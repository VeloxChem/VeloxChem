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

#include "NuclearPotentialGeom010PrimRecFH.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_fh(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_fh,
                                        const size_t              idx_npot_geom_010_0_ph,
                                        const size_t              idx_npot_geom_010_1_ph,
                                        const size_t              idx_npot_geom_010_0_dg,
                                        const size_t              idx_npot_geom_010_1_dg,
                                        const size_t              idx_npot_1_dh,
                                        const size_t              idx_npot_geom_010_0_dh,
                                        const size_t              idx_npot_geom_010_1_dh,
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

    // Set up components of auxiliary buffer : PH

    auto ta1_x_x_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph);

    auto ta1_x_x_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 1);

    auto ta1_x_x_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 2);

    auto ta1_x_x_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 3);

    auto ta1_x_x_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 4);

    auto ta1_x_x_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 5);

    auto ta1_x_x_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 6);

    auto ta1_x_x_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 7);

    auto ta1_x_x_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 8);

    auto ta1_x_x_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 9);

    auto ta1_x_x_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 10);

    auto ta1_x_x_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 11);

    auto ta1_x_x_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 12);

    auto ta1_x_x_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 13);

    auto ta1_x_x_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 14);

    auto ta1_x_x_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 15);

    auto ta1_x_x_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 16);

    auto ta1_x_x_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 17);

    auto ta1_x_x_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 18);

    auto ta1_x_x_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 19);

    auto ta1_x_x_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 20);

    auto ta1_x_y_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 21);

    auto ta1_x_y_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 22);

    auto ta1_x_y_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 23);

    auto ta1_x_y_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 24);

    auto ta1_x_y_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 25);

    auto ta1_x_y_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 26);

    auto ta1_x_y_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 27);

    auto ta1_x_y_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 28);

    auto ta1_x_y_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 29);

    auto ta1_x_y_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 30);

    auto ta1_x_y_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 31);

    auto ta1_x_y_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 32);

    auto ta1_x_y_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 33);

    auto ta1_x_y_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 34);

    auto ta1_x_y_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 35);

    auto ta1_x_y_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 36);

    auto ta1_x_y_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 37);

    auto ta1_x_y_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 38);

    auto ta1_x_y_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 39);

    auto ta1_x_y_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 40);

    auto ta1_x_y_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 41);

    auto ta1_x_z_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 42);

    auto ta1_x_z_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 43);

    auto ta1_x_z_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 44);

    auto ta1_x_z_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 45);

    auto ta1_x_z_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 46);

    auto ta1_x_z_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 47);

    auto ta1_x_z_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 48);

    auto ta1_x_z_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 49);

    auto ta1_x_z_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 50);

    auto ta1_x_z_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 51);

    auto ta1_x_z_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 52);

    auto ta1_x_z_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 53);

    auto ta1_x_z_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 54);

    auto ta1_x_z_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 55);

    auto ta1_x_z_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 56);

    auto ta1_x_z_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 57);

    auto ta1_x_z_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 58);

    auto ta1_x_z_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 59);

    auto ta1_x_z_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 60);

    auto ta1_x_z_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 61);

    auto ta1_x_z_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 62);

    auto ta1_y_x_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 63);

    auto ta1_y_x_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 64);

    auto ta1_y_x_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 65);

    auto ta1_y_x_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 66);

    auto ta1_y_x_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 67);

    auto ta1_y_x_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 68);

    auto ta1_y_x_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 69);

    auto ta1_y_x_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 70);

    auto ta1_y_x_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 71);

    auto ta1_y_x_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 72);

    auto ta1_y_x_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 73);

    auto ta1_y_x_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 74);

    auto ta1_y_x_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 75);

    auto ta1_y_x_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 76);

    auto ta1_y_x_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 77);

    auto ta1_y_x_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 78);

    auto ta1_y_x_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 79);

    auto ta1_y_x_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 80);

    auto ta1_y_x_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 81);

    auto ta1_y_x_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 82);

    auto ta1_y_x_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 83);

    auto ta1_y_y_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 84);

    auto ta1_y_y_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 85);

    auto ta1_y_y_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 86);

    auto ta1_y_y_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 87);

    auto ta1_y_y_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 88);

    auto ta1_y_y_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 89);

    auto ta1_y_y_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 90);

    auto ta1_y_y_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 91);

    auto ta1_y_y_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 92);

    auto ta1_y_y_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 93);

    auto ta1_y_y_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 94);

    auto ta1_y_y_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 95);

    auto ta1_y_y_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 96);

    auto ta1_y_y_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 97);

    auto ta1_y_y_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 98);

    auto ta1_y_y_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 99);

    auto ta1_y_y_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 100);

    auto ta1_y_y_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 101);

    auto ta1_y_y_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 102);

    auto ta1_y_y_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 103);

    auto ta1_y_y_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 104);

    auto ta1_y_z_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 105);

    auto ta1_y_z_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 106);

    auto ta1_y_z_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 107);

    auto ta1_y_z_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 108);

    auto ta1_y_z_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 109);

    auto ta1_y_z_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 110);

    auto ta1_y_z_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 111);

    auto ta1_y_z_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 112);

    auto ta1_y_z_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 113);

    auto ta1_y_z_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 114);

    auto ta1_y_z_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 115);

    auto ta1_y_z_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 116);

    auto ta1_y_z_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 117);

    auto ta1_y_z_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 118);

    auto ta1_y_z_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 119);

    auto ta1_y_z_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 120);

    auto ta1_y_z_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 121);

    auto ta1_y_z_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 122);

    auto ta1_y_z_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 123);

    auto ta1_y_z_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 124);

    auto ta1_y_z_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 125);

    auto ta1_z_x_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 126);

    auto ta1_z_x_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 127);

    auto ta1_z_x_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 128);

    auto ta1_z_x_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 129);

    auto ta1_z_x_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 130);

    auto ta1_z_x_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 131);

    auto ta1_z_x_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 132);

    auto ta1_z_x_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 133);

    auto ta1_z_x_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 134);

    auto ta1_z_x_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 135);

    auto ta1_z_x_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 136);

    auto ta1_z_x_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 137);

    auto ta1_z_x_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 138);

    auto ta1_z_x_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 139);

    auto ta1_z_x_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 140);

    auto ta1_z_x_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 141);

    auto ta1_z_x_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 142);

    auto ta1_z_x_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 143);

    auto ta1_z_x_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 144);

    auto ta1_z_x_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 145);

    auto ta1_z_x_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 146);

    auto ta1_z_y_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 147);

    auto ta1_z_y_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 148);

    auto ta1_z_y_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 149);

    auto ta1_z_y_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 150);

    auto ta1_z_y_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 151);

    auto ta1_z_y_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 152);

    auto ta1_z_y_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 153);

    auto ta1_z_y_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 154);

    auto ta1_z_y_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 155);

    auto ta1_z_y_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 156);

    auto ta1_z_y_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 157);

    auto ta1_z_y_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 158);

    auto ta1_z_y_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 159);

    auto ta1_z_y_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 160);

    auto ta1_z_y_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 161);

    auto ta1_z_y_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 162);

    auto ta1_z_y_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 163);

    auto ta1_z_y_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 164);

    auto ta1_z_y_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 165);

    auto ta1_z_y_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 166);

    auto ta1_z_y_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 167);

    auto ta1_z_z_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 168);

    auto ta1_z_z_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 169);

    auto ta1_z_z_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 170);

    auto ta1_z_z_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 171);

    auto ta1_z_z_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 172);

    auto ta1_z_z_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 173);

    auto ta1_z_z_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 174);

    auto ta1_z_z_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 175);

    auto ta1_z_z_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 176);

    auto ta1_z_z_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 177);

    auto ta1_z_z_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 178);

    auto ta1_z_z_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 179);

    auto ta1_z_z_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 180);

    auto ta1_z_z_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 181);

    auto ta1_z_z_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 182);

    auto ta1_z_z_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 183);

    auto ta1_z_z_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 184);

    auto ta1_z_z_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 185);

    auto ta1_z_z_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 186);

    auto ta1_z_z_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 187);

    auto ta1_z_z_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 188);

    // Set up components of auxiliary buffer : PH

    auto ta1_x_x_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph);

    auto ta1_x_x_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 1);

    auto ta1_x_x_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 2);

    auto ta1_x_x_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 3);

    auto ta1_x_x_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 4);

    auto ta1_x_x_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 5);

    auto ta1_x_x_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 6);

    auto ta1_x_x_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 7);

    auto ta1_x_x_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 8);

    auto ta1_x_x_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 9);

    auto ta1_x_x_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 10);

    auto ta1_x_x_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 11);

    auto ta1_x_x_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 12);

    auto ta1_x_x_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 13);

    auto ta1_x_x_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 14);

    auto ta1_x_x_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 15);

    auto ta1_x_x_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 16);

    auto ta1_x_x_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 17);

    auto ta1_x_x_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 18);

    auto ta1_x_x_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 19);

    auto ta1_x_x_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 20);

    auto ta1_x_y_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 21);

    auto ta1_x_y_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 22);

    auto ta1_x_y_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 23);

    auto ta1_x_y_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 24);

    auto ta1_x_y_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 25);

    auto ta1_x_y_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 26);

    auto ta1_x_y_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 27);

    auto ta1_x_y_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 28);

    auto ta1_x_y_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 29);

    auto ta1_x_y_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 30);

    auto ta1_x_y_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 31);

    auto ta1_x_y_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 32);

    auto ta1_x_y_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 33);

    auto ta1_x_y_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 34);

    auto ta1_x_y_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 35);

    auto ta1_x_y_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 36);

    auto ta1_x_y_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 37);

    auto ta1_x_y_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 38);

    auto ta1_x_y_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 39);

    auto ta1_x_y_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 40);

    auto ta1_x_y_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 41);

    auto ta1_x_z_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 42);

    auto ta1_x_z_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 43);

    auto ta1_x_z_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 44);

    auto ta1_x_z_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 45);

    auto ta1_x_z_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 46);

    auto ta1_x_z_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 47);

    auto ta1_x_z_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 48);

    auto ta1_x_z_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 49);

    auto ta1_x_z_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 50);

    auto ta1_x_z_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 51);

    auto ta1_x_z_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 52);

    auto ta1_x_z_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 53);

    auto ta1_x_z_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 54);

    auto ta1_x_z_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 55);

    auto ta1_x_z_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 56);

    auto ta1_x_z_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 57);

    auto ta1_x_z_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 58);

    auto ta1_x_z_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 59);

    auto ta1_x_z_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 60);

    auto ta1_x_z_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 61);

    auto ta1_x_z_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 62);

    auto ta1_y_x_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 63);

    auto ta1_y_x_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 64);

    auto ta1_y_x_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 65);

    auto ta1_y_x_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 66);

    auto ta1_y_x_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 67);

    auto ta1_y_x_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 68);

    auto ta1_y_x_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 69);

    auto ta1_y_x_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 70);

    auto ta1_y_x_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 71);

    auto ta1_y_x_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 72);

    auto ta1_y_x_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 73);

    auto ta1_y_x_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 74);

    auto ta1_y_x_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 75);

    auto ta1_y_x_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 76);

    auto ta1_y_x_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 77);

    auto ta1_y_x_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 78);

    auto ta1_y_x_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 79);

    auto ta1_y_x_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 80);

    auto ta1_y_x_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 81);

    auto ta1_y_x_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 82);

    auto ta1_y_x_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 83);

    auto ta1_y_y_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 84);

    auto ta1_y_y_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 85);

    auto ta1_y_y_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 86);

    auto ta1_y_y_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 87);

    auto ta1_y_y_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 88);

    auto ta1_y_y_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 89);

    auto ta1_y_y_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 90);

    auto ta1_y_y_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 91);

    auto ta1_y_y_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 92);

    auto ta1_y_y_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 93);

    auto ta1_y_y_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 94);

    auto ta1_y_y_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 95);

    auto ta1_y_y_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 96);

    auto ta1_y_y_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 97);

    auto ta1_y_y_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 98);

    auto ta1_y_y_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 99);

    auto ta1_y_y_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 100);

    auto ta1_y_y_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 101);

    auto ta1_y_y_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 102);

    auto ta1_y_y_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 103);

    auto ta1_y_y_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 104);

    auto ta1_y_z_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 105);

    auto ta1_y_z_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 106);

    auto ta1_y_z_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 107);

    auto ta1_y_z_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 108);

    auto ta1_y_z_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 109);

    auto ta1_y_z_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 110);

    auto ta1_y_z_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 111);

    auto ta1_y_z_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 112);

    auto ta1_y_z_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 113);

    auto ta1_y_z_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 114);

    auto ta1_y_z_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 115);

    auto ta1_y_z_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 116);

    auto ta1_y_z_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 117);

    auto ta1_y_z_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 118);

    auto ta1_y_z_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 119);

    auto ta1_y_z_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 120);

    auto ta1_y_z_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 121);

    auto ta1_y_z_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 122);

    auto ta1_y_z_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 123);

    auto ta1_y_z_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 124);

    auto ta1_y_z_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 125);

    auto ta1_z_x_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 126);

    auto ta1_z_x_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 127);

    auto ta1_z_x_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 128);

    auto ta1_z_x_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 129);

    auto ta1_z_x_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 130);

    auto ta1_z_x_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 131);

    auto ta1_z_x_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 132);

    auto ta1_z_x_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 133);

    auto ta1_z_x_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 134);

    auto ta1_z_x_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 135);

    auto ta1_z_x_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 136);

    auto ta1_z_x_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 137);

    auto ta1_z_x_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 138);

    auto ta1_z_x_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 139);

    auto ta1_z_x_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 140);

    auto ta1_z_x_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 141);

    auto ta1_z_x_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 142);

    auto ta1_z_x_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 143);

    auto ta1_z_x_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 144);

    auto ta1_z_x_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 145);

    auto ta1_z_x_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 146);

    auto ta1_z_y_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 147);

    auto ta1_z_y_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 148);

    auto ta1_z_y_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 149);

    auto ta1_z_y_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 150);

    auto ta1_z_y_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 151);

    auto ta1_z_y_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 152);

    auto ta1_z_y_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 153);

    auto ta1_z_y_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 154);

    auto ta1_z_y_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 155);

    auto ta1_z_y_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 156);

    auto ta1_z_y_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 157);

    auto ta1_z_y_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 158);

    auto ta1_z_y_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 159);

    auto ta1_z_y_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 160);

    auto ta1_z_y_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 161);

    auto ta1_z_y_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 162);

    auto ta1_z_y_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 163);

    auto ta1_z_y_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 164);

    auto ta1_z_y_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 165);

    auto ta1_z_y_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 166);

    auto ta1_z_y_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 167);

    auto ta1_z_z_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 168);

    auto ta1_z_z_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 169);

    auto ta1_z_z_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 170);

    auto ta1_z_z_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 171);

    auto ta1_z_z_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 172);

    auto ta1_z_z_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 173);

    auto ta1_z_z_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 174);

    auto ta1_z_z_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 175);

    auto ta1_z_z_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 176);

    auto ta1_z_z_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 177);

    auto ta1_z_z_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 178);

    auto ta1_z_z_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 179);

    auto ta1_z_z_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 180);

    auto ta1_z_z_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 181);

    auto ta1_z_z_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 182);

    auto ta1_z_z_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 183);

    auto ta1_z_z_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 184);

    auto ta1_z_z_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 185);

    auto ta1_z_z_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 186);

    auto ta1_z_z_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 187);

    auto ta1_z_z_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 188);

    // Set up components of auxiliary buffer : DG

    auto ta1_x_xx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg);

    auto ta1_x_xx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 1);

    auto ta1_x_xx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 2);

    auto ta1_x_xx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 3);

    auto ta1_x_xx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 4);

    auto ta1_x_xx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 5);

    auto ta1_x_xx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 6);

    auto ta1_x_xx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 7);

    auto ta1_x_xx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 8);

    auto ta1_x_xx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 9);

    auto ta1_x_xx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 10);

    auto ta1_x_xx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 11);

    auto ta1_x_xx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 12);

    auto ta1_x_xx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 13);

    auto ta1_x_xx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 14);

    auto ta1_x_xz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 32);

    auto ta1_x_xz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 34);

    auto ta1_x_xz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 35);

    auto ta1_x_xz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 37);

    auto ta1_x_xz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 38);

    auto ta1_x_xz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 39);

    auto ta1_x_yy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 45);

    auto ta1_x_yy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 46);

    auto ta1_x_yy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 47);

    auto ta1_x_yy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 48);

    auto ta1_x_yy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 49);

    auto ta1_x_yy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 50);

    auto ta1_x_yy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 51);

    auto ta1_x_yy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 52);

    auto ta1_x_yy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 53);

    auto ta1_x_yy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 54);

    auto ta1_x_yy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 55);

    auto ta1_x_yy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 56);

    auto ta1_x_yy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 57);

    auto ta1_x_yy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 58);

    auto ta1_x_yy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 59);

    auto ta1_x_zz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 75);

    auto ta1_x_zz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 76);

    auto ta1_x_zz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 77);

    auto ta1_x_zz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 78);

    auto ta1_x_zz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 79);

    auto ta1_x_zz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 80);

    auto ta1_x_zz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 81);

    auto ta1_x_zz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 82);

    auto ta1_x_zz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 83);

    auto ta1_x_zz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 84);

    auto ta1_x_zz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 85);

    auto ta1_x_zz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 86);

    auto ta1_x_zz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 87);

    auto ta1_x_zz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 88);

    auto ta1_x_zz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 89);

    auto ta1_y_xx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 90);

    auto ta1_y_xx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 91);

    auto ta1_y_xx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 92);

    auto ta1_y_xx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 93);

    auto ta1_y_xx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 94);

    auto ta1_y_xx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 95);

    auto ta1_y_xx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 96);

    auto ta1_y_xx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 97);

    auto ta1_y_xx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 98);

    auto ta1_y_xx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 99);

    auto ta1_y_xx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 100);

    auto ta1_y_xx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 101);

    auto ta1_y_xx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 102);

    auto ta1_y_xx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 103);

    auto ta1_y_xx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 104);

    auto ta1_y_yy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 135);

    auto ta1_y_yy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 136);

    auto ta1_y_yy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 137);

    auto ta1_y_yy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 138);

    auto ta1_y_yy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 139);

    auto ta1_y_yy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 140);

    auto ta1_y_yy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 141);

    auto ta1_y_yy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 142);

    auto ta1_y_yy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 143);

    auto ta1_y_yy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 144);

    auto ta1_y_yy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 145);

    auto ta1_y_yy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 146);

    auto ta1_y_yy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 147);

    auto ta1_y_yy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 148);

    auto ta1_y_yy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 149);

    auto ta1_y_yz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 154);

    auto ta1_y_yz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 157);

    auto ta1_y_yz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 158);

    auto ta1_y_yz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 161);

    auto ta1_y_yz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 162);

    auto ta1_y_yz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 163);

    auto ta1_y_zz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 165);

    auto ta1_y_zz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 166);

    auto ta1_y_zz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 167);

    auto ta1_y_zz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 168);

    auto ta1_y_zz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 169);

    auto ta1_y_zz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 170);

    auto ta1_y_zz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 171);

    auto ta1_y_zz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 172);

    auto ta1_y_zz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 173);

    auto ta1_y_zz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 174);

    auto ta1_y_zz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 175);

    auto ta1_y_zz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 176);

    auto ta1_y_zz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 177);

    auto ta1_y_zz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 178);

    auto ta1_y_zz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 179);

    auto ta1_z_xx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 180);

    auto ta1_z_xx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 181);

    auto ta1_z_xx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 182);

    auto ta1_z_xx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 183);

    auto ta1_z_xx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 184);

    auto ta1_z_xx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 185);

    auto ta1_z_xx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 186);

    auto ta1_z_xx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 187);

    auto ta1_z_xx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 188);

    auto ta1_z_xx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 189);

    auto ta1_z_xx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 190);

    auto ta1_z_xx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 191);

    auto ta1_z_xx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 192);

    auto ta1_z_xx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 193);

    auto ta1_z_xx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 194);

    auto ta1_z_yy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 225);

    auto ta1_z_yy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 226);

    auto ta1_z_yy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 227);

    auto ta1_z_yy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 228);

    auto ta1_z_yy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 229);

    auto ta1_z_yy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 230);

    auto ta1_z_yy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 231);

    auto ta1_z_yy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 232);

    auto ta1_z_yy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 233);

    auto ta1_z_yy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 234);

    auto ta1_z_yy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 235);

    auto ta1_z_yy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 236);

    auto ta1_z_yy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 237);

    auto ta1_z_yy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 238);

    auto ta1_z_yy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 239);

    auto ta1_z_yz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 244);

    auto ta1_z_yz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 247);

    auto ta1_z_yz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 248);

    auto ta1_z_yz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 251);

    auto ta1_z_yz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 252);

    auto ta1_z_yz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 253);

    auto ta1_z_zz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 255);

    auto ta1_z_zz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 256);

    auto ta1_z_zz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 257);

    auto ta1_z_zz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 258);

    auto ta1_z_zz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 259);

    auto ta1_z_zz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 260);

    auto ta1_z_zz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 261);

    auto ta1_z_zz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 262);

    auto ta1_z_zz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 263);

    auto ta1_z_zz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 264);

    auto ta1_z_zz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 265);

    auto ta1_z_zz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 266);

    auto ta1_z_zz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 267);

    auto ta1_z_zz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 268);

    auto ta1_z_zz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 269);

    // Set up components of auxiliary buffer : DG

    auto ta1_x_xx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg);

    auto ta1_x_xx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 1);

    auto ta1_x_xx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 2);

    auto ta1_x_xx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 3);

    auto ta1_x_xx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 4);

    auto ta1_x_xx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 5);

    auto ta1_x_xx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 6);

    auto ta1_x_xx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 7);

    auto ta1_x_xx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 8);

    auto ta1_x_xx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 9);

    auto ta1_x_xx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 10);

    auto ta1_x_xx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 11);

    auto ta1_x_xx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 12);

    auto ta1_x_xx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 13);

    auto ta1_x_xx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 14);

    auto ta1_x_xz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 32);

    auto ta1_x_xz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 34);

    auto ta1_x_xz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 35);

    auto ta1_x_xz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 37);

    auto ta1_x_xz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 38);

    auto ta1_x_xz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 39);

    auto ta1_x_yy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 45);

    auto ta1_x_yy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 46);

    auto ta1_x_yy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 47);

    auto ta1_x_yy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 48);

    auto ta1_x_yy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 49);

    auto ta1_x_yy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 50);

    auto ta1_x_yy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 51);

    auto ta1_x_yy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 52);

    auto ta1_x_yy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 53);

    auto ta1_x_yy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 54);

    auto ta1_x_yy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 55);

    auto ta1_x_yy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 56);

    auto ta1_x_yy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 57);

    auto ta1_x_yy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 58);

    auto ta1_x_yy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 59);

    auto ta1_x_zz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 75);

    auto ta1_x_zz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 76);

    auto ta1_x_zz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 77);

    auto ta1_x_zz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 78);

    auto ta1_x_zz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 79);

    auto ta1_x_zz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 80);

    auto ta1_x_zz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 81);

    auto ta1_x_zz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 82);

    auto ta1_x_zz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 83);

    auto ta1_x_zz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 84);

    auto ta1_x_zz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 85);

    auto ta1_x_zz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 86);

    auto ta1_x_zz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 87);

    auto ta1_x_zz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 88);

    auto ta1_x_zz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 89);

    auto ta1_y_xx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 90);

    auto ta1_y_xx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 91);

    auto ta1_y_xx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 92);

    auto ta1_y_xx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 93);

    auto ta1_y_xx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 94);

    auto ta1_y_xx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 95);

    auto ta1_y_xx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 96);

    auto ta1_y_xx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 97);

    auto ta1_y_xx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 98);

    auto ta1_y_xx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 99);

    auto ta1_y_xx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 100);

    auto ta1_y_xx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 101);

    auto ta1_y_xx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 102);

    auto ta1_y_xx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 103);

    auto ta1_y_xx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 104);

    auto ta1_y_yy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 135);

    auto ta1_y_yy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 136);

    auto ta1_y_yy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 137);

    auto ta1_y_yy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 138);

    auto ta1_y_yy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 139);

    auto ta1_y_yy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 140);

    auto ta1_y_yy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 141);

    auto ta1_y_yy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 142);

    auto ta1_y_yy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 143);

    auto ta1_y_yy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 144);

    auto ta1_y_yy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 145);

    auto ta1_y_yy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 146);

    auto ta1_y_yy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 147);

    auto ta1_y_yy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 148);

    auto ta1_y_yy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 149);

    auto ta1_y_yz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 154);

    auto ta1_y_yz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 157);

    auto ta1_y_yz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 158);

    auto ta1_y_yz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 161);

    auto ta1_y_yz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 162);

    auto ta1_y_yz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 163);

    auto ta1_y_zz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 165);

    auto ta1_y_zz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 166);

    auto ta1_y_zz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 167);

    auto ta1_y_zz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 168);

    auto ta1_y_zz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 169);

    auto ta1_y_zz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 170);

    auto ta1_y_zz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 171);

    auto ta1_y_zz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 172);

    auto ta1_y_zz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 173);

    auto ta1_y_zz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 174);

    auto ta1_y_zz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 175);

    auto ta1_y_zz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 176);

    auto ta1_y_zz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 177);

    auto ta1_y_zz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 178);

    auto ta1_y_zz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 179);

    auto ta1_z_xx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 180);

    auto ta1_z_xx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 181);

    auto ta1_z_xx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 182);

    auto ta1_z_xx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 183);

    auto ta1_z_xx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 184);

    auto ta1_z_xx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 185);

    auto ta1_z_xx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 186);

    auto ta1_z_xx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 187);

    auto ta1_z_xx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 188);

    auto ta1_z_xx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 189);

    auto ta1_z_xx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 190);

    auto ta1_z_xx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 191);

    auto ta1_z_xx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 192);

    auto ta1_z_xx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 193);

    auto ta1_z_xx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 194);

    auto ta1_z_yy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 225);

    auto ta1_z_yy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 226);

    auto ta1_z_yy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 227);

    auto ta1_z_yy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 228);

    auto ta1_z_yy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 229);

    auto ta1_z_yy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 230);

    auto ta1_z_yy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 231);

    auto ta1_z_yy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 232);

    auto ta1_z_yy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 233);

    auto ta1_z_yy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 234);

    auto ta1_z_yy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 235);

    auto ta1_z_yy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 236);

    auto ta1_z_yy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 237);

    auto ta1_z_yy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 238);

    auto ta1_z_yy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 239);

    auto ta1_z_yz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 244);

    auto ta1_z_yz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 247);

    auto ta1_z_yz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 248);

    auto ta1_z_yz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 251);

    auto ta1_z_yz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 252);

    auto ta1_z_yz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 253);

    auto ta1_z_zz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 255);

    auto ta1_z_zz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 256);

    auto ta1_z_zz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 257);

    auto ta1_z_zz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 258);

    auto ta1_z_zz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 259);

    auto ta1_z_zz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 260);

    auto ta1_z_zz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 261);

    auto ta1_z_zz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 262);

    auto ta1_z_zz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 263);

    auto ta1_z_zz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 264);

    auto ta1_z_zz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 265);

    auto ta1_z_zz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 266);

    auto ta1_z_zz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 267);

    auto ta1_z_zz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 268);

    auto ta1_z_zz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 269);

    // Set up components of auxiliary buffer : DH

    auto ta_xx_xxxxx_1 = pbuffer.data(idx_npot_1_dh);

    auto ta_xx_xxxxy_1 = pbuffer.data(idx_npot_1_dh + 1);

    auto ta_xx_xxxxz_1 = pbuffer.data(idx_npot_1_dh + 2);

    auto ta_xx_xxxyy_1 = pbuffer.data(idx_npot_1_dh + 3);

    auto ta_xx_xxxyz_1 = pbuffer.data(idx_npot_1_dh + 4);

    auto ta_xx_xxxzz_1 = pbuffer.data(idx_npot_1_dh + 5);

    auto ta_xx_xxyyy_1 = pbuffer.data(idx_npot_1_dh + 6);

    auto ta_xx_xxyyz_1 = pbuffer.data(idx_npot_1_dh + 7);

    auto ta_xx_xxyzz_1 = pbuffer.data(idx_npot_1_dh + 8);

    auto ta_xx_xxzzz_1 = pbuffer.data(idx_npot_1_dh + 9);

    auto ta_xx_xyyyy_1 = pbuffer.data(idx_npot_1_dh + 10);

    auto ta_xx_xyyyz_1 = pbuffer.data(idx_npot_1_dh + 11);

    auto ta_xx_xyyzz_1 = pbuffer.data(idx_npot_1_dh + 12);

    auto ta_xx_xyzzz_1 = pbuffer.data(idx_npot_1_dh + 13);

    auto ta_xx_xzzzz_1 = pbuffer.data(idx_npot_1_dh + 14);

    auto ta_xx_yyyyy_1 = pbuffer.data(idx_npot_1_dh + 15);

    auto ta_xx_yyyyz_1 = pbuffer.data(idx_npot_1_dh + 16);

    auto ta_xx_yyyzz_1 = pbuffer.data(idx_npot_1_dh + 17);

    auto ta_xx_yyzzz_1 = pbuffer.data(idx_npot_1_dh + 18);

    auto ta_xx_yzzzz_1 = pbuffer.data(idx_npot_1_dh + 19);

    auto ta_xx_zzzzz_1 = pbuffer.data(idx_npot_1_dh + 20);

    auto ta_xy_xxxxy_1 = pbuffer.data(idx_npot_1_dh + 22);

    auto ta_xy_xxxyy_1 = pbuffer.data(idx_npot_1_dh + 24);

    auto ta_xy_xxyyy_1 = pbuffer.data(idx_npot_1_dh + 27);

    auto ta_xy_xyyyy_1 = pbuffer.data(idx_npot_1_dh + 31);

    auto ta_xz_xxxxz_1 = pbuffer.data(idx_npot_1_dh + 44);

    auto ta_xz_xxxzz_1 = pbuffer.data(idx_npot_1_dh + 47);

    auto ta_xz_xxzzz_1 = pbuffer.data(idx_npot_1_dh + 51);

    auto ta_xz_xzzzz_1 = pbuffer.data(idx_npot_1_dh + 56);

    auto ta_yy_xxxxx_1 = pbuffer.data(idx_npot_1_dh + 63);

    auto ta_yy_xxxxy_1 = pbuffer.data(idx_npot_1_dh + 64);

    auto ta_yy_xxxxz_1 = pbuffer.data(idx_npot_1_dh + 65);

    auto ta_yy_xxxyy_1 = pbuffer.data(idx_npot_1_dh + 66);

    auto ta_yy_xxxyz_1 = pbuffer.data(idx_npot_1_dh + 67);

    auto ta_yy_xxxzz_1 = pbuffer.data(idx_npot_1_dh + 68);

    auto ta_yy_xxyyy_1 = pbuffer.data(idx_npot_1_dh + 69);

    auto ta_yy_xxyyz_1 = pbuffer.data(idx_npot_1_dh + 70);

    auto ta_yy_xxyzz_1 = pbuffer.data(idx_npot_1_dh + 71);

    auto ta_yy_xxzzz_1 = pbuffer.data(idx_npot_1_dh + 72);

    auto ta_yy_xyyyy_1 = pbuffer.data(idx_npot_1_dh + 73);

    auto ta_yy_xyyyz_1 = pbuffer.data(idx_npot_1_dh + 74);

    auto ta_yy_xyyzz_1 = pbuffer.data(idx_npot_1_dh + 75);

    auto ta_yy_xyzzz_1 = pbuffer.data(idx_npot_1_dh + 76);

    auto ta_yy_xzzzz_1 = pbuffer.data(idx_npot_1_dh + 77);

    auto ta_yy_yyyyy_1 = pbuffer.data(idx_npot_1_dh + 78);

    auto ta_yy_yyyyz_1 = pbuffer.data(idx_npot_1_dh + 79);

    auto ta_yy_yyyzz_1 = pbuffer.data(idx_npot_1_dh + 80);

    auto ta_yy_yyzzz_1 = pbuffer.data(idx_npot_1_dh + 81);

    auto ta_yy_yzzzz_1 = pbuffer.data(idx_npot_1_dh + 82);

    auto ta_yy_zzzzz_1 = pbuffer.data(idx_npot_1_dh + 83);

    auto ta_yz_yyyyz_1 = pbuffer.data(idx_npot_1_dh + 100);

    auto ta_yz_yyyzz_1 = pbuffer.data(idx_npot_1_dh + 101);

    auto ta_yz_yyzzz_1 = pbuffer.data(idx_npot_1_dh + 102);

    auto ta_yz_yzzzz_1 = pbuffer.data(idx_npot_1_dh + 103);

    auto ta_zz_xxxxx_1 = pbuffer.data(idx_npot_1_dh + 105);

    auto ta_zz_xxxxy_1 = pbuffer.data(idx_npot_1_dh + 106);

    auto ta_zz_xxxxz_1 = pbuffer.data(idx_npot_1_dh + 107);

    auto ta_zz_xxxyy_1 = pbuffer.data(idx_npot_1_dh + 108);

    auto ta_zz_xxxyz_1 = pbuffer.data(idx_npot_1_dh + 109);

    auto ta_zz_xxxzz_1 = pbuffer.data(idx_npot_1_dh + 110);

    auto ta_zz_xxyyy_1 = pbuffer.data(idx_npot_1_dh + 111);

    auto ta_zz_xxyyz_1 = pbuffer.data(idx_npot_1_dh + 112);

    auto ta_zz_xxyzz_1 = pbuffer.data(idx_npot_1_dh + 113);

    auto ta_zz_xxzzz_1 = pbuffer.data(idx_npot_1_dh + 114);

    auto ta_zz_xyyyy_1 = pbuffer.data(idx_npot_1_dh + 115);

    auto ta_zz_xyyyz_1 = pbuffer.data(idx_npot_1_dh + 116);

    auto ta_zz_xyyzz_1 = pbuffer.data(idx_npot_1_dh + 117);

    auto ta_zz_xyzzz_1 = pbuffer.data(idx_npot_1_dh + 118);

    auto ta_zz_xzzzz_1 = pbuffer.data(idx_npot_1_dh + 119);

    auto ta_zz_yyyyy_1 = pbuffer.data(idx_npot_1_dh + 120);

    auto ta_zz_yyyyz_1 = pbuffer.data(idx_npot_1_dh + 121);

    auto ta_zz_yyyzz_1 = pbuffer.data(idx_npot_1_dh + 122);

    auto ta_zz_yyzzz_1 = pbuffer.data(idx_npot_1_dh + 123);

    auto ta_zz_yzzzz_1 = pbuffer.data(idx_npot_1_dh + 124);

    auto ta_zz_zzzzz_1 = pbuffer.data(idx_npot_1_dh + 125);

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

    auto ta1_x_xy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 22);

    auto ta1_x_xy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 23);

    auto ta1_x_xy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 24);

    auto ta1_x_xy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 26);

    auto ta1_x_xy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 27);

    auto ta1_x_xy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 30);

    auto ta1_x_xy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 31);

    auto ta1_x_xy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 35);

    auto ta1_x_xy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 36);

    auto ta1_x_xz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 42);

    auto ta1_x_xz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 43);

    auto ta1_x_xz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 44);

    auto ta1_x_xz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 45);

    auto ta1_x_xz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 46);

    auto ta1_x_xz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 47);

    auto ta1_x_xz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 48);

    auto ta1_x_xz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 49);

    auto ta1_x_xz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 50);

    auto ta1_x_xz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 51);

    auto ta1_x_xz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 52);

    auto ta1_x_xz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 53);

    auto ta1_x_xz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 54);

    auto ta1_x_xz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 55);

    auto ta1_x_xz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 56);

    auto ta1_x_xz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 62);

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

    auto ta1_x_yz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 100);

    auto ta1_x_yz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 101);

    auto ta1_x_yz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 102);

    auto ta1_x_yz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 103);

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

    auto ta1_y_xy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 147);

    auto ta1_y_xy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 148);

    auto ta1_y_xy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 150);

    auto ta1_y_xy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 153);

    auto ta1_y_xy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 157);

    auto ta1_y_xy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 162);

    auto ta1_y_xy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 163);

    auto ta1_y_xy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 164);

    auto ta1_y_xy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 165);

    auto ta1_y_xy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 166);

    auto ta1_y_xz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 170);

    auto ta1_y_xz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 173);

    auto ta1_y_xz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 177);

    auto ta1_y_xz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 182);

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

    auto ta1_y_yz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 214);

    auto ta1_y_yz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 216);

    auto ta1_y_yz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 217);

    auto ta1_y_yz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 218);

    auto ta1_y_yz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 220);

    auto ta1_y_yz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 221);

    auto ta1_y_yz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 222);

    auto ta1_y_yz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 223);

    auto ta1_y_yz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 225);

    auto ta1_y_yz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 226);

    auto ta1_y_yz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 227);

    auto ta1_y_yz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 228);

    auto ta1_y_yz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 229);

    auto ta1_y_yz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 230);

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

    auto ta1_z_xy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 274);

    auto ta1_z_xy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 276);

    auto ta1_z_xy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 279);

    auto ta1_z_xy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 283);

    auto ta1_z_xy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 288);

    auto ta1_z_xy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 289);

    auto ta1_z_xy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 290);

    auto ta1_z_xy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 291);

    auto ta1_z_xy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 292);

    auto ta1_z_xz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 294);

    auto ta1_z_xz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 296);

    auto ta1_z_xz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 299);

    auto ta1_z_xz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 303);

    auto ta1_z_xz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 308);

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

    auto ta1_z_yz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 340);

    auto ta1_z_yz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 341);

    auto ta1_z_yz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 343);

    auto ta1_z_yz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 344);

    auto ta1_z_yz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 345);

    auto ta1_z_yz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 347);

    auto ta1_z_yz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 348);

    auto ta1_z_yz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 349);

    auto ta1_z_yz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 350);

    auto ta1_z_yz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 351);

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

    auto ta1_x_xy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 22);

    auto ta1_x_xy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 23);

    auto ta1_x_xy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 24);

    auto ta1_x_xy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 26);

    auto ta1_x_xy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 27);

    auto ta1_x_xy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 30);

    auto ta1_x_xy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 31);

    auto ta1_x_xy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 35);

    auto ta1_x_xy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 36);

    auto ta1_x_xz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 42);

    auto ta1_x_xz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 43);

    auto ta1_x_xz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 44);

    auto ta1_x_xz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 45);

    auto ta1_x_xz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 46);

    auto ta1_x_xz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 47);

    auto ta1_x_xz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 48);

    auto ta1_x_xz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 49);

    auto ta1_x_xz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 50);

    auto ta1_x_xz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 51);

    auto ta1_x_xz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 52);

    auto ta1_x_xz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 53);

    auto ta1_x_xz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 54);

    auto ta1_x_xz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 55);

    auto ta1_x_xz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 56);

    auto ta1_x_xz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 62);

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

    auto ta1_x_yz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 100);

    auto ta1_x_yz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 101);

    auto ta1_x_yz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 102);

    auto ta1_x_yz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 103);

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

    auto ta1_y_xy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 147);

    auto ta1_y_xy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 148);

    auto ta1_y_xy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 150);

    auto ta1_y_xy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 153);

    auto ta1_y_xy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 157);

    auto ta1_y_xy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 162);

    auto ta1_y_xy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 163);

    auto ta1_y_xy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 164);

    auto ta1_y_xy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 165);

    auto ta1_y_xy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 166);

    auto ta1_y_xz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 170);

    auto ta1_y_xz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 173);

    auto ta1_y_xz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 177);

    auto ta1_y_xz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 182);

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

    auto ta1_y_yz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 214);

    auto ta1_y_yz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 216);

    auto ta1_y_yz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 217);

    auto ta1_y_yz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 218);

    auto ta1_y_yz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 220);

    auto ta1_y_yz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 221);

    auto ta1_y_yz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 222);

    auto ta1_y_yz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 223);

    auto ta1_y_yz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 225);

    auto ta1_y_yz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 226);

    auto ta1_y_yz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 227);

    auto ta1_y_yz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 228);

    auto ta1_y_yz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 229);

    auto ta1_y_yz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 230);

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

    auto ta1_z_xy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 274);

    auto ta1_z_xy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 276);

    auto ta1_z_xy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 279);

    auto ta1_z_xy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 283);

    auto ta1_z_xy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 288);

    auto ta1_z_xy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 289);

    auto ta1_z_xy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 290);

    auto ta1_z_xy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 291);

    auto ta1_z_xy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 292);

    auto ta1_z_xz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 294);

    auto ta1_z_xz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 296);

    auto ta1_z_xz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 299);

    auto ta1_z_xz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 303);

    auto ta1_z_xz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 308);

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

    auto ta1_z_yz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 340);

    auto ta1_z_yz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 341);

    auto ta1_z_yz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 343);

    auto ta1_z_yz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 344);

    auto ta1_z_yz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 345);

    auto ta1_z_yz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 347);

    auto ta1_z_yz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 348);

    auto ta1_z_yz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 349);

    auto ta1_z_yz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 350);

    auto ta1_z_yz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 351);

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

    // Set up 0-21 components of targeted buffer : FH

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

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_x_x_xxxxx_0,   \
                             ta1_x_x_xxxxx_1,   \
                             ta1_x_x_xxxxy_0,   \
                             ta1_x_x_xxxxy_1,   \
                             ta1_x_x_xxxxz_0,   \
                             ta1_x_x_xxxxz_1,   \
                             ta1_x_x_xxxyy_0,   \
                             ta1_x_x_xxxyy_1,   \
                             ta1_x_x_xxxyz_0,   \
                             ta1_x_x_xxxyz_1,   \
                             ta1_x_x_xxxzz_0,   \
                             ta1_x_x_xxxzz_1,   \
                             ta1_x_x_xxyyy_0,   \
                             ta1_x_x_xxyyy_1,   \
                             ta1_x_x_xxyyz_0,   \
                             ta1_x_x_xxyyz_1,   \
                             ta1_x_x_xxyzz_0,   \
                             ta1_x_x_xxyzz_1,   \
                             ta1_x_x_xxzzz_0,   \
                             ta1_x_x_xxzzz_1,   \
                             ta1_x_x_xyyyy_0,   \
                             ta1_x_x_xyyyy_1,   \
                             ta1_x_x_xyyyz_0,   \
                             ta1_x_x_xyyyz_1,   \
                             ta1_x_x_xyyzz_0,   \
                             ta1_x_x_xyyzz_1,   \
                             ta1_x_x_xyzzz_0,   \
                             ta1_x_x_xyzzz_1,   \
                             ta1_x_x_xzzzz_0,   \
                             ta1_x_x_xzzzz_1,   \
                             ta1_x_x_yyyyy_0,   \
                             ta1_x_x_yyyyy_1,   \
                             ta1_x_x_yyyyz_0,   \
                             ta1_x_x_yyyyz_1,   \
                             ta1_x_x_yyyzz_0,   \
                             ta1_x_x_yyyzz_1,   \
                             ta1_x_x_yyzzz_0,   \
                             ta1_x_x_yyzzz_1,   \
                             ta1_x_x_yzzzz_0,   \
                             ta1_x_x_yzzzz_1,   \
                             ta1_x_x_zzzzz_0,   \
                             ta1_x_x_zzzzz_1,   \
                             ta1_x_xx_xxxx_0,   \
                             ta1_x_xx_xxxx_1,   \
                             ta1_x_xx_xxxxx_0,  \
                             ta1_x_xx_xxxxx_1,  \
                             ta1_x_xx_xxxxy_0,  \
                             ta1_x_xx_xxxxy_1,  \
                             ta1_x_xx_xxxxz_0,  \
                             ta1_x_xx_xxxxz_1,  \
                             ta1_x_xx_xxxy_0,   \
                             ta1_x_xx_xxxy_1,   \
                             ta1_x_xx_xxxyy_0,  \
                             ta1_x_xx_xxxyy_1,  \
                             ta1_x_xx_xxxyz_0,  \
                             ta1_x_xx_xxxyz_1,  \
                             ta1_x_xx_xxxz_0,   \
                             ta1_x_xx_xxxz_1,   \
                             ta1_x_xx_xxxzz_0,  \
                             ta1_x_xx_xxxzz_1,  \
                             ta1_x_xx_xxyy_0,   \
                             ta1_x_xx_xxyy_1,   \
                             ta1_x_xx_xxyyy_0,  \
                             ta1_x_xx_xxyyy_1,  \
                             ta1_x_xx_xxyyz_0,  \
                             ta1_x_xx_xxyyz_1,  \
                             ta1_x_xx_xxyz_0,   \
                             ta1_x_xx_xxyz_1,   \
                             ta1_x_xx_xxyzz_0,  \
                             ta1_x_xx_xxyzz_1,  \
                             ta1_x_xx_xxzz_0,   \
                             ta1_x_xx_xxzz_1,   \
                             ta1_x_xx_xxzzz_0,  \
                             ta1_x_xx_xxzzz_1,  \
                             ta1_x_xx_xyyy_0,   \
                             ta1_x_xx_xyyy_1,   \
                             ta1_x_xx_xyyyy_0,  \
                             ta1_x_xx_xyyyy_1,  \
                             ta1_x_xx_xyyyz_0,  \
                             ta1_x_xx_xyyyz_1,  \
                             ta1_x_xx_xyyz_0,   \
                             ta1_x_xx_xyyz_1,   \
                             ta1_x_xx_xyyzz_0,  \
                             ta1_x_xx_xyyzz_1,  \
                             ta1_x_xx_xyzz_0,   \
                             ta1_x_xx_xyzz_1,   \
                             ta1_x_xx_xyzzz_0,  \
                             ta1_x_xx_xyzzz_1,  \
                             ta1_x_xx_xzzz_0,   \
                             ta1_x_xx_xzzz_1,   \
                             ta1_x_xx_xzzzz_0,  \
                             ta1_x_xx_xzzzz_1,  \
                             ta1_x_xx_yyyy_0,   \
                             ta1_x_xx_yyyy_1,   \
                             ta1_x_xx_yyyyy_0,  \
                             ta1_x_xx_yyyyy_1,  \
                             ta1_x_xx_yyyyz_0,  \
                             ta1_x_xx_yyyyz_1,  \
                             ta1_x_xx_yyyz_0,   \
                             ta1_x_xx_yyyz_1,   \
                             ta1_x_xx_yyyzz_0,  \
                             ta1_x_xx_yyyzz_1,  \
                             ta1_x_xx_yyzz_0,   \
                             ta1_x_xx_yyzz_1,   \
                             ta1_x_xx_yyzzz_0,  \
                             ta1_x_xx_yyzzz_1,  \
                             ta1_x_xx_yzzz_0,   \
                             ta1_x_xx_yzzz_1,   \
                             ta1_x_xx_yzzzz_0,  \
                             ta1_x_xx_yzzzz_1,  \
                             ta1_x_xx_zzzz_0,   \
                             ta1_x_xx_zzzz_1,   \
                             ta1_x_xx_zzzzz_0,  \
                             ta1_x_xx_zzzzz_1,  \
                             ta1_x_xxx_xxxxx_0, \
                             ta1_x_xxx_xxxxy_0, \
                             ta1_x_xxx_xxxxz_0, \
                             ta1_x_xxx_xxxyy_0, \
                             ta1_x_xxx_xxxyz_0, \
                             ta1_x_xxx_xxxzz_0, \
                             ta1_x_xxx_xxyyy_0, \
                             ta1_x_xxx_xxyyz_0, \
                             ta1_x_xxx_xxyzz_0, \
                             ta1_x_xxx_xxzzz_0, \
                             ta1_x_xxx_xyyyy_0, \
                             ta1_x_xxx_xyyyz_0, \
                             ta1_x_xxx_xyyzz_0, \
                             ta1_x_xxx_xyzzz_0, \
                             ta1_x_xxx_xzzzz_0, \
                             ta1_x_xxx_yyyyy_0, \
                             ta1_x_xxx_yyyyz_0, \
                             ta1_x_xxx_yyyzz_0, \
                             ta1_x_xxx_yyzzz_0, \
                             ta1_x_xxx_yzzzz_0, \
                             ta1_x_xxx_zzzzz_0, \
                             ta_xx_xxxxx_1,     \
                             ta_xx_xxxxy_1,     \
                             ta_xx_xxxxz_1,     \
                             ta_xx_xxxyy_1,     \
                             ta_xx_xxxyz_1,     \
                             ta_xx_xxxzz_1,     \
                             ta_xx_xxyyy_1,     \
                             ta_xx_xxyyz_1,     \
                             ta_xx_xxyzz_1,     \
                             ta_xx_xxzzz_1,     \
                             ta_xx_xyyyy_1,     \
                             ta_xx_xyyyz_1,     \
                             ta_xx_xyyzz_1,     \
                             ta_xx_xyzzz_1,     \
                             ta_xx_xzzzz_1,     \
                             ta_xx_yyyyy_1,     \
                             ta_xx_yyyyz_1,     \
                             ta_xx_yyyzz_1,     \
                             ta_xx_yyzzz_1,     \
                             ta_xx_yzzzz_1,     \
                             ta_xx_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxx_xxxxx_0[i] = 2.0 * ta1_x_x_xxxxx_0[i] * fe_0 - 2.0 * ta1_x_x_xxxxx_1[i] * fe_0 + 5.0 * ta1_x_xx_xxxx_0[i] * fe_0 -
                               5.0 * ta1_x_xx_xxxx_1[i] * fe_0 + ta_xx_xxxxx_1[i] + ta1_x_xx_xxxxx_0[i] * pa_x[i] - ta1_x_xx_xxxxx_1[i] * pc_x[i];

        ta1_x_xxx_xxxxy_0[i] = 2.0 * ta1_x_x_xxxxy_0[i] * fe_0 - 2.0 * ta1_x_x_xxxxy_1[i] * fe_0 + 4.0 * ta1_x_xx_xxxy_0[i] * fe_0 -
                               4.0 * ta1_x_xx_xxxy_1[i] * fe_0 + ta_xx_xxxxy_1[i] + ta1_x_xx_xxxxy_0[i] * pa_x[i] - ta1_x_xx_xxxxy_1[i] * pc_x[i];

        ta1_x_xxx_xxxxz_0[i] = 2.0 * ta1_x_x_xxxxz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxxz_1[i] * fe_0 + 4.0 * ta1_x_xx_xxxz_0[i] * fe_0 -
                               4.0 * ta1_x_xx_xxxz_1[i] * fe_0 + ta_xx_xxxxz_1[i] + ta1_x_xx_xxxxz_0[i] * pa_x[i] - ta1_x_xx_xxxxz_1[i] * pc_x[i];

        ta1_x_xxx_xxxyy_0[i] = 2.0 * ta1_x_x_xxxyy_0[i] * fe_0 - 2.0 * ta1_x_x_xxxyy_1[i] * fe_0 + 3.0 * ta1_x_xx_xxyy_0[i] * fe_0 -
                               3.0 * ta1_x_xx_xxyy_1[i] * fe_0 + ta_xx_xxxyy_1[i] + ta1_x_xx_xxxyy_0[i] * pa_x[i] - ta1_x_xx_xxxyy_1[i] * pc_x[i];

        ta1_x_xxx_xxxyz_0[i] = 2.0 * ta1_x_x_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxyz_1[i] * fe_0 + 3.0 * ta1_x_xx_xxyz_0[i] * fe_0 -
                               3.0 * ta1_x_xx_xxyz_1[i] * fe_0 + ta_xx_xxxyz_1[i] + ta1_x_xx_xxxyz_0[i] * pa_x[i] - ta1_x_xx_xxxyz_1[i] * pc_x[i];

        ta1_x_xxx_xxxzz_0[i] = 2.0 * ta1_x_x_xxxzz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxzz_1[i] * fe_0 + 3.0 * ta1_x_xx_xxzz_0[i] * fe_0 -
                               3.0 * ta1_x_xx_xxzz_1[i] * fe_0 + ta_xx_xxxzz_1[i] + ta1_x_xx_xxxzz_0[i] * pa_x[i] - ta1_x_xx_xxxzz_1[i] * pc_x[i];

        ta1_x_xxx_xxyyy_0[i] = 2.0 * ta1_x_x_xxyyy_0[i] * fe_0 - 2.0 * ta1_x_x_xxyyy_1[i] * fe_0 + 2.0 * ta1_x_xx_xyyy_0[i] * fe_0 -
                               2.0 * ta1_x_xx_xyyy_1[i] * fe_0 + ta_xx_xxyyy_1[i] + ta1_x_xx_xxyyy_0[i] * pa_x[i] - ta1_x_xx_xxyyy_1[i] * pc_x[i];

        ta1_x_xxx_xxyyz_0[i] = 2.0 * ta1_x_x_xxyyz_0[i] * fe_0 - 2.0 * ta1_x_x_xxyyz_1[i] * fe_0 + 2.0 * ta1_x_xx_xyyz_0[i] * fe_0 -
                               2.0 * ta1_x_xx_xyyz_1[i] * fe_0 + ta_xx_xxyyz_1[i] + ta1_x_xx_xxyyz_0[i] * pa_x[i] - ta1_x_xx_xxyyz_1[i] * pc_x[i];

        ta1_x_xxx_xxyzz_0[i] = 2.0 * ta1_x_x_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_x_xxyzz_1[i] * fe_0 + 2.0 * ta1_x_xx_xyzz_0[i] * fe_0 -
                               2.0 * ta1_x_xx_xyzz_1[i] * fe_0 + ta_xx_xxyzz_1[i] + ta1_x_xx_xxyzz_0[i] * pa_x[i] - ta1_x_xx_xxyzz_1[i] * pc_x[i];

        ta1_x_xxx_xxzzz_0[i] = 2.0 * ta1_x_x_xxzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xxzzz_1[i] * fe_0 + 2.0 * ta1_x_xx_xzzz_0[i] * fe_0 -
                               2.0 * ta1_x_xx_xzzz_1[i] * fe_0 + ta_xx_xxzzz_1[i] + ta1_x_xx_xxzzz_0[i] * pa_x[i] - ta1_x_xx_xxzzz_1[i] * pc_x[i];

        ta1_x_xxx_xyyyy_0[i] = 2.0 * ta1_x_x_xyyyy_0[i] * fe_0 - 2.0 * ta1_x_x_xyyyy_1[i] * fe_0 + ta1_x_xx_yyyy_0[i] * fe_0 -
                               ta1_x_xx_yyyy_1[i] * fe_0 + ta_xx_xyyyy_1[i] + ta1_x_xx_xyyyy_0[i] * pa_x[i] - ta1_x_xx_xyyyy_1[i] * pc_x[i];

        ta1_x_xxx_xyyyz_0[i] = 2.0 * ta1_x_x_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_x_xyyyz_1[i] * fe_0 + ta1_x_xx_yyyz_0[i] * fe_0 -
                               ta1_x_xx_yyyz_1[i] * fe_0 + ta_xx_xyyyz_1[i] + ta1_x_xx_xyyyz_0[i] * pa_x[i] - ta1_x_xx_xyyyz_1[i] * pc_x[i];

        ta1_x_xxx_xyyzz_0[i] = 2.0 * ta1_x_x_xyyzz_0[i] * fe_0 - 2.0 * ta1_x_x_xyyzz_1[i] * fe_0 + ta1_x_xx_yyzz_0[i] * fe_0 -
                               ta1_x_xx_yyzz_1[i] * fe_0 + ta_xx_xyyzz_1[i] + ta1_x_xx_xyyzz_0[i] * pa_x[i] - ta1_x_xx_xyyzz_1[i] * pc_x[i];

        ta1_x_xxx_xyzzz_0[i] = 2.0 * ta1_x_x_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xyzzz_1[i] * fe_0 + ta1_x_xx_yzzz_0[i] * fe_0 -
                               ta1_x_xx_yzzz_1[i] * fe_0 + ta_xx_xyzzz_1[i] + ta1_x_xx_xyzzz_0[i] * pa_x[i] - ta1_x_xx_xyzzz_1[i] * pc_x[i];

        ta1_x_xxx_xzzzz_0[i] = 2.0 * ta1_x_x_xzzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xzzzz_1[i] * fe_0 + ta1_x_xx_zzzz_0[i] * fe_0 -
                               ta1_x_xx_zzzz_1[i] * fe_0 + ta_xx_xzzzz_1[i] + ta1_x_xx_xzzzz_0[i] * pa_x[i] - ta1_x_xx_xzzzz_1[i] * pc_x[i];

        ta1_x_xxx_yyyyy_0[i] = 2.0 * ta1_x_x_yyyyy_0[i] * fe_0 - 2.0 * ta1_x_x_yyyyy_1[i] * fe_0 + ta_xx_yyyyy_1[i] + ta1_x_xx_yyyyy_0[i] * pa_x[i] -
                               ta1_x_xx_yyyyy_1[i] * pc_x[i];

        ta1_x_xxx_yyyyz_0[i] = 2.0 * ta1_x_x_yyyyz_0[i] * fe_0 - 2.0 * ta1_x_x_yyyyz_1[i] * fe_0 + ta_xx_yyyyz_1[i] + ta1_x_xx_yyyyz_0[i] * pa_x[i] -
                               ta1_x_xx_yyyyz_1[i] * pc_x[i];

        ta1_x_xxx_yyyzz_0[i] = 2.0 * ta1_x_x_yyyzz_0[i] * fe_0 - 2.0 * ta1_x_x_yyyzz_1[i] * fe_0 + ta_xx_yyyzz_1[i] + ta1_x_xx_yyyzz_0[i] * pa_x[i] -
                               ta1_x_xx_yyyzz_1[i] * pc_x[i];

        ta1_x_xxx_yyzzz_0[i] = 2.0 * ta1_x_x_yyzzz_0[i] * fe_0 - 2.0 * ta1_x_x_yyzzz_1[i] * fe_0 + ta_xx_yyzzz_1[i] + ta1_x_xx_yyzzz_0[i] * pa_x[i] -
                               ta1_x_xx_yyzzz_1[i] * pc_x[i];

        ta1_x_xxx_yzzzz_0[i] = 2.0 * ta1_x_x_yzzzz_0[i] * fe_0 - 2.0 * ta1_x_x_yzzzz_1[i] * fe_0 + ta_xx_yzzzz_1[i] + ta1_x_xx_yzzzz_0[i] * pa_x[i] -
                               ta1_x_xx_yzzzz_1[i] * pc_x[i];

        ta1_x_xxx_zzzzz_0[i] = 2.0 * ta1_x_x_zzzzz_0[i] * fe_0 - 2.0 * ta1_x_x_zzzzz_1[i] * fe_0 + ta_xx_zzzzz_1[i] + ta1_x_xx_zzzzz_0[i] * pa_x[i] -
                               ta1_x_xx_zzzzz_1[i] * pc_x[i];
    }

    // Set up 21-42 components of targeted buffer : FH

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

    auto ta1_x_xxy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 37);

    auto ta1_x_xxy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 38);

    auto ta1_x_xxy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 39);

    auto ta1_x_xxy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 40);

    auto ta1_x_xxy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 41);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_x_xx_xxxx_0,   \
                             ta1_x_xx_xxxx_1,   \
                             ta1_x_xx_xxxxx_0,  \
                             ta1_x_xx_xxxxx_1,  \
                             ta1_x_xx_xxxxy_0,  \
                             ta1_x_xx_xxxxy_1,  \
                             ta1_x_xx_xxxxz_0,  \
                             ta1_x_xx_xxxxz_1,  \
                             ta1_x_xx_xxxy_0,   \
                             ta1_x_xx_xxxy_1,   \
                             ta1_x_xx_xxxyy_0,  \
                             ta1_x_xx_xxxyy_1,  \
                             ta1_x_xx_xxxyz_0,  \
                             ta1_x_xx_xxxyz_1,  \
                             ta1_x_xx_xxxz_0,   \
                             ta1_x_xx_xxxz_1,   \
                             ta1_x_xx_xxxzz_0,  \
                             ta1_x_xx_xxxzz_1,  \
                             ta1_x_xx_xxyy_0,   \
                             ta1_x_xx_xxyy_1,   \
                             ta1_x_xx_xxyyy_0,  \
                             ta1_x_xx_xxyyy_1,  \
                             ta1_x_xx_xxyyz_0,  \
                             ta1_x_xx_xxyyz_1,  \
                             ta1_x_xx_xxyz_0,   \
                             ta1_x_xx_xxyz_1,   \
                             ta1_x_xx_xxyzz_0,  \
                             ta1_x_xx_xxyzz_1,  \
                             ta1_x_xx_xxzz_0,   \
                             ta1_x_xx_xxzz_1,   \
                             ta1_x_xx_xxzzz_0,  \
                             ta1_x_xx_xxzzz_1,  \
                             ta1_x_xx_xyyy_0,   \
                             ta1_x_xx_xyyy_1,   \
                             ta1_x_xx_xyyyy_0,  \
                             ta1_x_xx_xyyyy_1,  \
                             ta1_x_xx_xyyyz_0,  \
                             ta1_x_xx_xyyyz_1,  \
                             ta1_x_xx_xyyz_0,   \
                             ta1_x_xx_xyyz_1,   \
                             ta1_x_xx_xyyzz_0,  \
                             ta1_x_xx_xyyzz_1,  \
                             ta1_x_xx_xyzz_0,   \
                             ta1_x_xx_xyzz_1,   \
                             ta1_x_xx_xyzzz_0,  \
                             ta1_x_xx_xyzzz_1,  \
                             ta1_x_xx_xzzz_0,   \
                             ta1_x_xx_xzzz_1,   \
                             ta1_x_xx_xzzzz_0,  \
                             ta1_x_xx_xzzzz_1,  \
                             ta1_x_xx_yyyy_0,   \
                             ta1_x_xx_yyyy_1,   \
                             ta1_x_xx_yyyyy_0,  \
                             ta1_x_xx_yyyyy_1,  \
                             ta1_x_xx_yyyyz_0,  \
                             ta1_x_xx_yyyyz_1,  \
                             ta1_x_xx_yyyz_0,   \
                             ta1_x_xx_yyyz_1,   \
                             ta1_x_xx_yyyzz_0,  \
                             ta1_x_xx_yyyzz_1,  \
                             ta1_x_xx_yyzz_0,   \
                             ta1_x_xx_yyzz_1,   \
                             ta1_x_xx_yyzzz_0,  \
                             ta1_x_xx_yyzzz_1,  \
                             ta1_x_xx_yzzz_0,   \
                             ta1_x_xx_yzzz_1,   \
                             ta1_x_xx_yzzzz_0,  \
                             ta1_x_xx_yzzzz_1,  \
                             ta1_x_xx_zzzz_0,   \
                             ta1_x_xx_zzzz_1,   \
                             ta1_x_xx_zzzzz_0,  \
                             ta1_x_xx_zzzzz_1,  \
                             ta1_x_xxy_xxxxx_0, \
                             ta1_x_xxy_xxxxy_0, \
                             ta1_x_xxy_xxxxz_0, \
                             ta1_x_xxy_xxxyy_0, \
                             ta1_x_xxy_xxxyz_0, \
                             ta1_x_xxy_xxxzz_0, \
                             ta1_x_xxy_xxyyy_0, \
                             ta1_x_xxy_xxyyz_0, \
                             ta1_x_xxy_xxyzz_0, \
                             ta1_x_xxy_xxzzz_0, \
                             ta1_x_xxy_xyyyy_0, \
                             ta1_x_xxy_xyyyz_0, \
                             ta1_x_xxy_xyyzz_0, \
                             ta1_x_xxy_xyzzz_0, \
                             ta1_x_xxy_xzzzz_0, \
                             ta1_x_xxy_yyyyy_0, \
                             ta1_x_xxy_yyyyz_0, \
                             ta1_x_xxy_yyyzz_0, \
                             ta1_x_xxy_yyzzz_0, \
                             ta1_x_xxy_yzzzz_0, \
                             ta1_x_xxy_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxy_xxxxx_0[i] = ta1_x_xx_xxxxx_0[i] * pa_y[i] - ta1_x_xx_xxxxx_1[i] * pc_y[i];

        ta1_x_xxy_xxxxy_0[i] = ta1_x_xx_xxxx_0[i] * fe_0 - ta1_x_xx_xxxx_1[i] * fe_0 + ta1_x_xx_xxxxy_0[i] * pa_y[i] - ta1_x_xx_xxxxy_1[i] * pc_y[i];

        ta1_x_xxy_xxxxz_0[i] = ta1_x_xx_xxxxz_0[i] * pa_y[i] - ta1_x_xx_xxxxz_1[i] * pc_y[i];

        ta1_x_xxy_xxxyy_0[i] =
            2.0 * ta1_x_xx_xxxy_0[i] * fe_0 - 2.0 * ta1_x_xx_xxxy_1[i] * fe_0 + ta1_x_xx_xxxyy_0[i] * pa_y[i] - ta1_x_xx_xxxyy_1[i] * pc_y[i];

        ta1_x_xxy_xxxyz_0[i] = ta1_x_xx_xxxz_0[i] * fe_0 - ta1_x_xx_xxxz_1[i] * fe_0 + ta1_x_xx_xxxyz_0[i] * pa_y[i] - ta1_x_xx_xxxyz_1[i] * pc_y[i];

        ta1_x_xxy_xxxzz_0[i] = ta1_x_xx_xxxzz_0[i] * pa_y[i] - ta1_x_xx_xxxzz_1[i] * pc_y[i];

        ta1_x_xxy_xxyyy_0[i] =
            3.0 * ta1_x_xx_xxyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxyy_1[i] * fe_0 + ta1_x_xx_xxyyy_0[i] * pa_y[i] - ta1_x_xx_xxyyy_1[i] * pc_y[i];

        ta1_x_xxy_xxyyz_0[i] =
            2.0 * ta1_x_xx_xxyz_0[i] * fe_0 - 2.0 * ta1_x_xx_xxyz_1[i] * fe_0 + ta1_x_xx_xxyyz_0[i] * pa_y[i] - ta1_x_xx_xxyyz_1[i] * pc_y[i];

        ta1_x_xxy_xxyzz_0[i] = ta1_x_xx_xxzz_0[i] * fe_0 - ta1_x_xx_xxzz_1[i] * fe_0 + ta1_x_xx_xxyzz_0[i] * pa_y[i] - ta1_x_xx_xxyzz_1[i] * pc_y[i];

        ta1_x_xxy_xxzzz_0[i] = ta1_x_xx_xxzzz_0[i] * pa_y[i] - ta1_x_xx_xxzzz_1[i] * pc_y[i];

        ta1_x_xxy_xyyyy_0[i] =
            4.0 * ta1_x_xx_xyyy_0[i] * fe_0 - 4.0 * ta1_x_xx_xyyy_1[i] * fe_0 + ta1_x_xx_xyyyy_0[i] * pa_y[i] - ta1_x_xx_xyyyy_1[i] * pc_y[i];

        ta1_x_xxy_xyyyz_0[i] =
            3.0 * ta1_x_xx_xyyz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyyz_1[i] * fe_0 + ta1_x_xx_xyyyz_0[i] * pa_y[i] - ta1_x_xx_xyyyz_1[i] * pc_y[i];

        ta1_x_xxy_xyyzz_0[i] =
            2.0 * ta1_x_xx_xyzz_0[i] * fe_0 - 2.0 * ta1_x_xx_xyzz_1[i] * fe_0 + ta1_x_xx_xyyzz_0[i] * pa_y[i] - ta1_x_xx_xyyzz_1[i] * pc_y[i];

        ta1_x_xxy_xyzzz_0[i] = ta1_x_xx_xzzz_0[i] * fe_0 - ta1_x_xx_xzzz_1[i] * fe_0 + ta1_x_xx_xyzzz_0[i] * pa_y[i] - ta1_x_xx_xyzzz_1[i] * pc_y[i];

        ta1_x_xxy_xzzzz_0[i] = ta1_x_xx_xzzzz_0[i] * pa_y[i] - ta1_x_xx_xzzzz_1[i] * pc_y[i];

        ta1_x_xxy_yyyyy_0[i] =
            5.0 * ta1_x_xx_yyyy_0[i] * fe_0 - 5.0 * ta1_x_xx_yyyy_1[i] * fe_0 + ta1_x_xx_yyyyy_0[i] * pa_y[i] - ta1_x_xx_yyyyy_1[i] * pc_y[i];

        ta1_x_xxy_yyyyz_0[i] =
            4.0 * ta1_x_xx_yyyz_0[i] * fe_0 - 4.0 * ta1_x_xx_yyyz_1[i] * fe_0 + ta1_x_xx_yyyyz_0[i] * pa_y[i] - ta1_x_xx_yyyyz_1[i] * pc_y[i];

        ta1_x_xxy_yyyzz_0[i] =
            3.0 * ta1_x_xx_yyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyzz_1[i] * fe_0 + ta1_x_xx_yyyzz_0[i] * pa_y[i] - ta1_x_xx_yyyzz_1[i] * pc_y[i];

        ta1_x_xxy_yyzzz_0[i] =
            2.0 * ta1_x_xx_yzzz_0[i] * fe_0 - 2.0 * ta1_x_xx_yzzz_1[i] * fe_0 + ta1_x_xx_yyzzz_0[i] * pa_y[i] - ta1_x_xx_yyzzz_1[i] * pc_y[i];

        ta1_x_xxy_yzzzz_0[i] = ta1_x_xx_zzzz_0[i] * fe_0 - ta1_x_xx_zzzz_1[i] * fe_0 + ta1_x_xx_yzzzz_0[i] * pa_y[i] - ta1_x_xx_yzzzz_1[i] * pc_y[i];

        ta1_x_xxy_zzzzz_0[i] = ta1_x_xx_zzzzz_0[i] * pa_y[i] - ta1_x_xx_zzzzz_1[i] * pc_y[i];
    }

    // Set up 42-63 components of targeted buffer : FH

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

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_x_xx_xxxx_0,   \
                             ta1_x_xx_xxxx_1,   \
                             ta1_x_xx_xxxxx_0,  \
                             ta1_x_xx_xxxxx_1,  \
                             ta1_x_xx_xxxxy_0,  \
                             ta1_x_xx_xxxxy_1,  \
                             ta1_x_xx_xxxxz_0,  \
                             ta1_x_xx_xxxxz_1,  \
                             ta1_x_xx_xxxy_0,   \
                             ta1_x_xx_xxxy_1,   \
                             ta1_x_xx_xxxyy_0,  \
                             ta1_x_xx_xxxyy_1,  \
                             ta1_x_xx_xxxyz_0,  \
                             ta1_x_xx_xxxyz_1,  \
                             ta1_x_xx_xxxz_0,   \
                             ta1_x_xx_xxxz_1,   \
                             ta1_x_xx_xxxzz_0,  \
                             ta1_x_xx_xxxzz_1,  \
                             ta1_x_xx_xxyy_0,   \
                             ta1_x_xx_xxyy_1,   \
                             ta1_x_xx_xxyyy_0,  \
                             ta1_x_xx_xxyyy_1,  \
                             ta1_x_xx_xxyyz_0,  \
                             ta1_x_xx_xxyyz_1,  \
                             ta1_x_xx_xxyz_0,   \
                             ta1_x_xx_xxyz_1,   \
                             ta1_x_xx_xxyzz_0,  \
                             ta1_x_xx_xxyzz_1,  \
                             ta1_x_xx_xxzz_0,   \
                             ta1_x_xx_xxzz_1,   \
                             ta1_x_xx_xxzzz_0,  \
                             ta1_x_xx_xxzzz_1,  \
                             ta1_x_xx_xyyy_0,   \
                             ta1_x_xx_xyyy_1,   \
                             ta1_x_xx_xyyyy_0,  \
                             ta1_x_xx_xyyyy_1,  \
                             ta1_x_xx_xyyyz_0,  \
                             ta1_x_xx_xyyyz_1,  \
                             ta1_x_xx_xyyz_0,   \
                             ta1_x_xx_xyyz_1,   \
                             ta1_x_xx_xyyzz_0,  \
                             ta1_x_xx_xyyzz_1,  \
                             ta1_x_xx_xyzz_0,   \
                             ta1_x_xx_xyzz_1,   \
                             ta1_x_xx_xyzzz_0,  \
                             ta1_x_xx_xyzzz_1,  \
                             ta1_x_xx_xzzz_0,   \
                             ta1_x_xx_xzzz_1,   \
                             ta1_x_xx_xzzzz_0,  \
                             ta1_x_xx_xzzzz_1,  \
                             ta1_x_xx_yyyy_0,   \
                             ta1_x_xx_yyyy_1,   \
                             ta1_x_xx_yyyyy_0,  \
                             ta1_x_xx_yyyyy_1,  \
                             ta1_x_xx_yyyyz_0,  \
                             ta1_x_xx_yyyyz_1,  \
                             ta1_x_xx_yyyz_0,   \
                             ta1_x_xx_yyyz_1,   \
                             ta1_x_xx_yyyzz_0,  \
                             ta1_x_xx_yyyzz_1,  \
                             ta1_x_xx_yyzz_0,   \
                             ta1_x_xx_yyzz_1,   \
                             ta1_x_xx_yyzzz_0,  \
                             ta1_x_xx_yyzzz_1,  \
                             ta1_x_xx_yzzz_0,   \
                             ta1_x_xx_yzzz_1,   \
                             ta1_x_xx_yzzzz_0,  \
                             ta1_x_xx_yzzzz_1,  \
                             ta1_x_xx_zzzz_0,   \
                             ta1_x_xx_zzzz_1,   \
                             ta1_x_xx_zzzzz_0,  \
                             ta1_x_xx_zzzzz_1,  \
                             ta1_x_xxz_xxxxx_0, \
                             ta1_x_xxz_xxxxy_0, \
                             ta1_x_xxz_xxxxz_0, \
                             ta1_x_xxz_xxxyy_0, \
                             ta1_x_xxz_xxxyz_0, \
                             ta1_x_xxz_xxxzz_0, \
                             ta1_x_xxz_xxyyy_0, \
                             ta1_x_xxz_xxyyz_0, \
                             ta1_x_xxz_xxyzz_0, \
                             ta1_x_xxz_xxzzz_0, \
                             ta1_x_xxz_xyyyy_0, \
                             ta1_x_xxz_xyyyz_0, \
                             ta1_x_xxz_xyyzz_0, \
                             ta1_x_xxz_xyzzz_0, \
                             ta1_x_xxz_xzzzz_0, \
                             ta1_x_xxz_yyyyy_0, \
                             ta1_x_xxz_yyyyz_0, \
                             ta1_x_xxz_yyyzz_0, \
                             ta1_x_xxz_yyzzz_0, \
                             ta1_x_xxz_yzzzz_0, \
                             ta1_x_xxz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxz_xxxxx_0[i] = ta1_x_xx_xxxxx_0[i] * pa_z[i] - ta1_x_xx_xxxxx_1[i] * pc_z[i];

        ta1_x_xxz_xxxxy_0[i] = ta1_x_xx_xxxxy_0[i] * pa_z[i] - ta1_x_xx_xxxxy_1[i] * pc_z[i];

        ta1_x_xxz_xxxxz_0[i] = ta1_x_xx_xxxx_0[i] * fe_0 - ta1_x_xx_xxxx_1[i] * fe_0 + ta1_x_xx_xxxxz_0[i] * pa_z[i] - ta1_x_xx_xxxxz_1[i] * pc_z[i];

        ta1_x_xxz_xxxyy_0[i] = ta1_x_xx_xxxyy_0[i] * pa_z[i] - ta1_x_xx_xxxyy_1[i] * pc_z[i];

        ta1_x_xxz_xxxyz_0[i] = ta1_x_xx_xxxy_0[i] * fe_0 - ta1_x_xx_xxxy_1[i] * fe_0 + ta1_x_xx_xxxyz_0[i] * pa_z[i] - ta1_x_xx_xxxyz_1[i] * pc_z[i];

        ta1_x_xxz_xxxzz_0[i] =
            2.0 * ta1_x_xx_xxxz_0[i] * fe_0 - 2.0 * ta1_x_xx_xxxz_1[i] * fe_0 + ta1_x_xx_xxxzz_0[i] * pa_z[i] - ta1_x_xx_xxxzz_1[i] * pc_z[i];

        ta1_x_xxz_xxyyy_0[i] = ta1_x_xx_xxyyy_0[i] * pa_z[i] - ta1_x_xx_xxyyy_1[i] * pc_z[i];

        ta1_x_xxz_xxyyz_0[i] = ta1_x_xx_xxyy_0[i] * fe_0 - ta1_x_xx_xxyy_1[i] * fe_0 + ta1_x_xx_xxyyz_0[i] * pa_z[i] - ta1_x_xx_xxyyz_1[i] * pc_z[i];

        ta1_x_xxz_xxyzz_0[i] =
            2.0 * ta1_x_xx_xxyz_0[i] * fe_0 - 2.0 * ta1_x_xx_xxyz_1[i] * fe_0 + ta1_x_xx_xxyzz_0[i] * pa_z[i] - ta1_x_xx_xxyzz_1[i] * pc_z[i];

        ta1_x_xxz_xxzzz_0[i] =
            3.0 * ta1_x_xx_xxzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxzz_1[i] * fe_0 + ta1_x_xx_xxzzz_0[i] * pa_z[i] - ta1_x_xx_xxzzz_1[i] * pc_z[i];

        ta1_x_xxz_xyyyy_0[i] = ta1_x_xx_xyyyy_0[i] * pa_z[i] - ta1_x_xx_xyyyy_1[i] * pc_z[i];

        ta1_x_xxz_xyyyz_0[i] = ta1_x_xx_xyyy_0[i] * fe_0 - ta1_x_xx_xyyy_1[i] * fe_0 + ta1_x_xx_xyyyz_0[i] * pa_z[i] - ta1_x_xx_xyyyz_1[i] * pc_z[i];

        ta1_x_xxz_xyyzz_0[i] =
            2.0 * ta1_x_xx_xyyz_0[i] * fe_0 - 2.0 * ta1_x_xx_xyyz_1[i] * fe_0 + ta1_x_xx_xyyzz_0[i] * pa_z[i] - ta1_x_xx_xyyzz_1[i] * pc_z[i];

        ta1_x_xxz_xyzzz_0[i] =
            3.0 * ta1_x_xx_xyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyzz_1[i] * fe_0 + ta1_x_xx_xyzzz_0[i] * pa_z[i] - ta1_x_xx_xyzzz_1[i] * pc_z[i];

        ta1_x_xxz_xzzzz_0[i] =
            4.0 * ta1_x_xx_xzzz_0[i] * fe_0 - 4.0 * ta1_x_xx_xzzz_1[i] * fe_0 + ta1_x_xx_xzzzz_0[i] * pa_z[i] - ta1_x_xx_xzzzz_1[i] * pc_z[i];

        ta1_x_xxz_yyyyy_0[i] = ta1_x_xx_yyyyy_0[i] * pa_z[i] - ta1_x_xx_yyyyy_1[i] * pc_z[i];

        ta1_x_xxz_yyyyz_0[i] = ta1_x_xx_yyyy_0[i] * fe_0 - ta1_x_xx_yyyy_1[i] * fe_0 + ta1_x_xx_yyyyz_0[i] * pa_z[i] - ta1_x_xx_yyyyz_1[i] * pc_z[i];

        ta1_x_xxz_yyyzz_0[i] =
            2.0 * ta1_x_xx_yyyz_0[i] * fe_0 - 2.0 * ta1_x_xx_yyyz_1[i] * fe_0 + ta1_x_xx_yyyzz_0[i] * pa_z[i] - ta1_x_xx_yyyzz_1[i] * pc_z[i];

        ta1_x_xxz_yyzzz_0[i] =
            3.0 * ta1_x_xx_yyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyzz_1[i] * fe_0 + ta1_x_xx_yyzzz_0[i] * pa_z[i] - ta1_x_xx_yyzzz_1[i] * pc_z[i];

        ta1_x_xxz_yzzzz_0[i] =
            4.0 * ta1_x_xx_yzzz_0[i] * fe_0 - 4.0 * ta1_x_xx_yzzz_1[i] * fe_0 + ta1_x_xx_yzzzz_0[i] * pa_z[i] - ta1_x_xx_yzzzz_1[i] * pc_z[i];

        ta1_x_xxz_zzzzz_0[i] =
            5.0 * ta1_x_xx_zzzz_0[i] * fe_0 - 5.0 * ta1_x_xx_zzzz_1[i] * fe_0 + ta1_x_xx_zzzzz_0[i] * pa_z[i] - ta1_x_xx_zzzzz_1[i] * pc_z[i];
    }

    // Set up 63-84 components of targeted buffer : FH

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

    auto ta1_x_xyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 83);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_x_x_xxxxx_0,   \
                             ta1_x_x_xxxxx_1,   \
                             ta1_x_x_xxxxz_0,   \
                             ta1_x_x_xxxxz_1,   \
                             ta1_x_x_xxxzz_0,   \
                             ta1_x_x_xxxzz_1,   \
                             ta1_x_x_xxzzz_0,   \
                             ta1_x_x_xxzzz_1,   \
                             ta1_x_x_xzzzz_0,   \
                             ta1_x_x_xzzzz_1,   \
                             ta1_x_xy_xxxxx_0,  \
                             ta1_x_xy_xxxxx_1,  \
                             ta1_x_xy_xxxxz_0,  \
                             ta1_x_xy_xxxxz_1,  \
                             ta1_x_xy_xxxzz_0,  \
                             ta1_x_xy_xxxzz_1,  \
                             ta1_x_xy_xxzzz_0,  \
                             ta1_x_xy_xxzzz_1,  \
                             ta1_x_xy_xzzzz_0,  \
                             ta1_x_xy_xzzzz_1,  \
                             ta1_x_xyy_xxxxx_0, \
                             ta1_x_xyy_xxxxy_0, \
                             ta1_x_xyy_xxxxz_0, \
                             ta1_x_xyy_xxxyy_0, \
                             ta1_x_xyy_xxxyz_0, \
                             ta1_x_xyy_xxxzz_0, \
                             ta1_x_xyy_xxyyy_0, \
                             ta1_x_xyy_xxyyz_0, \
                             ta1_x_xyy_xxyzz_0, \
                             ta1_x_xyy_xxzzz_0, \
                             ta1_x_xyy_xyyyy_0, \
                             ta1_x_xyy_xyyyz_0, \
                             ta1_x_xyy_xyyzz_0, \
                             ta1_x_xyy_xyzzz_0, \
                             ta1_x_xyy_xzzzz_0, \
                             ta1_x_xyy_yyyyy_0, \
                             ta1_x_xyy_yyyyz_0, \
                             ta1_x_xyy_yyyzz_0, \
                             ta1_x_xyy_yyzzz_0, \
                             ta1_x_xyy_yzzzz_0, \
                             ta1_x_xyy_zzzzz_0, \
                             ta1_x_yy_xxxxy_0,  \
                             ta1_x_yy_xxxxy_1,  \
                             ta1_x_yy_xxxy_0,   \
                             ta1_x_yy_xxxy_1,   \
                             ta1_x_yy_xxxyy_0,  \
                             ta1_x_yy_xxxyy_1,  \
                             ta1_x_yy_xxxyz_0,  \
                             ta1_x_yy_xxxyz_1,  \
                             ta1_x_yy_xxyy_0,   \
                             ta1_x_yy_xxyy_1,   \
                             ta1_x_yy_xxyyy_0,  \
                             ta1_x_yy_xxyyy_1,  \
                             ta1_x_yy_xxyyz_0,  \
                             ta1_x_yy_xxyyz_1,  \
                             ta1_x_yy_xxyz_0,   \
                             ta1_x_yy_xxyz_1,   \
                             ta1_x_yy_xxyzz_0,  \
                             ta1_x_yy_xxyzz_1,  \
                             ta1_x_yy_xyyy_0,   \
                             ta1_x_yy_xyyy_1,   \
                             ta1_x_yy_xyyyy_0,  \
                             ta1_x_yy_xyyyy_1,  \
                             ta1_x_yy_xyyyz_0,  \
                             ta1_x_yy_xyyyz_1,  \
                             ta1_x_yy_xyyz_0,   \
                             ta1_x_yy_xyyz_1,   \
                             ta1_x_yy_xyyzz_0,  \
                             ta1_x_yy_xyyzz_1,  \
                             ta1_x_yy_xyzz_0,   \
                             ta1_x_yy_xyzz_1,   \
                             ta1_x_yy_xyzzz_0,  \
                             ta1_x_yy_xyzzz_1,  \
                             ta1_x_yy_yyyy_0,   \
                             ta1_x_yy_yyyy_1,   \
                             ta1_x_yy_yyyyy_0,  \
                             ta1_x_yy_yyyyy_1,  \
                             ta1_x_yy_yyyyz_0,  \
                             ta1_x_yy_yyyyz_1,  \
                             ta1_x_yy_yyyz_0,   \
                             ta1_x_yy_yyyz_1,   \
                             ta1_x_yy_yyyzz_0,  \
                             ta1_x_yy_yyyzz_1,  \
                             ta1_x_yy_yyzz_0,   \
                             ta1_x_yy_yyzz_1,   \
                             ta1_x_yy_yyzzz_0,  \
                             ta1_x_yy_yyzzz_1,  \
                             ta1_x_yy_yzzz_0,   \
                             ta1_x_yy_yzzz_1,   \
                             ta1_x_yy_yzzzz_0,  \
                             ta1_x_yy_yzzzz_1,  \
                             ta1_x_yy_zzzzz_0,  \
                             ta1_x_yy_zzzzz_1,  \
                             ta_yy_xxxxy_1,     \
                             ta_yy_xxxyy_1,     \
                             ta_yy_xxxyz_1,     \
                             ta_yy_xxyyy_1,     \
                             ta_yy_xxyyz_1,     \
                             ta_yy_xxyzz_1,     \
                             ta_yy_xyyyy_1,     \
                             ta_yy_xyyyz_1,     \
                             ta_yy_xyyzz_1,     \
                             ta_yy_xyzzz_1,     \
                             ta_yy_yyyyy_1,     \
                             ta_yy_yyyyz_1,     \
                             ta_yy_yyyzz_1,     \
                             ta_yy_yyzzz_1,     \
                             ta_yy_yzzzz_1,     \
                             ta_yy_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyy_xxxxx_0[i] = ta1_x_x_xxxxx_0[i] * fe_0 - ta1_x_x_xxxxx_1[i] * fe_0 + ta1_x_xy_xxxxx_0[i] * pa_y[i] - ta1_x_xy_xxxxx_1[i] * pc_y[i];

        ta1_x_xyy_xxxxy_0[i] = 4.0 * ta1_x_yy_xxxy_0[i] * fe_0 - 4.0 * ta1_x_yy_xxxy_1[i] * fe_0 + ta_yy_xxxxy_1[i] + ta1_x_yy_xxxxy_0[i] * pa_x[i] -
                               ta1_x_yy_xxxxy_1[i] * pc_x[i];

        ta1_x_xyy_xxxxz_0[i] = ta1_x_x_xxxxz_0[i] * fe_0 - ta1_x_x_xxxxz_1[i] * fe_0 + ta1_x_xy_xxxxz_0[i] * pa_y[i] - ta1_x_xy_xxxxz_1[i] * pc_y[i];

        ta1_x_xyy_xxxyy_0[i] = 3.0 * ta1_x_yy_xxyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyy_1[i] * fe_0 + ta_yy_xxxyy_1[i] + ta1_x_yy_xxxyy_0[i] * pa_x[i] -
                               ta1_x_yy_xxxyy_1[i] * pc_x[i];

        ta1_x_xyy_xxxyz_0[i] = 3.0 * ta1_x_yy_xxyz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyz_1[i] * fe_0 + ta_yy_xxxyz_1[i] + ta1_x_yy_xxxyz_0[i] * pa_x[i] -
                               ta1_x_yy_xxxyz_1[i] * pc_x[i];

        ta1_x_xyy_xxxzz_0[i] = ta1_x_x_xxxzz_0[i] * fe_0 - ta1_x_x_xxxzz_1[i] * fe_0 + ta1_x_xy_xxxzz_0[i] * pa_y[i] - ta1_x_xy_xxxzz_1[i] * pc_y[i];

        ta1_x_xyy_xxyyy_0[i] = 2.0 * ta1_x_yy_xyyy_0[i] * fe_0 - 2.0 * ta1_x_yy_xyyy_1[i] * fe_0 + ta_yy_xxyyy_1[i] + ta1_x_yy_xxyyy_0[i] * pa_x[i] -
                               ta1_x_yy_xxyyy_1[i] * pc_x[i];

        ta1_x_xyy_xxyyz_0[i] = 2.0 * ta1_x_yy_xyyz_0[i] * fe_0 - 2.0 * ta1_x_yy_xyyz_1[i] * fe_0 + ta_yy_xxyyz_1[i] + ta1_x_yy_xxyyz_0[i] * pa_x[i] -
                               ta1_x_yy_xxyyz_1[i] * pc_x[i];

        ta1_x_xyy_xxyzz_0[i] = 2.0 * ta1_x_yy_xyzz_0[i] * fe_0 - 2.0 * ta1_x_yy_xyzz_1[i] * fe_0 + ta_yy_xxyzz_1[i] + ta1_x_yy_xxyzz_0[i] * pa_x[i] -
                               ta1_x_yy_xxyzz_1[i] * pc_x[i];

        ta1_x_xyy_xxzzz_0[i] = ta1_x_x_xxzzz_0[i] * fe_0 - ta1_x_x_xxzzz_1[i] * fe_0 + ta1_x_xy_xxzzz_0[i] * pa_y[i] - ta1_x_xy_xxzzz_1[i] * pc_y[i];

        ta1_x_xyy_xyyyy_0[i] =
            ta1_x_yy_yyyy_0[i] * fe_0 - ta1_x_yy_yyyy_1[i] * fe_0 + ta_yy_xyyyy_1[i] + ta1_x_yy_xyyyy_0[i] * pa_x[i] - ta1_x_yy_xyyyy_1[i] * pc_x[i];

        ta1_x_xyy_xyyyz_0[i] =
            ta1_x_yy_yyyz_0[i] * fe_0 - ta1_x_yy_yyyz_1[i] * fe_0 + ta_yy_xyyyz_1[i] + ta1_x_yy_xyyyz_0[i] * pa_x[i] - ta1_x_yy_xyyyz_1[i] * pc_x[i];

        ta1_x_xyy_xyyzz_0[i] =
            ta1_x_yy_yyzz_0[i] * fe_0 - ta1_x_yy_yyzz_1[i] * fe_0 + ta_yy_xyyzz_1[i] + ta1_x_yy_xyyzz_0[i] * pa_x[i] - ta1_x_yy_xyyzz_1[i] * pc_x[i];

        ta1_x_xyy_xyzzz_0[i] =
            ta1_x_yy_yzzz_0[i] * fe_0 - ta1_x_yy_yzzz_1[i] * fe_0 + ta_yy_xyzzz_1[i] + ta1_x_yy_xyzzz_0[i] * pa_x[i] - ta1_x_yy_xyzzz_1[i] * pc_x[i];

        ta1_x_xyy_xzzzz_0[i] = ta1_x_x_xzzzz_0[i] * fe_0 - ta1_x_x_xzzzz_1[i] * fe_0 + ta1_x_xy_xzzzz_0[i] * pa_y[i] - ta1_x_xy_xzzzz_1[i] * pc_y[i];

        ta1_x_xyy_yyyyy_0[i] = ta_yy_yyyyy_1[i] + ta1_x_yy_yyyyy_0[i] * pa_x[i] - ta1_x_yy_yyyyy_1[i] * pc_x[i];

        ta1_x_xyy_yyyyz_0[i] = ta_yy_yyyyz_1[i] + ta1_x_yy_yyyyz_0[i] * pa_x[i] - ta1_x_yy_yyyyz_1[i] * pc_x[i];

        ta1_x_xyy_yyyzz_0[i] = ta_yy_yyyzz_1[i] + ta1_x_yy_yyyzz_0[i] * pa_x[i] - ta1_x_yy_yyyzz_1[i] * pc_x[i];

        ta1_x_xyy_yyzzz_0[i] = ta_yy_yyzzz_1[i] + ta1_x_yy_yyzzz_0[i] * pa_x[i] - ta1_x_yy_yyzzz_1[i] * pc_x[i];

        ta1_x_xyy_yzzzz_0[i] = ta_yy_yzzzz_1[i] + ta1_x_yy_yzzzz_0[i] * pa_x[i] - ta1_x_yy_yzzzz_1[i] * pc_x[i];

        ta1_x_xyy_zzzzz_0[i] = ta_yy_zzzzz_1[i] + ta1_x_yy_zzzzz_0[i] * pa_x[i] - ta1_x_yy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 84-105 components of targeted buffer : FH

    auto ta1_x_xyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 84);

    auto ta1_x_xyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 85);

    auto ta1_x_xyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 86);

    auto ta1_x_xyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 87);

    auto ta1_x_xyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 88);

    auto ta1_x_xyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 89);

    auto ta1_x_xyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 90);

    auto ta1_x_xyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 91);

    auto ta1_x_xyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 92);

    auto ta1_x_xyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 93);

    auto ta1_x_xyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 94);

    auto ta1_x_xyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 95);

    auto ta1_x_xyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 96);

    auto ta1_x_xyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 97);

    auto ta1_x_xyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 98);

    auto ta1_x_xyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 99);

    auto ta1_x_xyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 100);

    auto ta1_x_xyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 101);

    auto ta1_x_xyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 102);

    auto ta1_x_xyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 103);

    auto ta1_x_xyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 104);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_xy_xxxxy_0,  \
                             ta1_x_xy_xxxxy_1,  \
                             ta1_x_xy_xxxyy_0,  \
                             ta1_x_xy_xxxyy_1,  \
                             ta1_x_xy_xxyyy_0,  \
                             ta1_x_xy_xxyyy_1,  \
                             ta1_x_xy_xyyyy_0,  \
                             ta1_x_xy_xyyyy_1,  \
                             ta1_x_xy_yyyyy_0,  \
                             ta1_x_xy_yyyyy_1,  \
                             ta1_x_xyz_xxxxx_0, \
                             ta1_x_xyz_xxxxy_0, \
                             ta1_x_xyz_xxxxz_0, \
                             ta1_x_xyz_xxxyy_0, \
                             ta1_x_xyz_xxxyz_0, \
                             ta1_x_xyz_xxxzz_0, \
                             ta1_x_xyz_xxyyy_0, \
                             ta1_x_xyz_xxyyz_0, \
                             ta1_x_xyz_xxyzz_0, \
                             ta1_x_xyz_xxzzz_0, \
                             ta1_x_xyz_xyyyy_0, \
                             ta1_x_xyz_xyyyz_0, \
                             ta1_x_xyz_xyyzz_0, \
                             ta1_x_xyz_xyzzz_0, \
                             ta1_x_xyz_xzzzz_0, \
                             ta1_x_xyz_yyyyy_0, \
                             ta1_x_xyz_yyyyz_0, \
                             ta1_x_xyz_yyyzz_0, \
                             ta1_x_xyz_yyzzz_0, \
                             ta1_x_xyz_yzzzz_0, \
                             ta1_x_xyz_zzzzz_0, \
                             ta1_x_xz_xxxxx_0,  \
                             ta1_x_xz_xxxxx_1,  \
                             ta1_x_xz_xxxxz_0,  \
                             ta1_x_xz_xxxxz_1,  \
                             ta1_x_xz_xxxyz_0,  \
                             ta1_x_xz_xxxyz_1,  \
                             ta1_x_xz_xxxz_0,   \
                             ta1_x_xz_xxxz_1,   \
                             ta1_x_xz_xxxzz_0,  \
                             ta1_x_xz_xxxzz_1,  \
                             ta1_x_xz_xxyyz_0,  \
                             ta1_x_xz_xxyyz_1,  \
                             ta1_x_xz_xxyz_0,   \
                             ta1_x_xz_xxyz_1,   \
                             ta1_x_xz_xxyzz_0,  \
                             ta1_x_xz_xxyzz_1,  \
                             ta1_x_xz_xxzz_0,   \
                             ta1_x_xz_xxzz_1,   \
                             ta1_x_xz_xxzzz_0,  \
                             ta1_x_xz_xxzzz_1,  \
                             ta1_x_xz_xyyyz_0,  \
                             ta1_x_xz_xyyyz_1,  \
                             ta1_x_xz_xyyz_0,   \
                             ta1_x_xz_xyyz_1,   \
                             ta1_x_xz_xyyzz_0,  \
                             ta1_x_xz_xyyzz_1,  \
                             ta1_x_xz_xyzz_0,   \
                             ta1_x_xz_xyzz_1,   \
                             ta1_x_xz_xyzzz_0,  \
                             ta1_x_xz_xyzzz_1,  \
                             ta1_x_xz_xzzz_0,   \
                             ta1_x_xz_xzzz_1,   \
                             ta1_x_xz_xzzzz_0,  \
                             ta1_x_xz_xzzzz_1,  \
                             ta1_x_xz_zzzzz_0,  \
                             ta1_x_xz_zzzzz_1,  \
                             ta1_x_yz_yyyyz_0,  \
                             ta1_x_yz_yyyyz_1,  \
                             ta1_x_yz_yyyzz_0,  \
                             ta1_x_yz_yyyzz_1,  \
                             ta1_x_yz_yyzzz_0,  \
                             ta1_x_yz_yyzzz_1,  \
                             ta1_x_yz_yzzzz_0,  \
                             ta1_x_yz_yzzzz_1,  \
                             ta_yz_yyyyz_1,     \
                             ta_yz_yyyzz_1,     \
                             ta_yz_yyzzz_1,     \
                             ta_yz_yzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyz_xxxxx_0[i] = ta1_x_xz_xxxxx_0[i] * pa_y[i] - ta1_x_xz_xxxxx_1[i] * pc_y[i];

        ta1_x_xyz_xxxxy_0[i] = ta1_x_xy_xxxxy_0[i] * pa_z[i] - ta1_x_xy_xxxxy_1[i] * pc_z[i];

        ta1_x_xyz_xxxxz_0[i] = ta1_x_xz_xxxxz_0[i] * pa_y[i] - ta1_x_xz_xxxxz_1[i] * pc_y[i];

        ta1_x_xyz_xxxyy_0[i] = ta1_x_xy_xxxyy_0[i] * pa_z[i] - ta1_x_xy_xxxyy_1[i] * pc_z[i];

        ta1_x_xyz_xxxyz_0[i] = ta1_x_xz_xxxz_0[i] * fe_0 - ta1_x_xz_xxxz_1[i] * fe_0 + ta1_x_xz_xxxyz_0[i] * pa_y[i] - ta1_x_xz_xxxyz_1[i] * pc_y[i];

        ta1_x_xyz_xxxzz_0[i] = ta1_x_xz_xxxzz_0[i] * pa_y[i] - ta1_x_xz_xxxzz_1[i] * pc_y[i];

        ta1_x_xyz_xxyyy_0[i] = ta1_x_xy_xxyyy_0[i] * pa_z[i] - ta1_x_xy_xxyyy_1[i] * pc_z[i];

        ta1_x_xyz_xxyyz_0[i] =
            2.0 * ta1_x_xz_xxyz_0[i] * fe_0 - 2.0 * ta1_x_xz_xxyz_1[i] * fe_0 + ta1_x_xz_xxyyz_0[i] * pa_y[i] - ta1_x_xz_xxyyz_1[i] * pc_y[i];

        ta1_x_xyz_xxyzz_0[i] = ta1_x_xz_xxzz_0[i] * fe_0 - ta1_x_xz_xxzz_1[i] * fe_0 + ta1_x_xz_xxyzz_0[i] * pa_y[i] - ta1_x_xz_xxyzz_1[i] * pc_y[i];

        ta1_x_xyz_xxzzz_0[i] = ta1_x_xz_xxzzz_0[i] * pa_y[i] - ta1_x_xz_xxzzz_1[i] * pc_y[i];

        ta1_x_xyz_xyyyy_0[i] = ta1_x_xy_xyyyy_0[i] * pa_z[i] - ta1_x_xy_xyyyy_1[i] * pc_z[i];

        ta1_x_xyz_xyyyz_0[i] =
            3.0 * ta1_x_xz_xyyz_0[i] * fe_0 - 3.0 * ta1_x_xz_xyyz_1[i] * fe_0 + ta1_x_xz_xyyyz_0[i] * pa_y[i] - ta1_x_xz_xyyyz_1[i] * pc_y[i];

        ta1_x_xyz_xyyzz_0[i] =
            2.0 * ta1_x_xz_xyzz_0[i] * fe_0 - 2.0 * ta1_x_xz_xyzz_1[i] * fe_0 + ta1_x_xz_xyyzz_0[i] * pa_y[i] - ta1_x_xz_xyyzz_1[i] * pc_y[i];

        ta1_x_xyz_xyzzz_0[i] = ta1_x_xz_xzzz_0[i] * fe_0 - ta1_x_xz_xzzz_1[i] * fe_0 + ta1_x_xz_xyzzz_0[i] * pa_y[i] - ta1_x_xz_xyzzz_1[i] * pc_y[i];

        ta1_x_xyz_xzzzz_0[i] = ta1_x_xz_xzzzz_0[i] * pa_y[i] - ta1_x_xz_xzzzz_1[i] * pc_y[i];

        ta1_x_xyz_yyyyy_0[i] = ta1_x_xy_yyyyy_0[i] * pa_z[i] - ta1_x_xy_yyyyy_1[i] * pc_z[i];

        ta1_x_xyz_yyyyz_0[i] = ta_yz_yyyyz_1[i] + ta1_x_yz_yyyyz_0[i] * pa_x[i] - ta1_x_yz_yyyyz_1[i] * pc_x[i];

        ta1_x_xyz_yyyzz_0[i] = ta_yz_yyyzz_1[i] + ta1_x_yz_yyyzz_0[i] * pa_x[i] - ta1_x_yz_yyyzz_1[i] * pc_x[i];

        ta1_x_xyz_yyzzz_0[i] = ta_yz_yyzzz_1[i] + ta1_x_yz_yyzzz_0[i] * pa_x[i] - ta1_x_yz_yyzzz_1[i] * pc_x[i];

        ta1_x_xyz_yzzzz_0[i] = ta_yz_yzzzz_1[i] + ta1_x_yz_yzzzz_0[i] * pa_x[i] - ta1_x_yz_yzzzz_1[i] * pc_x[i];

        ta1_x_xyz_zzzzz_0[i] = ta1_x_xz_zzzzz_0[i] * pa_y[i] - ta1_x_xz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 105-126 components of targeted buffer : FH

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

    auto ta1_x_xzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 120);

    auto ta1_x_xzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 121);

    auto ta1_x_xzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 122);

    auto ta1_x_xzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 123);

    auto ta1_x_xzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 124);

    auto ta1_x_xzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 125);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_x_x_xxxxx_0,   \
                             ta1_x_x_xxxxx_1,   \
                             ta1_x_x_xxxxy_0,   \
                             ta1_x_x_xxxxy_1,   \
                             ta1_x_x_xxxyy_0,   \
                             ta1_x_x_xxxyy_1,   \
                             ta1_x_x_xxyyy_0,   \
                             ta1_x_x_xxyyy_1,   \
                             ta1_x_x_xyyyy_0,   \
                             ta1_x_x_xyyyy_1,   \
                             ta1_x_xz_xxxxx_0,  \
                             ta1_x_xz_xxxxx_1,  \
                             ta1_x_xz_xxxxy_0,  \
                             ta1_x_xz_xxxxy_1,  \
                             ta1_x_xz_xxxyy_0,  \
                             ta1_x_xz_xxxyy_1,  \
                             ta1_x_xz_xxyyy_0,  \
                             ta1_x_xz_xxyyy_1,  \
                             ta1_x_xz_xyyyy_0,  \
                             ta1_x_xz_xyyyy_1,  \
                             ta1_x_xzz_xxxxx_0, \
                             ta1_x_xzz_xxxxy_0, \
                             ta1_x_xzz_xxxxz_0, \
                             ta1_x_xzz_xxxyy_0, \
                             ta1_x_xzz_xxxyz_0, \
                             ta1_x_xzz_xxxzz_0, \
                             ta1_x_xzz_xxyyy_0, \
                             ta1_x_xzz_xxyyz_0, \
                             ta1_x_xzz_xxyzz_0, \
                             ta1_x_xzz_xxzzz_0, \
                             ta1_x_xzz_xyyyy_0, \
                             ta1_x_xzz_xyyyz_0, \
                             ta1_x_xzz_xyyzz_0, \
                             ta1_x_xzz_xyzzz_0, \
                             ta1_x_xzz_xzzzz_0, \
                             ta1_x_xzz_yyyyy_0, \
                             ta1_x_xzz_yyyyz_0, \
                             ta1_x_xzz_yyyzz_0, \
                             ta1_x_xzz_yyzzz_0, \
                             ta1_x_xzz_yzzzz_0, \
                             ta1_x_xzz_zzzzz_0, \
                             ta1_x_zz_xxxxz_0,  \
                             ta1_x_zz_xxxxz_1,  \
                             ta1_x_zz_xxxyz_0,  \
                             ta1_x_zz_xxxyz_1,  \
                             ta1_x_zz_xxxz_0,   \
                             ta1_x_zz_xxxz_1,   \
                             ta1_x_zz_xxxzz_0,  \
                             ta1_x_zz_xxxzz_1,  \
                             ta1_x_zz_xxyyz_0,  \
                             ta1_x_zz_xxyyz_1,  \
                             ta1_x_zz_xxyz_0,   \
                             ta1_x_zz_xxyz_1,   \
                             ta1_x_zz_xxyzz_0,  \
                             ta1_x_zz_xxyzz_1,  \
                             ta1_x_zz_xxzz_0,   \
                             ta1_x_zz_xxzz_1,   \
                             ta1_x_zz_xxzzz_0,  \
                             ta1_x_zz_xxzzz_1,  \
                             ta1_x_zz_xyyyz_0,  \
                             ta1_x_zz_xyyyz_1,  \
                             ta1_x_zz_xyyz_0,   \
                             ta1_x_zz_xyyz_1,   \
                             ta1_x_zz_xyyzz_0,  \
                             ta1_x_zz_xyyzz_1,  \
                             ta1_x_zz_xyzz_0,   \
                             ta1_x_zz_xyzz_1,   \
                             ta1_x_zz_xyzzz_0,  \
                             ta1_x_zz_xyzzz_1,  \
                             ta1_x_zz_xzzz_0,   \
                             ta1_x_zz_xzzz_1,   \
                             ta1_x_zz_xzzzz_0,  \
                             ta1_x_zz_xzzzz_1,  \
                             ta1_x_zz_yyyyy_0,  \
                             ta1_x_zz_yyyyy_1,  \
                             ta1_x_zz_yyyyz_0,  \
                             ta1_x_zz_yyyyz_1,  \
                             ta1_x_zz_yyyz_0,   \
                             ta1_x_zz_yyyz_1,   \
                             ta1_x_zz_yyyzz_0,  \
                             ta1_x_zz_yyyzz_1,  \
                             ta1_x_zz_yyzz_0,   \
                             ta1_x_zz_yyzz_1,   \
                             ta1_x_zz_yyzzz_0,  \
                             ta1_x_zz_yyzzz_1,  \
                             ta1_x_zz_yzzz_0,   \
                             ta1_x_zz_yzzz_1,   \
                             ta1_x_zz_yzzzz_0,  \
                             ta1_x_zz_yzzzz_1,  \
                             ta1_x_zz_zzzz_0,   \
                             ta1_x_zz_zzzz_1,   \
                             ta1_x_zz_zzzzz_0,  \
                             ta1_x_zz_zzzzz_1,  \
                             ta_zz_xxxxz_1,     \
                             ta_zz_xxxyz_1,     \
                             ta_zz_xxxzz_1,     \
                             ta_zz_xxyyz_1,     \
                             ta_zz_xxyzz_1,     \
                             ta_zz_xxzzz_1,     \
                             ta_zz_xyyyz_1,     \
                             ta_zz_xyyzz_1,     \
                             ta_zz_xyzzz_1,     \
                             ta_zz_xzzzz_1,     \
                             ta_zz_yyyyy_1,     \
                             ta_zz_yyyyz_1,     \
                             ta_zz_yyyzz_1,     \
                             ta_zz_yyzzz_1,     \
                             ta_zz_yzzzz_1,     \
                             ta_zz_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzz_xxxxx_0[i] = ta1_x_x_xxxxx_0[i] * fe_0 - ta1_x_x_xxxxx_1[i] * fe_0 + ta1_x_xz_xxxxx_0[i] * pa_z[i] - ta1_x_xz_xxxxx_1[i] * pc_z[i];

        ta1_x_xzz_xxxxy_0[i] = ta1_x_x_xxxxy_0[i] * fe_0 - ta1_x_x_xxxxy_1[i] * fe_0 + ta1_x_xz_xxxxy_0[i] * pa_z[i] - ta1_x_xz_xxxxy_1[i] * pc_z[i];

        ta1_x_xzz_xxxxz_0[i] = 4.0 * ta1_x_zz_xxxz_0[i] * fe_0 - 4.0 * ta1_x_zz_xxxz_1[i] * fe_0 + ta_zz_xxxxz_1[i] + ta1_x_zz_xxxxz_0[i] * pa_x[i] -
                               ta1_x_zz_xxxxz_1[i] * pc_x[i];

        ta1_x_xzz_xxxyy_0[i] = ta1_x_x_xxxyy_0[i] * fe_0 - ta1_x_x_xxxyy_1[i] * fe_0 + ta1_x_xz_xxxyy_0[i] * pa_z[i] - ta1_x_xz_xxxyy_1[i] * pc_z[i];

        ta1_x_xzz_xxxyz_0[i] = 3.0 * ta1_x_zz_xxyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyz_1[i] * fe_0 + ta_zz_xxxyz_1[i] + ta1_x_zz_xxxyz_0[i] * pa_x[i] -
                               ta1_x_zz_xxxyz_1[i] * pc_x[i];

        ta1_x_xzz_xxxzz_0[i] = 3.0 * ta1_x_zz_xxzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxzz_1[i] * fe_0 + ta_zz_xxxzz_1[i] + ta1_x_zz_xxxzz_0[i] * pa_x[i] -
                               ta1_x_zz_xxxzz_1[i] * pc_x[i];

        ta1_x_xzz_xxyyy_0[i] = ta1_x_x_xxyyy_0[i] * fe_0 - ta1_x_x_xxyyy_1[i] * fe_0 + ta1_x_xz_xxyyy_0[i] * pa_z[i] - ta1_x_xz_xxyyy_1[i] * pc_z[i];

        ta1_x_xzz_xxyyz_0[i] = 2.0 * ta1_x_zz_xyyz_0[i] * fe_0 - 2.0 * ta1_x_zz_xyyz_1[i] * fe_0 + ta_zz_xxyyz_1[i] + ta1_x_zz_xxyyz_0[i] * pa_x[i] -
                               ta1_x_zz_xxyyz_1[i] * pc_x[i];

        ta1_x_xzz_xxyzz_0[i] = 2.0 * ta1_x_zz_xyzz_0[i] * fe_0 - 2.0 * ta1_x_zz_xyzz_1[i] * fe_0 + ta_zz_xxyzz_1[i] + ta1_x_zz_xxyzz_0[i] * pa_x[i] -
                               ta1_x_zz_xxyzz_1[i] * pc_x[i];

        ta1_x_xzz_xxzzz_0[i] = 2.0 * ta1_x_zz_xzzz_0[i] * fe_0 - 2.0 * ta1_x_zz_xzzz_1[i] * fe_0 + ta_zz_xxzzz_1[i] + ta1_x_zz_xxzzz_0[i] * pa_x[i] -
                               ta1_x_zz_xxzzz_1[i] * pc_x[i];

        ta1_x_xzz_xyyyy_0[i] = ta1_x_x_xyyyy_0[i] * fe_0 - ta1_x_x_xyyyy_1[i] * fe_0 + ta1_x_xz_xyyyy_0[i] * pa_z[i] - ta1_x_xz_xyyyy_1[i] * pc_z[i];

        ta1_x_xzz_xyyyz_0[i] =
            ta1_x_zz_yyyz_0[i] * fe_0 - ta1_x_zz_yyyz_1[i] * fe_0 + ta_zz_xyyyz_1[i] + ta1_x_zz_xyyyz_0[i] * pa_x[i] - ta1_x_zz_xyyyz_1[i] * pc_x[i];

        ta1_x_xzz_xyyzz_0[i] =
            ta1_x_zz_yyzz_0[i] * fe_0 - ta1_x_zz_yyzz_1[i] * fe_0 + ta_zz_xyyzz_1[i] + ta1_x_zz_xyyzz_0[i] * pa_x[i] - ta1_x_zz_xyyzz_1[i] * pc_x[i];

        ta1_x_xzz_xyzzz_0[i] =
            ta1_x_zz_yzzz_0[i] * fe_0 - ta1_x_zz_yzzz_1[i] * fe_0 + ta_zz_xyzzz_1[i] + ta1_x_zz_xyzzz_0[i] * pa_x[i] - ta1_x_zz_xyzzz_1[i] * pc_x[i];

        ta1_x_xzz_xzzzz_0[i] =
            ta1_x_zz_zzzz_0[i] * fe_0 - ta1_x_zz_zzzz_1[i] * fe_0 + ta_zz_xzzzz_1[i] + ta1_x_zz_xzzzz_0[i] * pa_x[i] - ta1_x_zz_xzzzz_1[i] * pc_x[i];

        ta1_x_xzz_yyyyy_0[i] = ta_zz_yyyyy_1[i] + ta1_x_zz_yyyyy_0[i] * pa_x[i] - ta1_x_zz_yyyyy_1[i] * pc_x[i];

        ta1_x_xzz_yyyyz_0[i] = ta_zz_yyyyz_1[i] + ta1_x_zz_yyyyz_0[i] * pa_x[i] - ta1_x_zz_yyyyz_1[i] * pc_x[i];

        ta1_x_xzz_yyyzz_0[i] = ta_zz_yyyzz_1[i] + ta1_x_zz_yyyzz_0[i] * pa_x[i] - ta1_x_zz_yyyzz_1[i] * pc_x[i];

        ta1_x_xzz_yyzzz_0[i] = ta_zz_yyzzz_1[i] + ta1_x_zz_yyzzz_0[i] * pa_x[i] - ta1_x_zz_yyzzz_1[i] * pc_x[i];

        ta1_x_xzz_yzzzz_0[i] = ta_zz_yzzzz_1[i] + ta1_x_zz_yzzzz_0[i] * pa_x[i] - ta1_x_zz_yzzzz_1[i] * pc_x[i];

        ta1_x_xzz_zzzzz_0[i] = ta_zz_zzzzz_1[i] + ta1_x_zz_zzzzz_0[i] * pa_x[i] - ta1_x_zz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 126-147 components of targeted buffer : FH

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

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_x_y_xxxxx_0,   \
                             ta1_x_y_xxxxx_1,   \
                             ta1_x_y_xxxxy_0,   \
                             ta1_x_y_xxxxy_1,   \
                             ta1_x_y_xxxxz_0,   \
                             ta1_x_y_xxxxz_1,   \
                             ta1_x_y_xxxyy_0,   \
                             ta1_x_y_xxxyy_1,   \
                             ta1_x_y_xxxyz_0,   \
                             ta1_x_y_xxxyz_1,   \
                             ta1_x_y_xxxzz_0,   \
                             ta1_x_y_xxxzz_1,   \
                             ta1_x_y_xxyyy_0,   \
                             ta1_x_y_xxyyy_1,   \
                             ta1_x_y_xxyyz_0,   \
                             ta1_x_y_xxyyz_1,   \
                             ta1_x_y_xxyzz_0,   \
                             ta1_x_y_xxyzz_1,   \
                             ta1_x_y_xxzzz_0,   \
                             ta1_x_y_xxzzz_1,   \
                             ta1_x_y_xyyyy_0,   \
                             ta1_x_y_xyyyy_1,   \
                             ta1_x_y_xyyyz_0,   \
                             ta1_x_y_xyyyz_1,   \
                             ta1_x_y_xyyzz_0,   \
                             ta1_x_y_xyyzz_1,   \
                             ta1_x_y_xyzzz_0,   \
                             ta1_x_y_xyzzz_1,   \
                             ta1_x_y_xzzzz_0,   \
                             ta1_x_y_xzzzz_1,   \
                             ta1_x_y_yyyyy_0,   \
                             ta1_x_y_yyyyy_1,   \
                             ta1_x_y_yyyyz_0,   \
                             ta1_x_y_yyyyz_1,   \
                             ta1_x_y_yyyzz_0,   \
                             ta1_x_y_yyyzz_1,   \
                             ta1_x_y_yyzzz_0,   \
                             ta1_x_y_yyzzz_1,   \
                             ta1_x_y_yzzzz_0,   \
                             ta1_x_y_yzzzz_1,   \
                             ta1_x_y_zzzzz_0,   \
                             ta1_x_y_zzzzz_1,   \
                             ta1_x_yy_xxxx_0,   \
                             ta1_x_yy_xxxx_1,   \
                             ta1_x_yy_xxxxx_0,  \
                             ta1_x_yy_xxxxx_1,  \
                             ta1_x_yy_xxxxy_0,  \
                             ta1_x_yy_xxxxy_1,  \
                             ta1_x_yy_xxxxz_0,  \
                             ta1_x_yy_xxxxz_1,  \
                             ta1_x_yy_xxxy_0,   \
                             ta1_x_yy_xxxy_1,   \
                             ta1_x_yy_xxxyy_0,  \
                             ta1_x_yy_xxxyy_1,  \
                             ta1_x_yy_xxxyz_0,  \
                             ta1_x_yy_xxxyz_1,  \
                             ta1_x_yy_xxxz_0,   \
                             ta1_x_yy_xxxz_1,   \
                             ta1_x_yy_xxxzz_0,  \
                             ta1_x_yy_xxxzz_1,  \
                             ta1_x_yy_xxyy_0,   \
                             ta1_x_yy_xxyy_1,   \
                             ta1_x_yy_xxyyy_0,  \
                             ta1_x_yy_xxyyy_1,  \
                             ta1_x_yy_xxyyz_0,  \
                             ta1_x_yy_xxyyz_1,  \
                             ta1_x_yy_xxyz_0,   \
                             ta1_x_yy_xxyz_1,   \
                             ta1_x_yy_xxyzz_0,  \
                             ta1_x_yy_xxyzz_1,  \
                             ta1_x_yy_xxzz_0,   \
                             ta1_x_yy_xxzz_1,   \
                             ta1_x_yy_xxzzz_0,  \
                             ta1_x_yy_xxzzz_1,  \
                             ta1_x_yy_xyyy_0,   \
                             ta1_x_yy_xyyy_1,   \
                             ta1_x_yy_xyyyy_0,  \
                             ta1_x_yy_xyyyy_1,  \
                             ta1_x_yy_xyyyz_0,  \
                             ta1_x_yy_xyyyz_1,  \
                             ta1_x_yy_xyyz_0,   \
                             ta1_x_yy_xyyz_1,   \
                             ta1_x_yy_xyyzz_0,  \
                             ta1_x_yy_xyyzz_1,  \
                             ta1_x_yy_xyzz_0,   \
                             ta1_x_yy_xyzz_1,   \
                             ta1_x_yy_xyzzz_0,  \
                             ta1_x_yy_xyzzz_1,  \
                             ta1_x_yy_xzzz_0,   \
                             ta1_x_yy_xzzz_1,   \
                             ta1_x_yy_xzzzz_0,  \
                             ta1_x_yy_xzzzz_1,  \
                             ta1_x_yy_yyyy_0,   \
                             ta1_x_yy_yyyy_1,   \
                             ta1_x_yy_yyyyy_0,  \
                             ta1_x_yy_yyyyy_1,  \
                             ta1_x_yy_yyyyz_0,  \
                             ta1_x_yy_yyyyz_1,  \
                             ta1_x_yy_yyyz_0,   \
                             ta1_x_yy_yyyz_1,   \
                             ta1_x_yy_yyyzz_0,  \
                             ta1_x_yy_yyyzz_1,  \
                             ta1_x_yy_yyzz_0,   \
                             ta1_x_yy_yyzz_1,   \
                             ta1_x_yy_yyzzz_0,  \
                             ta1_x_yy_yyzzz_1,  \
                             ta1_x_yy_yzzz_0,   \
                             ta1_x_yy_yzzz_1,   \
                             ta1_x_yy_yzzzz_0,  \
                             ta1_x_yy_yzzzz_1,  \
                             ta1_x_yy_zzzz_0,   \
                             ta1_x_yy_zzzz_1,   \
                             ta1_x_yy_zzzzz_0,  \
                             ta1_x_yy_zzzzz_1,  \
                             ta1_x_yyy_xxxxx_0, \
                             ta1_x_yyy_xxxxy_0, \
                             ta1_x_yyy_xxxxz_0, \
                             ta1_x_yyy_xxxyy_0, \
                             ta1_x_yyy_xxxyz_0, \
                             ta1_x_yyy_xxxzz_0, \
                             ta1_x_yyy_xxyyy_0, \
                             ta1_x_yyy_xxyyz_0, \
                             ta1_x_yyy_xxyzz_0, \
                             ta1_x_yyy_xxzzz_0, \
                             ta1_x_yyy_xyyyy_0, \
                             ta1_x_yyy_xyyyz_0, \
                             ta1_x_yyy_xyyzz_0, \
                             ta1_x_yyy_xyzzz_0, \
                             ta1_x_yyy_xzzzz_0, \
                             ta1_x_yyy_yyyyy_0, \
                             ta1_x_yyy_yyyyz_0, \
                             ta1_x_yyy_yyyzz_0, \
                             ta1_x_yyy_yyzzz_0, \
                             ta1_x_yyy_yzzzz_0, \
                             ta1_x_yyy_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyy_xxxxx_0[i] =
            2.0 * ta1_x_y_xxxxx_0[i] * fe_0 - 2.0 * ta1_x_y_xxxxx_1[i] * fe_0 + ta1_x_yy_xxxxx_0[i] * pa_y[i] - ta1_x_yy_xxxxx_1[i] * pc_y[i];

        ta1_x_yyy_xxxxy_0[i] = 2.0 * ta1_x_y_xxxxy_0[i] * fe_0 - 2.0 * ta1_x_y_xxxxy_1[i] * fe_0 + ta1_x_yy_xxxx_0[i] * fe_0 -
                               ta1_x_yy_xxxx_1[i] * fe_0 + ta1_x_yy_xxxxy_0[i] * pa_y[i] - ta1_x_yy_xxxxy_1[i] * pc_y[i];

        ta1_x_yyy_xxxxz_0[i] =
            2.0 * ta1_x_y_xxxxz_0[i] * fe_0 - 2.0 * ta1_x_y_xxxxz_1[i] * fe_0 + ta1_x_yy_xxxxz_0[i] * pa_y[i] - ta1_x_yy_xxxxz_1[i] * pc_y[i];

        ta1_x_yyy_xxxyy_0[i] = 2.0 * ta1_x_y_xxxyy_0[i] * fe_0 - 2.0 * ta1_x_y_xxxyy_1[i] * fe_0 + 2.0 * ta1_x_yy_xxxy_0[i] * fe_0 -
                               2.0 * ta1_x_yy_xxxy_1[i] * fe_0 + ta1_x_yy_xxxyy_0[i] * pa_y[i] - ta1_x_yy_xxxyy_1[i] * pc_y[i];

        ta1_x_yyy_xxxyz_0[i] = 2.0 * ta1_x_y_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_y_xxxyz_1[i] * fe_0 + ta1_x_yy_xxxz_0[i] * fe_0 -
                               ta1_x_yy_xxxz_1[i] * fe_0 + ta1_x_yy_xxxyz_0[i] * pa_y[i] - ta1_x_yy_xxxyz_1[i] * pc_y[i];

        ta1_x_yyy_xxxzz_0[i] =
            2.0 * ta1_x_y_xxxzz_0[i] * fe_0 - 2.0 * ta1_x_y_xxxzz_1[i] * fe_0 + ta1_x_yy_xxxzz_0[i] * pa_y[i] - ta1_x_yy_xxxzz_1[i] * pc_y[i];

        ta1_x_yyy_xxyyy_0[i] = 2.0 * ta1_x_y_xxyyy_0[i] * fe_0 - 2.0 * ta1_x_y_xxyyy_1[i] * fe_0 + 3.0 * ta1_x_yy_xxyy_0[i] * fe_0 -
                               3.0 * ta1_x_yy_xxyy_1[i] * fe_0 + ta1_x_yy_xxyyy_0[i] * pa_y[i] - ta1_x_yy_xxyyy_1[i] * pc_y[i];

        ta1_x_yyy_xxyyz_0[i] = 2.0 * ta1_x_y_xxyyz_0[i] * fe_0 - 2.0 * ta1_x_y_xxyyz_1[i] * fe_0 + 2.0 * ta1_x_yy_xxyz_0[i] * fe_0 -
                               2.0 * ta1_x_yy_xxyz_1[i] * fe_0 + ta1_x_yy_xxyyz_0[i] * pa_y[i] - ta1_x_yy_xxyyz_1[i] * pc_y[i];

        ta1_x_yyy_xxyzz_0[i] = 2.0 * ta1_x_y_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_y_xxyzz_1[i] * fe_0 + ta1_x_yy_xxzz_0[i] * fe_0 -
                               ta1_x_yy_xxzz_1[i] * fe_0 + ta1_x_yy_xxyzz_0[i] * pa_y[i] - ta1_x_yy_xxyzz_1[i] * pc_y[i];

        ta1_x_yyy_xxzzz_0[i] =
            2.0 * ta1_x_y_xxzzz_0[i] * fe_0 - 2.0 * ta1_x_y_xxzzz_1[i] * fe_0 + ta1_x_yy_xxzzz_0[i] * pa_y[i] - ta1_x_yy_xxzzz_1[i] * pc_y[i];

        ta1_x_yyy_xyyyy_0[i] = 2.0 * ta1_x_y_xyyyy_0[i] * fe_0 - 2.0 * ta1_x_y_xyyyy_1[i] * fe_0 + 4.0 * ta1_x_yy_xyyy_0[i] * fe_0 -
                               4.0 * ta1_x_yy_xyyy_1[i] * fe_0 + ta1_x_yy_xyyyy_0[i] * pa_y[i] - ta1_x_yy_xyyyy_1[i] * pc_y[i];

        ta1_x_yyy_xyyyz_0[i] = 2.0 * ta1_x_y_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_y_xyyyz_1[i] * fe_0 + 3.0 * ta1_x_yy_xyyz_0[i] * fe_0 -
                               3.0 * ta1_x_yy_xyyz_1[i] * fe_0 + ta1_x_yy_xyyyz_0[i] * pa_y[i] - ta1_x_yy_xyyyz_1[i] * pc_y[i];

        ta1_x_yyy_xyyzz_0[i] = 2.0 * ta1_x_y_xyyzz_0[i] * fe_0 - 2.0 * ta1_x_y_xyyzz_1[i] * fe_0 + 2.0 * ta1_x_yy_xyzz_0[i] * fe_0 -
                               2.0 * ta1_x_yy_xyzz_1[i] * fe_0 + ta1_x_yy_xyyzz_0[i] * pa_y[i] - ta1_x_yy_xyyzz_1[i] * pc_y[i];

        ta1_x_yyy_xyzzz_0[i] = 2.0 * ta1_x_y_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_y_xyzzz_1[i] * fe_0 + ta1_x_yy_xzzz_0[i] * fe_0 -
                               ta1_x_yy_xzzz_1[i] * fe_0 + ta1_x_yy_xyzzz_0[i] * pa_y[i] - ta1_x_yy_xyzzz_1[i] * pc_y[i];

        ta1_x_yyy_xzzzz_0[i] =
            2.0 * ta1_x_y_xzzzz_0[i] * fe_0 - 2.0 * ta1_x_y_xzzzz_1[i] * fe_0 + ta1_x_yy_xzzzz_0[i] * pa_y[i] - ta1_x_yy_xzzzz_1[i] * pc_y[i];

        ta1_x_yyy_yyyyy_0[i] = 2.0 * ta1_x_y_yyyyy_0[i] * fe_0 - 2.0 * ta1_x_y_yyyyy_1[i] * fe_0 + 5.0 * ta1_x_yy_yyyy_0[i] * fe_0 -
                               5.0 * ta1_x_yy_yyyy_1[i] * fe_0 + ta1_x_yy_yyyyy_0[i] * pa_y[i] - ta1_x_yy_yyyyy_1[i] * pc_y[i];

        ta1_x_yyy_yyyyz_0[i] = 2.0 * ta1_x_y_yyyyz_0[i] * fe_0 - 2.0 * ta1_x_y_yyyyz_1[i] * fe_0 + 4.0 * ta1_x_yy_yyyz_0[i] * fe_0 -
                               4.0 * ta1_x_yy_yyyz_1[i] * fe_0 + ta1_x_yy_yyyyz_0[i] * pa_y[i] - ta1_x_yy_yyyyz_1[i] * pc_y[i];

        ta1_x_yyy_yyyzz_0[i] = 2.0 * ta1_x_y_yyyzz_0[i] * fe_0 - 2.0 * ta1_x_y_yyyzz_1[i] * fe_0 + 3.0 * ta1_x_yy_yyzz_0[i] * fe_0 -
                               3.0 * ta1_x_yy_yyzz_1[i] * fe_0 + ta1_x_yy_yyyzz_0[i] * pa_y[i] - ta1_x_yy_yyyzz_1[i] * pc_y[i];

        ta1_x_yyy_yyzzz_0[i] = 2.0 * ta1_x_y_yyzzz_0[i] * fe_0 - 2.0 * ta1_x_y_yyzzz_1[i] * fe_0 + 2.0 * ta1_x_yy_yzzz_0[i] * fe_0 -
                               2.0 * ta1_x_yy_yzzz_1[i] * fe_0 + ta1_x_yy_yyzzz_0[i] * pa_y[i] - ta1_x_yy_yyzzz_1[i] * pc_y[i];

        ta1_x_yyy_yzzzz_0[i] = 2.0 * ta1_x_y_yzzzz_0[i] * fe_0 - 2.0 * ta1_x_y_yzzzz_1[i] * fe_0 + ta1_x_yy_zzzz_0[i] * fe_0 -
                               ta1_x_yy_zzzz_1[i] * fe_0 + ta1_x_yy_yzzzz_0[i] * pa_y[i] - ta1_x_yy_yzzzz_1[i] * pc_y[i];

        ta1_x_yyy_zzzzz_0[i] =
            2.0 * ta1_x_y_zzzzz_0[i] * fe_0 - 2.0 * ta1_x_y_zzzzz_1[i] * fe_0 + ta1_x_yy_zzzzz_0[i] * pa_y[i] - ta1_x_yy_zzzzz_1[i] * pc_y[i];
    }

    // Set up 147-168 components of targeted buffer : FH

    auto ta1_x_yyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 147);

    auto ta1_x_yyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 148);

    auto ta1_x_yyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 149);

    auto ta1_x_yyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 150);

    auto ta1_x_yyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 151);

    auto ta1_x_yyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 152);

    auto ta1_x_yyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 153);

    auto ta1_x_yyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 154);

    auto ta1_x_yyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 155);

    auto ta1_x_yyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 156);

    auto ta1_x_yyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 157);

    auto ta1_x_yyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 158);

    auto ta1_x_yyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 159);

    auto ta1_x_yyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 160);

    auto ta1_x_yyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 161);

    auto ta1_x_yyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 162);

    auto ta1_x_yyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 163);

    auto ta1_x_yyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 164);

    auto ta1_x_yyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 165);

    auto ta1_x_yyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 166);

    auto ta1_x_yyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 167);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_yy_xxxxx_0,  \
                             ta1_x_yy_xxxxx_1,  \
                             ta1_x_yy_xxxxy_0,  \
                             ta1_x_yy_xxxxy_1,  \
                             ta1_x_yy_xxxy_0,   \
                             ta1_x_yy_xxxy_1,   \
                             ta1_x_yy_xxxyy_0,  \
                             ta1_x_yy_xxxyy_1,  \
                             ta1_x_yy_xxxyz_0,  \
                             ta1_x_yy_xxxyz_1,  \
                             ta1_x_yy_xxyy_0,   \
                             ta1_x_yy_xxyy_1,   \
                             ta1_x_yy_xxyyy_0,  \
                             ta1_x_yy_xxyyy_1,  \
                             ta1_x_yy_xxyyz_0,  \
                             ta1_x_yy_xxyyz_1,  \
                             ta1_x_yy_xxyz_0,   \
                             ta1_x_yy_xxyz_1,   \
                             ta1_x_yy_xxyzz_0,  \
                             ta1_x_yy_xxyzz_1,  \
                             ta1_x_yy_xyyy_0,   \
                             ta1_x_yy_xyyy_1,   \
                             ta1_x_yy_xyyyy_0,  \
                             ta1_x_yy_xyyyy_1,  \
                             ta1_x_yy_xyyyz_0,  \
                             ta1_x_yy_xyyyz_1,  \
                             ta1_x_yy_xyyz_0,   \
                             ta1_x_yy_xyyz_1,   \
                             ta1_x_yy_xyyzz_0,  \
                             ta1_x_yy_xyyzz_1,  \
                             ta1_x_yy_xyzz_0,   \
                             ta1_x_yy_xyzz_1,   \
                             ta1_x_yy_xyzzz_0,  \
                             ta1_x_yy_xyzzz_1,  \
                             ta1_x_yy_yyyy_0,   \
                             ta1_x_yy_yyyy_1,   \
                             ta1_x_yy_yyyyy_0,  \
                             ta1_x_yy_yyyyy_1,  \
                             ta1_x_yy_yyyyz_0,  \
                             ta1_x_yy_yyyyz_1,  \
                             ta1_x_yy_yyyz_0,   \
                             ta1_x_yy_yyyz_1,   \
                             ta1_x_yy_yyyzz_0,  \
                             ta1_x_yy_yyyzz_1,  \
                             ta1_x_yy_yyzz_0,   \
                             ta1_x_yy_yyzz_1,   \
                             ta1_x_yy_yyzzz_0,  \
                             ta1_x_yy_yyzzz_1,  \
                             ta1_x_yy_yzzz_0,   \
                             ta1_x_yy_yzzz_1,   \
                             ta1_x_yy_yzzzz_0,  \
                             ta1_x_yy_yzzzz_1,  \
                             ta1_x_yyz_xxxxx_0, \
                             ta1_x_yyz_xxxxy_0, \
                             ta1_x_yyz_xxxxz_0, \
                             ta1_x_yyz_xxxyy_0, \
                             ta1_x_yyz_xxxyz_0, \
                             ta1_x_yyz_xxxzz_0, \
                             ta1_x_yyz_xxyyy_0, \
                             ta1_x_yyz_xxyyz_0, \
                             ta1_x_yyz_xxyzz_0, \
                             ta1_x_yyz_xxzzz_0, \
                             ta1_x_yyz_xyyyy_0, \
                             ta1_x_yyz_xyyyz_0, \
                             ta1_x_yyz_xyyzz_0, \
                             ta1_x_yyz_xyzzz_0, \
                             ta1_x_yyz_xzzzz_0, \
                             ta1_x_yyz_yyyyy_0, \
                             ta1_x_yyz_yyyyz_0, \
                             ta1_x_yyz_yyyzz_0, \
                             ta1_x_yyz_yyzzz_0, \
                             ta1_x_yyz_yzzzz_0, \
                             ta1_x_yyz_zzzzz_0, \
                             ta1_x_yz_xxxxz_0,  \
                             ta1_x_yz_xxxxz_1,  \
                             ta1_x_yz_xxxzz_0,  \
                             ta1_x_yz_xxxzz_1,  \
                             ta1_x_yz_xxzzz_0,  \
                             ta1_x_yz_xxzzz_1,  \
                             ta1_x_yz_xzzzz_0,  \
                             ta1_x_yz_xzzzz_1,  \
                             ta1_x_yz_zzzzz_0,  \
                             ta1_x_yz_zzzzz_1,  \
                             ta1_x_z_xxxxz_0,   \
                             ta1_x_z_xxxxz_1,   \
                             ta1_x_z_xxxzz_0,   \
                             ta1_x_z_xxxzz_1,   \
                             ta1_x_z_xxzzz_0,   \
                             ta1_x_z_xxzzz_1,   \
                             ta1_x_z_xzzzz_0,   \
                             ta1_x_z_xzzzz_1,   \
                             ta1_x_z_zzzzz_0,   \
                             ta1_x_z_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyz_xxxxx_0[i] = ta1_x_yy_xxxxx_0[i] * pa_z[i] - ta1_x_yy_xxxxx_1[i] * pc_z[i];

        ta1_x_yyz_xxxxy_0[i] = ta1_x_yy_xxxxy_0[i] * pa_z[i] - ta1_x_yy_xxxxy_1[i] * pc_z[i];

        ta1_x_yyz_xxxxz_0[i] = ta1_x_z_xxxxz_0[i] * fe_0 - ta1_x_z_xxxxz_1[i] * fe_0 + ta1_x_yz_xxxxz_0[i] * pa_y[i] - ta1_x_yz_xxxxz_1[i] * pc_y[i];

        ta1_x_yyz_xxxyy_0[i] = ta1_x_yy_xxxyy_0[i] * pa_z[i] - ta1_x_yy_xxxyy_1[i] * pc_z[i];

        ta1_x_yyz_xxxyz_0[i] = ta1_x_yy_xxxy_0[i] * fe_0 - ta1_x_yy_xxxy_1[i] * fe_0 + ta1_x_yy_xxxyz_0[i] * pa_z[i] - ta1_x_yy_xxxyz_1[i] * pc_z[i];

        ta1_x_yyz_xxxzz_0[i] = ta1_x_z_xxxzz_0[i] * fe_0 - ta1_x_z_xxxzz_1[i] * fe_0 + ta1_x_yz_xxxzz_0[i] * pa_y[i] - ta1_x_yz_xxxzz_1[i] * pc_y[i];

        ta1_x_yyz_xxyyy_0[i] = ta1_x_yy_xxyyy_0[i] * pa_z[i] - ta1_x_yy_xxyyy_1[i] * pc_z[i];

        ta1_x_yyz_xxyyz_0[i] = ta1_x_yy_xxyy_0[i] * fe_0 - ta1_x_yy_xxyy_1[i] * fe_0 + ta1_x_yy_xxyyz_0[i] * pa_z[i] - ta1_x_yy_xxyyz_1[i] * pc_z[i];

        ta1_x_yyz_xxyzz_0[i] =
            2.0 * ta1_x_yy_xxyz_0[i] * fe_0 - 2.0 * ta1_x_yy_xxyz_1[i] * fe_0 + ta1_x_yy_xxyzz_0[i] * pa_z[i] - ta1_x_yy_xxyzz_1[i] * pc_z[i];

        ta1_x_yyz_xxzzz_0[i] = ta1_x_z_xxzzz_0[i] * fe_0 - ta1_x_z_xxzzz_1[i] * fe_0 + ta1_x_yz_xxzzz_0[i] * pa_y[i] - ta1_x_yz_xxzzz_1[i] * pc_y[i];

        ta1_x_yyz_xyyyy_0[i] = ta1_x_yy_xyyyy_0[i] * pa_z[i] - ta1_x_yy_xyyyy_1[i] * pc_z[i];

        ta1_x_yyz_xyyyz_0[i] = ta1_x_yy_xyyy_0[i] * fe_0 - ta1_x_yy_xyyy_1[i] * fe_0 + ta1_x_yy_xyyyz_0[i] * pa_z[i] - ta1_x_yy_xyyyz_1[i] * pc_z[i];

        ta1_x_yyz_xyyzz_0[i] =
            2.0 * ta1_x_yy_xyyz_0[i] * fe_0 - 2.0 * ta1_x_yy_xyyz_1[i] * fe_0 + ta1_x_yy_xyyzz_0[i] * pa_z[i] - ta1_x_yy_xyyzz_1[i] * pc_z[i];

        ta1_x_yyz_xyzzz_0[i] =
            3.0 * ta1_x_yy_xyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xyzz_1[i] * fe_0 + ta1_x_yy_xyzzz_0[i] * pa_z[i] - ta1_x_yy_xyzzz_1[i] * pc_z[i];

        ta1_x_yyz_xzzzz_0[i] = ta1_x_z_xzzzz_0[i] * fe_0 - ta1_x_z_xzzzz_1[i] * fe_0 + ta1_x_yz_xzzzz_0[i] * pa_y[i] - ta1_x_yz_xzzzz_1[i] * pc_y[i];

        ta1_x_yyz_yyyyy_0[i] = ta1_x_yy_yyyyy_0[i] * pa_z[i] - ta1_x_yy_yyyyy_1[i] * pc_z[i];

        ta1_x_yyz_yyyyz_0[i] = ta1_x_yy_yyyy_0[i] * fe_0 - ta1_x_yy_yyyy_1[i] * fe_0 + ta1_x_yy_yyyyz_0[i] * pa_z[i] - ta1_x_yy_yyyyz_1[i] * pc_z[i];

        ta1_x_yyz_yyyzz_0[i] =
            2.0 * ta1_x_yy_yyyz_0[i] * fe_0 - 2.0 * ta1_x_yy_yyyz_1[i] * fe_0 + ta1_x_yy_yyyzz_0[i] * pa_z[i] - ta1_x_yy_yyyzz_1[i] * pc_z[i];

        ta1_x_yyz_yyzzz_0[i] =
            3.0 * ta1_x_yy_yyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyzz_1[i] * fe_0 + ta1_x_yy_yyzzz_0[i] * pa_z[i] - ta1_x_yy_yyzzz_1[i] * pc_z[i];

        ta1_x_yyz_yzzzz_0[i] =
            4.0 * ta1_x_yy_yzzz_0[i] * fe_0 - 4.0 * ta1_x_yy_yzzz_1[i] * fe_0 + ta1_x_yy_yzzzz_0[i] * pa_z[i] - ta1_x_yy_yzzzz_1[i] * pc_z[i];

        ta1_x_yyz_zzzzz_0[i] = ta1_x_z_zzzzz_0[i] * fe_0 - ta1_x_z_zzzzz_1[i] * fe_0 + ta1_x_yz_zzzzz_0[i] * pa_y[i] - ta1_x_yz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 168-189 components of targeted buffer : FH

    auto ta1_x_yzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 168);

    auto ta1_x_yzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 169);

    auto ta1_x_yzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 170);

    auto ta1_x_yzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 171);

    auto ta1_x_yzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 172);

    auto ta1_x_yzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 173);

    auto ta1_x_yzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 174);

    auto ta1_x_yzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 175);

    auto ta1_x_yzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 176);

    auto ta1_x_yzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 177);

    auto ta1_x_yzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 178);

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

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_x_yzz_xxxxx_0, \
                             ta1_x_yzz_xxxxy_0, \
                             ta1_x_yzz_xxxxz_0, \
                             ta1_x_yzz_xxxyy_0, \
                             ta1_x_yzz_xxxyz_0, \
                             ta1_x_yzz_xxxzz_0, \
                             ta1_x_yzz_xxyyy_0, \
                             ta1_x_yzz_xxyyz_0, \
                             ta1_x_yzz_xxyzz_0, \
                             ta1_x_yzz_xxzzz_0, \
                             ta1_x_yzz_xyyyy_0, \
                             ta1_x_yzz_xyyyz_0, \
                             ta1_x_yzz_xyyzz_0, \
                             ta1_x_yzz_xyzzz_0, \
                             ta1_x_yzz_xzzzz_0, \
                             ta1_x_yzz_yyyyy_0, \
                             ta1_x_yzz_yyyyz_0, \
                             ta1_x_yzz_yyyzz_0, \
                             ta1_x_yzz_yyzzz_0, \
                             ta1_x_yzz_yzzzz_0, \
                             ta1_x_yzz_zzzzz_0, \
                             ta1_x_zz_xxxx_0,   \
                             ta1_x_zz_xxxx_1,   \
                             ta1_x_zz_xxxxx_0,  \
                             ta1_x_zz_xxxxx_1,  \
                             ta1_x_zz_xxxxy_0,  \
                             ta1_x_zz_xxxxy_1,  \
                             ta1_x_zz_xxxxz_0,  \
                             ta1_x_zz_xxxxz_1,  \
                             ta1_x_zz_xxxy_0,   \
                             ta1_x_zz_xxxy_1,   \
                             ta1_x_zz_xxxyy_0,  \
                             ta1_x_zz_xxxyy_1,  \
                             ta1_x_zz_xxxyz_0,  \
                             ta1_x_zz_xxxyz_1,  \
                             ta1_x_zz_xxxz_0,   \
                             ta1_x_zz_xxxz_1,   \
                             ta1_x_zz_xxxzz_0,  \
                             ta1_x_zz_xxxzz_1,  \
                             ta1_x_zz_xxyy_0,   \
                             ta1_x_zz_xxyy_1,   \
                             ta1_x_zz_xxyyy_0,  \
                             ta1_x_zz_xxyyy_1,  \
                             ta1_x_zz_xxyyz_0,  \
                             ta1_x_zz_xxyyz_1,  \
                             ta1_x_zz_xxyz_0,   \
                             ta1_x_zz_xxyz_1,   \
                             ta1_x_zz_xxyzz_0,  \
                             ta1_x_zz_xxyzz_1,  \
                             ta1_x_zz_xxzz_0,   \
                             ta1_x_zz_xxzz_1,   \
                             ta1_x_zz_xxzzz_0,  \
                             ta1_x_zz_xxzzz_1,  \
                             ta1_x_zz_xyyy_0,   \
                             ta1_x_zz_xyyy_1,   \
                             ta1_x_zz_xyyyy_0,  \
                             ta1_x_zz_xyyyy_1,  \
                             ta1_x_zz_xyyyz_0,  \
                             ta1_x_zz_xyyyz_1,  \
                             ta1_x_zz_xyyz_0,   \
                             ta1_x_zz_xyyz_1,   \
                             ta1_x_zz_xyyzz_0,  \
                             ta1_x_zz_xyyzz_1,  \
                             ta1_x_zz_xyzz_0,   \
                             ta1_x_zz_xyzz_1,   \
                             ta1_x_zz_xyzzz_0,  \
                             ta1_x_zz_xyzzz_1,  \
                             ta1_x_zz_xzzz_0,   \
                             ta1_x_zz_xzzz_1,   \
                             ta1_x_zz_xzzzz_0,  \
                             ta1_x_zz_xzzzz_1,  \
                             ta1_x_zz_yyyy_0,   \
                             ta1_x_zz_yyyy_1,   \
                             ta1_x_zz_yyyyy_0,  \
                             ta1_x_zz_yyyyy_1,  \
                             ta1_x_zz_yyyyz_0,  \
                             ta1_x_zz_yyyyz_1,  \
                             ta1_x_zz_yyyz_0,   \
                             ta1_x_zz_yyyz_1,   \
                             ta1_x_zz_yyyzz_0,  \
                             ta1_x_zz_yyyzz_1,  \
                             ta1_x_zz_yyzz_0,   \
                             ta1_x_zz_yyzz_1,   \
                             ta1_x_zz_yyzzz_0,  \
                             ta1_x_zz_yyzzz_1,  \
                             ta1_x_zz_yzzz_0,   \
                             ta1_x_zz_yzzz_1,   \
                             ta1_x_zz_yzzzz_0,  \
                             ta1_x_zz_yzzzz_1,  \
                             ta1_x_zz_zzzz_0,   \
                             ta1_x_zz_zzzz_1,   \
                             ta1_x_zz_zzzzz_0,  \
                             ta1_x_zz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzz_xxxxx_0[i] = ta1_x_zz_xxxxx_0[i] * pa_y[i] - ta1_x_zz_xxxxx_1[i] * pc_y[i];

        ta1_x_yzz_xxxxy_0[i] = ta1_x_zz_xxxx_0[i] * fe_0 - ta1_x_zz_xxxx_1[i] * fe_0 + ta1_x_zz_xxxxy_0[i] * pa_y[i] - ta1_x_zz_xxxxy_1[i] * pc_y[i];

        ta1_x_yzz_xxxxz_0[i] = ta1_x_zz_xxxxz_0[i] * pa_y[i] - ta1_x_zz_xxxxz_1[i] * pc_y[i];

        ta1_x_yzz_xxxyy_0[i] =
            2.0 * ta1_x_zz_xxxy_0[i] * fe_0 - 2.0 * ta1_x_zz_xxxy_1[i] * fe_0 + ta1_x_zz_xxxyy_0[i] * pa_y[i] - ta1_x_zz_xxxyy_1[i] * pc_y[i];

        ta1_x_yzz_xxxyz_0[i] = ta1_x_zz_xxxz_0[i] * fe_0 - ta1_x_zz_xxxz_1[i] * fe_0 + ta1_x_zz_xxxyz_0[i] * pa_y[i] - ta1_x_zz_xxxyz_1[i] * pc_y[i];

        ta1_x_yzz_xxxzz_0[i] = ta1_x_zz_xxxzz_0[i] * pa_y[i] - ta1_x_zz_xxxzz_1[i] * pc_y[i];

        ta1_x_yzz_xxyyy_0[i] =
            3.0 * ta1_x_zz_xxyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyy_1[i] * fe_0 + ta1_x_zz_xxyyy_0[i] * pa_y[i] - ta1_x_zz_xxyyy_1[i] * pc_y[i];

        ta1_x_yzz_xxyyz_0[i] =
            2.0 * ta1_x_zz_xxyz_0[i] * fe_0 - 2.0 * ta1_x_zz_xxyz_1[i] * fe_0 + ta1_x_zz_xxyyz_0[i] * pa_y[i] - ta1_x_zz_xxyyz_1[i] * pc_y[i];

        ta1_x_yzz_xxyzz_0[i] = ta1_x_zz_xxzz_0[i] * fe_0 - ta1_x_zz_xxzz_1[i] * fe_0 + ta1_x_zz_xxyzz_0[i] * pa_y[i] - ta1_x_zz_xxyzz_1[i] * pc_y[i];

        ta1_x_yzz_xxzzz_0[i] = ta1_x_zz_xxzzz_0[i] * pa_y[i] - ta1_x_zz_xxzzz_1[i] * pc_y[i];

        ta1_x_yzz_xyyyy_0[i] =
            4.0 * ta1_x_zz_xyyy_0[i] * fe_0 - 4.0 * ta1_x_zz_xyyy_1[i] * fe_0 + ta1_x_zz_xyyyy_0[i] * pa_y[i] - ta1_x_zz_xyyyy_1[i] * pc_y[i];

        ta1_x_yzz_xyyyz_0[i] =
            3.0 * ta1_x_zz_xyyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xyyz_1[i] * fe_0 + ta1_x_zz_xyyyz_0[i] * pa_y[i] - ta1_x_zz_xyyyz_1[i] * pc_y[i];

        ta1_x_yzz_xyyzz_0[i] =
            2.0 * ta1_x_zz_xyzz_0[i] * fe_0 - 2.0 * ta1_x_zz_xyzz_1[i] * fe_0 + ta1_x_zz_xyyzz_0[i] * pa_y[i] - ta1_x_zz_xyyzz_1[i] * pc_y[i];

        ta1_x_yzz_xyzzz_0[i] = ta1_x_zz_xzzz_0[i] * fe_0 - ta1_x_zz_xzzz_1[i] * fe_0 + ta1_x_zz_xyzzz_0[i] * pa_y[i] - ta1_x_zz_xyzzz_1[i] * pc_y[i];

        ta1_x_yzz_xzzzz_0[i] = ta1_x_zz_xzzzz_0[i] * pa_y[i] - ta1_x_zz_xzzzz_1[i] * pc_y[i];

        ta1_x_yzz_yyyyy_0[i] =
            5.0 * ta1_x_zz_yyyy_0[i] * fe_0 - 5.0 * ta1_x_zz_yyyy_1[i] * fe_0 + ta1_x_zz_yyyyy_0[i] * pa_y[i] - ta1_x_zz_yyyyy_1[i] * pc_y[i];

        ta1_x_yzz_yyyyz_0[i] =
            4.0 * ta1_x_zz_yyyz_0[i] * fe_0 - 4.0 * ta1_x_zz_yyyz_1[i] * fe_0 + ta1_x_zz_yyyyz_0[i] * pa_y[i] - ta1_x_zz_yyyyz_1[i] * pc_y[i];

        ta1_x_yzz_yyyzz_0[i] =
            3.0 * ta1_x_zz_yyzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyzz_1[i] * fe_0 + ta1_x_zz_yyyzz_0[i] * pa_y[i] - ta1_x_zz_yyyzz_1[i] * pc_y[i];

        ta1_x_yzz_yyzzz_0[i] =
            2.0 * ta1_x_zz_yzzz_0[i] * fe_0 - 2.0 * ta1_x_zz_yzzz_1[i] * fe_0 + ta1_x_zz_yyzzz_0[i] * pa_y[i] - ta1_x_zz_yyzzz_1[i] * pc_y[i];

        ta1_x_yzz_yzzzz_0[i] = ta1_x_zz_zzzz_0[i] * fe_0 - ta1_x_zz_zzzz_1[i] * fe_0 + ta1_x_zz_yzzzz_0[i] * pa_y[i] - ta1_x_zz_yzzzz_1[i] * pc_y[i];

        ta1_x_yzz_zzzzz_0[i] = ta1_x_zz_zzzzz_0[i] * pa_y[i] - ta1_x_zz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 189-210 components of targeted buffer : FH

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

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_x_z_xxxxx_0,   \
                             ta1_x_z_xxxxx_1,   \
                             ta1_x_z_xxxxy_0,   \
                             ta1_x_z_xxxxy_1,   \
                             ta1_x_z_xxxxz_0,   \
                             ta1_x_z_xxxxz_1,   \
                             ta1_x_z_xxxyy_0,   \
                             ta1_x_z_xxxyy_1,   \
                             ta1_x_z_xxxyz_0,   \
                             ta1_x_z_xxxyz_1,   \
                             ta1_x_z_xxxzz_0,   \
                             ta1_x_z_xxxzz_1,   \
                             ta1_x_z_xxyyy_0,   \
                             ta1_x_z_xxyyy_1,   \
                             ta1_x_z_xxyyz_0,   \
                             ta1_x_z_xxyyz_1,   \
                             ta1_x_z_xxyzz_0,   \
                             ta1_x_z_xxyzz_1,   \
                             ta1_x_z_xxzzz_0,   \
                             ta1_x_z_xxzzz_1,   \
                             ta1_x_z_xyyyy_0,   \
                             ta1_x_z_xyyyy_1,   \
                             ta1_x_z_xyyyz_0,   \
                             ta1_x_z_xyyyz_1,   \
                             ta1_x_z_xyyzz_0,   \
                             ta1_x_z_xyyzz_1,   \
                             ta1_x_z_xyzzz_0,   \
                             ta1_x_z_xyzzz_1,   \
                             ta1_x_z_xzzzz_0,   \
                             ta1_x_z_xzzzz_1,   \
                             ta1_x_z_yyyyy_0,   \
                             ta1_x_z_yyyyy_1,   \
                             ta1_x_z_yyyyz_0,   \
                             ta1_x_z_yyyyz_1,   \
                             ta1_x_z_yyyzz_0,   \
                             ta1_x_z_yyyzz_1,   \
                             ta1_x_z_yyzzz_0,   \
                             ta1_x_z_yyzzz_1,   \
                             ta1_x_z_yzzzz_0,   \
                             ta1_x_z_yzzzz_1,   \
                             ta1_x_z_zzzzz_0,   \
                             ta1_x_z_zzzzz_1,   \
                             ta1_x_zz_xxxx_0,   \
                             ta1_x_zz_xxxx_1,   \
                             ta1_x_zz_xxxxx_0,  \
                             ta1_x_zz_xxxxx_1,  \
                             ta1_x_zz_xxxxy_0,  \
                             ta1_x_zz_xxxxy_1,  \
                             ta1_x_zz_xxxxz_0,  \
                             ta1_x_zz_xxxxz_1,  \
                             ta1_x_zz_xxxy_0,   \
                             ta1_x_zz_xxxy_1,   \
                             ta1_x_zz_xxxyy_0,  \
                             ta1_x_zz_xxxyy_1,  \
                             ta1_x_zz_xxxyz_0,  \
                             ta1_x_zz_xxxyz_1,  \
                             ta1_x_zz_xxxz_0,   \
                             ta1_x_zz_xxxz_1,   \
                             ta1_x_zz_xxxzz_0,  \
                             ta1_x_zz_xxxzz_1,  \
                             ta1_x_zz_xxyy_0,   \
                             ta1_x_zz_xxyy_1,   \
                             ta1_x_zz_xxyyy_0,  \
                             ta1_x_zz_xxyyy_1,  \
                             ta1_x_zz_xxyyz_0,  \
                             ta1_x_zz_xxyyz_1,  \
                             ta1_x_zz_xxyz_0,   \
                             ta1_x_zz_xxyz_1,   \
                             ta1_x_zz_xxyzz_0,  \
                             ta1_x_zz_xxyzz_1,  \
                             ta1_x_zz_xxzz_0,   \
                             ta1_x_zz_xxzz_1,   \
                             ta1_x_zz_xxzzz_0,  \
                             ta1_x_zz_xxzzz_1,  \
                             ta1_x_zz_xyyy_0,   \
                             ta1_x_zz_xyyy_1,   \
                             ta1_x_zz_xyyyy_0,  \
                             ta1_x_zz_xyyyy_1,  \
                             ta1_x_zz_xyyyz_0,  \
                             ta1_x_zz_xyyyz_1,  \
                             ta1_x_zz_xyyz_0,   \
                             ta1_x_zz_xyyz_1,   \
                             ta1_x_zz_xyyzz_0,  \
                             ta1_x_zz_xyyzz_1,  \
                             ta1_x_zz_xyzz_0,   \
                             ta1_x_zz_xyzz_1,   \
                             ta1_x_zz_xyzzz_0,  \
                             ta1_x_zz_xyzzz_1,  \
                             ta1_x_zz_xzzz_0,   \
                             ta1_x_zz_xzzz_1,   \
                             ta1_x_zz_xzzzz_0,  \
                             ta1_x_zz_xzzzz_1,  \
                             ta1_x_zz_yyyy_0,   \
                             ta1_x_zz_yyyy_1,   \
                             ta1_x_zz_yyyyy_0,  \
                             ta1_x_zz_yyyyy_1,  \
                             ta1_x_zz_yyyyz_0,  \
                             ta1_x_zz_yyyyz_1,  \
                             ta1_x_zz_yyyz_0,   \
                             ta1_x_zz_yyyz_1,   \
                             ta1_x_zz_yyyzz_0,  \
                             ta1_x_zz_yyyzz_1,  \
                             ta1_x_zz_yyzz_0,   \
                             ta1_x_zz_yyzz_1,   \
                             ta1_x_zz_yyzzz_0,  \
                             ta1_x_zz_yyzzz_1,  \
                             ta1_x_zz_yzzz_0,   \
                             ta1_x_zz_yzzz_1,   \
                             ta1_x_zz_yzzzz_0,  \
                             ta1_x_zz_yzzzz_1,  \
                             ta1_x_zz_zzzz_0,   \
                             ta1_x_zz_zzzz_1,   \
                             ta1_x_zz_zzzzz_0,  \
                             ta1_x_zz_zzzzz_1,  \
                             ta1_x_zzz_xxxxx_0, \
                             ta1_x_zzz_xxxxy_0, \
                             ta1_x_zzz_xxxxz_0, \
                             ta1_x_zzz_xxxyy_0, \
                             ta1_x_zzz_xxxyz_0, \
                             ta1_x_zzz_xxxzz_0, \
                             ta1_x_zzz_xxyyy_0, \
                             ta1_x_zzz_xxyyz_0, \
                             ta1_x_zzz_xxyzz_0, \
                             ta1_x_zzz_xxzzz_0, \
                             ta1_x_zzz_xyyyy_0, \
                             ta1_x_zzz_xyyyz_0, \
                             ta1_x_zzz_xyyzz_0, \
                             ta1_x_zzz_xyzzz_0, \
                             ta1_x_zzz_xzzzz_0, \
                             ta1_x_zzz_yyyyy_0, \
                             ta1_x_zzz_yyyyz_0, \
                             ta1_x_zzz_yyyzz_0, \
                             ta1_x_zzz_yyzzz_0, \
                             ta1_x_zzz_yzzzz_0, \
                             ta1_x_zzz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzz_xxxxx_0[i] =
            2.0 * ta1_x_z_xxxxx_0[i] * fe_0 - 2.0 * ta1_x_z_xxxxx_1[i] * fe_0 + ta1_x_zz_xxxxx_0[i] * pa_z[i] - ta1_x_zz_xxxxx_1[i] * pc_z[i];

        ta1_x_zzz_xxxxy_0[i] =
            2.0 * ta1_x_z_xxxxy_0[i] * fe_0 - 2.0 * ta1_x_z_xxxxy_1[i] * fe_0 + ta1_x_zz_xxxxy_0[i] * pa_z[i] - ta1_x_zz_xxxxy_1[i] * pc_z[i];

        ta1_x_zzz_xxxxz_0[i] = 2.0 * ta1_x_z_xxxxz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxxz_1[i] * fe_0 + ta1_x_zz_xxxx_0[i] * fe_0 -
                               ta1_x_zz_xxxx_1[i] * fe_0 + ta1_x_zz_xxxxz_0[i] * pa_z[i] - ta1_x_zz_xxxxz_1[i] * pc_z[i];

        ta1_x_zzz_xxxyy_0[i] =
            2.0 * ta1_x_z_xxxyy_0[i] * fe_0 - 2.0 * ta1_x_z_xxxyy_1[i] * fe_0 + ta1_x_zz_xxxyy_0[i] * pa_z[i] - ta1_x_zz_xxxyy_1[i] * pc_z[i];

        ta1_x_zzz_xxxyz_0[i] = 2.0 * ta1_x_z_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxyz_1[i] * fe_0 + ta1_x_zz_xxxy_0[i] * fe_0 -
                               ta1_x_zz_xxxy_1[i] * fe_0 + ta1_x_zz_xxxyz_0[i] * pa_z[i] - ta1_x_zz_xxxyz_1[i] * pc_z[i];

        ta1_x_zzz_xxxzz_0[i] = 2.0 * ta1_x_z_xxxzz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxzz_1[i] * fe_0 + 2.0 * ta1_x_zz_xxxz_0[i] * fe_0 -
                               2.0 * ta1_x_zz_xxxz_1[i] * fe_0 + ta1_x_zz_xxxzz_0[i] * pa_z[i] - ta1_x_zz_xxxzz_1[i] * pc_z[i];

        ta1_x_zzz_xxyyy_0[i] =
            2.0 * ta1_x_z_xxyyy_0[i] * fe_0 - 2.0 * ta1_x_z_xxyyy_1[i] * fe_0 + ta1_x_zz_xxyyy_0[i] * pa_z[i] - ta1_x_zz_xxyyy_1[i] * pc_z[i];

        ta1_x_zzz_xxyyz_0[i] = 2.0 * ta1_x_z_xxyyz_0[i] * fe_0 - 2.0 * ta1_x_z_xxyyz_1[i] * fe_0 + ta1_x_zz_xxyy_0[i] * fe_0 -
                               ta1_x_zz_xxyy_1[i] * fe_0 + ta1_x_zz_xxyyz_0[i] * pa_z[i] - ta1_x_zz_xxyyz_1[i] * pc_z[i];

        ta1_x_zzz_xxyzz_0[i] = 2.0 * ta1_x_z_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_z_xxyzz_1[i] * fe_0 + 2.0 * ta1_x_zz_xxyz_0[i] * fe_0 -
                               2.0 * ta1_x_zz_xxyz_1[i] * fe_0 + ta1_x_zz_xxyzz_0[i] * pa_z[i] - ta1_x_zz_xxyzz_1[i] * pc_z[i];

        ta1_x_zzz_xxzzz_0[i] = 2.0 * ta1_x_z_xxzzz_0[i] * fe_0 - 2.0 * ta1_x_z_xxzzz_1[i] * fe_0 + 3.0 * ta1_x_zz_xxzz_0[i] * fe_0 -
                               3.0 * ta1_x_zz_xxzz_1[i] * fe_0 + ta1_x_zz_xxzzz_0[i] * pa_z[i] - ta1_x_zz_xxzzz_1[i] * pc_z[i];

        ta1_x_zzz_xyyyy_0[i] =
            2.0 * ta1_x_z_xyyyy_0[i] * fe_0 - 2.0 * ta1_x_z_xyyyy_1[i] * fe_0 + ta1_x_zz_xyyyy_0[i] * pa_z[i] - ta1_x_zz_xyyyy_1[i] * pc_z[i];

        ta1_x_zzz_xyyyz_0[i] = 2.0 * ta1_x_z_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_z_xyyyz_1[i] * fe_0 + ta1_x_zz_xyyy_0[i] * fe_0 -
                               ta1_x_zz_xyyy_1[i] * fe_0 + ta1_x_zz_xyyyz_0[i] * pa_z[i] - ta1_x_zz_xyyyz_1[i] * pc_z[i];

        ta1_x_zzz_xyyzz_0[i] = 2.0 * ta1_x_z_xyyzz_0[i] * fe_0 - 2.0 * ta1_x_z_xyyzz_1[i] * fe_0 + 2.0 * ta1_x_zz_xyyz_0[i] * fe_0 -
                               2.0 * ta1_x_zz_xyyz_1[i] * fe_0 + ta1_x_zz_xyyzz_0[i] * pa_z[i] - ta1_x_zz_xyyzz_1[i] * pc_z[i];

        ta1_x_zzz_xyzzz_0[i] = 2.0 * ta1_x_z_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_z_xyzzz_1[i] * fe_0 + 3.0 * ta1_x_zz_xyzz_0[i] * fe_0 -
                               3.0 * ta1_x_zz_xyzz_1[i] * fe_0 + ta1_x_zz_xyzzz_0[i] * pa_z[i] - ta1_x_zz_xyzzz_1[i] * pc_z[i];

        ta1_x_zzz_xzzzz_0[i] = 2.0 * ta1_x_z_xzzzz_0[i] * fe_0 - 2.0 * ta1_x_z_xzzzz_1[i] * fe_0 + 4.0 * ta1_x_zz_xzzz_0[i] * fe_0 -
                               4.0 * ta1_x_zz_xzzz_1[i] * fe_0 + ta1_x_zz_xzzzz_0[i] * pa_z[i] - ta1_x_zz_xzzzz_1[i] * pc_z[i];

        ta1_x_zzz_yyyyy_0[i] =
            2.0 * ta1_x_z_yyyyy_0[i] * fe_0 - 2.0 * ta1_x_z_yyyyy_1[i] * fe_0 + ta1_x_zz_yyyyy_0[i] * pa_z[i] - ta1_x_zz_yyyyy_1[i] * pc_z[i];

        ta1_x_zzz_yyyyz_0[i] = 2.0 * ta1_x_z_yyyyz_0[i] * fe_0 - 2.0 * ta1_x_z_yyyyz_1[i] * fe_0 + ta1_x_zz_yyyy_0[i] * fe_0 -
                               ta1_x_zz_yyyy_1[i] * fe_0 + ta1_x_zz_yyyyz_0[i] * pa_z[i] - ta1_x_zz_yyyyz_1[i] * pc_z[i];

        ta1_x_zzz_yyyzz_0[i] = 2.0 * ta1_x_z_yyyzz_0[i] * fe_0 - 2.0 * ta1_x_z_yyyzz_1[i] * fe_0 + 2.0 * ta1_x_zz_yyyz_0[i] * fe_0 -
                               2.0 * ta1_x_zz_yyyz_1[i] * fe_0 + ta1_x_zz_yyyzz_0[i] * pa_z[i] - ta1_x_zz_yyyzz_1[i] * pc_z[i];

        ta1_x_zzz_yyzzz_0[i] = 2.0 * ta1_x_z_yyzzz_0[i] * fe_0 - 2.0 * ta1_x_z_yyzzz_1[i] * fe_0 + 3.0 * ta1_x_zz_yyzz_0[i] * fe_0 -
                               3.0 * ta1_x_zz_yyzz_1[i] * fe_0 + ta1_x_zz_yyzzz_0[i] * pa_z[i] - ta1_x_zz_yyzzz_1[i] * pc_z[i];

        ta1_x_zzz_yzzzz_0[i] = 2.0 * ta1_x_z_yzzzz_0[i] * fe_0 - 2.0 * ta1_x_z_yzzzz_1[i] * fe_0 + 4.0 * ta1_x_zz_yzzz_0[i] * fe_0 -
                               4.0 * ta1_x_zz_yzzz_1[i] * fe_0 + ta1_x_zz_yzzzz_0[i] * pa_z[i] - ta1_x_zz_yzzzz_1[i] * pc_z[i];

        ta1_x_zzz_zzzzz_0[i] = 2.0 * ta1_x_z_zzzzz_0[i] * fe_0 - 2.0 * ta1_x_z_zzzzz_1[i] * fe_0 + 5.0 * ta1_x_zz_zzzz_0[i] * fe_0 -
                               5.0 * ta1_x_zz_zzzz_1[i] * fe_0 + ta1_x_zz_zzzzz_0[i] * pa_z[i] - ta1_x_zz_zzzzz_1[i] * pc_z[i];
    }

    // Set up 210-231 components of targeted buffer : FH

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

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_y_x_xxxxx_0,   \
                             ta1_y_x_xxxxx_1,   \
                             ta1_y_x_xxxxy_0,   \
                             ta1_y_x_xxxxy_1,   \
                             ta1_y_x_xxxxz_0,   \
                             ta1_y_x_xxxxz_1,   \
                             ta1_y_x_xxxyy_0,   \
                             ta1_y_x_xxxyy_1,   \
                             ta1_y_x_xxxyz_0,   \
                             ta1_y_x_xxxyz_1,   \
                             ta1_y_x_xxxzz_0,   \
                             ta1_y_x_xxxzz_1,   \
                             ta1_y_x_xxyyy_0,   \
                             ta1_y_x_xxyyy_1,   \
                             ta1_y_x_xxyyz_0,   \
                             ta1_y_x_xxyyz_1,   \
                             ta1_y_x_xxyzz_0,   \
                             ta1_y_x_xxyzz_1,   \
                             ta1_y_x_xxzzz_0,   \
                             ta1_y_x_xxzzz_1,   \
                             ta1_y_x_xyyyy_0,   \
                             ta1_y_x_xyyyy_1,   \
                             ta1_y_x_xyyyz_0,   \
                             ta1_y_x_xyyyz_1,   \
                             ta1_y_x_xyyzz_0,   \
                             ta1_y_x_xyyzz_1,   \
                             ta1_y_x_xyzzz_0,   \
                             ta1_y_x_xyzzz_1,   \
                             ta1_y_x_xzzzz_0,   \
                             ta1_y_x_xzzzz_1,   \
                             ta1_y_x_yyyyy_0,   \
                             ta1_y_x_yyyyy_1,   \
                             ta1_y_x_yyyyz_0,   \
                             ta1_y_x_yyyyz_1,   \
                             ta1_y_x_yyyzz_0,   \
                             ta1_y_x_yyyzz_1,   \
                             ta1_y_x_yyzzz_0,   \
                             ta1_y_x_yyzzz_1,   \
                             ta1_y_x_yzzzz_0,   \
                             ta1_y_x_yzzzz_1,   \
                             ta1_y_x_zzzzz_0,   \
                             ta1_y_x_zzzzz_1,   \
                             ta1_y_xx_xxxx_0,   \
                             ta1_y_xx_xxxx_1,   \
                             ta1_y_xx_xxxxx_0,  \
                             ta1_y_xx_xxxxx_1,  \
                             ta1_y_xx_xxxxy_0,  \
                             ta1_y_xx_xxxxy_1,  \
                             ta1_y_xx_xxxxz_0,  \
                             ta1_y_xx_xxxxz_1,  \
                             ta1_y_xx_xxxy_0,   \
                             ta1_y_xx_xxxy_1,   \
                             ta1_y_xx_xxxyy_0,  \
                             ta1_y_xx_xxxyy_1,  \
                             ta1_y_xx_xxxyz_0,  \
                             ta1_y_xx_xxxyz_1,  \
                             ta1_y_xx_xxxz_0,   \
                             ta1_y_xx_xxxz_1,   \
                             ta1_y_xx_xxxzz_0,  \
                             ta1_y_xx_xxxzz_1,  \
                             ta1_y_xx_xxyy_0,   \
                             ta1_y_xx_xxyy_1,   \
                             ta1_y_xx_xxyyy_0,  \
                             ta1_y_xx_xxyyy_1,  \
                             ta1_y_xx_xxyyz_0,  \
                             ta1_y_xx_xxyyz_1,  \
                             ta1_y_xx_xxyz_0,   \
                             ta1_y_xx_xxyz_1,   \
                             ta1_y_xx_xxyzz_0,  \
                             ta1_y_xx_xxyzz_1,  \
                             ta1_y_xx_xxzz_0,   \
                             ta1_y_xx_xxzz_1,   \
                             ta1_y_xx_xxzzz_0,  \
                             ta1_y_xx_xxzzz_1,  \
                             ta1_y_xx_xyyy_0,   \
                             ta1_y_xx_xyyy_1,   \
                             ta1_y_xx_xyyyy_0,  \
                             ta1_y_xx_xyyyy_1,  \
                             ta1_y_xx_xyyyz_0,  \
                             ta1_y_xx_xyyyz_1,  \
                             ta1_y_xx_xyyz_0,   \
                             ta1_y_xx_xyyz_1,   \
                             ta1_y_xx_xyyzz_0,  \
                             ta1_y_xx_xyyzz_1,  \
                             ta1_y_xx_xyzz_0,   \
                             ta1_y_xx_xyzz_1,   \
                             ta1_y_xx_xyzzz_0,  \
                             ta1_y_xx_xyzzz_1,  \
                             ta1_y_xx_xzzz_0,   \
                             ta1_y_xx_xzzz_1,   \
                             ta1_y_xx_xzzzz_0,  \
                             ta1_y_xx_xzzzz_1,  \
                             ta1_y_xx_yyyy_0,   \
                             ta1_y_xx_yyyy_1,   \
                             ta1_y_xx_yyyyy_0,  \
                             ta1_y_xx_yyyyy_1,  \
                             ta1_y_xx_yyyyz_0,  \
                             ta1_y_xx_yyyyz_1,  \
                             ta1_y_xx_yyyz_0,   \
                             ta1_y_xx_yyyz_1,   \
                             ta1_y_xx_yyyzz_0,  \
                             ta1_y_xx_yyyzz_1,  \
                             ta1_y_xx_yyzz_0,   \
                             ta1_y_xx_yyzz_1,   \
                             ta1_y_xx_yyzzz_0,  \
                             ta1_y_xx_yyzzz_1,  \
                             ta1_y_xx_yzzz_0,   \
                             ta1_y_xx_yzzz_1,   \
                             ta1_y_xx_yzzzz_0,  \
                             ta1_y_xx_yzzzz_1,  \
                             ta1_y_xx_zzzz_0,   \
                             ta1_y_xx_zzzz_1,   \
                             ta1_y_xx_zzzzz_0,  \
                             ta1_y_xx_zzzzz_1,  \
                             ta1_y_xxx_xxxxx_0, \
                             ta1_y_xxx_xxxxy_0, \
                             ta1_y_xxx_xxxxz_0, \
                             ta1_y_xxx_xxxyy_0, \
                             ta1_y_xxx_xxxyz_0, \
                             ta1_y_xxx_xxxzz_0, \
                             ta1_y_xxx_xxyyy_0, \
                             ta1_y_xxx_xxyyz_0, \
                             ta1_y_xxx_xxyzz_0, \
                             ta1_y_xxx_xxzzz_0, \
                             ta1_y_xxx_xyyyy_0, \
                             ta1_y_xxx_xyyyz_0, \
                             ta1_y_xxx_xyyzz_0, \
                             ta1_y_xxx_xyzzz_0, \
                             ta1_y_xxx_xzzzz_0, \
                             ta1_y_xxx_yyyyy_0, \
                             ta1_y_xxx_yyyyz_0, \
                             ta1_y_xxx_yyyzz_0, \
                             ta1_y_xxx_yyzzz_0, \
                             ta1_y_xxx_yzzzz_0, \
                             ta1_y_xxx_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxx_xxxxx_0[i] = 2.0 * ta1_y_x_xxxxx_0[i] * fe_0 - 2.0 * ta1_y_x_xxxxx_1[i] * fe_0 + 5.0 * ta1_y_xx_xxxx_0[i] * fe_0 -
                               5.0 * ta1_y_xx_xxxx_1[i] * fe_0 + ta1_y_xx_xxxxx_0[i] * pa_x[i] - ta1_y_xx_xxxxx_1[i] * pc_x[i];

        ta1_y_xxx_xxxxy_0[i] = 2.0 * ta1_y_x_xxxxy_0[i] * fe_0 - 2.0 * ta1_y_x_xxxxy_1[i] * fe_0 + 4.0 * ta1_y_xx_xxxy_0[i] * fe_0 -
                               4.0 * ta1_y_xx_xxxy_1[i] * fe_0 + ta1_y_xx_xxxxy_0[i] * pa_x[i] - ta1_y_xx_xxxxy_1[i] * pc_x[i];

        ta1_y_xxx_xxxxz_0[i] = 2.0 * ta1_y_x_xxxxz_0[i] * fe_0 - 2.0 * ta1_y_x_xxxxz_1[i] * fe_0 + 4.0 * ta1_y_xx_xxxz_0[i] * fe_0 -
                               4.0 * ta1_y_xx_xxxz_1[i] * fe_0 + ta1_y_xx_xxxxz_0[i] * pa_x[i] - ta1_y_xx_xxxxz_1[i] * pc_x[i];

        ta1_y_xxx_xxxyy_0[i] = 2.0 * ta1_y_x_xxxyy_0[i] * fe_0 - 2.0 * ta1_y_x_xxxyy_1[i] * fe_0 + 3.0 * ta1_y_xx_xxyy_0[i] * fe_0 -
                               3.0 * ta1_y_xx_xxyy_1[i] * fe_0 + ta1_y_xx_xxxyy_0[i] * pa_x[i] - ta1_y_xx_xxxyy_1[i] * pc_x[i];

        ta1_y_xxx_xxxyz_0[i] = 2.0 * ta1_y_x_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_x_xxxyz_1[i] * fe_0 + 3.0 * ta1_y_xx_xxyz_0[i] * fe_0 -
                               3.0 * ta1_y_xx_xxyz_1[i] * fe_0 + ta1_y_xx_xxxyz_0[i] * pa_x[i] - ta1_y_xx_xxxyz_1[i] * pc_x[i];

        ta1_y_xxx_xxxzz_0[i] = 2.0 * ta1_y_x_xxxzz_0[i] * fe_0 - 2.0 * ta1_y_x_xxxzz_1[i] * fe_0 + 3.0 * ta1_y_xx_xxzz_0[i] * fe_0 -
                               3.0 * ta1_y_xx_xxzz_1[i] * fe_0 + ta1_y_xx_xxxzz_0[i] * pa_x[i] - ta1_y_xx_xxxzz_1[i] * pc_x[i];

        ta1_y_xxx_xxyyy_0[i] = 2.0 * ta1_y_x_xxyyy_0[i] * fe_0 - 2.0 * ta1_y_x_xxyyy_1[i] * fe_0 + 2.0 * ta1_y_xx_xyyy_0[i] * fe_0 -
                               2.0 * ta1_y_xx_xyyy_1[i] * fe_0 + ta1_y_xx_xxyyy_0[i] * pa_x[i] - ta1_y_xx_xxyyy_1[i] * pc_x[i];

        ta1_y_xxx_xxyyz_0[i] = 2.0 * ta1_y_x_xxyyz_0[i] * fe_0 - 2.0 * ta1_y_x_xxyyz_1[i] * fe_0 + 2.0 * ta1_y_xx_xyyz_0[i] * fe_0 -
                               2.0 * ta1_y_xx_xyyz_1[i] * fe_0 + ta1_y_xx_xxyyz_0[i] * pa_x[i] - ta1_y_xx_xxyyz_1[i] * pc_x[i];

        ta1_y_xxx_xxyzz_0[i] = 2.0 * ta1_y_x_xxyzz_0[i] * fe_0 - 2.0 * ta1_y_x_xxyzz_1[i] * fe_0 + 2.0 * ta1_y_xx_xyzz_0[i] * fe_0 -
                               2.0 * ta1_y_xx_xyzz_1[i] * fe_0 + ta1_y_xx_xxyzz_0[i] * pa_x[i] - ta1_y_xx_xxyzz_1[i] * pc_x[i];

        ta1_y_xxx_xxzzz_0[i] = 2.0 * ta1_y_x_xxzzz_0[i] * fe_0 - 2.0 * ta1_y_x_xxzzz_1[i] * fe_0 + 2.0 * ta1_y_xx_xzzz_0[i] * fe_0 -
                               2.0 * ta1_y_xx_xzzz_1[i] * fe_0 + ta1_y_xx_xxzzz_0[i] * pa_x[i] - ta1_y_xx_xxzzz_1[i] * pc_x[i];

        ta1_y_xxx_xyyyy_0[i] = 2.0 * ta1_y_x_xyyyy_0[i] * fe_0 - 2.0 * ta1_y_x_xyyyy_1[i] * fe_0 + ta1_y_xx_yyyy_0[i] * fe_0 -
                               ta1_y_xx_yyyy_1[i] * fe_0 + ta1_y_xx_xyyyy_0[i] * pa_x[i] - ta1_y_xx_xyyyy_1[i] * pc_x[i];

        ta1_y_xxx_xyyyz_0[i] = 2.0 * ta1_y_x_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_x_xyyyz_1[i] * fe_0 + ta1_y_xx_yyyz_0[i] * fe_0 -
                               ta1_y_xx_yyyz_1[i] * fe_0 + ta1_y_xx_xyyyz_0[i] * pa_x[i] - ta1_y_xx_xyyyz_1[i] * pc_x[i];

        ta1_y_xxx_xyyzz_0[i] = 2.0 * ta1_y_x_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_x_xyyzz_1[i] * fe_0 + ta1_y_xx_yyzz_0[i] * fe_0 -
                               ta1_y_xx_yyzz_1[i] * fe_0 + ta1_y_xx_xyyzz_0[i] * pa_x[i] - ta1_y_xx_xyyzz_1[i] * pc_x[i];

        ta1_y_xxx_xyzzz_0[i] = 2.0 * ta1_y_x_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_x_xyzzz_1[i] * fe_0 + ta1_y_xx_yzzz_0[i] * fe_0 -
                               ta1_y_xx_yzzz_1[i] * fe_0 + ta1_y_xx_xyzzz_0[i] * pa_x[i] - ta1_y_xx_xyzzz_1[i] * pc_x[i];

        ta1_y_xxx_xzzzz_0[i] = 2.0 * ta1_y_x_xzzzz_0[i] * fe_0 - 2.0 * ta1_y_x_xzzzz_1[i] * fe_0 + ta1_y_xx_zzzz_0[i] * fe_0 -
                               ta1_y_xx_zzzz_1[i] * fe_0 + ta1_y_xx_xzzzz_0[i] * pa_x[i] - ta1_y_xx_xzzzz_1[i] * pc_x[i];

        ta1_y_xxx_yyyyy_0[i] =
            2.0 * ta1_y_x_yyyyy_0[i] * fe_0 - 2.0 * ta1_y_x_yyyyy_1[i] * fe_0 + ta1_y_xx_yyyyy_0[i] * pa_x[i] - ta1_y_xx_yyyyy_1[i] * pc_x[i];

        ta1_y_xxx_yyyyz_0[i] =
            2.0 * ta1_y_x_yyyyz_0[i] * fe_0 - 2.0 * ta1_y_x_yyyyz_1[i] * fe_0 + ta1_y_xx_yyyyz_0[i] * pa_x[i] - ta1_y_xx_yyyyz_1[i] * pc_x[i];

        ta1_y_xxx_yyyzz_0[i] =
            2.0 * ta1_y_x_yyyzz_0[i] * fe_0 - 2.0 * ta1_y_x_yyyzz_1[i] * fe_0 + ta1_y_xx_yyyzz_0[i] * pa_x[i] - ta1_y_xx_yyyzz_1[i] * pc_x[i];

        ta1_y_xxx_yyzzz_0[i] =
            2.0 * ta1_y_x_yyzzz_0[i] * fe_0 - 2.0 * ta1_y_x_yyzzz_1[i] * fe_0 + ta1_y_xx_yyzzz_0[i] * pa_x[i] - ta1_y_xx_yyzzz_1[i] * pc_x[i];

        ta1_y_xxx_yzzzz_0[i] =
            2.0 * ta1_y_x_yzzzz_0[i] * fe_0 - 2.0 * ta1_y_x_yzzzz_1[i] * fe_0 + ta1_y_xx_yzzzz_0[i] * pa_x[i] - ta1_y_xx_yzzzz_1[i] * pc_x[i];

        ta1_y_xxx_zzzzz_0[i] =
            2.0 * ta1_y_x_zzzzz_0[i] * fe_0 - 2.0 * ta1_y_x_zzzzz_1[i] * fe_0 + ta1_y_xx_zzzzz_0[i] * pa_x[i] - ta1_y_xx_zzzzz_1[i] * pc_x[i];
    }

    // Set up 231-252 components of targeted buffer : FH

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

    auto ta1_y_xxy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 251);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_y_xx_xxxx_0,   \
                             ta1_y_xx_xxxx_1,   \
                             ta1_y_xx_xxxxx_0,  \
                             ta1_y_xx_xxxxx_1,  \
                             ta1_y_xx_xxxxy_0,  \
                             ta1_y_xx_xxxxy_1,  \
                             ta1_y_xx_xxxxz_0,  \
                             ta1_y_xx_xxxxz_1,  \
                             ta1_y_xx_xxxy_0,   \
                             ta1_y_xx_xxxy_1,   \
                             ta1_y_xx_xxxyy_0,  \
                             ta1_y_xx_xxxyy_1,  \
                             ta1_y_xx_xxxyz_0,  \
                             ta1_y_xx_xxxyz_1,  \
                             ta1_y_xx_xxxz_0,   \
                             ta1_y_xx_xxxz_1,   \
                             ta1_y_xx_xxxzz_0,  \
                             ta1_y_xx_xxxzz_1,  \
                             ta1_y_xx_xxyy_0,   \
                             ta1_y_xx_xxyy_1,   \
                             ta1_y_xx_xxyyy_0,  \
                             ta1_y_xx_xxyyy_1,  \
                             ta1_y_xx_xxyyz_0,  \
                             ta1_y_xx_xxyyz_1,  \
                             ta1_y_xx_xxyz_0,   \
                             ta1_y_xx_xxyz_1,   \
                             ta1_y_xx_xxyzz_0,  \
                             ta1_y_xx_xxyzz_1,  \
                             ta1_y_xx_xxzz_0,   \
                             ta1_y_xx_xxzz_1,   \
                             ta1_y_xx_xxzzz_0,  \
                             ta1_y_xx_xxzzz_1,  \
                             ta1_y_xx_xyyy_0,   \
                             ta1_y_xx_xyyy_1,   \
                             ta1_y_xx_xyyyy_0,  \
                             ta1_y_xx_xyyyy_1,  \
                             ta1_y_xx_xyyyz_0,  \
                             ta1_y_xx_xyyyz_1,  \
                             ta1_y_xx_xyyz_0,   \
                             ta1_y_xx_xyyz_1,   \
                             ta1_y_xx_xyyzz_0,  \
                             ta1_y_xx_xyyzz_1,  \
                             ta1_y_xx_xyzz_0,   \
                             ta1_y_xx_xyzz_1,   \
                             ta1_y_xx_xyzzz_0,  \
                             ta1_y_xx_xyzzz_1,  \
                             ta1_y_xx_xzzz_0,   \
                             ta1_y_xx_xzzz_1,   \
                             ta1_y_xx_xzzzz_0,  \
                             ta1_y_xx_xzzzz_1,  \
                             ta1_y_xx_zzzzz_0,  \
                             ta1_y_xx_zzzzz_1,  \
                             ta1_y_xxy_xxxxx_0, \
                             ta1_y_xxy_xxxxy_0, \
                             ta1_y_xxy_xxxxz_0, \
                             ta1_y_xxy_xxxyy_0, \
                             ta1_y_xxy_xxxyz_0, \
                             ta1_y_xxy_xxxzz_0, \
                             ta1_y_xxy_xxyyy_0, \
                             ta1_y_xxy_xxyyz_0, \
                             ta1_y_xxy_xxyzz_0, \
                             ta1_y_xxy_xxzzz_0, \
                             ta1_y_xxy_xyyyy_0, \
                             ta1_y_xxy_xyyyz_0, \
                             ta1_y_xxy_xyyzz_0, \
                             ta1_y_xxy_xyzzz_0, \
                             ta1_y_xxy_xzzzz_0, \
                             ta1_y_xxy_yyyyy_0, \
                             ta1_y_xxy_yyyyz_0, \
                             ta1_y_xxy_yyyzz_0, \
                             ta1_y_xxy_yyzzz_0, \
                             ta1_y_xxy_yzzzz_0, \
                             ta1_y_xxy_zzzzz_0, \
                             ta1_y_xy_yyyyy_0,  \
                             ta1_y_xy_yyyyy_1,  \
                             ta1_y_xy_yyyyz_0,  \
                             ta1_y_xy_yyyyz_1,  \
                             ta1_y_xy_yyyzz_0,  \
                             ta1_y_xy_yyyzz_1,  \
                             ta1_y_xy_yyzzz_0,  \
                             ta1_y_xy_yyzzz_1,  \
                             ta1_y_xy_yzzzz_0,  \
                             ta1_y_xy_yzzzz_1,  \
                             ta1_y_y_yyyyy_0,   \
                             ta1_y_y_yyyyy_1,   \
                             ta1_y_y_yyyyz_0,   \
                             ta1_y_y_yyyyz_1,   \
                             ta1_y_y_yyyzz_0,   \
                             ta1_y_y_yyyzz_1,   \
                             ta1_y_y_yyzzz_0,   \
                             ta1_y_y_yyzzz_1,   \
                             ta1_y_y_yzzzz_0,   \
                             ta1_y_y_yzzzz_1,   \
                             ta_xx_xxxxx_1,     \
                             ta_xx_xxxxy_1,     \
                             ta_xx_xxxxz_1,     \
                             ta_xx_xxxyy_1,     \
                             ta_xx_xxxyz_1,     \
                             ta_xx_xxxzz_1,     \
                             ta_xx_xxyyy_1,     \
                             ta_xx_xxyyz_1,     \
                             ta_xx_xxyzz_1,     \
                             ta_xx_xxzzz_1,     \
                             ta_xx_xyyyy_1,     \
                             ta_xx_xyyyz_1,     \
                             ta_xx_xyyzz_1,     \
                             ta_xx_xyzzz_1,     \
                             ta_xx_xzzzz_1,     \
                             ta_xx_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxy_xxxxx_0[i] = ta_xx_xxxxx_1[i] + ta1_y_xx_xxxxx_0[i] * pa_y[i] - ta1_y_xx_xxxxx_1[i] * pc_y[i];

        ta1_y_xxy_xxxxy_0[i] =
            ta1_y_xx_xxxx_0[i] * fe_0 - ta1_y_xx_xxxx_1[i] * fe_0 + ta_xx_xxxxy_1[i] + ta1_y_xx_xxxxy_0[i] * pa_y[i] - ta1_y_xx_xxxxy_1[i] * pc_y[i];

        ta1_y_xxy_xxxxz_0[i] = ta_xx_xxxxz_1[i] + ta1_y_xx_xxxxz_0[i] * pa_y[i] - ta1_y_xx_xxxxz_1[i] * pc_y[i];

        ta1_y_xxy_xxxyy_0[i] = 2.0 * ta1_y_xx_xxxy_0[i] * fe_0 - 2.0 * ta1_y_xx_xxxy_1[i] * fe_0 + ta_xx_xxxyy_1[i] + ta1_y_xx_xxxyy_0[i] * pa_y[i] -
                               ta1_y_xx_xxxyy_1[i] * pc_y[i];

        ta1_y_xxy_xxxyz_0[i] =
            ta1_y_xx_xxxz_0[i] * fe_0 - ta1_y_xx_xxxz_1[i] * fe_0 + ta_xx_xxxyz_1[i] + ta1_y_xx_xxxyz_0[i] * pa_y[i] - ta1_y_xx_xxxyz_1[i] * pc_y[i];

        ta1_y_xxy_xxxzz_0[i] = ta_xx_xxxzz_1[i] + ta1_y_xx_xxxzz_0[i] * pa_y[i] - ta1_y_xx_xxxzz_1[i] * pc_y[i];

        ta1_y_xxy_xxyyy_0[i] = 3.0 * ta1_y_xx_xxyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxyy_1[i] * fe_0 + ta_xx_xxyyy_1[i] + ta1_y_xx_xxyyy_0[i] * pa_y[i] -
                               ta1_y_xx_xxyyy_1[i] * pc_y[i];

        ta1_y_xxy_xxyyz_0[i] = 2.0 * ta1_y_xx_xxyz_0[i] * fe_0 - 2.0 * ta1_y_xx_xxyz_1[i] * fe_0 + ta_xx_xxyyz_1[i] + ta1_y_xx_xxyyz_0[i] * pa_y[i] -
                               ta1_y_xx_xxyyz_1[i] * pc_y[i];

        ta1_y_xxy_xxyzz_0[i] =
            ta1_y_xx_xxzz_0[i] * fe_0 - ta1_y_xx_xxzz_1[i] * fe_0 + ta_xx_xxyzz_1[i] + ta1_y_xx_xxyzz_0[i] * pa_y[i] - ta1_y_xx_xxyzz_1[i] * pc_y[i];

        ta1_y_xxy_xxzzz_0[i] = ta_xx_xxzzz_1[i] + ta1_y_xx_xxzzz_0[i] * pa_y[i] - ta1_y_xx_xxzzz_1[i] * pc_y[i];

        ta1_y_xxy_xyyyy_0[i] = 4.0 * ta1_y_xx_xyyy_0[i] * fe_0 - 4.0 * ta1_y_xx_xyyy_1[i] * fe_0 + ta_xx_xyyyy_1[i] + ta1_y_xx_xyyyy_0[i] * pa_y[i] -
                               ta1_y_xx_xyyyy_1[i] * pc_y[i];

        ta1_y_xxy_xyyyz_0[i] = 3.0 * ta1_y_xx_xyyz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyyz_1[i] * fe_0 + ta_xx_xyyyz_1[i] + ta1_y_xx_xyyyz_0[i] * pa_y[i] -
                               ta1_y_xx_xyyyz_1[i] * pc_y[i];

        ta1_y_xxy_xyyzz_0[i] = 2.0 * ta1_y_xx_xyzz_0[i] * fe_0 - 2.0 * ta1_y_xx_xyzz_1[i] * fe_0 + ta_xx_xyyzz_1[i] + ta1_y_xx_xyyzz_0[i] * pa_y[i] -
                               ta1_y_xx_xyyzz_1[i] * pc_y[i];

        ta1_y_xxy_xyzzz_0[i] =
            ta1_y_xx_xzzz_0[i] * fe_0 - ta1_y_xx_xzzz_1[i] * fe_0 + ta_xx_xyzzz_1[i] + ta1_y_xx_xyzzz_0[i] * pa_y[i] - ta1_y_xx_xyzzz_1[i] * pc_y[i];

        ta1_y_xxy_xzzzz_0[i] = ta_xx_xzzzz_1[i] + ta1_y_xx_xzzzz_0[i] * pa_y[i] - ta1_y_xx_xzzzz_1[i] * pc_y[i];

        ta1_y_xxy_yyyyy_0[i] = ta1_y_y_yyyyy_0[i] * fe_0 - ta1_y_y_yyyyy_1[i] * fe_0 + ta1_y_xy_yyyyy_0[i] * pa_x[i] - ta1_y_xy_yyyyy_1[i] * pc_x[i];

        ta1_y_xxy_yyyyz_0[i] = ta1_y_y_yyyyz_0[i] * fe_0 - ta1_y_y_yyyyz_1[i] * fe_0 + ta1_y_xy_yyyyz_0[i] * pa_x[i] - ta1_y_xy_yyyyz_1[i] * pc_x[i];

        ta1_y_xxy_yyyzz_0[i] = ta1_y_y_yyyzz_0[i] * fe_0 - ta1_y_y_yyyzz_1[i] * fe_0 + ta1_y_xy_yyyzz_0[i] * pa_x[i] - ta1_y_xy_yyyzz_1[i] * pc_x[i];

        ta1_y_xxy_yyzzz_0[i] = ta1_y_y_yyzzz_0[i] * fe_0 - ta1_y_y_yyzzz_1[i] * fe_0 + ta1_y_xy_yyzzz_0[i] * pa_x[i] - ta1_y_xy_yyzzz_1[i] * pc_x[i];

        ta1_y_xxy_yzzzz_0[i] = ta1_y_y_yzzzz_0[i] * fe_0 - ta1_y_y_yzzzz_1[i] * fe_0 + ta1_y_xy_yzzzz_0[i] * pa_x[i] - ta1_y_xy_yzzzz_1[i] * pc_x[i];

        ta1_y_xxy_zzzzz_0[i] = ta_xx_zzzzz_1[i] + ta1_y_xx_zzzzz_0[i] * pa_y[i] - ta1_y_xx_zzzzz_1[i] * pc_y[i];
    }

    // Set up 252-273 components of targeted buffer : FH

    auto ta1_y_xxz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 252);

    auto ta1_y_xxz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 253);

    auto ta1_y_xxz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 254);

    auto ta1_y_xxz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 255);

    auto ta1_y_xxz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 256);

    auto ta1_y_xxz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 257);

    auto ta1_y_xxz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 258);

    auto ta1_y_xxz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 259);

    auto ta1_y_xxz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 260);

    auto ta1_y_xxz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 261);

    auto ta1_y_xxz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 262);

    auto ta1_y_xxz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 263);

    auto ta1_y_xxz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 264);

    auto ta1_y_xxz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 265);

    auto ta1_y_xxz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 266);

    auto ta1_y_xxz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 267);

    auto ta1_y_xxz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 268);

    auto ta1_y_xxz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 269);

    auto ta1_y_xxz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 270);

    auto ta1_y_xxz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 271);

    auto ta1_y_xxz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 272);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_y_xx_xxxx_0,   \
                             ta1_y_xx_xxxx_1,   \
                             ta1_y_xx_xxxxx_0,  \
                             ta1_y_xx_xxxxx_1,  \
                             ta1_y_xx_xxxxy_0,  \
                             ta1_y_xx_xxxxy_1,  \
                             ta1_y_xx_xxxxz_0,  \
                             ta1_y_xx_xxxxz_1,  \
                             ta1_y_xx_xxxy_0,   \
                             ta1_y_xx_xxxy_1,   \
                             ta1_y_xx_xxxyy_0,  \
                             ta1_y_xx_xxxyy_1,  \
                             ta1_y_xx_xxxyz_0,  \
                             ta1_y_xx_xxxyz_1,  \
                             ta1_y_xx_xxxz_0,   \
                             ta1_y_xx_xxxz_1,   \
                             ta1_y_xx_xxxzz_0,  \
                             ta1_y_xx_xxxzz_1,  \
                             ta1_y_xx_xxyy_0,   \
                             ta1_y_xx_xxyy_1,   \
                             ta1_y_xx_xxyyy_0,  \
                             ta1_y_xx_xxyyy_1,  \
                             ta1_y_xx_xxyyz_0,  \
                             ta1_y_xx_xxyyz_1,  \
                             ta1_y_xx_xxyz_0,   \
                             ta1_y_xx_xxyz_1,   \
                             ta1_y_xx_xxyzz_0,  \
                             ta1_y_xx_xxyzz_1,  \
                             ta1_y_xx_xxzz_0,   \
                             ta1_y_xx_xxzz_1,   \
                             ta1_y_xx_xxzzz_0,  \
                             ta1_y_xx_xxzzz_1,  \
                             ta1_y_xx_xyyy_0,   \
                             ta1_y_xx_xyyy_1,   \
                             ta1_y_xx_xyyyy_0,  \
                             ta1_y_xx_xyyyy_1,  \
                             ta1_y_xx_xyyyz_0,  \
                             ta1_y_xx_xyyyz_1,  \
                             ta1_y_xx_xyyz_0,   \
                             ta1_y_xx_xyyz_1,   \
                             ta1_y_xx_xyyzz_0,  \
                             ta1_y_xx_xyyzz_1,  \
                             ta1_y_xx_xyzz_0,   \
                             ta1_y_xx_xyzz_1,   \
                             ta1_y_xx_xyzzz_0,  \
                             ta1_y_xx_xyzzz_1,  \
                             ta1_y_xx_xzzz_0,   \
                             ta1_y_xx_xzzz_1,   \
                             ta1_y_xx_xzzzz_0,  \
                             ta1_y_xx_xzzzz_1,  \
                             ta1_y_xx_yyyyy_0,  \
                             ta1_y_xx_yyyyy_1,  \
                             ta1_y_xxz_xxxxx_0, \
                             ta1_y_xxz_xxxxy_0, \
                             ta1_y_xxz_xxxxz_0, \
                             ta1_y_xxz_xxxyy_0, \
                             ta1_y_xxz_xxxyz_0, \
                             ta1_y_xxz_xxxzz_0, \
                             ta1_y_xxz_xxyyy_0, \
                             ta1_y_xxz_xxyyz_0, \
                             ta1_y_xxz_xxyzz_0, \
                             ta1_y_xxz_xxzzz_0, \
                             ta1_y_xxz_xyyyy_0, \
                             ta1_y_xxz_xyyyz_0, \
                             ta1_y_xxz_xyyzz_0, \
                             ta1_y_xxz_xyzzz_0, \
                             ta1_y_xxz_xzzzz_0, \
                             ta1_y_xxz_yyyyy_0, \
                             ta1_y_xxz_yyyyz_0, \
                             ta1_y_xxz_yyyzz_0, \
                             ta1_y_xxz_yyzzz_0, \
                             ta1_y_xxz_yzzzz_0, \
                             ta1_y_xxz_zzzzz_0, \
                             ta1_y_xz_yyyyz_0,  \
                             ta1_y_xz_yyyyz_1,  \
                             ta1_y_xz_yyyzz_0,  \
                             ta1_y_xz_yyyzz_1,  \
                             ta1_y_xz_yyzzz_0,  \
                             ta1_y_xz_yyzzz_1,  \
                             ta1_y_xz_yzzzz_0,  \
                             ta1_y_xz_yzzzz_1,  \
                             ta1_y_xz_zzzzz_0,  \
                             ta1_y_xz_zzzzz_1,  \
                             ta1_y_z_yyyyz_0,   \
                             ta1_y_z_yyyyz_1,   \
                             ta1_y_z_yyyzz_0,   \
                             ta1_y_z_yyyzz_1,   \
                             ta1_y_z_yyzzz_0,   \
                             ta1_y_z_yyzzz_1,   \
                             ta1_y_z_yzzzz_0,   \
                             ta1_y_z_yzzzz_1,   \
                             ta1_y_z_zzzzz_0,   \
                             ta1_y_z_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxz_xxxxx_0[i] = ta1_y_xx_xxxxx_0[i] * pa_z[i] - ta1_y_xx_xxxxx_1[i] * pc_z[i];

        ta1_y_xxz_xxxxy_0[i] = ta1_y_xx_xxxxy_0[i] * pa_z[i] - ta1_y_xx_xxxxy_1[i] * pc_z[i];

        ta1_y_xxz_xxxxz_0[i] = ta1_y_xx_xxxx_0[i] * fe_0 - ta1_y_xx_xxxx_1[i] * fe_0 + ta1_y_xx_xxxxz_0[i] * pa_z[i] - ta1_y_xx_xxxxz_1[i] * pc_z[i];

        ta1_y_xxz_xxxyy_0[i] = ta1_y_xx_xxxyy_0[i] * pa_z[i] - ta1_y_xx_xxxyy_1[i] * pc_z[i];

        ta1_y_xxz_xxxyz_0[i] = ta1_y_xx_xxxy_0[i] * fe_0 - ta1_y_xx_xxxy_1[i] * fe_0 + ta1_y_xx_xxxyz_0[i] * pa_z[i] - ta1_y_xx_xxxyz_1[i] * pc_z[i];

        ta1_y_xxz_xxxzz_0[i] =
            2.0 * ta1_y_xx_xxxz_0[i] * fe_0 - 2.0 * ta1_y_xx_xxxz_1[i] * fe_0 + ta1_y_xx_xxxzz_0[i] * pa_z[i] - ta1_y_xx_xxxzz_1[i] * pc_z[i];

        ta1_y_xxz_xxyyy_0[i] = ta1_y_xx_xxyyy_0[i] * pa_z[i] - ta1_y_xx_xxyyy_1[i] * pc_z[i];

        ta1_y_xxz_xxyyz_0[i] = ta1_y_xx_xxyy_0[i] * fe_0 - ta1_y_xx_xxyy_1[i] * fe_0 + ta1_y_xx_xxyyz_0[i] * pa_z[i] - ta1_y_xx_xxyyz_1[i] * pc_z[i];

        ta1_y_xxz_xxyzz_0[i] =
            2.0 * ta1_y_xx_xxyz_0[i] * fe_0 - 2.0 * ta1_y_xx_xxyz_1[i] * fe_0 + ta1_y_xx_xxyzz_0[i] * pa_z[i] - ta1_y_xx_xxyzz_1[i] * pc_z[i];

        ta1_y_xxz_xxzzz_0[i] =
            3.0 * ta1_y_xx_xxzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxzz_1[i] * fe_0 + ta1_y_xx_xxzzz_0[i] * pa_z[i] - ta1_y_xx_xxzzz_1[i] * pc_z[i];

        ta1_y_xxz_xyyyy_0[i] = ta1_y_xx_xyyyy_0[i] * pa_z[i] - ta1_y_xx_xyyyy_1[i] * pc_z[i];

        ta1_y_xxz_xyyyz_0[i] = ta1_y_xx_xyyy_0[i] * fe_0 - ta1_y_xx_xyyy_1[i] * fe_0 + ta1_y_xx_xyyyz_0[i] * pa_z[i] - ta1_y_xx_xyyyz_1[i] * pc_z[i];

        ta1_y_xxz_xyyzz_0[i] =
            2.0 * ta1_y_xx_xyyz_0[i] * fe_0 - 2.0 * ta1_y_xx_xyyz_1[i] * fe_0 + ta1_y_xx_xyyzz_0[i] * pa_z[i] - ta1_y_xx_xyyzz_1[i] * pc_z[i];

        ta1_y_xxz_xyzzz_0[i] =
            3.0 * ta1_y_xx_xyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyzz_1[i] * fe_0 + ta1_y_xx_xyzzz_0[i] * pa_z[i] - ta1_y_xx_xyzzz_1[i] * pc_z[i];

        ta1_y_xxz_xzzzz_0[i] =
            4.0 * ta1_y_xx_xzzz_0[i] * fe_0 - 4.0 * ta1_y_xx_xzzz_1[i] * fe_0 + ta1_y_xx_xzzzz_0[i] * pa_z[i] - ta1_y_xx_xzzzz_1[i] * pc_z[i];

        ta1_y_xxz_yyyyy_0[i] = ta1_y_xx_yyyyy_0[i] * pa_z[i] - ta1_y_xx_yyyyy_1[i] * pc_z[i];

        ta1_y_xxz_yyyyz_0[i] = ta1_y_z_yyyyz_0[i] * fe_0 - ta1_y_z_yyyyz_1[i] * fe_0 + ta1_y_xz_yyyyz_0[i] * pa_x[i] - ta1_y_xz_yyyyz_1[i] * pc_x[i];

        ta1_y_xxz_yyyzz_0[i] = ta1_y_z_yyyzz_0[i] * fe_0 - ta1_y_z_yyyzz_1[i] * fe_0 + ta1_y_xz_yyyzz_0[i] * pa_x[i] - ta1_y_xz_yyyzz_1[i] * pc_x[i];

        ta1_y_xxz_yyzzz_0[i] = ta1_y_z_yyzzz_0[i] * fe_0 - ta1_y_z_yyzzz_1[i] * fe_0 + ta1_y_xz_yyzzz_0[i] * pa_x[i] - ta1_y_xz_yyzzz_1[i] * pc_x[i];

        ta1_y_xxz_yzzzz_0[i] = ta1_y_z_yzzzz_0[i] * fe_0 - ta1_y_z_yzzzz_1[i] * fe_0 + ta1_y_xz_yzzzz_0[i] * pa_x[i] - ta1_y_xz_yzzzz_1[i] * pc_x[i];

        ta1_y_xxz_zzzzz_0[i] = ta1_y_z_zzzzz_0[i] * fe_0 - ta1_y_z_zzzzz_1[i] * fe_0 + ta1_y_xz_zzzzz_0[i] * pa_x[i] - ta1_y_xz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 273-294 components of targeted buffer : FH

    auto ta1_y_xyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 273);

    auto ta1_y_xyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 274);

    auto ta1_y_xyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 275);

    auto ta1_y_xyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 276);

    auto ta1_y_xyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 277);

    auto ta1_y_xyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 278);

    auto ta1_y_xyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 279);

    auto ta1_y_xyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 280);

    auto ta1_y_xyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 281);

    auto ta1_y_xyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 282);

    auto ta1_y_xyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 283);

    auto ta1_y_xyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 284);

    auto ta1_y_xyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 285);

    auto ta1_y_xyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 286);

    auto ta1_y_xyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 287);

    auto ta1_y_xyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 288);

    auto ta1_y_xyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 289);

    auto ta1_y_xyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 290);

    auto ta1_y_xyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 291);

    auto ta1_y_xyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 292);

    auto ta1_y_xyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 293);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_y_xyy_xxxxx_0, \
                             ta1_y_xyy_xxxxy_0, \
                             ta1_y_xyy_xxxxz_0, \
                             ta1_y_xyy_xxxyy_0, \
                             ta1_y_xyy_xxxyz_0, \
                             ta1_y_xyy_xxxzz_0, \
                             ta1_y_xyy_xxyyy_0, \
                             ta1_y_xyy_xxyyz_0, \
                             ta1_y_xyy_xxyzz_0, \
                             ta1_y_xyy_xxzzz_0, \
                             ta1_y_xyy_xyyyy_0, \
                             ta1_y_xyy_xyyyz_0, \
                             ta1_y_xyy_xyyzz_0, \
                             ta1_y_xyy_xyzzz_0, \
                             ta1_y_xyy_xzzzz_0, \
                             ta1_y_xyy_yyyyy_0, \
                             ta1_y_xyy_yyyyz_0, \
                             ta1_y_xyy_yyyzz_0, \
                             ta1_y_xyy_yyzzz_0, \
                             ta1_y_xyy_yzzzz_0, \
                             ta1_y_xyy_zzzzz_0, \
                             ta1_y_yy_xxxx_0,   \
                             ta1_y_yy_xxxx_1,   \
                             ta1_y_yy_xxxxx_0,  \
                             ta1_y_yy_xxxxx_1,  \
                             ta1_y_yy_xxxxy_0,  \
                             ta1_y_yy_xxxxy_1,  \
                             ta1_y_yy_xxxxz_0,  \
                             ta1_y_yy_xxxxz_1,  \
                             ta1_y_yy_xxxy_0,   \
                             ta1_y_yy_xxxy_1,   \
                             ta1_y_yy_xxxyy_0,  \
                             ta1_y_yy_xxxyy_1,  \
                             ta1_y_yy_xxxyz_0,  \
                             ta1_y_yy_xxxyz_1,  \
                             ta1_y_yy_xxxz_0,   \
                             ta1_y_yy_xxxz_1,   \
                             ta1_y_yy_xxxzz_0,  \
                             ta1_y_yy_xxxzz_1,  \
                             ta1_y_yy_xxyy_0,   \
                             ta1_y_yy_xxyy_1,   \
                             ta1_y_yy_xxyyy_0,  \
                             ta1_y_yy_xxyyy_1,  \
                             ta1_y_yy_xxyyz_0,  \
                             ta1_y_yy_xxyyz_1,  \
                             ta1_y_yy_xxyz_0,   \
                             ta1_y_yy_xxyz_1,   \
                             ta1_y_yy_xxyzz_0,  \
                             ta1_y_yy_xxyzz_1,  \
                             ta1_y_yy_xxzz_0,   \
                             ta1_y_yy_xxzz_1,   \
                             ta1_y_yy_xxzzz_0,  \
                             ta1_y_yy_xxzzz_1,  \
                             ta1_y_yy_xyyy_0,   \
                             ta1_y_yy_xyyy_1,   \
                             ta1_y_yy_xyyyy_0,  \
                             ta1_y_yy_xyyyy_1,  \
                             ta1_y_yy_xyyyz_0,  \
                             ta1_y_yy_xyyyz_1,  \
                             ta1_y_yy_xyyz_0,   \
                             ta1_y_yy_xyyz_1,   \
                             ta1_y_yy_xyyzz_0,  \
                             ta1_y_yy_xyyzz_1,  \
                             ta1_y_yy_xyzz_0,   \
                             ta1_y_yy_xyzz_1,   \
                             ta1_y_yy_xyzzz_0,  \
                             ta1_y_yy_xyzzz_1,  \
                             ta1_y_yy_xzzz_0,   \
                             ta1_y_yy_xzzz_1,   \
                             ta1_y_yy_xzzzz_0,  \
                             ta1_y_yy_xzzzz_1,  \
                             ta1_y_yy_yyyy_0,   \
                             ta1_y_yy_yyyy_1,   \
                             ta1_y_yy_yyyyy_0,  \
                             ta1_y_yy_yyyyy_1,  \
                             ta1_y_yy_yyyyz_0,  \
                             ta1_y_yy_yyyyz_1,  \
                             ta1_y_yy_yyyz_0,   \
                             ta1_y_yy_yyyz_1,   \
                             ta1_y_yy_yyyzz_0,  \
                             ta1_y_yy_yyyzz_1,  \
                             ta1_y_yy_yyzz_0,   \
                             ta1_y_yy_yyzz_1,   \
                             ta1_y_yy_yyzzz_0,  \
                             ta1_y_yy_yyzzz_1,  \
                             ta1_y_yy_yzzz_0,   \
                             ta1_y_yy_yzzz_1,   \
                             ta1_y_yy_yzzzz_0,  \
                             ta1_y_yy_yzzzz_1,  \
                             ta1_y_yy_zzzz_0,   \
                             ta1_y_yy_zzzz_1,   \
                             ta1_y_yy_zzzzz_0,  \
                             ta1_y_yy_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyy_xxxxx_0[i] =
            5.0 * ta1_y_yy_xxxx_0[i] * fe_0 - 5.0 * ta1_y_yy_xxxx_1[i] * fe_0 + ta1_y_yy_xxxxx_0[i] * pa_x[i] - ta1_y_yy_xxxxx_1[i] * pc_x[i];

        ta1_y_xyy_xxxxy_0[i] =
            4.0 * ta1_y_yy_xxxy_0[i] * fe_0 - 4.0 * ta1_y_yy_xxxy_1[i] * fe_0 + ta1_y_yy_xxxxy_0[i] * pa_x[i] - ta1_y_yy_xxxxy_1[i] * pc_x[i];

        ta1_y_xyy_xxxxz_0[i] =
            4.0 * ta1_y_yy_xxxz_0[i] * fe_0 - 4.0 * ta1_y_yy_xxxz_1[i] * fe_0 + ta1_y_yy_xxxxz_0[i] * pa_x[i] - ta1_y_yy_xxxxz_1[i] * pc_x[i];

        ta1_y_xyy_xxxyy_0[i] =
            3.0 * ta1_y_yy_xxyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyy_1[i] * fe_0 + ta1_y_yy_xxxyy_0[i] * pa_x[i] - ta1_y_yy_xxxyy_1[i] * pc_x[i];

        ta1_y_xyy_xxxyz_0[i] =
            3.0 * ta1_y_yy_xxyz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyz_1[i] * fe_0 + ta1_y_yy_xxxyz_0[i] * pa_x[i] - ta1_y_yy_xxxyz_1[i] * pc_x[i];

        ta1_y_xyy_xxxzz_0[i] =
            3.0 * ta1_y_yy_xxzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxzz_1[i] * fe_0 + ta1_y_yy_xxxzz_0[i] * pa_x[i] - ta1_y_yy_xxxzz_1[i] * pc_x[i];

        ta1_y_xyy_xxyyy_0[i] =
            2.0 * ta1_y_yy_xyyy_0[i] * fe_0 - 2.0 * ta1_y_yy_xyyy_1[i] * fe_0 + ta1_y_yy_xxyyy_0[i] * pa_x[i] - ta1_y_yy_xxyyy_1[i] * pc_x[i];

        ta1_y_xyy_xxyyz_0[i] =
            2.0 * ta1_y_yy_xyyz_0[i] * fe_0 - 2.0 * ta1_y_yy_xyyz_1[i] * fe_0 + ta1_y_yy_xxyyz_0[i] * pa_x[i] - ta1_y_yy_xxyyz_1[i] * pc_x[i];

        ta1_y_xyy_xxyzz_0[i] =
            2.0 * ta1_y_yy_xyzz_0[i] * fe_0 - 2.0 * ta1_y_yy_xyzz_1[i] * fe_0 + ta1_y_yy_xxyzz_0[i] * pa_x[i] - ta1_y_yy_xxyzz_1[i] * pc_x[i];

        ta1_y_xyy_xxzzz_0[i] =
            2.0 * ta1_y_yy_xzzz_0[i] * fe_0 - 2.0 * ta1_y_yy_xzzz_1[i] * fe_0 + ta1_y_yy_xxzzz_0[i] * pa_x[i] - ta1_y_yy_xxzzz_1[i] * pc_x[i];

        ta1_y_xyy_xyyyy_0[i] = ta1_y_yy_yyyy_0[i] * fe_0 - ta1_y_yy_yyyy_1[i] * fe_0 + ta1_y_yy_xyyyy_0[i] * pa_x[i] - ta1_y_yy_xyyyy_1[i] * pc_x[i];

        ta1_y_xyy_xyyyz_0[i] = ta1_y_yy_yyyz_0[i] * fe_0 - ta1_y_yy_yyyz_1[i] * fe_0 + ta1_y_yy_xyyyz_0[i] * pa_x[i] - ta1_y_yy_xyyyz_1[i] * pc_x[i];

        ta1_y_xyy_xyyzz_0[i] = ta1_y_yy_yyzz_0[i] * fe_0 - ta1_y_yy_yyzz_1[i] * fe_0 + ta1_y_yy_xyyzz_0[i] * pa_x[i] - ta1_y_yy_xyyzz_1[i] * pc_x[i];

        ta1_y_xyy_xyzzz_0[i] = ta1_y_yy_yzzz_0[i] * fe_0 - ta1_y_yy_yzzz_1[i] * fe_0 + ta1_y_yy_xyzzz_0[i] * pa_x[i] - ta1_y_yy_xyzzz_1[i] * pc_x[i];

        ta1_y_xyy_xzzzz_0[i] = ta1_y_yy_zzzz_0[i] * fe_0 - ta1_y_yy_zzzz_1[i] * fe_0 + ta1_y_yy_xzzzz_0[i] * pa_x[i] - ta1_y_yy_xzzzz_1[i] * pc_x[i];

        ta1_y_xyy_yyyyy_0[i] = ta1_y_yy_yyyyy_0[i] * pa_x[i] - ta1_y_yy_yyyyy_1[i] * pc_x[i];

        ta1_y_xyy_yyyyz_0[i] = ta1_y_yy_yyyyz_0[i] * pa_x[i] - ta1_y_yy_yyyyz_1[i] * pc_x[i];

        ta1_y_xyy_yyyzz_0[i] = ta1_y_yy_yyyzz_0[i] * pa_x[i] - ta1_y_yy_yyyzz_1[i] * pc_x[i];

        ta1_y_xyy_yyzzz_0[i] = ta1_y_yy_yyzzz_0[i] * pa_x[i] - ta1_y_yy_yyzzz_1[i] * pc_x[i];

        ta1_y_xyy_yzzzz_0[i] = ta1_y_yy_yzzzz_0[i] * pa_x[i] - ta1_y_yy_yzzzz_1[i] * pc_x[i];

        ta1_y_xyy_zzzzz_0[i] = ta1_y_yy_zzzzz_0[i] * pa_x[i] - ta1_y_yy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 294-315 components of targeted buffer : FH

    auto ta1_y_xyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 294);

    auto ta1_y_xyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 295);

    auto ta1_y_xyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 296);

    auto ta1_y_xyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 297);

    auto ta1_y_xyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 298);

    auto ta1_y_xyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 299);

    auto ta1_y_xyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 300);

    auto ta1_y_xyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 301);

    auto ta1_y_xyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 302);

    auto ta1_y_xyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 303);

    auto ta1_y_xyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 304);

    auto ta1_y_xyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 305);

    auto ta1_y_xyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 306);

    auto ta1_y_xyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 307);

    auto ta1_y_xyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 308);

    auto ta1_y_xyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 309);

    auto ta1_y_xyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 310);

    auto ta1_y_xyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 311);

    auto ta1_y_xyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 312);

    auto ta1_y_xyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 313);

    auto ta1_y_xyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 314);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_xy_xxxxx_0,  \
                             ta1_y_xy_xxxxx_1,  \
                             ta1_y_xy_xxxxy_0,  \
                             ta1_y_xy_xxxxy_1,  \
                             ta1_y_xy_xxxyy_0,  \
                             ta1_y_xy_xxxyy_1,  \
                             ta1_y_xy_xxyyy_0,  \
                             ta1_y_xy_xxyyy_1,  \
                             ta1_y_xy_xyyyy_0,  \
                             ta1_y_xy_xyyyy_1,  \
                             ta1_y_xyz_xxxxx_0, \
                             ta1_y_xyz_xxxxy_0, \
                             ta1_y_xyz_xxxxz_0, \
                             ta1_y_xyz_xxxyy_0, \
                             ta1_y_xyz_xxxyz_0, \
                             ta1_y_xyz_xxxzz_0, \
                             ta1_y_xyz_xxyyy_0, \
                             ta1_y_xyz_xxyyz_0, \
                             ta1_y_xyz_xxyzz_0, \
                             ta1_y_xyz_xxzzz_0, \
                             ta1_y_xyz_xyyyy_0, \
                             ta1_y_xyz_xyyyz_0, \
                             ta1_y_xyz_xyyzz_0, \
                             ta1_y_xyz_xyzzz_0, \
                             ta1_y_xyz_xzzzz_0, \
                             ta1_y_xyz_yyyyy_0, \
                             ta1_y_xyz_yyyyz_0, \
                             ta1_y_xyz_yyyzz_0, \
                             ta1_y_xyz_yyzzz_0, \
                             ta1_y_xyz_yzzzz_0, \
                             ta1_y_xyz_zzzzz_0, \
                             ta1_y_xz_xxxxz_0,  \
                             ta1_y_xz_xxxxz_1,  \
                             ta1_y_xz_xxxzz_0,  \
                             ta1_y_xz_xxxzz_1,  \
                             ta1_y_xz_xxzzz_0,  \
                             ta1_y_xz_xxzzz_1,  \
                             ta1_y_xz_xzzzz_0,  \
                             ta1_y_xz_xzzzz_1,  \
                             ta1_y_yz_xxxyz_0,  \
                             ta1_y_yz_xxxyz_1,  \
                             ta1_y_yz_xxyyz_0,  \
                             ta1_y_yz_xxyyz_1,  \
                             ta1_y_yz_xxyz_0,   \
                             ta1_y_yz_xxyz_1,   \
                             ta1_y_yz_xxyzz_0,  \
                             ta1_y_yz_xxyzz_1,  \
                             ta1_y_yz_xyyyz_0,  \
                             ta1_y_yz_xyyyz_1,  \
                             ta1_y_yz_xyyz_0,   \
                             ta1_y_yz_xyyz_1,   \
                             ta1_y_yz_xyyzz_0,  \
                             ta1_y_yz_xyyzz_1,  \
                             ta1_y_yz_xyzz_0,   \
                             ta1_y_yz_xyzz_1,   \
                             ta1_y_yz_xyzzz_0,  \
                             ta1_y_yz_xyzzz_1,  \
                             ta1_y_yz_yyyyy_0,  \
                             ta1_y_yz_yyyyy_1,  \
                             ta1_y_yz_yyyyz_0,  \
                             ta1_y_yz_yyyyz_1,  \
                             ta1_y_yz_yyyz_0,   \
                             ta1_y_yz_yyyz_1,   \
                             ta1_y_yz_yyyzz_0,  \
                             ta1_y_yz_yyyzz_1,  \
                             ta1_y_yz_yyzz_0,   \
                             ta1_y_yz_yyzz_1,   \
                             ta1_y_yz_yyzzz_0,  \
                             ta1_y_yz_yyzzz_1,  \
                             ta1_y_yz_yzzz_0,   \
                             ta1_y_yz_yzzz_1,   \
                             ta1_y_yz_yzzzz_0,  \
                             ta1_y_yz_yzzzz_1,  \
                             ta1_y_yz_zzzzz_0,  \
                             ta1_y_yz_zzzzz_1,  \
                             ta_xz_xxxxz_1,     \
                             ta_xz_xxxzz_1,     \
                             ta_xz_xxzzz_1,     \
                             ta_xz_xzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyz_xxxxx_0[i] = ta1_y_xy_xxxxx_0[i] * pa_z[i] - ta1_y_xy_xxxxx_1[i] * pc_z[i];

        ta1_y_xyz_xxxxy_0[i] = ta1_y_xy_xxxxy_0[i] * pa_z[i] - ta1_y_xy_xxxxy_1[i] * pc_z[i];

        ta1_y_xyz_xxxxz_0[i] = ta_xz_xxxxz_1[i] + ta1_y_xz_xxxxz_0[i] * pa_y[i] - ta1_y_xz_xxxxz_1[i] * pc_y[i];

        ta1_y_xyz_xxxyy_0[i] = ta1_y_xy_xxxyy_0[i] * pa_z[i] - ta1_y_xy_xxxyy_1[i] * pc_z[i];

        ta1_y_xyz_xxxyz_0[i] =
            3.0 * ta1_y_yz_xxyz_0[i] * fe_0 - 3.0 * ta1_y_yz_xxyz_1[i] * fe_0 + ta1_y_yz_xxxyz_0[i] * pa_x[i] - ta1_y_yz_xxxyz_1[i] * pc_x[i];

        ta1_y_xyz_xxxzz_0[i] = ta_xz_xxxzz_1[i] + ta1_y_xz_xxxzz_0[i] * pa_y[i] - ta1_y_xz_xxxzz_1[i] * pc_y[i];

        ta1_y_xyz_xxyyy_0[i] = ta1_y_xy_xxyyy_0[i] * pa_z[i] - ta1_y_xy_xxyyy_1[i] * pc_z[i];

        ta1_y_xyz_xxyyz_0[i] =
            2.0 * ta1_y_yz_xyyz_0[i] * fe_0 - 2.0 * ta1_y_yz_xyyz_1[i] * fe_0 + ta1_y_yz_xxyyz_0[i] * pa_x[i] - ta1_y_yz_xxyyz_1[i] * pc_x[i];

        ta1_y_xyz_xxyzz_0[i] =
            2.0 * ta1_y_yz_xyzz_0[i] * fe_0 - 2.0 * ta1_y_yz_xyzz_1[i] * fe_0 + ta1_y_yz_xxyzz_0[i] * pa_x[i] - ta1_y_yz_xxyzz_1[i] * pc_x[i];

        ta1_y_xyz_xxzzz_0[i] = ta_xz_xxzzz_1[i] + ta1_y_xz_xxzzz_0[i] * pa_y[i] - ta1_y_xz_xxzzz_1[i] * pc_y[i];

        ta1_y_xyz_xyyyy_0[i] = ta1_y_xy_xyyyy_0[i] * pa_z[i] - ta1_y_xy_xyyyy_1[i] * pc_z[i];

        ta1_y_xyz_xyyyz_0[i] = ta1_y_yz_yyyz_0[i] * fe_0 - ta1_y_yz_yyyz_1[i] * fe_0 + ta1_y_yz_xyyyz_0[i] * pa_x[i] - ta1_y_yz_xyyyz_1[i] * pc_x[i];

        ta1_y_xyz_xyyzz_0[i] = ta1_y_yz_yyzz_0[i] * fe_0 - ta1_y_yz_yyzz_1[i] * fe_0 + ta1_y_yz_xyyzz_0[i] * pa_x[i] - ta1_y_yz_xyyzz_1[i] * pc_x[i];

        ta1_y_xyz_xyzzz_0[i] = ta1_y_yz_yzzz_0[i] * fe_0 - ta1_y_yz_yzzz_1[i] * fe_0 + ta1_y_yz_xyzzz_0[i] * pa_x[i] - ta1_y_yz_xyzzz_1[i] * pc_x[i];

        ta1_y_xyz_xzzzz_0[i] = ta_xz_xzzzz_1[i] + ta1_y_xz_xzzzz_0[i] * pa_y[i] - ta1_y_xz_xzzzz_1[i] * pc_y[i];

        ta1_y_xyz_yyyyy_0[i] = ta1_y_yz_yyyyy_0[i] * pa_x[i] - ta1_y_yz_yyyyy_1[i] * pc_x[i];

        ta1_y_xyz_yyyyz_0[i] = ta1_y_yz_yyyyz_0[i] * pa_x[i] - ta1_y_yz_yyyyz_1[i] * pc_x[i];

        ta1_y_xyz_yyyzz_0[i] = ta1_y_yz_yyyzz_0[i] * pa_x[i] - ta1_y_yz_yyyzz_1[i] * pc_x[i];

        ta1_y_xyz_yyzzz_0[i] = ta1_y_yz_yyzzz_0[i] * pa_x[i] - ta1_y_yz_yyzzz_1[i] * pc_x[i];

        ta1_y_xyz_yzzzz_0[i] = ta1_y_yz_yzzzz_0[i] * pa_x[i] - ta1_y_yz_yzzzz_1[i] * pc_x[i];

        ta1_y_xyz_zzzzz_0[i] = ta1_y_yz_zzzzz_0[i] * pa_x[i] - ta1_y_yz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 315-336 components of targeted buffer : FH

    auto ta1_y_xzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 315);

    auto ta1_y_xzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 316);

    auto ta1_y_xzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 317);

    auto ta1_y_xzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 318);

    auto ta1_y_xzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 319);

    auto ta1_y_xzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 320);

    auto ta1_y_xzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 321);

    auto ta1_y_xzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 322);

    auto ta1_y_xzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 323);

    auto ta1_y_xzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 324);

    auto ta1_y_xzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 325);

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

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_y_xzz_xxxxx_0, \
                             ta1_y_xzz_xxxxy_0, \
                             ta1_y_xzz_xxxxz_0, \
                             ta1_y_xzz_xxxyy_0, \
                             ta1_y_xzz_xxxyz_0, \
                             ta1_y_xzz_xxxzz_0, \
                             ta1_y_xzz_xxyyy_0, \
                             ta1_y_xzz_xxyyz_0, \
                             ta1_y_xzz_xxyzz_0, \
                             ta1_y_xzz_xxzzz_0, \
                             ta1_y_xzz_xyyyy_0, \
                             ta1_y_xzz_xyyyz_0, \
                             ta1_y_xzz_xyyzz_0, \
                             ta1_y_xzz_xyzzz_0, \
                             ta1_y_xzz_xzzzz_0, \
                             ta1_y_xzz_yyyyy_0, \
                             ta1_y_xzz_yyyyz_0, \
                             ta1_y_xzz_yyyzz_0, \
                             ta1_y_xzz_yyzzz_0, \
                             ta1_y_xzz_yzzzz_0, \
                             ta1_y_xzz_zzzzz_0, \
                             ta1_y_zz_xxxx_0,   \
                             ta1_y_zz_xxxx_1,   \
                             ta1_y_zz_xxxxx_0,  \
                             ta1_y_zz_xxxxx_1,  \
                             ta1_y_zz_xxxxy_0,  \
                             ta1_y_zz_xxxxy_1,  \
                             ta1_y_zz_xxxxz_0,  \
                             ta1_y_zz_xxxxz_1,  \
                             ta1_y_zz_xxxy_0,   \
                             ta1_y_zz_xxxy_1,   \
                             ta1_y_zz_xxxyy_0,  \
                             ta1_y_zz_xxxyy_1,  \
                             ta1_y_zz_xxxyz_0,  \
                             ta1_y_zz_xxxyz_1,  \
                             ta1_y_zz_xxxz_0,   \
                             ta1_y_zz_xxxz_1,   \
                             ta1_y_zz_xxxzz_0,  \
                             ta1_y_zz_xxxzz_1,  \
                             ta1_y_zz_xxyy_0,   \
                             ta1_y_zz_xxyy_1,   \
                             ta1_y_zz_xxyyy_0,  \
                             ta1_y_zz_xxyyy_1,  \
                             ta1_y_zz_xxyyz_0,  \
                             ta1_y_zz_xxyyz_1,  \
                             ta1_y_zz_xxyz_0,   \
                             ta1_y_zz_xxyz_1,   \
                             ta1_y_zz_xxyzz_0,  \
                             ta1_y_zz_xxyzz_1,  \
                             ta1_y_zz_xxzz_0,   \
                             ta1_y_zz_xxzz_1,   \
                             ta1_y_zz_xxzzz_0,  \
                             ta1_y_zz_xxzzz_1,  \
                             ta1_y_zz_xyyy_0,   \
                             ta1_y_zz_xyyy_1,   \
                             ta1_y_zz_xyyyy_0,  \
                             ta1_y_zz_xyyyy_1,  \
                             ta1_y_zz_xyyyz_0,  \
                             ta1_y_zz_xyyyz_1,  \
                             ta1_y_zz_xyyz_0,   \
                             ta1_y_zz_xyyz_1,   \
                             ta1_y_zz_xyyzz_0,  \
                             ta1_y_zz_xyyzz_1,  \
                             ta1_y_zz_xyzz_0,   \
                             ta1_y_zz_xyzz_1,   \
                             ta1_y_zz_xyzzz_0,  \
                             ta1_y_zz_xyzzz_1,  \
                             ta1_y_zz_xzzz_0,   \
                             ta1_y_zz_xzzz_1,   \
                             ta1_y_zz_xzzzz_0,  \
                             ta1_y_zz_xzzzz_1,  \
                             ta1_y_zz_yyyy_0,   \
                             ta1_y_zz_yyyy_1,   \
                             ta1_y_zz_yyyyy_0,  \
                             ta1_y_zz_yyyyy_1,  \
                             ta1_y_zz_yyyyz_0,  \
                             ta1_y_zz_yyyyz_1,  \
                             ta1_y_zz_yyyz_0,   \
                             ta1_y_zz_yyyz_1,   \
                             ta1_y_zz_yyyzz_0,  \
                             ta1_y_zz_yyyzz_1,  \
                             ta1_y_zz_yyzz_0,   \
                             ta1_y_zz_yyzz_1,   \
                             ta1_y_zz_yyzzz_0,  \
                             ta1_y_zz_yyzzz_1,  \
                             ta1_y_zz_yzzz_0,   \
                             ta1_y_zz_yzzz_1,   \
                             ta1_y_zz_yzzzz_0,  \
                             ta1_y_zz_yzzzz_1,  \
                             ta1_y_zz_zzzz_0,   \
                             ta1_y_zz_zzzz_1,   \
                             ta1_y_zz_zzzzz_0,  \
                             ta1_y_zz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzz_xxxxx_0[i] =
            5.0 * ta1_y_zz_xxxx_0[i] * fe_0 - 5.0 * ta1_y_zz_xxxx_1[i] * fe_0 + ta1_y_zz_xxxxx_0[i] * pa_x[i] - ta1_y_zz_xxxxx_1[i] * pc_x[i];

        ta1_y_xzz_xxxxy_0[i] =
            4.0 * ta1_y_zz_xxxy_0[i] * fe_0 - 4.0 * ta1_y_zz_xxxy_1[i] * fe_0 + ta1_y_zz_xxxxy_0[i] * pa_x[i] - ta1_y_zz_xxxxy_1[i] * pc_x[i];

        ta1_y_xzz_xxxxz_0[i] =
            4.0 * ta1_y_zz_xxxz_0[i] * fe_0 - 4.0 * ta1_y_zz_xxxz_1[i] * fe_0 + ta1_y_zz_xxxxz_0[i] * pa_x[i] - ta1_y_zz_xxxxz_1[i] * pc_x[i];

        ta1_y_xzz_xxxyy_0[i] =
            3.0 * ta1_y_zz_xxyy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyy_1[i] * fe_0 + ta1_y_zz_xxxyy_0[i] * pa_x[i] - ta1_y_zz_xxxyy_1[i] * pc_x[i];

        ta1_y_xzz_xxxyz_0[i] =
            3.0 * ta1_y_zz_xxyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyz_1[i] * fe_0 + ta1_y_zz_xxxyz_0[i] * pa_x[i] - ta1_y_zz_xxxyz_1[i] * pc_x[i];

        ta1_y_xzz_xxxzz_0[i] =
            3.0 * ta1_y_zz_xxzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxzz_1[i] * fe_0 + ta1_y_zz_xxxzz_0[i] * pa_x[i] - ta1_y_zz_xxxzz_1[i] * pc_x[i];

        ta1_y_xzz_xxyyy_0[i] =
            2.0 * ta1_y_zz_xyyy_0[i] * fe_0 - 2.0 * ta1_y_zz_xyyy_1[i] * fe_0 + ta1_y_zz_xxyyy_0[i] * pa_x[i] - ta1_y_zz_xxyyy_1[i] * pc_x[i];

        ta1_y_xzz_xxyyz_0[i] =
            2.0 * ta1_y_zz_xyyz_0[i] * fe_0 - 2.0 * ta1_y_zz_xyyz_1[i] * fe_0 + ta1_y_zz_xxyyz_0[i] * pa_x[i] - ta1_y_zz_xxyyz_1[i] * pc_x[i];

        ta1_y_xzz_xxyzz_0[i] =
            2.0 * ta1_y_zz_xyzz_0[i] * fe_0 - 2.0 * ta1_y_zz_xyzz_1[i] * fe_0 + ta1_y_zz_xxyzz_0[i] * pa_x[i] - ta1_y_zz_xxyzz_1[i] * pc_x[i];

        ta1_y_xzz_xxzzz_0[i] =
            2.0 * ta1_y_zz_xzzz_0[i] * fe_0 - 2.0 * ta1_y_zz_xzzz_1[i] * fe_0 + ta1_y_zz_xxzzz_0[i] * pa_x[i] - ta1_y_zz_xxzzz_1[i] * pc_x[i];

        ta1_y_xzz_xyyyy_0[i] = ta1_y_zz_yyyy_0[i] * fe_0 - ta1_y_zz_yyyy_1[i] * fe_0 + ta1_y_zz_xyyyy_0[i] * pa_x[i] - ta1_y_zz_xyyyy_1[i] * pc_x[i];

        ta1_y_xzz_xyyyz_0[i] = ta1_y_zz_yyyz_0[i] * fe_0 - ta1_y_zz_yyyz_1[i] * fe_0 + ta1_y_zz_xyyyz_0[i] * pa_x[i] - ta1_y_zz_xyyyz_1[i] * pc_x[i];

        ta1_y_xzz_xyyzz_0[i] = ta1_y_zz_yyzz_0[i] * fe_0 - ta1_y_zz_yyzz_1[i] * fe_0 + ta1_y_zz_xyyzz_0[i] * pa_x[i] - ta1_y_zz_xyyzz_1[i] * pc_x[i];

        ta1_y_xzz_xyzzz_0[i] = ta1_y_zz_yzzz_0[i] * fe_0 - ta1_y_zz_yzzz_1[i] * fe_0 + ta1_y_zz_xyzzz_0[i] * pa_x[i] - ta1_y_zz_xyzzz_1[i] * pc_x[i];

        ta1_y_xzz_xzzzz_0[i] = ta1_y_zz_zzzz_0[i] * fe_0 - ta1_y_zz_zzzz_1[i] * fe_0 + ta1_y_zz_xzzzz_0[i] * pa_x[i] - ta1_y_zz_xzzzz_1[i] * pc_x[i];

        ta1_y_xzz_yyyyy_0[i] = ta1_y_zz_yyyyy_0[i] * pa_x[i] - ta1_y_zz_yyyyy_1[i] * pc_x[i];

        ta1_y_xzz_yyyyz_0[i] = ta1_y_zz_yyyyz_0[i] * pa_x[i] - ta1_y_zz_yyyyz_1[i] * pc_x[i];

        ta1_y_xzz_yyyzz_0[i] = ta1_y_zz_yyyzz_0[i] * pa_x[i] - ta1_y_zz_yyyzz_1[i] * pc_x[i];

        ta1_y_xzz_yyzzz_0[i] = ta1_y_zz_yyzzz_0[i] * pa_x[i] - ta1_y_zz_yyzzz_1[i] * pc_x[i];

        ta1_y_xzz_yzzzz_0[i] = ta1_y_zz_yzzzz_0[i] * pa_x[i] - ta1_y_zz_yzzzz_1[i] * pc_x[i];

        ta1_y_xzz_zzzzz_0[i] = ta1_y_zz_zzzzz_0[i] * pa_x[i] - ta1_y_zz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 336-357 components of targeted buffer : FH

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

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_y_y_xxxxx_0,   \
                             ta1_y_y_xxxxx_1,   \
                             ta1_y_y_xxxxy_0,   \
                             ta1_y_y_xxxxy_1,   \
                             ta1_y_y_xxxxz_0,   \
                             ta1_y_y_xxxxz_1,   \
                             ta1_y_y_xxxyy_0,   \
                             ta1_y_y_xxxyy_1,   \
                             ta1_y_y_xxxyz_0,   \
                             ta1_y_y_xxxyz_1,   \
                             ta1_y_y_xxxzz_0,   \
                             ta1_y_y_xxxzz_1,   \
                             ta1_y_y_xxyyy_0,   \
                             ta1_y_y_xxyyy_1,   \
                             ta1_y_y_xxyyz_0,   \
                             ta1_y_y_xxyyz_1,   \
                             ta1_y_y_xxyzz_0,   \
                             ta1_y_y_xxyzz_1,   \
                             ta1_y_y_xxzzz_0,   \
                             ta1_y_y_xxzzz_1,   \
                             ta1_y_y_xyyyy_0,   \
                             ta1_y_y_xyyyy_1,   \
                             ta1_y_y_xyyyz_0,   \
                             ta1_y_y_xyyyz_1,   \
                             ta1_y_y_xyyzz_0,   \
                             ta1_y_y_xyyzz_1,   \
                             ta1_y_y_xyzzz_0,   \
                             ta1_y_y_xyzzz_1,   \
                             ta1_y_y_xzzzz_0,   \
                             ta1_y_y_xzzzz_1,   \
                             ta1_y_y_yyyyy_0,   \
                             ta1_y_y_yyyyy_1,   \
                             ta1_y_y_yyyyz_0,   \
                             ta1_y_y_yyyyz_1,   \
                             ta1_y_y_yyyzz_0,   \
                             ta1_y_y_yyyzz_1,   \
                             ta1_y_y_yyzzz_0,   \
                             ta1_y_y_yyzzz_1,   \
                             ta1_y_y_yzzzz_0,   \
                             ta1_y_y_yzzzz_1,   \
                             ta1_y_y_zzzzz_0,   \
                             ta1_y_y_zzzzz_1,   \
                             ta1_y_yy_xxxx_0,   \
                             ta1_y_yy_xxxx_1,   \
                             ta1_y_yy_xxxxx_0,  \
                             ta1_y_yy_xxxxx_1,  \
                             ta1_y_yy_xxxxy_0,  \
                             ta1_y_yy_xxxxy_1,  \
                             ta1_y_yy_xxxxz_0,  \
                             ta1_y_yy_xxxxz_1,  \
                             ta1_y_yy_xxxy_0,   \
                             ta1_y_yy_xxxy_1,   \
                             ta1_y_yy_xxxyy_0,  \
                             ta1_y_yy_xxxyy_1,  \
                             ta1_y_yy_xxxyz_0,  \
                             ta1_y_yy_xxxyz_1,  \
                             ta1_y_yy_xxxz_0,   \
                             ta1_y_yy_xxxz_1,   \
                             ta1_y_yy_xxxzz_0,  \
                             ta1_y_yy_xxxzz_1,  \
                             ta1_y_yy_xxyy_0,   \
                             ta1_y_yy_xxyy_1,   \
                             ta1_y_yy_xxyyy_0,  \
                             ta1_y_yy_xxyyy_1,  \
                             ta1_y_yy_xxyyz_0,  \
                             ta1_y_yy_xxyyz_1,  \
                             ta1_y_yy_xxyz_0,   \
                             ta1_y_yy_xxyz_1,   \
                             ta1_y_yy_xxyzz_0,  \
                             ta1_y_yy_xxyzz_1,  \
                             ta1_y_yy_xxzz_0,   \
                             ta1_y_yy_xxzz_1,   \
                             ta1_y_yy_xxzzz_0,  \
                             ta1_y_yy_xxzzz_1,  \
                             ta1_y_yy_xyyy_0,   \
                             ta1_y_yy_xyyy_1,   \
                             ta1_y_yy_xyyyy_0,  \
                             ta1_y_yy_xyyyy_1,  \
                             ta1_y_yy_xyyyz_0,  \
                             ta1_y_yy_xyyyz_1,  \
                             ta1_y_yy_xyyz_0,   \
                             ta1_y_yy_xyyz_1,   \
                             ta1_y_yy_xyyzz_0,  \
                             ta1_y_yy_xyyzz_1,  \
                             ta1_y_yy_xyzz_0,   \
                             ta1_y_yy_xyzz_1,   \
                             ta1_y_yy_xyzzz_0,  \
                             ta1_y_yy_xyzzz_1,  \
                             ta1_y_yy_xzzz_0,   \
                             ta1_y_yy_xzzz_1,   \
                             ta1_y_yy_xzzzz_0,  \
                             ta1_y_yy_xzzzz_1,  \
                             ta1_y_yy_yyyy_0,   \
                             ta1_y_yy_yyyy_1,   \
                             ta1_y_yy_yyyyy_0,  \
                             ta1_y_yy_yyyyy_1,  \
                             ta1_y_yy_yyyyz_0,  \
                             ta1_y_yy_yyyyz_1,  \
                             ta1_y_yy_yyyz_0,   \
                             ta1_y_yy_yyyz_1,   \
                             ta1_y_yy_yyyzz_0,  \
                             ta1_y_yy_yyyzz_1,  \
                             ta1_y_yy_yyzz_0,   \
                             ta1_y_yy_yyzz_1,   \
                             ta1_y_yy_yyzzz_0,  \
                             ta1_y_yy_yyzzz_1,  \
                             ta1_y_yy_yzzz_0,   \
                             ta1_y_yy_yzzz_1,   \
                             ta1_y_yy_yzzzz_0,  \
                             ta1_y_yy_yzzzz_1,  \
                             ta1_y_yy_zzzz_0,   \
                             ta1_y_yy_zzzz_1,   \
                             ta1_y_yy_zzzzz_0,  \
                             ta1_y_yy_zzzzz_1,  \
                             ta1_y_yyy_xxxxx_0, \
                             ta1_y_yyy_xxxxy_0, \
                             ta1_y_yyy_xxxxz_0, \
                             ta1_y_yyy_xxxyy_0, \
                             ta1_y_yyy_xxxyz_0, \
                             ta1_y_yyy_xxxzz_0, \
                             ta1_y_yyy_xxyyy_0, \
                             ta1_y_yyy_xxyyz_0, \
                             ta1_y_yyy_xxyzz_0, \
                             ta1_y_yyy_xxzzz_0, \
                             ta1_y_yyy_xyyyy_0, \
                             ta1_y_yyy_xyyyz_0, \
                             ta1_y_yyy_xyyzz_0, \
                             ta1_y_yyy_xyzzz_0, \
                             ta1_y_yyy_xzzzz_0, \
                             ta1_y_yyy_yyyyy_0, \
                             ta1_y_yyy_yyyyz_0, \
                             ta1_y_yyy_yyyzz_0, \
                             ta1_y_yyy_yyzzz_0, \
                             ta1_y_yyy_yzzzz_0, \
                             ta1_y_yyy_zzzzz_0, \
                             ta_yy_xxxxx_1,     \
                             ta_yy_xxxxy_1,     \
                             ta_yy_xxxxz_1,     \
                             ta_yy_xxxyy_1,     \
                             ta_yy_xxxyz_1,     \
                             ta_yy_xxxzz_1,     \
                             ta_yy_xxyyy_1,     \
                             ta_yy_xxyyz_1,     \
                             ta_yy_xxyzz_1,     \
                             ta_yy_xxzzz_1,     \
                             ta_yy_xyyyy_1,     \
                             ta_yy_xyyyz_1,     \
                             ta_yy_xyyzz_1,     \
                             ta_yy_xyzzz_1,     \
                             ta_yy_xzzzz_1,     \
                             ta_yy_yyyyy_1,     \
                             ta_yy_yyyyz_1,     \
                             ta_yy_yyyzz_1,     \
                             ta_yy_yyzzz_1,     \
                             ta_yy_yzzzz_1,     \
                             ta_yy_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyy_xxxxx_0[i] = 2.0 * ta1_y_y_xxxxx_0[i] * fe_0 - 2.0 * ta1_y_y_xxxxx_1[i] * fe_0 + ta_yy_xxxxx_1[i] + ta1_y_yy_xxxxx_0[i] * pa_y[i] -
                               ta1_y_yy_xxxxx_1[i] * pc_y[i];

        ta1_y_yyy_xxxxy_0[i] = 2.0 * ta1_y_y_xxxxy_0[i] * fe_0 - 2.0 * ta1_y_y_xxxxy_1[i] * fe_0 + ta1_y_yy_xxxx_0[i] * fe_0 -
                               ta1_y_yy_xxxx_1[i] * fe_0 + ta_yy_xxxxy_1[i] + ta1_y_yy_xxxxy_0[i] * pa_y[i] - ta1_y_yy_xxxxy_1[i] * pc_y[i];

        ta1_y_yyy_xxxxz_0[i] = 2.0 * ta1_y_y_xxxxz_0[i] * fe_0 - 2.0 * ta1_y_y_xxxxz_1[i] * fe_0 + ta_yy_xxxxz_1[i] + ta1_y_yy_xxxxz_0[i] * pa_y[i] -
                               ta1_y_yy_xxxxz_1[i] * pc_y[i];

        ta1_y_yyy_xxxyy_0[i] = 2.0 * ta1_y_y_xxxyy_0[i] * fe_0 - 2.0 * ta1_y_y_xxxyy_1[i] * fe_0 + 2.0 * ta1_y_yy_xxxy_0[i] * fe_0 -
                               2.0 * ta1_y_yy_xxxy_1[i] * fe_0 + ta_yy_xxxyy_1[i] + ta1_y_yy_xxxyy_0[i] * pa_y[i] - ta1_y_yy_xxxyy_1[i] * pc_y[i];

        ta1_y_yyy_xxxyz_0[i] = 2.0 * ta1_y_y_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_y_xxxyz_1[i] * fe_0 + ta1_y_yy_xxxz_0[i] * fe_0 -
                               ta1_y_yy_xxxz_1[i] * fe_0 + ta_yy_xxxyz_1[i] + ta1_y_yy_xxxyz_0[i] * pa_y[i] - ta1_y_yy_xxxyz_1[i] * pc_y[i];

        ta1_y_yyy_xxxzz_0[i] = 2.0 * ta1_y_y_xxxzz_0[i] * fe_0 - 2.0 * ta1_y_y_xxxzz_1[i] * fe_0 + ta_yy_xxxzz_1[i] + ta1_y_yy_xxxzz_0[i] * pa_y[i] -
                               ta1_y_yy_xxxzz_1[i] * pc_y[i];

        ta1_y_yyy_xxyyy_0[i] = 2.0 * ta1_y_y_xxyyy_0[i] * fe_0 - 2.0 * ta1_y_y_xxyyy_1[i] * fe_0 + 3.0 * ta1_y_yy_xxyy_0[i] * fe_0 -
                               3.0 * ta1_y_yy_xxyy_1[i] * fe_0 + ta_yy_xxyyy_1[i] + ta1_y_yy_xxyyy_0[i] * pa_y[i] - ta1_y_yy_xxyyy_1[i] * pc_y[i];

        ta1_y_yyy_xxyyz_0[i] = 2.0 * ta1_y_y_xxyyz_0[i] * fe_0 - 2.0 * ta1_y_y_xxyyz_1[i] * fe_0 + 2.0 * ta1_y_yy_xxyz_0[i] * fe_0 -
                               2.0 * ta1_y_yy_xxyz_1[i] * fe_0 + ta_yy_xxyyz_1[i] + ta1_y_yy_xxyyz_0[i] * pa_y[i] - ta1_y_yy_xxyyz_1[i] * pc_y[i];

        ta1_y_yyy_xxyzz_0[i] = 2.0 * ta1_y_y_xxyzz_0[i] * fe_0 - 2.0 * ta1_y_y_xxyzz_1[i] * fe_0 + ta1_y_yy_xxzz_0[i] * fe_0 -
                               ta1_y_yy_xxzz_1[i] * fe_0 + ta_yy_xxyzz_1[i] + ta1_y_yy_xxyzz_0[i] * pa_y[i] - ta1_y_yy_xxyzz_1[i] * pc_y[i];

        ta1_y_yyy_xxzzz_0[i] = 2.0 * ta1_y_y_xxzzz_0[i] * fe_0 - 2.0 * ta1_y_y_xxzzz_1[i] * fe_0 + ta_yy_xxzzz_1[i] + ta1_y_yy_xxzzz_0[i] * pa_y[i] -
                               ta1_y_yy_xxzzz_1[i] * pc_y[i];

        ta1_y_yyy_xyyyy_0[i] = 2.0 * ta1_y_y_xyyyy_0[i] * fe_0 - 2.0 * ta1_y_y_xyyyy_1[i] * fe_0 + 4.0 * ta1_y_yy_xyyy_0[i] * fe_0 -
                               4.0 * ta1_y_yy_xyyy_1[i] * fe_0 + ta_yy_xyyyy_1[i] + ta1_y_yy_xyyyy_0[i] * pa_y[i] - ta1_y_yy_xyyyy_1[i] * pc_y[i];

        ta1_y_yyy_xyyyz_0[i] = 2.0 * ta1_y_y_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_y_xyyyz_1[i] * fe_0 + 3.0 * ta1_y_yy_xyyz_0[i] * fe_0 -
                               3.0 * ta1_y_yy_xyyz_1[i] * fe_0 + ta_yy_xyyyz_1[i] + ta1_y_yy_xyyyz_0[i] * pa_y[i] - ta1_y_yy_xyyyz_1[i] * pc_y[i];

        ta1_y_yyy_xyyzz_0[i] = 2.0 * ta1_y_y_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_y_xyyzz_1[i] * fe_0 + 2.0 * ta1_y_yy_xyzz_0[i] * fe_0 -
                               2.0 * ta1_y_yy_xyzz_1[i] * fe_0 + ta_yy_xyyzz_1[i] + ta1_y_yy_xyyzz_0[i] * pa_y[i] - ta1_y_yy_xyyzz_1[i] * pc_y[i];

        ta1_y_yyy_xyzzz_0[i] = 2.0 * ta1_y_y_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_y_xyzzz_1[i] * fe_0 + ta1_y_yy_xzzz_0[i] * fe_0 -
                               ta1_y_yy_xzzz_1[i] * fe_0 + ta_yy_xyzzz_1[i] + ta1_y_yy_xyzzz_0[i] * pa_y[i] - ta1_y_yy_xyzzz_1[i] * pc_y[i];

        ta1_y_yyy_xzzzz_0[i] = 2.0 * ta1_y_y_xzzzz_0[i] * fe_0 - 2.0 * ta1_y_y_xzzzz_1[i] * fe_0 + ta_yy_xzzzz_1[i] + ta1_y_yy_xzzzz_0[i] * pa_y[i] -
                               ta1_y_yy_xzzzz_1[i] * pc_y[i];

        ta1_y_yyy_yyyyy_0[i] = 2.0 * ta1_y_y_yyyyy_0[i] * fe_0 - 2.0 * ta1_y_y_yyyyy_1[i] * fe_0 + 5.0 * ta1_y_yy_yyyy_0[i] * fe_0 -
                               5.0 * ta1_y_yy_yyyy_1[i] * fe_0 + ta_yy_yyyyy_1[i] + ta1_y_yy_yyyyy_0[i] * pa_y[i] - ta1_y_yy_yyyyy_1[i] * pc_y[i];

        ta1_y_yyy_yyyyz_0[i] = 2.0 * ta1_y_y_yyyyz_0[i] * fe_0 - 2.0 * ta1_y_y_yyyyz_1[i] * fe_0 + 4.0 * ta1_y_yy_yyyz_0[i] * fe_0 -
                               4.0 * ta1_y_yy_yyyz_1[i] * fe_0 + ta_yy_yyyyz_1[i] + ta1_y_yy_yyyyz_0[i] * pa_y[i] - ta1_y_yy_yyyyz_1[i] * pc_y[i];

        ta1_y_yyy_yyyzz_0[i] = 2.0 * ta1_y_y_yyyzz_0[i] * fe_0 - 2.0 * ta1_y_y_yyyzz_1[i] * fe_0 + 3.0 * ta1_y_yy_yyzz_0[i] * fe_0 -
                               3.0 * ta1_y_yy_yyzz_1[i] * fe_0 + ta_yy_yyyzz_1[i] + ta1_y_yy_yyyzz_0[i] * pa_y[i] - ta1_y_yy_yyyzz_1[i] * pc_y[i];

        ta1_y_yyy_yyzzz_0[i] = 2.0 * ta1_y_y_yyzzz_0[i] * fe_0 - 2.0 * ta1_y_y_yyzzz_1[i] * fe_0 + 2.0 * ta1_y_yy_yzzz_0[i] * fe_0 -
                               2.0 * ta1_y_yy_yzzz_1[i] * fe_0 + ta_yy_yyzzz_1[i] + ta1_y_yy_yyzzz_0[i] * pa_y[i] - ta1_y_yy_yyzzz_1[i] * pc_y[i];

        ta1_y_yyy_yzzzz_0[i] = 2.0 * ta1_y_y_yzzzz_0[i] * fe_0 - 2.0 * ta1_y_y_yzzzz_1[i] * fe_0 + ta1_y_yy_zzzz_0[i] * fe_0 -
                               ta1_y_yy_zzzz_1[i] * fe_0 + ta_yy_yzzzz_1[i] + ta1_y_yy_yzzzz_0[i] * pa_y[i] - ta1_y_yy_yzzzz_1[i] * pc_y[i];

        ta1_y_yyy_zzzzz_0[i] = 2.0 * ta1_y_y_zzzzz_0[i] * fe_0 - 2.0 * ta1_y_y_zzzzz_1[i] * fe_0 + ta_yy_zzzzz_1[i] + ta1_y_yy_zzzzz_0[i] * pa_y[i] -
                               ta1_y_yy_zzzzz_1[i] * pc_y[i];
    }

    // Set up 357-378 components of targeted buffer : FH

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

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_y_yy_xxxx_0,   \
                             ta1_y_yy_xxxx_1,   \
                             ta1_y_yy_xxxxx_0,  \
                             ta1_y_yy_xxxxx_1,  \
                             ta1_y_yy_xxxxy_0,  \
                             ta1_y_yy_xxxxy_1,  \
                             ta1_y_yy_xxxxz_0,  \
                             ta1_y_yy_xxxxz_1,  \
                             ta1_y_yy_xxxy_0,   \
                             ta1_y_yy_xxxy_1,   \
                             ta1_y_yy_xxxyy_0,  \
                             ta1_y_yy_xxxyy_1,  \
                             ta1_y_yy_xxxyz_0,  \
                             ta1_y_yy_xxxyz_1,  \
                             ta1_y_yy_xxxz_0,   \
                             ta1_y_yy_xxxz_1,   \
                             ta1_y_yy_xxxzz_0,  \
                             ta1_y_yy_xxxzz_1,  \
                             ta1_y_yy_xxyy_0,   \
                             ta1_y_yy_xxyy_1,   \
                             ta1_y_yy_xxyyy_0,  \
                             ta1_y_yy_xxyyy_1,  \
                             ta1_y_yy_xxyyz_0,  \
                             ta1_y_yy_xxyyz_1,  \
                             ta1_y_yy_xxyz_0,   \
                             ta1_y_yy_xxyz_1,   \
                             ta1_y_yy_xxyzz_0,  \
                             ta1_y_yy_xxyzz_1,  \
                             ta1_y_yy_xxzz_0,   \
                             ta1_y_yy_xxzz_1,   \
                             ta1_y_yy_xxzzz_0,  \
                             ta1_y_yy_xxzzz_1,  \
                             ta1_y_yy_xyyy_0,   \
                             ta1_y_yy_xyyy_1,   \
                             ta1_y_yy_xyyyy_0,  \
                             ta1_y_yy_xyyyy_1,  \
                             ta1_y_yy_xyyyz_0,  \
                             ta1_y_yy_xyyyz_1,  \
                             ta1_y_yy_xyyz_0,   \
                             ta1_y_yy_xyyz_1,   \
                             ta1_y_yy_xyyzz_0,  \
                             ta1_y_yy_xyyzz_1,  \
                             ta1_y_yy_xyzz_0,   \
                             ta1_y_yy_xyzz_1,   \
                             ta1_y_yy_xyzzz_0,  \
                             ta1_y_yy_xyzzz_1,  \
                             ta1_y_yy_xzzz_0,   \
                             ta1_y_yy_xzzz_1,   \
                             ta1_y_yy_xzzzz_0,  \
                             ta1_y_yy_xzzzz_1,  \
                             ta1_y_yy_yyyy_0,   \
                             ta1_y_yy_yyyy_1,   \
                             ta1_y_yy_yyyyy_0,  \
                             ta1_y_yy_yyyyy_1,  \
                             ta1_y_yy_yyyyz_0,  \
                             ta1_y_yy_yyyyz_1,  \
                             ta1_y_yy_yyyz_0,   \
                             ta1_y_yy_yyyz_1,   \
                             ta1_y_yy_yyyzz_0,  \
                             ta1_y_yy_yyyzz_1,  \
                             ta1_y_yy_yyzz_0,   \
                             ta1_y_yy_yyzz_1,   \
                             ta1_y_yy_yyzzz_0,  \
                             ta1_y_yy_yyzzz_1,  \
                             ta1_y_yy_yzzz_0,   \
                             ta1_y_yy_yzzz_1,   \
                             ta1_y_yy_yzzzz_0,  \
                             ta1_y_yy_yzzzz_1,  \
                             ta1_y_yy_zzzz_0,   \
                             ta1_y_yy_zzzz_1,   \
                             ta1_y_yy_zzzzz_0,  \
                             ta1_y_yy_zzzzz_1,  \
                             ta1_y_yyz_xxxxx_0, \
                             ta1_y_yyz_xxxxy_0, \
                             ta1_y_yyz_xxxxz_0, \
                             ta1_y_yyz_xxxyy_0, \
                             ta1_y_yyz_xxxyz_0, \
                             ta1_y_yyz_xxxzz_0, \
                             ta1_y_yyz_xxyyy_0, \
                             ta1_y_yyz_xxyyz_0, \
                             ta1_y_yyz_xxyzz_0, \
                             ta1_y_yyz_xxzzz_0, \
                             ta1_y_yyz_xyyyy_0, \
                             ta1_y_yyz_xyyyz_0, \
                             ta1_y_yyz_xyyzz_0, \
                             ta1_y_yyz_xyzzz_0, \
                             ta1_y_yyz_xzzzz_0, \
                             ta1_y_yyz_yyyyy_0, \
                             ta1_y_yyz_yyyyz_0, \
                             ta1_y_yyz_yyyzz_0, \
                             ta1_y_yyz_yyzzz_0, \
                             ta1_y_yyz_yzzzz_0, \
                             ta1_y_yyz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyz_xxxxx_0[i] = ta1_y_yy_xxxxx_0[i] * pa_z[i] - ta1_y_yy_xxxxx_1[i] * pc_z[i];

        ta1_y_yyz_xxxxy_0[i] = ta1_y_yy_xxxxy_0[i] * pa_z[i] - ta1_y_yy_xxxxy_1[i] * pc_z[i];

        ta1_y_yyz_xxxxz_0[i] = ta1_y_yy_xxxx_0[i] * fe_0 - ta1_y_yy_xxxx_1[i] * fe_0 + ta1_y_yy_xxxxz_0[i] * pa_z[i] - ta1_y_yy_xxxxz_1[i] * pc_z[i];

        ta1_y_yyz_xxxyy_0[i] = ta1_y_yy_xxxyy_0[i] * pa_z[i] - ta1_y_yy_xxxyy_1[i] * pc_z[i];

        ta1_y_yyz_xxxyz_0[i] = ta1_y_yy_xxxy_0[i] * fe_0 - ta1_y_yy_xxxy_1[i] * fe_0 + ta1_y_yy_xxxyz_0[i] * pa_z[i] - ta1_y_yy_xxxyz_1[i] * pc_z[i];

        ta1_y_yyz_xxxzz_0[i] =
            2.0 * ta1_y_yy_xxxz_0[i] * fe_0 - 2.0 * ta1_y_yy_xxxz_1[i] * fe_0 + ta1_y_yy_xxxzz_0[i] * pa_z[i] - ta1_y_yy_xxxzz_1[i] * pc_z[i];

        ta1_y_yyz_xxyyy_0[i] = ta1_y_yy_xxyyy_0[i] * pa_z[i] - ta1_y_yy_xxyyy_1[i] * pc_z[i];

        ta1_y_yyz_xxyyz_0[i] = ta1_y_yy_xxyy_0[i] * fe_0 - ta1_y_yy_xxyy_1[i] * fe_0 + ta1_y_yy_xxyyz_0[i] * pa_z[i] - ta1_y_yy_xxyyz_1[i] * pc_z[i];

        ta1_y_yyz_xxyzz_0[i] =
            2.0 * ta1_y_yy_xxyz_0[i] * fe_0 - 2.0 * ta1_y_yy_xxyz_1[i] * fe_0 + ta1_y_yy_xxyzz_0[i] * pa_z[i] - ta1_y_yy_xxyzz_1[i] * pc_z[i];

        ta1_y_yyz_xxzzz_0[i] =
            3.0 * ta1_y_yy_xxzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxzz_1[i] * fe_0 + ta1_y_yy_xxzzz_0[i] * pa_z[i] - ta1_y_yy_xxzzz_1[i] * pc_z[i];

        ta1_y_yyz_xyyyy_0[i] = ta1_y_yy_xyyyy_0[i] * pa_z[i] - ta1_y_yy_xyyyy_1[i] * pc_z[i];

        ta1_y_yyz_xyyyz_0[i] = ta1_y_yy_xyyy_0[i] * fe_0 - ta1_y_yy_xyyy_1[i] * fe_0 + ta1_y_yy_xyyyz_0[i] * pa_z[i] - ta1_y_yy_xyyyz_1[i] * pc_z[i];

        ta1_y_yyz_xyyzz_0[i] =
            2.0 * ta1_y_yy_xyyz_0[i] * fe_0 - 2.0 * ta1_y_yy_xyyz_1[i] * fe_0 + ta1_y_yy_xyyzz_0[i] * pa_z[i] - ta1_y_yy_xyyzz_1[i] * pc_z[i];

        ta1_y_yyz_xyzzz_0[i] =
            3.0 * ta1_y_yy_xyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xyzz_1[i] * fe_0 + ta1_y_yy_xyzzz_0[i] * pa_z[i] - ta1_y_yy_xyzzz_1[i] * pc_z[i];

        ta1_y_yyz_xzzzz_0[i] =
            4.0 * ta1_y_yy_xzzz_0[i] * fe_0 - 4.0 * ta1_y_yy_xzzz_1[i] * fe_0 + ta1_y_yy_xzzzz_0[i] * pa_z[i] - ta1_y_yy_xzzzz_1[i] * pc_z[i];

        ta1_y_yyz_yyyyy_0[i] = ta1_y_yy_yyyyy_0[i] * pa_z[i] - ta1_y_yy_yyyyy_1[i] * pc_z[i];

        ta1_y_yyz_yyyyz_0[i] = ta1_y_yy_yyyy_0[i] * fe_0 - ta1_y_yy_yyyy_1[i] * fe_0 + ta1_y_yy_yyyyz_0[i] * pa_z[i] - ta1_y_yy_yyyyz_1[i] * pc_z[i];

        ta1_y_yyz_yyyzz_0[i] =
            2.0 * ta1_y_yy_yyyz_0[i] * fe_0 - 2.0 * ta1_y_yy_yyyz_1[i] * fe_0 + ta1_y_yy_yyyzz_0[i] * pa_z[i] - ta1_y_yy_yyyzz_1[i] * pc_z[i];

        ta1_y_yyz_yyzzz_0[i] =
            3.0 * ta1_y_yy_yyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyzz_1[i] * fe_0 + ta1_y_yy_yyzzz_0[i] * pa_z[i] - ta1_y_yy_yyzzz_1[i] * pc_z[i];

        ta1_y_yyz_yzzzz_0[i] =
            4.0 * ta1_y_yy_yzzz_0[i] * fe_0 - 4.0 * ta1_y_yy_yzzz_1[i] * fe_0 + ta1_y_yy_yzzzz_0[i] * pa_z[i] - ta1_y_yy_yzzzz_1[i] * pc_z[i];

        ta1_y_yyz_zzzzz_0[i] =
            5.0 * ta1_y_yy_zzzz_0[i] * fe_0 - 5.0 * ta1_y_yy_zzzz_1[i] * fe_0 + ta1_y_yy_zzzzz_0[i] * pa_z[i] - ta1_y_yy_zzzzz_1[i] * pc_z[i];
    }

    // Set up 378-399 components of targeted buffer : FH

    auto ta1_y_yzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 378);

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

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_y_xxxxy_0,   \
                             ta1_y_y_xxxxy_1,   \
                             ta1_y_y_xxxyy_0,   \
                             ta1_y_y_xxxyy_1,   \
                             ta1_y_y_xxyyy_0,   \
                             ta1_y_y_xxyyy_1,   \
                             ta1_y_y_xyyyy_0,   \
                             ta1_y_y_xyyyy_1,   \
                             ta1_y_y_yyyyy_0,   \
                             ta1_y_y_yyyyy_1,   \
                             ta1_y_yz_xxxxy_0,  \
                             ta1_y_yz_xxxxy_1,  \
                             ta1_y_yz_xxxyy_0,  \
                             ta1_y_yz_xxxyy_1,  \
                             ta1_y_yz_xxyyy_0,  \
                             ta1_y_yz_xxyyy_1,  \
                             ta1_y_yz_xyyyy_0,  \
                             ta1_y_yz_xyyyy_1,  \
                             ta1_y_yz_yyyyy_0,  \
                             ta1_y_yz_yyyyy_1,  \
                             ta1_y_yzz_xxxxx_0, \
                             ta1_y_yzz_xxxxy_0, \
                             ta1_y_yzz_xxxxz_0, \
                             ta1_y_yzz_xxxyy_0, \
                             ta1_y_yzz_xxxyz_0, \
                             ta1_y_yzz_xxxzz_0, \
                             ta1_y_yzz_xxyyy_0, \
                             ta1_y_yzz_xxyyz_0, \
                             ta1_y_yzz_xxyzz_0, \
                             ta1_y_yzz_xxzzz_0, \
                             ta1_y_yzz_xyyyy_0, \
                             ta1_y_yzz_xyyyz_0, \
                             ta1_y_yzz_xyyzz_0, \
                             ta1_y_yzz_xyzzz_0, \
                             ta1_y_yzz_xzzzz_0, \
                             ta1_y_yzz_yyyyy_0, \
                             ta1_y_yzz_yyyyz_0, \
                             ta1_y_yzz_yyyzz_0, \
                             ta1_y_yzz_yyzzz_0, \
                             ta1_y_yzz_yzzzz_0, \
                             ta1_y_yzz_zzzzz_0, \
                             ta1_y_zz_xxxxx_0,  \
                             ta1_y_zz_xxxxx_1,  \
                             ta1_y_zz_xxxxz_0,  \
                             ta1_y_zz_xxxxz_1,  \
                             ta1_y_zz_xxxyz_0,  \
                             ta1_y_zz_xxxyz_1,  \
                             ta1_y_zz_xxxz_0,   \
                             ta1_y_zz_xxxz_1,   \
                             ta1_y_zz_xxxzz_0,  \
                             ta1_y_zz_xxxzz_1,  \
                             ta1_y_zz_xxyyz_0,  \
                             ta1_y_zz_xxyyz_1,  \
                             ta1_y_zz_xxyz_0,   \
                             ta1_y_zz_xxyz_1,   \
                             ta1_y_zz_xxyzz_0,  \
                             ta1_y_zz_xxyzz_1,  \
                             ta1_y_zz_xxzz_0,   \
                             ta1_y_zz_xxzz_1,   \
                             ta1_y_zz_xxzzz_0,  \
                             ta1_y_zz_xxzzz_1,  \
                             ta1_y_zz_xyyyz_0,  \
                             ta1_y_zz_xyyyz_1,  \
                             ta1_y_zz_xyyz_0,   \
                             ta1_y_zz_xyyz_1,   \
                             ta1_y_zz_xyyzz_0,  \
                             ta1_y_zz_xyyzz_1,  \
                             ta1_y_zz_xyzz_0,   \
                             ta1_y_zz_xyzz_1,   \
                             ta1_y_zz_xyzzz_0,  \
                             ta1_y_zz_xyzzz_1,  \
                             ta1_y_zz_xzzz_0,   \
                             ta1_y_zz_xzzz_1,   \
                             ta1_y_zz_xzzzz_0,  \
                             ta1_y_zz_xzzzz_1,  \
                             ta1_y_zz_yyyyz_0,  \
                             ta1_y_zz_yyyyz_1,  \
                             ta1_y_zz_yyyz_0,   \
                             ta1_y_zz_yyyz_1,   \
                             ta1_y_zz_yyyzz_0,  \
                             ta1_y_zz_yyyzz_1,  \
                             ta1_y_zz_yyzz_0,   \
                             ta1_y_zz_yyzz_1,   \
                             ta1_y_zz_yyzzz_0,  \
                             ta1_y_zz_yyzzz_1,  \
                             ta1_y_zz_yzzz_0,   \
                             ta1_y_zz_yzzz_1,   \
                             ta1_y_zz_yzzzz_0,  \
                             ta1_y_zz_yzzzz_1,  \
                             ta1_y_zz_zzzz_0,   \
                             ta1_y_zz_zzzz_1,   \
                             ta1_y_zz_zzzzz_0,  \
                             ta1_y_zz_zzzzz_1,  \
                             ta_zz_xxxxx_1,     \
                             ta_zz_xxxxz_1,     \
                             ta_zz_xxxyz_1,     \
                             ta_zz_xxxzz_1,     \
                             ta_zz_xxyyz_1,     \
                             ta_zz_xxyzz_1,     \
                             ta_zz_xxzzz_1,     \
                             ta_zz_xyyyz_1,     \
                             ta_zz_xyyzz_1,     \
                             ta_zz_xyzzz_1,     \
                             ta_zz_xzzzz_1,     \
                             ta_zz_yyyyz_1,     \
                             ta_zz_yyyzz_1,     \
                             ta_zz_yyzzz_1,     \
                             ta_zz_yzzzz_1,     \
                             ta_zz_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzz_xxxxx_0[i] = ta_zz_xxxxx_1[i] + ta1_y_zz_xxxxx_0[i] * pa_y[i] - ta1_y_zz_xxxxx_1[i] * pc_y[i];

        ta1_y_yzz_xxxxy_0[i] = ta1_y_y_xxxxy_0[i] * fe_0 - ta1_y_y_xxxxy_1[i] * fe_0 + ta1_y_yz_xxxxy_0[i] * pa_z[i] - ta1_y_yz_xxxxy_1[i] * pc_z[i];

        ta1_y_yzz_xxxxz_0[i] = ta_zz_xxxxz_1[i] + ta1_y_zz_xxxxz_0[i] * pa_y[i] - ta1_y_zz_xxxxz_1[i] * pc_y[i];

        ta1_y_yzz_xxxyy_0[i] = ta1_y_y_xxxyy_0[i] * fe_0 - ta1_y_y_xxxyy_1[i] * fe_0 + ta1_y_yz_xxxyy_0[i] * pa_z[i] - ta1_y_yz_xxxyy_1[i] * pc_z[i];

        ta1_y_yzz_xxxyz_0[i] =
            ta1_y_zz_xxxz_0[i] * fe_0 - ta1_y_zz_xxxz_1[i] * fe_0 + ta_zz_xxxyz_1[i] + ta1_y_zz_xxxyz_0[i] * pa_y[i] - ta1_y_zz_xxxyz_1[i] * pc_y[i];

        ta1_y_yzz_xxxzz_0[i] = ta_zz_xxxzz_1[i] + ta1_y_zz_xxxzz_0[i] * pa_y[i] - ta1_y_zz_xxxzz_1[i] * pc_y[i];

        ta1_y_yzz_xxyyy_0[i] = ta1_y_y_xxyyy_0[i] * fe_0 - ta1_y_y_xxyyy_1[i] * fe_0 + ta1_y_yz_xxyyy_0[i] * pa_z[i] - ta1_y_yz_xxyyy_1[i] * pc_z[i];

        ta1_y_yzz_xxyyz_0[i] = 2.0 * ta1_y_zz_xxyz_0[i] * fe_0 - 2.0 * ta1_y_zz_xxyz_1[i] * fe_0 + ta_zz_xxyyz_1[i] + ta1_y_zz_xxyyz_0[i] * pa_y[i] -
                               ta1_y_zz_xxyyz_1[i] * pc_y[i];

        ta1_y_yzz_xxyzz_0[i] =
            ta1_y_zz_xxzz_0[i] * fe_0 - ta1_y_zz_xxzz_1[i] * fe_0 + ta_zz_xxyzz_1[i] + ta1_y_zz_xxyzz_0[i] * pa_y[i] - ta1_y_zz_xxyzz_1[i] * pc_y[i];

        ta1_y_yzz_xxzzz_0[i] = ta_zz_xxzzz_1[i] + ta1_y_zz_xxzzz_0[i] * pa_y[i] - ta1_y_zz_xxzzz_1[i] * pc_y[i];

        ta1_y_yzz_xyyyy_0[i] = ta1_y_y_xyyyy_0[i] * fe_0 - ta1_y_y_xyyyy_1[i] * fe_0 + ta1_y_yz_xyyyy_0[i] * pa_z[i] - ta1_y_yz_xyyyy_1[i] * pc_z[i];

        ta1_y_yzz_xyyyz_0[i] = 3.0 * ta1_y_zz_xyyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xyyz_1[i] * fe_0 + ta_zz_xyyyz_1[i] + ta1_y_zz_xyyyz_0[i] * pa_y[i] -
                               ta1_y_zz_xyyyz_1[i] * pc_y[i];

        ta1_y_yzz_xyyzz_0[i] = 2.0 * ta1_y_zz_xyzz_0[i] * fe_0 - 2.0 * ta1_y_zz_xyzz_1[i] * fe_0 + ta_zz_xyyzz_1[i] + ta1_y_zz_xyyzz_0[i] * pa_y[i] -
                               ta1_y_zz_xyyzz_1[i] * pc_y[i];

        ta1_y_yzz_xyzzz_0[i] =
            ta1_y_zz_xzzz_0[i] * fe_0 - ta1_y_zz_xzzz_1[i] * fe_0 + ta_zz_xyzzz_1[i] + ta1_y_zz_xyzzz_0[i] * pa_y[i] - ta1_y_zz_xyzzz_1[i] * pc_y[i];

        ta1_y_yzz_xzzzz_0[i] = ta_zz_xzzzz_1[i] + ta1_y_zz_xzzzz_0[i] * pa_y[i] - ta1_y_zz_xzzzz_1[i] * pc_y[i];

        ta1_y_yzz_yyyyy_0[i] = ta1_y_y_yyyyy_0[i] * fe_0 - ta1_y_y_yyyyy_1[i] * fe_0 + ta1_y_yz_yyyyy_0[i] * pa_z[i] - ta1_y_yz_yyyyy_1[i] * pc_z[i];

        ta1_y_yzz_yyyyz_0[i] = 4.0 * ta1_y_zz_yyyz_0[i] * fe_0 - 4.0 * ta1_y_zz_yyyz_1[i] * fe_0 + ta_zz_yyyyz_1[i] + ta1_y_zz_yyyyz_0[i] * pa_y[i] -
                               ta1_y_zz_yyyyz_1[i] * pc_y[i];

        ta1_y_yzz_yyyzz_0[i] = 3.0 * ta1_y_zz_yyzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyzz_1[i] * fe_0 + ta_zz_yyyzz_1[i] + ta1_y_zz_yyyzz_0[i] * pa_y[i] -
                               ta1_y_zz_yyyzz_1[i] * pc_y[i];

        ta1_y_yzz_yyzzz_0[i] = 2.0 * ta1_y_zz_yzzz_0[i] * fe_0 - 2.0 * ta1_y_zz_yzzz_1[i] * fe_0 + ta_zz_yyzzz_1[i] + ta1_y_zz_yyzzz_0[i] * pa_y[i] -
                               ta1_y_zz_yyzzz_1[i] * pc_y[i];

        ta1_y_yzz_yzzzz_0[i] =
            ta1_y_zz_zzzz_0[i] * fe_0 - ta1_y_zz_zzzz_1[i] * fe_0 + ta_zz_yzzzz_1[i] + ta1_y_zz_yzzzz_0[i] * pa_y[i] - ta1_y_zz_yzzzz_1[i] * pc_y[i];

        ta1_y_yzz_zzzzz_0[i] = ta_zz_zzzzz_1[i] + ta1_y_zz_zzzzz_0[i] * pa_y[i] - ta1_y_zz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 399-420 components of targeted buffer : FH

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

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_y_z_xxxxx_0,   \
                             ta1_y_z_xxxxx_1,   \
                             ta1_y_z_xxxxy_0,   \
                             ta1_y_z_xxxxy_1,   \
                             ta1_y_z_xxxxz_0,   \
                             ta1_y_z_xxxxz_1,   \
                             ta1_y_z_xxxyy_0,   \
                             ta1_y_z_xxxyy_1,   \
                             ta1_y_z_xxxyz_0,   \
                             ta1_y_z_xxxyz_1,   \
                             ta1_y_z_xxxzz_0,   \
                             ta1_y_z_xxxzz_1,   \
                             ta1_y_z_xxyyy_0,   \
                             ta1_y_z_xxyyy_1,   \
                             ta1_y_z_xxyyz_0,   \
                             ta1_y_z_xxyyz_1,   \
                             ta1_y_z_xxyzz_0,   \
                             ta1_y_z_xxyzz_1,   \
                             ta1_y_z_xxzzz_0,   \
                             ta1_y_z_xxzzz_1,   \
                             ta1_y_z_xyyyy_0,   \
                             ta1_y_z_xyyyy_1,   \
                             ta1_y_z_xyyyz_0,   \
                             ta1_y_z_xyyyz_1,   \
                             ta1_y_z_xyyzz_0,   \
                             ta1_y_z_xyyzz_1,   \
                             ta1_y_z_xyzzz_0,   \
                             ta1_y_z_xyzzz_1,   \
                             ta1_y_z_xzzzz_0,   \
                             ta1_y_z_xzzzz_1,   \
                             ta1_y_z_yyyyy_0,   \
                             ta1_y_z_yyyyy_1,   \
                             ta1_y_z_yyyyz_0,   \
                             ta1_y_z_yyyyz_1,   \
                             ta1_y_z_yyyzz_0,   \
                             ta1_y_z_yyyzz_1,   \
                             ta1_y_z_yyzzz_0,   \
                             ta1_y_z_yyzzz_1,   \
                             ta1_y_z_yzzzz_0,   \
                             ta1_y_z_yzzzz_1,   \
                             ta1_y_z_zzzzz_0,   \
                             ta1_y_z_zzzzz_1,   \
                             ta1_y_zz_xxxx_0,   \
                             ta1_y_zz_xxxx_1,   \
                             ta1_y_zz_xxxxx_0,  \
                             ta1_y_zz_xxxxx_1,  \
                             ta1_y_zz_xxxxy_0,  \
                             ta1_y_zz_xxxxy_1,  \
                             ta1_y_zz_xxxxz_0,  \
                             ta1_y_zz_xxxxz_1,  \
                             ta1_y_zz_xxxy_0,   \
                             ta1_y_zz_xxxy_1,   \
                             ta1_y_zz_xxxyy_0,  \
                             ta1_y_zz_xxxyy_1,  \
                             ta1_y_zz_xxxyz_0,  \
                             ta1_y_zz_xxxyz_1,  \
                             ta1_y_zz_xxxz_0,   \
                             ta1_y_zz_xxxz_1,   \
                             ta1_y_zz_xxxzz_0,  \
                             ta1_y_zz_xxxzz_1,  \
                             ta1_y_zz_xxyy_0,   \
                             ta1_y_zz_xxyy_1,   \
                             ta1_y_zz_xxyyy_0,  \
                             ta1_y_zz_xxyyy_1,  \
                             ta1_y_zz_xxyyz_0,  \
                             ta1_y_zz_xxyyz_1,  \
                             ta1_y_zz_xxyz_0,   \
                             ta1_y_zz_xxyz_1,   \
                             ta1_y_zz_xxyzz_0,  \
                             ta1_y_zz_xxyzz_1,  \
                             ta1_y_zz_xxzz_0,   \
                             ta1_y_zz_xxzz_1,   \
                             ta1_y_zz_xxzzz_0,  \
                             ta1_y_zz_xxzzz_1,  \
                             ta1_y_zz_xyyy_0,   \
                             ta1_y_zz_xyyy_1,   \
                             ta1_y_zz_xyyyy_0,  \
                             ta1_y_zz_xyyyy_1,  \
                             ta1_y_zz_xyyyz_0,  \
                             ta1_y_zz_xyyyz_1,  \
                             ta1_y_zz_xyyz_0,   \
                             ta1_y_zz_xyyz_1,   \
                             ta1_y_zz_xyyzz_0,  \
                             ta1_y_zz_xyyzz_1,  \
                             ta1_y_zz_xyzz_0,   \
                             ta1_y_zz_xyzz_1,   \
                             ta1_y_zz_xyzzz_0,  \
                             ta1_y_zz_xyzzz_1,  \
                             ta1_y_zz_xzzz_0,   \
                             ta1_y_zz_xzzz_1,   \
                             ta1_y_zz_xzzzz_0,  \
                             ta1_y_zz_xzzzz_1,  \
                             ta1_y_zz_yyyy_0,   \
                             ta1_y_zz_yyyy_1,   \
                             ta1_y_zz_yyyyy_0,  \
                             ta1_y_zz_yyyyy_1,  \
                             ta1_y_zz_yyyyz_0,  \
                             ta1_y_zz_yyyyz_1,  \
                             ta1_y_zz_yyyz_0,   \
                             ta1_y_zz_yyyz_1,   \
                             ta1_y_zz_yyyzz_0,  \
                             ta1_y_zz_yyyzz_1,  \
                             ta1_y_zz_yyzz_0,   \
                             ta1_y_zz_yyzz_1,   \
                             ta1_y_zz_yyzzz_0,  \
                             ta1_y_zz_yyzzz_1,  \
                             ta1_y_zz_yzzz_0,   \
                             ta1_y_zz_yzzz_1,   \
                             ta1_y_zz_yzzzz_0,  \
                             ta1_y_zz_yzzzz_1,  \
                             ta1_y_zz_zzzz_0,   \
                             ta1_y_zz_zzzz_1,   \
                             ta1_y_zz_zzzzz_0,  \
                             ta1_y_zz_zzzzz_1,  \
                             ta1_y_zzz_xxxxx_0, \
                             ta1_y_zzz_xxxxy_0, \
                             ta1_y_zzz_xxxxz_0, \
                             ta1_y_zzz_xxxyy_0, \
                             ta1_y_zzz_xxxyz_0, \
                             ta1_y_zzz_xxxzz_0, \
                             ta1_y_zzz_xxyyy_0, \
                             ta1_y_zzz_xxyyz_0, \
                             ta1_y_zzz_xxyzz_0, \
                             ta1_y_zzz_xxzzz_0, \
                             ta1_y_zzz_xyyyy_0, \
                             ta1_y_zzz_xyyyz_0, \
                             ta1_y_zzz_xyyzz_0, \
                             ta1_y_zzz_xyzzz_0, \
                             ta1_y_zzz_xzzzz_0, \
                             ta1_y_zzz_yyyyy_0, \
                             ta1_y_zzz_yyyyz_0, \
                             ta1_y_zzz_yyyzz_0, \
                             ta1_y_zzz_yyzzz_0, \
                             ta1_y_zzz_yzzzz_0, \
                             ta1_y_zzz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzz_xxxxx_0[i] =
            2.0 * ta1_y_z_xxxxx_0[i] * fe_0 - 2.0 * ta1_y_z_xxxxx_1[i] * fe_0 + ta1_y_zz_xxxxx_0[i] * pa_z[i] - ta1_y_zz_xxxxx_1[i] * pc_z[i];

        ta1_y_zzz_xxxxy_0[i] =
            2.0 * ta1_y_z_xxxxy_0[i] * fe_0 - 2.0 * ta1_y_z_xxxxy_1[i] * fe_0 + ta1_y_zz_xxxxy_0[i] * pa_z[i] - ta1_y_zz_xxxxy_1[i] * pc_z[i];

        ta1_y_zzz_xxxxz_0[i] = 2.0 * ta1_y_z_xxxxz_0[i] * fe_0 - 2.0 * ta1_y_z_xxxxz_1[i] * fe_0 + ta1_y_zz_xxxx_0[i] * fe_0 -
                               ta1_y_zz_xxxx_1[i] * fe_0 + ta1_y_zz_xxxxz_0[i] * pa_z[i] - ta1_y_zz_xxxxz_1[i] * pc_z[i];

        ta1_y_zzz_xxxyy_0[i] =
            2.0 * ta1_y_z_xxxyy_0[i] * fe_0 - 2.0 * ta1_y_z_xxxyy_1[i] * fe_0 + ta1_y_zz_xxxyy_0[i] * pa_z[i] - ta1_y_zz_xxxyy_1[i] * pc_z[i];

        ta1_y_zzz_xxxyz_0[i] = 2.0 * ta1_y_z_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_z_xxxyz_1[i] * fe_0 + ta1_y_zz_xxxy_0[i] * fe_0 -
                               ta1_y_zz_xxxy_1[i] * fe_0 + ta1_y_zz_xxxyz_0[i] * pa_z[i] - ta1_y_zz_xxxyz_1[i] * pc_z[i];

        ta1_y_zzz_xxxzz_0[i] = 2.0 * ta1_y_z_xxxzz_0[i] * fe_0 - 2.0 * ta1_y_z_xxxzz_1[i] * fe_0 + 2.0 * ta1_y_zz_xxxz_0[i] * fe_0 -
                               2.0 * ta1_y_zz_xxxz_1[i] * fe_0 + ta1_y_zz_xxxzz_0[i] * pa_z[i] - ta1_y_zz_xxxzz_1[i] * pc_z[i];

        ta1_y_zzz_xxyyy_0[i] =
            2.0 * ta1_y_z_xxyyy_0[i] * fe_0 - 2.0 * ta1_y_z_xxyyy_1[i] * fe_0 + ta1_y_zz_xxyyy_0[i] * pa_z[i] - ta1_y_zz_xxyyy_1[i] * pc_z[i];

        ta1_y_zzz_xxyyz_0[i] = 2.0 * ta1_y_z_xxyyz_0[i] * fe_0 - 2.0 * ta1_y_z_xxyyz_1[i] * fe_0 + ta1_y_zz_xxyy_0[i] * fe_0 -
                               ta1_y_zz_xxyy_1[i] * fe_0 + ta1_y_zz_xxyyz_0[i] * pa_z[i] - ta1_y_zz_xxyyz_1[i] * pc_z[i];

        ta1_y_zzz_xxyzz_0[i] = 2.0 * ta1_y_z_xxyzz_0[i] * fe_0 - 2.0 * ta1_y_z_xxyzz_1[i] * fe_0 + 2.0 * ta1_y_zz_xxyz_0[i] * fe_0 -
                               2.0 * ta1_y_zz_xxyz_1[i] * fe_0 + ta1_y_zz_xxyzz_0[i] * pa_z[i] - ta1_y_zz_xxyzz_1[i] * pc_z[i];

        ta1_y_zzz_xxzzz_0[i] = 2.0 * ta1_y_z_xxzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xxzzz_1[i] * fe_0 + 3.0 * ta1_y_zz_xxzz_0[i] * fe_0 -
                               3.0 * ta1_y_zz_xxzz_1[i] * fe_0 + ta1_y_zz_xxzzz_0[i] * pa_z[i] - ta1_y_zz_xxzzz_1[i] * pc_z[i];

        ta1_y_zzz_xyyyy_0[i] =
            2.0 * ta1_y_z_xyyyy_0[i] * fe_0 - 2.0 * ta1_y_z_xyyyy_1[i] * fe_0 + ta1_y_zz_xyyyy_0[i] * pa_z[i] - ta1_y_zz_xyyyy_1[i] * pc_z[i];

        ta1_y_zzz_xyyyz_0[i] = 2.0 * ta1_y_z_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_z_xyyyz_1[i] * fe_0 + ta1_y_zz_xyyy_0[i] * fe_0 -
                               ta1_y_zz_xyyy_1[i] * fe_0 + ta1_y_zz_xyyyz_0[i] * pa_z[i] - ta1_y_zz_xyyyz_1[i] * pc_z[i];

        ta1_y_zzz_xyyzz_0[i] = 2.0 * ta1_y_z_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_z_xyyzz_1[i] * fe_0 + 2.0 * ta1_y_zz_xyyz_0[i] * fe_0 -
                               2.0 * ta1_y_zz_xyyz_1[i] * fe_0 + ta1_y_zz_xyyzz_0[i] * pa_z[i] - ta1_y_zz_xyyzz_1[i] * pc_z[i];

        ta1_y_zzz_xyzzz_0[i] = 2.0 * ta1_y_z_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xyzzz_1[i] * fe_0 + 3.0 * ta1_y_zz_xyzz_0[i] * fe_0 -
                               3.0 * ta1_y_zz_xyzz_1[i] * fe_0 + ta1_y_zz_xyzzz_0[i] * pa_z[i] - ta1_y_zz_xyzzz_1[i] * pc_z[i];

        ta1_y_zzz_xzzzz_0[i] = 2.0 * ta1_y_z_xzzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xzzzz_1[i] * fe_0 + 4.0 * ta1_y_zz_xzzz_0[i] * fe_0 -
                               4.0 * ta1_y_zz_xzzz_1[i] * fe_0 + ta1_y_zz_xzzzz_0[i] * pa_z[i] - ta1_y_zz_xzzzz_1[i] * pc_z[i];

        ta1_y_zzz_yyyyy_0[i] =
            2.0 * ta1_y_z_yyyyy_0[i] * fe_0 - 2.0 * ta1_y_z_yyyyy_1[i] * fe_0 + ta1_y_zz_yyyyy_0[i] * pa_z[i] - ta1_y_zz_yyyyy_1[i] * pc_z[i];

        ta1_y_zzz_yyyyz_0[i] = 2.0 * ta1_y_z_yyyyz_0[i] * fe_0 - 2.0 * ta1_y_z_yyyyz_1[i] * fe_0 + ta1_y_zz_yyyy_0[i] * fe_0 -
                               ta1_y_zz_yyyy_1[i] * fe_0 + ta1_y_zz_yyyyz_0[i] * pa_z[i] - ta1_y_zz_yyyyz_1[i] * pc_z[i];

        ta1_y_zzz_yyyzz_0[i] = 2.0 * ta1_y_z_yyyzz_0[i] * fe_0 - 2.0 * ta1_y_z_yyyzz_1[i] * fe_0 + 2.0 * ta1_y_zz_yyyz_0[i] * fe_0 -
                               2.0 * ta1_y_zz_yyyz_1[i] * fe_0 + ta1_y_zz_yyyzz_0[i] * pa_z[i] - ta1_y_zz_yyyzz_1[i] * pc_z[i];

        ta1_y_zzz_yyzzz_0[i] = 2.0 * ta1_y_z_yyzzz_0[i] * fe_0 - 2.0 * ta1_y_z_yyzzz_1[i] * fe_0 + 3.0 * ta1_y_zz_yyzz_0[i] * fe_0 -
                               3.0 * ta1_y_zz_yyzz_1[i] * fe_0 + ta1_y_zz_yyzzz_0[i] * pa_z[i] - ta1_y_zz_yyzzz_1[i] * pc_z[i];

        ta1_y_zzz_yzzzz_0[i] = 2.0 * ta1_y_z_yzzzz_0[i] * fe_0 - 2.0 * ta1_y_z_yzzzz_1[i] * fe_0 + 4.0 * ta1_y_zz_yzzz_0[i] * fe_0 -
                               4.0 * ta1_y_zz_yzzz_1[i] * fe_0 + ta1_y_zz_yzzzz_0[i] * pa_z[i] - ta1_y_zz_yzzzz_1[i] * pc_z[i];

        ta1_y_zzz_zzzzz_0[i] = 2.0 * ta1_y_z_zzzzz_0[i] * fe_0 - 2.0 * ta1_y_z_zzzzz_1[i] * fe_0 + 5.0 * ta1_y_zz_zzzz_0[i] * fe_0 -
                               5.0 * ta1_y_zz_zzzz_1[i] * fe_0 + ta1_y_zz_zzzzz_0[i] * pa_z[i] - ta1_y_zz_zzzzz_1[i] * pc_z[i];
    }

    // Set up 420-441 components of targeted buffer : FH

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

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_z_x_xxxxx_0,   \
                             ta1_z_x_xxxxx_1,   \
                             ta1_z_x_xxxxy_0,   \
                             ta1_z_x_xxxxy_1,   \
                             ta1_z_x_xxxxz_0,   \
                             ta1_z_x_xxxxz_1,   \
                             ta1_z_x_xxxyy_0,   \
                             ta1_z_x_xxxyy_1,   \
                             ta1_z_x_xxxyz_0,   \
                             ta1_z_x_xxxyz_1,   \
                             ta1_z_x_xxxzz_0,   \
                             ta1_z_x_xxxzz_1,   \
                             ta1_z_x_xxyyy_0,   \
                             ta1_z_x_xxyyy_1,   \
                             ta1_z_x_xxyyz_0,   \
                             ta1_z_x_xxyyz_1,   \
                             ta1_z_x_xxyzz_0,   \
                             ta1_z_x_xxyzz_1,   \
                             ta1_z_x_xxzzz_0,   \
                             ta1_z_x_xxzzz_1,   \
                             ta1_z_x_xyyyy_0,   \
                             ta1_z_x_xyyyy_1,   \
                             ta1_z_x_xyyyz_0,   \
                             ta1_z_x_xyyyz_1,   \
                             ta1_z_x_xyyzz_0,   \
                             ta1_z_x_xyyzz_1,   \
                             ta1_z_x_xyzzz_0,   \
                             ta1_z_x_xyzzz_1,   \
                             ta1_z_x_xzzzz_0,   \
                             ta1_z_x_xzzzz_1,   \
                             ta1_z_x_yyyyy_0,   \
                             ta1_z_x_yyyyy_1,   \
                             ta1_z_x_yyyyz_0,   \
                             ta1_z_x_yyyyz_1,   \
                             ta1_z_x_yyyzz_0,   \
                             ta1_z_x_yyyzz_1,   \
                             ta1_z_x_yyzzz_0,   \
                             ta1_z_x_yyzzz_1,   \
                             ta1_z_x_yzzzz_0,   \
                             ta1_z_x_yzzzz_1,   \
                             ta1_z_x_zzzzz_0,   \
                             ta1_z_x_zzzzz_1,   \
                             ta1_z_xx_xxxx_0,   \
                             ta1_z_xx_xxxx_1,   \
                             ta1_z_xx_xxxxx_0,  \
                             ta1_z_xx_xxxxx_1,  \
                             ta1_z_xx_xxxxy_0,  \
                             ta1_z_xx_xxxxy_1,  \
                             ta1_z_xx_xxxxz_0,  \
                             ta1_z_xx_xxxxz_1,  \
                             ta1_z_xx_xxxy_0,   \
                             ta1_z_xx_xxxy_1,   \
                             ta1_z_xx_xxxyy_0,  \
                             ta1_z_xx_xxxyy_1,  \
                             ta1_z_xx_xxxyz_0,  \
                             ta1_z_xx_xxxyz_1,  \
                             ta1_z_xx_xxxz_0,   \
                             ta1_z_xx_xxxz_1,   \
                             ta1_z_xx_xxxzz_0,  \
                             ta1_z_xx_xxxzz_1,  \
                             ta1_z_xx_xxyy_0,   \
                             ta1_z_xx_xxyy_1,   \
                             ta1_z_xx_xxyyy_0,  \
                             ta1_z_xx_xxyyy_1,  \
                             ta1_z_xx_xxyyz_0,  \
                             ta1_z_xx_xxyyz_1,  \
                             ta1_z_xx_xxyz_0,   \
                             ta1_z_xx_xxyz_1,   \
                             ta1_z_xx_xxyzz_0,  \
                             ta1_z_xx_xxyzz_1,  \
                             ta1_z_xx_xxzz_0,   \
                             ta1_z_xx_xxzz_1,   \
                             ta1_z_xx_xxzzz_0,  \
                             ta1_z_xx_xxzzz_1,  \
                             ta1_z_xx_xyyy_0,   \
                             ta1_z_xx_xyyy_1,   \
                             ta1_z_xx_xyyyy_0,  \
                             ta1_z_xx_xyyyy_1,  \
                             ta1_z_xx_xyyyz_0,  \
                             ta1_z_xx_xyyyz_1,  \
                             ta1_z_xx_xyyz_0,   \
                             ta1_z_xx_xyyz_1,   \
                             ta1_z_xx_xyyzz_0,  \
                             ta1_z_xx_xyyzz_1,  \
                             ta1_z_xx_xyzz_0,   \
                             ta1_z_xx_xyzz_1,   \
                             ta1_z_xx_xyzzz_0,  \
                             ta1_z_xx_xyzzz_1,  \
                             ta1_z_xx_xzzz_0,   \
                             ta1_z_xx_xzzz_1,   \
                             ta1_z_xx_xzzzz_0,  \
                             ta1_z_xx_xzzzz_1,  \
                             ta1_z_xx_yyyy_0,   \
                             ta1_z_xx_yyyy_1,   \
                             ta1_z_xx_yyyyy_0,  \
                             ta1_z_xx_yyyyy_1,  \
                             ta1_z_xx_yyyyz_0,  \
                             ta1_z_xx_yyyyz_1,  \
                             ta1_z_xx_yyyz_0,   \
                             ta1_z_xx_yyyz_1,   \
                             ta1_z_xx_yyyzz_0,  \
                             ta1_z_xx_yyyzz_1,  \
                             ta1_z_xx_yyzz_0,   \
                             ta1_z_xx_yyzz_1,   \
                             ta1_z_xx_yyzzz_0,  \
                             ta1_z_xx_yyzzz_1,  \
                             ta1_z_xx_yzzz_0,   \
                             ta1_z_xx_yzzz_1,   \
                             ta1_z_xx_yzzzz_0,  \
                             ta1_z_xx_yzzzz_1,  \
                             ta1_z_xx_zzzz_0,   \
                             ta1_z_xx_zzzz_1,   \
                             ta1_z_xx_zzzzz_0,  \
                             ta1_z_xx_zzzzz_1,  \
                             ta1_z_xxx_xxxxx_0, \
                             ta1_z_xxx_xxxxy_0, \
                             ta1_z_xxx_xxxxz_0, \
                             ta1_z_xxx_xxxyy_0, \
                             ta1_z_xxx_xxxyz_0, \
                             ta1_z_xxx_xxxzz_0, \
                             ta1_z_xxx_xxyyy_0, \
                             ta1_z_xxx_xxyyz_0, \
                             ta1_z_xxx_xxyzz_0, \
                             ta1_z_xxx_xxzzz_0, \
                             ta1_z_xxx_xyyyy_0, \
                             ta1_z_xxx_xyyyz_0, \
                             ta1_z_xxx_xyyzz_0, \
                             ta1_z_xxx_xyzzz_0, \
                             ta1_z_xxx_xzzzz_0, \
                             ta1_z_xxx_yyyyy_0, \
                             ta1_z_xxx_yyyyz_0, \
                             ta1_z_xxx_yyyzz_0, \
                             ta1_z_xxx_yyzzz_0, \
                             ta1_z_xxx_yzzzz_0, \
                             ta1_z_xxx_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxx_xxxxx_0[i] = 2.0 * ta1_z_x_xxxxx_0[i] * fe_0 - 2.0 * ta1_z_x_xxxxx_1[i] * fe_0 + 5.0 * ta1_z_xx_xxxx_0[i] * fe_0 -
                               5.0 * ta1_z_xx_xxxx_1[i] * fe_0 + ta1_z_xx_xxxxx_0[i] * pa_x[i] - ta1_z_xx_xxxxx_1[i] * pc_x[i];

        ta1_z_xxx_xxxxy_0[i] = 2.0 * ta1_z_x_xxxxy_0[i] * fe_0 - 2.0 * ta1_z_x_xxxxy_1[i] * fe_0 + 4.0 * ta1_z_xx_xxxy_0[i] * fe_0 -
                               4.0 * ta1_z_xx_xxxy_1[i] * fe_0 + ta1_z_xx_xxxxy_0[i] * pa_x[i] - ta1_z_xx_xxxxy_1[i] * pc_x[i];

        ta1_z_xxx_xxxxz_0[i] = 2.0 * ta1_z_x_xxxxz_0[i] * fe_0 - 2.0 * ta1_z_x_xxxxz_1[i] * fe_0 + 4.0 * ta1_z_xx_xxxz_0[i] * fe_0 -
                               4.0 * ta1_z_xx_xxxz_1[i] * fe_0 + ta1_z_xx_xxxxz_0[i] * pa_x[i] - ta1_z_xx_xxxxz_1[i] * pc_x[i];

        ta1_z_xxx_xxxyy_0[i] = 2.0 * ta1_z_x_xxxyy_0[i] * fe_0 - 2.0 * ta1_z_x_xxxyy_1[i] * fe_0 + 3.0 * ta1_z_xx_xxyy_0[i] * fe_0 -
                               3.0 * ta1_z_xx_xxyy_1[i] * fe_0 + ta1_z_xx_xxxyy_0[i] * pa_x[i] - ta1_z_xx_xxxyy_1[i] * pc_x[i];

        ta1_z_xxx_xxxyz_0[i] = 2.0 * ta1_z_x_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_x_xxxyz_1[i] * fe_0 + 3.0 * ta1_z_xx_xxyz_0[i] * fe_0 -
                               3.0 * ta1_z_xx_xxyz_1[i] * fe_0 + ta1_z_xx_xxxyz_0[i] * pa_x[i] - ta1_z_xx_xxxyz_1[i] * pc_x[i];

        ta1_z_xxx_xxxzz_0[i] = 2.0 * ta1_z_x_xxxzz_0[i] * fe_0 - 2.0 * ta1_z_x_xxxzz_1[i] * fe_0 + 3.0 * ta1_z_xx_xxzz_0[i] * fe_0 -
                               3.0 * ta1_z_xx_xxzz_1[i] * fe_0 + ta1_z_xx_xxxzz_0[i] * pa_x[i] - ta1_z_xx_xxxzz_1[i] * pc_x[i];

        ta1_z_xxx_xxyyy_0[i] = 2.0 * ta1_z_x_xxyyy_0[i] * fe_0 - 2.0 * ta1_z_x_xxyyy_1[i] * fe_0 + 2.0 * ta1_z_xx_xyyy_0[i] * fe_0 -
                               2.0 * ta1_z_xx_xyyy_1[i] * fe_0 + ta1_z_xx_xxyyy_0[i] * pa_x[i] - ta1_z_xx_xxyyy_1[i] * pc_x[i];

        ta1_z_xxx_xxyyz_0[i] = 2.0 * ta1_z_x_xxyyz_0[i] * fe_0 - 2.0 * ta1_z_x_xxyyz_1[i] * fe_0 + 2.0 * ta1_z_xx_xyyz_0[i] * fe_0 -
                               2.0 * ta1_z_xx_xyyz_1[i] * fe_0 + ta1_z_xx_xxyyz_0[i] * pa_x[i] - ta1_z_xx_xxyyz_1[i] * pc_x[i];

        ta1_z_xxx_xxyzz_0[i] = 2.0 * ta1_z_x_xxyzz_0[i] * fe_0 - 2.0 * ta1_z_x_xxyzz_1[i] * fe_0 + 2.0 * ta1_z_xx_xyzz_0[i] * fe_0 -
                               2.0 * ta1_z_xx_xyzz_1[i] * fe_0 + ta1_z_xx_xxyzz_0[i] * pa_x[i] - ta1_z_xx_xxyzz_1[i] * pc_x[i];

        ta1_z_xxx_xxzzz_0[i] = 2.0 * ta1_z_x_xxzzz_0[i] * fe_0 - 2.0 * ta1_z_x_xxzzz_1[i] * fe_0 + 2.0 * ta1_z_xx_xzzz_0[i] * fe_0 -
                               2.0 * ta1_z_xx_xzzz_1[i] * fe_0 + ta1_z_xx_xxzzz_0[i] * pa_x[i] - ta1_z_xx_xxzzz_1[i] * pc_x[i];

        ta1_z_xxx_xyyyy_0[i] = 2.0 * ta1_z_x_xyyyy_0[i] * fe_0 - 2.0 * ta1_z_x_xyyyy_1[i] * fe_0 + ta1_z_xx_yyyy_0[i] * fe_0 -
                               ta1_z_xx_yyyy_1[i] * fe_0 + ta1_z_xx_xyyyy_0[i] * pa_x[i] - ta1_z_xx_xyyyy_1[i] * pc_x[i];

        ta1_z_xxx_xyyyz_0[i] = 2.0 * ta1_z_x_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_x_xyyyz_1[i] * fe_0 + ta1_z_xx_yyyz_0[i] * fe_0 -
                               ta1_z_xx_yyyz_1[i] * fe_0 + ta1_z_xx_xyyyz_0[i] * pa_x[i] - ta1_z_xx_xyyyz_1[i] * pc_x[i];

        ta1_z_xxx_xyyzz_0[i] = 2.0 * ta1_z_x_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_x_xyyzz_1[i] * fe_0 + ta1_z_xx_yyzz_0[i] * fe_0 -
                               ta1_z_xx_yyzz_1[i] * fe_0 + ta1_z_xx_xyyzz_0[i] * pa_x[i] - ta1_z_xx_xyyzz_1[i] * pc_x[i];

        ta1_z_xxx_xyzzz_0[i] = 2.0 * ta1_z_x_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_x_xyzzz_1[i] * fe_0 + ta1_z_xx_yzzz_0[i] * fe_0 -
                               ta1_z_xx_yzzz_1[i] * fe_0 + ta1_z_xx_xyzzz_0[i] * pa_x[i] - ta1_z_xx_xyzzz_1[i] * pc_x[i];

        ta1_z_xxx_xzzzz_0[i] = 2.0 * ta1_z_x_xzzzz_0[i] * fe_0 - 2.0 * ta1_z_x_xzzzz_1[i] * fe_0 + ta1_z_xx_zzzz_0[i] * fe_0 -
                               ta1_z_xx_zzzz_1[i] * fe_0 + ta1_z_xx_xzzzz_0[i] * pa_x[i] - ta1_z_xx_xzzzz_1[i] * pc_x[i];

        ta1_z_xxx_yyyyy_0[i] =
            2.0 * ta1_z_x_yyyyy_0[i] * fe_0 - 2.0 * ta1_z_x_yyyyy_1[i] * fe_0 + ta1_z_xx_yyyyy_0[i] * pa_x[i] - ta1_z_xx_yyyyy_1[i] * pc_x[i];

        ta1_z_xxx_yyyyz_0[i] =
            2.0 * ta1_z_x_yyyyz_0[i] * fe_0 - 2.0 * ta1_z_x_yyyyz_1[i] * fe_0 + ta1_z_xx_yyyyz_0[i] * pa_x[i] - ta1_z_xx_yyyyz_1[i] * pc_x[i];

        ta1_z_xxx_yyyzz_0[i] =
            2.0 * ta1_z_x_yyyzz_0[i] * fe_0 - 2.0 * ta1_z_x_yyyzz_1[i] * fe_0 + ta1_z_xx_yyyzz_0[i] * pa_x[i] - ta1_z_xx_yyyzz_1[i] * pc_x[i];

        ta1_z_xxx_yyzzz_0[i] =
            2.0 * ta1_z_x_yyzzz_0[i] * fe_0 - 2.0 * ta1_z_x_yyzzz_1[i] * fe_0 + ta1_z_xx_yyzzz_0[i] * pa_x[i] - ta1_z_xx_yyzzz_1[i] * pc_x[i];

        ta1_z_xxx_yzzzz_0[i] =
            2.0 * ta1_z_x_yzzzz_0[i] * fe_0 - 2.0 * ta1_z_x_yzzzz_1[i] * fe_0 + ta1_z_xx_yzzzz_0[i] * pa_x[i] - ta1_z_xx_yzzzz_1[i] * pc_x[i];

        ta1_z_xxx_zzzzz_0[i] =
            2.0 * ta1_z_x_zzzzz_0[i] * fe_0 - 2.0 * ta1_z_x_zzzzz_1[i] * fe_0 + ta1_z_xx_zzzzz_0[i] * pa_x[i] - ta1_z_xx_zzzzz_1[i] * pc_x[i];
    }

    // Set up 441-462 components of targeted buffer : FH

    auto ta1_z_xxy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 441);

    auto ta1_z_xxy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 442);

    auto ta1_z_xxy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 443);

    auto ta1_z_xxy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 444);

    auto ta1_z_xxy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 445);

    auto ta1_z_xxy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 446);

    auto ta1_z_xxy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 447);

    auto ta1_z_xxy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 448);

    auto ta1_z_xxy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 449);

    auto ta1_z_xxy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 450);

    auto ta1_z_xxy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 451);

    auto ta1_z_xxy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 452);

    auto ta1_z_xxy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 453);

    auto ta1_z_xxy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 454);

    auto ta1_z_xxy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 455);

    auto ta1_z_xxy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 456);

    auto ta1_z_xxy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 457);

    auto ta1_z_xxy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 458);

    auto ta1_z_xxy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 459);

    auto ta1_z_xxy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 460);

    auto ta1_z_xxy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 461);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_z_xx_xxxx_0,   \
                             ta1_z_xx_xxxx_1,   \
                             ta1_z_xx_xxxxx_0,  \
                             ta1_z_xx_xxxxx_1,  \
                             ta1_z_xx_xxxxy_0,  \
                             ta1_z_xx_xxxxy_1,  \
                             ta1_z_xx_xxxxz_0,  \
                             ta1_z_xx_xxxxz_1,  \
                             ta1_z_xx_xxxy_0,   \
                             ta1_z_xx_xxxy_1,   \
                             ta1_z_xx_xxxyy_0,  \
                             ta1_z_xx_xxxyy_1,  \
                             ta1_z_xx_xxxyz_0,  \
                             ta1_z_xx_xxxyz_1,  \
                             ta1_z_xx_xxxz_0,   \
                             ta1_z_xx_xxxz_1,   \
                             ta1_z_xx_xxxzz_0,  \
                             ta1_z_xx_xxxzz_1,  \
                             ta1_z_xx_xxyy_0,   \
                             ta1_z_xx_xxyy_1,   \
                             ta1_z_xx_xxyyy_0,  \
                             ta1_z_xx_xxyyy_1,  \
                             ta1_z_xx_xxyyz_0,  \
                             ta1_z_xx_xxyyz_1,  \
                             ta1_z_xx_xxyz_0,   \
                             ta1_z_xx_xxyz_1,   \
                             ta1_z_xx_xxyzz_0,  \
                             ta1_z_xx_xxyzz_1,  \
                             ta1_z_xx_xxzz_0,   \
                             ta1_z_xx_xxzz_1,   \
                             ta1_z_xx_xxzzz_0,  \
                             ta1_z_xx_xxzzz_1,  \
                             ta1_z_xx_xyyy_0,   \
                             ta1_z_xx_xyyy_1,   \
                             ta1_z_xx_xyyyy_0,  \
                             ta1_z_xx_xyyyy_1,  \
                             ta1_z_xx_xyyyz_0,  \
                             ta1_z_xx_xyyyz_1,  \
                             ta1_z_xx_xyyz_0,   \
                             ta1_z_xx_xyyz_1,   \
                             ta1_z_xx_xyyzz_0,  \
                             ta1_z_xx_xyyzz_1,  \
                             ta1_z_xx_xyzz_0,   \
                             ta1_z_xx_xyzz_1,   \
                             ta1_z_xx_xyzzz_0,  \
                             ta1_z_xx_xyzzz_1,  \
                             ta1_z_xx_xzzz_0,   \
                             ta1_z_xx_xzzz_1,   \
                             ta1_z_xx_xzzzz_0,  \
                             ta1_z_xx_xzzzz_1,  \
                             ta1_z_xx_zzzzz_0,  \
                             ta1_z_xx_zzzzz_1,  \
                             ta1_z_xxy_xxxxx_0, \
                             ta1_z_xxy_xxxxy_0, \
                             ta1_z_xxy_xxxxz_0, \
                             ta1_z_xxy_xxxyy_0, \
                             ta1_z_xxy_xxxyz_0, \
                             ta1_z_xxy_xxxzz_0, \
                             ta1_z_xxy_xxyyy_0, \
                             ta1_z_xxy_xxyyz_0, \
                             ta1_z_xxy_xxyzz_0, \
                             ta1_z_xxy_xxzzz_0, \
                             ta1_z_xxy_xyyyy_0, \
                             ta1_z_xxy_xyyyz_0, \
                             ta1_z_xxy_xyyzz_0, \
                             ta1_z_xxy_xyzzz_0, \
                             ta1_z_xxy_xzzzz_0, \
                             ta1_z_xxy_yyyyy_0, \
                             ta1_z_xxy_yyyyz_0, \
                             ta1_z_xxy_yyyzz_0, \
                             ta1_z_xxy_yyzzz_0, \
                             ta1_z_xxy_yzzzz_0, \
                             ta1_z_xxy_zzzzz_0, \
                             ta1_z_xy_yyyyy_0,  \
                             ta1_z_xy_yyyyy_1,  \
                             ta1_z_xy_yyyyz_0,  \
                             ta1_z_xy_yyyyz_1,  \
                             ta1_z_xy_yyyzz_0,  \
                             ta1_z_xy_yyyzz_1,  \
                             ta1_z_xy_yyzzz_0,  \
                             ta1_z_xy_yyzzz_1,  \
                             ta1_z_xy_yzzzz_0,  \
                             ta1_z_xy_yzzzz_1,  \
                             ta1_z_y_yyyyy_0,   \
                             ta1_z_y_yyyyy_1,   \
                             ta1_z_y_yyyyz_0,   \
                             ta1_z_y_yyyyz_1,   \
                             ta1_z_y_yyyzz_0,   \
                             ta1_z_y_yyyzz_1,   \
                             ta1_z_y_yyzzz_0,   \
                             ta1_z_y_yyzzz_1,   \
                             ta1_z_y_yzzzz_0,   \
                             ta1_z_y_yzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxy_xxxxx_0[i] = ta1_z_xx_xxxxx_0[i] * pa_y[i] - ta1_z_xx_xxxxx_1[i] * pc_y[i];

        ta1_z_xxy_xxxxy_0[i] = ta1_z_xx_xxxx_0[i] * fe_0 - ta1_z_xx_xxxx_1[i] * fe_0 + ta1_z_xx_xxxxy_0[i] * pa_y[i] - ta1_z_xx_xxxxy_1[i] * pc_y[i];

        ta1_z_xxy_xxxxz_0[i] = ta1_z_xx_xxxxz_0[i] * pa_y[i] - ta1_z_xx_xxxxz_1[i] * pc_y[i];

        ta1_z_xxy_xxxyy_0[i] =
            2.0 * ta1_z_xx_xxxy_0[i] * fe_0 - 2.0 * ta1_z_xx_xxxy_1[i] * fe_0 + ta1_z_xx_xxxyy_0[i] * pa_y[i] - ta1_z_xx_xxxyy_1[i] * pc_y[i];

        ta1_z_xxy_xxxyz_0[i] = ta1_z_xx_xxxz_0[i] * fe_0 - ta1_z_xx_xxxz_1[i] * fe_0 + ta1_z_xx_xxxyz_0[i] * pa_y[i] - ta1_z_xx_xxxyz_1[i] * pc_y[i];

        ta1_z_xxy_xxxzz_0[i] = ta1_z_xx_xxxzz_0[i] * pa_y[i] - ta1_z_xx_xxxzz_1[i] * pc_y[i];

        ta1_z_xxy_xxyyy_0[i] =
            3.0 * ta1_z_xx_xxyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxyy_1[i] * fe_0 + ta1_z_xx_xxyyy_0[i] * pa_y[i] - ta1_z_xx_xxyyy_1[i] * pc_y[i];

        ta1_z_xxy_xxyyz_0[i] =
            2.0 * ta1_z_xx_xxyz_0[i] * fe_0 - 2.0 * ta1_z_xx_xxyz_1[i] * fe_0 + ta1_z_xx_xxyyz_0[i] * pa_y[i] - ta1_z_xx_xxyyz_1[i] * pc_y[i];

        ta1_z_xxy_xxyzz_0[i] = ta1_z_xx_xxzz_0[i] * fe_0 - ta1_z_xx_xxzz_1[i] * fe_0 + ta1_z_xx_xxyzz_0[i] * pa_y[i] - ta1_z_xx_xxyzz_1[i] * pc_y[i];

        ta1_z_xxy_xxzzz_0[i] = ta1_z_xx_xxzzz_0[i] * pa_y[i] - ta1_z_xx_xxzzz_1[i] * pc_y[i];

        ta1_z_xxy_xyyyy_0[i] =
            4.0 * ta1_z_xx_xyyy_0[i] * fe_0 - 4.0 * ta1_z_xx_xyyy_1[i] * fe_0 + ta1_z_xx_xyyyy_0[i] * pa_y[i] - ta1_z_xx_xyyyy_1[i] * pc_y[i];

        ta1_z_xxy_xyyyz_0[i] =
            3.0 * ta1_z_xx_xyyz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyyz_1[i] * fe_0 + ta1_z_xx_xyyyz_0[i] * pa_y[i] - ta1_z_xx_xyyyz_1[i] * pc_y[i];

        ta1_z_xxy_xyyzz_0[i] =
            2.0 * ta1_z_xx_xyzz_0[i] * fe_0 - 2.0 * ta1_z_xx_xyzz_1[i] * fe_0 + ta1_z_xx_xyyzz_0[i] * pa_y[i] - ta1_z_xx_xyyzz_1[i] * pc_y[i];

        ta1_z_xxy_xyzzz_0[i] = ta1_z_xx_xzzz_0[i] * fe_0 - ta1_z_xx_xzzz_1[i] * fe_0 + ta1_z_xx_xyzzz_0[i] * pa_y[i] - ta1_z_xx_xyzzz_1[i] * pc_y[i];

        ta1_z_xxy_xzzzz_0[i] = ta1_z_xx_xzzzz_0[i] * pa_y[i] - ta1_z_xx_xzzzz_1[i] * pc_y[i];

        ta1_z_xxy_yyyyy_0[i] = ta1_z_y_yyyyy_0[i] * fe_0 - ta1_z_y_yyyyy_1[i] * fe_0 + ta1_z_xy_yyyyy_0[i] * pa_x[i] - ta1_z_xy_yyyyy_1[i] * pc_x[i];

        ta1_z_xxy_yyyyz_0[i] = ta1_z_y_yyyyz_0[i] * fe_0 - ta1_z_y_yyyyz_1[i] * fe_0 + ta1_z_xy_yyyyz_0[i] * pa_x[i] - ta1_z_xy_yyyyz_1[i] * pc_x[i];

        ta1_z_xxy_yyyzz_0[i] = ta1_z_y_yyyzz_0[i] * fe_0 - ta1_z_y_yyyzz_1[i] * fe_0 + ta1_z_xy_yyyzz_0[i] * pa_x[i] - ta1_z_xy_yyyzz_1[i] * pc_x[i];

        ta1_z_xxy_yyzzz_0[i] = ta1_z_y_yyzzz_0[i] * fe_0 - ta1_z_y_yyzzz_1[i] * fe_0 + ta1_z_xy_yyzzz_0[i] * pa_x[i] - ta1_z_xy_yyzzz_1[i] * pc_x[i];

        ta1_z_xxy_yzzzz_0[i] = ta1_z_y_yzzzz_0[i] * fe_0 - ta1_z_y_yzzzz_1[i] * fe_0 + ta1_z_xy_yzzzz_0[i] * pa_x[i] - ta1_z_xy_yzzzz_1[i] * pc_x[i];

        ta1_z_xxy_zzzzz_0[i] = ta1_z_xx_zzzzz_0[i] * pa_y[i] - ta1_z_xx_zzzzz_1[i] * pc_y[i];
    }

    // Set up 462-483 components of targeted buffer : FH

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

    auto ta1_z_xxz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 477);

    auto ta1_z_xxz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 478);

    auto ta1_z_xxz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 479);

    auto ta1_z_xxz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 480);

    auto ta1_z_xxz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 481);

    auto ta1_z_xxz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 482);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_z_xx_xxxx_0,   \
                             ta1_z_xx_xxxx_1,   \
                             ta1_z_xx_xxxxx_0,  \
                             ta1_z_xx_xxxxx_1,  \
                             ta1_z_xx_xxxxy_0,  \
                             ta1_z_xx_xxxxy_1,  \
                             ta1_z_xx_xxxxz_0,  \
                             ta1_z_xx_xxxxz_1,  \
                             ta1_z_xx_xxxy_0,   \
                             ta1_z_xx_xxxy_1,   \
                             ta1_z_xx_xxxyy_0,  \
                             ta1_z_xx_xxxyy_1,  \
                             ta1_z_xx_xxxyz_0,  \
                             ta1_z_xx_xxxyz_1,  \
                             ta1_z_xx_xxxz_0,   \
                             ta1_z_xx_xxxz_1,   \
                             ta1_z_xx_xxxzz_0,  \
                             ta1_z_xx_xxxzz_1,  \
                             ta1_z_xx_xxyy_0,   \
                             ta1_z_xx_xxyy_1,   \
                             ta1_z_xx_xxyyy_0,  \
                             ta1_z_xx_xxyyy_1,  \
                             ta1_z_xx_xxyyz_0,  \
                             ta1_z_xx_xxyyz_1,  \
                             ta1_z_xx_xxyz_0,   \
                             ta1_z_xx_xxyz_1,   \
                             ta1_z_xx_xxyzz_0,  \
                             ta1_z_xx_xxyzz_1,  \
                             ta1_z_xx_xxzz_0,   \
                             ta1_z_xx_xxzz_1,   \
                             ta1_z_xx_xxzzz_0,  \
                             ta1_z_xx_xxzzz_1,  \
                             ta1_z_xx_xyyy_0,   \
                             ta1_z_xx_xyyy_1,   \
                             ta1_z_xx_xyyyy_0,  \
                             ta1_z_xx_xyyyy_1,  \
                             ta1_z_xx_xyyyz_0,  \
                             ta1_z_xx_xyyyz_1,  \
                             ta1_z_xx_xyyz_0,   \
                             ta1_z_xx_xyyz_1,   \
                             ta1_z_xx_xyyzz_0,  \
                             ta1_z_xx_xyyzz_1,  \
                             ta1_z_xx_xyzz_0,   \
                             ta1_z_xx_xyzz_1,   \
                             ta1_z_xx_xyzzz_0,  \
                             ta1_z_xx_xyzzz_1,  \
                             ta1_z_xx_xzzz_0,   \
                             ta1_z_xx_xzzz_1,   \
                             ta1_z_xx_xzzzz_0,  \
                             ta1_z_xx_xzzzz_1,  \
                             ta1_z_xx_yyyyy_0,  \
                             ta1_z_xx_yyyyy_1,  \
                             ta1_z_xxz_xxxxx_0, \
                             ta1_z_xxz_xxxxy_0, \
                             ta1_z_xxz_xxxxz_0, \
                             ta1_z_xxz_xxxyy_0, \
                             ta1_z_xxz_xxxyz_0, \
                             ta1_z_xxz_xxxzz_0, \
                             ta1_z_xxz_xxyyy_0, \
                             ta1_z_xxz_xxyyz_0, \
                             ta1_z_xxz_xxyzz_0, \
                             ta1_z_xxz_xxzzz_0, \
                             ta1_z_xxz_xyyyy_0, \
                             ta1_z_xxz_xyyyz_0, \
                             ta1_z_xxz_xyyzz_0, \
                             ta1_z_xxz_xyzzz_0, \
                             ta1_z_xxz_xzzzz_0, \
                             ta1_z_xxz_yyyyy_0, \
                             ta1_z_xxz_yyyyz_0, \
                             ta1_z_xxz_yyyzz_0, \
                             ta1_z_xxz_yyzzz_0, \
                             ta1_z_xxz_yzzzz_0, \
                             ta1_z_xxz_zzzzz_0, \
                             ta1_z_xz_yyyyz_0,  \
                             ta1_z_xz_yyyyz_1,  \
                             ta1_z_xz_yyyzz_0,  \
                             ta1_z_xz_yyyzz_1,  \
                             ta1_z_xz_yyzzz_0,  \
                             ta1_z_xz_yyzzz_1,  \
                             ta1_z_xz_yzzzz_0,  \
                             ta1_z_xz_yzzzz_1,  \
                             ta1_z_xz_zzzzz_0,  \
                             ta1_z_xz_zzzzz_1,  \
                             ta1_z_z_yyyyz_0,   \
                             ta1_z_z_yyyyz_1,   \
                             ta1_z_z_yyyzz_0,   \
                             ta1_z_z_yyyzz_1,   \
                             ta1_z_z_yyzzz_0,   \
                             ta1_z_z_yyzzz_1,   \
                             ta1_z_z_yzzzz_0,   \
                             ta1_z_z_yzzzz_1,   \
                             ta1_z_z_zzzzz_0,   \
                             ta1_z_z_zzzzz_1,   \
                             ta_xx_xxxxx_1,     \
                             ta_xx_xxxxy_1,     \
                             ta_xx_xxxxz_1,     \
                             ta_xx_xxxyy_1,     \
                             ta_xx_xxxyz_1,     \
                             ta_xx_xxxzz_1,     \
                             ta_xx_xxyyy_1,     \
                             ta_xx_xxyyz_1,     \
                             ta_xx_xxyzz_1,     \
                             ta_xx_xxzzz_1,     \
                             ta_xx_xyyyy_1,     \
                             ta_xx_xyyyz_1,     \
                             ta_xx_xyyzz_1,     \
                             ta_xx_xyzzz_1,     \
                             ta_xx_xzzzz_1,     \
                             ta_xx_yyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxz_xxxxx_0[i] = ta_xx_xxxxx_1[i] + ta1_z_xx_xxxxx_0[i] * pa_z[i] - ta1_z_xx_xxxxx_1[i] * pc_z[i];

        ta1_z_xxz_xxxxy_0[i] = ta_xx_xxxxy_1[i] + ta1_z_xx_xxxxy_0[i] * pa_z[i] - ta1_z_xx_xxxxy_1[i] * pc_z[i];

        ta1_z_xxz_xxxxz_0[i] =
            ta1_z_xx_xxxx_0[i] * fe_0 - ta1_z_xx_xxxx_1[i] * fe_0 + ta_xx_xxxxz_1[i] + ta1_z_xx_xxxxz_0[i] * pa_z[i] - ta1_z_xx_xxxxz_1[i] * pc_z[i];

        ta1_z_xxz_xxxyy_0[i] = ta_xx_xxxyy_1[i] + ta1_z_xx_xxxyy_0[i] * pa_z[i] - ta1_z_xx_xxxyy_1[i] * pc_z[i];

        ta1_z_xxz_xxxyz_0[i] =
            ta1_z_xx_xxxy_0[i] * fe_0 - ta1_z_xx_xxxy_1[i] * fe_0 + ta_xx_xxxyz_1[i] + ta1_z_xx_xxxyz_0[i] * pa_z[i] - ta1_z_xx_xxxyz_1[i] * pc_z[i];

        ta1_z_xxz_xxxzz_0[i] = 2.0 * ta1_z_xx_xxxz_0[i] * fe_0 - 2.0 * ta1_z_xx_xxxz_1[i] * fe_0 + ta_xx_xxxzz_1[i] + ta1_z_xx_xxxzz_0[i] * pa_z[i] -
                               ta1_z_xx_xxxzz_1[i] * pc_z[i];

        ta1_z_xxz_xxyyy_0[i] = ta_xx_xxyyy_1[i] + ta1_z_xx_xxyyy_0[i] * pa_z[i] - ta1_z_xx_xxyyy_1[i] * pc_z[i];

        ta1_z_xxz_xxyyz_0[i] =
            ta1_z_xx_xxyy_0[i] * fe_0 - ta1_z_xx_xxyy_1[i] * fe_0 + ta_xx_xxyyz_1[i] + ta1_z_xx_xxyyz_0[i] * pa_z[i] - ta1_z_xx_xxyyz_1[i] * pc_z[i];

        ta1_z_xxz_xxyzz_0[i] = 2.0 * ta1_z_xx_xxyz_0[i] * fe_0 - 2.0 * ta1_z_xx_xxyz_1[i] * fe_0 + ta_xx_xxyzz_1[i] + ta1_z_xx_xxyzz_0[i] * pa_z[i] -
                               ta1_z_xx_xxyzz_1[i] * pc_z[i];

        ta1_z_xxz_xxzzz_0[i] = 3.0 * ta1_z_xx_xxzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxzz_1[i] * fe_0 + ta_xx_xxzzz_1[i] + ta1_z_xx_xxzzz_0[i] * pa_z[i] -
                               ta1_z_xx_xxzzz_1[i] * pc_z[i];

        ta1_z_xxz_xyyyy_0[i] = ta_xx_xyyyy_1[i] + ta1_z_xx_xyyyy_0[i] * pa_z[i] - ta1_z_xx_xyyyy_1[i] * pc_z[i];

        ta1_z_xxz_xyyyz_0[i] =
            ta1_z_xx_xyyy_0[i] * fe_0 - ta1_z_xx_xyyy_1[i] * fe_0 + ta_xx_xyyyz_1[i] + ta1_z_xx_xyyyz_0[i] * pa_z[i] - ta1_z_xx_xyyyz_1[i] * pc_z[i];

        ta1_z_xxz_xyyzz_0[i] = 2.0 * ta1_z_xx_xyyz_0[i] * fe_0 - 2.0 * ta1_z_xx_xyyz_1[i] * fe_0 + ta_xx_xyyzz_1[i] + ta1_z_xx_xyyzz_0[i] * pa_z[i] -
                               ta1_z_xx_xyyzz_1[i] * pc_z[i];

        ta1_z_xxz_xyzzz_0[i] = 3.0 * ta1_z_xx_xyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyzz_1[i] * fe_0 + ta_xx_xyzzz_1[i] + ta1_z_xx_xyzzz_0[i] * pa_z[i] -
                               ta1_z_xx_xyzzz_1[i] * pc_z[i];

        ta1_z_xxz_xzzzz_0[i] = 4.0 * ta1_z_xx_xzzz_0[i] * fe_0 - 4.0 * ta1_z_xx_xzzz_1[i] * fe_0 + ta_xx_xzzzz_1[i] + ta1_z_xx_xzzzz_0[i] * pa_z[i] -
                               ta1_z_xx_xzzzz_1[i] * pc_z[i];

        ta1_z_xxz_yyyyy_0[i] = ta_xx_yyyyy_1[i] + ta1_z_xx_yyyyy_0[i] * pa_z[i] - ta1_z_xx_yyyyy_1[i] * pc_z[i];

        ta1_z_xxz_yyyyz_0[i] = ta1_z_z_yyyyz_0[i] * fe_0 - ta1_z_z_yyyyz_1[i] * fe_0 + ta1_z_xz_yyyyz_0[i] * pa_x[i] - ta1_z_xz_yyyyz_1[i] * pc_x[i];

        ta1_z_xxz_yyyzz_0[i] = ta1_z_z_yyyzz_0[i] * fe_0 - ta1_z_z_yyyzz_1[i] * fe_0 + ta1_z_xz_yyyzz_0[i] * pa_x[i] - ta1_z_xz_yyyzz_1[i] * pc_x[i];

        ta1_z_xxz_yyzzz_0[i] = ta1_z_z_yyzzz_0[i] * fe_0 - ta1_z_z_yyzzz_1[i] * fe_0 + ta1_z_xz_yyzzz_0[i] * pa_x[i] - ta1_z_xz_yyzzz_1[i] * pc_x[i];

        ta1_z_xxz_yzzzz_0[i] = ta1_z_z_yzzzz_0[i] * fe_0 - ta1_z_z_yzzzz_1[i] * fe_0 + ta1_z_xz_yzzzz_0[i] * pa_x[i] - ta1_z_xz_yzzzz_1[i] * pc_x[i];

        ta1_z_xxz_zzzzz_0[i] = ta1_z_z_zzzzz_0[i] * fe_0 - ta1_z_z_zzzzz_1[i] * fe_0 + ta1_z_xz_zzzzz_0[i] * pa_x[i] - ta1_z_xz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 483-504 components of targeted buffer : FH

    auto ta1_z_xyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 483);

    auto ta1_z_xyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 484);

    auto ta1_z_xyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 485);

    auto ta1_z_xyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 486);

    auto ta1_z_xyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 487);

    auto ta1_z_xyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 488);

    auto ta1_z_xyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 489);

    auto ta1_z_xyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 490);

    auto ta1_z_xyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 491);

    auto ta1_z_xyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 492);

    auto ta1_z_xyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 493);

    auto ta1_z_xyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 494);

    auto ta1_z_xyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 495);

    auto ta1_z_xyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 496);

    auto ta1_z_xyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 497);

    auto ta1_z_xyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 498);

    auto ta1_z_xyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 499);

    auto ta1_z_xyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 500);

    auto ta1_z_xyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 501);

    auto ta1_z_xyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 502);

    auto ta1_z_xyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 503);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_z_xyy_xxxxx_0, \
                             ta1_z_xyy_xxxxy_0, \
                             ta1_z_xyy_xxxxz_0, \
                             ta1_z_xyy_xxxyy_0, \
                             ta1_z_xyy_xxxyz_0, \
                             ta1_z_xyy_xxxzz_0, \
                             ta1_z_xyy_xxyyy_0, \
                             ta1_z_xyy_xxyyz_0, \
                             ta1_z_xyy_xxyzz_0, \
                             ta1_z_xyy_xxzzz_0, \
                             ta1_z_xyy_xyyyy_0, \
                             ta1_z_xyy_xyyyz_0, \
                             ta1_z_xyy_xyyzz_0, \
                             ta1_z_xyy_xyzzz_0, \
                             ta1_z_xyy_xzzzz_0, \
                             ta1_z_xyy_yyyyy_0, \
                             ta1_z_xyy_yyyyz_0, \
                             ta1_z_xyy_yyyzz_0, \
                             ta1_z_xyy_yyzzz_0, \
                             ta1_z_xyy_yzzzz_0, \
                             ta1_z_xyy_zzzzz_0, \
                             ta1_z_yy_xxxx_0,   \
                             ta1_z_yy_xxxx_1,   \
                             ta1_z_yy_xxxxx_0,  \
                             ta1_z_yy_xxxxx_1,  \
                             ta1_z_yy_xxxxy_0,  \
                             ta1_z_yy_xxxxy_1,  \
                             ta1_z_yy_xxxxz_0,  \
                             ta1_z_yy_xxxxz_1,  \
                             ta1_z_yy_xxxy_0,   \
                             ta1_z_yy_xxxy_1,   \
                             ta1_z_yy_xxxyy_0,  \
                             ta1_z_yy_xxxyy_1,  \
                             ta1_z_yy_xxxyz_0,  \
                             ta1_z_yy_xxxyz_1,  \
                             ta1_z_yy_xxxz_0,   \
                             ta1_z_yy_xxxz_1,   \
                             ta1_z_yy_xxxzz_0,  \
                             ta1_z_yy_xxxzz_1,  \
                             ta1_z_yy_xxyy_0,   \
                             ta1_z_yy_xxyy_1,   \
                             ta1_z_yy_xxyyy_0,  \
                             ta1_z_yy_xxyyy_1,  \
                             ta1_z_yy_xxyyz_0,  \
                             ta1_z_yy_xxyyz_1,  \
                             ta1_z_yy_xxyz_0,   \
                             ta1_z_yy_xxyz_1,   \
                             ta1_z_yy_xxyzz_0,  \
                             ta1_z_yy_xxyzz_1,  \
                             ta1_z_yy_xxzz_0,   \
                             ta1_z_yy_xxzz_1,   \
                             ta1_z_yy_xxzzz_0,  \
                             ta1_z_yy_xxzzz_1,  \
                             ta1_z_yy_xyyy_0,   \
                             ta1_z_yy_xyyy_1,   \
                             ta1_z_yy_xyyyy_0,  \
                             ta1_z_yy_xyyyy_1,  \
                             ta1_z_yy_xyyyz_0,  \
                             ta1_z_yy_xyyyz_1,  \
                             ta1_z_yy_xyyz_0,   \
                             ta1_z_yy_xyyz_1,   \
                             ta1_z_yy_xyyzz_0,  \
                             ta1_z_yy_xyyzz_1,  \
                             ta1_z_yy_xyzz_0,   \
                             ta1_z_yy_xyzz_1,   \
                             ta1_z_yy_xyzzz_0,  \
                             ta1_z_yy_xyzzz_1,  \
                             ta1_z_yy_xzzz_0,   \
                             ta1_z_yy_xzzz_1,   \
                             ta1_z_yy_xzzzz_0,  \
                             ta1_z_yy_xzzzz_1,  \
                             ta1_z_yy_yyyy_0,   \
                             ta1_z_yy_yyyy_1,   \
                             ta1_z_yy_yyyyy_0,  \
                             ta1_z_yy_yyyyy_1,  \
                             ta1_z_yy_yyyyz_0,  \
                             ta1_z_yy_yyyyz_1,  \
                             ta1_z_yy_yyyz_0,   \
                             ta1_z_yy_yyyz_1,   \
                             ta1_z_yy_yyyzz_0,  \
                             ta1_z_yy_yyyzz_1,  \
                             ta1_z_yy_yyzz_0,   \
                             ta1_z_yy_yyzz_1,   \
                             ta1_z_yy_yyzzz_0,  \
                             ta1_z_yy_yyzzz_1,  \
                             ta1_z_yy_yzzz_0,   \
                             ta1_z_yy_yzzz_1,   \
                             ta1_z_yy_yzzzz_0,  \
                             ta1_z_yy_yzzzz_1,  \
                             ta1_z_yy_zzzz_0,   \
                             ta1_z_yy_zzzz_1,   \
                             ta1_z_yy_zzzzz_0,  \
                             ta1_z_yy_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyy_xxxxx_0[i] =
            5.0 * ta1_z_yy_xxxx_0[i] * fe_0 - 5.0 * ta1_z_yy_xxxx_1[i] * fe_0 + ta1_z_yy_xxxxx_0[i] * pa_x[i] - ta1_z_yy_xxxxx_1[i] * pc_x[i];

        ta1_z_xyy_xxxxy_0[i] =
            4.0 * ta1_z_yy_xxxy_0[i] * fe_0 - 4.0 * ta1_z_yy_xxxy_1[i] * fe_0 + ta1_z_yy_xxxxy_0[i] * pa_x[i] - ta1_z_yy_xxxxy_1[i] * pc_x[i];

        ta1_z_xyy_xxxxz_0[i] =
            4.0 * ta1_z_yy_xxxz_0[i] * fe_0 - 4.0 * ta1_z_yy_xxxz_1[i] * fe_0 + ta1_z_yy_xxxxz_0[i] * pa_x[i] - ta1_z_yy_xxxxz_1[i] * pc_x[i];

        ta1_z_xyy_xxxyy_0[i] =
            3.0 * ta1_z_yy_xxyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyy_1[i] * fe_0 + ta1_z_yy_xxxyy_0[i] * pa_x[i] - ta1_z_yy_xxxyy_1[i] * pc_x[i];

        ta1_z_xyy_xxxyz_0[i] =
            3.0 * ta1_z_yy_xxyz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyz_1[i] * fe_0 + ta1_z_yy_xxxyz_0[i] * pa_x[i] - ta1_z_yy_xxxyz_1[i] * pc_x[i];

        ta1_z_xyy_xxxzz_0[i] =
            3.0 * ta1_z_yy_xxzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxzz_1[i] * fe_0 + ta1_z_yy_xxxzz_0[i] * pa_x[i] - ta1_z_yy_xxxzz_1[i] * pc_x[i];

        ta1_z_xyy_xxyyy_0[i] =
            2.0 * ta1_z_yy_xyyy_0[i] * fe_0 - 2.0 * ta1_z_yy_xyyy_1[i] * fe_0 + ta1_z_yy_xxyyy_0[i] * pa_x[i] - ta1_z_yy_xxyyy_1[i] * pc_x[i];

        ta1_z_xyy_xxyyz_0[i] =
            2.0 * ta1_z_yy_xyyz_0[i] * fe_0 - 2.0 * ta1_z_yy_xyyz_1[i] * fe_0 + ta1_z_yy_xxyyz_0[i] * pa_x[i] - ta1_z_yy_xxyyz_1[i] * pc_x[i];

        ta1_z_xyy_xxyzz_0[i] =
            2.0 * ta1_z_yy_xyzz_0[i] * fe_0 - 2.0 * ta1_z_yy_xyzz_1[i] * fe_0 + ta1_z_yy_xxyzz_0[i] * pa_x[i] - ta1_z_yy_xxyzz_1[i] * pc_x[i];

        ta1_z_xyy_xxzzz_0[i] =
            2.0 * ta1_z_yy_xzzz_0[i] * fe_0 - 2.0 * ta1_z_yy_xzzz_1[i] * fe_0 + ta1_z_yy_xxzzz_0[i] * pa_x[i] - ta1_z_yy_xxzzz_1[i] * pc_x[i];

        ta1_z_xyy_xyyyy_0[i] = ta1_z_yy_yyyy_0[i] * fe_0 - ta1_z_yy_yyyy_1[i] * fe_0 + ta1_z_yy_xyyyy_0[i] * pa_x[i] - ta1_z_yy_xyyyy_1[i] * pc_x[i];

        ta1_z_xyy_xyyyz_0[i] = ta1_z_yy_yyyz_0[i] * fe_0 - ta1_z_yy_yyyz_1[i] * fe_0 + ta1_z_yy_xyyyz_0[i] * pa_x[i] - ta1_z_yy_xyyyz_1[i] * pc_x[i];

        ta1_z_xyy_xyyzz_0[i] = ta1_z_yy_yyzz_0[i] * fe_0 - ta1_z_yy_yyzz_1[i] * fe_0 + ta1_z_yy_xyyzz_0[i] * pa_x[i] - ta1_z_yy_xyyzz_1[i] * pc_x[i];

        ta1_z_xyy_xyzzz_0[i] = ta1_z_yy_yzzz_0[i] * fe_0 - ta1_z_yy_yzzz_1[i] * fe_0 + ta1_z_yy_xyzzz_0[i] * pa_x[i] - ta1_z_yy_xyzzz_1[i] * pc_x[i];

        ta1_z_xyy_xzzzz_0[i] = ta1_z_yy_zzzz_0[i] * fe_0 - ta1_z_yy_zzzz_1[i] * fe_0 + ta1_z_yy_xzzzz_0[i] * pa_x[i] - ta1_z_yy_xzzzz_1[i] * pc_x[i];

        ta1_z_xyy_yyyyy_0[i] = ta1_z_yy_yyyyy_0[i] * pa_x[i] - ta1_z_yy_yyyyy_1[i] * pc_x[i];

        ta1_z_xyy_yyyyz_0[i] = ta1_z_yy_yyyyz_0[i] * pa_x[i] - ta1_z_yy_yyyyz_1[i] * pc_x[i];

        ta1_z_xyy_yyyzz_0[i] = ta1_z_yy_yyyzz_0[i] * pa_x[i] - ta1_z_yy_yyyzz_1[i] * pc_x[i];

        ta1_z_xyy_yyzzz_0[i] = ta1_z_yy_yyzzz_0[i] * pa_x[i] - ta1_z_yy_yyzzz_1[i] * pc_x[i];

        ta1_z_xyy_yzzzz_0[i] = ta1_z_yy_yzzzz_0[i] * pa_x[i] - ta1_z_yy_yzzzz_1[i] * pc_x[i];

        ta1_z_xyy_zzzzz_0[i] = ta1_z_yy_zzzzz_0[i] * pa_x[i] - ta1_z_yy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 504-525 components of targeted buffer : FH

    auto ta1_z_xyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 504);

    auto ta1_z_xyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 505);

    auto ta1_z_xyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 506);

    auto ta1_z_xyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 507);

    auto ta1_z_xyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 508);

    auto ta1_z_xyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 509);

    auto ta1_z_xyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 510);

    auto ta1_z_xyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 511);

    auto ta1_z_xyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 512);

    auto ta1_z_xyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 513);

    auto ta1_z_xyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 514);

    auto ta1_z_xyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 515);

    auto ta1_z_xyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 516);

    auto ta1_z_xyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 517);

    auto ta1_z_xyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 518);

    auto ta1_z_xyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 519);

    auto ta1_z_xyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 520);

    auto ta1_z_xyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 521);

    auto ta1_z_xyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 522);

    auto ta1_z_xyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 523);

    auto ta1_z_xyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 524);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_xy_xxxxy_0,  \
                             ta1_z_xy_xxxxy_1,  \
                             ta1_z_xy_xxxyy_0,  \
                             ta1_z_xy_xxxyy_1,  \
                             ta1_z_xy_xxyyy_0,  \
                             ta1_z_xy_xxyyy_1,  \
                             ta1_z_xy_xyyyy_0,  \
                             ta1_z_xy_xyyyy_1,  \
                             ta1_z_xyz_xxxxx_0, \
                             ta1_z_xyz_xxxxy_0, \
                             ta1_z_xyz_xxxxz_0, \
                             ta1_z_xyz_xxxyy_0, \
                             ta1_z_xyz_xxxyz_0, \
                             ta1_z_xyz_xxxzz_0, \
                             ta1_z_xyz_xxyyy_0, \
                             ta1_z_xyz_xxyyz_0, \
                             ta1_z_xyz_xxyzz_0, \
                             ta1_z_xyz_xxzzz_0, \
                             ta1_z_xyz_xyyyy_0, \
                             ta1_z_xyz_xyyyz_0, \
                             ta1_z_xyz_xyyzz_0, \
                             ta1_z_xyz_xyzzz_0, \
                             ta1_z_xyz_xzzzz_0, \
                             ta1_z_xyz_yyyyy_0, \
                             ta1_z_xyz_yyyyz_0, \
                             ta1_z_xyz_yyyzz_0, \
                             ta1_z_xyz_yyzzz_0, \
                             ta1_z_xyz_yzzzz_0, \
                             ta1_z_xyz_zzzzz_0, \
                             ta1_z_xz_xxxxx_0,  \
                             ta1_z_xz_xxxxx_1,  \
                             ta1_z_xz_xxxxz_0,  \
                             ta1_z_xz_xxxxz_1,  \
                             ta1_z_xz_xxxzz_0,  \
                             ta1_z_xz_xxxzz_1,  \
                             ta1_z_xz_xxzzz_0,  \
                             ta1_z_xz_xxzzz_1,  \
                             ta1_z_xz_xzzzz_0,  \
                             ta1_z_xz_xzzzz_1,  \
                             ta1_z_yz_xxxyz_0,  \
                             ta1_z_yz_xxxyz_1,  \
                             ta1_z_yz_xxyyz_0,  \
                             ta1_z_yz_xxyyz_1,  \
                             ta1_z_yz_xxyz_0,   \
                             ta1_z_yz_xxyz_1,   \
                             ta1_z_yz_xxyzz_0,  \
                             ta1_z_yz_xxyzz_1,  \
                             ta1_z_yz_xyyyz_0,  \
                             ta1_z_yz_xyyyz_1,  \
                             ta1_z_yz_xyyz_0,   \
                             ta1_z_yz_xyyz_1,   \
                             ta1_z_yz_xyyzz_0,  \
                             ta1_z_yz_xyyzz_1,  \
                             ta1_z_yz_xyzz_0,   \
                             ta1_z_yz_xyzz_1,   \
                             ta1_z_yz_xyzzz_0,  \
                             ta1_z_yz_xyzzz_1,  \
                             ta1_z_yz_yyyyy_0,  \
                             ta1_z_yz_yyyyy_1,  \
                             ta1_z_yz_yyyyz_0,  \
                             ta1_z_yz_yyyyz_1,  \
                             ta1_z_yz_yyyz_0,   \
                             ta1_z_yz_yyyz_1,   \
                             ta1_z_yz_yyyzz_0,  \
                             ta1_z_yz_yyyzz_1,  \
                             ta1_z_yz_yyzz_0,   \
                             ta1_z_yz_yyzz_1,   \
                             ta1_z_yz_yyzzz_0,  \
                             ta1_z_yz_yyzzz_1,  \
                             ta1_z_yz_yzzz_0,   \
                             ta1_z_yz_yzzz_1,   \
                             ta1_z_yz_yzzzz_0,  \
                             ta1_z_yz_yzzzz_1,  \
                             ta1_z_yz_zzzzz_0,  \
                             ta1_z_yz_zzzzz_1,  \
                             ta_xy_xxxxy_1,     \
                             ta_xy_xxxyy_1,     \
                             ta_xy_xxyyy_1,     \
                             ta_xy_xyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyz_xxxxx_0[i] = ta1_z_xz_xxxxx_0[i] * pa_y[i] - ta1_z_xz_xxxxx_1[i] * pc_y[i];

        ta1_z_xyz_xxxxy_0[i] = ta_xy_xxxxy_1[i] + ta1_z_xy_xxxxy_0[i] * pa_z[i] - ta1_z_xy_xxxxy_1[i] * pc_z[i];

        ta1_z_xyz_xxxxz_0[i] = ta1_z_xz_xxxxz_0[i] * pa_y[i] - ta1_z_xz_xxxxz_1[i] * pc_y[i];

        ta1_z_xyz_xxxyy_0[i] = ta_xy_xxxyy_1[i] + ta1_z_xy_xxxyy_0[i] * pa_z[i] - ta1_z_xy_xxxyy_1[i] * pc_z[i];

        ta1_z_xyz_xxxyz_0[i] =
            3.0 * ta1_z_yz_xxyz_0[i] * fe_0 - 3.0 * ta1_z_yz_xxyz_1[i] * fe_0 + ta1_z_yz_xxxyz_0[i] * pa_x[i] - ta1_z_yz_xxxyz_1[i] * pc_x[i];

        ta1_z_xyz_xxxzz_0[i] = ta1_z_xz_xxxzz_0[i] * pa_y[i] - ta1_z_xz_xxxzz_1[i] * pc_y[i];

        ta1_z_xyz_xxyyy_0[i] = ta_xy_xxyyy_1[i] + ta1_z_xy_xxyyy_0[i] * pa_z[i] - ta1_z_xy_xxyyy_1[i] * pc_z[i];

        ta1_z_xyz_xxyyz_0[i] =
            2.0 * ta1_z_yz_xyyz_0[i] * fe_0 - 2.0 * ta1_z_yz_xyyz_1[i] * fe_0 + ta1_z_yz_xxyyz_0[i] * pa_x[i] - ta1_z_yz_xxyyz_1[i] * pc_x[i];

        ta1_z_xyz_xxyzz_0[i] =
            2.0 * ta1_z_yz_xyzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xyzz_1[i] * fe_0 + ta1_z_yz_xxyzz_0[i] * pa_x[i] - ta1_z_yz_xxyzz_1[i] * pc_x[i];

        ta1_z_xyz_xxzzz_0[i] = ta1_z_xz_xxzzz_0[i] * pa_y[i] - ta1_z_xz_xxzzz_1[i] * pc_y[i];

        ta1_z_xyz_xyyyy_0[i] = ta_xy_xyyyy_1[i] + ta1_z_xy_xyyyy_0[i] * pa_z[i] - ta1_z_xy_xyyyy_1[i] * pc_z[i];

        ta1_z_xyz_xyyyz_0[i] = ta1_z_yz_yyyz_0[i] * fe_0 - ta1_z_yz_yyyz_1[i] * fe_0 + ta1_z_yz_xyyyz_0[i] * pa_x[i] - ta1_z_yz_xyyyz_1[i] * pc_x[i];

        ta1_z_xyz_xyyzz_0[i] = ta1_z_yz_yyzz_0[i] * fe_0 - ta1_z_yz_yyzz_1[i] * fe_0 + ta1_z_yz_xyyzz_0[i] * pa_x[i] - ta1_z_yz_xyyzz_1[i] * pc_x[i];

        ta1_z_xyz_xyzzz_0[i] = ta1_z_yz_yzzz_0[i] * fe_0 - ta1_z_yz_yzzz_1[i] * fe_0 + ta1_z_yz_xyzzz_0[i] * pa_x[i] - ta1_z_yz_xyzzz_1[i] * pc_x[i];

        ta1_z_xyz_xzzzz_0[i] = ta1_z_xz_xzzzz_0[i] * pa_y[i] - ta1_z_xz_xzzzz_1[i] * pc_y[i];

        ta1_z_xyz_yyyyy_0[i] = ta1_z_yz_yyyyy_0[i] * pa_x[i] - ta1_z_yz_yyyyy_1[i] * pc_x[i];

        ta1_z_xyz_yyyyz_0[i] = ta1_z_yz_yyyyz_0[i] * pa_x[i] - ta1_z_yz_yyyyz_1[i] * pc_x[i];

        ta1_z_xyz_yyyzz_0[i] = ta1_z_yz_yyyzz_0[i] * pa_x[i] - ta1_z_yz_yyyzz_1[i] * pc_x[i];

        ta1_z_xyz_yyzzz_0[i] = ta1_z_yz_yyzzz_0[i] * pa_x[i] - ta1_z_yz_yyzzz_1[i] * pc_x[i];

        ta1_z_xyz_yzzzz_0[i] = ta1_z_yz_yzzzz_0[i] * pa_x[i] - ta1_z_yz_yzzzz_1[i] * pc_x[i];

        ta1_z_xyz_zzzzz_0[i] = ta1_z_yz_zzzzz_0[i] * pa_x[i] - ta1_z_yz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 525-546 components of targeted buffer : FH

    auto ta1_z_xzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 525);

    auto ta1_z_xzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 526);

    auto ta1_z_xzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 527);

    auto ta1_z_xzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 528);

    auto ta1_z_xzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 529);

    auto ta1_z_xzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 530);

    auto ta1_z_xzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 531);

    auto ta1_z_xzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 532);

    auto ta1_z_xzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 533);

    auto ta1_z_xzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 534);

    auto ta1_z_xzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 535);

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

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_z_xzz_xxxxx_0, \
                             ta1_z_xzz_xxxxy_0, \
                             ta1_z_xzz_xxxxz_0, \
                             ta1_z_xzz_xxxyy_0, \
                             ta1_z_xzz_xxxyz_0, \
                             ta1_z_xzz_xxxzz_0, \
                             ta1_z_xzz_xxyyy_0, \
                             ta1_z_xzz_xxyyz_0, \
                             ta1_z_xzz_xxyzz_0, \
                             ta1_z_xzz_xxzzz_0, \
                             ta1_z_xzz_xyyyy_0, \
                             ta1_z_xzz_xyyyz_0, \
                             ta1_z_xzz_xyyzz_0, \
                             ta1_z_xzz_xyzzz_0, \
                             ta1_z_xzz_xzzzz_0, \
                             ta1_z_xzz_yyyyy_0, \
                             ta1_z_xzz_yyyyz_0, \
                             ta1_z_xzz_yyyzz_0, \
                             ta1_z_xzz_yyzzz_0, \
                             ta1_z_xzz_yzzzz_0, \
                             ta1_z_xzz_zzzzz_0, \
                             ta1_z_zz_xxxx_0,   \
                             ta1_z_zz_xxxx_1,   \
                             ta1_z_zz_xxxxx_0,  \
                             ta1_z_zz_xxxxx_1,  \
                             ta1_z_zz_xxxxy_0,  \
                             ta1_z_zz_xxxxy_1,  \
                             ta1_z_zz_xxxxz_0,  \
                             ta1_z_zz_xxxxz_1,  \
                             ta1_z_zz_xxxy_0,   \
                             ta1_z_zz_xxxy_1,   \
                             ta1_z_zz_xxxyy_0,  \
                             ta1_z_zz_xxxyy_1,  \
                             ta1_z_zz_xxxyz_0,  \
                             ta1_z_zz_xxxyz_1,  \
                             ta1_z_zz_xxxz_0,   \
                             ta1_z_zz_xxxz_1,   \
                             ta1_z_zz_xxxzz_0,  \
                             ta1_z_zz_xxxzz_1,  \
                             ta1_z_zz_xxyy_0,   \
                             ta1_z_zz_xxyy_1,   \
                             ta1_z_zz_xxyyy_0,  \
                             ta1_z_zz_xxyyy_1,  \
                             ta1_z_zz_xxyyz_0,  \
                             ta1_z_zz_xxyyz_1,  \
                             ta1_z_zz_xxyz_0,   \
                             ta1_z_zz_xxyz_1,   \
                             ta1_z_zz_xxyzz_0,  \
                             ta1_z_zz_xxyzz_1,  \
                             ta1_z_zz_xxzz_0,   \
                             ta1_z_zz_xxzz_1,   \
                             ta1_z_zz_xxzzz_0,  \
                             ta1_z_zz_xxzzz_1,  \
                             ta1_z_zz_xyyy_0,   \
                             ta1_z_zz_xyyy_1,   \
                             ta1_z_zz_xyyyy_0,  \
                             ta1_z_zz_xyyyy_1,  \
                             ta1_z_zz_xyyyz_0,  \
                             ta1_z_zz_xyyyz_1,  \
                             ta1_z_zz_xyyz_0,   \
                             ta1_z_zz_xyyz_1,   \
                             ta1_z_zz_xyyzz_0,  \
                             ta1_z_zz_xyyzz_1,  \
                             ta1_z_zz_xyzz_0,   \
                             ta1_z_zz_xyzz_1,   \
                             ta1_z_zz_xyzzz_0,  \
                             ta1_z_zz_xyzzz_1,  \
                             ta1_z_zz_xzzz_0,   \
                             ta1_z_zz_xzzz_1,   \
                             ta1_z_zz_xzzzz_0,  \
                             ta1_z_zz_xzzzz_1,  \
                             ta1_z_zz_yyyy_0,   \
                             ta1_z_zz_yyyy_1,   \
                             ta1_z_zz_yyyyy_0,  \
                             ta1_z_zz_yyyyy_1,  \
                             ta1_z_zz_yyyyz_0,  \
                             ta1_z_zz_yyyyz_1,  \
                             ta1_z_zz_yyyz_0,   \
                             ta1_z_zz_yyyz_1,   \
                             ta1_z_zz_yyyzz_0,  \
                             ta1_z_zz_yyyzz_1,  \
                             ta1_z_zz_yyzz_0,   \
                             ta1_z_zz_yyzz_1,   \
                             ta1_z_zz_yyzzz_0,  \
                             ta1_z_zz_yyzzz_1,  \
                             ta1_z_zz_yzzz_0,   \
                             ta1_z_zz_yzzz_1,   \
                             ta1_z_zz_yzzzz_0,  \
                             ta1_z_zz_yzzzz_1,  \
                             ta1_z_zz_zzzz_0,   \
                             ta1_z_zz_zzzz_1,   \
                             ta1_z_zz_zzzzz_0,  \
                             ta1_z_zz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzz_xxxxx_0[i] =
            5.0 * ta1_z_zz_xxxx_0[i] * fe_0 - 5.0 * ta1_z_zz_xxxx_1[i] * fe_0 + ta1_z_zz_xxxxx_0[i] * pa_x[i] - ta1_z_zz_xxxxx_1[i] * pc_x[i];

        ta1_z_xzz_xxxxy_0[i] =
            4.0 * ta1_z_zz_xxxy_0[i] * fe_0 - 4.0 * ta1_z_zz_xxxy_1[i] * fe_0 + ta1_z_zz_xxxxy_0[i] * pa_x[i] - ta1_z_zz_xxxxy_1[i] * pc_x[i];

        ta1_z_xzz_xxxxz_0[i] =
            4.0 * ta1_z_zz_xxxz_0[i] * fe_0 - 4.0 * ta1_z_zz_xxxz_1[i] * fe_0 + ta1_z_zz_xxxxz_0[i] * pa_x[i] - ta1_z_zz_xxxxz_1[i] * pc_x[i];

        ta1_z_xzz_xxxyy_0[i] =
            3.0 * ta1_z_zz_xxyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyy_1[i] * fe_0 + ta1_z_zz_xxxyy_0[i] * pa_x[i] - ta1_z_zz_xxxyy_1[i] * pc_x[i];

        ta1_z_xzz_xxxyz_0[i] =
            3.0 * ta1_z_zz_xxyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyz_1[i] * fe_0 + ta1_z_zz_xxxyz_0[i] * pa_x[i] - ta1_z_zz_xxxyz_1[i] * pc_x[i];

        ta1_z_xzz_xxxzz_0[i] =
            3.0 * ta1_z_zz_xxzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxzz_1[i] * fe_0 + ta1_z_zz_xxxzz_0[i] * pa_x[i] - ta1_z_zz_xxxzz_1[i] * pc_x[i];

        ta1_z_xzz_xxyyy_0[i] =
            2.0 * ta1_z_zz_xyyy_0[i] * fe_0 - 2.0 * ta1_z_zz_xyyy_1[i] * fe_0 + ta1_z_zz_xxyyy_0[i] * pa_x[i] - ta1_z_zz_xxyyy_1[i] * pc_x[i];

        ta1_z_xzz_xxyyz_0[i] =
            2.0 * ta1_z_zz_xyyz_0[i] * fe_0 - 2.0 * ta1_z_zz_xyyz_1[i] * fe_0 + ta1_z_zz_xxyyz_0[i] * pa_x[i] - ta1_z_zz_xxyyz_1[i] * pc_x[i];

        ta1_z_xzz_xxyzz_0[i] =
            2.0 * ta1_z_zz_xyzz_0[i] * fe_0 - 2.0 * ta1_z_zz_xyzz_1[i] * fe_0 + ta1_z_zz_xxyzz_0[i] * pa_x[i] - ta1_z_zz_xxyzz_1[i] * pc_x[i];

        ta1_z_xzz_xxzzz_0[i] =
            2.0 * ta1_z_zz_xzzz_0[i] * fe_0 - 2.0 * ta1_z_zz_xzzz_1[i] * fe_0 + ta1_z_zz_xxzzz_0[i] * pa_x[i] - ta1_z_zz_xxzzz_1[i] * pc_x[i];

        ta1_z_xzz_xyyyy_0[i] = ta1_z_zz_yyyy_0[i] * fe_0 - ta1_z_zz_yyyy_1[i] * fe_0 + ta1_z_zz_xyyyy_0[i] * pa_x[i] - ta1_z_zz_xyyyy_1[i] * pc_x[i];

        ta1_z_xzz_xyyyz_0[i] = ta1_z_zz_yyyz_0[i] * fe_0 - ta1_z_zz_yyyz_1[i] * fe_0 + ta1_z_zz_xyyyz_0[i] * pa_x[i] - ta1_z_zz_xyyyz_1[i] * pc_x[i];

        ta1_z_xzz_xyyzz_0[i] = ta1_z_zz_yyzz_0[i] * fe_0 - ta1_z_zz_yyzz_1[i] * fe_0 + ta1_z_zz_xyyzz_0[i] * pa_x[i] - ta1_z_zz_xyyzz_1[i] * pc_x[i];

        ta1_z_xzz_xyzzz_0[i] = ta1_z_zz_yzzz_0[i] * fe_0 - ta1_z_zz_yzzz_1[i] * fe_0 + ta1_z_zz_xyzzz_0[i] * pa_x[i] - ta1_z_zz_xyzzz_1[i] * pc_x[i];

        ta1_z_xzz_xzzzz_0[i] = ta1_z_zz_zzzz_0[i] * fe_0 - ta1_z_zz_zzzz_1[i] * fe_0 + ta1_z_zz_xzzzz_0[i] * pa_x[i] - ta1_z_zz_xzzzz_1[i] * pc_x[i];

        ta1_z_xzz_yyyyy_0[i] = ta1_z_zz_yyyyy_0[i] * pa_x[i] - ta1_z_zz_yyyyy_1[i] * pc_x[i];

        ta1_z_xzz_yyyyz_0[i] = ta1_z_zz_yyyyz_0[i] * pa_x[i] - ta1_z_zz_yyyyz_1[i] * pc_x[i];

        ta1_z_xzz_yyyzz_0[i] = ta1_z_zz_yyyzz_0[i] * pa_x[i] - ta1_z_zz_yyyzz_1[i] * pc_x[i];

        ta1_z_xzz_yyzzz_0[i] = ta1_z_zz_yyzzz_0[i] * pa_x[i] - ta1_z_zz_yyzzz_1[i] * pc_x[i];

        ta1_z_xzz_yzzzz_0[i] = ta1_z_zz_yzzzz_0[i] * pa_x[i] - ta1_z_zz_yzzzz_1[i] * pc_x[i];

        ta1_z_xzz_zzzzz_0[i] = ta1_z_zz_zzzzz_0[i] * pa_x[i] - ta1_z_zz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 546-567 components of targeted buffer : FH

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

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_z_y_xxxxx_0,   \
                             ta1_z_y_xxxxx_1,   \
                             ta1_z_y_xxxxy_0,   \
                             ta1_z_y_xxxxy_1,   \
                             ta1_z_y_xxxxz_0,   \
                             ta1_z_y_xxxxz_1,   \
                             ta1_z_y_xxxyy_0,   \
                             ta1_z_y_xxxyy_1,   \
                             ta1_z_y_xxxyz_0,   \
                             ta1_z_y_xxxyz_1,   \
                             ta1_z_y_xxxzz_0,   \
                             ta1_z_y_xxxzz_1,   \
                             ta1_z_y_xxyyy_0,   \
                             ta1_z_y_xxyyy_1,   \
                             ta1_z_y_xxyyz_0,   \
                             ta1_z_y_xxyyz_1,   \
                             ta1_z_y_xxyzz_0,   \
                             ta1_z_y_xxyzz_1,   \
                             ta1_z_y_xxzzz_0,   \
                             ta1_z_y_xxzzz_1,   \
                             ta1_z_y_xyyyy_0,   \
                             ta1_z_y_xyyyy_1,   \
                             ta1_z_y_xyyyz_0,   \
                             ta1_z_y_xyyyz_1,   \
                             ta1_z_y_xyyzz_0,   \
                             ta1_z_y_xyyzz_1,   \
                             ta1_z_y_xyzzz_0,   \
                             ta1_z_y_xyzzz_1,   \
                             ta1_z_y_xzzzz_0,   \
                             ta1_z_y_xzzzz_1,   \
                             ta1_z_y_yyyyy_0,   \
                             ta1_z_y_yyyyy_1,   \
                             ta1_z_y_yyyyz_0,   \
                             ta1_z_y_yyyyz_1,   \
                             ta1_z_y_yyyzz_0,   \
                             ta1_z_y_yyyzz_1,   \
                             ta1_z_y_yyzzz_0,   \
                             ta1_z_y_yyzzz_1,   \
                             ta1_z_y_yzzzz_0,   \
                             ta1_z_y_yzzzz_1,   \
                             ta1_z_y_zzzzz_0,   \
                             ta1_z_y_zzzzz_1,   \
                             ta1_z_yy_xxxx_0,   \
                             ta1_z_yy_xxxx_1,   \
                             ta1_z_yy_xxxxx_0,  \
                             ta1_z_yy_xxxxx_1,  \
                             ta1_z_yy_xxxxy_0,  \
                             ta1_z_yy_xxxxy_1,  \
                             ta1_z_yy_xxxxz_0,  \
                             ta1_z_yy_xxxxz_1,  \
                             ta1_z_yy_xxxy_0,   \
                             ta1_z_yy_xxxy_1,   \
                             ta1_z_yy_xxxyy_0,  \
                             ta1_z_yy_xxxyy_1,  \
                             ta1_z_yy_xxxyz_0,  \
                             ta1_z_yy_xxxyz_1,  \
                             ta1_z_yy_xxxz_0,   \
                             ta1_z_yy_xxxz_1,   \
                             ta1_z_yy_xxxzz_0,  \
                             ta1_z_yy_xxxzz_1,  \
                             ta1_z_yy_xxyy_0,   \
                             ta1_z_yy_xxyy_1,   \
                             ta1_z_yy_xxyyy_0,  \
                             ta1_z_yy_xxyyy_1,  \
                             ta1_z_yy_xxyyz_0,  \
                             ta1_z_yy_xxyyz_1,  \
                             ta1_z_yy_xxyz_0,   \
                             ta1_z_yy_xxyz_1,   \
                             ta1_z_yy_xxyzz_0,  \
                             ta1_z_yy_xxyzz_1,  \
                             ta1_z_yy_xxzz_0,   \
                             ta1_z_yy_xxzz_1,   \
                             ta1_z_yy_xxzzz_0,  \
                             ta1_z_yy_xxzzz_1,  \
                             ta1_z_yy_xyyy_0,   \
                             ta1_z_yy_xyyy_1,   \
                             ta1_z_yy_xyyyy_0,  \
                             ta1_z_yy_xyyyy_1,  \
                             ta1_z_yy_xyyyz_0,  \
                             ta1_z_yy_xyyyz_1,  \
                             ta1_z_yy_xyyz_0,   \
                             ta1_z_yy_xyyz_1,   \
                             ta1_z_yy_xyyzz_0,  \
                             ta1_z_yy_xyyzz_1,  \
                             ta1_z_yy_xyzz_0,   \
                             ta1_z_yy_xyzz_1,   \
                             ta1_z_yy_xyzzz_0,  \
                             ta1_z_yy_xyzzz_1,  \
                             ta1_z_yy_xzzz_0,   \
                             ta1_z_yy_xzzz_1,   \
                             ta1_z_yy_xzzzz_0,  \
                             ta1_z_yy_xzzzz_1,  \
                             ta1_z_yy_yyyy_0,   \
                             ta1_z_yy_yyyy_1,   \
                             ta1_z_yy_yyyyy_0,  \
                             ta1_z_yy_yyyyy_1,  \
                             ta1_z_yy_yyyyz_0,  \
                             ta1_z_yy_yyyyz_1,  \
                             ta1_z_yy_yyyz_0,   \
                             ta1_z_yy_yyyz_1,   \
                             ta1_z_yy_yyyzz_0,  \
                             ta1_z_yy_yyyzz_1,  \
                             ta1_z_yy_yyzz_0,   \
                             ta1_z_yy_yyzz_1,   \
                             ta1_z_yy_yyzzz_0,  \
                             ta1_z_yy_yyzzz_1,  \
                             ta1_z_yy_yzzz_0,   \
                             ta1_z_yy_yzzz_1,   \
                             ta1_z_yy_yzzzz_0,  \
                             ta1_z_yy_yzzzz_1,  \
                             ta1_z_yy_zzzz_0,   \
                             ta1_z_yy_zzzz_1,   \
                             ta1_z_yy_zzzzz_0,  \
                             ta1_z_yy_zzzzz_1,  \
                             ta1_z_yyy_xxxxx_0, \
                             ta1_z_yyy_xxxxy_0, \
                             ta1_z_yyy_xxxxz_0, \
                             ta1_z_yyy_xxxyy_0, \
                             ta1_z_yyy_xxxyz_0, \
                             ta1_z_yyy_xxxzz_0, \
                             ta1_z_yyy_xxyyy_0, \
                             ta1_z_yyy_xxyyz_0, \
                             ta1_z_yyy_xxyzz_0, \
                             ta1_z_yyy_xxzzz_0, \
                             ta1_z_yyy_xyyyy_0, \
                             ta1_z_yyy_xyyyz_0, \
                             ta1_z_yyy_xyyzz_0, \
                             ta1_z_yyy_xyzzz_0, \
                             ta1_z_yyy_xzzzz_0, \
                             ta1_z_yyy_yyyyy_0, \
                             ta1_z_yyy_yyyyz_0, \
                             ta1_z_yyy_yyyzz_0, \
                             ta1_z_yyy_yyzzz_0, \
                             ta1_z_yyy_yzzzz_0, \
                             ta1_z_yyy_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyy_xxxxx_0[i] =
            2.0 * ta1_z_y_xxxxx_0[i] * fe_0 - 2.0 * ta1_z_y_xxxxx_1[i] * fe_0 + ta1_z_yy_xxxxx_0[i] * pa_y[i] - ta1_z_yy_xxxxx_1[i] * pc_y[i];

        ta1_z_yyy_xxxxy_0[i] = 2.0 * ta1_z_y_xxxxy_0[i] * fe_0 - 2.0 * ta1_z_y_xxxxy_1[i] * fe_0 + ta1_z_yy_xxxx_0[i] * fe_0 -
                               ta1_z_yy_xxxx_1[i] * fe_0 + ta1_z_yy_xxxxy_0[i] * pa_y[i] - ta1_z_yy_xxxxy_1[i] * pc_y[i];

        ta1_z_yyy_xxxxz_0[i] =
            2.0 * ta1_z_y_xxxxz_0[i] * fe_0 - 2.0 * ta1_z_y_xxxxz_1[i] * fe_0 + ta1_z_yy_xxxxz_0[i] * pa_y[i] - ta1_z_yy_xxxxz_1[i] * pc_y[i];

        ta1_z_yyy_xxxyy_0[i] = 2.0 * ta1_z_y_xxxyy_0[i] * fe_0 - 2.0 * ta1_z_y_xxxyy_1[i] * fe_0 + 2.0 * ta1_z_yy_xxxy_0[i] * fe_0 -
                               2.0 * ta1_z_yy_xxxy_1[i] * fe_0 + ta1_z_yy_xxxyy_0[i] * pa_y[i] - ta1_z_yy_xxxyy_1[i] * pc_y[i];

        ta1_z_yyy_xxxyz_0[i] = 2.0 * ta1_z_y_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_y_xxxyz_1[i] * fe_0 + ta1_z_yy_xxxz_0[i] * fe_0 -
                               ta1_z_yy_xxxz_1[i] * fe_0 + ta1_z_yy_xxxyz_0[i] * pa_y[i] - ta1_z_yy_xxxyz_1[i] * pc_y[i];

        ta1_z_yyy_xxxzz_0[i] =
            2.0 * ta1_z_y_xxxzz_0[i] * fe_0 - 2.0 * ta1_z_y_xxxzz_1[i] * fe_0 + ta1_z_yy_xxxzz_0[i] * pa_y[i] - ta1_z_yy_xxxzz_1[i] * pc_y[i];

        ta1_z_yyy_xxyyy_0[i] = 2.0 * ta1_z_y_xxyyy_0[i] * fe_0 - 2.0 * ta1_z_y_xxyyy_1[i] * fe_0 + 3.0 * ta1_z_yy_xxyy_0[i] * fe_0 -
                               3.0 * ta1_z_yy_xxyy_1[i] * fe_0 + ta1_z_yy_xxyyy_0[i] * pa_y[i] - ta1_z_yy_xxyyy_1[i] * pc_y[i];

        ta1_z_yyy_xxyyz_0[i] = 2.0 * ta1_z_y_xxyyz_0[i] * fe_0 - 2.0 * ta1_z_y_xxyyz_1[i] * fe_0 + 2.0 * ta1_z_yy_xxyz_0[i] * fe_0 -
                               2.0 * ta1_z_yy_xxyz_1[i] * fe_0 + ta1_z_yy_xxyyz_0[i] * pa_y[i] - ta1_z_yy_xxyyz_1[i] * pc_y[i];

        ta1_z_yyy_xxyzz_0[i] = 2.0 * ta1_z_y_xxyzz_0[i] * fe_0 - 2.0 * ta1_z_y_xxyzz_1[i] * fe_0 + ta1_z_yy_xxzz_0[i] * fe_0 -
                               ta1_z_yy_xxzz_1[i] * fe_0 + ta1_z_yy_xxyzz_0[i] * pa_y[i] - ta1_z_yy_xxyzz_1[i] * pc_y[i];

        ta1_z_yyy_xxzzz_0[i] =
            2.0 * ta1_z_y_xxzzz_0[i] * fe_0 - 2.0 * ta1_z_y_xxzzz_1[i] * fe_0 + ta1_z_yy_xxzzz_0[i] * pa_y[i] - ta1_z_yy_xxzzz_1[i] * pc_y[i];

        ta1_z_yyy_xyyyy_0[i] = 2.0 * ta1_z_y_xyyyy_0[i] * fe_0 - 2.0 * ta1_z_y_xyyyy_1[i] * fe_0 + 4.0 * ta1_z_yy_xyyy_0[i] * fe_0 -
                               4.0 * ta1_z_yy_xyyy_1[i] * fe_0 + ta1_z_yy_xyyyy_0[i] * pa_y[i] - ta1_z_yy_xyyyy_1[i] * pc_y[i];

        ta1_z_yyy_xyyyz_0[i] = 2.0 * ta1_z_y_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_y_xyyyz_1[i] * fe_0 + 3.0 * ta1_z_yy_xyyz_0[i] * fe_0 -
                               3.0 * ta1_z_yy_xyyz_1[i] * fe_0 + ta1_z_yy_xyyyz_0[i] * pa_y[i] - ta1_z_yy_xyyyz_1[i] * pc_y[i];

        ta1_z_yyy_xyyzz_0[i] = 2.0 * ta1_z_y_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_y_xyyzz_1[i] * fe_0 + 2.0 * ta1_z_yy_xyzz_0[i] * fe_0 -
                               2.0 * ta1_z_yy_xyzz_1[i] * fe_0 + ta1_z_yy_xyyzz_0[i] * pa_y[i] - ta1_z_yy_xyyzz_1[i] * pc_y[i];

        ta1_z_yyy_xyzzz_0[i] = 2.0 * ta1_z_y_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_y_xyzzz_1[i] * fe_0 + ta1_z_yy_xzzz_0[i] * fe_0 -
                               ta1_z_yy_xzzz_1[i] * fe_0 + ta1_z_yy_xyzzz_0[i] * pa_y[i] - ta1_z_yy_xyzzz_1[i] * pc_y[i];

        ta1_z_yyy_xzzzz_0[i] =
            2.0 * ta1_z_y_xzzzz_0[i] * fe_0 - 2.0 * ta1_z_y_xzzzz_1[i] * fe_0 + ta1_z_yy_xzzzz_0[i] * pa_y[i] - ta1_z_yy_xzzzz_1[i] * pc_y[i];

        ta1_z_yyy_yyyyy_0[i] = 2.0 * ta1_z_y_yyyyy_0[i] * fe_0 - 2.0 * ta1_z_y_yyyyy_1[i] * fe_0 + 5.0 * ta1_z_yy_yyyy_0[i] * fe_0 -
                               5.0 * ta1_z_yy_yyyy_1[i] * fe_0 + ta1_z_yy_yyyyy_0[i] * pa_y[i] - ta1_z_yy_yyyyy_1[i] * pc_y[i];

        ta1_z_yyy_yyyyz_0[i] = 2.0 * ta1_z_y_yyyyz_0[i] * fe_0 - 2.0 * ta1_z_y_yyyyz_1[i] * fe_0 + 4.0 * ta1_z_yy_yyyz_0[i] * fe_0 -
                               4.0 * ta1_z_yy_yyyz_1[i] * fe_0 + ta1_z_yy_yyyyz_0[i] * pa_y[i] - ta1_z_yy_yyyyz_1[i] * pc_y[i];

        ta1_z_yyy_yyyzz_0[i] = 2.0 * ta1_z_y_yyyzz_0[i] * fe_0 - 2.0 * ta1_z_y_yyyzz_1[i] * fe_0 + 3.0 * ta1_z_yy_yyzz_0[i] * fe_0 -
                               3.0 * ta1_z_yy_yyzz_1[i] * fe_0 + ta1_z_yy_yyyzz_0[i] * pa_y[i] - ta1_z_yy_yyyzz_1[i] * pc_y[i];

        ta1_z_yyy_yyzzz_0[i] = 2.0 * ta1_z_y_yyzzz_0[i] * fe_0 - 2.0 * ta1_z_y_yyzzz_1[i] * fe_0 + 2.0 * ta1_z_yy_yzzz_0[i] * fe_0 -
                               2.0 * ta1_z_yy_yzzz_1[i] * fe_0 + ta1_z_yy_yyzzz_0[i] * pa_y[i] - ta1_z_yy_yyzzz_1[i] * pc_y[i];

        ta1_z_yyy_yzzzz_0[i] = 2.0 * ta1_z_y_yzzzz_0[i] * fe_0 - 2.0 * ta1_z_y_yzzzz_1[i] * fe_0 + ta1_z_yy_zzzz_0[i] * fe_0 -
                               ta1_z_yy_zzzz_1[i] * fe_0 + ta1_z_yy_yzzzz_0[i] * pa_y[i] - ta1_z_yy_yzzzz_1[i] * pc_y[i];

        ta1_z_yyy_zzzzz_0[i] =
            2.0 * ta1_z_y_zzzzz_0[i] * fe_0 - 2.0 * ta1_z_y_zzzzz_1[i] * fe_0 + ta1_z_yy_zzzzz_0[i] * pa_y[i] - ta1_z_yy_zzzzz_1[i] * pc_y[i];
    }

    // Set up 567-588 components of targeted buffer : FH

    auto ta1_z_yyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 567);

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

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_yy_xxxxx_0,  \
                             ta1_z_yy_xxxxx_1,  \
                             ta1_z_yy_xxxxy_0,  \
                             ta1_z_yy_xxxxy_1,  \
                             ta1_z_yy_xxxy_0,   \
                             ta1_z_yy_xxxy_1,   \
                             ta1_z_yy_xxxyy_0,  \
                             ta1_z_yy_xxxyy_1,  \
                             ta1_z_yy_xxxyz_0,  \
                             ta1_z_yy_xxxyz_1,  \
                             ta1_z_yy_xxyy_0,   \
                             ta1_z_yy_xxyy_1,   \
                             ta1_z_yy_xxyyy_0,  \
                             ta1_z_yy_xxyyy_1,  \
                             ta1_z_yy_xxyyz_0,  \
                             ta1_z_yy_xxyyz_1,  \
                             ta1_z_yy_xxyz_0,   \
                             ta1_z_yy_xxyz_1,   \
                             ta1_z_yy_xxyzz_0,  \
                             ta1_z_yy_xxyzz_1,  \
                             ta1_z_yy_xyyy_0,   \
                             ta1_z_yy_xyyy_1,   \
                             ta1_z_yy_xyyyy_0,  \
                             ta1_z_yy_xyyyy_1,  \
                             ta1_z_yy_xyyyz_0,  \
                             ta1_z_yy_xyyyz_1,  \
                             ta1_z_yy_xyyz_0,   \
                             ta1_z_yy_xyyz_1,   \
                             ta1_z_yy_xyyzz_0,  \
                             ta1_z_yy_xyyzz_1,  \
                             ta1_z_yy_xyzz_0,   \
                             ta1_z_yy_xyzz_1,   \
                             ta1_z_yy_xyzzz_0,  \
                             ta1_z_yy_xyzzz_1,  \
                             ta1_z_yy_yyyy_0,   \
                             ta1_z_yy_yyyy_1,   \
                             ta1_z_yy_yyyyy_0,  \
                             ta1_z_yy_yyyyy_1,  \
                             ta1_z_yy_yyyyz_0,  \
                             ta1_z_yy_yyyyz_1,  \
                             ta1_z_yy_yyyz_0,   \
                             ta1_z_yy_yyyz_1,   \
                             ta1_z_yy_yyyzz_0,  \
                             ta1_z_yy_yyyzz_1,  \
                             ta1_z_yy_yyzz_0,   \
                             ta1_z_yy_yyzz_1,   \
                             ta1_z_yy_yyzzz_0,  \
                             ta1_z_yy_yyzzz_1,  \
                             ta1_z_yy_yzzz_0,   \
                             ta1_z_yy_yzzz_1,   \
                             ta1_z_yy_yzzzz_0,  \
                             ta1_z_yy_yzzzz_1,  \
                             ta1_z_yyz_xxxxx_0, \
                             ta1_z_yyz_xxxxy_0, \
                             ta1_z_yyz_xxxxz_0, \
                             ta1_z_yyz_xxxyy_0, \
                             ta1_z_yyz_xxxyz_0, \
                             ta1_z_yyz_xxxzz_0, \
                             ta1_z_yyz_xxyyy_0, \
                             ta1_z_yyz_xxyyz_0, \
                             ta1_z_yyz_xxyzz_0, \
                             ta1_z_yyz_xxzzz_0, \
                             ta1_z_yyz_xyyyy_0, \
                             ta1_z_yyz_xyyyz_0, \
                             ta1_z_yyz_xyyzz_0, \
                             ta1_z_yyz_xyzzz_0, \
                             ta1_z_yyz_xzzzz_0, \
                             ta1_z_yyz_yyyyy_0, \
                             ta1_z_yyz_yyyyz_0, \
                             ta1_z_yyz_yyyzz_0, \
                             ta1_z_yyz_yyzzz_0, \
                             ta1_z_yyz_yzzzz_0, \
                             ta1_z_yyz_zzzzz_0, \
                             ta1_z_yz_xxxxz_0,  \
                             ta1_z_yz_xxxxz_1,  \
                             ta1_z_yz_xxxzz_0,  \
                             ta1_z_yz_xxxzz_1,  \
                             ta1_z_yz_xxzzz_0,  \
                             ta1_z_yz_xxzzz_1,  \
                             ta1_z_yz_xzzzz_0,  \
                             ta1_z_yz_xzzzz_1,  \
                             ta1_z_yz_zzzzz_0,  \
                             ta1_z_yz_zzzzz_1,  \
                             ta1_z_z_xxxxz_0,   \
                             ta1_z_z_xxxxz_1,   \
                             ta1_z_z_xxxzz_0,   \
                             ta1_z_z_xxxzz_1,   \
                             ta1_z_z_xxzzz_0,   \
                             ta1_z_z_xxzzz_1,   \
                             ta1_z_z_xzzzz_0,   \
                             ta1_z_z_xzzzz_1,   \
                             ta1_z_z_zzzzz_0,   \
                             ta1_z_z_zzzzz_1,   \
                             ta_yy_xxxxx_1,     \
                             ta_yy_xxxxy_1,     \
                             ta_yy_xxxyy_1,     \
                             ta_yy_xxxyz_1,     \
                             ta_yy_xxyyy_1,     \
                             ta_yy_xxyyz_1,     \
                             ta_yy_xxyzz_1,     \
                             ta_yy_xyyyy_1,     \
                             ta_yy_xyyyz_1,     \
                             ta_yy_xyyzz_1,     \
                             ta_yy_xyzzz_1,     \
                             ta_yy_yyyyy_1,     \
                             ta_yy_yyyyz_1,     \
                             ta_yy_yyyzz_1,     \
                             ta_yy_yyzzz_1,     \
                             ta_yy_yzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyz_xxxxx_0[i] = ta_yy_xxxxx_1[i] + ta1_z_yy_xxxxx_0[i] * pa_z[i] - ta1_z_yy_xxxxx_1[i] * pc_z[i];

        ta1_z_yyz_xxxxy_0[i] = ta_yy_xxxxy_1[i] + ta1_z_yy_xxxxy_0[i] * pa_z[i] - ta1_z_yy_xxxxy_1[i] * pc_z[i];

        ta1_z_yyz_xxxxz_0[i] = ta1_z_z_xxxxz_0[i] * fe_0 - ta1_z_z_xxxxz_1[i] * fe_0 + ta1_z_yz_xxxxz_0[i] * pa_y[i] - ta1_z_yz_xxxxz_1[i] * pc_y[i];

        ta1_z_yyz_xxxyy_0[i] = ta_yy_xxxyy_1[i] + ta1_z_yy_xxxyy_0[i] * pa_z[i] - ta1_z_yy_xxxyy_1[i] * pc_z[i];

        ta1_z_yyz_xxxyz_0[i] =
            ta1_z_yy_xxxy_0[i] * fe_0 - ta1_z_yy_xxxy_1[i] * fe_0 + ta_yy_xxxyz_1[i] + ta1_z_yy_xxxyz_0[i] * pa_z[i] - ta1_z_yy_xxxyz_1[i] * pc_z[i];

        ta1_z_yyz_xxxzz_0[i] = ta1_z_z_xxxzz_0[i] * fe_0 - ta1_z_z_xxxzz_1[i] * fe_0 + ta1_z_yz_xxxzz_0[i] * pa_y[i] - ta1_z_yz_xxxzz_1[i] * pc_y[i];

        ta1_z_yyz_xxyyy_0[i] = ta_yy_xxyyy_1[i] + ta1_z_yy_xxyyy_0[i] * pa_z[i] - ta1_z_yy_xxyyy_1[i] * pc_z[i];

        ta1_z_yyz_xxyyz_0[i] =
            ta1_z_yy_xxyy_0[i] * fe_0 - ta1_z_yy_xxyy_1[i] * fe_0 + ta_yy_xxyyz_1[i] + ta1_z_yy_xxyyz_0[i] * pa_z[i] - ta1_z_yy_xxyyz_1[i] * pc_z[i];

        ta1_z_yyz_xxyzz_0[i] = 2.0 * ta1_z_yy_xxyz_0[i] * fe_0 - 2.0 * ta1_z_yy_xxyz_1[i] * fe_0 + ta_yy_xxyzz_1[i] + ta1_z_yy_xxyzz_0[i] * pa_z[i] -
                               ta1_z_yy_xxyzz_1[i] * pc_z[i];

        ta1_z_yyz_xxzzz_0[i] = ta1_z_z_xxzzz_0[i] * fe_0 - ta1_z_z_xxzzz_1[i] * fe_0 + ta1_z_yz_xxzzz_0[i] * pa_y[i] - ta1_z_yz_xxzzz_1[i] * pc_y[i];

        ta1_z_yyz_xyyyy_0[i] = ta_yy_xyyyy_1[i] + ta1_z_yy_xyyyy_0[i] * pa_z[i] - ta1_z_yy_xyyyy_1[i] * pc_z[i];

        ta1_z_yyz_xyyyz_0[i] =
            ta1_z_yy_xyyy_0[i] * fe_0 - ta1_z_yy_xyyy_1[i] * fe_0 + ta_yy_xyyyz_1[i] + ta1_z_yy_xyyyz_0[i] * pa_z[i] - ta1_z_yy_xyyyz_1[i] * pc_z[i];

        ta1_z_yyz_xyyzz_0[i] = 2.0 * ta1_z_yy_xyyz_0[i] * fe_0 - 2.0 * ta1_z_yy_xyyz_1[i] * fe_0 + ta_yy_xyyzz_1[i] + ta1_z_yy_xyyzz_0[i] * pa_z[i] -
                               ta1_z_yy_xyyzz_1[i] * pc_z[i];

        ta1_z_yyz_xyzzz_0[i] = 3.0 * ta1_z_yy_xyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xyzz_1[i] * fe_0 + ta_yy_xyzzz_1[i] + ta1_z_yy_xyzzz_0[i] * pa_z[i] -
                               ta1_z_yy_xyzzz_1[i] * pc_z[i];

        ta1_z_yyz_xzzzz_0[i] = ta1_z_z_xzzzz_0[i] * fe_0 - ta1_z_z_xzzzz_1[i] * fe_0 + ta1_z_yz_xzzzz_0[i] * pa_y[i] - ta1_z_yz_xzzzz_1[i] * pc_y[i];

        ta1_z_yyz_yyyyy_0[i] = ta_yy_yyyyy_1[i] + ta1_z_yy_yyyyy_0[i] * pa_z[i] - ta1_z_yy_yyyyy_1[i] * pc_z[i];

        ta1_z_yyz_yyyyz_0[i] =
            ta1_z_yy_yyyy_0[i] * fe_0 - ta1_z_yy_yyyy_1[i] * fe_0 + ta_yy_yyyyz_1[i] + ta1_z_yy_yyyyz_0[i] * pa_z[i] - ta1_z_yy_yyyyz_1[i] * pc_z[i];

        ta1_z_yyz_yyyzz_0[i] = 2.0 * ta1_z_yy_yyyz_0[i] * fe_0 - 2.0 * ta1_z_yy_yyyz_1[i] * fe_0 + ta_yy_yyyzz_1[i] + ta1_z_yy_yyyzz_0[i] * pa_z[i] -
                               ta1_z_yy_yyyzz_1[i] * pc_z[i];

        ta1_z_yyz_yyzzz_0[i] = 3.0 * ta1_z_yy_yyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyzz_1[i] * fe_0 + ta_yy_yyzzz_1[i] + ta1_z_yy_yyzzz_0[i] * pa_z[i] -
                               ta1_z_yy_yyzzz_1[i] * pc_z[i];

        ta1_z_yyz_yzzzz_0[i] = 4.0 * ta1_z_yy_yzzz_0[i] * fe_0 - 4.0 * ta1_z_yy_yzzz_1[i] * fe_0 + ta_yy_yzzzz_1[i] + ta1_z_yy_yzzzz_0[i] * pa_z[i] -
                               ta1_z_yy_yzzzz_1[i] * pc_z[i];

        ta1_z_yyz_zzzzz_0[i] = ta1_z_z_zzzzz_0[i] * fe_0 - ta1_z_z_zzzzz_1[i] * fe_0 + ta1_z_yz_zzzzz_0[i] * pa_y[i] - ta1_z_yz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 588-609 components of targeted buffer : FH

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

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_z_yzz_xxxxx_0, \
                             ta1_z_yzz_xxxxy_0, \
                             ta1_z_yzz_xxxxz_0, \
                             ta1_z_yzz_xxxyy_0, \
                             ta1_z_yzz_xxxyz_0, \
                             ta1_z_yzz_xxxzz_0, \
                             ta1_z_yzz_xxyyy_0, \
                             ta1_z_yzz_xxyyz_0, \
                             ta1_z_yzz_xxyzz_0, \
                             ta1_z_yzz_xxzzz_0, \
                             ta1_z_yzz_xyyyy_0, \
                             ta1_z_yzz_xyyyz_0, \
                             ta1_z_yzz_xyyzz_0, \
                             ta1_z_yzz_xyzzz_0, \
                             ta1_z_yzz_xzzzz_0, \
                             ta1_z_yzz_yyyyy_0, \
                             ta1_z_yzz_yyyyz_0, \
                             ta1_z_yzz_yyyzz_0, \
                             ta1_z_yzz_yyzzz_0, \
                             ta1_z_yzz_yzzzz_0, \
                             ta1_z_yzz_zzzzz_0, \
                             ta1_z_zz_xxxx_0,   \
                             ta1_z_zz_xxxx_1,   \
                             ta1_z_zz_xxxxx_0,  \
                             ta1_z_zz_xxxxx_1,  \
                             ta1_z_zz_xxxxy_0,  \
                             ta1_z_zz_xxxxy_1,  \
                             ta1_z_zz_xxxxz_0,  \
                             ta1_z_zz_xxxxz_1,  \
                             ta1_z_zz_xxxy_0,   \
                             ta1_z_zz_xxxy_1,   \
                             ta1_z_zz_xxxyy_0,  \
                             ta1_z_zz_xxxyy_1,  \
                             ta1_z_zz_xxxyz_0,  \
                             ta1_z_zz_xxxyz_1,  \
                             ta1_z_zz_xxxz_0,   \
                             ta1_z_zz_xxxz_1,   \
                             ta1_z_zz_xxxzz_0,  \
                             ta1_z_zz_xxxzz_1,  \
                             ta1_z_zz_xxyy_0,   \
                             ta1_z_zz_xxyy_1,   \
                             ta1_z_zz_xxyyy_0,  \
                             ta1_z_zz_xxyyy_1,  \
                             ta1_z_zz_xxyyz_0,  \
                             ta1_z_zz_xxyyz_1,  \
                             ta1_z_zz_xxyz_0,   \
                             ta1_z_zz_xxyz_1,   \
                             ta1_z_zz_xxyzz_0,  \
                             ta1_z_zz_xxyzz_1,  \
                             ta1_z_zz_xxzz_0,   \
                             ta1_z_zz_xxzz_1,   \
                             ta1_z_zz_xxzzz_0,  \
                             ta1_z_zz_xxzzz_1,  \
                             ta1_z_zz_xyyy_0,   \
                             ta1_z_zz_xyyy_1,   \
                             ta1_z_zz_xyyyy_0,  \
                             ta1_z_zz_xyyyy_1,  \
                             ta1_z_zz_xyyyz_0,  \
                             ta1_z_zz_xyyyz_1,  \
                             ta1_z_zz_xyyz_0,   \
                             ta1_z_zz_xyyz_1,   \
                             ta1_z_zz_xyyzz_0,  \
                             ta1_z_zz_xyyzz_1,  \
                             ta1_z_zz_xyzz_0,   \
                             ta1_z_zz_xyzz_1,   \
                             ta1_z_zz_xyzzz_0,  \
                             ta1_z_zz_xyzzz_1,  \
                             ta1_z_zz_xzzz_0,   \
                             ta1_z_zz_xzzz_1,   \
                             ta1_z_zz_xzzzz_0,  \
                             ta1_z_zz_xzzzz_1,  \
                             ta1_z_zz_yyyy_0,   \
                             ta1_z_zz_yyyy_1,   \
                             ta1_z_zz_yyyyy_0,  \
                             ta1_z_zz_yyyyy_1,  \
                             ta1_z_zz_yyyyz_0,  \
                             ta1_z_zz_yyyyz_1,  \
                             ta1_z_zz_yyyz_0,   \
                             ta1_z_zz_yyyz_1,   \
                             ta1_z_zz_yyyzz_0,  \
                             ta1_z_zz_yyyzz_1,  \
                             ta1_z_zz_yyzz_0,   \
                             ta1_z_zz_yyzz_1,   \
                             ta1_z_zz_yyzzz_0,  \
                             ta1_z_zz_yyzzz_1,  \
                             ta1_z_zz_yzzz_0,   \
                             ta1_z_zz_yzzz_1,   \
                             ta1_z_zz_yzzzz_0,  \
                             ta1_z_zz_yzzzz_1,  \
                             ta1_z_zz_zzzz_0,   \
                             ta1_z_zz_zzzz_1,   \
                             ta1_z_zz_zzzzz_0,  \
                             ta1_z_zz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzz_xxxxx_0[i] = ta1_z_zz_xxxxx_0[i] * pa_y[i] - ta1_z_zz_xxxxx_1[i] * pc_y[i];

        ta1_z_yzz_xxxxy_0[i] = ta1_z_zz_xxxx_0[i] * fe_0 - ta1_z_zz_xxxx_1[i] * fe_0 + ta1_z_zz_xxxxy_0[i] * pa_y[i] - ta1_z_zz_xxxxy_1[i] * pc_y[i];

        ta1_z_yzz_xxxxz_0[i] = ta1_z_zz_xxxxz_0[i] * pa_y[i] - ta1_z_zz_xxxxz_1[i] * pc_y[i];

        ta1_z_yzz_xxxyy_0[i] =
            2.0 * ta1_z_zz_xxxy_0[i] * fe_0 - 2.0 * ta1_z_zz_xxxy_1[i] * fe_0 + ta1_z_zz_xxxyy_0[i] * pa_y[i] - ta1_z_zz_xxxyy_1[i] * pc_y[i];

        ta1_z_yzz_xxxyz_0[i] = ta1_z_zz_xxxz_0[i] * fe_0 - ta1_z_zz_xxxz_1[i] * fe_0 + ta1_z_zz_xxxyz_0[i] * pa_y[i] - ta1_z_zz_xxxyz_1[i] * pc_y[i];

        ta1_z_yzz_xxxzz_0[i] = ta1_z_zz_xxxzz_0[i] * pa_y[i] - ta1_z_zz_xxxzz_1[i] * pc_y[i];

        ta1_z_yzz_xxyyy_0[i] =
            3.0 * ta1_z_zz_xxyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyy_1[i] * fe_0 + ta1_z_zz_xxyyy_0[i] * pa_y[i] - ta1_z_zz_xxyyy_1[i] * pc_y[i];

        ta1_z_yzz_xxyyz_0[i] =
            2.0 * ta1_z_zz_xxyz_0[i] * fe_0 - 2.0 * ta1_z_zz_xxyz_1[i] * fe_0 + ta1_z_zz_xxyyz_0[i] * pa_y[i] - ta1_z_zz_xxyyz_1[i] * pc_y[i];

        ta1_z_yzz_xxyzz_0[i] = ta1_z_zz_xxzz_0[i] * fe_0 - ta1_z_zz_xxzz_1[i] * fe_0 + ta1_z_zz_xxyzz_0[i] * pa_y[i] - ta1_z_zz_xxyzz_1[i] * pc_y[i];

        ta1_z_yzz_xxzzz_0[i] = ta1_z_zz_xxzzz_0[i] * pa_y[i] - ta1_z_zz_xxzzz_1[i] * pc_y[i];

        ta1_z_yzz_xyyyy_0[i] =
            4.0 * ta1_z_zz_xyyy_0[i] * fe_0 - 4.0 * ta1_z_zz_xyyy_1[i] * fe_0 + ta1_z_zz_xyyyy_0[i] * pa_y[i] - ta1_z_zz_xyyyy_1[i] * pc_y[i];

        ta1_z_yzz_xyyyz_0[i] =
            3.0 * ta1_z_zz_xyyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xyyz_1[i] * fe_0 + ta1_z_zz_xyyyz_0[i] * pa_y[i] - ta1_z_zz_xyyyz_1[i] * pc_y[i];

        ta1_z_yzz_xyyzz_0[i] =
            2.0 * ta1_z_zz_xyzz_0[i] * fe_0 - 2.0 * ta1_z_zz_xyzz_1[i] * fe_0 + ta1_z_zz_xyyzz_0[i] * pa_y[i] - ta1_z_zz_xyyzz_1[i] * pc_y[i];

        ta1_z_yzz_xyzzz_0[i] = ta1_z_zz_xzzz_0[i] * fe_0 - ta1_z_zz_xzzz_1[i] * fe_0 + ta1_z_zz_xyzzz_0[i] * pa_y[i] - ta1_z_zz_xyzzz_1[i] * pc_y[i];

        ta1_z_yzz_xzzzz_0[i] = ta1_z_zz_xzzzz_0[i] * pa_y[i] - ta1_z_zz_xzzzz_1[i] * pc_y[i];

        ta1_z_yzz_yyyyy_0[i] =
            5.0 * ta1_z_zz_yyyy_0[i] * fe_0 - 5.0 * ta1_z_zz_yyyy_1[i] * fe_0 + ta1_z_zz_yyyyy_0[i] * pa_y[i] - ta1_z_zz_yyyyy_1[i] * pc_y[i];

        ta1_z_yzz_yyyyz_0[i] =
            4.0 * ta1_z_zz_yyyz_0[i] * fe_0 - 4.0 * ta1_z_zz_yyyz_1[i] * fe_0 + ta1_z_zz_yyyyz_0[i] * pa_y[i] - ta1_z_zz_yyyyz_1[i] * pc_y[i];

        ta1_z_yzz_yyyzz_0[i] =
            3.0 * ta1_z_zz_yyzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyzz_1[i] * fe_0 + ta1_z_zz_yyyzz_0[i] * pa_y[i] - ta1_z_zz_yyyzz_1[i] * pc_y[i];

        ta1_z_yzz_yyzzz_0[i] =
            2.0 * ta1_z_zz_yzzz_0[i] * fe_0 - 2.0 * ta1_z_zz_yzzz_1[i] * fe_0 + ta1_z_zz_yyzzz_0[i] * pa_y[i] - ta1_z_zz_yyzzz_1[i] * pc_y[i];

        ta1_z_yzz_yzzzz_0[i] = ta1_z_zz_zzzz_0[i] * fe_0 - ta1_z_zz_zzzz_1[i] * fe_0 + ta1_z_zz_yzzzz_0[i] * pa_y[i] - ta1_z_zz_yzzzz_1[i] * pc_y[i];

        ta1_z_yzz_zzzzz_0[i] = ta1_z_zz_zzzzz_0[i] * pa_y[i] - ta1_z_zz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 609-630 components of targeted buffer : FH

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

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_z_z_xxxxx_0,   \
                             ta1_z_z_xxxxx_1,   \
                             ta1_z_z_xxxxy_0,   \
                             ta1_z_z_xxxxy_1,   \
                             ta1_z_z_xxxxz_0,   \
                             ta1_z_z_xxxxz_1,   \
                             ta1_z_z_xxxyy_0,   \
                             ta1_z_z_xxxyy_1,   \
                             ta1_z_z_xxxyz_0,   \
                             ta1_z_z_xxxyz_1,   \
                             ta1_z_z_xxxzz_0,   \
                             ta1_z_z_xxxzz_1,   \
                             ta1_z_z_xxyyy_0,   \
                             ta1_z_z_xxyyy_1,   \
                             ta1_z_z_xxyyz_0,   \
                             ta1_z_z_xxyyz_1,   \
                             ta1_z_z_xxyzz_0,   \
                             ta1_z_z_xxyzz_1,   \
                             ta1_z_z_xxzzz_0,   \
                             ta1_z_z_xxzzz_1,   \
                             ta1_z_z_xyyyy_0,   \
                             ta1_z_z_xyyyy_1,   \
                             ta1_z_z_xyyyz_0,   \
                             ta1_z_z_xyyyz_1,   \
                             ta1_z_z_xyyzz_0,   \
                             ta1_z_z_xyyzz_1,   \
                             ta1_z_z_xyzzz_0,   \
                             ta1_z_z_xyzzz_1,   \
                             ta1_z_z_xzzzz_0,   \
                             ta1_z_z_xzzzz_1,   \
                             ta1_z_z_yyyyy_0,   \
                             ta1_z_z_yyyyy_1,   \
                             ta1_z_z_yyyyz_0,   \
                             ta1_z_z_yyyyz_1,   \
                             ta1_z_z_yyyzz_0,   \
                             ta1_z_z_yyyzz_1,   \
                             ta1_z_z_yyzzz_0,   \
                             ta1_z_z_yyzzz_1,   \
                             ta1_z_z_yzzzz_0,   \
                             ta1_z_z_yzzzz_1,   \
                             ta1_z_z_zzzzz_0,   \
                             ta1_z_z_zzzzz_1,   \
                             ta1_z_zz_xxxx_0,   \
                             ta1_z_zz_xxxx_1,   \
                             ta1_z_zz_xxxxx_0,  \
                             ta1_z_zz_xxxxx_1,  \
                             ta1_z_zz_xxxxy_0,  \
                             ta1_z_zz_xxxxy_1,  \
                             ta1_z_zz_xxxxz_0,  \
                             ta1_z_zz_xxxxz_1,  \
                             ta1_z_zz_xxxy_0,   \
                             ta1_z_zz_xxxy_1,   \
                             ta1_z_zz_xxxyy_0,  \
                             ta1_z_zz_xxxyy_1,  \
                             ta1_z_zz_xxxyz_0,  \
                             ta1_z_zz_xxxyz_1,  \
                             ta1_z_zz_xxxz_0,   \
                             ta1_z_zz_xxxz_1,   \
                             ta1_z_zz_xxxzz_0,  \
                             ta1_z_zz_xxxzz_1,  \
                             ta1_z_zz_xxyy_0,   \
                             ta1_z_zz_xxyy_1,   \
                             ta1_z_zz_xxyyy_0,  \
                             ta1_z_zz_xxyyy_1,  \
                             ta1_z_zz_xxyyz_0,  \
                             ta1_z_zz_xxyyz_1,  \
                             ta1_z_zz_xxyz_0,   \
                             ta1_z_zz_xxyz_1,   \
                             ta1_z_zz_xxyzz_0,  \
                             ta1_z_zz_xxyzz_1,  \
                             ta1_z_zz_xxzz_0,   \
                             ta1_z_zz_xxzz_1,   \
                             ta1_z_zz_xxzzz_0,  \
                             ta1_z_zz_xxzzz_1,  \
                             ta1_z_zz_xyyy_0,   \
                             ta1_z_zz_xyyy_1,   \
                             ta1_z_zz_xyyyy_0,  \
                             ta1_z_zz_xyyyy_1,  \
                             ta1_z_zz_xyyyz_0,  \
                             ta1_z_zz_xyyyz_1,  \
                             ta1_z_zz_xyyz_0,   \
                             ta1_z_zz_xyyz_1,   \
                             ta1_z_zz_xyyzz_0,  \
                             ta1_z_zz_xyyzz_1,  \
                             ta1_z_zz_xyzz_0,   \
                             ta1_z_zz_xyzz_1,   \
                             ta1_z_zz_xyzzz_0,  \
                             ta1_z_zz_xyzzz_1,  \
                             ta1_z_zz_xzzz_0,   \
                             ta1_z_zz_xzzz_1,   \
                             ta1_z_zz_xzzzz_0,  \
                             ta1_z_zz_xzzzz_1,  \
                             ta1_z_zz_yyyy_0,   \
                             ta1_z_zz_yyyy_1,   \
                             ta1_z_zz_yyyyy_0,  \
                             ta1_z_zz_yyyyy_1,  \
                             ta1_z_zz_yyyyz_0,  \
                             ta1_z_zz_yyyyz_1,  \
                             ta1_z_zz_yyyz_0,   \
                             ta1_z_zz_yyyz_1,   \
                             ta1_z_zz_yyyzz_0,  \
                             ta1_z_zz_yyyzz_1,  \
                             ta1_z_zz_yyzz_0,   \
                             ta1_z_zz_yyzz_1,   \
                             ta1_z_zz_yyzzz_0,  \
                             ta1_z_zz_yyzzz_1,  \
                             ta1_z_zz_yzzz_0,   \
                             ta1_z_zz_yzzz_1,   \
                             ta1_z_zz_yzzzz_0,  \
                             ta1_z_zz_yzzzz_1,  \
                             ta1_z_zz_zzzz_0,   \
                             ta1_z_zz_zzzz_1,   \
                             ta1_z_zz_zzzzz_0,  \
                             ta1_z_zz_zzzzz_1,  \
                             ta1_z_zzz_xxxxx_0, \
                             ta1_z_zzz_xxxxy_0, \
                             ta1_z_zzz_xxxxz_0, \
                             ta1_z_zzz_xxxyy_0, \
                             ta1_z_zzz_xxxyz_0, \
                             ta1_z_zzz_xxxzz_0, \
                             ta1_z_zzz_xxyyy_0, \
                             ta1_z_zzz_xxyyz_0, \
                             ta1_z_zzz_xxyzz_0, \
                             ta1_z_zzz_xxzzz_0, \
                             ta1_z_zzz_xyyyy_0, \
                             ta1_z_zzz_xyyyz_0, \
                             ta1_z_zzz_xyyzz_0, \
                             ta1_z_zzz_xyzzz_0, \
                             ta1_z_zzz_xzzzz_0, \
                             ta1_z_zzz_yyyyy_0, \
                             ta1_z_zzz_yyyyz_0, \
                             ta1_z_zzz_yyyzz_0, \
                             ta1_z_zzz_yyzzz_0, \
                             ta1_z_zzz_yzzzz_0, \
                             ta1_z_zzz_zzzzz_0, \
                             ta_zz_xxxxx_1,     \
                             ta_zz_xxxxy_1,     \
                             ta_zz_xxxxz_1,     \
                             ta_zz_xxxyy_1,     \
                             ta_zz_xxxyz_1,     \
                             ta_zz_xxxzz_1,     \
                             ta_zz_xxyyy_1,     \
                             ta_zz_xxyyz_1,     \
                             ta_zz_xxyzz_1,     \
                             ta_zz_xxzzz_1,     \
                             ta_zz_xyyyy_1,     \
                             ta_zz_xyyyz_1,     \
                             ta_zz_xyyzz_1,     \
                             ta_zz_xyzzz_1,     \
                             ta_zz_xzzzz_1,     \
                             ta_zz_yyyyy_1,     \
                             ta_zz_yyyyz_1,     \
                             ta_zz_yyyzz_1,     \
                             ta_zz_yyzzz_1,     \
                             ta_zz_yzzzz_1,     \
                             ta_zz_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzz_xxxxx_0[i] = 2.0 * ta1_z_z_xxxxx_0[i] * fe_0 - 2.0 * ta1_z_z_xxxxx_1[i] * fe_0 + ta_zz_xxxxx_1[i] + ta1_z_zz_xxxxx_0[i] * pa_z[i] -
                               ta1_z_zz_xxxxx_1[i] * pc_z[i];

        ta1_z_zzz_xxxxy_0[i] = 2.0 * ta1_z_z_xxxxy_0[i] * fe_0 - 2.0 * ta1_z_z_xxxxy_1[i] * fe_0 + ta_zz_xxxxy_1[i] + ta1_z_zz_xxxxy_0[i] * pa_z[i] -
                               ta1_z_zz_xxxxy_1[i] * pc_z[i];

        ta1_z_zzz_xxxxz_0[i] = 2.0 * ta1_z_z_xxxxz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxxz_1[i] * fe_0 + ta1_z_zz_xxxx_0[i] * fe_0 -
                               ta1_z_zz_xxxx_1[i] * fe_0 + ta_zz_xxxxz_1[i] + ta1_z_zz_xxxxz_0[i] * pa_z[i] - ta1_z_zz_xxxxz_1[i] * pc_z[i];

        ta1_z_zzz_xxxyy_0[i] = 2.0 * ta1_z_z_xxxyy_0[i] * fe_0 - 2.0 * ta1_z_z_xxxyy_1[i] * fe_0 + ta_zz_xxxyy_1[i] + ta1_z_zz_xxxyy_0[i] * pa_z[i] -
                               ta1_z_zz_xxxyy_1[i] * pc_z[i];

        ta1_z_zzz_xxxyz_0[i] = 2.0 * ta1_z_z_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxyz_1[i] * fe_0 + ta1_z_zz_xxxy_0[i] * fe_0 -
                               ta1_z_zz_xxxy_1[i] * fe_0 + ta_zz_xxxyz_1[i] + ta1_z_zz_xxxyz_0[i] * pa_z[i] - ta1_z_zz_xxxyz_1[i] * pc_z[i];

        ta1_z_zzz_xxxzz_0[i] = 2.0 * ta1_z_z_xxxzz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxzz_1[i] * fe_0 + 2.0 * ta1_z_zz_xxxz_0[i] * fe_0 -
                               2.0 * ta1_z_zz_xxxz_1[i] * fe_0 + ta_zz_xxxzz_1[i] + ta1_z_zz_xxxzz_0[i] * pa_z[i] - ta1_z_zz_xxxzz_1[i] * pc_z[i];

        ta1_z_zzz_xxyyy_0[i] = 2.0 * ta1_z_z_xxyyy_0[i] * fe_0 - 2.0 * ta1_z_z_xxyyy_1[i] * fe_0 + ta_zz_xxyyy_1[i] + ta1_z_zz_xxyyy_0[i] * pa_z[i] -
                               ta1_z_zz_xxyyy_1[i] * pc_z[i];

        ta1_z_zzz_xxyyz_0[i] = 2.0 * ta1_z_z_xxyyz_0[i] * fe_0 - 2.0 * ta1_z_z_xxyyz_1[i] * fe_0 + ta1_z_zz_xxyy_0[i] * fe_0 -
                               ta1_z_zz_xxyy_1[i] * fe_0 + ta_zz_xxyyz_1[i] + ta1_z_zz_xxyyz_0[i] * pa_z[i] - ta1_z_zz_xxyyz_1[i] * pc_z[i];

        ta1_z_zzz_xxyzz_0[i] = 2.0 * ta1_z_z_xxyzz_0[i] * fe_0 - 2.0 * ta1_z_z_xxyzz_1[i] * fe_0 + 2.0 * ta1_z_zz_xxyz_0[i] * fe_0 -
                               2.0 * ta1_z_zz_xxyz_1[i] * fe_0 + ta_zz_xxyzz_1[i] + ta1_z_zz_xxyzz_0[i] * pa_z[i] - ta1_z_zz_xxyzz_1[i] * pc_z[i];

        ta1_z_zzz_xxzzz_0[i] = 2.0 * ta1_z_z_xxzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xxzzz_1[i] * fe_0 + 3.0 * ta1_z_zz_xxzz_0[i] * fe_0 -
                               3.0 * ta1_z_zz_xxzz_1[i] * fe_0 + ta_zz_xxzzz_1[i] + ta1_z_zz_xxzzz_0[i] * pa_z[i] - ta1_z_zz_xxzzz_1[i] * pc_z[i];

        ta1_z_zzz_xyyyy_0[i] = 2.0 * ta1_z_z_xyyyy_0[i] * fe_0 - 2.0 * ta1_z_z_xyyyy_1[i] * fe_0 + ta_zz_xyyyy_1[i] + ta1_z_zz_xyyyy_0[i] * pa_z[i] -
                               ta1_z_zz_xyyyy_1[i] * pc_z[i];

        ta1_z_zzz_xyyyz_0[i] = 2.0 * ta1_z_z_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_z_xyyyz_1[i] * fe_0 + ta1_z_zz_xyyy_0[i] * fe_0 -
                               ta1_z_zz_xyyy_1[i] * fe_0 + ta_zz_xyyyz_1[i] + ta1_z_zz_xyyyz_0[i] * pa_z[i] - ta1_z_zz_xyyyz_1[i] * pc_z[i];

        ta1_z_zzz_xyyzz_0[i] = 2.0 * ta1_z_z_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_z_xyyzz_1[i] * fe_0 + 2.0 * ta1_z_zz_xyyz_0[i] * fe_0 -
                               2.0 * ta1_z_zz_xyyz_1[i] * fe_0 + ta_zz_xyyzz_1[i] + ta1_z_zz_xyyzz_0[i] * pa_z[i] - ta1_z_zz_xyyzz_1[i] * pc_z[i];

        ta1_z_zzz_xyzzz_0[i] = 2.0 * ta1_z_z_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xyzzz_1[i] * fe_0 + 3.0 * ta1_z_zz_xyzz_0[i] * fe_0 -
                               3.0 * ta1_z_zz_xyzz_1[i] * fe_0 + ta_zz_xyzzz_1[i] + ta1_z_zz_xyzzz_0[i] * pa_z[i] - ta1_z_zz_xyzzz_1[i] * pc_z[i];

        ta1_z_zzz_xzzzz_0[i] = 2.0 * ta1_z_z_xzzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xzzzz_1[i] * fe_0 + 4.0 * ta1_z_zz_xzzz_0[i] * fe_0 -
                               4.0 * ta1_z_zz_xzzz_1[i] * fe_0 + ta_zz_xzzzz_1[i] + ta1_z_zz_xzzzz_0[i] * pa_z[i] - ta1_z_zz_xzzzz_1[i] * pc_z[i];

        ta1_z_zzz_yyyyy_0[i] = 2.0 * ta1_z_z_yyyyy_0[i] * fe_0 - 2.0 * ta1_z_z_yyyyy_1[i] * fe_0 + ta_zz_yyyyy_1[i] + ta1_z_zz_yyyyy_0[i] * pa_z[i] -
                               ta1_z_zz_yyyyy_1[i] * pc_z[i];

        ta1_z_zzz_yyyyz_0[i] = 2.0 * ta1_z_z_yyyyz_0[i] * fe_0 - 2.0 * ta1_z_z_yyyyz_1[i] * fe_0 + ta1_z_zz_yyyy_0[i] * fe_0 -
                               ta1_z_zz_yyyy_1[i] * fe_0 + ta_zz_yyyyz_1[i] + ta1_z_zz_yyyyz_0[i] * pa_z[i] - ta1_z_zz_yyyyz_1[i] * pc_z[i];

        ta1_z_zzz_yyyzz_0[i] = 2.0 * ta1_z_z_yyyzz_0[i] * fe_0 - 2.0 * ta1_z_z_yyyzz_1[i] * fe_0 + 2.0 * ta1_z_zz_yyyz_0[i] * fe_0 -
                               2.0 * ta1_z_zz_yyyz_1[i] * fe_0 + ta_zz_yyyzz_1[i] + ta1_z_zz_yyyzz_0[i] * pa_z[i] - ta1_z_zz_yyyzz_1[i] * pc_z[i];

        ta1_z_zzz_yyzzz_0[i] = 2.0 * ta1_z_z_yyzzz_0[i] * fe_0 - 2.0 * ta1_z_z_yyzzz_1[i] * fe_0 + 3.0 * ta1_z_zz_yyzz_0[i] * fe_0 -
                               3.0 * ta1_z_zz_yyzz_1[i] * fe_0 + ta_zz_yyzzz_1[i] + ta1_z_zz_yyzzz_0[i] * pa_z[i] - ta1_z_zz_yyzzz_1[i] * pc_z[i];

        ta1_z_zzz_yzzzz_0[i] = 2.0 * ta1_z_z_yzzzz_0[i] * fe_0 - 2.0 * ta1_z_z_yzzzz_1[i] * fe_0 + 4.0 * ta1_z_zz_yzzz_0[i] * fe_0 -
                               4.0 * ta1_z_zz_yzzz_1[i] * fe_0 + ta_zz_yzzzz_1[i] + ta1_z_zz_yzzzz_0[i] * pa_z[i] - ta1_z_zz_yzzzz_1[i] * pc_z[i];

        ta1_z_zzz_zzzzz_0[i] = 2.0 * ta1_z_z_zzzzz_0[i] * fe_0 - 2.0 * ta1_z_z_zzzzz_1[i] * fe_0 + 5.0 * ta1_z_zz_zzzz_0[i] * fe_0 -
                               5.0 * ta1_z_zz_zzzz_1[i] * fe_0 + ta_zz_zzzzz_1[i] + ta1_z_zz_zzzzz_0[i] * pa_z[i] - ta1_z_zz_zzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
