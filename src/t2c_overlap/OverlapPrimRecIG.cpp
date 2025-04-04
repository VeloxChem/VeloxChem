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

#include "OverlapPrimRecIG.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_ig(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_ig,
                     const size_t idx_ovl_gg,
                     const size_t idx_ovl_hf,
                     const size_t idx_ovl_hg,
                     const CSimdArray<double>& factors,
                     const size_t idx_rpa,
                     const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

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

    auto ts_xxxy_xxxx = pbuffer.data(idx_ovl_gg + 15);

    auto ts_xxxy_xxxz = pbuffer.data(idx_ovl_gg + 17);

    auto ts_xxxy_xxzz = pbuffer.data(idx_ovl_gg + 20);

    auto ts_xxxy_xzzz = pbuffer.data(idx_ovl_gg + 24);

    auto ts_xxxy_yyyy = pbuffer.data(idx_ovl_gg + 25);

    auto ts_xxxy_yyyz = pbuffer.data(idx_ovl_gg + 26);

    auto ts_xxxy_yyzz = pbuffer.data(idx_ovl_gg + 27);

    auto ts_xxxy_yzzz = pbuffer.data(idx_ovl_gg + 28);

    auto ts_xxxz_xxxx = pbuffer.data(idx_ovl_gg + 30);

    auto ts_xxxz_xxxy = pbuffer.data(idx_ovl_gg + 31);

    auto ts_xxxz_xxxz = pbuffer.data(idx_ovl_gg + 32);

    auto ts_xxxz_xxyy = pbuffer.data(idx_ovl_gg + 33);

    auto ts_xxxz_xxzz = pbuffer.data(idx_ovl_gg + 35);

    auto ts_xxxz_xyyy = pbuffer.data(idx_ovl_gg + 36);

    auto ts_xxxz_xzzz = pbuffer.data(idx_ovl_gg + 39);

    auto ts_xxxz_yyyz = pbuffer.data(idx_ovl_gg + 41);

    auto ts_xxxz_yyzz = pbuffer.data(idx_ovl_gg + 42);

    auto ts_xxxz_yzzz = pbuffer.data(idx_ovl_gg + 43);

    auto ts_xxxz_zzzz = pbuffer.data(idx_ovl_gg + 44);

    auto ts_xxyy_xxxx = pbuffer.data(idx_ovl_gg + 45);

    auto ts_xxyy_xxxy = pbuffer.data(idx_ovl_gg + 46);

    auto ts_xxyy_xxxz = pbuffer.data(idx_ovl_gg + 47);

    auto ts_xxyy_xxyy = pbuffer.data(idx_ovl_gg + 48);

    auto ts_xxyy_xxyz = pbuffer.data(idx_ovl_gg + 49);

    auto ts_xxyy_xxzz = pbuffer.data(idx_ovl_gg + 50);

    auto ts_xxyy_xyyy = pbuffer.data(idx_ovl_gg + 51);

    auto ts_xxyy_xyyz = pbuffer.data(idx_ovl_gg + 52);

    auto ts_xxyy_xyzz = pbuffer.data(idx_ovl_gg + 53);

    auto ts_xxyy_xzzz = pbuffer.data(idx_ovl_gg + 54);

    auto ts_xxyy_yyyy = pbuffer.data(idx_ovl_gg + 55);

    auto ts_xxyy_yyyz = pbuffer.data(idx_ovl_gg + 56);

    auto ts_xxyy_yyzz = pbuffer.data(idx_ovl_gg + 57);

    auto ts_xxyy_yzzz = pbuffer.data(idx_ovl_gg + 58);

    auto ts_xxyy_zzzz = pbuffer.data(idx_ovl_gg + 59);

    auto ts_xxyz_xxxz = pbuffer.data(idx_ovl_gg + 62);

    auto ts_xxyz_xxzz = pbuffer.data(idx_ovl_gg + 65);

    auto ts_xxyz_xzzz = pbuffer.data(idx_ovl_gg + 69);

    auto ts_xxyz_yyyz = pbuffer.data(idx_ovl_gg + 71);

    auto ts_xxyz_yyzz = pbuffer.data(idx_ovl_gg + 72);

    auto ts_xxyz_yzzz = pbuffer.data(idx_ovl_gg + 73);

    auto ts_xxzz_xxxx = pbuffer.data(idx_ovl_gg + 75);

    auto ts_xxzz_xxxy = pbuffer.data(idx_ovl_gg + 76);

    auto ts_xxzz_xxxz = pbuffer.data(idx_ovl_gg + 77);

    auto ts_xxzz_xxyy = pbuffer.data(idx_ovl_gg + 78);

    auto ts_xxzz_xxyz = pbuffer.data(idx_ovl_gg + 79);

    auto ts_xxzz_xxzz = pbuffer.data(idx_ovl_gg + 80);

    auto ts_xxzz_xyyy = pbuffer.data(idx_ovl_gg + 81);

    auto ts_xxzz_xyyz = pbuffer.data(idx_ovl_gg + 82);

    auto ts_xxzz_xyzz = pbuffer.data(idx_ovl_gg + 83);

    auto ts_xxzz_xzzz = pbuffer.data(idx_ovl_gg + 84);

    auto ts_xxzz_yyyy = pbuffer.data(idx_ovl_gg + 85);

    auto ts_xxzz_yyyz = pbuffer.data(idx_ovl_gg + 86);

    auto ts_xxzz_yyzz = pbuffer.data(idx_ovl_gg + 87);

    auto ts_xxzz_yzzz = pbuffer.data(idx_ovl_gg + 88);

    auto ts_xxzz_zzzz = pbuffer.data(idx_ovl_gg + 89);

    auto ts_xyyy_xxxy = pbuffer.data(idx_ovl_gg + 91);

    auto ts_xyyy_xxyy = pbuffer.data(idx_ovl_gg + 93);

    auto ts_xyyy_xxyz = pbuffer.data(idx_ovl_gg + 94);

    auto ts_xyyy_xyyy = pbuffer.data(idx_ovl_gg + 96);

    auto ts_xyyy_xyyz = pbuffer.data(idx_ovl_gg + 97);

    auto ts_xyyy_xyzz = pbuffer.data(idx_ovl_gg + 98);

    auto ts_xyyy_yyyy = pbuffer.data(idx_ovl_gg + 100);

    auto ts_xyyy_yyyz = pbuffer.data(idx_ovl_gg + 101);

    auto ts_xyyy_yyzz = pbuffer.data(idx_ovl_gg + 102);

    auto ts_xyyy_yzzz = pbuffer.data(idx_ovl_gg + 103);

    auto ts_xyyy_zzzz = pbuffer.data(idx_ovl_gg + 104);

    auto ts_xyyz_yyyz = pbuffer.data(idx_ovl_gg + 116);

    auto ts_xyyz_yyzz = pbuffer.data(idx_ovl_gg + 117);

    auto ts_xyyz_yzzz = pbuffer.data(idx_ovl_gg + 118);

    auto ts_xyyz_zzzz = pbuffer.data(idx_ovl_gg + 119);

    auto ts_xyzz_yyyy = pbuffer.data(idx_ovl_gg + 130);

    auto ts_xyzz_yyyz = pbuffer.data(idx_ovl_gg + 131);

    auto ts_xyzz_yyzz = pbuffer.data(idx_ovl_gg + 132);

    auto ts_xyzz_yzzz = pbuffer.data(idx_ovl_gg + 133);

    auto ts_xzzz_xxxz = pbuffer.data(idx_ovl_gg + 137);

    auto ts_xzzz_xxyz = pbuffer.data(idx_ovl_gg + 139);

    auto ts_xzzz_xxzz = pbuffer.data(idx_ovl_gg + 140);

    auto ts_xzzz_xyyz = pbuffer.data(idx_ovl_gg + 142);

    auto ts_xzzz_xyzz = pbuffer.data(idx_ovl_gg + 143);

    auto ts_xzzz_xzzz = pbuffer.data(idx_ovl_gg + 144);

    auto ts_xzzz_yyyy = pbuffer.data(idx_ovl_gg + 145);

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

    auto ts_yyyz_xxxy = pbuffer.data(idx_ovl_gg + 166);

    auto ts_yyyz_xxxz = pbuffer.data(idx_ovl_gg + 167);

    auto ts_yyyz_xxyy = pbuffer.data(idx_ovl_gg + 168);

    auto ts_yyyz_xxzz = pbuffer.data(idx_ovl_gg + 170);

    auto ts_yyyz_xyyy = pbuffer.data(idx_ovl_gg + 171);

    auto ts_yyyz_xzzz = pbuffer.data(idx_ovl_gg + 174);

    auto ts_yyyz_yyyy = pbuffer.data(idx_ovl_gg + 175);

    auto ts_yyyz_yyyz = pbuffer.data(idx_ovl_gg + 176);

    auto ts_yyyz_yyzz = pbuffer.data(idx_ovl_gg + 177);

    auto ts_yyyz_yzzz = pbuffer.data(idx_ovl_gg + 178);

    auto ts_yyyz_zzzz = pbuffer.data(idx_ovl_gg + 179);

    auto ts_yyzz_xxxx = pbuffer.data(idx_ovl_gg + 180);

    auto ts_yyzz_xxxy = pbuffer.data(idx_ovl_gg + 181);

    auto ts_yyzz_xxxz = pbuffer.data(idx_ovl_gg + 182);

    auto ts_yyzz_xxyy = pbuffer.data(idx_ovl_gg + 183);

    auto ts_yyzz_xxyz = pbuffer.data(idx_ovl_gg + 184);

    auto ts_yyzz_xxzz = pbuffer.data(idx_ovl_gg + 185);

    auto ts_yyzz_xyyy = pbuffer.data(idx_ovl_gg + 186);

    auto ts_yyzz_xyyz = pbuffer.data(idx_ovl_gg + 187);

    auto ts_yyzz_xyzz = pbuffer.data(idx_ovl_gg + 188);

    auto ts_yyzz_xzzz = pbuffer.data(idx_ovl_gg + 189);

    auto ts_yyzz_yyyy = pbuffer.data(idx_ovl_gg + 190);

    auto ts_yyzz_yyyz = pbuffer.data(idx_ovl_gg + 191);

    auto ts_yyzz_yyzz = pbuffer.data(idx_ovl_gg + 192);

    auto ts_yyzz_yzzz = pbuffer.data(idx_ovl_gg + 193);

    auto ts_yyzz_zzzz = pbuffer.data(idx_ovl_gg + 194);

    auto ts_yzzz_xxxx = pbuffer.data(idx_ovl_gg + 195);

    auto ts_yzzz_xxxz = pbuffer.data(idx_ovl_gg + 197);

    auto ts_yzzz_xxyz = pbuffer.data(idx_ovl_gg + 199);

    auto ts_yzzz_xxzz = pbuffer.data(idx_ovl_gg + 200);

    auto ts_yzzz_xyyz = pbuffer.data(idx_ovl_gg + 202);

    auto ts_yzzz_xyzz = pbuffer.data(idx_ovl_gg + 203);

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

    // Set up components of auxiliary buffer : HF

    auto ts_xxxxx_xxx = pbuffer.data(idx_ovl_hf);

    auto ts_xxxxx_xxy = pbuffer.data(idx_ovl_hf + 1);

    auto ts_xxxxx_xxz = pbuffer.data(idx_ovl_hf + 2);

    auto ts_xxxxx_xyy = pbuffer.data(idx_ovl_hf + 3);

    auto ts_xxxxx_xyz = pbuffer.data(idx_ovl_hf + 4);

    auto ts_xxxxx_xzz = pbuffer.data(idx_ovl_hf + 5);

    auto ts_xxxxx_yyy = pbuffer.data(idx_ovl_hf + 6);

    auto ts_xxxxx_yyz = pbuffer.data(idx_ovl_hf + 7);

    auto ts_xxxxx_yzz = pbuffer.data(idx_ovl_hf + 8);

    auto ts_xxxxx_zzz = pbuffer.data(idx_ovl_hf + 9);

    auto ts_xxxxz_xxz = pbuffer.data(idx_ovl_hf + 22);

    auto ts_xxxxz_xyz = pbuffer.data(idx_ovl_hf + 24);

    auto ts_xxxxz_xzz = pbuffer.data(idx_ovl_hf + 25);

    auto ts_xxxyy_xxy = pbuffer.data(idx_ovl_hf + 31);

    auto ts_xxxyy_xyy = pbuffer.data(idx_ovl_hf + 33);

    auto ts_xxxyy_xyz = pbuffer.data(idx_ovl_hf + 34);

    auto ts_xxxyy_yyy = pbuffer.data(idx_ovl_hf + 36);

    auto ts_xxxyy_yyz = pbuffer.data(idx_ovl_hf + 37);

    auto ts_xxxyy_yzz = pbuffer.data(idx_ovl_hf + 38);

    auto ts_xxxzz_xxx = pbuffer.data(idx_ovl_hf + 50);

    auto ts_xxxzz_xxy = pbuffer.data(idx_ovl_hf + 51);

    auto ts_xxxzz_xxz = pbuffer.data(idx_ovl_hf + 52);

    auto ts_xxxzz_xyy = pbuffer.data(idx_ovl_hf + 53);

    auto ts_xxxzz_xyz = pbuffer.data(idx_ovl_hf + 54);

    auto ts_xxxzz_xzz = pbuffer.data(idx_ovl_hf + 55);

    auto ts_xxxzz_yyz = pbuffer.data(idx_ovl_hf + 57);

    auto ts_xxxzz_yzz = pbuffer.data(idx_ovl_hf + 58);

    auto ts_xxxzz_zzz = pbuffer.data(idx_ovl_hf + 59);

    auto ts_xxyyy_xxy = pbuffer.data(idx_ovl_hf + 61);

    auto ts_xxyyy_xyy = pbuffer.data(idx_ovl_hf + 63);

    auto ts_xxyyy_xyz = pbuffer.data(idx_ovl_hf + 64);

    auto ts_xxyyy_yyy = pbuffer.data(idx_ovl_hf + 66);

    auto ts_xxyyy_yyz = pbuffer.data(idx_ovl_hf + 67);

    auto ts_xxyyy_yzz = pbuffer.data(idx_ovl_hf + 68);

    auto ts_xxzzz_xxx = pbuffer.data(idx_ovl_hf + 90);

    auto ts_xxzzz_xxy = pbuffer.data(idx_ovl_hf + 91);

    auto ts_xxzzz_xxz = pbuffer.data(idx_ovl_hf + 92);

    auto ts_xxzzz_xyy = pbuffer.data(idx_ovl_hf + 93);

    auto ts_xxzzz_xyz = pbuffer.data(idx_ovl_hf + 94);

    auto ts_xxzzz_xzz = pbuffer.data(idx_ovl_hf + 95);

    auto ts_xxzzz_yyz = pbuffer.data(idx_ovl_hf + 97);

    auto ts_xxzzz_yzz = pbuffer.data(idx_ovl_hf + 98);

    auto ts_xxzzz_zzz = pbuffer.data(idx_ovl_hf + 99);

    auto ts_xyyyy_xxy = pbuffer.data(idx_ovl_hf + 101);

    auto ts_xyyyy_xyy = pbuffer.data(idx_ovl_hf + 103);

    auto ts_xyyyy_xyz = pbuffer.data(idx_ovl_hf + 104);

    auto ts_xyyyy_yyy = pbuffer.data(idx_ovl_hf + 106);

    auto ts_xyyyy_yyz = pbuffer.data(idx_ovl_hf + 107);

    auto ts_xyyyy_yzz = pbuffer.data(idx_ovl_hf + 108);

    auto ts_xyyzz_xyz = pbuffer.data(idx_ovl_hf + 124);

    auto ts_xyyzz_yyz = pbuffer.data(idx_ovl_hf + 127);

    auto ts_xyyzz_yzz = pbuffer.data(idx_ovl_hf + 128);

    auto ts_xzzzz_xxz = pbuffer.data(idx_ovl_hf + 142);

    auto ts_xzzzz_xyz = pbuffer.data(idx_ovl_hf + 144);

    auto ts_xzzzz_xzz = pbuffer.data(idx_ovl_hf + 145);

    auto ts_xzzzz_yyz = pbuffer.data(idx_ovl_hf + 147);

    auto ts_xzzzz_yzz = pbuffer.data(idx_ovl_hf + 148);

    auto ts_xzzzz_zzz = pbuffer.data(idx_ovl_hf + 149);

    auto ts_yyyyy_xxx = pbuffer.data(idx_ovl_hf + 150);

    auto ts_yyyyy_xxy = pbuffer.data(idx_ovl_hf + 151);

    auto ts_yyyyy_xxz = pbuffer.data(idx_ovl_hf + 152);

    auto ts_yyyyy_xyy = pbuffer.data(idx_ovl_hf + 153);

    auto ts_yyyyy_xyz = pbuffer.data(idx_ovl_hf + 154);

    auto ts_yyyyy_xzz = pbuffer.data(idx_ovl_hf + 155);

    auto ts_yyyyy_yyy = pbuffer.data(idx_ovl_hf + 156);

    auto ts_yyyyy_yyz = pbuffer.data(idx_ovl_hf + 157);

    auto ts_yyyyy_yzz = pbuffer.data(idx_ovl_hf + 158);

    auto ts_yyyyy_zzz = pbuffer.data(idx_ovl_hf + 159);

    auto ts_yyyyz_xxz = pbuffer.data(idx_ovl_hf + 162);

    auto ts_yyyyz_xyz = pbuffer.data(idx_ovl_hf + 164);

    auto ts_yyyyz_xzz = pbuffer.data(idx_ovl_hf + 165);

    auto ts_yyyyz_yyz = pbuffer.data(idx_ovl_hf + 167);

    auto ts_yyyyz_yzz = pbuffer.data(idx_ovl_hf + 168);

    auto ts_yyyyz_zzz = pbuffer.data(idx_ovl_hf + 169);

    auto ts_yyyzz_xxx = pbuffer.data(idx_ovl_hf + 170);

    auto ts_yyyzz_xxy = pbuffer.data(idx_ovl_hf + 171);

    auto ts_yyyzz_xxz = pbuffer.data(idx_ovl_hf + 172);

    auto ts_yyyzz_xyy = pbuffer.data(idx_ovl_hf + 173);

    auto ts_yyyzz_xyz = pbuffer.data(idx_ovl_hf + 174);

    auto ts_yyyzz_xzz = pbuffer.data(idx_ovl_hf + 175);

    auto ts_yyyzz_yyy = pbuffer.data(idx_ovl_hf + 176);

    auto ts_yyyzz_yyz = pbuffer.data(idx_ovl_hf + 177);

    auto ts_yyyzz_yzz = pbuffer.data(idx_ovl_hf + 178);

    auto ts_yyyzz_zzz = pbuffer.data(idx_ovl_hf + 179);

    auto ts_yyzzz_xxx = pbuffer.data(idx_ovl_hf + 180);

    auto ts_yyzzz_xxy = pbuffer.data(idx_ovl_hf + 181);

    auto ts_yyzzz_xxz = pbuffer.data(idx_ovl_hf + 182);

    auto ts_yyzzz_xyy = pbuffer.data(idx_ovl_hf + 183);

    auto ts_yyzzz_xyz = pbuffer.data(idx_ovl_hf + 184);

    auto ts_yyzzz_xzz = pbuffer.data(idx_ovl_hf + 185);

    auto ts_yyzzz_yyy = pbuffer.data(idx_ovl_hf + 186);

    auto ts_yyzzz_yyz = pbuffer.data(idx_ovl_hf + 187);

    auto ts_yyzzz_yzz = pbuffer.data(idx_ovl_hf + 188);

    auto ts_yyzzz_zzz = pbuffer.data(idx_ovl_hf + 189);

    auto ts_yzzzz_xxy = pbuffer.data(idx_ovl_hf + 191);

    auto ts_yzzzz_xxz = pbuffer.data(idx_ovl_hf + 192);

    auto ts_yzzzz_xyy = pbuffer.data(idx_ovl_hf + 193);

    auto ts_yzzzz_xyz = pbuffer.data(idx_ovl_hf + 194);

    auto ts_yzzzz_xzz = pbuffer.data(idx_ovl_hf + 195);

    auto ts_yzzzz_yyy = pbuffer.data(idx_ovl_hf + 196);

    auto ts_yzzzz_yyz = pbuffer.data(idx_ovl_hf + 197);

    auto ts_yzzzz_yzz = pbuffer.data(idx_ovl_hf + 198);

    auto ts_yzzzz_zzz = pbuffer.data(idx_ovl_hf + 199);

    auto ts_zzzzz_xxx = pbuffer.data(idx_ovl_hf + 200);

    auto ts_zzzzz_xxy = pbuffer.data(idx_ovl_hf + 201);

    auto ts_zzzzz_xxz = pbuffer.data(idx_ovl_hf + 202);

    auto ts_zzzzz_xyy = pbuffer.data(idx_ovl_hf + 203);

    auto ts_zzzzz_xyz = pbuffer.data(idx_ovl_hf + 204);

    auto ts_zzzzz_xzz = pbuffer.data(idx_ovl_hf + 205);

    auto ts_zzzzz_yyy = pbuffer.data(idx_ovl_hf + 206);

    auto ts_zzzzz_yyz = pbuffer.data(idx_ovl_hf + 207);

    auto ts_zzzzz_yzz = pbuffer.data(idx_ovl_hf + 208);

    auto ts_zzzzz_zzz = pbuffer.data(idx_ovl_hf + 209);

    // Set up components of auxiliary buffer : HG

    auto ts_xxxxx_xxxx = pbuffer.data(idx_ovl_hg);

    auto ts_xxxxx_xxxy = pbuffer.data(idx_ovl_hg + 1);

    auto ts_xxxxx_xxxz = pbuffer.data(idx_ovl_hg + 2);

    auto ts_xxxxx_xxyy = pbuffer.data(idx_ovl_hg + 3);

    auto ts_xxxxx_xxyz = pbuffer.data(idx_ovl_hg + 4);

    auto ts_xxxxx_xxzz = pbuffer.data(idx_ovl_hg + 5);

    auto ts_xxxxx_xyyy = pbuffer.data(idx_ovl_hg + 6);

    auto ts_xxxxx_xyyz = pbuffer.data(idx_ovl_hg + 7);

    auto ts_xxxxx_xyzz = pbuffer.data(idx_ovl_hg + 8);

    auto ts_xxxxx_xzzz = pbuffer.data(idx_ovl_hg + 9);

    auto ts_xxxxx_yyyy = pbuffer.data(idx_ovl_hg + 10);

    auto ts_xxxxx_yyyz = pbuffer.data(idx_ovl_hg + 11);

    auto ts_xxxxx_yyzz = pbuffer.data(idx_ovl_hg + 12);

    auto ts_xxxxx_yzzz = pbuffer.data(idx_ovl_hg + 13);

    auto ts_xxxxx_zzzz = pbuffer.data(idx_ovl_hg + 14);

    auto ts_xxxxy_xxxx = pbuffer.data(idx_ovl_hg + 15);

    auto ts_xxxxy_xxxy = pbuffer.data(idx_ovl_hg + 16);

    auto ts_xxxxy_xxxz = pbuffer.data(idx_ovl_hg + 17);

    auto ts_xxxxy_xxyy = pbuffer.data(idx_ovl_hg + 18);

    auto ts_xxxxy_xxzz = pbuffer.data(idx_ovl_hg + 20);

    auto ts_xxxxy_xyyy = pbuffer.data(idx_ovl_hg + 21);

    auto ts_xxxxy_xzzz = pbuffer.data(idx_ovl_hg + 24);

    auto ts_xxxxy_yyyy = pbuffer.data(idx_ovl_hg + 25);

    auto ts_xxxxy_yyyz = pbuffer.data(idx_ovl_hg + 26);

    auto ts_xxxxy_yyzz = pbuffer.data(idx_ovl_hg + 27);

    auto ts_xxxxy_yzzz = pbuffer.data(idx_ovl_hg + 28);

    auto ts_xxxxz_xxxx = pbuffer.data(idx_ovl_hg + 30);

    auto ts_xxxxz_xxxy = pbuffer.data(idx_ovl_hg + 31);

    auto ts_xxxxz_xxxz = pbuffer.data(idx_ovl_hg + 32);

    auto ts_xxxxz_xxyy = pbuffer.data(idx_ovl_hg + 33);

    auto ts_xxxxz_xxyz = pbuffer.data(idx_ovl_hg + 34);

    auto ts_xxxxz_xxzz = pbuffer.data(idx_ovl_hg + 35);

    auto ts_xxxxz_xyyy = pbuffer.data(idx_ovl_hg + 36);

    auto ts_xxxxz_xyyz = pbuffer.data(idx_ovl_hg + 37);

    auto ts_xxxxz_xyzz = pbuffer.data(idx_ovl_hg + 38);

    auto ts_xxxxz_xzzz = pbuffer.data(idx_ovl_hg + 39);

    auto ts_xxxxz_yyyz = pbuffer.data(idx_ovl_hg + 41);

    auto ts_xxxxz_yyzz = pbuffer.data(idx_ovl_hg + 42);

    auto ts_xxxxz_yzzz = pbuffer.data(idx_ovl_hg + 43);

    auto ts_xxxxz_zzzz = pbuffer.data(idx_ovl_hg + 44);

    auto ts_xxxyy_xxxx = pbuffer.data(idx_ovl_hg + 45);

    auto ts_xxxyy_xxxy = pbuffer.data(idx_ovl_hg + 46);

    auto ts_xxxyy_xxxz = pbuffer.data(idx_ovl_hg + 47);

    auto ts_xxxyy_xxyy = pbuffer.data(idx_ovl_hg + 48);

    auto ts_xxxyy_xxyz = pbuffer.data(idx_ovl_hg + 49);

    auto ts_xxxyy_xxzz = pbuffer.data(idx_ovl_hg + 50);

    auto ts_xxxyy_xyyy = pbuffer.data(idx_ovl_hg + 51);

    auto ts_xxxyy_xyyz = pbuffer.data(idx_ovl_hg + 52);

    auto ts_xxxyy_xyzz = pbuffer.data(idx_ovl_hg + 53);

    auto ts_xxxyy_xzzz = pbuffer.data(idx_ovl_hg + 54);

    auto ts_xxxyy_yyyy = pbuffer.data(idx_ovl_hg + 55);

    auto ts_xxxyy_yyyz = pbuffer.data(idx_ovl_hg + 56);

    auto ts_xxxyy_yyzz = pbuffer.data(idx_ovl_hg + 57);

    auto ts_xxxyy_yzzz = pbuffer.data(idx_ovl_hg + 58);

    auto ts_xxxyy_zzzz = pbuffer.data(idx_ovl_hg + 59);

    auto ts_xxxyz_xxxz = pbuffer.data(idx_ovl_hg + 62);

    auto ts_xxxyz_xxzz = pbuffer.data(idx_ovl_hg + 65);

    auto ts_xxxyz_xzzz = pbuffer.data(idx_ovl_hg + 69);

    auto ts_xxxyz_yyyz = pbuffer.data(idx_ovl_hg + 71);

    auto ts_xxxyz_yyzz = pbuffer.data(idx_ovl_hg + 72);

    auto ts_xxxyz_yzzz = pbuffer.data(idx_ovl_hg + 73);

    auto ts_xxxzz_xxxx = pbuffer.data(idx_ovl_hg + 75);

    auto ts_xxxzz_xxxy = pbuffer.data(idx_ovl_hg + 76);

    auto ts_xxxzz_xxxz = pbuffer.data(idx_ovl_hg + 77);

    auto ts_xxxzz_xxyy = pbuffer.data(idx_ovl_hg + 78);

    auto ts_xxxzz_xxyz = pbuffer.data(idx_ovl_hg + 79);

    auto ts_xxxzz_xxzz = pbuffer.data(idx_ovl_hg + 80);

    auto ts_xxxzz_xyyy = pbuffer.data(idx_ovl_hg + 81);

    auto ts_xxxzz_xyyz = pbuffer.data(idx_ovl_hg + 82);

    auto ts_xxxzz_xyzz = pbuffer.data(idx_ovl_hg + 83);

    auto ts_xxxzz_xzzz = pbuffer.data(idx_ovl_hg + 84);

    auto ts_xxxzz_yyyy = pbuffer.data(idx_ovl_hg + 85);

    auto ts_xxxzz_yyyz = pbuffer.data(idx_ovl_hg + 86);

    auto ts_xxxzz_yyzz = pbuffer.data(idx_ovl_hg + 87);

    auto ts_xxxzz_yzzz = pbuffer.data(idx_ovl_hg + 88);

    auto ts_xxxzz_zzzz = pbuffer.data(idx_ovl_hg + 89);

    auto ts_xxyyy_xxxx = pbuffer.data(idx_ovl_hg + 90);

    auto ts_xxyyy_xxxy = pbuffer.data(idx_ovl_hg + 91);

    auto ts_xxyyy_xxxz = pbuffer.data(idx_ovl_hg + 92);

    auto ts_xxyyy_xxyy = pbuffer.data(idx_ovl_hg + 93);

    auto ts_xxyyy_xxyz = pbuffer.data(idx_ovl_hg + 94);

    auto ts_xxyyy_xxzz = pbuffer.data(idx_ovl_hg + 95);

    auto ts_xxyyy_xyyy = pbuffer.data(idx_ovl_hg + 96);

    auto ts_xxyyy_xyyz = pbuffer.data(idx_ovl_hg + 97);

    auto ts_xxyyy_xyzz = pbuffer.data(idx_ovl_hg + 98);

    auto ts_xxyyy_xzzz = pbuffer.data(idx_ovl_hg + 99);

    auto ts_xxyyy_yyyy = pbuffer.data(idx_ovl_hg + 100);

    auto ts_xxyyy_yyyz = pbuffer.data(idx_ovl_hg + 101);

    auto ts_xxyyy_yyzz = pbuffer.data(idx_ovl_hg + 102);

    auto ts_xxyyy_yzzz = pbuffer.data(idx_ovl_hg + 103);

    auto ts_xxyyy_zzzz = pbuffer.data(idx_ovl_hg + 104);

    auto ts_xxyyz_xxxy = pbuffer.data(idx_ovl_hg + 106);

    auto ts_xxyyz_xxxz = pbuffer.data(idx_ovl_hg + 107);

    auto ts_xxyyz_xxyy = pbuffer.data(idx_ovl_hg + 108);

    auto ts_xxyyz_xxzz = pbuffer.data(idx_ovl_hg + 110);

    auto ts_xxyyz_xyyy = pbuffer.data(idx_ovl_hg + 111);

    auto ts_xxyyz_xzzz = pbuffer.data(idx_ovl_hg + 114);

    auto ts_xxyyz_yyyz = pbuffer.data(idx_ovl_hg + 116);

    auto ts_xxyyz_yyzz = pbuffer.data(idx_ovl_hg + 117);

    auto ts_xxyyz_yzzz = pbuffer.data(idx_ovl_hg + 118);

    auto ts_xxyyz_zzzz = pbuffer.data(idx_ovl_hg + 119);

    auto ts_xxyzz_xxxx = pbuffer.data(idx_ovl_hg + 120);

    auto ts_xxyzz_xxxz = pbuffer.data(idx_ovl_hg + 122);

    auto ts_xxyzz_xxzz = pbuffer.data(idx_ovl_hg + 125);

    auto ts_xxyzz_xzzz = pbuffer.data(idx_ovl_hg + 129);

    auto ts_xxyzz_yyyy = pbuffer.data(idx_ovl_hg + 130);

    auto ts_xxyzz_yyyz = pbuffer.data(idx_ovl_hg + 131);

    auto ts_xxyzz_yyzz = pbuffer.data(idx_ovl_hg + 132);

    auto ts_xxyzz_yzzz = pbuffer.data(idx_ovl_hg + 133);

    auto ts_xxzzz_xxxx = pbuffer.data(idx_ovl_hg + 135);

    auto ts_xxzzz_xxxy = pbuffer.data(idx_ovl_hg + 136);

    auto ts_xxzzz_xxxz = pbuffer.data(idx_ovl_hg + 137);

    auto ts_xxzzz_xxyy = pbuffer.data(idx_ovl_hg + 138);

    auto ts_xxzzz_xxyz = pbuffer.data(idx_ovl_hg + 139);

    auto ts_xxzzz_xxzz = pbuffer.data(idx_ovl_hg + 140);

    auto ts_xxzzz_xyyy = pbuffer.data(idx_ovl_hg + 141);

    auto ts_xxzzz_xyyz = pbuffer.data(idx_ovl_hg + 142);

    auto ts_xxzzz_xyzz = pbuffer.data(idx_ovl_hg + 143);

    auto ts_xxzzz_xzzz = pbuffer.data(idx_ovl_hg + 144);

    auto ts_xxzzz_yyyy = pbuffer.data(idx_ovl_hg + 145);

    auto ts_xxzzz_yyyz = pbuffer.data(idx_ovl_hg + 146);

    auto ts_xxzzz_yyzz = pbuffer.data(idx_ovl_hg + 147);

    auto ts_xxzzz_yzzz = pbuffer.data(idx_ovl_hg + 148);

    auto ts_xxzzz_zzzz = pbuffer.data(idx_ovl_hg + 149);

    auto ts_xyyyy_xxxx = pbuffer.data(idx_ovl_hg + 150);

    auto ts_xyyyy_xxxy = pbuffer.data(idx_ovl_hg + 151);

    auto ts_xyyyy_xxyy = pbuffer.data(idx_ovl_hg + 153);

    auto ts_xyyyy_xxyz = pbuffer.data(idx_ovl_hg + 154);

    auto ts_xyyyy_xyyy = pbuffer.data(idx_ovl_hg + 156);

    auto ts_xyyyy_xyyz = pbuffer.data(idx_ovl_hg + 157);

    auto ts_xyyyy_xyzz = pbuffer.data(idx_ovl_hg + 158);

    auto ts_xyyyy_yyyy = pbuffer.data(idx_ovl_hg + 160);

    auto ts_xyyyy_yyyz = pbuffer.data(idx_ovl_hg + 161);

    auto ts_xyyyy_yyzz = pbuffer.data(idx_ovl_hg + 162);

    auto ts_xyyyy_yzzz = pbuffer.data(idx_ovl_hg + 163);

    auto ts_xyyyy_zzzz = pbuffer.data(idx_ovl_hg + 164);

    auto ts_xyyyz_yyyz = pbuffer.data(idx_ovl_hg + 176);

    auto ts_xyyyz_yyzz = pbuffer.data(idx_ovl_hg + 177);

    auto ts_xyyyz_yzzz = pbuffer.data(idx_ovl_hg + 178);

    auto ts_xyyyz_zzzz = pbuffer.data(idx_ovl_hg + 179);

    auto ts_xyyzz_xxyz = pbuffer.data(idx_ovl_hg + 184);

    auto ts_xyyzz_xyyz = pbuffer.data(idx_ovl_hg + 187);

    auto ts_xyyzz_xyzz = pbuffer.data(idx_ovl_hg + 188);

    auto ts_xyyzz_yyyy = pbuffer.data(idx_ovl_hg + 190);

    auto ts_xyyzz_yyyz = pbuffer.data(idx_ovl_hg + 191);

    auto ts_xyyzz_yyzz = pbuffer.data(idx_ovl_hg + 192);

    auto ts_xyyzz_yzzz = pbuffer.data(idx_ovl_hg + 193);

    auto ts_xyyzz_zzzz = pbuffer.data(idx_ovl_hg + 194);

    auto ts_xyzzz_yyyy = pbuffer.data(idx_ovl_hg + 205);

    auto ts_xyzzz_yyyz = pbuffer.data(idx_ovl_hg + 206);

    auto ts_xyzzz_yyzz = pbuffer.data(idx_ovl_hg + 207);

    auto ts_xyzzz_yzzz = pbuffer.data(idx_ovl_hg + 208);

    auto ts_xzzzz_xxxx = pbuffer.data(idx_ovl_hg + 210);

    auto ts_xzzzz_xxxz = pbuffer.data(idx_ovl_hg + 212);

    auto ts_xzzzz_xxyz = pbuffer.data(idx_ovl_hg + 214);

    auto ts_xzzzz_xxzz = pbuffer.data(idx_ovl_hg + 215);

    auto ts_xzzzz_xyyz = pbuffer.data(idx_ovl_hg + 217);

    auto ts_xzzzz_xyzz = pbuffer.data(idx_ovl_hg + 218);

    auto ts_xzzzz_xzzz = pbuffer.data(idx_ovl_hg + 219);

    auto ts_xzzzz_yyyy = pbuffer.data(idx_ovl_hg + 220);

    auto ts_xzzzz_yyyz = pbuffer.data(idx_ovl_hg + 221);

    auto ts_xzzzz_yyzz = pbuffer.data(idx_ovl_hg + 222);

    auto ts_xzzzz_yzzz = pbuffer.data(idx_ovl_hg + 223);

    auto ts_xzzzz_zzzz = pbuffer.data(idx_ovl_hg + 224);

    auto ts_yyyyy_xxxx = pbuffer.data(idx_ovl_hg + 225);

    auto ts_yyyyy_xxxy = pbuffer.data(idx_ovl_hg + 226);

    auto ts_yyyyy_xxxz = pbuffer.data(idx_ovl_hg + 227);

    auto ts_yyyyy_xxyy = pbuffer.data(idx_ovl_hg + 228);

    auto ts_yyyyy_xxyz = pbuffer.data(idx_ovl_hg + 229);

    auto ts_yyyyy_xxzz = pbuffer.data(idx_ovl_hg + 230);

    auto ts_yyyyy_xyyy = pbuffer.data(idx_ovl_hg + 231);

    auto ts_yyyyy_xyyz = pbuffer.data(idx_ovl_hg + 232);

    auto ts_yyyyy_xyzz = pbuffer.data(idx_ovl_hg + 233);

    auto ts_yyyyy_xzzz = pbuffer.data(idx_ovl_hg + 234);

    auto ts_yyyyy_yyyy = pbuffer.data(idx_ovl_hg + 235);

    auto ts_yyyyy_yyyz = pbuffer.data(idx_ovl_hg + 236);

    auto ts_yyyyy_yyzz = pbuffer.data(idx_ovl_hg + 237);

    auto ts_yyyyy_yzzz = pbuffer.data(idx_ovl_hg + 238);

    auto ts_yyyyy_zzzz = pbuffer.data(idx_ovl_hg + 239);

    auto ts_yyyyz_xxxy = pbuffer.data(idx_ovl_hg + 241);

    auto ts_yyyyz_xxxz = pbuffer.data(idx_ovl_hg + 242);

    auto ts_yyyyz_xxyy = pbuffer.data(idx_ovl_hg + 243);

    auto ts_yyyyz_xxyz = pbuffer.data(idx_ovl_hg + 244);

    auto ts_yyyyz_xxzz = pbuffer.data(idx_ovl_hg + 245);

    auto ts_yyyyz_xyyy = pbuffer.data(idx_ovl_hg + 246);

    auto ts_yyyyz_xyyz = pbuffer.data(idx_ovl_hg + 247);

    auto ts_yyyyz_xyzz = pbuffer.data(idx_ovl_hg + 248);

    auto ts_yyyyz_xzzz = pbuffer.data(idx_ovl_hg + 249);

    auto ts_yyyyz_yyyy = pbuffer.data(idx_ovl_hg + 250);

    auto ts_yyyyz_yyyz = pbuffer.data(idx_ovl_hg + 251);

    auto ts_yyyyz_yyzz = pbuffer.data(idx_ovl_hg + 252);

    auto ts_yyyyz_yzzz = pbuffer.data(idx_ovl_hg + 253);

    auto ts_yyyyz_zzzz = pbuffer.data(idx_ovl_hg + 254);

    auto ts_yyyzz_xxxx = pbuffer.data(idx_ovl_hg + 255);

    auto ts_yyyzz_xxxy = pbuffer.data(idx_ovl_hg + 256);

    auto ts_yyyzz_xxxz = pbuffer.data(idx_ovl_hg + 257);

    auto ts_yyyzz_xxyy = pbuffer.data(idx_ovl_hg + 258);

    auto ts_yyyzz_xxyz = pbuffer.data(idx_ovl_hg + 259);

    auto ts_yyyzz_xxzz = pbuffer.data(idx_ovl_hg + 260);

    auto ts_yyyzz_xyyy = pbuffer.data(idx_ovl_hg + 261);

    auto ts_yyyzz_xyyz = pbuffer.data(idx_ovl_hg + 262);

    auto ts_yyyzz_xyzz = pbuffer.data(idx_ovl_hg + 263);

    auto ts_yyyzz_xzzz = pbuffer.data(idx_ovl_hg + 264);

    auto ts_yyyzz_yyyy = pbuffer.data(idx_ovl_hg + 265);

    auto ts_yyyzz_yyyz = pbuffer.data(idx_ovl_hg + 266);

    auto ts_yyyzz_yyzz = pbuffer.data(idx_ovl_hg + 267);

    auto ts_yyyzz_yzzz = pbuffer.data(idx_ovl_hg + 268);

    auto ts_yyyzz_zzzz = pbuffer.data(idx_ovl_hg + 269);

    auto ts_yyzzz_xxxx = pbuffer.data(idx_ovl_hg + 270);

    auto ts_yyzzz_xxxy = pbuffer.data(idx_ovl_hg + 271);

    auto ts_yyzzz_xxxz = pbuffer.data(idx_ovl_hg + 272);

    auto ts_yyzzz_xxyy = pbuffer.data(idx_ovl_hg + 273);

    auto ts_yyzzz_xxyz = pbuffer.data(idx_ovl_hg + 274);

    auto ts_yyzzz_xxzz = pbuffer.data(idx_ovl_hg + 275);

    auto ts_yyzzz_xyyy = pbuffer.data(idx_ovl_hg + 276);

    auto ts_yyzzz_xyyz = pbuffer.data(idx_ovl_hg + 277);

    auto ts_yyzzz_xyzz = pbuffer.data(idx_ovl_hg + 278);

    auto ts_yyzzz_xzzz = pbuffer.data(idx_ovl_hg + 279);

    auto ts_yyzzz_yyyy = pbuffer.data(idx_ovl_hg + 280);

    auto ts_yyzzz_yyyz = pbuffer.data(idx_ovl_hg + 281);

    auto ts_yyzzz_yyzz = pbuffer.data(idx_ovl_hg + 282);

    auto ts_yyzzz_yzzz = pbuffer.data(idx_ovl_hg + 283);

    auto ts_yyzzz_zzzz = pbuffer.data(idx_ovl_hg + 284);

    auto ts_yzzzz_xxxx = pbuffer.data(idx_ovl_hg + 285);

    auto ts_yzzzz_xxxy = pbuffer.data(idx_ovl_hg + 286);

    auto ts_yzzzz_xxxz = pbuffer.data(idx_ovl_hg + 287);

    auto ts_yzzzz_xxyy = pbuffer.data(idx_ovl_hg + 288);

    auto ts_yzzzz_xxyz = pbuffer.data(idx_ovl_hg + 289);

    auto ts_yzzzz_xxzz = pbuffer.data(idx_ovl_hg + 290);

    auto ts_yzzzz_xyyy = pbuffer.data(idx_ovl_hg + 291);

    auto ts_yzzzz_xyyz = pbuffer.data(idx_ovl_hg + 292);

    auto ts_yzzzz_xyzz = pbuffer.data(idx_ovl_hg + 293);

    auto ts_yzzzz_xzzz = pbuffer.data(idx_ovl_hg + 294);

    auto ts_yzzzz_yyyy = pbuffer.data(idx_ovl_hg + 295);

    auto ts_yzzzz_yyyz = pbuffer.data(idx_ovl_hg + 296);

    auto ts_yzzzz_yyzz = pbuffer.data(idx_ovl_hg + 297);

    auto ts_yzzzz_yzzz = pbuffer.data(idx_ovl_hg + 298);

    auto ts_yzzzz_zzzz = pbuffer.data(idx_ovl_hg + 299);

    auto ts_zzzzz_xxxx = pbuffer.data(idx_ovl_hg + 300);

    auto ts_zzzzz_xxxy = pbuffer.data(idx_ovl_hg + 301);

    auto ts_zzzzz_xxxz = pbuffer.data(idx_ovl_hg + 302);

    auto ts_zzzzz_xxyy = pbuffer.data(idx_ovl_hg + 303);

    auto ts_zzzzz_xxyz = pbuffer.data(idx_ovl_hg + 304);

    auto ts_zzzzz_xxzz = pbuffer.data(idx_ovl_hg + 305);

    auto ts_zzzzz_xyyy = pbuffer.data(idx_ovl_hg + 306);

    auto ts_zzzzz_xyyz = pbuffer.data(idx_ovl_hg + 307);

    auto ts_zzzzz_xyzz = pbuffer.data(idx_ovl_hg + 308);

    auto ts_zzzzz_xzzz = pbuffer.data(idx_ovl_hg + 309);

    auto ts_zzzzz_yyyy = pbuffer.data(idx_ovl_hg + 310);

    auto ts_zzzzz_yyyz = pbuffer.data(idx_ovl_hg + 311);

    auto ts_zzzzz_yyzz = pbuffer.data(idx_ovl_hg + 312);

    auto ts_zzzzz_yzzz = pbuffer.data(idx_ovl_hg + 313);

    auto ts_zzzzz_zzzz = pbuffer.data(idx_ovl_hg + 314);

    // Set up 0-15 components of targeted buffer : IG

    auto ts_xxxxxx_xxxx = pbuffer.data(idx_ovl_ig);

    auto ts_xxxxxx_xxxy = pbuffer.data(idx_ovl_ig + 1);

    auto ts_xxxxxx_xxxz = pbuffer.data(idx_ovl_ig + 2);

    auto ts_xxxxxx_xxyy = pbuffer.data(idx_ovl_ig + 3);

    auto ts_xxxxxx_xxyz = pbuffer.data(idx_ovl_ig + 4);

    auto ts_xxxxxx_xxzz = pbuffer.data(idx_ovl_ig + 5);

    auto ts_xxxxxx_xyyy = pbuffer.data(idx_ovl_ig + 6);

    auto ts_xxxxxx_xyyz = pbuffer.data(idx_ovl_ig + 7);

    auto ts_xxxxxx_xyzz = pbuffer.data(idx_ovl_ig + 8);

    auto ts_xxxxxx_xzzz = pbuffer.data(idx_ovl_ig + 9);

    auto ts_xxxxxx_yyyy = pbuffer.data(idx_ovl_ig + 10);

    auto ts_xxxxxx_yyyz = pbuffer.data(idx_ovl_ig + 11);

    auto ts_xxxxxx_yyzz = pbuffer.data(idx_ovl_ig + 12);

    auto ts_xxxxxx_yzzz = pbuffer.data(idx_ovl_ig + 13);

    auto ts_xxxxxx_zzzz = pbuffer.data(idx_ovl_ig + 14);

    #pragma omp simd aligned(pa_x, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxxz, ts_xxxx_xxyy, ts_xxxx_xxyz, ts_xxxx_xxzz, ts_xxxx_xyyy, ts_xxxx_xyyz, ts_xxxx_xyzz, ts_xxxx_xzzz, ts_xxxx_yyyy, ts_xxxx_yyyz, ts_xxxx_yyzz, ts_xxxx_yzzz, ts_xxxx_zzzz, ts_xxxxx_xxx, ts_xxxxx_xxxx, ts_xxxxx_xxxy, ts_xxxxx_xxxz, ts_xxxxx_xxy, ts_xxxxx_xxyy, ts_xxxxx_xxyz, ts_xxxxx_xxz, ts_xxxxx_xxzz, ts_xxxxx_xyy, ts_xxxxx_xyyy, ts_xxxxx_xyyz, ts_xxxxx_xyz, ts_xxxxx_xyzz, ts_xxxxx_xzz, ts_xxxxx_xzzz, ts_xxxxx_yyy, ts_xxxxx_yyyy, ts_xxxxx_yyyz, ts_xxxxx_yyz, ts_xxxxx_yyzz, ts_xxxxx_yzz, ts_xxxxx_yzzz, ts_xxxxx_zzz, ts_xxxxx_zzzz, ts_xxxxxx_xxxx, ts_xxxxxx_xxxy, ts_xxxxxx_xxxz, ts_xxxxxx_xxyy, ts_xxxxxx_xxyz, ts_xxxxxx_xxzz, ts_xxxxxx_xyyy, ts_xxxxxx_xyyz, ts_xxxxxx_xyzz, ts_xxxxxx_xzzz, ts_xxxxxx_yyyy, ts_xxxxxx_yyyz, ts_xxxxxx_yyzz, ts_xxxxxx_yzzz, ts_xxxxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxx_xxxx[i] = 5.0 * ts_xxxx_xxxx[i] * fe_0 + 4.0 * ts_xxxxx_xxx[i] * fe_0 + ts_xxxxx_xxxx[i] * pa_x[i];

        ts_xxxxxx_xxxy[i] = 5.0 * ts_xxxx_xxxy[i] * fe_0 + 3.0 * ts_xxxxx_xxy[i] * fe_0 + ts_xxxxx_xxxy[i] * pa_x[i];

        ts_xxxxxx_xxxz[i] = 5.0 * ts_xxxx_xxxz[i] * fe_0 + 3.0 * ts_xxxxx_xxz[i] * fe_0 + ts_xxxxx_xxxz[i] * pa_x[i];

        ts_xxxxxx_xxyy[i] = 5.0 * ts_xxxx_xxyy[i] * fe_0 + 2.0 * ts_xxxxx_xyy[i] * fe_0 + ts_xxxxx_xxyy[i] * pa_x[i];

        ts_xxxxxx_xxyz[i] = 5.0 * ts_xxxx_xxyz[i] * fe_0 + 2.0 * ts_xxxxx_xyz[i] * fe_0 + ts_xxxxx_xxyz[i] * pa_x[i];

        ts_xxxxxx_xxzz[i] = 5.0 * ts_xxxx_xxzz[i] * fe_0 + 2.0 * ts_xxxxx_xzz[i] * fe_0 + ts_xxxxx_xxzz[i] * pa_x[i];

        ts_xxxxxx_xyyy[i] = 5.0 * ts_xxxx_xyyy[i] * fe_0 + ts_xxxxx_yyy[i] * fe_0 + ts_xxxxx_xyyy[i] * pa_x[i];

        ts_xxxxxx_xyyz[i] = 5.0 * ts_xxxx_xyyz[i] * fe_0 + ts_xxxxx_yyz[i] * fe_0 + ts_xxxxx_xyyz[i] * pa_x[i];

        ts_xxxxxx_xyzz[i] = 5.0 * ts_xxxx_xyzz[i] * fe_0 + ts_xxxxx_yzz[i] * fe_0 + ts_xxxxx_xyzz[i] * pa_x[i];

        ts_xxxxxx_xzzz[i] = 5.0 * ts_xxxx_xzzz[i] * fe_0 + ts_xxxxx_zzz[i] * fe_0 + ts_xxxxx_xzzz[i] * pa_x[i];

        ts_xxxxxx_yyyy[i] = 5.0 * ts_xxxx_yyyy[i] * fe_0 + ts_xxxxx_yyyy[i] * pa_x[i];

        ts_xxxxxx_yyyz[i] = 5.0 * ts_xxxx_yyyz[i] * fe_0 + ts_xxxxx_yyyz[i] * pa_x[i];

        ts_xxxxxx_yyzz[i] = 5.0 * ts_xxxx_yyzz[i] * fe_0 + ts_xxxxx_yyzz[i] * pa_x[i];

        ts_xxxxxx_yzzz[i] = 5.0 * ts_xxxx_yzzz[i] * fe_0 + ts_xxxxx_yzzz[i] * pa_x[i];

        ts_xxxxxx_zzzz[i] = 5.0 * ts_xxxx_zzzz[i] * fe_0 + ts_xxxxx_zzzz[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : IG

    auto ts_xxxxxy_xxxx = pbuffer.data(idx_ovl_ig + 15);

    auto ts_xxxxxy_xxxy = pbuffer.data(idx_ovl_ig + 16);

    auto ts_xxxxxy_xxxz = pbuffer.data(idx_ovl_ig + 17);

    auto ts_xxxxxy_xxyy = pbuffer.data(idx_ovl_ig + 18);

    auto ts_xxxxxy_xxyz = pbuffer.data(idx_ovl_ig + 19);

    auto ts_xxxxxy_xxzz = pbuffer.data(idx_ovl_ig + 20);

    auto ts_xxxxxy_xyyy = pbuffer.data(idx_ovl_ig + 21);

    auto ts_xxxxxy_xyyz = pbuffer.data(idx_ovl_ig + 22);

    auto ts_xxxxxy_xyzz = pbuffer.data(idx_ovl_ig + 23);

    auto ts_xxxxxy_xzzz = pbuffer.data(idx_ovl_ig + 24);

    auto ts_xxxxxy_yyyy = pbuffer.data(idx_ovl_ig + 25);

    auto ts_xxxxxy_yyyz = pbuffer.data(idx_ovl_ig + 26);

    auto ts_xxxxxy_yyzz = pbuffer.data(idx_ovl_ig + 27);

    auto ts_xxxxxy_yzzz = pbuffer.data(idx_ovl_ig + 28);

    auto ts_xxxxxy_zzzz = pbuffer.data(idx_ovl_ig + 29);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxxx_xxx, ts_xxxxx_xxxx, ts_xxxxx_xxxy, ts_xxxxx_xxxz, ts_xxxxx_xxy, ts_xxxxx_xxyy, ts_xxxxx_xxyz, ts_xxxxx_xxz, ts_xxxxx_xxzz, ts_xxxxx_xyy, ts_xxxxx_xyyy, ts_xxxxx_xyyz, ts_xxxxx_xyz, ts_xxxxx_xyzz, ts_xxxxx_xzz, ts_xxxxx_xzzz, ts_xxxxx_zzzz, ts_xxxxxy_xxxx, ts_xxxxxy_xxxy, ts_xxxxxy_xxxz, ts_xxxxxy_xxyy, ts_xxxxxy_xxyz, ts_xxxxxy_xxzz, ts_xxxxxy_xyyy, ts_xxxxxy_xyyz, ts_xxxxxy_xyzz, ts_xxxxxy_xzzz, ts_xxxxxy_yyyy, ts_xxxxxy_yyyz, ts_xxxxxy_yyzz, ts_xxxxxy_yzzz, ts_xxxxxy_zzzz, ts_xxxxy_yyyy, ts_xxxxy_yyyz, ts_xxxxy_yyzz, ts_xxxxy_yzzz, ts_xxxy_yyyy, ts_xxxy_yyyz, ts_xxxy_yyzz, ts_xxxy_yzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxy_xxxx[i] = ts_xxxxx_xxxx[i] * pa_y[i];

        ts_xxxxxy_xxxy[i] = ts_xxxxx_xxx[i] * fe_0 + ts_xxxxx_xxxy[i] * pa_y[i];

        ts_xxxxxy_xxxz[i] = ts_xxxxx_xxxz[i] * pa_y[i];

        ts_xxxxxy_xxyy[i] = 2.0 * ts_xxxxx_xxy[i] * fe_0 + ts_xxxxx_xxyy[i] * pa_y[i];

        ts_xxxxxy_xxyz[i] = ts_xxxxx_xxz[i] * fe_0 + ts_xxxxx_xxyz[i] * pa_y[i];

        ts_xxxxxy_xxzz[i] = ts_xxxxx_xxzz[i] * pa_y[i];

        ts_xxxxxy_xyyy[i] = 3.0 * ts_xxxxx_xyy[i] * fe_0 + ts_xxxxx_xyyy[i] * pa_y[i];

        ts_xxxxxy_xyyz[i] = 2.0 * ts_xxxxx_xyz[i] * fe_0 + ts_xxxxx_xyyz[i] * pa_y[i];

        ts_xxxxxy_xyzz[i] = ts_xxxxx_xzz[i] * fe_0 + ts_xxxxx_xyzz[i] * pa_y[i];

        ts_xxxxxy_xzzz[i] = ts_xxxxx_xzzz[i] * pa_y[i];

        ts_xxxxxy_yyyy[i] = 4.0 * ts_xxxy_yyyy[i] * fe_0 + ts_xxxxy_yyyy[i] * pa_x[i];

        ts_xxxxxy_yyyz[i] = 4.0 * ts_xxxy_yyyz[i] * fe_0 + ts_xxxxy_yyyz[i] * pa_x[i];

        ts_xxxxxy_yyzz[i] = 4.0 * ts_xxxy_yyzz[i] * fe_0 + ts_xxxxy_yyzz[i] * pa_x[i];

        ts_xxxxxy_yzzz[i] = 4.0 * ts_xxxy_yzzz[i] * fe_0 + ts_xxxxy_yzzz[i] * pa_x[i];

        ts_xxxxxy_zzzz[i] = ts_xxxxx_zzzz[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : IG

    auto ts_xxxxxz_xxxx = pbuffer.data(idx_ovl_ig + 30);

    auto ts_xxxxxz_xxxy = pbuffer.data(idx_ovl_ig + 31);

    auto ts_xxxxxz_xxxz = pbuffer.data(idx_ovl_ig + 32);

    auto ts_xxxxxz_xxyy = pbuffer.data(idx_ovl_ig + 33);

    auto ts_xxxxxz_xxyz = pbuffer.data(idx_ovl_ig + 34);

    auto ts_xxxxxz_xxzz = pbuffer.data(idx_ovl_ig + 35);

    auto ts_xxxxxz_xyyy = pbuffer.data(idx_ovl_ig + 36);

    auto ts_xxxxxz_xyyz = pbuffer.data(idx_ovl_ig + 37);

    auto ts_xxxxxz_xyzz = pbuffer.data(idx_ovl_ig + 38);

    auto ts_xxxxxz_xzzz = pbuffer.data(idx_ovl_ig + 39);

    auto ts_xxxxxz_yyyy = pbuffer.data(idx_ovl_ig + 40);

    auto ts_xxxxxz_yyyz = pbuffer.data(idx_ovl_ig + 41);

    auto ts_xxxxxz_yyzz = pbuffer.data(idx_ovl_ig + 42);

    auto ts_xxxxxz_yzzz = pbuffer.data(idx_ovl_ig + 43);

    auto ts_xxxxxz_zzzz = pbuffer.data(idx_ovl_ig + 44);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxxx_xxx, ts_xxxxx_xxxx, ts_xxxxx_xxxy, ts_xxxxx_xxxz, ts_xxxxx_xxy, ts_xxxxx_xxyy, ts_xxxxx_xxyz, ts_xxxxx_xxz, ts_xxxxx_xxzz, ts_xxxxx_xyy, ts_xxxxx_xyyy, ts_xxxxx_xyyz, ts_xxxxx_xyz, ts_xxxxx_xyzz, ts_xxxxx_xzz, ts_xxxxx_xzzz, ts_xxxxx_yyyy, ts_xxxxxz_xxxx, ts_xxxxxz_xxxy, ts_xxxxxz_xxxz, ts_xxxxxz_xxyy, ts_xxxxxz_xxyz, ts_xxxxxz_xxzz, ts_xxxxxz_xyyy, ts_xxxxxz_xyyz, ts_xxxxxz_xyzz, ts_xxxxxz_xzzz, ts_xxxxxz_yyyy, ts_xxxxxz_yyyz, ts_xxxxxz_yyzz, ts_xxxxxz_yzzz, ts_xxxxxz_zzzz, ts_xxxxz_yyyz, ts_xxxxz_yyzz, ts_xxxxz_yzzz, ts_xxxxz_zzzz, ts_xxxz_yyyz, ts_xxxz_yyzz, ts_xxxz_yzzz, ts_xxxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxz_xxxx[i] = ts_xxxxx_xxxx[i] * pa_z[i];

        ts_xxxxxz_xxxy[i] = ts_xxxxx_xxxy[i] * pa_z[i];

        ts_xxxxxz_xxxz[i] = ts_xxxxx_xxx[i] * fe_0 + ts_xxxxx_xxxz[i] * pa_z[i];

        ts_xxxxxz_xxyy[i] = ts_xxxxx_xxyy[i] * pa_z[i];

        ts_xxxxxz_xxyz[i] = ts_xxxxx_xxy[i] * fe_0 + ts_xxxxx_xxyz[i] * pa_z[i];

        ts_xxxxxz_xxzz[i] = 2.0 * ts_xxxxx_xxz[i] * fe_0 + ts_xxxxx_xxzz[i] * pa_z[i];

        ts_xxxxxz_xyyy[i] = ts_xxxxx_xyyy[i] * pa_z[i];

        ts_xxxxxz_xyyz[i] = ts_xxxxx_xyy[i] * fe_0 + ts_xxxxx_xyyz[i] * pa_z[i];

        ts_xxxxxz_xyzz[i] = 2.0 * ts_xxxxx_xyz[i] * fe_0 + ts_xxxxx_xyzz[i] * pa_z[i];

        ts_xxxxxz_xzzz[i] = 3.0 * ts_xxxxx_xzz[i] * fe_0 + ts_xxxxx_xzzz[i] * pa_z[i];

        ts_xxxxxz_yyyy[i] = ts_xxxxx_yyyy[i] * pa_z[i];

        ts_xxxxxz_yyyz[i] = 4.0 * ts_xxxz_yyyz[i] * fe_0 + ts_xxxxz_yyyz[i] * pa_x[i];

        ts_xxxxxz_yyzz[i] = 4.0 * ts_xxxz_yyzz[i] * fe_0 + ts_xxxxz_yyzz[i] * pa_x[i];

        ts_xxxxxz_yzzz[i] = 4.0 * ts_xxxz_yzzz[i] * fe_0 + ts_xxxxz_yzzz[i] * pa_x[i];

        ts_xxxxxz_zzzz[i] = 4.0 * ts_xxxz_zzzz[i] * fe_0 + ts_xxxxz_zzzz[i] * pa_x[i];
    }

    // Set up 45-60 components of targeted buffer : IG

    auto ts_xxxxyy_xxxx = pbuffer.data(idx_ovl_ig + 45);

    auto ts_xxxxyy_xxxy = pbuffer.data(idx_ovl_ig + 46);

    auto ts_xxxxyy_xxxz = pbuffer.data(idx_ovl_ig + 47);

    auto ts_xxxxyy_xxyy = pbuffer.data(idx_ovl_ig + 48);

    auto ts_xxxxyy_xxyz = pbuffer.data(idx_ovl_ig + 49);

    auto ts_xxxxyy_xxzz = pbuffer.data(idx_ovl_ig + 50);

    auto ts_xxxxyy_xyyy = pbuffer.data(idx_ovl_ig + 51);

    auto ts_xxxxyy_xyyz = pbuffer.data(idx_ovl_ig + 52);

    auto ts_xxxxyy_xyzz = pbuffer.data(idx_ovl_ig + 53);

    auto ts_xxxxyy_xzzz = pbuffer.data(idx_ovl_ig + 54);

    auto ts_xxxxyy_yyyy = pbuffer.data(idx_ovl_ig + 55);

    auto ts_xxxxyy_yyyz = pbuffer.data(idx_ovl_ig + 56);

    auto ts_xxxxyy_yyzz = pbuffer.data(idx_ovl_ig + 57);

    auto ts_xxxxyy_yzzz = pbuffer.data(idx_ovl_ig + 58);

    auto ts_xxxxyy_zzzz = pbuffer.data(idx_ovl_ig + 59);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxx_xxxx, ts_xxxx_xxxz, ts_xxxx_xxzz, ts_xxxx_xzzz, ts_xxxxy_xxxx, ts_xxxxy_xxxz, ts_xxxxy_xxzz, ts_xxxxy_xzzz, ts_xxxxyy_xxxx, ts_xxxxyy_xxxy, ts_xxxxyy_xxxz, ts_xxxxyy_xxyy, ts_xxxxyy_xxyz, ts_xxxxyy_xxzz, ts_xxxxyy_xyyy, ts_xxxxyy_xyyz, ts_xxxxyy_xyzz, ts_xxxxyy_xzzz, ts_xxxxyy_yyyy, ts_xxxxyy_yyyz, ts_xxxxyy_yyzz, ts_xxxxyy_yzzz, ts_xxxxyy_zzzz, ts_xxxyy_xxxy, ts_xxxyy_xxy, ts_xxxyy_xxyy, ts_xxxyy_xxyz, ts_xxxyy_xyy, ts_xxxyy_xyyy, ts_xxxyy_xyyz, ts_xxxyy_xyz, ts_xxxyy_xyzz, ts_xxxyy_yyy, ts_xxxyy_yyyy, ts_xxxyy_yyyz, ts_xxxyy_yyz, ts_xxxyy_yyzz, ts_xxxyy_yzz, ts_xxxyy_yzzz, ts_xxxyy_zzzz, ts_xxyy_xxxy, ts_xxyy_xxyy, ts_xxyy_xxyz, ts_xxyy_xyyy, ts_xxyy_xyyz, ts_xxyy_xyzz, ts_xxyy_yyyy, ts_xxyy_yyyz, ts_xxyy_yyzz, ts_xxyy_yzzz, ts_xxyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxyy_xxxx[i] = ts_xxxx_xxxx[i] * fe_0 + ts_xxxxy_xxxx[i] * pa_y[i];

        ts_xxxxyy_xxxy[i] = 3.0 * ts_xxyy_xxxy[i] * fe_0 + 3.0 * ts_xxxyy_xxy[i] * fe_0 + ts_xxxyy_xxxy[i] * pa_x[i];

        ts_xxxxyy_xxxz[i] = ts_xxxx_xxxz[i] * fe_0 + ts_xxxxy_xxxz[i] * pa_y[i];

        ts_xxxxyy_xxyy[i] = 3.0 * ts_xxyy_xxyy[i] * fe_0 + 2.0 * ts_xxxyy_xyy[i] * fe_0 + ts_xxxyy_xxyy[i] * pa_x[i];

        ts_xxxxyy_xxyz[i] = 3.0 * ts_xxyy_xxyz[i] * fe_0 + 2.0 * ts_xxxyy_xyz[i] * fe_0 + ts_xxxyy_xxyz[i] * pa_x[i];

        ts_xxxxyy_xxzz[i] = ts_xxxx_xxzz[i] * fe_0 + ts_xxxxy_xxzz[i] * pa_y[i];

        ts_xxxxyy_xyyy[i] = 3.0 * ts_xxyy_xyyy[i] * fe_0 + ts_xxxyy_yyy[i] * fe_0 + ts_xxxyy_xyyy[i] * pa_x[i];

        ts_xxxxyy_xyyz[i] = 3.0 * ts_xxyy_xyyz[i] * fe_0 + ts_xxxyy_yyz[i] * fe_0 + ts_xxxyy_xyyz[i] * pa_x[i];

        ts_xxxxyy_xyzz[i] = 3.0 * ts_xxyy_xyzz[i] * fe_0 + ts_xxxyy_yzz[i] * fe_0 + ts_xxxyy_xyzz[i] * pa_x[i];

        ts_xxxxyy_xzzz[i] = ts_xxxx_xzzz[i] * fe_0 + ts_xxxxy_xzzz[i] * pa_y[i];

        ts_xxxxyy_yyyy[i] = 3.0 * ts_xxyy_yyyy[i] * fe_0 + ts_xxxyy_yyyy[i] * pa_x[i];

        ts_xxxxyy_yyyz[i] = 3.0 * ts_xxyy_yyyz[i] * fe_0 + ts_xxxyy_yyyz[i] * pa_x[i];

        ts_xxxxyy_yyzz[i] = 3.0 * ts_xxyy_yyzz[i] * fe_0 + ts_xxxyy_yyzz[i] * pa_x[i];

        ts_xxxxyy_yzzz[i] = 3.0 * ts_xxyy_yzzz[i] * fe_0 + ts_xxxyy_yzzz[i] * pa_x[i];

        ts_xxxxyy_zzzz[i] = 3.0 * ts_xxyy_zzzz[i] * fe_0 + ts_xxxyy_zzzz[i] * pa_x[i];
    }

    // Set up 60-75 components of targeted buffer : IG

    auto ts_xxxxyz_xxxx = pbuffer.data(idx_ovl_ig + 60);

    auto ts_xxxxyz_xxxy = pbuffer.data(idx_ovl_ig + 61);

    auto ts_xxxxyz_xxxz = pbuffer.data(idx_ovl_ig + 62);

    auto ts_xxxxyz_xxyy = pbuffer.data(idx_ovl_ig + 63);

    auto ts_xxxxyz_xxyz = pbuffer.data(idx_ovl_ig + 64);

    auto ts_xxxxyz_xxzz = pbuffer.data(idx_ovl_ig + 65);

    auto ts_xxxxyz_xyyy = pbuffer.data(idx_ovl_ig + 66);

    auto ts_xxxxyz_xyyz = pbuffer.data(idx_ovl_ig + 67);

    auto ts_xxxxyz_xyzz = pbuffer.data(idx_ovl_ig + 68);

    auto ts_xxxxyz_xzzz = pbuffer.data(idx_ovl_ig + 69);

    auto ts_xxxxyz_yyyy = pbuffer.data(idx_ovl_ig + 70);

    auto ts_xxxxyz_yyyz = pbuffer.data(idx_ovl_ig + 71);

    auto ts_xxxxyz_yyzz = pbuffer.data(idx_ovl_ig + 72);

    auto ts_xxxxyz_yzzz = pbuffer.data(idx_ovl_ig + 73);

    auto ts_xxxxyz_zzzz = pbuffer.data(idx_ovl_ig + 74);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxxy_xxxy, ts_xxxxy_xxyy, ts_xxxxy_xyyy, ts_xxxxy_yyyy, ts_xxxxyz_xxxx, ts_xxxxyz_xxxy, ts_xxxxyz_xxxz, ts_xxxxyz_xxyy, ts_xxxxyz_xxyz, ts_xxxxyz_xxzz, ts_xxxxyz_xyyy, ts_xxxxyz_xyyz, ts_xxxxyz_xyzz, ts_xxxxyz_xzzz, ts_xxxxyz_yyyy, ts_xxxxyz_yyyz, ts_xxxxyz_yyzz, ts_xxxxyz_yzzz, ts_xxxxyz_zzzz, ts_xxxxz_xxxx, ts_xxxxz_xxxz, ts_xxxxz_xxyz, ts_xxxxz_xxz, ts_xxxxz_xxzz, ts_xxxxz_xyyz, ts_xxxxz_xyz, ts_xxxxz_xyzz, ts_xxxxz_xzz, ts_xxxxz_xzzz, ts_xxxxz_zzzz, ts_xxxyz_yyyz, ts_xxxyz_yyzz, ts_xxxyz_yzzz, ts_xxyz_yyyz, ts_xxyz_yyzz, ts_xxyz_yzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxyz_xxxx[i] = ts_xxxxz_xxxx[i] * pa_y[i];

        ts_xxxxyz_xxxy[i] = ts_xxxxy_xxxy[i] * pa_z[i];

        ts_xxxxyz_xxxz[i] = ts_xxxxz_xxxz[i] * pa_y[i];

        ts_xxxxyz_xxyy[i] = ts_xxxxy_xxyy[i] * pa_z[i];

        ts_xxxxyz_xxyz[i] = ts_xxxxz_xxz[i] * fe_0 + ts_xxxxz_xxyz[i] * pa_y[i];

        ts_xxxxyz_xxzz[i] = ts_xxxxz_xxzz[i] * pa_y[i];

        ts_xxxxyz_xyyy[i] = ts_xxxxy_xyyy[i] * pa_z[i];

        ts_xxxxyz_xyyz[i] = 2.0 * ts_xxxxz_xyz[i] * fe_0 + ts_xxxxz_xyyz[i] * pa_y[i];

        ts_xxxxyz_xyzz[i] = ts_xxxxz_xzz[i] * fe_0 + ts_xxxxz_xyzz[i] * pa_y[i];

        ts_xxxxyz_xzzz[i] = ts_xxxxz_xzzz[i] * pa_y[i];

        ts_xxxxyz_yyyy[i] = ts_xxxxy_yyyy[i] * pa_z[i];

        ts_xxxxyz_yyyz[i] = 3.0 * ts_xxyz_yyyz[i] * fe_0 + ts_xxxyz_yyyz[i] * pa_x[i];

        ts_xxxxyz_yyzz[i] = 3.0 * ts_xxyz_yyzz[i] * fe_0 + ts_xxxyz_yyzz[i] * pa_x[i];

        ts_xxxxyz_yzzz[i] = 3.0 * ts_xxyz_yzzz[i] * fe_0 + ts_xxxyz_yzzz[i] * pa_x[i];

        ts_xxxxyz_zzzz[i] = ts_xxxxz_zzzz[i] * pa_y[i];
    }

    // Set up 75-90 components of targeted buffer : IG

    auto ts_xxxxzz_xxxx = pbuffer.data(idx_ovl_ig + 75);

    auto ts_xxxxzz_xxxy = pbuffer.data(idx_ovl_ig + 76);

    auto ts_xxxxzz_xxxz = pbuffer.data(idx_ovl_ig + 77);

    auto ts_xxxxzz_xxyy = pbuffer.data(idx_ovl_ig + 78);

    auto ts_xxxxzz_xxyz = pbuffer.data(idx_ovl_ig + 79);

    auto ts_xxxxzz_xxzz = pbuffer.data(idx_ovl_ig + 80);

    auto ts_xxxxzz_xyyy = pbuffer.data(idx_ovl_ig + 81);

    auto ts_xxxxzz_xyyz = pbuffer.data(idx_ovl_ig + 82);

    auto ts_xxxxzz_xyzz = pbuffer.data(idx_ovl_ig + 83);

    auto ts_xxxxzz_xzzz = pbuffer.data(idx_ovl_ig + 84);

    auto ts_xxxxzz_yyyy = pbuffer.data(idx_ovl_ig + 85);

    auto ts_xxxxzz_yyyz = pbuffer.data(idx_ovl_ig + 86);

    auto ts_xxxxzz_yyzz = pbuffer.data(idx_ovl_ig + 87);

    auto ts_xxxxzz_yzzz = pbuffer.data(idx_ovl_ig + 88);

    auto ts_xxxxzz_zzzz = pbuffer.data(idx_ovl_ig + 89);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxyy, ts_xxxx_xyyy, ts_xxxxz_xxxx, ts_xxxxz_xxxy, ts_xxxxz_xxyy, ts_xxxxz_xyyy, ts_xxxxzz_xxxx, ts_xxxxzz_xxxy, ts_xxxxzz_xxxz, ts_xxxxzz_xxyy, ts_xxxxzz_xxyz, ts_xxxxzz_xxzz, ts_xxxxzz_xyyy, ts_xxxxzz_xyyz, ts_xxxxzz_xyzz, ts_xxxxzz_xzzz, ts_xxxxzz_yyyy, ts_xxxxzz_yyyz, ts_xxxxzz_yyzz, ts_xxxxzz_yzzz, ts_xxxxzz_zzzz, ts_xxxzz_xxxz, ts_xxxzz_xxyz, ts_xxxzz_xxz, ts_xxxzz_xxzz, ts_xxxzz_xyyz, ts_xxxzz_xyz, ts_xxxzz_xyzz, ts_xxxzz_xzz, ts_xxxzz_xzzz, ts_xxxzz_yyyy, ts_xxxzz_yyyz, ts_xxxzz_yyz, ts_xxxzz_yyzz, ts_xxxzz_yzz, ts_xxxzz_yzzz, ts_xxxzz_zzz, ts_xxxzz_zzzz, ts_xxzz_xxxz, ts_xxzz_xxyz, ts_xxzz_xxzz, ts_xxzz_xyyz, ts_xxzz_xyzz, ts_xxzz_xzzz, ts_xxzz_yyyy, ts_xxzz_yyyz, ts_xxzz_yyzz, ts_xxzz_yzzz, ts_xxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxzz_xxxx[i] = ts_xxxx_xxxx[i] * fe_0 + ts_xxxxz_xxxx[i] * pa_z[i];

        ts_xxxxzz_xxxy[i] = ts_xxxx_xxxy[i] * fe_0 + ts_xxxxz_xxxy[i] * pa_z[i];

        ts_xxxxzz_xxxz[i] = 3.0 * ts_xxzz_xxxz[i] * fe_0 + 3.0 * ts_xxxzz_xxz[i] * fe_0 + ts_xxxzz_xxxz[i] * pa_x[i];

        ts_xxxxzz_xxyy[i] = ts_xxxx_xxyy[i] * fe_0 + ts_xxxxz_xxyy[i] * pa_z[i];

        ts_xxxxzz_xxyz[i] = 3.0 * ts_xxzz_xxyz[i] * fe_0 + 2.0 * ts_xxxzz_xyz[i] * fe_0 + ts_xxxzz_xxyz[i] * pa_x[i];

        ts_xxxxzz_xxzz[i] = 3.0 * ts_xxzz_xxzz[i] * fe_0 + 2.0 * ts_xxxzz_xzz[i] * fe_0 + ts_xxxzz_xxzz[i] * pa_x[i];

        ts_xxxxzz_xyyy[i] = ts_xxxx_xyyy[i] * fe_0 + ts_xxxxz_xyyy[i] * pa_z[i];

        ts_xxxxzz_xyyz[i] = 3.0 * ts_xxzz_xyyz[i] * fe_0 + ts_xxxzz_yyz[i] * fe_0 + ts_xxxzz_xyyz[i] * pa_x[i];

        ts_xxxxzz_xyzz[i] = 3.0 * ts_xxzz_xyzz[i] * fe_0 + ts_xxxzz_yzz[i] * fe_0 + ts_xxxzz_xyzz[i] * pa_x[i];

        ts_xxxxzz_xzzz[i] = 3.0 * ts_xxzz_xzzz[i] * fe_0 + ts_xxxzz_zzz[i] * fe_0 + ts_xxxzz_xzzz[i] * pa_x[i];

        ts_xxxxzz_yyyy[i] = 3.0 * ts_xxzz_yyyy[i] * fe_0 + ts_xxxzz_yyyy[i] * pa_x[i];

        ts_xxxxzz_yyyz[i] = 3.0 * ts_xxzz_yyyz[i] * fe_0 + ts_xxxzz_yyyz[i] * pa_x[i];

        ts_xxxxzz_yyzz[i] = 3.0 * ts_xxzz_yyzz[i] * fe_0 + ts_xxxzz_yyzz[i] * pa_x[i];

        ts_xxxxzz_yzzz[i] = 3.0 * ts_xxzz_yzzz[i] * fe_0 + ts_xxxzz_yzzz[i] * pa_x[i];

        ts_xxxxzz_zzzz[i] = 3.0 * ts_xxzz_zzzz[i] * fe_0 + ts_xxxzz_zzzz[i] * pa_x[i];
    }

    // Set up 90-105 components of targeted buffer : IG

    auto ts_xxxyyy_xxxx = pbuffer.data(idx_ovl_ig + 90);

    auto ts_xxxyyy_xxxy = pbuffer.data(idx_ovl_ig + 91);

    auto ts_xxxyyy_xxxz = pbuffer.data(idx_ovl_ig + 92);

    auto ts_xxxyyy_xxyy = pbuffer.data(idx_ovl_ig + 93);

    auto ts_xxxyyy_xxyz = pbuffer.data(idx_ovl_ig + 94);

    auto ts_xxxyyy_xxzz = pbuffer.data(idx_ovl_ig + 95);

    auto ts_xxxyyy_xyyy = pbuffer.data(idx_ovl_ig + 96);

    auto ts_xxxyyy_xyyz = pbuffer.data(idx_ovl_ig + 97);

    auto ts_xxxyyy_xyzz = pbuffer.data(idx_ovl_ig + 98);

    auto ts_xxxyyy_xzzz = pbuffer.data(idx_ovl_ig + 99);

    auto ts_xxxyyy_yyyy = pbuffer.data(idx_ovl_ig + 100);

    auto ts_xxxyyy_yyyz = pbuffer.data(idx_ovl_ig + 101);

    auto ts_xxxyyy_yyzz = pbuffer.data(idx_ovl_ig + 102);

    auto ts_xxxyyy_yzzz = pbuffer.data(idx_ovl_ig + 103);

    auto ts_xxxyyy_zzzz = pbuffer.data(idx_ovl_ig + 104);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxy_xxxx, ts_xxxy_xxxz, ts_xxxy_xxzz, ts_xxxy_xzzz, ts_xxxyy_xxxx, ts_xxxyy_xxxz, ts_xxxyy_xxzz, ts_xxxyy_xzzz, ts_xxxyyy_xxxx, ts_xxxyyy_xxxy, ts_xxxyyy_xxxz, ts_xxxyyy_xxyy, ts_xxxyyy_xxyz, ts_xxxyyy_xxzz, ts_xxxyyy_xyyy, ts_xxxyyy_xyyz, ts_xxxyyy_xyzz, ts_xxxyyy_xzzz, ts_xxxyyy_yyyy, ts_xxxyyy_yyyz, ts_xxxyyy_yyzz, ts_xxxyyy_yzzz, ts_xxxyyy_zzzz, ts_xxyyy_xxxy, ts_xxyyy_xxy, ts_xxyyy_xxyy, ts_xxyyy_xxyz, ts_xxyyy_xyy, ts_xxyyy_xyyy, ts_xxyyy_xyyz, ts_xxyyy_xyz, ts_xxyyy_xyzz, ts_xxyyy_yyy, ts_xxyyy_yyyy, ts_xxyyy_yyyz, ts_xxyyy_yyz, ts_xxyyy_yyzz, ts_xxyyy_yzz, ts_xxyyy_yzzz, ts_xxyyy_zzzz, ts_xyyy_xxxy, ts_xyyy_xxyy, ts_xyyy_xxyz, ts_xyyy_xyyy, ts_xyyy_xyyz, ts_xyyy_xyzz, ts_xyyy_yyyy, ts_xyyy_yyyz, ts_xyyy_yyzz, ts_xyyy_yzzz, ts_xyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyyy_xxxx[i] = 2.0 * ts_xxxy_xxxx[i] * fe_0 + ts_xxxyy_xxxx[i] * pa_y[i];

        ts_xxxyyy_xxxy[i] = 2.0 * ts_xyyy_xxxy[i] * fe_0 + 3.0 * ts_xxyyy_xxy[i] * fe_0 + ts_xxyyy_xxxy[i] * pa_x[i];

        ts_xxxyyy_xxxz[i] = 2.0 * ts_xxxy_xxxz[i] * fe_0 + ts_xxxyy_xxxz[i] * pa_y[i];

        ts_xxxyyy_xxyy[i] = 2.0 * ts_xyyy_xxyy[i] * fe_0 + 2.0 * ts_xxyyy_xyy[i] * fe_0 + ts_xxyyy_xxyy[i] * pa_x[i];

        ts_xxxyyy_xxyz[i] = 2.0 * ts_xyyy_xxyz[i] * fe_0 + 2.0 * ts_xxyyy_xyz[i] * fe_0 + ts_xxyyy_xxyz[i] * pa_x[i];

        ts_xxxyyy_xxzz[i] = 2.0 * ts_xxxy_xxzz[i] * fe_0 + ts_xxxyy_xxzz[i] * pa_y[i];

        ts_xxxyyy_xyyy[i] = 2.0 * ts_xyyy_xyyy[i] * fe_0 + ts_xxyyy_yyy[i] * fe_0 + ts_xxyyy_xyyy[i] * pa_x[i];

        ts_xxxyyy_xyyz[i] = 2.0 * ts_xyyy_xyyz[i] * fe_0 + ts_xxyyy_yyz[i] * fe_0 + ts_xxyyy_xyyz[i] * pa_x[i];

        ts_xxxyyy_xyzz[i] = 2.0 * ts_xyyy_xyzz[i] * fe_0 + ts_xxyyy_yzz[i] * fe_0 + ts_xxyyy_xyzz[i] * pa_x[i];

        ts_xxxyyy_xzzz[i] = 2.0 * ts_xxxy_xzzz[i] * fe_0 + ts_xxxyy_xzzz[i] * pa_y[i];

        ts_xxxyyy_yyyy[i] = 2.0 * ts_xyyy_yyyy[i] * fe_0 + ts_xxyyy_yyyy[i] * pa_x[i];

        ts_xxxyyy_yyyz[i] = 2.0 * ts_xyyy_yyyz[i] * fe_0 + ts_xxyyy_yyyz[i] * pa_x[i];

        ts_xxxyyy_yyzz[i] = 2.0 * ts_xyyy_yyzz[i] * fe_0 + ts_xxyyy_yyzz[i] * pa_x[i];

        ts_xxxyyy_yzzz[i] = 2.0 * ts_xyyy_yzzz[i] * fe_0 + ts_xxyyy_yzzz[i] * pa_x[i];

        ts_xxxyyy_zzzz[i] = 2.0 * ts_xyyy_zzzz[i] * fe_0 + ts_xxyyy_zzzz[i] * pa_x[i];
    }

    // Set up 105-120 components of targeted buffer : IG

    auto ts_xxxyyz_xxxx = pbuffer.data(idx_ovl_ig + 105);

    auto ts_xxxyyz_xxxy = pbuffer.data(idx_ovl_ig + 106);

    auto ts_xxxyyz_xxxz = pbuffer.data(idx_ovl_ig + 107);

    auto ts_xxxyyz_xxyy = pbuffer.data(idx_ovl_ig + 108);

    auto ts_xxxyyz_xxyz = pbuffer.data(idx_ovl_ig + 109);

    auto ts_xxxyyz_xxzz = pbuffer.data(idx_ovl_ig + 110);

    auto ts_xxxyyz_xyyy = pbuffer.data(idx_ovl_ig + 111);

    auto ts_xxxyyz_xyyz = pbuffer.data(idx_ovl_ig + 112);

    auto ts_xxxyyz_xyzz = pbuffer.data(idx_ovl_ig + 113);

    auto ts_xxxyyz_xzzz = pbuffer.data(idx_ovl_ig + 114);

    auto ts_xxxyyz_yyyy = pbuffer.data(idx_ovl_ig + 115);

    auto ts_xxxyyz_yyyz = pbuffer.data(idx_ovl_ig + 116);

    auto ts_xxxyyz_yyzz = pbuffer.data(idx_ovl_ig + 117);

    auto ts_xxxyyz_yzzz = pbuffer.data(idx_ovl_ig + 118);

    auto ts_xxxyyz_zzzz = pbuffer.data(idx_ovl_ig + 119);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxyy_xxxx, ts_xxxyy_xxxy, ts_xxxyy_xxy, ts_xxxyy_xxyy, ts_xxxyy_xxyz, ts_xxxyy_xyy, ts_xxxyy_xyyy, ts_xxxyy_xyyz, ts_xxxyy_xyz, ts_xxxyy_xyzz, ts_xxxyy_yyyy, ts_xxxyyz_xxxx, ts_xxxyyz_xxxy, ts_xxxyyz_xxxz, ts_xxxyyz_xxyy, ts_xxxyyz_xxyz, ts_xxxyyz_xxzz, ts_xxxyyz_xyyy, ts_xxxyyz_xyyz, ts_xxxyyz_xyzz, ts_xxxyyz_xzzz, ts_xxxyyz_yyyy, ts_xxxyyz_yyyz, ts_xxxyyz_yyzz, ts_xxxyyz_yzzz, ts_xxxyyz_zzzz, ts_xxxyz_xxxz, ts_xxxyz_xxzz, ts_xxxyz_xzzz, ts_xxxz_xxxz, ts_xxxz_xxzz, ts_xxxz_xzzz, ts_xxyyz_yyyz, ts_xxyyz_yyzz, ts_xxyyz_yzzz, ts_xxyyz_zzzz, ts_xyyz_yyyz, ts_xyyz_yyzz, ts_xyyz_yzzz, ts_xyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyyz_xxxx[i] = ts_xxxyy_xxxx[i] * pa_z[i];

        ts_xxxyyz_xxxy[i] = ts_xxxyy_xxxy[i] * pa_z[i];

        ts_xxxyyz_xxxz[i] = ts_xxxz_xxxz[i] * fe_0 + ts_xxxyz_xxxz[i] * pa_y[i];

        ts_xxxyyz_xxyy[i] = ts_xxxyy_xxyy[i] * pa_z[i];

        ts_xxxyyz_xxyz[i] = ts_xxxyy_xxy[i] * fe_0 + ts_xxxyy_xxyz[i] * pa_z[i];

        ts_xxxyyz_xxzz[i] = ts_xxxz_xxzz[i] * fe_0 + ts_xxxyz_xxzz[i] * pa_y[i];

        ts_xxxyyz_xyyy[i] = ts_xxxyy_xyyy[i] * pa_z[i];

        ts_xxxyyz_xyyz[i] = ts_xxxyy_xyy[i] * fe_0 + ts_xxxyy_xyyz[i] * pa_z[i];

        ts_xxxyyz_xyzz[i] = 2.0 * ts_xxxyy_xyz[i] * fe_0 + ts_xxxyy_xyzz[i] * pa_z[i];

        ts_xxxyyz_xzzz[i] = ts_xxxz_xzzz[i] * fe_0 + ts_xxxyz_xzzz[i] * pa_y[i];

        ts_xxxyyz_yyyy[i] = ts_xxxyy_yyyy[i] * pa_z[i];

        ts_xxxyyz_yyyz[i] = 2.0 * ts_xyyz_yyyz[i] * fe_0 + ts_xxyyz_yyyz[i] * pa_x[i];

        ts_xxxyyz_yyzz[i] = 2.0 * ts_xyyz_yyzz[i] * fe_0 + ts_xxyyz_yyzz[i] * pa_x[i];

        ts_xxxyyz_yzzz[i] = 2.0 * ts_xyyz_yzzz[i] * fe_0 + ts_xxyyz_yzzz[i] * pa_x[i];

        ts_xxxyyz_zzzz[i] = 2.0 * ts_xyyz_zzzz[i] * fe_0 + ts_xxyyz_zzzz[i] * pa_x[i];
    }

    // Set up 120-135 components of targeted buffer : IG

    auto ts_xxxyzz_xxxx = pbuffer.data(idx_ovl_ig + 120);

    auto ts_xxxyzz_xxxy = pbuffer.data(idx_ovl_ig + 121);

    auto ts_xxxyzz_xxxz = pbuffer.data(idx_ovl_ig + 122);

    auto ts_xxxyzz_xxyy = pbuffer.data(idx_ovl_ig + 123);

    auto ts_xxxyzz_xxyz = pbuffer.data(idx_ovl_ig + 124);

    auto ts_xxxyzz_xxzz = pbuffer.data(idx_ovl_ig + 125);

    auto ts_xxxyzz_xyyy = pbuffer.data(idx_ovl_ig + 126);

    auto ts_xxxyzz_xyyz = pbuffer.data(idx_ovl_ig + 127);

    auto ts_xxxyzz_xyzz = pbuffer.data(idx_ovl_ig + 128);

    auto ts_xxxyzz_xzzz = pbuffer.data(idx_ovl_ig + 129);

    auto ts_xxxyzz_yyyy = pbuffer.data(idx_ovl_ig + 130);

    auto ts_xxxyzz_yyyz = pbuffer.data(idx_ovl_ig + 131);

    auto ts_xxxyzz_yyzz = pbuffer.data(idx_ovl_ig + 132);

    auto ts_xxxyzz_yzzz = pbuffer.data(idx_ovl_ig + 133);

    auto ts_xxxyzz_zzzz = pbuffer.data(idx_ovl_ig + 134);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxyzz_xxxx, ts_xxxyzz_xxxy, ts_xxxyzz_xxxz, ts_xxxyzz_xxyy, ts_xxxyzz_xxyz, ts_xxxyzz_xxzz, ts_xxxyzz_xyyy, ts_xxxyzz_xyyz, ts_xxxyzz_xyzz, ts_xxxyzz_xzzz, ts_xxxyzz_yyyy, ts_xxxyzz_yyyz, ts_xxxyzz_yyzz, ts_xxxyzz_yzzz, ts_xxxyzz_zzzz, ts_xxxzz_xxx, ts_xxxzz_xxxx, ts_xxxzz_xxxy, ts_xxxzz_xxxz, ts_xxxzz_xxy, ts_xxxzz_xxyy, ts_xxxzz_xxyz, ts_xxxzz_xxz, ts_xxxzz_xxzz, ts_xxxzz_xyy, ts_xxxzz_xyyy, ts_xxxzz_xyyz, ts_xxxzz_xyz, ts_xxxzz_xyzz, ts_xxxzz_xzz, ts_xxxzz_xzzz, ts_xxxzz_zzzz, ts_xxyzz_yyyy, ts_xxyzz_yyyz, ts_xxyzz_yyzz, ts_xxyzz_yzzz, ts_xyzz_yyyy, ts_xyzz_yyyz, ts_xyzz_yyzz, ts_xyzz_yzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyzz_xxxx[i] = ts_xxxzz_xxxx[i] * pa_y[i];

        ts_xxxyzz_xxxy[i] = ts_xxxzz_xxx[i] * fe_0 + ts_xxxzz_xxxy[i] * pa_y[i];

        ts_xxxyzz_xxxz[i] = ts_xxxzz_xxxz[i] * pa_y[i];

        ts_xxxyzz_xxyy[i] = 2.0 * ts_xxxzz_xxy[i] * fe_0 + ts_xxxzz_xxyy[i] * pa_y[i];

        ts_xxxyzz_xxyz[i] = ts_xxxzz_xxz[i] * fe_0 + ts_xxxzz_xxyz[i] * pa_y[i];

        ts_xxxyzz_xxzz[i] = ts_xxxzz_xxzz[i] * pa_y[i];

        ts_xxxyzz_xyyy[i] = 3.0 * ts_xxxzz_xyy[i] * fe_0 + ts_xxxzz_xyyy[i] * pa_y[i];

        ts_xxxyzz_xyyz[i] = 2.0 * ts_xxxzz_xyz[i] * fe_0 + ts_xxxzz_xyyz[i] * pa_y[i];

        ts_xxxyzz_xyzz[i] = ts_xxxzz_xzz[i] * fe_0 + ts_xxxzz_xyzz[i] * pa_y[i];

        ts_xxxyzz_xzzz[i] = ts_xxxzz_xzzz[i] * pa_y[i];

        ts_xxxyzz_yyyy[i] = 2.0 * ts_xyzz_yyyy[i] * fe_0 + ts_xxyzz_yyyy[i] * pa_x[i];

        ts_xxxyzz_yyyz[i] = 2.0 * ts_xyzz_yyyz[i] * fe_0 + ts_xxyzz_yyyz[i] * pa_x[i];

        ts_xxxyzz_yyzz[i] = 2.0 * ts_xyzz_yyzz[i] * fe_0 + ts_xxyzz_yyzz[i] * pa_x[i];

        ts_xxxyzz_yzzz[i] = 2.0 * ts_xyzz_yzzz[i] * fe_0 + ts_xxyzz_yzzz[i] * pa_x[i];

        ts_xxxyzz_zzzz[i] = ts_xxxzz_zzzz[i] * pa_y[i];
    }

    // Set up 135-150 components of targeted buffer : IG

    auto ts_xxxzzz_xxxx = pbuffer.data(idx_ovl_ig + 135);

    auto ts_xxxzzz_xxxy = pbuffer.data(idx_ovl_ig + 136);

    auto ts_xxxzzz_xxxz = pbuffer.data(idx_ovl_ig + 137);

    auto ts_xxxzzz_xxyy = pbuffer.data(idx_ovl_ig + 138);

    auto ts_xxxzzz_xxyz = pbuffer.data(idx_ovl_ig + 139);

    auto ts_xxxzzz_xxzz = pbuffer.data(idx_ovl_ig + 140);

    auto ts_xxxzzz_xyyy = pbuffer.data(idx_ovl_ig + 141);

    auto ts_xxxzzz_xyyz = pbuffer.data(idx_ovl_ig + 142);

    auto ts_xxxzzz_xyzz = pbuffer.data(idx_ovl_ig + 143);

    auto ts_xxxzzz_xzzz = pbuffer.data(idx_ovl_ig + 144);

    auto ts_xxxzzz_yyyy = pbuffer.data(idx_ovl_ig + 145);

    auto ts_xxxzzz_yyyz = pbuffer.data(idx_ovl_ig + 146);

    auto ts_xxxzzz_yyzz = pbuffer.data(idx_ovl_ig + 147);

    auto ts_xxxzzz_yzzz = pbuffer.data(idx_ovl_ig + 148);

    auto ts_xxxzzz_zzzz = pbuffer.data(idx_ovl_ig + 149);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxz_xxxx, ts_xxxz_xxxy, ts_xxxz_xxyy, ts_xxxz_xyyy, ts_xxxzz_xxxx, ts_xxxzz_xxxy, ts_xxxzz_xxyy, ts_xxxzz_xyyy, ts_xxxzzz_xxxx, ts_xxxzzz_xxxy, ts_xxxzzz_xxxz, ts_xxxzzz_xxyy, ts_xxxzzz_xxyz, ts_xxxzzz_xxzz, ts_xxxzzz_xyyy, ts_xxxzzz_xyyz, ts_xxxzzz_xyzz, ts_xxxzzz_xzzz, ts_xxxzzz_yyyy, ts_xxxzzz_yyyz, ts_xxxzzz_yyzz, ts_xxxzzz_yzzz, ts_xxxzzz_zzzz, ts_xxzzz_xxxz, ts_xxzzz_xxyz, ts_xxzzz_xxz, ts_xxzzz_xxzz, ts_xxzzz_xyyz, ts_xxzzz_xyz, ts_xxzzz_xyzz, ts_xxzzz_xzz, ts_xxzzz_xzzz, ts_xxzzz_yyyy, ts_xxzzz_yyyz, ts_xxzzz_yyz, ts_xxzzz_yyzz, ts_xxzzz_yzz, ts_xxzzz_yzzz, ts_xxzzz_zzz, ts_xxzzz_zzzz, ts_xzzz_xxxz, ts_xzzz_xxyz, ts_xzzz_xxzz, ts_xzzz_xyyz, ts_xzzz_xyzz, ts_xzzz_xzzz, ts_xzzz_yyyy, ts_xzzz_yyyz, ts_xzzz_yyzz, ts_xzzz_yzzz, ts_xzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxzzz_xxxx[i] = 2.0 * ts_xxxz_xxxx[i] * fe_0 + ts_xxxzz_xxxx[i] * pa_z[i];

        ts_xxxzzz_xxxy[i] = 2.0 * ts_xxxz_xxxy[i] * fe_0 + ts_xxxzz_xxxy[i] * pa_z[i];

        ts_xxxzzz_xxxz[i] = 2.0 * ts_xzzz_xxxz[i] * fe_0 + 3.0 * ts_xxzzz_xxz[i] * fe_0 + ts_xxzzz_xxxz[i] * pa_x[i];

        ts_xxxzzz_xxyy[i] = 2.0 * ts_xxxz_xxyy[i] * fe_0 + ts_xxxzz_xxyy[i] * pa_z[i];

        ts_xxxzzz_xxyz[i] = 2.0 * ts_xzzz_xxyz[i] * fe_0 + 2.0 * ts_xxzzz_xyz[i] * fe_0 + ts_xxzzz_xxyz[i] * pa_x[i];

        ts_xxxzzz_xxzz[i] = 2.0 * ts_xzzz_xxzz[i] * fe_0 + 2.0 * ts_xxzzz_xzz[i] * fe_0 + ts_xxzzz_xxzz[i] * pa_x[i];

        ts_xxxzzz_xyyy[i] = 2.0 * ts_xxxz_xyyy[i] * fe_0 + ts_xxxzz_xyyy[i] * pa_z[i];

        ts_xxxzzz_xyyz[i] = 2.0 * ts_xzzz_xyyz[i] * fe_0 + ts_xxzzz_yyz[i] * fe_0 + ts_xxzzz_xyyz[i] * pa_x[i];

        ts_xxxzzz_xyzz[i] = 2.0 * ts_xzzz_xyzz[i] * fe_0 + ts_xxzzz_yzz[i] * fe_0 + ts_xxzzz_xyzz[i] * pa_x[i];

        ts_xxxzzz_xzzz[i] = 2.0 * ts_xzzz_xzzz[i] * fe_0 + ts_xxzzz_zzz[i] * fe_0 + ts_xxzzz_xzzz[i] * pa_x[i];

        ts_xxxzzz_yyyy[i] = 2.0 * ts_xzzz_yyyy[i] * fe_0 + ts_xxzzz_yyyy[i] * pa_x[i];

        ts_xxxzzz_yyyz[i] = 2.0 * ts_xzzz_yyyz[i] * fe_0 + ts_xxzzz_yyyz[i] * pa_x[i];

        ts_xxxzzz_yyzz[i] = 2.0 * ts_xzzz_yyzz[i] * fe_0 + ts_xxzzz_yyzz[i] * pa_x[i];

        ts_xxxzzz_yzzz[i] = 2.0 * ts_xzzz_yzzz[i] * fe_0 + ts_xxzzz_yzzz[i] * pa_x[i];

        ts_xxxzzz_zzzz[i] = 2.0 * ts_xzzz_zzzz[i] * fe_0 + ts_xxzzz_zzzz[i] * pa_x[i];
    }

    // Set up 150-165 components of targeted buffer : IG

    auto ts_xxyyyy_xxxx = pbuffer.data(idx_ovl_ig + 150);

    auto ts_xxyyyy_xxxy = pbuffer.data(idx_ovl_ig + 151);

    auto ts_xxyyyy_xxxz = pbuffer.data(idx_ovl_ig + 152);

    auto ts_xxyyyy_xxyy = pbuffer.data(idx_ovl_ig + 153);

    auto ts_xxyyyy_xxyz = pbuffer.data(idx_ovl_ig + 154);

    auto ts_xxyyyy_xxzz = pbuffer.data(idx_ovl_ig + 155);

    auto ts_xxyyyy_xyyy = pbuffer.data(idx_ovl_ig + 156);

    auto ts_xxyyyy_xyyz = pbuffer.data(idx_ovl_ig + 157);

    auto ts_xxyyyy_xyzz = pbuffer.data(idx_ovl_ig + 158);

    auto ts_xxyyyy_xzzz = pbuffer.data(idx_ovl_ig + 159);

    auto ts_xxyyyy_yyyy = pbuffer.data(idx_ovl_ig + 160);

    auto ts_xxyyyy_yyyz = pbuffer.data(idx_ovl_ig + 161);

    auto ts_xxyyyy_yyzz = pbuffer.data(idx_ovl_ig + 162);

    auto ts_xxyyyy_yzzz = pbuffer.data(idx_ovl_ig + 163);

    auto ts_xxyyyy_zzzz = pbuffer.data(idx_ovl_ig + 164);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyy_xxxx, ts_xxyy_xxxz, ts_xxyy_xxzz, ts_xxyy_xzzz, ts_xxyyy_xxxx, ts_xxyyy_xxxz, ts_xxyyy_xxzz, ts_xxyyy_xzzz, ts_xxyyyy_xxxx, ts_xxyyyy_xxxy, ts_xxyyyy_xxxz, ts_xxyyyy_xxyy, ts_xxyyyy_xxyz, ts_xxyyyy_xxzz, ts_xxyyyy_xyyy, ts_xxyyyy_xyyz, ts_xxyyyy_xyzz, ts_xxyyyy_xzzz, ts_xxyyyy_yyyy, ts_xxyyyy_yyyz, ts_xxyyyy_yyzz, ts_xxyyyy_yzzz, ts_xxyyyy_zzzz, ts_xyyyy_xxxy, ts_xyyyy_xxy, ts_xyyyy_xxyy, ts_xyyyy_xxyz, ts_xyyyy_xyy, ts_xyyyy_xyyy, ts_xyyyy_xyyz, ts_xyyyy_xyz, ts_xyyyy_xyzz, ts_xyyyy_yyy, ts_xyyyy_yyyy, ts_xyyyy_yyyz, ts_xyyyy_yyz, ts_xyyyy_yyzz, ts_xyyyy_yzz, ts_xyyyy_yzzz, ts_xyyyy_zzzz, ts_yyyy_xxxy, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyzz, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyzz, ts_yyyy_yzzz, ts_yyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyyy_xxxx[i] = 3.0 * ts_xxyy_xxxx[i] * fe_0 + ts_xxyyy_xxxx[i] * pa_y[i];

        ts_xxyyyy_xxxy[i] = ts_yyyy_xxxy[i] * fe_0 + 3.0 * ts_xyyyy_xxy[i] * fe_0 + ts_xyyyy_xxxy[i] * pa_x[i];

        ts_xxyyyy_xxxz[i] = 3.0 * ts_xxyy_xxxz[i] * fe_0 + ts_xxyyy_xxxz[i] * pa_y[i];

        ts_xxyyyy_xxyy[i] = ts_yyyy_xxyy[i] * fe_0 + 2.0 * ts_xyyyy_xyy[i] * fe_0 + ts_xyyyy_xxyy[i] * pa_x[i];

        ts_xxyyyy_xxyz[i] = ts_yyyy_xxyz[i] * fe_0 + 2.0 * ts_xyyyy_xyz[i] * fe_0 + ts_xyyyy_xxyz[i] * pa_x[i];

        ts_xxyyyy_xxzz[i] = 3.0 * ts_xxyy_xxzz[i] * fe_0 + ts_xxyyy_xxzz[i] * pa_y[i];

        ts_xxyyyy_xyyy[i] = ts_yyyy_xyyy[i] * fe_0 + ts_xyyyy_yyy[i] * fe_0 + ts_xyyyy_xyyy[i] * pa_x[i];

        ts_xxyyyy_xyyz[i] = ts_yyyy_xyyz[i] * fe_0 + ts_xyyyy_yyz[i] * fe_0 + ts_xyyyy_xyyz[i] * pa_x[i];

        ts_xxyyyy_xyzz[i] = ts_yyyy_xyzz[i] * fe_0 + ts_xyyyy_yzz[i] * fe_0 + ts_xyyyy_xyzz[i] * pa_x[i];

        ts_xxyyyy_xzzz[i] = 3.0 * ts_xxyy_xzzz[i] * fe_0 + ts_xxyyy_xzzz[i] * pa_y[i];

        ts_xxyyyy_yyyy[i] = ts_yyyy_yyyy[i] * fe_0 + ts_xyyyy_yyyy[i] * pa_x[i];

        ts_xxyyyy_yyyz[i] = ts_yyyy_yyyz[i] * fe_0 + ts_xyyyy_yyyz[i] * pa_x[i];

        ts_xxyyyy_yyzz[i] = ts_yyyy_yyzz[i] * fe_0 + ts_xyyyy_yyzz[i] * pa_x[i];

        ts_xxyyyy_yzzz[i] = ts_yyyy_yzzz[i] * fe_0 + ts_xyyyy_yzzz[i] * pa_x[i];

        ts_xxyyyy_zzzz[i] = ts_yyyy_zzzz[i] * fe_0 + ts_xyyyy_zzzz[i] * pa_x[i];
    }

    // Set up 165-180 components of targeted buffer : IG

    auto ts_xxyyyz_xxxx = pbuffer.data(idx_ovl_ig + 165);

    auto ts_xxyyyz_xxxy = pbuffer.data(idx_ovl_ig + 166);

    auto ts_xxyyyz_xxxz = pbuffer.data(idx_ovl_ig + 167);

    auto ts_xxyyyz_xxyy = pbuffer.data(idx_ovl_ig + 168);

    auto ts_xxyyyz_xxyz = pbuffer.data(idx_ovl_ig + 169);

    auto ts_xxyyyz_xxzz = pbuffer.data(idx_ovl_ig + 170);

    auto ts_xxyyyz_xyyy = pbuffer.data(idx_ovl_ig + 171);

    auto ts_xxyyyz_xyyz = pbuffer.data(idx_ovl_ig + 172);

    auto ts_xxyyyz_xyzz = pbuffer.data(idx_ovl_ig + 173);

    auto ts_xxyyyz_xzzz = pbuffer.data(idx_ovl_ig + 174);

    auto ts_xxyyyz_yyyy = pbuffer.data(idx_ovl_ig + 175);

    auto ts_xxyyyz_yyyz = pbuffer.data(idx_ovl_ig + 176);

    auto ts_xxyyyz_yyzz = pbuffer.data(idx_ovl_ig + 177);

    auto ts_xxyyyz_yzzz = pbuffer.data(idx_ovl_ig + 178);

    auto ts_xxyyyz_zzzz = pbuffer.data(idx_ovl_ig + 179);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyyy_xxxx, ts_xxyyy_xxxy, ts_xxyyy_xxy, ts_xxyyy_xxyy, ts_xxyyy_xxyz, ts_xxyyy_xyy, ts_xxyyy_xyyy, ts_xxyyy_xyyz, ts_xxyyy_xyz, ts_xxyyy_xyzz, ts_xxyyy_yyyy, ts_xxyyyz_xxxx, ts_xxyyyz_xxxy, ts_xxyyyz_xxxz, ts_xxyyyz_xxyy, ts_xxyyyz_xxyz, ts_xxyyyz_xxzz, ts_xxyyyz_xyyy, ts_xxyyyz_xyyz, ts_xxyyyz_xyzz, ts_xxyyyz_xzzz, ts_xxyyyz_yyyy, ts_xxyyyz_yyyz, ts_xxyyyz_yyzz, ts_xxyyyz_yzzz, ts_xxyyyz_zzzz, ts_xxyyz_xxxz, ts_xxyyz_xxzz, ts_xxyyz_xzzz, ts_xxyz_xxxz, ts_xxyz_xxzz, ts_xxyz_xzzz, ts_xyyyz_yyyz, ts_xyyyz_yyzz, ts_xyyyz_yzzz, ts_xyyyz_zzzz, ts_yyyz_yyyz, ts_yyyz_yyzz, ts_yyyz_yzzz, ts_yyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyyz_xxxx[i] = ts_xxyyy_xxxx[i] * pa_z[i];

        ts_xxyyyz_xxxy[i] = ts_xxyyy_xxxy[i] * pa_z[i];

        ts_xxyyyz_xxxz[i] = 2.0 * ts_xxyz_xxxz[i] * fe_0 + ts_xxyyz_xxxz[i] * pa_y[i];

        ts_xxyyyz_xxyy[i] = ts_xxyyy_xxyy[i] * pa_z[i];

        ts_xxyyyz_xxyz[i] = ts_xxyyy_xxy[i] * fe_0 + ts_xxyyy_xxyz[i] * pa_z[i];

        ts_xxyyyz_xxzz[i] = 2.0 * ts_xxyz_xxzz[i] * fe_0 + ts_xxyyz_xxzz[i] * pa_y[i];

        ts_xxyyyz_xyyy[i] = ts_xxyyy_xyyy[i] * pa_z[i];

        ts_xxyyyz_xyyz[i] = ts_xxyyy_xyy[i] * fe_0 + ts_xxyyy_xyyz[i] * pa_z[i];

        ts_xxyyyz_xyzz[i] = 2.0 * ts_xxyyy_xyz[i] * fe_0 + ts_xxyyy_xyzz[i] * pa_z[i];

        ts_xxyyyz_xzzz[i] = 2.0 * ts_xxyz_xzzz[i] * fe_0 + ts_xxyyz_xzzz[i] * pa_y[i];

        ts_xxyyyz_yyyy[i] = ts_xxyyy_yyyy[i] * pa_z[i];

        ts_xxyyyz_yyyz[i] = ts_yyyz_yyyz[i] * fe_0 + ts_xyyyz_yyyz[i] * pa_x[i];

        ts_xxyyyz_yyzz[i] = ts_yyyz_yyzz[i] * fe_0 + ts_xyyyz_yyzz[i] * pa_x[i];

        ts_xxyyyz_yzzz[i] = ts_yyyz_yzzz[i] * fe_0 + ts_xyyyz_yzzz[i] * pa_x[i];

        ts_xxyyyz_zzzz[i] = ts_yyyz_zzzz[i] * fe_0 + ts_xyyyz_zzzz[i] * pa_x[i];
    }

    // Set up 180-195 components of targeted buffer : IG

    auto ts_xxyyzz_xxxx = pbuffer.data(idx_ovl_ig + 180);

    auto ts_xxyyzz_xxxy = pbuffer.data(idx_ovl_ig + 181);

    auto ts_xxyyzz_xxxz = pbuffer.data(idx_ovl_ig + 182);

    auto ts_xxyyzz_xxyy = pbuffer.data(idx_ovl_ig + 183);

    auto ts_xxyyzz_xxyz = pbuffer.data(idx_ovl_ig + 184);

    auto ts_xxyyzz_xxzz = pbuffer.data(idx_ovl_ig + 185);

    auto ts_xxyyzz_xyyy = pbuffer.data(idx_ovl_ig + 186);

    auto ts_xxyyzz_xyyz = pbuffer.data(idx_ovl_ig + 187);

    auto ts_xxyyzz_xyzz = pbuffer.data(idx_ovl_ig + 188);

    auto ts_xxyyzz_xzzz = pbuffer.data(idx_ovl_ig + 189);

    auto ts_xxyyzz_yyyy = pbuffer.data(idx_ovl_ig + 190);

    auto ts_xxyyzz_yyyz = pbuffer.data(idx_ovl_ig + 191);

    auto ts_xxyyzz_yyzz = pbuffer.data(idx_ovl_ig + 192);

    auto ts_xxyyzz_yzzz = pbuffer.data(idx_ovl_ig + 193);

    auto ts_xxyyzz_zzzz = pbuffer.data(idx_ovl_ig + 194);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyy_xxxy, ts_xxyy_xxyy, ts_xxyy_xyyy, ts_xxyyz_xxxy, ts_xxyyz_xxyy, ts_xxyyz_xyyy, ts_xxyyzz_xxxx, ts_xxyyzz_xxxy, ts_xxyyzz_xxxz, ts_xxyyzz_xxyy, ts_xxyyzz_xxyz, ts_xxyyzz_xxzz, ts_xxyyzz_xyyy, ts_xxyyzz_xyyz, ts_xxyyzz_xyzz, ts_xxyyzz_xzzz, ts_xxyyzz_yyyy, ts_xxyyzz_yyyz, ts_xxyyzz_yyzz, ts_xxyyzz_yzzz, ts_xxyyzz_zzzz, ts_xxyzz_xxxx, ts_xxyzz_xxxz, ts_xxyzz_xxzz, ts_xxyzz_xzzz, ts_xxzz_xxxx, ts_xxzz_xxxz, ts_xxzz_xxzz, ts_xxzz_xzzz, ts_xyyzz_xxyz, ts_xyyzz_xyyz, ts_xyyzz_xyz, ts_xyyzz_xyzz, ts_xyyzz_yyyy, ts_xyyzz_yyyz, ts_xyyzz_yyz, ts_xyyzz_yyzz, ts_xyyzz_yzz, ts_xyyzz_yzzz, ts_xyyzz_zzzz, ts_yyzz_xxyz, ts_yyzz_xyyz, ts_yyzz_xyzz, ts_yyzz_yyyy, ts_yyzz_yyyz, ts_yyzz_yyzz, ts_yyzz_yzzz, ts_yyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyzz_xxxx[i] = ts_xxzz_xxxx[i] * fe_0 + ts_xxyzz_xxxx[i] * pa_y[i];

        ts_xxyyzz_xxxy[i] = ts_xxyy_xxxy[i] * fe_0 + ts_xxyyz_xxxy[i] * pa_z[i];

        ts_xxyyzz_xxxz[i] = ts_xxzz_xxxz[i] * fe_0 + ts_xxyzz_xxxz[i] * pa_y[i];

        ts_xxyyzz_xxyy[i] = ts_xxyy_xxyy[i] * fe_0 + ts_xxyyz_xxyy[i] * pa_z[i];

        ts_xxyyzz_xxyz[i] = ts_yyzz_xxyz[i] * fe_0 + 2.0 * ts_xyyzz_xyz[i] * fe_0 + ts_xyyzz_xxyz[i] * pa_x[i];

        ts_xxyyzz_xxzz[i] = ts_xxzz_xxzz[i] * fe_0 + ts_xxyzz_xxzz[i] * pa_y[i];

        ts_xxyyzz_xyyy[i] = ts_xxyy_xyyy[i] * fe_0 + ts_xxyyz_xyyy[i] * pa_z[i];

        ts_xxyyzz_xyyz[i] = ts_yyzz_xyyz[i] * fe_0 + ts_xyyzz_yyz[i] * fe_0 + ts_xyyzz_xyyz[i] * pa_x[i];

        ts_xxyyzz_xyzz[i] = ts_yyzz_xyzz[i] * fe_0 + ts_xyyzz_yzz[i] * fe_0 + ts_xyyzz_xyzz[i] * pa_x[i];

        ts_xxyyzz_xzzz[i] = ts_xxzz_xzzz[i] * fe_0 + ts_xxyzz_xzzz[i] * pa_y[i];

        ts_xxyyzz_yyyy[i] = ts_yyzz_yyyy[i] * fe_0 + ts_xyyzz_yyyy[i] * pa_x[i];

        ts_xxyyzz_yyyz[i] = ts_yyzz_yyyz[i] * fe_0 + ts_xyyzz_yyyz[i] * pa_x[i];

        ts_xxyyzz_yyzz[i] = ts_yyzz_yyzz[i] * fe_0 + ts_xyyzz_yyzz[i] * pa_x[i];

        ts_xxyyzz_yzzz[i] = ts_yyzz_yzzz[i] * fe_0 + ts_xyyzz_yzzz[i] * pa_x[i];

        ts_xxyyzz_zzzz[i] = ts_yyzz_zzzz[i] * fe_0 + ts_xyyzz_zzzz[i] * pa_x[i];
    }

    // Set up 195-210 components of targeted buffer : IG

    auto ts_xxyzzz_xxxx = pbuffer.data(idx_ovl_ig + 195);

    auto ts_xxyzzz_xxxy = pbuffer.data(idx_ovl_ig + 196);

    auto ts_xxyzzz_xxxz = pbuffer.data(idx_ovl_ig + 197);

    auto ts_xxyzzz_xxyy = pbuffer.data(idx_ovl_ig + 198);

    auto ts_xxyzzz_xxyz = pbuffer.data(idx_ovl_ig + 199);

    auto ts_xxyzzz_xxzz = pbuffer.data(idx_ovl_ig + 200);

    auto ts_xxyzzz_xyyy = pbuffer.data(idx_ovl_ig + 201);

    auto ts_xxyzzz_xyyz = pbuffer.data(idx_ovl_ig + 202);

    auto ts_xxyzzz_xyzz = pbuffer.data(idx_ovl_ig + 203);

    auto ts_xxyzzz_xzzz = pbuffer.data(idx_ovl_ig + 204);

    auto ts_xxyzzz_yyyy = pbuffer.data(idx_ovl_ig + 205);

    auto ts_xxyzzz_yyyz = pbuffer.data(idx_ovl_ig + 206);

    auto ts_xxyzzz_yyzz = pbuffer.data(idx_ovl_ig + 207);

    auto ts_xxyzzz_yzzz = pbuffer.data(idx_ovl_ig + 208);

    auto ts_xxyzzz_zzzz = pbuffer.data(idx_ovl_ig + 209);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyzzz_xxxx, ts_xxyzzz_xxxy, ts_xxyzzz_xxxz, ts_xxyzzz_xxyy, ts_xxyzzz_xxyz, ts_xxyzzz_xxzz, ts_xxyzzz_xyyy, ts_xxyzzz_xyyz, ts_xxyzzz_xyzz, ts_xxyzzz_xzzz, ts_xxyzzz_yyyy, ts_xxyzzz_yyyz, ts_xxyzzz_yyzz, ts_xxyzzz_yzzz, ts_xxyzzz_zzzz, ts_xxzzz_xxx, ts_xxzzz_xxxx, ts_xxzzz_xxxy, ts_xxzzz_xxxz, ts_xxzzz_xxy, ts_xxzzz_xxyy, ts_xxzzz_xxyz, ts_xxzzz_xxz, ts_xxzzz_xxzz, ts_xxzzz_xyy, ts_xxzzz_xyyy, ts_xxzzz_xyyz, ts_xxzzz_xyz, ts_xxzzz_xyzz, ts_xxzzz_xzz, ts_xxzzz_xzzz, ts_xxzzz_zzzz, ts_xyzzz_yyyy, ts_xyzzz_yyyz, ts_xyzzz_yyzz, ts_xyzzz_yzzz, ts_yzzz_yyyy, ts_yzzz_yyyz, ts_yzzz_yyzz, ts_yzzz_yzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyzzz_xxxx[i] = ts_xxzzz_xxxx[i] * pa_y[i];

        ts_xxyzzz_xxxy[i] = ts_xxzzz_xxx[i] * fe_0 + ts_xxzzz_xxxy[i] * pa_y[i];

        ts_xxyzzz_xxxz[i] = ts_xxzzz_xxxz[i] * pa_y[i];

        ts_xxyzzz_xxyy[i] = 2.0 * ts_xxzzz_xxy[i] * fe_0 + ts_xxzzz_xxyy[i] * pa_y[i];

        ts_xxyzzz_xxyz[i] = ts_xxzzz_xxz[i] * fe_0 + ts_xxzzz_xxyz[i] * pa_y[i];

        ts_xxyzzz_xxzz[i] = ts_xxzzz_xxzz[i] * pa_y[i];

        ts_xxyzzz_xyyy[i] = 3.0 * ts_xxzzz_xyy[i] * fe_0 + ts_xxzzz_xyyy[i] * pa_y[i];

        ts_xxyzzz_xyyz[i] = 2.0 * ts_xxzzz_xyz[i] * fe_0 + ts_xxzzz_xyyz[i] * pa_y[i];

        ts_xxyzzz_xyzz[i] = ts_xxzzz_xzz[i] * fe_0 + ts_xxzzz_xyzz[i] * pa_y[i];

        ts_xxyzzz_xzzz[i] = ts_xxzzz_xzzz[i] * pa_y[i];

        ts_xxyzzz_yyyy[i] = ts_yzzz_yyyy[i] * fe_0 + ts_xyzzz_yyyy[i] * pa_x[i];

        ts_xxyzzz_yyyz[i] = ts_yzzz_yyyz[i] * fe_0 + ts_xyzzz_yyyz[i] * pa_x[i];

        ts_xxyzzz_yyzz[i] = ts_yzzz_yyzz[i] * fe_0 + ts_xyzzz_yyzz[i] * pa_x[i];

        ts_xxyzzz_yzzz[i] = ts_yzzz_yzzz[i] * fe_0 + ts_xyzzz_yzzz[i] * pa_x[i];

        ts_xxyzzz_zzzz[i] = ts_xxzzz_zzzz[i] * pa_y[i];
    }

    // Set up 210-225 components of targeted buffer : IG

    auto ts_xxzzzz_xxxx = pbuffer.data(idx_ovl_ig + 210);

    auto ts_xxzzzz_xxxy = pbuffer.data(idx_ovl_ig + 211);

    auto ts_xxzzzz_xxxz = pbuffer.data(idx_ovl_ig + 212);

    auto ts_xxzzzz_xxyy = pbuffer.data(idx_ovl_ig + 213);

    auto ts_xxzzzz_xxyz = pbuffer.data(idx_ovl_ig + 214);

    auto ts_xxzzzz_xxzz = pbuffer.data(idx_ovl_ig + 215);

    auto ts_xxzzzz_xyyy = pbuffer.data(idx_ovl_ig + 216);

    auto ts_xxzzzz_xyyz = pbuffer.data(idx_ovl_ig + 217);

    auto ts_xxzzzz_xyzz = pbuffer.data(idx_ovl_ig + 218);

    auto ts_xxzzzz_xzzz = pbuffer.data(idx_ovl_ig + 219);

    auto ts_xxzzzz_yyyy = pbuffer.data(idx_ovl_ig + 220);

    auto ts_xxzzzz_yyyz = pbuffer.data(idx_ovl_ig + 221);

    auto ts_xxzzzz_yyzz = pbuffer.data(idx_ovl_ig + 222);

    auto ts_xxzzzz_yzzz = pbuffer.data(idx_ovl_ig + 223);

    auto ts_xxzzzz_zzzz = pbuffer.data(idx_ovl_ig + 224);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxzz_xxxx, ts_xxzz_xxxy, ts_xxzz_xxyy, ts_xxzz_xyyy, ts_xxzzz_xxxx, ts_xxzzz_xxxy, ts_xxzzz_xxyy, ts_xxzzz_xyyy, ts_xxzzzz_xxxx, ts_xxzzzz_xxxy, ts_xxzzzz_xxxz, ts_xxzzzz_xxyy, ts_xxzzzz_xxyz, ts_xxzzzz_xxzz, ts_xxzzzz_xyyy, ts_xxzzzz_xyyz, ts_xxzzzz_xyzz, ts_xxzzzz_xzzz, ts_xxzzzz_yyyy, ts_xxzzzz_yyyz, ts_xxzzzz_yyzz, ts_xxzzzz_yzzz, ts_xxzzzz_zzzz, ts_xzzzz_xxxz, ts_xzzzz_xxyz, ts_xzzzz_xxz, ts_xzzzz_xxzz, ts_xzzzz_xyyz, ts_xzzzz_xyz, ts_xzzzz_xyzz, ts_xzzzz_xzz, ts_xzzzz_xzzz, ts_xzzzz_yyyy, ts_xzzzz_yyyz, ts_xzzzz_yyz, ts_xzzzz_yyzz, ts_xzzzz_yzz, ts_xzzzz_yzzz, ts_xzzzz_zzz, ts_xzzzz_zzzz, ts_zzzz_xxxz, ts_zzzz_xxyz, ts_zzzz_xxzz, ts_zzzz_xyyz, ts_zzzz_xyzz, ts_zzzz_xzzz, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyzz, ts_zzzz_yzzz, ts_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzzzz_xxxx[i] = 3.0 * ts_xxzz_xxxx[i] * fe_0 + ts_xxzzz_xxxx[i] * pa_z[i];

        ts_xxzzzz_xxxy[i] = 3.0 * ts_xxzz_xxxy[i] * fe_0 + ts_xxzzz_xxxy[i] * pa_z[i];

        ts_xxzzzz_xxxz[i] = ts_zzzz_xxxz[i] * fe_0 + 3.0 * ts_xzzzz_xxz[i] * fe_0 + ts_xzzzz_xxxz[i] * pa_x[i];

        ts_xxzzzz_xxyy[i] = 3.0 * ts_xxzz_xxyy[i] * fe_0 + ts_xxzzz_xxyy[i] * pa_z[i];

        ts_xxzzzz_xxyz[i] = ts_zzzz_xxyz[i] * fe_0 + 2.0 * ts_xzzzz_xyz[i] * fe_0 + ts_xzzzz_xxyz[i] * pa_x[i];

        ts_xxzzzz_xxzz[i] = ts_zzzz_xxzz[i] * fe_0 + 2.0 * ts_xzzzz_xzz[i] * fe_0 + ts_xzzzz_xxzz[i] * pa_x[i];

        ts_xxzzzz_xyyy[i] = 3.0 * ts_xxzz_xyyy[i] * fe_0 + ts_xxzzz_xyyy[i] * pa_z[i];

        ts_xxzzzz_xyyz[i] = ts_zzzz_xyyz[i] * fe_0 + ts_xzzzz_yyz[i] * fe_0 + ts_xzzzz_xyyz[i] * pa_x[i];

        ts_xxzzzz_xyzz[i] = ts_zzzz_xyzz[i] * fe_0 + ts_xzzzz_yzz[i] * fe_0 + ts_xzzzz_xyzz[i] * pa_x[i];

        ts_xxzzzz_xzzz[i] = ts_zzzz_xzzz[i] * fe_0 + ts_xzzzz_zzz[i] * fe_0 + ts_xzzzz_xzzz[i] * pa_x[i];

        ts_xxzzzz_yyyy[i] = ts_zzzz_yyyy[i] * fe_0 + ts_xzzzz_yyyy[i] * pa_x[i];

        ts_xxzzzz_yyyz[i] = ts_zzzz_yyyz[i] * fe_0 + ts_xzzzz_yyyz[i] * pa_x[i];

        ts_xxzzzz_yyzz[i] = ts_zzzz_yyzz[i] * fe_0 + ts_xzzzz_yyzz[i] * pa_x[i];

        ts_xxzzzz_yzzz[i] = ts_zzzz_yzzz[i] * fe_0 + ts_xzzzz_yzzz[i] * pa_x[i];

        ts_xxzzzz_zzzz[i] = ts_zzzz_zzzz[i] * fe_0 + ts_xzzzz_zzzz[i] * pa_x[i];
    }

    // Set up 225-240 components of targeted buffer : IG

    auto ts_xyyyyy_xxxx = pbuffer.data(idx_ovl_ig + 225);

    auto ts_xyyyyy_xxxy = pbuffer.data(idx_ovl_ig + 226);

    auto ts_xyyyyy_xxxz = pbuffer.data(idx_ovl_ig + 227);

    auto ts_xyyyyy_xxyy = pbuffer.data(idx_ovl_ig + 228);

    auto ts_xyyyyy_xxyz = pbuffer.data(idx_ovl_ig + 229);

    auto ts_xyyyyy_xxzz = pbuffer.data(idx_ovl_ig + 230);

    auto ts_xyyyyy_xyyy = pbuffer.data(idx_ovl_ig + 231);

    auto ts_xyyyyy_xyyz = pbuffer.data(idx_ovl_ig + 232);

    auto ts_xyyyyy_xyzz = pbuffer.data(idx_ovl_ig + 233);

    auto ts_xyyyyy_xzzz = pbuffer.data(idx_ovl_ig + 234);

    auto ts_xyyyyy_yyyy = pbuffer.data(idx_ovl_ig + 235);

    auto ts_xyyyyy_yyyz = pbuffer.data(idx_ovl_ig + 236);

    auto ts_xyyyyy_yyzz = pbuffer.data(idx_ovl_ig + 237);

    auto ts_xyyyyy_yzzz = pbuffer.data(idx_ovl_ig + 238);

    auto ts_xyyyyy_zzzz = pbuffer.data(idx_ovl_ig + 239);

    #pragma omp simd aligned(pa_x, ts_xyyyyy_xxxx, ts_xyyyyy_xxxy, ts_xyyyyy_xxxz, ts_xyyyyy_xxyy, ts_xyyyyy_xxyz, ts_xyyyyy_xxzz, ts_xyyyyy_xyyy, ts_xyyyyy_xyyz, ts_xyyyyy_xyzz, ts_xyyyyy_xzzz, ts_xyyyyy_yyyy, ts_xyyyyy_yyyz, ts_xyyyyy_yyzz, ts_xyyyyy_yzzz, ts_xyyyyy_zzzz, ts_yyyyy_xxx, ts_yyyyy_xxxx, ts_yyyyy_xxxy, ts_yyyyy_xxxz, ts_yyyyy_xxy, ts_yyyyy_xxyy, ts_yyyyy_xxyz, ts_yyyyy_xxz, ts_yyyyy_xxzz, ts_yyyyy_xyy, ts_yyyyy_xyyy, ts_yyyyy_xyyz, ts_yyyyy_xyz, ts_yyyyy_xyzz, ts_yyyyy_xzz, ts_yyyyy_xzzz, ts_yyyyy_yyy, ts_yyyyy_yyyy, ts_yyyyy_yyyz, ts_yyyyy_yyz, ts_yyyyy_yyzz, ts_yyyyy_yzz, ts_yyyyy_yzzz, ts_yyyyy_zzz, ts_yyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyyy_xxxx[i] = 4.0 * ts_yyyyy_xxx[i] * fe_0 + ts_yyyyy_xxxx[i] * pa_x[i];

        ts_xyyyyy_xxxy[i] = 3.0 * ts_yyyyy_xxy[i] * fe_0 + ts_yyyyy_xxxy[i] * pa_x[i];

        ts_xyyyyy_xxxz[i] = 3.0 * ts_yyyyy_xxz[i] * fe_0 + ts_yyyyy_xxxz[i] * pa_x[i];

        ts_xyyyyy_xxyy[i] = 2.0 * ts_yyyyy_xyy[i] * fe_0 + ts_yyyyy_xxyy[i] * pa_x[i];

        ts_xyyyyy_xxyz[i] = 2.0 * ts_yyyyy_xyz[i] * fe_0 + ts_yyyyy_xxyz[i] * pa_x[i];

        ts_xyyyyy_xxzz[i] = 2.0 * ts_yyyyy_xzz[i] * fe_0 + ts_yyyyy_xxzz[i] * pa_x[i];

        ts_xyyyyy_xyyy[i] = ts_yyyyy_yyy[i] * fe_0 + ts_yyyyy_xyyy[i] * pa_x[i];

        ts_xyyyyy_xyyz[i] = ts_yyyyy_yyz[i] * fe_0 + ts_yyyyy_xyyz[i] * pa_x[i];

        ts_xyyyyy_xyzz[i] = ts_yyyyy_yzz[i] * fe_0 + ts_yyyyy_xyzz[i] * pa_x[i];

        ts_xyyyyy_xzzz[i] = ts_yyyyy_zzz[i] * fe_0 + ts_yyyyy_xzzz[i] * pa_x[i];

        ts_xyyyyy_yyyy[i] = ts_yyyyy_yyyy[i] * pa_x[i];

        ts_xyyyyy_yyyz[i] = ts_yyyyy_yyyz[i] * pa_x[i];

        ts_xyyyyy_yyzz[i] = ts_yyyyy_yyzz[i] * pa_x[i];

        ts_xyyyyy_yzzz[i] = ts_yyyyy_yzzz[i] * pa_x[i];

        ts_xyyyyy_zzzz[i] = ts_yyyyy_zzzz[i] * pa_x[i];
    }

    // Set up 240-255 components of targeted buffer : IG

    auto ts_xyyyyz_xxxx = pbuffer.data(idx_ovl_ig + 240);

    auto ts_xyyyyz_xxxy = pbuffer.data(idx_ovl_ig + 241);

    auto ts_xyyyyz_xxxz = pbuffer.data(idx_ovl_ig + 242);

    auto ts_xyyyyz_xxyy = pbuffer.data(idx_ovl_ig + 243);

    auto ts_xyyyyz_xxyz = pbuffer.data(idx_ovl_ig + 244);

    auto ts_xyyyyz_xxzz = pbuffer.data(idx_ovl_ig + 245);

    auto ts_xyyyyz_xyyy = pbuffer.data(idx_ovl_ig + 246);

    auto ts_xyyyyz_xyyz = pbuffer.data(idx_ovl_ig + 247);

    auto ts_xyyyyz_xyzz = pbuffer.data(idx_ovl_ig + 248);

    auto ts_xyyyyz_xzzz = pbuffer.data(idx_ovl_ig + 249);

    auto ts_xyyyyz_yyyy = pbuffer.data(idx_ovl_ig + 250);

    auto ts_xyyyyz_yyyz = pbuffer.data(idx_ovl_ig + 251);

    auto ts_xyyyyz_yyzz = pbuffer.data(idx_ovl_ig + 252);

    auto ts_xyyyyz_yzzz = pbuffer.data(idx_ovl_ig + 253);

    auto ts_xyyyyz_zzzz = pbuffer.data(idx_ovl_ig + 254);

    #pragma omp simd aligned(pa_x, pa_z, ts_xyyyy_xxxx, ts_xyyyy_xxxy, ts_xyyyy_xxyy, ts_xyyyy_xyyy, ts_xyyyyz_xxxx, ts_xyyyyz_xxxy, ts_xyyyyz_xxxz, ts_xyyyyz_xxyy, ts_xyyyyz_xxyz, ts_xyyyyz_xxzz, ts_xyyyyz_xyyy, ts_xyyyyz_xyyz, ts_xyyyyz_xyzz, ts_xyyyyz_xzzz, ts_xyyyyz_yyyy, ts_xyyyyz_yyyz, ts_xyyyyz_yyzz, ts_xyyyyz_yzzz, ts_xyyyyz_zzzz, ts_yyyyz_xxxz, ts_yyyyz_xxyz, ts_yyyyz_xxz, ts_yyyyz_xxzz, ts_yyyyz_xyyz, ts_yyyyz_xyz, ts_yyyyz_xyzz, ts_yyyyz_xzz, ts_yyyyz_xzzz, ts_yyyyz_yyyy, ts_yyyyz_yyyz, ts_yyyyz_yyz, ts_yyyyz_yyzz, ts_yyyyz_yzz, ts_yyyyz_yzzz, ts_yyyyz_zzz, ts_yyyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyyz_xxxx[i] = ts_xyyyy_xxxx[i] * pa_z[i];

        ts_xyyyyz_xxxy[i] = ts_xyyyy_xxxy[i] * pa_z[i];

        ts_xyyyyz_xxxz[i] = 3.0 * ts_yyyyz_xxz[i] * fe_0 + ts_yyyyz_xxxz[i] * pa_x[i];

        ts_xyyyyz_xxyy[i] = ts_xyyyy_xxyy[i] * pa_z[i];

        ts_xyyyyz_xxyz[i] = 2.0 * ts_yyyyz_xyz[i] * fe_0 + ts_yyyyz_xxyz[i] * pa_x[i];

        ts_xyyyyz_xxzz[i] = 2.0 * ts_yyyyz_xzz[i] * fe_0 + ts_yyyyz_xxzz[i] * pa_x[i];

        ts_xyyyyz_xyyy[i] = ts_xyyyy_xyyy[i] * pa_z[i];

        ts_xyyyyz_xyyz[i] = ts_yyyyz_yyz[i] * fe_0 + ts_yyyyz_xyyz[i] * pa_x[i];

        ts_xyyyyz_xyzz[i] = ts_yyyyz_yzz[i] * fe_0 + ts_yyyyz_xyzz[i] * pa_x[i];

        ts_xyyyyz_xzzz[i] = ts_yyyyz_zzz[i] * fe_0 + ts_yyyyz_xzzz[i] * pa_x[i];

        ts_xyyyyz_yyyy[i] = ts_yyyyz_yyyy[i] * pa_x[i];

        ts_xyyyyz_yyyz[i] = ts_yyyyz_yyyz[i] * pa_x[i];

        ts_xyyyyz_yyzz[i] = ts_yyyyz_yyzz[i] * pa_x[i];

        ts_xyyyyz_yzzz[i] = ts_yyyyz_yzzz[i] * pa_x[i];

        ts_xyyyyz_zzzz[i] = ts_yyyyz_zzzz[i] * pa_x[i];
    }

    // Set up 255-270 components of targeted buffer : IG

    auto ts_xyyyzz_xxxx = pbuffer.data(idx_ovl_ig + 255);

    auto ts_xyyyzz_xxxy = pbuffer.data(idx_ovl_ig + 256);

    auto ts_xyyyzz_xxxz = pbuffer.data(idx_ovl_ig + 257);

    auto ts_xyyyzz_xxyy = pbuffer.data(idx_ovl_ig + 258);

    auto ts_xyyyzz_xxyz = pbuffer.data(idx_ovl_ig + 259);

    auto ts_xyyyzz_xxzz = pbuffer.data(idx_ovl_ig + 260);

    auto ts_xyyyzz_xyyy = pbuffer.data(idx_ovl_ig + 261);

    auto ts_xyyyzz_xyyz = pbuffer.data(idx_ovl_ig + 262);

    auto ts_xyyyzz_xyzz = pbuffer.data(idx_ovl_ig + 263);

    auto ts_xyyyzz_xzzz = pbuffer.data(idx_ovl_ig + 264);

    auto ts_xyyyzz_yyyy = pbuffer.data(idx_ovl_ig + 265);

    auto ts_xyyyzz_yyyz = pbuffer.data(idx_ovl_ig + 266);

    auto ts_xyyyzz_yyzz = pbuffer.data(idx_ovl_ig + 267);

    auto ts_xyyyzz_yzzz = pbuffer.data(idx_ovl_ig + 268);

    auto ts_xyyyzz_zzzz = pbuffer.data(idx_ovl_ig + 269);

    #pragma omp simd aligned(pa_x, ts_xyyyzz_xxxx, ts_xyyyzz_xxxy, ts_xyyyzz_xxxz, ts_xyyyzz_xxyy, ts_xyyyzz_xxyz, ts_xyyyzz_xxzz, ts_xyyyzz_xyyy, ts_xyyyzz_xyyz, ts_xyyyzz_xyzz, ts_xyyyzz_xzzz, ts_xyyyzz_yyyy, ts_xyyyzz_yyyz, ts_xyyyzz_yyzz, ts_xyyyzz_yzzz, ts_xyyyzz_zzzz, ts_yyyzz_xxx, ts_yyyzz_xxxx, ts_yyyzz_xxxy, ts_yyyzz_xxxz, ts_yyyzz_xxy, ts_yyyzz_xxyy, ts_yyyzz_xxyz, ts_yyyzz_xxz, ts_yyyzz_xxzz, ts_yyyzz_xyy, ts_yyyzz_xyyy, ts_yyyzz_xyyz, ts_yyyzz_xyz, ts_yyyzz_xyzz, ts_yyyzz_xzz, ts_yyyzz_xzzz, ts_yyyzz_yyy, ts_yyyzz_yyyy, ts_yyyzz_yyyz, ts_yyyzz_yyz, ts_yyyzz_yyzz, ts_yyyzz_yzz, ts_yyyzz_yzzz, ts_yyyzz_zzz, ts_yyyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyzz_xxxx[i] = 4.0 * ts_yyyzz_xxx[i] * fe_0 + ts_yyyzz_xxxx[i] * pa_x[i];

        ts_xyyyzz_xxxy[i] = 3.0 * ts_yyyzz_xxy[i] * fe_0 + ts_yyyzz_xxxy[i] * pa_x[i];

        ts_xyyyzz_xxxz[i] = 3.0 * ts_yyyzz_xxz[i] * fe_0 + ts_yyyzz_xxxz[i] * pa_x[i];

        ts_xyyyzz_xxyy[i] = 2.0 * ts_yyyzz_xyy[i] * fe_0 + ts_yyyzz_xxyy[i] * pa_x[i];

        ts_xyyyzz_xxyz[i] = 2.0 * ts_yyyzz_xyz[i] * fe_0 + ts_yyyzz_xxyz[i] * pa_x[i];

        ts_xyyyzz_xxzz[i] = 2.0 * ts_yyyzz_xzz[i] * fe_0 + ts_yyyzz_xxzz[i] * pa_x[i];

        ts_xyyyzz_xyyy[i] = ts_yyyzz_yyy[i] * fe_0 + ts_yyyzz_xyyy[i] * pa_x[i];

        ts_xyyyzz_xyyz[i] = ts_yyyzz_yyz[i] * fe_0 + ts_yyyzz_xyyz[i] * pa_x[i];

        ts_xyyyzz_xyzz[i] = ts_yyyzz_yzz[i] * fe_0 + ts_yyyzz_xyzz[i] * pa_x[i];

        ts_xyyyzz_xzzz[i] = ts_yyyzz_zzz[i] * fe_0 + ts_yyyzz_xzzz[i] * pa_x[i];

        ts_xyyyzz_yyyy[i] = ts_yyyzz_yyyy[i] * pa_x[i];

        ts_xyyyzz_yyyz[i] = ts_yyyzz_yyyz[i] * pa_x[i];

        ts_xyyyzz_yyzz[i] = ts_yyyzz_yyzz[i] * pa_x[i];

        ts_xyyyzz_yzzz[i] = ts_yyyzz_yzzz[i] * pa_x[i];

        ts_xyyyzz_zzzz[i] = ts_yyyzz_zzzz[i] * pa_x[i];
    }

    // Set up 270-285 components of targeted buffer : IG

    auto ts_xyyzzz_xxxx = pbuffer.data(idx_ovl_ig + 270);

    auto ts_xyyzzz_xxxy = pbuffer.data(idx_ovl_ig + 271);

    auto ts_xyyzzz_xxxz = pbuffer.data(idx_ovl_ig + 272);

    auto ts_xyyzzz_xxyy = pbuffer.data(idx_ovl_ig + 273);

    auto ts_xyyzzz_xxyz = pbuffer.data(idx_ovl_ig + 274);

    auto ts_xyyzzz_xxzz = pbuffer.data(idx_ovl_ig + 275);

    auto ts_xyyzzz_xyyy = pbuffer.data(idx_ovl_ig + 276);

    auto ts_xyyzzz_xyyz = pbuffer.data(idx_ovl_ig + 277);

    auto ts_xyyzzz_xyzz = pbuffer.data(idx_ovl_ig + 278);

    auto ts_xyyzzz_xzzz = pbuffer.data(idx_ovl_ig + 279);

    auto ts_xyyzzz_yyyy = pbuffer.data(idx_ovl_ig + 280);

    auto ts_xyyzzz_yyyz = pbuffer.data(idx_ovl_ig + 281);

    auto ts_xyyzzz_yyzz = pbuffer.data(idx_ovl_ig + 282);

    auto ts_xyyzzz_yzzz = pbuffer.data(idx_ovl_ig + 283);

    auto ts_xyyzzz_zzzz = pbuffer.data(idx_ovl_ig + 284);

    #pragma omp simd aligned(pa_x, ts_xyyzzz_xxxx, ts_xyyzzz_xxxy, ts_xyyzzz_xxxz, ts_xyyzzz_xxyy, ts_xyyzzz_xxyz, ts_xyyzzz_xxzz, ts_xyyzzz_xyyy, ts_xyyzzz_xyyz, ts_xyyzzz_xyzz, ts_xyyzzz_xzzz, ts_xyyzzz_yyyy, ts_xyyzzz_yyyz, ts_xyyzzz_yyzz, ts_xyyzzz_yzzz, ts_xyyzzz_zzzz, ts_yyzzz_xxx, ts_yyzzz_xxxx, ts_yyzzz_xxxy, ts_yyzzz_xxxz, ts_yyzzz_xxy, ts_yyzzz_xxyy, ts_yyzzz_xxyz, ts_yyzzz_xxz, ts_yyzzz_xxzz, ts_yyzzz_xyy, ts_yyzzz_xyyy, ts_yyzzz_xyyz, ts_yyzzz_xyz, ts_yyzzz_xyzz, ts_yyzzz_xzz, ts_yyzzz_xzzz, ts_yyzzz_yyy, ts_yyzzz_yyyy, ts_yyzzz_yyyz, ts_yyzzz_yyz, ts_yyzzz_yyzz, ts_yyzzz_yzz, ts_yyzzz_yzzz, ts_yyzzz_zzz, ts_yyzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyzzz_xxxx[i] = 4.0 * ts_yyzzz_xxx[i] * fe_0 + ts_yyzzz_xxxx[i] * pa_x[i];

        ts_xyyzzz_xxxy[i] = 3.0 * ts_yyzzz_xxy[i] * fe_0 + ts_yyzzz_xxxy[i] * pa_x[i];

        ts_xyyzzz_xxxz[i] = 3.0 * ts_yyzzz_xxz[i] * fe_0 + ts_yyzzz_xxxz[i] * pa_x[i];

        ts_xyyzzz_xxyy[i] = 2.0 * ts_yyzzz_xyy[i] * fe_0 + ts_yyzzz_xxyy[i] * pa_x[i];

        ts_xyyzzz_xxyz[i] = 2.0 * ts_yyzzz_xyz[i] * fe_0 + ts_yyzzz_xxyz[i] * pa_x[i];

        ts_xyyzzz_xxzz[i] = 2.0 * ts_yyzzz_xzz[i] * fe_0 + ts_yyzzz_xxzz[i] * pa_x[i];

        ts_xyyzzz_xyyy[i] = ts_yyzzz_yyy[i] * fe_0 + ts_yyzzz_xyyy[i] * pa_x[i];

        ts_xyyzzz_xyyz[i] = ts_yyzzz_yyz[i] * fe_0 + ts_yyzzz_xyyz[i] * pa_x[i];

        ts_xyyzzz_xyzz[i] = ts_yyzzz_yzz[i] * fe_0 + ts_yyzzz_xyzz[i] * pa_x[i];

        ts_xyyzzz_xzzz[i] = ts_yyzzz_zzz[i] * fe_0 + ts_yyzzz_xzzz[i] * pa_x[i];

        ts_xyyzzz_yyyy[i] = ts_yyzzz_yyyy[i] * pa_x[i];

        ts_xyyzzz_yyyz[i] = ts_yyzzz_yyyz[i] * pa_x[i];

        ts_xyyzzz_yyzz[i] = ts_yyzzz_yyzz[i] * pa_x[i];

        ts_xyyzzz_yzzz[i] = ts_yyzzz_yzzz[i] * pa_x[i];

        ts_xyyzzz_zzzz[i] = ts_yyzzz_zzzz[i] * pa_x[i];
    }

    // Set up 285-300 components of targeted buffer : IG

    auto ts_xyzzzz_xxxx = pbuffer.data(idx_ovl_ig + 285);

    auto ts_xyzzzz_xxxy = pbuffer.data(idx_ovl_ig + 286);

    auto ts_xyzzzz_xxxz = pbuffer.data(idx_ovl_ig + 287);

    auto ts_xyzzzz_xxyy = pbuffer.data(idx_ovl_ig + 288);

    auto ts_xyzzzz_xxyz = pbuffer.data(idx_ovl_ig + 289);

    auto ts_xyzzzz_xxzz = pbuffer.data(idx_ovl_ig + 290);

    auto ts_xyzzzz_xyyy = pbuffer.data(idx_ovl_ig + 291);

    auto ts_xyzzzz_xyyz = pbuffer.data(idx_ovl_ig + 292);

    auto ts_xyzzzz_xyzz = pbuffer.data(idx_ovl_ig + 293);

    auto ts_xyzzzz_xzzz = pbuffer.data(idx_ovl_ig + 294);

    auto ts_xyzzzz_yyyy = pbuffer.data(idx_ovl_ig + 295);

    auto ts_xyzzzz_yyyz = pbuffer.data(idx_ovl_ig + 296);

    auto ts_xyzzzz_yyzz = pbuffer.data(idx_ovl_ig + 297);

    auto ts_xyzzzz_yzzz = pbuffer.data(idx_ovl_ig + 298);

    auto ts_xyzzzz_zzzz = pbuffer.data(idx_ovl_ig + 299);

    #pragma omp simd aligned(pa_x, pa_y, ts_xyzzzz_xxxx, ts_xyzzzz_xxxy, ts_xyzzzz_xxxz, ts_xyzzzz_xxyy, ts_xyzzzz_xxyz, ts_xyzzzz_xxzz, ts_xyzzzz_xyyy, ts_xyzzzz_xyyz, ts_xyzzzz_xyzz, ts_xyzzzz_xzzz, ts_xyzzzz_yyyy, ts_xyzzzz_yyyz, ts_xyzzzz_yyzz, ts_xyzzzz_yzzz, ts_xyzzzz_zzzz, ts_xzzzz_xxxx, ts_xzzzz_xxxz, ts_xzzzz_xxzz, ts_xzzzz_xzzz, ts_yzzzz_xxxy, ts_yzzzz_xxy, ts_yzzzz_xxyy, ts_yzzzz_xxyz, ts_yzzzz_xyy, ts_yzzzz_xyyy, ts_yzzzz_xyyz, ts_yzzzz_xyz, ts_yzzzz_xyzz, ts_yzzzz_yyy, ts_yzzzz_yyyy, ts_yzzzz_yyyz, ts_yzzzz_yyz, ts_yzzzz_yyzz, ts_yzzzz_yzz, ts_yzzzz_yzzz, ts_yzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyzzzz_xxxx[i] = ts_xzzzz_xxxx[i] * pa_y[i];

        ts_xyzzzz_xxxy[i] = 3.0 * ts_yzzzz_xxy[i] * fe_0 + ts_yzzzz_xxxy[i] * pa_x[i];

        ts_xyzzzz_xxxz[i] = ts_xzzzz_xxxz[i] * pa_y[i];

        ts_xyzzzz_xxyy[i] = 2.0 * ts_yzzzz_xyy[i] * fe_0 + ts_yzzzz_xxyy[i] * pa_x[i];

        ts_xyzzzz_xxyz[i] = 2.0 * ts_yzzzz_xyz[i] * fe_0 + ts_yzzzz_xxyz[i] * pa_x[i];

        ts_xyzzzz_xxzz[i] = ts_xzzzz_xxzz[i] * pa_y[i];

        ts_xyzzzz_xyyy[i] = ts_yzzzz_yyy[i] * fe_0 + ts_yzzzz_xyyy[i] * pa_x[i];

        ts_xyzzzz_xyyz[i] = ts_yzzzz_yyz[i] * fe_0 + ts_yzzzz_xyyz[i] * pa_x[i];

        ts_xyzzzz_xyzz[i] = ts_yzzzz_yzz[i] * fe_0 + ts_yzzzz_xyzz[i] * pa_x[i];

        ts_xyzzzz_xzzz[i] = ts_xzzzz_xzzz[i] * pa_y[i];

        ts_xyzzzz_yyyy[i] = ts_yzzzz_yyyy[i] * pa_x[i];

        ts_xyzzzz_yyyz[i] = ts_yzzzz_yyyz[i] * pa_x[i];

        ts_xyzzzz_yyzz[i] = ts_yzzzz_yyzz[i] * pa_x[i];

        ts_xyzzzz_yzzz[i] = ts_yzzzz_yzzz[i] * pa_x[i];

        ts_xyzzzz_zzzz[i] = ts_yzzzz_zzzz[i] * pa_x[i];
    }

    // Set up 300-315 components of targeted buffer : IG

    auto ts_xzzzzz_xxxx = pbuffer.data(idx_ovl_ig + 300);

    auto ts_xzzzzz_xxxy = pbuffer.data(idx_ovl_ig + 301);

    auto ts_xzzzzz_xxxz = pbuffer.data(idx_ovl_ig + 302);

    auto ts_xzzzzz_xxyy = pbuffer.data(idx_ovl_ig + 303);

    auto ts_xzzzzz_xxyz = pbuffer.data(idx_ovl_ig + 304);

    auto ts_xzzzzz_xxzz = pbuffer.data(idx_ovl_ig + 305);

    auto ts_xzzzzz_xyyy = pbuffer.data(idx_ovl_ig + 306);

    auto ts_xzzzzz_xyyz = pbuffer.data(idx_ovl_ig + 307);

    auto ts_xzzzzz_xyzz = pbuffer.data(idx_ovl_ig + 308);

    auto ts_xzzzzz_xzzz = pbuffer.data(idx_ovl_ig + 309);

    auto ts_xzzzzz_yyyy = pbuffer.data(idx_ovl_ig + 310);

    auto ts_xzzzzz_yyyz = pbuffer.data(idx_ovl_ig + 311);

    auto ts_xzzzzz_yyzz = pbuffer.data(idx_ovl_ig + 312);

    auto ts_xzzzzz_yzzz = pbuffer.data(idx_ovl_ig + 313);

    auto ts_xzzzzz_zzzz = pbuffer.data(idx_ovl_ig + 314);

    #pragma omp simd aligned(pa_x, ts_xzzzzz_xxxx, ts_xzzzzz_xxxy, ts_xzzzzz_xxxz, ts_xzzzzz_xxyy, ts_xzzzzz_xxyz, ts_xzzzzz_xxzz, ts_xzzzzz_xyyy, ts_xzzzzz_xyyz, ts_xzzzzz_xyzz, ts_xzzzzz_xzzz, ts_xzzzzz_yyyy, ts_xzzzzz_yyyz, ts_xzzzzz_yyzz, ts_xzzzzz_yzzz, ts_xzzzzz_zzzz, ts_zzzzz_xxx, ts_zzzzz_xxxx, ts_zzzzz_xxxy, ts_zzzzz_xxxz, ts_zzzzz_xxy, ts_zzzzz_xxyy, ts_zzzzz_xxyz, ts_zzzzz_xxz, ts_zzzzz_xxzz, ts_zzzzz_xyy, ts_zzzzz_xyyy, ts_zzzzz_xyyz, ts_zzzzz_xyz, ts_zzzzz_xyzz, ts_zzzzz_xzz, ts_zzzzz_xzzz, ts_zzzzz_yyy, ts_zzzzz_yyyy, ts_zzzzz_yyyz, ts_zzzzz_yyz, ts_zzzzz_yyzz, ts_zzzzz_yzz, ts_zzzzz_yzzz, ts_zzzzz_zzz, ts_zzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzzzz_xxxx[i] = 4.0 * ts_zzzzz_xxx[i] * fe_0 + ts_zzzzz_xxxx[i] * pa_x[i];

        ts_xzzzzz_xxxy[i] = 3.0 * ts_zzzzz_xxy[i] * fe_0 + ts_zzzzz_xxxy[i] * pa_x[i];

        ts_xzzzzz_xxxz[i] = 3.0 * ts_zzzzz_xxz[i] * fe_0 + ts_zzzzz_xxxz[i] * pa_x[i];

        ts_xzzzzz_xxyy[i] = 2.0 * ts_zzzzz_xyy[i] * fe_0 + ts_zzzzz_xxyy[i] * pa_x[i];

        ts_xzzzzz_xxyz[i] = 2.0 * ts_zzzzz_xyz[i] * fe_0 + ts_zzzzz_xxyz[i] * pa_x[i];

        ts_xzzzzz_xxzz[i] = 2.0 * ts_zzzzz_xzz[i] * fe_0 + ts_zzzzz_xxzz[i] * pa_x[i];

        ts_xzzzzz_xyyy[i] = ts_zzzzz_yyy[i] * fe_0 + ts_zzzzz_xyyy[i] * pa_x[i];

        ts_xzzzzz_xyyz[i] = ts_zzzzz_yyz[i] * fe_0 + ts_zzzzz_xyyz[i] * pa_x[i];

        ts_xzzzzz_xyzz[i] = ts_zzzzz_yzz[i] * fe_0 + ts_zzzzz_xyzz[i] * pa_x[i];

        ts_xzzzzz_xzzz[i] = ts_zzzzz_zzz[i] * fe_0 + ts_zzzzz_xzzz[i] * pa_x[i];

        ts_xzzzzz_yyyy[i] = ts_zzzzz_yyyy[i] * pa_x[i];

        ts_xzzzzz_yyyz[i] = ts_zzzzz_yyyz[i] * pa_x[i];

        ts_xzzzzz_yyzz[i] = ts_zzzzz_yyzz[i] * pa_x[i];

        ts_xzzzzz_yzzz[i] = ts_zzzzz_yzzz[i] * pa_x[i];

        ts_xzzzzz_zzzz[i] = ts_zzzzz_zzzz[i] * pa_x[i];
    }

    // Set up 315-330 components of targeted buffer : IG

    auto ts_yyyyyy_xxxx = pbuffer.data(idx_ovl_ig + 315);

    auto ts_yyyyyy_xxxy = pbuffer.data(idx_ovl_ig + 316);

    auto ts_yyyyyy_xxxz = pbuffer.data(idx_ovl_ig + 317);

    auto ts_yyyyyy_xxyy = pbuffer.data(idx_ovl_ig + 318);

    auto ts_yyyyyy_xxyz = pbuffer.data(idx_ovl_ig + 319);

    auto ts_yyyyyy_xxzz = pbuffer.data(idx_ovl_ig + 320);

    auto ts_yyyyyy_xyyy = pbuffer.data(idx_ovl_ig + 321);

    auto ts_yyyyyy_xyyz = pbuffer.data(idx_ovl_ig + 322);

    auto ts_yyyyyy_xyzz = pbuffer.data(idx_ovl_ig + 323);

    auto ts_yyyyyy_xzzz = pbuffer.data(idx_ovl_ig + 324);

    auto ts_yyyyyy_yyyy = pbuffer.data(idx_ovl_ig + 325);

    auto ts_yyyyyy_yyyz = pbuffer.data(idx_ovl_ig + 326);

    auto ts_yyyyyy_yyzz = pbuffer.data(idx_ovl_ig + 327);

    auto ts_yyyyyy_yzzz = pbuffer.data(idx_ovl_ig + 328);

    auto ts_yyyyyy_zzzz = pbuffer.data(idx_ovl_ig + 329);

    #pragma omp simd aligned(pa_y, ts_yyyy_xxxx, ts_yyyy_xxxy, ts_yyyy_xxxz, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xxzz, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyzz, ts_yyyy_xzzz, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyzz, ts_yyyy_yzzz, ts_yyyy_zzzz, ts_yyyyy_xxx, ts_yyyyy_xxxx, ts_yyyyy_xxxy, ts_yyyyy_xxxz, ts_yyyyy_xxy, ts_yyyyy_xxyy, ts_yyyyy_xxyz, ts_yyyyy_xxz, ts_yyyyy_xxzz, ts_yyyyy_xyy, ts_yyyyy_xyyy, ts_yyyyy_xyyz, ts_yyyyy_xyz, ts_yyyyy_xyzz, ts_yyyyy_xzz, ts_yyyyy_xzzz, ts_yyyyy_yyy, ts_yyyyy_yyyy, ts_yyyyy_yyyz, ts_yyyyy_yyz, ts_yyyyy_yyzz, ts_yyyyy_yzz, ts_yyyyy_yzzz, ts_yyyyy_zzz, ts_yyyyy_zzzz, ts_yyyyyy_xxxx, ts_yyyyyy_xxxy, ts_yyyyyy_xxxz, ts_yyyyyy_xxyy, ts_yyyyyy_xxyz, ts_yyyyyy_xxzz, ts_yyyyyy_xyyy, ts_yyyyyy_xyyz, ts_yyyyyy_xyzz, ts_yyyyyy_xzzz, ts_yyyyyy_yyyy, ts_yyyyyy_yyyz, ts_yyyyyy_yyzz, ts_yyyyyy_yzzz, ts_yyyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyyy_xxxx[i] = 5.0 * ts_yyyy_xxxx[i] * fe_0 + ts_yyyyy_xxxx[i] * pa_y[i];

        ts_yyyyyy_xxxy[i] = 5.0 * ts_yyyy_xxxy[i] * fe_0 + ts_yyyyy_xxx[i] * fe_0 + ts_yyyyy_xxxy[i] * pa_y[i];

        ts_yyyyyy_xxxz[i] = 5.0 * ts_yyyy_xxxz[i] * fe_0 + ts_yyyyy_xxxz[i] * pa_y[i];

        ts_yyyyyy_xxyy[i] = 5.0 * ts_yyyy_xxyy[i] * fe_0 + 2.0 * ts_yyyyy_xxy[i] * fe_0 + ts_yyyyy_xxyy[i] * pa_y[i];

        ts_yyyyyy_xxyz[i] = 5.0 * ts_yyyy_xxyz[i] * fe_0 + ts_yyyyy_xxz[i] * fe_0 + ts_yyyyy_xxyz[i] * pa_y[i];

        ts_yyyyyy_xxzz[i] = 5.0 * ts_yyyy_xxzz[i] * fe_0 + ts_yyyyy_xxzz[i] * pa_y[i];

        ts_yyyyyy_xyyy[i] = 5.0 * ts_yyyy_xyyy[i] * fe_0 + 3.0 * ts_yyyyy_xyy[i] * fe_0 + ts_yyyyy_xyyy[i] * pa_y[i];

        ts_yyyyyy_xyyz[i] = 5.0 * ts_yyyy_xyyz[i] * fe_0 + 2.0 * ts_yyyyy_xyz[i] * fe_0 + ts_yyyyy_xyyz[i] * pa_y[i];

        ts_yyyyyy_xyzz[i] = 5.0 * ts_yyyy_xyzz[i] * fe_0 + ts_yyyyy_xzz[i] * fe_0 + ts_yyyyy_xyzz[i] * pa_y[i];

        ts_yyyyyy_xzzz[i] = 5.0 * ts_yyyy_xzzz[i] * fe_0 + ts_yyyyy_xzzz[i] * pa_y[i];

        ts_yyyyyy_yyyy[i] = 5.0 * ts_yyyy_yyyy[i] * fe_0 + 4.0 * ts_yyyyy_yyy[i] * fe_0 + ts_yyyyy_yyyy[i] * pa_y[i];

        ts_yyyyyy_yyyz[i] = 5.0 * ts_yyyy_yyyz[i] * fe_0 + 3.0 * ts_yyyyy_yyz[i] * fe_0 + ts_yyyyy_yyyz[i] * pa_y[i];

        ts_yyyyyy_yyzz[i] = 5.0 * ts_yyyy_yyzz[i] * fe_0 + 2.0 * ts_yyyyy_yzz[i] * fe_0 + ts_yyyyy_yyzz[i] * pa_y[i];

        ts_yyyyyy_yzzz[i] = 5.0 * ts_yyyy_yzzz[i] * fe_0 + ts_yyyyy_zzz[i] * fe_0 + ts_yyyyy_yzzz[i] * pa_y[i];

        ts_yyyyyy_zzzz[i] = 5.0 * ts_yyyy_zzzz[i] * fe_0 + ts_yyyyy_zzzz[i] * pa_y[i];
    }

    // Set up 330-345 components of targeted buffer : IG

    auto ts_yyyyyz_xxxx = pbuffer.data(idx_ovl_ig + 330);

    auto ts_yyyyyz_xxxy = pbuffer.data(idx_ovl_ig + 331);

    auto ts_yyyyyz_xxxz = pbuffer.data(idx_ovl_ig + 332);

    auto ts_yyyyyz_xxyy = pbuffer.data(idx_ovl_ig + 333);

    auto ts_yyyyyz_xxyz = pbuffer.data(idx_ovl_ig + 334);

    auto ts_yyyyyz_xxzz = pbuffer.data(idx_ovl_ig + 335);

    auto ts_yyyyyz_xyyy = pbuffer.data(idx_ovl_ig + 336);

    auto ts_yyyyyz_xyyz = pbuffer.data(idx_ovl_ig + 337);

    auto ts_yyyyyz_xyzz = pbuffer.data(idx_ovl_ig + 338);

    auto ts_yyyyyz_xzzz = pbuffer.data(idx_ovl_ig + 339);

    auto ts_yyyyyz_yyyy = pbuffer.data(idx_ovl_ig + 340);

    auto ts_yyyyyz_yyyz = pbuffer.data(idx_ovl_ig + 341);

    auto ts_yyyyyz_yyzz = pbuffer.data(idx_ovl_ig + 342);

    auto ts_yyyyyz_yzzz = pbuffer.data(idx_ovl_ig + 343);

    auto ts_yyyyyz_zzzz = pbuffer.data(idx_ovl_ig + 344);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyyy_xxxx, ts_yyyyy_xxxy, ts_yyyyy_xxy, ts_yyyyy_xxyy, ts_yyyyy_xxyz, ts_yyyyy_xyy, ts_yyyyy_xyyy, ts_yyyyy_xyyz, ts_yyyyy_xyz, ts_yyyyy_xyzz, ts_yyyyy_yyy, ts_yyyyy_yyyy, ts_yyyyy_yyyz, ts_yyyyy_yyz, ts_yyyyy_yyzz, ts_yyyyy_yzz, ts_yyyyy_yzzz, ts_yyyyyz_xxxx, ts_yyyyyz_xxxy, ts_yyyyyz_xxxz, ts_yyyyyz_xxyy, ts_yyyyyz_xxyz, ts_yyyyyz_xxzz, ts_yyyyyz_xyyy, ts_yyyyyz_xyyz, ts_yyyyyz_xyzz, ts_yyyyyz_xzzz, ts_yyyyyz_yyyy, ts_yyyyyz_yyyz, ts_yyyyyz_yyzz, ts_yyyyyz_yzzz, ts_yyyyyz_zzzz, ts_yyyyz_xxxz, ts_yyyyz_xxzz, ts_yyyyz_xzzz, ts_yyyyz_zzzz, ts_yyyz_xxxz, ts_yyyz_xxzz, ts_yyyz_xzzz, ts_yyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyyz_xxxx[i] = ts_yyyyy_xxxx[i] * pa_z[i];

        ts_yyyyyz_xxxy[i] = ts_yyyyy_xxxy[i] * pa_z[i];

        ts_yyyyyz_xxxz[i] = 4.0 * ts_yyyz_xxxz[i] * fe_0 + ts_yyyyz_xxxz[i] * pa_y[i];

        ts_yyyyyz_xxyy[i] = ts_yyyyy_xxyy[i] * pa_z[i];

        ts_yyyyyz_xxyz[i] = ts_yyyyy_xxy[i] * fe_0 + ts_yyyyy_xxyz[i] * pa_z[i];

        ts_yyyyyz_xxzz[i] = 4.0 * ts_yyyz_xxzz[i] * fe_0 + ts_yyyyz_xxzz[i] * pa_y[i];

        ts_yyyyyz_xyyy[i] = ts_yyyyy_xyyy[i] * pa_z[i];

        ts_yyyyyz_xyyz[i] = ts_yyyyy_xyy[i] * fe_0 + ts_yyyyy_xyyz[i] * pa_z[i];

        ts_yyyyyz_xyzz[i] = 2.0 * ts_yyyyy_xyz[i] * fe_0 + ts_yyyyy_xyzz[i] * pa_z[i];

        ts_yyyyyz_xzzz[i] = 4.0 * ts_yyyz_xzzz[i] * fe_0 + ts_yyyyz_xzzz[i] * pa_y[i];

        ts_yyyyyz_yyyy[i] = ts_yyyyy_yyyy[i] * pa_z[i];

        ts_yyyyyz_yyyz[i] = ts_yyyyy_yyy[i] * fe_0 + ts_yyyyy_yyyz[i] * pa_z[i];

        ts_yyyyyz_yyzz[i] = 2.0 * ts_yyyyy_yyz[i] * fe_0 + ts_yyyyy_yyzz[i] * pa_z[i];

        ts_yyyyyz_yzzz[i] = 3.0 * ts_yyyyy_yzz[i] * fe_0 + ts_yyyyy_yzzz[i] * pa_z[i];

        ts_yyyyyz_zzzz[i] = 4.0 * ts_yyyz_zzzz[i] * fe_0 + ts_yyyyz_zzzz[i] * pa_y[i];
    }

    // Set up 345-360 components of targeted buffer : IG

    auto ts_yyyyzz_xxxx = pbuffer.data(idx_ovl_ig + 345);

    auto ts_yyyyzz_xxxy = pbuffer.data(idx_ovl_ig + 346);

    auto ts_yyyyzz_xxxz = pbuffer.data(idx_ovl_ig + 347);

    auto ts_yyyyzz_xxyy = pbuffer.data(idx_ovl_ig + 348);

    auto ts_yyyyzz_xxyz = pbuffer.data(idx_ovl_ig + 349);

    auto ts_yyyyzz_xxzz = pbuffer.data(idx_ovl_ig + 350);

    auto ts_yyyyzz_xyyy = pbuffer.data(idx_ovl_ig + 351);

    auto ts_yyyyzz_xyyz = pbuffer.data(idx_ovl_ig + 352);

    auto ts_yyyyzz_xyzz = pbuffer.data(idx_ovl_ig + 353);

    auto ts_yyyyzz_xzzz = pbuffer.data(idx_ovl_ig + 354);

    auto ts_yyyyzz_yyyy = pbuffer.data(idx_ovl_ig + 355);

    auto ts_yyyyzz_yyyz = pbuffer.data(idx_ovl_ig + 356);

    auto ts_yyyyzz_yyzz = pbuffer.data(idx_ovl_ig + 357);

    auto ts_yyyyzz_yzzz = pbuffer.data(idx_ovl_ig + 358);

    auto ts_yyyyzz_zzzz = pbuffer.data(idx_ovl_ig + 359);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyy_xxxy, ts_yyyy_xxyy, ts_yyyy_xyyy, ts_yyyy_yyyy, ts_yyyyz_xxxy, ts_yyyyz_xxyy, ts_yyyyz_xyyy, ts_yyyyz_yyyy, ts_yyyyzz_xxxx, ts_yyyyzz_xxxy, ts_yyyyzz_xxxz, ts_yyyyzz_xxyy, ts_yyyyzz_xxyz, ts_yyyyzz_xxzz, ts_yyyyzz_xyyy, ts_yyyyzz_xyyz, ts_yyyyzz_xyzz, ts_yyyyzz_xzzz, ts_yyyyzz_yyyy, ts_yyyyzz_yyyz, ts_yyyyzz_yyzz, ts_yyyyzz_yzzz, ts_yyyyzz_zzzz, ts_yyyzz_xxxx, ts_yyyzz_xxxz, ts_yyyzz_xxyz, ts_yyyzz_xxz, ts_yyyzz_xxzz, ts_yyyzz_xyyz, ts_yyyzz_xyz, ts_yyyzz_xyzz, ts_yyyzz_xzz, ts_yyyzz_xzzz, ts_yyyzz_yyyz, ts_yyyzz_yyz, ts_yyyzz_yyzz, ts_yyyzz_yzz, ts_yyyzz_yzzz, ts_yyyzz_zzz, ts_yyyzz_zzzz, ts_yyzz_xxxx, ts_yyzz_xxxz, ts_yyzz_xxyz, ts_yyzz_xxzz, ts_yyzz_xyyz, ts_yyzz_xyzz, ts_yyzz_xzzz, ts_yyzz_yyyz, ts_yyzz_yyzz, ts_yyzz_yzzz, ts_yyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyzz_xxxx[i] = 3.0 * ts_yyzz_xxxx[i] * fe_0 + ts_yyyzz_xxxx[i] * pa_y[i];

        ts_yyyyzz_xxxy[i] = ts_yyyy_xxxy[i] * fe_0 + ts_yyyyz_xxxy[i] * pa_z[i];

        ts_yyyyzz_xxxz[i] = 3.0 * ts_yyzz_xxxz[i] * fe_0 + ts_yyyzz_xxxz[i] * pa_y[i];

        ts_yyyyzz_xxyy[i] = ts_yyyy_xxyy[i] * fe_0 + ts_yyyyz_xxyy[i] * pa_z[i];

        ts_yyyyzz_xxyz[i] = 3.0 * ts_yyzz_xxyz[i] * fe_0 + ts_yyyzz_xxz[i] * fe_0 + ts_yyyzz_xxyz[i] * pa_y[i];

        ts_yyyyzz_xxzz[i] = 3.0 * ts_yyzz_xxzz[i] * fe_0 + ts_yyyzz_xxzz[i] * pa_y[i];

        ts_yyyyzz_xyyy[i] = ts_yyyy_xyyy[i] * fe_0 + ts_yyyyz_xyyy[i] * pa_z[i];

        ts_yyyyzz_xyyz[i] = 3.0 * ts_yyzz_xyyz[i] * fe_0 + 2.0 * ts_yyyzz_xyz[i] * fe_0 + ts_yyyzz_xyyz[i] * pa_y[i];

        ts_yyyyzz_xyzz[i] = 3.0 * ts_yyzz_xyzz[i] * fe_0 + ts_yyyzz_xzz[i] * fe_0 + ts_yyyzz_xyzz[i] * pa_y[i];

        ts_yyyyzz_xzzz[i] = 3.0 * ts_yyzz_xzzz[i] * fe_0 + ts_yyyzz_xzzz[i] * pa_y[i];

        ts_yyyyzz_yyyy[i] = ts_yyyy_yyyy[i] * fe_0 + ts_yyyyz_yyyy[i] * pa_z[i];

        ts_yyyyzz_yyyz[i] = 3.0 * ts_yyzz_yyyz[i] * fe_0 + 3.0 * ts_yyyzz_yyz[i] * fe_0 + ts_yyyzz_yyyz[i] * pa_y[i];

        ts_yyyyzz_yyzz[i] = 3.0 * ts_yyzz_yyzz[i] * fe_0 + 2.0 * ts_yyyzz_yzz[i] * fe_0 + ts_yyyzz_yyzz[i] * pa_y[i];

        ts_yyyyzz_yzzz[i] = 3.0 * ts_yyzz_yzzz[i] * fe_0 + ts_yyyzz_zzz[i] * fe_0 + ts_yyyzz_yzzz[i] * pa_y[i];

        ts_yyyyzz_zzzz[i] = 3.0 * ts_yyzz_zzzz[i] * fe_0 + ts_yyyzz_zzzz[i] * pa_y[i];
    }

    // Set up 360-375 components of targeted buffer : IG

    auto ts_yyyzzz_xxxx = pbuffer.data(idx_ovl_ig + 360);

    auto ts_yyyzzz_xxxy = pbuffer.data(idx_ovl_ig + 361);

    auto ts_yyyzzz_xxxz = pbuffer.data(idx_ovl_ig + 362);

    auto ts_yyyzzz_xxyy = pbuffer.data(idx_ovl_ig + 363);

    auto ts_yyyzzz_xxyz = pbuffer.data(idx_ovl_ig + 364);

    auto ts_yyyzzz_xxzz = pbuffer.data(idx_ovl_ig + 365);

    auto ts_yyyzzz_xyyy = pbuffer.data(idx_ovl_ig + 366);

    auto ts_yyyzzz_xyyz = pbuffer.data(idx_ovl_ig + 367);

    auto ts_yyyzzz_xyzz = pbuffer.data(idx_ovl_ig + 368);

    auto ts_yyyzzz_xzzz = pbuffer.data(idx_ovl_ig + 369);

    auto ts_yyyzzz_yyyy = pbuffer.data(idx_ovl_ig + 370);

    auto ts_yyyzzz_yyyz = pbuffer.data(idx_ovl_ig + 371);

    auto ts_yyyzzz_yyzz = pbuffer.data(idx_ovl_ig + 372);

    auto ts_yyyzzz_yzzz = pbuffer.data(idx_ovl_ig + 373);

    auto ts_yyyzzz_zzzz = pbuffer.data(idx_ovl_ig + 374);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyz_xxxy, ts_yyyz_xxyy, ts_yyyz_xyyy, ts_yyyz_yyyy, ts_yyyzz_xxxy, ts_yyyzz_xxyy, ts_yyyzz_xyyy, ts_yyyzz_yyyy, ts_yyyzzz_xxxx, ts_yyyzzz_xxxy, ts_yyyzzz_xxxz, ts_yyyzzz_xxyy, ts_yyyzzz_xxyz, ts_yyyzzz_xxzz, ts_yyyzzz_xyyy, ts_yyyzzz_xyyz, ts_yyyzzz_xyzz, ts_yyyzzz_xzzz, ts_yyyzzz_yyyy, ts_yyyzzz_yyyz, ts_yyyzzz_yyzz, ts_yyyzzz_yzzz, ts_yyyzzz_zzzz, ts_yyzzz_xxxx, ts_yyzzz_xxxz, ts_yyzzz_xxyz, ts_yyzzz_xxz, ts_yyzzz_xxzz, ts_yyzzz_xyyz, ts_yyzzz_xyz, ts_yyzzz_xyzz, ts_yyzzz_xzz, ts_yyzzz_xzzz, ts_yyzzz_yyyz, ts_yyzzz_yyz, ts_yyzzz_yyzz, ts_yyzzz_yzz, ts_yyzzz_yzzz, ts_yyzzz_zzz, ts_yyzzz_zzzz, ts_yzzz_xxxx, ts_yzzz_xxxz, ts_yzzz_xxyz, ts_yzzz_xxzz, ts_yzzz_xyyz, ts_yzzz_xyzz, ts_yzzz_xzzz, ts_yzzz_yyyz, ts_yzzz_yyzz, ts_yzzz_yzzz, ts_yzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyzzz_xxxx[i] = 2.0 * ts_yzzz_xxxx[i] * fe_0 + ts_yyzzz_xxxx[i] * pa_y[i];

        ts_yyyzzz_xxxy[i] = 2.0 * ts_yyyz_xxxy[i] * fe_0 + ts_yyyzz_xxxy[i] * pa_z[i];

        ts_yyyzzz_xxxz[i] = 2.0 * ts_yzzz_xxxz[i] * fe_0 + ts_yyzzz_xxxz[i] * pa_y[i];

        ts_yyyzzz_xxyy[i] = 2.0 * ts_yyyz_xxyy[i] * fe_0 + ts_yyyzz_xxyy[i] * pa_z[i];

        ts_yyyzzz_xxyz[i] = 2.0 * ts_yzzz_xxyz[i] * fe_0 + ts_yyzzz_xxz[i] * fe_0 + ts_yyzzz_xxyz[i] * pa_y[i];

        ts_yyyzzz_xxzz[i] = 2.0 * ts_yzzz_xxzz[i] * fe_0 + ts_yyzzz_xxzz[i] * pa_y[i];

        ts_yyyzzz_xyyy[i] = 2.0 * ts_yyyz_xyyy[i] * fe_0 + ts_yyyzz_xyyy[i] * pa_z[i];

        ts_yyyzzz_xyyz[i] = 2.0 * ts_yzzz_xyyz[i] * fe_0 + 2.0 * ts_yyzzz_xyz[i] * fe_0 + ts_yyzzz_xyyz[i] * pa_y[i];

        ts_yyyzzz_xyzz[i] = 2.0 * ts_yzzz_xyzz[i] * fe_0 + ts_yyzzz_xzz[i] * fe_0 + ts_yyzzz_xyzz[i] * pa_y[i];

        ts_yyyzzz_xzzz[i] = 2.0 * ts_yzzz_xzzz[i] * fe_0 + ts_yyzzz_xzzz[i] * pa_y[i];

        ts_yyyzzz_yyyy[i] = 2.0 * ts_yyyz_yyyy[i] * fe_0 + ts_yyyzz_yyyy[i] * pa_z[i];

        ts_yyyzzz_yyyz[i] = 2.0 * ts_yzzz_yyyz[i] * fe_0 + 3.0 * ts_yyzzz_yyz[i] * fe_0 + ts_yyzzz_yyyz[i] * pa_y[i];

        ts_yyyzzz_yyzz[i] = 2.0 * ts_yzzz_yyzz[i] * fe_0 + 2.0 * ts_yyzzz_yzz[i] * fe_0 + ts_yyzzz_yyzz[i] * pa_y[i];

        ts_yyyzzz_yzzz[i] = 2.0 * ts_yzzz_yzzz[i] * fe_0 + ts_yyzzz_zzz[i] * fe_0 + ts_yyzzz_yzzz[i] * pa_y[i];

        ts_yyyzzz_zzzz[i] = 2.0 * ts_yzzz_zzzz[i] * fe_0 + ts_yyzzz_zzzz[i] * pa_y[i];
    }

    // Set up 375-390 components of targeted buffer : IG

    auto ts_yyzzzz_xxxx = pbuffer.data(idx_ovl_ig + 375);

    auto ts_yyzzzz_xxxy = pbuffer.data(idx_ovl_ig + 376);

    auto ts_yyzzzz_xxxz = pbuffer.data(idx_ovl_ig + 377);

    auto ts_yyzzzz_xxyy = pbuffer.data(idx_ovl_ig + 378);

    auto ts_yyzzzz_xxyz = pbuffer.data(idx_ovl_ig + 379);

    auto ts_yyzzzz_xxzz = pbuffer.data(idx_ovl_ig + 380);

    auto ts_yyzzzz_xyyy = pbuffer.data(idx_ovl_ig + 381);

    auto ts_yyzzzz_xyyz = pbuffer.data(idx_ovl_ig + 382);

    auto ts_yyzzzz_xyzz = pbuffer.data(idx_ovl_ig + 383);

    auto ts_yyzzzz_xzzz = pbuffer.data(idx_ovl_ig + 384);

    auto ts_yyzzzz_yyyy = pbuffer.data(idx_ovl_ig + 385);

    auto ts_yyzzzz_yyyz = pbuffer.data(idx_ovl_ig + 386);

    auto ts_yyzzzz_yyzz = pbuffer.data(idx_ovl_ig + 387);

    auto ts_yyzzzz_yzzz = pbuffer.data(idx_ovl_ig + 388);

    auto ts_yyzzzz_zzzz = pbuffer.data(idx_ovl_ig + 389);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyzz_xxxy, ts_yyzz_xxyy, ts_yyzz_xyyy, ts_yyzz_yyyy, ts_yyzzz_xxxy, ts_yyzzz_xxyy, ts_yyzzz_xyyy, ts_yyzzz_yyyy, ts_yyzzzz_xxxx, ts_yyzzzz_xxxy, ts_yyzzzz_xxxz, ts_yyzzzz_xxyy, ts_yyzzzz_xxyz, ts_yyzzzz_xxzz, ts_yyzzzz_xyyy, ts_yyzzzz_xyyz, ts_yyzzzz_xyzz, ts_yyzzzz_xzzz, ts_yyzzzz_yyyy, ts_yyzzzz_yyyz, ts_yyzzzz_yyzz, ts_yyzzzz_yzzz, ts_yyzzzz_zzzz, ts_yzzzz_xxxx, ts_yzzzz_xxxz, ts_yzzzz_xxyz, ts_yzzzz_xxz, ts_yzzzz_xxzz, ts_yzzzz_xyyz, ts_yzzzz_xyz, ts_yzzzz_xyzz, ts_yzzzz_xzz, ts_yzzzz_xzzz, ts_yzzzz_yyyz, ts_yzzzz_yyz, ts_yzzzz_yyzz, ts_yzzzz_yzz, ts_yzzzz_yzzz, ts_yzzzz_zzz, ts_yzzzz_zzzz, ts_zzzz_xxxx, ts_zzzz_xxxz, ts_zzzz_xxyz, ts_zzzz_xxzz, ts_zzzz_xyyz, ts_zzzz_xyzz, ts_zzzz_xzzz, ts_zzzz_yyyz, ts_zzzz_yyzz, ts_zzzz_yzzz, ts_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzzzz_xxxx[i] = ts_zzzz_xxxx[i] * fe_0 + ts_yzzzz_xxxx[i] * pa_y[i];

        ts_yyzzzz_xxxy[i] = 3.0 * ts_yyzz_xxxy[i] * fe_0 + ts_yyzzz_xxxy[i] * pa_z[i];

        ts_yyzzzz_xxxz[i] = ts_zzzz_xxxz[i] * fe_0 + ts_yzzzz_xxxz[i] * pa_y[i];

        ts_yyzzzz_xxyy[i] = 3.0 * ts_yyzz_xxyy[i] * fe_0 + ts_yyzzz_xxyy[i] * pa_z[i];

        ts_yyzzzz_xxyz[i] = ts_zzzz_xxyz[i] * fe_0 + ts_yzzzz_xxz[i] * fe_0 + ts_yzzzz_xxyz[i] * pa_y[i];

        ts_yyzzzz_xxzz[i] = ts_zzzz_xxzz[i] * fe_0 + ts_yzzzz_xxzz[i] * pa_y[i];

        ts_yyzzzz_xyyy[i] = 3.0 * ts_yyzz_xyyy[i] * fe_0 + ts_yyzzz_xyyy[i] * pa_z[i];

        ts_yyzzzz_xyyz[i] = ts_zzzz_xyyz[i] * fe_0 + 2.0 * ts_yzzzz_xyz[i] * fe_0 + ts_yzzzz_xyyz[i] * pa_y[i];

        ts_yyzzzz_xyzz[i] = ts_zzzz_xyzz[i] * fe_0 + ts_yzzzz_xzz[i] * fe_0 + ts_yzzzz_xyzz[i] * pa_y[i];

        ts_yyzzzz_xzzz[i] = ts_zzzz_xzzz[i] * fe_0 + ts_yzzzz_xzzz[i] * pa_y[i];

        ts_yyzzzz_yyyy[i] = 3.0 * ts_yyzz_yyyy[i] * fe_0 + ts_yyzzz_yyyy[i] * pa_z[i];

        ts_yyzzzz_yyyz[i] = ts_zzzz_yyyz[i] * fe_0 + 3.0 * ts_yzzzz_yyz[i] * fe_0 + ts_yzzzz_yyyz[i] * pa_y[i];

        ts_yyzzzz_yyzz[i] = ts_zzzz_yyzz[i] * fe_0 + 2.0 * ts_yzzzz_yzz[i] * fe_0 + ts_yzzzz_yyzz[i] * pa_y[i];

        ts_yyzzzz_yzzz[i] = ts_zzzz_yzzz[i] * fe_0 + ts_yzzzz_zzz[i] * fe_0 + ts_yzzzz_yzzz[i] * pa_y[i];

        ts_yyzzzz_zzzz[i] = ts_zzzz_zzzz[i] * fe_0 + ts_yzzzz_zzzz[i] * pa_y[i];
    }

    // Set up 390-405 components of targeted buffer : IG

    auto ts_yzzzzz_xxxx = pbuffer.data(idx_ovl_ig + 390);

    auto ts_yzzzzz_xxxy = pbuffer.data(idx_ovl_ig + 391);

    auto ts_yzzzzz_xxxz = pbuffer.data(idx_ovl_ig + 392);

    auto ts_yzzzzz_xxyy = pbuffer.data(idx_ovl_ig + 393);

    auto ts_yzzzzz_xxyz = pbuffer.data(idx_ovl_ig + 394);

    auto ts_yzzzzz_xxzz = pbuffer.data(idx_ovl_ig + 395);

    auto ts_yzzzzz_xyyy = pbuffer.data(idx_ovl_ig + 396);

    auto ts_yzzzzz_xyyz = pbuffer.data(idx_ovl_ig + 397);

    auto ts_yzzzzz_xyzz = pbuffer.data(idx_ovl_ig + 398);

    auto ts_yzzzzz_xzzz = pbuffer.data(idx_ovl_ig + 399);

    auto ts_yzzzzz_yyyy = pbuffer.data(idx_ovl_ig + 400);

    auto ts_yzzzzz_yyyz = pbuffer.data(idx_ovl_ig + 401);

    auto ts_yzzzzz_yyzz = pbuffer.data(idx_ovl_ig + 402);

    auto ts_yzzzzz_yzzz = pbuffer.data(idx_ovl_ig + 403);

    auto ts_yzzzzz_zzzz = pbuffer.data(idx_ovl_ig + 404);

    #pragma omp simd aligned(pa_y, ts_yzzzzz_xxxx, ts_yzzzzz_xxxy, ts_yzzzzz_xxxz, ts_yzzzzz_xxyy, ts_yzzzzz_xxyz, ts_yzzzzz_xxzz, ts_yzzzzz_xyyy, ts_yzzzzz_xyyz, ts_yzzzzz_xyzz, ts_yzzzzz_xzzz, ts_yzzzzz_yyyy, ts_yzzzzz_yyyz, ts_yzzzzz_yyzz, ts_yzzzzz_yzzz, ts_yzzzzz_zzzz, ts_zzzzz_xxx, ts_zzzzz_xxxx, ts_zzzzz_xxxy, ts_zzzzz_xxxz, ts_zzzzz_xxy, ts_zzzzz_xxyy, ts_zzzzz_xxyz, ts_zzzzz_xxz, ts_zzzzz_xxzz, ts_zzzzz_xyy, ts_zzzzz_xyyy, ts_zzzzz_xyyz, ts_zzzzz_xyz, ts_zzzzz_xyzz, ts_zzzzz_xzz, ts_zzzzz_xzzz, ts_zzzzz_yyy, ts_zzzzz_yyyy, ts_zzzzz_yyyz, ts_zzzzz_yyz, ts_zzzzz_yyzz, ts_zzzzz_yzz, ts_zzzzz_yzzz, ts_zzzzz_zzz, ts_zzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzzzz_xxxx[i] = ts_zzzzz_xxxx[i] * pa_y[i];

        ts_yzzzzz_xxxy[i] = ts_zzzzz_xxx[i] * fe_0 + ts_zzzzz_xxxy[i] * pa_y[i];

        ts_yzzzzz_xxxz[i] = ts_zzzzz_xxxz[i] * pa_y[i];

        ts_yzzzzz_xxyy[i] = 2.0 * ts_zzzzz_xxy[i] * fe_0 + ts_zzzzz_xxyy[i] * pa_y[i];

        ts_yzzzzz_xxyz[i] = ts_zzzzz_xxz[i] * fe_0 + ts_zzzzz_xxyz[i] * pa_y[i];

        ts_yzzzzz_xxzz[i] = ts_zzzzz_xxzz[i] * pa_y[i];

        ts_yzzzzz_xyyy[i] = 3.0 * ts_zzzzz_xyy[i] * fe_0 + ts_zzzzz_xyyy[i] * pa_y[i];

        ts_yzzzzz_xyyz[i] = 2.0 * ts_zzzzz_xyz[i] * fe_0 + ts_zzzzz_xyyz[i] * pa_y[i];

        ts_yzzzzz_xyzz[i] = ts_zzzzz_xzz[i] * fe_0 + ts_zzzzz_xyzz[i] * pa_y[i];

        ts_yzzzzz_xzzz[i] = ts_zzzzz_xzzz[i] * pa_y[i];

        ts_yzzzzz_yyyy[i] = 4.0 * ts_zzzzz_yyy[i] * fe_0 + ts_zzzzz_yyyy[i] * pa_y[i];

        ts_yzzzzz_yyyz[i] = 3.0 * ts_zzzzz_yyz[i] * fe_0 + ts_zzzzz_yyyz[i] * pa_y[i];

        ts_yzzzzz_yyzz[i] = 2.0 * ts_zzzzz_yzz[i] * fe_0 + ts_zzzzz_yyzz[i] * pa_y[i];

        ts_yzzzzz_yzzz[i] = ts_zzzzz_zzz[i] * fe_0 + ts_zzzzz_yzzz[i] * pa_y[i];

        ts_yzzzzz_zzzz[i] = ts_zzzzz_zzzz[i] * pa_y[i];
    }

    // Set up 405-420 components of targeted buffer : IG

    auto ts_zzzzzz_xxxx = pbuffer.data(idx_ovl_ig + 405);

    auto ts_zzzzzz_xxxy = pbuffer.data(idx_ovl_ig + 406);

    auto ts_zzzzzz_xxxz = pbuffer.data(idx_ovl_ig + 407);

    auto ts_zzzzzz_xxyy = pbuffer.data(idx_ovl_ig + 408);

    auto ts_zzzzzz_xxyz = pbuffer.data(idx_ovl_ig + 409);

    auto ts_zzzzzz_xxzz = pbuffer.data(idx_ovl_ig + 410);

    auto ts_zzzzzz_xyyy = pbuffer.data(idx_ovl_ig + 411);

    auto ts_zzzzzz_xyyz = pbuffer.data(idx_ovl_ig + 412);

    auto ts_zzzzzz_xyzz = pbuffer.data(idx_ovl_ig + 413);

    auto ts_zzzzzz_xzzz = pbuffer.data(idx_ovl_ig + 414);

    auto ts_zzzzzz_yyyy = pbuffer.data(idx_ovl_ig + 415);

    auto ts_zzzzzz_yyyz = pbuffer.data(idx_ovl_ig + 416);

    auto ts_zzzzzz_yyzz = pbuffer.data(idx_ovl_ig + 417);

    auto ts_zzzzzz_yzzz = pbuffer.data(idx_ovl_ig + 418);

    auto ts_zzzzzz_zzzz = pbuffer.data(idx_ovl_ig + 419);

    #pragma omp simd aligned(pa_z, ts_zzzz_xxxx, ts_zzzz_xxxy, ts_zzzz_xxxz, ts_zzzz_xxyy, ts_zzzz_xxyz, ts_zzzz_xxzz, ts_zzzz_xyyy, ts_zzzz_xyyz, ts_zzzz_xyzz, ts_zzzz_xzzz, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyzz, ts_zzzz_yzzz, ts_zzzz_zzzz, ts_zzzzz_xxx, ts_zzzzz_xxxx, ts_zzzzz_xxxy, ts_zzzzz_xxxz, ts_zzzzz_xxy, ts_zzzzz_xxyy, ts_zzzzz_xxyz, ts_zzzzz_xxz, ts_zzzzz_xxzz, ts_zzzzz_xyy, ts_zzzzz_xyyy, ts_zzzzz_xyyz, ts_zzzzz_xyz, ts_zzzzz_xyzz, ts_zzzzz_xzz, ts_zzzzz_xzzz, ts_zzzzz_yyy, ts_zzzzz_yyyy, ts_zzzzz_yyyz, ts_zzzzz_yyz, ts_zzzzz_yyzz, ts_zzzzz_yzz, ts_zzzzz_yzzz, ts_zzzzz_zzz, ts_zzzzz_zzzz, ts_zzzzzz_xxxx, ts_zzzzzz_xxxy, ts_zzzzzz_xxxz, ts_zzzzzz_xxyy, ts_zzzzzz_xxyz, ts_zzzzzz_xxzz, ts_zzzzzz_xyyy, ts_zzzzzz_xyyz, ts_zzzzzz_xyzz, ts_zzzzzz_xzzz, ts_zzzzzz_yyyy, ts_zzzzzz_yyyz, ts_zzzzzz_yyzz, ts_zzzzzz_yzzz, ts_zzzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzzzz_xxxx[i] = 5.0 * ts_zzzz_xxxx[i] * fe_0 + ts_zzzzz_xxxx[i] * pa_z[i];

        ts_zzzzzz_xxxy[i] = 5.0 * ts_zzzz_xxxy[i] * fe_0 + ts_zzzzz_xxxy[i] * pa_z[i];

        ts_zzzzzz_xxxz[i] = 5.0 * ts_zzzz_xxxz[i] * fe_0 + ts_zzzzz_xxx[i] * fe_0 + ts_zzzzz_xxxz[i] * pa_z[i];

        ts_zzzzzz_xxyy[i] = 5.0 * ts_zzzz_xxyy[i] * fe_0 + ts_zzzzz_xxyy[i] * pa_z[i];

        ts_zzzzzz_xxyz[i] = 5.0 * ts_zzzz_xxyz[i] * fe_0 + ts_zzzzz_xxy[i] * fe_0 + ts_zzzzz_xxyz[i] * pa_z[i];

        ts_zzzzzz_xxzz[i] = 5.0 * ts_zzzz_xxzz[i] * fe_0 + 2.0 * ts_zzzzz_xxz[i] * fe_0 + ts_zzzzz_xxzz[i] * pa_z[i];

        ts_zzzzzz_xyyy[i] = 5.0 * ts_zzzz_xyyy[i] * fe_0 + ts_zzzzz_xyyy[i] * pa_z[i];

        ts_zzzzzz_xyyz[i] = 5.0 * ts_zzzz_xyyz[i] * fe_0 + ts_zzzzz_xyy[i] * fe_0 + ts_zzzzz_xyyz[i] * pa_z[i];

        ts_zzzzzz_xyzz[i] = 5.0 * ts_zzzz_xyzz[i] * fe_0 + 2.0 * ts_zzzzz_xyz[i] * fe_0 + ts_zzzzz_xyzz[i] * pa_z[i];

        ts_zzzzzz_xzzz[i] = 5.0 * ts_zzzz_xzzz[i] * fe_0 + 3.0 * ts_zzzzz_xzz[i] * fe_0 + ts_zzzzz_xzzz[i] * pa_z[i];

        ts_zzzzzz_yyyy[i] = 5.0 * ts_zzzz_yyyy[i] * fe_0 + ts_zzzzz_yyyy[i] * pa_z[i];

        ts_zzzzzz_yyyz[i] = 5.0 * ts_zzzz_yyyz[i] * fe_0 + ts_zzzzz_yyy[i] * fe_0 + ts_zzzzz_yyyz[i] * pa_z[i];

        ts_zzzzzz_yyzz[i] = 5.0 * ts_zzzz_yyzz[i] * fe_0 + 2.0 * ts_zzzzz_yyz[i] * fe_0 + ts_zzzzz_yyzz[i] * pa_z[i];

        ts_zzzzzz_yzzz[i] = 5.0 * ts_zzzz_yzzz[i] * fe_0 + 3.0 * ts_zzzzz_yzz[i] * fe_0 + ts_zzzzz_yzzz[i] * pa_z[i];

        ts_zzzzzz_zzzz[i] = 5.0 * ts_zzzz_zzzz[i] * fe_0 + 4.0 * ts_zzzzz_zzz[i] * fe_0 + ts_zzzzz_zzzz[i] * pa_z[i];
    }

}

} // ovlrec namespace

