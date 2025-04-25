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

#include "OverlapPrimRecIH.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_ih(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_ih,
                     const size_t idx_ovl_gh,
                     const size_t idx_ovl_hg,
                     const size_t idx_ovl_hh,
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

    // Set up components of auxiliary buffer : GH

    auto ts_xxxx_xxxxx = pbuffer.data(idx_ovl_gh);

    auto ts_xxxx_xxxxy = pbuffer.data(idx_ovl_gh + 1);

    auto ts_xxxx_xxxxz = pbuffer.data(idx_ovl_gh + 2);

    auto ts_xxxx_xxxyy = pbuffer.data(idx_ovl_gh + 3);

    auto ts_xxxx_xxxyz = pbuffer.data(idx_ovl_gh + 4);

    auto ts_xxxx_xxxzz = pbuffer.data(idx_ovl_gh + 5);

    auto ts_xxxx_xxyyy = pbuffer.data(idx_ovl_gh + 6);

    auto ts_xxxx_xxyyz = pbuffer.data(idx_ovl_gh + 7);

    auto ts_xxxx_xxyzz = pbuffer.data(idx_ovl_gh + 8);

    auto ts_xxxx_xxzzz = pbuffer.data(idx_ovl_gh + 9);

    auto ts_xxxx_xyyyy = pbuffer.data(idx_ovl_gh + 10);

    auto ts_xxxx_xyyyz = pbuffer.data(idx_ovl_gh + 11);

    auto ts_xxxx_xyyzz = pbuffer.data(idx_ovl_gh + 12);

    auto ts_xxxx_xyzzz = pbuffer.data(idx_ovl_gh + 13);

    auto ts_xxxx_xzzzz = pbuffer.data(idx_ovl_gh + 14);

    auto ts_xxxx_yyyyy = pbuffer.data(idx_ovl_gh + 15);

    auto ts_xxxx_yyyyz = pbuffer.data(idx_ovl_gh + 16);

    auto ts_xxxx_yyyzz = pbuffer.data(idx_ovl_gh + 17);

    auto ts_xxxx_yyzzz = pbuffer.data(idx_ovl_gh + 18);

    auto ts_xxxx_yzzzz = pbuffer.data(idx_ovl_gh + 19);

    auto ts_xxxx_zzzzz = pbuffer.data(idx_ovl_gh + 20);

    auto ts_xxxy_xxxxx = pbuffer.data(idx_ovl_gh + 21);

    auto ts_xxxy_xxxxz = pbuffer.data(idx_ovl_gh + 23);

    auto ts_xxxy_xxxzz = pbuffer.data(idx_ovl_gh + 26);

    auto ts_xxxy_xxzzz = pbuffer.data(idx_ovl_gh + 30);

    auto ts_xxxy_xzzzz = pbuffer.data(idx_ovl_gh + 35);

    auto ts_xxxy_yyyyy = pbuffer.data(idx_ovl_gh + 36);

    auto ts_xxxy_yyyyz = pbuffer.data(idx_ovl_gh + 37);

    auto ts_xxxy_yyyzz = pbuffer.data(idx_ovl_gh + 38);

    auto ts_xxxy_yyzzz = pbuffer.data(idx_ovl_gh + 39);

    auto ts_xxxy_yzzzz = pbuffer.data(idx_ovl_gh + 40);

    auto ts_xxxz_xxxxx = pbuffer.data(idx_ovl_gh + 42);

    auto ts_xxxz_xxxxy = pbuffer.data(idx_ovl_gh + 43);

    auto ts_xxxz_xxxxz = pbuffer.data(idx_ovl_gh + 44);

    auto ts_xxxz_xxxyy = pbuffer.data(idx_ovl_gh + 45);

    auto ts_xxxz_xxxzz = pbuffer.data(idx_ovl_gh + 47);

    auto ts_xxxz_xxyyy = pbuffer.data(idx_ovl_gh + 48);

    auto ts_xxxz_xxzzz = pbuffer.data(idx_ovl_gh + 51);

    auto ts_xxxz_xyyyy = pbuffer.data(idx_ovl_gh + 52);

    auto ts_xxxz_xzzzz = pbuffer.data(idx_ovl_gh + 56);

    auto ts_xxxz_yyyyz = pbuffer.data(idx_ovl_gh + 58);

    auto ts_xxxz_yyyzz = pbuffer.data(idx_ovl_gh + 59);

    auto ts_xxxz_yyzzz = pbuffer.data(idx_ovl_gh + 60);

    auto ts_xxxz_yzzzz = pbuffer.data(idx_ovl_gh + 61);

    auto ts_xxxz_zzzzz = pbuffer.data(idx_ovl_gh + 62);

    auto ts_xxyy_xxxxx = pbuffer.data(idx_ovl_gh + 63);

    auto ts_xxyy_xxxxy = pbuffer.data(idx_ovl_gh + 64);

    auto ts_xxyy_xxxxz = pbuffer.data(idx_ovl_gh + 65);

    auto ts_xxyy_xxxyy = pbuffer.data(idx_ovl_gh + 66);

    auto ts_xxyy_xxxyz = pbuffer.data(idx_ovl_gh + 67);

    auto ts_xxyy_xxxzz = pbuffer.data(idx_ovl_gh + 68);

    auto ts_xxyy_xxyyy = pbuffer.data(idx_ovl_gh + 69);

    auto ts_xxyy_xxyyz = pbuffer.data(idx_ovl_gh + 70);

    auto ts_xxyy_xxyzz = pbuffer.data(idx_ovl_gh + 71);

    auto ts_xxyy_xxzzz = pbuffer.data(idx_ovl_gh + 72);

    auto ts_xxyy_xyyyy = pbuffer.data(idx_ovl_gh + 73);

    auto ts_xxyy_xyyyz = pbuffer.data(idx_ovl_gh + 74);

    auto ts_xxyy_xyyzz = pbuffer.data(idx_ovl_gh + 75);

    auto ts_xxyy_xyzzz = pbuffer.data(idx_ovl_gh + 76);

    auto ts_xxyy_xzzzz = pbuffer.data(idx_ovl_gh + 77);

    auto ts_xxyy_yyyyy = pbuffer.data(idx_ovl_gh + 78);

    auto ts_xxyy_yyyyz = pbuffer.data(idx_ovl_gh + 79);

    auto ts_xxyy_yyyzz = pbuffer.data(idx_ovl_gh + 80);

    auto ts_xxyy_yyzzz = pbuffer.data(idx_ovl_gh + 81);

    auto ts_xxyy_yzzzz = pbuffer.data(idx_ovl_gh + 82);

    auto ts_xxyy_zzzzz = pbuffer.data(idx_ovl_gh + 83);

    auto ts_xxyz_xxxxz = pbuffer.data(idx_ovl_gh + 86);

    auto ts_xxyz_xxxzz = pbuffer.data(idx_ovl_gh + 89);

    auto ts_xxyz_xxzzz = pbuffer.data(idx_ovl_gh + 93);

    auto ts_xxyz_xzzzz = pbuffer.data(idx_ovl_gh + 98);

    auto ts_xxyz_yyyyz = pbuffer.data(idx_ovl_gh + 100);

    auto ts_xxyz_yyyzz = pbuffer.data(idx_ovl_gh + 101);

    auto ts_xxyz_yyzzz = pbuffer.data(idx_ovl_gh + 102);

    auto ts_xxyz_yzzzz = pbuffer.data(idx_ovl_gh + 103);

    auto ts_xxzz_xxxxx = pbuffer.data(idx_ovl_gh + 105);

    auto ts_xxzz_xxxxy = pbuffer.data(idx_ovl_gh + 106);

    auto ts_xxzz_xxxxz = pbuffer.data(idx_ovl_gh + 107);

    auto ts_xxzz_xxxyy = pbuffer.data(idx_ovl_gh + 108);

    auto ts_xxzz_xxxyz = pbuffer.data(idx_ovl_gh + 109);

    auto ts_xxzz_xxxzz = pbuffer.data(idx_ovl_gh + 110);

    auto ts_xxzz_xxyyy = pbuffer.data(idx_ovl_gh + 111);

    auto ts_xxzz_xxyyz = pbuffer.data(idx_ovl_gh + 112);

    auto ts_xxzz_xxyzz = pbuffer.data(idx_ovl_gh + 113);

    auto ts_xxzz_xxzzz = pbuffer.data(idx_ovl_gh + 114);

    auto ts_xxzz_xyyyy = pbuffer.data(idx_ovl_gh + 115);

    auto ts_xxzz_xyyyz = pbuffer.data(idx_ovl_gh + 116);

    auto ts_xxzz_xyyzz = pbuffer.data(idx_ovl_gh + 117);

    auto ts_xxzz_xyzzz = pbuffer.data(idx_ovl_gh + 118);

    auto ts_xxzz_xzzzz = pbuffer.data(idx_ovl_gh + 119);

    auto ts_xxzz_yyyyy = pbuffer.data(idx_ovl_gh + 120);

    auto ts_xxzz_yyyyz = pbuffer.data(idx_ovl_gh + 121);

    auto ts_xxzz_yyyzz = pbuffer.data(idx_ovl_gh + 122);

    auto ts_xxzz_yyzzz = pbuffer.data(idx_ovl_gh + 123);

    auto ts_xxzz_yzzzz = pbuffer.data(idx_ovl_gh + 124);

    auto ts_xxzz_zzzzz = pbuffer.data(idx_ovl_gh + 125);

    auto ts_xyyy_xxxxy = pbuffer.data(idx_ovl_gh + 127);

    auto ts_xyyy_xxxyy = pbuffer.data(idx_ovl_gh + 129);

    auto ts_xyyy_xxxyz = pbuffer.data(idx_ovl_gh + 130);

    auto ts_xyyy_xxyyy = pbuffer.data(idx_ovl_gh + 132);

    auto ts_xyyy_xxyyz = pbuffer.data(idx_ovl_gh + 133);

    auto ts_xyyy_xxyzz = pbuffer.data(idx_ovl_gh + 134);

    auto ts_xyyy_xyyyy = pbuffer.data(idx_ovl_gh + 136);

    auto ts_xyyy_xyyyz = pbuffer.data(idx_ovl_gh + 137);

    auto ts_xyyy_xyyzz = pbuffer.data(idx_ovl_gh + 138);

    auto ts_xyyy_xyzzz = pbuffer.data(idx_ovl_gh + 139);

    auto ts_xyyy_yyyyy = pbuffer.data(idx_ovl_gh + 141);

    auto ts_xyyy_yyyyz = pbuffer.data(idx_ovl_gh + 142);

    auto ts_xyyy_yyyzz = pbuffer.data(idx_ovl_gh + 143);

    auto ts_xyyy_yyzzz = pbuffer.data(idx_ovl_gh + 144);

    auto ts_xyyy_yzzzz = pbuffer.data(idx_ovl_gh + 145);

    auto ts_xyyy_zzzzz = pbuffer.data(idx_ovl_gh + 146);

    auto ts_xyyz_yyyyz = pbuffer.data(idx_ovl_gh + 163);

    auto ts_xyyz_yyyzz = pbuffer.data(idx_ovl_gh + 164);

    auto ts_xyyz_yyzzz = pbuffer.data(idx_ovl_gh + 165);

    auto ts_xyyz_yzzzz = pbuffer.data(idx_ovl_gh + 166);

    auto ts_xyyz_zzzzz = pbuffer.data(idx_ovl_gh + 167);

    auto ts_xyzz_yyyyy = pbuffer.data(idx_ovl_gh + 183);

    auto ts_xyzz_yyyyz = pbuffer.data(idx_ovl_gh + 184);

    auto ts_xyzz_yyyzz = pbuffer.data(idx_ovl_gh + 185);

    auto ts_xyzz_yyzzz = pbuffer.data(idx_ovl_gh + 186);

    auto ts_xyzz_yzzzz = pbuffer.data(idx_ovl_gh + 187);

    auto ts_xzzz_xxxxz = pbuffer.data(idx_ovl_gh + 191);

    auto ts_xzzz_xxxyz = pbuffer.data(idx_ovl_gh + 193);

    auto ts_xzzz_xxxzz = pbuffer.data(idx_ovl_gh + 194);

    auto ts_xzzz_xxyyz = pbuffer.data(idx_ovl_gh + 196);

    auto ts_xzzz_xxyzz = pbuffer.data(idx_ovl_gh + 197);

    auto ts_xzzz_xxzzz = pbuffer.data(idx_ovl_gh + 198);

    auto ts_xzzz_xyyyz = pbuffer.data(idx_ovl_gh + 200);

    auto ts_xzzz_xyyzz = pbuffer.data(idx_ovl_gh + 201);

    auto ts_xzzz_xyzzz = pbuffer.data(idx_ovl_gh + 202);

    auto ts_xzzz_xzzzz = pbuffer.data(idx_ovl_gh + 203);

    auto ts_xzzz_yyyyy = pbuffer.data(idx_ovl_gh + 204);

    auto ts_xzzz_yyyyz = pbuffer.data(idx_ovl_gh + 205);

    auto ts_xzzz_yyyzz = pbuffer.data(idx_ovl_gh + 206);

    auto ts_xzzz_yyzzz = pbuffer.data(idx_ovl_gh + 207);

    auto ts_xzzz_yzzzz = pbuffer.data(idx_ovl_gh + 208);

    auto ts_xzzz_zzzzz = pbuffer.data(idx_ovl_gh + 209);

    auto ts_yyyy_xxxxx = pbuffer.data(idx_ovl_gh + 210);

    auto ts_yyyy_xxxxy = pbuffer.data(idx_ovl_gh + 211);

    auto ts_yyyy_xxxxz = pbuffer.data(idx_ovl_gh + 212);

    auto ts_yyyy_xxxyy = pbuffer.data(idx_ovl_gh + 213);

    auto ts_yyyy_xxxyz = pbuffer.data(idx_ovl_gh + 214);

    auto ts_yyyy_xxxzz = pbuffer.data(idx_ovl_gh + 215);

    auto ts_yyyy_xxyyy = pbuffer.data(idx_ovl_gh + 216);

    auto ts_yyyy_xxyyz = pbuffer.data(idx_ovl_gh + 217);

    auto ts_yyyy_xxyzz = pbuffer.data(idx_ovl_gh + 218);

    auto ts_yyyy_xxzzz = pbuffer.data(idx_ovl_gh + 219);

    auto ts_yyyy_xyyyy = pbuffer.data(idx_ovl_gh + 220);

    auto ts_yyyy_xyyyz = pbuffer.data(idx_ovl_gh + 221);

    auto ts_yyyy_xyyzz = pbuffer.data(idx_ovl_gh + 222);

    auto ts_yyyy_xyzzz = pbuffer.data(idx_ovl_gh + 223);

    auto ts_yyyy_xzzzz = pbuffer.data(idx_ovl_gh + 224);

    auto ts_yyyy_yyyyy = pbuffer.data(idx_ovl_gh + 225);

    auto ts_yyyy_yyyyz = pbuffer.data(idx_ovl_gh + 226);

    auto ts_yyyy_yyyzz = pbuffer.data(idx_ovl_gh + 227);

    auto ts_yyyy_yyzzz = pbuffer.data(idx_ovl_gh + 228);

    auto ts_yyyy_yzzzz = pbuffer.data(idx_ovl_gh + 229);

    auto ts_yyyy_zzzzz = pbuffer.data(idx_ovl_gh + 230);

    auto ts_yyyz_xxxxy = pbuffer.data(idx_ovl_gh + 232);

    auto ts_yyyz_xxxxz = pbuffer.data(idx_ovl_gh + 233);

    auto ts_yyyz_xxxyy = pbuffer.data(idx_ovl_gh + 234);

    auto ts_yyyz_xxxzz = pbuffer.data(idx_ovl_gh + 236);

    auto ts_yyyz_xxyyy = pbuffer.data(idx_ovl_gh + 237);

    auto ts_yyyz_xxzzz = pbuffer.data(idx_ovl_gh + 240);

    auto ts_yyyz_xyyyy = pbuffer.data(idx_ovl_gh + 241);

    auto ts_yyyz_xzzzz = pbuffer.data(idx_ovl_gh + 245);

    auto ts_yyyz_yyyyy = pbuffer.data(idx_ovl_gh + 246);

    auto ts_yyyz_yyyyz = pbuffer.data(idx_ovl_gh + 247);

    auto ts_yyyz_yyyzz = pbuffer.data(idx_ovl_gh + 248);

    auto ts_yyyz_yyzzz = pbuffer.data(idx_ovl_gh + 249);

    auto ts_yyyz_yzzzz = pbuffer.data(idx_ovl_gh + 250);

    auto ts_yyyz_zzzzz = pbuffer.data(idx_ovl_gh + 251);

    auto ts_yyzz_xxxxx = pbuffer.data(idx_ovl_gh + 252);

    auto ts_yyzz_xxxxy = pbuffer.data(idx_ovl_gh + 253);

    auto ts_yyzz_xxxxz = pbuffer.data(idx_ovl_gh + 254);

    auto ts_yyzz_xxxyy = pbuffer.data(idx_ovl_gh + 255);

    auto ts_yyzz_xxxyz = pbuffer.data(idx_ovl_gh + 256);

    auto ts_yyzz_xxxzz = pbuffer.data(idx_ovl_gh + 257);

    auto ts_yyzz_xxyyy = pbuffer.data(idx_ovl_gh + 258);

    auto ts_yyzz_xxyyz = pbuffer.data(idx_ovl_gh + 259);

    auto ts_yyzz_xxyzz = pbuffer.data(idx_ovl_gh + 260);

    auto ts_yyzz_xxzzz = pbuffer.data(idx_ovl_gh + 261);

    auto ts_yyzz_xyyyy = pbuffer.data(idx_ovl_gh + 262);

    auto ts_yyzz_xyyyz = pbuffer.data(idx_ovl_gh + 263);

    auto ts_yyzz_xyyzz = pbuffer.data(idx_ovl_gh + 264);

    auto ts_yyzz_xyzzz = pbuffer.data(idx_ovl_gh + 265);

    auto ts_yyzz_xzzzz = pbuffer.data(idx_ovl_gh + 266);

    auto ts_yyzz_yyyyy = pbuffer.data(idx_ovl_gh + 267);

    auto ts_yyzz_yyyyz = pbuffer.data(idx_ovl_gh + 268);

    auto ts_yyzz_yyyzz = pbuffer.data(idx_ovl_gh + 269);

    auto ts_yyzz_yyzzz = pbuffer.data(idx_ovl_gh + 270);

    auto ts_yyzz_yzzzz = pbuffer.data(idx_ovl_gh + 271);

    auto ts_yyzz_zzzzz = pbuffer.data(idx_ovl_gh + 272);

    auto ts_yzzz_xxxxx = pbuffer.data(idx_ovl_gh + 273);

    auto ts_yzzz_xxxxz = pbuffer.data(idx_ovl_gh + 275);

    auto ts_yzzz_xxxyz = pbuffer.data(idx_ovl_gh + 277);

    auto ts_yzzz_xxxzz = pbuffer.data(idx_ovl_gh + 278);

    auto ts_yzzz_xxyyz = pbuffer.data(idx_ovl_gh + 280);

    auto ts_yzzz_xxyzz = pbuffer.data(idx_ovl_gh + 281);

    auto ts_yzzz_xxzzz = pbuffer.data(idx_ovl_gh + 282);

    auto ts_yzzz_xyyyz = pbuffer.data(idx_ovl_gh + 284);

    auto ts_yzzz_xyyzz = pbuffer.data(idx_ovl_gh + 285);

    auto ts_yzzz_xyzzz = pbuffer.data(idx_ovl_gh + 286);

    auto ts_yzzz_xzzzz = pbuffer.data(idx_ovl_gh + 287);

    auto ts_yzzz_yyyyy = pbuffer.data(idx_ovl_gh + 288);

    auto ts_yzzz_yyyyz = pbuffer.data(idx_ovl_gh + 289);

    auto ts_yzzz_yyyzz = pbuffer.data(idx_ovl_gh + 290);

    auto ts_yzzz_yyzzz = pbuffer.data(idx_ovl_gh + 291);

    auto ts_yzzz_yzzzz = pbuffer.data(idx_ovl_gh + 292);

    auto ts_yzzz_zzzzz = pbuffer.data(idx_ovl_gh + 293);

    auto ts_zzzz_xxxxx = pbuffer.data(idx_ovl_gh + 294);

    auto ts_zzzz_xxxxy = pbuffer.data(idx_ovl_gh + 295);

    auto ts_zzzz_xxxxz = pbuffer.data(idx_ovl_gh + 296);

    auto ts_zzzz_xxxyy = pbuffer.data(idx_ovl_gh + 297);

    auto ts_zzzz_xxxyz = pbuffer.data(idx_ovl_gh + 298);

    auto ts_zzzz_xxxzz = pbuffer.data(idx_ovl_gh + 299);

    auto ts_zzzz_xxyyy = pbuffer.data(idx_ovl_gh + 300);

    auto ts_zzzz_xxyyz = pbuffer.data(idx_ovl_gh + 301);

    auto ts_zzzz_xxyzz = pbuffer.data(idx_ovl_gh + 302);

    auto ts_zzzz_xxzzz = pbuffer.data(idx_ovl_gh + 303);

    auto ts_zzzz_xyyyy = pbuffer.data(idx_ovl_gh + 304);

    auto ts_zzzz_xyyyz = pbuffer.data(idx_ovl_gh + 305);

    auto ts_zzzz_xyyzz = pbuffer.data(idx_ovl_gh + 306);

    auto ts_zzzz_xyzzz = pbuffer.data(idx_ovl_gh + 307);

    auto ts_zzzz_xzzzz = pbuffer.data(idx_ovl_gh + 308);

    auto ts_zzzz_yyyyy = pbuffer.data(idx_ovl_gh + 309);

    auto ts_zzzz_yyyyz = pbuffer.data(idx_ovl_gh + 310);

    auto ts_zzzz_yyyzz = pbuffer.data(idx_ovl_gh + 311);

    auto ts_zzzz_yyzzz = pbuffer.data(idx_ovl_gh + 312);

    auto ts_zzzz_yzzzz = pbuffer.data(idx_ovl_gh + 313);

    auto ts_zzzz_zzzzz = pbuffer.data(idx_ovl_gh + 314);

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

    auto ts_xxxxz_xxxz = pbuffer.data(idx_ovl_hg + 32);

    auto ts_xxxxz_xxyz = pbuffer.data(idx_ovl_hg + 34);

    auto ts_xxxxz_xxzz = pbuffer.data(idx_ovl_hg + 35);

    auto ts_xxxxz_xyyz = pbuffer.data(idx_ovl_hg + 37);

    auto ts_xxxxz_xyzz = pbuffer.data(idx_ovl_hg + 38);

    auto ts_xxxxz_xzzz = pbuffer.data(idx_ovl_hg + 39);

    auto ts_xxxyy_xxxy = pbuffer.data(idx_ovl_hg + 46);

    auto ts_xxxyy_xxyy = pbuffer.data(idx_ovl_hg + 48);

    auto ts_xxxyy_xxyz = pbuffer.data(idx_ovl_hg + 49);

    auto ts_xxxyy_xyyy = pbuffer.data(idx_ovl_hg + 51);

    auto ts_xxxyy_xyyz = pbuffer.data(idx_ovl_hg + 52);

    auto ts_xxxyy_xyzz = pbuffer.data(idx_ovl_hg + 53);

    auto ts_xxxyy_yyyy = pbuffer.data(idx_ovl_hg + 55);

    auto ts_xxxyy_yyyz = pbuffer.data(idx_ovl_hg + 56);

    auto ts_xxxyy_yyzz = pbuffer.data(idx_ovl_hg + 57);

    auto ts_xxxyy_yzzz = pbuffer.data(idx_ovl_hg + 58);

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

    auto ts_xxxzz_yyyz = pbuffer.data(idx_ovl_hg + 86);

    auto ts_xxxzz_yyzz = pbuffer.data(idx_ovl_hg + 87);

    auto ts_xxxzz_yzzz = pbuffer.data(idx_ovl_hg + 88);

    auto ts_xxxzz_zzzz = pbuffer.data(idx_ovl_hg + 89);

    auto ts_xxyyy_xxxy = pbuffer.data(idx_ovl_hg + 91);

    auto ts_xxyyy_xxyy = pbuffer.data(idx_ovl_hg + 93);

    auto ts_xxyyy_xxyz = pbuffer.data(idx_ovl_hg + 94);

    auto ts_xxyyy_xyyy = pbuffer.data(idx_ovl_hg + 96);

    auto ts_xxyyy_xyyz = pbuffer.data(idx_ovl_hg + 97);

    auto ts_xxyyy_xyzz = pbuffer.data(idx_ovl_hg + 98);

    auto ts_xxyyy_yyyy = pbuffer.data(idx_ovl_hg + 100);

    auto ts_xxyyy_yyyz = pbuffer.data(idx_ovl_hg + 101);

    auto ts_xxyyy_yyzz = pbuffer.data(idx_ovl_hg + 102);

    auto ts_xxyyy_yzzz = pbuffer.data(idx_ovl_hg + 103);

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

    auto ts_xxzzz_yyyz = pbuffer.data(idx_ovl_hg + 146);

    auto ts_xxzzz_yyzz = pbuffer.data(idx_ovl_hg + 147);

    auto ts_xxzzz_yzzz = pbuffer.data(idx_ovl_hg + 148);

    auto ts_xxzzz_zzzz = pbuffer.data(idx_ovl_hg + 149);

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

    auto ts_xyyzz_xxyz = pbuffer.data(idx_ovl_hg + 184);

    auto ts_xyyzz_xyyz = pbuffer.data(idx_ovl_hg + 187);

    auto ts_xyyzz_xyzz = pbuffer.data(idx_ovl_hg + 188);

    auto ts_xyyzz_yyyz = pbuffer.data(idx_ovl_hg + 191);

    auto ts_xyyzz_yyzz = pbuffer.data(idx_ovl_hg + 192);

    auto ts_xyyzz_yzzz = pbuffer.data(idx_ovl_hg + 193);

    auto ts_xzzzz_xxxz = pbuffer.data(idx_ovl_hg + 212);

    auto ts_xzzzz_xxyz = pbuffer.data(idx_ovl_hg + 214);

    auto ts_xzzzz_xxzz = pbuffer.data(idx_ovl_hg + 215);

    auto ts_xzzzz_xyyz = pbuffer.data(idx_ovl_hg + 217);

    auto ts_xzzzz_xyzz = pbuffer.data(idx_ovl_hg + 218);

    auto ts_xzzzz_xzzz = pbuffer.data(idx_ovl_hg + 219);

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

    auto ts_yyyyz_xxxz = pbuffer.data(idx_ovl_hg + 242);

    auto ts_yyyyz_xxyz = pbuffer.data(idx_ovl_hg + 244);

    auto ts_yyyyz_xxzz = pbuffer.data(idx_ovl_hg + 245);

    auto ts_yyyyz_xyyz = pbuffer.data(idx_ovl_hg + 247);

    auto ts_yyyyz_xyzz = pbuffer.data(idx_ovl_hg + 248);

    auto ts_yyyyz_xzzz = pbuffer.data(idx_ovl_hg + 249);

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

    // Set up components of auxiliary buffer : HH

    auto ts_xxxxx_xxxxx = pbuffer.data(idx_ovl_hh);

    auto ts_xxxxx_xxxxy = pbuffer.data(idx_ovl_hh + 1);

    auto ts_xxxxx_xxxxz = pbuffer.data(idx_ovl_hh + 2);

    auto ts_xxxxx_xxxyy = pbuffer.data(idx_ovl_hh + 3);

    auto ts_xxxxx_xxxyz = pbuffer.data(idx_ovl_hh + 4);

    auto ts_xxxxx_xxxzz = pbuffer.data(idx_ovl_hh + 5);

    auto ts_xxxxx_xxyyy = pbuffer.data(idx_ovl_hh + 6);

    auto ts_xxxxx_xxyyz = pbuffer.data(idx_ovl_hh + 7);

    auto ts_xxxxx_xxyzz = pbuffer.data(idx_ovl_hh + 8);

    auto ts_xxxxx_xxzzz = pbuffer.data(idx_ovl_hh + 9);

    auto ts_xxxxx_xyyyy = pbuffer.data(idx_ovl_hh + 10);

    auto ts_xxxxx_xyyyz = pbuffer.data(idx_ovl_hh + 11);

    auto ts_xxxxx_xyyzz = pbuffer.data(idx_ovl_hh + 12);

    auto ts_xxxxx_xyzzz = pbuffer.data(idx_ovl_hh + 13);

    auto ts_xxxxx_xzzzz = pbuffer.data(idx_ovl_hh + 14);

    auto ts_xxxxx_yyyyy = pbuffer.data(idx_ovl_hh + 15);

    auto ts_xxxxx_yyyyz = pbuffer.data(idx_ovl_hh + 16);

    auto ts_xxxxx_yyyzz = pbuffer.data(idx_ovl_hh + 17);

    auto ts_xxxxx_yyzzz = pbuffer.data(idx_ovl_hh + 18);

    auto ts_xxxxx_yzzzz = pbuffer.data(idx_ovl_hh + 19);

    auto ts_xxxxx_zzzzz = pbuffer.data(idx_ovl_hh + 20);

    auto ts_xxxxy_xxxxx = pbuffer.data(idx_ovl_hh + 21);

    auto ts_xxxxy_xxxxy = pbuffer.data(idx_ovl_hh + 22);

    auto ts_xxxxy_xxxxz = pbuffer.data(idx_ovl_hh + 23);

    auto ts_xxxxy_xxxyy = pbuffer.data(idx_ovl_hh + 24);

    auto ts_xxxxy_xxxzz = pbuffer.data(idx_ovl_hh + 26);

    auto ts_xxxxy_xxyyy = pbuffer.data(idx_ovl_hh + 27);

    auto ts_xxxxy_xxzzz = pbuffer.data(idx_ovl_hh + 30);

    auto ts_xxxxy_xyyyy = pbuffer.data(idx_ovl_hh + 31);

    auto ts_xxxxy_xzzzz = pbuffer.data(idx_ovl_hh + 35);

    auto ts_xxxxy_yyyyy = pbuffer.data(idx_ovl_hh + 36);

    auto ts_xxxxy_yyyyz = pbuffer.data(idx_ovl_hh + 37);

    auto ts_xxxxy_yyyzz = pbuffer.data(idx_ovl_hh + 38);

    auto ts_xxxxy_yyzzz = pbuffer.data(idx_ovl_hh + 39);

    auto ts_xxxxy_yzzzz = pbuffer.data(idx_ovl_hh + 40);

    auto ts_xxxxz_xxxxx = pbuffer.data(idx_ovl_hh + 42);

    auto ts_xxxxz_xxxxy = pbuffer.data(idx_ovl_hh + 43);

    auto ts_xxxxz_xxxxz = pbuffer.data(idx_ovl_hh + 44);

    auto ts_xxxxz_xxxyy = pbuffer.data(idx_ovl_hh + 45);

    auto ts_xxxxz_xxxyz = pbuffer.data(idx_ovl_hh + 46);

    auto ts_xxxxz_xxxzz = pbuffer.data(idx_ovl_hh + 47);

    auto ts_xxxxz_xxyyy = pbuffer.data(idx_ovl_hh + 48);

    auto ts_xxxxz_xxyyz = pbuffer.data(idx_ovl_hh + 49);

    auto ts_xxxxz_xxyzz = pbuffer.data(idx_ovl_hh + 50);

    auto ts_xxxxz_xxzzz = pbuffer.data(idx_ovl_hh + 51);

    auto ts_xxxxz_xyyyy = pbuffer.data(idx_ovl_hh + 52);

    auto ts_xxxxz_xyyyz = pbuffer.data(idx_ovl_hh + 53);

    auto ts_xxxxz_xyyzz = pbuffer.data(idx_ovl_hh + 54);

    auto ts_xxxxz_xyzzz = pbuffer.data(idx_ovl_hh + 55);

    auto ts_xxxxz_xzzzz = pbuffer.data(idx_ovl_hh + 56);

    auto ts_xxxxz_yyyyz = pbuffer.data(idx_ovl_hh + 58);

    auto ts_xxxxz_yyyzz = pbuffer.data(idx_ovl_hh + 59);

    auto ts_xxxxz_yyzzz = pbuffer.data(idx_ovl_hh + 60);

    auto ts_xxxxz_yzzzz = pbuffer.data(idx_ovl_hh + 61);

    auto ts_xxxxz_zzzzz = pbuffer.data(idx_ovl_hh + 62);

    auto ts_xxxyy_xxxxx = pbuffer.data(idx_ovl_hh + 63);

    auto ts_xxxyy_xxxxy = pbuffer.data(idx_ovl_hh + 64);

    auto ts_xxxyy_xxxxz = pbuffer.data(idx_ovl_hh + 65);

    auto ts_xxxyy_xxxyy = pbuffer.data(idx_ovl_hh + 66);

    auto ts_xxxyy_xxxyz = pbuffer.data(idx_ovl_hh + 67);

    auto ts_xxxyy_xxxzz = pbuffer.data(idx_ovl_hh + 68);

    auto ts_xxxyy_xxyyy = pbuffer.data(idx_ovl_hh + 69);

    auto ts_xxxyy_xxyyz = pbuffer.data(idx_ovl_hh + 70);

    auto ts_xxxyy_xxyzz = pbuffer.data(idx_ovl_hh + 71);

    auto ts_xxxyy_xxzzz = pbuffer.data(idx_ovl_hh + 72);

    auto ts_xxxyy_xyyyy = pbuffer.data(idx_ovl_hh + 73);

    auto ts_xxxyy_xyyyz = pbuffer.data(idx_ovl_hh + 74);

    auto ts_xxxyy_xyyzz = pbuffer.data(idx_ovl_hh + 75);

    auto ts_xxxyy_xyzzz = pbuffer.data(idx_ovl_hh + 76);

    auto ts_xxxyy_xzzzz = pbuffer.data(idx_ovl_hh + 77);

    auto ts_xxxyy_yyyyy = pbuffer.data(idx_ovl_hh + 78);

    auto ts_xxxyy_yyyyz = pbuffer.data(idx_ovl_hh + 79);

    auto ts_xxxyy_yyyzz = pbuffer.data(idx_ovl_hh + 80);

    auto ts_xxxyy_yyzzz = pbuffer.data(idx_ovl_hh + 81);

    auto ts_xxxyy_yzzzz = pbuffer.data(idx_ovl_hh + 82);

    auto ts_xxxyy_zzzzz = pbuffer.data(idx_ovl_hh + 83);

    auto ts_xxxyz_xxxxz = pbuffer.data(idx_ovl_hh + 86);

    auto ts_xxxyz_xxxzz = pbuffer.data(idx_ovl_hh + 89);

    auto ts_xxxyz_xxzzz = pbuffer.data(idx_ovl_hh + 93);

    auto ts_xxxyz_xzzzz = pbuffer.data(idx_ovl_hh + 98);

    auto ts_xxxyz_yyyyz = pbuffer.data(idx_ovl_hh + 100);

    auto ts_xxxyz_yyyzz = pbuffer.data(idx_ovl_hh + 101);

    auto ts_xxxyz_yyzzz = pbuffer.data(idx_ovl_hh + 102);

    auto ts_xxxyz_yzzzz = pbuffer.data(idx_ovl_hh + 103);

    auto ts_xxxzz_xxxxx = pbuffer.data(idx_ovl_hh + 105);

    auto ts_xxxzz_xxxxy = pbuffer.data(idx_ovl_hh + 106);

    auto ts_xxxzz_xxxxz = pbuffer.data(idx_ovl_hh + 107);

    auto ts_xxxzz_xxxyy = pbuffer.data(idx_ovl_hh + 108);

    auto ts_xxxzz_xxxyz = pbuffer.data(idx_ovl_hh + 109);

    auto ts_xxxzz_xxxzz = pbuffer.data(idx_ovl_hh + 110);

    auto ts_xxxzz_xxyyy = pbuffer.data(idx_ovl_hh + 111);

    auto ts_xxxzz_xxyyz = pbuffer.data(idx_ovl_hh + 112);

    auto ts_xxxzz_xxyzz = pbuffer.data(idx_ovl_hh + 113);

    auto ts_xxxzz_xxzzz = pbuffer.data(idx_ovl_hh + 114);

    auto ts_xxxzz_xyyyy = pbuffer.data(idx_ovl_hh + 115);

    auto ts_xxxzz_xyyyz = pbuffer.data(idx_ovl_hh + 116);

    auto ts_xxxzz_xyyzz = pbuffer.data(idx_ovl_hh + 117);

    auto ts_xxxzz_xyzzz = pbuffer.data(idx_ovl_hh + 118);

    auto ts_xxxzz_xzzzz = pbuffer.data(idx_ovl_hh + 119);

    auto ts_xxxzz_yyyyy = pbuffer.data(idx_ovl_hh + 120);

    auto ts_xxxzz_yyyyz = pbuffer.data(idx_ovl_hh + 121);

    auto ts_xxxzz_yyyzz = pbuffer.data(idx_ovl_hh + 122);

    auto ts_xxxzz_yyzzz = pbuffer.data(idx_ovl_hh + 123);

    auto ts_xxxzz_yzzzz = pbuffer.data(idx_ovl_hh + 124);

    auto ts_xxxzz_zzzzz = pbuffer.data(idx_ovl_hh + 125);

    auto ts_xxyyy_xxxxx = pbuffer.data(idx_ovl_hh + 126);

    auto ts_xxyyy_xxxxy = pbuffer.data(idx_ovl_hh + 127);

    auto ts_xxyyy_xxxxz = pbuffer.data(idx_ovl_hh + 128);

    auto ts_xxyyy_xxxyy = pbuffer.data(idx_ovl_hh + 129);

    auto ts_xxyyy_xxxyz = pbuffer.data(idx_ovl_hh + 130);

    auto ts_xxyyy_xxxzz = pbuffer.data(idx_ovl_hh + 131);

    auto ts_xxyyy_xxyyy = pbuffer.data(idx_ovl_hh + 132);

    auto ts_xxyyy_xxyyz = pbuffer.data(idx_ovl_hh + 133);

    auto ts_xxyyy_xxyzz = pbuffer.data(idx_ovl_hh + 134);

    auto ts_xxyyy_xxzzz = pbuffer.data(idx_ovl_hh + 135);

    auto ts_xxyyy_xyyyy = pbuffer.data(idx_ovl_hh + 136);

    auto ts_xxyyy_xyyyz = pbuffer.data(idx_ovl_hh + 137);

    auto ts_xxyyy_xyyzz = pbuffer.data(idx_ovl_hh + 138);

    auto ts_xxyyy_xyzzz = pbuffer.data(idx_ovl_hh + 139);

    auto ts_xxyyy_xzzzz = pbuffer.data(idx_ovl_hh + 140);

    auto ts_xxyyy_yyyyy = pbuffer.data(idx_ovl_hh + 141);

    auto ts_xxyyy_yyyyz = pbuffer.data(idx_ovl_hh + 142);

    auto ts_xxyyy_yyyzz = pbuffer.data(idx_ovl_hh + 143);

    auto ts_xxyyy_yyzzz = pbuffer.data(idx_ovl_hh + 144);

    auto ts_xxyyy_yzzzz = pbuffer.data(idx_ovl_hh + 145);

    auto ts_xxyyy_zzzzz = pbuffer.data(idx_ovl_hh + 146);

    auto ts_xxyyz_xxxxy = pbuffer.data(idx_ovl_hh + 148);

    auto ts_xxyyz_xxxxz = pbuffer.data(idx_ovl_hh + 149);

    auto ts_xxyyz_xxxyy = pbuffer.data(idx_ovl_hh + 150);

    auto ts_xxyyz_xxxzz = pbuffer.data(idx_ovl_hh + 152);

    auto ts_xxyyz_xxyyy = pbuffer.data(idx_ovl_hh + 153);

    auto ts_xxyyz_xxzzz = pbuffer.data(idx_ovl_hh + 156);

    auto ts_xxyyz_xyyyy = pbuffer.data(idx_ovl_hh + 157);

    auto ts_xxyyz_xzzzz = pbuffer.data(idx_ovl_hh + 161);

    auto ts_xxyyz_yyyyz = pbuffer.data(idx_ovl_hh + 163);

    auto ts_xxyyz_yyyzz = pbuffer.data(idx_ovl_hh + 164);

    auto ts_xxyyz_yyzzz = pbuffer.data(idx_ovl_hh + 165);

    auto ts_xxyyz_yzzzz = pbuffer.data(idx_ovl_hh + 166);

    auto ts_xxyyz_zzzzz = pbuffer.data(idx_ovl_hh + 167);

    auto ts_xxyzz_xxxxx = pbuffer.data(idx_ovl_hh + 168);

    auto ts_xxyzz_xxxxz = pbuffer.data(idx_ovl_hh + 170);

    auto ts_xxyzz_xxxzz = pbuffer.data(idx_ovl_hh + 173);

    auto ts_xxyzz_xxzzz = pbuffer.data(idx_ovl_hh + 177);

    auto ts_xxyzz_xzzzz = pbuffer.data(idx_ovl_hh + 182);

    auto ts_xxyzz_yyyyy = pbuffer.data(idx_ovl_hh + 183);

    auto ts_xxyzz_yyyyz = pbuffer.data(idx_ovl_hh + 184);

    auto ts_xxyzz_yyyzz = pbuffer.data(idx_ovl_hh + 185);

    auto ts_xxyzz_yyzzz = pbuffer.data(idx_ovl_hh + 186);

    auto ts_xxyzz_yzzzz = pbuffer.data(idx_ovl_hh + 187);

    auto ts_xxzzz_xxxxx = pbuffer.data(idx_ovl_hh + 189);

    auto ts_xxzzz_xxxxy = pbuffer.data(idx_ovl_hh + 190);

    auto ts_xxzzz_xxxxz = pbuffer.data(idx_ovl_hh + 191);

    auto ts_xxzzz_xxxyy = pbuffer.data(idx_ovl_hh + 192);

    auto ts_xxzzz_xxxyz = pbuffer.data(idx_ovl_hh + 193);

    auto ts_xxzzz_xxxzz = pbuffer.data(idx_ovl_hh + 194);

    auto ts_xxzzz_xxyyy = pbuffer.data(idx_ovl_hh + 195);

    auto ts_xxzzz_xxyyz = pbuffer.data(idx_ovl_hh + 196);

    auto ts_xxzzz_xxyzz = pbuffer.data(idx_ovl_hh + 197);

    auto ts_xxzzz_xxzzz = pbuffer.data(idx_ovl_hh + 198);

    auto ts_xxzzz_xyyyy = pbuffer.data(idx_ovl_hh + 199);

    auto ts_xxzzz_xyyyz = pbuffer.data(idx_ovl_hh + 200);

    auto ts_xxzzz_xyyzz = pbuffer.data(idx_ovl_hh + 201);

    auto ts_xxzzz_xyzzz = pbuffer.data(idx_ovl_hh + 202);

    auto ts_xxzzz_xzzzz = pbuffer.data(idx_ovl_hh + 203);

    auto ts_xxzzz_yyyyy = pbuffer.data(idx_ovl_hh + 204);

    auto ts_xxzzz_yyyyz = pbuffer.data(idx_ovl_hh + 205);

    auto ts_xxzzz_yyyzz = pbuffer.data(idx_ovl_hh + 206);

    auto ts_xxzzz_yyzzz = pbuffer.data(idx_ovl_hh + 207);

    auto ts_xxzzz_yzzzz = pbuffer.data(idx_ovl_hh + 208);

    auto ts_xxzzz_zzzzz = pbuffer.data(idx_ovl_hh + 209);

    auto ts_xyyyy_xxxxx = pbuffer.data(idx_ovl_hh + 210);

    auto ts_xyyyy_xxxxy = pbuffer.data(idx_ovl_hh + 211);

    auto ts_xyyyy_xxxyy = pbuffer.data(idx_ovl_hh + 213);

    auto ts_xyyyy_xxxyz = pbuffer.data(idx_ovl_hh + 214);

    auto ts_xyyyy_xxyyy = pbuffer.data(idx_ovl_hh + 216);

    auto ts_xyyyy_xxyyz = pbuffer.data(idx_ovl_hh + 217);

    auto ts_xyyyy_xxyzz = pbuffer.data(idx_ovl_hh + 218);

    auto ts_xyyyy_xyyyy = pbuffer.data(idx_ovl_hh + 220);

    auto ts_xyyyy_xyyyz = pbuffer.data(idx_ovl_hh + 221);

    auto ts_xyyyy_xyyzz = pbuffer.data(idx_ovl_hh + 222);

    auto ts_xyyyy_xyzzz = pbuffer.data(idx_ovl_hh + 223);

    auto ts_xyyyy_yyyyy = pbuffer.data(idx_ovl_hh + 225);

    auto ts_xyyyy_yyyyz = pbuffer.data(idx_ovl_hh + 226);

    auto ts_xyyyy_yyyzz = pbuffer.data(idx_ovl_hh + 227);

    auto ts_xyyyy_yyzzz = pbuffer.data(idx_ovl_hh + 228);

    auto ts_xyyyy_yzzzz = pbuffer.data(idx_ovl_hh + 229);

    auto ts_xyyyy_zzzzz = pbuffer.data(idx_ovl_hh + 230);

    auto ts_xyyyz_yyyyz = pbuffer.data(idx_ovl_hh + 247);

    auto ts_xyyyz_yyyzz = pbuffer.data(idx_ovl_hh + 248);

    auto ts_xyyyz_yyzzz = pbuffer.data(idx_ovl_hh + 249);

    auto ts_xyyyz_yzzzz = pbuffer.data(idx_ovl_hh + 250);

    auto ts_xyyyz_zzzzz = pbuffer.data(idx_ovl_hh + 251);

    auto ts_xyyzz_xxxyz = pbuffer.data(idx_ovl_hh + 256);

    auto ts_xyyzz_xxyyz = pbuffer.data(idx_ovl_hh + 259);

    auto ts_xyyzz_xxyzz = pbuffer.data(idx_ovl_hh + 260);

    auto ts_xyyzz_xyyyz = pbuffer.data(idx_ovl_hh + 263);

    auto ts_xyyzz_xyyzz = pbuffer.data(idx_ovl_hh + 264);

    auto ts_xyyzz_xyzzz = pbuffer.data(idx_ovl_hh + 265);

    auto ts_xyyzz_yyyyy = pbuffer.data(idx_ovl_hh + 267);

    auto ts_xyyzz_yyyyz = pbuffer.data(idx_ovl_hh + 268);

    auto ts_xyyzz_yyyzz = pbuffer.data(idx_ovl_hh + 269);

    auto ts_xyyzz_yyzzz = pbuffer.data(idx_ovl_hh + 270);

    auto ts_xyyzz_yzzzz = pbuffer.data(idx_ovl_hh + 271);

    auto ts_xyyzz_zzzzz = pbuffer.data(idx_ovl_hh + 272);

    auto ts_xyzzz_yyyyy = pbuffer.data(idx_ovl_hh + 288);

    auto ts_xyzzz_yyyyz = pbuffer.data(idx_ovl_hh + 289);

    auto ts_xyzzz_yyyzz = pbuffer.data(idx_ovl_hh + 290);

    auto ts_xyzzz_yyzzz = pbuffer.data(idx_ovl_hh + 291);

    auto ts_xyzzz_yzzzz = pbuffer.data(idx_ovl_hh + 292);

    auto ts_xzzzz_xxxxx = pbuffer.data(idx_ovl_hh + 294);

    auto ts_xzzzz_xxxxz = pbuffer.data(idx_ovl_hh + 296);

    auto ts_xzzzz_xxxyz = pbuffer.data(idx_ovl_hh + 298);

    auto ts_xzzzz_xxxzz = pbuffer.data(idx_ovl_hh + 299);

    auto ts_xzzzz_xxyyz = pbuffer.data(idx_ovl_hh + 301);

    auto ts_xzzzz_xxyzz = pbuffer.data(idx_ovl_hh + 302);

    auto ts_xzzzz_xxzzz = pbuffer.data(idx_ovl_hh + 303);

    auto ts_xzzzz_xyyyz = pbuffer.data(idx_ovl_hh + 305);

    auto ts_xzzzz_xyyzz = pbuffer.data(idx_ovl_hh + 306);

    auto ts_xzzzz_xyzzz = pbuffer.data(idx_ovl_hh + 307);

    auto ts_xzzzz_xzzzz = pbuffer.data(idx_ovl_hh + 308);

    auto ts_xzzzz_yyyyy = pbuffer.data(idx_ovl_hh + 309);

    auto ts_xzzzz_yyyyz = pbuffer.data(idx_ovl_hh + 310);

    auto ts_xzzzz_yyyzz = pbuffer.data(idx_ovl_hh + 311);

    auto ts_xzzzz_yyzzz = pbuffer.data(idx_ovl_hh + 312);

    auto ts_xzzzz_yzzzz = pbuffer.data(idx_ovl_hh + 313);

    auto ts_xzzzz_zzzzz = pbuffer.data(idx_ovl_hh + 314);

    auto ts_yyyyy_xxxxx = pbuffer.data(idx_ovl_hh + 315);

    auto ts_yyyyy_xxxxy = pbuffer.data(idx_ovl_hh + 316);

    auto ts_yyyyy_xxxxz = pbuffer.data(idx_ovl_hh + 317);

    auto ts_yyyyy_xxxyy = pbuffer.data(idx_ovl_hh + 318);

    auto ts_yyyyy_xxxyz = pbuffer.data(idx_ovl_hh + 319);

    auto ts_yyyyy_xxxzz = pbuffer.data(idx_ovl_hh + 320);

    auto ts_yyyyy_xxyyy = pbuffer.data(idx_ovl_hh + 321);

    auto ts_yyyyy_xxyyz = pbuffer.data(idx_ovl_hh + 322);

    auto ts_yyyyy_xxyzz = pbuffer.data(idx_ovl_hh + 323);

    auto ts_yyyyy_xxzzz = pbuffer.data(idx_ovl_hh + 324);

    auto ts_yyyyy_xyyyy = pbuffer.data(idx_ovl_hh + 325);

    auto ts_yyyyy_xyyyz = pbuffer.data(idx_ovl_hh + 326);

    auto ts_yyyyy_xyyzz = pbuffer.data(idx_ovl_hh + 327);

    auto ts_yyyyy_xyzzz = pbuffer.data(idx_ovl_hh + 328);

    auto ts_yyyyy_xzzzz = pbuffer.data(idx_ovl_hh + 329);

    auto ts_yyyyy_yyyyy = pbuffer.data(idx_ovl_hh + 330);

    auto ts_yyyyy_yyyyz = pbuffer.data(idx_ovl_hh + 331);

    auto ts_yyyyy_yyyzz = pbuffer.data(idx_ovl_hh + 332);

    auto ts_yyyyy_yyzzz = pbuffer.data(idx_ovl_hh + 333);

    auto ts_yyyyy_yzzzz = pbuffer.data(idx_ovl_hh + 334);

    auto ts_yyyyy_zzzzz = pbuffer.data(idx_ovl_hh + 335);

    auto ts_yyyyz_xxxxy = pbuffer.data(idx_ovl_hh + 337);

    auto ts_yyyyz_xxxxz = pbuffer.data(idx_ovl_hh + 338);

    auto ts_yyyyz_xxxyy = pbuffer.data(idx_ovl_hh + 339);

    auto ts_yyyyz_xxxyz = pbuffer.data(idx_ovl_hh + 340);

    auto ts_yyyyz_xxxzz = pbuffer.data(idx_ovl_hh + 341);

    auto ts_yyyyz_xxyyy = pbuffer.data(idx_ovl_hh + 342);

    auto ts_yyyyz_xxyyz = pbuffer.data(idx_ovl_hh + 343);

    auto ts_yyyyz_xxyzz = pbuffer.data(idx_ovl_hh + 344);

    auto ts_yyyyz_xxzzz = pbuffer.data(idx_ovl_hh + 345);

    auto ts_yyyyz_xyyyy = pbuffer.data(idx_ovl_hh + 346);

    auto ts_yyyyz_xyyyz = pbuffer.data(idx_ovl_hh + 347);

    auto ts_yyyyz_xyyzz = pbuffer.data(idx_ovl_hh + 348);

    auto ts_yyyyz_xyzzz = pbuffer.data(idx_ovl_hh + 349);

    auto ts_yyyyz_xzzzz = pbuffer.data(idx_ovl_hh + 350);

    auto ts_yyyyz_yyyyy = pbuffer.data(idx_ovl_hh + 351);

    auto ts_yyyyz_yyyyz = pbuffer.data(idx_ovl_hh + 352);

    auto ts_yyyyz_yyyzz = pbuffer.data(idx_ovl_hh + 353);

    auto ts_yyyyz_yyzzz = pbuffer.data(idx_ovl_hh + 354);

    auto ts_yyyyz_yzzzz = pbuffer.data(idx_ovl_hh + 355);

    auto ts_yyyyz_zzzzz = pbuffer.data(idx_ovl_hh + 356);

    auto ts_yyyzz_xxxxx = pbuffer.data(idx_ovl_hh + 357);

    auto ts_yyyzz_xxxxy = pbuffer.data(idx_ovl_hh + 358);

    auto ts_yyyzz_xxxxz = pbuffer.data(idx_ovl_hh + 359);

    auto ts_yyyzz_xxxyy = pbuffer.data(idx_ovl_hh + 360);

    auto ts_yyyzz_xxxyz = pbuffer.data(idx_ovl_hh + 361);

    auto ts_yyyzz_xxxzz = pbuffer.data(idx_ovl_hh + 362);

    auto ts_yyyzz_xxyyy = pbuffer.data(idx_ovl_hh + 363);

    auto ts_yyyzz_xxyyz = pbuffer.data(idx_ovl_hh + 364);

    auto ts_yyyzz_xxyzz = pbuffer.data(idx_ovl_hh + 365);

    auto ts_yyyzz_xxzzz = pbuffer.data(idx_ovl_hh + 366);

    auto ts_yyyzz_xyyyy = pbuffer.data(idx_ovl_hh + 367);

    auto ts_yyyzz_xyyyz = pbuffer.data(idx_ovl_hh + 368);

    auto ts_yyyzz_xyyzz = pbuffer.data(idx_ovl_hh + 369);

    auto ts_yyyzz_xyzzz = pbuffer.data(idx_ovl_hh + 370);

    auto ts_yyyzz_xzzzz = pbuffer.data(idx_ovl_hh + 371);

    auto ts_yyyzz_yyyyy = pbuffer.data(idx_ovl_hh + 372);

    auto ts_yyyzz_yyyyz = pbuffer.data(idx_ovl_hh + 373);

    auto ts_yyyzz_yyyzz = pbuffer.data(idx_ovl_hh + 374);

    auto ts_yyyzz_yyzzz = pbuffer.data(idx_ovl_hh + 375);

    auto ts_yyyzz_yzzzz = pbuffer.data(idx_ovl_hh + 376);

    auto ts_yyyzz_zzzzz = pbuffer.data(idx_ovl_hh + 377);

    auto ts_yyzzz_xxxxx = pbuffer.data(idx_ovl_hh + 378);

    auto ts_yyzzz_xxxxy = pbuffer.data(idx_ovl_hh + 379);

    auto ts_yyzzz_xxxxz = pbuffer.data(idx_ovl_hh + 380);

    auto ts_yyzzz_xxxyy = pbuffer.data(idx_ovl_hh + 381);

    auto ts_yyzzz_xxxyz = pbuffer.data(idx_ovl_hh + 382);

    auto ts_yyzzz_xxxzz = pbuffer.data(idx_ovl_hh + 383);

    auto ts_yyzzz_xxyyy = pbuffer.data(idx_ovl_hh + 384);

    auto ts_yyzzz_xxyyz = pbuffer.data(idx_ovl_hh + 385);

    auto ts_yyzzz_xxyzz = pbuffer.data(idx_ovl_hh + 386);

    auto ts_yyzzz_xxzzz = pbuffer.data(idx_ovl_hh + 387);

    auto ts_yyzzz_xyyyy = pbuffer.data(idx_ovl_hh + 388);

    auto ts_yyzzz_xyyyz = pbuffer.data(idx_ovl_hh + 389);

    auto ts_yyzzz_xyyzz = pbuffer.data(idx_ovl_hh + 390);

    auto ts_yyzzz_xyzzz = pbuffer.data(idx_ovl_hh + 391);

    auto ts_yyzzz_xzzzz = pbuffer.data(idx_ovl_hh + 392);

    auto ts_yyzzz_yyyyy = pbuffer.data(idx_ovl_hh + 393);

    auto ts_yyzzz_yyyyz = pbuffer.data(idx_ovl_hh + 394);

    auto ts_yyzzz_yyyzz = pbuffer.data(idx_ovl_hh + 395);

    auto ts_yyzzz_yyzzz = pbuffer.data(idx_ovl_hh + 396);

    auto ts_yyzzz_yzzzz = pbuffer.data(idx_ovl_hh + 397);

    auto ts_yyzzz_zzzzz = pbuffer.data(idx_ovl_hh + 398);

    auto ts_yzzzz_xxxxx = pbuffer.data(idx_ovl_hh + 399);

    auto ts_yzzzz_xxxxy = pbuffer.data(idx_ovl_hh + 400);

    auto ts_yzzzz_xxxxz = pbuffer.data(idx_ovl_hh + 401);

    auto ts_yzzzz_xxxyy = pbuffer.data(idx_ovl_hh + 402);

    auto ts_yzzzz_xxxyz = pbuffer.data(idx_ovl_hh + 403);

    auto ts_yzzzz_xxxzz = pbuffer.data(idx_ovl_hh + 404);

    auto ts_yzzzz_xxyyy = pbuffer.data(idx_ovl_hh + 405);

    auto ts_yzzzz_xxyyz = pbuffer.data(idx_ovl_hh + 406);

    auto ts_yzzzz_xxyzz = pbuffer.data(idx_ovl_hh + 407);

    auto ts_yzzzz_xxzzz = pbuffer.data(idx_ovl_hh + 408);

    auto ts_yzzzz_xyyyy = pbuffer.data(idx_ovl_hh + 409);

    auto ts_yzzzz_xyyyz = pbuffer.data(idx_ovl_hh + 410);

    auto ts_yzzzz_xyyzz = pbuffer.data(idx_ovl_hh + 411);

    auto ts_yzzzz_xyzzz = pbuffer.data(idx_ovl_hh + 412);

    auto ts_yzzzz_xzzzz = pbuffer.data(idx_ovl_hh + 413);

    auto ts_yzzzz_yyyyy = pbuffer.data(idx_ovl_hh + 414);

    auto ts_yzzzz_yyyyz = pbuffer.data(idx_ovl_hh + 415);

    auto ts_yzzzz_yyyzz = pbuffer.data(idx_ovl_hh + 416);

    auto ts_yzzzz_yyzzz = pbuffer.data(idx_ovl_hh + 417);

    auto ts_yzzzz_yzzzz = pbuffer.data(idx_ovl_hh + 418);

    auto ts_yzzzz_zzzzz = pbuffer.data(idx_ovl_hh + 419);

    auto ts_zzzzz_xxxxx = pbuffer.data(idx_ovl_hh + 420);

    auto ts_zzzzz_xxxxy = pbuffer.data(idx_ovl_hh + 421);

    auto ts_zzzzz_xxxxz = pbuffer.data(idx_ovl_hh + 422);

    auto ts_zzzzz_xxxyy = pbuffer.data(idx_ovl_hh + 423);

    auto ts_zzzzz_xxxyz = pbuffer.data(idx_ovl_hh + 424);

    auto ts_zzzzz_xxxzz = pbuffer.data(idx_ovl_hh + 425);

    auto ts_zzzzz_xxyyy = pbuffer.data(idx_ovl_hh + 426);

    auto ts_zzzzz_xxyyz = pbuffer.data(idx_ovl_hh + 427);

    auto ts_zzzzz_xxyzz = pbuffer.data(idx_ovl_hh + 428);

    auto ts_zzzzz_xxzzz = pbuffer.data(idx_ovl_hh + 429);

    auto ts_zzzzz_xyyyy = pbuffer.data(idx_ovl_hh + 430);

    auto ts_zzzzz_xyyyz = pbuffer.data(idx_ovl_hh + 431);

    auto ts_zzzzz_xyyzz = pbuffer.data(idx_ovl_hh + 432);

    auto ts_zzzzz_xyzzz = pbuffer.data(idx_ovl_hh + 433);

    auto ts_zzzzz_xzzzz = pbuffer.data(idx_ovl_hh + 434);

    auto ts_zzzzz_yyyyy = pbuffer.data(idx_ovl_hh + 435);

    auto ts_zzzzz_yyyyz = pbuffer.data(idx_ovl_hh + 436);

    auto ts_zzzzz_yyyzz = pbuffer.data(idx_ovl_hh + 437);

    auto ts_zzzzz_yyzzz = pbuffer.data(idx_ovl_hh + 438);

    auto ts_zzzzz_yzzzz = pbuffer.data(idx_ovl_hh + 439);

    auto ts_zzzzz_zzzzz = pbuffer.data(idx_ovl_hh + 440);

    // Set up 0-21 components of targeted buffer : IH

    auto ts_xxxxxx_xxxxx = pbuffer.data(idx_ovl_ih);

    auto ts_xxxxxx_xxxxy = pbuffer.data(idx_ovl_ih + 1);

    auto ts_xxxxxx_xxxxz = pbuffer.data(idx_ovl_ih + 2);

    auto ts_xxxxxx_xxxyy = pbuffer.data(idx_ovl_ih + 3);

    auto ts_xxxxxx_xxxyz = pbuffer.data(idx_ovl_ih + 4);

    auto ts_xxxxxx_xxxzz = pbuffer.data(idx_ovl_ih + 5);

    auto ts_xxxxxx_xxyyy = pbuffer.data(idx_ovl_ih + 6);

    auto ts_xxxxxx_xxyyz = pbuffer.data(idx_ovl_ih + 7);

    auto ts_xxxxxx_xxyzz = pbuffer.data(idx_ovl_ih + 8);

    auto ts_xxxxxx_xxzzz = pbuffer.data(idx_ovl_ih + 9);

    auto ts_xxxxxx_xyyyy = pbuffer.data(idx_ovl_ih + 10);

    auto ts_xxxxxx_xyyyz = pbuffer.data(idx_ovl_ih + 11);

    auto ts_xxxxxx_xyyzz = pbuffer.data(idx_ovl_ih + 12);

    auto ts_xxxxxx_xyzzz = pbuffer.data(idx_ovl_ih + 13);

    auto ts_xxxxxx_xzzzz = pbuffer.data(idx_ovl_ih + 14);

    auto ts_xxxxxx_yyyyy = pbuffer.data(idx_ovl_ih + 15);

    auto ts_xxxxxx_yyyyz = pbuffer.data(idx_ovl_ih + 16);

    auto ts_xxxxxx_yyyzz = pbuffer.data(idx_ovl_ih + 17);

    auto ts_xxxxxx_yyzzz = pbuffer.data(idx_ovl_ih + 18);

    auto ts_xxxxxx_yzzzz = pbuffer.data(idx_ovl_ih + 19);

    auto ts_xxxxxx_zzzzz = pbuffer.data(idx_ovl_ih + 20);

    #pragma omp simd aligned(pa_x, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxxz, ts_xxxx_xxxyy, ts_xxxx_xxxyz, ts_xxxx_xxxzz, ts_xxxx_xxyyy, ts_xxxx_xxyyz, ts_xxxx_xxyzz, ts_xxxx_xxzzz, ts_xxxx_xyyyy, ts_xxxx_xyyyz, ts_xxxx_xyyzz, ts_xxxx_xyzzz, ts_xxxx_xzzzz, ts_xxxx_yyyyy, ts_xxxx_yyyyz, ts_xxxx_yyyzz, ts_xxxx_yyzzz, ts_xxxx_yzzzz, ts_xxxx_zzzzz, ts_xxxxx_xxxx, ts_xxxxx_xxxxx, ts_xxxxx_xxxxy, ts_xxxxx_xxxxz, ts_xxxxx_xxxy, ts_xxxxx_xxxyy, ts_xxxxx_xxxyz, ts_xxxxx_xxxz, ts_xxxxx_xxxzz, ts_xxxxx_xxyy, ts_xxxxx_xxyyy, ts_xxxxx_xxyyz, ts_xxxxx_xxyz, ts_xxxxx_xxyzz, ts_xxxxx_xxzz, ts_xxxxx_xxzzz, ts_xxxxx_xyyy, ts_xxxxx_xyyyy, ts_xxxxx_xyyyz, ts_xxxxx_xyyz, ts_xxxxx_xyyzz, ts_xxxxx_xyzz, ts_xxxxx_xyzzz, ts_xxxxx_xzzz, ts_xxxxx_xzzzz, ts_xxxxx_yyyy, ts_xxxxx_yyyyy, ts_xxxxx_yyyyz, ts_xxxxx_yyyz, ts_xxxxx_yyyzz, ts_xxxxx_yyzz, ts_xxxxx_yyzzz, ts_xxxxx_yzzz, ts_xxxxx_yzzzz, ts_xxxxx_zzzz, ts_xxxxx_zzzzz, ts_xxxxxx_xxxxx, ts_xxxxxx_xxxxy, ts_xxxxxx_xxxxz, ts_xxxxxx_xxxyy, ts_xxxxxx_xxxyz, ts_xxxxxx_xxxzz, ts_xxxxxx_xxyyy, ts_xxxxxx_xxyyz, ts_xxxxxx_xxyzz, ts_xxxxxx_xxzzz, ts_xxxxxx_xyyyy, ts_xxxxxx_xyyyz, ts_xxxxxx_xyyzz, ts_xxxxxx_xyzzz, ts_xxxxxx_xzzzz, ts_xxxxxx_yyyyy, ts_xxxxxx_yyyyz, ts_xxxxxx_yyyzz, ts_xxxxxx_yyzzz, ts_xxxxxx_yzzzz, ts_xxxxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxx_xxxxx[i] = 5.0 * ts_xxxx_xxxxx[i] * fe_0 + 5.0 * ts_xxxxx_xxxx[i] * fe_0 + ts_xxxxx_xxxxx[i] * pa_x[i];

        ts_xxxxxx_xxxxy[i] = 5.0 * ts_xxxx_xxxxy[i] * fe_0 + 4.0 * ts_xxxxx_xxxy[i] * fe_0 + ts_xxxxx_xxxxy[i] * pa_x[i];

        ts_xxxxxx_xxxxz[i] = 5.0 * ts_xxxx_xxxxz[i] * fe_0 + 4.0 * ts_xxxxx_xxxz[i] * fe_0 + ts_xxxxx_xxxxz[i] * pa_x[i];

        ts_xxxxxx_xxxyy[i] = 5.0 * ts_xxxx_xxxyy[i] * fe_0 + 3.0 * ts_xxxxx_xxyy[i] * fe_0 + ts_xxxxx_xxxyy[i] * pa_x[i];

        ts_xxxxxx_xxxyz[i] = 5.0 * ts_xxxx_xxxyz[i] * fe_0 + 3.0 * ts_xxxxx_xxyz[i] * fe_0 + ts_xxxxx_xxxyz[i] * pa_x[i];

        ts_xxxxxx_xxxzz[i] = 5.0 * ts_xxxx_xxxzz[i] * fe_0 + 3.0 * ts_xxxxx_xxzz[i] * fe_0 + ts_xxxxx_xxxzz[i] * pa_x[i];

        ts_xxxxxx_xxyyy[i] = 5.0 * ts_xxxx_xxyyy[i] * fe_0 + 2.0 * ts_xxxxx_xyyy[i] * fe_0 + ts_xxxxx_xxyyy[i] * pa_x[i];

        ts_xxxxxx_xxyyz[i] = 5.0 * ts_xxxx_xxyyz[i] * fe_0 + 2.0 * ts_xxxxx_xyyz[i] * fe_0 + ts_xxxxx_xxyyz[i] * pa_x[i];

        ts_xxxxxx_xxyzz[i] = 5.0 * ts_xxxx_xxyzz[i] * fe_0 + 2.0 * ts_xxxxx_xyzz[i] * fe_0 + ts_xxxxx_xxyzz[i] * pa_x[i];

        ts_xxxxxx_xxzzz[i] = 5.0 * ts_xxxx_xxzzz[i] * fe_0 + 2.0 * ts_xxxxx_xzzz[i] * fe_0 + ts_xxxxx_xxzzz[i] * pa_x[i];

        ts_xxxxxx_xyyyy[i] = 5.0 * ts_xxxx_xyyyy[i] * fe_0 + ts_xxxxx_yyyy[i] * fe_0 + ts_xxxxx_xyyyy[i] * pa_x[i];

        ts_xxxxxx_xyyyz[i] = 5.0 * ts_xxxx_xyyyz[i] * fe_0 + ts_xxxxx_yyyz[i] * fe_0 + ts_xxxxx_xyyyz[i] * pa_x[i];

        ts_xxxxxx_xyyzz[i] = 5.0 * ts_xxxx_xyyzz[i] * fe_0 + ts_xxxxx_yyzz[i] * fe_0 + ts_xxxxx_xyyzz[i] * pa_x[i];

        ts_xxxxxx_xyzzz[i] = 5.0 * ts_xxxx_xyzzz[i] * fe_0 + ts_xxxxx_yzzz[i] * fe_0 + ts_xxxxx_xyzzz[i] * pa_x[i];

        ts_xxxxxx_xzzzz[i] = 5.0 * ts_xxxx_xzzzz[i] * fe_0 + ts_xxxxx_zzzz[i] * fe_0 + ts_xxxxx_xzzzz[i] * pa_x[i];

        ts_xxxxxx_yyyyy[i] = 5.0 * ts_xxxx_yyyyy[i] * fe_0 + ts_xxxxx_yyyyy[i] * pa_x[i];

        ts_xxxxxx_yyyyz[i] = 5.0 * ts_xxxx_yyyyz[i] * fe_0 + ts_xxxxx_yyyyz[i] * pa_x[i];

        ts_xxxxxx_yyyzz[i] = 5.0 * ts_xxxx_yyyzz[i] * fe_0 + ts_xxxxx_yyyzz[i] * pa_x[i];

        ts_xxxxxx_yyzzz[i] = 5.0 * ts_xxxx_yyzzz[i] * fe_0 + ts_xxxxx_yyzzz[i] * pa_x[i];

        ts_xxxxxx_yzzzz[i] = 5.0 * ts_xxxx_yzzzz[i] * fe_0 + ts_xxxxx_yzzzz[i] * pa_x[i];

        ts_xxxxxx_zzzzz[i] = 5.0 * ts_xxxx_zzzzz[i] * fe_0 + ts_xxxxx_zzzzz[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : IH

    auto ts_xxxxxy_xxxxx = pbuffer.data(idx_ovl_ih + 21);

    auto ts_xxxxxy_xxxxy = pbuffer.data(idx_ovl_ih + 22);

    auto ts_xxxxxy_xxxxz = pbuffer.data(idx_ovl_ih + 23);

    auto ts_xxxxxy_xxxyy = pbuffer.data(idx_ovl_ih + 24);

    auto ts_xxxxxy_xxxyz = pbuffer.data(idx_ovl_ih + 25);

    auto ts_xxxxxy_xxxzz = pbuffer.data(idx_ovl_ih + 26);

    auto ts_xxxxxy_xxyyy = pbuffer.data(idx_ovl_ih + 27);

    auto ts_xxxxxy_xxyyz = pbuffer.data(idx_ovl_ih + 28);

    auto ts_xxxxxy_xxyzz = pbuffer.data(idx_ovl_ih + 29);

    auto ts_xxxxxy_xxzzz = pbuffer.data(idx_ovl_ih + 30);

    auto ts_xxxxxy_xyyyy = pbuffer.data(idx_ovl_ih + 31);

    auto ts_xxxxxy_xyyyz = pbuffer.data(idx_ovl_ih + 32);

    auto ts_xxxxxy_xyyzz = pbuffer.data(idx_ovl_ih + 33);

    auto ts_xxxxxy_xyzzz = pbuffer.data(idx_ovl_ih + 34);

    auto ts_xxxxxy_xzzzz = pbuffer.data(idx_ovl_ih + 35);

    auto ts_xxxxxy_yyyyy = pbuffer.data(idx_ovl_ih + 36);

    auto ts_xxxxxy_yyyyz = pbuffer.data(idx_ovl_ih + 37);

    auto ts_xxxxxy_yyyzz = pbuffer.data(idx_ovl_ih + 38);

    auto ts_xxxxxy_yyzzz = pbuffer.data(idx_ovl_ih + 39);

    auto ts_xxxxxy_yzzzz = pbuffer.data(idx_ovl_ih + 40);

    auto ts_xxxxxy_zzzzz = pbuffer.data(idx_ovl_ih + 41);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxxx_xxxx, ts_xxxxx_xxxxx, ts_xxxxx_xxxxy, ts_xxxxx_xxxxz, ts_xxxxx_xxxy, ts_xxxxx_xxxyy, ts_xxxxx_xxxyz, ts_xxxxx_xxxz, ts_xxxxx_xxxzz, ts_xxxxx_xxyy, ts_xxxxx_xxyyy, ts_xxxxx_xxyyz, ts_xxxxx_xxyz, ts_xxxxx_xxyzz, ts_xxxxx_xxzz, ts_xxxxx_xxzzz, ts_xxxxx_xyyy, ts_xxxxx_xyyyy, ts_xxxxx_xyyyz, ts_xxxxx_xyyz, ts_xxxxx_xyyzz, ts_xxxxx_xyzz, ts_xxxxx_xyzzz, ts_xxxxx_xzzz, ts_xxxxx_xzzzz, ts_xxxxx_zzzzz, ts_xxxxxy_xxxxx, ts_xxxxxy_xxxxy, ts_xxxxxy_xxxxz, ts_xxxxxy_xxxyy, ts_xxxxxy_xxxyz, ts_xxxxxy_xxxzz, ts_xxxxxy_xxyyy, ts_xxxxxy_xxyyz, ts_xxxxxy_xxyzz, ts_xxxxxy_xxzzz, ts_xxxxxy_xyyyy, ts_xxxxxy_xyyyz, ts_xxxxxy_xyyzz, ts_xxxxxy_xyzzz, ts_xxxxxy_xzzzz, ts_xxxxxy_yyyyy, ts_xxxxxy_yyyyz, ts_xxxxxy_yyyzz, ts_xxxxxy_yyzzz, ts_xxxxxy_yzzzz, ts_xxxxxy_zzzzz, ts_xxxxy_yyyyy, ts_xxxxy_yyyyz, ts_xxxxy_yyyzz, ts_xxxxy_yyzzz, ts_xxxxy_yzzzz, ts_xxxy_yyyyy, ts_xxxy_yyyyz, ts_xxxy_yyyzz, ts_xxxy_yyzzz, ts_xxxy_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxy_xxxxx[i] = ts_xxxxx_xxxxx[i] * pa_y[i];

        ts_xxxxxy_xxxxy[i] = ts_xxxxx_xxxx[i] * fe_0 + ts_xxxxx_xxxxy[i] * pa_y[i];

        ts_xxxxxy_xxxxz[i] = ts_xxxxx_xxxxz[i] * pa_y[i];

        ts_xxxxxy_xxxyy[i] = 2.0 * ts_xxxxx_xxxy[i] * fe_0 + ts_xxxxx_xxxyy[i] * pa_y[i];

        ts_xxxxxy_xxxyz[i] = ts_xxxxx_xxxz[i] * fe_0 + ts_xxxxx_xxxyz[i] * pa_y[i];

        ts_xxxxxy_xxxzz[i] = ts_xxxxx_xxxzz[i] * pa_y[i];

        ts_xxxxxy_xxyyy[i] = 3.0 * ts_xxxxx_xxyy[i] * fe_0 + ts_xxxxx_xxyyy[i] * pa_y[i];

        ts_xxxxxy_xxyyz[i] = 2.0 * ts_xxxxx_xxyz[i] * fe_0 + ts_xxxxx_xxyyz[i] * pa_y[i];

        ts_xxxxxy_xxyzz[i] = ts_xxxxx_xxzz[i] * fe_0 + ts_xxxxx_xxyzz[i] * pa_y[i];

        ts_xxxxxy_xxzzz[i] = ts_xxxxx_xxzzz[i] * pa_y[i];

        ts_xxxxxy_xyyyy[i] = 4.0 * ts_xxxxx_xyyy[i] * fe_0 + ts_xxxxx_xyyyy[i] * pa_y[i];

        ts_xxxxxy_xyyyz[i] = 3.0 * ts_xxxxx_xyyz[i] * fe_0 + ts_xxxxx_xyyyz[i] * pa_y[i];

        ts_xxxxxy_xyyzz[i] = 2.0 * ts_xxxxx_xyzz[i] * fe_0 + ts_xxxxx_xyyzz[i] * pa_y[i];

        ts_xxxxxy_xyzzz[i] = ts_xxxxx_xzzz[i] * fe_0 + ts_xxxxx_xyzzz[i] * pa_y[i];

        ts_xxxxxy_xzzzz[i] = ts_xxxxx_xzzzz[i] * pa_y[i];

        ts_xxxxxy_yyyyy[i] = 4.0 * ts_xxxy_yyyyy[i] * fe_0 + ts_xxxxy_yyyyy[i] * pa_x[i];

        ts_xxxxxy_yyyyz[i] = 4.0 * ts_xxxy_yyyyz[i] * fe_0 + ts_xxxxy_yyyyz[i] * pa_x[i];

        ts_xxxxxy_yyyzz[i] = 4.0 * ts_xxxy_yyyzz[i] * fe_0 + ts_xxxxy_yyyzz[i] * pa_x[i];

        ts_xxxxxy_yyzzz[i] = 4.0 * ts_xxxy_yyzzz[i] * fe_0 + ts_xxxxy_yyzzz[i] * pa_x[i];

        ts_xxxxxy_yzzzz[i] = 4.0 * ts_xxxy_yzzzz[i] * fe_0 + ts_xxxxy_yzzzz[i] * pa_x[i];

        ts_xxxxxy_zzzzz[i] = ts_xxxxx_zzzzz[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : IH

    auto ts_xxxxxz_xxxxx = pbuffer.data(idx_ovl_ih + 42);

    auto ts_xxxxxz_xxxxy = pbuffer.data(idx_ovl_ih + 43);

    auto ts_xxxxxz_xxxxz = pbuffer.data(idx_ovl_ih + 44);

    auto ts_xxxxxz_xxxyy = pbuffer.data(idx_ovl_ih + 45);

    auto ts_xxxxxz_xxxyz = pbuffer.data(idx_ovl_ih + 46);

    auto ts_xxxxxz_xxxzz = pbuffer.data(idx_ovl_ih + 47);

    auto ts_xxxxxz_xxyyy = pbuffer.data(idx_ovl_ih + 48);

    auto ts_xxxxxz_xxyyz = pbuffer.data(idx_ovl_ih + 49);

    auto ts_xxxxxz_xxyzz = pbuffer.data(idx_ovl_ih + 50);

    auto ts_xxxxxz_xxzzz = pbuffer.data(idx_ovl_ih + 51);

    auto ts_xxxxxz_xyyyy = pbuffer.data(idx_ovl_ih + 52);

    auto ts_xxxxxz_xyyyz = pbuffer.data(idx_ovl_ih + 53);

    auto ts_xxxxxz_xyyzz = pbuffer.data(idx_ovl_ih + 54);

    auto ts_xxxxxz_xyzzz = pbuffer.data(idx_ovl_ih + 55);

    auto ts_xxxxxz_xzzzz = pbuffer.data(idx_ovl_ih + 56);

    auto ts_xxxxxz_yyyyy = pbuffer.data(idx_ovl_ih + 57);

    auto ts_xxxxxz_yyyyz = pbuffer.data(idx_ovl_ih + 58);

    auto ts_xxxxxz_yyyzz = pbuffer.data(idx_ovl_ih + 59);

    auto ts_xxxxxz_yyzzz = pbuffer.data(idx_ovl_ih + 60);

    auto ts_xxxxxz_yzzzz = pbuffer.data(idx_ovl_ih + 61);

    auto ts_xxxxxz_zzzzz = pbuffer.data(idx_ovl_ih + 62);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxxx_xxxx, ts_xxxxx_xxxxx, ts_xxxxx_xxxxy, ts_xxxxx_xxxxz, ts_xxxxx_xxxy, ts_xxxxx_xxxyy, ts_xxxxx_xxxyz, ts_xxxxx_xxxz, ts_xxxxx_xxxzz, ts_xxxxx_xxyy, ts_xxxxx_xxyyy, ts_xxxxx_xxyyz, ts_xxxxx_xxyz, ts_xxxxx_xxyzz, ts_xxxxx_xxzz, ts_xxxxx_xxzzz, ts_xxxxx_xyyy, ts_xxxxx_xyyyy, ts_xxxxx_xyyyz, ts_xxxxx_xyyz, ts_xxxxx_xyyzz, ts_xxxxx_xyzz, ts_xxxxx_xyzzz, ts_xxxxx_xzzz, ts_xxxxx_xzzzz, ts_xxxxx_yyyyy, ts_xxxxxz_xxxxx, ts_xxxxxz_xxxxy, ts_xxxxxz_xxxxz, ts_xxxxxz_xxxyy, ts_xxxxxz_xxxyz, ts_xxxxxz_xxxzz, ts_xxxxxz_xxyyy, ts_xxxxxz_xxyyz, ts_xxxxxz_xxyzz, ts_xxxxxz_xxzzz, ts_xxxxxz_xyyyy, ts_xxxxxz_xyyyz, ts_xxxxxz_xyyzz, ts_xxxxxz_xyzzz, ts_xxxxxz_xzzzz, ts_xxxxxz_yyyyy, ts_xxxxxz_yyyyz, ts_xxxxxz_yyyzz, ts_xxxxxz_yyzzz, ts_xxxxxz_yzzzz, ts_xxxxxz_zzzzz, ts_xxxxz_yyyyz, ts_xxxxz_yyyzz, ts_xxxxz_yyzzz, ts_xxxxz_yzzzz, ts_xxxxz_zzzzz, ts_xxxz_yyyyz, ts_xxxz_yyyzz, ts_xxxz_yyzzz, ts_xxxz_yzzzz, ts_xxxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxz_xxxxx[i] = ts_xxxxx_xxxxx[i] * pa_z[i];

        ts_xxxxxz_xxxxy[i] = ts_xxxxx_xxxxy[i] * pa_z[i];

        ts_xxxxxz_xxxxz[i] = ts_xxxxx_xxxx[i] * fe_0 + ts_xxxxx_xxxxz[i] * pa_z[i];

        ts_xxxxxz_xxxyy[i] = ts_xxxxx_xxxyy[i] * pa_z[i];

        ts_xxxxxz_xxxyz[i] = ts_xxxxx_xxxy[i] * fe_0 + ts_xxxxx_xxxyz[i] * pa_z[i];

        ts_xxxxxz_xxxzz[i] = 2.0 * ts_xxxxx_xxxz[i] * fe_0 + ts_xxxxx_xxxzz[i] * pa_z[i];

        ts_xxxxxz_xxyyy[i] = ts_xxxxx_xxyyy[i] * pa_z[i];

        ts_xxxxxz_xxyyz[i] = ts_xxxxx_xxyy[i] * fe_0 + ts_xxxxx_xxyyz[i] * pa_z[i];

        ts_xxxxxz_xxyzz[i] = 2.0 * ts_xxxxx_xxyz[i] * fe_0 + ts_xxxxx_xxyzz[i] * pa_z[i];

        ts_xxxxxz_xxzzz[i] = 3.0 * ts_xxxxx_xxzz[i] * fe_0 + ts_xxxxx_xxzzz[i] * pa_z[i];

        ts_xxxxxz_xyyyy[i] = ts_xxxxx_xyyyy[i] * pa_z[i];

        ts_xxxxxz_xyyyz[i] = ts_xxxxx_xyyy[i] * fe_0 + ts_xxxxx_xyyyz[i] * pa_z[i];

        ts_xxxxxz_xyyzz[i] = 2.0 * ts_xxxxx_xyyz[i] * fe_0 + ts_xxxxx_xyyzz[i] * pa_z[i];

        ts_xxxxxz_xyzzz[i] = 3.0 * ts_xxxxx_xyzz[i] * fe_0 + ts_xxxxx_xyzzz[i] * pa_z[i];

        ts_xxxxxz_xzzzz[i] = 4.0 * ts_xxxxx_xzzz[i] * fe_0 + ts_xxxxx_xzzzz[i] * pa_z[i];

        ts_xxxxxz_yyyyy[i] = ts_xxxxx_yyyyy[i] * pa_z[i];

        ts_xxxxxz_yyyyz[i] = 4.0 * ts_xxxz_yyyyz[i] * fe_0 + ts_xxxxz_yyyyz[i] * pa_x[i];

        ts_xxxxxz_yyyzz[i] = 4.0 * ts_xxxz_yyyzz[i] * fe_0 + ts_xxxxz_yyyzz[i] * pa_x[i];

        ts_xxxxxz_yyzzz[i] = 4.0 * ts_xxxz_yyzzz[i] * fe_0 + ts_xxxxz_yyzzz[i] * pa_x[i];

        ts_xxxxxz_yzzzz[i] = 4.0 * ts_xxxz_yzzzz[i] * fe_0 + ts_xxxxz_yzzzz[i] * pa_x[i];

        ts_xxxxxz_zzzzz[i] = 4.0 * ts_xxxz_zzzzz[i] * fe_0 + ts_xxxxz_zzzzz[i] * pa_x[i];
    }

    // Set up 63-84 components of targeted buffer : IH

    auto ts_xxxxyy_xxxxx = pbuffer.data(idx_ovl_ih + 63);

    auto ts_xxxxyy_xxxxy = pbuffer.data(idx_ovl_ih + 64);

    auto ts_xxxxyy_xxxxz = pbuffer.data(idx_ovl_ih + 65);

    auto ts_xxxxyy_xxxyy = pbuffer.data(idx_ovl_ih + 66);

    auto ts_xxxxyy_xxxyz = pbuffer.data(idx_ovl_ih + 67);

    auto ts_xxxxyy_xxxzz = pbuffer.data(idx_ovl_ih + 68);

    auto ts_xxxxyy_xxyyy = pbuffer.data(idx_ovl_ih + 69);

    auto ts_xxxxyy_xxyyz = pbuffer.data(idx_ovl_ih + 70);

    auto ts_xxxxyy_xxyzz = pbuffer.data(idx_ovl_ih + 71);

    auto ts_xxxxyy_xxzzz = pbuffer.data(idx_ovl_ih + 72);

    auto ts_xxxxyy_xyyyy = pbuffer.data(idx_ovl_ih + 73);

    auto ts_xxxxyy_xyyyz = pbuffer.data(idx_ovl_ih + 74);

    auto ts_xxxxyy_xyyzz = pbuffer.data(idx_ovl_ih + 75);

    auto ts_xxxxyy_xyzzz = pbuffer.data(idx_ovl_ih + 76);

    auto ts_xxxxyy_xzzzz = pbuffer.data(idx_ovl_ih + 77);

    auto ts_xxxxyy_yyyyy = pbuffer.data(idx_ovl_ih + 78);

    auto ts_xxxxyy_yyyyz = pbuffer.data(idx_ovl_ih + 79);

    auto ts_xxxxyy_yyyzz = pbuffer.data(idx_ovl_ih + 80);

    auto ts_xxxxyy_yyzzz = pbuffer.data(idx_ovl_ih + 81);

    auto ts_xxxxyy_yzzzz = pbuffer.data(idx_ovl_ih + 82);

    auto ts_xxxxyy_zzzzz = pbuffer.data(idx_ovl_ih + 83);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxx_xxxxx, ts_xxxx_xxxxz, ts_xxxx_xxxzz, ts_xxxx_xxzzz, ts_xxxx_xzzzz, ts_xxxxy_xxxxx, ts_xxxxy_xxxxz, ts_xxxxy_xxxzz, ts_xxxxy_xxzzz, ts_xxxxy_xzzzz, ts_xxxxyy_xxxxx, ts_xxxxyy_xxxxy, ts_xxxxyy_xxxxz, ts_xxxxyy_xxxyy, ts_xxxxyy_xxxyz, ts_xxxxyy_xxxzz, ts_xxxxyy_xxyyy, ts_xxxxyy_xxyyz, ts_xxxxyy_xxyzz, ts_xxxxyy_xxzzz, ts_xxxxyy_xyyyy, ts_xxxxyy_xyyyz, ts_xxxxyy_xyyzz, ts_xxxxyy_xyzzz, ts_xxxxyy_xzzzz, ts_xxxxyy_yyyyy, ts_xxxxyy_yyyyz, ts_xxxxyy_yyyzz, ts_xxxxyy_yyzzz, ts_xxxxyy_yzzzz, ts_xxxxyy_zzzzz, ts_xxxyy_xxxxy, ts_xxxyy_xxxy, ts_xxxyy_xxxyy, ts_xxxyy_xxxyz, ts_xxxyy_xxyy, ts_xxxyy_xxyyy, ts_xxxyy_xxyyz, ts_xxxyy_xxyz, ts_xxxyy_xxyzz, ts_xxxyy_xyyy, ts_xxxyy_xyyyy, ts_xxxyy_xyyyz, ts_xxxyy_xyyz, ts_xxxyy_xyyzz, ts_xxxyy_xyzz, ts_xxxyy_xyzzz, ts_xxxyy_yyyy, ts_xxxyy_yyyyy, ts_xxxyy_yyyyz, ts_xxxyy_yyyz, ts_xxxyy_yyyzz, ts_xxxyy_yyzz, ts_xxxyy_yyzzz, ts_xxxyy_yzzz, ts_xxxyy_yzzzz, ts_xxxyy_zzzzz, ts_xxyy_xxxxy, ts_xxyy_xxxyy, ts_xxyy_xxxyz, ts_xxyy_xxyyy, ts_xxyy_xxyyz, ts_xxyy_xxyzz, ts_xxyy_xyyyy, ts_xxyy_xyyyz, ts_xxyy_xyyzz, ts_xxyy_xyzzz, ts_xxyy_yyyyy, ts_xxyy_yyyyz, ts_xxyy_yyyzz, ts_xxyy_yyzzz, ts_xxyy_yzzzz, ts_xxyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxyy_xxxxx[i] = ts_xxxx_xxxxx[i] * fe_0 + ts_xxxxy_xxxxx[i] * pa_y[i];

        ts_xxxxyy_xxxxy[i] = 3.0 * ts_xxyy_xxxxy[i] * fe_0 + 4.0 * ts_xxxyy_xxxy[i] * fe_0 + ts_xxxyy_xxxxy[i] * pa_x[i];

        ts_xxxxyy_xxxxz[i] = ts_xxxx_xxxxz[i] * fe_0 + ts_xxxxy_xxxxz[i] * pa_y[i];

        ts_xxxxyy_xxxyy[i] = 3.0 * ts_xxyy_xxxyy[i] * fe_0 + 3.0 * ts_xxxyy_xxyy[i] * fe_0 + ts_xxxyy_xxxyy[i] * pa_x[i];

        ts_xxxxyy_xxxyz[i] = 3.0 * ts_xxyy_xxxyz[i] * fe_0 + 3.0 * ts_xxxyy_xxyz[i] * fe_0 + ts_xxxyy_xxxyz[i] * pa_x[i];

        ts_xxxxyy_xxxzz[i] = ts_xxxx_xxxzz[i] * fe_0 + ts_xxxxy_xxxzz[i] * pa_y[i];

        ts_xxxxyy_xxyyy[i] = 3.0 * ts_xxyy_xxyyy[i] * fe_0 + 2.0 * ts_xxxyy_xyyy[i] * fe_0 + ts_xxxyy_xxyyy[i] * pa_x[i];

        ts_xxxxyy_xxyyz[i] = 3.0 * ts_xxyy_xxyyz[i] * fe_0 + 2.0 * ts_xxxyy_xyyz[i] * fe_0 + ts_xxxyy_xxyyz[i] * pa_x[i];

        ts_xxxxyy_xxyzz[i] = 3.0 * ts_xxyy_xxyzz[i] * fe_0 + 2.0 * ts_xxxyy_xyzz[i] * fe_0 + ts_xxxyy_xxyzz[i] * pa_x[i];

        ts_xxxxyy_xxzzz[i] = ts_xxxx_xxzzz[i] * fe_0 + ts_xxxxy_xxzzz[i] * pa_y[i];

        ts_xxxxyy_xyyyy[i] = 3.0 * ts_xxyy_xyyyy[i] * fe_0 + ts_xxxyy_yyyy[i] * fe_0 + ts_xxxyy_xyyyy[i] * pa_x[i];

        ts_xxxxyy_xyyyz[i] = 3.0 * ts_xxyy_xyyyz[i] * fe_0 + ts_xxxyy_yyyz[i] * fe_0 + ts_xxxyy_xyyyz[i] * pa_x[i];

        ts_xxxxyy_xyyzz[i] = 3.0 * ts_xxyy_xyyzz[i] * fe_0 + ts_xxxyy_yyzz[i] * fe_0 + ts_xxxyy_xyyzz[i] * pa_x[i];

        ts_xxxxyy_xyzzz[i] = 3.0 * ts_xxyy_xyzzz[i] * fe_0 + ts_xxxyy_yzzz[i] * fe_0 + ts_xxxyy_xyzzz[i] * pa_x[i];

        ts_xxxxyy_xzzzz[i] = ts_xxxx_xzzzz[i] * fe_0 + ts_xxxxy_xzzzz[i] * pa_y[i];

        ts_xxxxyy_yyyyy[i] = 3.0 * ts_xxyy_yyyyy[i] * fe_0 + ts_xxxyy_yyyyy[i] * pa_x[i];

        ts_xxxxyy_yyyyz[i] = 3.0 * ts_xxyy_yyyyz[i] * fe_0 + ts_xxxyy_yyyyz[i] * pa_x[i];

        ts_xxxxyy_yyyzz[i] = 3.0 * ts_xxyy_yyyzz[i] * fe_0 + ts_xxxyy_yyyzz[i] * pa_x[i];

        ts_xxxxyy_yyzzz[i] = 3.0 * ts_xxyy_yyzzz[i] * fe_0 + ts_xxxyy_yyzzz[i] * pa_x[i];

        ts_xxxxyy_yzzzz[i] = 3.0 * ts_xxyy_yzzzz[i] * fe_0 + ts_xxxyy_yzzzz[i] * pa_x[i];

        ts_xxxxyy_zzzzz[i] = 3.0 * ts_xxyy_zzzzz[i] * fe_0 + ts_xxxyy_zzzzz[i] * pa_x[i];
    }

    // Set up 84-105 components of targeted buffer : IH

    auto ts_xxxxyz_xxxxx = pbuffer.data(idx_ovl_ih + 84);

    auto ts_xxxxyz_xxxxy = pbuffer.data(idx_ovl_ih + 85);

    auto ts_xxxxyz_xxxxz = pbuffer.data(idx_ovl_ih + 86);

    auto ts_xxxxyz_xxxyy = pbuffer.data(idx_ovl_ih + 87);

    auto ts_xxxxyz_xxxyz = pbuffer.data(idx_ovl_ih + 88);

    auto ts_xxxxyz_xxxzz = pbuffer.data(idx_ovl_ih + 89);

    auto ts_xxxxyz_xxyyy = pbuffer.data(idx_ovl_ih + 90);

    auto ts_xxxxyz_xxyyz = pbuffer.data(idx_ovl_ih + 91);

    auto ts_xxxxyz_xxyzz = pbuffer.data(idx_ovl_ih + 92);

    auto ts_xxxxyz_xxzzz = pbuffer.data(idx_ovl_ih + 93);

    auto ts_xxxxyz_xyyyy = pbuffer.data(idx_ovl_ih + 94);

    auto ts_xxxxyz_xyyyz = pbuffer.data(idx_ovl_ih + 95);

    auto ts_xxxxyz_xyyzz = pbuffer.data(idx_ovl_ih + 96);

    auto ts_xxxxyz_xyzzz = pbuffer.data(idx_ovl_ih + 97);

    auto ts_xxxxyz_xzzzz = pbuffer.data(idx_ovl_ih + 98);

    auto ts_xxxxyz_yyyyy = pbuffer.data(idx_ovl_ih + 99);

    auto ts_xxxxyz_yyyyz = pbuffer.data(idx_ovl_ih + 100);

    auto ts_xxxxyz_yyyzz = pbuffer.data(idx_ovl_ih + 101);

    auto ts_xxxxyz_yyzzz = pbuffer.data(idx_ovl_ih + 102);

    auto ts_xxxxyz_yzzzz = pbuffer.data(idx_ovl_ih + 103);

    auto ts_xxxxyz_zzzzz = pbuffer.data(idx_ovl_ih + 104);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxxy_xxxxy, ts_xxxxy_xxxyy, ts_xxxxy_xxyyy, ts_xxxxy_xyyyy, ts_xxxxy_yyyyy, ts_xxxxyz_xxxxx, ts_xxxxyz_xxxxy, ts_xxxxyz_xxxxz, ts_xxxxyz_xxxyy, ts_xxxxyz_xxxyz, ts_xxxxyz_xxxzz, ts_xxxxyz_xxyyy, ts_xxxxyz_xxyyz, ts_xxxxyz_xxyzz, ts_xxxxyz_xxzzz, ts_xxxxyz_xyyyy, ts_xxxxyz_xyyyz, ts_xxxxyz_xyyzz, ts_xxxxyz_xyzzz, ts_xxxxyz_xzzzz, ts_xxxxyz_yyyyy, ts_xxxxyz_yyyyz, ts_xxxxyz_yyyzz, ts_xxxxyz_yyzzz, ts_xxxxyz_yzzzz, ts_xxxxyz_zzzzz, ts_xxxxz_xxxxx, ts_xxxxz_xxxxz, ts_xxxxz_xxxyz, ts_xxxxz_xxxz, ts_xxxxz_xxxzz, ts_xxxxz_xxyyz, ts_xxxxz_xxyz, ts_xxxxz_xxyzz, ts_xxxxz_xxzz, ts_xxxxz_xxzzz, ts_xxxxz_xyyyz, ts_xxxxz_xyyz, ts_xxxxz_xyyzz, ts_xxxxz_xyzz, ts_xxxxz_xyzzz, ts_xxxxz_xzzz, ts_xxxxz_xzzzz, ts_xxxxz_zzzzz, ts_xxxyz_yyyyz, ts_xxxyz_yyyzz, ts_xxxyz_yyzzz, ts_xxxyz_yzzzz, ts_xxyz_yyyyz, ts_xxyz_yyyzz, ts_xxyz_yyzzz, ts_xxyz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxyz_xxxxx[i] = ts_xxxxz_xxxxx[i] * pa_y[i];

        ts_xxxxyz_xxxxy[i] = ts_xxxxy_xxxxy[i] * pa_z[i];

        ts_xxxxyz_xxxxz[i] = ts_xxxxz_xxxxz[i] * pa_y[i];

        ts_xxxxyz_xxxyy[i] = ts_xxxxy_xxxyy[i] * pa_z[i];

        ts_xxxxyz_xxxyz[i] = ts_xxxxz_xxxz[i] * fe_0 + ts_xxxxz_xxxyz[i] * pa_y[i];

        ts_xxxxyz_xxxzz[i] = ts_xxxxz_xxxzz[i] * pa_y[i];

        ts_xxxxyz_xxyyy[i] = ts_xxxxy_xxyyy[i] * pa_z[i];

        ts_xxxxyz_xxyyz[i] = 2.0 * ts_xxxxz_xxyz[i] * fe_0 + ts_xxxxz_xxyyz[i] * pa_y[i];

        ts_xxxxyz_xxyzz[i] = ts_xxxxz_xxzz[i] * fe_0 + ts_xxxxz_xxyzz[i] * pa_y[i];

        ts_xxxxyz_xxzzz[i] = ts_xxxxz_xxzzz[i] * pa_y[i];

        ts_xxxxyz_xyyyy[i] = ts_xxxxy_xyyyy[i] * pa_z[i];

        ts_xxxxyz_xyyyz[i] = 3.0 * ts_xxxxz_xyyz[i] * fe_0 + ts_xxxxz_xyyyz[i] * pa_y[i];

        ts_xxxxyz_xyyzz[i] = 2.0 * ts_xxxxz_xyzz[i] * fe_0 + ts_xxxxz_xyyzz[i] * pa_y[i];

        ts_xxxxyz_xyzzz[i] = ts_xxxxz_xzzz[i] * fe_0 + ts_xxxxz_xyzzz[i] * pa_y[i];

        ts_xxxxyz_xzzzz[i] = ts_xxxxz_xzzzz[i] * pa_y[i];

        ts_xxxxyz_yyyyy[i] = ts_xxxxy_yyyyy[i] * pa_z[i];

        ts_xxxxyz_yyyyz[i] = 3.0 * ts_xxyz_yyyyz[i] * fe_0 + ts_xxxyz_yyyyz[i] * pa_x[i];

        ts_xxxxyz_yyyzz[i] = 3.0 * ts_xxyz_yyyzz[i] * fe_0 + ts_xxxyz_yyyzz[i] * pa_x[i];

        ts_xxxxyz_yyzzz[i] = 3.0 * ts_xxyz_yyzzz[i] * fe_0 + ts_xxxyz_yyzzz[i] * pa_x[i];

        ts_xxxxyz_yzzzz[i] = 3.0 * ts_xxyz_yzzzz[i] * fe_0 + ts_xxxyz_yzzzz[i] * pa_x[i];

        ts_xxxxyz_zzzzz[i] = ts_xxxxz_zzzzz[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : IH

    auto ts_xxxxzz_xxxxx = pbuffer.data(idx_ovl_ih + 105);

    auto ts_xxxxzz_xxxxy = pbuffer.data(idx_ovl_ih + 106);

    auto ts_xxxxzz_xxxxz = pbuffer.data(idx_ovl_ih + 107);

    auto ts_xxxxzz_xxxyy = pbuffer.data(idx_ovl_ih + 108);

    auto ts_xxxxzz_xxxyz = pbuffer.data(idx_ovl_ih + 109);

    auto ts_xxxxzz_xxxzz = pbuffer.data(idx_ovl_ih + 110);

    auto ts_xxxxzz_xxyyy = pbuffer.data(idx_ovl_ih + 111);

    auto ts_xxxxzz_xxyyz = pbuffer.data(idx_ovl_ih + 112);

    auto ts_xxxxzz_xxyzz = pbuffer.data(idx_ovl_ih + 113);

    auto ts_xxxxzz_xxzzz = pbuffer.data(idx_ovl_ih + 114);

    auto ts_xxxxzz_xyyyy = pbuffer.data(idx_ovl_ih + 115);

    auto ts_xxxxzz_xyyyz = pbuffer.data(idx_ovl_ih + 116);

    auto ts_xxxxzz_xyyzz = pbuffer.data(idx_ovl_ih + 117);

    auto ts_xxxxzz_xyzzz = pbuffer.data(idx_ovl_ih + 118);

    auto ts_xxxxzz_xzzzz = pbuffer.data(idx_ovl_ih + 119);

    auto ts_xxxxzz_yyyyy = pbuffer.data(idx_ovl_ih + 120);

    auto ts_xxxxzz_yyyyz = pbuffer.data(idx_ovl_ih + 121);

    auto ts_xxxxzz_yyyzz = pbuffer.data(idx_ovl_ih + 122);

    auto ts_xxxxzz_yyzzz = pbuffer.data(idx_ovl_ih + 123);

    auto ts_xxxxzz_yzzzz = pbuffer.data(idx_ovl_ih + 124);

    auto ts_xxxxzz_zzzzz = pbuffer.data(idx_ovl_ih + 125);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxyy, ts_xxxx_xxyyy, ts_xxxx_xyyyy, ts_xxxxz_xxxxx, ts_xxxxz_xxxxy, ts_xxxxz_xxxyy, ts_xxxxz_xxyyy, ts_xxxxz_xyyyy, ts_xxxxzz_xxxxx, ts_xxxxzz_xxxxy, ts_xxxxzz_xxxxz, ts_xxxxzz_xxxyy, ts_xxxxzz_xxxyz, ts_xxxxzz_xxxzz, ts_xxxxzz_xxyyy, ts_xxxxzz_xxyyz, ts_xxxxzz_xxyzz, ts_xxxxzz_xxzzz, ts_xxxxzz_xyyyy, ts_xxxxzz_xyyyz, ts_xxxxzz_xyyzz, ts_xxxxzz_xyzzz, ts_xxxxzz_xzzzz, ts_xxxxzz_yyyyy, ts_xxxxzz_yyyyz, ts_xxxxzz_yyyzz, ts_xxxxzz_yyzzz, ts_xxxxzz_yzzzz, ts_xxxxzz_zzzzz, ts_xxxzz_xxxxz, ts_xxxzz_xxxyz, ts_xxxzz_xxxz, ts_xxxzz_xxxzz, ts_xxxzz_xxyyz, ts_xxxzz_xxyz, ts_xxxzz_xxyzz, ts_xxxzz_xxzz, ts_xxxzz_xxzzz, ts_xxxzz_xyyyz, ts_xxxzz_xyyz, ts_xxxzz_xyyzz, ts_xxxzz_xyzz, ts_xxxzz_xyzzz, ts_xxxzz_xzzz, ts_xxxzz_xzzzz, ts_xxxzz_yyyyy, ts_xxxzz_yyyyz, ts_xxxzz_yyyz, ts_xxxzz_yyyzz, ts_xxxzz_yyzz, ts_xxxzz_yyzzz, ts_xxxzz_yzzz, ts_xxxzz_yzzzz, ts_xxxzz_zzzz, ts_xxxzz_zzzzz, ts_xxzz_xxxxz, ts_xxzz_xxxyz, ts_xxzz_xxxzz, ts_xxzz_xxyyz, ts_xxzz_xxyzz, ts_xxzz_xxzzz, ts_xxzz_xyyyz, ts_xxzz_xyyzz, ts_xxzz_xyzzz, ts_xxzz_xzzzz, ts_xxzz_yyyyy, ts_xxzz_yyyyz, ts_xxzz_yyyzz, ts_xxzz_yyzzz, ts_xxzz_yzzzz, ts_xxzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxzz_xxxxx[i] = ts_xxxx_xxxxx[i] * fe_0 + ts_xxxxz_xxxxx[i] * pa_z[i];

        ts_xxxxzz_xxxxy[i] = ts_xxxx_xxxxy[i] * fe_0 + ts_xxxxz_xxxxy[i] * pa_z[i];

        ts_xxxxzz_xxxxz[i] = 3.0 * ts_xxzz_xxxxz[i] * fe_0 + 4.0 * ts_xxxzz_xxxz[i] * fe_0 + ts_xxxzz_xxxxz[i] * pa_x[i];

        ts_xxxxzz_xxxyy[i] = ts_xxxx_xxxyy[i] * fe_0 + ts_xxxxz_xxxyy[i] * pa_z[i];

        ts_xxxxzz_xxxyz[i] = 3.0 * ts_xxzz_xxxyz[i] * fe_0 + 3.0 * ts_xxxzz_xxyz[i] * fe_0 + ts_xxxzz_xxxyz[i] * pa_x[i];

        ts_xxxxzz_xxxzz[i] = 3.0 * ts_xxzz_xxxzz[i] * fe_0 + 3.0 * ts_xxxzz_xxzz[i] * fe_0 + ts_xxxzz_xxxzz[i] * pa_x[i];

        ts_xxxxzz_xxyyy[i] = ts_xxxx_xxyyy[i] * fe_0 + ts_xxxxz_xxyyy[i] * pa_z[i];

        ts_xxxxzz_xxyyz[i] = 3.0 * ts_xxzz_xxyyz[i] * fe_0 + 2.0 * ts_xxxzz_xyyz[i] * fe_0 + ts_xxxzz_xxyyz[i] * pa_x[i];

        ts_xxxxzz_xxyzz[i] = 3.0 * ts_xxzz_xxyzz[i] * fe_0 + 2.0 * ts_xxxzz_xyzz[i] * fe_0 + ts_xxxzz_xxyzz[i] * pa_x[i];

        ts_xxxxzz_xxzzz[i] = 3.0 * ts_xxzz_xxzzz[i] * fe_0 + 2.0 * ts_xxxzz_xzzz[i] * fe_0 + ts_xxxzz_xxzzz[i] * pa_x[i];

        ts_xxxxzz_xyyyy[i] = ts_xxxx_xyyyy[i] * fe_0 + ts_xxxxz_xyyyy[i] * pa_z[i];

        ts_xxxxzz_xyyyz[i] = 3.0 * ts_xxzz_xyyyz[i] * fe_0 + ts_xxxzz_yyyz[i] * fe_0 + ts_xxxzz_xyyyz[i] * pa_x[i];

        ts_xxxxzz_xyyzz[i] = 3.0 * ts_xxzz_xyyzz[i] * fe_0 + ts_xxxzz_yyzz[i] * fe_0 + ts_xxxzz_xyyzz[i] * pa_x[i];

        ts_xxxxzz_xyzzz[i] = 3.0 * ts_xxzz_xyzzz[i] * fe_0 + ts_xxxzz_yzzz[i] * fe_0 + ts_xxxzz_xyzzz[i] * pa_x[i];

        ts_xxxxzz_xzzzz[i] = 3.0 * ts_xxzz_xzzzz[i] * fe_0 + ts_xxxzz_zzzz[i] * fe_0 + ts_xxxzz_xzzzz[i] * pa_x[i];

        ts_xxxxzz_yyyyy[i] = 3.0 * ts_xxzz_yyyyy[i] * fe_0 + ts_xxxzz_yyyyy[i] * pa_x[i];

        ts_xxxxzz_yyyyz[i] = 3.0 * ts_xxzz_yyyyz[i] * fe_0 + ts_xxxzz_yyyyz[i] * pa_x[i];

        ts_xxxxzz_yyyzz[i] = 3.0 * ts_xxzz_yyyzz[i] * fe_0 + ts_xxxzz_yyyzz[i] * pa_x[i];

        ts_xxxxzz_yyzzz[i] = 3.0 * ts_xxzz_yyzzz[i] * fe_0 + ts_xxxzz_yyzzz[i] * pa_x[i];

        ts_xxxxzz_yzzzz[i] = 3.0 * ts_xxzz_yzzzz[i] * fe_0 + ts_xxxzz_yzzzz[i] * pa_x[i];

        ts_xxxxzz_zzzzz[i] = 3.0 * ts_xxzz_zzzzz[i] * fe_0 + ts_xxxzz_zzzzz[i] * pa_x[i];
    }

    // Set up 126-147 components of targeted buffer : IH

    auto ts_xxxyyy_xxxxx = pbuffer.data(idx_ovl_ih + 126);

    auto ts_xxxyyy_xxxxy = pbuffer.data(idx_ovl_ih + 127);

    auto ts_xxxyyy_xxxxz = pbuffer.data(idx_ovl_ih + 128);

    auto ts_xxxyyy_xxxyy = pbuffer.data(idx_ovl_ih + 129);

    auto ts_xxxyyy_xxxyz = pbuffer.data(idx_ovl_ih + 130);

    auto ts_xxxyyy_xxxzz = pbuffer.data(idx_ovl_ih + 131);

    auto ts_xxxyyy_xxyyy = pbuffer.data(idx_ovl_ih + 132);

    auto ts_xxxyyy_xxyyz = pbuffer.data(idx_ovl_ih + 133);

    auto ts_xxxyyy_xxyzz = pbuffer.data(idx_ovl_ih + 134);

    auto ts_xxxyyy_xxzzz = pbuffer.data(idx_ovl_ih + 135);

    auto ts_xxxyyy_xyyyy = pbuffer.data(idx_ovl_ih + 136);

    auto ts_xxxyyy_xyyyz = pbuffer.data(idx_ovl_ih + 137);

    auto ts_xxxyyy_xyyzz = pbuffer.data(idx_ovl_ih + 138);

    auto ts_xxxyyy_xyzzz = pbuffer.data(idx_ovl_ih + 139);

    auto ts_xxxyyy_xzzzz = pbuffer.data(idx_ovl_ih + 140);

    auto ts_xxxyyy_yyyyy = pbuffer.data(idx_ovl_ih + 141);

    auto ts_xxxyyy_yyyyz = pbuffer.data(idx_ovl_ih + 142);

    auto ts_xxxyyy_yyyzz = pbuffer.data(idx_ovl_ih + 143);

    auto ts_xxxyyy_yyzzz = pbuffer.data(idx_ovl_ih + 144);

    auto ts_xxxyyy_yzzzz = pbuffer.data(idx_ovl_ih + 145);

    auto ts_xxxyyy_zzzzz = pbuffer.data(idx_ovl_ih + 146);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxy_xxxxx, ts_xxxy_xxxxz, ts_xxxy_xxxzz, ts_xxxy_xxzzz, ts_xxxy_xzzzz, ts_xxxyy_xxxxx, ts_xxxyy_xxxxz, ts_xxxyy_xxxzz, ts_xxxyy_xxzzz, ts_xxxyy_xzzzz, ts_xxxyyy_xxxxx, ts_xxxyyy_xxxxy, ts_xxxyyy_xxxxz, ts_xxxyyy_xxxyy, ts_xxxyyy_xxxyz, ts_xxxyyy_xxxzz, ts_xxxyyy_xxyyy, ts_xxxyyy_xxyyz, ts_xxxyyy_xxyzz, ts_xxxyyy_xxzzz, ts_xxxyyy_xyyyy, ts_xxxyyy_xyyyz, ts_xxxyyy_xyyzz, ts_xxxyyy_xyzzz, ts_xxxyyy_xzzzz, ts_xxxyyy_yyyyy, ts_xxxyyy_yyyyz, ts_xxxyyy_yyyzz, ts_xxxyyy_yyzzz, ts_xxxyyy_yzzzz, ts_xxxyyy_zzzzz, ts_xxyyy_xxxxy, ts_xxyyy_xxxy, ts_xxyyy_xxxyy, ts_xxyyy_xxxyz, ts_xxyyy_xxyy, ts_xxyyy_xxyyy, ts_xxyyy_xxyyz, ts_xxyyy_xxyz, ts_xxyyy_xxyzz, ts_xxyyy_xyyy, ts_xxyyy_xyyyy, ts_xxyyy_xyyyz, ts_xxyyy_xyyz, ts_xxyyy_xyyzz, ts_xxyyy_xyzz, ts_xxyyy_xyzzz, ts_xxyyy_yyyy, ts_xxyyy_yyyyy, ts_xxyyy_yyyyz, ts_xxyyy_yyyz, ts_xxyyy_yyyzz, ts_xxyyy_yyzz, ts_xxyyy_yyzzz, ts_xxyyy_yzzz, ts_xxyyy_yzzzz, ts_xxyyy_zzzzz, ts_xyyy_xxxxy, ts_xyyy_xxxyy, ts_xyyy_xxxyz, ts_xyyy_xxyyy, ts_xyyy_xxyyz, ts_xyyy_xxyzz, ts_xyyy_xyyyy, ts_xyyy_xyyyz, ts_xyyy_xyyzz, ts_xyyy_xyzzz, ts_xyyy_yyyyy, ts_xyyy_yyyyz, ts_xyyy_yyyzz, ts_xyyy_yyzzz, ts_xyyy_yzzzz, ts_xyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyyy_xxxxx[i] = 2.0 * ts_xxxy_xxxxx[i] * fe_0 + ts_xxxyy_xxxxx[i] * pa_y[i];

        ts_xxxyyy_xxxxy[i] = 2.0 * ts_xyyy_xxxxy[i] * fe_0 + 4.0 * ts_xxyyy_xxxy[i] * fe_0 + ts_xxyyy_xxxxy[i] * pa_x[i];

        ts_xxxyyy_xxxxz[i] = 2.0 * ts_xxxy_xxxxz[i] * fe_0 + ts_xxxyy_xxxxz[i] * pa_y[i];

        ts_xxxyyy_xxxyy[i] = 2.0 * ts_xyyy_xxxyy[i] * fe_0 + 3.0 * ts_xxyyy_xxyy[i] * fe_0 + ts_xxyyy_xxxyy[i] * pa_x[i];

        ts_xxxyyy_xxxyz[i] = 2.0 * ts_xyyy_xxxyz[i] * fe_0 + 3.0 * ts_xxyyy_xxyz[i] * fe_0 + ts_xxyyy_xxxyz[i] * pa_x[i];

        ts_xxxyyy_xxxzz[i] = 2.0 * ts_xxxy_xxxzz[i] * fe_0 + ts_xxxyy_xxxzz[i] * pa_y[i];

        ts_xxxyyy_xxyyy[i] = 2.0 * ts_xyyy_xxyyy[i] * fe_0 + 2.0 * ts_xxyyy_xyyy[i] * fe_0 + ts_xxyyy_xxyyy[i] * pa_x[i];

        ts_xxxyyy_xxyyz[i] = 2.0 * ts_xyyy_xxyyz[i] * fe_0 + 2.0 * ts_xxyyy_xyyz[i] * fe_0 + ts_xxyyy_xxyyz[i] * pa_x[i];

        ts_xxxyyy_xxyzz[i] = 2.0 * ts_xyyy_xxyzz[i] * fe_0 + 2.0 * ts_xxyyy_xyzz[i] * fe_0 + ts_xxyyy_xxyzz[i] * pa_x[i];

        ts_xxxyyy_xxzzz[i] = 2.0 * ts_xxxy_xxzzz[i] * fe_0 + ts_xxxyy_xxzzz[i] * pa_y[i];

        ts_xxxyyy_xyyyy[i] = 2.0 * ts_xyyy_xyyyy[i] * fe_0 + ts_xxyyy_yyyy[i] * fe_0 + ts_xxyyy_xyyyy[i] * pa_x[i];

        ts_xxxyyy_xyyyz[i] = 2.0 * ts_xyyy_xyyyz[i] * fe_0 + ts_xxyyy_yyyz[i] * fe_0 + ts_xxyyy_xyyyz[i] * pa_x[i];

        ts_xxxyyy_xyyzz[i] = 2.0 * ts_xyyy_xyyzz[i] * fe_0 + ts_xxyyy_yyzz[i] * fe_0 + ts_xxyyy_xyyzz[i] * pa_x[i];

        ts_xxxyyy_xyzzz[i] = 2.0 * ts_xyyy_xyzzz[i] * fe_0 + ts_xxyyy_yzzz[i] * fe_0 + ts_xxyyy_xyzzz[i] * pa_x[i];

        ts_xxxyyy_xzzzz[i] = 2.0 * ts_xxxy_xzzzz[i] * fe_0 + ts_xxxyy_xzzzz[i] * pa_y[i];

        ts_xxxyyy_yyyyy[i] = 2.0 * ts_xyyy_yyyyy[i] * fe_0 + ts_xxyyy_yyyyy[i] * pa_x[i];

        ts_xxxyyy_yyyyz[i] = 2.0 * ts_xyyy_yyyyz[i] * fe_0 + ts_xxyyy_yyyyz[i] * pa_x[i];

        ts_xxxyyy_yyyzz[i] = 2.0 * ts_xyyy_yyyzz[i] * fe_0 + ts_xxyyy_yyyzz[i] * pa_x[i];

        ts_xxxyyy_yyzzz[i] = 2.0 * ts_xyyy_yyzzz[i] * fe_0 + ts_xxyyy_yyzzz[i] * pa_x[i];

        ts_xxxyyy_yzzzz[i] = 2.0 * ts_xyyy_yzzzz[i] * fe_0 + ts_xxyyy_yzzzz[i] * pa_x[i];

        ts_xxxyyy_zzzzz[i] = 2.0 * ts_xyyy_zzzzz[i] * fe_0 + ts_xxyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 147-168 components of targeted buffer : IH

    auto ts_xxxyyz_xxxxx = pbuffer.data(idx_ovl_ih + 147);

    auto ts_xxxyyz_xxxxy = pbuffer.data(idx_ovl_ih + 148);

    auto ts_xxxyyz_xxxxz = pbuffer.data(idx_ovl_ih + 149);

    auto ts_xxxyyz_xxxyy = pbuffer.data(idx_ovl_ih + 150);

    auto ts_xxxyyz_xxxyz = pbuffer.data(idx_ovl_ih + 151);

    auto ts_xxxyyz_xxxzz = pbuffer.data(idx_ovl_ih + 152);

    auto ts_xxxyyz_xxyyy = pbuffer.data(idx_ovl_ih + 153);

    auto ts_xxxyyz_xxyyz = pbuffer.data(idx_ovl_ih + 154);

    auto ts_xxxyyz_xxyzz = pbuffer.data(idx_ovl_ih + 155);

    auto ts_xxxyyz_xxzzz = pbuffer.data(idx_ovl_ih + 156);

    auto ts_xxxyyz_xyyyy = pbuffer.data(idx_ovl_ih + 157);

    auto ts_xxxyyz_xyyyz = pbuffer.data(idx_ovl_ih + 158);

    auto ts_xxxyyz_xyyzz = pbuffer.data(idx_ovl_ih + 159);

    auto ts_xxxyyz_xyzzz = pbuffer.data(idx_ovl_ih + 160);

    auto ts_xxxyyz_xzzzz = pbuffer.data(idx_ovl_ih + 161);

    auto ts_xxxyyz_yyyyy = pbuffer.data(idx_ovl_ih + 162);

    auto ts_xxxyyz_yyyyz = pbuffer.data(idx_ovl_ih + 163);

    auto ts_xxxyyz_yyyzz = pbuffer.data(idx_ovl_ih + 164);

    auto ts_xxxyyz_yyzzz = pbuffer.data(idx_ovl_ih + 165);

    auto ts_xxxyyz_yzzzz = pbuffer.data(idx_ovl_ih + 166);

    auto ts_xxxyyz_zzzzz = pbuffer.data(idx_ovl_ih + 167);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxyy_xxxxx, ts_xxxyy_xxxxy, ts_xxxyy_xxxy, ts_xxxyy_xxxyy, ts_xxxyy_xxxyz, ts_xxxyy_xxyy, ts_xxxyy_xxyyy, ts_xxxyy_xxyyz, ts_xxxyy_xxyz, ts_xxxyy_xxyzz, ts_xxxyy_xyyy, ts_xxxyy_xyyyy, ts_xxxyy_xyyyz, ts_xxxyy_xyyz, ts_xxxyy_xyyzz, ts_xxxyy_xyzz, ts_xxxyy_xyzzz, ts_xxxyy_yyyyy, ts_xxxyyz_xxxxx, ts_xxxyyz_xxxxy, ts_xxxyyz_xxxxz, ts_xxxyyz_xxxyy, ts_xxxyyz_xxxyz, ts_xxxyyz_xxxzz, ts_xxxyyz_xxyyy, ts_xxxyyz_xxyyz, ts_xxxyyz_xxyzz, ts_xxxyyz_xxzzz, ts_xxxyyz_xyyyy, ts_xxxyyz_xyyyz, ts_xxxyyz_xyyzz, ts_xxxyyz_xyzzz, ts_xxxyyz_xzzzz, ts_xxxyyz_yyyyy, ts_xxxyyz_yyyyz, ts_xxxyyz_yyyzz, ts_xxxyyz_yyzzz, ts_xxxyyz_yzzzz, ts_xxxyyz_zzzzz, ts_xxxyz_xxxxz, ts_xxxyz_xxxzz, ts_xxxyz_xxzzz, ts_xxxyz_xzzzz, ts_xxxz_xxxxz, ts_xxxz_xxxzz, ts_xxxz_xxzzz, ts_xxxz_xzzzz, ts_xxyyz_yyyyz, ts_xxyyz_yyyzz, ts_xxyyz_yyzzz, ts_xxyyz_yzzzz, ts_xxyyz_zzzzz, ts_xyyz_yyyyz, ts_xyyz_yyyzz, ts_xyyz_yyzzz, ts_xyyz_yzzzz, ts_xyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyyz_xxxxx[i] = ts_xxxyy_xxxxx[i] * pa_z[i];

        ts_xxxyyz_xxxxy[i] = ts_xxxyy_xxxxy[i] * pa_z[i];

        ts_xxxyyz_xxxxz[i] = ts_xxxz_xxxxz[i] * fe_0 + ts_xxxyz_xxxxz[i] * pa_y[i];

        ts_xxxyyz_xxxyy[i] = ts_xxxyy_xxxyy[i] * pa_z[i];

        ts_xxxyyz_xxxyz[i] = ts_xxxyy_xxxy[i] * fe_0 + ts_xxxyy_xxxyz[i] * pa_z[i];

        ts_xxxyyz_xxxzz[i] = ts_xxxz_xxxzz[i] * fe_0 + ts_xxxyz_xxxzz[i] * pa_y[i];

        ts_xxxyyz_xxyyy[i] = ts_xxxyy_xxyyy[i] * pa_z[i];

        ts_xxxyyz_xxyyz[i] = ts_xxxyy_xxyy[i] * fe_0 + ts_xxxyy_xxyyz[i] * pa_z[i];

        ts_xxxyyz_xxyzz[i] = 2.0 * ts_xxxyy_xxyz[i] * fe_0 + ts_xxxyy_xxyzz[i] * pa_z[i];

        ts_xxxyyz_xxzzz[i] = ts_xxxz_xxzzz[i] * fe_0 + ts_xxxyz_xxzzz[i] * pa_y[i];

        ts_xxxyyz_xyyyy[i] = ts_xxxyy_xyyyy[i] * pa_z[i];

        ts_xxxyyz_xyyyz[i] = ts_xxxyy_xyyy[i] * fe_0 + ts_xxxyy_xyyyz[i] * pa_z[i];

        ts_xxxyyz_xyyzz[i] = 2.0 * ts_xxxyy_xyyz[i] * fe_0 + ts_xxxyy_xyyzz[i] * pa_z[i];

        ts_xxxyyz_xyzzz[i] = 3.0 * ts_xxxyy_xyzz[i] * fe_0 + ts_xxxyy_xyzzz[i] * pa_z[i];

        ts_xxxyyz_xzzzz[i] = ts_xxxz_xzzzz[i] * fe_0 + ts_xxxyz_xzzzz[i] * pa_y[i];

        ts_xxxyyz_yyyyy[i] = ts_xxxyy_yyyyy[i] * pa_z[i];

        ts_xxxyyz_yyyyz[i] = 2.0 * ts_xyyz_yyyyz[i] * fe_0 + ts_xxyyz_yyyyz[i] * pa_x[i];

        ts_xxxyyz_yyyzz[i] = 2.0 * ts_xyyz_yyyzz[i] * fe_0 + ts_xxyyz_yyyzz[i] * pa_x[i];

        ts_xxxyyz_yyzzz[i] = 2.0 * ts_xyyz_yyzzz[i] * fe_0 + ts_xxyyz_yyzzz[i] * pa_x[i];

        ts_xxxyyz_yzzzz[i] = 2.0 * ts_xyyz_yzzzz[i] * fe_0 + ts_xxyyz_yzzzz[i] * pa_x[i];

        ts_xxxyyz_zzzzz[i] = 2.0 * ts_xyyz_zzzzz[i] * fe_0 + ts_xxyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 168-189 components of targeted buffer : IH

    auto ts_xxxyzz_xxxxx = pbuffer.data(idx_ovl_ih + 168);

    auto ts_xxxyzz_xxxxy = pbuffer.data(idx_ovl_ih + 169);

    auto ts_xxxyzz_xxxxz = pbuffer.data(idx_ovl_ih + 170);

    auto ts_xxxyzz_xxxyy = pbuffer.data(idx_ovl_ih + 171);

    auto ts_xxxyzz_xxxyz = pbuffer.data(idx_ovl_ih + 172);

    auto ts_xxxyzz_xxxzz = pbuffer.data(idx_ovl_ih + 173);

    auto ts_xxxyzz_xxyyy = pbuffer.data(idx_ovl_ih + 174);

    auto ts_xxxyzz_xxyyz = pbuffer.data(idx_ovl_ih + 175);

    auto ts_xxxyzz_xxyzz = pbuffer.data(idx_ovl_ih + 176);

    auto ts_xxxyzz_xxzzz = pbuffer.data(idx_ovl_ih + 177);

    auto ts_xxxyzz_xyyyy = pbuffer.data(idx_ovl_ih + 178);

    auto ts_xxxyzz_xyyyz = pbuffer.data(idx_ovl_ih + 179);

    auto ts_xxxyzz_xyyzz = pbuffer.data(idx_ovl_ih + 180);

    auto ts_xxxyzz_xyzzz = pbuffer.data(idx_ovl_ih + 181);

    auto ts_xxxyzz_xzzzz = pbuffer.data(idx_ovl_ih + 182);

    auto ts_xxxyzz_yyyyy = pbuffer.data(idx_ovl_ih + 183);

    auto ts_xxxyzz_yyyyz = pbuffer.data(idx_ovl_ih + 184);

    auto ts_xxxyzz_yyyzz = pbuffer.data(idx_ovl_ih + 185);

    auto ts_xxxyzz_yyzzz = pbuffer.data(idx_ovl_ih + 186);

    auto ts_xxxyzz_yzzzz = pbuffer.data(idx_ovl_ih + 187);

    auto ts_xxxyzz_zzzzz = pbuffer.data(idx_ovl_ih + 188);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxyzz_xxxxx, ts_xxxyzz_xxxxy, ts_xxxyzz_xxxxz, ts_xxxyzz_xxxyy, ts_xxxyzz_xxxyz, ts_xxxyzz_xxxzz, ts_xxxyzz_xxyyy, ts_xxxyzz_xxyyz, ts_xxxyzz_xxyzz, ts_xxxyzz_xxzzz, ts_xxxyzz_xyyyy, ts_xxxyzz_xyyyz, ts_xxxyzz_xyyzz, ts_xxxyzz_xyzzz, ts_xxxyzz_xzzzz, ts_xxxyzz_yyyyy, ts_xxxyzz_yyyyz, ts_xxxyzz_yyyzz, ts_xxxyzz_yyzzz, ts_xxxyzz_yzzzz, ts_xxxyzz_zzzzz, ts_xxxzz_xxxx, ts_xxxzz_xxxxx, ts_xxxzz_xxxxy, ts_xxxzz_xxxxz, ts_xxxzz_xxxy, ts_xxxzz_xxxyy, ts_xxxzz_xxxyz, ts_xxxzz_xxxz, ts_xxxzz_xxxzz, ts_xxxzz_xxyy, ts_xxxzz_xxyyy, ts_xxxzz_xxyyz, ts_xxxzz_xxyz, ts_xxxzz_xxyzz, ts_xxxzz_xxzz, ts_xxxzz_xxzzz, ts_xxxzz_xyyy, ts_xxxzz_xyyyy, ts_xxxzz_xyyyz, ts_xxxzz_xyyz, ts_xxxzz_xyyzz, ts_xxxzz_xyzz, ts_xxxzz_xyzzz, ts_xxxzz_xzzz, ts_xxxzz_xzzzz, ts_xxxzz_zzzzz, ts_xxyzz_yyyyy, ts_xxyzz_yyyyz, ts_xxyzz_yyyzz, ts_xxyzz_yyzzz, ts_xxyzz_yzzzz, ts_xyzz_yyyyy, ts_xyzz_yyyyz, ts_xyzz_yyyzz, ts_xyzz_yyzzz, ts_xyzz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyzz_xxxxx[i] = ts_xxxzz_xxxxx[i] * pa_y[i];

        ts_xxxyzz_xxxxy[i] = ts_xxxzz_xxxx[i] * fe_0 + ts_xxxzz_xxxxy[i] * pa_y[i];

        ts_xxxyzz_xxxxz[i] = ts_xxxzz_xxxxz[i] * pa_y[i];

        ts_xxxyzz_xxxyy[i] = 2.0 * ts_xxxzz_xxxy[i] * fe_0 + ts_xxxzz_xxxyy[i] * pa_y[i];

        ts_xxxyzz_xxxyz[i] = ts_xxxzz_xxxz[i] * fe_0 + ts_xxxzz_xxxyz[i] * pa_y[i];

        ts_xxxyzz_xxxzz[i] = ts_xxxzz_xxxzz[i] * pa_y[i];

        ts_xxxyzz_xxyyy[i] = 3.0 * ts_xxxzz_xxyy[i] * fe_0 + ts_xxxzz_xxyyy[i] * pa_y[i];

        ts_xxxyzz_xxyyz[i] = 2.0 * ts_xxxzz_xxyz[i] * fe_0 + ts_xxxzz_xxyyz[i] * pa_y[i];

        ts_xxxyzz_xxyzz[i] = ts_xxxzz_xxzz[i] * fe_0 + ts_xxxzz_xxyzz[i] * pa_y[i];

        ts_xxxyzz_xxzzz[i] = ts_xxxzz_xxzzz[i] * pa_y[i];

        ts_xxxyzz_xyyyy[i] = 4.0 * ts_xxxzz_xyyy[i] * fe_0 + ts_xxxzz_xyyyy[i] * pa_y[i];

        ts_xxxyzz_xyyyz[i] = 3.0 * ts_xxxzz_xyyz[i] * fe_0 + ts_xxxzz_xyyyz[i] * pa_y[i];

        ts_xxxyzz_xyyzz[i] = 2.0 * ts_xxxzz_xyzz[i] * fe_0 + ts_xxxzz_xyyzz[i] * pa_y[i];

        ts_xxxyzz_xyzzz[i] = ts_xxxzz_xzzz[i] * fe_0 + ts_xxxzz_xyzzz[i] * pa_y[i];

        ts_xxxyzz_xzzzz[i] = ts_xxxzz_xzzzz[i] * pa_y[i];

        ts_xxxyzz_yyyyy[i] = 2.0 * ts_xyzz_yyyyy[i] * fe_0 + ts_xxyzz_yyyyy[i] * pa_x[i];

        ts_xxxyzz_yyyyz[i] = 2.0 * ts_xyzz_yyyyz[i] * fe_0 + ts_xxyzz_yyyyz[i] * pa_x[i];

        ts_xxxyzz_yyyzz[i] = 2.0 * ts_xyzz_yyyzz[i] * fe_0 + ts_xxyzz_yyyzz[i] * pa_x[i];

        ts_xxxyzz_yyzzz[i] = 2.0 * ts_xyzz_yyzzz[i] * fe_0 + ts_xxyzz_yyzzz[i] * pa_x[i];

        ts_xxxyzz_yzzzz[i] = 2.0 * ts_xyzz_yzzzz[i] * fe_0 + ts_xxyzz_yzzzz[i] * pa_x[i];

        ts_xxxyzz_zzzzz[i] = ts_xxxzz_zzzzz[i] * pa_y[i];
    }

    // Set up 189-210 components of targeted buffer : IH

    auto ts_xxxzzz_xxxxx = pbuffer.data(idx_ovl_ih + 189);

    auto ts_xxxzzz_xxxxy = pbuffer.data(idx_ovl_ih + 190);

    auto ts_xxxzzz_xxxxz = pbuffer.data(idx_ovl_ih + 191);

    auto ts_xxxzzz_xxxyy = pbuffer.data(idx_ovl_ih + 192);

    auto ts_xxxzzz_xxxyz = pbuffer.data(idx_ovl_ih + 193);

    auto ts_xxxzzz_xxxzz = pbuffer.data(idx_ovl_ih + 194);

    auto ts_xxxzzz_xxyyy = pbuffer.data(idx_ovl_ih + 195);

    auto ts_xxxzzz_xxyyz = pbuffer.data(idx_ovl_ih + 196);

    auto ts_xxxzzz_xxyzz = pbuffer.data(idx_ovl_ih + 197);

    auto ts_xxxzzz_xxzzz = pbuffer.data(idx_ovl_ih + 198);

    auto ts_xxxzzz_xyyyy = pbuffer.data(idx_ovl_ih + 199);

    auto ts_xxxzzz_xyyyz = pbuffer.data(idx_ovl_ih + 200);

    auto ts_xxxzzz_xyyzz = pbuffer.data(idx_ovl_ih + 201);

    auto ts_xxxzzz_xyzzz = pbuffer.data(idx_ovl_ih + 202);

    auto ts_xxxzzz_xzzzz = pbuffer.data(idx_ovl_ih + 203);

    auto ts_xxxzzz_yyyyy = pbuffer.data(idx_ovl_ih + 204);

    auto ts_xxxzzz_yyyyz = pbuffer.data(idx_ovl_ih + 205);

    auto ts_xxxzzz_yyyzz = pbuffer.data(idx_ovl_ih + 206);

    auto ts_xxxzzz_yyzzz = pbuffer.data(idx_ovl_ih + 207);

    auto ts_xxxzzz_yzzzz = pbuffer.data(idx_ovl_ih + 208);

    auto ts_xxxzzz_zzzzz = pbuffer.data(idx_ovl_ih + 209);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxz_xxxxx, ts_xxxz_xxxxy, ts_xxxz_xxxyy, ts_xxxz_xxyyy, ts_xxxz_xyyyy, ts_xxxzz_xxxxx, ts_xxxzz_xxxxy, ts_xxxzz_xxxyy, ts_xxxzz_xxyyy, ts_xxxzz_xyyyy, ts_xxxzzz_xxxxx, ts_xxxzzz_xxxxy, ts_xxxzzz_xxxxz, ts_xxxzzz_xxxyy, ts_xxxzzz_xxxyz, ts_xxxzzz_xxxzz, ts_xxxzzz_xxyyy, ts_xxxzzz_xxyyz, ts_xxxzzz_xxyzz, ts_xxxzzz_xxzzz, ts_xxxzzz_xyyyy, ts_xxxzzz_xyyyz, ts_xxxzzz_xyyzz, ts_xxxzzz_xyzzz, ts_xxxzzz_xzzzz, ts_xxxzzz_yyyyy, ts_xxxzzz_yyyyz, ts_xxxzzz_yyyzz, ts_xxxzzz_yyzzz, ts_xxxzzz_yzzzz, ts_xxxzzz_zzzzz, ts_xxzzz_xxxxz, ts_xxzzz_xxxyz, ts_xxzzz_xxxz, ts_xxzzz_xxxzz, ts_xxzzz_xxyyz, ts_xxzzz_xxyz, ts_xxzzz_xxyzz, ts_xxzzz_xxzz, ts_xxzzz_xxzzz, ts_xxzzz_xyyyz, ts_xxzzz_xyyz, ts_xxzzz_xyyzz, ts_xxzzz_xyzz, ts_xxzzz_xyzzz, ts_xxzzz_xzzz, ts_xxzzz_xzzzz, ts_xxzzz_yyyyy, ts_xxzzz_yyyyz, ts_xxzzz_yyyz, ts_xxzzz_yyyzz, ts_xxzzz_yyzz, ts_xxzzz_yyzzz, ts_xxzzz_yzzz, ts_xxzzz_yzzzz, ts_xxzzz_zzzz, ts_xxzzz_zzzzz, ts_xzzz_xxxxz, ts_xzzz_xxxyz, ts_xzzz_xxxzz, ts_xzzz_xxyyz, ts_xzzz_xxyzz, ts_xzzz_xxzzz, ts_xzzz_xyyyz, ts_xzzz_xyyzz, ts_xzzz_xyzzz, ts_xzzz_xzzzz, ts_xzzz_yyyyy, ts_xzzz_yyyyz, ts_xzzz_yyyzz, ts_xzzz_yyzzz, ts_xzzz_yzzzz, ts_xzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxzzz_xxxxx[i] = 2.0 * ts_xxxz_xxxxx[i] * fe_0 + ts_xxxzz_xxxxx[i] * pa_z[i];

        ts_xxxzzz_xxxxy[i] = 2.0 * ts_xxxz_xxxxy[i] * fe_0 + ts_xxxzz_xxxxy[i] * pa_z[i];

        ts_xxxzzz_xxxxz[i] = 2.0 * ts_xzzz_xxxxz[i] * fe_0 + 4.0 * ts_xxzzz_xxxz[i] * fe_0 + ts_xxzzz_xxxxz[i] * pa_x[i];

        ts_xxxzzz_xxxyy[i] = 2.0 * ts_xxxz_xxxyy[i] * fe_0 + ts_xxxzz_xxxyy[i] * pa_z[i];

        ts_xxxzzz_xxxyz[i] = 2.0 * ts_xzzz_xxxyz[i] * fe_0 + 3.0 * ts_xxzzz_xxyz[i] * fe_0 + ts_xxzzz_xxxyz[i] * pa_x[i];

        ts_xxxzzz_xxxzz[i] = 2.0 * ts_xzzz_xxxzz[i] * fe_0 + 3.0 * ts_xxzzz_xxzz[i] * fe_0 + ts_xxzzz_xxxzz[i] * pa_x[i];

        ts_xxxzzz_xxyyy[i] = 2.0 * ts_xxxz_xxyyy[i] * fe_0 + ts_xxxzz_xxyyy[i] * pa_z[i];

        ts_xxxzzz_xxyyz[i] = 2.0 * ts_xzzz_xxyyz[i] * fe_0 + 2.0 * ts_xxzzz_xyyz[i] * fe_0 + ts_xxzzz_xxyyz[i] * pa_x[i];

        ts_xxxzzz_xxyzz[i] = 2.0 * ts_xzzz_xxyzz[i] * fe_0 + 2.0 * ts_xxzzz_xyzz[i] * fe_0 + ts_xxzzz_xxyzz[i] * pa_x[i];

        ts_xxxzzz_xxzzz[i] = 2.0 * ts_xzzz_xxzzz[i] * fe_0 + 2.0 * ts_xxzzz_xzzz[i] * fe_0 + ts_xxzzz_xxzzz[i] * pa_x[i];

        ts_xxxzzz_xyyyy[i] = 2.0 * ts_xxxz_xyyyy[i] * fe_0 + ts_xxxzz_xyyyy[i] * pa_z[i];

        ts_xxxzzz_xyyyz[i] = 2.0 * ts_xzzz_xyyyz[i] * fe_0 + ts_xxzzz_yyyz[i] * fe_0 + ts_xxzzz_xyyyz[i] * pa_x[i];

        ts_xxxzzz_xyyzz[i] = 2.0 * ts_xzzz_xyyzz[i] * fe_0 + ts_xxzzz_yyzz[i] * fe_0 + ts_xxzzz_xyyzz[i] * pa_x[i];

        ts_xxxzzz_xyzzz[i] = 2.0 * ts_xzzz_xyzzz[i] * fe_0 + ts_xxzzz_yzzz[i] * fe_0 + ts_xxzzz_xyzzz[i] * pa_x[i];

        ts_xxxzzz_xzzzz[i] = 2.0 * ts_xzzz_xzzzz[i] * fe_0 + ts_xxzzz_zzzz[i] * fe_0 + ts_xxzzz_xzzzz[i] * pa_x[i];

        ts_xxxzzz_yyyyy[i] = 2.0 * ts_xzzz_yyyyy[i] * fe_0 + ts_xxzzz_yyyyy[i] * pa_x[i];

        ts_xxxzzz_yyyyz[i] = 2.0 * ts_xzzz_yyyyz[i] * fe_0 + ts_xxzzz_yyyyz[i] * pa_x[i];

        ts_xxxzzz_yyyzz[i] = 2.0 * ts_xzzz_yyyzz[i] * fe_0 + ts_xxzzz_yyyzz[i] * pa_x[i];

        ts_xxxzzz_yyzzz[i] = 2.0 * ts_xzzz_yyzzz[i] * fe_0 + ts_xxzzz_yyzzz[i] * pa_x[i];

        ts_xxxzzz_yzzzz[i] = 2.0 * ts_xzzz_yzzzz[i] * fe_0 + ts_xxzzz_yzzzz[i] * pa_x[i];

        ts_xxxzzz_zzzzz[i] = 2.0 * ts_xzzz_zzzzz[i] * fe_0 + ts_xxzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 210-231 components of targeted buffer : IH

    auto ts_xxyyyy_xxxxx = pbuffer.data(idx_ovl_ih + 210);

    auto ts_xxyyyy_xxxxy = pbuffer.data(idx_ovl_ih + 211);

    auto ts_xxyyyy_xxxxz = pbuffer.data(idx_ovl_ih + 212);

    auto ts_xxyyyy_xxxyy = pbuffer.data(idx_ovl_ih + 213);

    auto ts_xxyyyy_xxxyz = pbuffer.data(idx_ovl_ih + 214);

    auto ts_xxyyyy_xxxzz = pbuffer.data(idx_ovl_ih + 215);

    auto ts_xxyyyy_xxyyy = pbuffer.data(idx_ovl_ih + 216);

    auto ts_xxyyyy_xxyyz = pbuffer.data(idx_ovl_ih + 217);

    auto ts_xxyyyy_xxyzz = pbuffer.data(idx_ovl_ih + 218);

    auto ts_xxyyyy_xxzzz = pbuffer.data(idx_ovl_ih + 219);

    auto ts_xxyyyy_xyyyy = pbuffer.data(idx_ovl_ih + 220);

    auto ts_xxyyyy_xyyyz = pbuffer.data(idx_ovl_ih + 221);

    auto ts_xxyyyy_xyyzz = pbuffer.data(idx_ovl_ih + 222);

    auto ts_xxyyyy_xyzzz = pbuffer.data(idx_ovl_ih + 223);

    auto ts_xxyyyy_xzzzz = pbuffer.data(idx_ovl_ih + 224);

    auto ts_xxyyyy_yyyyy = pbuffer.data(idx_ovl_ih + 225);

    auto ts_xxyyyy_yyyyz = pbuffer.data(idx_ovl_ih + 226);

    auto ts_xxyyyy_yyyzz = pbuffer.data(idx_ovl_ih + 227);

    auto ts_xxyyyy_yyzzz = pbuffer.data(idx_ovl_ih + 228);

    auto ts_xxyyyy_yzzzz = pbuffer.data(idx_ovl_ih + 229);

    auto ts_xxyyyy_zzzzz = pbuffer.data(idx_ovl_ih + 230);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyy_xxxxx, ts_xxyy_xxxxz, ts_xxyy_xxxzz, ts_xxyy_xxzzz, ts_xxyy_xzzzz, ts_xxyyy_xxxxx, ts_xxyyy_xxxxz, ts_xxyyy_xxxzz, ts_xxyyy_xxzzz, ts_xxyyy_xzzzz, ts_xxyyyy_xxxxx, ts_xxyyyy_xxxxy, ts_xxyyyy_xxxxz, ts_xxyyyy_xxxyy, ts_xxyyyy_xxxyz, ts_xxyyyy_xxxzz, ts_xxyyyy_xxyyy, ts_xxyyyy_xxyyz, ts_xxyyyy_xxyzz, ts_xxyyyy_xxzzz, ts_xxyyyy_xyyyy, ts_xxyyyy_xyyyz, ts_xxyyyy_xyyzz, ts_xxyyyy_xyzzz, ts_xxyyyy_xzzzz, ts_xxyyyy_yyyyy, ts_xxyyyy_yyyyz, ts_xxyyyy_yyyzz, ts_xxyyyy_yyzzz, ts_xxyyyy_yzzzz, ts_xxyyyy_zzzzz, ts_xyyyy_xxxxy, ts_xyyyy_xxxy, ts_xyyyy_xxxyy, ts_xyyyy_xxxyz, ts_xyyyy_xxyy, ts_xyyyy_xxyyy, ts_xyyyy_xxyyz, ts_xyyyy_xxyz, ts_xyyyy_xxyzz, ts_xyyyy_xyyy, ts_xyyyy_xyyyy, ts_xyyyy_xyyyz, ts_xyyyy_xyyz, ts_xyyyy_xyyzz, ts_xyyyy_xyzz, ts_xyyyy_xyzzz, ts_xyyyy_yyyy, ts_xyyyy_yyyyy, ts_xyyyy_yyyyz, ts_xyyyy_yyyz, ts_xyyyy_yyyzz, ts_xyyyy_yyzz, ts_xyyyy_yyzzz, ts_xyyyy_yzzz, ts_xyyyy_yzzzz, ts_xyyyy_zzzzz, ts_yyyy_xxxxy, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyzz, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyzz, ts_yyyy_xyzzz, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyzz, ts_yyyy_yyzzz, ts_yyyy_yzzzz, ts_yyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyyy_xxxxx[i] = 3.0 * ts_xxyy_xxxxx[i] * fe_0 + ts_xxyyy_xxxxx[i] * pa_y[i];

        ts_xxyyyy_xxxxy[i] = ts_yyyy_xxxxy[i] * fe_0 + 4.0 * ts_xyyyy_xxxy[i] * fe_0 + ts_xyyyy_xxxxy[i] * pa_x[i];

        ts_xxyyyy_xxxxz[i] = 3.0 * ts_xxyy_xxxxz[i] * fe_0 + ts_xxyyy_xxxxz[i] * pa_y[i];

        ts_xxyyyy_xxxyy[i] = ts_yyyy_xxxyy[i] * fe_0 + 3.0 * ts_xyyyy_xxyy[i] * fe_0 + ts_xyyyy_xxxyy[i] * pa_x[i];

        ts_xxyyyy_xxxyz[i] = ts_yyyy_xxxyz[i] * fe_0 + 3.0 * ts_xyyyy_xxyz[i] * fe_0 + ts_xyyyy_xxxyz[i] * pa_x[i];

        ts_xxyyyy_xxxzz[i] = 3.0 * ts_xxyy_xxxzz[i] * fe_0 + ts_xxyyy_xxxzz[i] * pa_y[i];

        ts_xxyyyy_xxyyy[i] = ts_yyyy_xxyyy[i] * fe_0 + 2.0 * ts_xyyyy_xyyy[i] * fe_0 + ts_xyyyy_xxyyy[i] * pa_x[i];

        ts_xxyyyy_xxyyz[i] = ts_yyyy_xxyyz[i] * fe_0 + 2.0 * ts_xyyyy_xyyz[i] * fe_0 + ts_xyyyy_xxyyz[i] * pa_x[i];

        ts_xxyyyy_xxyzz[i] = ts_yyyy_xxyzz[i] * fe_0 + 2.0 * ts_xyyyy_xyzz[i] * fe_0 + ts_xyyyy_xxyzz[i] * pa_x[i];

        ts_xxyyyy_xxzzz[i] = 3.0 * ts_xxyy_xxzzz[i] * fe_0 + ts_xxyyy_xxzzz[i] * pa_y[i];

        ts_xxyyyy_xyyyy[i] = ts_yyyy_xyyyy[i] * fe_0 + ts_xyyyy_yyyy[i] * fe_0 + ts_xyyyy_xyyyy[i] * pa_x[i];

        ts_xxyyyy_xyyyz[i] = ts_yyyy_xyyyz[i] * fe_0 + ts_xyyyy_yyyz[i] * fe_0 + ts_xyyyy_xyyyz[i] * pa_x[i];

        ts_xxyyyy_xyyzz[i] = ts_yyyy_xyyzz[i] * fe_0 + ts_xyyyy_yyzz[i] * fe_0 + ts_xyyyy_xyyzz[i] * pa_x[i];

        ts_xxyyyy_xyzzz[i] = ts_yyyy_xyzzz[i] * fe_0 + ts_xyyyy_yzzz[i] * fe_0 + ts_xyyyy_xyzzz[i] * pa_x[i];

        ts_xxyyyy_xzzzz[i] = 3.0 * ts_xxyy_xzzzz[i] * fe_0 + ts_xxyyy_xzzzz[i] * pa_y[i];

        ts_xxyyyy_yyyyy[i] = ts_yyyy_yyyyy[i] * fe_0 + ts_xyyyy_yyyyy[i] * pa_x[i];

        ts_xxyyyy_yyyyz[i] = ts_yyyy_yyyyz[i] * fe_0 + ts_xyyyy_yyyyz[i] * pa_x[i];

        ts_xxyyyy_yyyzz[i] = ts_yyyy_yyyzz[i] * fe_0 + ts_xyyyy_yyyzz[i] * pa_x[i];

        ts_xxyyyy_yyzzz[i] = ts_yyyy_yyzzz[i] * fe_0 + ts_xyyyy_yyzzz[i] * pa_x[i];

        ts_xxyyyy_yzzzz[i] = ts_yyyy_yzzzz[i] * fe_0 + ts_xyyyy_yzzzz[i] * pa_x[i];

        ts_xxyyyy_zzzzz[i] = ts_yyyy_zzzzz[i] * fe_0 + ts_xyyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 231-252 components of targeted buffer : IH

    auto ts_xxyyyz_xxxxx = pbuffer.data(idx_ovl_ih + 231);

    auto ts_xxyyyz_xxxxy = pbuffer.data(idx_ovl_ih + 232);

    auto ts_xxyyyz_xxxxz = pbuffer.data(idx_ovl_ih + 233);

    auto ts_xxyyyz_xxxyy = pbuffer.data(idx_ovl_ih + 234);

    auto ts_xxyyyz_xxxyz = pbuffer.data(idx_ovl_ih + 235);

    auto ts_xxyyyz_xxxzz = pbuffer.data(idx_ovl_ih + 236);

    auto ts_xxyyyz_xxyyy = pbuffer.data(idx_ovl_ih + 237);

    auto ts_xxyyyz_xxyyz = pbuffer.data(idx_ovl_ih + 238);

    auto ts_xxyyyz_xxyzz = pbuffer.data(idx_ovl_ih + 239);

    auto ts_xxyyyz_xxzzz = pbuffer.data(idx_ovl_ih + 240);

    auto ts_xxyyyz_xyyyy = pbuffer.data(idx_ovl_ih + 241);

    auto ts_xxyyyz_xyyyz = pbuffer.data(idx_ovl_ih + 242);

    auto ts_xxyyyz_xyyzz = pbuffer.data(idx_ovl_ih + 243);

    auto ts_xxyyyz_xyzzz = pbuffer.data(idx_ovl_ih + 244);

    auto ts_xxyyyz_xzzzz = pbuffer.data(idx_ovl_ih + 245);

    auto ts_xxyyyz_yyyyy = pbuffer.data(idx_ovl_ih + 246);

    auto ts_xxyyyz_yyyyz = pbuffer.data(idx_ovl_ih + 247);

    auto ts_xxyyyz_yyyzz = pbuffer.data(idx_ovl_ih + 248);

    auto ts_xxyyyz_yyzzz = pbuffer.data(idx_ovl_ih + 249);

    auto ts_xxyyyz_yzzzz = pbuffer.data(idx_ovl_ih + 250);

    auto ts_xxyyyz_zzzzz = pbuffer.data(idx_ovl_ih + 251);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyyy_xxxxx, ts_xxyyy_xxxxy, ts_xxyyy_xxxy, ts_xxyyy_xxxyy, ts_xxyyy_xxxyz, ts_xxyyy_xxyy, ts_xxyyy_xxyyy, ts_xxyyy_xxyyz, ts_xxyyy_xxyz, ts_xxyyy_xxyzz, ts_xxyyy_xyyy, ts_xxyyy_xyyyy, ts_xxyyy_xyyyz, ts_xxyyy_xyyz, ts_xxyyy_xyyzz, ts_xxyyy_xyzz, ts_xxyyy_xyzzz, ts_xxyyy_yyyyy, ts_xxyyyz_xxxxx, ts_xxyyyz_xxxxy, ts_xxyyyz_xxxxz, ts_xxyyyz_xxxyy, ts_xxyyyz_xxxyz, ts_xxyyyz_xxxzz, ts_xxyyyz_xxyyy, ts_xxyyyz_xxyyz, ts_xxyyyz_xxyzz, ts_xxyyyz_xxzzz, ts_xxyyyz_xyyyy, ts_xxyyyz_xyyyz, ts_xxyyyz_xyyzz, ts_xxyyyz_xyzzz, ts_xxyyyz_xzzzz, ts_xxyyyz_yyyyy, ts_xxyyyz_yyyyz, ts_xxyyyz_yyyzz, ts_xxyyyz_yyzzz, ts_xxyyyz_yzzzz, ts_xxyyyz_zzzzz, ts_xxyyz_xxxxz, ts_xxyyz_xxxzz, ts_xxyyz_xxzzz, ts_xxyyz_xzzzz, ts_xxyz_xxxxz, ts_xxyz_xxxzz, ts_xxyz_xxzzz, ts_xxyz_xzzzz, ts_xyyyz_yyyyz, ts_xyyyz_yyyzz, ts_xyyyz_yyzzz, ts_xyyyz_yzzzz, ts_xyyyz_zzzzz, ts_yyyz_yyyyz, ts_yyyz_yyyzz, ts_yyyz_yyzzz, ts_yyyz_yzzzz, ts_yyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyyz_xxxxx[i] = ts_xxyyy_xxxxx[i] * pa_z[i];

        ts_xxyyyz_xxxxy[i] = ts_xxyyy_xxxxy[i] * pa_z[i];

        ts_xxyyyz_xxxxz[i] = 2.0 * ts_xxyz_xxxxz[i] * fe_0 + ts_xxyyz_xxxxz[i] * pa_y[i];

        ts_xxyyyz_xxxyy[i] = ts_xxyyy_xxxyy[i] * pa_z[i];

        ts_xxyyyz_xxxyz[i] = ts_xxyyy_xxxy[i] * fe_0 + ts_xxyyy_xxxyz[i] * pa_z[i];

        ts_xxyyyz_xxxzz[i] = 2.0 * ts_xxyz_xxxzz[i] * fe_0 + ts_xxyyz_xxxzz[i] * pa_y[i];

        ts_xxyyyz_xxyyy[i] = ts_xxyyy_xxyyy[i] * pa_z[i];

        ts_xxyyyz_xxyyz[i] = ts_xxyyy_xxyy[i] * fe_0 + ts_xxyyy_xxyyz[i] * pa_z[i];

        ts_xxyyyz_xxyzz[i] = 2.0 * ts_xxyyy_xxyz[i] * fe_0 + ts_xxyyy_xxyzz[i] * pa_z[i];

        ts_xxyyyz_xxzzz[i] = 2.0 * ts_xxyz_xxzzz[i] * fe_0 + ts_xxyyz_xxzzz[i] * pa_y[i];

        ts_xxyyyz_xyyyy[i] = ts_xxyyy_xyyyy[i] * pa_z[i];

        ts_xxyyyz_xyyyz[i] = ts_xxyyy_xyyy[i] * fe_0 + ts_xxyyy_xyyyz[i] * pa_z[i];

        ts_xxyyyz_xyyzz[i] = 2.0 * ts_xxyyy_xyyz[i] * fe_0 + ts_xxyyy_xyyzz[i] * pa_z[i];

        ts_xxyyyz_xyzzz[i] = 3.0 * ts_xxyyy_xyzz[i] * fe_0 + ts_xxyyy_xyzzz[i] * pa_z[i];

        ts_xxyyyz_xzzzz[i] = 2.0 * ts_xxyz_xzzzz[i] * fe_0 + ts_xxyyz_xzzzz[i] * pa_y[i];

        ts_xxyyyz_yyyyy[i] = ts_xxyyy_yyyyy[i] * pa_z[i];

        ts_xxyyyz_yyyyz[i] = ts_yyyz_yyyyz[i] * fe_0 + ts_xyyyz_yyyyz[i] * pa_x[i];

        ts_xxyyyz_yyyzz[i] = ts_yyyz_yyyzz[i] * fe_0 + ts_xyyyz_yyyzz[i] * pa_x[i];

        ts_xxyyyz_yyzzz[i] = ts_yyyz_yyzzz[i] * fe_0 + ts_xyyyz_yyzzz[i] * pa_x[i];

        ts_xxyyyz_yzzzz[i] = ts_yyyz_yzzzz[i] * fe_0 + ts_xyyyz_yzzzz[i] * pa_x[i];

        ts_xxyyyz_zzzzz[i] = ts_yyyz_zzzzz[i] * fe_0 + ts_xyyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 252-273 components of targeted buffer : IH

    auto ts_xxyyzz_xxxxx = pbuffer.data(idx_ovl_ih + 252);

    auto ts_xxyyzz_xxxxy = pbuffer.data(idx_ovl_ih + 253);

    auto ts_xxyyzz_xxxxz = pbuffer.data(idx_ovl_ih + 254);

    auto ts_xxyyzz_xxxyy = pbuffer.data(idx_ovl_ih + 255);

    auto ts_xxyyzz_xxxyz = pbuffer.data(idx_ovl_ih + 256);

    auto ts_xxyyzz_xxxzz = pbuffer.data(idx_ovl_ih + 257);

    auto ts_xxyyzz_xxyyy = pbuffer.data(idx_ovl_ih + 258);

    auto ts_xxyyzz_xxyyz = pbuffer.data(idx_ovl_ih + 259);

    auto ts_xxyyzz_xxyzz = pbuffer.data(idx_ovl_ih + 260);

    auto ts_xxyyzz_xxzzz = pbuffer.data(idx_ovl_ih + 261);

    auto ts_xxyyzz_xyyyy = pbuffer.data(idx_ovl_ih + 262);

    auto ts_xxyyzz_xyyyz = pbuffer.data(idx_ovl_ih + 263);

    auto ts_xxyyzz_xyyzz = pbuffer.data(idx_ovl_ih + 264);

    auto ts_xxyyzz_xyzzz = pbuffer.data(idx_ovl_ih + 265);

    auto ts_xxyyzz_xzzzz = pbuffer.data(idx_ovl_ih + 266);

    auto ts_xxyyzz_yyyyy = pbuffer.data(idx_ovl_ih + 267);

    auto ts_xxyyzz_yyyyz = pbuffer.data(idx_ovl_ih + 268);

    auto ts_xxyyzz_yyyzz = pbuffer.data(idx_ovl_ih + 269);

    auto ts_xxyyzz_yyzzz = pbuffer.data(idx_ovl_ih + 270);

    auto ts_xxyyzz_yzzzz = pbuffer.data(idx_ovl_ih + 271);

    auto ts_xxyyzz_zzzzz = pbuffer.data(idx_ovl_ih + 272);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyy_xxxxy, ts_xxyy_xxxyy, ts_xxyy_xxyyy, ts_xxyy_xyyyy, ts_xxyyz_xxxxy, ts_xxyyz_xxxyy, ts_xxyyz_xxyyy, ts_xxyyz_xyyyy, ts_xxyyzz_xxxxx, ts_xxyyzz_xxxxy, ts_xxyyzz_xxxxz, ts_xxyyzz_xxxyy, ts_xxyyzz_xxxyz, ts_xxyyzz_xxxzz, ts_xxyyzz_xxyyy, ts_xxyyzz_xxyyz, ts_xxyyzz_xxyzz, ts_xxyyzz_xxzzz, ts_xxyyzz_xyyyy, ts_xxyyzz_xyyyz, ts_xxyyzz_xyyzz, ts_xxyyzz_xyzzz, ts_xxyyzz_xzzzz, ts_xxyyzz_yyyyy, ts_xxyyzz_yyyyz, ts_xxyyzz_yyyzz, ts_xxyyzz_yyzzz, ts_xxyyzz_yzzzz, ts_xxyyzz_zzzzz, ts_xxyzz_xxxxx, ts_xxyzz_xxxxz, ts_xxyzz_xxxzz, ts_xxyzz_xxzzz, ts_xxyzz_xzzzz, ts_xxzz_xxxxx, ts_xxzz_xxxxz, ts_xxzz_xxxzz, ts_xxzz_xxzzz, ts_xxzz_xzzzz, ts_xyyzz_xxxyz, ts_xyyzz_xxyyz, ts_xyyzz_xxyz, ts_xyyzz_xxyzz, ts_xyyzz_xyyyz, ts_xyyzz_xyyz, ts_xyyzz_xyyzz, ts_xyyzz_xyzz, ts_xyyzz_xyzzz, ts_xyyzz_yyyyy, ts_xyyzz_yyyyz, ts_xyyzz_yyyz, ts_xyyzz_yyyzz, ts_xyyzz_yyzz, ts_xyyzz_yyzzz, ts_xyyzz_yzzz, ts_xyyzz_yzzzz, ts_xyyzz_zzzzz, ts_yyzz_xxxyz, ts_yyzz_xxyyz, ts_yyzz_xxyzz, ts_yyzz_xyyyz, ts_yyzz_xyyzz, ts_yyzz_xyzzz, ts_yyzz_yyyyy, ts_yyzz_yyyyz, ts_yyzz_yyyzz, ts_yyzz_yyzzz, ts_yyzz_yzzzz, ts_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyzz_xxxxx[i] = ts_xxzz_xxxxx[i] * fe_0 + ts_xxyzz_xxxxx[i] * pa_y[i];

        ts_xxyyzz_xxxxy[i] = ts_xxyy_xxxxy[i] * fe_0 + ts_xxyyz_xxxxy[i] * pa_z[i];

        ts_xxyyzz_xxxxz[i] = ts_xxzz_xxxxz[i] * fe_0 + ts_xxyzz_xxxxz[i] * pa_y[i];

        ts_xxyyzz_xxxyy[i] = ts_xxyy_xxxyy[i] * fe_0 + ts_xxyyz_xxxyy[i] * pa_z[i];

        ts_xxyyzz_xxxyz[i] = ts_yyzz_xxxyz[i] * fe_0 + 3.0 * ts_xyyzz_xxyz[i] * fe_0 + ts_xyyzz_xxxyz[i] * pa_x[i];

        ts_xxyyzz_xxxzz[i] = ts_xxzz_xxxzz[i] * fe_0 + ts_xxyzz_xxxzz[i] * pa_y[i];

        ts_xxyyzz_xxyyy[i] = ts_xxyy_xxyyy[i] * fe_0 + ts_xxyyz_xxyyy[i] * pa_z[i];

        ts_xxyyzz_xxyyz[i] = ts_yyzz_xxyyz[i] * fe_0 + 2.0 * ts_xyyzz_xyyz[i] * fe_0 + ts_xyyzz_xxyyz[i] * pa_x[i];

        ts_xxyyzz_xxyzz[i] = ts_yyzz_xxyzz[i] * fe_0 + 2.0 * ts_xyyzz_xyzz[i] * fe_0 + ts_xyyzz_xxyzz[i] * pa_x[i];

        ts_xxyyzz_xxzzz[i] = ts_xxzz_xxzzz[i] * fe_0 + ts_xxyzz_xxzzz[i] * pa_y[i];

        ts_xxyyzz_xyyyy[i] = ts_xxyy_xyyyy[i] * fe_0 + ts_xxyyz_xyyyy[i] * pa_z[i];

        ts_xxyyzz_xyyyz[i] = ts_yyzz_xyyyz[i] * fe_0 + ts_xyyzz_yyyz[i] * fe_0 + ts_xyyzz_xyyyz[i] * pa_x[i];

        ts_xxyyzz_xyyzz[i] = ts_yyzz_xyyzz[i] * fe_0 + ts_xyyzz_yyzz[i] * fe_0 + ts_xyyzz_xyyzz[i] * pa_x[i];

        ts_xxyyzz_xyzzz[i] = ts_yyzz_xyzzz[i] * fe_0 + ts_xyyzz_yzzz[i] * fe_0 + ts_xyyzz_xyzzz[i] * pa_x[i];

        ts_xxyyzz_xzzzz[i] = ts_xxzz_xzzzz[i] * fe_0 + ts_xxyzz_xzzzz[i] * pa_y[i];

        ts_xxyyzz_yyyyy[i] = ts_yyzz_yyyyy[i] * fe_0 + ts_xyyzz_yyyyy[i] * pa_x[i];

        ts_xxyyzz_yyyyz[i] = ts_yyzz_yyyyz[i] * fe_0 + ts_xyyzz_yyyyz[i] * pa_x[i];

        ts_xxyyzz_yyyzz[i] = ts_yyzz_yyyzz[i] * fe_0 + ts_xyyzz_yyyzz[i] * pa_x[i];

        ts_xxyyzz_yyzzz[i] = ts_yyzz_yyzzz[i] * fe_0 + ts_xyyzz_yyzzz[i] * pa_x[i];

        ts_xxyyzz_yzzzz[i] = ts_yyzz_yzzzz[i] * fe_0 + ts_xyyzz_yzzzz[i] * pa_x[i];

        ts_xxyyzz_zzzzz[i] = ts_yyzz_zzzzz[i] * fe_0 + ts_xyyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 273-294 components of targeted buffer : IH

    auto ts_xxyzzz_xxxxx = pbuffer.data(idx_ovl_ih + 273);

    auto ts_xxyzzz_xxxxy = pbuffer.data(idx_ovl_ih + 274);

    auto ts_xxyzzz_xxxxz = pbuffer.data(idx_ovl_ih + 275);

    auto ts_xxyzzz_xxxyy = pbuffer.data(idx_ovl_ih + 276);

    auto ts_xxyzzz_xxxyz = pbuffer.data(idx_ovl_ih + 277);

    auto ts_xxyzzz_xxxzz = pbuffer.data(idx_ovl_ih + 278);

    auto ts_xxyzzz_xxyyy = pbuffer.data(idx_ovl_ih + 279);

    auto ts_xxyzzz_xxyyz = pbuffer.data(idx_ovl_ih + 280);

    auto ts_xxyzzz_xxyzz = pbuffer.data(idx_ovl_ih + 281);

    auto ts_xxyzzz_xxzzz = pbuffer.data(idx_ovl_ih + 282);

    auto ts_xxyzzz_xyyyy = pbuffer.data(idx_ovl_ih + 283);

    auto ts_xxyzzz_xyyyz = pbuffer.data(idx_ovl_ih + 284);

    auto ts_xxyzzz_xyyzz = pbuffer.data(idx_ovl_ih + 285);

    auto ts_xxyzzz_xyzzz = pbuffer.data(idx_ovl_ih + 286);

    auto ts_xxyzzz_xzzzz = pbuffer.data(idx_ovl_ih + 287);

    auto ts_xxyzzz_yyyyy = pbuffer.data(idx_ovl_ih + 288);

    auto ts_xxyzzz_yyyyz = pbuffer.data(idx_ovl_ih + 289);

    auto ts_xxyzzz_yyyzz = pbuffer.data(idx_ovl_ih + 290);

    auto ts_xxyzzz_yyzzz = pbuffer.data(idx_ovl_ih + 291);

    auto ts_xxyzzz_yzzzz = pbuffer.data(idx_ovl_ih + 292);

    auto ts_xxyzzz_zzzzz = pbuffer.data(idx_ovl_ih + 293);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyzzz_xxxxx, ts_xxyzzz_xxxxy, ts_xxyzzz_xxxxz, ts_xxyzzz_xxxyy, ts_xxyzzz_xxxyz, ts_xxyzzz_xxxzz, ts_xxyzzz_xxyyy, ts_xxyzzz_xxyyz, ts_xxyzzz_xxyzz, ts_xxyzzz_xxzzz, ts_xxyzzz_xyyyy, ts_xxyzzz_xyyyz, ts_xxyzzz_xyyzz, ts_xxyzzz_xyzzz, ts_xxyzzz_xzzzz, ts_xxyzzz_yyyyy, ts_xxyzzz_yyyyz, ts_xxyzzz_yyyzz, ts_xxyzzz_yyzzz, ts_xxyzzz_yzzzz, ts_xxyzzz_zzzzz, ts_xxzzz_xxxx, ts_xxzzz_xxxxx, ts_xxzzz_xxxxy, ts_xxzzz_xxxxz, ts_xxzzz_xxxy, ts_xxzzz_xxxyy, ts_xxzzz_xxxyz, ts_xxzzz_xxxz, ts_xxzzz_xxxzz, ts_xxzzz_xxyy, ts_xxzzz_xxyyy, ts_xxzzz_xxyyz, ts_xxzzz_xxyz, ts_xxzzz_xxyzz, ts_xxzzz_xxzz, ts_xxzzz_xxzzz, ts_xxzzz_xyyy, ts_xxzzz_xyyyy, ts_xxzzz_xyyyz, ts_xxzzz_xyyz, ts_xxzzz_xyyzz, ts_xxzzz_xyzz, ts_xxzzz_xyzzz, ts_xxzzz_xzzz, ts_xxzzz_xzzzz, ts_xxzzz_zzzzz, ts_xyzzz_yyyyy, ts_xyzzz_yyyyz, ts_xyzzz_yyyzz, ts_xyzzz_yyzzz, ts_xyzzz_yzzzz, ts_yzzz_yyyyy, ts_yzzz_yyyyz, ts_yzzz_yyyzz, ts_yzzz_yyzzz, ts_yzzz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyzzz_xxxxx[i] = ts_xxzzz_xxxxx[i] * pa_y[i];

        ts_xxyzzz_xxxxy[i] = ts_xxzzz_xxxx[i] * fe_0 + ts_xxzzz_xxxxy[i] * pa_y[i];

        ts_xxyzzz_xxxxz[i] = ts_xxzzz_xxxxz[i] * pa_y[i];

        ts_xxyzzz_xxxyy[i] = 2.0 * ts_xxzzz_xxxy[i] * fe_0 + ts_xxzzz_xxxyy[i] * pa_y[i];

        ts_xxyzzz_xxxyz[i] = ts_xxzzz_xxxz[i] * fe_0 + ts_xxzzz_xxxyz[i] * pa_y[i];

        ts_xxyzzz_xxxzz[i] = ts_xxzzz_xxxzz[i] * pa_y[i];

        ts_xxyzzz_xxyyy[i] = 3.0 * ts_xxzzz_xxyy[i] * fe_0 + ts_xxzzz_xxyyy[i] * pa_y[i];

        ts_xxyzzz_xxyyz[i] = 2.0 * ts_xxzzz_xxyz[i] * fe_0 + ts_xxzzz_xxyyz[i] * pa_y[i];

        ts_xxyzzz_xxyzz[i] = ts_xxzzz_xxzz[i] * fe_0 + ts_xxzzz_xxyzz[i] * pa_y[i];

        ts_xxyzzz_xxzzz[i] = ts_xxzzz_xxzzz[i] * pa_y[i];

        ts_xxyzzz_xyyyy[i] = 4.0 * ts_xxzzz_xyyy[i] * fe_0 + ts_xxzzz_xyyyy[i] * pa_y[i];

        ts_xxyzzz_xyyyz[i] = 3.0 * ts_xxzzz_xyyz[i] * fe_0 + ts_xxzzz_xyyyz[i] * pa_y[i];

        ts_xxyzzz_xyyzz[i] = 2.0 * ts_xxzzz_xyzz[i] * fe_0 + ts_xxzzz_xyyzz[i] * pa_y[i];

        ts_xxyzzz_xyzzz[i] = ts_xxzzz_xzzz[i] * fe_0 + ts_xxzzz_xyzzz[i] * pa_y[i];

        ts_xxyzzz_xzzzz[i] = ts_xxzzz_xzzzz[i] * pa_y[i];

        ts_xxyzzz_yyyyy[i] = ts_yzzz_yyyyy[i] * fe_0 + ts_xyzzz_yyyyy[i] * pa_x[i];

        ts_xxyzzz_yyyyz[i] = ts_yzzz_yyyyz[i] * fe_0 + ts_xyzzz_yyyyz[i] * pa_x[i];

        ts_xxyzzz_yyyzz[i] = ts_yzzz_yyyzz[i] * fe_0 + ts_xyzzz_yyyzz[i] * pa_x[i];

        ts_xxyzzz_yyzzz[i] = ts_yzzz_yyzzz[i] * fe_0 + ts_xyzzz_yyzzz[i] * pa_x[i];

        ts_xxyzzz_yzzzz[i] = ts_yzzz_yzzzz[i] * fe_0 + ts_xyzzz_yzzzz[i] * pa_x[i];

        ts_xxyzzz_zzzzz[i] = ts_xxzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 294-315 components of targeted buffer : IH

    auto ts_xxzzzz_xxxxx = pbuffer.data(idx_ovl_ih + 294);

    auto ts_xxzzzz_xxxxy = pbuffer.data(idx_ovl_ih + 295);

    auto ts_xxzzzz_xxxxz = pbuffer.data(idx_ovl_ih + 296);

    auto ts_xxzzzz_xxxyy = pbuffer.data(idx_ovl_ih + 297);

    auto ts_xxzzzz_xxxyz = pbuffer.data(idx_ovl_ih + 298);

    auto ts_xxzzzz_xxxzz = pbuffer.data(idx_ovl_ih + 299);

    auto ts_xxzzzz_xxyyy = pbuffer.data(idx_ovl_ih + 300);

    auto ts_xxzzzz_xxyyz = pbuffer.data(idx_ovl_ih + 301);

    auto ts_xxzzzz_xxyzz = pbuffer.data(idx_ovl_ih + 302);

    auto ts_xxzzzz_xxzzz = pbuffer.data(idx_ovl_ih + 303);

    auto ts_xxzzzz_xyyyy = pbuffer.data(idx_ovl_ih + 304);

    auto ts_xxzzzz_xyyyz = pbuffer.data(idx_ovl_ih + 305);

    auto ts_xxzzzz_xyyzz = pbuffer.data(idx_ovl_ih + 306);

    auto ts_xxzzzz_xyzzz = pbuffer.data(idx_ovl_ih + 307);

    auto ts_xxzzzz_xzzzz = pbuffer.data(idx_ovl_ih + 308);

    auto ts_xxzzzz_yyyyy = pbuffer.data(idx_ovl_ih + 309);

    auto ts_xxzzzz_yyyyz = pbuffer.data(idx_ovl_ih + 310);

    auto ts_xxzzzz_yyyzz = pbuffer.data(idx_ovl_ih + 311);

    auto ts_xxzzzz_yyzzz = pbuffer.data(idx_ovl_ih + 312);

    auto ts_xxzzzz_yzzzz = pbuffer.data(idx_ovl_ih + 313);

    auto ts_xxzzzz_zzzzz = pbuffer.data(idx_ovl_ih + 314);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxzz_xxxxx, ts_xxzz_xxxxy, ts_xxzz_xxxyy, ts_xxzz_xxyyy, ts_xxzz_xyyyy, ts_xxzzz_xxxxx, ts_xxzzz_xxxxy, ts_xxzzz_xxxyy, ts_xxzzz_xxyyy, ts_xxzzz_xyyyy, ts_xxzzzz_xxxxx, ts_xxzzzz_xxxxy, ts_xxzzzz_xxxxz, ts_xxzzzz_xxxyy, ts_xxzzzz_xxxyz, ts_xxzzzz_xxxzz, ts_xxzzzz_xxyyy, ts_xxzzzz_xxyyz, ts_xxzzzz_xxyzz, ts_xxzzzz_xxzzz, ts_xxzzzz_xyyyy, ts_xxzzzz_xyyyz, ts_xxzzzz_xyyzz, ts_xxzzzz_xyzzz, ts_xxzzzz_xzzzz, ts_xxzzzz_yyyyy, ts_xxzzzz_yyyyz, ts_xxzzzz_yyyzz, ts_xxzzzz_yyzzz, ts_xxzzzz_yzzzz, ts_xxzzzz_zzzzz, ts_xzzzz_xxxxz, ts_xzzzz_xxxyz, ts_xzzzz_xxxz, ts_xzzzz_xxxzz, ts_xzzzz_xxyyz, ts_xzzzz_xxyz, ts_xzzzz_xxyzz, ts_xzzzz_xxzz, ts_xzzzz_xxzzz, ts_xzzzz_xyyyz, ts_xzzzz_xyyz, ts_xzzzz_xyyzz, ts_xzzzz_xyzz, ts_xzzzz_xyzzz, ts_xzzzz_xzzz, ts_xzzzz_xzzzz, ts_xzzzz_yyyyy, ts_xzzzz_yyyyz, ts_xzzzz_yyyz, ts_xzzzz_yyyzz, ts_xzzzz_yyzz, ts_xzzzz_yyzzz, ts_xzzzz_yzzz, ts_xzzzz_yzzzz, ts_xzzzz_zzzz, ts_xzzzz_zzzzz, ts_zzzz_xxxxz, ts_zzzz_xxxyz, ts_zzzz_xxxzz, ts_zzzz_xxyyz, ts_zzzz_xxyzz, ts_zzzz_xxzzz, ts_zzzz_xyyyz, ts_zzzz_xyyzz, ts_zzzz_xyzzz, ts_zzzz_xzzzz, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyzz, ts_zzzz_yyzzz, ts_zzzz_yzzzz, ts_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzzzz_xxxxx[i] = 3.0 * ts_xxzz_xxxxx[i] * fe_0 + ts_xxzzz_xxxxx[i] * pa_z[i];

        ts_xxzzzz_xxxxy[i] = 3.0 * ts_xxzz_xxxxy[i] * fe_0 + ts_xxzzz_xxxxy[i] * pa_z[i];

        ts_xxzzzz_xxxxz[i] = ts_zzzz_xxxxz[i] * fe_0 + 4.0 * ts_xzzzz_xxxz[i] * fe_0 + ts_xzzzz_xxxxz[i] * pa_x[i];

        ts_xxzzzz_xxxyy[i] = 3.0 * ts_xxzz_xxxyy[i] * fe_0 + ts_xxzzz_xxxyy[i] * pa_z[i];

        ts_xxzzzz_xxxyz[i] = ts_zzzz_xxxyz[i] * fe_0 + 3.0 * ts_xzzzz_xxyz[i] * fe_0 + ts_xzzzz_xxxyz[i] * pa_x[i];

        ts_xxzzzz_xxxzz[i] = ts_zzzz_xxxzz[i] * fe_0 + 3.0 * ts_xzzzz_xxzz[i] * fe_0 + ts_xzzzz_xxxzz[i] * pa_x[i];

        ts_xxzzzz_xxyyy[i] = 3.0 * ts_xxzz_xxyyy[i] * fe_0 + ts_xxzzz_xxyyy[i] * pa_z[i];

        ts_xxzzzz_xxyyz[i] = ts_zzzz_xxyyz[i] * fe_0 + 2.0 * ts_xzzzz_xyyz[i] * fe_0 + ts_xzzzz_xxyyz[i] * pa_x[i];

        ts_xxzzzz_xxyzz[i] = ts_zzzz_xxyzz[i] * fe_0 + 2.0 * ts_xzzzz_xyzz[i] * fe_0 + ts_xzzzz_xxyzz[i] * pa_x[i];

        ts_xxzzzz_xxzzz[i] = ts_zzzz_xxzzz[i] * fe_0 + 2.0 * ts_xzzzz_xzzz[i] * fe_0 + ts_xzzzz_xxzzz[i] * pa_x[i];

        ts_xxzzzz_xyyyy[i] = 3.0 * ts_xxzz_xyyyy[i] * fe_0 + ts_xxzzz_xyyyy[i] * pa_z[i];

        ts_xxzzzz_xyyyz[i] = ts_zzzz_xyyyz[i] * fe_0 + ts_xzzzz_yyyz[i] * fe_0 + ts_xzzzz_xyyyz[i] * pa_x[i];

        ts_xxzzzz_xyyzz[i] = ts_zzzz_xyyzz[i] * fe_0 + ts_xzzzz_yyzz[i] * fe_0 + ts_xzzzz_xyyzz[i] * pa_x[i];

        ts_xxzzzz_xyzzz[i] = ts_zzzz_xyzzz[i] * fe_0 + ts_xzzzz_yzzz[i] * fe_0 + ts_xzzzz_xyzzz[i] * pa_x[i];

        ts_xxzzzz_xzzzz[i] = ts_zzzz_xzzzz[i] * fe_0 + ts_xzzzz_zzzz[i] * fe_0 + ts_xzzzz_xzzzz[i] * pa_x[i];

        ts_xxzzzz_yyyyy[i] = ts_zzzz_yyyyy[i] * fe_0 + ts_xzzzz_yyyyy[i] * pa_x[i];

        ts_xxzzzz_yyyyz[i] = ts_zzzz_yyyyz[i] * fe_0 + ts_xzzzz_yyyyz[i] * pa_x[i];

        ts_xxzzzz_yyyzz[i] = ts_zzzz_yyyzz[i] * fe_0 + ts_xzzzz_yyyzz[i] * pa_x[i];

        ts_xxzzzz_yyzzz[i] = ts_zzzz_yyzzz[i] * fe_0 + ts_xzzzz_yyzzz[i] * pa_x[i];

        ts_xxzzzz_yzzzz[i] = ts_zzzz_yzzzz[i] * fe_0 + ts_xzzzz_yzzzz[i] * pa_x[i];

        ts_xxzzzz_zzzzz[i] = ts_zzzz_zzzzz[i] * fe_0 + ts_xzzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 315-336 components of targeted buffer : IH

    auto ts_xyyyyy_xxxxx = pbuffer.data(idx_ovl_ih + 315);

    auto ts_xyyyyy_xxxxy = pbuffer.data(idx_ovl_ih + 316);

    auto ts_xyyyyy_xxxxz = pbuffer.data(idx_ovl_ih + 317);

    auto ts_xyyyyy_xxxyy = pbuffer.data(idx_ovl_ih + 318);

    auto ts_xyyyyy_xxxyz = pbuffer.data(idx_ovl_ih + 319);

    auto ts_xyyyyy_xxxzz = pbuffer.data(idx_ovl_ih + 320);

    auto ts_xyyyyy_xxyyy = pbuffer.data(idx_ovl_ih + 321);

    auto ts_xyyyyy_xxyyz = pbuffer.data(idx_ovl_ih + 322);

    auto ts_xyyyyy_xxyzz = pbuffer.data(idx_ovl_ih + 323);

    auto ts_xyyyyy_xxzzz = pbuffer.data(idx_ovl_ih + 324);

    auto ts_xyyyyy_xyyyy = pbuffer.data(idx_ovl_ih + 325);

    auto ts_xyyyyy_xyyyz = pbuffer.data(idx_ovl_ih + 326);

    auto ts_xyyyyy_xyyzz = pbuffer.data(idx_ovl_ih + 327);

    auto ts_xyyyyy_xyzzz = pbuffer.data(idx_ovl_ih + 328);

    auto ts_xyyyyy_xzzzz = pbuffer.data(idx_ovl_ih + 329);

    auto ts_xyyyyy_yyyyy = pbuffer.data(idx_ovl_ih + 330);

    auto ts_xyyyyy_yyyyz = pbuffer.data(idx_ovl_ih + 331);

    auto ts_xyyyyy_yyyzz = pbuffer.data(idx_ovl_ih + 332);

    auto ts_xyyyyy_yyzzz = pbuffer.data(idx_ovl_ih + 333);

    auto ts_xyyyyy_yzzzz = pbuffer.data(idx_ovl_ih + 334);

    auto ts_xyyyyy_zzzzz = pbuffer.data(idx_ovl_ih + 335);

    #pragma omp simd aligned(pa_x, ts_xyyyyy_xxxxx, ts_xyyyyy_xxxxy, ts_xyyyyy_xxxxz, ts_xyyyyy_xxxyy, ts_xyyyyy_xxxyz, ts_xyyyyy_xxxzz, ts_xyyyyy_xxyyy, ts_xyyyyy_xxyyz, ts_xyyyyy_xxyzz, ts_xyyyyy_xxzzz, ts_xyyyyy_xyyyy, ts_xyyyyy_xyyyz, ts_xyyyyy_xyyzz, ts_xyyyyy_xyzzz, ts_xyyyyy_xzzzz, ts_xyyyyy_yyyyy, ts_xyyyyy_yyyyz, ts_xyyyyy_yyyzz, ts_xyyyyy_yyzzz, ts_xyyyyy_yzzzz, ts_xyyyyy_zzzzz, ts_yyyyy_xxxx, ts_yyyyy_xxxxx, ts_yyyyy_xxxxy, ts_yyyyy_xxxxz, ts_yyyyy_xxxy, ts_yyyyy_xxxyy, ts_yyyyy_xxxyz, ts_yyyyy_xxxz, ts_yyyyy_xxxzz, ts_yyyyy_xxyy, ts_yyyyy_xxyyy, ts_yyyyy_xxyyz, ts_yyyyy_xxyz, ts_yyyyy_xxyzz, ts_yyyyy_xxzz, ts_yyyyy_xxzzz, ts_yyyyy_xyyy, ts_yyyyy_xyyyy, ts_yyyyy_xyyyz, ts_yyyyy_xyyz, ts_yyyyy_xyyzz, ts_yyyyy_xyzz, ts_yyyyy_xyzzz, ts_yyyyy_xzzz, ts_yyyyy_xzzzz, ts_yyyyy_yyyy, ts_yyyyy_yyyyy, ts_yyyyy_yyyyz, ts_yyyyy_yyyz, ts_yyyyy_yyyzz, ts_yyyyy_yyzz, ts_yyyyy_yyzzz, ts_yyyyy_yzzz, ts_yyyyy_yzzzz, ts_yyyyy_zzzz, ts_yyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyyy_xxxxx[i] = 5.0 * ts_yyyyy_xxxx[i] * fe_0 + ts_yyyyy_xxxxx[i] * pa_x[i];

        ts_xyyyyy_xxxxy[i] = 4.0 * ts_yyyyy_xxxy[i] * fe_0 + ts_yyyyy_xxxxy[i] * pa_x[i];

        ts_xyyyyy_xxxxz[i] = 4.0 * ts_yyyyy_xxxz[i] * fe_0 + ts_yyyyy_xxxxz[i] * pa_x[i];

        ts_xyyyyy_xxxyy[i] = 3.0 * ts_yyyyy_xxyy[i] * fe_0 + ts_yyyyy_xxxyy[i] * pa_x[i];

        ts_xyyyyy_xxxyz[i] = 3.0 * ts_yyyyy_xxyz[i] * fe_0 + ts_yyyyy_xxxyz[i] * pa_x[i];

        ts_xyyyyy_xxxzz[i] = 3.0 * ts_yyyyy_xxzz[i] * fe_0 + ts_yyyyy_xxxzz[i] * pa_x[i];

        ts_xyyyyy_xxyyy[i] = 2.0 * ts_yyyyy_xyyy[i] * fe_0 + ts_yyyyy_xxyyy[i] * pa_x[i];

        ts_xyyyyy_xxyyz[i] = 2.0 * ts_yyyyy_xyyz[i] * fe_0 + ts_yyyyy_xxyyz[i] * pa_x[i];

        ts_xyyyyy_xxyzz[i] = 2.0 * ts_yyyyy_xyzz[i] * fe_0 + ts_yyyyy_xxyzz[i] * pa_x[i];

        ts_xyyyyy_xxzzz[i] = 2.0 * ts_yyyyy_xzzz[i] * fe_0 + ts_yyyyy_xxzzz[i] * pa_x[i];

        ts_xyyyyy_xyyyy[i] = ts_yyyyy_yyyy[i] * fe_0 + ts_yyyyy_xyyyy[i] * pa_x[i];

        ts_xyyyyy_xyyyz[i] = ts_yyyyy_yyyz[i] * fe_0 + ts_yyyyy_xyyyz[i] * pa_x[i];

        ts_xyyyyy_xyyzz[i] = ts_yyyyy_yyzz[i] * fe_0 + ts_yyyyy_xyyzz[i] * pa_x[i];

        ts_xyyyyy_xyzzz[i] = ts_yyyyy_yzzz[i] * fe_0 + ts_yyyyy_xyzzz[i] * pa_x[i];

        ts_xyyyyy_xzzzz[i] = ts_yyyyy_zzzz[i] * fe_0 + ts_yyyyy_xzzzz[i] * pa_x[i];

        ts_xyyyyy_yyyyy[i] = ts_yyyyy_yyyyy[i] * pa_x[i];

        ts_xyyyyy_yyyyz[i] = ts_yyyyy_yyyyz[i] * pa_x[i];

        ts_xyyyyy_yyyzz[i] = ts_yyyyy_yyyzz[i] * pa_x[i];

        ts_xyyyyy_yyzzz[i] = ts_yyyyy_yyzzz[i] * pa_x[i];

        ts_xyyyyy_yzzzz[i] = ts_yyyyy_yzzzz[i] * pa_x[i];

        ts_xyyyyy_zzzzz[i] = ts_yyyyy_zzzzz[i] * pa_x[i];
    }

    // Set up 336-357 components of targeted buffer : IH

    auto ts_xyyyyz_xxxxx = pbuffer.data(idx_ovl_ih + 336);

    auto ts_xyyyyz_xxxxy = pbuffer.data(idx_ovl_ih + 337);

    auto ts_xyyyyz_xxxxz = pbuffer.data(idx_ovl_ih + 338);

    auto ts_xyyyyz_xxxyy = pbuffer.data(idx_ovl_ih + 339);

    auto ts_xyyyyz_xxxyz = pbuffer.data(idx_ovl_ih + 340);

    auto ts_xyyyyz_xxxzz = pbuffer.data(idx_ovl_ih + 341);

    auto ts_xyyyyz_xxyyy = pbuffer.data(idx_ovl_ih + 342);

    auto ts_xyyyyz_xxyyz = pbuffer.data(idx_ovl_ih + 343);

    auto ts_xyyyyz_xxyzz = pbuffer.data(idx_ovl_ih + 344);

    auto ts_xyyyyz_xxzzz = pbuffer.data(idx_ovl_ih + 345);

    auto ts_xyyyyz_xyyyy = pbuffer.data(idx_ovl_ih + 346);

    auto ts_xyyyyz_xyyyz = pbuffer.data(idx_ovl_ih + 347);

    auto ts_xyyyyz_xyyzz = pbuffer.data(idx_ovl_ih + 348);

    auto ts_xyyyyz_xyzzz = pbuffer.data(idx_ovl_ih + 349);

    auto ts_xyyyyz_xzzzz = pbuffer.data(idx_ovl_ih + 350);

    auto ts_xyyyyz_yyyyy = pbuffer.data(idx_ovl_ih + 351);

    auto ts_xyyyyz_yyyyz = pbuffer.data(idx_ovl_ih + 352);

    auto ts_xyyyyz_yyyzz = pbuffer.data(idx_ovl_ih + 353);

    auto ts_xyyyyz_yyzzz = pbuffer.data(idx_ovl_ih + 354);

    auto ts_xyyyyz_yzzzz = pbuffer.data(idx_ovl_ih + 355);

    auto ts_xyyyyz_zzzzz = pbuffer.data(idx_ovl_ih + 356);

    #pragma omp simd aligned(pa_x, pa_z, ts_xyyyy_xxxxx, ts_xyyyy_xxxxy, ts_xyyyy_xxxyy, ts_xyyyy_xxyyy, ts_xyyyy_xyyyy, ts_xyyyyz_xxxxx, ts_xyyyyz_xxxxy, ts_xyyyyz_xxxxz, ts_xyyyyz_xxxyy, ts_xyyyyz_xxxyz, ts_xyyyyz_xxxzz, ts_xyyyyz_xxyyy, ts_xyyyyz_xxyyz, ts_xyyyyz_xxyzz, ts_xyyyyz_xxzzz, ts_xyyyyz_xyyyy, ts_xyyyyz_xyyyz, ts_xyyyyz_xyyzz, ts_xyyyyz_xyzzz, ts_xyyyyz_xzzzz, ts_xyyyyz_yyyyy, ts_xyyyyz_yyyyz, ts_xyyyyz_yyyzz, ts_xyyyyz_yyzzz, ts_xyyyyz_yzzzz, ts_xyyyyz_zzzzz, ts_yyyyz_xxxxz, ts_yyyyz_xxxyz, ts_yyyyz_xxxz, ts_yyyyz_xxxzz, ts_yyyyz_xxyyz, ts_yyyyz_xxyz, ts_yyyyz_xxyzz, ts_yyyyz_xxzz, ts_yyyyz_xxzzz, ts_yyyyz_xyyyz, ts_yyyyz_xyyz, ts_yyyyz_xyyzz, ts_yyyyz_xyzz, ts_yyyyz_xyzzz, ts_yyyyz_xzzz, ts_yyyyz_xzzzz, ts_yyyyz_yyyyy, ts_yyyyz_yyyyz, ts_yyyyz_yyyz, ts_yyyyz_yyyzz, ts_yyyyz_yyzz, ts_yyyyz_yyzzz, ts_yyyyz_yzzz, ts_yyyyz_yzzzz, ts_yyyyz_zzzz, ts_yyyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyyz_xxxxx[i] = ts_xyyyy_xxxxx[i] * pa_z[i];

        ts_xyyyyz_xxxxy[i] = ts_xyyyy_xxxxy[i] * pa_z[i];

        ts_xyyyyz_xxxxz[i] = 4.0 * ts_yyyyz_xxxz[i] * fe_0 + ts_yyyyz_xxxxz[i] * pa_x[i];

        ts_xyyyyz_xxxyy[i] = ts_xyyyy_xxxyy[i] * pa_z[i];

        ts_xyyyyz_xxxyz[i] = 3.0 * ts_yyyyz_xxyz[i] * fe_0 + ts_yyyyz_xxxyz[i] * pa_x[i];

        ts_xyyyyz_xxxzz[i] = 3.0 * ts_yyyyz_xxzz[i] * fe_0 + ts_yyyyz_xxxzz[i] * pa_x[i];

        ts_xyyyyz_xxyyy[i] = ts_xyyyy_xxyyy[i] * pa_z[i];

        ts_xyyyyz_xxyyz[i] = 2.0 * ts_yyyyz_xyyz[i] * fe_0 + ts_yyyyz_xxyyz[i] * pa_x[i];

        ts_xyyyyz_xxyzz[i] = 2.0 * ts_yyyyz_xyzz[i] * fe_0 + ts_yyyyz_xxyzz[i] * pa_x[i];

        ts_xyyyyz_xxzzz[i] = 2.0 * ts_yyyyz_xzzz[i] * fe_0 + ts_yyyyz_xxzzz[i] * pa_x[i];

        ts_xyyyyz_xyyyy[i] = ts_xyyyy_xyyyy[i] * pa_z[i];

        ts_xyyyyz_xyyyz[i] = ts_yyyyz_yyyz[i] * fe_0 + ts_yyyyz_xyyyz[i] * pa_x[i];

        ts_xyyyyz_xyyzz[i] = ts_yyyyz_yyzz[i] * fe_0 + ts_yyyyz_xyyzz[i] * pa_x[i];

        ts_xyyyyz_xyzzz[i] = ts_yyyyz_yzzz[i] * fe_0 + ts_yyyyz_xyzzz[i] * pa_x[i];

        ts_xyyyyz_xzzzz[i] = ts_yyyyz_zzzz[i] * fe_0 + ts_yyyyz_xzzzz[i] * pa_x[i];

        ts_xyyyyz_yyyyy[i] = ts_yyyyz_yyyyy[i] * pa_x[i];

        ts_xyyyyz_yyyyz[i] = ts_yyyyz_yyyyz[i] * pa_x[i];

        ts_xyyyyz_yyyzz[i] = ts_yyyyz_yyyzz[i] * pa_x[i];

        ts_xyyyyz_yyzzz[i] = ts_yyyyz_yyzzz[i] * pa_x[i];

        ts_xyyyyz_yzzzz[i] = ts_yyyyz_yzzzz[i] * pa_x[i];

        ts_xyyyyz_zzzzz[i] = ts_yyyyz_zzzzz[i] * pa_x[i];
    }

    // Set up 357-378 components of targeted buffer : IH

    auto ts_xyyyzz_xxxxx = pbuffer.data(idx_ovl_ih + 357);

    auto ts_xyyyzz_xxxxy = pbuffer.data(idx_ovl_ih + 358);

    auto ts_xyyyzz_xxxxz = pbuffer.data(idx_ovl_ih + 359);

    auto ts_xyyyzz_xxxyy = pbuffer.data(idx_ovl_ih + 360);

    auto ts_xyyyzz_xxxyz = pbuffer.data(idx_ovl_ih + 361);

    auto ts_xyyyzz_xxxzz = pbuffer.data(idx_ovl_ih + 362);

    auto ts_xyyyzz_xxyyy = pbuffer.data(idx_ovl_ih + 363);

    auto ts_xyyyzz_xxyyz = pbuffer.data(idx_ovl_ih + 364);

    auto ts_xyyyzz_xxyzz = pbuffer.data(idx_ovl_ih + 365);

    auto ts_xyyyzz_xxzzz = pbuffer.data(idx_ovl_ih + 366);

    auto ts_xyyyzz_xyyyy = pbuffer.data(idx_ovl_ih + 367);

    auto ts_xyyyzz_xyyyz = pbuffer.data(idx_ovl_ih + 368);

    auto ts_xyyyzz_xyyzz = pbuffer.data(idx_ovl_ih + 369);

    auto ts_xyyyzz_xyzzz = pbuffer.data(idx_ovl_ih + 370);

    auto ts_xyyyzz_xzzzz = pbuffer.data(idx_ovl_ih + 371);

    auto ts_xyyyzz_yyyyy = pbuffer.data(idx_ovl_ih + 372);

    auto ts_xyyyzz_yyyyz = pbuffer.data(idx_ovl_ih + 373);

    auto ts_xyyyzz_yyyzz = pbuffer.data(idx_ovl_ih + 374);

    auto ts_xyyyzz_yyzzz = pbuffer.data(idx_ovl_ih + 375);

    auto ts_xyyyzz_yzzzz = pbuffer.data(idx_ovl_ih + 376);

    auto ts_xyyyzz_zzzzz = pbuffer.data(idx_ovl_ih + 377);

    #pragma omp simd aligned(pa_x, ts_xyyyzz_xxxxx, ts_xyyyzz_xxxxy, ts_xyyyzz_xxxxz, ts_xyyyzz_xxxyy, ts_xyyyzz_xxxyz, ts_xyyyzz_xxxzz, ts_xyyyzz_xxyyy, ts_xyyyzz_xxyyz, ts_xyyyzz_xxyzz, ts_xyyyzz_xxzzz, ts_xyyyzz_xyyyy, ts_xyyyzz_xyyyz, ts_xyyyzz_xyyzz, ts_xyyyzz_xyzzz, ts_xyyyzz_xzzzz, ts_xyyyzz_yyyyy, ts_xyyyzz_yyyyz, ts_xyyyzz_yyyzz, ts_xyyyzz_yyzzz, ts_xyyyzz_yzzzz, ts_xyyyzz_zzzzz, ts_yyyzz_xxxx, ts_yyyzz_xxxxx, ts_yyyzz_xxxxy, ts_yyyzz_xxxxz, ts_yyyzz_xxxy, ts_yyyzz_xxxyy, ts_yyyzz_xxxyz, ts_yyyzz_xxxz, ts_yyyzz_xxxzz, ts_yyyzz_xxyy, ts_yyyzz_xxyyy, ts_yyyzz_xxyyz, ts_yyyzz_xxyz, ts_yyyzz_xxyzz, ts_yyyzz_xxzz, ts_yyyzz_xxzzz, ts_yyyzz_xyyy, ts_yyyzz_xyyyy, ts_yyyzz_xyyyz, ts_yyyzz_xyyz, ts_yyyzz_xyyzz, ts_yyyzz_xyzz, ts_yyyzz_xyzzz, ts_yyyzz_xzzz, ts_yyyzz_xzzzz, ts_yyyzz_yyyy, ts_yyyzz_yyyyy, ts_yyyzz_yyyyz, ts_yyyzz_yyyz, ts_yyyzz_yyyzz, ts_yyyzz_yyzz, ts_yyyzz_yyzzz, ts_yyyzz_yzzz, ts_yyyzz_yzzzz, ts_yyyzz_zzzz, ts_yyyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyzz_xxxxx[i] = 5.0 * ts_yyyzz_xxxx[i] * fe_0 + ts_yyyzz_xxxxx[i] * pa_x[i];

        ts_xyyyzz_xxxxy[i] = 4.0 * ts_yyyzz_xxxy[i] * fe_0 + ts_yyyzz_xxxxy[i] * pa_x[i];

        ts_xyyyzz_xxxxz[i] = 4.0 * ts_yyyzz_xxxz[i] * fe_0 + ts_yyyzz_xxxxz[i] * pa_x[i];

        ts_xyyyzz_xxxyy[i] = 3.0 * ts_yyyzz_xxyy[i] * fe_0 + ts_yyyzz_xxxyy[i] * pa_x[i];

        ts_xyyyzz_xxxyz[i] = 3.0 * ts_yyyzz_xxyz[i] * fe_0 + ts_yyyzz_xxxyz[i] * pa_x[i];

        ts_xyyyzz_xxxzz[i] = 3.0 * ts_yyyzz_xxzz[i] * fe_0 + ts_yyyzz_xxxzz[i] * pa_x[i];

        ts_xyyyzz_xxyyy[i] = 2.0 * ts_yyyzz_xyyy[i] * fe_0 + ts_yyyzz_xxyyy[i] * pa_x[i];

        ts_xyyyzz_xxyyz[i] = 2.0 * ts_yyyzz_xyyz[i] * fe_0 + ts_yyyzz_xxyyz[i] * pa_x[i];

        ts_xyyyzz_xxyzz[i] = 2.0 * ts_yyyzz_xyzz[i] * fe_0 + ts_yyyzz_xxyzz[i] * pa_x[i];

        ts_xyyyzz_xxzzz[i] = 2.0 * ts_yyyzz_xzzz[i] * fe_0 + ts_yyyzz_xxzzz[i] * pa_x[i];

        ts_xyyyzz_xyyyy[i] = ts_yyyzz_yyyy[i] * fe_0 + ts_yyyzz_xyyyy[i] * pa_x[i];

        ts_xyyyzz_xyyyz[i] = ts_yyyzz_yyyz[i] * fe_0 + ts_yyyzz_xyyyz[i] * pa_x[i];

        ts_xyyyzz_xyyzz[i] = ts_yyyzz_yyzz[i] * fe_0 + ts_yyyzz_xyyzz[i] * pa_x[i];

        ts_xyyyzz_xyzzz[i] = ts_yyyzz_yzzz[i] * fe_0 + ts_yyyzz_xyzzz[i] * pa_x[i];

        ts_xyyyzz_xzzzz[i] = ts_yyyzz_zzzz[i] * fe_0 + ts_yyyzz_xzzzz[i] * pa_x[i];

        ts_xyyyzz_yyyyy[i] = ts_yyyzz_yyyyy[i] * pa_x[i];

        ts_xyyyzz_yyyyz[i] = ts_yyyzz_yyyyz[i] * pa_x[i];

        ts_xyyyzz_yyyzz[i] = ts_yyyzz_yyyzz[i] * pa_x[i];

        ts_xyyyzz_yyzzz[i] = ts_yyyzz_yyzzz[i] * pa_x[i];

        ts_xyyyzz_yzzzz[i] = ts_yyyzz_yzzzz[i] * pa_x[i];

        ts_xyyyzz_zzzzz[i] = ts_yyyzz_zzzzz[i] * pa_x[i];
    }

    // Set up 378-399 components of targeted buffer : IH

    auto ts_xyyzzz_xxxxx = pbuffer.data(idx_ovl_ih + 378);

    auto ts_xyyzzz_xxxxy = pbuffer.data(idx_ovl_ih + 379);

    auto ts_xyyzzz_xxxxz = pbuffer.data(idx_ovl_ih + 380);

    auto ts_xyyzzz_xxxyy = pbuffer.data(idx_ovl_ih + 381);

    auto ts_xyyzzz_xxxyz = pbuffer.data(idx_ovl_ih + 382);

    auto ts_xyyzzz_xxxzz = pbuffer.data(idx_ovl_ih + 383);

    auto ts_xyyzzz_xxyyy = pbuffer.data(idx_ovl_ih + 384);

    auto ts_xyyzzz_xxyyz = pbuffer.data(idx_ovl_ih + 385);

    auto ts_xyyzzz_xxyzz = pbuffer.data(idx_ovl_ih + 386);

    auto ts_xyyzzz_xxzzz = pbuffer.data(idx_ovl_ih + 387);

    auto ts_xyyzzz_xyyyy = pbuffer.data(idx_ovl_ih + 388);

    auto ts_xyyzzz_xyyyz = pbuffer.data(idx_ovl_ih + 389);

    auto ts_xyyzzz_xyyzz = pbuffer.data(idx_ovl_ih + 390);

    auto ts_xyyzzz_xyzzz = pbuffer.data(idx_ovl_ih + 391);

    auto ts_xyyzzz_xzzzz = pbuffer.data(idx_ovl_ih + 392);

    auto ts_xyyzzz_yyyyy = pbuffer.data(idx_ovl_ih + 393);

    auto ts_xyyzzz_yyyyz = pbuffer.data(idx_ovl_ih + 394);

    auto ts_xyyzzz_yyyzz = pbuffer.data(idx_ovl_ih + 395);

    auto ts_xyyzzz_yyzzz = pbuffer.data(idx_ovl_ih + 396);

    auto ts_xyyzzz_yzzzz = pbuffer.data(idx_ovl_ih + 397);

    auto ts_xyyzzz_zzzzz = pbuffer.data(idx_ovl_ih + 398);

    #pragma omp simd aligned(pa_x, ts_xyyzzz_xxxxx, ts_xyyzzz_xxxxy, ts_xyyzzz_xxxxz, ts_xyyzzz_xxxyy, ts_xyyzzz_xxxyz, ts_xyyzzz_xxxzz, ts_xyyzzz_xxyyy, ts_xyyzzz_xxyyz, ts_xyyzzz_xxyzz, ts_xyyzzz_xxzzz, ts_xyyzzz_xyyyy, ts_xyyzzz_xyyyz, ts_xyyzzz_xyyzz, ts_xyyzzz_xyzzz, ts_xyyzzz_xzzzz, ts_xyyzzz_yyyyy, ts_xyyzzz_yyyyz, ts_xyyzzz_yyyzz, ts_xyyzzz_yyzzz, ts_xyyzzz_yzzzz, ts_xyyzzz_zzzzz, ts_yyzzz_xxxx, ts_yyzzz_xxxxx, ts_yyzzz_xxxxy, ts_yyzzz_xxxxz, ts_yyzzz_xxxy, ts_yyzzz_xxxyy, ts_yyzzz_xxxyz, ts_yyzzz_xxxz, ts_yyzzz_xxxzz, ts_yyzzz_xxyy, ts_yyzzz_xxyyy, ts_yyzzz_xxyyz, ts_yyzzz_xxyz, ts_yyzzz_xxyzz, ts_yyzzz_xxzz, ts_yyzzz_xxzzz, ts_yyzzz_xyyy, ts_yyzzz_xyyyy, ts_yyzzz_xyyyz, ts_yyzzz_xyyz, ts_yyzzz_xyyzz, ts_yyzzz_xyzz, ts_yyzzz_xyzzz, ts_yyzzz_xzzz, ts_yyzzz_xzzzz, ts_yyzzz_yyyy, ts_yyzzz_yyyyy, ts_yyzzz_yyyyz, ts_yyzzz_yyyz, ts_yyzzz_yyyzz, ts_yyzzz_yyzz, ts_yyzzz_yyzzz, ts_yyzzz_yzzz, ts_yyzzz_yzzzz, ts_yyzzz_zzzz, ts_yyzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyzzz_xxxxx[i] = 5.0 * ts_yyzzz_xxxx[i] * fe_0 + ts_yyzzz_xxxxx[i] * pa_x[i];

        ts_xyyzzz_xxxxy[i] = 4.0 * ts_yyzzz_xxxy[i] * fe_0 + ts_yyzzz_xxxxy[i] * pa_x[i];

        ts_xyyzzz_xxxxz[i] = 4.0 * ts_yyzzz_xxxz[i] * fe_0 + ts_yyzzz_xxxxz[i] * pa_x[i];

        ts_xyyzzz_xxxyy[i] = 3.0 * ts_yyzzz_xxyy[i] * fe_0 + ts_yyzzz_xxxyy[i] * pa_x[i];

        ts_xyyzzz_xxxyz[i] = 3.0 * ts_yyzzz_xxyz[i] * fe_0 + ts_yyzzz_xxxyz[i] * pa_x[i];

        ts_xyyzzz_xxxzz[i] = 3.0 * ts_yyzzz_xxzz[i] * fe_0 + ts_yyzzz_xxxzz[i] * pa_x[i];

        ts_xyyzzz_xxyyy[i] = 2.0 * ts_yyzzz_xyyy[i] * fe_0 + ts_yyzzz_xxyyy[i] * pa_x[i];

        ts_xyyzzz_xxyyz[i] = 2.0 * ts_yyzzz_xyyz[i] * fe_0 + ts_yyzzz_xxyyz[i] * pa_x[i];

        ts_xyyzzz_xxyzz[i] = 2.0 * ts_yyzzz_xyzz[i] * fe_0 + ts_yyzzz_xxyzz[i] * pa_x[i];

        ts_xyyzzz_xxzzz[i] = 2.0 * ts_yyzzz_xzzz[i] * fe_0 + ts_yyzzz_xxzzz[i] * pa_x[i];

        ts_xyyzzz_xyyyy[i] = ts_yyzzz_yyyy[i] * fe_0 + ts_yyzzz_xyyyy[i] * pa_x[i];

        ts_xyyzzz_xyyyz[i] = ts_yyzzz_yyyz[i] * fe_0 + ts_yyzzz_xyyyz[i] * pa_x[i];

        ts_xyyzzz_xyyzz[i] = ts_yyzzz_yyzz[i] * fe_0 + ts_yyzzz_xyyzz[i] * pa_x[i];

        ts_xyyzzz_xyzzz[i] = ts_yyzzz_yzzz[i] * fe_0 + ts_yyzzz_xyzzz[i] * pa_x[i];

        ts_xyyzzz_xzzzz[i] = ts_yyzzz_zzzz[i] * fe_0 + ts_yyzzz_xzzzz[i] * pa_x[i];

        ts_xyyzzz_yyyyy[i] = ts_yyzzz_yyyyy[i] * pa_x[i];

        ts_xyyzzz_yyyyz[i] = ts_yyzzz_yyyyz[i] * pa_x[i];

        ts_xyyzzz_yyyzz[i] = ts_yyzzz_yyyzz[i] * pa_x[i];

        ts_xyyzzz_yyzzz[i] = ts_yyzzz_yyzzz[i] * pa_x[i];

        ts_xyyzzz_yzzzz[i] = ts_yyzzz_yzzzz[i] * pa_x[i];

        ts_xyyzzz_zzzzz[i] = ts_yyzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 399-420 components of targeted buffer : IH

    auto ts_xyzzzz_xxxxx = pbuffer.data(idx_ovl_ih + 399);

    auto ts_xyzzzz_xxxxy = pbuffer.data(idx_ovl_ih + 400);

    auto ts_xyzzzz_xxxxz = pbuffer.data(idx_ovl_ih + 401);

    auto ts_xyzzzz_xxxyy = pbuffer.data(idx_ovl_ih + 402);

    auto ts_xyzzzz_xxxyz = pbuffer.data(idx_ovl_ih + 403);

    auto ts_xyzzzz_xxxzz = pbuffer.data(idx_ovl_ih + 404);

    auto ts_xyzzzz_xxyyy = pbuffer.data(idx_ovl_ih + 405);

    auto ts_xyzzzz_xxyyz = pbuffer.data(idx_ovl_ih + 406);

    auto ts_xyzzzz_xxyzz = pbuffer.data(idx_ovl_ih + 407);

    auto ts_xyzzzz_xxzzz = pbuffer.data(idx_ovl_ih + 408);

    auto ts_xyzzzz_xyyyy = pbuffer.data(idx_ovl_ih + 409);

    auto ts_xyzzzz_xyyyz = pbuffer.data(idx_ovl_ih + 410);

    auto ts_xyzzzz_xyyzz = pbuffer.data(idx_ovl_ih + 411);

    auto ts_xyzzzz_xyzzz = pbuffer.data(idx_ovl_ih + 412);

    auto ts_xyzzzz_xzzzz = pbuffer.data(idx_ovl_ih + 413);

    auto ts_xyzzzz_yyyyy = pbuffer.data(idx_ovl_ih + 414);

    auto ts_xyzzzz_yyyyz = pbuffer.data(idx_ovl_ih + 415);

    auto ts_xyzzzz_yyyzz = pbuffer.data(idx_ovl_ih + 416);

    auto ts_xyzzzz_yyzzz = pbuffer.data(idx_ovl_ih + 417);

    auto ts_xyzzzz_yzzzz = pbuffer.data(idx_ovl_ih + 418);

    auto ts_xyzzzz_zzzzz = pbuffer.data(idx_ovl_ih + 419);

    #pragma omp simd aligned(pa_x, pa_y, ts_xyzzzz_xxxxx, ts_xyzzzz_xxxxy, ts_xyzzzz_xxxxz, ts_xyzzzz_xxxyy, ts_xyzzzz_xxxyz, ts_xyzzzz_xxxzz, ts_xyzzzz_xxyyy, ts_xyzzzz_xxyyz, ts_xyzzzz_xxyzz, ts_xyzzzz_xxzzz, ts_xyzzzz_xyyyy, ts_xyzzzz_xyyyz, ts_xyzzzz_xyyzz, ts_xyzzzz_xyzzz, ts_xyzzzz_xzzzz, ts_xyzzzz_yyyyy, ts_xyzzzz_yyyyz, ts_xyzzzz_yyyzz, ts_xyzzzz_yyzzz, ts_xyzzzz_yzzzz, ts_xyzzzz_zzzzz, ts_xzzzz_xxxxx, ts_xzzzz_xxxxz, ts_xzzzz_xxxzz, ts_xzzzz_xxzzz, ts_xzzzz_xzzzz, ts_yzzzz_xxxxy, ts_yzzzz_xxxy, ts_yzzzz_xxxyy, ts_yzzzz_xxxyz, ts_yzzzz_xxyy, ts_yzzzz_xxyyy, ts_yzzzz_xxyyz, ts_yzzzz_xxyz, ts_yzzzz_xxyzz, ts_yzzzz_xyyy, ts_yzzzz_xyyyy, ts_yzzzz_xyyyz, ts_yzzzz_xyyz, ts_yzzzz_xyyzz, ts_yzzzz_xyzz, ts_yzzzz_xyzzz, ts_yzzzz_yyyy, ts_yzzzz_yyyyy, ts_yzzzz_yyyyz, ts_yzzzz_yyyz, ts_yzzzz_yyyzz, ts_yzzzz_yyzz, ts_yzzzz_yyzzz, ts_yzzzz_yzzz, ts_yzzzz_yzzzz, ts_yzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyzzzz_xxxxx[i] = ts_xzzzz_xxxxx[i] * pa_y[i];

        ts_xyzzzz_xxxxy[i] = 4.0 * ts_yzzzz_xxxy[i] * fe_0 + ts_yzzzz_xxxxy[i] * pa_x[i];

        ts_xyzzzz_xxxxz[i] = ts_xzzzz_xxxxz[i] * pa_y[i];

        ts_xyzzzz_xxxyy[i] = 3.0 * ts_yzzzz_xxyy[i] * fe_0 + ts_yzzzz_xxxyy[i] * pa_x[i];

        ts_xyzzzz_xxxyz[i] = 3.0 * ts_yzzzz_xxyz[i] * fe_0 + ts_yzzzz_xxxyz[i] * pa_x[i];

        ts_xyzzzz_xxxzz[i] = ts_xzzzz_xxxzz[i] * pa_y[i];

        ts_xyzzzz_xxyyy[i] = 2.0 * ts_yzzzz_xyyy[i] * fe_0 + ts_yzzzz_xxyyy[i] * pa_x[i];

        ts_xyzzzz_xxyyz[i] = 2.0 * ts_yzzzz_xyyz[i] * fe_0 + ts_yzzzz_xxyyz[i] * pa_x[i];

        ts_xyzzzz_xxyzz[i] = 2.0 * ts_yzzzz_xyzz[i] * fe_0 + ts_yzzzz_xxyzz[i] * pa_x[i];

        ts_xyzzzz_xxzzz[i] = ts_xzzzz_xxzzz[i] * pa_y[i];

        ts_xyzzzz_xyyyy[i] = ts_yzzzz_yyyy[i] * fe_0 + ts_yzzzz_xyyyy[i] * pa_x[i];

        ts_xyzzzz_xyyyz[i] = ts_yzzzz_yyyz[i] * fe_0 + ts_yzzzz_xyyyz[i] * pa_x[i];

        ts_xyzzzz_xyyzz[i] = ts_yzzzz_yyzz[i] * fe_0 + ts_yzzzz_xyyzz[i] * pa_x[i];

        ts_xyzzzz_xyzzz[i] = ts_yzzzz_yzzz[i] * fe_0 + ts_yzzzz_xyzzz[i] * pa_x[i];

        ts_xyzzzz_xzzzz[i] = ts_xzzzz_xzzzz[i] * pa_y[i];

        ts_xyzzzz_yyyyy[i] = ts_yzzzz_yyyyy[i] * pa_x[i];

        ts_xyzzzz_yyyyz[i] = ts_yzzzz_yyyyz[i] * pa_x[i];

        ts_xyzzzz_yyyzz[i] = ts_yzzzz_yyyzz[i] * pa_x[i];

        ts_xyzzzz_yyzzz[i] = ts_yzzzz_yyzzz[i] * pa_x[i];

        ts_xyzzzz_yzzzz[i] = ts_yzzzz_yzzzz[i] * pa_x[i];

        ts_xyzzzz_zzzzz[i] = ts_yzzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 420-441 components of targeted buffer : IH

    auto ts_xzzzzz_xxxxx = pbuffer.data(idx_ovl_ih + 420);

    auto ts_xzzzzz_xxxxy = pbuffer.data(idx_ovl_ih + 421);

    auto ts_xzzzzz_xxxxz = pbuffer.data(idx_ovl_ih + 422);

    auto ts_xzzzzz_xxxyy = pbuffer.data(idx_ovl_ih + 423);

    auto ts_xzzzzz_xxxyz = pbuffer.data(idx_ovl_ih + 424);

    auto ts_xzzzzz_xxxzz = pbuffer.data(idx_ovl_ih + 425);

    auto ts_xzzzzz_xxyyy = pbuffer.data(idx_ovl_ih + 426);

    auto ts_xzzzzz_xxyyz = pbuffer.data(idx_ovl_ih + 427);

    auto ts_xzzzzz_xxyzz = pbuffer.data(idx_ovl_ih + 428);

    auto ts_xzzzzz_xxzzz = pbuffer.data(idx_ovl_ih + 429);

    auto ts_xzzzzz_xyyyy = pbuffer.data(idx_ovl_ih + 430);

    auto ts_xzzzzz_xyyyz = pbuffer.data(idx_ovl_ih + 431);

    auto ts_xzzzzz_xyyzz = pbuffer.data(idx_ovl_ih + 432);

    auto ts_xzzzzz_xyzzz = pbuffer.data(idx_ovl_ih + 433);

    auto ts_xzzzzz_xzzzz = pbuffer.data(idx_ovl_ih + 434);

    auto ts_xzzzzz_yyyyy = pbuffer.data(idx_ovl_ih + 435);

    auto ts_xzzzzz_yyyyz = pbuffer.data(idx_ovl_ih + 436);

    auto ts_xzzzzz_yyyzz = pbuffer.data(idx_ovl_ih + 437);

    auto ts_xzzzzz_yyzzz = pbuffer.data(idx_ovl_ih + 438);

    auto ts_xzzzzz_yzzzz = pbuffer.data(idx_ovl_ih + 439);

    auto ts_xzzzzz_zzzzz = pbuffer.data(idx_ovl_ih + 440);

    #pragma omp simd aligned(pa_x, ts_xzzzzz_xxxxx, ts_xzzzzz_xxxxy, ts_xzzzzz_xxxxz, ts_xzzzzz_xxxyy, ts_xzzzzz_xxxyz, ts_xzzzzz_xxxzz, ts_xzzzzz_xxyyy, ts_xzzzzz_xxyyz, ts_xzzzzz_xxyzz, ts_xzzzzz_xxzzz, ts_xzzzzz_xyyyy, ts_xzzzzz_xyyyz, ts_xzzzzz_xyyzz, ts_xzzzzz_xyzzz, ts_xzzzzz_xzzzz, ts_xzzzzz_yyyyy, ts_xzzzzz_yyyyz, ts_xzzzzz_yyyzz, ts_xzzzzz_yyzzz, ts_xzzzzz_yzzzz, ts_xzzzzz_zzzzz, ts_zzzzz_xxxx, ts_zzzzz_xxxxx, ts_zzzzz_xxxxy, ts_zzzzz_xxxxz, ts_zzzzz_xxxy, ts_zzzzz_xxxyy, ts_zzzzz_xxxyz, ts_zzzzz_xxxz, ts_zzzzz_xxxzz, ts_zzzzz_xxyy, ts_zzzzz_xxyyy, ts_zzzzz_xxyyz, ts_zzzzz_xxyz, ts_zzzzz_xxyzz, ts_zzzzz_xxzz, ts_zzzzz_xxzzz, ts_zzzzz_xyyy, ts_zzzzz_xyyyy, ts_zzzzz_xyyyz, ts_zzzzz_xyyz, ts_zzzzz_xyyzz, ts_zzzzz_xyzz, ts_zzzzz_xyzzz, ts_zzzzz_xzzz, ts_zzzzz_xzzzz, ts_zzzzz_yyyy, ts_zzzzz_yyyyy, ts_zzzzz_yyyyz, ts_zzzzz_yyyz, ts_zzzzz_yyyzz, ts_zzzzz_yyzz, ts_zzzzz_yyzzz, ts_zzzzz_yzzz, ts_zzzzz_yzzzz, ts_zzzzz_zzzz, ts_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzzzz_xxxxx[i] = 5.0 * ts_zzzzz_xxxx[i] * fe_0 + ts_zzzzz_xxxxx[i] * pa_x[i];

        ts_xzzzzz_xxxxy[i] = 4.0 * ts_zzzzz_xxxy[i] * fe_0 + ts_zzzzz_xxxxy[i] * pa_x[i];

        ts_xzzzzz_xxxxz[i] = 4.0 * ts_zzzzz_xxxz[i] * fe_0 + ts_zzzzz_xxxxz[i] * pa_x[i];

        ts_xzzzzz_xxxyy[i] = 3.0 * ts_zzzzz_xxyy[i] * fe_0 + ts_zzzzz_xxxyy[i] * pa_x[i];

        ts_xzzzzz_xxxyz[i] = 3.0 * ts_zzzzz_xxyz[i] * fe_0 + ts_zzzzz_xxxyz[i] * pa_x[i];

        ts_xzzzzz_xxxzz[i] = 3.0 * ts_zzzzz_xxzz[i] * fe_0 + ts_zzzzz_xxxzz[i] * pa_x[i];

        ts_xzzzzz_xxyyy[i] = 2.0 * ts_zzzzz_xyyy[i] * fe_0 + ts_zzzzz_xxyyy[i] * pa_x[i];

        ts_xzzzzz_xxyyz[i] = 2.0 * ts_zzzzz_xyyz[i] * fe_0 + ts_zzzzz_xxyyz[i] * pa_x[i];

        ts_xzzzzz_xxyzz[i] = 2.0 * ts_zzzzz_xyzz[i] * fe_0 + ts_zzzzz_xxyzz[i] * pa_x[i];

        ts_xzzzzz_xxzzz[i] = 2.0 * ts_zzzzz_xzzz[i] * fe_0 + ts_zzzzz_xxzzz[i] * pa_x[i];

        ts_xzzzzz_xyyyy[i] = ts_zzzzz_yyyy[i] * fe_0 + ts_zzzzz_xyyyy[i] * pa_x[i];

        ts_xzzzzz_xyyyz[i] = ts_zzzzz_yyyz[i] * fe_0 + ts_zzzzz_xyyyz[i] * pa_x[i];

        ts_xzzzzz_xyyzz[i] = ts_zzzzz_yyzz[i] * fe_0 + ts_zzzzz_xyyzz[i] * pa_x[i];

        ts_xzzzzz_xyzzz[i] = ts_zzzzz_yzzz[i] * fe_0 + ts_zzzzz_xyzzz[i] * pa_x[i];

        ts_xzzzzz_xzzzz[i] = ts_zzzzz_zzzz[i] * fe_0 + ts_zzzzz_xzzzz[i] * pa_x[i];

        ts_xzzzzz_yyyyy[i] = ts_zzzzz_yyyyy[i] * pa_x[i];

        ts_xzzzzz_yyyyz[i] = ts_zzzzz_yyyyz[i] * pa_x[i];

        ts_xzzzzz_yyyzz[i] = ts_zzzzz_yyyzz[i] * pa_x[i];

        ts_xzzzzz_yyzzz[i] = ts_zzzzz_yyzzz[i] * pa_x[i];

        ts_xzzzzz_yzzzz[i] = ts_zzzzz_yzzzz[i] * pa_x[i];

        ts_xzzzzz_zzzzz[i] = ts_zzzzz_zzzzz[i] * pa_x[i];
    }

    // Set up 441-462 components of targeted buffer : IH

    auto ts_yyyyyy_xxxxx = pbuffer.data(idx_ovl_ih + 441);

    auto ts_yyyyyy_xxxxy = pbuffer.data(idx_ovl_ih + 442);

    auto ts_yyyyyy_xxxxz = pbuffer.data(idx_ovl_ih + 443);

    auto ts_yyyyyy_xxxyy = pbuffer.data(idx_ovl_ih + 444);

    auto ts_yyyyyy_xxxyz = pbuffer.data(idx_ovl_ih + 445);

    auto ts_yyyyyy_xxxzz = pbuffer.data(idx_ovl_ih + 446);

    auto ts_yyyyyy_xxyyy = pbuffer.data(idx_ovl_ih + 447);

    auto ts_yyyyyy_xxyyz = pbuffer.data(idx_ovl_ih + 448);

    auto ts_yyyyyy_xxyzz = pbuffer.data(idx_ovl_ih + 449);

    auto ts_yyyyyy_xxzzz = pbuffer.data(idx_ovl_ih + 450);

    auto ts_yyyyyy_xyyyy = pbuffer.data(idx_ovl_ih + 451);

    auto ts_yyyyyy_xyyyz = pbuffer.data(idx_ovl_ih + 452);

    auto ts_yyyyyy_xyyzz = pbuffer.data(idx_ovl_ih + 453);

    auto ts_yyyyyy_xyzzz = pbuffer.data(idx_ovl_ih + 454);

    auto ts_yyyyyy_xzzzz = pbuffer.data(idx_ovl_ih + 455);

    auto ts_yyyyyy_yyyyy = pbuffer.data(idx_ovl_ih + 456);

    auto ts_yyyyyy_yyyyz = pbuffer.data(idx_ovl_ih + 457);

    auto ts_yyyyyy_yyyzz = pbuffer.data(idx_ovl_ih + 458);

    auto ts_yyyyyy_yyzzz = pbuffer.data(idx_ovl_ih + 459);

    auto ts_yyyyyy_yzzzz = pbuffer.data(idx_ovl_ih + 460);

    auto ts_yyyyyy_zzzzz = pbuffer.data(idx_ovl_ih + 461);

    #pragma omp simd aligned(pa_y, ts_yyyy_xxxxx, ts_yyyy_xxxxy, ts_yyyy_xxxxz, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxxzz, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyzz, ts_yyyy_xxzzz, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyzz, ts_yyyy_xyzzz, ts_yyyy_xzzzz, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyzz, ts_yyyy_yyzzz, ts_yyyy_yzzzz, ts_yyyy_zzzzz, ts_yyyyy_xxxx, ts_yyyyy_xxxxx, ts_yyyyy_xxxxy, ts_yyyyy_xxxxz, ts_yyyyy_xxxy, ts_yyyyy_xxxyy, ts_yyyyy_xxxyz, ts_yyyyy_xxxz, ts_yyyyy_xxxzz, ts_yyyyy_xxyy, ts_yyyyy_xxyyy, ts_yyyyy_xxyyz, ts_yyyyy_xxyz, ts_yyyyy_xxyzz, ts_yyyyy_xxzz, ts_yyyyy_xxzzz, ts_yyyyy_xyyy, ts_yyyyy_xyyyy, ts_yyyyy_xyyyz, ts_yyyyy_xyyz, ts_yyyyy_xyyzz, ts_yyyyy_xyzz, ts_yyyyy_xyzzz, ts_yyyyy_xzzz, ts_yyyyy_xzzzz, ts_yyyyy_yyyy, ts_yyyyy_yyyyy, ts_yyyyy_yyyyz, ts_yyyyy_yyyz, ts_yyyyy_yyyzz, ts_yyyyy_yyzz, ts_yyyyy_yyzzz, ts_yyyyy_yzzz, ts_yyyyy_yzzzz, ts_yyyyy_zzzz, ts_yyyyy_zzzzz, ts_yyyyyy_xxxxx, ts_yyyyyy_xxxxy, ts_yyyyyy_xxxxz, ts_yyyyyy_xxxyy, ts_yyyyyy_xxxyz, ts_yyyyyy_xxxzz, ts_yyyyyy_xxyyy, ts_yyyyyy_xxyyz, ts_yyyyyy_xxyzz, ts_yyyyyy_xxzzz, ts_yyyyyy_xyyyy, ts_yyyyyy_xyyyz, ts_yyyyyy_xyyzz, ts_yyyyyy_xyzzz, ts_yyyyyy_xzzzz, ts_yyyyyy_yyyyy, ts_yyyyyy_yyyyz, ts_yyyyyy_yyyzz, ts_yyyyyy_yyzzz, ts_yyyyyy_yzzzz, ts_yyyyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyyy_xxxxx[i] = 5.0 * ts_yyyy_xxxxx[i] * fe_0 + ts_yyyyy_xxxxx[i] * pa_y[i];

        ts_yyyyyy_xxxxy[i] = 5.0 * ts_yyyy_xxxxy[i] * fe_0 + ts_yyyyy_xxxx[i] * fe_0 + ts_yyyyy_xxxxy[i] * pa_y[i];

        ts_yyyyyy_xxxxz[i] = 5.0 * ts_yyyy_xxxxz[i] * fe_0 + ts_yyyyy_xxxxz[i] * pa_y[i];

        ts_yyyyyy_xxxyy[i] = 5.0 * ts_yyyy_xxxyy[i] * fe_0 + 2.0 * ts_yyyyy_xxxy[i] * fe_0 + ts_yyyyy_xxxyy[i] * pa_y[i];

        ts_yyyyyy_xxxyz[i] = 5.0 * ts_yyyy_xxxyz[i] * fe_0 + ts_yyyyy_xxxz[i] * fe_0 + ts_yyyyy_xxxyz[i] * pa_y[i];

        ts_yyyyyy_xxxzz[i] = 5.0 * ts_yyyy_xxxzz[i] * fe_0 + ts_yyyyy_xxxzz[i] * pa_y[i];

        ts_yyyyyy_xxyyy[i] = 5.0 * ts_yyyy_xxyyy[i] * fe_0 + 3.0 * ts_yyyyy_xxyy[i] * fe_0 + ts_yyyyy_xxyyy[i] * pa_y[i];

        ts_yyyyyy_xxyyz[i] = 5.0 * ts_yyyy_xxyyz[i] * fe_0 + 2.0 * ts_yyyyy_xxyz[i] * fe_0 + ts_yyyyy_xxyyz[i] * pa_y[i];

        ts_yyyyyy_xxyzz[i] = 5.0 * ts_yyyy_xxyzz[i] * fe_0 + ts_yyyyy_xxzz[i] * fe_0 + ts_yyyyy_xxyzz[i] * pa_y[i];

        ts_yyyyyy_xxzzz[i] = 5.0 * ts_yyyy_xxzzz[i] * fe_0 + ts_yyyyy_xxzzz[i] * pa_y[i];

        ts_yyyyyy_xyyyy[i] = 5.0 * ts_yyyy_xyyyy[i] * fe_0 + 4.0 * ts_yyyyy_xyyy[i] * fe_0 + ts_yyyyy_xyyyy[i] * pa_y[i];

        ts_yyyyyy_xyyyz[i] = 5.0 * ts_yyyy_xyyyz[i] * fe_0 + 3.0 * ts_yyyyy_xyyz[i] * fe_0 + ts_yyyyy_xyyyz[i] * pa_y[i];

        ts_yyyyyy_xyyzz[i] = 5.0 * ts_yyyy_xyyzz[i] * fe_0 + 2.0 * ts_yyyyy_xyzz[i] * fe_0 + ts_yyyyy_xyyzz[i] * pa_y[i];

        ts_yyyyyy_xyzzz[i] = 5.0 * ts_yyyy_xyzzz[i] * fe_0 + ts_yyyyy_xzzz[i] * fe_0 + ts_yyyyy_xyzzz[i] * pa_y[i];

        ts_yyyyyy_xzzzz[i] = 5.0 * ts_yyyy_xzzzz[i] * fe_0 + ts_yyyyy_xzzzz[i] * pa_y[i];

        ts_yyyyyy_yyyyy[i] = 5.0 * ts_yyyy_yyyyy[i] * fe_0 + 5.0 * ts_yyyyy_yyyy[i] * fe_0 + ts_yyyyy_yyyyy[i] * pa_y[i];

        ts_yyyyyy_yyyyz[i] = 5.0 * ts_yyyy_yyyyz[i] * fe_0 + 4.0 * ts_yyyyy_yyyz[i] * fe_0 + ts_yyyyy_yyyyz[i] * pa_y[i];

        ts_yyyyyy_yyyzz[i] = 5.0 * ts_yyyy_yyyzz[i] * fe_0 + 3.0 * ts_yyyyy_yyzz[i] * fe_0 + ts_yyyyy_yyyzz[i] * pa_y[i];

        ts_yyyyyy_yyzzz[i] = 5.0 * ts_yyyy_yyzzz[i] * fe_0 + 2.0 * ts_yyyyy_yzzz[i] * fe_0 + ts_yyyyy_yyzzz[i] * pa_y[i];

        ts_yyyyyy_yzzzz[i] = 5.0 * ts_yyyy_yzzzz[i] * fe_0 + ts_yyyyy_zzzz[i] * fe_0 + ts_yyyyy_yzzzz[i] * pa_y[i];

        ts_yyyyyy_zzzzz[i] = 5.0 * ts_yyyy_zzzzz[i] * fe_0 + ts_yyyyy_zzzzz[i] * pa_y[i];
    }

    // Set up 462-483 components of targeted buffer : IH

    auto ts_yyyyyz_xxxxx = pbuffer.data(idx_ovl_ih + 462);

    auto ts_yyyyyz_xxxxy = pbuffer.data(idx_ovl_ih + 463);

    auto ts_yyyyyz_xxxxz = pbuffer.data(idx_ovl_ih + 464);

    auto ts_yyyyyz_xxxyy = pbuffer.data(idx_ovl_ih + 465);

    auto ts_yyyyyz_xxxyz = pbuffer.data(idx_ovl_ih + 466);

    auto ts_yyyyyz_xxxzz = pbuffer.data(idx_ovl_ih + 467);

    auto ts_yyyyyz_xxyyy = pbuffer.data(idx_ovl_ih + 468);

    auto ts_yyyyyz_xxyyz = pbuffer.data(idx_ovl_ih + 469);

    auto ts_yyyyyz_xxyzz = pbuffer.data(idx_ovl_ih + 470);

    auto ts_yyyyyz_xxzzz = pbuffer.data(idx_ovl_ih + 471);

    auto ts_yyyyyz_xyyyy = pbuffer.data(idx_ovl_ih + 472);

    auto ts_yyyyyz_xyyyz = pbuffer.data(idx_ovl_ih + 473);

    auto ts_yyyyyz_xyyzz = pbuffer.data(idx_ovl_ih + 474);

    auto ts_yyyyyz_xyzzz = pbuffer.data(idx_ovl_ih + 475);

    auto ts_yyyyyz_xzzzz = pbuffer.data(idx_ovl_ih + 476);

    auto ts_yyyyyz_yyyyy = pbuffer.data(idx_ovl_ih + 477);

    auto ts_yyyyyz_yyyyz = pbuffer.data(idx_ovl_ih + 478);

    auto ts_yyyyyz_yyyzz = pbuffer.data(idx_ovl_ih + 479);

    auto ts_yyyyyz_yyzzz = pbuffer.data(idx_ovl_ih + 480);

    auto ts_yyyyyz_yzzzz = pbuffer.data(idx_ovl_ih + 481);

    auto ts_yyyyyz_zzzzz = pbuffer.data(idx_ovl_ih + 482);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyyy_xxxxx, ts_yyyyy_xxxxy, ts_yyyyy_xxxy, ts_yyyyy_xxxyy, ts_yyyyy_xxxyz, ts_yyyyy_xxyy, ts_yyyyy_xxyyy, ts_yyyyy_xxyyz, ts_yyyyy_xxyz, ts_yyyyy_xxyzz, ts_yyyyy_xyyy, ts_yyyyy_xyyyy, ts_yyyyy_xyyyz, ts_yyyyy_xyyz, ts_yyyyy_xyyzz, ts_yyyyy_xyzz, ts_yyyyy_xyzzz, ts_yyyyy_yyyy, ts_yyyyy_yyyyy, ts_yyyyy_yyyyz, ts_yyyyy_yyyz, ts_yyyyy_yyyzz, ts_yyyyy_yyzz, ts_yyyyy_yyzzz, ts_yyyyy_yzzz, ts_yyyyy_yzzzz, ts_yyyyyz_xxxxx, ts_yyyyyz_xxxxy, ts_yyyyyz_xxxxz, ts_yyyyyz_xxxyy, ts_yyyyyz_xxxyz, ts_yyyyyz_xxxzz, ts_yyyyyz_xxyyy, ts_yyyyyz_xxyyz, ts_yyyyyz_xxyzz, ts_yyyyyz_xxzzz, ts_yyyyyz_xyyyy, ts_yyyyyz_xyyyz, ts_yyyyyz_xyyzz, ts_yyyyyz_xyzzz, ts_yyyyyz_xzzzz, ts_yyyyyz_yyyyy, ts_yyyyyz_yyyyz, ts_yyyyyz_yyyzz, ts_yyyyyz_yyzzz, ts_yyyyyz_yzzzz, ts_yyyyyz_zzzzz, ts_yyyyz_xxxxz, ts_yyyyz_xxxzz, ts_yyyyz_xxzzz, ts_yyyyz_xzzzz, ts_yyyyz_zzzzz, ts_yyyz_xxxxz, ts_yyyz_xxxzz, ts_yyyz_xxzzz, ts_yyyz_xzzzz, ts_yyyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyyz_xxxxx[i] = ts_yyyyy_xxxxx[i] * pa_z[i];

        ts_yyyyyz_xxxxy[i] = ts_yyyyy_xxxxy[i] * pa_z[i];

        ts_yyyyyz_xxxxz[i] = 4.0 * ts_yyyz_xxxxz[i] * fe_0 + ts_yyyyz_xxxxz[i] * pa_y[i];

        ts_yyyyyz_xxxyy[i] = ts_yyyyy_xxxyy[i] * pa_z[i];

        ts_yyyyyz_xxxyz[i] = ts_yyyyy_xxxy[i] * fe_0 + ts_yyyyy_xxxyz[i] * pa_z[i];

        ts_yyyyyz_xxxzz[i] = 4.0 * ts_yyyz_xxxzz[i] * fe_0 + ts_yyyyz_xxxzz[i] * pa_y[i];

        ts_yyyyyz_xxyyy[i] = ts_yyyyy_xxyyy[i] * pa_z[i];

        ts_yyyyyz_xxyyz[i] = ts_yyyyy_xxyy[i] * fe_0 + ts_yyyyy_xxyyz[i] * pa_z[i];

        ts_yyyyyz_xxyzz[i] = 2.0 * ts_yyyyy_xxyz[i] * fe_0 + ts_yyyyy_xxyzz[i] * pa_z[i];

        ts_yyyyyz_xxzzz[i] = 4.0 * ts_yyyz_xxzzz[i] * fe_0 + ts_yyyyz_xxzzz[i] * pa_y[i];

        ts_yyyyyz_xyyyy[i] = ts_yyyyy_xyyyy[i] * pa_z[i];

        ts_yyyyyz_xyyyz[i] = ts_yyyyy_xyyy[i] * fe_0 + ts_yyyyy_xyyyz[i] * pa_z[i];

        ts_yyyyyz_xyyzz[i] = 2.0 * ts_yyyyy_xyyz[i] * fe_0 + ts_yyyyy_xyyzz[i] * pa_z[i];

        ts_yyyyyz_xyzzz[i] = 3.0 * ts_yyyyy_xyzz[i] * fe_0 + ts_yyyyy_xyzzz[i] * pa_z[i];

        ts_yyyyyz_xzzzz[i] = 4.0 * ts_yyyz_xzzzz[i] * fe_0 + ts_yyyyz_xzzzz[i] * pa_y[i];

        ts_yyyyyz_yyyyy[i] = ts_yyyyy_yyyyy[i] * pa_z[i];

        ts_yyyyyz_yyyyz[i] = ts_yyyyy_yyyy[i] * fe_0 + ts_yyyyy_yyyyz[i] * pa_z[i];

        ts_yyyyyz_yyyzz[i] = 2.0 * ts_yyyyy_yyyz[i] * fe_0 + ts_yyyyy_yyyzz[i] * pa_z[i];

        ts_yyyyyz_yyzzz[i] = 3.0 * ts_yyyyy_yyzz[i] * fe_0 + ts_yyyyy_yyzzz[i] * pa_z[i];

        ts_yyyyyz_yzzzz[i] = 4.0 * ts_yyyyy_yzzz[i] * fe_0 + ts_yyyyy_yzzzz[i] * pa_z[i];

        ts_yyyyyz_zzzzz[i] = 4.0 * ts_yyyz_zzzzz[i] * fe_0 + ts_yyyyz_zzzzz[i] * pa_y[i];
    }

    // Set up 483-504 components of targeted buffer : IH

    auto ts_yyyyzz_xxxxx = pbuffer.data(idx_ovl_ih + 483);

    auto ts_yyyyzz_xxxxy = pbuffer.data(idx_ovl_ih + 484);

    auto ts_yyyyzz_xxxxz = pbuffer.data(idx_ovl_ih + 485);

    auto ts_yyyyzz_xxxyy = pbuffer.data(idx_ovl_ih + 486);

    auto ts_yyyyzz_xxxyz = pbuffer.data(idx_ovl_ih + 487);

    auto ts_yyyyzz_xxxzz = pbuffer.data(idx_ovl_ih + 488);

    auto ts_yyyyzz_xxyyy = pbuffer.data(idx_ovl_ih + 489);

    auto ts_yyyyzz_xxyyz = pbuffer.data(idx_ovl_ih + 490);

    auto ts_yyyyzz_xxyzz = pbuffer.data(idx_ovl_ih + 491);

    auto ts_yyyyzz_xxzzz = pbuffer.data(idx_ovl_ih + 492);

    auto ts_yyyyzz_xyyyy = pbuffer.data(idx_ovl_ih + 493);

    auto ts_yyyyzz_xyyyz = pbuffer.data(idx_ovl_ih + 494);

    auto ts_yyyyzz_xyyzz = pbuffer.data(idx_ovl_ih + 495);

    auto ts_yyyyzz_xyzzz = pbuffer.data(idx_ovl_ih + 496);

    auto ts_yyyyzz_xzzzz = pbuffer.data(idx_ovl_ih + 497);

    auto ts_yyyyzz_yyyyy = pbuffer.data(idx_ovl_ih + 498);

    auto ts_yyyyzz_yyyyz = pbuffer.data(idx_ovl_ih + 499);

    auto ts_yyyyzz_yyyzz = pbuffer.data(idx_ovl_ih + 500);

    auto ts_yyyyzz_yyzzz = pbuffer.data(idx_ovl_ih + 501);

    auto ts_yyyyzz_yzzzz = pbuffer.data(idx_ovl_ih + 502);

    auto ts_yyyyzz_zzzzz = pbuffer.data(idx_ovl_ih + 503);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyy_xxxxy, ts_yyyy_xxxyy, ts_yyyy_xxyyy, ts_yyyy_xyyyy, ts_yyyy_yyyyy, ts_yyyyz_xxxxy, ts_yyyyz_xxxyy, ts_yyyyz_xxyyy, ts_yyyyz_xyyyy, ts_yyyyz_yyyyy, ts_yyyyzz_xxxxx, ts_yyyyzz_xxxxy, ts_yyyyzz_xxxxz, ts_yyyyzz_xxxyy, ts_yyyyzz_xxxyz, ts_yyyyzz_xxxzz, ts_yyyyzz_xxyyy, ts_yyyyzz_xxyyz, ts_yyyyzz_xxyzz, ts_yyyyzz_xxzzz, ts_yyyyzz_xyyyy, ts_yyyyzz_xyyyz, ts_yyyyzz_xyyzz, ts_yyyyzz_xyzzz, ts_yyyyzz_xzzzz, ts_yyyyzz_yyyyy, ts_yyyyzz_yyyyz, ts_yyyyzz_yyyzz, ts_yyyyzz_yyzzz, ts_yyyyzz_yzzzz, ts_yyyyzz_zzzzz, ts_yyyzz_xxxxx, ts_yyyzz_xxxxz, ts_yyyzz_xxxyz, ts_yyyzz_xxxz, ts_yyyzz_xxxzz, ts_yyyzz_xxyyz, ts_yyyzz_xxyz, ts_yyyzz_xxyzz, ts_yyyzz_xxzz, ts_yyyzz_xxzzz, ts_yyyzz_xyyyz, ts_yyyzz_xyyz, ts_yyyzz_xyyzz, ts_yyyzz_xyzz, ts_yyyzz_xyzzz, ts_yyyzz_xzzz, ts_yyyzz_xzzzz, ts_yyyzz_yyyyz, ts_yyyzz_yyyz, ts_yyyzz_yyyzz, ts_yyyzz_yyzz, ts_yyyzz_yyzzz, ts_yyyzz_yzzz, ts_yyyzz_yzzzz, ts_yyyzz_zzzz, ts_yyyzz_zzzzz, ts_yyzz_xxxxx, ts_yyzz_xxxxz, ts_yyzz_xxxyz, ts_yyzz_xxxzz, ts_yyzz_xxyyz, ts_yyzz_xxyzz, ts_yyzz_xxzzz, ts_yyzz_xyyyz, ts_yyzz_xyyzz, ts_yyzz_xyzzz, ts_yyzz_xzzzz, ts_yyzz_yyyyz, ts_yyzz_yyyzz, ts_yyzz_yyzzz, ts_yyzz_yzzzz, ts_yyzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyzz_xxxxx[i] = 3.0 * ts_yyzz_xxxxx[i] * fe_0 + ts_yyyzz_xxxxx[i] * pa_y[i];

        ts_yyyyzz_xxxxy[i] = ts_yyyy_xxxxy[i] * fe_0 + ts_yyyyz_xxxxy[i] * pa_z[i];

        ts_yyyyzz_xxxxz[i] = 3.0 * ts_yyzz_xxxxz[i] * fe_0 + ts_yyyzz_xxxxz[i] * pa_y[i];

        ts_yyyyzz_xxxyy[i] = ts_yyyy_xxxyy[i] * fe_0 + ts_yyyyz_xxxyy[i] * pa_z[i];

        ts_yyyyzz_xxxyz[i] = 3.0 * ts_yyzz_xxxyz[i] * fe_0 + ts_yyyzz_xxxz[i] * fe_0 + ts_yyyzz_xxxyz[i] * pa_y[i];

        ts_yyyyzz_xxxzz[i] = 3.0 * ts_yyzz_xxxzz[i] * fe_0 + ts_yyyzz_xxxzz[i] * pa_y[i];

        ts_yyyyzz_xxyyy[i] = ts_yyyy_xxyyy[i] * fe_0 + ts_yyyyz_xxyyy[i] * pa_z[i];

        ts_yyyyzz_xxyyz[i] = 3.0 * ts_yyzz_xxyyz[i] * fe_0 + 2.0 * ts_yyyzz_xxyz[i] * fe_0 + ts_yyyzz_xxyyz[i] * pa_y[i];

        ts_yyyyzz_xxyzz[i] = 3.0 * ts_yyzz_xxyzz[i] * fe_0 + ts_yyyzz_xxzz[i] * fe_0 + ts_yyyzz_xxyzz[i] * pa_y[i];

        ts_yyyyzz_xxzzz[i] = 3.0 * ts_yyzz_xxzzz[i] * fe_0 + ts_yyyzz_xxzzz[i] * pa_y[i];

        ts_yyyyzz_xyyyy[i] = ts_yyyy_xyyyy[i] * fe_0 + ts_yyyyz_xyyyy[i] * pa_z[i];

        ts_yyyyzz_xyyyz[i] = 3.0 * ts_yyzz_xyyyz[i] * fe_0 + 3.0 * ts_yyyzz_xyyz[i] * fe_0 + ts_yyyzz_xyyyz[i] * pa_y[i];

        ts_yyyyzz_xyyzz[i] = 3.0 * ts_yyzz_xyyzz[i] * fe_0 + 2.0 * ts_yyyzz_xyzz[i] * fe_0 + ts_yyyzz_xyyzz[i] * pa_y[i];

        ts_yyyyzz_xyzzz[i] = 3.0 * ts_yyzz_xyzzz[i] * fe_0 + ts_yyyzz_xzzz[i] * fe_0 + ts_yyyzz_xyzzz[i] * pa_y[i];

        ts_yyyyzz_xzzzz[i] = 3.0 * ts_yyzz_xzzzz[i] * fe_0 + ts_yyyzz_xzzzz[i] * pa_y[i];

        ts_yyyyzz_yyyyy[i] = ts_yyyy_yyyyy[i] * fe_0 + ts_yyyyz_yyyyy[i] * pa_z[i];

        ts_yyyyzz_yyyyz[i] = 3.0 * ts_yyzz_yyyyz[i] * fe_0 + 4.0 * ts_yyyzz_yyyz[i] * fe_0 + ts_yyyzz_yyyyz[i] * pa_y[i];

        ts_yyyyzz_yyyzz[i] = 3.0 * ts_yyzz_yyyzz[i] * fe_0 + 3.0 * ts_yyyzz_yyzz[i] * fe_0 + ts_yyyzz_yyyzz[i] * pa_y[i];

        ts_yyyyzz_yyzzz[i] = 3.0 * ts_yyzz_yyzzz[i] * fe_0 + 2.0 * ts_yyyzz_yzzz[i] * fe_0 + ts_yyyzz_yyzzz[i] * pa_y[i];

        ts_yyyyzz_yzzzz[i] = 3.0 * ts_yyzz_yzzzz[i] * fe_0 + ts_yyyzz_zzzz[i] * fe_0 + ts_yyyzz_yzzzz[i] * pa_y[i];

        ts_yyyyzz_zzzzz[i] = 3.0 * ts_yyzz_zzzzz[i] * fe_0 + ts_yyyzz_zzzzz[i] * pa_y[i];
    }

    // Set up 504-525 components of targeted buffer : IH

    auto ts_yyyzzz_xxxxx = pbuffer.data(idx_ovl_ih + 504);

    auto ts_yyyzzz_xxxxy = pbuffer.data(idx_ovl_ih + 505);

    auto ts_yyyzzz_xxxxz = pbuffer.data(idx_ovl_ih + 506);

    auto ts_yyyzzz_xxxyy = pbuffer.data(idx_ovl_ih + 507);

    auto ts_yyyzzz_xxxyz = pbuffer.data(idx_ovl_ih + 508);

    auto ts_yyyzzz_xxxzz = pbuffer.data(idx_ovl_ih + 509);

    auto ts_yyyzzz_xxyyy = pbuffer.data(idx_ovl_ih + 510);

    auto ts_yyyzzz_xxyyz = pbuffer.data(idx_ovl_ih + 511);

    auto ts_yyyzzz_xxyzz = pbuffer.data(idx_ovl_ih + 512);

    auto ts_yyyzzz_xxzzz = pbuffer.data(idx_ovl_ih + 513);

    auto ts_yyyzzz_xyyyy = pbuffer.data(idx_ovl_ih + 514);

    auto ts_yyyzzz_xyyyz = pbuffer.data(idx_ovl_ih + 515);

    auto ts_yyyzzz_xyyzz = pbuffer.data(idx_ovl_ih + 516);

    auto ts_yyyzzz_xyzzz = pbuffer.data(idx_ovl_ih + 517);

    auto ts_yyyzzz_xzzzz = pbuffer.data(idx_ovl_ih + 518);

    auto ts_yyyzzz_yyyyy = pbuffer.data(idx_ovl_ih + 519);

    auto ts_yyyzzz_yyyyz = pbuffer.data(idx_ovl_ih + 520);

    auto ts_yyyzzz_yyyzz = pbuffer.data(idx_ovl_ih + 521);

    auto ts_yyyzzz_yyzzz = pbuffer.data(idx_ovl_ih + 522);

    auto ts_yyyzzz_yzzzz = pbuffer.data(idx_ovl_ih + 523);

    auto ts_yyyzzz_zzzzz = pbuffer.data(idx_ovl_ih + 524);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyz_xxxxy, ts_yyyz_xxxyy, ts_yyyz_xxyyy, ts_yyyz_xyyyy, ts_yyyz_yyyyy, ts_yyyzz_xxxxy, ts_yyyzz_xxxyy, ts_yyyzz_xxyyy, ts_yyyzz_xyyyy, ts_yyyzz_yyyyy, ts_yyyzzz_xxxxx, ts_yyyzzz_xxxxy, ts_yyyzzz_xxxxz, ts_yyyzzz_xxxyy, ts_yyyzzz_xxxyz, ts_yyyzzz_xxxzz, ts_yyyzzz_xxyyy, ts_yyyzzz_xxyyz, ts_yyyzzz_xxyzz, ts_yyyzzz_xxzzz, ts_yyyzzz_xyyyy, ts_yyyzzz_xyyyz, ts_yyyzzz_xyyzz, ts_yyyzzz_xyzzz, ts_yyyzzz_xzzzz, ts_yyyzzz_yyyyy, ts_yyyzzz_yyyyz, ts_yyyzzz_yyyzz, ts_yyyzzz_yyzzz, ts_yyyzzz_yzzzz, ts_yyyzzz_zzzzz, ts_yyzzz_xxxxx, ts_yyzzz_xxxxz, ts_yyzzz_xxxyz, ts_yyzzz_xxxz, ts_yyzzz_xxxzz, ts_yyzzz_xxyyz, ts_yyzzz_xxyz, ts_yyzzz_xxyzz, ts_yyzzz_xxzz, ts_yyzzz_xxzzz, ts_yyzzz_xyyyz, ts_yyzzz_xyyz, ts_yyzzz_xyyzz, ts_yyzzz_xyzz, ts_yyzzz_xyzzz, ts_yyzzz_xzzz, ts_yyzzz_xzzzz, ts_yyzzz_yyyyz, ts_yyzzz_yyyz, ts_yyzzz_yyyzz, ts_yyzzz_yyzz, ts_yyzzz_yyzzz, ts_yyzzz_yzzz, ts_yyzzz_yzzzz, ts_yyzzz_zzzz, ts_yyzzz_zzzzz, ts_yzzz_xxxxx, ts_yzzz_xxxxz, ts_yzzz_xxxyz, ts_yzzz_xxxzz, ts_yzzz_xxyyz, ts_yzzz_xxyzz, ts_yzzz_xxzzz, ts_yzzz_xyyyz, ts_yzzz_xyyzz, ts_yzzz_xyzzz, ts_yzzz_xzzzz, ts_yzzz_yyyyz, ts_yzzz_yyyzz, ts_yzzz_yyzzz, ts_yzzz_yzzzz, ts_yzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyzzz_xxxxx[i] = 2.0 * ts_yzzz_xxxxx[i] * fe_0 + ts_yyzzz_xxxxx[i] * pa_y[i];

        ts_yyyzzz_xxxxy[i] = 2.0 * ts_yyyz_xxxxy[i] * fe_0 + ts_yyyzz_xxxxy[i] * pa_z[i];

        ts_yyyzzz_xxxxz[i] = 2.0 * ts_yzzz_xxxxz[i] * fe_0 + ts_yyzzz_xxxxz[i] * pa_y[i];

        ts_yyyzzz_xxxyy[i] = 2.0 * ts_yyyz_xxxyy[i] * fe_0 + ts_yyyzz_xxxyy[i] * pa_z[i];

        ts_yyyzzz_xxxyz[i] = 2.0 * ts_yzzz_xxxyz[i] * fe_0 + ts_yyzzz_xxxz[i] * fe_0 + ts_yyzzz_xxxyz[i] * pa_y[i];

        ts_yyyzzz_xxxzz[i] = 2.0 * ts_yzzz_xxxzz[i] * fe_0 + ts_yyzzz_xxxzz[i] * pa_y[i];

        ts_yyyzzz_xxyyy[i] = 2.0 * ts_yyyz_xxyyy[i] * fe_0 + ts_yyyzz_xxyyy[i] * pa_z[i];

        ts_yyyzzz_xxyyz[i] = 2.0 * ts_yzzz_xxyyz[i] * fe_0 + 2.0 * ts_yyzzz_xxyz[i] * fe_0 + ts_yyzzz_xxyyz[i] * pa_y[i];

        ts_yyyzzz_xxyzz[i] = 2.0 * ts_yzzz_xxyzz[i] * fe_0 + ts_yyzzz_xxzz[i] * fe_0 + ts_yyzzz_xxyzz[i] * pa_y[i];

        ts_yyyzzz_xxzzz[i] = 2.0 * ts_yzzz_xxzzz[i] * fe_0 + ts_yyzzz_xxzzz[i] * pa_y[i];

        ts_yyyzzz_xyyyy[i] = 2.0 * ts_yyyz_xyyyy[i] * fe_0 + ts_yyyzz_xyyyy[i] * pa_z[i];

        ts_yyyzzz_xyyyz[i] = 2.0 * ts_yzzz_xyyyz[i] * fe_0 + 3.0 * ts_yyzzz_xyyz[i] * fe_0 + ts_yyzzz_xyyyz[i] * pa_y[i];

        ts_yyyzzz_xyyzz[i] = 2.0 * ts_yzzz_xyyzz[i] * fe_0 + 2.0 * ts_yyzzz_xyzz[i] * fe_0 + ts_yyzzz_xyyzz[i] * pa_y[i];

        ts_yyyzzz_xyzzz[i] = 2.0 * ts_yzzz_xyzzz[i] * fe_0 + ts_yyzzz_xzzz[i] * fe_0 + ts_yyzzz_xyzzz[i] * pa_y[i];

        ts_yyyzzz_xzzzz[i] = 2.0 * ts_yzzz_xzzzz[i] * fe_0 + ts_yyzzz_xzzzz[i] * pa_y[i];

        ts_yyyzzz_yyyyy[i] = 2.0 * ts_yyyz_yyyyy[i] * fe_0 + ts_yyyzz_yyyyy[i] * pa_z[i];

        ts_yyyzzz_yyyyz[i] = 2.0 * ts_yzzz_yyyyz[i] * fe_0 + 4.0 * ts_yyzzz_yyyz[i] * fe_0 + ts_yyzzz_yyyyz[i] * pa_y[i];

        ts_yyyzzz_yyyzz[i] = 2.0 * ts_yzzz_yyyzz[i] * fe_0 + 3.0 * ts_yyzzz_yyzz[i] * fe_0 + ts_yyzzz_yyyzz[i] * pa_y[i];

        ts_yyyzzz_yyzzz[i] = 2.0 * ts_yzzz_yyzzz[i] * fe_0 + 2.0 * ts_yyzzz_yzzz[i] * fe_0 + ts_yyzzz_yyzzz[i] * pa_y[i];

        ts_yyyzzz_yzzzz[i] = 2.0 * ts_yzzz_yzzzz[i] * fe_0 + ts_yyzzz_zzzz[i] * fe_0 + ts_yyzzz_yzzzz[i] * pa_y[i];

        ts_yyyzzz_zzzzz[i] = 2.0 * ts_yzzz_zzzzz[i] * fe_0 + ts_yyzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 525-546 components of targeted buffer : IH

    auto ts_yyzzzz_xxxxx = pbuffer.data(idx_ovl_ih + 525);

    auto ts_yyzzzz_xxxxy = pbuffer.data(idx_ovl_ih + 526);

    auto ts_yyzzzz_xxxxz = pbuffer.data(idx_ovl_ih + 527);

    auto ts_yyzzzz_xxxyy = pbuffer.data(idx_ovl_ih + 528);

    auto ts_yyzzzz_xxxyz = pbuffer.data(idx_ovl_ih + 529);

    auto ts_yyzzzz_xxxzz = pbuffer.data(idx_ovl_ih + 530);

    auto ts_yyzzzz_xxyyy = pbuffer.data(idx_ovl_ih + 531);

    auto ts_yyzzzz_xxyyz = pbuffer.data(idx_ovl_ih + 532);

    auto ts_yyzzzz_xxyzz = pbuffer.data(idx_ovl_ih + 533);

    auto ts_yyzzzz_xxzzz = pbuffer.data(idx_ovl_ih + 534);

    auto ts_yyzzzz_xyyyy = pbuffer.data(idx_ovl_ih + 535);

    auto ts_yyzzzz_xyyyz = pbuffer.data(idx_ovl_ih + 536);

    auto ts_yyzzzz_xyyzz = pbuffer.data(idx_ovl_ih + 537);

    auto ts_yyzzzz_xyzzz = pbuffer.data(idx_ovl_ih + 538);

    auto ts_yyzzzz_xzzzz = pbuffer.data(idx_ovl_ih + 539);

    auto ts_yyzzzz_yyyyy = pbuffer.data(idx_ovl_ih + 540);

    auto ts_yyzzzz_yyyyz = pbuffer.data(idx_ovl_ih + 541);

    auto ts_yyzzzz_yyyzz = pbuffer.data(idx_ovl_ih + 542);

    auto ts_yyzzzz_yyzzz = pbuffer.data(idx_ovl_ih + 543);

    auto ts_yyzzzz_yzzzz = pbuffer.data(idx_ovl_ih + 544);

    auto ts_yyzzzz_zzzzz = pbuffer.data(idx_ovl_ih + 545);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyzz_xxxxy, ts_yyzz_xxxyy, ts_yyzz_xxyyy, ts_yyzz_xyyyy, ts_yyzz_yyyyy, ts_yyzzz_xxxxy, ts_yyzzz_xxxyy, ts_yyzzz_xxyyy, ts_yyzzz_xyyyy, ts_yyzzz_yyyyy, ts_yyzzzz_xxxxx, ts_yyzzzz_xxxxy, ts_yyzzzz_xxxxz, ts_yyzzzz_xxxyy, ts_yyzzzz_xxxyz, ts_yyzzzz_xxxzz, ts_yyzzzz_xxyyy, ts_yyzzzz_xxyyz, ts_yyzzzz_xxyzz, ts_yyzzzz_xxzzz, ts_yyzzzz_xyyyy, ts_yyzzzz_xyyyz, ts_yyzzzz_xyyzz, ts_yyzzzz_xyzzz, ts_yyzzzz_xzzzz, ts_yyzzzz_yyyyy, ts_yyzzzz_yyyyz, ts_yyzzzz_yyyzz, ts_yyzzzz_yyzzz, ts_yyzzzz_yzzzz, ts_yyzzzz_zzzzz, ts_yzzzz_xxxxx, ts_yzzzz_xxxxz, ts_yzzzz_xxxyz, ts_yzzzz_xxxz, ts_yzzzz_xxxzz, ts_yzzzz_xxyyz, ts_yzzzz_xxyz, ts_yzzzz_xxyzz, ts_yzzzz_xxzz, ts_yzzzz_xxzzz, ts_yzzzz_xyyyz, ts_yzzzz_xyyz, ts_yzzzz_xyyzz, ts_yzzzz_xyzz, ts_yzzzz_xyzzz, ts_yzzzz_xzzz, ts_yzzzz_xzzzz, ts_yzzzz_yyyyz, ts_yzzzz_yyyz, ts_yzzzz_yyyzz, ts_yzzzz_yyzz, ts_yzzzz_yyzzz, ts_yzzzz_yzzz, ts_yzzzz_yzzzz, ts_yzzzz_zzzz, ts_yzzzz_zzzzz, ts_zzzz_xxxxx, ts_zzzz_xxxxz, ts_zzzz_xxxyz, ts_zzzz_xxxzz, ts_zzzz_xxyyz, ts_zzzz_xxyzz, ts_zzzz_xxzzz, ts_zzzz_xyyyz, ts_zzzz_xyyzz, ts_zzzz_xyzzz, ts_zzzz_xzzzz, ts_zzzz_yyyyz, ts_zzzz_yyyzz, ts_zzzz_yyzzz, ts_zzzz_yzzzz, ts_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzzzz_xxxxx[i] = ts_zzzz_xxxxx[i] * fe_0 + ts_yzzzz_xxxxx[i] * pa_y[i];

        ts_yyzzzz_xxxxy[i] = 3.0 * ts_yyzz_xxxxy[i] * fe_0 + ts_yyzzz_xxxxy[i] * pa_z[i];

        ts_yyzzzz_xxxxz[i] = ts_zzzz_xxxxz[i] * fe_0 + ts_yzzzz_xxxxz[i] * pa_y[i];

        ts_yyzzzz_xxxyy[i] = 3.0 * ts_yyzz_xxxyy[i] * fe_0 + ts_yyzzz_xxxyy[i] * pa_z[i];

        ts_yyzzzz_xxxyz[i] = ts_zzzz_xxxyz[i] * fe_0 + ts_yzzzz_xxxz[i] * fe_0 + ts_yzzzz_xxxyz[i] * pa_y[i];

        ts_yyzzzz_xxxzz[i] = ts_zzzz_xxxzz[i] * fe_0 + ts_yzzzz_xxxzz[i] * pa_y[i];

        ts_yyzzzz_xxyyy[i] = 3.0 * ts_yyzz_xxyyy[i] * fe_0 + ts_yyzzz_xxyyy[i] * pa_z[i];

        ts_yyzzzz_xxyyz[i] = ts_zzzz_xxyyz[i] * fe_0 + 2.0 * ts_yzzzz_xxyz[i] * fe_0 + ts_yzzzz_xxyyz[i] * pa_y[i];

        ts_yyzzzz_xxyzz[i] = ts_zzzz_xxyzz[i] * fe_0 + ts_yzzzz_xxzz[i] * fe_0 + ts_yzzzz_xxyzz[i] * pa_y[i];

        ts_yyzzzz_xxzzz[i] = ts_zzzz_xxzzz[i] * fe_0 + ts_yzzzz_xxzzz[i] * pa_y[i];

        ts_yyzzzz_xyyyy[i] = 3.0 * ts_yyzz_xyyyy[i] * fe_0 + ts_yyzzz_xyyyy[i] * pa_z[i];

        ts_yyzzzz_xyyyz[i] = ts_zzzz_xyyyz[i] * fe_0 + 3.0 * ts_yzzzz_xyyz[i] * fe_0 + ts_yzzzz_xyyyz[i] * pa_y[i];

        ts_yyzzzz_xyyzz[i] = ts_zzzz_xyyzz[i] * fe_0 + 2.0 * ts_yzzzz_xyzz[i] * fe_0 + ts_yzzzz_xyyzz[i] * pa_y[i];

        ts_yyzzzz_xyzzz[i] = ts_zzzz_xyzzz[i] * fe_0 + ts_yzzzz_xzzz[i] * fe_0 + ts_yzzzz_xyzzz[i] * pa_y[i];

        ts_yyzzzz_xzzzz[i] = ts_zzzz_xzzzz[i] * fe_0 + ts_yzzzz_xzzzz[i] * pa_y[i];

        ts_yyzzzz_yyyyy[i] = 3.0 * ts_yyzz_yyyyy[i] * fe_0 + ts_yyzzz_yyyyy[i] * pa_z[i];

        ts_yyzzzz_yyyyz[i] = ts_zzzz_yyyyz[i] * fe_0 + 4.0 * ts_yzzzz_yyyz[i] * fe_0 + ts_yzzzz_yyyyz[i] * pa_y[i];

        ts_yyzzzz_yyyzz[i] = ts_zzzz_yyyzz[i] * fe_0 + 3.0 * ts_yzzzz_yyzz[i] * fe_0 + ts_yzzzz_yyyzz[i] * pa_y[i];

        ts_yyzzzz_yyzzz[i] = ts_zzzz_yyzzz[i] * fe_0 + 2.0 * ts_yzzzz_yzzz[i] * fe_0 + ts_yzzzz_yyzzz[i] * pa_y[i];

        ts_yyzzzz_yzzzz[i] = ts_zzzz_yzzzz[i] * fe_0 + ts_yzzzz_zzzz[i] * fe_0 + ts_yzzzz_yzzzz[i] * pa_y[i];

        ts_yyzzzz_zzzzz[i] = ts_zzzz_zzzzz[i] * fe_0 + ts_yzzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 546-567 components of targeted buffer : IH

    auto ts_yzzzzz_xxxxx = pbuffer.data(idx_ovl_ih + 546);

    auto ts_yzzzzz_xxxxy = pbuffer.data(idx_ovl_ih + 547);

    auto ts_yzzzzz_xxxxz = pbuffer.data(idx_ovl_ih + 548);

    auto ts_yzzzzz_xxxyy = pbuffer.data(idx_ovl_ih + 549);

    auto ts_yzzzzz_xxxyz = pbuffer.data(idx_ovl_ih + 550);

    auto ts_yzzzzz_xxxzz = pbuffer.data(idx_ovl_ih + 551);

    auto ts_yzzzzz_xxyyy = pbuffer.data(idx_ovl_ih + 552);

    auto ts_yzzzzz_xxyyz = pbuffer.data(idx_ovl_ih + 553);

    auto ts_yzzzzz_xxyzz = pbuffer.data(idx_ovl_ih + 554);

    auto ts_yzzzzz_xxzzz = pbuffer.data(idx_ovl_ih + 555);

    auto ts_yzzzzz_xyyyy = pbuffer.data(idx_ovl_ih + 556);

    auto ts_yzzzzz_xyyyz = pbuffer.data(idx_ovl_ih + 557);

    auto ts_yzzzzz_xyyzz = pbuffer.data(idx_ovl_ih + 558);

    auto ts_yzzzzz_xyzzz = pbuffer.data(idx_ovl_ih + 559);

    auto ts_yzzzzz_xzzzz = pbuffer.data(idx_ovl_ih + 560);

    auto ts_yzzzzz_yyyyy = pbuffer.data(idx_ovl_ih + 561);

    auto ts_yzzzzz_yyyyz = pbuffer.data(idx_ovl_ih + 562);

    auto ts_yzzzzz_yyyzz = pbuffer.data(idx_ovl_ih + 563);

    auto ts_yzzzzz_yyzzz = pbuffer.data(idx_ovl_ih + 564);

    auto ts_yzzzzz_yzzzz = pbuffer.data(idx_ovl_ih + 565);

    auto ts_yzzzzz_zzzzz = pbuffer.data(idx_ovl_ih + 566);

    #pragma omp simd aligned(pa_y, ts_yzzzzz_xxxxx, ts_yzzzzz_xxxxy, ts_yzzzzz_xxxxz, ts_yzzzzz_xxxyy, ts_yzzzzz_xxxyz, ts_yzzzzz_xxxzz, ts_yzzzzz_xxyyy, ts_yzzzzz_xxyyz, ts_yzzzzz_xxyzz, ts_yzzzzz_xxzzz, ts_yzzzzz_xyyyy, ts_yzzzzz_xyyyz, ts_yzzzzz_xyyzz, ts_yzzzzz_xyzzz, ts_yzzzzz_xzzzz, ts_yzzzzz_yyyyy, ts_yzzzzz_yyyyz, ts_yzzzzz_yyyzz, ts_yzzzzz_yyzzz, ts_yzzzzz_yzzzz, ts_yzzzzz_zzzzz, ts_zzzzz_xxxx, ts_zzzzz_xxxxx, ts_zzzzz_xxxxy, ts_zzzzz_xxxxz, ts_zzzzz_xxxy, ts_zzzzz_xxxyy, ts_zzzzz_xxxyz, ts_zzzzz_xxxz, ts_zzzzz_xxxzz, ts_zzzzz_xxyy, ts_zzzzz_xxyyy, ts_zzzzz_xxyyz, ts_zzzzz_xxyz, ts_zzzzz_xxyzz, ts_zzzzz_xxzz, ts_zzzzz_xxzzz, ts_zzzzz_xyyy, ts_zzzzz_xyyyy, ts_zzzzz_xyyyz, ts_zzzzz_xyyz, ts_zzzzz_xyyzz, ts_zzzzz_xyzz, ts_zzzzz_xyzzz, ts_zzzzz_xzzz, ts_zzzzz_xzzzz, ts_zzzzz_yyyy, ts_zzzzz_yyyyy, ts_zzzzz_yyyyz, ts_zzzzz_yyyz, ts_zzzzz_yyyzz, ts_zzzzz_yyzz, ts_zzzzz_yyzzz, ts_zzzzz_yzzz, ts_zzzzz_yzzzz, ts_zzzzz_zzzz, ts_zzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzzzz_xxxxx[i] = ts_zzzzz_xxxxx[i] * pa_y[i];

        ts_yzzzzz_xxxxy[i] = ts_zzzzz_xxxx[i] * fe_0 + ts_zzzzz_xxxxy[i] * pa_y[i];

        ts_yzzzzz_xxxxz[i] = ts_zzzzz_xxxxz[i] * pa_y[i];

        ts_yzzzzz_xxxyy[i] = 2.0 * ts_zzzzz_xxxy[i] * fe_0 + ts_zzzzz_xxxyy[i] * pa_y[i];

        ts_yzzzzz_xxxyz[i] = ts_zzzzz_xxxz[i] * fe_0 + ts_zzzzz_xxxyz[i] * pa_y[i];

        ts_yzzzzz_xxxzz[i] = ts_zzzzz_xxxzz[i] * pa_y[i];

        ts_yzzzzz_xxyyy[i] = 3.0 * ts_zzzzz_xxyy[i] * fe_0 + ts_zzzzz_xxyyy[i] * pa_y[i];

        ts_yzzzzz_xxyyz[i] = 2.0 * ts_zzzzz_xxyz[i] * fe_0 + ts_zzzzz_xxyyz[i] * pa_y[i];

        ts_yzzzzz_xxyzz[i] = ts_zzzzz_xxzz[i] * fe_0 + ts_zzzzz_xxyzz[i] * pa_y[i];

        ts_yzzzzz_xxzzz[i] = ts_zzzzz_xxzzz[i] * pa_y[i];

        ts_yzzzzz_xyyyy[i] = 4.0 * ts_zzzzz_xyyy[i] * fe_0 + ts_zzzzz_xyyyy[i] * pa_y[i];

        ts_yzzzzz_xyyyz[i] = 3.0 * ts_zzzzz_xyyz[i] * fe_0 + ts_zzzzz_xyyyz[i] * pa_y[i];

        ts_yzzzzz_xyyzz[i] = 2.0 * ts_zzzzz_xyzz[i] * fe_0 + ts_zzzzz_xyyzz[i] * pa_y[i];

        ts_yzzzzz_xyzzz[i] = ts_zzzzz_xzzz[i] * fe_0 + ts_zzzzz_xyzzz[i] * pa_y[i];

        ts_yzzzzz_xzzzz[i] = ts_zzzzz_xzzzz[i] * pa_y[i];

        ts_yzzzzz_yyyyy[i] = 5.0 * ts_zzzzz_yyyy[i] * fe_0 + ts_zzzzz_yyyyy[i] * pa_y[i];

        ts_yzzzzz_yyyyz[i] = 4.0 * ts_zzzzz_yyyz[i] * fe_0 + ts_zzzzz_yyyyz[i] * pa_y[i];

        ts_yzzzzz_yyyzz[i] = 3.0 * ts_zzzzz_yyzz[i] * fe_0 + ts_zzzzz_yyyzz[i] * pa_y[i];

        ts_yzzzzz_yyzzz[i] = 2.0 * ts_zzzzz_yzzz[i] * fe_0 + ts_zzzzz_yyzzz[i] * pa_y[i];

        ts_yzzzzz_yzzzz[i] = ts_zzzzz_zzzz[i] * fe_0 + ts_zzzzz_yzzzz[i] * pa_y[i];

        ts_yzzzzz_zzzzz[i] = ts_zzzzz_zzzzz[i] * pa_y[i];
    }

    // Set up 567-588 components of targeted buffer : IH

    auto ts_zzzzzz_xxxxx = pbuffer.data(idx_ovl_ih + 567);

    auto ts_zzzzzz_xxxxy = pbuffer.data(idx_ovl_ih + 568);

    auto ts_zzzzzz_xxxxz = pbuffer.data(idx_ovl_ih + 569);

    auto ts_zzzzzz_xxxyy = pbuffer.data(idx_ovl_ih + 570);

    auto ts_zzzzzz_xxxyz = pbuffer.data(idx_ovl_ih + 571);

    auto ts_zzzzzz_xxxzz = pbuffer.data(idx_ovl_ih + 572);

    auto ts_zzzzzz_xxyyy = pbuffer.data(idx_ovl_ih + 573);

    auto ts_zzzzzz_xxyyz = pbuffer.data(idx_ovl_ih + 574);

    auto ts_zzzzzz_xxyzz = pbuffer.data(idx_ovl_ih + 575);

    auto ts_zzzzzz_xxzzz = pbuffer.data(idx_ovl_ih + 576);

    auto ts_zzzzzz_xyyyy = pbuffer.data(idx_ovl_ih + 577);

    auto ts_zzzzzz_xyyyz = pbuffer.data(idx_ovl_ih + 578);

    auto ts_zzzzzz_xyyzz = pbuffer.data(idx_ovl_ih + 579);

    auto ts_zzzzzz_xyzzz = pbuffer.data(idx_ovl_ih + 580);

    auto ts_zzzzzz_xzzzz = pbuffer.data(idx_ovl_ih + 581);

    auto ts_zzzzzz_yyyyy = pbuffer.data(idx_ovl_ih + 582);

    auto ts_zzzzzz_yyyyz = pbuffer.data(idx_ovl_ih + 583);

    auto ts_zzzzzz_yyyzz = pbuffer.data(idx_ovl_ih + 584);

    auto ts_zzzzzz_yyzzz = pbuffer.data(idx_ovl_ih + 585);

    auto ts_zzzzzz_yzzzz = pbuffer.data(idx_ovl_ih + 586);

    auto ts_zzzzzz_zzzzz = pbuffer.data(idx_ovl_ih + 587);

    #pragma omp simd aligned(pa_z, ts_zzzz_xxxxx, ts_zzzz_xxxxy, ts_zzzz_xxxxz, ts_zzzz_xxxyy, ts_zzzz_xxxyz, ts_zzzz_xxxzz, ts_zzzz_xxyyy, ts_zzzz_xxyyz, ts_zzzz_xxyzz, ts_zzzz_xxzzz, ts_zzzz_xyyyy, ts_zzzz_xyyyz, ts_zzzz_xyyzz, ts_zzzz_xyzzz, ts_zzzz_xzzzz, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyzz, ts_zzzz_yyzzz, ts_zzzz_yzzzz, ts_zzzz_zzzzz, ts_zzzzz_xxxx, ts_zzzzz_xxxxx, ts_zzzzz_xxxxy, ts_zzzzz_xxxxz, ts_zzzzz_xxxy, ts_zzzzz_xxxyy, ts_zzzzz_xxxyz, ts_zzzzz_xxxz, ts_zzzzz_xxxzz, ts_zzzzz_xxyy, ts_zzzzz_xxyyy, ts_zzzzz_xxyyz, ts_zzzzz_xxyz, ts_zzzzz_xxyzz, ts_zzzzz_xxzz, ts_zzzzz_xxzzz, ts_zzzzz_xyyy, ts_zzzzz_xyyyy, ts_zzzzz_xyyyz, ts_zzzzz_xyyz, ts_zzzzz_xyyzz, ts_zzzzz_xyzz, ts_zzzzz_xyzzz, ts_zzzzz_xzzz, ts_zzzzz_xzzzz, ts_zzzzz_yyyy, ts_zzzzz_yyyyy, ts_zzzzz_yyyyz, ts_zzzzz_yyyz, ts_zzzzz_yyyzz, ts_zzzzz_yyzz, ts_zzzzz_yyzzz, ts_zzzzz_yzzz, ts_zzzzz_yzzzz, ts_zzzzz_zzzz, ts_zzzzz_zzzzz, ts_zzzzzz_xxxxx, ts_zzzzzz_xxxxy, ts_zzzzzz_xxxxz, ts_zzzzzz_xxxyy, ts_zzzzzz_xxxyz, ts_zzzzzz_xxxzz, ts_zzzzzz_xxyyy, ts_zzzzzz_xxyyz, ts_zzzzzz_xxyzz, ts_zzzzzz_xxzzz, ts_zzzzzz_xyyyy, ts_zzzzzz_xyyyz, ts_zzzzzz_xyyzz, ts_zzzzzz_xyzzz, ts_zzzzzz_xzzzz, ts_zzzzzz_yyyyy, ts_zzzzzz_yyyyz, ts_zzzzzz_yyyzz, ts_zzzzzz_yyzzz, ts_zzzzzz_yzzzz, ts_zzzzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzzzz_xxxxx[i] = 5.0 * ts_zzzz_xxxxx[i] * fe_0 + ts_zzzzz_xxxxx[i] * pa_z[i];

        ts_zzzzzz_xxxxy[i] = 5.0 * ts_zzzz_xxxxy[i] * fe_0 + ts_zzzzz_xxxxy[i] * pa_z[i];

        ts_zzzzzz_xxxxz[i] = 5.0 * ts_zzzz_xxxxz[i] * fe_0 + ts_zzzzz_xxxx[i] * fe_0 + ts_zzzzz_xxxxz[i] * pa_z[i];

        ts_zzzzzz_xxxyy[i] = 5.0 * ts_zzzz_xxxyy[i] * fe_0 + ts_zzzzz_xxxyy[i] * pa_z[i];

        ts_zzzzzz_xxxyz[i] = 5.0 * ts_zzzz_xxxyz[i] * fe_0 + ts_zzzzz_xxxy[i] * fe_0 + ts_zzzzz_xxxyz[i] * pa_z[i];

        ts_zzzzzz_xxxzz[i] = 5.0 * ts_zzzz_xxxzz[i] * fe_0 + 2.0 * ts_zzzzz_xxxz[i] * fe_0 + ts_zzzzz_xxxzz[i] * pa_z[i];

        ts_zzzzzz_xxyyy[i] = 5.0 * ts_zzzz_xxyyy[i] * fe_0 + ts_zzzzz_xxyyy[i] * pa_z[i];

        ts_zzzzzz_xxyyz[i] = 5.0 * ts_zzzz_xxyyz[i] * fe_0 + ts_zzzzz_xxyy[i] * fe_0 + ts_zzzzz_xxyyz[i] * pa_z[i];

        ts_zzzzzz_xxyzz[i] = 5.0 * ts_zzzz_xxyzz[i] * fe_0 + 2.0 * ts_zzzzz_xxyz[i] * fe_0 + ts_zzzzz_xxyzz[i] * pa_z[i];

        ts_zzzzzz_xxzzz[i] = 5.0 * ts_zzzz_xxzzz[i] * fe_0 + 3.0 * ts_zzzzz_xxzz[i] * fe_0 + ts_zzzzz_xxzzz[i] * pa_z[i];

        ts_zzzzzz_xyyyy[i] = 5.0 * ts_zzzz_xyyyy[i] * fe_0 + ts_zzzzz_xyyyy[i] * pa_z[i];

        ts_zzzzzz_xyyyz[i] = 5.0 * ts_zzzz_xyyyz[i] * fe_0 + ts_zzzzz_xyyy[i] * fe_0 + ts_zzzzz_xyyyz[i] * pa_z[i];

        ts_zzzzzz_xyyzz[i] = 5.0 * ts_zzzz_xyyzz[i] * fe_0 + 2.0 * ts_zzzzz_xyyz[i] * fe_0 + ts_zzzzz_xyyzz[i] * pa_z[i];

        ts_zzzzzz_xyzzz[i] = 5.0 * ts_zzzz_xyzzz[i] * fe_0 + 3.0 * ts_zzzzz_xyzz[i] * fe_0 + ts_zzzzz_xyzzz[i] * pa_z[i];

        ts_zzzzzz_xzzzz[i] = 5.0 * ts_zzzz_xzzzz[i] * fe_0 + 4.0 * ts_zzzzz_xzzz[i] * fe_0 + ts_zzzzz_xzzzz[i] * pa_z[i];

        ts_zzzzzz_yyyyy[i] = 5.0 * ts_zzzz_yyyyy[i] * fe_0 + ts_zzzzz_yyyyy[i] * pa_z[i];

        ts_zzzzzz_yyyyz[i] = 5.0 * ts_zzzz_yyyyz[i] * fe_0 + ts_zzzzz_yyyy[i] * fe_0 + ts_zzzzz_yyyyz[i] * pa_z[i];

        ts_zzzzzz_yyyzz[i] = 5.0 * ts_zzzz_yyyzz[i] * fe_0 + 2.0 * ts_zzzzz_yyyz[i] * fe_0 + ts_zzzzz_yyyzz[i] * pa_z[i];

        ts_zzzzzz_yyzzz[i] = 5.0 * ts_zzzz_yyzzz[i] * fe_0 + 3.0 * ts_zzzzz_yyzz[i] * fe_0 + ts_zzzzz_yyzzz[i] * pa_z[i];

        ts_zzzzzz_yzzzz[i] = 5.0 * ts_zzzz_yzzzz[i] * fe_0 + 4.0 * ts_zzzzz_yzzz[i] * fe_0 + ts_zzzzz_yzzzz[i] * pa_z[i];

        ts_zzzzzz_zzzzz[i] = 5.0 * ts_zzzz_zzzzz[i] * fe_0 + 5.0 * ts_zzzzz_zzzz[i] * fe_0 + ts_zzzzz_zzzzz[i] * pa_z[i];
    }

}

} // ovlrec namespace

