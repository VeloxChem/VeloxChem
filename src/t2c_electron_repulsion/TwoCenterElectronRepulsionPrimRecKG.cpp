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

#include "TwoCenterElectronRepulsionPrimRecKG.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_kg(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_kg,
                                const size_t idx_eri_0_hg,
                                const size_t idx_eri_1_hg,
                                const size_t idx_eri_1_if,
                                const size_t idx_eri_1_ig,
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

    // Set up components of auxiliary buffer : HG

    auto g_xxxxx_xxxx_0 = pbuffer.data(idx_eri_0_hg);

    auto g_xxxxx_xxxy_0 = pbuffer.data(idx_eri_0_hg + 1);

    auto g_xxxxx_xxxz_0 = pbuffer.data(idx_eri_0_hg + 2);

    auto g_xxxxx_xxyy_0 = pbuffer.data(idx_eri_0_hg + 3);

    auto g_xxxxx_xxyz_0 = pbuffer.data(idx_eri_0_hg + 4);

    auto g_xxxxx_xxzz_0 = pbuffer.data(idx_eri_0_hg + 5);

    auto g_xxxxx_xyyy_0 = pbuffer.data(idx_eri_0_hg + 6);

    auto g_xxxxx_xyyz_0 = pbuffer.data(idx_eri_0_hg + 7);

    auto g_xxxxx_xyzz_0 = pbuffer.data(idx_eri_0_hg + 8);

    auto g_xxxxx_xzzz_0 = pbuffer.data(idx_eri_0_hg + 9);

    auto g_xxxxx_yyyy_0 = pbuffer.data(idx_eri_0_hg + 10);

    auto g_xxxxx_yyyz_0 = pbuffer.data(idx_eri_0_hg + 11);

    auto g_xxxxx_yyzz_0 = pbuffer.data(idx_eri_0_hg + 12);

    auto g_xxxxx_yzzz_0 = pbuffer.data(idx_eri_0_hg + 13);

    auto g_xxxxx_zzzz_0 = pbuffer.data(idx_eri_0_hg + 14);

    auto g_xxxxy_xxxx_0 = pbuffer.data(idx_eri_0_hg + 15);

    auto g_xxxxy_xxxz_0 = pbuffer.data(idx_eri_0_hg + 17);

    auto g_xxxxy_xxzz_0 = pbuffer.data(idx_eri_0_hg + 20);

    auto g_xxxxy_xzzz_0 = pbuffer.data(idx_eri_0_hg + 24);

    auto g_xxxxz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 30);

    auto g_xxxxz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 31);

    auto g_xxxxz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 33);

    auto g_xxxxz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 36);

    auto g_xxxyy_xxxx_0 = pbuffer.data(idx_eri_0_hg + 45);

    auto g_xxxyy_xxxy_0 = pbuffer.data(idx_eri_0_hg + 46);

    auto g_xxxyy_xxxz_0 = pbuffer.data(idx_eri_0_hg + 47);

    auto g_xxxyy_xxyy_0 = pbuffer.data(idx_eri_0_hg + 48);

    auto g_xxxyy_xxyz_0 = pbuffer.data(idx_eri_0_hg + 49);

    auto g_xxxyy_xxzz_0 = pbuffer.data(idx_eri_0_hg + 50);

    auto g_xxxyy_xyyy_0 = pbuffer.data(idx_eri_0_hg + 51);

    auto g_xxxyy_xyyz_0 = pbuffer.data(idx_eri_0_hg + 52);

    auto g_xxxyy_xyzz_0 = pbuffer.data(idx_eri_0_hg + 53);

    auto g_xxxyy_xzzz_0 = pbuffer.data(idx_eri_0_hg + 54);

    auto g_xxxyy_yyyy_0 = pbuffer.data(idx_eri_0_hg + 55);

    auto g_xxxyy_yyyz_0 = pbuffer.data(idx_eri_0_hg + 56);

    auto g_xxxyy_yyzz_0 = pbuffer.data(idx_eri_0_hg + 57);

    auto g_xxxyy_yzzz_0 = pbuffer.data(idx_eri_0_hg + 58);

    auto g_xxxyy_zzzz_0 = pbuffer.data(idx_eri_0_hg + 59);

    auto g_xxxzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 75);

    auto g_xxxzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 76);

    auto g_xxxzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 77);

    auto g_xxxzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 78);

    auto g_xxxzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 79);

    auto g_xxxzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 80);

    auto g_xxxzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 81);

    auto g_xxxzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 82);

    auto g_xxxzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 83);

    auto g_xxxzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 84);

    auto g_xxxzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 85);

    auto g_xxxzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 86);

    auto g_xxxzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 87);

    auto g_xxxzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 88);

    auto g_xxxzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 89);

    auto g_xxyyy_xxxx_0 = pbuffer.data(idx_eri_0_hg + 90);

    auto g_xxyyy_xxxy_0 = pbuffer.data(idx_eri_0_hg + 91);

    auto g_xxyyy_xxxz_0 = pbuffer.data(idx_eri_0_hg + 92);

    auto g_xxyyy_xxyy_0 = pbuffer.data(idx_eri_0_hg + 93);

    auto g_xxyyy_xxyz_0 = pbuffer.data(idx_eri_0_hg + 94);

    auto g_xxyyy_xxzz_0 = pbuffer.data(idx_eri_0_hg + 95);

    auto g_xxyyy_xyyy_0 = pbuffer.data(idx_eri_0_hg + 96);

    auto g_xxyyy_xyyz_0 = pbuffer.data(idx_eri_0_hg + 97);

    auto g_xxyyy_xyzz_0 = pbuffer.data(idx_eri_0_hg + 98);

    auto g_xxyyy_xzzz_0 = pbuffer.data(idx_eri_0_hg + 99);

    auto g_xxyyy_yyyy_0 = pbuffer.data(idx_eri_0_hg + 100);

    auto g_xxyyy_yyyz_0 = pbuffer.data(idx_eri_0_hg + 101);

    auto g_xxyyy_yyzz_0 = pbuffer.data(idx_eri_0_hg + 102);

    auto g_xxyyy_yzzz_0 = pbuffer.data(idx_eri_0_hg + 103);

    auto g_xxyyy_zzzz_0 = pbuffer.data(idx_eri_0_hg + 104);

    auto g_xxyyz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 106);

    auto g_xxyyz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 108);

    auto g_xxyyz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 111);

    auto g_xxyzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 120);

    auto g_xxyzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 122);

    auto g_xxyzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 125);

    auto g_xxyzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 129);

    auto g_xxzzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 135);

    auto g_xxzzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 136);

    auto g_xxzzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 137);

    auto g_xxzzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 138);

    auto g_xxzzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 139);

    auto g_xxzzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 140);

    auto g_xxzzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 141);

    auto g_xxzzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 142);

    auto g_xxzzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 143);

    auto g_xxzzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 144);

    auto g_xxzzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 145);

    auto g_xxzzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 146);

    auto g_xxzzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 147);

    auto g_xxzzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 148);

    auto g_xxzzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 149);

    auto g_xyyyy_xxxy_0 = pbuffer.data(idx_eri_0_hg + 151);

    auto g_xyyyy_xxyy_0 = pbuffer.data(idx_eri_0_hg + 153);

    auto g_xyyyy_xxyz_0 = pbuffer.data(idx_eri_0_hg + 154);

    auto g_xyyyy_xyyy_0 = pbuffer.data(idx_eri_0_hg + 156);

    auto g_xyyyy_xyyz_0 = pbuffer.data(idx_eri_0_hg + 157);

    auto g_xyyyy_xyzz_0 = pbuffer.data(idx_eri_0_hg + 158);

    auto g_xyyyy_yyyy_0 = pbuffer.data(idx_eri_0_hg + 160);

    auto g_xyyyy_yyyz_0 = pbuffer.data(idx_eri_0_hg + 161);

    auto g_xyyyy_yyzz_0 = pbuffer.data(idx_eri_0_hg + 162);

    auto g_xyyyy_yzzz_0 = pbuffer.data(idx_eri_0_hg + 163);

    auto g_xyyyy_zzzz_0 = pbuffer.data(idx_eri_0_hg + 164);

    auto g_xyyzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 184);

    auto g_xyyzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 187);

    auto g_xyyzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 188);

    auto g_xyyzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 190);

    auto g_xyyzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 191);

    auto g_xyyzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 192);

    auto g_xyyzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 193);

    auto g_xyyzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 194);

    auto g_xzzzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 212);

    auto g_xzzzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 214);

    auto g_xzzzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 215);

    auto g_xzzzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 217);

    auto g_xzzzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 218);

    auto g_xzzzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 219);

    auto g_xzzzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 220);

    auto g_xzzzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 221);

    auto g_xzzzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 222);

    auto g_xzzzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 223);

    auto g_xzzzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 224);

    auto g_yyyyy_xxxx_0 = pbuffer.data(idx_eri_0_hg + 225);

    auto g_yyyyy_xxxy_0 = pbuffer.data(idx_eri_0_hg + 226);

    auto g_yyyyy_xxxz_0 = pbuffer.data(idx_eri_0_hg + 227);

    auto g_yyyyy_xxyy_0 = pbuffer.data(idx_eri_0_hg + 228);

    auto g_yyyyy_xxyz_0 = pbuffer.data(idx_eri_0_hg + 229);

    auto g_yyyyy_xxzz_0 = pbuffer.data(idx_eri_0_hg + 230);

    auto g_yyyyy_xyyy_0 = pbuffer.data(idx_eri_0_hg + 231);

    auto g_yyyyy_xyyz_0 = pbuffer.data(idx_eri_0_hg + 232);

    auto g_yyyyy_xyzz_0 = pbuffer.data(idx_eri_0_hg + 233);

    auto g_yyyyy_xzzz_0 = pbuffer.data(idx_eri_0_hg + 234);

    auto g_yyyyy_yyyy_0 = pbuffer.data(idx_eri_0_hg + 235);

    auto g_yyyyy_yyyz_0 = pbuffer.data(idx_eri_0_hg + 236);

    auto g_yyyyy_yyzz_0 = pbuffer.data(idx_eri_0_hg + 237);

    auto g_yyyyy_yzzz_0 = pbuffer.data(idx_eri_0_hg + 238);

    auto g_yyyyy_zzzz_0 = pbuffer.data(idx_eri_0_hg + 239);

    auto g_yyyyz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 241);

    auto g_yyyyz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 243);

    auto g_yyyyz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 246);

    auto g_yyyyz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 250);

    auto g_yyyzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 255);

    auto g_yyyzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 256);

    auto g_yyyzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 257);

    auto g_yyyzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 258);

    auto g_yyyzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 259);

    auto g_yyyzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 260);

    auto g_yyyzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 261);

    auto g_yyyzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 262);

    auto g_yyyzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 263);

    auto g_yyyzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 264);

    auto g_yyyzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 265);

    auto g_yyyzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 266);

    auto g_yyyzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 267);

    auto g_yyyzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 268);

    auto g_yyyzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 269);

    auto g_yyzzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 270);

    auto g_yyzzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 271);

    auto g_yyzzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 272);

    auto g_yyzzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 273);

    auto g_yyzzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 274);

    auto g_yyzzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 275);

    auto g_yyzzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 276);

    auto g_yyzzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 277);

    auto g_yyzzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 278);

    auto g_yyzzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 279);

    auto g_yyzzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 280);

    auto g_yyzzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 281);

    auto g_yyzzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 282);

    auto g_yyzzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 283);

    auto g_yyzzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 284);

    auto g_yzzzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 285);

    auto g_yzzzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 287);

    auto g_yzzzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 289);

    auto g_yzzzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 290);

    auto g_yzzzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 292);

    auto g_yzzzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 293);

    auto g_yzzzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 294);

    auto g_yzzzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 296);

    auto g_yzzzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 297);

    auto g_yzzzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 298);

    auto g_yzzzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 299);

    auto g_zzzzz_xxxx_0 = pbuffer.data(idx_eri_0_hg + 300);

    auto g_zzzzz_xxxy_0 = pbuffer.data(idx_eri_0_hg + 301);

    auto g_zzzzz_xxxz_0 = pbuffer.data(idx_eri_0_hg + 302);

    auto g_zzzzz_xxyy_0 = pbuffer.data(idx_eri_0_hg + 303);

    auto g_zzzzz_xxyz_0 = pbuffer.data(idx_eri_0_hg + 304);

    auto g_zzzzz_xxzz_0 = pbuffer.data(idx_eri_0_hg + 305);

    auto g_zzzzz_xyyy_0 = pbuffer.data(idx_eri_0_hg + 306);

    auto g_zzzzz_xyyz_0 = pbuffer.data(idx_eri_0_hg + 307);

    auto g_zzzzz_xyzz_0 = pbuffer.data(idx_eri_0_hg + 308);

    auto g_zzzzz_xzzz_0 = pbuffer.data(idx_eri_0_hg + 309);

    auto g_zzzzz_yyyy_0 = pbuffer.data(idx_eri_0_hg + 310);

    auto g_zzzzz_yyyz_0 = pbuffer.data(idx_eri_0_hg + 311);

    auto g_zzzzz_yyzz_0 = pbuffer.data(idx_eri_0_hg + 312);

    auto g_zzzzz_yzzz_0 = pbuffer.data(idx_eri_0_hg + 313);

    auto g_zzzzz_zzzz_0 = pbuffer.data(idx_eri_0_hg + 314);

    // Set up components of auxiliary buffer : HG

    auto g_xxxxx_xxxx_1 = pbuffer.data(idx_eri_1_hg);

    auto g_xxxxx_xxxy_1 = pbuffer.data(idx_eri_1_hg + 1);

    auto g_xxxxx_xxxz_1 = pbuffer.data(idx_eri_1_hg + 2);

    auto g_xxxxx_xxyy_1 = pbuffer.data(idx_eri_1_hg + 3);

    auto g_xxxxx_xxyz_1 = pbuffer.data(idx_eri_1_hg + 4);

    auto g_xxxxx_xxzz_1 = pbuffer.data(idx_eri_1_hg + 5);

    auto g_xxxxx_xyyy_1 = pbuffer.data(idx_eri_1_hg + 6);

    auto g_xxxxx_xyyz_1 = pbuffer.data(idx_eri_1_hg + 7);

    auto g_xxxxx_xyzz_1 = pbuffer.data(idx_eri_1_hg + 8);

    auto g_xxxxx_xzzz_1 = pbuffer.data(idx_eri_1_hg + 9);

    auto g_xxxxx_yyyy_1 = pbuffer.data(idx_eri_1_hg + 10);

    auto g_xxxxx_yyyz_1 = pbuffer.data(idx_eri_1_hg + 11);

    auto g_xxxxx_yyzz_1 = pbuffer.data(idx_eri_1_hg + 12);

    auto g_xxxxx_yzzz_1 = pbuffer.data(idx_eri_1_hg + 13);

    auto g_xxxxx_zzzz_1 = pbuffer.data(idx_eri_1_hg + 14);

    auto g_xxxxy_xxxx_1 = pbuffer.data(idx_eri_1_hg + 15);

    auto g_xxxxy_xxxz_1 = pbuffer.data(idx_eri_1_hg + 17);

    auto g_xxxxy_xxzz_1 = pbuffer.data(idx_eri_1_hg + 20);

    auto g_xxxxy_xzzz_1 = pbuffer.data(idx_eri_1_hg + 24);

    auto g_xxxxz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 30);

    auto g_xxxxz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 31);

    auto g_xxxxz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 33);

    auto g_xxxxz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 36);

    auto g_xxxyy_xxxx_1 = pbuffer.data(idx_eri_1_hg + 45);

    auto g_xxxyy_xxxy_1 = pbuffer.data(idx_eri_1_hg + 46);

    auto g_xxxyy_xxxz_1 = pbuffer.data(idx_eri_1_hg + 47);

    auto g_xxxyy_xxyy_1 = pbuffer.data(idx_eri_1_hg + 48);

    auto g_xxxyy_xxyz_1 = pbuffer.data(idx_eri_1_hg + 49);

    auto g_xxxyy_xxzz_1 = pbuffer.data(idx_eri_1_hg + 50);

    auto g_xxxyy_xyyy_1 = pbuffer.data(idx_eri_1_hg + 51);

    auto g_xxxyy_xyyz_1 = pbuffer.data(idx_eri_1_hg + 52);

    auto g_xxxyy_xyzz_1 = pbuffer.data(idx_eri_1_hg + 53);

    auto g_xxxyy_xzzz_1 = pbuffer.data(idx_eri_1_hg + 54);

    auto g_xxxyy_yyyy_1 = pbuffer.data(idx_eri_1_hg + 55);

    auto g_xxxyy_yyyz_1 = pbuffer.data(idx_eri_1_hg + 56);

    auto g_xxxyy_yyzz_1 = pbuffer.data(idx_eri_1_hg + 57);

    auto g_xxxyy_yzzz_1 = pbuffer.data(idx_eri_1_hg + 58);

    auto g_xxxyy_zzzz_1 = pbuffer.data(idx_eri_1_hg + 59);

    auto g_xxxzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 75);

    auto g_xxxzz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 76);

    auto g_xxxzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 77);

    auto g_xxxzz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 78);

    auto g_xxxzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 79);

    auto g_xxxzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 80);

    auto g_xxxzz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 81);

    auto g_xxxzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 82);

    auto g_xxxzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 83);

    auto g_xxxzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 84);

    auto g_xxxzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 85);

    auto g_xxxzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 86);

    auto g_xxxzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 87);

    auto g_xxxzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 88);

    auto g_xxxzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 89);

    auto g_xxyyy_xxxx_1 = pbuffer.data(idx_eri_1_hg + 90);

    auto g_xxyyy_xxxy_1 = pbuffer.data(idx_eri_1_hg + 91);

    auto g_xxyyy_xxxz_1 = pbuffer.data(idx_eri_1_hg + 92);

    auto g_xxyyy_xxyy_1 = pbuffer.data(idx_eri_1_hg + 93);

    auto g_xxyyy_xxyz_1 = pbuffer.data(idx_eri_1_hg + 94);

    auto g_xxyyy_xxzz_1 = pbuffer.data(idx_eri_1_hg + 95);

    auto g_xxyyy_xyyy_1 = pbuffer.data(idx_eri_1_hg + 96);

    auto g_xxyyy_xyyz_1 = pbuffer.data(idx_eri_1_hg + 97);

    auto g_xxyyy_xyzz_1 = pbuffer.data(idx_eri_1_hg + 98);

    auto g_xxyyy_xzzz_1 = pbuffer.data(idx_eri_1_hg + 99);

    auto g_xxyyy_yyyy_1 = pbuffer.data(idx_eri_1_hg + 100);

    auto g_xxyyy_yyyz_1 = pbuffer.data(idx_eri_1_hg + 101);

    auto g_xxyyy_yyzz_1 = pbuffer.data(idx_eri_1_hg + 102);

    auto g_xxyyy_yzzz_1 = pbuffer.data(idx_eri_1_hg + 103);

    auto g_xxyyy_zzzz_1 = pbuffer.data(idx_eri_1_hg + 104);

    auto g_xxyyz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 106);

    auto g_xxyyz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 108);

    auto g_xxyyz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 111);

    auto g_xxyzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 120);

    auto g_xxyzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 122);

    auto g_xxyzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 125);

    auto g_xxyzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 129);

    auto g_xxzzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 135);

    auto g_xxzzz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 136);

    auto g_xxzzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 137);

    auto g_xxzzz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 138);

    auto g_xxzzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 139);

    auto g_xxzzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 140);

    auto g_xxzzz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 141);

    auto g_xxzzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 142);

    auto g_xxzzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 143);

    auto g_xxzzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 144);

    auto g_xxzzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 145);

    auto g_xxzzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 146);

    auto g_xxzzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 147);

    auto g_xxzzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 148);

    auto g_xxzzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 149);

    auto g_xyyyy_xxxy_1 = pbuffer.data(idx_eri_1_hg + 151);

    auto g_xyyyy_xxyy_1 = pbuffer.data(idx_eri_1_hg + 153);

    auto g_xyyyy_xxyz_1 = pbuffer.data(idx_eri_1_hg + 154);

    auto g_xyyyy_xyyy_1 = pbuffer.data(idx_eri_1_hg + 156);

    auto g_xyyyy_xyyz_1 = pbuffer.data(idx_eri_1_hg + 157);

    auto g_xyyyy_xyzz_1 = pbuffer.data(idx_eri_1_hg + 158);

    auto g_xyyyy_yyyy_1 = pbuffer.data(idx_eri_1_hg + 160);

    auto g_xyyyy_yyyz_1 = pbuffer.data(idx_eri_1_hg + 161);

    auto g_xyyyy_yyzz_1 = pbuffer.data(idx_eri_1_hg + 162);

    auto g_xyyyy_yzzz_1 = pbuffer.data(idx_eri_1_hg + 163);

    auto g_xyyyy_zzzz_1 = pbuffer.data(idx_eri_1_hg + 164);

    auto g_xyyzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 184);

    auto g_xyyzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 187);

    auto g_xyyzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 188);

    auto g_xyyzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 190);

    auto g_xyyzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 191);

    auto g_xyyzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 192);

    auto g_xyyzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 193);

    auto g_xyyzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 194);

    auto g_xzzzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 212);

    auto g_xzzzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 214);

    auto g_xzzzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 215);

    auto g_xzzzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 217);

    auto g_xzzzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 218);

    auto g_xzzzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 219);

    auto g_xzzzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 220);

    auto g_xzzzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 221);

    auto g_xzzzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 222);

    auto g_xzzzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 223);

    auto g_xzzzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 224);

    auto g_yyyyy_xxxx_1 = pbuffer.data(idx_eri_1_hg + 225);

    auto g_yyyyy_xxxy_1 = pbuffer.data(idx_eri_1_hg + 226);

    auto g_yyyyy_xxxz_1 = pbuffer.data(idx_eri_1_hg + 227);

    auto g_yyyyy_xxyy_1 = pbuffer.data(idx_eri_1_hg + 228);

    auto g_yyyyy_xxyz_1 = pbuffer.data(idx_eri_1_hg + 229);

    auto g_yyyyy_xxzz_1 = pbuffer.data(idx_eri_1_hg + 230);

    auto g_yyyyy_xyyy_1 = pbuffer.data(idx_eri_1_hg + 231);

    auto g_yyyyy_xyyz_1 = pbuffer.data(idx_eri_1_hg + 232);

    auto g_yyyyy_xyzz_1 = pbuffer.data(idx_eri_1_hg + 233);

    auto g_yyyyy_xzzz_1 = pbuffer.data(idx_eri_1_hg + 234);

    auto g_yyyyy_yyyy_1 = pbuffer.data(idx_eri_1_hg + 235);

    auto g_yyyyy_yyyz_1 = pbuffer.data(idx_eri_1_hg + 236);

    auto g_yyyyy_yyzz_1 = pbuffer.data(idx_eri_1_hg + 237);

    auto g_yyyyy_yzzz_1 = pbuffer.data(idx_eri_1_hg + 238);

    auto g_yyyyy_zzzz_1 = pbuffer.data(idx_eri_1_hg + 239);

    auto g_yyyyz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 241);

    auto g_yyyyz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 243);

    auto g_yyyyz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 246);

    auto g_yyyyz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 250);

    auto g_yyyzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 255);

    auto g_yyyzz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 256);

    auto g_yyyzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 257);

    auto g_yyyzz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 258);

    auto g_yyyzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 259);

    auto g_yyyzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 260);

    auto g_yyyzz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 261);

    auto g_yyyzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 262);

    auto g_yyyzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 263);

    auto g_yyyzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 264);

    auto g_yyyzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 265);

    auto g_yyyzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 266);

    auto g_yyyzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 267);

    auto g_yyyzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 268);

    auto g_yyyzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 269);

    auto g_yyzzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 270);

    auto g_yyzzz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 271);

    auto g_yyzzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 272);

    auto g_yyzzz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 273);

    auto g_yyzzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 274);

    auto g_yyzzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 275);

    auto g_yyzzz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 276);

    auto g_yyzzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 277);

    auto g_yyzzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 278);

    auto g_yyzzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 279);

    auto g_yyzzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 280);

    auto g_yyzzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 281);

    auto g_yyzzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 282);

    auto g_yyzzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 283);

    auto g_yyzzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 284);

    auto g_yzzzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 285);

    auto g_yzzzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 287);

    auto g_yzzzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 289);

    auto g_yzzzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 290);

    auto g_yzzzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 292);

    auto g_yzzzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 293);

    auto g_yzzzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 294);

    auto g_yzzzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 296);

    auto g_yzzzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 297);

    auto g_yzzzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 298);

    auto g_yzzzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 299);

    auto g_zzzzz_xxxx_1 = pbuffer.data(idx_eri_1_hg + 300);

    auto g_zzzzz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 301);

    auto g_zzzzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 302);

    auto g_zzzzz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 303);

    auto g_zzzzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 304);

    auto g_zzzzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 305);

    auto g_zzzzz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 306);

    auto g_zzzzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 307);

    auto g_zzzzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 308);

    auto g_zzzzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 309);

    auto g_zzzzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 310);

    auto g_zzzzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 311);

    auto g_zzzzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 312);

    auto g_zzzzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 313);

    auto g_zzzzz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 314);

    // Set up components of auxiliary buffer : IF

    auto g_xxxxxx_xxx_1 = pbuffer.data(idx_eri_1_if);

    auto g_xxxxxx_xxy_1 = pbuffer.data(idx_eri_1_if + 1);

    auto g_xxxxxx_xxz_1 = pbuffer.data(idx_eri_1_if + 2);

    auto g_xxxxxx_xyy_1 = pbuffer.data(idx_eri_1_if + 3);

    auto g_xxxxxx_xyz_1 = pbuffer.data(idx_eri_1_if + 4);

    auto g_xxxxxx_xzz_1 = pbuffer.data(idx_eri_1_if + 5);

    auto g_xxxxxx_yyy_1 = pbuffer.data(idx_eri_1_if + 6);

    auto g_xxxxxx_yyz_1 = pbuffer.data(idx_eri_1_if + 7);

    auto g_xxxxxx_yzz_1 = pbuffer.data(idx_eri_1_if + 8);

    auto g_xxxxxx_zzz_1 = pbuffer.data(idx_eri_1_if + 9);

    auto g_xxxxxz_xxz_1 = pbuffer.data(idx_eri_1_if + 22);

    auto g_xxxxxz_xyz_1 = pbuffer.data(idx_eri_1_if + 24);

    auto g_xxxxxz_xzz_1 = pbuffer.data(idx_eri_1_if + 25);

    auto g_xxxxxz_yyz_1 = pbuffer.data(idx_eri_1_if + 27);

    auto g_xxxxxz_yzz_1 = pbuffer.data(idx_eri_1_if + 28);

    auto g_xxxxxz_zzz_1 = pbuffer.data(idx_eri_1_if + 29);

    auto g_xxxxyy_xxx_1 = pbuffer.data(idx_eri_1_if + 30);

    auto g_xxxxyy_xxy_1 = pbuffer.data(idx_eri_1_if + 31);

    auto g_xxxxyy_xxz_1 = pbuffer.data(idx_eri_1_if + 32);

    auto g_xxxxyy_xyy_1 = pbuffer.data(idx_eri_1_if + 33);

    auto g_xxxxyy_xyz_1 = pbuffer.data(idx_eri_1_if + 34);

    auto g_xxxxyy_xzz_1 = pbuffer.data(idx_eri_1_if + 35);

    auto g_xxxxyy_yyy_1 = pbuffer.data(idx_eri_1_if + 36);

    auto g_xxxxyy_yyz_1 = pbuffer.data(idx_eri_1_if + 37);

    auto g_xxxxyy_yzz_1 = pbuffer.data(idx_eri_1_if + 38);

    auto g_xxxxyy_zzz_1 = pbuffer.data(idx_eri_1_if + 39);

    auto g_xxxxzz_xxx_1 = pbuffer.data(idx_eri_1_if + 50);

    auto g_xxxxzz_xxy_1 = pbuffer.data(idx_eri_1_if + 51);

    auto g_xxxxzz_xxz_1 = pbuffer.data(idx_eri_1_if + 52);

    auto g_xxxxzz_xyy_1 = pbuffer.data(idx_eri_1_if + 53);

    auto g_xxxxzz_xyz_1 = pbuffer.data(idx_eri_1_if + 54);

    auto g_xxxxzz_xzz_1 = pbuffer.data(idx_eri_1_if + 55);

    auto g_xxxxzz_yyy_1 = pbuffer.data(idx_eri_1_if + 56);

    auto g_xxxxzz_yyz_1 = pbuffer.data(idx_eri_1_if + 57);

    auto g_xxxxzz_yzz_1 = pbuffer.data(idx_eri_1_if + 58);

    auto g_xxxxzz_zzz_1 = pbuffer.data(idx_eri_1_if + 59);

    auto g_xxxyyy_xxx_1 = pbuffer.data(idx_eri_1_if + 60);

    auto g_xxxyyy_xxy_1 = pbuffer.data(idx_eri_1_if + 61);

    auto g_xxxyyy_xxz_1 = pbuffer.data(idx_eri_1_if + 62);

    auto g_xxxyyy_xyy_1 = pbuffer.data(idx_eri_1_if + 63);

    auto g_xxxyyy_xyz_1 = pbuffer.data(idx_eri_1_if + 64);

    auto g_xxxyyy_xzz_1 = pbuffer.data(idx_eri_1_if + 65);

    auto g_xxxyyy_yyy_1 = pbuffer.data(idx_eri_1_if + 66);

    auto g_xxxyyy_yyz_1 = pbuffer.data(idx_eri_1_if + 67);

    auto g_xxxyyy_yzz_1 = pbuffer.data(idx_eri_1_if + 68);

    auto g_xxxyyy_zzz_1 = pbuffer.data(idx_eri_1_if + 69);

    auto g_xxxzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 90);

    auto g_xxxzzz_xxy_1 = pbuffer.data(idx_eri_1_if + 91);

    auto g_xxxzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 92);

    auto g_xxxzzz_xyy_1 = pbuffer.data(idx_eri_1_if + 93);

    auto g_xxxzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 94);

    auto g_xxxzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 95);

    auto g_xxxzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 96);

    auto g_xxxzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 97);

    auto g_xxxzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 98);

    auto g_xxxzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 99);

    auto g_xxyyyy_xxx_1 = pbuffer.data(idx_eri_1_if + 100);

    auto g_xxyyyy_xxy_1 = pbuffer.data(idx_eri_1_if + 101);

    auto g_xxyyyy_xxz_1 = pbuffer.data(idx_eri_1_if + 102);

    auto g_xxyyyy_xyy_1 = pbuffer.data(idx_eri_1_if + 103);

    auto g_xxyyyy_xyz_1 = pbuffer.data(idx_eri_1_if + 104);

    auto g_xxyyyy_xzz_1 = pbuffer.data(idx_eri_1_if + 105);

    auto g_xxyyyy_yyy_1 = pbuffer.data(idx_eri_1_if + 106);

    auto g_xxyyyy_yyz_1 = pbuffer.data(idx_eri_1_if + 107);

    auto g_xxyyyy_yzz_1 = pbuffer.data(idx_eri_1_if + 108);

    auto g_xxyyyy_zzz_1 = pbuffer.data(idx_eri_1_if + 109);

    auto g_xxyyzz_xyz_1 = pbuffer.data(idx_eri_1_if + 124);

    auto g_xxyyzz_yyz_1 = pbuffer.data(idx_eri_1_if + 127);

    auto g_xxyyzz_yzz_1 = pbuffer.data(idx_eri_1_if + 128);

    auto g_xxzzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 140);

    auto g_xxzzzz_xxy_1 = pbuffer.data(idx_eri_1_if + 141);

    auto g_xxzzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 142);

    auto g_xxzzzz_xyy_1 = pbuffer.data(idx_eri_1_if + 143);

    auto g_xxzzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 144);

    auto g_xxzzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 145);

    auto g_xxzzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 146);

    auto g_xxzzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 147);

    auto g_xxzzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 148);

    auto g_xxzzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 149);

    auto g_xyyyyy_xxy_1 = pbuffer.data(idx_eri_1_if + 151);

    auto g_xyyyyy_xyy_1 = pbuffer.data(idx_eri_1_if + 153);

    auto g_xyyyyy_xyz_1 = pbuffer.data(idx_eri_1_if + 154);

    auto g_xyyyyy_yyy_1 = pbuffer.data(idx_eri_1_if + 156);

    auto g_xyyyyy_yyz_1 = pbuffer.data(idx_eri_1_if + 157);

    auto g_xyyyyy_yzz_1 = pbuffer.data(idx_eri_1_if + 158);

    auto g_xyyyzz_xyz_1 = pbuffer.data(idx_eri_1_if + 174);

    auto g_xyyyzz_yyz_1 = pbuffer.data(idx_eri_1_if + 177);

    auto g_xyyyzz_yzz_1 = pbuffer.data(idx_eri_1_if + 178);

    auto g_xyyzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 184);

    auto g_xyyzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 187);

    auto g_xyyzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 188);

    auto g_xzzzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 202);

    auto g_xzzzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 204);

    auto g_xzzzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 205);

    auto g_xzzzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 207);

    auto g_xzzzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 208);

    auto g_xzzzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 209);

    auto g_yyyyyy_xxx_1 = pbuffer.data(idx_eri_1_if + 210);

    auto g_yyyyyy_xxy_1 = pbuffer.data(idx_eri_1_if + 211);

    auto g_yyyyyy_xxz_1 = pbuffer.data(idx_eri_1_if + 212);

    auto g_yyyyyy_xyy_1 = pbuffer.data(idx_eri_1_if + 213);

    auto g_yyyyyy_xyz_1 = pbuffer.data(idx_eri_1_if + 214);

    auto g_yyyyyy_xzz_1 = pbuffer.data(idx_eri_1_if + 215);

    auto g_yyyyyy_yyy_1 = pbuffer.data(idx_eri_1_if + 216);

    auto g_yyyyyy_yyz_1 = pbuffer.data(idx_eri_1_if + 217);

    auto g_yyyyyy_yzz_1 = pbuffer.data(idx_eri_1_if + 218);

    auto g_yyyyyy_zzz_1 = pbuffer.data(idx_eri_1_if + 219);

    auto g_yyyyyz_xxz_1 = pbuffer.data(idx_eri_1_if + 222);

    auto g_yyyyyz_xyz_1 = pbuffer.data(idx_eri_1_if + 224);

    auto g_yyyyyz_xzz_1 = pbuffer.data(idx_eri_1_if + 225);

    auto g_yyyyyz_yyz_1 = pbuffer.data(idx_eri_1_if + 227);

    auto g_yyyyyz_yzz_1 = pbuffer.data(idx_eri_1_if + 228);

    auto g_yyyyyz_zzz_1 = pbuffer.data(idx_eri_1_if + 229);

    auto g_yyyyzz_xxx_1 = pbuffer.data(idx_eri_1_if + 230);

    auto g_yyyyzz_xxy_1 = pbuffer.data(idx_eri_1_if + 231);

    auto g_yyyyzz_xxz_1 = pbuffer.data(idx_eri_1_if + 232);

    auto g_yyyyzz_xyy_1 = pbuffer.data(idx_eri_1_if + 233);

    auto g_yyyyzz_xyz_1 = pbuffer.data(idx_eri_1_if + 234);

    auto g_yyyyzz_xzz_1 = pbuffer.data(idx_eri_1_if + 235);

    auto g_yyyyzz_yyy_1 = pbuffer.data(idx_eri_1_if + 236);

    auto g_yyyyzz_yyz_1 = pbuffer.data(idx_eri_1_if + 237);

    auto g_yyyyzz_yzz_1 = pbuffer.data(idx_eri_1_if + 238);

    auto g_yyyyzz_zzz_1 = pbuffer.data(idx_eri_1_if + 239);

    auto g_yyyzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 240);

    auto g_yyyzzz_xxy_1 = pbuffer.data(idx_eri_1_if + 241);

    auto g_yyyzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 242);

    auto g_yyyzzz_xyy_1 = pbuffer.data(idx_eri_1_if + 243);

    auto g_yyyzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 244);

    auto g_yyyzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 245);

    auto g_yyyzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 246);

    auto g_yyyzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 247);

    auto g_yyyzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 248);

    auto g_yyyzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 249);

    auto g_yyzzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 250);

    auto g_yyzzzz_xxy_1 = pbuffer.data(idx_eri_1_if + 251);

    auto g_yyzzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 252);

    auto g_yyzzzz_xyy_1 = pbuffer.data(idx_eri_1_if + 253);

    auto g_yyzzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 254);

    auto g_yyzzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 255);

    auto g_yyzzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 256);

    auto g_yyzzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 257);

    auto g_yyzzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 258);

    auto g_yyzzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 259);

    auto g_yzzzzz_xxy_1 = pbuffer.data(idx_eri_1_if + 261);

    auto g_yzzzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 262);

    auto g_yzzzzz_xyy_1 = pbuffer.data(idx_eri_1_if + 263);

    auto g_yzzzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 264);

    auto g_yzzzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 265);

    auto g_yzzzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 266);

    auto g_yzzzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 267);

    auto g_yzzzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 268);

    auto g_yzzzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 269);

    auto g_zzzzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 270);

    auto g_zzzzzz_xxy_1 = pbuffer.data(idx_eri_1_if + 271);

    auto g_zzzzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 272);

    auto g_zzzzzz_xyy_1 = pbuffer.data(idx_eri_1_if + 273);

    auto g_zzzzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 274);

    auto g_zzzzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 275);

    auto g_zzzzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 276);

    auto g_zzzzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 277);

    auto g_zzzzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 278);

    auto g_zzzzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 279);

    // Set up components of auxiliary buffer : IG

    auto g_xxxxxx_xxxx_1 = pbuffer.data(idx_eri_1_ig);

    auto g_xxxxxx_xxxy_1 = pbuffer.data(idx_eri_1_ig + 1);

    auto g_xxxxxx_xxxz_1 = pbuffer.data(idx_eri_1_ig + 2);

    auto g_xxxxxx_xxyy_1 = pbuffer.data(idx_eri_1_ig + 3);

    auto g_xxxxxx_xxyz_1 = pbuffer.data(idx_eri_1_ig + 4);

    auto g_xxxxxx_xxzz_1 = pbuffer.data(idx_eri_1_ig + 5);

    auto g_xxxxxx_xyyy_1 = pbuffer.data(idx_eri_1_ig + 6);

    auto g_xxxxxx_xyyz_1 = pbuffer.data(idx_eri_1_ig + 7);

    auto g_xxxxxx_xyzz_1 = pbuffer.data(idx_eri_1_ig + 8);

    auto g_xxxxxx_xzzz_1 = pbuffer.data(idx_eri_1_ig + 9);

    auto g_xxxxxx_yyyy_1 = pbuffer.data(idx_eri_1_ig + 10);

    auto g_xxxxxx_yyyz_1 = pbuffer.data(idx_eri_1_ig + 11);

    auto g_xxxxxx_yyzz_1 = pbuffer.data(idx_eri_1_ig + 12);

    auto g_xxxxxx_yzzz_1 = pbuffer.data(idx_eri_1_ig + 13);

    auto g_xxxxxx_zzzz_1 = pbuffer.data(idx_eri_1_ig + 14);

    auto g_xxxxxy_xxxx_1 = pbuffer.data(idx_eri_1_ig + 15);

    auto g_xxxxxy_xxxy_1 = pbuffer.data(idx_eri_1_ig + 16);

    auto g_xxxxxy_xxxz_1 = pbuffer.data(idx_eri_1_ig + 17);

    auto g_xxxxxy_xxyy_1 = pbuffer.data(idx_eri_1_ig + 18);

    auto g_xxxxxy_xxzz_1 = pbuffer.data(idx_eri_1_ig + 20);

    auto g_xxxxxy_xyyy_1 = pbuffer.data(idx_eri_1_ig + 21);

    auto g_xxxxxy_xzzz_1 = pbuffer.data(idx_eri_1_ig + 24);

    auto g_xxxxxy_yyyy_1 = pbuffer.data(idx_eri_1_ig + 25);

    auto g_xxxxxz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 30);

    auto g_xxxxxz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 31);

    auto g_xxxxxz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 32);

    auto g_xxxxxz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 33);

    auto g_xxxxxz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 34);

    auto g_xxxxxz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 35);

    auto g_xxxxxz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 36);

    auto g_xxxxxz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 37);

    auto g_xxxxxz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 38);

    auto g_xxxxxz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 39);

    auto g_xxxxxz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 41);

    auto g_xxxxxz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 42);

    auto g_xxxxxz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 43);

    auto g_xxxxxz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 44);

    auto g_xxxxyy_xxxx_1 = pbuffer.data(idx_eri_1_ig + 45);

    auto g_xxxxyy_xxxy_1 = pbuffer.data(idx_eri_1_ig + 46);

    auto g_xxxxyy_xxxz_1 = pbuffer.data(idx_eri_1_ig + 47);

    auto g_xxxxyy_xxyy_1 = pbuffer.data(idx_eri_1_ig + 48);

    auto g_xxxxyy_xxyz_1 = pbuffer.data(idx_eri_1_ig + 49);

    auto g_xxxxyy_xxzz_1 = pbuffer.data(idx_eri_1_ig + 50);

    auto g_xxxxyy_xyyy_1 = pbuffer.data(idx_eri_1_ig + 51);

    auto g_xxxxyy_xyyz_1 = pbuffer.data(idx_eri_1_ig + 52);

    auto g_xxxxyy_xyzz_1 = pbuffer.data(idx_eri_1_ig + 53);

    auto g_xxxxyy_xzzz_1 = pbuffer.data(idx_eri_1_ig + 54);

    auto g_xxxxyy_yyyy_1 = pbuffer.data(idx_eri_1_ig + 55);

    auto g_xxxxyy_yyyz_1 = pbuffer.data(idx_eri_1_ig + 56);

    auto g_xxxxyy_yyzz_1 = pbuffer.data(idx_eri_1_ig + 57);

    auto g_xxxxyy_yzzz_1 = pbuffer.data(idx_eri_1_ig + 58);

    auto g_xxxxyy_zzzz_1 = pbuffer.data(idx_eri_1_ig + 59);

    auto g_xxxxzz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 75);

    auto g_xxxxzz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 76);

    auto g_xxxxzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 77);

    auto g_xxxxzz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 78);

    auto g_xxxxzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 79);

    auto g_xxxxzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 80);

    auto g_xxxxzz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 81);

    auto g_xxxxzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 82);

    auto g_xxxxzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 83);

    auto g_xxxxzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 84);

    auto g_xxxxzz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 85);

    auto g_xxxxzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 86);

    auto g_xxxxzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 87);

    auto g_xxxxzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 88);

    auto g_xxxxzz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 89);

    auto g_xxxyyy_xxxx_1 = pbuffer.data(idx_eri_1_ig + 90);

    auto g_xxxyyy_xxxy_1 = pbuffer.data(idx_eri_1_ig + 91);

    auto g_xxxyyy_xxxz_1 = pbuffer.data(idx_eri_1_ig + 92);

    auto g_xxxyyy_xxyy_1 = pbuffer.data(idx_eri_1_ig + 93);

    auto g_xxxyyy_xxyz_1 = pbuffer.data(idx_eri_1_ig + 94);

    auto g_xxxyyy_xxzz_1 = pbuffer.data(idx_eri_1_ig + 95);

    auto g_xxxyyy_xyyy_1 = pbuffer.data(idx_eri_1_ig + 96);

    auto g_xxxyyy_xyyz_1 = pbuffer.data(idx_eri_1_ig + 97);

    auto g_xxxyyy_xyzz_1 = pbuffer.data(idx_eri_1_ig + 98);

    auto g_xxxyyy_xzzz_1 = pbuffer.data(idx_eri_1_ig + 99);

    auto g_xxxyyy_yyyy_1 = pbuffer.data(idx_eri_1_ig + 100);

    auto g_xxxyyy_yyyz_1 = pbuffer.data(idx_eri_1_ig + 101);

    auto g_xxxyyy_yyzz_1 = pbuffer.data(idx_eri_1_ig + 102);

    auto g_xxxyyy_yzzz_1 = pbuffer.data(idx_eri_1_ig + 103);

    auto g_xxxyyy_zzzz_1 = pbuffer.data(idx_eri_1_ig + 104);

    auto g_xxxyyz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 106);

    auto g_xxxyyz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 108);

    auto g_xxxyyz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 111);

    auto g_xxxyzz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 120);

    auto g_xxxyzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 122);

    auto g_xxxyzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 125);

    auto g_xxxyzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 129);

    auto g_xxxzzz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 135);

    auto g_xxxzzz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 136);

    auto g_xxxzzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 137);

    auto g_xxxzzz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 138);

    auto g_xxxzzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 139);

    auto g_xxxzzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 140);

    auto g_xxxzzz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 141);

    auto g_xxxzzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 142);

    auto g_xxxzzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 143);

    auto g_xxxzzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 144);

    auto g_xxxzzz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 145);

    auto g_xxxzzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 146);

    auto g_xxxzzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 147);

    auto g_xxxzzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 148);

    auto g_xxxzzz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 149);

    auto g_xxyyyy_xxxx_1 = pbuffer.data(idx_eri_1_ig + 150);

    auto g_xxyyyy_xxxy_1 = pbuffer.data(idx_eri_1_ig + 151);

    auto g_xxyyyy_xxxz_1 = pbuffer.data(idx_eri_1_ig + 152);

    auto g_xxyyyy_xxyy_1 = pbuffer.data(idx_eri_1_ig + 153);

    auto g_xxyyyy_xxyz_1 = pbuffer.data(idx_eri_1_ig + 154);

    auto g_xxyyyy_xxzz_1 = pbuffer.data(idx_eri_1_ig + 155);

    auto g_xxyyyy_xyyy_1 = pbuffer.data(idx_eri_1_ig + 156);

    auto g_xxyyyy_xyyz_1 = pbuffer.data(idx_eri_1_ig + 157);

    auto g_xxyyyy_xyzz_1 = pbuffer.data(idx_eri_1_ig + 158);

    auto g_xxyyyy_xzzz_1 = pbuffer.data(idx_eri_1_ig + 159);

    auto g_xxyyyy_yyyy_1 = pbuffer.data(idx_eri_1_ig + 160);

    auto g_xxyyyy_yyyz_1 = pbuffer.data(idx_eri_1_ig + 161);

    auto g_xxyyyy_yyzz_1 = pbuffer.data(idx_eri_1_ig + 162);

    auto g_xxyyyy_yzzz_1 = pbuffer.data(idx_eri_1_ig + 163);

    auto g_xxyyyy_zzzz_1 = pbuffer.data(idx_eri_1_ig + 164);

    auto g_xxyyyz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 166);

    auto g_xxyyyz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 168);

    auto g_xxyyyz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 171);

    auto g_xxyyzz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 180);

    auto g_xxyyzz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 181);

    auto g_xxyyzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 182);

    auto g_xxyyzz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 183);

    auto g_xxyyzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 184);

    auto g_xxyyzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 185);

    auto g_xxyyzz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 186);

    auto g_xxyyzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 187);

    auto g_xxyyzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 188);

    auto g_xxyyzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 189);

    auto g_xxyyzz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 190);

    auto g_xxyyzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 191);

    auto g_xxyyzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 192);

    auto g_xxyyzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 193);

    auto g_xxyyzz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 194);

    auto g_xxyzzz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 195);

    auto g_xxyzzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 197);

    auto g_xxyzzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 200);

    auto g_xxyzzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 204);

    auto g_xxzzzz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 210);

    auto g_xxzzzz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 211);

    auto g_xxzzzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 212);

    auto g_xxzzzz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 213);

    auto g_xxzzzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 214);

    auto g_xxzzzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 215);

    auto g_xxzzzz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 216);

    auto g_xxzzzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 217);

    auto g_xxzzzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 218);

    auto g_xxzzzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 219);

    auto g_xxzzzz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 220);

    auto g_xxzzzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 221);

    auto g_xxzzzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 222);

    auto g_xxzzzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 223);

    auto g_xxzzzz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 224);

    auto g_xyyyyy_xxxx_1 = pbuffer.data(idx_eri_1_ig + 225);

    auto g_xyyyyy_xxxy_1 = pbuffer.data(idx_eri_1_ig + 226);

    auto g_xyyyyy_xxyy_1 = pbuffer.data(idx_eri_1_ig + 228);

    auto g_xyyyyy_xxyz_1 = pbuffer.data(idx_eri_1_ig + 229);

    auto g_xyyyyy_xyyy_1 = pbuffer.data(idx_eri_1_ig + 231);

    auto g_xyyyyy_xyyz_1 = pbuffer.data(idx_eri_1_ig + 232);

    auto g_xyyyyy_xyzz_1 = pbuffer.data(idx_eri_1_ig + 233);

    auto g_xyyyyy_yyyy_1 = pbuffer.data(idx_eri_1_ig + 235);

    auto g_xyyyyy_yyyz_1 = pbuffer.data(idx_eri_1_ig + 236);

    auto g_xyyyyy_yyzz_1 = pbuffer.data(idx_eri_1_ig + 237);

    auto g_xyyyyy_yzzz_1 = pbuffer.data(idx_eri_1_ig + 238);

    auto g_xyyyyy_zzzz_1 = pbuffer.data(idx_eri_1_ig + 239);

    auto g_xyyyzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 259);

    auto g_xyyyzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 262);

    auto g_xyyyzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 263);

    auto g_xyyyzz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 265);

    auto g_xyyyzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 266);

    auto g_xyyyzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 267);

    auto g_xyyyzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 268);

    auto g_xyyyzz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 269);

    auto g_xyyzzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 274);

    auto g_xyyzzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 277);

    auto g_xyyzzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 278);

    auto g_xyyzzz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 280);

    auto g_xyyzzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 281);

    auto g_xyyzzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 282);

    auto g_xyyzzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 283);

    auto g_xyyzzz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 284);

    auto g_xzzzzz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 300);

    auto g_xzzzzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 302);

    auto g_xzzzzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 304);

    auto g_xzzzzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 305);

    auto g_xzzzzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 307);

    auto g_xzzzzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 308);

    auto g_xzzzzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 309);

    auto g_xzzzzz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 310);

    auto g_xzzzzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 311);

    auto g_xzzzzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 312);

    auto g_xzzzzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 313);

    auto g_xzzzzz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 314);

    auto g_yyyyyy_xxxx_1 = pbuffer.data(idx_eri_1_ig + 315);

    auto g_yyyyyy_xxxy_1 = pbuffer.data(idx_eri_1_ig + 316);

    auto g_yyyyyy_xxxz_1 = pbuffer.data(idx_eri_1_ig + 317);

    auto g_yyyyyy_xxyy_1 = pbuffer.data(idx_eri_1_ig + 318);

    auto g_yyyyyy_xxyz_1 = pbuffer.data(idx_eri_1_ig + 319);

    auto g_yyyyyy_xxzz_1 = pbuffer.data(idx_eri_1_ig + 320);

    auto g_yyyyyy_xyyy_1 = pbuffer.data(idx_eri_1_ig + 321);

    auto g_yyyyyy_xyyz_1 = pbuffer.data(idx_eri_1_ig + 322);

    auto g_yyyyyy_xyzz_1 = pbuffer.data(idx_eri_1_ig + 323);

    auto g_yyyyyy_xzzz_1 = pbuffer.data(idx_eri_1_ig + 324);

    auto g_yyyyyy_yyyy_1 = pbuffer.data(idx_eri_1_ig + 325);

    auto g_yyyyyy_yyyz_1 = pbuffer.data(idx_eri_1_ig + 326);

    auto g_yyyyyy_yyzz_1 = pbuffer.data(idx_eri_1_ig + 327);

    auto g_yyyyyy_yzzz_1 = pbuffer.data(idx_eri_1_ig + 328);

    auto g_yyyyyy_zzzz_1 = pbuffer.data(idx_eri_1_ig + 329);

    auto g_yyyyyz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 331);

    auto g_yyyyyz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 332);

    auto g_yyyyyz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 333);

    auto g_yyyyyz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 334);

    auto g_yyyyyz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 335);

    auto g_yyyyyz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 336);

    auto g_yyyyyz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 337);

    auto g_yyyyyz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 338);

    auto g_yyyyyz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 339);

    auto g_yyyyyz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 340);

    auto g_yyyyyz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 341);

    auto g_yyyyyz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 342);

    auto g_yyyyyz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 343);

    auto g_yyyyyz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 344);

    auto g_yyyyzz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 345);

    auto g_yyyyzz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 346);

    auto g_yyyyzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 347);

    auto g_yyyyzz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 348);

    auto g_yyyyzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 349);

    auto g_yyyyzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 350);

    auto g_yyyyzz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 351);

    auto g_yyyyzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 352);

    auto g_yyyyzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 353);

    auto g_yyyyzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 354);

    auto g_yyyyzz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 355);

    auto g_yyyyzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 356);

    auto g_yyyyzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 357);

    auto g_yyyyzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 358);

    auto g_yyyyzz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 359);

    auto g_yyyzzz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 360);

    auto g_yyyzzz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 361);

    auto g_yyyzzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 362);

    auto g_yyyzzz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 363);

    auto g_yyyzzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 364);

    auto g_yyyzzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 365);

    auto g_yyyzzz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 366);

    auto g_yyyzzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 367);

    auto g_yyyzzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 368);

    auto g_yyyzzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 369);

    auto g_yyyzzz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 370);

    auto g_yyyzzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 371);

    auto g_yyyzzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 372);

    auto g_yyyzzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 373);

    auto g_yyyzzz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 374);

    auto g_yyzzzz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 375);

    auto g_yyzzzz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 376);

    auto g_yyzzzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 377);

    auto g_yyzzzz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 378);

    auto g_yyzzzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 379);

    auto g_yyzzzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 380);

    auto g_yyzzzz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 381);

    auto g_yyzzzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 382);

    auto g_yyzzzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 383);

    auto g_yyzzzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 384);

    auto g_yyzzzz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 385);

    auto g_yyzzzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 386);

    auto g_yyzzzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 387);

    auto g_yyzzzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 388);

    auto g_yyzzzz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 389);

    auto g_yzzzzz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 390);

    auto g_yzzzzz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 391);

    auto g_yzzzzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 392);

    auto g_yzzzzz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 393);

    auto g_yzzzzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 394);

    auto g_yzzzzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 395);

    auto g_yzzzzz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 396);

    auto g_yzzzzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 397);

    auto g_yzzzzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 398);

    auto g_yzzzzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 399);

    auto g_yzzzzz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 400);

    auto g_yzzzzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 401);

    auto g_yzzzzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 402);

    auto g_yzzzzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 403);

    auto g_yzzzzz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 404);

    auto g_zzzzzz_xxxx_1 = pbuffer.data(idx_eri_1_ig + 405);

    auto g_zzzzzz_xxxy_1 = pbuffer.data(idx_eri_1_ig + 406);

    auto g_zzzzzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 407);

    auto g_zzzzzz_xxyy_1 = pbuffer.data(idx_eri_1_ig + 408);

    auto g_zzzzzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 409);

    auto g_zzzzzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 410);

    auto g_zzzzzz_xyyy_1 = pbuffer.data(idx_eri_1_ig + 411);

    auto g_zzzzzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 412);

    auto g_zzzzzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 413);

    auto g_zzzzzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 414);

    auto g_zzzzzz_yyyy_1 = pbuffer.data(idx_eri_1_ig + 415);

    auto g_zzzzzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 416);

    auto g_zzzzzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 417);

    auto g_zzzzzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 418);

    auto g_zzzzzz_zzzz_1 = pbuffer.data(idx_eri_1_ig + 419);

    // Set up 0-15 components of targeted buffer : KG

    auto g_xxxxxxx_xxxx_0 = pbuffer.data(idx_eri_0_kg);

    auto g_xxxxxxx_xxxy_0 = pbuffer.data(idx_eri_0_kg + 1);

    auto g_xxxxxxx_xxxz_0 = pbuffer.data(idx_eri_0_kg + 2);

    auto g_xxxxxxx_xxyy_0 = pbuffer.data(idx_eri_0_kg + 3);

    auto g_xxxxxxx_xxyz_0 = pbuffer.data(idx_eri_0_kg + 4);

    auto g_xxxxxxx_xxzz_0 = pbuffer.data(idx_eri_0_kg + 5);

    auto g_xxxxxxx_xyyy_0 = pbuffer.data(idx_eri_0_kg + 6);

    auto g_xxxxxxx_xyyz_0 = pbuffer.data(idx_eri_0_kg + 7);

    auto g_xxxxxxx_xyzz_0 = pbuffer.data(idx_eri_0_kg + 8);

    auto g_xxxxxxx_xzzz_0 = pbuffer.data(idx_eri_0_kg + 9);

    auto g_xxxxxxx_yyyy_0 = pbuffer.data(idx_eri_0_kg + 10);

    auto g_xxxxxxx_yyyz_0 = pbuffer.data(idx_eri_0_kg + 11);

    auto g_xxxxxxx_yyzz_0 = pbuffer.data(idx_eri_0_kg + 12);

    auto g_xxxxxxx_yzzz_0 = pbuffer.data(idx_eri_0_kg + 13);

    auto g_xxxxxxx_zzzz_0 = pbuffer.data(idx_eri_0_kg + 14);

    #pragma omp simd aligned(g_xxxxx_xxxx_0, g_xxxxx_xxxx_1, g_xxxxx_xxxy_0, g_xxxxx_xxxy_1, g_xxxxx_xxxz_0, g_xxxxx_xxxz_1, g_xxxxx_xxyy_0, g_xxxxx_xxyy_1, g_xxxxx_xxyz_0, g_xxxxx_xxyz_1, g_xxxxx_xxzz_0, g_xxxxx_xxzz_1, g_xxxxx_xyyy_0, g_xxxxx_xyyy_1, g_xxxxx_xyyz_0, g_xxxxx_xyyz_1, g_xxxxx_xyzz_0, g_xxxxx_xyzz_1, g_xxxxx_xzzz_0, g_xxxxx_xzzz_1, g_xxxxx_yyyy_0, g_xxxxx_yyyy_1, g_xxxxx_yyyz_0, g_xxxxx_yyyz_1, g_xxxxx_yyzz_0, g_xxxxx_yyzz_1, g_xxxxx_yzzz_0, g_xxxxx_yzzz_1, g_xxxxx_zzzz_0, g_xxxxx_zzzz_1, g_xxxxxx_xxx_1, g_xxxxxx_xxxx_1, g_xxxxxx_xxxy_1, g_xxxxxx_xxxz_1, g_xxxxxx_xxy_1, g_xxxxxx_xxyy_1, g_xxxxxx_xxyz_1, g_xxxxxx_xxz_1, g_xxxxxx_xxzz_1, g_xxxxxx_xyy_1, g_xxxxxx_xyyy_1, g_xxxxxx_xyyz_1, g_xxxxxx_xyz_1, g_xxxxxx_xyzz_1, g_xxxxxx_xzz_1, g_xxxxxx_xzzz_1, g_xxxxxx_yyy_1, g_xxxxxx_yyyy_1, g_xxxxxx_yyyz_1, g_xxxxxx_yyz_1, g_xxxxxx_yyzz_1, g_xxxxxx_yzz_1, g_xxxxxx_yzzz_1, g_xxxxxx_zzz_1, g_xxxxxx_zzzz_1, g_xxxxxxx_xxxx_0, g_xxxxxxx_xxxy_0, g_xxxxxxx_xxxz_0, g_xxxxxxx_xxyy_0, g_xxxxxxx_xxyz_0, g_xxxxxxx_xxzz_0, g_xxxxxxx_xyyy_0, g_xxxxxxx_xyyz_0, g_xxxxxxx_xyzz_0, g_xxxxxxx_xzzz_0, g_xxxxxxx_yyyy_0, g_xxxxxxx_yyyz_0, g_xxxxxxx_yyzz_0, g_xxxxxxx_yzzz_0, g_xxxxxxx_zzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxxx_xxxx_0[i] = 6.0 * g_xxxxx_xxxx_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxx_1[i] * fz_be_0 + 4.0 * g_xxxxxx_xxx_1[i] * fe_0 + g_xxxxxx_xxxx_1[i] * pa_x[i];

        g_xxxxxxx_xxxy_0[i] = 6.0 * g_xxxxx_xxxy_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxy_1[i] * fz_be_0 + 3.0 * g_xxxxxx_xxy_1[i] * fe_0 + g_xxxxxx_xxxy_1[i] * pa_x[i];

        g_xxxxxxx_xxxz_0[i] = 6.0 * g_xxxxx_xxxz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_xxz_1[i] * fe_0 + g_xxxxxx_xxxz_1[i] * pa_x[i];

        g_xxxxxxx_xxyy_0[i] = 6.0 * g_xxxxx_xxyy_0[i] * fbe_0 - 6.0 * g_xxxxx_xxyy_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xyy_1[i] * fe_0 + g_xxxxxx_xxyy_1[i] * pa_x[i];

        g_xxxxxxx_xxyz_0[i] = 6.0 * g_xxxxx_xxyz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xyz_1[i] * fe_0 + g_xxxxxx_xxyz_1[i] * pa_x[i];

        g_xxxxxxx_xxzz_0[i] = 6.0 * g_xxxxx_xxzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xzz_1[i] * fe_0 + g_xxxxxx_xxzz_1[i] * pa_x[i];

        g_xxxxxxx_xyyy_0[i] = 6.0 * g_xxxxx_xyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_xyyy_1[i] * fz_be_0 + g_xxxxxx_yyy_1[i] * fe_0 + g_xxxxxx_xyyy_1[i] * pa_x[i];

        g_xxxxxxx_xyyz_0[i] = 6.0 * g_xxxxx_xyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_xyyz_1[i] * fz_be_0 + g_xxxxxx_yyz_1[i] * fe_0 + g_xxxxxx_xyyz_1[i] * pa_x[i];

        g_xxxxxxx_xyzz_0[i] = 6.0 * g_xxxxx_xyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xyzz_1[i] * fz_be_0 + g_xxxxxx_yzz_1[i] * fe_0 + g_xxxxxx_xyzz_1[i] * pa_x[i];

        g_xxxxxxx_xzzz_0[i] = 6.0 * g_xxxxx_xzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xzzz_1[i] * fz_be_0 + g_xxxxxx_zzz_1[i] * fe_0 + g_xxxxxx_xzzz_1[i] * pa_x[i];

        g_xxxxxxx_yyyy_0[i] = 6.0 * g_xxxxx_yyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_yyyy_1[i] * fz_be_0 + g_xxxxxx_yyyy_1[i] * pa_x[i];

        g_xxxxxxx_yyyz_0[i] = 6.0 * g_xxxxx_yyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_yyyz_1[i] * fz_be_0 + g_xxxxxx_yyyz_1[i] * pa_x[i];

        g_xxxxxxx_yyzz_0[i] = 6.0 * g_xxxxx_yyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_yyzz_1[i] * fz_be_0 + g_xxxxxx_yyzz_1[i] * pa_x[i];

        g_xxxxxxx_yzzz_0[i] = 6.0 * g_xxxxx_yzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_yzzz_1[i] * fz_be_0 + g_xxxxxx_yzzz_1[i] * pa_x[i];

        g_xxxxxxx_zzzz_0[i] = 6.0 * g_xxxxx_zzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_zzzz_1[i] * fz_be_0 + g_xxxxxx_zzzz_1[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : KG

    auto g_xxxxxxy_xxxx_0 = pbuffer.data(idx_eri_0_kg + 15);

    auto g_xxxxxxy_xxxy_0 = pbuffer.data(idx_eri_0_kg + 16);

    auto g_xxxxxxy_xxxz_0 = pbuffer.data(idx_eri_0_kg + 17);

    auto g_xxxxxxy_xxyy_0 = pbuffer.data(idx_eri_0_kg + 18);

    auto g_xxxxxxy_xxyz_0 = pbuffer.data(idx_eri_0_kg + 19);

    auto g_xxxxxxy_xxzz_0 = pbuffer.data(idx_eri_0_kg + 20);

    auto g_xxxxxxy_xyyy_0 = pbuffer.data(idx_eri_0_kg + 21);

    auto g_xxxxxxy_xyyz_0 = pbuffer.data(idx_eri_0_kg + 22);

    auto g_xxxxxxy_xyzz_0 = pbuffer.data(idx_eri_0_kg + 23);

    auto g_xxxxxxy_xzzz_0 = pbuffer.data(idx_eri_0_kg + 24);

    auto g_xxxxxxy_yyyy_0 = pbuffer.data(idx_eri_0_kg + 25);

    auto g_xxxxxxy_yyyz_0 = pbuffer.data(idx_eri_0_kg + 26);

    auto g_xxxxxxy_yyzz_0 = pbuffer.data(idx_eri_0_kg + 27);

    auto g_xxxxxxy_yzzz_0 = pbuffer.data(idx_eri_0_kg + 28);

    auto g_xxxxxxy_zzzz_0 = pbuffer.data(idx_eri_0_kg + 29);

    #pragma omp simd aligned(g_xxxxxx_xxx_1, g_xxxxxx_xxxx_1, g_xxxxxx_xxxy_1, g_xxxxxx_xxxz_1, g_xxxxxx_xxy_1, g_xxxxxx_xxyy_1, g_xxxxxx_xxyz_1, g_xxxxxx_xxz_1, g_xxxxxx_xxzz_1, g_xxxxxx_xyy_1, g_xxxxxx_xyyy_1, g_xxxxxx_xyyz_1, g_xxxxxx_xyz_1, g_xxxxxx_xyzz_1, g_xxxxxx_xzz_1, g_xxxxxx_xzzz_1, g_xxxxxx_yyy_1, g_xxxxxx_yyyy_1, g_xxxxxx_yyyz_1, g_xxxxxx_yyz_1, g_xxxxxx_yyzz_1, g_xxxxxx_yzz_1, g_xxxxxx_yzzz_1, g_xxxxxx_zzz_1, g_xxxxxx_zzzz_1, g_xxxxxxy_xxxx_0, g_xxxxxxy_xxxy_0, g_xxxxxxy_xxxz_0, g_xxxxxxy_xxyy_0, g_xxxxxxy_xxyz_0, g_xxxxxxy_xxzz_0, g_xxxxxxy_xyyy_0, g_xxxxxxy_xyyz_0, g_xxxxxxy_xyzz_0, g_xxxxxxy_xzzz_0, g_xxxxxxy_yyyy_0, g_xxxxxxy_yyyz_0, g_xxxxxxy_yyzz_0, g_xxxxxxy_yzzz_0, g_xxxxxxy_zzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxxy_xxxx_0[i] = g_xxxxxx_xxxx_1[i] * pa_y[i];

        g_xxxxxxy_xxxy_0[i] = g_xxxxxx_xxx_1[i] * fe_0 + g_xxxxxx_xxxy_1[i] * pa_y[i];

        g_xxxxxxy_xxxz_0[i] = g_xxxxxx_xxxz_1[i] * pa_y[i];

        g_xxxxxxy_xxyy_0[i] = 2.0 * g_xxxxxx_xxy_1[i] * fe_0 + g_xxxxxx_xxyy_1[i] * pa_y[i];

        g_xxxxxxy_xxyz_0[i] = g_xxxxxx_xxz_1[i] * fe_0 + g_xxxxxx_xxyz_1[i] * pa_y[i];

        g_xxxxxxy_xxzz_0[i] = g_xxxxxx_xxzz_1[i] * pa_y[i];

        g_xxxxxxy_xyyy_0[i] = 3.0 * g_xxxxxx_xyy_1[i] * fe_0 + g_xxxxxx_xyyy_1[i] * pa_y[i];

        g_xxxxxxy_xyyz_0[i] = 2.0 * g_xxxxxx_xyz_1[i] * fe_0 + g_xxxxxx_xyyz_1[i] * pa_y[i];

        g_xxxxxxy_xyzz_0[i] = g_xxxxxx_xzz_1[i] * fe_0 + g_xxxxxx_xyzz_1[i] * pa_y[i];

        g_xxxxxxy_xzzz_0[i] = g_xxxxxx_xzzz_1[i] * pa_y[i];

        g_xxxxxxy_yyyy_0[i] = 4.0 * g_xxxxxx_yyy_1[i] * fe_0 + g_xxxxxx_yyyy_1[i] * pa_y[i];

        g_xxxxxxy_yyyz_0[i] = 3.0 * g_xxxxxx_yyz_1[i] * fe_0 + g_xxxxxx_yyyz_1[i] * pa_y[i];

        g_xxxxxxy_yyzz_0[i] = 2.0 * g_xxxxxx_yzz_1[i] * fe_0 + g_xxxxxx_yyzz_1[i] * pa_y[i];

        g_xxxxxxy_yzzz_0[i] = g_xxxxxx_zzz_1[i] * fe_0 + g_xxxxxx_yzzz_1[i] * pa_y[i];

        g_xxxxxxy_zzzz_0[i] = g_xxxxxx_zzzz_1[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : KG

    auto g_xxxxxxz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 30);

    auto g_xxxxxxz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 31);

    auto g_xxxxxxz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 32);

    auto g_xxxxxxz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 33);

    auto g_xxxxxxz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 34);

    auto g_xxxxxxz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 35);

    auto g_xxxxxxz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 36);

    auto g_xxxxxxz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 37);

    auto g_xxxxxxz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 38);

    auto g_xxxxxxz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 39);

    auto g_xxxxxxz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 40);

    auto g_xxxxxxz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 41);

    auto g_xxxxxxz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 42);

    auto g_xxxxxxz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 43);

    auto g_xxxxxxz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 44);

    #pragma omp simd aligned(g_xxxxxx_xxx_1, g_xxxxxx_xxxx_1, g_xxxxxx_xxxy_1, g_xxxxxx_xxxz_1, g_xxxxxx_xxy_1, g_xxxxxx_xxyy_1, g_xxxxxx_xxyz_1, g_xxxxxx_xxz_1, g_xxxxxx_xxzz_1, g_xxxxxx_xyy_1, g_xxxxxx_xyyy_1, g_xxxxxx_xyyz_1, g_xxxxxx_xyz_1, g_xxxxxx_xyzz_1, g_xxxxxx_xzz_1, g_xxxxxx_xzzz_1, g_xxxxxx_yyy_1, g_xxxxxx_yyyy_1, g_xxxxxx_yyyz_1, g_xxxxxx_yyz_1, g_xxxxxx_yyzz_1, g_xxxxxx_yzz_1, g_xxxxxx_yzzz_1, g_xxxxxx_zzz_1, g_xxxxxx_zzzz_1, g_xxxxxxz_xxxx_0, g_xxxxxxz_xxxy_0, g_xxxxxxz_xxxz_0, g_xxxxxxz_xxyy_0, g_xxxxxxz_xxyz_0, g_xxxxxxz_xxzz_0, g_xxxxxxz_xyyy_0, g_xxxxxxz_xyyz_0, g_xxxxxxz_xyzz_0, g_xxxxxxz_xzzz_0, g_xxxxxxz_yyyy_0, g_xxxxxxz_yyyz_0, g_xxxxxxz_yyzz_0, g_xxxxxxz_yzzz_0, g_xxxxxxz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxxz_xxxx_0[i] = g_xxxxxx_xxxx_1[i] * pa_z[i];

        g_xxxxxxz_xxxy_0[i] = g_xxxxxx_xxxy_1[i] * pa_z[i];

        g_xxxxxxz_xxxz_0[i] = g_xxxxxx_xxx_1[i] * fe_0 + g_xxxxxx_xxxz_1[i] * pa_z[i];

        g_xxxxxxz_xxyy_0[i] = g_xxxxxx_xxyy_1[i] * pa_z[i];

        g_xxxxxxz_xxyz_0[i] = g_xxxxxx_xxy_1[i] * fe_0 + g_xxxxxx_xxyz_1[i] * pa_z[i];

        g_xxxxxxz_xxzz_0[i] = 2.0 * g_xxxxxx_xxz_1[i] * fe_0 + g_xxxxxx_xxzz_1[i] * pa_z[i];

        g_xxxxxxz_xyyy_0[i] = g_xxxxxx_xyyy_1[i] * pa_z[i];

        g_xxxxxxz_xyyz_0[i] = g_xxxxxx_xyy_1[i] * fe_0 + g_xxxxxx_xyyz_1[i] * pa_z[i];

        g_xxxxxxz_xyzz_0[i] = 2.0 * g_xxxxxx_xyz_1[i] * fe_0 + g_xxxxxx_xyzz_1[i] * pa_z[i];

        g_xxxxxxz_xzzz_0[i] = 3.0 * g_xxxxxx_xzz_1[i] * fe_0 + g_xxxxxx_xzzz_1[i] * pa_z[i];

        g_xxxxxxz_yyyy_0[i] = g_xxxxxx_yyyy_1[i] * pa_z[i];

        g_xxxxxxz_yyyz_0[i] = g_xxxxxx_yyy_1[i] * fe_0 + g_xxxxxx_yyyz_1[i] * pa_z[i];

        g_xxxxxxz_yyzz_0[i] = 2.0 * g_xxxxxx_yyz_1[i] * fe_0 + g_xxxxxx_yyzz_1[i] * pa_z[i];

        g_xxxxxxz_yzzz_0[i] = 3.0 * g_xxxxxx_yzz_1[i] * fe_0 + g_xxxxxx_yzzz_1[i] * pa_z[i];

        g_xxxxxxz_zzzz_0[i] = 4.0 * g_xxxxxx_zzz_1[i] * fe_0 + g_xxxxxx_zzzz_1[i] * pa_z[i];
    }

    // Set up 45-60 components of targeted buffer : KG

    auto g_xxxxxyy_xxxx_0 = pbuffer.data(idx_eri_0_kg + 45);

    auto g_xxxxxyy_xxxy_0 = pbuffer.data(idx_eri_0_kg + 46);

    auto g_xxxxxyy_xxxz_0 = pbuffer.data(idx_eri_0_kg + 47);

    auto g_xxxxxyy_xxyy_0 = pbuffer.data(idx_eri_0_kg + 48);

    auto g_xxxxxyy_xxyz_0 = pbuffer.data(idx_eri_0_kg + 49);

    auto g_xxxxxyy_xxzz_0 = pbuffer.data(idx_eri_0_kg + 50);

    auto g_xxxxxyy_xyyy_0 = pbuffer.data(idx_eri_0_kg + 51);

    auto g_xxxxxyy_xyyz_0 = pbuffer.data(idx_eri_0_kg + 52);

    auto g_xxxxxyy_xyzz_0 = pbuffer.data(idx_eri_0_kg + 53);

    auto g_xxxxxyy_xzzz_0 = pbuffer.data(idx_eri_0_kg + 54);

    auto g_xxxxxyy_yyyy_0 = pbuffer.data(idx_eri_0_kg + 55);

    auto g_xxxxxyy_yyyz_0 = pbuffer.data(idx_eri_0_kg + 56);

    auto g_xxxxxyy_yyzz_0 = pbuffer.data(idx_eri_0_kg + 57);

    auto g_xxxxxyy_yzzz_0 = pbuffer.data(idx_eri_0_kg + 58);

    auto g_xxxxxyy_zzzz_0 = pbuffer.data(idx_eri_0_kg + 59);

    #pragma omp simd aligned(g_xxxxx_xxxx_0, g_xxxxx_xxxx_1, g_xxxxx_xxxz_0, g_xxxxx_xxxz_1, g_xxxxx_xxzz_0, g_xxxxx_xxzz_1, g_xxxxx_xzzz_0, g_xxxxx_xzzz_1, g_xxxxxy_xxxx_1, g_xxxxxy_xxxz_1, g_xxxxxy_xxzz_1, g_xxxxxy_xzzz_1, g_xxxxxyy_xxxx_0, g_xxxxxyy_xxxy_0, g_xxxxxyy_xxxz_0, g_xxxxxyy_xxyy_0, g_xxxxxyy_xxyz_0, g_xxxxxyy_xxzz_0, g_xxxxxyy_xyyy_0, g_xxxxxyy_xyyz_0, g_xxxxxyy_xyzz_0, g_xxxxxyy_xzzz_0, g_xxxxxyy_yyyy_0, g_xxxxxyy_yyyz_0, g_xxxxxyy_yyzz_0, g_xxxxxyy_yzzz_0, g_xxxxxyy_zzzz_0, g_xxxxyy_xxxy_1, g_xxxxyy_xxy_1, g_xxxxyy_xxyy_1, g_xxxxyy_xxyz_1, g_xxxxyy_xyy_1, g_xxxxyy_xyyy_1, g_xxxxyy_xyyz_1, g_xxxxyy_xyz_1, g_xxxxyy_xyzz_1, g_xxxxyy_yyy_1, g_xxxxyy_yyyy_1, g_xxxxyy_yyyz_1, g_xxxxyy_yyz_1, g_xxxxyy_yyzz_1, g_xxxxyy_yzz_1, g_xxxxyy_yzzz_1, g_xxxxyy_zzzz_1, g_xxxyy_xxxy_0, g_xxxyy_xxxy_1, g_xxxyy_xxyy_0, g_xxxyy_xxyy_1, g_xxxyy_xxyz_0, g_xxxyy_xxyz_1, g_xxxyy_xyyy_0, g_xxxyy_xyyy_1, g_xxxyy_xyyz_0, g_xxxyy_xyyz_1, g_xxxyy_xyzz_0, g_xxxyy_xyzz_1, g_xxxyy_yyyy_0, g_xxxyy_yyyy_1, g_xxxyy_yyyz_0, g_xxxyy_yyyz_1, g_xxxyy_yyzz_0, g_xxxyy_yyzz_1, g_xxxyy_yzzz_0, g_xxxyy_yzzz_1, g_xxxyy_zzzz_0, g_xxxyy_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxyy_xxxx_0[i] = g_xxxxx_xxxx_0[i] * fbe_0 - g_xxxxx_xxxx_1[i] * fz_be_0 + g_xxxxxy_xxxx_1[i] * pa_y[i];

        g_xxxxxyy_xxxy_0[i] = 4.0 * g_xxxyy_xxxy_0[i] * fbe_0 - 4.0 * g_xxxyy_xxxy_1[i] * fz_be_0 + 3.0 * g_xxxxyy_xxy_1[i] * fe_0 + g_xxxxyy_xxxy_1[i] * pa_x[i];

        g_xxxxxyy_xxxz_0[i] = g_xxxxx_xxxz_0[i] * fbe_0 - g_xxxxx_xxxz_1[i] * fz_be_0 + g_xxxxxy_xxxz_1[i] * pa_y[i];

        g_xxxxxyy_xxyy_0[i] = 4.0 * g_xxxyy_xxyy_0[i] * fbe_0 - 4.0 * g_xxxyy_xxyy_1[i] * fz_be_0 + 2.0 * g_xxxxyy_xyy_1[i] * fe_0 + g_xxxxyy_xxyy_1[i] * pa_x[i];

        g_xxxxxyy_xxyz_0[i] = 4.0 * g_xxxyy_xxyz_0[i] * fbe_0 - 4.0 * g_xxxyy_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_xyz_1[i] * fe_0 + g_xxxxyy_xxyz_1[i] * pa_x[i];

        g_xxxxxyy_xxzz_0[i] = g_xxxxx_xxzz_0[i] * fbe_0 - g_xxxxx_xxzz_1[i] * fz_be_0 + g_xxxxxy_xxzz_1[i] * pa_y[i];

        g_xxxxxyy_xyyy_0[i] = 4.0 * g_xxxyy_xyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_xyyy_1[i] * fz_be_0 + g_xxxxyy_yyy_1[i] * fe_0 + g_xxxxyy_xyyy_1[i] * pa_x[i];

        g_xxxxxyy_xyyz_0[i] = 4.0 * g_xxxyy_xyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_xyyz_1[i] * fz_be_0 + g_xxxxyy_yyz_1[i] * fe_0 + g_xxxxyy_xyyz_1[i] * pa_x[i];

        g_xxxxxyy_xyzz_0[i] = 4.0 * g_xxxyy_xyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_xyzz_1[i] * fz_be_0 + g_xxxxyy_yzz_1[i] * fe_0 + g_xxxxyy_xyzz_1[i] * pa_x[i];

        g_xxxxxyy_xzzz_0[i] = g_xxxxx_xzzz_0[i] * fbe_0 - g_xxxxx_xzzz_1[i] * fz_be_0 + g_xxxxxy_xzzz_1[i] * pa_y[i];

        g_xxxxxyy_yyyy_0[i] = 4.0 * g_xxxyy_yyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_yyyy_1[i] * fz_be_0 + g_xxxxyy_yyyy_1[i] * pa_x[i];

        g_xxxxxyy_yyyz_0[i] = 4.0 * g_xxxyy_yyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_yyyz_1[i] * fz_be_0 + g_xxxxyy_yyyz_1[i] * pa_x[i];

        g_xxxxxyy_yyzz_0[i] = 4.0 * g_xxxyy_yyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_yyzz_1[i] * fz_be_0 + g_xxxxyy_yyzz_1[i] * pa_x[i];

        g_xxxxxyy_yzzz_0[i] = 4.0 * g_xxxyy_yzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_yzzz_1[i] * fz_be_0 + g_xxxxyy_yzzz_1[i] * pa_x[i];

        g_xxxxxyy_zzzz_0[i] = 4.0 * g_xxxyy_zzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_zzzz_1[i] * fz_be_0 + g_xxxxyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 60-75 components of targeted buffer : KG

    auto g_xxxxxyz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 60);

    auto g_xxxxxyz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 61);

    auto g_xxxxxyz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 62);

    auto g_xxxxxyz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 63);

    auto g_xxxxxyz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 64);

    auto g_xxxxxyz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 65);

    auto g_xxxxxyz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 66);

    auto g_xxxxxyz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 67);

    auto g_xxxxxyz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 68);

    auto g_xxxxxyz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 69);

    auto g_xxxxxyz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 70);

    auto g_xxxxxyz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 71);

    auto g_xxxxxyz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 72);

    auto g_xxxxxyz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 73);

    auto g_xxxxxyz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 74);

    #pragma omp simd aligned(g_xxxxxy_xxxy_1, g_xxxxxy_xxyy_1, g_xxxxxy_xyyy_1, g_xxxxxy_yyyy_1, g_xxxxxyz_xxxx_0, g_xxxxxyz_xxxy_0, g_xxxxxyz_xxxz_0, g_xxxxxyz_xxyy_0, g_xxxxxyz_xxyz_0, g_xxxxxyz_xxzz_0, g_xxxxxyz_xyyy_0, g_xxxxxyz_xyyz_0, g_xxxxxyz_xyzz_0, g_xxxxxyz_xzzz_0, g_xxxxxyz_yyyy_0, g_xxxxxyz_yyyz_0, g_xxxxxyz_yyzz_0, g_xxxxxyz_yzzz_0, g_xxxxxyz_zzzz_0, g_xxxxxz_xxxx_1, g_xxxxxz_xxxz_1, g_xxxxxz_xxyz_1, g_xxxxxz_xxz_1, g_xxxxxz_xxzz_1, g_xxxxxz_xyyz_1, g_xxxxxz_xyz_1, g_xxxxxz_xyzz_1, g_xxxxxz_xzz_1, g_xxxxxz_xzzz_1, g_xxxxxz_yyyz_1, g_xxxxxz_yyz_1, g_xxxxxz_yyzz_1, g_xxxxxz_yzz_1, g_xxxxxz_yzzz_1, g_xxxxxz_zzz_1, g_xxxxxz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxyz_xxxx_0[i] = g_xxxxxz_xxxx_1[i] * pa_y[i];

        g_xxxxxyz_xxxy_0[i] = g_xxxxxy_xxxy_1[i] * pa_z[i];

        g_xxxxxyz_xxxz_0[i] = g_xxxxxz_xxxz_1[i] * pa_y[i];

        g_xxxxxyz_xxyy_0[i] = g_xxxxxy_xxyy_1[i] * pa_z[i];

        g_xxxxxyz_xxyz_0[i] = g_xxxxxz_xxz_1[i] * fe_0 + g_xxxxxz_xxyz_1[i] * pa_y[i];

        g_xxxxxyz_xxzz_0[i] = g_xxxxxz_xxzz_1[i] * pa_y[i];

        g_xxxxxyz_xyyy_0[i] = g_xxxxxy_xyyy_1[i] * pa_z[i];

        g_xxxxxyz_xyyz_0[i] = 2.0 * g_xxxxxz_xyz_1[i] * fe_0 + g_xxxxxz_xyyz_1[i] * pa_y[i];

        g_xxxxxyz_xyzz_0[i] = g_xxxxxz_xzz_1[i] * fe_0 + g_xxxxxz_xyzz_1[i] * pa_y[i];

        g_xxxxxyz_xzzz_0[i] = g_xxxxxz_xzzz_1[i] * pa_y[i];

        g_xxxxxyz_yyyy_0[i] = g_xxxxxy_yyyy_1[i] * pa_z[i];

        g_xxxxxyz_yyyz_0[i] = 3.0 * g_xxxxxz_yyz_1[i] * fe_0 + g_xxxxxz_yyyz_1[i] * pa_y[i];

        g_xxxxxyz_yyzz_0[i] = 2.0 * g_xxxxxz_yzz_1[i] * fe_0 + g_xxxxxz_yyzz_1[i] * pa_y[i];

        g_xxxxxyz_yzzz_0[i] = g_xxxxxz_zzz_1[i] * fe_0 + g_xxxxxz_yzzz_1[i] * pa_y[i];

        g_xxxxxyz_zzzz_0[i] = g_xxxxxz_zzzz_1[i] * pa_y[i];
    }

    // Set up 75-90 components of targeted buffer : KG

    auto g_xxxxxzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 75);

    auto g_xxxxxzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 76);

    auto g_xxxxxzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 77);

    auto g_xxxxxzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 78);

    auto g_xxxxxzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 79);

    auto g_xxxxxzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 80);

    auto g_xxxxxzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 81);

    auto g_xxxxxzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 82);

    auto g_xxxxxzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 83);

    auto g_xxxxxzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 84);

    auto g_xxxxxzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 85);

    auto g_xxxxxzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 86);

    auto g_xxxxxzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 87);

    auto g_xxxxxzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 88);

    auto g_xxxxxzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 89);

    #pragma omp simd aligned(g_xxxxx_xxxx_0, g_xxxxx_xxxx_1, g_xxxxx_xxxy_0, g_xxxxx_xxxy_1, g_xxxxx_xxyy_0, g_xxxxx_xxyy_1, g_xxxxx_xyyy_0, g_xxxxx_xyyy_1, g_xxxxxz_xxxx_1, g_xxxxxz_xxxy_1, g_xxxxxz_xxyy_1, g_xxxxxz_xyyy_1, g_xxxxxzz_xxxx_0, g_xxxxxzz_xxxy_0, g_xxxxxzz_xxxz_0, g_xxxxxzz_xxyy_0, g_xxxxxzz_xxyz_0, g_xxxxxzz_xxzz_0, g_xxxxxzz_xyyy_0, g_xxxxxzz_xyyz_0, g_xxxxxzz_xyzz_0, g_xxxxxzz_xzzz_0, g_xxxxxzz_yyyy_0, g_xxxxxzz_yyyz_0, g_xxxxxzz_yyzz_0, g_xxxxxzz_yzzz_0, g_xxxxxzz_zzzz_0, g_xxxxzz_xxxz_1, g_xxxxzz_xxyz_1, g_xxxxzz_xxz_1, g_xxxxzz_xxzz_1, g_xxxxzz_xyyz_1, g_xxxxzz_xyz_1, g_xxxxzz_xyzz_1, g_xxxxzz_xzz_1, g_xxxxzz_xzzz_1, g_xxxxzz_yyyy_1, g_xxxxzz_yyyz_1, g_xxxxzz_yyz_1, g_xxxxzz_yyzz_1, g_xxxxzz_yzz_1, g_xxxxzz_yzzz_1, g_xxxxzz_zzz_1, g_xxxxzz_zzzz_1, g_xxxzz_xxxz_0, g_xxxzz_xxxz_1, g_xxxzz_xxyz_0, g_xxxzz_xxyz_1, g_xxxzz_xxzz_0, g_xxxzz_xxzz_1, g_xxxzz_xyyz_0, g_xxxzz_xyyz_1, g_xxxzz_xyzz_0, g_xxxzz_xyzz_1, g_xxxzz_xzzz_0, g_xxxzz_xzzz_1, g_xxxzz_yyyy_0, g_xxxzz_yyyy_1, g_xxxzz_yyyz_0, g_xxxzz_yyyz_1, g_xxxzz_yyzz_0, g_xxxzz_yyzz_1, g_xxxzz_yzzz_0, g_xxxzz_yzzz_1, g_xxxzz_zzzz_0, g_xxxzz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxzz_xxxx_0[i] = g_xxxxx_xxxx_0[i] * fbe_0 - g_xxxxx_xxxx_1[i] * fz_be_0 + g_xxxxxz_xxxx_1[i] * pa_z[i];

        g_xxxxxzz_xxxy_0[i] = g_xxxxx_xxxy_0[i] * fbe_0 - g_xxxxx_xxxy_1[i] * fz_be_0 + g_xxxxxz_xxxy_1[i] * pa_z[i];

        g_xxxxxzz_xxxz_0[i] = 4.0 * g_xxxzz_xxxz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxxz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_xxz_1[i] * fe_0 + g_xxxxzz_xxxz_1[i] * pa_x[i];

        g_xxxxxzz_xxyy_0[i] = g_xxxxx_xxyy_0[i] * fbe_0 - g_xxxxx_xxyy_1[i] * fz_be_0 + g_xxxxxz_xxyy_1[i] * pa_z[i];

        g_xxxxxzz_xxyz_0[i] = 4.0 * g_xxxzz_xxyz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_xyz_1[i] * fe_0 + g_xxxxzz_xxyz_1[i] * pa_x[i];

        g_xxxxxzz_xxzz_0[i] = 4.0 * g_xxxzz_xxzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_xzz_1[i] * fe_0 + g_xxxxzz_xxzz_1[i] * pa_x[i];

        g_xxxxxzz_xyyy_0[i] = g_xxxxx_xyyy_0[i] * fbe_0 - g_xxxxx_xyyy_1[i] * fz_be_0 + g_xxxxxz_xyyy_1[i] * pa_z[i];

        g_xxxxxzz_xyyz_0[i] = 4.0 * g_xxxzz_xyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_xyyz_1[i] * fz_be_0 + g_xxxxzz_yyz_1[i] * fe_0 + g_xxxxzz_xyyz_1[i] * pa_x[i];

        g_xxxxxzz_xyzz_0[i] = 4.0 * g_xxxzz_xyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xyzz_1[i] * fz_be_0 + g_xxxxzz_yzz_1[i] * fe_0 + g_xxxxzz_xyzz_1[i] * pa_x[i];

        g_xxxxxzz_xzzz_0[i] = 4.0 * g_xxxzz_xzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xzzz_1[i] * fz_be_0 + g_xxxxzz_zzz_1[i] * fe_0 + g_xxxxzz_xzzz_1[i] * pa_x[i];

        g_xxxxxzz_yyyy_0[i] = 4.0 * g_xxxzz_yyyy_0[i] * fbe_0 - 4.0 * g_xxxzz_yyyy_1[i] * fz_be_0 + g_xxxxzz_yyyy_1[i] * pa_x[i];

        g_xxxxxzz_yyyz_0[i] = 4.0 * g_xxxzz_yyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_yyyz_1[i] * fz_be_0 + g_xxxxzz_yyyz_1[i] * pa_x[i];

        g_xxxxxzz_yyzz_0[i] = 4.0 * g_xxxzz_yyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_yyzz_1[i] * fz_be_0 + g_xxxxzz_yyzz_1[i] * pa_x[i];

        g_xxxxxzz_yzzz_0[i] = 4.0 * g_xxxzz_yzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_yzzz_1[i] * fz_be_0 + g_xxxxzz_yzzz_1[i] * pa_x[i];

        g_xxxxxzz_zzzz_0[i] = 4.0 * g_xxxzz_zzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_zzzz_1[i] * fz_be_0 + g_xxxxzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 90-105 components of targeted buffer : KG

    auto g_xxxxyyy_xxxx_0 = pbuffer.data(idx_eri_0_kg + 90);

    auto g_xxxxyyy_xxxy_0 = pbuffer.data(idx_eri_0_kg + 91);

    auto g_xxxxyyy_xxxz_0 = pbuffer.data(idx_eri_0_kg + 92);

    auto g_xxxxyyy_xxyy_0 = pbuffer.data(idx_eri_0_kg + 93);

    auto g_xxxxyyy_xxyz_0 = pbuffer.data(idx_eri_0_kg + 94);

    auto g_xxxxyyy_xxzz_0 = pbuffer.data(idx_eri_0_kg + 95);

    auto g_xxxxyyy_xyyy_0 = pbuffer.data(idx_eri_0_kg + 96);

    auto g_xxxxyyy_xyyz_0 = pbuffer.data(idx_eri_0_kg + 97);

    auto g_xxxxyyy_xyzz_0 = pbuffer.data(idx_eri_0_kg + 98);

    auto g_xxxxyyy_xzzz_0 = pbuffer.data(idx_eri_0_kg + 99);

    auto g_xxxxyyy_yyyy_0 = pbuffer.data(idx_eri_0_kg + 100);

    auto g_xxxxyyy_yyyz_0 = pbuffer.data(idx_eri_0_kg + 101);

    auto g_xxxxyyy_yyzz_0 = pbuffer.data(idx_eri_0_kg + 102);

    auto g_xxxxyyy_yzzz_0 = pbuffer.data(idx_eri_0_kg + 103);

    auto g_xxxxyyy_zzzz_0 = pbuffer.data(idx_eri_0_kg + 104);

    #pragma omp simd aligned(g_xxxxy_xxxx_0, g_xxxxy_xxxx_1, g_xxxxy_xxxz_0, g_xxxxy_xxxz_1, g_xxxxy_xxzz_0, g_xxxxy_xxzz_1, g_xxxxy_xzzz_0, g_xxxxy_xzzz_1, g_xxxxyy_xxxx_1, g_xxxxyy_xxxz_1, g_xxxxyy_xxzz_1, g_xxxxyy_xzzz_1, g_xxxxyyy_xxxx_0, g_xxxxyyy_xxxy_0, g_xxxxyyy_xxxz_0, g_xxxxyyy_xxyy_0, g_xxxxyyy_xxyz_0, g_xxxxyyy_xxzz_0, g_xxxxyyy_xyyy_0, g_xxxxyyy_xyyz_0, g_xxxxyyy_xyzz_0, g_xxxxyyy_xzzz_0, g_xxxxyyy_yyyy_0, g_xxxxyyy_yyyz_0, g_xxxxyyy_yyzz_0, g_xxxxyyy_yzzz_0, g_xxxxyyy_zzzz_0, g_xxxyyy_xxxy_1, g_xxxyyy_xxy_1, g_xxxyyy_xxyy_1, g_xxxyyy_xxyz_1, g_xxxyyy_xyy_1, g_xxxyyy_xyyy_1, g_xxxyyy_xyyz_1, g_xxxyyy_xyz_1, g_xxxyyy_xyzz_1, g_xxxyyy_yyy_1, g_xxxyyy_yyyy_1, g_xxxyyy_yyyz_1, g_xxxyyy_yyz_1, g_xxxyyy_yyzz_1, g_xxxyyy_yzz_1, g_xxxyyy_yzzz_1, g_xxxyyy_zzzz_1, g_xxyyy_xxxy_0, g_xxyyy_xxxy_1, g_xxyyy_xxyy_0, g_xxyyy_xxyy_1, g_xxyyy_xxyz_0, g_xxyyy_xxyz_1, g_xxyyy_xyyy_0, g_xxyyy_xyyy_1, g_xxyyy_xyyz_0, g_xxyyy_xyyz_1, g_xxyyy_xyzz_0, g_xxyyy_xyzz_1, g_xxyyy_yyyy_0, g_xxyyy_yyyy_1, g_xxyyy_yyyz_0, g_xxyyy_yyyz_1, g_xxyyy_yyzz_0, g_xxyyy_yyzz_1, g_xxyyy_yzzz_0, g_xxyyy_yzzz_1, g_xxyyy_zzzz_0, g_xxyyy_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxyyy_xxxx_0[i] = 2.0 * g_xxxxy_xxxx_0[i] * fbe_0 - 2.0 * g_xxxxy_xxxx_1[i] * fz_be_0 + g_xxxxyy_xxxx_1[i] * pa_y[i];

        g_xxxxyyy_xxxy_0[i] = 3.0 * g_xxyyy_xxxy_0[i] * fbe_0 - 3.0 * g_xxyyy_xxxy_1[i] * fz_be_0 + 3.0 * g_xxxyyy_xxy_1[i] * fe_0 + g_xxxyyy_xxxy_1[i] * pa_x[i];

        g_xxxxyyy_xxxz_0[i] = 2.0 * g_xxxxy_xxxz_0[i] * fbe_0 - 2.0 * g_xxxxy_xxxz_1[i] * fz_be_0 + g_xxxxyy_xxxz_1[i] * pa_y[i];

        g_xxxxyyy_xxyy_0[i] = 3.0 * g_xxyyy_xxyy_0[i] * fbe_0 - 3.0 * g_xxyyy_xxyy_1[i] * fz_be_0 + 2.0 * g_xxxyyy_xyy_1[i] * fe_0 + g_xxxyyy_xxyy_1[i] * pa_x[i];

        g_xxxxyyy_xxyz_0[i] = 3.0 * g_xxyyy_xxyz_0[i] * fbe_0 - 3.0 * g_xxyyy_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_xyz_1[i] * fe_0 + g_xxxyyy_xxyz_1[i] * pa_x[i];

        g_xxxxyyy_xxzz_0[i] = 2.0 * g_xxxxy_xxzz_0[i] * fbe_0 - 2.0 * g_xxxxy_xxzz_1[i] * fz_be_0 + g_xxxxyy_xxzz_1[i] * pa_y[i];

        g_xxxxyyy_xyyy_0[i] = 3.0 * g_xxyyy_xyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_xyyy_1[i] * fz_be_0 + g_xxxyyy_yyy_1[i] * fe_0 + g_xxxyyy_xyyy_1[i] * pa_x[i];

        g_xxxxyyy_xyyz_0[i] = 3.0 * g_xxyyy_xyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_xyyz_1[i] * fz_be_0 + g_xxxyyy_yyz_1[i] * fe_0 + g_xxxyyy_xyyz_1[i] * pa_x[i];

        g_xxxxyyy_xyzz_0[i] = 3.0 * g_xxyyy_xyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_xyzz_1[i] * fz_be_0 + g_xxxyyy_yzz_1[i] * fe_0 + g_xxxyyy_xyzz_1[i] * pa_x[i];

        g_xxxxyyy_xzzz_0[i] = 2.0 * g_xxxxy_xzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_xzzz_1[i] * fz_be_0 + g_xxxxyy_xzzz_1[i] * pa_y[i];

        g_xxxxyyy_yyyy_0[i] = 3.0 * g_xxyyy_yyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_yyyy_1[i] * fz_be_0 + g_xxxyyy_yyyy_1[i] * pa_x[i];

        g_xxxxyyy_yyyz_0[i] = 3.0 * g_xxyyy_yyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_yyyz_1[i] * fz_be_0 + g_xxxyyy_yyyz_1[i] * pa_x[i];

        g_xxxxyyy_yyzz_0[i] = 3.0 * g_xxyyy_yyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_yyzz_1[i] * fz_be_0 + g_xxxyyy_yyzz_1[i] * pa_x[i];

        g_xxxxyyy_yzzz_0[i] = 3.0 * g_xxyyy_yzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_yzzz_1[i] * fz_be_0 + g_xxxyyy_yzzz_1[i] * pa_x[i];

        g_xxxxyyy_zzzz_0[i] = 3.0 * g_xxyyy_zzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_zzzz_1[i] * fz_be_0 + g_xxxyyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 105-120 components of targeted buffer : KG

    auto g_xxxxyyz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 105);

    auto g_xxxxyyz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 106);

    auto g_xxxxyyz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 107);

    auto g_xxxxyyz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 108);

    auto g_xxxxyyz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 109);

    auto g_xxxxyyz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 110);

    auto g_xxxxyyz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 111);

    auto g_xxxxyyz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 112);

    auto g_xxxxyyz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 113);

    auto g_xxxxyyz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 114);

    auto g_xxxxyyz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 115);

    auto g_xxxxyyz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 116);

    auto g_xxxxyyz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 117);

    auto g_xxxxyyz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 118);

    auto g_xxxxyyz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 119);

    #pragma omp simd aligned(g_xxxxyy_xxx_1, g_xxxxyy_xxxx_1, g_xxxxyy_xxxy_1, g_xxxxyy_xxxz_1, g_xxxxyy_xxy_1, g_xxxxyy_xxyy_1, g_xxxxyy_xxyz_1, g_xxxxyy_xxz_1, g_xxxxyy_xxzz_1, g_xxxxyy_xyy_1, g_xxxxyy_xyyy_1, g_xxxxyy_xyyz_1, g_xxxxyy_xyz_1, g_xxxxyy_xyzz_1, g_xxxxyy_xzz_1, g_xxxxyy_xzzz_1, g_xxxxyy_yyy_1, g_xxxxyy_yyyy_1, g_xxxxyy_yyyz_1, g_xxxxyy_yyz_1, g_xxxxyy_yyzz_1, g_xxxxyy_yzz_1, g_xxxxyy_yzzz_1, g_xxxxyy_zzz_1, g_xxxxyy_zzzz_1, g_xxxxyyz_xxxx_0, g_xxxxyyz_xxxy_0, g_xxxxyyz_xxxz_0, g_xxxxyyz_xxyy_0, g_xxxxyyz_xxyz_0, g_xxxxyyz_xxzz_0, g_xxxxyyz_xyyy_0, g_xxxxyyz_xyyz_0, g_xxxxyyz_xyzz_0, g_xxxxyyz_xzzz_0, g_xxxxyyz_yyyy_0, g_xxxxyyz_yyyz_0, g_xxxxyyz_yyzz_0, g_xxxxyyz_yzzz_0, g_xxxxyyz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyyz_xxxx_0[i] = g_xxxxyy_xxxx_1[i] * pa_z[i];

        g_xxxxyyz_xxxy_0[i] = g_xxxxyy_xxxy_1[i] * pa_z[i];

        g_xxxxyyz_xxxz_0[i] = g_xxxxyy_xxx_1[i] * fe_0 + g_xxxxyy_xxxz_1[i] * pa_z[i];

        g_xxxxyyz_xxyy_0[i] = g_xxxxyy_xxyy_1[i] * pa_z[i];

        g_xxxxyyz_xxyz_0[i] = g_xxxxyy_xxy_1[i] * fe_0 + g_xxxxyy_xxyz_1[i] * pa_z[i];

        g_xxxxyyz_xxzz_0[i] = 2.0 * g_xxxxyy_xxz_1[i] * fe_0 + g_xxxxyy_xxzz_1[i] * pa_z[i];

        g_xxxxyyz_xyyy_0[i] = g_xxxxyy_xyyy_1[i] * pa_z[i];

        g_xxxxyyz_xyyz_0[i] = g_xxxxyy_xyy_1[i] * fe_0 + g_xxxxyy_xyyz_1[i] * pa_z[i];

        g_xxxxyyz_xyzz_0[i] = 2.0 * g_xxxxyy_xyz_1[i] * fe_0 + g_xxxxyy_xyzz_1[i] * pa_z[i];

        g_xxxxyyz_xzzz_0[i] = 3.0 * g_xxxxyy_xzz_1[i] * fe_0 + g_xxxxyy_xzzz_1[i] * pa_z[i];

        g_xxxxyyz_yyyy_0[i] = g_xxxxyy_yyyy_1[i] * pa_z[i];

        g_xxxxyyz_yyyz_0[i] = g_xxxxyy_yyy_1[i] * fe_0 + g_xxxxyy_yyyz_1[i] * pa_z[i];

        g_xxxxyyz_yyzz_0[i] = 2.0 * g_xxxxyy_yyz_1[i] * fe_0 + g_xxxxyy_yyzz_1[i] * pa_z[i];

        g_xxxxyyz_yzzz_0[i] = 3.0 * g_xxxxyy_yzz_1[i] * fe_0 + g_xxxxyy_yzzz_1[i] * pa_z[i];

        g_xxxxyyz_zzzz_0[i] = 4.0 * g_xxxxyy_zzz_1[i] * fe_0 + g_xxxxyy_zzzz_1[i] * pa_z[i];
    }

    // Set up 120-135 components of targeted buffer : KG

    auto g_xxxxyzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 120);

    auto g_xxxxyzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 121);

    auto g_xxxxyzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 122);

    auto g_xxxxyzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 123);

    auto g_xxxxyzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 124);

    auto g_xxxxyzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 125);

    auto g_xxxxyzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 126);

    auto g_xxxxyzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 127);

    auto g_xxxxyzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 128);

    auto g_xxxxyzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 129);

    auto g_xxxxyzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 130);

    auto g_xxxxyzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 131);

    auto g_xxxxyzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 132);

    auto g_xxxxyzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 133);

    auto g_xxxxyzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 134);

    #pragma omp simd aligned(g_xxxxyzz_xxxx_0, g_xxxxyzz_xxxy_0, g_xxxxyzz_xxxz_0, g_xxxxyzz_xxyy_0, g_xxxxyzz_xxyz_0, g_xxxxyzz_xxzz_0, g_xxxxyzz_xyyy_0, g_xxxxyzz_xyyz_0, g_xxxxyzz_xyzz_0, g_xxxxyzz_xzzz_0, g_xxxxyzz_yyyy_0, g_xxxxyzz_yyyz_0, g_xxxxyzz_yyzz_0, g_xxxxyzz_yzzz_0, g_xxxxyzz_zzzz_0, g_xxxxzz_xxx_1, g_xxxxzz_xxxx_1, g_xxxxzz_xxxy_1, g_xxxxzz_xxxz_1, g_xxxxzz_xxy_1, g_xxxxzz_xxyy_1, g_xxxxzz_xxyz_1, g_xxxxzz_xxz_1, g_xxxxzz_xxzz_1, g_xxxxzz_xyy_1, g_xxxxzz_xyyy_1, g_xxxxzz_xyyz_1, g_xxxxzz_xyz_1, g_xxxxzz_xyzz_1, g_xxxxzz_xzz_1, g_xxxxzz_xzzz_1, g_xxxxzz_yyy_1, g_xxxxzz_yyyy_1, g_xxxxzz_yyyz_1, g_xxxxzz_yyz_1, g_xxxxzz_yyzz_1, g_xxxxzz_yzz_1, g_xxxxzz_yzzz_1, g_xxxxzz_zzz_1, g_xxxxzz_zzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyzz_xxxx_0[i] = g_xxxxzz_xxxx_1[i] * pa_y[i];

        g_xxxxyzz_xxxy_0[i] = g_xxxxzz_xxx_1[i] * fe_0 + g_xxxxzz_xxxy_1[i] * pa_y[i];

        g_xxxxyzz_xxxz_0[i] = g_xxxxzz_xxxz_1[i] * pa_y[i];

        g_xxxxyzz_xxyy_0[i] = 2.0 * g_xxxxzz_xxy_1[i] * fe_0 + g_xxxxzz_xxyy_1[i] * pa_y[i];

        g_xxxxyzz_xxyz_0[i] = g_xxxxzz_xxz_1[i] * fe_0 + g_xxxxzz_xxyz_1[i] * pa_y[i];

        g_xxxxyzz_xxzz_0[i] = g_xxxxzz_xxzz_1[i] * pa_y[i];

        g_xxxxyzz_xyyy_0[i] = 3.0 * g_xxxxzz_xyy_1[i] * fe_0 + g_xxxxzz_xyyy_1[i] * pa_y[i];

        g_xxxxyzz_xyyz_0[i] = 2.0 * g_xxxxzz_xyz_1[i] * fe_0 + g_xxxxzz_xyyz_1[i] * pa_y[i];

        g_xxxxyzz_xyzz_0[i] = g_xxxxzz_xzz_1[i] * fe_0 + g_xxxxzz_xyzz_1[i] * pa_y[i];

        g_xxxxyzz_xzzz_0[i] = g_xxxxzz_xzzz_1[i] * pa_y[i];

        g_xxxxyzz_yyyy_0[i] = 4.0 * g_xxxxzz_yyy_1[i] * fe_0 + g_xxxxzz_yyyy_1[i] * pa_y[i];

        g_xxxxyzz_yyyz_0[i] = 3.0 * g_xxxxzz_yyz_1[i] * fe_0 + g_xxxxzz_yyyz_1[i] * pa_y[i];

        g_xxxxyzz_yyzz_0[i] = 2.0 * g_xxxxzz_yzz_1[i] * fe_0 + g_xxxxzz_yyzz_1[i] * pa_y[i];

        g_xxxxyzz_yzzz_0[i] = g_xxxxzz_zzz_1[i] * fe_0 + g_xxxxzz_yzzz_1[i] * pa_y[i];

        g_xxxxyzz_zzzz_0[i] = g_xxxxzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 135-150 components of targeted buffer : KG

    auto g_xxxxzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 135);

    auto g_xxxxzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 136);

    auto g_xxxxzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 137);

    auto g_xxxxzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 138);

    auto g_xxxxzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 139);

    auto g_xxxxzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 140);

    auto g_xxxxzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 141);

    auto g_xxxxzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 142);

    auto g_xxxxzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 143);

    auto g_xxxxzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 144);

    auto g_xxxxzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 145);

    auto g_xxxxzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 146);

    auto g_xxxxzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 147);

    auto g_xxxxzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 148);

    auto g_xxxxzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 149);

    #pragma omp simd aligned(g_xxxxz_xxxx_0, g_xxxxz_xxxx_1, g_xxxxz_xxxy_0, g_xxxxz_xxxy_1, g_xxxxz_xxyy_0, g_xxxxz_xxyy_1, g_xxxxz_xyyy_0, g_xxxxz_xyyy_1, g_xxxxzz_xxxx_1, g_xxxxzz_xxxy_1, g_xxxxzz_xxyy_1, g_xxxxzz_xyyy_1, g_xxxxzzz_xxxx_0, g_xxxxzzz_xxxy_0, g_xxxxzzz_xxxz_0, g_xxxxzzz_xxyy_0, g_xxxxzzz_xxyz_0, g_xxxxzzz_xxzz_0, g_xxxxzzz_xyyy_0, g_xxxxzzz_xyyz_0, g_xxxxzzz_xyzz_0, g_xxxxzzz_xzzz_0, g_xxxxzzz_yyyy_0, g_xxxxzzz_yyyz_0, g_xxxxzzz_yyzz_0, g_xxxxzzz_yzzz_0, g_xxxxzzz_zzzz_0, g_xxxzzz_xxxz_1, g_xxxzzz_xxyz_1, g_xxxzzz_xxz_1, g_xxxzzz_xxzz_1, g_xxxzzz_xyyz_1, g_xxxzzz_xyz_1, g_xxxzzz_xyzz_1, g_xxxzzz_xzz_1, g_xxxzzz_xzzz_1, g_xxxzzz_yyyy_1, g_xxxzzz_yyyz_1, g_xxxzzz_yyz_1, g_xxxzzz_yyzz_1, g_xxxzzz_yzz_1, g_xxxzzz_yzzz_1, g_xxxzzz_zzz_1, g_xxxzzz_zzzz_1, g_xxzzz_xxxz_0, g_xxzzz_xxxz_1, g_xxzzz_xxyz_0, g_xxzzz_xxyz_1, g_xxzzz_xxzz_0, g_xxzzz_xxzz_1, g_xxzzz_xyyz_0, g_xxzzz_xyyz_1, g_xxzzz_xyzz_0, g_xxzzz_xyzz_1, g_xxzzz_xzzz_0, g_xxzzz_xzzz_1, g_xxzzz_yyyy_0, g_xxzzz_yyyy_1, g_xxzzz_yyyz_0, g_xxzzz_yyyz_1, g_xxzzz_yyzz_0, g_xxzzz_yyzz_1, g_xxzzz_yzzz_0, g_xxzzz_yzzz_1, g_xxzzz_zzzz_0, g_xxzzz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxzzz_xxxx_0[i] = 2.0 * g_xxxxz_xxxx_0[i] * fbe_0 - 2.0 * g_xxxxz_xxxx_1[i] * fz_be_0 + g_xxxxzz_xxxx_1[i] * pa_z[i];

        g_xxxxzzz_xxxy_0[i] = 2.0 * g_xxxxz_xxxy_0[i] * fbe_0 - 2.0 * g_xxxxz_xxxy_1[i] * fz_be_0 + g_xxxxzz_xxxy_1[i] * pa_z[i];

        g_xxxxzzz_xxxz_0[i] = 3.0 * g_xxzzz_xxxz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxxz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_xxz_1[i] * fe_0 + g_xxxzzz_xxxz_1[i] * pa_x[i];

        g_xxxxzzz_xxyy_0[i] = 2.0 * g_xxxxz_xxyy_0[i] * fbe_0 - 2.0 * g_xxxxz_xxyy_1[i] * fz_be_0 + g_xxxxzz_xxyy_1[i] * pa_z[i];

        g_xxxxzzz_xxyz_0[i] = 3.0 * g_xxzzz_xxyz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_xyz_1[i] * fe_0 + g_xxxzzz_xxyz_1[i] * pa_x[i];

        g_xxxxzzz_xxzz_0[i] = 3.0 * g_xxzzz_xxzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_xzz_1[i] * fe_0 + g_xxxzzz_xxzz_1[i] * pa_x[i];

        g_xxxxzzz_xyyy_0[i] = 2.0 * g_xxxxz_xyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_xyyy_1[i] * fz_be_0 + g_xxxxzz_xyyy_1[i] * pa_z[i];

        g_xxxxzzz_xyyz_0[i] = 3.0 * g_xxzzz_xyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_xyyz_1[i] * fz_be_0 + g_xxxzzz_yyz_1[i] * fe_0 + g_xxxzzz_xyyz_1[i] * pa_x[i];

        g_xxxxzzz_xyzz_0[i] = 3.0 * g_xxzzz_xyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xyzz_1[i] * fz_be_0 + g_xxxzzz_yzz_1[i] * fe_0 + g_xxxzzz_xyzz_1[i] * pa_x[i];

        g_xxxxzzz_xzzz_0[i] = 3.0 * g_xxzzz_xzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xzzz_1[i] * fz_be_0 + g_xxxzzz_zzz_1[i] * fe_0 + g_xxxzzz_xzzz_1[i] * pa_x[i];

        g_xxxxzzz_yyyy_0[i] = 3.0 * g_xxzzz_yyyy_0[i] * fbe_0 - 3.0 * g_xxzzz_yyyy_1[i] * fz_be_0 + g_xxxzzz_yyyy_1[i] * pa_x[i];

        g_xxxxzzz_yyyz_0[i] = 3.0 * g_xxzzz_yyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_yyyz_1[i] * fz_be_0 + g_xxxzzz_yyyz_1[i] * pa_x[i];

        g_xxxxzzz_yyzz_0[i] = 3.0 * g_xxzzz_yyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_yyzz_1[i] * fz_be_0 + g_xxxzzz_yyzz_1[i] * pa_x[i];

        g_xxxxzzz_yzzz_0[i] = 3.0 * g_xxzzz_yzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_yzzz_1[i] * fz_be_0 + g_xxxzzz_yzzz_1[i] * pa_x[i];

        g_xxxxzzz_zzzz_0[i] = 3.0 * g_xxzzz_zzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_zzzz_1[i] * fz_be_0 + g_xxxzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 150-165 components of targeted buffer : KG

    auto g_xxxyyyy_xxxx_0 = pbuffer.data(idx_eri_0_kg + 150);

    auto g_xxxyyyy_xxxy_0 = pbuffer.data(idx_eri_0_kg + 151);

    auto g_xxxyyyy_xxxz_0 = pbuffer.data(idx_eri_0_kg + 152);

    auto g_xxxyyyy_xxyy_0 = pbuffer.data(idx_eri_0_kg + 153);

    auto g_xxxyyyy_xxyz_0 = pbuffer.data(idx_eri_0_kg + 154);

    auto g_xxxyyyy_xxzz_0 = pbuffer.data(idx_eri_0_kg + 155);

    auto g_xxxyyyy_xyyy_0 = pbuffer.data(idx_eri_0_kg + 156);

    auto g_xxxyyyy_xyyz_0 = pbuffer.data(idx_eri_0_kg + 157);

    auto g_xxxyyyy_xyzz_0 = pbuffer.data(idx_eri_0_kg + 158);

    auto g_xxxyyyy_xzzz_0 = pbuffer.data(idx_eri_0_kg + 159);

    auto g_xxxyyyy_yyyy_0 = pbuffer.data(idx_eri_0_kg + 160);

    auto g_xxxyyyy_yyyz_0 = pbuffer.data(idx_eri_0_kg + 161);

    auto g_xxxyyyy_yyzz_0 = pbuffer.data(idx_eri_0_kg + 162);

    auto g_xxxyyyy_yzzz_0 = pbuffer.data(idx_eri_0_kg + 163);

    auto g_xxxyyyy_zzzz_0 = pbuffer.data(idx_eri_0_kg + 164);

    #pragma omp simd aligned(g_xxxyy_xxxx_0, g_xxxyy_xxxx_1, g_xxxyy_xxxz_0, g_xxxyy_xxxz_1, g_xxxyy_xxzz_0, g_xxxyy_xxzz_1, g_xxxyy_xzzz_0, g_xxxyy_xzzz_1, g_xxxyyy_xxxx_1, g_xxxyyy_xxxz_1, g_xxxyyy_xxzz_1, g_xxxyyy_xzzz_1, g_xxxyyyy_xxxx_0, g_xxxyyyy_xxxy_0, g_xxxyyyy_xxxz_0, g_xxxyyyy_xxyy_0, g_xxxyyyy_xxyz_0, g_xxxyyyy_xxzz_0, g_xxxyyyy_xyyy_0, g_xxxyyyy_xyyz_0, g_xxxyyyy_xyzz_0, g_xxxyyyy_xzzz_0, g_xxxyyyy_yyyy_0, g_xxxyyyy_yyyz_0, g_xxxyyyy_yyzz_0, g_xxxyyyy_yzzz_0, g_xxxyyyy_zzzz_0, g_xxyyyy_xxxy_1, g_xxyyyy_xxy_1, g_xxyyyy_xxyy_1, g_xxyyyy_xxyz_1, g_xxyyyy_xyy_1, g_xxyyyy_xyyy_1, g_xxyyyy_xyyz_1, g_xxyyyy_xyz_1, g_xxyyyy_xyzz_1, g_xxyyyy_yyy_1, g_xxyyyy_yyyy_1, g_xxyyyy_yyyz_1, g_xxyyyy_yyz_1, g_xxyyyy_yyzz_1, g_xxyyyy_yzz_1, g_xxyyyy_yzzz_1, g_xxyyyy_zzzz_1, g_xyyyy_xxxy_0, g_xyyyy_xxxy_1, g_xyyyy_xxyy_0, g_xyyyy_xxyy_1, g_xyyyy_xxyz_0, g_xyyyy_xxyz_1, g_xyyyy_xyyy_0, g_xyyyy_xyyy_1, g_xyyyy_xyyz_0, g_xyyyy_xyyz_1, g_xyyyy_xyzz_0, g_xyyyy_xyzz_1, g_xyyyy_yyyy_0, g_xyyyy_yyyy_1, g_xyyyy_yyyz_0, g_xyyyy_yyyz_1, g_xyyyy_yyzz_0, g_xyyyy_yyzz_1, g_xyyyy_yzzz_0, g_xyyyy_yzzz_1, g_xyyyy_zzzz_0, g_xyyyy_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyyy_xxxx_0[i] = 3.0 * g_xxxyy_xxxx_0[i] * fbe_0 - 3.0 * g_xxxyy_xxxx_1[i] * fz_be_0 + g_xxxyyy_xxxx_1[i] * pa_y[i];

        g_xxxyyyy_xxxy_0[i] = 2.0 * g_xyyyy_xxxy_0[i] * fbe_0 - 2.0 * g_xyyyy_xxxy_1[i] * fz_be_0 + 3.0 * g_xxyyyy_xxy_1[i] * fe_0 + g_xxyyyy_xxxy_1[i] * pa_x[i];

        g_xxxyyyy_xxxz_0[i] = 3.0 * g_xxxyy_xxxz_0[i] * fbe_0 - 3.0 * g_xxxyy_xxxz_1[i] * fz_be_0 + g_xxxyyy_xxxz_1[i] * pa_y[i];

        g_xxxyyyy_xxyy_0[i] = 2.0 * g_xyyyy_xxyy_0[i] * fbe_0 - 2.0 * g_xyyyy_xxyy_1[i] * fz_be_0 + 2.0 * g_xxyyyy_xyy_1[i] * fe_0 + g_xxyyyy_xxyy_1[i] * pa_x[i];

        g_xxxyyyy_xxyz_0[i] = 2.0 * g_xyyyy_xxyz_0[i] * fbe_0 - 2.0 * g_xyyyy_xxyz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_xyz_1[i] * fe_0 + g_xxyyyy_xxyz_1[i] * pa_x[i];

        g_xxxyyyy_xxzz_0[i] = 3.0 * g_xxxyy_xxzz_0[i] * fbe_0 - 3.0 * g_xxxyy_xxzz_1[i] * fz_be_0 + g_xxxyyy_xxzz_1[i] * pa_y[i];

        g_xxxyyyy_xyyy_0[i] = 2.0 * g_xyyyy_xyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_xyyy_1[i] * fz_be_0 + g_xxyyyy_yyy_1[i] * fe_0 + g_xxyyyy_xyyy_1[i] * pa_x[i];

        g_xxxyyyy_xyyz_0[i] = 2.0 * g_xyyyy_xyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_xyyz_1[i] * fz_be_0 + g_xxyyyy_yyz_1[i] * fe_0 + g_xxyyyy_xyyz_1[i] * pa_x[i];

        g_xxxyyyy_xyzz_0[i] = 2.0 * g_xyyyy_xyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_xyzz_1[i] * fz_be_0 + g_xxyyyy_yzz_1[i] * fe_0 + g_xxyyyy_xyzz_1[i] * pa_x[i];

        g_xxxyyyy_xzzz_0[i] = 3.0 * g_xxxyy_xzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_xzzz_1[i] * fz_be_0 + g_xxxyyy_xzzz_1[i] * pa_y[i];

        g_xxxyyyy_yyyy_0[i] = 2.0 * g_xyyyy_yyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_yyyy_1[i] * fz_be_0 + g_xxyyyy_yyyy_1[i] * pa_x[i];

        g_xxxyyyy_yyyz_0[i] = 2.0 * g_xyyyy_yyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_yyyz_1[i] * fz_be_0 + g_xxyyyy_yyyz_1[i] * pa_x[i];

        g_xxxyyyy_yyzz_0[i] = 2.0 * g_xyyyy_yyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_yyzz_1[i] * fz_be_0 + g_xxyyyy_yyzz_1[i] * pa_x[i];

        g_xxxyyyy_yzzz_0[i] = 2.0 * g_xyyyy_yzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_yzzz_1[i] * fz_be_0 + g_xxyyyy_yzzz_1[i] * pa_x[i];

        g_xxxyyyy_zzzz_0[i] = 2.0 * g_xyyyy_zzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_zzzz_1[i] * fz_be_0 + g_xxyyyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 165-180 components of targeted buffer : KG

    auto g_xxxyyyz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 165);

    auto g_xxxyyyz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 166);

    auto g_xxxyyyz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 167);

    auto g_xxxyyyz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 168);

    auto g_xxxyyyz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 169);

    auto g_xxxyyyz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 170);

    auto g_xxxyyyz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 171);

    auto g_xxxyyyz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 172);

    auto g_xxxyyyz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 173);

    auto g_xxxyyyz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 174);

    auto g_xxxyyyz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 175);

    auto g_xxxyyyz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 176);

    auto g_xxxyyyz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 177);

    auto g_xxxyyyz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 178);

    auto g_xxxyyyz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 179);

    #pragma omp simd aligned(g_xxxyyy_xxx_1, g_xxxyyy_xxxx_1, g_xxxyyy_xxxy_1, g_xxxyyy_xxxz_1, g_xxxyyy_xxy_1, g_xxxyyy_xxyy_1, g_xxxyyy_xxyz_1, g_xxxyyy_xxz_1, g_xxxyyy_xxzz_1, g_xxxyyy_xyy_1, g_xxxyyy_xyyy_1, g_xxxyyy_xyyz_1, g_xxxyyy_xyz_1, g_xxxyyy_xyzz_1, g_xxxyyy_xzz_1, g_xxxyyy_xzzz_1, g_xxxyyy_yyy_1, g_xxxyyy_yyyy_1, g_xxxyyy_yyyz_1, g_xxxyyy_yyz_1, g_xxxyyy_yyzz_1, g_xxxyyy_yzz_1, g_xxxyyy_yzzz_1, g_xxxyyy_zzz_1, g_xxxyyy_zzzz_1, g_xxxyyyz_xxxx_0, g_xxxyyyz_xxxy_0, g_xxxyyyz_xxxz_0, g_xxxyyyz_xxyy_0, g_xxxyyyz_xxyz_0, g_xxxyyyz_xxzz_0, g_xxxyyyz_xyyy_0, g_xxxyyyz_xyyz_0, g_xxxyyyz_xyzz_0, g_xxxyyyz_xzzz_0, g_xxxyyyz_yyyy_0, g_xxxyyyz_yyyz_0, g_xxxyyyz_yyzz_0, g_xxxyyyz_yzzz_0, g_xxxyyyz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyyyz_xxxx_0[i] = g_xxxyyy_xxxx_1[i] * pa_z[i];

        g_xxxyyyz_xxxy_0[i] = g_xxxyyy_xxxy_1[i] * pa_z[i];

        g_xxxyyyz_xxxz_0[i] = g_xxxyyy_xxx_1[i] * fe_0 + g_xxxyyy_xxxz_1[i] * pa_z[i];

        g_xxxyyyz_xxyy_0[i] = g_xxxyyy_xxyy_1[i] * pa_z[i];

        g_xxxyyyz_xxyz_0[i] = g_xxxyyy_xxy_1[i] * fe_0 + g_xxxyyy_xxyz_1[i] * pa_z[i];

        g_xxxyyyz_xxzz_0[i] = 2.0 * g_xxxyyy_xxz_1[i] * fe_0 + g_xxxyyy_xxzz_1[i] * pa_z[i];

        g_xxxyyyz_xyyy_0[i] = g_xxxyyy_xyyy_1[i] * pa_z[i];

        g_xxxyyyz_xyyz_0[i] = g_xxxyyy_xyy_1[i] * fe_0 + g_xxxyyy_xyyz_1[i] * pa_z[i];

        g_xxxyyyz_xyzz_0[i] = 2.0 * g_xxxyyy_xyz_1[i] * fe_0 + g_xxxyyy_xyzz_1[i] * pa_z[i];

        g_xxxyyyz_xzzz_0[i] = 3.0 * g_xxxyyy_xzz_1[i] * fe_0 + g_xxxyyy_xzzz_1[i] * pa_z[i];

        g_xxxyyyz_yyyy_0[i] = g_xxxyyy_yyyy_1[i] * pa_z[i];

        g_xxxyyyz_yyyz_0[i] = g_xxxyyy_yyy_1[i] * fe_0 + g_xxxyyy_yyyz_1[i] * pa_z[i];

        g_xxxyyyz_yyzz_0[i] = 2.0 * g_xxxyyy_yyz_1[i] * fe_0 + g_xxxyyy_yyzz_1[i] * pa_z[i];

        g_xxxyyyz_yzzz_0[i] = 3.0 * g_xxxyyy_yzz_1[i] * fe_0 + g_xxxyyy_yzzz_1[i] * pa_z[i];

        g_xxxyyyz_zzzz_0[i] = 4.0 * g_xxxyyy_zzz_1[i] * fe_0 + g_xxxyyy_zzzz_1[i] * pa_z[i];
    }

    // Set up 180-195 components of targeted buffer : KG

    auto g_xxxyyzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 180);

    auto g_xxxyyzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 181);

    auto g_xxxyyzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 182);

    auto g_xxxyyzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 183);

    auto g_xxxyyzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 184);

    auto g_xxxyyzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 185);

    auto g_xxxyyzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 186);

    auto g_xxxyyzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 187);

    auto g_xxxyyzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 188);

    auto g_xxxyyzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 189);

    auto g_xxxyyzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 190);

    auto g_xxxyyzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 191);

    auto g_xxxyyzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 192);

    auto g_xxxyyzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 193);

    auto g_xxxyyzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 194);

    #pragma omp simd aligned(g_xxxyy_xxxy_0, g_xxxyy_xxxy_1, g_xxxyy_xxyy_0, g_xxxyy_xxyy_1, g_xxxyy_xyyy_0, g_xxxyy_xyyy_1, g_xxxyyz_xxxy_1, g_xxxyyz_xxyy_1, g_xxxyyz_xyyy_1, g_xxxyyzz_xxxx_0, g_xxxyyzz_xxxy_0, g_xxxyyzz_xxxz_0, g_xxxyyzz_xxyy_0, g_xxxyyzz_xxyz_0, g_xxxyyzz_xxzz_0, g_xxxyyzz_xyyy_0, g_xxxyyzz_xyyz_0, g_xxxyyzz_xyzz_0, g_xxxyyzz_xzzz_0, g_xxxyyzz_yyyy_0, g_xxxyyzz_yyyz_0, g_xxxyyzz_yyzz_0, g_xxxyyzz_yzzz_0, g_xxxyyzz_zzzz_0, g_xxxyzz_xxxx_1, g_xxxyzz_xxxz_1, g_xxxyzz_xxzz_1, g_xxxyzz_xzzz_1, g_xxxzz_xxxx_0, g_xxxzz_xxxx_1, g_xxxzz_xxxz_0, g_xxxzz_xxxz_1, g_xxxzz_xxzz_0, g_xxxzz_xxzz_1, g_xxxzz_xzzz_0, g_xxxzz_xzzz_1, g_xxyyzz_xxyz_1, g_xxyyzz_xyyz_1, g_xxyyzz_xyz_1, g_xxyyzz_xyzz_1, g_xxyyzz_yyyy_1, g_xxyyzz_yyyz_1, g_xxyyzz_yyz_1, g_xxyyzz_yyzz_1, g_xxyyzz_yzz_1, g_xxyyzz_yzzz_1, g_xxyyzz_zzzz_1, g_xyyzz_xxyz_0, g_xyyzz_xxyz_1, g_xyyzz_xyyz_0, g_xyyzz_xyyz_1, g_xyyzz_xyzz_0, g_xyyzz_xyzz_1, g_xyyzz_yyyy_0, g_xyyzz_yyyy_1, g_xyyzz_yyyz_0, g_xyyzz_yyyz_1, g_xyyzz_yyzz_0, g_xyyzz_yyzz_1, g_xyyzz_yzzz_0, g_xyyzz_yzzz_1, g_xyyzz_zzzz_0, g_xyyzz_zzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyzz_xxxx_0[i] = g_xxxzz_xxxx_0[i] * fbe_0 - g_xxxzz_xxxx_1[i] * fz_be_0 + g_xxxyzz_xxxx_1[i] * pa_y[i];

        g_xxxyyzz_xxxy_0[i] = g_xxxyy_xxxy_0[i] * fbe_0 - g_xxxyy_xxxy_1[i] * fz_be_0 + g_xxxyyz_xxxy_1[i] * pa_z[i];

        g_xxxyyzz_xxxz_0[i] = g_xxxzz_xxxz_0[i] * fbe_0 - g_xxxzz_xxxz_1[i] * fz_be_0 + g_xxxyzz_xxxz_1[i] * pa_y[i];

        g_xxxyyzz_xxyy_0[i] = g_xxxyy_xxyy_0[i] * fbe_0 - g_xxxyy_xxyy_1[i] * fz_be_0 + g_xxxyyz_xxyy_1[i] * pa_z[i];

        g_xxxyyzz_xxyz_0[i] = 2.0 * g_xyyzz_xxyz_0[i] * fbe_0 - 2.0 * g_xyyzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_xyz_1[i] * fe_0 + g_xxyyzz_xxyz_1[i] * pa_x[i];

        g_xxxyyzz_xxzz_0[i] = g_xxxzz_xxzz_0[i] * fbe_0 - g_xxxzz_xxzz_1[i] * fz_be_0 + g_xxxyzz_xxzz_1[i] * pa_y[i];

        g_xxxyyzz_xyyy_0[i] = g_xxxyy_xyyy_0[i] * fbe_0 - g_xxxyy_xyyy_1[i] * fz_be_0 + g_xxxyyz_xyyy_1[i] * pa_z[i];

        g_xxxyyzz_xyyz_0[i] = 2.0 * g_xyyzz_xyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_xyyz_1[i] * fz_be_0 + g_xxyyzz_yyz_1[i] * fe_0 + g_xxyyzz_xyyz_1[i] * pa_x[i];

        g_xxxyyzz_xyzz_0[i] = 2.0 * g_xyyzz_xyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_xyzz_1[i] * fz_be_0 + g_xxyyzz_yzz_1[i] * fe_0 + g_xxyyzz_xyzz_1[i] * pa_x[i];

        g_xxxyyzz_xzzz_0[i] = g_xxxzz_xzzz_0[i] * fbe_0 - g_xxxzz_xzzz_1[i] * fz_be_0 + g_xxxyzz_xzzz_1[i] * pa_y[i];

        g_xxxyyzz_yyyy_0[i] = 2.0 * g_xyyzz_yyyy_0[i] * fbe_0 - 2.0 * g_xyyzz_yyyy_1[i] * fz_be_0 + g_xxyyzz_yyyy_1[i] * pa_x[i];

        g_xxxyyzz_yyyz_0[i] = 2.0 * g_xyyzz_yyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_yyyz_1[i] * fz_be_0 + g_xxyyzz_yyyz_1[i] * pa_x[i];

        g_xxxyyzz_yyzz_0[i] = 2.0 * g_xyyzz_yyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_yyzz_1[i] * fz_be_0 + g_xxyyzz_yyzz_1[i] * pa_x[i];

        g_xxxyyzz_yzzz_0[i] = 2.0 * g_xyyzz_yzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_yzzz_1[i] * fz_be_0 + g_xxyyzz_yzzz_1[i] * pa_x[i];

        g_xxxyyzz_zzzz_0[i] = 2.0 * g_xyyzz_zzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_zzzz_1[i] * fz_be_0 + g_xxyyzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 195-210 components of targeted buffer : KG

    auto g_xxxyzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 195);

    auto g_xxxyzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 196);

    auto g_xxxyzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 197);

    auto g_xxxyzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 198);

    auto g_xxxyzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 199);

    auto g_xxxyzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 200);

    auto g_xxxyzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 201);

    auto g_xxxyzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 202);

    auto g_xxxyzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 203);

    auto g_xxxyzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 204);

    auto g_xxxyzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 205);

    auto g_xxxyzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 206);

    auto g_xxxyzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 207);

    auto g_xxxyzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 208);

    auto g_xxxyzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 209);

    #pragma omp simd aligned(g_xxxyzzz_xxxx_0, g_xxxyzzz_xxxy_0, g_xxxyzzz_xxxz_0, g_xxxyzzz_xxyy_0, g_xxxyzzz_xxyz_0, g_xxxyzzz_xxzz_0, g_xxxyzzz_xyyy_0, g_xxxyzzz_xyyz_0, g_xxxyzzz_xyzz_0, g_xxxyzzz_xzzz_0, g_xxxyzzz_yyyy_0, g_xxxyzzz_yyyz_0, g_xxxyzzz_yyzz_0, g_xxxyzzz_yzzz_0, g_xxxyzzz_zzzz_0, g_xxxzzz_xxx_1, g_xxxzzz_xxxx_1, g_xxxzzz_xxxy_1, g_xxxzzz_xxxz_1, g_xxxzzz_xxy_1, g_xxxzzz_xxyy_1, g_xxxzzz_xxyz_1, g_xxxzzz_xxz_1, g_xxxzzz_xxzz_1, g_xxxzzz_xyy_1, g_xxxzzz_xyyy_1, g_xxxzzz_xyyz_1, g_xxxzzz_xyz_1, g_xxxzzz_xyzz_1, g_xxxzzz_xzz_1, g_xxxzzz_xzzz_1, g_xxxzzz_yyy_1, g_xxxzzz_yyyy_1, g_xxxzzz_yyyz_1, g_xxxzzz_yyz_1, g_xxxzzz_yyzz_1, g_xxxzzz_yzz_1, g_xxxzzz_yzzz_1, g_xxxzzz_zzz_1, g_xxxzzz_zzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyzzz_xxxx_0[i] = g_xxxzzz_xxxx_1[i] * pa_y[i];

        g_xxxyzzz_xxxy_0[i] = g_xxxzzz_xxx_1[i] * fe_0 + g_xxxzzz_xxxy_1[i] * pa_y[i];

        g_xxxyzzz_xxxz_0[i] = g_xxxzzz_xxxz_1[i] * pa_y[i];

        g_xxxyzzz_xxyy_0[i] = 2.0 * g_xxxzzz_xxy_1[i] * fe_0 + g_xxxzzz_xxyy_1[i] * pa_y[i];

        g_xxxyzzz_xxyz_0[i] = g_xxxzzz_xxz_1[i] * fe_0 + g_xxxzzz_xxyz_1[i] * pa_y[i];

        g_xxxyzzz_xxzz_0[i] = g_xxxzzz_xxzz_1[i] * pa_y[i];

        g_xxxyzzz_xyyy_0[i] = 3.0 * g_xxxzzz_xyy_1[i] * fe_0 + g_xxxzzz_xyyy_1[i] * pa_y[i];

        g_xxxyzzz_xyyz_0[i] = 2.0 * g_xxxzzz_xyz_1[i] * fe_0 + g_xxxzzz_xyyz_1[i] * pa_y[i];

        g_xxxyzzz_xyzz_0[i] = g_xxxzzz_xzz_1[i] * fe_0 + g_xxxzzz_xyzz_1[i] * pa_y[i];

        g_xxxyzzz_xzzz_0[i] = g_xxxzzz_xzzz_1[i] * pa_y[i];

        g_xxxyzzz_yyyy_0[i] = 4.0 * g_xxxzzz_yyy_1[i] * fe_0 + g_xxxzzz_yyyy_1[i] * pa_y[i];

        g_xxxyzzz_yyyz_0[i] = 3.0 * g_xxxzzz_yyz_1[i] * fe_0 + g_xxxzzz_yyyz_1[i] * pa_y[i];

        g_xxxyzzz_yyzz_0[i] = 2.0 * g_xxxzzz_yzz_1[i] * fe_0 + g_xxxzzz_yyzz_1[i] * pa_y[i];

        g_xxxyzzz_yzzz_0[i] = g_xxxzzz_zzz_1[i] * fe_0 + g_xxxzzz_yzzz_1[i] * pa_y[i];

        g_xxxyzzz_zzzz_0[i] = g_xxxzzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 210-225 components of targeted buffer : KG

    auto g_xxxzzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 210);

    auto g_xxxzzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 211);

    auto g_xxxzzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 212);

    auto g_xxxzzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 213);

    auto g_xxxzzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 214);

    auto g_xxxzzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 215);

    auto g_xxxzzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 216);

    auto g_xxxzzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 217);

    auto g_xxxzzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 218);

    auto g_xxxzzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 219);

    auto g_xxxzzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 220);

    auto g_xxxzzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 221);

    auto g_xxxzzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 222);

    auto g_xxxzzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 223);

    auto g_xxxzzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 224);

    #pragma omp simd aligned(g_xxxzz_xxxx_0, g_xxxzz_xxxx_1, g_xxxzz_xxxy_0, g_xxxzz_xxxy_1, g_xxxzz_xxyy_0, g_xxxzz_xxyy_1, g_xxxzz_xyyy_0, g_xxxzz_xyyy_1, g_xxxzzz_xxxx_1, g_xxxzzz_xxxy_1, g_xxxzzz_xxyy_1, g_xxxzzz_xyyy_1, g_xxxzzzz_xxxx_0, g_xxxzzzz_xxxy_0, g_xxxzzzz_xxxz_0, g_xxxzzzz_xxyy_0, g_xxxzzzz_xxyz_0, g_xxxzzzz_xxzz_0, g_xxxzzzz_xyyy_0, g_xxxzzzz_xyyz_0, g_xxxzzzz_xyzz_0, g_xxxzzzz_xzzz_0, g_xxxzzzz_yyyy_0, g_xxxzzzz_yyyz_0, g_xxxzzzz_yyzz_0, g_xxxzzzz_yzzz_0, g_xxxzzzz_zzzz_0, g_xxzzzz_xxxz_1, g_xxzzzz_xxyz_1, g_xxzzzz_xxz_1, g_xxzzzz_xxzz_1, g_xxzzzz_xyyz_1, g_xxzzzz_xyz_1, g_xxzzzz_xyzz_1, g_xxzzzz_xzz_1, g_xxzzzz_xzzz_1, g_xxzzzz_yyyy_1, g_xxzzzz_yyyz_1, g_xxzzzz_yyz_1, g_xxzzzz_yyzz_1, g_xxzzzz_yzz_1, g_xxzzzz_yzzz_1, g_xxzzzz_zzz_1, g_xxzzzz_zzzz_1, g_xzzzz_xxxz_0, g_xzzzz_xxxz_1, g_xzzzz_xxyz_0, g_xzzzz_xxyz_1, g_xzzzz_xxzz_0, g_xzzzz_xxzz_1, g_xzzzz_xyyz_0, g_xzzzz_xyyz_1, g_xzzzz_xyzz_0, g_xzzzz_xyzz_1, g_xzzzz_xzzz_0, g_xzzzz_xzzz_1, g_xzzzz_yyyy_0, g_xzzzz_yyyy_1, g_xzzzz_yyyz_0, g_xzzzz_yyyz_1, g_xzzzz_yyzz_0, g_xzzzz_yyzz_1, g_xzzzz_yzzz_0, g_xzzzz_yzzz_1, g_xzzzz_zzzz_0, g_xzzzz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzzzz_xxxx_0[i] = 3.0 * g_xxxzz_xxxx_0[i] * fbe_0 - 3.0 * g_xxxzz_xxxx_1[i] * fz_be_0 + g_xxxzzz_xxxx_1[i] * pa_z[i];

        g_xxxzzzz_xxxy_0[i] = 3.0 * g_xxxzz_xxxy_0[i] * fbe_0 - 3.0 * g_xxxzz_xxxy_1[i] * fz_be_0 + g_xxxzzz_xxxy_1[i] * pa_z[i];

        g_xxxzzzz_xxxz_0[i] = 2.0 * g_xzzzz_xxxz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxxz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_xxz_1[i] * fe_0 + g_xxzzzz_xxxz_1[i] * pa_x[i];

        g_xxxzzzz_xxyy_0[i] = 3.0 * g_xxxzz_xxyy_0[i] * fbe_0 - 3.0 * g_xxxzz_xxyy_1[i] * fz_be_0 + g_xxxzzz_xxyy_1[i] * pa_z[i];

        g_xxxzzzz_xxyz_0[i] = 2.0 * g_xzzzz_xxyz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_xyz_1[i] * fe_0 + g_xxzzzz_xxyz_1[i] * pa_x[i];

        g_xxxzzzz_xxzz_0[i] = 2.0 * g_xzzzz_xxzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_xzz_1[i] * fe_0 + g_xxzzzz_xxzz_1[i] * pa_x[i];

        g_xxxzzzz_xyyy_0[i] = 3.0 * g_xxxzz_xyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_xyyy_1[i] * fz_be_0 + g_xxxzzz_xyyy_1[i] * pa_z[i];

        g_xxxzzzz_xyyz_0[i] = 2.0 * g_xzzzz_xyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_xyyz_1[i] * fz_be_0 + g_xxzzzz_yyz_1[i] * fe_0 + g_xxzzzz_xyyz_1[i] * pa_x[i];

        g_xxxzzzz_xyzz_0[i] = 2.0 * g_xzzzz_xyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xyzz_1[i] * fz_be_0 + g_xxzzzz_yzz_1[i] * fe_0 + g_xxzzzz_xyzz_1[i] * pa_x[i];

        g_xxxzzzz_xzzz_0[i] = 2.0 * g_xzzzz_xzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xzzz_1[i] * fz_be_0 + g_xxzzzz_zzz_1[i] * fe_0 + g_xxzzzz_xzzz_1[i] * pa_x[i];

        g_xxxzzzz_yyyy_0[i] = 2.0 * g_xzzzz_yyyy_0[i] * fbe_0 - 2.0 * g_xzzzz_yyyy_1[i] * fz_be_0 + g_xxzzzz_yyyy_1[i] * pa_x[i];

        g_xxxzzzz_yyyz_0[i] = 2.0 * g_xzzzz_yyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_yyyz_1[i] * fz_be_0 + g_xxzzzz_yyyz_1[i] * pa_x[i];

        g_xxxzzzz_yyzz_0[i] = 2.0 * g_xzzzz_yyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_yyzz_1[i] * fz_be_0 + g_xxzzzz_yyzz_1[i] * pa_x[i];

        g_xxxzzzz_yzzz_0[i] = 2.0 * g_xzzzz_yzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_yzzz_1[i] * fz_be_0 + g_xxzzzz_yzzz_1[i] * pa_x[i];

        g_xxxzzzz_zzzz_0[i] = 2.0 * g_xzzzz_zzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_zzzz_1[i] * fz_be_0 + g_xxzzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 225-240 components of targeted buffer : KG

    auto g_xxyyyyy_xxxx_0 = pbuffer.data(idx_eri_0_kg + 225);

    auto g_xxyyyyy_xxxy_0 = pbuffer.data(idx_eri_0_kg + 226);

    auto g_xxyyyyy_xxxz_0 = pbuffer.data(idx_eri_0_kg + 227);

    auto g_xxyyyyy_xxyy_0 = pbuffer.data(idx_eri_0_kg + 228);

    auto g_xxyyyyy_xxyz_0 = pbuffer.data(idx_eri_0_kg + 229);

    auto g_xxyyyyy_xxzz_0 = pbuffer.data(idx_eri_0_kg + 230);

    auto g_xxyyyyy_xyyy_0 = pbuffer.data(idx_eri_0_kg + 231);

    auto g_xxyyyyy_xyyz_0 = pbuffer.data(idx_eri_0_kg + 232);

    auto g_xxyyyyy_xyzz_0 = pbuffer.data(idx_eri_0_kg + 233);

    auto g_xxyyyyy_xzzz_0 = pbuffer.data(idx_eri_0_kg + 234);

    auto g_xxyyyyy_yyyy_0 = pbuffer.data(idx_eri_0_kg + 235);

    auto g_xxyyyyy_yyyz_0 = pbuffer.data(idx_eri_0_kg + 236);

    auto g_xxyyyyy_yyzz_0 = pbuffer.data(idx_eri_0_kg + 237);

    auto g_xxyyyyy_yzzz_0 = pbuffer.data(idx_eri_0_kg + 238);

    auto g_xxyyyyy_zzzz_0 = pbuffer.data(idx_eri_0_kg + 239);

    #pragma omp simd aligned(g_xxyyy_xxxx_0, g_xxyyy_xxxx_1, g_xxyyy_xxxz_0, g_xxyyy_xxxz_1, g_xxyyy_xxzz_0, g_xxyyy_xxzz_1, g_xxyyy_xzzz_0, g_xxyyy_xzzz_1, g_xxyyyy_xxxx_1, g_xxyyyy_xxxz_1, g_xxyyyy_xxzz_1, g_xxyyyy_xzzz_1, g_xxyyyyy_xxxx_0, g_xxyyyyy_xxxy_0, g_xxyyyyy_xxxz_0, g_xxyyyyy_xxyy_0, g_xxyyyyy_xxyz_0, g_xxyyyyy_xxzz_0, g_xxyyyyy_xyyy_0, g_xxyyyyy_xyyz_0, g_xxyyyyy_xyzz_0, g_xxyyyyy_xzzz_0, g_xxyyyyy_yyyy_0, g_xxyyyyy_yyyz_0, g_xxyyyyy_yyzz_0, g_xxyyyyy_yzzz_0, g_xxyyyyy_zzzz_0, g_xyyyyy_xxxy_1, g_xyyyyy_xxy_1, g_xyyyyy_xxyy_1, g_xyyyyy_xxyz_1, g_xyyyyy_xyy_1, g_xyyyyy_xyyy_1, g_xyyyyy_xyyz_1, g_xyyyyy_xyz_1, g_xyyyyy_xyzz_1, g_xyyyyy_yyy_1, g_xyyyyy_yyyy_1, g_xyyyyy_yyyz_1, g_xyyyyy_yyz_1, g_xyyyyy_yyzz_1, g_xyyyyy_yzz_1, g_xyyyyy_yzzz_1, g_xyyyyy_zzzz_1, g_yyyyy_xxxy_0, g_yyyyy_xxxy_1, g_yyyyy_xxyy_0, g_yyyyy_xxyy_1, g_yyyyy_xxyz_0, g_yyyyy_xxyz_1, g_yyyyy_xyyy_0, g_yyyyy_xyyy_1, g_yyyyy_xyyz_0, g_yyyyy_xyyz_1, g_yyyyy_xyzz_0, g_yyyyy_xyzz_1, g_yyyyy_yyyy_0, g_yyyyy_yyyy_1, g_yyyyy_yyyz_0, g_yyyyy_yyyz_1, g_yyyyy_yyzz_0, g_yyyyy_yyzz_1, g_yyyyy_yzzz_0, g_yyyyy_yzzz_1, g_yyyyy_zzzz_0, g_yyyyy_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyyy_xxxx_0[i] = 4.0 * g_xxyyy_xxxx_0[i] * fbe_0 - 4.0 * g_xxyyy_xxxx_1[i] * fz_be_0 + g_xxyyyy_xxxx_1[i] * pa_y[i];

        g_xxyyyyy_xxxy_0[i] = g_yyyyy_xxxy_0[i] * fbe_0 - g_yyyyy_xxxy_1[i] * fz_be_0 + 3.0 * g_xyyyyy_xxy_1[i] * fe_0 + g_xyyyyy_xxxy_1[i] * pa_x[i];

        g_xxyyyyy_xxxz_0[i] = 4.0 * g_xxyyy_xxxz_0[i] * fbe_0 - 4.0 * g_xxyyy_xxxz_1[i] * fz_be_0 + g_xxyyyy_xxxz_1[i] * pa_y[i];

        g_xxyyyyy_xxyy_0[i] = g_yyyyy_xxyy_0[i] * fbe_0 - g_yyyyy_xxyy_1[i] * fz_be_0 + 2.0 * g_xyyyyy_xyy_1[i] * fe_0 + g_xyyyyy_xxyy_1[i] * pa_x[i];

        g_xxyyyyy_xxyz_0[i] = g_yyyyy_xxyz_0[i] * fbe_0 - g_yyyyy_xxyz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_xyz_1[i] * fe_0 + g_xyyyyy_xxyz_1[i] * pa_x[i];

        g_xxyyyyy_xxzz_0[i] = 4.0 * g_xxyyy_xxzz_0[i] * fbe_0 - 4.0 * g_xxyyy_xxzz_1[i] * fz_be_0 + g_xxyyyy_xxzz_1[i] * pa_y[i];

        g_xxyyyyy_xyyy_0[i] = g_yyyyy_xyyy_0[i] * fbe_0 - g_yyyyy_xyyy_1[i] * fz_be_0 + g_xyyyyy_yyy_1[i] * fe_0 + g_xyyyyy_xyyy_1[i] * pa_x[i];

        g_xxyyyyy_xyyz_0[i] = g_yyyyy_xyyz_0[i] * fbe_0 - g_yyyyy_xyyz_1[i] * fz_be_0 + g_xyyyyy_yyz_1[i] * fe_0 + g_xyyyyy_xyyz_1[i] * pa_x[i];

        g_xxyyyyy_xyzz_0[i] = g_yyyyy_xyzz_0[i] * fbe_0 - g_yyyyy_xyzz_1[i] * fz_be_0 + g_xyyyyy_yzz_1[i] * fe_0 + g_xyyyyy_xyzz_1[i] * pa_x[i];

        g_xxyyyyy_xzzz_0[i] = 4.0 * g_xxyyy_xzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_xzzz_1[i] * fz_be_0 + g_xxyyyy_xzzz_1[i] * pa_y[i];

        g_xxyyyyy_yyyy_0[i] = g_yyyyy_yyyy_0[i] * fbe_0 - g_yyyyy_yyyy_1[i] * fz_be_0 + g_xyyyyy_yyyy_1[i] * pa_x[i];

        g_xxyyyyy_yyyz_0[i] = g_yyyyy_yyyz_0[i] * fbe_0 - g_yyyyy_yyyz_1[i] * fz_be_0 + g_xyyyyy_yyyz_1[i] * pa_x[i];

        g_xxyyyyy_yyzz_0[i] = g_yyyyy_yyzz_0[i] * fbe_0 - g_yyyyy_yyzz_1[i] * fz_be_0 + g_xyyyyy_yyzz_1[i] * pa_x[i];

        g_xxyyyyy_yzzz_0[i] = g_yyyyy_yzzz_0[i] * fbe_0 - g_yyyyy_yzzz_1[i] * fz_be_0 + g_xyyyyy_yzzz_1[i] * pa_x[i];

        g_xxyyyyy_zzzz_0[i] = g_yyyyy_zzzz_0[i] * fbe_0 - g_yyyyy_zzzz_1[i] * fz_be_0 + g_xyyyyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 240-255 components of targeted buffer : KG

    auto g_xxyyyyz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 240);

    auto g_xxyyyyz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 241);

    auto g_xxyyyyz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 242);

    auto g_xxyyyyz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 243);

    auto g_xxyyyyz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 244);

    auto g_xxyyyyz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 245);

    auto g_xxyyyyz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 246);

    auto g_xxyyyyz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 247);

    auto g_xxyyyyz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 248);

    auto g_xxyyyyz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 249);

    auto g_xxyyyyz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 250);

    auto g_xxyyyyz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 251);

    auto g_xxyyyyz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 252);

    auto g_xxyyyyz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 253);

    auto g_xxyyyyz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 254);

    #pragma omp simd aligned(g_xxyyyy_xxx_1, g_xxyyyy_xxxx_1, g_xxyyyy_xxxy_1, g_xxyyyy_xxxz_1, g_xxyyyy_xxy_1, g_xxyyyy_xxyy_1, g_xxyyyy_xxyz_1, g_xxyyyy_xxz_1, g_xxyyyy_xxzz_1, g_xxyyyy_xyy_1, g_xxyyyy_xyyy_1, g_xxyyyy_xyyz_1, g_xxyyyy_xyz_1, g_xxyyyy_xyzz_1, g_xxyyyy_xzz_1, g_xxyyyy_xzzz_1, g_xxyyyy_yyy_1, g_xxyyyy_yyyy_1, g_xxyyyy_yyyz_1, g_xxyyyy_yyz_1, g_xxyyyy_yyzz_1, g_xxyyyy_yzz_1, g_xxyyyy_yzzz_1, g_xxyyyy_zzz_1, g_xxyyyy_zzzz_1, g_xxyyyyz_xxxx_0, g_xxyyyyz_xxxy_0, g_xxyyyyz_xxxz_0, g_xxyyyyz_xxyy_0, g_xxyyyyz_xxyz_0, g_xxyyyyz_xxzz_0, g_xxyyyyz_xyyy_0, g_xxyyyyz_xyyz_0, g_xxyyyyz_xyzz_0, g_xxyyyyz_xzzz_0, g_xxyyyyz_yyyy_0, g_xxyyyyz_yyyz_0, g_xxyyyyz_yyzz_0, g_xxyyyyz_yzzz_0, g_xxyyyyz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyyyz_xxxx_0[i] = g_xxyyyy_xxxx_1[i] * pa_z[i];

        g_xxyyyyz_xxxy_0[i] = g_xxyyyy_xxxy_1[i] * pa_z[i];

        g_xxyyyyz_xxxz_0[i] = g_xxyyyy_xxx_1[i] * fe_0 + g_xxyyyy_xxxz_1[i] * pa_z[i];

        g_xxyyyyz_xxyy_0[i] = g_xxyyyy_xxyy_1[i] * pa_z[i];

        g_xxyyyyz_xxyz_0[i] = g_xxyyyy_xxy_1[i] * fe_0 + g_xxyyyy_xxyz_1[i] * pa_z[i];

        g_xxyyyyz_xxzz_0[i] = 2.0 * g_xxyyyy_xxz_1[i] * fe_0 + g_xxyyyy_xxzz_1[i] * pa_z[i];

        g_xxyyyyz_xyyy_0[i] = g_xxyyyy_xyyy_1[i] * pa_z[i];

        g_xxyyyyz_xyyz_0[i] = g_xxyyyy_xyy_1[i] * fe_0 + g_xxyyyy_xyyz_1[i] * pa_z[i];

        g_xxyyyyz_xyzz_0[i] = 2.0 * g_xxyyyy_xyz_1[i] * fe_0 + g_xxyyyy_xyzz_1[i] * pa_z[i];

        g_xxyyyyz_xzzz_0[i] = 3.0 * g_xxyyyy_xzz_1[i] * fe_0 + g_xxyyyy_xzzz_1[i] * pa_z[i];

        g_xxyyyyz_yyyy_0[i] = g_xxyyyy_yyyy_1[i] * pa_z[i];

        g_xxyyyyz_yyyz_0[i] = g_xxyyyy_yyy_1[i] * fe_0 + g_xxyyyy_yyyz_1[i] * pa_z[i];

        g_xxyyyyz_yyzz_0[i] = 2.0 * g_xxyyyy_yyz_1[i] * fe_0 + g_xxyyyy_yyzz_1[i] * pa_z[i];

        g_xxyyyyz_yzzz_0[i] = 3.0 * g_xxyyyy_yzz_1[i] * fe_0 + g_xxyyyy_yzzz_1[i] * pa_z[i];

        g_xxyyyyz_zzzz_0[i] = 4.0 * g_xxyyyy_zzz_1[i] * fe_0 + g_xxyyyy_zzzz_1[i] * pa_z[i];
    }

    // Set up 255-270 components of targeted buffer : KG

    auto g_xxyyyzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 255);

    auto g_xxyyyzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 256);

    auto g_xxyyyzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 257);

    auto g_xxyyyzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 258);

    auto g_xxyyyzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 259);

    auto g_xxyyyzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 260);

    auto g_xxyyyzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 261);

    auto g_xxyyyzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 262);

    auto g_xxyyyzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 263);

    auto g_xxyyyzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 264);

    auto g_xxyyyzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 265);

    auto g_xxyyyzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 266);

    auto g_xxyyyzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 267);

    auto g_xxyyyzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 268);

    auto g_xxyyyzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 269);

    #pragma omp simd aligned(g_xxyyy_xxxy_0, g_xxyyy_xxxy_1, g_xxyyy_xxyy_0, g_xxyyy_xxyy_1, g_xxyyy_xyyy_0, g_xxyyy_xyyy_1, g_xxyyyz_xxxy_1, g_xxyyyz_xxyy_1, g_xxyyyz_xyyy_1, g_xxyyyzz_xxxx_0, g_xxyyyzz_xxxy_0, g_xxyyyzz_xxxz_0, g_xxyyyzz_xxyy_0, g_xxyyyzz_xxyz_0, g_xxyyyzz_xxzz_0, g_xxyyyzz_xyyy_0, g_xxyyyzz_xyyz_0, g_xxyyyzz_xyzz_0, g_xxyyyzz_xzzz_0, g_xxyyyzz_yyyy_0, g_xxyyyzz_yyyz_0, g_xxyyyzz_yyzz_0, g_xxyyyzz_yzzz_0, g_xxyyyzz_zzzz_0, g_xxyyzz_xxxx_1, g_xxyyzz_xxxz_1, g_xxyyzz_xxzz_1, g_xxyyzz_xzzz_1, g_xxyzz_xxxx_0, g_xxyzz_xxxx_1, g_xxyzz_xxxz_0, g_xxyzz_xxxz_1, g_xxyzz_xxzz_0, g_xxyzz_xxzz_1, g_xxyzz_xzzz_0, g_xxyzz_xzzz_1, g_xyyyzz_xxyz_1, g_xyyyzz_xyyz_1, g_xyyyzz_xyz_1, g_xyyyzz_xyzz_1, g_xyyyzz_yyyy_1, g_xyyyzz_yyyz_1, g_xyyyzz_yyz_1, g_xyyyzz_yyzz_1, g_xyyyzz_yzz_1, g_xyyyzz_yzzz_1, g_xyyyzz_zzzz_1, g_yyyzz_xxyz_0, g_yyyzz_xxyz_1, g_yyyzz_xyyz_0, g_yyyzz_xyyz_1, g_yyyzz_xyzz_0, g_yyyzz_xyzz_1, g_yyyzz_yyyy_0, g_yyyzz_yyyy_1, g_yyyzz_yyyz_0, g_yyyzz_yyyz_1, g_yyyzz_yyzz_0, g_yyyzz_yyzz_1, g_yyyzz_yzzz_0, g_yyyzz_yzzz_1, g_yyyzz_zzzz_0, g_yyyzz_zzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyzz_xxxx_0[i] = 2.0 * g_xxyzz_xxxx_0[i] * fbe_0 - 2.0 * g_xxyzz_xxxx_1[i] * fz_be_0 + g_xxyyzz_xxxx_1[i] * pa_y[i];

        g_xxyyyzz_xxxy_0[i] = g_xxyyy_xxxy_0[i] * fbe_0 - g_xxyyy_xxxy_1[i] * fz_be_0 + g_xxyyyz_xxxy_1[i] * pa_z[i];

        g_xxyyyzz_xxxz_0[i] = 2.0 * g_xxyzz_xxxz_0[i] * fbe_0 - 2.0 * g_xxyzz_xxxz_1[i] * fz_be_0 + g_xxyyzz_xxxz_1[i] * pa_y[i];

        g_xxyyyzz_xxyy_0[i] = g_xxyyy_xxyy_0[i] * fbe_0 - g_xxyyy_xxyy_1[i] * fz_be_0 + g_xxyyyz_xxyy_1[i] * pa_z[i];

        g_xxyyyzz_xxyz_0[i] = g_yyyzz_xxyz_0[i] * fbe_0 - g_yyyzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_xyz_1[i] * fe_0 + g_xyyyzz_xxyz_1[i] * pa_x[i];

        g_xxyyyzz_xxzz_0[i] = 2.0 * g_xxyzz_xxzz_0[i] * fbe_0 - 2.0 * g_xxyzz_xxzz_1[i] * fz_be_0 + g_xxyyzz_xxzz_1[i] * pa_y[i];

        g_xxyyyzz_xyyy_0[i] = g_xxyyy_xyyy_0[i] * fbe_0 - g_xxyyy_xyyy_1[i] * fz_be_0 + g_xxyyyz_xyyy_1[i] * pa_z[i];

        g_xxyyyzz_xyyz_0[i] = g_yyyzz_xyyz_0[i] * fbe_0 - g_yyyzz_xyyz_1[i] * fz_be_0 + g_xyyyzz_yyz_1[i] * fe_0 + g_xyyyzz_xyyz_1[i] * pa_x[i];

        g_xxyyyzz_xyzz_0[i] = g_yyyzz_xyzz_0[i] * fbe_0 - g_yyyzz_xyzz_1[i] * fz_be_0 + g_xyyyzz_yzz_1[i] * fe_0 + g_xyyyzz_xyzz_1[i] * pa_x[i];

        g_xxyyyzz_xzzz_0[i] = 2.0 * g_xxyzz_xzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_xzzz_1[i] * fz_be_0 + g_xxyyzz_xzzz_1[i] * pa_y[i];

        g_xxyyyzz_yyyy_0[i] = g_yyyzz_yyyy_0[i] * fbe_0 - g_yyyzz_yyyy_1[i] * fz_be_0 + g_xyyyzz_yyyy_1[i] * pa_x[i];

        g_xxyyyzz_yyyz_0[i] = g_yyyzz_yyyz_0[i] * fbe_0 - g_yyyzz_yyyz_1[i] * fz_be_0 + g_xyyyzz_yyyz_1[i] * pa_x[i];

        g_xxyyyzz_yyzz_0[i] = g_yyyzz_yyzz_0[i] * fbe_0 - g_yyyzz_yyzz_1[i] * fz_be_0 + g_xyyyzz_yyzz_1[i] * pa_x[i];

        g_xxyyyzz_yzzz_0[i] = g_yyyzz_yzzz_0[i] * fbe_0 - g_yyyzz_yzzz_1[i] * fz_be_0 + g_xyyyzz_yzzz_1[i] * pa_x[i];

        g_xxyyyzz_zzzz_0[i] = g_yyyzz_zzzz_0[i] * fbe_0 - g_yyyzz_zzzz_1[i] * fz_be_0 + g_xyyyzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 270-285 components of targeted buffer : KG

    auto g_xxyyzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 270);

    auto g_xxyyzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 271);

    auto g_xxyyzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 272);

    auto g_xxyyzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 273);

    auto g_xxyyzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 274);

    auto g_xxyyzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 275);

    auto g_xxyyzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 276);

    auto g_xxyyzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 277);

    auto g_xxyyzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 278);

    auto g_xxyyzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 279);

    auto g_xxyyzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 280);

    auto g_xxyyzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 281);

    auto g_xxyyzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 282);

    auto g_xxyyzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 283);

    auto g_xxyyzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 284);

    #pragma omp simd aligned(g_xxyyz_xxxy_0, g_xxyyz_xxxy_1, g_xxyyz_xxyy_0, g_xxyyz_xxyy_1, g_xxyyz_xyyy_0, g_xxyyz_xyyy_1, g_xxyyzz_xxxy_1, g_xxyyzz_xxyy_1, g_xxyyzz_xyyy_1, g_xxyyzzz_xxxx_0, g_xxyyzzz_xxxy_0, g_xxyyzzz_xxxz_0, g_xxyyzzz_xxyy_0, g_xxyyzzz_xxyz_0, g_xxyyzzz_xxzz_0, g_xxyyzzz_xyyy_0, g_xxyyzzz_xyyz_0, g_xxyyzzz_xyzz_0, g_xxyyzzz_xzzz_0, g_xxyyzzz_yyyy_0, g_xxyyzzz_yyyz_0, g_xxyyzzz_yyzz_0, g_xxyyzzz_yzzz_0, g_xxyyzzz_zzzz_0, g_xxyzzz_xxxx_1, g_xxyzzz_xxxz_1, g_xxyzzz_xxzz_1, g_xxyzzz_xzzz_1, g_xxzzz_xxxx_0, g_xxzzz_xxxx_1, g_xxzzz_xxxz_0, g_xxzzz_xxxz_1, g_xxzzz_xxzz_0, g_xxzzz_xxzz_1, g_xxzzz_xzzz_0, g_xxzzz_xzzz_1, g_xyyzzz_xxyz_1, g_xyyzzz_xyyz_1, g_xyyzzz_xyz_1, g_xyyzzz_xyzz_1, g_xyyzzz_yyyy_1, g_xyyzzz_yyyz_1, g_xyyzzz_yyz_1, g_xyyzzz_yyzz_1, g_xyyzzz_yzz_1, g_xyyzzz_yzzz_1, g_xyyzzz_zzzz_1, g_yyzzz_xxyz_0, g_yyzzz_xxyz_1, g_yyzzz_xyyz_0, g_yyzzz_xyyz_1, g_yyzzz_xyzz_0, g_yyzzz_xyzz_1, g_yyzzz_yyyy_0, g_yyzzz_yyyy_1, g_yyzzz_yyyz_0, g_yyzzz_yyyz_1, g_yyzzz_yyzz_0, g_yyzzz_yyzz_1, g_yyzzz_yzzz_0, g_yyzzz_yzzz_1, g_yyzzz_zzzz_0, g_yyzzz_zzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyzzz_xxxx_0[i] = g_xxzzz_xxxx_0[i] * fbe_0 - g_xxzzz_xxxx_1[i] * fz_be_0 + g_xxyzzz_xxxx_1[i] * pa_y[i];

        g_xxyyzzz_xxxy_0[i] = 2.0 * g_xxyyz_xxxy_0[i] * fbe_0 - 2.0 * g_xxyyz_xxxy_1[i] * fz_be_0 + g_xxyyzz_xxxy_1[i] * pa_z[i];

        g_xxyyzzz_xxxz_0[i] = g_xxzzz_xxxz_0[i] * fbe_0 - g_xxzzz_xxxz_1[i] * fz_be_0 + g_xxyzzz_xxxz_1[i] * pa_y[i];

        g_xxyyzzz_xxyy_0[i] = 2.0 * g_xxyyz_xxyy_0[i] * fbe_0 - 2.0 * g_xxyyz_xxyy_1[i] * fz_be_0 + g_xxyyzz_xxyy_1[i] * pa_z[i];

        g_xxyyzzz_xxyz_0[i] = g_yyzzz_xxyz_0[i] * fbe_0 - g_yyzzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_xyz_1[i] * fe_0 + g_xyyzzz_xxyz_1[i] * pa_x[i];

        g_xxyyzzz_xxzz_0[i] = g_xxzzz_xxzz_0[i] * fbe_0 - g_xxzzz_xxzz_1[i] * fz_be_0 + g_xxyzzz_xxzz_1[i] * pa_y[i];

        g_xxyyzzz_xyyy_0[i] = 2.0 * g_xxyyz_xyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_xyyy_1[i] * fz_be_0 + g_xxyyzz_xyyy_1[i] * pa_z[i];

        g_xxyyzzz_xyyz_0[i] = g_yyzzz_xyyz_0[i] * fbe_0 - g_yyzzz_xyyz_1[i] * fz_be_0 + g_xyyzzz_yyz_1[i] * fe_0 + g_xyyzzz_xyyz_1[i] * pa_x[i];

        g_xxyyzzz_xyzz_0[i] = g_yyzzz_xyzz_0[i] * fbe_0 - g_yyzzz_xyzz_1[i] * fz_be_0 + g_xyyzzz_yzz_1[i] * fe_0 + g_xyyzzz_xyzz_1[i] * pa_x[i];

        g_xxyyzzz_xzzz_0[i] = g_xxzzz_xzzz_0[i] * fbe_0 - g_xxzzz_xzzz_1[i] * fz_be_0 + g_xxyzzz_xzzz_1[i] * pa_y[i];

        g_xxyyzzz_yyyy_0[i] = g_yyzzz_yyyy_0[i] * fbe_0 - g_yyzzz_yyyy_1[i] * fz_be_0 + g_xyyzzz_yyyy_1[i] * pa_x[i];

        g_xxyyzzz_yyyz_0[i] = g_yyzzz_yyyz_0[i] * fbe_0 - g_yyzzz_yyyz_1[i] * fz_be_0 + g_xyyzzz_yyyz_1[i] * pa_x[i];

        g_xxyyzzz_yyzz_0[i] = g_yyzzz_yyzz_0[i] * fbe_0 - g_yyzzz_yyzz_1[i] * fz_be_0 + g_xyyzzz_yyzz_1[i] * pa_x[i];

        g_xxyyzzz_yzzz_0[i] = g_yyzzz_yzzz_0[i] * fbe_0 - g_yyzzz_yzzz_1[i] * fz_be_0 + g_xyyzzz_yzzz_1[i] * pa_x[i];

        g_xxyyzzz_zzzz_0[i] = g_yyzzz_zzzz_0[i] * fbe_0 - g_yyzzz_zzzz_1[i] * fz_be_0 + g_xyyzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 285-300 components of targeted buffer : KG

    auto g_xxyzzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 285);

    auto g_xxyzzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 286);

    auto g_xxyzzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 287);

    auto g_xxyzzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 288);

    auto g_xxyzzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 289);

    auto g_xxyzzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 290);

    auto g_xxyzzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 291);

    auto g_xxyzzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 292);

    auto g_xxyzzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 293);

    auto g_xxyzzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 294);

    auto g_xxyzzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 295);

    auto g_xxyzzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 296);

    auto g_xxyzzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 297);

    auto g_xxyzzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 298);

    auto g_xxyzzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 299);

    #pragma omp simd aligned(g_xxyzzzz_xxxx_0, g_xxyzzzz_xxxy_0, g_xxyzzzz_xxxz_0, g_xxyzzzz_xxyy_0, g_xxyzzzz_xxyz_0, g_xxyzzzz_xxzz_0, g_xxyzzzz_xyyy_0, g_xxyzzzz_xyyz_0, g_xxyzzzz_xyzz_0, g_xxyzzzz_xzzz_0, g_xxyzzzz_yyyy_0, g_xxyzzzz_yyyz_0, g_xxyzzzz_yyzz_0, g_xxyzzzz_yzzz_0, g_xxyzzzz_zzzz_0, g_xxzzzz_xxx_1, g_xxzzzz_xxxx_1, g_xxzzzz_xxxy_1, g_xxzzzz_xxxz_1, g_xxzzzz_xxy_1, g_xxzzzz_xxyy_1, g_xxzzzz_xxyz_1, g_xxzzzz_xxz_1, g_xxzzzz_xxzz_1, g_xxzzzz_xyy_1, g_xxzzzz_xyyy_1, g_xxzzzz_xyyz_1, g_xxzzzz_xyz_1, g_xxzzzz_xyzz_1, g_xxzzzz_xzz_1, g_xxzzzz_xzzz_1, g_xxzzzz_yyy_1, g_xxzzzz_yyyy_1, g_xxzzzz_yyyz_1, g_xxzzzz_yyz_1, g_xxzzzz_yyzz_1, g_xxzzzz_yzz_1, g_xxzzzz_yzzz_1, g_xxzzzz_zzz_1, g_xxzzzz_zzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzzzz_xxxx_0[i] = g_xxzzzz_xxxx_1[i] * pa_y[i];

        g_xxyzzzz_xxxy_0[i] = g_xxzzzz_xxx_1[i] * fe_0 + g_xxzzzz_xxxy_1[i] * pa_y[i];

        g_xxyzzzz_xxxz_0[i] = g_xxzzzz_xxxz_1[i] * pa_y[i];

        g_xxyzzzz_xxyy_0[i] = 2.0 * g_xxzzzz_xxy_1[i] * fe_0 + g_xxzzzz_xxyy_1[i] * pa_y[i];

        g_xxyzzzz_xxyz_0[i] = g_xxzzzz_xxz_1[i] * fe_0 + g_xxzzzz_xxyz_1[i] * pa_y[i];

        g_xxyzzzz_xxzz_0[i] = g_xxzzzz_xxzz_1[i] * pa_y[i];

        g_xxyzzzz_xyyy_0[i] = 3.0 * g_xxzzzz_xyy_1[i] * fe_0 + g_xxzzzz_xyyy_1[i] * pa_y[i];

        g_xxyzzzz_xyyz_0[i] = 2.0 * g_xxzzzz_xyz_1[i] * fe_0 + g_xxzzzz_xyyz_1[i] * pa_y[i];

        g_xxyzzzz_xyzz_0[i] = g_xxzzzz_xzz_1[i] * fe_0 + g_xxzzzz_xyzz_1[i] * pa_y[i];

        g_xxyzzzz_xzzz_0[i] = g_xxzzzz_xzzz_1[i] * pa_y[i];

        g_xxyzzzz_yyyy_0[i] = 4.0 * g_xxzzzz_yyy_1[i] * fe_0 + g_xxzzzz_yyyy_1[i] * pa_y[i];

        g_xxyzzzz_yyyz_0[i] = 3.0 * g_xxzzzz_yyz_1[i] * fe_0 + g_xxzzzz_yyyz_1[i] * pa_y[i];

        g_xxyzzzz_yyzz_0[i] = 2.0 * g_xxzzzz_yzz_1[i] * fe_0 + g_xxzzzz_yyzz_1[i] * pa_y[i];

        g_xxyzzzz_yzzz_0[i] = g_xxzzzz_zzz_1[i] * fe_0 + g_xxzzzz_yzzz_1[i] * pa_y[i];

        g_xxyzzzz_zzzz_0[i] = g_xxzzzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 300-315 components of targeted buffer : KG

    auto g_xxzzzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 300);

    auto g_xxzzzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 301);

    auto g_xxzzzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 302);

    auto g_xxzzzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 303);

    auto g_xxzzzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 304);

    auto g_xxzzzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 305);

    auto g_xxzzzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 306);

    auto g_xxzzzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 307);

    auto g_xxzzzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 308);

    auto g_xxzzzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 309);

    auto g_xxzzzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 310);

    auto g_xxzzzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 311);

    auto g_xxzzzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 312);

    auto g_xxzzzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 313);

    auto g_xxzzzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 314);

    #pragma omp simd aligned(g_xxzzz_xxxx_0, g_xxzzz_xxxx_1, g_xxzzz_xxxy_0, g_xxzzz_xxxy_1, g_xxzzz_xxyy_0, g_xxzzz_xxyy_1, g_xxzzz_xyyy_0, g_xxzzz_xyyy_1, g_xxzzzz_xxxx_1, g_xxzzzz_xxxy_1, g_xxzzzz_xxyy_1, g_xxzzzz_xyyy_1, g_xxzzzzz_xxxx_0, g_xxzzzzz_xxxy_0, g_xxzzzzz_xxxz_0, g_xxzzzzz_xxyy_0, g_xxzzzzz_xxyz_0, g_xxzzzzz_xxzz_0, g_xxzzzzz_xyyy_0, g_xxzzzzz_xyyz_0, g_xxzzzzz_xyzz_0, g_xxzzzzz_xzzz_0, g_xxzzzzz_yyyy_0, g_xxzzzzz_yyyz_0, g_xxzzzzz_yyzz_0, g_xxzzzzz_yzzz_0, g_xxzzzzz_zzzz_0, g_xzzzzz_xxxz_1, g_xzzzzz_xxyz_1, g_xzzzzz_xxz_1, g_xzzzzz_xxzz_1, g_xzzzzz_xyyz_1, g_xzzzzz_xyz_1, g_xzzzzz_xyzz_1, g_xzzzzz_xzz_1, g_xzzzzz_xzzz_1, g_xzzzzz_yyyy_1, g_xzzzzz_yyyz_1, g_xzzzzz_yyz_1, g_xzzzzz_yyzz_1, g_xzzzzz_yzz_1, g_xzzzzz_yzzz_1, g_xzzzzz_zzz_1, g_xzzzzz_zzzz_1, g_zzzzz_xxxz_0, g_zzzzz_xxxz_1, g_zzzzz_xxyz_0, g_zzzzz_xxyz_1, g_zzzzz_xxzz_0, g_zzzzz_xxzz_1, g_zzzzz_xyyz_0, g_zzzzz_xyyz_1, g_zzzzz_xyzz_0, g_zzzzz_xyzz_1, g_zzzzz_xzzz_0, g_zzzzz_xzzz_1, g_zzzzz_yyyy_0, g_zzzzz_yyyy_1, g_zzzzz_yyyz_0, g_zzzzz_yyyz_1, g_zzzzz_yyzz_0, g_zzzzz_yyzz_1, g_zzzzz_yzzz_0, g_zzzzz_yzzz_1, g_zzzzz_zzzz_0, g_zzzzz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzzzz_xxxx_0[i] = 4.0 * g_xxzzz_xxxx_0[i] * fbe_0 - 4.0 * g_xxzzz_xxxx_1[i] * fz_be_0 + g_xxzzzz_xxxx_1[i] * pa_z[i];

        g_xxzzzzz_xxxy_0[i] = 4.0 * g_xxzzz_xxxy_0[i] * fbe_0 - 4.0 * g_xxzzz_xxxy_1[i] * fz_be_0 + g_xxzzzz_xxxy_1[i] * pa_z[i];

        g_xxzzzzz_xxxz_0[i] = g_zzzzz_xxxz_0[i] * fbe_0 - g_zzzzz_xxxz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_xxz_1[i] * fe_0 + g_xzzzzz_xxxz_1[i] * pa_x[i];

        g_xxzzzzz_xxyy_0[i] = 4.0 * g_xxzzz_xxyy_0[i] * fbe_0 - 4.0 * g_xxzzz_xxyy_1[i] * fz_be_0 + g_xxzzzz_xxyy_1[i] * pa_z[i];

        g_xxzzzzz_xxyz_0[i] = g_zzzzz_xxyz_0[i] * fbe_0 - g_zzzzz_xxyz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_xyz_1[i] * fe_0 + g_xzzzzz_xxyz_1[i] * pa_x[i];

        g_xxzzzzz_xxzz_0[i] = g_zzzzz_xxzz_0[i] * fbe_0 - g_zzzzz_xxzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_xzz_1[i] * fe_0 + g_xzzzzz_xxzz_1[i] * pa_x[i];

        g_xxzzzzz_xyyy_0[i] = 4.0 * g_xxzzz_xyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_xyyy_1[i] * fz_be_0 + g_xxzzzz_xyyy_1[i] * pa_z[i];

        g_xxzzzzz_xyyz_0[i] = g_zzzzz_xyyz_0[i] * fbe_0 - g_zzzzz_xyyz_1[i] * fz_be_0 + g_xzzzzz_yyz_1[i] * fe_0 + g_xzzzzz_xyyz_1[i] * pa_x[i];

        g_xxzzzzz_xyzz_0[i] = g_zzzzz_xyzz_0[i] * fbe_0 - g_zzzzz_xyzz_1[i] * fz_be_0 + g_xzzzzz_yzz_1[i] * fe_0 + g_xzzzzz_xyzz_1[i] * pa_x[i];

        g_xxzzzzz_xzzz_0[i] = g_zzzzz_xzzz_0[i] * fbe_0 - g_zzzzz_xzzz_1[i] * fz_be_0 + g_xzzzzz_zzz_1[i] * fe_0 + g_xzzzzz_xzzz_1[i] * pa_x[i];

        g_xxzzzzz_yyyy_0[i] = g_zzzzz_yyyy_0[i] * fbe_0 - g_zzzzz_yyyy_1[i] * fz_be_0 + g_xzzzzz_yyyy_1[i] * pa_x[i];

        g_xxzzzzz_yyyz_0[i] = g_zzzzz_yyyz_0[i] * fbe_0 - g_zzzzz_yyyz_1[i] * fz_be_0 + g_xzzzzz_yyyz_1[i] * pa_x[i];

        g_xxzzzzz_yyzz_0[i] = g_zzzzz_yyzz_0[i] * fbe_0 - g_zzzzz_yyzz_1[i] * fz_be_0 + g_xzzzzz_yyzz_1[i] * pa_x[i];

        g_xxzzzzz_yzzz_0[i] = g_zzzzz_yzzz_0[i] * fbe_0 - g_zzzzz_yzzz_1[i] * fz_be_0 + g_xzzzzz_yzzz_1[i] * pa_x[i];

        g_xxzzzzz_zzzz_0[i] = g_zzzzz_zzzz_0[i] * fbe_0 - g_zzzzz_zzzz_1[i] * fz_be_0 + g_xzzzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 315-330 components of targeted buffer : KG

    auto g_xyyyyyy_xxxx_0 = pbuffer.data(idx_eri_0_kg + 315);

    auto g_xyyyyyy_xxxy_0 = pbuffer.data(idx_eri_0_kg + 316);

    auto g_xyyyyyy_xxxz_0 = pbuffer.data(idx_eri_0_kg + 317);

    auto g_xyyyyyy_xxyy_0 = pbuffer.data(idx_eri_0_kg + 318);

    auto g_xyyyyyy_xxyz_0 = pbuffer.data(idx_eri_0_kg + 319);

    auto g_xyyyyyy_xxzz_0 = pbuffer.data(idx_eri_0_kg + 320);

    auto g_xyyyyyy_xyyy_0 = pbuffer.data(idx_eri_0_kg + 321);

    auto g_xyyyyyy_xyyz_0 = pbuffer.data(idx_eri_0_kg + 322);

    auto g_xyyyyyy_xyzz_0 = pbuffer.data(idx_eri_0_kg + 323);

    auto g_xyyyyyy_xzzz_0 = pbuffer.data(idx_eri_0_kg + 324);

    auto g_xyyyyyy_yyyy_0 = pbuffer.data(idx_eri_0_kg + 325);

    auto g_xyyyyyy_yyyz_0 = pbuffer.data(idx_eri_0_kg + 326);

    auto g_xyyyyyy_yyzz_0 = pbuffer.data(idx_eri_0_kg + 327);

    auto g_xyyyyyy_yzzz_0 = pbuffer.data(idx_eri_0_kg + 328);

    auto g_xyyyyyy_zzzz_0 = pbuffer.data(idx_eri_0_kg + 329);

    #pragma omp simd aligned(g_xyyyyyy_xxxx_0, g_xyyyyyy_xxxy_0, g_xyyyyyy_xxxz_0, g_xyyyyyy_xxyy_0, g_xyyyyyy_xxyz_0, g_xyyyyyy_xxzz_0, g_xyyyyyy_xyyy_0, g_xyyyyyy_xyyz_0, g_xyyyyyy_xyzz_0, g_xyyyyyy_xzzz_0, g_xyyyyyy_yyyy_0, g_xyyyyyy_yyyz_0, g_xyyyyyy_yyzz_0, g_xyyyyyy_yzzz_0, g_xyyyyyy_zzzz_0, g_yyyyyy_xxx_1, g_yyyyyy_xxxx_1, g_yyyyyy_xxxy_1, g_yyyyyy_xxxz_1, g_yyyyyy_xxy_1, g_yyyyyy_xxyy_1, g_yyyyyy_xxyz_1, g_yyyyyy_xxz_1, g_yyyyyy_xxzz_1, g_yyyyyy_xyy_1, g_yyyyyy_xyyy_1, g_yyyyyy_xyyz_1, g_yyyyyy_xyz_1, g_yyyyyy_xyzz_1, g_yyyyyy_xzz_1, g_yyyyyy_xzzz_1, g_yyyyyy_yyy_1, g_yyyyyy_yyyy_1, g_yyyyyy_yyyz_1, g_yyyyyy_yyz_1, g_yyyyyy_yyzz_1, g_yyyyyy_yzz_1, g_yyyyyy_yzzz_1, g_yyyyyy_zzz_1, g_yyyyyy_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyyy_xxxx_0[i] = 4.0 * g_yyyyyy_xxx_1[i] * fe_0 + g_yyyyyy_xxxx_1[i] * pa_x[i];

        g_xyyyyyy_xxxy_0[i] = 3.0 * g_yyyyyy_xxy_1[i] * fe_0 + g_yyyyyy_xxxy_1[i] * pa_x[i];

        g_xyyyyyy_xxxz_0[i] = 3.0 * g_yyyyyy_xxz_1[i] * fe_0 + g_yyyyyy_xxxz_1[i] * pa_x[i];

        g_xyyyyyy_xxyy_0[i] = 2.0 * g_yyyyyy_xyy_1[i] * fe_0 + g_yyyyyy_xxyy_1[i] * pa_x[i];

        g_xyyyyyy_xxyz_0[i] = 2.0 * g_yyyyyy_xyz_1[i] * fe_0 + g_yyyyyy_xxyz_1[i] * pa_x[i];

        g_xyyyyyy_xxzz_0[i] = 2.0 * g_yyyyyy_xzz_1[i] * fe_0 + g_yyyyyy_xxzz_1[i] * pa_x[i];

        g_xyyyyyy_xyyy_0[i] = g_yyyyyy_yyy_1[i] * fe_0 + g_yyyyyy_xyyy_1[i] * pa_x[i];

        g_xyyyyyy_xyyz_0[i] = g_yyyyyy_yyz_1[i] * fe_0 + g_yyyyyy_xyyz_1[i] * pa_x[i];

        g_xyyyyyy_xyzz_0[i] = g_yyyyyy_yzz_1[i] * fe_0 + g_yyyyyy_xyzz_1[i] * pa_x[i];

        g_xyyyyyy_xzzz_0[i] = g_yyyyyy_zzz_1[i] * fe_0 + g_yyyyyy_xzzz_1[i] * pa_x[i];

        g_xyyyyyy_yyyy_0[i] = g_yyyyyy_yyyy_1[i] * pa_x[i];

        g_xyyyyyy_yyyz_0[i] = g_yyyyyy_yyyz_1[i] * pa_x[i];

        g_xyyyyyy_yyzz_0[i] = g_yyyyyy_yyzz_1[i] * pa_x[i];

        g_xyyyyyy_yzzz_0[i] = g_yyyyyy_yzzz_1[i] * pa_x[i];

        g_xyyyyyy_zzzz_0[i] = g_yyyyyy_zzzz_1[i] * pa_x[i];
    }

    // Set up 330-345 components of targeted buffer : KG

    auto g_xyyyyyz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 330);

    auto g_xyyyyyz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 331);

    auto g_xyyyyyz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 332);

    auto g_xyyyyyz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 333);

    auto g_xyyyyyz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 334);

    auto g_xyyyyyz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 335);

    auto g_xyyyyyz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 336);

    auto g_xyyyyyz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 337);

    auto g_xyyyyyz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 338);

    auto g_xyyyyyz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 339);

    auto g_xyyyyyz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 340);

    auto g_xyyyyyz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 341);

    auto g_xyyyyyz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 342);

    auto g_xyyyyyz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 343);

    auto g_xyyyyyz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 344);

    #pragma omp simd aligned(g_xyyyyy_xxxx_1, g_xyyyyy_xxxy_1, g_xyyyyy_xxyy_1, g_xyyyyy_xyyy_1, g_xyyyyyz_xxxx_0, g_xyyyyyz_xxxy_0, g_xyyyyyz_xxxz_0, g_xyyyyyz_xxyy_0, g_xyyyyyz_xxyz_0, g_xyyyyyz_xxzz_0, g_xyyyyyz_xyyy_0, g_xyyyyyz_xyyz_0, g_xyyyyyz_xyzz_0, g_xyyyyyz_xzzz_0, g_xyyyyyz_yyyy_0, g_xyyyyyz_yyyz_0, g_xyyyyyz_yyzz_0, g_xyyyyyz_yzzz_0, g_xyyyyyz_zzzz_0, g_yyyyyz_xxxz_1, g_yyyyyz_xxyz_1, g_yyyyyz_xxz_1, g_yyyyyz_xxzz_1, g_yyyyyz_xyyz_1, g_yyyyyz_xyz_1, g_yyyyyz_xyzz_1, g_yyyyyz_xzz_1, g_yyyyyz_xzzz_1, g_yyyyyz_yyyy_1, g_yyyyyz_yyyz_1, g_yyyyyz_yyz_1, g_yyyyyz_yyzz_1, g_yyyyyz_yzz_1, g_yyyyyz_yzzz_1, g_yyyyyz_zzz_1, g_yyyyyz_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyyz_xxxx_0[i] = g_xyyyyy_xxxx_1[i] * pa_z[i];

        g_xyyyyyz_xxxy_0[i] = g_xyyyyy_xxxy_1[i] * pa_z[i];

        g_xyyyyyz_xxxz_0[i] = 3.0 * g_yyyyyz_xxz_1[i] * fe_0 + g_yyyyyz_xxxz_1[i] * pa_x[i];

        g_xyyyyyz_xxyy_0[i] = g_xyyyyy_xxyy_1[i] * pa_z[i];

        g_xyyyyyz_xxyz_0[i] = 2.0 * g_yyyyyz_xyz_1[i] * fe_0 + g_yyyyyz_xxyz_1[i] * pa_x[i];

        g_xyyyyyz_xxzz_0[i] = 2.0 * g_yyyyyz_xzz_1[i] * fe_0 + g_yyyyyz_xxzz_1[i] * pa_x[i];

        g_xyyyyyz_xyyy_0[i] = g_xyyyyy_xyyy_1[i] * pa_z[i];

        g_xyyyyyz_xyyz_0[i] = g_yyyyyz_yyz_1[i] * fe_0 + g_yyyyyz_xyyz_1[i] * pa_x[i];

        g_xyyyyyz_xyzz_0[i] = g_yyyyyz_yzz_1[i] * fe_0 + g_yyyyyz_xyzz_1[i] * pa_x[i];

        g_xyyyyyz_xzzz_0[i] = g_yyyyyz_zzz_1[i] * fe_0 + g_yyyyyz_xzzz_1[i] * pa_x[i];

        g_xyyyyyz_yyyy_0[i] = g_yyyyyz_yyyy_1[i] * pa_x[i];

        g_xyyyyyz_yyyz_0[i] = g_yyyyyz_yyyz_1[i] * pa_x[i];

        g_xyyyyyz_yyzz_0[i] = g_yyyyyz_yyzz_1[i] * pa_x[i];

        g_xyyyyyz_yzzz_0[i] = g_yyyyyz_yzzz_1[i] * pa_x[i];

        g_xyyyyyz_zzzz_0[i] = g_yyyyyz_zzzz_1[i] * pa_x[i];
    }

    // Set up 345-360 components of targeted buffer : KG

    auto g_xyyyyzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 345);

    auto g_xyyyyzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 346);

    auto g_xyyyyzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 347);

    auto g_xyyyyzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 348);

    auto g_xyyyyzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 349);

    auto g_xyyyyzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 350);

    auto g_xyyyyzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 351);

    auto g_xyyyyzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 352);

    auto g_xyyyyzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 353);

    auto g_xyyyyzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 354);

    auto g_xyyyyzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 355);

    auto g_xyyyyzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 356);

    auto g_xyyyyzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 357);

    auto g_xyyyyzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 358);

    auto g_xyyyyzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 359);

    #pragma omp simd aligned(g_xyyyyzz_xxxx_0, g_xyyyyzz_xxxy_0, g_xyyyyzz_xxxz_0, g_xyyyyzz_xxyy_0, g_xyyyyzz_xxyz_0, g_xyyyyzz_xxzz_0, g_xyyyyzz_xyyy_0, g_xyyyyzz_xyyz_0, g_xyyyyzz_xyzz_0, g_xyyyyzz_xzzz_0, g_xyyyyzz_yyyy_0, g_xyyyyzz_yyyz_0, g_xyyyyzz_yyzz_0, g_xyyyyzz_yzzz_0, g_xyyyyzz_zzzz_0, g_yyyyzz_xxx_1, g_yyyyzz_xxxx_1, g_yyyyzz_xxxy_1, g_yyyyzz_xxxz_1, g_yyyyzz_xxy_1, g_yyyyzz_xxyy_1, g_yyyyzz_xxyz_1, g_yyyyzz_xxz_1, g_yyyyzz_xxzz_1, g_yyyyzz_xyy_1, g_yyyyzz_xyyy_1, g_yyyyzz_xyyz_1, g_yyyyzz_xyz_1, g_yyyyzz_xyzz_1, g_yyyyzz_xzz_1, g_yyyyzz_xzzz_1, g_yyyyzz_yyy_1, g_yyyyzz_yyyy_1, g_yyyyzz_yyyz_1, g_yyyyzz_yyz_1, g_yyyyzz_yyzz_1, g_yyyyzz_yzz_1, g_yyyyzz_yzzz_1, g_yyyyzz_zzz_1, g_yyyyzz_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyzz_xxxx_0[i] = 4.0 * g_yyyyzz_xxx_1[i] * fe_0 + g_yyyyzz_xxxx_1[i] * pa_x[i];

        g_xyyyyzz_xxxy_0[i] = 3.0 * g_yyyyzz_xxy_1[i] * fe_0 + g_yyyyzz_xxxy_1[i] * pa_x[i];

        g_xyyyyzz_xxxz_0[i] = 3.0 * g_yyyyzz_xxz_1[i] * fe_0 + g_yyyyzz_xxxz_1[i] * pa_x[i];

        g_xyyyyzz_xxyy_0[i] = 2.0 * g_yyyyzz_xyy_1[i] * fe_0 + g_yyyyzz_xxyy_1[i] * pa_x[i];

        g_xyyyyzz_xxyz_0[i] = 2.0 * g_yyyyzz_xyz_1[i] * fe_0 + g_yyyyzz_xxyz_1[i] * pa_x[i];

        g_xyyyyzz_xxzz_0[i] = 2.0 * g_yyyyzz_xzz_1[i] * fe_0 + g_yyyyzz_xxzz_1[i] * pa_x[i];

        g_xyyyyzz_xyyy_0[i] = g_yyyyzz_yyy_1[i] * fe_0 + g_yyyyzz_xyyy_1[i] * pa_x[i];

        g_xyyyyzz_xyyz_0[i] = g_yyyyzz_yyz_1[i] * fe_0 + g_yyyyzz_xyyz_1[i] * pa_x[i];

        g_xyyyyzz_xyzz_0[i] = g_yyyyzz_yzz_1[i] * fe_0 + g_yyyyzz_xyzz_1[i] * pa_x[i];

        g_xyyyyzz_xzzz_0[i] = g_yyyyzz_zzz_1[i] * fe_0 + g_yyyyzz_xzzz_1[i] * pa_x[i];

        g_xyyyyzz_yyyy_0[i] = g_yyyyzz_yyyy_1[i] * pa_x[i];

        g_xyyyyzz_yyyz_0[i] = g_yyyyzz_yyyz_1[i] * pa_x[i];

        g_xyyyyzz_yyzz_0[i] = g_yyyyzz_yyzz_1[i] * pa_x[i];

        g_xyyyyzz_yzzz_0[i] = g_yyyyzz_yzzz_1[i] * pa_x[i];

        g_xyyyyzz_zzzz_0[i] = g_yyyyzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 360-375 components of targeted buffer : KG

    auto g_xyyyzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 360);

    auto g_xyyyzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 361);

    auto g_xyyyzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 362);

    auto g_xyyyzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 363);

    auto g_xyyyzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 364);

    auto g_xyyyzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 365);

    auto g_xyyyzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 366);

    auto g_xyyyzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 367);

    auto g_xyyyzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 368);

    auto g_xyyyzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 369);

    auto g_xyyyzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 370);

    auto g_xyyyzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 371);

    auto g_xyyyzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 372);

    auto g_xyyyzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 373);

    auto g_xyyyzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 374);

    #pragma omp simd aligned(g_xyyyzzz_xxxx_0, g_xyyyzzz_xxxy_0, g_xyyyzzz_xxxz_0, g_xyyyzzz_xxyy_0, g_xyyyzzz_xxyz_0, g_xyyyzzz_xxzz_0, g_xyyyzzz_xyyy_0, g_xyyyzzz_xyyz_0, g_xyyyzzz_xyzz_0, g_xyyyzzz_xzzz_0, g_xyyyzzz_yyyy_0, g_xyyyzzz_yyyz_0, g_xyyyzzz_yyzz_0, g_xyyyzzz_yzzz_0, g_xyyyzzz_zzzz_0, g_yyyzzz_xxx_1, g_yyyzzz_xxxx_1, g_yyyzzz_xxxy_1, g_yyyzzz_xxxz_1, g_yyyzzz_xxy_1, g_yyyzzz_xxyy_1, g_yyyzzz_xxyz_1, g_yyyzzz_xxz_1, g_yyyzzz_xxzz_1, g_yyyzzz_xyy_1, g_yyyzzz_xyyy_1, g_yyyzzz_xyyz_1, g_yyyzzz_xyz_1, g_yyyzzz_xyzz_1, g_yyyzzz_xzz_1, g_yyyzzz_xzzz_1, g_yyyzzz_yyy_1, g_yyyzzz_yyyy_1, g_yyyzzz_yyyz_1, g_yyyzzz_yyz_1, g_yyyzzz_yyzz_1, g_yyyzzz_yzz_1, g_yyyzzz_yzzz_1, g_yyyzzz_zzz_1, g_yyyzzz_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyzzz_xxxx_0[i] = 4.0 * g_yyyzzz_xxx_1[i] * fe_0 + g_yyyzzz_xxxx_1[i] * pa_x[i];

        g_xyyyzzz_xxxy_0[i] = 3.0 * g_yyyzzz_xxy_1[i] * fe_0 + g_yyyzzz_xxxy_1[i] * pa_x[i];

        g_xyyyzzz_xxxz_0[i] = 3.0 * g_yyyzzz_xxz_1[i] * fe_0 + g_yyyzzz_xxxz_1[i] * pa_x[i];

        g_xyyyzzz_xxyy_0[i] = 2.0 * g_yyyzzz_xyy_1[i] * fe_0 + g_yyyzzz_xxyy_1[i] * pa_x[i];

        g_xyyyzzz_xxyz_0[i] = 2.0 * g_yyyzzz_xyz_1[i] * fe_0 + g_yyyzzz_xxyz_1[i] * pa_x[i];

        g_xyyyzzz_xxzz_0[i] = 2.0 * g_yyyzzz_xzz_1[i] * fe_0 + g_yyyzzz_xxzz_1[i] * pa_x[i];

        g_xyyyzzz_xyyy_0[i] = g_yyyzzz_yyy_1[i] * fe_0 + g_yyyzzz_xyyy_1[i] * pa_x[i];

        g_xyyyzzz_xyyz_0[i] = g_yyyzzz_yyz_1[i] * fe_0 + g_yyyzzz_xyyz_1[i] * pa_x[i];

        g_xyyyzzz_xyzz_0[i] = g_yyyzzz_yzz_1[i] * fe_0 + g_yyyzzz_xyzz_1[i] * pa_x[i];

        g_xyyyzzz_xzzz_0[i] = g_yyyzzz_zzz_1[i] * fe_0 + g_yyyzzz_xzzz_1[i] * pa_x[i];

        g_xyyyzzz_yyyy_0[i] = g_yyyzzz_yyyy_1[i] * pa_x[i];

        g_xyyyzzz_yyyz_0[i] = g_yyyzzz_yyyz_1[i] * pa_x[i];

        g_xyyyzzz_yyzz_0[i] = g_yyyzzz_yyzz_1[i] * pa_x[i];

        g_xyyyzzz_yzzz_0[i] = g_yyyzzz_yzzz_1[i] * pa_x[i];

        g_xyyyzzz_zzzz_0[i] = g_yyyzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 375-390 components of targeted buffer : KG

    auto g_xyyzzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 375);

    auto g_xyyzzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 376);

    auto g_xyyzzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 377);

    auto g_xyyzzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 378);

    auto g_xyyzzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 379);

    auto g_xyyzzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 380);

    auto g_xyyzzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 381);

    auto g_xyyzzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 382);

    auto g_xyyzzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 383);

    auto g_xyyzzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 384);

    auto g_xyyzzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 385);

    auto g_xyyzzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 386);

    auto g_xyyzzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 387);

    auto g_xyyzzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 388);

    auto g_xyyzzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 389);

    #pragma omp simd aligned(g_xyyzzzz_xxxx_0, g_xyyzzzz_xxxy_0, g_xyyzzzz_xxxz_0, g_xyyzzzz_xxyy_0, g_xyyzzzz_xxyz_0, g_xyyzzzz_xxzz_0, g_xyyzzzz_xyyy_0, g_xyyzzzz_xyyz_0, g_xyyzzzz_xyzz_0, g_xyyzzzz_xzzz_0, g_xyyzzzz_yyyy_0, g_xyyzzzz_yyyz_0, g_xyyzzzz_yyzz_0, g_xyyzzzz_yzzz_0, g_xyyzzzz_zzzz_0, g_yyzzzz_xxx_1, g_yyzzzz_xxxx_1, g_yyzzzz_xxxy_1, g_yyzzzz_xxxz_1, g_yyzzzz_xxy_1, g_yyzzzz_xxyy_1, g_yyzzzz_xxyz_1, g_yyzzzz_xxz_1, g_yyzzzz_xxzz_1, g_yyzzzz_xyy_1, g_yyzzzz_xyyy_1, g_yyzzzz_xyyz_1, g_yyzzzz_xyz_1, g_yyzzzz_xyzz_1, g_yyzzzz_xzz_1, g_yyzzzz_xzzz_1, g_yyzzzz_yyy_1, g_yyzzzz_yyyy_1, g_yyzzzz_yyyz_1, g_yyzzzz_yyz_1, g_yyzzzz_yyzz_1, g_yyzzzz_yzz_1, g_yyzzzz_yzzz_1, g_yyzzzz_zzz_1, g_yyzzzz_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzzzz_xxxx_0[i] = 4.0 * g_yyzzzz_xxx_1[i] * fe_0 + g_yyzzzz_xxxx_1[i] * pa_x[i];

        g_xyyzzzz_xxxy_0[i] = 3.0 * g_yyzzzz_xxy_1[i] * fe_0 + g_yyzzzz_xxxy_1[i] * pa_x[i];

        g_xyyzzzz_xxxz_0[i] = 3.0 * g_yyzzzz_xxz_1[i] * fe_0 + g_yyzzzz_xxxz_1[i] * pa_x[i];

        g_xyyzzzz_xxyy_0[i] = 2.0 * g_yyzzzz_xyy_1[i] * fe_0 + g_yyzzzz_xxyy_1[i] * pa_x[i];

        g_xyyzzzz_xxyz_0[i] = 2.0 * g_yyzzzz_xyz_1[i] * fe_0 + g_yyzzzz_xxyz_1[i] * pa_x[i];

        g_xyyzzzz_xxzz_0[i] = 2.0 * g_yyzzzz_xzz_1[i] * fe_0 + g_yyzzzz_xxzz_1[i] * pa_x[i];

        g_xyyzzzz_xyyy_0[i] = g_yyzzzz_yyy_1[i] * fe_0 + g_yyzzzz_xyyy_1[i] * pa_x[i];

        g_xyyzzzz_xyyz_0[i] = g_yyzzzz_yyz_1[i] * fe_0 + g_yyzzzz_xyyz_1[i] * pa_x[i];

        g_xyyzzzz_xyzz_0[i] = g_yyzzzz_yzz_1[i] * fe_0 + g_yyzzzz_xyzz_1[i] * pa_x[i];

        g_xyyzzzz_xzzz_0[i] = g_yyzzzz_zzz_1[i] * fe_0 + g_yyzzzz_xzzz_1[i] * pa_x[i];

        g_xyyzzzz_yyyy_0[i] = g_yyzzzz_yyyy_1[i] * pa_x[i];

        g_xyyzzzz_yyyz_0[i] = g_yyzzzz_yyyz_1[i] * pa_x[i];

        g_xyyzzzz_yyzz_0[i] = g_yyzzzz_yyzz_1[i] * pa_x[i];

        g_xyyzzzz_yzzz_0[i] = g_yyzzzz_yzzz_1[i] * pa_x[i];

        g_xyyzzzz_zzzz_0[i] = g_yyzzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 390-405 components of targeted buffer : KG

    auto g_xyzzzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 390);

    auto g_xyzzzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 391);

    auto g_xyzzzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 392);

    auto g_xyzzzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 393);

    auto g_xyzzzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 394);

    auto g_xyzzzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 395);

    auto g_xyzzzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 396);

    auto g_xyzzzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 397);

    auto g_xyzzzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 398);

    auto g_xyzzzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 399);

    auto g_xyzzzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 400);

    auto g_xyzzzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 401);

    auto g_xyzzzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 402);

    auto g_xyzzzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 403);

    auto g_xyzzzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 404);

    #pragma omp simd aligned(g_xyzzzzz_xxxx_0, g_xyzzzzz_xxxy_0, g_xyzzzzz_xxxz_0, g_xyzzzzz_xxyy_0, g_xyzzzzz_xxyz_0, g_xyzzzzz_xxzz_0, g_xyzzzzz_xyyy_0, g_xyzzzzz_xyyz_0, g_xyzzzzz_xyzz_0, g_xyzzzzz_xzzz_0, g_xyzzzzz_yyyy_0, g_xyzzzzz_yyyz_0, g_xyzzzzz_yyzz_0, g_xyzzzzz_yzzz_0, g_xyzzzzz_zzzz_0, g_xzzzzz_xxxx_1, g_xzzzzz_xxxz_1, g_xzzzzz_xxzz_1, g_xzzzzz_xzzz_1, g_yzzzzz_xxxy_1, g_yzzzzz_xxy_1, g_yzzzzz_xxyy_1, g_yzzzzz_xxyz_1, g_yzzzzz_xyy_1, g_yzzzzz_xyyy_1, g_yzzzzz_xyyz_1, g_yzzzzz_xyz_1, g_yzzzzz_xyzz_1, g_yzzzzz_yyy_1, g_yzzzzz_yyyy_1, g_yzzzzz_yyyz_1, g_yzzzzz_yyz_1, g_yzzzzz_yyzz_1, g_yzzzzz_yzz_1, g_yzzzzz_yzzz_1, g_yzzzzz_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzzzz_xxxx_0[i] = g_xzzzzz_xxxx_1[i] * pa_y[i];

        g_xyzzzzz_xxxy_0[i] = 3.0 * g_yzzzzz_xxy_1[i] * fe_0 + g_yzzzzz_xxxy_1[i] * pa_x[i];

        g_xyzzzzz_xxxz_0[i] = g_xzzzzz_xxxz_1[i] * pa_y[i];

        g_xyzzzzz_xxyy_0[i] = 2.0 * g_yzzzzz_xyy_1[i] * fe_0 + g_yzzzzz_xxyy_1[i] * pa_x[i];

        g_xyzzzzz_xxyz_0[i] = 2.0 * g_yzzzzz_xyz_1[i] * fe_0 + g_yzzzzz_xxyz_1[i] * pa_x[i];

        g_xyzzzzz_xxzz_0[i] = g_xzzzzz_xxzz_1[i] * pa_y[i];

        g_xyzzzzz_xyyy_0[i] = g_yzzzzz_yyy_1[i] * fe_0 + g_yzzzzz_xyyy_1[i] * pa_x[i];

        g_xyzzzzz_xyyz_0[i] = g_yzzzzz_yyz_1[i] * fe_0 + g_yzzzzz_xyyz_1[i] * pa_x[i];

        g_xyzzzzz_xyzz_0[i] = g_yzzzzz_yzz_1[i] * fe_0 + g_yzzzzz_xyzz_1[i] * pa_x[i];

        g_xyzzzzz_xzzz_0[i] = g_xzzzzz_xzzz_1[i] * pa_y[i];

        g_xyzzzzz_yyyy_0[i] = g_yzzzzz_yyyy_1[i] * pa_x[i];

        g_xyzzzzz_yyyz_0[i] = g_yzzzzz_yyyz_1[i] * pa_x[i];

        g_xyzzzzz_yyzz_0[i] = g_yzzzzz_yyzz_1[i] * pa_x[i];

        g_xyzzzzz_yzzz_0[i] = g_yzzzzz_yzzz_1[i] * pa_x[i];

        g_xyzzzzz_zzzz_0[i] = g_yzzzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 405-420 components of targeted buffer : KG

    auto g_xzzzzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 405);

    auto g_xzzzzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 406);

    auto g_xzzzzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 407);

    auto g_xzzzzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 408);

    auto g_xzzzzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 409);

    auto g_xzzzzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 410);

    auto g_xzzzzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 411);

    auto g_xzzzzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 412);

    auto g_xzzzzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 413);

    auto g_xzzzzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 414);

    auto g_xzzzzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 415);

    auto g_xzzzzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 416);

    auto g_xzzzzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 417);

    auto g_xzzzzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 418);

    auto g_xzzzzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 419);

    #pragma omp simd aligned(g_xzzzzzz_xxxx_0, g_xzzzzzz_xxxy_0, g_xzzzzzz_xxxz_0, g_xzzzzzz_xxyy_0, g_xzzzzzz_xxyz_0, g_xzzzzzz_xxzz_0, g_xzzzzzz_xyyy_0, g_xzzzzzz_xyyz_0, g_xzzzzzz_xyzz_0, g_xzzzzzz_xzzz_0, g_xzzzzzz_yyyy_0, g_xzzzzzz_yyyz_0, g_xzzzzzz_yyzz_0, g_xzzzzzz_yzzz_0, g_xzzzzzz_zzzz_0, g_zzzzzz_xxx_1, g_zzzzzz_xxxx_1, g_zzzzzz_xxxy_1, g_zzzzzz_xxxz_1, g_zzzzzz_xxy_1, g_zzzzzz_xxyy_1, g_zzzzzz_xxyz_1, g_zzzzzz_xxz_1, g_zzzzzz_xxzz_1, g_zzzzzz_xyy_1, g_zzzzzz_xyyy_1, g_zzzzzz_xyyz_1, g_zzzzzz_xyz_1, g_zzzzzz_xyzz_1, g_zzzzzz_xzz_1, g_zzzzzz_xzzz_1, g_zzzzzz_yyy_1, g_zzzzzz_yyyy_1, g_zzzzzz_yyyz_1, g_zzzzzz_yyz_1, g_zzzzzz_yyzz_1, g_zzzzzz_yzz_1, g_zzzzzz_yzzz_1, g_zzzzzz_zzz_1, g_zzzzzz_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzzzz_xxxx_0[i] = 4.0 * g_zzzzzz_xxx_1[i] * fe_0 + g_zzzzzz_xxxx_1[i] * pa_x[i];

        g_xzzzzzz_xxxy_0[i] = 3.0 * g_zzzzzz_xxy_1[i] * fe_0 + g_zzzzzz_xxxy_1[i] * pa_x[i];

        g_xzzzzzz_xxxz_0[i] = 3.0 * g_zzzzzz_xxz_1[i] * fe_0 + g_zzzzzz_xxxz_1[i] * pa_x[i];

        g_xzzzzzz_xxyy_0[i] = 2.0 * g_zzzzzz_xyy_1[i] * fe_0 + g_zzzzzz_xxyy_1[i] * pa_x[i];

        g_xzzzzzz_xxyz_0[i] = 2.0 * g_zzzzzz_xyz_1[i] * fe_0 + g_zzzzzz_xxyz_1[i] * pa_x[i];

        g_xzzzzzz_xxzz_0[i] = 2.0 * g_zzzzzz_xzz_1[i] * fe_0 + g_zzzzzz_xxzz_1[i] * pa_x[i];

        g_xzzzzzz_xyyy_0[i] = g_zzzzzz_yyy_1[i] * fe_0 + g_zzzzzz_xyyy_1[i] * pa_x[i];

        g_xzzzzzz_xyyz_0[i] = g_zzzzzz_yyz_1[i] * fe_0 + g_zzzzzz_xyyz_1[i] * pa_x[i];

        g_xzzzzzz_xyzz_0[i] = g_zzzzzz_yzz_1[i] * fe_0 + g_zzzzzz_xyzz_1[i] * pa_x[i];

        g_xzzzzzz_xzzz_0[i] = g_zzzzzz_zzz_1[i] * fe_0 + g_zzzzzz_xzzz_1[i] * pa_x[i];

        g_xzzzzzz_yyyy_0[i] = g_zzzzzz_yyyy_1[i] * pa_x[i];

        g_xzzzzzz_yyyz_0[i] = g_zzzzzz_yyyz_1[i] * pa_x[i];

        g_xzzzzzz_yyzz_0[i] = g_zzzzzz_yyzz_1[i] * pa_x[i];

        g_xzzzzzz_yzzz_0[i] = g_zzzzzz_yzzz_1[i] * pa_x[i];

        g_xzzzzzz_zzzz_0[i] = g_zzzzzz_zzzz_1[i] * pa_x[i];
    }

    // Set up 420-435 components of targeted buffer : KG

    auto g_yyyyyyy_xxxx_0 = pbuffer.data(idx_eri_0_kg + 420);

    auto g_yyyyyyy_xxxy_0 = pbuffer.data(idx_eri_0_kg + 421);

    auto g_yyyyyyy_xxxz_0 = pbuffer.data(idx_eri_0_kg + 422);

    auto g_yyyyyyy_xxyy_0 = pbuffer.data(idx_eri_0_kg + 423);

    auto g_yyyyyyy_xxyz_0 = pbuffer.data(idx_eri_0_kg + 424);

    auto g_yyyyyyy_xxzz_0 = pbuffer.data(idx_eri_0_kg + 425);

    auto g_yyyyyyy_xyyy_0 = pbuffer.data(idx_eri_0_kg + 426);

    auto g_yyyyyyy_xyyz_0 = pbuffer.data(idx_eri_0_kg + 427);

    auto g_yyyyyyy_xyzz_0 = pbuffer.data(idx_eri_0_kg + 428);

    auto g_yyyyyyy_xzzz_0 = pbuffer.data(idx_eri_0_kg + 429);

    auto g_yyyyyyy_yyyy_0 = pbuffer.data(idx_eri_0_kg + 430);

    auto g_yyyyyyy_yyyz_0 = pbuffer.data(idx_eri_0_kg + 431);

    auto g_yyyyyyy_yyzz_0 = pbuffer.data(idx_eri_0_kg + 432);

    auto g_yyyyyyy_yzzz_0 = pbuffer.data(idx_eri_0_kg + 433);

    auto g_yyyyyyy_zzzz_0 = pbuffer.data(idx_eri_0_kg + 434);

    #pragma omp simd aligned(g_yyyyy_xxxx_0, g_yyyyy_xxxx_1, g_yyyyy_xxxy_0, g_yyyyy_xxxy_1, g_yyyyy_xxxz_0, g_yyyyy_xxxz_1, g_yyyyy_xxyy_0, g_yyyyy_xxyy_1, g_yyyyy_xxyz_0, g_yyyyy_xxyz_1, g_yyyyy_xxzz_0, g_yyyyy_xxzz_1, g_yyyyy_xyyy_0, g_yyyyy_xyyy_1, g_yyyyy_xyyz_0, g_yyyyy_xyyz_1, g_yyyyy_xyzz_0, g_yyyyy_xyzz_1, g_yyyyy_xzzz_0, g_yyyyy_xzzz_1, g_yyyyy_yyyy_0, g_yyyyy_yyyy_1, g_yyyyy_yyyz_0, g_yyyyy_yyyz_1, g_yyyyy_yyzz_0, g_yyyyy_yyzz_1, g_yyyyy_yzzz_0, g_yyyyy_yzzz_1, g_yyyyy_zzzz_0, g_yyyyy_zzzz_1, g_yyyyyy_xxx_1, g_yyyyyy_xxxx_1, g_yyyyyy_xxxy_1, g_yyyyyy_xxxz_1, g_yyyyyy_xxy_1, g_yyyyyy_xxyy_1, g_yyyyyy_xxyz_1, g_yyyyyy_xxz_1, g_yyyyyy_xxzz_1, g_yyyyyy_xyy_1, g_yyyyyy_xyyy_1, g_yyyyyy_xyyz_1, g_yyyyyy_xyz_1, g_yyyyyy_xyzz_1, g_yyyyyy_xzz_1, g_yyyyyy_xzzz_1, g_yyyyyy_yyy_1, g_yyyyyy_yyyy_1, g_yyyyyy_yyyz_1, g_yyyyyy_yyz_1, g_yyyyyy_yyzz_1, g_yyyyyy_yzz_1, g_yyyyyy_yzzz_1, g_yyyyyy_zzz_1, g_yyyyyy_zzzz_1, g_yyyyyyy_xxxx_0, g_yyyyyyy_xxxy_0, g_yyyyyyy_xxxz_0, g_yyyyyyy_xxyy_0, g_yyyyyyy_xxyz_0, g_yyyyyyy_xxzz_0, g_yyyyyyy_xyyy_0, g_yyyyyyy_xyyz_0, g_yyyyyyy_xyzz_0, g_yyyyyyy_xzzz_0, g_yyyyyyy_yyyy_0, g_yyyyyyy_yyyz_0, g_yyyyyyy_yyzz_0, g_yyyyyyy_yzzz_0, g_yyyyyyy_zzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyyy_xxxx_0[i] = 6.0 * g_yyyyy_xxxx_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxx_1[i] * fz_be_0 + g_yyyyyy_xxxx_1[i] * pa_y[i];

        g_yyyyyyy_xxxy_0[i] = 6.0 * g_yyyyy_xxxy_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxy_1[i] * fz_be_0 + g_yyyyyy_xxx_1[i] * fe_0 + g_yyyyyy_xxxy_1[i] * pa_y[i];

        g_yyyyyyy_xxxz_0[i] = 6.0 * g_yyyyy_xxxz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxz_1[i] * fz_be_0 + g_yyyyyy_xxxz_1[i] * pa_y[i];

        g_yyyyyyy_xxyy_0[i] = 6.0 * g_yyyyy_xxyy_0[i] * fbe_0 - 6.0 * g_yyyyy_xxyy_1[i] * fz_be_0 + 2.0 * g_yyyyyy_xxy_1[i] * fe_0 + g_yyyyyy_xxyy_1[i] * pa_y[i];

        g_yyyyyyy_xxyz_0[i] = 6.0 * g_yyyyy_xxyz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxyz_1[i] * fz_be_0 + g_yyyyyy_xxz_1[i] * fe_0 + g_yyyyyy_xxyz_1[i] * pa_y[i];

        g_yyyyyyy_xxzz_0[i] = 6.0 * g_yyyyy_xxzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxzz_1[i] * fz_be_0 + g_yyyyyy_xxzz_1[i] * pa_y[i];

        g_yyyyyyy_xyyy_0[i] = 6.0 * g_yyyyy_xyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_xyyy_1[i] * fz_be_0 + 3.0 * g_yyyyyy_xyy_1[i] * fe_0 + g_yyyyyy_xyyy_1[i] * pa_y[i];

        g_yyyyyyy_xyyz_0[i] = 6.0 * g_yyyyy_xyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_xyyz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_xyz_1[i] * fe_0 + g_yyyyyy_xyyz_1[i] * pa_y[i];

        g_yyyyyyy_xyzz_0[i] = 6.0 * g_yyyyy_xyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xyzz_1[i] * fz_be_0 + g_yyyyyy_xzz_1[i] * fe_0 + g_yyyyyy_xyzz_1[i] * pa_y[i];

        g_yyyyyyy_xzzz_0[i] = 6.0 * g_yyyyy_xzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xzzz_1[i] * fz_be_0 + g_yyyyyy_xzzz_1[i] * pa_y[i];

        g_yyyyyyy_yyyy_0[i] = 6.0 * g_yyyyy_yyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_yyyy_1[i] * fz_be_0 + 4.0 * g_yyyyyy_yyy_1[i] * fe_0 + g_yyyyyy_yyyy_1[i] * pa_y[i];

        g_yyyyyyy_yyyz_0[i] = 6.0 * g_yyyyy_yyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_yyyz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_yyz_1[i] * fe_0 + g_yyyyyy_yyyz_1[i] * pa_y[i];

        g_yyyyyyy_yyzz_0[i] = 6.0 * g_yyyyy_yyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_yyzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_yzz_1[i] * fe_0 + g_yyyyyy_yyzz_1[i] * pa_y[i];

        g_yyyyyyy_yzzz_0[i] = 6.0 * g_yyyyy_yzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_yzzz_1[i] * fz_be_0 + g_yyyyyy_zzz_1[i] * fe_0 + g_yyyyyy_yzzz_1[i] * pa_y[i];

        g_yyyyyyy_zzzz_0[i] = 6.0 * g_yyyyy_zzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_zzzz_1[i] * fz_be_0 + g_yyyyyy_zzzz_1[i] * pa_y[i];
    }

    // Set up 435-450 components of targeted buffer : KG

    auto g_yyyyyyz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 435);

    auto g_yyyyyyz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 436);

    auto g_yyyyyyz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 437);

    auto g_yyyyyyz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 438);

    auto g_yyyyyyz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 439);

    auto g_yyyyyyz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 440);

    auto g_yyyyyyz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 441);

    auto g_yyyyyyz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 442);

    auto g_yyyyyyz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 443);

    auto g_yyyyyyz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 444);

    auto g_yyyyyyz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 445);

    auto g_yyyyyyz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 446);

    auto g_yyyyyyz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 447);

    auto g_yyyyyyz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 448);

    auto g_yyyyyyz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 449);

    #pragma omp simd aligned(g_yyyyyy_xxx_1, g_yyyyyy_xxxx_1, g_yyyyyy_xxxy_1, g_yyyyyy_xxxz_1, g_yyyyyy_xxy_1, g_yyyyyy_xxyy_1, g_yyyyyy_xxyz_1, g_yyyyyy_xxz_1, g_yyyyyy_xxzz_1, g_yyyyyy_xyy_1, g_yyyyyy_xyyy_1, g_yyyyyy_xyyz_1, g_yyyyyy_xyz_1, g_yyyyyy_xyzz_1, g_yyyyyy_xzz_1, g_yyyyyy_xzzz_1, g_yyyyyy_yyy_1, g_yyyyyy_yyyy_1, g_yyyyyy_yyyz_1, g_yyyyyy_yyz_1, g_yyyyyy_yyzz_1, g_yyyyyy_yzz_1, g_yyyyyy_yzzz_1, g_yyyyyy_zzz_1, g_yyyyyy_zzzz_1, g_yyyyyyz_xxxx_0, g_yyyyyyz_xxxy_0, g_yyyyyyz_xxxz_0, g_yyyyyyz_xxyy_0, g_yyyyyyz_xxyz_0, g_yyyyyyz_xxzz_0, g_yyyyyyz_xyyy_0, g_yyyyyyz_xyyz_0, g_yyyyyyz_xyzz_0, g_yyyyyyz_xzzz_0, g_yyyyyyz_yyyy_0, g_yyyyyyz_yyyz_0, g_yyyyyyz_yyzz_0, g_yyyyyyz_yzzz_0, g_yyyyyyz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyyyz_xxxx_0[i] = g_yyyyyy_xxxx_1[i] * pa_z[i];

        g_yyyyyyz_xxxy_0[i] = g_yyyyyy_xxxy_1[i] * pa_z[i];

        g_yyyyyyz_xxxz_0[i] = g_yyyyyy_xxx_1[i] * fe_0 + g_yyyyyy_xxxz_1[i] * pa_z[i];

        g_yyyyyyz_xxyy_0[i] = g_yyyyyy_xxyy_1[i] * pa_z[i];

        g_yyyyyyz_xxyz_0[i] = g_yyyyyy_xxy_1[i] * fe_0 + g_yyyyyy_xxyz_1[i] * pa_z[i];

        g_yyyyyyz_xxzz_0[i] = 2.0 * g_yyyyyy_xxz_1[i] * fe_0 + g_yyyyyy_xxzz_1[i] * pa_z[i];

        g_yyyyyyz_xyyy_0[i] = g_yyyyyy_xyyy_1[i] * pa_z[i];

        g_yyyyyyz_xyyz_0[i] = g_yyyyyy_xyy_1[i] * fe_0 + g_yyyyyy_xyyz_1[i] * pa_z[i];

        g_yyyyyyz_xyzz_0[i] = 2.0 * g_yyyyyy_xyz_1[i] * fe_0 + g_yyyyyy_xyzz_1[i] * pa_z[i];

        g_yyyyyyz_xzzz_0[i] = 3.0 * g_yyyyyy_xzz_1[i] * fe_0 + g_yyyyyy_xzzz_1[i] * pa_z[i];

        g_yyyyyyz_yyyy_0[i] = g_yyyyyy_yyyy_1[i] * pa_z[i];

        g_yyyyyyz_yyyz_0[i] = g_yyyyyy_yyy_1[i] * fe_0 + g_yyyyyy_yyyz_1[i] * pa_z[i];

        g_yyyyyyz_yyzz_0[i] = 2.0 * g_yyyyyy_yyz_1[i] * fe_0 + g_yyyyyy_yyzz_1[i] * pa_z[i];

        g_yyyyyyz_yzzz_0[i] = 3.0 * g_yyyyyy_yzz_1[i] * fe_0 + g_yyyyyy_yzzz_1[i] * pa_z[i];

        g_yyyyyyz_zzzz_0[i] = 4.0 * g_yyyyyy_zzz_1[i] * fe_0 + g_yyyyyy_zzzz_1[i] * pa_z[i];
    }

    // Set up 450-465 components of targeted buffer : KG

    auto g_yyyyyzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 450);

    auto g_yyyyyzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 451);

    auto g_yyyyyzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 452);

    auto g_yyyyyzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 453);

    auto g_yyyyyzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 454);

    auto g_yyyyyzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 455);

    auto g_yyyyyzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 456);

    auto g_yyyyyzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 457);

    auto g_yyyyyzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 458);

    auto g_yyyyyzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 459);

    auto g_yyyyyzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 460);

    auto g_yyyyyzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 461);

    auto g_yyyyyzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 462);

    auto g_yyyyyzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 463);

    auto g_yyyyyzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 464);

    #pragma omp simd aligned(g_yyyyy_xxxy_0, g_yyyyy_xxxy_1, g_yyyyy_xxyy_0, g_yyyyy_xxyy_1, g_yyyyy_xyyy_0, g_yyyyy_xyyy_1, g_yyyyy_yyyy_0, g_yyyyy_yyyy_1, g_yyyyyz_xxxy_1, g_yyyyyz_xxyy_1, g_yyyyyz_xyyy_1, g_yyyyyz_yyyy_1, g_yyyyyzz_xxxx_0, g_yyyyyzz_xxxy_0, g_yyyyyzz_xxxz_0, g_yyyyyzz_xxyy_0, g_yyyyyzz_xxyz_0, g_yyyyyzz_xxzz_0, g_yyyyyzz_xyyy_0, g_yyyyyzz_xyyz_0, g_yyyyyzz_xyzz_0, g_yyyyyzz_xzzz_0, g_yyyyyzz_yyyy_0, g_yyyyyzz_yyyz_0, g_yyyyyzz_yyzz_0, g_yyyyyzz_yzzz_0, g_yyyyyzz_zzzz_0, g_yyyyzz_xxxx_1, g_yyyyzz_xxxz_1, g_yyyyzz_xxyz_1, g_yyyyzz_xxz_1, g_yyyyzz_xxzz_1, g_yyyyzz_xyyz_1, g_yyyyzz_xyz_1, g_yyyyzz_xyzz_1, g_yyyyzz_xzz_1, g_yyyyzz_xzzz_1, g_yyyyzz_yyyz_1, g_yyyyzz_yyz_1, g_yyyyzz_yyzz_1, g_yyyyzz_yzz_1, g_yyyyzz_yzzz_1, g_yyyyzz_zzz_1, g_yyyyzz_zzzz_1, g_yyyzz_xxxx_0, g_yyyzz_xxxx_1, g_yyyzz_xxxz_0, g_yyyzz_xxxz_1, g_yyyzz_xxyz_0, g_yyyzz_xxyz_1, g_yyyzz_xxzz_0, g_yyyzz_xxzz_1, g_yyyzz_xyyz_0, g_yyyzz_xyyz_1, g_yyyzz_xyzz_0, g_yyyzz_xyzz_1, g_yyyzz_xzzz_0, g_yyyzz_xzzz_1, g_yyyzz_yyyz_0, g_yyyzz_yyyz_1, g_yyyzz_yyzz_0, g_yyyzz_yyzz_1, g_yyyzz_yzzz_0, g_yyyzz_yzzz_1, g_yyyzz_zzzz_0, g_yyyzz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyzz_xxxx_0[i] = 4.0 * g_yyyzz_xxxx_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxx_1[i] * fz_be_0 + g_yyyyzz_xxxx_1[i] * pa_y[i];

        g_yyyyyzz_xxxy_0[i] = g_yyyyy_xxxy_0[i] * fbe_0 - g_yyyyy_xxxy_1[i] * fz_be_0 + g_yyyyyz_xxxy_1[i] * pa_z[i];

        g_yyyyyzz_xxxz_0[i] = 4.0 * g_yyyzz_xxxz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxz_1[i] * fz_be_0 + g_yyyyzz_xxxz_1[i] * pa_y[i];

        g_yyyyyzz_xxyy_0[i] = g_yyyyy_xxyy_0[i] * fbe_0 - g_yyyyy_xxyy_1[i] * fz_be_0 + g_yyyyyz_xxyy_1[i] * pa_z[i];

        g_yyyyyzz_xxyz_0[i] = 4.0 * g_yyyzz_xxyz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxyz_1[i] * fz_be_0 + g_yyyyzz_xxz_1[i] * fe_0 + g_yyyyzz_xxyz_1[i] * pa_y[i];

        g_yyyyyzz_xxzz_0[i] = 4.0 * g_yyyzz_xxzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxzz_1[i] * fz_be_0 + g_yyyyzz_xxzz_1[i] * pa_y[i];

        g_yyyyyzz_xyyy_0[i] = g_yyyyy_xyyy_0[i] * fbe_0 - g_yyyyy_xyyy_1[i] * fz_be_0 + g_yyyyyz_xyyy_1[i] * pa_z[i];

        g_yyyyyzz_xyyz_0[i] = 4.0 * g_yyyzz_xyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_xyyz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_xyz_1[i] * fe_0 + g_yyyyzz_xyyz_1[i] * pa_y[i];

        g_yyyyyzz_xyzz_0[i] = 4.0 * g_yyyzz_xyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xyzz_1[i] * fz_be_0 + g_yyyyzz_xzz_1[i] * fe_0 + g_yyyyzz_xyzz_1[i] * pa_y[i];

        g_yyyyyzz_xzzz_0[i] = 4.0 * g_yyyzz_xzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xzzz_1[i] * fz_be_0 + g_yyyyzz_xzzz_1[i] * pa_y[i];

        g_yyyyyzz_yyyy_0[i] = g_yyyyy_yyyy_0[i] * fbe_0 - g_yyyyy_yyyy_1[i] * fz_be_0 + g_yyyyyz_yyyy_1[i] * pa_z[i];

        g_yyyyyzz_yyyz_0[i] = 4.0 * g_yyyzz_yyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_yyyz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_yyz_1[i] * fe_0 + g_yyyyzz_yyyz_1[i] * pa_y[i];

        g_yyyyyzz_yyzz_0[i] = 4.0 * g_yyyzz_yyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_yyzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_yzz_1[i] * fe_0 + g_yyyyzz_yyzz_1[i] * pa_y[i];

        g_yyyyyzz_yzzz_0[i] = 4.0 * g_yyyzz_yzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_yzzz_1[i] * fz_be_0 + g_yyyyzz_zzz_1[i] * fe_0 + g_yyyyzz_yzzz_1[i] * pa_y[i];

        g_yyyyyzz_zzzz_0[i] = 4.0 * g_yyyzz_zzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_zzzz_1[i] * fz_be_0 + g_yyyyzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 465-480 components of targeted buffer : KG

    auto g_yyyyzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 465);

    auto g_yyyyzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 466);

    auto g_yyyyzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 467);

    auto g_yyyyzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 468);

    auto g_yyyyzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 469);

    auto g_yyyyzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 470);

    auto g_yyyyzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 471);

    auto g_yyyyzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 472);

    auto g_yyyyzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 473);

    auto g_yyyyzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 474);

    auto g_yyyyzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 475);

    auto g_yyyyzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 476);

    auto g_yyyyzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 477);

    auto g_yyyyzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 478);

    auto g_yyyyzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 479);

    #pragma omp simd aligned(g_yyyyz_xxxy_0, g_yyyyz_xxxy_1, g_yyyyz_xxyy_0, g_yyyyz_xxyy_1, g_yyyyz_xyyy_0, g_yyyyz_xyyy_1, g_yyyyz_yyyy_0, g_yyyyz_yyyy_1, g_yyyyzz_xxxy_1, g_yyyyzz_xxyy_1, g_yyyyzz_xyyy_1, g_yyyyzz_yyyy_1, g_yyyyzzz_xxxx_0, g_yyyyzzz_xxxy_0, g_yyyyzzz_xxxz_0, g_yyyyzzz_xxyy_0, g_yyyyzzz_xxyz_0, g_yyyyzzz_xxzz_0, g_yyyyzzz_xyyy_0, g_yyyyzzz_xyyz_0, g_yyyyzzz_xyzz_0, g_yyyyzzz_xzzz_0, g_yyyyzzz_yyyy_0, g_yyyyzzz_yyyz_0, g_yyyyzzz_yyzz_0, g_yyyyzzz_yzzz_0, g_yyyyzzz_zzzz_0, g_yyyzzz_xxxx_1, g_yyyzzz_xxxz_1, g_yyyzzz_xxyz_1, g_yyyzzz_xxz_1, g_yyyzzz_xxzz_1, g_yyyzzz_xyyz_1, g_yyyzzz_xyz_1, g_yyyzzz_xyzz_1, g_yyyzzz_xzz_1, g_yyyzzz_xzzz_1, g_yyyzzz_yyyz_1, g_yyyzzz_yyz_1, g_yyyzzz_yyzz_1, g_yyyzzz_yzz_1, g_yyyzzz_yzzz_1, g_yyyzzz_zzz_1, g_yyyzzz_zzzz_1, g_yyzzz_xxxx_0, g_yyzzz_xxxx_1, g_yyzzz_xxxz_0, g_yyzzz_xxxz_1, g_yyzzz_xxyz_0, g_yyzzz_xxyz_1, g_yyzzz_xxzz_0, g_yyzzz_xxzz_1, g_yyzzz_xyyz_0, g_yyzzz_xyyz_1, g_yyzzz_xyzz_0, g_yyzzz_xyzz_1, g_yyzzz_xzzz_0, g_yyzzz_xzzz_1, g_yyzzz_yyyz_0, g_yyzzz_yyyz_1, g_yyzzz_yyzz_0, g_yyzzz_yyzz_1, g_yyzzz_yzzz_0, g_yyzzz_yzzz_1, g_yyzzz_zzzz_0, g_yyzzz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyzzz_xxxx_0[i] = 3.0 * g_yyzzz_xxxx_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxx_1[i] * fz_be_0 + g_yyyzzz_xxxx_1[i] * pa_y[i];

        g_yyyyzzz_xxxy_0[i] = 2.0 * g_yyyyz_xxxy_0[i] * fbe_0 - 2.0 * g_yyyyz_xxxy_1[i] * fz_be_0 + g_yyyyzz_xxxy_1[i] * pa_z[i];

        g_yyyyzzz_xxxz_0[i] = 3.0 * g_yyzzz_xxxz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxz_1[i] * fz_be_0 + g_yyyzzz_xxxz_1[i] * pa_y[i];

        g_yyyyzzz_xxyy_0[i] = 2.0 * g_yyyyz_xxyy_0[i] * fbe_0 - 2.0 * g_yyyyz_xxyy_1[i] * fz_be_0 + g_yyyyzz_xxyy_1[i] * pa_z[i];

        g_yyyyzzz_xxyz_0[i] = 3.0 * g_yyzzz_xxyz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxyz_1[i] * fz_be_0 + g_yyyzzz_xxz_1[i] * fe_0 + g_yyyzzz_xxyz_1[i] * pa_y[i];

        g_yyyyzzz_xxzz_0[i] = 3.0 * g_yyzzz_xxzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxzz_1[i] * fz_be_0 + g_yyyzzz_xxzz_1[i] * pa_y[i];

        g_yyyyzzz_xyyy_0[i] = 2.0 * g_yyyyz_xyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_xyyy_1[i] * fz_be_0 + g_yyyyzz_xyyy_1[i] * pa_z[i];

        g_yyyyzzz_xyyz_0[i] = 3.0 * g_yyzzz_xyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_xyyz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_xyz_1[i] * fe_0 + g_yyyzzz_xyyz_1[i] * pa_y[i];

        g_yyyyzzz_xyzz_0[i] = 3.0 * g_yyzzz_xyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xyzz_1[i] * fz_be_0 + g_yyyzzz_xzz_1[i] * fe_0 + g_yyyzzz_xyzz_1[i] * pa_y[i];

        g_yyyyzzz_xzzz_0[i] = 3.0 * g_yyzzz_xzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xzzz_1[i] * fz_be_0 + g_yyyzzz_xzzz_1[i] * pa_y[i];

        g_yyyyzzz_yyyy_0[i] = 2.0 * g_yyyyz_yyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_yyyy_1[i] * fz_be_0 + g_yyyyzz_yyyy_1[i] * pa_z[i];

        g_yyyyzzz_yyyz_0[i] = 3.0 * g_yyzzz_yyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_yyyz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_yyz_1[i] * fe_0 + g_yyyzzz_yyyz_1[i] * pa_y[i];

        g_yyyyzzz_yyzz_0[i] = 3.0 * g_yyzzz_yyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_yyzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_yzz_1[i] * fe_0 + g_yyyzzz_yyzz_1[i] * pa_y[i];

        g_yyyyzzz_yzzz_0[i] = 3.0 * g_yyzzz_yzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_yzzz_1[i] * fz_be_0 + g_yyyzzz_zzz_1[i] * fe_0 + g_yyyzzz_yzzz_1[i] * pa_y[i];

        g_yyyyzzz_zzzz_0[i] = 3.0 * g_yyzzz_zzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_zzzz_1[i] * fz_be_0 + g_yyyzzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 480-495 components of targeted buffer : KG

    auto g_yyyzzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 480);

    auto g_yyyzzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 481);

    auto g_yyyzzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 482);

    auto g_yyyzzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 483);

    auto g_yyyzzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 484);

    auto g_yyyzzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 485);

    auto g_yyyzzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 486);

    auto g_yyyzzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 487);

    auto g_yyyzzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 488);

    auto g_yyyzzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 489);

    auto g_yyyzzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 490);

    auto g_yyyzzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 491);

    auto g_yyyzzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 492);

    auto g_yyyzzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 493);

    auto g_yyyzzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 494);

    #pragma omp simd aligned(g_yyyzz_xxxy_0, g_yyyzz_xxxy_1, g_yyyzz_xxyy_0, g_yyyzz_xxyy_1, g_yyyzz_xyyy_0, g_yyyzz_xyyy_1, g_yyyzz_yyyy_0, g_yyyzz_yyyy_1, g_yyyzzz_xxxy_1, g_yyyzzz_xxyy_1, g_yyyzzz_xyyy_1, g_yyyzzz_yyyy_1, g_yyyzzzz_xxxx_0, g_yyyzzzz_xxxy_0, g_yyyzzzz_xxxz_0, g_yyyzzzz_xxyy_0, g_yyyzzzz_xxyz_0, g_yyyzzzz_xxzz_0, g_yyyzzzz_xyyy_0, g_yyyzzzz_xyyz_0, g_yyyzzzz_xyzz_0, g_yyyzzzz_xzzz_0, g_yyyzzzz_yyyy_0, g_yyyzzzz_yyyz_0, g_yyyzzzz_yyzz_0, g_yyyzzzz_yzzz_0, g_yyyzzzz_zzzz_0, g_yyzzzz_xxxx_1, g_yyzzzz_xxxz_1, g_yyzzzz_xxyz_1, g_yyzzzz_xxz_1, g_yyzzzz_xxzz_1, g_yyzzzz_xyyz_1, g_yyzzzz_xyz_1, g_yyzzzz_xyzz_1, g_yyzzzz_xzz_1, g_yyzzzz_xzzz_1, g_yyzzzz_yyyz_1, g_yyzzzz_yyz_1, g_yyzzzz_yyzz_1, g_yyzzzz_yzz_1, g_yyzzzz_yzzz_1, g_yyzzzz_zzz_1, g_yyzzzz_zzzz_1, g_yzzzz_xxxx_0, g_yzzzz_xxxx_1, g_yzzzz_xxxz_0, g_yzzzz_xxxz_1, g_yzzzz_xxyz_0, g_yzzzz_xxyz_1, g_yzzzz_xxzz_0, g_yzzzz_xxzz_1, g_yzzzz_xyyz_0, g_yzzzz_xyyz_1, g_yzzzz_xyzz_0, g_yzzzz_xyzz_1, g_yzzzz_xzzz_0, g_yzzzz_xzzz_1, g_yzzzz_yyyz_0, g_yzzzz_yyyz_1, g_yzzzz_yyzz_0, g_yzzzz_yyzz_1, g_yzzzz_yzzz_0, g_yzzzz_yzzz_1, g_yzzzz_zzzz_0, g_yzzzz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzzzz_xxxx_0[i] = 2.0 * g_yzzzz_xxxx_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxx_1[i] * fz_be_0 + g_yyzzzz_xxxx_1[i] * pa_y[i];

        g_yyyzzzz_xxxy_0[i] = 3.0 * g_yyyzz_xxxy_0[i] * fbe_0 - 3.0 * g_yyyzz_xxxy_1[i] * fz_be_0 + g_yyyzzz_xxxy_1[i] * pa_z[i];

        g_yyyzzzz_xxxz_0[i] = 2.0 * g_yzzzz_xxxz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxz_1[i] * fz_be_0 + g_yyzzzz_xxxz_1[i] * pa_y[i];

        g_yyyzzzz_xxyy_0[i] = 3.0 * g_yyyzz_xxyy_0[i] * fbe_0 - 3.0 * g_yyyzz_xxyy_1[i] * fz_be_0 + g_yyyzzz_xxyy_1[i] * pa_z[i];

        g_yyyzzzz_xxyz_0[i] = 2.0 * g_yzzzz_xxyz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxyz_1[i] * fz_be_0 + g_yyzzzz_xxz_1[i] * fe_0 + g_yyzzzz_xxyz_1[i] * pa_y[i];

        g_yyyzzzz_xxzz_0[i] = 2.0 * g_yzzzz_xxzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxzz_1[i] * fz_be_0 + g_yyzzzz_xxzz_1[i] * pa_y[i];

        g_yyyzzzz_xyyy_0[i] = 3.0 * g_yyyzz_xyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_xyyy_1[i] * fz_be_0 + g_yyyzzz_xyyy_1[i] * pa_z[i];

        g_yyyzzzz_xyyz_0[i] = 2.0 * g_yzzzz_xyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_xyyz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_xyz_1[i] * fe_0 + g_yyzzzz_xyyz_1[i] * pa_y[i];

        g_yyyzzzz_xyzz_0[i] = 2.0 * g_yzzzz_xyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xyzz_1[i] * fz_be_0 + g_yyzzzz_xzz_1[i] * fe_0 + g_yyzzzz_xyzz_1[i] * pa_y[i];

        g_yyyzzzz_xzzz_0[i] = 2.0 * g_yzzzz_xzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xzzz_1[i] * fz_be_0 + g_yyzzzz_xzzz_1[i] * pa_y[i];

        g_yyyzzzz_yyyy_0[i] = 3.0 * g_yyyzz_yyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_yyyy_1[i] * fz_be_0 + g_yyyzzz_yyyy_1[i] * pa_z[i];

        g_yyyzzzz_yyyz_0[i] = 2.0 * g_yzzzz_yyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_yyyz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_yyz_1[i] * fe_0 + g_yyzzzz_yyyz_1[i] * pa_y[i];

        g_yyyzzzz_yyzz_0[i] = 2.0 * g_yzzzz_yyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_yyzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_yzz_1[i] * fe_0 + g_yyzzzz_yyzz_1[i] * pa_y[i];

        g_yyyzzzz_yzzz_0[i] = 2.0 * g_yzzzz_yzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_yzzz_1[i] * fz_be_0 + g_yyzzzz_zzz_1[i] * fe_0 + g_yyzzzz_yzzz_1[i] * pa_y[i];

        g_yyyzzzz_zzzz_0[i] = 2.0 * g_yzzzz_zzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_zzzz_1[i] * fz_be_0 + g_yyzzzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 495-510 components of targeted buffer : KG

    auto g_yyzzzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 495);

    auto g_yyzzzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 496);

    auto g_yyzzzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 497);

    auto g_yyzzzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 498);

    auto g_yyzzzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 499);

    auto g_yyzzzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 500);

    auto g_yyzzzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 501);

    auto g_yyzzzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 502);

    auto g_yyzzzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 503);

    auto g_yyzzzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 504);

    auto g_yyzzzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 505);

    auto g_yyzzzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 506);

    auto g_yyzzzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 507);

    auto g_yyzzzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 508);

    auto g_yyzzzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 509);

    #pragma omp simd aligned(g_yyzzz_xxxy_0, g_yyzzz_xxxy_1, g_yyzzz_xxyy_0, g_yyzzz_xxyy_1, g_yyzzz_xyyy_0, g_yyzzz_xyyy_1, g_yyzzz_yyyy_0, g_yyzzz_yyyy_1, g_yyzzzz_xxxy_1, g_yyzzzz_xxyy_1, g_yyzzzz_xyyy_1, g_yyzzzz_yyyy_1, g_yyzzzzz_xxxx_0, g_yyzzzzz_xxxy_0, g_yyzzzzz_xxxz_0, g_yyzzzzz_xxyy_0, g_yyzzzzz_xxyz_0, g_yyzzzzz_xxzz_0, g_yyzzzzz_xyyy_0, g_yyzzzzz_xyyz_0, g_yyzzzzz_xyzz_0, g_yyzzzzz_xzzz_0, g_yyzzzzz_yyyy_0, g_yyzzzzz_yyyz_0, g_yyzzzzz_yyzz_0, g_yyzzzzz_yzzz_0, g_yyzzzzz_zzzz_0, g_yzzzzz_xxxx_1, g_yzzzzz_xxxz_1, g_yzzzzz_xxyz_1, g_yzzzzz_xxz_1, g_yzzzzz_xxzz_1, g_yzzzzz_xyyz_1, g_yzzzzz_xyz_1, g_yzzzzz_xyzz_1, g_yzzzzz_xzz_1, g_yzzzzz_xzzz_1, g_yzzzzz_yyyz_1, g_yzzzzz_yyz_1, g_yzzzzz_yyzz_1, g_yzzzzz_yzz_1, g_yzzzzz_yzzz_1, g_yzzzzz_zzz_1, g_yzzzzz_zzzz_1, g_zzzzz_xxxx_0, g_zzzzz_xxxx_1, g_zzzzz_xxxz_0, g_zzzzz_xxxz_1, g_zzzzz_xxyz_0, g_zzzzz_xxyz_1, g_zzzzz_xxzz_0, g_zzzzz_xxzz_1, g_zzzzz_xyyz_0, g_zzzzz_xyyz_1, g_zzzzz_xyzz_0, g_zzzzz_xyzz_1, g_zzzzz_xzzz_0, g_zzzzz_xzzz_1, g_zzzzz_yyyz_0, g_zzzzz_yyyz_1, g_zzzzz_yyzz_0, g_zzzzz_yyzz_1, g_zzzzz_yzzz_0, g_zzzzz_yzzz_1, g_zzzzz_zzzz_0, g_zzzzz_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzzzz_xxxx_0[i] = g_zzzzz_xxxx_0[i] * fbe_0 - g_zzzzz_xxxx_1[i] * fz_be_0 + g_yzzzzz_xxxx_1[i] * pa_y[i];

        g_yyzzzzz_xxxy_0[i] = 4.0 * g_yyzzz_xxxy_0[i] * fbe_0 - 4.0 * g_yyzzz_xxxy_1[i] * fz_be_0 + g_yyzzzz_xxxy_1[i] * pa_z[i];

        g_yyzzzzz_xxxz_0[i] = g_zzzzz_xxxz_0[i] * fbe_0 - g_zzzzz_xxxz_1[i] * fz_be_0 + g_yzzzzz_xxxz_1[i] * pa_y[i];

        g_yyzzzzz_xxyy_0[i] = 4.0 * g_yyzzz_xxyy_0[i] * fbe_0 - 4.0 * g_yyzzz_xxyy_1[i] * fz_be_0 + g_yyzzzz_xxyy_1[i] * pa_z[i];

        g_yyzzzzz_xxyz_0[i] = g_zzzzz_xxyz_0[i] * fbe_0 - g_zzzzz_xxyz_1[i] * fz_be_0 + g_yzzzzz_xxz_1[i] * fe_0 + g_yzzzzz_xxyz_1[i] * pa_y[i];

        g_yyzzzzz_xxzz_0[i] = g_zzzzz_xxzz_0[i] * fbe_0 - g_zzzzz_xxzz_1[i] * fz_be_0 + g_yzzzzz_xxzz_1[i] * pa_y[i];

        g_yyzzzzz_xyyy_0[i] = 4.0 * g_yyzzz_xyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_xyyy_1[i] * fz_be_0 + g_yyzzzz_xyyy_1[i] * pa_z[i];

        g_yyzzzzz_xyyz_0[i] = g_zzzzz_xyyz_0[i] * fbe_0 - g_zzzzz_xyyz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_xyz_1[i] * fe_0 + g_yzzzzz_xyyz_1[i] * pa_y[i];

        g_yyzzzzz_xyzz_0[i] = g_zzzzz_xyzz_0[i] * fbe_0 - g_zzzzz_xyzz_1[i] * fz_be_0 + g_yzzzzz_xzz_1[i] * fe_0 + g_yzzzzz_xyzz_1[i] * pa_y[i];

        g_yyzzzzz_xzzz_0[i] = g_zzzzz_xzzz_0[i] * fbe_0 - g_zzzzz_xzzz_1[i] * fz_be_0 + g_yzzzzz_xzzz_1[i] * pa_y[i];

        g_yyzzzzz_yyyy_0[i] = 4.0 * g_yyzzz_yyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_yyyy_1[i] * fz_be_0 + g_yyzzzz_yyyy_1[i] * pa_z[i];

        g_yyzzzzz_yyyz_0[i] = g_zzzzz_yyyz_0[i] * fbe_0 - g_zzzzz_yyyz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_yyz_1[i] * fe_0 + g_yzzzzz_yyyz_1[i] * pa_y[i];

        g_yyzzzzz_yyzz_0[i] = g_zzzzz_yyzz_0[i] * fbe_0 - g_zzzzz_yyzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_yzz_1[i] * fe_0 + g_yzzzzz_yyzz_1[i] * pa_y[i];

        g_yyzzzzz_yzzz_0[i] = g_zzzzz_yzzz_0[i] * fbe_0 - g_zzzzz_yzzz_1[i] * fz_be_0 + g_yzzzzz_zzz_1[i] * fe_0 + g_yzzzzz_yzzz_1[i] * pa_y[i];

        g_yyzzzzz_zzzz_0[i] = g_zzzzz_zzzz_0[i] * fbe_0 - g_zzzzz_zzzz_1[i] * fz_be_0 + g_yzzzzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 510-525 components of targeted buffer : KG

    auto g_yzzzzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 510);

    auto g_yzzzzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 511);

    auto g_yzzzzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 512);

    auto g_yzzzzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 513);

    auto g_yzzzzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 514);

    auto g_yzzzzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 515);

    auto g_yzzzzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 516);

    auto g_yzzzzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 517);

    auto g_yzzzzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 518);

    auto g_yzzzzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 519);

    auto g_yzzzzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 520);

    auto g_yzzzzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 521);

    auto g_yzzzzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 522);

    auto g_yzzzzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 523);

    auto g_yzzzzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 524);

    #pragma omp simd aligned(g_yzzzzzz_xxxx_0, g_yzzzzzz_xxxy_0, g_yzzzzzz_xxxz_0, g_yzzzzzz_xxyy_0, g_yzzzzzz_xxyz_0, g_yzzzzzz_xxzz_0, g_yzzzzzz_xyyy_0, g_yzzzzzz_xyyz_0, g_yzzzzzz_xyzz_0, g_yzzzzzz_xzzz_0, g_yzzzzzz_yyyy_0, g_yzzzzzz_yyyz_0, g_yzzzzzz_yyzz_0, g_yzzzzzz_yzzz_0, g_yzzzzzz_zzzz_0, g_zzzzzz_xxx_1, g_zzzzzz_xxxx_1, g_zzzzzz_xxxy_1, g_zzzzzz_xxxz_1, g_zzzzzz_xxy_1, g_zzzzzz_xxyy_1, g_zzzzzz_xxyz_1, g_zzzzzz_xxz_1, g_zzzzzz_xxzz_1, g_zzzzzz_xyy_1, g_zzzzzz_xyyy_1, g_zzzzzz_xyyz_1, g_zzzzzz_xyz_1, g_zzzzzz_xyzz_1, g_zzzzzz_xzz_1, g_zzzzzz_xzzz_1, g_zzzzzz_yyy_1, g_zzzzzz_yyyy_1, g_zzzzzz_yyyz_1, g_zzzzzz_yyz_1, g_zzzzzz_yyzz_1, g_zzzzzz_yzz_1, g_zzzzzz_yzzz_1, g_zzzzzz_zzz_1, g_zzzzzz_zzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzzzz_xxxx_0[i] = g_zzzzzz_xxxx_1[i] * pa_y[i];

        g_yzzzzzz_xxxy_0[i] = g_zzzzzz_xxx_1[i] * fe_0 + g_zzzzzz_xxxy_1[i] * pa_y[i];

        g_yzzzzzz_xxxz_0[i] = g_zzzzzz_xxxz_1[i] * pa_y[i];

        g_yzzzzzz_xxyy_0[i] = 2.0 * g_zzzzzz_xxy_1[i] * fe_0 + g_zzzzzz_xxyy_1[i] * pa_y[i];

        g_yzzzzzz_xxyz_0[i] = g_zzzzzz_xxz_1[i] * fe_0 + g_zzzzzz_xxyz_1[i] * pa_y[i];

        g_yzzzzzz_xxzz_0[i] = g_zzzzzz_xxzz_1[i] * pa_y[i];

        g_yzzzzzz_xyyy_0[i] = 3.0 * g_zzzzzz_xyy_1[i] * fe_0 + g_zzzzzz_xyyy_1[i] * pa_y[i];

        g_yzzzzzz_xyyz_0[i] = 2.0 * g_zzzzzz_xyz_1[i] * fe_0 + g_zzzzzz_xyyz_1[i] * pa_y[i];

        g_yzzzzzz_xyzz_0[i] = g_zzzzzz_xzz_1[i] * fe_0 + g_zzzzzz_xyzz_1[i] * pa_y[i];

        g_yzzzzzz_xzzz_0[i] = g_zzzzzz_xzzz_1[i] * pa_y[i];

        g_yzzzzzz_yyyy_0[i] = 4.0 * g_zzzzzz_yyy_1[i] * fe_0 + g_zzzzzz_yyyy_1[i] * pa_y[i];

        g_yzzzzzz_yyyz_0[i] = 3.0 * g_zzzzzz_yyz_1[i] * fe_0 + g_zzzzzz_yyyz_1[i] * pa_y[i];

        g_yzzzzzz_yyzz_0[i] = 2.0 * g_zzzzzz_yzz_1[i] * fe_0 + g_zzzzzz_yyzz_1[i] * pa_y[i];

        g_yzzzzzz_yzzz_0[i] = g_zzzzzz_zzz_1[i] * fe_0 + g_zzzzzz_yzzz_1[i] * pa_y[i];

        g_yzzzzzz_zzzz_0[i] = g_zzzzzz_zzzz_1[i] * pa_y[i];
    }

    // Set up 525-540 components of targeted buffer : KG

    auto g_zzzzzzz_xxxx_0 = pbuffer.data(idx_eri_0_kg + 525);

    auto g_zzzzzzz_xxxy_0 = pbuffer.data(idx_eri_0_kg + 526);

    auto g_zzzzzzz_xxxz_0 = pbuffer.data(idx_eri_0_kg + 527);

    auto g_zzzzzzz_xxyy_0 = pbuffer.data(idx_eri_0_kg + 528);

    auto g_zzzzzzz_xxyz_0 = pbuffer.data(idx_eri_0_kg + 529);

    auto g_zzzzzzz_xxzz_0 = pbuffer.data(idx_eri_0_kg + 530);

    auto g_zzzzzzz_xyyy_0 = pbuffer.data(idx_eri_0_kg + 531);

    auto g_zzzzzzz_xyyz_0 = pbuffer.data(idx_eri_0_kg + 532);

    auto g_zzzzzzz_xyzz_0 = pbuffer.data(idx_eri_0_kg + 533);

    auto g_zzzzzzz_xzzz_0 = pbuffer.data(idx_eri_0_kg + 534);

    auto g_zzzzzzz_yyyy_0 = pbuffer.data(idx_eri_0_kg + 535);

    auto g_zzzzzzz_yyyz_0 = pbuffer.data(idx_eri_0_kg + 536);

    auto g_zzzzzzz_yyzz_0 = pbuffer.data(idx_eri_0_kg + 537);

    auto g_zzzzzzz_yzzz_0 = pbuffer.data(idx_eri_0_kg + 538);

    auto g_zzzzzzz_zzzz_0 = pbuffer.data(idx_eri_0_kg + 539);

    #pragma omp simd aligned(g_zzzzz_xxxx_0, g_zzzzz_xxxx_1, g_zzzzz_xxxy_0, g_zzzzz_xxxy_1, g_zzzzz_xxxz_0, g_zzzzz_xxxz_1, g_zzzzz_xxyy_0, g_zzzzz_xxyy_1, g_zzzzz_xxyz_0, g_zzzzz_xxyz_1, g_zzzzz_xxzz_0, g_zzzzz_xxzz_1, g_zzzzz_xyyy_0, g_zzzzz_xyyy_1, g_zzzzz_xyyz_0, g_zzzzz_xyyz_1, g_zzzzz_xyzz_0, g_zzzzz_xyzz_1, g_zzzzz_xzzz_0, g_zzzzz_xzzz_1, g_zzzzz_yyyy_0, g_zzzzz_yyyy_1, g_zzzzz_yyyz_0, g_zzzzz_yyyz_1, g_zzzzz_yyzz_0, g_zzzzz_yyzz_1, g_zzzzz_yzzz_0, g_zzzzz_yzzz_1, g_zzzzz_zzzz_0, g_zzzzz_zzzz_1, g_zzzzzz_xxx_1, g_zzzzzz_xxxx_1, g_zzzzzz_xxxy_1, g_zzzzzz_xxxz_1, g_zzzzzz_xxy_1, g_zzzzzz_xxyy_1, g_zzzzzz_xxyz_1, g_zzzzzz_xxz_1, g_zzzzzz_xxzz_1, g_zzzzzz_xyy_1, g_zzzzzz_xyyy_1, g_zzzzzz_xyyz_1, g_zzzzzz_xyz_1, g_zzzzzz_xyzz_1, g_zzzzzz_xzz_1, g_zzzzzz_xzzz_1, g_zzzzzz_yyy_1, g_zzzzzz_yyyy_1, g_zzzzzz_yyyz_1, g_zzzzzz_yyz_1, g_zzzzzz_yyzz_1, g_zzzzzz_yzz_1, g_zzzzzz_yzzz_1, g_zzzzzz_zzz_1, g_zzzzzz_zzzz_1, g_zzzzzzz_xxxx_0, g_zzzzzzz_xxxy_0, g_zzzzzzz_xxxz_0, g_zzzzzzz_xxyy_0, g_zzzzzzz_xxyz_0, g_zzzzzzz_xxzz_0, g_zzzzzzz_xyyy_0, g_zzzzzzz_xyyz_0, g_zzzzzzz_xyzz_0, g_zzzzzzz_xzzz_0, g_zzzzzzz_yyyy_0, g_zzzzzzz_yyyz_0, g_zzzzzzz_yyzz_0, g_zzzzzzz_yzzz_0, g_zzzzzzz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzzzz_xxxx_0[i] = 6.0 * g_zzzzz_xxxx_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxx_1[i] * fz_be_0 + g_zzzzzz_xxxx_1[i] * pa_z[i];

        g_zzzzzzz_xxxy_0[i] = 6.0 * g_zzzzz_xxxy_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxy_1[i] * fz_be_0 + g_zzzzzz_xxxy_1[i] * pa_z[i];

        g_zzzzzzz_xxxz_0[i] = 6.0 * g_zzzzz_xxxz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxz_1[i] * fz_be_0 + g_zzzzzz_xxx_1[i] * fe_0 + g_zzzzzz_xxxz_1[i] * pa_z[i];

        g_zzzzzzz_xxyy_0[i] = 6.0 * g_zzzzz_xxyy_0[i] * fbe_0 - 6.0 * g_zzzzz_xxyy_1[i] * fz_be_0 + g_zzzzzz_xxyy_1[i] * pa_z[i];

        g_zzzzzzz_xxyz_0[i] = 6.0 * g_zzzzz_xxyz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxyz_1[i] * fz_be_0 + g_zzzzzz_xxy_1[i] * fe_0 + g_zzzzzz_xxyz_1[i] * pa_z[i];

        g_zzzzzzz_xxzz_0[i] = 6.0 * g_zzzzz_xxzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_xxz_1[i] * fe_0 + g_zzzzzz_xxzz_1[i] * pa_z[i];

        g_zzzzzzz_xyyy_0[i] = 6.0 * g_zzzzz_xyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_xyyy_1[i] * fz_be_0 + g_zzzzzz_xyyy_1[i] * pa_z[i];

        g_zzzzzzz_xyyz_0[i] = 6.0 * g_zzzzz_xyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_xyyz_1[i] * fz_be_0 + g_zzzzzz_xyy_1[i] * fe_0 + g_zzzzzz_xyyz_1[i] * pa_z[i];

        g_zzzzzzz_xyzz_0[i] = 6.0 * g_zzzzz_xyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_xyz_1[i] * fe_0 + g_zzzzzz_xyzz_1[i] * pa_z[i];

        g_zzzzzzz_xzzz_0[i] = 6.0 * g_zzzzz_xzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_xzz_1[i] * fe_0 + g_zzzzzz_xzzz_1[i] * pa_z[i];

        g_zzzzzzz_yyyy_0[i] = 6.0 * g_zzzzz_yyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_yyyy_1[i] * fz_be_0 + g_zzzzzz_yyyy_1[i] * pa_z[i];

        g_zzzzzzz_yyyz_0[i] = 6.0 * g_zzzzz_yyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_yyyz_1[i] * fz_be_0 + g_zzzzzz_yyy_1[i] * fe_0 + g_zzzzzz_yyyz_1[i] * pa_z[i];

        g_zzzzzzz_yyzz_0[i] = 6.0 * g_zzzzz_yyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_yyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_yyz_1[i] * fe_0 + g_zzzzzz_yyzz_1[i] * pa_z[i];

        g_zzzzzzz_yzzz_0[i] = 6.0 * g_zzzzz_yzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_yzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_yzz_1[i] * fe_0 + g_zzzzzz_yzzz_1[i] * pa_z[i];

        g_zzzzzzz_zzzz_0[i] = 6.0 * g_zzzzz_zzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_zzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_zzz_1[i] * fe_0 + g_zzzzzz_zzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

