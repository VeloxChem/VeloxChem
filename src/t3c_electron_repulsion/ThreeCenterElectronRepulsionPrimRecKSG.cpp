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

#include "ThreeCenterElectronRepulsionPrimRecKSG.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ksg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksg,
                                 size_t idx_eri_0_hsg,
                                 size_t idx_eri_1_hsg,
                                 size_t idx_eri_1_isf,
                                 size_t idx_eri_1_isg,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : HSG

    auto g_xxxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg);

    auto g_xxxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 1);

    auto g_xxxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 2);

    auto g_xxxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 3);

    auto g_xxxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 4);

    auto g_xxxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 5);

    auto g_xxxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 6);

    auto g_xxxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 7);

    auto g_xxxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 8);

    auto g_xxxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 9);

    auto g_xxxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 10);

    auto g_xxxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 11);

    auto g_xxxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 12);

    auto g_xxxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 13);

    auto g_xxxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 14);

    auto g_xxxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 15);

    auto g_xxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 17);

    auto g_xxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 20);

    auto g_xxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 24);

    auto g_xxxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 30);

    auto g_xxxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 31);

    auto g_xxxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 33);

    auto g_xxxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 36);

    auto g_xxxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 45);

    auto g_xxxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 46);

    auto g_xxxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 47);

    auto g_xxxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 48);

    auto g_xxxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 49);

    auto g_xxxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 50);

    auto g_xxxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 51);

    auto g_xxxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 52);

    auto g_xxxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 53);

    auto g_xxxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 54);

    auto g_xxxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 55);

    auto g_xxxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 56);

    auto g_xxxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 57);

    auto g_xxxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 58);

    auto g_xxxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 59);

    auto g_xxxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 75);

    auto g_xxxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 76);

    auto g_xxxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 77);

    auto g_xxxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 78);

    auto g_xxxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 79);

    auto g_xxxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 80);

    auto g_xxxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 81);

    auto g_xxxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 82);

    auto g_xxxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 83);

    auto g_xxxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 84);

    auto g_xxxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 85);

    auto g_xxxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 86);

    auto g_xxxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 87);

    auto g_xxxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 88);

    auto g_xxxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 89);

    auto g_xxyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 90);

    auto g_xxyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 91);

    auto g_xxyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 92);

    auto g_xxyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 93);

    auto g_xxyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 94);

    auto g_xxyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 95);

    auto g_xxyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 96);

    auto g_xxyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 97);

    auto g_xxyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 98);

    auto g_xxyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 99);

    auto g_xxyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 100);

    auto g_xxyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 101);

    auto g_xxyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 102);

    auto g_xxyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 103);

    auto g_xxyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 104);

    auto g_xxyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 106);

    auto g_xxyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 108);

    auto g_xxyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 111);

    auto g_xxyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 120);

    auto g_xxyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 122);

    auto g_xxyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 125);

    auto g_xxyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 129);

    auto g_xxzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 135);

    auto g_xxzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 136);

    auto g_xxzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 137);

    auto g_xxzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 138);

    auto g_xxzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 139);

    auto g_xxzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 140);

    auto g_xxzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 141);

    auto g_xxzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 142);

    auto g_xxzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 143);

    auto g_xxzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 144);

    auto g_xxzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 145);

    auto g_xxzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 146);

    auto g_xxzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 147);

    auto g_xxzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 148);

    auto g_xxzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 149);

    auto g_xyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 151);

    auto g_xyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 153);

    auto g_xyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 154);

    auto g_xyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 156);

    auto g_xyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 157);

    auto g_xyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 158);

    auto g_xyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 160);

    auto g_xyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 161);

    auto g_xyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 162);

    auto g_xyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 163);

    auto g_xyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 164);

    auto g_xyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 184);

    auto g_xyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 187);

    auto g_xyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 188);

    auto g_xyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 190);

    auto g_xyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 191);

    auto g_xyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 192);

    auto g_xyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 193);

    auto g_xyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 194);

    auto g_xzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 212);

    auto g_xzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 214);

    auto g_xzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 215);

    auto g_xzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 217);

    auto g_xzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 218);

    auto g_xzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 219);

    auto g_xzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 220);

    auto g_xzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 221);

    auto g_xzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 222);

    auto g_xzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 223);

    auto g_xzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 224);

    auto g_yyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 225);

    auto g_yyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 226);

    auto g_yyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 227);

    auto g_yyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 228);

    auto g_yyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 229);

    auto g_yyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 230);

    auto g_yyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 231);

    auto g_yyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 232);

    auto g_yyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 233);

    auto g_yyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 234);

    auto g_yyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 235);

    auto g_yyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 236);

    auto g_yyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 237);

    auto g_yyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 238);

    auto g_yyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 239);

    auto g_yyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 241);

    auto g_yyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 243);

    auto g_yyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 246);

    auto g_yyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 250);

    auto g_yyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 255);

    auto g_yyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 256);

    auto g_yyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 257);

    auto g_yyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 258);

    auto g_yyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 259);

    auto g_yyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 260);

    auto g_yyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 261);

    auto g_yyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 262);

    auto g_yyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 263);

    auto g_yyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 264);

    auto g_yyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 265);

    auto g_yyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 266);

    auto g_yyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 267);

    auto g_yyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 268);

    auto g_yyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 269);

    auto g_yyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 270);

    auto g_yyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 271);

    auto g_yyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 272);

    auto g_yyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 273);

    auto g_yyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 274);

    auto g_yyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 275);

    auto g_yyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 276);

    auto g_yyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 277);

    auto g_yyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 278);

    auto g_yyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 279);

    auto g_yyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 280);

    auto g_yyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 281);

    auto g_yyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 282);

    auto g_yyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 283);

    auto g_yyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 284);

    auto g_yzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 285);

    auto g_yzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 287);

    auto g_yzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 289);

    auto g_yzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 290);

    auto g_yzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 292);

    auto g_yzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 293);

    auto g_yzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 294);

    auto g_yzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 296);

    auto g_yzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 297);

    auto g_yzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 298);

    auto g_yzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 299);

    auto g_zzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_hsg + 300);

    auto g_zzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_hsg + 301);

    auto g_zzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_hsg + 302);

    auto g_zzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_hsg + 303);

    auto g_zzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_hsg + 304);

    auto g_zzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_hsg + 305);

    auto g_zzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_hsg + 306);

    auto g_zzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_hsg + 307);

    auto g_zzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_hsg + 308);

    auto g_zzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_hsg + 309);

    auto g_zzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_hsg + 310);

    auto g_zzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_hsg + 311);

    auto g_zzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_hsg + 312);

    auto g_zzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_hsg + 313);

    auto g_zzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_hsg + 314);

    /// Set up components of auxilary buffer : HSG

    auto g_xxxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg);

    auto g_xxxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 1);

    auto g_xxxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 2);

    auto g_xxxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 3);

    auto g_xxxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 4);

    auto g_xxxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 5);

    auto g_xxxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 6);

    auto g_xxxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 7);

    auto g_xxxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 8);

    auto g_xxxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 9);

    auto g_xxxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 10);

    auto g_xxxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 11);

    auto g_xxxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 12);

    auto g_xxxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 13);

    auto g_xxxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 14);

    auto g_xxxxy_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 15);

    auto g_xxxxy_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 17);

    auto g_xxxxy_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 20);

    auto g_xxxxy_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 24);

    auto g_xxxxz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 30);

    auto g_xxxxz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 31);

    auto g_xxxxz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 33);

    auto g_xxxxz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 36);

    auto g_xxxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 45);

    auto g_xxxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 46);

    auto g_xxxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 47);

    auto g_xxxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 48);

    auto g_xxxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 49);

    auto g_xxxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 50);

    auto g_xxxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 51);

    auto g_xxxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 52);

    auto g_xxxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 53);

    auto g_xxxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 54);

    auto g_xxxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 55);

    auto g_xxxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 56);

    auto g_xxxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 57);

    auto g_xxxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 58);

    auto g_xxxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 59);

    auto g_xxxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 75);

    auto g_xxxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 76);

    auto g_xxxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 77);

    auto g_xxxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 78);

    auto g_xxxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 79);

    auto g_xxxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 80);

    auto g_xxxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 81);

    auto g_xxxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 82);

    auto g_xxxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 83);

    auto g_xxxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 84);

    auto g_xxxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 85);

    auto g_xxxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 86);

    auto g_xxxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 87);

    auto g_xxxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 88);

    auto g_xxxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 89);

    auto g_xxyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 90);

    auto g_xxyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 91);

    auto g_xxyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 92);

    auto g_xxyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 93);

    auto g_xxyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 94);

    auto g_xxyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 95);

    auto g_xxyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 96);

    auto g_xxyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 97);

    auto g_xxyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 98);

    auto g_xxyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 99);

    auto g_xxyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 100);

    auto g_xxyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 101);

    auto g_xxyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 102);

    auto g_xxyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 103);

    auto g_xxyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 104);

    auto g_xxyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 106);

    auto g_xxyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 108);

    auto g_xxyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 111);

    auto g_xxyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 120);

    auto g_xxyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 122);

    auto g_xxyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 125);

    auto g_xxyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 129);

    auto g_xxzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 135);

    auto g_xxzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 136);

    auto g_xxzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 137);

    auto g_xxzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 138);

    auto g_xxzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 139);

    auto g_xxzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 140);

    auto g_xxzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 141);

    auto g_xxzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 142);

    auto g_xxzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 143);

    auto g_xxzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 144);

    auto g_xxzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 145);

    auto g_xxzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 146);

    auto g_xxzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 147);

    auto g_xxzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 148);

    auto g_xxzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 149);

    auto g_xyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 151);

    auto g_xyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 153);

    auto g_xyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 154);

    auto g_xyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 156);

    auto g_xyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 157);

    auto g_xyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 158);

    auto g_xyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 160);

    auto g_xyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 161);

    auto g_xyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 162);

    auto g_xyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 163);

    auto g_xyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 164);

    auto g_xyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 184);

    auto g_xyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 187);

    auto g_xyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 188);

    auto g_xyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 190);

    auto g_xyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 191);

    auto g_xyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 192);

    auto g_xyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 193);

    auto g_xyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 194);

    auto g_xzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 212);

    auto g_xzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 214);

    auto g_xzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 215);

    auto g_xzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 217);

    auto g_xzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 218);

    auto g_xzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 219);

    auto g_xzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 220);

    auto g_xzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 221);

    auto g_xzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 222);

    auto g_xzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 223);

    auto g_xzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 224);

    auto g_yyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 225);

    auto g_yyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 226);

    auto g_yyyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 227);

    auto g_yyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 228);

    auto g_yyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 229);

    auto g_yyyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 230);

    auto g_yyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 231);

    auto g_yyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 232);

    auto g_yyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 233);

    auto g_yyyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 234);

    auto g_yyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 235);

    auto g_yyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 236);

    auto g_yyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 237);

    auto g_yyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 238);

    auto g_yyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 239);

    auto g_yyyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 241);

    auto g_yyyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 243);

    auto g_yyyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 246);

    auto g_yyyyz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 250);

    auto g_yyyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 255);

    auto g_yyyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 256);

    auto g_yyyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 257);

    auto g_yyyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 258);

    auto g_yyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 259);

    auto g_yyyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 260);

    auto g_yyyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 261);

    auto g_yyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 262);

    auto g_yyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 263);

    auto g_yyyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 264);

    auto g_yyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 265);

    auto g_yyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 266);

    auto g_yyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 267);

    auto g_yyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 268);

    auto g_yyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 269);

    auto g_yyzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 270);

    auto g_yyzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 271);

    auto g_yyzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 272);

    auto g_yyzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 273);

    auto g_yyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 274);

    auto g_yyzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 275);

    auto g_yyzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 276);

    auto g_yyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 277);

    auto g_yyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 278);

    auto g_yyzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 279);

    auto g_yyzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 280);

    auto g_yyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 281);

    auto g_yyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 282);

    auto g_yyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 283);

    auto g_yyzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 284);

    auto g_yzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 285);

    auto g_yzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 287);

    auto g_yzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 289);

    auto g_yzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 290);

    auto g_yzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 292);

    auto g_yzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 293);

    auto g_yzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 294);

    auto g_yzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 296);

    auto g_yzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 297);

    auto g_yzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 298);

    auto g_yzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 299);

    auto g_zzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_hsg + 300);

    auto g_zzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_hsg + 301);

    auto g_zzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_hsg + 302);

    auto g_zzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_hsg + 303);

    auto g_zzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_hsg + 304);

    auto g_zzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_hsg + 305);

    auto g_zzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_hsg + 306);

    auto g_zzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_hsg + 307);

    auto g_zzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_hsg + 308);

    auto g_zzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_hsg + 309);

    auto g_zzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_hsg + 310);

    auto g_zzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_hsg + 311);

    auto g_zzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_hsg + 312);

    auto g_zzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_hsg + 313);

    auto g_zzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_hsg + 314);

    /// Set up components of auxilary buffer : ISF

    auto g_xxxxxx_0_xxx_1 = pbuffer.data(idx_eri_1_isf);

    auto g_xxxxxx_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 1);

    auto g_xxxxxx_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 2);

    auto g_xxxxxx_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 3);

    auto g_xxxxxx_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 4);

    auto g_xxxxxx_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 5);

    auto g_xxxxxx_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 6);

    auto g_xxxxxx_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 7);

    auto g_xxxxxx_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 8);

    auto g_xxxxxx_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 9);

    auto g_xxxxxz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 22);

    auto g_xxxxxz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 24);

    auto g_xxxxxz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 25);

    auto g_xxxxxz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 27);

    auto g_xxxxxz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 28);

    auto g_xxxxxz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 29);

    auto g_xxxxyy_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 30);

    auto g_xxxxyy_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 31);

    auto g_xxxxyy_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 32);

    auto g_xxxxyy_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 33);

    auto g_xxxxyy_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 34);

    auto g_xxxxyy_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 35);

    auto g_xxxxyy_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 36);

    auto g_xxxxyy_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 37);

    auto g_xxxxyy_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 38);

    auto g_xxxxyy_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 39);

    auto g_xxxxzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 50);

    auto g_xxxxzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 51);

    auto g_xxxxzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 52);

    auto g_xxxxzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 53);

    auto g_xxxxzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 54);

    auto g_xxxxzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 55);

    auto g_xxxxzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 56);

    auto g_xxxxzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 57);

    auto g_xxxxzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 58);

    auto g_xxxxzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 59);

    auto g_xxxyyy_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 60);

    auto g_xxxyyy_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 61);

    auto g_xxxyyy_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 62);

    auto g_xxxyyy_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 63);

    auto g_xxxyyy_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 64);

    auto g_xxxyyy_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 65);

    auto g_xxxyyy_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 66);

    auto g_xxxyyy_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 67);

    auto g_xxxyyy_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 68);

    auto g_xxxyyy_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 69);

    auto g_xxxzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 90);

    auto g_xxxzzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 91);

    auto g_xxxzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 92);

    auto g_xxxzzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 93);

    auto g_xxxzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 94);

    auto g_xxxzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 95);

    auto g_xxxzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 96);

    auto g_xxxzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 97);

    auto g_xxxzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 98);

    auto g_xxxzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 99);

    auto g_xxyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 100);

    auto g_xxyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 101);

    auto g_xxyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 102);

    auto g_xxyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 103);

    auto g_xxyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 104);

    auto g_xxyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 105);

    auto g_xxyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 106);

    auto g_xxyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 107);

    auto g_xxyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 108);

    auto g_xxyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 109);

    auto g_xxyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 124);

    auto g_xxyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 127);

    auto g_xxyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 128);

    auto g_xxzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 140);

    auto g_xxzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 141);

    auto g_xxzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 142);

    auto g_xxzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 143);

    auto g_xxzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 144);

    auto g_xxzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 145);

    auto g_xxzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 146);

    auto g_xxzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 147);

    auto g_xxzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 148);

    auto g_xxzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 149);

    auto g_xyyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 151);

    auto g_xyyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 153);

    auto g_xyyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 154);

    auto g_xyyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 156);

    auto g_xyyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 157);

    auto g_xyyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 158);

    auto g_xyyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 174);

    auto g_xyyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 177);

    auto g_xyyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 178);

    auto g_xyyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 184);

    auto g_xyyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 187);

    auto g_xyyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 188);

    auto g_xzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 202);

    auto g_xzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 204);

    auto g_xzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 205);

    auto g_xzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 207);

    auto g_xzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 208);

    auto g_xzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 209);

    auto g_yyyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 210);

    auto g_yyyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 211);

    auto g_yyyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 212);

    auto g_yyyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 213);

    auto g_yyyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 214);

    auto g_yyyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 215);

    auto g_yyyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 216);

    auto g_yyyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 217);

    auto g_yyyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 218);

    auto g_yyyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 219);

    auto g_yyyyyz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 222);

    auto g_yyyyyz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 224);

    auto g_yyyyyz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 225);

    auto g_yyyyyz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 227);

    auto g_yyyyyz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 228);

    auto g_yyyyyz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 229);

    auto g_yyyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 230);

    auto g_yyyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 231);

    auto g_yyyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 232);

    auto g_yyyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 233);

    auto g_yyyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 234);

    auto g_yyyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 235);

    auto g_yyyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 236);

    auto g_yyyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 237);

    auto g_yyyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 238);

    auto g_yyyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 239);

    auto g_yyyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 240);

    auto g_yyyzzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 241);

    auto g_yyyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 242);

    auto g_yyyzzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 243);

    auto g_yyyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 244);

    auto g_yyyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 245);

    auto g_yyyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 246);

    auto g_yyyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 247);

    auto g_yyyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 248);

    auto g_yyyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 249);

    auto g_yyzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 250);

    auto g_yyzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 251);

    auto g_yyzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 252);

    auto g_yyzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 253);

    auto g_yyzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 254);

    auto g_yyzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 255);

    auto g_yyzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 256);

    auto g_yyzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 257);

    auto g_yyzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 258);

    auto g_yyzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 259);

    auto g_yzzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 261);

    auto g_yzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 262);

    auto g_yzzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 263);

    auto g_yzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 264);

    auto g_yzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 265);

    auto g_yzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 266);

    auto g_yzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 267);

    auto g_yzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 268);

    auto g_yzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 269);

    auto g_zzzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_isf + 270);

    auto g_zzzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_isf + 271);

    auto g_zzzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_isf + 272);

    auto g_zzzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_isf + 273);

    auto g_zzzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_isf + 274);

    auto g_zzzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_isf + 275);

    auto g_zzzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_isf + 276);

    auto g_zzzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_isf + 277);

    auto g_zzzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_isf + 278);

    auto g_zzzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_isf + 279);

    /// Set up components of auxilary buffer : ISG

    auto g_xxxxxx_0_xxxx_1 = pbuffer.data(idx_eri_1_isg);

    auto g_xxxxxx_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 1);

    auto g_xxxxxx_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 2);

    auto g_xxxxxx_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 3);

    auto g_xxxxxx_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 4);

    auto g_xxxxxx_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 5);

    auto g_xxxxxx_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 6);

    auto g_xxxxxx_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 7);

    auto g_xxxxxx_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 8);

    auto g_xxxxxx_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 9);

    auto g_xxxxxx_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 10);

    auto g_xxxxxx_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 11);

    auto g_xxxxxx_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 12);

    auto g_xxxxxx_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 13);

    auto g_xxxxxx_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 14);

    auto g_xxxxxy_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 15);

    auto g_xxxxxy_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 16);

    auto g_xxxxxy_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 17);

    auto g_xxxxxy_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 18);

    auto g_xxxxxy_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 20);

    auto g_xxxxxy_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 21);

    auto g_xxxxxy_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 24);

    auto g_xxxxxy_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 25);

    auto g_xxxxxz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 30);

    auto g_xxxxxz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 31);

    auto g_xxxxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 32);

    auto g_xxxxxz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 33);

    auto g_xxxxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 34);

    auto g_xxxxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 35);

    auto g_xxxxxz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 36);

    auto g_xxxxxz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 37);

    auto g_xxxxxz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 38);

    auto g_xxxxxz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 39);

    auto g_xxxxxz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 41);

    auto g_xxxxxz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 42);

    auto g_xxxxxz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 43);

    auto g_xxxxxz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 44);

    auto g_xxxxyy_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 45);

    auto g_xxxxyy_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 46);

    auto g_xxxxyy_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 47);

    auto g_xxxxyy_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 48);

    auto g_xxxxyy_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 49);

    auto g_xxxxyy_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 50);

    auto g_xxxxyy_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 51);

    auto g_xxxxyy_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 52);

    auto g_xxxxyy_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 53);

    auto g_xxxxyy_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 54);

    auto g_xxxxyy_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 55);

    auto g_xxxxyy_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 56);

    auto g_xxxxyy_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 57);

    auto g_xxxxyy_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 58);

    auto g_xxxxyy_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 59);

    auto g_xxxxzz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 75);

    auto g_xxxxzz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 76);

    auto g_xxxxzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 77);

    auto g_xxxxzz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 78);

    auto g_xxxxzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 79);

    auto g_xxxxzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 80);

    auto g_xxxxzz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 81);

    auto g_xxxxzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 82);

    auto g_xxxxzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 83);

    auto g_xxxxzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 84);

    auto g_xxxxzz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 85);

    auto g_xxxxzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 86);

    auto g_xxxxzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 87);

    auto g_xxxxzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 88);

    auto g_xxxxzz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 89);

    auto g_xxxyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 90);

    auto g_xxxyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 91);

    auto g_xxxyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 92);

    auto g_xxxyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 93);

    auto g_xxxyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 94);

    auto g_xxxyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 95);

    auto g_xxxyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 96);

    auto g_xxxyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 97);

    auto g_xxxyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 98);

    auto g_xxxyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 99);

    auto g_xxxyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 100);

    auto g_xxxyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 101);

    auto g_xxxyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 102);

    auto g_xxxyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 103);

    auto g_xxxyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 104);

    auto g_xxxyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 106);

    auto g_xxxyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 108);

    auto g_xxxyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 111);

    auto g_xxxyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 120);

    auto g_xxxyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 122);

    auto g_xxxyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 125);

    auto g_xxxyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 129);

    auto g_xxxzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 135);

    auto g_xxxzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 136);

    auto g_xxxzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 137);

    auto g_xxxzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 138);

    auto g_xxxzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 139);

    auto g_xxxzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 140);

    auto g_xxxzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 141);

    auto g_xxxzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 142);

    auto g_xxxzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 143);

    auto g_xxxzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 144);

    auto g_xxxzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 145);

    auto g_xxxzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 146);

    auto g_xxxzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 147);

    auto g_xxxzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 148);

    auto g_xxxzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 149);

    auto g_xxyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 150);

    auto g_xxyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 151);

    auto g_xxyyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 152);

    auto g_xxyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 153);

    auto g_xxyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 154);

    auto g_xxyyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 155);

    auto g_xxyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 156);

    auto g_xxyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 157);

    auto g_xxyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 158);

    auto g_xxyyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 159);

    auto g_xxyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 160);

    auto g_xxyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 161);

    auto g_xxyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 162);

    auto g_xxyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 163);

    auto g_xxyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 164);

    auto g_xxyyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 166);

    auto g_xxyyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 168);

    auto g_xxyyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 171);

    auto g_xxyyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 180);

    auto g_xxyyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 181);

    auto g_xxyyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 182);

    auto g_xxyyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 183);

    auto g_xxyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 184);

    auto g_xxyyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 185);

    auto g_xxyyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 186);

    auto g_xxyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 187);

    auto g_xxyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 188);

    auto g_xxyyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 189);

    auto g_xxyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 190);

    auto g_xxyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 191);

    auto g_xxyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 192);

    auto g_xxyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 193);

    auto g_xxyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 194);

    auto g_xxyzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 195);

    auto g_xxyzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 197);

    auto g_xxyzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 200);

    auto g_xxyzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 204);

    auto g_xxzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 210);

    auto g_xxzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 211);

    auto g_xxzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 212);

    auto g_xxzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 213);

    auto g_xxzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 214);

    auto g_xxzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 215);

    auto g_xxzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 216);

    auto g_xxzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 217);

    auto g_xxzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 218);

    auto g_xxzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 219);

    auto g_xxzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 220);

    auto g_xxzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 221);

    auto g_xxzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 222);

    auto g_xxzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 223);

    auto g_xxzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 224);

    auto g_xyyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 225);

    auto g_xyyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 226);

    auto g_xyyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 228);

    auto g_xyyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 229);

    auto g_xyyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 231);

    auto g_xyyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 232);

    auto g_xyyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 233);

    auto g_xyyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 235);

    auto g_xyyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 236);

    auto g_xyyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 237);

    auto g_xyyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 238);

    auto g_xyyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 239);

    auto g_xyyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 259);

    auto g_xyyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 262);

    auto g_xyyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 263);

    auto g_xyyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 265);

    auto g_xyyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 266);

    auto g_xyyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 267);

    auto g_xyyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 268);

    auto g_xyyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 269);

    auto g_xyyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 274);

    auto g_xyyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 277);

    auto g_xyyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 278);

    auto g_xyyzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 280);

    auto g_xyyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 281);

    auto g_xyyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 282);

    auto g_xyyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 283);

    auto g_xyyzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 284);

    auto g_xzzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 300);

    auto g_xzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 302);

    auto g_xzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 304);

    auto g_xzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 305);

    auto g_xzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 307);

    auto g_xzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 308);

    auto g_xzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 309);

    auto g_xzzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 310);

    auto g_xzzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 311);

    auto g_xzzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 312);

    auto g_xzzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 313);

    auto g_xzzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 314);

    auto g_yyyyyy_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 315);

    auto g_yyyyyy_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 316);

    auto g_yyyyyy_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 317);

    auto g_yyyyyy_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 318);

    auto g_yyyyyy_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 319);

    auto g_yyyyyy_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 320);

    auto g_yyyyyy_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 321);

    auto g_yyyyyy_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 322);

    auto g_yyyyyy_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 323);

    auto g_yyyyyy_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 324);

    auto g_yyyyyy_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 325);

    auto g_yyyyyy_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 326);

    auto g_yyyyyy_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 327);

    auto g_yyyyyy_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 328);

    auto g_yyyyyy_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 329);

    auto g_yyyyyz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 331);

    auto g_yyyyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 332);

    auto g_yyyyyz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 333);

    auto g_yyyyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 334);

    auto g_yyyyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 335);

    auto g_yyyyyz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 336);

    auto g_yyyyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 337);

    auto g_yyyyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 338);

    auto g_yyyyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 339);

    auto g_yyyyyz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 340);

    auto g_yyyyyz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 341);

    auto g_yyyyyz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 342);

    auto g_yyyyyz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 343);

    auto g_yyyyyz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 344);

    auto g_yyyyzz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 345);

    auto g_yyyyzz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 346);

    auto g_yyyyzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 347);

    auto g_yyyyzz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 348);

    auto g_yyyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 349);

    auto g_yyyyzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 350);

    auto g_yyyyzz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 351);

    auto g_yyyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 352);

    auto g_yyyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 353);

    auto g_yyyyzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 354);

    auto g_yyyyzz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 355);

    auto g_yyyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 356);

    auto g_yyyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 357);

    auto g_yyyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 358);

    auto g_yyyyzz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 359);

    auto g_yyyzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 360);

    auto g_yyyzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 361);

    auto g_yyyzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 362);

    auto g_yyyzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 363);

    auto g_yyyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 364);

    auto g_yyyzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 365);

    auto g_yyyzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 366);

    auto g_yyyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 367);

    auto g_yyyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 368);

    auto g_yyyzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 369);

    auto g_yyyzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 370);

    auto g_yyyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 371);

    auto g_yyyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 372);

    auto g_yyyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 373);

    auto g_yyyzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 374);

    auto g_yyzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 375);

    auto g_yyzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 376);

    auto g_yyzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 377);

    auto g_yyzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 378);

    auto g_yyzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 379);

    auto g_yyzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 380);

    auto g_yyzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 381);

    auto g_yyzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 382);

    auto g_yyzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 383);

    auto g_yyzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 384);

    auto g_yyzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 385);

    auto g_yyzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 386);

    auto g_yyzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 387);

    auto g_yyzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 388);

    auto g_yyzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 389);

    auto g_yzzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 390);

    auto g_yzzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 391);

    auto g_yzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 392);

    auto g_yzzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 393);

    auto g_yzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 394);

    auto g_yzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 395);

    auto g_yzzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 396);

    auto g_yzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 397);

    auto g_yzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 398);

    auto g_yzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 399);

    auto g_yzzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 400);

    auto g_yzzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 401);

    auto g_yzzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 402);

    auto g_yzzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 403);

    auto g_yzzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 404);

    auto g_zzzzzz_0_xxxx_1 = pbuffer.data(idx_eri_1_isg + 405);

    auto g_zzzzzz_0_xxxy_1 = pbuffer.data(idx_eri_1_isg + 406);

    auto g_zzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 407);

    auto g_zzzzzz_0_xxyy_1 = pbuffer.data(idx_eri_1_isg + 408);

    auto g_zzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 409);

    auto g_zzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 410);

    auto g_zzzzzz_0_xyyy_1 = pbuffer.data(idx_eri_1_isg + 411);

    auto g_zzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 412);

    auto g_zzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 413);

    auto g_zzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 414);

    auto g_zzzzzz_0_yyyy_1 = pbuffer.data(idx_eri_1_isg + 415);

    auto g_zzzzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 416);

    auto g_zzzzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 417);

    auto g_zzzzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 418);

    auto g_zzzzzz_0_zzzz_1 = pbuffer.data(idx_eri_1_isg + 419);

    /// Set up 0-15 components of targeted buffer : KSG

    auto g_xxxxxxx_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg);

    auto g_xxxxxxx_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 1);

    auto g_xxxxxxx_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 2);

    auto g_xxxxxxx_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 3);

    auto g_xxxxxxx_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 4);

    auto g_xxxxxxx_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 5);

    auto g_xxxxxxx_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 6);

    auto g_xxxxxxx_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 7);

    auto g_xxxxxxx_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 8);

    auto g_xxxxxxx_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 9);

    auto g_xxxxxxx_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 10);

    auto g_xxxxxxx_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 11);

    auto g_xxxxxxx_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 12);

    auto g_xxxxxxx_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 13);

    auto g_xxxxxxx_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 14);

    #pragma omp simd aligned(g_xxxxx_0_xxxx_0, g_xxxxx_0_xxxx_1, g_xxxxx_0_xxxy_0, g_xxxxx_0_xxxy_1, g_xxxxx_0_xxxz_0, g_xxxxx_0_xxxz_1, g_xxxxx_0_xxyy_0, g_xxxxx_0_xxyy_1, g_xxxxx_0_xxyz_0, g_xxxxx_0_xxyz_1, g_xxxxx_0_xxzz_0, g_xxxxx_0_xxzz_1, g_xxxxx_0_xyyy_0, g_xxxxx_0_xyyy_1, g_xxxxx_0_xyyz_0, g_xxxxx_0_xyyz_1, g_xxxxx_0_xyzz_0, g_xxxxx_0_xyzz_1, g_xxxxx_0_xzzz_0, g_xxxxx_0_xzzz_1, g_xxxxx_0_yyyy_0, g_xxxxx_0_yyyy_1, g_xxxxx_0_yyyz_0, g_xxxxx_0_yyyz_1, g_xxxxx_0_yyzz_0, g_xxxxx_0_yyzz_1, g_xxxxx_0_yzzz_0, g_xxxxx_0_yzzz_1, g_xxxxx_0_zzzz_0, g_xxxxx_0_zzzz_1, g_xxxxxx_0_xxx_1, g_xxxxxx_0_xxxx_1, g_xxxxxx_0_xxxy_1, g_xxxxxx_0_xxxz_1, g_xxxxxx_0_xxy_1, g_xxxxxx_0_xxyy_1, g_xxxxxx_0_xxyz_1, g_xxxxxx_0_xxz_1, g_xxxxxx_0_xxzz_1, g_xxxxxx_0_xyy_1, g_xxxxxx_0_xyyy_1, g_xxxxxx_0_xyyz_1, g_xxxxxx_0_xyz_1, g_xxxxxx_0_xyzz_1, g_xxxxxx_0_xzz_1, g_xxxxxx_0_xzzz_1, g_xxxxxx_0_yyy_1, g_xxxxxx_0_yyyy_1, g_xxxxxx_0_yyyz_1, g_xxxxxx_0_yyz_1, g_xxxxxx_0_yyzz_1, g_xxxxxx_0_yzz_1, g_xxxxxx_0_yzzz_1, g_xxxxxx_0_zzz_1, g_xxxxxx_0_zzzz_1, g_xxxxxxx_0_xxxx_0, g_xxxxxxx_0_xxxy_0, g_xxxxxxx_0_xxxz_0, g_xxxxxxx_0_xxyy_0, g_xxxxxxx_0_xxyz_0, g_xxxxxxx_0_xxzz_0, g_xxxxxxx_0_xyyy_0, g_xxxxxxx_0_xyyz_0, g_xxxxxxx_0_xyzz_0, g_xxxxxxx_0_xzzz_0, g_xxxxxxx_0_yyyy_0, g_xxxxxxx_0_yyyz_0, g_xxxxxxx_0_yyzz_0, g_xxxxxxx_0_yzzz_0, g_xxxxxxx_0_zzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxxx_0_xxxx_0[i] = 6.0 * g_xxxxx_0_xxxx_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxx_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxx_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxy_0[i] = 6.0 * g_xxxxx_0_xxxy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxz_0[i] = 6.0 * g_xxxxx_0_xxxz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyy_0[i] = 6.0 * g_xxxxx_0_xxyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyz_0[i] = 6.0 * g_xxxxx_0_xxyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxzz_0[i] = 6.0 * g_xxxxx_0_xxzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyy_0[i] = 6.0 * g_xxxxx_0_xyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyy_1[i] * fz_be_0 + g_xxxxxx_0_yyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyz_0[i] = 6.0 * g_xxxxx_0_xyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyz_1[i] * fz_be_0 + g_xxxxxx_0_yyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyzz_0[i] = 6.0 * g_xxxxx_0_xyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyzz_1[i] * fz_be_0 + g_xxxxxx_0_yzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xzzz_0[i] = 6.0 * g_xxxxx_0_xzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xzzz_1[i] * fz_be_0 + g_xxxxxx_0_zzz_1[i] * fi_acd_0 + g_xxxxxx_0_xzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyy_0[i] = 6.0 * g_xxxxx_0_yyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyy_1[i] * fz_be_0 + g_xxxxxx_0_yyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyz_0[i] = 6.0 * g_xxxxx_0_yyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyz_1[i] * fz_be_0 + g_xxxxxx_0_yyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyzz_0[i] = 6.0 * g_xxxxx_0_yyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyzz_1[i] * fz_be_0 + g_xxxxxx_0_yyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yzzz_0[i] = 6.0 * g_xxxxx_0_yzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yzzz_1[i] * fz_be_0 + g_xxxxxx_0_yzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_zzzz_0[i] = 6.0 * g_xxxxx_0_zzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_zzzz_1[i] * fz_be_0 + g_xxxxxx_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 15-30 components of targeted buffer : KSG

    auto g_xxxxxxy_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 15);

    auto g_xxxxxxy_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 16);

    auto g_xxxxxxy_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 17);

    auto g_xxxxxxy_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 18);

    auto g_xxxxxxy_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 19);

    auto g_xxxxxxy_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 20);

    auto g_xxxxxxy_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 21);

    auto g_xxxxxxy_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 22);

    auto g_xxxxxxy_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 23);

    auto g_xxxxxxy_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 24);

    auto g_xxxxxxy_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 25);

    auto g_xxxxxxy_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 26);

    auto g_xxxxxxy_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 27);

    auto g_xxxxxxy_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 28);

    auto g_xxxxxxy_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 29);

    #pragma omp simd aligned(g_xxxxxx_0_xxx_1, g_xxxxxx_0_xxxx_1, g_xxxxxx_0_xxxy_1, g_xxxxxx_0_xxxz_1, g_xxxxxx_0_xxy_1, g_xxxxxx_0_xxyy_1, g_xxxxxx_0_xxyz_1, g_xxxxxx_0_xxz_1, g_xxxxxx_0_xxzz_1, g_xxxxxx_0_xyy_1, g_xxxxxx_0_xyyy_1, g_xxxxxx_0_xyyz_1, g_xxxxxx_0_xyz_1, g_xxxxxx_0_xyzz_1, g_xxxxxx_0_xzz_1, g_xxxxxx_0_xzzz_1, g_xxxxxx_0_yyy_1, g_xxxxxx_0_yyyy_1, g_xxxxxx_0_yyyz_1, g_xxxxxx_0_yyz_1, g_xxxxxx_0_yyzz_1, g_xxxxxx_0_yzz_1, g_xxxxxx_0_yzzz_1, g_xxxxxx_0_zzz_1, g_xxxxxx_0_zzzz_1, g_xxxxxxy_0_xxxx_0, g_xxxxxxy_0_xxxy_0, g_xxxxxxy_0_xxxz_0, g_xxxxxxy_0_xxyy_0, g_xxxxxxy_0_xxyz_0, g_xxxxxxy_0_xxzz_0, g_xxxxxxy_0_xyyy_0, g_xxxxxxy_0_xyyz_0, g_xxxxxxy_0_xyzz_0, g_xxxxxxy_0_xzzz_0, g_xxxxxxy_0_yyyy_0, g_xxxxxxy_0_yyyz_0, g_xxxxxxy_0_yyzz_0, g_xxxxxxy_0_yzzz_0, g_xxxxxxy_0_zzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxy_0_xxxx_0[i] = g_xxxxxx_0_xxxx_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxy_0[i] = g_xxxxxx_0_xxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxz_0[i] = g_xxxxxx_0_xxxz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyy_0[i] = 2.0 * g_xxxxxx_0_xxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyz_0[i] = g_xxxxxx_0_xxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxzz_0[i] = g_xxxxxx_0_xxzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyy_0[i] = 3.0 * g_xxxxxx_0_xyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyz_0[i] = 2.0 * g_xxxxxx_0_xyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyzz_0[i] = g_xxxxxx_0_xzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xzzz_0[i] = g_xxxxxx_0_xzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyy_0[i] = 4.0 * g_xxxxxx_0_yyy_1[i] * fi_acd_0 + g_xxxxxx_0_yyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyz_0[i] = 3.0 * g_xxxxxx_0_yyz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyzz_0[i] = 2.0 * g_xxxxxx_0_yzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yzzz_0[i] = g_xxxxxx_0_zzz_1[i] * fi_acd_0 + g_xxxxxx_0_yzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_zzzz_0[i] = g_xxxxxx_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 30-45 components of targeted buffer : KSG

    auto g_xxxxxxz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 30);

    auto g_xxxxxxz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 31);

    auto g_xxxxxxz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 32);

    auto g_xxxxxxz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 33);

    auto g_xxxxxxz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 34);

    auto g_xxxxxxz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 35);

    auto g_xxxxxxz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 36);

    auto g_xxxxxxz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 37);

    auto g_xxxxxxz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 38);

    auto g_xxxxxxz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 39);

    auto g_xxxxxxz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 40);

    auto g_xxxxxxz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 41);

    auto g_xxxxxxz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 42);

    auto g_xxxxxxz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 43);

    auto g_xxxxxxz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 44);

    #pragma omp simd aligned(g_xxxxxx_0_xxx_1, g_xxxxxx_0_xxxx_1, g_xxxxxx_0_xxxy_1, g_xxxxxx_0_xxxz_1, g_xxxxxx_0_xxy_1, g_xxxxxx_0_xxyy_1, g_xxxxxx_0_xxyz_1, g_xxxxxx_0_xxz_1, g_xxxxxx_0_xxzz_1, g_xxxxxx_0_xyy_1, g_xxxxxx_0_xyyy_1, g_xxxxxx_0_xyyz_1, g_xxxxxx_0_xyz_1, g_xxxxxx_0_xyzz_1, g_xxxxxx_0_xzz_1, g_xxxxxx_0_xzzz_1, g_xxxxxx_0_yyy_1, g_xxxxxx_0_yyyy_1, g_xxxxxx_0_yyyz_1, g_xxxxxx_0_yyz_1, g_xxxxxx_0_yyzz_1, g_xxxxxx_0_yzz_1, g_xxxxxx_0_yzzz_1, g_xxxxxx_0_zzz_1, g_xxxxxx_0_zzzz_1, g_xxxxxxz_0_xxxx_0, g_xxxxxxz_0_xxxy_0, g_xxxxxxz_0_xxxz_0, g_xxxxxxz_0_xxyy_0, g_xxxxxxz_0_xxyz_0, g_xxxxxxz_0_xxzz_0, g_xxxxxxz_0_xyyy_0, g_xxxxxxz_0_xyyz_0, g_xxxxxxz_0_xyzz_0, g_xxxxxxz_0_xzzz_0, g_xxxxxxz_0_yyyy_0, g_xxxxxxz_0_yyyz_0, g_xxxxxxz_0_yyzz_0, g_xxxxxxz_0_yzzz_0, g_xxxxxxz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxz_0_xxxx_0[i] = g_xxxxxx_0_xxxx_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxy_0[i] = g_xxxxxx_0_xxxy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxz_0[i] = g_xxxxxx_0_xxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyy_0[i] = g_xxxxxx_0_xxyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyz_0[i] = g_xxxxxx_0_xxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxzz_0[i] = 2.0 * g_xxxxxx_0_xxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyy_0[i] = g_xxxxxx_0_xyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyz_0[i] = g_xxxxxx_0_xyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyzz_0[i] = 2.0 * g_xxxxxx_0_xyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xzzz_0[i] = 3.0 * g_xxxxxx_0_xzz_1[i] * fi_acd_0 + g_xxxxxx_0_xzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyy_0[i] = g_xxxxxx_0_yyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyz_0[i] = g_xxxxxx_0_yyy_1[i] * fi_acd_0 + g_xxxxxx_0_yyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyzz_0[i] = 2.0 * g_xxxxxx_0_yyz_1[i] * fi_acd_0 + g_xxxxxx_0_yyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yzzz_0[i] = 3.0 * g_xxxxxx_0_yzz_1[i] * fi_acd_0 + g_xxxxxx_0_yzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_zzzz_0[i] = 4.0 * g_xxxxxx_0_zzz_1[i] * fi_acd_0 + g_xxxxxx_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 45-60 components of targeted buffer : KSG

    auto g_xxxxxyy_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 45);

    auto g_xxxxxyy_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 46);

    auto g_xxxxxyy_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 47);

    auto g_xxxxxyy_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 48);

    auto g_xxxxxyy_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 49);

    auto g_xxxxxyy_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 50);

    auto g_xxxxxyy_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 51);

    auto g_xxxxxyy_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 52);

    auto g_xxxxxyy_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 53);

    auto g_xxxxxyy_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 54);

    auto g_xxxxxyy_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 55);

    auto g_xxxxxyy_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 56);

    auto g_xxxxxyy_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 57);

    auto g_xxxxxyy_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 58);

    auto g_xxxxxyy_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 59);

    #pragma omp simd aligned(g_xxxxx_0_xxxx_0, g_xxxxx_0_xxxx_1, g_xxxxx_0_xxxz_0, g_xxxxx_0_xxxz_1, g_xxxxx_0_xxzz_0, g_xxxxx_0_xxzz_1, g_xxxxx_0_xzzz_0, g_xxxxx_0_xzzz_1, g_xxxxxy_0_xxxx_1, g_xxxxxy_0_xxxz_1, g_xxxxxy_0_xxzz_1, g_xxxxxy_0_xzzz_1, g_xxxxxyy_0_xxxx_0, g_xxxxxyy_0_xxxy_0, g_xxxxxyy_0_xxxz_0, g_xxxxxyy_0_xxyy_0, g_xxxxxyy_0_xxyz_0, g_xxxxxyy_0_xxzz_0, g_xxxxxyy_0_xyyy_0, g_xxxxxyy_0_xyyz_0, g_xxxxxyy_0_xyzz_0, g_xxxxxyy_0_xzzz_0, g_xxxxxyy_0_yyyy_0, g_xxxxxyy_0_yyyz_0, g_xxxxxyy_0_yyzz_0, g_xxxxxyy_0_yzzz_0, g_xxxxxyy_0_zzzz_0, g_xxxxyy_0_xxxy_1, g_xxxxyy_0_xxy_1, g_xxxxyy_0_xxyy_1, g_xxxxyy_0_xxyz_1, g_xxxxyy_0_xyy_1, g_xxxxyy_0_xyyy_1, g_xxxxyy_0_xyyz_1, g_xxxxyy_0_xyz_1, g_xxxxyy_0_xyzz_1, g_xxxxyy_0_yyy_1, g_xxxxyy_0_yyyy_1, g_xxxxyy_0_yyyz_1, g_xxxxyy_0_yyz_1, g_xxxxyy_0_yyzz_1, g_xxxxyy_0_yzz_1, g_xxxxyy_0_yzzz_1, g_xxxxyy_0_zzzz_1, g_xxxyy_0_xxxy_0, g_xxxyy_0_xxxy_1, g_xxxyy_0_xxyy_0, g_xxxyy_0_xxyy_1, g_xxxyy_0_xxyz_0, g_xxxyy_0_xxyz_1, g_xxxyy_0_xyyy_0, g_xxxyy_0_xyyy_1, g_xxxyy_0_xyyz_0, g_xxxyy_0_xyyz_1, g_xxxyy_0_xyzz_0, g_xxxyy_0_xyzz_1, g_xxxyy_0_yyyy_0, g_xxxyy_0_yyyy_1, g_xxxyy_0_yyyz_0, g_xxxyy_0_yyyz_1, g_xxxyy_0_yyzz_0, g_xxxyy_0_yyzz_1, g_xxxyy_0_yzzz_0, g_xxxyy_0_yzzz_1, g_xxxyy_0_zzzz_0, g_xxxyy_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxyy_0_xxxx_0[i] = g_xxxxx_0_xxxx_0[i] * fbe_0 - g_xxxxx_0_xxxx_1[i] * fz_be_0 + g_xxxxxy_0_xxxx_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxy_0[i] = 4.0 * g_xxxyy_0_xxxy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxz_0[i] = g_xxxxx_0_xxxz_0[i] * fbe_0 - g_xxxxx_0_xxxz_1[i] * fz_be_0 + g_xxxxxy_0_xxxz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxyy_0[i] = 4.0 * g_xxxyy_0_xxyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyz_0[i] = 4.0 * g_xxxyy_0_xxyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxzz_0[i] = g_xxxxx_0_xxzz_0[i] * fbe_0 - g_xxxxx_0_xxzz_1[i] * fz_be_0 + g_xxxxxy_0_xxzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xyyy_0[i] = 4.0 * g_xxxyy_0_xyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyy_1[i] * fz_be_0 + g_xxxxyy_0_yyy_1[i] * fi_acd_0 + g_xxxxyy_0_xyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyz_0[i] = 4.0 * g_xxxyy_0_xyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyz_1[i] * fz_be_0 + g_xxxxyy_0_yyz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyzz_0[i] = 4.0 * g_xxxyy_0_xyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyzz_1[i] * fz_be_0 + g_xxxxyy_0_yzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xzzz_0[i] = g_xxxxx_0_xzzz_0[i] * fbe_0 - g_xxxxx_0_xzzz_1[i] * fz_be_0 + g_xxxxxy_0_xzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_yyyy_0[i] = 4.0 * g_xxxyy_0_yyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyy_1[i] * fz_be_0 + g_xxxxyy_0_yyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyz_0[i] = 4.0 * g_xxxyy_0_yyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyz_1[i] * fz_be_0 + g_xxxxyy_0_yyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyzz_0[i] = 4.0 * g_xxxyy_0_yyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyzz_1[i] * fz_be_0 + g_xxxxyy_0_yyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yzzz_0[i] = 4.0 * g_xxxyy_0_yzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yzzz_1[i] * fz_be_0 + g_xxxxyy_0_yzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_zzzz_0[i] = 4.0 * g_xxxyy_0_zzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_zzzz_1[i] * fz_be_0 + g_xxxxyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 60-75 components of targeted buffer : KSG

    auto g_xxxxxyz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 60);

    auto g_xxxxxyz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 61);

    auto g_xxxxxyz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 62);

    auto g_xxxxxyz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 63);

    auto g_xxxxxyz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 64);

    auto g_xxxxxyz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 65);

    auto g_xxxxxyz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 66);

    auto g_xxxxxyz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 67);

    auto g_xxxxxyz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 68);

    auto g_xxxxxyz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 69);

    auto g_xxxxxyz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 70);

    auto g_xxxxxyz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 71);

    auto g_xxxxxyz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 72);

    auto g_xxxxxyz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 73);

    auto g_xxxxxyz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 74);

    #pragma omp simd aligned(g_xxxxxy_0_xxxy_1, g_xxxxxy_0_xxyy_1, g_xxxxxy_0_xyyy_1, g_xxxxxy_0_yyyy_1, g_xxxxxyz_0_xxxx_0, g_xxxxxyz_0_xxxy_0, g_xxxxxyz_0_xxxz_0, g_xxxxxyz_0_xxyy_0, g_xxxxxyz_0_xxyz_0, g_xxxxxyz_0_xxzz_0, g_xxxxxyz_0_xyyy_0, g_xxxxxyz_0_xyyz_0, g_xxxxxyz_0_xyzz_0, g_xxxxxyz_0_xzzz_0, g_xxxxxyz_0_yyyy_0, g_xxxxxyz_0_yyyz_0, g_xxxxxyz_0_yyzz_0, g_xxxxxyz_0_yzzz_0, g_xxxxxyz_0_zzzz_0, g_xxxxxz_0_xxxx_1, g_xxxxxz_0_xxxz_1, g_xxxxxz_0_xxyz_1, g_xxxxxz_0_xxz_1, g_xxxxxz_0_xxzz_1, g_xxxxxz_0_xyyz_1, g_xxxxxz_0_xyz_1, g_xxxxxz_0_xyzz_1, g_xxxxxz_0_xzz_1, g_xxxxxz_0_xzzz_1, g_xxxxxz_0_yyyz_1, g_xxxxxz_0_yyz_1, g_xxxxxz_0_yyzz_1, g_xxxxxz_0_yzz_1, g_xxxxxz_0_yzzz_1, g_xxxxxz_0_zzz_1, g_xxxxxz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxyz_0_xxxx_0[i] = g_xxxxxz_0_xxxx_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxy_0[i] = g_xxxxxy_0_xxxy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxz_0[i] = g_xxxxxz_0_xxxz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyy_0[i] = g_xxxxxy_0_xxyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxyz_0[i] = g_xxxxxz_0_xxz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxzz_0[i] = g_xxxxxz_0_xxzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyy_0[i] = g_xxxxxy_0_xyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xyyz_0[i] = 2.0 * g_xxxxxz_0_xyz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyzz_0[i] = g_xxxxxz_0_xzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xzzz_0[i] = g_xxxxxz_0_xzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyy_0[i] = g_xxxxxy_0_yyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_yyyz_0[i] = 3.0 * g_xxxxxz_0_yyz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyzz_0[i] = 2.0 * g_xxxxxz_0_yzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yzzz_0[i] = g_xxxxxz_0_zzz_1[i] * fi_acd_0 + g_xxxxxz_0_yzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_zzzz_0[i] = g_xxxxxz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 75-90 components of targeted buffer : KSG

    auto g_xxxxxzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 75);

    auto g_xxxxxzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 76);

    auto g_xxxxxzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 77);

    auto g_xxxxxzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 78);

    auto g_xxxxxzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 79);

    auto g_xxxxxzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 80);

    auto g_xxxxxzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 81);

    auto g_xxxxxzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 82);

    auto g_xxxxxzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 83);

    auto g_xxxxxzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 84);

    auto g_xxxxxzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 85);

    auto g_xxxxxzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 86);

    auto g_xxxxxzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 87);

    auto g_xxxxxzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 88);

    auto g_xxxxxzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 89);

    #pragma omp simd aligned(g_xxxxx_0_xxxx_0, g_xxxxx_0_xxxx_1, g_xxxxx_0_xxxy_0, g_xxxxx_0_xxxy_1, g_xxxxx_0_xxyy_0, g_xxxxx_0_xxyy_1, g_xxxxx_0_xyyy_0, g_xxxxx_0_xyyy_1, g_xxxxxz_0_xxxx_1, g_xxxxxz_0_xxxy_1, g_xxxxxz_0_xxyy_1, g_xxxxxz_0_xyyy_1, g_xxxxxzz_0_xxxx_0, g_xxxxxzz_0_xxxy_0, g_xxxxxzz_0_xxxz_0, g_xxxxxzz_0_xxyy_0, g_xxxxxzz_0_xxyz_0, g_xxxxxzz_0_xxzz_0, g_xxxxxzz_0_xyyy_0, g_xxxxxzz_0_xyyz_0, g_xxxxxzz_0_xyzz_0, g_xxxxxzz_0_xzzz_0, g_xxxxxzz_0_yyyy_0, g_xxxxxzz_0_yyyz_0, g_xxxxxzz_0_yyzz_0, g_xxxxxzz_0_yzzz_0, g_xxxxxzz_0_zzzz_0, g_xxxxzz_0_xxxz_1, g_xxxxzz_0_xxyz_1, g_xxxxzz_0_xxz_1, g_xxxxzz_0_xxzz_1, g_xxxxzz_0_xyyz_1, g_xxxxzz_0_xyz_1, g_xxxxzz_0_xyzz_1, g_xxxxzz_0_xzz_1, g_xxxxzz_0_xzzz_1, g_xxxxzz_0_yyyy_1, g_xxxxzz_0_yyyz_1, g_xxxxzz_0_yyz_1, g_xxxxzz_0_yyzz_1, g_xxxxzz_0_yzz_1, g_xxxxzz_0_yzzz_1, g_xxxxzz_0_zzz_1, g_xxxxzz_0_zzzz_1, g_xxxzz_0_xxxz_0, g_xxxzz_0_xxxz_1, g_xxxzz_0_xxyz_0, g_xxxzz_0_xxyz_1, g_xxxzz_0_xxzz_0, g_xxxzz_0_xxzz_1, g_xxxzz_0_xyyz_0, g_xxxzz_0_xyyz_1, g_xxxzz_0_xyzz_0, g_xxxzz_0_xyzz_1, g_xxxzz_0_xzzz_0, g_xxxzz_0_xzzz_1, g_xxxzz_0_yyyy_0, g_xxxzz_0_yyyy_1, g_xxxzz_0_yyyz_0, g_xxxzz_0_yyyz_1, g_xxxzz_0_yyzz_0, g_xxxzz_0_yyzz_1, g_xxxzz_0_yzzz_0, g_xxxzz_0_yzzz_1, g_xxxzz_0_zzzz_0, g_xxxzz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxzz_0_xxxx_0[i] = g_xxxxx_0_xxxx_0[i] * fbe_0 - g_xxxxx_0_xxxx_1[i] * fz_be_0 + g_xxxxxz_0_xxxx_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxy_0[i] = g_xxxxx_0_xxxy_0[i] * fbe_0 - g_xxxxx_0_xxxy_1[i] * fz_be_0 + g_xxxxxz_0_xxxy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxz_0[i] = 4.0 * g_xxxzz_0_xxxz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyy_0[i] = g_xxxxx_0_xxyy_0[i] * fbe_0 - g_xxxxx_0_xxyy_1[i] * fz_be_0 + g_xxxxxz_0_xxyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxyz_0[i] = 4.0 * g_xxxzz_0_xxyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxzz_0[i] = 4.0 * g_xxxzz_0_xxzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyy_0[i] = g_xxxxx_0_xyyy_0[i] * fbe_0 - g_xxxxx_0_xyyy_1[i] * fz_be_0 + g_xxxxxz_0_xyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xyyz_0[i] = 4.0 * g_xxxzz_0_xyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyz_1[i] * fz_be_0 + g_xxxxzz_0_yyz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyzz_0[i] = 4.0 * g_xxxzz_0_xyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyzz_1[i] * fz_be_0 + g_xxxxzz_0_yzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xzzz_0[i] = 4.0 * g_xxxzz_0_xzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xzzz_1[i] * fz_be_0 + g_xxxxzz_0_zzz_1[i] * fi_acd_0 + g_xxxxzz_0_xzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyy_0[i] = 4.0 * g_xxxzz_0_yyyy_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyy_1[i] * fz_be_0 + g_xxxxzz_0_yyyy_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyz_0[i] = 4.0 * g_xxxzz_0_yyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyz_1[i] * fz_be_0 + g_xxxxzz_0_yyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyzz_0[i] = 4.0 * g_xxxzz_0_yyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyzz_1[i] * fz_be_0 + g_xxxxzz_0_yyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yzzz_0[i] = 4.0 * g_xxxzz_0_yzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yzzz_1[i] * fz_be_0 + g_xxxxzz_0_yzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_zzzz_0[i] = 4.0 * g_xxxzz_0_zzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_zzzz_1[i] * fz_be_0 + g_xxxxzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 90-105 components of targeted buffer : KSG

    auto g_xxxxyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 90);

    auto g_xxxxyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 91);

    auto g_xxxxyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 92);

    auto g_xxxxyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 93);

    auto g_xxxxyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 94);

    auto g_xxxxyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 95);

    auto g_xxxxyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 96);

    auto g_xxxxyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 97);

    auto g_xxxxyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 98);

    auto g_xxxxyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 99);

    auto g_xxxxyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 100);

    auto g_xxxxyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 101);

    auto g_xxxxyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 102);

    auto g_xxxxyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 103);

    auto g_xxxxyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 104);

    #pragma omp simd aligned(g_xxxxy_0_xxxx_0, g_xxxxy_0_xxxx_1, g_xxxxy_0_xxxz_0, g_xxxxy_0_xxxz_1, g_xxxxy_0_xxzz_0, g_xxxxy_0_xxzz_1, g_xxxxy_0_xzzz_0, g_xxxxy_0_xzzz_1, g_xxxxyy_0_xxxx_1, g_xxxxyy_0_xxxz_1, g_xxxxyy_0_xxzz_1, g_xxxxyy_0_xzzz_1, g_xxxxyyy_0_xxxx_0, g_xxxxyyy_0_xxxy_0, g_xxxxyyy_0_xxxz_0, g_xxxxyyy_0_xxyy_0, g_xxxxyyy_0_xxyz_0, g_xxxxyyy_0_xxzz_0, g_xxxxyyy_0_xyyy_0, g_xxxxyyy_0_xyyz_0, g_xxxxyyy_0_xyzz_0, g_xxxxyyy_0_xzzz_0, g_xxxxyyy_0_yyyy_0, g_xxxxyyy_0_yyyz_0, g_xxxxyyy_0_yyzz_0, g_xxxxyyy_0_yzzz_0, g_xxxxyyy_0_zzzz_0, g_xxxyyy_0_xxxy_1, g_xxxyyy_0_xxy_1, g_xxxyyy_0_xxyy_1, g_xxxyyy_0_xxyz_1, g_xxxyyy_0_xyy_1, g_xxxyyy_0_xyyy_1, g_xxxyyy_0_xyyz_1, g_xxxyyy_0_xyz_1, g_xxxyyy_0_xyzz_1, g_xxxyyy_0_yyy_1, g_xxxyyy_0_yyyy_1, g_xxxyyy_0_yyyz_1, g_xxxyyy_0_yyz_1, g_xxxyyy_0_yyzz_1, g_xxxyyy_0_yzz_1, g_xxxyyy_0_yzzz_1, g_xxxyyy_0_zzzz_1, g_xxyyy_0_xxxy_0, g_xxyyy_0_xxxy_1, g_xxyyy_0_xxyy_0, g_xxyyy_0_xxyy_1, g_xxyyy_0_xxyz_0, g_xxyyy_0_xxyz_1, g_xxyyy_0_xyyy_0, g_xxyyy_0_xyyy_1, g_xxyyy_0_xyyz_0, g_xxyyy_0_xyyz_1, g_xxyyy_0_xyzz_0, g_xxyyy_0_xyzz_1, g_xxyyy_0_yyyy_0, g_xxyyy_0_yyyy_1, g_xxyyy_0_yyyz_0, g_xxyyy_0_yyyz_1, g_xxyyy_0_yyzz_0, g_xxyyy_0_yyzz_1, g_xxyyy_0_yzzz_0, g_xxyyy_0_yzzz_1, g_xxyyy_0_zzzz_0, g_xxyyy_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyyy_0_xxxx_0[i] = 2.0 * g_xxxxy_0_xxxx_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxx_1[i] * fz_be_0 + g_xxxxyy_0_xxxx_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxy_0[i] = 3.0 * g_xxyyy_0_xxxy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxz_0[i] = 2.0 * g_xxxxy_0_xxxz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxz_1[i] * fz_be_0 + g_xxxxyy_0_xxxz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxyy_0[i] = 3.0 * g_xxyyy_0_xxyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyz_0[i] = 3.0 * g_xxyyy_0_xxyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxzz_0[i] = 2.0 * g_xxxxy_0_xxzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxzz_1[i] * fz_be_0 + g_xxxxyy_0_xxzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xyyy_0[i] = 3.0 * g_xxyyy_0_xyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyy_1[i] * fz_be_0 + g_xxxyyy_0_yyy_1[i] * fi_acd_0 + g_xxxyyy_0_xyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyz_0[i] = 3.0 * g_xxyyy_0_xyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyz_1[i] * fz_be_0 + g_xxxyyy_0_yyz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyzz_0[i] = 3.0 * g_xxyyy_0_xyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyzz_1[i] * fz_be_0 + g_xxxyyy_0_yzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xzzz_0[i] = 2.0 * g_xxxxy_0_xzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xzzz_1[i] * fz_be_0 + g_xxxxyy_0_xzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_yyyy_0[i] = 3.0 * g_xxyyy_0_yyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyy_1[i] * fz_be_0 + g_xxxyyy_0_yyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyz_0[i] = 3.0 * g_xxyyy_0_yyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyz_1[i] * fz_be_0 + g_xxxyyy_0_yyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyzz_0[i] = 3.0 * g_xxyyy_0_yyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyzz_1[i] * fz_be_0 + g_xxxyyy_0_yyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yzzz_0[i] = 3.0 * g_xxyyy_0_yzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yzzz_1[i] * fz_be_0 + g_xxxyyy_0_yzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_zzzz_0[i] = 3.0 * g_xxyyy_0_zzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_zzzz_1[i] * fz_be_0 + g_xxxyyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 105-120 components of targeted buffer : KSG

    auto g_xxxxyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 105);

    auto g_xxxxyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 106);

    auto g_xxxxyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 107);

    auto g_xxxxyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 108);

    auto g_xxxxyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 109);

    auto g_xxxxyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 110);

    auto g_xxxxyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 111);

    auto g_xxxxyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 112);

    auto g_xxxxyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 113);

    auto g_xxxxyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 114);

    auto g_xxxxyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 115);

    auto g_xxxxyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 116);

    auto g_xxxxyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 117);

    auto g_xxxxyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 118);

    auto g_xxxxyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 119);

    #pragma omp simd aligned(g_xxxxyy_0_xxx_1, g_xxxxyy_0_xxxx_1, g_xxxxyy_0_xxxy_1, g_xxxxyy_0_xxxz_1, g_xxxxyy_0_xxy_1, g_xxxxyy_0_xxyy_1, g_xxxxyy_0_xxyz_1, g_xxxxyy_0_xxz_1, g_xxxxyy_0_xxzz_1, g_xxxxyy_0_xyy_1, g_xxxxyy_0_xyyy_1, g_xxxxyy_0_xyyz_1, g_xxxxyy_0_xyz_1, g_xxxxyy_0_xyzz_1, g_xxxxyy_0_xzz_1, g_xxxxyy_0_xzzz_1, g_xxxxyy_0_yyy_1, g_xxxxyy_0_yyyy_1, g_xxxxyy_0_yyyz_1, g_xxxxyy_0_yyz_1, g_xxxxyy_0_yyzz_1, g_xxxxyy_0_yzz_1, g_xxxxyy_0_yzzz_1, g_xxxxyy_0_zzz_1, g_xxxxyy_0_zzzz_1, g_xxxxyyz_0_xxxx_0, g_xxxxyyz_0_xxxy_0, g_xxxxyyz_0_xxxz_0, g_xxxxyyz_0_xxyy_0, g_xxxxyyz_0_xxyz_0, g_xxxxyyz_0_xxzz_0, g_xxxxyyz_0_xyyy_0, g_xxxxyyz_0_xyyz_0, g_xxxxyyz_0_xyzz_0, g_xxxxyyz_0_xzzz_0, g_xxxxyyz_0_yyyy_0, g_xxxxyyz_0_yyyz_0, g_xxxxyyz_0_yyzz_0, g_xxxxyyz_0_yzzz_0, g_xxxxyyz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyyz_0_xxxx_0[i] = g_xxxxyy_0_xxxx_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxy_0[i] = g_xxxxyy_0_xxxy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxz_0[i] = g_xxxxyy_0_xxx_1[i] * fi_acd_0 + g_xxxxyy_0_xxxz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyy_0[i] = g_xxxxyy_0_xxyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyz_0[i] = g_xxxxyy_0_xxy_1[i] * fi_acd_0 + g_xxxxyy_0_xxyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxzz_0[i] = 2.0 * g_xxxxyy_0_xxz_1[i] * fi_acd_0 + g_xxxxyy_0_xxzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyy_0[i] = g_xxxxyy_0_xyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyz_0[i] = g_xxxxyy_0_xyy_1[i] * fi_acd_0 + g_xxxxyy_0_xyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyzz_0[i] = 2.0 * g_xxxxyy_0_xyz_1[i] * fi_acd_0 + g_xxxxyy_0_xyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xzzz_0[i] = 3.0 * g_xxxxyy_0_xzz_1[i] * fi_acd_0 + g_xxxxyy_0_xzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyy_0[i] = g_xxxxyy_0_yyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyz_0[i] = g_xxxxyy_0_yyy_1[i] * fi_acd_0 + g_xxxxyy_0_yyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyzz_0[i] = 2.0 * g_xxxxyy_0_yyz_1[i] * fi_acd_0 + g_xxxxyy_0_yyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yzzz_0[i] = 3.0 * g_xxxxyy_0_yzz_1[i] * fi_acd_0 + g_xxxxyy_0_yzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_zzzz_0[i] = 4.0 * g_xxxxyy_0_zzz_1[i] * fi_acd_0 + g_xxxxyy_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 120-135 components of targeted buffer : KSG

    auto g_xxxxyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 120);

    auto g_xxxxyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 121);

    auto g_xxxxyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 122);

    auto g_xxxxyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 123);

    auto g_xxxxyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 124);

    auto g_xxxxyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 125);

    auto g_xxxxyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 126);

    auto g_xxxxyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 127);

    auto g_xxxxyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 128);

    auto g_xxxxyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 129);

    auto g_xxxxyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 130);

    auto g_xxxxyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 131);

    auto g_xxxxyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 132);

    auto g_xxxxyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 133);

    auto g_xxxxyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 134);

    #pragma omp simd aligned(g_xxxxyzz_0_xxxx_0, g_xxxxyzz_0_xxxy_0, g_xxxxyzz_0_xxxz_0, g_xxxxyzz_0_xxyy_0, g_xxxxyzz_0_xxyz_0, g_xxxxyzz_0_xxzz_0, g_xxxxyzz_0_xyyy_0, g_xxxxyzz_0_xyyz_0, g_xxxxyzz_0_xyzz_0, g_xxxxyzz_0_xzzz_0, g_xxxxyzz_0_yyyy_0, g_xxxxyzz_0_yyyz_0, g_xxxxyzz_0_yyzz_0, g_xxxxyzz_0_yzzz_0, g_xxxxyzz_0_zzzz_0, g_xxxxzz_0_xxx_1, g_xxxxzz_0_xxxx_1, g_xxxxzz_0_xxxy_1, g_xxxxzz_0_xxxz_1, g_xxxxzz_0_xxy_1, g_xxxxzz_0_xxyy_1, g_xxxxzz_0_xxyz_1, g_xxxxzz_0_xxz_1, g_xxxxzz_0_xxzz_1, g_xxxxzz_0_xyy_1, g_xxxxzz_0_xyyy_1, g_xxxxzz_0_xyyz_1, g_xxxxzz_0_xyz_1, g_xxxxzz_0_xyzz_1, g_xxxxzz_0_xzz_1, g_xxxxzz_0_xzzz_1, g_xxxxzz_0_yyy_1, g_xxxxzz_0_yyyy_1, g_xxxxzz_0_yyyz_1, g_xxxxzz_0_yyz_1, g_xxxxzz_0_yyzz_1, g_xxxxzz_0_yzz_1, g_xxxxzz_0_yzzz_1, g_xxxxzz_0_zzz_1, g_xxxxzz_0_zzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyzz_0_xxxx_0[i] = g_xxxxzz_0_xxxx_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxy_0[i] = g_xxxxzz_0_xxx_1[i] * fi_acd_0 + g_xxxxzz_0_xxxy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxz_0[i] = g_xxxxzz_0_xxxz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyy_0[i] = 2.0 * g_xxxxzz_0_xxy_1[i] * fi_acd_0 + g_xxxxzz_0_xxyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyz_0[i] = g_xxxxzz_0_xxz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxzz_0[i] = g_xxxxzz_0_xxzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyy_0[i] = 3.0 * g_xxxxzz_0_xyy_1[i] * fi_acd_0 + g_xxxxzz_0_xyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyz_0[i] = 2.0 * g_xxxxzz_0_xyz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyzz_0[i] = g_xxxxzz_0_xzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xzzz_0[i] = g_xxxxzz_0_xzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyy_0[i] = 4.0 * g_xxxxzz_0_yyy_1[i] * fi_acd_0 + g_xxxxzz_0_yyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyz_0[i] = 3.0 * g_xxxxzz_0_yyz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyzz_0[i] = 2.0 * g_xxxxzz_0_yzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yzzz_0[i] = g_xxxxzz_0_zzz_1[i] * fi_acd_0 + g_xxxxzz_0_yzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_zzzz_0[i] = g_xxxxzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 135-150 components of targeted buffer : KSG

    auto g_xxxxzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 135);

    auto g_xxxxzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 136);

    auto g_xxxxzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 137);

    auto g_xxxxzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 138);

    auto g_xxxxzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 139);

    auto g_xxxxzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 140);

    auto g_xxxxzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 141);

    auto g_xxxxzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 142);

    auto g_xxxxzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 143);

    auto g_xxxxzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 144);

    auto g_xxxxzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 145);

    auto g_xxxxzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 146);

    auto g_xxxxzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 147);

    auto g_xxxxzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 148);

    auto g_xxxxzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 149);

    #pragma omp simd aligned(g_xxxxz_0_xxxx_0, g_xxxxz_0_xxxx_1, g_xxxxz_0_xxxy_0, g_xxxxz_0_xxxy_1, g_xxxxz_0_xxyy_0, g_xxxxz_0_xxyy_1, g_xxxxz_0_xyyy_0, g_xxxxz_0_xyyy_1, g_xxxxzz_0_xxxx_1, g_xxxxzz_0_xxxy_1, g_xxxxzz_0_xxyy_1, g_xxxxzz_0_xyyy_1, g_xxxxzzz_0_xxxx_0, g_xxxxzzz_0_xxxy_0, g_xxxxzzz_0_xxxz_0, g_xxxxzzz_0_xxyy_0, g_xxxxzzz_0_xxyz_0, g_xxxxzzz_0_xxzz_0, g_xxxxzzz_0_xyyy_0, g_xxxxzzz_0_xyyz_0, g_xxxxzzz_0_xyzz_0, g_xxxxzzz_0_xzzz_0, g_xxxxzzz_0_yyyy_0, g_xxxxzzz_0_yyyz_0, g_xxxxzzz_0_yyzz_0, g_xxxxzzz_0_yzzz_0, g_xxxxzzz_0_zzzz_0, g_xxxzzz_0_xxxz_1, g_xxxzzz_0_xxyz_1, g_xxxzzz_0_xxz_1, g_xxxzzz_0_xxzz_1, g_xxxzzz_0_xyyz_1, g_xxxzzz_0_xyz_1, g_xxxzzz_0_xyzz_1, g_xxxzzz_0_xzz_1, g_xxxzzz_0_xzzz_1, g_xxxzzz_0_yyyy_1, g_xxxzzz_0_yyyz_1, g_xxxzzz_0_yyz_1, g_xxxzzz_0_yyzz_1, g_xxxzzz_0_yzz_1, g_xxxzzz_0_yzzz_1, g_xxxzzz_0_zzz_1, g_xxxzzz_0_zzzz_1, g_xxzzz_0_xxxz_0, g_xxzzz_0_xxxz_1, g_xxzzz_0_xxyz_0, g_xxzzz_0_xxyz_1, g_xxzzz_0_xxzz_0, g_xxzzz_0_xxzz_1, g_xxzzz_0_xyyz_0, g_xxzzz_0_xyyz_1, g_xxzzz_0_xyzz_0, g_xxzzz_0_xyzz_1, g_xxzzz_0_xzzz_0, g_xxzzz_0_xzzz_1, g_xxzzz_0_yyyy_0, g_xxzzz_0_yyyy_1, g_xxzzz_0_yyyz_0, g_xxzzz_0_yyyz_1, g_xxzzz_0_yyzz_0, g_xxzzz_0_yyzz_1, g_xxzzz_0_yzzz_0, g_xxzzz_0_yzzz_1, g_xxzzz_0_zzzz_0, g_xxzzz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzzz_0_xxxx_0[i] = 2.0 * g_xxxxz_0_xxxx_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxx_1[i] * fz_be_0 + g_xxxxzz_0_xxxx_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxy_0[i] = 2.0 * g_xxxxz_0_xxxy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxy_1[i] * fz_be_0 + g_xxxxzz_0_xxxy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxz_0[i] = 3.0 * g_xxzzz_0_xxxz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyy_0[i] = 2.0 * g_xxxxz_0_xxyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxyy_1[i] * fz_be_0 + g_xxxxzz_0_xxyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxyz_0[i] = 3.0 * g_xxzzz_0_xxyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxzz_0[i] = 3.0 * g_xxzzz_0_xxzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyy_0[i] = 2.0 * g_xxxxz_0_xyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xyyy_1[i] * fz_be_0 + g_xxxxzz_0_xyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xyyz_0[i] = 3.0 * g_xxzzz_0_xyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyz_1[i] * fz_be_0 + g_xxxzzz_0_yyz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyzz_0[i] = 3.0 * g_xxzzz_0_xyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyzz_1[i] * fz_be_0 + g_xxxzzz_0_yzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xzzz_0[i] = 3.0 * g_xxzzz_0_xzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xzzz_1[i] * fz_be_0 + g_xxxzzz_0_zzz_1[i] * fi_acd_0 + g_xxxzzz_0_xzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyy_0[i] = 3.0 * g_xxzzz_0_yyyy_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyy_1[i] * fz_be_0 + g_xxxzzz_0_yyyy_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyz_0[i] = 3.0 * g_xxzzz_0_yyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyz_1[i] * fz_be_0 + g_xxxzzz_0_yyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyzz_0[i] = 3.0 * g_xxzzz_0_yyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyzz_1[i] * fz_be_0 + g_xxxzzz_0_yyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yzzz_0[i] = 3.0 * g_xxzzz_0_yzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yzzz_1[i] * fz_be_0 + g_xxxzzz_0_yzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_zzzz_0[i] = 3.0 * g_xxzzz_0_zzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_zzzz_1[i] * fz_be_0 + g_xxxzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 150-165 components of targeted buffer : KSG

    auto g_xxxyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 150);

    auto g_xxxyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 151);

    auto g_xxxyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 152);

    auto g_xxxyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 153);

    auto g_xxxyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 154);

    auto g_xxxyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 155);

    auto g_xxxyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 156);

    auto g_xxxyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 157);

    auto g_xxxyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 158);

    auto g_xxxyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 159);

    auto g_xxxyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 160);

    auto g_xxxyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 161);

    auto g_xxxyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 162);

    auto g_xxxyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 163);

    auto g_xxxyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 164);

    #pragma omp simd aligned(g_xxxyy_0_xxxx_0, g_xxxyy_0_xxxx_1, g_xxxyy_0_xxxz_0, g_xxxyy_0_xxxz_1, g_xxxyy_0_xxzz_0, g_xxxyy_0_xxzz_1, g_xxxyy_0_xzzz_0, g_xxxyy_0_xzzz_1, g_xxxyyy_0_xxxx_1, g_xxxyyy_0_xxxz_1, g_xxxyyy_0_xxzz_1, g_xxxyyy_0_xzzz_1, g_xxxyyyy_0_xxxx_0, g_xxxyyyy_0_xxxy_0, g_xxxyyyy_0_xxxz_0, g_xxxyyyy_0_xxyy_0, g_xxxyyyy_0_xxyz_0, g_xxxyyyy_0_xxzz_0, g_xxxyyyy_0_xyyy_0, g_xxxyyyy_0_xyyz_0, g_xxxyyyy_0_xyzz_0, g_xxxyyyy_0_xzzz_0, g_xxxyyyy_0_yyyy_0, g_xxxyyyy_0_yyyz_0, g_xxxyyyy_0_yyzz_0, g_xxxyyyy_0_yzzz_0, g_xxxyyyy_0_zzzz_0, g_xxyyyy_0_xxxy_1, g_xxyyyy_0_xxy_1, g_xxyyyy_0_xxyy_1, g_xxyyyy_0_xxyz_1, g_xxyyyy_0_xyy_1, g_xxyyyy_0_xyyy_1, g_xxyyyy_0_xyyz_1, g_xxyyyy_0_xyz_1, g_xxyyyy_0_xyzz_1, g_xxyyyy_0_yyy_1, g_xxyyyy_0_yyyy_1, g_xxyyyy_0_yyyz_1, g_xxyyyy_0_yyz_1, g_xxyyyy_0_yyzz_1, g_xxyyyy_0_yzz_1, g_xxyyyy_0_yzzz_1, g_xxyyyy_0_zzzz_1, g_xyyyy_0_xxxy_0, g_xyyyy_0_xxxy_1, g_xyyyy_0_xxyy_0, g_xyyyy_0_xxyy_1, g_xyyyy_0_xxyz_0, g_xyyyy_0_xxyz_1, g_xyyyy_0_xyyy_0, g_xyyyy_0_xyyy_1, g_xyyyy_0_xyyz_0, g_xyyyy_0_xyyz_1, g_xyyyy_0_xyzz_0, g_xyyyy_0_xyzz_1, g_xyyyy_0_yyyy_0, g_xyyyy_0_yyyy_1, g_xyyyy_0_yyyz_0, g_xyyyy_0_yyyz_1, g_xyyyy_0_yyzz_0, g_xyyyy_0_yyzz_1, g_xyyyy_0_yzzz_0, g_xyyyy_0_yzzz_1, g_xyyyy_0_zzzz_0, g_xyyyy_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyyy_0_xxxx_0[i] = 3.0 * g_xxxyy_0_xxxx_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxx_1[i] * fz_be_0 + g_xxxyyy_0_xxxx_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxy_0[i] = 2.0 * g_xyyyy_0_xxxy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxz_0[i] = 3.0 * g_xxxyy_0_xxxz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxz_1[i] * fz_be_0 + g_xxxyyy_0_xxxz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxyy_0[i] = 2.0 * g_xyyyy_0_xxyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyz_0[i] = 2.0 * g_xyyyy_0_xxyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxzz_0[i] = 3.0 * g_xxxyy_0_xxzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxzz_1[i] * fz_be_0 + g_xxxyyy_0_xxzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xyyy_0[i] = 2.0 * g_xyyyy_0_xyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyy_1[i] * fz_be_0 + g_xxyyyy_0_yyy_1[i] * fi_acd_0 + g_xxyyyy_0_xyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyz_0[i] = 2.0 * g_xyyyy_0_xyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyz_1[i] * fz_be_0 + g_xxyyyy_0_yyz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyzz_0[i] = 2.0 * g_xyyyy_0_xyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyzz_1[i] * fz_be_0 + g_xxyyyy_0_yzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xzzz_0[i] = 3.0 * g_xxxyy_0_xzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xzzz_1[i] * fz_be_0 + g_xxxyyy_0_xzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_yyyy_0[i] = 2.0 * g_xyyyy_0_yyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyy_1[i] * fz_be_0 + g_xxyyyy_0_yyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyz_0[i] = 2.0 * g_xyyyy_0_yyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyz_1[i] * fz_be_0 + g_xxyyyy_0_yyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyzz_0[i] = 2.0 * g_xyyyy_0_yyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyzz_1[i] * fz_be_0 + g_xxyyyy_0_yyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yzzz_0[i] = 2.0 * g_xyyyy_0_yzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yzzz_1[i] * fz_be_0 + g_xxyyyy_0_yzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_zzzz_0[i] = 2.0 * g_xyyyy_0_zzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_zzzz_1[i] * fz_be_0 + g_xxyyyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 165-180 components of targeted buffer : KSG

    auto g_xxxyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 165);

    auto g_xxxyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 166);

    auto g_xxxyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 167);

    auto g_xxxyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 168);

    auto g_xxxyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 169);

    auto g_xxxyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 170);

    auto g_xxxyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 171);

    auto g_xxxyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 172);

    auto g_xxxyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 173);

    auto g_xxxyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 174);

    auto g_xxxyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 175);

    auto g_xxxyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 176);

    auto g_xxxyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 177);

    auto g_xxxyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 178);

    auto g_xxxyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 179);

    #pragma omp simd aligned(g_xxxyyy_0_xxx_1, g_xxxyyy_0_xxxx_1, g_xxxyyy_0_xxxy_1, g_xxxyyy_0_xxxz_1, g_xxxyyy_0_xxy_1, g_xxxyyy_0_xxyy_1, g_xxxyyy_0_xxyz_1, g_xxxyyy_0_xxz_1, g_xxxyyy_0_xxzz_1, g_xxxyyy_0_xyy_1, g_xxxyyy_0_xyyy_1, g_xxxyyy_0_xyyz_1, g_xxxyyy_0_xyz_1, g_xxxyyy_0_xyzz_1, g_xxxyyy_0_xzz_1, g_xxxyyy_0_xzzz_1, g_xxxyyy_0_yyy_1, g_xxxyyy_0_yyyy_1, g_xxxyyy_0_yyyz_1, g_xxxyyy_0_yyz_1, g_xxxyyy_0_yyzz_1, g_xxxyyy_0_yzz_1, g_xxxyyy_0_yzzz_1, g_xxxyyy_0_zzz_1, g_xxxyyy_0_zzzz_1, g_xxxyyyz_0_xxxx_0, g_xxxyyyz_0_xxxy_0, g_xxxyyyz_0_xxxz_0, g_xxxyyyz_0_xxyy_0, g_xxxyyyz_0_xxyz_0, g_xxxyyyz_0_xxzz_0, g_xxxyyyz_0_xyyy_0, g_xxxyyyz_0_xyyz_0, g_xxxyyyz_0_xyzz_0, g_xxxyyyz_0_xzzz_0, g_xxxyyyz_0_yyyy_0, g_xxxyyyz_0_yyyz_0, g_xxxyyyz_0_yyzz_0, g_xxxyyyz_0_yzzz_0, g_xxxyyyz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyyz_0_xxxx_0[i] = g_xxxyyy_0_xxxx_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxy_0[i] = g_xxxyyy_0_xxxy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxz_0[i] = g_xxxyyy_0_xxx_1[i] * fi_acd_0 + g_xxxyyy_0_xxxz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyy_0[i] = g_xxxyyy_0_xxyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyz_0[i] = g_xxxyyy_0_xxy_1[i] * fi_acd_0 + g_xxxyyy_0_xxyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxzz_0[i] = 2.0 * g_xxxyyy_0_xxz_1[i] * fi_acd_0 + g_xxxyyy_0_xxzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyy_0[i] = g_xxxyyy_0_xyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyz_0[i] = g_xxxyyy_0_xyy_1[i] * fi_acd_0 + g_xxxyyy_0_xyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyzz_0[i] = 2.0 * g_xxxyyy_0_xyz_1[i] * fi_acd_0 + g_xxxyyy_0_xyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xzzz_0[i] = 3.0 * g_xxxyyy_0_xzz_1[i] * fi_acd_0 + g_xxxyyy_0_xzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyy_0[i] = g_xxxyyy_0_yyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyz_0[i] = g_xxxyyy_0_yyy_1[i] * fi_acd_0 + g_xxxyyy_0_yyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyzz_0[i] = 2.0 * g_xxxyyy_0_yyz_1[i] * fi_acd_0 + g_xxxyyy_0_yyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yzzz_0[i] = 3.0 * g_xxxyyy_0_yzz_1[i] * fi_acd_0 + g_xxxyyy_0_yzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_zzzz_0[i] = 4.0 * g_xxxyyy_0_zzz_1[i] * fi_acd_0 + g_xxxyyy_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 180-195 components of targeted buffer : KSG

    auto g_xxxyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 180);

    auto g_xxxyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 181);

    auto g_xxxyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 182);

    auto g_xxxyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 183);

    auto g_xxxyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 184);

    auto g_xxxyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 185);

    auto g_xxxyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 186);

    auto g_xxxyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 187);

    auto g_xxxyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 188);

    auto g_xxxyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 189);

    auto g_xxxyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 190);

    auto g_xxxyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 191);

    auto g_xxxyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 192);

    auto g_xxxyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 193);

    auto g_xxxyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 194);

    #pragma omp simd aligned(g_xxxyy_0_xxxy_0, g_xxxyy_0_xxxy_1, g_xxxyy_0_xxyy_0, g_xxxyy_0_xxyy_1, g_xxxyy_0_xyyy_0, g_xxxyy_0_xyyy_1, g_xxxyyz_0_xxxy_1, g_xxxyyz_0_xxyy_1, g_xxxyyz_0_xyyy_1, g_xxxyyzz_0_xxxx_0, g_xxxyyzz_0_xxxy_0, g_xxxyyzz_0_xxxz_0, g_xxxyyzz_0_xxyy_0, g_xxxyyzz_0_xxyz_0, g_xxxyyzz_0_xxzz_0, g_xxxyyzz_0_xyyy_0, g_xxxyyzz_0_xyyz_0, g_xxxyyzz_0_xyzz_0, g_xxxyyzz_0_xzzz_0, g_xxxyyzz_0_yyyy_0, g_xxxyyzz_0_yyyz_0, g_xxxyyzz_0_yyzz_0, g_xxxyyzz_0_yzzz_0, g_xxxyyzz_0_zzzz_0, g_xxxyzz_0_xxxx_1, g_xxxyzz_0_xxxz_1, g_xxxyzz_0_xxzz_1, g_xxxyzz_0_xzzz_1, g_xxxzz_0_xxxx_0, g_xxxzz_0_xxxx_1, g_xxxzz_0_xxxz_0, g_xxxzz_0_xxxz_1, g_xxxzz_0_xxzz_0, g_xxxzz_0_xxzz_1, g_xxxzz_0_xzzz_0, g_xxxzz_0_xzzz_1, g_xxyyzz_0_xxyz_1, g_xxyyzz_0_xyyz_1, g_xxyyzz_0_xyz_1, g_xxyyzz_0_xyzz_1, g_xxyyzz_0_yyyy_1, g_xxyyzz_0_yyyz_1, g_xxyyzz_0_yyz_1, g_xxyyzz_0_yyzz_1, g_xxyyzz_0_yzz_1, g_xxyyzz_0_yzzz_1, g_xxyyzz_0_zzzz_1, g_xyyzz_0_xxyz_0, g_xyyzz_0_xxyz_1, g_xyyzz_0_xyyz_0, g_xyyzz_0_xyyz_1, g_xyyzz_0_xyzz_0, g_xyyzz_0_xyzz_1, g_xyyzz_0_yyyy_0, g_xyyzz_0_yyyy_1, g_xyyzz_0_yyyz_0, g_xyyzz_0_yyyz_1, g_xyyzz_0_yyzz_0, g_xyyzz_0_yyzz_1, g_xyyzz_0_yzzz_0, g_xyyzz_0_yzzz_1, g_xyyzz_0_zzzz_0, g_xyyzz_0_zzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyzz_0_xxxx_0[i] = g_xxxzz_0_xxxx_0[i] * fbe_0 - g_xxxzz_0_xxxx_1[i] * fz_be_0 + g_xxxyzz_0_xxxx_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxy_0[i] = g_xxxyy_0_xxxy_0[i] * fbe_0 - g_xxxyy_0_xxxy_1[i] * fz_be_0 + g_xxxyyz_0_xxxy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxz_0[i] = g_xxxzz_0_xxxz_0[i] * fbe_0 - g_xxxzz_0_xxxz_1[i] * fz_be_0 + g_xxxyzz_0_xxxz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxyy_0[i] = g_xxxyy_0_xxyy_0[i] * fbe_0 - g_xxxyy_0_xxyy_1[i] * fz_be_0 + g_xxxyyz_0_xxyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxyz_0[i] = 2.0 * g_xyyzz_0_xxyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxzz_0[i] = g_xxxzz_0_xxzz_0[i] * fbe_0 - g_xxxzz_0_xxzz_1[i] * fz_be_0 + g_xxxyzz_0_xxzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xyyy_0[i] = g_xxxyy_0_xyyy_0[i] * fbe_0 - g_xxxyy_0_xyyy_1[i] * fz_be_0 + g_xxxyyz_0_xyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xyyz_0[i] = 2.0 * g_xyyzz_0_xyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyz_1[i] * fz_be_0 + g_xxyyzz_0_yyz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyzz_0[i] = 2.0 * g_xyyzz_0_xyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyzz_1[i] * fz_be_0 + g_xxyyzz_0_yzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xzzz_0[i] = g_xxxzz_0_xzzz_0[i] * fbe_0 - g_xxxzz_0_xzzz_1[i] * fz_be_0 + g_xxxyzz_0_xzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_yyyy_0[i] = 2.0 * g_xyyzz_0_yyyy_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyy_1[i] * fz_be_0 + g_xxyyzz_0_yyyy_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyz_0[i] = 2.0 * g_xyyzz_0_yyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyz_1[i] * fz_be_0 + g_xxyyzz_0_yyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyzz_0[i] = 2.0 * g_xyyzz_0_yyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyzz_1[i] * fz_be_0 + g_xxyyzz_0_yyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yzzz_0[i] = 2.0 * g_xyyzz_0_yzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yzzz_1[i] * fz_be_0 + g_xxyyzz_0_yzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_zzzz_0[i] = 2.0 * g_xyyzz_0_zzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_zzzz_1[i] * fz_be_0 + g_xxyyzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 195-210 components of targeted buffer : KSG

    auto g_xxxyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 195);

    auto g_xxxyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 196);

    auto g_xxxyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 197);

    auto g_xxxyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 198);

    auto g_xxxyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 199);

    auto g_xxxyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 200);

    auto g_xxxyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 201);

    auto g_xxxyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 202);

    auto g_xxxyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 203);

    auto g_xxxyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 204);

    auto g_xxxyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 205);

    auto g_xxxyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 206);

    auto g_xxxyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 207);

    auto g_xxxyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 208);

    auto g_xxxyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 209);

    #pragma omp simd aligned(g_xxxyzzz_0_xxxx_0, g_xxxyzzz_0_xxxy_0, g_xxxyzzz_0_xxxz_0, g_xxxyzzz_0_xxyy_0, g_xxxyzzz_0_xxyz_0, g_xxxyzzz_0_xxzz_0, g_xxxyzzz_0_xyyy_0, g_xxxyzzz_0_xyyz_0, g_xxxyzzz_0_xyzz_0, g_xxxyzzz_0_xzzz_0, g_xxxyzzz_0_yyyy_0, g_xxxyzzz_0_yyyz_0, g_xxxyzzz_0_yyzz_0, g_xxxyzzz_0_yzzz_0, g_xxxyzzz_0_zzzz_0, g_xxxzzz_0_xxx_1, g_xxxzzz_0_xxxx_1, g_xxxzzz_0_xxxy_1, g_xxxzzz_0_xxxz_1, g_xxxzzz_0_xxy_1, g_xxxzzz_0_xxyy_1, g_xxxzzz_0_xxyz_1, g_xxxzzz_0_xxz_1, g_xxxzzz_0_xxzz_1, g_xxxzzz_0_xyy_1, g_xxxzzz_0_xyyy_1, g_xxxzzz_0_xyyz_1, g_xxxzzz_0_xyz_1, g_xxxzzz_0_xyzz_1, g_xxxzzz_0_xzz_1, g_xxxzzz_0_xzzz_1, g_xxxzzz_0_yyy_1, g_xxxzzz_0_yyyy_1, g_xxxzzz_0_yyyz_1, g_xxxzzz_0_yyz_1, g_xxxzzz_0_yyzz_1, g_xxxzzz_0_yzz_1, g_xxxzzz_0_yzzz_1, g_xxxzzz_0_zzz_1, g_xxxzzz_0_zzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzzz_0_xxxx_0[i] = g_xxxzzz_0_xxxx_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxy_0[i] = g_xxxzzz_0_xxx_1[i] * fi_acd_0 + g_xxxzzz_0_xxxy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxz_0[i] = g_xxxzzz_0_xxxz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyy_0[i] = 2.0 * g_xxxzzz_0_xxy_1[i] * fi_acd_0 + g_xxxzzz_0_xxyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyz_0[i] = g_xxxzzz_0_xxz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxzz_0[i] = g_xxxzzz_0_xxzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyy_0[i] = 3.0 * g_xxxzzz_0_xyy_1[i] * fi_acd_0 + g_xxxzzz_0_xyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyz_0[i] = 2.0 * g_xxxzzz_0_xyz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyzz_0[i] = g_xxxzzz_0_xzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xzzz_0[i] = g_xxxzzz_0_xzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyy_0[i] = 4.0 * g_xxxzzz_0_yyy_1[i] * fi_acd_0 + g_xxxzzz_0_yyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyz_0[i] = 3.0 * g_xxxzzz_0_yyz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyzz_0[i] = 2.0 * g_xxxzzz_0_yzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yzzz_0[i] = g_xxxzzz_0_zzz_1[i] * fi_acd_0 + g_xxxzzz_0_yzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_zzzz_0[i] = g_xxxzzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 210-225 components of targeted buffer : KSG

    auto g_xxxzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 210);

    auto g_xxxzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 211);

    auto g_xxxzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 212);

    auto g_xxxzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 213);

    auto g_xxxzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 214);

    auto g_xxxzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 215);

    auto g_xxxzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 216);

    auto g_xxxzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 217);

    auto g_xxxzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 218);

    auto g_xxxzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 219);

    auto g_xxxzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 220);

    auto g_xxxzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 221);

    auto g_xxxzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 222);

    auto g_xxxzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 223);

    auto g_xxxzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 224);

    #pragma omp simd aligned(g_xxxzz_0_xxxx_0, g_xxxzz_0_xxxx_1, g_xxxzz_0_xxxy_0, g_xxxzz_0_xxxy_1, g_xxxzz_0_xxyy_0, g_xxxzz_0_xxyy_1, g_xxxzz_0_xyyy_0, g_xxxzz_0_xyyy_1, g_xxxzzz_0_xxxx_1, g_xxxzzz_0_xxxy_1, g_xxxzzz_0_xxyy_1, g_xxxzzz_0_xyyy_1, g_xxxzzzz_0_xxxx_0, g_xxxzzzz_0_xxxy_0, g_xxxzzzz_0_xxxz_0, g_xxxzzzz_0_xxyy_0, g_xxxzzzz_0_xxyz_0, g_xxxzzzz_0_xxzz_0, g_xxxzzzz_0_xyyy_0, g_xxxzzzz_0_xyyz_0, g_xxxzzzz_0_xyzz_0, g_xxxzzzz_0_xzzz_0, g_xxxzzzz_0_yyyy_0, g_xxxzzzz_0_yyyz_0, g_xxxzzzz_0_yyzz_0, g_xxxzzzz_0_yzzz_0, g_xxxzzzz_0_zzzz_0, g_xxzzzz_0_xxxz_1, g_xxzzzz_0_xxyz_1, g_xxzzzz_0_xxz_1, g_xxzzzz_0_xxzz_1, g_xxzzzz_0_xyyz_1, g_xxzzzz_0_xyz_1, g_xxzzzz_0_xyzz_1, g_xxzzzz_0_xzz_1, g_xxzzzz_0_xzzz_1, g_xxzzzz_0_yyyy_1, g_xxzzzz_0_yyyz_1, g_xxzzzz_0_yyz_1, g_xxzzzz_0_yyzz_1, g_xxzzzz_0_yzz_1, g_xxzzzz_0_yzzz_1, g_xxzzzz_0_zzz_1, g_xxzzzz_0_zzzz_1, g_xzzzz_0_xxxz_0, g_xzzzz_0_xxxz_1, g_xzzzz_0_xxyz_0, g_xzzzz_0_xxyz_1, g_xzzzz_0_xxzz_0, g_xzzzz_0_xxzz_1, g_xzzzz_0_xyyz_0, g_xzzzz_0_xyyz_1, g_xzzzz_0_xyzz_0, g_xzzzz_0_xyzz_1, g_xzzzz_0_xzzz_0, g_xzzzz_0_xzzz_1, g_xzzzz_0_yyyy_0, g_xzzzz_0_yyyy_1, g_xzzzz_0_yyyz_0, g_xzzzz_0_yyyz_1, g_xzzzz_0_yyzz_0, g_xzzzz_0_yyzz_1, g_xzzzz_0_yzzz_0, g_xzzzz_0_yzzz_1, g_xzzzz_0_zzzz_0, g_xzzzz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzzz_0_xxxx_0[i] = 3.0 * g_xxxzz_0_xxxx_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxx_1[i] * fz_be_0 + g_xxxzzz_0_xxxx_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxy_0[i] = 3.0 * g_xxxzz_0_xxxy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxy_1[i] * fz_be_0 + g_xxxzzz_0_xxxy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxz_0[i] = 2.0 * g_xzzzz_0_xxxz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyy_0[i] = 3.0 * g_xxxzz_0_xxyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxyy_1[i] * fz_be_0 + g_xxxzzz_0_xxyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxyz_0[i] = 2.0 * g_xzzzz_0_xxyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxzz_0[i] = 2.0 * g_xzzzz_0_xxzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyy_0[i] = 3.0 * g_xxxzz_0_xyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xyyy_1[i] * fz_be_0 + g_xxxzzz_0_xyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xyyz_0[i] = 2.0 * g_xzzzz_0_xyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyz_1[i] * fz_be_0 + g_xxzzzz_0_yyz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyzz_0[i] = 2.0 * g_xzzzz_0_xyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyzz_1[i] * fz_be_0 + g_xxzzzz_0_yzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xzzz_0[i] = 2.0 * g_xzzzz_0_xzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xzzz_1[i] * fz_be_0 + g_xxzzzz_0_zzz_1[i] * fi_acd_0 + g_xxzzzz_0_xzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyy_0[i] = 2.0 * g_xzzzz_0_yyyy_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyy_1[i] * fz_be_0 + g_xxzzzz_0_yyyy_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyz_0[i] = 2.0 * g_xzzzz_0_yyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyz_1[i] * fz_be_0 + g_xxzzzz_0_yyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyzz_0[i] = 2.0 * g_xzzzz_0_yyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyzz_1[i] * fz_be_0 + g_xxzzzz_0_yyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yzzz_0[i] = 2.0 * g_xzzzz_0_yzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yzzz_1[i] * fz_be_0 + g_xxzzzz_0_yzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_zzzz_0[i] = 2.0 * g_xzzzz_0_zzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_zzzz_1[i] * fz_be_0 + g_xxzzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 225-240 components of targeted buffer : KSG

    auto g_xxyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 225);

    auto g_xxyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 226);

    auto g_xxyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 227);

    auto g_xxyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 228);

    auto g_xxyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 229);

    auto g_xxyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 230);

    auto g_xxyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 231);

    auto g_xxyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 232);

    auto g_xxyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 233);

    auto g_xxyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 234);

    auto g_xxyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 235);

    auto g_xxyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 236);

    auto g_xxyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 237);

    auto g_xxyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 238);

    auto g_xxyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 239);

    #pragma omp simd aligned(g_xxyyy_0_xxxx_0, g_xxyyy_0_xxxx_1, g_xxyyy_0_xxxz_0, g_xxyyy_0_xxxz_1, g_xxyyy_0_xxzz_0, g_xxyyy_0_xxzz_1, g_xxyyy_0_xzzz_0, g_xxyyy_0_xzzz_1, g_xxyyyy_0_xxxx_1, g_xxyyyy_0_xxxz_1, g_xxyyyy_0_xxzz_1, g_xxyyyy_0_xzzz_1, g_xxyyyyy_0_xxxx_0, g_xxyyyyy_0_xxxy_0, g_xxyyyyy_0_xxxz_0, g_xxyyyyy_0_xxyy_0, g_xxyyyyy_0_xxyz_0, g_xxyyyyy_0_xxzz_0, g_xxyyyyy_0_xyyy_0, g_xxyyyyy_0_xyyz_0, g_xxyyyyy_0_xyzz_0, g_xxyyyyy_0_xzzz_0, g_xxyyyyy_0_yyyy_0, g_xxyyyyy_0_yyyz_0, g_xxyyyyy_0_yyzz_0, g_xxyyyyy_0_yzzz_0, g_xxyyyyy_0_zzzz_0, g_xyyyyy_0_xxxy_1, g_xyyyyy_0_xxy_1, g_xyyyyy_0_xxyy_1, g_xyyyyy_0_xxyz_1, g_xyyyyy_0_xyy_1, g_xyyyyy_0_xyyy_1, g_xyyyyy_0_xyyz_1, g_xyyyyy_0_xyz_1, g_xyyyyy_0_xyzz_1, g_xyyyyy_0_yyy_1, g_xyyyyy_0_yyyy_1, g_xyyyyy_0_yyyz_1, g_xyyyyy_0_yyz_1, g_xyyyyy_0_yyzz_1, g_xyyyyy_0_yzz_1, g_xyyyyy_0_yzzz_1, g_xyyyyy_0_zzzz_1, g_yyyyy_0_xxxy_0, g_yyyyy_0_xxxy_1, g_yyyyy_0_xxyy_0, g_yyyyy_0_xxyy_1, g_yyyyy_0_xxyz_0, g_yyyyy_0_xxyz_1, g_yyyyy_0_xyyy_0, g_yyyyy_0_xyyy_1, g_yyyyy_0_xyyz_0, g_yyyyy_0_xyyz_1, g_yyyyy_0_xyzz_0, g_yyyyy_0_xyzz_1, g_yyyyy_0_yyyy_0, g_yyyyy_0_yyyy_1, g_yyyyy_0_yyyz_0, g_yyyyy_0_yyyz_1, g_yyyyy_0_yyzz_0, g_yyyyy_0_yyzz_1, g_yyyyy_0_yzzz_0, g_yyyyy_0_yzzz_1, g_yyyyy_0_zzzz_0, g_yyyyy_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyyy_0_xxxx_0[i] = 4.0 * g_xxyyy_0_xxxx_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxx_1[i] * fz_be_0 + g_xxyyyy_0_xxxx_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxy_0[i] = g_yyyyy_0_xxxy_0[i] * fbe_0 - g_yyyyy_0_xxxy_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxz_0[i] = 4.0 * g_xxyyy_0_xxxz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxz_1[i] * fz_be_0 + g_xxyyyy_0_xxxz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxyy_0[i] = g_yyyyy_0_xxyy_0[i] * fbe_0 - g_yyyyy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyz_0[i] = g_yyyyy_0_xxyz_0[i] * fbe_0 - g_yyyyy_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxzz_0[i] = 4.0 * g_xxyyy_0_xxzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxzz_1[i] * fz_be_0 + g_xxyyyy_0_xxzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xyyy_0[i] = g_yyyyy_0_xyyy_0[i] * fbe_0 - g_yyyyy_0_xyyy_1[i] * fz_be_0 + g_xyyyyy_0_yyy_1[i] * fi_acd_0 + g_xyyyyy_0_xyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyz_0[i] = g_yyyyy_0_xyyz_0[i] * fbe_0 - g_yyyyy_0_xyyz_1[i] * fz_be_0 + g_xyyyyy_0_yyz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyzz_0[i] = g_yyyyy_0_xyzz_0[i] * fbe_0 - g_yyyyy_0_xyzz_1[i] * fz_be_0 + g_xyyyyy_0_yzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xzzz_0[i] = 4.0 * g_xxyyy_0_xzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xzzz_1[i] * fz_be_0 + g_xxyyyy_0_xzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_yyyy_0[i] = g_yyyyy_0_yyyy_0[i] * fbe_0 - g_yyyyy_0_yyyy_1[i] * fz_be_0 + g_xyyyyy_0_yyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyz_0[i] = g_yyyyy_0_yyyz_0[i] * fbe_0 - g_yyyyy_0_yyyz_1[i] * fz_be_0 + g_xyyyyy_0_yyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyzz_0[i] = g_yyyyy_0_yyzz_0[i] * fbe_0 - g_yyyyy_0_yyzz_1[i] * fz_be_0 + g_xyyyyy_0_yyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yzzz_0[i] = g_yyyyy_0_yzzz_0[i] * fbe_0 - g_yyyyy_0_yzzz_1[i] * fz_be_0 + g_xyyyyy_0_yzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_zzzz_0[i] = g_yyyyy_0_zzzz_0[i] * fbe_0 - g_yyyyy_0_zzzz_1[i] * fz_be_0 + g_xyyyyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 240-255 components of targeted buffer : KSG

    auto g_xxyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 240);

    auto g_xxyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 241);

    auto g_xxyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 242);

    auto g_xxyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 243);

    auto g_xxyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 244);

    auto g_xxyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 245);

    auto g_xxyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 246);

    auto g_xxyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 247);

    auto g_xxyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 248);

    auto g_xxyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 249);

    auto g_xxyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 250);

    auto g_xxyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 251);

    auto g_xxyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 252);

    auto g_xxyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 253);

    auto g_xxyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 254);

    #pragma omp simd aligned(g_xxyyyy_0_xxx_1, g_xxyyyy_0_xxxx_1, g_xxyyyy_0_xxxy_1, g_xxyyyy_0_xxxz_1, g_xxyyyy_0_xxy_1, g_xxyyyy_0_xxyy_1, g_xxyyyy_0_xxyz_1, g_xxyyyy_0_xxz_1, g_xxyyyy_0_xxzz_1, g_xxyyyy_0_xyy_1, g_xxyyyy_0_xyyy_1, g_xxyyyy_0_xyyz_1, g_xxyyyy_0_xyz_1, g_xxyyyy_0_xyzz_1, g_xxyyyy_0_xzz_1, g_xxyyyy_0_xzzz_1, g_xxyyyy_0_yyy_1, g_xxyyyy_0_yyyy_1, g_xxyyyy_0_yyyz_1, g_xxyyyy_0_yyz_1, g_xxyyyy_0_yyzz_1, g_xxyyyy_0_yzz_1, g_xxyyyy_0_yzzz_1, g_xxyyyy_0_zzz_1, g_xxyyyy_0_zzzz_1, g_xxyyyyz_0_xxxx_0, g_xxyyyyz_0_xxxy_0, g_xxyyyyz_0_xxxz_0, g_xxyyyyz_0_xxyy_0, g_xxyyyyz_0_xxyz_0, g_xxyyyyz_0_xxzz_0, g_xxyyyyz_0_xyyy_0, g_xxyyyyz_0_xyyz_0, g_xxyyyyz_0_xyzz_0, g_xxyyyyz_0_xzzz_0, g_xxyyyyz_0_yyyy_0, g_xxyyyyz_0_yyyz_0, g_xxyyyyz_0_yyzz_0, g_xxyyyyz_0_yzzz_0, g_xxyyyyz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyyz_0_xxxx_0[i] = g_xxyyyy_0_xxxx_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxy_0[i] = g_xxyyyy_0_xxxy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxz_0[i] = g_xxyyyy_0_xxx_1[i] * fi_acd_0 + g_xxyyyy_0_xxxz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyy_0[i] = g_xxyyyy_0_xxyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyz_0[i] = g_xxyyyy_0_xxy_1[i] * fi_acd_0 + g_xxyyyy_0_xxyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxzz_0[i] = 2.0 * g_xxyyyy_0_xxz_1[i] * fi_acd_0 + g_xxyyyy_0_xxzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyy_0[i] = g_xxyyyy_0_xyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyz_0[i] = g_xxyyyy_0_xyy_1[i] * fi_acd_0 + g_xxyyyy_0_xyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyzz_0[i] = 2.0 * g_xxyyyy_0_xyz_1[i] * fi_acd_0 + g_xxyyyy_0_xyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xzzz_0[i] = 3.0 * g_xxyyyy_0_xzz_1[i] * fi_acd_0 + g_xxyyyy_0_xzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyy_0[i] = g_xxyyyy_0_yyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyz_0[i] = g_xxyyyy_0_yyy_1[i] * fi_acd_0 + g_xxyyyy_0_yyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyzz_0[i] = 2.0 * g_xxyyyy_0_yyz_1[i] * fi_acd_0 + g_xxyyyy_0_yyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yzzz_0[i] = 3.0 * g_xxyyyy_0_yzz_1[i] * fi_acd_0 + g_xxyyyy_0_yzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_zzzz_0[i] = 4.0 * g_xxyyyy_0_zzz_1[i] * fi_acd_0 + g_xxyyyy_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 255-270 components of targeted buffer : KSG

    auto g_xxyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 255);

    auto g_xxyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 256);

    auto g_xxyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 257);

    auto g_xxyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 258);

    auto g_xxyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 259);

    auto g_xxyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 260);

    auto g_xxyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 261);

    auto g_xxyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 262);

    auto g_xxyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 263);

    auto g_xxyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 264);

    auto g_xxyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 265);

    auto g_xxyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 266);

    auto g_xxyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 267);

    auto g_xxyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 268);

    auto g_xxyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 269);

    #pragma omp simd aligned(g_xxyyy_0_xxxy_0, g_xxyyy_0_xxxy_1, g_xxyyy_0_xxyy_0, g_xxyyy_0_xxyy_1, g_xxyyy_0_xyyy_0, g_xxyyy_0_xyyy_1, g_xxyyyz_0_xxxy_1, g_xxyyyz_0_xxyy_1, g_xxyyyz_0_xyyy_1, g_xxyyyzz_0_xxxx_0, g_xxyyyzz_0_xxxy_0, g_xxyyyzz_0_xxxz_0, g_xxyyyzz_0_xxyy_0, g_xxyyyzz_0_xxyz_0, g_xxyyyzz_0_xxzz_0, g_xxyyyzz_0_xyyy_0, g_xxyyyzz_0_xyyz_0, g_xxyyyzz_0_xyzz_0, g_xxyyyzz_0_xzzz_0, g_xxyyyzz_0_yyyy_0, g_xxyyyzz_0_yyyz_0, g_xxyyyzz_0_yyzz_0, g_xxyyyzz_0_yzzz_0, g_xxyyyzz_0_zzzz_0, g_xxyyzz_0_xxxx_1, g_xxyyzz_0_xxxz_1, g_xxyyzz_0_xxzz_1, g_xxyyzz_0_xzzz_1, g_xxyzz_0_xxxx_0, g_xxyzz_0_xxxx_1, g_xxyzz_0_xxxz_0, g_xxyzz_0_xxxz_1, g_xxyzz_0_xxzz_0, g_xxyzz_0_xxzz_1, g_xxyzz_0_xzzz_0, g_xxyzz_0_xzzz_1, g_xyyyzz_0_xxyz_1, g_xyyyzz_0_xyyz_1, g_xyyyzz_0_xyz_1, g_xyyyzz_0_xyzz_1, g_xyyyzz_0_yyyy_1, g_xyyyzz_0_yyyz_1, g_xyyyzz_0_yyz_1, g_xyyyzz_0_yyzz_1, g_xyyyzz_0_yzz_1, g_xyyyzz_0_yzzz_1, g_xyyyzz_0_zzzz_1, g_yyyzz_0_xxyz_0, g_yyyzz_0_xxyz_1, g_yyyzz_0_xyyz_0, g_yyyzz_0_xyyz_1, g_yyyzz_0_xyzz_0, g_yyyzz_0_xyzz_1, g_yyyzz_0_yyyy_0, g_yyyzz_0_yyyy_1, g_yyyzz_0_yyyz_0, g_yyyzz_0_yyyz_1, g_yyyzz_0_yyzz_0, g_yyyzz_0_yyzz_1, g_yyyzz_0_yzzz_0, g_yyyzz_0_yzzz_1, g_yyyzz_0_zzzz_0, g_yyyzz_0_zzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyzz_0_xxxx_0[i] = 2.0 * g_xxyzz_0_xxxx_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxx_1[i] * fz_be_0 + g_xxyyzz_0_xxxx_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxy_0[i] = g_xxyyy_0_xxxy_0[i] * fbe_0 - g_xxyyy_0_xxxy_1[i] * fz_be_0 + g_xxyyyz_0_xxxy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxz_0[i] = 2.0 * g_xxyzz_0_xxxz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxz_1[i] * fz_be_0 + g_xxyyzz_0_xxxz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxyy_0[i] = g_xxyyy_0_xxyy_0[i] * fbe_0 - g_xxyyy_0_xxyy_1[i] * fz_be_0 + g_xxyyyz_0_xxyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxyz_0[i] = g_yyyzz_0_xxyz_0[i] * fbe_0 - g_yyyzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxzz_0[i] = 2.0 * g_xxyzz_0_xxzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxzz_1[i] * fz_be_0 + g_xxyyzz_0_xxzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xyyy_0[i] = g_xxyyy_0_xyyy_0[i] * fbe_0 - g_xxyyy_0_xyyy_1[i] * fz_be_0 + g_xxyyyz_0_xyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xyyz_0[i] = g_yyyzz_0_xyyz_0[i] * fbe_0 - g_yyyzz_0_xyyz_1[i] * fz_be_0 + g_xyyyzz_0_yyz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyzz_0[i] = g_yyyzz_0_xyzz_0[i] * fbe_0 - g_yyyzz_0_xyzz_1[i] * fz_be_0 + g_xyyyzz_0_yzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xzzz_0[i] = 2.0 * g_xxyzz_0_xzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xzzz_1[i] * fz_be_0 + g_xxyyzz_0_xzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_yyyy_0[i] = g_yyyzz_0_yyyy_0[i] * fbe_0 - g_yyyzz_0_yyyy_1[i] * fz_be_0 + g_xyyyzz_0_yyyy_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyz_0[i] = g_yyyzz_0_yyyz_0[i] * fbe_0 - g_yyyzz_0_yyyz_1[i] * fz_be_0 + g_xyyyzz_0_yyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyzz_0[i] = g_yyyzz_0_yyzz_0[i] * fbe_0 - g_yyyzz_0_yyzz_1[i] * fz_be_0 + g_xyyyzz_0_yyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yzzz_0[i] = g_yyyzz_0_yzzz_0[i] * fbe_0 - g_yyyzz_0_yzzz_1[i] * fz_be_0 + g_xyyyzz_0_yzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_zzzz_0[i] = g_yyyzz_0_zzzz_0[i] * fbe_0 - g_yyyzz_0_zzzz_1[i] * fz_be_0 + g_xyyyzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 270-285 components of targeted buffer : KSG

    auto g_xxyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 270);

    auto g_xxyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 271);

    auto g_xxyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 272);

    auto g_xxyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 273);

    auto g_xxyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 274);

    auto g_xxyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 275);

    auto g_xxyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 276);

    auto g_xxyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 277);

    auto g_xxyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 278);

    auto g_xxyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 279);

    auto g_xxyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 280);

    auto g_xxyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 281);

    auto g_xxyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 282);

    auto g_xxyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 283);

    auto g_xxyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 284);

    #pragma omp simd aligned(g_xxyyz_0_xxxy_0, g_xxyyz_0_xxxy_1, g_xxyyz_0_xxyy_0, g_xxyyz_0_xxyy_1, g_xxyyz_0_xyyy_0, g_xxyyz_0_xyyy_1, g_xxyyzz_0_xxxy_1, g_xxyyzz_0_xxyy_1, g_xxyyzz_0_xyyy_1, g_xxyyzzz_0_xxxx_0, g_xxyyzzz_0_xxxy_0, g_xxyyzzz_0_xxxz_0, g_xxyyzzz_0_xxyy_0, g_xxyyzzz_0_xxyz_0, g_xxyyzzz_0_xxzz_0, g_xxyyzzz_0_xyyy_0, g_xxyyzzz_0_xyyz_0, g_xxyyzzz_0_xyzz_0, g_xxyyzzz_0_xzzz_0, g_xxyyzzz_0_yyyy_0, g_xxyyzzz_0_yyyz_0, g_xxyyzzz_0_yyzz_0, g_xxyyzzz_0_yzzz_0, g_xxyyzzz_0_zzzz_0, g_xxyzzz_0_xxxx_1, g_xxyzzz_0_xxxz_1, g_xxyzzz_0_xxzz_1, g_xxyzzz_0_xzzz_1, g_xxzzz_0_xxxx_0, g_xxzzz_0_xxxx_1, g_xxzzz_0_xxxz_0, g_xxzzz_0_xxxz_1, g_xxzzz_0_xxzz_0, g_xxzzz_0_xxzz_1, g_xxzzz_0_xzzz_0, g_xxzzz_0_xzzz_1, g_xyyzzz_0_xxyz_1, g_xyyzzz_0_xyyz_1, g_xyyzzz_0_xyz_1, g_xyyzzz_0_xyzz_1, g_xyyzzz_0_yyyy_1, g_xyyzzz_0_yyyz_1, g_xyyzzz_0_yyz_1, g_xyyzzz_0_yyzz_1, g_xyyzzz_0_yzz_1, g_xyyzzz_0_yzzz_1, g_xyyzzz_0_zzzz_1, g_yyzzz_0_xxyz_0, g_yyzzz_0_xxyz_1, g_yyzzz_0_xyyz_0, g_yyzzz_0_xyyz_1, g_yyzzz_0_xyzz_0, g_yyzzz_0_xyzz_1, g_yyzzz_0_yyyy_0, g_yyzzz_0_yyyy_1, g_yyzzz_0_yyyz_0, g_yyzzz_0_yyyz_1, g_yyzzz_0_yyzz_0, g_yyzzz_0_yyzz_1, g_yyzzz_0_yzzz_0, g_yyzzz_0_yzzz_1, g_yyzzz_0_zzzz_0, g_yyzzz_0_zzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyzzz_0_xxxx_0[i] = g_xxzzz_0_xxxx_0[i] * fbe_0 - g_xxzzz_0_xxxx_1[i] * fz_be_0 + g_xxyzzz_0_xxxx_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxy_0[i] = 2.0 * g_xxyyz_0_xxxy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxy_1[i] * fz_be_0 + g_xxyyzz_0_xxxy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxz_0[i] = g_xxzzz_0_xxxz_0[i] * fbe_0 - g_xxzzz_0_xxxz_1[i] * fz_be_0 + g_xxyzzz_0_xxxz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxyy_0[i] = 2.0 * g_xxyyz_0_xxyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxyy_1[i] * fz_be_0 + g_xxyyzz_0_xxyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxyz_0[i] = g_yyzzz_0_xxyz_0[i] * fbe_0 - g_yyzzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxzz_0[i] = g_xxzzz_0_xxzz_0[i] * fbe_0 - g_xxzzz_0_xxzz_1[i] * fz_be_0 + g_xxyzzz_0_xxzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xyyy_0[i] = 2.0 * g_xxyyz_0_xyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xyyy_1[i] * fz_be_0 + g_xxyyzz_0_xyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xyyz_0[i] = g_yyzzz_0_xyyz_0[i] * fbe_0 - g_yyzzz_0_xyyz_1[i] * fz_be_0 + g_xyyzzz_0_yyz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyzz_0[i] = g_yyzzz_0_xyzz_0[i] * fbe_0 - g_yyzzz_0_xyzz_1[i] * fz_be_0 + g_xyyzzz_0_yzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xzzz_0[i] = g_xxzzz_0_xzzz_0[i] * fbe_0 - g_xxzzz_0_xzzz_1[i] * fz_be_0 + g_xxyzzz_0_xzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_yyyy_0[i] = g_yyzzz_0_yyyy_0[i] * fbe_0 - g_yyzzz_0_yyyy_1[i] * fz_be_0 + g_xyyzzz_0_yyyy_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyz_0[i] = g_yyzzz_0_yyyz_0[i] * fbe_0 - g_yyzzz_0_yyyz_1[i] * fz_be_0 + g_xyyzzz_0_yyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyzz_0[i] = g_yyzzz_0_yyzz_0[i] * fbe_0 - g_yyzzz_0_yyzz_1[i] * fz_be_0 + g_xyyzzz_0_yyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yzzz_0[i] = g_yyzzz_0_yzzz_0[i] * fbe_0 - g_yyzzz_0_yzzz_1[i] * fz_be_0 + g_xyyzzz_0_yzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_zzzz_0[i] = g_yyzzz_0_zzzz_0[i] * fbe_0 - g_yyzzz_0_zzzz_1[i] * fz_be_0 + g_xyyzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 285-300 components of targeted buffer : KSG

    auto g_xxyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 285);

    auto g_xxyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 286);

    auto g_xxyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 287);

    auto g_xxyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 288);

    auto g_xxyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 289);

    auto g_xxyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 290);

    auto g_xxyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 291);

    auto g_xxyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 292);

    auto g_xxyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 293);

    auto g_xxyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 294);

    auto g_xxyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 295);

    auto g_xxyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 296);

    auto g_xxyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 297);

    auto g_xxyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 298);

    auto g_xxyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 299);

    #pragma omp simd aligned(g_xxyzzzz_0_xxxx_0, g_xxyzzzz_0_xxxy_0, g_xxyzzzz_0_xxxz_0, g_xxyzzzz_0_xxyy_0, g_xxyzzzz_0_xxyz_0, g_xxyzzzz_0_xxzz_0, g_xxyzzzz_0_xyyy_0, g_xxyzzzz_0_xyyz_0, g_xxyzzzz_0_xyzz_0, g_xxyzzzz_0_xzzz_0, g_xxyzzzz_0_yyyy_0, g_xxyzzzz_0_yyyz_0, g_xxyzzzz_0_yyzz_0, g_xxyzzzz_0_yzzz_0, g_xxyzzzz_0_zzzz_0, g_xxzzzz_0_xxx_1, g_xxzzzz_0_xxxx_1, g_xxzzzz_0_xxxy_1, g_xxzzzz_0_xxxz_1, g_xxzzzz_0_xxy_1, g_xxzzzz_0_xxyy_1, g_xxzzzz_0_xxyz_1, g_xxzzzz_0_xxz_1, g_xxzzzz_0_xxzz_1, g_xxzzzz_0_xyy_1, g_xxzzzz_0_xyyy_1, g_xxzzzz_0_xyyz_1, g_xxzzzz_0_xyz_1, g_xxzzzz_0_xyzz_1, g_xxzzzz_0_xzz_1, g_xxzzzz_0_xzzz_1, g_xxzzzz_0_yyy_1, g_xxzzzz_0_yyyy_1, g_xxzzzz_0_yyyz_1, g_xxzzzz_0_yyz_1, g_xxzzzz_0_yyzz_1, g_xxzzzz_0_yzz_1, g_xxzzzz_0_yzzz_1, g_xxzzzz_0_zzz_1, g_xxzzzz_0_zzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzzz_0_xxxx_0[i] = g_xxzzzz_0_xxxx_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxy_0[i] = g_xxzzzz_0_xxx_1[i] * fi_acd_0 + g_xxzzzz_0_xxxy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxz_0[i] = g_xxzzzz_0_xxxz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyy_0[i] = 2.0 * g_xxzzzz_0_xxy_1[i] * fi_acd_0 + g_xxzzzz_0_xxyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyz_0[i] = g_xxzzzz_0_xxz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxzz_0[i] = g_xxzzzz_0_xxzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyy_0[i] = 3.0 * g_xxzzzz_0_xyy_1[i] * fi_acd_0 + g_xxzzzz_0_xyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyz_0[i] = 2.0 * g_xxzzzz_0_xyz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyzz_0[i] = g_xxzzzz_0_xzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xzzz_0[i] = g_xxzzzz_0_xzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyy_0[i] = 4.0 * g_xxzzzz_0_yyy_1[i] * fi_acd_0 + g_xxzzzz_0_yyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyz_0[i] = 3.0 * g_xxzzzz_0_yyz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyzz_0[i] = 2.0 * g_xxzzzz_0_yzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yzzz_0[i] = g_xxzzzz_0_zzz_1[i] * fi_acd_0 + g_xxzzzz_0_yzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_zzzz_0[i] = g_xxzzzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 300-315 components of targeted buffer : KSG

    auto g_xxzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 300);

    auto g_xxzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 301);

    auto g_xxzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 302);

    auto g_xxzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 303);

    auto g_xxzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 304);

    auto g_xxzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 305);

    auto g_xxzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 306);

    auto g_xxzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 307);

    auto g_xxzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 308);

    auto g_xxzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 309);

    auto g_xxzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 310);

    auto g_xxzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 311);

    auto g_xxzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 312);

    auto g_xxzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 313);

    auto g_xxzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 314);

    #pragma omp simd aligned(g_xxzzz_0_xxxx_0, g_xxzzz_0_xxxx_1, g_xxzzz_0_xxxy_0, g_xxzzz_0_xxxy_1, g_xxzzz_0_xxyy_0, g_xxzzz_0_xxyy_1, g_xxzzz_0_xyyy_0, g_xxzzz_0_xyyy_1, g_xxzzzz_0_xxxx_1, g_xxzzzz_0_xxxy_1, g_xxzzzz_0_xxyy_1, g_xxzzzz_0_xyyy_1, g_xxzzzzz_0_xxxx_0, g_xxzzzzz_0_xxxy_0, g_xxzzzzz_0_xxxz_0, g_xxzzzzz_0_xxyy_0, g_xxzzzzz_0_xxyz_0, g_xxzzzzz_0_xxzz_0, g_xxzzzzz_0_xyyy_0, g_xxzzzzz_0_xyyz_0, g_xxzzzzz_0_xyzz_0, g_xxzzzzz_0_xzzz_0, g_xxzzzzz_0_yyyy_0, g_xxzzzzz_0_yyyz_0, g_xxzzzzz_0_yyzz_0, g_xxzzzzz_0_yzzz_0, g_xxzzzzz_0_zzzz_0, g_xzzzzz_0_xxxz_1, g_xzzzzz_0_xxyz_1, g_xzzzzz_0_xxz_1, g_xzzzzz_0_xxzz_1, g_xzzzzz_0_xyyz_1, g_xzzzzz_0_xyz_1, g_xzzzzz_0_xyzz_1, g_xzzzzz_0_xzz_1, g_xzzzzz_0_xzzz_1, g_xzzzzz_0_yyyy_1, g_xzzzzz_0_yyyz_1, g_xzzzzz_0_yyz_1, g_xzzzzz_0_yyzz_1, g_xzzzzz_0_yzz_1, g_xzzzzz_0_yzzz_1, g_xzzzzz_0_zzz_1, g_xzzzzz_0_zzzz_1, g_zzzzz_0_xxxz_0, g_zzzzz_0_xxxz_1, g_zzzzz_0_xxyz_0, g_zzzzz_0_xxyz_1, g_zzzzz_0_xxzz_0, g_zzzzz_0_xxzz_1, g_zzzzz_0_xyyz_0, g_zzzzz_0_xyyz_1, g_zzzzz_0_xyzz_0, g_zzzzz_0_xyzz_1, g_zzzzz_0_xzzz_0, g_zzzzz_0_xzzz_1, g_zzzzz_0_yyyy_0, g_zzzzz_0_yyyy_1, g_zzzzz_0_yyyz_0, g_zzzzz_0_yyyz_1, g_zzzzz_0_yyzz_0, g_zzzzz_0_yyzz_1, g_zzzzz_0_yzzz_0, g_zzzzz_0_yzzz_1, g_zzzzz_0_zzzz_0, g_zzzzz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzzz_0_xxxx_0[i] = 4.0 * g_xxzzz_0_xxxx_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxx_1[i] * fz_be_0 + g_xxzzzz_0_xxxx_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxy_0[i] = 4.0 * g_xxzzz_0_xxxy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxy_1[i] * fz_be_0 + g_xxzzzz_0_xxxy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxz_0[i] = g_zzzzz_0_xxxz_0[i] * fbe_0 - g_zzzzz_0_xxxz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyy_0[i] = 4.0 * g_xxzzz_0_xxyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxyy_1[i] * fz_be_0 + g_xxzzzz_0_xxyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxyz_0[i] = g_zzzzz_0_xxyz_0[i] * fbe_0 - g_zzzzz_0_xxyz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxzz_0[i] = g_zzzzz_0_xxzz_0[i] * fbe_0 - g_zzzzz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyy_0[i] = 4.0 * g_xxzzz_0_xyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xyyy_1[i] * fz_be_0 + g_xxzzzz_0_xyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xyyz_0[i] = g_zzzzz_0_xyyz_0[i] * fbe_0 - g_zzzzz_0_xyyz_1[i] * fz_be_0 + g_xzzzzz_0_yyz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyzz_0[i] = g_zzzzz_0_xyzz_0[i] * fbe_0 - g_zzzzz_0_xyzz_1[i] * fz_be_0 + g_xzzzzz_0_yzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xzzz_0[i] = g_zzzzz_0_xzzz_0[i] * fbe_0 - g_zzzzz_0_xzzz_1[i] * fz_be_0 + g_xzzzzz_0_zzz_1[i] * fi_acd_0 + g_xzzzzz_0_xzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyy_0[i] = g_zzzzz_0_yyyy_0[i] * fbe_0 - g_zzzzz_0_yyyy_1[i] * fz_be_0 + g_xzzzzz_0_yyyy_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyz_0[i] = g_zzzzz_0_yyyz_0[i] * fbe_0 - g_zzzzz_0_yyyz_1[i] * fz_be_0 + g_xzzzzz_0_yyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyzz_0[i] = g_zzzzz_0_yyzz_0[i] * fbe_0 - g_zzzzz_0_yyzz_1[i] * fz_be_0 + g_xzzzzz_0_yyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yzzz_0[i] = g_zzzzz_0_yzzz_0[i] * fbe_0 - g_zzzzz_0_yzzz_1[i] * fz_be_0 + g_xzzzzz_0_yzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_zzzz_0[i] = g_zzzzz_0_zzzz_0[i] * fbe_0 - g_zzzzz_0_zzzz_1[i] * fz_be_0 + g_xzzzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 315-330 components of targeted buffer : KSG

    auto g_xyyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 315);

    auto g_xyyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 316);

    auto g_xyyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 317);

    auto g_xyyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 318);

    auto g_xyyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 319);

    auto g_xyyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 320);

    auto g_xyyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 321);

    auto g_xyyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 322);

    auto g_xyyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 323);

    auto g_xyyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 324);

    auto g_xyyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 325);

    auto g_xyyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 326);

    auto g_xyyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 327);

    auto g_xyyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 328);

    auto g_xyyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 329);

    #pragma omp simd aligned(g_xyyyyyy_0_xxxx_0, g_xyyyyyy_0_xxxy_0, g_xyyyyyy_0_xxxz_0, g_xyyyyyy_0_xxyy_0, g_xyyyyyy_0_xxyz_0, g_xyyyyyy_0_xxzz_0, g_xyyyyyy_0_xyyy_0, g_xyyyyyy_0_xyyz_0, g_xyyyyyy_0_xyzz_0, g_xyyyyyy_0_xzzz_0, g_xyyyyyy_0_yyyy_0, g_xyyyyyy_0_yyyz_0, g_xyyyyyy_0_yyzz_0, g_xyyyyyy_0_yzzz_0, g_xyyyyyy_0_zzzz_0, g_yyyyyy_0_xxx_1, g_yyyyyy_0_xxxx_1, g_yyyyyy_0_xxxy_1, g_yyyyyy_0_xxxz_1, g_yyyyyy_0_xxy_1, g_yyyyyy_0_xxyy_1, g_yyyyyy_0_xxyz_1, g_yyyyyy_0_xxz_1, g_yyyyyy_0_xxzz_1, g_yyyyyy_0_xyy_1, g_yyyyyy_0_xyyy_1, g_yyyyyy_0_xyyz_1, g_yyyyyy_0_xyz_1, g_yyyyyy_0_xyzz_1, g_yyyyyy_0_xzz_1, g_yyyyyy_0_xzzz_1, g_yyyyyy_0_yyy_1, g_yyyyyy_0_yyyy_1, g_yyyyyy_0_yyyz_1, g_yyyyyy_0_yyz_1, g_yyyyyy_0_yyzz_1, g_yyyyyy_0_yzz_1, g_yyyyyy_0_yzzz_1, g_yyyyyy_0_zzz_1, g_yyyyyy_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyy_0_xxxx_0[i] = 4.0 * g_yyyyyy_0_xxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxx_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxy_0[i] = 3.0 * g_yyyyyy_0_xxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxz_0[i] = 3.0 * g_yyyyyy_0_xxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyy_0[i] = 2.0 * g_yyyyyy_0_xyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyz_0[i] = 2.0 * g_yyyyyy_0_xyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxzz_0[i] = 2.0 * g_yyyyyy_0_xzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyy_0[i] = g_yyyyyy_0_yyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyz_0[i] = g_yyyyyy_0_yyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyzz_0[i] = g_yyyyyy_0_yzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xzzz_0[i] = g_yyyyyy_0_zzz_1[i] * fi_acd_0 + g_yyyyyy_0_xzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyy_0[i] = g_yyyyyy_0_yyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyz_0[i] = g_yyyyyy_0_yyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyzz_0[i] = g_yyyyyy_0_yyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yzzz_0[i] = g_yyyyyy_0_yzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_zzzz_0[i] = g_yyyyyy_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 330-345 components of targeted buffer : KSG

    auto g_xyyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 330);

    auto g_xyyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 331);

    auto g_xyyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 332);

    auto g_xyyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 333);

    auto g_xyyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 334);

    auto g_xyyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 335);

    auto g_xyyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 336);

    auto g_xyyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 337);

    auto g_xyyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 338);

    auto g_xyyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 339);

    auto g_xyyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 340);

    auto g_xyyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 341);

    auto g_xyyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 342);

    auto g_xyyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 343);

    auto g_xyyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 344);

    #pragma omp simd aligned(g_xyyyyy_0_xxxx_1, g_xyyyyy_0_xxxy_1, g_xyyyyy_0_xxyy_1, g_xyyyyy_0_xyyy_1, g_xyyyyyz_0_xxxx_0, g_xyyyyyz_0_xxxy_0, g_xyyyyyz_0_xxxz_0, g_xyyyyyz_0_xxyy_0, g_xyyyyyz_0_xxyz_0, g_xyyyyyz_0_xxzz_0, g_xyyyyyz_0_xyyy_0, g_xyyyyyz_0_xyyz_0, g_xyyyyyz_0_xyzz_0, g_xyyyyyz_0_xzzz_0, g_xyyyyyz_0_yyyy_0, g_xyyyyyz_0_yyyz_0, g_xyyyyyz_0_yyzz_0, g_xyyyyyz_0_yzzz_0, g_xyyyyyz_0_zzzz_0, g_yyyyyz_0_xxxz_1, g_yyyyyz_0_xxyz_1, g_yyyyyz_0_xxz_1, g_yyyyyz_0_xxzz_1, g_yyyyyz_0_xyyz_1, g_yyyyyz_0_xyz_1, g_yyyyyz_0_xyzz_1, g_yyyyyz_0_xzz_1, g_yyyyyz_0_xzzz_1, g_yyyyyz_0_yyyy_1, g_yyyyyz_0_yyyz_1, g_yyyyyz_0_yyz_1, g_yyyyyz_0_yyzz_1, g_yyyyyz_0_yzz_1, g_yyyyyz_0_yzzz_1, g_yyyyyz_0_zzz_1, g_yyyyyz_0_zzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyz_0_xxxx_0[i] = g_xyyyyy_0_xxxx_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxy_0[i] = g_xyyyyy_0_xxxy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxz_0[i] = 3.0 * g_yyyyyz_0_xxz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyy_0[i] = g_xyyyyy_0_xxyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxyz_0[i] = 2.0 * g_yyyyyz_0_xyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxzz_0[i] = 2.0 * g_yyyyyz_0_xzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyy_0[i] = g_xyyyyy_0_xyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xyyz_0[i] = g_yyyyyz_0_yyz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyzz_0[i] = g_yyyyyz_0_yzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xzzz_0[i] = g_yyyyyz_0_zzz_1[i] * fi_acd_0 + g_yyyyyz_0_xzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyy_0[i] = g_yyyyyz_0_yyyy_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyz_0[i] = g_yyyyyz_0_yyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyzz_0[i] = g_yyyyyz_0_yyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yzzz_0[i] = g_yyyyyz_0_yzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_zzzz_0[i] = g_yyyyyz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 345-360 components of targeted buffer : KSG

    auto g_xyyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 345);

    auto g_xyyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 346);

    auto g_xyyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 347);

    auto g_xyyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 348);

    auto g_xyyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 349);

    auto g_xyyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 350);

    auto g_xyyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 351);

    auto g_xyyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 352);

    auto g_xyyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 353);

    auto g_xyyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 354);

    auto g_xyyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 355);

    auto g_xyyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 356);

    auto g_xyyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 357);

    auto g_xyyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 358);

    auto g_xyyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 359);

    #pragma omp simd aligned(g_xyyyyzz_0_xxxx_0, g_xyyyyzz_0_xxxy_0, g_xyyyyzz_0_xxxz_0, g_xyyyyzz_0_xxyy_0, g_xyyyyzz_0_xxyz_0, g_xyyyyzz_0_xxzz_0, g_xyyyyzz_0_xyyy_0, g_xyyyyzz_0_xyyz_0, g_xyyyyzz_0_xyzz_0, g_xyyyyzz_0_xzzz_0, g_xyyyyzz_0_yyyy_0, g_xyyyyzz_0_yyyz_0, g_xyyyyzz_0_yyzz_0, g_xyyyyzz_0_yzzz_0, g_xyyyyzz_0_zzzz_0, g_yyyyzz_0_xxx_1, g_yyyyzz_0_xxxx_1, g_yyyyzz_0_xxxy_1, g_yyyyzz_0_xxxz_1, g_yyyyzz_0_xxy_1, g_yyyyzz_0_xxyy_1, g_yyyyzz_0_xxyz_1, g_yyyyzz_0_xxz_1, g_yyyyzz_0_xxzz_1, g_yyyyzz_0_xyy_1, g_yyyyzz_0_xyyy_1, g_yyyyzz_0_xyyz_1, g_yyyyzz_0_xyz_1, g_yyyyzz_0_xyzz_1, g_yyyyzz_0_xzz_1, g_yyyyzz_0_xzzz_1, g_yyyyzz_0_yyy_1, g_yyyyzz_0_yyyy_1, g_yyyyzz_0_yyyz_1, g_yyyyzz_0_yyz_1, g_yyyyzz_0_yyzz_1, g_yyyyzz_0_yzz_1, g_yyyyzz_0_yzzz_1, g_yyyyzz_0_zzz_1, g_yyyyzz_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyzz_0_xxxx_0[i] = 4.0 * g_yyyyzz_0_xxx_1[i] * fi_acd_0 + g_yyyyzz_0_xxxx_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxy_0[i] = 3.0 * g_yyyyzz_0_xxy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxz_0[i] = 3.0 * g_yyyyzz_0_xxz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyy_0[i] = 2.0 * g_yyyyzz_0_xyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyz_0[i] = 2.0 * g_yyyyzz_0_xyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxzz_0[i] = 2.0 * g_yyyyzz_0_xzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyy_0[i] = g_yyyyzz_0_yyy_1[i] * fi_acd_0 + g_yyyyzz_0_xyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyz_0[i] = g_yyyyzz_0_yyz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyzz_0[i] = g_yyyyzz_0_yzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xzzz_0[i] = g_yyyyzz_0_zzz_1[i] * fi_acd_0 + g_yyyyzz_0_xzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyy_0[i] = g_yyyyzz_0_yyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyz_0[i] = g_yyyyzz_0_yyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyzz_0[i] = g_yyyyzz_0_yyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yzzz_0[i] = g_yyyyzz_0_yzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_zzzz_0[i] = g_yyyyzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 360-375 components of targeted buffer : KSG

    auto g_xyyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 360);

    auto g_xyyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 361);

    auto g_xyyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 362);

    auto g_xyyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 363);

    auto g_xyyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 364);

    auto g_xyyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 365);

    auto g_xyyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 366);

    auto g_xyyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 367);

    auto g_xyyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 368);

    auto g_xyyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 369);

    auto g_xyyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 370);

    auto g_xyyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 371);

    auto g_xyyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 372);

    auto g_xyyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 373);

    auto g_xyyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 374);

    #pragma omp simd aligned(g_xyyyzzz_0_xxxx_0, g_xyyyzzz_0_xxxy_0, g_xyyyzzz_0_xxxz_0, g_xyyyzzz_0_xxyy_0, g_xyyyzzz_0_xxyz_0, g_xyyyzzz_0_xxzz_0, g_xyyyzzz_0_xyyy_0, g_xyyyzzz_0_xyyz_0, g_xyyyzzz_0_xyzz_0, g_xyyyzzz_0_xzzz_0, g_xyyyzzz_0_yyyy_0, g_xyyyzzz_0_yyyz_0, g_xyyyzzz_0_yyzz_0, g_xyyyzzz_0_yzzz_0, g_xyyyzzz_0_zzzz_0, g_yyyzzz_0_xxx_1, g_yyyzzz_0_xxxx_1, g_yyyzzz_0_xxxy_1, g_yyyzzz_0_xxxz_1, g_yyyzzz_0_xxy_1, g_yyyzzz_0_xxyy_1, g_yyyzzz_0_xxyz_1, g_yyyzzz_0_xxz_1, g_yyyzzz_0_xxzz_1, g_yyyzzz_0_xyy_1, g_yyyzzz_0_xyyy_1, g_yyyzzz_0_xyyz_1, g_yyyzzz_0_xyz_1, g_yyyzzz_0_xyzz_1, g_yyyzzz_0_xzz_1, g_yyyzzz_0_xzzz_1, g_yyyzzz_0_yyy_1, g_yyyzzz_0_yyyy_1, g_yyyzzz_0_yyyz_1, g_yyyzzz_0_yyz_1, g_yyyzzz_0_yyzz_1, g_yyyzzz_0_yzz_1, g_yyyzzz_0_yzzz_1, g_yyyzzz_0_zzz_1, g_yyyzzz_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzzz_0_xxxx_0[i] = 4.0 * g_yyyzzz_0_xxx_1[i] * fi_acd_0 + g_yyyzzz_0_xxxx_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxy_0[i] = 3.0 * g_yyyzzz_0_xxy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxz_0[i] = 3.0 * g_yyyzzz_0_xxz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyy_0[i] = 2.0 * g_yyyzzz_0_xyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyz_0[i] = 2.0 * g_yyyzzz_0_xyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxzz_0[i] = 2.0 * g_yyyzzz_0_xzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyy_0[i] = g_yyyzzz_0_yyy_1[i] * fi_acd_0 + g_yyyzzz_0_xyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyz_0[i] = g_yyyzzz_0_yyz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyzz_0[i] = g_yyyzzz_0_yzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xzzz_0[i] = g_yyyzzz_0_zzz_1[i] * fi_acd_0 + g_yyyzzz_0_xzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyy_0[i] = g_yyyzzz_0_yyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyz_0[i] = g_yyyzzz_0_yyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyzz_0[i] = g_yyyzzz_0_yyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yzzz_0[i] = g_yyyzzz_0_yzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_zzzz_0[i] = g_yyyzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 375-390 components of targeted buffer : KSG

    auto g_xyyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 375);

    auto g_xyyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 376);

    auto g_xyyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 377);

    auto g_xyyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 378);

    auto g_xyyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 379);

    auto g_xyyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 380);

    auto g_xyyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 381);

    auto g_xyyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 382);

    auto g_xyyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 383);

    auto g_xyyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 384);

    auto g_xyyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 385);

    auto g_xyyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 386);

    auto g_xyyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 387);

    auto g_xyyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 388);

    auto g_xyyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 389);

    #pragma omp simd aligned(g_xyyzzzz_0_xxxx_0, g_xyyzzzz_0_xxxy_0, g_xyyzzzz_0_xxxz_0, g_xyyzzzz_0_xxyy_0, g_xyyzzzz_0_xxyz_0, g_xyyzzzz_0_xxzz_0, g_xyyzzzz_0_xyyy_0, g_xyyzzzz_0_xyyz_0, g_xyyzzzz_0_xyzz_0, g_xyyzzzz_0_xzzz_0, g_xyyzzzz_0_yyyy_0, g_xyyzzzz_0_yyyz_0, g_xyyzzzz_0_yyzz_0, g_xyyzzzz_0_yzzz_0, g_xyyzzzz_0_zzzz_0, g_yyzzzz_0_xxx_1, g_yyzzzz_0_xxxx_1, g_yyzzzz_0_xxxy_1, g_yyzzzz_0_xxxz_1, g_yyzzzz_0_xxy_1, g_yyzzzz_0_xxyy_1, g_yyzzzz_0_xxyz_1, g_yyzzzz_0_xxz_1, g_yyzzzz_0_xxzz_1, g_yyzzzz_0_xyy_1, g_yyzzzz_0_xyyy_1, g_yyzzzz_0_xyyz_1, g_yyzzzz_0_xyz_1, g_yyzzzz_0_xyzz_1, g_yyzzzz_0_xzz_1, g_yyzzzz_0_xzzz_1, g_yyzzzz_0_yyy_1, g_yyzzzz_0_yyyy_1, g_yyzzzz_0_yyyz_1, g_yyzzzz_0_yyz_1, g_yyzzzz_0_yyzz_1, g_yyzzzz_0_yzz_1, g_yyzzzz_0_yzzz_1, g_yyzzzz_0_zzz_1, g_yyzzzz_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzzz_0_xxxx_0[i] = 4.0 * g_yyzzzz_0_xxx_1[i] * fi_acd_0 + g_yyzzzz_0_xxxx_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxy_0[i] = 3.0 * g_yyzzzz_0_xxy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxz_0[i] = 3.0 * g_yyzzzz_0_xxz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyy_0[i] = 2.0 * g_yyzzzz_0_xyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyz_0[i] = 2.0 * g_yyzzzz_0_xyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxzz_0[i] = 2.0 * g_yyzzzz_0_xzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyy_0[i] = g_yyzzzz_0_yyy_1[i] * fi_acd_0 + g_yyzzzz_0_xyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyz_0[i] = g_yyzzzz_0_yyz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyzz_0[i] = g_yyzzzz_0_yzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xzzz_0[i] = g_yyzzzz_0_zzz_1[i] * fi_acd_0 + g_yyzzzz_0_xzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyy_0[i] = g_yyzzzz_0_yyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyz_0[i] = g_yyzzzz_0_yyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyzz_0[i] = g_yyzzzz_0_yyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yzzz_0[i] = g_yyzzzz_0_yzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_zzzz_0[i] = g_yyzzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 390-405 components of targeted buffer : KSG

    auto g_xyzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 390);

    auto g_xyzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 391);

    auto g_xyzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 392);

    auto g_xyzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 393);

    auto g_xyzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 394);

    auto g_xyzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 395);

    auto g_xyzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 396);

    auto g_xyzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 397);

    auto g_xyzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 398);

    auto g_xyzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 399);

    auto g_xyzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 400);

    auto g_xyzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 401);

    auto g_xyzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 402);

    auto g_xyzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 403);

    auto g_xyzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 404);

    #pragma omp simd aligned(g_xyzzzzz_0_xxxx_0, g_xyzzzzz_0_xxxy_0, g_xyzzzzz_0_xxxz_0, g_xyzzzzz_0_xxyy_0, g_xyzzzzz_0_xxyz_0, g_xyzzzzz_0_xxzz_0, g_xyzzzzz_0_xyyy_0, g_xyzzzzz_0_xyyz_0, g_xyzzzzz_0_xyzz_0, g_xyzzzzz_0_xzzz_0, g_xyzzzzz_0_yyyy_0, g_xyzzzzz_0_yyyz_0, g_xyzzzzz_0_yyzz_0, g_xyzzzzz_0_yzzz_0, g_xyzzzzz_0_zzzz_0, g_xzzzzz_0_xxxx_1, g_xzzzzz_0_xxxz_1, g_xzzzzz_0_xxzz_1, g_xzzzzz_0_xzzz_1, g_yzzzzz_0_xxxy_1, g_yzzzzz_0_xxy_1, g_yzzzzz_0_xxyy_1, g_yzzzzz_0_xxyz_1, g_yzzzzz_0_xyy_1, g_yzzzzz_0_xyyy_1, g_yzzzzz_0_xyyz_1, g_yzzzzz_0_xyz_1, g_yzzzzz_0_xyzz_1, g_yzzzzz_0_yyy_1, g_yzzzzz_0_yyyy_1, g_yzzzzz_0_yyyz_1, g_yzzzzz_0_yyz_1, g_yzzzzz_0_yyzz_1, g_yzzzzz_0_yzz_1, g_yzzzzz_0_yzzz_1, g_yzzzzz_0_zzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzzz_0_xxxx_0[i] = g_xzzzzz_0_xxxx_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxy_0[i] = 3.0 * g_yzzzzz_0_xxy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxz_0[i] = g_xzzzzz_0_xxxz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxyy_0[i] = 2.0 * g_yzzzzz_0_xyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyz_0[i] = 2.0 * g_yzzzzz_0_xyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxzz_0[i] = g_xzzzzz_0_xxzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xyyy_0[i] = g_yzzzzz_0_yyy_1[i] * fi_acd_0 + g_yzzzzz_0_xyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyz_0[i] = g_yzzzzz_0_yyz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyzz_0[i] = g_yzzzzz_0_yzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xzzz_0[i] = g_xzzzzz_0_xzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_yyyy_0[i] = g_yzzzzz_0_yyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyz_0[i] = g_yzzzzz_0_yyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyzz_0[i] = g_yzzzzz_0_yyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yzzz_0[i] = g_yzzzzz_0_yzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_zzzz_0[i] = g_yzzzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 405-420 components of targeted buffer : KSG

    auto g_xzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 405);

    auto g_xzzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 406);

    auto g_xzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 407);

    auto g_xzzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 408);

    auto g_xzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 409);

    auto g_xzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 410);

    auto g_xzzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 411);

    auto g_xzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 412);

    auto g_xzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 413);

    auto g_xzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 414);

    auto g_xzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 415);

    auto g_xzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 416);

    auto g_xzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 417);

    auto g_xzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 418);

    auto g_xzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 419);

    #pragma omp simd aligned(g_xzzzzzz_0_xxxx_0, g_xzzzzzz_0_xxxy_0, g_xzzzzzz_0_xxxz_0, g_xzzzzzz_0_xxyy_0, g_xzzzzzz_0_xxyz_0, g_xzzzzzz_0_xxzz_0, g_xzzzzzz_0_xyyy_0, g_xzzzzzz_0_xyyz_0, g_xzzzzzz_0_xyzz_0, g_xzzzzzz_0_xzzz_0, g_xzzzzzz_0_yyyy_0, g_xzzzzzz_0_yyyz_0, g_xzzzzzz_0_yyzz_0, g_xzzzzzz_0_yzzz_0, g_xzzzzzz_0_zzzz_0, g_zzzzzz_0_xxx_1, g_zzzzzz_0_xxxx_1, g_zzzzzz_0_xxxy_1, g_zzzzzz_0_xxxz_1, g_zzzzzz_0_xxy_1, g_zzzzzz_0_xxyy_1, g_zzzzzz_0_xxyz_1, g_zzzzzz_0_xxz_1, g_zzzzzz_0_xxzz_1, g_zzzzzz_0_xyy_1, g_zzzzzz_0_xyyy_1, g_zzzzzz_0_xyyz_1, g_zzzzzz_0_xyz_1, g_zzzzzz_0_xyzz_1, g_zzzzzz_0_xzz_1, g_zzzzzz_0_xzzz_1, g_zzzzzz_0_yyy_1, g_zzzzzz_0_yyyy_1, g_zzzzzz_0_yyyz_1, g_zzzzzz_0_yyz_1, g_zzzzzz_0_yyzz_1, g_zzzzzz_0_yzz_1, g_zzzzzz_0_yzzz_1, g_zzzzzz_0_zzz_1, g_zzzzzz_0_zzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzzz_0_xxxx_0[i] = 4.0 * g_zzzzzz_0_xxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxx_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxy_0[i] = 3.0 * g_zzzzzz_0_xxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxz_0[i] = 3.0 * g_zzzzzz_0_xxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyy_0[i] = 2.0 * g_zzzzzz_0_xyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyz_0[i] = 2.0 * g_zzzzzz_0_xyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxzz_0[i] = 2.0 * g_zzzzzz_0_xzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyy_0[i] = g_zzzzzz_0_yyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyz_0[i] = g_zzzzzz_0_yyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyzz_0[i] = g_zzzzzz_0_yzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xzzz_0[i] = g_zzzzzz_0_zzz_1[i] * fi_acd_0 + g_zzzzzz_0_xzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyy_0[i] = g_zzzzzz_0_yyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyz_0[i] = g_zzzzzz_0_yyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyzz_0[i] = g_zzzzzz_0_yyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yzzz_0[i] = g_zzzzzz_0_yzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_zzzz_0[i] = g_zzzzzz_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 420-435 components of targeted buffer : KSG

    auto g_yyyyyyy_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 420);

    auto g_yyyyyyy_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 421);

    auto g_yyyyyyy_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 422);

    auto g_yyyyyyy_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 423);

    auto g_yyyyyyy_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 424);

    auto g_yyyyyyy_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 425);

    auto g_yyyyyyy_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 426);

    auto g_yyyyyyy_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 427);

    auto g_yyyyyyy_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 428);

    auto g_yyyyyyy_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 429);

    auto g_yyyyyyy_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 430);

    auto g_yyyyyyy_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 431);

    auto g_yyyyyyy_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 432);

    auto g_yyyyyyy_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 433);

    auto g_yyyyyyy_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 434);

    #pragma omp simd aligned(g_yyyyy_0_xxxx_0, g_yyyyy_0_xxxx_1, g_yyyyy_0_xxxy_0, g_yyyyy_0_xxxy_1, g_yyyyy_0_xxxz_0, g_yyyyy_0_xxxz_1, g_yyyyy_0_xxyy_0, g_yyyyy_0_xxyy_1, g_yyyyy_0_xxyz_0, g_yyyyy_0_xxyz_1, g_yyyyy_0_xxzz_0, g_yyyyy_0_xxzz_1, g_yyyyy_0_xyyy_0, g_yyyyy_0_xyyy_1, g_yyyyy_0_xyyz_0, g_yyyyy_0_xyyz_1, g_yyyyy_0_xyzz_0, g_yyyyy_0_xyzz_1, g_yyyyy_0_xzzz_0, g_yyyyy_0_xzzz_1, g_yyyyy_0_yyyy_0, g_yyyyy_0_yyyy_1, g_yyyyy_0_yyyz_0, g_yyyyy_0_yyyz_1, g_yyyyy_0_yyzz_0, g_yyyyy_0_yyzz_1, g_yyyyy_0_yzzz_0, g_yyyyy_0_yzzz_1, g_yyyyy_0_zzzz_0, g_yyyyy_0_zzzz_1, g_yyyyyy_0_xxx_1, g_yyyyyy_0_xxxx_1, g_yyyyyy_0_xxxy_1, g_yyyyyy_0_xxxz_1, g_yyyyyy_0_xxy_1, g_yyyyyy_0_xxyy_1, g_yyyyyy_0_xxyz_1, g_yyyyyy_0_xxz_1, g_yyyyyy_0_xxzz_1, g_yyyyyy_0_xyy_1, g_yyyyyy_0_xyyy_1, g_yyyyyy_0_xyyz_1, g_yyyyyy_0_xyz_1, g_yyyyyy_0_xyzz_1, g_yyyyyy_0_xzz_1, g_yyyyyy_0_xzzz_1, g_yyyyyy_0_yyy_1, g_yyyyyy_0_yyyy_1, g_yyyyyy_0_yyyz_1, g_yyyyyy_0_yyz_1, g_yyyyyy_0_yyzz_1, g_yyyyyy_0_yzz_1, g_yyyyyy_0_yzzz_1, g_yyyyyy_0_zzz_1, g_yyyyyy_0_zzzz_1, g_yyyyyyy_0_xxxx_0, g_yyyyyyy_0_xxxy_0, g_yyyyyyy_0_xxxz_0, g_yyyyyyy_0_xxyy_0, g_yyyyyyy_0_xxyz_0, g_yyyyyyy_0_xxzz_0, g_yyyyyyy_0_xyyy_0, g_yyyyyyy_0_xyyz_0, g_yyyyyyy_0_xyzz_0, g_yyyyyyy_0_xzzz_0, g_yyyyyyy_0_yyyy_0, g_yyyyyyy_0_yyyz_0, g_yyyyyyy_0_yyzz_0, g_yyyyyyy_0_yzzz_0, g_yyyyyyy_0_zzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyyy_0_xxxx_0[i] = 6.0 * g_yyyyy_0_xxxx_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxx_1[i] * fz_be_0 + g_yyyyyy_0_xxxx_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxy_0[i] = 6.0 * g_yyyyy_0_xxxy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxy_1[i] * fz_be_0 + g_yyyyyy_0_xxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxz_0[i] = 6.0 * g_yyyyy_0_xxxz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxz_1[i] * fz_be_0 + g_yyyyyy_0_xxxz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyy_0[i] = 6.0 * g_yyyyy_0_xxyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyy_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyz_0[i] = 6.0 * g_yyyyy_0_xxyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyz_1[i] * fz_be_0 + g_yyyyyy_0_xxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxzz_0[i] = 6.0 * g_yyyyy_0_xxzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxzz_1[i] * fz_be_0 + g_yyyyyy_0_xxzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyy_0[i] = 6.0 * g_yyyyy_0_xyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyy_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyz_0[i] = 6.0 * g_yyyyy_0_xyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyzz_0[i] = 6.0 * g_yyyyy_0_xyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyzz_1[i] * fz_be_0 + g_yyyyyy_0_xzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xzzz_0[i] = 6.0 * g_yyyyy_0_xzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xzzz_1[i] * fz_be_0 + g_yyyyyy_0_xzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyy_0[i] = 6.0 * g_yyyyy_0_yyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyy_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_yyy_1[i] * fi_acd_0 + g_yyyyyy_0_yyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyz_0[i] = 6.0 * g_yyyyy_0_yyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_yyz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyzz_0[i] = 6.0 * g_yyyyy_0_yyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_yzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yzzz_0[i] = 6.0 * g_yyyyy_0_yzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yzzz_1[i] * fz_be_0 + g_yyyyyy_0_zzz_1[i] * fi_acd_0 + g_yyyyyy_0_yzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_zzzz_0[i] = 6.0 * g_yyyyy_0_zzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_zzzz_1[i] * fz_be_0 + g_yyyyyy_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 435-450 components of targeted buffer : KSG

    auto g_yyyyyyz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 435);

    auto g_yyyyyyz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 436);

    auto g_yyyyyyz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 437);

    auto g_yyyyyyz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 438);

    auto g_yyyyyyz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 439);

    auto g_yyyyyyz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 440);

    auto g_yyyyyyz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 441);

    auto g_yyyyyyz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 442);

    auto g_yyyyyyz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 443);

    auto g_yyyyyyz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 444);

    auto g_yyyyyyz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 445);

    auto g_yyyyyyz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 446);

    auto g_yyyyyyz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 447);

    auto g_yyyyyyz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 448);

    auto g_yyyyyyz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 449);

    #pragma omp simd aligned(g_yyyyyy_0_xxx_1, g_yyyyyy_0_xxxx_1, g_yyyyyy_0_xxxy_1, g_yyyyyy_0_xxxz_1, g_yyyyyy_0_xxy_1, g_yyyyyy_0_xxyy_1, g_yyyyyy_0_xxyz_1, g_yyyyyy_0_xxz_1, g_yyyyyy_0_xxzz_1, g_yyyyyy_0_xyy_1, g_yyyyyy_0_xyyy_1, g_yyyyyy_0_xyyz_1, g_yyyyyy_0_xyz_1, g_yyyyyy_0_xyzz_1, g_yyyyyy_0_xzz_1, g_yyyyyy_0_xzzz_1, g_yyyyyy_0_yyy_1, g_yyyyyy_0_yyyy_1, g_yyyyyy_0_yyyz_1, g_yyyyyy_0_yyz_1, g_yyyyyy_0_yyzz_1, g_yyyyyy_0_yzz_1, g_yyyyyy_0_yzzz_1, g_yyyyyy_0_zzz_1, g_yyyyyy_0_zzzz_1, g_yyyyyyz_0_xxxx_0, g_yyyyyyz_0_xxxy_0, g_yyyyyyz_0_xxxz_0, g_yyyyyyz_0_xxyy_0, g_yyyyyyz_0_xxyz_0, g_yyyyyyz_0_xxzz_0, g_yyyyyyz_0_xyyy_0, g_yyyyyyz_0_xyyz_0, g_yyyyyyz_0_xyzz_0, g_yyyyyyz_0_xzzz_0, g_yyyyyyz_0_yyyy_0, g_yyyyyyz_0_yyyz_0, g_yyyyyyz_0_yyzz_0, g_yyyyyyz_0_yzzz_0, g_yyyyyyz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyyz_0_xxxx_0[i] = g_yyyyyy_0_xxxx_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxy_0[i] = g_yyyyyy_0_xxxy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxz_0[i] = g_yyyyyy_0_xxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyy_0[i] = g_yyyyyy_0_xxyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyz_0[i] = g_yyyyyy_0_xxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxzz_0[i] = 2.0 * g_yyyyyy_0_xxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyy_0[i] = g_yyyyyy_0_xyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyz_0[i] = g_yyyyyy_0_xyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyzz_0[i] = 2.0 * g_yyyyyy_0_xyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xzzz_0[i] = 3.0 * g_yyyyyy_0_xzz_1[i] * fi_acd_0 + g_yyyyyy_0_xzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyy_0[i] = g_yyyyyy_0_yyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyz_0[i] = g_yyyyyy_0_yyy_1[i] * fi_acd_0 + g_yyyyyy_0_yyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyzz_0[i] = 2.0 * g_yyyyyy_0_yyz_1[i] * fi_acd_0 + g_yyyyyy_0_yyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yzzz_0[i] = 3.0 * g_yyyyyy_0_yzz_1[i] * fi_acd_0 + g_yyyyyy_0_yzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_zzzz_0[i] = 4.0 * g_yyyyyy_0_zzz_1[i] * fi_acd_0 + g_yyyyyy_0_zzzz_1[i] * wa_z[i];
    }

    /// Set up 450-465 components of targeted buffer : KSG

    auto g_yyyyyzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 450);

    auto g_yyyyyzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 451);

    auto g_yyyyyzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 452);

    auto g_yyyyyzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 453);

    auto g_yyyyyzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 454);

    auto g_yyyyyzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 455);

    auto g_yyyyyzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 456);

    auto g_yyyyyzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 457);

    auto g_yyyyyzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 458);

    auto g_yyyyyzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 459);

    auto g_yyyyyzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 460);

    auto g_yyyyyzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 461);

    auto g_yyyyyzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 462);

    auto g_yyyyyzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 463);

    auto g_yyyyyzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 464);

    #pragma omp simd aligned(g_yyyyy_0_xxxy_0, g_yyyyy_0_xxxy_1, g_yyyyy_0_xxyy_0, g_yyyyy_0_xxyy_1, g_yyyyy_0_xyyy_0, g_yyyyy_0_xyyy_1, g_yyyyy_0_yyyy_0, g_yyyyy_0_yyyy_1, g_yyyyyz_0_xxxy_1, g_yyyyyz_0_xxyy_1, g_yyyyyz_0_xyyy_1, g_yyyyyz_0_yyyy_1, g_yyyyyzz_0_xxxx_0, g_yyyyyzz_0_xxxy_0, g_yyyyyzz_0_xxxz_0, g_yyyyyzz_0_xxyy_0, g_yyyyyzz_0_xxyz_0, g_yyyyyzz_0_xxzz_0, g_yyyyyzz_0_xyyy_0, g_yyyyyzz_0_xyyz_0, g_yyyyyzz_0_xyzz_0, g_yyyyyzz_0_xzzz_0, g_yyyyyzz_0_yyyy_0, g_yyyyyzz_0_yyyz_0, g_yyyyyzz_0_yyzz_0, g_yyyyyzz_0_yzzz_0, g_yyyyyzz_0_zzzz_0, g_yyyyzz_0_xxxx_1, g_yyyyzz_0_xxxz_1, g_yyyyzz_0_xxyz_1, g_yyyyzz_0_xxz_1, g_yyyyzz_0_xxzz_1, g_yyyyzz_0_xyyz_1, g_yyyyzz_0_xyz_1, g_yyyyzz_0_xyzz_1, g_yyyyzz_0_xzz_1, g_yyyyzz_0_xzzz_1, g_yyyyzz_0_yyyz_1, g_yyyyzz_0_yyz_1, g_yyyyzz_0_yyzz_1, g_yyyyzz_0_yzz_1, g_yyyyzz_0_yzzz_1, g_yyyyzz_0_zzz_1, g_yyyyzz_0_zzzz_1, g_yyyzz_0_xxxx_0, g_yyyzz_0_xxxx_1, g_yyyzz_0_xxxz_0, g_yyyzz_0_xxxz_1, g_yyyzz_0_xxyz_0, g_yyyzz_0_xxyz_1, g_yyyzz_0_xxzz_0, g_yyyzz_0_xxzz_1, g_yyyzz_0_xyyz_0, g_yyyzz_0_xyyz_1, g_yyyzz_0_xyzz_0, g_yyyzz_0_xyzz_1, g_yyyzz_0_xzzz_0, g_yyyzz_0_xzzz_1, g_yyyzz_0_yyyz_0, g_yyyzz_0_yyyz_1, g_yyyzz_0_yyzz_0, g_yyyzz_0_yyzz_1, g_yyyzz_0_yzzz_0, g_yyyzz_0_yzzz_1, g_yyyzz_0_zzzz_0, g_yyyzz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyzz_0_xxxx_0[i] = 4.0 * g_yyyzz_0_xxxx_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxx_1[i] * fz_be_0 + g_yyyyzz_0_xxxx_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxy_0[i] = g_yyyyy_0_xxxy_0[i] * fbe_0 - g_yyyyy_0_xxxy_1[i] * fz_be_0 + g_yyyyyz_0_xxxy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxz_0[i] = 4.0 * g_yyyzz_0_xxxz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxz_1[i] * fz_be_0 + g_yyyyzz_0_xxxz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyy_0[i] = g_yyyyy_0_xxyy_0[i] * fbe_0 - g_yyyyy_0_xxyy_1[i] * fz_be_0 + g_yyyyyz_0_xxyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxyz_0[i] = 4.0 * g_yyyzz_0_xxyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyz_1[i] * fz_be_0 + g_yyyyzz_0_xxz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxzz_0[i] = 4.0 * g_yyyzz_0_xxzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxzz_1[i] * fz_be_0 + g_yyyyzz_0_xxzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyy_0[i] = g_yyyyy_0_xyyy_0[i] * fbe_0 - g_yyyyy_0_xyyy_1[i] * fz_be_0 + g_yyyyyz_0_xyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xyyz_0[i] = 4.0 * g_yyyzz_0_xyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xyz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyzz_0[i] = 4.0 * g_yyyzz_0_xyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyzz_1[i] * fz_be_0 + g_yyyyzz_0_xzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xzzz_0[i] = 4.0 * g_yyyzz_0_xzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xzzz_1[i] * fz_be_0 + g_yyyyzz_0_xzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyy_0[i] = g_yyyyy_0_yyyy_0[i] * fbe_0 - g_yyyyy_0_yyyy_1[i] * fz_be_0 + g_yyyyyz_0_yyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_yyyz_0[i] = 4.0 * g_yyyzz_0_yyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_yyz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyzz_0[i] = 4.0 * g_yyyzz_0_yyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_yzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yzzz_0[i] = 4.0 * g_yyyzz_0_yzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yzzz_1[i] * fz_be_0 + g_yyyyzz_0_zzz_1[i] * fi_acd_0 + g_yyyyzz_0_yzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_zzzz_0[i] = 4.0 * g_yyyzz_0_zzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_zzzz_1[i] * fz_be_0 + g_yyyyzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 465-480 components of targeted buffer : KSG

    auto g_yyyyzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 465);

    auto g_yyyyzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 466);

    auto g_yyyyzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 467);

    auto g_yyyyzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 468);

    auto g_yyyyzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 469);

    auto g_yyyyzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 470);

    auto g_yyyyzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 471);

    auto g_yyyyzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 472);

    auto g_yyyyzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 473);

    auto g_yyyyzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 474);

    auto g_yyyyzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 475);

    auto g_yyyyzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 476);

    auto g_yyyyzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 477);

    auto g_yyyyzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 478);

    auto g_yyyyzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 479);

    #pragma omp simd aligned(g_yyyyz_0_xxxy_0, g_yyyyz_0_xxxy_1, g_yyyyz_0_xxyy_0, g_yyyyz_0_xxyy_1, g_yyyyz_0_xyyy_0, g_yyyyz_0_xyyy_1, g_yyyyz_0_yyyy_0, g_yyyyz_0_yyyy_1, g_yyyyzz_0_xxxy_1, g_yyyyzz_0_xxyy_1, g_yyyyzz_0_xyyy_1, g_yyyyzz_0_yyyy_1, g_yyyyzzz_0_xxxx_0, g_yyyyzzz_0_xxxy_0, g_yyyyzzz_0_xxxz_0, g_yyyyzzz_0_xxyy_0, g_yyyyzzz_0_xxyz_0, g_yyyyzzz_0_xxzz_0, g_yyyyzzz_0_xyyy_0, g_yyyyzzz_0_xyyz_0, g_yyyyzzz_0_xyzz_0, g_yyyyzzz_0_xzzz_0, g_yyyyzzz_0_yyyy_0, g_yyyyzzz_0_yyyz_0, g_yyyyzzz_0_yyzz_0, g_yyyyzzz_0_yzzz_0, g_yyyyzzz_0_zzzz_0, g_yyyzzz_0_xxxx_1, g_yyyzzz_0_xxxz_1, g_yyyzzz_0_xxyz_1, g_yyyzzz_0_xxz_1, g_yyyzzz_0_xxzz_1, g_yyyzzz_0_xyyz_1, g_yyyzzz_0_xyz_1, g_yyyzzz_0_xyzz_1, g_yyyzzz_0_xzz_1, g_yyyzzz_0_xzzz_1, g_yyyzzz_0_yyyz_1, g_yyyzzz_0_yyz_1, g_yyyzzz_0_yyzz_1, g_yyyzzz_0_yzz_1, g_yyyzzz_0_yzzz_1, g_yyyzzz_0_zzz_1, g_yyyzzz_0_zzzz_1, g_yyzzz_0_xxxx_0, g_yyzzz_0_xxxx_1, g_yyzzz_0_xxxz_0, g_yyzzz_0_xxxz_1, g_yyzzz_0_xxyz_0, g_yyzzz_0_xxyz_1, g_yyzzz_0_xxzz_0, g_yyzzz_0_xxzz_1, g_yyzzz_0_xyyz_0, g_yyzzz_0_xyyz_1, g_yyzzz_0_xyzz_0, g_yyzzz_0_xyzz_1, g_yyzzz_0_xzzz_0, g_yyzzz_0_xzzz_1, g_yyzzz_0_yyyz_0, g_yyzzz_0_yyyz_1, g_yyzzz_0_yyzz_0, g_yyzzz_0_yyzz_1, g_yyzzz_0_yzzz_0, g_yyzzz_0_yzzz_1, g_yyzzz_0_zzzz_0, g_yyzzz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzzz_0_xxxx_0[i] = 3.0 * g_yyzzz_0_xxxx_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxx_1[i] * fz_be_0 + g_yyyzzz_0_xxxx_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxy_0[i] = 2.0 * g_yyyyz_0_xxxy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxy_1[i] * fz_be_0 + g_yyyyzz_0_xxxy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxz_0[i] = 3.0 * g_yyzzz_0_xxxz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxz_1[i] * fz_be_0 + g_yyyzzz_0_xxxz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyy_0[i] = 2.0 * g_yyyyz_0_xxyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxyy_1[i] * fz_be_0 + g_yyyyzz_0_xxyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxyz_0[i] = 3.0 * g_yyzzz_0_xxyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyz_1[i] * fz_be_0 + g_yyyzzz_0_xxz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxzz_0[i] = 3.0 * g_yyzzz_0_xxzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxzz_1[i] * fz_be_0 + g_yyyzzz_0_xxzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyy_0[i] = 2.0 * g_yyyyz_0_xyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xyyy_1[i] * fz_be_0 + g_yyyyzz_0_xyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xyyz_0[i] = 3.0 * g_yyzzz_0_xyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xyz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyzz_0[i] = 3.0 * g_yyzzz_0_xyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyzz_1[i] * fz_be_0 + g_yyyzzz_0_xzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xzzz_0[i] = 3.0 * g_yyzzz_0_xzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xzzz_1[i] * fz_be_0 + g_yyyzzz_0_xzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyy_0[i] = 2.0 * g_yyyyz_0_yyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_yyyy_1[i] * fz_be_0 + g_yyyyzz_0_yyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_yyyz_0[i] = 3.0 * g_yyzzz_0_yyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_yyz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyzz_0[i] = 3.0 * g_yyzzz_0_yyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_yzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yzzz_0[i] = 3.0 * g_yyzzz_0_yzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yzzz_1[i] * fz_be_0 + g_yyyzzz_0_zzz_1[i] * fi_acd_0 + g_yyyzzz_0_yzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_zzzz_0[i] = 3.0 * g_yyzzz_0_zzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_zzzz_1[i] * fz_be_0 + g_yyyzzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 480-495 components of targeted buffer : KSG

    auto g_yyyzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 480);

    auto g_yyyzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 481);

    auto g_yyyzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 482);

    auto g_yyyzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 483);

    auto g_yyyzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 484);

    auto g_yyyzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 485);

    auto g_yyyzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 486);

    auto g_yyyzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 487);

    auto g_yyyzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 488);

    auto g_yyyzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 489);

    auto g_yyyzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 490);

    auto g_yyyzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 491);

    auto g_yyyzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 492);

    auto g_yyyzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 493);

    auto g_yyyzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 494);

    #pragma omp simd aligned(g_yyyzz_0_xxxy_0, g_yyyzz_0_xxxy_1, g_yyyzz_0_xxyy_0, g_yyyzz_0_xxyy_1, g_yyyzz_0_xyyy_0, g_yyyzz_0_xyyy_1, g_yyyzz_0_yyyy_0, g_yyyzz_0_yyyy_1, g_yyyzzz_0_xxxy_1, g_yyyzzz_0_xxyy_1, g_yyyzzz_0_xyyy_1, g_yyyzzz_0_yyyy_1, g_yyyzzzz_0_xxxx_0, g_yyyzzzz_0_xxxy_0, g_yyyzzzz_0_xxxz_0, g_yyyzzzz_0_xxyy_0, g_yyyzzzz_0_xxyz_0, g_yyyzzzz_0_xxzz_0, g_yyyzzzz_0_xyyy_0, g_yyyzzzz_0_xyyz_0, g_yyyzzzz_0_xyzz_0, g_yyyzzzz_0_xzzz_0, g_yyyzzzz_0_yyyy_0, g_yyyzzzz_0_yyyz_0, g_yyyzzzz_0_yyzz_0, g_yyyzzzz_0_yzzz_0, g_yyyzzzz_0_zzzz_0, g_yyzzzz_0_xxxx_1, g_yyzzzz_0_xxxz_1, g_yyzzzz_0_xxyz_1, g_yyzzzz_0_xxz_1, g_yyzzzz_0_xxzz_1, g_yyzzzz_0_xyyz_1, g_yyzzzz_0_xyz_1, g_yyzzzz_0_xyzz_1, g_yyzzzz_0_xzz_1, g_yyzzzz_0_xzzz_1, g_yyzzzz_0_yyyz_1, g_yyzzzz_0_yyz_1, g_yyzzzz_0_yyzz_1, g_yyzzzz_0_yzz_1, g_yyzzzz_0_yzzz_1, g_yyzzzz_0_zzz_1, g_yyzzzz_0_zzzz_1, g_yzzzz_0_xxxx_0, g_yzzzz_0_xxxx_1, g_yzzzz_0_xxxz_0, g_yzzzz_0_xxxz_1, g_yzzzz_0_xxyz_0, g_yzzzz_0_xxyz_1, g_yzzzz_0_xxzz_0, g_yzzzz_0_xxzz_1, g_yzzzz_0_xyyz_0, g_yzzzz_0_xyyz_1, g_yzzzz_0_xyzz_0, g_yzzzz_0_xyzz_1, g_yzzzz_0_xzzz_0, g_yzzzz_0_xzzz_1, g_yzzzz_0_yyyz_0, g_yzzzz_0_yyyz_1, g_yzzzz_0_yyzz_0, g_yzzzz_0_yyzz_1, g_yzzzz_0_yzzz_0, g_yzzzz_0_yzzz_1, g_yzzzz_0_zzzz_0, g_yzzzz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzzz_0_xxxx_0[i] = 2.0 * g_yzzzz_0_xxxx_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxx_1[i] * fz_be_0 + g_yyzzzz_0_xxxx_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxy_0[i] = 3.0 * g_yyyzz_0_xxxy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxy_1[i] * fz_be_0 + g_yyyzzz_0_xxxy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxz_0[i] = 2.0 * g_yzzzz_0_xxxz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxz_1[i] * fz_be_0 + g_yyzzzz_0_xxxz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyy_0[i] = 3.0 * g_yyyzz_0_xxyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxyy_1[i] * fz_be_0 + g_yyyzzz_0_xxyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxyz_0[i] = 2.0 * g_yzzzz_0_xxyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyz_1[i] * fz_be_0 + g_yyzzzz_0_xxz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxzz_0[i] = 2.0 * g_yzzzz_0_xxzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxzz_1[i] * fz_be_0 + g_yyzzzz_0_xxzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyy_0[i] = 3.0 * g_yyyzz_0_xyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xyyy_1[i] * fz_be_0 + g_yyyzzz_0_xyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xyyz_0[i] = 2.0 * g_yzzzz_0_xyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xyz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyzz_0[i] = 2.0 * g_yzzzz_0_xyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyzz_1[i] * fz_be_0 + g_yyzzzz_0_xzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xzzz_0[i] = 2.0 * g_yzzzz_0_xzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xzzz_1[i] * fz_be_0 + g_yyzzzz_0_xzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyy_0[i] = 3.0 * g_yyyzz_0_yyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_yyyy_1[i] * fz_be_0 + g_yyyzzz_0_yyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_yyyz_0[i] = 2.0 * g_yzzzz_0_yyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_yyz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyzz_0[i] = 2.0 * g_yzzzz_0_yyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_yzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yzzz_0[i] = 2.0 * g_yzzzz_0_yzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yzzz_1[i] * fz_be_0 + g_yyzzzz_0_zzz_1[i] * fi_acd_0 + g_yyzzzz_0_yzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_zzzz_0[i] = 2.0 * g_yzzzz_0_zzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_zzzz_1[i] * fz_be_0 + g_yyzzzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 495-510 components of targeted buffer : KSG

    auto g_yyzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 495);

    auto g_yyzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 496);

    auto g_yyzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 497);

    auto g_yyzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 498);

    auto g_yyzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 499);

    auto g_yyzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 500);

    auto g_yyzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 501);

    auto g_yyzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 502);

    auto g_yyzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 503);

    auto g_yyzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 504);

    auto g_yyzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 505);

    auto g_yyzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 506);

    auto g_yyzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 507);

    auto g_yyzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 508);

    auto g_yyzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 509);

    #pragma omp simd aligned(g_yyzzz_0_xxxy_0, g_yyzzz_0_xxxy_1, g_yyzzz_0_xxyy_0, g_yyzzz_0_xxyy_1, g_yyzzz_0_xyyy_0, g_yyzzz_0_xyyy_1, g_yyzzz_0_yyyy_0, g_yyzzz_0_yyyy_1, g_yyzzzz_0_xxxy_1, g_yyzzzz_0_xxyy_1, g_yyzzzz_0_xyyy_1, g_yyzzzz_0_yyyy_1, g_yyzzzzz_0_xxxx_0, g_yyzzzzz_0_xxxy_0, g_yyzzzzz_0_xxxz_0, g_yyzzzzz_0_xxyy_0, g_yyzzzzz_0_xxyz_0, g_yyzzzzz_0_xxzz_0, g_yyzzzzz_0_xyyy_0, g_yyzzzzz_0_xyyz_0, g_yyzzzzz_0_xyzz_0, g_yyzzzzz_0_xzzz_0, g_yyzzzzz_0_yyyy_0, g_yyzzzzz_0_yyyz_0, g_yyzzzzz_0_yyzz_0, g_yyzzzzz_0_yzzz_0, g_yyzzzzz_0_zzzz_0, g_yzzzzz_0_xxxx_1, g_yzzzzz_0_xxxz_1, g_yzzzzz_0_xxyz_1, g_yzzzzz_0_xxz_1, g_yzzzzz_0_xxzz_1, g_yzzzzz_0_xyyz_1, g_yzzzzz_0_xyz_1, g_yzzzzz_0_xyzz_1, g_yzzzzz_0_xzz_1, g_yzzzzz_0_xzzz_1, g_yzzzzz_0_yyyz_1, g_yzzzzz_0_yyz_1, g_yzzzzz_0_yyzz_1, g_yzzzzz_0_yzz_1, g_yzzzzz_0_yzzz_1, g_yzzzzz_0_zzz_1, g_yzzzzz_0_zzzz_1, g_zzzzz_0_xxxx_0, g_zzzzz_0_xxxx_1, g_zzzzz_0_xxxz_0, g_zzzzz_0_xxxz_1, g_zzzzz_0_xxyz_0, g_zzzzz_0_xxyz_1, g_zzzzz_0_xxzz_0, g_zzzzz_0_xxzz_1, g_zzzzz_0_xyyz_0, g_zzzzz_0_xyyz_1, g_zzzzz_0_xyzz_0, g_zzzzz_0_xyzz_1, g_zzzzz_0_xzzz_0, g_zzzzz_0_xzzz_1, g_zzzzz_0_yyyz_0, g_zzzzz_0_yyyz_1, g_zzzzz_0_yyzz_0, g_zzzzz_0_yyzz_1, g_zzzzz_0_yzzz_0, g_zzzzz_0_yzzz_1, g_zzzzz_0_zzzz_0, g_zzzzz_0_zzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzzz_0_xxxx_0[i] = g_zzzzz_0_xxxx_0[i] * fbe_0 - g_zzzzz_0_xxxx_1[i] * fz_be_0 + g_yzzzzz_0_xxxx_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxy_0[i] = 4.0 * g_yyzzz_0_xxxy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxy_1[i] * fz_be_0 + g_yyzzzz_0_xxxy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxz_0[i] = g_zzzzz_0_xxxz_0[i] * fbe_0 - g_zzzzz_0_xxxz_1[i] * fz_be_0 + g_yzzzzz_0_xxxz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyy_0[i] = 4.0 * g_yyzzz_0_xxyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxyy_1[i] * fz_be_0 + g_yyzzzz_0_xxyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxyz_0[i] = g_zzzzz_0_xxyz_0[i] * fbe_0 - g_zzzzz_0_xxyz_1[i] * fz_be_0 + g_yzzzzz_0_xxz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxzz_0[i] = g_zzzzz_0_xxzz_0[i] * fbe_0 - g_zzzzz_0_xxzz_1[i] * fz_be_0 + g_yzzzzz_0_xxzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyy_0[i] = 4.0 * g_yyzzz_0_xyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xyyy_1[i] * fz_be_0 + g_yyzzzz_0_xyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xyyz_0[i] = g_zzzzz_0_xyyz_0[i] * fbe_0 - g_zzzzz_0_xyyz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xyz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyzz_0[i] = g_zzzzz_0_xyzz_0[i] * fbe_0 - g_zzzzz_0_xyzz_1[i] * fz_be_0 + g_yzzzzz_0_xzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xzzz_0[i] = g_zzzzz_0_xzzz_0[i] * fbe_0 - g_zzzzz_0_xzzz_1[i] * fz_be_0 + g_yzzzzz_0_xzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyy_0[i] = 4.0 * g_yyzzz_0_yyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_yyyy_1[i] * fz_be_0 + g_yyzzzz_0_yyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_yyyz_0[i] = g_zzzzz_0_yyyz_0[i] * fbe_0 - g_zzzzz_0_yyyz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_yyz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyzz_0[i] = g_zzzzz_0_yyzz_0[i] * fbe_0 - g_zzzzz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_yzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yzzz_0[i] = g_zzzzz_0_yzzz_0[i] * fbe_0 - g_zzzzz_0_yzzz_1[i] * fz_be_0 + g_yzzzzz_0_zzz_1[i] * fi_acd_0 + g_yzzzzz_0_yzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_zzzz_0[i] = g_zzzzz_0_zzzz_0[i] * fbe_0 - g_zzzzz_0_zzzz_1[i] * fz_be_0 + g_yzzzzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 510-525 components of targeted buffer : KSG

    auto g_yzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 510);

    auto g_yzzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 511);

    auto g_yzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 512);

    auto g_yzzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 513);

    auto g_yzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 514);

    auto g_yzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 515);

    auto g_yzzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 516);

    auto g_yzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 517);

    auto g_yzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 518);

    auto g_yzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 519);

    auto g_yzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 520);

    auto g_yzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 521);

    auto g_yzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 522);

    auto g_yzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 523);

    auto g_yzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 524);

    #pragma omp simd aligned(g_yzzzzzz_0_xxxx_0, g_yzzzzzz_0_xxxy_0, g_yzzzzzz_0_xxxz_0, g_yzzzzzz_0_xxyy_0, g_yzzzzzz_0_xxyz_0, g_yzzzzzz_0_xxzz_0, g_yzzzzzz_0_xyyy_0, g_yzzzzzz_0_xyyz_0, g_yzzzzzz_0_xyzz_0, g_yzzzzzz_0_xzzz_0, g_yzzzzzz_0_yyyy_0, g_yzzzzzz_0_yyyz_0, g_yzzzzzz_0_yyzz_0, g_yzzzzzz_0_yzzz_0, g_yzzzzzz_0_zzzz_0, g_zzzzzz_0_xxx_1, g_zzzzzz_0_xxxx_1, g_zzzzzz_0_xxxy_1, g_zzzzzz_0_xxxz_1, g_zzzzzz_0_xxy_1, g_zzzzzz_0_xxyy_1, g_zzzzzz_0_xxyz_1, g_zzzzzz_0_xxz_1, g_zzzzzz_0_xxzz_1, g_zzzzzz_0_xyy_1, g_zzzzzz_0_xyyy_1, g_zzzzzz_0_xyyz_1, g_zzzzzz_0_xyz_1, g_zzzzzz_0_xyzz_1, g_zzzzzz_0_xzz_1, g_zzzzzz_0_xzzz_1, g_zzzzzz_0_yyy_1, g_zzzzzz_0_yyyy_1, g_zzzzzz_0_yyyz_1, g_zzzzzz_0_yyz_1, g_zzzzzz_0_yyzz_1, g_zzzzzz_0_yzz_1, g_zzzzzz_0_yzzz_1, g_zzzzzz_0_zzz_1, g_zzzzzz_0_zzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzzz_0_xxxx_0[i] = g_zzzzzz_0_xxxx_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxy_0[i] = g_zzzzzz_0_xxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxz_0[i] = g_zzzzzz_0_xxxz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyy_0[i] = 2.0 * g_zzzzzz_0_xxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyz_0[i] = g_zzzzzz_0_xxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxzz_0[i] = g_zzzzzz_0_xxzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyy_0[i] = 3.0 * g_zzzzzz_0_xyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyz_0[i] = 2.0 * g_zzzzzz_0_xyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyzz_0[i] = g_zzzzzz_0_xzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xzzz_0[i] = g_zzzzzz_0_xzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyy_0[i] = 4.0 * g_zzzzzz_0_yyy_1[i] * fi_acd_0 + g_zzzzzz_0_yyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyz_0[i] = 3.0 * g_zzzzzz_0_yyz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyzz_0[i] = 2.0 * g_zzzzzz_0_yzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yzzz_0[i] = g_zzzzzz_0_zzz_1[i] * fi_acd_0 + g_zzzzzz_0_yzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_zzzz_0[i] = g_zzzzzz_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 525-540 components of targeted buffer : KSG

    auto g_zzzzzzz_0_xxxx_0 = pbuffer.data(idx_eri_0_ksg + 525);

    auto g_zzzzzzz_0_xxxy_0 = pbuffer.data(idx_eri_0_ksg + 526);

    auto g_zzzzzzz_0_xxxz_0 = pbuffer.data(idx_eri_0_ksg + 527);

    auto g_zzzzzzz_0_xxyy_0 = pbuffer.data(idx_eri_0_ksg + 528);

    auto g_zzzzzzz_0_xxyz_0 = pbuffer.data(idx_eri_0_ksg + 529);

    auto g_zzzzzzz_0_xxzz_0 = pbuffer.data(idx_eri_0_ksg + 530);

    auto g_zzzzzzz_0_xyyy_0 = pbuffer.data(idx_eri_0_ksg + 531);

    auto g_zzzzzzz_0_xyyz_0 = pbuffer.data(idx_eri_0_ksg + 532);

    auto g_zzzzzzz_0_xyzz_0 = pbuffer.data(idx_eri_0_ksg + 533);

    auto g_zzzzzzz_0_xzzz_0 = pbuffer.data(idx_eri_0_ksg + 534);

    auto g_zzzzzzz_0_yyyy_0 = pbuffer.data(idx_eri_0_ksg + 535);

    auto g_zzzzzzz_0_yyyz_0 = pbuffer.data(idx_eri_0_ksg + 536);

    auto g_zzzzzzz_0_yyzz_0 = pbuffer.data(idx_eri_0_ksg + 537);

    auto g_zzzzzzz_0_yzzz_0 = pbuffer.data(idx_eri_0_ksg + 538);

    auto g_zzzzzzz_0_zzzz_0 = pbuffer.data(idx_eri_0_ksg + 539);

    #pragma omp simd aligned(g_zzzzz_0_xxxx_0, g_zzzzz_0_xxxx_1, g_zzzzz_0_xxxy_0, g_zzzzz_0_xxxy_1, g_zzzzz_0_xxxz_0, g_zzzzz_0_xxxz_1, g_zzzzz_0_xxyy_0, g_zzzzz_0_xxyy_1, g_zzzzz_0_xxyz_0, g_zzzzz_0_xxyz_1, g_zzzzz_0_xxzz_0, g_zzzzz_0_xxzz_1, g_zzzzz_0_xyyy_0, g_zzzzz_0_xyyy_1, g_zzzzz_0_xyyz_0, g_zzzzz_0_xyyz_1, g_zzzzz_0_xyzz_0, g_zzzzz_0_xyzz_1, g_zzzzz_0_xzzz_0, g_zzzzz_0_xzzz_1, g_zzzzz_0_yyyy_0, g_zzzzz_0_yyyy_1, g_zzzzz_0_yyyz_0, g_zzzzz_0_yyyz_1, g_zzzzz_0_yyzz_0, g_zzzzz_0_yyzz_1, g_zzzzz_0_yzzz_0, g_zzzzz_0_yzzz_1, g_zzzzz_0_zzzz_0, g_zzzzz_0_zzzz_1, g_zzzzzz_0_xxx_1, g_zzzzzz_0_xxxx_1, g_zzzzzz_0_xxxy_1, g_zzzzzz_0_xxxz_1, g_zzzzzz_0_xxy_1, g_zzzzzz_0_xxyy_1, g_zzzzzz_0_xxyz_1, g_zzzzzz_0_xxz_1, g_zzzzzz_0_xxzz_1, g_zzzzzz_0_xyy_1, g_zzzzzz_0_xyyy_1, g_zzzzzz_0_xyyz_1, g_zzzzzz_0_xyz_1, g_zzzzzz_0_xyzz_1, g_zzzzzz_0_xzz_1, g_zzzzzz_0_xzzz_1, g_zzzzzz_0_yyy_1, g_zzzzzz_0_yyyy_1, g_zzzzzz_0_yyyz_1, g_zzzzzz_0_yyz_1, g_zzzzzz_0_yyzz_1, g_zzzzzz_0_yzz_1, g_zzzzzz_0_yzzz_1, g_zzzzzz_0_zzz_1, g_zzzzzz_0_zzzz_1, g_zzzzzzz_0_xxxx_0, g_zzzzzzz_0_xxxy_0, g_zzzzzzz_0_xxxz_0, g_zzzzzzz_0_xxyy_0, g_zzzzzzz_0_xxyz_0, g_zzzzzzz_0_xxzz_0, g_zzzzzzz_0_xyyy_0, g_zzzzzzz_0_xyyz_0, g_zzzzzzz_0_xyzz_0, g_zzzzzzz_0_xzzz_0, g_zzzzzzz_0_yyyy_0, g_zzzzzzz_0_yyyz_0, g_zzzzzzz_0_yyzz_0, g_zzzzzzz_0_yzzz_0, g_zzzzzzz_0_zzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzzz_0_xxxx_0[i] = 6.0 * g_zzzzz_0_xxxx_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxx_1[i] * fz_be_0 + g_zzzzzz_0_xxxx_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxy_0[i] = 6.0 * g_zzzzz_0_xxxy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxy_1[i] * fz_be_0 + g_zzzzzz_0_xxxy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxz_0[i] = 6.0 * g_zzzzz_0_xxxz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxz_1[i] * fz_be_0 + g_zzzzzz_0_xxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyy_0[i] = 6.0 * g_zzzzz_0_xxyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyy_1[i] * fz_be_0 + g_zzzzzz_0_xxyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyz_0[i] = 6.0 * g_zzzzz_0_xxyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyz_1[i] * fz_be_0 + g_zzzzzz_0_xxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxzz_0[i] = 6.0 * g_zzzzz_0_xxzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyy_0[i] = 6.0 * g_zzzzz_0_xyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyy_1[i] * fz_be_0 + g_zzzzzz_0_xyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyz_0[i] = 6.0 * g_zzzzz_0_xyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyz_1[i] * fz_be_0 + g_zzzzzz_0_xyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyzz_0[i] = 6.0 * g_zzzzz_0_xyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xzzz_0[i] = 6.0 * g_zzzzz_0_xzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xzz_1[i] * fi_acd_0 + g_zzzzzz_0_xzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyy_0[i] = 6.0 * g_zzzzz_0_yyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyy_1[i] * fz_be_0 + g_zzzzzz_0_yyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyz_0[i] = 6.0 * g_zzzzz_0_yyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyz_1[i] * fz_be_0 + g_zzzzzz_0_yyy_1[i] * fi_acd_0 + g_zzzzzz_0_yyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyzz_0[i] = 6.0 * g_zzzzz_0_yyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_yyz_1[i] * fi_acd_0 + g_zzzzzz_0_yyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yzzz_0[i] = 6.0 * g_zzzzz_0_yzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_yzz_1[i] * fi_acd_0 + g_zzzzzz_0_yzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_zzzz_0[i] = 6.0 * g_zzzzz_0_zzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_zzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_zzz_1[i] * fi_acd_0 + g_zzzzzz_0_zzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

