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

#include "TwoCenterElectronRepulsionPrimRecIH.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ih(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ih,
                                const size_t idx_eri_0_gh,
                                const size_t idx_eri_1_gh,
                                const size_t idx_eri_1_hg,
                                const size_t idx_eri_1_hh,
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

    auto g_xxxx_xxxxx_0 = pbuffer.data(idx_eri_0_gh);

    auto g_xxxx_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 1);

    auto g_xxxx_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 2);

    auto g_xxxx_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 3);

    auto g_xxxx_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 4);

    auto g_xxxx_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 5);

    auto g_xxxx_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 6);

    auto g_xxxx_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 7);

    auto g_xxxx_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 8);

    auto g_xxxx_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 9);

    auto g_xxxx_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 10);

    auto g_xxxx_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 11);

    auto g_xxxx_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 12);

    auto g_xxxx_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 13);

    auto g_xxxx_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 14);

    auto g_xxxx_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 15);

    auto g_xxxx_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 16);

    auto g_xxxx_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 17);

    auto g_xxxx_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 18);

    auto g_xxxx_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 19);

    auto g_xxxx_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 20);

    auto g_xxxy_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 21);

    auto g_xxxy_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 23);

    auto g_xxxy_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 26);

    auto g_xxxy_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 30);

    auto g_xxxy_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 35);

    auto g_xxxz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 42);

    auto g_xxxz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 43);

    auto g_xxxz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 45);

    auto g_xxxz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 48);

    auto g_xxxz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 52);

    auto g_xxyy_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 63);

    auto g_xxyy_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 64);

    auto g_xxyy_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 65);

    auto g_xxyy_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 66);

    auto g_xxyy_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 67);

    auto g_xxyy_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 68);

    auto g_xxyy_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 69);

    auto g_xxyy_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 70);

    auto g_xxyy_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 71);

    auto g_xxyy_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 72);

    auto g_xxyy_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 73);

    auto g_xxyy_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 74);

    auto g_xxyy_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 75);

    auto g_xxyy_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 76);

    auto g_xxyy_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 77);

    auto g_xxyy_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 78);

    auto g_xxyy_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 79);

    auto g_xxyy_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 80);

    auto g_xxyy_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 81);

    auto g_xxyy_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 82);

    auto g_xxyy_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 83);

    auto g_xxzz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 105);

    auto g_xxzz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 106);

    auto g_xxzz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 107);

    auto g_xxzz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 108);

    auto g_xxzz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 109);

    auto g_xxzz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 110);

    auto g_xxzz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 111);

    auto g_xxzz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 112);

    auto g_xxzz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 113);

    auto g_xxzz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 114);

    auto g_xxzz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 115);

    auto g_xxzz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 116);

    auto g_xxzz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 117);

    auto g_xxzz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 118);

    auto g_xxzz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 119);

    auto g_xxzz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 120);

    auto g_xxzz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 121);

    auto g_xxzz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 122);

    auto g_xxzz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 123);

    auto g_xxzz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 124);

    auto g_xxzz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 125);

    auto g_xyyy_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 127);

    auto g_xyyy_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 129);

    auto g_xyyy_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 130);

    auto g_xyyy_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 132);

    auto g_xyyy_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 133);

    auto g_xyyy_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 134);

    auto g_xyyy_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 136);

    auto g_xyyy_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 137);

    auto g_xyyy_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 138);

    auto g_xyyy_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 139);

    auto g_xyyy_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 141);

    auto g_xyyy_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 142);

    auto g_xyyy_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 143);

    auto g_xyyy_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 144);

    auto g_xyyy_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 145);

    auto g_xyyy_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 146);

    auto g_xzzz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 191);

    auto g_xzzz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 193);

    auto g_xzzz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 194);

    auto g_xzzz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 196);

    auto g_xzzz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 197);

    auto g_xzzz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 198);

    auto g_xzzz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 200);

    auto g_xzzz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 201);

    auto g_xzzz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 202);

    auto g_xzzz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 203);

    auto g_xzzz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 204);

    auto g_xzzz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 205);

    auto g_xzzz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 206);

    auto g_xzzz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 207);

    auto g_xzzz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 208);

    auto g_xzzz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 209);

    auto g_yyyy_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 210);

    auto g_yyyy_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 211);

    auto g_yyyy_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 212);

    auto g_yyyy_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 213);

    auto g_yyyy_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 214);

    auto g_yyyy_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 215);

    auto g_yyyy_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 216);

    auto g_yyyy_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 217);

    auto g_yyyy_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 218);

    auto g_yyyy_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 219);

    auto g_yyyy_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 220);

    auto g_yyyy_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 221);

    auto g_yyyy_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 222);

    auto g_yyyy_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 223);

    auto g_yyyy_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 224);

    auto g_yyyy_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 225);

    auto g_yyyy_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 226);

    auto g_yyyy_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 227);

    auto g_yyyy_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 228);

    auto g_yyyy_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 229);

    auto g_yyyy_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 230);

    auto g_yyyz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 232);

    auto g_yyyz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 234);

    auto g_yyyz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 237);

    auto g_yyyz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 241);

    auto g_yyyz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 246);

    auto g_yyzz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 252);

    auto g_yyzz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 253);

    auto g_yyzz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 254);

    auto g_yyzz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 255);

    auto g_yyzz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 256);

    auto g_yyzz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 257);

    auto g_yyzz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 258);

    auto g_yyzz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 259);

    auto g_yyzz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 260);

    auto g_yyzz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 261);

    auto g_yyzz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 262);

    auto g_yyzz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 263);

    auto g_yyzz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 264);

    auto g_yyzz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 265);

    auto g_yyzz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 266);

    auto g_yyzz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 267);

    auto g_yyzz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 268);

    auto g_yyzz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 269);

    auto g_yyzz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 270);

    auto g_yyzz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 271);

    auto g_yyzz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 272);

    auto g_yzzz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 273);

    auto g_yzzz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 275);

    auto g_yzzz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 277);

    auto g_yzzz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 278);

    auto g_yzzz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 280);

    auto g_yzzz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 281);

    auto g_yzzz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 282);

    auto g_yzzz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 284);

    auto g_yzzz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 285);

    auto g_yzzz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 286);

    auto g_yzzz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 287);

    auto g_yzzz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 289);

    auto g_yzzz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 290);

    auto g_yzzz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 291);

    auto g_yzzz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 292);

    auto g_yzzz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 293);

    auto g_zzzz_xxxxx_0 = pbuffer.data(idx_eri_0_gh + 294);

    auto g_zzzz_xxxxy_0 = pbuffer.data(idx_eri_0_gh + 295);

    auto g_zzzz_xxxxz_0 = pbuffer.data(idx_eri_0_gh + 296);

    auto g_zzzz_xxxyy_0 = pbuffer.data(idx_eri_0_gh + 297);

    auto g_zzzz_xxxyz_0 = pbuffer.data(idx_eri_0_gh + 298);

    auto g_zzzz_xxxzz_0 = pbuffer.data(idx_eri_0_gh + 299);

    auto g_zzzz_xxyyy_0 = pbuffer.data(idx_eri_0_gh + 300);

    auto g_zzzz_xxyyz_0 = pbuffer.data(idx_eri_0_gh + 301);

    auto g_zzzz_xxyzz_0 = pbuffer.data(idx_eri_0_gh + 302);

    auto g_zzzz_xxzzz_0 = pbuffer.data(idx_eri_0_gh + 303);

    auto g_zzzz_xyyyy_0 = pbuffer.data(idx_eri_0_gh + 304);

    auto g_zzzz_xyyyz_0 = pbuffer.data(idx_eri_0_gh + 305);

    auto g_zzzz_xyyzz_0 = pbuffer.data(idx_eri_0_gh + 306);

    auto g_zzzz_xyzzz_0 = pbuffer.data(idx_eri_0_gh + 307);

    auto g_zzzz_xzzzz_0 = pbuffer.data(idx_eri_0_gh + 308);

    auto g_zzzz_yyyyy_0 = pbuffer.data(idx_eri_0_gh + 309);

    auto g_zzzz_yyyyz_0 = pbuffer.data(idx_eri_0_gh + 310);

    auto g_zzzz_yyyzz_0 = pbuffer.data(idx_eri_0_gh + 311);

    auto g_zzzz_yyzzz_0 = pbuffer.data(idx_eri_0_gh + 312);

    auto g_zzzz_yzzzz_0 = pbuffer.data(idx_eri_0_gh + 313);

    auto g_zzzz_zzzzz_0 = pbuffer.data(idx_eri_0_gh + 314);

    // Set up components of auxiliary buffer : GH

    auto g_xxxx_xxxxx_1 = pbuffer.data(idx_eri_1_gh);

    auto g_xxxx_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 1);

    auto g_xxxx_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 2);

    auto g_xxxx_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 3);

    auto g_xxxx_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 4);

    auto g_xxxx_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 5);

    auto g_xxxx_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 6);

    auto g_xxxx_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 7);

    auto g_xxxx_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 8);

    auto g_xxxx_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 9);

    auto g_xxxx_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 10);

    auto g_xxxx_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 11);

    auto g_xxxx_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 12);

    auto g_xxxx_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 13);

    auto g_xxxx_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 14);

    auto g_xxxx_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 15);

    auto g_xxxx_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 16);

    auto g_xxxx_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 17);

    auto g_xxxx_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 18);

    auto g_xxxx_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 19);

    auto g_xxxx_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 20);

    auto g_xxxy_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 21);

    auto g_xxxy_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 23);

    auto g_xxxy_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 26);

    auto g_xxxy_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 30);

    auto g_xxxy_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 35);

    auto g_xxxz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 42);

    auto g_xxxz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 43);

    auto g_xxxz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 45);

    auto g_xxxz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 48);

    auto g_xxxz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 52);

    auto g_xxyy_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 63);

    auto g_xxyy_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 64);

    auto g_xxyy_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 65);

    auto g_xxyy_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 66);

    auto g_xxyy_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 67);

    auto g_xxyy_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 68);

    auto g_xxyy_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 69);

    auto g_xxyy_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 70);

    auto g_xxyy_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 71);

    auto g_xxyy_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 72);

    auto g_xxyy_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 73);

    auto g_xxyy_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 74);

    auto g_xxyy_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 75);

    auto g_xxyy_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 76);

    auto g_xxyy_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 77);

    auto g_xxyy_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 78);

    auto g_xxyy_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 79);

    auto g_xxyy_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 80);

    auto g_xxyy_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 81);

    auto g_xxyy_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 82);

    auto g_xxyy_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 83);

    auto g_xxzz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 105);

    auto g_xxzz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 106);

    auto g_xxzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 107);

    auto g_xxzz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 108);

    auto g_xxzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 109);

    auto g_xxzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 110);

    auto g_xxzz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 111);

    auto g_xxzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 112);

    auto g_xxzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 113);

    auto g_xxzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 114);

    auto g_xxzz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 115);

    auto g_xxzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 116);

    auto g_xxzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 117);

    auto g_xxzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 118);

    auto g_xxzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 119);

    auto g_xxzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 120);

    auto g_xxzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 121);

    auto g_xxzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 122);

    auto g_xxzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 123);

    auto g_xxzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 124);

    auto g_xxzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 125);

    auto g_xyyy_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 127);

    auto g_xyyy_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 129);

    auto g_xyyy_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 130);

    auto g_xyyy_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 132);

    auto g_xyyy_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 133);

    auto g_xyyy_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 134);

    auto g_xyyy_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 136);

    auto g_xyyy_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 137);

    auto g_xyyy_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 138);

    auto g_xyyy_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 139);

    auto g_xyyy_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 141);

    auto g_xyyy_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 142);

    auto g_xyyy_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 143);

    auto g_xyyy_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 144);

    auto g_xyyy_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 145);

    auto g_xyyy_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 146);

    auto g_xzzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 191);

    auto g_xzzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 193);

    auto g_xzzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 194);

    auto g_xzzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 196);

    auto g_xzzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 197);

    auto g_xzzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 198);

    auto g_xzzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 200);

    auto g_xzzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 201);

    auto g_xzzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 202);

    auto g_xzzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 203);

    auto g_xzzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 204);

    auto g_xzzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 205);

    auto g_xzzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 206);

    auto g_xzzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 207);

    auto g_xzzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 208);

    auto g_xzzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 209);

    auto g_yyyy_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 210);

    auto g_yyyy_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 211);

    auto g_yyyy_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 212);

    auto g_yyyy_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 213);

    auto g_yyyy_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 214);

    auto g_yyyy_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 215);

    auto g_yyyy_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 216);

    auto g_yyyy_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 217);

    auto g_yyyy_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 218);

    auto g_yyyy_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 219);

    auto g_yyyy_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 220);

    auto g_yyyy_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 221);

    auto g_yyyy_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 222);

    auto g_yyyy_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 223);

    auto g_yyyy_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 224);

    auto g_yyyy_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 225);

    auto g_yyyy_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 226);

    auto g_yyyy_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 227);

    auto g_yyyy_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 228);

    auto g_yyyy_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 229);

    auto g_yyyy_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 230);

    auto g_yyyz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 232);

    auto g_yyyz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 234);

    auto g_yyyz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 237);

    auto g_yyyz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 241);

    auto g_yyyz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 246);

    auto g_yyzz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 252);

    auto g_yyzz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 253);

    auto g_yyzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 254);

    auto g_yyzz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 255);

    auto g_yyzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 256);

    auto g_yyzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 257);

    auto g_yyzz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 258);

    auto g_yyzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 259);

    auto g_yyzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 260);

    auto g_yyzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 261);

    auto g_yyzz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 262);

    auto g_yyzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 263);

    auto g_yyzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 264);

    auto g_yyzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 265);

    auto g_yyzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 266);

    auto g_yyzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 267);

    auto g_yyzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 268);

    auto g_yyzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 269);

    auto g_yyzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 270);

    auto g_yyzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 271);

    auto g_yyzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 272);

    auto g_yzzz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 273);

    auto g_yzzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 275);

    auto g_yzzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 277);

    auto g_yzzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 278);

    auto g_yzzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 280);

    auto g_yzzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 281);

    auto g_yzzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 282);

    auto g_yzzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 284);

    auto g_yzzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 285);

    auto g_yzzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 286);

    auto g_yzzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 287);

    auto g_yzzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 289);

    auto g_yzzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 290);

    auto g_yzzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 291);

    auto g_yzzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 292);

    auto g_yzzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 293);

    auto g_zzzz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 294);

    auto g_zzzz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 295);

    auto g_zzzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 296);

    auto g_zzzz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 297);

    auto g_zzzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 298);

    auto g_zzzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 299);

    auto g_zzzz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 300);

    auto g_zzzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 301);

    auto g_zzzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 302);

    auto g_zzzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 303);

    auto g_zzzz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 304);

    auto g_zzzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 305);

    auto g_zzzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 306);

    auto g_zzzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 307);

    auto g_zzzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 308);

    auto g_zzzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 309);

    auto g_zzzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 310);

    auto g_zzzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 311);

    auto g_zzzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 312);

    auto g_zzzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 313);

    auto g_zzzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 314);

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

    auto g_xxxxz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 32);

    auto g_xxxxz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 34);

    auto g_xxxxz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 35);

    auto g_xxxxz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 37);

    auto g_xxxxz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 38);

    auto g_xxxxz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 39);

    auto g_xxxxz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 41);

    auto g_xxxxz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 42);

    auto g_xxxxz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 43);

    auto g_xxxxz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 44);

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

    auto g_xyyzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 184);

    auto g_xyyzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 187);

    auto g_xyyzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 188);

    auto g_xyyzz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 191);

    auto g_xyyzz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 192);

    auto g_xyyzz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 193);

    auto g_xzzzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 212);

    auto g_xzzzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 214);

    auto g_xzzzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 215);

    auto g_xzzzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 217);

    auto g_xzzzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 218);

    auto g_xzzzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 219);

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

    auto g_yyyyz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 242);

    auto g_yyyyz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 244);

    auto g_yyyyz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 245);

    auto g_yyyyz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 247);

    auto g_yyyyz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 248);

    auto g_yyyyz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 249);

    auto g_yyyyz_yyyz_1 = pbuffer.data(idx_eri_1_hg + 251);

    auto g_yyyyz_yyzz_1 = pbuffer.data(idx_eri_1_hg + 252);

    auto g_yyyyz_yzzz_1 = pbuffer.data(idx_eri_1_hg + 253);

    auto g_yyyyz_zzzz_1 = pbuffer.data(idx_eri_1_hg + 254);

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

    auto g_yzzzz_xxxy_1 = pbuffer.data(idx_eri_1_hg + 286);

    auto g_yzzzz_xxxz_1 = pbuffer.data(idx_eri_1_hg + 287);

    auto g_yzzzz_xxyy_1 = pbuffer.data(idx_eri_1_hg + 288);

    auto g_yzzzz_xxyz_1 = pbuffer.data(idx_eri_1_hg + 289);

    auto g_yzzzz_xxzz_1 = pbuffer.data(idx_eri_1_hg + 290);

    auto g_yzzzz_xyyy_1 = pbuffer.data(idx_eri_1_hg + 291);

    auto g_yzzzz_xyyz_1 = pbuffer.data(idx_eri_1_hg + 292);

    auto g_yzzzz_xyzz_1 = pbuffer.data(idx_eri_1_hg + 293);

    auto g_yzzzz_xzzz_1 = pbuffer.data(idx_eri_1_hg + 294);

    auto g_yzzzz_yyyy_1 = pbuffer.data(idx_eri_1_hg + 295);

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

    // Set up components of auxiliary buffer : HH

    auto g_xxxxx_xxxxx_1 = pbuffer.data(idx_eri_1_hh);

    auto g_xxxxx_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 1);

    auto g_xxxxx_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 2);

    auto g_xxxxx_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 3);

    auto g_xxxxx_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 4);

    auto g_xxxxx_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 5);

    auto g_xxxxx_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 6);

    auto g_xxxxx_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 7);

    auto g_xxxxx_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 8);

    auto g_xxxxx_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 9);

    auto g_xxxxx_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 10);

    auto g_xxxxx_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 11);

    auto g_xxxxx_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 12);

    auto g_xxxxx_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 13);

    auto g_xxxxx_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 14);

    auto g_xxxxx_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 15);

    auto g_xxxxx_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 16);

    auto g_xxxxx_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 17);

    auto g_xxxxx_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 18);

    auto g_xxxxx_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 19);

    auto g_xxxxx_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 20);

    auto g_xxxxy_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 21);

    auto g_xxxxy_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 22);

    auto g_xxxxy_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 23);

    auto g_xxxxy_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 24);

    auto g_xxxxy_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 26);

    auto g_xxxxy_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 27);

    auto g_xxxxy_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 30);

    auto g_xxxxy_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 31);

    auto g_xxxxy_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 35);

    auto g_xxxxy_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 36);

    auto g_xxxxz_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 42);

    auto g_xxxxz_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 43);

    auto g_xxxxz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 44);

    auto g_xxxxz_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 45);

    auto g_xxxxz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 46);

    auto g_xxxxz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 47);

    auto g_xxxxz_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 48);

    auto g_xxxxz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 49);

    auto g_xxxxz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 50);

    auto g_xxxxz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 51);

    auto g_xxxxz_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 52);

    auto g_xxxxz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 53);

    auto g_xxxxz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 54);

    auto g_xxxxz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 55);

    auto g_xxxxz_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 56);

    auto g_xxxxz_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 58);

    auto g_xxxxz_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 59);

    auto g_xxxxz_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 60);

    auto g_xxxxz_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 61);

    auto g_xxxxz_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 62);

    auto g_xxxyy_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 63);

    auto g_xxxyy_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 64);

    auto g_xxxyy_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 65);

    auto g_xxxyy_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 66);

    auto g_xxxyy_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 67);

    auto g_xxxyy_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 68);

    auto g_xxxyy_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 69);

    auto g_xxxyy_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 70);

    auto g_xxxyy_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 71);

    auto g_xxxyy_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 72);

    auto g_xxxyy_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 73);

    auto g_xxxyy_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 74);

    auto g_xxxyy_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 75);

    auto g_xxxyy_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 76);

    auto g_xxxyy_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 77);

    auto g_xxxyy_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 78);

    auto g_xxxyy_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 79);

    auto g_xxxyy_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 80);

    auto g_xxxyy_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 81);

    auto g_xxxyy_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 82);

    auto g_xxxyy_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 83);

    auto g_xxxzz_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 105);

    auto g_xxxzz_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 106);

    auto g_xxxzz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 107);

    auto g_xxxzz_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 108);

    auto g_xxxzz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 109);

    auto g_xxxzz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 110);

    auto g_xxxzz_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 111);

    auto g_xxxzz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 112);

    auto g_xxxzz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 113);

    auto g_xxxzz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 114);

    auto g_xxxzz_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 115);

    auto g_xxxzz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 116);

    auto g_xxxzz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 117);

    auto g_xxxzz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 118);

    auto g_xxxzz_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 119);

    auto g_xxxzz_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 120);

    auto g_xxxzz_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 121);

    auto g_xxxzz_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 122);

    auto g_xxxzz_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 123);

    auto g_xxxzz_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 124);

    auto g_xxxzz_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 125);

    auto g_xxyyy_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 126);

    auto g_xxyyy_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 127);

    auto g_xxyyy_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 128);

    auto g_xxyyy_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 129);

    auto g_xxyyy_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 130);

    auto g_xxyyy_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 131);

    auto g_xxyyy_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 132);

    auto g_xxyyy_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 133);

    auto g_xxyyy_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 134);

    auto g_xxyyy_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 135);

    auto g_xxyyy_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 136);

    auto g_xxyyy_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 137);

    auto g_xxyyy_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 138);

    auto g_xxyyy_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 139);

    auto g_xxyyy_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 140);

    auto g_xxyyy_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 141);

    auto g_xxyyy_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 142);

    auto g_xxyyy_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 143);

    auto g_xxyyy_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 144);

    auto g_xxyyy_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 145);

    auto g_xxyyy_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 146);

    auto g_xxyyz_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 148);

    auto g_xxyyz_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 150);

    auto g_xxyyz_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 153);

    auto g_xxyyz_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 157);

    auto g_xxyzz_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 168);

    auto g_xxyzz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 170);

    auto g_xxyzz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 173);

    auto g_xxyzz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 177);

    auto g_xxyzz_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 182);

    auto g_xxzzz_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 189);

    auto g_xxzzz_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 190);

    auto g_xxzzz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 191);

    auto g_xxzzz_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 192);

    auto g_xxzzz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 193);

    auto g_xxzzz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 194);

    auto g_xxzzz_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 195);

    auto g_xxzzz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 196);

    auto g_xxzzz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 197);

    auto g_xxzzz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 198);

    auto g_xxzzz_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 199);

    auto g_xxzzz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 200);

    auto g_xxzzz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 201);

    auto g_xxzzz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 202);

    auto g_xxzzz_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 203);

    auto g_xxzzz_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 204);

    auto g_xxzzz_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 205);

    auto g_xxzzz_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 206);

    auto g_xxzzz_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 207);

    auto g_xxzzz_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 208);

    auto g_xxzzz_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 209);

    auto g_xyyyy_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 210);

    auto g_xyyyy_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 211);

    auto g_xyyyy_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 213);

    auto g_xyyyy_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 214);

    auto g_xyyyy_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 216);

    auto g_xyyyy_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 217);

    auto g_xyyyy_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 218);

    auto g_xyyyy_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 220);

    auto g_xyyyy_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 221);

    auto g_xyyyy_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 222);

    auto g_xyyyy_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 223);

    auto g_xyyyy_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 225);

    auto g_xyyyy_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 226);

    auto g_xyyyy_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 227);

    auto g_xyyyy_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 228);

    auto g_xyyyy_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 229);

    auto g_xyyyy_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 230);

    auto g_xyyzz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 256);

    auto g_xyyzz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 259);

    auto g_xyyzz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 260);

    auto g_xyyzz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 263);

    auto g_xyyzz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 264);

    auto g_xyyzz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 265);

    auto g_xyyzz_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 267);

    auto g_xyyzz_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 268);

    auto g_xyyzz_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 269);

    auto g_xyyzz_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 270);

    auto g_xyyzz_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 271);

    auto g_xyyzz_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 272);

    auto g_xzzzz_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 294);

    auto g_xzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 296);

    auto g_xzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 298);

    auto g_xzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 299);

    auto g_xzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 301);

    auto g_xzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 302);

    auto g_xzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 303);

    auto g_xzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 305);

    auto g_xzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 306);

    auto g_xzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 307);

    auto g_xzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 308);

    auto g_xzzzz_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 309);

    auto g_xzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 310);

    auto g_xzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 311);

    auto g_xzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 312);

    auto g_xzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 313);

    auto g_xzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 314);

    auto g_yyyyy_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 315);

    auto g_yyyyy_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 316);

    auto g_yyyyy_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 317);

    auto g_yyyyy_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 318);

    auto g_yyyyy_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 319);

    auto g_yyyyy_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 320);

    auto g_yyyyy_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 321);

    auto g_yyyyy_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 322);

    auto g_yyyyy_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 323);

    auto g_yyyyy_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 324);

    auto g_yyyyy_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 325);

    auto g_yyyyy_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 326);

    auto g_yyyyy_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 327);

    auto g_yyyyy_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 328);

    auto g_yyyyy_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 329);

    auto g_yyyyy_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 330);

    auto g_yyyyy_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 331);

    auto g_yyyyy_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 332);

    auto g_yyyyy_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 333);

    auto g_yyyyy_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 334);

    auto g_yyyyy_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 335);

    auto g_yyyyz_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 337);

    auto g_yyyyz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 338);

    auto g_yyyyz_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 339);

    auto g_yyyyz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 340);

    auto g_yyyyz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 341);

    auto g_yyyyz_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 342);

    auto g_yyyyz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 343);

    auto g_yyyyz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 344);

    auto g_yyyyz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 345);

    auto g_yyyyz_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 346);

    auto g_yyyyz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 347);

    auto g_yyyyz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 348);

    auto g_yyyyz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 349);

    auto g_yyyyz_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 350);

    auto g_yyyyz_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 351);

    auto g_yyyyz_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 352);

    auto g_yyyyz_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 353);

    auto g_yyyyz_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 354);

    auto g_yyyyz_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 355);

    auto g_yyyyz_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 356);

    auto g_yyyzz_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 357);

    auto g_yyyzz_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 358);

    auto g_yyyzz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 359);

    auto g_yyyzz_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 360);

    auto g_yyyzz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 361);

    auto g_yyyzz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 362);

    auto g_yyyzz_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 363);

    auto g_yyyzz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 364);

    auto g_yyyzz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 365);

    auto g_yyyzz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 366);

    auto g_yyyzz_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 367);

    auto g_yyyzz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 368);

    auto g_yyyzz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 369);

    auto g_yyyzz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 370);

    auto g_yyyzz_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 371);

    auto g_yyyzz_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 372);

    auto g_yyyzz_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 373);

    auto g_yyyzz_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 374);

    auto g_yyyzz_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 375);

    auto g_yyyzz_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 376);

    auto g_yyyzz_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 377);

    auto g_yyzzz_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 378);

    auto g_yyzzz_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 379);

    auto g_yyzzz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 380);

    auto g_yyzzz_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 381);

    auto g_yyzzz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 382);

    auto g_yyzzz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 383);

    auto g_yyzzz_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 384);

    auto g_yyzzz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 385);

    auto g_yyzzz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 386);

    auto g_yyzzz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 387);

    auto g_yyzzz_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 388);

    auto g_yyzzz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 389);

    auto g_yyzzz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 390);

    auto g_yyzzz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 391);

    auto g_yyzzz_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 392);

    auto g_yyzzz_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 393);

    auto g_yyzzz_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 394);

    auto g_yyzzz_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 395);

    auto g_yyzzz_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 396);

    auto g_yyzzz_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 397);

    auto g_yyzzz_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 398);

    auto g_yzzzz_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 399);

    auto g_yzzzz_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 400);

    auto g_yzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 401);

    auto g_yzzzz_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 402);

    auto g_yzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 403);

    auto g_yzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 404);

    auto g_yzzzz_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 405);

    auto g_yzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 406);

    auto g_yzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 407);

    auto g_yzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 408);

    auto g_yzzzz_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 409);

    auto g_yzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 410);

    auto g_yzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 411);

    auto g_yzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 412);

    auto g_yzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 413);

    auto g_yzzzz_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 414);

    auto g_yzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 415);

    auto g_yzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 416);

    auto g_yzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 417);

    auto g_yzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 418);

    auto g_yzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 419);

    auto g_zzzzz_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 420);

    auto g_zzzzz_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 421);

    auto g_zzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 422);

    auto g_zzzzz_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 423);

    auto g_zzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 424);

    auto g_zzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 425);

    auto g_zzzzz_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 426);

    auto g_zzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 427);

    auto g_zzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 428);

    auto g_zzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 429);

    auto g_zzzzz_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 430);

    auto g_zzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 431);

    auto g_zzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 432);

    auto g_zzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 433);

    auto g_zzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 434);

    auto g_zzzzz_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 435);

    auto g_zzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 436);

    auto g_zzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 437);

    auto g_zzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 438);

    auto g_zzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 439);

    auto g_zzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_hh + 440);

    // Set up 0-21 components of targeted buffer : IH

    auto g_xxxxxx_xxxxx_0 = pbuffer.data(idx_eri_0_ih);

    auto g_xxxxxx_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 1);

    auto g_xxxxxx_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 2);

    auto g_xxxxxx_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 3);

    auto g_xxxxxx_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 4);

    auto g_xxxxxx_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 5);

    auto g_xxxxxx_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 6);

    auto g_xxxxxx_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 7);

    auto g_xxxxxx_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 8);

    auto g_xxxxxx_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 9);

    auto g_xxxxxx_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 10);

    auto g_xxxxxx_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 11);

    auto g_xxxxxx_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 12);

    auto g_xxxxxx_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 13);

    auto g_xxxxxx_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 14);

    auto g_xxxxxx_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 15);

    auto g_xxxxxx_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 16);

    auto g_xxxxxx_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 17);

    auto g_xxxxxx_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 18);

    auto g_xxxxxx_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 19);

    auto g_xxxxxx_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 20);

    #pragma omp simd aligned(g_xxxx_xxxxx_0, g_xxxx_xxxxx_1, g_xxxx_xxxxy_0, g_xxxx_xxxxy_1, g_xxxx_xxxxz_0, g_xxxx_xxxxz_1, g_xxxx_xxxyy_0, g_xxxx_xxxyy_1, g_xxxx_xxxyz_0, g_xxxx_xxxyz_1, g_xxxx_xxxzz_0, g_xxxx_xxxzz_1, g_xxxx_xxyyy_0, g_xxxx_xxyyy_1, g_xxxx_xxyyz_0, g_xxxx_xxyyz_1, g_xxxx_xxyzz_0, g_xxxx_xxyzz_1, g_xxxx_xxzzz_0, g_xxxx_xxzzz_1, g_xxxx_xyyyy_0, g_xxxx_xyyyy_1, g_xxxx_xyyyz_0, g_xxxx_xyyyz_1, g_xxxx_xyyzz_0, g_xxxx_xyyzz_1, g_xxxx_xyzzz_0, g_xxxx_xyzzz_1, g_xxxx_xzzzz_0, g_xxxx_xzzzz_1, g_xxxx_yyyyy_0, g_xxxx_yyyyy_1, g_xxxx_yyyyz_0, g_xxxx_yyyyz_1, g_xxxx_yyyzz_0, g_xxxx_yyyzz_1, g_xxxx_yyzzz_0, g_xxxx_yyzzz_1, g_xxxx_yzzzz_0, g_xxxx_yzzzz_1, g_xxxx_zzzzz_0, g_xxxx_zzzzz_1, g_xxxxx_xxxx_1, g_xxxxx_xxxxx_1, g_xxxxx_xxxxy_1, g_xxxxx_xxxxz_1, g_xxxxx_xxxy_1, g_xxxxx_xxxyy_1, g_xxxxx_xxxyz_1, g_xxxxx_xxxz_1, g_xxxxx_xxxzz_1, g_xxxxx_xxyy_1, g_xxxxx_xxyyy_1, g_xxxxx_xxyyz_1, g_xxxxx_xxyz_1, g_xxxxx_xxyzz_1, g_xxxxx_xxzz_1, g_xxxxx_xxzzz_1, g_xxxxx_xyyy_1, g_xxxxx_xyyyy_1, g_xxxxx_xyyyz_1, g_xxxxx_xyyz_1, g_xxxxx_xyyzz_1, g_xxxxx_xyzz_1, g_xxxxx_xyzzz_1, g_xxxxx_xzzz_1, g_xxxxx_xzzzz_1, g_xxxxx_yyyy_1, g_xxxxx_yyyyy_1, g_xxxxx_yyyyz_1, g_xxxxx_yyyz_1, g_xxxxx_yyyzz_1, g_xxxxx_yyzz_1, g_xxxxx_yyzzz_1, g_xxxxx_yzzz_1, g_xxxxx_yzzzz_1, g_xxxxx_zzzz_1, g_xxxxx_zzzzz_1, g_xxxxxx_xxxxx_0, g_xxxxxx_xxxxy_0, g_xxxxxx_xxxxz_0, g_xxxxxx_xxxyy_0, g_xxxxxx_xxxyz_0, g_xxxxxx_xxxzz_0, g_xxxxxx_xxyyy_0, g_xxxxxx_xxyyz_0, g_xxxxxx_xxyzz_0, g_xxxxxx_xxzzz_0, g_xxxxxx_xyyyy_0, g_xxxxxx_xyyyz_0, g_xxxxxx_xyyzz_0, g_xxxxxx_xyzzz_0, g_xxxxxx_xzzzz_0, g_xxxxxx_yyyyy_0, g_xxxxxx_yyyyz_0, g_xxxxxx_yyyzz_0, g_xxxxxx_yyzzz_0, g_xxxxxx_yzzzz_0, g_xxxxxx_zzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxx_xxxxx_0[i] = 5.0 * g_xxxx_xxxxx_0[i] * fbe_0 - 5.0 * g_xxxx_xxxxx_1[i] * fz_be_0 + 5.0 * g_xxxxx_xxxx_1[i] * fe_0 + g_xxxxx_xxxxx_1[i] * pa_x[i];

        g_xxxxxx_xxxxy_0[i] = 5.0 * g_xxxx_xxxxy_0[i] * fbe_0 - 5.0 * g_xxxx_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxxxx_xxxy_1[i] * fe_0 + g_xxxxx_xxxxy_1[i] * pa_x[i];

        g_xxxxxx_xxxxz_0[i] = 5.0 * g_xxxx_xxxxz_0[i] * fbe_0 - 5.0 * g_xxxx_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxxxx_xxxz_1[i] * fe_0 + g_xxxxx_xxxxz_1[i] * pa_x[i];

        g_xxxxxx_xxxyy_0[i] = 5.0 * g_xxxx_xxxyy_0[i] * fbe_0 - 5.0 * g_xxxx_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxxxx_xxyy_1[i] * fe_0 + g_xxxxx_xxxyy_1[i] * pa_x[i];

        g_xxxxxx_xxxyz_0[i] = 5.0 * g_xxxx_xxxyz_0[i] * fbe_0 - 5.0 * g_xxxx_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxxx_xxyz_1[i] * fe_0 + g_xxxxx_xxxyz_1[i] * pa_x[i];

        g_xxxxxx_xxxzz_0[i] = 5.0 * g_xxxx_xxxzz_0[i] * fbe_0 - 5.0 * g_xxxx_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_xxzz_1[i] * fe_0 + g_xxxxx_xxxzz_1[i] * pa_x[i];

        g_xxxxxx_xxyyy_0[i] = 5.0 * g_xxxx_xxyyy_0[i] * fbe_0 - 5.0 * g_xxxx_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxxxx_xyyy_1[i] * fe_0 + g_xxxxx_xxyyy_1[i] * pa_x[i];

        g_xxxxxx_xxyyz_0[i] = 5.0 * g_xxxx_xxyyz_0[i] * fbe_0 - 5.0 * g_xxxx_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxxx_xyyz_1[i] * fe_0 + g_xxxxx_xxyyz_1[i] * pa_x[i];

        g_xxxxxx_xxyzz_0[i] = 5.0 * g_xxxx_xxyzz_0[i] * fbe_0 - 5.0 * g_xxxx_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_xyzz_1[i] * fe_0 + g_xxxxx_xxyzz_1[i] * pa_x[i];

        g_xxxxxx_xxzzz_0[i] = 5.0 * g_xxxx_xxzzz_0[i] * fbe_0 - 5.0 * g_xxxx_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_xzzz_1[i] * fe_0 + g_xxxxx_xxzzz_1[i] * pa_x[i];

        g_xxxxxx_xyyyy_0[i] = 5.0 * g_xxxx_xyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_xyyyy_1[i] * fz_be_0 + g_xxxxx_yyyy_1[i] * fe_0 + g_xxxxx_xyyyy_1[i] * pa_x[i];

        g_xxxxxx_xyyyz_0[i] = 5.0 * g_xxxx_xyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_xyyyz_1[i] * fz_be_0 + g_xxxxx_yyyz_1[i] * fe_0 + g_xxxxx_xyyyz_1[i] * pa_x[i];

        g_xxxxxx_xyyzz_0[i] = 5.0 * g_xxxx_xyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_xyyzz_1[i] * fz_be_0 + g_xxxxx_yyzz_1[i] * fe_0 + g_xxxxx_xyyzz_1[i] * pa_x[i];

        g_xxxxxx_xyzzz_0[i] = 5.0 * g_xxxx_xyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_xyzzz_1[i] * fz_be_0 + g_xxxxx_yzzz_1[i] * fe_0 + g_xxxxx_xyzzz_1[i] * pa_x[i];

        g_xxxxxx_xzzzz_0[i] = 5.0 * g_xxxx_xzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_xzzzz_1[i] * fz_be_0 + g_xxxxx_zzzz_1[i] * fe_0 + g_xxxxx_xzzzz_1[i] * pa_x[i];

        g_xxxxxx_yyyyy_0[i] = 5.0 * g_xxxx_yyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_yyyyy_1[i] * fz_be_0 + g_xxxxx_yyyyy_1[i] * pa_x[i];

        g_xxxxxx_yyyyz_0[i] = 5.0 * g_xxxx_yyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_yyyyz_1[i] * fz_be_0 + g_xxxxx_yyyyz_1[i] * pa_x[i];

        g_xxxxxx_yyyzz_0[i] = 5.0 * g_xxxx_yyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_yyyzz_1[i] * fz_be_0 + g_xxxxx_yyyzz_1[i] * pa_x[i];

        g_xxxxxx_yyzzz_0[i] = 5.0 * g_xxxx_yyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_yyzzz_1[i] * fz_be_0 + g_xxxxx_yyzzz_1[i] * pa_x[i];

        g_xxxxxx_yzzzz_0[i] = 5.0 * g_xxxx_yzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_yzzzz_1[i] * fz_be_0 + g_xxxxx_yzzzz_1[i] * pa_x[i];

        g_xxxxxx_zzzzz_0[i] = 5.0 * g_xxxx_zzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_zzzzz_1[i] * fz_be_0 + g_xxxxx_zzzzz_1[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : IH

    auto g_xxxxxy_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 21);

    auto g_xxxxxy_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 22);

    auto g_xxxxxy_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 23);

    auto g_xxxxxy_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 24);

    auto g_xxxxxy_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 25);

    auto g_xxxxxy_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 26);

    auto g_xxxxxy_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 27);

    auto g_xxxxxy_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 28);

    auto g_xxxxxy_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 29);

    auto g_xxxxxy_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 30);

    auto g_xxxxxy_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 31);

    auto g_xxxxxy_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 32);

    auto g_xxxxxy_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 33);

    auto g_xxxxxy_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 34);

    auto g_xxxxxy_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 35);

    auto g_xxxxxy_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 36);

    auto g_xxxxxy_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 37);

    auto g_xxxxxy_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 38);

    auto g_xxxxxy_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 39);

    auto g_xxxxxy_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 40);

    auto g_xxxxxy_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 41);

    #pragma omp simd aligned(g_xxxxx_xxxx_1, g_xxxxx_xxxxx_1, g_xxxxx_xxxxy_1, g_xxxxx_xxxxz_1, g_xxxxx_xxxy_1, g_xxxxx_xxxyy_1, g_xxxxx_xxxyz_1, g_xxxxx_xxxz_1, g_xxxxx_xxxzz_1, g_xxxxx_xxyy_1, g_xxxxx_xxyyy_1, g_xxxxx_xxyyz_1, g_xxxxx_xxyz_1, g_xxxxx_xxyzz_1, g_xxxxx_xxzz_1, g_xxxxx_xxzzz_1, g_xxxxx_xyyy_1, g_xxxxx_xyyyy_1, g_xxxxx_xyyyz_1, g_xxxxx_xyyz_1, g_xxxxx_xyyzz_1, g_xxxxx_xyzz_1, g_xxxxx_xyzzz_1, g_xxxxx_xzzz_1, g_xxxxx_xzzzz_1, g_xxxxx_yyyy_1, g_xxxxx_yyyyy_1, g_xxxxx_yyyyz_1, g_xxxxx_yyyz_1, g_xxxxx_yyyzz_1, g_xxxxx_yyzz_1, g_xxxxx_yyzzz_1, g_xxxxx_yzzz_1, g_xxxxx_yzzzz_1, g_xxxxx_zzzz_1, g_xxxxx_zzzzz_1, g_xxxxxy_xxxxx_0, g_xxxxxy_xxxxy_0, g_xxxxxy_xxxxz_0, g_xxxxxy_xxxyy_0, g_xxxxxy_xxxyz_0, g_xxxxxy_xxxzz_0, g_xxxxxy_xxyyy_0, g_xxxxxy_xxyyz_0, g_xxxxxy_xxyzz_0, g_xxxxxy_xxzzz_0, g_xxxxxy_xyyyy_0, g_xxxxxy_xyyyz_0, g_xxxxxy_xyyzz_0, g_xxxxxy_xyzzz_0, g_xxxxxy_xzzzz_0, g_xxxxxy_yyyyy_0, g_xxxxxy_yyyyz_0, g_xxxxxy_yyyzz_0, g_xxxxxy_yyzzz_0, g_xxxxxy_yzzzz_0, g_xxxxxy_zzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxy_xxxxx_0[i] = g_xxxxx_xxxxx_1[i] * pa_y[i];

        g_xxxxxy_xxxxy_0[i] = g_xxxxx_xxxx_1[i] * fe_0 + g_xxxxx_xxxxy_1[i] * pa_y[i];

        g_xxxxxy_xxxxz_0[i] = g_xxxxx_xxxxz_1[i] * pa_y[i];

        g_xxxxxy_xxxyy_0[i] = 2.0 * g_xxxxx_xxxy_1[i] * fe_0 + g_xxxxx_xxxyy_1[i] * pa_y[i];

        g_xxxxxy_xxxyz_0[i] = g_xxxxx_xxxz_1[i] * fe_0 + g_xxxxx_xxxyz_1[i] * pa_y[i];

        g_xxxxxy_xxxzz_0[i] = g_xxxxx_xxxzz_1[i] * pa_y[i];

        g_xxxxxy_xxyyy_0[i] = 3.0 * g_xxxxx_xxyy_1[i] * fe_0 + g_xxxxx_xxyyy_1[i] * pa_y[i];

        g_xxxxxy_xxyyz_0[i] = 2.0 * g_xxxxx_xxyz_1[i] * fe_0 + g_xxxxx_xxyyz_1[i] * pa_y[i];

        g_xxxxxy_xxyzz_0[i] = g_xxxxx_xxzz_1[i] * fe_0 + g_xxxxx_xxyzz_1[i] * pa_y[i];

        g_xxxxxy_xxzzz_0[i] = g_xxxxx_xxzzz_1[i] * pa_y[i];

        g_xxxxxy_xyyyy_0[i] = 4.0 * g_xxxxx_xyyy_1[i] * fe_0 + g_xxxxx_xyyyy_1[i] * pa_y[i];

        g_xxxxxy_xyyyz_0[i] = 3.0 * g_xxxxx_xyyz_1[i] * fe_0 + g_xxxxx_xyyyz_1[i] * pa_y[i];

        g_xxxxxy_xyyzz_0[i] = 2.0 * g_xxxxx_xyzz_1[i] * fe_0 + g_xxxxx_xyyzz_1[i] * pa_y[i];

        g_xxxxxy_xyzzz_0[i] = g_xxxxx_xzzz_1[i] * fe_0 + g_xxxxx_xyzzz_1[i] * pa_y[i];

        g_xxxxxy_xzzzz_0[i] = g_xxxxx_xzzzz_1[i] * pa_y[i];

        g_xxxxxy_yyyyy_0[i] = 5.0 * g_xxxxx_yyyy_1[i] * fe_0 + g_xxxxx_yyyyy_1[i] * pa_y[i];

        g_xxxxxy_yyyyz_0[i] = 4.0 * g_xxxxx_yyyz_1[i] * fe_0 + g_xxxxx_yyyyz_1[i] * pa_y[i];

        g_xxxxxy_yyyzz_0[i] = 3.0 * g_xxxxx_yyzz_1[i] * fe_0 + g_xxxxx_yyyzz_1[i] * pa_y[i];

        g_xxxxxy_yyzzz_0[i] = 2.0 * g_xxxxx_yzzz_1[i] * fe_0 + g_xxxxx_yyzzz_1[i] * pa_y[i];

        g_xxxxxy_yzzzz_0[i] = g_xxxxx_zzzz_1[i] * fe_0 + g_xxxxx_yzzzz_1[i] * pa_y[i];

        g_xxxxxy_zzzzz_0[i] = g_xxxxx_zzzzz_1[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : IH

    auto g_xxxxxz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 42);

    auto g_xxxxxz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 43);

    auto g_xxxxxz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 44);

    auto g_xxxxxz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 45);

    auto g_xxxxxz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 46);

    auto g_xxxxxz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 47);

    auto g_xxxxxz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 48);

    auto g_xxxxxz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 49);

    auto g_xxxxxz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 50);

    auto g_xxxxxz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 51);

    auto g_xxxxxz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 52);

    auto g_xxxxxz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 53);

    auto g_xxxxxz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 54);

    auto g_xxxxxz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 55);

    auto g_xxxxxz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 56);

    auto g_xxxxxz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 57);

    auto g_xxxxxz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 58);

    auto g_xxxxxz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 59);

    auto g_xxxxxz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 60);

    auto g_xxxxxz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 61);

    auto g_xxxxxz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 62);

    #pragma omp simd aligned(g_xxxxx_xxxx_1, g_xxxxx_xxxxx_1, g_xxxxx_xxxxy_1, g_xxxxx_xxxxz_1, g_xxxxx_xxxy_1, g_xxxxx_xxxyy_1, g_xxxxx_xxxyz_1, g_xxxxx_xxxz_1, g_xxxxx_xxxzz_1, g_xxxxx_xxyy_1, g_xxxxx_xxyyy_1, g_xxxxx_xxyyz_1, g_xxxxx_xxyz_1, g_xxxxx_xxyzz_1, g_xxxxx_xxzz_1, g_xxxxx_xxzzz_1, g_xxxxx_xyyy_1, g_xxxxx_xyyyy_1, g_xxxxx_xyyyz_1, g_xxxxx_xyyz_1, g_xxxxx_xyyzz_1, g_xxxxx_xyzz_1, g_xxxxx_xyzzz_1, g_xxxxx_xzzz_1, g_xxxxx_xzzzz_1, g_xxxxx_yyyy_1, g_xxxxx_yyyyy_1, g_xxxxx_yyyyz_1, g_xxxxx_yyyz_1, g_xxxxx_yyyzz_1, g_xxxxx_yyzz_1, g_xxxxx_yyzzz_1, g_xxxxx_yzzz_1, g_xxxxx_yzzzz_1, g_xxxxx_zzzz_1, g_xxxxx_zzzzz_1, g_xxxxxz_xxxxx_0, g_xxxxxz_xxxxy_0, g_xxxxxz_xxxxz_0, g_xxxxxz_xxxyy_0, g_xxxxxz_xxxyz_0, g_xxxxxz_xxxzz_0, g_xxxxxz_xxyyy_0, g_xxxxxz_xxyyz_0, g_xxxxxz_xxyzz_0, g_xxxxxz_xxzzz_0, g_xxxxxz_xyyyy_0, g_xxxxxz_xyyyz_0, g_xxxxxz_xyyzz_0, g_xxxxxz_xyzzz_0, g_xxxxxz_xzzzz_0, g_xxxxxz_yyyyy_0, g_xxxxxz_yyyyz_0, g_xxxxxz_yyyzz_0, g_xxxxxz_yyzzz_0, g_xxxxxz_yzzzz_0, g_xxxxxz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxz_xxxxx_0[i] = g_xxxxx_xxxxx_1[i] * pa_z[i];

        g_xxxxxz_xxxxy_0[i] = g_xxxxx_xxxxy_1[i] * pa_z[i];

        g_xxxxxz_xxxxz_0[i] = g_xxxxx_xxxx_1[i] * fe_0 + g_xxxxx_xxxxz_1[i] * pa_z[i];

        g_xxxxxz_xxxyy_0[i] = g_xxxxx_xxxyy_1[i] * pa_z[i];

        g_xxxxxz_xxxyz_0[i] = g_xxxxx_xxxy_1[i] * fe_0 + g_xxxxx_xxxyz_1[i] * pa_z[i];

        g_xxxxxz_xxxzz_0[i] = 2.0 * g_xxxxx_xxxz_1[i] * fe_0 + g_xxxxx_xxxzz_1[i] * pa_z[i];

        g_xxxxxz_xxyyy_0[i] = g_xxxxx_xxyyy_1[i] * pa_z[i];

        g_xxxxxz_xxyyz_0[i] = g_xxxxx_xxyy_1[i] * fe_0 + g_xxxxx_xxyyz_1[i] * pa_z[i];

        g_xxxxxz_xxyzz_0[i] = 2.0 * g_xxxxx_xxyz_1[i] * fe_0 + g_xxxxx_xxyzz_1[i] * pa_z[i];

        g_xxxxxz_xxzzz_0[i] = 3.0 * g_xxxxx_xxzz_1[i] * fe_0 + g_xxxxx_xxzzz_1[i] * pa_z[i];

        g_xxxxxz_xyyyy_0[i] = g_xxxxx_xyyyy_1[i] * pa_z[i];

        g_xxxxxz_xyyyz_0[i] = g_xxxxx_xyyy_1[i] * fe_0 + g_xxxxx_xyyyz_1[i] * pa_z[i];

        g_xxxxxz_xyyzz_0[i] = 2.0 * g_xxxxx_xyyz_1[i] * fe_0 + g_xxxxx_xyyzz_1[i] * pa_z[i];

        g_xxxxxz_xyzzz_0[i] = 3.0 * g_xxxxx_xyzz_1[i] * fe_0 + g_xxxxx_xyzzz_1[i] * pa_z[i];

        g_xxxxxz_xzzzz_0[i] = 4.0 * g_xxxxx_xzzz_1[i] * fe_0 + g_xxxxx_xzzzz_1[i] * pa_z[i];

        g_xxxxxz_yyyyy_0[i] = g_xxxxx_yyyyy_1[i] * pa_z[i];

        g_xxxxxz_yyyyz_0[i] = g_xxxxx_yyyy_1[i] * fe_0 + g_xxxxx_yyyyz_1[i] * pa_z[i];

        g_xxxxxz_yyyzz_0[i] = 2.0 * g_xxxxx_yyyz_1[i] * fe_0 + g_xxxxx_yyyzz_1[i] * pa_z[i];

        g_xxxxxz_yyzzz_0[i] = 3.0 * g_xxxxx_yyzz_1[i] * fe_0 + g_xxxxx_yyzzz_1[i] * pa_z[i];

        g_xxxxxz_yzzzz_0[i] = 4.0 * g_xxxxx_yzzz_1[i] * fe_0 + g_xxxxx_yzzzz_1[i] * pa_z[i];

        g_xxxxxz_zzzzz_0[i] = 5.0 * g_xxxxx_zzzz_1[i] * fe_0 + g_xxxxx_zzzzz_1[i] * pa_z[i];
    }

    // Set up 63-84 components of targeted buffer : IH

    auto g_xxxxyy_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 63);

    auto g_xxxxyy_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 64);

    auto g_xxxxyy_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 65);

    auto g_xxxxyy_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 66);

    auto g_xxxxyy_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 67);

    auto g_xxxxyy_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 68);

    auto g_xxxxyy_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 69);

    auto g_xxxxyy_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 70);

    auto g_xxxxyy_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 71);

    auto g_xxxxyy_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 72);

    auto g_xxxxyy_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 73);

    auto g_xxxxyy_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 74);

    auto g_xxxxyy_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 75);

    auto g_xxxxyy_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 76);

    auto g_xxxxyy_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 77);

    auto g_xxxxyy_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 78);

    auto g_xxxxyy_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 79);

    auto g_xxxxyy_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 80);

    auto g_xxxxyy_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 81);

    auto g_xxxxyy_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 82);

    auto g_xxxxyy_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 83);

    #pragma omp simd aligned(g_xxxx_xxxxx_0, g_xxxx_xxxxx_1, g_xxxx_xxxxz_0, g_xxxx_xxxxz_1, g_xxxx_xxxzz_0, g_xxxx_xxxzz_1, g_xxxx_xxzzz_0, g_xxxx_xxzzz_1, g_xxxx_xzzzz_0, g_xxxx_xzzzz_1, g_xxxxy_xxxxx_1, g_xxxxy_xxxxz_1, g_xxxxy_xxxzz_1, g_xxxxy_xxzzz_1, g_xxxxy_xzzzz_1, g_xxxxyy_xxxxx_0, g_xxxxyy_xxxxy_0, g_xxxxyy_xxxxz_0, g_xxxxyy_xxxyy_0, g_xxxxyy_xxxyz_0, g_xxxxyy_xxxzz_0, g_xxxxyy_xxyyy_0, g_xxxxyy_xxyyz_0, g_xxxxyy_xxyzz_0, g_xxxxyy_xxzzz_0, g_xxxxyy_xyyyy_0, g_xxxxyy_xyyyz_0, g_xxxxyy_xyyzz_0, g_xxxxyy_xyzzz_0, g_xxxxyy_xzzzz_0, g_xxxxyy_yyyyy_0, g_xxxxyy_yyyyz_0, g_xxxxyy_yyyzz_0, g_xxxxyy_yyzzz_0, g_xxxxyy_yzzzz_0, g_xxxxyy_zzzzz_0, g_xxxyy_xxxxy_1, g_xxxyy_xxxy_1, g_xxxyy_xxxyy_1, g_xxxyy_xxxyz_1, g_xxxyy_xxyy_1, g_xxxyy_xxyyy_1, g_xxxyy_xxyyz_1, g_xxxyy_xxyz_1, g_xxxyy_xxyzz_1, g_xxxyy_xyyy_1, g_xxxyy_xyyyy_1, g_xxxyy_xyyyz_1, g_xxxyy_xyyz_1, g_xxxyy_xyyzz_1, g_xxxyy_xyzz_1, g_xxxyy_xyzzz_1, g_xxxyy_yyyy_1, g_xxxyy_yyyyy_1, g_xxxyy_yyyyz_1, g_xxxyy_yyyz_1, g_xxxyy_yyyzz_1, g_xxxyy_yyzz_1, g_xxxyy_yyzzz_1, g_xxxyy_yzzz_1, g_xxxyy_yzzzz_1, g_xxxyy_zzzzz_1, g_xxyy_xxxxy_0, g_xxyy_xxxxy_1, g_xxyy_xxxyy_0, g_xxyy_xxxyy_1, g_xxyy_xxxyz_0, g_xxyy_xxxyz_1, g_xxyy_xxyyy_0, g_xxyy_xxyyy_1, g_xxyy_xxyyz_0, g_xxyy_xxyyz_1, g_xxyy_xxyzz_0, g_xxyy_xxyzz_1, g_xxyy_xyyyy_0, g_xxyy_xyyyy_1, g_xxyy_xyyyz_0, g_xxyy_xyyyz_1, g_xxyy_xyyzz_0, g_xxyy_xyyzz_1, g_xxyy_xyzzz_0, g_xxyy_xyzzz_1, g_xxyy_yyyyy_0, g_xxyy_yyyyy_1, g_xxyy_yyyyz_0, g_xxyy_yyyyz_1, g_xxyy_yyyzz_0, g_xxyy_yyyzz_1, g_xxyy_yyzzz_0, g_xxyy_yyzzz_1, g_xxyy_yzzzz_0, g_xxyy_yzzzz_1, g_xxyy_zzzzz_0, g_xxyy_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxyy_xxxxx_0[i] = g_xxxx_xxxxx_0[i] * fbe_0 - g_xxxx_xxxxx_1[i] * fz_be_0 + g_xxxxy_xxxxx_1[i] * pa_y[i];

        g_xxxxyy_xxxxy_0[i] = 3.0 * g_xxyy_xxxxy_0[i] * fbe_0 - 3.0 * g_xxyy_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxxyy_xxxy_1[i] * fe_0 + g_xxxyy_xxxxy_1[i] * pa_x[i];

        g_xxxxyy_xxxxz_0[i] = g_xxxx_xxxxz_0[i] * fbe_0 - g_xxxx_xxxxz_1[i] * fz_be_0 + g_xxxxy_xxxxz_1[i] * pa_y[i];

        g_xxxxyy_xxxyy_0[i] = 3.0 * g_xxyy_xxxyy_0[i] * fbe_0 - 3.0 * g_xxyy_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxxyy_xxyy_1[i] * fe_0 + g_xxxyy_xxxyy_1[i] * pa_x[i];

        g_xxxxyy_xxxyz_0[i] = 3.0 * g_xxyy_xxxyz_0[i] * fbe_0 - 3.0 * g_xxyy_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxyy_xxyz_1[i] * fe_0 + g_xxxyy_xxxyz_1[i] * pa_x[i];

        g_xxxxyy_xxxzz_0[i] = g_xxxx_xxxzz_0[i] * fbe_0 - g_xxxx_xxxzz_1[i] * fz_be_0 + g_xxxxy_xxxzz_1[i] * pa_y[i];

        g_xxxxyy_xxyyy_0[i] = 3.0 * g_xxyy_xxyyy_0[i] * fbe_0 - 3.0 * g_xxyy_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxxyy_xyyy_1[i] * fe_0 + g_xxxyy_xxyyy_1[i] * pa_x[i];

        g_xxxxyy_xxyyz_0[i] = 3.0 * g_xxyy_xxyyz_0[i] * fbe_0 - 3.0 * g_xxyy_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxyy_xyyz_1[i] * fe_0 + g_xxxyy_xxyyz_1[i] * pa_x[i];

        g_xxxxyy_xxyzz_0[i] = 3.0 * g_xxyy_xxyzz_0[i] * fbe_0 - 3.0 * g_xxyy_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_xyzz_1[i] * fe_0 + g_xxxyy_xxyzz_1[i] * pa_x[i];

        g_xxxxyy_xxzzz_0[i] = g_xxxx_xxzzz_0[i] * fbe_0 - g_xxxx_xxzzz_1[i] * fz_be_0 + g_xxxxy_xxzzz_1[i] * pa_y[i];

        g_xxxxyy_xyyyy_0[i] = 3.0 * g_xxyy_xyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_xyyyy_1[i] * fz_be_0 + g_xxxyy_yyyy_1[i] * fe_0 + g_xxxyy_xyyyy_1[i] * pa_x[i];

        g_xxxxyy_xyyyz_0[i] = 3.0 * g_xxyy_xyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_xyyyz_1[i] * fz_be_0 + g_xxxyy_yyyz_1[i] * fe_0 + g_xxxyy_xyyyz_1[i] * pa_x[i];

        g_xxxxyy_xyyzz_0[i] = 3.0 * g_xxyy_xyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_xyyzz_1[i] * fz_be_0 + g_xxxyy_yyzz_1[i] * fe_0 + g_xxxyy_xyyzz_1[i] * pa_x[i];

        g_xxxxyy_xyzzz_0[i] = 3.0 * g_xxyy_xyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_xyzzz_1[i] * fz_be_0 + g_xxxyy_yzzz_1[i] * fe_0 + g_xxxyy_xyzzz_1[i] * pa_x[i];

        g_xxxxyy_xzzzz_0[i] = g_xxxx_xzzzz_0[i] * fbe_0 - g_xxxx_xzzzz_1[i] * fz_be_0 + g_xxxxy_xzzzz_1[i] * pa_y[i];

        g_xxxxyy_yyyyy_0[i] = 3.0 * g_xxyy_yyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_yyyyy_1[i] * fz_be_0 + g_xxxyy_yyyyy_1[i] * pa_x[i];

        g_xxxxyy_yyyyz_0[i] = 3.0 * g_xxyy_yyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_yyyyz_1[i] * fz_be_0 + g_xxxyy_yyyyz_1[i] * pa_x[i];

        g_xxxxyy_yyyzz_0[i] = 3.0 * g_xxyy_yyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_yyyzz_1[i] * fz_be_0 + g_xxxyy_yyyzz_1[i] * pa_x[i];

        g_xxxxyy_yyzzz_0[i] = 3.0 * g_xxyy_yyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_yyzzz_1[i] * fz_be_0 + g_xxxyy_yyzzz_1[i] * pa_x[i];

        g_xxxxyy_yzzzz_0[i] = 3.0 * g_xxyy_yzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_yzzzz_1[i] * fz_be_0 + g_xxxyy_yzzzz_1[i] * pa_x[i];

        g_xxxxyy_zzzzz_0[i] = 3.0 * g_xxyy_zzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_zzzzz_1[i] * fz_be_0 + g_xxxyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 84-105 components of targeted buffer : IH

    auto g_xxxxyz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 84);

    auto g_xxxxyz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 85);

    auto g_xxxxyz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 86);

    auto g_xxxxyz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 87);

    auto g_xxxxyz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 88);

    auto g_xxxxyz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 89);

    auto g_xxxxyz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 90);

    auto g_xxxxyz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 91);

    auto g_xxxxyz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 92);

    auto g_xxxxyz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 93);

    auto g_xxxxyz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 94);

    auto g_xxxxyz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 95);

    auto g_xxxxyz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 96);

    auto g_xxxxyz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 97);

    auto g_xxxxyz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 98);

    auto g_xxxxyz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 99);

    auto g_xxxxyz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 100);

    auto g_xxxxyz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 101);

    auto g_xxxxyz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 102);

    auto g_xxxxyz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 103);

    auto g_xxxxyz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 104);

    #pragma omp simd aligned(g_xxxxy_xxxxy_1, g_xxxxy_xxxyy_1, g_xxxxy_xxyyy_1, g_xxxxy_xyyyy_1, g_xxxxy_yyyyy_1, g_xxxxyz_xxxxx_0, g_xxxxyz_xxxxy_0, g_xxxxyz_xxxxz_0, g_xxxxyz_xxxyy_0, g_xxxxyz_xxxyz_0, g_xxxxyz_xxxzz_0, g_xxxxyz_xxyyy_0, g_xxxxyz_xxyyz_0, g_xxxxyz_xxyzz_0, g_xxxxyz_xxzzz_0, g_xxxxyz_xyyyy_0, g_xxxxyz_xyyyz_0, g_xxxxyz_xyyzz_0, g_xxxxyz_xyzzz_0, g_xxxxyz_xzzzz_0, g_xxxxyz_yyyyy_0, g_xxxxyz_yyyyz_0, g_xxxxyz_yyyzz_0, g_xxxxyz_yyzzz_0, g_xxxxyz_yzzzz_0, g_xxxxyz_zzzzz_0, g_xxxxz_xxxxx_1, g_xxxxz_xxxxz_1, g_xxxxz_xxxyz_1, g_xxxxz_xxxz_1, g_xxxxz_xxxzz_1, g_xxxxz_xxyyz_1, g_xxxxz_xxyz_1, g_xxxxz_xxyzz_1, g_xxxxz_xxzz_1, g_xxxxz_xxzzz_1, g_xxxxz_xyyyz_1, g_xxxxz_xyyz_1, g_xxxxz_xyyzz_1, g_xxxxz_xyzz_1, g_xxxxz_xyzzz_1, g_xxxxz_xzzz_1, g_xxxxz_xzzzz_1, g_xxxxz_yyyyz_1, g_xxxxz_yyyz_1, g_xxxxz_yyyzz_1, g_xxxxz_yyzz_1, g_xxxxz_yyzzz_1, g_xxxxz_yzzz_1, g_xxxxz_yzzzz_1, g_xxxxz_zzzz_1, g_xxxxz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyz_xxxxx_0[i] = g_xxxxz_xxxxx_1[i] * pa_y[i];

        g_xxxxyz_xxxxy_0[i] = g_xxxxy_xxxxy_1[i] * pa_z[i];

        g_xxxxyz_xxxxz_0[i] = g_xxxxz_xxxxz_1[i] * pa_y[i];

        g_xxxxyz_xxxyy_0[i] = g_xxxxy_xxxyy_1[i] * pa_z[i];

        g_xxxxyz_xxxyz_0[i] = g_xxxxz_xxxz_1[i] * fe_0 + g_xxxxz_xxxyz_1[i] * pa_y[i];

        g_xxxxyz_xxxzz_0[i] = g_xxxxz_xxxzz_1[i] * pa_y[i];

        g_xxxxyz_xxyyy_0[i] = g_xxxxy_xxyyy_1[i] * pa_z[i];

        g_xxxxyz_xxyyz_0[i] = 2.0 * g_xxxxz_xxyz_1[i] * fe_0 + g_xxxxz_xxyyz_1[i] * pa_y[i];

        g_xxxxyz_xxyzz_0[i] = g_xxxxz_xxzz_1[i] * fe_0 + g_xxxxz_xxyzz_1[i] * pa_y[i];

        g_xxxxyz_xxzzz_0[i] = g_xxxxz_xxzzz_1[i] * pa_y[i];

        g_xxxxyz_xyyyy_0[i] = g_xxxxy_xyyyy_1[i] * pa_z[i];

        g_xxxxyz_xyyyz_0[i] = 3.0 * g_xxxxz_xyyz_1[i] * fe_0 + g_xxxxz_xyyyz_1[i] * pa_y[i];

        g_xxxxyz_xyyzz_0[i] = 2.0 * g_xxxxz_xyzz_1[i] * fe_0 + g_xxxxz_xyyzz_1[i] * pa_y[i];

        g_xxxxyz_xyzzz_0[i] = g_xxxxz_xzzz_1[i] * fe_0 + g_xxxxz_xyzzz_1[i] * pa_y[i];

        g_xxxxyz_xzzzz_0[i] = g_xxxxz_xzzzz_1[i] * pa_y[i];

        g_xxxxyz_yyyyy_0[i] = g_xxxxy_yyyyy_1[i] * pa_z[i];

        g_xxxxyz_yyyyz_0[i] = 4.0 * g_xxxxz_yyyz_1[i] * fe_0 + g_xxxxz_yyyyz_1[i] * pa_y[i];

        g_xxxxyz_yyyzz_0[i] = 3.0 * g_xxxxz_yyzz_1[i] * fe_0 + g_xxxxz_yyyzz_1[i] * pa_y[i];

        g_xxxxyz_yyzzz_0[i] = 2.0 * g_xxxxz_yzzz_1[i] * fe_0 + g_xxxxz_yyzzz_1[i] * pa_y[i];

        g_xxxxyz_yzzzz_0[i] = g_xxxxz_zzzz_1[i] * fe_0 + g_xxxxz_yzzzz_1[i] * pa_y[i];

        g_xxxxyz_zzzzz_0[i] = g_xxxxz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : IH

    auto g_xxxxzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 105);

    auto g_xxxxzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 106);

    auto g_xxxxzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 107);

    auto g_xxxxzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 108);

    auto g_xxxxzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 109);

    auto g_xxxxzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 110);

    auto g_xxxxzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 111);

    auto g_xxxxzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 112);

    auto g_xxxxzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 113);

    auto g_xxxxzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 114);

    auto g_xxxxzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 115);

    auto g_xxxxzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 116);

    auto g_xxxxzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 117);

    auto g_xxxxzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 118);

    auto g_xxxxzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 119);

    auto g_xxxxzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 120);

    auto g_xxxxzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 121);

    auto g_xxxxzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 122);

    auto g_xxxxzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 123);

    auto g_xxxxzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 124);

    auto g_xxxxzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 125);

    #pragma omp simd aligned(g_xxxx_xxxxx_0, g_xxxx_xxxxx_1, g_xxxx_xxxxy_0, g_xxxx_xxxxy_1, g_xxxx_xxxyy_0, g_xxxx_xxxyy_1, g_xxxx_xxyyy_0, g_xxxx_xxyyy_1, g_xxxx_xyyyy_0, g_xxxx_xyyyy_1, g_xxxxz_xxxxx_1, g_xxxxz_xxxxy_1, g_xxxxz_xxxyy_1, g_xxxxz_xxyyy_1, g_xxxxz_xyyyy_1, g_xxxxzz_xxxxx_0, g_xxxxzz_xxxxy_0, g_xxxxzz_xxxxz_0, g_xxxxzz_xxxyy_0, g_xxxxzz_xxxyz_0, g_xxxxzz_xxxzz_0, g_xxxxzz_xxyyy_0, g_xxxxzz_xxyyz_0, g_xxxxzz_xxyzz_0, g_xxxxzz_xxzzz_0, g_xxxxzz_xyyyy_0, g_xxxxzz_xyyyz_0, g_xxxxzz_xyyzz_0, g_xxxxzz_xyzzz_0, g_xxxxzz_xzzzz_0, g_xxxxzz_yyyyy_0, g_xxxxzz_yyyyz_0, g_xxxxzz_yyyzz_0, g_xxxxzz_yyzzz_0, g_xxxxzz_yzzzz_0, g_xxxxzz_zzzzz_0, g_xxxzz_xxxxz_1, g_xxxzz_xxxyz_1, g_xxxzz_xxxz_1, g_xxxzz_xxxzz_1, g_xxxzz_xxyyz_1, g_xxxzz_xxyz_1, g_xxxzz_xxyzz_1, g_xxxzz_xxzz_1, g_xxxzz_xxzzz_1, g_xxxzz_xyyyz_1, g_xxxzz_xyyz_1, g_xxxzz_xyyzz_1, g_xxxzz_xyzz_1, g_xxxzz_xyzzz_1, g_xxxzz_xzzz_1, g_xxxzz_xzzzz_1, g_xxxzz_yyyyy_1, g_xxxzz_yyyyz_1, g_xxxzz_yyyz_1, g_xxxzz_yyyzz_1, g_xxxzz_yyzz_1, g_xxxzz_yyzzz_1, g_xxxzz_yzzz_1, g_xxxzz_yzzzz_1, g_xxxzz_zzzz_1, g_xxxzz_zzzzz_1, g_xxzz_xxxxz_0, g_xxzz_xxxxz_1, g_xxzz_xxxyz_0, g_xxzz_xxxyz_1, g_xxzz_xxxzz_0, g_xxzz_xxxzz_1, g_xxzz_xxyyz_0, g_xxzz_xxyyz_1, g_xxzz_xxyzz_0, g_xxzz_xxyzz_1, g_xxzz_xxzzz_0, g_xxzz_xxzzz_1, g_xxzz_xyyyz_0, g_xxzz_xyyyz_1, g_xxzz_xyyzz_0, g_xxzz_xyyzz_1, g_xxzz_xyzzz_0, g_xxzz_xyzzz_1, g_xxzz_xzzzz_0, g_xxzz_xzzzz_1, g_xxzz_yyyyy_0, g_xxzz_yyyyy_1, g_xxzz_yyyyz_0, g_xxzz_yyyyz_1, g_xxzz_yyyzz_0, g_xxzz_yyyzz_1, g_xxzz_yyzzz_0, g_xxzz_yyzzz_1, g_xxzz_yzzzz_0, g_xxzz_yzzzz_1, g_xxzz_zzzzz_0, g_xxzz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxzz_xxxxx_0[i] = g_xxxx_xxxxx_0[i] * fbe_0 - g_xxxx_xxxxx_1[i] * fz_be_0 + g_xxxxz_xxxxx_1[i] * pa_z[i];

        g_xxxxzz_xxxxy_0[i] = g_xxxx_xxxxy_0[i] * fbe_0 - g_xxxx_xxxxy_1[i] * fz_be_0 + g_xxxxz_xxxxy_1[i] * pa_z[i];

        g_xxxxzz_xxxxz_0[i] = 3.0 * g_xxzz_xxxxz_0[i] * fbe_0 - 3.0 * g_xxzz_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxxzz_xxxz_1[i] * fe_0 + g_xxxzz_xxxxz_1[i] * pa_x[i];

        g_xxxxzz_xxxyy_0[i] = g_xxxx_xxxyy_0[i] * fbe_0 - g_xxxx_xxxyy_1[i] * fz_be_0 + g_xxxxz_xxxyy_1[i] * pa_z[i];

        g_xxxxzz_xxxyz_0[i] = 3.0 * g_xxzz_xxxyz_0[i] * fbe_0 - 3.0 * g_xxzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxzz_xxyz_1[i] * fe_0 + g_xxxzz_xxxyz_1[i] * pa_x[i];

        g_xxxxzz_xxxzz_0[i] = 3.0 * g_xxzz_xxxzz_0[i] * fbe_0 - 3.0 * g_xxzz_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_xxzz_1[i] * fe_0 + g_xxxzz_xxxzz_1[i] * pa_x[i];

        g_xxxxzz_xxyyy_0[i] = g_xxxx_xxyyy_0[i] * fbe_0 - g_xxxx_xxyyy_1[i] * fz_be_0 + g_xxxxz_xxyyy_1[i] * pa_z[i];

        g_xxxxzz_xxyyz_0[i] = 3.0 * g_xxzz_xxyyz_0[i] * fbe_0 - 3.0 * g_xxzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxzz_xyyz_1[i] * fe_0 + g_xxxzz_xxyyz_1[i] * pa_x[i];

        g_xxxxzz_xxyzz_0[i] = 3.0 * g_xxzz_xxyzz_0[i] * fbe_0 - 3.0 * g_xxzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_xyzz_1[i] * fe_0 + g_xxxzz_xxyzz_1[i] * pa_x[i];

        g_xxxxzz_xxzzz_0[i] = 3.0 * g_xxzz_xxzzz_0[i] * fbe_0 - 3.0 * g_xxzz_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_xzzz_1[i] * fe_0 + g_xxxzz_xxzzz_1[i] * pa_x[i];

        g_xxxxzz_xyyyy_0[i] = g_xxxx_xyyyy_0[i] * fbe_0 - g_xxxx_xyyyy_1[i] * fz_be_0 + g_xxxxz_xyyyy_1[i] * pa_z[i];

        g_xxxxzz_xyyyz_0[i] = 3.0 * g_xxzz_xyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_xyyyz_1[i] * fz_be_0 + g_xxxzz_yyyz_1[i] * fe_0 + g_xxxzz_xyyyz_1[i] * pa_x[i];

        g_xxxxzz_xyyzz_0[i] = 3.0 * g_xxzz_xyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_xyyzz_1[i] * fz_be_0 + g_xxxzz_yyzz_1[i] * fe_0 + g_xxxzz_xyyzz_1[i] * pa_x[i];

        g_xxxxzz_xyzzz_0[i] = 3.0 * g_xxzz_xyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_xyzzz_1[i] * fz_be_0 + g_xxxzz_yzzz_1[i] * fe_0 + g_xxxzz_xyzzz_1[i] * pa_x[i];

        g_xxxxzz_xzzzz_0[i] = 3.0 * g_xxzz_xzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_xzzzz_1[i] * fz_be_0 + g_xxxzz_zzzz_1[i] * fe_0 + g_xxxzz_xzzzz_1[i] * pa_x[i];

        g_xxxxzz_yyyyy_0[i] = 3.0 * g_xxzz_yyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_yyyyy_1[i] * fz_be_0 + g_xxxzz_yyyyy_1[i] * pa_x[i];

        g_xxxxzz_yyyyz_0[i] = 3.0 * g_xxzz_yyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_yyyyz_1[i] * fz_be_0 + g_xxxzz_yyyyz_1[i] * pa_x[i];

        g_xxxxzz_yyyzz_0[i] = 3.0 * g_xxzz_yyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_yyyzz_1[i] * fz_be_0 + g_xxxzz_yyyzz_1[i] * pa_x[i];

        g_xxxxzz_yyzzz_0[i] = 3.0 * g_xxzz_yyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_yyzzz_1[i] * fz_be_0 + g_xxxzz_yyzzz_1[i] * pa_x[i];

        g_xxxxzz_yzzzz_0[i] = 3.0 * g_xxzz_yzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_yzzzz_1[i] * fz_be_0 + g_xxxzz_yzzzz_1[i] * pa_x[i];

        g_xxxxzz_zzzzz_0[i] = 3.0 * g_xxzz_zzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_zzzzz_1[i] * fz_be_0 + g_xxxzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 126-147 components of targeted buffer : IH

    auto g_xxxyyy_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 126);

    auto g_xxxyyy_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 127);

    auto g_xxxyyy_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 128);

    auto g_xxxyyy_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 129);

    auto g_xxxyyy_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 130);

    auto g_xxxyyy_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 131);

    auto g_xxxyyy_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 132);

    auto g_xxxyyy_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 133);

    auto g_xxxyyy_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 134);

    auto g_xxxyyy_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 135);

    auto g_xxxyyy_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 136);

    auto g_xxxyyy_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 137);

    auto g_xxxyyy_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 138);

    auto g_xxxyyy_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 139);

    auto g_xxxyyy_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 140);

    auto g_xxxyyy_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 141);

    auto g_xxxyyy_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 142);

    auto g_xxxyyy_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 143);

    auto g_xxxyyy_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 144);

    auto g_xxxyyy_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 145);

    auto g_xxxyyy_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 146);

    #pragma omp simd aligned(g_xxxy_xxxxx_0, g_xxxy_xxxxx_1, g_xxxy_xxxxz_0, g_xxxy_xxxxz_1, g_xxxy_xxxzz_0, g_xxxy_xxxzz_1, g_xxxy_xxzzz_0, g_xxxy_xxzzz_1, g_xxxy_xzzzz_0, g_xxxy_xzzzz_1, g_xxxyy_xxxxx_1, g_xxxyy_xxxxz_1, g_xxxyy_xxxzz_1, g_xxxyy_xxzzz_1, g_xxxyy_xzzzz_1, g_xxxyyy_xxxxx_0, g_xxxyyy_xxxxy_0, g_xxxyyy_xxxxz_0, g_xxxyyy_xxxyy_0, g_xxxyyy_xxxyz_0, g_xxxyyy_xxxzz_0, g_xxxyyy_xxyyy_0, g_xxxyyy_xxyyz_0, g_xxxyyy_xxyzz_0, g_xxxyyy_xxzzz_0, g_xxxyyy_xyyyy_0, g_xxxyyy_xyyyz_0, g_xxxyyy_xyyzz_0, g_xxxyyy_xyzzz_0, g_xxxyyy_xzzzz_0, g_xxxyyy_yyyyy_0, g_xxxyyy_yyyyz_0, g_xxxyyy_yyyzz_0, g_xxxyyy_yyzzz_0, g_xxxyyy_yzzzz_0, g_xxxyyy_zzzzz_0, g_xxyyy_xxxxy_1, g_xxyyy_xxxy_1, g_xxyyy_xxxyy_1, g_xxyyy_xxxyz_1, g_xxyyy_xxyy_1, g_xxyyy_xxyyy_1, g_xxyyy_xxyyz_1, g_xxyyy_xxyz_1, g_xxyyy_xxyzz_1, g_xxyyy_xyyy_1, g_xxyyy_xyyyy_1, g_xxyyy_xyyyz_1, g_xxyyy_xyyz_1, g_xxyyy_xyyzz_1, g_xxyyy_xyzz_1, g_xxyyy_xyzzz_1, g_xxyyy_yyyy_1, g_xxyyy_yyyyy_1, g_xxyyy_yyyyz_1, g_xxyyy_yyyz_1, g_xxyyy_yyyzz_1, g_xxyyy_yyzz_1, g_xxyyy_yyzzz_1, g_xxyyy_yzzz_1, g_xxyyy_yzzzz_1, g_xxyyy_zzzzz_1, g_xyyy_xxxxy_0, g_xyyy_xxxxy_1, g_xyyy_xxxyy_0, g_xyyy_xxxyy_1, g_xyyy_xxxyz_0, g_xyyy_xxxyz_1, g_xyyy_xxyyy_0, g_xyyy_xxyyy_1, g_xyyy_xxyyz_0, g_xyyy_xxyyz_1, g_xyyy_xxyzz_0, g_xyyy_xxyzz_1, g_xyyy_xyyyy_0, g_xyyy_xyyyy_1, g_xyyy_xyyyz_0, g_xyyy_xyyyz_1, g_xyyy_xyyzz_0, g_xyyy_xyyzz_1, g_xyyy_xyzzz_0, g_xyyy_xyzzz_1, g_xyyy_yyyyy_0, g_xyyy_yyyyy_1, g_xyyy_yyyyz_0, g_xyyy_yyyyz_1, g_xyyy_yyyzz_0, g_xyyy_yyyzz_1, g_xyyy_yyzzz_0, g_xyyy_yyzzz_1, g_xyyy_yzzzz_0, g_xyyy_yzzzz_1, g_xyyy_zzzzz_0, g_xyyy_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyy_xxxxx_0[i] = 2.0 * g_xxxy_xxxxx_0[i] * fbe_0 - 2.0 * g_xxxy_xxxxx_1[i] * fz_be_0 + g_xxxyy_xxxxx_1[i] * pa_y[i];

        g_xxxyyy_xxxxy_0[i] = 2.0 * g_xyyy_xxxxy_0[i] * fbe_0 - 2.0 * g_xyyy_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxyyy_xxxy_1[i] * fe_0 + g_xxyyy_xxxxy_1[i] * pa_x[i];

        g_xxxyyy_xxxxz_0[i] = 2.0 * g_xxxy_xxxxz_0[i] * fbe_0 - 2.0 * g_xxxy_xxxxz_1[i] * fz_be_0 + g_xxxyy_xxxxz_1[i] * pa_y[i];

        g_xxxyyy_xxxyy_0[i] = 2.0 * g_xyyy_xxxyy_0[i] * fbe_0 - 2.0 * g_xyyy_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxyyy_xxyy_1[i] * fe_0 + g_xxyyy_xxxyy_1[i] * pa_x[i];

        g_xxxyyy_xxxyz_0[i] = 2.0 * g_xyyy_xxxyz_0[i] * fbe_0 - 2.0 * g_xyyy_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxyyy_xxyz_1[i] * fe_0 + g_xxyyy_xxxyz_1[i] * pa_x[i];

        g_xxxyyy_xxxzz_0[i] = 2.0 * g_xxxy_xxxzz_0[i] * fbe_0 - 2.0 * g_xxxy_xxxzz_1[i] * fz_be_0 + g_xxxyy_xxxzz_1[i] * pa_y[i];

        g_xxxyyy_xxyyy_0[i] = 2.0 * g_xyyy_xxyyy_0[i] * fbe_0 - 2.0 * g_xyyy_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxyyy_xyyy_1[i] * fe_0 + g_xxyyy_xxyyy_1[i] * pa_x[i];

        g_xxxyyy_xxyyz_0[i] = 2.0 * g_xyyy_xxyyz_0[i] * fbe_0 - 2.0 * g_xyyy_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxyyy_xyyz_1[i] * fe_0 + g_xxyyy_xxyyz_1[i] * pa_x[i];

        g_xxxyyy_xxyzz_0[i] = 2.0 * g_xyyy_xxyzz_0[i] * fbe_0 - 2.0 * g_xyyy_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_xyzz_1[i] * fe_0 + g_xxyyy_xxyzz_1[i] * pa_x[i];

        g_xxxyyy_xxzzz_0[i] = 2.0 * g_xxxy_xxzzz_0[i] * fbe_0 - 2.0 * g_xxxy_xxzzz_1[i] * fz_be_0 + g_xxxyy_xxzzz_1[i] * pa_y[i];

        g_xxxyyy_xyyyy_0[i] = 2.0 * g_xyyy_xyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_xyyyy_1[i] * fz_be_0 + g_xxyyy_yyyy_1[i] * fe_0 + g_xxyyy_xyyyy_1[i] * pa_x[i];

        g_xxxyyy_xyyyz_0[i] = 2.0 * g_xyyy_xyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_xyyyz_1[i] * fz_be_0 + g_xxyyy_yyyz_1[i] * fe_0 + g_xxyyy_xyyyz_1[i] * pa_x[i];

        g_xxxyyy_xyyzz_0[i] = 2.0 * g_xyyy_xyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_xyyzz_1[i] * fz_be_0 + g_xxyyy_yyzz_1[i] * fe_0 + g_xxyyy_xyyzz_1[i] * pa_x[i];

        g_xxxyyy_xyzzz_0[i] = 2.0 * g_xyyy_xyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_xyzzz_1[i] * fz_be_0 + g_xxyyy_yzzz_1[i] * fe_0 + g_xxyyy_xyzzz_1[i] * pa_x[i];

        g_xxxyyy_xzzzz_0[i] = 2.0 * g_xxxy_xzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_xzzzz_1[i] * fz_be_0 + g_xxxyy_xzzzz_1[i] * pa_y[i];

        g_xxxyyy_yyyyy_0[i] = 2.0 * g_xyyy_yyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_yyyyy_1[i] * fz_be_0 + g_xxyyy_yyyyy_1[i] * pa_x[i];

        g_xxxyyy_yyyyz_0[i] = 2.0 * g_xyyy_yyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_yyyyz_1[i] * fz_be_0 + g_xxyyy_yyyyz_1[i] * pa_x[i];

        g_xxxyyy_yyyzz_0[i] = 2.0 * g_xyyy_yyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_yyyzz_1[i] * fz_be_0 + g_xxyyy_yyyzz_1[i] * pa_x[i];

        g_xxxyyy_yyzzz_0[i] = 2.0 * g_xyyy_yyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_yyzzz_1[i] * fz_be_0 + g_xxyyy_yyzzz_1[i] * pa_x[i];

        g_xxxyyy_yzzzz_0[i] = 2.0 * g_xyyy_yzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_yzzzz_1[i] * fz_be_0 + g_xxyyy_yzzzz_1[i] * pa_x[i];

        g_xxxyyy_zzzzz_0[i] = 2.0 * g_xyyy_zzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_zzzzz_1[i] * fz_be_0 + g_xxyyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 147-168 components of targeted buffer : IH

    auto g_xxxyyz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 147);

    auto g_xxxyyz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 148);

    auto g_xxxyyz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 149);

    auto g_xxxyyz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 150);

    auto g_xxxyyz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 151);

    auto g_xxxyyz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 152);

    auto g_xxxyyz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 153);

    auto g_xxxyyz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 154);

    auto g_xxxyyz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 155);

    auto g_xxxyyz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 156);

    auto g_xxxyyz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 157);

    auto g_xxxyyz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 158);

    auto g_xxxyyz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 159);

    auto g_xxxyyz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 160);

    auto g_xxxyyz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 161);

    auto g_xxxyyz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 162);

    auto g_xxxyyz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 163);

    auto g_xxxyyz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 164);

    auto g_xxxyyz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 165);

    auto g_xxxyyz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 166);

    auto g_xxxyyz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 167);

    #pragma omp simd aligned(g_xxxyy_xxxx_1, g_xxxyy_xxxxx_1, g_xxxyy_xxxxy_1, g_xxxyy_xxxxz_1, g_xxxyy_xxxy_1, g_xxxyy_xxxyy_1, g_xxxyy_xxxyz_1, g_xxxyy_xxxz_1, g_xxxyy_xxxzz_1, g_xxxyy_xxyy_1, g_xxxyy_xxyyy_1, g_xxxyy_xxyyz_1, g_xxxyy_xxyz_1, g_xxxyy_xxyzz_1, g_xxxyy_xxzz_1, g_xxxyy_xxzzz_1, g_xxxyy_xyyy_1, g_xxxyy_xyyyy_1, g_xxxyy_xyyyz_1, g_xxxyy_xyyz_1, g_xxxyy_xyyzz_1, g_xxxyy_xyzz_1, g_xxxyy_xyzzz_1, g_xxxyy_xzzz_1, g_xxxyy_xzzzz_1, g_xxxyy_yyyy_1, g_xxxyy_yyyyy_1, g_xxxyy_yyyyz_1, g_xxxyy_yyyz_1, g_xxxyy_yyyzz_1, g_xxxyy_yyzz_1, g_xxxyy_yyzzz_1, g_xxxyy_yzzz_1, g_xxxyy_yzzzz_1, g_xxxyy_zzzz_1, g_xxxyy_zzzzz_1, g_xxxyyz_xxxxx_0, g_xxxyyz_xxxxy_0, g_xxxyyz_xxxxz_0, g_xxxyyz_xxxyy_0, g_xxxyyz_xxxyz_0, g_xxxyyz_xxxzz_0, g_xxxyyz_xxyyy_0, g_xxxyyz_xxyyz_0, g_xxxyyz_xxyzz_0, g_xxxyyz_xxzzz_0, g_xxxyyz_xyyyy_0, g_xxxyyz_xyyyz_0, g_xxxyyz_xyyzz_0, g_xxxyyz_xyzzz_0, g_xxxyyz_xzzzz_0, g_xxxyyz_yyyyy_0, g_xxxyyz_yyyyz_0, g_xxxyyz_yyyzz_0, g_xxxyyz_yyzzz_0, g_xxxyyz_yzzzz_0, g_xxxyyz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyyz_xxxxx_0[i] = g_xxxyy_xxxxx_1[i] * pa_z[i];

        g_xxxyyz_xxxxy_0[i] = g_xxxyy_xxxxy_1[i] * pa_z[i];

        g_xxxyyz_xxxxz_0[i] = g_xxxyy_xxxx_1[i] * fe_0 + g_xxxyy_xxxxz_1[i] * pa_z[i];

        g_xxxyyz_xxxyy_0[i] = g_xxxyy_xxxyy_1[i] * pa_z[i];

        g_xxxyyz_xxxyz_0[i] = g_xxxyy_xxxy_1[i] * fe_0 + g_xxxyy_xxxyz_1[i] * pa_z[i];

        g_xxxyyz_xxxzz_0[i] = 2.0 * g_xxxyy_xxxz_1[i] * fe_0 + g_xxxyy_xxxzz_1[i] * pa_z[i];

        g_xxxyyz_xxyyy_0[i] = g_xxxyy_xxyyy_1[i] * pa_z[i];

        g_xxxyyz_xxyyz_0[i] = g_xxxyy_xxyy_1[i] * fe_0 + g_xxxyy_xxyyz_1[i] * pa_z[i];

        g_xxxyyz_xxyzz_0[i] = 2.0 * g_xxxyy_xxyz_1[i] * fe_0 + g_xxxyy_xxyzz_1[i] * pa_z[i];

        g_xxxyyz_xxzzz_0[i] = 3.0 * g_xxxyy_xxzz_1[i] * fe_0 + g_xxxyy_xxzzz_1[i] * pa_z[i];

        g_xxxyyz_xyyyy_0[i] = g_xxxyy_xyyyy_1[i] * pa_z[i];

        g_xxxyyz_xyyyz_0[i] = g_xxxyy_xyyy_1[i] * fe_0 + g_xxxyy_xyyyz_1[i] * pa_z[i];

        g_xxxyyz_xyyzz_0[i] = 2.0 * g_xxxyy_xyyz_1[i] * fe_0 + g_xxxyy_xyyzz_1[i] * pa_z[i];

        g_xxxyyz_xyzzz_0[i] = 3.0 * g_xxxyy_xyzz_1[i] * fe_0 + g_xxxyy_xyzzz_1[i] * pa_z[i];

        g_xxxyyz_xzzzz_0[i] = 4.0 * g_xxxyy_xzzz_1[i] * fe_0 + g_xxxyy_xzzzz_1[i] * pa_z[i];

        g_xxxyyz_yyyyy_0[i] = g_xxxyy_yyyyy_1[i] * pa_z[i];

        g_xxxyyz_yyyyz_0[i] = g_xxxyy_yyyy_1[i] * fe_0 + g_xxxyy_yyyyz_1[i] * pa_z[i];

        g_xxxyyz_yyyzz_0[i] = 2.0 * g_xxxyy_yyyz_1[i] * fe_0 + g_xxxyy_yyyzz_1[i] * pa_z[i];

        g_xxxyyz_yyzzz_0[i] = 3.0 * g_xxxyy_yyzz_1[i] * fe_0 + g_xxxyy_yyzzz_1[i] * pa_z[i];

        g_xxxyyz_yzzzz_0[i] = 4.0 * g_xxxyy_yzzz_1[i] * fe_0 + g_xxxyy_yzzzz_1[i] * pa_z[i];

        g_xxxyyz_zzzzz_0[i] = 5.0 * g_xxxyy_zzzz_1[i] * fe_0 + g_xxxyy_zzzzz_1[i] * pa_z[i];
    }

    // Set up 168-189 components of targeted buffer : IH

    auto g_xxxyzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 168);

    auto g_xxxyzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 169);

    auto g_xxxyzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 170);

    auto g_xxxyzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 171);

    auto g_xxxyzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 172);

    auto g_xxxyzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 173);

    auto g_xxxyzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 174);

    auto g_xxxyzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 175);

    auto g_xxxyzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 176);

    auto g_xxxyzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 177);

    auto g_xxxyzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 178);

    auto g_xxxyzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 179);

    auto g_xxxyzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 180);

    auto g_xxxyzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 181);

    auto g_xxxyzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 182);

    auto g_xxxyzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 183);

    auto g_xxxyzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 184);

    auto g_xxxyzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 185);

    auto g_xxxyzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 186);

    auto g_xxxyzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 187);

    auto g_xxxyzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 188);

    #pragma omp simd aligned(g_xxxyzz_xxxxx_0, g_xxxyzz_xxxxy_0, g_xxxyzz_xxxxz_0, g_xxxyzz_xxxyy_0, g_xxxyzz_xxxyz_0, g_xxxyzz_xxxzz_0, g_xxxyzz_xxyyy_0, g_xxxyzz_xxyyz_0, g_xxxyzz_xxyzz_0, g_xxxyzz_xxzzz_0, g_xxxyzz_xyyyy_0, g_xxxyzz_xyyyz_0, g_xxxyzz_xyyzz_0, g_xxxyzz_xyzzz_0, g_xxxyzz_xzzzz_0, g_xxxyzz_yyyyy_0, g_xxxyzz_yyyyz_0, g_xxxyzz_yyyzz_0, g_xxxyzz_yyzzz_0, g_xxxyzz_yzzzz_0, g_xxxyzz_zzzzz_0, g_xxxzz_xxxx_1, g_xxxzz_xxxxx_1, g_xxxzz_xxxxy_1, g_xxxzz_xxxxz_1, g_xxxzz_xxxy_1, g_xxxzz_xxxyy_1, g_xxxzz_xxxyz_1, g_xxxzz_xxxz_1, g_xxxzz_xxxzz_1, g_xxxzz_xxyy_1, g_xxxzz_xxyyy_1, g_xxxzz_xxyyz_1, g_xxxzz_xxyz_1, g_xxxzz_xxyzz_1, g_xxxzz_xxzz_1, g_xxxzz_xxzzz_1, g_xxxzz_xyyy_1, g_xxxzz_xyyyy_1, g_xxxzz_xyyyz_1, g_xxxzz_xyyz_1, g_xxxzz_xyyzz_1, g_xxxzz_xyzz_1, g_xxxzz_xyzzz_1, g_xxxzz_xzzz_1, g_xxxzz_xzzzz_1, g_xxxzz_yyyy_1, g_xxxzz_yyyyy_1, g_xxxzz_yyyyz_1, g_xxxzz_yyyz_1, g_xxxzz_yyyzz_1, g_xxxzz_yyzz_1, g_xxxzz_yyzzz_1, g_xxxzz_yzzz_1, g_xxxzz_yzzzz_1, g_xxxzz_zzzz_1, g_xxxzz_zzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyzz_xxxxx_0[i] = g_xxxzz_xxxxx_1[i] * pa_y[i];

        g_xxxyzz_xxxxy_0[i] = g_xxxzz_xxxx_1[i] * fe_0 + g_xxxzz_xxxxy_1[i] * pa_y[i];

        g_xxxyzz_xxxxz_0[i] = g_xxxzz_xxxxz_1[i] * pa_y[i];

        g_xxxyzz_xxxyy_0[i] = 2.0 * g_xxxzz_xxxy_1[i] * fe_0 + g_xxxzz_xxxyy_1[i] * pa_y[i];

        g_xxxyzz_xxxyz_0[i] = g_xxxzz_xxxz_1[i] * fe_0 + g_xxxzz_xxxyz_1[i] * pa_y[i];

        g_xxxyzz_xxxzz_0[i] = g_xxxzz_xxxzz_1[i] * pa_y[i];

        g_xxxyzz_xxyyy_0[i] = 3.0 * g_xxxzz_xxyy_1[i] * fe_0 + g_xxxzz_xxyyy_1[i] * pa_y[i];

        g_xxxyzz_xxyyz_0[i] = 2.0 * g_xxxzz_xxyz_1[i] * fe_0 + g_xxxzz_xxyyz_1[i] * pa_y[i];

        g_xxxyzz_xxyzz_0[i] = g_xxxzz_xxzz_1[i] * fe_0 + g_xxxzz_xxyzz_1[i] * pa_y[i];

        g_xxxyzz_xxzzz_0[i] = g_xxxzz_xxzzz_1[i] * pa_y[i];

        g_xxxyzz_xyyyy_0[i] = 4.0 * g_xxxzz_xyyy_1[i] * fe_0 + g_xxxzz_xyyyy_1[i] * pa_y[i];

        g_xxxyzz_xyyyz_0[i] = 3.0 * g_xxxzz_xyyz_1[i] * fe_0 + g_xxxzz_xyyyz_1[i] * pa_y[i];

        g_xxxyzz_xyyzz_0[i] = 2.0 * g_xxxzz_xyzz_1[i] * fe_0 + g_xxxzz_xyyzz_1[i] * pa_y[i];

        g_xxxyzz_xyzzz_0[i] = g_xxxzz_xzzz_1[i] * fe_0 + g_xxxzz_xyzzz_1[i] * pa_y[i];

        g_xxxyzz_xzzzz_0[i] = g_xxxzz_xzzzz_1[i] * pa_y[i];

        g_xxxyzz_yyyyy_0[i] = 5.0 * g_xxxzz_yyyy_1[i] * fe_0 + g_xxxzz_yyyyy_1[i] * pa_y[i];

        g_xxxyzz_yyyyz_0[i] = 4.0 * g_xxxzz_yyyz_1[i] * fe_0 + g_xxxzz_yyyyz_1[i] * pa_y[i];

        g_xxxyzz_yyyzz_0[i] = 3.0 * g_xxxzz_yyzz_1[i] * fe_0 + g_xxxzz_yyyzz_1[i] * pa_y[i];

        g_xxxyzz_yyzzz_0[i] = 2.0 * g_xxxzz_yzzz_1[i] * fe_0 + g_xxxzz_yyzzz_1[i] * pa_y[i];

        g_xxxyzz_yzzzz_0[i] = g_xxxzz_zzzz_1[i] * fe_0 + g_xxxzz_yzzzz_1[i] * pa_y[i];

        g_xxxyzz_zzzzz_0[i] = g_xxxzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 189-210 components of targeted buffer : IH

    auto g_xxxzzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 189);

    auto g_xxxzzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 190);

    auto g_xxxzzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 191);

    auto g_xxxzzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 192);

    auto g_xxxzzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 193);

    auto g_xxxzzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 194);

    auto g_xxxzzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 195);

    auto g_xxxzzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 196);

    auto g_xxxzzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 197);

    auto g_xxxzzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 198);

    auto g_xxxzzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 199);

    auto g_xxxzzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 200);

    auto g_xxxzzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 201);

    auto g_xxxzzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 202);

    auto g_xxxzzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 203);

    auto g_xxxzzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 204);

    auto g_xxxzzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 205);

    auto g_xxxzzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 206);

    auto g_xxxzzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 207);

    auto g_xxxzzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 208);

    auto g_xxxzzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 209);

    #pragma omp simd aligned(g_xxxz_xxxxx_0, g_xxxz_xxxxx_1, g_xxxz_xxxxy_0, g_xxxz_xxxxy_1, g_xxxz_xxxyy_0, g_xxxz_xxxyy_1, g_xxxz_xxyyy_0, g_xxxz_xxyyy_1, g_xxxz_xyyyy_0, g_xxxz_xyyyy_1, g_xxxzz_xxxxx_1, g_xxxzz_xxxxy_1, g_xxxzz_xxxyy_1, g_xxxzz_xxyyy_1, g_xxxzz_xyyyy_1, g_xxxzzz_xxxxx_0, g_xxxzzz_xxxxy_0, g_xxxzzz_xxxxz_0, g_xxxzzz_xxxyy_0, g_xxxzzz_xxxyz_0, g_xxxzzz_xxxzz_0, g_xxxzzz_xxyyy_0, g_xxxzzz_xxyyz_0, g_xxxzzz_xxyzz_0, g_xxxzzz_xxzzz_0, g_xxxzzz_xyyyy_0, g_xxxzzz_xyyyz_0, g_xxxzzz_xyyzz_0, g_xxxzzz_xyzzz_0, g_xxxzzz_xzzzz_0, g_xxxzzz_yyyyy_0, g_xxxzzz_yyyyz_0, g_xxxzzz_yyyzz_0, g_xxxzzz_yyzzz_0, g_xxxzzz_yzzzz_0, g_xxxzzz_zzzzz_0, g_xxzzz_xxxxz_1, g_xxzzz_xxxyz_1, g_xxzzz_xxxz_1, g_xxzzz_xxxzz_1, g_xxzzz_xxyyz_1, g_xxzzz_xxyz_1, g_xxzzz_xxyzz_1, g_xxzzz_xxzz_1, g_xxzzz_xxzzz_1, g_xxzzz_xyyyz_1, g_xxzzz_xyyz_1, g_xxzzz_xyyzz_1, g_xxzzz_xyzz_1, g_xxzzz_xyzzz_1, g_xxzzz_xzzz_1, g_xxzzz_xzzzz_1, g_xxzzz_yyyyy_1, g_xxzzz_yyyyz_1, g_xxzzz_yyyz_1, g_xxzzz_yyyzz_1, g_xxzzz_yyzz_1, g_xxzzz_yyzzz_1, g_xxzzz_yzzz_1, g_xxzzz_yzzzz_1, g_xxzzz_zzzz_1, g_xxzzz_zzzzz_1, g_xzzz_xxxxz_0, g_xzzz_xxxxz_1, g_xzzz_xxxyz_0, g_xzzz_xxxyz_1, g_xzzz_xxxzz_0, g_xzzz_xxxzz_1, g_xzzz_xxyyz_0, g_xzzz_xxyyz_1, g_xzzz_xxyzz_0, g_xzzz_xxyzz_1, g_xzzz_xxzzz_0, g_xzzz_xxzzz_1, g_xzzz_xyyyz_0, g_xzzz_xyyyz_1, g_xzzz_xyyzz_0, g_xzzz_xyyzz_1, g_xzzz_xyzzz_0, g_xzzz_xyzzz_1, g_xzzz_xzzzz_0, g_xzzz_xzzzz_1, g_xzzz_yyyyy_0, g_xzzz_yyyyy_1, g_xzzz_yyyyz_0, g_xzzz_yyyyz_1, g_xzzz_yyyzz_0, g_xzzz_yyyzz_1, g_xzzz_yyzzz_0, g_xzzz_yyzzz_1, g_xzzz_yzzzz_0, g_xzzz_yzzzz_1, g_xzzz_zzzzz_0, g_xzzz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzzz_xxxxx_0[i] = 2.0 * g_xxxz_xxxxx_0[i] * fbe_0 - 2.0 * g_xxxz_xxxxx_1[i] * fz_be_0 + g_xxxzz_xxxxx_1[i] * pa_z[i];

        g_xxxzzz_xxxxy_0[i] = 2.0 * g_xxxz_xxxxy_0[i] * fbe_0 - 2.0 * g_xxxz_xxxxy_1[i] * fz_be_0 + g_xxxzz_xxxxy_1[i] * pa_z[i];

        g_xxxzzz_xxxxz_0[i] = 2.0 * g_xzzz_xxxxz_0[i] * fbe_0 - 2.0 * g_xzzz_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxzzz_xxxz_1[i] * fe_0 + g_xxzzz_xxxxz_1[i] * pa_x[i];

        g_xxxzzz_xxxyy_0[i] = 2.0 * g_xxxz_xxxyy_0[i] * fbe_0 - 2.0 * g_xxxz_xxxyy_1[i] * fz_be_0 + g_xxxzz_xxxyy_1[i] * pa_z[i];

        g_xxxzzz_xxxyz_0[i] = 2.0 * g_xzzz_xxxyz_0[i] * fbe_0 - 2.0 * g_xzzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxzzz_xxyz_1[i] * fe_0 + g_xxzzz_xxxyz_1[i] * pa_x[i];

        g_xxxzzz_xxxzz_0[i] = 2.0 * g_xzzz_xxxzz_0[i] * fbe_0 - 2.0 * g_xzzz_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_xxzz_1[i] * fe_0 + g_xxzzz_xxxzz_1[i] * pa_x[i];

        g_xxxzzz_xxyyy_0[i] = 2.0 * g_xxxz_xxyyy_0[i] * fbe_0 - 2.0 * g_xxxz_xxyyy_1[i] * fz_be_0 + g_xxxzz_xxyyy_1[i] * pa_z[i];

        g_xxxzzz_xxyyz_0[i] = 2.0 * g_xzzz_xxyyz_0[i] * fbe_0 - 2.0 * g_xzzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxzzz_xyyz_1[i] * fe_0 + g_xxzzz_xxyyz_1[i] * pa_x[i];

        g_xxxzzz_xxyzz_0[i] = 2.0 * g_xzzz_xxyzz_0[i] * fbe_0 - 2.0 * g_xzzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_xyzz_1[i] * fe_0 + g_xxzzz_xxyzz_1[i] * pa_x[i];

        g_xxxzzz_xxzzz_0[i] = 2.0 * g_xzzz_xxzzz_0[i] * fbe_0 - 2.0 * g_xzzz_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_xzzz_1[i] * fe_0 + g_xxzzz_xxzzz_1[i] * pa_x[i];

        g_xxxzzz_xyyyy_0[i] = 2.0 * g_xxxz_xyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_xyyyy_1[i] * fz_be_0 + g_xxxzz_xyyyy_1[i] * pa_z[i];

        g_xxxzzz_xyyyz_0[i] = 2.0 * g_xzzz_xyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_xyyyz_1[i] * fz_be_0 + g_xxzzz_yyyz_1[i] * fe_0 + g_xxzzz_xyyyz_1[i] * pa_x[i];

        g_xxxzzz_xyyzz_0[i] = 2.0 * g_xzzz_xyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_xyyzz_1[i] * fz_be_0 + g_xxzzz_yyzz_1[i] * fe_0 + g_xxzzz_xyyzz_1[i] * pa_x[i];

        g_xxxzzz_xyzzz_0[i] = 2.0 * g_xzzz_xyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_xyzzz_1[i] * fz_be_0 + g_xxzzz_yzzz_1[i] * fe_0 + g_xxzzz_xyzzz_1[i] * pa_x[i];

        g_xxxzzz_xzzzz_0[i] = 2.0 * g_xzzz_xzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_xzzzz_1[i] * fz_be_0 + g_xxzzz_zzzz_1[i] * fe_0 + g_xxzzz_xzzzz_1[i] * pa_x[i];

        g_xxxzzz_yyyyy_0[i] = 2.0 * g_xzzz_yyyyy_0[i] * fbe_0 - 2.0 * g_xzzz_yyyyy_1[i] * fz_be_0 + g_xxzzz_yyyyy_1[i] * pa_x[i];

        g_xxxzzz_yyyyz_0[i] = 2.0 * g_xzzz_yyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_yyyyz_1[i] * fz_be_0 + g_xxzzz_yyyyz_1[i] * pa_x[i];

        g_xxxzzz_yyyzz_0[i] = 2.0 * g_xzzz_yyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_yyyzz_1[i] * fz_be_0 + g_xxzzz_yyyzz_1[i] * pa_x[i];

        g_xxxzzz_yyzzz_0[i] = 2.0 * g_xzzz_yyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_yyzzz_1[i] * fz_be_0 + g_xxzzz_yyzzz_1[i] * pa_x[i];

        g_xxxzzz_yzzzz_0[i] = 2.0 * g_xzzz_yzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_yzzzz_1[i] * fz_be_0 + g_xxzzz_yzzzz_1[i] * pa_x[i];

        g_xxxzzz_zzzzz_0[i] = 2.0 * g_xzzz_zzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_zzzzz_1[i] * fz_be_0 + g_xxzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 210-231 components of targeted buffer : IH

    auto g_xxyyyy_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 210);

    auto g_xxyyyy_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 211);

    auto g_xxyyyy_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 212);

    auto g_xxyyyy_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 213);

    auto g_xxyyyy_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 214);

    auto g_xxyyyy_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 215);

    auto g_xxyyyy_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 216);

    auto g_xxyyyy_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 217);

    auto g_xxyyyy_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 218);

    auto g_xxyyyy_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 219);

    auto g_xxyyyy_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 220);

    auto g_xxyyyy_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 221);

    auto g_xxyyyy_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 222);

    auto g_xxyyyy_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 223);

    auto g_xxyyyy_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 224);

    auto g_xxyyyy_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 225);

    auto g_xxyyyy_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 226);

    auto g_xxyyyy_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 227);

    auto g_xxyyyy_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 228);

    auto g_xxyyyy_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 229);

    auto g_xxyyyy_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 230);

    #pragma omp simd aligned(g_xxyy_xxxxx_0, g_xxyy_xxxxx_1, g_xxyy_xxxxz_0, g_xxyy_xxxxz_1, g_xxyy_xxxzz_0, g_xxyy_xxxzz_1, g_xxyy_xxzzz_0, g_xxyy_xxzzz_1, g_xxyy_xzzzz_0, g_xxyy_xzzzz_1, g_xxyyy_xxxxx_1, g_xxyyy_xxxxz_1, g_xxyyy_xxxzz_1, g_xxyyy_xxzzz_1, g_xxyyy_xzzzz_1, g_xxyyyy_xxxxx_0, g_xxyyyy_xxxxy_0, g_xxyyyy_xxxxz_0, g_xxyyyy_xxxyy_0, g_xxyyyy_xxxyz_0, g_xxyyyy_xxxzz_0, g_xxyyyy_xxyyy_0, g_xxyyyy_xxyyz_0, g_xxyyyy_xxyzz_0, g_xxyyyy_xxzzz_0, g_xxyyyy_xyyyy_0, g_xxyyyy_xyyyz_0, g_xxyyyy_xyyzz_0, g_xxyyyy_xyzzz_0, g_xxyyyy_xzzzz_0, g_xxyyyy_yyyyy_0, g_xxyyyy_yyyyz_0, g_xxyyyy_yyyzz_0, g_xxyyyy_yyzzz_0, g_xxyyyy_yzzzz_0, g_xxyyyy_zzzzz_0, g_xyyyy_xxxxy_1, g_xyyyy_xxxy_1, g_xyyyy_xxxyy_1, g_xyyyy_xxxyz_1, g_xyyyy_xxyy_1, g_xyyyy_xxyyy_1, g_xyyyy_xxyyz_1, g_xyyyy_xxyz_1, g_xyyyy_xxyzz_1, g_xyyyy_xyyy_1, g_xyyyy_xyyyy_1, g_xyyyy_xyyyz_1, g_xyyyy_xyyz_1, g_xyyyy_xyyzz_1, g_xyyyy_xyzz_1, g_xyyyy_xyzzz_1, g_xyyyy_yyyy_1, g_xyyyy_yyyyy_1, g_xyyyy_yyyyz_1, g_xyyyy_yyyz_1, g_xyyyy_yyyzz_1, g_xyyyy_yyzz_1, g_xyyyy_yyzzz_1, g_xyyyy_yzzz_1, g_xyyyy_yzzzz_1, g_xyyyy_zzzzz_1, g_yyyy_xxxxy_0, g_yyyy_xxxxy_1, g_yyyy_xxxyy_0, g_yyyy_xxxyy_1, g_yyyy_xxxyz_0, g_yyyy_xxxyz_1, g_yyyy_xxyyy_0, g_yyyy_xxyyy_1, g_yyyy_xxyyz_0, g_yyyy_xxyyz_1, g_yyyy_xxyzz_0, g_yyyy_xxyzz_1, g_yyyy_xyyyy_0, g_yyyy_xyyyy_1, g_yyyy_xyyyz_0, g_yyyy_xyyyz_1, g_yyyy_xyyzz_0, g_yyyy_xyyzz_1, g_yyyy_xyzzz_0, g_yyyy_xyzzz_1, g_yyyy_yyyyy_0, g_yyyy_yyyyy_1, g_yyyy_yyyyz_0, g_yyyy_yyyyz_1, g_yyyy_yyyzz_0, g_yyyy_yyyzz_1, g_yyyy_yyzzz_0, g_yyyy_yyzzz_1, g_yyyy_yzzzz_0, g_yyyy_yzzzz_1, g_yyyy_zzzzz_0, g_yyyy_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyy_xxxxx_0[i] = 3.0 * g_xxyy_xxxxx_0[i] * fbe_0 - 3.0 * g_xxyy_xxxxx_1[i] * fz_be_0 + g_xxyyy_xxxxx_1[i] * pa_y[i];

        g_xxyyyy_xxxxy_0[i] = g_yyyy_xxxxy_0[i] * fbe_0 - g_yyyy_xxxxy_1[i] * fz_be_0 + 4.0 * g_xyyyy_xxxy_1[i] * fe_0 + g_xyyyy_xxxxy_1[i] * pa_x[i];

        g_xxyyyy_xxxxz_0[i] = 3.0 * g_xxyy_xxxxz_0[i] * fbe_0 - 3.0 * g_xxyy_xxxxz_1[i] * fz_be_0 + g_xxyyy_xxxxz_1[i] * pa_y[i];

        g_xxyyyy_xxxyy_0[i] = g_yyyy_xxxyy_0[i] * fbe_0 - g_yyyy_xxxyy_1[i] * fz_be_0 + 3.0 * g_xyyyy_xxyy_1[i] * fe_0 + g_xyyyy_xxxyy_1[i] * pa_x[i];

        g_xxyyyy_xxxyz_0[i] = g_yyyy_xxxyz_0[i] * fbe_0 - g_yyyy_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyyyy_xxyz_1[i] * fe_0 + g_xyyyy_xxxyz_1[i] * pa_x[i];

        g_xxyyyy_xxxzz_0[i] = 3.0 * g_xxyy_xxxzz_0[i] * fbe_0 - 3.0 * g_xxyy_xxxzz_1[i] * fz_be_0 + g_xxyyy_xxxzz_1[i] * pa_y[i];

        g_xxyyyy_xxyyy_0[i] = g_yyyy_xxyyy_0[i] * fbe_0 - g_yyyy_xxyyy_1[i] * fz_be_0 + 2.0 * g_xyyyy_xyyy_1[i] * fe_0 + g_xyyyy_xxyyy_1[i] * pa_x[i];

        g_xxyyyy_xxyyz_0[i] = g_yyyy_xxyyz_0[i] * fbe_0 - g_yyyy_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyyyy_xyyz_1[i] * fe_0 + g_xyyyy_xxyyz_1[i] * pa_x[i];

        g_xxyyyy_xxyzz_0[i] = g_yyyy_xxyzz_0[i] * fbe_0 - g_yyyy_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_xyzz_1[i] * fe_0 + g_xyyyy_xxyzz_1[i] * pa_x[i];

        g_xxyyyy_xxzzz_0[i] = 3.0 * g_xxyy_xxzzz_0[i] * fbe_0 - 3.0 * g_xxyy_xxzzz_1[i] * fz_be_0 + g_xxyyy_xxzzz_1[i] * pa_y[i];

        g_xxyyyy_xyyyy_0[i] = g_yyyy_xyyyy_0[i] * fbe_0 - g_yyyy_xyyyy_1[i] * fz_be_0 + g_xyyyy_yyyy_1[i] * fe_0 + g_xyyyy_xyyyy_1[i] * pa_x[i];

        g_xxyyyy_xyyyz_0[i] = g_yyyy_xyyyz_0[i] * fbe_0 - g_yyyy_xyyyz_1[i] * fz_be_0 + g_xyyyy_yyyz_1[i] * fe_0 + g_xyyyy_xyyyz_1[i] * pa_x[i];

        g_xxyyyy_xyyzz_0[i] = g_yyyy_xyyzz_0[i] * fbe_0 - g_yyyy_xyyzz_1[i] * fz_be_0 + g_xyyyy_yyzz_1[i] * fe_0 + g_xyyyy_xyyzz_1[i] * pa_x[i];

        g_xxyyyy_xyzzz_0[i] = g_yyyy_xyzzz_0[i] * fbe_0 - g_yyyy_xyzzz_1[i] * fz_be_0 + g_xyyyy_yzzz_1[i] * fe_0 + g_xyyyy_xyzzz_1[i] * pa_x[i];

        g_xxyyyy_xzzzz_0[i] = 3.0 * g_xxyy_xzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_xzzzz_1[i] * fz_be_0 + g_xxyyy_xzzzz_1[i] * pa_y[i];

        g_xxyyyy_yyyyy_0[i] = g_yyyy_yyyyy_0[i] * fbe_0 - g_yyyy_yyyyy_1[i] * fz_be_0 + g_xyyyy_yyyyy_1[i] * pa_x[i];

        g_xxyyyy_yyyyz_0[i] = g_yyyy_yyyyz_0[i] * fbe_0 - g_yyyy_yyyyz_1[i] * fz_be_0 + g_xyyyy_yyyyz_1[i] * pa_x[i];

        g_xxyyyy_yyyzz_0[i] = g_yyyy_yyyzz_0[i] * fbe_0 - g_yyyy_yyyzz_1[i] * fz_be_0 + g_xyyyy_yyyzz_1[i] * pa_x[i];

        g_xxyyyy_yyzzz_0[i] = g_yyyy_yyzzz_0[i] * fbe_0 - g_yyyy_yyzzz_1[i] * fz_be_0 + g_xyyyy_yyzzz_1[i] * pa_x[i];

        g_xxyyyy_yzzzz_0[i] = g_yyyy_yzzzz_0[i] * fbe_0 - g_yyyy_yzzzz_1[i] * fz_be_0 + g_xyyyy_yzzzz_1[i] * pa_x[i];

        g_xxyyyy_zzzzz_0[i] = g_yyyy_zzzzz_0[i] * fbe_0 - g_yyyy_zzzzz_1[i] * fz_be_0 + g_xyyyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 231-252 components of targeted buffer : IH

    auto g_xxyyyz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 231);

    auto g_xxyyyz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 232);

    auto g_xxyyyz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 233);

    auto g_xxyyyz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 234);

    auto g_xxyyyz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 235);

    auto g_xxyyyz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 236);

    auto g_xxyyyz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 237);

    auto g_xxyyyz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 238);

    auto g_xxyyyz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 239);

    auto g_xxyyyz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 240);

    auto g_xxyyyz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 241);

    auto g_xxyyyz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 242);

    auto g_xxyyyz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 243);

    auto g_xxyyyz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 244);

    auto g_xxyyyz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 245);

    auto g_xxyyyz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 246);

    auto g_xxyyyz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 247);

    auto g_xxyyyz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 248);

    auto g_xxyyyz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 249);

    auto g_xxyyyz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 250);

    auto g_xxyyyz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 251);

    #pragma omp simd aligned(g_xxyyy_xxxx_1, g_xxyyy_xxxxx_1, g_xxyyy_xxxxy_1, g_xxyyy_xxxxz_1, g_xxyyy_xxxy_1, g_xxyyy_xxxyy_1, g_xxyyy_xxxyz_1, g_xxyyy_xxxz_1, g_xxyyy_xxxzz_1, g_xxyyy_xxyy_1, g_xxyyy_xxyyy_1, g_xxyyy_xxyyz_1, g_xxyyy_xxyz_1, g_xxyyy_xxyzz_1, g_xxyyy_xxzz_1, g_xxyyy_xxzzz_1, g_xxyyy_xyyy_1, g_xxyyy_xyyyy_1, g_xxyyy_xyyyz_1, g_xxyyy_xyyz_1, g_xxyyy_xyyzz_1, g_xxyyy_xyzz_1, g_xxyyy_xyzzz_1, g_xxyyy_xzzz_1, g_xxyyy_xzzzz_1, g_xxyyy_yyyy_1, g_xxyyy_yyyyy_1, g_xxyyy_yyyyz_1, g_xxyyy_yyyz_1, g_xxyyy_yyyzz_1, g_xxyyy_yyzz_1, g_xxyyy_yyzzz_1, g_xxyyy_yzzz_1, g_xxyyy_yzzzz_1, g_xxyyy_zzzz_1, g_xxyyy_zzzzz_1, g_xxyyyz_xxxxx_0, g_xxyyyz_xxxxy_0, g_xxyyyz_xxxxz_0, g_xxyyyz_xxxyy_0, g_xxyyyz_xxxyz_0, g_xxyyyz_xxxzz_0, g_xxyyyz_xxyyy_0, g_xxyyyz_xxyyz_0, g_xxyyyz_xxyzz_0, g_xxyyyz_xxzzz_0, g_xxyyyz_xyyyy_0, g_xxyyyz_xyyyz_0, g_xxyyyz_xyyzz_0, g_xxyyyz_xyzzz_0, g_xxyyyz_xzzzz_0, g_xxyyyz_yyyyy_0, g_xxyyyz_yyyyz_0, g_xxyyyz_yyyzz_0, g_xxyyyz_yyzzz_0, g_xxyyyz_yzzzz_0, g_xxyyyz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyyz_xxxxx_0[i] = g_xxyyy_xxxxx_1[i] * pa_z[i];

        g_xxyyyz_xxxxy_0[i] = g_xxyyy_xxxxy_1[i] * pa_z[i];

        g_xxyyyz_xxxxz_0[i] = g_xxyyy_xxxx_1[i] * fe_0 + g_xxyyy_xxxxz_1[i] * pa_z[i];

        g_xxyyyz_xxxyy_0[i] = g_xxyyy_xxxyy_1[i] * pa_z[i];

        g_xxyyyz_xxxyz_0[i] = g_xxyyy_xxxy_1[i] * fe_0 + g_xxyyy_xxxyz_1[i] * pa_z[i];

        g_xxyyyz_xxxzz_0[i] = 2.0 * g_xxyyy_xxxz_1[i] * fe_0 + g_xxyyy_xxxzz_1[i] * pa_z[i];

        g_xxyyyz_xxyyy_0[i] = g_xxyyy_xxyyy_1[i] * pa_z[i];

        g_xxyyyz_xxyyz_0[i] = g_xxyyy_xxyy_1[i] * fe_0 + g_xxyyy_xxyyz_1[i] * pa_z[i];

        g_xxyyyz_xxyzz_0[i] = 2.0 * g_xxyyy_xxyz_1[i] * fe_0 + g_xxyyy_xxyzz_1[i] * pa_z[i];

        g_xxyyyz_xxzzz_0[i] = 3.0 * g_xxyyy_xxzz_1[i] * fe_0 + g_xxyyy_xxzzz_1[i] * pa_z[i];

        g_xxyyyz_xyyyy_0[i] = g_xxyyy_xyyyy_1[i] * pa_z[i];

        g_xxyyyz_xyyyz_0[i] = g_xxyyy_xyyy_1[i] * fe_0 + g_xxyyy_xyyyz_1[i] * pa_z[i];

        g_xxyyyz_xyyzz_0[i] = 2.0 * g_xxyyy_xyyz_1[i] * fe_0 + g_xxyyy_xyyzz_1[i] * pa_z[i];

        g_xxyyyz_xyzzz_0[i] = 3.0 * g_xxyyy_xyzz_1[i] * fe_0 + g_xxyyy_xyzzz_1[i] * pa_z[i];

        g_xxyyyz_xzzzz_0[i] = 4.0 * g_xxyyy_xzzz_1[i] * fe_0 + g_xxyyy_xzzzz_1[i] * pa_z[i];

        g_xxyyyz_yyyyy_0[i] = g_xxyyy_yyyyy_1[i] * pa_z[i];

        g_xxyyyz_yyyyz_0[i] = g_xxyyy_yyyy_1[i] * fe_0 + g_xxyyy_yyyyz_1[i] * pa_z[i];

        g_xxyyyz_yyyzz_0[i] = 2.0 * g_xxyyy_yyyz_1[i] * fe_0 + g_xxyyy_yyyzz_1[i] * pa_z[i];

        g_xxyyyz_yyzzz_0[i] = 3.0 * g_xxyyy_yyzz_1[i] * fe_0 + g_xxyyy_yyzzz_1[i] * pa_z[i];

        g_xxyyyz_yzzzz_0[i] = 4.0 * g_xxyyy_yzzz_1[i] * fe_0 + g_xxyyy_yzzzz_1[i] * pa_z[i];

        g_xxyyyz_zzzzz_0[i] = 5.0 * g_xxyyy_zzzz_1[i] * fe_0 + g_xxyyy_zzzzz_1[i] * pa_z[i];
    }

    // Set up 252-273 components of targeted buffer : IH

    auto g_xxyyzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 252);

    auto g_xxyyzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 253);

    auto g_xxyyzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 254);

    auto g_xxyyzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 255);

    auto g_xxyyzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 256);

    auto g_xxyyzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 257);

    auto g_xxyyzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 258);

    auto g_xxyyzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 259);

    auto g_xxyyzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 260);

    auto g_xxyyzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 261);

    auto g_xxyyzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 262);

    auto g_xxyyzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 263);

    auto g_xxyyzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 264);

    auto g_xxyyzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 265);

    auto g_xxyyzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 266);

    auto g_xxyyzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 267);

    auto g_xxyyzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 268);

    auto g_xxyyzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 269);

    auto g_xxyyzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 270);

    auto g_xxyyzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 271);

    auto g_xxyyzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 272);

    #pragma omp simd aligned(g_xxyy_xxxxy_0, g_xxyy_xxxxy_1, g_xxyy_xxxyy_0, g_xxyy_xxxyy_1, g_xxyy_xxyyy_0, g_xxyy_xxyyy_1, g_xxyy_xyyyy_0, g_xxyy_xyyyy_1, g_xxyyz_xxxxy_1, g_xxyyz_xxxyy_1, g_xxyyz_xxyyy_1, g_xxyyz_xyyyy_1, g_xxyyzz_xxxxx_0, g_xxyyzz_xxxxy_0, g_xxyyzz_xxxxz_0, g_xxyyzz_xxxyy_0, g_xxyyzz_xxxyz_0, g_xxyyzz_xxxzz_0, g_xxyyzz_xxyyy_0, g_xxyyzz_xxyyz_0, g_xxyyzz_xxyzz_0, g_xxyyzz_xxzzz_0, g_xxyyzz_xyyyy_0, g_xxyyzz_xyyyz_0, g_xxyyzz_xyyzz_0, g_xxyyzz_xyzzz_0, g_xxyyzz_xzzzz_0, g_xxyyzz_yyyyy_0, g_xxyyzz_yyyyz_0, g_xxyyzz_yyyzz_0, g_xxyyzz_yyzzz_0, g_xxyyzz_yzzzz_0, g_xxyyzz_zzzzz_0, g_xxyzz_xxxxx_1, g_xxyzz_xxxxz_1, g_xxyzz_xxxzz_1, g_xxyzz_xxzzz_1, g_xxyzz_xzzzz_1, g_xxzz_xxxxx_0, g_xxzz_xxxxx_1, g_xxzz_xxxxz_0, g_xxzz_xxxxz_1, g_xxzz_xxxzz_0, g_xxzz_xxxzz_1, g_xxzz_xxzzz_0, g_xxzz_xxzzz_1, g_xxzz_xzzzz_0, g_xxzz_xzzzz_1, g_xyyzz_xxxyz_1, g_xyyzz_xxyyz_1, g_xyyzz_xxyz_1, g_xyyzz_xxyzz_1, g_xyyzz_xyyyz_1, g_xyyzz_xyyz_1, g_xyyzz_xyyzz_1, g_xyyzz_xyzz_1, g_xyyzz_xyzzz_1, g_xyyzz_yyyyy_1, g_xyyzz_yyyyz_1, g_xyyzz_yyyz_1, g_xyyzz_yyyzz_1, g_xyyzz_yyzz_1, g_xyyzz_yyzzz_1, g_xyyzz_yzzz_1, g_xyyzz_yzzzz_1, g_xyyzz_zzzzz_1, g_yyzz_xxxyz_0, g_yyzz_xxxyz_1, g_yyzz_xxyyz_0, g_yyzz_xxyyz_1, g_yyzz_xxyzz_0, g_yyzz_xxyzz_1, g_yyzz_xyyyz_0, g_yyzz_xyyyz_1, g_yyzz_xyyzz_0, g_yyzz_xyyzz_1, g_yyzz_xyzzz_0, g_yyzz_xyzzz_1, g_yyzz_yyyyy_0, g_yyzz_yyyyy_1, g_yyzz_yyyyz_0, g_yyzz_yyyyz_1, g_yyzz_yyyzz_0, g_yyzz_yyyzz_1, g_yyzz_yyzzz_0, g_yyzz_yyzzz_1, g_yyzz_yzzzz_0, g_yyzz_yzzzz_1, g_yyzz_zzzzz_0, g_yyzz_zzzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyzz_xxxxx_0[i] = g_xxzz_xxxxx_0[i] * fbe_0 - g_xxzz_xxxxx_1[i] * fz_be_0 + g_xxyzz_xxxxx_1[i] * pa_y[i];

        g_xxyyzz_xxxxy_0[i] = g_xxyy_xxxxy_0[i] * fbe_0 - g_xxyy_xxxxy_1[i] * fz_be_0 + g_xxyyz_xxxxy_1[i] * pa_z[i];

        g_xxyyzz_xxxxz_0[i] = g_xxzz_xxxxz_0[i] * fbe_0 - g_xxzz_xxxxz_1[i] * fz_be_0 + g_xxyzz_xxxxz_1[i] * pa_y[i];

        g_xxyyzz_xxxyy_0[i] = g_xxyy_xxxyy_0[i] * fbe_0 - g_xxyy_xxxyy_1[i] * fz_be_0 + g_xxyyz_xxxyy_1[i] * pa_z[i];

        g_xxyyzz_xxxyz_0[i] = g_yyzz_xxxyz_0[i] * fbe_0 - g_yyzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyyzz_xxyz_1[i] * fe_0 + g_xyyzz_xxxyz_1[i] * pa_x[i];

        g_xxyyzz_xxxzz_0[i] = g_xxzz_xxxzz_0[i] * fbe_0 - g_xxzz_xxxzz_1[i] * fz_be_0 + g_xxyzz_xxxzz_1[i] * pa_y[i];

        g_xxyyzz_xxyyy_0[i] = g_xxyy_xxyyy_0[i] * fbe_0 - g_xxyy_xxyyy_1[i] * fz_be_0 + g_xxyyz_xxyyy_1[i] * pa_z[i];

        g_xxyyzz_xxyyz_0[i] = g_yyzz_xxyyz_0[i] * fbe_0 - g_yyzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyyzz_xyyz_1[i] * fe_0 + g_xyyzz_xxyyz_1[i] * pa_x[i];

        g_xxyyzz_xxyzz_0[i] = g_yyzz_xxyzz_0[i] * fbe_0 - g_yyzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_xyzz_1[i] * fe_0 + g_xyyzz_xxyzz_1[i] * pa_x[i];

        g_xxyyzz_xxzzz_0[i] = g_xxzz_xxzzz_0[i] * fbe_0 - g_xxzz_xxzzz_1[i] * fz_be_0 + g_xxyzz_xxzzz_1[i] * pa_y[i];

        g_xxyyzz_xyyyy_0[i] = g_xxyy_xyyyy_0[i] * fbe_0 - g_xxyy_xyyyy_1[i] * fz_be_0 + g_xxyyz_xyyyy_1[i] * pa_z[i];

        g_xxyyzz_xyyyz_0[i] = g_yyzz_xyyyz_0[i] * fbe_0 - g_yyzz_xyyyz_1[i] * fz_be_0 + g_xyyzz_yyyz_1[i] * fe_0 + g_xyyzz_xyyyz_1[i] * pa_x[i];

        g_xxyyzz_xyyzz_0[i] = g_yyzz_xyyzz_0[i] * fbe_0 - g_yyzz_xyyzz_1[i] * fz_be_0 + g_xyyzz_yyzz_1[i] * fe_0 + g_xyyzz_xyyzz_1[i] * pa_x[i];

        g_xxyyzz_xyzzz_0[i] = g_yyzz_xyzzz_0[i] * fbe_0 - g_yyzz_xyzzz_1[i] * fz_be_0 + g_xyyzz_yzzz_1[i] * fe_0 + g_xyyzz_xyzzz_1[i] * pa_x[i];

        g_xxyyzz_xzzzz_0[i] = g_xxzz_xzzzz_0[i] * fbe_0 - g_xxzz_xzzzz_1[i] * fz_be_0 + g_xxyzz_xzzzz_1[i] * pa_y[i];

        g_xxyyzz_yyyyy_0[i] = g_yyzz_yyyyy_0[i] * fbe_0 - g_yyzz_yyyyy_1[i] * fz_be_0 + g_xyyzz_yyyyy_1[i] * pa_x[i];

        g_xxyyzz_yyyyz_0[i] = g_yyzz_yyyyz_0[i] * fbe_0 - g_yyzz_yyyyz_1[i] * fz_be_0 + g_xyyzz_yyyyz_1[i] * pa_x[i];

        g_xxyyzz_yyyzz_0[i] = g_yyzz_yyyzz_0[i] * fbe_0 - g_yyzz_yyyzz_1[i] * fz_be_0 + g_xyyzz_yyyzz_1[i] * pa_x[i];

        g_xxyyzz_yyzzz_0[i] = g_yyzz_yyzzz_0[i] * fbe_0 - g_yyzz_yyzzz_1[i] * fz_be_0 + g_xyyzz_yyzzz_1[i] * pa_x[i];

        g_xxyyzz_yzzzz_0[i] = g_yyzz_yzzzz_0[i] * fbe_0 - g_yyzz_yzzzz_1[i] * fz_be_0 + g_xyyzz_yzzzz_1[i] * pa_x[i];

        g_xxyyzz_zzzzz_0[i] = g_yyzz_zzzzz_0[i] * fbe_0 - g_yyzz_zzzzz_1[i] * fz_be_0 + g_xyyzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 273-294 components of targeted buffer : IH

    auto g_xxyzzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 273);

    auto g_xxyzzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 274);

    auto g_xxyzzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 275);

    auto g_xxyzzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 276);

    auto g_xxyzzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 277);

    auto g_xxyzzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 278);

    auto g_xxyzzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 279);

    auto g_xxyzzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 280);

    auto g_xxyzzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 281);

    auto g_xxyzzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 282);

    auto g_xxyzzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 283);

    auto g_xxyzzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 284);

    auto g_xxyzzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 285);

    auto g_xxyzzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 286);

    auto g_xxyzzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 287);

    auto g_xxyzzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 288);

    auto g_xxyzzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 289);

    auto g_xxyzzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 290);

    auto g_xxyzzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 291);

    auto g_xxyzzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 292);

    auto g_xxyzzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 293);

    #pragma omp simd aligned(g_xxyzzz_xxxxx_0, g_xxyzzz_xxxxy_0, g_xxyzzz_xxxxz_0, g_xxyzzz_xxxyy_0, g_xxyzzz_xxxyz_0, g_xxyzzz_xxxzz_0, g_xxyzzz_xxyyy_0, g_xxyzzz_xxyyz_0, g_xxyzzz_xxyzz_0, g_xxyzzz_xxzzz_0, g_xxyzzz_xyyyy_0, g_xxyzzz_xyyyz_0, g_xxyzzz_xyyzz_0, g_xxyzzz_xyzzz_0, g_xxyzzz_xzzzz_0, g_xxyzzz_yyyyy_0, g_xxyzzz_yyyyz_0, g_xxyzzz_yyyzz_0, g_xxyzzz_yyzzz_0, g_xxyzzz_yzzzz_0, g_xxyzzz_zzzzz_0, g_xxzzz_xxxx_1, g_xxzzz_xxxxx_1, g_xxzzz_xxxxy_1, g_xxzzz_xxxxz_1, g_xxzzz_xxxy_1, g_xxzzz_xxxyy_1, g_xxzzz_xxxyz_1, g_xxzzz_xxxz_1, g_xxzzz_xxxzz_1, g_xxzzz_xxyy_1, g_xxzzz_xxyyy_1, g_xxzzz_xxyyz_1, g_xxzzz_xxyz_1, g_xxzzz_xxyzz_1, g_xxzzz_xxzz_1, g_xxzzz_xxzzz_1, g_xxzzz_xyyy_1, g_xxzzz_xyyyy_1, g_xxzzz_xyyyz_1, g_xxzzz_xyyz_1, g_xxzzz_xyyzz_1, g_xxzzz_xyzz_1, g_xxzzz_xyzzz_1, g_xxzzz_xzzz_1, g_xxzzz_xzzzz_1, g_xxzzz_yyyy_1, g_xxzzz_yyyyy_1, g_xxzzz_yyyyz_1, g_xxzzz_yyyz_1, g_xxzzz_yyyzz_1, g_xxzzz_yyzz_1, g_xxzzz_yyzzz_1, g_xxzzz_yzzz_1, g_xxzzz_yzzzz_1, g_xxzzz_zzzz_1, g_xxzzz_zzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzzz_xxxxx_0[i] = g_xxzzz_xxxxx_1[i] * pa_y[i];

        g_xxyzzz_xxxxy_0[i] = g_xxzzz_xxxx_1[i] * fe_0 + g_xxzzz_xxxxy_1[i] * pa_y[i];

        g_xxyzzz_xxxxz_0[i] = g_xxzzz_xxxxz_1[i] * pa_y[i];

        g_xxyzzz_xxxyy_0[i] = 2.0 * g_xxzzz_xxxy_1[i] * fe_0 + g_xxzzz_xxxyy_1[i] * pa_y[i];

        g_xxyzzz_xxxyz_0[i] = g_xxzzz_xxxz_1[i] * fe_0 + g_xxzzz_xxxyz_1[i] * pa_y[i];

        g_xxyzzz_xxxzz_0[i] = g_xxzzz_xxxzz_1[i] * pa_y[i];

        g_xxyzzz_xxyyy_0[i] = 3.0 * g_xxzzz_xxyy_1[i] * fe_0 + g_xxzzz_xxyyy_1[i] * pa_y[i];

        g_xxyzzz_xxyyz_0[i] = 2.0 * g_xxzzz_xxyz_1[i] * fe_0 + g_xxzzz_xxyyz_1[i] * pa_y[i];

        g_xxyzzz_xxyzz_0[i] = g_xxzzz_xxzz_1[i] * fe_0 + g_xxzzz_xxyzz_1[i] * pa_y[i];

        g_xxyzzz_xxzzz_0[i] = g_xxzzz_xxzzz_1[i] * pa_y[i];

        g_xxyzzz_xyyyy_0[i] = 4.0 * g_xxzzz_xyyy_1[i] * fe_0 + g_xxzzz_xyyyy_1[i] * pa_y[i];

        g_xxyzzz_xyyyz_0[i] = 3.0 * g_xxzzz_xyyz_1[i] * fe_0 + g_xxzzz_xyyyz_1[i] * pa_y[i];

        g_xxyzzz_xyyzz_0[i] = 2.0 * g_xxzzz_xyzz_1[i] * fe_0 + g_xxzzz_xyyzz_1[i] * pa_y[i];

        g_xxyzzz_xyzzz_0[i] = g_xxzzz_xzzz_1[i] * fe_0 + g_xxzzz_xyzzz_1[i] * pa_y[i];

        g_xxyzzz_xzzzz_0[i] = g_xxzzz_xzzzz_1[i] * pa_y[i];

        g_xxyzzz_yyyyy_0[i] = 5.0 * g_xxzzz_yyyy_1[i] * fe_0 + g_xxzzz_yyyyy_1[i] * pa_y[i];

        g_xxyzzz_yyyyz_0[i] = 4.0 * g_xxzzz_yyyz_1[i] * fe_0 + g_xxzzz_yyyyz_1[i] * pa_y[i];

        g_xxyzzz_yyyzz_0[i] = 3.0 * g_xxzzz_yyzz_1[i] * fe_0 + g_xxzzz_yyyzz_1[i] * pa_y[i];

        g_xxyzzz_yyzzz_0[i] = 2.0 * g_xxzzz_yzzz_1[i] * fe_0 + g_xxzzz_yyzzz_1[i] * pa_y[i];

        g_xxyzzz_yzzzz_0[i] = g_xxzzz_zzzz_1[i] * fe_0 + g_xxzzz_yzzzz_1[i] * pa_y[i];

        g_xxyzzz_zzzzz_0[i] = g_xxzzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 294-315 components of targeted buffer : IH

    auto g_xxzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 294);

    auto g_xxzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 295);

    auto g_xxzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 296);

    auto g_xxzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 297);

    auto g_xxzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 298);

    auto g_xxzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 299);

    auto g_xxzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 300);

    auto g_xxzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 301);

    auto g_xxzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 302);

    auto g_xxzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 303);

    auto g_xxzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 304);

    auto g_xxzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 305);

    auto g_xxzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 306);

    auto g_xxzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 307);

    auto g_xxzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 308);

    auto g_xxzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 309);

    auto g_xxzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 310);

    auto g_xxzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 311);

    auto g_xxzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 312);

    auto g_xxzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 313);

    auto g_xxzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 314);

    #pragma omp simd aligned(g_xxzz_xxxxx_0, g_xxzz_xxxxx_1, g_xxzz_xxxxy_0, g_xxzz_xxxxy_1, g_xxzz_xxxyy_0, g_xxzz_xxxyy_1, g_xxzz_xxyyy_0, g_xxzz_xxyyy_1, g_xxzz_xyyyy_0, g_xxzz_xyyyy_1, g_xxzzz_xxxxx_1, g_xxzzz_xxxxy_1, g_xxzzz_xxxyy_1, g_xxzzz_xxyyy_1, g_xxzzz_xyyyy_1, g_xxzzzz_xxxxx_0, g_xxzzzz_xxxxy_0, g_xxzzzz_xxxxz_0, g_xxzzzz_xxxyy_0, g_xxzzzz_xxxyz_0, g_xxzzzz_xxxzz_0, g_xxzzzz_xxyyy_0, g_xxzzzz_xxyyz_0, g_xxzzzz_xxyzz_0, g_xxzzzz_xxzzz_0, g_xxzzzz_xyyyy_0, g_xxzzzz_xyyyz_0, g_xxzzzz_xyyzz_0, g_xxzzzz_xyzzz_0, g_xxzzzz_xzzzz_0, g_xxzzzz_yyyyy_0, g_xxzzzz_yyyyz_0, g_xxzzzz_yyyzz_0, g_xxzzzz_yyzzz_0, g_xxzzzz_yzzzz_0, g_xxzzzz_zzzzz_0, g_xzzzz_xxxxz_1, g_xzzzz_xxxyz_1, g_xzzzz_xxxz_1, g_xzzzz_xxxzz_1, g_xzzzz_xxyyz_1, g_xzzzz_xxyz_1, g_xzzzz_xxyzz_1, g_xzzzz_xxzz_1, g_xzzzz_xxzzz_1, g_xzzzz_xyyyz_1, g_xzzzz_xyyz_1, g_xzzzz_xyyzz_1, g_xzzzz_xyzz_1, g_xzzzz_xyzzz_1, g_xzzzz_xzzz_1, g_xzzzz_xzzzz_1, g_xzzzz_yyyyy_1, g_xzzzz_yyyyz_1, g_xzzzz_yyyz_1, g_xzzzz_yyyzz_1, g_xzzzz_yyzz_1, g_xzzzz_yyzzz_1, g_xzzzz_yzzz_1, g_xzzzz_yzzzz_1, g_xzzzz_zzzz_1, g_xzzzz_zzzzz_1, g_zzzz_xxxxz_0, g_zzzz_xxxxz_1, g_zzzz_xxxyz_0, g_zzzz_xxxyz_1, g_zzzz_xxxzz_0, g_zzzz_xxxzz_1, g_zzzz_xxyyz_0, g_zzzz_xxyyz_1, g_zzzz_xxyzz_0, g_zzzz_xxyzz_1, g_zzzz_xxzzz_0, g_zzzz_xxzzz_1, g_zzzz_xyyyz_0, g_zzzz_xyyyz_1, g_zzzz_xyyzz_0, g_zzzz_xyyzz_1, g_zzzz_xyzzz_0, g_zzzz_xyzzz_1, g_zzzz_xzzzz_0, g_zzzz_xzzzz_1, g_zzzz_yyyyy_0, g_zzzz_yyyyy_1, g_zzzz_yyyyz_0, g_zzzz_yyyyz_1, g_zzzz_yyyzz_0, g_zzzz_yyyzz_1, g_zzzz_yyzzz_0, g_zzzz_yyzzz_1, g_zzzz_yzzzz_0, g_zzzz_yzzzz_1, g_zzzz_zzzzz_0, g_zzzz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzzz_xxxxx_0[i] = 3.0 * g_xxzz_xxxxx_0[i] * fbe_0 - 3.0 * g_xxzz_xxxxx_1[i] * fz_be_0 + g_xxzzz_xxxxx_1[i] * pa_z[i];

        g_xxzzzz_xxxxy_0[i] = 3.0 * g_xxzz_xxxxy_0[i] * fbe_0 - 3.0 * g_xxzz_xxxxy_1[i] * fz_be_0 + g_xxzzz_xxxxy_1[i] * pa_z[i];

        g_xxzzzz_xxxxz_0[i] = g_zzzz_xxxxz_0[i] * fbe_0 - g_zzzz_xxxxz_1[i] * fz_be_0 + 4.0 * g_xzzzz_xxxz_1[i] * fe_0 + g_xzzzz_xxxxz_1[i] * pa_x[i];

        g_xxzzzz_xxxyy_0[i] = 3.0 * g_xxzz_xxxyy_0[i] * fbe_0 - 3.0 * g_xxzz_xxxyy_1[i] * fz_be_0 + g_xxzzz_xxxyy_1[i] * pa_z[i];

        g_xxzzzz_xxxyz_0[i] = g_zzzz_xxxyz_0[i] * fbe_0 - g_zzzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xzzzz_xxyz_1[i] * fe_0 + g_xzzzz_xxxyz_1[i] * pa_x[i];

        g_xxzzzz_xxxzz_0[i] = g_zzzz_xxxzz_0[i] * fbe_0 - g_zzzz_xxxzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_xxzz_1[i] * fe_0 + g_xzzzz_xxxzz_1[i] * pa_x[i];

        g_xxzzzz_xxyyy_0[i] = 3.0 * g_xxzz_xxyyy_0[i] * fbe_0 - 3.0 * g_xxzz_xxyyy_1[i] * fz_be_0 + g_xxzzz_xxyyy_1[i] * pa_z[i];

        g_xxzzzz_xxyyz_0[i] = g_zzzz_xxyyz_0[i] * fbe_0 - g_zzzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xzzzz_xyyz_1[i] * fe_0 + g_xzzzz_xxyyz_1[i] * pa_x[i];

        g_xxzzzz_xxyzz_0[i] = g_zzzz_xxyzz_0[i] * fbe_0 - g_zzzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_xyzz_1[i] * fe_0 + g_xzzzz_xxyzz_1[i] * pa_x[i];

        g_xxzzzz_xxzzz_0[i] = g_zzzz_xxzzz_0[i] * fbe_0 - g_zzzz_xxzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_xzzz_1[i] * fe_0 + g_xzzzz_xxzzz_1[i] * pa_x[i];

        g_xxzzzz_xyyyy_0[i] = 3.0 * g_xxzz_xyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_xyyyy_1[i] * fz_be_0 + g_xxzzz_xyyyy_1[i] * pa_z[i];

        g_xxzzzz_xyyyz_0[i] = g_zzzz_xyyyz_0[i] * fbe_0 - g_zzzz_xyyyz_1[i] * fz_be_0 + g_xzzzz_yyyz_1[i] * fe_0 + g_xzzzz_xyyyz_1[i] * pa_x[i];

        g_xxzzzz_xyyzz_0[i] = g_zzzz_xyyzz_0[i] * fbe_0 - g_zzzz_xyyzz_1[i] * fz_be_0 + g_xzzzz_yyzz_1[i] * fe_0 + g_xzzzz_xyyzz_1[i] * pa_x[i];

        g_xxzzzz_xyzzz_0[i] = g_zzzz_xyzzz_0[i] * fbe_0 - g_zzzz_xyzzz_1[i] * fz_be_0 + g_xzzzz_yzzz_1[i] * fe_0 + g_xzzzz_xyzzz_1[i] * pa_x[i];

        g_xxzzzz_xzzzz_0[i] = g_zzzz_xzzzz_0[i] * fbe_0 - g_zzzz_xzzzz_1[i] * fz_be_0 + g_xzzzz_zzzz_1[i] * fe_0 + g_xzzzz_xzzzz_1[i] * pa_x[i];

        g_xxzzzz_yyyyy_0[i] = g_zzzz_yyyyy_0[i] * fbe_0 - g_zzzz_yyyyy_1[i] * fz_be_0 + g_xzzzz_yyyyy_1[i] * pa_x[i];

        g_xxzzzz_yyyyz_0[i] = g_zzzz_yyyyz_0[i] * fbe_0 - g_zzzz_yyyyz_1[i] * fz_be_0 + g_xzzzz_yyyyz_1[i] * pa_x[i];

        g_xxzzzz_yyyzz_0[i] = g_zzzz_yyyzz_0[i] * fbe_0 - g_zzzz_yyyzz_1[i] * fz_be_0 + g_xzzzz_yyyzz_1[i] * pa_x[i];

        g_xxzzzz_yyzzz_0[i] = g_zzzz_yyzzz_0[i] * fbe_0 - g_zzzz_yyzzz_1[i] * fz_be_0 + g_xzzzz_yyzzz_1[i] * pa_x[i];

        g_xxzzzz_yzzzz_0[i] = g_zzzz_yzzzz_0[i] * fbe_0 - g_zzzz_yzzzz_1[i] * fz_be_0 + g_xzzzz_yzzzz_1[i] * pa_x[i];

        g_xxzzzz_zzzzz_0[i] = g_zzzz_zzzzz_0[i] * fbe_0 - g_zzzz_zzzzz_1[i] * fz_be_0 + g_xzzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 315-336 components of targeted buffer : IH

    auto g_xyyyyy_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 315);

    auto g_xyyyyy_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 316);

    auto g_xyyyyy_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 317);

    auto g_xyyyyy_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 318);

    auto g_xyyyyy_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 319);

    auto g_xyyyyy_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 320);

    auto g_xyyyyy_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 321);

    auto g_xyyyyy_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 322);

    auto g_xyyyyy_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 323);

    auto g_xyyyyy_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 324);

    auto g_xyyyyy_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 325);

    auto g_xyyyyy_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 326);

    auto g_xyyyyy_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 327);

    auto g_xyyyyy_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 328);

    auto g_xyyyyy_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 329);

    auto g_xyyyyy_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 330);

    auto g_xyyyyy_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 331);

    auto g_xyyyyy_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 332);

    auto g_xyyyyy_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 333);

    auto g_xyyyyy_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 334);

    auto g_xyyyyy_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 335);

    #pragma omp simd aligned(g_xyyyyy_xxxxx_0, g_xyyyyy_xxxxy_0, g_xyyyyy_xxxxz_0, g_xyyyyy_xxxyy_0, g_xyyyyy_xxxyz_0, g_xyyyyy_xxxzz_0, g_xyyyyy_xxyyy_0, g_xyyyyy_xxyyz_0, g_xyyyyy_xxyzz_0, g_xyyyyy_xxzzz_0, g_xyyyyy_xyyyy_0, g_xyyyyy_xyyyz_0, g_xyyyyy_xyyzz_0, g_xyyyyy_xyzzz_0, g_xyyyyy_xzzzz_0, g_xyyyyy_yyyyy_0, g_xyyyyy_yyyyz_0, g_xyyyyy_yyyzz_0, g_xyyyyy_yyzzz_0, g_xyyyyy_yzzzz_0, g_xyyyyy_zzzzz_0, g_yyyyy_xxxx_1, g_yyyyy_xxxxx_1, g_yyyyy_xxxxy_1, g_yyyyy_xxxxz_1, g_yyyyy_xxxy_1, g_yyyyy_xxxyy_1, g_yyyyy_xxxyz_1, g_yyyyy_xxxz_1, g_yyyyy_xxxzz_1, g_yyyyy_xxyy_1, g_yyyyy_xxyyy_1, g_yyyyy_xxyyz_1, g_yyyyy_xxyz_1, g_yyyyy_xxyzz_1, g_yyyyy_xxzz_1, g_yyyyy_xxzzz_1, g_yyyyy_xyyy_1, g_yyyyy_xyyyy_1, g_yyyyy_xyyyz_1, g_yyyyy_xyyz_1, g_yyyyy_xyyzz_1, g_yyyyy_xyzz_1, g_yyyyy_xyzzz_1, g_yyyyy_xzzz_1, g_yyyyy_xzzzz_1, g_yyyyy_yyyy_1, g_yyyyy_yyyyy_1, g_yyyyy_yyyyz_1, g_yyyyy_yyyz_1, g_yyyyy_yyyzz_1, g_yyyyy_yyzz_1, g_yyyyy_yyzzz_1, g_yyyyy_yzzz_1, g_yyyyy_yzzzz_1, g_yyyyy_zzzz_1, g_yyyyy_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyy_xxxxx_0[i] = 5.0 * g_yyyyy_xxxx_1[i] * fe_0 + g_yyyyy_xxxxx_1[i] * pa_x[i];

        g_xyyyyy_xxxxy_0[i] = 4.0 * g_yyyyy_xxxy_1[i] * fe_0 + g_yyyyy_xxxxy_1[i] * pa_x[i];

        g_xyyyyy_xxxxz_0[i] = 4.0 * g_yyyyy_xxxz_1[i] * fe_0 + g_yyyyy_xxxxz_1[i] * pa_x[i];

        g_xyyyyy_xxxyy_0[i] = 3.0 * g_yyyyy_xxyy_1[i] * fe_0 + g_yyyyy_xxxyy_1[i] * pa_x[i];

        g_xyyyyy_xxxyz_0[i] = 3.0 * g_yyyyy_xxyz_1[i] * fe_0 + g_yyyyy_xxxyz_1[i] * pa_x[i];

        g_xyyyyy_xxxzz_0[i] = 3.0 * g_yyyyy_xxzz_1[i] * fe_0 + g_yyyyy_xxxzz_1[i] * pa_x[i];

        g_xyyyyy_xxyyy_0[i] = 2.0 * g_yyyyy_xyyy_1[i] * fe_0 + g_yyyyy_xxyyy_1[i] * pa_x[i];

        g_xyyyyy_xxyyz_0[i] = 2.0 * g_yyyyy_xyyz_1[i] * fe_0 + g_yyyyy_xxyyz_1[i] * pa_x[i];

        g_xyyyyy_xxyzz_0[i] = 2.0 * g_yyyyy_xyzz_1[i] * fe_0 + g_yyyyy_xxyzz_1[i] * pa_x[i];

        g_xyyyyy_xxzzz_0[i] = 2.0 * g_yyyyy_xzzz_1[i] * fe_0 + g_yyyyy_xxzzz_1[i] * pa_x[i];

        g_xyyyyy_xyyyy_0[i] = g_yyyyy_yyyy_1[i] * fe_0 + g_yyyyy_xyyyy_1[i] * pa_x[i];

        g_xyyyyy_xyyyz_0[i] = g_yyyyy_yyyz_1[i] * fe_0 + g_yyyyy_xyyyz_1[i] * pa_x[i];

        g_xyyyyy_xyyzz_0[i] = g_yyyyy_yyzz_1[i] * fe_0 + g_yyyyy_xyyzz_1[i] * pa_x[i];

        g_xyyyyy_xyzzz_0[i] = g_yyyyy_yzzz_1[i] * fe_0 + g_yyyyy_xyzzz_1[i] * pa_x[i];

        g_xyyyyy_xzzzz_0[i] = g_yyyyy_zzzz_1[i] * fe_0 + g_yyyyy_xzzzz_1[i] * pa_x[i];

        g_xyyyyy_yyyyy_0[i] = g_yyyyy_yyyyy_1[i] * pa_x[i];

        g_xyyyyy_yyyyz_0[i] = g_yyyyy_yyyyz_1[i] * pa_x[i];

        g_xyyyyy_yyyzz_0[i] = g_yyyyy_yyyzz_1[i] * pa_x[i];

        g_xyyyyy_yyzzz_0[i] = g_yyyyy_yyzzz_1[i] * pa_x[i];

        g_xyyyyy_yzzzz_0[i] = g_yyyyy_yzzzz_1[i] * pa_x[i];

        g_xyyyyy_zzzzz_0[i] = g_yyyyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 336-357 components of targeted buffer : IH

    auto g_xyyyyz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 336);

    auto g_xyyyyz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 337);

    auto g_xyyyyz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 338);

    auto g_xyyyyz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 339);

    auto g_xyyyyz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 340);

    auto g_xyyyyz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 341);

    auto g_xyyyyz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 342);

    auto g_xyyyyz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 343);

    auto g_xyyyyz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 344);

    auto g_xyyyyz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 345);

    auto g_xyyyyz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 346);

    auto g_xyyyyz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 347);

    auto g_xyyyyz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 348);

    auto g_xyyyyz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 349);

    auto g_xyyyyz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 350);

    auto g_xyyyyz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 351);

    auto g_xyyyyz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 352);

    auto g_xyyyyz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 353);

    auto g_xyyyyz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 354);

    auto g_xyyyyz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 355);

    auto g_xyyyyz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 356);

    #pragma omp simd aligned(g_xyyyy_xxxxx_1, g_xyyyy_xxxxy_1, g_xyyyy_xxxyy_1, g_xyyyy_xxyyy_1, g_xyyyy_xyyyy_1, g_xyyyyz_xxxxx_0, g_xyyyyz_xxxxy_0, g_xyyyyz_xxxxz_0, g_xyyyyz_xxxyy_0, g_xyyyyz_xxxyz_0, g_xyyyyz_xxxzz_0, g_xyyyyz_xxyyy_0, g_xyyyyz_xxyyz_0, g_xyyyyz_xxyzz_0, g_xyyyyz_xxzzz_0, g_xyyyyz_xyyyy_0, g_xyyyyz_xyyyz_0, g_xyyyyz_xyyzz_0, g_xyyyyz_xyzzz_0, g_xyyyyz_xzzzz_0, g_xyyyyz_yyyyy_0, g_xyyyyz_yyyyz_0, g_xyyyyz_yyyzz_0, g_xyyyyz_yyzzz_0, g_xyyyyz_yzzzz_0, g_xyyyyz_zzzzz_0, g_yyyyz_xxxxz_1, g_yyyyz_xxxyz_1, g_yyyyz_xxxz_1, g_yyyyz_xxxzz_1, g_yyyyz_xxyyz_1, g_yyyyz_xxyz_1, g_yyyyz_xxyzz_1, g_yyyyz_xxzz_1, g_yyyyz_xxzzz_1, g_yyyyz_xyyyz_1, g_yyyyz_xyyz_1, g_yyyyz_xyyzz_1, g_yyyyz_xyzz_1, g_yyyyz_xyzzz_1, g_yyyyz_xzzz_1, g_yyyyz_xzzzz_1, g_yyyyz_yyyyy_1, g_yyyyz_yyyyz_1, g_yyyyz_yyyz_1, g_yyyyz_yyyzz_1, g_yyyyz_yyzz_1, g_yyyyz_yyzzz_1, g_yyyyz_yzzz_1, g_yyyyz_yzzzz_1, g_yyyyz_zzzz_1, g_yyyyz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyz_xxxxx_0[i] = g_xyyyy_xxxxx_1[i] * pa_z[i];

        g_xyyyyz_xxxxy_0[i] = g_xyyyy_xxxxy_1[i] * pa_z[i];

        g_xyyyyz_xxxxz_0[i] = 4.0 * g_yyyyz_xxxz_1[i] * fe_0 + g_yyyyz_xxxxz_1[i] * pa_x[i];

        g_xyyyyz_xxxyy_0[i] = g_xyyyy_xxxyy_1[i] * pa_z[i];

        g_xyyyyz_xxxyz_0[i] = 3.0 * g_yyyyz_xxyz_1[i] * fe_0 + g_yyyyz_xxxyz_1[i] * pa_x[i];

        g_xyyyyz_xxxzz_0[i] = 3.0 * g_yyyyz_xxzz_1[i] * fe_0 + g_yyyyz_xxxzz_1[i] * pa_x[i];

        g_xyyyyz_xxyyy_0[i] = g_xyyyy_xxyyy_1[i] * pa_z[i];

        g_xyyyyz_xxyyz_0[i] = 2.0 * g_yyyyz_xyyz_1[i] * fe_0 + g_yyyyz_xxyyz_1[i] * pa_x[i];

        g_xyyyyz_xxyzz_0[i] = 2.0 * g_yyyyz_xyzz_1[i] * fe_0 + g_yyyyz_xxyzz_1[i] * pa_x[i];

        g_xyyyyz_xxzzz_0[i] = 2.0 * g_yyyyz_xzzz_1[i] * fe_0 + g_yyyyz_xxzzz_1[i] * pa_x[i];

        g_xyyyyz_xyyyy_0[i] = g_xyyyy_xyyyy_1[i] * pa_z[i];

        g_xyyyyz_xyyyz_0[i] = g_yyyyz_yyyz_1[i] * fe_0 + g_yyyyz_xyyyz_1[i] * pa_x[i];

        g_xyyyyz_xyyzz_0[i] = g_yyyyz_yyzz_1[i] * fe_0 + g_yyyyz_xyyzz_1[i] * pa_x[i];

        g_xyyyyz_xyzzz_0[i] = g_yyyyz_yzzz_1[i] * fe_0 + g_yyyyz_xyzzz_1[i] * pa_x[i];

        g_xyyyyz_xzzzz_0[i] = g_yyyyz_zzzz_1[i] * fe_0 + g_yyyyz_xzzzz_1[i] * pa_x[i];

        g_xyyyyz_yyyyy_0[i] = g_yyyyz_yyyyy_1[i] * pa_x[i];

        g_xyyyyz_yyyyz_0[i] = g_yyyyz_yyyyz_1[i] * pa_x[i];

        g_xyyyyz_yyyzz_0[i] = g_yyyyz_yyyzz_1[i] * pa_x[i];

        g_xyyyyz_yyzzz_0[i] = g_yyyyz_yyzzz_1[i] * pa_x[i];

        g_xyyyyz_yzzzz_0[i] = g_yyyyz_yzzzz_1[i] * pa_x[i];

        g_xyyyyz_zzzzz_0[i] = g_yyyyz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 357-378 components of targeted buffer : IH

    auto g_xyyyzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 357);

    auto g_xyyyzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 358);

    auto g_xyyyzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 359);

    auto g_xyyyzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 360);

    auto g_xyyyzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 361);

    auto g_xyyyzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 362);

    auto g_xyyyzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 363);

    auto g_xyyyzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 364);

    auto g_xyyyzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 365);

    auto g_xyyyzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 366);

    auto g_xyyyzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 367);

    auto g_xyyyzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 368);

    auto g_xyyyzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 369);

    auto g_xyyyzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 370);

    auto g_xyyyzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 371);

    auto g_xyyyzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 372);

    auto g_xyyyzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 373);

    auto g_xyyyzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 374);

    auto g_xyyyzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 375);

    auto g_xyyyzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 376);

    auto g_xyyyzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 377);

    #pragma omp simd aligned(g_xyyyzz_xxxxx_0, g_xyyyzz_xxxxy_0, g_xyyyzz_xxxxz_0, g_xyyyzz_xxxyy_0, g_xyyyzz_xxxyz_0, g_xyyyzz_xxxzz_0, g_xyyyzz_xxyyy_0, g_xyyyzz_xxyyz_0, g_xyyyzz_xxyzz_0, g_xyyyzz_xxzzz_0, g_xyyyzz_xyyyy_0, g_xyyyzz_xyyyz_0, g_xyyyzz_xyyzz_0, g_xyyyzz_xyzzz_0, g_xyyyzz_xzzzz_0, g_xyyyzz_yyyyy_0, g_xyyyzz_yyyyz_0, g_xyyyzz_yyyzz_0, g_xyyyzz_yyzzz_0, g_xyyyzz_yzzzz_0, g_xyyyzz_zzzzz_0, g_yyyzz_xxxx_1, g_yyyzz_xxxxx_1, g_yyyzz_xxxxy_1, g_yyyzz_xxxxz_1, g_yyyzz_xxxy_1, g_yyyzz_xxxyy_1, g_yyyzz_xxxyz_1, g_yyyzz_xxxz_1, g_yyyzz_xxxzz_1, g_yyyzz_xxyy_1, g_yyyzz_xxyyy_1, g_yyyzz_xxyyz_1, g_yyyzz_xxyz_1, g_yyyzz_xxyzz_1, g_yyyzz_xxzz_1, g_yyyzz_xxzzz_1, g_yyyzz_xyyy_1, g_yyyzz_xyyyy_1, g_yyyzz_xyyyz_1, g_yyyzz_xyyz_1, g_yyyzz_xyyzz_1, g_yyyzz_xyzz_1, g_yyyzz_xyzzz_1, g_yyyzz_xzzz_1, g_yyyzz_xzzzz_1, g_yyyzz_yyyy_1, g_yyyzz_yyyyy_1, g_yyyzz_yyyyz_1, g_yyyzz_yyyz_1, g_yyyzz_yyyzz_1, g_yyyzz_yyzz_1, g_yyyzz_yyzzz_1, g_yyyzz_yzzz_1, g_yyyzz_yzzzz_1, g_yyyzz_zzzz_1, g_yyyzz_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyzz_xxxxx_0[i] = 5.0 * g_yyyzz_xxxx_1[i] * fe_0 + g_yyyzz_xxxxx_1[i] * pa_x[i];

        g_xyyyzz_xxxxy_0[i] = 4.0 * g_yyyzz_xxxy_1[i] * fe_0 + g_yyyzz_xxxxy_1[i] * pa_x[i];

        g_xyyyzz_xxxxz_0[i] = 4.0 * g_yyyzz_xxxz_1[i] * fe_0 + g_yyyzz_xxxxz_1[i] * pa_x[i];

        g_xyyyzz_xxxyy_0[i] = 3.0 * g_yyyzz_xxyy_1[i] * fe_0 + g_yyyzz_xxxyy_1[i] * pa_x[i];

        g_xyyyzz_xxxyz_0[i] = 3.0 * g_yyyzz_xxyz_1[i] * fe_0 + g_yyyzz_xxxyz_1[i] * pa_x[i];

        g_xyyyzz_xxxzz_0[i] = 3.0 * g_yyyzz_xxzz_1[i] * fe_0 + g_yyyzz_xxxzz_1[i] * pa_x[i];

        g_xyyyzz_xxyyy_0[i] = 2.0 * g_yyyzz_xyyy_1[i] * fe_0 + g_yyyzz_xxyyy_1[i] * pa_x[i];

        g_xyyyzz_xxyyz_0[i] = 2.0 * g_yyyzz_xyyz_1[i] * fe_0 + g_yyyzz_xxyyz_1[i] * pa_x[i];

        g_xyyyzz_xxyzz_0[i] = 2.0 * g_yyyzz_xyzz_1[i] * fe_0 + g_yyyzz_xxyzz_1[i] * pa_x[i];

        g_xyyyzz_xxzzz_0[i] = 2.0 * g_yyyzz_xzzz_1[i] * fe_0 + g_yyyzz_xxzzz_1[i] * pa_x[i];

        g_xyyyzz_xyyyy_0[i] = g_yyyzz_yyyy_1[i] * fe_0 + g_yyyzz_xyyyy_1[i] * pa_x[i];

        g_xyyyzz_xyyyz_0[i] = g_yyyzz_yyyz_1[i] * fe_0 + g_yyyzz_xyyyz_1[i] * pa_x[i];

        g_xyyyzz_xyyzz_0[i] = g_yyyzz_yyzz_1[i] * fe_0 + g_yyyzz_xyyzz_1[i] * pa_x[i];

        g_xyyyzz_xyzzz_0[i] = g_yyyzz_yzzz_1[i] * fe_0 + g_yyyzz_xyzzz_1[i] * pa_x[i];

        g_xyyyzz_xzzzz_0[i] = g_yyyzz_zzzz_1[i] * fe_0 + g_yyyzz_xzzzz_1[i] * pa_x[i];

        g_xyyyzz_yyyyy_0[i] = g_yyyzz_yyyyy_1[i] * pa_x[i];

        g_xyyyzz_yyyyz_0[i] = g_yyyzz_yyyyz_1[i] * pa_x[i];

        g_xyyyzz_yyyzz_0[i] = g_yyyzz_yyyzz_1[i] * pa_x[i];

        g_xyyyzz_yyzzz_0[i] = g_yyyzz_yyzzz_1[i] * pa_x[i];

        g_xyyyzz_yzzzz_0[i] = g_yyyzz_yzzzz_1[i] * pa_x[i];

        g_xyyyzz_zzzzz_0[i] = g_yyyzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 378-399 components of targeted buffer : IH

    auto g_xyyzzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 378);

    auto g_xyyzzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 379);

    auto g_xyyzzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 380);

    auto g_xyyzzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 381);

    auto g_xyyzzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 382);

    auto g_xyyzzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 383);

    auto g_xyyzzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 384);

    auto g_xyyzzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 385);

    auto g_xyyzzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 386);

    auto g_xyyzzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 387);

    auto g_xyyzzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 388);

    auto g_xyyzzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 389);

    auto g_xyyzzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 390);

    auto g_xyyzzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 391);

    auto g_xyyzzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 392);

    auto g_xyyzzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 393);

    auto g_xyyzzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 394);

    auto g_xyyzzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 395);

    auto g_xyyzzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 396);

    auto g_xyyzzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 397);

    auto g_xyyzzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 398);

    #pragma omp simd aligned(g_xyyzzz_xxxxx_0, g_xyyzzz_xxxxy_0, g_xyyzzz_xxxxz_0, g_xyyzzz_xxxyy_0, g_xyyzzz_xxxyz_0, g_xyyzzz_xxxzz_0, g_xyyzzz_xxyyy_0, g_xyyzzz_xxyyz_0, g_xyyzzz_xxyzz_0, g_xyyzzz_xxzzz_0, g_xyyzzz_xyyyy_0, g_xyyzzz_xyyyz_0, g_xyyzzz_xyyzz_0, g_xyyzzz_xyzzz_0, g_xyyzzz_xzzzz_0, g_xyyzzz_yyyyy_0, g_xyyzzz_yyyyz_0, g_xyyzzz_yyyzz_0, g_xyyzzz_yyzzz_0, g_xyyzzz_yzzzz_0, g_xyyzzz_zzzzz_0, g_yyzzz_xxxx_1, g_yyzzz_xxxxx_1, g_yyzzz_xxxxy_1, g_yyzzz_xxxxz_1, g_yyzzz_xxxy_1, g_yyzzz_xxxyy_1, g_yyzzz_xxxyz_1, g_yyzzz_xxxz_1, g_yyzzz_xxxzz_1, g_yyzzz_xxyy_1, g_yyzzz_xxyyy_1, g_yyzzz_xxyyz_1, g_yyzzz_xxyz_1, g_yyzzz_xxyzz_1, g_yyzzz_xxzz_1, g_yyzzz_xxzzz_1, g_yyzzz_xyyy_1, g_yyzzz_xyyyy_1, g_yyzzz_xyyyz_1, g_yyzzz_xyyz_1, g_yyzzz_xyyzz_1, g_yyzzz_xyzz_1, g_yyzzz_xyzzz_1, g_yyzzz_xzzz_1, g_yyzzz_xzzzz_1, g_yyzzz_yyyy_1, g_yyzzz_yyyyy_1, g_yyzzz_yyyyz_1, g_yyzzz_yyyz_1, g_yyzzz_yyyzz_1, g_yyzzz_yyzz_1, g_yyzzz_yyzzz_1, g_yyzzz_yzzz_1, g_yyzzz_yzzzz_1, g_yyzzz_zzzz_1, g_yyzzz_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzzz_xxxxx_0[i] = 5.0 * g_yyzzz_xxxx_1[i] * fe_0 + g_yyzzz_xxxxx_1[i] * pa_x[i];

        g_xyyzzz_xxxxy_0[i] = 4.0 * g_yyzzz_xxxy_1[i] * fe_0 + g_yyzzz_xxxxy_1[i] * pa_x[i];

        g_xyyzzz_xxxxz_0[i] = 4.0 * g_yyzzz_xxxz_1[i] * fe_0 + g_yyzzz_xxxxz_1[i] * pa_x[i];

        g_xyyzzz_xxxyy_0[i] = 3.0 * g_yyzzz_xxyy_1[i] * fe_0 + g_yyzzz_xxxyy_1[i] * pa_x[i];

        g_xyyzzz_xxxyz_0[i] = 3.0 * g_yyzzz_xxyz_1[i] * fe_0 + g_yyzzz_xxxyz_1[i] * pa_x[i];

        g_xyyzzz_xxxzz_0[i] = 3.0 * g_yyzzz_xxzz_1[i] * fe_0 + g_yyzzz_xxxzz_1[i] * pa_x[i];

        g_xyyzzz_xxyyy_0[i] = 2.0 * g_yyzzz_xyyy_1[i] * fe_0 + g_yyzzz_xxyyy_1[i] * pa_x[i];

        g_xyyzzz_xxyyz_0[i] = 2.0 * g_yyzzz_xyyz_1[i] * fe_0 + g_yyzzz_xxyyz_1[i] * pa_x[i];

        g_xyyzzz_xxyzz_0[i] = 2.0 * g_yyzzz_xyzz_1[i] * fe_0 + g_yyzzz_xxyzz_1[i] * pa_x[i];

        g_xyyzzz_xxzzz_0[i] = 2.0 * g_yyzzz_xzzz_1[i] * fe_0 + g_yyzzz_xxzzz_1[i] * pa_x[i];

        g_xyyzzz_xyyyy_0[i] = g_yyzzz_yyyy_1[i] * fe_0 + g_yyzzz_xyyyy_1[i] * pa_x[i];

        g_xyyzzz_xyyyz_0[i] = g_yyzzz_yyyz_1[i] * fe_0 + g_yyzzz_xyyyz_1[i] * pa_x[i];

        g_xyyzzz_xyyzz_0[i] = g_yyzzz_yyzz_1[i] * fe_0 + g_yyzzz_xyyzz_1[i] * pa_x[i];

        g_xyyzzz_xyzzz_0[i] = g_yyzzz_yzzz_1[i] * fe_0 + g_yyzzz_xyzzz_1[i] * pa_x[i];

        g_xyyzzz_xzzzz_0[i] = g_yyzzz_zzzz_1[i] * fe_0 + g_yyzzz_xzzzz_1[i] * pa_x[i];

        g_xyyzzz_yyyyy_0[i] = g_yyzzz_yyyyy_1[i] * pa_x[i];

        g_xyyzzz_yyyyz_0[i] = g_yyzzz_yyyyz_1[i] * pa_x[i];

        g_xyyzzz_yyyzz_0[i] = g_yyzzz_yyyzz_1[i] * pa_x[i];

        g_xyyzzz_yyzzz_0[i] = g_yyzzz_yyzzz_1[i] * pa_x[i];

        g_xyyzzz_yzzzz_0[i] = g_yyzzz_yzzzz_1[i] * pa_x[i];

        g_xyyzzz_zzzzz_0[i] = g_yyzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 399-420 components of targeted buffer : IH

    auto g_xyzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 399);

    auto g_xyzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 400);

    auto g_xyzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 401);

    auto g_xyzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 402);

    auto g_xyzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 403);

    auto g_xyzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 404);

    auto g_xyzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 405);

    auto g_xyzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 406);

    auto g_xyzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 407);

    auto g_xyzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 408);

    auto g_xyzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 409);

    auto g_xyzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 410);

    auto g_xyzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 411);

    auto g_xyzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 412);

    auto g_xyzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 413);

    auto g_xyzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 414);

    auto g_xyzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 415);

    auto g_xyzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 416);

    auto g_xyzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 417);

    auto g_xyzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 418);

    auto g_xyzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 419);

    #pragma omp simd aligned(g_xyzzzz_xxxxx_0, g_xyzzzz_xxxxy_0, g_xyzzzz_xxxxz_0, g_xyzzzz_xxxyy_0, g_xyzzzz_xxxyz_0, g_xyzzzz_xxxzz_0, g_xyzzzz_xxyyy_0, g_xyzzzz_xxyyz_0, g_xyzzzz_xxyzz_0, g_xyzzzz_xxzzz_0, g_xyzzzz_xyyyy_0, g_xyzzzz_xyyyz_0, g_xyzzzz_xyyzz_0, g_xyzzzz_xyzzz_0, g_xyzzzz_xzzzz_0, g_xyzzzz_yyyyy_0, g_xyzzzz_yyyyz_0, g_xyzzzz_yyyzz_0, g_xyzzzz_yyzzz_0, g_xyzzzz_yzzzz_0, g_xyzzzz_zzzzz_0, g_xzzzz_xxxxx_1, g_xzzzz_xxxxz_1, g_xzzzz_xxxzz_1, g_xzzzz_xxzzz_1, g_xzzzz_xzzzz_1, g_yzzzz_xxxxy_1, g_yzzzz_xxxy_1, g_yzzzz_xxxyy_1, g_yzzzz_xxxyz_1, g_yzzzz_xxyy_1, g_yzzzz_xxyyy_1, g_yzzzz_xxyyz_1, g_yzzzz_xxyz_1, g_yzzzz_xxyzz_1, g_yzzzz_xyyy_1, g_yzzzz_xyyyy_1, g_yzzzz_xyyyz_1, g_yzzzz_xyyz_1, g_yzzzz_xyyzz_1, g_yzzzz_xyzz_1, g_yzzzz_xyzzz_1, g_yzzzz_yyyy_1, g_yzzzz_yyyyy_1, g_yzzzz_yyyyz_1, g_yzzzz_yyyz_1, g_yzzzz_yyyzz_1, g_yzzzz_yyzz_1, g_yzzzz_yyzzz_1, g_yzzzz_yzzz_1, g_yzzzz_yzzzz_1, g_yzzzz_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzzz_xxxxx_0[i] = g_xzzzz_xxxxx_1[i] * pa_y[i];

        g_xyzzzz_xxxxy_0[i] = 4.0 * g_yzzzz_xxxy_1[i] * fe_0 + g_yzzzz_xxxxy_1[i] * pa_x[i];

        g_xyzzzz_xxxxz_0[i] = g_xzzzz_xxxxz_1[i] * pa_y[i];

        g_xyzzzz_xxxyy_0[i] = 3.0 * g_yzzzz_xxyy_1[i] * fe_0 + g_yzzzz_xxxyy_1[i] * pa_x[i];

        g_xyzzzz_xxxyz_0[i] = 3.0 * g_yzzzz_xxyz_1[i] * fe_0 + g_yzzzz_xxxyz_1[i] * pa_x[i];

        g_xyzzzz_xxxzz_0[i] = g_xzzzz_xxxzz_1[i] * pa_y[i];

        g_xyzzzz_xxyyy_0[i] = 2.0 * g_yzzzz_xyyy_1[i] * fe_0 + g_yzzzz_xxyyy_1[i] * pa_x[i];

        g_xyzzzz_xxyyz_0[i] = 2.0 * g_yzzzz_xyyz_1[i] * fe_0 + g_yzzzz_xxyyz_1[i] * pa_x[i];

        g_xyzzzz_xxyzz_0[i] = 2.0 * g_yzzzz_xyzz_1[i] * fe_0 + g_yzzzz_xxyzz_1[i] * pa_x[i];

        g_xyzzzz_xxzzz_0[i] = g_xzzzz_xxzzz_1[i] * pa_y[i];

        g_xyzzzz_xyyyy_0[i] = g_yzzzz_yyyy_1[i] * fe_0 + g_yzzzz_xyyyy_1[i] * pa_x[i];

        g_xyzzzz_xyyyz_0[i] = g_yzzzz_yyyz_1[i] * fe_0 + g_yzzzz_xyyyz_1[i] * pa_x[i];

        g_xyzzzz_xyyzz_0[i] = g_yzzzz_yyzz_1[i] * fe_0 + g_yzzzz_xyyzz_1[i] * pa_x[i];

        g_xyzzzz_xyzzz_0[i] = g_yzzzz_yzzz_1[i] * fe_0 + g_yzzzz_xyzzz_1[i] * pa_x[i];

        g_xyzzzz_xzzzz_0[i] = g_xzzzz_xzzzz_1[i] * pa_y[i];

        g_xyzzzz_yyyyy_0[i] = g_yzzzz_yyyyy_1[i] * pa_x[i];

        g_xyzzzz_yyyyz_0[i] = g_yzzzz_yyyyz_1[i] * pa_x[i];

        g_xyzzzz_yyyzz_0[i] = g_yzzzz_yyyzz_1[i] * pa_x[i];

        g_xyzzzz_yyzzz_0[i] = g_yzzzz_yyzzz_1[i] * pa_x[i];

        g_xyzzzz_yzzzz_0[i] = g_yzzzz_yzzzz_1[i] * pa_x[i];

        g_xyzzzz_zzzzz_0[i] = g_yzzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 420-441 components of targeted buffer : IH

    auto g_xzzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 420);

    auto g_xzzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 421);

    auto g_xzzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 422);

    auto g_xzzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 423);

    auto g_xzzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 424);

    auto g_xzzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 425);

    auto g_xzzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 426);

    auto g_xzzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 427);

    auto g_xzzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 428);

    auto g_xzzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 429);

    auto g_xzzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 430);

    auto g_xzzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 431);

    auto g_xzzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 432);

    auto g_xzzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 433);

    auto g_xzzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 434);

    auto g_xzzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 435);

    auto g_xzzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 436);

    auto g_xzzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 437);

    auto g_xzzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 438);

    auto g_xzzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 439);

    auto g_xzzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 440);

    #pragma omp simd aligned(g_xzzzzz_xxxxx_0, g_xzzzzz_xxxxy_0, g_xzzzzz_xxxxz_0, g_xzzzzz_xxxyy_0, g_xzzzzz_xxxyz_0, g_xzzzzz_xxxzz_0, g_xzzzzz_xxyyy_0, g_xzzzzz_xxyyz_0, g_xzzzzz_xxyzz_0, g_xzzzzz_xxzzz_0, g_xzzzzz_xyyyy_0, g_xzzzzz_xyyyz_0, g_xzzzzz_xyyzz_0, g_xzzzzz_xyzzz_0, g_xzzzzz_xzzzz_0, g_xzzzzz_yyyyy_0, g_xzzzzz_yyyyz_0, g_xzzzzz_yyyzz_0, g_xzzzzz_yyzzz_0, g_xzzzzz_yzzzz_0, g_xzzzzz_zzzzz_0, g_zzzzz_xxxx_1, g_zzzzz_xxxxx_1, g_zzzzz_xxxxy_1, g_zzzzz_xxxxz_1, g_zzzzz_xxxy_1, g_zzzzz_xxxyy_1, g_zzzzz_xxxyz_1, g_zzzzz_xxxz_1, g_zzzzz_xxxzz_1, g_zzzzz_xxyy_1, g_zzzzz_xxyyy_1, g_zzzzz_xxyyz_1, g_zzzzz_xxyz_1, g_zzzzz_xxyzz_1, g_zzzzz_xxzz_1, g_zzzzz_xxzzz_1, g_zzzzz_xyyy_1, g_zzzzz_xyyyy_1, g_zzzzz_xyyyz_1, g_zzzzz_xyyz_1, g_zzzzz_xyyzz_1, g_zzzzz_xyzz_1, g_zzzzz_xyzzz_1, g_zzzzz_xzzz_1, g_zzzzz_xzzzz_1, g_zzzzz_yyyy_1, g_zzzzz_yyyyy_1, g_zzzzz_yyyyz_1, g_zzzzz_yyyz_1, g_zzzzz_yyyzz_1, g_zzzzz_yyzz_1, g_zzzzz_yyzzz_1, g_zzzzz_yzzz_1, g_zzzzz_yzzzz_1, g_zzzzz_zzzz_1, g_zzzzz_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzzz_xxxxx_0[i] = 5.0 * g_zzzzz_xxxx_1[i] * fe_0 + g_zzzzz_xxxxx_1[i] * pa_x[i];

        g_xzzzzz_xxxxy_0[i] = 4.0 * g_zzzzz_xxxy_1[i] * fe_0 + g_zzzzz_xxxxy_1[i] * pa_x[i];

        g_xzzzzz_xxxxz_0[i] = 4.0 * g_zzzzz_xxxz_1[i] * fe_0 + g_zzzzz_xxxxz_1[i] * pa_x[i];

        g_xzzzzz_xxxyy_0[i] = 3.0 * g_zzzzz_xxyy_1[i] * fe_0 + g_zzzzz_xxxyy_1[i] * pa_x[i];

        g_xzzzzz_xxxyz_0[i] = 3.0 * g_zzzzz_xxyz_1[i] * fe_0 + g_zzzzz_xxxyz_1[i] * pa_x[i];

        g_xzzzzz_xxxzz_0[i] = 3.0 * g_zzzzz_xxzz_1[i] * fe_0 + g_zzzzz_xxxzz_1[i] * pa_x[i];

        g_xzzzzz_xxyyy_0[i] = 2.0 * g_zzzzz_xyyy_1[i] * fe_0 + g_zzzzz_xxyyy_1[i] * pa_x[i];

        g_xzzzzz_xxyyz_0[i] = 2.0 * g_zzzzz_xyyz_1[i] * fe_0 + g_zzzzz_xxyyz_1[i] * pa_x[i];

        g_xzzzzz_xxyzz_0[i] = 2.0 * g_zzzzz_xyzz_1[i] * fe_0 + g_zzzzz_xxyzz_1[i] * pa_x[i];

        g_xzzzzz_xxzzz_0[i] = 2.0 * g_zzzzz_xzzz_1[i] * fe_0 + g_zzzzz_xxzzz_1[i] * pa_x[i];

        g_xzzzzz_xyyyy_0[i] = g_zzzzz_yyyy_1[i] * fe_0 + g_zzzzz_xyyyy_1[i] * pa_x[i];

        g_xzzzzz_xyyyz_0[i] = g_zzzzz_yyyz_1[i] * fe_0 + g_zzzzz_xyyyz_1[i] * pa_x[i];

        g_xzzzzz_xyyzz_0[i] = g_zzzzz_yyzz_1[i] * fe_0 + g_zzzzz_xyyzz_1[i] * pa_x[i];

        g_xzzzzz_xyzzz_0[i] = g_zzzzz_yzzz_1[i] * fe_0 + g_zzzzz_xyzzz_1[i] * pa_x[i];

        g_xzzzzz_xzzzz_0[i] = g_zzzzz_zzzz_1[i] * fe_0 + g_zzzzz_xzzzz_1[i] * pa_x[i];

        g_xzzzzz_yyyyy_0[i] = g_zzzzz_yyyyy_1[i] * pa_x[i];

        g_xzzzzz_yyyyz_0[i] = g_zzzzz_yyyyz_1[i] * pa_x[i];

        g_xzzzzz_yyyzz_0[i] = g_zzzzz_yyyzz_1[i] * pa_x[i];

        g_xzzzzz_yyzzz_0[i] = g_zzzzz_yyzzz_1[i] * pa_x[i];

        g_xzzzzz_yzzzz_0[i] = g_zzzzz_yzzzz_1[i] * pa_x[i];

        g_xzzzzz_zzzzz_0[i] = g_zzzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 441-462 components of targeted buffer : IH

    auto g_yyyyyy_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 441);

    auto g_yyyyyy_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 442);

    auto g_yyyyyy_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 443);

    auto g_yyyyyy_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 444);

    auto g_yyyyyy_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 445);

    auto g_yyyyyy_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 446);

    auto g_yyyyyy_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 447);

    auto g_yyyyyy_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 448);

    auto g_yyyyyy_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 449);

    auto g_yyyyyy_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 450);

    auto g_yyyyyy_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 451);

    auto g_yyyyyy_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 452);

    auto g_yyyyyy_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 453);

    auto g_yyyyyy_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 454);

    auto g_yyyyyy_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 455);

    auto g_yyyyyy_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 456);

    auto g_yyyyyy_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 457);

    auto g_yyyyyy_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 458);

    auto g_yyyyyy_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 459);

    auto g_yyyyyy_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 460);

    auto g_yyyyyy_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 461);

    #pragma omp simd aligned(g_yyyy_xxxxx_0, g_yyyy_xxxxx_1, g_yyyy_xxxxy_0, g_yyyy_xxxxy_1, g_yyyy_xxxxz_0, g_yyyy_xxxxz_1, g_yyyy_xxxyy_0, g_yyyy_xxxyy_1, g_yyyy_xxxyz_0, g_yyyy_xxxyz_1, g_yyyy_xxxzz_0, g_yyyy_xxxzz_1, g_yyyy_xxyyy_0, g_yyyy_xxyyy_1, g_yyyy_xxyyz_0, g_yyyy_xxyyz_1, g_yyyy_xxyzz_0, g_yyyy_xxyzz_1, g_yyyy_xxzzz_0, g_yyyy_xxzzz_1, g_yyyy_xyyyy_0, g_yyyy_xyyyy_1, g_yyyy_xyyyz_0, g_yyyy_xyyyz_1, g_yyyy_xyyzz_0, g_yyyy_xyyzz_1, g_yyyy_xyzzz_0, g_yyyy_xyzzz_1, g_yyyy_xzzzz_0, g_yyyy_xzzzz_1, g_yyyy_yyyyy_0, g_yyyy_yyyyy_1, g_yyyy_yyyyz_0, g_yyyy_yyyyz_1, g_yyyy_yyyzz_0, g_yyyy_yyyzz_1, g_yyyy_yyzzz_0, g_yyyy_yyzzz_1, g_yyyy_yzzzz_0, g_yyyy_yzzzz_1, g_yyyy_zzzzz_0, g_yyyy_zzzzz_1, g_yyyyy_xxxx_1, g_yyyyy_xxxxx_1, g_yyyyy_xxxxy_1, g_yyyyy_xxxxz_1, g_yyyyy_xxxy_1, g_yyyyy_xxxyy_1, g_yyyyy_xxxyz_1, g_yyyyy_xxxz_1, g_yyyyy_xxxzz_1, g_yyyyy_xxyy_1, g_yyyyy_xxyyy_1, g_yyyyy_xxyyz_1, g_yyyyy_xxyz_1, g_yyyyy_xxyzz_1, g_yyyyy_xxzz_1, g_yyyyy_xxzzz_1, g_yyyyy_xyyy_1, g_yyyyy_xyyyy_1, g_yyyyy_xyyyz_1, g_yyyyy_xyyz_1, g_yyyyy_xyyzz_1, g_yyyyy_xyzz_1, g_yyyyy_xyzzz_1, g_yyyyy_xzzz_1, g_yyyyy_xzzzz_1, g_yyyyy_yyyy_1, g_yyyyy_yyyyy_1, g_yyyyy_yyyyz_1, g_yyyyy_yyyz_1, g_yyyyy_yyyzz_1, g_yyyyy_yyzz_1, g_yyyyy_yyzzz_1, g_yyyyy_yzzz_1, g_yyyyy_yzzzz_1, g_yyyyy_zzzz_1, g_yyyyy_zzzzz_1, g_yyyyyy_xxxxx_0, g_yyyyyy_xxxxy_0, g_yyyyyy_xxxxz_0, g_yyyyyy_xxxyy_0, g_yyyyyy_xxxyz_0, g_yyyyyy_xxxzz_0, g_yyyyyy_xxyyy_0, g_yyyyyy_xxyyz_0, g_yyyyyy_xxyzz_0, g_yyyyyy_xxzzz_0, g_yyyyyy_xyyyy_0, g_yyyyyy_xyyyz_0, g_yyyyyy_xyyzz_0, g_yyyyyy_xyzzz_0, g_yyyyyy_xzzzz_0, g_yyyyyy_yyyyy_0, g_yyyyyy_yyyyz_0, g_yyyyyy_yyyzz_0, g_yyyyyy_yyzzz_0, g_yyyyyy_yzzzz_0, g_yyyyyy_zzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyy_xxxxx_0[i] = 5.0 * g_yyyy_xxxxx_0[i] * fbe_0 - 5.0 * g_yyyy_xxxxx_1[i] * fz_be_0 + g_yyyyy_xxxxx_1[i] * pa_y[i];

        g_yyyyyy_xxxxy_0[i] = 5.0 * g_yyyy_xxxxy_0[i] * fbe_0 - 5.0 * g_yyyy_xxxxy_1[i] * fz_be_0 + g_yyyyy_xxxx_1[i] * fe_0 + g_yyyyy_xxxxy_1[i] * pa_y[i];

        g_yyyyyy_xxxxz_0[i] = 5.0 * g_yyyy_xxxxz_0[i] * fbe_0 - 5.0 * g_yyyy_xxxxz_1[i] * fz_be_0 + g_yyyyy_xxxxz_1[i] * pa_y[i];

        g_yyyyyy_xxxyy_0[i] = 5.0 * g_yyyy_xxxyy_0[i] * fbe_0 - 5.0 * g_yyyy_xxxyy_1[i] * fz_be_0 + 2.0 * g_yyyyy_xxxy_1[i] * fe_0 + g_yyyyy_xxxyy_1[i] * pa_y[i];

        g_yyyyyy_xxxyz_0[i] = 5.0 * g_yyyy_xxxyz_0[i] * fbe_0 - 5.0 * g_yyyy_xxxyz_1[i] * fz_be_0 + g_yyyyy_xxxz_1[i] * fe_0 + g_yyyyy_xxxyz_1[i] * pa_y[i];

        g_yyyyyy_xxxzz_0[i] = 5.0 * g_yyyy_xxxzz_0[i] * fbe_0 - 5.0 * g_yyyy_xxxzz_1[i] * fz_be_0 + g_yyyyy_xxxzz_1[i] * pa_y[i];

        g_yyyyyy_xxyyy_0[i] = 5.0 * g_yyyy_xxyyy_0[i] * fbe_0 - 5.0 * g_yyyy_xxyyy_1[i] * fz_be_0 + 3.0 * g_yyyyy_xxyy_1[i] * fe_0 + g_yyyyy_xxyyy_1[i] * pa_y[i];

        g_yyyyyy_xxyyz_0[i] = 5.0 * g_yyyy_xxyyz_0[i] * fbe_0 - 5.0 * g_yyyy_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyy_xxyz_1[i] * fe_0 + g_yyyyy_xxyyz_1[i] * pa_y[i];

        g_yyyyyy_xxyzz_0[i] = 5.0 * g_yyyy_xxyzz_0[i] * fbe_0 - 5.0 * g_yyyy_xxyzz_1[i] * fz_be_0 + g_yyyyy_xxzz_1[i] * fe_0 + g_yyyyy_xxyzz_1[i] * pa_y[i];

        g_yyyyyy_xxzzz_0[i] = 5.0 * g_yyyy_xxzzz_0[i] * fbe_0 - 5.0 * g_yyyy_xxzzz_1[i] * fz_be_0 + g_yyyyy_xxzzz_1[i] * pa_y[i];

        g_yyyyyy_xyyyy_0[i] = 5.0 * g_yyyy_xyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_xyyyy_1[i] * fz_be_0 + 4.0 * g_yyyyy_xyyy_1[i] * fe_0 + g_yyyyy_xyyyy_1[i] * pa_y[i];

        g_yyyyyy_xyyyz_0[i] = 5.0 * g_yyyy_xyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyy_xyyz_1[i] * fe_0 + g_yyyyy_xyyyz_1[i] * pa_y[i];

        g_yyyyyy_xyyzz_0[i] = 5.0 * g_yyyy_xyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_xyzz_1[i] * fe_0 + g_yyyyy_xyyzz_1[i] * pa_y[i];

        g_yyyyyy_xyzzz_0[i] = 5.0 * g_yyyy_xyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_xyzzz_1[i] * fz_be_0 + g_yyyyy_xzzz_1[i] * fe_0 + g_yyyyy_xyzzz_1[i] * pa_y[i];

        g_yyyyyy_xzzzz_0[i] = 5.0 * g_yyyy_xzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_xzzzz_1[i] * fz_be_0 + g_yyyyy_xzzzz_1[i] * pa_y[i];

        g_yyyyyy_yyyyy_0[i] = 5.0 * g_yyyy_yyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_yyyyy_1[i] * fz_be_0 + 5.0 * g_yyyyy_yyyy_1[i] * fe_0 + g_yyyyy_yyyyy_1[i] * pa_y[i];

        g_yyyyyy_yyyyz_0[i] = 5.0 * g_yyyy_yyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyy_yyyz_1[i] * fe_0 + g_yyyyy_yyyyz_1[i] * pa_y[i];

        g_yyyyyy_yyyzz_0[i] = 5.0 * g_yyyy_yyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_yyzz_1[i] * fe_0 + g_yyyyy_yyyzz_1[i] * pa_y[i];

        g_yyyyyy_yyzzz_0[i] = 5.0 * g_yyyy_yyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_yzzz_1[i] * fe_0 + g_yyyyy_yyzzz_1[i] * pa_y[i];

        g_yyyyyy_yzzzz_0[i] = 5.0 * g_yyyy_yzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_yzzzz_1[i] * fz_be_0 + g_yyyyy_zzzz_1[i] * fe_0 + g_yyyyy_yzzzz_1[i] * pa_y[i];

        g_yyyyyy_zzzzz_0[i] = 5.0 * g_yyyy_zzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_zzzzz_1[i] * fz_be_0 + g_yyyyy_zzzzz_1[i] * pa_y[i];
    }

    // Set up 462-483 components of targeted buffer : IH

    auto g_yyyyyz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 462);

    auto g_yyyyyz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 463);

    auto g_yyyyyz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 464);

    auto g_yyyyyz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 465);

    auto g_yyyyyz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 466);

    auto g_yyyyyz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 467);

    auto g_yyyyyz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 468);

    auto g_yyyyyz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 469);

    auto g_yyyyyz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 470);

    auto g_yyyyyz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 471);

    auto g_yyyyyz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 472);

    auto g_yyyyyz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 473);

    auto g_yyyyyz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 474);

    auto g_yyyyyz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 475);

    auto g_yyyyyz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 476);

    auto g_yyyyyz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 477);

    auto g_yyyyyz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 478);

    auto g_yyyyyz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 479);

    auto g_yyyyyz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 480);

    auto g_yyyyyz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 481);

    auto g_yyyyyz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 482);

    #pragma omp simd aligned(g_yyyyy_xxxx_1, g_yyyyy_xxxxx_1, g_yyyyy_xxxxy_1, g_yyyyy_xxxxz_1, g_yyyyy_xxxy_1, g_yyyyy_xxxyy_1, g_yyyyy_xxxyz_1, g_yyyyy_xxxz_1, g_yyyyy_xxxzz_1, g_yyyyy_xxyy_1, g_yyyyy_xxyyy_1, g_yyyyy_xxyyz_1, g_yyyyy_xxyz_1, g_yyyyy_xxyzz_1, g_yyyyy_xxzz_1, g_yyyyy_xxzzz_1, g_yyyyy_xyyy_1, g_yyyyy_xyyyy_1, g_yyyyy_xyyyz_1, g_yyyyy_xyyz_1, g_yyyyy_xyyzz_1, g_yyyyy_xyzz_1, g_yyyyy_xyzzz_1, g_yyyyy_xzzz_1, g_yyyyy_xzzzz_1, g_yyyyy_yyyy_1, g_yyyyy_yyyyy_1, g_yyyyy_yyyyz_1, g_yyyyy_yyyz_1, g_yyyyy_yyyzz_1, g_yyyyy_yyzz_1, g_yyyyy_yyzzz_1, g_yyyyy_yzzz_1, g_yyyyy_yzzzz_1, g_yyyyy_zzzz_1, g_yyyyy_zzzzz_1, g_yyyyyz_xxxxx_0, g_yyyyyz_xxxxy_0, g_yyyyyz_xxxxz_0, g_yyyyyz_xxxyy_0, g_yyyyyz_xxxyz_0, g_yyyyyz_xxxzz_0, g_yyyyyz_xxyyy_0, g_yyyyyz_xxyyz_0, g_yyyyyz_xxyzz_0, g_yyyyyz_xxzzz_0, g_yyyyyz_xyyyy_0, g_yyyyyz_xyyyz_0, g_yyyyyz_xyyzz_0, g_yyyyyz_xyzzz_0, g_yyyyyz_xzzzz_0, g_yyyyyz_yyyyy_0, g_yyyyyz_yyyyz_0, g_yyyyyz_yyyzz_0, g_yyyyyz_yyzzz_0, g_yyyyyz_yzzzz_0, g_yyyyyz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyyz_xxxxx_0[i] = g_yyyyy_xxxxx_1[i] * pa_z[i];

        g_yyyyyz_xxxxy_0[i] = g_yyyyy_xxxxy_1[i] * pa_z[i];

        g_yyyyyz_xxxxz_0[i] = g_yyyyy_xxxx_1[i] * fe_0 + g_yyyyy_xxxxz_1[i] * pa_z[i];

        g_yyyyyz_xxxyy_0[i] = g_yyyyy_xxxyy_1[i] * pa_z[i];

        g_yyyyyz_xxxyz_0[i] = g_yyyyy_xxxy_1[i] * fe_0 + g_yyyyy_xxxyz_1[i] * pa_z[i];

        g_yyyyyz_xxxzz_0[i] = 2.0 * g_yyyyy_xxxz_1[i] * fe_0 + g_yyyyy_xxxzz_1[i] * pa_z[i];

        g_yyyyyz_xxyyy_0[i] = g_yyyyy_xxyyy_1[i] * pa_z[i];

        g_yyyyyz_xxyyz_0[i] = g_yyyyy_xxyy_1[i] * fe_0 + g_yyyyy_xxyyz_1[i] * pa_z[i];

        g_yyyyyz_xxyzz_0[i] = 2.0 * g_yyyyy_xxyz_1[i] * fe_0 + g_yyyyy_xxyzz_1[i] * pa_z[i];

        g_yyyyyz_xxzzz_0[i] = 3.0 * g_yyyyy_xxzz_1[i] * fe_0 + g_yyyyy_xxzzz_1[i] * pa_z[i];

        g_yyyyyz_xyyyy_0[i] = g_yyyyy_xyyyy_1[i] * pa_z[i];

        g_yyyyyz_xyyyz_0[i] = g_yyyyy_xyyy_1[i] * fe_0 + g_yyyyy_xyyyz_1[i] * pa_z[i];

        g_yyyyyz_xyyzz_0[i] = 2.0 * g_yyyyy_xyyz_1[i] * fe_0 + g_yyyyy_xyyzz_1[i] * pa_z[i];

        g_yyyyyz_xyzzz_0[i] = 3.0 * g_yyyyy_xyzz_1[i] * fe_0 + g_yyyyy_xyzzz_1[i] * pa_z[i];

        g_yyyyyz_xzzzz_0[i] = 4.0 * g_yyyyy_xzzz_1[i] * fe_0 + g_yyyyy_xzzzz_1[i] * pa_z[i];

        g_yyyyyz_yyyyy_0[i] = g_yyyyy_yyyyy_1[i] * pa_z[i];

        g_yyyyyz_yyyyz_0[i] = g_yyyyy_yyyy_1[i] * fe_0 + g_yyyyy_yyyyz_1[i] * pa_z[i];

        g_yyyyyz_yyyzz_0[i] = 2.0 * g_yyyyy_yyyz_1[i] * fe_0 + g_yyyyy_yyyzz_1[i] * pa_z[i];

        g_yyyyyz_yyzzz_0[i] = 3.0 * g_yyyyy_yyzz_1[i] * fe_0 + g_yyyyy_yyzzz_1[i] * pa_z[i];

        g_yyyyyz_yzzzz_0[i] = 4.0 * g_yyyyy_yzzz_1[i] * fe_0 + g_yyyyy_yzzzz_1[i] * pa_z[i];

        g_yyyyyz_zzzzz_0[i] = 5.0 * g_yyyyy_zzzz_1[i] * fe_0 + g_yyyyy_zzzzz_1[i] * pa_z[i];
    }

    // Set up 483-504 components of targeted buffer : IH

    auto g_yyyyzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 483);

    auto g_yyyyzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 484);

    auto g_yyyyzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 485);

    auto g_yyyyzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 486);

    auto g_yyyyzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 487);

    auto g_yyyyzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 488);

    auto g_yyyyzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 489);

    auto g_yyyyzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 490);

    auto g_yyyyzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 491);

    auto g_yyyyzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 492);

    auto g_yyyyzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 493);

    auto g_yyyyzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 494);

    auto g_yyyyzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 495);

    auto g_yyyyzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 496);

    auto g_yyyyzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 497);

    auto g_yyyyzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 498);

    auto g_yyyyzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 499);

    auto g_yyyyzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 500);

    auto g_yyyyzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 501);

    auto g_yyyyzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 502);

    auto g_yyyyzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 503);

    #pragma omp simd aligned(g_yyyy_xxxxy_0, g_yyyy_xxxxy_1, g_yyyy_xxxyy_0, g_yyyy_xxxyy_1, g_yyyy_xxyyy_0, g_yyyy_xxyyy_1, g_yyyy_xyyyy_0, g_yyyy_xyyyy_1, g_yyyy_yyyyy_0, g_yyyy_yyyyy_1, g_yyyyz_xxxxy_1, g_yyyyz_xxxyy_1, g_yyyyz_xxyyy_1, g_yyyyz_xyyyy_1, g_yyyyz_yyyyy_1, g_yyyyzz_xxxxx_0, g_yyyyzz_xxxxy_0, g_yyyyzz_xxxxz_0, g_yyyyzz_xxxyy_0, g_yyyyzz_xxxyz_0, g_yyyyzz_xxxzz_0, g_yyyyzz_xxyyy_0, g_yyyyzz_xxyyz_0, g_yyyyzz_xxyzz_0, g_yyyyzz_xxzzz_0, g_yyyyzz_xyyyy_0, g_yyyyzz_xyyyz_0, g_yyyyzz_xyyzz_0, g_yyyyzz_xyzzz_0, g_yyyyzz_xzzzz_0, g_yyyyzz_yyyyy_0, g_yyyyzz_yyyyz_0, g_yyyyzz_yyyzz_0, g_yyyyzz_yyzzz_0, g_yyyyzz_yzzzz_0, g_yyyyzz_zzzzz_0, g_yyyzz_xxxxx_1, g_yyyzz_xxxxz_1, g_yyyzz_xxxyz_1, g_yyyzz_xxxz_1, g_yyyzz_xxxzz_1, g_yyyzz_xxyyz_1, g_yyyzz_xxyz_1, g_yyyzz_xxyzz_1, g_yyyzz_xxzz_1, g_yyyzz_xxzzz_1, g_yyyzz_xyyyz_1, g_yyyzz_xyyz_1, g_yyyzz_xyyzz_1, g_yyyzz_xyzz_1, g_yyyzz_xyzzz_1, g_yyyzz_xzzz_1, g_yyyzz_xzzzz_1, g_yyyzz_yyyyz_1, g_yyyzz_yyyz_1, g_yyyzz_yyyzz_1, g_yyyzz_yyzz_1, g_yyyzz_yyzzz_1, g_yyyzz_yzzz_1, g_yyyzz_yzzzz_1, g_yyyzz_zzzz_1, g_yyyzz_zzzzz_1, g_yyzz_xxxxx_0, g_yyzz_xxxxx_1, g_yyzz_xxxxz_0, g_yyzz_xxxxz_1, g_yyzz_xxxyz_0, g_yyzz_xxxyz_1, g_yyzz_xxxzz_0, g_yyzz_xxxzz_1, g_yyzz_xxyyz_0, g_yyzz_xxyyz_1, g_yyzz_xxyzz_0, g_yyzz_xxyzz_1, g_yyzz_xxzzz_0, g_yyzz_xxzzz_1, g_yyzz_xyyyz_0, g_yyzz_xyyyz_1, g_yyzz_xyyzz_0, g_yyzz_xyyzz_1, g_yyzz_xyzzz_0, g_yyzz_xyzzz_1, g_yyzz_xzzzz_0, g_yyzz_xzzzz_1, g_yyzz_yyyyz_0, g_yyzz_yyyyz_1, g_yyzz_yyyzz_0, g_yyzz_yyyzz_1, g_yyzz_yyzzz_0, g_yyzz_yyzzz_1, g_yyzz_yzzzz_0, g_yyzz_yzzzz_1, g_yyzz_zzzzz_0, g_yyzz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyzz_xxxxx_0[i] = 3.0 * g_yyzz_xxxxx_0[i] * fbe_0 - 3.0 * g_yyzz_xxxxx_1[i] * fz_be_0 + g_yyyzz_xxxxx_1[i] * pa_y[i];

        g_yyyyzz_xxxxy_0[i] = g_yyyy_xxxxy_0[i] * fbe_0 - g_yyyy_xxxxy_1[i] * fz_be_0 + g_yyyyz_xxxxy_1[i] * pa_z[i];

        g_yyyyzz_xxxxz_0[i] = 3.0 * g_yyzz_xxxxz_0[i] * fbe_0 - 3.0 * g_yyzz_xxxxz_1[i] * fz_be_0 + g_yyyzz_xxxxz_1[i] * pa_y[i];

        g_yyyyzz_xxxyy_0[i] = g_yyyy_xxxyy_0[i] * fbe_0 - g_yyyy_xxxyy_1[i] * fz_be_0 + g_yyyyz_xxxyy_1[i] * pa_z[i];

        g_yyyyzz_xxxyz_0[i] = 3.0 * g_yyzz_xxxyz_0[i] * fbe_0 - 3.0 * g_yyzz_xxxyz_1[i] * fz_be_0 + g_yyyzz_xxxz_1[i] * fe_0 + g_yyyzz_xxxyz_1[i] * pa_y[i];

        g_yyyyzz_xxxzz_0[i] = 3.0 * g_yyzz_xxxzz_0[i] * fbe_0 - 3.0 * g_yyzz_xxxzz_1[i] * fz_be_0 + g_yyyzz_xxxzz_1[i] * pa_y[i];

        g_yyyyzz_xxyyy_0[i] = g_yyyy_xxyyy_0[i] * fbe_0 - g_yyyy_xxyyy_1[i] * fz_be_0 + g_yyyyz_xxyyy_1[i] * pa_z[i];

        g_yyyyzz_xxyyz_0[i] = 3.0 * g_yyzz_xxyyz_0[i] * fbe_0 - 3.0 * g_yyzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyyzz_xxyz_1[i] * fe_0 + g_yyyzz_xxyyz_1[i] * pa_y[i];

        g_yyyyzz_xxyzz_0[i] = 3.0 * g_yyzz_xxyzz_0[i] * fbe_0 - 3.0 * g_yyzz_xxyzz_1[i] * fz_be_0 + g_yyyzz_xxzz_1[i] * fe_0 + g_yyyzz_xxyzz_1[i] * pa_y[i];

        g_yyyyzz_xxzzz_0[i] = 3.0 * g_yyzz_xxzzz_0[i] * fbe_0 - 3.0 * g_yyzz_xxzzz_1[i] * fz_be_0 + g_yyyzz_xxzzz_1[i] * pa_y[i];

        g_yyyyzz_xyyyy_0[i] = g_yyyy_xyyyy_0[i] * fbe_0 - g_yyyy_xyyyy_1[i] * fz_be_0 + g_yyyyz_xyyyy_1[i] * pa_z[i];

        g_yyyyzz_xyyyz_0[i] = 3.0 * g_yyzz_xyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyyzz_xyyz_1[i] * fe_0 + g_yyyzz_xyyyz_1[i] * pa_y[i];

        g_yyyyzz_xyyzz_0[i] = 3.0 * g_yyzz_xyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_xyzz_1[i] * fe_0 + g_yyyzz_xyyzz_1[i] * pa_y[i];

        g_yyyyzz_xyzzz_0[i] = 3.0 * g_yyzz_xyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_xyzzz_1[i] * fz_be_0 + g_yyyzz_xzzz_1[i] * fe_0 + g_yyyzz_xyzzz_1[i] * pa_y[i];

        g_yyyyzz_xzzzz_0[i] = 3.0 * g_yyzz_xzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_xzzzz_1[i] * fz_be_0 + g_yyyzz_xzzzz_1[i] * pa_y[i];

        g_yyyyzz_yyyyy_0[i] = g_yyyy_yyyyy_0[i] * fbe_0 - g_yyyy_yyyyy_1[i] * fz_be_0 + g_yyyyz_yyyyy_1[i] * pa_z[i];

        g_yyyyzz_yyyyz_0[i] = 3.0 * g_yyzz_yyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyyzz_yyyz_1[i] * fe_0 + g_yyyzz_yyyyz_1[i] * pa_y[i];

        g_yyyyzz_yyyzz_0[i] = 3.0 * g_yyzz_yyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_yyzz_1[i] * fe_0 + g_yyyzz_yyyzz_1[i] * pa_y[i];

        g_yyyyzz_yyzzz_0[i] = 3.0 * g_yyzz_yyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_yzzz_1[i] * fe_0 + g_yyyzz_yyzzz_1[i] * pa_y[i];

        g_yyyyzz_yzzzz_0[i] = 3.0 * g_yyzz_yzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_yzzzz_1[i] * fz_be_0 + g_yyyzz_zzzz_1[i] * fe_0 + g_yyyzz_yzzzz_1[i] * pa_y[i];

        g_yyyyzz_zzzzz_0[i] = 3.0 * g_yyzz_zzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_zzzzz_1[i] * fz_be_0 + g_yyyzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 504-525 components of targeted buffer : IH

    auto g_yyyzzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 504);

    auto g_yyyzzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 505);

    auto g_yyyzzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 506);

    auto g_yyyzzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 507);

    auto g_yyyzzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 508);

    auto g_yyyzzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 509);

    auto g_yyyzzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 510);

    auto g_yyyzzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 511);

    auto g_yyyzzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 512);

    auto g_yyyzzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 513);

    auto g_yyyzzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 514);

    auto g_yyyzzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 515);

    auto g_yyyzzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 516);

    auto g_yyyzzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 517);

    auto g_yyyzzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 518);

    auto g_yyyzzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 519);

    auto g_yyyzzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 520);

    auto g_yyyzzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 521);

    auto g_yyyzzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 522);

    auto g_yyyzzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 523);

    auto g_yyyzzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 524);

    #pragma omp simd aligned(g_yyyz_xxxxy_0, g_yyyz_xxxxy_1, g_yyyz_xxxyy_0, g_yyyz_xxxyy_1, g_yyyz_xxyyy_0, g_yyyz_xxyyy_1, g_yyyz_xyyyy_0, g_yyyz_xyyyy_1, g_yyyz_yyyyy_0, g_yyyz_yyyyy_1, g_yyyzz_xxxxy_1, g_yyyzz_xxxyy_1, g_yyyzz_xxyyy_1, g_yyyzz_xyyyy_1, g_yyyzz_yyyyy_1, g_yyyzzz_xxxxx_0, g_yyyzzz_xxxxy_0, g_yyyzzz_xxxxz_0, g_yyyzzz_xxxyy_0, g_yyyzzz_xxxyz_0, g_yyyzzz_xxxzz_0, g_yyyzzz_xxyyy_0, g_yyyzzz_xxyyz_0, g_yyyzzz_xxyzz_0, g_yyyzzz_xxzzz_0, g_yyyzzz_xyyyy_0, g_yyyzzz_xyyyz_0, g_yyyzzz_xyyzz_0, g_yyyzzz_xyzzz_0, g_yyyzzz_xzzzz_0, g_yyyzzz_yyyyy_0, g_yyyzzz_yyyyz_0, g_yyyzzz_yyyzz_0, g_yyyzzz_yyzzz_0, g_yyyzzz_yzzzz_0, g_yyyzzz_zzzzz_0, g_yyzzz_xxxxx_1, g_yyzzz_xxxxz_1, g_yyzzz_xxxyz_1, g_yyzzz_xxxz_1, g_yyzzz_xxxzz_1, g_yyzzz_xxyyz_1, g_yyzzz_xxyz_1, g_yyzzz_xxyzz_1, g_yyzzz_xxzz_1, g_yyzzz_xxzzz_1, g_yyzzz_xyyyz_1, g_yyzzz_xyyz_1, g_yyzzz_xyyzz_1, g_yyzzz_xyzz_1, g_yyzzz_xyzzz_1, g_yyzzz_xzzz_1, g_yyzzz_xzzzz_1, g_yyzzz_yyyyz_1, g_yyzzz_yyyz_1, g_yyzzz_yyyzz_1, g_yyzzz_yyzz_1, g_yyzzz_yyzzz_1, g_yyzzz_yzzz_1, g_yyzzz_yzzzz_1, g_yyzzz_zzzz_1, g_yyzzz_zzzzz_1, g_yzzz_xxxxx_0, g_yzzz_xxxxx_1, g_yzzz_xxxxz_0, g_yzzz_xxxxz_1, g_yzzz_xxxyz_0, g_yzzz_xxxyz_1, g_yzzz_xxxzz_0, g_yzzz_xxxzz_1, g_yzzz_xxyyz_0, g_yzzz_xxyyz_1, g_yzzz_xxyzz_0, g_yzzz_xxyzz_1, g_yzzz_xxzzz_0, g_yzzz_xxzzz_1, g_yzzz_xyyyz_0, g_yzzz_xyyyz_1, g_yzzz_xyyzz_0, g_yzzz_xyyzz_1, g_yzzz_xyzzz_0, g_yzzz_xyzzz_1, g_yzzz_xzzzz_0, g_yzzz_xzzzz_1, g_yzzz_yyyyz_0, g_yzzz_yyyyz_1, g_yzzz_yyyzz_0, g_yzzz_yyyzz_1, g_yzzz_yyzzz_0, g_yzzz_yyzzz_1, g_yzzz_yzzzz_0, g_yzzz_yzzzz_1, g_yzzz_zzzzz_0, g_yzzz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzzz_xxxxx_0[i] = 2.0 * g_yzzz_xxxxx_0[i] * fbe_0 - 2.0 * g_yzzz_xxxxx_1[i] * fz_be_0 + g_yyzzz_xxxxx_1[i] * pa_y[i];

        g_yyyzzz_xxxxy_0[i] = 2.0 * g_yyyz_xxxxy_0[i] * fbe_0 - 2.0 * g_yyyz_xxxxy_1[i] * fz_be_0 + g_yyyzz_xxxxy_1[i] * pa_z[i];

        g_yyyzzz_xxxxz_0[i] = 2.0 * g_yzzz_xxxxz_0[i] * fbe_0 - 2.0 * g_yzzz_xxxxz_1[i] * fz_be_0 + g_yyzzz_xxxxz_1[i] * pa_y[i];

        g_yyyzzz_xxxyy_0[i] = 2.0 * g_yyyz_xxxyy_0[i] * fbe_0 - 2.0 * g_yyyz_xxxyy_1[i] * fz_be_0 + g_yyyzz_xxxyy_1[i] * pa_z[i];

        g_yyyzzz_xxxyz_0[i] = 2.0 * g_yzzz_xxxyz_0[i] * fbe_0 - 2.0 * g_yzzz_xxxyz_1[i] * fz_be_0 + g_yyzzz_xxxz_1[i] * fe_0 + g_yyzzz_xxxyz_1[i] * pa_y[i];

        g_yyyzzz_xxxzz_0[i] = 2.0 * g_yzzz_xxxzz_0[i] * fbe_0 - 2.0 * g_yzzz_xxxzz_1[i] * fz_be_0 + g_yyzzz_xxxzz_1[i] * pa_y[i];

        g_yyyzzz_xxyyy_0[i] = 2.0 * g_yyyz_xxyyy_0[i] * fbe_0 - 2.0 * g_yyyz_xxyyy_1[i] * fz_be_0 + g_yyyzz_xxyyy_1[i] * pa_z[i];

        g_yyyzzz_xxyyz_0[i] = 2.0 * g_yzzz_xxyyz_0[i] * fbe_0 - 2.0 * g_yzzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyzzz_xxyz_1[i] * fe_0 + g_yyzzz_xxyyz_1[i] * pa_y[i];

        g_yyyzzz_xxyzz_0[i] = 2.0 * g_yzzz_xxyzz_0[i] * fbe_0 - 2.0 * g_yzzz_xxyzz_1[i] * fz_be_0 + g_yyzzz_xxzz_1[i] * fe_0 + g_yyzzz_xxyzz_1[i] * pa_y[i];

        g_yyyzzz_xxzzz_0[i] = 2.0 * g_yzzz_xxzzz_0[i] * fbe_0 - 2.0 * g_yzzz_xxzzz_1[i] * fz_be_0 + g_yyzzz_xxzzz_1[i] * pa_y[i];

        g_yyyzzz_xyyyy_0[i] = 2.0 * g_yyyz_xyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_xyyyy_1[i] * fz_be_0 + g_yyyzz_xyyyy_1[i] * pa_z[i];

        g_yyyzzz_xyyyz_0[i] = 2.0 * g_yzzz_xyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyzzz_xyyz_1[i] * fe_0 + g_yyzzz_xyyyz_1[i] * pa_y[i];

        g_yyyzzz_xyyzz_0[i] = 2.0 * g_yzzz_xyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_xyzz_1[i] * fe_0 + g_yyzzz_xyyzz_1[i] * pa_y[i];

        g_yyyzzz_xyzzz_0[i] = 2.0 * g_yzzz_xyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_xyzzz_1[i] * fz_be_0 + g_yyzzz_xzzz_1[i] * fe_0 + g_yyzzz_xyzzz_1[i] * pa_y[i];

        g_yyyzzz_xzzzz_0[i] = 2.0 * g_yzzz_xzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_xzzzz_1[i] * fz_be_0 + g_yyzzz_xzzzz_1[i] * pa_y[i];

        g_yyyzzz_yyyyy_0[i] = 2.0 * g_yyyz_yyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_yyyyy_1[i] * fz_be_0 + g_yyyzz_yyyyy_1[i] * pa_z[i];

        g_yyyzzz_yyyyz_0[i] = 2.0 * g_yzzz_yyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyzzz_yyyz_1[i] * fe_0 + g_yyzzz_yyyyz_1[i] * pa_y[i];

        g_yyyzzz_yyyzz_0[i] = 2.0 * g_yzzz_yyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_yyzz_1[i] * fe_0 + g_yyzzz_yyyzz_1[i] * pa_y[i];

        g_yyyzzz_yyzzz_0[i] = 2.0 * g_yzzz_yyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_yzzz_1[i] * fe_0 + g_yyzzz_yyzzz_1[i] * pa_y[i];

        g_yyyzzz_yzzzz_0[i] = 2.0 * g_yzzz_yzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_yzzzz_1[i] * fz_be_0 + g_yyzzz_zzzz_1[i] * fe_0 + g_yyzzz_yzzzz_1[i] * pa_y[i];

        g_yyyzzz_zzzzz_0[i] = 2.0 * g_yzzz_zzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_zzzzz_1[i] * fz_be_0 + g_yyzzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 525-546 components of targeted buffer : IH

    auto g_yyzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 525);

    auto g_yyzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 526);

    auto g_yyzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 527);

    auto g_yyzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 528);

    auto g_yyzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 529);

    auto g_yyzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 530);

    auto g_yyzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 531);

    auto g_yyzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 532);

    auto g_yyzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 533);

    auto g_yyzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 534);

    auto g_yyzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 535);

    auto g_yyzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 536);

    auto g_yyzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 537);

    auto g_yyzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 538);

    auto g_yyzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 539);

    auto g_yyzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 540);

    auto g_yyzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 541);

    auto g_yyzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 542);

    auto g_yyzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 543);

    auto g_yyzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 544);

    auto g_yyzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 545);

    #pragma omp simd aligned(g_yyzz_xxxxy_0, g_yyzz_xxxxy_1, g_yyzz_xxxyy_0, g_yyzz_xxxyy_1, g_yyzz_xxyyy_0, g_yyzz_xxyyy_1, g_yyzz_xyyyy_0, g_yyzz_xyyyy_1, g_yyzz_yyyyy_0, g_yyzz_yyyyy_1, g_yyzzz_xxxxy_1, g_yyzzz_xxxyy_1, g_yyzzz_xxyyy_1, g_yyzzz_xyyyy_1, g_yyzzz_yyyyy_1, g_yyzzzz_xxxxx_0, g_yyzzzz_xxxxy_0, g_yyzzzz_xxxxz_0, g_yyzzzz_xxxyy_0, g_yyzzzz_xxxyz_0, g_yyzzzz_xxxzz_0, g_yyzzzz_xxyyy_0, g_yyzzzz_xxyyz_0, g_yyzzzz_xxyzz_0, g_yyzzzz_xxzzz_0, g_yyzzzz_xyyyy_0, g_yyzzzz_xyyyz_0, g_yyzzzz_xyyzz_0, g_yyzzzz_xyzzz_0, g_yyzzzz_xzzzz_0, g_yyzzzz_yyyyy_0, g_yyzzzz_yyyyz_0, g_yyzzzz_yyyzz_0, g_yyzzzz_yyzzz_0, g_yyzzzz_yzzzz_0, g_yyzzzz_zzzzz_0, g_yzzzz_xxxxx_1, g_yzzzz_xxxxz_1, g_yzzzz_xxxyz_1, g_yzzzz_xxxz_1, g_yzzzz_xxxzz_1, g_yzzzz_xxyyz_1, g_yzzzz_xxyz_1, g_yzzzz_xxyzz_1, g_yzzzz_xxzz_1, g_yzzzz_xxzzz_1, g_yzzzz_xyyyz_1, g_yzzzz_xyyz_1, g_yzzzz_xyyzz_1, g_yzzzz_xyzz_1, g_yzzzz_xyzzz_1, g_yzzzz_xzzz_1, g_yzzzz_xzzzz_1, g_yzzzz_yyyyz_1, g_yzzzz_yyyz_1, g_yzzzz_yyyzz_1, g_yzzzz_yyzz_1, g_yzzzz_yyzzz_1, g_yzzzz_yzzz_1, g_yzzzz_yzzzz_1, g_yzzzz_zzzz_1, g_yzzzz_zzzzz_1, g_zzzz_xxxxx_0, g_zzzz_xxxxx_1, g_zzzz_xxxxz_0, g_zzzz_xxxxz_1, g_zzzz_xxxyz_0, g_zzzz_xxxyz_1, g_zzzz_xxxzz_0, g_zzzz_xxxzz_1, g_zzzz_xxyyz_0, g_zzzz_xxyyz_1, g_zzzz_xxyzz_0, g_zzzz_xxyzz_1, g_zzzz_xxzzz_0, g_zzzz_xxzzz_1, g_zzzz_xyyyz_0, g_zzzz_xyyyz_1, g_zzzz_xyyzz_0, g_zzzz_xyyzz_1, g_zzzz_xyzzz_0, g_zzzz_xyzzz_1, g_zzzz_xzzzz_0, g_zzzz_xzzzz_1, g_zzzz_yyyyz_0, g_zzzz_yyyyz_1, g_zzzz_yyyzz_0, g_zzzz_yyyzz_1, g_zzzz_yyzzz_0, g_zzzz_yyzzz_1, g_zzzz_yzzzz_0, g_zzzz_yzzzz_1, g_zzzz_zzzzz_0, g_zzzz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzzz_xxxxx_0[i] = g_zzzz_xxxxx_0[i] * fbe_0 - g_zzzz_xxxxx_1[i] * fz_be_0 + g_yzzzz_xxxxx_1[i] * pa_y[i];

        g_yyzzzz_xxxxy_0[i] = 3.0 * g_yyzz_xxxxy_0[i] * fbe_0 - 3.0 * g_yyzz_xxxxy_1[i] * fz_be_0 + g_yyzzz_xxxxy_1[i] * pa_z[i];

        g_yyzzzz_xxxxz_0[i] = g_zzzz_xxxxz_0[i] * fbe_0 - g_zzzz_xxxxz_1[i] * fz_be_0 + g_yzzzz_xxxxz_1[i] * pa_y[i];

        g_yyzzzz_xxxyy_0[i] = 3.0 * g_yyzz_xxxyy_0[i] * fbe_0 - 3.0 * g_yyzz_xxxyy_1[i] * fz_be_0 + g_yyzzz_xxxyy_1[i] * pa_z[i];

        g_yyzzzz_xxxyz_0[i] = g_zzzz_xxxyz_0[i] * fbe_0 - g_zzzz_xxxyz_1[i] * fz_be_0 + g_yzzzz_xxxz_1[i] * fe_0 + g_yzzzz_xxxyz_1[i] * pa_y[i];

        g_yyzzzz_xxxzz_0[i] = g_zzzz_xxxzz_0[i] * fbe_0 - g_zzzz_xxxzz_1[i] * fz_be_0 + g_yzzzz_xxxzz_1[i] * pa_y[i];

        g_yyzzzz_xxyyy_0[i] = 3.0 * g_yyzz_xxyyy_0[i] * fbe_0 - 3.0 * g_yyzz_xxyyy_1[i] * fz_be_0 + g_yyzzz_xxyyy_1[i] * pa_z[i];

        g_yyzzzz_xxyyz_0[i] = g_zzzz_xxyyz_0[i] * fbe_0 - g_zzzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_yzzzz_xxyz_1[i] * fe_0 + g_yzzzz_xxyyz_1[i] * pa_y[i];

        g_yyzzzz_xxyzz_0[i] = g_zzzz_xxyzz_0[i] * fbe_0 - g_zzzz_xxyzz_1[i] * fz_be_0 + g_yzzzz_xxzz_1[i] * fe_0 + g_yzzzz_xxyzz_1[i] * pa_y[i];

        g_yyzzzz_xxzzz_0[i] = g_zzzz_xxzzz_0[i] * fbe_0 - g_zzzz_xxzzz_1[i] * fz_be_0 + g_yzzzz_xxzzz_1[i] * pa_y[i];

        g_yyzzzz_xyyyy_0[i] = 3.0 * g_yyzz_xyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_xyyyy_1[i] * fz_be_0 + g_yyzzz_xyyyy_1[i] * pa_z[i];

        g_yyzzzz_xyyyz_0[i] = g_zzzz_xyyyz_0[i] * fbe_0 - g_zzzz_xyyyz_1[i] * fz_be_0 + 3.0 * g_yzzzz_xyyz_1[i] * fe_0 + g_yzzzz_xyyyz_1[i] * pa_y[i];

        g_yyzzzz_xyyzz_0[i] = g_zzzz_xyyzz_0[i] * fbe_0 - g_zzzz_xyyzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_xyzz_1[i] * fe_0 + g_yzzzz_xyyzz_1[i] * pa_y[i];

        g_yyzzzz_xyzzz_0[i] = g_zzzz_xyzzz_0[i] * fbe_0 - g_zzzz_xyzzz_1[i] * fz_be_0 + g_yzzzz_xzzz_1[i] * fe_0 + g_yzzzz_xyzzz_1[i] * pa_y[i];

        g_yyzzzz_xzzzz_0[i] = g_zzzz_xzzzz_0[i] * fbe_0 - g_zzzz_xzzzz_1[i] * fz_be_0 + g_yzzzz_xzzzz_1[i] * pa_y[i];

        g_yyzzzz_yyyyy_0[i] = 3.0 * g_yyzz_yyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_yyyyy_1[i] * fz_be_0 + g_yyzzz_yyyyy_1[i] * pa_z[i];

        g_yyzzzz_yyyyz_0[i] = g_zzzz_yyyyz_0[i] * fbe_0 - g_zzzz_yyyyz_1[i] * fz_be_0 + 4.0 * g_yzzzz_yyyz_1[i] * fe_0 + g_yzzzz_yyyyz_1[i] * pa_y[i];

        g_yyzzzz_yyyzz_0[i] = g_zzzz_yyyzz_0[i] * fbe_0 - g_zzzz_yyyzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_yyzz_1[i] * fe_0 + g_yzzzz_yyyzz_1[i] * pa_y[i];

        g_yyzzzz_yyzzz_0[i] = g_zzzz_yyzzz_0[i] * fbe_0 - g_zzzz_yyzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_yzzz_1[i] * fe_0 + g_yzzzz_yyzzz_1[i] * pa_y[i];

        g_yyzzzz_yzzzz_0[i] = g_zzzz_yzzzz_0[i] * fbe_0 - g_zzzz_yzzzz_1[i] * fz_be_0 + g_yzzzz_zzzz_1[i] * fe_0 + g_yzzzz_yzzzz_1[i] * pa_y[i];

        g_yyzzzz_zzzzz_0[i] = g_zzzz_zzzzz_0[i] * fbe_0 - g_zzzz_zzzzz_1[i] * fz_be_0 + g_yzzzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 546-567 components of targeted buffer : IH

    auto g_yzzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 546);

    auto g_yzzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 547);

    auto g_yzzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 548);

    auto g_yzzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 549);

    auto g_yzzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 550);

    auto g_yzzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 551);

    auto g_yzzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 552);

    auto g_yzzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 553);

    auto g_yzzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 554);

    auto g_yzzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 555);

    auto g_yzzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 556);

    auto g_yzzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 557);

    auto g_yzzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 558);

    auto g_yzzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 559);

    auto g_yzzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 560);

    auto g_yzzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 561);

    auto g_yzzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 562);

    auto g_yzzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 563);

    auto g_yzzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 564);

    auto g_yzzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 565);

    auto g_yzzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 566);

    #pragma omp simd aligned(g_yzzzzz_xxxxx_0, g_yzzzzz_xxxxy_0, g_yzzzzz_xxxxz_0, g_yzzzzz_xxxyy_0, g_yzzzzz_xxxyz_0, g_yzzzzz_xxxzz_0, g_yzzzzz_xxyyy_0, g_yzzzzz_xxyyz_0, g_yzzzzz_xxyzz_0, g_yzzzzz_xxzzz_0, g_yzzzzz_xyyyy_0, g_yzzzzz_xyyyz_0, g_yzzzzz_xyyzz_0, g_yzzzzz_xyzzz_0, g_yzzzzz_xzzzz_0, g_yzzzzz_yyyyy_0, g_yzzzzz_yyyyz_0, g_yzzzzz_yyyzz_0, g_yzzzzz_yyzzz_0, g_yzzzzz_yzzzz_0, g_yzzzzz_zzzzz_0, g_zzzzz_xxxx_1, g_zzzzz_xxxxx_1, g_zzzzz_xxxxy_1, g_zzzzz_xxxxz_1, g_zzzzz_xxxy_1, g_zzzzz_xxxyy_1, g_zzzzz_xxxyz_1, g_zzzzz_xxxz_1, g_zzzzz_xxxzz_1, g_zzzzz_xxyy_1, g_zzzzz_xxyyy_1, g_zzzzz_xxyyz_1, g_zzzzz_xxyz_1, g_zzzzz_xxyzz_1, g_zzzzz_xxzz_1, g_zzzzz_xxzzz_1, g_zzzzz_xyyy_1, g_zzzzz_xyyyy_1, g_zzzzz_xyyyz_1, g_zzzzz_xyyz_1, g_zzzzz_xyyzz_1, g_zzzzz_xyzz_1, g_zzzzz_xyzzz_1, g_zzzzz_xzzz_1, g_zzzzz_xzzzz_1, g_zzzzz_yyyy_1, g_zzzzz_yyyyy_1, g_zzzzz_yyyyz_1, g_zzzzz_yyyz_1, g_zzzzz_yyyzz_1, g_zzzzz_yyzz_1, g_zzzzz_yyzzz_1, g_zzzzz_yzzz_1, g_zzzzz_yzzzz_1, g_zzzzz_zzzz_1, g_zzzzz_zzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzzz_xxxxx_0[i] = g_zzzzz_xxxxx_1[i] * pa_y[i];

        g_yzzzzz_xxxxy_0[i] = g_zzzzz_xxxx_1[i] * fe_0 + g_zzzzz_xxxxy_1[i] * pa_y[i];

        g_yzzzzz_xxxxz_0[i] = g_zzzzz_xxxxz_1[i] * pa_y[i];

        g_yzzzzz_xxxyy_0[i] = 2.0 * g_zzzzz_xxxy_1[i] * fe_0 + g_zzzzz_xxxyy_1[i] * pa_y[i];

        g_yzzzzz_xxxyz_0[i] = g_zzzzz_xxxz_1[i] * fe_0 + g_zzzzz_xxxyz_1[i] * pa_y[i];

        g_yzzzzz_xxxzz_0[i] = g_zzzzz_xxxzz_1[i] * pa_y[i];

        g_yzzzzz_xxyyy_0[i] = 3.0 * g_zzzzz_xxyy_1[i] * fe_0 + g_zzzzz_xxyyy_1[i] * pa_y[i];

        g_yzzzzz_xxyyz_0[i] = 2.0 * g_zzzzz_xxyz_1[i] * fe_0 + g_zzzzz_xxyyz_1[i] * pa_y[i];

        g_yzzzzz_xxyzz_0[i] = g_zzzzz_xxzz_1[i] * fe_0 + g_zzzzz_xxyzz_1[i] * pa_y[i];

        g_yzzzzz_xxzzz_0[i] = g_zzzzz_xxzzz_1[i] * pa_y[i];

        g_yzzzzz_xyyyy_0[i] = 4.0 * g_zzzzz_xyyy_1[i] * fe_0 + g_zzzzz_xyyyy_1[i] * pa_y[i];

        g_yzzzzz_xyyyz_0[i] = 3.0 * g_zzzzz_xyyz_1[i] * fe_0 + g_zzzzz_xyyyz_1[i] * pa_y[i];

        g_yzzzzz_xyyzz_0[i] = 2.0 * g_zzzzz_xyzz_1[i] * fe_0 + g_zzzzz_xyyzz_1[i] * pa_y[i];

        g_yzzzzz_xyzzz_0[i] = g_zzzzz_xzzz_1[i] * fe_0 + g_zzzzz_xyzzz_1[i] * pa_y[i];

        g_yzzzzz_xzzzz_0[i] = g_zzzzz_xzzzz_1[i] * pa_y[i];

        g_yzzzzz_yyyyy_0[i] = 5.0 * g_zzzzz_yyyy_1[i] * fe_0 + g_zzzzz_yyyyy_1[i] * pa_y[i];

        g_yzzzzz_yyyyz_0[i] = 4.0 * g_zzzzz_yyyz_1[i] * fe_0 + g_zzzzz_yyyyz_1[i] * pa_y[i];

        g_yzzzzz_yyyzz_0[i] = 3.0 * g_zzzzz_yyzz_1[i] * fe_0 + g_zzzzz_yyyzz_1[i] * pa_y[i];

        g_yzzzzz_yyzzz_0[i] = 2.0 * g_zzzzz_yzzz_1[i] * fe_0 + g_zzzzz_yyzzz_1[i] * pa_y[i];

        g_yzzzzz_yzzzz_0[i] = g_zzzzz_zzzz_1[i] * fe_0 + g_zzzzz_yzzzz_1[i] * pa_y[i];

        g_yzzzzz_zzzzz_0[i] = g_zzzzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 567-588 components of targeted buffer : IH

    auto g_zzzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_ih + 567);

    auto g_zzzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_ih + 568);

    auto g_zzzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_ih + 569);

    auto g_zzzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_ih + 570);

    auto g_zzzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_ih + 571);

    auto g_zzzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_ih + 572);

    auto g_zzzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_ih + 573);

    auto g_zzzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_ih + 574);

    auto g_zzzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_ih + 575);

    auto g_zzzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_ih + 576);

    auto g_zzzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_ih + 577);

    auto g_zzzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_ih + 578);

    auto g_zzzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_ih + 579);

    auto g_zzzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_ih + 580);

    auto g_zzzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_ih + 581);

    auto g_zzzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_ih + 582);

    auto g_zzzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_ih + 583);

    auto g_zzzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_ih + 584);

    auto g_zzzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_ih + 585);

    auto g_zzzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_ih + 586);

    auto g_zzzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_ih + 587);

    #pragma omp simd aligned(g_zzzz_xxxxx_0, g_zzzz_xxxxx_1, g_zzzz_xxxxy_0, g_zzzz_xxxxy_1, g_zzzz_xxxxz_0, g_zzzz_xxxxz_1, g_zzzz_xxxyy_0, g_zzzz_xxxyy_1, g_zzzz_xxxyz_0, g_zzzz_xxxyz_1, g_zzzz_xxxzz_0, g_zzzz_xxxzz_1, g_zzzz_xxyyy_0, g_zzzz_xxyyy_1, g_zzzz_xxyyz_0, g_zzzz_xxyyz_1, g_zzzz_xxyzz_0, g_zzzz_xxyzz_1, g_zzzz_xxzzz_0, g_zzzz_xxzzz_1, g_zzzz_xyyyy_0, g_zzzz_xyyyy_1, g_zzzz_xyyyz_0, g_zzzz_xyyyz_1, g_zzzz_xyyzz_0, g_zzzz_xyyzz_1, g_zzzz_xyzzz_0, g_zzzz_xyzzz_1, g_zzzz_xzzzz_0, g_zzzz_xzzzz_1, g_zzzz_yyyyy_0, g_zzzz_yyyyy_1, g_zzzz_yyyyz_0, g_zzzz_yyyyz_1, g_zzzz_yyyzz_0, g_zzzz_yyyzz_1, g_zzzz_yyzzz_0, g_zzzz_yyzzz_1, g_zzzz_yzzzz_0, g_zzzz_yzzzz_1, g_zzzz_zzzzz_0, g_zzzz_zzzzz_1, g_zzzzz_xxxx_1, g_zzzzz_xxxxx_1, g_zzzzz_xxxxy_1, g_zzzzz_xxxxz_1, g_zzzzz_xxxy_1, g_zzzzz_xxxyy_1, g_zzzzz_xxxyz_1, g_zzzzz_xxxz_1, g_zzzzz_xxxzz_1, g_zzzzz_xxyy_1, g_zzzzz_xxyyy_1, g_zzzzz_xxyyz_1, g_zzzzz_xxyz_1, g_zzzzz_xxyzz_1, g_zzzzz_xxzz_1, g_zzzzz_xxzzz_1, g_zzzzz_xyyy_1, g_zzzzz_xyyyy_1, g_zzzzz_xyyyz_1, g_zzzzz_xyyz_1, g_zzzzz_xyyzz_1, g_zzzzz_xyzz_1, g_zzzzz_xyzzz_1, g_zzzzz_xzzz_1, g_zzzzz_xzzzz_1, g_zzzzz_yyyy_1, g_zzzzz_yyyyy_1, g_zzzzz_yyyyz_1, g_zzzzz_yyyz_1, g_zzzzz_yyyzz_1, g_zzzzz_yyzz_1, g_zzzzz_yyzzz_1, g_zzzzz_yzzz_1, g_zzzzz_yzzzz_1, g_zzzzz_zzzz_1, g_zzzzz_zzzzz_1, g_zzzzzz_xxxxx_0, g_zzzzzz_xxxxy_0, g_zzzzzz_xxxxz_0, g_zzzzzz_xxxyy_0, g_zzzzzz_xxxyz_0, g_zzzzzz_xxxzz_0, g_zzzzzz_xxyyy_0, g_zzzzzz_xxyyz_0, g_zzzzzz_xxyzz_0, g_zzzzzz_xxzzz_0, g_zzzzzz_xyyyy_0, g_zzzzzz_xyyyz_0, g_zzzzzz_xyyzz_0, g_zzzzzz_xyzzz_0, g_zzzzzz_xzzzz_0, g_zzzzzz_yyyyy_0, g_zzzzzz_yyyyz_0, g_zzzzzz_yyyzz_0, g_zzzzzz_yyzzz_0, g_zzzzzz_yzzzz_0, g_zzzzzz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzzz_xxxxx_0[i] = 5.0 * g_zzzz_xxxxx_0[i] * fbe_0 - 5.0 * g_zzzz_xxxxx_1[i] * fz_be_0 + g_zzzzz_xxxxx_1[i] * pa_z[i];

        g_zzzzzz_xxxxy_0[i] = 5.0 * g_zzzz_xxxxy_0[i] * fbe_0 - 5.0 * g_zzzz_xxxxy_1[i] * fz_be_0 + g_zzzzz_xxxxy_1[i] * pa_z[i];

        g_zzzzzz_xxxxz_0[i] = 5.0 * g_zzzz_xxxxz_0[i] * fbe_0 - 5.0 * g_zzzz_xxxxz_1[i] * fz_be_0 + g_zzzzz_xxxx_1[i] * fe_0 + g_zzzzz_xxxxz_1[i] * pa_z[i];

        g_zzzzzz_xxxyy_0[i] = 5.0 * g_zzzz_xxxyy_0[i] * fbe_0 - 5.0 * g_zzzz_xxxyy_1[i] * fz_be_0 + g_zzzzz_xxxyy_1[i] * pa_z[i];

        g_zzzzzz_xxxyz_0[i] = 5.0 * g_zzzz_xxxyz_0[i] * fbe_0 - 5.0 * g_zzzz_xxxyz_1[i] * fz_be_0 + g_zzzzz_xxxy_1[i] * fe_0 + g_zzzzz_xxxyz_1[i] * pa_z[i];

        g_zzzzzz_xxxzz_0[i] = 5.0 * g_zzzz_xxxzz_0[i] * fbe_0 - 5.0 * g_zzzz_xxxzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_xxxz_1[i] * fe_0 + g_zzzzz_xxxzz_1[i] * pa_z[i];

        g_zzzzzz_xxyyy_0[i] = 5.0 * g_zzzz_xxyyy_0[i] * fbe_0 - 5.0 * g_zzzz_xxyyy_1[i] * fz_be_0 + g_zzzzz_xxyyy_1[i] * pa_z[i];

        g_zzzzzz_xxyyz_0[i] = 5.0 * g_zzzz_xxyyz_0[i] * fbe_0 - 5.0 * g_zzzz_xxyyz_1[i] * fz_be_0 + g_zzzzz_xxyy_1[i] * fe_0 + g_zzzzz_xxyyz_1[i] * pa_z[i];

        g_zzzzzz_xxyzz_0[i] = 5.0 * g_zzzz_xxyzz_0[i] * fbe_0 - 5.0 * g_zzzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_xxyz_1[i] * fe_0 + g_zzzzz_xxyzz_1[i] * pa_z[i];

        g_zzzzzz_xxzzz_0[i] = 5.0 * g_zzzz_xxzzz_0[i] * fbe_0 - 5.0 * g_zzzz_xxzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_xxzz_1[i] * fe_0 + g_zzzzz_xxzzz_1[i] * pa_z[i];

        g_zzzzzz_xyyyy_0[i] = 5.0 * g_zzzz_xyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_xyyyy_1[i] * fz_be_0 + g_zzzzz_xyyyy_1[i] * pa_z[i];

        g_zzzzzz_xyyyz_0[i] = 5.0 * g_zzzz_xyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_xyyyz_1[i] * fz_be_0 + g_zzzzz_xyyy_1[i] * fe_0 + g_zzzzz_xyyyz_1[i] * pa_z[i];

        g_zzzzzz_xyyzz_0[i] = 5.0 * g_zzzz_xyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_xyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_xyyz_1[i] * fe_0 + g_zzzzz_xyyzz_1[i] * pa_z[i];

        g_zzzzzz_xyzzz_0[i] = 5.0 * g_zzzz_xyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_xyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_xyzz_1[i] * fe_0 + g_zzzzz_xyzzz_1[i] * pa_z[i];

        g_zzzzzz_xzzzz_0[i] = 5.0 * g_zzzz_xzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_xzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_xzzz_1[i] * fe_0 + g_zzzzz_xzzzz_1[i] * pa_z[i];

        g_zzzzzz_yyyyy_0[i] = 5.0 * g_zzzz_yyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_yyyyy_1[i] * fz_be_0 + g_zzzzz_yyyyy_1[i] * pa_z[i];

        g_zzzzzz_yyyyz_0[i] = 5.0 * g_zzzz_yyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_yyyyz_1[i] * fz_be_0 + g_zzzzz_yyyy_1[i] * fe_0 + g_zzzzz_yyyyz_1[i] * pa_z[i];

        g_zzzzzz_yyyzz_0[i] = 5.0 * g_zzzz_yyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_yyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_yyyz_1[i] * fe_0 + g_zzzzz_yyyzz_1[i] * pa_z[i];

        g_zzzzzz_yyzzz_0[i] = 5.0 * g_zzzz_yyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_yyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_yyzz_1[i] * fe_0 + g_zzzzz_yyzzz_1[i] * pa_z[i];

        g_zzzzzz_yzzzz_0[i] = 5.0 * g_zzzz_yzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_yzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_yzzz_1[i] * fe_0 + g_zzzzz_yzzzz_1[i] * pa_z[i];

        g_zzzzzz_zzzzz_0[i] = 5.0 * g_zzzz_zzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_zzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_zzzz_1[i] * fe_0 + g_zzzzz_zzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

