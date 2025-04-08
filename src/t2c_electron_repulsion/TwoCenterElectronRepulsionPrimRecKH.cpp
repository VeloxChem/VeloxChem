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

#include "TwoCenterElectronRepulsionPrimRecKH.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_kh(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_kh,
                                const size_t idx_eri_0_hh,
                                const size_t idx_eri_1_hh,
                                const size_t idx_eri_1_ig,
                                const size_t idx_eri_1_ih,
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

    // Set up components of auxiliary buffer : HH

    auto g_xxxxx_xxxxx_0 = pbuffer.data(idx_eri_0_hh);

    auto g_xxxxx_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 1);

    auto g_xxxxx_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 2);

    auto g_xxxxx_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 3);

    auto g_xxxxx_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 4);

    auto g_xxxxx_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 5);

    auto g_xxxxx_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 6);

    auto g_xxxxx_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 7);

    auto g_xxxxx_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 8);

    auto g_xxxxx_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 9);

    auto g_xxxxx_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 10);

    auto g_xxxxx_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 11);

    auto g_xxxxx_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 12);

    auto g_xxxxx_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 13);

    auto g_xxxxx_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 14);

    auto g_xxxxx_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 15);

    auto g_xxxxx_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 16);

    auto g_xxxxx_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 17);

    auto g_xxxxx_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 18);

    auto g_xxxxx_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 19);

    auto g_xxxxx_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 20);

    auto g_xxxxy_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 21);

    auto g_xxxxy_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 23);

    auto g_xxxxy_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 26);

    auto g_xxxxy_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 30);

    auto g_xxxxy_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 35);

    auto g_xxxxz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 42);

    auto g_xxxxz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 43);

    auto g_xxxxz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 45);

    auto g_xxxxz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 48);

    auto g_xxxxz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 52);

    auto g_xxxyy_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 63);

    auto g_xxxyy_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 64);

    auto g_xxxyy_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 65);

    auto g_xxxyy_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 66);

    auto g_xxxyy_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 67);

    auto g_xxxyy_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 68);

    auto g_xxxyy_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 69);

    auto g_xxxyy_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 70);

    auto g_xxxyy_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 71);

    auto g_xxxyy_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 72);

    auto g_xxxyy_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 73);

    auto g_xxxyy_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 74);

    auto g_xxxyy_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 75);

    auto g_xxxyy_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 76);

    auto g_xxxyy_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 77);

    auto g_xxxyy_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 78);

    auto g_xxxyy_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 79);

    auto g_xxxyy_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 80);

    auto g_xxxyy_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 81);

    auto g_xxxyy_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 82);

    auto g_xxxyy_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 83);

    auto g_xxxzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 105);

    auto g_xxxzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 106);

    auto g_xxxzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 107);

    auto g_xxxzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 108);

    auto g_xxxzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 109);

    auto g_xxxzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 110);

    auto g_xxxzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 111);

    auto g_xxxzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 112);

    auto g_xxxzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 113);

    auto g_xxxzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 114);

    auto g_xxxzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 115);

    auto g_xxxzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 116);

    auto g_xxxzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 117);

    auto g_xxxzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 118);

    auto g_xxxzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 119);

    auto g_xxxzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 120);

    auto g_xxxzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 121);

    auto g_xxxzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 122);

    auto g_xxxzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 123);

    auto g_xxxzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 124);

    auto g_xxxzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 125);

    auto g_xxyyy_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 126);

    auto g_xxyyy_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 127);

    auto g_xxyyy_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 128);

    auto g_xxyyy_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 129);

    auto g_xxyyy_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 130);

    auto g_xxyyy_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 131);

    auto g_xxyyy_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 132);

    auto g_xxyyy_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 133);

    auto g_xxyyy_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 134);

    auto g_xxyyy_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 135);

    auto g_xxyyy_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 136);

    auto g_xxyyy_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 137);

    auto g_xxyyy_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 138);

    auto g_xxyyy_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 139);

    auto g_xxyyy_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 140);

    auto g_xxyyy_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 141);

    auto g_xxyyy_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 142);

    auto g_xxyyy_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 143);

    auto g_xxyyy_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 144);

    auto g_xxyyy_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 145);

    auto g_xxyyy_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 146);

    auto g_xxyyz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 148);

    auto g_xxyyz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 150);

    auto g_xxyyz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 153);

    auto g_xxyyz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 157);

    auto g_xxyzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 168);

    auto g_xxyzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 170);

    auto g_xxyzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 173);

    auto g_xxyzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 177);

    auto g_xxyzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 182);

    auto g_xxzzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 189);

    auto g_xxzzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 190);

    auto g_xxzzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 191);

    auto g_xxzzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 192);

    auto g_xxzzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 193);

    auto g_xxzzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 194);

    auto g_xxzzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 195);

    auto g_xxzzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 196);

    auto g_xxzzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 197);

    auto g_xxzzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 198);

    auto g_xxzzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 199);

    auto g_xxzzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 200);

    auto g_xxzzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 201);

    auto g_xxzzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 202);

    auto g_xxzzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 203);

    auto g_xxzzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 204);

    auto g_xxzzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 205);

    auto g_xxzzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 206);

    auto g_xxzzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 207);

    auto g_xxzzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 208);

    auto g_xxzzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 209);

    auto g_xyyyy_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 211);

    auto g_xyyyy_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 213);

    auto g_xyyyy_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 214);

    auto g_xyyyy_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 216);

    auto g_xyyyy_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 217);

    auto g_xyyyy_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 218);

    auto g_xyyyy_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 220);

    auto g_xyyyy_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 221);

    auto g_xyyyy_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 222);

    auto g_xyyyy_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 223);

    auto g_xyyyy_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 225);

    auto g_xyyyy_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 226);

    auto g_xyyyy_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 227);

    auto g_xyyyy_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 228);

    auto g_xyyyy_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 229);

    auto g_xyyyy_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 230);

    auto g_xyyzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 256);

    auto g_xyyzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 259);

    auto g_xyyzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 260);

    auto g_xyyzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 263);

    auto g_xyyzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 264);

    auto g_xyyzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 265);

    auto g_xyyzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 267);

    auto g_xyyzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 268);

    auto g_xyyzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 269);

    auto g_xyyzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 270);

    auto g_xyyzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 271);

    auto g_xyyzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 272);

    auto g_xzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 296);

    auto g_xzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 298);

    auto g_xzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 299);

    auto g_xzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 301);

    auto g_xzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 302);

    auto g_xzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 303);

    auto g_xzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 305);

    auto g_xzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 306);

    auto g_xzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 307);

    auto g_xzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 308);

    auto g_xzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 309);

    auto g_xzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 310);

    auto g_xzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 311);

    auto g_xzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 312);

    auto g_xzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 313);

    auto g_xzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 314);

    auto g_yyyyy_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 315);

    auto g_yyyyy_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 316);

    auto g_yyyyy_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 317);

    auto g_yyyyy_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 318);

    auto g_yyyyy_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 319);

    auto g_yyyyy_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 320);

    auto g_yyyyy_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 321);

    auto g_yyyyy_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 322);

    auto g_yyyyy_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 323);

    auto g_yyyyy_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 324);

    auto g_yyyyy_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 325);

    auto g_yyyyy_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 326);

    auto g_yyyyy_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 327);

    auto g_yyyyy_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 328);

    auto g_yyyyy_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 329);

    auto g_yyyyy_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 330);

    auto g_yyyyy_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 331);

    auto g_yyyyy_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 332);

    auto g_yyyyy_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 333);

    auto g_yyyyy_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 334);

    auto g_yyyyy_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 335);

    auto g_yyyyz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 337);

    auto g_yyyyz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 339);

    auto g_yyyyz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 342);

    auto g_yyyyz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 346);

    auto g_yyyyz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 351);

    auto g_yyyzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 357);

    auto g_yyyzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 358);

    auto g_yyyzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 359);

    auto g_yyyzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 360);

    auto g_yyyzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 361);

    auto g_yyyzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 362);

    auto g_yyyzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 363);

    auto g_yyyzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 364);

    auto g_yyyzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 365);

    auto g_yyyzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 366);

    auto g_yyyzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 367);

    auto g_yyyzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 368);

    auto g_yyyzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 369);

    auto g_yyyzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 370);

    auto g_yyyzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 371);

    auto g_yyyzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 372);

    auto g_yyyzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 373);

    auto g_yyyzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 374);

    auto g_yyyzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 375);

    auto g_yyyzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 376);

    auto g_yyyzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 377);

    auto g_yyzzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 378);

    auto g_yyzzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 379);

    auto g_yyzzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 380);

    auto g_yyzzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 381);

    auto g_yyzzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 382);

    auto g_yyzzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 383);

    auto g_yyzzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 384);

    auto g_yyzzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 385);

    auto g_yyzzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 386);

    auto g_yyzzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 387);

    auto g_yyzzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 388);

    auto g_yyzzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 389);

    auto g_yyzzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 390);

    auto g_yyzzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 391);

    auto g_yyzzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 392);

    auto g_yyzzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 393);

    auto g_yyzzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 394);

    auto g_yyzzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 395);

    auto g_yyzzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 396);

    auto g_yyzzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 397);

    auto g_yyzzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 398);

    auto g_yzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 399);

    auto g_yzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 401);

    auto g_yzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 403);

    auto g_yzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 404);

    auto g_yzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 406);

    auto g_yzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 407);

    auto g_yzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 408);

    auto g_yzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 410);

    auto g_yzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 411);

    auto g_yzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 412);

    auto g_yzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 413);

    auto g_yzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 415);

    auto g_yzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 416);

    auto g_yzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 417);

    auto g_yzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 418);

    auto g_yzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 419);

    auto g_zzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_hh + 420);

    auto g_zzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_hh + 421);

    auto g_zzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_hh + 422);

    auto g_zzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_hh + 423);

    auto g_zzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_hh + 424);

    auto g_zzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_hh + 425);

    auto g_zzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_hh + 426);

    auto g_zzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_hh + 427);

    auto g_zzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_hh + 428);

    auto g_zzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_hh + 429);

    auto g_zzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_hh + 430);

    auto g_zzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_hh + 431);

    auto g_zzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_hh + 432);

    auto g_zzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_hh + 433);

    auto g_zzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_hh + 434);

    auto g_zzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_hh + 435);

    auto g_zzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_hh + 436);

    auto g_zzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_hh + 437);

    auto g_zzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_hh + 438);

    auto g_zzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_hh + 439);

    auto g_zzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_hh + 440);

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

    auto g_xxxxy_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 23);

    auto g_xxxxy_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 26);

    auto g_xxxxy_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 30);

    auto g_xxxxy_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 35);

    auto g_xxxxz_xxxxx_1 = pbuffer.data(idx_eri_1_hh + 42);

    auto g_xxxxz_xxxxy_1 = pbuffer.data(idx_eri_1_hh + 43);

    auto g_xxxxz_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 45);

    auto g_xxxxz_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 48);

    auto g_xxxxz_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 52);

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

    auto g_yyyyz_xxxyy_1 = pbuffer.data(idx_eri_1_hh + 339);

    auto g_yyyyz_xxyyy_1 = pbuffer.data(idx_eri_1_hh + 342);

    auto g_yyyyz_xyyyy_1 = pbuffer.data(idx_eri_1_hh + 346);

    auto g_yyyyz_yyyyy_1 = pbuffer.data(idx_eri_1_hh + 351);

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

    auto g_yzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 401);

    auto g_yzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 403);

    auto g_yzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 404);

    auto g_yzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 406);

    auto g_yzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 407);

    auto g_yzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 408);

    auto g_yzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 410);

    auto g_yzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 411);

    auto g_yzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 412);

    auto g_yzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 413);

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

    auto g_xxxxxz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 32);

    auto g_xxxxxz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 34);

    auto g_xxxxxz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 35);

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

    auto g_xxyyzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 184);

    auto g_xxyyzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 187);

    auto g_xxyyzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 188);

    auto g_xxyyzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 191);

    auto g_xxyyzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 192);

    auto g_xxyyzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 193);

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

    auto g_xyyyzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 259);

    auto g_xyyyzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 262);

    auto g_xyyyzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 263);

    auto g_xyyyzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 266);

    auto g_xyyyzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 267);

    auto g_xyyyzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 268);

    auto g_xyyzzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 274);

    auto g_xyyzzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 277);

    auto g_xyyzzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 278);

    auto g_xyyzzz_yyyz_1 = pbuffer.data(idx_eri_1_ig + 281);

    auto g_xyyzzz_yyzz_1 = pbuffer.data(idx_eri_1_ig + 282);

    auto g_xyyzzz_yzzz_1 = pbuffer.data(idx_eri_1_ig + 283);

    auto g_xzzzzz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 302);

    auto g_xzzzzz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 304);

    auto g_xzzzzz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 305);

    auto g_xzzzzz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 307);

    auto g_xzzzzz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 308);

    auto g_xzzzzz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 309);

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

    auto g_yyyyyz_xxxz_1 = pbuffer.data(idx_eri_1_ig + 332);

    auto g_yyyyyz_xxyz_1 = pbuffer.data(idx_eri_1_ig + 334);

    auto g_yyyyyz_xxzz_1 = pbuffer.data(idx_eri_1_ig + 335);

    auto g_yyyyyz_xyyz_1 = pbuffer.data(idx_eri_1_ig + 337);

    auto g_yyyyyz_xyzz_1 = pbuffer.data(idx_eri_1_ig + 338);

    auto g_yyyyyz_xzzz_1 = pbuffer.data(idx_eri_1_ig + 339);

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

    // Set up components of auxiliary buffer : IH

    auto g_xxxxxx_xxxxx_1 = pbuffer.data(idx_eri_1_ih);

    auto g_xxxxxx_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 1);

    auto g_xxxxxx_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 2);

    auto g_xxxxxx_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 3);

    auto g_xxxxxx_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 4);

    auto g_xxxxxx_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 5);

    auto g_xxxxxx_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 6);

    auto g_xxxxxx_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 7);

    auto g_xxxxxx_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 8);

    auto g_xxxxxx_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 9);

    auto g_xxxxxx_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 10);

    auto g_xxxxxx_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 11);

    auto g_xxxxxx_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 12);

    auto g_xxxxxx_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 13);

    auto g_xxxxxx_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 14);

    auto g_xxxxxx_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 15);

    auto g_xxxxxx_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 16);

    auto g_xxxxxx_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 17);

    auto g_xxxxxx_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 18);

    auto g_xxxxxx_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 19);

    auto g_xxxxxx_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 20);

    auto g_xxxxxy_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 21);

    auto g_xxxxxy_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 22);

    auto g_xxxxxy_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 23);

    auto g_xxxxxy_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 24);

    auto g_xxxxxy_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 26);

    auto g_xxxxxy_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 27);

    auto g_xxxxxy_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 30);

    auto g_xxxxxy_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 31);

    auto g_xxxxxy_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 35);

    auto g_xxxxxy_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 36);

    auto g_xxxxxz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 42);

    auto g_xxxxxz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 43);

    auto g_xxxxxz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 44);

    auto g_xxxxxz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 45);

    auto g_xxxxxz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 46);

    auto g_xxxxxz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 47);

    auto g_xxxxxz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 48);

    auto g_xxxxxz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 49);

    auto g_xxxxxz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 50);

    auto g_xxxxxz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 51);

    auto g_xxxxxz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 52);

    auto g_xxxxxz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 53);

    auto g_xxxxxz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 54);

    auto g_xxxxxz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 55);

    auto g_xxxxxz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 56);

    auto g_xxxxxz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 58);

    auto g_xxxxxz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 59);

    auto g_xxxxxz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 60);

    auto g_xxxxxz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 61);

    auto g_xxxxxz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 62);

    auto g_xxxxyy_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 63);

    auto g_xxxxyy_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 64);

    auto g_xxxxyy_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 65);

    auto g_xxxxyy_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 66);

    auto g_xxxxyy_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 67);

    auto g_xxxxyy_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 68);

    auto g_xxxxyy_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 69);

    auto g_xxxxyy_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 70);

    auto g_xxxxyy_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 71);

    auto g_xxxxyy_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 72);

    auto g_xxxxyy_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 73);

    auto g_xxxxyy_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 74);

    auto g_xxxxyy_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 75);

    auto g_xxxxyy_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 76);

    auto g_xxxxyy_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 77);

    auto g_xxxxyy_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 78);

    auto g_xxxxyy_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 79);

    auto g_xxxxyy_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 80);

    auto g_xxxxyy_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 81);

    auto g_xxxxyy_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 82);

    auto g_xxxxyy_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 83);

    auto g_xxxxzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 105);

    auto g_xxxxzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 106);

    auto g_xxxxzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 107);

    auto g_xxxxzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 108);

    auto g_xxxxzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 109);

    auto g_xxxxzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 110);

    auto g_xxxxzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 111);

    auto g_xxxxzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 112);

    auto g_xxxxzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 113);

    auto g_xxxxzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 114);

    auto g_xxxxzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 115);

    auto g_xxxxzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 116);

    auto g_xxxxzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 117);

    auto g_xxxxzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 118);

    auto g_xxxxzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 119);

    auto g_xxxxzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 120);

    auto g_xxxxzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 121);

    auto g_xxxxzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 122);

    auto g_xxxxzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 123);

    auto g_xxxxzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 124);

    auto g_xxxxzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 125);

    auto g_xxxyyy_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 126);

    auto g_xxxyyy_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 127);

    auto g_xxxyyy_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 128);

    auto g_xxxyyy_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 129);

    auto g_xxxyyy_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 130);

    auto g_xxxyyy_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 131);

    auto g_xxxyyy_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 132);

    auto g_xxxyyy_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 133);

    auto g_xxxyyy_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 134);

    auto g_xxxyyy_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 135);

    auto g_xxxyyy_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 136);

    auto g_xxxyyy_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 137);

    auto g_xxxyyy_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 138);

    auto g_xxxyyy_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 139);

    auto g_xxxyyy_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 140);

    auto g_xxxyyy_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 141);

    auto g_xxxyyy_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 142);

    auto g_xxxyyy_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 143);

    auto g_xxxyyy_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 144);

    auto g_xxxyyy_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 145);

    auto g_xxxyyy_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 146);

    auto g_xxxyyz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 148);

    auto g_xxxyyz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 150);

    auto g_xxxyyz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 153);

    auto g_xxxyyz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 157);

    auto g_xxxyzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 168);

    auto g_xxxyzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 170);

    auto g_xxxyzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 173);

    auto g_xxxyzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 177);

    auto g_xxxyzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 182);

    auto g_xxxzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 189);

    auto g_xxxzzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 190);

    auto g_xxxzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 191);

    auto g_xxxzzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 192);

    auto g_xxxzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 193);

    auto g_xxxzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 194);

    auto g_xxxzzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 195);

    auto g_xxxzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 196);

    auto g_xxxzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 197);

    auto g_xxxzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 198);

    auto g_xxxzzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 199);

    auto g_xxxzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 200);

    auto g_xxxzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 201);

    auto g_xxxzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 202);

    auto g_xxxzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 203);

    auto g_xxxzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 204);

    auto g_xxxzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 205);

    auto g_xxxzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 206);

    auto g_xxxzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 207);

    auto g_xxxzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 208);

    auto g_xxxzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 209);

    auto g_xxyyyy_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 210);

    auto g_xxyyyy_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 211);

    auto g_xxyyyy_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 212);

    auto g_xxyyyy_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 213);

    auto g_xxyyyy_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 214);

    auto g_xxyyyy_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 215);

    auto g_xxyyyy_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 216);

    auto g_xxyyyy_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 217);

    auto g_xxyyyy_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 218);

    auto g_xxyyyy_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 219);

    auto g_xxyyyy_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 220);

    auto g_xxyyyy_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 221);

    auto g_xxyyyy_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 222);

    auto g_xxyyyy_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 223);

    auto g_xxyyyy_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 224);

    auto g_xxyyyy_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 225);

    auto g_xxyyyy_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 226);

    auto g_xxyyyy_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 227);

    auto g_xxyyyy_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 228);

    auto g_xxyyyy_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 229);

    auto g_xxyyyy_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 230);

    auto g_xxyyyz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 232);

    auto g_xxyyyz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 234);

    auto g_xxyyyz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 237);

    auto g_xxyyyz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 241);

    auto g_xxyyzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 252);

    auto g_xxyyzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 253);

    auto g_xxyyzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 254);

    auto g_xxyyzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 255);

    auto g_xxyyzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 256);

    auto g_xxyyzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 257);

    auto g_xxyyzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 258);

    auto g_xxyyzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 259);

    auto g_xxyyzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 260);

    auto g_xxyyzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 261);

    auto g_xxyyzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 262);

    auto g_xxyyzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 263);

    auto g_xxyyzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 264);

    auto g_xxyyzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 265);

    auto g_xxyyzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 266);

    auto g_xxyyzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 267);

    auto g_xxyyzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 268);

    auto g_xxyyzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 269);

    auto g_xxyyzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 270);

    auto g_xxyyzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 271);

    auto g_xxyyzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 272);

    auto g_xxyzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 273);

    auto g_xxyzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 275);

    auto g_xxyzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 278);

    auto g_xxyzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 282);

    auto g_xxyzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 287);

    auto g_xxzzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 294);

    auto g_xxzzzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 295);

    auto g_xxzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 296);

    auto g_xxzzzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 297);

    auto g_xxzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 298);

    auto g_xxzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 299);

    auto g_xxzzzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 300);

    auto g_xxzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 301);

    auto g_xxzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 302);

    auto g_xxzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 303);

    auto g_xxzzzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 304);

    auto g_xxzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 305);

    auto g_xxzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 306);

    auto g_xxzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 307);

    auto g_xxzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 308);

    auto g_xxzzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 309);

    auto g_xxzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 310);

    auto g_xxzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 311);

    auto g_xxzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 312);

    auto g_xxzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 313);

    auto g_xxzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 314);

    auto g_xyyyyy_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 315);

    auto g_xyyyyy_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 316);

    auto g_xyyyyy_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 318);

    auto g_xyyyyy_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 319);

    auto g_xyyyyy_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 321);

    auto g_xyyyyy_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 322);

    auto g_xyyyyy_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 323);

    auto g_xyyyyy_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 325);

    auto g_xyyyyy_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 326);

    auto g_xyyyyy_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 327);

    auto g_xyyyyy_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 328);

    auto g_xyyyyy_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 330);

    auto g_xyyyyy_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 331);

    auto g_xyyyyy_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 332);

    auto g_xyyyyy_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 333);

    auto g_xyyyyy_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 334);

    auto g_xyyyyy_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 335);

    auto g_xyyyzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 361);

    auto g_xyyyzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 364);

    auto g_xyyyzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 365);

    auto g_xyyyzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 368);

    auto g_xyyyzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 369);

    auto g_xyyyzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 370);

    auto g_xyyyzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 372);

    auto g_xyyyzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 373);

    auto g_xyyyzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 374);

    auto g_xyyyzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 375);

    auto g_xyyyzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 376);

    auto g_xyyyzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 377);

    auto g_xyyzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 382);

    auto g_xyyzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 385);

    auto g_xyyzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 386);

    auto g_xyyzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 389);

    auto g_xyyzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 390);

    auto g_xyyzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 391);

    auto g_xyyzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 393);

    auto g_xyyzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 394);

    auto g_xyyzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 395);

    auto g_xyyzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 396);

    auto g_xyyzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 397);

    auto g_xyyzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 398);

    auto g_xzzzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 420);

    auto g_xzzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 422);

    auto g_xzzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 424);

    auto g_xzzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 425);

    auto g_xzzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 427);

    auto g_xzzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 428);

    auto g_xzzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 429);

    auto g_xzzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 431);

    auto g_xzzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 432);

    auto g_xzzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 433);

    auto g_xzzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 434);

    auto g_xzzzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 435);

    auto g_xzzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 436);

    auto g_xzzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 437);

    auto g_xzzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 438);

    auto g_xzzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 439);

    auto g_xzzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 440);

    auto g_yyyyyy_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 441);

    auto g_yyyyyy_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 442);

    auto g_yyyyyy_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 443);

    auto g_yyyyyy_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 444);

    auto g_yyyyyy_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 445);

    auto g_yyyyyy_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 446);

    auto g_yyyyyy_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 447);

    auto g_yyyyyy_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 448);

    auto g_yyyyyy_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 449);

    auto g_yyyyyy_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 450);

    auto g_yyyyyy_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 451);

    auto g_yyyyyy_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 452);

    auto g_yyyyyy_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 453);

    auto g_yyyyyy_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 454);

    auto g_yyyyyy_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 455);

    auto g_yyyyyy_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 456);

    auto g_yyyyyy_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 457);

    auto g_yyyyyy_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 458);

    auto g_yyyyyy_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 459);

    auto g_yyyyyy_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 460);

    auto g_yyyyyy_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 461);

    auto g_yyyyyz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 463);

    auto g_yyyyyz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 464);

    auto g_yyyyyz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 465);

    auto g_yyyyyz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 466);

    auto g_yyyyyz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 467);

    auto g_yyyyyz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 468);

    auto g_yyyyyz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 469);

    auto g_yyyyyz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 470);

    auto g_yyyyyz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 471);

    auto g_yyyyyz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 472);

    auto g_yyyyyz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 473);

    auto g_yyyyyz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 474);

    auto g_yyyyyz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 475);

    auto g_yyyyyz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 476);

    auto g_yyyyyz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 477);

    auto g_yyyyyz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 478);

    auto g_yyyyyz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 479);

    auto g_yyyyyz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 480);

    auto g_yyyyyz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 481);

    auto g_yyyyyz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 482);

    auto g_yyyyzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 483);

    auto g_yyyyzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 484);

    auto g_yyyyzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 485);

    auto g_yyyyzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 486);

    auto g_yyyyzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 487);

    auto g_yyyyzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 488);

    auto g_yyyyzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 489);

    auto g_yyyyzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 490);

    auto g_yyyyzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 491);

    auto g_yyyyzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 492);

    auto g_yyyyzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 493);

    auto g_yyyyzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 494);

    auto g_yyyyzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 495);

    auto g_yyyyzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 496);

    auto g_yyyyzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 497);

    auto g_yyyyzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 498);

    auto g_yyyyzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 499);

    auto g_yyyyzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 500);

    auto g_yyyyzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 501);

    auto g_yyyyzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 502);

    auto g_yyyyzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 503);

    auto g_yyyzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 504);

    auto g_yyyzzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 505);

    auto g_yyyzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 506);

    auto g_yyyzzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 507);

    auto g_yyyzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 508);

    auto g_yyyzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 509);

    auto g_yyyzzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 510);

    auto g_yyyzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 511);

    auto g_yyyzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 512);

    auto g_yyyzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 513);

    auto g_yyyzzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 514);

    auto g_yyyzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 515);

    auto g_yyyzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 516);

    auto g_yyyzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 517);

    auto g_yyyzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 518);

    auto g_yyyzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 519);

    auto g_yyyzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 520);

    auto g_yyyzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 521);

    auto g_yyyzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 522);

    auto g_yyyzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 523);

    auto g_yyyzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 524);

    auto g_yyzzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 525);

    auto g_yyzzzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 526);

    auto g_yyzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 527);

    auto g_yyzzzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 528);

    auto g_yyzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 529);

    auto g_yyzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 530);

    auto g_yyzzzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 531);

    auto g_yyzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 532);

    auto g_yyzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 533);

    auto g_yyzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 534);

    auto g_yyzzzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 535);

    auto g_yyzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 536);

    auto g_yyzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 537);

    auto g_yyzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 538);

    auto g_yyzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 539);

    auto g_yyzzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 540);

    auto g_yyzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 541);

    auto g_yyzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 542);

    auto g_yyzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 543);

    auto g_yyzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 544);

    auto g_yyzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 545);

    auto g_yzzzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 546);

    auto g_yzzzzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 547);

    auto g_yzzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 548);

    auto g_yzzzzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 549);

    auto g_yzzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 550);

    auto g_yzzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 551);

    auto g_yzzzzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 552);

    auto g_yzzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 553);

    auto g_yzzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 554);

    auto g_yzzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 555);

    auto g_yzzzzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 556);

    auto g_yzzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 557);

    auto g_yzzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 558);

    auto g_yzzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 559);

    auto g_yzzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 560);

    auto g_yzzzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 561);

    auto g_yzzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 562);

    auto g_yzzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 563);

    auto g_yzzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 564);

    auto g_yzzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 565);

    auto g_yzzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 566);

    auto g_zzzzzz_xxxxx_1 = pbuffer.data(idx_eri_1_ih + 567);

    auto g_zzzzzz_xxxxy_1 = pbuffer.data(idx_eri_1_ih + 568);

    auto g_zzzzzz_xxxxz_1 = pbuffer.data(idx_eri_1_ih + 569);

    auto g_zzzzzz_xxxyy_1 = pbuffer.data(idx_eri_1_ih + 570);

    auto g_zzzzzz_xxxyz_1 = pbuffer.data(idx_eri_1_ih + 571);

    auto g_zzzzzz_xxxzz_1 = pbuffer.data(idx_eri_1_ih + 572);

    auto g_zzzzzz_xxyyy_1 = pbuffer.data(idx_eri_1_ih + 573);

    auto g_zzzzzz_xxyyz_1 = pbuffer.data(idx_eri_1_ih + 574);

    auto g_zzzzzz_xxyzz_1 = pbuffer.data(idx_eri_1_ih + 575);

    auto g_zzzzzz_xxzzz_1 = pbuffer.data(idx_eri_1_ih + 576);

    auto g_zzzzzz_xyyyy_1 = pbuffer.data(idx_eri_1_ih + 577);

    auto g_zzzzzz_xyyyz_1 = pbuffer.data(idx_eri_1_ih + 578);

    auto g_zzzzzz_xyyzz_1 = pbuffer.data(idx_eri_1_ih + 579);

    auto g_zzzzzz_xyzzz_1 = pbuffer.data(idx_eri_1_ih + 580);

    auto g_zzzzzz_xzzzz_1 = pbuffer.data(idx_eri_1_ih + 581);

    auto g_zzzzzz_yyyyy_1 = pbuffer.data(idx_eri_1_ih + 582);

    auto g_zzzzzz_yyyyz_1 = pbuffer.data(idx_eri_1_ih + 583);

    auto g_zzzzzz_yyyzz_1 = pbuffer.data(idx_eri_1_ih + 584);

    auto g_zzzzzz_yyzzz_1 = pbuffer.data(idx_eri_1_ih + 585);

    auto g_zzzzzz_yzzzz_1 = pbuffer.data(idx_eri_1_ih + 586);

    auto g_zzzzzz_zzzzz_1 = pbuffer.data(idx_eri_1_ih + 587);

    // Set up 0-21 components of targeted buffer : KH

    auto g_xxxxxxx_xxxxx_0 = pbuffer.data(idx_eri_0_kh);

    auto g_xxxxxxx_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 1);

    auto g_xxxxxxx_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 2);

    auto g_xxxxxxx_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 3);

    auto g_xxxxxxx_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 4);

    auto g_xxxxxxx_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 5);

    auto g_xxxxxxx_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 6);

    auto g_xxxxxxx_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 7);

    auto g_xxxxxxx_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 8);

    auto g_xxxxxxx_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 9);

    auto g_xxxxxxx_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 10);

    auto g_xxxxxxx_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 11);

    auto g_xxxxxxx_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 12);

    auto g_xxxxxxx_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 13);

    auto g_xxxxxxx_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 14);

    auto g_xxxxxxx_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 15);

    auto g_xxxxxxx_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 16);

    auto g_xxxxxxx_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 17);

    auto g_xxxxxxx_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 18);

    auto g_xxxxxxx_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 19);

    auto g_xxxxxxx_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 20);

    #pragma omp simd aligned(g_xxxxx_xxxxx_0, g_xxxxx_xxxxx_1, g_xxxxx_xxxxy_0, g_xxxxx_xxxxy_1, g_xxxxx_xxxxz_0, g_xxxxx_xxxxz_1, g_xxxxx_xxxyy_0, g_xxxxx_xxxyy_1, g_xxxxx_xxxyz_0, g_xxxxx_xxxyz_1, g_xxxxx_xxxzz_0, g_xxxxx_xxxzz_1, g_xxxxx_xxyyy_0, g_xxxxx_xxyyy_1, g_xxxxx_xxyyz_0, g_xxxxx_xxyyz_1, g_xxxxx_xxyzz_0, g_xxxxx_xxyzz_1, g_xxxxx_xxzzz_0, g_xxxxx_xxzzz_1, g_xxxxx_xyyyy_0, g_xxxxx_xyyyy_1, g_xxxxx_xyyyz_0, g_xxxxx_xyyyz_1, g_xxxxx_xyyzz_0, g_xxxxx_xyyzz_1, g_xxxxx_xyzzz_0, g_xxxxx_xyzzz_1, g_xxxxx_xzzzz_0, g_xxxxx_xzzzz_1, g_xxxxx_yyyyy_0, g_xxxxx_yyyyy_1, g_xxxxx_yyyyz_0, g_xxxxx_yyyyz_1, g_xxxxx_yyyzz_0, g_xxxxx_yyyzz_1, g_xxxxx_yyzzz_0, g_xxxxx_yyzzz_1, g_xxxxx_yzzzz_0, g_xxxxx_yzzzz_1, g_xxxxx_zzzzz_0, g_xxxxx_zzzzz_1, g_xxxxxx_xxxx_1, g_xxxxxx_xxxxx_1, g_xxxxxx_xxxxy_1, g_xxxxxx_xxxxz_1, g_xxxxxx_xxxy_1, g_xxxxxx_xxxyy_1, g_xxxxxx_xxxyz_1, g_xxxxxx_xxxz_1, g_xxxxxx_xxxzz_1, g_xxxxxx_xxyy_1, g_xxxxxx_xxyyy_1, g_xxxxxx_xxyyz_1, g_xxxxxx_xxyz_1, g_xxxxxx_xxyzz_1, g_xxxxxx_xxzz_1, g_xxxxxx_xxzzz_1, g_xxxxxx_xyyy_1, g_xxxxxx_xyyyy_1, g_xxxxxx_xyyyz_1, g_xxxxxx_xyyz_1, g_xxxxxx_xyyzz_1, g_xxxxxx_xyzz_1, g_xxxxxx_xyzzz_1, g_xxxxxx_xzzz_1, g_xxxxxx_xzzzz_1, g_xxxxxx_yyyy_1, g_xxxxxx_yyyyy_1, g_xxxxxx_yyyyz_1, g_xxxxxx_yyyz_1, g_xxxxxx_yyyzz_1, g_xxxxxx_yyzz_1, g_xxxxxx_yyzzz_1, g_xxxxxx_yzzz_1, g_xxxxxx_yzzzz_1, g_xxxxxx_zzzz_1, g_xxxxxx_zzzzz_1, g_xxxxxxx_xxxxx_0, g_xxxxxxx_xxxxy_0, g_xxxxxxx_xxxxz_0, g_xxxxxxx_xxxyy_0, g_xxxxxxx_xxxyz_0, g_xxxxxxx_xxxzz_0, g_xxxxxxx_xxyyy_0, g_xxxxxxx_xxyyz_0, g_xxxxxxx_xxyzz_0, g_xxxxxxx_xxzzz_0, g_xxxxxxx_xyyyy_0, g_xxxxxxx_xyyyz_0, g_xxxxxxx_xyyzz_0, g_xxxxxxx_xyzzz_0, g_xxxxxxx_xzzzz_0, g_xxxxxxx_yyyyy_0, g_xxxxxxx_yyyyz_0, g_xxxxxxx_yyyzz_0, g_xxxxxxx_yyzzz_0, g_xxxxxxx_yzzzz_0, g_xxxxxxx_zzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxxx_xxxxx_0[i] = 6.0 * g_xxxxx_xxxxx_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxxx_1[i] * fz_be_0 + 5.0 * g_xxxxxx_xxxx_1[i] * fe_0 + g_xxxxxx_xxxxx_1[i] * pa_x[i];

        g_xxxxxxx_xxxxy_0[i] = 6.0 * g_xxxxx_xxxxy_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxxxxx_xxxy_1[i] * fe_0 + g_xxxxxx_xxxxy_1[i] * pa_x[i];

        g_xxxxxxx_xxxxz_0[i] = 6.0 * g_xxxxx_xxxxz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_xxxz_1[i] * fe_0 + g_xxxxxx_xxxxz_1[i] * pa_x[i];

        g_xxxxxxx_xxxyy_0[i] = 6.0 * g_xxxxx_xxxyy_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxxxxx_xxyy_1[i] * fe_0 + g_xxxxxx_xxxyy_1[i] * pa_x[i];

        g_xxxxxxx_xxxyz_0[i] = 6.0 * g_xxxxx_xxxyz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_xxyz_1[i] * fe_0 + g_xxxxxx_xxxyz_1[i] * pa_x[i];

        g_xxxxxxx_xxxzz_0[i] = 6.0 * g_xxxxx_xxxzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_xxzz_1[i] * fe_0 + g_xxxxxx_xxxzz_1[i] * pa_x[i];

        g_xxxxxxx_xxyyy_0[i] = 6.0 * g_xxxxx_xxyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xyyy_1[i] * fe_0 + g_xxxxxx_xxyyy_1[i] * pa_x[i];

        g_xxxxxxx_xxyyz_0[i] = 6.0 * g_xxxxx_xxyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xyyz_1[i] * fe_0 + g_xxxxxx_xxyyz_1[i] * pa_x[i];

        g_xxxxxxx_xxyzz_0[i] = 6.0 * g_xxxxx_xxyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xyzz_1[i] * fe_0 + g_xxxxxx_xxyzz_1[i] * pa_x[i];

        g_xxxxxxx_xxzzz_0[i] = 6.0 * g_xxxxx_xxzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xzzz_1[i] * fe_0 + g_xxxxxx_xxzzz_1[i] * pa_x[i];

        g_xxxxxxx_xyyyy_0[i] = 6.0 * g_xxxxx_xyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_xyyyy_1[i] * fz_be_0 + g_xxxxxx_yyyy_1[i] * fe_0 + g_xxxxxx_xyyyy_1[i] * pa_x[i];

        g_xxxxxxx_xyyyz_0[i] = 6.0 * g_xxxxx_xyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_xyyyz_1[i] * fz_be_0 + g_xxxxxx_yyyz_1[i] * fe_0 + g_xxxxxx_xyyyz_1[i] * pa_x[i];

        g_xxxxxxx_xyyzz_0[i] = 6.0 * g_xxxxx_xyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xyyzz_1[i] * fz_be_0 + g_xxxxxx_yyzz_1[i] * fe_0 + g_xxxxxx_xyyzz_1[i] * pa_x[i];

        g_xxxxxxx_xyzzz_0[i] = 6.0 * g_xxxxx_xyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xyzzz_1[i] * fz_be_0 + g_xxxxxx_yzzz_1[i] * fe_0 + g_xxxxxx_xyzzz_1[i] * pa_x[i];

        g_xxxxxxx_xzzzz_0[i] = 6.0 * g_xxxxx_xzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xzzzz_1[i] * fz_be_0 + g_xxxxxx_zzzz_1[i] * fe_0 + g_xxxxxx_xzzzz_1[i] * pa_x[i];

        g_xxxxxxx_yyyyy_0[i] = 6.0 * g_xxxxx_yyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_yyyyy_1[i] * fz_be_0 + g_xxxxxx_yyyyy_1[i] * pa_x[i];

        g_xxxxxxx_yyyyz_0[i] = 6.0 * g_xxxxx_yyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_yyyyz_1[i] * fz_be_0 + g_xxxxxx_yyyyz_1[i] * pa_x[i];

        g_xxxxxxx_yyyzz_0[i] = 6.0 * g_xxxxx_yyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_yyyzz_1[i] * fz_be_0 + g_xxxxxx_yyyzz_1[i] * pa_x[i];

        g_xxxxxxx_yyzzz_0[i] = 6.0 * g_xxxxx_yyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_yyzzz_1[i] * fz_be_0 + g_xxxxxx_yyzzz_1[i] * pa_x[i];

        g_xxxxxxx_yzzzz_0[i] = 6.0 * g_xxxxx_yzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_yzzzz_1[i] * fz_be_0 + g_xxxxxx_yzzzz_1[i] * pa_x[i];

        g_xxxxxxx_zzzzz_0[i] = 6.0 * g_xxxxx_zzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_zzzzz_1[i] * fz_be_0 + g_xxxxxx_zzzzz_1[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : KH

    auto g_xxxxxxy_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 21);

    auto g_xxxxxxy_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 22);

    auto g_xxxxxxy_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 23);

    auto g_xxxxxxy_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 24);

    auto g_xxxxxxy_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 25);

    auto g_xxxxxxy_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 26);

    auto g_xxxxxxy_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 27);

    auto g_xxxxxxy_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 28);

    auto g_xxxxxxy_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 29);

    auto g_xxxxxxy_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 30);

    auto g_xxxxxxy_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 31);

    auto g_xxxxxxy_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 32);

    auto g_xxxxxxy_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 33);

    auto g_xxxxxxy_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 34);

    auto g_xxxxxxy_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 35);

    auto g_xxxxxxy_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 36);

    auto g_xxxxxxy_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 37);

    auto g_xxxxxxy_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 38);

    auto g_xxxxxxy_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 39);

    auto g_xxxxxxy_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 40);

    auto g_xxxxxxy_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 41);

    #pragma omp simd aligned(g_xxxxxx_xxxx_1, g_xxxxxx_xxxxx_1, g_xxxxxx_xxxxy_1, g_xxxxxx_xxxxz_1, g_xxxxxx_xxxy_1, g_xxxxxx_xxxyy_1, g_xxxxxx_xxxyz_1, g_xxxxxx_xxxz_1, g_xxxxxx_xxxzz_1, g_xxxxxx_xxyy_1, g_xxxxxx_xxyyy_1, g_xxxxxx_xxyyz_1, g_xxxxxx_xxyz_1, g_xxxxxx_xxyzz_1, g_xxxxxx_xxzz_1, g_xxxxxx_xxzzz_1, g_xxxxxx_xyyy_1, g_xxxxxx_xyyyy_1, g_xxxxxx_xyyyz_1, g_xxxxxx_xyyz_1, g_xxxxxx_xyyzz_1, g_xxxxxx_xyzz_1, g_xxxxxx_xyzzz_1, g_xxxxxx_xzzz_1, g_xxxxxx_xzzzz_1, g_xxxxxx_yyyy_1, g_xxxxxx_yyyyy_1, g_xxxxxx_yyyyz_1, g_xxxxxx_yyyz_1, g_xxxxxx_yyyzz_1, g_xxxxxx_yyzz_1, g_xxxxxx_yyzzz_1, g_xxxxxx_yzzz_1, g_xxxxxx_yzzzz_1, g_xxxxxx_zzzz_1, g_xxxxxx_zzzzz_1, g_xxxxxxy_xxxxx_0, g_xxxxxxy_xxxxy_0, g_xxxxxxy_xxxxz_0, g_xxxxxxy_xxxyy_0, g_xxxxxxy_xxxyz_0, g_xxxxxxy_xxxzz_0, g_xxxxxxy_xxyyy_0, g_xxxxxxy_xxyyz_0, g_xxxxxxy_xxyzz_0, g_xxxxxxy_xxzzz_0, g_xxxxxxy_xyyyy_0, g_xxxxxxy_xyyyz_0, g_xxxxxxy_xyyzz_0, g_xxxxxxy_xyzzz_0, g_xxxxxxy_xzzzz_0, g_xxxxxxy_yyyyy_0, g_xxxxxxy_yyyyz_0, g_xxxxxxy_yyyzz_0, g_xxxxxxy_yyzzz_0, g_xxxxxxy_yzzzz_0, g_xxxxxxy_zzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxxy_xxxxx_0[i] = g_xxxxxx_xxxxx_1[i] * pa_y[i];

        g_xxxxxxy_xxxxy_0[i] = g_xxxxxx_xxxx_1[i] * fe_0 + g_xxxxxx_xxxxy_1[i] * pa_y[i];

        g_xxxxxxy_xxxxz_0[i] = g_xxxxxx_xxxxz_1[i] * pa_y[i];

        g_xxxxxxy_xxxyy_0[i] = 2.0 * g_xxxxxx_xxxy_1[i] * fe_0 + g_xxxxxx_xxxyy_1[i] * pa_y[i];

        g_xxxxxxy_xxxyz_0[i] = g_xxxxxx_xxxz_1[i] * fe_0 + g_xxxxxx_xxxyz_1[i] * pa_y[i];

        g_xxxxxxy_xxxzz_0[i] = g_xxxxxx_xxxzz_1[i] * pa_y[i];

        g_xxxxxxy_xxyyy_0[i] = 3.0 * g_xxxxxx_xxyy_1[i] * fe_0 + g_xxxxxx_xxyyy_1[i] * pa_y[i];

        g_xxxxxxy_xxyyz_0[i] = 2.0 * g_xxxxxx_xxyz_1[i] * fe_0 + g_xxxxxx_xxyyz_1[i] * pa_y[i];

        g_xxxxxxy_xxyzz_0[i] = g_xxxxxx_xxzz_1[i] * fe_0 + g_xxxxxx_xxyzz_1[i] * pa_y[i];

        g_xxxxxxy_xxzzz_0[i] = g_xxxxxx_xxzzz_1[i] * pa_y[i];

        g_xxxxxxy_xyyyy_0[i] = 4.0 * g_xxxxxx_xyyy_1[i] * fe_0 + g_xxxxxx_xyyyy_1[i] * pa_y[i];

        g_xxxxxxy_xyyyz_0[i] = 3.0 * g_xxxxxx_xyyz_1[i] * fe_0 + g_xxxxxx_xyyyz_1[i] * pa_y[i];

        g_xxxxxxy_xyyzz_0[i] = 2.0 * g_xxxxxx_xyzz_1[i] * fe_0 + g_xxxxxx_xyyzz_1[i] * pa_y[i];

        g_xxxxxxy_xyzzz_0[i] = g_xxxxxx_xzzz_1[i] * fe_0 + g_xxxxxx_xyzzz_1[i] * pa_y[i];

        g_xxxxxxy_xzzzz_0[i] = g_xxxxxx_xzzzz_1[i] * pa_y[i];

        g_xxxxxxy_yyyyy_0[i] = 5.0 * g_xxxxxx_yyyy_1[i] * fe_0 + g_xxxxxx_yyyyy_1[i] * pa_y[i];

        g_xxxxxxy_yyyyz_0[i] = 4.0 * g_xxxxxx_yyyz_1[i] * fe_0 + g_xxxxxx_yyyyz_1[i] * pa_y[i];

        g_xxxxxxy_yyyzz_0[i] = 3.0 * g_xxxxxx_yyzz_1[i] * fe_0 + g_xxxxxx_yyyzz_1[i] * pa_y[i];

        g_xxxxxxy_yyzzz_0[i] = 2.0 * g_xxxxxx_yzzz_1[i] * fe_0 + g_xxxxxx_yyzzz_1[i] * pa_y[i];

        g_xxxxxxy_yzzzz_0[i] = g_xxxxxx_zzzz_1[i] * fe_0 + g_xxxxxx_yzzzz_1[i] * pa_y[i];

        g_xxxxxxy_zzzzz_0[i] = g_xxxxxx_zzzzz_1[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : KH

    auto g_xxxxxxz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 42);

    auto g_xxxxxxz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 43);

    auto g_xxxxxxz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 44);

    auto g_xxxxxxz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 45);

    auto g_xxxxxxz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 46);

    auto g_xxxxxxz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 47);

    auto g_xxxxxxz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 48);

    auto g_xxxxxxz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 49);

    auto g_xxxxxxz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 50);

    auto g_xxxxxxz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 51);

    auto g_xxxxxxz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 52);

    auto g_xxxxxxz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 53);

    auto g_xxxxxxz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 54);

    auto g_xxxxxxz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 55);

    auto g_xxxxxxz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 56);

    auto g_xxxxxxz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 57);

    auto g_xxxxxxz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 58);

    auto g_xxxxxxz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 59);

    auto g_xxxxxxz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 60);

    auto g_xxxxxxz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 61);

    auto g_xxxxxxz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 62);

    #pragma omp simd aligned(g_xxxxxx_xxxx_1, g_xxxxxx_xxxxx_1, g_xxxxxx_xxxxy_1, g_xxxxxx_xxxxz_1, g_xxxxxx_xxxy_1, g_xxxxxx_xxxyy_1, g_xxxxxx_xxxyz_1, g_xxxxxx_xxxz_1, g_xxxxxx_xxxzz_1, g_xxxxxx_xxyy_1, g_xxxxxx_xxyyy_1, g_xxxxxx_xxyyz_1, g_xxxxxx_xxyz_1, g_xxxxxx_xxyzz_1, g_xxxxxx_xxzz_1, g_xxxxxx_xxzzz_1, g_xxxxxx_xyyy_1, g_xxxxxx_xyyyy_1, g_xxxxxx_xyyyz_1, g_xxxxxx_xyyz_1, g_xxxxxx_xyyzz_1, g_xxxxxx_xyzz_1, g_xxxxxx_xyzzz_1, g_xxxxxx_xzzz_1, g_xxxxxx_xzzzz_1, g_xxxxxx_yyyy_1, g_xxxxxx_yyyyy_1, g_xxxxxx_yyyyz_1, g_xxxxxx_yyyz_1, g_xxxxxx_yyyzz_1, g_xxxxxx_yyzz_1, g_xxxxxx_yyzzz_1, g_xxxxxx_yzzz_1, g_xxxxxx_yzzzz_1, g_xxxxxx_zzzz_1, g_xxxxxx_zzzzz_1, g_xxxxxxz_xxxxx_0, g_xxxxxxz_xxxxy_0, g_xxxxxxz_xxxxz_0, g_xxxxxxz_xxxyy_0, g_xxxxxxz_xxxyz_0, g_xxxxxxz_xxxzz_0, g_xxxxxxz_xxyyy_0, g_xxxxxxz_xxyyz_0, g_xxxxxxz_xxyzz_0, g_xxxxxxz_xxzzz_0, g_xxxxxxz_xyyyy_0, g_xxxxxxz_xyyyz_0, g_xxxxxxz_xyyzz_0, g_xxxxxxz_xyzzz_0, g_xxxxxxz_xzzzz_0, g_xxxxxxz_yyyyy_0, g_xxxxxxz_yyyyz_0, g_xxxxxxz_yyyzz_0, g_xxxxxxz_yyzzz_0, g_xxxxxxz_yzzzz_0, g_xxxxxxz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxxz_xxxxx_0[i] = g_xxxxxx_xxxxx_1[i] * pa_z[i];

        g_xxxxxxz_xxxxy_0[i] = g_xxxxxx_xxxxy_1[i] * pa_z[i];

        g_xxxxxxz_xxxxz_0[i] = g_xxxxxx_xxxx_1[i] * fe_0 + g_xxxxxx_xxxxz_1[i] * pa_z[i];

        g_xxxxxxz_xxxyy_0[i] = g_xxxxxx_xxxyy_1[i] * pa_z[i];

        g_xxxxxxz_xxxyz_0[i] = g_xxxxxx_xxxy_1[i] * fe_0 + g_xxxxxx_xxxyz_1[i] * pa_z[i];

        g_xxxxxxz_xxxzz_0[i] = 2.0 * g_xxxxxx_xxxz_1[i] * fe_0 + g_xxxxxx_xxxzz_1[i] * pa_z[i];

        g_xxxxxxz_xxyyy_0[i] = g_xxxxxx_xxyyy_1[i] * pa_z[i];

        g_xxxxxxz_xxyyz_0[i] = g_xxxxxx_xxyy_1[i] * fe_0 + g_xxxxxx_xxyyz_1[i] * pa_z[i];

        g_xxxxxxz_xxyzz_0[i] = 2.0 * g_xxxxxx_xxyz_1[i] * fe_0 + g_xxxxxx_xxyzz_1[i] * pa_z[i];

        g_xxxxxxz_xxzzz_0[i] = 3.0 * g_xxxxxx_xxzz_1[i] * fe_0 + g_xxxxxx_xxzzz_1[i] * pa_z[i];

        g_xxxxxxz_xyyyy_0[i] = g_xxxxxx_xyyyy_1[i] * pa_z[i];

        g_xxxxxxz_xyyyz_0[i] = g_xxxxxx_xyyy_1[i] * fe_0 + g_xxxxxx_xyyyz_1[i] * pa_z[i];

        g_xxxxxxz_xyyzz_0[i] = 2.0 * g_xxxxxx_xyyz_1[i] * fe_0 + g_xxxxxx_xyyzz_1[i] * pa_z[i];

        g_xxxxxxz_xyzzz_0[i] = 3.0 * g_xxxxxx_xyzz_1[i] * fe_0 + g_xxxxxx_xyzzz_1[i] * pa_z[i];

        g_xxxxxxz_xzzzz_0[i] = 4.0 * g_xxxxxx_xzzz_1[i] * fe_0 + g_xxxxxx_xzzzz_1[i] * pa_z[i];

        g_xxxxxxz_yyyyy_0[i] = g_xxxxxx_yyyyy_1[i] * pa_z[i];

        g_xxxxxxz_yyyyz_0[i] = g_xxxxxx_yyyy_1[i] * fe_0 + g_xxxxxx_yyyyz_1[i] * pa_z[i];

        g_xxxxxxz_yyyzz_0[i] = 2.0 * g_xxxxxx_yyyz_1[i] * fe_0 + g_xxxxxx_yyyzz_1[i] * pa_z[i];

        g_xxxxxxz_yyzzz_0[i] = 3.0 * g_xxxxxx_yyzz_1[i] * fe_0 + g_xxxxxx_yyzzz_1[i] * pa_z[i];

        g_xxxxxxz_yzzzz_0[i] = 4.0 * g_xxxxxx_yzzz_1[i] * fe_0 + g_xxxxxx_yzzzz_1[i] * pa_z[i];

        g_xxxxxxz_zzzzz_0[i] = 5.0 * g_xxxxxx_zzzz_1[i] * fe_0 + g_xxxxxx_zzzzz_1[i] * pa_z[i];
    }

    // Set up 63-84 components of targeted buffer : KH

    auto g_xxxxxyy_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 63);

    auto g_xxxxxyy_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 64);

    auto g_xxxxxyy_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 65);

    auto g_xxxxxyy_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 66);

    auto g_xxxxxyy_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 67);

    auto g_xxxxxyy_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 68);

    auto g_xxxxxyy_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 69);

    auto g_xxxxxyy_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 70);

    auto g_xxxxxyy_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 71);

    auto g_xxxxxyy_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 72);

    auto g_xxxxxyy_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 73);

    auto g_xxxxxyy_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 74);

    auto g_xxxxxyy_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 75);

    auto g_xxxxxyy_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 76);

    auto g_xxxxxyy_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 77);

    auto g_xxxxxyy_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 78);

    auto g_xxxxxyy_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 79);

    auto g_xxxxxyy_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 80);

    auto g_xxxxxyy_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 81);

    auto g_xxxxxyy_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 82);

    auto g_xxxxxyy_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 83);

    #pragma omp simd aligned(g_xxxxx_xxxxx_0, g_xxxxx_xxxxx_1, g_xxxxx_xxxxz_0, g_xxxxx_xxxxz_1, g_xxxxx_xxxzz_0, g_xxxxx_xxxzz_1, g_xxxxx_xxzzz_0, g_xxxxx_xxzzz_1, g_xxxxx_xzzzz_0, g_xxxxx_xzzzz_1, g_xxxxxy_xxxxx_1, g_xxxxxy_xxxxz_1, g_xxxxxy_xxxzz_1, g_xxxxxy_xxzzz_1, g_xxxxxy_xzzzz_1, g_xxxxxyy_xxxxx_0, g_xxxxxyy_xxxxy_0, g_xxxxxyy_xxxxz_0, g_xxxxxyy_xxxyy_0, g_xxxxxyy_xxxyz_0, g_xxxxxyy_xxxzz_0, g_xxxxxyy_xxyyy_0, g_xxxxxyy_xxyyz_0, g_xxxxxyy_xxyzz_0, g_xxxxxyy_xxzzz_0, g_xxxxxyy_xyyyy_0, g_xxxxxyy_xyyyz_0, g_xxxxxyy_xyyzz_0, g_xxxxxyy_xyzzz_0, g_xxxxxyy_xzzzz_0, g_xxxxxyy_yyyyy_0, g_xxxxxyy_yyyyz_0, g_xxxxxyy_yyyzz_0, g_xxxxxyy_yyzzz_0, g_xxxxxyy_yzzzz_0, g_xxxxxyy_zzzzz_0, g_xxxxyy_xxxxy_1, g_xxxxyy_xxxy_1, g_xxxxyy_xxxyy_1, g_xxxxyy_xxxyz_1, g_xxxxyy_xxyy_1, g_xxxxyy_xxyyy_1, g_xxxxyy_xxyyz_1, g_xxxxyy_xxyz_1, g_xxxxyy_xxyzz_1, g_xxxxyy_xyyy_1, g_xxxxyy_xyyyy_1, g_xxxxyy_xyyyz_1, g_xxxxyy_xyyz_1, g_xxxxyy_xyyzz_1, g_xxxxyy_xyzz_1, g_xxxxyy_xyzzz_1, g_xxxxyy_yyyy_1, g_xxxxyy_yyyyy_1, g_xxxxyy_yyyyz_1, g_xxxxyy_yyyz_1, g_xxxxyy_yyyzz_1, g_xxxxyy_yyzz_1, g_xxxxyy_yyzzz_1, g_xxxxyy_yzzz_1, g_xxxxyy_yzzzz_1, g_xxxxyy_zzzzz_1, g_xxxyy_xxxxy_0, g_xxxyy_xxxxy_1, g_xxxyy_xxxyy_0, g_xxxyy_xxxyy_1, g_xxxyy_xxxyz_0, g_xxxyy_xxxyz_1, g_xxxyy_xxyyy_0, g_xxxyy_xxyyy_1, g_xxxyy_xxyyz_0, g_xxxyy_xxyyz_1, g_xxxyy_xxyzz_0, g_xxxyy_xxyzz_1, g_xxxyy_xyyyy_0, g_xxxyy_xyyyy_1, g_xxxyy_xyyyz_0, g_xxxyy_xyyyz_1, g_xxxyy_xyyzz_0, g_xxxyy_xyyzz_1, g_xxxyy_xyzzz_0, g_xxxyy_xyzzz_1, g_xxxyy_yyyyy_0, g_xxxyy_yyyyy_1, g_xxxyy_yyyyz_0, g_xxxyy_yyyyz_1, g_xxxyy_yyyzz_0, g_xxxyy_yyyzz_1, g_xxxyy_yyzzz_0, g_xxxyy_yyzzz_1, g_xxxyy_yzzzz_0, g_xxxyy_yzzzz_1, g_xxxyy_zzzzz_0, g_xxxyy_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxyy_xxxxx_0[i] = g_xxxxx_xxxxx_0[i] * fbe_0 - g_xxxxx_xxxxx_1[i] * fz_be_0 + g_xxxxxy_xxxxx_1[i] * pa_y[i];

        g_xxxxxyy_xxxxy_0[i] = 4.0 * g_xxxyy_xxxxy_0[i] * fbe_0 - 4.0 * g_xxxyy_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxxxyy_xxxy_1[i] * fe_0 + g_xxxxyy_xxxxy_1[i] * pa_x[i];

        g_xxxxxyy_xxxxz_0[i] = g_xxxxx_xxxxz_0[i] * fbe_0 - g_xxxxx_xxxxz_1[i] * fz_be_0 + g_xxxxxy_xxxxz_1[i] * pa_y[i];

        g_xxxxxyy_xxxyy_0[i] = 4.0 * g_xxxyy_xxxyy_0[i] * fbe_0 - 4.0 * g_xxxyy_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxxxyy_xxyy_1[i] * fe_0 + g_xxxxyy_xxxyy_1[i] * pa_x[i];

        g_xxxxxyy_xxxyz_0[i] = 4.0 * g_xxxyy_xxxyz_0[i] * fbe_0 - 4.0 * g_xxxyy_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_xxyz_1[i] * fe_0 + g_xxxxyy_xxxyz_1[i] * pa_x[i];

        g_xxxxxyy_xxxzz_0[i] = g_xxxxx_xxxzz_0[i] * fbe_0 - g_xxxxx_xxxzz_1[i] * fz_be_0 + g_xxxxxy_xxxzz_1[i] * pa_y[i];

        g_xxxxxyy_xxyyy_0[i] = 4.0 * g_xxxyy_xxyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxxxyy_xyyy_1[i] * fe_0 + g_xxxxyy_xxyyy_1[i] * pa_x[i];

        g_xxxxxyy_xxyyz_0[i] = 4.0 * g_xxxyy_xxyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_xyyz_1[i] * fe_0 + g_xxxxyy_xxyyz_1[i] * pa_x[i];

        g_xxxxxyy_xxyzz_0[i] = 4.0 * g_xxxyy_xxyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_xyzz_1[i] * fe_0 + g_xxxxyy_xxyzz_1[i] * pa_x[i];

        g_xxxxxyy_xxzzz_0[i] = g_xxxxx_xxzzz_0[i] * fbe_0 - g_xxxxx_xxzzz_1[i] * fz_be_0 + g_xxxxxy_xxzzz_1[i] * pa_y[i];

        g_xxxxxyy_xyyyy_0[i] = 4.0 * g_xxxyy_xyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_xyyyy_1[i] * fz_be_0 + g_xxxxyy_yyyy_1[i] * fe_0 + g_xxxxyy_xyyyy_1[i] * pa_x[i];

        g_xxxxxyy_xyyyz_0[i] = 4.0 * g_xxxyy_xyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_xyyyz_1[i] * fz_be_0 + g_xxxxyy_yyyz_1[i] * fe_0 + g_xxxxyy_xyyyz_1[i] * pa_x[i];

        g_xxxxxyy_xyyzz_0[i] = 4.0 * g_xxxyy_xyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_xyyzz_1[i] * fz_be_0 + g_xxxxyy_yyzz_1[i] * fe_0 + g_xxxxyy_xyyzz_1[i] * pa_x[i];

        g_xxxxxyy_xyzzz_0[i] = 4.0 * g_xxxyy_xyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_xyzzz_1[i] * fz_be_0 + g_xxxxyy_yzzz_1[i] * fe_0 + g_xxxxyy_xyzzz_1[i] * pa_x[i];

        g_xxxxxyy_xzzzz_0[i] = g_xxxxx_xzzzz_0[i] * fbe_0 - g_xxxxx_xzzzz_1[i] * fz_be_0 + g_xxxxxy_xzzzz_1[i] * pa_y[i];

        g_xxxxxyy_yyyyy_0[i] = 4.0 * g_xxxyy_yyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_yyyyy_1[i] * fz_be_0 + g_xxxxyy_yyyyy_1[i] * pa_x[i];

        g_xxxxxyy_yyyyz_0[i] = 4.0 * g_xxxyy_yyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_yyyyz_1[i] * fz_be_0 + g_xxxxyy_yyyyz_1[i] * pa_x[i];

        g_xxxxxyy_yyyzz_0[i] = 4.0 * g_xxxyy_yyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_yyyzz_1[i] * fz_be_0 + g_xxxxyy_yyyzz_1[i] * pa_x[i];

        g_xxxxxyy_yyzzz_0[i] = 4.0 * g_xxxyy_yyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_yyzzz_1[i] * fz_be_0 + g_xxxxyy_yyzzz_1[i] * pa_x[i];

        g_xxxxxyy_yzzzz_0[i] = 4.0 * g_xxxyy_yzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_yzzzz_1[i] * fz_be_0 + g_xxxxyy_yzzzz_1[i] * pa_x[i];

        g_xxxxxyy_zzzzz_0[i] = 4.0 * g_xxxyy_zzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_zzzzz_1[i] * fz_be_0 + g_xxxxyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 84-105 components of targeted buffer : KH

    auto g_xxxxxyz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 84);

    auto g_xxxxxyz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 85);

    auto g_xxxxxyz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 86);

    auto g_xxxxxyz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 87);

    auto g_xxxxxyz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 88);

    auto g_xxxxxyz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 89);

    auto g_xxxxxyz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 90);

    auto g_xxxxxyz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 91);

    auto g_xxxxxyz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 92);

    auto g_xxxxxyz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 93);

    auto g_xxxxxyz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 94);

    auto g_xxxxxyz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 95);

    auto g_xxxxxyz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 96);

    auto g_xxxxxyz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 97);

    auto g_xxxxxyz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 98);

    auto g_xxxxxyz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 99);

    auto g_xxxxxyz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 100);

    auto g_xxxxxyz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 101);

    auto g_xxxxxyz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 102);

    auto g_xxxxxyz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 103);

    auto g_xxxxxyz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 104);

    #pragma omp simd aligned(g_xxxxxy_xxxxy_1, g_xxxxxy_xxxyy_1, g_xxxxxy_xxyyy_1, g_xxxxxy_xyyyy_1, g_xxxxxy_yyyyy_1, g_xxxxxyz_xxxxx_0, g_xxxxxyz_xxxxy_0, g_xxxxxyz_xxxxz_0, g_xxxxxyz_xxxyy_0, g_xxxxxyz_xxxyz_0, g_xxxxxyz_xxxzz_0, g_xxxxxyz_xxyyy_0, g_xxxxxyz_xxyyz_0, g_xxxxxyz_xxyzz_0, g_xxxxxyz_xxzzz_0, g_xxxxxyz_xyyyy_0, g_xxxxxyz_xyyyz_0, g_xxxxxyz_xyyzz_0, g_xxxxxyz_xyzzz_0, g_xxxxxyz_xzzzz_0, g_xxxxxyz_yyyyy_0, g_xxxxxyz_yyyyz_0, g_xxxxxyz_yyyzz_0, g_xxxxxyz_yyzzz_0, g_xxxxxyz_yzzzz_0, g_xxxxxyz_zzzzz_0, g_xxxxxz_xxxxx_1, g_xxxxxz_xxxxz_1, g_xxxxxz_xxxyz_1, g_xxxxxz_xxxz_1, g_xxxxxz_xxxzz_1, g_xxxxxz_xxyyz_1, g_xxxxxz_xxyz_1, g_xxxxxz_xxyzz_1, g_xxxxxz_xxzz_1, g_xxxxxz_xxzzz_1, g_xxxxxz_xyyyz_1, g_xxxxxz_xyyz_1, g_xxxxxz_xyyzz_1, g_xxxxxz_xyzz_1, g_xxxxxz_xyzzz_1, g_xxxxxz_xzzz_1, g_xxxxxz_xzzzz_1, g_xxxxxz_yyyyz_1, g_xxxxxz_yyyz_1, g_xxxxxz_yyyzz_1, g_xxxxxz_yyzz_1, g_xxxxxz_yyzzz_1, g_xxxxxz_yzzz_1, g_xxxxxz_yzzzz_1, g_xxxxxz_zzzz_1, g_xxxxxz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxyz_xxxxx_0[i] = g_xxxxxz_xxxxx_1[i] * pa_y[i];

        g_xxxxxyz_xxxxy_0[i] = g_xxxxxy_xxxxy_1[i] * pa_z[i];

        g_xxxxxyz_xxxxz_0[i] = g_xxxxxz_xxxxz_1[i] * pa_y[i];

        g_xxxxxyz_xxxyy_0[i] = g_xxxxxy_xxxyy_1[i] * pa_z[i];

        g_xxxxxyz_xxxyz_0[i] = g_xxxxxz_xxxz_1[i] * fe_0 + g_xxxxxz_xxxyz_1[i] * pa_y[i];

        g_xxxxxyz_xxxzz_0[i] = g_xxxxxz_xxxzz_1[i] * pa_y[i];

        g_xxxxxyz_xxyyy_0[i] = g_xxxxxy_xxyyy_1[i] * pa_z[i];

        g_xxxxxyz_xxyyz_0[i] = 2.0 * g_xxxxxz_xxyz_1[i] * fe_0 + g_xxxxxz_xxyyz_1[i] * pa_y[i];

        g_xxxxxyz_xxyzz_0[i] = g_xxxxxz_xxzz_1[i] * fe_0 + g_xxxxxz_xxyzz_1[i] * pa_y[i];

        g_xxxxxyz_xxzzz_0[i] = g_xxxxxz_xxzzz_1[i] * pa_y[i];

        g_xxxxxyz_xyyyy_0[i] = g_xxxxxy_xyyyy_1[i] * pa_z[i];

        g_xxxxxyz_xyyyz_0[i] = 3.0 * g_xxxxxz_xyyz_1[i] * fe_0 + g_xxxxxz_xyyyz_1[i] * pa_y[i];

        g_xxxxxyz_xyyzz_0[i] = 2.0 * g_xxxxxz_xyzz_1[i] * fe_0 + g_xxxxxz_xyyzz_1[i] * pa_y[i];

        g_xxxxxyz_xyzzz_0[i] = g_xxxxxz_xzzz_1[i] * fe_0 + g_xxxxxz_xyzzz_1[i] * pa_y[i];

        g_xxxxxyz_xzzzz_0[i] = g_xxxxxz_xzzzz_1[i] * pa_y[i];

        g_xxxxxyz_yyyyy_0[i] = g_xxxxxy_yyyyy_1[i] * pa_z[i];

        g_xxxxxyz_yyyyz_0[i] = 4.0 * g_xxxxxz_yyyz_1[i] * fe_0 + g_xxxxxz_yyyyz_1[i] * pa_y[i];

        g_xxxxxyz_yyyzz_0[i] = 3.0 * g_xxxxxz_yyzz_1[i] * fe_0 + g_xxxxxz_yyyzz_1[i] * pa_y[i];

        g_xxxxxyz_yyzzz_0[i] = 2.0 * g_xxxxxz_yzzz_1[i] * fe_0 + g_xxxxxz_yyzzz_1[i] * pa_y[i];

        g_xxxxxyz_yzzzz_0[i] = g_xxxxxz_zzzz_1[i] * fe_0 + g_xxxxxz_yzzzz_1[i] * pa_y[i];

        g_xxxxxyz_zzzzz_0[i] = g_xxxxxz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : KH

    auto g_xxxxxzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 105);

    auto g_xxxxxzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 106);

    auto g_xxxxxzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 107);

    auto g_xxxxxzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 108);

    auto g_xxxxxzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 109);

    auto g_xxxxxzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 110);

    auto g_xxxxxzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 111);

    auto g_xxxxxzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 112);

    auto g_xxxxxzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 113);

    auto g_xxxxxzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 114);

    auto g_xxxxxzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 115);

    auto g_xxxxxzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 116);

    auto g_xxxxxzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 117);

    auto g_xxxxxzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 118);

    auto g_xxxxxzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 119);

    auto g_xxxxxzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 120);

    auto g_xxxxxzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 121);

    auto g_xxxxxzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 122);

    auto g_xxxxxzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 123);

    auto g_xxxxxzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 124);

    auto g_xxxxxzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 125);

    #pragma omp simd aligned(g_xxxxx_xxxxx_0, g_xxxxx_xxxxx_1, g_xxxxx_xxxxy_0, g_xxxxx_xxxxy_1, g_xxxxx_xxxyy_0, g_xxxxx_xxxyy_1, g_xxxxx_xxyyy_0, g_xxxxx_xxyyy_1, g_xxxxx_xyyyy_0, g_xxxxx_xyyyy_1, g_xxxxxz_xxxxx_1, g_xxxxxz_xxxxy_1, g_xxxxxz_xxxyy_1, g_xxxxxz_xxyyy_1, g_xxxxxz_xyyyy_1, g_xxxxxzz_xxxxx_0, g_xxxxxzz_xxxxy_0, g_xxxxxzz_xxxxz_0, g_xxxxxzz_xxxyy_0, g_xxxxxzz_xxxyz_0, g_xxxxxzz_xxxzz_0, g_xxxxxzz_xxyyy_0, g_xxxxxzz_xxyyz_0, g_xxxxxzz_xxyzz_0, g_xxxxxzz_xxzzz_0, g_xxxxxzz_xyyyy_0, g_xxxxxzz_xyyyz_0, g_xxxxxzz_xyyzz_0, g_xxxxxzz_xyzzz_0, g_xxxxxzz_xzzzz_0, g_xxxxxzz_yyyyy_0, g_xxxxxzz_yyyyz_0, g_xxxxxzz_yyyzz_0, g_xxxxxzz_yyzzz_0, g_xxxxxzz_yzzzz_0, g_xxxxxzz_zzzzz_0, g_xxxxzz_xxxxz_1, g_xxxxzz_xxxyz_1, g_xxxxzz_xxxz_1, g_xxxxzz_xxxzz_1, g_xxxxzz_xxyyz_1, g_xxxxzz_xxyz_1, g_xxxxzz_xxyzz_1, g_xxxxzz_xxzz_1, g_xxxxzz_xxzzz_1, g_xxxxzz_xyyyz_1, g_xxxxzz_xyyz_1, g_xxxxzz_xyyzz_1, g_xxxxzz_xyzz_1, g_xxxxzz_xyzzz_1, g_xxxxzz_xzzz_1, g_xxxxzz_xzzzz_1, g_xxxxzz_yyyyy_1, g_xxxxzz_yyyyz_1, g_xxxxzz_yyyz_1, g_xxxxzz_yyyzz_1, g_xxxxzz_yyzz_1, g_xxxxzz_yyzzz_1, g_xxxxzz_yzzz_1, g_xxxxzz_yzzzz_1, g_xxxxzz_zzzz_1, g_xxxxzz_zzzzz_1, g_xxxzz_xxxxz_0, g_xxxzz_xxxxz_1, g_xxxzz_xxxyz_0, g_xxxzz_xxxyz_1, g_xxxzz_xxxzz_0, g_xxxzz_xxxzz_1, g_xxxzz_xxyyz_0, g_xxxzz_xxyyz_1, g_xxxzz_xxyzz_0, g_xxxzz_xxyzz_1, g_xxxzz_xxzzz_0, g_xxxzz_xxzzz_1, g_xxxzz_xyyyz_0, g_xxxzz_xyyyz_1, g_xxxzz_xyyzz_0, g_xxxzz_xyyzz_1, g_xxxzz_xyzzz_0, g_xxxzz_xyzzz_1, g_xxxzz_xzzzz_0, g_xxxzz_xzzzz_1, g_xxxzz_yyyyy_0, g_xxxzz_yyyyy_1, g_xxxzz_yyyyz_0, g_xxxzz_yyyyz_1, g_xxxzz_yyyzz_0, g_xxxzz_yyyzz_1, g_xxxzz_yyzzz_0, g_xxxzz_yyzzz_1, g_xxxzz_yzzzz_0, g_xxxzz_yzzzz_1, g_xxxzz_zzzzz_0, g_xxxzz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxzz_xxxxx_0[i] = g_xxxxx_xxxxx_0[i] * fbe_0 - g_xxxxx_xxxxx_1[i] * fz_be_0 + g_xxxxxz_xxxxx_1[i] * pa_z[i];

        g_xxxxxzz_xxxxy_0[i] = g_xxxxx_xxxxy_0[i] * fbe_0 - g_xxxxx_xxxxy_1[i] * fz_be_0 + g_xxxxxz_xxxxy_1[i] * pa_z[i];

        g_xxxxxzz_xxxxz_0[i] = 4.0 * g_xxxzz_xxxxz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_xxxz_1[i] * fe_0 + g_xxxxzz_xxxxz_1[i] * pa_x[i];

        g_xxxxxzz_xxxyy_0[i] = g_xxxxx_xxxyy_0[i] * fbe_0 - g_xxxxx_xxxyy_1[i] * fz_be_0 + g_xxxxxz_xxxyy_1[i] * pa_z[i];

        g_xxxxxzz_xxxyz_0[i] = 4.0 * g_xxxzz_xxxyz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_xxyz_1[i] * fe_0 + g_xxxxzz_xxxyz_1[i] * pa_x[i];

        g_xxxxxzz_xxxzz_0[i] = 4.0 * g_xxxzz_xxxzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_xxzz_1[i] * fe_0 + g_xxxxzz_xxxzz_1[i] * pa_x[i];

        g_xxxxxzz_xxyyy_0[i] = g_xxxxx_xxyyy_0[i] * fbe_0 - g_xxxxx_xxyyy_1[i] * fz_be_0 + g_xxxxxz_xxyyy_1[i] * pa_z[i];

        g_xxxxxzz_xxyyz_0[i] = 4.0 * g_xxxzz_xxyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_xyyz_1[i] * fe_0 + g_xxxxzz_xxyyz_1[i] * pa_x[i];

        g_xxxxxzz_xxyzz_0[i] = 4.0 * g_xxxzz_xxyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_xyzz_1[i] * fe_0 + g_xxxxzz_xxyzz_1[i] * pa_x[i];

        g_xxxxxzz_xxzzz_0[i] = 4.0 * g_xxxzz_xxzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_xzzz_1[i] * fe_0 + g_xxxxzz_xxzzz_1[i] * pa_x[i];

        g_xxxxxzz_xyyyy_0[i] = g_xxxxx_xyyyy_0[i] * fbe_0 - g_xxxxx_xyyyy_1[i] * fz_be_0 + g_xxxxxz_xyyyy_1[i] * pa_z[i];

        g_xxxxxzz_xyyyz_0[i] = 4.0 * g_xxxzz_xyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_xyyyz_1[i] * fz_be_0 + g_xxxxzz_yyyz_1[i] * fe_0 + g_xxxxzz_xyyyz_1[i] * pa_x[i];

        g_xxxxxzz_xyyzz_0[i] = 4.0 * g_xxxzz_xyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xyyzz_1[i] * fz_be_0 + g_xxxxzz_yyzz_1[i] * fe_0 + g_xxxxzz_xyyzz_1[i] * pa_x[i];

        g_xxxxxzz_xyzzz_0[i] = 4.0 * g_xxxzz_xyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xyzzz_1[i] * fz_be_0 + g_xxxxzz_yzzz_1[i] * fe_0 + g_xxxxzz_xyzzz_1[i] * pa_x[i];

        g_xxxxxzz_xzzzz_0[i] = 4.0 * g_xxxzz_xzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xzzzz_1[i] * fz_be_0 + g_xxxxzz_zzzz_1[i] * fe_0 + g_xxxxzz_xzzzz_1[i] * pa_x[i];

        g_xxxxxzz_yyyyy_0[i] = 4.0 * g_xxxzz_yyyyy_0[i] * fbe_0 - 4.0 * g_xxxzz_yyyyy_1[i] * fz_be_0 + g_xxxxzz_yyyyy_1[i] * pa_x[i];

        g_xxxxxzz_yyyyz_0[i] = 4.0 * g_xxxzz_yyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_yyyyz_1[i] * fz_be_0 + g_xxxxzz_yyyyz_1[i] * pa_x[i];

        g_xxxxxzz_yyyzz_0[i] = 4.0 * g_xxxzz_yyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_yyyzz_1[i] * fz_be_0 + g_xxxxzz_yyyzz_1[i] * pa_x[i];

        g_xxxxxzz_yyzzz_0[i] = 4.0 * g_xxxzz_yyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_yyzzz_1[i] * fz_be_0 + g_xxxxzz_yyzzz_1[i] * pa_x[i];

        g_xxxxxzz_yzzzz_0[i] = 4.0 * g_xxxzz_yzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_yzzzz_1[i] * fz_be_0 + g_xxxxzz_yzzzz_1[i] * pa_x[i];

        g_xxxxxzz_zzzzz_0[i] = 4.0 * g_xxxzz_zzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_zzzzz_1[i] * fz_be_0 + g_xxxxzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 126-147 components of targeted buffer : KH

    auto g_xxxxyyy_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 126);

    auto g_xxxxyyy_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 127);

    auto g_xxxxyyy_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 128);

    auto g_xxxxyyy_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 129);

    auto g_xxxxyyy_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 130);

    auto g_xxxxyyy_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 131);

    auto g_xxxxyyy_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 132);

    auto g_xxxxyyy_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 133);

    auto g_xxxxyyy_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 134);

    auto g_xxxxyyy_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 135);

    auto g_xxxxyyy_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 136);

    auto g_xxxxyyy_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 137);

    auto g_xxxxyyy_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 138);

    auto g_xxxxyyy_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 139);

    auto g_xxxxyyy_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 140);

    auto g_xxxxyyy_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 141);

    auto g_xxxxyyy_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 142);

    auto g_xxxxyyy_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 143);

    auto g_xxxxyyy_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 144);

    auto g_xxxxyyy_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 145);

    auto g_xxxxyyy_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 146);

    #pragma omp simd aligned(g_xxxxy_xxxxx_0, g_xxxxy_xxxxx_1, g_xxxxy_xxxxz_0, g_xxxxy_xxxxz_1, g_xxxxy_xxxzz_0, g_xxxxy_xxxzz_1, g_xxxxy_xxzzz_0, g_xxxxy_xxzzz_1, g_xxxxy_xzzzz_0, g_xxxxy_xzzzz_1, g_xxxxyy_xxxxx_1, g_xxxxyy_xxxxz_1, g_xxxxyy_xxxzz_1, g_xxxxyy_xxzzz_1, g_xxxxyy_xzzzz_1, g_xxxxyyy_xxxxx_0, g_xxxxyyy_xxxxy_0, g_xxxxyyy_xxxxz_0, g_xxxxyyy_xxxyy_0, g_xxxxyyy_xxxyz_0, g_xxxxyyy_xxxzz_0, g_xxxxyyy_xxyyy_0, g_xxxxyyy_xxyyz_0, g_xxxxyyy_xxyzz_0, g_xxxxyyy_xxzzz_0, g_xxxxyyy_xyyyy_0, g_xxxxyyy_xyyyz_0, g_xxxxyyy_xyyzz_0, g_xxxxyyy_xyzzz_0, g_xxxxyyy_xzzzz_0, g_xxxxyyy_yyyyy_0, g_xxxxyyy_yyyyz_0, g_xxxxyyy_yyyzz_0, g_xxxxyyy_yyzzz_0, g_xxxxyyy_yzzzz_0, g_xxxxyyy_zzzzz_0, g_xxxyyy_xxxxy_1, g_xxxyyy_xxxy_1, g_xxxyyy_xxxyy_1, g_xxxyyy_xxxyz_1, g_xxxyyy_xxyy_1, g_xxxyyy_xxyyy_1, g_xxxyyy_xxyyz_1, g_xxxyyy_xxyz_1, g_xxxyyy_xxyzz_1, g_xxxyyy_xyyy_1, g_xxxyyy_xyyyy_1, g_xxxyyy_xyyyz_1, g_xxxyyy_xyyz_1, g_xxxyyy_xyyzz_1, g_xxxyyy_xyzz_1, g_xxxyyy_xyzzz_1, g_xxxyyy_yyyy_1, g_xxxyyy_yyyyy_1, g_xxxyyy_yyyyz_1, g_xxxyyy_yyyz_1, g_xxxyyy_yyyzz_1, g_xxxyyy_yyzz_1, g_xxxyyy_yyzzz_1, g_xxxyyy_yzzz_1, g_xxxyyy_yzzzz_1, g_xxxyyy_zzzzz_1, g_xxyyy_xxxxy_0, g_xxyyy_xxxxy_1, g_xxyyy_xxxyy_0, g_xxyyy_xxxyy_1, g_xxyyy_xxxyz_0, g_xxyyy_xxxyz_1, g_xxyyy_xxyyy_0, g_xxyyy_xxyyy_1, g_xxyyy_xxyyz_0, g_xxyyy_xxyyz_1, g_xxyyy_xxyzz_0, g_xxyyy_xxyzz_1, g_xxyyy_xyyyy_0, g_xxyyy_xyyyy_1, g_xxyyy_xyyyz_0, g_xxyyy_xyyyz_1, g_xxyyy_xyyzz_0, g_xxyyy_xyyzz_1, g_xxyyy_xyzzz_0, g_xxyyy_xyzzz_1, g_xxyyy_yyyyy_0, g_xxyyy_yyyyy_1, g_xxyyy_yyyyz_0, g_xxyyy_yyyyz_1, g_xxyyy_yyyzz_0, g_xxyyy_yyyzz_1, g_xxyyy_yyzzz_0, g_xxyyy_yyzzz_1, g_xxyyy_yzzzz_0, g_xxyyy_yzzzz_1, g_xxyyy_zzzzz_0, g_xxyyy_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxyyy_xxxxx_0[i] = 2.0 * g_xxxxy_xxxxx_0[i] * fbe_0 - 2.0 * g_xxxxy_xxxxx_1[i] * fz_be_0 + g_xxxxyy_xxxxx_1[i] * pa_y[i];

        g_xxxxyyy_xxxxy_0[i] = 3.0 * g_xxyyy_xxxxy_0[i] * fbe_0 - 3.0 * g_xxyyy_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxxyyy_xxxy_1[i] * fe_0 + g_xxxyyy_xxxxy_1[i] * pa_x[i];

        g_xxxxyyy_xxxxz_0[i] = 2.0 * g_xxxxy_xxxxz_0[i] * fbe_0 - 2.0 * g_xxxxy_xxxxz_1[i] * fz_be_0 + g_xxxxyy_xxxxz_1[i] * pa_y[i];

        g_xxxxyyy_xxxyy_0[i] = 3.0 * g_xxyyy_xxxyy_0[i] * fbe_0 - 3.0 * g_xxyyy_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxxyyy_xxyy_1[i] * fe_0 + g_xxxyyy_xxxyy_1[i] * pa_x[i];

        g_xxxxyyy_xxxyz_0[i] = 3.0 * g_xxyyy_xxxyz_0[i] * fbe_0 - 3.0 * g_xxyyy_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_xxyz_1[i] * fe_0 + g_xxxyyy_xxxyz_1[i] * pa_x[i];

        g_xxxxyyy_xxxzz_0[i] = 2.0 * g_xxxxy_xxxzz_0[i] * fbe_0 - 2.0 * g_xxxxy_xxxzz_1[i] * fz_be_0 + g_xxxxyy_xxxzz_1[i] * pa_y[i];

        g_xxxxyyy_xxyyy_0[i] = 3.0 * g_xxyyy_xxyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxxyyy_xyyy_1[i] * fe_0 + g_xxxyyy_xxyyy_1[i] * pa_x[i];

        g_xxxxyyy_xxyyz_0[i] = 3.0 * g_xxyyy_xxyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_xyyz_1[i] * fe_0 + g_xxxyyy_xxyyz_1[i] * pa_x[i];

        g_xxxxyyy_xxyzz_0[i] = 3.0 * g_xxyyy_xxyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_xyzz_1[i] * fe_0 + g_xxxyyy_xxyzz_1[i] * pa_x[i];

        g_xxxxyyy_xxzzz_0[i] = 2.0 * g_xxxxy_xxzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_xxzzz_1[i] * fz_be_0 + g_xxxxyy_xxzzz_1[i] * pa_y[i];

        g_xxxxyyy_xyyyy_0[i] = 3.0 * g_xxyyy_xyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_xyyyy_1[i] * fz_be_0 + g_xxxyyy_yyyy_1[i] * fe_0 + g_xxxyyy_xyyyy_1[i] * pa_x[i];

        g_xxxxyyy_xyyyz_0[i] = 3.0 * g_xxyyy_xyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_xyyyz_1[i] * fz_be_0 + g_xxxyyy_yyyz_1[i] * fe_0 + g_xxxyyy_xyyyz_1[i] * pa_x[i];

        g_xxxxyyy_xyyzz_0[i] = 3.0 * g_xxyyy_xyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_xyyzz_1[i] * fz_be_0 + g_xxxyyy_yyzz_1[i] * fe_0 + g_xxxyyy_xyyzz_1[i] * pa_x[i];

        g_xxxxyyy_xyzzz_0[i] = 3.0 * g_xxyyy_xyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_xyzzz_1[i] * fz_be_0 + g_xxxyyy_yzzz_1[i] * fe_0 + g_xxxyyy_xyzzz_1[i] * pa_x[i];

        g_xxxxyyy_xzzzz_0[i] = 2.0 * g_xxxxy_xzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_xzzzz_1[i] * fz_be_0 + g_xxxxyy_xzzzz_1[i] * pa_y[i];

        g_xxxxyyy_yyyyy_0[i] = 3.0 * g_xxyyy_yyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_yyyyy_1[i] * fz_be_0 + g_xxxyyy_yyyyy_1[i] * pa_x[i];

        g_xxxxyyy_yyyyz_0[i] = 3.0 * g_xxyyy_yyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_yyyyz_1[i] * fz_be_0 + g_xxxyyy_yyyyz_1[i] * pa_x[i];

        g_xxxxyyy_yyyzz_0[i] = 3.0 * g_xxyyy_yyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_yyyzz_1[i] * fz_be_0 + g_xxxyyy_yyyzz_1[i] * pa_x[i];

        g_xxxxyyy_yyzzz_0[i] = 3.0 * g_xxyyy_yyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_yyzzz_1[i] * fz_be_0 + g_xxxyyy_yyzzz_1[i] * pa_x[i];

        g_xxxxyyy_yzzzz_0[i] = 3.0 * g_xxyyy_yzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_yzzzz_1[i] * fz_be_0 + g_xxxyyy_yzzzz_1[i] * pa_x[i];

        g_xxxxyyy_zzzzz_0[i] = 3.0 * g_xxyyy_zzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_zzzzz_1[i] * fz_be_0 + g_xxxyyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 147-168 components of targeted buffer : KH

    auto g_xxxxyyz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 147);

    auto g_xxxxyyz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 148);

    auto g_xxxxyyz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 149);

    auto g_xxxxyyz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 150);

    auto g_xxxxyyz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 151);

    auto g_xxxxyyz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 152);

    auto g_xxxxyyz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 153);

    auto g_xxxxyyz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 154);

    auto g_xxxxyyz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 155);

    auto g_xxxxyyz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 156);

    auto g_xxxxyyz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 157);

    auto g_xxxxyyz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 158);

    auto g_xxxxyyz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 159);

    auto g_xxxxyyz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 160);

    auto g_xxxxyyz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 161);

    auto g_xxxxyyz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 162);

    auto g_xxxxyyz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 163);

    auto g_xxxxyyz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 164);

    auto g_xxxxyyz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 165);

    auto g_xxxxyyz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 166);

    auto g_xxxxyyz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 167);

    #pragma omp simd aligned(g_xxxxyy_xxxx_1, g_xxxxyy_xxxxx_1, g_xxxxyy_xxxxy_1, g_xxxxyy_xxxxz_1, g_xxxxyy_xxxy_1, g_xxxxyy_xxxyy_1, g_xxxxyy_xxxyz_1, g_xxxxyy_xxxz_1, g_xxxxyy_xxxzz_1, g_xxxxyy_xxyy_1, g_xxxxyy_xxyyy_1, g_xxxxyy_xxyyz_1, g_xxxxyy_xxyz_1, g_xxxxyy_xxyzz_1, g_xxxxyy_xxzz_1, g_xxxxyy_xxzzz_1, g_xxxxyy_xyyy_1, g_xxxxyy_xyyyy_1, g_xxxxyy_xyyyz_1, g_xxxxyy_xyyz_1, g_xxxxyy_xyyzz_1, g_xxxxyy_xyzz_1, g_xxxxyy_xyzzz_1, g_xxxxyy_xzzz_1, g_xxxxyy_xzzzz_1, g_xxxxyy_yyyy_1, g_xxxxyy_yyyyy_1, g_xxxxyy_yyyyz_1, g_xxxxyy_yyyz_1, g_xxxxyy_yyyzz_1, g_xxxxyy_yyzz_1, g_xxxxyy_yyzzz_1, g_xxxxyy_yzzz_1, g_xxxxyy_yzzzz_1, g_xxxxyy_zzzz_1, g_xxxxyy_zzzzz_1, g_xxxxyyz_xxxxx_0, g_xxxxyyz_xxxxy_0, g_xxxxyyz_xxxxz_0, g_xxxxyyz_xxxyy_0, g_xxxxyyz_xxxyz_0, g_xxxxyyz_xxxzz_0, g_xxxxyyz_xxyyy_0, g_xxxxyyz_xxyyz_0, g_xxxxyyz_xxyzz_0, g_xxxxyyz_xxzzz_0, g_xxxxyyz_xyyyy_0, g_xxxxyyz_xyyyz_0, g_xxxxyyz_xyyzz_0, g_xxxxyyz_xyzzz_0, g_xxxxyyz_xzzzz_0, g_xxxxyyz_yyyyy_0, g_xxxxyyz_yyyyz_0, g_xxxxyyz_yyyzz_0, g_xxxxyyz_yyzzz_0, g_xxxxyyz_yzzzz_0, g_xxxxyyz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyyz_xxxxx_0[i] = g_xxxxyy_xxxxx_1[i] * pa_z[i];

        g_xxxxyyz_xxxxy_0[i] = g_xxxxyy_xxxxy_1[i] * pa_z[i];

        g_xxxxyyz_xxxxz_0[i] = g_xxxxyy_xxxx_1[i] * fe_0 + g_xxxxyy_xxxxz_1[i] * pa_z[i];

        g_xxxxyyz_xxxyy_0[i] = g_xxxxyy_xxxyy_1[i] * pa_z[i];

        g_xxxxyyz_xxxyz_0[i] = g_xxxxyy_xxxy_1[i] * fe_0 + g_xxxxyy_xxxyz_1[i] * pa_z[i];

        g_xxxxyyz_xxxzz_0[i] = 2.0 * g_xxxxyy_xxxz_1[i] * fe_0 + g_xxxxyy_xxxzz_1[i] * pa_z[i];

        g_xxxxyyz_xxyyy_0[i] = g_xxxxyy_xxyyy_1[i] * pa_z[i];

        g_xxxxyyz_xxyyz_0[i] = g_xxxxyy_xxyy_1[i] * fe_0 + g_xxxxyy_xxyyz_1[i] * pa_z[i];

        g_xxxxyyz_xxyzz_0[i] = 2.0 * g_xxxxyy_xxyz_1[i] * fe_0 + g_xxxxyy_xxyzz_1[i] * pa_z[i];

        g_xxxxyyz_xxzzz_0[i] = 3.0 * g_xxxxyy_xxzz_1[i] * fe_0 + g_xxxxyy_xxzzz_1[i] * pa_z[i];

        g_xxxxyyz_xyyyy_0[i] = g_xxxxyy_xyyyy_1[i] * pa_z[i];

        g_xxxxyyz_xyyyz_0[i] = g_xxxxyy_xyyy_1[i] * fe_0 + g_xxxxyy_xyyyz_1[i] * pa_z[i];

        g_xxxxyyz_xyyzz_0[i] = 2.0 * g_xxxxyy_xyyz_1[i] * fe_0 + g_xxxxyy_xyyzz_1[i] * pa_z[i];

        g_xxxxyyz_xyzzz_0[i] = 3.0 * g_xxxxyy_xyzz_1[i] * fe_0 + g_xxxxyy_xyzzz_1[i] * pa_z[i];

        g_xxxxyyz_xzzzz_0[i] = 4.0 * g_xxxxyy_xzzz_1[i] * fe_0 + g_xxxxyy_xzzzz_1[i] * pa_z[i];

        g_xxxxyyz_yyyyy_0[i] = g_xxxxyy_yyyyy_1[i] * pa_z[i];

        g_xxxxyyz_yyyyz_0[i] = g_xxxxyy_yyyy_1[i] * fe_0 + g_xxxxyy_yyyyz_1[i] * pa_z[i];

        g_xxxxyyz_yyyzz_0[i] = 2.0 * g_xxxxyy_yyyz_1[i] * fe_0 + g_xxxxyy_yyyzz_1[i] * pa_z[i];

        g_xxxxyyz_yyzzz_0[i] = 3.0 * g_xxxxyy_yyzz_1[i] * fe_0 + g_xxxxyy_yyzzz_1[i] * pa_z[i];

        g_xxxxyyz_yzzzz_0[i] = 4.0 * g_xxxxyy_yzzz_1[i] * fe_0 + g_xxxxyy_yzzzz_1[i] * pa_z[i];

        g_xxxxyyz_zzzzz_0[i] = 5.0 * g_xxxxyy_zzzz_1[i] * fe_0 + g_xxxxyy_zzzzz_1[i] * pa_z[i];
    }

    // Set up 168-189 components of targeted buffer : KH

    auto g_xxxxyzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 168);

    auto g_xxxxyzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 169);

    auto g_xxxxyzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 170);

    auto g_xxxxyzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 171);

    auto g_xxxxyzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 172);

    auto g_xxxxyzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 173);

    auto g_xxxxyzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 174);

    auto g_xxxxyzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 175);

    auto g_xxxxyzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 176);

    auto g_xxxxyzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 177);

    auto g_xxxxyzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 178);

    auto g_xxxxyzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 179);

    auto g_xxxxyzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 180);

    auto g_xxxxyzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 181);

    auto g_xxxxyzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 182);

    auto g_xxxxyzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 183);

    auto g_xxxxyzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 184);

    auto g_xxxxyzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 185);

    auto g_xxxxyzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 186);

    auto g_xxxxyzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 187);

    auto g_xxxxyzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 188);

    #pragma omp simd aligned(g_xxxxyzz_xxxxx_0, g_xxxxyzz_xxxxy_0, g_xxxxyzz_xxxxz_0, g_xxxxyzz_xxxyy_0, g_xxxxyzz_xxxyz_0, g_xxxxyzz_xxxzz_0, g_xxxxyzz_xxyyy_0, g_xxxxyzz_xxyyz_0, g_xxxxyzz_xxyzz_0, g_xxxxyzz_xxzzz_0, g_xxxxyzz_xyyyy_0, g_xxxxyzz_xyyyz_0, g_xxxxyzz_xyyzz_0, g_xxxxyzz_xyzzz_0, g_xxxxyzz_xzzzz_0, g_xxxxyzz_yyyyy_0, g_xxxxyzz_yyyyz_0, g_xxxxyzz_yyyzz_0, g_xxxxyzz_yyzzz_0, g_xxxxyzz_yzzzz_0, g_xxxxyzz_zzzzz_0, g_xxxxzz_xxxx_1, g_xxxxzz_xxxxx_1, g_xxxxzz_xxxxy_1, g_xxxxzz_xxxxz_1, g_xxxxzz_xxxy_1, g_xxxxzz_xxxyy_1, g_xxxxzz_xxxyz_1, g_xxxxzz_xxxz_1, g_xxxxzz_xxxzz_1, g_xxxxzz_xxyy_1, g_xxxxzz_xxyyy_1, g_xxxxzz_xxyyz_1, g_xxxxzz_xxyz_1, g_xxxxzz_xxyzz_1, g_xxxxzz_xxzz_1, g_xxxxzz_xxzzz_1, g_xxxxzz_xyyy_1, g_xxxxzz_xyyyy_1, g_xxxxzz_xyyyz_1, g_xxxxzz_xyyz_1, g_xxxxzz_xyyzz_1, g_xxxxzz_xyzz_1, g_xxxxzz_xyzzz_1, g_xxxxzz_xzzz_1, g_xxxxzz_xzzzz_1, g_xxxxzz_yyyy_1, g_xxxxzz_yyyyy_1, g_xxxxzz_yyyyz_1, g_xxxxzz_yyyz_1, g_xxxxzz_yyyzz_1, g_xxxxzz_yyzz_1, g_xxxxzz_yyzzz_1, g_xxxxzz_yzzz_1, g_xxxxzz_yzzzz_1, g_xxxxzz_zzzz_1, g_xxxxzz_zzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyzz_xxxxx_0[i] = g_xxxxzz_xxxxx_1[i] * pa_y[i];

        g_xxxxyzz_xxxxy_0[i] = g_xxxxzz_xxxx_1[i] * fe_0 + g_xxxxzz_xxxxy_1[i] * pa_y[i];

        g_xxxxyzz_xxxxz_0[i] = g_xxxxzz_xxxxz_1[i] * pa_y[i];

        g_xxxxyzz_xxxyy_0[i] = 2.0 * g_xxxxzz_xxxy_1[i] * fe_0 + g_xxxxzz_xxxyy_1[i] * pa_y[i];

        g_xxxxyzz_xxxyz_0[i] = g_xxxxzz_xxxz_1[i] * fe_0 + g_xxxxzz_xxxyz_1[i] * pa_y[i];

        g_xxxxyzz_xxxzz_0[i] = g_xxxxzz_xxxzz_1[i] * pa_y[i];

        g_xxxxyzz_xxyyy_0[i] = 3.0 * g_xxxxzz_xxyy_1[i] * fe_0 + g_xxxxzz_xxyyy_1[i] * pa_y[i];

        g_xxxxyzz_xxyyz_0[i] = 2.0 * g_xxxxzz_xxyz_1[i] * fe_0 + g_xxxxzz_xxyyz_1[i] * pa_y[i];

        g_xxxxyzz_xxyzz_0[i] = g_xxxxzz_xxzz_1[i] * fe_0 + g_xxxxzz_xxyzz_1[i] * pa_y[i];

        g_xxxxyzz_xxzzz_0[i] = g_xxxxzz_xxzzz_1[i] * pa_y[i];

        g_xxxxyzz_xyyyy_0[i] = 4.0 * g_xxxxzz_xyyy_1[i] * fe_0 + g_xxxxzz_xyyyy_1[i] * pa_y[i];

        g_xxxxyzz_xyyyz_0[i] = 3.0 * g_xxxxzz_xyyz_1[i] * fe_0 + g_xxxxzz_xyyyz_1[i] * pa_y[i];

        g_xxxxyzz_xyyzz_0[i] = 2.0 * g_xxxxzz_xyzz_1[i] * fe_0 + g_xxxxzz_xyyzz_1[i] * pa_y[i];

        g_xxxxyzz_xyzzz_0[i] = g_xxxxzz_xzzz_1[i] * fe_0 + g_xxxxzz_xyzzz_1[i] * pa_y[i];

        g_xxxxyzz_xzzzz_0[i] = g_xxxxzz_xzzzz_1[i] * pa_y[i];

        g_xxxxyzz_yyyyy_0[i] = 5.0 * g_xxxxzz_yyyy_1[i] * fe_0 + g_xxxxzz_yyyyy_1[i] * pa_y[i];

        g_xxxxyzz_yyyyz_0[i] = 4.0 * g_xxxxzz_yyyz_1[i] * fe_0 + g_xxxxzz_yyyyz_1[i] * pa_y[i];

        g_xxxxyzz_yyyzz_0[i] = 3.0 * g_xxxxzz_yyzz_1[i] * fe_0 + g_xxxxzz_yyyzz_1[i] * pa_y[i];

        g_xxxxyzz_yyzzz_0[i] = 2.0 * g_xxxxzz_yzzz_1[i] * fe_0 + g_xxxxzz_yyzzz_1[i] * pa_y[i];

        g_xxxxyzz_yzzzz_0[i] = g_xxxxzz_zzzz_1[i] * fe_0 + g_xxxxzz_yzzzz_1[i] * pa_y[i];

        g_xxxxyzz_zzzzz_0[i] = g_xxxxzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 189-210 components of targeted buffer : KH

    auto g_xxxxzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 189);

    auto g_xxxxzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 190);

    auto g_xxxxzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 191);

    auto g_xxxxzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 192);

    auto g_xxxxzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 193);

    auto g_xxxxzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 194);

    auto g_xxxxzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 195);

    auto g_xxxxzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 196);

    auto g_xxxxzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 197);

    auto g_xxxxzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 198);

    auto g_xxxxzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 199);

    auto g_xxxxzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 200);

    auto g_xxxxzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 201);

    auto g_xxxxzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 202);

    auto g_xxxxzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 203);

    auto g_xxxxzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 204);

    auto g_xxxxzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 205);

    auto g_xxxxzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 206);

    auto g_xxxxzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 207);

    auto g_xxxxzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 208);

    auto g_xxxxzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 209);

    #pragma omp simd aligned(g_xxxxz_xxxxx_0, g_xxxxz_xxxxx_1, g_xxxxz_xxxxy_0, g_xxxxz_xxxxy_1, g_xxxxz_xxxyy_0, g_xxxxz_xxxyy_1, g_xxxxz_xxyyy_0, g_xxxxz_xxyyy_1, g_xxxxz_xyyyy_0, g_xxxxz_xyyyy_1, g_xxxxzz_xxxxx_1, g_xxxxzz_xxxxy_1, g_xxxxzz_xxxyy_1, g_xxxxzz_xxyyy_1, g_xxxxzz_xyyyy_1, g_xxxxzzz_xxxxx_0, g_xxxxzzz_xxxxy_0, g_xxxxzzz_xxxxz_0, g_xxxxzzz_xxxyy_0, g_xxxxzzz_xxxyz_0, g_xxxxzzz_xxxzz_0, g_xxxxzzz_xxyyy_0, g_xxxxzzz_xxyyz_0, g_xxxxzzz_xxyzz_0, g_xxxxzzz_xxzzz_0, g_xxxxzzz_xyyyy_0, g_xxxxzzz_xyyyz_0, g_xxxxzzz_xyyzz_0, g_xxxxzzz_xyzzz_0, g_xxxxzzz_xzzzz_0, g_xxxxzzz_yyyyy_0, g_xxxxzzz_yyyyz_0, g_xxxxzzz_yyyzz_0, g_xxxxzzz_yyzzz_0, g_xxxxzzz_yzzzz_0, g_xxxxzzz_zzzzz_0, g_xxxzzz_xxxxz_1, g_xxxzzz_xxxyz_1, g_xxxzzz_xxxz_1, g_xxxzzz_xxxzz_1, g_xxxzzz_xxyyz_1, g_xxxzzz_xxyz_1, g_xxxzzz_xxyzz_1, g_xxxzzz_xxzz_1, g_xxxzzz_xxzzz_1, g_xxxzzz_xyyyz_1, g_xxxzzz_xyyz_1, g_xxxzzz_xyyzz_1, g_xxxzzz_xyzz_1, g_xxxzzz_xyzzz_1, g_xxxzzz_xzzz_1, g_xxxzzz_xzzzz_1, g_xxxzzz_yyyyy_1, g_xxxzzz_yyyyz_1, g_xxxzzz_yyyz_1, g_xxxzzz_yyyzz_1, g_xxxzzz_yyzz_1, g_xxxzzz_yyzzz_1, g_xxxzzz_yzzz_1, g_xxxzzz_yzzzz_1, g_xxxzzz_zzzz_1, g_xxxzzz_zzzzz_1, g_xxzzz_xxxxz_0, g_xxzzz_xxxxz_1, g_xxzzz_xxxyz_0, g_xxzzz_xxxyz_1, g_xxzzz_xxxzz_0, g_xxzzz_xxxzz_1, g_xxzzz_xxyyz_0, g_xxzzz_xxyyz_1, g_xxzzz_xxyzz_0, g_xxzzz_xxyzz_1, g_xxzzz_xxzzz_0, g_xxzzz_xxzzz_1, g_xxzzz_xyyyz_0, g_xxzzz_xyyyz_1, g_xxzzz_xyyzz_0, g_xxzzz_xyyzz_1, g_xxzzz_xyzzz_0, g_xxzzz_xyzzz_1, g_xxzzz_xzzzz_0, g_xxzzz_xzzzz_1, g_xxzzz_yyyyy_0, g_xxzzz_yyyyy_1, g_xxzzz_yyyyz_0, g_xxzzz_yyyyz_1, g_xxzzz_yyyzz_0, g_xxzzz_yyyzz_1, g_xxzzz_yyzzz_0, g_xxzzz_yyzzz_1, g_xxzzz_yzzzz_0, g_xxzzz_yzzzz_1, g_xxzzz_zzzzz_0, g_xxzzz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxzzz_xxxxx_0[i] = 2.0 * g_xxxxz_xxxxx_0[i] * fbe_0 - 2.0 * g_xxxxz_xxxxx_1[i] * fz_be_0 + g_xxxxzz_xxxxx_1[i] * pa_z[i];

        g_xxxxzzz_xxxxy_0[i] = 2.0 * g_xxxxz_xxxxy_0[i] * fbe_0 - 2.0 * g_xxxxz_xxxxy_1[i] * fz_be_0 + g_xxxxzz_xxxxy_1[i] * pa_z[i];

        g_xxxxzzz_xxxxz_0[i] = 3.0 * g_xxzzz_xxxxz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_xxxz_1[i] * fe_0 + g_xxxzzz_xxxxz_1[i] * pa_x[i];

        g_xxxxzzz_xxxyy_0[i] = 2.0 * g_xxxxz_xxxyy_0[i] * fbe_0 - 2.0 * g_xxxxz_xxxyy_1[i] * fz_be_0 + g_xxxxzz_xxxyy_1[i] * pa_z[i];

        g_xxxxzzz_xxxyz_0[i] = 3.0 * g_xxzzz_xxxyz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_xxyz_1[i] * fe_0 + g_xxxzzz_xxxyz_1[i] * pa_x[i];

        g_xxxxzzz_xxxzz_0[i] = 3.0 * g_xxzzz_xxxzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_xxzz_1[i] * fe_0 + g_xxxzzz_xxxzz_1[i] * pa_x[i];

        g_xxxxzzz_xxyyy_0[i] = 2.0 * g_xxxxz_xxyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_xxyyy_1[i] * fz_be_0 + g_xxxxzz_xxyyy_1[i] * pa_z[i];

        g_xxxxzzz_xxyyz_0[i] = 3.0 * g_xxzzz_xxyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_xyyz_1[i] * fe_0 + g_xxxzzz_xxyyz_1[i] * pa_x[i];

        g_xxxxzzz_xxyzz_0[i] = 3.0 * g_xxzzz_xxyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_xyzz_1[i] * fe_0 + g_xxxzzz_xxyzz_1[i] * pa_x[i];

        g_xxxxzzz_xxzzz_0[i] = 3.0 * g_xxzzz_xxzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_xzzz_1[i] * fe_0 + g_xxxzzz_xxzzz_1[i] * pa_x[i];

        g_xxxxzzz_xyyyy_0[i] = 2.0 * g_xxxxz_xyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_xyyyy_1[i] * fz_be_0 + g_xxxxzz_xyyyy_1[i] * pa_z[i];

        g_xxxxzzz_xyyyz_0[i] = 3.0 * g_xxzzz_xyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_xyyyz_1[i] * fz_be_0 + g_xxxzzz_yyyz_1[i] * fe_0 + g_xxxzzz_xyyyz_1[i] * pa_x[i];

        g_xxxxzzz_xyyzz_0[i] = 3.0 * g_xxzzz_xyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xyyzz_1[i] * fz_be_0 + g_xxxzzz_yyzz_1[i] * fe_0 + g_xxxzzz_xyyzz_1[i] * pa_x[i];

        g_xxxxzzz_xyzzz_0[i] = 3.0 * g_xxzzz_xyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xyzzz_1[i] * fz_be_0 + g_xxxzzz_yzzz_1[i] * fe_0 + g_xxxzzz_xyzzz_1[i] * pa_x[i];

        g_xxxxzzz_xzzzz_0[i] = 3.0 * g_xxzzz_xzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xzzzz_1[i] * fz_be_0 + g_xxxzzz_zzzz_1[i] * fe_0 + g_xxxzzz_xzzzz_1[i] * pa_x[i];

        g_xxxxzzz_yyyyy_0[i] = 3.0 * g_xxzzz_yyyyy_0[i] * fbe_0 - 3.0 * g_xxzzz_yyyyy_1[i] * fz_be_0 + g_xxxzzz_yyyyy_1[i] * pa_x[i];

        g_xxxxzzz_yyyyz_0[i] = 3.0 * g_xxzzz_yyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_yyyyz_1[i] * fz_be_0 + g_xxxzzz_yyyyz_1[i] * pa_x[i];

        g_xxxxzzz_yyyzz_0[i] = 3.0 * g_xxzzz_yyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_yyyzz_1[i] * fz_be_0 + g_xxxzzz_yyyzz_1[i] * pa_x[i];

        g_xxxxzzz_yyzzz_0[i] = 3.0 * g_xxzzz_yyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_yyzzz_1[i] * fz_be_0 + g_xxxzzz_yyzzz_1[i] * pa_x[i];

        g_xxxxzzz_yzzzz_0[i] = 3.0 * g_xxzzz_yzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_yzzzz_1[i] * fz_be_0 + g_xxxzzz_yzzzz_1[i] * pa_x[i];

        g_xxxxzzz_zzzzz_0[i] = 3.0 * g_xxzzz_zzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_zzzzz_1[i] * fz_be_0 + g_xxxzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 210-231 components of targeted buffer : KH

    auto g_xxxyyyy_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 210);

    auto g_xxxyyyy_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 211);

    auto g_xxxyyyy_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 212);

    auto g_xxxyyyy_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 213);

    auto g_xxxyyyy_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 214);

    auto g_xxxyyyy_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 215);

    auto g_xxxyyyy_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 216);

    auto g_xxxyyyy_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 217);

    auto g_xxxyyyy_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 218);

    auto g_xxxyyyy_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 219);

    auto g_xxxyyyy_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 220);

    auto g_xxxyyyy_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 221);

    auto g_xxxyyyy_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 222);

    auto g_xxxyyyy_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 223);

    auto g_xxxyyyy_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 224);

    auto g_xxxyyyy_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 225);

    auto g_xxxyyyy_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 226);

    auto g_xxxyyyy_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 227);

    auto g_xxxyyyy_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 228);

    auto g_xxxyyyy_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 229);

    auto g_xxxyyyy_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 230);

    #pragma omp simd aligned(g_xxxyy_xxxxx_0, g_xxxyy_xxxxx_1, g_xxxyy_xxxxz_0, g_xxxyy_xxxxz_1, g_xxxyy_xxxzz_0, g_xxxyy_xxxzz_1, g_xxxyy_xxzzz_0, g_xxxyy_xxzzz_1, g_xxxyy_xzzzz_0, g_xxxyy_xzzzz_1, g_xxxyyy_xxxxx_1, g_xxxyyy_xxxxz_1, g_xxxyyy_xxxzz_1, g_xxxyyy_xxzzz_1, g_xxxyyy_xzzzz_1, g_xxxyyyy_xxxxx_0, g_xxxyyyy_xxxxy_0, g_xxxyyyy_xxxxz_0, g_xxxyyyy_xxxyy_0, g_xxxyyyy_xxxyz_0, g_xxxyyyy_xxxzz_0, g_xxxyyyy_xxyyy_0, g_xxxyyyy_xxyyz_0, g_xxxyyyy_xxyzz_0, g_xxxyyyy_xxzzz_0, g_xxxyyyy_xyyyy_0, g_xxxyyyy_xyyyz_0, g_xxxyyyy_xyyzz_0, g_xxxyyyy_xyzzz_0, g_xxxyyyy_xzzzz_0, g_xxxyyyy_yyyyy_0, g_xxxyyyy_yyyyz_0, g_xxxyyyy_yyyzz_0, g_xxxyyyy_yyzzz_0, g_xxxyyyy_yzzzz_0, g_xxxyyyy_zzzzz_0, g_xxyyyy_xxxxy_1, g_xxyyyy_xxxy_1, g_xxyyyy_xxxyy_1, g_xxyyyy_xxxyz_1, g_xxyyyy_xxyy_1, g_xxyyyy_xxyyy_1, g_xxyyyy_xxyyz_1, g_xxyyyy_xxyz_1, g_xxyyyy_xxyzz_1, g_xxyyyy_xyyy_1, g_xxyyyy_xyyyy_1, g_xxyyyy_xyyyz_1, g_xxyyyy_xyyz_1, g_xxyyyy_xyyzz_1, g_xxyyyy_xyzz_1, g_xxyyyy_xyzzz_1, g_xxyyyy_yyyy_1, g_xxyyyy_yyyyy_1, g_xxyyyy_yyyyz_1, g_xxyyyy_yyyz_1, g_xxyyyy_yyyzz_1, g_xxyyyy_yyzz_1, g_xxyyyy_yyzzz_1, g_xxyyyy_yzzz_1, g_xxyyyy_yzzzz_1, g_xxyyyy_zzzzz_1, g_xyyyy_xxxxy_0, g_xyyyy_xxxxy_1, g_xyyyy_xxxyy_0, g_xyyyy_xxxyy_1, g_xyyyy_xxxyz_0, g_xyyyy_xxxyz_1, g_xyyyy_xxyyy_0, g_xyyyy_xxyyy_1, g_xyyyy_xxyyz_0, g_xyyyy_xxyyz_1, g_xyyyy_xxyzz_0, g_xyyyy_xxyzz_1, g_xyyyy_xyyyy_0, g_xyyyy_xyyyy_1, g_xyyyy_xyyyz_0, g_xyyyy_xyyyz_1, g_xyyyy_xyyzz_0, g_xyyyy_xyyzz_1, g_xyyyy_xyzzz_0, g_xyyyy_xyzzz_1, g_xyyyy_yyyyy_0, g_xyyyy_yyyyy_1, g_xyyyy_yyyyz_0, g_xyyyy_yyyyz_1, g_xyyyy_yyyzz_0, g_xyyyy_yyyzz_1, g_xyyyy_yyzzz_0, g_xyyyy_yyzzz_1, g_xyyyy_yzzzz_0, g_xyyyy_yzzzz_1, g_xyyyy_zzzzz_0, g_xyyyy_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyyy_xxxxx_0[i] = 3.0 * g_xxxyy_xxxxx_0[i] * fbe_0 - 3.0 * g_xxxyy_xxxxx_1[i] * fz_be_0 + g_xxxyyy_xxxxx_1[i] * pa_y[i];

        g_xxxyyyy_xxxxy_0[i] = 2.0 * g_xyyyy_xxxxy_0[i] * fbe_0 - 2.0 * g_xyyyy_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxyyyy_xxxy_1[i] * fe_0 + g_xxyyyy_xxxxy_1[i] * pa_x[i];

        g_xxxyyyy_xxxxz_0[i] = 3.0 * g_xxxyy_xxxxz_0[i] * fbe_0 - 3.0 * g_xxxyy_xxxxz_1[i] * fz_be_0 + g_xxxyyy_xxxxz_1[i] * pa_y[i];

        g_xxxyyyy_xxxyy_0[i] = 2.0 * g_xyyyy_xxxyy_0[i] * fbe_0 - 2.0 * g_xyyyy_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxyyyy_xxyy_1[i] * fe_0 + g_xxyyyy_xxxyy_1[i] * pa_x[i];

        g_xxxyyyy_xxxyz_0[i] = 2.0 * g_xyyyy_xxxyz_0[i] * fbe_0 - 2.0 * g_xyyyy_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_xxyz_1[i] * fe_0 + g_xxyyyy_xxxyz_1[i] * pa_x[i];

        g_xxxyyyy_xxxzz_0[i] = 3.0 * g_xxxyy_xxxzz_0[i] * fbe_0 - 3.0 * g_xxxyy_xxxzz_1[i] * fz_be_0 + g_xxxyyy_xxxzz_1[i] * pa_y[i];

        g_xxxyyyy_xxyyy_0[i] = 2.0 * g_xyyyy_xxyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxyyyy_xyyy_1[i] * fe_0 + g_xxyyyy_xxyyy_1[i] * pa_x[i];

        g_xxxyyyy_xxyyz_0[i] = 2.0 * g_xyyyy_xxyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_xyyz_1[i] * fe_0 + g_xxyyyy_xxyyz_1[i] * pa_x[i];

        g_xxxyyyy_xxyzz_0[i] = 2.0 * g_xyyyy_xxyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_xyzz_1[i] * fe_0 + g_xxyyyy_xxyzz_1[i] * pa_x[i];

        g_xxxyyyy_xxzzz_0[i] = 3.0 * g_xxxyy_xxzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_xxzzz_1[i] * fz_be_0 + g_xxxyyy_xxzzz_1[i] * pa_y[i];

        g_xxxyyyy_xyyyy_0[i] = 2.0 * g_xyyyy_xyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_xyyyy_1[i] * fz_be_0 + g_xxyyyy_yyyy_1[i] * fe_0 + g_xxyyyy_xyyyy_1[i] * pa_x[i];

        g_xxxyyyy_xyyyz_0[i] = 2.0 * g_xyyyy_xyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_xyyyz_1[i] * fz_be_0 + g_xxyyyy_yyyz_1[i] * fe_0 + g_xxyyyy_xyyyz_1[i] * pa_x[i];

        g_xxxyyyy_xyyzz_0[i] = 2.0 * g_xyyyy_xyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_xyyzz_1[i] * fz_be_0 + g_xxyyyy_yyzz_1[i] * fe_0 + g_xxyyyy_xyyzz_1[i] * pa_x[i];

        g_xxxyyyy_xyzzz_0[i] = 2.0 * g_xyyyy_xyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_xyzzz_1[i] * fz_be_0 + g_xxyyyy_yzzz_1[i] * fe_0 + g_xxyyyy_xyzzz_1[i] * pa_x[i];

        g_xxxyyyy_xzzzz_0[i] = 3.0 * g_xxxyy_xzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_xzzzz_1[i] * fz_be_0 + g_xxxyyy_xzzzz_1[i] * pa_y[i];

        g_xxxyyyy_yyyyy_0[i] = 2.0 * g_xyyyy_yyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_yyyyy_1[i] * fz_be_0 + g_xxyyyy_yyyyy_1[i] * pa_x[i];

        g_xxxyyyy_yyyyz_0[i] = 2.0 * g_xyyyy_yyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_yyyyz_1[i] * fz_be_0 + g_xxyyyy_yyyyz_1[i] * pa_x[i];

        g_xxxyyyy_yyyzz_0[i] = 2.0 * g_xyyyy_yyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_yyyzz_1[i] * fz_be_0 + g_xxyyyy_yyyzz_1[i] * pa_x[i];

        g_xxxyyyy_yyzzz_0[i] = 2.0 * g_xyyyy_yyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_yyzzz_1[i] * fz_be_0 + g_xxyyyy_yyzzz_1[i] * pa_x[i];

        g_xxxyyyy_yzzzz_0[i] = 2.0 * g_xyyyy_yzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_yzzzz_1[i] * fz_be_0 + g_xxyyyy_yzzzz_1[i] * pa_x[i];

        g_xxxyyyy_zzzzz_0[i] = 2.0 * g_xyyyy_zzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_zzzzz_1[i] * fz_be_0 + g_xxyyyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 231-252 components of targeted buffer : KH

    auto g_xxxyyyz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 231);

    auto g_xxxyyyz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 232);

    auto g_xxxyyyz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 233);

    auto g_xxxyyyz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 234);

    auto g_xxxyyyz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 235);

    auto g_xxxyyyz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 236);

    auto g_xxxyyyz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 237);

    auto g_xxxyyyz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 238);

    auto g_xxxyyyz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 239);

    auto g_xxxyyyz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 240);

    auto g_xxxyyyz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 241);

    auto g_xxxyyyz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 242);

    auto g_xxxyyyz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 243);

    auto g_xxxyyyz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 244);

    auto g_xxxyyyz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 245);

    auto g_xxxyyyz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 246);

    auto g_xxxyyyz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 247);

    auto g_xxxyyyz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 248);

    auto g_xxxyyyz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 249);

    auto g_xxxyyyz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 250);

    auto g_xxxyyyz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 251);

    #pragma omp simd aligned(g_xxxyyy_xxxx_1, g_xxxyyy_xxxxx_1, g_xxxyyy_xxxxy_1, g_xxxyyy_xxxxz_1, g_xxxyyy_xxxy_1, g_xxxyyy_xxxyy_1, g_xxxyyy_xxxyz_1, g_xxxyyy_xxxz_1, g_xxxyyy_xxxzz_1, g_xxxyyy_xxyy_1, g_xxxyyy_xxyyy_1, g_xxxyyy_xxyyz_1, g_xxxyyy_xxyz_1, g_xxxyyy_xxyzz_1, g_xxxyyy_xxzz_1, g_xxxyyy_xxzzz_1, g_xxxyyy_xyyy_1, g_xxxyyy_xyyyy_1, g_xxxyyy_xyyyz_1, g_xxxyyy_xyyz_1, g_xxxyyy_xyyzz_1, g_xxxyyy_xyzz_1, g_xxxyyy_xyzzz_1, g_xxxyyy_xzzz_1, g_xxxyyy_xzzzz_1, g_xxxyyy_yyyy_1, g_xxxyyy_yyyyy_1, g_xxxyyy_yyyyz_1, g_xxxyyy_yyyz_1, g_xxxyyy_yyyzz_1, g_xxxyyy_yyzz_1, g_xxxyyy_yyzzz_1, g_xxxyyy_yzzz_1, g_xxxyyy_yzzzz_1, g_xxxyyy_zzzz_1, g_xxxyyy_zzzzz_1, g_xxxyyyz_xxxxx_0, g_xxxyyyz_xxxxy_0, g_xxxyyyz_xxxxz_0, g_xxxyyyz_xxxyy_0, g_xxxyyyz_xxxyz_0, g_xxxyyyz_xxxzz_0, g_xxxyyyz_xxyyy_0, g_xxxyyyz_xxyyz_0, g_xxxyyyz_xxyzz_0, g_xxxyyyz_xxzzz_0, g_xxxyyyz_xyyyy_0, g_xxxyyyz_xyyyz_0, g_xxxyyyz_xyyzz_0, g_xxxyyyz_xyzzz_0, g_xxxyyyz_xzzzz_0, g_xxxyyyz_yyyyy_0, g_xxxyyyz_yyyyz_0, g_xxxyyyz_yyyzz_0, g_xxxyyyz_yyzzz_0, g_xxxyyyz_yzzzz_0, g_xxxyyyz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyyyz_xxxxx_0[i] = g_xxxyyy_xxxxx_1[i] * pa_z[i];

        g_xxxyyyz_xxxxy_0[i] = g_xxxyyy_xxxxy_1[i] * pa_z[i];

        g_xxxyyyz_xxxxz_0[i] = g_xxxyyy_xxxx_1[i] * fe_0 + g_xxxyyy_xxxxz_1[i] * pa_z[i];

        g_xxxyyyz_xxxyy_0[i] = g_xxxyyy_xxxyy_1[i] * pa_z[i];

        g_xxxyyyz_xxxyz_0[i] = g_xxxyyy_xxxy_1[i] * fe_0 + g_xxxyyy_xxxyz_1[i] * pa_z[i];

        g_xxxyyyz_xxxzz_0[i] = 2.0 * g_xxxyyy_xxxz_1[i] * fe_0 + g_xxxyyy_xxxzz_1[i] * pa_z[i];

        g_xxxyyyz_xxyyy_0[i] = g_xxxyyy_xxyyy_1[i] * pa_z[i];

        g_xxxyyyz_xxyyz_0[i] = g_xxxyyy_xxyy_1[i] * fe_0 + g_xxxyyy_xxyyz_1[i] * pa_z[i];

        g_xxxyyyz_xxyzz_0[i] = 2.0 * g_xxxyyy_xxyz_1[i] * fe_0 + g_xxxyyy_xxyzz_1[i] * pa_z[i];

        g_xxxyyyz_xxzzz_0[i] = 3.0 * g_xxxyyy_xxzz_1[i] * fe_0 + g_xxxyyy_xxzzz_1[i] * pa_z[i];

        g_xxxyyyz_xyyyy_0[i] = g_xxxyyy_xyyyy_1[i] * pa_z[i];

        g_xxxyyyz_xyyyz_0[i] = g_xxxyyy_xyyy_1[i] * fe_0 + g_xxxyyy_xyyyz_1[i] * pa_z[i];

        g_xxxyyyz_xyyzz_0[i] = 2.0 * g_xxxyyy_xyyz_1[i] * fe_0 + g_xxxyyy_xyyzz_1[i] * pa_z[i];

        g_xxxyyyz_xyzzz_0[i] = 3.0 * g_xxxyyy_xyzz_1[i] * fe_0 + g_xxxyyy_xyzzz_1[i] * pa_z[i];

        g_xxxyyyz_xzzzz_0[i] = 4.0 * g_xxxyyy_xzzz_1[i] * fe_0 + g_xxxyyy_xzzzz_1[i] * pa_z[i];

        g_xxxyyyz_yyyyy_0[i] = g_xxxyyy_yyyyy_1[i] * pa_z[i];

        g_xxxyyyz_yyyyz_0[i] = g_xxxyyy_yyyy_1[i] * fe_0 + g_xxxyyy_yyyyz_1[i] * pa_z[i];

        g_xxxyyyz_yyyzz_0[i] = 2.0 * g_xxxyyy_yyyz_1[i] * fe_0 + g_xxxyyy_yyyzz_1[i] * pa_z[i];

        g_xxxyyyz_yyzzz_0[i] = 3.0 * g_xxxyyy_yyzz_1[i] * fe_0 + g_xxxyyy_yyzzz_1[i] * pa_z[i];

        g_xxxyyyz_yzzzz_0[i] = 4.0 * g_xxxyyy_yzzz_1[i] * fe_0 + g_xxxyyy_yzzzz_1[i] * pa_z[i];

        g_xxxyyyz_zzzzz_0[i] = 5.0 * g_xxxyyy_zzzz_1[i] * fe_0 + g_xxxyyy_zzzzz_1[i] * pa_z[i];
    }

    // Set up 252-273 components of targeted buffer : KH

    auto g_xxxyyzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 252);

    auto g_xxxyyzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 253);

    auto g_xxxyyzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 254);

    auto g_xxxyyzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 255);

    auto g_xxxyyzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 256);

    auto g_xxxyyzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 257);

    auto g_xxxyyzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 258);

    auto g_xxxyyzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 259);

    auto g_xxxyyzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 260);

    auto g_xxxyyzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 261);

    auto g_xxxyyzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 262);

    auto g_xxxyyzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 263);

    auto g_xxxyyzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 264);

    auto g_xxxyyzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 265);

    auto g_xxxyyzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 266);

    auto g_xxxyyzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 267);

    auto g_xxxyyzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 268);

    auto g_xxxyyzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 269);

    auto g_xxxyyzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 270);

    auto g_xxxyyzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 271);

    auto g_xxxyyzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 272);

    #pragma omp simd aligned(g_xxxyy_xxxxy_0, g_xxxyy_xxxxy_1, g_xxxyy_xxxyy_0, g_xxxyy_xxxyy_1, g_xxxyy_xxyyy_0, g_xxxyy_xxyyy_1, g_xxxyy_xyyyy_0, g_xxxyy_xyyyy_1, g_xxxyyz_xxxxy_1, g_xxxyyz_xxxyy_1, g_xxxyyz_xxyyy_1, g_xxxyyz_xyyyy_1, g_xxxyyzz_xxxxx_0, g_xxxyyzz_xxxxy_0, g_xxxyyzz_xxxxz_0, g_xxxyyzz_xxxyy_0, g_xxxyyzz_xxxyz_0, g_xxxyyzz_xxxzz_0, g_xxxyyzz_xxyyy_0, g_xxxyyzz_xxyyz_0, g_xxxyyzz_xxyzz_0, g_xxxyyzz_xxzzz_0, g_xxxyyzz_xyyyy_0, g_xxxyyzz_xyyyz_0, g_xxxyyzz_xyyzz_0, g_xxxyyzz_xyzzz_0, g_xxxyyzz_xzzzz_0, g_xxxyyzz_yyyyy_0, g_xxxyyzz_yyyyz_0, g_xxxyyzz_yyyzz_0, g_xxxyyzz_yyzzz_0, g_xxxyyzz_yzzzz_0, g_xxxyyzz_zzzzz_0, g_xxxyzz_xxxxx_1, g_xxxyzz_xxxxz_1, g_xxxyzz_xxxzz_1, g_xxxyzz_xxzzz_1, g_xxxyzz_xzzzz_1, g_xxxzz_xxxxx_0, g_xxxzz_xxxxx_1, g_xxxzz_xxxxz_0, g_xxxzz_xxxxz_1, g_xxxzz_xxxzz_0, g_xxxzz_xxxzz_1, g_xxxzz_xxzzz_0, g_xxxzz_xxzzz_1, g_xxxzz_xzzzz_0, g_xxxzz_xzzzz_1, g_xxyyzz_xxxyz_1, g_xxyyzz_xxyyz_1, g_xxyyzz_xxyz_1, g_xxyyzz_xxyzz_1, g_xxyyzz_xyyyz_1, g_xxyyzz_xyyz_1, g_xxyyzz_xyyzz_1, g_xxyyzz_xyzz_1, g_xxyyzz_xyzzz_1, g_xxyyzz_yyyyy_1, g_xxyyzz_yyyyz_1, g_xxyyzz_yyyz_1, g_xxyyzz_yyyzz_1, g_xxyyzz_yyzz_1, g_xxyyzz_yyzzz_1, g_xxyyzz_yzzz_1, g_xxyyzz_yzzzz_1, g_xxyyzz_zzzzz_1, g_xyyzz_xxxyz_0, g_xyyzz_xxxyz_1, g_xyyzz_xxyyz_0, g_xyyzz_xxyyz_1, g_xyyzz_xxyzz_0, g_xyyzz_xxyzz_1, g_xyyzz_xyyyz_0, g_xyyzz_xyyyz_1, g_xyyzz_xyyzz_0, g_xyyzz_xyyzz_1, g_xyyzz_xyzzz_0, g_xyyzz_xyzzz_1, g_xyyzz_yyyyy_0, g_xyyzz_yyyyy_1, g_xyyzz_yyyyz_0, g_xyyzz_yyyyz_1, g_xyyzz_yyyzz_0, g_xyyzz_yyyzz_1, g_xyyzz_yyzzz_0, g_xyyzz_yyzzz_1, g_xyyzz_yzzzz_0, g_xyyzz_yzzzz_1, g_xyyzz_zzzzz_0, g_xyyzz_zzzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyzz_xxxxx_0[i] = g_xxxzz_xxxxx_0[i] * fbe_0 - g_xxxzz_xxxxx_1[i] * fz_be_0 + g_xxxyzz_xxxxx_1[i] * pa_y[i];

        g_xxxyyzz_xxxxy_0[i] = g_xxxyy_xxxxy_0[i] * fbe_0 - g_xxxyy_xxxxy_1[i] * fz_be_0 + g_xxxyyz_xxxxy_1[i] * pa_z[i];

        g_xxxyyzz_xxxxz_0[i] = g_xxxzz_xxxxz_0[i] * fbe_0 - g_xxxzz_xxxxz_1[i] * fz_be_0 + g_xxxyzz_xxxxz_1[i] * pa_y[i];

        g_xxxyyzz_xxxyy_0[i] = g_xxxyy_xxxyy_0[i] * fbe_0 - g_xxxyy_xxxyy_1[i] * fz_be_0 + g_xxxyyz_xxxyy_1[i] * pa_z[i];

        g_xxxyyzz_xxxyz_0[i] = 2.0 * g_xyyzz_xxxyz_0[i] * fbe_0 - 2.0 * g_xyyzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_xxyz_1[i] * fe_0 + g_xxyyzz_xxxyz_1[i] * pa_x[i];

        g_xxxyyzz_xxxzz_0[i] = g_xxxzz_xxxzz_0[i] * fbe_0 - g_xxxzz_xxxzz_1[i] * fz_be_0 + g_xxxyzz_xxxzz_1[i] * pa_y[i];

        g_xxxyyzz_xxyyy_0[i] = g_xxxyy_xxyyy_0[i] * fbe_0 - g_xxxyy_xxyyy_1[i] * fz_be_0 + g_xxxyyz_xxyyy_1[i] * pa_z[i];

        g_xxxyyzz_xxyyz_0[i] = 2.0 * g_xyyzz_xxyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_xyyz_1[i] * fe_0 + g_xxyyzz_xxyyz_1[i] * pa_x[i];

        g_xxxyyzz_xxyzz_0[i] = 2.0 * g_xyyzz_xxyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_xyzz_1[i] * fe_0 + g_xxyyzz_xxyzz_1[i] * pa_x[i];

        g_xxxyyzz_xxzzz_0[i] = g_xxxzz_xxzzz_0[i] * fbe_0 - g_xxxzz_xxzzz_1[i] * fz_be_0 + g_xxxyzz_xxzzz_1[i] * pa_y[i];

        g_xxxyyzz_xyyyy_0[i] = g_xxxyy_xyyyy_0[i] * fbe_0 - g_xxxyy_xyyyy_1[i] * fz_be_0 + g_xxxyyz_xyyyy_1[i] * pa_z[i];

        g_xxxyyzz_xyyyz_0[i] = 2.0 * g_xyyzz_xyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_xyyyz_1[i] * fz_be_0 + g_xxyyzz_yyyz_1[i] * fe_0 + g_xxyyzz_xyyyz_1[i] * pa_x[i];

        g_xxxyyzz_xyyzz_0[i] = 2.0 * g_xyyzz_xyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_xyyzz_1[i] * fz_be_0 + g_xxyyzz_yyzz_1[i] * fe_0 + g_xxyyzz_xyyzz_1[i] * pa_x[i];

        g_xxxyyzz_xyzzz_0[i] = 2.0 * g_xyyzz_xyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_xyzzz_1[i] * fz_be_0 + g_xxyyzz_yzzz_1[i] * fe_0 + g_xxyyzz_xyzzz_1[i] * pa_x[i];

        g_xxxyyzz_xzzzz_0[i] = g_xxxzz_xzzzz_0[i] * fbe_0 - g_xxxzz_xzzzz_1[i] * fz_be_0 + g_xxxyzz_xzzzz_1[i] * pa_y[i];

        g_xxxyyzz_yyyyy_0[i] = 2.0 * g_xyyzz_yyyyy_0[i] * fbe_0 - 2.0 * g_xyyzz_yyyyy_1[i] * fz_be_0 + g_xxyyzz_yyyyy_1[i] * pa_x[i];

        g_xxxyyzz_yyyyz_0[i] = 2.0 * g_xyyzz_yyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_yyyyz_1[i] * fz_be_0 + g_xxyyzz_yyyyz_1[i] * pa_x[i];

        g_xxxyyzz_yyyzz_0[i] = 2.0 * g_xyyzz_yyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_yyyzz_1[i] * fz_be_0 + g_xxyyzz_yyyzz_1[i] * pa_x[i];

        g_xxxyyzz_yyzzz_0[i] = 2.0 * g_xyyzz_yyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_yyzzz_1[i] * fz_be_0 + g_xxyyzz_yyzzz_1[i] * pa_x[i];

        g_xxxyyzz_yzzzz_0[i] = 2.0 * g_xyyzz_yzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_yzzzz_1[i] * fz_be_0 + g_xxyyzz_yzzzz_1[i] * pa_x[i];

        g_xxxyyzz_zzzzz_0[i] = 2.0 * g_xyyzz_zzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_zzzzz_1[i] * fz_be_0 + g_xxyyzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 273-294 components of targeted buffer : KH

    auto g_xxxyzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 273);

    auto g_xxxyzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 274);

    auto g_xxxyzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 275);

    auto g_xxxyzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 276);

    auto g_xxxyzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 277);

    auto g_xxxyzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 278);

    auto g_xxxyzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 279);

    auto g_xxxyzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 280);

    auto g_xxxyzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 281);

    auto g_xxxyzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 282);

    auto g_xxxyzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 283);

    auto g_xxxyzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 284);

    auto g_xxxyzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 285);

    auto g_xxxyzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 286);

    auto g_xxxyzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 287);

    auto g_xxxyzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 288);

    auto g_xxxyzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 289);

    auto g_xxxyzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 290);

    auto g_xxxyzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 291);

    auto g_xxxyzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 292);

    auto g_xxxyzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 293);

    #pragma omp simd aligned(g_xxxyzzz_xxxxx_0, g_xxxyzzz_xxxxy_0, g_xxxyzzz_xxxxz_0, g_xxxyzzz_xxxyy_0, g_xxxyzzz_xxxyz_0, g_xxxyzzz_xxxzz_0, g_xxxyzzz_xxyyy_0, g_xxxyzzz_xxyyz_0, g_xxxyzzz_xxyzz_0, g_xxxyzzz_xxzzz_0, g_xxxyzzz_xyyyy_0, g_xxxyzzz_xyyyz_0, g_xxxyzzz_xyyzz_0, g_xxxyzzz_xyzzz_0, g_xxxyzzz_xzzzz_0, g_xxxyzzz_yyyyy_0, g_xxxyzzz_yyyyz_0, g_xxxyzzz_yyyzz_0, g_xxxyzzz_yyzzz_0, g_xxxyzzz_yzzzz_0, g_xxxyzzz_zzzzz_0, g_xxxzzz_xxxx_1, g_xxxzzz_xxxxx_1, g_xxxzzz_xxxxy_1, g_xxxzzz_xxxxz_1, g_xxxzzz_xxxy_1, g_xxxzzz_xxxyy_1, g_xxxzzz_xxxyz_1, g_xxxzzz_xxxz_1, g_xxxzzz_xxxzz_1, g_xxxzzz_xxyy_1, g_xxxzzz_xxyyy_1, g_xxxzzz_xxyyz_1, g_xxxzzz_xxyz_1, g_xxxzzz_xxyzz_1, g_xxxzzz_xxzz_1, g_xxxzzz_xxzzz_1, g_xxxzzz_xyyy_1, g_xxxzzz_xyyyy_1, g_xxxzzz_xyyyz_1, g_xxxzzz_xyyz_1, g_xxxzzz_xyyzz_1, g_xxxzzz_xyzz_1, g_xxxzzz_xyzzz_1, g_xxxzzz_xzzz_1, g_xxxzzz_xzzzz_1, g_xxxzzz_yyyy_1, g_xxxzzz_yyyyy_1, g_xxxzzz_yyyyz_1, g_xxxzzz_yyyz_1, g_xxxzzz_yyyzz_1, g_xxxzzz_yyzz_1, g_xxxzzz_yyzzz_1, g_xxxzzz_yzzz_1, g_xxxzzz_yzzzz_1, g_xxxzzz_zzzz_1, g_xxxzzz_zzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyzzz_xxxxx_0[i] = g_xxxzzz_xxxxx_1[i] * pa_y[i];

        g_xxxyzzz_xxxxy_0[i] = g_xxxzzz_xxxx_1[i] * fe_0 + g_xxxzzz_xxxxy_1[i] * pa_y[i];

        g_xxxyzzz_xxxxz_0[i] = g_xxxzzz_xxxxz_1[i] * pa_y[i];

        g_xxxyzzz_xxxyy_0[i] = 2.0 * g_xxxzzz_xxxy_1[i] * fe_0 + g_xxxzzz_xxxyy_1[i] * pa_y[i];

        g_xxxyzzz_xxxyz_0[i] = g_xxxzzz_xxxz_1[i] * fe_0 + g_xxxzzz_xxxyz_1[i] * pa_y[i];

        g_xxxyzzz_xxxzz_0[i] = g_xxxzzz_xxxzz_1[i] * pa_y[i];

        g_xxxyzzz_xxyyy_0[i] = 3.0 * g_xxxzzz_xxyy_1[i] * fe_0 + g_xxxzzz_xxyyy_1[i] * pa_y[i];

        g_xxxyzzz_xxyyz_0[i] = 2.0 * g_xxxzzz_xxyz_1[i] * fe_0 + g_xxxzzz_xxyyz_1[i] * pa_y[i];

        g_xxxyzzz_xxyzz_0[i] = g_xxxzzz_xxzz_1[i] * fe_0 + g_xxxzzz_xxyzz_1[i] * pa_y[i];

        g_xxxyzzz_xxzzz_0[i] = g_xxxzzz_xxzzz_1[i] * pa_y[i];

        g_xxxyzzz_xyyyy_0[i] = 4.0 * g_xxxzzz_xyyy_1[i] * fe_0 + g_xxxzzz_xyyyy_1[i] * pa_y[i];

        g_xxxyzzz_xyyyz_0[i] = 3.0 * g_xxxzzz_xyyz_1[i] * fe_0 + g_xxxzzz_xyyyz_1[i] * pa_y[i];

        g_xxxyzzz_xyyzz_0[i] = 2.0 * g_xxxzzz_xyzz_1[i] * fe_0 + g_xxxzzz_xyyzz_1[i] * pa_y[i];

        g_xxxyzzz_xyzzz_0[i] = g_xxxzzz_xzzz_1[i] * fe_0 + g_xxxzzz_xyzzz_1[i] * pa_y[i];

        g_xxxyzzz_xzzzz_0[i] = g_xxxzzz_xzzzz_1[i] * pa_y[i];

        g_xxxyzzz_yyyyy_0[i] = 5.0 * g_xxxzzz_yyyy_1[i] * fe_0 + g_xxxzzz_yyyyy_1[i] * pa_y[i];

        g_xxxyzzz_yyyyz_0[i] = 4.0 * g_xxxzzz_yyyz_1[i] * fe_0 + g_xxxzzz_yyyyz_1[i] * pa_y[i];

        g_xxxyzzz_yyyzz_0[i] = 3.0 * g_xxxzzz_yyzz_1[i] * fe_0 + g_xxxzzz_yyyzz_1[i] * pa_y[i];

        g_xxxyzzz_yyzzz_0[i] = 2.0 * g_xxxzzz_yzzz_1[i] * fe_0 + g_xxxzzz_yyzzz_1[i] * pa_y[i];

        g_xxxyzzz_yzzzz_0[i] = g_xxxzzz_zzzz_1[i] * fe_0 + g_xxxzzz_yzzzz_1[i] * pa_y[i];

        g_xxxyzzz_zzzzz_0[i] = g_xxxzzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 294-315 components of targeted buffer : KH

    auto g_xxxzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 294);

    auto g_xxxzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 295);

    auto g_xxxzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 296);

    auto g_xxxzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 297);

    auto g_xxxzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 298);

    auto g_xxxzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 299);

    auto g_xxxzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 300);

    auto g_xxxzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 301);

    auto g_xxxzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 302);

    auto g_xxxzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 303);

    auto g_xxxzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 304);

    auto g_xxxzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 305);

    auto g_xxxzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 306);

    auto g_xxxzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 307);

    auto g_xxxzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 308);

    auto g_xxxzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 309);

    auto g_xxxzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 310);

    auto g_xxxzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 311);

    auto g_xxxzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 312);

    auto g_xxxzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 313);

    auto g_xxxzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 314);

    #pragma omp simd aligned(g_xxxzz_xxxxx_0, g_xxxzz_xxxxx_1, g_xxxzz_xxxxy_0, g_xxxzz_xxxxy_1, g_xxxzz_xxxyy_0, g_xxxzz_xxxyy_1, g_xxxzz_xxyyy_0, g_xxxzz_xxyyy_1, g_xxxzz_xyyyy_0, g_xxxzz_xyyyy_1, g_xxxzzz_xxxxx_1, g_xxxzzz_xxxxy_1, g_xxxzzz_xxxyy_1, g_xxxzzz_xxyyy_1, g_xxxzzz_xyyyy_1, g_xxxzzzz_xxxxx_0, g_xxxzzzz_xxxxy_0, g_xxxzzzz_xxxxz_0, g_xxxzzzz_xxxyy_0, g_xxxzzzz_xxxyz_0, g_xxxzzzz_xxxzz_0, g_xxxzzzz_xxyyy_0, g_xxxzzzz_xxyyz_0, g_xxxzzzz_xxyzz_0, g_xxxzzzz_xxzzz_0, g_xxxzzzz_xyyyy_0, g_xxxzzzz_xyyyz_0, g_xxxzzzz_xyyzz_0, g_xxxzzzz_xyzzz_0, g_xxxzzzz_xzzzz_0, g_xxxzzzz_yyyyy_0, g_xxxzzzz_yyyyz_0, g_xxxzzzz_yyyzz_0, g_xxxzzzz_yyzzz_0, g_xxxzzzz_yzzzz_0, g_xxxzzzz_zzzzz_0, g_xxzzzz_xxxxz_1, g_xxzzzz_xxxyz_1, g_xxzzzz_xxxz_1, g_xxzzzz_xxxzz_1, g_xxzzzz_xxyyz_1, g_xxzzzz_xxyz_1, g_xxzzzz_xxyzz_1, g_xxzzzz_xxzz_1, g_xxzzzz_xxzzz_1, g_xxzzzz_xyyyz_1, g_xxzzzz_xyyz_1, g_xxzzzz_xyyzz_1, g_xxzzzz_xyzz_1, g_xxzzzz_xyzzz_1, g_xxzzzz_xzzz_1, g_xxzzzz_xzzzz_1, g_xxzzzz_yyyyy_1, g_xxzzzz_yyyyz_1, g_xxzzzz_yyyz_1, g_xxzzzz_yyyzz_1, g_xxzzzz_yyzz_1, g_xxzzzz_yyzzz_1, g_xxzzzz_yzzz_1, g_xxzzzz_yzzzz_1, g_xxzzzz_zzzz_1, g_xxzzzz_zzzzz_1, g_xzzzz_xxxxz_0, g_xzzzz_xxxxz_1, g_xzzzz_xxxyz_0, g_xzzzz_xxxyz_1, g_xzzzz_xxxzz_0, g_xzzzz_xxxzz_1, g_xzzzz_xxyyz_0, g_xzzzz_xxyyz_1, g_xzzzz_xxyzz_0, g_xzzzz_xxyzz_1, g_xzzzz_xxzzz_0, g_xzzzz_xxzzz_1, g_xzzzz_xyyyz_0, g_xzzzz_xyyyz_1, g_xzzzz_xyyzz_0, g_xzzzz_xyyzz_1, g_xzzzz_xyzzz_0, g_xzzzz_xyzzz_1, g_xzzzz_xzzzz_0, g_xzzzz_xzzzz_1, g_xzzzz_yyyyy_0, g_xzzzz_yyyyy_1, g_xzzzz_yyyyz_0, g_xzzzz_yyyyz_1, g_xzzzz_yyyzz_0, g_xzzzz_yyyzz_1, g_xzzzz_yyzzz_0, g_xzzzz_yyzzz_1, g_xzzzz_yzzzz_0, g_xzzzz_yzzzz_1, g_xzzzz_zzzzz_0, g_xzzzz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzzzz_xxxxx_0[i] = 3.0 * g_xxxzz_xxxxx_0[i] * fbe_0 - 3.0 * g_xxxzz_xxxxx_1[i] * fz_be_0 + g_xxxzzz_xxxxx_1[i] * pa_z[i];

        g_xxxzzzz_xxxxy_0[i] = 3.0 * g_xxxzz_xxxxy_0[i] * fbe_0 - 3.0 * g_xxxzz_xxxxy_1[i] * fz_be_0 + g_xxxzzz_xxxxy_1[i] * pa_z[i];

        g_xxxzzzz_xxxxz_0[i] = 2.0 * g_xzzzz_xxxxz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_xxxz_1[i] * fe_0 + g_xxzzzz_xxxxz_1[i] * pa_x[i];

        g_xxxzzzz_xxxyy_0[i] = 3.0 * g_xxxzz_xxxyy_0[i] * fbe_0 - 3.0 * g_xxxzz_xxxyy_1[i] * fz_be_0 + g_xxxzzz_xxxyy_1[i] * pa_z[i];

        g_xxxzzzz_xxxyz_0[i] = 2.0 * g_xzzzz_xxxyz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_xxyz_1[i] * fe_0 + g_xxzzzz_xxxyz_1[i] * pa_x[i];

        g_xxxzzzz_xxxzz_0[i] = 2.0 * g_xzzzz_xxxzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_xxzz_1[i] * fe_0 + g_xxzzzz_xxxzz_1[i] * pa_x[i];

        g_xxxzzzz_xxyyy_0[i] = 3.0 * g_xxxzz_xxyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_xxyyy_1[i] * fz_be_0 + g_xxxzzz_xxyyy_1[i] * pa_z[i];

        g_xxxzzzz_xxyyz_0[i] = 2.0 * g_xzzzz_xxyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_xyyz_1[i] * fe_0 + g_xxzzzz_xxyyz_1[i] * pa_x[i];

        g_xxxzzzz_xxyzz_0[i] = 2.0 * g_xzzzz_xxyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_xyzz_1[i] * fe_0 + g_xxzzzz_xxyzz_1[i] * pa_x[i];

        g_xxxzzzz_xxzzz_0[i] = 2.0 * g_xzzzz_xxzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_xzzz_1[i] * fe_0 + g_xxzzzz_xxzzz_1[i] * pa_x[i];

        g_xxxzzzz_xyyyy_0[i] = 3.0 * g_xxxzz_xyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_xyyyy_1[i] * fz_be_0 + g_xxxzzz_xyyyy_1[i] * pa_z[i];

        g_xxxzzzz_xyyyz_0[i] = 2.0 * g_xzzzz_xyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_xyyyz_1[i] * fz_be_0 + g_xxzzzz_yyyz_1[i] * fe_0 + g_xxzzzz_xyyyz_1[i] * pa_x[i];

        g_xxxzzzz_xyyzz_0[i] = 2.0 * g_xzzzz_xyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xyyzz_1[i] * fz_be_0 + g_xxzzzz_yyzz_1[i] * fe_0 + g_xxzzzz_xyyzz_1[i] * pa_x[i];

        g_xxxzzzz_xyzzz_0[i] = 2.0 * g_xzzzz_xyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xyzzz_1[i] * fz_be_0 + g_xxzzzz_yzzz_1[i] * fe_0 + g_xxzzzz_xyzzz_1[i] * pa_x[i];

        g_xxxzzzz_xzzzz_0[i] = 2.0 * g_xzzzz_xzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xzzzz_1[i] * fz_be_0 + g_xxzzzz_zzzz_1[i] * fe_0 + g_xxzzzz_xzzzz_1[i] * pa_x[i];

        g_xxxzzzz_yyyyy_0[i] = 2.0 * g_xzzzz_yyyyy_0[i] * fbe_0 - 2.0 * g_xzzzz_yyyyy_1[i] * fz_be_0 + g_xxzzzz_yyyyy_1[i] * pa_x[i];

        g_xxxzzzz_yyyyz_0[i] = 2.0 * g_xzzzz_yyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_yyyyz_1[i] * fz_be_0 + g_xxzzzz_yyyyz_1[i] * pa_x[i];

        g_xxxzzzz_yyyzz_0[i] = 2.0 * g_xzzzz_yyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_yyyzz_1[i] * fz_be_0 + g_xxzzzz_yyyzz_1[i] * pa_x[i];

        g_xxxzzzz_yyzzz_0[i] = 2.0 * g_xzzzz_yyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_yyzzz_1[i] * fz_be_0 + g_xxzzzz_yyzzz_1[i] * pa_x[i];

        g_xxxzzzz_yzzzz_0[i] = 2.0 * g_xzzzz_yzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_yzzzz_1[i] * fz_be_0 + g_xxzzzz_yzzzz_1[i] * pa_x[i];

        g_xxxzzzz_zzzzz_0[i] = 2.0 * g_xzzzz_zzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_zzzzz_1[i] * fz_be_0 + g_xxzzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 315-336 components of targeted buffer : KH

    auto g_xxyyyyy_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 315);

    auto g_xxyyyyy_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 316);

    auto g_xxyyyyy_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 317);

    auto g_xxyyyyy_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 318);

    auto g_xxyyyyy_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 319);

    auto g_xxyyyyy_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 320);

    auto g_xxyyyyy_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 321);

    auto g_xxyyyyy_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 322);

    auto g_xxyyyyy_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 323);

    auto g_xxyyyyy_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 324);

    auto g_xxyyyyy_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 325);

    auto g_xxyyyyy_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 326);

    auto g_xxyyyyy_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 327);

    auto g_xxyyyyy_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 328);

    auto g_xxyyyyy_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 329);

    auto g_xxyyyyy_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 330);

    auto g_xxyyyyy_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 331);

    auto g_xxyyyyy_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 332);

    auto g_xxyyyyy_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 333);

    auto g_xxyyyyy_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 334);

    auto g_xxyyyyy_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 335);

    #pragma omp simd aligned(g_xxyyy_xxxxx_0, g_xxyyy_xxxxx_1, g_xxyyy_xxxxz_0, g_xxyyy_xxxxz_1, g_xxyyy_xxxzz_0, g_xxyyy_xxxzz_1, g_xxyyy_xxzzz_0, g_xxyyy_xxzzz_1, g_xxyyy_xzzzz_0, g_xxyyy_xzzzz_1, g_xxyyyy_xxxxx_1, g_xxyyyy_xxxxz_1, g_xxyyyy_xxxzz_1, g_xxyyyy_xxzzz_1, g_xxyyyy_xzzzz_1, g_xxyyyyy_xxxxx_0, g_xxyyyyy_xxxxy_0, g_xxyyyyy_xxxxz_0, g_xxyyyyy_xxxyy_0, g_xxyyyyy_xxxyz_0, g_xxyyyyy_xxxzz_0, g_xxyyyyy_xxyyy_0, g_xxyyyyy_xxyyz_0, g_xxyyyyy_xxyzz_0, g_xxyyyyy_xxzzz_0, g_xxyyyyy_xyyyy_0, g_xxyyyyy_xyyyz_0, g_xxyyyyy_xyyzz_0, g_xxyyyyy_xyzzz_0, g_xxyyyyy_xzzzz_0, g_xxyyyyy_yyyyy_0, g_xxyyyyy_yyyyz_0, g_xxyyyyy_yyyzz_0, g_xxyyyyy_yyzzz_0, g_xxyyyyy_yzzzz_0, g_xxyyyyy_zzzzz_0, g_xyyyyy_xxxxy_1, g_xyyyyy_xxxy_1, g_xyyyyy_xxxyy_1, g_xyyyyy_xxxyz_1, g_xyyyyy_xxyy_1, g_xyyyyy_xxyyy_1, g_xyyyyy_xxyyz_1, g_xyyyyy_xxyz_1, g_xyyyyy_xxyzz_1, g_xyyyyy_xyyy_1, g_xyyyyy_xyyyy_1, g_xyyyyy_xyyyz_1, g_xyyyyy_xyyz_1, g_xyyyyy_xyyzz_1, g_xyyyyy_xyzz_1, g_xyyyyy_xyzzz_1, g_xyyyyy_yyyy_1, g_xyyyyy_yyyyy_1, g_xyyyyy_yyyyz_1, g_xyyyyy_yyyz_1, g_xyyyyy_yyyzz_1, g_xyyyyy_yyzz_1, g_xyyyyy_yyzzz_1, g_xyyyyy_yzzz_1, g_xyyyyy_yzzzz_1, g_xyyyyy_zzzzz_1, g_yyyyy_xxxxy_0, g_yyyyy_xxxxy_1, g_yyyyy_xxxyy_0, g_yyyyy_xxxyy_1, g_yyyyy_xxxyz_0, g_yyyyy_xxxyz_1, g_yyyyy_xxyyy_0, g_yyyyy_xxyyy_1, g_yyyyy_xxyyz_0, g_yyyyy_xxyyz_1, g_yyyyy_xxyzz_0, g_yyyyy_xxyzz_1, g_yyyyy_xyyyy_0, g_yyyyy_xyyyy_1, g_yyyyy_xyyyz_0, g_yyyyy_xyyyz_1, g_yyyyy_xyyzz_0, g_yyyyy_xyyzz_1, g_yyyyy_xyzzz_0, g_yyyyy_xyzzz_1, g_yyyyy_yyyyy_0, g_yyyyy_yyyyy_1, g_yyyyy_yyyyz_0, g_yyyyy_yyyyz_1, g_yyyyy_yyyzz_0, g_yyyyy_yyyzz_1, g_yyyyy_yyzzz_0, g_yyyyy_yyzzz_1, g_yyyyy_yzzzz_0, g_yyyyy_yzzzz_1, g_yyyyy_zzzzz_0, g_yyyyy_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyyy_xxxxx_0[i] = 4.0 * g_xxyyy_xxxxx_0[i] * fbe_0 - 4.0 * g_xxyyy_xxxxx_1[i] * fz_be_0 + g_xxyyyy_xxxxx_1[i] * pa_y[i];

        g_xxyyyyy_xxxxy_0[i] = g_yyyyy_xxxxy_0[i] * fbe_0 - g_yyyyy_xxxxy_1[i] * fz_be_0 + 4.0 * g_xyyyyy_xxxy_1[i] * fe_0 + g_xyyyyy_xxxxy_1[i] * pa_x[i];

        g_xxyyyyy_xxxxz_0[i] = 4.0 * g_xxyyy_xxxxz_0[i] * fbe_0 - 4.0 * g_xxyyy_xxxxz_1[i] * fz_be_0 + g_xxyyyy_xxxxz_1[i] * pa_y[i];

        g_xxyyyyy_xxxyy_0[i] = g_yyyyy_xxxyy_0[i] * fbe_0 - g_yyyyy_xxxyy_1[i] * fz_be_0 + 3.0 * g_xyyyyy_xxyy_1[i] * fe_0 + g_xyyyyy_xxxyy_1[i] * pa_x[i];

        g_xxyyyyy_xxxyz_0[i] = g_yyyyy_xxxyz_0[i] * fbe_0 - g_yyyyy_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_xxyz_1[i] * fe_0 + g_xyyyyy_xxxyz_1[i] * pa_x[i];

        g_xxyyyyy_xxxzz_0[i] = 4.0 * g_xxyyy_xxxzz_0[i] * fbe_0 - 4.0 * g_xxyyy_xxxzz_1[i] * fz_be_0 + g_xxyyyy_xxxzz_1[i] * pa_y[i];

        g_xxyyyyy_xxyyy_0[i] = g_yyyyy_xxyyy_0[i] * fbe_0 - g_yyyyy_xxyyy_1[i] * fz_be_0 + 2.0 * g_xyyyyy_xyyy_1[i] * fe_0 + g_xyyyyy_xxyyy_1[i] * pa_x[i];

        g_xxyyyyy_xxyyz_0[i] = g_yyyyy_xxyyz_0[i] * fbe_0 - g_yyyyy_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_xyyz_1[i] * fe_0 + g_xyyyyy_xxyyz_1[i] * pa_x[i];

        g_xxyyyyy_xxyzz_0[i] = g_yyyyy_xxyzz_0[i] * fbe_0 - g_yyyyy_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_xyzz_1[i] * fe_0 + g_xyyyyy_xxyzz_1[i] * pa_x[i];

        g_xxyyyyy_xxzzz_0[i] = 4.0 * g_xxyyy_xxzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_xxzzz_1[i] * fz_be_0 + g_xxyyyy_xxzzz_1[i] * pa_y[i];

        g_xxyyyyy_xyyyy_0[i] = g_yyyyy_xyyyy_0[i] * fbe_0 - g_yyyyy_xyyyy_1[i] * fz_be_0 + g_xyyyyy_yyyy_1[i] * fe_0 + g_xyyyyy_xyyyy_1[i] * pa_x[i];

        g_xxyyyyy_xyyyz_0[i] = g_yyyyy_xyyyz_0[i] * fbe_0 - g_yyyyy_xyyyz_1[i] * fz_be_0 + g_xyyyyy_yyyz_1[i] * fe_0 + g_xyyyyy_xyyyz_1[i] * pa_x[i];

        g_xxyyyyy_xyyzz_0[i] = g_yyyyy_xyyzz_0[i] * fbe_0 - g_yyyyy_xyyzz_1[i] * fz_be_0 + g_xyyyyy_yyzz_1[i] * fe_0 + g_xyyyyy_xyyzz_1[i] * pa_x[i];

        g_xxyyyyy_xyzzz_0[i] = g_yyyyy_xyzzz_0[i] * fbe_0 - g_yyyyy_xyzzz_1[i] * fz_be_0 + g_xyyyyy_yzzz_1[i] * fe_0 + g_xyyyyy_xyzzz_1[i] * pa_x[i];

        g_xxyyyyy_xzzzz_0[i] = 4.0 * g_xxyyy_xzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_xzzzz_1[i] * fz_be_0 + g_xxyyyy_xzzzz_1[i] * pa_y[i];

        g_xxyyyyy_yyyyy_0[i] = g_yyyyy_yyyyy_0[i] * fbe_0 - g_yyyyy_yyyyy_1[i] * fz_be_0 + g_xyyyyy_yyyyy_1[i] * pa_x[i];

        g_xxyyyyy_yyyyz_0[i] = g_yyyyy_yyyyz_0[i] * fbe_0 - g_yyyyy_yyyyz_1[i] * fz_be_0 + g_xyyyyy_yyyyz_1[i] * pa_x[i];

        g_xxyyyyy_yyyzz_0[i] = g_yyyyy_yyyzz_0[i] * fbe_0 - g_yyyyy_yyyzz_1[i] * fz_be_0 + g_xyyyyy_yyyzz_1[i] * pa_x[i];

        g_xxyyyyy_yyzzz_0[i] = g_yyyyy_yyzzz_0[i] * fbe_0 - g_yyyyy_yyzzz_1[i] * fz_be_0 + g_xyyyyy_yyzzz_1[i] * pa_x[i];

        g_xxyyyyy_yzzzz_0[i] = g_yyyyy_yzzzz_0[i] * fbe_0 - g_yyyyy_yzzzz_1[i] * fz_be_0 + g_xyyyyy_yzzzz_1[i] * pa_x[i];

        g_xxyyyyy_zzzzz_0[i] = g_yyyyy_zzzzz_0[i] * fbe_0 - g_yyyyy_zzzzz_1[i] * fz_be_0 + g_xyyyyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 336-357 components of targeted buffer : KH

    auto g_xxyyyyz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 336);

    auto g_xxyyyyz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 337);

    auto g_xxyyyyz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 338);

    auto g_xxyyyyz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 339);

    auto g_xxyyyyz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 340);

    auto g_xxyyyyz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 341);

    auto g_xxyyyyz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 342);

    auto g_xxyyyyz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 343);

    auto g_xxyyyyz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 344);

    auto g_xxyyyyz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 345);

    auto g_xxyyyyz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 346);

    auto g_xxyyyyz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 347);

    auto g_xxyyyyz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 348);

    auto g_xxyyyyz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 349);

    auto g_xxyyyyz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 350);

    auto g_xxyyyyz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 351);

    auto g_xxyyyyz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 352);

    auto g_xxyyyyz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 353);

    auto g_xxyyyyz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 354);

    auto g_xxyyyyz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 355);

    auto g_xxyyyyz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 356);

    #pragma omp simd aligned(g_xxyyyy_xxxx_1, g_xxyyyy_xxxxx_1, g_xxyyyy_xxxxy_1, g_xxyyyy_xxxxz_1, g_xxyyyy_xxxy_1, g_xxyyyy_xxxyy_1, g_xxyyyy_xxxyz_1, g_xxyyyy_xxxz_1, g_xxyyyy_xxxzz_1, g_xxyyyy_xxyy_1, g_xxyyyy_xxyyy_1, g_xxyyyy_xxyyz_1, g_xxyyyy_xxyz_1, g_xxyyyy_xxyzz_1, g_xxyyyy_xxzz_1, g_xxyyyy_xxzzz_1, g_xxyyyy_xyyy_1, g_xxyyyy_xyyyy_1, g_xxyyyy_xyyyz_1, g_xxyyyy_xyyz_1, g_xxyyyy_xyyzz_1, g_xxyyyy_xyzz_1, g_xxyyyy_xyzzz_1, g_xxyyyy_xzzz_1, g_xxyyyy_xzzzz_1, g_xxyyyy_yyyy_1, g_xxyyyy_yyyyy_1, g_xxyyyy_yyyyz_1, g_xxyyyy_yyyz_1, g_xxyyyy_yyyzz_1, g_xxyyyy_yyzz_1, g_xxyyyy_yyzzz_1, g_xxyyyy_yzzz_1, g_xxyyyy_yzzzz_1, g_xxyyyy_zzzz_1, g_xxyyyy_zzzzz_1, g_xxyyyyz_xxxxx_0, g_xxyyyyz_xxxxy_0, g_xxyyyyz_xxxxz_0, g_xxyyyyz_xxxyy_0, g_xxyyyyz_xxxyz_0, g_xxyyyyz_xxxzz_0, g_xxyyyyz_xxyyy_0, g_xxyyyyz_xxyyz_0, g_xxyyyyz_xxyzz_0, g_xxyyyyz_xxzzz_0, g_xxyyyyz_xyyyy_0, g_xxyyyyz_xyyyz_0, g_xxyyyyz_xyyzz_0, g_xxyyyyz_xyzzz_0, g_xxyyyyz_xzzzz_0, g_xxyyyyz_yyyyy_0, g_xxyyyyz_yyyyz_0, g_xxyyyyz_yyyzz_0, g_xxyyyyz_yyzzz_0, g_xxyyyyz_yzzzz_0, g_xxyyyyz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyyyz_xxxxx_0[i] = g_xxyyyy_xxxxx_1[i] * pa_z[i];

        g_xxyyyyz_xxxxy_0[i] = g_xxyyyy_xxxxy_1[i] * pa_z[i];

        g_xxyyyyz_xxxxz_0[i] = g_xxyyyy_xxxx_1[i] * fe_0 + g_xxyyyy_xxxxz_1[i] * pa_z[i];

        g_xxyyyyz_xxxyy_0[i] = g_xxyyyy_xxxyy_1[i] * pa_z[i];

        g_xxyyyyz_xxxyz_0[i] = g_xxyyyy_xxxy_1[i] * fe_0 + g_xxyyyy_xxxyz_1[i] * pa_z[i];

        g_xxyyyyz_xxxzz_0[i] = 2.0 * g_xxyyyy_xxxz_1[i] * fe_0 + g_xxyyyy_xxxzz_1[i] * pa_z[i];

        g_xxyyyyz_xxyyy_0[i] = g_xxyyyy_xxyyy_1[i] * pa_z[i];

        g_xxyyyyz_xxyyz_0[i] = g_xxyyyy_xxyy_1[i] * fe_0 + g_xxyyyy_xxyyz_1[i] * pa_z[i];

        g_xxyyyyz_xxyzz_0[i] = 2.0 * g_xxyyyy_xxyz_1[i] * fe_0 + g_xxyyyy_xxyzz_1[i] * pa_z[i];

        g_xxyyyyz_xxzzz_0[i] = 3.0 * g_xxyyyy_xxzz_1[i] * fe_0 + g_xxyyyy_xxzzz_1[i] * pa_z[i];

        g_xxyyyyz_xyyyy_0[i] = g_xxyyyy_xyyyy_1[i] * pa_z[i];

        g_xxyyyyz_xyyyz_0[i] = g_xxyyyy_xyyy_1[i] * fe_0 + g_xxyyyy_xyyyz_1[i] * pa_z[i];

        g_xxyyyyz_xyyzz_0[i] = 2.0 * g_xxyyyy_xyyz_1[i] * fe_0 + g_xxyyyy_xyyzz_1[i] * pa_z[i];

        g_xxyyyyz_xyzzz_0[i] = 3.0 * g_xxyyyy_xyzz_1[i] * fe_0 + g_xxyyyy_xyzzz_1[i] * pa_z[i];

        g_xxyyyyz_xzzzz_0[i] = 4.0 * g_xxyyyy_xzzz_1[i] * fe_0 + g_xxyyyy_xzzzz_1[i] * pa_z[i];

        g_xxyyyyz_yyyyy_0[i] = g_xxyyyy_yyyyy_1[i] * pa_z[i];

        g_xxyyyyz_yyyyz_0[i] = g_xxyyyy_yyyy_1[i] * fe_0 + g_xxyyyy_yyyyz_1[i] * pa_z[i];

        g_xxyyyyz_yyyzz_0[i] = 2.0 * g_xxyyyy_yyyz_1[i] * fe_0 + g_xxyyyy_yyyzz_1[i] * pa_z[i];

        g_xxyyyyz_yyzzz_0[i] = 3.0 * g_xxyyyy_yyzz_1[i] * fe_0 + g_xxyyyy_yyzzz_1[i] * pa_z[i];

        g_xxyyyyz_yzzzz_0[i] = 4.0 * g_xxyyyy_yzzz_1[i] * fe_0 + g_xxyyyy_yzzzz_1[i] * pa_z[i];

        g_xxyyyyz_zzzzz_0[i] = 5.0 * g_xxyyyy_zzzz_1[i] * fe_0 + g_xxyyyy_zzzzz_1[i] * pa_z[i];
    }

    // Set up 357-378 components of targeted buffer : KH

    auto g_xxyyyzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 357);

    auto g_xxyyyzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 358);

    auto g_xxyyyzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 359);

    auto g_xxyyyzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 360);

    auto g_xxyyyzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 361);

    auto g_xxyyyzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 362);

    auto g_xxyyyzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 363);

    auto g_xxyyyzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 364);

    auto g_xxyyyzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 365);

    auto g_xxyyyzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 366);

    auto g_xxyyyzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 367);

    auto g_xxyyyzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 368);

    auto g_xxyyyzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 369);

    auto g_xxyyyzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 370);

    auto g_xxyyyzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 371);

    auto g_xxyyyzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 372);

    auto g_xxyyyzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 373);

    auto g_xxyyyzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 374);

    auto g_xxyyyzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 375);

    auto g_xxyyyzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 376);

    auto g_xxyyyzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 377);

    #pragma omp simd aligned(g_xxyyy_xxxxy_0, g_xxyyy_xxxxy_1, g_xxyyy_xxxyy_0, g_xxyyy_xxxyy_1, g_xxyyy_xxyyy_0, g_xxyyy_xxyyy_1, g_xxyyy_xyyyy_0, g_xxyyy_xyyyy_1, g_xxyyyz_xxxxy_1, g_xxyyyz_xxxyy_1, g_xxyyyz_xxyyy_1, g_xxyyyz_xyyyy_1, g_xxyyyzz_xxxxx_0, g_xxyyyzz_xxxxy_0, g_xxyyyzz_xxxxz_0, g_xxyyyzz_xxxyy_0, g_xxyyyzz_xxxyz_0, g_xxyyyzz_xxxzz_0, g_xxyyyzz_xxyyy_0, g_xxyyyzz_xxyyz_0, g_xxyyyzz_xxyzz_0, g_xxyyyzz_xxzzz_0, g_xxyyyzz_xyyyy_0, g_xxyyyzz_xyyyz_0, g_xxyyyzz_xyyzz_0, g_xxyyyzz_xyzzz_0, g_xxyyyzz_xzzzz_0, g_xxyyyzz_yyyyy_0, g_xxyyyzz_yyyyz_0, g_xxyyyzz_yyyzz_0, g_xxyyyzz_yyzzz_0, g_xxyyyzz_yzzzz_0, g_xxyyyzz_zzzzz_0, g_xxyyzz_xxxxx_1, g_xxyyzz_xxxxz_1, g_xxyyzz_xxxzz_1, g_xxyyzz_xxzzz_1, g_xxyyzz_xzzzz_1, g_xxyzz_xxxxx_0, g_xxyzz_xxxxx_1, g_xxyzz_xxxxz_0, g_xxyzz_xxxxz_1, g_xxyzz_xxxzz_0, g_xxyzz_xxxzz_1, g_xxyzz_xxzzz_0, g_xxyzz_xxzzz_1, g_xxyzz_xzzzz_0, g_xxyzz_xzzzz_1, g_xyyyzz_xxxyz_1, g_xyyyzz_xxyyz_1, g_xyyyzz_xxyz_1, g_xyyyzz_xxyzz_1, g_xyyyzz_xyyyz_1, g_xyyyzz_xyyz_1, g_xyyyzz_xyyzz_1, g_xyyyzz_xyzz_1, g_xyyyzz_xyzzz_1, g_xyyyzz_yyyyy_1, g_xyyyzz_yyyyz_1, g_xyyyzz_yyyz_1, g_xyyyzz_yyyzz_1, g_xyyyzz_yyzz_1, g_xyyyzz_yyzzz_1, g_xyyyzz_yzzz_1, g_xyyyzz_yzzzz_1, g_xyyyzz_zzzzz_1, g_yyyzz_xxxyz_0, g_yyyzz_xxxyz_1, g_yyyzz_xxyyz_0, g_yyyzz_xxyyz_1, g_yyyzz_xxyzz_0, g_yyyzz_xxyzz_1, g_yyyzz_xyyyz_0, g_yyyzz_xyyyz_1, g_yyyzz_xyyzz_0, g_yyyzz_xyyzz_1, g_yyyzz_xyzzz_0, g_yyyzz_xyzzz_1, g_yyyzz_yyyyy_0, g_yyyzz_yyyyy_1, g_yyyzz_yyyyz_0, g_yyyzz_yyyyz_1, g_yyyzz_yyyzz_0, g_yyyzz_yyyzz_1, g_yyyzz_yyzzz_0, g_yyyzz_yyzzz_1, g_yyyzz_yzzzz_0, g_yyyzz_yzzzz_1, g_yyyzz_zzzzz_0, g_yyyzz_zzzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyzz_xxxxx_0[i] = 2.0 * g_xxyzz_xxxxx_0[i] * fbe_0 - 2.0 * g_xxyzz_xxxxx_1[i] * fz_be_0 + g_xxyyzz_xxxxx_1[i] * pa_y[i];

        g_xxyyyzz_xxxxy_0[i] = g_xxyyy_xxxxy_0[i] * fbe_0 - g_xxyyy_xxxxy_1[i] * fz_be_0 + g_xxyyyz_xxxxy_1[i] * pa_z[i];

        g_xxyyyzz_xxxxz_0[i] = 2.0 * g_xxyzz_xxxxz_0[i] * fbe_0 - 2.0 * g_xxyzz_xxxxz_1[i] * fz_be_0 + g_xxyyzz_xxxxz_1[i] * pa_y[i];

        g_xxyyyzz_xxxyy_0[i] = g_xxyyy_xxxyy_0[i] * fbe_0 - g_xxyyy_xxxyy_1[i] * fz_be_0 + g_xxyyyz_xxxyy_1[i] * pa_z[i];

        g_xxyyyzz_xxxyz_0[i] = g_yyyzz_xxxyz_0[i] * fbe_0 - g_yyyzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_xxyz_1[i] * fe_0 + g_xyyyzz_xxxyz_1[i] * pa_x[i];

        g_xxyyyzz_xxxzz_0[i] = 2.0 * g_xxyzz_xxxzz_0[i] * fbe_0 - 2.0 * g_xxyzz_xxxzz_1[i] * fz_be_0 + g_xxyyzz_xxxzz_1[i] * pa_y[i];

        g_xxyyyzz_xxyyy_0[i] = g_xxyyy_xxyyy_0[i] * fbe_0 - g_xxyyy_xxyyy_1[i] * fz_be_0 + g_xxyyyz_xxyyy_1[i] * pa_z[i];

        g_xxyyyzz_xxyyz_0[i] = g_yyyzz_xxyyz_0[i] * fbe_0 - g_yyyzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_xyyz_1[i] * fe_0 + g_xyyyzz_xxyyz_1[i] * pa_x[i];

        g_xxyyyzz_xxyzz_0[i] = g_yyyzz_xxyzz_0[i] * fbe_0 - g_yyyzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_xyzz_1[i] * fe_0 + g_xyyyzz_xxyzz_1[i] * pa_x[i];

        g_xxyyyzz_xxzzz_0[i] = 2.0 * g_xxyzz_xxzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_xxzzz_1[i] * fz_be_0 + g_xxyyzz_xxzzz_1[i] * pa_y[i];

        g_xxyyyzz_xyyyy_0[i] = g_xxyyy_xyyyy_0[i] * fbe_0 - g_xxyyy_xyyyy_1[i] * fz_be_0 + g_xxyyyz_xyyyy_1[i] * pa_z[i];

        g_xxyyyzz_xyyyz_0[i] = g_yyyzz_xyyyz_0[i] * fbe_0 - g_yyyzz_xyyyz_1[i] * fz_be_0 + g_xyyyzz_yyyz_1[i] * fe_0 + g_xyyyzz_xyyyz_1[i] * pa_x[i];

        g_xxyyyzz_xyyzz_0[i] = g_yyyzz_xyyzz_0[i] * fbe_0 - g_yyyzz_xyyzz_1[i] * fz_be_0 + g_xyyyzz_yyzz_1[i] * fe_0 + g_xyyyzz_xyyzz_1[i] * pa_x[i];

        g_xxyyyzz_xyzzz_0[i] = g_yyyzz_xyzzz_0[i] * fbe_0 - g_yyyzz_xyzzz_1[i] * fz_be_0 + g_xyyyzz_yzzz_1[i] * fe_0 + g_xyyyzz_xyzzz_1[i] * pa_x[i];

        g_xxyyyzz_xzzzz_0[i] = 2.0 * g_xxyzz_xzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_xzzzz_1[i] * fz_be_0 + g_xxyyzz_xzzzz_1[i] * pa_y[i];

        g_xxyyyzz_yyyyy_0[i] = g_yyyzz_yyyyy_0[i] * fbe_0 - g_yyyzz_yyyyy_1[i] * fz_be_0 + g_xyyyzz_yyyyy_1[i] * pa_x[i];

        g_xxyyyzz_yyyyz_0[i] = g_yyyzz_yyyyz_0[i] * fbe_0 - g_yyyzz_yyyyz_1[i] * fz_be_0 + g_xyyyzz_yyyyz_1[i] * pa_x[i];

        g_xxyyyzz_yyyzz_0[i] = g_yyyzz_yyyzz_0[i] * fbe_0 - g_yyyzz_yyyzz_1[i] * fz_be_0 + g_xyyyzz_yyyzz_1[i] * pa_x[i];

        g_xxyyyzz_yyzzz_0[i] = g_yyyzz_yyzzz_0[i] * fbe_0 - g_yyyzz_yyzzz_1[i] * fz_be_0 + g_xyyyzz_yyzzz_1[i] * pa_x[i];

        g_xxyyyzz_yzzzz_0[i] = g_yyyzz_yzzzz_0[i] * fbe_0 - g_yyyzz_yzzzz_1[i] * fz_be_0 + g_xyyyzz_yzzzz_1[i] * pa_x[i];

        g_xxyyyzz_zzzzz_0[i] = g_yyyzz_zzzzz_0[i] * fbe_0 - g_yyyzz_zzzzz_1[i] * fz_be_0 + g_xyyyzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 378-399 components of targeted buffer : KH

    auto g_xxyyzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 378);

    auto g_xxyyzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 379);

    auto g_xxyyzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 380);

    auto g_xxyyzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 381);

    auto g_xxyyzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 382);

    auto g_xxyyzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 383);

    auto g_xxyyzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 384);

    auto g_xxyyzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 385);

    auto g_xxyyzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 386);

    auto g_xxyyzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 387);

    auto g_xxyyzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 388);

    auto g_xxyyzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 389);

    auto g_xxyyzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 390);

    auto g_xxyyzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 391);

    auto g_xxyyzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 392);

    auto g_xxyyzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 393);

    auto g_xxyyzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 394);

    auto g_xxyyzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 395);

    auto g_xxyyzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 396);

    auto g_xxyyzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 397);

    auto g_xxyyzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 398);

    #pragma omp simd aligned(g_xxyyz_xxxxy_0, g_xxyyz_xxxxy_1, g_xxyyz_xxxyy_0, g_xxyyz_xxxyy_1, g_xxyyz_xxyyy_0, g_xxyyz_xxyyy_1, g_xxyyz_xyyyy_0, g_xxyyz_xyyyy_1, g_xxyyzz_xxxxy_1, g_xxyyzz_xxxyy_1, g_xxyyzz_xxyyy_1, g_xxyyzz_xyyyy_1, g_xxyyzzz_xxxxx_0, g_xxyyzzz_xxxxy_0, g_xxyyzzz_xxxxz_0, g_xxyyzzz_xxxyy_0, g_xxyyzzz_xxxyz_0, g_xxyyzzz_xxxzz_0, g_xxyyzzz_xxyyy_0, g_xxyyzzz_xxyyz_0, g_xxyyzzz_xxyzz_0, g_xxyyzzz_xxzzz_0, g_xxyyzzz_xyyyy_0, g_xxyyzzz_xyyyz_0, g_xxyyzzz_xyyzz_0, g_xxyyzzz_xyzzz_0, g_xxyyzzz_xzzzz_0, g_xxyyzzz_yyyyy_0, g_xxyyzzz_yyyyz_0, g_xxyyzzz_yyyzz_0, g_xxyyzzz_yyzzz_0, g_xxyyzzz_yzzzz_0, g_xxyyzzz_zzzzz_0, g_xxyzzz_xxxxx_1, g_xxyzzz_xxxxz_1, g_xxyzzz_xxxzz_1, g_xxyzzz_xxzzz_1, g_xxyzzz_xzzzz_1, g_xxzzz_xxxxx_0, g_xxzzz_xxxxx_1, g_xxzzz_xxxxz_0, g_xxzzz_xxxxz_1, g_xxzzz_xxxzz_0, g_xxzzz_xxxzz_1, g_xxzzz_xxzzz_0, g_xxzzz_xxzzz_1, g_xxzzz_xzzzz_0, g_xxzzz_xzzzz_1, g_xyyzzz_xxxyz_1, g_xyyzzz_xxyyz_1, g_xyyzzz_xxyz_1, g_xyyzzz_xxyzz_1, g_xyyzzz_xyyyz_1, g_xyyzzz_xyyz_1, g_xyyzzz_xyyzz_1, g_xyyzzz_xyzz_1, g_xyyzzz_xyzzz_1, g_xyyzzz_yyyyy_1, g_xyyzzz_yyyyz_1, g_xyyzzz_yyyz_1, g_xyyzzz_yyyzz_1, g_xyyzzz_yyzz_1, g_xyyzzz_yyzzz_1, g_xyyzzz_yzzz_1, g_xyyzzz_yzzzz_1, g_xyyzzz_zzzzz_1, g_yyzzz_xxxyz_0, g_yyzzz_xxxyz_1, g_yyzzz_xxyyz_0, g_yyzzz_xxyyz_1, g_yyzzz_xxyzz_0, g_yyzzz_xxyzz_1, g_yyzzz_xyyyz_0, g_yyzzz_xyyyz_1, g_yyzzz_xyyzz_0, g_yyzzz_xyyzz_1, g_yyzzz_xyzzz_0, g_yyzzz_xyzzz_1, g_yyzzz_yyyyy_0, g_yyzzz_yyyyy_1, g_yyzzz_yyyyz_0, g_yyzzz_yyyyz_1, g_yyzzz_yyyzz_0, g_yyzzz_yyyzz_1, g_yyzzz_yyzzz_0, g_yyzzz_yyzzz_1, g_yyzzz_yzzzz_0, g_yyzzz_yzzzz_1, g_yyzzz_zzzzz_0, g_yyzzz_zzzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyzzz_xxxxx_0[i] = g_xxzzz_xxxxx_0[i] * fbe_0 - g_xxzzz_xxxxx_1[i] * fz_be_0 + g_xxyzzz_xxxxx_1[i] * pa_y[i];

        g_xxyyzzz_xxxxy_0[i] = 2.0 * g_xxyyz_xxxxy_0[i] * fbe_0 - 2.0 * g_xxyyz_xxxxy_1[i] * fz_be_0 + g_xxyyzz_xxxxy_1[i] * pa_z[i];

        g_xxyyzzz_xxxxz_0[i] = g_xxzzz_xxxxz_0[i] * fbe_0 - g_xxzzz_xxxxz_1[i] * fz_be_0 + g_xxyzzz_xxxxz_1[i] * pa_y[i];

        g_xxyyzzz_xxxyy_0[i] = 2.0 * g_xxyyz_xxxyy_0[i] * fbe_0 - 2.0 * g_xxyyz_xxxyy_1[i] * fz_be_0 + g_xxyyzz_xxxyy_1[i] * pa_z[i];

        g_xxyyzzz_xxxyz_0[i] = g_yyzzz_xxxyz_0[i] * fbe_0 - g_yyzzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_xxyz_1[i] * fe_0 + g_xyyzzz_xxxyz_1[i] * pa_x[i];

        g_xxyyzzz_xxxzz_0[i] = g_xxzzz_xxxzz_0[i] * fbe_0 - g_xxzzz_xxxzz_1[i] * fz_be_0 + g_xxyzzz_xxxzz_1[i] * pa_y[i];

        g_xxyyzzz_xxyyy_0[i] = 2.0 * g_xxyyz_xxyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_xxyyy_1[i] * fz_be_0 + g_xxyyzz_xxyyy_1[i] * pa_z[i];

        g_xxyyzzz_xxyyz_0[i] = g_yyzzz_xxyyz_0[i] * fbe_0 - g_yyzzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_xyyz_1[i] * fe_0 + g_xyyzzz_xxyyz_1[i] * pa_x[i];

        g_xxyyzzz_xxyzz_0[i] = g_yyzzz_xxyzz_0[i] * fbe_0 - g_yyzzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_xyzz_1[i] * fe_0 + g_xyyzzz_xxyzz_1[i] * pa_x[i];

        g_xxyyzzz_xxzzz_0[i] = g_xxzzz_xxzzz_0[i] * fbe_0 - g_xxzzz_xxzzz_1[i] * fz_be_0 + g_xxyzzz_xxzzz_1[i] * pa_y[i];

        g_xxyyzzz_xyyyy_0[i] = 2.0 * g_xxyyz_xyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_xyyyy_1[i] * fz_be_0 + g_xxyyzz_xyyyy_1[i] * pa_z[i];

        g_xxyyzzz_xyyyz_0[i] = g_yyzzz_xyyyz_0[i] * fbe_0 - g_yyzzz_xyyyz_1[i] * fz_be_0 + g_xyyzzz_yyyz_1[i] * fe_0 + g_xyyzzz_xyyyz_1[i] * pa_x[i];

        g_xxyyzzz_xyyzz_0[i] = g_yyzzz_xyyzz_0[i] * fbe_0 - g_yyzzz_xyyzz_1[i] * fz_be_0 + g_xyyzzz_yyzz_1[i] * fe_0 + g_xyyzzz_xyyzz_1[i] * pa_x[i];

        g_xxyyzzz_xyzzz_0[i] = g_yyzzz_xyzzz_0[i] * fbe_0 - g_yyzzz_xyzzz_1[i] * fz_be_0 + g_xyyzzz_yzzz_1[i] * fe_0 + g_xyyzzz_xyzzz_1[i] * pa_x[i];

        g_xxyyzzz_xzzzz_0[i] = g_xxzzz_xzzzz_0[i] * fbe_0 - g_xxzzz_xzzzz_1[i] * fz_be_0 + g_xxyzzz_xzzzz_1[i] * pa_y[i];

        g_xxyyzzz_yyyyy_0[i] = g_yyzzz_yyyyy_0[i] * fbe_0 - g_yyzzz_yyyyy_1[i] * fz_be_0 + g_xyyzzz_yyyyy_1[i] * pa_x[i];

        g_xxyyzzz_yyyyz_0[i] = g_yyzzz_yyyyz_0[i] * fbe_0 - g_yyzzz_yyyyz_1[i] * fz_be_0 + g_xyyzzz_yyyyz_1[i] * pa_x[i];

        g_xxyyzzz_yyyzz_0[i] = g_yyzzz_yyyzz_0[i] * fbe_0 - g_yyzzz_yyyzz_1[i] * fz_be_0 + g_xyyzzz_yyyzz_1[i] * pa_x[i];

        g_xxyyzzz_yyzzz_0[i] = g_yyzzz_yyzzz_0[i] * fbe_0 - g_yyzzz_yyzzz_1[i] * fz_be_0 + g_xyyzzz_yyzzz_1[i] * pa_x[i];

        g_xxyyzzz_yzzzz_0[i] = g_yyzzz_yzzzz_0[i] * fbe_0 - g_yyzzz_yzzzz_1[i] * fz_be_0 + g_xyyzzz_yzzzz_1[i] * pa_x[i];

        g_xxyyzzz_zzzzz_0[i] = g_yyzzz_zzzzz_0[i] * fbe_0 - g_yyzzz_zzzzz_1[i] * fz_be_0 + g_xyyzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 399-420 components of targeted buffer : KH

    auto g_xxyzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 399);

    auto g_xxyzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 400);

    auto g_xxyzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 401);

    auto g_xxyzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 402);

    auto g_xxyzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 403);

    auto g_xxyzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 404);

    auto g_xxyzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 405);

    auto g_xxyzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 406);

    auto g_xxyzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 407);

    auto g_xxyzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 408);

    auto g_xxyzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 409);

    auto g_xxyzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 410);

    auto g_xxyzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 411);

    auto g_xxyzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 412);

    auto g_xxyzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 413);

    auto g_xxyzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 414);

    auto g_xxyzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 415);

    auto g_xxyzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 416);

    auto g_xxyzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 417);

    auto g_xxyzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 418);

    auto g_xxyzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 419);

    #pragma omp simd aligned(g_xxyzzzz_xxxxx_0, g_xxyzzzz_xxxxy_0, g_xxyzzzz_xxxxz_0, g_xxyzzzz_xxxyy_0, g_xxyzzzz_xxxyz_0, g_xxyzzzz_xxxzz_0, g_xxyzzzz_xxyyy_0, g_xxyzzzz_xxyyz_0, g_xxyzzzz_xxyzz_0, g_xxyzzzz_xxzzz_0, g_xxyzzzz_xyyyy_0, g_xxyzzzz_xyyyz_0, g_xxyzzzz_xyyzz_0, g_xxyzzzz_xyzzz_0, g_xxyzzzz_xzzzz_0, g_xxyzzzz_yyyyy_0, g_xxyzzzz_yyyyz_0, g_xxyzzzz_yyyzz_0, g_xxyzzzz_yyzzz_0, g_xxyzzzz_yzzzz_0, g_xxyzzzz_zzzzz_0, g_xxzzzz_xxxx_1, g_xxzzzz_xxxxx_1, g_xxzzzz_xxxxy_1, g_xxzzzz_xxxxz_1, g_xxzzzz_xxxy_1, g_xxzzzz_xxxyy_1, g_xxzzzz_xxxyz_1, g_xxzzzz_xxxz_1, g_xxzzzz_xxxzz_1, g_xxzzzz_xxyy_1, g_xxzzzz_xxyyy_1, g_xxzzzz_xxyyz_1, g_xxzzzz_xxyz_1, g_xxzzzz_xxyzz_1, g_xxzzzz_xxzz_1, g_xxzzzz_xxzzz_1, g_xxzzzz_xyyy_1, g_xxzzzz_xyyyy_1, g_xxzzzz_xyyyz_1, g_xxzzzz_xyyz_1, g_xxzzzz_xyyzz_1, g_xxzzzz_xyzz_1, g_xxzzzz_xyzzz_1, g_xxzzzz_xzzz_1, g_xxzzzz_xzzzz_1, g_xxzzzz_yyyy_1, g_xxzzzz_yyyyy_1, g_xxzzzz_yyyyz_1, g_xxzzzz_yyyz_1, g_xxzzzz_yyyzz_1, g_xxzzzz_yyzz_1, g_xxzzzz_yyzzz_1, g_xxzzzz_yzzz_1, g_xxzzzz_yzzzz_1, g_xxzzzz_zzzz_1, g_xxzzzz_zzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzzzz_xxxxx_0[i] = g_xxzzzz_xxxxx_1[i] * pa_y[i];

        g_xxyzzzz_xxxxy_0[i] = g_xxzzzz_xxxx_1[i] * fe_0 + g_xxzzzz_xxxxy_1[i] * pa_y[i];

        g_xxyzzzz_xxxxz_0[i] = g_xxzzzz_xxxxz_1[i] * pa_y[i];

        g_xxyzzzz_xxxyy_0[i] = 2.0 * g_xxzzzz_xxxy_1[i] * fe_0 + g_xxzzzz_xxxyy_1[i] * pa_y[i];

        g_xxyzzzz_xxxyz_0[i] = g_xxzzzz_xxxz_1[i] * fe_0 + g_xxzzzz_xxxyz_1[i] * pa_y[i];

        g_xxyzzzz_xxxzz_0[i] = g_xxzzzz_xxxzz_1[i] * pa_y[i];

        g_xxyzzzz_xxyyy_0[i] = 3.0 * g_xxzzzz_xxyy_1[i] * fe_0 + g_xxzzzz_xxyyy_1[i] * pa_y[i];

        g_xxyzzzz_xxyyz_0[i] = 2.0 * g_xxzzzz_xxyz_1[i] * fe_0 + g_xxzzzz_xxyyz_1[i] * pa_y[i];

        g_xxyzzzz_xxyzz_0[i] = g_xxzzzz_xxzz_1[i] * fe_0 + g_xxzzzz_xxyzz_1[i] * pa_y[i];

        g_xxyzzzz_xxzzz_0[i] = g_xxzzzz_xxzzz_1[i] * pa_y[i];

        g_xxyzzzz_xyyyy_0[i] = 4.0 * g_xxzzzz_xyyy_1[i] * fe_0 + g_xxzzzz_xyyyy_1[i] * pa_y[i];

        g_xxyzzzz_xyyyz_0[i] = 3.0 * g_xxzzzz_xyyz_1[i] * fe_0 + g_xxzzzz_xyyyz_1[i] * pa_y[i];

        g_xxyzzzz_xyyzz_0[i] = 2.0 * g_xxzzzz_xyzz_1[i] * fe_0 + g_xxzzzz_xyyzz_1[i] * pa_y[i];

        g_xxyzzzz_xyzzz_0[i] = g_xxzzzz_xzzz_1[i] * fe_0 + g_xxzzzz_xyzzz_1[i] * pa_y[i];

        g_xxyzzzz_xzzzz_0[i] = g_xxzzzz_xzzzz_1[i] * pa_y[i];

        g_xxyzzzz_yyyyy_0[i] = 5.0 * g_xxzzzz_yyyy_1[i] * fe_0 + g_xxzzzz_yyyyy_1[i] * pa_y[i];

        g_xxyzzzz_yyyyz_0[i] = 4.0 * g_xxzzzz_yyyz_1[i] * fe_0 + g_xxzzzz_yyyyz_1[i] * pa_y[i];

        g_xxyzzzz_yyyzz_0[i] = 3.0 * g_xxzzzz_yyzz_1[i] * fe_0 + g_xxzzzz_yyyzz_1[i] * pa_y[i];

        g_xxyzzzz_yyzzz_0[i] = 2.0 * g_xxzzzz_yzzz_1[i] * fe_0 + g_xxzzzz_yyzzz_1[i] * pa_y[i];

        g_xxyzzzz_yzzzz_0[i] = g_xxzzzz_zzzz_1[i] * fe_0 + g_xxzzzz_yzzzz_1[i] * pa_y[i];

        g_xxyzzzz_zzzzz_0[i] = g_xxzzzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 420-441 components of targeted buffer : KH

    auto g_xxzzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 420);

    auto g_xxzzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 421);

    auto g_xxzzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 422);

    auto g_xxzzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 423);

    auto g_xxzzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 424);

    auto g_xxzzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 425);

    auto g_xxzzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 426);

    auto g_xxzzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 427);

    auto g_xxzzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 428);

    auto g_xxzzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 429);

    auto g_xxzzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 430);

    auto g_xxzzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 431);

    auto g_xxzzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 432);

    auto g_xxzzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 433);

    auto g_xxzzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 434);

    auto g_xxzzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 435);

    auto g_xxzzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 436);

    auto g_xxzzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 437);

    auto g_xxzzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 438);

    auto g_xxzzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 439);

    auto g_xxzzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 440);

    #pragma omp simd aligned(g_xxzzz_xxxxx_0, g_xxzzz_xxxxx_1, g_xxzzz_xxxxy_0, g_xxzzz_xxxxy_1, g_xxzzz_xxxyy_0, g_xxzzz_xxxyy_1, g_xxzzz_xxyyy_0, g_xxzzz_xxyyy_1, g_xxzzz_xyyyy_0, g_xxzzz_xyyyy_1, g_xxzzzz_xxxxx_1, g_xxzzzz_xxxxy_1, g_xxzzzz_xxxyy_1, g_xxzzzz_xxyyy_1, g_xxzzzz_xyyyy_1, g_xxzzzzz_xxxxx_0, g_xxzzzzz_xxxxy_0, g_xxzzzzz_xxxxz_0, g_xxzzzzz_xxxyy_0, g_xxzzzzz_xxxyz_0, g_xxzzzzz_xxxzz_0, g_xxzzzzz_xxyyy_0, g_xxzzzzz_xxyyz_0, g_xxzzzzz_xxyzz_0, g_xxzzzzz_xxzzz_0, g_xxzzzzz_xyyyy_0, g_xxzzzzz_xyyyz_0, g_xxzzzzz_xyyzz_0, g_xxzzzzz_xyzzz_0, g_xxzzzzz_xzzzz_0, g_xxzzzzz_yyyyy_0, g_xxzzzzz_yyyyz_0, g_xxzzzzz_yyyzz_0, g_xxzzzzz_yyzzz_0, g_xxzzzzz_yzzzz_0, g_xxzzzzz_zzzzz_0, g_xzzzzz_xxxxz_1, g_xzzzzz_xxxyz_1, g_xzzzzz_xxxz_1, g_xzzzzz_xxxzz_1, g_xzzzzz_xxyyz_1, g_xzzzzz_xxyz_1, g_xzzzzz_xxyzz_1, g_xzzzzz_xxzz_1, g_xzzzzz_xxzzz_1, g_xzzzzz_xyyyz_1, g_xzzzzz_xyyz_1, g_xzzzzz_xyyzz_1, g_xzzzzz_xyzz_1, g_xzzzzz_xyzzz_1, g_xzzzzz_xzzz_1, g_xzzzzz_xzzzz_1, g_xzzzzz_yyyyy_1, g_xzzzzz_yyyyz_1, g_xzzzzz_yyyz_1, g_xzzzzz_yyyzz_1, g_xzzzzz_yyzz_1, g_xzzzzz_yyzzz_1, g_xzzzzz_yzzz_1, g_xzzzzz_yzzzz_1, g_xzzzzz_zzzz_1, g_xzzzzz_zzzzz_1, g_zzzzz_xxxxz_0, g_zzzzz_xxxxz_1, g_zzzzz_xxxyz_0, g_zzzzz_xxxyz_1, g_zzzzz_xxxzz_0, g_zzzzz_xxxzz_1, g_zzzzz_xxyyz_0, g_zzzzz_xxyyz_1, g_zzzzz_xxyzz_0, g_zzzzz_xxyzz_1, g_zzzzz_xxzzz_0, g_zzzzz_xxzzz_1, g_zzzzz_xyyyz_0, g_zzzzz_xyyyz_1, g_zzzzz_xyyzz_0, g_zzzzz_xyyzz_1, g_zzzzz_xyzzz_0, g_zzzzz_xyzzz_1, g_zzzzz_xzzzz_0, g_zzzzz_xzzzz_1, g_zzzzz_yyyyy_0, g_zzzzz_yyyyy_1, g_zzzzz_yyyyz_0, g_zzzzz_yyyyz_1, g_zzzzz_yyyzz_0, g_zzzzz_yyyzz_1, g_zzzzz_yyzzz_0, g_zzzzz_yyzzz_1, g_zzzzz_yzzzz_0, g_zzzzz_yzzzz_1, g_zzzzz_zzzzz_0, g_zzzzz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzzzz_xxxxx_0[i] = 4.0 * g_xxzzz_xxxxx_0[i] * fbe_0 - 4.0 * g_xxzzz_xxxxx_1[i] * fz_be_0 + g_xxzzzz_xxxxx_1[i] * pa_z[i];

        g_xxzzzzz_xxxxy_0[i] = 4.0 * g_xxzzz_xxxxy_0[i] * fbe_0 - 4.0 * g_xxzzz_xxxxy_1[i] * fz_be_0 + g_xxzzzz_xxxxy_1[i] * pa_z[i];

        g_xxzzzzz_xxxxz_0[i] = g_zzzzz_xxxxz_0[i] * fbe_0 - g_zzzzz_xxxxz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_xxxz_1[i] * fe_0 + g_xzzzzz_xxxxz_1[i] * pa_x[i];

        g_xxzzzzz_xxxyy_0[i] = 4.0 * g_xxzzz_xxxyy_0[i] * fbe_0 - 4.0 * g_xxzzz_xxxyy_1[i] * fz_be_0 + g_xxzzzz_xxxyy_1[i] * pa_z[i];

        g_xxzzzzz_xxxyz_0[i] = g_zzzzz_xxxyz_0[i] * fbe_0 - g_zzzzz_xxxyz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_xxyz_1[i] * fe_0 + g_xzzzzz_xxxyz_1[i] * pa_x[i];

        g_xxzzzzz_xxxzz_0[i] = g_zzzzz_xxxzz_0[i] * fbe_0 - g_zzzzz_xxxzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_xxzz_1[i] * fe_0 + g_xzzzzz_xxxzz_1[i] * pa_x[i];

        g_xxzzzzz_xxyyy_0[i] = 4.0 * g_xxzzz_xxyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_xxyyy_1[i] * fz_be_0 + g_xxzzzz_xxyyy_1[i] * pa_z[i];

        g_xxzzzzz_xxyyz_0[i] = g_zzzzz_xxyyz_0[i] * fbe_0 - g_zzzzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_xyyz_1[i] * fe_0 + g_xzzzzz_xxyyz_1[i] * pa_x[i];

        g_xxzzzzz_xxyzz_0[i] = g_zzzzz_xxyzz_0[i] * fbe_0 - g_zzzzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_xyzz_1[i] * fe_0 + g_xzzzzz_xxyzz_1[i] * pa_x[i];

        g_xxzzzzz_xxzzz_0[i] = g_zzzzz_xxzzz_0[i] * fbe_0 - g_zzzzz_xxzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_xzzz_1[i] * fe_0 + g_xzzzzz_xxzzz_1[i] * pa_x[i];

        g_xxzzzzz_xyyyy_0[i] = 4.0 * g_xxzzz_xyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_xyyyy_1[i] * fz_be_0 + g_xxzzzz_xyyyy_1[i] * pa_z[i];

        g_xxzzzzz_xyyyz_0[i] = g_zzzzz_xyyyz_0[i] * fbe_0 - g_zzzzz_xyyyz_1[i] * fz_be_0 + g_xzzzzz_yyyz_1[i] * fe_0 + g_xzzzzz_xyyyz_1[i] * pa_x[i];

        g_xxzzzzz_xyyzz_0[i] = g_zzzzz_xyyzz_0[i] * fbe_0 - g_zzzzz_xyyzz_1[i] * fz_be_0 + g_xzzzzz_yyzz_1[i] * fe_0 + g_xzzzzz_xyyzz_1[i] * pa_x[i];

        g_xxzzzzz_xyzzz_0[i] = g_zzzzz_xyzzz_0[i] * fbe_0 - g_zzzzz_xyzzz_1[i] * fz_be_0 + g_xzzzzz_yzzz_1[i] * fe_0 + g_xzzzzz_xyzzz_1[i] * pa_x[i];

        g_xxzzzzz_xzzzz_0[i] = g_zzzzz_xzzzz_0[i] * fbe_0 - g_zzzzz_xzzzz_1[i] * fz_be_0 + g_xzzzzz_zzzz_1[i] * fe_0 + g_xzzzzz_xzzzz_1[i] * pa_x[i];

        g_xxzzzzz_yyyyy_0[i] = g_zzzzz_yyyyy_0[i] * fbe_0 - g_zzzzz_yyyyy_1[i] * fz_be_0 + g_xzzzzz_yyyyy_1[i] * pa_x[i];

        g_xxzzzzz_yyyyz_0[i] = g_zzzzz_yyyyz_0[i] * fbe_0 - g_zzzzz_yyyyz_1[i] * fz_be_0 + g_xzzzzz_yyyyz_1[i] * pa_x[i];

        g_xxzzzzz_yyyzz_0[i] = g_zzzzz_yyyzz_0[i] * fbe_0 - g_zzzzz_yyyzz_1[i] * fz_be_0 + g_xzzzzz_yyyzz_1[i] * pa_x[i];

        g_xxzzzzz_yyzzz_0[i] = g_zzzzz_yyzzz_0[i] * fbe_0 - g_zzzzz_yyzzz_1[i] * fz_be_0 + g_xzzzzz_yyzzz_1[i] * pa_x[i];

        g_xxzzzzz_yzzzz_0[i] = g_zzzzz_yzzzz_0[i] * fbe_0 - g_zzzzz_yzzzz_1[i] * fz_be_0 + g_xzzzzz_yzzzz_1[i] * pa_x[i];

        g_xxzzzzz_zzzzz_0[i] = g_zzzzz_zzzzz_0[i] * fbe_0 - g_zzzzz_zzzzz_1[i] * fz_be_0 + g_xzzzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 441-462 components of targeted buffer : KH

    auto g_xyyyyyy_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 441);

    auto g_xyyyyyy_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 442);

    auto g_xyyyyyy_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 443);

    auto g_xyyyyyy_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 444);

    auto g_xyyyyyy_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 445);

    auto g_xyyyyyy_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 446);

    auto g_xyyyyyy_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 447);

    auto g_xyyyyyy_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 448);

    auto g_xyyyyyy_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 449);

    auto g_xyyyyyy_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 450);

    auto g_xyyyyyy_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 451);

    auto g_xyyyyyy_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 452);

    auto g_xyyyyyy_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 453);

    auto g_xyyyyyy_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 454);

    auto g_xyyyyyy_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 455);

    auto g_xyyyyyy_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 456);

    auto g_xyyyyyy_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 457);

    auto g_xyyyyyy_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 458);

    auto g_xyyyyyy_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 459);

    auto g_xyyyyyy_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 460);

    auto g_xyyyyyy_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 461);

    #pragma omp simd aligned(g_xyyyyyy_xxxxx_0, g_xyyyyyy_xxxxy_0, g_xyyyyyy_xxxxz_0, g_xyyyyyy_xxxyy_0, g_xyyyyyy_xxxyz_0, g_xyyyyyy_xxxzz_0, g_xyyyyyy_xxyyy_0, g_xyyyyyy_xxyyz_0, g_xyyyyyy_xxyzz_0, g_xyyyyyy_xxzzz_0, g_xyyyyyy_xyyyy_0, g_xyyyyyy_xyyyz_0, g_xyyyyyy_xyyzz_0, g_xyyyyyy_xyzzz_0, g_xyyyyyy_xzzzz_0, g_xyyyyyy_yyyyy_0, g_xyyyyyy_yyyyz_0, g_xyyyyyy_yyyzz_0, g_xyyyyyy_yyzzz_0, g_xyyyyyy_yzzzz_0, g_xyyyyyy_zzzzz_0, g_yyyyyy_xxxx_1, g_yyyyyy_xxxxx_1, g_yyyyyy_xxxxy_1, g_yyyyyy_xxxxz_1, g_yyyyyy_xxxy_1, g_yyyyyy_xxxyy_1, g_yyyyyy_xxxyz_1, g_yyyyyy_xxxz_1, g_yyyyyy_xxxzz_1, g_yyyyyy_xxyy_1, g_yyyyyy_xxyyy_1, g_yyyyyy_xxyyz_1, g_yyyyyy_xxyz_1, g_yyyyyy_xxyzz_1, g_yyyyyy_xxzz_1, g_yyyyyy_xxzzz_1, g_yyyyyy_xyyy_1, g_yyyyyy_xyyyy_1, g_yyyyyy_xyyyz_1, g_yyyyyy_xyyz_1, g_yyyyyy_xyyzz_1, g_yyyyyy_xyzz_1, g_yyyyyy_xyzzz_1, g_yyyyyy_xzzz_1, g_yyyyyy_xzzzz_1, g_yyyyyy_yyyy_1, g_yyyyyy_yyyyy_1, g_yyyyyy_yyyyz_1, g_yyyyyy_yyyz_1, g_yyyyyy_yyyzz_1, g_yyyyyy_yyzz_1, g_yyyyyy_yyzzz_1, g_yyyyyy_yzzz_1, g_yyyyyy_yzzzz_1, g_yyyyyy_zzzz_1, g_yyyyyy_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyyy_xxxxx_0[i] = 5.0 * g_yyyyyy_xxxx_1[i] * fe_0 + g_yyyyyy_xxxxx_1[i] * pa_x[i];

        g_xyyyyyy_xxxxy_0[i] = 4.0 * g_yyyyyy_xxxy_1[i] * fe_0 + g_yyyyyy_xxxxy_1[i] * pa_x[i];

        g_xyyyyyy_xxxxz_0[i] = 4.0 * g_yyyyyy_xxxz_1[i] * fe_0 + g_yyyyyy_xxxxz_1[i] * pa_x[i];

        g_xyyyyyy_xxxyy_0[i] = 3.0 * g_yyyyyy_xxyy_1[i] * fe_0 + g_yyyyyy_xxxyy_1[i] * pa_x[i];

        g_xyyyyyy_xxxyz_0[i] = 3.0 * g_yyyyyy_xxyz_1[i] * fe_0 + g_yyyyyy_xxxyz_1[i] * pa_x[i];

        g_xyyyyyy_xxxzz_0[i] = 3.0 * g_yyyyyy_xxzz_1[i] * fe_0 + g_yyyyyy_xxxzz_1[i] * pa_x[i];

        g_xyyyyyy_xxyyy_0[i] = 2.0 * g_yyyyyy_xyyy_1[i] * fe_0 + g_yyyyyy_xxyyy_1[i] * pa_x[i];

        g_xyyyyyy_xxyyz_0[i] = 2.0 * g_yyyyyy_xyyz_1[i] * fe_0 + g_yyyyyy_xxyyz_1[i] * pa_x[i];

        g_xyyyyyy_xxyzz_0[i] = 2.0 * g_yyyyyy_xyzz_1[i] * fe_0 + g_yyyyyy_xxyzz_1[i] * pa_x[i];

        g_xyyyyyy_xxzzz_0[i] = 2.0 * g_yyyyyy_xzzz_1[i] * fe_0 + g_yyyyyy_xxzzz_1[i] * pa_x[i];

        g_xyyyyyy_xyyyy_0[i] = g_yyyyyy_yyyy_1[i] * fe_0 + g_yyyyyy_xyyyy_1[i] * pa_x[i];

        g_xyyyyyy_xyyyz_0[i] = g_yyyyyy_yyyz_1[i] * fe_0 + g_yyyyyy_xyyyz_1[i] * pa_x[i];

        g_xyyyyyy_xyyzz_0[i] = g_yyyyyy_yyzz_1[i] * fe_0 + g_yyyyyy_xyyzz_1[i] * pa_x[i];

        g_xyyyyyy_xyzzz_0[i] = g_yyyyyy_yzzz_1[i] * fe_0 + g_yyyyyy_xyzzz_1[i] * pa_x[i];

        g_xyyyyyy_xzzzz_0[i] = g_yyyyyy_zzzz_1[i] * fe_0 + g_yyyyyy_xzzzz_1[i] * pa_x[i];

        g_xyyyyyy_yyyyy_0[i] = g_yyyyyy_yyyyy_1[i] * pa_x[i];

        g_xyyyyyy_yyyyz_0[i] = g_yyyyyy_yyyyz_1[i] * pa_x[i];

        g_xyyyyyy_yyyzz_0[i] = g_yyyyyy_yyyzz_1[i] * pa_x[i];

        g_xyyyyyy_yyzzz_0[i] = g_yyyyyy_yyzzz_1[i] * pa_x[i];

        g_xyyyyyy_yzzzz_0[i] = g_yyyyyy_yzzzz_1[i] * pa_x[i];

        g_xyyyyyy_zzzzz_0[i] = g_yyyyyy_zzzzz_1[i] * pa_x[i];
    }

    // Set up 462-483 components of targeted buffer : KH

    auto g_xyyyyyz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 462);

    auto g_xyyyyyz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 463);

    auto g_xyyyyyz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 464);

    auto g_xyyyyyz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 465);

    auto g_xyyyyyz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 466);

    auto g_xyyyyyz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 467);

    auto g_xyyyyyz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 468);

    auto g_xyyyyyz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 469);

    auto g_xyyyyyz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 470);

    auto g_xyyyyyz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 471);

    auto g_xyyyyyz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 472);

    auto g_xyyyyyz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 473);

    auto g_xyyyyyz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 474);

    auto g_xyyyyyz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 475);

    auto g_xyyyyyz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 476);

    auto g_xyyyyyz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 477);

    auto g_xyyyyyz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 478);

    auto g_xyyyyyz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 479);

    auto g_xyyyyyz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 480);

    auto g_xyyyyyz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 481);

    auto g_xyyyyyz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 482);

    #pragma omp simd aligned(g_xyyyyy_xxxxx_1, g_xyyyyy_xxxxy_1, g_xyyyyy_xxxyy_1, g_xyyyyy_xxyyy_1, g_xyyyyy_xyyyy_1, g_xyyyyyz_xxxxx_0, g_xyyyyyz_xxxxy_0, g_xyyyyyz_xxxxz_0, g_xyyyyyz_xxxyy_0, g_xyyyyyz_xxxyz_0, g_xyyyyyz_xxxzz_0, g_xyyyyyz_xxyyy_0, g_xyyyyyz_xxyyz_0, g_xyyyyyz_xxyzz_0, g_xyyyyyz_xxzzz_0, g_xyyyyyz_xyyyy_0, g_xyyyyyz_xyyyz_0, g_xyyyyyz_xyyzz_0, g_xyyyyyz_xyzzz_0, g_xyyyyyz_xzzzz_0, g_xyyyyyz_yyyyy_0, g_xyyyyyz_yyyyz_0, g_xyyyyyz_yyyzz_0, g_xyyyyyz_yyzzz_0, g_xyyyyyz_yzzzz_0, g_xyyyyyz_zzzzz_0, g_yyyyyz_xxxxz_1, g_yyyyyz_xxxyz_1, g_yyyyyz_xxxz_1, g_yyyyyz_xxxzz_1, g_yyyyyz_xxyyz_1, g_yyyyyz_xxyz_1, g_yyyyyz_xxyzz_1, g_yyyyyz_xxzz_1, g_yyyyyz_xxzzz_1, g_yyyyyz_xyyyz_1, g_yyyyyz_xyyz_1, g_yyyyyz_xyyzz_1, g_yyyyyz_xyzz_1, g_yyyyyz_xyzzz_1, g_yyyyyz_xzzz_1, g_yyyyyz_xzzzz_1, g_yyyyyz_yyyyy_1, g_yyyyyz_yyyyz_1, g_yyyyyz_yyyz_1, g_yyyyyz_yyyzz_1, g_yyyyyz_yyzz_1, g_yyyyyz_yyzzz_1, g_yyyyyz_yzzz_1, g_yyyyyz_yzzzz_1, g_yyyyyz_zzzz_1, g_yyyyyz_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyyz_xxxxx_0[i] = g_xyyyyy_xxxxx_1[i] * pa_z[i];

        g_xyyyyyz_xxxxy_0[i] = g_xyyyyy_xxxxy_1[i] * pa_z[i];

        g_xyyyyyz_xxxxz_0[i] = 4.0 * g_yyyyyz_xxxz_1[i] * fe_0 + g_yyyyyz_xxxxz_1[i] * pa_x[i];

        g_xyyyyyz_xxxyy_0[i] = g_xyyyyy_xxxyy_1[i] * pa_z[i];

        g_xyyyyyz_xxxyz_0[i] = 3.0 * g_yyyyyz_xxyz_1[i] * fe_0 + g_yyyyyz_xxxyz_1[i] * pa_x[i];

        g_xyyyyyz_xxxzz_0[i] = 3.0 * g_yyyyyz_xxzz_1[i] * fe_0 + g_yyyyyz_xxxzz_1[i] * pa_x[i];

        g_xyyyyyz_xxyyy_0[i] = g_xyyyyy_xxyyy_1[i] * pa_z[i];

        g_xyyyyyz_xxyyz_0[i] = 2.0 * g_yyyyyz_xyyz_1[i] * fe_0 + g_yyyyyz_xxyyz_1[i] * pa_x[i];

        g_xyyyyyz_xxyzz_0[i] = 2.0 * g_yyyyyz_xyzz_1[i] * fe_0 + g_yyyyyz_xxyzz_1[i] * pa_x[i];

        g_xyyyyyz_xxzzz_0[i] = 2.0 * g_yyyyyz_xzzz_1[i] * fe_0 + g_yyyyyz_xxzzz_1[i] * pa_x[i];

        g_xyyyyyz_xyyyy_0[i] = g_xyyyyy_xyyyy_1[i] * pa_z[i];

        g_xyyyyyz_xyyyz_0[i] = g_yyyyyz_yyyz_1[i] * fe_0 + g_yyyyyz_xyyyz_1[i] * pa_x[i];

        g_xyyyyyz_xyyzz_0[i] = g_yyyyyz_yyzz_1[i] * fe_0 + g_yyyyyz_xyyzz_1[i] * pa_x[i];

        g_xyyyyyz_xyzzz_0[i] = g_yyyyyz_yzzz_1[i] * fe_0 + g_yyyyyz_xyzzz_1[i] * pa_x[i];

        g_xyyyyyz_xzzzz_0[i] = g_yyyyyz_zzzz_1[i] * fe_0 + g_yyyyyz_xzzzz_1[i] * pa_x[i];

        g_xyyyyyz_yyyyy_0[i] = g_yyyyyz_yyyyy_1[i] * pa_x[i];

        g_xyyyyyz_yyyyz_0[i] = g_yyyyyz_yyyyz_1[i] * pa_x[i];

        g_xyyyyyz_yyyzz_0[i] = g_yyyyyz_yyyzz_1[i] * pa_x[i];

        g_xyyyyyz_yyzzz_0[i] = g_yyyyyz_yyzzz_1[i] * pa_x[i];

        g_xyyyyyz_yzzzz_0[i] = g_yyyyyz_yzzzz_1[i] * pa_x[i];

        g_xyyyyyz_zzzzz_0[i] = g_yyyyyz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 483-504 components of targeted buffer : KH

    auto g_xyyyyzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 483);

    auto g_xyyyyzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 484);

    auto g_xyyyyzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 485);

    auto g_xyyyyzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 486);

    auto g_xyyyyzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 487);

    auto g_xyyyyzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 488);

    auto g_xyyyyzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 489);

    auto g_xyyyyzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 490);

    auto g_xyyyyzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 491);

    auto g_xyyyyzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 492);

    auto g_xyyyyzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 493);

    auto g_xyyyyzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 494);

    auto g_xyyyyzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 495);

    auto g_xyyyyzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 496);

    auto g_xyyyyzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 497);

    auto g_xyyyyzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 498);

    auto g_xyyyyzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 499);

    auto g_xyyyyzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 500);

    auto g_xyyyyzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 501);

    auto g_xyyyyzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 502);

    auto g_xyyyyzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 503);

    #pragma omp simd aligned(g_xyyyyzz_xxxxx_0, g_xyyyyzz_xxxxy_0, g_xyyyyzz_xxxxz_0, g_xyyyyzz_xxxyy_0, g_xyyyyzz_xxxyz_0, g_xyyyyzz_xxxzz_0, g_xyyyyzz_xxyyy_0, g_xyyyyzz_xxyyz_0, g_xyyyyzz_xxyzz_0, g_xyyyyzz_xxzzz_0, g_xyyyyzz_xyyyy_0, g_xyyyyzz_xyyyz_0, g_xyyyyzz_xyyzz_0, g_xyyyyzz_xyzzz_0, g_xyyyyzz_xzzzz_0, g_xyyyyzz_yyyyy_0, g_xyyyyzz_yyyyz_0, g_xyyyyzz_yyyzz_0, g_xyyyyzz_yyzzz_0, g_xyyyyzz_yzzzz_0, g_xyyyyzz_zzzzz_0, g_yyyyzz_xxxx_1, g_yyyyzz_xxxxx_1, g_yyyyzz_xxxxy_1, g_yyyyzz_xxxxz_1, g_yyyyzz_xxxy_1, g_yyyyzz_xxxyy_1, g_yyyyzz_xxxyz_1, g_yyyyzz_xxxz_1, g_yyyyzz_xxxzz_1, g_yyyyzz_xxyy_1, g_yyyyzz_xxyyy_1, g_yyyyzz_xxyyz_1, g_yyyyzz_xxyz_1, g_yyyyzz_xxyzz_1, g_yyyyzz_xxzz_1, g_yyyyzz_xxzzz_1, g_yyyyzz_xyyy_1, g_yyyyzz_xyyyy_1, g_yyyyzz_xyyyz_1, g_yyyyzz_xyyz_1, g_yyyyzz_xyyzz_1, g_yyyyzz_xyzz_1, g_yyyyzz_xyzzz_1, g_yyyyzz_xzzz_1, g_yyyyzz_xzzzz_1, g_yyyyzz_yyyy_1, g_yyyyzz_yyyyy_1, g_yyyyzz_yyyyz_1, g_yyyyzz_yyyz_1, g_yyyyzz_yyyzz_1, g_yyyyzz_yyzz_1, g_yyyyzz_yyzzz_1, g_yyyyzz_yzzz_1, g_yyyyzz_yzzzz_1, g_yyyyzz_zzzz_1, g_yyyyzz_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyzz_xxxxx_0[i] = 5.0 * g_yyyyzz_xxxx_1[i] * fe_0 + g_yyyyzz_xxxxx_1[i] * pa_x[i];

        g_xyyyyzz_xxxxy_0[i] = 4.0 * g_yyyyzz_xxxy_1[i] * fe_0 + g_yyyyzz_xxxxy_1[i] * pa_x[i];

        g_xyyyyzz_xxxxz_0[i] = 4.0 * g_yyyyzz_xxxz_1[i] * fe_0 + g_yyyyzz_xxxxz_1[i] * pa_x[i];

        g_xyyyyzz_xxxyy_0[i] = 3.0 * g_yyyyzz_xxyy_1[i] * fe_0 + g_yyyyzz_xxxyy_1[i] * pa_x[i];

        g_xyyyyzz_xxxyz_0[i] = 3.0 * g_yyyyzz_xxyz_1[i] * fe_0 + g_yyyyzz_xxxyz_1[i] * pa_x[i];

        g_xyyyyzz_xxxzz_0[i] = 3.0 * g_yyyyzz_xxzz_1[i] * fe_0 + g_yyyyzz_xxxzz_1[i] * pa_x[i];

        g_xyyyyzz_xxyyy_0[i] = 2.0 * g_yyyyzz_xyyy_1[i] * fe_0 + g_yyyyzz_xxyyy_1[i] * pa_x[i];

        g_xyyyyzz_xxyyz_0[i] = 2.0 * g_yyyyzz_xyyz_1[i] * fe_0 + g_yyyyzz_xxyyz_1[i] * pa_x[i];

        g_xyyyyzz_xxyzz_0[i] = 2.0 * g_yyyyzz_xyzz_1[i] * fe_0 + g_yyyyzz_xxyzz_1[i] * pa_x[i];

        g_xyyyyzz_xxzzz_0[i] = 2.0 * g_yyyyzz_xzzz_1[i] * fe_0 + g_yyyyzz_xxzzz_1[i] * pa_x[i];

        g_xyyyyzz_xyyyy_0[i] = g_yyyyzz_yyyy_1[i] * fe_0 + g_yyyyzz_xyyyy_1[i] * pa_x[i];

        g_xyyyyzz_xyyyz_0[i] = g_yyyyzz_yyyz_1[i] * fe_0 + g_yyyyzz_xyyyz_1[i] * pa_x[i];

        g_xyyyyzz_xyyzz_0[i] = g_yyyyzz_yyzz_1[i] * fe_0 + g_yyyyzz_xyyzz_1[i] * pa_x[i];

        g_xyyyyzz_xyzzz_0[i] = g_yyyyzz_yzzz_1[i] * fe_0 + g_yyyyzz_xyzzz_1[i] * pa_x[i];

        g_xyyyyzz_xzzzz_0[i] = g_yyyyzz_zzzz_1[i] * fe_0 + g_yyyyzz_xzzzz_1[i] * pa_x[i];

        g_xyyyyzz_yyyyy_0[i] = g_yyyyzz_yyyyy_1[i] * pa_x[i];

        g_xyyyyzz_yyyyz_0[i] = g_yyyyzz_yyyyz_1[i] * pa_x[i];

        g_xyyyyzz_yyyzz_0[i] = g_yyyyzz_yyyzz_1[i] * pa_x[i];

        g_xyyyyzz_yyzzz_0[i] = g_yyyyzz_yyzzz_1[i] * pa_x[i];

        g_xyyyyzz_yzzzz_0[i] = g_yyyyzz_yzzzz_1[i] * pa_x[i];

        g_xyyyyzz_zzzzz_0[i] = g_yyyyzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 504-525 components of targeted buffer : KH

    auto g_xyyyzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 504);

    auto g_xyyyzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 505);

    auto g_xyyyzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 506);

    auto g_xyyyzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 507);

    auto g_xyyyzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 508);

    auto g_xyyyzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 509);

    auto g_xyyyzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 510);

    auto g_xyyyzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 511);

    auto g_xyyyzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 512);

    auto g_xyyyzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 513);

    auto g_xyyyzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 514);

    auto g_xyyyzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 515);

    auto g_xyyyzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 516);

    auto g_xyyyzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 517);

    auto g_xyyyzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 518);

    auto g_xyyyzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 519);

    auto g_xyyyzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 520);

    auto g_xyyyzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 521);

    auto g_xyyyzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 522);

    auto g_xyyyzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 523);

    auto g_xyyyzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 524);

    #pragma omp simd aligned(g_xyyyzzz_xxxxx_0, g_xyyyzzz_xxxxy_0, g_xyyyzzz_xxxxz_0, g_xyyyzzz_xxxyy_0, g_xyyyzzz_xxxyz_0, g_xyyyzzz_xxxzz_0, g_xyyyzzz_xxyyy_0, g_xyyyzzz_xxyyz_0, g_xyyyzzz_xxyzz_0, g_xyyyzzz_xxzzz_0, g_xyyyzzz_xyyyy_0, g_xyyyzzz_xyyyz_0, g_xyyyzzz_xyyzz_0, g_xyyyzzz_xyzzz_0, g_xyyyzzz_xzzzz_0, g_xyyyzzz_yyyyy_0, g_xyyyzzz_yyyyz_0, g_xyyyzzz_yyyzz_0, g_xyyyzzz_yyzzz_0, g_xyyyzzz_yzzzz_0, g_xyyyzzz_zzzzz_0, g_yyyzzz_xxxx_1, g_yyyzzz_xxxxx_1, g_yyyzzz_xxxxy_1, g_yyyzzz_xxxxz_1, g_yyyzzz_xxxy_1, g_yyyzzz_xxxyy_1, g_yyyzzz_xxxyz_1, g_yyyzzz_xxxz_1, g_yyyzzz_xxxzz_1, g_yyyzzz_xxyy_1, g_yyyzzz_xxyyy_1, g_yyyzzz_xxyyz_1, g_yyyzzz_xxyz_1, g_yyyzzz_xxyzz_1, g_yyyzzz_xxzz_1, g_yyyzzz_xxzzz_1, g_yyyzzz_xyyy_1, g_yyyzzz_xyyyy_1, g_yyyzzz_xyyyz_1, g_yyyzzz_xyyz_1, g_yyyzzz_xyyzz_1, g_yyyzzz_xyzz_1, g_yyyzzz_xyzzz_1, g_yyyzzz_xzzz_1, g_yyyzzz_xzzzz_1, g_yyyzzz_yyyy_1, g_yyyzzz_yyyyy_1, g_yyyzzz_yyyyz_1, g_yyyzzz_yyyz_1, g_yyyzzz_yyyzz_1, g_yyyzzz_yyzz_1, g_yyyzzz_yyzzz_1, g_yyyzzz_yzzz_1, g_yyyzzz_yzzzz_1, g_yyyzzz_zzzz_1, g_yyyzzz_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyzzz_xxxxx_0[i] = 5.0 * g_yyyzzz_xxxx_1[i] * fe_0 + g_yyyzzz_xxxxx_1[i] * pa_x[i];

        g_xyyyzzz_xxxxy_0[i] = 4.0 * g_yyyzzz_xxxy_1[i] * fe_0 + g_yyyzzz_xxxxy_1[i] * pa_x[i];

        g_xyyyzzz_xxxxz_0[i] = 4.0 * g_yyyzzz_xxxz_1[i] * fe_0 + g_yyyzzz_xxxxz_1[i] * pa_x[i];

        g_xyyyzzz_xxxyy_0[i] = 3.0 * g_yyyzzz_xxyy_1[i] * fe_0 + g_yyyzzz_xxxyy_1[i] * pa_x[i];

        g_xyyyzzz_xxxyz_0[i] = 3.0 * g_yyyzzz_xxyz_1[i] * fe_0 + g_yyyzzz_xxxyz_1[i] * pa_x[i];

        g_xyyyzzz_xxxzz_0[i] = 3.0 * g_yyyzzz_xxzz_1[i] * fe_0 + g_yyyzzz_xxxzz_1[i] * pa_x[i];

        g_xyyyzzz_xxyyy_0[i] = 2.0 * g_yyyzzz_xyyy_1[i] * fe_0 + g_yyyzzz_xxyyy_1[i] * pa_x[i];

        g_xyyyzzz_xxyyz_0[i] = 2.0 * g_yyyzzz_xyyz_1[i] * fe_0 + g_yyyzzz_xxyyz_1[i] * pa_x[i];

        g_xyyyzzz_xxyzz_0[i] = 2.0 * g_yyyzzz_xyzz_1[i] * fe_0 + g_yyyzzz_xxyzz_1[i] * pa_x[i];

        g_xyyyzzz_xxzzz_0[i] = 2.0 * g_yyyzzz_xzzz_1[i] * fe_0 + g_yyyzzz_xxzzz_1[i] * pa_x[i];

        g_xyyyzzz_xyyyy_0[i] = g_yyyzzz_yyyy_1[i] * fe_0 + g_yyyzzz_xyyyy_1[i] * pa_x[i];

        g_xyyyzzz_xyyyz_0[i] = g_yyyzzz_yyyz_1[i] * fe_0 + g_yyyzzz_xyyyz_1[i] * pa_x[i];

        g_xyyyzzz_xyyzz_0[i] = g_yyyzzz_yyzz_1[i] * fe_0 + g_yyyzzz_xyyzz_1[i] * pa_x[i];

        g_xyyyzzz_xyzzz_0[i] = g_yyyzzz_yzzz_1[i] * fe_0 + g_yyyzzz_xyzzz_1[i] * pa_x[i];

        g_xyyyzzz_xzzzz_0[i] = g_yyyzzz_zzzz_1[i] * fe_0 + g_yyyzzz_xzzzz_1[i] * pa_x[i];

        g_xyyyzzz_yyyyy_0[i] = g_yyyzzz_yyyyy_1[i] * pa_x[i];

        g_xyyyzzz_yyyyz_0[i] = g_yyyzzz_yyyyz_1[i] * pa_x[i];

        g_xyyyzzz_yyyzz_0[i] = g_yyyzzz_yyyzz_1[i] * pa_x[i];

        g_xyyyzzz_yyzzz_0[i] = g_yyyzzz_yyzzz_1[i] * pa_x[i];

        g_xyyyzzz_yzzzz_0[i] = g_yyyzzz_yzzzz_1[i] * pa_x[i];

        g_xyyyzzz_zzzzz_0[i] = g_yyyzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 525-546 components of targeted buffer : KH

    auto g_xyyzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 525);

    auto g_xyyzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 526);

    auto g_xyyzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 527);

    auto g_xyyzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 528);

    auto g_xyyzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 529);

    auto g_xyyzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 530);

    auto g_xyyzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 531);

    auto g_xyyzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 532);

    auto g_xyyzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 533);

    auto g_xyyzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 534);

    auto g_xyyzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 535);

    auto g_xyyzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 536);

    auto g_xyyzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 537);

    auto g_xyyzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 538);

    auto g_xyyzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 539);

    auto g_xyyzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 540);

    auto g_xyyzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 541);

    auto g_xyyzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 542);

    auto g_xyyzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 543);

    auto g_xyyzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 544);

    auto g_xyyzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 545);

    #pragma omp simd aligned(g_xyyzzzz_xxxxx_0, g_xyyzzzz_xxxxy_0, g_xyyzzzz_xxxxz_0, g_xyyzzzz_xxxyy_0, g_xyyzzzz_xxxyz_0, g_xyyzzzz_xxxzz_0, g_xyyzzzz_xxyyy_0, g_xyyzzzz_xxyyz_0, g_xyyzzzz_xxyzz_0, g_xyyzzzz_xxzzz_0, g_xyyzzzz_xyyyy_0, g_xyyzzzz_xyyyz_0, g_xyyzzzz_xyyzz_0, g_xyyzzzz_xyzzz_0, g_xyyzzzz_xzzzz_0, g_xyyzzzz_yyyyy_0, g_xyyzzzz_yyyyz_0, g_xyyzzzz_yyyzz_0, g_xyyzzzz_yyzzz_0, g_xyyzzzz_yzzzz_0, g_xyyzzzz_zzzzz_0, g_yyzzzz_xxxx_1, g_yyzzzz_xxxxx_1, g_yyzzzz_xxxxy_1, g_yyzzzz_xxxxz_1, g_yyzzzz_xxxy_1, g_yyzzzz_xxxyy_1, g_yyzzzz_xxxyz_1, g_yyzzzz_xxxz_1, g_yyzzzz_xxxzz_1, g_yyzzzz_xxyy_1, g_yyzzzz_xxyyy_1, g_yyzzzz_xxyyz_1, g_yyzzzz_xxyz_1, g_yyzzzz_xxyzz_1, g_yyzzzz_xxzz_1, g_yyzzzz_xxzzz_1, g_yyzzzz_xyyy_1, g_yyzzzz_xyyyy_1, g_yyzzzz_xyyyz_1, g_yyzzzz_xyyz_1, g_yyzzzz_xyyzz_1, g_yyzzzz_xyzz_1, g_yyzzzz_xyzzz_1, g_yyzzzz_xzzz_1, g_yyzzzz_xzzzz_1, g_yyzzzz_yyyy_1, g_yyzzzz_yyyyy_1, g_yyzzzz_yyyyz_1, g_yyzzzz_yyyz_1, g_yyzzzz_yyyzz_1, g_yyzzzz_yyzz_1, g_yyzzzz_yyzzz_1, g_yyzzzz_yzzz_1, g_yyzzzz_yzzzz_1, g_yyzzzz_zzzz_1, g_yyzzzz_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzzzz_xxxxx_0[i] = 5.0 * g_yyzzzz_xxxx_1[i] * fe_0 + g_yyzzzz_xxxxx_1[i] * pa_x[i];

        g_xyyzzzz_xxxxy_0[i] = 4.0 * g_yyzzzz_xxxy_1[i] * fe_0 + g_yyzzzz_xxxxy_1[i] * pa_x[i];

        g_xyyzzzz_xxxxz_0[i] = 4.0 * g_yyzzzz_xxxz_1[i] * fe_0 + g_yyzzzz_xxxxz_1[i] * pa_x[i];

        g_xyyzzzz_xxxyy_0[i] = 3.0 * g_yyzzzz_xxyy_1[i] * fe_0 + g_yyzzzz_xxxyy_1[i] * pa_x[i];

        g_xyyzzzz_xxxyz_0[i] = 3.0 * g_yyzzzz_xxyz_1[i] * fe_0 + g_yyzzzz_xxxyz_1[i] * pa_x[i];

        g_xyyzzzz_xxxzz_0[i] = 3.0 * g_yyzzzz_xxzz_1[i] * fe_0 + g_yyzzzz_xxxzz_1[i] * pa_x[i];

        g_xyyzzzz_xxyyy_0[i] = 2.0 * g_yyzzzz_xyyy_1[i] * fe_0 + g_yyzzzz_xxyyy_1[i] * pa_x[i];

        g_xyyzzzz_xxyyz_0[i] = 2.0 * g_yyzzzz_xyyz_1[i] * fe_0 + g_yyzzzz_xxyyz_1[i] * pa_x[i];

        g_xyyzzzz_xxyzz_0[i] = 2.0 * g_yyzzzz_xyzz_1[i] * fe_0 + g_yyzzzz_xxyzz_1[i] * pa_x[i];

        g_xyyzzzz_xxzzz_0[i] = 2.0 * g_yyzzzz_xzzz_1[i] * fe_0 + g_yyzzzz_xxzzz_1[i] * pa_x[i];

        g_xyyzzzz_xyyyy_0[i] = g_yyzzzz_yyyy_1[i] * fe_0 + g_yyzzzz_xyyyy_1[i] * pa_x[i];

        g_xyyzzzz_xyyyz_0[i] = g_yyzzzz_yyyz_1[i] * fe_0 + g_yyzzzz_xyyyz_1[i] * pa_x[i];

        g_xyyzzzz_xyyzz_0[i] = g_yyzzzz_yyzz_1[i] * fe_0 + g_yyzzzz_xyyzz_1[i] * pa_x[i];

        g_xyyzzzz_xyzzz_0[i] = g_yyzzzz_yzzz_1[i] * fe_0 + g_yyzzzz_xyzzz_1[i] * pa_x[i];

        g_xyyzzzz_xzzzz_0[i] = g_yyzzzz_zzzz_1[i] * fe_0 + g_yyzzzz_xzzzz_1[i] * pa_x[i];

        g_xyyzzzz_yyyyy_0[i] = g_yyzzzz_yyyyy_1[i] * pa_x[i];

        g_xyyzzzz_yyyyz_0[i] = g_yyzzzz_yyyyz_1[i] * pa_x[i];

        g_xyyzzzz_yyyzz_0[i] = g_yyzzzz_yyyzz_1[i] * pa_x[i];

        g_xyyzzzz_yyzzz_0[i] = g_yyzzzz_yyzzz_1[i] * pa_x[i];

        g_xyyzzzz_yzzzz_0[i] = g_yyzzzz_yzzzz_1[i] * pa_x[i];

        g_xyyzzzz_zzzzz_0[i] = g_yyzzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 546-567 components of targeted buffer : KH

    auto g_xyzzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 546);

    auto g_xyzzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 547);

    auto g_xyzzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 548);

    auto g_xyzzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 549);

    auto g_xyzzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 550);

    auto g_xyzzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 551);

    auto g_xyzzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 552);

    auto g_xyzzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 553);

    auto g_xyzzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 554);

    auto g_xyzzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 555);

    auto g_xyzzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 556);

    auto g_xyzzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 557);

    auto g_xyzzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 558);

    auto g_xyzzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 559);

    auto g_xyzzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 560);

    auto g_xyzzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 561);

    auto g_xyzzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 562);

    auto g_xyzzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 563);

    auto g_xyzzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 564);

    auto g_xyzzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 565);

    auto g_xyzzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 566);

    #pragma omp simd aligned(g_xyzzzzz_xxxxx_0, g_xyzzzzz_xxxxy_0, g_xyzzzzz_xxxxz_0, g_xyzzzzz_xxxyy_0, g_xyzzzzz_xxxyz_0, g_xyzzzzz_xxxzz_0, g_xyzzzzz_xxyyy_0, g_xyzzzzz_xxyyz_0, g_xyzzzzz_xxyzz_0, g_xyzzzzz_xxzzz_0, g_xyzzzzz_xyyyy_0, g_xyzzzzz_xyyyz_0, g_xyzzzzz_xyyzz_0, g_xyzzzzz_xyzzz_0, g_xyzzzzz_xzzzz_0, g_xyzzzzz_yyyyy_0, g_xyzzzzz_yyyyz_0, g_xyzzzzz_yyyzz_0, g_xyzzzzz_yyzzz_0, g_xyzzzzz_yzzzz_0, g_xyzzzzz_zzzzz_0, g_xzzzzz_xxxxx_1, g_xzzzzz_xxxxz_1, g_xzzzzz_xxxzz_1, g_xzzzzz_xxzzz_1, g_xzzzzz_xzzzz_1, g_yzzzzz_xxxxy_1, g_yzzzzz_xxxy_1, g_yzzzzz_xxxyy_1, g_yzzzzz_xxxyz_1, g_yzzzzz_xxyy_1, g_yzzzzz_xxyyy_1, g_yzzzzz_xxyyz_1, g_yzzzzz_xxyz_1, g_yzzzzz_xxyzz_1, g_yzzzzz_xyyy_1, g_yzzzzz_xyyyy_1, g_yzzzzz_xyyyz_1, g_yzzzzz_xyyz_1, g_yzzzzz_xyyzz_1, g_yzzzzz_xyzz_1, g_yzzzzz_xyzzz_1, g_yzzzzz_yyyy_1, g_yzzzzz_yyyyy_1, g_yzzzzz_yyyyz_1, g_yzzzzz_yyyz_1, g_yzzzzz_yyyzz_1, g_yzzzzz_yyzz_1, g_yzzzzz_yyzzz_1, g_yzzzzz_yzzz_1, g_yzzzzz_yzzzz_1, g_yzzzzz_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzzzz_xxxxx_0[i] = g_xzzzzz_xxxxx_1[i] * pa_y[i];

        g_xyzzzzz_xxxxy_0[i] = 4.0 * g_yzzzzz_xxxy_1[i] * fe_0 + g_yzzzzz_xxxxy_1[i] * pa_x[i];

        g_xyzzzzz_xxxxz_0[i] = g_xzzzzz_xxxxz_1[i] * pa_y[i];

        g_xyzzzzz_xxxyy_0[i] = 3.0 * g_yzzzzz_xxyy_1[i] * fe_0 + g_yzzzzz_xxxyy_1[i] * pa_x[i];

        g_xyzzzzz_xxxyz_0[i] = 3.0 * g_yzzzzz_xxyz_1[i] * fe_0 + g_yzzzzz_xxxyz_1[i] * pa_x[i];

        g_xyzzzzz_xxxzz_0[i] = g_xzzzzz_xxxzz_1[i] * pa_y[i];

        g_xyzzzzz_xxyyy_0[i] = 2.0 * g_yzzzzz_xyyy_1[i] * fe_0 + g_yzzzzz_xxyyy_1[i] * pa_x[i];

        g_xyzzzzz_xxyyz_0[i] = 2.0 * g_yzzzzz_xyyz_1[i] * fe_0 + g_yzzzzz_xxyyz_1[i] * pa_x[i];

        g_xyzzzzz_xxyzz_0[i] = 2.0 * g_yzzzzz_xyzz_1[i] * fe_0 + g_yzzzzz_xxyzz_1[i] * pa_x[i];

        g_xyzzzzz_xxzzz_0[i] = g_xzzzzz_xxzzz_1[i] * pa_y[i];

        g_xyzzzzz_xyyyy_0[i] = g_yzzzzz_yyyy_1[i] * fe_0 + g_yzzzzz_xyyyy_1[i] * pa_x[i];

        g_xyzzzzz_xyyyz_0[i] = g_yzzzzz_yyyz_1[i] * fe_0 + g_yzzzzz_xyyyz_1[i] * pa_x[i];

        g_xyzzzzz_xyyzz_0[i] = g_yzzzzz_yyzz_1[i] * fe_0 + g_yzzzzz_xyyzz_1[i] * pa_x[i];

        g_xyzzzzz_xyzzz_0[i] = g_yzzzzz_yzzz_1[i] * fe_0 + g_yzzzzz_xyzzz_1[i] * pa_x[i];

        g_xyzzzzz_xzzzz_0[i] = g_xzzzzz_xzzzz_1[i] * pa_y[i];

        g_xyzzzzz_yyyyy_0[i] = g_yzzzzz_yyyyy_1[i] * pa_x[i];

        g_xyzzzzz_yyyyz_0[i] = g_yzzzzz_yyyyz_1[i] * pa_x[i];

        g_xyzzzzz_yyyzz_0[i] = g_yzzzzz_yyyzz_1[i] * pa_x[i];

        g_xyzzzzz_yyzzz_0[i] = g_yzzzzz_yyzzz_1[i] * pa_x[i];

        g_xyzzzzz_yzzzz_0[i] = g_yzzzzz_yzzzz_1[i] * pa_x[i];

        g_xyzzzzz_zzzzz_0[i] = g_yzzzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 567-588 components of targeted buffer : KH

    auto g_xzzzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 567);

    auto g_xzzzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 568);

    auto g_xzzzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 569);

    auto g_xzzzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 570);

    auto g_xzzzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 571);

    auto g_xzzzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 572);

    auto g_xzzzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 573);

    auto g_xzzzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 574);

    auto g_xzzzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 575);

    auto g_xzzzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 576);

    auto g_xzzzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 577);

    auto g_xzzzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 578);

    auto g_xzzzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 579);

    auto g_xzzzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 580);

    auto g_xzzzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 581);

    auto g_xzzzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 582);

    auto g_xzzzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 583);

    auto g_xzzzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 584);

    auto g_xzzzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 585);

    auto g_xzzzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 586);

    auto g_xzzzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 587);

    #pragma omp simd aligned(g_xzzzzzz_xxxxx_0, g_xzzzzzz_xxxxy_0, g_xzzzzzz_xxxxz_0, g_xzzzzzz_xxxyy_0, g_xzzzzzz_xxxyz_0, g_xzzzzzz_xxxzz_0, g_xzzzzzz_xxyyy_0, g_xzzzzzz_xxyyz_0, g_xzzzzzz_xxyzz_0, g_xzzzzzz_xxzzz_0, g_xzzzzzz_xyyyy_0, g_xzzzzzz_xyyyz_0, g_xzzzzzz_xyyzz_0, g_xzzzzzz_xyzzz_0, g_xzzzzzz_xzzzz_0, g_xzzzzzz_yyyyy_0, g_xzzzzzz_yyyyz_0, g_xzzzzzz_yyyzz_0, g_xzzzzzz_yyzzz_0, g_xzzzzzz_yzzzz_0, g_xzzzzzz_zzzzz_0, g_zzzzzz_xxxx_1, g_zzzzzz_xxxxx_1, g_zzzzzz_xxxxy_1, g_zzzzzz_xxxxz_1, g_zzzzzz_xxxy_1, g_zzzzzz_xxxyy_1, g_zzzzzz_xxxyz_1, g_zzzzzz_xxxz_1, g_zzzzzz_xxxzz_1, g_zzzzzz_xxyy_1, g_zzzzzz_xxyyy_1, g_zzzzzz_xxyyz_1, g_zzzzzz_xxyz_1, g_zzzzzz_xxyzz_1, g_zzzzzz_xxzz_1, g_zzzzzz_xxzzz_1, g_zzzzzz_xyyy_1, g_zzzzzz_xyyyy_1, g_zzzzzz_xyyyz_1, g_zzzzzz_xyyz_1, g_zzzzzz_xyyzz_1, g_zzzzzz_xyzz_1, g_zzzzzz_xyzzz_1, g_zzzzzz_xzzz_1, g_zzzzzz_xzzzz_1, g_zzzzzz_yyyy_1, g_zzzzzz_yyyyy_1, g_zzzzzz_yyyyz_1, g_zzzzzz_yyyz_1, g_zzzzzz_yyyzz_1, g_zzzzzz_yyzz_1, g_zzzzzz_yyzzz_1, g_zzzzzz_yzzz_1, g_zzzzzz_yzzzz_1, g_zzzzzz_zzzz_1, g_zzzzzz_zzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzzzz_xxxxx_0[i] = 5.0 * g_zzzzzz_xxxx_1[i] * fe_0 + g_zzzzzz_xxxxx_1[i] * pa_x[i];

        g_xzzzzzz_xxxxy_0[i] = 4.0 * g_zzzzzz_xxxy_1[i] * fe_0 + g_zzzzzz_xxxxy_1[i] * pa_x[i];

        g_xzzzzzz_xxxxz_0[i] = 4.0 * g_zzzzzz_xxxz_1[i] * fe_0 + g_zzzzzz_xxxxz_1[i] * pa_x[i];

        g_xzzzzzz_xxxyy_0[i] = 3.0 * g_zzzzzz_xxyy_1[i] * fe_0 + g_zzzzzz_xxxyy_1[i] * pa_x[i];

        g_xzzzzzz_xxxyz_0[i] = 3.0 * g_zzzzzz_xxyz_1[i] * fe_0 + g_zzzzzz_xxxyz_1[i] * pa_x[i];

        g_xzzzzzz_xxxzz_0[i] = 3.0 * g_zzzzzz_xxzz_1[i] * fe_0 + g_zzzzzz_xxxzz_1[i] * pa_x[i];

        g_xzzzzzz_xxyyy_0[i] = 2.0 * g_zzzzzz_xyyy_1[i] * fe_0 + g_zzzzzz_xxyyy_1[i] * pa_x[i];

        g_xzzzzzz_xxyyz_0[i] = 2.0 * g_zzzzzz_xyyz_1[i] * fe_0 + g_zzzzzz_xxyyz_1[i] * pa_x[i];

        g_xzzzzzz_xxyzz_0[i] = 2.0 * g_zzzzzz_xyzz_1[i] * fe_0 + g_zzzzzz_xxyzz_1[i] * pa_x[i];

        g_xzzzzzz_xxzzz_0[i] = 2.0 * g_zzzzzz_xzzz_1[i] * fe_0 + g_zzzzzz_xxzzz_1[i] * pa_x[i];

        g_xzzzzzz_xyyyy_0[i] = g_zzzzzz_yyyy_1[i] * fe_0 + g_zzzzzz_xyyyy_1[i] * pa_x[i];

        g_xzzzzzz_xyyyz_0[i] = g_zzzzzz_yyyz_1[i] * fe_0 + g_zzzzzz_xyyyz_1[i] * pa_x[i];

        g_xzzzzzz_xyyzz_0[i] = g_zzzzzz_yyzz_1[i] * fe_0 + g_zzzzzz_xyyzz_1[i] * pa_x[i];

        g_xzzzzzz_xyzzz_0[i] = g_zzzzzz_yzzz_1[i] * fe_0 + g_zzzzzz_xyzzz_1[i] * pa_x[i];

        g_xzzzzzz_xzzzz_0[i] = g_zzzzzz_zzzz_1[i] * fe_0 + g_zzzzzz_xzzzz_1[i] * pa_x[i];

        g_xzzzzzz_yyyyy_0[i] = g_zzzzzz_yyyyy_1[i] * pa_x[i];

        g_xzzzzzz_yyyyz_0[i] = g_zzzzzz_yyyyz_1[i] * pa_x[i];

        g_xzzzzzz_yyyzz_0[i] = g_zzzzzz_yyyzz_1[i] * pa_x[i];

        g_xzzzzzz_yyzzz_0[i] = g_zzzzzz_yyzzz_1[i] * pa_x[i];

        g_xzzzzzz_yzzzz_0[i] = g_zzzzzz_yzzzz_1[i] * pa_x[i];

        g_xzzzzzz_zzzzz_0[i] = g_zzzzzz_zzzzz_1[i] * pa_x[i];
    }

    // Set up 588-609 components of targeted buffer : KH

    auto g_yyyyyyy_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 588);

    auto g_yyyyyyy_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 589);

    auto g_yyyyyyy_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 590);

    auto g_yyyyyyy_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 591);

    auto g_yyyyyyy_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 592);

    auto g_yyyyyyy_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 593);

    auto g_yyyyyyy_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 594);

    auto g_yyyyyyy_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 595);

    auto g_yyyyyyy_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 596);

    auto g_yyyyyyy_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 597);

    auto g_yyyyyyy_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 598);

    auto g_yyyyyyy_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 599);

    auto g_yyyyyyy_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 600);

    auto g_yyyyyyy_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 601);

    auto g_yyyyyyy_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 602);

    auto g_yyyyyyy_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 603);

    auto g_yyyyyyy_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 604);

    auto g_yyyyyyy_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 605);

    auto g_yyyyyyy_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 606);

    auto g_yyyyyyy_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 607);

    auto g_yyyyyyy_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 608);

    #pragma omp simd aligned(g_yyyyy_xxxxx_0, g_yyyyy_xxxxx_1, g_yyyyy_xxxxy_0, g_yyyyy_xxxxy_1, g_yyyyy_xxxxz_0, g_yyyyy_xxxxz_1, g_yyyyy_xxxyy_0, g_yyyyy_xxxyy_1, g_yyyyy_xxxyz_0, g_yyyyy_xxxyz_1, g_yyyyy_xxxzz_0, g_yyyyy_xxxzz_1, g_yyyyy_xxyyy_0, g_yyyyy_xxyyy_1, g_yyyyy_xxyyz_0, g_yyyyy_xxyyz_1, g_yyyyy_xxyzz_0, g_yyyyy_xxyzz_1, g_yyyyy_xxzzz_0, g_yyyyy_xxzzz_1, g_yyyyy_xyyyy_0, g_yyyyy_xyyyy_1, g_yyyyy_xyyyz_0, g_yyyyy_xyyyz_1, g_yyyyy_xyyzz_0, g_yyyyy_xyyzz_1, g_yyyyy_xyzzz_0, g_yyyyy_xyzzz_1, g_yyyyy_xzzzz_0, g_yyyyy_xzzzz_1, g_yyyyy_yyyyy_0, g_yyyyy_yyyyy_1, g_yyyyy_yyyyz_0, g_yyyyy_yyyyz_1, g_yyyyy_yyyzz_0, g_yyyyy_yyyzz_1, g_yyyyy_yyzzz_0, g_yyyyy_yyzzz_1, g_yyyyy_yzzzz_0, g_yyyyy_yzzzz_1, g_yyyyy_zzzzz_0, g_yyyyy_zzzzz_1, g_yyyyyy_xxxx_1, g_yyyyyy_xxxxx_1, g_yyyyyy_xxxxy_1, g_yyyyyy_xxxxz_1, g_yyyyyy_xxxy_1, g_yyyyyy_xxxyy_1, g_yyyyyy_xxxyz_1, g_yyyyyy_xxxz_1, g_yyyyyy_xxxzz_1, g_yyyyyy_xxyy_1, g_yyyyyy_xxyyy_1, g_yyyyyy_xxyyz_1, g_yyyyyy_xxyz_1, g_yyyyyy_xxyzz_1, g_yyyyyy_xxzz_1, g_yyyyyy_xxzzz_1, g_yyyyyy_xyyy_1, g_yyyyyy_xyyyy_1, g_yyyyyy_xyyyz_1, g_yyyyyy_xyyz_1, g_yyyyyy_xyyzz_1, g_yyyyyy_xyzz_1, g_yyyyyy_xyzzz_1, g_yyyyyy_xzzz_1, g_yyyyyy_xzzzz_1, g_yyyyyy_yyyy_1, g_yyyyyy_yyyyy_1, g_yyyyyy_yyyyz_1, g_yyyyyy_yyyz_1, g_yyyyyy_yyyzz_1, g_yyyyyy_yyzz_1, g_yyyyyy_yyzzz_1, g_yyyyyy_yzzz_1, g_yyyyyy_yzzzz_1, g_yyyyyy_zzzz_1, g_yyyyyy_zzzzz_1, g_yyyyyyy_xxxxx_0, g_yyyyyyy_xxxxy_0, g_yyyyyyy_xxxxz_0, g_yyyyyyy_xxxyy_0, g_yyyyyyy_xxxyz_0, g_yyyyyyy_xxxzz_0, g_yyyyyyy_xxyyy_0, g_yyyyyyy_xxyyz_0, g_yyyyyyy_xxyzz_0, g_yyyyyyy_xxzzz_0, g_yyyyyyy_xyyyy_0, g_yyyyyyy_xyyyz_0, g_yyyyyyy_xyyzz_0, g_yyyyyyy_xyzzz_0, g_yyyyyyy_xzzzz_0, g_yyyyyyy_yyyyy_0, g_yyyyyyy_yyyyz_0, g_yyyyyyy_yyyzz_0, g_yyyyyyy_yyzzz_0, g_yyyyyyy_yzzzz_0, g_yyyyyyy_zzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyyy_xxxxx_0[i] = 6.0 * g_yyyyy_xxxxx_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxxx_1[i] * fz_be_0 + g_yyyyyy_xxxxx_1[i] * pa_y[i];

        g_yyyyyyy_xxxxy_0[i] = 6.0 * g_yyyyy_xxxxy_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxxy_1[i] * fz_be_0 + g_yyyyyy_xxxx_1[i] * fe_0 + g_yyyyyy_xxxxy_1[i] * pa_y[i];

        g_yyyyyyy_xxxxz_0[i] = 6.0 * g_yyyyy_xxxxz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxxz_1[i] * fz_be_0 + g_yyyyyy_xxxxz_1[i] * pa_y[i];

        g_yyyyyyy_xxxyy_0[i] = 6.0 * g_yyyyy_xxxyy_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxyy_1[i] * fz_be_0 + 2.0 * g_yyyyyy_xxxy_1[i] * fe_0 + g_yyyyyy_xxxyy_1[i] * pa_y[i];

        g_yyyyyyy_xxxyz_0[i] = 6.0 * g_yyyyy_xxxyz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxyz_1[i] * fz_be_0 + g_yyyyyy_xxxz_1[i] * fe_0 + g_yyyyyy_xxxyz_1[i] * pa_y[i];

        g_yyyyyyy_xxxzz_0[i] = 6.0 * g_yyyyy_xxxzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxxzz_1[i] * fz_be_0 + g_yyyyyy_xxxzz_1[i] * pa_y[i];

        g_yyyyyyy_xxyyy_0[i] = 6.0 * g_yyyyy_xxyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_xxyyy_1[i] * fz_be_0 + 3.0 * g_yyyyyy_xxyy_1[i] * fe_0 + g_yyyyyy_xxyyy_1[i] * pa_y[i];

        g_yyyyyyy_xxyyz_0[i] = 6.0 * g_yyyyy_xxyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_xxyz_1[i] * fe_0 + g_yyyyyy_xxyyz_1[i] * pa_y[i];

        g_yyyyyyy_xxyzz_0[i] = 6.0 * g_yyyyy_xxyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxyzz_1[i] * fz_be_0 + g_yyyyyy_xxzz_1[i] * fe_0 + g_yyyyyy_xxyzz_1[i] * pa_y[i];

        g_yyyyyyy_xxzzz_0[i] = 6.0 * g_yyyyy_xxzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxzzz_1[i] * fz_be_0 + g_yyyyyy_xxzzz_1[i] * pa_y[i];

        g_yyyyyyy_xyyyy_0[i] = 6.0 * g_yyyyy_xyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_xyyyy_1[i] * fz_be_0 + 4.0 * g_yyyyyy_xyyy_1[i] * fe_0 + g_yyyyyy_xyyyy_1[i] * pa_y[i];

        g_yyyyyyy_xyyyz_0[i] = 6.0 * g_yyyyy_xyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_xyyz_1[i] * fe_0 + g_yyyyyy_xyyyz_1[i] * pa_y[i];

        g_yyyyyyy_xyyzz_0[i] = 6.0 * g_yyyyy_xyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_xyzz_1[i] * fe_0 + g_yyyyyy_xyyzz_1[i] * pa_y[i];

        g_yyyyyyy_xyzzz_0[i] = 6.0 * g_yyyyy_xyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xyzzz_1[i] * fz_be_0 + g_yyyyyy_xzzz_1[i] * fe_0 + g_yyyyyy_xyzzz_1[i] * pa_y[i];

        g_yyyyyyy_xzzzz_0[i] = 6.0 * g_yyyyy_xzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xzzzz_1[i] * fz_be_0 + g_yyyyyy_xzzzz_1[i] * pa_y[i];

        g_yyyyyyy_yyyyy_0[i] = 6.0 * g_yyyyy_yyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_yyyyy_1[i] * fz_be_0 + 5.0 * g_yyyyyy_yyyy_1[i] * fe_0 + g_yyyyyy_yyyyy_1[i] * pa_y[i];

        g_yyyyyyy_yyyyz_0[i] = 6.0 * g_yyyyy_yyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_yyyz_1[i] * fe_0 + g_yyyyyy_yyyyz_1[i] * pa_y[i];

        g_yyyyyyy_yyyzz_0[i] = 6.0 * g_yyyyy_yyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_yyzz_1[i] * fe_0 + g_yyyyyy_yyyzz_1[i] * pa_y[i];

        g_yyyyyyy_yyzzz_0[i] = 6.0 * g_yyyyy_yyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_yzzz_1[i] * fe_0 + g_yyyyyy_yyzzz_1[i] * pa_y[i];

        g_yyyyyyy_yzzzz_0[i] = 6.0 * g_yyyyy_yzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_yzzzz_1[i] * fz_be_0 + g_yyyyyy_zzzz_1[i] * fe_0 + g_yyyyyy_yzzzz_1[i] * pa_y[i];

        g_yyyyyyy_zzzzz_0[i] = 6.0 * g_yyyyy_zzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_zzzzz_1[i] * fz_be_0 + g_yyyyyy_zzzzz_1[i] * pa_y[i];
    }

    // Set up 609-630 components of targeted buffer : KH

    auto g_yyyyyyz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 609);

    auto g_yyyyyyz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 610);

    auto g_yyyyyyz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 611);

    auto g_yyyyyyz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 612);

    auto g_yyyyyyz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 613);

    auto g_yyyyyyz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 614);

    auto g_yyyyyyz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 615);

    auto g_yyyyyyz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 616);

    auto g_yyyyyyz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 617);

    auto g_yyyyyyz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 618);

    auto g_yyyyyyz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 619);

    auto g_yyyyyyz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 620);

    auto g_yyyyyyz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 621);

    auto g_yyyyyyz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 622);

    auto g_yyyyyyz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 623);

    auto g_yyyyyyz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 624);

    auto g_yyyyyyz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 625);

    auto g_yyyyyyz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 626);

    auto g_yyyyyyz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 627);

    auto g_yyyyyyz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 628);

    auto g_yyyyyyz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 629);

    #pragma omp simd aligned(g_yyyyyy_xxxx_1, g_yyyyyy_xxxxx_1, g_yyyyyy_xxxxy_1, g_yyyyyy_xxxxz_1, g_yyyyyy_xxxy_1, g_yyyyyy_xxxyy_1, g_yyyyyy_xxxyz_1, g_yyyyyy_xxxz_1, g_yyyyyy_xxxzz_1, g_yyyyyy_xxyy_1, g_yyyyyy_xxyyy_1, g_yyyyyy_xxyyz_1, g_yyyyyy_xxyz_1, g_yyyyyy_xxyzz_1, g_yyyyyy_xxzz_1, g_yyyyyy_xxzzz_1, g_yyyyyy_xyyy_1, g_yyyyyy_xyyyy_1, g_yyyyyy_xyyyz_1, g_yyyyyy_xyyz_1, g_yyyyyy_xyyzz_1, g_yyyyyy_xyzz_1, g_yyyyyy_xyzzz_1, g_yyyyyy_xzzz_1, g_yyyyyy_xzzzz_1, g_yyyyyy_yyyy_1, g_yyyyyy_yyyyy_1, g_yyyyyy_yyyyz_1, g_yyyyyy_yyyz_1, g_yyyyyy_yyyzz_1, g_yyyyyy_yyzz_1, g_yyyyyy_yyzzz_1, g_yyyyyy_yzzz_1, g_yyyyyy_yzzzz_1, g_yyyyyy_zzzz_1, g_yyyyyy_zzzzz_1, g_yyyyyyz_xxxxx_0, g_yyyyyyz_xxxxy_0, g_yyyyyyz_xxxxz_0, g_yyyyyyz_xxxyy_0, g_yyyyyyz_xxxyz_0, g_yyyyyyz_xxxzz_0, g_yyyyyyz_xxyyy_0, g_yyyyyyz_xxyyz_0, g_yyyyyyz_xxyzz_0, g_yyyyyyz_xxzzz_0, g_yyyyyyz_xyyyy_0, g_yyyyyyz_xyyyz_0, g_yyyyyyz_xyyzz_0, g_yyyyyyz_xyzzz_0, g_yyyyyyz_xzzzz_0, g_yyyyyyz_yyyyy_0, g_yyyyyyz_yyyyz_0, g_yyyyyyz_yyyzz_0, g_yyyyyyz_yyzzz_0, g_yyyyyyz_yzzzz_0, g_yyyyyyz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyyyz_xxxxx_0[i] = g_yyyyyy_xxxxx_1[i] * pa_z[i];

        g_yyyyyyz_xxxxy_0[i] = g_yyyyyy_xxxxy_1[i] * pa_z[i];

        g_yyyyyyz_xxxxz_0[i] = g_yyyyyy_xxxx_1[i] * fe_0 + g_yyyyyy_xxxxz_1[i] * pa_z[i];

        g_yyyyyyz_xxxyy_0[i] = g_yyyyyy_xxxyy_1[i] * pa_z[i];

        g_yyyyyyz_xxxyz_0[i] = g_yyyyyy_xxxy_1[i] * fe_0 + g_yyyyyy_xxxyz_1[i] * pa_z[i];

        g_yyyyyyz_xxxzz_0[i] = 2.0 * g_yyyyyy_xxxz_1[i] * fe_0 + g_yyyyyy_xxxzz_1[i] * pa_z[i];

        g_yyyyyyz_xxyyy_0[i] = g_yyyyyy_xxyyy_1[i] * pa_z[i];

        g_yyyyyyz_xxyyz_0[i] = g_yyyyyy_xxyy_1[i] * fe_0 + g_yyyyyy_xxyyz_1[i] * pa_z[i];

        g_yyyyyyz_xxyzz_0[i] = 2.0 * g_yyyyyy_xxyz_1[i] * fe_0 + g_yyyyyy_xxyzz_1[i] * pa_z[i];

        g_yyyyyyz_xxzzz_0[i] = 3.0 * g_yyyyyy_xxzz_1[i] * fe_0 + g_yyyyyy_xxzzz_1[i] * pa_z[i];

        g_yyyyyyz_xyyyy_0[i] = g_yyyyyy_xyyyy_1[i] * pa_z[i];

        g_yyyyyyz_xyyyz_0[i] = g_yyyyyy_xyyy_1[i] * fe_0 + g_yyyyyy_xyyyz_1[i] * pa_z[i];

        g_yyyyyyz_xyyzz_0[i] = 2.0 * g_yyyyyy_xyyz_1[i] * fe_0 + g_yyyyyy_xyyzz_1[i] * pa_z[i];

        g_yyyyyyz_xyzzz_0[i] = 3.0 * g_yyyyyy_xyzz_1[i] * fe_0 + g_yyyyyy_xyzzz_1[i] * pa_z[i];

        g_yyyyyyz_xzzzz_0[i] = 4.0 * g_yyyyyy_xzzz_1[i] * fe_0 + g_yyyyyy_xzzzz_1[i] * pa_z[i];

        g_yyyyyyz_yyyyy_0[i] = g_yyyyyy_yyyyy_1[i] * pa_z[i];

        g_yyyyyyz_yyyyz_0[i] = g_yyyyyy_yyyy_1[i] * fe_0 + g_yyyyyy_yyyyz_1[i] * pa_z[i];

        g_yyyyyyz_yyyzz_0[i] = 2.0 * g_yyyyyy_yyyz_1[i] * fe_0 + g_yyyyyy_yyyzz_1[i] * pa_z[i];

        g_yyyyyyz_yyzzz_0[i] = 3.0 * g_yyyyyy_yyzz_1[i] * fe_0 + g_yyyyyy_yyzzz_1[i] * pa_z[i];

        g_yyyyyyz_yzzzz_0[i] = 4.0 * g_yyyyyy_yzzz_1[i] * fe_0 + g_yyyyyy_yzzzz_1[i] * pa_z[i];

        g_yyyyyyz_zzzzz_0[i] = 5.0 * g_yyyyyy_zzzz_1[i] * fe_0 + g_yyyyyy_zzzzz_1[i] * pa_z[i];
    }

    // Set up 630-651 components of targeted buffer : KH

    auto g_yyyyyzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 630);

    auto g_yyyyyzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 631);

    auto g_yyyyyzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 632);

    auto g_yyyyyzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 633);

    auto g_yyyyyzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 634);

    auto g_yyyyyzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 635);

    auto g_yyyyyzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 636);

    auto g_yyyyyzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 637);

    auto g_yyyyyzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 638);

    auto g_yyyyyzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 639);

    auto g_yyyyyzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 640);

    auto g_yyyyyzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 641);

    auto g_yyyyyzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 642);

    auto g_yyyyyzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 643);

    auto g_yyyyyzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 644);

    auto g_yyyyyzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 645);

    auto g_yyyyyzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 646);

    auto g_yyyyyzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 647);

    auto g_yyyyyzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 648);

    auto g_yyyyyzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 649);

    auto g_yyyyyzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 650);

    #pragma omp simd aligned(g_yyyyy_xxxxy_0, g_yyyyy_xxxxy_1, g_yyyyy_xxxyy_0, g_yyyyy_xxxyy_1, g_yyyyy_xxyyy_0, g_yyyyy_xxyyy_1, g_yyyyy_xyyyy_0, g_yyyyy_xyyyy_1, g_yyyyy_yyyyy_0, g_yyyyy_yyyyy_1, g_yyyyyz_xxxxy_1, g_yyyyyz_xxxyy_1, g_yyyyyz_xxyyy_1, g_yyyyyz_xyyyy_1, g_yyyyyz_yyyyy_1, g_yyyyyzz_xxxxx_0, g_yyyyyzz_xxxxy_0, g_yyyyyzz_xxxxz_0, g_yyyyyzz_xxxyy_0, g_yyyyyzz_xxxyz_0, g_yyyyyzz_xxxzz_0, g_yyyyyzz_xxyyy_0, g_yyyyyzz_xxyyz_0, g_yyyyyzz_xxyzz_0, g_yyyyyzz_xxzzz_0, g_yyyyyzz_xyyyy_0, g_yyyyyzz_xyyyz_0, g_yyyyyzz_xyyzz_0, g_yyyyyzz_xyzzz_0, g_yyyyyzz_xzzzz_0, g_yyyyyzz_yyyyy_0, g_yyyyyzz_yyyyz_0, g_yyyyyzz_yyyzz_0, g_yyyyyzz_yyzzz_0, g_yyyyyzz_yzzzz_0, g_yyyyyzz_zzzzz_0, g_yyyyzz_xxxxx_1, g_yyyyzz_xxxxz_1, g_yyyyzz_xxxyz_1, g_yyyyzz_xxxz_1, g_yyyyzz_xxxzz_1, g_yyyyzz_xxyyz_1, g_yyyyzz_xxyz_1, g_yyyyzz_xxyzz_1, g_yyyyzz_xxzz_1, g_yyyyzz_xxzzz_1, g_yyyyzz_xyyyz_1, g_yyyyzz_xyyz_1, g_yyyyzz_xyyzz_1, g_yyyyzz_xyzz_1, g_yyyyzz_xyzzz_1, g_yyyyzz_xzzz_1, g_yyyyzz_xzzzz_1, g_yyyyzz_yyyyz_1, g_yyyyzz_yyyz_1, g_yyyyzz_yyyzz_1, g_yyyyzz_yyzz_1, g_yyyyzz_yyzzz_1, g_yyyyzz_yzzz_1, g_yyyyzz_yzzzz_1, g_yyyyzz_zzzz_1, g_yyyyzz_zzzzz_1, g_yyyzz_xxxxx_0, g_yyyzz_xxxxx_1, g_yyyzz_xxxxz_0, g_yyyzz_xxxxz_1, g_yyyzz_xxxyz_0, g_yyyzz_xxxyz_1, g_yyyzz_xxxzz_0, g_yyyzz_xxxzz_1, g_yyyzz_xxyyz_0, g_yyyzz_xxyyz_1, g_yyyzz_xxyzz_0, g_yyyzz_xxyzz_1, g_yyyzz_xxzzz_0, g_yyyzz_xxzzz_1, g_yyyzz_xyyyz_0, g_yyyzz_xyyyz_1, g_yyyzz_xyyzz_0, g_yyyzz_xyyzz_1, g_yyyzz_xyzzz_0, g_yyyzz_xyzzz_1, g_yyyzz_xzzzz_0, g_yyyzz_xzzzz_1, g_yyyzz_yyyyz_0, g_yyyzz_yyyyz_1, g_yyyzz_yyyzz_0, g_yyyzz_yyyzz_1, g_yyyzz_yyzzz_0, g_yyyzz_yyzzz_1, g_yyyzz_yzzzz_0, g_yyyzz_yzzzz_1, g_yyyzz_zzzzz_0, g_yyyzz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyzz_xxxxx_0[i] = 4.0 * g_yyyzz_xxxxx_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxxx_1[i] * fz_be_0 + g_yyyyzz_xxxxx_1[i] * pa_y[i];

        g_yyyyyzz_xxxxy_0[i] = g_yyyyy_xxxxy_0[i] * fbe_0 - g_yyyyy_xxxxy_1[i] * fz_be_0 + g_yyyyyz_xxxxy_1[i] * pa_z[i];

        g_yyyyyzz_xxxxz_0[i] = 4.0 * g_yyyzz_xxxxz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxxz_1[i] * fz_be_0 + g_yyyyzz_xxxxz_1[i] * pa_y[i];

        g_yyyyyzz_xxxyy_0[i] = g_yyyyy_xxxyy_0[i] * fbe_0 - g_yyyyy_xxxyy_1[i] * fz_be_0 + g_yyyyyz_xxxyy_1[i] * pa_z[i];

        g_yyyyyzz_xxxyz_0[i] = 4.0 * g_yyyzz_xxxyz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxyz_1[i] * fz_be_0 + g_yyyyzz_xxxz_1[i] * fe_0 + g_yyyyzz_xxxyz_1[i] * pa_y[i];

        g_yyyyyzz_xxxzz_0[i] = 4.0 * g_yyyzz_xxxzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxxzz_1[i] * fz_be_0 + g_yyyyzz_xxxzz_1[i] * pa_y[i];

        g_yyyyyzz_xxyyy_0[i] = g_yyyyy_xxyyy_0[i] * fbe_0 - g_yyyyy_xxyyy_1[i] * fz_be_0 + g_yyyyyz_xxyyy_1[i] * pa_z[i];

        g_yyyyyzz_xxyyz_0[i] = 4.0 * g_yyyzz_xxyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_xxyz_1[i] * fe_0 + g_yyyyzz_xxyyz_1[i] * pa_y[i];

        g_yyyyyzz_xxyzz_0[i] = 4.0 * g_yyyzz_xxyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxyzz_1[i] * fz_be_0 + g_yyyyzz_xxzz_1[i] * fe_0 + g_yyyyzz_xxyzz_1[i] * pa_y[i];

        g_yyyyyzz_xxzzz_0[i] = 4.0 * g_yyyzz_xxzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxzzz_1[i] * fz_be_0 + g_yyyyzz_xxzzz_1[i] * pa_y[i];

        g_yyyyyzz_xyyyy_0[i] = g_yyyyy_xyyyy_0[i] * fbe_0 - g_yyyyy_xyyyy_1[i] * fz_be_0 + g_yyyyyz_xyyyy_1[i] * pa_z[i];

        g_yyyyyzz_xyyyz_0[i] = 4.0 * g_yyyzz_xyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_xyyz_1[i] * fe_0 + g_yyyyzz_xyyyz_1[i] * pa_y[i];

        g_yyyyyzz_xyyzz_0[i] = 4.0 * g_yyyzz_xyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_xyzz_1[i] * fe_0 + g_yyyyzz_xyyzz_1[i] * pa_y[i];

        g_yyyyyzz_xyzzz_0[i] = 4.0 * g_yyyzz_xyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xyzzz_1[i] * fz_be_0 + g_yyyyzz_xzzz_1[i] * fe_0 + g_yyyyzz_xyzzz_1[i] * pa_y[i];

        g_yyyyyzz_xzzzz_0[i] = 4.0 * g_yyyzz_xzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xzzzz_1[i] * fz_be_0 + g_yyyyzz_xzzzz_1[i] * pa_y[i];

        g_yyyyyzz_yyyyy_0[i] = g_yyyyy_yyyyy_0[i] * fbe_0 - g_yyyyy_yyyyy_1[i] * fz_be_0 + g_yyyyyz_yyyyy_1[i] * pa_z[i];

        g_yyyyyzz_yyyyz_0[i] = 4.0 * g_yyyzz_yyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_yyyz_1[i] * fe_0 + g_yyyyzz_yyyyz_1[i] * pa_y[i];

        g_yyyyyzz_yyyzz_0[i] = 4.0 * g_yyyzz_yyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_yyzz_1[i] * fe_0 + g_yyyyzz_yyyzz_1[i] * pa_y[i];

        g_yyyyyzz_yyzzz_0[i] = 4.0 * g_yyyzz_yyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_yzzz_1[i] * fe_0 + g_yyyyzz_yyzzz_1[i] * pa_y[i];

        g_yyyyyzz_yzzzz_0[i] = 4.0 * g_yyyzz_yzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_yzzzz_1[i] * fz_be_0 + g_yyyyzz_zzzz_1[i] * fe_0 + g_yyyyzz_yzzzz_1[i] * pa_y[i];

        g_yyyyyzz_zzzzz_0[i] = 4.0 * g_yyyzz_zzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_zzzzz_1[i] * fz_be_0 + g_yyyyzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 651-672 components of targeted buffer : KH

    auto g_yyyyzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 651);

    auto g_yyyyzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 652);

    auto g_yyyyzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 653);

    auto g_yyyyzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 654);

    auto g_yyyyzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 655);

    auto g_yyyyzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 656);

    auto g_yyyyzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 657);

    auto g_yyyyzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 658);

    auto g_yyyyzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 659);

    auto g_yyyyzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 660);

    auto g_yyyyzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 661);

    auto g_yyyyzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 662);

    auto g_yyyyzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 663);

    auto g_yyyyzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 664);

    auto g_yyyyzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 665);

    auto g_yyyyzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 666);

    auto g_yyyyzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 667);

    auto g_yyyyzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 668);

    auto g_yyyyzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 669);

    auto g_yyyyzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 670);

    auto g_yyyyzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 671);

    #pragma omp simd aligned(g_yyyyz_xxxxy_0, g_yyyyz_xxxxy_1, g_yyyyz_xxxyy_0, g_yyyyz_xxxyy_1, g_yyyyz_xxyyy_0, g_yyyyz_xxyyy_1, g_yyyyz_xyyyy_0, g_yyyyz_xyyyy_1, g_yyyyz_yyyyy_0, g_yyyyz_yyyyy_1, g_yyyyzz_xxxxy_1, g_yyyyzz_xxxyy_1, g_yyyyzz_xxyyy_1, g_yyyyzz_xyyyy_1, g_yyyyzz_yyyyy_1, g_yyyyzzz_xxxxx_0, g_yyyyzzz_xxxxy_0, g_yyyyzzz_xxxxz_0, g_yyyyzzz_xxxyy_0, g_yyyyzzz_xxxyz_0, g_yyyyzzz_xxxzz_0, g_yyyyzzz_xxyyy_0, g_yyyyzzz_xxyyz_0, g_yyyyzzz_xxyzz_0, g_yyyyzzz_xxzzz_0, g_yyyyzzz_xyyyy_0, g_yyyyzzz_xyyyz_0, g_yyyyzzz_xyyzz_0, g_yyyyzzz_xyzzz_0, g_yyyyzzz_xzzzz_0, g_yyyyzzz_yyyyy_0, g_yyyyzzz_yyyyz_0, g_yyyyzzz_yyyzz_0, g_yyyyzzz_yyzzz_0, g_yyyyzzz_yzzzz_0, g_yyyyzzz_zzzzz_0, g_yyyzzz_xxxxx_1, g_yyyzzz_xxxxz_1, g_yyyzzz_xxxyz_1, g_yyyzzz_xxxz_1, g_yyyzzz_xxxzz_1, g_yyyzzz_xxyyz_1, g_yyyzzz_xxyz_1, g_yyyzzz_xxyzz_1, g_yyyzzz_xxzz_1, g_yyyzzz_xxzzz_1, g_yyyzzz_xyyyz_1, g_yyyzzz_xyyz_1, g_yyyzzz_xyyzz_1, g_yyyzzz_xyzz_1, g_yyyzzz_xyzzz_1, g_yyyzzz_xzzz_1, g_yyyzzz_xzzzz_1, g_yyyzzz_yyyyz_1, g_yyyzzz_yyyz_1, g_yyyzzz_yyyzz_1, g_yyyzzz_yyzz_1, g_yyyzzz_yyzzz_1, g_yyyzzz_yzzz_1, g_yyyzzz_yzzzz_1, g_yyyzzz_zzzz_1, g_yyyzzz_zzzzz_1, g_yyzzz_xxxxx_0, g_yyzzz_xxxxx_1, g_yyzzz_xxxxz_0, g_yyzzz_xxxxz_1, g_yyzzz_xxxyz_0, g_yyzzz_xxxyz_1, g_yyzzz_xxxzz_0, g_yyzzz_xxxzz_1, g_yyzzz_xxyyz_0, g_yyzzz_xxyyz_1, g_yyzzz_xxyzz_0, g_yyzzz_xxyzz_1, g_yyzzz_xxzzz_0, g_yyzzz_xxzzz_1, g_yyzzz_xyyyz_0, g_yyzzz_xyyyz_1, g_yyzzz_xyyzz_0, g_yyzzz_xyyzz_1, g_yyzzz_xyzzz_0, g_yyzzz_xyzzz_1, g_yyzzz_xzzzz_0, g_yyzzz_xzzzz_1, g_yyzzz_yyyyz_0, g_yyzzz_yyyyz_1, g_yyzzz_yyyzz_0, g_yyzzz_yyyzz_1, g_yyzzz_yyzzz_0, g_yyzzz_yyzzz_1, g_yyzzz_yzzzz_0, g_yyzzz_yzzzz_1, g_yyzzz_zzzzz_0, g_yyzzz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyzzz_xxxxx_0[i] = 3.0 * g_yyzzz_xxxxx_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxxx_1[i] * fz_be_0 + g_yyyzzz_xxxxx_1[i] * pa_y[i];

        g_yyyyzzz_xxxxy_0[i] = 2.0 * g_yyyyz_xxxxy_0[i] * fbe_0 - 2.0 * g_yyyyz_xxxxy_1[i] * fz_be_0 + g_yyyyzz_xxxxy_1[i] * pa_z[i];

        g_yyyyzzz_xxxxz_0[i] = 3.0 * g_yyzzz_xxxxz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxxz_1[i] * fz_be_0 + g_yyyzzz_xxxxz_1[i] * pa_y[i];

        g_yyyyzzz_xxxyy_0[i] = 2.0 * g_yyyyz_xxxyy_0[i] * fbe_0 - 2.0 * g_yyyyz_xxxyy_1[i] * fz_be_0 + g_yyyyzz_xxxyy_1[i] * pa_z[i];

        g_yyyyzzz_xxxyz_0[i] = 3.0 * g_yyzzz_xxxyz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxyz_1[i] * fz_be_0 + g_yyyzzz_xxxz_1[i] * fe_0 + g_yyyzzz_xxxyz_1[i] * pa_y[i];

        g_yyyyzzz_xxxzz_0[i] = 3.0 * g_yyzzz_xxxzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxxzz_1[i] * fz_be_0 + g_yyyzzz_xxxzz_1[i] * pa_y[i];

        g_yyyyzzz_xxyyy_0[i] = 2.0 * g_yyyyz_xxyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_xxyyy_1[i] * fz_be_0 + g_yyyyzz_xxyyy_1[i] * pa_z[i];

        g_yyyyzzz_xxyyz_0[i] = 3.0 * g_yyzzz_xxyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_xxyz_1[i] * fe_0 + g_yyyzzz_xxyyz_1[i] * pa_y[i];

        g_yyyyzzz_xxyzz_0[i] = 3.0 * g_yyzzz_xxyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxyzz_1[i] * fz_be_0 + g_yyyzzz_xxzz_1[i] * fe_0 + g_yyyzzz_xxyzz_1[i] * pa_y[i];

        g_yyyyzzz_xxzzz_0[i] = 3.0 * g_yyzzz_xxzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxzzz_1[i] * fz_be_0 + g_yyyzzz_xxzzz_1[i] * pa_y[i];

        g_yyyyzzz_xyyyy_0[i] = 2.0 * g_yyyyz_xyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_xyyyy_1[i] * fz_be_0 + g_yyyyzz_xyyyy_1[i] * pa_z[i];

        g_yyyyzzz_xyyyz_0[i] = 3.0 * g_yyzzz_xyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_xyyz_1[i] * fe_0 + g_yyyzzz_xyyyz_1[i] * pa_y[i];

        g_yyyyzzz_xyyzz_0[i] = 3.0 * g_yyzzz_xyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_xyzz_1[i] * fe_0 + g_yyyzzz_xyyzz_1[i] * pa_y[i];

        g_yyyyzzz_xyzzz_0[i] = 3.0 * g_yyzzz_xyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xyzzz_1[i] * fz_be_0 + g_yyyzzz_xzzz_1[i] * fe_0 + g_yyyzzz_xyzzz_1[i] * pa_y[i];

        g_yyyyzzz_xzzzz_0[i] = 3.0 * g_yyzzz_xzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xzzzz_1[i] * fz_be_0 + g_yyyzzz_xzzzz_1[i] * pa_y[i];

        g_yyyyzzz_yyyyy_0[i] = 2.0 * g_yyyyz_yyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_yyyyy_1[i] * fz_be_0 + g_yyyyzz_yyyyy_1[i] * pa_z[i];

        g_yyyyzzz_yyyyz_0[i] = 3.0 * g_yyzzz_yyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_yyyz_1[i] * fe_0 + g_yyyzzz_yyyyz_1[i] * pa_y[i];

        g_yyyyzzz_yyyzz_0[i] = 3.0 * g_yyzzz_yyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_yyzz_1[i] * fe_0 + g_yyyzzz_yyyzz_1[i] * pa_y[i];

        g_yyyyzzz_yyzzz_0[i] = 3.0 * g_yyzzz_yyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_yzzz_1[i] * fe_0 + g_yyyzzz_yyzzz_1[i] * pa_y[i];

        g_yyyyzzz_yzzzz_0[i] = 3.0 * g_yyzzz_yzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_yzzzz_1[i] * fz_be_0 + g_yyyzzz_zzzz_1[i] * fe_0 + g_yyyzzz_yzzzz_1[i] * pa_y[i];

        g_yyyyzzz_zzzzz_0[i] = 3.0 * g_yyzzz_zzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_zzzzz_1[i] * fz_be_0 + g_yyyzzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 672-693 components of targeted buffer : KH

    auto g_yyyzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 672);

    auto g_yyyzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 673);

    auto g_yyyzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 674);

    auto g_yyyzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 675);

    auto g_yyyzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 676);

    auto g_yyyzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 677);

    auto g_yyyzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 678);

    auto g_yyyzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 679);

    auto g_yyyzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 680);

    auto g_yyyzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 681);

    auto g_yyyzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 682);

    auto g_yyyzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 683);

    auto g_yyyzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 684);

    auto g_yyyzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 685);

    auto g_yyyzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 686);

    auto g_yyyzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 687);

    auto g_yyyzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 688);

    auto g_yyyzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 689);

    auto g_yyyzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 690);

    auto g_yyyzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 691);

    auto g_yyyzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 692);

    #pragma omp simd aligned(g_yyyzz_xxxxy_0, g_yyyzz_xxxxy_1, g_yyyzz_xxxyy_0, g_yyyzz_xxxyy_1, g_yyyzz_xxyyy_0, g_yyyzz_xxyyy_1, g_yyyzz_xyyyy_0, g_yyyzz_xyyyy_1, g_yyyzz_yyyyy_0, g_yyyzz_yyyyy_1, g_yyyzzz_xxxxy_1, g_yyyzzz_xxxyy_1, g_yyyzzz_xxyyy_1, g_yyyzzz_xyyyy_1, g_yyyzzz_yyyyy_1, g_yyyzzzz_xxxxx_0, g_yyyzzzz_xxxxy_0, g_yyyzzzz_xxxxz_0, g_yyyzzzz_xxxyy_0, g_yyyzzzz_xxxyz_0, g_yyyzzzz_xxxzz_0, g_yyyzzzz_xxyyy_0, g_yyyzzzz_xxyyz_0, g_yyyzzzz_xxyzz_0, g_yyyzzzz_xxzzz_0, g_yyyzzzz_xyyyy_0, g_yyyzzzz_xyyyz_0, g_yyyzzzz_xyyzz_0, g_yyyzzzz_xyzzz_0, g_yyyzzzz_xzzzz_0, g_yyyzzzz_yyyyy_0, g_yyyzzzz_yyyyz_0, g_yyyzzzz_yyyzz_0, g_yyyzzzz_yyzzz_0, g_yyyzzzz_yzzzz_0, g_yyyzzzz_zzzzz_0, g_yyzzzz_xxxxx_1, g_yyzzzz_xxxxz_1, g_yyzzzz_xxxyz_1, g_yyzzzz_xxxz_1, g_yyzzzz_xxxzz_1, g_yyzzzz_xxyyz_1, g_yyzzzz_xxyz_1, g_yyzzzz_xxyzz_1, g_yyzzzz_xxzz_1, g_yyzzzz_xxzzz_1, g_yyzzzz_xyyyz_1, g_yyzzzz_xyyz_1, g_yyzzzz_xyyzz_1, g_yyzzzz_xyzz_1, g_yyzzzz_xyzzz_1, g_yyzzzz_xzzz_1, g_yyzzzz_xzzzz_1, g_yyzzzz_yyyyz_1, g_yyzzzz_yyyz_1, g_yyzzzz_yyyzz_1, g_yyzzzz_yyzz_1, g_yyzzzz_yyzzz_1, g_yyzzzz_yzzz_1, g_yyzzzz_yzzzz_1, g_yyzzzz_zzzz_1, g_yyzzzz_zzzzz_1, g_yzzzz_xxxxx_0, g_yzzzz_xxxxx_1, g_yzzzz_xxxxz_0, g_yzzzz_xxxxz_1, g_yzzzz_xxxyz_0, g_yzzzz_xxxyz_1, g_yzzzz_xxxzz_0, g_yzzzz_xxxzz_1, g_yzzzz_xxyyz_0, g_yzzzz_xxyyz_1, g_yzzzz_xxyzz_0, g_yzzzz_xxyzz_1, g_yzzzz_xxzzz_0, g_yzzzz_xxzzz_1, g_yzzzz_xyyyz_0, g_yzzzz_xyyyz_1, g_yzzzz_xyyzz_0, g_yzzzz_xyyzz_1, g_yzzzz_xyzzz_0, g_yzzzz_xyzzz_1, g_yzzzz_xzzzz_0, g_yzzzz_xzzzz_1, g_yzzzz_yyyyz_0, g_yzzzz_yyyyz_1, g_yzzzz_yyyzz_0, g_yzzzz_yyyzz_1, g_yzzzz_yyzzz_0, g_yzzzz_yyzzz_1, g_yzzzz_yzzzz_0, g_yzzzz_yzzzz_1, g_yzzzz_zzzzz_0, g_yzzzz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzzzz_xxxxx_0[i] = 2.0 * g_yzzzz_xxxxx_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxxx_1[i] * fz_be_0 + g_yyzzzz_xxxxx_1[i] * pa_y[i];

        g_yyyzzzz_xxxxy_0[i] = 3.0 * g_yyyzz_xxxxy_0[i] * fbe_0 - 3.0 * g_yyyzz_xxxxy_1[i] * fz_be_0 + g_yyyzzz_xxxxy_1[i] * pa_z[i];

        g_yyyzzzz_xxxxz_0[i] = 2.0 * g_yzzzz_xxxxz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxxz_1[i] * fz_be_0 + g_yyzzzz_xxxxz_1[i] * pa_y[i];

        g_yyyzzzz_xxxyy_0[i] = 3.0 * g_yyyzz_xxxyy_0[i] * fbe_0 - 3.0 * g_yyyzz_xxxyy_1[i] * fz_be_0 + g_yyyzzz_xxxyy_1[i] * pa_z[i];

        g_yyyzzzz_xxxyz_0[i] = 2.0 * g_yzzzz_xxxyz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxyz_1[i] * fz_be_0 + g_yyzzzz_xxxz_1[i] * fe_0 + g_yyzzzz_xxxyz_1[i] * pa_y[i];

        g_yyyzzzz_xxxzz_0[i] = 2.0 * g_yzzzz_xxxzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxxzz_1[i] * fz_be_0 + g_yyzzzz_xxxzz_1[i] * pa_y[i];

        g_yyyzzzz_xxyyy_0[i] = 3.0 * g_yyyzz_xxyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_xxyyy_1[i] * fz_be_0 + g_yyyzzz_xxyyy_1[i] * pa_z[i];

        g_yyyzzzz_xxyyz_0[i] = 2.0 * g_yzzzz_xxyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_xxyz_1[i] * fe_0 + g_yyzzzz_xxyyz_1[i] * pa_y[i];

        g_yyyzzzz_xxyzz_0[i] = 2.0 * g_yzzzz_xxyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxyzz_1[i] * fz_be_0 + g_yyzzzz_xxzz_1[i] * fe_0 + g_yyzzzz_xxyzz_1[i] * pa_y[i];

        g_yyyzzzz_xxzzz_0[i] = 2.0 * g_yzzzz_xxzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxzzz_1[i] * fz_be_0 + g_yyzzzz_xxzzz_1[i] * pa_y[i];

        g_yyyzzzz_xyyyy_0[i] = 3.0 * g_yyyzz_xyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_xyyyy_1[i] * fz_be_0 + g_yyyzzz_xyyyy_1[i] * pa_z[i];

        g_yyyzzzz_xyyyz_0[i] = 2.0 * g_yzzzz_xyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_xyyz_1[i] * fe_0 + g_yyzzzz_xyyyz_1[i] * pa_y[i];

        g_yyyzzzz_xyyzz_0[i] = 2.0 * g_yzzzz_xyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_xyzz_1[i] * fe_0 + g_yyzzzz_xyyzz_1[i] * pa_y[i];

        g_yyyzzzz_xyzzz_0[i] = 2.0 * g_yzzzz_xyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xyzzz_1[i] * fz_be_0 + g_yyzzzz_xzzz_1[i] * fe_0 + g_yyzzzz_xyzzz_1[i] * pa_y[i];

        g_yyyzzzz_xzzzz_0[i] = 2.0 * g_yzzzz_xzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xzzzz_1[i] * fz_be_0 + g_yyzzzz_xzzzz_1[i] * pa_y[i];

        g_yyyzzzz_yyyyy_0[i] = 3.0 * g_yyyzz_yyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_yyyyy_1[i] * fz_be_0 + g_yyyzzz_yyyyy_1[i] * pa_z[i];

        g_yyyzzzz_yyyyz_0[i] = 2.0 * g_yzzzz_yyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_yyyz_1[i] * fe_0 + g_yyzzzz_yyyyz_1[i] * pa_y[i];

        g_yyyzzzz_yyyzz_0[i] = 2.0 * g_yzzzz_yyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_yyzz_1[i] * fe_0 + g_yyzzzz_yyyzz_1[i] * pa_y[i];

        g_yyyzzzz_yyzzz_0[i] = 2.0 * g_yzzzz_yyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_yzzz_1[i] * fe_0 + g_yyzzzz_yyzzz_1[i] * pa_y[i];

        g_yyyzzzz_yzzzz_0[i] = 2.0 * g_yzzzz_yzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_yzzzz_1[i] * fz_be_0 + g_yyzzzz_zzzz_1[i] * fe_0 + g_yyzzzz_yzzzz_1[i] * pa_y[i];

        g_yyyzzzz_zzzzz_0[i] = 2.0 * g_yzzzz_zzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_zzzzz_1[i] * fz_be_0 + g_yyzzzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 693-714 components of targeted buffer : KH

    auto g_yyzzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 693);

    auto g_yyzzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 694);

    auto g_yyzzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 695);

    auto g_yyzzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 696);

    auto g_yyzzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 697);

    auto g_yyzzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 698);

    auto g_yyzzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 699);

    auto g_yyzzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 700);

    auto g_yyzzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 701);

    auto g_yyzzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 702);

    auto g_yyzzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 703);

    auto g_yyzzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 704);

    auto g_yyzzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 705);

    auto g_yyzzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 706);

    auto g_yyzzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 707);

    auto g_yyzzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 708);

    auto g_yyzzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 709);

    auto g_yyzzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 710);

    auto g_yyzzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 711);

    auto g_yyzzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 712);

    auto g_yyzzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 713);

    #pragma omp simd aligned(g_yyzzz_xxxxy_0, g_yyzzz_xxxxy_1, g_yyzzz_xxxyy_0, g_yyzzz_xxxyy_1, g_yyzzz_xxyyy_0, g_yyzzz_xxyyy_1, g_yyzzz_xyyyy_0, g_yyzzz_xyyyy_1, g_yyzzz_yyyyy_0, g_yyzzz_yyyyy_1, g_yyzzzz_xxxxy_1, g_yyzzzz_xxxyy_1, g_yyzzzz_xxyyy_1, g_yyzzzz_xyyyy_1, g_yyzzzz_yyyyy_1, g_yyzzzzz_xxxxx_0, g_yyzzzzz_xxxxy_0, g_yyzzzzz_xxxxz_0, g_yyzzzzz_xxxyy_0, g_yyzzzzz_xxxyz_0, g_yyzzzzz_xxxzz_0, g_yyzzzzz_xxyyy_0, g_yyzzzzz_xxyyz_0, g_yyzzzzz_xxyzz_0, g_yyzzzzz_xxzzz_0, g_yyzzzzz_xyyyy_0, g_yyzzzzz_xyyyz_0, g_yyzzzzz_xyyzz_0, g_yyzzzzz_xyzzz_0, g_yyzzzzz_xzzzz_0, g_yyzzzzz_yyyyy_0, g_yyzzzzz_yyyyz_0, g_yyzzzzz_yyyzz_0, g_yyzzzzz_yyzzz_0, g_yyzzzzz_yzzzz_0, g_yyzzzzz_zzzzz_0, g_yzzzzz_xxxxx_1, g_yzzzzz_xxxxz_1, g_yzzzzz_xxxyz_1, g_yzzzzz_xxxz_1, g_yzzzzz_xxxzz_1, g_yzzzzz_xxyyz_1, g_yzzzzz_xxyz_1, g_yzzzzz_xxyzz_1, g_yzzzzz_xxzz_1, g_yzzzzz_xxzzz_1, g_yzzzzz_xyyyz_1, g_yzzzzz_xyyz_1, g_yzzzzz_xyyzz_1, g_yzzzzz_xyzz_1, g_yzzzzz_xyzzz_1, g_yzzzzz_xzzz_1, g_yzzzzz_xzzzz_1, g_yzzzzz_yyyyz_1, g_yzzzzz_yyyz_1, g_yzzzzz_yyyzz_1, g_yzzzzz_yyzz_1, g_yzzzzz_yyzzz_1, g_yzzzzz_yzzz_1, g_yzzzzz_yzzzz_1, g_yzzzzz_zzzz_1, g_yzzzzz_zzzzz_1, g_zzzzz_xxxxx_0, g_zzzzz_xxxxx_1, g_zzzzz_xxxxz_0, g_zzzzz_xxxxz_1, g_zzzzz_xxxyz_0, g_zzzzz_xxxyz_1, g_zzzzz_xxxzz_0, g_zzzzz_xxxzz_1, g_zzzzz_xxyyz_0, g_zzzzz_xxyyz_1, g_zzzzz_xxyzz_0, g_zzzzz_xxyzz_1, g_zzzzz_xxzzz_0, g_zzzzz_xxzzz_1, g_zzzzz_xyyyz_0, g_zzzzz_xyyyz_1, g_zzzzz_xyyzz_0, g_zzzzz_xyyzz_1, g_zzzzz_xyzzz_0, g_zzzzz_xyzzz_1, g_zzzzz_xzzzz_0, g_zzzzz_xzzzz_1, g_zzzzz_yyyyz_0, g_zzzzz_yyyyz_1, g_zzzzz_yyyzz_0, g_zzzzz_yyyzz_1, g_zzzzz_yyzzz_0, g_zzzzz_yyzzz_1, g_zzzzz_yzzzz_0, g_zzzzz_yzzzz_1, g_zzzzz_zzzzz_0, g_zzzzz_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzzzz_xxxxx_0[i] = g_zzzzz_xxxxx_0[i] * fbe_0 - g_zzzzz_xxxxx_1[i] * fz_be_0 + g_yzzzzz_xxxxx_1[i] * pa_y[i];

        g_yyzzzzz_xxxxy_0[i] = 4.0 * g_yyzzz_xxxxy_0[i] * fbe_0 - 4.0 * g_yyzzz_xxxxy_1[i] * fz_be_0 + g_yyzzzz_xxxxy_1[i] * pa_z[i];

        g_yyzzzzz_xxxxz_0[i] = g_zzzzz_xxxxz_0[i] * fbe_0 - g_zzzzz_xxxxz_1[i] * fz_be_0 + g_yzzzzz_xxxxz_1[i] * pa_y[i];

        g_yyzzzzz_xxxyy_0[i] = 4.0 * g_yyzzz_xxxyy_0[i] * fbe_0 - 4.0 * g_yyzzz_xxxyy_1[i] * fz_be_0 + g_yyzzzz_xxxyy_1[i] * pa_z[i];

        g_yyzzzzz_xxxyz_0[i] = g_zzzzz_xxxyz_0[i] * fbe_0 - g_zzzzz_xxxyz_1[i] * fz_be_0 + g_yzzzzz_xxxz_1[i] * fe_0 + g_yzzzzz_xxxyz_1[i] * pa_y[i];

        g_yyzzzzz_xxxzz_0[i] = g_zzzzz_xxxzz_0[i] * fbe_0 - g_zzzzz_xxxzz_1[i] * fz_be_0 + g_yzzzzz_xxxzz_1[i] * pa_y[i];

        g_yyzzzzz_xxyyy_0[i] = 4.0 * g_yyzzz_xxyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_xxyyy_1[i] * fz_be_0 + g_yyzzzz_xxyyy_1[i] * pa_z[i];

        g_yyzzzzz_xxyyz_0[i] = g_zzzzz_xxyyz_0[i] * fbe_0 - g_zzzzz_xxyyz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_xxyz_1[i] * fe_0 + g_yzzzzz_xxyyz_1[i] * pa_y[i];

        g_yyzzzzz_xxyzz_0[i] = g_zzzzz_xxyzz_0[i] * fbe_0 - g_zzzzz_xxyzz_1[i] * fz_be_0 + g_yzzzzz_xxzz_1[i] * fe_0 + g_yzzzzz_xxyzz_1[i] * pa_y[i];

        g_yyzzzzz_xxzzz_0[i] = g_zzzzz_xxzzz_0[i] * fbe_0 - g_zzzzz_xxzzz_1[i] * fz_be_0 + g_yzzzzz_xxzzz_1[i] * pa_y[i];

        g_yyzzzzz_xyyyy_0[i] = 4.0 * g_yyzzz_xyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_xyyyy_1[i] * fz_be_0 + g_yyzzzz_xyyyy_1[i] * pa_z[i];

        g_yyzzzzz_xyyyz_0[i] = g_zzzzz_xyyyz_0[i] * fbe_0 - g_zzzzz_xyyyz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_xyyz_1[i] * fe_0 + g_yzzzzz_xyyyz_1[i] * pa_y[i];

        g_yyzzzzz_xyyzz_0[i] = g_zzzzz_xyyzz_0[i] * fbe_0 - g_zzzzz_xyyzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_xyzz_1[i] * fe_0 + g_yzzzzz_xyyzz_1[i] * pa_y[i];

        g_yyzzzzz_xyzzz_0[i] = g_zzzzz_xyzzz_0[i] * fbe_0 - g_zzzzz_xyzzz_1[i] * fz_be_0 + g_yzzzzz_xzzz_1[i] * fe_0 + g_yzzzzz_xyzzz_1[i] * pa_y[i];

        g_yyzzzzz_xzzzz_0[i] = g_zzzzz_xzzzz_0[i] * fbe_0 - g_zzzzz_xzzzz_1[i] * fz_be_0 + g_yzzzzz_xzzzz_1[i] * pa_y[i];

        g_yyzzzzz_yyyyy_0[i] = 4.0 * g_yyzzz_yyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_yyyyy_1[i] * fz_be_0 + g_yyzzzz_yyyyy_1[i] * pa_z[i];

        g_yyzzzzz_yyyyz_0[i] = g_zzzzz_yyyyz_0[i] * fbe_0 - g_zzzzz_yyyyz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_yyyz_1[i] * fe_0 + g_yzzzzz_yyyyz_1[i] * pa_y[i];

        g_yyzzzzz_yyyzz_0[i] = g_zzzzz_yyyzz_0[i] * fbe_0 - g_zzzzz_yyyzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_yyzz_1[i] * fe_0 + g_yzzzzz_yyyzz_1[i] * pa_y[i];

        g_yyzzzzz_yyzzz_0[i] = g_zzzzz_yyzzz_0[i] * fbe_0 - g_zzzzz_yyzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_yzzz_1[i] * fe_0 + g_yzzzzz_yyzzz_1[i] * pa_y[i];

        g_yyzzzzz_yzzzz_0[i] = g_zzzzz_yzzzz_0[i] * fbe_0 - g_zzzzz_yzzzz_1[i] * fz_be_0 + g_yzzzzz_zzzz_1[i] * fe_0 + g_yzzzzz_yzzzz_1[i] * pa_y[i];

        g_yyzzzzz_zzzzz_0[i] = g_zzzzz_zzzzz_0[i] * fbe_0 - g_zzzzz_zzzzz_1[i] * fz_be_0 + g_yzzzzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 714-735 components of targeted buffer : KH

    auto g_yzzzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 714);

    auto g_yzzzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 715);

    auto g_yzzzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 716);

    auto g_yzzzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 717);

    auto g_yzzzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 718);

    auto g_yzzzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 719);

    auto g_yzzzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 720);

    auto g_yzzzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 721);

    auto g_yzzzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 722);

    auto g_yzzzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 723);

    auto g_yzzzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 724);

    auto g_yzzzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 725);

    auto g_yzzzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 726);

    auto g_yzzzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 727);

    auto g_yzzzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 728);

    auto g_yzzzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 729);

    auto g_yzzzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 730);

    auto g_yzzzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 731);

    auto g_yzzzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 732);

    auto g_yzzzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 733);

    auto g_yzzzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 734);

    #pragma omp simd aligned(g_yzzzzzz_xxxxx_0, g_yzzzzzz_xxxxy_0, g_yzzzzzz_xxxxz_0, g_yzzzzzz_xxxyy_0, g_yzzzzzz_xxxyz_0, g_yzzzzzz_xxxzz_0, g_yzzzzzz_xxyyy_0, g_yzzzzzz_xxyyz_0, g_yzzzzzz_xxyzz_0, g_yzzzzzz_xxzzz_0, g_yzzzzzz_xyyyy_0, g_yzzzzzz_xyyyz_0, g_yzzzzzz_xyyzz_0, g_yzzzzzz_xyzzz_0, g_yzzzzzz_xzzzz_0, g_yzzzzzz_yyyyy_0, g_yzzzzzz_yyyyz_0, g_yzzzzzz_yyyzz_0, g_yzzzzzz_yyzzz_0, g_yzzzzzz_yzzzz_0, g_yzzzzzz_zzzzz_0, g_zzzzzz_xxxx_1, g_zzzzzz_xxxxx_1, g_zzzzzz_xxxxy_1, g_zzzzzz_xxxxz_1, g_zzzzzz_xxxy_1, g_zzzzzz_xxxyy_1, g_zzzzzz_xxxyz_1, g_zzzzzz_xxxz_1, g_zzzzzz_xxxzz_1, g_zzzzzz_xxyy_1, g_zzzzzz_xxyyy_1, g_zzzzzz_xxyyz_1, g_zzzzzz_xxyz_1, g_zzzzzz_xxyzz_1, g_zzzzzz_xxzz_1, g_zzzzzz_xxzzz_1, g_zzzzzz_xyyy_1, g_zzzzzz_xyyyy_1, g_zzzzzz_xyyyz_1, g_zzzzzz_xyyz_1, g_zzzzzz_xyyzz_1, g_zzzzzz_xyzz_1, g_zzzzzz_xyzzz_1, g_zzzzzz_xzzz_1, g_zzzzzz_xzzzz_1, g_zzzzzz_yyyy_1, g_zzzzzz_yyyyy_1, g_zzzzzz_yyyyz_1, g_zzzzzz_yyyz_1, g_zzzzzz_yyyzz_1, g_zzzzzz_yyzz_1, g_zzzzzz_yyzzz_1, g_zzzzzz_yzzz_1, g_zzzzzz_yzzzz_1, g_zzzzzz_zzzz_1, g_zzzzzz_zzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzzzz_xxxxx_0[i] = g_zzzzzz_xxxxx_1[i] * pa_y[i];

        g_yzzzzzz_xxxxy_0[i] = g_zzzzzz_xxxx_1[i] * fe_0 + g_zzzzzz_xxxxy_1[i] * pa_y[i];

        g_yzzzzzz_xxxxz_0[i] = g_zzzzzz_xxxxz_1[i] * pa_y[i];

        g_yzzzzzz_xxxyy_0[i] = 2.0 * g_zzzzzz_xxxy_1[i] * fe_0 + g_zzzzzz_xxxyy_1[i] * pa_y[i];

        g_yzzzzzz_xxxyz_0[i] = g_zzzzzz_xxxz_1[i] * fe_0 + g_zzzzzz_xxxyz_1[i] * pa_y[i];

        g_yzzzzzz_xxxzz_0[i] = g_zzzzzz_xxxzz_1[i] * pa_y[i];

        g_yzzzzzz_xxyyy_0[i] = 3.0 * g_zzzzzz_xxyy_1[i] * fe_0 + g_zzzzzz_xxyyy_1[i] * pa_y[i];

        g_yzzzzzz_xxyyz_0[i] = 2.0 * g_zzzzzz_xxyz_1[i] * fe_0 + g_zzzzzz_xxyyz_1[i] * pa_y[i];

        g_yzzzzzz_xxyzz_0[i] = g_zzzzzz_xxzz_1[i] * fe_0 + g_zzzzzz_xxyzz_1[i] * pa_y[i];

        g_yzzzzzz_xxzzz_0[i] = g_zzzzzz_xxzzz_1[i] * pa_y[i];

        g_yzzzzzz_xyyyy_0[i] = 4.0 * g_zzzzzz_xyyy_1[i] * fe_0 + g_zzzzzz_xyyyy_1[i] * pa_y[i];

        g_yzzzzzz_xyyyz_0[i] = 3.0 * g_zzzzzz_xyyz_1[i] * fe_0 + g_zzzzzz_xyyyz_1[i] * pa_y[i];

        g_yzzzzzz_xyyzz_0[i] = 2.0 * g_zzzzzz_xyzz_1[i] * fe_0 + g_zzzzzz_xyyzz_1[i] * pa_y[i];

        g_yzzzzzz_xyzzz_0[i] = g_zzzzzz_xzzz_1[i] * fe_0 + g_zzzzzz_xyzzz_1[i] * pa_y[i];

        g_yzzzzzz_xzzzz_0[i] = g_zzzzzz_xzzzz_1[i] * pa_y[i];

        g_yzzzzzz_yyyyy_0[i] = 5.0 * g_zzzzzz_yyyy_1[i] * fe_0 + g_zzzzzz_yyyyy_1[i] * pa_y[i];

        g_yzzzzzz_yyyyz_0[i] = 4.0 * g_zzzzzz_yyyz_1[i] * fe_0 + g_zzzzzz_yyyyz_1[i] * pa_y[i];

        g_yzzzzzz_yyyzz_0[i] = 3.0 * g_zzzzzz_yyzz_1[i] * fe_0 + g_zzzzzz_yyyzz_1[i] * pa_y[i];

        g_yzzzzzz_yyzzz_0[i] = 2.0 * g_zzzzzz_yzzz_1[i] * fe_0 + g_zzzzzz_yyzzz_1[i] * pa_y[i];

        g_yzzzzzz_yzzzz_0[i] = g_zzzzzz_zzzz_1[i] * fe_0 + g_zzzzzz_yzzzz_1[i] * pa_y[i];

        g_yzzzzzz_zzzzz_0[i] = g_zzzzzz_zzzzz_1[i] * pa_y[i];
    }

    // Set up 735-756 components of targeted buffer : KH

    auto g_zzzzzzz_xxxxx_0 = pbuffer.data(idx_eri_0_kh + 735);

    auto g_zzzzzzz_xxxxy_0 = pbuffer.data(idx_eri_0_kh + 736);

    auto g_zzzzzzz_xxxxz_0 = pbuffer.data(idx_eri_0_kh + 737);

    auto g_zzzzzzz_xxxyy_0 = pbuffer.data(idx_eri_0_kh + 738);

    auto g_zzzzzzz_xxxyz_0 = pbuffer.data(idx_eri_0_kh + 739);

    auto g_zzzzzzz_xxxzz_0 = pbuffer.data(idx_eri_0_kh + 740);

    auto g_zzzzzzz_xxyyy_0 = pbuffer.data(idx_eri_0_kh + 741);

    auto g_zzzzzzz_xxyyz_0 = pbuffer.data(idx_eri_0_kh + 742);

    auto g_zzzzzzz_xxyzz_0 = pbuffer.data(idx_eri_0_kh + 743);

    auto g_zzzzzzz_xxzzz_0 = pbuffer.data(idx_eri_0_kh + 744);

    auto g_zzzzzzz_xyyyy_0 = pbuffer.data(idx_eri_0_kh + 745);

    auto g_zzzzzzz_xyyyz_0 = pbuffer.data(idx_eri_0_kh + 746);

    auto g_zzzzzzz_xyyzz_0 = pbuffer.data(idx_eri_0_kh + 747);

    auto g_zzzzzzz_xyzzz_0 = pbuffer.data(idx_eri_0_kh + 748);

    auto g_zzzzzzz_xzzzz_0 = pbuffer.data(idx_eri_0_kh + 749);

    auto g_zzzzzzz_yyyyy_0 = pbuffer.data(idx_eri_0_kh + 750);

    auto g_zzzzzzz_yyyyz_0 = pbuffer.data(idx_eri_0_kh + 751);

    auto g_zzzzzzz_yyyzz_0 = pbuffer.data(idx_eri_0_kh + 752);

    auto g_zzzzzzz_yyzzz_0 = pbuffer.data(idx_eri_0_kh + 753);

    auto g_zzzzzzz_yzzzz_0 = pbuffer.data(idx_eri_0_kh + 754);

    auto g_zzzzzzz_zzzzz_0 = pbuffer.data(idx_eri_0_kh + 755);

    #pragma omp simd aligned(g_zzzzz_xxxxx_0, g_zzzzz_xxxxx_1, g_zzzzz_xxxxy_0, g_zzzzz_xxxxy_1, g_zzzzz_xxxxz_0, g_zzzzz_xxxxz_1, g_zzzzz_xxxyy_0, g_zzzzz_xxxyy_1, g_zzzzz_xxxyz_0, g_zzzzz_xxxyz_1, g_zzzzz_xxxzz_0, g_zzzzz_xxxzz_1, g_zzzzz_xxyyy_0, g_zzzzz_xxyyy_1, g_zzzzz_xxyyz_0, g_zzzzz_xxyyz_1, g_zzzzz_xxyzz_0, g_zzzzz_xxyzz_1, g_zzzzz_xxzzz_0, g_zzzzz_xxzzz_1, g_zzzzz_xyyyy_0, g_zzzzz_xyyyy_1, g_zzzzz_xyyyz_0, g_zzzzz_xyyyz_1, g_zzzzz_xyyzz_0, g_zzzzz_xyyzz_1, g_zzzzz_xyzzz_0, g_zzzzz_xyzzz_1, g_zzzzz_xzzzz_0, g_zzzzz_xzzzz_1, g_zzzzz_yyyyy_0, g_zzzzz_yyyyy_1, g_zzzzz_yyyyz_0, g_zzzzz_yyyyz_1, g_zzzzz_yyyzz_0, g_zzzzz_yyyzz_1, g_zzzzz_yyzzz_0, g_zzzzz_yyzzz_1, g_zzzzz_yzzzz_0, g_zzzzz_yzzzz_1, g_zzzzz_zzzzz_0, g_zzzzz_zzzzz_1, g_zzzzzz_xxxx_1, g_zzzzzz_xxxxx_1, g_zzzzzz_xxxxy_1, g_zzzzzz_xxxxz_1, g_zzzzzz_xxxy_1, g_zzzzzz_xxxyy_1, g_zzzzzz_xxxyz_1, g_zzzzzz_xxxz_1, g_zzzzzz_xxxzz_1, g_zzzzzz_xxyy_1, g_zzzzzz_xxyyy_1, g_zzzzzz_xxyyz_1, g_zzzzzz_xxyz_1, g_zzzzzz_xxyzz_1, g_zzzzzz_xxzz_1, g_zzzzzz_xxzzz_1, g_zzzzzz_xyyy_1, g_zzzzzz_xyyyy_1, g_zzzzzz_xyyyz_1, g_zzzzzz_xyyz_1, g_zzzzzz_xyyzz_1, g_zzzzzz_xyzz_1, g_zzzzzz_xyzzz_1, g_zzzzzz_xzzz_1, g_zzzzzz_xzzzz_1, g_zzzzzz_yyyy_1, g_zzzzzz_yyyyy_1, g_zzzzzz_yyyyz_1, g_zzzzzz_yyyz_1, g_zzzzzz_yyyzz_1, g_zzzzzz_yyzz_1, g_zzzzzz_yyzzz_1, g_zzzzzz_yzzz_1, g_zzzzzz_yzzzz_1, g_zzzzzz_zzzz_1, g_zzzzzz_zzzzz_1, g_zzzzzzz_xxxxx_0, g_zzzzzzz_xxxxy_0, g_zzzzzzz_xxxxz_0, g_zzzzzzz_xxxyy_0, g_zzzzzzz_xxxyz_0, g_zzzzzzz_xxxzz_0, g_zzzzzzz_xxyyy_0, g_zzzzzzz_xxyyz_0, g_zzzzzzz_xxyzz_0, g_zzzzzzz_xxzzz_0, g_zzzzzzz_xyyyy_0, g_zzzzzzz_xyyyz_0, g_zzzzzzz_xyyzz_0, g_zzzzzzz_xyzzz_0, g_zzzzzzz_xzzzz_0, g_zzzzzzz_yyyyy_0, g_zzzzzzz_yyyyz_0, g_zzzzzzz_yyyzz_0, g_zzzzzzz_yyzzz_0, g_zzzzzzz_yzzzz_0, g_zzzzzzz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzzzz_xxxxx_0[i] = 6.0 * g_zzzzz_xxxxx_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxxx_1[i] * fz_be_0 + g_zzzzzz_xxxxx_1[i] * pa_z[i];

        g_zzzzzzz_xxxxy_0[i] = 6.0 * g_zzzzz_xxxxy_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxxy_1[i] * fz_be_0 + g_zzzzzz_xxxxy_1[i] * pa_z[i];

        g_zzzzzzz_xxxxz_0[i] = 6.0 * g_zzzzz_xxxxz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxxz_1[i] * fz_be_0 + g_zzzzzz_xxxx_1[i] * fe_0 + g_zzzzzz_xxxxz_1[i] * pa_z[i];

        g_zzzzzzz_xxxyy_0[i] = 6.0 * g_zzzzz_xxxyy_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxyy_1[i] * fz_be_0 + g_zzzzzz_xxxyy_1[i] * pa_z[i];

        g_zzzzzzz_xxxyz_0[i] = 6.0 * g_zzzzz_xxxyz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxyz_1[i] * fz_be_0 + g_zzzzzz_xxxy_1[i] * fe_0 + g_zzzzzz_xxxyz_1[i] * pa_z[i];

        g_zzzzzzz_xxxzz_0[i] = 6.0 * g_zzzzz_xxxzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxxzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_xxxz_1[i] * fe_0 + g_zzzzzz_xxxzz_1[i] * pa_z[i];

        g_zzzzzzz_xxyyy_0[i] = 6.0 * g_zzzzz_xxyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_xxyyy_1[i] * fz_be_0 + g_zzzzzz_xxyyy_1[i] * pa_z[i];

        g_zzzzzzz_xxyyz_0[i] = 6.0 * g_zzzzz_xxyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxyyz_1[i] * fz_be_0 + g_zzzzzz_xxyy_1[i] * fe_0 + g_zzzzzz_xxyyz_1[i] * pa_z[i];

        g_zzzzzzz_xxyzz_0[i] = 6.0 * g_zzzzz_xxyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_xxyz_1[i] * fe_0 + g_zzzzzz_xxyzz_1[i] * pa_z[i];

        g_zzzzzzz_xxzzz_0[i] = 6.0 * g_zzzzz_xxzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_xxzz_1[i] * fe_0 + g_zzzzzz_xxzzz_1[i] * pa_z[i];

        g_zzzzzzz_xyyyy_0[i] = 6.0 * g_zzzzz_xyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_xyyyy_1[i] * fz_be_0 + g_zzzzzz_xyyyy_1[i] * pa_z[i];

        g_zzzzzzz_xyyyz_0[i] = 6.0 * g_zzzzz_xyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_xyyyz_1[i] * fz_be_0 + g_zzzzzz_xyyy_1[i] * fe_0 + g_zzzzzz_xyyyz_1[i] * pa_z[i];

        g_zzzzzzz_xyyzz_0[i] = 6.0 * g_zzzzz_xyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_xyyz_1[i] * fe_0 + g_zzzzzz_xyyzz_1[i] * pa_z[i];

        g_zzzzzzz_xyzzz_0[i] = 6.0 * g_zzzzz_xyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_xyzz_1[i] * fe_0 + g_zzzzzz_xyzzz_1[i] * pa_z[i];

        g_zzzzzzz_xzzzz_0[i] = 6.0 * g_zzzzz_xzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_xzzz_1[i] * fe_0 + g_zzzzzz_xzzzz_1[i] * pa_z[i];

        g_zzzzzzz_yyyyy_0[i] = 6.0 * g_zzzzz_yyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_yyyyy_1[i] * fz_be_0 + g_zzzzzz_yyyyy_1[i] * pa_z[i];

        g_zzzzzzz_yyyyz_0[i] = 6.0 * g_zzzzz_yyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_yyyyz_1[i] * fz_be_0 + g_zzzzzz_yyyy_1[i] * fe_0 + g_zzzzzz_yyyyz_1[i] * pa_z[i];

        g_zzzzzzz_yyyzz_0[i] = 6.0 * g_zzzzz_yyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_yyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_yyyz_1[i] * fe_0 + g_zzzzzz_yyyzz_1[i] * pa_z[i];

        g_zzzzzzz_yyzzz_0[i] = 6.0 * g_zzzzz_yyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_yyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_yyzz_1[i] * fe_0 + g_zzzzzz_yyzzz_1[i] * pa_z[i];

        g_zzzzzzz_yzzzz_0[i] = 6.0 * g_zzzzz_yzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_yzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_yzzz_1[i] * fe_0 + g_zzzzzz_yzzzz_1[i] * pa_z[i];

        g_zzzzzzz_zzzzz_0[i] = 6.0 * g_zzzzz_zzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_zzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_zzzz_1[i] * fe_0 + g_zzzzzz_zzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

