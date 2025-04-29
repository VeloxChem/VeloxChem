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

#include "ThreeCenterElectronRepulsionPrimRecKSH.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ksh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksh,
                                 size_t idx_eri_0_hsh,
                                 size_t idx_eri_1_hsh,
                                 size_t idx_eri_1_isg,
                                 size_t idx_eri_1_ish,
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

    /// Set up components of auxilary buffer : HSH

    auto g_xxxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh);

    auto g_xxxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 1);

    auto g_xxxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 2);

    auto g_xxxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 3);

    auto g_xxxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 4);

    auto g_xxxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 5);

    auto g_xxxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 6);

    auto g_xxxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 7);

    auto g_xxxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 8);

    auto g_xxxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 9);

    auto g_xxxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 10);

    auto g_xxxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 11);

    auto g_xxxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 12);

    auto g_xxxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 13);

    auto g_xxxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 14);

    auto g_xxxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 15);

    auto g_xxxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 16);

    auto g_xxxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 17);

    auto g_xxxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 18);

    auto g_xxxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 19);

    auto g_xxxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 20);

    auto g_xxxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 21);

    auto g_xxxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 23);

    auto g_xxxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 26);

    auto g_xxxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 30);

    auto g_xxxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 35);

    auto g_xxxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 42);

    auto g_xxxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 43);

    auto g_xxxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 45);

    auto g_xxxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 48);

    auto g_xxxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 52);

    auto g_xxxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 63);

    auto g_xxxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 64);

    auto g_xxxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 65);

    auto g_xxxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 66);

    auto g_xxxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 67);

    auto g_xxxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 68);

    auto g_xxxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 69);

    auto g_xxxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 70);

    auto g_xxxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 71);

    auto g_xxxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 72);

    auto g_xxxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 73);

    auto g_xxxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 74);

    auto g_xxxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 75);

    auto g_xxxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 76);

    auto g_xxxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 77);

    auto g_xxxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 78);

    auto g_xxxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 79);

    auto g_xxxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 80);

    auto g_xxxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 81);

    auto g_xxxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 82);

    auto g_xxxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 83);

    auto g_xxxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 105);

    auto g_xxxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 106);

    auto g_xxxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 107);

    auto g_xxxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 108);

    auto g_xxxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 109);

    auto g_xxxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 110);

    auto g_xxxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 111);

    auto g_xxxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 112);

    auto g_xxxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 113);

    auto g_xxxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 114);

    auto g_xxxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 115);

    auto g_xxxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 116);

    auto g_xxxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 117);

    auto g_xxxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 118);

    auto g_xxxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 119);

    auto g_xxxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 120);

    auto g_xxxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 121);

    auto g_xxxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 122);

    auto g_xxxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 123);

    auto g_xxxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 124);

    auto g_xxxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 125);

    auto g_xxyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 126);

    auto g_xxyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 127);

    auto g_xxyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 128);

    auto g_xxyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 129);

    auto g_xxyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 130);

    auto g_xxyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 131);

    auto g_xxyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 132);

    auto g_xxyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 133);

    auto g_xxyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 134);

    auto g_xxyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 135);

    auto g_xxyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 136);

    auto g_xxyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 137);

    auto g_xxyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 138);

    auto g_xxyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 139);

    auto g_xxyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 140);

    auto g_xxyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 141);

    auto g_xxyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 142);

    auto g_xxyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 143);

    auto g_xxyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 144);

    auto g_xxyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 145);

    auto g_xxyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 146);

    auto g_xxyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 148);

    auto g_xxyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 150);

    auto g_xxyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 153);

    auto g_xxyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 157);

    auto g_xxyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 168);

    auto g_xxyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 170);

    auto g_xxyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 173);

    auto g_xxyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 177);

    auto g_xxyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 182);

    auto g_xxzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 189);

    auto g_xxzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 190);

    auto g_xxzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 191);

    auto g_xxzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 192);

    auto g_xxzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 193);

    auto g_xxzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 194);

    auto g_xxzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 195);

    auto g_xxzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 196);

    auto g_xxzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 197);

    auto g_xxzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 198);

    auto g_xxzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 199);

    auto g_xxzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 200);

    auto g_xxzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 201);

    auto g_xxzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 202);

    auto g_xxzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 203);

    auto g_xxzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 204);

    auto g_xxzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 205);

    auto g_xxzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 206);

    auto g_xxzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 207);

    auto g_xxzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 208);

    auto g_xxzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 209);

    auto g_xyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 211);

    auto g_xyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 213);

    auto g_xyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 214);

    auto g_xyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 216);

    auto g_xyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 217);

    auto g_xyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 218);

    auto g_xyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 220);

    auto g_xyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 221);

    auto g_xyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 222);

    auto g_xyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 223);

    auto g_xyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 225);

    auto g_xyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 226);

    auto g_xyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 227);

    auto g_xyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 228);

    auto g_xyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 229);

    auto g_xyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 230);

    auto g_xyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 256);

    auto g_xyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 259);

    auto g_xyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 260);

    auto g_xyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 263);

    auto g_xyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 264);

    auto g_xyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 265);

    auto g_xyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 267);

    auto g_xyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 268);

    auto g_xyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 269);

    auto g_xyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 270);

    auto g_xyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 271);

    auto g_xyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 272);

    auto g_xzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 296);

    auto g_xzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 298);

    auto g_xzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 299);

    auto g_xzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 301);

    auto g_xzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 302);

    auto g_xzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 303);

    auto g_xzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 305);

    auto g_xzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 306);

    auto g_xzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 307);

    auto g_xzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 308);

    auto g_xzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 309);

    auto g_xzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 310);

    auto g_xzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 311);

    auto g_xzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 312);

    auto g_xzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 313);

    auto g_xzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 314);

    auto g_yyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 315);

    auto g_yyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 316);

    auto g_yyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 317);

    auto g_yyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 318);

    auto g_yyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 319);

    auto g_yyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 320);

    auto g_yyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 321);

    auto g_yyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 322);

    auto g_yyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 323);

    auto g_yyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 324);

    auto g_yyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 325);

    auto g_yyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 326);

    auto g_yyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 327);

    auto g_yyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 328);

    auto g_yyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 329);

    auto g_yyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 330);

    auto g_yyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 331);

    auto g_yyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 332);

    auto g_yyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 333);

    auto g_yyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 334);

    auto g_yyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 335);

    auto g_yyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 337);

    auto g_yyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 339);

    auto g_yyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 342);

    auto g_yyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 346);

    auto g_yyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 351);

    auto g_yyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 357);

    auto g_yyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 358);

    auto g_yyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 359);

    auto g_yyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 360);

    auto g_yyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 361);

    auto g_yyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 362);

    auto g_yyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 363);

    auto g_yyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 364);

    auto g_yyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 365);

    auto g_yyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 366);

    auto g_yyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 367);

    auto g_yyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 368);

    auto g_yyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 369);

    auto g_yyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 370);

    auto g_yyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 371);

    auto g_yyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 372);

    auto g_yyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 373);

    auto g_yyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 374);

    auto g_yyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 375);

    auto g_yyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 376);

    auto g_yyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 377);

    auto g_yyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 378);

    auto g_yyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 379);

    auto g_yyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 380);

    auto g_yyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 381);

    auto g_yyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 382);

    auto g_yyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 383);

    auto g_yyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 384);

    auto g_yyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 385);

    auto g_yyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 386);

    auto g_yyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 387);

    auto g_yyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 388);

    auto g_yyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 389);

    auto g_yyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 390);

    auto g_yyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 391);

    auto g_yyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 392);

    auto g_yyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 393);

    auto g_yyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 394);

    auto g_yyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 395);

    auto g_yyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 396);

    auto g_yyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 397);

    auto g_yyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 398);

    auto g_yzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 399);

    auto g_yzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 401);

    auto g_yzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 403);

    auto g_yzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 404);

    auto g_yzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 406);

    auto g_yzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 407);

    auto g_yzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 408);

    auto g_yzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 410);

    auto g_yzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 411);

    auto g_yzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 412);

    auto g_yzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 413);

    auto g_yzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 415);

    auto g_yzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 416);

    auto g_yzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 417);

    auto g_yzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 418);

    auto g_yzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 419);

    auto g_zzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_hsh + 420);

    auto g_zzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_hsh + 421);

    auto g_zzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_hsh + 422);

    auto g_zzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_hsh + 423);

    auto g_zzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_hsh + 424);

    auto g_zzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_hsh + 425);

    auto g_zzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_hsh + 426);

    auto g_zzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_hsh + 427);

    auto g_zzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_hsh + 428);

    auto g_zzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_hsh + 429);

    auto g_zzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_hsh + 430);

    auto g_zzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_hsh + 431);

    auto g_zzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_hsh + 432);

    auto g_zzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_hsh + 433);

    auto g_zzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_hsh + 434);

    auto g_zzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_hsh + 435);

    auto g_zzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_hsh + 436);

    auto g_zzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_hsh + 437);

    auto g_zzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_hsh + 438);

    auto g_zzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_hsh + 439);

    auto g_zzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_hsh + 440);

    /// Set up components of auxilary buffer : HSH

    auto g_xxxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh);

    auto g_xxxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 1);

    auto g_xxxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 2);

    auto g_xxxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 3);

    auto g_xxxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 4);

    auto g_xxxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 5);

    auto g_xxxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 6);

    auto g_xxxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 7);

    auto g_xxxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 8);

    auto g_xxxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 9);

    auto g_xxxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 10);

    auto g_xxxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 11);

    auto g_xxxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 12);

    auto g_xxxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 13);

    auto g_xxxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 14);

    auto g_xxxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 15);

    auto g_xxxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 16);

    auto g_xxxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 17);

    auto g_xxxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 18);

    auto g_xxxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 19);

    auto g_xxxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 20);

    auto g_xxxxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 21);

    auto g_xxxxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 23);

    auto g_xxxxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 26);

    auto g_xxxxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 30);

    auto g_xxxxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 35);

    auto g_xxxxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 42);

    auto g_xxxxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 43);

    auto g_xxxxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 45);

    auto g_xxxxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 48);

    auto g_xxxxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 52);

    auto g_xxxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 63);

    auto g_xxxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 64);

    auto g_xxxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 65);

    auto g_xxxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 66);

    auto g_xxxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 67);

    auto g_xxxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 68);

    auto g_xxxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 69);

    auto g_xxxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 70);

    auto g_xxxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 71);

    auto g_xxxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 72);

    auto g_xxxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 73);

    auto g_xxxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 74);

    auto g_xxxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 75);

    auto g_xxxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 76);

    auto g_xxxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 77);

    auto g_xxxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 78);

    auto g_xxxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 79);

    auto g_xxxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 80);

    auto g_xxxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 81);

    auto g_xxxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 82);

    auto g_xxxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 83);

    auto g_xxxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 105);

    auto g_xxxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 106);

    auto g_xxxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 107);

    auto g_xxxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 108);

    auto g_xxxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 109);

    auto g_xxxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 110);

    auto g_xxxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 111);

    auto g_xxxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 112);

    auto g_xxxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 113);

    auto g_xxxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 114);

    auto g_xxxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 115);

    auto g_xxxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 116);

    auto g_xxxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 117);

    auto g_xxxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 118);

    auto g_xxxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 119);

    auto g_xxxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 120);

    auto g_xxxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 121);

    auto g_xxxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 122);

    auto g_xxxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 123);

    auto g_xxxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 124);

    auto g_xxxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 125);

    auto g_xxyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 126);

    auto g_xxyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 127);

    auto g_xxyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 128);

    auto g_xxyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 129);

    auto g_xxyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 130);

    auto g_xxyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 131);

    auto g_xxyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 132);

    auto g_xxyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 133);

    auto g_xxyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 134);

    auto g_xxyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 135);

    auto g_xxyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 136);

    auto g_xxyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 137);

    auto g_xxyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 138);

    auto g_xxyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 139);

    auto g_xxyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 140);

    auto g_xxyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 141);

    auto g_xxyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 142);

    auto g_xxyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 143);

    auto g_xxyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 144);

    auto g_xxyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 145);

    auto g_xxyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 146);

    auto g_xxyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 148);

    auto g_xxyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 150);

    auto g_xxyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 153);

    auto g_xxyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 157);

    auto g_xxyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 168);

    auto g_xxyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 170);

    auto g_xxyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 173);

    auto g_xxyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 177);

    auto g_xxyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 182);

    auto g_xxzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 189);

    auto g_xxzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 190);

    auto g_xxzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 191);

    auto g_xxzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 192);

    auto g_xxzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 193);

    auto g_xxzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 194);

    auto g_xxzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 195);

    auto g_xxzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 196);

    auto g_xxzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 197);

    auto g_xxzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 198);

    auto g_xxzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 199);

    auto g_xxzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 200);

    auto g_xxzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 201);

    auto g_xxzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 202);

    auto g_xxzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 203);

    auto g_xxzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 204);

    auto g_xxzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 205);

    auto g_xxzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 206);

    auto g_xxzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 207);

    auto g_xxzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 208);

    auto g_xxzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 209);

    auto g_xyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 211);

    auto g_xyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 213);

    auto g_xyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 214);

    auto g_xyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 216);

    auto g_xyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 217);

    auto g_xyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 218);

    auto g_xyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 220);

    auto g_xyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 221);

    auto g_xyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 222);

    auto g_xyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 223);

    auto g_xyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 225);

    auto g_xyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 226);

    auto g_xyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 227);

    auto g_xyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 228);

    auto g_xyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 229);

    auto g_xyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 230);

    auto g_xyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 256);

    auto g_xyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 259);

    auto g_xyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 260);

    auto g_xyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 263);

    auto g_xyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 264);

    auto g_xyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 265);

    auto g_xyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 267);

    auto g_xyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 268);

    auto g_xyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 269);

    auto g_xyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 270);

    auto g_xyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 271);

    auto g_xyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 272);

    auto g_xzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 296);

    auto g_xzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 298);

    auto g_xzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 299);

    auto g_xzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 301);

    auto g_xzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 302);

    auto g_xzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 303);

    auto g_xzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 305);

    auto g_xzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 306);

    auto g_xzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 307);

    auto g_xzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 308);

    auto g_xzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 309);

    auto g_xzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 310);

    auto g_xzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 311);

    auto g_xzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 312);

    auto g_xzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 313);

    auto g_xzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 314);

    auto g_yyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 315);

    auto g_yyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 316);

    auto g_yyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 317);

    auto g_yyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 318);

    auto g_yyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 319);

    auto g_yyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 320);

    auto g_yyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 321);

    auto g_yyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 322);

    auto g_yyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 323);

    auto g_yyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 324);

    auto g_yyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 325);

    auto g_yyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 326);

    auto g_yyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 327);

    auto g_yyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 328);

    auto g_yyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 329);

    auto g_yyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 330);

    auto g_yyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 331);

    auto g_yyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 332);

    auto g_yyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 333);

    auto g_yyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 334);

    auto g_yyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 335);

    auto g_yyyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 337);

    auto g_yyyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 339);

    auto g_yyyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 342);

    auto g_yyyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 346);

    auto g_yyyyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 351);

    auto g_yyyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 357);

    auto g_yyyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 358);

    auto g_yyyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 359);

    auto g_yyyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 360);

    auto g_yyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 361);

    auto g_yyyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 362);

    auto g_yyyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 363);

    auto g_yyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 364);

    auto g_yyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 365);

    auto g_yyyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 366);

    auto g_yyyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 367);

    auto g_yyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 368);

    auto g_yyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 369);

    auto g_yyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 370);

    auto g_yyyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 371);

    auto g_yyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 372);

    auto g_yyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 373);

    auto g_yyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 374);

    auto g_yyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 375);

    auto g_yyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 376);

    auto g_yyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 377);

    auto g_yyzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 378);

    auto g_yyzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 379);

    auto g_yyzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 380);

    auto g_yyzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 381);

    auto g_yyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 382);

    auto g_yyzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 383);

    auto g_yyzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 384);

    auto g_yyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 385);

    auto g_yyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 386);

    auto g_yyzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 387);

    auto g_yyzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 388);

    auto g_yyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 389);

    auto g_yyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 390);

    auto g_yyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 391);

    auto g_yyzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 392);

    auto g_yyzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 393);

    auto g_yyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 394);

    auto g_yyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 395);

    auto g_yyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 396);

    auto g_yyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 397);

    auto g_yyzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 398);

    auto g_yzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 399);

    auto g_yzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 401);

    auto g_yzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 403);

    auto g_yzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 404);

    auto g_yzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 406);

    auto g_yzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 407);

    auto g_yzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 408);

    auto g_yzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 410);

    auto g_yzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 411);

    auto g_yzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 412);

    auto g_yzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 413);

    auto g_yzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 415);

    auto g_yzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 416);

    auto g_yzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 417);

    auto g_yzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 418);

    auto g_yzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 419);

    auto g_zzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_hsh + 420);

    auto g_zzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 421);

    auto g_zzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 422);

    auto g_zzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 423);

    auto g_zzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 424);

    auto g_zzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 425);

    auto g_zzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 426);

    auto g_zzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 427);

    auto g_zzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 428);

    auto g_zzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 429);

    auto g_zzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 430);

    auto g_zzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 431);

    auto g_zzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 432);

    auto g_zzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 433);

    auto g_zzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 434);

    auto g_zzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 435);

    auto g_zzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 436);

    auto g_zzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 437);

    auto g_zzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 438);

    auto g_zzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 439);

    auto g_zzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 440);

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

    auto g_xxxxxz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 32);

    auto g_xxxxxz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 34);

    auto g_xxxxxz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 35);

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

    auto g_xxyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 184);

    auto g_xxyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 187);

    auto g_xxyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 188);

    auto g_xxyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 191);

    auto g_xxyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 192);

    auto g_xxyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 193);

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

    auto g_xyyyzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 259);

    auto g_xyyyzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 262);

    auto g_xyyyzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 263);

    auto g_xyyyzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 266);

    auto g_xyyyzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 267);

    auto g_xyyyzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 268);

    auto g_xyyzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 274);

    auto g_xyyzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 277);

    auto g_xyyzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 278);

    auto g_xyyzzz_0_yyyz_1 = pbuffer.data(idx_eri_1_isg + 281);

    auto g_xyyzzz_0_yyzz_1 = pbuffer.data(idx_eri_1_isg + 282);

    auto g_xyyzzz_0_yzzz_1 = pbuffer.data(idx_eri_1_isg + 283);

    auto g_xzzzzz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 302);

    auto g_xzzzzz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 304);

    auto g_xzzzzz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 305);

    auto g_xzzzzz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 307);

    auto g_xzzzzz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 308);

    auto g_xzzzzz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 309);

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

    auto g_yyyyyz_0_xxxz_1 = pbuffer.data(idx_eri_1_isg + 332);

    auto g_yyyyyz_0_xxyz_1 = pbuffer.data(idx_eri_1_isg + 334);

    auto g_yyyyyz_0_xxzz_1 = pbuffer.data(idx_eri_1_isg + 335);

    auto g_yyyyyz_0_xyyz_1 = pbuffer.data(idx_eri_1_isg + 337);

    auto g_yyyyyz_0_xyzz_1 = pbuffer.data(idx_eri_1_isg + 338);

    auto g_yyyyyz_0_xzzz_1 = pbuffer.data(idx_eri_1_isg + 339);

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

    /// Set up components of auxilary buffer : ISH

    auto g_xxxxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish);

    auto g_xxxxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 1);

    auto g_xxxxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 2);

    auto g_xxxxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 3);

    auto g_xxxxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 4);

    auto g_xxxxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 5);

    auto g_xxxxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 6);

    auto g_xxxxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 7);

    auto g_xxxxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 8);

    auto g_xxxxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 9);

    auto g_xxxxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 10);

    auto g_xxxxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 11);

    auto g_xxxxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 12);

    auto g_xxxxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 13);

    auto g_xxxxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 14);

    auto g_xxxxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 15);

    auto g_xxxxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 16);

    auto g_xxxxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 17);

    auto g_xxxxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 18);

    auto g_xxxxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 19);

    auto g_xxxxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 20);

    auto g_xxxxxy_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 21);

    auto g_xxxxxy_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 22);

    auto g_xxxxxy_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 23);

    auto g_xxxxxy_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 24);

    auto g_xxxxxy_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 26);

    auto g_xxxxxy_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 27);

    auto g_xxxxxy_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 30);

    auto g_xxxxxy_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 31);

    auto g_xxxxxy_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 35);

    auto g_xxxxxy_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 36);

    auto g_xxxxxz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 42);

    auto g_xxxxxz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 43);

    auto g_xxxxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 44);

    auto g_xxxxxz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 45);

    auto g_xxxxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 46);

    auto g_xxxxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 47);

    auto g_xxxxxz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 48);

    auto g_xxxxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 49);

    auto g_xxxxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 50);

    auto g_xxxxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 51);

    auto g_xxxxxz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 52);

    auto g_xxxxxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 53);

    auto g_xxxxxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 54);

    auto g_xxxxxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 55);

    auto g_xxxxxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 56);

    auto g_xxxxxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 58);

    auto g_xxxxxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 59);

    auto g_xxxxxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 60);

    auto g_xxxxxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 61);

    auto g_xxxxxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 62);

    auto g_xxxxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 63);

    auto g_xxxxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 64);

    auto g_xxxxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 65);

    auto g_xxxxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 66);

    auto g_xxxxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 67);

    auto g_xxxxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 68);

    auto g_xxxxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 69);

    auto g_xxxxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 70);

    auto g_xxxxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 71);

    auto g_xxxxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 72);

    auto g_xxxxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 73);

    auto g_xxxxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 74);

    auto g_xxxxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 75);

    auto g_xxxxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 76);

    auto g_xxxxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 77);

    auto g_xxxxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 78);

    auto g_xxxxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 79);

    auto g_xxxxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 80);

    auto g_xxxxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 81);

    auto g_xxxxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 82);

    auto g_xxxxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 83);

    auto g_xxxxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 105);

    auto g_xxxxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 106);

    auto g_xxxxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 107);

    auto g_xxxxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 108);

    auto g_xxxxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 109);

    auto g_xxxxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 110);

    auto g_xxxxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 111);

    auto g_xxxxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 112);

    auto g_xxxxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 113);

    auto g_xxxxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 114);

    auto g_xxxxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 115);

    auto g_xxxxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 116);

    auto g_xxxxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 117);

    auto g_xxxxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 118);

    auto g_xxxxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 119);

    auto g_xxxxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 120);

    auto g_xxxxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 121);

    auto g_xxxxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 122);

    auto g_xxxxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 123);

    auto g_xxxxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 124);

    auto g_xxxxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 125);

    auto g_xxxyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 126);

    auto g_xxxyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 127);

    auto g_xxxyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 128);

    auto g_xxxyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 129);

    auto g_xxxyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 130);

    auto g_xxxyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 131);

    auto g_xxxyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 132);

    auto g_xxxyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 133);

    auto g_xxxyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 134);

    auto g_xxxyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 135);

    auto g_xxxyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 136);

    auto g_xxxyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 137);

    auto g_xxxyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 138);

    auto g_xxxyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 139);

    auto g_xxxyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 140);

    auto g_xxxyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 141);

    auto g_xxxyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 142);

    auto g_xxxyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 143);

    auto g_xxxyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 144);

    auto g_xxxyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 145);

    auto g_xxxyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 146);

    auto g_xxxyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 148);

    auto g_xxxyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 150);

    auto g_xxxyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 153);

    auto g_xxxyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 157);

    auto g_xxxyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 168);

    auto g_xxxyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 170);

    auto g_xxxyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 173);

    auto g_xxxyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 177);

    auto g_xxxyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 182);

    auto g_xxxzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 189);

    auto g_xxxzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 190);

    auto g_xxxzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 191);

    auto g_xxxzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 192);

    auto g_xxxzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 193);

    auto g_xxxzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 194);

    auto g_xxxzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 195);

    auto g_xxxzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 196);

    auto g_xxxzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 197);

    auto g_xxxzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 198);

    auto g_xxxzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 199);

    auto g_xxxzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 200);

    auto g_xxxzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 201);

    auto g_xxxzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 202);

    auto g_xxxzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 203);

    auto g_xxxzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 204);

    auto g_xxxzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 205);

    auto g_xxxzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 206);

    auto g_xxxzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 207);

    auto g_xxxzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 208);

    auto g_xxxzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 209);

    auto g_xxyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 210);

    auto g_xxyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 211);

    auto g_xxyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 212);

    auto g_xxyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 213);

    auto g_xxyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 214);

    auto g_xxyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 215);

    auto g_xxyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 216);

    auto g_xxyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 217);

    auto g_xxyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 218);

    auto g_xxyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 219);

    auto g_xxyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 220);

    auto g_xxyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 221);

    auto g_xxyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 222);

    auto g_xxyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 223);

    auto g_xxyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 224);

    auto g_xxyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 225);

    auto g_xxyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 226);

    auto g_xxyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 227);

    auto g_xxyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 228);

    auto g_xxyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 229);

    auto g_xxyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 230);

    auto g_xxyyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 232);

    auto g_xxyyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 234);

    auto g_xxyyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 237);

    auto g_xxyyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 241);

    auto g_xxyyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 252);

    auto g_xxyyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 253);

    auto g_xxyyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 254);

    auto g_xxyyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 255);

    auto g_xxyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 256);

    auto g_xxyyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 257);

    auto g_xxyyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 258);

    auto g_xxyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 259);

    auto g_xxyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 260);

    auto g_xxyyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 261);

    auto g_xxyyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 262);

    auto g_xxyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 263);

    auto g_xxyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 264);

    auto g_xxyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 265);

    auto g_xxyyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 266);

    auto g_xxyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 267);

    auto g_xxyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 268);

    auto g_xxyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 269);

    auto g_xxyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 270);

    auto g_xxyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 271);

    auto g_xxyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 272);

    auto g_xxyzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 273);

    auto g_xxyzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 275);

    auto g_xxyzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 278);

    auto g_xxyzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 282);

    auto g_xxyzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 287);

    auto g_xxzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 294);

    auto g_xxzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 295);

    auto g_xxzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 296);

    auto g_xxzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 297);

    auto g_xxzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 298);

    auto g_xxzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 299);

    auto g_xxzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 300);

    auto g_xxzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 301);

    auto g_xxzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 302);

    auto g_xxzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 303);

    auto g_xxzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 304);

    auto g_xxzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 305);

    auto g_xxzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 306);

    auto g_xxzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 307);

    auto g_xxzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 308);

    auto g_xxzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 309);

    auto g_xxzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 310);

    auto g_xxzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 311);

    auto g_xxzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 312);

    auto g_xxzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 313);

    auto g_xxzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 314);

    auto g_xyyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 315);

    auto g_xyyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 316);

    auto g_xyyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 318);

    auto g_xyyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 319);

    auto g_xyyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 321);

    auto g_xyyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 322);

    auto g_xyyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 323);

    auto g_xyyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 325);

    auto g_xyyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 326);

    auto g_xyyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 327);

    auto g_xyyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 328);

    auto g_xyyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 330);

    auto g_xyyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 331);

    auto g_xyyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 332);

    auto g_xyyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 333);

    auto g_xyyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 334);

    auto g_xyyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 335);

    auto g_xyyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 361);

    auto g_xyyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 364);

    auto g_xyyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 365);

    auto g_xyyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 368);

    auto g_xyyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 369);

    auto g_xyyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 370);

    auto g_xyyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 372);

    auto g_xyyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 373);

    auto g_xyyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 374);

    auto g_xyyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 375);

    auto g_xyyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 376);

    auto g_xyyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 377);

    auto g_xyyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 382);

    auto g_xyyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 385);

    auto g_xyyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 386);

    auto g_xyyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 389);

    auto g_xyyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 390);

    auto g_xyyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 391);

    auto g_xyyzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 393);

    auto g_xyyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 394);

    auto g_xyyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 395);

    auto g_xyyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 396);

    auto g_xyyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 397);

    auto g_xyyzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 398);

    auto g_xzzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 420);

    auto g_xzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 422);

    auto g_xzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 424);

    auto g_xzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 425);

    auto g_xzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 427);

    auto g_xzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 428);

    auto g_xzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 429);

    auto g_xzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 431);

    auto g_xzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 432);

    auto g_xzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 433);

    auto g_xzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 434);

    auto g_xzzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 435);

    auto g_xzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 436);

    auto g_xzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 437);

    auto g_xzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 438);

    auto g_xzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 439);

    auto g_xzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 440);

    auto g_yyyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 441);

    auto g_yyyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 442);

    auto g_yyyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 443);

    auto g_yyyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 444);

    auto g_yyyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 445);

    auto g_yyyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 446);

    auto g_yyyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 447);

    auto g_yyyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 448);

    auto g_yyyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 449);

    auto g_yyyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 450);

    auto g_yyyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 451);

    auto g_yyyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 452);

    auto g_yyyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 453);

    auto g_yyyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 454);

    auto g_yyyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 455);

    auto g_yyyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 456);

    auto g_yyyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 457);

    auto g_yyyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 458);

    auto g_yyyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 459);

    auto g_yyyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 460);

    auto g_yyyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 461);

    auto g_yyyyyz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 463);

    auto g_yyyyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 464);

    auto g_yyyyyz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 465);

    auto g_yyyyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 466);

    auto g_yyyyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 467);

    auto g_yyyyyz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 468);

    auto g_yyyyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 469);

    auto g_yyyyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 470);

    auto g_yyyyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 471);

    auto g_yyyyyz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 472);

    auto g_yyyyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 473);

    auto g_yyyyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 474);

    auto g_yyyyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 475);

    auto g_yyyyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 476);

    auto g_yyyyyz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 477);

    auto g_yyyyyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 478);

    auto g_yyyyyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 479);

    auto g_yyyyyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 480);

    auto g_yyyyyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 481);

    auto g_yyyyyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 482);

    auto g_yyyyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 483);

    auto g_yyyyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 484);

    auto g_yyyyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 485);

    auto g_yyyyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 486);

    auto g_yyyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 487);

    auto g_yyyyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 488);

    auto g_yyyyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 489);

    auto g_yyyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 490);

    auto g_yyyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 491);

    auto g_yyyyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 492);

    auto g_yyyyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 493);

    auto g_yyyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 494);

    auto g_yyyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 495);

    auto g_yyyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 496);

    auto g_yyyyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 497);

    auto g_yyyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 498);

    auto g_yyyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 499);

    auto g_yyyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 500);

    auto g_yyyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 501);

    auto g_yyyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 502);

    auto g_yyyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 503);

    auto g_yyyzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 504);

    auto g_yyyzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 505);

    auto g_yyyzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 506);

    auto g_yyyzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 507);

    auto g_yyyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 508);

    auto g_yyyzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 509);

    auto g_yyyzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 510);

    auto g_yyyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 511);

    auto g_yyyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 512);

    auto g_yyyzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 513);

    auto g_yyyzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 514);

    auto g_yyyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 515);

    auto g_yyyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 516);

    auto g_yyyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 517);

    auto g_yyyzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 518);

    auto g_yyyzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 519);

    auto g_yyyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 520);

    auto g_yyyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 521);

    auto g_yyyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 522);

    auto g_yyyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 523);

    auto g_yyyzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 524);

    auto g_yyzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 525);

    auto g_yyzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 526);

    auto g_yyzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 527);

    auto g_yyzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 528);

    auto g_yyzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 529);

    auto g_yyzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 530);

    auto g_yyzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 531);

    auto g_yyzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 532);

    auto g_yyzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 533);

    auto g_yyzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 534);

    auto g_yyzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 535);

    auto g_yyzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 536);

    auto g_yyzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 537);

    auto g_yyzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 538);

    auto g_yyzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 539);

    auto g_yyzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 540);

    auto g_yyzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 541);

    auto g_yyzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 542);

    auto g_yyzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 543);

    auto g_yyzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 544);

    auto g_yyzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 545);

    auto g_yzzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 546);

    auto g_yzzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 547);

    auto g_yzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 548);

    auto g_yzzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 549);

    auto g_yzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 550);

    auto g_yzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 551);

    auto g_yzzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 552);

    auto g_yzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 553);

    auto g_yzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 554);

    auto g_yzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 555);

    auto g_yzzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 556);

    auto g_yzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 557);

    auto g_yzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 558);

    auto g_yzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 559);

    auto g_yzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 560);

    auto g_yzzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 561);

    auto g_yzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 562);

    auto g_yzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 563);

    auto g_yzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 564);

    auto g_yzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 565);

    auto g_yzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 566);

    auto g_zzzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_ish + 567);

    auto g_zzzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_ish + 568);

    auto g_zzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 569);

    auto g_zzzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_ish + 570);

    auto g_zzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 571);

    auto g_zzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 572);

    auto g_zzzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_ish + 573);

    auto g_zzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 574);

    auto g_zzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 575);

    auto g_zzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 576);

    auto g_zzzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_ish + 577);

    auto g_zzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 578);

    auto g_zzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 579);

    auto g_zzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 580);

    auto g_zzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 581);

    auto g_zzzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_ish + 582);

    auto g_zzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 583);

    auto g_zzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 584);

    auto g_zzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 585);

    auto g_zzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 586);

    auto g_zzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_ish + 587);

    /// Set up 0-21 components of targeted buffer : KSH

    auto g_xxxxxxx_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh);

    auto g_xxxxxxx_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 1);

    auto g_xxxxxxx_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 2);

    auto g_xxxxxxx_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 3);

    auto g_xxxxxxx_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 4);

    auto g_xxxxxxx_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 5);

    auto g_xxxxxxx_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 6);

    auto g_xxxxxxx_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 7);

    auto g_xxxxxxx_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 8);

    auto g_xxxxxxx_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 9);

    auto g_xxxxxxx_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 10);

    auto g_xxxxxxx_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 11);

    auto g_xxxxxxx_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 12);

    auto g_xxxxxxx_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 13);

    auto g_xxxxxxx_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 14);

    auto g_xxxxxxx_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 15);

    auto g_xxxxxxx_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 16);

    auto g_xxxxxxx_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 17);

    auto g_xxxxxxx_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 18);

    auto g_xxxxxxx_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 19);

    auto g_xxxxxxx_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 20);

    #pragma omp simd aligned(g_xxxxx_0_xxxxx_0, g_xxxxx_0_xxxxx_1, g_xxxxx_0_xxxxy_0, g_xxxxx_0_xxxxy_1, g_xxxxx_0_xxxxz_0, g_xxxxx_0_xxxxz_1, g_xxxxx_0_xxxyy_0, g_xxxxx_0_xxxyy_1, g_xxxxx_0_xxxyz_0, g_xxxxx_0_xxxyz_1, g_xxxxx_0_xxxzz_0, g_xxxxx_0_xxxzz_1, g_xxxxx_0_xxyyy_0, g_xxxxx_0_xxyyy_1, g_xxxxx_0_xxyyz_0, g_xxxxx_0_xxyyz_1, g_xxxxx_0_xxyzz_0, g_xxxxx_0_xxyzz_1, g_xxxxx_0_xxzzz_0, g_xxxxx_0_xxzzz_1, g_xxxxx_0_xyyyy_0, g_xxxxx_0_xyyyy_1, g_xxxxx_0_xyyyz_0, g_xxxxx_0_xyyyz_1, g_xxxxx_0_xyyzz_0, g_xxxxx_0_xyyzz_1, g_xxxxx_0_xyzzz_0, g_xxxxx_0_xyzzz_1, g_xxxxx_0_xzzzz_0, g_xxxxx_0_xzzzz_1, g_xxxxx_0_yyyyy_0, g_xxxxx_0_yyyyy_1, g_xxxxx_0_yyyyz_0, g_xxxxx_0_yyyyz_1, g_xxxxx_0_yyyzz_0, g_xxxxx_0_yyyzz_1, g_xxxxx_0_yyzzz_0, g_xxxxx_0_yyzzz_1, g_xxxxx_0_yzzzz_0, g_xxxxx_0_yzzzz_1, g_xxxxx_0_zzzzz_0, g_xxxxx_0_zzzzz_1, g_xxxxxx_0_xxxx_1, g_xxxxxx_0_xxxxx_1, g_xxxxxx_0_xxxxy_1, g_xxxxxx_0_xxxxz_1, g_xxxxxx_0_xxxy_1, g_xxxxxx_0_xxxyy_1, g_xxxxxx_0_xxxyz_1, g_xxxxxx_0_xxxz_1, g_xxxxxx_0_xxxzz_1, g_xxxxxx_0_xxyy_1, g_xxxxxx_0_xxyyy_1, g_xxxxxx_0_xxyyz_1, g_xxxxxx_0_xxyz_1, g_xxxxxx_0_xxyzz_1, g_xxxxxx_0_xxzz_1, g_xxxxxx_0_xxzzz_1, g_xxxxxx_0_xyyy_1, g_xxxxxx_0_xyyyy_1, g_xxxxxx_0_xyyyz_1, g_xxxxxx_0_xyyz_1, g_xxxxxx_0_xyyzz_1, g_xxxxxx_0_xyzz_1, g_xxxxxx_0_xyzzz_1, g_xxxxxx_0_xzzz_1, g_xxxxxx_0_xzzzz_1, g_xxxxxx_0_yyyy_1, g_xxxxxx_0_yyyyy_1, g_xxxxxx_0_yyyyz_1, g_xxxxxx_0_yyyz_1, g_xxxxxx_0_yyyzz_1, g_xxxxxx_0_yyzz_1, g_xxxxxx_0_yyzzz_1, g_xxxxxx_0_yzzz_1, g_xxxxxx_0_yzzzz_1, g_xxxxxx_0_zzzz_1, g_xxxxxx_0_zzzzz_1, g_xxxxxxx_0_xxxxx_0, g_xxxxxxx_0_xxxxy_0, g_xxxxxxx_0_xxxxz_0, g_xxxxxxx_0_xxxyy_0, g_xxxxxxx_0_xxxyz_0, g_xxxxxxx_0_xxxzz_0, g_xxxxxxx_0_xxyyy_0, g_xxxxxxx_0_xxyyz_0, g_xxxxxxx_0_xxyzz_0, g_xxxxxxx_0_xxzzz_0, g_xxxxxxx_0_xyyyy_0, g_xxxxxxx_0_xyyyz_0, g_xxxxxxx_0_xyyzz_0, g_xxxxxxx_0_xyzzz_0, g_xxxxxxx_0_xzzzz_0, g_xxxxxxx_0_yyyyy_0, g_xxxxxxx_0_yyyyz_0, g_xxxxxxx_0_yyyzz_0, g_xxxxxxx_0_yyzzz_0, g_xxxxxxx_0_yzzzz_0, g_xxxxxxx_0_zzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxxx_0_xxxxx_0[i] = 6.0 * g_xxxxx_0_xxxxx_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxx_1[i] * fz_be_0 + 5.0 * g_xxxxxx_0_xxxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxx_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxy_0[i] = 6.0 * g_xxxxx_0_xxxxy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxz_0[i] = 6.0 * g_xxxxx_0_xxxxz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyy_0[i] = 6.0 * g_xxxxx_0_xxxyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyz_0[i] = 6.0 * g_xxxxx_0_xxxyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxzz_0[i] = 6.0 * g_xxxxx_0_xxxzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyy_0[i] = 6.0 * g_xxxxx_0_xxyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyz_0[i] = 6.0 * g_xxxxx_0_xxyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyzz_0[i] = 6.0 * g_xxxxx_0_xxyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxzzz_0[i] = 6.0 * g_xxxxx_0_xxzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyy_0[i] = 6.0 * g_xxxxx_0_xyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyy_1[i] * fz_be_0 + g_xxxxxx_0_yyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyz_0[i] = 6.0 * g_xxxxx_0_xyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyz_1[i] * fz_be_0 + g_xxxxxx_0_yyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyzz_0[i] = 6.0 * g_xxxxx_0_xyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyzz_1[i] * fz_be_0 + g_xxxxxx_0_yyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyzzz_0[i] = 6.0 * g_xxxxx_0_xyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyzzz_1[i] * fz_be_0 + g_xxxxxx_0_yzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xzzzz_0[i] = 6.0 * g_xxxxx_0_xzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xzzzz_1[i] * fz_be_0 + g_xxxxxx_0_zzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyy_0[i] = 6.0 * g_xxxxx_0_yyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyy_1[i] * fz_be_0 + g_xxxxxx_0_yyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyz_0[i] = 6.0 * g_xxxxx_0_yyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyzz_0[i] = 6.0 * g_xxxxx_0_yyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyzzz_0[i] = 6.0 * g_xxxxx_0_yyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yzzzz_0[i] = 6.0 * g_xxxxx_0_yzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_zzzzz_0[i] = 6.0 * g_xxxxx_0_zzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_zzzzz_1[i] * fz_be_0 + g_xxxxxx_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 21-42 components of targeted buffer : KSH

    auto g_xxxxxxy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 21);

    auto g_xxxxxxy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 22);

    auto g_xxxxxxy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 23);

    auto g_xxxxxxy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 24);

    auto g_xxxxxxy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 25);

    auto g_xxxxxxy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 26);

    auto g_xxxxxxy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 27);

    auto g_xxxxxxy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 28);

    auto g_xxxxxxy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 29);

    auto g_xxxxxxy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 30);

    auto g_xxxxxxy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 31);

    auto g_xxxxxxy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 32);

    auto g_xxxxxxy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 33);

    auto g_xxxxxxy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 34);

    auto g_xxxxxxy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 35);

    auto g_xxxxxxy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 36);

    auto g_xxxxxxy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 37);

    auto g_xxxxxxy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 38);

    auto g_xxxxxxy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 39);

    auto g_xxxxxxy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 40);

    auto g_xxxxxxy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 41);

    #pragma omp simd aligned(g_xxxxxx_0_xxxx_1, g_xxxxxx_0_xxxxx_1, g_xxxxxx_0_xxxxy_1, g_xxxxxx_0_xxxxz_1, g_xxxxxx_0_xxxy_1, g_xxxxxx_0_xxxyy_1, g_xxxxxx_0_xxxyz_1, g_xxxxxx_0_xxxz_1, g_xxxxxx_0_xxxzz_1, g_xxxxxx_0_xxyy_1, g_xxxxxx_0_xxyyy_1, g_xxxxxx_0_xxyyz_1, g_xxxxxx_0_xxyz_1, g_xxxxxx_0_xxyzz_1, g_xxxxxx_0_xxzz_1, g_xxxxxx_0_xxzzz_1, g_xxxxxx_0_xyyy_1, g_xxxxxx_0_xyyyy_1, g_xxxxxx_0_xyyyz_1, g_xxxxxx_0_xyyz_1, g_xxxxxx_0_xyyzz_1, g_xxxxxx_0_xyzz_1, g_xxxxxx_0_xyzzz_1, g_xxxxxx_0_xzzz_1, g_xxxxxx_0_xzzzz_1, g_xxxxxx_0_yyyy_1, g_xxxxxx_0_yyyyy_1, g_xxxxxx_0_yyyyz_1, g_xxxxxx_0_yyyz_1, g_xxxxxx_0_yyyzz_1, g_xxxxxx_0_yyzz_1, g_xxxxxx_0_yyzzz_1, g_xxxxxx_0_yzzz_1, g_xxxxxx_0_yzzzz_1, g_xxxxxx_0_zzzz_1, g_xxxxxx_0_zzzzz_1, g_xxxxxxy_0_xxxxx_0, g_xxxxxxy_0_xxxxy_0, g_xxxxxxy_0_xxxxz_0, g_xxxxxxy_0_xxxyy_0, g_xxxxxxy_0_xxxyz_0, g_xxxxxxy_0_xxxzz_0, g_xxxxxxy_0_xxyyy_0, g_xxxxxxy_0_xxyyz_0, g_xxxxxxy_0_xxyzz_0, g_xxxxxxy_0_xxzzz_0, g_xxxxxxy_0_xyyyy_0, g_xxxxxxy_0_xyyyz_0, g_xxxxxxy_0_xyyzz_0, g_xxxxxxy_0_xyzzz_0, g_xxxxxxy_0_xzzzz_0, g_xxxxxxy_0_yyyyy_0, g_xxxxxxy_0_yyyyz_0, g_xxxxxxy_0_yyyzz_0, g_xxxxxxy_0_yyzzz_0, g_xxxxxxy_0_yzzzz_0, g_xxxxxxy_0_zzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxy_0_xxxxx_0[i] = g_xxxxxx_0_xxxxx_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxy_0[i] = g_xxxxxx_0_xxxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxz_0[i] = g_xxxxxx_0_xxxxz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyy_0[i] = 2.0 * g_xxxxxx_0_xxxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyz_0[i] = g_xxxxxx_0_xxxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxzz_0[i] = g_xxxxxx_0_xxxzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyy_0[i] = 3.0 * g_xxxxxx_0_xxyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyz_0[i] = 2.0 * g_xxxxxx_0_xxyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyzz_0[i] = g_xxxxxx_0_xxzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxzzz_0[i] = g_xxxxxx_0_xxzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyy_0[i] = 4.0 * g_xxxxxx_0_xyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyz_0[i] = 3.0 * g_xxxxxx_0_xyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyzz_0[i] = 2.0 * g_xxxxxx_0_xyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyzzz_0[i] = g_xxxxxx_0_xzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xzzzz_0[i] = g_xxxxxx_0_xzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyy_0[i] = 5.0 * g_xxxxxx_0_yyyy_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyz_0[i] = 4.0 * g_xxxxxx_0_yyyz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyzz_0[i] = 3.0 * g_xxxxxx_0_yyzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyzzz_0[i] = 2.0 * g_xxxxxx_0_yzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yzzzz_0[i] = g_xxxxxx_0_zzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_zzzzz_0[i] = g_xxxxxx_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 42-63 components of targeted buffer : KSH

    auto g_xxxxxxz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 42);

    auto g_xxxxxxz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 43);

    auto g_xxxxxxz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 44);

    auto g_xxxxxxz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 45);

    auto g_xxxxxxz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 46);

    auto g_xxxxxxz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 47);

    auto g_xxxxxxz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 48);

    auto g_xxxxxxz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 49);

    auto g_xxxxxxz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 50);

    auto g_xxxxxxz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 51);

    auto g_xxxxxxz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 52);

    auto g_xxxxxxz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 53);

    auto g_xxxxxxz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 54);

    auto g_xxxxxxz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 55);

    auto g_xxxxxxz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 56);

    auto g_xxxxxxz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 57);

    auto g_xxxxxxz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 58);

    auto g_xxxxxxz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 59);

    auto g_xxxxxxz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 60);

    auto g_xxxxxxz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 61);

    auto g_xxxxxxz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 62);

    #pragma omp simd aligned(g_xxxxxx_0_xxxx_1, g_xxxxxx_0_xxxxx_1, g_xxxxxx_0_xxxxy_1, g_xxxxxx_0_xxxxz_1, g_xxxxxx_0_xxxy_1, g_xxxxxx_0_xxxyy_1, g_xxxxxx_0_xxxyz_1, g_xxxxxx_0_xxxz_1, g_xxxxxx_0_xxxzz_1, g_xxxxxx_0_xxyy_1, g_xxxxxx_0_xxyyy_1, g_xxxxxx_0_xxyyz_1, g_xxxxxx_0_xxyz_1, g_xxxxxx_0_xxyzz_1, g_xxxxxx_0_xxzz_1, g_xxxxxx_0_xxzzz_1, g_xxxxxx_0_xyyy_1, g_xxxxxx_0_xyyyy_1, g_xxxxxx_0_xyyyz_1, g_xxxxxx_0_xyyz_1, g_xxxxxx_0_xyyzz_1, g_xxxxxx_0_xyzz_1, g_xxxxxx_0_xyzzz_1, g_xxxxxx_0_xzzz_1, g_xxxxxx_0_xzzzz_1, g_xxxxxx_0_yyyy_1, g_xxxxxx_0_yyyyy_1, g_xxxxxx_0_yyyyz_1, g_xxxxxx_0_yyyz_1, g_xxxxxx_0_yyyzz_1, g_xxxxxx_0_yyzz_1, g_xxxxxx_0_yyzzz_1, g_xxxxxx_0_yzzz_1, g_xxxxxx_0_yzzzz_1, g_xxxxxx_0_zzzz_1, g_xxxxxx_0_zzzzz_1, g_xxxxxxz_0_xxxxx_0, g_xxxxxxz_0_xxxxy_0, g_xxxxxxz_0_xxxxz_0, g_xxxxxxz_0_xxxyy_0, g_xxxxxxz_0_xxxyz_0, g_xxxxxxz_0_xxxzz_0, g_xxxxxxz_0_xxyyy_0, g_xxxxxxz_0_xxyyz_0, g_xxxxxxz_0_xxyzz_0, g_xxxxxxz_0_xxzzz_0, g_xxxxxxz_0_xyyyy_0, g_xxxxxxz_0_xyyyz_0, g_xxxxxxz_0_xyyzz_0, g_xxxxxxz_0_xyzzz_0, g_xxxxxxz_0_xzzzz_0, g_xxxxxxz_0_yyyyy_0, g_xxxxxxz_0_yyyyz_0, g_xxxxxxz_0_yyyzz_0, g_xxxxxxz_0_yyzzz_0, g_xxxxxxz_0_yzzzz_0, g_xxxxxxz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxz_0_xxxxx_0[i] = g_xxxxxx_0_xxxxx_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxy_0[i] = g_xxxxxx_0_xxxxy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxz_0[i] = g_xxxxxx_0_xxxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyy_0[i] = g_xxxxxx_0_xxxyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyz_0[i] = g_xxxxxx_0_xxxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxzz_0[i] = 2.0 * g_xxxxxx_0_xxxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyy_0[i] = g_xxxxxx_0_xxyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyz_0[i] = g_xxxxxx_0_xxyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyzz_0[i] = 2.0 * g_xxxxxx_0_xxyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxzzz_0[i] = 3.0 * g_xxxxxx_0_xxzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyy_0[i] = g_xxxxxx_0_xyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyz_0[i] = g_xxxxxx_0_xyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyzz_0[i] = 2.0 * g_xxxxxx_0_xyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyzzz_0[i] = 3.0 * g_xxxxxx_0_xyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xzzzz_0[i] = 4.0 * g_xxxxxx_0_xzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyy_0[i] = g_xxxxxx_0_yyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyz_0[i] = g_xxxxxx_0_yyyy_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyzz_0[i] = 2.0 * g_xxxxxx_0_yyyz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyzzz_0[i] = 3.0 * g_xxxxxx_0_yyzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yzzzz_0[i] = 4.0 * g_xxxxxx_0_yzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_zzzzz_0[i] = 5.0 * g_xxxxxx_0_zzzz_1[i] * fi_acd_0 + g_xxxxxx_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 63-84 components of targeted buffer : KSH

    auto g_xxxxxyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 63);

    auto g_xxxxxyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 64);

    auto g_xxxxxyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 65);

    auto g_xxxxxyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 66);

    auto g_xxxxxyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 67);

    auto g_xxxxxyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 68);

    auto g_xxxxxyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 69);

    auto g_xxxxxyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 70);

    auto g_xxxxxyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 71);

    auto g_xxxxxyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 72);

    auto g_xxxxxyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 73);

    auto g_xxxxxyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 74);

    auto g_xxxxxyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 75);

    auto g_xxxxxyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 76);

    auto g_xxxxxyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 77);

    auto g_xxxxxyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 78);

    auto g_xxxxxyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 79);

    auto g_xxxxxyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 80);

    auto g_xxxxxyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 81);

    auto g_xxxxxyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 82);

    auto g_xxxxxyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 83);

    #pragma omp simd aligned(g_xxxxx_0_xxxxx_0, g_xxxxx_0_xxxxx_1, g_xxxxx_0_xxxxz_0, g_xxxxx_0_xxxxz_1, g_xxxxx_0_xxxzz_0, g_xxxxx_0_xxxzz_1, g_xxxxx_0_xxzzz_0, g_xxxxx_0_xxzzz_1, g_xxxxx_0_xzzzz_0, g_xxxxx_0_xzzzz_1, g_xxxxxy_0_xxxxx_1, g_xxxxxy_0_xxxxz_1, g_xxxxxy_0_xxxzz_1, g_xxxxxy_0_xxzzz_1, g_xxxxxy_0_xzzzz_1, g_xxxxxyy_0_xxxxx_0, g_xxxxxyy_0_xxxxy_0, g_xxxxxyy_0_xxxxz_0, g_xxxxxyy_0_xxxyy_0, g_xxxxxyy_0_xxxyz_0, g_xxxxxyy_0_xxxzz_0, g_xxxxxyy_0_xxyyy_0, g_xxxxxyy_0_xxyyz_0, g_xxxxxyy_0_xxyzz_0, g_xxxxxyy_0_xxzzz_0, g_xxxxxyy_0_xyyyy_0, g_xxxxxyy_0_xyyyz_0, g_xxxxxyy_0_xyyzz_0, g_xxxxxyy_0_xyzzz_0, g_xxxxxyy_0_xzzzz_0, g_xxxxxyy_0_yyyyy_0, g_xxxxxyy_0_yyyyz_0, g_xxxxxyy_0_yyyzz_0, g_xxxxxyy_0_yyzzz_0, g_xxxxxyy_0_yzzzz_0, g_xxxxxyy_0_zzzzz_0, g_xxxxyy_0_xxxxy_1, g_xxxxyy_0_xxxy_1, g_xxxxyy_0_xxxyy_1, g_xxxxyy_0_xxxyz_1, g_xxxxyy_0_xxyy_1, g_xxxxyy_0_xxyyy_1, g_xxxxyy_0_xxyyz_1, g_xxxxyy_0_xxyz_1, g_xxxxyy_0_xxyzz_1, g_xxxxyy_0_xyyy_1, g_xxxxyy_0_xyyyy_1, g_xxxxyy_0_xyyyz_1, g_xxxxyy_0_xyyz_1, g_xxxxyy_0_xyyzz_1, g_xxxxyy_0_xyzz_1, g_xxxxyy_0_xyzzz_1, g_xxxxyy_0_yyyy_1, g_xxxxyy_0_yyyyy_1, g_xxxxyy_0_yyyyz_1, g_xxxxyy_0_yyyz_1, g_xxxxyy_0_yyyzz_1, g_xxxxyy_0_yyzz_1, g_xxxxyy_0_yyzzz_1, g_xxxxyy_0_yzzz_1, g_xxxxyy_0_yzzzz_1, g_xxxxyy_0_zzzzz_1, g_xxxyy_0_xxxxy_0, g_xxxyy_0_xxxxy_1, g_xxxyy_0_xxxyy_0, g_xxxyy_0_xxxyy_1, g_xxxyy_0_xxxyz_0, g_xxxyy_0_xxxyz_1, g_xxxyy_0_xxyyy_0, g_xxxyy_0_xxyyy_1, g_xxxyy_0_xxyyz_0, g_xxxyy_0_xxyyz_1, g_xxxyy_0_xxyzz_0, g_xxxyy_0_xxyzz_1, g_xxxyy_0_xyyyy_0, g_xxxyy_0_xyyyy_1, g_xxxyy_0_xyyyz_0, g_xxxyy_0_xyyyz_1, g_xxxyy_0_xyyzz_0, g_xxxyy_0_xyyzz_1, g_xxxyy_0_xyzzz_0, g_xxxyy_0_xyzzz_1, g_xxxyy_0_yyyyy_0, g_xxxyy_0_yyyyy_1, g_xxxyy_0_yyyyz_0, g_xxxyy_0_yyyyz_1, g_xxxyy_0_yyyzz_0, g_xxxyy_0_yyyzz_1, g_xxxyy_0_yyzzz_0, g_xxxyy_0_yyzzz_1, g_xxxyy_0_yzzzz_0, g_xxxyy_0_yzzzz_1, g_xxxyy_0_zzzzz_0, g_xxxyy_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxyy_0_xxxxx_0[i] = g_xxxxx_0_xxxxx_0[i] * fbe_0 - g_xxxxx_0_xxxxx_1[i] * fz_be_0 + g_xxxxxy_0_xxxxx_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxxy_0[i] = 4.0 * g_xxxyy_0_xxxxy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxxxyy_0_xxxy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxz_0[i] = g_xxxxx_0_xxxxz_0[i] * fbe_0 - g_xxxxx_0_xxxxz_1[i] * fz_be_0 + g_xxxxxy_0_xxxxz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxyy_0[i] = 4.0 * g_xxxyy_0_xxxyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxyz_0[i] = 4.0 * g_xxxyy_0_xxxyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxzz_0[i] = g_xxxxx_0_xxxzz_0[i] * fbe_0 - g_xxxxx_0_xxxzz_1[i] * fz_be_0 + g_xxxxxy_0_xxxzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxyyy_0[i] = 4.0 * g_xxxyy_0_xxyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyyz_0[i] = 4.0 * g_xxxyy_0_xxyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyzz_0[i] = 4.0 * g_xxxyy_0_xxyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxzzz_0[i] = g_xxxxx_0_xxzzz_0[i] * fbe_0 - g_xxxxx_0_xxzzz_1[i] * fz_be_0 + g_xxxxxy_0_xxzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xyyyy_0[i] = 4.0 * g_xxxyy_0_xyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyy_1[i] * fz_be_0 + g_xxxxyy_0_yyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyyz_0[i] = 4.0 * g_xxxyy_0_xyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyz_1[i] * fz_be_0 + g_xxxxyy_0_yyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyzz_0[i] = 4.0 * g_xxxyy_0_xyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyzz_1[i] * fz_be_0 + g_xxxxyy_0_yyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyzzz_0[i] = 4.0 * g_xxxyy_0_xyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyzzz_1[i] * fz_be_0 + g_xxxxyy_0_yzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xzzzz_0[i] = g_xxxxx_0_xzzzz_0[i] * fbe_0 - g_xxxxx_0_xzzzz_1[i] * fz_be_0 + g_xxxxxy_0_xzzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_yyyyy_0[i] = 4.0 * g_xxxyy_0_yyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyy_1[i] * fz_be_0 + g_xxxxyy_0_yyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyyz_0[i] = 4.0 * g_xxxyy_0_yyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyzz_0[i] = 4.0 * g_xxxyy_0_yyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyzzz_0[i] = 4.0 * g_xxxyy_0_yyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yzzzz_0[i] = 4.0 * g_xxxyy_0_yzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_zzzzz_0[i] = 4.0 * g_xxxyy_0_zzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_zzzzz_1[i] * fz_be_0 + g_xxxxyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 84-105 components of targeted buffer : KSH

    auto g_xxxxxyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 84);

    auto g_xxxxxyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 85);

    auto g_xxxxxyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 86);

    auto g_xxxxxyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 87);

    auto g_xxxxxyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 88);

    auto g_xxxxxyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 89);

    auto g_xxxxxyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 90);

    auto g_xxxxxyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 91);

    auto g_xxxxxyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 92);

    auto g_xxxxxyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 93);

    auto g_xxxxxyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 94);

    auto g_xxxxxyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 95);

    auto g_xxxxxyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 96);

    auto g_xxxxxyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 97);

    auto g_xxxxxyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 98);

    auto g_xxxxxyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 99);

    auto g_xxxxxyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 100);

    auto g_xxxxxyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 101);

    auto g_xxxxxyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 102);

    auto g_xxxxxyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 103);

    auto g_xxxxxyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 104);

    #pragma omp simd aligned(g_xxxxxy_0_xxxxy_1, g_xxxxxy_0_xxxyy_1, g_xxxxxy_0_xxyyy_1, g_xxxxxy_0_xyyyy_1, g_xxxxxy_0_yyyyy_1, g_xxxxxyz_0_xxxxx_0, g_xxxxxyz_0_xxxxy_0, g_xxxxxyz_0_xxxxz_0, g_xxxxxyz_0_xxxyy_0, g_xxxxxyz_0_xxxyz_0, g_xxxxxyz_0_xxxzz_0, g_xxxxxyz_0_xxyyy_0, g_xxxxxyz_0_xxyyz_0, g_xxxxxyz_0_xxyzz_0, g_xxxxxyz_0_xxzzz_0, g_xxxxxyz_0_xyyyy_0, g_xxxxxyz_0_xyyyz_0, g_xxxxxyz_0_xyyzz_0, g_xxxxxyz_0_xyzzz_0, g_xxxxxyz_0_xzzzz_0, g_xxxxxyz_0_yyyyy_0, g_xxxxxyz_0_yyyyz_0, g_xxxxxyz_0_yyyzz_0, g_xxxxxyz_0_yyzzz_0, g_xxxxxyz_0_yzzzz_0, g_xxxxxyz_0_zzzzz_0, g_xxxxxz_0_xxxxx_1, g_xxxxxz_0_xxxxz_1, g_xxxxxz_0_xxxyz_1, g_xxxxxz_0_xxxz_1, g_xxxxxz_0_xxxzz_1, g_xxxxxz_0_xxyyz_1, g_xxxxxz_0_xxyz_1, g_xxxxxz_0_xxyzz_1, g_xxxxxz_0_xxzz_1, g_xxxxxz_0_xxzzz_1, g_xxxxxz_0_xyyyz_1, g_xxxxxz_0_xyyz_1, g_xxxxxz_0_xyyzz_1, g_xxxxxz_0_xyzz_1, g_xxxxxz_0_xyzzz_1, g_xxxxxz_0_xzzz_1, g_xxxxxz_0_xzzzz_1, g_xxxxxz_0_yyyyz_1, g_xxxxxz_0_yyyz_1, g_xxxxxz_0_yyyzz_1, g_xxxxxz_0_yyzz_1, g_xxxxxz_0_yyzzz_1, g_xxxxxz_0_yzzz_1, g_xxxxxz_0_yzzzz_1, g_xxxxxz_0_zzzz_1, g_xxxxxz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxyz_0_xxxxx_0[i] = g_xxxxxz_0_xxxxx_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxy_0[i] = g_xxxxxy_0_xxxxy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxxz_0[i] = g_xxxxxz_0_xxxxz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxyy_0[i] = g_xxxxxy_0_xxxyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxyz_0[i] = g_xxxxxz_0_xxxz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxzz_0[i] = g_xxxxxz_0_xxxzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyyy_0[i] = g_xxxxxy_0_xxyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxyyz_0[i] = 2.0 * g_xxxxxz_0_xxyz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyzz_0[i] = g_xxxxxz_0_xxzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxzzz_0[i] = g_xxxxxz_0_xxzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyyy_0[i] = g_xxxxxy_0_xyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xyyyz_0[i] = 3.0 * g_xxxxxz_0_xyyz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyzz_0[i] = 2.0 * g_xxxxxz_0_xyzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyzzz_0[i] = g_xxxxxz_0_xzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xzzzz_0[i] = g_xxxxxz_0_xzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyyy_0[i] = g_xxxxxy_0_yyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_yyyyz_0[i] = 4.0 * g_xxxxxz_0_yyyz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyzz_0[i] = 3.0 * g_xxxxxz_0_yyzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyzzz_0[i] = 2.0 * g_xxxxxz_0_yzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yzzzz_0[i] = g_xxxxxz_0_zzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_zzzzz_0[i] = g_xxxxxz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 105-126 components of targeted buffer : KSH

    auto g_xxxxxzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 105);

    auto g_xxxxxzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 106);

    auto g_xxxxxzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 107);

    auto g_xxxxxzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 108);

    auto g_xxxxxzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 109);

    auto g_xxxxxzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 110);

    auto g_xxxxxzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 111);

    auto g_xxxxxzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 112);

    auto g_xxxxxzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 113);

    auto g_xxxxxzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 114);

    auto g_xxxxxzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 115);

    auto g_xxxxxzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 116);

    auto g_xxxxxzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 117);

    auto g_xxxxxzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 118);

    auto g_xxxxxzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 119);

    auto g_xxxxxzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 120);

    auto g_xxxxxzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 121);

    auto g_xxxxxzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 122);

    auto g_xxxxxzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 123);

    auto g_xxxxxzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 124);

    auto g_xxxxxzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 125);

    #pragma omp simd aligned(g_xxxxx_0_xxxxx_0, g_xxxxx_0_xxxxx_1, g_xxxxx_0_xxxxy_0, g_xxxxx_0_xxxxy_1, g_xxxxx_0_xxxyy_0, g_xxxxx_0_xxxyy_1, g_xxxxx_0_xxyyy_0, g_xxxxx_0_xxyyy_1, g_xxxxx_0_xyyyy_0, g_xxxxx_0_xyyyy_1, g_xxxxxz_0_xxxxx_1, g_xxxxxz_0_xxxxy_1, g_xxxxxz_0_xxxyy_1, g_xxxxxz_0_xxyyy_1, g_xxxxxz_0_xyyyy_1, g_xxxxxzz_0_xxxxx_0, g_xxxxxzz_0_xxxxy_0, g_xxxxxzz_0_xxxxz_0, g_xxxxxzz_0_xxxyy_0, g_xxxxxzz_0_xxxyz_0, g_xxxxxzz_0_xxxzz_0, g_xxxxxzz_0_xxyyy_0, g_xxxxxzz_0_xxyyz_0, g_xxxxxzz_0_xxyzz_0, g_xxxxxzz_0_xxzzz_0, g_xxxxxzz_0_xyyyy_0, g_xxxxxzz_0_xyyyz_0, g_xxxxxzz_0_xyyzz_0, g_xxxxxzz_0_xyzzz_0, g_xxxxxzz_0_xzzzz_0, g_xxxxxzz_0_yyyyy_0, g_xxxxxzz_0_yyyyz_0, g_xxxxxzz_0_yyyzz_0, g_xxxxxzz_0_yyzzz_0, g_xxxxxzz_0_yzzzz_0, g_xxxxxzz_0_zzzzz_0, g_xxxxzz_0_xxxxz_1, g_xxxxzz_0_xxxyz_1, g_xxxxzz_0_xxxz_1, g_xxxxzz_0_xxxzz_1, g_xxxxzz_0_xxyyz_1, g_xxxxzz_0_xxyz_1, g_xxxxzz_0_xxyzz_1, g_xxxxzz_0_xxzz_1, g_xxxxzz_0_xxzzz_1, g_xxxxzz_0_xyyyz_1, g_xxxxzz_0_xyyz_1, g_xxxxzz_0_xyyzz_1, g_xxxxzz_0_xyzz_1, g_xxxxzz_0_xyzzz_1, g_xxxxzz_0_xzzz_1, g_xxxxzz_0_xzzzz_1, g_xxxxzz_0_yyyyy_1, g_xxxxzz_0_yyyyz_1, g_xxxxzz_0_yyyz_1, g_xxxxzz_0_yyyzz_1, g_xxxxzz_0_yyzz_1, g_xxxxzz_0_yyzzz_1, g_xxxxzz_0_yzzz_1, g_xxxxzz_0_yzzzz_1, g_xxxxzz_0_zzzz_1, g_xxxxzz_0_zzzzz_1, g_xxxzz_0_xxxxz_0, g_xxxzz_0_xxxxz_1, g_xxxzz_0_xxxyz_0, g_xxxzz_0_xxxyz_1, g_xxxzz_0_xxxzz_0, g_xxxzz_0_xxxzz_1, g_xxxzz_0_xxyyz_0, g_xxxzz_0_xxyyz_1, g_xxxzz_0_xxyzz_0, g_xxxzz_0_xxyzz_1, g_xxxzz_0_xxzzz_0, g_xxxzz_0_xxzzz_1, g_xxxzz_0_xyyyz_0, g_xxxzz_0_xyyyz_1, g_xxxzz_0_xyyzz_0, g_xxxzz_0_xyyzz_1, g_xxxzz_0_xyzzz_0, g_xxxzz_0_xyzzz_1, g_xxxzz_0_xzzzz_0, g_xxxzz_0_xzzzz_1, g_xxxzz_0_yyyyy_0, g_xxxzz_0_yyyyy_1, g_xxxzz_0_yyyyz_0, g_xxxzz_0_yyyyz_1, g_xxxzz_0_yyyzz_0, g_xxxzz_0_yyyzz_1, g_xxxzz_0_yyzzz_0, g_xxxzz_0_yyzzz_1, g_xxxzz_0_yzzzz_0, g_xxxzz_0_yzzzz_1, g_xxxzz_0_zzzzz_0, g_xxxzz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxzz_0_xxxxx_0[i] = g_xxxxx_0_xxxxx_0[i] * fbe_0 - g_xxxxx_0_xxxxx_1[i] * fz_be_0 + g_xxxxxz_0_xxxxx_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxy_0[i] = g_xxxxx_0_xxxxy_0[i] * fbe_0 - g_xxxxx_0_xxxxy_1[i] * fz_be_0 + g_xxxxxz_0_xxxxy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxz_0[i] = 4.0 * g_xxxzz_0_xxxxz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_0_xxxz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxyy_0[i] = g_xxxxx_0_xxxyy_0[i] * fbe_0 - g_xxxxx_0_xxxyy_1[i] * fz_be_0 + g_xxxxxz_0_xxxyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxyz_0[i] = 4.0 * g_xxxzz_0_xxxyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxzz_0[i] = 4.0 * g_xxxzz_0_xxxzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyyy_0[i] = g_xxxxx_0_xxyyy_0[i] * fbe_0 - g_xxxxx_0_xxyyy_1[i] * fz_be_0 + g_xxxxxz_0_xxyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxyyz_0[i] = 4.0 * g_xxxzz_0_xxyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyzz_0[i] = 4.0 * g_xxxzz_0_xxyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxzzz_0[i] = 4.0 * g_xxxzz_0_xxzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyyy_0[i] = g_xxxxx_0_xyyyy_0[i] * fbe_0 - g_xxxxx_0_xyyyy_1[i] * fz_be_0 + g_xxxxxz_0_xyyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xyyyz_0[i] = 4.0 * g_xxxzz_0_xyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyyz_1[i] * fz_be_0 + g_xxxxzz_0_yyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyzz_0[i] = 4.0 * g_xxxzz_0_xyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyzz_1[i] * fz_be_0 + g_xxxxzz_0_yyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyzzz_0[i] = 4.0 * g_xxxzz_0_xyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyzzz_1[i] * fz_be_0 + g_xxxxzz_0_yzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xzzzz_0[i] = 4.0 * g_xxxzz_0_xzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xzzzz_1[i] * fz_be_0 + g_xxxxzz_0_zzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyy_0[i] = 4.0 * g_xxxzz_0_yyyyy_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyy_1[i] * fz_be_0 + g_xxxxzz_0_yyyyy_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyz_0[i] = 4.0 * g_xxxzz_0_yyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyzz_0[i] = 4.0 * g_xxxzz_0_yyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyzzz_0[i] = 4.0 * g_xxxzz_0_yyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yzzzz_0[i] = 4.0 * g_xxxzz_0_yzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_zzzzz_0[i] = 4.0 * g_xxxzz_0_zzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_zzzzz_1[i] * fz_be_0 + g_xxxxzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 126-147 components of targeted buffer : KSH

    auto g_xxxxyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 126);

    auto g_xxxxyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 127);

    auto g_xxxxyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 128);

    auto g_xxxxyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 129);

    auto g_xxxxyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 130);

    auto g_xxxxyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 131);

    auto g_xxxxyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 132);

    auto g_xxxxyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 133);

    auto g_xxxxyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 134);

    auto g_xxxxyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 135);

    auto g_xxxxyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 136);

    auto g_xxxxyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 137);

    auto g_xxxxyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 138);

    auto g_xxxxyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 139);

    auto g_xxxxyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 140);

    auto g_xxxxyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 141);

    auto g_xxxxyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 142);

    auto g_xxxxyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 143);

    auto g_xxxxyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 144);

    auto g_xxxxyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 145);

    auto g_xxxxyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 146);

    #pragma omp simd aligned(g_xxxxy_0_xxxxx_0, g_xxxxy_0_xxxxx_1, g_xxxxy_0_xxxxz_0, g_xxxxy_0_xxxxz_1, g_xxxxy_0_xxxzz_0, g_xxxxy_0_xxxzz_1, g_xxxxy_0_xxzzz_0, g_xxxxy_0_xxzzz_1, g_xxxxy_0_xzzzz_0, g_xxxxy_0_xzzzz_1, g_xxxxyy_0_xxxxx_1, g_xxxxyy_0_xxxxz_1, g_xxxxyy_0_xxxzz_1, g_xxxxyy_0_xxzzz_1, g_xxxxyy_0_xzzzz_1, g_xxxxyyy_0_xxxxx_0, g_xxxxyyy_0_xxxxy_0, g_xxxxyyy_0_xxxxz_0, g_xxxxyyy_0_xxxyy_0, g_xxxxyyy_0_xxxyz_0, g_xxxxyyy_0_xxxzz_0, g_xxxxyyy_0_xxyyy_0, g_xxxxyyy_0_xxyyz_0, g_xxxxyyy_0_xxyzz_0, g_xxxxyyy_0_xxzzz_0, g_xxxxyyy_0_xyyyy_0, g_xxxxyyy_0_xyyyz_0, g_xxxxyyy_0_xyyzz_0, g_xxxxyyy_0_xyzzz_0, g_xxxxyyy_0_xzzzz_0, g_xxxxyyy_0_yyyyy_0, g_xxxxyyy_0_yyyyz_0, g_xxxxyyy_0_yyyzz_0, g_xxxxyyy_0_yyzzz_0, g_xxxxyyy_0_yzzzz_0, g_xxxxyyy_0_zzzzz_0, g_xxxyyy_0_xxxxy_1, g_xxxyyy_0_xxxy_1, g_xxxyyy_0_xxxyy_1, g_xxxyyy_0_xxxyz_1, g_xxxyyy_0_xxyy_1, g_xxxyyy_0_xxyyy_1, g_xxxyyy_0_xxyyz_1, g_xxxyyy_0_xxyz_1, g_xxxyyy_0_xxyzz_1, g_xxxyyy_0_xyyy_1, g_xxxyyy_0_xyyyy_1, g_xxxyyy_0_xyyyz_1, g_xxxyyy_0_xyyz_1, g_xxxyyy_0_xyyzz_1, g_xxxyyy_0_xyzz_1, g_xxxyyy_0_xyzzz_1, g_xxxyyy_0_yyyy_1, g_xxxyyy_0_yyyyy_1, g_xxxyyy_0_yyyyz_1, g_xxxyyy_0_yyyz_1, g_xxxyyy_0_yyyzz_1, g_xxxyyy_0_yyzz_1, g_xxxyyy_0_yyzzz_1, g_xxxyyy_0_yzzz_1, g_xxxyyy_0_yzzzz_1, g_xxxyyy_0_zzzzz_1, g_xxyyy_0_xxxxy_0, g_xxyyy_0_xxxxy_1, g_xxyyy_0_xxxyy_0, g_xxyyy_0_xxxyy_1, g_xxyyy_0_xxxyz_0, g_xxyyy_0_xxxyz_1, g_xxyyy_0_xxyyy_0, g_xxyyy_0_xxyyy_1, g_xxyyy_0_xxyyz_0, g_xxyyy_0_xxyyz_1, g_xxyyy_0_xxyzz_0, g_xxyyy_0_xxyzz_1, g_xxyyy_0_xyyyy_0, g_xxyyy_0_xyyyy_1, g_xxyyy_0_xyyyz_0, g_xxyyy_0_xyyyz_1, g_xxyyy_0_xyyzz_0, g_xxyyy_0_xyyzz_1, g_xxyyy_0_xyzzz_0, g_xxyyy_0_xyzzz_1, g_xxyyy_0_yyyyy_0, g_xxyyy_0_yyyyy_1, g_xxyyy_0_yyyyz_0, g_xxyyy_0_yyyyz_1, g_xxyyy_0_yyyzz_0, g_xxyyy_0_yyyzz_1, g_xxyyy_0_yyzzz_0, g_xxyyy_0_yyzzz_1, g_xxyyy_0_yzzzz_0, g_xxyyy_0_yzzzz_1, g_xxyyy_0_zzzzz_0, g_xxyyy_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyyy_0_xxxxx_0[i] = 2.0 * g_xxxxy_0_xxxxx_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxx_1[i] * fz_be_0 + g_xxxxyy_0_xxxxx_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxxy_0[i] = 3.0 * g_xxyyy_0_xxxxy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxxyyy_0_xxxy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxz_0[i] = 2.0 * g_xxxxy_0_xxxxz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxz_1[i] * fz_be_0 + g_xxxxyy_0_xxxxz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxyy_0[i] = 3.0 * g_xxyyy_0_xxxyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxyz_0[i] = 3.0 * g_xxyyy_0_xxxyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxzz_0[i] = 2.0 * g_xxxxy_0_xxxzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxzz_1[i] * fz_be_0 + g_xxxxyy_0_xxxzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxyyy_0[i] = 3.0 * g_xxyyy_0_xxyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyyz_0[i] = 3.0 * g_xxyyy_0_xxyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyzz_0[i] = 3.0 * g_xxyyy_0_xxyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxzzz_0[i] = 2.0 * g_xxxxy_0_xxzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxzzz_1[i] * fz_be_0 + g_xxxxyy_0_xxzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xyyyy_0[i] = 3.0 * g_xxyyy_0_xyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyy_1[i] * fz_be_0 + g_xxxyyy_0_yyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyyz_0[i] = 3.0 * g_xxyyy_0_xyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyz_1[i] * fz_be_0 + g_xxxyyy_0_yyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyzz_0[i] = 3.0 * g_xxyyy_0_xyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyzz_1[i] * fz_be_0 + g_xxxyyy_0_yyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyzzz_0[i] = 3.0 * g_xxyyy_0_xyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyzzz_1[i] * fz_be_0 + g_xxxyyy_0_yzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xzzzz_0[i] = 2.0 * g_xxxxy_0_xzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xzzzz_1[i] * fz_be_0 + g_xxxxyy_0_xzzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_yyyyy_0[i] = 3.0 * g_xxyyy_0_yyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyy_1[i] * fz_be_0 + g_xxxyyy_0_yyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyyz_0[i] = 3.0 * g_xxyyy_0_yyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyzz_0[i] = 3.0 * g_xxyyy_0_yyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyzzz_0[i] = 3.0 * g_xxyyy_0_yyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yzzzz_0[i] = 3.0 * g_xxyyy_0_yzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_zzzzz_0[i] = 3.0 * g_xxyyy_0_zzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_zzzzz_1[i] * fz_be_0 + g_xxxyyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 147-168 components of targeted buffer : KSH

    auto g_xxxxyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 147);

    auto g_xxxxyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 148);

    auto g_xxxxyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 149);

    auto g_xxxxyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 150);

    auto g_xxxxyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 151);

    auto g_xxxxyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 152);

    auto g_xxxxyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 153);

    auto g_xxxxyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 154);

    auto g_xxxxyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 155);

    auto g_xxxxyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 156);

    auto g_xxxxyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 157);

    auto g_xxxxyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 158);

    auto g_xxxxyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 159);

    auto g_xxxxyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 160);

    auto g_xxxxyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 161);

    auto g_xxxxyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 162);

    auto g_xxxxyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 163);

    auto g_xxxxyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 164);

    auto g_xxxxyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 165);

    auto g_xxxxyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 166);

    auto g_xxxxyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 167);

    #pragma omp simd aligned(g_xxxxyy_0_xxxx_1, g_xxxxyy_0_xxxxx_1, g_xxxxyy_0_xxxxy_1, g_xxxxyy_0_xxxxz_1, g_xxxxyy_0_xxxy_1, g_xxxxyy_0_xxxyy_1, g_xxxxyy_0_xxxyz_1, g_xxxxyy_0_xxxz_1, g_xxxxyy_0_xxxzz_1, g_xxxxyy_0_xxyy_1, g_xxxxyy_0_xxyyy_1, g_xxxxyy_0_xxyyz_1, g_xxxxyy_0_xxyz_1, g_xxxxyy_0_xxyzz_1, g_xxxxyy_0_xxzz_1, g_xxxxyy_0_xxzzz_1, g_xxxxyy_0_xyyy_1, g_xxxxyy_0_xyyyy_1, g_xxxxyy_0_xyyyz_1, g_xxxxyy_0_xyyz_1, g_xxxxyy_0_xyyzz_1, g_xxxxyy_0_xyzz_1, g_xxxxyy_0_xyzzz_1, g_xxxxyy_0_xzzz_1, g_xxxxyy_0_xzzzz_1, g_xxxxyy_0_yyyy_1, g_xxxxyy_0_yyyyy_1, g_xxxxyy_0_yyyyz_1, g_xxxxyy_0_yyyz_1, g_xxxxyy_0_yyyzz_1, g_xxxxyy_0_yyzz_1, g_xxxxyy_0_yyzzz_1, g_xxxxyy_0_yzzz_1, g_xxxxyy_0_yzzzz_1, g_xxxxyy_0_zzzz_1, g_xxxxyy_0_zzzzz_1, g_xxxxyyz_0_xxxxx_0, g_xxxxyyz_0_xxxxy_0, g_xxxxyyz_0_xxxxz_0, g_xxxxyyz_0_xxxyy_0, g_xxxxyyz_0_xxxyz_0, g_xxxxyyz_0_xxxzz_0, g_xxxxyyz_0_xxyyy_0, g_xxxxyyz_0_xxyyz_0, g_xxxxyyz_0_xxyzz_0, g_xxxxyyz_0_xxzzz_0, g_xxxxyyz_0_xyyyy_0, g_xxxxyyz_0_xyyyz_0, g_xxxxyyz_0_xyyzz_0, g_xxxxyyz_0_xyzzz_0, g_xxxxyyz_0_xzzzz_0, g_xxxxyyz_0_yyyyy_0, g_xxxxyyz_0_yyyyz_0, g_xxxxyyz_0_yyyzz_0, g_xxxxyyz_0_yyzzz_0, g_xxxxyyz_0_yzzzz_0, g_xxxxyyz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyyz_0_xxxxx_0[i] = g_xxxxyy_0_xxxxx_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxy_0[i] = g_xxxxyy_0_xxxxy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxz_0[i] = g_xxxxyy_0_xxxx_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyy_0[i] = g_xxxxyy_0_xxxyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyz_0[i] = g_xxxxyy_0_xxxy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxzz_0[i] = 2.0 * g_xxxxyy_0_xxxz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyy_0[i] = g_xxxxyy_0_xxyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyz_0[i] = g_xxxxyy_0_xxyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyzz_0[i] = 2.0 * g_xxxxyy_0_xxyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxzzz_0[i] = 3.0 * g_xxxxyy_0_xxzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyy_0[i] = g_xxxxyy_0_xyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyz_0[i] = g_xxxxyy_0_xyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyzz_0[i] = 2.0 * g_xxxxyy_0_xyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyzzz_0[i] = 3.0 * g_xxxxyy_0_xyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xzzzz_0[i] = 4.0 * g_xxxxyy_0_xzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyy_0[i] = g_xxxxyy_0_yyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyz_0[i] = g_xxxxyy_0_yyyy_1[i] * fi_acd_0 + g_xxxxyy_0_yyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyzz_0[i] = 2.0 * g_xxxxyy_0_yyyz_1[i] * fi_acd_0 + g_xxxxyy_0_yyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyzzz_0[i] = 3.0 * g_xxxxyy_0_yyzz_1[i] * fi_acd_0 + g_xxxxyy_0_yyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yzzzz_0[i] = 4.0 * g_xxxxyy_0_yzzz_1[i] * fi_acd_0 + g_xxxxyy_0_yzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_zzzzz_0[i] = 5.0 * g_xxxxyy_0_zzzz_1[i] * fi_acd_0 + g_xxxxyy_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 168-189 components of targeted buffer : KSH

    auto g_xxxxyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 168);

    auto g_xxxxyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 169);

    auto g_xxxxyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 170);

    auto g_xxxxyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 171);

    auto g_xxxxyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 172);

    auto g_xxxxyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 173);

    auto g_xxxxyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 174);

    auto g_xxxxyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 175);

    auto g_xxxxyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 176);

    auto g_xxxxyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 177);

    auto g_xxxxyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 178);

    auto g_xxxxyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 179);

    auto g_xxxxyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 180);

    auto g_xxxxyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 181);

    auto g_xxxxyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 182);

    auto g_xxxxyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 183);

    auto g_xxxxyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 184);

    auto g_xxxxyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 185);

    auto g_xxxxyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 186);

    auto g_xxxxyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 187);

    auto g_xxxxyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 188);

    #pragma omp simd aligned(g_xxxxyzz_0_xxxxx_0, g_xxxxyzz_0_xxxxy_0, g_xxxxyzz_0_xxxxz_0, g_xxxxyzz_0_xxxyy_0, g_xxxxyzz_0_xxxyz_0, g_xxxxyzz_0_xxxzz_0, g_xxxxyzz_0_xxyyy_0, g_xxxxyzz_0_xxyyz_0, g_xxxxyzz_0_xxyzz_0, g_xxxxyzz_0_xxzzz_0, g_xxxxyzz_0_xyyyy_0, g_xxxxyzz_0_xyyyz_0, g_xxxxyzz_0_xyyzz_0, g_xxxxyzz_0_xyzzz_0, g_xxxxyzz_0_xzzzz_0, g_xxxxyzz_0_yyyyy_0, g_xxxxyzz_0_yyyyz_0, g_xxxxyzz_0_yyyzz_0, g_xxxxyzz_0_yyzzz_0, g_xxxxyzz_0_yzzzz_0, g_xxxxyzz_0_zzzzz_0, g_xxxxzz_0_xxxx_1, g_xxxxzz_0_xxxxx_1, g_xxxxzz_0_xxxxy_1, g_xxxxzz_0_xxxxz_1, g_xxxxzz_0_xxxy_1, g_xxxxzz_0_xxxyy_1, g_xxxxzz_0_xxxyz_1, g_xxxxzz_0_xxxz_1, g_xxxxzz_0_xxxzz_1, g_xxxxzz_0_xxyy_1, g_xxxxzz_0_xxyyy_1, g_xxxxzz_0_xxyyz_1, g_xxxxzz_0_xxyz_1, g_xxxxzz_0_xxyzz_1, g_xxxxzz_0_xxzz_1, g_xxxxzz_0_xxzzz_1, g_xxxxzz_0_xyyy_1, g_xxxxzz_0_xyyyy_1, g_xxxxzz_0_xyyyz_1, g_xxxxzz_0_xyyz_1, g_xxxxzz_0_xyyzz_1, g_xxxxzz_0_xyzz_1, g_xxxxzz_0_xyzzz_1, g_xxxxzz_0_xzzz_1, g_xxxxzz_0_xzzzz_1, g_xxxxzz_0_yyyy_1, g_xxxxzz_0_yyyyy_1, g_xxxxzz_0_yyyyz_1, g_xxxxzz_0_yyyz_1, g_xxxxzz_0_yyyzz_1, g_xxxxzz_0_yyzz_1, g_xxxxzz_0_yyzzz_1, g_xxxxzz_0_yzzz_1, g_xxxxzz_0_yzzzz_1, g_xxxxzz_0_zzzz_1, g_xxxxzz_0_zzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyzz_0_xxxxx_0[i] = g_xxxxzz_0_xxxxx_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxy_0[i] = g_xxxxzz_0_xxxx_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxz_0[i] = g_xxxxzz_0_xxxxz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyy_0[i] = 2.0 * g_xxxxzz_0_xxxy_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyz_0[i] = g_xxxxzz_0_xxxz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxzz_0[i] = g_xxxxzz_0_xxxzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyy_0[i] = 3.0 * g_xxxxzz_0_xxyy_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyz_0[i] = 2.0 * g_xxxxzz_0_xxyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyzz_0[i] = g_xxxxzz_0_xxzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxzzz_0[i] = g_xxxxzz_0_xxzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyy_0[i] = 4.0 * g_xxxxzz_0_xyyy_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyz_0[i] = 3.0 * g_xxxxzz_0_xyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyzz_0[i] = 2.0 * g_xxxxzz_0_xyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyzzz_0[i] = g_xxxxzz_0_xzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xzzzz_0[i] = g_xxxxzz_0_xzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyy_0[i] = 5.0 * g_xxxxzz_0_yyyy_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyz_0[i] = 4.0 * g_xxxxzz_0_yyyz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyzz_0[i] = 3.0 * g_xxxxzz_0_yyzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyzzz_0[i] = 2.0 * g_xxxxzz_0_yzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yzzzz_0[i] = g_xxxxzz_0_zzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_zzzzz_0[i] = g_xxxxzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 189-210 components of targeted buffer : KSH

    auto g_xxxxzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 189);

    auto g_xxxxzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 190);

    auto g_xxxxzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 191);

    auto g_xxxxzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 192);

    auto g_xxxxzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 193);

    auto g_xxxxzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 194);

    auto g_xxxxzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 195);

    auto g_xxxxzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 196);

    auto g_xxxxzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 197);

    auto g_xxxxzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 198);

    auto g_xxxxzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 199);

    auto g_xxxxzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 200);

    auto g_xxxxzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 201);

    auto g_xxxxzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 202);

    auto g_xxxxzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 203);

    auto g_xxxxzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 204);

    auto g_xxxxzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 205);

    auto g_xxxxzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 206);

    auto g_xxxxzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 207);

    auto g_xxxxzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 208);

    auto g_xxxxzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 209);

    #pragma omp simd aligned(g_xxxxz_0_xxxxx_0, g_xxxxz_0_xxxxx_1, g_xxxxz_0_xxxxy_0, g_xxxxz_0_xxxxy_1, g_xxxxz_0_xxxyy_0, g_xxxxz_0_xxxyy_1, g_xxxxz_0_xxyyy_0, g_xxxxz_0_xxyyy_1, g_xxxxz_0_xyyyy_0, g_xxxxz_0_xyyyy_1, g_xxxxzz_0_xxxxx_1, g_xxxxzz_0_xxxxy_1, g_xxxxzz_0_xxxyy_1, g_xxxxzz_0_xxyyy_1, g_xxxxzz_0_xyyyy_1, g_xxxxzzz_0_xxxxx_0, g_xxxxzzz_0_xxxxy_0, g_xxxxzzz_0_xxxxz_0, g_xxxxzzz_0_xxxyy_0, g_xxxxzzz_0_xxxyz_0, g_xxxxzzz_0_xxxzz_0, g_xxxxzzz_0_xxyyy_0, g_xxxxzzz_0_xxyyz_0, g_xxxxzzz_0_xxyzz_0, g_xxxxzzz_0_xxzzz_0, g_xxxxzzz_0_xyyyy_0, g_xxxxzzz_0_xyyyz_0, g_xxxxzzz_0_xyyzz_0, g_xxxxzzz_0_xyzzz_0, g_xxxxzzz_0_xzzzz_0, g_xxxxzzz_0_yyyyy_0, g_xxxxzzz_0_yyyyz_0, g_xxxxzzz_0_yyyzz_0, g_xxxxzzz_0_yyzzz_0, g_xxxxzzz_0_yzzzz_0, g_xxxxzzz_0_zzzzz_0, g_xxxzzz_0_xxxxz_1, g_xxxzzz_0_xxxyz_1, g_xxxzzz_0_xxxz_1, g_xxxzzz_0_xxxzz_1, g_xxxzzz_0_xxyyz_1, g_xxxzzz_0_xxyz_1, g_xxxzzz_0_xxyzz_1, g_xxxzzz_0_xxzz_1, g_xxxzzz_0_xxzzz_1, g_xxxzzz_0_xyyyz_1, g_xxxzzz_0_xyyz_1, g_xxxzzz_0_xyyzz_1, g_xxxzzz_0_xyzz_1, g_xxxzzz_0_xyzzz_1, g_xxxzzz_0_xzzz_1, g_xxxzzz_0_xzzzz_1, g_xxxzzz_0_yyyyy_1, g_xxxzzz_0_yyyyz_1, g_xxxzzz_0_yyyz_1, g_xxxzzz_0_yyyzz_1, g_xxxzzz_0_yyzz_1, g_xxxzzz_0_yyzzz_1, g_xxxzzz_0_yzzz_1, g_xxxzzz_0_yzzzz_1, g_xxxzzz_0_zzzz_1, g_xxxzzz_0_zzzzz_1, g_xxzzz_0_xxxxz_0, g_xxzzz_0_xxxxz_1, g_xxzzz_0_xxxyz_0, g_xxzzz_0_xxxyz_1, g_xxzzz_0_xxxzz_0, g_xxzzz_0_xxxzz_1, g_xxzzz_0_xxyyz_0, g_xxzzz_0_xxyyz_1, g_xxzzz_0_xxyzz_0, g_xxzzz_0_xxyzz_1, g_xxzzz_0_xxzzz_0, g_xxzzz_0_xxzzz_1, g_xxzzz_0_xyyyz_0, g_xxzzz_0_xyyyz_1, g_xxzzz_0_xyyzz_0, g_xxzzz_0_xyyzz_1, g_xxzzz_0_xyzzz_0, g_xxzzz_0_xyzzz_1, g_xxzzz_0_xzzzz_0, g_xxzzz_0_xzzzz_1, g_xxzzz_0_yyyyy_0, g_xxzzz_0_yyyyy_1, g_xxzzz_0_yyyyz_0, g_xxzzz_0_yyyyz_1, g_xxzzz_0_yyyzz_0, g_xxzzz_0_yyyzz_1, g_xxzzz_0_yyzzz_0, g_xxzzz_0_yyzzz_1, g_xxzzz_0_yzzzz_0, g_xxzzz_0_yzzzz_1, g_xxzzz_0_zzzzz_0, g_xxzzz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzzz_0_xxxxx_0[i] = 2.0 * g_xxxxz_0_xxxxx_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxx_1[i] * fz_be_0 + g_xxxxzz_0_xxxxx_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxy_0[i] = 2.0 * g_xxxxz_0_xxxxy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxy_1[i] * fz_be_0 + g_xxxxzz_0_xxxxy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxz_0[i] = 3.0 * g_xxzzz_0_xxxxz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_0_xxxz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxyy_0[i] = 2.0 * g_xxxxz_0_xxxyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxyy_1[i] * fz_be_0 + g_xxxxzz_0_xxxyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxyz_0[i] = 3.0 * g_xxzzz_0_xxxyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxzz_0[i] = 3.0 * g_xxzzz_0_xxxzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyyy_0[i] = 2.0 * g_xxxxz_0_xxyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxyyy_1[i] * fz_be_0 + g_xxxxzz_0_xxyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxyyz_0[i] = 3.0 * g_xxzzz_0_xxyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyzz_0[i] = 3.0 * g_xxzzz_0_xxyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxzzz_0[i] = 3.0 * g_xxzzz_0_xxzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyyy_0[i] = 2.0 * g_xxxxz_0_xyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xyyyy_1[i] * fz_be_0 + g_xxxxzz_0_xyyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xyyyz_0[i] = 3.0 * g_xxzzz_0_xyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyyz_1[i] * fz_be_0 + g_xxxzzz_0_yyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyzz_0[i] = 3.0 * g_xxzzz_0_xyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyzz_1[i] * fz_be_0 + g_xxxzzz_0_yyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyzzz_0[i] = 3.0 * g_xxzzz_0_xyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyzzz_1[i] * fz_be_0 + g_xxxzzz_0_yzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xzzzz_0[i] = 3.0 * g_xxzzz_0_xzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xzzzz_1[i] * fz_be_0 + g_xxxzzz_0_zzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyy_0[i] = 3.0 * g_xxzzz_0_yyyyy_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyy_1[i] * fz_be_0 + g_xxxzzz_0_yyyyy_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyz_0[i] = 3.0 * g_xxzzz_0_yyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyzz_0[i] = 3.0 * g_xxzzz_0_yyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyzzz_0[i] = 3.0 * g_xxzzz_0_yyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yzzzz_0[i] = 3.0 * g_xxzzz_0_yzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_zzzzz_0[i] = 3.0 * g_xxzzz_0_zzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_zzzzz_1[i] * fz_be_0 + g_xxxzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 210-231 components of targeted buffer : KSH

    auto g_xxxyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 210);

    auto g_xxxyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 211);

    auto g_xxxyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 212);

    auto g_xxxyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 213);

    auto g_xxxyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 214);

    auto g_xxxyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 215);

    auto g_xxxyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 216);

    auto g_xxxyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 217);

    auto g_xxxyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 218);

    auto g_xxxyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 219);

    auto g_xxxyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 220);

    auto g_xxxyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 221);

    auto g_xxxyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 222);

    auto g_xxxyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 223);

    auto g_xxxyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 224);

    auto g_xxxyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 225);

    auto g_xxxyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 226);

    auto g_xxxyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 227);

    auto g_xxxyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 228);

    auto g_xxxyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 229);

    auto g_xxxyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 230);

    #pragma omp simd aligned(g_xxxyy_0_xxxxx_0, g_xxxyy_0_xxxxx_1, g_xxxyy_0_xxxxz_0, g_xxxyy_0_xxxxz_1, g_xxxyy_0_xxxzz_0, g_xxxyy_0_xxxzz_1, g_xxxyy_0_xxzzz_0, g_xxxyy_0_xxzzz_1, g_xxxyy_0_xzzzz_0, g_xxxyy_0_xzzzz_1, g_xxxyyy_0_xxxxx_1, g_xxxyyy_0_xxxxz_1, g_xxxyyy_0_xxxzz_1, g_xxxyyy_0_xxzzz_1, g_xxxyyy_0_xzzzz_1, g_xxxyyyy_0_xxxxx_0, g_xxxyyyy_0_xxxxy_0, g_xxxyyyy_0_xxxxz_0, g_xxxyyyy_0_xxxyy_0, g_xxxyyyy_0_xxxyz_0, g_xxxyyyy_0_xxxzz_0, g_xxxyyyy_0_xxyyy_0, g_xxxyyyy_0_xxyyz_0, g_xxxyyyy_0_xxyzz_0, g_xxxyyyy_0_xxzzz_0, g_xxxyyyy_0_xyyyy_0, g_xxxyyyy_0_xyyyz_0, g_xxxyyyy_0_xyyzz_0, g_xxxyyyy_0_xyzzz_0, g_xxxyyyy_0_xzzzz_0, g_xxxyyyy_0_yyyyy_0, g_xxxyyyy_0_yyyyz_0, g_xxxyyyy_0_yyyzz_0, g_xxxyyyy_0_yyzzz_0, g_xxxyyyy_0_yzzzz_0, g_xxxyyyy_0_zzzzz_0, g_xxyyyy_0_xxxxy_1, g_xxyyyy_0_xxxy_1, g_xxyyyy_0_xxxyy_1, g_xxyyyy_0_xxxyz_1, g_xxyyyy_0_xxyy_1, g_xxyyyy_0_xxyyy_1, g_xxyyyy_0_xxyyz_1, g_xxyyyy_0_xxyz_1, g_xxyyyy_0_xxyzz_1, g_xxyyyy_0_xyyy_1, g_xxyyyy_0_xyyyy_1, g_xxyyyy_0_xyyyz_1, g_xxyyyy_0_xyyz_1, g_xxyyyy_0_xyyzz_1, g_xxyyyy_0_xyzz_1, g_xxyyyy_0_xyzzz_1, g_xxyyyy_0_yyyy_1, g_xxyyyy_0_yyyyy_1, g_xxyyyy_0_yyyyz_1, g_xxyyyy_0_yyyz_1, g_xxyyyy_0_yyyzz_1, g_xxyyyy_0_yyzz_1, g_xxyyyy_0_yyzzz_1, g_xxyyyy_0_yzzz_1, g_xxyyyy_0_yzzzz_1, g_xxyyyy_0_zzzzz_1, g_xyyyy_0_xxxxy_0, g_xyyyy_0_xxxxy_1, g_xyyyy_0_xxxyy_0, g_xyyyy_0_xxxyy_1, g_xyyyy_0_xxxyz_0, g_xyyyy_0_xxxyz_1, g_xyyyy_0_xxyyy_0, g_xyyyy_0_xxyyy_1, g_xyyyy_0_xxyyz_0, g_xyyyy_0_xxyyz_1, g_xyyyy_0_xxyzz_0, g_xyyyy_0_xxyzz_1, g_xyyyy_0_xyyyy_0, g_xyyyy_0_xyyyy_1, g_xyyyy_0_xyyyz_0, g_xyyyy_0_xyyyz_1, g_xyyyy_0_xyyzz_0, g_xyyyy_0_xyyzz_1, g_xyyyy_0_xyzzz_0, g_xyyyy_0_xyzzz_1, g_xyyyy_0_yyyyy_0, g_xyyyy_0_yyyyy_1, g_xyyyy_0_yyyyz_0, g_xyyyy_0_yyyyz_1, g_xyyyy_0_yyyzz_0, g_xyyyy_0_yyyzz_1, g_xyyyy_0_yyzzz_0, g_xyyyy_0_yyzzz_1, g_xyyyy_0_yzzzz_0, g_xyyyy_0_yzzzz_1, g_xyyyy_0_zzzzz_0, g_xyyyy_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyyy_0_xxxxx_0[i] = 3.0 * g_xxxyy_0_xxxxx_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxx_1[i] * fz_be_0 + g_xxxyyy_0_xxxxx_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxxy_0[i] = 2.0 * g_xyyyy_0_xxxxy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xxyyyy_0_xxxy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxz_0[i] = 3.0 * g_xxxyy_0_xxxxz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxz_1[i] * fz_be_0 + g_xxxyyy_0_xxxxz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxyy_0[i] = 2.0 * g_xyyyy_0_xxxyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxyz_0[i] = 2.0 * g_xyyyy_0_xxxyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxzz_0[i] = 3.0 * g_xxxyy_0_xxxzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxzz_1[i] * fz_be_0 + g_xxxyyy_0_xxxzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxyyy_0[i] = 2.0 * g_xyyyy_0_xxyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyyz_0[i] = 2.0 * g_xyyyy_0_xxyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyzz_0[i] = 2.0 * g_xyyyy_0_xxyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxzzz_0[i] = 3.0 * g_xxxyy_0_xxzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxzzz_1[i] * fz_be_0 + g_xxxyyy_0_xxzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xyyyy_0[i] = 2.0 * g_xyyyy_0_xyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyy_1[i] * fz_be_0 + g_xxyyyy_0_yyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyyz_0[i] = 2.0 * g_xyyyy_0_xyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyz_1[i] * fz_be_0 + g_xxyyyy_0_yyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyzz_0[i] = 2.0 * g_xyyyy_0_xyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyzz_1[i] * fz_be_0 + g_xxyyyy_0_yyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyzzz_0[i] = 2.0 * g_xyyyy_0_xyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyzzz_1[i] * fz_be_0 + g_xxyyyy_0_yzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xzzzz_0[i] = 3.0 * g_xxxyy_0_xzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xzzzz_1[i] * fz_be_0 + g_xxxyyy_0_xzzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_yyyyy_0[i] = 2.0 * g_xyyyy_0_yyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyy_1[i] * fz_be_0 + g_xxyyyy_0_yyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyyz_0[i] = 2.0 * g_xyyyy_0_yyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyzz_0[i] = 2.0 * g_xyyyy_0_yyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyzzz_0[i] = 2.0 * g_xyyyy_0_yyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yzzzz_0[i] = 2.0 * g_xyyyy_0_yzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_zzzzz_0[i] = 2.0 * g_xyyyy_0_zzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_zzzzz_1[i] * fz_be_0 + g_xxyyyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 231-252 components of targeted buffer : KSH

    auto g_xxxyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 231);

    auto g_xxxyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 232);

    auto g_xxxyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 233);

    auto g_xxxyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 234);

    auto g_xxxyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 235);

    auto g_xxxyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 236);

    auto g_xxxyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 237);

    auto g_xxxyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 238);

    auto g_xxxyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 239);

    auto g_xxxyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 240);

    auto g_xxxyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 241);

    auto g_xxxyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 242);

    auto g_xxxyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 243);

    auto g_xxxyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 244);

    auto g_xxxyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 245);

    auto g_xxxyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 246);

    auto g_xxxyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 247);

    auto g_xxxyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 248);

    auto g_xxxyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 249);

    auto g_xxxyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 250);

    auto g_xxxyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 251);

    #pragma omp simd aligned(g_xxxyyy_0_xxxx_1, g_xxxyyy_0_xxxxx_1, g_xxxyyy_0_xxxxy_1, g_xxxyyy_0_xxxxz_1, g_xxxyyy_0_xxxy_1, g_xxxyyy_0_xxxyy_1, g_xxxyyy_0_xxxyz_1, g_xxxyyy_0_xxxz_1, g_xxxyyy_0_xxxzz_1, g_xxxyyy_0_xxyy_1, g_xxxyyy_0_xxyyy_1, g_xxxyyy_0_xxyyz_1, g_xxxyyy_0_xxyz_1, g_xxxyyy_0_xxyzz_1, g_xxxyyy_0_xxzz_1, g_xxxyyy_0_xxzzz_1, g_xxxyyy_0_xyyy_1, g_xxxyyy_0_xyyyy_1, g_xxxyyy_0_xyyyz_1, g_xxxyyy_0_xyyz_1, g_xxxyyy_0_xyyzz_1, g_xxxyyy_0_xyzz_1, g_xxxyyy_0_xyzzz_1, g_xxxyyy_0_xzzz_1, g_xxxyyy_0_xzzzz_1, g_xxxyyy_0_yyyy_1, g_xxxyyy_0_yyyyy_1, g_xxxyyy_0_yyyyz_1, g_xxxyyy_0_yyyz_1, g_xxxyyy_0_yyyzz_1, g_xxxyyy_0_yyzz_1, g_xxxyyy_0_yyzzz_1, g_xxxyyy_0_yzzz_1, g_xxxyyy_0_yzzzz_1, g_xxxyyy_0_zzzz_1, g_xxxyyy_0_zzzzz_1, g_xxxyyyz_0_xxxxx_0, g_xxxyyyz_0_xxxxy_0, g_xxxyyyz_0_xxxxz_0, g_xxxyyyz_0_xxxyy_0, g_xxxyyyz_0_xxxyz_0, g_xxxyyyz_0_xxxzz_0, g_xxxyyyz_0_xxyyy_0, g_xxxyyyz_0_xxyyz_0, g_xxxyyyz_0_xxyzz_0, g_xxxyyyz_0_xxzzz_0, g_xxxyyyz_0_xyyyy_0, g_xxxyyyz_0_xyyyz_0, g_xxxyyyz_0_xyyzz_0, g_xxxyyyz_0_xyzzz_0, g_xxxyyyz_0_xzzzz_0, g_xxxyyyz_0_yyyyy_0, g_xxxyyyz_0_yyyyz_0, g_xxxyyyz_0_yyyzz_0, g_xxxyyyz_0_yyzzz_0, g_xxxyyyz_0_yzzzz_0, g_xxxyyyz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyyz_0_xxxxx_0[i] = g_xxxyyy_0_xxxxx_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxy_0[i] = g_xxxyyy_0_xxxxy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxz_0[i] = g_xxxyyy_0_xxxx_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyy_0[i] = g_xxxyyy_0_xxxyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyz_0[i] = g_xxxyyy_0_xxxy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxzz_0[i] = 2.0 * g_xxxyyy_0_xxxz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyy_0[i] = g_xxxyyy_0_xxyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyz_0[i] = g_xxxyyy_0_xxyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyzz_0[i] = 2.0 * g_xxxyyy_0_xxyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxzzz_0[i] = 3.0 * g_xxxyyy_0_xxzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyy_0[i] = g_xxxyyy_0_xyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyz_0[i] = g_xxxyyy_0_xyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyzz_0[i] = 2.0 * g_xxxyyy_0_xyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyzzz_0[i] = 3.0 * g_xxxyyy_0_xyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xzzzz_0[i] = 4.0 * g_xxxyyy_0_xzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyy_0[i] = g_xxxyyy_0_yyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyz_0[i] = g_xxxyyy_0_yyyy_1[i] * fi_acd_0 + g_xxxyyy_0_yyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyzz_0[i] = 2.0 * g_xxxyyy_0_yyyz_1[i] * fi_acd_0 + g_xxxyyy_0_yyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyzzz_0[i] = 3.0 * g_xxxyyy_0_yyzz_1[i] * fi_acd_0 + g_xxxyyy_0_yyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yzzzz_0[i] = 4.0 * g_xxxyyy_0_yzzz_1[i] * fi_acd_0 + g_xxxyyy_0_yzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_zzzzz_0[i] = 5.0 * g_xxxyyy_0_zzzz_1[i] * fi_acd_0 + g_xxxyyy_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 252-273 components of targeted buffer : KSH

    auto g_xxxyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 252);

    auto g_xxxyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 253);

    auto g_xxxyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 254);

    auto g_xxxyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 255);

    auto g_xxxyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 256);

    auto g_xxxyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 257);

    auto g_xxxyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 258);

    auto g_xxxyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 259);

    auto g_xxxyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 260);

    auto g_xxxyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 261);

    auto g_xxxyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 262);

    auto g_xxxyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 263);

    auto g_xxxyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 264);

    auto g_xxxyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 265);

    auto g_xxxyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 266);

    auto g_xxxyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 267);

    auto g_xxxyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 268);

    auto g_xxxyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 269);

    auto g_xxxyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 270);

    auto g_xxxyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 271);

    auto g_xxxyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 272);

    #pragma omp simd aligned(g_xxxyy_0_xxxxy_0, g_xxxyy_0_xxxxy_1, g_xxxyy_0_xxxyy_0, g_xxxyy_0_xxxyy_1, g_xxxyy_0_xxyyy_0, g_xxxyy_0_xxyyy_1, g_xxxyy_0_xyyyy_0, g_xxxyy_0_xyyyy_1, g_xxxyyz_0_xxxxy_1, g_xxxyyz_0_xxxyy_1, g_xxxyyz_0_xxyyy_1, g_xxxyyz_0_xyyyy_1, g_xxxyyzz_0_xxxxx_0, g_xxxyyzz_0_xxxxy_0, g_xxxyyzz_0_xxxxz_0, g_xxxyyzz_0_xxxyy_0, g_xxxyyzz_0_xxxyz_0, g_xxxyyzz_0_xxxzz_0, g_xxxyyzz_0_xxyyy_0, g_xxxyyzz_0_xxyyz_0, g_xxxyyzz_0_xxyzz_0, g_xxxyyzz_0_xxzzz_0, g_xxxyyzz_0_xyyyy_0, g_xxxyyzz_0_xyyyz_0, g_xxxyyzz_0_xyyzz_0, g_xxxyyzz_0_xyzzz_0, g_xxxyyzz_0_xzzzz_0, g_xxxyyzz_0_yyyyy_0, g_xxxyyzz_0_yyyyz_0, g_xxxyyzz_0_yyyzz_0, g_xxxyyzz_0_yyzzz_0, g_xxxyyzz_0_yzzzz_0, g_xxxyyzz_0_zzzzz_0, g_xxxyzz_0_xxxxx_1, g_xxxyzz_0_xxxxz_1, g_xxxyzz_0_xxxzz_1, g_xxxyzz_0_xxzzz_1, g_xxxyzz_0_xzzzz_1, g_xxxzz_0_xxxxx_0, g_xxxzz_0_xxxxx_1, g_xxxzz_0_xxxxz_0, g_xxxzz_0_xxxxz_1, g_xxxzz_0_xxxzz_0, g_xxxzz_0_xxxzz_1, g_xxxzz_0_xxzzz_0, g_xxxzz_0_xxzzz_1, g_xxxzz_0_xzzzz_0, g_xxxzz_0_xzzzz_1, g_xxyyzz_0_xxxyz_1, g_xxyyzz_0_xxyyz_1, g_xxyyzz_0_xxyz_1, g_xxyyzz_0_xxyzz_1, g_xxyyzz_0_xyyyz_1, g_xxyyzz_0_xyyz_1, g_xxyyzz_0_xyyzz_1, g_xxyyzz_0_xyzz_1, g_xxyyzz_0_xyzzz_1, g_xxyyzz_0_yyyyy_1, g_xxyyzz_0_yyyyz_1, g_xxyyzz_0_yyyz_1, g_xxyyzz_0_yyyzz_1, g_xxyyzz_0_yyzz_1, g_xxyyzz_0_yyzzz_1, g_xxyyzz_0_yzzz_1, g_xxyyzz_0_yzzzz_1, g_xxyyzz_0_zzzzz_1, g_xyyzz_0_xxxyz_0, g_xyyzz_0_xxxyz_1, g_xyyzz_0_xxyyz_0, g_xyyzz_0_xxyyz_1, g_xyyzz_0_xxyzz_0, g_xyyzz_0_xxyzz_1, g_xyyzz_0_xyyyz_0, g_xyyzz_0_xyyyz_1, g_xyyzz_0_xyyzz_0, g_xyyzz_0_xyyzz_1, g_xyyzz_0_xyzzz_0, g_xyyzz_0_xyzzz_1, g_xyyzz_0_yyyyy_0, g_xyyzz_0_yyyyy_1, g_xyyzz_0_yyyyz_0, g_xyyzz_0_yyyyz_1, g_xyyzz_0_yyyzz_0, g_xyyzz_0_yyyzz_1, g_xyyzz_0_yyzzz_0, g_xyyzz_0_yyzzz_1, g_xyyzz_0_yzzzz_0, g_xyyzz_0_yzzzz_1, g_xyyzz_0_zzzzz_0, g_xyyzz_0_zzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyzz_0_xxxxx_0[i] = g_xxxzz_0_xxxxx_0[i] * fbe_0 - g_xxxzz_0_xxxxx_1[i] * fz_be_0 + g_xxxyzz_0_xxxxx_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxxy_0[i] = g_xxxyy_0_xxxxy_0[i] * fbe_0 - g_xxxyy_0_xxxxy_1[i] * fz_be_0 + g_xxxyyz_0_xxxxy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxxz_0[i] = g_xxxzz_0_xxxxz_0[i] * fbe_0 - g_xxxzz_0_xxxxz_1[i] * fz_be_0 + g_xxxyzz_0_xxxxz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxyy_0[i] = g_xxxyy_0_xxxyy_0[i] * fbe_0 - g_xxxyy_0_xxxyy_1[i] * fz_be_0 + g_xxxyyz_0_xxxyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxyz_0[i] = 2.0 * g_xyyzz_0_xxxyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_0_xxyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxzz_0[i] = g_xxxzz_0_xxxzz_0[i] * fbe_0 - g_xxxzz_0_xxxzz_1[i] * fz_be_0 + g_xxxyzz_0_xxxzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxyyy_0[i] = g_xxxyy_0_xxyyy_0[i] * fbe_0 - g_xxxyy_0_xxyyy_1[i] * fz_be_0 + g_xxxyyz_0_xxyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxyyz_0[i] = 2.0 * g_xyyzz_0_xxyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxyzz_0[i] = 2.0 * g_xyyzz_0_xxyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxzzz_0[i] = g_xxxzz_0_xxzzz_0[i] * fbe_0 - g_xxxzz_0_xxzzz_1[i] * fz_be_0 + g_xxxyzz_0_xxzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xyyyy_0[i] = g_xxxyy_0_xyyyy_0[i] * fbe_0 - g_xxxyy_0_xyyyy_1[i] * fz_be_0 + g_xxxyyz_0_xyyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xyyyz_0[i] = 2.0 * g_xyyzz_0_xyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyyz_1[i] * fz_be_0 + g_xxyyzz_0_yyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyyzz_0[i] = 2.0 * g_xyyzz_0_xyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyzz_1[i] * fz_be_0 + g_xxyyzz_0_yyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyzzz_0[i] = 2.0 * g_xyyzz_0_xyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyzzz_1[i] * fz_be_0 + g_xxyyzz_0_yzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xzzzz_0[i] = g_xxxzz_0_xzzzz_0[i] * fbe_0 - g_xxxzz_0_xzzzz_1[i] * fz_be_0 + g_xxxyzz_0_xzzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_yyyyy_0[i] = 2.0 * g_xyyzz_0_yyyyy_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyy_1[i] * fz_be_0 + g_xxyyzz_0_yyyyy_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyyz_0[i] = 2.0 * g_xyyzz_0_yyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyzz_0[i] = 2.0 * g_xyyzz_0_yyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyzzz_0[i] = 2.0 * g_xyyzz_0_yyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yzzzz_0[i] = 2.0 * g_xyyzz_0_yzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_zzzzz_0[i] = 2.0 * g_xyyzz_0_zzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_zzzzz_1[i] * fz_be_0 + g_xxyyzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 273-294 components of targeted buffer : KSH

    auto g_xxxyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 273);

    auto g_xxxyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 274);

    auto g_xxxyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 275);

    auto g_xxxyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 276);

    auto g_xxxyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 277);

    auto g_xxxyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 278);

    auto g_xxxyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 279);

    auto g_xxxyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 280);

    auto g_xxxyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 281);

    auto g_xxxyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 282);

    auto g_xxxyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 283);

    auto g_xxxyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 284);

    auto g_xxxyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 285);

    auto g_xxxyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 286);

    auto g_xxxyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 287);

    auto g_xxxyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 288);

    auto g_xxxyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 289);

    auto g_xxxyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 290);

    auto g_xxxyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 291);

    auto g_xxxyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 292);

    auto g_xxxyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 293);

    #pragma omp simd aligned(g_xxxyzzz_0_xxxxx_0, g_xxxyzzz_0_xxxxy_0, g_xxxyzzz_0_xxxxz_0, g_xxxyzzz_0_xxxyy_0, g_xxxyzzz_0_xxxyz_0, g_xxxyzzz_0_xxxzz_0, g_xxxyzzz_0_xxyyy_0, g_xxxyzzz_0_xxyyz_0, g_xxxyzzz_0_xxyzz_0, g_xxxyzzz_0_xxzzz_0, g_xxxyzzz_0_xyyyy_0, g_xxxyzzz_0_xyyyz_0, g_xxxyzzz_0_xyyzz_0, g_xxxyzzz_0_xyzzz_0, g_xxxyzzz_0_xzzzz_0, g_xxxyzzz_0_yyyyy_0, g_xxxyzzz_0_yyyyz_0, g_xxxyzzz_0_yyyzz_0, g_xxxyzzz_0_yyzzz_0, g_xxxyzzz_0_yzzzz_0, g_xxxyzzz_0_zzzzz_0, g_xxxzzz_0_xxxx_1, g_xxxzzz_0_xxxxx_1, g_xxxzzz_0_xxxxy_1, g_xxxzzz_0_xxxxz_1, g_xxxzzz_0_xxxy_1, g_xxxzzz_0_xxxyy_1, g_xxxzzz_0_xxxyz_1, g_xxxzzz_0_xxxz_1, g_xxxzzz_0_xxxzz_1, g_xxxzzz_0_xxyy_1, g_xxxzzz_0_xxyyy_1, g_xxxzzz_0_xxyyz_1, g_xxxzzz_0_xxyz_1, g_xxxzzz_0_xxyzz_1, g_xxxzzz_0_xxzz_1, g_xxxzzz_0_xxzzz_1, g_xxxzzz_0_xyyy_1, g_xxxzzz_0_xyyyy_1, g_xxxzzz_0_xyyyz_1, g_xxxzzz_0_xyyz_1, g_xxxzzz_0_xyyzz_1, g_xxxzzz_0_xyzz_1, g_xxxzzz_0_xyzzz_1, g_xxxzzz_0_xzzz_1, g_xxxzzz_0_xzzzz_1, g_xxxzzz_0_yyyy_1, g_xxxzzz_0_yyyyy_1, g_xxxzzz_0_yyyyz_1, g_xxxzzz_0_yyyz_1, g_xxxzzz_0_yyyzz_1, g_xxxzzz_0_yyzz_1, g_xxxzzz_0_yyzzz_1, g_xxxzzz_0_yzzz_1, g_xxxzzz_0_yzzzz_1, g_xxxzzz_0_zzzz_1, g_xxxzzz_0_zzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzzz_0_xxxxx_0[i] = g_xxxzzz_0_xxxxx_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxy_0[i] = g_xxxzzz_0_xxxx_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxz_0[i] = g_xxxzzz_0_xxxxz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyy_0[i] = 2.0 * g_xxxzzz_0_xxxy_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyz_0[i] = g_xxxzzz_0_xxxz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxzz_0[i] = g_xxxzzz_0_xxxzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyy_0[i] = 3.0 * g_xxxzzz_0_xxyy_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyz_0[i] = 2.0 * g_xxxzzz_0_xxyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyzz_0[i] = g_xxxzzz_0_xxzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxzzz_0[i] = g_xxxzzz_0_xxzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyy_0[i] = 4.0 * g_xxxzzz_0_xyyy_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyz_0[i] = 3.0 * g_xxxzzz_0_xyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyzz_0[i] = 2.0 * g_xxxzzz_0_xyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyzzz_0[i] = g_xxxzzz_0_xzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xzzzz_0[i] = g_xxxzzz_0_xzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyy_0[i] = 5.0 * g_xxxzzz_0_yyyy_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyz_0[i] = 4.0 * g_xxxzzz_0_yyyz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyzz_0[i] = 3.0 * g_xxxzzz_0_yyzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyzzz_0[i] = 2.0 * g_xxxzzz_0_yzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yzzzz_0[i] = g_xxxzzz_0_zzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_zzzzz_0[i] = g_xxxzzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 294-315 components of targeted buffer : KSH

    auto g_xxxzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 294);

    auto g_xxxzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 295);

    auto g_xxxzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 296);

    auto g_xxxzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 297);

    auto g_xxxzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 298);

    auto g_xxxzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 299);

    auto g_xxxzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 300);

    auto g_xxxzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 301);

    auto g_xxxzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 302);

    auto g_xxxzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 303);

    auto g_xxxzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 304);

    auto g_xxxzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 305);

    auto g_xxxzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 306);

    auto g_xxxzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 307);

    auto g_xxxzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 308);

    auto g_xxxzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 309);

    auto g_xxxzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 310);

    auto g_xxxzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 311);

    auto g_xxxzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 312);

    auto g_xxxzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 313);

    auto g_xxxzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 314);

    #pragma omp simd aligned(g_xxxzz_0_xxxxx_0, g_xxxzz_0_xxxxx_1, g_xxxzz_0_xxxxy_0, g_xxxzz_0_xxxxy_1, g_xxxzz_0_xxxyy_0, g_xxxzz_0_xxxyy_1, g_xxxzz_0_xxyyy_0, g_xxxzz_0_xxyyy_1, g_xxxzz_0_xyyyy_0, g_xxxzz_0_xyyyy_1, g_xxxzzz_0_xxxxx_1, g_xxxzzz_0_xxxxy_1, g_xxxzzz_0_xxxyy_1, g_xxxzzz_0_xxyyy_1, g_xxxzzz_0_xyyyy_1, g_xxxzzzz_0_xxxxx_0, g_xxxzzzz_0_xxxxy_0, g_xxxzzzz_0_xxxxz_0, g_xxxzzzz_0_xxxyy_0, g_xxxzzzz_0_xxxyz_0, g_xxxzzzz_0_xxxzz_0, g_xxxzzzz_0_xxyyy_0, g_xxxzzzz_0_xxyyz_0, g_xxxzzzz_0_xxyzz_0, g_xxxzzzz_0_xxzzz_0, g_xxxzzzz_0_xyyyy_0, g_xxxzzzz_0_xyyyz_0, g_xxxzzzz_0_xyyzz_0, g_xxxzzzz_0_xyzzz_0, g_xxxzzzz_0_xzzzz_0, g_xxxzzzz_0_yyyyy_0, g_xxxzzzz_0_yyyyz_0, g_xxxzzzz_0_yyyzz_0, g_xxxzzzz_0_yyzzz_0, g_xxxzzzz_0_yzzzz_0, g_xxxzzzz_0_zzzzz_0, g_xxzzzz_0_xxxxz_1, g_xxzzzz_0_xxxyz_1, g_xxzzzz_0_xxxz_1, g_xxzzzz_0_xxxzz_1, g_xxzzzz_0_xxyyz_1, g_xxzzzz_0_xxyz_1, g_xxzzzz_0_xxyzz_1, g_xxzzzz_0_xxzz_1, g_xxzzzz_0_xxzzz_1, g_xxzzzz_0_xyyyz_1, g_xxzzzz_0_xyyz_1, g_xxzzzz_0_xyyzz_1, g_xxzzzz_0_xyzz_1, g_xxzzzz_0_xyzzz_1, g_xxzzzz_0_xzzz_1, g_xxzzzz_0_xzzzz_1, g_xxzzzz_0_yyyyy_1, g_xxzzzz_0_yyyyz_1, g_xxzzzz_0_yyyz_1, g_xxzzzz_0_yyyzz_1, g_xxzzzz_0_yyzz_1, g_xxzzzz_0_yyzzz_1, g_xxzzzz_0_yzzz_1, g_xxzzzz_0_yzzzz_1, g_xxzzzz_0_zzzz_1, g_xxzzzz_0_zzzzz_1, g_xzzzz_0_xxxxz_0, g_xzzzz_0_xxxxz_1, g_xzzzz_0_xxxyz_0, g_xzzzz_0_xxxyz_1, g_xzzzz_0_xxxzz_0, g_xzzzz_0_xxxzz_1, g_xzzzz_0_xxyyz_0, g_xzzzz_0_xxyyz_1, g_xzzzz_0_xxyzz_0, g_xzzzz_0_xxyzz_1, g_xzzzz_0_xxzzz_0, g_xzzzz_0_xxzzz_1, g_xzzzz_0_xyyyz_0, g_xzzzz_0_xyyyz_1, g_xzzzz_0_xyyzz_0, g_xzzzz_0_xyyzz_1, g_xzzzz_0_xyzzz_0, g_xzzzz_0_xyzzz_1, g_xzzzz_0_xzzzz_0, g_xzzzz_0_xzzzz_1, g_xzzzz_0_yyyyy_0, g_xzzzz_0_yyyyy_1, g_xzzzz_0_yyyyz_0, g_xzzzz_0_yyyyz_1, g_xzzzz_0_yyyzz_0, g_xzzzz_0_yyyzz_1, g_xzzzz_0_yyzzz_0, g_xzzzz_0_yyzzz_1, g_xzzzz_0_yzzzz_0, g_xzzzz_0_yzzzz_1, g_xzzzz_0_zzzzz_0, g_xzzzz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzzz_0_xxxxx_0[i] = 3.0 * g_xxxzz_0_xxxxx_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxx_1[i] * fz_be_0 + g_xxxzzz_0_xxxxx_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxy_0[i] = 3.0 * g_xxxzz_0_xxxxy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxy_1[i] * fz_be_0 + g_xxxzzz_0_xxxxy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxz_0[i] = 2.0 * g_xzzzz_0_xxxxz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_0_xxxz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxyy_0[i] = 3.0 * g_xxxzz_0_xxxyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxyy_1[i] * fz_be_0 + g_xxxzzz_0_xxxyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxyz_0[i] = 2.0 * g_xzzzz_0_xxxyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxzz_0[i] = 2.0 * g_xzzzz_0_xxxzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyyy_0[i] = 3.0 * g_xxxzz_0_xxyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxyyy_1[i] * fz_be_0 + g_xxxzzz_0_xxyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxyyz_0[i] = 2.0 * g_xzzzz_0_xxyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyzz_0[i] = 2.0 * g_xzzzz_0_xxyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxzzz_0[i] = 2.0 * g_xzzzz_0_xxzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyyy_0[i] = 3.0 * g_xxxzz_0_xyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xyyyy_1[i] * fz_be_0 + g_xxxzzz_0_xyyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xyyyz_0[i] = 2.0 * g_xzzzz_0_xyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyyz_1[i] * fz_be_0 + g_xxzzzz_0_yyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyzz_0[i] = 2.0 * g_xzzzz_0_xyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyzz_1[i] * fz_be_0 + g_xxzzzz_0_yyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyzzz_0[i] = 2.0 * g_xzzzz_0_xyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyzzz_1[i] * fz_be_0 + g_xxzzzz_0_yzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xzzzz_0[i] = 2.0 * g_xzzzz_0_xzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xzzzz_1[i] * fz_be_0 + g_xxzzzz_0_zzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyy_0[i] = 2.0 * g_xzzzz_0_yyyyy_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyy_1[i] * fz_be_0 + g_xxzzzz_0_yyyyy_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyz_0[i] = 2.0 * g_xzzzz_0_yyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyzz_0[i] = 2.0 * g_xzzzz_0_yyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyzzz_0[i] = 2.0 * g_xzzzz_0_yyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yzzzz_0[i] = 2.0 * g_xzzzz_0_yzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_zzzzz_0[i] = 2.0 * g_xzzzz_0_zzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_zzzzz_1[i] * fz_be_0 + g_xxzzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 315-336 components of targeted buffer : KSH

    auto g_xxyyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 315);

    auto g_xxyyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 316);

    auto g_xxyyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 317);

    auto g_xxyyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 318);

    auto g_xxyyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 319);

    auto g_xxyyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 320);

    auto g_xxyyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 321);

    auto g_xxyyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 322);

    auto g_xxyyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 323);

    auto g_xxyyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 324);

    auto g_xxyyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 325);

    auto g_xxyyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 326);

    auto g_xxyyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 327);

    auto g_xxyyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 328);

    auto g_xxyyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 329);

    auto g_xxyyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 330);

    auto g_xxyyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 331);

    auto g_xxyyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 332);

    auto g_xxyyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 333);

    auto g_xxyyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 334);

    auto g_xxyyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 335);

    #pragma omp simd aligned(g_xxyyy_0_xxxxx_0, g_xxyyy_0_xxxxx_1, g_xxyyy_0_xxxxz_0, g_xxyyy_0_xxxxz_1, g_xxyyy_0_xxxzz_0, g_xxyyy_0_xxxzz_1, g_xxyyy_0_xxzzz_0, g_xxyyy_0_xxzzz_1, g_xxyyy_0_xzzzz_0, g_xxyyy_0_xzzzz_1, g_xxyyyy_0_xxxxx_1, g_xxyyyy_0_xxxxz_1, g_xxyyyy_0_xxxzz_1, g_xxyyyy_0_xxzzz_1, g_xxyyyy_0_xzzzz_1, g_xxyyyyy_0_xxxxx_0, g_xxyyyyy_0_xxxxy_0, g_xxyyyyy_0_xxxxz_0, g_xxyyyyy_0_xxxyy_0, g_xxyyyyy_0_xxxyz_0, g_xxyyyyy_0_xxxzz_0, g_xxyyyyy_0_xxyyy_0, g_xxyyyyy_0_xxyyz_0, g_xxyyyyy_0_xxyzz_0, g_xxyyyyy_0_xxzzz_0, g_xxyyyyy_0_xyyyy_0, g_xxyyyyy_0_xyyyz_0, g_xxyyyyy_0_xyyzz_0, g_xxyyyyy_0_xyzzz_0, g_xxyyyyy_0_xzzzz_0, g_xxyyyyy_0_yyyyy_0, g_xxyyyyy_0_yyyyz_0, g_xxyyyyy_0_yyyzz_0, g_xxyyyyy_0_yyzzz_0, g_xxyyyyy_0_yzzzz_0, g_xxyyyyy_0_zzzzz_0, g_xyyyyy_0_xxxxy_1, g_xyyyyy_0_xxxy_1, g_xyyyyy_0_xxxyy_1, g_xyyyyy_0_xxxyz_1, g_xyyyyy_0_xxyy_1, g_xyyyyy_0_xxyyy_1, g_xyyyyy_0_xxyyz_1, g_xyyyyy_0_xxyz_1, g_xyyyyy_0_xxyzz_1, g_xyyyyy_0_xyyy_1, g_xyyyyy_0_xyyyy_1, g_xyyyyy_0_xyyyz_1, g_xyyyyy_0_xyyz_1, g_xyyyyy_0_xyyzz_1, g_xyyyyy_0_xyzz_1, g_xyyyyy_0_xyzzz_1, g_xyyyyy_0_yyyy_1, g_xyyyyy_0_yyyyy_1, g_xyyyyy_0_yyyyz_1, g_xyyyyy_0_yyyz_1, g_xyyyyy_0_yyyzz_1, g_xyyyyy_0_yyzz_1, g_xyyyyy_0_yyzzz_1, g_xyyyyy_0_yzzz_1, g_xyyyyy_0_yzzzz_1, g_xyyyyy_0_zzzzz_1, g_yyyyy_0_xxxxy_0, g_yyyyy_0_xxxxy_1, g_yyyyy_0_xxxyy_0, g_yyyyy_0_xxxyy_1, g_yyyyy_0_xxxyz_0, g_yyyyy_0_xxxyz_1, g_yyyyy_0_xxyyy_0, g_yyyyy_0_xxyyy_1, g_yyyyy_0_xxyyz_0, g_yyyyy_0_xxyyz_1, g_yyyyy_0_xxyzz_0, g_yyyyy_0_xxyzz_1, g_yyyyy_0_xyyyy_0, g_yyyyy_0_xyyyy_1, g_yyyyy_0_xyyyz_0, g_yyyyy_0_xyyyz_1, g_yyyyy_0_xyyzz_0, g_yyyyy_0_xyyzz_1, g_yyyyy_0_xyzzz_0, g_yyyyy_0_xyzzz_1, g_yyyyy_0_yyyyy_0, g_yyyyy_0_yyyyy_1, g_yyyyy_0_yyyyz_0, g_yyyyy_0_yyyyz_1, g_yyyyy_0_yyyzz_0, g_yyyyy_0_yyyzz_1, g_yyyyy_0_yyzzz_0, g_yyyyy_0_yyzzz_1, g_yyyyy_0_yzzzz_0, g_yyyyy_0_yzzzz_1, g_yyyyy_0_zzzzz_0, g_yyyyy_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyyy_0_xxxxx_0[i] = 4.0 * g_xxyyy_0_xxxxx_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxx_1[i] * fz_be_0 + g_xxyyyy_0_xxxxx_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxxy_0[i] = g_yyyyy_0_xxxxy_0[i] * fbe_0 - g_yyyyy_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_xyyyyy_0_xxxy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxz_0[i] = 4.0 * g_xxyyy_0_xxxxz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxz_1[i] * fz_be_0 + g_xxyyyy_0_xxxxz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxyy_0[i] = g_yyyyy_0_xxxyy_0[i] * fbe_0 - g_yyyyy_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxyz_0[i] = g_yyyyy_0_xxxyz_0[i] * fbe_0 - g_yyyyy_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxzz_0[i] = 4.0 * g_xxyyy_0_xxxzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxzz_1[i] * fz_be_0 + g_xxyyyy_0_xxxzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxyyy_0[i] = g_yyyyy_0_xxyyy_0[i] * fbe_0 - g_yyyyy_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyyz_0[i] = g_yyyyy_0_xxyyz_0[i] * fbe_0 - g_yyyyy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyzz_0[i] = g_yyyyy_0_xxyzz_0[i] * fbe_0 - g_yyyyy_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxzzz_0[i] = 4.0 * g_xxyyy_0_xxzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxzzz_1[i] * fz_be_0 + g_xxyyyy_0_xxzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xyyyy_0[i] = g_yyyyy_0_xyyyy_0[i] * fbe_0 - g_yyyyy_0_xyyyy_1[i] * fz_be_0 + g_xyyyyy_0_yyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyyz_0[i] = g_yyyyy_0_xyyyz_0[i] * fbe_0 - g_yyyyy_0_xyyyz_1[i] * fz_be_0 + g_xyyyyy_0_yyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyzz_0[i] = g_yyyyy_0_xyyzz_0[i] * fbe_0 - g_yyyyy_0_xyyzz_1[i] * fz_be_0 + g_xyyyyy_0_yyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyzzz_0[i] = g_yyyyy_0_xyzzz_0[i] * fbe_0 - g_yyyyy_0_xyzzz_1[i] * fz_be_0 + g_xyyyyy_0_yzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xzzzz_0[i] = 4.0 * g_xxyyy_0_xzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xzzzz_1[i] * fz_be_0 + g_xxyyyy_0_xzzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_yyyyy_0[i] = g_yyyyy_0_yyyyy_0[i] * fbe_0 - g_yyyyy_0_yyyyy_1[i] * fz_be_0 + g_xyyyyy_0_yyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyyz_0[i] = g_yyyyy_0_yyyyz_0[i] * fbe_0 - g_yyyyy_0_yyyyz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyzz_0[i] = g_yyyyy_0_yyyzz_0[i] * fbe_0 - g_yyyyy_0_yyyzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyzzz_0[i] = g_yyyyy_0_yyzzz_0[i] * fbe_0 - g_yyyyy_0_yyzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yzzzz_0[i] = g_yyyyy_0_yzzzz_0[i] * fbe_0 - g_yyyyy_0_yzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_zzzzz_0[i] = g_yyyyy_0_zzzzz_0[i] * fbe_0 - g_yyyyy_0_zzzzz_1[i] * fz_be_0 + g_xyyyyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 336-357 components of targeted buffer : KSH

    auto g_xxyyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 336);

    auto g_xxyyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 337);

    auto g_xxyyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 338);

    auto g_xxyyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 339);

    auto g_xxyyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 340);

    auto g_xxyyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 341);

    auto g_xxyyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 342);

    auto g_xxyyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 343);

    auto g_xxyyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 344);

    auto g_xxyyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 345);

    auto g_xxyyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 346);

    auto g_xxyyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 347);

    auto g_xxyyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 348);

    auto g_xxyyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 349);

    auto g_xxyyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 350);

    auto g_xxyyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 351);

    auto g_xxyyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 352);

    auto g_xxyyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 353);

    auto g_xxyyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 354);

    auto g_xxyyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 355);

    auto g_xxyyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 356);

    #pragma omp simd aligned(g_xxyyyy_0_xxxx_1, g_xxyyyy_0_xxxxx_1, g_xxyyyy_0_xxxxy_1, g_xxyyyy_0_xxxxz_1, g_xxyyyy_0_xxxy_1, g_xxyyyy_0_xxxyy_1, g_xxyyyy_0_xxxyz_1, g_xxyyyy_0_xxxz_1, g_xxyyyy_0_xxxzz_1, g_xxyyyy_0_xxyy_1, g_xxyyyy_0_xxyyy_1, g_xxyyyy_0_xxyyz_1, g_xxyyyy_0_xxyz_1, g_xxyyyy_0_xxyzz_1, g_xxyyyy_0_xxzz_1, g_xxyyyy_0_xxzzz_1, g_xxyyyy_0_xyyy_1, g_xxyyyy_0_xyyyy_1, g_xxyyyy_0_xyyyz_1, g_xxyyyy_0_xyyz_1, g_xxyyyy_0_xyyzz_1, g_xxyyyy_0_xyzz_1, g_xxyyyy_0_xyzzz_1, g_xxyyyy_0_xzzz_1, g_xxyyyy_0_xzzzz_1, g_xxyyyy_0_yyyy_1, g_xxyyyy_0_yyyyy_1, g_xxyyyy_0_yyyyz_1, g_xxyyyy_0_yyyz_1, g_xxyyyy_0_yyyzz_1, g_xxyyyy_0_yyzz_1, g_xxyyyy_0_yyzzz_1, g_xxyyyy_0_yzzz_1, g_xxyyyy_0_yzzzz_1, g_xxyyyy_0_zzzz_1, g_xxyyyy_0_zzzzz_1, g_xxyyyyz_0_xxxxx_0, g_xxyyyyz_0_xxxxy_0, g_xxyyyyz_0_xxxxz_0, g_xxyyyyz_0_xxxyy_0, g_xxyyyyz_0_xxxyz_0, g_xxyyyyz_0_xxxzz_0, g_xxyyyyz_0_xxyyy_0, g_xxyyyyz_0_xxyyz_0, g_xxyyyyz_0_xxyzz_0, g_xxyyyyz_0_xxzzz_0, g_xxyyyyz_0_xyyyy_0, g_xxyyyyz_0_xyyyz_0, g_xxyyyyz_0_xyyzz_0, g_xxyyyyz_0_xyzzz_0, g_xxyyyyz_0_xzzzz_0, g_xxyyyyz_0_yyyyy_0, g_xxyyyyz_0_yyyyz_0, g_xxyyyyz_0_yyyzz_0, g_xxyyyyz_0_yyzzz_0, g_xxyyyyz_0_yzzzz_0, g_xxyyyyz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyyz_0_xxxxx_0[i] = g_xxyyyy_0_xxxxx_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxy_0[i] = g_xxyyyy_0_xxxxy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxz_0[i] = g_xxyyyy_0_xxxx_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyy_0[i] = g_xxyyyy_0_xxxyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyz_0[i] = g_xxyyyy_0_xxxy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxzz_0[i] = 2.0 * g_xxyyyy_0_xxxz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyy_0[i] = g_xxyyyy_0_xxyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyz_0[i] = g_xxyyyy_0_xxyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyzz_0[i] = 2.0 * g_xxyyyy_0_xxyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxzzz_0[i] = 3.0 * g_xxyyyy_0_xxzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyy_0[i] = g_xxyyyy_0_xyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyz_0[i] = g_xxyyyy_0_xyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyzz_0[i] = 2.0 * g_xxyyyy_0_xyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyzzz_0[i] = 3.0 * g_xxyyyy_0_xyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xzzzz_0[i] = 4.0 * g_xxyyyy_0_xzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyy_0[i] = g_xxyyyy_0_yyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyz_0[i] = g_xxyyyy_0_yyyy_1[i] * fi_acd_0 + g_xxyyyy_0_yyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyzz_0[i] = 2.0 * g_xxyyyy_0_yyyz_1[i] * fi_acd_0 + g_xxyyyy_0_yyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyzzz_0[i] = 3.0 * g_xxyyyy_0_yyzz_1[i] * fi_acd_0 + g_xxyyyy_0_yyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yzzzz_0[i] = 4.0 * g_xxyyyy_0_yzzz_1[i] * fi_acd_0 + g_xxyyyy_0_yzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_zzzzz_0[i] = 5.0 * g_xxyyyy_0_zzzz_1[i] * fi_acd_0 + g_xxyyyy_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 357-378 components of targeted buffer : KSH

    auto g_xxyyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 357);

    auto g_xxyyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 358);

    auto g_xxyyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 359);

    auto g_xxyyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 360);

    auto g_xxyyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 361);

    auto g_xxyyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 362);

    auto g_xxyyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 363);

    auto g_xxyyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 364);

    auto g_xxyyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 365);

    auto g_xxyyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 366);

    auto g_xxyyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 367);

    auto g_xxyyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 368);

    auto g_xxyyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 369);

    auto g_xxyyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 370);

    auto g_xxyyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 371);

    auto g_xxyyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 372);

    auto g_xxyyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 373);

    auto g_xxyyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 374);

    auto g_xxyyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 375);

    auto g_xxyyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 376);

    auto g_xxyyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 377);

    #pragma omp simd aligned(g_xxyyy_0_xxxxy_0, g_xxyyy_0_xxxxy_1, g_xxyyy_0_xxxyy_0, g_xxyyy_0_xxxyy_1, g_xxyyy_0_xxyyy_0, g_xxyyy_0_xxyyy_1, g_xxyyy_0_xyyyy_0, g_xxyyy_0_xyyyy_1, g_xxyyyz_0_xxxxy_1, g_xxyyyz_0_xxxyy_1, g_xxyyyz_0_xxyyy_1, g_xxyyyz_0_xyyyy_1, g_xxyyyzz_0_xxxxx_0, g_xxyyyzz_0_xxxxy_0, g_xxyyyzz_0_xxxxz_0, g_xxyyyzz_0_xxxyy_0, g_xxyyyzz_0_xxxyz_0, g_xxyyyzz_0_xxxzz_0, g_xxyyyzz_0_xxyyy_0, g_xxyyyzz_0_xxyyz_0, g_xxyyyzz_0_xxyzz_0, g_xxyyyzz_0_xxzzz_0, g_xxyyyzz_0_xyyyy_0, g_xxyyyzz_0_xyyyz_0, g_xxyyyzz_0_xyyzz_0, g_xxyyyzz_0_xyzzz_0, g_xxyyyzz_0_xzzzz_0, g_xxyyyzz_0_yyyyy_0, g_xxyyyzz_0_yyyyz_0, g_xxyyyzz_0_yyyzz_0, g_xxyyyzz_0_yyzzz_0, g_xxyyyzz_0_yzzzz_0, g_xxyyyzz_0_zzzzz_0, g_xxyyzz_0_xxxxx_1, g_xxyyzz_0_xxxxz_1, g_xxyyzz_0_xxxzz_1, g_xxyyzz_0_xxzzz_1, g_xxyyzz_0_xzzzz_1, g_xxyzz_0_xxxxx_0, g_xxyzz_0_xxxxx_1, g_xxyzz_0_xxxxz_0, g_xxyzz_0_xxxxz_1, g_xxyzz_0_xxxzz_0, g_xxyzz_0_xxxzz_1, g_xxyzz_0_xxzzz_0, g_xxyzz_0_xxzzz_1, g_xxyzz_0_xzzzz_0, g_xxyzz_0_xzzzz_1, g_xyyyzz_0_xxxyz_1, g_xyyyzz_0_xxyyz_1, g_xyyyzz_0_xxyz_1, g_xyyyzz_0_xxyzz_1, g_xyyyzz_0_xyyyz_1, g_xyyyzz_0_xyyz_1, g_xyyyzz_0_xyyzz_1, g_xyyyzz_0_xyzz_1, g_xyyyzz_0_xyzzz_1, g_xyyyzz_0_yyyyy_1, g_xyyyzz_0_yyyyz_1, g_xyyyzz_0_yyyz_1, g_xyyyzz_0_yyyzz_1, g_xyyyzz_0_yyzz_1, g_xyyyzz_0_yyzzz_1, g_xyyyzz_0_yzzz_1, g_xyyyzz_0_yzzzz_1, g_xyyyzz_0_zzzzz_1, g_yyyzz_0_xxxyz_0, g_yyyzz_0_xxxyz_1, g_yyyzz_0_xxyyz_0, g_yyyzz_0_xxyyz_1, g_yyyzz_0_xxyzz_0, g_yyyzz_0_xxyzz_1, g_yyyzz_0_xyyyz_0, g_yyyzz_0_xyyyz_1, g_yyyzz_0_xyyzz_0, g_yyyzz_0_xyyzz_1, g_yyyzz_0_xyzzz_0, g_yyyzz_0_xyzzz_1, g_yyyzz_0_yyyyy_0, g_yyyzz_0_yyyyy_1, g_yyyzz_0_yyyyz_0, g_yyyzz_0_yyyyz_1, g_yyyzz_0_yyyzz_0, g_yyyzz_0_yyyzz_1, g_yyyzz_0_yyzzz_0, g_yyyzz_0_yyzzz_1, g_yyyzz_0_yzzzz_0, g_yyyzz_0_yzzzz_1, g_yyyzz_0_zzzzz_0, g_yyyzz_0_zzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyzz_0_xxxxx_0[i] = 2.0 * g_xxyzz_0_xxxxx_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxx_1[i] * fz_be_0 + g_xxyyzz_0_xxxxx_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxxy_0[i] = g_xxyyy_0_xxxxy_0[i] * fbe_0 - g_xxyyy_0_xxxxy_1[i] * fz_be_0 + g_xxyyyz_0_xxxxy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxxz_0[i] = 2.0 * g_xxyzz_0_xxxxz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxz_1[i] * fz_be_0 + g_xxyyzz_0_xxxxz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxyy_0[i] = g_xxyyy_0_xxxyy_0[i] * fbe_0 - g_xxyyy_0_xxxyy_1[i] * fz_be_0 + g_xxyyyz_0_xxxyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxyz_0[i] = g_yyyzz_0_xxxyz_0[i] * fbe_0 - g_yyyzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_0_xxyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxzz_0[i] = 2.0 * g_xxyzz_0_xxxzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxzz_1[i] * fz_be_0 + g_xxyyzz_0_xxxzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxyyy_0[i] = g_xxyyy_0_xxyyy_0[i] * fbe_0 - g_xxyyy_0_xxyyy_1[i] * fz_be_0 + g_xxyyyz_0_xxyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxyyz_0[i] = g_yyyzz_0_xxyyz_0[i] * fbe_0 - g_yyyzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxyzz_0[i] = g_yyyzz_0_xxyzz_0[i] * fbe_0 - g_yyyzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxzzz_0[i] = 2.0 * g_xxyzz_0_xxzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxzzz_1[i] * fz_be_0 + g_xxyyzz_0_xxzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xyyyy_0[i] = g_xxyyy_0_xyyyy_0[i] * fbe_0 - g_xxyyy_0_xyyyy_1[i] * fz_be_0 + g_xxyyyz_0_xyyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xyyyz_0[i] = g_yyyzz_0_xyyyz_0[i] * fbe_0 - g_yyyzz_0_xyyyz_1[i] * fz_be_0 + g_xyyyzz_0_yyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyyzz_0[i] = g_yyyzz_0_xyyzz_0[i] * fbe_0 - g_yyyzz_0_xyyzz_1[i] * fz_be_0 + g_xyyyzz_0_yyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyzzz_0[i] = g_yyyzz_0_xyzzz_0[i] * fbe_0 - g_yyyzz_0_xyzzz_1[i] * fz_be_0 + g_xyyyzz_0_yzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xzzzz_0[i] = 2.0 * g_xxyzz_0_xzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xzzzz_1[i] * fz_be_0 + g_xxyyzz_0_xzzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_yyyyy_0[i] = g_yyyzz_0_yyyyy_0[i] * fbe_0 - g_yyyzz_0_yyyyy_1[i] * fz_be_0 + g_xyyyzz_0_yyyyy_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyyz_0[i] = g_yyyzz_0_yyyyz_0[i] * fbe_0 - g_yyyzz_0_yyyyz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyzz_0[i] = g_yyyzz_0_yyyzz_0[i] * fbe_0 - g_yyyzz_0_yyyzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyzzz_0[i] = g_yyyzz_0_yyzzz_0[i] * fbe_0 - g_yyyzz_0_yyzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yzzzz_0[i] = g_yyyzz_0_yzzzz_0[i] * fbe_0 - g_yyyzz_0_yzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_zzzzz_0[i] = g_yyyzz_0_zzzzz_0[i] * fbe_0 - g_yyyzz_0_zzzzz_1[i] * fz_be_0 + g_xyyyzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 378-399 components of targeted buffer : KSH

    auto g_xxyyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 378);

    auto g_xxyyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 379);

    auto g_xxyyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 380);

    auto g_xxyyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 381);

    auto g_xxyyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 382);

    auto g_xxyyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 383);

    auto g_xxyyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 384);

    auto g_xxyyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 385);

    auto g_xxyyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 386);

    auto g_xxyyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 387);

    auto g_xxyyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 388);

    auto g_xxyyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 389);

    auto g_xxyyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 390);

    auto g_xxyyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 391);

    auto g_xxyyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 392);

    auto g_xxyyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 393);

    auto g_xxyyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 394);

    auto g_xxyyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 395);

    auto g_xxyyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 396);

    auto g_xxyyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 397);

    auto g_xxyyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 398);

    #pragma omp simd aligned(g_xxyyz_0_xxxxy_0, g_xxyyz_0_xxxxy_1, g_xxyyz_0_xxxyy_0, g_xxyyz_0_xxxyy_1, g_xxyyz_0_xxyyy_0, g_xxyyz_0_xxyyy_1, g_xxyyz_0_xyyyy_0, g_xxyyz_0_xyyyy_1, g_xxyyzz_0_xxxxy_1, g_xxyyzz_0_xxxyy_1, g_xxyyzz_0_xxyyy_1, g_xxyyzz_0_xyyyy_1, g_xxyyzzz_0_xxxxx_0, g_xxyyzzz_0_xxxxy_0, g_xxyyzzz_0_xxxxz_0, g_xxyyzzz_0_xxxyy_0, g_xxyyzzz_0_xxxyz_0, g_xxyyzzz_0_xxxzz_0, g_xxyyzzz_0_xxyyy_0, g_xxyyzzz_0_xxyyz_0, g_xxyyzzz_0_xxyzz_0, g_xxyyzzz_0_xxzzz_0, g_xxyyzzz_0_xyyyy_0, g_xxyyzzz_0_xyyyz_0, g_xxyyzzz_0_xyyzz_0, g_xxyyzzz_0_xyzzz_0, g_xxyyzzz_0_xzzzz_0, g_xxyyzzz_0_yyyyy_0, g_xxyyzzz_0_yyyyz_0, g_xxyyzzz_0_yyyzz_0, g_xxyyzzz_0_yyzzz_0, g_xxyyzzz_0_yzzzz_0, g_xxyyzzz_0_zzzzz_0, g_xxyzzz_0_xxxxx_1, g_xxyzzz_0_xxxxz_1, g_xxyzzz_0_xxxzz_1, g_xxyzzz_0_xxzzz_1, g_xxyzzz_0_xzzzz_1, g_xxzzz_0_xxxxx_0, g_xxzzz_0_xxxxx_1, g_xxzzz_0_xxxxz_0, g_xxzzz_0_xxxxz_1, g_xxzzz_0_xxxzz_0, g_xxzzz_0_xxxzz_1, g_xxzzz_0_xxzzz_0, g_xxzzz_0_xxzzz_1, g_xxzzz_0_xzzzz_0, g_xxzzz_0_xzzzz_1, g_xyyzzz_0_xxxyz_1, g_xyyzzz_0_xxyyz_1, g_xyyzzz_0_xxyz_1, g_xyyzzz_0_xxyzz_1, g_xyyzzz_0_xyyyz_1, g_xyyzzz_0_xyyz_1, g_xyyzzz_0_xyyzz_1, g_xyyzzz_0_xyzz_1, g_xyyzzz_0_xyzzz_1, g_xyyzzz_0_yyyyy_1, g_xyyzzz_0_yyyyz_1, g_xyyzzz_0_yyyz_1, g_xyyzzz_0_yyyzz_1, g_xyyzzz_0_yyzz_1, g_xyyzzz_0_yyzzz_1, g_xyyzzz_0_yzzz_1, g_xyyzzz_0_yzzzz_1, g_xyyzzz_0_zzzzz_1, g_yyzzz_0_xxxyz_0, g_yyzzz_0_xxxyz_1, g_yyzzz_0_xxyyz_0, g_yyzzz_0_xxyyz_1, g_yyzzz_0_xxyzz_0, g_yyzzz_0_xxyzz_1, g_yyzzz_0_xyyyz_0, g_yyzzz_0_xyyyz_1, g_yyzzz_0_xyyzz_0, g_yyzzz_0_xyyzz_1, g_yyzzz_0_xyzzz_0, g_yyzzz_0_xyzzz_1, g_yyzzz_0_yyyyy_0, g_yyzzz_0_yyyyy_1, g_yyzzz_0_yyyyz_0, g_yyzzz_0_yyyyz_1, g_yyzzz_0_yyyzz_0, g_yyzzz_0_yyyzz_1, g_yyzzz_0_yyzzz_0, g_yyzzz_0_yyzzz_1, g_yyzzz_0_yzzzz_0, g_yyzzz_0_yzzzz_1, g_yyzzz_0_zzzzz_0, g_yyzzz_0_zzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyzzz_0_xxxxx_0[i] = g_xxzzz_0_xxxxx_0[i] * fbe_0 - g_xxzzz_0_xxxxx_1[i] * fz_be_0 + g_xxyzzz_0_xxxxx_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxxy_0[i] = 2.0 * g_xxyyz_0_xxxxy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxxy_1[i] * fz_be_0 + g_xxyyzz_0_xxxxy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxxz_0[i] = g_xxzzz_0_xxxxz_0[i] * fbe_0 - g_xxzzz_0_xxxxz_1[i] * fz_be_0 + g_xxyzzz_0_xxxxz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxyy_0[i] = 2.0 * g_xxyyz_0_xxxyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxyy_1[i] * fz_be_0 + g_xxyyzz_0_xxxyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxyz_0[i] = g_yyzzz_0_xxxyz_0[i] * fbe_0 - g_yyzzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_0_xxyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxzz_0[i] = g_xxzzz_0_xxxzz_0[i] * fbe_0 - g_xxzzz_0_xxxzz_1[i] * fz_be_0 + g_xxyzzz_0_xxxzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxyyy_0[i] = 2.0 * g_xxyyz_0_xxyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxyyy_1[i] * fz_be_0 + g_xxyyzz_0_xxyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxyyz_0[i] = g_yyzzz_0_xxyyz_0[i] * fbe_0 - g_yyzzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxyzz_0[i] = g_yyzzz_0_xxyzz_0[i] * fbe_0 - g_yyzzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxzzz_0[i] = g_xxzzz_0_xxzzz_0[i] * fbe_0 - g_xxzzz_0_xxzzz_1[i] * fz_be_0 + g_xxyzzz_0_xxzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xyyyy_0[i] = 2.0 * g_xxyyz_0_xyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xyyyy_1[i] * fz_be_0 + g_xxyyzz_0_xyyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xyyyz_0[i] = g_yyzzz_0_xyyyz_0[i] * fbe_0 - g_yyzzz_0_xyyyz_1[i] * fz_be_0 + g_xyyzzz_0_yyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyyzz_0[i] = g_yyzzz_0_xyyzz_0[i] * fbe_0 - g_yyzzz_0_xyyzz_1[i] * fz_be_0 + g_xyyzzz_0_yyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyzzz_0[i] = g_yyzzz_0_xyzzz_0[i] * fbe_0 - g_yyzzz_0_xyzzz_1[i] * fz_be_0 + g_xyyzzz_0_yzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xzzzz_0[i] = g_xxzzz_0_xzzzz_0[i] * fbe_0 - g_xxzzz_0_xzzzz_1[i] * fz_be_0 + g_xxyzzz_0_xzzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_yyyyy_0[i] = g_yyzzz_0_yyyyy_0[i] * fbe_0 - g_yyzzz_0_yyyyy_1[i] * fz_be_0 + g_xyyzzz_0_yyyyy_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyyz_0[i] = g_yyzzz_0_yyyyz_0[i] * fbe_0 - g_yyzzz_0_yyyyz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyzz_0[i] = g_yyzzz_0_yyyzz_0[i] * fbe_0 - g_yyzzz_0_yyyzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyzzz_0[i] = g_yyzzz_0_yyzzz_0[i] * fbe_0 - g_yyzzz_0_yyzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yzzzz_0[i] = g_yyzzz_0_yzzzz_0[i] * fbe_0 - g_yyzzz_0_yzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_zzzzz_0[i] = g_yyzzz_0_zzzzz_0[i] * fbe_0 - g_yyzzz_0_zzzzz_1[i] * fz_be_0 + g_xyyzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 399-420 components of targeted buffer : KSH

    auto g_xxyzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 399);

    auto g_xxyzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 400);

    auto g_xxyzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 401);

    auto g_xxyzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 402);

    auto g_xxyzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 403);

    auto g_xxyzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 404);

    auto g_xxyzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 405);

    auto g_xxyzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 406);

    auto g_xxyzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 407);

    auto g_xxyzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 408);

    auto g_xxyzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 409);

    auto g_xxyzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 410);

    auto g_xxyzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 411);

    auto g_xxyzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 412);

    auto g_xxyzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 413);

    auto g_xxyzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 414);

    auto g_xxyzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 415);

    auto g_xxyzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 416);

    auto g_xxyzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 417);

    auto g_xxyzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 418);

    auto g_xxyzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 419);

    #pragma omp simd aligned(g_xxyzzzz_0_xxxxx_0, g_xxyzzzz_0_xxxxy_0, g_xxyzzzz_0_xxxxz_0, g_xxyzzzz_0_xxxyy_0, g_xxyzzzz_0_xxxyz_0, g_xxyzzzz_0_xxxzz_0, g_xxyzzzz_0_xxyyy_0, g_xxyzzzz_0_xxyyz_0, g_xxyzzzz_0_xxyzz_0, g_xxyzzzz_0_xxzzz_0, g_xxyzzzz_0_xyyyy_0, g_xxyzzzz_0_xyyyz_0, g_xxyzzzz_0_xyyzz_0, g_xxyzzzz_0_xyzzz_0, g_xxyzzzz_0_xzzzz_0, g_xxyzzzz_0_yyyyy_0, g_xxyzzzz_0_yyyyz_0, g_xxyzzzz_0_yyyzz_0, g_xxyzzzz_0_yyzzz_0, g_xxyzzzz_0_yzzzz_0, g_xxyzzzz_0_zzzzz_0, g_xxzzzz_0_xxxx_1, g_xxzzzz_0_xxxxx_1, g_xxzzzz_0_xxxxy_1, g_xxzzzz_0_xxxxz_1, g_xxzzzz_0_xxxy_1, g_xxzzzz_0_xxxyy_1, g_xxzzzz_0_xxxyz_1, g_xxzzzz_0_xxxz_1, g_xxzzzz_0_xxxzz_1, g_xxzzzz_0_xxyy_1, g_xxzzzz_0_xxyyy_1, g_xxzzzz_0_xxyyz_1, g_xxzzzz_0_xxyz_1, g_xxzzzz_0_xxyzz_1, g_xxzzzz_0_xxzz_1, g_xxzzzz_0_xxzzz_1, g_xxzzzz_0_xyyy_1, g_xxzzzz_0_xyyyy_1, g_xxzzzz_0_xyyyz_1, g_xxzzzz_0_xyyz_1, g_xxzzzz_0_xyyzz_1, g_xxzzzz_0_xyzz_1, g_xxzzzz_0_xyzzz_1, g_xxzzzz_0_xzzz_1, g_xxzzzz_0_xzzzz_1, g_xxzzzz_0_yyyy_1, g_xxzzzz_0_yyyyy_1, g_xxzzzz_0_yyyyz_1, g_xxzzzz_0_yyyz_1, g_xxzzzz_0_yyyzz_1, g_xxzzzz_0_yyzz_1, g_xxzzzz_0_yyzzz_1, g_xxzzzz_0_yzzz_1, g_xxzzzz_0_yzzzz_1, g_xxzzzz_0_zzzz_1, g_xxzzzz_0_zzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzzz_0_xxxxx_0[i] = g_xxzzzz_0_xxxxx_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxy_0[i] = g_xxzzzz_0_xxxx_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxz_0[i] = g_xxzzzz_0_xxxxz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyy_0[i] = 2.0 * g_xxzzzz_0_xxxy_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyz_0[i] = g_xxzzzz_0_xxxz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxzz_0[i] = g_xxzzzz_0_xxxzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyy_0[i] = 3.0 * g_xxzzzz_0_xxyy_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyz_0[i] = 2.0 * g_xxzzzz_0_xxyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyzz_0[i] = g_xxzzzz_0_xxzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxzzz_0[i] = g_xxzzzz_0_xxzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyy_0[i] = 4.0 * g_xxzzzz_0_xyyy_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyz_0[i] = 3.0 * g_xxzzzz_0_xyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyzz_0[i] = 2.0 * g_xxzzzz_0_xyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyzzz_0[i] = g_xxzzzz_0_xzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xzzzz_0[i] = g_xxzzzz_0_xzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyy_0[i] = 5.0 * g_xxzzzz_0_yyyy_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyz_0[i] = 4.0 * g_xxzzzz_0_yyyz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyzz_0[i] = 3.0 * g_xxzzzz_0_yyzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyzzz_0[i] = 2.0 * g_xxzzzz_0_yzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yzzzz_0[i] = g_xxzzzz_0_zzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_zzzzz_0[i] = g_xxzzzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 420-441 components of targeted buffer : KSH

    auto g_xxzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 420);

    auto g_xxzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 421);

    auto g_xxzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 422);

    auto g_xxzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 423);

    auto g_xxzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 424);

    auto g_xxzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 425);

    auto g_xxzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 426);

    auto g_xxzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 427);

    auto g_xxzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 428);

    auto g_xxzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 429);

    auto g_xxzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 430);

    auto g_xxzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 431);

    auto g_xxzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 432);

    auto g_xxzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 433);

    auto g_xxzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 434);

    auto g_xxzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 435);

    auto g_xxzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 436);

    auto g_xxzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 437);

    auto g_xxzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 438);

    auto g_xxzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 439);

    auto g_xxzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 440);

    #pragma omp simd aligned(g_xxzzz_0_xxxxx_0, g_xxzzz_0_xxxxx_1, g_xxzzz_0_xxxxy_0, g_xxzzz_0_xxxxy_1, g_xxzzz_0_xxxyy_0, g_xxzzz_0_xxxyy_1, g_xxzzz_0_xxyyy_0, g_xxzzz_0_xxyyy_1, g_xxzzz_0_xyyyy_0, g_xxzzz_0_xyyyy_1, g_xxzzzz_0_xxxxx_1, g_xxzzzz_0_xxxxy_1, g_xxzzzz_0_xxxyy_1, g_xxzzzz_0_xxyyy_1, g_xxzzzz_0_xyyyy_1, g_xxzzzzz_0_xxxxx_0, g_xxzzzzz_0_xxxxy_0, g_xxzzzzz_0_xxxxz_0, g_xxzzzzz_0_xxxyy_0, g_xxzzzzz_0_xxxyz_0, g_xxzzzzz_0_xxxzz_0, g_xxzzzzz_0_xxyyy_0, g_xxzzzzz_0_xxyyz_0, g_xxzzzzz_0_xxyzz_0, g_xxzzzzz_0_xxzzz_0, g_xxzzzzz_0_xyyyy_0, g_xxzzzzz_0_xyyyz_0, g_xxzzzzz_0_xyyzz_0, g_xxzzzzz_0_xyzzz_0, g_xxzzzzz_0_xzzzz_0, g_xxzzzzz_0_yyyyy_0, g_xxzzzzz_0_yyyyz_0, g_xxzzzzz_0_yyyzz_0, g_xxzzzzz_0_yyzzz_0, g_xxzzzzz_0_yzzzz_0, g_xxzzzzz_0_zzzzz_0, g_xzzzzz_0_xxxxz_1, g_xzzzzz_0_xxxyz_1, g_xzzzzz_0_xxxz_1, g_xzzzzz_0_xxxzz_1, g_xzzzzz_0_xxyyz_1, g_xzzzzz_0_xxyz_1, g_xzzzzz_0_xxyzz_1, g_xzzzzz_0_xxzz_1, g_xzzzzz_0_xxzzz_1, g_xzzzzz_0_xyyyz_1, g_xzzzzz_0_xyyz_1, g_xzzzzz_0_xyyzz_1, g_xzzzzz_0_xyzz_1, g_xzzzzz_0_xyzzz_1, g_xzzzzz_0_xzzz_1, g_xzzzzz_0_xzzzz_1, g_xzzzzz_0_yyyyy_1, g_xzzzzz_0_yyyyz_1, g_xzzzzz_0_yyyz_1, g_xzzzzz_0_yyyzz_1, g_xzzzzz_0_yyzz_1, g_xzzzzz_0_yyzzz_1, g_xzzzzz_0_yzzz_1, g_xzzzzz_0_yzzzz_1, g_xzzzzz_0_zzzz_1, g_xzzzzz_0_zzzzz_1, g_zzzzz_0_xxxxz_0, g_zzzzz_0_xxxxz_1, g_zzzzz_0_xxxyz_0, g_zzzzz_0_xxxyz_1, g_zzzzz_0_xxxzz_0, g_zzzzz_0_xxxzz_1, g_zzzzz_0_xxyyz_0, g_zzzzz_0_xxyyz_1, g_zzzzz_0_xxyzz_0, g_zzzzz_0_xxyzz_1, g_zzzzz_0_xxzzz_0, g_zzzzz_0_xxzzz_1, g_zzzzz_0_xyyyz_0, g_zzzzz_0_xyyyz_1, g_zzzzz_0_xyyzz_0, g_zzzzz_0_xyyzz_1, g_zzzzz_0_xyzzz_0, g_zzzzz_0_xyzzz_1, g_zzzzz_0_xzzzz_0, g_zzzzz_0_xzzzz_1, g_zzzzz_0_yyyyy_0, g_zzzzz_0_yyyyy_1, g_zzzzz_0_yyyyz_0, g_zzzzz_0_yyyyz_1, g_zzzzz_0_yyyzz_0, g_zzzzz_0_yyyzz_1, g_zzzzz_0_yyzzz_0, g_zzzzz_0_yyzzz_1, g_zzzzz_0_yzzzz_0, g_zzzzz_0_yzzzz_1, g_zzzzz_0_zzzzz_0, g_zzzzz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzzz_0_xxxxx_0[i] = 4.0 * g_xxzzz_0_xxxxx_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxx_1[i] * fz_be_0 + g_xxzzzz_0_xxxxx_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxy_0[i] = 4.0 * g_xxzzz_0_xxxxy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxy_1[i] * fz_be_0 + g_xxzzzz_0_xxxxy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxz_0[i] = g_zzzzz_0_xxxxz_0[i] * fbe_0 - g_zzzzz_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_0_xxxz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxyy_0[i] = 4.0 * g_xxzzz_0_xxxyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxyy_1[i] * fz_be_0 + g_xxzzzz_0_xxxyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxyz_0[i] = g_zzzzz_0_xxxyz_0[i] * fbe_0 - g_zzzzz_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxzz_0[i] = g_zzzzz_0_xxxzz_0[i] * fbe_0 - g_zzzzz_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyyy_0[i] = 4.0 * g_xxzzz_0_xxyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxyyy_1[i] * fz_be_0 + g_xxzzzz_0_xxyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxyyz_0[i] = g_zzzzz_0_xxyyz_0[i] * fbe_0 - g_zzzzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyzz_0[i] = g_zzzzz_0_xxyzz_0[i] * fbe_0 - g_zzzzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxzzz_0[i] = g_zzzzz_0_xxzzz_0[i] * fbe_0 - g_zzzzz_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyyy_0[i] = 4.0 * g_xxzzz_0_xyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xyyyy_1[i] * fz_be_0 + g_xxzzzz_0_xyyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xyyyz_0[i] = g_zzzzz_0_xyyyz_0[i] * fbe_0 - g_zzzzz_0_xyyyz_1[i] * fz_be_0 + g_xzzzzz_0_yyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyzz_0[i] = g_zzzzz_0_xyyzz_0[i] * fbe_0 - g_zzzzz_0_xyyzz_1[i] * fz_be_0 + g_xzzzzz_0_yyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyzzz_0[i] = g_zzzzz_0_xyzzz_0[i] * fbe_0 - g_zzzzz_0_xyzzz_1[i] * fz_be_0 + g_xzzzzz_0_yzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xzzzz_0[i] = g_zzzzz_0_xzzzz_0[i] * fbe_0 - g_zzzzz_0_xzzzz_1[i] * fz_be_0 + g_xzzzzz_0_zzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyy_0[i] = g_zzzzz_0_yyyyy_0[i] * fbe_0 - g_zzzzz_0_yyyyy_1[i] * fz_be_0 + g_xzzzzz_0_yyyyy_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyz_0[i] = g_zzzzz_0_yyyyz_0[i] * fbe_0 - g_zzzzz_0_yyyyz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyzz_0[i] = g_zzzzz_0_yyyzz_0[i] * fbe_0 - g_zzzzz_0_yyyzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyzzz_0[i] = g_zzzzz_0_yyzzz_0[i] * fbe_0 - g_zzzzz_0_yyzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yzzzz_0[i] = g_zzzzz_0_yzzzz_0[i] * fbe_0 - g_zzzzz_0_yzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_zzzzz_0[i] = g_zzzzz_0_zzzzz_0[i] * fbe_0 - g_zzzzz_0_zzzzz_1[i] * fz_be_0 + g_xzzzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 441-462 components of targeted buffer : KSH

    auto g_xyyyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 441);

    auto g_xyyyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 442);

    auto g_xyyyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 443);

    auto g_xyyyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 444);

    auto g_xyyyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 445);

    auto g_xyyyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 446);

    auto g_xyyyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 447);

    auto g_xyyyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 448);

    auto g_xyyyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 449);

    auto g_xyyyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 450);

    auto g_xyyyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 451);

    auto g_xyyyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 452);

    auto g_xyyyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 453);

    auto g_xyyyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 454);

    auto g_xyyyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 455);

    auto g_xyyyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 456);

    auto g_xyyyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 457);

    auto g_xyyyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 458);

    auto g_xyyyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 459);

    auto g_xyyyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 460);

    auto g_xyyyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 461);

    #pragma omp simd aligned(g_xyyyyyy_0_xxxxx_0, g_xyyyyyy_0_xxxxy_0, g_xyyyyyy_0_xxxxz_0, g_xyyyyyy_0_xxxyy_0, g_xyyyyyy_0_xxxyz_0, g_xyyyyyy_0_xxxzz_0, g_xyyyyyy_0_xxyyy_0, g_xyyyyyy_0_xxyyz_0, g_xyyyyyy_0_xxyzz_0, g_xyyyyyy_0_xxzzz_0, g_xyyyyyy_0_xyyyy_0, g_xyyyyyy_0_xyyyz_0, g_xyyyyyy_0_xyyzz_0, g_xyyyyyy_0_xyzzz_0, g_xyyyyyy_0_xzzzz_0, g_xyyyyyy_0_yyyyy_0, g_xyyyyyy_0_yyyyz_0, g_xyyyyyy_0_yyyzz_0, g_xyyyyyy_0_yyzzz_0, g_xyyyyyy_0_yzzzz_0, g_xyyyyyy_0_zzzzz_0, g_yyyyyy_0_xxxx_1, g_yyyyyy_0_xxxxx_1, g_yyyyyy_0_xxxxy_1, g_yyyyyy_0_xxxxz_1, g_yyyyyy_0_xxxy_1, g_yyyyyy_0_xxxyy_1, g_yyyyyy_0_xxxyz_1, g_yyyyyy_0_xxxz_1, g_yyyyyy_0_xxxzz_1, g_yyyyyy_0_xxyy_1, g_yyyyyy_0_xxyyy_1, g_yyyyyy_0_xxyyz_1, g_yyyyyy_0_xxyz_1, g_yyyyyy_0_xxyzz_1, g_yyyyyy_0_xxzz_1, g_yyyyyy_0_xxzzz_1, g_yyyyyy_0_xyyy_1, g_yyyyyy_0_xyyyy_1, g_yyyyyy_0_xyyyz_1, g_yyyyyy_0_xyyz_1, g_yyyyyy_0_xyyzz_1, g_yyyyyy_0_xyzz_1, g_yyyyyy_0_xyzzz_1, g_yyyyyy_0_xzzz_1, g_yyyyyy_0_xzzzz_1, g_yyyyyy_0_yyyy_1, g_yyyyyy_0_yyyyy_1, g_yyyyyy_0_yyyyz_1, g_yyyyyy_0_yyyz_1, g_yyyyyy_0_yyyzz_1, g_yyyyyy_0_yyzz_1, g_yyyyyy_0_yyzzz_1, g_yyyyyy_0_yzzz_1, g_yyyyyy_0_yzzzz_1, g_yyyyyy_0_zzzz_1, g_yyyyyy_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyy_0_xxxxx_0[i] = 5.0 * g_yyyyyy_0_xxxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxx_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxy_0[i] = 4.0 * g_yyyyyy_0_xxxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxz_0[i] = 4.0 * g_yyyyyy_0_xxxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyy_0[i] = 3.0 * g_yyyyyy_0_xxyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyz_0[i] = 3.0 * g_yyyyyy_0_xxyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxzz_0[i] = 3.0 * g_yyyyyy_0_xxzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyy_0[i] = 2.0 * g_yyyyyy_0_xyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyz_0[i] = 2.0 * g_yyyyyy_0_xyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyzz_0[i] = 2.0 * g_yyyyyy_0_xyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxzzz_0[i] = 2.0 * g_yyyyyy_0_xzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyy_0[i] = g_yyyyyy_0_yyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyz_0[i] = g_yyyyyy_0_yyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyzz_0[i] = g_yyyyyy_0_yyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyzzz_0[i] = g_yyyyyy_0_yzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xzzzz_0[i] = g_yyyyyy_0_zzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyy_0[i] = g_yyyyyy_0_yyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyz_0[i] = g_yyyyyy_0_yyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyzz_0[i] = g_yyyyyy_0_yyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyzzz_0[i] = g_yyyyyy_0_yyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yzzzz_0[i] = g_yyyyyy_0_yzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_zzzzz_0[i] = g_yyyyyy_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 462-483 components of targeted buffer : KSH

    auto g_xyyyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 462);

    auto g_xyyyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 463);

    auto g_xyyyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 464);

    auto g_xyyyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 465);

    auto g_xyyyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 466);

    auto g_xyyyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 467);

    auto g_xyyyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 468);

    auto g_xyyyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 469);

    auto g_xyyyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 470);

    auto g_xyyyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 471);

    auto g_xyyyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 472);

    auto g_xyyyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 473);

    auto g_xyyyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 474);

    auto g_xyyyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 475);

    auto g_xyyyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 476);

    auto g_xyyyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 477);

    auto g_xyyyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 478);

    auto g_xyyyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 479);

    auto g_xyyyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 480);

    auto g_xyyyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 481);

    auto g_xyyyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 482);

    #pragma omp simd aligned(g_xyyyyy_0_xxxxx_1, g_xyyyyy_0_xxxxy_1, g_xyyyyy_0_xxxyy_1, g_xyyyyy_0_xxyyy_1, g_xyyyyy_0_xyyyy_1, g_xyyyyyz_0_xxxxx_0, g_xyyyyyz_0_xxxxy_0, g_xyyyyyz_0_xxxxz_0, g_xyyyyyz_0_xxxyy_0, g_xyyyyyz_0_xxxyz_0, g_xyyyyyz_0_xxxzz_0, g_xyyyyyz_0_xxyyy_0, g_xyyyyyz_0_xxyyz_0, g_xyyyyyz_0_xxyzz_0, g_xyyyyyz_0_xxzzz_0, g_xyyyyyz_0_xyyyy_0, g_xyyyyyz_0_xyyyz_0, g_xyyyyyz_0_xyyzz_0, g_xyyyyyz_0_xyzzz_0, g_xyyyyyz_0_xzzzz_0, g_xyyyyyz_0_yyyyy_0, g_xyyyyyz_0_yyyyz_0, g_xyyyyyz_0_yyyzz_0, g_xyyyyyz_0_yyzzz_0, g_xyyyyyz_0_yzzzz_0, g_xyyyyyz_0_zzzzz_0, g_yyyyyz_0_xxxxz_1, g_yyyyyz_0_xxxyz_1, g_yyyyyz_0_xxxz_1, g_yyyyyz_0_xxxzz_1, g_yyyyyz_0_xxyyz_1, g_yyyyyz_0_xxyz_1, g_yyyyyz_0_xxyzz_1, g_yyyyyz_0_xxzz_1, g_yyyyyz_0_xxzzz_1, g_yyyyyz_0_xyyyz_1, g_yyyyyz_0_xyyz_1, g_yyyyyz_0_xyyzz_1, g_yyyyyz_0_xyzz_1, g_yyyyyz_0_xyzzz_1, g_yyyyyz_0_xzzz_1, g_yyyyyz_0_xzzzz_1, g_yyyyyz_0_yyyyy_1, g_yyyyyz_0_yyyyz_1, g_yyyyyz_0_yyyz_1, g_yyyyyz_0_yyyzz_1, g_yyyyyz_0_yyzz_1, g_yyyyyz_0_yyzzz_1, g_yyyyyz_0_yzzz_1, g_yyyyyz_0_yzzzz_1, g_yyyyyz_0_zzzz_1, g_yyyyyz_0_zzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyz_0_xxxxx_0[i] = g_xyyyyy_0_xxxxx_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxy_0[i] = g_xyyyyy_0_xxxxy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxz_0[i] = 4.0 * g_yyyyyz_0_xxxz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxyy_0[i] = g_xyyyyy_0_xxxyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxyz_0[i] = 3.0 * g_yyyyyz_0_xxyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxzz_0[i] = 3.0 * g_yyyyyz_0_xxzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyyy_0[i] = g_xyyyyy_0_xxyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxyyz_0[i] = 2.0 * g_yyyyyz_0_xyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyzz_0[i] = 2.0 * g_yyyyyz_0_xyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxzzz_0[i] = 2.0 * g_yyyyyz_0_xzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyyy_0[i] = g_xyyyyy_0_xyyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xyyyz_0[i] = g_yyyyyz_0_yyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyzz_0[i] = g_yyyyyz_0_yyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyzzz_0[i] = g_yyyyyz_0_yzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xzzzz_0[i] = g_yyyyyz_0_zzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyy_0[i] = g_yyyyyz_0_yyyyy_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyz_0[i] = g_yyyyyz_0_yyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyzz_0[i] = g_yyyyyz_0_yyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyzzz_0[i] = g_yyyyyz_0_yyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yzzzz_0[i] = g_yyyyyz_0_yzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_zzzzz_0[i] = g_yyyyyz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 483-504 components of targeted buffer : KSH

    auto g_xyyyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 483);

    auto g_xyyyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 484);

    auto g_xyyyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 485);

    auto g_xyyyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 486);

    auto g_xyyyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 487);

    auto g_xyyyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 488);

    auto g_xyyyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 489);

    auto g_xyyyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 490);

    auto g_xyyyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 491);

    auto g_xyyyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 492);

    auto g_xyyyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 493);

    auto g_xyyyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 494);

    auto g_xyyyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 495);

    auto g_xyyyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 496);

    auto g_xyyyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 497);

    auto g_xyyyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 498);

    auto g_xyyyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 499);

    auto g_xyyyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 500);

    auto g_xyyyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 501);

    auto g_xyyyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 502);

    auto g_xyyyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 503);

    #pragma omp simd aligned(g_xyyyyzz_0_xxxxx_0, g_xyyyyzz_0_xxxxy_0, g_xyyyyzz_0_xxxxz_0, g_xyyyyzz_0_xxxyy_0, g_xyyyyzz_0_xxxyz_0, g_xyyyyzz_0_xxxzz_0, g_xyyyyzz_0_xxyyy_0, g_xyyyyzz_0_xxyyz_0, g_xyyyyzz_0_xxyzz_0, g_xyyyyzz_0_xxzzz_0, g_xyyyyzz_0_xyyyy_0, g_xyyyyzz_0_xyyyz_0, g_xyyyyzz_0_xyyzz_0, g_xyyyyzz_0_xyzzz_0, g_xyyyyzz_0_xzzzz_0, g_xyyyyzz_0_yyyyy_0, g_xyyyyzz_0_yyyyz_0, g_xyyyyzz_0_yyyzz_0, g_xyyyyzz_0_yyzzz_0, g_xyyyyzz_0_yzzzz_0, g_xyyyyzz_0_zzzzz_0, g_yyyyzz_0_xxxx_1, g_yyyyzz_0_xxxxx_1, g_yyyyzz_0_xxxxy_1, g_yyyyzz_0_xxxxz_1, g_yyyyzz_0_xxxy_1, g_yyyyzz_0_xxxyy_1, g_yyyyzz_0_xxxyz_1, g_yyyyzz_0_xxxz_1, g_yyyyzz_0_xxxzz_1, g_yyyyzz_0_xxyy_1, g_yyyyzz_0_xxyyy_1, g_yyyyzz_0_xxyyz_1, g_yyyyzz_0_xxyz_1, g_yyyyzz_0_xxyzz_1, g_yyyyzz_0_xxzz_1, g_yyyyzz_0_xxzzz_1, g_yyyyzz_0_xyyy_1, g_yyyyzz_0_xyyyy_1, g_yyyyzz_0_xyyyz_1, g_yyyyzz_0_xyyz_1, g_yyyyzz_0_xyyzz_1, g_yyyyzz_0_xyzz_1, g_yyyyzz_0_xyzzz_1, g_yyyyzz_0_xzzz_1, g_yyyyzz_0_xzzzz_1, g_yyyyzz_0_yyyy_1, g_yyyyzz_0_yyyyy_1, g_yyyyzz_0_yyyyz_1, g_yyyyzz_0_yyyz_1, g_yyyyzz_0_yyyzz_1, g_yyyyzz_0_yyzz_1, g_yyyyzz_0_yyzzz_1, g_yyyyzz_0_yzzz_1, g_yyyyzz_0_yzzzz_1, g_yyyyzz_0_zzzz_1, g_yyyyzz_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyzz_0_xxxxx_0[i] = 5.0 * g_yyyyzz_0_xxxx_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxx_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxy_0[i] = 4.0 * g_yyyyzz_0_xxxy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxz_0[i] = 4.0 * g_yyyyzz_0_xxxz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyy_0[i] = 3.0 * g_yyyyzz_0_xxyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyz_0[i] = 3.0 * g_yyyyzz_0_xxyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxzz_0[i] = 3.0 * g_yyyyzz_0_xxzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyy_0[i] = 2.0 * g_yyyyzz_0_xyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyz_0[i] = 2.0 * g_yyyyzz_0_xyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyzz_0[i] = 2.0 * g_yyyyzz_0_xyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxzzz_0[i] = 2.0 * g_yyyyzz_0_xzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyy_0[i] = g_yyyyzz_0_yyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyz_0[i] = g_yyyyzz_0_yyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyzz_0[i] = g_yyyyzz_0_yyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyzzz_0[i] = g_yyyyzz_0_yzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xzzzz_0[i] = g_yyyyzz_0_zzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyy_0[i] = g_yyyyzz_0_yyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyz_0[i] = g_yyyyzz_0_yyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyzz_0[i] = g_yyyyzz_0_yyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyzzz_0[i] = g_yyyyzz_0_yyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yzzzz_0[i] = g_yyyyzz_0_yzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_zzzzz_0[i] = g_yyyyzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 504-525 components of targeted buffer : KSH

    auto g_xyyyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 504);

    auto g_xyyyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 505);

    auto g_xyyyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 506);

    auto g_xyyyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 507);

    auto g_xyyyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 508);

    auto g_xyyyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 509);

    auto g_xyyyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 510);

    auto g_xyyyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 511);

    auto g_xyyyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 512);

    auto g_xyyyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 513);

    auto g_xyyyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 514);

    auto g_xyyyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 515);

    auto g_xyyyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 516);

    auto g_xyyyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 517);

    auto g_xyyyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 518);

    auto g_xyyyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 519);

    auto g_xyyyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 520);

    auto g_xyyyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 521);

    auto g_xyyyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 522);

    auto g_xyyyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 523);

    auto g_xyyyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 524);

    #pragma omp simd aligned(g_xyyyzzz_0_xxxxx_0, g_xyyyzzz_0_xxxxy_0, g_xyyyzzz_0_xxxxz_0, g_xyyyzzz_0_xxxyy_0, g_xyyyzzz_0_xxxyz_0, g_xyyyzzz_0_xxxzz_0, g_xyyyzzz_0_xxyyy_0, g_xyyyzzz_0_xxyyz_0, g_xyyyzzz_0_xxyzz_0, g_xyyyzzz_0_xxzzz_0, g_xyyyzzz_0_xyyyy_0, g_xyyyzzz_0_xyyyz_0, g_xyyyzzz_0_xyyzz_0, g_xyyyzzz_0_xyzzz_0, g_xyyyzzz_0_xzzzz_0, g_xyyyzzz_0_yyyyy_0, g_xyyyzzz_0_yyyyz_0, g_xyyyzzz_0_yyyzz_0, g_xyyyzzz_0_yyzzz_0, g_xyyyzzz_0_yzzzz_0, g_xyyyzzz_0_zzzzz_0, g_yyyzzz_0_xxxx_1, g_yyyzzz_0_xxxxx_1, g_yyyzzz_0_xxxxy_1, g_yyyzzz_0_xxxxz_1, g_yyyzzz_0_xxxy_1, g_yyyzzz_0_xxxyy_1, g_yyyzzz_0_xxxyz_1, g_yyyzzz_0_xxxz_1, g_yyyzzz_0_xxxzz_1, g_yyyzzz_0_xxyy_1, g_yyyzzz_0_xxyyy_1, g_yyyzzz_0_xxyyz_1, g_yyyzzz_0_xxyz_1, g_yyyzzz_0_xxyzz_1, g_yyyzzz_0_xxzz_1, g_yyyzzz_0_xxzzz_1, g_yyyzzz_0_xyyy_1, g_yyyzzz_0_xyyyy_1, g_yyyzzz_0_xyyyz_1, g_yyyzzz_0_xyyz_1, g_yyyzzz_0_xyyzz_1, g_yyyzzz_0_xyzz_1, g_yyyzzz_0_xyzzz_1, g_yyyzzz_0_xzzz_1, g_yyyzzz_0_xzzzz_1, g_yyyzzz_0_yyyy_1, g_yyyzzz_0_yyyyy_1, g_yyyzzz_0_yyyyz_1, g_yyyzzz_0_yyyz_1, g_yyyzzz_0_yyyzz_1, g_yyyzzz_0_yyzz_1, g_yyyzzz_0_yyzzz_1, g_yyyzzz_0_yzzz_1, g_yyyzzz_0_yzzzz_1, g_yyyzzz_0_zzzz_1, g_yyyzzz_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzzz_0_xxxxx_0[i] = 5.0 * g_yyyzzz_0_xxxx_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxx_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxy_0[i] = 4.0 * g_yyyzzz_0_xxxy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxz_0[i] = 4.0 * g_yyyzzz_0_xxxz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyy_0[i] = 3.0 * g_yyyzzz_0_xxyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyz_0[i] = 3.0 * g_yyyzzz_0_xxyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxzz_0[i] = 3.0 * g_yyyzzz_0_xxzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyy_0[i] = 2.0 * g_yyyzzz_0_xyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyz_0[i] = 2.0 * g_yyyzzz_0_xyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyzz_0[i] = 2.0 * g_yyyzzz_0_xyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxzzz_0[i] = 2.0 * g_yyyzzz_0_xzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyy_0[i] = g_yyyzzz_0_yyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyz_0[i] = g_yyyzzz_0_yyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyzz_0[i] = g_yyyzzz_0_yyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyzzz_0[i] = g_yyyzzz_0_yzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xzzzz_0[i] = g_yyyzzz_0_zzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyy_0[i] = g_yyyzzz_0_yyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyz_0[i] = g_yyyzzz_0_yyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyzz_0[i] = g_yyyzzz_0_yyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyzzz_0[i] = g_yyyzzz_0_yyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yzzzz_0[i] = g_yyyzzz_0_yzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_zzzzz_0[i] = g_yyyzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 525-546 components of targeted buffer : KSH

    auto g_xyyzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 525);

    auto g_xyyzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 526);

    auto g_xyyzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 527);

    auto g_xyyzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 528);

    auto g_xyyzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 529);

    auto g_xyyzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 530);

    auto g_xyyzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 531);

    auto g_xyyzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 532);

    auto g_xyyzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 533);

    auto g_xyyzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 534);

    auto g_xyyzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 535);

    auto g_xyyzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 536);

    auto g_xyyzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 537);

    auto g_xyyzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 538);

    auto g_xyyzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 539);

    auto g_xyyzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 540);

    auto g_xyyzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 541);

    auto g_xyyzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 542);

    auto g_xyyzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 543);

    auto g_xyyzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 544);

    auto g_xyyzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 545);

    #pragma omp simd aligned(g_xyyzzzz_0_xxxxx_0, g_xyyzzzz_0_xxxxy_0, g_xyyzzzz_0_xxxxz_0, g_xyyzzzz_0_xxxyy_0, g_xyyzzzz_0_xxxyz_0, g_xyyzzzz_0_xxxzz_0, g_xyyzzzz_0_xxyyy_0, g_xyyzzzz_0_xxyyz_0, g_xyyzzzz_0_xxyzz_0, g_xyyzzzz_0_xxzzz_0, g_xyyzzzz_0_xyyyy_0, g_xyyzzzz_0_xyyyz_0, g_xyyzzzz_0_xyyzz_0, g_xyyzzzz_0_xyzzz_0, g_xyyzzzz_0_xzzzz_0, g_xyyzzzz_0_yyyyy_0, g_xyyzzzz_0_yyyyz_0, g_xyyzzzz_0_yyyzz_0, g_xyyzzzz_0_yyzzz_0, g_xyyzzzz_0_yzzzz_0, g_xyyzzzz_0_zzzzz_0, g_yyzzzz_0_xxxx_1, g_yyzzzz_0_xxxxx_1, g_yyzzzz_0_xxxxy_1, g_yyzzzz_0_xxxxz_1, g_yyzzzz_0_xxxy_1, g_yyzzzz_0_xxxyy_1, g_yyzzzz_0_xxxyz_1, g_yyzzzz_0_xxxz_1, g_yyzzzz_0_xxxzz_1, g_yyzzzz_0_xxyy_1, g_yyzzzz_0_xxyyy_1, g_yyzzzz_0_xxyyz_1, g_yyzzzz_0_xxyz_1, g_yyzzzz_0_xxyzz_1, g_yyzzzz_0_xxzz_1, g_yyzzzz_0_xxzzz_1, g_yyzzzz_0_xyyy_1, g_yyzzzz_0_xyyyy_1, g_yyzzzz_0_xyyyz_1, g_yyzzzz_0_xyyz_1, g_yyzzzz_0_xyyzz_1, g_yyzzzz_0_xyzz_1, g_yyzzzz_0_xyzzz_1, g_yyzzzz_0_xzzz_1, g_yyzzzz_0_xzzzz_1, g_yyzzzz_0_yyyy_1, g_yyzzzz_0_yyyyy_1, g_yyzzzz_0_yyyyz_1, g_yyzzzz_0_yyyz_1, g_yyzzzz_0_yyyzz_1, g_yyzzzz_0_yyzz_1, g_yyzzzz_0_yyzzz_1, g_yyzzzz_0_yzzz_1, g_yyzzzz_0_yzzzz_1, g_yyzzzz_0_zzzz_1, g_yyzzzz_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzzz_0_xxxxx_0[i] = 5.0 * g_yyzzzz_0_xxxx_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxx_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxy_0[i] = 4.0 * g_yyzzzz_0_xxxy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxz_0[i] = 4.0 * g_yyzzzz_0_xxxz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyy_0[i] = 3.0 * g_yyzzzz_0_xxyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyz_0[i] = 3.0 * g_yyzzzz_0_xxyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxzz_0[i] = 3.0 * g_yyzzzz_0_xxzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyy_0[i] = 2.0 * g_yyzzzz_0_xyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyz_0[i] = 2.0 * g_yyzzzz_0_xyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyzz_0[i] = 2.0 * g_yyzzzz_0_xyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxzzz_0[i] = 2.0 * g_yyzzzz_0_xzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyy_0[i] = g_yyzzzz_0_yyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyz_0[i] = g_yyzzzz_0_yyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyzz_0[i] = g_yyzzzz_0_yyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyzzz_0[i] = g_yyzzzz_0_yzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xzzzz_0[i] = g_yyzzzz_0_zzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyy_0[i] = g_yyzzzz_0_yyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyz_0[i] = g_yyzzzz_0_yyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyzz_0[i] = g_yyzzzz_0_yyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyzzz_0[i] = g_yyzzzz_0_yyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yzzzz_0[i] = g_yyzzzz_0_yzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_zzzzz_0[i] = g_yyzzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 546-567 components of targeted buffer : KSH

    auto g_xyzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 546);

    auto g_xyzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 547);

    auto g_xyzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 548);

    auto g_xyzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 549);

    auto g_xyzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 550);

    auto g_xyzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 551);

    auto g_xyzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 552);

    auto g_xyzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 553);

    auto g_xyzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 554);

    auto g_xyzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 555);

    auto g_xyzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 556);

    auto g_xyzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 557);

    auto g_xyzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 558);

    auto g_xyzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 559);

    auto g_xyzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 560);

    auto g_xyzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 561);

    auto g_xyzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 562);

    auto g_xyzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 563);

    auto g_xyzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 564);

    auto g_xyzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 565);

    auto g_xyzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 566);

    #pragma omp simd aligned(g_xyzzzzz_0_xxxxx_0, g_xyzzzzz_0_xxxxy_0, g_xyzzzzz_0_xxxxz_0, g_xyzzzzz_0_xxxyy_0, g_xyzzzzz_0_xxxyz_0, g_xyzzzzz_0_xxxzz_0, g_xyzzzzz_0_xxyyy_0, g_xyzzzzz_0_xxyyz_0, g_xyzzzzz_0_xxyzz_0, g_xyzzzzz_0_xxzzz_0, g_xyzzzzz_0_xyyyy_0, g_xyzzzzz_0_xyyyz_0, g_xyzzzzz_0_xyyzz_0, g_xyzzzzz_0_xyzzz_0, g_xyzzzzz_0_xzzzz_0, g_xyzzzzz_0_yyyyy_0, g_xyzzzzz_0_yyyyz_0, g_xyzzzzz_0_yyyzz_0, g_xyzzzzz_0_yyzzz_0, g_xyzzzzz_0_yzzzz_0, g_xyzzzzz_0_zzzzz_0, g_xzzzzz_0_xxxxx_1, g_xzzzzz_0_xxxxz_1, g_xzzzzz_0_xxxzz_1, g_xzzzzz_0_xxzzz_1, g_xzzzzz_0_xzzzz_1, g_yzzzzz_0_xxxxy_1, g_yzzzzz_0_xxxy_1, g_yzzzzz_0_xxxyy_1, g_yzzzzz_0_xxxyz_1, g_yzzzzz_0_xxyy_1, g_yzzzzz_0_xxyyy_1, g_yzzzzz_0_xxyyz_1, g_yzzzzz_0_xxyz_1, g_yzzzzz_0_xxyzz_1, g_yzzzzz_0_xyyy_1, g_yzzzzz_0_xyyyy_1, g_yzzzzz_0_xyyyz_1, g_yzzzzz_0_xyyz_1, g_yzzzzz_0_xyyzz_1, g_yzzzzz_0_xyzz_1, g_yzzzzz_0_xyzzz_1, g_yzzzzz_0_yyyy_1, g_yzzzzz_0_yyyyy_1, g_yzzzzz_0_yyyyz_1, g_yzzzzz_0_yyyz_1, g_yzzzzz_0_yyyzz_1, g_yzzzzz_0_yyzz_1, g_yzzzzz_0_yyzzz_1, g_yzzzzz_0_yzzz_1, g_yzzzzz_0_yzzzz_1, g_yzzzzz_0_zzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzzz_0_xxxxx_0[i] = g_xzzzzz_0_xxxxx_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxxy_0[i] = 4.0 * g_yzzzzz_0_xxxy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxz_0[i] = g_xzzzzz_0_xxxxz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxyy_0[i] = 3.0 * g_yzzzzz_0_xxyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxyz_0[i] = 3.0 * g_yzzzzz_0_xxyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxzz_0[i] = g_xzzzzz_0_xxxzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxyyy_0[i] = 2.0 * g_yzzzzz_0_xyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyyz_0[i] = 2.0 * g_yzzzzz_0_xyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyzz_0[i] = 2.0 * g_yzzzzz_0_xyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxzzz_0[i] = g_xzzzzz_0_xxzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xyyyy_0[i] = g_yzzzzz_0_yyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyyz_0[i] = g_yzzzzz_0_yyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyzz_0[i] = g_yzzzzz_0_yyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyzzz_0[i] = g_yzzzzz_0_yzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xzzzz_0[i] = g_xzzzzz_0_xzzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_yyyyy_0[i] = g_yzzzzz_0_yyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyyz_0[i] = g_yzzzzz_0_yyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyzz_0[i] = g_yzzzzz_0_yyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyzzz_0[i] = g_yzzzzz_0_yyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yzzzz_0[i] = g_yzzzzz_0_yzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_zzzzz_0[i] = g_yzzzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 567-588 components of targeted buffer : KSH

    auto g_xzzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 567);

    auto g_xzzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 568);

    auto g_xzzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 569);

    auto g_xzzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 570);

    auto g_xzzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 571);

    auto g_xzzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 572);

    auto g_xzzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 573);

    auto g_xzzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 574);

    auto g_xzzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 575);

    auto g_xzzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 576);

    auto g_xzzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 577);

    auto g_xzzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 578);

    auto g_xzzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 579);

    auto g_xzzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 580);

    auto g_xzzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 581);

    auto g_xzzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 582);

    auto g_xzzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 583);

    auto g_xzzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 584);

    auto g_xzzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 585);

    auto g_xzzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 586);

    auto g_xzzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 587);

    #pragma omp simd aligned(g_xzzzzzz_0_xxxxx_0, g_xzzzzzz_0_xxxxy_0, g_xzzzzzz_0_xxxxz_0, g_xzzzzzz_0_xxxyy_0, g_xzzzzzz_0_xxxyz_0, g_xzzzzzz_0_xxxzz_0, g_xzzzzzz_0_xxyyy_0, g_xzzzzzz_0_xxyyz_0, g_xzzzzzz_0_xxyzz_0, g_xzzzzzz_0_xxzzz_0, g_xzzzzzz_0_xyyyy_0, g_xzzzzzz_0_xyyyz_0, g_xzzzzzz_0_xyyzz_0, g_xzzzzzz_0_xyzzz_0, g_xzzzzzz_0_xzzzz_0, g_xzzzzzz_0_yyyyy_0, g_xzzzzzz_0_yyyyz_0, g_xzzzzzz_0_yyyzz_0, g_xzzzzzz_0_yyzzz_0, g_xzzzzzz_0_yzzzz_0, g_xzzzzzz_0_zzzzz_0, g_zzzzzz_0_xxxx_1, g_zzzzzz_0_xxxxx_1, g_zzzzzz_0_xxxxy_1, g_zzzzzz_0_xxxxz_1, g_zzzzzz_0_xxxy_1, g_zzzzzz_0_xxxyy_1, g_zzzzzz_0_xxxyz_1, g_zzzzzz_0_xxxz_1, g_zzzzzz_0_xxxzz_1, g_zzzzzz_0_xxyy_1, g_zzzzzz_0_xxyyy_1, g_zzzzzz_0_xxyyz_1, g_zzzzzz_0_xxyz_1, g_zzzzzz_0_xxyzz_1, g_zzzzzz_0_xxzz_1, g_zzzzzz_0_xxzzz_1, g_zzzzzz_0_xyyy_1, g_zzzzzz_0_xyyyy_1, g_zzzzzz_0_xyyyz_1, g_zzzzzz_0_xyyz_1, g_zzzzzz_0_xyyzz_1, g_zzzzzz_0_xyzz_1, g_zzzzzz_0_xyzzz_1, g_zzzzzz_0_xzzz_1, g_zzzzzz_0_xzzzz_1, g_zzzzzz_0_yyyy_1, g_zzzzzz_0_yyyyy_1, g_zzzzzz_0_yyyyz_1, g_zzzzzz_0_yyyz_1, g_zzzzzz_0_yyyzz_1, g_zzzzzz_0_yyzz_1, g_zzzzzz_0_yyzzz_1, g_zzzzzz_0_yzzz_1, g_zzzzzz_0_yzzzz_1, g_zzzzzz_0_zzzz_1, g_zzzzzz_0_zzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzzz_0_xxxxx_0[i] = 5.0 * g_zzzzzz_0_xxxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxx_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxy_0[i] = 4.0 * g_zzzzzz_0_xxxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxz_0[i] = 4.0 * g_zzzzzz_0_xxxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyy_0[i] = 3.0 * g_zzzzzz_0_xxyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyz_0[i] = 3.0 * g_zzzzzz_0_xxyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxzz_0[i] = 3.0 * g_zzzzzz_0_xxzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyy_0[i] = 2.0 * g_zzzzzz_0_xyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyz_0[i] = 2.0 * g_zzzzzz_0_xyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyzz_0[i] = 2.0 * g_zzzzzz_0_xyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxzzz_0[i] = 2.0 * g_zzzzzz_0_xzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyy_0[i] = g_zzzzzz_0_yyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyz_0[i] = g_zzzzzz_0_yyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyzz_0[i] = g_zzzzzz_0_yyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyzzz_0[i] = g_zzzzzz_0_yzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xzzzz_0[i] = g_zzzzzz_0_zzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyy_0[i] = g_zzzzzz_0_yyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyz_0[i] = g_zzzzzz_0_yyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyzz_0[i] = g_zzzzzz_0_yyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyzzz_0[i] = g_zzzzzz_0_yyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yzzzz_0[i] = g_zzzzzz_0_yzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_zzzzz_0[i] = g_zzzzzz_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 588-609 components of targeted buffer : KSH

    auto g_yyyyyyy_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 588);

    auto g_yyyyyyy_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 589);

    auto g_yyyyyyy_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 590);

    auto g_yyyyyyy_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 591);

    auto g_yyyyyyy_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 592);

    auto g_yyyyyyy_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 593);

    auto g_yyyyyyy_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 594);

    auto g_yyyyyyy_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 595);

    auto g_yyyyyyy_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 596);

    auto g_yyyyyyy_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 597);

    auto g_yyyyyyy_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 598);

    auto g_yyyyyyy_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 599);

    auto g_yyyyyyy_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 600);

    auto g_yyyyyyy_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 601);

    auto g_yyyyyyy_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 602);

    auto g_yyyyyyy_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 603);

    auto g_yyyyyyy_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 604);

    auto g_yyyyyyy_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 605);

    auto g_yyyyyyy_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 606);

    auto g_yyyyyyy_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 607);

    auto g_yyyyyyy_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 608);

    #pragma omp simd aligned(g_yyyyy_0_xxxxx_0, g_yyyyy_0_xxxxx_1, g_yyyyy_0_xxxxy_0, g_yyyyy_0_xxxxy_1, g_yyyyy_0_xxxxz_0, g_yyyyy_0_xxxxz_1, g_yyyyy_0_xxxyy_0, g_yyyyy_0_xxxyy_1, g_yyyyy_0_xxxyz_0, g_yyyyy_0_xxxyz_1, g_yyyyy_0_xxxzz_0, g_yyyyy_0_xxxzz_1, g_yyyyy_0_xxyyy_0, g_yyyyy_0_xxyyy_1, g_yyyyy_0_xxyyz_0, g_yyyyy_0_xxyyz_1, g_yyyyy_0_xxyzz_0, g_yyyyy_0_xxyzz_1, g_yyyyy_0_xxzzz_0, g_yyyyy_0_xxzzz_1, g_yyyyy_0_xyyyy_0, g_yyyyy_0_xyyyy_1, g_yyyyy_0_xyyyz_0, g_yyyyy_0_xyyyz_1, g_yyyyy_0_xyyzz_0, g_yyyyy_0_xyyzz_1, g_yyyyy_0_xyzzz_0, g_yyyyy_0_xyzzz_1, g_yyyyy_0_xzzzz_0, g_yyyyy_0_xzzzz_1, g_yyyyy_0_yyyyy_0, g_yyyyy_0_yyyyy_1, g_yyyyy_0_yyyyz_0, g_yyyyy_0_yyyyz_1, g_yyyyy_0_yyyzz_0, g_yyyyy_0_yyyzz_1, g_yyyyy_0_yyzzz_0, g_yyyyy_0_yyzzz_1, g_yyyyy_0_yzzzz_0, g_yyyyy_0_yzzzz_1, g_yyyyy_0_zzzzz_0, g_yyyyy_0_zzzzz_1, g_yyyyyy_0_xxxx_1, g_yyyyyy_0_xxxxx_1, g_yyyyyy_0_xxxxy_1, g_yyyyyy_0_xxxxz_1, g_yyyyyy_0_xxxy_1, g_yyyyyy_0_xxxyy_1, g_yyyyyy_0_xxxyz_1, g_yyyyyy_0_xxxz_1, g_yyyyyy_0_xxxzz_1, g_yyyyyy_0_xxyy_1, g_yyyyyy_0_xxyyy_1, g_yyyyyy_0_xxyyz_1, g_yyyyyy_0_xxyz_1, g_yyyyyy_0_xxyzz_1, g_yyyyyy_0_xxzz_1, g_yyyyyy_0_xxzzz_1, g_yyyyyy_0_xyyy_1, g_yyyyyy_0_xyyyy_1, g_yyyyyy_0_xyyyz_1, g_yyyyyy_0_xyyz_1, g_yyyyyy_0_xyyzz_1, g_yyyyyy_0_xyzz_1, g_yyyyyy_0_xyzzz_1, g_yyyyyy_0_xzzz_1, g_yyyyyy_0_xzzzz_1, g_yyyyyy_0_yyyy_1, g_yyyyyy_0_yyyyy_1, g_yyyyyy_0_yyyyz_1, g_yyyyyy_0_yyyz_1, g_yyyyyy_0_yyyzz_1, g_yyyyyy_0_yyzz_1, g_yyyyyy_0_yyzzz_1, g_yyyyyy_0_yzzz_1, g_yyyyyy_0_yzzzz_1, g_yyyyyy_0_zzzz_1, g_yyyyyy_0_zzzzz_1, g_yyyyyyy_0_xxxxx_0, g_yyyyyyy_0_xxxxy_0, g_yyyyyyy_0_xxxxz_0, g_yyyyyyy_0_xxxyy_0, g_yyyyyyy_0_xxxyz_0, g_yyyyyyy_0_xxxzz_0, g_yyyyyyy_0_xxyyy_0, g_yyyyyyy_0_xxyyz_0, g_yyyyyyy_0_xxyzz_0, g_yyyyyyy_0_xxzzz_0, g_yyyyyyy_0_xyyyy_0, g_yyyyyyy_0_xyyyz_0, g_yyyyyyy_0_xyyzz_0, g_yyyyyyy_0_xyzzz_0, g_yyyyyyy_0_xzzzz_0, g_yyyyyyy_0_yyyyy_0, g_yyyyyyy_0_yyyyz_0, g_yyyyyyy_0_yyyzz_0, g_yyyyyyy_0_yyzzz_0, g_yyyyyyy_0_yzzzz_0, g_yyyyyyy_0_zzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyyy_0_xxxxx_0[i] = 6.0 * g_yyyyy_0_xxxxx_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxx_1[i] * fz_be_0 + g_yyyyyy_0_xxxxx_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxy_0[i] = 6.0 * g_yyyyy_0_xxxxy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxy_1[i] * fz_be_0 + g_yyyyyy_0_xxxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxz_0[i] = 6.0 * g_yyyyy_0_xxxxz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyy_0[i] = 6.0 * g_yyyyy_0_xxxyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyy_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyz_0[i] = 6.0 * g_yyyyy_0_xxxyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyz_1[i] * fz_be_0 + g_yyyyyy_0_xxxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxzz_0[i] = 6.0 * g_yyyyy_0_xxxzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyy_0[i] = 6.0 * g_yyyyy_0_xxyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyy_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xxyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyz_0[i] = 6.0 * g_yyyyy_0_xxyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyzz_0[i] = 6.0 * g_yyyyy_0_xxyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyzz_1[i] * fz_be_0 + g_yyyyyy_0_xxzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxzzz_0[i] = 6.0 * g_yyyyy_0_xxzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyy_0[i] = 6.0 * g_yyyyy_0_xyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyy_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_xyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyz_0[i] = 6.0 * g_yyyyy_0_xyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyzz_0[i] = 6.0 * g_yyyyy_0_xyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyzzz_0[i] = 6.0 * g_yyyyy_0_xyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyzzz_1[i] * fz_be_0 + g_yyyyyy_0_xzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xzzzz_0[i] = 6.0 * g_yyyyy_0_xzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyy_0[i] = 6.0 * g_yyyyy_0_yyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyy_1[i] * fz_be_0 + 5.0 * g_yyyyyy_0_yyyy_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyz_0[i] = 6.0 * g_yyyyy_0_yyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_yyyz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyzz_0[i] = 6.0 * g_yyyyy_0_yyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_yyzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyzzz_0[i] = 6.0 * g_yyyyy_0_yyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_yzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yzzzz_0[i] = 6.0 * g_yyyyy_0_yzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yzzzz_1[i] * fz_be_0 + g_yyyyyy_0_zzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_zzzzz_0[i] = 6.0 * g_yyyyy_0_zzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_zzzzz_1[i] * fz_be_0 + g_yyyyyy_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 609-630 components of targeted buffer : KSH

    auto g_yyyyyyz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 609);

    auto g_yyyyyyz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 610);

    auto g_yyyyyyz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 611);

    auto g_yyyyyyz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 612);

    auto g_yyyyyyz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 613);

    auto g_yyyyyyz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 614);

    auto g_yyyyyyz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 615);

    auto g_yyyyyyz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 616);

    auto g_yyyyyyz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 617);

    auto g_yyyyyyz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 618);

    auto g_yyyyyyz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 619);

    auto g_yyyyyyz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 620);

    auto g_yyyyyyz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 621);

    auto g_yyyyyyz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 622);

    auto g_yyyyyyz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 623);

    auto g_yyyyyyz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 624);

    auto g_yyyyyyz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 625);

    auto g_yyyyyyz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 626);

    auto g_yyyyyyz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 627);

    auto g_yyyyyyz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 628);

    auto g_yyyyyyz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 629);

    #pragma omp simd aligned(g_yyyyyy_0_xxxx_1, g_yyyyyy_0_xxxxx_1, g_yyyyyy_0_xxxxy_1, g_yyyyyy_0_xxxxz_1, g_yyyyyy_0_xxxy_1, g_yyyyyy_0_xxxyy_1, g_yyyyyy_0_xxxyz_1, g_yyyyyy_0_xxxz_1, g_yyyyyy_0_xxxzz_1, g_yyyyyy_0_xxyy_1, g_yyyyyy_0_xxyyy_1, g_yyyyyy_0_xxyyz_1, g_yyyyyy_0_xxyz_1, g_yyyyyy_0_xxyzz_1, g_yyyyyy_0_xxzz_1, g_yyyyyy_0_xxzzz_1, g_yyyyyy_0_xyyy_1, g_yyyyyy_0_xyyyy_1, g_yyyyyy_0_xyyyz_1, g_yyyyyy_0_xyyz_1, g_yyyyyy_0_xyyzz_1, g_yyyyyy_0_xyzz_1, g_yyyyyy_0_xyzzz_1, g_yyyyyy_0_xzzz_1, g_yyyyyy_0_xzzzz_1, g_yyyyyy_0_yyyy_1, g_yyyyyy_0_yyyyy_1, g_yyyyyy_0_yyyyz_1, g_yyyyyy_0_yyyz_1, g_yyyyyy_0_yyyzz_1, g_yyyyyy_0_yyzz_1, g_yyyyyy_0_yyzzz_1, g_yyyyyy_0_yzzz_1, g_yyyyyy_0_yzzzz_1, g_yyyyyy_0_zzzz_1, g_yyyyyy_0_zzzzz_1, g_yyyyyyz_0_xxxxx_0, g_yyyyyyz_0_xxxxy_0, g_yyyyyyz_0_xxxxz_0, g_yyyyyyz_0_xxxyy_0, g_yyyyyyz_0_xxxyz_0, g_yyyyyyz_0_xxxzz_0, g_yyyyyyz_0_xxyyy_0, g_yyyyyyz_0_xxyyz_0, g_yyyyyyz_0_xxyzz_0, g_yyyyyyz_0_xxzzz_0, g_yyyyyyz_0_xyyyy_0, g_yyyyyyz_0_xyyyz_0, g_yyyyyyz_0_xyyzz_0, g_yyyyyyz_0_xyzzz_0, g_yyyyyyz_0_xzzzz_0, g_yyyyyyz_0_yyyyy_0, g_yyyyyyz_0_yyyyz_0, g_yyyyyyz_0_yyyzz_0, g_yyyyyyz_0_yyzzz_0, g_yyyyyyz_0_yzzzz_0, g_yyyyyyz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyyz_0_xxxxx_0[i] = g_yyyyyy_0_xxxxx_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxy_0[i] = g_yyyyyy_0_xxxxy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxz_0[i] = g_yyyyyy_0_xxxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyy_0[i] = g_yyyyyy_0_xxxyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyz_0[i] = g_yyyyyy_0_xxxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxzz_0[i] = 2.0 * g_yyyyyy_0_xxxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyy_0[i] = g_yyyyyy_0_xxyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyz_0[i] = g_yyyyyy_0_xxyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyzz_0[i] = 2.0 * g_yyyyyy_0_xxyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxzzz_0[i] = 3.0 * g_yyyyyy_0_xxzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyy_0[i] = g_yyyyyy_0_xyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyz_0[i] = g_yyyyyy_0_xyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyzz_0[i] = 2.0 * g_yyyyyy_0_xyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyzzz_0[i] = 3.0 * g_yyyyyy_0_xyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xzzzz_0[i] = 4.0 * g_yyyyyy_0_xzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyy_0[i] = g_yyyyyy_0_yyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyz_0[i] = g_yyyyyy_0_yyyy_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyzz_0[i] = 2.0 * g_yyyyyy_0_yyyz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyzzz_0[i] = 3.0 * g_yyyyyy_0_yyzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yzzzz_0[i] = 4.0 * g_yyyyyy_0_yzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_zzzzz_0[i] = 5.0 * g_yyyyyy_0_zzzz_1[i] * fi_acd_0 + g_yyyyyy_0_zzzzz_1[i] * wa_z[i];
    }

    /// Set up 630-651 components of targeted buffer : KSH

    auto g_yyyyyzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 630);

    auto g_yyyyyzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 631);

    auto g_yyyyyzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 632);

    auto g_yyyyyzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 633);

    auto g_yyyyyzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 634);

    auto g_yyyyyzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 635);

    auto g_yyyyyzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 636);

    auto g_yyyyyzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 637);

    auto g_yyyyyzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 638);

    auto g_yyyyyzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 639);

    auto g_yyyyyzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 640);

    auto g_yyyyyzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 641);

    auto g_yyyyyzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 642);

    auto g_yyyyyzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 643);

    auto g_yyyyyzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 644);

    auto g_yyyyyzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 645);

    auto g_yyyyyzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 646);

    auto g_yyyyyzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 647);

    auto g_yyyyyzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 648);

    auto g_yyyyyzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 649);

    auto g_yyyyyzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 650);

    #pragma omp simd aligned(g_yyyyy_0_xxxxy_0, g_yyyyy_0_xxxxy_1, g_yyyyy_0_xxxyy_0, g_yyyyy_0_xxxyy_1, g_yyyyy_0_xxyyy_0, g_yyyyy_0_xxyyy_1, g_yyyyy_0_xyyyy_0, g_yyyyy_0_xyyyy_1, g_yyyyy_0_yyyyy_0, g_yyyyy_0_yyyyy_1, g_yyyyyz_0_xxxxy_1, g_yyyyyz_0_xxxyy_1, g_yyyyyz_0_xxyyy_1, g_yyyyyz_0_xyyyy_1, g_yyyyyz_0_yyyyy_1, g_yyyyyzz_0_xxxxx_0, g_yyyyyzz_0_xxxxy_0, g_yyyyyzz_0_xxxxz_0, g_yyyyyzz_0_xxxyy_0, g_yyyyyzz_0_xxxyz_0, g_yyyyyzz_0_xxxzz_0, g_yyyyyzz_0_xxyyy_0, g_yyyyyzz_0_xxyyz_0, g_yyyyyzz_0_xxyzz_0, g_yyyyyzz_0_xxzzz_0, g_yyyyyzz_0_xyyyy_0, g_yyyyyzz_0_xyyyz_0, g_yyyyyzz_0_xyyzz_0, g_yyyyyzz_0_xyzzz_0, g_yyyyyzz_0_xzzzz_0, g_yyyyyzz_0_yyyyy_0, g_yyyyyzz_0_yyyyz_0, g_yyyyyzz_0_yyyzz_0, g_yyyyyzz_0_yyzzz_0, g_yyyyyzz_0_yzzzz_0, g_yyyyyzz_0_zzzzz_0, g_yyyyzz_0_xxxxx_1, g_yyyyzz_0_xxxxz_1, g_yyyyzz_0_xxxyz_1, g_yyyyzz_0_xxxz_1, g_yyyyzz_0_xxxzz_1, g_yyyyzz_0_xxyyz_1, g_yyyyzz_0_xxyz_1, g_yyyyzz_0_xxyzz_1, g_yyyyzz_0_xxzz_1, g_yyyyzz_0_xxzzz_1, g_yyyyzz_0_xyyyz_1, g_yyyyzz_0_xyyz_1, g_yyyyzz_0_xyyzz_1, g_yyyyzz_0_xyzz_1, g_yyyyzz_0_xyzzz_1, g_yyyyzz_0_xzzz_1, g_yyyyzz_0_xzzzz_1, g_yyyyzz_0_yyyyz_1, g_yyyyzz_0_yyyz_1, g_yyyyzz_0_yyyzz_1, g_yyyyzz_0_yyzz_1, g_yyyyzz_0_yyzzz_1, g_yyyyzz_0_yzzz_1, g_yyyyzz_0_yzzzz_1, g_yyyyzz_0_zzzz_1, g_yyyyzz_0_zzzzz_1, g_yyyzz_0_xxxxx_0, g_yyyzz_0_xxxxx_1, g_yyyzz_0_xxxxz_0, g_yyyzz_0_xxxxz_1, g_yyyzz_0_xxxyz_0, g_yyyzz_0_xxxyz_1, g_yyyzz_0_xxxzz_0, g_yyyzz_0_xxxzz_1, g_yyyzz_0_xxyyz_0, g_yyyzz_0_xxyyz_1, g_yyyzz_0_xxyzz_0, g_yyyzz_0_xxyzz_1, g_yyyzz_0_xxzzz_0, g_yyyzz_0_xxzzz_1, g_yyyzz_0_xyyyz_0, g_yyyzz_0_xyyyz_1, g_yyyzz_0_xyyzz_0, g_yyyzz_0_xyyzz_1, g_yyyzz_0_xyzzz_0, g_yyyzz_0_xyzzz_1, g_yyyzz_0_xzzzz_0, g_yyyzz_0_xzzzz_1, g_yyyzz_0_yyyyz_0, g_yyyzz_0_yyyyz_1, g_yyyzz_0_yyyzz_0, g_yyyzz_0_yyyzz_1, g_yyyzz_0_yyzzz_0, g_yyyzz_0_yyzzz_1, g_yyyzz_0_yzzzz_0, g_yyyzz_0_yzzzz_1, g_yyyzz_0_zzzzz_0, g_yyyzz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyzz_0_xxxxx_0[i] = 4.0 * g_yyyzz_0_xxxxx_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxx_1[i] * fz_be_0 + g_yyyyzz_0_xxxxx_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxy_0[i] = g_yyyyy_0_xxxxy_0[i] * fbe_0 - g_yyyyy_0_xxxxy_1[i] * fz_be_0 + g_yyyyyz_0_xxxxy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxxz_0[i] = 4.0 * g_yyyzz_0_xxxxz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxyy_0[i] = g_yyyyy_0_xxxyy_0[i] * fbe_0 - g_yyyyy_0_xxxyy_1[i] * fz_be_0 + g_yyyyyz_0_xxxyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxyz_0[i] = 4.0 * g_yyyzz_0_xxxyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxyz_1[i] * fz_be_0 + g_yyyyzz_0_xxxz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxzz_0[i] = 4.0 * g_yyyzz_0_xxxzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyyy_0[i] = g_yyyyy_0_xxyyy_0[i] * fbe_0 - g_yyyyy_0_xxyyy_1[i] * fz_be_0 + g_yyyyyz_0_xxyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxyyz_0[i] = 4.0 * g_yyyzz_0_xxyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xxyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyzz_0[i] = 4.0 * g_yyyzz_0_xxyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyzz_1[i] * fz_be_0 + g_yyyyzz_0_xxzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxzzz_0[i] = 4.0 * g_yyyzz_0_xxzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyyy_0[i] = g_yyyyy_0_xyyyy_0[i] * fbe_0 - g_yyyyy_0_xyyyy_1[i] * fz_be_0 + g_yyyyyz_0_xyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xyyyz_0[i] = 4.0 * g_yyyzz_0_xyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_xyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyzz_0[i] = 4.0 * g_yyyzz_0_xyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyzzz_0[i] = 4.0 * g_yyyzz_0_xyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyzzz_1[i] * fz_be_0 + g_yyyyzz_0_xzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xzzzz_0[i] = 4.0 * g_yyyzz_0_xzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyyy_0[i] = g_yyyyy_0_yyyyy_0[i] * fbe_0 - g_yyyyy_0_yyyyy_1[i] * fz_be_0 + g_yyyyyz_0_yyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_yyyyz_0[i] = 4.0 * g_yyyzz_0_yyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_0_yyyz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyzz_0[i] = 4.0 * g_yyyzz_0_yyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_yyzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyzzz_0[i] = 4.0 * g_yyyzz_0_yyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_yzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yzzzz_0[i] = 4.0 * g_yyyzz_0_yzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yzzzz_1[i] * fz_be_0 + g_yyyyzz_0_zzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_zzzzz_0[i] = 4.0 * g_yyyzz_0_zzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_zzzzz_1[i] * fz_be_0 + g_yyyyzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 651-672 components of targeted buffer : KSH

    auto g_yyyyzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 651);

    auto g_yyyyzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 652);

    auto g_yyyyzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 653);

    auto g_yyyyzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 654);

    auto g_yyyyzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 655);

    auto g_yyyyzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 656);

    auto g_yyyyzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 657);

    auto g_yyyyzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 658);

    auto g_yyyyzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 659);

    auto g_yyyyzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 660);

    auto g_yyyyzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 661);

    auto g_yyyyzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 662);

    auto g_yyyyzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 663);

    auto g_yyyyzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 664);

    auto g_yyyyzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 665);

    auto g_yyyyzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 666);

    auto g_yyyyzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 667);

    auto g_yyyyzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 668);

    auto g_yyyyzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 669);

    auto g_yyyyzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 670);

    auto g_yyyyzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 671);

    #pragma omp simd aligned(g_yyyyz_0_xxxxy_0, g_yyyyz_0_xxxxy_1, g_yyyyz_0_xxxyy_0, g_yyyyz_0_xxxyy_1, g_yyyyz_0_xxyyy_0, g_yyyyz_0_xxyyy_1, g_yyyyz_0_xyyyy_0, g_yyyyz_0_xyyyy_1, g_yyyyz_0_yyyyy_0, g_yyyyz_0_yyyyy_1, g_yyyyzz_0_xxxxy_1, g_yyyyzz_0_xxxyy_1, g_yyyyzz_0_xxyyy_1, g_yyyyzz_0_xyyyy_1, g_yyyyzz_0_yyyyy_1, g_yyyyzzz_0_xxxxx_0, g_yyyyzzz_0_xxxxy_0, g_yyyyzzz_0_xxxxz_0, g_yyyyzzz_0_xxxyy_0, g_yyyyzzz_0_xxxyz_0, g_yyyyzzz_0_xxxzz_0, g_yyyyzzz_0_xxyyy_0, g_yyyyzzz_0_xxyyz_0, g_yyyyzzz_0_xxyzz_0, g_yyyyzzz_0_xxzzz_0, g_yyyyzzz_0_xyyyy_0, g_yyyyzzz_0_xyyyz_0, g_yyyyzzz_0_xyyzz_0, g_yyyyzzz_0_xyzzz_0, g_yyyyzzz_0_xzzzz_0, g_yyyyzzz_0_yyyyy_0, g_yyyyzzz_0_yyyyz_0, g_yyyyzzz_0_yyyzz_0, g_yyyyzzz_0_yyzzz_0, g_yyyyzzz_0_yzzzz_0, g_yyyyzzz_0_zzzzz_0, g_yyyzzz_0_xxxxx_1, g_yyyzzz_0_xxxxz_1, g_yyyzzz_0_xxxyz_1, g_yyyzzz_0_xxxz_1, g_yyyzzz_0_xxxzz_1, g_yyyzzz_0_xxyyz_1, g_yyyzzz_0_xxyz_1, g_yyyzzz_0_xxyzz_1, g_yyyzzz_0_xxzz_1, g_yyyzzz_0_xxzzz_1, g_yyyzzz_0_xyyyz_1, g_yyyzzz_0_xyyz_1, g_yyyzzz_0_xyyzz_1, g_yyyzzz_0_xyzz_1, g_yyyzzz_0_xyzzz_1, g_yyyzzz_0_xzzz_1, g_yyyzzz_0_xzzzz_1, g_yyyzzz_0_yyyyz_1, g_yyyzzz_0_yyyz_1, g_yyyzzz_0_yyyzz_1, g_yyyzzz_0_yyzz_1, g_yyyzzz_0_yyzzz_1, g_yyyzzz_0_yzzz_1, g_yyyzzz_0_yzzzz_1, g_yyyzzz_0_zzzz_1, g_yyyzzz_0_zzzzz_1, g_yyzzz_0_xxxxx_0, g_yyzzz_0_xxxxx_1, g_yyzzz_0_xxxxz_0, g_yyzzz_0_xxxxz_1, g_yyzzz_0_xxxyz_0, g_yyzzz_0_xxxyz_1, g_yyzzz_0_xxxzz_0, g_yyzzz_0_xxxzz_1, g_yyzzz_0_xxyyz_0, g_yyzzz_0_xxyyz_1, g_yyzzz_0_xxyzz_0, g_yyzzz_0_xxyzz_1, g_yyzzz_0_xxzzz_0, g_yyzzz_0_xxzzz_1, g_yyzzz_0_xyyyz_0, g_yyzzz_0_xyyyz_1, g_yyzzz_0_xyyzz_0, g_yyzzz_0_xyyzz_1, g_yyzzz_0_xyzzz_0, g_yyzzz_0_xyzzz_1, g_yyzzz_0_xzzzz_0, g_yyzzz_0_xzzzz_1, g_yyzzz_0_yyyyz_0, g_yyzzz_0_yyyyz_1, g_yyzzz_0_yyyzz_0, g_yyzzz_0_yyyzz_1, g_yyzzz_0_yyzzz_0, g_yyzzz_0_yyzzz_1, g_yyzzz_0_yzzzz_0, g_yyzzz_0_yzzzz_1, g_yyzzz_0_zzzzz_0, g_yyzzz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzzz_0_xxxxx_0[i] = 3.0 * g_yyzzz_0_xxxxx_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxx_1[i] * fz_be_0 + g_yyyzzz_0_xxxxx_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxy_0[i] = 2.0 * g_yyyyz_0_xxxxy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxxy_1[i] * fz_be_0 + g_yyyyzz_0_xxxxy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxxz_0[i] = 3.0 * g_yyzzz_0_xxxxz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxyy_0[i] = 2.0 * g_yyyyz_0_xxxyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxyy_1[i] * fz_be_0 + g_yyyyzz_0_xxxyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxyz_0[i] = 3.0 * g_yyzzz_0_xxxyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxyz_1[i] * fz_be_0 + g_yyyzzz_0_xxxz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxzz_0[i] = 3.0 * g_yyzzz_0_xxxzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyyy_0[i] = 2.0 * g_yyyyz_0_xxyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxyyy_1[i] * fz_be_0 + g_yyyyzz_0_xxyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxyyz_0[i] = 3.0 * g_yyzzz_0_xxyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xxyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyzz_0[i] = 3.0 * g_yyzzz_0_xxyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyzz_1[i] * fz_be_0 + g_yyyzzz_0_xxzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxzzz_0[i] = 3.0 * g_yyzzz_0_xxzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyyy_0[i] = 2.0 * g_yyyyz_0_xyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xyyyy_1[i] * fz_be_0 + g_yyyyzz_0_xyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xyyyz_0[i] = 3.0 * g_yyzzz_0_xyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_xyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyzz_0[i] = 3.0 * g_yyzzz_0_xyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyzzz_0[i] = 3.0 * g_yyzzz_0_xyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyzzz_1[i] * fz_be_0 + g_yyyzzz_0_xzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xzzzz_0[i] = 3.0 * g_yyzzz_0_xzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyyy_0[i] = 2.0 * g_yyyyz_0_yyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_yyyyy_1[i] * fz_be_0 + g_yyyyzz_0_yyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_yyyyz_0[i] = 3.0 * g_yyzzz_0_yyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_0_yyyz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyzz_0[i] = 3.0 * g_yyzzz_0_yyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_yyzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyzzz_0[i] = 3.0 * g_yyzzz_0_yyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_yzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yzzzz_0[i] = 3.0 * g_yyzzz_0_yzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yzzzz_1[i] * fz_be_0 + g_yyyzzz_0_zzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_zzzzz_0[i] = 3.0 * g_yyzzz_0_zzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_zzzzz_1[i] * fz_be_0 + g_yyyzzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 672-693 components of targeted buffer : KSH

    auto g_yyyzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 672);

    auto g_yyyzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 673);

    auto g_yyyzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 674);

    auto g_yyyzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 675);

    auto g_yyyzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 676);

    auto g_yyyzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 677);

    auto g_yyyzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 678);

    auto g_yyyzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 679);

    auto g_yyyzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 680);

    auto g_yyyzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 681);

    auto g_yyyzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 682);

    auto g_yyyzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 683);

    auto g_yyyzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 684);

    auto g_yyyzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 685);

    auto g_yyyzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 686);

    auto g_yyyzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 687);

    auto g_yyyzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 688);

    auto g_yyyzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 689);

    auto g_yyyzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 690);

    auto g_yyyzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 691);

    auto g_yyyzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 692);

    #pragma omp simd aligned(g_yyyzz_0_xxxxy_0, g_yyyzz_0_xxxxy_1, g_yyyzz_0_xxxyy_0, g_yyyzz_0_xxxyy_1, g_yyyzz_0_xxyyy_0, g_yyyzz_0_xxyyy_1, g_yyyzz_0_xyyyy_0, g_yyyzz_0_xyyyy_1, g_yyyzz_0_yyyyy_0, g_yyyzz_0_yyyyy_1, g_yyyzzz_0_xxxxy_1, g_yyyzzz_0_xxxyy_1, g_yyyzzz_0_xxyyy_1, g_yyyzzz_0_xyyyy_1, g_yyyzzz_0_yyyyy_1, g_yyyzzzz_0_xxxxx_0, g_yyyzzzz_0_xxxxy_0, g_yyyzzzz_0_xxxxz_0, g_yyyzzzz_0_xxxyy_0, g_yyyzzzz_0_xxxyz_0, g_yyyzzzz_0_xxxzz_0, g_yyyzzzz_0_xxyyy_0, g_yyyzzzz_0_xxyyz_0, g_yyyzzzz_0_xxyzz_0, g_yyyzzzz_0_xxzzz_0, g_yyyzzzz_0_xyyyy_0, g_yyyzzzz_0_xyyyz_0, g_yyyzzzz_0_xyyzz_0, g_yyyzzzz_0_xyzzz_0, g_yyyzzzz_0_xzzzz_0, g_yyyzzzz_0_yyyyy_0, g_yyyzzzz_0_yyyyz_0, g_yyyzzzz_0_yyyzz_0, g_yyyzzzz_0_yyzzz_0, g_yyyzzzz_0_yzzzz_0, g_yyyzzzz_0_zzzzz_0, g_yyzzzz_0_xxxxx_1, g_yyzzzz_0_xxxxz_1, g_yyzzzz_0_xxxyz_1, g_yyzzzz_0_xxxz_1, g_yyzzzz_0_xxxzz_1, g_yyzzzz_0_xxyyz_1, g_yyzzzz_0_xxyz_1, g_yyzzzz_0_xxyzz_1, g_yyzzzz_0_xxzz_1, g_yyzzzz_0_xxzzz_1, g_yyzzzz_0_xyyyz_1, g_yyzzzz_0_xyyz_1, g_yyzzzz_0_xyyzz_1, g_yyzzzz_0_xyzz_1, g_yyzzzz_0_xyzzz_1, g_yyzzzz_0_xzzz_1, g_yyzzzz_0_xzzzz_1, g_yyzzzz_0_yyyyz_1, g_yyzzzz_0_yyyz_1, g_yyzzzz_0_yyyzz_1, g_yyzzzz_0_yyzz_1, g_yyzzzz_0_yyzzz_1, g_yyzzzz_0_yzzz_1, g_yyzzzz_0_yzzzz_1, g_yyzzzz_0_zzzz_1, g_yyzzzz_0_zzzzz_1, g_yzzzz_0_xxxxx_0, g_yzzzz_0_xxxxx_1, g_yzzzz_0_xxxxz_0, g_yzzzz_0_xxxxz_1, g_yzzzz_0_xxxyz_0, g_yzzzz_0_xxxyz_1, g_yzzzz_0_xxxzz_0, g_yzzzz_0_xxxzz_1, g_yzzzz_0_xxyyz_0, g_yzzzz_0_xxyyz_1, g_yzzzz_0_xxyzz_0, g_yzzzz_0_xxyzz_1, g_yzzzz_0_xxzzz_0, g_yzzzz_0_xxzzz_1, g_yzzzz_0_xyyyz_0, g_yzzzz_0_xyyyz_1, g_yzzzz_0_xyyzz_0, g_yzzzz_0_xyyzz_1, g_yzzzz_0_xyzzz_0, g_yzzzz_0_xyzzz_1, g_yzzzz_0_xzzzz_0, g_yzzzz_0_xzzzz_1, g_yzzzz_0_yyyyz_0, g_yzzzz_0_yyyyz_1, g_yzzzz_0_yyyzz_0, g_yzzzz_0_yyyzz_1, g_yzzzz_0_yyzzz_0, g_yzzzz_0_yyzzz_1, g_yzzzz_0_yzzzz_0, g_yzzzz_0_yzzzz_1, g_yzzzz_0_zzzzz_0, g_yzzzz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzzz_0_xxxxx_0[i] = 2.0 * g_yzzzz_0_xxxxx_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxx_1[i] * fz_be_0 + g_yyzzzz_0_xxxxx_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxy_0[i] = 3.0 * g_yyyzz_0_xxxxy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxxy_1[i] * fz_be_0 + g_yyyzzz_0_xxxxy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxxz_0[i] = 2.0 * g_yzzzz_0_xxxxz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxyy_0[i] = 3.0 * g_yyyzz_0_xxxyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxyy_1[i] * fz_be_0 + g_yyyzzz_0_xxxyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxyz_0[i] = 2.0 * g_yzzzz_0_xxxyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxyz_1[i] * fz_be_0 + g_yyzzzz_0_xxxz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxzz_0[i] = 2.0 * g_yzzzz_0_xxxzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyyy_0[i] = 3.0 * g_yyyzz_0_xxyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxyyy_1[i] * fz_be_0 + g_yyyzzz_0_xxyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxyyz_0[i] = 2.0 * g_yzzzz_0_xxyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xxyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyzz_0[i] = 2.0 * g_yzzzz_0_xxyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyzz_1[i] * fz_be_0 + g_yyzzzz_0_xxzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxzzz_0[i] = 2.0 * g_yzzzz_0_xxzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyyy_0[i] = 3.0 * g_yyyzz_0_xyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xyyyy_1[i] * fz_be_0 + g_yyyzzz_0_xyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xyyyz_0[i] = 2.0 * g_yzzzz_0_xyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_xyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyzz_0[i] = 2.0 * g_yzzzz_0_xyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyzzz_0[i] = 2.0 * g_yzzzz_0_xyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyzzz_1[i] * fz_be_0 + g_yyzzzz_0_xzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xzzzz_0[i] = 2.0 * g_yzzzz_0_xzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyyy_0[i] = 3.0 * g_yyyzz_0_yyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_yyyyy_1[i] * fz_be_0 + g_yyyzzz_0_yyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_yyyyz_0[i] = 2.0 * g_yzzzz_0_yyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_0_yyyz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyzz_0[i] = 2.0 * g_yzzzz_0_yyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_yyzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyzzz_0[i] = 2.0 * g_yzzzz_0_yyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_yzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yzzzz_0[i] = 2.0 * g_yzzzz_0_yzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yzzzz_1[i] * fz_be_0 + g_yyzzzz_0_zzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_zzzzz_0[i] = 2.0 * g_yzzzz_0_zzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_zzzzz_1[i] * fz_be_0 + g_yyzzzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 693-714 components of targeted buffer : KSH

    auto g_yyzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 693);

    auto g_yyzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 694);

    auto g_yyzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 695);

    auto g_yyzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 696);

    auto g_yyzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 697);

    auto g_yyzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 698);

    auto g_yyzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 699);

    auto g_yyzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 700);

    auto g_yyzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 701);

    auto g_yyzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 702);

    auto g_yyzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 703);

    auto g_yyzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 704);

    auto g_yyzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 705);

    auto g_yyzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 706);

    auto g_yyzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 707);

    auto g_yyzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 708);

    auto g_yyzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 709);

    auto g_yyzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 710);

    auto g_yyzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 711);

    auto g_yyzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 712);

    auto g_yyzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 713);

    #pragma omp simd aligned(g_yyzzz_0_xxxxy_0, g_yyzzz_0_xxxxy_1, g_yyzzz_0_xxxyy_0, g_yyzzz_0_xxxyy_1, g_yyzzz_0_xxyyy_0, g_yyzzz_0_xxyyy_1, g_yyzzz_0_xyyyy_0, g_yyzzz_0_xyyyy_1, g_yyzzz_0_yyyyy_0, g_yyzzz_0_yyyyy_1, g_yyzzzz_0_xxxxy_1, g_yyzzzz_0_xxxyy_1, g_yyzzzz_0_xxyyy_1, g_yyzzzz_0_xyyyy_1, g_yyzzzz_0_yyyyy_1, g_yyzzzzz_0_xxxxx_0, g_yyzzzzz_0_xxxxy_0, g_yyzzzzz_0_xxxxz_0, g_yyzzzzz_0_xxxyy_0, g_yyzzzzz_0_xxxyz_0, g_yyzzzzz_0_xxxzz_0, g_yyzzzzz_0_xxyyy_0, g_yyzzzzz_0_xxyyz_0, g_yyzzzzz_0_xxyzz_0, g_yyzzzzz_0_xxzzz_0, g_yyzzzzz_0_xyyyy_0, g_yyzzzzz_0_xyyyz_0, g_yyzzzzz_0_xyyzz_0, g_yyzzzzz_0_xyzzz_0, g_yyzzzzz_0_xzzzz_0, g_yyzzzzz_0_yyyyy_0, g_yyzzzzz_0_yyyyz_0, g_yyzzzzz_0_yyyzz_0, g_yyzzzzz_0_yyzzz_0, g_yyzzzzz_0_yzzzz_0, g_yyzzzzz_0_zzzzz_0, g_yzzzzz_0_xxxxx_1, g_yzzzzz_0_xxxxz_1, g_yzzzzz_0_xxxyz_1, g_yzzzzz_0_xxxz_1, g_yzzzzz_0_xxxzz_1, g_yzzzzz_0_xxyyz_1, g_yzzzzz_0_xxyz_1, g_yzzzzz_0_xxyzz_1, g_yzzzzz_0_xxzz_1, g_yzzzzz_0_xxzzz_1, g_yzzzzz_0_xyyyz_1, g_yzzzzz_0_xyyz_1, g_yzzzzz_0_xyyzz_1, g_yzzzzz_0_xyzz_1, g_yzzzzz_0_xyzzz_1, g_yzzzzz_0_xzzz_1, g_yzzzzz_0_xzzzz_1, g_yzzzzz_0_yyyyz_1, g_yzzzzz_0_yyyz_1, g_yzzzzz_0_yyyzz_1, g_yzzzzz_0_yyzz_1, g_yzzzzz_0_yyzzz_1, g_yzzzzz_0_yzzz_1, g_yzzzzz_0_yzzzz_1, g_yzzzzz_0_zzzz_1, g_yzzzzz_0_zzzzz_1, g_zzzzz_0_xxxxx_0, g_zzzzz_0_xxxxx_1, g_zzzzz_0_xxxxz_0, g_zzzzz_0_xxxxz_1, g_zzzzz_0_xxxyz_0, g_zzzzz_0_xxxyz_1, g_zzzzz_0_xxxzz_0, g_zzzzz_0_xxxzz_1, g_zzzzz_0_xxyyz_0, g_zzzzz_0_xxyyz_1, g_zzzzz_0_xxyzz_0, g_zzzzz_0_xxyzz_1, g_zzzzz_0_xxzzz_0, g_zzzzz_0_xxzzz_1, g_zzzzz_0_xyyyz_0, g_zzzzz_0_xyyyz_1, g_zzzzz_0_xyyzz_0, g_zzzzz_0_xyyzz_1, g_zzzzz_0_xyzzz_0, g_zzzzz_0_xyzzz_1, g_zzzzz_0_xzzzz_0, g_zzzzz_0_xzzzz_1, g_zzzzz_0_yyyyz_0, g_zzzzz_0_yyyyz_1, g_zzzzz_0_yyyzz_0, g_zzzzz_0_yyyzz_1, g_zzzzz_0_yyzzz_0, g_zzzzz_0_yyzzz_1, g_zzzzz_0_yzzzz_0, g_zzzzz_0_yzzzz_1, g_zzzzz_0_zzzzz_0, g_zzzzz_0_zzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzzz_0_xxxxx_0[i] = g_zzzzz_0_xxxxx_0[i] * fbe_0 - g_zzzzz_0_xxxxx_1[i] * fz_be_0 + g_yzzzzz_0_xxxxx_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxy_0[i] = 4.0 * g_yyzzz_0_xxxxy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxxy_1[i] * fz_be_0 + g_yyzzzz_0_xxxxy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxxz_0[i] = g_zzzzz_0_xxxxz_0[i] * fbe_0 - g_zzzzz_0_xxxxz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxyy_0[i] = 4.0 * g_yyzzz_0_xxxyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxyy_1[i] * fz_be_0 + g_yyzzzz_0_xxxyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxyz_0[i] = g_zzzzz_0_xxxyz_0[i] * fbe_0 - g_zzzzz_0_xxxyz_1[i] * fz_be_0 + g_yzzzzz_0_xxxz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxzz_0[i] = g_zzzzz_0_xxxzz_0[i] * fbe_0 - g_zzzzz_0_xxxzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyyy_0[i] = 4.0 * g_yyzzz_0_xxyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxyyy_1[i] * fz_be_0 + g_yyzzzz_0_xxyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxyyz_0[i] = g_zzzzz_0_xxyyz_0[i] * fbe_0 - g_zzzzz_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xxyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyzz_0[i] = g_zzzzz_0_xxyzz_0[i] * fbe_0 - g_zzzzz_0_xxyzz_1[i] * fz_be_0 + g_yzzzzz_0_xxzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxzzz_0[i] = g_zzzzz_0_xxzzz_0[i] * fbe_0 - g_zzzzz_0_xxzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyyy_0[i] = 4.0 * g_yyzzz_0_xyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xyyyy_1[i] * fz_be_0 + g_yyzzzz_0_xyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xyyyz_0[i] = g_zzzzz_0_xyyyz_0[i] * fbe_0 - g_zzzzz_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_xyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyzz_0[i] = g_zzzzz_0_xyyzz_0[i] * fbe_0 - g_zzzzz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyzzz_0[i] = g_zzzzz_0_xyzzz_0[i] * fbe_0 - g_zzzzz_0_xyzzz_1[i] * fz_be_0 + g_yzzzzz_0_xzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xzzzz_0[i] = g_zzzzz_0_xzzzz_0[i] * fbe_0 - g_zzzzz_0_xzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyyy_0[i] = 4.0 * g_yyzzz_0_yyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_yyyyy_1[i] * fz_be_0 + g_yyzzzz_0_yyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_yyyyz_0[i] = g_zzzzz_0_yyyyz_0[i] * fbe_0 - g_zzzzz_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_0_yyyz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyzz_0[i] = g_zzzzz_0_yyyzz_0[i] * fbe_0 - g_zzzzz_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_yyzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyzzz_0[i] = g_zzzzz_0_yyzzz_0[i] * fbe_0 - g_zzzzz_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_yzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yzzzz_0[i] = g_zzzzz_0_yzzzz_0[i] * fbe_0 - g_zzzzz_0_yzzzz_1[i] * fz_be_0 + g_yzzzzz_0_zzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_zzzzz_0[i] = g_zzzzz_0_zzzzz_0[i] * fbe_0 - g_zzzzz_0_zzzzz_1[i] * fz_be_0 + g_yzzzzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 714-735 components of targeted buffer : KSH

    auto g_yzzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 714);

    auto g_yzzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 715);

    auto g_yzzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 716);

    auto g_yzzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 717);

    auto g_yzzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 718);

    auto g_yzzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 719);

    auto g_yzzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 720);

    auto g_yzzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 721);

    auto g_yzzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 722);

    auto g_yzzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 723);

    auto g_yzzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 724);

    auto g_yzzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 725);

    auto g_yzzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 726);

    auto g_yzzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 727);

    auto g_yzzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 728);

    auto g_yzzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 729);

    auto g_yzzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 730);

    auto g_yzzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 731);

    auto g_yzzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 732);

    auto g_yzzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 733);

    auto g_yzzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 734);

    #pragma omp simd aligned(g_yzzzzzz_0_xxxxx_0, g_yzzzzzz_0_xxxxy_0, g_yzzzzzz_0_xxxxz_0, g_yzzzzzz_0_xxxyy_0, g_yzzzzzz_0_xxxyz_0, g_yzzzzzz_0_xxxzz_0, g_yzzzzzz_0_xxyyy_0, g_yzzzzzz_0_xxyyz_0, g_yzzzzzz_0_xxyzz_0, g_yzzzzzz_0_xxzzz_0, g_yzzzzzz_0_xyyyy_0, g_yzzzzzz_0_xyyyz_0, g_yzzzzzz_0_xyyzz_0, g_yzzzzzz_0_xyzzz_0, g_yzzzzzz_0_xzzzz_0, g_yzzzzzz_0_yyyyy_0, g_yzzzzzz_0_yyyyz_0, g_yzzzzzz_0_yyyzz_0, g_yzzzzzz_0_yyzzz_0, g_yzzzzzz_0_yzzzz_0, g_yzzzzzz_0_zzzzz_0, g_zzzzzz_0_xxxx_1, g_zzzzzz_0_xxxxx_1, g_zzzzzz_0_xxxxy_1, g_zzzzzz_0_xxxxz_1, g_zzzzzz_0_xxxy_1, g_zzzzzz_0_xxxyy_1, g_zzzzzz_0_xxxyz_1, g_zzzzzz_0_xxxz_1, g_zzzzzz_0_xxxzz_1, g_zzzzzz_0_xxyy_1, g_zzzzzz_0_xxyyy_1, g_zzzzzz_0_xxyyz_1, g_zzzzzz_0_xxyz_1, g_zzzzzz_0_xxyzz_1, g_zzzzzz_0_xxzz_1, g_zzzzzz_0_xxzzz_1, g_zzzzzz_0_xyyy_1, g_zzzzzz_0_xyyyy_1, g_zzzzzz_0_xyyyz_1, g_zzzzzz_0_xyyz_1, g_zzzzzz_0_xyyzz_1, g_zzzzzz_0_xyzz_1, g_zzzzzz_0_xyzzz_1, g_zzzzzz_0_xzzz_1, g_zzzzzz_0_xzzzz_1, g_zzzzzz_0_yyyy_1, g_zzzzzz_0_yyyyy_1, g_zzzzzz_0_yyyyz_1, g_zzzzzz_0_yyyz_1, g_zzzzzz_0_yyyzz_1, g_zzzzzz_0_yyzz_1, g_zzzzzz_0_yyzzz_1, g_zzzzzz_0_yzzz_1, g_zzzzzz_0_yzzzz_1, g_zzzzzz_0_zzzz_1, g_zzzzzz_0_zzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzzz_0_xxxxx_0[i] = g_zzzzzz_0_xxxxx_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxy_0[i] = g_zzzzzz_0_xxxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxz_0[i] = g_zzzzzz_0_xxxxz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyy_0[i] = 2.0 * g_zzzzzz_0_xxxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyz_0[i] = g_zzzzzz_0_xxxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxzz_0[i] = g_zzzzzz_0_xxxzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyy_0[i] = 3.0 * g_zzzzzz_0_xxyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyz_0[i] = 2.0 * g_zzzzzz_0_xxyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyzz_0[i] = g_zzzzzz_0_xxzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxzzz_0[i] = g_zzzzzz_0_xxzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyy_0[i] = 4.0 * g_zzzzzz_0_xyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyz_0[i] = 3.0 * g_zzzzzz_0_xyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyzz_0[i] = 2.0 * g_zzzzzz_0_xyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyzzz_0[i] = g_zzzzzz_0_xzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xzzzz_0[i] = g_zzzzzz_0_xzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyy_0[i] = 5.0 * g_zzzzzz_0_yyyy_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyz_0[i] = 4.0 * g_zzzzzz_0_yyyz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyzz_0[i] = 3.0 * g_zzzzzz_0_yyzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyzzz_0[i] = 2.0 * g_zzzzzz_0_yzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yzzzz_0[i] = g_zzzzzz_0_zzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_zzzzz_0[i] = g_zzzzzz_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 735-756 components of targeted buffer : KSH

    auto g_zzzzzzz_0_xxxxx_0 = pbuffer.data(idx_eri_0_ksh + 735);

    auto g_zzzzzzz_0_xxxxy_0 = pbuffer.data(idx_eri_0_ksh + 736);

    auto g_zzzzzzz_0_xxxxz_0 = pbuffer.data(idx_eri_0_ksh + 737);

    auto g_zzzzzzz_0_xxxyy_0 = pbuffer.data(idx_eri_0_ksh + 738);

    auto g_zzzzzzz_0_xxxyz_0 = pbuffer.data(idx_eri_0_ksh + 739);

    auto g_zzzzzzz_0_xxxzz_0 = pbuffer.data(idx_eri_0_ksh + 740);

    auto g_zzzzzzz_0_xxyyy_0 = pbuffer.data(idx_eri_0_ksh + 741);

    auto g_zzzzzzz_0_xxyyz_0 = pbuffer.data(idx_eri_0_ksh + 742);

    auto g_zzzzzzz_0_xxyzz_0 = pbuffer.data(idx_eri_0_ksh + 743);

    auto g_zzzzzzz_0_xxzzz_0 = pbuffer.data(idx_eri_0_ksh + 744);

    auto g_zzzzzzz_0_xyyyy_0 = pbuffer.data(idx_eri_0_ksh + 745);

    auto g_zzzzzzz_0_xyyyz_0 = pbuffer.data(idx_eri_0_ksh + 746);

    auto g_zzzzzzz_0_xyyzz_0 = pbuffer.data(idx_eri_0_ksh + 747);

    auto g_zzzzzzz_0_xyzzz_0 = pbuffer.data(idx_eri_0_ksh + 748);

    auto g_zzzzzzz_0_xzzzz_0 = pbuffer.data(idx_eri_0_ksh + 749);

    auto g_zzzzzzz_0_yyyyy_0 = pbuffer.data(idx_eri_0_ksh + 750);

    auto g_zzzzzzz_0_yyyyz_0 = pbuffer.data(idx_eri_0_ksh + 751);

    auto g_zzzzzzz_0_yyyzz_0 = pbuffer.data(idx_eri_0_ksh + 752);

    auto g_zzzzzzz_0_yyzzz_0 = pbuffer.data(idx_eri_0_ksh + 753);

    auto g_zzzzzzz_0_yzzzz_0 = pbuffer.data(idx_eri_0_ksh + 754);

    auto g_zzzzzzz_0_zzzzz_0 = pbuffer.data(idx_eri_0_ksh + 755);

    #pragma omp simd aligned(g_zzzzz_0_xxxxx_0, g_zzzzz_0_xxxxx_1, g_zzzzz_0_xxxxy_0, g_zzzzz_0_xxxxy_1, g_zzzzz_0_xxxxz_0, g_zzzzz_0_xxxxz_1, g_zzzzz_0_xxxyy_0, g_zzzzz_0_xxxyy_1, g_zzzzz_0_xxxyz_0, g_zzzzz_0_xxxyz_1, g_zzzzz_0_xxxzz_0, g_zzzzz_0_xxxzz_1, g_zzzzz_0_xxyyy_0, g_zzzzz_0_xxyyy_1, g_zzzzz_0_xxyyz_0, g_zzzzz_0_xxyyz_1, g_zzzzz_0_xxyzz_0, g_zzzzz_0_xxyzz_1, g_zzzzz_0_xxzzz_0, g_zzzzz_0_xxzzz_1, g_zzzzz_0_xyyyy_0, g_zzzzz_0_xyyyy_1, g_zzzzz_0_xyyyz_0, g_zzzzz_0_xyyyz_1, g_zzzzz_0_xyyzz_0, g_zzzzz_0_xyyzz_1, g_zzzzz_0_xyzzz_0, g_zzzzz_0_xyzzz_1, g_zzzzz_0_xzzzz_0, g_zzzzz_0_xzzzz_1, g_zzzzz_0_yyyyy_0, g_zzzzz_0_yyyyy_1, g_zzzzz_0_yyyyz_0, g_zzzzz_0_yyyyz_1, g_zzzzz_0_yyyzz_0, g_zzzzz_0_yyyzz_1, g_zzzzz_0_yyzzz_0, g_zzzzz_0_yyzzz_1, g_zzzzz_0_yzzzz_0, g_zzzzz_0_yzzzz_1, g_zzzzz_0_zzzzz_0, g_zzzzz_0_zzzzz_1, g_zzzzzz_0_xxxx_1, g_zzzzzz_0_xxxxx_1, g_zzzzzz_0_xxxxy_1, g_zzzzzz_0_xxxxz_1, g_zzzzzz_0_xxxy_1, g_zzzzzz_0_xxxyy_1, g_zzzzzz_0_xxxyz_1, g_zzzzzz_0_xxxz_1, g_zzzzzz_0_xxxzz_1, g_zzzzzz_0_xxyy_1, g_zzzzzz_0_xxyyy_1, g_zzzzzz_0_xxyyz_1, g_zzzzzz_0_xxyz_1, g_zzzzzz_0_xxyzz_1, g_zzzzzz_0_xxzz_1, g_zzzzzz_0_xxzzz_1, g_zzzzzz_0_xyyy_1, g_zzzzzz_0_xyyyy_1, g_zzzzzz_0_xyyyz_1, g_zzzzzz_0_xyyz_1, g_zzzzzz_0_xyyzz_1, g_zzzzzz_0_xyzz_1, g_zzzzzz_0_xyzzz_1, g_zzzzzz_0_xzzz_1, g_zzzzzz_0_xzzzz_1, g_zzzzzz_0_yyyy_1, g_zzzzzz_0_yyyyy_1, g_zzzzzz_0_yyyyz_1, g_zzzzzz_0_yyyz_1, g_zzzzzz_0_yyyzz_1, g_zzzzzz_0_yyzz_1, g_zzzzzz_0_yyzzz_1, g_zzzzzz_0_yzzz_1, g_zzzzzz_0_yzzzz_1, g_zzzzzz_0_zzzz_1, g_zzzzzz_0_zzzzz_1, g_zzzzzzz_0_xxxxx_0, g_zzzzzzz_0_xxxxy_0, g_zzzzzzz_0_xxxxz_0, g_zzzzzzz_0_xxxyy_0, g_zzzzzzz_0_xxxyz_0, g_zzzzzzz_0_xxxzz_0, g_zzzzzzz_0_xxyyy_0, g_zzzzzzz_0_xxyyz_0, g_zzzzzzz_0_xxyzz_0, g_zzzzzzz_0_xxzzz_0, g_zzzzzzz_0_xyyyy_0, g_zzzzzzz_0_xyyyz_0, g_zzzzzzz_0_xyyzz_0, g_zzzzzzz_0_xyzzz_0, g_zzzzzzz_0_xzzzz_0, g_zzzzzzz_0_yyyyy_0, g_zzzzzzz_0_yyyyz_0, g_zzzzzzz_0_yyyzz_0, g_zzzzzzz_0_yyzzz_0, g_zzzzzzz_0_yzzzz_0, g_zzzzzzz_0_zzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzzz_0_xxxxx_0[i] = 6.0 * g_zzzzz_0_xxxxx_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxx_1[i] * fz_be_0 + g_zzzzzz_0_xxxxx_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxy_0[i] = 6.0 * g_zzzzz_0_xxxxy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxy_1[i] * fz_be_0 + g_zzzzzz_0_xxxxy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxz_0[i] = 6.0 * g_zzzzz_0_xxxxz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxz_1[i] * fz_be_0 + g_zzzzzz_0_xxxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyy_0[i] = 6.0 * g_zzzzz_0_xxxyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyy_1[i] * fz_be_0 + g_zzzzzz_0_xxxyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyz_0[i] = 6.0 * g_zzzzz_0_xxxyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyz_1[i] * fz_be_0 + g_zzzzzz_0_xxxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxzz_0[i] = 6.0 * g_zzzzz_0_xxxzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyy_0[i] = 6.0 * g_zzzzz_0_xxyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyy_1[i] * fz_be_0 + g_zzzzzz_0_xxyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyz_0[i] = 6.0 * g_zzzzz_0_xxyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyz_1[i] * fz_be_0 + g_zzzzzz_0_xxyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyzz_0[i] = 6.0 * g_zzzzz_0_xxyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxzzz_0[i] = 6.0 * g_zzzzz_0_xxzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xxzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyy_0[i] = 6.0 * g_zzzzz_0_xyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyy_1[i] * fz_be_0 + g_zzzzzz_0_xyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyz_0[i] = 6.0 * g_zzzzz_0_xyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyz_1[i] * fz_be_0 + g_zzzzzz_0_xyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyzz_0[i] = 6.0 * g_zzzzz_0_xyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyzzz_0[i] = 6.0 * g_zzzzz_0_xyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xzzzz_0[i] = 6.0 * g_zzzzz_0_xzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_xzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyy_0[i] = 6.0 * g_zzzzz_0_yyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyy_1[i] * fz_be_0 + g_zzzzzz_0_yyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyz_0[i] = 6.0 * g_zzzzz_0_yyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyz_1[i] * fz_be_0 + g_zzzzzz_0_yyyy_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyzz_0[i] = 6.0 * g_zzzzz_0_yyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_yyyz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyzzz_0[i] = 6.0 * g_zzzzz_0_yyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_yyzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yzzzz_0[i] = 6.0 * g_zzzzz_0_yzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_yzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_zzzzz_0[i] = 6.0 * g_zzzzz_0_zzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_zzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_0_zzzz_1[i] * fi_acd_0 + g_zzzzzz_0_zzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

