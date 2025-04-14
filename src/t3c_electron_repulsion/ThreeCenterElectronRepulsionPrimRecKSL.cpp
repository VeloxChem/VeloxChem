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

#include "ThreeCenterElectronRepulsionPrimRecKSL.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ksl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksl,
                                 size_t idx_eri_0_hsl,
                                 size_t idx_eri_1_hsl,
                                 size_t idx_eri_1_isk,
                                 size_t idx_eri_1_isl,
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

    /// Set up components of auxilary buffer : HSL

    auto g_xxxxx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl);

    auto g_xxxxx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 1);

    auto g_xxxxx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 2);

    auto g_xxxxx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 3);

    auto g_xxxxx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 4);

    auto g_xxxxx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 5);

    auto g_xxxxx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 6);

    auto g_xxxxx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 7);

    auto g_xxxxx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 8);

    auto g_xxxxx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 9);

    auto g_xxxxx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 10);

    auto g_xxxxx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 11);

    auto g_xxxxx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 12);

    auto g_xxxxx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 13);

    auto g_xxxxx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 14);

    auto g_xxxxx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 15);

    auto g_xxxxx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 16);

    auto g_xxxxx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 17);

    auto g_xxxxx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 18);

    auto g_xxxxx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 19);

    auto g_xxxxx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 20);

    auto g_xxxxx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 21);

    auto g_xxxxx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 22);

    auto g_xxxxx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 23);

    auto g_xxxxx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 24);

    auto g_xxxxx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 25);

    auto g_xxxxx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 26);

    auto g_xxxxx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 27);

    auto g_xxxxx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 28);

    auto g_xxxxx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 29);

    auto g_xxxxx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 30);

    auto g_xxxxx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 31);

    auto g_xxxxx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 32);

    auto g_xxxxx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 33);

    auto g_xxxxx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 34);

    auto g_xxxxx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 35);

    auto g_xxxxx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 36);

    auto g_xxxxx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 37);

    auto g_xxxxx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 38);

    auto g_xxxxx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 39);

    auto g_xxxxx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 40);

    auto g_xxxxx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 41);

    auto g_xxxxx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 42);

    auto g_xxxxx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 43);

    auto g_xxxxx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 44);

    auto g_xxxxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 45);

    auto g_xxxxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 47);

    auto g_xxxxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 50);

    auto g_xxxxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 54);

    auto g_xxxxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 59);

    auto g_xxxxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 65);

    auto g_xxxxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 72);

    auto g_xxxxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 80);

    auto g_xxxxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 90);

    auto g_xxxxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 91);

    auto g_xxxxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 93);

    auto g_xxxxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 96);

    auto g_xxxxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 100);

    auto g_xxxxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 105);

    auto g_xxxxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 111);

    auto g_xxxxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 118);

    auto g_xxxyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 135);

    auto g_xxxyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 136);

    auto g_xxxyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 137);

    auto g_xxxyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 138);

    auto g_xxxyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 139);

    auto g_xxxyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 140);

    auto g_xxxyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 141);

    auto g_xxxyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 142);

    auto g_xxxyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 143);

    auto g_xxxyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 144);

    auto g_xxxyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 145);

    auto g_xxxyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 146);

    auto g_xxxyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 147);

    auto g_xxxyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 148);

    auto g_xxxyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 149);

    auto g_xxxyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 150);

    auto g_xxxyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 151);

    auto g_xxxyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 152);

    auto g_xxxyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 153);

    auto g_xxxyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 154);

    auto g_xxxyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 155);

    auto g_xxxyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 156);

    auto g_xxxyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 157);

    auto g_xxxyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 158);

    auto g_xxxyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 159);

    auto g_xxxyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 160);

    auto g_xxxyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 161);

    auto g_xxxyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 162);

    auto g_xxxyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 163);

    auto g_xxxyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 164);

    auto g_xxxyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 165);

    auto g_xxxyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 166);

    auto g_xxxyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 167);

    auto g_xxxyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 168);

    auto g_xxxyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 169);

    auto g_xxxyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 170);

    auto g_xxxyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 171);

    auto g_xxxyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 172);

    auto g_xxxyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 173);

    auto g_xxxyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 174);

    auto g_xxxyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 175);

    auto g_xxxyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 176);

    auto g_xxxyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 177);

    auto g_xxxyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 178);

    auto g_xxxyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 179);

    auto g_xxxzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 225);

    auto g_xxxzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 226);

    auto g_xxxzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 227);

    auto g_xxxzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 228);

    auto g_xxxzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 229);

    auto g_xxxzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 230);

    auto g_xxxzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 231);

    auto g_xxxzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 232);

    auto g_xxxzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 233);

    auto g_xxxzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 234);

    auto g_xxxzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 235);

    auto g_xxxzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 236);

    auto g_xxxzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 237);

    auto g_xxxzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 238);

    auto g_xxxzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 239);

    auto g_xxxzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 240);

    auto g_xxxzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 241);

    auto g_xxxzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 242);

    auto g_xxxzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 243);

    auto g_xxxzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 244);

    auto g_xxxzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 245);

    auto g_xxxzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 246);

    auto g_xxxzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 247);

    auto g_xxxzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 248);

    auto g_xxxzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 249);

    auto g_xxxzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 250);

    auto g_xxxzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 251);

    auto g_xxxzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 252);

    auto g_xxxzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 253);

    auto g_xxxzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 254);

    auto g_xxxzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 255);

    auto g_xxxzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 256);

    auto g_xxxzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 257);

    auto g_xxxzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 258);

    auto g_xxxzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 259);

    auto g_xxxzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 260);

    auto g_xxxzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 261);

    auto g_xxxzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 262);

    auto g_xxxzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 263);

    auto g_xxxzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 264);

    auto g_xxxzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 265);

    auto g_xxxzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 266);

    auto g_xxxzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 267);

    auto g_xxxzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 268);

    auto g_xxxzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 269);

    auto g_xxyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 270);

    auto g_xxyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 271);

    auto g_xxyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 272);

    auto g_xxyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 273);

    auto g_xxyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 274);

    auto g_xxyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 275);

    auto g_xxyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 276);

    auto g_xxyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 277);

    auto g_xxyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 278);

    auto g_xxyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 279);

    auto g_xxyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 280);

    auto g_xxyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 281);

    auto g_xxyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 282);

    auto g_xxyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 283);

    auto g_xxyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 284);

    auto g_xxyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 285);

    auto g_xxyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 286);

    auto g_xxyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 287);

    auto g_xxyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 288);

    auto g_xxyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 289);

    auto g_xxyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 290);

    auto g_xxyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 291);

    auto g_xxyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 292);

    auto g_xxyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 293);

    auto g_xxyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 294);

    auto g_xxyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 295);

    auto g_xxyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 296);

    auto g_xxyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 297);

    auto g_xxyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 298);

    auto g_xxyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 299);

    auto g_xxyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 300);

    auto g_xxyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 301);

    auto g_xxyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 302);

    auto g_xxyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 303);

    auto g_xxyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 304);

    auto g_xxyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 305);

    auto g_xxyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 306);

    auto g_xxyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 307);

    auto g_xxyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 308);

    auto g_xxyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 309);

    auto g_xxyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 310);

    auto g_xxyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 311);

    auto g_xxyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 312);

    auto g_xxyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 313);

    auto g_xxyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 314);

    auto g_xxyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 316);

    auto g_xxyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 318);

    auto g_xxyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 321);

    auto g_xxyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 325);

    auto g_xxyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 330);

    auto g_xxyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 336);

    auto g_xxyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 343);

    auto g_xxyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 360);

    auto g_xxyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 362);

    auto g_xxyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 365);

    auto g_xxyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 369);

    auto g_xxyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 374);

    auto g_xxyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 380);

    auto g_xxyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 387);

    auto g_xxyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 395);

    auto g_xxzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 405);

    auto g_xxzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 406);

    auto g_xxzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 407);

    auto g_xxzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 408);

    auto g_xxzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 409);

    auto g_xxzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 410);

    auto g_xxzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 411);

    auto g_xxzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 412);

    auto g_xxzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 413);

    auto g_xxzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 414);

    auto g_xxzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 415);

    auto g_xxzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 416);

    auto g_xxzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 417);

    auto g_xxzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 418);

    auto g_xxzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 419);

    auto g_xxzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 420);

    auto g_xxzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 421);

    auto g_xxzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 422);

    auto g_xxzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 423);

    auto g_xxzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 424);

    auto g_xxzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 425);

    auto g_xxzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 426);

    auto g_xxzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 427);

    auto g_xxzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 428);

    auto g_xxzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 429);

    auto g_xxzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 430);

    auto g_xxzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 431);

    auto g_xxzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 432);

    auto g_xxzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 433);

    auto g_xxzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 434);

    auto g_xxzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 435);

    auto g_xxzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 436);

    auto g_xxzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 437);

    auto g_xxzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 438);

    auto g_xxzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 439);

    auto g_xxzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 440);

    auto g_xxzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 441);

    auto g_xxzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 442);

    auto g_xxzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 443);

    auto g_xxzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 444);

    auto g_xxzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 445);

    auto g_xxzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 446);

    auto g_xxzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 447);

    auto g_xxzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 448);

    auto g_xxzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 449);

    auto g_xyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 451);

    auto g_xyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 453);

    auto g_xyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 454);

    auto g_xyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 456);

    auto g_xyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 457);

    auto g_xyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 458);

    auto g_xyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 460);

    auto g_xyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 461);

    auto g_xyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 462);

    auto g_xyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 463);

    auto g_xyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 465);

    auto g_xyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 466);

    auto g_xyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 467);

    auto g_xyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 468);

    auto g_xyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 469);

    auto g_xyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 471);

    auto g_xyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 472);

    auto g_xyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 473);

    auto g_xyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 474);

    auto g_xyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 475);

    auto g_xyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 476);

    auto g_xyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 478);

    auto g_xyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 479);

    auto g_xyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 480);

    auto g_xyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 481);

    auto g_xyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 482);

    auto g_xyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 483);

    auto g_xyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 484);

    auto g_xyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 486);

    auto g_xyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 487);

    auto g_xyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 488);

    auto g_xyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 489);

    auto g_xyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 490);

    auto g_xyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 491);

    auto g_xyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 492);

    auto g_xyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 493);

    auto g_xyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 494);

    auto g_xyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 544);

    auto g_xyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 547);

    auto g_xyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 548);

    auto g_xyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 551);

    auto g_xyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 552);

    auto g_xyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 553);

    auto g_xyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 556);

    auto g_xyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 557);

    auto g_xyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 558);

    auto g_xyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 559);

    auto g_xyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 562);

    auto g_xyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 563);

    auto g_xyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 564);

    auto g_xyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 565);

    auto g_xyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 566);

    auto g_xyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 569);

    auto g_xyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 570);

    auto g_xyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 571);

    auto g_xyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 572);

    auto g_xyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 573);

    auto g_xyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 574);

    auto g_xyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 576);

    auto g_xyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 577);

    auto g_xyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 578);

    auto g_xyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 579);

    auto g_xyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 580);

    auto g_xyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 581);

    auto g_xyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 582);

    auto g_xyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 583);

    auto g_xyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 584);

    auto g_xzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 632);

    auto g_xzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 634);

    auto g_xzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 635);

    auto g_xzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 637);

    auto g_xzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 638);

    auto g_xzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 639);

    auto g_xzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 641);

    auto g_xzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 642);

    auto g_xzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 643);

    auto g_xzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 644);

    auto g_xzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 646);

    auto g_xzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 647);

    auto g_xzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 648);

    auto g_xzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 649);

    auto g_xzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 650);

    auto g_xzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 652);

    auto g_xzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 653);

    auto g_xzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 654);

    auto g_xzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 655);

    auto g_xzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 656);

    auto g_xzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 657);

    auto g_xzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 659);

    auto g_xzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 660);

    auto g_xzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 661);

    auto g_xzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 662);

    auto g_xzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 663);

    auto g_xzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 664);

    auto g_xzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 665);

    auto g_xzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 666);

    auto g_xzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 667);

    auto g_xzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 668);

    auto g_xzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 669);

    auto g_xzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 670);

    auto g_xzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 671);

    auto g_xzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 672);

    auto g_xzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 673);

    auto g_xzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 674);

    auto g_yyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 675);

    auto g_yyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 676);

    auto g_yyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 677);

    auto g_yyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 678);

    auto g_yyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 679);

    auto g_yyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 680);

    auto g_yyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 681);

    auto g_yyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 682);

    auto g_yyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 683);

    auto g_yyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 684);

    auto g_yyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 685);

    auto g_yyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 686);

    auto g_yyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 687);

    auto g_yyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 688);

    auto g_yyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 689);

    auto g_yyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 690);

    auto g_yyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 691);

    auto g_yyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 692);

    auto g_yyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 693);

    auto g_yyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 694);

    auto g_yyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 695);

    auto g_yyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 696);

    auto g_yyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 697);

    auto g_yyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 698);

    auto g_yyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 699);

    auto g_yyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 700);

    auto g_yyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 701);

    auto g_yyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 702);

    auto g_yyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 703);

    auto g_yyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 704);

    auto g_yyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 705);

    auto g_yyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 706);

    auto g_yyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 707);

    auto g_yyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 708);

    auto g_yyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 709);

    auto g_yyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 710);

    auto g_yyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 711);

    auto g_yyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 712);

    auto g_yyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 713);

    auto g_yyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 714);

    auto g_yyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 715);

    auto g_yyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 716);

    auto g_yyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 717);

    auto g_yyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 718);

    auto g_yyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 719);

    auto g_yyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 721);

    auto g_yyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 723);

    auto g_yyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 726);

    auto g_yyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 730);

    auto g_yyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 735);

    auto g_yyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 741);

    auto g_yyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 748);

    auto g_yyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 756);

    auto g_yyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 765);

    auto g_yyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 766);

    auto g_yyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 767);

    auto g_yyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 768);

    auto g_yyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 769);

    auto g_yyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 770);

    auto g_yyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 771);

    auto g_yyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 772);

    auto g_yyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 773);

    auto g_yyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 774);

    auto g_yyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 775);

    auto g_yyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 776);

    auto g_yyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 777);

    auto g_yyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 778);

    auto g_yyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 779);

    auto g_yyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 780);

    auto g_yyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 781);

    auto g_yyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 782);

    auto g_yyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 783);

    auto g_yyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 784);

    auto g_yyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 785);

    auto g_yyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 786);

    auto g_yyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 787);

    auto g_yyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 788);

    auto g_yyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 789);

    auto g_yyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 790);

    auto g_yyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 791);

    auto g_yyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 792);

    auto g_yyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 793);

    auto g_yyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 794);

    auto g_yyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 795);

    auto g_yyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 796);

    auto g_yyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 797);

    auto g_yyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 798);

    auto g_yyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 799);

    auto g_yyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 800);

    auto g_yyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 801);

    auto g_yyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 802);

    auto g_yyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 803);

    auto g_yyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 804);

    auto g_yyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 805);

    auto g_yyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 806);

    auto g_yyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 807);

    auto g_yyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 808);

    auto g_yyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 809);

    auto g_yyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 810);

    auto g_yyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 811);

    auto g_yyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 812);

    auto g_yyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 813);

    auto g_yyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 814);

    auto g_yyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 815);

    auto g_yyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 816);

    auto g_yyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 817);

    auto g_yyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 818);

    auto g_yyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 819);

    auto g_yyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 820);

    auto g_yyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 821);

    auto g_yyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 822);

    auto g_yyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 823);

    auto g_yyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 824);

    auto g_yyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 825);

    auto g_yyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 826);

    auto g_yyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 827);

    auto g_yyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 828);

    auto g_yyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 829);

    auto g_yyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 830);

    auto g_yyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 831);

    auto g_yyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 832);

    auto g_yyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 833);

    auto g_yyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 834);

    auto g_yyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 835);

    auto g_yyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 836);

    auto g_yyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 837);

    auto g_yyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 838);

    auto g_yyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 839);

    auto g_yyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 840);

    auto g_yyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 841);

    auto g_yyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 842);

    auto g_yyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 843);

    auto g_yyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 844);

    auto g_yyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 845);

    auto g_yyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 846);

    auto g_yyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 847);

    auto g_yyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 848);

    auto g_yyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 849);

    auto g_yyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 850);

    auto g_yyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 851);

    auto g_yyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 852);

    auto g_yyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 853);

    auto g_yyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 854);

    auto g_yzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 855);

    auto g_yzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 857);

    auto g_yzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 859);

    auto g_yzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 860);

    auto g_yzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 862);

    auto g_yzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 863);

    auto g_yzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 864);

    auto g_yzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 866);

    auto g_yzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 867);

    auto g_yzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 868);

    auto g_yzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 869);

    auto g_yzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 871);

    auto g_yzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 872);

    auto g_yzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 873);

    auto g_yzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 874);

    auto g_yzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 875);

    auto g_yzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 877);

    auto g_yzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 878);

    auto g_yzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 879);

    auto g_yzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 880);

    auto g_yzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 881);

    auto g_yzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 882);

    auto g_yzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 884);

    auto g_yzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 885);

    auto g_yzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 886);

    auto g_yzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 887);

    auto g_yzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 888);

    auto g_yzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 889);

    auto g_yzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 890);

    auto g_yzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 892);

    auto g_yzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 893);

    auto g_yzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 894);

    auto g_yzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 895);

    auto g_yzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 896);

    auto g_yzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 897);

    auto g_yzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 898);

    auto g_yzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 899);

    auto g_zzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 900);

    auto g_zzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 901);

    auto g_zzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 902);

    auto g_zzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 903);

    auto g_zzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 904);

    auto g_zzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 905);

    auto g_zzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 906);

    auto g_zzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 907);

    auto g_zzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 908);

    auto g_zzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 909);

    auto g_zzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 910);

    auto g_zzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 911);

    auto g_zzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 912);

    auto g_zzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 913);

    auto g_zzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 914);

    auto g_zzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 915);

    auto g_zzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 916);

    auto g_zzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 917);

    auto g_zzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 918);

    auto g_zzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 919);

    auto g_zzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 920);

    auto g_zzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 921);

    auto g_zzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 922);

    auto g_zzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 923);

    auto g_zzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 924);

    auto g_zzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 925);

    auto g_zzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 926);

    auto g_zzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 927);

    auto g_zzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 928);

    auto g_zzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 929);

    auto g_zzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 930);

    auto g_zzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 931);

    auto g_zzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 932);

    auto g_zzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 933);

    auto g_zzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 934);

    auto g_zzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 935);

    auto g_zzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 936);

    auto g_zzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 937);

    auto g_zzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 938);

    auto g_zzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 939);

    auto g_zzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 940);

    auto g_zzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 941);

    auto g_zzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 942);

    auto g_zzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 943);

    auto g_zzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 944);

    /// Set up components of auxilary buffer : HSL

    auto g_xxxxx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl);

    auto g_xxxxx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 1);

    auto g_xxxxx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 2);

    auto g_xxxxx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 3);

    auto g_xxxxx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 4);

    auto g_xxxxx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 5);

    auto g_xxxxx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 6);

    auto g_xxxxx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 7);

    auto g_xxxxx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 8);

    auto g_xxxxx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 9);

    auto g_xxxxx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 10);

    auto g_xxxxx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 11);

    auto g_xxxxx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 12);

    auto g_xxxxx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 13);

    auto g_xxxxx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 14);

    auto g_xxxxx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 15);

    auto g_xxxxx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 16);

    auto g_xxxxx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 17);

    auto g_xxxxx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 18);

    auto g_xxxxx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 19);

    auto g_xxxxx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 20);

    auto g_xxxxx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 21);

    auto g_xxxxx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 22);

    auto g_xxxxx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 23);

    auto g_xxxxx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 24);

    auto g_xxxxx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 25);

    auto g_xxxxx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 26);

    auto g_xxxxx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 27);

    auto g_xxxxx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 28);

    auto g_xxxxx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 29);

    auto g_xxxxx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 30);

    auto g_xxxxx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 31);

    auto g_xxxxx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 32);

    auto g_xxxxx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 33);

    auto g_xxxxx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 34);

    auto g_xxxxx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 35);

    auto g_xxxxx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 36);

    auto g_xxxxx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 37);

    auto g_xxxxx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 38);

    auto g_xxxxx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 39);

    auto g_xxxxx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 40);

    auto g_xxxxx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 41);

    auto g_xxxxx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 42);

    auto g_xxxxx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 43);

    auto g_xxxxx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 44);

    auto g_xxxxy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 45);

    auto g_xxxxy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 47);

    auto g_xxxxy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 50);

    auto g_xxxxy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 54);

    auto g_xxxxy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 59);

    auto g_xxxxy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 65);

    auto g_xxxxy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 72);

    auto g_xxxxy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 80);

    auto g_xxxxz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 90);

    auto g_xxxxz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 91);

    auto g_xxxxz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 93);

    auto g_xxxxz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 96);

    auto g_xxxxz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 100);

    auto g_xxxxz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 105);

    auto g_xxxxz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 111);

    auto g_xxxxz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 118);

    auto g_xxxyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 135);

    auto g_xxxyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 136);

    auto g_xxxyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 137);

    auto g_xxxyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 138);

    auto g_xxxyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 139);

    auto g_xxxyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 140);

    auto g_xxxyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 141);

    auto g_xxxyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 142);

    auto g_xxxyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 143);

    auto g_xxxyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 144);

    auto g_xxxyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 145);

    auto g_xxxyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 146);

    auto g_xxxyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 147);

    auto g_xxxyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 148);

    auto g_xxxyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 149);

    auto g_xxxyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 150);

    auto g_xxxyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 151);

    auto g_xxxyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 152);

    auto g_xxxyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 153);

    auto g_xxxyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 154);

    auto g_xxxyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 155);

    auto g_xxxyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 156);

    auto g_xxxyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 157);

    auto g_xxxyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 158);

    auto g_xxxyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 159);

    auto g_xxxyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 160);

    auto g_xxxyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 161);

    auto g_xxxyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 162);

    auto g_xxxyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 163);

    auto g_xxxyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 164);

    auto g_xxxyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 165);

    auto g_xxxyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 166);

    auto g_xxxyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 167);

    auto g_xxxyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 168);

    auto g_xxxyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 169);

    auto g_xxxyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 170);

    auto g_xxxyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 171);

    auto g_xxxyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 172);

    auto g_xxxyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 173);

    auto g_xxxyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 174);

    auto g_xxxyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 175);

    auto g_xxxyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 176);

    auto g_xxxyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 177);

    auto g_xxxyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 178);

    auto g_xxxyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 179);

    auto g_xxxzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 225);

    auto g_xxxzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 226);

    auto g_xxxzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 227);

    auto g_xxxzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 228);

    auto g_xxxzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 229);

    auto g_xxxzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 230);

    auto g_xxxzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 231);

    auto g_xxxzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 232);

    auto g_xxxzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 233);

    auto g_xxxzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 234);

    auto g_xxxzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 235);

    auto g_xxxzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 236);

    auto g_xxxzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 237);

    auto g_xxxzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 238);

    auto g_xxxzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 239);

    auto g_xxxzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 240);

    auto g_xxxzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 241);

    auto g_xxxzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 242);

    auto g_xxxzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 243);

    auto g_xxxzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 244);

    auto g_xxxzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 245);

    auto g_xxxzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 246);

    auto g_xxxzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 247);

    auto g_xxxzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 248);

    auto g_xxxzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 249);

    auto g_xxxzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 250);

    auto g_xxxzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 251);

    auto g_xxxzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 252);

    auto g_xxxzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 253);

    auto g_xxxzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 254);

    auto g_xxxzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 255);

    auto g_xxxzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 256);

    auto g_xxxzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 257);

    auto g_xxxzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 258);

    auto g_xxxzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 259);

    auto g_xxxzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 260);

    auto g_xxxzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 261);

    auto g_xxxzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 262);

    auto g_xxxzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 263);

    auto g_xxxzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 264);

    auto g_xxxzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 265);

    auto g_xxxzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 266);

    auto g_xxxzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 267);

    auto g_xxxzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 268);

    auto g_xxxzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 269);

    auto g_xxyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 270);

    auto g_xxyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 271);

    auto g_xxyyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 272);

    auto g_xxyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 273);

    auto g_xxyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 274);

    auto g_xxyyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 275);

    auto g_xxyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 276);

    auto g_xxyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 277);

    auto g_xxyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 278);

    auto g_xxyyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 279);

    auto g_xxyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 280);

    auto g_xxyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 281);

    auto g_xxyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 282);

    auto g_xxyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 283);

    auto g_xxyyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 284);

    auto g_xxyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 285);

    auto g_xxyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 286);

    auto g_xxyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 287);

    auto g_xxyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 288);

    auto g_xxyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 289);

    auto g_xxyyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 290);

    auto g_xxyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 291);

    auto g_xxyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 292);

    auto g_xxyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 293);

    auto g_xxyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 294);

    auto g_xxyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 295);

    auto g_xxyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 296);

    auto g_xxyyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 297);

    auto g_xxyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 298);

    auto g_xxyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 299);

    auto g_xxyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 300);

    auto g_xxyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 301);

    auto g_xxyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 302);

    auto g_xxyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 303);

    auto g_xxyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 304);

    auto g_xxyyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 305);

    auto g_xxyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 306);

    auto g_xxyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 307);

    auto g_xxyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 308);

    auto g_xxyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 309);

    auto g_xxyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 310);

    auto g_xxyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 311);

    auto g_xxyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 312);

    auto g_xxyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 313);

    auto g_xxyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 314);

    auto g_xxyyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 316);

    auto g_xxyyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 318);

    auto g_xxyyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 321);

    auto g_xxyyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 325);

    auto g_xxyyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 330);

    auto g_xxyyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 336);

    auto g_xxyyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 343);

    auto g_xxyzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 360);

    auto g_xxyzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 362);

    auto g_xxyzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 365);

    auto g_xxyzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 369);

    auto g_xxyzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 374);

    auto g_xxyzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 380);

    auto g_xxyzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 387);

    auto g_xxyzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 395);

    auto g_xxzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 405);

    auto g_xxzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 406);

    auto g_xxzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 407);

    auto g_xxzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 408);

    auto g_xxzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 409);

    auto g_xxzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 410);

    auto g_xxzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 411);

    auto g_xxzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 412);

    auto g_xxzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 413);

    auto g_xxzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 414);

    auto g_xxzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 415);

    auto g_xxzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 416);

    auto g_xxzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 417);

    auto g_xxzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 418);

    auto g_xxzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 419);

    auto g_xxzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 420);

    auto g_xxzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 421);

    auto g_xxzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 422);

    auto g_xxzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 423);

    auto g_xxzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 424);

    auto g_xxzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 425);

    auto g_xxzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 426);

    auto g_xxzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 427);

    auto g_xxzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 428);

    auto g_xxzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 429);

    auto g_xxzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 430);

    auto g_xxzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 431);

    auto g_xxzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 432);

    auto g_xxzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 433);

    auto g_xxzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 434);

    auto g_xxzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 435);

    auto g_xxzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 436);

    auto g_xxzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 437);

    auto g_xxzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 438);

    auto g_xxzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 439);

    auto g_xxzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 440);

    auto g_xxzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 441);

    auto g_xxzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 442);

    auto g_xxzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 443);

    auto g_xxzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 444);

    auto g_xxzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 445);

    auto g_xxzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 446);

    auto g_xxzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 447);

    auto g_xxzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 448);

    auto g_xxzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 449);

    auto g_xyyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 451);

    auto g_xyyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 453);

    auto g_xyyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 454);

    auto g_xyyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 456);

    auto g_xyyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 457);

    auto g_xyyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 458);

    auto g_xyyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 460);

    auto g_xyyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 461);

    auto g_xyyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 462);

    auto g_xyyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 463);

    auto g_xyyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 465);

    auto g_xyyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 466);

    auto g_xyyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 467);

    auto g_xyyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 468);

    auto g_xyyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 469);

    auto g_xyyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 471);

    auto g_xyyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 472);

    auto g_xyyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 473);

    auto g_xyyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 474);

    auto g_xyyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 475);

    auto g_xyyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 476);

    auto g_xyyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 478);

    auto g_xyyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 479);

    auto g_xyyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 480);

    auto g_xyyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 481);

    auto g_xyyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 482);

    auto g_xyyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 483);

    auto g_xyyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 484);

    auto g_xyyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 486);

    auto g_xyyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 487);

    auto g_xyyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 488);

    auto g_xyyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 489);

    auto g_xyyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 490);

    auto g_xyyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 491);

    auto g_xyyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 492);

    auto g_xyyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 493);

    auto g_xyyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 494);

    auto g_xyyzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 544);

    auto g_xyyzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 547);

    auto g_xyyzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 548);

    auto g_xyyzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 551);

    auto g_xyyzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 552);

    auto g_xyyzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 553);

    auto g_xyyzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 556);

    auto g_xyyzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 557);

    auto g_xyyzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 558);

    auto g_xyyzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 559);

    auto g_xyyzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 562);

    auto g_xyyzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 563);

    auto g_xyyzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 564);

    auto g_xyyzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 565);

    auto g_xyyzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 566);

    auto g_xyyzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 569);

    auto g_xyyzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 570);

    auto g_xyyzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 571);

    auto g_xyyzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 572);

    auto g_xyyzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 573);

    auto g_xyyzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 574);

    auto g_xyyzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 576);

    auto g_xyyzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 577);

    auto g_xyyzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 578);

    auto g_xyyzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 579);

    auto g_xyyzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 580);

    auto g_xyyzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 581);

    auto g_xyyzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 582);

    auto g_xyyzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 583);

    auto g_xyyzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 584);

    auto g_xzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 632);

    auto g_xzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 634);

    auto g_xzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 635);

    auto g_xzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 637);

    auto g_xzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 638);

    auto g_xzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 639);

    auto g_xzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 641);

    auto g_xzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 642);

    auto g_xzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 643);

    auto g_xzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 644);

    auto g_xzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 646);

    auto g_xzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 647);

    auto g_xzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 648);

    auto g_xzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 649);

    auto g_xzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 650);

    auto g_xzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 652);

    auto g_xzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 653);

    auto g_xzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 654);

    auto g_xzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 655);

    auto g_xzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 656);

    auto g_xzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 657);

    auto g_xzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 659);

    auto g_xzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 660);

    auto g_xzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 661);

    auto g_xzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 662);

    auto g_xzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 663);

    auto g_xzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 664);

    auto g_xzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 665);

    auto g_xzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 666);

    auto g_xzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 667);

    auto g_xzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 668);

    auto g_xzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 669);

    auto g_xzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 670);

    auto g_xzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 671);

    auto g_xzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 672);

    auto g_xzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 673);

    auto g_xzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 674);

    auto g_yyyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 675);

    auto g_yyyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 676);

    auto g_yyyyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 677);

    auto g_yyyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 678);

    auto g_yyyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 679);

    auto g_yyyyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 680);

    auto g_yyyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 681);

    auto g_yyyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 682);

    auto g_yyyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 683);

    auto g_yyyyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 684);

    auto g_yyyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 685);

    auto g_yyyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 686);

    auto g_yyyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 687);

    auto g_yyyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 688);

    auto g_yyyyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 689);

    auto g_yyyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 690);

    auto g_yyyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 691);

    auto g_yyyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 692);

    auto g_yyyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 693);

    auto g_yyyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 694);

    auto g_yyyyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 695);

    auto g_yyyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 696);

    auto g_yyyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 697);

    auto g_yyyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 698);

    auto g_yyyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 699);

    auto g_yyyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 700);

    auto g_yyyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 701);

    auto g_yyyyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 702);

    auto g_yyyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 703);

    auto g_yyyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 704);

    auto g_yyyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 705);

    auto g_yyyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 706);

    auto g_yyyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 707);

    auto g_yyyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 708);

    auto g_yyyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 709);

    auto g_yyyyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 710);

    auto g_yyyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 711);

    auto g_yyyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 712);

    auto g_yyyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 713);

    auto g_yyyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 714);

    auto g_yyyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 715);

    auto g_yyyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 716);

    auto g_yyyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 717);

    auto g_yyyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 718);

    auto g_yyyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 719);

    auto g_yyyyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 721);

    auto g_yyyyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 723);

    auto g_yyyyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 726);

    auto g_yyyyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 730);

    auto g_yyyyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 735);

    auto g_yyyyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 741);

    auto g_yyyyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 748);

    auto g_yyyyz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 756);

    auto g_yyyzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 765);

    auto g_yyyzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 766);

    auto g_yyyzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 767);

    auto g_yyyzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 768);

    auto g_yyyzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 769);

    auto g_yyyzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 770);

    auto g_yyyzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 771);

    auto g_yyyzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 772);

    auto g_yyyzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 773);

    auto g_yyyzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 774);

    auto g_yyyzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 775);

    auto g_yyyzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 776);

    auto g_yyyzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 777);

    auto g_yyyzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 778);

    auto g_yyyzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 779);

    auto g_yyyzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 780);

    auto g_yyyzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 781);

    auto g_yyyzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 782);

    auto g_yyyzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 783);

    auto g_yyyzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 784);

    auto g_yyyzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 785);

    auto g_yyyzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 786);

    auto g_yyyzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 787);

    auto g_yyyzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 788);

    auto g_yyyzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 789);

    auto g_yyyzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 790);

    auto g_yyyzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 791);

    auto g_yyyzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 792);

    auto g_yyyzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 793);

    auto g_yyyzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 794);

    auto g_yyyzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 795);

    auto g_yyyzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 796);

    auto g_yyyzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 797);

    auto g_yyyzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 798);

    auto g_yyyzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 799);

    auto g_yyyzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 800);

    auto g_yyyzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 801);

    auto g_yyyzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 802);

    auto g_yyyzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 803);

    auto g_yyyzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 804);

    auto g_yyyzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 805);

    auto g_yyyzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 806);

    auto g_yyyzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 807);

    auto g_yyyzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 808);

    auto g_yyyzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 809);

    auto g_yyzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 810);

    auto g_yyzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 811);

    auto g_yyzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 812);

    auto g_yyzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 813);

    auto g_yyzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 814);

    auto g_yyzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 815);

    auto g_yyzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 816);

    auto g_yyzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 817);

    auto g_yyzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 818);

    auto g_yyzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 819);

    auto g_yyzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 820);

    auto g_yyzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 821);

    auto g_yyzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 822);

    auto g_yyzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 823);

    auto g_yyzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 824);

    auto g_yyzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 825);

    auto g_yyzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 826);

    auto g_yyzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 827);

    auto g_yyzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 828);

    auto g_yyzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 829);

    auto g_yyzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 830);

    auto g_yyzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 831);

    auto g_yyzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 832);

    auto g_yyzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 833);

    auto g_yyzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 834);

    auto g_yyzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 835);

    auto g_yyzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 836);

    auto g_yyzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 837);

    auto g_yyzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 838);

    auto g_yyzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 839);

    auto g_yyzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 840);

    auto g_yyzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 841);

    auto g_yyzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 842);

    auto g_yyzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 843);

    auto g_yyzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 844);

    auto g_yyzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 845);

    auto g_yyzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 846);

    auto g_yyzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 847);

    auto g_yyzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 848);

    auto g_yyzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 849);

    auto g_yyzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 850);

    auto g_yyzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 851);

    auto g_yyzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 852);

    auto g_yyzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 853);

    auto g_yyzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 854);

    auto g_yzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 855);

    auto g_yzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 857);

    auto g_yzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 859);

    auto g_yzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 860);

    auto g_yzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 862);

    auto g_yzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 863);

    auto g_yzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 864);

    auto g_yzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 866);

    auto g_yzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 867);

    auto g_yzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 868);

    auto g_yzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 869);

    auto g_yzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 871);

    auto g_yzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 872);

    auto g_yzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 873);

    auto g_yzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 874);

    auto g_yzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 875);

    auto g_yzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 877);

    auto g_yzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 878);

    auto g_yzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 879);

    auto g_yzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 880);

    auto g_yzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 881);

    auto g_yzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 882);

    auto g_yzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 884);

    auto g_yzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 885);

    auto g_yzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 886);

    auto g_yzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 887);

    auto g_yzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 888);

    auto g_yzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 889);

    auto g_yzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 890);

    auto g_yzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 892);

    auto g_yzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 893);

    auto g_yzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 894);

    auto g_yzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 895);

    auto g_yzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 896);

    auto g_yzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 897);

    auto g_yzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 898);

    auto g_yzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 899);

    auto g_zzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_hsl + 900);

    auto g_zzzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_hsl + 901);

    auto g_zzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_hsl + 902);

    auto g_zzzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_hsl + 903);

    auto g_zzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_hsl + 904);

    auto g_zzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_hsl + 905);

    auto g_zzzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_hsl + 906);

    auto g_zzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_hsl + 907);

    auto g_zzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_hsl + 908);

    auto g_zzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_hsl + 909);

    auto g_zzzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_hsl + 910);

    auto g_zzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_hsl + 911);

    auto g_zzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_hsl + 912);

    auto g_zzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_hsl + 913);

    auto g_zzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_hsl + 914);

    auto g_zzzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 915);

    auto g_zzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 916);

    auto g_zzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 917);

    auto g_zzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 918);

    auto g_zzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 919);

    auto g_zzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 920);

    auto g_zzzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 921);

    auto g_zzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 922);

    auto g_zzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 923);

    auto g_zzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 924);

    auto g_zzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 925);

    auto g_zzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 926);

    auto g_zzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 927);

    auto g_zzzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 928);

    auto g_zzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 929);

    auto g_zzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 930);

    auto g_zzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 931);

    auto g_zzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 932);

    auto g_zzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 933);

    auto g_zzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 934);

    auto g_zzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 935);

    auto g_zzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_hsl + 936);

    auto g_zzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_hsl + 937);

    auto g_zzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_hsl + 938);

    auto g_zzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_hsl + 939);

    auto g_zzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_hsl + 940);

    auto g_zzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 941);

    auto g_zzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 942);

    auto g_zzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 943);

    auto g_zzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_hsl + 944);

    /// Set up components of auxilary buffer : ISK

    auto g_xxxxxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk);

    auto g_xxxxxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 1);

    auto g_xxxxxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 2);

    auto g_xxxxxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 3);

    auto g_xxxxxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 4);

    auto g_xxxxxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 5);

    auto g_xxxxxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 6);

    auto g_xxxxxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 7);

    auto g_xxxxxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 8);

    auto g_xxxxxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 9);

    auto g_xxxxxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 10);

    auto g_xxxxxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 11);

    auto g_xxxxxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 12);

    auto g_xxxxxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 13);

    auto g_xxxxxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 14);

    auto g_xxxxxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 15);

    auto g_xxxxxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 16);

    auto g_xxxxxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 17);

    auto g_xxxxxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 18);

    auto g_xxxxxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 19);

    auto g_xxxxxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 20);

    auto g_xxxxxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 21);

    auto g_xxxxxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 22);

    auto g_xxxxxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 23);

    auto g_xxxxxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 24);

    auto g_xxxxxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 25);

    auto g_xxxxxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 26);

    auto g_xxxxxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 27);

    auto g_xxxxxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 28);

    auto g_xxxxxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 29);

    auto g_xxxxxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 30);

    auto g_xxxxxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 31);

    auto g_xxxxxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 32);

    auto g_xxxxxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 33);

    auto g_xxxxxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 34);

    auto g_xxxxxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 35);

    auto g_xxxxxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 74);

    auto g_xxxxxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 76);

    auto g_xxxxxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 77);

    auto g_xxxxxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 79);

    auto g_xxxxxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 80);

    auto g_xxxxxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 81);

    auto g_xxxxxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 83);

    auto g_xxxxxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 84);

    auto g_xxxxxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 85);

    auto g_xxxxxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 86);

    auto g_xxxxxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 88);

    auto g_xxxxxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 89);

    auto g_xxxxxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 90);

    auto g_xxxxxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 91);

    auto g_xxxxxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 92);

    auto g_xxxxxz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 94);

    auto g_xxxxxz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 95);

    auto g_xxxxxz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 96);

    auto g_xxxxxz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 97);

    auto g_xxxxxz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 98);

    auto g_xxxxxz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 99);

    auto g_xxxxxz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 101);

    auto g_xxxxxz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 102);

    auto g_xxxxxz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 103);

    auto g_xxxxxz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 104);

    auto g_xxxxxz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 105);

    auto g_xxxxxz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 106);

    auto g_xxxxxz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 107);

    auto g_xxxxyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 108);

    auto g_xxxxyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 109);

    auto g_xxxxyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 110);

    auto g_xxxxyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 111);

    auto g_xxxxyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 112);

    auto g_xxxxyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 113);

    auto g_xxxxyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 114);

    auto g_xxxxyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 115);

    auto g_xxxxyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 116);

    auto g_xxxxyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 117);

    auto g_xxxxyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 118);

    auto g_xxxxyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 119);

    auto g_xxxxyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 120);

    auto g_xxxxyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 121);

    auto g_xxxxyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 122);

    auto g_xxxxyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 123);

    auto g_xxxxyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 124);

    auto g_xxxxyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 125);

    auto g_xxxxyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 126);

    auto g_xxxxyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 127);

    auto g_xxxxyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 128);

    auto g_xxxxyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 129);

    auto g_xxxxyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 130);

    auto g_xxxxyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 131);

    auto g_xxxxyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 132);

    auto g_xxxxyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 133);

    auto g_xxxxyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 134);

    auto g_xxxxyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 135);

    auto g_xxxxyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 136);

    auto g_xxxxyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 137);

    auto g_xxxxyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 138);

    auto g_xxxxyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 139);

    auto g_xxxxyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 140);

    auto g_xxxxyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 141);

    auto g_xxxxyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 142);

    auto g_xxxxyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 143);

    auto g_xxxxzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 180);

    auto g_xxxxzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 181);

    auto g_xxxxzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 182);

    auto g_xxxxzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 183);

    auto g_xxxxzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 184);

    auto g_xxxxzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 185);

    auto g_xxxxzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 186);

    auto g_xxxxzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 187);

    auto g_xxxxzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 188);

    auto g_xxxxzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 189);

    auto g_xxxxzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 190);

    auto g_xxxxzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 191);

    auto g_xxxxzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 192);

    auto g_xxxxzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 193);

    auto g_xxxxzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 194);

    auto g_xxxxzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 195);

    auto g_xxxxzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 196);

    auto g_xxxxzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 197);

    auto g_xxxxzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 198);

    auto g_xxxxzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 199);

    auto g_xxxxzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 200);

    auto g_xxxxzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 201);

    auto g_xxxxzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 202);

    auto g_xxxxzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 203);

    auto g_xxxxzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 204);

    auto g_xxxxzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 205);

    auto g_xxxxzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 206);

    auto g_xxxxzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 207);

    auto g_xxxxzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 208);

    auto g_xxxxzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 209);

    auto g_xxxxzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 210);

    auto g_xxxxzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 211);

    auto g_xxxxzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 212);

    auto g_xxxxzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 213);

    auto g_xxxxzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 214);

    auto g_xxxxzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 215);

    auto g_xxxyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 216);

    auto g_xxxyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 217);

    auto g_xxxyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 218);

    auto g_xxxyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 219);

    auto g_xxxyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 220);

    auto g_xxxyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 221);

    auto g_xxxyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 222);

    auto g_xxxyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 223);

    auto g_xxxyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 224);

    auto g_xxxyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 225);

    auto g_xxxyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 226);

    auto g_xxxyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 227);

    auto g_xxxyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 228);

    auto g_xxxyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 229);

    auto g_xxxyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 230);

    auto g_xxxyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 231);

    auto g_xxxyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 232);

    auto g_xxxyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 233);

    auto g_xxxyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 234);

    auto g_xxxyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 235);

    auto g_xxxyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 236);

    auto g_xxxyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 237);

    auto g_xxxyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 238);

    auto g_xxxyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 239);

    auto g_xxxyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 240);

    auto g_xxxyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 241);

    auto g_xxxyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 242);

    auto g_xxxyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 243);

    auto g_xxxyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 244);

    auto g_xxxyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 245);

    auto g_xxxyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 246);

    auto g_xxxyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 247);

    auto g_xxxyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 248);

    auto g_xxxyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 249);

    auto g_xxxyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 250);

    auto g_xxxyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 251);

    auto g_xxxzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 324);

    auto g_xxxzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 325);

    auto g_xxxzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 326);

    auto g_xxxzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 327);

    auto g_xxxzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 328);

    auto g_xxxzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 329);

    auto g_xxxzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 330);

    auto g_xxxzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 331);

    auto g_xxxzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 332);

    auto g_xxxzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 333);

    auto g_xxxzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 334);

    auto g_xxxzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 335);

    auto g_xxxzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 336);

    auto g_xxxzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 337);

    auto g_xxxzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 338);

    auto g_xxxzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 339);

    auto g_xxxzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 340);

    auto g_xxxzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 341);

    auto g_xxxzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 342);

    auto g_xxxzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 343);

    auto g_xxxzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 344);

    auto g_xxxzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 345);

    auto g_xxxzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 346);

    auto g_xxxzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 347);

    auto g_xxxzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 348);

    auto g_xxxzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 349);

    auto g_xxxzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 350);

    auto g_xxxzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 351);

    auto g_xxxzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 352);

    auto g_xxxzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 353);

    auto g_xxxzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 354);

    auto g_xxxzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 355);

    auto g_xxxzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 356);

    auto g_xxxzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 357);

    auto g_xxxzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 358);

    auto g_xxxzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 359);

    auto g_xxyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 360);

    auto g_xxyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 361);

    auto g_xxyyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 362);

    auto g_xxyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 363);

    auto g_xxyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 364);

    auto g_xxyyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 365);

    auto g_xxyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 366);

    auto g_xxyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 367);

    auto g_xxyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 368);

    auto g_xxyyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 369);

    auto g_xxyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 370);

    auto g_xxyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 371);

    auto g_xxyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 372);

    auto g_xxyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 373);

    auto g_xxyyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 374);

    auto g_xxyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 375);

    auto g_xxyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 376);

    auto g_xxyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 377);

    auto g_xxyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 378);

    auto g_xxyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 379);

    auto g_xxyyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 380);

    auto g_xxyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 381);

    auto g_xxyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 382);

    auto g_xxyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 383);

    auto g_xxyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 384);

    auto g_xxyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 385);

    auto g_xxyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 386);

    auto g_xxyyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 387);

    auto g_xxyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 388);

    auto g_xxyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 389);

    auto g_xxyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 390);

    auto g_xxyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 391);

    auto g_xxyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 392);

    auto g_xxyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 393);

    auto g_xxyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 394);

    auto g_xxyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 395);

    auto g_xxyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 436);

    auto g_xxyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 439);

    auto g_xxyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 440);

    auto g_xxyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 443);

    auto g_xxyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 444);

    auto g_xxyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 445);

    auto g_xxyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 448);

    auto g_xxyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 449);

    auto g_xxyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 450);

    auto g_xxyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 451);

    auto g_xxyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 454);

    auto g_xxyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 455);

    auto g_xxyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 456);

    auto g_xxyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 457);

    auto g_xxyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 458);

    auto g_xxyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 461);

    auto g_xxyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 462);

    auto g_xxyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 463);

    auto g_xxyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 464);

    auto g_xxyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 465);

    auto g_xxyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 466);

    auto g_xxzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 504);

    auto g_xxzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 505);

    auto g_xxzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 506);

    auto g_xxzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 507);

    auto g_xxzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 508);

    auto g_xxzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 509);

    auto g_xxzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 510);

    auto g_xxzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 511);

    auto g_xxzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 512);

    auto g_xxzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 513);

    auto g_xxzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 514);

    auto g_xxzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 515);

    auto g_xxzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 516);

    auto g_xxzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 517);

    auto g_xxzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 518);

    auto g_xxzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 519);

    auto g_xxzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 520);

    auto g_xxzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 521);

    auto g_xxzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 522);

    auto g_xxzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 523);

    auto g_xxzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 524);

    auto g_xxzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 525);

    auto g_xxzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 526);

    auto g_xxzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 527);

    auto g_xxzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 528);

    auto g_xxzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 529);

    auto g_xxzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 530);

    auto g_xxzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 531);

    auto g_xxzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 532);

    auto g_xxzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 533);

    auto g_xxzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 534);

    auto g_xxzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 535);

    auto g_xxzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 536);

    auto g_xxzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 537);

    auto g_xxzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 538);

    auto g_xxzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 539);

    auto g_xyyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 541);

    auto g_xyyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 543);

    auto g_xyyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 544);

    auto g_xyyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 546);

    auto g_xyyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 547);

    auto g_xyyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 548);

    auto g_xyyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 550);

    auto g_xyyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 551);

    auto g_xyyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 552);

    auto g_xyyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 553);

    auto g_xyyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 555);

    auto g_xyyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 556);

    auto g_xyyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 557);

    auto g_xyyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 558);

    auto g_xyyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 559);

    auto g_xyyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 561);

    auto g_xyyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 562);

    auto g_xyyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 563);

    auto g_xyyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 564);

    auto g_xyyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 565);

    auto g_xyyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 566);

    auto g_xyyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 568);

    auto g_xyyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 569);

    auto g_xyyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 570);

    auto g_xyyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 571);

    auto g_xyyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 572);

    auto g_xyyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 573);

    auto g_xyyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 574);

    auto g_xyyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 616);

    auto g_xyyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 619);

    auto g_xyyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 620);

    auto g_xyyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 623);

    auto g_xyyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 624);

    auto g_xyyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 625);

    auto g_xyyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 628);

    auto g_xyyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 629);

    auto g_xyyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 630);

    auto g_xyyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 631);

    auto g_xyyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 634);

    auto g_xyyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 635);

    auto g_xyyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 636);

    auto g_xyyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 637);

    auto g_xyyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 638);

    auto g_xyyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 641);

    auto g_xyyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 642);

    auto g_xyyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 643);

    auto g_xyyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 644);

    auto g_xyyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 645);

    auto g_xyyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 646);

    auto g_xyyzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 652);

    auto g_xyyzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 655);

    auto g_xyyzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 656);

    auto g_xyyzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 659);

    auto g_xyyzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 660);

    auto g_xyyzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 661);

    auto g_xyyzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 664);

    auto g_xyyzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 665);

    auto g_xyyzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 666);

    auto g_xyyzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 667);

    auto g_xyyzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 670);

    auto g_xyyzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 671);

    auto g_xyyzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 672);

    auto g_xyyzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 673);

    auto g_xyyzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 674);

    auto g_xyyzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 677);

    auto g_xyyzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 678);

    auto g_xyyzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 679);

    auto g_xyyzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 680);

    auto g_xyyzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 681);

    auto g_xyyzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 682);

    auto g_xzzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 722);

    auto g_xzzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 724);

    auto g_xzzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 725);

    auto g_xzzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 727);

    auto g_xzzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 728);

    auto g_xzzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 729);

    auto g_xzzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 731);

    auto g_xzzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 732);

    auto g_xzzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 733);

    auto g_xzzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 734);

    auto g_xzzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 736);

    auto g_xzzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 737);

    auto g_xzzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 738);

    auto g_xzzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 739);

    auto g_xzzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 740);

    auto g_xzzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 742);

    auto g_xzzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 743);

    auto g_xzzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 744);

    auto g_xzzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 745);

    auto g_xzzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 746);

    auto g_xzzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 747);

    auto g_xzzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 749);

    auto g_xzzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 750);

    auto g_xzzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 751);

    auto g_xzzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 752);

    auto g_xzzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 753);

    auto g_xzzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 754);

    auto g_xzzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 755);

    auto g_yyyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 756);

    auto g_yyyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 757);

    auto g_yyyyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 758);

    auto g_yyyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 759);

    auto g_yyyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 760);

    auto g_yyyyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 761);

    auto g_yyyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 762);

    auto g_yyyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 763);

    auto g_yyyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 764);

    auto g_yyyyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 765);

    auto g_yyyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 766);

    auto g_yyyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 767);

    auto g_yyyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 768);

    auto g_yyyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 769);

    auto g_yyyyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 770);

    auto g_yyyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 771);

    auto g_yyyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 772);

    auto g_yyyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 773);

    auto g_yyyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 774);

    auto g_yyyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 775);

    auto g_yyyyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 776);

    auto g_yyyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 777);

    auto g_yyyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 778);

    auto g_yyyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 779);

    auto g_yyyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 780);

    auto g_yyyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 781);

    auto g_yyyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 782);

    auto g_yyyyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 783);

    auto g_yyyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 784);

    auto g_yyyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 785);

    auto g_yyyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 786);

    auto g_yyyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 787);

    auto g_yyyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 788);

    auto g_yyyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 789);

    auto g_yyyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 790);

    auto g_yyyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 791);

    auto g_yyyyyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 794);

    auto g_yyyyyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 796);

    auto g_yyyyyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 797);

    auto g_yyyyyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 799);

    auto g_yyyyyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 800);

    auto g_yyyyyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 801);

    auto g_yyyyyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 803);

    auto g_yyyyyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 804);

    auto g_yyyyyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 805);

    auto g_yyyyyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 806);

    auto g_yyyyyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 808);

    auto g_yyyyyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 809);

    auto g_yyyyyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 810);

    auto g_yyyyyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 811);

    auto g_yyyyyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 812);

    auto g_yyyyyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 814);

    auto g_yyyyyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 815);

    auto g_yyyyyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 816);

    auto g_yyyyyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 817);

    auto g_yyyyyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 818);

    auto g_yyyyyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 819);

    auto g_yyyyyz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 821);

    auto g_yyyyyz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 822);

    auto g_yyyyyz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 823);

    auto g_yyyyyz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 824);

    auto g_yyyyyz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 825);

    auto g_yyyyyz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 826);

    auto g_yyyyyz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 827);

    auto g_yyyyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 828);

    auto g_yyyyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 829);

    auto g_yyyyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 830);

    auto g_yyyyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 831);

    auto g_yyyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 832);

    auto g_yyyyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 833);

    auto g_yyyyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 834);

    auto g_yyyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 835);

    auto g_yyyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 836);

    auto g_yyyyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 837);

    auto g_yyyyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 838);

    auto g_yyyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 839);

    auto g_yyyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 840);

    auto g_yyyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 841);

    auto g_yyyyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 842);

    auto g_yyyyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 843);

    auto g_yyyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 844);

    auto g_yyyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 845);

    auto g_yyyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 846);

    auto g_yyyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 847);

    auto g_yyyyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 848);

    auto g_yyyyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 849);

    auto g_yyyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 850);

    auto g_yyyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 851);

    auto g_yyyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 852);

    auto g_yyyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 853);

    auto g_yyyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 854);

    auto g_yyyyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 855);

    auto g_yyyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 856);

    auto g_yyyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 857);

    auto g_yyyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 858);

    auto g_yyyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 859);

    auto g_yyyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 860);

    auto g_yyyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 861);

    auto g_yyyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 862);

    auto g_yyyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 863);

    auto g_yyyzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 864);

    auto g_yyyzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 865);

    auto g_yyyzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 866);

    auto g_yyyzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 867);

    auto g_yyyzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 868);

    auto g_yyyzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 869);

    auto g_yyyzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 870);

    auto g_yyyzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 871);

    auto g_yyyzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 872);

    auto g_yyyzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 873);

    auto g_yyyzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 874);

    auto g_yyyzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 875);

    auto g_yyyzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 876);

    auto g_yyyzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 877);

    auto g_yyyzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 878);

    auto g_yyyzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 879);

    auto g_yyyzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 880);

    auto g_yyyzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 881);

    auto g_yyyzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 882);

    auto g_yyyzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 883);

    auto g_yyyzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 884);

    auto g_yyyzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 885);

    auto g_yyyzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 886);

    auto g_yyyzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 887);

    auto g_yyyzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 888);

    auto g_yyyzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 889);

    auto g_yyyzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 890);

    auto g_yyyzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 891);

    auto g_yyyzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 892);

    auto g_yyyzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 893);

    auto g_yyyzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 894);

    auto g_yyyzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 895);

    auto g_yyyzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 896);

    auto g_yyyzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 897);

    auto g_yyyzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 898);

    auto g_yyyzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 899);

    auto g_yyzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 900);

    auto g_yyzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 901);

    auto g_yyzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 902);

    auto g_yyzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 903);

    auto g_yyzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 904);

    auto g_yyzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 905);

    auto g_yyzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 906);

    auto g_yyzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 907);

    auto g_yyzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 908);

    auto g_yyzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 909);

    auto g_yyzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 910);

    auto g_yyzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 911);

    auto g_yyzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 912);

    auto g_yyzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 913);

    auto g_yyzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 914);

    auto g_yyzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 915);

    auto g_yyzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 916);

    auto g_yyzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 917);

    auto g_yyzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 918);

    auto g_yyzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 919);

    auto g_yyzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 920);

    auto g_yyzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 921);

    auto g_yyzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 922);

    auto g_yyzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 923);

    auto g_yyzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 924);

    auto g_yyzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 925);

    auto g_yyzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 926);

    auto g_yyzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 927);

    auto g_yyzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 928);

    auto g_yyzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 929);

    auto g_yyzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 930);

    auto g_yyzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 931);

    auto g_yyzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 932);

    auto g_yyzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 933);

    auto g_yyzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 934);

    auto g_yyzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 935);

    auto g_yzzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 937);

    auto g_yzzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 938);

    auto g_yzzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 939);

    auto g_yzzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 940);

    auto g_yzzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 941);

    auto g_yzzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 942);

    auto g_yzzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 943);

    auto g_yzzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 944);

    auto g_yzzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 945);

    auto g_yzzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 946);

    auto g_yzzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 947);

    auto g_yzzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 948);

    auto g_yzzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 949);

    auto g_yzzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 950);

    auto g_yzzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 951);

    auto g_yzzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 952);

    auto g_yzzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 953);

    auto g_yzzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 954);

    auto g_yzzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 955);

    auto g_yzzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 956);

    auto g_yzzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 957);

    auto g_yzzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 958);

    auto g_yzzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 959);

    auto g_yzzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 960);

    auto g_yzzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 961);

    auto g_yzzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 962);

    auto g_yzzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 963);

    auto g_yzzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 964);

    auto g_yzzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 965);

    auto g_yzzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 966);

    auto g_yzzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 967);

    auto g_yzzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 968);

    auto g_yzzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 969);

    auto g_yzzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 970);

    auto g_yzzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 971);

    auto g_zzzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 972);

    auto g_zzzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 973);

    auto g_zzzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 974);

    auto g_zzzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 975);

    auto g_zzzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 976);

    auto g_zzzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 977);

    auto g_zzzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 978);

    auto g_zzzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 979);

    auto g_zzzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 980);

    auto g_zzzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 981);

    auto g_zzzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 982);

    auto g_zzzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 983);

    auto g_zzzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 984);

    auto g_zzzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 985);

    auto g_zzzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 986);

    auto g_zzzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 987);

    auto g_zzzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 988);

    auto g_zzzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 989);

    auto g_zzzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 990);

    auto g_zzzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 991);

    auto g_zzzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 992);

    auto g_zzzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 993);

    auto g_zzzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 994);

    auto g_zzzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 995);

    auto g_zzzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 996);

    auto g_zzzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 997);

    auto g_zzzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 998);

    auto g_zzzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 999);

    auto g_zzzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 1000);

    auto g_zzzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 1001);

    auto g_zzzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 1002);

    auto g_zzzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 1003);

    auto g_zzzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 1004);

    auto g_zzzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 1005);

    auto g_zzzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 1006);

    auto g_zzzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 1007);

    /// Set up components of auxilary buffer : ISL

    auto g_xxxxxx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl);

    auto g_xxxxxx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 1);

    auto g_xxxxxx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 2);

    auto g_xxxxxx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 3);

    auto g_xxxxxx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 4);

    auto g_xxxxxx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 5);

    auto g_xxxxxx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 6);

    auto g_xxxxxx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 7);

    auto g_xxxxxx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 8);

    auto g_xxxxxx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 9);

    auto g_xxxxxx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 10);

    auto g_xxxxxx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 11);

    auto g_xxxxxx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 12);

    auto g_xxxxxx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 13);

    auto g_xxxxxx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 14);

    auto g_xxxxxx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 15);

    auto g_xxxxxx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 16);

    auto g_xxxxxx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 17);

    auto g_xxxxxx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 18);

    auto g_xxxxxx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 19);

    auto g_xxxxxx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 20);

    auto g_xxxxxx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 21);

    auto g_xxxxxx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 22);

    auto g_xxxxxx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 23);

    auto g_xxxxxx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 24);

    auto g_xxxxxx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 25);

    auto g_xxxxxx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 26);

    auto g_xxxxxx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 27);

    auto g_xxxxxx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 28);

    auto g_xxxxxx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 29);

    auto g_xxxxxx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 30);

    auto g_xxxxxx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 31);

    auto g_xxxxxx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 32);

    auto g_xxxxxx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 33);

    auto g_xxxxxx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 34);

    auto g_xxxxxx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 35);

    auto g_xxxxxx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 36);

    auto g_xxxxxx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 37);

    auto g_xxxxxx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 38);

    auto g_xxxxxx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 39);

    auto g_xxxxxx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 40);

    auto g_xxxxxx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 41);

    auto g_xxxxxx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 42);

    auto g_xxxxxx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 43);

    auto g_xxxxxx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 44);

    auto g_xxxxxy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 45);

    auto g_xxxxxy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 46);

    auto g_xxxxxy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 47);

    auto g_xxxxxy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 48);

    auto g_xxxxxy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 50);

    auto g_xxxxxy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 51);

    auto g_xxxxxy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 54);

    auto g_xxxxxy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 55);

    auto g_xxxxxy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 59);

    auto g_xxxxxy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 60);

    auto g_xxxxxy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 65);

    auto g_xxxxxy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 66);

    auto g_xxxxxy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 72);

    auto g_xxxxxy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 73);

    auto g_xxxxxy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 80);

    auto g_xxxxxy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 81);

    auto g_xxxxxz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 90);

    auto g_xxxxxz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 91);

    auto g_xxxxxz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 92);

    auto g_xxxxxz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 93);

    auto g_xxxxxz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 94);

    auto g_xxxxxz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 95);

    auto g_xxxxxz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 96);

    auto g_xxxxxz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 97);

    auto g_xxxxxz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 98);

    auto g_xxxxxz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 99);

    auto g_xxxxxz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 100);

    auto g_xxxxxz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 101);

    auto g_xxxxxz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 102);

    auto g_xxxxxz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 103);

    auto g_xxxxxz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 104);

    auto g_xxxxxz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 105);

    auto g_xxxxxz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 106);

    auto g_xxxxxz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 107);

    auto g_xxxxxz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 108);

    auto g_xxxxxz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 109);

    auto g_xxxxxz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 110);

    auto g_xxxxxz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 111);

    auto g_xxxxxz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 112);

    auto g_xxxxxz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 113);

    auto g_xxxxxz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 114);

    auto g_xxxxxz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 115);

    auto g_xxxxxz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 116);

    auto g_xxxxxz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 117);

    auto g_xxxxxz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 118);

    auto g_xxxxxz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 119);

    auto g_xxxxxz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 120);

    auto g_xxxxxz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 121);

    auto g_xxxxxz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 122);

    auto g_xxxxxz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 123);

    auto g_xxxxxz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 124);

    auto g_xxxxxz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 125);

    auto g_xxxxxz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 127);

    auto g_xxxxxz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 128);

    auto g_xxxxxz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 129);

    auto g_xxxxxz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 130);

    auto g_xxxxxz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 131);

    auto g_xxxxxz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 132);

    auto g_xxxxxz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 133);

    auto g_xxxxxz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 134);

    auto g_xxxxyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 135);

    auto g_xxxxyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 136);

    auto g_xxxxyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 137);

    auto g_xxxxyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 138);

    auto g_xxxxyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 139);

    auto g_xxxxyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 140);

    auto g_xxxxyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 141);

    auto g_xxxxyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 142);

    auto g_xxxxyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 143);

    auto g_xxxxyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 144);

    auto g_xxxxyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 145);

    auto g_xxxxyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 146);

    auto g_xxxxyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 147);

    auto g_xxxxyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 148);

    auto g_xxxxyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 149);

    auto g_xxxxyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 150);

    auto g_xxxxyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 151);

    auto g_xxxxyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 152);

    auto g_xxxxyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 153);

    auto g_xxxxyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 154);

    auto g_xxxxyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 155);

    auto g_xxxxyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 156);

    auto g_xxxxyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 157);

    auto g_xxxxyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 158);

    auto g_xxxxyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 159);

    auto g_xxxxyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 160);

    auto g_xxxxyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 161);

    auto g_xxxxyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 162);

    auto g_xxxxyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 163);

    auto g_xxxxyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 164);

    auto g_xxxxyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 165);

    auto g_xxxxyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 166);

    auto g_xxxxyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 167);

    auto g_xxxxyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 168);

    auto g_xxxxyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 169);

    auto g_xxxxyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 170);

    auto g_xxxxyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 171);

    auto g_xxxxyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 172);

    auto g_xxxxyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 173);

    auto g_xxxxyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 174);

    auto g_xxxxyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 175);

    auto g_xxxxyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 176);

    auto g_xxxxyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 177);

    auto g_xxxxyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 178);

    auto g_xxxxyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 179);

    auto g_xxxxzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 225);

    auto g_xxxxzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 226);

    auto g_xxxxzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 227);

    auto g_xxxxzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 228);

    auto g_xxxxzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 229);

    auto g_xxxxzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 230);

    auto g_xxxxzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 231);

    auto g_xxxxzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 232);

    auto g_xxxxzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 233);

    auto g_xxxxzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 234);

    auto g_xxxxzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 235);

    auto g_xxxxzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 236);

    auto g_xxxxzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 237);

    auto g_xxxxzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 238);

    auto g_xxxxzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 239);

    auto g_xxxxzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 240);

    auto g_xxxxzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 241);

    auto g_xxxxzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 242);

    auto g_xxxxzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 243);

    auto g_xxxxzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 244);

    auto g_xxxxzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 245);

    auto g_xxxxzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 246);

    auto g_xxxxzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 247);

    auto g_xxxxzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 248);

    auto g_xxxxzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 249);

    auto g_xxxxzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 250);

    auto g_xxxxzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 251);

    auto g_xxxxzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 252);

    auto g_xxxxzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 253);

    auto g_xxxxzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 254);

    auto g_xxxxzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 255);

    auto g_xxxxzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 256);

    auto g_xxxxzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 257);

    auto g_xxxxzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 258);

    auto g_xxxxzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 259);

    auto g_xxxxzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 260);

    auto g_xxxxzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 261);

    auto g_xxxxzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 262);

    auto g_xxxxzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 263);

    auto g_xxxxzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 264);

    auto g_xxxxzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 265);

    auto g_xxxxzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 266);

    auto g_xxxxzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 267);

    auto g_xxxxzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 268);

    auto g_xxxxzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 269);

    auto g_xxxyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 270);

    auto g_xxxyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 271);

    auto g_xxxyyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 272);

    auto g_xxxyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 273);

    auto g_xxxyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 274);

    auto g_xxxyyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 275);

    auto g_xxxyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 276);

    auto g_xxxyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 277);

    auto g_xxxyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 278);

    auto g_xxxyyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 279);

    auto g_xxxyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 280);

    auto g_xxxyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 281);

    auto g_xxxyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 282);

    auto g_xxxyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 283);

    auto g_xxxyyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 284);

    auto g_xxxyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 285);

    auto g_xxxyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 286);

    auto g_xxxyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 287);

    auto g_xxxyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 288);

    auto g_xxxyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 289);

    auto g_xxxyyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 290);

    auto g_xxxyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 291);

    auto g_xxxyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 292);

    auto g_xxxyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 293);

    auto g_xxxyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 294);

    auto g_xxxyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 295);

    auto g_xxxyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 296);

    auto g_xxxyyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 297);

    auto g_xxxyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 298);

    auto g_xxxyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 299);

    auto g_xxxyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 300);

    auto g_xxxyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 301);

    auto g_xxxyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 302);

    auto g_xxxyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 303);

    auto g_xxxyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 304);

    auto g_xxxyyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 305);

    auto g_xxxyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 306);

    auto g_xxxyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 307);

    auto g_xxxyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 308);

    auto g_xxxyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 309);

    auto g_xxxyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 310);

    auto g_xxxyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 311);

    auto g_xxxyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 312);

    auto g_xxxyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 313);

    auto g_xxxyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 314);

    auto g_xxxyyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 316);

    auto g_xxxyyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 318);

    auto g_xxxyyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 321);

    auto g_xxxyyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 325);

    auto g_xxxyyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 330);

    auto g_xxxyyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 336);

    auto g_xxxyyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 343);

    auto g_xxxyzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 360);

    auto g_xxxyzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 362);

    auto g_xxxyzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 365);

    auto g_xxxyzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 369);

    auto g_xxxyzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 374);

    auto g_xxxyzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 380);

    auto g_xxxyzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 387);

    auto g_xxxyzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 395);

    auto g_xxxzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 405);

    auto g_xxxzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 406);

    auto g_xxxzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 407);

    auto g_xxxzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 408);

    auto g_xxxzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 409);

    auto g_xxxzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 410);

    auto g_xxxzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 411);

    auto g_xxxzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 412);

    auto g_xxxzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 413);

    auto g_xxxzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 414);

    auto g_xxxzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 415);

    auto g_xxxzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 416);

    auto g_xxxzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 417);

    auto g_xxxzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 418);

    auto g_xxxzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 419);

    auto g_xxxzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 420);

    auto g_xxxzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 421);

    auto g_xxxzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 422);

    auto g_xxxzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 423);

    auto g_xxxzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 424);

    auto g_xxxzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 425);

    auto g_xxxzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 426);

    auto g_xxxzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 427);

    auto g_xxxzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 428);

    auto g_xxxzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 429);

    auto g_xxxzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 430);

    auto g_xxxzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 431);

    auto g_xxxzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 432);

    auto g_xxxzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 433);

    auto g_xxxzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 434);

    auto g_xxxzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 435);

    auto g_xxxzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 436);

    auto g_xxxzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 437);

    auto g_xxxzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 438);

    auto g_xxxzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 439);

    auto g_xxxzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 440);

    auto g_xxxzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 441);

    auto g_xxxzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 442);

    auto g_xxxzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 443);

    auto g_xxxzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 444);

    auto g_xxxzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 445);

    auto g_xxxzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 446);

    auto g_xxxzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 447);

    auto g_xxxzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 448);

    auto g_xxxzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 449);

    auto g_xxyyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 450);

    auto g_xxyyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 451);

    auto g_xxyyyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 452);

    auto g_xxyyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 453);

    auto g_xxyyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 454);

    auto g_xxyyyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 455);

    auto g_xxyyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 456);

    auto g_xxyyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 457);

    auto g_xxyyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 458);

    auto g_xxyyyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 459);

    auto g_xxyyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 460);

    auto g_xxyyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 461);

    auto g_xxyyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 462);

    auto g_xxyyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 463);

    auto g_xxyyyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 464);

    auto g_xxyyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 465);

    auto g_xxyyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 466);

    auto g_xxyyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 467);

    auto g_xxyyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 468);

    auto g_xxyyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 469);

    auto g_xxyyyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 470);

    auto g_xxyyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 471);

    auto g_xxyyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 472);

    auto g_xxyyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 473);

    auto g_xxyyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 474);

    auto g_xxyyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 475);

    auto g_xxyyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 476);

    auto g_xxyyyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 477);

    auto g_xxyyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 478);

    auto g_xxyyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 479);

    auto g_xxyyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 480);

    auto g_xxyyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 481);

    auto g_xxyyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 482);

    auto g_xxyyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 483);

    auto g_xxyyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 484);

    auto g_xxyyyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 485);

    auto g_xxyyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 486);

    auto g_xxyyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 487);

    auto g_xxyyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 488);

    auto g_xxyyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 489);

    auto g_xxyyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 490);

    auto g_xxyyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 491);

    auto g_xxyyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 492);

    auto g_xxyyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 493);

    auto g_xxyyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 494);

    auto g_xxyyyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 496);

    auto g_xxyyyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 498);

    auto g_xxyyyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 501);

    auto g_xxyyyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 505);

    auto g_xxyyyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 510);

    auto g_xxyyyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 516);

    auto g_xxyyyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 523);

    auto g_xxyyzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 540);

    auto g_xxyyzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 541);

    auto g_xxyyzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 542);

    auto g_xxyyzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 543);

    auto g_xxyyzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 544);

    auto g_xxyyzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 545);

    auto g_xxyyzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 546);

    auto g_xxyyzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 547);

    auto g_xxyyzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 548);

    auto g_xxyyzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 549);

    auto g_xxyyzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 550);

    auto g_xxyyzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 551);

    auto g_xxyyzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 552);

    auto g_xxyyzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 553);

    auto g_xxyyzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 554);

    auto g_xxyyzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 555);

    auto g_xxyyzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 556);

    auto g_xxyyzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 557);

    auto g_xxyyzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 558);

    auto g_xxyyzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 559);

    auto g_xxyyzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 560);

    auto g_xxyyzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 561);

    auto g_xxyyzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 562);

    auto g_xxyyzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 563);

    auto g_xxyyzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 564);

    auto g_xxyyzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 565);

    auto g_xxyyzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 566);

    auto g_xxyyzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 567);

    auto g_xxyyzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 568);

    auto g_xxyyzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 569);

    auto g_xxyyzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 570);

    auto g_xxyyzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 571);

    auto g_xxyyzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 572);

    auto g_xxyyzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 573);

    auto g_xxyyzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 574);

    auto g_xxyyzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 575);

    auto g_xxyyzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 576);

    auto g_xxyyzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 577);

    auto g_xxyyzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 578);

    auto g_xxyyzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 579);

    auto g_xxyyzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 580);

    auto g_xxyyzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 581);

    auto g_xxyyzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 582);

    auto g_xxyyzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 583);

    auto g_xxyyzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 584);

    auto g_xxyzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 585);

    auto g_xxyzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 587);

    auto g_xxyzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 590);

    auto g_xxyzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 594);

    auto g_xxyzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 599);

    auto g_xxyzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 605);

    auto g_xxyzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 612);

    auto g_xxyzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 620);

    auto g_xxzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 630);

    auto g_xxzzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 631);

    auto g_xxzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 632);

    auto g_xxzzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 633);

    auto g_xxzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 634);

    auto g_xxzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 635);

    auto g_xxzzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 636);

    auto g_xxzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 637);

    auto g_xxzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 638);

    auto g_xxzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 639);

    auto g_xxzzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 640);

    auto g_xxzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 641);

    auto g_xxzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 642);

    auto g_xxzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 643);

    auto g_xxzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 644);

    auto g_xxzzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 645);

    auto g_xxzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 646);

    auto g_xxzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 647);

    auto g_xxzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 648);

    auto g_xxzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 649);

    auto g_xxzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 650);

    auto g_xxzzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 651);

    auto g_xxzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 652);

    auto g_xxzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 653);

    auto g_xxzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 654);

    auto g_xxzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 655);

    auto g_xxzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 656);

    auto g_xxzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 657);

    auto g_xxzzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 658);

    auto g_xxzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 659);

    auto g_xxzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 660);

    auto g_xxzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 661);

    auto g_xxzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 662);

    auto g_xxzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 663);

    auto g_xxzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 664);

    auto g_xxzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 665);

    auto g_xxzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 666);

    auto g_xxzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 667);

    auto g_xxzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 668);

    auto g_xxzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 669);

    auto g_xxzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 670);

    auto g_xxzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 671);

    auto g_xxzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 672);

    auto g_xxzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 673);

    auto g_xxzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 674);

    auto g_xyyyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 675);

    auto g_xyyyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 676);

    auto g_xyyyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 678);

    auto g_xyyyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 679);

    auto g_xyyyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 681);

    auto g_xyyyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 682);

    auto g_xyyyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 683);

    auto g_xyyyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 685);

    auto g_xyyyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 686);

    auto g_xyyyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 687);

    auto g_xyyyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 688);

    auto g_xyyyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 690);

    auto g_xyyyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 691);

    auto g_xyyyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 692);

    auto g_xyyyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 693);

    auto g_xyyyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 694);

    auto g_xyyyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 696);

    auto g_xyyyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 697);

    auto g_xyyyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 698);

    auto g_xyyyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 699);

    auto g_xyyyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 700);

    auto g_xyyyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 701);

    auto g_xyyyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 703);

    auto g_xyyyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 704);

    auto g_xyyyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 705);

    auto g_xyyyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 706);

    auto g_xyyyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 707);

    auto g_xyyyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 708);

    auto g_xyyyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 709);

    auto g_xyyyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 711);

    auto g_xyyyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 712);

    auto g_xyyyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 713);

    auto g_xyyyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 714);

    auto g_xyyyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 715);

    auto g_xyyyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 716);

    auto g_xyyyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 717);

    auto g_xyyyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 718);

    auto g_xyyyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 719);

    auto g_xyyyzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 769);

    auto g_xyyyzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 772);

    auto g_xyyyzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 773);

    auto g_xyyyzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 776);

    auto g_xyyyzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 777);

    auto g_xyyyzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 778);

    auto g_xyyyzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 781);

    auto g_xyyyzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 782);

    auto g_xyyyzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 783);

    auto g_xyyyzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 784);

    auto g_xyyyzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 787);

    auto g_xyyyzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 788);

    auto g_xyyyzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 789);

    auto g_xyyyzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 790);

    auto g_xyyyzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 791);

    auto g_xyyyzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 794);

    auto g_xyyyzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 795);

    auto g_xyyyzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 796);

    auto g_xyyyzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 797);

    auto g_xyyyzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 798);

    auto g_xyyyzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 799);

    auto g_xyyyzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 801);

    auto g_xyyyzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 802);

    auto g_xyyyzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 803);

    auto g_xyyyzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 804);

    auto g_xyyyzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 805);

    auto g_xyyyzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 806);

    auto g_xyyyzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 807);

    auto g_xyyyzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 808);

    auto g_xyyyzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 809);

    auto g_xyyzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 814);

    auto g_xyyzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 817);

    auto g_xyyzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 818);

    auto g_xyyzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 821);

    auto g_xyyzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 822);

    auto g_xyyzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 823);

    auto g_xyyzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 826);

    auto g_xyyzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 827);

    auto g_xyyzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 828);

    auto g_xyyzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 829);

    auto g_xyyzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 832);

    auto g_xyyzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 833);

    auto g_xyyzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 834);

    auto g_xyyzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 835);

    auto g_xyyzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 836);

    auto g_xyyzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 839);

    auto g_xyyzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 840);

    auto g_xyyzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 841);

    auto g_xyyzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 842);

    auto g_xyyzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 843);

    auto g_xyyzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 844);

    auto g_xyyzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 846);

    auto g_xyyzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 847);

    auto g_xyyzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 848);

    auto g_xyyzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 849);

    auto g_xyyzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 850);

    auto g_xyyzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 851);

    auto g_xyyzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 852);

    auto g_xyyzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 853);

    auto g_xyyzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 854);

    auto g_xzzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 900);

    auto g_xzzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 902);

    auto g_xzzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 904);

    auto g_xzzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 905);

    auto g_xzzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 907);

    auto g_xzzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 908);

    auto g_xzzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 909);

    auto g_xzzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 911);

    auto g_xzzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 912);

    auto g_xzzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 913);

    auto g_xzzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 914);

    auto g_xzzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 916);

    auto g_xzzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 917);

    auto g_xzzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 918);

    auto g_xzzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 919);

    auto g_xzzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 920);

    auto g_xzzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 922);

    auto g_xzzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 923);

    auto g_xzzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 924);

    auto g_xzzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 925);

    auto g_xzzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 926);

    auto g_xzzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 927);

    auto g_xzzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 929);

    auto g_xzzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 930);

    auto g_xzzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 931);

    auto g_xzzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 932);

    auto g_xzzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 933);

    auto g_xzzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 934);

    auto g_xzzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 935);

    auto g_xzzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 936);

    auto g_xzzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 937);

    auto g_xzzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 938);

    auto g_xzzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 939);

    auto g_xzzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 940);

    auto g_xzzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 941);

    auto g_xzzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 942);

    auto g_xzzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 943);

    auto g_xzzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 944);

    auto g_yyyyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 945);

    auto g_yyyyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 946);

    auto g_yyyyyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 947);

    auto g_yyyyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 948);

    auto g_yyyyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 949);

    auto g_yyyyyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 950);

    auto g_yyyyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 951);

    auto g_yyyyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 952);

    auto g_yyyyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 953);

    auto g_yyyyyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 954);

    auto g_yyyyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 955);

    auto g_yyyyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 956);

    auto g_yyyyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 957);

    auto g_yyyyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 958);

    auto g_yyyyyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 959);

    auto g_yyyyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 960);

    auto g_yyyyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 961);

    auto g_yyyyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 962);

    auto g_yyyyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 963);

    auto g_yyyyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 964);

    auto g_yyyyyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 965);

    auto g_yyyyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 966);

    auto g_yyyyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 967);

    auto g_yyyyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 968);

    auto g_yyyyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 969);

    auto g_yyyyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 970);

    auto g_yyyyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 971);

    auto g_yyyyyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 972);

    auto g_yyyyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 973);

    auto g_yyyyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 974);

    auto g_yyyyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 975);

    auto g_yyyyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 976);

    auto g_yyyyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 977);

    auto g_yyyyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 978);

    auto g_yyyyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 979);

    auto g_yyyyyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 980);

    auto g_yyyyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 981);

    auto g_yyyyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 982);

    auto g_yyyyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 983);

    auto g_yyyyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 984);

    auto g_yyyyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 985);

    auto g_yyyyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 986);

    auto g_yyyyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 987);

    auto g_yyyyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 988);

    auto g_yyyyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 989);

    auto g_yyyyyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 991);

    auto g_yyyyyz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 992);

    auto g_yyyyyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 993);

    auto g_yyyyyz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 994);

    auto g_yyyyyz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 995);

    auto g_yyyyyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 996);

    auto g_yyyyyz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 997);

    auto g_yyyyyz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 998);

    auto g_yyyyyz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 999);

    auto g_yyyyyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 1000);

    auto g_yyyyyz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 1001);

    auto g_yyyyyz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 1002);

    auto g_yyyyyz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 1003);

    auto g_yyyyyz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 1004);

    auto g_yyyyyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1005);

    auto g_yyyyyz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1006);

    auto g_yyyyyz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1007);

    auto g_yyyyyz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1008);

    auto g_yyyyyz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1009);

    auto g_yyyyyz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1010);

    auto g_yyyyyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1011);

    auto g_yyyyyz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1012);

    auto g_yyyyyz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1013);

    auto g_yyyyyz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1014);

    auto g_yyyyyz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1015);

    auto g_yyyyyz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1016);

    auto g_yyyyyz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1017);

    auto g_yyyyyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1018);

    auto g_yyyyyz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1019);

    auto g_yyyyyz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1020);

    auto g_yyyyyz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1021);

    auto g_yyyyyz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1022);

    auto g_yyyyyz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1023);

    auto g_yyyyyz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1024);

    auto g_yyyyyz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1025);

    auto g_yyyyyz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1026);

    auto g_yyyyyz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1027);

    auto g_yyyyyz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1028);

    auto g_yyyyyz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1029);

    auto g_yyyyyz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1030);

    auto g_yyyyyz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1031);

    auto g_yyyyyz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1032);

    auto g_yyyyyz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1033);

    auto g_yyyyyz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1034);

    auto g_yyyyzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 1035);

    auto g_yyyyzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 1036);

    auto g_yyyyzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 1037);

    auto g_yyyyzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 1038);

    auto g_yyyyzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 1039);

    auto g_yyyyzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 1040);

    auto g_yyyyzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 1041);

    auto g_yyyyzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 1042);

    auto g_yyyyzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 1043);

    auto g_yyyyzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 1044);

    auto g_yyyyzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 1045);

    auto g_yyyyzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 1046);

    auto g_yyyyzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 1047);

    auto g_yyyyzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 1048);

    auto g_yyyyzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 1049);

    auto g_yyyyzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1050);

    auto g_yyyyzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1051);

    auto g_yyyyzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1052);

    auto g_yyyyzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1053);

    auto g_yyyyzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1054);

    auto g_yyyyzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1055);

    auto g_yyyyzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1056);

    auto g_yyyyzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1057);

    auto g_yyyyzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1058);

    auto g_yyyyzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1059);

    auto g_yyyyzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1060);

    auto g_yyyyzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1061);

    auto g_yyyyzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1062);

    auto g_yyyyzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1063);

    auto g_yyyyzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1064);

    auto g_yyyyzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1065);

    auto g_yyyyzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1066);

    auto g_yyyyzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1067);

    auto g_yyyyzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1068);

    auto g_yyyyzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1069);

    auto g_yyyyzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1070);

    auto g_yyyyzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1071);

    auto g_yyyyzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1072);

    auto g_yyyyzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1073);

    auto g_yyyyzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1074);

    auto g_yyyyzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1075);

    auto g_yyyyzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1076);

    auto g_yyyyzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1077);

    auto g_yyyyzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1078);

    auto g_yyyyzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1079);

    auto g_yyyzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 1080);

    auto g_yyyzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 1081);

    auto g_yyyzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 1082);

    auto g_yyyzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 1083);

    auto g_yyyzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 1084);

    auto g_yyyzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 1085);

    auto g_yyyzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 1086);

    auto g_yyyzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 1087);

    auto g_yyyzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 1088);

    auto g_yyyzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 1089);

    auto g_yyyzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 1090);

    auto g_yyyzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 1091);

    auto g_yyyzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 1092);

    auto g_yyyzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 1093);

    auto g_yyyzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 1094);

    auto g_yyyzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1095);

    auto g_yyyzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1096);

    auto g_yyyzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1097);

    auto g_yyyzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1098);

    auto g_yyyzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1099);

    auto g_yyyzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1100);

    auto g_yyyzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1101);

    auto g_yyyzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1102);

    auto g_yyyzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1103);

    auto g_yyyzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1104);

    auto g_yyyzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1105);

    auto g_yyyzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1106);

    auto g_yyyzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1107);

    auto g_yyyzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1108);

    auto g_yyyzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1109);

    auto g_yyyzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1110);

    auto g_yyyzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1111);

    auto g_yyyzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1112);

    auto g_yyyzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1113);

    auto g_yyyzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1114);

    auto g_yyyzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1115);

    auto g_yyyzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1116);

    auto g_yyyzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1117);

    auto g_yyyzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1118);

    auto g_yyyzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1119);

    auto g_yyyzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1120);

    auto g_yyyzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1121);

    auto g_yyyzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1122);

    auto g_yyyzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1123);

    auto g_yyyzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1124);

    auto g_yyzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 1125);

    auto g_yyzzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 1126);

    auto g_yyzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 1127);

    auto g_yyzzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 1128);

    auto g_yyzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 1129);

    auto g_yyzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 1130);

    auto g_yyzzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 1131);

    auto g_yyzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 1132);

    auto g_yyzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 1133);

    auto g_yyzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 1134);

    auto g_yyzzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 1135);

    auto g_yyzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 1136);

    auto g_yyzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 1137);

    auto g_yyzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 1138);

    auto g_yyzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 1139);

    auto g_yyzzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1140);

    auto g_yyzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1141);

    auto g_yyzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1142);

    auto g_yyzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1143);

    auto g_yyzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1144);

    auto g_yyzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1145);

    auto g_yyzzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1146);

    auto g_yyzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1147);

    auto g_yyzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1148);

    auto g_yyzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1149);

    auto g_yyzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1150);

    auto g_yyzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1151);

    auto g_yyzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1152);

    auto g_yyzzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1153);

    auto g_yyzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1154);

    auto g_yyzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1155);

    auto g_yyzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1156);

    auto g_yyzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1157);

    auto g_yyzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1158);

    auto g_yyzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1159);

    auto g_yyzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1160);

    auto g_yyzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1161);

    auto g_yyzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1162);

    auto g_yyzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1163);

    auto g_yyzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1164);

    auto g_yyzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1165);

    auto g_yyzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1166);

    auto g_yyzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1167);

    auto g_yyzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1168);

    auto g_yyzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1169);

    auto g_yzzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 1170);

    auto g_yzzzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 1171);

    auto g_yzzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 1172);

    auto g_yzzzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 1173);

    auto g_yzzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 1174);

    auto g_yzzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 1175);

    auto g_yzzzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 1176);

    auto g_yzzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 1177);

    auto g_yzzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 1178);

    auto g_yzzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 1179);

    auto g_yzzzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 1180);

    auto g_yzzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 1181);

    auto g_yzzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 1182);

    auto g_yzzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 1183);

    auto g_yzzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 1184);

    auto g_yzzzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1185);

    auto g_yzzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1186);

    auto g_yzzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1187);

    auto g_yzzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1188);

    auto g_yzzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1189);

    auto g_yzzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1190);

    auto g_yzzzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1191);

    auto g_yzzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1192);

    auto g_yzzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1193);

    auto g_yzzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1194);

    auto g_yzzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1195);

    auto g_yzzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1196);

    auto g_yzzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1197);

    auto g_yzzzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1198);

    auto g_yzzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1199);

    auto g_yzzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1200);

    auto g_yzzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1201);

    auto g_yzzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1202);

    auto g_yzzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1203);

    auto g_yzzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1204);

    auto g_yzzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1205);

    auto g_yzzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1206);

    auto g_yzzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1207);

    auto g_yzzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1208);

    auto g_yzzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1209);

    auto g_yzzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1210);

    auto g_yzzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1211);

    auto g_yzzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1212);

    auto g_yzzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1213);

    auto g_yzzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1214);

    auto g_zzzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_isl + 1215);

    auto g_zzzzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_isl + 1216);

    auto g_zzzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_isl + 1217);

    auto g_zzzzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_isl + 1218);

    auto g_zzzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_isl + 1219);

    auto g_zzzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_isl + 1220);

    auto g_zzzzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_isl + 1221);

    auto g_zzzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_isl + 1222);

    auto g_zzzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_isl + 1223);

    auto g_zzzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_isl + 1224);

    auto g_zzzzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_isl + 1225);

    auto g_zzzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_isl + 1226);

    auto g_zzzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_isl + 1227);

    auto g_zzzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_isl + 1228);

    auto g_zzzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_isl + 1229);

    auto g_zzzzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1230);

    auto g_zzzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1231);

    auto g_zzzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1232);

    auto g_zzzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1233);

    auto g_zzzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1234);

    auto g_zzzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1235);

    auto g_zzzzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1236);

    auto g_zzzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1237);

    auto g_zzzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1238);

    auto g_zzzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1239);

    auto g_zzzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1240);

    auto g_zzzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1241);

    auto g_zzzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1242);

    auto g_zzzzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1243);

    auto g_zzzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1244);

    auto g_zzzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1245);

    auto g_zzzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1246);

    auto g_zzzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1247);

    auto g_zzzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1248);

    auto g_zzzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1249);

    auto g_zzzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1250);

    auto g_zzzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_isl + 1251);

    auto g_zzzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_isl + 1252);

    auto g_zzzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_isl + 1253);

    auto g_zzzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_isl + 1254);

    auto g_zzzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_isl + 1255);

    auto g_zzzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1256);

    auto g_zzzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1257);

    auto g_zzzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1258);

    auto g_zzzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_isl + 1259);

    /// Set up 0-45 components of targeted buffer : KSL

    auto g_xxxxxxx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl);

    auto g_xxxxxxx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1);

    auto g_xxxxxxx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 2);

    auto g_xxxxxxx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 3);

    auto g_xxxxxxx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 4);

    auto g_xxxxxxx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 5);

    auto g_xxxxxxx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 6);

    auto g_xxxxxxx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 7);

    auto g_xxxxxxx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 8);

    auto g_xxxxxxx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 9);

    auto g_xxxxxxx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 10);

    auto g_xxxxxxx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 11);

    auto g_xxxxxxx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 12);

    auto g_xxxxxxx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 13);

    auto g_xxxxxxx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 14);

    auto g_xxxxxxx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 15);

    auto g_xxxxxxx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 16);

    auto g_xxxxxxx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 17);

    auto g_xxxxxxx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 18);

    auto g_xxxxxxx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 19);

    auto g_xxxxxxx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 20);

    auto g_xxxxxxx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 21);

    auto g_xxxxxxx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 22);

    auto g_xxxxxxx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 23);

    auto g_xxxxxxx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 24);

    auto g_xxxxxxx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 25);

    auto g_xxxxxxx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 26);

    auto g_xxxxxxx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 27);

    auto g_xxxxxxx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 28);

    auto g_xxxxxxx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 29);

    auto g_xxxxxxx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 30);

    auto g_xxxxxxx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 31);

    auto g_xxxxxxx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 32);

    auto g_xxxxxxx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 33);

    auto g_xxxxxxx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 34);

    auto g_xxxxxxx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 35);

    auto g_xxxxxxx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 36);

    auto g_xxxxxxx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 37);

    auto g_xxxxxxx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 38);

    auto g_xxxxxxx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 39);

    auto g_xxxxxxx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 40);

    auto g_xxxxxxx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 41);

    auto g_xxxxxxx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 42);

    auto g_xxxxxxx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 43);

    auto g_xxxxxxx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 44);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxxxx_0, g_xxxxx_0_xxxxxxxx_1, g_xxxxx_0_xxxxxxxy_0, g_xxxxx_0_xxxxxxxy_1, g_xxxxx_0_xxxxxxxz_0, g_xxxxx_0_xxxxxxxz_1, g_xxxxx_0_xxxxxxyy_0, g_xxxxx_0_xxxxxxyy_1, g_xxxxx_0_xxxxxxyz_0, g_xxxxx_0_xxxxxxyz_1, g_xxxxx_0_xxxxxxzz_0, g_xxxxx_0_xxxxxxzz_1, g_xxxxx_0_xxxxxyyy_0, g_xxxxx_0_xxxxxyyy_1, g_xxxxx_0_xxxxxyyz_0, g_xxxxx_0_xxxxxyyz_1, g_xxxxx_0_xxxxxyzz_0, g_xxxxx_0_xxxxxyzz_1, g_xxxxx_0_xxxxxzzz_0, g_xxxxx_0_xxxxxzzz_1, g_xxxxx_0_xxxxyyyy_0, g_xxxxx_0_xxxxyyyy_1, g_xxxxx_0_xxxxyyyz_0, g_xxxxx_0_xxxxyyyz_1, g_xxxxx_0_xxxxyyzz_0, g_xxxxx_0_xxxxyyzz_1, g_xxxxx_0_xxxxyzzz_0, g_xxxxx_0_xxxxyzzz_1, g_xxxxx_0_xxxxzzzz_0, g_xxxxx_0_xxxxzzzz_1, g_xxxxx_0_xxxyyyyy_0, g_xxxxx_0_xxxyyyyy_1, g_xxxxx_0_xxxyyyyz_0, g_xxxxx_0_xxxyyyyz_1, g_xxxxx_0_xxxyyyzz_0, g_xxxxx_0_xxxyyyzz_1, g_xxxxx_0_xxxyyzzz_0, g_xxxxx_0_xxxyyzzz_1, g_xxxxx_0_xxxyzzzz_0, g_xxxxx_0_xxxyzzzz_1, g_xxxxx_0_xxxzzzzz_0, g_xxxxx_0_xxxzzzzz_1, g_xxxxx_0_xxyyyyyy_0, g_xxxxx_0_xxyyyyyy_1, g_xxxxx_0_xxyyyyyz_0, g_xxxxx_0_xxyyyyyz_1, g_xxxxx_0_xxyyyyzz_0, g_xxxxx_0_xxyyyyzz_1, g_xxxxx_0_xxyyyzzz_0, g_xxxxx_0_xxyyyzzz_1, g_xxxxx_0_xxyyzzzz_0, g_xxxxx_0_xxyyzzzz_1, g_xxxxx_0_xxyzzzzz_0, g_xxxxx_0_xxyzzzzz_1, g_xxxxx_0_xxzzzzzz_0, g_xxxxx_0_xxzzzzzz_1, g_xxxxx_0_xyyyyyyy_0, g_xxxxx_0_xyyyyyyy_1, g_xxxxx_0_xyyyyyyz_0, g_xxxxx_0_xyyyyyyz_1, g_xxxxx_0_xyyyyyzz_0, g_xxxxx_0_xyyyyyzz_1, g_xxxxx_0_xyyyyzzz_0, g_xxxxx_0_xyyyyzzz_1, g_xxxxx_0_xyyyzzzz_0, g_xxxxx_0_xyyyzzzz_1, g_xxxxx_0_xyyzzzzz_0, g_xxxxx_0_xyyzzzzz_1, g_xxxxx_0_xyzzzzzz_0, g_xxxxx_0_xyzzzzzz_1, g_xxxxx_0_xzzzzzzz_0, g_xxxxx_0_xzzzzzzz_1, g_xxxxx_0_yyyyyyyy_0, g_xxxxx_0_yyyyyyyy_1, g_xxxxx_0_yyyyyyyz_0, g_xxxxx_0_yyyyyyyz_1, g_xxxxx_0_yyyyyyzz_0, g_xxxxx_0_yyyyyyzz_1, g_xxxxx_0_yyyyyzzz_0, g_xxxxx_0_yyyyyzzz_1, g_xxxxx_0_yyyyzzzz_0, g_xxxxx_0_yyyyzzzz_1, g_xxxxx_0_yyyzzzzz_0, g_xxxxx_0_yyyzzzzz_1, g_xxxxx_0_yyzzzzzz_0, g_xxxxx_0_yyzzzzzz_1, g_xxxxx_0_yzzzzzzz_0, g_xxxxx_0_yzzzzzzz_1, g_xxxxx_0_zzzzzzzz_0, g_xxxxx_0_zzzzzzzz_1, g_xxxxxx_0_xxxxxxx_1, g_xxxxxx_0_xxxxxxxx_1, g_xxxxxx_0_xxxxxxxy_1, g_xxxxxx_0_xxxxxxxz_1, g_xxxxxx_0_xxxxxxy_1, g_xxxxxx_0_xxxxxxyy_1, g_xxxxxx_0_xxxxxxyz_1, g_xxxxxx_0_xxxxxxz_1, g_xxxxxx_0_xxxxxxzz_1, g_xxxxxx_0_xxxxxyy_1, g_xxxxxx_0_xxxxxyyy_1, g_xxxxxx_0_xxxxxyyz_1, g_xxxxxx_0_xxxxxyz_1, g_xxxxxx_0_xxxxxyzz_1, g_xxxxxx_0_xxxxxzz_1, g_xxxxxx_0_xxxxxzzz_1, g_xxxxxx_0_xxxxyyy_1, g_xxxxxx_0_xxxxyyyy_1, g_xxxxxx_0_xxxxyyyz_1, g_xxxxxx_0_xxxxyyz_1, g_xxxxxx_0_xxxxyyzz_1, g_xxxxxx_0_xxxxyzz_1, g_xxxxxx_0_xxxxyzzz_1, g_xxxxxx_0_xxxxzzz_1, g_xxxxxx_0_xxxxzzzz_1, g_xxxxxx_0_xxxyyyy_1, g_xxxxxx_0_xxxyyyyy_1, g_xxxxxx_0_xxxyyyyz_1, g_xxxxxx_0_xxxyyyz_1, g_xxxxxx_0_xxxyyyzz_1, g_xxxxxx_0_xxxyyzz_1, g_xxxxxx_0_xxxyyzzz_1, g_xxxxxx_0_xxxyzzz_1, g_xxxxxx_0_xxxyzzzz_1, g_xxxxxx_0_xxxzzzz_1, g_xxxxxx_0_xxxzzzzz_1, g_xxxxxx_0_xxyyyyy_1, g_xxxxxx_0_xxyyyyyy_1, g_xxxxxx_0_xxyyyyyz_1, g_xxxxxx_0_xxyyyyz_1, g_xxxxxx_0_xxyyyyzz_1, g_xxxxxx_0_xxyyyzz_1, g_xxxxxx_0_xxyyyzzz_1, g_xxxxxx_0_xxyyzzz_1, g_xxxxxx_0_xxyyzzzz_1, g_xxxxxx_0_xxyzzzz_1, g_xxxxxx_0_xxyzzzzz_1, g_xxxxxx_0_xxzzzzz_1, g_xxxxxx_0_xxzzzzzz_1, g_xxxxxx_0_xyyyyyy_1, g_xxxxxx_0_xyyyyyyy_1, g_xxxxxx_0_xyyyyyyz_1, g_xxxxxx_0_xyyyyyz_1, g_xxxxxx_0_xyyyyyzz_1, g_xxxxxx_0_xyyyyzz_1, g_xxxxxx_0_xyyyyzzz_1, g_xxxxxx_0_xyyyzzz_1, g_xxxxxx_0_xyyyzzzz_1, g_xxxxxx_0_xyyzzzz_1, g_xxxxxx_0_xyyzzzzz_1, g_xxxxxx_0_xyzzzzz_1, g_xxxxxx_0_xyzzzzzz_1, g_xxxxxx_0_xzzzzzz_1, g_xxxxxx_0_xzzzzzzz_1, g_xxxxxx_0_yyyyyyy_1, g_xxxxxx_0_yyyyyyyy_1, g_xxxxxx_0_yyyyyyyz_1, g_xxxxxx_0_yyyyyyz_1, g_xxxxxx_0_yyyyyyzz_1, g_xxxxxx_0_yyyyyzz_1, g_xxxxxx_0_yyyyyzzz_1, g_xxxxxx_0_yyyyzzz_1, g_xxxxxx_0_yyyyzzzz_1, g_xxxxxx_0_yyyzzzz_1, g_xxxxxx_0_yyyzzzzz_1, g_xxxxxx_0_yyzzzzz_1, g_xxxxxx_0_yyzzzzzz_1, g_xxxxxx_0_yzzzzzz_1, g_xxxxxx_0_yzzzzzzz_1, g_xxxxxx_0_zzzzzzz_1, g_xxxxxx_0_zzzzzzzz_1, g_xxxxxxx_0_xxxxxxxx_0, g_xxxxxxx_0_xxxxxxxy_0, g_xxxxxxx_0_xxxxxxxz_0, g_xxxxxxx_0_xxxxxxyy_0, g_xxxxxxx_0_xxxxxxyz_0, g_xxxxxxx_0_xxxxxxzz_0, g_xxxxxxx_0_xxxxxyyy_0, g_xxxxxxx_0_xxxxxyyz_0, g_xxxxxxx_0_xxxxxyzz_0, g_xxxxxxx_0_xxxxxzzz_0, g_xxxxxxx_0_xxxxyyyy_0, g_xxxxxxx_0_xxxxyyyz_0, g_xxxxxxx_0_xxxxyyzz_0, g_xxxxxxx_0_xxxxyzzz_0, g_xxxxxxx_0_xxxxzzzz_0, g_xxxxxxx_0_xxxyyyyy_0, g_xxxxxxx_0_xxxyyyyz_0, g_xxxxxxx_0_xxxyyyzz_0, g_xxxxxxx_0_xxxyyzzz_0, g_xxxxxxx_0_xxxyzzzz_0, g_xxxxxxx_0_xxxzzzzz_0, g_xxxxxxx_0_xxyyyyyy_0, g_xxxxxxx_0_xxyyyyyz_0, g_xxxxxxx_0_xxyyyyzz_0, g_xxxxxxx_0_xxyyyzzz_0, g_xxxxxxx_0_xxyyzzzz_0, g_xxxxxxx_0_xxyzzzzz_0, g_xxxxxxx_0_xxzzzzzz_0, g_xxxxxxx_0_xyyyyyyy_0, g_xxxxxxx_0_xyyyyyyz_0, g_xxxxxxx_0_xyyyyyzz_0, g_xxxxxxx_0_xyyyyzzz_0, g_xxxxxxx_0_xyyyzzzz_0, g_xxxxxxx_0_xyyzzzzz_0, g_xxxxxxx_0_xyzzzzzz_0, g_xxxxxxx_0_xzzzzzzz_0, g_xxxxxxx_0_yyyyyyyy_0, g_xxxxxxx_0_yyyyyyyz_0, g_xxxxxxx_0_yyyyyyzz_0, g_xxxxxxx_0_yyyyyzzz_0, g_xxxxxxx_0_yyyyzzzz_0, g_xxxxxxx_0_yyyzzzzz_0, g_xxxxxxx_0_yyzzzzzz_0, g_xxxxxxx_0_yzzzzzzz_0, g_xxxxxxx_0_zzzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxxx_0_xxxxxxxx_0[i] = 6.0 * g_xxxxx_0_xxxxxxxx_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxxxx_1[i] * fz_be_0 + 8.0 * g_xxxxxx_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxxx_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxxxy_0[i] = 6.0 * g_xxxxx_0_xxxxxxxy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xxxxxx_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxxxz_0[i] = 6.0 * g_xxxxx_0_xxxxxxxz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xxxxxx_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxxyy_0[i] = 6.0 * g_xxxxx_0_xxxxxxyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xxxxxx_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxxyz_0[i] = 6.0 * g_xxxxx_0_xxxxxxyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxxxxx_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxxzz_0[i] = 6.0 * g_xxxxx_0_xxxxxxzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xxxxxx_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxyyy_0[i] = 6.0 * g_xxxxx_0_xxxxxyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xxxxxx_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxyyz_0[i] = 6.0 * g_xxxxx_0_xxxxxyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxxxxx_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxyzz_0[i] = 6.0 * g_xxxxx_0_xxxxxyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxxxxx_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxzzz_0[i] = 6.0 * g_xxxxx_0_xxxxxzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xxxxxx_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxyyyy_0[i] = 6.0 * g_xxxxx_0_xxxxyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxyyyz_0[i] = 6.0 * g_xxxxx_0_xxxxyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxyyzz_0[i] = 6.0 * g_xxxxx_0_xxxxyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxyzzz_0[i] = 6.0 * g_xxxxx_0_xxxxyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxzzzz_0[i] = 6.0 * g_xxxxx_0_xxxxzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyyyyy_0[i] = 6.0 * g_xxxxx_0_xxxyyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyyyyz_0[i] = 6.0 * g_xxxxx_0_xxxyyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyyyzz_0[i] = 6.0 * g_xxxxx_0_xxxyyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyyzzz_0[i] = 6.0 * g_xxxxx_0_xxxyyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyzzzz_0[i] = 6.0 * g_xxxxx_0_xxxyzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxzzzzz_0[i] = 6.0 * g_xxxxx_0_xxxzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyyyyy_0[i] = 6.0 * g_xxxxx_0_xxyyyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyyyyz_0[i] = 6.0 * g_xxxxx_0_xxyyyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyyyzz_0[i] = 6.0 * g_xxxxx_0_xxyyyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyyzzz_0[i] = 6.0 * g_xxxxx_0_xxyyyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyzzzz_0[i] = 6.0 * g_xxxxx_0_xxyyzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyzzzzz_0[i] = 6.0 * g_xxxxx_0_xxyzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxzzzzzz_0[i] = 6.0 * g_xxxxx_0_xxzzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyyyyy_0[i] = 6.0 * g_xxxxx_0_xyyyyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyyyyz_0[i] = 6.0 * g_xxxxx_0_xyyyyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyyyyz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyyyzz_0[i] = 6.0 * g_xxxxx_0_xyyyyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyyyzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyyzzz_0[i] = 6.0 * g_xxxxx_0_xyyyyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyyzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyzzzz_0[i] = 6.0 * g_xxxxx_0_xyyyzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyzzzzz_0[i] = 6.0 * g_xxxxx_0_xyyzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyzzzzzz_0[i] = 6.0 * g_xxxxx_0_xyzzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyzzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xzzzzzzz_0[i] = 6.0 * g_xxxxx_0_xzzzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyyyyy_0[i] = 6.0 * g_xxxxx_0_yyyyyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyyyyy_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyyyyz_0[i] = 6.0 * g_xxxxx_0_yyyyyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyyyyz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyyyzz_0[i] = 6.0 * g_xxxxx_0_yyyyyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyyyzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyyzzz_0[i] = 6.0 * g_xxxxx_0_yyyyyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyyzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyzzzz_0[i] = 6.0 * g_xxxxx_0_yyyyzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyzzzzz_0[i] = 6.0 * g_xxxxx_0_yyyzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyzzzzzz_0[i] = 6.0 * g_xxxxx_0_yyzzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyzzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yzzzzzzz_0[i] = 6.0 * g_xxxxx_0_yzzzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yzzzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_zzzzzzzz_0[i] = 6.0 * g_xxxxx_0_zzzzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_zzzzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 45-90 components of targeted buffer : KSL

    auto g_xxxxxxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 45);

    auto g_xxxxxxy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 46);

    auto g_xxxxxxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 47);

    auto g_xxxxxxy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 48);

    auto g_xxxxxxy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 49);

    auto g_xxxxxxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 50);

    auto g_xxxxxxy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 51);

    auto g_xxxxxxy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 52);

    auto g_xxxxxxy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 53);

    auto g_xxxxxxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 54);

    auto g_xxxxxxy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 55);

    auto g_xxxxxxy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 56);

    auto g_xxxxxxy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 57);

    auto g_xxxxxxy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 58);

    auto g_xxxxxxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 59);

    auto g_xxxxxxy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 60);

    auto g_xxxxxxy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 61);

    auto g_xxxxxxy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 62);

    auto g_xxxxxxy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 63);

    auto g_xxxxxxy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 64);

    auto g_xxxxxxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 65);

    auto g_xxxxxxy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 66);

    auto g_xxxxxxy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 67);

    auto g_xxxxxxy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 68);

    auto g_xxxxxxy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 69);

    auto g_xxxxxxy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 70);

    auto g_xxxxxxy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 71);

    auto g_xxxxxxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 72);

    auto g_xxxxxxy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 73);

    auto g_xxxxxxy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 74);

    auto g_xxxxxxy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 75);

    auto g_xxxxxxy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 76);

    auto g_xxxxxxy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 77);

    auto g_xxxxxxy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 78);

    auto g_xxxxxxy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 79);

    auto g_xxxxxxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 80);

    auto g_xxxxxxy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 81);

    auto g_xxxxxxy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 82);

    auto g_xxxxxxy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 83);

    auto g_xxxxxxy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 84);

    auto g_xxxxxxy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 85);

    auto g_xxxxxxy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 86);

    auto g_xxxxxxy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 87);

    auto g_xxxxxxy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 88);

    auto g_xxxxxxy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 89);

    #pragma omp simd aligned(g_xxxxxx_0_xxxxxxx_1, g_xxxxxx_0_xxxxxxxx_1, g_xxxxxx_0_xxxxxxxy_1, g_xxxxxx_0_xxxxxxxz_1, g_xxxxxx_0_xxxxxxy_1, g_xxxxxx_0_xxxxxxyy_1, g_xxxxxx_0_xxxxxxyz_1, g_xxxxxx_0_xxxxxxz_1, g_xxxxxx_0_xxxxxxzz_1, g_xxxxxx_0_xxxxxyy_1, g_xxxxxx_0_xxxxxyyy_1, g_xxxxxx_0_xxxxxyyz_1, g_xxxxxx_0_xxxxxyz_1, g_xxxxxx_0_xxxxxyzz_1, g_xxxxxx_0_xxxxxzz_1, g_xxxxxx_0_xxxxxzzz_1, g_xxxxxx_0_xxxxyyy_1, g_xxxxxx_0_xxxxyyyy_1, g_xxxxxx_0_xxxxyyyz_1, g_xxxxxx_0_xxxxyyz_1, g_xxxxxx_0_xxxxyyzz_1, g_xxxxxx_0_xxxxyzz_1, g_xxxxxx_0_xxxxyzzz_1, g_xxxxxx_0_xxxxzzz_1, g_xxxxxx_0_xxxxzzzz_1, g_xxxxxx_0_xxxyyyy_1, g_xxxxxx_0_xxxyyyyy_1, g_xxxxxx_0_xxxyyyyz_1, g_xxxxxx_0_xxxyyyz_1, g_xxxxxx_0_xxxyyyzz_1, g_xxxxxx_0_xxxyyzz_1, g_xxxxxx_0_xxxyyzzz_1, g_xxxxxx_0_xxxyzzz_1, g_xxxxxx_0_xxxyzzzz_1, g_xxxxxx_0_xxxzzzz_1, g_xxxxxx_0_xxxzzzzz_1, g_xxxxxx_0_xxyyyyy_1, g_xxxxxx_0_xxyyyyyy_1, g_xxxxxx_0_xxyyyyyz_1, g_xxxxxx_0_xxyyyyz_1, g_xxxxxx_0_xxyyyyzz_1, g_xxxxxx_0_xxyyyzz_1, g_xxxxxx_0_xxyyyzzz_1, g_xxxxxx_0_xxyyzzz_1, g_xxxxxx_0_xxyyzzzz_1, g_xxxxxx_0_xxyzzzz_1, g_xxxxxx_0_xxyzzzzz_1, g_xxxxxx_0_xxzzzzz_1, g_xxxxxx_0_xxzzzzzz_1, g_xxxxxx_0_xyyyyyy_1, g_xxxxxx_0_xyyyyyyy_1, g_xxxxxx_0_xyyyyyyz_1, g_xxxxxx_0_xyyyyyz_1, g_xxxxxx_0_xyyyyyzz_1, g_xxxxxx_0_xyyyyzz_1, g_xxxxxx_0_xyyyyzzz_1, g_xxxxxx_0_xyyyzzz_1, g_xxxxxx_0_xyyyzzzz_1, g_xxxxxx_0_xyyzzzz_1, g_xxxxxx_0_xyyzzzzz_1, g_xxxxxx_0_xyzzzzz_1, g_xxxxxx_0_xyzzzzzz_1, g_xxxxxx_0_xzzzzzz_1, g_xxxxxx_0_xzzzzzzz_1, g_xxxxxx_0_yyyyyyy_1, g_xxxxxx_0_yyyyyyyy_1, g_xxxxxx_0_yyyyyyyz_1, g_xxxxxx_0_yyyyyyz_1, g_xxxxxx_0_yyyyyyzz_1, g_xxxxxx_0_yyyyyzz_1, g_xxxxxx_0_yyyyyzzz_1, g_xxxxxx_0_yyyyzzz_1, g_xxxxxx_0_yyyyzzzz_1, g_xxxxxx_0_yyyzzzz_1, g_xxxxxx_0_yyyzzzzz_1, g_xxxxxx_0_yyzzzzz_1, g_xxxxxx_0_yyzzzzzz_1, g_xxxxxx_0_yzzzzzz_1, g_xxxxxx_0_yzzzzzzz_1, g_xxxxxx_0_zzzzzzz_1, g_xxxxxx_0_zzzzzzzz_1, g_xxxxxxy_0_xxxxxxxx_0, g_xxxxxxy_0_xxxxxxxy_0, g_xxxxxxy_0_xxxxxxxz_0, g_xxxxxxy_0_xxxxxxyy_0, g_xxxxxxy_0_xxxxxxyz_0, g_xxxxxxy_0_xxxxxxzz_0, g_xxxxxxy_0_xxxxxyyy_0, g_xxxxxxy_0_xxxxxyyz_0, g_xxxxxxy_0_xxxxxyzz_0, g_xxxxxxy_0_xxxxxzzz_0, g_xxxxxxy_0_xxxxyyyy_0, g_xxxxxxy_0_xxxxyyyz_0, g_xxxxxxy_0_xxxxyyzz_0, g_xxxxxxy_0_xxxxyzzz_0, g_xxxxxxy_0_xxxxzzzz_0, g_xxxxxxy_0_xxxyyyyy_0, g_xxxxxxy_0_xxxyyyyz_0, g_xxxxxxy_0_xxxyyyzz_0, g_xxxxxxy_0_xxxyyzzz_0, g_xxxxxxy_0_xxxyzzzz_0, g_xxxxxxy_0_xxxzzzzz_0, g_xxxxxxy_0_xxyyyyyy_0, g_xxxxxxy_0_xxyyyyyz_0, g_xxxxxxy_0_xxyyyyzz_0, g_xxxxxxy_0_xxyyyzzz_0, g_xxxxxxy_0_xxyyzzzz_0, g_xxxxxxy_0_xxyzzzzz_0, g_xxxxxxy_0_xxzzzzzz_0, g_xxxxxxy_0_xyyyyyyy_0, g_xxxxxxy_0_xyyyyyyz_0, g_xxxxxxy_0_xyyyyyzz_0, g_xxxxxxy_0_xyyyyzzz_0, g_xxxxxxy_0_xyyyzzzz_0, g_xxxxxxy_0_xyyzzzzz_0, g_xxxxxxy_0_xyzzzzzz_0, g_xxxxxxy_0_xzzzzzzz_0, g_xxxxxxy_0_yyyyyyyy_0, g_xxxxxxy_0_yyyyyyyz_0, g_xxxxxxy_0_yyyyyyzz_0, g_xxxxxxy_0_yyyyyzzz_0, g_xxxxxxy_0_yyyyzzzz_0, g_xxxxxxy_0_yyyzzzzz_0, g_xxxxxxy_0_yyzzzzzz_0, g_xxxxxxy_0_yzzzzzzz_0, g_xxxxxxy_0_zzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxy_0_xxxxxxxx_0[i] = g_xxxxxx_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxxxy_0[i] = g_xxxxxx_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxxy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxxxz_0[i] = g_xxxxxx_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxxyy_0[i] = 2.0 * g_xxxxxx_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxxyz_0[i] = g_xxxxxx_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxxzz_0[i] = g_xxxxxx_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxyyy_0[i] = 3.0 * g_xxxxxx_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxyyz_0[i] = 2.0 * g_xxxxxx_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxyzz_0[i] = g_xxxxxx_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxzzz_0[i] = g_xxxxxx_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxyyyy_0[i] = 4.0 * g_xxxxxx_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxyyyz_0[i] = 3.0 * g_xxxxxx_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxyyzz_0[i] = 2.0 * g_xxxxxx_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxyzzz_0[i] = g_xxxxxx_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxzzzz_0[i] = g_xxxxxx_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyyyyy_0[i] = 5.0 * g_xxxxxx_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyyyyz_0[i] = 4.0 * g_xxxxxx_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyyyzz_0[i] = 3.0 * g_xxxxxx_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyyzzz_0[i] = 2.0 * g_xxxxxx_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyzzzz_0[i] = g_xxxxxx_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxzzzzz_0[i] = g_xxxxxx_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyyyyy_0[i] = 6.0 * g_xxxxxx_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyyyyz_0[i] = 5.0 * g_xxxxxx_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyyyzz_0[i] = 4.0 * g_xxxxxx_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyyzzz_0[i] = 3.0 * g_xxxxxx_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyzzzz_0[i] = 2.0 * g_xxxxxx_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyzzzzz_0[i] = g_xxxxxx_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxzzzzzz_0[i] = g_xxxxxx_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyyyyy_0[i] = 7.0 * g_xxxxxx_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyyyyz_0[i] = 6.0 * g_xxxxxx_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyyyzz_0[i] = 5.0 * g_xxxxxx_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyyzzz_0[i] = 4.0 * g_xxxxxx_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyzzzz_0[i] = 3.0 * g_xxxxxx_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyzzzzz_0[i] = 2.0 * g_xxxxxx_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyzzzzzz_0[i] = g_xxxxxx_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xzzzzzzz_0[i] = g_xxxxxx_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyyyyy_0[i] = 8.0 * g_xxxxxx_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyyyyz_0[i] = 7.0 * g_xxxxxx_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyyyzz_0[i] = 6.0 * g_xxxxxx_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyyzzz_0[i] = 5.0 * g_xxxxxx_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyzzzz_0[i] = 4.0 * g_xxxxxx_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyzzzzz_0[i] = 3.0 * g_xxxxxx_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyzzzzzz_0[i] = 2.0 * g_xxxxxx_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yzzzzzzz_0[i] = g_xxxxxx_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_zzzzzzzz_0[i] = g_xxxxxx_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 90-135 components of targeted buffer : KSL

    auto g_xxxxxxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 90);

    auto g_xxxxxxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 91);

    auto g_xxxxxxz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 92);

    auto g_xxxxxxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 93);

    auto g_xxxxxxz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 94);

    auto g_xxxxxxz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 95);

    auto g_xxxxxxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 96);

    auto g_xxxxxxz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 97);

    auto g_xxxxxxz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 98);

    auto g_xxxxxxz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 99);

    auto g_xxxxxxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 100);

    auto g_xxxxxxz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 101);

    auto g_xxxxxxz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 102);

    auto g_xxxxxxz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 103);

    auto g_xxxxxxz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 104);

    auto g_xxxxxxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 105);

    auto g_xxxxxxz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 106);

    auto g_xxxxxxz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 107);

    auto g_xxxxxxz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 108);

    auto g_xxxxxxz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 109);

    auto g_xxxxxxz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 110);

    auto g_xxxxxxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 111);

    auto g_xxxxxxz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 112);

    auto g_xxxxxxz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 113);

    auto g_xxxxxxz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 114);

    auto g_xxxxxxz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 115);

    auto g_xxxxxxz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 116);

    auto g_xxxxxxz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 117);

    auto g_xxxxxxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 118);

    auto g_xxxxxxz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 119);

    auto g_xxxxxxz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 120);

    auto g_xxxxxxz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 121);

    auto g_xxxxxxz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 122);

    auto g_xxxxxxz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 123);

    auto g_xxxxxxz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 124);

    auto g_xxxxxxz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 125);

    auto g_xxxxxxz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 126);

    auto g_xxxxxxz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 127);

    auto g_xxxxxxz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 128);

    auto g_xxxxxxz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 129);

    auto g_xxxxxxz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 130);

    auto g_xxxxxxz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 131);

    auto g_xxxxxxz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 132);

    auto g_xxxxxxz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 133);

    auto g_xxxxxxz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 134);

    #pragma omp simd aligned(g_xxxxxx_0_xxxxxxx_1, g_xxxxxx_0_xxxxxxxx_1, g_xxxxxx_0_xxxxxxxy_1, g_xxxxxx_0_xxxxxxxz_1, g_xxxxxx_0_xxxxxxy_1, g_xxxxxx_0_xxxxxxyy_1, g_xxxxxx_0_xxxxxxyz_1, g_xxxxxx_0_xxxxxxz_1, g_xxxxxx_0_xxxxxxzz_1, g_xxxxxx_0_xxxxxyy_1, g_xxxxxx_0_xxxxxyyy_1, g_xxxxxx_0_xxxxxyyz_1, g_xxxxxx_0_xxxxxyz_1, g_xxxxxx_0_xxxxxyzz_1, g_xxxxxx_0_xxxxxzz_1, g_xxxxxx_0_xxxxxzzz_1, g_xxxxxx_0_xxxxyyy_1, g_xxxxxx_0_xxxxyyyy_1, g_xxxxxx_0_xxxxyyyz_1, g_xxxxxx_0_xxxxyyz_1, g_xxxxxx_0_xxxxyyzz_1, g_xxxxxx_0_xxxxyzz_1, g_xxxxxx_0_xxxxyzzz_1, g_xxxxxx_0_xxxxzzz_1, g_xxxxxx_0_xxxxzzzz_1, g_xxxxxx_0_xxxyyyy_1, g_xxxxxx_0_xxxyyyyy_1, g_xxxxxx_0_xxxyyyyz_1, g_xxxxxx_0_xxxyyyz_1, g_xxxxxx_0_xxxyyyzz_1, g_xxxxxx_0_xxxyyzz_1, g_xxxxxx_0_xxxyyzzz_1, g_xxxxxx_0_xxxyzzz_1, g_xxxxxx_0_xxxyzzzz_1, g_xxxxxx_0_xxxzzzz_1, g_xxxxxx_0_xxxzzzzz_1, g_xxxxxx_0_xxyyyyy_1, g_xxxxxx_0_xxyyyyyy_1, g_xxxxxx_0_xxyyyyyz_1, g_xxxxxx_0_xxyyyyz_1, g_xxxxxx_0_xxyyyyzz_1, g_xxxxxx_0_xxyyyzz_1, g_xxxxxx_0_xxyyyzzz_1, g_xxxxxx_0_xxyyzzz_1, g_xxxxxx_0_xxyyzzzz_1, g_xxxxxx_0_xxyzzzz_1, g_xxxxxx_0_xxyzzzzz_1, g_xxxxxx_0_xxzzzzz_1, g_xxxxxx_0_xxzzzzzz_1, g_xxxxxx_0_xyyyyyy_1, g_xxxxxx_0_xyyyyyyy_1, g_xxxxxx_0_xyyyyyyz_1, g_xxxxxx_0_xyyyyyz_1, g_xxxxxx_0_xyyyyyzz_1, g_xxxxxx_0_xyyyyzz_1, g_xxxxxx_0_xyyyyzzz_1, g_xxxxxx_0_xyyyzzz_1, g_xxxxxx_0_xyyyzzzz_1, g_xxxxxx_0_xyyzzzz_1, g_xxxxxx_0_xyyzzzzz_1, g_xxxxxx_0_xyzzzzz_1, g_xxxxxx_0_xyzzzzzz_1, g_xxxxxx_0_xzzzzzz_1, g_xxxxxx_0_xzzzzzzz_1, g_xxxxxx_0_yyyyyyy_1, g_xxxxxx_0_yyyyyyyy_1, g_xxxxxx_0_yyyyyyyz_1, g_xxxxxx_0_yyyyyyz_1, g_xxxxxx_0_yyyyyyzz_1, g_xxxxxx_0_yyyyyzz_1, g_xxxxxx_0_yyyyyzzz_1, g_xxxxxx_0_yyyyzzz_1, g_xxxxxx_0_yyyyzzzz_1, g_xxxxxx_0_yyyzzzz_1, g_xxxxxx_0_yyyzzzzz_1, g_xxxxxx_0_yyzzzzz_1, g_xxxxxx_0_yyzzzzzz_1, g_xxxxxx_0_yzzzzzz_1, g_xxxxxx_0_yzzzzzzz_1, g_xxxxxx_0_zzzzzzz_1, g_xxxxxx_0_zzzzzzzz_1, g_xxxxxxz_0_xxxxxxxx_0, g_xxxxxxz_0_xxxxxxxy_0, g_xxxxxxz_0_xxxxxxxz_0, g_xxxxxxz_0_xxxxxxyy_0, g_xxxxxxz_0_xxxxxxyz_0, g_xxxxxxz_0_xxxxxxzz_0, g_xxxxxxz_0_xxxxxyyy_0, g_xxxxxxz_0_xxxxxyyz_0, g_xxxxxxz_0_xxxxxyzz_0, g_xxxxxxz_0_xxxxxzzz_0, g_xxxxxxz_0_xxxxyyyy_0, g_xxxxxxz_0_xxxxyyyz_0, g_xxxxxxz_0_xxxxyyzz_0, g_xxxxxxz_0_xxxxyzzz_0, g_xxxxxxz_0_xxxxzzzz_0, g_xxxxxxz_0_xxxyyyyy_0, g_xxxxxxz_0_xxxyyyyz_0, g_xxxxxxz_0_xxxyyyzz_0, g_xxxxxxz_0_xxxyyzzz_0, g_xxxxxxz_0_xxxyzzzz_0, g_xxxxxxz_0_xxxzzzzz_0, g_xxxxxxz_0_xxyyyyyy_0, g_xxxxxxz_0_xxyyyyyz_0, g_xxxxxxz_0_xxyyyyzz_0, g_xxxxxxz_0_xxyyyzzz_0, g_xxxxxxz_0_xxyyzzzz_0, g_xxxxxxz_0_xxyzzzzz_0, g_xxxxxxz_0_xxzzzzzz_0, g_xxxxxxz_0_xyyyyyyy_0, g_xxxxxxz_0_xyyyyyyz_0, g_xxxxxxz_0_xyyyyyzz_0, g_xxxxxxz_0_xyyyyzzz_0, g_xxxxxxz_0_xyyyzzzz_0, g_xxxxxxz_0_xyyzzzzz_0, g_xxxxxxz_0_xyzzzzzz_0, g_xxxxxxz_0_xzzzzzzz_0, g_xxxxxxz_0_yyyyyyyy_0, g_xxxxxxz_0_yyyyyyyz_0, g_xxxxxxz_0_yyyyyyzz_0, g_xxxxxxz_0_yyyyyzzz_0, g_xxxxxxz_0_yyyyzzzz_0, g_xxxxxxz_0_yyyzzzzz_0, g_xxxxxxz_0_yyzzzzzz_0, g_xxxxxxz_0_yzzzzzzz_0, g_xxxxxxz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxz_0_xxxxxxxx_0[i] = g_xxxxxx_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxxxy_0[i] = g_xxxxxx_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxxxz_0[i] = g_xxxxxx_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxxz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxxyy_0[i] = g_xxxxxx_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxxyz_0[i] = g_xxxxxx_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxxzz_0[i] = 2.0 * g_xxxxxx_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxyyy_0[i] = g_xxxxxx_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxyyz_0[i] = g_xxxxxx_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxyzz_0[i] = 2.0 * g_xxxxxx_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxzzz_0[i] = 3.0 * g_xxxxxx_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxyyyy_0[i] = g_xxxxxx_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxyyyz_0[i] = g_xxxxxx_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxyyzz_0[i] = 2.0 * g_xxxxxx_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxyzzz_0[i] = 3.0 * g_xxxxxx_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxzzzz_0[i] = 4.0 * g_xxxxxx_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyyyyy_0[i] = g_xxxxxx_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyyyyz_0[i] = g_xxxxxx_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyyyzz_0[i] = 2.0 * g_xxxxxx_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyyzzz_0[i] = 3.0 * g_xxxxxx_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyzzzz_0[i] = 4.0 * g_xxxxxx_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxzzzzz_0[i] = 5.0 * g_xxxxxx_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyyyyy_0[i] = g_xxxxxx_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyyyyz_0[i] = g_xxxxxx_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyyyzz_0[i] = 2.0 * g_xxxxxx_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyyzzz_0[i] = 3.0 * g_xxxxxx_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyzzzz_0[i] = 4.0 * g_xxxxxx_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyzzzzz_0[i] = 5.0 * g_xxxxxx_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxzzzzzz_0[i] = 6.0 * g_xxxxxx_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxzzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyyyyy_0[i] = g_xxxxxx_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyyyyz_0[i] = g_xxxxxx_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyyyzz_0[i] = 2.0 * g_xxxxxx_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyyzzz_0[i] = 3.0 * g_xxxxxx_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyzzzz_0[i] = 4.0 * g_xxxxxx_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyzzzzz_0[i] = 5.0 * g_xxxxxx_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyzzzzzz_0[i] = 6.0 * g_xxxxxx_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xzzzzzzz_0[i] = 7.0 * g_xxxxxx_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xzzzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyyyyy_0[i] = g_xxxxxx_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyyyyz_0[i] = g_xxxxxx_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyyyzz_0[i] = 2.0 * g_xxxxxx_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyyzzz_0[i] = 3.0 * g_xxxxxx_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyzzzz_0[i] = 4.0 * g_xxxxxx_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyzzzzz_0[i] = 5.0 * g_xxxxxx_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyzzzzzz_0[i] = 6.0 * g_xxxxxx_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyzzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yzzzzzzz_0[i] = 7.0 * g_xxxxxx_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yzzzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_zzzzzzzz_0[i] = 8.0 * g_xxxxxx_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 135-180 components of targeted buffer : KSL

    auto g_xxxxxyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 135);

    auto g_xxxxxyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 136);

    auto g_xxxxxyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 137);

    auto g_xxxxxyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 138);

    auto g_xxxxxyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 139);

    auto g_xxxxxyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 140);

    auto g_xxxxxyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 141);

    auto g_xxxxxyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 142);

    auto g_xxxxxyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 143);

    auto g_xxxxxyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 144);

    auto g_xxxxxyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 145);

    auto g_xxxxxyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 146);

    auto g_xxxxxyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 147);

    auto g_xxxxxyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 148);

    auto g_xxxxxyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 149);

    auto g_xxxxxyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 150);

    auto g_xxxxxyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 151);

    auto g_xxxxxyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 152);

    auto g_xxxxxyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 153);

    auto g_xxxxxyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 154);

    auto g_xxxxxyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 155);

    auto g_xxxxxyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 156);

    auto g_xxxxxyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 157);

    auto g_xxxxxyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 158);

    auto g_xxxxxyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 159);

    auto g_xxxxxyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 160);

    auto g_xxxxxyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 161);

    auto g_xxxxxyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 162);

    auto g_xxxxxyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 163);

    auto g_xxxxxyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 164);

    auto g_xxxxxyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 165);

    auto g_xxxxxyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 166);

    auto g_xxxxxyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 167);

    auto g_xxxxxyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 168);

    auto g_xxxxxyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 169);

    auto g_xxxxxyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 170);

    auto g_xxxxxyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 171);

    auto g_xxxxxyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 172);

    auto g_xxxxxyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 173);

    auto g_xxxxxyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 174);

    auto g_xxxxxyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 175);

    auto g_xxxxxyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 176);

    auto g_xxxxxyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 177);

    auto g_xxxxxyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 178);

    auto g_xxxxxyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 179);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxxxx_0, g_xxxxx_0_xxxxxxxx_1, g_xxxxx_0_xxxxxxxz_0, g_xxxxx_0_xxxxxxxz_1, g_xxxxx_0_xxxxxxzz_0, g_xxxxx_0_xxxxxxzz_1, g_xxxxx_0_xxxxxzzz_0, g_xxxxx_0_xxxxxzzz_1, g_xxxxx_0_xxxxzzzz_0, g_xxxxx_0_xxxxzzzz_1, g_xxxxx_0_xxxzzzzz_0, g_xxxxx_0_xxxzzzzz_1, g_xxxxx_0_xxzzzzzz_0, g_xxxxx_0_xxzzzzzz_1, g_xxxxx_0_xzzzzzzz_0, g_xxxxx_0_xzzzzzzz_1, g_xxxxxy_0_xxxxxxxx_1, g_xxxxxy_0_xxxxxxxz_1, g_xxxxxy_0_xxxxxxzz_1, g_xxxxxy_0_xxxxxzzz_1, g_xxxxxy_0_xxxxzzzz_1, g_xxxxxy_0_xxxzzzzz_1, g_xxxxxy_0_xxzzzzzz_1, g_xxxxxy_0_xzzzzzzz_1, g_xxxxxyy_0_xxxxxxxx_0, g_xxxxxyy_0_xxxxxxxy_0, g_xxxxxyy_0_xxxxxxxz_0, g_xxxxxyy_0_xxxxxxyy_0, g_xxxxxyy_0_xxxxxxyz_0, g_xxxxxyy_0_xxxxxxzz_0, g_xxxxxyy_0_xxxxxyyy_0, g_xxxxxyy_0_xxxxxyyz_0, g_xxxxxyy_0_xxxxxyzz_0, g_xxxxxyy_0_xxxxxzzz_0, g_xxxxxyy_0_xxxxyyyy_0, g_xxxxxyy_0_xxxxyyyz_0, g_xxxxxyy_0_xxxxyyzz_0, g_xxxxxyy_0_xxxxyzzz_0, g_xxxxxyy_0_xxxxzzzz_0, g_xxxxxyy_0_xxxyyyyy_0, g_xxxxxyy_0_xxxyyyyz_0, g_xxxxxyy_0_xxxyyyzz_0, g_xxxxxyy_0_xxxyyzzz_0, g_xxxxxyy_0_xxxyzzzz_0, g_xxxxxyy_0_xxxzzzzz_0, g_xxxxxyy_0_xxyyyyyy_0, g_xxxxxyy_0_xxyyyyyz_0, g_xxxxxyy_0_xxyyyyzz_0, g_xxxxxyy_0_xxyyyzzz_0, g_xxxxxyy_0_xxyyzzzz_0, g_xxxxxyy_0_xxyzzzzz_0, g_xxxxxyy_0_xxzzzzzz_0, g_xxxxxyy_0_xyyyyyyy_0, g_xxxxxyy_0_xyyyyyyz_0, g_xxxxxyy_0_xyyyyyzz_0, g_xxxxxyy_0_xyyyyzzz_0, g_xxxxxyy_0_xyyyzzzz_0, g_xxxxxyy_0_xyyzzzzz_0, g_xxxxxyy_0_xyzzzzzz_0, g_xxxxxyy_0_xzzzzzzz_0, g_xxxxxyy_0_yyyyyyyy_0, g_xxxxxyy_0_yyyyyyyz_0, g_xxxxxyy_0_yyyyyyzz_0, g_xxxxxyy_0_yyyyyzzz_0, g_xxxxxyy_0_yyyyzzzz_0, g_xxxxxyy_0_yyyzzzzz_0, g_xxxxxyy_0_yyzzzzzz_0, g_xxxxxyy_0_yzzzzzzz_0, g_xxxxxyy_0_zzzzzzzz_0, g_xxxxyy_0_xxxxxxxy_1, g_xxxxyy_0_xxxxxxy_1, g_xxxxyy_0_xxxxxxyy_1, g_xxxxyy_0_xxxxxxyz_1, g_xxxxyy_0_xxxxxyy_1, g_xxxxyy_0_xxxxxyyy_1, g_xxxxyy_0_xxxxxyyz_1, g_xxxxyy_0_xxxxxyz_1, g_xxxxyy_0_xxxxxyzz_1, g_xxxxyy_0_xxxxyyy_1, g_xxxxyy_0_xxxxyyyy_1, g_xxxxyy_0_xxxxyyyz_1, g_xxxxyy_0_xxxxyyz_1, g_xxxxyy_0_xxxxyyzz_1, g_xxxxyy_0_xxxxyzz_1, g_xxxxyy_0_xxxxyzzz_1, g_xxxxyy_0_xxxyyyy_1, g_xxxxyy_0_xxxyyyyy_1, g_xxxxyy_0_xxxyyyyz_1, g_xxxxyy_0_xxxyyyz_1, g_xxxxyy_0_xxxyyyzz_1, g_xxxxyy_0_xxxyyzz_1, g_xxxxyy_0_xxxyyzzz_1, g_xxxxyy_0_xxxyzzz_1, g_xxxxyy_0_xxxyzzzz_1, g_xxxxyy_0_xxyyyyy_1, g_xxxxyy_0_xxyyyyyy_1, g_xxxxyy_0_xxyyyyyz_1, g_xxxxyy_0_xxyyyyz_1, g_xxxxyy_0_xxyyyyzz_1, g_xxxxyy_0_xxyyyzz_1, g_xxxxyy_0_xxyyyzzz_1, g_xxxxyy_0_xxyyzzz_1, g_xxxxyy_0_xxyyzzzz_1, g_xxxxyy_0_xxyzzzz_1, g_xxxxyy_0_xxyzzzzz_1, g_xxxxyy_0_xyyyyyy_1, g_xxxxyy_0_xyyyyyyy_1, g_xxxxyy_0_xyyyyyyz_1, g_xxxxyy_0_xyyyyyz_1, g_xxxxyy_0_xyyyyyzz_1, g_xxxxyy_0_xyyyyzz_1, g_xxxxyy_0_xyyyyzzz_1, g_xxxxyy_0_xyyyzzz_1, g_xxxxyy_0_xyyyzzzz_1, g_xxxxyy_0_xyyzzzz_1, g_xxxxyy_0_xyyzzzzz_1, g_xxxxyy_0_xyzzzzz_1, g_xxxxyy_0_xyzzzzzz_1, g_xxxxyy_0_yyyyyyy_1, g_xxxxyy_0_yyyyyyyy_1, g_xxxxyy_0_yyyyyyyz_1, g_xxxxyy_0_yyyyyyz_1, g_xxxxyy_0_yyyyyyzz_1, g_xxxxyy_0_yyyyyzz_1, g_xxxxyy_0_yyyyyzzz_1, g_xxxxyy_0_yyyyzzz_1, g_xxxxyy_0_yyyyzzzz_1, g_xxxxyy_0_yyyzzzz_1, g_xxxxyy_0_yyyzzzzz_1, g_xxxxyy_0_yyzzzzz_1, g_xxxxyy_0_yyzzzzzz_1, g_xxxxyy_0_yzzzzzz_1, g_xxxxyy_0_yzzzzzzz_1, g_xxxxyy_0_zzzzzzzz_1, g_xxxyy_0_xxxxxxxy_0, g_xxxyy_0_xxxxxxxy_1, g_xxxyy_0_xxxxxxyy_0, g_xxxyy_0_xxxxxxyy_1, g_xxxyy_0_xxxxxxyz_0, g_xxxyy_0_xxxxxxyz_1, g_xxxyy_0_xxxxxyyy_0, g_xxxyy_0_xxxxxyyy_1, g_xxxyy_0_xxxxxyyz_0, g_xxxyy_0_xxxxxyyz_1, g_xxxyy_0_xxxxxyzz_0, g_xxxyy_0_xxxxxyzz_1, g_xxxyy_0_xxxxyyyy_0, g_xxxyy_0_xxxxyyyy_1, g_xxxyy_0_xxxxyyyz_0, g_xxxyy_0_xxxxyyyz_1, g_xxxyy_0_xxxxyyzz_0, g_xxxyy_0_xxxxyyzz_1, g_xxxyy_0_xxxxyzzz_0, g_xxxyy_0_xxxxyzzz_1, g_xxxyy_0_xxxyyyyy_0, g_xxxyy_0_xxxyyyyy_1, g_xxxyy_0_xxxyyyyz_0, g_xxxyy_0_xxxyyyyz_1, g_xxxyy_0_xxxyyyzz_0, g_xxxyy_0_xxxyyyzz_1, g_xxxyy_0_xxxyyzzz_0, g_xxxyy_0_xxxyyzzz_1, g_xxxyy_0_xxxyzzzz_0, g_xxxyy_0_xxxyzzzz_1, g_xxxyy_0_xxyyyyyy_0, g_xxxyy_0_xxyyyyyy_1, g_xxxyy_0_xxyyyyyz_0, g_xxxyy_0_xxyyyyyz_1, g_xxxyy_0_xxyyyyzz_0, g_xxxyy_0_xxyyyyzz_1, g_xxxyy_0_xxyyyzzz_0, g_xxxyy_0_xxyyyzzz_1, g_xxxyy_0_xxyyzzzz_0, g_xxxyy_0_xxyyzzzz_1, g_xxxyy_0_xxyzzzzz_0, g_xxxyy_0_xxyzzzzz_1, g_xxxyy_0_xyyyyyyy_0, g_xxxyy_0_xyyyyyyy_1, g_xxxyy_0_xyyyyyyz_0, g_xxxyy_0_xyyyyyyz_1, g_xxxyy_0_xyyyyyzz_0, g_xxxyy_0_xyyyyyzz_1, g_xxxyy_0_xyyyyzzz_0, g_xxxyy_0_xyyyyzzz_1, g_xxxyy_0_xyyyzzzz_0, g_xxxyy_0_xyyyzzzz_1, g_xxxyy_0_xyyzzzzz_0, g_xxxyy_0_xyyzzzzz_1, g_xxxyy_0_xyzzzzzz_0, g_xxxyy_0_xyzzzzzz_1, g_xxxyy_0_yyyyyyyy_0, g_xxxyy_0_yyyyyyyy_1, g_xxxyy_0_yyyyyyyz_0, g_xxxyy_0_yyyyyyyz_1, g_xxxyy_0_yyyyyyzz_0, g_xxxyy_0_yyyyyyzz_1, g_xxxyy_0_yyyyyzzz_0, g_xxxyy_0_yyyyyzzz_1, g_xxxyy_0_yyyyzzzz_0, g_xxxyy_0_yyyyzzzz_1, g_xxxyy_0_yyyzzzzz_0, g_xxxyy_0_yyyzzzzz_1, g_xxxyy_0_yyzzzzzz_0, g_xxxyy_0_yyzzzzzz_1, g_xxxyy_0_yzzzzzzz_0, g_xxxyy_0_yzzzzzzz_1, g_xxxyy_0_zzzzzzzz_0, g_xxxyy_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxyy_0_xxxxxxxx_0[i] = g_xxxxx_0_xxxxxxxx_0[i] * fbe_0 - g_xxxxx_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxxxy_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxxxxxy_0[i] = 4.0 * g_xxxyy_0_xxxxxxxy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xxxxyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxxxxz_0[i] = g_xxxxx_0_xxxxxxxz_0[i] * fbe_0 - g_xxxxx_0_xxxxxxxz_1[i] * fz_be_0 + g_xxxxxy_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxxxxyy_0[i] = 4.0 * g_xxxyy_0_xxxxxxyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xxxxyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxxxyz_0[i] = 4.0 * g_xxxyy_0_xxxxxxyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxxxyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxxxzz_0[i] = g_xxxxx_0_xxxxxxzz_0[i] * fbe_0 - g_xxxxx_0_xxxxxxzz_1[i] * fz_be_0 + g_xxxxxy_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxxxyyy_0[i] = 4.0 * g_xxxyy_0_xxxxxyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xxxxyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxxyyz_0[i] = 4.0 * g_xxxyy_0_xxxxxyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxxxyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxxyzz_0[i] = 4.0 * g_xxxyy_0_xxxxxyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxxxyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxxzzz_0[i] = g_xxxxx_0_xxxxxzzz_0[i] * fbe_0 - g_xxxxx_0_xxxxxzzz_1[i] * fz_be_0 + g_xxxxxy_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxxyyyy_0[i] = 4.0 * g_xxxyy_0_xxxxyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xxxxyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxyyyz_0[i] = 4.0 * g_xxxyy_0_xxxxyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxxxyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxyyzz_0[i] = 4.0 * g_xxxyy_0_xxxxyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxxxyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxyzzz_0[i] = 4.0 * g_xxxyy_0_xxxxyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxxxyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxzzzz_0[i] = g_xxxxx_0_xxxxzzzz_0[i] * fbe_0 - g_xxxxx_0_xxxxzzzz_1[i] * fz_be_0 + g_xxxxxy_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxyyyyy_0[i] = 4.0 * g_xxxyy_0_xxxyyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxyyyyz_0[i] = 4.0 * g_xxxyy_0_xxxyyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxyyyzz_0[i] = 4.0 * g_xxxyy_0_xxxyyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxyyzzz_0[i] = 4.0 * g_xxxyy_0_xxxyyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxyzzzz_0[i] = 4.0 * g_xxxyy_0_xxxyzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxzzzzz_0[i] = g_xxxxx_0_xxxzzzzz_0[i] * fbe_0 - g_xxxxx_0_xxxzzzzz_1[i] * fz_be_0 + g_xxxxxy_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxyyyyyy_0[i] = 4.0 * g_xxxyy_0_xxyyyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyyyyyz_0[i] = 4.0 * g_xxxyy_0_xxyyyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyyyyzz_0[i] = 4.0 * g_xxxyy_0_xxyyyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyyyzzz_0[i] = 4.0 * g_xxxyy_0_xxyyyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyyzzzz_0[i] = 4.0 * g_xxxyy_0_xxyyzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyzzzzz_0[i] = 4.0 * g_xxxyy_0_xxyzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxzzzzzz_0[i] = g_xxxxx_0_xxzzzzzz_0[i] * fbe_0 - g_xxxxx_0_xxzzzzzz_1[i] * fz_be_0 + g_xxxxxy_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xyyyyyyy_0[i] = 4.0 * g_xxxyy_0_xyyyyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyyyyyz_0[i] = 4.0 * g_xxxyy_0_xyyyyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyyyyz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyyyyzz_0[i] = 4.0 * g_xxxyy_0_xyyyyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyyyzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyyyzzz_0[i] = 4.0 * g_xxxyy_0_xyyyyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyyzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyyzzzz_0[i] = 4.0 * g_xxxyy_0_xyyyzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyzzzzz_0[i] = 4.0 * g_xxxyy_0_xyyzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyzzzzzz_0[i] = 4.0 * g_xxxyy_0_xyzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyzzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xzzzzzzz_0[i] = g_xxxxx_0_xzzzzzzz_0[i] * fbe_0 - g_xxxxx_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxxxy_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_yyyyyyyy_0[i] = 4.0 * g_xxxyy_0_yyyyyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyyyyy_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyyyyyz_0[i] = 4.0 * g_xxxyy_0_yyyyyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyyyyz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyyyyzz_0[i] = 4.0 * g_xxxyy_0_yyyyyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyyyzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyyyzzz_0[i] = 4.0 * g_xxxyy_0_yyyyyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyyzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyyzzzz_0[i] = 4.0 * g_xxxyy_0_yyyyzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyzzzzz_0[i] = 4.0 * g_xxxyy_0_yyyzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyzzzzzz_0[i] = 4.0 * g_xxxyy_0_yyzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyzzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yzzzzzzz_0[i] = 4.0 * g_xxxyy_0_yzzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yzzzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_zzzzzzzz_0[i] = 4.0 * g_xxxyy_0_zzzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_zzzzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 180-225 components of targeted buffer : KSL

    auto g_xxxxxyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 180);

    auto g_xxxxxyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 181);

    auto g_xxxxxyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 182);

    auto g_xxxxxyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 183);

    auto g_xxxxxyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 184);

    auto g_xxxxxyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 185);

    auto g_xxxxxyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 186);

    auto g_xxxxxyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 187);

    auto g_xxxxxyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 188);

    auto g_xxxxxyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 189);

    auto g_xxxxxyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 190);

    auto g_xxxxxyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 191);

    auto g_xxxxxyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 192);

    auto g_xxxxxyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 193);

    auto g_xxxxxyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 194);

    auto g_xxxxxyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 195);

    auto g_xxxxxyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 196);

    auto g_xxxxxyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 197);

    auto g_xxxxxyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 198);

    auto g_xxxxxyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 199);

    auto g_xxxxxyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 200);

    auto g_xxxxxyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 201);

    auto g_xxxxxyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 202);

    auto g_xxxxxyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 203);

    auto g_xxxxxyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 204);

    auto g_xxxxxyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 205);

    auto g_xxxxxyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 206);

    auto g_xxxxxyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 207);

    auto g_xxxxxyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 208);

    auto g_xxxxxyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 209);

    auto g_xxxxxyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 210);

    auto g_xxxxxyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 211);

    auto g_xxxxxyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 212);

    auto g_xxxxxyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 213);

    auto g_xxxxxyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 214);

    auto g_xxxxxyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 215);

    auto g_xxxxxyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 216);

    auto g_xxxxxyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 217);

    auto g_xxxxxyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 218);

    auto g_xxxxxyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 219);

    auto g_xxxxxyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 220);

    auto g_xxxxxyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 221);

    auto g_xxxxxyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 222);

    auto g_xxxxxyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 223);

    auto g_xxxxxyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 224);

    #pragma omp simd aligned(g_xxxxxy_0_xxxxxxxy_1, g_xxxxxy_0_xxxxxxyy_1, g_xxxxxy_0_xxxxxyyy_1, g_xxxxxy_0_xxxxyyyy_1, g_xxxxxy_0_xxxyyyyy_1, g_xxxxxy_0_xxyyyyyy_1, g_xxxxxy_0_xyyyyyyy_1, g_xxxxxy_0_yyyyyyyy_1, g_xxxxxyz_0_xxxxxxxx_0, g_xxxxxyz_0_xxxxxxxy_0, g_xxxxxyz_0_xxxxxxxz_0, g_xxxxxyz_0_xxxxxxyy_0, g_xxxxxyz_0_xxxxxxyz_0, g_xxxxxyz_0_xxxxxxzz_0, g_xxxxxyz_0_xxxxxyyy_0, g_xxxxxyz_0_xxxxxyyz_0, g_xxxxxyz_0_xxxxxyzz_0, g_xxxxxyz_0_xxxxxzzz_0, g_xxxxxyz_0_xxxxyyyy_0, g_xxxxxyz_0_xxxxyyyz_0, g_xxxxxyz_0_xxxxyyzz_0, g_xxxxxyz_0_xxxxyzzz_0, g_xxxxxyz_0_xxxxzzzz_0, g_xxxxxyz_0_xxxyyyyy_0, g_xxxxxyz_0_xxxyyyyz_0, g_xxxxxyz_0_xxxyyyzz_0, g_xxxxxyz_0_xxxyyzzz_0, g_xxxxxyz_0_xxxyzzzz_0, g_xxxxxyz_0_xxxzzzzz_0, g_xxxxxyz_0_xxyyyyyy_0, g_xxxxxyz_0_xxyyyyyz_0, g_xxxxxyz_0_xxyyyyzz_0, g_xxxxxyz_0_xxyyyzzz_0, g_xxxxxyz_0_xxyyzzzz_0, g_xxxxxyz_0_xxyzzzzz_0, g_xxxxxyz_0_xxzzzzzz_0, g_xxxxxyz_0_xyyyyyyy_0, g_xxxxxyz_0_xyyyyyyz_0, g_xxxxxyz_0_xyyyyyzz_0, g_xxxxxyz_0_xyyyyzzz_0, g_xxxxxyz_0_xyyyzzzz_0, g_xxxxxyz_0_xyyzzzzz_0, g_xxxxxyz_0_xyzzzzzz_0, g_xxxxxyz_0_xzzzzzzz_0, g_xxxxxyz_0_yyyyyyyy_0, g_xxxxxyz_0_yyyyyyyz_0, g_xxxxxyz_0_yyyyyyzz_0, g_xxxxxyz_0_yyyyyzzz_0, g_xxxxxyz_0_yyyyzzzz_0, g_xxxxxyz_0_yyyzzzzz_0, g_xxxxxyz_0_yyzzzzzz_0, g_xxxxxyz_0_yzzzzzzz_0, g_xxxxxyz_0_zzzzzzzz_0, g_xxxxxz_0_xxxxxxxx_1, g_xxxxxz_0_xxxxxxxz_1, g_xxxxxz_0_xxxxxxyz_1, g_xxxxxz_0_xxxxxxz_1, g_xxxxxz_0_xxxxxxzz_1, g_xxxxxz_0_xxxxxyyz_1, g_xxxxxz_0_xxxxxyz_1, g_xxxxxz_0_xxxxxyzz_1, g_xxxxxz_0_xxxxxzz_1, g_xxxxxz_0_xxxxxzzz_1, g_xxxxxz_0_xxxxyyyz_1, g_xxxxxz_0_xxxxyyz_1, g_xxxxxz_0_xxxxyyzz_1, g_xxxxxz_0_xxxxyzz_1, g_xxxxxz_0_xxxxyzzz_1, g_xxxxxz_0_xxxxzzz_1, g_xxxxxz_0_xxxxzzzz_1, g_xxxxxz_0_xxxyyyyz_1, g_xxxxxz_0_xxxyyyz_1, g_xxxxxz_0_xxxyyyzz_1, g_xxxxxz_0_xxxyyzz_1, g_xxxxxz_0_xxxyyzzz_1, g_xxxxxz_0_xxxyzzz_1, g_xxxxxz_0_xxxyzzzz_1, g_xxxxxz_0_xxxzzzz_1, g_xxxxxz_0_xxxzzzzz_1, g_xxxxxz_0_xxyyyyyz_1, g_xxxxxz_0_xxyyyyz_1, g_xxxxxz_0_xxyyyyzz_1, g_xxxxxz_0_xxyyyzz_1, g_xxxxxz_0_xxyyyzzz_1, g_xxxxxz_0_xxyyzzz_1, g_xxxxxz_0_xxyyzzzz_1, g_xxxxxz_0_xxyzzzz_1, g_xxxxxz_0_xxyzzzzz_1, g_xxxxxz_0_xxzzzzz_1, g_xxxxxz_0_xxzzzzzz_1, g_xxxxxz_0_xyyyyyyz_1, g_xxxxxz_0_xyyyyyz_1, g_xxxxxz_0_xyyyyyzz_1, g_xxxxxz_0_xyyyyzz_1, g_xxxxxz_0_xyyyyzzz_1, g_xxxxxz_0_xyyyzzz_1, g_xxxxxz_0_xyyyzzzz_1, g_xxxxxz_0_xyyzzzz_1, g_xxxxxz_0_xyyzzzzz_1, g_xxxxxz_0_xyzzzzz_1, g_xxxxxz_0_xyzzzzzz_1, g_xxxxxz_0_xzzzzzz_1, g_xxxxxz_0_xzzzzzzz_1, g_xxxxxz_0_yyyyyyyz_1, g_xxxxxz_0_yyyyyyz_1, g_xxxxxz_0_yyyyyyzz_1, g_xxxxxz_0_yyyyyzz_1, g_xxxxxz_0_yyyyyzzz_1, g_xxxxxz_0_yyyyzzz_1, g_xxxxxz_0_yyyyzzzz_1, g_xxxxxz_0_yyyzzzz_1, g_xxxxxz_0_yyyzzzzz_1, g_xxxxxz_0_yyzzzzz_1, g_xxxxxz_0_yyzzzzzz_1, g_xxxxxz_0_yzzzzzz_1, g_xxxxxz_0_yzzzzzzz_1, g_xxxxxz_0_zzzzzzz_1, g_xxxxxz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxyz_0_xxxxxxxx_0[i] = g_xxxxxz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxxxxy_0[i] = g_xxxxxy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxxxxxz_0[i] = g_xxxxxz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxxxyy_0[i] = g_xxxxxy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxxxxyz_0[i] = g_xxxxxz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxxxzz_0[i] = g_xxxxxz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxxyyy_0[i] = g_xxxxxy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxxxyyz_0[i] = 2.0 * g_xxxxxz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxxyzz_0[i] = g_xxxxxz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxxzzz_0[i] = g_xxxxxz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxyyyy_0[i] = g_xxxxxy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxxyyyz_0[i] = 3.0 * g_xxxxxz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxyyzz_0[i] = 2.0 * g_xxxxxz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxyzzz_0[i] = g_xxxxxz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxzzzz_0[i] = g_xxxxxz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxyyyyy_0[i] = g_xxxxxy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxyyyyz_0[i] = 4.0 * g_xxxxxz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxyyyzz_0[i] = 3.0 * g_xxxxxz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxyyzzz_0[i] = 2.0 * g_xxxxxz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxyzzzz_0[i] = g_xxxxxz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxzzzzz_0[i] = g_xxxxxz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyyyyyy_0[i] = g_xxxxxy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxyyyyyz_0[i] = 5.0 * g_xxxxxz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyyyyzz_0[i] = 4.0 * g_xxxxxz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyyyzzz_0[i] = 3.0 * g_xxxxxz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyyzzzz_0[i] = 2.0 * g_xxxxxz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyzzzzz_0[i] = g_xxxxxz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxzzzzzz_0[i] = g_xxxxxz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyyyyyy_0[i] = g_xxxxxy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xyyyyyyz_0[i] = 6.0 * g_xxxxxz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyyyyzz_0[i] = 5.0 * g_xxxxxz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyyyzzz_0[i] = 4.0 * g_xxxxxz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyyzzzz_0[i] = 3.0 * g_xxxxxz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyzzzzz_0[i] = 2.0 * g_xxxxxz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyzzzzzz_0[i] = g_xxxxxz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xzzzzzzz_0[i] = g_xxxxxz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyyyyyy_0[i] = g_xxxxxy_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_yyyyyyyz_0[i] = 7.0 * g_xxxxxz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyyyyzz_0[i] = 6.0 * g_xxxxxz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyyyzzz_0[i] = 5.0 * g_xxxxxz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyyzzzz_0[i] = 4.0 * g_xxxxxz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyzzzzz_0[i] = 3.0 * g_xxxxxz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyzzzzzz_0[i] = 2.0 * g_xxxxxz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yzzzzzzz_0[i] = g_xxxxxz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_zzzzzzzz_0[i] = g_xxxxxz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 225-270 components of targeted buffer : KSL

    auto g_xxxxxzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 225);

    auto g_xxxxxzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 226);

    auto g_xxxxxzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 227);

    auto g_xxxxxzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 228);

    auto g_xxxxxzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 229);

    auto g_xxxxxzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 230);

    auto g_xxxxxzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 231);

    auto g_xxxxxzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 232);

    auto g_xxxxxzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 233);

    auto g_xxxxxzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 234);

    auto g_xxxxxzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 235);

    auto g_xxxxxzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 236);

    auto g_xxxxxzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 237);

    auto g_xxxxxzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 238);

    auto g_xxxxxzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 239);

    auto g_xxxxxzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 240);

    auto g_xxxxxzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 241);

    auto g_xxxxxzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 242);

    auto g_xxxxxzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 243);

    auto g_xxxxxzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 244);

    auto g_xxxxxzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 245);

    auto g_xxxxxzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 246);

    auto g_xxxxxzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 247);

    auto g_xxxxxzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 248);

    auto g_xxxxxzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 249);

    auto g_xxxxxzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 250);

    auto g_xxxxxzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 251);

    auto g_xxxxxzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 252);

    auto g_xxxxxzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 253);

    auto g_xxxxxzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 254);

    auto g_xxxxxzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 255);

    auto g_xxxxxzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 256);

    auto g_xxxxxzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 257);

    auto g_xxxxxzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 258);

    auto g_xxxxxzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 259);

    auto g_xxxxxzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 260);

    auto g_xxxxxzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 261);

    auto g_xxxxxzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 262);

    auto g_xxxxxzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 263);

    auto g_xxxxxzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 264);

    auto g_xxxxxzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 265);

    auto g_xxxxxzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 266);

    auto g_xxxxxzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 267);

    auto g_xxxxxzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 268);

    auto g_xxxxxzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 269);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxxxx_0, g_xxxxx_0_xxxxxxxx_1, g_xxxxx_0_xxxxxxxy_0, g_xxxxx_0_xxxxxxxy_1, g_xxxxx_0_xxxxxxyy_0, g_xxxxx_0_xxxxxxyy_1, g_xxxxx_0_xxxxxyyy_0, g_xxxxx_0_xxxxxyyy_1, g_xxxxx_0_xxxxyyyy_0, g_xxxxx_0_xxxxyyyy_1, g_xxxxx_0_xxxyyyyy_0, g_xxxxx_0_xxxyyyyy_1, g_xxxxx_0_xxyyyyyy_0, g_xxxxx_0_xxyyyyyy_1, g_xxxxx_0_xyyyyyyy_0, g_xxxxx_0_xyyyyyyy_1, g_xxxxxz_0_xxxxxxxx_1, g_xxxxxz_0_xxxxxxxy_1, g_xxxxxz_0_xxxxxxyy_1, g_xxxxxz_0_xxxxxyyy_1, g_xxxxxz_0_xxxxyyyy_1, g_xxxxxz_0_xxxyyyyy_1, g_xxxxxz_0_xxyyyyyy_1, g_xxxxxz_0_xyyyyyyy_1, g_xxxxxzz_0_xxxxxxxx_0, g_xxxxxzz_0_xxxxxxxy_0, g_xxxxxzz_0_xxxxxxxz_0, g_xxxxxzz_0_xxxxxxyy_0, g_xxxxxzz_0_xxxxxxyz_0, g_xxxxxzz_0_xxxxxxzz_0, g_xxxxxzz_0_xxxxxyyy_0, g_xxxxxzz_0_xxxxxyyz_0, g_xxxxxzz_0_xxxxxyzz_0, g_xxxxxzz_0_xxxxxzzz_0, g_xxxxxzz_0_xxxxyyyy_0, g_xxxxxzz_0_xxxxyyyz_0, g_xxxxxzz_0_xxxxyyzz_0, g_xxxxxzz_0_xxxxyzzz_0, g_xxxxxzz_0_xxxxzzzz_0, g_xxxxxzz_0_xxxyyyyy_0, g_xxxxxzz_0_xxxyyyyz_0, g_xxxxxzz_0_xxxyyyzz_0, g_xxxxxzz_0_xxxyyzzz_0, g_xxxxxzz_0_xxxyzzzz_0, g_xxxxxzz_0_xxxzzzzz_0, g_xxxxxzz_0_xxyyyyyy_0, g_xxxxxzz_0_xxyyyyyz_0, g_xxxxxzz_0_xxyyyyzz_0, g_xxxxxzz_0_xxyyyzzz_0, g_xxxxxzz_0_xxyyzzzz_0, g_xxxxxzz_0_xxyzzzzz_0, g_xxxxxzz_0_xxzzzzzz_0, g_xxxxxzz_0_xyyyyyyy_0, g_xxxxxzz_0_xyyyyyyz_0, g_xxxxxzz_0_xyyyyyzz_0, g_xxxxxzz_0_xyyyyzzz_0, g_xxxxxzz_0_xyyyzzzz_0, g_xxxxxzz_0_xyyzzzzz_0, g_xxxxxzz_0_xyzzzzzz_0, g_xxxxxzz_0_xzzzzzzz_0, g_xxxxxzz_0_yyyyyyyy_0, g_xxxxxzz_0_yyyyyyyz_0, g_xxxxxzz_0_yyyyyyzz_0, g_xxxxxzz_0_yyyyyzzz_0, g_xxxxxzz_0_yyyyzzzz_0, g_xxxxxzz_0_yyyzzzzz_0, g_xxxxxzz_0_yyzzzzzz_0, g_xxxxxzz_0_yzzzzzzz_0, g_xxxxxzz_0_zzzzzzzz_0, g_xxxxzz_0_xxxxxxxz_1, g_xxxxzz_0_xxxxxxyz_1, g_xxxxzz_0_xxxxxxz_1, g_xxxxzz_0_xxxxxxzz_1, g_xxxxzz_0_xxxxxyyz_1, g_xxxxzz_0_xxxxxyz_1, g_xxxxzz_0_xxxxxyzz_1, g_xxxxzz_0_xxxxxzz_1, g_xxxxzz_0_xxxxxzzz_1, g_xxxxzz_0_xxxxyyyz_1, g_xxxxzz_0_xxxxyyz_1, g_xxxxzz_0_xxxxyyzz_1, g_xxxxzz_0_xxxxyzz_1, g_xxxxzz_0_xxxxyzzz_1, g_xxxxzz_0_xxxxzzz_1, g_xxxxzz_0_xxxxzzzz_1, g_xxxxzz_0_xxxyyyyz_1, g_xxxxzz_0_xxxyyyz_1, g_xxxxzz_0_xxxyyyzz_1, g_xxxxzz_0_xxxyyzz_1, g_xxxxzz_0_xxxyyzzz_1, g_xxxxzz_0_xxxyzzz_1, g_xxxxzz_0_xxxyzzzz_1, g_xxxxzz_0_xxxzzzz_1, g_xxxxzz_0_xxxzzzzz_1, g_xxxxzz_0_xxyyyyyz_1, g_xxxxzz_0_xxyyyyz_1, g_xxxxzz_0_xxyyyyzz_1, g_xxxxzz_0_xxyyyzz_1, g_xxxxzz_0_xxyyyzzz_1, g_xxxxzz_0_xxyyzzz_1, g_xxxxzz_0_xxyyzzzz_1, g_xxxxzz_0_xxyzzzz_1, g_xxxxzz_0_xxyzzzzz_1, g_xxxxzz_0_xxzzzzz_1, g_xxxxzz_0_xxzzzzzz_1, g_xxxxzz_0_xyyyyyyz_1, g_xxxxzz_0_xyyyyyz_1, g_xxxxzz_0_xyyyyyzz_1, g_xxxxzz_0_xyyyyzz_1, g_xxxxzz_0_xyyyyzzz_1, g_xxxxzz_0_xyyyzzz_1, g_xxxxzz_0_xyyyzzzz_1, g_xxxxzz_0_xyyzzzz_1, g_xxxxzz_0_xyyzzzzz_1, g_xxxxzz_0_xyzzzzz_1, g_xxxxzz_0_xyzzzzzz_1, g_xxxxzz_0_xzzzzzz_1, g_xxxxzz_0_xzzzzzzz_1, g_xxxxzz_0_yyyyyyyy_1, g_xxxxzz_0_yyyyyyyz_1, g_xxxxzz_0_yyyyyyz_1, g_xxxxzz_0_yyyyyyzz_1, g_xxxxzz_0_yyyyyzz_1, g_xxxxzz_0_yyyyyzzz_1, g_xxxxzz_0_yyyyzzz_1, g_xxxxzz_0_yyyyzzzz_1, g_xxxxzz_0_yyyzzzz_1, g_xxxxzz_0_yyyzzzzz_1, g_xxxxzz_0_yyzzzzz_1, g_xxxxzz_0_yyzzzzzz_1, g_xxxxzz_0_yzzzzzz_1, g_xxxxzz_0_yzzzzzzz_1, g_xxxxzz_0_zzzzzzz_1, g_xxxxzz_0_zzzzzzzz_1, g_xxxzz_0_xxxxxxxz_0, g_xxxzz_0_xxxxxxxz_1, g_xxxzz_0_xxxxxxyz_0, g_xxxzz_0_xxxxxxyz_1, g_xxxzz_0_xxxxxxzz_0, g_xxxzz_0_xxxxxxzz_1, g_xxxzz_0_xxxxxyyz_0, g_xxxzz_0_xxxxxyyz_1, g_xxxzz_0_xxxxxyzz_0, g_xxxzz_0_xxxxxyzz_1, g_xxxzz_0_xxxxxzzz_0, g_xxxzz_0_xxxxxzzz_1, g_xxxzz_0_xxxxyyyz_0, g_xxxzz_0_xxxxyyyz_1, g_xxxzz_0_xxxxyyzz_0, g_xxxzz_0_xxxxyyzz_1, g_xxxzz_0_xxxxyzzz_0, g_xxxzz_0_xxxxyzzz_1, g_xxxzz_0_xxxxzzzz_0, g_xxxzz_0_xxxxzzzz_1, g_xxxzz_0_xxxyyyyz_0, g_xxxzz_0_xxxyyyyz_1, g_xxxzz_0_xxxyyyzz_0, g_xxxzz_0_xxxyyyzz_1, g_xxxzz_0_xxxyyzzz_0, g_xxxzz_0_xxxyyzzz_1, g_xxxzz_0_xxxyzzzz_0, g_xxxzz_0_xxxyzzzz_1, g_xxxzz_0_xxxzzzzz_0, g_xxxzz_0_xxxzzzzz_1, g_xxxzz_0_xxyyyyyz_0, g_xxxzz_0_xxyyyyyz_1, g_xxxzz_0_xxyyyyzz_0, g_xxxzz_0_xxyyyyzz_1, g_xxxzz_0_xxyyyzzz_0, g_xxxzz_0_xxyyyzzz_1, g_xxxzz_0_xxyyzzzz_0, g_xxxzz_0_xxyyzzzz_1, g_xxxzz_0_xxyzzzzz_0, g_xxxzz_0_xxyzzzzz_1, g_xxxzz_0_xxzzzzzz_0, g_xxxzz_0_xxzzzzzz_1, g_xxxzz_0_xyyyyyyz_0, g_xxxzz_0_xyyyyyyz_1, g_xxxzz_0_xyyyyyzz_0, g_xxxzz_0_xyyyyyzz_1, g_xxxzz_0_xyyyyzzz_0, g_xxxzz_0_xyyyyzzz_1, g_xxxzz_0_xyyyzzzz_0, g_xxxzz_0_xyyyzzzz_1, g_xxxzz_0_xyyzzzzz_0, g_xxxzz_0_xyyzzzzz_1, g_xxxzz_0_xyzzzzzz_0, g_xxxzz_0_xyzzzzzz_1, g_xxxzz_0_xzzzzzzz_0, g_xxxzz_0_xzzzzzzz_1, g_xxxzz_0_yyyyyyyy_0, g_xxxzz_0_yyyyyyyy_1, g_xxxzz_0_yyyyyyyz_0, g_xxxzz_0_yyyyyyyz_1, g_xxxzz_0_yyyyyyzz_0, g_xxxzz_0_yyyyyyzz_1, g_xxxzz_0_yyyyyzzz_0, g_xxxzz_0_yyyyyzzz_1, g_xxxzz_0_yyyyzzzz_0, g_xxxzz_0_yyyyzzzz_1, g_xxxzz_0_yyyzzzzz_0, g_xxxzz_0_yyyzzzzz_1, g_xxxzz_0_yyzzzzzz_0, g_xxxzz_0_yyzzzzzz_1, g_xxxzz_0_yzzzzzzz_0, g_xxxzz_0_yzzzzzzz_1, g_xxxzz_0_zzzzzzzz_0, g_xxxzz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxzz_0_xxxxxxxx_0[i] = g_xxxxx_0_xxxxxxxx_0[i] * fbe_0 - g_xxxxx_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxxxz_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxxxxy_0[i] = g_xxxxx_0_xxxxxxxy_0[i] * fbe_0 - g_xxxxx_0_xxxxxxxy_1[i] * fz_be_0 + g_xxxxxz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxxxxz_0[i] = 4.0 * g_xxxzz_0_xxxxxxxz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xxxxzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxxxyy_0[i] = g_xxxxx_0_xxxxxxyy_0[i] * fbe_0 - g_xxxxx_0_xxxxxxyy_1[i] * fz_be_0 + g_xxxxxz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxxxyz_0[i] = 4.0 * g_xxxzz_0_xxxxxxyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxxxzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxxxzz_0[i] = 4.0 * g_xxxzz_0_xxxxxxzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xxxxzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxxyyy_0[i] = g_xxxxx_0_xxxxxyyy_0[i] * fbe_0 - g_xxxxx_0_xxxxxyyy_1[i] * fz_be_0 + g_xxxxxz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxxyyz_0[i] = 4.0 * g_xxxzz_0_xxxxxyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxxxzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxxyzz_0[i] = 4.0 * g_xxxzz_0_xxxxxyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxxxzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxxzzz_0[i] = 4.0 * g_xxxzz_0_xxxxxzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xxxxzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxyyyy_0[i] = g_xxxxx_0_xxxxyyyy_0[i] * fbe_0 - g_xxxxx_0_xxxxyyyy_1[i] * fz_be_0 + g_xxxxxz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxyyyz_0[i] = 4.0 * g_xxxzz_0_xxxxyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxyyzz_0[i] = 4.0 * g_xxxzz_0_xxxxyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxyzzz_0[i] = 4.0 * g_xxxzz_0_xxxxyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxzzzz_0[i] = 4.0 * g_xxxzz_0_xxxxzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxyyyyy_0[i] = g_xxxxx_0_xxxyyyyy_0[i] * fbe_0 - g_xxxxx_0_xxxyyyyy_1[i] * fz_be_0 + g_xxxxxz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxyyyyz_0[i] = 4.0 * g_xxxzz_0_xxxyyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxyyyzz_0[i] = 4.0 * g_xxxzz_0_xxxyyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxyyzzz_0[i] = 4.0 * g_xxxzz_0_xxxyyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxyzzzz_0[i] = 4.0 * g_xxxzz_0_xxxyzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxzzzzz_0[i] = 4.0 * g_xxxzz_0_xxxzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyyyyyy_0[i] = g_xxxxx_0_xxyyyyyy_0[i] * fbe_0 - g_xxxxx_0_xxyyyyyy_1[i] * fz_be_0 + g_xxxxxz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxyyyyyz_0[i] = 4.0 * g_xxxzz_0_xxyyyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyyyyzz_0[i] = 4.0 * g_xxxzz_0_xxyyyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyyyzzz_0[i] = 4.0 * g_xxxzz_0_xxyyyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyyzzzz_0[i] = 4.0 * g_xxxzz_0_xxyyzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyzzzzz_0[i] = 4.0 * g_xxxzz_0_xxyzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxzzzzzz_0[i] = 4.0 * g_xxxzz_0_xxzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyyyyyy_0[i] = g_xxxxx_0_xyyyyyyy_0[i] * fbe_0 - g_xxxxx_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxxxz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xyyyyyyz_0[i] = 4.0 * g_xxxzz_0_xyyyyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyyyyzz_0[i] = 4.0 * g_xxxzz_0_xyyyyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyyyzzz_0[i] = 4.0 * g_xxxzz_0_xyyyyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyyzzzz_0[i] = 4.0 * g_xxxzz_0_xyyyzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyzzzzz_0[i] = 4.0 * g_xxxzz_0_xyyzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyzzzzzz_0[i] = 4.0 * g_xxxzz_0_xyzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xzzzzzzz_0[i] = 4.0 * g_xxxzz_0_xzzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyyyyy_0[i] = 4.0 * g_xxxzz_0_yyyyyyyy_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xxxxzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyyyyz_0[i] = 4.0 * g_xxxzz_0_yyyyyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyyyzz_0[i] = 4.0 * g_xxxzz_0_yyyyyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyyzzz_0[i] = 4.0 * g_xxxzz_0_yyyyyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyzzzz_0[i] = 4.0 * g_xxxzz_0_yyyyzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyzzzzz_0[i] = 4.0 * g_xxxzz_0_yyyzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyzzzzzz_0[i] = 4.0 * g_xxxzz_0_yyzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yzzzzzzz_0[i] = 4.0 * g_xxxzz_0_yzzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_zzzzzzzz_0[i] = 4.0 * g_xxxzz_0_zzzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 270-315 components of targeted buffer : KSL

    auto g_xxxxyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 270);

    auto g_xxxxyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 271);

    auto g_xxxxyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 272);

    auto g_xxxxyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 273);

    auto g_xxxxyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 274);

    auto g_xxxxyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 275);

    auto g_xxxxyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 276);

    auto g_xxxxyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 277);

    auto g_xxxxyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 278);

    auto g_xxxxyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 279);

    auto g_xxxxyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 280);

    auto g_xxxxyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 281);

    auto g_xxxxyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 282);

    auto g_xxxxyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 283);

    auto g_xxxxyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 284);

    auto g_xxxxyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 285);

    auto g_xxxxyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 286);

    auto g_xxxxyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 287);

    auto g_xxxxyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 288);

    auto g_xxxxyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 289);

    auto g_xxxxyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 290);

    auto g_xxxxyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 291);

    auto g_xxxxyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 292);

    auto g_xxxxyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 293);

    auto g_xxxxyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 294);

    auto g_xxxxyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 295);

    auto g_xxxxyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 296);

    auto g_xxxxyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 297);

    auto g_xxxxyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 298);

    auto g_xxxxyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 299);

    auto g_xxxxyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 300);

    auto g_xxxxyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 301);

    auto g_xxxxyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 302);

    auto g_xxxxyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 303);

    auto g_xxxxyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 304);

    auto g_xxxxyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 305);

    auto g_xxxxyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 306);

    auto g_xxxxyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 307);

    auto g_xxxxyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 308);

    auto g_xxxxyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 309);

    auto g_xxxxyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 310);

    auto g_xxxxyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 311);

    auto g_xxxxyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 312);

    auto g_xxxxyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 313);

    auto g_xxxxyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 314);

    #pragma omp simd aligned(g_xxxxy_0_xxxxxxxx_0, g_xxxxy_0_xxxxxxxx_1, g_xxxxy_0_xxxxxxxz_0, g_xxxxy_0_xxxxxxxz_1, g_xxxxy_0_xxxxxxzz_0, g_xxxxy_0_xxxxxxzz_1, g_xxxxy_0_xxxxxzzz_0, g_xxxxy_0_xxxxxzzz_1, g_xxxxy_0_xxxxzzzz_0, g_xxxxy_0_xxxxzzzz_1, g_xxxxy_0_xxxzzzzz_0, g_xxxxy_0_xxxzzzzz_1, g_xxxxy_0_xxzzzzzz_0, g_xxxxy_0_xxzzzzzz_1, g_xxxxy_0_xzzzzzzz_0, g_xxxxy_0_xzzzzzzz_1, g_xxxxyy_0_xxxxxxxx_1, g_xxxxyy_0_xxxxxxxz_1, g_xxxxyy_0_xxxxxxzz_1, g_xxxxyy_0_xxxxxzzz_1, g_xxxxyy_0_xxxxzzzz_1, g_xxxxyy_0_xxxzzzzz_1, g_xxxxyy_0_xxzzzzzz_1, g_xxxxyy_0_xzzzzzzz_1, g_xxxxyyy_0_xxxxxxxx_0, g_xxxxyyy_0_xxxxxxxy_0, g_xxxxyyy_0_xxxxxxxz_0, g_xxxxyyy_0_xxxxxxyy_0, g_xxxxyyy_0_xxxxxxyz_0, g_xxxxyyy_0_xxxxxxzz_0, g_xxxxyyy_0_xxxxxyyy_0, g_xxxxyyy_0_xxxxxyyz_0, g_xxxxyyy_0_xxxxxyzz_0, g_xxxxyyy_0_xxxxxzzz_0, g_xxxxyyy_0_xxxxyyyy_0, g_xxxxyyy_0_xxxxyyyz_0, g_xxxxyyy_0_xxxxyyzz_0, g_xxxxyyy_0_xxxxyzzz_0, g_xxxxyyy_0_xxxxzzzz_0, g_xxxxyyy_0_xxxyyyyy_0, g_xxxxyyy_0_xxxyyyyz_0, g_xxxxyyy_0_xxxyyyzz_0, g_xxxxyyy_0_xxxyyzzz_0, g_xxxxyyy_0_xxxyzzzz_0, g_xxxxyyy_0_xxxzzzzz_0, g_xxxxyyy_0_xxyyyyyy_0, g_xxxxyyy_0_xxyyyyyz_0, g_xxxxyyy_0_xxyyyyzz_0, g_xxxxyyy_0_xxyyyzzz_0, g_xxxxyyy_0_xxyyzzzz_0, g_xxxxyyy_0_xxyzzzzz_0, g_xxxxyyy_0_xxzzzzzz_0, g_xxxxyyy_0_xyyyyyyy_0, g_xxxxyyy_0_xyyyyyyz_0, g_xxxxyyy_0_xyyyyyzz_0, g_xxxxyyy_0_xyyyyzzz_0, g_xxxxyyy_0_xyyyzzzz_0, g_xxxxyyy_0_xyyzzzzz_0, g_xxxxyyy_0_xyzzzzzz_0, g_xxxxyyy_0_xzzzzzzz_0, g_xxxxyyy_0_yyyyyyyy_0, g_xxxxyyy_0_yyyyyyyz_0, g_xxxxyyy_0_yyyyyyzz_0, g_xxxxyyy_0_yyyyyzzz_0, g_xxxxyyy_0_yyyyzzzz_0, g_xxxxyyy_0_yyyzzzzz_0, g_xxxxyyy_0_yyzzzzzz_0, g_xxxxyyy_0_yzzzzzzz_0, g_xxxxyyy_0_zzzzzzzz_0, g_xxxyyy_0_xxxxxxxy_1, g_xxxyyy_0_xxxxxxy_1, g_xxxyyy_0_xxxxxxyy_1, g_xxxyyy_0_xxxxxxyz_1, g_xxxyyy_0_xxxxxyy_1, g_xxxyyy_0_xxxxxyyy_1, g_xxxyyy_0_xxxxxyyz_1, g_xxxyyy_0_xxxxxyz_1, g_xxxyyy_0_xxxxxyzz_1, g_xxxyyy_0_xxxxyyy_1, g_xxxyyy_0_xxxxyyyy_1, g_xxxyyy_0_xxxxyyyz_1, g_xxxyyy_0_xxxxyyz_1, g_xxxyyy_0_xxxxyyzz_1, g_xxxyyy_0_xxxxyzz_1, g_xxxyyy_0_xxxxyzzz_1, g_xxxyyy_0_xxxyyyy_1, g_xxxyyy_0_xxxyyyyy_1, g_xxxyyy_0_xxxyyyyz_1, g_xxxyyy_0_xxxyyyz_1, g_xxxyyy_0_xxxyyyzz_1, g_xxxyyy_0_xxxyyzz_1, g_xxxyyy_0_xxxyyzzz_1, g_xxxyyy_0_xxxyzzz_1, g_xxxyyy_0_xxxyzzzz_1, g_xxxyyy_0_xxyyyyy_1, g_xxxyyy_0_xxyyyyyy_1, g_xxxyyy_0_xxyyyyyz_1, g_xxxyyy_0_xxyyyyz_1, g_xxxyyy_0_xxyyyyzz_1, g_xxxyyy_0_xxyyyzz_1, g_xxxyyy_0_xxyyyzzz_1, g_xxxyyy_0_xxyyzzz_1, g_xxxyyy_0_xxyyzzzz_1, g_xxxyyy_0_xxyzzzz_1, g_xxxyyy_0_xxyzzzzz_1, g_xxxyyy_0_xyyyyyy_1, g_xxxyyy_0_xyyyyyyy_1, g_xxxyyy_0_xyyyyyyz_1, g_xxxyyy_0_xyyyyyz_1, g_xxxyyy_0_xyyyyyzz_1, g_xxxyyy_0_xyyyyzz_1, g_xxxyyy_0_xyyyyzzz_1, g_xxxyyy_0_xyyyzzz_1, g_xxxyyy_0_xyyyzzzz_1, g_xxxyyy_0_xyyzzzz_1, g_xxxyyy_0_xyyzzzzz_1, g_xxxyyy_0_xyzzzzz_1, g_xxxyyy_0_xyzzzzzz_1, g_xxxyyy_0_yyyyyyy_1, g_xxxyyy_0_yyyyyyyy_1, g_xxxyyy_0_yyyyyyyz_1, g_xxxyyy_0_yyyyyyz_1, g_xxxyyy_0_yyyyyyzz_1, g_xxxyyy_0_yyyyyzz_1, g_xxxyyy_0_yyyyyzzz_1, g_xxxyyy_0_yyyyzzz_1, g_xxxyyy_0_yyyyzzzz_1, g_xxxyyy_0_yyyzzzz_1, g_xxxyyy_0_yyyzzzzz_1, g_xxxyyy_0_yyzzzzz_1, g_xxxyyy_0_yyzzzzzz_1, g_xxxyyy_0_yzzzzzz_1, g_xxxyyy_0_yzzzzzzz_1, g_xxxyyy_0_zzzzzzzz_1, g_xxyyy_0_xxxxxxxy_0, g_xxyyy_0_xxxxxxxy_1, g_xxyyy_0_xxxxxxyy_0, g_xxyyy_0_xxxxxxyy_1, g_xxyyy_0_xxxxxxyz_0, g_xxyyy_0_xxxxxxyz_1, g_xxyyy_0_xxxxxyyy_0, g_xxyyy_0_xxxxxyyy_1, g_xxyyy_0_xxxxxyyz_0, g_xxyyy_0_xxxxxyyz_1, g_xxyyy_0_xxxxxyzz_0, g_xxyyy_0_xxxxxyzz_1, g_xxyyy_0_xxxxyyyy_0, g_xxyyy_0_xxxxyyyy_1, g_xxyyy_0_xxxxyyyz_0, g_xxyyy_0_xxxxyyyz_1, g_xxyyy_0_xxxxyyzz_0, g_xxyyy_0_xxxxyyzz_1, g_xxyyy_0_xxxxyzzz_0, g_xxyyy_0_xxxxyzzz_1, g_xxyyy_0_xxxyyyyy_0, g_xxyyy_0_xxxyyyyy_1, g_xxyyy_0_xxxyyyyz_0, g_xxyyy_0_xxxyyyyz_1, g_xxyyy_0_xxxyyyzz_0, g_xxyyy_0_xxxyyyzz_1, g_xxyyy_0_xxxyyzzz_0, g_xxyyy_0_xxxyyzzz_1, g_xxyyy_0_xxxyzzzz_0, g_xxyyy_0_xxxyzzzz_1, g_xxyyy_0_xxyyyyyy_0, g_xxyyy_0_xxyyyyyy_1, g_xxyyy_0_xxyyyyyz_0, g_xxyyy_0_xxyyyyyz_1, g_xxyyy_0_xxyyyyzz_0, g_xxyyy_0_xxyyyyzz_1, g_xxyyy_0_xxyyyzzz_0, g_xxyyy_0_xxyyyzzz_1, g_xxyyy_0_xxyyzzzz_0, g_xxyyy_0_xxyyzzzz_1, g_xxyyy_0_xxyzzzzz_0, g_xxyyy_0_xxyzzzzz_1, g_xxyyy_0_xyyyyyyy_0, g_xxyyy_0_xyyyyyyy_1, g_xxyyy_0_xyyyyyyz_0, g_xxyyy_0_xyyyyyyz_1, g_xxyyy_0_xyyyyyzz_0, g_xxyyy_0_xyyyyyzz_1, g_xxyyy_0_xyyyyzzz_0, g_xxyyy_0_xyyyyzzz_1, g_xxyyy_0_xyyyzzzz_0, g_xxyyy_0_xyyyzzzz_1, g_xxyyy_0_xyyzzzzz_0, g_xxyyy_0_xyyzzzzz_1, g_xxyyy_0_xyzzzzzz_0, g_xxyyy_0_xyzzzzzz_1, g_xxyyy_0_yyyyyyyy_0, g_xxyyy_0_yyyyyyyy_1, g_xxyyy_0_yyyyyyyz_0, g_xxyyy_0_yyyyyyyz_1, g_xxyyy_0_yyyyyyzz_0, g_xxyyy_0_yyyyyyzz_1, g_xxyyy_0_yyyyyzzz_0, g_xxyyy_0_yyyyyzzz_1, g_xxyyy_0_yyyyzzzz_0, g_xxyyy_0_yyyyzzzz_1, g_xxyyy_0_yyyzzzzz_0, g_xxyyy_0_yyyzzzzz_1, g_xxyyy_0_yyzzzzzz_0, g_xxyyy_0_yyzzzzzz_1, g_xxyyy_0_yzzzzzzz_0, g_xxyyy_0_yzzzzzzz_1, g_xxyyy_0_zzzzzzzz_0, g_xxyyy_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyyy_0_xxxxxxxx_0[i] = 2.0 * g_xxxxy_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxxyy_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxxxxxy_0[i] = 3.0 * g_xxyyy_0_xxxxxxxy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xxxyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxxxxz_0[i] = 2.0 * g_xxxxy_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxxxxz_1[i] * fz_be_0 + g_xxxxyy_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxxxxyy_0[i] = 3.0 * g_xxyyy_0_xxxxxxyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xxxyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxxxyz_0[i] = 3.0 * g_xxyyy_0_xxxxxxyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxxyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxxxzz_0[i] = 2.0 * g_xxxxy_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxxxzz_1[i] * fz_be_0 + g_xxxxyy_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxxxyyy_0[i] = 3.0 * g_xxyyy_0_xxxxxyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xxxyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxxyyz_0[i] = 3.0 * g_xxyyy_0_xxxxxyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxxyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxxyzz_0[i] = 3.0 * g_xxyyy_0_xxxxxyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxxyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxxzzz_0[i] = 2.0 * g_xxxxy_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxxzzz_1[i] * fz_be_0 + g_xxxxyy_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxxyyyy_0[i] = 3.0 * g_xxyyy_0_xxxxyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xxxyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxyyyz_0[i] = 3.0 * g_xxyyy_0_xxxxyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxxyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxyyzz_0[i] = 3.0 * g_xxyyy_0_xxxxyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxxyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxyzzz_0[i] = 3.0 * g_xxyyy_0_xxxxyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxxyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxzzzz_0[i] = 2.0 * g_xxxxy_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxzzzz_1[i] * fz_be_0 + g_xxxxyy_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxyyyyy_0[i] = 3.0 * g_xxyyy_0_xxxyyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxyyyyz_0[i] = 3.0 * g_xxyyy_0_xxxyyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxyyyzz_0[i] = 3.0 * g_xxyyy_0_xxxyyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxyyzzz_0[i] = 3.0 * g_xxyyy_0_xxxyyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxyzzzz_0[i] = 3.0 * g_xxyyy_0_xxxyzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxzzzzz_0[i] = 2.0 * g_xxxxy_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxyyyyyy_0[i] = 3.0 * g_xxyyy_0_xxyyyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyyyyyz_0[i] = 3.0 * g_xxyyy_0_xxyyyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyyyyzz_0[i] = 3.0 * g_xxyyy_0_xxyyyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyyyzzz_0[i] = 3.0 * g_xxyyy_0_xxyyyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyyzzzz_0[i] = 3.0 * g_xxyyy_0_xxyyzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyzzzzz_0[i] = 3.0 * g_xxyyy_0_xxyzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxzzzzzz_0[i] = 2.0 * g_xxxxy_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxzzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xyyyyyyy_0[i] = 3.0 * g_xxyyy_0_xyyyyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyyyyyz_0[i] = 3.0 * g_xxyyy_0_xyyyyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyyyyz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyyyyzz_0[i] = 3.0 * g_xxyyy_0_xyyyyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyyyzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyyyzzz_0[i] = 3.0 * g_xxyyy_0_xyyyyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyyzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyyzzzz_0[i] = 3.0 * g_xxyyy_0_xyyyzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyzzzzz_0[i] = 3.0 * g_xxyyy_0_xyyzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyzzzzzz_0[i] = 3.0 * g_xxyyy_0_xyzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyzzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xzzzzzzz_0[i] = 2.0 * g_xxxxy_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_yyyyyyyy_0[i] = 3.0 * g_xxyyy_0_yyyyyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyyyyy_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyyyyyz_0[i] = 3.0 * g_xxyyy_0_yyyyyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyyyyz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyyyyzz_0[i] = 3.0 * g_xxyyy_0_yyyyyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyyyzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyyyzzz_0[i] = 3.0 * g_xxyyy_0_yyyyyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyyzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyyzzzz_0[i] = 3.0 * g_xxyyy_0_yyyyzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyzzzzz_0[i] = 3.0 * g_xxyyy_0_yyyzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyzzzzzz_0[i] = 3.0 * g_xxyyy_0_yyzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyzzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yzzzzzzz_0[i] = 3.0 * g_xxyyy_0_yzzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yzzzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_zzzzzzzz_0[i] = 3.0 * g_xxyyy_0_zzzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_zzzzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 315-360 components of targeted buffer : KSL

    auto g_xxxxyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 315);

    auto g_xxxxyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 316);

    auto g_xxxxyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 317);

    auto g_xxxxyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 318);

    auto g_xxxxyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 319);

    auto g_xxxxyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 320);

    auto g_xxxxyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 321);

    auto g_xxxxyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 322);

    auto g_xxxxyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 323);

    auto g_xxxxyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 324);

    auto g_xxxxyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 325);

    auto g_xxxxyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 326);

    auto g_xxxxyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 327);

    auto g_xxxxyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 328);

    auto g_xxxxyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 329);

    auto g_xxxxyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 330);

    auto g_xxxxyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 331);

    auto g_xxxxyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 332);

    auto g_xxxxyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 333);

    auto g_xxxxyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 334);

    auto g_xxxxyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 335);

    auto g_xxxxyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 336);

    auto g_xxxxyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 337);

    auto g_xxxxyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 338);

    auto g_xxxxyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 339);

    auto g_xxxxyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 340);

    auto g_xxxxyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 341);

    auto g_xxxxyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 342);

    auto g_xxxxyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 343);

    auto g_xxxxyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 344);

    auto g_xxxxyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 345);

    auto g_xxxxyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 346);

    auto g_xxxxyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 347);

    auto g_xxxxyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 348);

    auto g_xxxxyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 349);

    auto g_xxxxyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 350);

    auto g_xxxxyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 351);

    auto g_xxxxyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 352);

    auto g_xxxxyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 353);

    auto g_xxxxyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 354);

    auto g_xxxxyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 355);

    auto g_xxxxyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 356);

    auto g_xxxxyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 357);

    auto g_xxxxyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 358);

    auto g_xxxxyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 359);

    #pragma omp simd aligned(g_xxxxyy_0_xxxxxxx_1, g_xxxxyy_0_xxxxxxxx_1, g_xxxxyy_0_xxxxxxxy_1, g_xxxxyy_0_xxxxxxxz_1, g_xxxxyy_0_xxxxxxy_1, g_xxxxyy_0_xxxxxxyy_1, g_xxxxyy_0_xxxxxxyz_1, g_xxxxyy_0_xxxxxxz_1, g_xxxxyy_0_xxxxxxzz_1, g_xxxxyy_0_xxxxxyy_1, g_xxxxyy_0_xxxxxyyy_1, g_xxxxyy_0_xxxxxyyz_1, g_xxxxyy_0_xxxxxyz_1, g_xxxxyy_0_xxxxxyzz_1, g_xxxxyy_0_xxxxxzz_1, g_xxxxyy_0_xxxxxzzz_1, g_xxxxyy_0_xxxxyyy_1, g_xxxxyy_0_xxxxyyyy_1, g_xxxxyy_0_xxxxyyyz_1, g_xxxxyy_0_xxxxyyz_1, g_xxxxyy_0_xxxxyyzz_1, g_xxxxyy_0_xxxxyzz_1, g_xxxxyy_0_xxxxyzzz_1, g_xxxxyy_0_xxxxzzz_1, g_xxxxyy_0_xxxxzzzz_1, g_xxxxyy_0_xxxyyyy_1, g_xxxxyy_0_xxxyyyyy_1, g_xxxxyy_0_xxxyyyyz_1, g_xxxxyy_0_xxxyyyz_1, g_xxxxyy_0_xxxyyyzz_1, g_xxxxyy_0_xxxyyzz_1, g_xxxxyy_0_xxxyyzzz_1, g_xxxxyy_0_xxxyzzz_1, g_xxxxyy_0_xxxyzzzz_1, g_xxxxyy_0_xxxzzzz_1, g_xxxxyy_0_xxxzzzzz_1, g_xxxxyy_0_xxyyyyy_1, g_xxxxyy_0_xxyyyyyy_1, g_xxxxyy_0_xxyyyyyz_1, g_xxxxyy_0_xxyyyyz_1, g_xxxxyy_0_xxyyyyzz_1, g_xxxxyy_0_xxyyyzz_1, g_xxxxyy_0_xxyyyzzz_1, g_xxxxyy_0_xxyyzzz_1, g_xxxxyy_0_xxyyzzzz_1, g_xxxxyy_0_xxyzzzz_1, g_xxxxyy_0_xxyzzzzz_1, g_xxxxyy_0_xxzzzzz_1, g_xxxxyy_0_xxzzzzzz_1, g_xxxxyy_0_xyyyyyy_1, g_xxxxyy_0_xyyyyyyy_1, g_xxxxyy_0_xyyyyyyz_1, g_xxxxyy_0_xyyyyyz_1, g_xxxxyy_0_xyyyyyzz_1, g_xxxxyy_0_xyyyyzz_1, g_xxxxyy_0_xyyyyzzz_1, g_xxxxyy_0_xyyyzzz_1, g_xxxxyy_0_xyyyzzzz_1, g_xxxxyy_0_xyyzzzz_1, g_xxxxyy_0_xyyzzzzz_1, g_xxxxyy_0_xyzzzzz_1, g_xxxxyy_0_xyzzzzzz_1, g_xxxxyy_0_xzzzzzz_1, g_xxxxyy_0_xzzzzzzz_1, g_xxxxyy_0_yyyyyyy_1, g_xxxxyy_0_yyyyyyyy_1, g_xxxxyy_0_yyyyyyyz_1, g_xxxxyy_0_yyyyyyz_1, g_xxxxyy_0_yyyyyyzz_1, g_xxxxyy_0_yyyyyzz_1, g_xxxxyy_0_yyyyyzzz_1, g_xxxxyy_0_yyyyzzz_1, g_xxxxyy_0_yyyyzzzz_1, g_xxxxyy_0_yyyzzzz_1, g_xxxxyy_0_yyyzzzzz_1, g_xxxxyy_0_yyzzzzz_1, g_xxxxyy_0_yyzzzzzz_1, g_xxxxyy_0_yzzzzzz_1, g_xxxxyy_0_yzzzzzzz_1, g_xxxxyy_0_zzzzzzz_1, g_xxxxyy_0_zzzzzzzz_1, g_xxxxyyz_0_xxxxxxxx_0, g_xxxxyyz_0_xxxxxxxy_0, g_xxxxyyz_0_xxxxxxxz_0, g_xxxxyyz_0_xxxxxxyy_0, g_xxxxyyz_0_xxxxxxyz_0, g_xxxxyyz_0_xxxxxxzz_0, g_xxxxyyz_0_xxxxxyyy_0, g_xxxxyyz_0_xxxxxyyz_0, g_xxxxyyz_0_xxxxxyzz_0, g_xxxxyyz_0_xxxxxzzz_0, g_xxxxyyz_0_xxxxyyyy_0, g_xxxxyyz_0_xxxxyyyz_0, g_xxxxyyz_0_xxxxyyzz_0, g_xxxxyyz_0_xxxxyzzz_0, g_xxxxyyz_0_xxxxzzzz_0, g_xxxxyyz_0_xxxyyyyy_0, g_xxxxyyz_0_xxxyyyyz_0, g_xxxxyyz_0_xxxyyyzz_0, g_xxxxyyz_0_xxxyyzzz_0, g_xxxxyyz_0_xxxyzzzz_0, g_xxxxyyz_0_xxxzzzzz_0, g_xxxxyyz_0_xxyyyyyy_0, g_xxxxyyz_0_xxyyyyyz_0, g_xxxxyyz_0_xxyyyyzz_0, g_xxxxyyz_0_xxyyyzzz_0, g_xxxxyyz_0_xxyyzzzz_0, g_xxxxyyz_0_xxyzzzzz_0, g_xxxxyyz_0_xxzzzzzz_0, g_xxxxyyz_0_xyyyyyyy_0, g_xxxxyyz_0_xyyyyyyz_0, g_xxxxyyz_0_xyyyyyzz_0, g_xxxxyyz_0_xyyyyzzz_0, g_xxxxyyz_0_xyyyzzzz_0, g_xxxxyyz_0_xyyzzzzz_0, g_xxxxyyz_0_xyzzzzzz_0, g_xxxxyyz_0_xzzzzzzz_0, g_xxxxyyz_0_yyyyyyyy_0, g_xxxxyyz_0_yyyyyyyz_0, g_xxxxyyz_0_yyyyyyzz_0, g_xxxxyyz_0_yyyyyzzz_0, g_xxxxyyz_0_yyyyzzzz_0, g_xxxxyyz_0_yyyzzzzz_0, g_xxxxyyz_0_yyzzzzzz_0, g_xxxxyyz_0_yzzzzzzz_0, g_xxxxyyz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyyz_0_xxxxxxxx_0[i] = g_xxxxyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxxxy_0[i] = g_xxxxyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxxxz_0[i] = g_xxxxyy_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxxxz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxxyy_0[i] = g_xxxxyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxxyz_0[i] = g_xxxxyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxxyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxxzz_0[i] = 2.0 * g_xxxxyy_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxxzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxyyy_0[i] = g_xxxxyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxyyz_0[i] = g_xxxxyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxyzz_0[i] = 2.0 * g_xxxxyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxzzz_0[i] = 3.0 * g_xxxxyy_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxyyyy_0[i] = g_xxxxyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxyyyz_0[i] = g_xxxxyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxyyzz_0[i] = 2.0 * g_xxxxyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxyzzz_0[i] = 3.0 * g_xxxxyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxzzzz_0[i] = 4.0 * g_xxxxyy_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyyyyy_0[i] = g_xxxxyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyyyyz_0[i] = g_xxxxyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyyyzz_0[i] = 2.0 * g_xxxxyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyyzzz_0[i] = 3.0 * g_xxxxyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyzzzz_0[i] = 4.0 * g_xxxxyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxzzzzz_0[i] = 5.0 * g_xxxxyy_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyyyyy_0[i] = g_xxxxyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyyyyz_0[i] = g_xxxxyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyyyzz_0[i] = 2.0 * g_xxxxyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyyzzz_0[i] = 3.0 * g_xxxxyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyzzzz_0[i] = 4.0 * g_xxxxyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyzzzzz_0[i] = 5.0 * g_xxxxyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxzzzzzz_0[i] = 6.0 * g_xxxxyy_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxzzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyyyyy_0[i] = g_xxxxyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyyyyz_0[i] = g_xxxxyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyyyzz_0[i] = 2.0 * g_xxxxyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyyzzz_0[i] = 3.0 * g_xxxxyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyzzzz_0[i] = 4.0 * g_xxxxyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyzzzzz_0[i] = 5.0 * g_xxxxyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyzzzzzz_0[i] = 6.0 * g_xxxxyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyzzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xzzzzzzz_0[i] = 7.0 * g_xxxxyy_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xzzzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyyyyy_0[i] = g_xxxxyy_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyyyyz_0[i] = g_xxxxyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_yyyyyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyyyzz_0[i] = 2.0 * g_xxxxyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_yyyyyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyyzzz_0[i] = 3.0 * g_xxxxyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_yyyyyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyzzzz_0[i] = 4.0 * g_xxxxyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_yyyyzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyzzzzz_0[i] = 5.0 * g_xxxxyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_yyyzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyzzzzzz_0[i] = 6.0 * g_xxxxyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_yyzzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yzzzzzzz_0[i] = 7.0 * g_xxxxyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_yzzzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_zzzzzzzz_0[i] = 8.0 * g_xxxxyy_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 360-405 components of targeted buffer : KSL

    auto g_xxxxyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 360);

    auto g_xxxxyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 361);

    auto g_xxxxyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 362);

    auto g_xxxxyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 363);

    auto g_xxxxyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 364);

    auto g_xxxxyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 365);

    auto g_xxxxyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 366);

    auto g_xxxxyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 367);

    auto g_xxxxyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 368);

    auto g_xxxxyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 369);

    auto g_xxxxyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 370);

    auto g_xxxxyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 371);

    auto g_xxxxyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 372);

    auto g_xxxxyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 373);

    auto g_xxxxyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 374);

    auto g_xxxxyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 375);

    auto g_xxxxyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 376);

    auto g_xxxxyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 377);

    auto g_xxxxyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 378);

    auto g_xxxxyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 379);

    auto g_xxxxyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 380);

    auto g_xxxxyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 381);

    auto g_xxxxyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 382);

    auto g_xxxxyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 383);

    auto g_xxxxyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 384);

    auto g_xxxxyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 385);

    auto g_xxxxyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 386);

    auto g_xxxxyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 387);

    auto g_xxxxyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 388);

    auto g_xxxxyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 389);

    auto g_xxxxyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 390);

    auto g_xxxxyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 391);

    auto g_xxxxyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 392);

    auto g_xxxxyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 393);

    auto g_xxxxyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 394);

    auto g_xxxxyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 395);

    auto g_xxxxyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 396);

    auto g_xxxxyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 397);

    auto g_xxxxyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 398);

    auto g_xxxxyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 399);

    auto g_xxxxyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 400);

    auto g_xxxxyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 401);

    auto g_xxxxyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 402);

    auto g_xxxxyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 403);

    auto g_xxxxyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 404);

    #pragma omp simd aligned(g_xxxxyzz_0_xxxxxxxx_0, g_xxxxyzz_0_xxxxxxxy_0, g_xxxxyzz_0_xxxxxxxz_0, g_xxxxyzz_0_xxxxxxyy_0, g_xxxxyzz_0_xxxxxxyz_0, g_xxxxyzz_0_xxxxxxzz_0, g_xxxxyzz_0_xxxxxyyy_0, g_xxxxyzz_0_xxxxxyyz_0, g_xxxxyzz_0_xxxxxyzz_0, g_xxxxyzz_0_xxxxxzzz_0, g_xxxxyzz_0_xxxxyyyy_0, g_xxxxyzz_0_xxxxyyyz_0, g_xxxxyzz_0_xxxxyyzz_0, g_xxxxyzz_0_xxxxyzzz_0, g_xxxxyzz_0_xxxxzzzz_0, g_xxxxyzz_0_xxxyyyyy_0, g_xxxxyzz_0_xxxyyyyz_0, g_xxxxyzz_0_xxxyyyzz_0, g_xxxxyzz_0_xxxyyzzz_0, g_xxxxyzz_0_xxxyzzzz_0, g_xxxxyzz_0_xxxzzzzz_0, g_xxxxyzz_0_xxyyyyyy_0, g_xxxxyzz_0_xxyyyyyz_0, g_xxxxyzz_0_xxyyyyzz_0, g_xxxxyzz_0_xxyyyzzz_0, g_xxxxyzz_0_xxyyzzzz_0, g_xxxxyzz_0_xxyzzzzz_0, g_xxxxyzz_0_xxzzzzzz_0, g_xxxxyzz_0_xyyyyyyy_0, g_xxxxyzz_0_xyyyyyyz_0, g_xxxxyzz_0_xyyyyyzz_0, g_xxxxyzz_0_xyyyyzzz_0, g_xxxxyzz_0_xyyyzzzz_0, g_xxxxyzz_0_xyyzzzzz_0, g_xxxxyzz_0_xyzzzzzz_0, g_xxxxyzz_0_xzzzzzzz_0, g_xxxxyzz_0_yyyyyyyy_0, g_xxxxyzz_0_yyyyyyyz_0, g_xxxxyzz_0_yyyyyyzz_0, g_xxxxyzz_0_yyyyyzzz_0, g_xxxxyzz_0_yyyyzzzz_0, g_xxxxyzz_0_yyyzzzzz_0, g_xxxxyzz_0_yyzzzzzz_0, g_xxxxyzz_0_yzzzzzzz_0, g_xxxxyzz_0_zzzzzzzz_0, g_xxxxzz_0_xxxxxxx_1, g_xxxxzz_0_xxxxxxxx_1, g_xxxxzz_0_xxxxxxxy_1, g_xxxxzz_0_xxxxxxxz_1, g_xxxxzz_0_xxxxxxy_1, g_xxxxzz_0_xxxxxxyy_1, g_xxxxzz_0_xxxxxxyz_1, g_xxxxzz_0_xxxxxxz_1, g_xxxxzz_0_xxxxxxzz_1, g_xxxxzz_0_xxxxxyy_1, g_xxxxzz_0_xxxxxyyy_1, g_xxxxzz_0_xxxxxyyz_1, g_xxxxzz_0_xxxxxyz_1, g_xxxxzz_0_xxxxxyzz_1, g_xxxxzz_0_xxxxxzz_1, g_xxxxzz_0_xxxxxzzz_1, g_xxxxzz_0_xxxxyyy_1, g_xxxxzz_0_xxxxyyyy_1, g_xxxxzz_0_xxxxyyyz_1, g_xxxxzz_0_xxxxyyz_1, g_xxxxzz_0_xxxxyyzz_1, g_xxxxzz_0_xxxxyzz_1, g_xxxxzz_0_xxxxyzzz_1, g_xxxxzz_0_xxxxzzz_1, g_xxxxzz_0_xxxxzzzz_1, g_xxxxzz_0_xxxyyyy_1, g_xxxxzz_0_xxxyyyyy_1, g_xxxxzz_0_xxxyyyyz_1, g_xxxxzz_0_xxxyyyz_1, g_xxxxzz_0_xxxyyyzz_1, g_xxxxzz_0_xxxyyzz_1, g_xxxxzz_0_xxxyyzzz_1, g_xxxxzz_0_xxxyzzz_1, g_xxxxzz_0_xxxyzzzz_1, g_xxxxzz_0_xxxzzzz_1, g_xxxxzz_0_xxxzzzzz_1, g_xxxxzz_0_xxyyyyy_1, g_xxxxzz_0_xxyyyyyy_1, g_xxxxzz_0_xxyyyyyz_1, g_xxxxzz_0_xxyyyyz_1, g_xxxxzz_0_xxyyyyzz_1, g_xxxxzz_0_xxyyyzz_1, g_xxxxzz_0_xxyyyzzz_1, g_xxxxzz_0_xxyyzzz_1, g_xxxxzz_0_xxyyzzzz_1, g_xxxxzz_0_xxyzzzz_1, g_xxxxzz_0_xxyzzzzz_1, g_xxxxzz_0_xxzzzzz_1, g_xxxxzz_0_xxzzzzzz_1, g_xxxxzz_0_xyyyyyy_1, g_xxxxzz_0_xyyyyyyy_1, g_xxxxzz_0_xyyyyyyz_1, g_xxxxzz_0_xyyyyyz_1, g_xxxxzz_0_xyyyyyzz_1, g_xxxxzz_0_xyyyyzz_1, g_xxxxzz_0_xyyyyzzz_1, g_xxxxzz_0_xyyyzzz_1, g_xxxxzz_0_xyyyzzzz_1, g_xxxxzz_0_xyyzzzz_1, g_xxxxzz_0_xyyzzzzz_1, g_xxxxzz_0_xyzzzzz_1, g_xxxxzz_0_xyzzzzzz_1, g_xxxxzz_0_xzzzzzz_1, g_xxxxzz_0_xzzzzzzz_1, g_xxxxzz_0_yyyyyyy_1, g_xxxxzz_0_yyyyyyyy_1, g_xxxxzz_0_yyyyyyyz_1, g_xxxxzz_0_yyyyyyz_1, g_xxxxzz_0_yyyyyyzz_1, g_xxxxzz_0_yyyyyzz_1, g_xxxxzz_0_yyyyyzzz_1, g_xxxxzz_0_yyyyzzz_1, g_xxxxzz_0_yyyyzzzz_1, g_xxxxzz_0_yyyzzzz_1, g_xxxxzz_0_yyyzzzzz_1, g_xxxxzz_0_yyzzzzz_1, g_xxxxzz_0_yyzzzzzz_1, g_xxxxzz_0_yzzzzzz_1, g_xxxxzz_0_yzzzzzzz_1, g_xxxxzz_0_zzzzzzz_1, g_xxxxzz_0_zzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyzz_0_xxxxxxxx_0[i] = g_xxxxzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxxxy_0[i] = g_xxxxzz_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxxxy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxxxz_0[i] = g_xxxxzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxxyy_0[i] = 2.0 * g_xxxxzz_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxxyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxxyz_0[i] = g_xxxxzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxxzz_0[i] = g_xxxxzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxyyy_0[i] = 3.0 * g_xxxxzz_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxyyz_0[i] = 2.0 * g_xxxxzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxyzz_0[i] = g_xxxxzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxzzz_0[i] = g_xxxxzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxyyyy_0[i] = 4.0 * g_xxxxzz_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxyyyz_0[i] = 3.0 * g_xxxxzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxyyzz_0[i] = 2.0 * g_xxxxzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxyzzz_0[i] = g_xxxxzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxzzzz_0[i] = g_xxxxzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyyyyy_0[i] = 5.0 * g_xxxxzz_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyyyyz_0[i] = 4.0 * g_xxxxzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyyyzz_0[i] = 3.0 * g_xxxxzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyyzzz_0[i] = 2.0 * g_xxxxzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyzzzz_0[i] = g_xxxxzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxzzzzz_0[i] = g_xxxxzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyyyyy_0[i] = 6.0 * g_xxxxzz_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyyyyz_0[i] = 5.0 * g_xxxxzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyyyzz_0[i] = 4.0 * g_xxxxzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyyzzz_0[i] = 3.0 * g_xxxxzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyzzzz_0[i] = 2.0 * g_xxxxzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyzzzzz_0[i] = g_xxxxzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxzzzzzz_0[i] = g_xxxxzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyyyyy_0[i] = 7.0 * g_xxxxzz_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyyyyz_0[i] = 6.0 * g_xxxxzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyyyzz_0[i] = 5.0 * g_xxxxzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyyzzz_0[i] = 4.0 * g_xxxxzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyzzzz_0[i] = 3.0 * g_xxxxzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyzzzzz_0[i] = 2.0 * g_xxxxzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyzzzzzz_0[i] = g_xxxxzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xzzzzzzz_0[i] = g_xxxxzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyyyyy_0[i] = 8.0 * g_xxxxzz_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyyyyz_0[i] = 7.0 * g_xxxxzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyyyzz_0[i] = 6.0 * g_xxxxzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyyzzz_0[i] = 5.0 * g_xxxxzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyzzzz_0[i] = 4.0 * g_xxxxzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyzzzzz_0[i] = 3.0 * g_xxxxzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyzzzzzz_0[i] = 2.0 * g_xxxxzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yzzzzzzz_0[i] = g_xxxxzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_zzzzzzzz_0[i] = g_xxxxzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 405-450 components of targeted buffer : KSL

    auto g_xxxxzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 405);

    auto g_xxxxzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 406);

    auto g_xxxxzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 407);

    auto g_xxxxzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 408);

    auto g_xxxxzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 409);

    auto g_xxxxzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 410);

    auto g_xxxxzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 411);

    auto g_xxxxzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 412);

    auto g_xxxxzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 413);

    auto g_xxxxzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 414);

    auto g_xxxxzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 415);

    auto g_xxxxzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 416);

    auto g_xxxxzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 417);

    auto g_xxxxzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 418);

    auto g_xxxxzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 419);

    auto g_xxxxzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 420);

    auto g_xxxxzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 421);

    auto g_xxxxzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 422);

    auto g_xxxxzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 423);

    auto g_xxxxzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 424);

    auto g_xxxxzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 425);

    auto g_xxxxzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 426);

    auto g_xxxxzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 427);

    auto g_xxxxzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 428);

    auto g_xxxxzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 429);

    auto g_xxxxzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 430);

    auto g_xxxxzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 431);

    auto g_xxxxzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 432);

    auto g_xxxxzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 433);

    auto g_xxxxzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 434);

    auto g_xxxxzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 435);

    auto g_xxxxzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 436);

    auto g_xxxxzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 437);

    auto g_xxxxzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 438);

    auto g_xxxxzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 439);

    auto g_xxxxzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 440);

    auto g_xxxxzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 441);

    auto g_xxxxzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 442);

    auto g_xxxxzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 443);

    auto g_xxxxzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 444);

    auto g_xxxxzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 445);

    auto g_xxxxzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 446);

    auto g_xxxxzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 447);

    auto g_xxxxzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 448);

    auto g_xxxxzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 449);

    #pragma omp simd aligned(g_xxxxz_0_xxxxxxxx_0, g_xxxxz_0_xxxxxxxx_1, g_xxxxz_0_xxxxxxxy_0, g_xxxxz_0_xxxxxxxy_1, g_xxxxz_0_xxxxxxyy_0, g_xxxxz_0_xxxxxxyy_1, g_xxxxz_0_xxxxxyyy_0, g_xxxxz_0_xxxxxyyy_1, g_xxxxz_0_xxxxyyyy_0, g_xxxxz_0_xxxxyyyy_1, g_xxxxz_0_xxxyyyyy_0, g_xxxxz_0_xxxyyyyy_1, g_xxxxz_0_xxyyyyyy_0, g_xxxxz_0_xxyyyyyy_1, g_xxxxz_0_xyyyyyyy_0, g_xxxxz_0_xyyyyyyy_1, g_xxxxzz_0_xxxxxxxx_1, g_xxxxzz_0_xxxxxxxy_1, g_xxxxzz_0_xxxxxxyy_1, g_xxxxzz_0_xxxxxyyy_1, g_xxxxzz_0_xxxxyyyy_1, g_xxxxzz_0_xxxyyyyy_1, g_xxxxzz_0_xxyyyyyy_1, g_xxxxzz_0_xyyyyyyy_1, g_xxxxzzz_0_xxxxxxxx_0, g_xxxxzzz_0_xxxxxxxy_0, g_xxxxzzz_0_xxxxxxxz_0, g_xxxxzzz_0_xxxxxxyy_0, g_xxxxzzz_0_xxxxxxyz_0, g_xxxxzzz_0_xxxxxxzz_0, g_xxxxzzz_0_xxxxxyyy_0, g_xxxxzzz_0_xxxxxyyz_0, g_xxxxzzz_0_xxxxxyzz_0, g_xxxxzzz_0_xxxxxzzz_0, g_xxxxzzz_0_xxxxyyyy_0, g_xxxxzzz_0_xxxxyyyz_0, g_xxxxzzz_0_xxxxyyzz_0, g_xxxxzzz_0_xxxxyzzz_0, g_xxxxzzz_0_xxxxzzzz_0, g_xxxxzzz_0_xxxyyyyy_0, g_xxxxzzz_0_xxxyyyyz_0, g_xxxxzzz_0_xxxyyyzz_0, g_xxxxzzz_0_xxxyyzzz_0, g_xxxxzzz_0_xxxyzzzz_0, g_xxxxzzz_0_xxxzzzzz_0, g_xxxxzzz_0_xxyyyyyy_0, g_xxxxzzz_0_xxyyyyyz_0, g_xxxxzzz_0_xxyyyyzz_0, g_xxxxzzz_0_xxyyyzzz_0, g_xxxxzzz_0_xxyyzzzz_0, g_xxxxzzz_0_xxyzzzzz_0, g_xxxxzzz_0_xxzzzzzz_0, g_xxxxzzz_0_xyyyyyyy_0, g_xxxxzzz_0_xyyyyyyz_0, g_xxxxzzz_0_xyyyyyzz_0, g_xxxxzzz_0_xyyyyzzz_0, g_xxxxzzz_0_xyyyzzzz_0, g_xxxxzzz_0_xyyzzzzz_0, g_xxxxzzz_0_xyzzzzzz_0, g_xxxxzzz_0_xzzzzzzz_0, g_xxxxzzz_0_yyyyyyyy_0, g_xxxxzzz_0_yyyyyyyz_0, g_xxxxzzz_0_yyyyyyzz_0, g_xxxxzzz_0_yyyyyzzz_0, g_xxxxzzz_0_yyyyzzzz_0, g_xxxxzzz_0_yyyzzzzz_0, g_xxxxzzz_0_yyzzzzzz_0, g_xxxxzzz_0_yzzzzzzz_0, g_xxxxzzz_0_zzzzzzzz_0, g_xxxzzz_0_xxxxxxxz_1, g_xxxzzz_0_xxxxxxyz_1, g_xxxzzz_0_xxxxxxz_1, g_xxxzzz_0_xxxxxxzz_1, g_xxxzzz_0_xxxxxyyz_1, g_xxxzzz_0_xxxxxyz_1, g_xxxzzz_0_xxxxxyzz_1, g_xxxzzz_0_xxxxxzz_1, g_xxxzzz_0_xxxxxzzz_1, g_xxxzzz_0_xxxxyyyz_1, g_xxxzzz_0_xxxxyyz_1, g_xxxzzz_0_xxxxyyzz_1, g_xxxzzz_0_xxxxyzz_1, g_xxxzzz_0_xxxxyzzz_1, g_xxxzzz_0_xxxxzzz_1, g_xxxzzz_0_xxxxzzzz_1, g_xxxzzz_0_xxxyyyyz_1, g_xxxzzz_0_xxxyyyz_1, g_xxxzzz_0_xxxyyyzz_1, g_xxxzzz_0_xxxyyzz_1, g_xxxzzz_0_xxxyyzzz_1, g_xxxzzz_0_xxxyzzz_1, g_xxxzzz_0_xxxyzzzz_1, g_xxxzzz_0_xxxzzzz_1, g_xxxzzz_0_xxxzzzzz_1, g_xxxzzz_0_xxyyyyyz_1, g_xxxzzz_0_xxyyyyz_1, g_xxxzzz_0_xxyyyyzz_1, g_xxxzzz_0_xxyyyzz_1, g_xxxzzz_0_xxyyyzzz_1, g_xxxzzz_0_xxyyzzz_1, g_xxxzzz_0_xxyyzzzz_1, g_xxxzzz_0_xxyzzzz_1, g_xxxzzz_0_xxyzzzzz_1, g_xxxzzz_0_xxzzzzz_1, g_xxxzzz_0_xxzzzzzz_1, g_xxxzzz_0_xyyyyyyz_1, g_xxxzzz_0_xyyyyyz_1, g_xxxzzz_0_xyyyyyzz_1, g_xxxzzz_0_xyyyyzz_1, g_xxxzzz_0_xyyyyzzz_1, g_xxxzzz_0_xyyyzzz_1, g_xxxzzz_0_xyyyzzzz_1, g_xxxzzz_0_xyyzzzz_1, g_xxxzzz_0_xyyzzzzz_1, g_xxxzzz_0_xyzzzzz_1, g_xxxzzz_0_xyzzzzzz_1, g_xxxzzz_0_xzzzzzz_1, g_xxxzzz_0_xzzzzzzz_1, g_xxxzzz_0_yyyyyyyy_1, g_xxxzzz_0_yyyyyyyz_1, g_xxxzzz_0_yyyyyyz_1, g_xxxzzz_0_yyyyyyzz_1, g_xxxzzz_0_yyyyyzz_1, g_xxxzzz_0_yyyyyzzz_1, g_xxxzzz_0_yyyyzzz_1, g_xxxzzz_0_yyyyzzzz_1, g_xxxzzz_0_yyyzzzz_1, g_xxxzzz_0_yyyzzzzz_1, g_xxxzzz_0_yyzzzzz_1, g_xxxzzz_0_yyzzzzzz_1, g_xxxzzz_0_yzzzzzz_1, g_xxxzzz_0_yzzzzzzz_1, g_xxxzzz_0_zzzzzzz_1, g_xxxzzz_0_zzzzzzzz_1, g_xxzzz_0_xxxxxxxz_0, g_xxzzz_0_xxxxxxxz_1, g_xxzzz_0_xxxxxxyz_0, g_xxzzz_0_xxxxxxyz_1, g_xxzzz_0_xxxxxxzz_0, g_xxzzz_0_xxxxxxzz_1, g_xxzzz_0_xxxxxyyz_0, g_xxzzz_0_xxxxxyyz_1, g_xxzzz_0_xxxxxyzz_0, g_xxzzz_0_xxxxxyzz_1, g_xxzzz_0_xxxxxzzz_0, g_xxzzz_0_xxxxxzzz_1, g_xxzzz_0_xxxxyyyz_0, g_xxzzz_0_xxxxyyyz_1, g_xxzzz_0_xxxxyyzz_0, g_xxzzz_0_xxxxyyzz_1, g_xxzzz_0_xxxxyzzz_0, g_xxzzz_0_xxxxyzzz_1, g_xxzzz_0_xxxxzzzz_0, g_xxzzz_0_xxxxzzzz_1, g_xxzzz_0_xxxyyyyz_0, g_xxzzz_0_xxxyyyyz_1, g_xxzzz_0_xxxyyyzz_0, g_xxzzz_0_xxxyyyzz_1, g_xxzzz_0_xxxyyzzz_0, g_xxzzz_0_xxxyyzzz_1, g_xxzzz_0_xxxyzzzz_0, g_xxzzz_0_xxxyzzzz_1, g_xxzzz_0_xxxzzzzz_0, g_xxzzz_0_xxxzzzzz_1, g_xxzzz_0_xxyyyyyz_0, g_xxzzz_0_xxyyyyyz_1, g_xxzzz_0_xxyyyyzz_0, g_xxzzz_0_xxyyyyzz_1, g_xxzzz_0_xxyyyzzz_0, g_xxzzz_0_xxyyyzzz_1, g_xxzzz_0_xxyyzzzz_0, g_xxzzz_0_xxyyzzzz_1, g_xxzzz_0_xxyzzzzz_0, g_xxzzz_0_xxyzzzzz_1, g_xxzzz_0_xxzzzzzz_0, g_xxzzz_0_xxzzzzzz_1, g_xxzzz_0_xyyyyyyz_0, g_xxzzz_0_xyyyyyyz_1, g_xxzzz_0_xyyyyyzz_0, g_xxzzz_0_xyyyyyzz_1, g_xxzzz_0_xyyyyzzz_0, g_xxzzz_0_xyyyyzzz_1, g_xxzzz_0_xyyyzzzz_0, g_xxzzz_0_xyyyzzzz_1, g_xxzzz_0_xyyzzzzz_0, g_xxzzz_0_xyyzzzzz_1, g_xxzzz_0_xyzzzzzz_0, g_xxzzz_0_xyzzzzzz_1, g_xxzzz_0_xzzzzzzz_0, g_xxzzz_0_xzzzzzzz_1, g_xxzzz_0_yyyyyyyy_0, g_xxzzz_0_yyyyyyyy_1, g_xxzzz_0_yyyyyyyz_0, g_xxzzz_0_yyyyyyyz_1, g_xxzzz_0_yyyyyyzz_0, g_xxzzz_0_yyyyyyzz_1, g_xxzzz_0_yyyyyzzz_0, g_xxzzz_0_yyyyyzzz_1, g_xxzzz_0_yyyyzzzz_0, g_xxzzz_0_yyyyzzzz_1, g_xxzzz_0_yyyzzzzz_0, g_xxzzz_0_yyyzzzzz_1, g_xxzzz_0_yyzzzzzz_0, g_xxzzz_0_yyzzzzzz_1, g_xxzzz_0_yzzzzzzz_0, g_xxzzz_0_yzzzzzzz_1, g_xxzzz_0_zzzzzzzz_0, g_xxzzz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzzz_0_xxxxxxxx_0[i] = 2.0 * g_xxxxz_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxxzz_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxxxxy_0[i] = 2.0 * g_xxxxz_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxxxxy_1[i] * fz_be_0 + g_xxxxzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxxxxz_0[i] = 3.0 * g_xxzzz_0_xxxxxxxz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xxxzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxxxyy_0[i] = 2.0 * g_xxxxz_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxxxyy_1[i] * fz_be_0 + g_xxxxzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxxxyz_0[i] = 3.0 * g_xxzzz_0_xxxxxxyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxxzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxxxzz_0[i] = 3.0 * g_xxzzz_0_xxxxxxzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xxxzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxxyyy_0[i] = 2.0 * g_xxxxz_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxxyyy_1[i] * fz_be_0 + g_xxxxzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxxyyz_0[i] = 3.0 * g_xxzzz_0_xxxxxyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxxzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxxyzz_0[i] = 3.0 * g_xxzzz_0_xxxxxyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxxzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxxzzz_0[i] = 3.0 * g_xxzzz_0_xxxxxzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xxxzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxyyyy_0[i] = 2.0 * g_xxxxz_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxyyyy_1[i] * fz_be_0 + g_xxxxzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxyyyz_0[i] = 3.0 * g_xxzzz_0_xxxxyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxyyzz_0[i] = 3.0 * g_xxzzz_0_xxxxyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxyzzz_0[i] = 3.0 * g_xxzzz_0_xxxxyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxzzzz_0[i] = 3.0 * g_xxzzz_0_xxxxzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxyyyyy_0[i] = 2.0 * g_xxxxz_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxyyyyy_1[i] * fz_be_0 + g_xxxxzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxyyyyz_0[i] = 3.0 * g_xxzzz_0_xxxyyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxyyyzz_0[i] = 3.0 * g_xxzzz_0_xxxyyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxyyzzz_0[i] = 3.0 * g_xxzzz_0_xxxyyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxyzzzz_0[i] = 3.0 * g_xxzzz_0_xxxyzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxzzzzz_0[i] = 3.0 * g_xxzzz_0_xxxzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyyyyyy_0[i] = 2.0 * g_xxxxz_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxyyyyyy_1[i] * fz_be_0 + g_xxxxzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxyyyyyz_0[i] = 3.0 * g_xxzzz_0_xxyyyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyyyyzz_0[i] = 3.0 * g_xxzzz_0_xxyyyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyyyzzz_0[i] = 3.0 * g_xxzzz_0_xxyyyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyyzzzz_0[i] = 3.0 * g_xxzzz_0_xxyyzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyzzzzz_0[i] = 3.0 * g_xxzzz_0_xxyzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxzzzzzz_0[i] = 3.0 * g_xxzzz_0_xxzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyyyyyy_0[i] = 2.0 * g_xxxxz_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxxzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xyyyyyyz_0[i] = 3.0 * g_xxzzz_0_xyyyyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyyyyzz_0[i] = 3.0 * g_xxzzz_0_xyyyyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyyyzzz_0[i] = 3.0 * g_xxzzz_0_xyyyyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyyzzzz_0[i] = 3.0 * g_xxzzz_0_xyyyzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyzzzzz_0[i] = 3.0 * g_xxzzz_0_xyyzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyzzzzzz_0[i] = 3.0 * g_xxzzz_0_xyzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xzzzzzzz_0[i] = 3.0 * g_xxzzz_0_xzzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyyyyy_0[i] = 3.0 * g_xxzzz_0_yyyyyyyy_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xxxzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyyyyz_0[i] = 3.0 * g_xxzzz_0_yyyyyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyyyzz_0[i] = 3.0 * g_xxzzz_0_yyyyyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyyzzz_0[i] = 3.0 * g_xxzzz_0_yyyyyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyzzzz_0[i] = 3.0 * g_xxzzz_0_yyyyzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyzzzzz_0[i] = 3.0 * g_xxzzz_0_yyyzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyzzzzzz_0[i] = 3.0 * g_xxzzz_0_yyzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yzzzzzzz_0[i] = 3.0 * g_xxzzz_0_yzzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_zzzzzzzz_0[i] = 3.0 * g_xxzzz_0_zzzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 450-495 components of targeted buffer : KSL

    auto g_xxxyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 450);

    auto g_xxxyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 451);

    auto g_xxxyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 452);

    auto g_xxxyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 453);

    auto g_xxxyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 454);

    auto g_xxxyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 455);

    auto g_xxxyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 456);

    auto g_xxxyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 457);

    auto g_xxxyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 458);

    auto g_xxxyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 459);

    auto g_xxxyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 460);

    auto g_xxxyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 461);

    auto g_xxxyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 462);

    auto g_xxxyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 463);

    auto g_xxxyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 464);

    auto g_xxxyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 465);

    auto g_xxxyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 466);

    auto g_xxxyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 467);

    auto g_xxxyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 468);

    auto g_xxxyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 469);

    auto g_xxxyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 470);

    auto g_xxxyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 471);

    auto g_xxxyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 472);

    auto g_xxxyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 473);

    auto g_xxxyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 474);

    auto g_xxxyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 475);

    auto g_xxxyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 476);

    auto g_xxxyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 477);

    auto g_xxxyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 478);

    auto g_xxxyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 479);

    auto g_xxxyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 480);

    auto g_xxxyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 481);

    auto g_xxxyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 482);

    auto g_xxxyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 483);

    auto g_xxxyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 484);

    auto g_xxxyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 485);

    auto g_xxxyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 486);

    auto g_xxxyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 487);

    auto g_xxxyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 488);

    auto g_xxxyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 489);

    auto g_xxxyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 490);

    auto g_xxxyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 491);

    auto g_xxxyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 492);

    auto g_xxxyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 493);

    auto g_xxxyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 494);

    #pragma omp simd aligned(g_xxxyy_0_xxxxxxxx_0, g_xxxyy_0_xxxxxxxx_1, g_xxxyy_0_xxxxxxxz_0, g_xxxyy_0_xxxxxxxz_1, g_xxxyy_0_xxxxxxzz_0, g_xxxyy_0_xxxxxxzz_1, g_xxxyy_0_xxxxxzzz_0, g_xxxyy_0_xxxxxzzz_1, g_xxxyy_0_xxxxzzzz_0, g_xxxyy_0_xxxxzzzz_1, g_xxxyy_0_xxxzzzzz_0, g_xxxyy_0_xxxzzzzz_1, g_xxxyy_0_xxzzzzzz_0, g_xxxyy_0_xxzzzzzz_1, g_xxxyy_0_xzzzzzzz_0, g_xxxyy_0_xzzzzzzz_1, g_xxxyyy_0_xxxxxxxx_1, g_xxxyyy_0_xxxxxxxz_1, g_xxxyyy_0_xxxxxxzz_1, g_xxxyyy_0_xxxxxzzz_1, g_xxxyyy_0_xxxxzzzz_1, g_xxxyyy_0_xxxzzzzz_1, g_xxxyyy_0_xxzzzzzz_1, g_xxxyyy_0_xzzzzzzz_1, g_xxxyyyy_0_xxxxxxxx_0, g_xxxyyyy_0_xxxxxxxy_0, g_xxxyyyy_0_xxxxxxxz_0, g_xxxyyyy_0_xxxxxxyy_0, g_xxxyyyy_0_xxxxxxyz_0, g_xxxyyyy_0_xxxxxxzz_0, g_xxxyyyy_0_xxxxxyyy_0, g_xxxyyyy_0_xxxxxyyz_0, g_xxxyyyy_0_xxxxxyzz_0, g_xxxyyyy_0_xxxxxzzz_0, g_xxxyyyy_0_xxxxyyyy_0, g_xxxyyyy_0_xxxxyyyz_0, g_xxxyyyy_0_xxxxyyzz_0, g_xxxyyyy_0_xxxxyzzz_0, g_xxxyyyy_0_xxxxzzzz_0, g_xxxyyyy_0_xxxyyyyy_0, g_xxxyyyy_0_xxxyyyyz_0, g_xxxyyyy_0_xxxyyyzz_0, g_xxxyyyy_0_xxxyyzzz_0, g_xxxyyyy_0_xxxyzzzz_0, g_xxxyyyy_0_xxxzzzzz_0, g_xxxyyyy_0_xxyyyyyy_0, g_xxxyyyy_0_xxyyyyyz_0, g_xxxyyyy_0_xxyyyyzz_0, g_xxxyyyy_0_xxyyyzzz_0, g_xxxyyyy_0_xxyyzzzz_0, g_xxxyyyy_0_xxyzzzzz_0, g_xxxyyyy_0_xxzzzzzz_0, g_xxxyyyy_0_xyyyyyyy_0, g_xxxyyyy_0_xyyyyyyz_0, g_xxxyyyy_0_xyyyyyzz_0, g_xxxyyyy_0_xyyyyzzz_0, g_xxxyyyy_0_xyyyzzzz_0, g_xxxyyyy_0_xyyzzzzz_0, g_xxxyyyy_0_xyzzzzzz_0, g_xxxyyyy_0_xzzzzzzz_0, g_xxxyyyy_0_yyyyyyyy_0, g_xxxyyyy_0_yyyyyyyz_0, g_xxxyyyy_0_yyyyyyzz_0, g_xxxyyyy_0_yyyyyzzz_0, g_xxxyyyy_0_yyyyzzzz_0, g_xxxyyyy_0_yyyzzzzz_0, g_xxxyyyy_0_yyzzzzzz_0, g_xxxyyyy_0_yzzzzzzz_0, g_xxxyyyy_0_zzzzzzzz_0, g_xxyyyy_0_xxxxxxxy_1, g_xxyyyy_0_xxxxxxy_1, g_xxyyyy_0_xxxxxxyy_1, g_xxyyyy_0_xxxxxxyz_1, g_xxyyyy_0_xxxxxyy_1, g_xxyyyy_0_xxxxxyyy_1, g_xxyyyy_0_xxxxxyyz_1, g_xxyyyy_0_xxxxxyz_1, g_xxyyyy_0_xxxxxyzz_1, g_xxyyyy_0_xxxxyyy_1, g_xxyyyy_0_xxxxyyyy_1, g_xxyyyy_0_xxxxyyyz_1, g_xxyyyy_0_xxxxyyz_1, g_xxyyyy_0_xxxxyyzz_1, g_xxyyyy_0_xxxxyzz_1, g_xxyyyy_0_xxxxyzzz_1, g_xxyyyy_0_xxxyyyy_1, g_xxyyyy_0_xxxyyyyy_1, g_xxyyyy_0_xxxyyyyz_1, g_xxyyyy_0_xxxyyyz_1, g_xxyyyy_0_xxxyyyzz_1, g_xxyyyy_0_xxxyyzz_1, g_xxyyyy_0_xxxyyzzz_1, g_xxyyyy_0_xxxyzzz_1, g_xxyyyy_0_xxxyzzzz_1, g_xxyyyy_0_xxyyyyy_1, g_xxyyyy_0_xxyyyyyy_1, g_xxyyyy_0_xxyyyyyz_1, g_xxyyyy_0_xxyyyyz_1, g_xxyyyy_0_xxyyyyzz_1, g_xxyyyy_0_xxyyyzz_1, g_xxyyyy_0_xxyyyzzz_1, g_xxyyyy_0_xxyyzzz_1, g_xxyyyy_0_xxyyzzzz_1, g_xxyyyy_0_xxyzzzz_1, g_xxyyyy_0_xxyzzzzz_1, g_xxyyyy_0_xyyyyyy_1, g_xxyyyy_0_xyyyyyyy_1, g_xxyyyy_0_xyyyyyyz_1, g_xxyyyy_0_xyyyyyz_1, g_xxyyyy_0_xyyyyyzz_1, g_xxyyyy_0_xyyyyzz_1, g_xxyyyy_0_xyyyyzzz_1, g_xxyyyy_0_xyyyzzz_1, g_xxyyyy_0_xyyyzzzz_1, g_xxyyyy_0_xyyzzzz_1, g_xxyyyy_0_xyyzzzzz_1, g_xxyyyy_0_xyzzzzz_1, g_xxyyyy_0_xyzzzzzz_1, g_xxyyyy_0_yyyyyyy_1, g_xxyyyy_0_yyyyyyyy_1, g_xxyyyy_0_yyyyyyyz_1, g_xxyyyy_0_yyyyyyz_1, g_xxyyyy_0_yyyyyyzz_1, g_xxyyyy_0_yyyyyzz_1, g_xxyyyy_0_yyyyyzzz_1, g_xxyyyy_0_yyyyzzz_1, g_xxyyyy_0_yyyyzzzz_1, g_xxyyyy_0_yyyzzzz_1, g_xxyyyy_0_yyyzzzzz_1, g_xxyyyy_0_yyzzzzz_1, g_xxyyyy_0_yyzzzzzz_1, g_xxyyyy_0_yzzzzzz_1, g_xxyyyy_0_yzzzzzzz_1, g_xxyyyy_0_zzzzzzzz_1, g_xyyyy_0_xxxxxxxy_0, g_xyyyy_0_xxxxxxxy_1, g_xyyyy_0_xxxxxxyy_0, g_xyyyy_0_xxxxxxyy_1, g_xyyyy_0_xxxxxxyz_0, g_xyyyy_0_xxxxxxyz_1, g_xyyyy_0_xxxxxyyy_0, g_xyyyy_0_xxxxxyyy_1, g_xyyyy_0_xxxxxyyz_0, g_xyyyy_0_xxxxxyyz_1, g_xyyyy_0_xxxxxyzz_0, g_xyyyy_0_xxxxxyzz_1, g_xyyyy_0_xxxxyyyy_0, g_xyyyy_0_xxxxyyyy_1, g_xyyyy_0_xxxxyyyz_0, g_xyyyy_0_xxxxyyyz_1, g_xyyyy_0_xxxxyyzz_0, g_xyyyy_0_xxxxyyzz_1, g_xyyyy_0_xxxxyzzz_0, g_xyyyy_0_xxxxyzzz_1, g_xyyyy_0_xxxyyyyy_0, g_xyyyy_0_xxxyyyyy_1, g_xyyyy_0_xxxyyyyz_0, g_xyyyy_0_xxxyyyyz_1, g_xyyyy_0_xxxyyyzz_0, g_xyyyy_0_xxxyyyzz_1, g_xyyyy_0_xxxyyzzz_0, g_xyyyy_0_xxxyyzzz_1, g_xyyyy_0_xxxyzzzz_0, g_xyyyy_0_xxxyzzzz_1, g_xyyyy_0_xxyyyyyy_0, g_xyyyy_0_xxyyyyyy_1, g_xyyyy_0_xxyyyyyz_0, g_xyyyy_0_xxyyyyyz_1, g_xyyyy_0_xxyyyyzz_0, g_xyyyy_0_xxyyyyzz_1, g_xyyyy_0_xxyyyzzz_0, g_xyyyy_0_xxyyyzzz_1, g_xyyyy_0_xxyyzzzz_0, g_xyyyy_0_xxyyzzzz_1, g_xyyyy_0_xxyzzzzz_0, g_xyyyy_0_xxyzzzzz_1, g_xyyyy_0_xyyyyyyy_0, g_xyyyy_0_xyyyyyyy_1, g_xyyyy_0_xyyyyyyz_0, g_xyyyy_0_xyyyyyyz_1, g_xyyyy_0_xyyyyyzz_0, g_xyyyy_0_xyyyyyzz_1, g_xyyyy_0_xyyyyzzz_0, g_xyyyy_0_xyyyyzzz_1, g_xyyyy_0_xyyyzzzz_0, g_xyyyy_0_xyyyzzzz_1, g_xyyyy_0_xyyzzzzz_0, g_xyyyy_0_xyyzzzzz_1, g_xyyyy_0_xyzzzzzz_0, g_xyyyy_0_xyzzzzzz_1, g_xyyyy_0_yyyyyyyy_0, g_xyyyy_0_yyyyyyyy_1, g_xyyyy_0_yyyyyyyz_0, g_xyyyy_0_yyyyyyyz_1, g_xyyyy_0_yyyyyyzz_0, g_xyyyy_0_yyyyyyzz_1, g_xyyyy_0_yyyyyzzz_0, g_xyyyy_0_yyyyyzzz_1, g_xyyyy_0_yyyyzzzz_0, g_xyyyy_0_yyyyzzzz_1, g_xyyyy_0_yyyzzzzz_0, g_xyyyy_0_yyyzzzzz_1, g_xyyyy_0_yyzzzzzz_0, g_xyyyy_0_yyzzzzzz_1, g_xyyyy_0_yzzzzzzz_0, g_xyyyy_0_yzzzzzzz_1, g_xyyyy_0_zzzzzzzz_0, g_xyyyy_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyyy_0_xxxxxxxx_0[i] = 3.0 * g_xxxyy_0_xxxxxxxx_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxyyy_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxxxxxy_0[i] = 2.0 * g_xyyyy_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xxyyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxxxxz_0[i] = 3.0 * g_xxxyy_0_xxxxxxxz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxxxxz_1[i] * fz_be_0 + g_xxxyyy_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxxxxyy_0[i] = 2.0 * g_xyyyy_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xxyyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxxxyz_0[i] = 2.0 * g_xyyyy_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxyyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxxxzz_0[i] = 3.0 * g_xxxyy_0_xxxxxxzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxxxzz_1[i] * fz_be_0 + g_xxxyyy_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxxxyyy_0[i] = 2.0 * g_xyyyy_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xxyyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxxyyz_0[i] = 2.0 * g_xyyyy_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxyyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxxyzz_0[i] = 2.0 * g_xyyyy_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxyyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxxzzz_0[i] = 3.0 * g_xxxyy_0_xxxxxzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxxzzz_1[i] * fz_be_0 + g_xxxyyy_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxxyyyy_0[i] = 2.0 * g_xyyyy_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xxyyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxyyyz_0[i] = 2.0 * g_xyyyy_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxyyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxyyzz_0[i] = 2.0 * g_xyyyy_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxyyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxyzzz_0[i] = 2.0 * g_xyyyy_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxyyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxzzzz_0[i] = 3.0 * g_xxxyy_0_xxxxzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxzzzz_1[i] * fz_be_0 + g_xxxyyy_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxyyyyy_0[i] = 2.0 * g_xyyyy_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxyyyyz_0[i] = 2.0 * g_xyyyy_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxyyyzz_0[i] = 2.0 * g_xyyyy_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxyyzzz_0[i] = 2.0 * g_xyyyy_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxyzzzz_0[i] = 2.0 * g_xyyyy_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxzzzzz_0[i] = 3.0 * g_xxxyy_0_xxxzzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxyyyyyy_0[i] = 2.0 * g_xyyyy_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyyyyyz_0[i] = 2.0 * g_xyyyy_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyyyyzz_0[i] = 2.0 * g_xyyyy_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyyyzzz_0[i] = 2.0 * g_xyyyy_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyyzzzz_0[i] = 2.0 * g_xyyyy_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyzzzzz_0[i] = 2.0 * g_xyyyy_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxzzzzzz_0[i] = 3.0 * g_xxxyy_0_xxzzzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxzzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xyyyyyyy_0[i] = 2.0 * g_xyyyy_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyyyyy_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyyyyyz_0[i] = 2.0 * g_xyyyy_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyyyyz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyyyyzz_0[i] = 2.0 * g_xyyyy_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyyyzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyyyzzz_0[i] = 2.0 * g_xyyyy_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyyzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyyzzzz_0[i] = 2.0 * g_xyyyy_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyzzzzz_0[i] = 2.0 * g_xyyyy_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyzzzzzz_0[i] = 2.0 * g_xyyyy_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyzzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xzzzzzzz_0[i] = 3.0 * g_xxxyy_0_xzzzzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_yyyyyyyy_0[i] = 2.0 * g_xyyyy_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyyyyy_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyyyyyz_0[i] = 2.0 * g_xyyyy_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyyyyz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyyyyzz_0[i] = 2.0 * g_xyyyy_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyyyzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyyyzzz_0[i] = 2.0 * g_xyyyy_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyyzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyyzzzz_0[i] = 2.0 * g_xyyyy_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyzzzzz_0[i] = 2.0 * g_xyyyy_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyzzzzzz_0[i] = 2.0 * g_xyyyy_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyzzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yzzzzzzz_0[i] = 2.0 * g_xyyyy_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yzzzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_zzzzzzzz_0[i] = 2.0 * g_xyyyy_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_zzzzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 495-540 components of targeted buffer : KSL

    auto g_xxxyyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 495);

    auto g_xxxyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 496);

    auto g_xxxyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 497);

    auto g_xxxyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 498);

    auto g_xxxyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 499);

    auto g_xxxyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 500);

    auto g_xxxyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 501);

    auto g_xxxyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 502);

    auto g_xxxyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 503);

    auto g_xxxyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 504);

    auto g_xxxyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 505);

    auto g_xxxyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 506);

    auto g_xxxyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 507);

    auto g_xxxyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 508);

    auto g_xxxyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 509);

    auto g_xxxyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 510);

    auto g_xxxyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 511);

    auto g_xxxyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 512);

    auto g_xxxyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 513);

    auto g_xxxyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 514);

    auto g_xxxyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 515);

    auto g_xxxyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 516);

    auto g_xxxyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 517);

    auto g_xxxyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 518);

    auto g_xxxyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 519);

    auto g_xxxyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 520);

    auto g_xxxyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 521);

    auto g_xxxyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 522);

    auto g_xxxyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 523);

    auto g_xxxyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 524);

    auto g_xxxyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 525);

    auto g_xxxyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 526);

    auto g_xxxyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 527);

    auto g_xxxyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 528);

    auto g_xxxyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 529);

    auto g_xxxyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 530);

    auto g_xxxyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 531);

    auto g_xxxyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 532);

    auto g_xxxyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 533);

    auto g_xxxyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 534);

    auto g_xxxyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 535);

    auto g_xxxyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 536);

    auto g_xxxyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 537);

    auto g_xxxyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 538);

    auto g_xxxyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 539);

    #pragma omp simd aligned(g_xxxyyy_0_xxxxxxx_1, g_xxxyyy_0_xxxxxxxx_1, g_xxxyyy_0_xxxxxxxy_1, g_xxxyyy_0_xxxxxxxz_1, g_xxxyyy_0_xxxxxxy_1, g_xxxyyy_0_xxxxxxyy_1, g_xxxyyy_0_xxxxxxyz_1, g_xxxyyy_0_xxxxxxz_1, g_xxxyyy_0_xxxxxxzz_1, g_xxxyyy_0_xxxxxyy_1, g_xxxyyy_0_xxxxxyyy_1, g_xxxyyy_0_xxxxxyyz_1, g_xxxyyy_0_xxxxxyz_1, g_xxxyyy_0_xxxxxyzz_1, g_xxxyyy_0_xxxxxzz_1, g_xxxyyy_0_xxxxxzzz_1, g_xxxyyy_0_xxxxyyy_1, g_xxxyyy_0_xxxxyyyy_1, g_xxxyyy_0_xxxxyyyz_1, g_xxxyyy_0_xxxxyyz_1, g_xxxyyy_0_xxxxyyzz_1, g_xxxyyy_0_xxxxyzz_1, g_xxxyyy_0_xxxxyzzz_1, g_xxxyyy_0_xxxxzzz_1, g_xxxyyy_0_xxxxzzzz_1, g_xxxyyy_0_xxxyyyy_1, g_xxxyyy_0_xxxyyyyy_1, g_xxxyyy_0_xxxyyyyz_1, g_xxxyyy_0_xxxyyyz_1, g_xxxyyy_0_xxxyyyzz_1, g_xxxyyy_0_xxxyyzz_1, g_xxxyyy_0_xxxyyzzz_1, g_xxxyyy_0_xxxyzzz_1, g_xxxyyy_0_xxxyzzzz_1, g_xxxyyy_0_xxxzzzz_1, g_xxxyyy_0_xxxzzzzz_1, g_xxxyyy_0_xxyyyyy_1, g_xxxyyy_0_xxyyyyyy_1, g_xxxyyy_0_xxyyyyyz_1, g_xxxyyy_0_xxyyyyz_1, g_xxxyyy_0_xxyyyyzz_1, g_xxxyyy_0_xxyyyzz_1, g_xxxyyy_0_xxyyyzzz_1, g_xxxyyy_0_xxyyzzz_1, g_xxxyyy_0_xxyyzzzz_1, g_xxxyyy_0_xxyzzzz_1, g_xxxyyy_0_xxyzzzzz_1, g_xxxyyy_0_xxzzzzz_1, g_xxxyyy_0_xxzzzzzz_1, g_xxxyyy_0_xyyyyyy_1, g_xxxyyy_0_xyyyyyyy_1, g_xxxyyy_0_xyyyyyyz_1, g_xxxyyy_0_xyyyyyz_1, g_xxxyyy_0_xyyyyyzz_1, g_xxxyyy_0_xyyyyzz_1, g_xxxyyy_0_xyyyyzzz_1, g_xxxyyy_0_xyyyzzz_1, g_xxxyyy_0_xyyyzzzz_1, g_xxxyyy_0_xyyzzzz_1, g_xxxyyy_0_xyyzzzzz_1, g_xxxyyy_0_xyzzzzz_1, g_xxxyyy_0_xyzzzzzz_1, g_xxxyyy_0_xzzzzzz_1, g_xxxyyy_0_xzzzzzzz_1, g_xxxyyy_0_yyyyyyy_1, g_xxxyyy_0_yyyyyyyy_1, g_xxxyyy_0_yyyyyyyz_1, g_xxxyyy_0_yyyyyyz_1, g_xxxyyy_0_yyyyyyzz_1, g_xxxyyy_0_yyyyyzz_1, g_xxxyyy_0_yyyyyzzz_1, g_xxxyyy_0_yyyyzzz_1, g_xxxyyy_0_yyyyzzzz_1, g_xxxyyy_0_yyyzzzz_1, g_xxxyyy_0_yyyzzzzz_1, g_xxxyyy_0_yyzzzzz_1, g_xxxyyy_0_yyzzzzzz_1, g_xxxyyy_0_yzzzzzz_1, g_xxxyyy_0_yzzzzzzz_1, g_xxxyyy_0_zzzzzzz_1, g_xxxyyy_0_zzzzzzzz_1, g_xxxyyyz_0_xxxxxxxx_0, g_xxxyyyz_0_xxxxxxxy_0, g_xxxyyyz_0_xxxxxxxz_0, g_xxxyyyz_0_xxxxxxyy_0, g_xxxyyyz_0_xxxxxxyz_0, g_xxxyyyz_0_xxxxxxzz_0, g_xxxyyyz_0_xxxxxyyy_0, g_xxxyyyz_0_xxxxxyyz_0, g_xxxyyyz_0_xxxxxyzz_0, g_xxxyyyz_0_xxxxxzzz_0, g_xxxyyyz_0_xxxxyyyy_0, g_xxxyyyz_0_xxxxyyyz_0, g_xxxyyyz_0_xxxxyyzz_0, g_xxxyyyz_0_xxxxyzzz_0, g_xxxyyyz_0_xxxxzzzz_0, g_xxxyyyz_0_xxxyyyyy_0, g_xxxyyyz_0_xxxyyyyz_0, g_xxxyyyz_0_xxxyyyzz_0, g_xxxyyyz_0_xxxyyzzz_0, g_xxxyyyz_0_xxxyzzzz_0, g_xxxyyyz_0_xxxzzzzz_0, g_xxxyyyz_0_xxyyyyyy_0, g_xxxyyyz_0_xxyyyyyz_0, g_xxxyyyz_0_xxyyyyzz_0, g_xxxyyyz_0_xxyyyzzz_0, g_xxxyyyz_0_xxyyzzzz_0, g_xxxyyyz_0_xxyzzzzz_0, g_xxxyyyz_0_xxzzzzzz_0, g_xxxyyyz_0_xyyyyyyy_0, g_xxxyyyz_0_xyyyyyyz_0, g_xxxyyyz_0_xyyyyyzz_0, g_xxxyyyz_0_xyyyyzzz_0, g_xxxyyyz_0_xyyyzzzz_0, g_xxxyyyz_0_xyyzzzzz_0, g_xxxyyyz_0_xyzzzzzz_0, g_xxxyyyz_0_xzzzzzzz_0, g_xxxyyyz_0_yyyyyyyy_0, g_xxxyyyz_0_yyyyyyyz_0, g_xxxyyyz_0_yyyyyyzz_0, g_xxxyyyz_0_yyyyyzzz_0, g_xxxyyyz_0_yyyyzzzz_0, g_xxxyyyz_0_yyyzzzzz_0, g_xxxyyyz_0_yyzzzzzz_0, g_xxxyyyz_0_yzzzzzzz_0, g_xxxyyyz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyyz_0_xxxxxxxx_0[i] = g_xxxyyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxxxy_0[i] = g_xxxyyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxxxz_0[i] = g_xxxyyy_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxxxz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxxyy_0[i] = g_xxxyyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxxyz_0[i] = g_xxxyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxxyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxxzz_0[i] = 2.0 * g_xxxyyy_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxxzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxyyy_0[i] = g_xxxyyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxyyz_0[i] = g_xxxyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxyzz_0[i] = 2.0 * g_xxxyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxzzz_0[i] = 3.0 * g_xxxyyy_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxyyyy_0[i] = g_xxxyyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxyyyz_0[i] = g_xxxyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxyyzz_0[i] = 2.0 * g_xxxyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxyzzz_0[i] = 3.0 * g_xxxyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxzzzz_0[i] = 4.0 * g_xxxyyy_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyyyyy_0[i] = g_xxxyyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyyyyz_0[i] = g_xxxyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyyyzz_0[i] = 2.0 * g_xxxyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyyzzz_0[i] = 3.0 * g_xxxyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyzzzz_0[i] = 4.0 * g_xxxyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxzzzzz_0[i] = 5.0 * g_xxxyyy_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyyyyy_0[i] = g_xxxyyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyyyyz_0[i] = g_xxxyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyyyzz_0[i] = 2.0 * g_xxxyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyyzzz_0[i] = 3.0 * g_xxxyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyzzzz_0[i] = 4.0 * g_xxxyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyzzzzz_0[i] = 5.0 * g_xxxyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxzzzzzz_0[i] = 6.0 * g_xxxyyy_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxzzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyyyyy_0[i] = g_xxxyyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyyyyz_0[i] = g_xxxyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyyyzz_0[i] = 2.0 * g_xxxyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyyzzz_0[i] = 3.0 * g_xxxyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyzzzz_0[i] = 4.0 * g_xxxyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyzzzzz_0[i] = 5.0 * g_xxxyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyzzzzzz_0[i] = 6.0 * g_xxxyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyzzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xzzzzzzz_0[i] = 7.0 * g_xxxyyy_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xzzzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyyyyy_0[i] = g_xxxyyy_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyyyyz_0[i] = g_xxxyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_yyyyyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyyyzz_0[i] = 2.0 * g_xxxyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_yyyyyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyyzzz_0[i] = 3.0 * g_xxxyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_yyyyyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyzzzz_0[i] = 4.0 * g_xxxyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_yyyyzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyzzzzz_0[i] = 5.0 * g_xxxyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_yyyzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyzzzzzz_0[i] = 6.0 * g_xxxyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_yyzzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yzzzzzzz_0[i] = 7.0 * g_xxxyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_yzzzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_zzzzzzzz_0[i] = 8.0 * g_xxxyyy_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 540-585 components of targeted buffer : KSL

    auto g_xxxyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 540);

    auto g_xxxyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 541);

    auto g_xxxyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 542);

    auto g_xxxyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 543);

    auto g_xxxyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 544);

    auto g_xxxyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 545);

    auto g_xxxyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 546);

    auto g_xxxyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 547);

    auto g_xxxyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 548);

    auto g_xxxyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 549);

    auto g_xxxyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 550);

    auto g_xxxyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 551);

    auto g_xxxyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 552);

    auto g_xxxyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 553);

    auto g_xxxyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 554);

    auto g_xxxyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 555);

    auto g_xxxyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 556);

    auto g_xxxyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 557);

    auto g_xxxyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 558);

    auto g_xxxyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 559);

    auto g_xxxyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 560);

    auto g_xxxyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 561);

    auto g_xxxyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 562);

    auto g_xxxyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 563);

    auto g_xxxyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 564);

    auto g_xxxyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 565);

    auto g_xxxyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 566);

    auto g_xxxyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 567);

    auto g_xxxyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 568);

    auto g_xxxyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 569);

    auto g_xxxyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 570);

    auto g_xxxyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 571);

    auto g_xxxyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 572);

    auto g_xxxyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 573);

    auto g_xxxyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 574);

    auto g_xxxyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 575);

    auto g_xxxyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 576);

    auto g_xxxyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 577);

    auto g_xxxyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 578);

    auto g_xxxyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 579);

    auto g_xxxyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 580);

    auto g_xxxyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 581);

    auto g_xxxyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 582);

    auto g_xxxyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 583);

    auto g_xxxyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 584);

    #pragma omp simd aligned(g_xxxyy_0_xxxxxxxy_0, g_xxxyy_0_xxxxxxxy_1, g_xxxyy_0_xxxxxxyy_0, g_xxxyy_0_xxxxxxyy_1, g_xxxyy_0_xxxxxyyy_0, g_xxxyy_0_xxxxxyyy_1, g_xxxyy_0_xxxxyyyy_0, g_xxxyy_0_xxxxyyyy_1, g_xxxyy_0_xxxyyyyy_0, g_xxxyy_0_xxxyyyyy_1, g_xxxyy_0_xxyyyyyy_0, g_xxxyy_0_xxyyyyyy_1, g_xxxyy_0_xyyyyyyy_0, g_xxxyy_0_xyyyyyyy_1, g_xxxyyz_0_xxxxxxxy_1, g_xxxyyz_0_xxxxxxyy_1, g_xxxyyz_0_xxxxxyyy_1, g_xxxyyz_0_xxxxyyyy_1, g_xxxyyz_0_xxxyyyyy_1, g_xxxyyz_0_xxyyyyyy_1, g_xxxyyz_0_xyyyyyyy_1, g_xxxyyzz_0_xxxxxxxx_0, g_xxxyyzz_0_xxxxxxxy_0, g_xxxyyzz_0_xxxxxxxz_0, g_xxxyyzz_0_xxxxxxyy_0, g_xxxyyzz_0_xxxxxxyz_0, g_xxxyyzz_0_xxxxxxzz_0, g_xxxyyzz_0_xxxxxyyy_0, g_xxxyyzz_0_xxxxxyyz_0, g_xxxyyzz_0_xxxxxyzz_0, g_xxxyyzz_0_xxxxxzzz_0, g_xxxyyzz_0_xxxxyyyy_0, g_xxxyyzz_0_xxxxyyyz_0, g_xxxyyzz_0_xxxxyyzz_0, g_xxxyyzz_0_xxxxyzzz_0, g_xxxyyzz_0_xxxxzzzz_0, g_xxxyyzz_0_xxxyyyyy_0, g_xxxyyzz_0_xxxyyyyz_0, g_xxxyyzz_0_xxxyyyzz_0, g_xxxyyzz_0_xxxyyzzz_0, g_xxxyyzz_0_xxxyzzzz_0, g_xxxyyzz_0_xxxzzzzz_0, g_xxxyyzz_0_xxyyyyyy_0, g_xxxyyzz_0_xxyyyyyz_0, g_xxxyyzz_0_xxyyyyzz_0, g_xxxyyzz_0_xxyyyzzz_0, g_xxxyyzz_0_xxyyzzzz_0, g_xxxyyzz_0_xxyzzzzz_0, g_xxxyyzz_0_xxzzzzzz_0, g_xxxyyzz_0_xyyyyyyy_0, g_xxxyyzz_0_xyyyyyyz_0, g_xxxyyzz_0_xyyyyyzz_0, g_xxxyyzz_0_xyyyyzzz_0, g_xxxyyzz_0_xyyyzzzz_0, g_xxxyyzz_0_xyyzzzzz_0, g_xxxyyzz_0_xyzzzzzz_0, g_xxxyyzz_0_xzzzzzzz_0, g_xxxyyzz_0_yyyyyyyy_0, g_xxxyyzz_0_yyyyyyyz_0, g_xxxyyzz_0_yyyyyyzz_0, g_xxxyyzz_0_yyyyyzzz_0, g_xxxyyzz_0_yyyyzzzz_0, g_xxxyyzz_0_yyyzzzzz_0, g_xxxyyzz_0_yyzzzzzz_0, g_xxxyyzz_0_yzzzzzzz_0, g_xxxyyzz_0_zzzzzzzz_0, g_xxxyzz_0_xxxxxxxx_1, g_xxxyzz_0_xxxxxxxz_1, g_xxxyzz_0_xxxxxxzz_1, g_xxxyzz_0_xxxxxzzz_1, g_xxxyzz_0_xxxxzzzz_1, g_xxxyzz_0_xxxzzzzz_1, g_xxxyzz_0_xxzzzzzz_1, g_xxxyzz_0_xzzzzzzz_1, g_xxxzz_0_xxxxxxxx_0, g_xxxzz_0_xxxxxxxx_1, g_xxxzz_0_xxxxxxxz_0, g_xxxzz_0_xxxxxxxz_1, g_xxxzz_0_xxxxxxzz_0, g_xxxzz_0_xxxxxxzz_1, g_xxxzz_0_xxxxxzzz_0, g_xxxzz_0_xxxxxzzz_1, g_xxxzz_0_xxxxzzzz_0, g_xxxzz_0_xxxxzzzz_1, g_xxxzz_0_xxxzzzzz_0, g_xxxzz_0_xxxzzzzz_1, g_xxxzz_0_xxzzzzzz_0, g_xxxzz_0_xxzzzzzz_1, g_xxxzz_0_xzzzzzzz_0, g_xxxzz_0_xzzzzzzz_1, g_xxyyzz_0_xxxxxxyz_1, g_xxyyzz_0_xxxxxyyz_1, g_xxyyzz_0_xxxxxyz_1, g_xxyyzz_0_xxxxxyzz_1, g_xxyyzz_0_xxxxyyyz_1, g_xxyyzz_0_xxxxyyz_1, g_xxyyzz_0_xxxxyyzz_1, g_xxyyzz_0_xxxxyzz_1, g_xxyyzz_0_xxxxyzzz_1, g_xxyyzz_0_xxxyyyyz_1, g_xxyyzz_0_xxxyyyz_1, g_xxyyzz_0_xxxyyyzz_1, g_xxyyzz_0_xxxyyzz_1, g_xxyyzz_0_xxxyyzzz_1, g_xxyyzz_0_xxxyzzz_1, g_xxyyzz_0_xxxyzzzz_1, g_xxyyzz_0_xxyyyyyz_1, g_xxyyzz_0_xxyyyyz_1, g_xxyyzz_0_xxyyyyzz_1, g_xxyyzz_0_xxyyyzz_1, g_xxyyzz_0_xxyyyzzz_1, g_xxyyzz_0_xxyyzzz_1, g_xxyyzz_0_xxyyzzzz_1, g_xxyyzz_0_xxyzzzz_1, g_xxyyzz_0_xxyzzzzz_1, g_xxyyzz_0_xyyyyyyz_1, g_xxyyzz_0_xyyyyyz_1, g_xxyyzz_0_xyyyyyzz_1, g_xxyyzz_0_xyyyyzz_1, g_xxyyzz_0_xyyyyzzz_1, g_xxyyzz_0_xyyyzzz_1, g_xxyyzz_0_xyyyzzzz_1, g_xxyyzz_0_xyyzzzz_1, g_xxyyzz_0_xyyzzzzz_1, g_xxyyzz_0_xyzzzzz_1, g_xxyyzz_0_xyzzzzzz_1, g_xxyyzz_0_yyyyyyyy_1, g_xxyyzz_0_yyyyyyyz_1, g_xxyyzz_0_yyyyyyz_1, g_xxyyzz_0_yyyyyyzz_1, g_xxyyzz_0_yyyyyzz_1, g_xxyyzz_0_yyyyyzzz_1, g_xxyyzz_0_yyyyzzz_1, g_xxyyzz_0_yyyyzzzz_1, g_xxyyzz_0_yyyzzzz_1, g_xxyyzz_0_yyyzzzzz_1, g_xxyyzz_0_yyzzzzz_1, g_xxyyzz_0_yyzzzzzz_1, g_xxyyzz_0_yzzzzzz_1, g_xxyyzz_0_yzzzzzzz_1, g_xxyyzz_0_zzzzzzzz_1, g_xyyzz_0_xxxxxxyz_0, g_xyyzz_0_xxxxxxyz_1, g_xyyzz_0_xxxxxyyz_0, g_xyyzz_0_xxxxxyyz_1, g_xyyzz_0_xxxxxyzz_0, g_xyyzz_0_xxxxxyzz_1, g_xyyzz_0_xxxxyyyz_0, g_xyyzz_0_xxxxyyyz_1, g_xyyzz_0_xxxxyyzz_0, g_xyyzz_0_xxxxyyzz_1, g_xyyzz_0_xxxxyzzz_0, g_xyyzz_0_xxxxyzzz_1, g_xyyzz_0_xxxyyyyz_0, g_xyyzz_0_xxxyyyyz_1, g_xyyzz_0_xxxyyyzz_0, g_xyyzz_0_xxxyyyzz_1, g_xyyzz_0_xxxyyzzz_0, g_xyyzz_0_xxxyyzzz_1, g_xyyzz_0_xxxyzzzz_0, g_xyyzz_0_xxxyzzzz_1, g_xyyzz_0_xxyyyyyz_0, g_xyyzz_0_xxyyyyyz_1, g_xyyzz_0_xxyyyyzz_0, g_xyyzz_0_xxyyyyzz_1, g_xyyzz_0_xxyyyzzz_0, g_xyyzz_0_xxyyyzzz_1, g_xyyzz_0_xxyyzzzz_0, g_xyyzz_0_xxyyzzzz_1, g_xyyzz_0_xxyzzzzz_0, g_xyyzz_0_xxyzzzzz_1, g_xyyzz_0_xyyyyyyz_0, g_xyyzz_0_xyyyyyyz_1, g_xyyzz_0_xyyyyyzz_0, g_xyyzz_0_xyyyyyzz_1, g_xyyzz_0_xyyyyzzz_0, g_xyyzz_0_xyyyyzzz_1, g_xyyzz_0_xyyyzzzz_0, g_xyyzz_0_xyyyzzzz_1, g_xyyzz_0_xyyzzzzz_0, g_xyyzz_0_xyyzzzzz_1, g_xyyzz_0_xyzzzzzz_0, g_xyyzz_0_xyzzzzzz_1, g_xyyzz_0_yyyyyyyy_0, g_xyyzz_0_yyyyyyyy_1, g_xyyzz_0_yyyyyyyz_0, g_xyyzz_0_yyyyyyyz_1, g_xyyzz_0_yyyyyyzz_0, g_xyyzz_0_yyyyyyzz_1, g_xyyzz_0_yyyyyzzz_0, g_xyyzz_0_yyyyyzzz_1, g_xyyzz_0_yyyyzzzz_0, g_xyyzz_0_yyyyzzzz_1, g_xyyzz_0_yyyzzzzz_0, g_xyyzz_0_yyyzzzzz_1, g_xyyzz_0_yyzzzzzz_0, g_xyyzz_0_yyzzzzzz_1, g_xyyzz_0_yzzzzzzz_0, g_xyyzz_0_yzzzzzzz_1, g_xyyzz_0_zzzzzzzz_0, g_xyyzz_0_zzzzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyzz_0_xxxxxxxx_0[i] = g_xxxzz_0_xxxxxxxx_0[i] * fbe_0 - g_xxxzz_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxyzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxxxxxy_0[i] = g_xxxyy_0_xxxxxxxy_0[i] * fbe_0 - g_xxxyy_0_xxxxxxxy_1[i] * fz_be_0 + g_xxxyyz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxxxxxz_0[i] = g_xxxzz_0_xxxxxxxz_0[i] * fbe_0 - g_xxxzz_0_xxxxxxxz_1[i] * fz_be_0 + g_xxxyzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxxxxyy_0[i] = g_xxxyy_0_xxxxxxyy_0[i] * fbe_0 - g_xxxyy_0_xxxxxxyy_1[i] * fz_be_0 + g_xxxyyz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxxxxyz_0[i] = 2.0 * g_xyyzz_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxyyzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxxxxzz_0[i] = g_xxxzz_0_xxxxxxzz_0[i] * fbe_0 - g_xxxzz_0_xxxxxxzz_1[i] * fz_be_0 + g_xxxyzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxxxyyy_0[i] = g_xxxyy_0_xxxxxyyy_0[i] * fbe_0 - g_xxxyy_0_xxxxxyyy_1[i] * fz_be_0 + g_xxxyyz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxxxyyz_0[i] = 2.0 * g_xyyzz_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxyyzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxxxyzz_0[i] = 2.0 * g_xyyzz_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxyyzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxxxzzz_0[i] = g_xxxzz_0_xxxxxzzz_0[i] * fbe_0 - g_xxxzz_0_xxxxxzzz_1[i] * fz_be_0 + g_xxxyzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxxyyyy_0[i] = g_xxxyy_0_xxxxyyyy_0[i] * fbe_0 - g_xxxyy_0_xxxxyyyy_1[i] * fz_be_0 + g_xxxyyz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxxyyyz_0[i] = 2.0 * g_xyyzz_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxyyzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxxyyzz_0[i] = 2.0 * g_xyyzz_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxyyzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxxyzzz_0[i] = 2.0 * g_xyyzz_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxyyzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxxzzzz_0[i] = g_xxxzz_0_xxxxzzzz_0[i] * fbe_0 - g_xxxzz_0_xxxxzzzz_1[i] * fz_be_0 + g_xxxyzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxyyyyy_0[i] = g_xxxyy_0_xxxyyyyy_0[i] * fbe_0 - g_xxxyy_0_xxxyyyyy_1[i] * fz_be_0 + g_xxxyyz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxyyyyz_0[i] = 2.0 * g_xyyzz_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxyyyzz_0[i] = 2.0 * g_xyyzz_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxyyzzz_0[i] = 2.0 * g_xyyzz_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxyzzzz_0[i] = 2.0 * g_xyyzz_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxzzzzz_0[i] = g_xxxzz_0_xxxzzzzz_0[i] * fbe_0 - g_xxxzz_0_xxxzzzzz_1[i] * fz_be_0 + g_xxxyzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxyyyyyy_0[i] = g_xxxyy_0_xxyyyyyy_0[i] * fbe_0 - g_xxxyy_0_xxyyyyyy_1[i] * fz_be_0 + g_xxxyyz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxyyyyyz_0[i] = 2.0 * g_xyyzz_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxyyyyzz_0[i] = 2.0 * g_xyyzz_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxyyyzzz_0[i] = 2.0 * g_xyyzz_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxyyzzzz_0[i] = 2.0 * g_xyyzz_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxyzzzzz_0[i] = 2.0 * g_xyyzz_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxzzzzzz_0[i] = g_xxxzz_0_xxzzzzzz_0[i] * fbe_0 - g_xxxzz_0_xxzzzzzz_1[i] * fz_be_0 + g_xxxyzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xyyyyyyy_0[i] = g_xxxyy_0_xyyyyyyy_0[i] * fbe_0 - g_xxxyy_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxyyz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xyyyyyyz_0[i] = 2.0 * g_xyyzz_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyyyyyzz_0[i] = 2.0 * g_xyyzz_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyyyyzzz_0[i] = 2.0 * g_xyyzz_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyyyzzzz_0[i] = 2.0 * g_xyyzz_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyyzzzzz_0[i] = 2.0 * g_xyyzz_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyzzzzzz_0[i] = 2.0 * g_xyyzz_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xzzzzzzz_0[i] = g_xxxzz_0_xzzzzzzz_0[i] * fbe_0 - g_xxxzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxyzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_yyyyyyyy_0[i] = 2.0 * g_xyyzz_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xxyyzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyyyyyz_0[i] = 2.0 * g_xyyzz_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyyyyzz_0[i] = 2.0 * g_xyyzz_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyyyzzz_0[i] = 2.0 * g_xyyzz_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyyzzzz_0[i] = 2.0 * g_xyyzz_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyzzzzz_0[i] = 2.0 * g_xyyzz_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyzzzzzz_0[i] = 2.0 * g_xyyzz_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yzzzzzzz_0[i] = 2.0 * g_xyyzz_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_zzzzzzzz_0[i] = 2.0 * g_xyyzz_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 585-630 components of targeted buffer : KSL

    auto g_xxxyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 585);

    auto g_xxxyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 586);

    auto g_xxxyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 587);

    auto g_xxxyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 588);

    auto g_xxxyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 589);

    auto g_xxxyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 590);

    auto g_xxxyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 591);

    auto g_xxxyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 592);

    auto g_xxxyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 593);

    auto g_xxxyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 594);

    auto g_xxxyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 595);

    auto g_xxxyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 596);

    auto g_xxxyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 597);

    auto g_xxxyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 598);

    auto g_xxxyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 599);

    auto g_xxxyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 600);

    auto g_xxxyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 601);

    auto g_xxxyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 602);

    auto g_xxxyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 603);

    auto g_xxxyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 604);

    auto g_xxxyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 605);

    auto g_xxxyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 606);

    auto g_xxxyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 607);

    auto g_xxxyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 608);

    auto g_xxxyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 609);

    auto g_xxxyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 610);

    auto g_xxxyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 611);

    auto g_xxxyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 612);

    auto g_xxxyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 613);

    auto g_xxxyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 614);

    auto g_xxxyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 615);

    auto g_xxxyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 616);

    auto g_xxxyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 617);

    auto g_xxxyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 618);

    auto g_xxxyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 619);

    auto g_xxxyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 620);

    auto g_xxxyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 621);

    auto g_xxxyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 622);

    auto g_xxxyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 623);

    auto g_xxxyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 624);

    auto g_xxxyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 625);

    auto g_xxxyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 626);

    auto g_xxxyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 627);

    auto g_xxxyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 628);

    auto g_xxxyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 629);

    #pragma omp simd aligned(g_xxxyzzz_0_xxxxxxxx_0, g_xxxyzzz_0_xxxxxxxy_0, g_xxxyzzz_0_xxxxxxxz_0, g_xxxyzzz_0_xxxxxxyy_0, g_xxxyzzz_0_xxxxxxyz_0, g_xxxyzzz_0_xxxxxxzz_0, g_xxxyzzz_0_xxxxxyyy_0, g_xxxyzzz_0_xxxxxyyz_0, g_xxxyzzz_0_xxxxxyzz_0, g_xxxyzzz_0_xxxxxzzz_0, g_xxxyzzz_0_xxxxyyyy_0, g_xxxyzzz_0_xxxxyyyz_0, g_xxxyzzz_0_xxxxyyzz_0, g_xxxyzzz_0_xxxxyzzz_0, g_xxxyzzz_0_xxxxzzzz_0, g_xxxyzzz_0_xxxyyyyy_0, g_xxxyzzz_0_xxxyyyyz_0, g_xxxyzzz_0_xxxyyyzz_0, g_xxxyzzz_0_xxxyyzzz_0, g_xxxyzzz_0_xxxyzzzz_0, g_xxxyzzz_0_xxxzzzzz_0, g_xxxyzzz_0_xxyyyyyy_0, g_xxxyzzz_0_xxyyyyyz_0, g_xxxyzzz_0_xxyyyyzz_0, g_xxxyzzz_0_xxyyyzzz_0, g_xxxyzzz_0_xxyyzzzz_0, g_xxxyzzz_0_xxyzzzzz_0, g_xxxyzzz_0_xxzzzzzz_0, g_xxxyzzz_0_xyyyyyyy_0, g_xxxyzzz_0_xyyyyyyz_0, g_xxxyzzz_0_xyyyyyzz_0, g_xxxyzzz_0_xyyyyzzz_0, g_xxxyzzz_0_xyyyzzzz_0, g_xxxyzzz_0_xyyzzzzz_0, g_xxxyzzz_0_xyzzzzzz_0, g_xxxyzzz_0_xzzzzzzz_0, g_xxxyzzz_0_yyyyyyyy_0, g_xxxyzzz_0_yyyyyyyz_0, g_xxxyzzz_0_yyyyyyzz_0, g_xxxyzzz_0_yyyyyzzz_0, g_xxxyzzz_0_yyyyzzzz_0, g_xxxyzzz_0_yyyzzzzz_0, g_xxxyzzz_0_yyzzzzzz_0, g_xxxyzzz_0_yzzzzzzz_0, g_xxxyzzz_0_zzzzzzzz_0, g_xxxzzz_0_xxxxxxx_1, g_xxxzzz_0_xxxxxxxx_1, g_xxxzzz_0_xxxxxxxy_1, g_xxxzzz_0_xxxxxxxz_1, g_xxxzzz_0_xxxxxxy_1, g_xxxzzz_0_xxxxxxyy_1, g_xxxzzz_0_xxxxxxyz_1, g_xxxzzz_0_xxxxxxz_1, g_xxxzzz_0_xxxxxxzz_1, g_xxxzzz_0_xxxxxyy_1, g_xxxzzz_0_xxxxxyyy_1, g_xxxzzz_0_xxxxxyyz_1, g_xxxzzz_0_xxxxxyz_1, g_xxxzzz_0_xxxxxyzz_1, g_xxxzzz_0_xxxxxzz_1, g_xxxzzz_0_xxxxxzzz_1, g_xxxzzz_0_xxxxyyy_1, g_xxxzzz_0_xxxxyyyy_1, g_xxxzzz_0_xxxxyyyz_1, g_xxxzzz_0_xxxxyyz_1, g_xxxzzz_0_xxxxyyzz_1, g_xxxzzz_0_xxxxyzz_1, g_xxxzzz_0_xxxxyzzz_1, g_xxxzzz_0_xxxxzzz_1, g_xxxzzz_0_xxxxzzzz_1, g_xxxzzz_0_xxxyyyy_1, g_xxxzzz_0_xxxyyyyy_1, g_xxxzzz_0_xxxyyyyz_1, g_xxxzzz_0_xxxyyyz_1, g_xxxzzz_0_xxxyyyzz_1, g_xxxzzz_0_xxxyyzz_1, g_xxxzzz_0_xxxyyzzz_1, g_xxxzzz_0_xxxyzzz_1, g_xxxzzz_0_xxxyzzzz_1, g_xxxzzz_0_xxxzzzz_1, g_xxxzzz_0_xxxzzzzz_1, g_xxxzzz_0_xxyyyyy_1, g_xxxzzz_0_xxyyyyyy_1, g_xxxzzz_0_xxyyyyyz_1, g_xxxzzz_0_xxyyyyz_1, g_xxxzzz_0_xxyyyyzz_1, g_xxxzzz_0_xxyyyzz_1, g_xxxzzz_0_xxyyyzzz_1, g_xxxzzz_0_xxyyzzz_1, g_xxxzzz_0_xxyyzzzz_1, g_xxxzzz_0_xxyzzzz_1, g_xxxzzz_0_xxyzzzzz_1, g_xxxzzz_0_xxzzzzz_1, g_xxxzzz_0_xxzzzzzz_1, g_xxxzzz_0_xyyyyyy_1, g_xxxzzz_0_xyyyyyyy_1, g_xxxzzz_0_xyyyyyyz_1, g_xxxzzz_0_xyyyyyz_1, g_xxxzzz_0_xyyyyyzz_1, g_xxxzzz_0_xyyyyzz_1, g_xxxzzz_0_xyyyyzzz_1, g_xxxzzz_0_xyyyzzz_1, g_xxxzzz_0_xyyyzzzz_1, g_xxxzzz_0_xyyzzzz_1, g_xxxzzz_0_xyyzzzzz_1, g_xxxzzz_0_xyzzzzz_1, g_xxxzzz_0_xyzzzzzz_1, g_xxxzzz_0_xzzzzzz_1, g_xxxzzz_0_xzzzzzzz_1, g_xxxzzz_0_yyyyyyy_1, g_xxxzzz_0_yyyyyyyy_1, g_xxxzzz_0_yyyyyyyz_1, g_xxxzzz_0_yyyyyyz_1, g_xxxzzz_0_yyyyyyzz_1, g_xxxzzz_0_yyyyyzz_1, g_xxxzzz_0_yyyyyzzz_1, g_xxxzzz_0_yyyyzzz_1, g_xxxzzz_0_yyyyzzzz_1, g_xxxzzz_0_yyyzzzz_1, g_xxxzzz_0_yyyzzzzz_1, g_xxxzzz_0_yyzzzzz_1, g_xxxzzz_0_yyzzzzzz_1, g_xxxzzz_0_yzzzzzz_1, g_xxxzzz_0_yzzzzzzz_1, g_xxxzzz_0_zzzzzzz_1, g_xxxzzz_0_zzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzzz_0_xxxxxxxx_0[i] = g_xxxzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxxxy_0[i] = g_xxxzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxxxy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxxxz_0[i] = g_xxxzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxxyy_0[i] = 2.0 * g_xxxzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxxyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxxyz_0[i] = g_xxxzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxxzz_0[i] = g_xxxzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxyyy_0[i] = 3.0 * g_xxxzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxyyz_0[i] = 2.0 * g_xxxzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxyzz_0[i] = g_xxxzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxzzz_0[i] = g_xxxzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxyyyy_0[i] = 4.0 * g_xxxzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxyyyz_0[i] = 3.0 * g_xxxzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxyyzz_0[i] = 2.0 * g_xxxzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxyzzz_0[i] = g_xxxzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxzzzz_0[i] = g_xxxzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyyyyy_0[i] = 5.0 * g_xxxzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyyyyz_0[i] = 4.0 * g_xxxzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyyyzz_0[i] = 3.0 * g_xxxzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyyzzz_0[i] = 2.0 * g_xxxzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyzzzz_0[i] = g_xxxzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxzzzzz_0[i] = g_xxxzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyyyyy_0[i] = 6.0 * g_xxxzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyyyyz_0[i] = 5.0 * g_xxxzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyyyzz_0[i] = 4.0 * g_xxxzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyyzzz_0[i] = 3.0 * g_xxxzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyzzzz_0[i] = 2.0 * g_xxxzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyzzzzz_0[i] = g_xxxzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxzzzzzz_0[i] = g_xxxzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyyyyy_0[i] = 7.0 * g_xxxzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyyyyz_0[i] = 6.0 * g_xxxzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyyyzz_0[i] = 5.0 * g_xxxzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyyzzz_0[i] = 4.0 * g_xxxzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyzzzz_0[i] = 3.0 * g_xxxzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyzzzzz_0[i] = 2.0 * g_xxxzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyzzzzzz_0[i] = g_xxxzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xzzzzzzz_0[i] = g_xxxzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyyyyy_0[i] = 8.0 * g_xxxzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyyyyz_0[i] = 7.0 * g_xxxzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyyyzz_0[i] = 6.0 * g_xxxzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyyzzz_0[i] = 5.0 * g_xxxzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyzzzz_0[i] = 4.0 * g_xxxzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyzzzzz_0[i] = 3.0 * g_xxxzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyzzzzzz_0[i] = 2.0 * g_xxxzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yzzzzzzz_0[i] = g_xxxzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_zzzzzzzz_0[i] = g_xxxzzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 630-675 components of targeted buffer : KSL

    auto g_xxxzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 630);

    auto g_xxxzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 631);

    auto g_xxxzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 632);

    auto g_xxxzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 633);

    auto g_xxxzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 634);

    auto g_xxxzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 635);

    auto g_xxxzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 636);

    auto g_xxxzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 637);

    auto g_xxxzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 638);

    auto g_xxxzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 639);

    auto g_xxxzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 640);

    auto g_xxxzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 641);

    auto g_xxxzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 642);

    auto g_xxxzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 643);

    auto g_xxxzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 644);

    auto g_xxxzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 645);

    auto g_xxxzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 646);

    auto g_xxxzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 647);

    auto g_xxxzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 648);

    auto g_xxxzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 649);

    auto g_xxxzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 650);

    auto g_xxxzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 651);

    auto g_xxxzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 652);

    auto g_xxxzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 653);

    auto g_xxxzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 654);

    auto g_xxxzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 655);

    auto g_xxxzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 656);

    auto g_xxxzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 657);

    auto g_xxxzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 658);

    auto g_xxxzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 659);

    auto g_xxxzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 660);

    auto g_xxxzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 661);

    auto g_xxxzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 662);

    auto g_xxxzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 663);

    auto g_xxxzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 664);

    auto g_xxxzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 665);

    auto g_xxxzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 666);

    auto g_xxxzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 667);

    auto g_xxxzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 668);

    auto g_xxxzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 669);

    auto g_xxxzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 670);

    auto g_xxxzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 671);

    auto g_xxxzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 672);

    auto g_xxxzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 673);

    auto g_xxxzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 674);

    #pragma omp simd aligned(g_xxxzz_0_xxxxxxxx_0, g_xxxzz_0_xxxxxxxx_1, g_xxxzz_0_xxxxxxxy_0, g_xxxzz_0_xxxxxxxy_1, g_xxxzz_0_xxxxxxyy_0, g_xxxzz_0_xxxxxxyy_1, g_xxxzz_0_xxxxxyyy_0, g_xxxzz_0_xxxxxyyy_1, g_xxxzz_0_xxxxyyyy_0, g_xxxzz_0_xxxxyyyy_1, g_xxxzz_0_xxxyyyyy_0, g_xxxzz_0_xxxyyyyy_1, g_xxxzz_0_xxyyyyyy_0, g_xxxzz_0_xxyyyyyy_1, g_xxxzz_0_xyyyyyyy_0, g_xxxzz_0_xyyyyyyy_1, g_xxxzzz_0_xxxxxxxx_1, g_xxxzzz_0_xxxxxxxy_1, g_xxxzzz_0_xxxxxxyy_1, g_xxxzzz_0_xxxxxyyy_1, g_xxxzzz_0_xxxxyyyy_1, g_xxxzzz_0_xxxyyyyy_1, g_xxxzzz_0_xxyyyyyy_1, g_xxxzzz_0_xyyyyyyy_1, g_xxxzzzz_0_xxxxxxxx_0, g_xxxzzzz_0_xxxxxxxy_0, g_xxxzzzz_0_xxxxxxxz_0, g_xxxzzzz_0_xxxxxxyy_0, g_xxxzzzz_0_xxxxxxyz_0, g_xxxzzzz_0_xxxxxxzz_0, g_xxxzzzz_0_xxxxxyyy_0, g_xxxzzzz_0_xxxxxyyz_0, g_xxxzzzz_0_xxxxxyzz_0, g_xxxzzzz_0_xxxxxzzz_0, g_xxxzzzz_0_xxxxyyyy_0, g_xxxzzzz_0_xxxxyyyz_0, g_xxxzzzz_0_xxxxyyzz_0, g_xxxzzzz_0_xxxxyzzz_0, g_xxxzzzz_0_xxxxzzzz_0, g_xxxzzzz_0_xxxyyyyy_0, g_xxxzzzz_0_xxxyyyyz_0, g_xxxzzzz_0_xxxyyyzz_0, g_xxxzzzz_0_xxxyyzzz_0, g_xxxzzzz_0_xxxyzzzz_0, g_xxxzzzz_0_xxxzzzzz_0, g_xxxzzzz_0_xxyyyyyy_0, g_xxxzzzz_0_xxyyyyyz_0, g_xxxzzzz_0_xxyyyyzz_0, g_xxxzzzz_0_xxyyyzzz_0, g_xxxzzzz_0_xxyyzzzz_0, g_xxxzzzz_0_xxyzzzzz_0, g_xxxzzzz_0_xxzzzzzz_0, g_xxxzzzz_0_xyyyyyyy_0, g_xxxzzzz_0_xyyyyyyz_0, g_xxxzzzz_0_xyyyyyzz_0, g_xxxzzzz_0_xyyyyzzz_0, g_xxxzzzz_0_xyyyzzzz_0, g_xxxzzzz_0_xyyzzzzz_0, g_xxxzzzz_0_xyzzzzzz_0, g_xxxzzzz_0_xzzzzzzz_0, g_xxxzzzz_0_yyyyyyyy_0, g_xxxzzzz_0_yyyyyyyz_0, g_xxxzzzz_0_yyyyyyzz_0, g_xxxzzzz_0_yyyyyzzz_0, g_xxxzzzz_0_yyyyzzzz_0, g_xxxzzzz_0_yyyzzzzz_0, g_xxxzzzz_0_yyzzzzzz_0, g_xxxzzzz_0_yzzzzzzz_0, g_xxxzzzz_0_zzzzzzzz_0, g_xxzzzz_0_xxxxxxxz_1, g_xxzzzz_0_xxxxxxyz_1, g_xxzzzz_0_xxxxxxz_1, g_xxzzzz_0_xxxxxxzz_1, g_xxzzzz_0_xxxxxyyz_1, g_xxzzzz_0_xxxxxyz_1, g_xxzzzz_0_xxxxxyzz_1, g_xxzzzz_0_xxxxxzz_1, g_xxzzzz_0_xxxxxzzz_1, g_xxzzzz_0_xxxxyyyz_1, g_xxzzzz_0_xxxxyyz_1, g_xxzzzz_0_xxxxyyzz_1, g_xxzzzz_0_xxxxyzz_1, g_xxzzzz_0_xxxxyzzz_1, g_xxzzzz_0_xxxxzzz_1, g_xxzzzz_0_xxxxzzzz_1, g_xxzzzz_0_xxxyyyyz_1, g_xxzzzz_0_xxxyyyz_1, g_xxzzzz_0_xxxyyyzz_1, g_xxzzzz_0_xxxyyzz_1, g_xxzzzz_0_xxxyyzzz_1, g_xxzzzz_0_xxxyzzz_1, g_xxzzzz_0_xxxyzzzz_1, g_xxzzzz_0_xxxzzzz_1, g_xxzzzz_0_xxxzzzzz_1, g_xxzzzz_0_xxyyyyyz_1, g_xxzzzz_0_xxyyyyz_1, g_xxzzzz_0_xxyyyyzz_1, g_xxzzzz_0_xxyyyzz_1, g_xxzzzz_0_xxyyyzzz_1, g_xxzzzz_0_xxyyzzz_1, g_xxzzzz_0_xxyyzzzz_1, g_xxzzzz_0_xxyzzzz_1, g_xxzzzz_0_xxyzzzzz_1, g_xxzzzz_0_xxzzzzz_1, g_xxzzzz_0_xxzzzzzz_1, g_xxzzzz_0_xyyyyyyz_1, g_xxzzzz_0_xyyyyyz_1, g_xxzzzz_0_xyyyyyzz_1, g_xxzzzz_0_xyyyyzz_1, g_xxzzzz_0_xyyyyzzz_1, g_xxzzzz_0_xyyyzzz_1, g_xxzzzz_0_xyyyzzzz_1, g_xxzzzz_0_xyyzzzz_1, g_xxzzzz_0_xyyzzzzz_1, g_xxzzzz_0_xyzzzzz_1, g_xxzzzz_0_xyzzzzzz_1, g_xxzzzz_0_xzzzzzz_1, g_xxzzzz_0_xzzzzzzz_1, g_xxzzzz_0_yyyyyyyy_1, g_xxzzzz_0_yyyyyyyz_1, g_xxzzzz_0_yyyyyyz_1, g_xxzzzz_0_yyyyyyzz_1, g_xxzzzz_0_yyyyyzz_1, g_xxzzzz_0_yyyyyzzz_1, g_xxzzzz_0_yyyyzzz_1, g_xxzzzz_0_yyyyzzzz_1, g_xxzzzz_0_yyyzzzz_1, g_xxzzzz_0_yyyzzzzz_1, g_xxzzzz_0_yyzzzzz_1, g_xxzzzz_0_yyzzzzzz_1, g_xxzzzz_0_yzzzzzz_1, g_xxzzzz_0_yzzzzzzz_1, g_xxzzzz_0_zzzzzzz_1, g_xxzzzz_0_zzzzzzzz_1, g_xzzzz_0_xxxxxxxz_0, g_xzzzz_0_xxxxxxxz_1, g_xzzzz_0_xxxxxxyz_0, g_xzzzz_0_xxxxxxyz_1, g_xzzzz_0_xxxxxxzz_0, g_xzzzz_0_xxxxxxzz_1, g_xzzzz_0_xxxxxyyz_0, g_xzzzz_0_xxxxxyyz_1, g_xzzzz_0_xxxxxyzz_0, g_xzzzz_0_xxxxxyzz_1, g_xzzzz_0_xxxxxzzz_0, g_xzzzz_0_xxxxxzzz_1, g_xzzzz_0_xxxxyyyz_0, g_xzzzz_0_xxxxyyyz_1, g_xzzzz_0_xxxxyyzz_0, g_xzzzz_0_xxxxyyzz_1, g_xzzzz_0_xxxxyzzz_0, g_xzzzz_0_xxxxyzzz_1, g_xzzzz_0_xxxxzzzz_0, g_xzzzz_0_xxxxzzzz_1, g_xzzzz_0_xxxyyyyz_0, g_xzzzz_0_xxxyyyyz_1, g_xzzzz_0_xxxyyyzz_0, g_xzzzz_0_xxxyyyzz_1, g_xzzzz_0_xxxyyzzz_0, g_xzzzz_0_xxxyyzzz_1, g_xzzzz_0_xxxyzzzz_0, g_xzzzz_0_xxxyzzzz_1, g_xzzzz_0_xxxzzzzz_0, g_xzzzz_0_xxxzzzzz_1, g_xzzzz_0_xxyyyyyz_0, g_xzzzz_0_xxyyyyyz_1, g_xzzzz_0_xxyyyyzz_0, g_xzzzz_0_xxyyyyzz_1, g_xzzzz_0_xxyyyzzz_0, g_xzzzz_0_xxyyyzzz_1, g_xzzzz_0_xxyyzzzz_0, g_xzzzz_0_xxyyzzzz_1, g_xzzzz_0_xxyzzzzz_0, g_xzzzz_0_xxyzzzzz_1, g_xzzzz_0_xxzzzzzz_0, g_xzzzz_0_xxzzzzzz_1, g_xzzzz_0_xyyyyyyz_0, g_xzzzz_0_xyyyyyyz_1, g_xzzzz_0_xyyyyyzz_0, g_xzzzz_0_xyyyyyzz_1, g_xzzzz_0_xyyyyzzz_0, g_xzzzz_0_xyyyyzzz_1, g_xzzzz_0_xyyyzzzz_0, g_xzzzz_0_xyyyzzzz_1, g_xzzzz_0_xyyzzzzz_0, g_xzzzz_0_xyyzzzzz_1, g_xzzzz_0_xyzzzzzz_0, g_xzzzz_0_xyzzzzzz_1, g_xzzzz_0_xzzzzzzz_0, g_xzzzz_0_xzzzzzzz_1, g_xzzzz_0_yyyyyyyy_0, g_xzzzz_0_yyyyyyyy_1, g_xzzzz_0_yyyyyyyz_0, g_xzzzz_0_yyyyyyyz_1, g_xzzzz_0_yyyyyyzz_0, g_xzzzz_0_yyyyyyzz_1, g_xzzzz_0_yyyyyzzz_0, g_xzzzz_0_yyyyyzzz_1, g_xzzzz_0_yyyyzzzz_0, g_xzzzz_0_yyyyzzzz_1, g_xzzzz_0_yyyzzzzz_0, g_xzzzz_0_yyyzzzzz_1, g_xzzzz_0_yyzzzzzz_0, g_xzzzz_0_yyzzzzzz_1, g_xzzzz_0_yzzzzzzz_0, g_xzzzz_0_yzzzzzzz_1, g_xzzzz_0_zzzzzzzz_0, g_xzzzz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzzz_0_xxxxxxxx_0[i] = 3.0 * g_xxxzz_0_xxxxxxxx_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxzzz_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxxxxy_0[i] = 3.0 * g_xxxzz_0_xxxxxxxy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxxxxy_1[i] * fz_be_0 + g_xxxzzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxxxxz_0[i] = 2.0 * g_xzzzz_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xxzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxxxyy_0[i] = 3.0 * g_xxxzz_0_xxxxxxyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxxxyy_1[i] * fz_be_0 + g_xxxzzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxxxyz_0[i] = 2.0 * g_xzzzz_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxxxzz_0[i] = 2.0 * g_xzzzz_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xxzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxxyyy_0[i] = 3.0 * g_xxxzz_0_xxxxxyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxxyyy_1[i] * fz_be_0 + g_xxxzzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxxyyz_0[i] = 2.0 * g_xzzzz_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxxyzz_0[i] = 2.0 * g_xzzzz_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxxzzz_0[i] = 2.0 * g_xzzzz_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xxzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxyyyy_0[i] = 3.0 * g_xxxzz_0_xxxxyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxyyyy_1[i] * fz_be_0 + g_xxxzzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxyyyz_0[i] = 2.0 * g_xzzzz_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxyyzz_0[i] = 2.0 * g_xzzzz_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxyzzz_0[i] = 2.0 * g_xzzzz_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxzzzz_0[i] = 2.0 * g_xzzzz_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxyyyyy_0[i] = 3.0 * g_xxxzz_0_xxxyyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxyyyyy_1[i] * fz_be_0 + g_xxxzzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxyyyyz_0[i] = 2.0 * g_xzzzz_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxyyyzz_0[i] = 2.0 * g_xzzzz_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxyyzzz_0[i] = 2.0 * g_xzzzz_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxyzzzz_0[i] = 2.0 * g_xzzzz_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxzzzzz_0[i] = 2.0 * g_xzzzz_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyyyyyy_0[i] = 3.0 * g_xxxzz_0_xxyyyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxyyyyyy_1[i] * fz_be_0 + g_xxxzzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxyyyyyz_0[i] = 2.0 * g_xzzzz_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyyyyzz_0[i] = 2.0 * g_xzzzz_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyyyzzz_0[i] = 2.0 * g_xzzzz_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyyzzzz_0[i] = 2.0 * g_xzzzz_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyzzzzz_0[i] = 2.0 * g_xzzzz_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxzzzzzz_0[i] = 2.0 * g_xzzzz_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyyyyyy_0[i] = 3.0 * g_xxxzz_0_xyyyyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxzzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xyyyyyyz_0[i] = 2.0 * g_xzzzz_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyyyyzz_0[i] = 2.0 * g_xzzzz_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyyyzzz_0[i] = 2.0 * g_xzzzz_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyyzzzz_0[i] = 2.0 * g_xzzzz_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyzzzzz_0[i] = 2.0 * g_xzzzz_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyzzzzzz_0[i] = 2.0 * g_xzzzz_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xzzzzzzz_0[i] = 2.0 * g_xzzzz_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyyyyy_0[i] = 2.0 * g_xzzzz_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xxzzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyyyyz_0[i] = 2.0 * g_xzzzz_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyyyzz_0[i] = 2.0 * g_xzzzz_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyyzzz_0[i] = 2.0 * g_xzzzz_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyzzzz_0[i] = 2.0 * g_xzzzz_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyzzzzz_0[i] = 2.0 * g_xzzzz_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyzzzzzz_0[i] = 2.0 * g_xzzzz_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yzzzzzzz_0[i] = 2.0 * g_xzzzz_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_zzzzzzzz_0[i] = 2.0 * g_xzzzz_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 675-720 components of targeted buffer : KSL

    auto g_xxyyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 675);

    auto g_xxyyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 676);

    auto g_xxyyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 677);

    auto g_xxyyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 678);

    auto g_xxyyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 679);

    auto g_xxyyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 680);

    auto g_xxyyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 681);

    auto g_xxyyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 682);

    auto g_xxyyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 683);

    auto g_xxyyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 684);

    auto g_xxyyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 685);

    auto g_xxyyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 686);

    auto g_xxyyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 687);

    auto g_xxyyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 688);

    auto g_xxyyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 689);

    auto g_xxyyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 690);

    auto g_xxyyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 691);

    auto g_xxyyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 692);

    auto g_xxyyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 693);

    auto g_xxyyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 694);

    auto g_xxyyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 695);

    auto g_xxyyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 696);

    auto g_xxyyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 697);

    auto g_xxyyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 698);

    auto g_xxyyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 699);

    auto g_xxyyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 700);

    auto g_xxyyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 701);

    auto g_xxyyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 702);

    auto g_xxyyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 703);

    auto g_xxyyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 704);

    auto g_xxyyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 705);

    auto g_xxyyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 706);

    auto g_xxyyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 707);

    auto g_xxyyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 708);

    auto g_xxyyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 709);

    auto g_xxyyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 710);

    auto g_xxyyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 711);

    auto g_xxyyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 712);

    auto g_xxyyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 713);

    auto g_xxyyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 714);

    auto g_xxyyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 715);

    auto g_xxyyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 716);

    auto g_xxyyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 717);

    auto g_xxyyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 718);

    auto g_xxyyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 719);

    #pragma omp simd aligned(g_xxyyy_0_xxxxxxxx_0, g_xxyyy_0_xxxxxxxx_1, g_xxyyy_0_xxxxxxxz_0, g_xxyyy_0_xxxxxxxz_1, g_xxyyy_0_xxxxxxzz_0, g_xxyyy_0_xxxxxxzz_1, g_xxyyy_0_xxxxxzzz_0, g_xxyyy_0_xxxxxzzz_1, g_xxyyy_0_xxxxzzzz_0, g_xxyyy_0_xxxxzzzz_1, g_xxyyy_0_xxxzzzzz_0, g_xxyyy_0_xxxzzzzz_1, g_xxyyy_0_xxzzzzzz_0, g_xxyyy_0_xxzzzzzz_1, g_xxyyy_0_xzzzzzzz_0, g_xxyyy_0_xzzzzzzz_1, g_xxyyyy_0_xxxxxxxx_1, g_xxyyyy_0_xxxxxxxz_1, g_xxyyyy_0_xxxxxxzz_1, g_xxyyyy_0_xxxxxzzz_1, g_xxyyyy_0_xxxxzzzz_1, g_xxyyyy_0_xxxzzzzz_1, g_xxyyyy_0_xxzzzzzz_1, g_xxyyyy_0_xzzzzzzz_1, g_xxyyyyy_0_xxxxxxxx_0, g_xxyyyyy_0_xxxxxxxy_0, g_xxyyyyy_0_xxxxxxxz_0, g_xxyyyyy_0_xxxxxxyy_0, g_xxyyyyy_0_xxxxxxyz_0, g_xxyyyyy_0_xxxxxxzz_0, g_xxyyyyy_0_xxxxxyyy_0, g_xxyyyyy_0_xxxxxyyz_0, g_xxyyyyy_0_xxxxxyzz_0, g_xxyyyyy_0_xxxxxzzz_0, g_xxyyyyy_0_xxxxyyyy_0, g_xxyyyyy_0_xxxxyyyz_0, g_xxyyyyy_0_xxxxyyzz_0, g_xxyyyyy_0_xxxxyzzz_0, g_xxyyyyy_0_xxxxzzzz_0, g_xxyyyyy_0_xxxyyyyy_0, g_xxyyyyy_0_xxxyyyyz_0, g_xxyyyyy_0_xxxyyyzz_0, g_xxyyyyy_0_xxxyyzzz_0, g_xxyyyyy_0_xxxyzzzz_0, g_xxyyyyy_0_xxxzzzzz_0, g_xxyyyyy_0_xxyyyyyy_0, g_xxyyyyy_0_xxyyyyyz_0, g_xxyyyyy_0_xxyyyyzz_0, g_xxyyyyy_0_xxyyyzzz_0, g_xxyyyyy_0_xxyyzzzz_0, g_xxyyyyy_0_xxyzzzzz_0, g_xxyyyyy_0_xxzzzzzz_0, g_xxyyyyy_0_xyyyyyyy_0, g_xxyyyyy_0_xyyyyyyz_0, g_xxyyyyy_0_xyyyyyzz_0, g_xxyyyyy_0_xyyyyzzz_0, g_xxyyyyy_0_xyyyzzzz_0, g_xxyyyyy_0_xyyzzzzz_0, g_xxyyyyy_0_xyzzzzzz_0, g_xxyyyyy_0_xzzzzzzz_0, g_xxyyyyy_0_yyyyyyyy_0, g_xxyyyyy_0_yyyyyyyz_0, g_xxyyyyy_0_yyyyyyzz_0, g_xxyyyyy_0_yyyyyzzz_0, g_xxyyyyy_0_yyyyzzzz_0, g_xxyyyyy_0_yyyzzzzz_0, g_xxyyyyy_0_yyzzzzzz_0, g_xxyyyyy_0_yzzzzzzz_0, g_xxyyyyy_0_zzzzzzzz_0, g_xyyyyy_0_xxxxxxxy_1, g_xyyyyy_0_xxxxxxy_1, g_xyyyyy_0_xxxxxxyy_1, g_xyyyyy_0_xxxxxxyz_1, g_xyyyyy_0_xxxxxyy_1, g_xyyyyy_0_xxxxxyyy_1, g_xyyyyy_0_xxxxxyyz_1, g_xyyyyy_0_xxxxxyz_1, g_xyyyyy_0_xxxxxyzz_1, g_xyyyyy_0_xxxxyyy_1, g_xyyyyy_0_xxxxyyyy_1, g_xyyyyy_0_xxxxyyyz_1, g_xyyyyy_0_xxxxyyz_1, g_xyyyyy_0_xxxxyyzz_1, g_xyyyyy_0_xxxxyzz_1, g_xyyyyy_0_xxxxyzzz_1, g_xyyyyy_0_xxxyyyy_1, g_xyyyyy_0_xxxyyyyy_1, g_xyyyyy_0_xxxyyyyz_1, g_xyyyyy_0_xxxyyyz_1, g_xyyyyy_0_xxxyyyzz_1, g_xyyyyy_0_xxxyyzz_1, g_xyyyyy_0_xxxyyzzz_1, g_xyyyyy_0_xxxyzzz_1, g_xyyyyy_0_xxxyzzzz_1, g_xyyyyy_0_xxyyyyy_1, g_xyyyyy_0_xxyyyyyy_1, g_xyyyyy_0_xxyyyyyz_1, g_xyyyyy_0_xxyyyyz_1, g_xyyyyy_0_xxyyyyzz_1, g_xyyyyy_0_xxyyyzz_1, g_xyyyyy_0_xxyyyzzz_1, g_xyyyyy_0_xxyyzzz_1, g_xyyyyy_0_xxyyzzzz_1, g_xyyyyy_0_xxyzzzz_1, g_xyyyyy_0_xxyzzzzz_1, g_xyyyyy_0_xyyyyyy_1, g_xyyyyy_0_xyyyyyyy_1, g_xyyyyy_0_xyyyyyyz_1, g_xyyyyy_0_xyyyyyz_1, g_xyyyyy_0_xyyyyyzz_1, g_xyyyyy_0_xyyyyzz_1, g_xyyyyy_0_xyyyyzzz_1, g_xyyyyy_0_xyyyzzz_1, g_xyyyyy_0_xyyyzzzz_1, g_xyyyyy_0_xyyzzzz_1, g_xyyyyy_0_xyyzzzzz_1, g_xyyyyy_0_xyzzzzz_1, g_xyyyyy_0_xyzzzzzz_1, g_xyyyyy_0_yyyyyyy_1, g_xyyyyy_0_yyyyyyyy_1, g_xyyyyy_0_yyyyyyyz_1, g_xyyyyy_0_yyyyyyz_1, g_xyyyyy_0_yyyyyyzz_1, g_xyyyyy_0_yyyyyzz_1, g_xyyyyy_0_yyyyyzzz_1, g_xyyyyy_0_yyyyzzz_1, g_xyyyyy_0_yyyyzzzz_1, g_xyyyyy_0_yyyzzzz_1, g_xyyyyy_0_yyyzzzzz_1, g_xyyyyy_0_yyzzzzz_1, g_xyyyyy_0_yyzzzzzz_1, g_xyyyyy_0_yzzzzzz_1, g_xyyyyy_0_yzzzzzzz_1, g_xyyyyy_0_zzzzzzzz_1, g_yyyyy_0_xxxxxxxy_0, g_yyyyy_0_xxxxxxxy_1, g_yyyyy_0_xxxxxxyy_0, g_yyyyy_0_xxxxxxyy_1, g_yyyyy_0_xxxxxxyz_0, g_yyyyy_0_xxxxxxyz_1, g_yyyyy_0_xxxxxyyy_0, g_yyyyy_0_xxxxxyyy_1, g_yyyyy_0_xxxxxyyz_0, g_yyyyy_0_xxxxxyyz_1, g_yyyyy_0_xxxxxyzz_0, g_yyyyy_0_xxxxxyzz_1, g_yyyyy_0_xxxxyyyy_0, g_yyyyy_0_xxxxyyyy_1, g_yyyyy_0_xxxxyyyz_0, g_yyyyy_0_xxxxyyyz_1, g_yyyyy_0_xxxxyyzz_0, g_yyyyy_0_xxxxyyzz_1, g_yyyyy_0_xxxxyzzz_0, g_yyyyy_0_xxxxyzzz_1, g_yyyyy_0_xxxyyyyy_0, g_yyyyy_0_xxxyyyyy_1, g_yyyyy_0_xxxyyyyz_0, g_yyyyy_0_xxxyyyyz_1, g_yyyyy_0_xxxyyyzz_0, g_yyyyy_0_xxxyyyzz_1, g_yyyyy_0_xxxyyzzz_0, g_yyyyy_0_xxxyyzzz_1, g_yyyyy_0_xxxyzzzz_0, g_yyyyy_0_xxxyzzzz_1, g_yyyyy_0_xxyyyyyy_0, g_yyyyy_0_xxyyyyyy_1, g_yyyyy_0_xxyyyyyz_0, g_yyyyy_0_xxyyyyyz_1, g_yyyyy_0_xxyyyyzz_0, g_yyyyy_0_xxyyyyzz_1, g_yyyyy_0_xxyyyzzz_0, g_yyyyy_0_xxyyyzzz_1, g_yyyyy_0_xxyyzzzz_0, g_yyyyy_0_xxyyzzzz_1, g_yyyyy_0_xxyzzzzz_0, g_yyyyy_0_xxyzzzzz_1, g_yyyyy_0_xyyyyyyy_0, g_yyyyy_0_xyyyyyyy_1, g_yyyyy_0_xyyyyyyz_0, g_yyyyy_0_xyyyyyyz_1, g_yyyyy_0_xyyyyyzz_0, g_yyyyy_0_xyyyyyzz_1, g_yyyyy_0_xyyyyzzz_0, g_yyyyy_0_xyyyyzzz_1, g_yyyyy_0_xyyyzzzz_0, g_yyyyy_0_xyyyzzzz_1, g_yyyyy_0_xyyzzzzz_0, g_yyyyy_0_xyyzzzzz_1, g_yyyyy_0_xyzzzzzz_0, g_yyyyy_0_xyzzzzzz_1, g_yyyyy_0_yyyyyyyy_0, g_yyyyy_0_yyyyyyyy_1, g_yyyyy_0_yyyyyyyz_0, g_yyyyy_0_yyyyyyyz_1, g_yyyyy_0_yyyyyyzz_0, g_yyyyy_0_yyyyyyzz_1, g_yyyyy_0_yyyyyzzz_0, g_yyyyy_0_yyyyyzzz_1, g_yyyyy_0_yyyyzzzz_0, g_yyyyy_0_yyyyzzzz_1, g_yyyyy_0_yyyzzzzz_0, g_yyyyy_0_yyyzzzzz_1, g_yyyyy_0_yyzzzzzz_0, g_yyyyy_0_yyzzzzzz_1, g_yyyyy_0_yzzzzzzz_0, g_yyyyy_0_yzzzzzzz_1, g_yyyyy_0_zzzzzzzz_0, g_yyyyy_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyyy_0_xxxxxxxx_0[i] = 4.0 * g_xxyyy_0_xxxxxxxx_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxxxxx_1[i] * fz_be_0 + g_xxyyyy_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxxxxxy_0[i] = g_yyyyy_0_xxxxxxxy_0[i] * fbe_0 - g_yyyyy_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xyyyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxxxxz_0[i] = 4.0 * g_xxyyy_0_xxxxxxxz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxxxxz_1[i] * fz_be_0 + g_xxyyyy_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxxxxyy_0[i] = g_yyyyy_0_xxxxxxyy_0[i] * fbe_0 - g_yyyyy_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xyyyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxxxyz_0[i] = g_yyyyy_0_xxxxxxyz_0[i] * fbe_0 - g_yyyyy_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xyyyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxxxzz_0[i] = 4.0 * g_xxyyy_0_xxxxxxzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxxxzz_1[i] * fz_be_0 + g_xxyyyy_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxxxyyy_0[i] = g_yyyyy_0_xxxxxyyy_0[i] * fbe_0 - g_yyyyy_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xyyyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxxyyz_0[i] = g_yyyyy_0_xxxxxyyz_0[i] * fbe_0 - g_yyyyy_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xyyyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxxyzz_0[i] = g_yyyyy_0_xxxxxyzz_0[i] * fbe_0 - g_yyyyy_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xyyyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxxzzz_0[i] = 4.0 * g_xxyyy_0_xxxxxzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxxzzz_1[i] * fz_be_0 + g_xxyyyy_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxxyyyy_0[i] = g_yyyyy_0_xxxxyyyy_0[i] * fbe_0 - g_yyyyy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xyyyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxyyyz_0[i] = g_yyyyy_0_xxxxyyyz_0[i] * fbe_0 - g_yyyyy_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xyyyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxyyzz_0[i] = g_yyyyy_0_xxxxyyzz_0[i] * fbe_0 - g_yyyyy_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xyyyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxyzzz_0[i] = g_yyyyy_0_xxxxyzzz_0[i] * fbe_0 - g_yyyyy_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xyyyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxzzzz_0[i] = 4.0 * g_xxyyy_0_xxxxzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxzzzz_1[i] * fz_be_0 + g_xxyyyy_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxyyyyy_0[i] = g_yyyyy_0_xxxyyyyy_0[i] * fbe_0 - g_yyyyy_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxyyyyz_0[i] = g_yyyyy_0_xxxyyyyz_0[i] * fbe_0 - g_yyyyy_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxyyyzz_0[i] = g_yyyyy_0_xxxyyyzz_0[i] * fbe_0 - g_yyyyy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxyyzzz_0[i] = g_yyyyy_0_xxxyyzzz_0[i] * fbe_0 - g_yyyyy_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxyzzzz_0[i] = g_yyyyy_0_xxxyzzzz_0[i] * fbe_0 - g_yyyyy_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxzzzzz_0[i] = 4.0 * g_xxyyy_0_xxxzzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxyyyyyy_0[i] = g_yyyyy_0_xxyyyyyy_0[i] * fbe_0 - g_yyyyy_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyyyyyz_0[i] = g_yyyyy_0_xxyyyyyz_0[i] * fbe_0 - g_yyyyy_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyyyyzz_0[i] = g_yyyyy_0_xxyyyyzz_0[i] * fbe_0 - g_yyyyy_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyyyzzz_0[i] = g_yyyyy_0_xxyyyzzz_0[i] * fbe_0 - g_yyyyy_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyyzzzz_0[i] = g_yyyyy_0_xxyyzzzz_0[i] * fbe_0 - g_yyyyy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyzzzzz_0[i] = g_yyyyy_0_xxyzzzzz_0[i] * fbe_0 - g_yyyyy_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxzzzzzz_0[i] = 4.0 * g_xxyyy_0_xxzzzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxzzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xyyyyyyy_0[i] = g_yyyyy_0_xyyyyyyy_0[i] * fbe_0 - g_yyyyy_0_xyyyyyyy_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyyyyyz_0[i] = g_yyyyy_0_xyyyyyyz_0[i] * fbe_0 - g_yyyyy_0_xyyyyyyz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyyyyzz_0[i] = g_yyyyy_0_xyyyyyzz_0[i] * fbe_0 - g_yyyyy_0_xyyyyyzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyyyzzz_0[i] = g_yyyyy_0_xyyyyzzz_0[i] * fbe_0 - g_yyyyy_0_xyyyyzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyyzzzz_0[i] = g_yyyyy_0_xyyyzzzz_0[i] * fbe_0 - g_yyyyy_0_xyyyzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyzzzzz_0[i] = g_yyyyy_0_xyyzzzzz_0[i] * fbe_0 - g_yyyyy_0_xyyzzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyzzzzzz_0[i] = g_yyyyy_0_xyzzzzzz_0[i] * fbe_0 - g_yyyyy_0_xyzzzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xzzzzzzz_0[i] = 4.0 * g_xxyyy_0_xzzzzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xzzzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_yyyyyyyy_0[i] = g_yyyyy_0_yyyyyyyy_0[i] * fbe_0 - g_yyyyy_0_yyyyyyyy_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyyyyyz_0[i] = g_yyyyy_0_yyyyyyyz_0[i] * fbe_0 - g_yyyyy_0_yyyyyyyz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyyyyzz_0[i] = g_yyyyy_0_yyyyyyzz_0[i] * fbe_0 - g_yyyyy_0_yyyyyyzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyyyzzz_0[i] = g_yyyyy_0_yyyyyzzz_0[i] * fbe_0 - g_yyyyy_0_yyyyyzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyyzzzz_0[i] = g_yyyyy_0_yyyyzzzz_0[i] * fbe_0 - g_yyyyy_0_yyyyzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyzzzzz_0[i] = g_yyyyy_0_yyyzzzzz_0[i] * fbe_0 - g_yyyyy_0_yyyzzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyzzzzzz_0[i] = g_yyyyy_0_yyzzzzzz_0[i] * fbe_0 - g_yyyyy_0_yyzzzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yzzzzzzz_0[i] = g_yyyyy_0_yzzzzzzz_0[i] * fbe_0 - g_yyyyy_0_yzzzzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_zzzzzzzz_0[i] = g_yyyyy_0_zzzzzzzz_0[i] * fbe_0 - g_yyyyy_0_zzzzzzzz_1[i] * fz_be_0 + g_xyyyyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 720-765 components of targeted buffer : KSL

    auto g_xxyyyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 720);

    auto g_xxyyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 721);

    auto g_xxyyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 722);

    auto g_xxyyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 723);

    auto g_xxyyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 724);

    auto g_xxyyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 725);

    auto g_xxyyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 726);

    auto g_xxyyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 727);

    auto g_xxyyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 728);

    auto g_xxyyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 729);

    auto g_xxyyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 730);

    auto g_xxyyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 731);

    auto g_xxyyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 732);

    auto g_xxyyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 733);

    auto g_xxyyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 734);

    auto g_xxyyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 735);

    auto g_xxyyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 736);

    auto g_xxyyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 737);

    auto g_xxyyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 738);

    auto g_xxyyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 739);

    auto g_xxyyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 740);

    auto g_xxyyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 741);

    auto g_xxyyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 742);

    auto g_xxyyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 743);

    auto g_xxyyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 744);

    auto g_xxyyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 745);

    auto g_xxyyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 746);

    auto g_xxyyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 747);

    auto g_xxyyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 748);

    auto g_xxyyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 749);

    auto g_xxyyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 750);

    auto g_xxyyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 751);

    auto g_xxyyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 752);

    auto g_xxyyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 753);

    auto g_xxyyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 754);

    auto g_xxyyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 755);

    auto g_xxyyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 756);

    auto g_xxyyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 757);

    auto g_xxyyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 758);

    auto g_xxyyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 759);

    auto g_xxyyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 760);

    auto g_xxyyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 761);

    auto g_xxyyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 762);

    auto g_xxyyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 763);

    auto g_xxyyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 764);

    #pragma omp simd aligned(g_xxyyyy_0_xxxxxxx_1, g_xxyyyy_0_xxxxxxxx_1, g_xxyyyy_0_xxxxxxxy_1, g_xxyyyy_0_xxxxxxxz_1, g_xxyyyy_0_xxxxxxy_1, g_xxyyyy_0_xxxxxxyy_1, g_xxyyyy_0_xxxxxxyz_1, g_xxyyyy_0_xxxxxxz_1, g_xxyyyy_0_xxxxxxzz_1, g_xxyyyy_0_xxxxxyy_1, g_xxyyyy_0_xxxxxyyy_1, g_xxyyyy_0_xxxxxyyz_1, g_xxyyyy_0_xxxxxyz_1, g_xxyyyy_0_xxxxxyzz_1, g_xxyyyy_0_xxxxxzz_1, g_xxyyyy_0_xxxxxzzz_1, g_xxyyyy_0_xxxxyyy_1, g_xxyyyy_0_xxxxyyyy_1, g_xxyyyy_0_xxxxyyyz_1, g_xxyyyy_0_xxxxyyz_1, g_xxyyyy_0_xxxxyyzz_1, g_xxyyyy_0_xxxxyzz_1, g_xxyyyy_0_xxxxyzzz_1, g_xxyyyy_0_xxxxzzz_1, g_xxyyyy_0_xxxxzzzz_1, g_xxyyyy_0_xxxyyyy_1, g_xxyyyy_0_xxxyyyyy_1, g_xxyyyy_0_xxxyyyyz_1, g_xxyyyy_0_xxxyyyz_1, g_xxyyyy_0_xxxyyyzz_1, g_xxyyyy_0_xxxyyzz_1, g_xxyyyy_0_xxxyyzzz_1, g_xxyyyy_0_xxxyzzz_1, g_xxyyyy_0_xxxyzzzz_1, g_xxyyyy_0_xxxzzzz_1, g_xxyyyy_0_xxxzzzzz_1, g_xxyyyy_0_xxyyyyy_1, g_xxyyyy_0_xxyyyyyy_1, g_xxyyyy_0_xxyyyyyz_1, g_xxyyyy_0_xxyyyyz_1, g_xxyyyy_0_xxyyyyzz_1, g_xxyyyy_0_xxyyyzz_1, g_xxyyyy_0_xxyyyzzz_1, g_xxyyyy_0_xxyyzzz_1, g_xxyyyy_0_xxyyzzzz_1, g_xxyyyy_0_xxyzzzz_1, g_xxyyyy_0_xxyzzzzz_1, g_xxyyyy_0_xxzzzzz_1, g_xxyyyy_0_xxzzzzzz_1, g_xxyyyy_0_xyyyyyy_1, g_xxyyyy_0_xyyyyyyy_1, g_xxyyyy_0_xyyyyyyz_1, g_xxyyyy_0_xyyyyyz_1, g_xxyyyy_0_xyyyyyzz_1, g_xxyyyy_0_xyyyyzz_1, g_xxyyyy_0_xyyyyzzz_1, g_xxyyyy_0_xyyyzzz_1, g_xxyyyy_0_xyyyzzzz_1, g_xxyyyy_0_xyyzzzz_1, g_xxyyyy_0_xyyzzzzz_1, g_xxyyyy_0_xyzzzzz_1, g_xxyyyy_0_xyzzzzzz_1, g_xxyyyy_0_xzzzzzz_1, g_xxyyyy_0_xzzzzzzz_1, g_xxyyyy_0_yyyyyyy_1, g_xxyyyy_0_yyyyyyyy_1, g_xxyyyy_0_yyyyyyyz_1, g_xxyyyy_0_yyyyyyz_1, g_xxyyyy_0_yyyyyyzz_1, g_xxyyyy_0_yyyyyzz_1, g_xxyyyy_0_yyyyyzzz_1, g_xxyyyy_0_yyyyzzz_1, g_xxyyyy_0_yyyyzzzz_1, g_xxyyyy_0_yyyzzzz_1, g_xxyyyy_0_yyyzzzzz_1, g_xxyyyy_0_yyzzzzz_1, g_xxyyyy_0_yyzzzzzz_1, g_xxyyyy_0_yzzzzzz_1, g_xxyyyy_0_yzzzzzzz_1, g_xxyyyy_0_zzzzzzz_1, g_xxyyyy_0_zzzzzzzz_1, g_xxyyyyz_0_xxxxxxxx_0, g_xxyyyyz_0_xxxxxxxy_0, g_xxyyyyz_0_xxxxxxxz_0, g_xxyyyyz_0_xxxxxxyy_0, g_xxyyyyz_0_xxxxxxyz_0, g_xxyyyyz_0_xxxxxxzz_0, g_xxyyyyz_0_xxxxxyyy_0, g_xxyyyyz_0_xxxxxyyz_0, g_xxyyyyz_0_xxxxxyzz_0, g_xxyyyyz_0_xxxxxzzz_0, g_xxyyyyz_0_xxxxyyyy_0, g_xxyyyyz_0_xxxxyyyz_0, g_xxyyyyz_0_xxxxyyzz_0, g_xxyyyyz_0_xxxxyzzz_0, g_xxyyyyz_0_xxxxzzzz_0, g_xxyyyyz_0_xxxyyyyy_0, g_xxyyyyz_0_xxxyyyyz_0, g_xxyyyyz_0_xxxyyyzz_0, g_xxyyyyz_0_xxxyyzzz_0, g_xxyyyyz_0_xxxyzzzz_0, g_xxyyyyz_0_xxxzzzzz_0, g_xxyyyyz_0_xxyyyyyy_0, g_xxyyyyz_0_xxyyyyyz_0, g_xxyyyyz_0_xxyyyyzz_0, g_xxyyyyz_0_xxyyyzzz_0, g_xxyyyyz_0_xxyyzzzz_0, g_xxyyyyz_0_xxyzzzzz_0, g_xxyyyyz_0_xxzzzzzz_0, g_xxyyyyz_0_xyyyyyyy_0, g_xxyyyyz_0_xyyyyyyz_0, g_xxyyyyz_0_xyyyyyzz_0, g_xxyyyyz_0_xyyyyzzz_0, g_xxyyyyz_0_xyyyzzzz_0, g_xxyyyyz_0_xyyzzzzz_0, g_xxyyyyz_0_xyzzzzzz_0, g_xxyyyyz_0_xzzzzzzz_0, g_xxyyyyz_0_yyyyyyyy_0, g_xxyyyyz_0_yyyyyyyz_0, g_xxyyyyz_0_yyyyyyzz_0, g_xxyyyyz_0_yyyyyzzz_0, g_xxyyyyz_0_yyyyzzzz_0, g_xxyyyyz_0_yyyzzzzz_0, g_xxyyyyz_0_yyzzzzzz_0, g_xxyyyyz_0_yzzzzzzz_0, g_xxyyyyz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyyz_0_xxxxxxxx_0[i] = g_xxyyyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxxxy_0[i] = g_xxyyyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxxxz_0[i] = g_xxyyyy_0_xxxxxxx_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxxxz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxxyy_0[i] = g_xxyyyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxxyz_0[i] = g_xxyyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxxyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxxzz_0[i] = 2.0 * g_xxyyyy_0_xxxxxxz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxxzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxyyy_0[i] = g_xxyyyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxyyz_0[i] = g_xxyyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxyzz_0[i] = 2.0 * g_xxyyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxzzz_0[i] = 3.0 * g_xxyyyy_0_xxxxxzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxyyyy_0[i] = g_xxyyyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxyyyz_0[i] = g_xxyyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxyyzz_0[i] = 2.0 * g_xxyyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxyzzz_0[i] = 3.0 * g_xxyyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxzzzz_0[i] = 4.0 * g_xxyyyy_0_xxxxzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyyyyy_0[i] = g_xxyyyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyyyyz_0[i] = g_xxyyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyyyzz_0[i] = 2.0 * g_xxyyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyyzzz_0[i] = 3.0 * g_xxyyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyzzzz_0[i] = 4.0 * g_xxyyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxzzzzz_0[i] = 5.0 * g_xxyyyy_0_xxxzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyyyyy_0[i] = g_xxyyyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyyyyz_0[i] = g_xxyyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyyyzz_0[i] = 2.0 * g_xxyyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyyzzz_0[i] = 3.0 * g_xxyyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyzzzz_0[i] = 4.0 * g_xxyyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyzzzzz_0[i] = 5.0 * g_xxyyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxzzzzzz_0[i] = 6.0 * g_xxyyyy_0_xxzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxzzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyyyyy_0[i] = g_xxyyyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyyyyz_0[i] = g_xxyyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyyyzz_0[i] = 2.0 * g_xxyyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyyzzz_0[i] = 3.0 * g_xxyyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyzzzz_0[i] = 4.0 * g_xxyyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyzzzzz_0[i] = 5.0 * g_xxyyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyzzzzzz_0[i] = 6.0 * g_xxyyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyzzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xzzzzzzz_0[i] = 7.0 * g_xxyyyy_0_xzzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xzzzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyyyyy_0[i] = g_xxyyyy_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyyyyz_0[i] = g_xxyyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_yyyyyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyyyzz_0[i] = 2.0 * g_xxyyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_yyyyyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyyzzz_0[i] = 3.0 * g_xxyyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_yyyyyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyzzzz_0[i] = 4.0 * g_xxyyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_yyyyzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyzzzzz_0[i] = 5.0 * g_xxyyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_yyyzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyzzzzzz_0[i] = 6.0 * g_xxyyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_yyzzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yzzzzzzz_0[i] = 7.0 * g_xxyyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_yzzzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_zzzzzzzz_0[i] = 8.0 * g_xxyyyy_0_zzzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 765-810 components of targeted buffer : KSL

    auto g_xxyyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 765);

    auto g_xxyyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 766);

    auto g_xxyyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 767);

    auto g_xxyyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 768);

    auto g_xxyyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 769);

    auto g_xxyyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 770);

    auto g_xxyyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 771);

    auto g_xxyyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 772);

    auto g_xxyyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 773);

    auto g_xxyyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 774);

    auto g_xxyyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 775);

    auto g_xxyyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 776);

    auto g_xxyyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 777);

    auto g_xxyyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 778);

    auto g_xxyyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 779);

    auto g_xxyyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 780);

    auto g_xxyyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 781);

    auto g_xxyyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 782);

    auto g_xxyyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 783);

    auto g_xxyyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 784);

    auto g_xxyyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 785);

    auto g_xxyyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 786);

    auto g_xxyyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 787);

    auto g_xxyyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 788);

    auto g_xxyyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 789);

    auto g_xxyyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 790);

    auto g_xxyyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 791);

    auto g_xxyyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 792);

    auto g_xxyyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 793);

    auto g_xxyyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 794);

    auto g_xxyyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 795);

    auto g_xxyyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 796);

    auto g_xxyyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 797);

    auto g_xxyyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 798);

    auto g_xxyyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 799);

    auto g_xxyyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 800);

    auto g_xxyyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 801);

    auto g_xxyyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 802);

    auto g_xxyyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 803);

    auto g_xxyyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 804);

    auto g_xxyyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 805);

    auto g_xxyyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 806);

    auto g_xxyyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 807);

    auto g_xxyyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 808);

    auto g_xxyyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 809);

    #pragma omp simd aligned(g_xxyyy_0_xxxxxxxy_0, g_xxyyy_0_xxxxxxxy_1, g_xxyyy_0_xxxxxxyy_0, g_xxyyy_0_xxxxxxyy_1, g_xxyyy_0_xxxxxyyy_0, g_xxyyy_0_xxxxxyyy_1, g_xxyyy_0_xxxxyyyy_0, g_xxyyy_0_xxxxyyyy_1, g_xxyyy_0_xxxyyyyy_0, g_xxyyy_0_xxxyyyyy_1, g_xxyyy_0_xxyyyyyy_0, g_xxyyy_0_xxyyyyyy_1, g_xxyyy_0_xyyyyyyy_0, g_xxyyy_0_xyyyyyyy_1, g_xxyyyz_0_xxxxxxxy_1, g_xxyyyz_0_xxxxxxyy_1, g_xxyyyz_0_xxxxxyyy_1, g_xxyyyz_0_xxxxyyyy_1, g_xxyyyz_0_xxxyyyyy_1, g_xxyyyz_0_xxyyyyyy_1, g_xxyyyz_0_xyyyyyyy_1, g_xxyyyzz_0_xxxxxxxx_0, g_xxyyyzz_0_xxxxxxxy_0, g_xxyyyzz_0_xxxxxxxz_0, g_xxyyyzz_0_xxxxxxyy_0, g_xxyyyzz_0_xxxxxxyz_0, g_xxyyyzz_0_xxxxxxzz_0, g_xxyyyzz_0_xxxxxyyy_0, g_xxyyyzz_0_xxxxxyyz_0, g_xxyyyzz_0_xxxxxyzz_0, g_xxyyyzz_0_xxxxxzzz_0, g_xxyyyzz_0_xxxxyyyy_0, g_xxyyyzz_0_xxxxyyyz_0, g_xxyyyzz_0_xxxxyyzz_0, g_xxyyyzz_0_xxxxyzzz_0, g_xxyyyzz_0_xxxxzzzz_0, g_xxyyyzz_0_xxxyyyyy_0, g_xxyyyzz_0_xxxyyyyz_0, g_xxyyyzz_0_xxxyyyzz_0, g_xxyyyzz_0_xxxyyzzz_0, g_xxyyyzz_0_xxxyzzzz_0, g_xxyyyzz_0_xxxzzzzz_0, g_xxyyyzz_0_xxyyyyyy_0, g_xxyyyzz_0_xxyyyyyz_0, g_xxyyyzz_0_xxyyyyzz_0, g_xxyyyzz_0_xxyyyzzz_0, g_xxyyyzz_0_xxyyzzzz_0, g_xxyyyzz_0_xxyzzzzz_0, g_xxyyyzz_0_xxzzzzzz_0, g_xxyyyzz_0_xyyyyyyy_0, g_xxyyyzz_0_xyyyyyyz_0, g_xxyyyzz_0_xyyyyyzz_0, g_xxyyyzz_0_xyyyyzzz_0, g_xxyyyzz_0_xyyyzzzz_0, g_xxyyyzz_0_xyyzzzzz_0, g_xxyyyzz_0_xyzzzzzz_0, g_xxyyyzz_0_xzzzzzzz_0, g_xxyyyzz_0_yyyyyyyy_0, g_xxyyyzz_0_yyyyyyyz_0, g_xxyyyzz_0_yyyyyyzz_0, g_xxyyyzz_0_yyyyyzzz_0, g_xxyyyzz_0_yyyyzzzz_0, g_xxyyyzz_0_yyyzzzzz_0, g_xxyyyzz_0_yyzzzzzz_0, g_xxyyyzz_0_yzzzzzzz_0, g_xxyyyzz_0_zzzzzzzz_0, g_xxyyzz_0_xxxxxxxx_1, g_xxyyzz_0_xxxxxxxz_1, g_xxyyzz_0_xxxxxxzz_1, g_xxyyzz_0_xxxxxzzz_1, g_xxyyzz_0_xxxxzzzz_1, g_xxyyzz_0_xxxzzzzz_1, g_xxyyzz_0_xxzzzzzz_1, g_xxyyzz_0_xzzzzzzz_1, g_xxyzz_0_xxxxxxxx_0, g_xxyzz_0_xxxxxxxx_1, g_xxyzz_0_xxxxxxxz_0, g_xxyzz_0_xxxxxxxz_1, g_xxyzz_0_xxxxxxzz_0, g_xxyzz_0_xxxxxxzz_1, g_xxyzz_0_xxxxxzzz_0, g_xxyzz_0_xxxxxzzz_1, g_xxyzz_0_xxxxzzzz_0, g_xxyzz_0_xxxxzzzz_1, g_xxyzz_0_xxxzzzzz_0, g_xxyzz_0_xxxzzzzz_1, g_xxyzz_0_xxzzzzzz_0, g_xxyzz_0_xxzzzzzz_1, g_xxyzz_0_xzzzzzzz_0, g_xxyzz_0_xzzzzzzz_1, g_xyyyzz_0_xxxxxxyz_1, g_xyyyzz_0_xxxxxyyz_1, g_xyyyzz_0_xxxxxyz_1, g_xyyyzz_0_xxxxxyzz_1, g_xyyyzz_0_xxxxyyyz_1, g_xyyyzz_0_xxxxyyz_1, g_xyyyzz_0_xxxxyyzz_1, g_xyyyzz_0_xxxxyzz_1, g_xyyyzz_0_xxxxyzzz_1, g_xyyyzz_0_xxxyyyyz_1, g_xyyyzz_0_xxxyyyz_1, g_xyyyzz_0_xxxyyyzz_1, g_xyyyzz_0_xxxyyzz_1, g_xyyyzz_0_xxxyyzzz_1, g_xyyyzz_0_xxxyzzz_1, g_xyyyzz_0_xxxyzzzz_1, g_xyyyzz_0_xxyyyyyz_1, g_xyyyzz_0_xxyyyyz_1, g_xyyyzz_0_xxyyyyzz_1, g_xyyyzz_0_xxyyyzz_1, g_xyyyzz_0_xxyyyzzz_1, g_xyyyzz_0_xxyyzzz_1, g_xyyyzz_0_xxyyzzzz_1, g_xyyyzz_0_xxyzzzz_1, g_xyyyzz_0_xxyzzzzz_1, g_xyyyzz_0_xyyyyyyz_1, g_xyyyzz_0_xyyyyyz_1, g_xyyyzz_0_xyyyyyzz_1, g_xyyyzz_0_xyyyyzz_1, g_xyyyzz_0_xyyyyzzz_1, g_xyyyzz_0_xyyyzzz_1, g_xyyyzz_0_xyyyzzzz_1, g_xyyyzz_0_xyyzzzz_1, g_xyyyzz_0_xyyzzzzz_1, g_xyyyzz_0_xyzzzzz_1, g_xyyyzz_0_xyzzzzzz_1, g_xyyyzz_0_yyyyyyyy_1, g_xyyyzz_0_yyyyyyyz_1, g_xyyyzz_0_yyyyyyz_1, g_xyyyzz_0_yyyyyyzz_1, g_xyyyzz_0_yyyyyzz_1, g_xyyyzz_0_yyyyyzzz_1, g_xyyyzz_0_yyyyzzz_1, g_xyyyzz_0_yyyyzzzz_1, g_xyyyzz_0_yyyzzzz_1, g_xyyyzz_0_yyyzzzzz_1, g_xyyyzz_0_yyzzzzz_1, g_xyyyzz_0_yyzzzzzz_1, g_xyyyzz_0_yzzzzzz_1, g_xyyyzz_0_yzzzzzzz_1, g_xyyyzz_0_zzzzzzzz_1, g_yyyzz_0_xxxxxxyz_0, g_yyyzz_0_xxxxxxyz_1, g_yyyzz_0_xxxxxyyz_0, g_yyyzz_0_xxxxxyyz_1, g_yyyzz_0_xxxxxyzz_0, g_yyyzz_0_xxxxxyzz_1, g_yyyzz_0_xxxxyyyz_0, g_yyyzz_0_xxxxyyyz_1, g_yyyzz_0_xxxxyyzz_0, g_yyyzz_0_xxxxyyzz_1, g_yyyzz_0_xxxxyzzz_0, g_yyyzz_0_xxxxyzzz_1, g_yyyzz_0_xxxyyyyz_0, g_yyyzz_0_xxxyyyyz_1, g_yyyzz_0_xxxyyyzz_0, g_yyyzz_0_xxxyyyzz_1, g_yyyzz_0_xxxyyzzz_0, g_yyyzz_0_xxxyyzzz_1, g_yyyzz_0_xxxyzzzz_0, g_yyyzz_0_xxxyzzzz_1, g_yyyzz_0_xxyyyyyz_0, g_yyyzz_0_xxyyyyyz_1, g_yyyzz_0_xxyyyyzz_0, g_yyyzz_0_xxyyyyzz_1, g_yyyzz_0_xxyyyzzz_0, g_yyyzz_0_xxyyyzzz_1, g_yyyzz_0_xxyyzzzz_0, g_yyyzz_0_xxyyzzzz_1, g_yyyzz_0_xxyzzzzz_0, g_yyyzz_0_xxyzzzzz_1, g_yyyzz_0_xyyyyyyz_0, g_yyyzz_0_xyyyyyyz_1, g_yyyzz_0_xyyyyyzz_0, g_yyyzz_0_xyyyyyzz_1, g_yyyzz_0_xyyyyzzz_0, g_yyyzz_0_xyyyyzzz_1, g_yyyzz_0_xyyyzzzz_0, g_yyyzz_0_xyyyzzzz_1, g_yyyzz_0_xyyzzzzz_0, g_yyyzz_0_xyyzzzzz_1, g_yyyzz_0_xyzzzzzz_0, g_yyyzz_0_xyzzzzzz_1, g_yyyzz_0_yyyyyyyy_0, g_yyyzz_0_yyyyyyyy_1, g_yyyzz_0_yyyyyyyz_0, g_yyyzz_0_yyyyyyyz_1, g_yyyzz_0_yyyyyyzz_0, g_yyyzz_0_yyyyyyzz_1, g_yyyzz_0_yyyyyzzz_0, g_yyyzz_0_yyyyyzzz_1, g_yyyzz_0_yyyyzzzz_0, g_yyyzz_0_yyyyzzzz_1, g_yyyzz_0_yyyzzzzz_0, g_yyyzz_0_yyyzzzzz_1, g_yyyzz_0_yyzzzzzz_0, g_yyyzz_0_yyzzzzzz_1, g_yyyzz_0_yzzzzzzz_0, g_yyyzz_0_yzzzzzzz_1, g_yyyzz_0_zzzzzzzz_0, g_yyyzz_0_zzzzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyzz_0_xxxxxxxx_0[i] = 2.0 * g_xxyzz_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxxxxx_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxxxxxy_0[i] = g_xxyyy_0_xxxxxxxy_0[i] * fbe_0 - g_xxyyy_0_xxxxxxxy_1[i] * fz_be_0 + g_xxyyyz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxxxxxz_0[i] = 2.0 * g_xxyzz_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxxxxz_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxxxxyy_0[i] = g_xxyyy_0_xxxxxxyy_0[i] * fbe_0 - g_xxyyy_0_xxxxxxyy_1[i] * fz_be_0 + g_xxyyyz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxxxxyz_0[i] = g_yyyzz_0_xxxxxxyz_0[i] * fbe_0 - g_yyyzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xyyyzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxxxxzz_0[i] = 2.0 * g_xxyzz_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxxxzz_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxxxyyy_0[i] = g_xxyyy_0_xxxxxyyy_0[i] * fbe_0 - g_xxyyy_0_xxxxxyyy_1[i] * fz_be_0 + g_xxyyyz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxxxyyz_0[i] = g_yyyzz_0_xxxxxyyz_0[i] * fbe_0 - g_yyyzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xyyyzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxxxyzz_0[i] = g_yyyzz_0_xxxxxyzz_0[i] * fbe_0 - g_yyyzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xyyyzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxxxzzz_0[i] = 2.0 * g_xxyzz_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxxzzz_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxxyyyy_0[i] = g_xxyyy_0_xxxxyyyy_0[i] * fbe_0 - g_xxyyy_0_xxxxyyyy_1[i] * fz_be_0 + g_xxyyyz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxxyyyz_0[i] = g_yyyzz_0_xxxxyyyz_0[i] * fbe_0 - g_yyyzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xyyyzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxxyyzz_0[i] = g_yyyzz_0_xxxxyyzz_0[i] * fbe_0 - g_yyyzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xyyyzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxxyzzz_0[i] = g_yyyzz_0_xxxxyzzz_0[i] * fbe_0 - g_yyyzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xyyyzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxxzzzz_0[i] = 2.0 * g_xxyzz_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxzzzz_1[i] * fz_be_0 + g_xxyyzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxyyyyy_0[i] = g_xxyyy_0_xxxyyyyy_0[i] * fbe_0 - g_xxyyy_0_xxxyyyyy_1[i] * fz_be_0 + g_xxyyyz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxyyyyz_0[i] = g_yyyzz_0_xxxyyyyz_0[i] * fbe_0 - g_yyyzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxyyyzz_0[i] = g_yyyzz_0_xxxyyyzz_0[i] * fbe_0 - g_yyyzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxyyzzz_0[i] = g_yyyzz_0_xxxyyzzz_0[i] * fbe_0 - g_yyyzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxyzzzz_0[i] = g_yyyzz_0_xxxyzzzz_0[i] * fbe_0 - g_yyyzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxzzzzz_0[i] = 2.0 * g_xxyzz_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxyyyyyy_0[i] = g_xxyyy_0_xxyyyyyy_0[i] * fbe_0 - g_xxyyy_0_xxyyyyyy_1[i] * fz_be_0 + g_xxyyyz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxyyyyyz_0[i] = g_yyyzz_0_xxyyyyyz_0[i] * fbe_0 - g_yyyzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxyyyyzz_0[i] = g_yyyzz_0_xxyyyyzz_0[i] * fbe_0 - g_yyyzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxyyyzzz_0[i] = g_yyyzz_0_xxyyyzzz_0[i] * fbe_0 - g_yyyzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxyyzzzz_0[i] = g_yyyzz_0_xxyyzzzz_0[i] * fbe_0 - g_yyyzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxyzzzzz_0[i] = g_yyyzz_0_xxyzzzzz_0[i] * fbe_0 - g_yyyzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxzzzzzz_0[i] = 2.0 * g_xxyzz_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxzzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xyyyyyyy_0[i] = g_xxyyy_0_xyyyyyyy_0[i] * fbe_0 - g_xxyyy_0_xyyyyyyy_1[i] * fz_be_0 + g_xxyyyz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xyyyyyyz_0[i] = g_yyyzz_0_xyyyyyyz_0[i] * fbe_0 - g_yyyzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyyyyyzz_0[i] = g_yyyzz_0_xyyyyyzz_0[i] * fbe_0 - g_yyyzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyyyyzzz_0[i] = g_yyyzz_0_xyyyyzzz_0[i] * fbe_0 - g_yyyzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyyyzzzz_0[i] = g_yyyzz_0_xyyyzzzz_0[i] * fbe_0 - g_yyyzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyyzzzzz_0[i] = g_yyyzz_0_xyyzzzzz_0[i] * fbe_0 - g_yyyzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyzzzzzz_0[i] = g_yyyzz_0_xyzzzzzz_0[i] * fbe_0 - g_yyyzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xzzzzzzz_0[i] = 2.0 * g_xxyzz_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_yyyyyyyy_0[i] = g_yyyzz_0_yyyyyyyy_0[i] * fbe_0 - g_yyyzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xyyyzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyyyyyz_0[i] = g_yyyzz_0_yyyyyyyz_0[i] * fbe_0 - g_yyyzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyyyyzz_0[i] = g_yyyzz_0_yyyyyyzz_0[i] * fbe_0 - g_yyyzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyyyzzz_0[i] = g_yyyzz_0_yyyyyzzz_0[i] * fbe_0 - g_yyyzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyyzzzz_0[i] = g_yyyzz_0_yyyyzzzz_0[i] * fbe_0 - g_yyyzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyzzzzz_0[i] = g_yyyzz_0_yyyzzzzz_0[i] * fbe_0 - g_yyyzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyzzzzzz_0[i] = g_yyyzz_0_yyzzzzzz_0[i] * fbe_0 - g_yyyzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yzzzzzzz_0[i] = g_yyyzz_0_yzzzzzzz_0[i] * fbe_0 - g_yyyzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_zzzzzzzz_0[i] = g_yyyzz_0_zzzzzzzz_0[i] * fbe_0 - g_yyyzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xyyyzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 810-855 components of targeted buffer : KSL

    auto g_xxyyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 810);

    auto g_xxyyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 811);

    auto g_xxyyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 812);

    auto g_xxyyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 813);

    auto g_xxyyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 814);

    auto g_xxyyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 815);

    auto g_xxyyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 816);

    auto g_xxyyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 817);

    auto g_xxyyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 818);

    auto g_xxyyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 819);

    auto g_xxyyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 820);

    auto g_xxyyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 821);

    auto g_xxyyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 822);

    auto g_xxyyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 823);

    auto g_xxyyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 824);

    auto g_xxyyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 825);

    auto g_xxyyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 826);

    auto g_xxyyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 827);

    auto g_xxyyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 828);

    auto g_xxyyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 829);

    auto g_xxyyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 830);

    auto g_xxyyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 831);

    auto g_xxyyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 832);

    auto g_xxyyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 833);

    auto g_xxyyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 834);

    auto g_xxyyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 835);

    auto g_xxyyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 836);

    auto g_xxyyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 837);

    auto g_xxyyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 838);

    auto g_xxyyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 839);

    auto g_xxyyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 840);

    auto g_xxyyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 841);

    auto g_xxyyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 842);

    auto g_xxyyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 843);

    auto g_xxyyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 844);

    auto g_xxyyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 845);

    auto g_xxyyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 846);

    auto g_xxyyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 847);

    auto g_xxyyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 848);

    auto g_xxyyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 849);

    auto g_xxyyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 850);

    auto g_xxyyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 851);

    auto g_xxyyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 852);

    auto g_xxyyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 853);

    auto g_xxyyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 854);

    #pragma omp simd aligned(g_xxyyz_0_xxxxxxxy_0, g_xxyyz_0_xxxxxxxy_1, g_xxyyz_0_xxxxxxyy_0, g_xxyyz_0_xxxxxxyy_1, g_xxyyz_0_xxxxxyyy_0, g_xxyyz_0_xxxxxyyy_1, g_xxyyz_0_xxxxyyyy_0, g_xxyyz_0_xxxxyyyy_1, g_xxyyz_0_xxxyyyyy_0, g_xxyyz_0_xxxyyyyy_1, g_xxyyz_0_xxyyyyyy_0, g_xxyyz_0_xxyyyyyy_1, g_xxyyz_0_xyyyyyyy_0, g_xxyyz_0_xyyyyyyy_1, g_xxyyzz_0_xxxxxxxy_1, g_xxyyzz_0_xxxxxxyy_1, g_xxyyzz_0_xxxxxyyy_1, g_xxyyzz_0_xxxxyyyy_1, g_xxyyzz_0_xxxyyyyy_1, g_xxyyzz_0_xxyyyyyy_1, g_xxyyzz_0_xyyyyyyy_1, g_xxyyzzz_0_xxxxxxxx_0, g_xxyyzzz_0_xxxxxxxy_0, g_xxyyzzz_0_xxxxxxxz_0, g_xxyyzzz_0_xxxxxxyy_0, g_xxyyzzz_0_xxxxxxyz_0, g_xxyyzzz_0_xxxxxxzz_0, g_xxyyzzz_0_xxxxxyyy_0, g_xxyyzzz_0_xxxxxyyz_0, g_xxyyzzz_0_xxxxxyzz_0, g_xxyyzzz_0_xxxxxzzz_0, g_xxyyzzz_0_xxxxyyyy_0, g_xxyyzzz_0_xxxxyyyz_0, g_xxyyzzz_0_xxxxyyzz_0, g_xxyyzzz_0_xxxxyzzz_0, g_xxyyzzz_0_xxxxzzzz_0, g_xxyyzzz_0_xxxyyyyy_0, g_xxyyzzz_0_xxxyyyyz_0, g_xxyyzzz_0_xxxyyyzz_0, g_xxyyzzz_0_xxxyyzzz_0, g_xxyyzzz_0_xxxyzzzz_0, g_xxyyzzz_0_xxxzzzzz_0, g_xxyyzzz_0_xxyyyyyy_0, g_xxyyzzz_0_xxyyyyyz_0, g_xxyyzzz_0_xxyyyyzz_0, g_xxyyzzz_0_xxyyyzzz_0, g_xxyyzzz_0_xxyyzzzz_0, g_xxyyzzz_0_xxyzzzzz_0, g_xxyyzzz_0_xxzzzzzz_0, g_xxyyzzz_0_xyyyyyyy_0, g_xxyyzzz_0_xyyyyyyz_0, g_xxyyzzz_0_xyyyyyzz_0, g_xxyyzzz_0_xyyyyzzz_0, g_xxyyzzz_0_xyyyzzzz_0, g_xxyyzzz_0_xyyzzzzz_0, g_xxyyzzz_0_xyzzzzzz_0, g_xxyyzzz_0_xzzzzzzz_0, g_xxyyzzz_0_yyyyyyyy_0, g_xxyyzzz_0_yyyyyyyz_0, g_xxyyzzz_0_yyyyyyzz_0, g_xxyyzzz_0_yyyyyzzz_0, g_xxyyzzz_0_yyyyzzzz_0, g_xxyyzzz_0_yyyzzzzz_0, g_xxyyzzz_0_yyzzzzzz_0, g_xxyyzzz_0_yzzzzzzz_0, g_xxyyzzz_0_zzzzzzzz_0, g_xxyzzz_0_xxxxxxxx_1, g_xxyzzz_0_xxxxxxxz_1, g_xxyzzz_0_xxxxxxzz_1, g_xxyzzz_0_xxxxxzzz_1, g_xxyzzz_0_xxxxzzzz_1, g_xxyzzz_0_xxxzzzzz_1, g_xxyzzz_0_xxzzzzzz_1, g_xxyzzz_0_xzzzzzzz_1, g_xxzzz_0_xxxxxxxx_0, g_xxzzz_0_xxxxxxxx_1, g_xxzzz_0_xxxxxxxz_0, g_xxzzz_0_xxxxxxxz_1, g_xxzzz_0_xxxxxxzz_0, g_xxzzz_0_xxxxxxzz_1, g_xxzzz_0_xxxxxzzz_0, g_xxzzz_0_xxxxxzzz_1, g_xxzzz_0_xxxxzzzz_0, g_xxzzz_0_xxxxzzzz_1, g_xxzzz_0_xxxzzzzz_0, g_xxzzz_0_xxxzzzzz_1, g_xxzzz_0_xxzzzzzz_0, g_xxzzz_0_xxzzzzzz_1, g_xxzzz_0_xzzzzzzz_0, g_xxzzz_0_xzzzzzzz_1, g_xyyzzz_0_xxxxxxyz_1, g_xyyzzz_0_xxxxxyyz_1, g_xyyzzz_0_xxxxxyz_1, g_xyyzzz_0_xxxxxyzz_1, g_xyyzzz_0_xxxxyyyz_1, g_xyyzzz_0_xxxxyyz_1, g_xyyzzz_0_xxxxyyzz_1, g_xyyzzz_0_xxxxyzz_1, g_xyyzzz_0_xxxxyzzz_1, g_xyyzzz_0_xxxyyyyz_1, g_xyyzzz_0_xxxyyyz_1, g_xyyzzz_0_xxxyyyzz_1, g_xyyzzz_0_xxxyyzz_1, g_xyyzzz_0_xxxyyzzz_1, g_xyyzzz_0_xxxyzzz_1, g_xyyzzz_0_xxxyzzzz_1, g_xyyzzz_0_xxyyyyyz_1, g_xyyzzz_0_xxyyyyz_1, g_xyyzzz_0_xxyyyyzz_1, g_xyyzzz_0_xxyyyzz_1, g_xyyzzz_0_xxyyyzzz_1, g_xyyzzz_0_xxyyzzz_1, g_xyyzzz_0_xxyyzzzz_1, g_xyyzzz_0_xxyzzzz_1, g_xyyzzz_0_xxyzzzzz_1, g_xyyzzz_0_xyyyyyyz_1, g_xyyzzz_0_xyyyyyz_1, g_xyyzzz_0_xyyyyyzz_1, g_xyyzzz_0_xyyyyzz_1, g_xyyzzz_0_xyyyyzzz_1, g_xyyzzz_0_xyyyzzz_1, g_xyyzzz_0_xyyyzzzz_1, g_xyyzzz_0_xyyzzzz_1, g_xyyzzz_0_xyyzzzzz_1, g_xyyzzz_0_xyzzzzz_1, g_xyyzzz_0_xyzzzzzz_1, g_xyyzzz_0_yyyyyyyy_1, g_xyyzzz_0_yyyyyyyz_1, g_xyyzzz_0_yyyyyyz_1, g_xyyzzz_0_yyyyyyzz_1, g_xyyzzz_0_yyyyyzz_1, g_xyyzzz_0_yyyyyzzz_1, g_xyyzzz_0_yyyyzzz_1, g_xyyzzz_0_yyyyzzzz_1, g_xyyzzz_0_yyyzzzz_1, g_xyyzzz_0_yyyzzzzz_1, g_xyyzzz_0_yyzzzzz_1, g_xyyzzz_0_yyzzzzzz_1, g_xyyzzz_0_yzzzzzz_1, g_xyyzzz_0_yzzzzzzz_1, g_xyyzzz_0_zzzzzzzz_1, g_yyzzz_0_xxxxxxyz_0, g_yyzzz_0_xxxxxxyz_1, g_yyzzz_0_xxxxxyyz_0, g_yyzzz_0_xxxxxyyz_1, g_yyzzz_0_xxxxxyzz_0, g_yyzzz_0_xxxxxyzz_1, g_yyzzz_0_xxxxyyyz_0, g_yyzzz_0_xxxxyyyz_1, g_yyzzz_0_xxxxyyzz_0, g_yyzzz_0_xxxxyyzz_1, g_yyzzz_0_xxxxyzzz_0, g_yyzzz_0_xxxxyzzz_1, g_yyzzz_0_xxxyyyyz_0, g_yyzzz_0_xxxyyyyz_1, g_yyzzz_0_xxxyyyzz_0, g_yyzzz_0_xxxyyyzz_1, g_yyzzz_0_xxxyyzzz_0, g_yyzzz_0_xxxyyzzz_1, g_yyzzz_0_xxxyzzzz_0, g_yyzzz_0_xxxyzzzz_1, g_yyzzz_0_xxyyyyyz_0, g_yyzzz_0_xxyyyyyz_1, g_yyzzz_0_xxyyyyzz_0, g_yyzzz_0_xxyyyyzz_1, g_yyzzz_0_xxyyyzzz_0, g_yyzzz_0_xxyyyzzz_1, g_yyzzz_0_xxyyzzzz_0, g_yyzzz_0_xxyyzzzz_1, g_yyzzz_0_xxyzzzzz_0, g_yyzzz_0_xxyzzzzz_1, g_yyzzz_0_xyyyyyyz_0, g_yyzzz_0_xyyyyyyz_1, g_yyzzz_0_xyyyyyzz_0, g_yyzzz_0_xyyyyyzz_1, g_yyzzz_0_xyyyyzzz_0, g_yyzzz_0_xyyyyzzz_1, g_yyzzz_0_xyyyzzzz_0, g_yyzzz_0_xyyyzzzz_1, g_yyzzz_0_xyyzzzzz_0, g_yyzzz_0_xyyzzzzz_1, g_yyzzz_0_xyzzzzzz_0, g_yyzzz_0_xyzzzzzz_1, g_yyzzz_0_yyyyyyyy_0, g_yyzzz_0_yyyyyyyy_1, g_yyzzz_0_yyyyyyyz_0, g_yyzzz_0_yyyyyyyz_1, g_yyzzz_0_yyyyyyzz_0, g_yyzzz_0_yyyyyyzz_1, g_yyzzz_0_yyyyyzzz_0, g_yyzzz_0_yyyyyzzz_1, g_yyzzz_0_yyyyzzzz_0, g_yyzzz_0_yyyyzzzz_1, g_yyzzz_0_yyyzzzzz_0, g_yyzzz_0_yyyzzzzz_1, g_yyzzz_0_yyzzzzzz_0, g_yyzzz_0_yyzzzzzz_1, g_yyzzz_0_yzzzzzzz_0, g_yyzzz_0_yzzzzzzz_1, g_yyzzz_0_zzzzzzzz_0, g_yyzzz_0_zzzzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyzzz_0_xxxxxxxx_0[i] = g_xxzzz_0_xxxxxxxx_0[i] * fbe_0 - g_xxzzz_0_xxxxxxxx_1[i] * fz_be_0 + g_xxyzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxxxxxy_0[i] = 2.0 * g_xxyyz_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxxxxxy_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxxxxxz_0[i] = g_xxzzz_0_xxxxxxxz_0[i] * fbe_0 - g_xxzzz_0_xxxxxxxz_1[i] * fz_be_0 + g_xxyzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxxxxyy_0[i] = 2.0 * g_xxyyz_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxxxxyy_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxxxxyz_0[i] = g_yyzzz_0_xxxxxxyz_0[i] * fbe_0 - g_yyzzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xyyzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxxxxzz_0[i] = g_xxzzz_0_xxxxxxzz_0[i] * fbe_0 - g_xxzzz_0_xxxxxxzz_1[i] * fz_be_0 + g_xxyzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxxxyyy_0[i] = 2.0 * g_xxyyz_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxxxyyy_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxxxyyz_0[i] = g_yyzzz_0_xxxxxyyz_0[i] * fbe_0 - g_yyzzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xyyzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxxxyzz_0[i] = g_yyzzz_0_xxxxxyzz_0[i] * fbe_0 - g_yyzzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xyyzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxxxzzz_0[i] = g_xxzzz_0_xxxxxzzz_0[i] * fbe_0 - g_xxzzz_0_xxxxxzzz_1[i] * fz_be_0 + g_xxyzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxxyyyy_0[i] = 2.0 * g_xxyyz_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxxyyyy_1[i] * fz_be_0 + g_xxyyzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxxyyyz_0[i] = g_yyzzz_0_xxxxyyyz_0[i] * fbe_0 - g_yyzzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xyyzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxxyyzz_0[i] = g_yyzzz_0_xxxxyyzz_0[i] * fbe_0 - g_yyzzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xyyzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxxyzzz_0[i] = g_yyzzz_0_xxxxyzzz_0[i] * fbe_0 - g_yyzzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xyyzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxxzzzz_0[i] = g_xxzzz_0_xxxxzzzz_0[i] * fbe_0 - g_xxzzz_0_xxxxzzzz_1[i] * fz_be_0 + g_xxyzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxyyyyy_0[i] = 2.0 * g_xxyyz_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxyyyyy_1[i] * fz_be_0 + g_xxyyzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxyyyyz_0[i] = g_yyzzz_0_xxxyyyyz_0[i] * fbe_0 - g_yyzzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxyyyzz_0[i] = g_yyzzz_0_xxxyyyzz_0[i] * fbe_0 - g_yyzzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxyyzzz_0[i] = g_yyzzz_0_xxxyyzzz_0[i] * fbe_0 - g_yyzzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxyzzzz_0[i] = g_yyzzz_0_xxxyzzzz_0[i] * fbe_0 - g_yyzzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxzzzzz_0[i] = g_xxzzz_0_xxxzzzzz_0[i] * fbe_0 - g_xxzzz_0_xxxzzzzz_1[i] * fz_be_0 + g_xxyzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxyyyyyy_0[i] = 2.0 * g_xxyyz_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxyyyyyy_1[i] * fz_be_0 + g_xxyyzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxyyyyyz_0[i] = g_yyzzz_0_xxyyyyyz_0[i] * fbe_0 - g_yyzzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxyyyyzz_0[i] = g_yyzzz_0_xxyyyyzz_0[i] * fbe_0 - g_yyzzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxyyyzzz_0[i] = g_yyzzz_0_xxyyyzzz_0[i] * fbe_0 - g_yyzzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxyyzzzz_0[i] = g_yyzzz_0_xxyyzzzz_0[i] * fbe_0 - g_yyzzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxyzzzzz_0[i] = g_yyzzz_0_xxyzzzzz_0[i] * fbe_0 - g_yyzzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxzzzzzz_0[i] = g_xxzzz_0_xxzzzzzz_0[i] * fbe_0 - g_xxzzz_0_xxzzzzzz_1[i] * fz_be_0 + g_xxyzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xyyyyyyy_0[i] = 2.0 * g_xxyyz_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xyyyyyyy_1[i] * fz_be_0 + g_xxyyzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xyyyyyyz_0[i] = g_yyzzz_0_xyyyyyyz_0[i] * fbe_0 - g_yyzzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyyyyyzz_0[i] = g_yyzzz_0_xyyyyyzz_0[i] * fbe_0 - g_yyzzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyyyyzzz_0[i] = g_yyzzz_0_xyyyyzzz_0[i] * fbe_0 - g_yyzzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyyyzzzz_0[i] = g_yyzzz_0_xyyyzzzz_0[i] * fbe_0 - g_yyzzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyyzzzzz_0[i] = g_yyzzz_0_xyyzzzzz_0[i] * fbe_0 - g_yyzzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyzzzzzz_0[i] = g_yyzzz_0_xyzzzzzz_0[i] * fbe_0 - g_yyzzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xzzzzzzz_0[i] = g_xxzzz_0_xzzzzzzz_0[i] * fbe_0 - g_xxzzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xxyzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_yyyyyyyy_0[i] = g_yyzzz_0_yyyyyyyy_0[i] * fbe_0 - g_yyzzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xyyzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyyyyyz_0[i] = g_yyzzz_0_yyyyyyyz_0[i] * fbe_0 - g_yyzzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyyyyzz_0[i] = g_yyzzz_0_yyyyyyzz_0[i] * fbe_0 - g_yyzzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyyyzzz_0[i] = g_yyzzz_0_yyyyyzzz_0[i] * fbe_0 - g_yyzzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyyzzzz_0[i] = g_yyzzz_0_yyyyzzzz_0[i] * fbe_0 - g_yyzzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyzzzzz_0[i] = g_yyzzz_0_yyyzzzzz_0[i] * fbe_0 - g_yyzzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyzzzzzz_0[i] = g_yyzzz_0_yyzzzzzz_0[i] * fbe_0 - g_yyzzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yzzzzzzz_0[i] = g_yyzzz_0_yzzzzzzz_0[i] * fbe_0 - g_yyzzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_zzzzzzzz_0[i] = g_yyzzz_0_zzzzzzzz_0[i] * fbe_0 - g_yyzzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xyyzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 855-900 components of targeted buffer : KSL

    auto g_xxyzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 855);

    auto g_xxyzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 856);

    auto g_xxyzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 857);

    auto g_xxyzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 858);

    auto g_xxyzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 859);

    auto g_xxyzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 860);

    auto g_xxyzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 861);

    auto g_xxyzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 862);

    auto g_xxyzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 863);

    auto g_xxyzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 864);

    auto g_xxyzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 865);

    auto g_xxyzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 866);

    auto g_xxyzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 867);

    auto g_xxyzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 868);

    auto g_xxyzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 869);

    auto g_xxyzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 870);

    auto g_xxyzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 871);

    auto g_xxyzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 872);

    auto g_xxyzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 873);

    auto g_xxyzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 874);

    auto g_xxyzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 875);

    auto g_xxyzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 876);

    auto g_xxyzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 877);

    auto g_xxyzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 878);

    auto g_xxyzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 879);

    auto g_xxyzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 880);

    auto g_xxyzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 881);

    auto g_xxyzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 882);

    auto g_xxyzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 883);

    auto g_xxyzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 884);

    auto g_xxyzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 885);

    auto g_xxyzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 886);

    auto g_xxyzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 887);

    auto g_xxyzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 888);

    auto g_xxyzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 889);

    auto g_xxyzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 890);

    auto g_xxyzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 891);

    auto g_xxyzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 892);

    auto g_xxyzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 893);

    auto g_xxyzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 894);

    auto g_xxyzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 895);

    auto g_xxyzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 896);

    auto g_xxyzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 897);

    auto g_xxyzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 898);

    auto g_xxyzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 899);

    #pragma omp simd aligned(g_xxyzzzz_0_xxxxxxxx_0, g_xxyzzzz_0_xxxxxxxy_0, g_xxyzzzz_0_xxxxxxxz_0, g_xxyzzzz_0_xxxxxxyy_0, g_xxyzzzz_0_xxxxxxyz_0, g_xxyzzzz_0_xxxxxxzz_0, g_xxyzzzz_0_xxxxxyyy_0, g_xxyzzzz_0_xxxxxyyz_0, g_xxyzzzz_0_xxxxxyzz_0, g_xxyzzzz_0_xxxxxzzz_0, g_xxyzzzz_0_xxxxyyyy_0, g_xxyzzzz_0_xxxxyyyz_0, g_xxyzzzz_0_xxxxyyzz_0, g_xxyzzzz_0_xxxxyzzz_0, g_xxyzzzz_0_xxxxzzzz_0, g_xxyzzzz_0_xxxyyyyy_0, g_xxyzzzz_0_xxxyyyyz_0, g_xxyzzzz_0_xxxyyyzz_0, g_xxyzzzz_0_xxxyyzzz_0, g_xxyzzzz_0_xxxyzzzz_0, g_xxyzzzz_0_xxxzzzzz_0, g_xxyzzzz_0_xxyyyyyy_0, g_xxyzzzz_0_xxyyyyyz_0, g_xxyzzzz_0_xxyyyyzz_0, g_xxyzzzz_0_xxyyyzzz_0, g_xxyzzzz_0_xxyyzzzz_0, g_xxyzzzz_0_xxyzzzzz_0, g_xxyzzzz_0_xxzzzzzz_0, g_xxyzzzz_0_xyyyyyyy_0, g_xxyzzzz_0_xyyyyyyz_0, g_xxyzzzz_0_xyyyyyzz_0, g_xxyzzzz_0_xyyyyzzz_0, g_xxyzzzz_0_xyyyzzzz_0, g_xxyzzzz_0_xyyzzzzz_0, g_xxyzzzz_0_xyzzzzzz_0, g_xxyzzzz_0_xzzzzzzz_0, g_xxyzzzz_0_yyyyyyyy_0, g_xxyzzzz_0_yyyyyyyz_0, g_xxyzzzz_0_yyyyyyzz_0, g_xxyzzzz_0_yyyyyzzz_0, g_xxyzzzz_0_yyyyzzzz_0, g_xxyzzzz_0_yyyzzzzz_0, g_xxyzzzz_0_yyzzzzzz_0, g_xxyzzzz_0_yzzzzzzz_0, g_xxyzzzz_0_zzzzzzzz_0, g_xxzzzz_0_xxxxxxx_1, g_xxzzzz_0_xxxxxxxx_1, g_xxzzzz_0_xxxxxxxy_1, g_xxzzzz_0_xxxxxxxz_1, g_xxzzzz_0_xxxxxxy_1, g_xxzzzz_0_xxxxxxyy_1, g_xxzzzz_0_xxxxxxyz_1, g_xxzzzz_0_xxxxxxz_1, g_xxzzzz_0_xxxxxxzz_1, g_xxzzzz_0_xxxxxyy_1, g_xxzzzz_0_xxxxxyyy_1, g_xxzzzz_0_xxxxxyyz_1, g_xxzzzz_0_xxxxxyz_1, g_xxzzzz_0_xxxxxyzz_1, g_xxzzzz_0_xxxxxzz_1, g_xxzzzz_0_xxxxxzzz_1, g_xxzzzz_0_xxxxyyy_1, g_xxzzzz_0_xxxxyyyy_1, g_xxzzzz_0_xxxxyyyz_1, g_xxzzzz_0_xxxxyyz_1, g_xxzzzz_0_xxxxyyzz_1, g_xxzzzz_0_xxxxyzz_1, g_xxzzzz_0_xxxxyzzz_1, g_xxzzzz_0_xxxxzzz_1, g_xxzzzz_0_xxxxzzzz_1, g_xxzzzz_0_xxxyyyy_1, g_xxzzzz_0_xxxyyyyy_1, g_xxzzzz_0_xxxyyyyz_1, g_xxzzzz_0_xxxyyyz_1, g_xxzzzz_0_xxxyyyzz_1, g_xxzzzz_0_xxxyyzz_1, g_xxzzzz_0_xxxyyzzz_1, g_xxzzzz_0_xxxyzzz_1, g_xxzzzz_0_xxxyzzzz_1, g_xxzzzz_0_xxxzzzz_1, g_xxzzzz_0_xxxzzzzz_1, g_xxzzzz_0_xxyyyyy_1, g_xxzzzz_0_xxyyyyyy_1, g_xxzzzz_0_xxyyyyyz_1, g_xxzzzz_0_xxyyyyz_1, g_xxzzzz_0_xxyyyyzz_1, g_xxzzzz_0_xxyyyzz_1, g_xxzzzz_0_xxyyyzzz_1, g_xxzzzz_0_xxyyzzz_1, g_xxzzzz_0_xxyyzzzz_1, g_xxzzzz_0_xxyzzzz_1, g_xxzzzz_0_xxyzzzzz_1, g_xxzzzz_0_xxzzzzz_1, g_xxzzzz_0_xxzzzzzz_1, g_xxzzzz_0_xyyyyyy_1, g_xxzzzz_0_xyyyyyyy_1, g_xxzzzz_0_xyyyyyyz_1, g_xxzzzz_0_xyyyyyz_1, g_xxzzzz_0_xyyyyyzz_1, g_xxzzzz_0_xyyyyzz_1, g_xxzzzz_0_xyyyyzzz_1, g_xxzzzz_0_xyyyzzz_1, g_xxzzzz_0_xyyyzzzz_1, g_xxzzzz_0_xyyzzzz_1, g_xxzzzz_0_xyyzzzzz_1, g_xxzzzz_0_xyzzzzz_1, g_xxzzzz_0_xyzzzzzz_1, g_xxzzzz_0_xzzzzzz_1, g_xxzzzz_0_xzzzzzzz_1, g_xxzzzz_0_yyyyyyy_1, g_xxzzzz_0_yyyyyyyy_1, g_xxzzzz_0_yyyyyyyz_1, g_xxzzzz_0_yyyyyyz_1, g_xxzzzz_0_yyyyyyzz_1, g_xxzzzz_0_yyyyyzz_1, g_xxzzzz_0_yyyyyzzz_1, g_xxzzzz_0_yyyyzzz_1, g_xxzzzz_0_yyyyzzzz_1, g_xxzzzz_0_yyyzzzz_1, g_xxzzzz_0_yyyzzzzz_1, g_xxzzzz_0_yyzzzzz_1, g_xxzzzz_0_yyzzzzzz_1, g_xxzzzz_0_yzzzzzz_1, g_xxzzzz_0_yzzzzzzz_1, g_xxzzzz_0_zzzzzzz_1, g_xxzzzz_0_zzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzzz_0_xxxxxxxx_0[i] = g_xxzzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxxxy_0[i] = g_xxzzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxxxy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxxxz_0[i] = g_xxzzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxxyy_0[i] = 2.0 * g_xxzzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxxyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxxyz_0[i] = g_xxzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxxzz_0[i] = g_xxzzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxyyy_0[i] = 3.0 * g_xxzzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxyyz_0[i] = 2.0 * g_xxzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxyzz_0[i] = g_xxzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxzzz_0[i] = g_xxzzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxyyyy_0[i] = 4.0 * g_xxzzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxyyyz_0[i] = 3.0 * g_xxzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxyyzz_0[i] = 2.0 * g_xxzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxyzzz_0[i] = g_xxzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxzzzz_0[i] = g_xxzzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyyyyy_0[i] = 5.0 * g_xxzzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyyyyz_0[i] = 4.0 * g_xxzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyyyzz_0[i] = 3.0 * g_xxzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyyzzz_0[i] = 2.0 * g_xxzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyzzzz_0[i] = g_xxzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxzzzzz_0[i] = g_xxzzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyyyyy_0[i] = 6.0 * g_xxzzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyyyyz_0[i] = 5.0 * g_xxzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyyyzz_0[i] = 4.0 * g_xxzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyyzzz_0[i] = 3.0 * g_xxzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyzzzz_0[i] = 2.0 * g_xxzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyzzzzz_0[i] = g_xxzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxzzzzzz_0[i] = g_xxzzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyyyyy_0[i] = 7.0 * g_xxzzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyyyyz_0[i] = 6.0 * g_xxzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyyyzz_0[i] = 5.0 * g_xxzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyyzzz_0[i] = 4.0 * g_xxzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyzzzz_0[i] = 3.0 * g_xxzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyzzzzz_0[i] = 2.0 * g_xxzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyzzzzzz_0[i] = g_xxzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xzzzzzzz_0[i] = g_xxzzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyyyyy_0[i] = 8.0 * g_xxzzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyyyyz_0[i] = 7.0 * g_xxzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyyyzz_0[i] = 6.0 * g_xxzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyyzzz_0[i] = 5.0 * g_xxzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyzzzz_0[i] = 4.0 * g_xxzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyzzzzz_0[i] = 3.0 * g_xxzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyzzzzzz_0[i] = 2.0 * g_xxzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yzzzzzzz_0[i] = g_xxzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_zzzzzzzz_0[i] = g_xxzzzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 900-945 components of targeted buffer : KSL

    auto g_xxzzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 900);

    auto g_xxzzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 901);

    auto g_xxzzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 902);

    auto g_xxzzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 903);

    auto g_xxzzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 904);

    auto g_xxzzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 905);

    auto g_xxzzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 906);

    auto g_xxzzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 907);

    auto g_xxzzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 908);

    auto g_xxzzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 909);

    auto g_xxzzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 910);

    auto g_xxzzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 911);

    auto g_xxzzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 912);

    auto g_xxzzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 913);

    auto g_xxzzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 914);

    auto g_xxzzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 915);

    auto g_xxzzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 916);

    auto g_xxzzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 917);

    auto g_xxzzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 918);

    auto g_xxzzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 919);

    auto g_xxzzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 920);

    auto g_xxzzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 921);

    auto g_xxzzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 922);

    auto g_xxzzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 923);

    auto g_xxzzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 924);

    auto g_xxzzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 925);

    auto g_xxzzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 926);

    auto g_xxzzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 927);

    auto g_xxzzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 928);

    auto g_xxzzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 929);

    auto g_xxzzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 930);

    auto g_xxzzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 931);

    auto g_xxzzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 932);

    auto g_xxzzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 933);

    auto g_xxzzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 934);

    auto g_xxzzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 935);

    auto g_xxzzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 936);

    auto g_xxzzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 937);

    auto g_xxzzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 938);

    auto g_xxzzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 939);

    auto g_xxzzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 940);

    auto g_xxzzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 941);

    auto g_xxzzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 942);

    auto g_xxzzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 943);

    auto g_xxzzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 944);

    #pragma omp simd aligned(g_xxzzz_0_xxxxxxxx_0, g_xxzzz_0_xxxxxxxx_1, g_xxzzz_0_xxxxxxxy_0, g_xxzzz_0_xxxxxxxy_1, g_xxzzz_0_xxxxxxyy_0, g_xxzzz_0_xxxxxxyy_1, g_xxzzz_0_xxxxxyyy_0, g_xxzzz_0_xxxxxyyy_1, g_xxzzz_0_xxxxyyyy_0, g_xxzzz_0_xxxxyyyy_1, g_xxzzz_0_xxxyyyyy_0, g_xxzzz_0_xxxyyyyy_1, g_xxzzz_0_xxyyyyyy_0, g_xxzzz_0_xxyyyyyy_1, g_xxzzz_0_xyyyyyyy_0, g_xxzzz_0_xyyyyyyy_1, g_xxzzzz_0_xxxxxxxx_1, g_xxzzzz_0_xxxxxxxy_1, g_xxzzzz_0_xxxxxxyy_1, g_xxzzzz_0_xxxxxyyy_1, g_xxzzzz_0_xxxxyyyy_1, g_xxzzzz_0_xxxyyyyy_1, g_xxzzzz_0_xxyyyyyy_1, g_xxzzzz_0_xyyyyyyy_1, g_xxzzzzz_0_xxxxxxxx_0, g_xxzzzzz_0_xxxxxxxy_0, g_xxzzzzz_0_xxxxxxxz_0, g_xxzzzzz_0_xxxxxxyy_0, g_xxzzzzz_0_xxxxxxyz_0, g_xxzzzzz_0_xxxxxxzz_0, g_xxzzzzz_0_xxxxxyyy_0, g_xxzzzzz_0_xxxxxyyz_0, g_xxzzzzz_0_xxxxxyzz_0, g_xxzzzzz_0_xxxxxzzz_0, g_xxzzzzz_0_xxxxyyyy_0, g_xxzzzzz_0_xxxxyyyz_0, g_xxzzzzz_0_xxxxyyzz_0, g_xxzzzzz_0_xxxxyzzz_0, g_xxzzzzz_0_xxxxzzzz_0, g_xxzzzzz_0_xxxyyyyy_0, g_xxzzzzz_0_xxxyyyyz_0, g_xxzzzzz_0_xxxyyyzz_0, g_xxzzzzz_0_xxxyyzzz_0, g_xxzzzzz_0_xxxyzzzz_0, g_xxzzzzz_0_xxxzzzzz_0, g_xxzzzzz_0_xxyyyyyy_0, g_xxzzzzz_0_xxyyyyyz_0, g_xxzzzzz_0_xxyyyyzz_0, g_xxzzzzz_0_xxyyyzzz_0, g_xxzzzzz_0_xxyyzzzz_0, g_xxzzzzz_0_xxyzzzzz_0, g_xxzzzzz_0_xxzzzzzz_0, g_xxzzzzz_0_xyyyyyyy_0, g_xxzzzzz_0_xyyyyyyz_0, g_xxzzzzz_0_xyyyyyzz_0, g_xxzzzzz_0_xyyyyzzz_0, g_xxzzzzz_0_xyyyzzzz_0, g_xxzzzzz_0_xyyzzzzz_0, g_xxzzzzz_0_xyzzzzzz_0, g_xxzzzzz_0_xzzzzzzz_0, g_xxzzzzz_0_yyyyyyyy_0, g_xxzzzzz_0_yyyyyyyz_0, g_xxzzzzz_0_yyyyyyzz_0, g_xxzzzzz_0_yyyyyzzz_0, g_xxzzzzz_0_yyyyzzzz_0, g_xxzzzzz_0_yyyzzzzz_0, g_xxzzzzz_0_yyzzzzzz_0, g_xxzzzzz_0_yzzzzzzz_0, g_xxzzzzz_0_zzzzzzzz_0, g_xzzzzz_0_xxxxxxxz_1, g_xzzzzz_0_xxxxxxyz_1, g_xzzzzz_0_xxxxxxz_1, g_xzzzzz_0_xxxxxxzz_1, g_xzzzzz_0_xxxxxyyz_1, g_xzzzzz_0_xxxxxyz_1, g_xzzzzz_0_xxxxxyzz_1, g_xzzzzz_0_xxxxxzz_1, g_xzzzzz_0_xxxxxzzz_1, g_xzzzzz_0_xxxxyyyz_1, g_xzzzzz_0_xxxxyyz_1, g_xzzzzz_0_xxxxyyzz_1, g_xzzzzz_0_xxxxyzz_1, g_xzzzzz_0_xxxxyzzz_1, g_xzzzzz_0_xxxxzzz_1, g_xzzzzz_0_xxxxzzzz_1, g_xzzzzz_0_xxxyyyyz_1, g_xzzzzz_0_xxxyyyz_1, g_xzzzzz_0_xxxyyyzz_1, g_xzzzzz_0_xxxyyzz_1, g_xzzzzz_0_xxxyyzzz_1, g_xzzzzz_0_xxxyzzz_1, g_xzzzzz_0_xxxyzzzz_1, g_xzzzzz_0_xxxzzzz_1, g_xzzzzz_0_xxxzzzzz_1, g_xzzzzz_0_xxyyyyyz_1, g_xzzzzz_0_xxyyyyz_1, g_xzzzzz_0_xxyyyyzz_1, g_xzzzzz_0_xxyyyzz_1, g_xzzzzz_0_xxyyyzzz_1, g_xzzzzz_0_xxyyzzz_1, g_xzzzzz_0_xxyyzzzz_1, g_xzzzzz_0_xxyzzzz_1, g_xzzzzz_0_xxyzzzzz_1, g_xzzzzz_0_xxzzzzz_1, g_xzzzzz_0_xxzzzzzz_1, g_xzzzzz_0_xyyyyyyz_1, g_xzzzzz_0_xyyyyyz_1, g_xzzzzz_0_xyyyyyzz_1, g_xzzzzz_0_xyyyyzz_1, g_xzzzzz_0_xyyyyzzz_1, g_xzzzzz_0_xyyyzzz_1, g_xzzzzz_0_xyyyzzzz_1, g_xzzzzz_0_xyyzzzz_1, g_xzzzzz_0_xyyzzzzz_1, g_xzzzzz_0_xyzzzzz_1, g_xzzzzz_0_xyzzzzzz_1, g_xzzzzz_0_xzzzzzz_1, g_xzzzzz_0_xzzzzzzz_1, g_xzzzzz_0_yyyyyyyy_1, g_xzzzzz_0_yyyyyyyz_1, g_xzzzzz_0_yyyyyyz_1, g_xzzzzz_0_yyyyyyzz_1, g_xzzzzz_0_yyyyyzz_1, g_xzzzzz_0_yyyyyzzz_1, g_xzzzzz_0_yyyyzzz_1, g_xzzzzz_0_yyyyzzzz_1, g_xzzzzz_0_yyyzzzz_1, g_xzzzzz_0_yyyzzzzz_1, g_xzzzzz_0_yyzzzzz_1, g_xzzzzz_0_yyzzzzzz_1, g_xzzzzz_0_yzzzzzz_1, g_xzzzzz_0_yzzzzzzz_1, g_xzzzzz_0_zzzzzzz_1, g_xzzzzz_0_zzzzzzzz_1, g_zzzzz_0_xxxxxxxz_0, g_zzzzz_0_xxxxxxxz_1, g_zzzzz_0_xxxxxxyz_0, g_zzzzz_0_xxxxxxyz_1, g_zzzzz_0_xxxxxxzz_0, g_zzzzz_0_xxxxxxzz_1, g_zzzzz_0_xxxxxyyz_0, g_zzzzz_0_xxxxxyyz_1, g_zzzzz_0_xxxxxyzz_0, g_zzzzz_0_xxxxxyzz_1, g_zzzzz_0_xxxxxzzz_0, g_zzzzz_0_xxxxxzzz_1, g_zzzzz_0_xxxxyyyz_0, g_zzzzz_0_xxxxyyyz_1, g_zzzzz_0_xxxxyyzz_0, g_zzzzz_0_xxxxyyzz_1, g_zzzzz_0_xxxxyzzz_0, g_zzzzz_0_xxxxyzzz_1, g_zzzzz_0_xxxxzzzz_0, g_zzzzz_0_xxxxzzzz_1, g_zzzzz_0_xxxyyyyz_0, g_zzzzz_0_xxxyyyyz_1, g_zzzzz_0_xxxyyyzz_0, g_zzzzz_0_xxxyyyzz_1, g_zzzzz_0_xxxyyzzz_0, g_zzzzz_0_xxxyyzzz_1, g_zzzzz_0_xxxyzzzz_0, g_zzzzz_0_xxxyzzzz_1, g_zzzzz_0_xxxzzzzz_0, g_zzzzz_0_xxxzzzzz_1, g_zzzzz_0_xxyyyyyz_0, g_zzzzz_0_xxyyyyyz_1, g_zzzzz_0_xxyyyyzz_0, g_zzzzz_0_xxyyyyzz_1, g_zzzzz_0_xxyyyzzz_0, g_zzzzz_0_xxyyyzzz_1, g_zzzzz_0_xxyyzzzz_0, g_zzzzz_0_xxyyzzzz_1, g_zzzzz_0_xxyzzzzz_0, g_zzzzz_0_xxyzzzzz_1, g_zzzzz_0_xxzzzzzz_0, g_zzzzz_0_xxzzzzzz_1, g_zzzzz_0_xyyyyyyz_0, g_zzzzz_0_xyyyyyyz_1, g_zzzzz_0_xyyyyyzz_0, g_zzzzz_0_xyyyyyzz_1, g_zzzzz_0_xyyyyzzz_0, g_zzzzz_0_xyyyyzzz_1, g_zzzzz_0_xyyyzzzz_0, g_zzzzz_0_xyyyzzzz_1, g_zzzzz_0_xyyzzzzz_0, g_zzzzz_0_xyyzzzzz_1, g_zzzzz_0_xyzzzzzz_0, g_zzzzz_0_xyzzzzzz_1, g_zzzzz_0_xzzzzzzz_0, g_zzzzz_0_xzzzzzzz_1, g_zzzzz_0_yyyyyyyy_0, g_zzzzz_0_yyyyyyyy_1, g_zzzzz_0_yyyyyyyz_0, g_zzzzz_0_yyyyyyyz_1, g_zzzzz_0_yyyyyyzz_0, g_zzzzz_0_yyyyyyzz_1, g_zzzzz_0_yyyyyzzz_0, g_zzzzz_0_yyyyyzzz_1, g_zzzzz_0_yyyyzzzz_0, g_zzzzz_0_yyyyzzzz_1, g_zzzzz_0_yyyzzzzz_0, g_zzzzz_0_yyyzzzzz_1, g_zzzzz_0_yyzzzzzz_0, g_zzzzz_0_yyzzzzzz_1, g_zzzzz_0_yzzzzzzz_0, g_zzzzz_0_yzzzzzzz_1, g_zzzzz_0_zzzzzzzz_0, g_zzzzz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzzz_0_xxxxxxxx_0[i] = 4.0 * g_xxzzz_0_xxxxxxxx_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxxxxx_1[i] * fz_be_0 + g_xxzzzz_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxxxxy_0[i] = 4.0 * g_xxzzz_0_xxxxxxxy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxxxxy_1[i] * fz_be_0 + g_xxzzzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxxxxz_0[i] = g_zzzzz_0_xxxxxxxz_0[i] * fbe_0 - g_zzzzz_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xzzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxxxyy_0[i] = 4.0 * g_xxzzz_0_xxxxxxyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxxxyy_1[i] * fz_be_0 + g_xxzzzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxxxyz_0[i] = g_zzzzz_0_xxxxxxyz_0[i] * fbe_0 - g_zzzzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xzzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxxxzz_0[i] = g_zzzzz_0_xxxxxxzz_0[i] * fbe_0 - g_zzzzz_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xzzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxxyyy_0[i] = 4.0 * g_xxzzz_0_xxxxxyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxxyyy_1[i] * fz_be_0 + g_xxzzzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxxyyz_0[i] = g_zzzzz_0_xxxxxyyz_0[i] * fbe_0 - g_zzzzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xzzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxxyzz_0[i] = g_zzzzz_0_xxxxxyzz_0[i] * fbe_0 - g_zzzzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xzzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxxzzz_0[i] = g_zzzzz_0_xxxxxzzz_0[i] * fbe_0 - g_zzzzz_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xzzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxyyyy_0[i] = 4.0 * g_xxzzz_0_xxxxyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxyyyy_1[i] * fz_be_0 + g_xxzzzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxyyyz_0[i] = g_zzzzz_0_xxxxyyyz_0[i] * fbe_0 - g_zzzzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxyyzz_0[i] = g_zzzzz_0_xxxxyyzz_0[i] * fbe_0 - g_zzzzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxyzzz_0[i] = g_zzzzz_0_xxxxyzzz_0[i] * fbe_0 - g_zzzzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxzzzz_0[i] = g_zzzzz_0_xxxxzzzz_0[i] * fbe_0 - g_zzzzz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxyyyyy_0[i] = 4.0 * g_xxzzz_0_xxxyyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxyyyyy_1[i] * fz_be_0 + g_xxzzzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxyyyyz_0[i] = g_zzzzz_0_xxxyyyyz_0[i] * fbe_0 - g_zzzzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxyyyzz_0[i] = g_zzzzz_0_xxxyyyzz_0[i] * fbe_0 - g_zzzzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxyyzzz_0[i] = g_zzzzz_0_xxxyyzzz_0[i] * fbe_0 - g_zzzzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxyzzzz_0[i] = g_zzzzz_0_xxxyzzzz_0[i] * fbe_0 - g_zzzzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxzzzzz_0[i] = g_zzzzz_0_xxxzzzzz_0[i] * fbe_0 - g_zzzzz_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyyyyyy_0[i] = 4.0 * g_xxzzz_0_xxyyyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxyyyyyy_1[i] * fz_be_0 + g_xxzzzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxyyyyyz_0[i] = g_zzzzz_0_xxyyyyyz_0[i] * fbe_0 - g_zzzzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyyyyzz_0[i] = g_zzzzz_0_xxyyyyzz_0[i] * fbe_0 - g_zzzzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyyyzzz_0[i] = g_zzzzz_0_xxyyyzzz_0[i] * fbe_0 - g_zzzzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyyzzzz_0[i] = g_zzzzz_0_xxyyzzzz_0[i] * fbe_0 - g_zzzzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyzzzzz_0[i] = g_zzzzz_0_xxyzzzzz_0[i] * fbe_0 - g_zzzzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxzzzzzz_0[i] = g_zzzzz_0_xxzzzzzz_0[i] * fbe_0 - g_zzzzz_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyyyyyy_0[i] = 4.0 * g_xxzzz_0_xyyyyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xyyyyyyy_1[i] * fz_be_0 + g_xxzzzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xyyyyyyz_0[i] = g_zzzzz_0_xyyyyyyz_0[i] * fbe_0 - g_zzzzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyyyyzz_0[i] = g_zzzzz_0_xyyyyyzz_0[i] * fbe_0 - g_zzzzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyyyzzz_0[i] = g_zzzzz_0_xyyyyzzz_0[i] * fbe_0 - g_zzzzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyyzzzz_0[i] = g_zzzzz_0_xyyyzzzz_0[i] * fbe_0 - g_zzzzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyzzzzz_0[i] = g_zzzzz_0_xyyzzzzz_0[i] * fbe_0 - g_zzzzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyzzzzzz_0[i] = g_zzzzz_0_xyzzzzzz_0[i] * fbe_0 - g_zzzzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xzzzzzzz_0[i] = g_zzzzz_0_xzzzzzzz_0[i] * fbe_0 - g_zzzzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyyyyy_0[i] = g_zzzzz_0_yyyyyyyy_0[i] * fbe_0 - g_zzzzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xzzzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyyyyz_0[i] = g_zzzzz_0_yyyyyyyz_0[i] * fbe_0 - g_zzzzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyyyzz_0[i] = g_zzzzz_0_yyyyyyzz_0[i] * fbe_0 - g_zzzzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyyzzz_0[i] = g_zzzzz_0_yyyyyzzz_0[i] * fbe_0 - g_zzzzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyzzzz_0[i] = g_zzzzz_0_yyyyzzzz_0[i] * fbe_0 - g_zzzzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyzzzzz_0[i] = g_zzzzz_0_yyyzzzzz_0[i] * fbe_0 - g_zzzzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyzzzzzz_0[i] = g_zzzzz_0_yyzzzzzz_0[i] * fbe_0 - g_zzzzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yzzzzzzz_0[i] = g_zzzzz_0_yzzzzzzz_0[i] * fbe_0 - g_zzzzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_zzzzzzzz_0[i] = g_zzzzz_0_zzzzzzzz_0[i] * fbe_0 - g_zzzzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 945-990 components of targeted buffer : KSL

    auto g_xyyyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 945);

    auto g_xyyyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 946);

    auto g_xyyyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 947);

    auto g_xyyyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 948);

    auto g_xyyyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 949);

    auto g_xyyyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 950);

    auto g_xyyyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 951);

    auto g_xyyyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 952);

    auto g_xyyyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 953);

    auto g_xyyyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 954);

    auto g_xyyyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 955);

    auto g_xyyyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 956);

    auto g_xyyyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 957);

    auto g_xyyyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 958);

    auto g_xyyyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 959);

    auto g_xyyyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 960);

    auto g_xyyyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 961);

    auto g_xyyyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 962);

    auto g_xyyyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 963);

    auto g_xyyyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 964);

    auto g_xyyyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 965);

    auto g_xyyyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 966);

    auto g_xyyyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 967);

    auto g_xyyyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 968);

    auto g_xyyyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 969);

    auto g_xyyyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 970);

    auto g_xyyyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 971);

    auto g_xyyyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 972);

    auto g_xyyyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 973);

    auto g_xyyyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 974);

    auto g_xyyyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 975);

    auto g_xyyyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 976);

    auto g_xyyyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 977);

    auto g_xyyyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 978);

    auto g_xyyyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 979);

    auto g_xyyyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 980);

    auto g_xyyyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 981);

    auto g_xyyyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 982);

    auto g_xyyyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 983);

    auto g_xyyyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 984);

    auto g_xyyyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 985);

    auto g_xyyyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 986);

    auto g_xyyyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 987);

    auto g_xyyyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 988);

    auto g_xyyyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 989);

    #pragma omp simd aligned(g_xyyyyyy_0_xxxxxxxx_0, g_xyyyyyy_0_xxxxxxxy_0, g_xyyyyyy_0_xxxxxxxz_0, g_xyyyyyy_0_xxxxxxyy_0, g_xyyyyyy_0_xxxxxxyz_0, g_xyyyyyy_0_xxxxxxzz_0, g_xyyyyyy_0_xxxxxyyy_0, g_xyyyyyy_0_xxxxxyyz_0, g_xyyyyyy_0_xxxxxyzz_0, g_xyyyyyy_0_xxxxxzzz_0, g_xyyyyyy_0_xxxxyyyy_0, g_xyyyyyy_0_xxxxyyyz_0, g_xyyyyyy_0_xxxxyyzz_0, g_xyyyyyy_0_xxxxyzzz_0, g_xyyyyyy_0_xxxxzzzz_0, g_xyyyyyy_0_xxxyyyyy_0, g_xyyyyyy_0_xxxyyyyz_0, g_xyyyyyy_0_xxxyyyzz_0, g_xyyyyyy_0_xxxyyzzz_0, g_xyyyyyy_0_xxxyzzzz_0, g_xyyyyyy_0_xxxzzzzz_0, g_xyyyyyy_0_xxyyyyyy_0, g_xyyyyyy_0_xxyyyyyz_0, g_xyyyyyy_0_xxyyyyzz_0, g_xyyyyyy_0_xxyyyzzz_0, g_xyyyyyy_0_xxyyzzzz_0, g_xyyyyyy_0_xxyzzzzz_0, g_xyyyyyy_0_xxzzzzzz_0, g_xyyyyyy_0_xyyyyyyy_0, g_xyyyyyy_0_xyyyyyyz_0, g_xyyyyyy_0_xyyyyyzz_0, g_xyyyyyy_0_xyyyyzzz_0, g_xyyyyyy_0_xyyyzzzz_0, g_xyyyyyy_0_xyyzzzzz_0, g_xyyyyyy_0_xyzzzzzz_0, g_xyyyyyy_0_xzzzzzzz_0, g_xyyyyyy_0_yyyyyyyy_0, g_xyyyyyy_0_yyyyyyyz_0, g_xyyyyyy_0_yyyyyyzz_0, g_xyyyyyy_0_yyyyyzzz_0, g_xyyyyyy_0_yyyyzzzz_0, g_xyyyyyy_0_yyyzzzzz_0, g_xyyyyyy_0_yyzzzzzz_0, g_xyyyyyy_0_yzzzzzzz_0, g_xyyyyyy_0_zzzzzzzz_0, g_yyyyyy_0_xxxxxxx_1, g_yyyyyy_0_xxxxxxxx_1, g_yyyyyy_0_xxxxxxxy_1, g_yyyyyy_0_xxxxxxxz_1, g_yyyyyy_0_xxxxxxy_1, g_yyyyyy_0_xxxxxxyy_1, g_yyyyyy_0_xxxxxxyz_1, g_yyyyyy_0_xxxxxxz_1, g_yyyyyy_0_xxxxxxzz_1, g_yyyyyy_0_xxxxxyy_1, g_yyyyyy_0_xxxxxyyy_1, g_yyyyyy_0_xxxxxyyz_1, g_yyyyyy_0_xxxxxyz_1, g_yyyyyy_0_xxxxxyzz_1, g_yyyyyy_0_xxxxxzz_1, g_yyyyyy_0_xxxxxzzz_1, g_yyyyyy_0_xxxxyyy_1, g_yyyyyy_0_xxxxyyyy_1, g_yyyyyy_0_xxxxyyyz_1, g_yyyyyy_0_xxxxyyz_1, g_yyyyyy_0_xxxxyyzz_1, g_yyyyyy_0_xxxxyzz_1, g_yyyyyy_0_xxxxyzzz_1, g_yyyyyy_0_xxxxzzz_1, g_yyyyyy_0_xxxxzzzz_1, g_yyyyyy_0_xxxyyyy_1, g_yyyyyy_0_xxxyyyyy_1, g_yyyyyy_0_xxxyyyyz_1, g_yyyyyy_0_xxxyyyz_1, g_yyyyyy_0_xxxyyyzz_1, g_yyyyyy_0_xxxyyzz_1, g_yyyyyy_0_xxxyyzzz_1, g_yyyyyy_0_xxxyzzz_1, g_yyyyyy_0_xxxyzzzz_1, g_yyyyyy_0_xxxzzzz_1, g_yyyyyy_0_xxxzzzzz_1, g_yyyyyy_0_xxyyyyy_1, g_yyyyyy_0_xxyyyyyy_1, g_yyyyyy_0_xxyyyyyz_1, g_yyyyyy_0_xxyyyyz_1, g_yyyyyy_0_xxyyyyzz_1, g_yyyyyy_0_xxyyyzz_1, g_yyyyyy_0_xxyyyzzz_1, g_yyyyyy_0_xxyyzzz_1, g_yyyyyy_0_xxyyzzzz_1, g_yyyyyy_0_xxyzzzz_1, g_yyyyyy_0_xxyzzzzz_1, g_yyyyyy_0_xxzzzzz_1, g_yyyyyy_0_xxzzzzzz_1, g_yyyyyy_0_xyyyyyy_1, g_yyyyyy_0_xyyyyyyy_1, g_yyyyyy_0_xyyyyyyz_1, g_yyyyyy_0_xyyyyyz_1, g_yyyyyy_0_xyyyyyzz_1, g_yyyyyy_0_xyyyyzz_1, g_yyyyyy_0_xyyyyzzz_1, g_yyyyyy_0_xyyyzzz_1, g_yyyyyy_0_xyyyzzzz_1, g_yyyyyy_0_xyyzzzz_1, g_yyyyyy_0_xyyzzzzz_1, g_yyyyyy_0_xyzzzzz_1, g_yyyyyy_0_xyzzzzzz_1, g_yyyyyy_0_xzzzzzz_1, g_yyyyyy_0_xzzzzzzz_1, g_yyyyyy_0_yyyyyyy_1, g_yyyyyy_0_yyyyyyyy_1, g_yyyyyy_0_yyyyyyyz_1, g_yyyyyy_0_yyyyyyz_1, g_yyyyyy_0_yyyyyyzz_1, g_yyyyyy_0_yyyyyzz_1, g_yyyyyy_0_yyyyyzzz_1, g_yyyyyy_0_yyyyzzz_1, g_yyyyyy_0_yyyyzzzz_1, g_yyyyyy_0_yyyzzzz_1, g_yyyyyy_0_yyyzzzzz_1, g_yyyyyy_0_yyzzzzz_1, g_yyyyyy_0_yyzzzzzz_1, g_yyyyyy_0_yzzzzzz_1, g_yyyyyy_0_yzzzzzzz_1, g_yyyyyy_0_zzzzzzz_1, g_yyyyyy_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyy_0_xxxxxxxx_0[i] = 8.0 * g_yyyyyy_0_xxxxxxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxxx_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxxxy_0[i] = 7.0 * g_yyyyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxxxz_0[i] = 7.0 * g_yyyyyy_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxxyy_0[i] = 6.0 * g_yyyyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxxyz_0[i] = 6.0 * g_yyyyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxxzz_0[i] = 6.0 * g_yyyyyy_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxyyy_0[i] = 5.0 * g_yyyyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxyyz_0[i] = 5.0 * g_yyyyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxyzz_0[i] = 5.0 * g_yyyyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxzzz_0[i] = 5.0 * g_yyyyyy_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxyyyy_0[i] = 4.0 * g_yyyyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxyyyz_0[i] = 4.0 * g_yyyyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxyyzz_0[i] = 4.0 * g_yyyyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxyzzz_0[i] = 4.0 * g_yyyyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxzzzz_0[i] = 4.0 * g_yyyyyy_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyyyyy_0[i] = 3.0 * g_yyyyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyyyyz_0[i] = 3.0 * g_yyyyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyyyzz_0[i] = 3.0 * g_yyyyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyyzzz_0[i] = 3.0 * g_yyyyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyzzzz_0[i] = 3.0 * g_yyyyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxzzzzz_0[i] = 3.0 * g_yyyyyy_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyyyyy_0[i] = 2.0 * g_yyyyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyyyyz_0[i] = 2.0 * g_yyyyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyyyzz_0[i] = 2.0 * g_yyyyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyyzzz_0[i] = 2.0 * g_yyyyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyzzzz_0[i] = 2.0 * g_yyyyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyzzzzz_0[i] = 2.0 * g_yyyyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxzzzzzz_0[i] = 2.0 * g_yyyyyy_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyyyyy_0[i] = g_yyyyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyyyyz_0[i] = g_yyyyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyyyzz_0[i] = g_yyyyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyyzzz_0[i] = g_yyyyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyzzzz_0[i] = g_yyyyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyzzzzz_0[i] = g_yyyyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyzzzzzz_0[i] = g_yyyyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xzzzzzzz_0[i] = g_yyyyyy_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyyyyy_0[i] = g_yyyyyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyyyyz_0[i] = g_yyyyyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyyyzz_0[i] = g_yyyyyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyyzzz_0[i] = g_yyyyyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyzzzz_0[i] = g_yyyyyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyzzzzz_0[i] = g_yyyyyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyzzzzzz_0[i] = g_yyyyyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yzzzzzzz_0[i] = g_yyyyyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_zzzzzzzz_0[i] = g_yyyyyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 990-1035 components of targeted buffer : KSL

    auto g_xyyyyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 990);

    auto g_xyyyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 991);

    auto g_xyyyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 992);

    auto g_xyyyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 993);

    auto g_xyyyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 994);

    auto g_xyyyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 995);

    auto g_xyyyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 996);

    auto g_xyyyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 997);

    auto g_xyyyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 998);

    auto g_xyyyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 999);

    auto g_xyyyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1000);

    auto g_xyyyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1001);

    auto g_xyyyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1002);

    auto g_xyyyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1003);

    auto g_xyyyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1004);

    auto g_xyyyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1005);

    auto g_xyyyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1006);

    auto g_xyyyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1007);

    auto g_xyyyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1008);

    auto g_xyyyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1009);

    auto g_xyyyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1010);

    auto g_xyyyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1011);

    auto g_xyyyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1012);

    auto g_xyyyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1013);

    auto g_xyyyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1014);

    auto g_xyyyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1015);

    auto g_xyyyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1016);

    auto g_xyyyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1017);

    auto g_xyyyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1018);

    auto g_xyyyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1019);

    auto g_xyyyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1020);

    auto g_xyyyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1021);

    auto g_xyyyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1022);

    auto g_xyyyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1023);

    auto g_xyyyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1024);

    auto g_xyyyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1025);

    auto g_xyyyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1026);

    auto g_xyyyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1027);

    auto g_xyyyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1028);

    auto g_xyyyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1029);

    auto g_xyyyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1030);

    auto g_xyyyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1031);

    auto g_xyyyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1032);

    auto g_xyyyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1033);

    auto g_xyyyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1034);

    #pragma omp simd aligned(g_xyyyyy_0_xxxxxxxx_1, g_xyyyyy_0_xxxxxxxy_1, g_xyyyyy_0_xxxxxxyy_1, g_xyyyyy_0_xxxxxyyy_1, g_xyyyyy_0_xxxxyyyy_1, g_xyyyyy_0_xxxyyyyy_1, g_xyyyyy_0_xxyyyyyy_1, g_xyyyyy_0_xyyyyyyy_1, g_xyyyyyz_0_xxxxxxxx_0, g_xyyyyyz_0_xxxxxxxy_0, g_xyyyyyz_0_xxxxxxxz_0, g_xyyyyyz_0_xxxxxxyy_0, g_xyyyyyz_0_xxxxxxyz_0, g_xyyyyyz_0_xxxxxxzz_0, g_xyyyyyz_0_xxxxxyyy_0, g_xyyyyyz_0_xxxxxyyz_0, g_xyyyyyz_0_xxxxxyzz_0, g_xyyyyyz_0_xxxxxzzz_0, g_xyyyyyz_0_xxxxyyyy_0, g_xyyyyyz_0_xxxxyyyz_0, g_xyyyyyz_0_xxxxyyzz_0, g_xyyyyyz_0_xxxxyzzz_0, g_xyyyyyz_0_xxxxzzzz_0, g_xyyyyyz_0_xxxyyyyy_0, g_xyyyyyz_0_xxxyyyyz_0, g_xyyyyyz_0_xxxyyyzz_0, g_xyyyyyz_0_xxxyyzzz_0, g_xyyyyyz_0_xxxyzzzz_0, g_xyyyyyz_0_xxxzzzzz_0, g_xyyyyyz_0_xxyyyyyy_0, g_xyyyyyz_0_xxyyyyyz_0, g_xyyyyyz_0_xxyyyyzz_0, g_xyyyyyz_0_xxyyyzzz_0, g_xyyyyyz_0_xxyyzzzz_0, g_xyyyyyz_0_xxyzzzzz_0, g_xyyyyyz_0_xxzzzzzz_0, g_xyyyyyz_0_xyyyyyyy_0, g_xyyyyyz_0_xyyyyyyz_0, g_xyyyyyz_0_xyyyyyzz_0, g_xyyyyyz_0_xyyyyzzz_0, g_xyyyyyz_0_xyyyzzzz_0, g_xyyyyyz_0_xyyzzzzz_0, g_xyyyyyz_0_xyzzzzzz_0, g_xyyyyyz_0_xzzzzzzz_0, g_xyyyyyz_0_yyyyyyyy_0, g_xyyyyyz_0_yyyyyyyz_0, g_xyyyyyz_0_yyyyyyzz_0, g_xyyyyyz_0_yyyyyzzz_0, g_xyyyyyz_0_yyyyzzzz_0, g_xyyyyyz_0_yyyzzzzz_0, g_xyyyyyz_0_yyzzzzzz_0, g_xyyyyyz_0_yzzzzzzz_0, g_xyyyyyz_0_zzzzzzzz_0, g_yyyyyz_0_xxxxxxxz_1, g_yyyyyz_0_xxxxxxyz_1, g_yyyyyz_0_xxxxxxz_1, g_yyyyyz_0_xxxxxxzz_1, g_yyyyyz_0_xxxxxyyz_1, g_yyyyyz_0_xxxxxyz_1, g_yyyyyz_0_xxxxxyzz_1, g_yyyyyz_0_xxxxxzz_1, g_yyyyyz_0_xxxxxzzz_1, g_yyyyyz_0_xxxxyyyz_1, g_yyyyyz_0_xxxxyyz_1, g_yyyyyz_0_xxxxyyzz_1, g_yyyyyz_0_xxxxyzz_1, g_yyyyyz_0_xxxxyzzz_1, g_yyyyyz_0_xxxxzzz_1, g_yyyyyz_0_xxxxzzzz_1, g_yyyyyz_0_xxxyyyyz_1, g_yyyyyz_0_xxxyyyz_1, g_yyyyyz_0_xxxyyyzz_1, g_yyyyyz_0_xxxyyzz_1, g_yyyyyz_0_xxxyyzzz_1, g_yyyyyz_0_xxxyzzz_1, g_yyyyyz_0_xxxyzzzz_1, g_yyyyyz_0_xxxzzzz_1, g_yyyyyz_0_xxxzzzzz_1, g_yyyyyz_0_xxyyyyyz_1, g_yyyyyz_0_xxyyyyz_1, g_yyyyyz_0_xxyyyyzz_1, g_yyyyyz_0_xxyyyzz_1, g_yyyyyz_0_xxyyyzzz_1, g_yyyyyz_0_xxyyzzz_1, g_yyyyyz_0_xxyyzzzz_1, g_yyyyyz_0_xxyzzzz_1, g_yyyyyz_0_xxyzzzzz_1, g_yyyyyz_0_xxzzzzz_1, g_yyyyyz_0_xxzzzzzz_1, g_yyyyyz_0_xyyyyyyz_1, g_yyyyyz_0_xyyyyyz_1, g_yyyyyz_0_xyyyyyzz_1, g_yyyyyz_0_xyyyyzz_1, g_yyyyyz_0_xyyyyzzz_1, g_yyyyyz_0_xyyyzzz_1, g_yyyyyz_0_xyyyzzzz_1, g_yyyyyz_0_xyyzzzz_1, g_yyyyyz_0_xyyzzzzz_1, g_yyyyyz_0_xyzzzzz_1, g_yyyyyz_0_xyzzzzzz_1, g_yyyyyz_0_xzzzzzz_1, g_yyyyyz_0_xzzzzzzz_1, g_yyyyyz_0_yyyyyyyy_1, g_yyyyyz_0_yyyyyyyz_1, g_yyyyyz_0_yyyyyyz_1, g_yyyyyz_0_yyyyyyzz_1, g_yyyyyz_0_yyyyyzz_1, g_yyyyyz_0_yyyyyzzz_1, g_yyyyyz_0_yyyyzzz_1, g_yyyyyz_0_yyyyzzzz_1, g_yyyyyz_0_yyyzzzz_1, g_yyyyyz_0_yyyzzzzz_1, g_yyyyyz_0_yyzzzzz_1, g_yyyyyz_0_yyzzzzzz_1, g_yyyyyz_0_yzzzzzz_1, g_yyyyyz_0_yzzzzzzz_1, g_yyyyyz_0_zzzzzzz_1, g_yyyyyz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyz_0_xxxxxxxx_0[i] = g_xyyyyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxxxxy_0[i] = g_xyyyyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxxxxz_0[i] = 7.0 * g_yyyyyz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxxxyy_0[i] = g_xyyyyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxxxyz_0[i] = 6.0 * g_yyyyyz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxxxzz_0[i] = 6.0 * g_yyyyyz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxxyyy_0[i] = g_xyyyyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxxyyz_0[i] = 5.0 * g_yyyyyz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxxyzz_0[i] = 5.0 * g_yyyyyz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxxzzz_0[i] = 5.0 * g_yyyyyz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxyyyy_0[i] = g_xyyyyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxyyyz_0[i] = 4.0 * g_yyyyyz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxyyzz_0[i] = 4.0 * g_yyyyyz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxyzzz_0[i] = 4.0 * g_yyyyyz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxzzzz_0[i] = 4.0 * g_yyyyyz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxyyyyy_0[i] = g_xyyyyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxyyyyz_0[i] = 3.0 * g_yyyyyz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxyyyzz_0[i] = 3.0 * g_yyyyyz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxyyzzz_0[i] = 3.0 * g_yyyyyz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxyzzzz_0[i] = 3.0 * g_yyyyyz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxzzzzz_0[i] = 3.0 * g_yyyyyz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyyyyyy_0[i] = g_xyyyyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxyyyyyz_0[i] = 2.0 * g_yyyyyz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyyyyzz_0[i] = 2.0 * g_yyyyyz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyyyzzz_0[i] = 2.0 * g_yyyyyz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyyzzzz_0[i] = 2.0 * g_yyyyyz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyzzzzz_0[i] = 2.0 * g_yyyyyz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxzzzzzz_0[i] = 2.0 * g_yyyyyz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyyyyyy_0[i] = g_xyyyyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xyyyyyyz_0[i] = g_yyyyyz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyyyyzz_0[i] = g_yyyyyz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyyyzzz_0[i] = g_yyyyyz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyyzzzz_0[i] = g_yyyyyz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyzzzzz_0[i] = g_yyyyyz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyzzzzzz_0[i] = g_yyyyyz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xzzzzzzz_0[i] = g_yyyyyz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyyyyy_0[i] = g_yyyyyz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyyyyz_0[i] = g_yyyyyz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyyyzz_0[i] = g_yyyyyz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyyzzz_0[i] = g_yyyyyz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyzzzz_0[i] = g_yyyyyz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyzzzzz_0[i] = g_yyyyyz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyzzzzzz_0[i] = g_yyyyyz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yzzzzzzz_0[i] = g_yyyyyz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_zzzzzzzz_0[i] = g_yyyyyz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 1035-1080 components of targeted buffer : KSL

    auto g_xyyyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1035);

    auto g_xyyyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1036);

    auto g_xyyyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1037);

    auto g_xyyyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1038);

    auto g_xyyyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1039);

    auto g_xyyyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1040);

    auto g_xyyyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1041);

    auto g_xyyyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1042);

    auto g_xyyyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1043);

    auto g_xyyyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1044);

    auto g_xyyyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1045);

    auto g_xyyyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1046);

    auto g_xyyyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1047);

    auto g_xyyyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1048);

    auto g_xyyyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1049);

    auto g_xyyyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1050);

    auto g_xyyyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1051);

    auto g_xyyyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1052);

    auto g_xyyyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1053);

    auto g_xyyyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1054);

    auto g_xyyyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1055);

    auto g_xyyyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1056);

    auto g_xyyyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1057);

    auto g_xyyyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1058);

    auto g_xyyyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1059);

    auto g_xyyyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1060);

    auto g_xyyyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1061);

    auto g_xyyyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1062);

    auto g_xyyyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1063);

    auto g_xyyyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1064);

    auto g_xyyyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1065);

    auto g_xyyyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1066);

    auto g_xyyyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1067);

    auto g_xyyyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1068);

    auto g_xyyyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1069);

    auto g_xyyyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1070);

    auto g_xyyyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1071);

    auto g_xyyyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1072);

    auto g_xyyyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1073);

    auto g_xyyyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1074);

    auto g_xyyyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1075);

    auto g_xyyyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1076);

    auto g_xyyyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1077);

    auto g_xyyyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1078);

    auto g_xyyyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1079);

    #pragma omp simd aligned(g_xyyyyzz_0_xxxxxxxx_0, g_xyyyyzz_0_xxxxxxxy_0, g_xyyyyzz_0_xxxxxxxz_0, g_xyyyyzz_0_xxxxxxyy_0, g_xyyyyzz_0_xxxxxxyz_0, g_xyyyyzz_0_xxxxxxzz_0, g_xyyyyzz_0_xxxxxyyy_0, g_xyyyyzz_0_xxxxxyyz_0, g_xyyyyzz_0_xxxxxyzz_0, g_xyyyyzz_0_xxxxxzzz_0, g_xyyyyzz_0_xxxxyyyy_0, g_xyyyyzz_0_xxxxyyyz_0, g_xyyyyzz_0_xxxxyyzz_0, g_xyyyyzz_0_xxxxyzzz_0, g_xyyyyzz_0_xxxxzzzz_0, g_xyyyyzz_0_xxxyyyyy_0, g_xyyyyzz_0_xxxyyyyz_0, g_xyyyyzz_0_xxxyyyzz_0, g_xyyyyzz_0_xxxyyzzz_0, g_xyyyyzz_0_xxxyzzzz_0, g_xyyyyzz_0_xxxzzzzz_0, g_xyyyyzz_0_xxyyyyyy_0, g_xyyyyzz_0_xxyyyyyz_0, g_xyyyyzz_0_xxyyyyzz_0, g_xyyyyzz_0_xxyyyzzz_0, g_xyyyyzz_0_xxyyzzzz_0, g_xyyyyzz_0_xxyzzzzz_0, g_xyyyyzz_0_xxzzzzzz_0, g_xyyyyzz_0_xyyyyyyy_0, g_xyyyyzz_0_xyyyyyyz_0, g_xyyyyzz_0_xyyyyyzz_0, g_xyyyyzz_0_xyyyyzzz_0, g_xyyyyzz_0_xyyyzzzz_0, g_xyyyyzz_0_xyyzzzzz_0, g_xyyyyzz_0_xyzzzzzz_0, g_xyyyyzz_0_xzzzzzzz_0, g_xyyyyzz_0_yyyyyyyy_0, g_xyyyyzz_0_yyyyyyyz_0, g_xyyyyzz_0_yyyyyyzz_0, g_xyyyyzz_0_yyyyyzzz_0, g_xyyyyzz_0_yyyyzzzz_0, g_xyyyyzz_0_yyyzzzzz_0, g_xyyyyzz_0_yyzzzzzz_0, g_xyyyyzz_0_yzzzzzzz_0, g_xyyyyzz_0_zzzzzzzz_0, g_yyyyzz_0_xxxxxxx_1, g_yyyyzz_0_xxxxxxxx_1, g_yyyyzz_0_xxxxxxxy_1, g_yyyyzz_0_xxxxxxxz_1, g_yyyyzz_0_xxxxxxy_1, g_yyyyzz_0_xxxxxxyy_1, g_yyyyzz_0_xxxxxxyz_1, g_yyyyzz_0_xxxxxxz_1, g_yyyyzz_0_xxxxxxzz_1, g_yyyyzz_0_xxxxxyy_1, g_yyyyzz_0_xxxxxyyy_1, g_yyyyzz_0_xxxxxyyz_1, g_yyyyzz_0_xxxxxyz_1, g_yyyyzz_0_xxxxxyzz_1, g_yyyyzz_0_xxxxxzz_1, g_yyyyzz_0_xxxxxzzz_1, g_yyyyzz_0_xxxxyyy_1, g_yyyyzz_0_xxxxyyyy_1, g_yyyyzz_0_xxxxyyyz_1, g_yyyyzz_0_xxxxyyz_1, g_yyyyzz_0_xxxxyyzz_1, g_yyyyzz_0_xxxxyzz_1, g_yyyyzz_0_xxxxyzzz_1, g_yyyyzz_0_xxxxzzz_1, g_yyyyzz_0_xxxxzzzz_1, g_yyyyzz_0_xxxyyyy_1, g_yyyyzz_0_xxxyyyyy_1, g_yyyyzz_0_xxxyyyyz_1, g_yyyyzz_0_xxxyyyz_1, g_yyyyzz_0_xxxyyyzz_1, g_yyyyzz_0_xxxyyzz_1, g_yyyyzz_0_xxxyyzzz_1, g_yyyyzz_0_xxxyzzz_1, g_yyyyzz_0_xxxyzzzz_1, g_yyyyzz_0_xxxzzzz_1, g_yyyyzz_0_xxxzzzzz_1, g_yyyyzz_0_xxyyyyy_1, g_yyyyzz_0_xxyyyyyy_1, g_yyyyzz_0_xxyyyyyz_1, g_yyyyzz_0_xxyyyyz_1, g_yyyyzz_0_xxyyyyzz_1, g_yyyyzz_0_xxyyyzz_1, g_yyyyzz_0_xxyyyzzz_1, g_yyyyzz_0_xxyyzzz_1, g_yyyyzz_0_xxyyzzzz_1, g_yyyyzz_0_xxyzzzz_1, g_yyyyzz_0_xxyzzzzz_1, g_yyyyzz_0_xxzzzzz_1, g_yyyyzz_0_xxzzzzzz_1, g_yyyyzz_0_xyyyyyy_1, g_yyyyzz_0_xyyyyyyy_1, g_yyyyzz_0_xyyyyyyz_1, g_yyyyzz_0_xyyyyyz_1, g_yyyyzz_0_xyyyyyzz_1, g_yyyyzz_0_xyyyyzz_1, g_yyyyzz_0_xyyyyzzz_1, g_yyyyzz_0_xyyyzzz_1, g_yyyyzz_0_xyyyzzzz_1, g_yyyyzz_0_xyyzzzz_1, g_yyyyzz_0_xyyzzzzz_1, g_yyyyzz_0_xyzzzzz_1, g_yyyyzz_0_xyzzzzzz_1, g_yyyyzz_0_xzzzzzz_1, g_yyyyzz_0_xzzzzzzz_1, g_yyyyzz_0_yyyyyyy_1, g_yyyyzz_0_yyyyyyyy_1, g_yyyyzz_0_yyyyyyyz_1, g_yyyyzz_0_yyyyyyz_1, g_yyyyzz_0_yyyyyyzz_1, g_yyyyzz_0_yyyyyzz_1, g_yyyyzz_0_yyyyyzzz_1, g_yyyyzz_0_yyyyzzz_1, g_yyyyzz_0_yyyyzzzz_1, g_yyyyzz_0_yyyzzzz_1, g_yyyyzz_0_yyyzzzzz_1, g_yyyyzz_0_yyzzzzz_1, g_yyyyzz_0_yyzzzzzz_1, g_yyyyzz_0_yzzzzzz_1, g_yyyyzz_0_yzzzzzzz_1, g_yyyyzz_0_zzzzzzz_1, g_yyyyzz_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyzz_0_xxxxxxxx_0[i] = 8.0 * g_yyyyzz_0_xxxxxxx_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxxxx_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxxxy_0[i] = 7.0 * g_yyyyzz_0_xxxxxxy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxxxz_0[i] = 7.0 * g_yyyyzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxxyy_0[i] = 6.0 * g_yyyyzz_0_xxxxxyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxxyz_0[i] = 6.0 * g_yyyyzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxxzz_0[i] = 6.0 * g_yyyyzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxyyy_0[i] = 5.0 * g_yyyyzz_0_xxxxyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxyyz_0[i] = 5.0 * g_yyyyzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxyzz_0[i] = 5.0 * g_yyyyzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxzzz_0[i] = 5.0 * g_yyyyzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxyyyy_0[i] = 4.0 * g_yyyyzz_0_xxxyyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxyyyz_0[i] = 4.0 * g_yyyyzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxyyzz_0[i] = 4.0 * g_yyyyzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxyzzz_0[i] = 4.0 * g_yyyyzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxzzzz_0[i] = 4.0 * g_yyyyzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyyyyy_0[i] = 3.0 * g_yyyyzz_0_xxyyyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyyyyz_0[i] = 3.0 * g_yyyyzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyyyzz_0[i] = 3.0 * g_yyyyzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyyzzz_0[i] = 3.0 * g_yyyyzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyzzzz_0[i] = 3.0 * g_yyyyzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxzzzzz_0[i] = 3.0 * g_yyyyzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyyyyy_0[i] = 2.0 * g_yyyyzz_0_xyyyyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyyyyz_0[i] = 2.0 * g_yyyyzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyyyzz_0[i] = 2.0 * g_yyyyzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyyzzz_0[i] = 2.0 * g_yyyyzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyzzzz_0[i] = 2.0 * g_yyyyzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyzzzzz_0[i] = 2.0 * g_yyyyzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxzzzzzz_0[i] = 2.0 * g_yyyyzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyyyyy_0[i] = g_yyyyzz_0_yyyyyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyyyyz_0[i] = g_yyyyzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyyyzz_0[i] = g_yyyyzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyyzzz_0[i] = g_yyyyzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyzzzz_0[i] = g_yyyyzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyzzzzz_0[i] = g_yyyyzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyzzzzzz_0[i] = g_yyyyzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xzzzzzzz_0[i] = g_yyyyzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyyyyy_0[i] = g_yyyyzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyyyyz_0[i] = g_yyyyzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyyyzz_0[i] = g_yyyyzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyyzzz_0[i] = g_yyyyzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyzzzz_0[i] = g_yyyyzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyzzzzz_0[i] = g_yyyyzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyzzzzzz_0[i] = g_yyyyzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yzzzzzzz_0[i] = g_yyyyzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_zzzzzzzz_0[i] = g_yyyyzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 1080-1125 components of targeted buffer : KSL

    auto g_xyyyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1080);

    auto g_xyyyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1081);

    auto g_xyyyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1082);

    auto g_xyyyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1083);

    auto g_xyyyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1084);

    auto g_xyyyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1085);

    auto g_xyyyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1086);

    auto g_xyyyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1087);

    auto g_xyyyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1088);

    auto g_xyyyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1089);

    auto g_xyyyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1090);

    auto g_xyyyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1091);

    auto g_xyyyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1092);

    auto g_xyyyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1093);

    auto g_xyyyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1094);

    auto g_xyyyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1095);

    auto g_xyyyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1096);

    auto g_xyyyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1097);

    auto g_xyyyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1098);

    auto g_xyyyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1099);

    auto g_xyyyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1100);

    auto g_xyyyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1101);

    auto g_xyyyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1102);

    auto g_xyyyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1103);

    auto g_xyyyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1104);

    auto g_xyyyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1105);

    auto g_xyyyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1106);

    auto g_xyyyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1107);

    auto g_xyyyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1108);

    auto g_xyyyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1109);

    auto g_xyyyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1110);

    auto g_xyyyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1111);

    auto g_xyyyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1112);

    auto g_xyyyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1113);

    auto g_xyyyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1114);

    auto g_xyyyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1115);

    auto g_xyyyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1116);

    auto g_xyyyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1117);

    auto g_xyyyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1118);

    auto g_xyyyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1119);

    auto g_xyyyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1120);

    auto g_xyyyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1121);

    auto g_xyyyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1122);

    auto g_xyyyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1123);

    auto g_xyyyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1124);

    #pragma omp simd aligned(g_xyyyzzz_0_xxxxxxxx_0, g_xyyyzzz_0_xxxxxxxy_0, g_xyyyzzz_0_xxxxxxxz_0, g_xyyyzzz_0_xxxxxxyy_0, g_xyyyzzz_0_xxxxxxyz_0, g_xyyyzzz_0_xxxxxxzz_0, g_xyyyzzz_0_xxxxxyyy_0, g_xyyyzzz_0_xxxxxyyz_0, g_xyyyzzz_0_xxxxxyzz_0, g_xyyyzzz_0_xxxxxzzz_0, g_xyyyzzz_0_xxxxyyyy_0, g_xyyyzzz_0_xxxxyyyz_0, g_xyyyzzz_0_xxxxyyzz_0, g_xyyyzzz_0_xxxxyzzz_0, g_xyyyzzz_0_xxxxzzzz_0, g_xyyyzzz_0_xxxyyyyy_0, g_xyyyzzz_0_xxxyyyyz_0, g_xyyyzzz_0_xxxyyyzz_0, g_xyyyzzz_0_xxxyyzzz_0, g_xyyyzzz_0_xxxyzzzz_0, g_xyyyzzz_0_xxxzzzzz_0, g_xyyyzzz_0_xxyyyyyy_0, g_xyyyzzz_0_xxyyyyyz_0, g_xyyyzzz_0_xxyyyyzz_0, g_xyyyzzz_0_xxyyyzzz_0, g_xyyyzzz_0_xxyyzzzz_0, g_xyyyzzz_0_xxyzzzzz_0, g_xyyyzzz_0_xxzzzzzz_0, g_xyyyzzz_0_xyyyyyyy_0, g_xyyyzzz_0_xyyyyyyz_0, g_xyyyzzz_0_xyyyyyzz_0, g_xyyyzzz_0_xyyyyzzz_0, g_xyyyzzz_0_xyyyzzzz_0, g_xyyyzzz_0_xyyzzzzz_0, g_xyyyzzz_0_xyzzzzzz_0, g_xyyyzzz_0_xzzzzzzz_0, g_xyyyzzz_0_yyyyyyyy_0, g_xyyyzzz_0_yyyyyyyz_0, g_xyyyzzz_0_yyyyyyzz_0, g_xyyyzzz_0_yyyyyzzz_0, g_xyyyzzz_0_yyyyzzzz_0, g_xyyyzzz_0_yyyzzzzz_0, g_xyyyzzz_0_yyzzzzzz_0, g_xyyyzzz_0_yzzzzzzz_0, g_xyyyzzz_0_zzzzzzzz_0, g_yyyzzz_0_xxxxxxx_1, g_yyyzzz_0_xxxxxxxx_1, g_yyyzzz_0_xxxxxxxy_1, g_yyyzzz_0_xxxxxxxz_1, g_yyyzzz_0_xxxxxxy_1, g_yyyzzz_0_xxxxxxyy_1, g_yyyzzz_0_xxxxxxyz_1, g_yyyzzz_0_xxxxxxz_1, g_yyyzzz_0_xxxxxxzz_1, g_yyyzzz_0_xxxxxyy_1, g_yyyzzz_0_xxxxxyyy_1, g_yyyzzz_0_xxxxxyyz_1, g_yyyzzz_0_xxxxxyz_1, g_yyyzzz_0_xxxxxyzz_1, g_yyyzzz_0_xxxxxzz_1, g_yyyzzz_0_xxxxxzzz_1, g_yyyzzz_0_xxxxyyy_1, g_yyyzzz_0_xxxxyyyy_1, g_yyyzzz_0_xxxxyyyz_1, g_yyyzzz_0_xxxxyyz_1, g_yyyzzz_0_xxxxyyzz_1, g_yyyzzz_0_xxxxyzz_1, g_yyyzzz_0_xxxxyzzz_1, g_yyyzzz_0_xxxxzzz_1, g_yyyzzz_0_xxxxzzzz_1, g_yyyzzz_0_xxxyyyy_1, g_yyyzzz_0_xxxyyyyy_1, g_yyyzzz_0_xxxyyyyz_1, g_yyyzzz_0_xxxyyyz_1, g_yyyzzz_0_xxxyyyzz_1, g_yyyzzz_0_xxxyyzz_1, g_yyyzzz_0_xxxyyzzz_1, g_yyyzzz_0_xxxyzzz_1, g_yyyzzz_0_xxxyzzzz_1, g_yyyzzz_0_xxxzzzz_1, g_yyyzzz_0_xxxzzzzz_1, g_yyyzzz_0_xxyyyyy_1, g_yyyzzz_0_xxyyyyyy_1, g_yyyzzz_0_xxyyyyyz_1, g_yyyzzz_0_xxyyyyz_1, g_yyyzzz_0_xxyyyyzz_1, g_yyyzzz_0_xxyyyzz_1, g_yyyzzz_0_xxyyyzzz_1, g_yyyzzz_0_xxyyzzz_1, g_yyyzzz_0_xxyyzzzz_1, g_yyyzzz_0_xxyzzzz_1, g_yyyzzz_0_xxyzzzzz_1, g_yyyzzz_0_xxzzzzz_1, g_yyyzzz_0_xxzzzzzz_1, g_yyyzzz_0_xyyyyyy_1, g_yyyzzz_0_xyyyyyyy_1, g_yyyzzz_0_xyyyyyyz_1, g_yyyzzz_0_xyyyyyz_1, g_yyyzzz_0_xyyyyyzz_1, g_yyyzzz_0_xyyyyzz_1, g_yyyzzz_0_xyyyyzzz_1, g_yyyzzz_0_xyyyzzz_1, g_yyyzzz_0_xyyyzzzz_1, g_yyyzzz_0_xyyzzzz_1, g_yyyzzz_0_xyyzzzzz_1, g_yyyzzz_0_xyzzzzz_1, g_yyyzzz_0_xyzzzzzz_1, g_yyyzzz_0_xzzzzzz_1, g_yyyzzz_0_xzzzzzzz_1, g_yyyzzz_0_yyyyyyy_1, g_yyyzzz_0_yyyyyyyy_1, g_yyyzzz_0_yyyyyyyz_1, g_yyyzzz_0_yyyyyyz_1, g_yyyzzz_0_yyyyyyzz_1, g_yyyzzz_0_yyyyyzz_1, g_yyyzzz_0_yyyyyzzz_1, g_yyyzzz_0_yyyyzzz_1, g_yyyzzz_0_yyyyzzzz_1, g_yyyzzz_0_yyyzzzz_1, g_yyyzzz_0_yyyzzzzz_1, g_yyyzzz_0_yyzzzzz_1, g_yyyzzz_0_yyzzzzzz_1, g_yyyzzz_0_yzzzzzz_1, g_yyyzzz_0_yzzzzzzz_1, g_yyyzzz_0_zzzzzzz_1, g_yyyzzz_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzzz_0_xxxxxxxx_0[i] = 8.0 * g_yyyzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxxxx_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxxxy_0[i] = 7.0 * g_yyyzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxxxz_0[i] = 7.0 * g_yyyzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxxyy_0[i] = 6.0 * g_yyyzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxxyz_0[i] = 6.0 * g_yyyzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxxzz_0[i] = 6.0 * g_yyyzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxyyy_0[i] = 5.0 * g_yyyzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxyyz_0[i] = 5.0 * g_yyyzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxyzz_0[i] = 5.0 * g_yyyzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxzzz_0[i] = 5.0 * g_yyyzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxyyyy_0[i] = 4.0 * g_yyyzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxyyyz_0[i] = 4.0 * g_yyyzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxyyzz_0[i] = 4.0 * g_yyyzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxyzzz_0[i] = 4.0 * g_yyyzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxzzzz_0[i] = 4.0 * g_yyyzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyyyyy_0[i] = 3.0 * g_yyyzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyyyyz_0[i] = 3.0 * g_yyyzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyyyzz_0[i] = 3.0 * g_yyyzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyyzzz_0[i] = 3.0 * g_yyyzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyzzzz_0[i] = 3.0 * g_yyyzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxzzzzz_0[i] = 3.0 * g_yyyzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyyyyy_0[i] = 2.0 * g_yyyzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyyyyz_0[i] = 2.0 * g_yyyzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyyyzz_0[i] = 2.0 * g_yyyzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyyzzz_0[i] = 2.0 * g_yyyzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyzzzz_0[i] = 2.0 * g_yyyzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyzzzzz_0[i] = 2.0 * g_yyyzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxzzzzzz_0[i] = 2.0 * g_yyyzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyyyyy_0[i] = g_yyyzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyyyyz_0[i] = g_yyyzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyyyzz_0[i] = g_yyyzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyyzzz_0[i] = g_yyyzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyzzzz_0[i] = g_yyyzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyzzzzz_0[i] = g_yyyzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyzzzzzz_0[i] = g_yyyzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xzzzzzzz_0[i] = g_yyyzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyyyyy_0[i] = g_yyyzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyyyyz_0[i] = g_yyyzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyyyzz_0[i] = g_yyyzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyyzzz_0[i] = g_yyyzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyzzzz_0[i] = g_yyyzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyzzzzz_0[i] = g_yyyzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyzzzzzz_0[i] = g_yyyzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yzzzzzzz_0[i] = g_yyyzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_zzzzzzzz_0[i] = g_yyyzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 1125-1170 components of targeted buffer : KSL

    auto g_xyyzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1125);

    auto g_xyyzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1126);

    auto g_xyyzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1127);

    auto g_xyyzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1128);

    auto g_xyyzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1129);

    auto g_xyyzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1130);

    auto g_xyyzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1131);

    auto g_xyyzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1132);

    auto g_xyyzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1133);

    auto g_xyyzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1134);

    auto g_xyyzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1135);

    auto g_xyyzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1136);

    auto g_xyyzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1137);

    auto g_xyyzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1138);

    auto g_xyyzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1139);

    auto g_xyyzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1140);

    auto g_xyyzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1141);

    auto g_xyyzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1142);

    auto g_xyyzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1143);

    auto g_xyyzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1144);

    auto g_xyyzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1145);

    auto g_xyyzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1146);

    auto g_xyyzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1147);

    auto g_xyyzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1148);

    auto g_xyyzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1149);

    auto g_xyyzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1150);

    auto g_xyyzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1151);

    auto g_xyyzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1152);

    auto g_xyyzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1153);

    auto g_xyyzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1154);

    auto g_xyyzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1155);

    auto g_xyyzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1156);

    auto g_xyyzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1157);

    auto g_xyyzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1158);

    auto g_xyyzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1159);

    auto g_xyyzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1160);

    auto g_xyyzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1161);

    auto g_xyyzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1162);

    auto g_xyyzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1163);

    auto g_xyyzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1164);

    auto g_xyyzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1165);

    auto g_xyyzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1166);

    auto g_xyyzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1167);

    auto g_xyyzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1168);

    auto g_xyyzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1169);

    #pragma omp simd aligned(g_xyyzzzz_0_xxxxxxxx_0, g_xyyzzzz_0_xxxxxxxy_0, g_xyyzzzz_0_xxxxxxxz_0, g_xyyzzzz_0_xxxxxxyy_0, g_xyyzzzz_0_xxxxxxyz_0, g_xyyzzzz_0_xxxxxxzz_0, g_xyyzzzz_0_xxxxxyyy_0, g_xyyzzzz_0_xxxxxyyz_0, g_xyyzzzz_0_xxxxxyzz_0, g_xyyzzzz_0_xxxxxzzz_0, g_xyyzzzz_0_xxxxyyyy_0, g_xyyzzzz_0_xxxxyyyz_0, g_xyyzzzz_0_xxxxyyzz_0, g_xyyzzzz_0_xxxxyzzz_0, g_xyyzzzz_0_xxxxzzzz_0, g_xyyzzzz_0_xxxyyyyy_0, g_xyyzzzz_0_xxxyyyyz_0, g_xyyzzzz_0_xxxyyyzz_0, g_xyyzzzz_0_xxxyyzzz_0, g_xyyzzzz_0_xxxyzzzz_0, g_xyyzzzz_0_xxxzzzzz_0, g_xyyzzzz_0_xxyyyyyy_0, g_xyyzzzz_0_xxyyyyyz_0, g_xyyzzzz_0_xxyyyyzz_0, g_xyyzzzz_0_xxyyyzzz_0, g_xyyzzzz_0_xxyyzzzz_0, g_xyyzzzz_0_xxyzzzzz_0, g_xyyzzzz_0_xxzzzzzz_0, g_xyyzzzz_0_xyyyyyyy_0, g_xyyzzzz_0_xyyyyyyz_0, g_xyyzzzz_0_xyyyyyzz_0, g_xyyzzzz_0_xyyyyzzz_0, g_xyyzzzz_0_xyyyzzzz_0, g_xyyzzzz_0_xyyzzzzz_0, g_xyyzzzz_0_xyzzzzzz_0, g_xyyzzzz_0_xzzzzzzz_0, g_xyyzzzz_0_yyyyyyyy_0, g_xyyzzzz_0_yyyyyyyz_0, g_xyyzzzz_0_yyyyyyzz_0, g_xyyzzzz_0_yyyyyzzz_0, g_xyyzzzz_0_yyyyzzzz_0, g_xyyzzzz_0_yyyzzzzz_0, g_xyyzzzz_0_yyzzzzzz_0, g_xyyzzzz_0_yzzzzzzz_0, g_xyyzzzz_0_zzzzzzzz_0, g_yyzzzz_0_xxxxxxx_1, g_yyzzzz_0_xxxxxxxx_1, g_yyzzzz_0_xxxxxxxy_1, g_yyzzzz_0_xxxxxxxz_1, g_yyzzzz_0_xxxxxxy_1, g_yyzzzz_0_xxxxxxyy_1, g_yyzzzz_0_xxxxxxyz_1, g_yyzzzz_0_xxxxxxz_1, g_yyzzzz_0_xxxxxxzz_1, g_yyzzzz_0_xxxxxyy_1, g_yyzzzz_0_xxxxxyyy_1, g_yyzzzz_0_xxxxxyyz_1, g_yyzzzz_0_xxxxxyz_1, g_yyzzzz_0_xxxxxyzz_1, g_yyzzzz_0_xxxxxzz_1, g_yyzzzz_0_xxxxxzzz_1, g_yyzzzz_0_xxxxyyy_1, g_yyzzzz_0_xxxxyyyy_1, g_yyzzzz_0_xxxxyyyz_1, g_yyzzzz_0_xxxxyyz_1, g_yyzzzz_0_xxxxyyzz_1, g_yyzzzz_0_xxxxyzz_1, g_yyzzzz_0_xxxxyzzz_1, g_yyzzzz_0_xxxxzzz_1, g_yyzzzz_0_xxxxzzzz_1, g_yyzzzz_0_xxxyyyy_1, g_yyzzzz_0_xxxyyyyy_1, g_yyzzzz_0_xxxyyyyz_1, g_yyzzzz_0_xxxyyyz_1, g_yyzzzz_0_xxxyyyzz_1, g_yyzzzz_0_xxxyyzz_1, g_yyzzzz_0_xxxyyzzz_1, g_yyzzzz_0_xxxyzzz_1, g_yyzzzz_0_xxxyzzzz_1, g_yyzzzz_0_xxxzzzz_1, g_yyzzzz_0_xxxzzzzz_1, g_yyzzzz_0_xxyyyyy_1, g_yyzzzz_0_xxyyyyyy_1, g_yyzzzz_0_xxyyyyyz_1, g_yyzzzz_0_xxyyyyz_1, g_yyzzzz_0_xxyyyyzz_1, g_yyzzzz_0_xxyyyzz_1, g_yyzzzz_0_xxyyyzzz_1, g_yyzzzz_0_xxyyzzz_1, g_yyzzzz_0_xxyyzzzz_1, g_yyzzzz_0_xxyzzzz_1, g_yyzzzz_0_xxyzzzzz_1, g_yyzzzz_0_xxzzzzz_1, g_yyzzzz_0_xxzzzzzz_1, g_yyzzzz_0_xyyyyyy_1, g_yyzzzz_0_xyyyyyyy_1, g_yyzzzz_0_xyyyyyyz_1, g_yyzzzz_0_xyyyyyz_1, g_yyzzzz_0_xyyyyyzz_1, g_yyzzzz_0_xyyyyzz_1, g_yyzzzz_0_xyyyyzzz_1, g_yyzzzz_0_xyyyzzz_1, g_yyzzzz_0_xyyyzzzz_1, g_yyzzzz_0_xyyzzzz_1, g_yyzzzz_0_xyyzzzzz_1, g_yyzzzz_0_xyzzzzz_1, g_yyzzzz_0_xyzzzzzz_1, g_yyzzzz_0_xzzzzzz_1, g_yyzzzz_0_xzzzzzzz_1, g_yyzzzz_0_yyyyyyy_1, g_yyzzzz_0_yyyyyyyy_1, g_yyzzzz_0_yyyyyyyz_1, g_yyzzzz_0_yyyyyyz_1, g_yyzzzz_0_yyyyyyzz_1, g_yyzzzz_0_yyyyyzz_1, g_yyzzzz_0_yyyyyzzz_1, g_yyzzzz_0_yyyyzzz_1, g_yyzzzz_0_yyyyzzzz_1, g_yyzzzz_0_yyyzzzz_1, g_yyzzzz_0_yyyzzzzz_1, g_yyzzzz_0_yyzzzzz_1, g_yyzzzz_0_yyzzzzzz_1, g_yyzzzz_0_yzzzzzz_1, g_yyzzzz_0_yzzzzzzz_1, g_yyzzzz_0_zzzzzzz_1, g_yyzzzz_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzzz_0_xxxxxxxx_0[i] = 8.0 * g_yyzzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxxxx_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxxxy_0[i] = 7.0 * g_yyzzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxxxz_0[i] = 7.0 * g_yyzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxxyy_0[i] = 6.0 * g_yyzzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxxyz_0[i] = 6.0 * g_yyzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxxzz_0[i] = 6.0 * g_yyzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxyyy_0[i] = 5.0 * g_yyzzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxyyz_0[i] = 5.0 * g_yyzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxyzz_0[i] = 5.0 * g_yyzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxzzz_0[i] = 5.0 * g_yyzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxyyyy_0[i] = 4.0 * g_yyzzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxyyyz_0[i] = 4.0 * g_yyzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxyyzz_0[i] = 4.0 * g_yyzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxyzzz_0[i] = 4.0 * g_yyzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxzzzz_0[i] = 4.0 * g_yyzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyyyyy_0[i] = 3.0 * g_yyzzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyyyyz_0[i] = 3.0 * g_yyzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyyyzz_0[i] = 3.0 * g_yyzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyyzzz_0[i] = 3.0 * g_yyzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyzzzz_0[i] = 3.0 * g_yyzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxzzzzz_0[i] = 3.0 * g_yyzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyyyyy_0[i] = 2.0 * g_yyzzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyyyyz_0[i] = 2.0 * g_yyzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyyyzz_0[i] = 2.0 * g_yyzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyyzzz_0[i] = 2.0 * g_yyzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyzzzz_0[i] = 2.0 * g_yyzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyzzzzz_0[i] = 2.0 * g_yyzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxzzzzzz_0[i] = 2.0 * g_yyzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyyyyy_0[i] = g_yyzzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyyyyz_0[i] = g_yyzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyyyzz_0[i] = g_yyzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyyzzz_0[i] = g_yyzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyzzzz_0[i] = g_yyzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyzzzzz_0[i] = g_yyzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyzzzzzz_0[i] = g_yyzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xzzzzzzz_0[i] = g_yyzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyyyyy_0[i] = g_yyzzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyyyyz_0[i] = g_yyzzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyyyzz_0[i] = g_yyzzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyyzzz_0[i] = g_yyzzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyzzzz_0[i] = g_yyzzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyzzzzz_0[i] = g_yyzzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyzzzzzz_0[i] = g_yyzzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yzzzzzzz_0[i] = g_yyzzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_zzzzzzzz_0[i] = g_yyzzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 1170-1215 components of targeted buffer : KSL

    auto g_xyzzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1170);

    auto g_xyzzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1171);

    auto g_xyzzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1172);

    auto g_xyzzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1173);

    auto g_xyzzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1174);

    auto g_xyzzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1175);

    auto g_xyzzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1176);

    auto g_xyzzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1177);

    auto g_xyzzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1178);

    auto g_xyzzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1179);

    auto g_xyzzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1180);

    auto g_xyzzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1181);

    auto g_xyzzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1182);

    auto g_xyzzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1183);

    auto g_xyzzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1184);

    auto g_xyzzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1185);

    auto g_xyzzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1186);

    auto g_xyzzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1187);

    auto g_xyzzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1188);

    auto g_xyzzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1189);

    auto g_xyzzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1190);

    auto g_xyzzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1191);

    auto g_xyzzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1192);

    auto g_xyzzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1193);

    auto g_xyzzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1194);

    auto g_xyzzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1195);

    auto g_xyzzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1196);

    auto g_xyzzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1197);

    auto g_xyzzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1198);

    auto g_xyzzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1199);

    auto g_xyzzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1200);

    auto g_xyzzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1201);

    auto g_xyzzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1202);

    auto g_xyzzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1203);

    auto g_xyzzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1204);

    auto g_xyzzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1205);

    auto g_xyzzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1206);

    auto g_xyzzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1207);

    auto g_xyzzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1208);

    auto g_xyzzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1209);

    auto g_xyzzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1210);

    auto g_xyzzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1211);

    auto g_xyzzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1212);

    auto g_xyzzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1213);

    auto g_xyzzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1214);

    #pragma omp simd aligned(g_xyzzzzz_0_xxxxxxxx_0, g_xyzzzzz_0_xxxxxxxy_0, g_xyzzzzz_0_xxxxxxxz_0, g_xyzzzzz_0_xxxxxxyy_0, g_xyzzzzz_0_xxxxxxyz_0, g_xyzzzzz_0_xxxxxxzz_0, g_xyzzzzz_0_xxxxxyyy_0, g_xyzzzzz_0_xxxxxyyz_0, g_xyzzzzz_0_xxxxxyzz_0, g_xyzzzzz_0_xxxxxzzz_0, g_xyzzzzz_0_xxxxyyyy_0, g_xyzzzzz_0_xxxxyyyz_0, g_xyzzzzz_0_xxxxyyzz_0, g_xyzzzzz_0_xxxxyzzz_0, g_xyzzzzz_0_xxxxzzzz_0, g_xyzzzzz_0_xxxyyyyy_0, g_xyzzzzz_0_xxxyyyyz_0, g_xyzzzzz_0_xxxyyyzz_0, g_xyzzzzz_0_xxxyyzzz_0, g_xyzzzzz_0_xxxyzzzz_0, g_xyzzzzz_0_xxxzzzzz_0, g_xyzzzzz_0_xxyyyyyy_0, g_xyzzzzz_0_xxyyyyyz_0, g_xyzzzzz_0_xxyyyyzz_0, g_xyzzzzz_0_xxyyyzzz_0, g_xyzzzzz_0_xxyyzzzz_0, g_xyzzzzz_0_xxyzzzzz_0, g_xyzzzzz_0_xxzzzzzz_0, g_xyzzzzz_0_xyyyyyyy_0, g_xyzzzzz_0_xyyyyyyz_0, g_xyzzzzz_0_xyyyyyzz_0, g_xyzzzzz_0_xyyyyzzz_0, g_xyzzzzz_0_xyyyzzzz_0, g_xyzzzzz_0_xyyzzzzz_0, g_xyzzzzz_0_xyzzzzzz_0, g_xyzzzzz_0_xzzzzzzz_0, g_xyzzzzz_0_yyyyyyyy_0, g_xyzzzzz_0_yyyyyyyz_0, g_xyzzzzz_0_yyyyyyzz_0, g_xyzzzzz_0_yyyyyzzz_0, g_xyzzzzz_0_yyyyzzzz_0, g_xyzzzzz_0_yyyzzzzz_0, g_xyzzzzz_0_yyzzzzzz_0, g_xyzzzzz_0_yzzzzzzz_0, g_xyzzzzz_0_zzzzzzzz_0, g_xzzzzz_0_xxxxxxxx_1, g_xzzzzz_0_xxxxxxxz_1, g_xzzzzz_0_xxxxxxzz_1, g_xzzzzz_0_xxxxxzzz_1, g_xzzzzz_0_xxxxzzzz_1, g_xzzzzz_0_xxxzzzzz_1, g_xzzzzz_0_xxzzzzzz_1, g_xzzzzz_0_xzzzzzzz_1, g_yzzzzz_0_xxxxxxxy_1, g_yzzzzz_0_xxxxxxy_1, g_yzzzzz_0_xxxxxxyy_1, g_yzzzzz_0_xxxxxxyz_1, g_yzzzzz_0_xxxxxyy_1, g_yzzzzz_0_xxxxxyyy_1, g_yzzzzz_0_xxxxxyyz_1, g_yzzzzz_0_xxxxxyz_1, g_yzzzzz_0_xxxxxyzz_1, g_yzzzzz_0_xxxxyyy_1, g_yzzzzz_0_xxxxyyyy_1, g_yzzzzz_0_xxxxyyyz_1, g_yzzzzz_0_xxxxyyz_1, g_yzzzzz_0_xxxxyyzz_1, g_yzzzzz_0_xxxxyzz_1, g_yzzzzz_0_xxxxyzzz_1, g_yzzzzz_0_xxxyyyy_1, g_yzzzzz_0_xxxyyyyy_1, g_yzzzzz_0_xxxyyyyz_1, g_yzzzzz_0_xxxyyyz_1, g_yzzzzz_0_xxxyyyzz_1, g_yzzzzz_0_xxxyyzz_1, g_yzzzzz_0_xxxyyzzz_1, g_yzzzzz_0_xxxyzzz_1, g_yzzzzz_0_xxxyzzzz_1, g_yzzzzz_0_xxyyyyy_1, g_yzzzzz_0_xxyyyyyy_1, g_yzzzzz_0_xxyyyyyz_1, g_yzzzzz_0_xxyyyyz_1, g_yzzzzz_0_xxyyyyzz_1, g_yzzzzz_0_xxyyyzz_1, g_yzzzzz_0_xxyyyzzz_1, g_yzzzzz_0_xxyyzzz_1, g_yzzzzz_0_xxyyzzzz_1, g_yzzzzz_0_xxyzzzz_1, g_yzzzzz_0_xxyzzzzz_1, g_yzzzzz_0_xyyyyyy_1, g_yzzzzz_0_xyyyyyyy_1, g_yzzzzz_0_xyyyyyyz_1, g_yzzzzz_0_xyyyyyz_1, g_yzzzzz_0_xyyyyyzz_1, g_yzzzzz_0_xyyyyzz_1, g_yzzzzz_0_xyyyyzzz_1, g_yzzzzz_0_xyyyzzz_1, g_yzzzzz_0_xyyyzzzz_1, g_yzzzzz_0_xyyzzzz_1, g_yzzzzz_0_xyyzzzzz_1, g_yzzzzz_0_xyzzzzz_1, g_yzzzzz_0_xyzzzzzz_1, g_yzzzzz_0_yyyyyyy_1, g_yzzzzz_0_yyyyyyyy_1, g_yzzzzz_0_yyyyyyyz_1, g_yzzzzz_0_yyyyyyz_1, g_yzzzzz_0_yyyyyyzz_1, g_yzzzzz_0_yyyyyzz_1, g_yzzzzz_0_yyyyyzzz_1, g_yzzzzz_0_yyyyzzz_1, g_yzzzzz_0_yyyyzzzz_1, g_yzzzzz_0_yyyzzzz_1, g_yzzzzz_0_yyyzzzzz_1, g_yzzzzz_0_yyzzzzz_1, g_yzzzzz_0_yyzzzzzz_1, g_yzzzzz_0_yzzzzzz_1, g_yzzzzz_0_yzzzzzzz_1, g_yzzzzz_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzzz_0_xxxxxxxx_0[i] = g_xzzzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxxxxxy_0[i] = 7.0 * g_yzzzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxxxxz_0[i] = g_xzzzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxxxxyy_0[i] = 6.0 * g_yzzzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxxxyz_0[i] = 6.0 * g_yzzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxxxzz_0[i] = g_xzzzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxxxyyy_0[i] = 5.0 * g_yzzzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxxyyz_0[i] = 5.0 * g_yzzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxxyzz_0[i] = 5.0 * g_yzzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxxzzz_0[i] = g_xzzzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxxyyyy_0[i] = 4.0 * g_yzzzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxyyyz_0[i] = 4.0 * g_yzzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxyyzz_0[i] = 4.0 * g_yzzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxyzzz_0[i] = 4.0 * g_yzzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxzzzz_0[i] = g_xzzzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxyyyyy_0[i] = 3.0 * g_yzzzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxyyyyz_0[i] = 3.0 * g_yzzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxyyyzz_0[i] = 3.0 * g_yzzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxyyzzz_0[i] = 3.0 * g_yzzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxyzzzz_0[i] = 3.0 * g_yzzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxzzzzz_0[i] = g_xzzzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxyyyyyy_0[i] = 2.0 * g_yzzzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyyyyyz_0[i] = 2.0 * g_yzzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyyyyzz_0[i] = 2.0 * g_yzzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyyyzzz_0[i] = 2.0 * g_yzzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyyzzzz_0[i] = 2.0 * g_yzzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyzzzzz_0[i] = 2.0 * g_yzzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxzzzzzz_0[i] = g_xzzzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xyyyyyyy_0[i] = g_yzzzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyyyyyz_0[i] = g_yzzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyyyyzz_0[i] = g_yzzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyyyzzz_0[i] = g_yzzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyyzzzz_0[i] = g_yzzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyzzzzz_0[i] = g_yzzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyzzzzzz_0[i] = g_yzzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xzzzzzzz_0[i] = g_xzzzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_yyyyyyyy_0[i] = g_yzzzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyyyyyz_0[i] = g_yzzzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyyyyzz_0[i] = g_yzzzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyyyzzz_0[i] = g_yzzzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyyzzzz_0[i] = g_yzzzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyzzzzz_0[i] = g_yzzzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyzzzzzz_0[i] = g_yzzzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yzzzzzzz_0[i] = g_yzzzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_zzzzzzzz_0[i] = g_yzzzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 1215-1260 components of targeted buffer : KSL

    auto g_xzzzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1215);

    auto g_xzzzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1216);

    auto g_xzzzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1217);

    auto g_xzzzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1218);

    auto g_xzzzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1219);

    auto g_xzzzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1220);

    auto g_xzzzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1221);

    auto g_xzzzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1222);

    auto g_xzzzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1223);

    auto g_xzzzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1224);

    auto g_xzzzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1225);

    auto g_xzzzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1226);

    auto g_xzzzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1227);

    auto g_xzzzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1228);

    auto g_xzzzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1229);

    auto g_xzzzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1230);

    auto g_xzzzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1231);

    auto g_xzzzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1232);

    auto g_xzzzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1233);

    auto g_xzzzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1234);

    auto g_xzzzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1235);

    auto g_xzzzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1236);

    auto g_xzzzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1237);

    auto g_xzzzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1238);

    auto g_xzzzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1239);

    auto g_xzzzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1240);

    auto g_xzzzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1241);

    auto g_xzzzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1242);

    auto g_xzzzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1243);

    auto g_xzzzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1244);

    auto g_xzzzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1245);

    auto g_xzzzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1246);

    auto g_xzzzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1247);

    auto g_xzzzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1248);

    auto g_xzzzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1249);

    auto g_xzzzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1250);

    auto g_xzzzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1251);

    auto g_xzzzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1252);

    auto g_xzzzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1253);

    auto g_xzzzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1254);

    auto g_xzzzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1255);

    auto g_xzzzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1256);

    auto g_xzzzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1257);

    auto g_xzzzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1258);

    auto g_xzzzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1259);

    #pragma omp simd aligned(g_xzzzzzz_0_xxxxxxxx_0, g_xzzzzzz_0_xxxxxxxy_0, g_xzzzzzz_0_xxxxxxxz_0, g_xzzzzzz_0_xxxxxxyy_0, g_xzzzzzz_0_xxxxxxyz_0, g_xzzzzzz_0_xxxxxxzz_0, g_xzzzzzz_0_xxxxxyyy_0, g_xzzzzzz_0_xxxxxyyz_0, g_xzzzzzz_0_xxxxxyzz_0, g_xzzzzzz_0_xxxxxzzz_0, g_xzzzzzz_0_xxxxyyyy_0, g_xzzzzzz_0_xxxxyyyz_0, g_xzzzzzz_0_xxxxyyzz_0, g_xzzzzzz_0_xxxxyzzz_0, g_xzzzzzz_0_xxxxzzzz_0, g_xzzzzzz_0_xxxyyyyy_0, g_xzzzzzz_0_xxxyyyyz_0, g_xzzzzzz_0_xxxyyyzz_0, g_xzzzzzz_0_xxxyyzzz_0, g_xzzzzzz_0_xxxyzzzz_0, g_xzzzzzz_0_xxxzzzzz_0, g_xzzzzzz_0_xxyyyyyy_0, g_xzzzzzz_0_xxyyyyyz_0, g_xzzzzzz_0_xxyyyyzz_0, g_xzzzzzz_0_xxyyyzzz_0, g_xzzzzzz_0_xxyyzzzz_0, g_xzzzzzz_0_xxyzzzzz_0, g_xzzzzzz_0_xxzzzzzz_0, g_xzzzzzz_0_xyyyyyyy_0, g_xzzzzzz_0_xyyyyyyz_0, g_xzzzzzz_0_xyyyyyzz_0, g_xzzzzzz_0_xyyyyzzz_0, g_xzzzzzz_0_xyyyzzzz_0, g_xzzzzzz_0_xyyzzzzz_0, g_xzzzzzz_0_xyzzzzzz_0, g_xzzzzzz_0_xzzzzzzz_0, g_xzzzzzz_0_yyyyyyyy_0, g_xzzzzzz_0_yyyyyyyz_0, g_xzzzzzz_0_yyyyyyzz_0, g_xzzzzzz_0_yyyyyzzz_0, g_xzzzzzz_0_yyyyzzzz_0, g_xzzzzzz_0_yyyzzzzz_0, g_xzzzzzz_0_yyzzzzzz_0, g_xzzzzzz_0_yzzzzzzz_0, g_xzzzzzz_0_zzzzzzzz_0, g_zzzzzz_0_xxxxxxx_1, g_zzzzzz_0_xxxxxxxx_1, g_zzzzzz_0_xxxxxxxy_1, g_zzzzzz_0_xxxxxxxz_1, g_zzzzzz_0_xxxxxxy_1, g_zzzzzz_0_xxxxxxyy_1, g_zzzzzz_0_xxxxxxyz_1, g_zzzzzz_0_xxxxxxz_1, g_zzzzzz_0_xxxxxxzz_1, g_zzzzzz_0_xxxxxyy_1, g_zzzzzz_0_xxxxxyyy_1, g_zzzzzz_0_xxxxxyyz_1, g_zzzzzz_0_xxxxxyz_1, g_zzzzzz_0_xxxxxyzz_1, g_zzzzzz_0_xxxxxzz_1, g_zzzzzz_0_xxxxxzzz_1, g_zzzzzz_0_xxxxyyy_1, g_zzzzzz_0_xxxxyyyy_1, g_zzzzzz_0_xxxxyyyz_1, g_zzzzzz_0_xxxxyyz_1, g_zzzzzz_0_xxxxyyzz_1, g_zzzzzz_0_xxxxyzz_1, g_zzzzzz_0_xxxxyzzz_1, g_zzzzzz_0_xxxxzzz_1, g_zzzzzz_0_xxxxzzzz_1, g_zzzzzz_0_xxxyyyy_1, g_zzzzzz_0_xxxyyyyy_1, g_zzzzzz_0_xxxyyyyz_1, g_zzzzzz_0_xxxyyyz_1, g_zzzzzz_0_xxxyyyzz_1, g_zzzzzz_0_xxxyyzz_1, g_zzzzzz_0_xxxyyzzz_1, g_zzzzzz_0_xxxyzzz_1, g_zzzzzz_0_xxxyzzzz_1, g_zzzzzz_0_xxxzzzz_1, g_zzzzzz_0_xxxzzzzz_1, g_zzzzzz_0_xxyyyyy_1, g_zzzzzz_0_xxyyyyyy_1, g_zzzzzz_0_xxyyyyyz_1, g_zzzzzz_0_xxyyyyz_1, g_zzzzzz_0_xxyyyyzz_1, g_zzzzzz_0_xxyyyzz_1, g_zzzzzz_0_xxyyyzzz_1, g_zzzzzz_0_xxyyzzz_1, g_zzzzzz_0_xxyyzzzz_1, g_zzzzzz_0_xxyzzzz_1, g_zzzzzz_0_xxyzzzzz_1, g_zzzzzz_0_xxzzzzz_1, g_zzzzzz_0_xxzzzzzz_1, g_zzzzzz_0_xyyyyyy_1, g_zzzzzz_0_xyyyyyyy_1, g_zzzzzz_0_xyyyyyyz_1, g_zzzzzz_0_xyyyyyz_1, g_zzzzzz_0_xyyyyyzz_1, g_zzzzzz_0_xyyyyzz_1, g_zzzzzz_0_xyyyyzzz_1, g_zzzzzz_0_xyyyzzz_1, g_zzzzzz_0_xyyyzzzz_1, g_zzzzzz_0_xyyzzzz_1, g_zzzzzz_0_xyyzzzzz_1, g_zzzzzz_0_xyzzzzz_1, g_zzzzzz_0_xyzzzzzz_1, g_zzzzzz_0_xzzzzzz_1, g_zzzzzz_0_xzzzzzzz_1, g_zzzzzz_0_yyyyyyy_1, g_zzzzzz_0_yyyyyyyy_1, g_zzzzzz_0_yyyyyyyz_1, g_zzzzzz_0_yyyyyyz_1, g_zzzzzz_0_yyyyyyzz_1, g_zzzzzz_0_yyyyyzz_1, g_zzzzzz_0_yyyyyzzz_1, g_zzzzzz_0_yyyyzzz_1, g_zzzzzz_0_yyyyzzzz_1, g_zzzzzz_0_yyyzzzz_1, g_zzzzzz_0_yyyzzzzz_1, g_zzzzzz_0_yyzzzzz_1, g_zzzzzz_0_yyzzzzzz_1, g_zzzzzz_0_yzzzzzz_1, g_zzzzzz_0_yzzzzzzz_1, g_zzzzzz_0_zzzzzzz_1, g_zzzzzz_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzzz_0_xxxxxxxx_0[i] = 8.0 * g_zzzzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxxx_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxxxy_0[i] = 7.0 * g_zzzzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxxxz_0[i] = 7.0 * g_zzzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxxyy_0[i] = 6.0 * g_zzzzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxxyz_0[i] = 6.0 * g_zzzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxxzz_0[i] = 6.0 * g_zzzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxyyy_0[i] = 5.0 * g_zzzzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxyyz_0[i] = 5.0 * g_zzzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxyzz_0[i] = 5.0 * g_zzzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxzzz_0[i] = 5.0 * g_zzzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxyyyy_0[i] = 4.0 * g_zzzzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxyyyz_0[i] = 4.0 * g_zzzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxyyzz_0[i] = 4.0 * g_zzzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxyzzz_0[i] = 4.0 * g_zzzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxzzzz_0[i] = 4.0 * g_zzzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyyyyy_0[i] = 3.0 * g_zzzzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyyyyz_0[i] = 3.0 * g_zzzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyyyzz_0[i] = 3.0 * g_zzzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyyzzz_0[i] = 3.0 * g_zzzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyzzzz_0[i] = 3.0 * g_zzzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxzzzzz_0[i] = 3.0 * g_zzzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyyyyy_0[i] = 2.0 * g_zzzzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyyyyz_0[i] = 2.0 * g_zzzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyyyzz_0[i] = 2.0 * g_zzzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyyzzz_0[i] = 2.0 * g_zzzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyzzzz_0[i] = 2.0 * g_zzzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyzzzzz_0[i] = 2.0 * g_zzzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxzzzzzz_0[i] = 2.0 * g_zzzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyyyyy_0[i] = g_zzzzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyyyyz_0[i] = g_zzzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyyyzz_0[i] = g_zzzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyyzzz_0[i] = g_zzzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyzzzz_0[i] = g_zzzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyzzzzz_0[i] = g_zzzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyzzzzzz_0[i] = g_zzzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xzzzzzzz_0[i] = g_zzzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyyyyy_0[i] = g_zzzzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyyyyz_0[i] = g_zzzzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyyyzz_0[i] = g_zzzzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyyzzz_0[i] = g_zzzzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyzzzz_0[i] = g_zzzzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyzzzzz_0[i] = g_zzzzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyzzzzzz_0[i] = g_zzzzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yzzzzzzz_0[i] = g_zzzzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_zzzzzzzz_0[i] = g_zzzzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 1260-1305 components of targeted buffer : KSL

    auto g_yyyyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1260);

    auto g_yyyyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1261);

    auto g_yyyyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1262);

    auto g_yyyyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1263);

    auto g_yyyyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1264);

    auto g_yyyyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1265);

    auto g_yyyyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1266);

    auto g_yyyyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1267);

    auto g_yyyyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1268);

    auto g_yyyyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1269);

    auto g_yyyyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1270);

    auto g_yyyyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1271);

    auto g_yyyyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1272);

    auto g_yyyyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1273);

    auto g_yyyyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1274);

    auto g_yyyyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1275);

    auto g_yyyyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1276);

    auto g_yyyyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1277);

    auto g_yyyyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1278);

    auto g_yyyyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1279);

    auto g_yyyyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1280);

    auto g_yyyyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1281);

    auto g_yyyyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1282);

    auto g_yyyyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1283);

    auto g_yyyyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1284);

    auto g_yyyyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1285);

    auto g_yyyyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1286);

    auto g_yyyyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1287);

    auto g_yyyyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1288);

    auto g_yyyyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1289);

    auto g_yyyyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1290);

    auto g_yyyyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1291);

    auto g_yyyyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1292);

    auto g_yyyyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1293);

    auto g_yyyyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1294);

    auto g_yyyyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1295);

    auto g_yyyyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1296);

    auto g_yyyyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1297);

    auto g_yyyyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1298);

    auto g_yyyyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1299);

    auto g_yyyyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1300);

    auto g_yyyyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1301);

    auto g_yyyyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1302);

    auto g_yyyyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1303);

    auto g_yyyyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1304);

    #pragma omp simd aligned(g_yyyyy_0_xxxxxxxx_0, g_yyyyy_0_xxxxxxxx_1, g_yyyyy_0_xxxxxxxy_0, g_yyyyy_0_xxxxxxxy_1, g_yyyyy_0_xxxxxxxz_0, g_yyyyy_0_xxxxxxxz_1, g_yyyyy_0_xxxxxxyy_0, g_yyyyy_0_xxxxxxyy_1, g_yyyyy_0_xxxxxxyz_0, g_yyyyy_0_xxxxxxyz_1, g_yyyyy_0_xxxxxxzz_0, g_yyyyy_0_xxxxxxzz_1, g_yyyyy_0_xxxxxyyy_0, g_yyyyy_0_xxxxxyyy_1, g_yyyyy_0_xxxxxyyz_0, g_yyyyy_0_xxxxxyyz_1, g_yyyyy_0_xxxxxyzz_0, g_yyyyy_0_xxxxxyzz_1, g_yyyyy_0_xxxxxzzz_0, g_yyyyy_0_xxxxxzzz_1, g_yyyyy_0_xxxxyyyy_0, g_yyyyy_0_xxxxyyyy_1, g_yyyyy_0_xxxxyyyz_0, g_yyyyy_0_xxxxyyyz_1, g_yyyyy_0_xxxxyyzz_0, g_yyyyy_0_xxxxyyzz_1, g_yyyyy_0_xxxxyzzz_0, g_yyyyy_0_xxxxyzzz_1, g_yyyyy_0_xxxxzzzz_0, g_yyyyy_0_xxxxzzzz_1, g_yyyyy_0_xxxyyyyy_0, g_yyyyy_0_xxxyyyyy_1, g_yyyyy_0_xxxyyyyz_0, g_yyyyy_0_xxxyyyyz_1, g_yyyyy_0_xxxyyyzz_0, g_yyyyy_0_xxxyyyzz_1, g_yyyyy_0_xxxyyzzz_0, g_yyyyy_0_xxxyyzzz_1, g_yyyyy_0_xxxyzzzz_0, g_yyyyy_0_xxxyzzzz_1, g_yyyyy_0_xxxzzzzz_0, g_yyyyy_0_xxxzzzzz_1, g_yyyyy_0_xxyyyyyy_0, g_yyyyy_0_xxyyyyyy_1, g_yyyyy_0_xxyyyyyz_0, g_yyyyy_0_xxyyyyyz_1, g_yyyyy_0_xxyyyyzz_0, g_yyyyy_0_xxyyyyzz_1, g_yyyyy_0_xxyyyzzz_0, g_yyyyy_0_xxyyyzzz_1, g_yyyyy_0_xxyyzzzz_0, g_yyyyy_0_xxyyzzzz_1, g_yyyyy_0_xxyzzzzz_0, g_yyyyy_0_xxyzzzzz_1, g_yyyyy_0_xxzzzzzz_0, g_yyyyy_0_xxzzzzzz_1, g_yyyyy_0_xyyyyyyy_0, g_yyyyy_0_xyyyyyyy_1, g_yyyyy_0_xyyyyyyz_0, g_yyyyy_0_xyyyyyyz_1, g_yyyyy_0_xyyyyyzz_0, g_yyyyy_0_xyyyyyzz_1, g_yyyyy_0_xyyyyzzz_0, g_yyyyy_0_xyyyyzzz_1, g_yyyyy_0_xyyyzzzz_0, g_yyyyy_0_xyyyzzzz_1, g_yyyyy_0_xyyzzzzz_0, g_yyyyy_0_xyyzzzzz_1, g_yyyyy_0_xyzzzzzz_0, g_yyyyy_0_xyzzzzzz_1, g_yyyyy_0_xzzzzzzz_0, g_yyyyy_0_xzzzzzzz_1, g_yyyyy_0_yyyyyyyy_0, g_yyyyy_0_yyyyyyyy_1, g_yyyyy_0_yyyyyyyz_0, g_yyyyy_0_yyyyyyyz_1, g_yyyyy_0_yyyyyyzz_0, g_yyyyy_0_yyyyyyzz_1, g_yyyyy_0_yyyyyzzz_0, g_yyyyy_0_yyyyyzzz_1, g_yyyyy_0_yyyyzzzz_0, g_yyyyy_0_yyyyzzzz_1, g_yyyyy_0_yyyzzzzz_0, g_yyyyy_0_yyyzzzzz_1, g_yyyyy_0_yyzzzzzz_0, g_yyyyy_0_yyzzzzzz_1, g_yyyyy_0_yzzzzzzz_0, g_yyyyy_0_yzzzzzzz_1, g_yyyyy_0_zzzzzzzz_0, g_yyyyy_0_zzzzzzzz_1, g_yyyyyy_0_xxxxxxx_1, g_yyyyyy_0_xxxxxxxx_1, g_yyyyyy_0_xxxxxxxy_1, g_yyyyyy_0_xxxxxxxz_1, g_yyyyyy_0_xxxxxxy_1, g_yyyyyy_0_xxxxxxyy_1, g_yyyyyy_0_xxxxxxyz_1, g_yyyyyy_0_xxxxxxz_1, g_yyyyyy_0_xxxxxxzz_1, g_yyyyyy_0_xxxxxyy_1, g_yyyyyy_0_xxxxxyyy_1, g_yyyyyy_0_xxxxxyyz_1, g_yyyyyy_0_xxxxxyz_1, g_yyyyyy_0_xxxxxyzz_1, g_yyyyyy_0_xxxxxzz_1, g_yyyyyy_0_xxxxxzzz_1, g_yyyyyy_0_xxxxyyy_1, g_yyyyyy_0_xxxxyyyy_1, g_yyyyyy_0_xxxxyyyz_1, g_yyyyyy_0_xxxxyyz_1, g_yyyyyy_0_xxxxyyzz_1, g_yyyyyy_0_xxxxyzz_1, g_yyyyyy_0_xxxxyzzz_1, g_yyyyyy_0_xxxxzzz_1, g_yyyyyy_0_xxxxzzzz_1, g_yyyyyy_0_xxxyyyy_1, g_yyyyyy_0_xxxyyyyy_1, g_yyyyyy_0_xxxyyyyz_1, g_yyyyyy_0_xxxyyyz_1, g_yyyyyy_0_xxxyyyzz_1, g_yyyyyy_0_xxxyyzz_1, g_yyyyyy_0_xxxyyzzz_1, g_yyyyyy_0_xxxyzzz_1, g_yyyyyy_0_xxxyzzzz_1, g_yyyyyy_0_xxxzzzz_1, g_yyyyyy_0_xxxzzzzz_1, g_yyyyyy_0_xxyyyyy_1, g_yyyyyy_0_xxyyyyyy_1, g_yyyyyy_0_xxyyyyyz_1, g_yyyyyy_0_xxyyyyz_1, g_yyyyyy_0_xxyyyyzz_1, g_yyyyyy_0_xxyyyzz_1, g_yyyyyy_0_xxyyyzzz_1, g_yyyyyy_0_xxyyzzz_1, g_yyyyyy_0_xxyyzzzz_1, g_yyyyyy_0_xxyzzzz_1, g_yyyyyy_0_xxyzzzzz_1, g_yyyyyy_0_xxzzzzz_1, g_yyyyyy_0_xxzzzzzz_1, g_yyyyyy_0_xyyyyyy_1, g_yyyyyy_0_xyyyyyyy_1, g_yyyyyy_0_xyyyyyyz_1, g_yyyyyy_0_xyyyyyz_1, g_yyyyyy_0_xyyyyyzz_1, g_yyyyyy_0_xyyyyzz_1, g_yyyyyy_0_xyyyyzzz_1, g_yyyyyy_0_xyyyzzz_1, g_yyyyyy_0_xyyyzzzz_1, g_yyyyyy_0_xyyzzzz_1, g_yyyyyy_0_xyyzzzzz_1, g_yyyyyy_0_xyzzzzz_1, g_yyyyyy_0_xyzzzzzz_1, g_yyyyyy_0_xzzzzzz_1, g_yyyyyy_0_xzzzzzzz_1, g_yyyyyy_0_yyyyyyy_1, g_yyyyyy_0_yyyyyyyy_1, g_yyyyyy_0_yyyyyyyz_1, g_yyyyyy_0_yyyyyyz_1, g_yyyyyy_0_yyyyyyzz_1, g_yyyyyy_0_yyyyyzz_1, g_yyyyyy_0_yyyyyzzz_1, g_yyyyyy_0_yyyyzzz_1, g_yyyyyy_0_yyyyzzzz_1, g_yyyyyy_0_yyyzzzz_1, g_yyyyyy_0_yyyzzzzz_1, g_yyyyyy_0_yyzzzzz_1, g_yyyyyy_0_yyzzzzzz_1, g_yyyyyy_0_yzzzzzz_1, g_yyyyyy_0_yzzzzzzz_1, g_yyyyyy_0_zzzzzzz_1, g_yyyyyy_0_zzzzzzzz_1, g_yyyyyyy_0_xxxxxxxx_0, g_yyyyyyy_0_xxxxxxxy_0, g_yyyyyyy_0_xxxxxxxz_0, g_yyyyyyy_0_xxxxxxyy_0, g_yyyyyyy_0_xxxxxxyz_0, g_yyyyyyy_0_xxxxxxzz_0, g_yyyyyyy_0_xxxxxyyy_0, g_yyyyyyy_0_xxxxxyyz_0, g_yyyyyyy_0_xxxxxyzz_0, g_yyyyyyy_0_xxxxxzzz_0, g_yyyyyyy_0_xxxxyyyy_0, g_yyyyyyy_0_xxxxyyyz_0, g_yyyyyyy_0_xxxxyyzz_0, g_yyyyyyy_0_xxxxyzzz_0, g_yyyyyyy_0_xxxxzzzz_0, g_yyyyyyy_0_xxxyyyyy_0, g_yyyyyyy_0_xxxyyyyz_0, g_yyyyyyy_0_xxxyyyzz_0, g_yyyyyyy_0_xxxyyzzz_0, g_yyyyyyy_0_xxxyzzzz_0, g_yyyyyyy_0_xxxzzzzz_0, g_yyyyyyy_0_xxyyyyyy_0, g_yyyyyyy_0_xxyyyyyz_0, g_yyyyyyy_0_xxyyyyzz_0, g_yyyyyyy_0_xxyyyzzz_0, g_yyyyyyy_0_xxyyzzzz_0, g_yyyyyyy_0_xxyzzzzz_0, g_yyyyyyy_0_xxzzzzzz_0, g_yyyyyyy_0_xyyyyyyy_0, g_yyyyyyy_0_xyyyyyyz_0, g_yyyyyyy_0_xyyyyyzz_0, g_yyyyyyy_0_xyyyyzzz_0, g_yyyyyyy_0_xyyyzzzz_0, g_yyyyyyy_0_xyyzzzzz_0, g_yyyyyyy_0_xyzzzzzz_0, g_yyyyyyy_0_xzzzzzzz_0, g_yyyyyyy_0_yyyyyyyy_0, g_yyyyyyy_0_yyyyyyyz_0, g_yyyyyyy_0_yyyyyyzz_0, g_yyyyyyy_0_yyyyyzzz_0, g_yyyyyyy_0_yyyyzzzz_0, g_yyyyyyy_0_yyyzzzzz_0, g_yyyyyyy_0_yyzzzzzz_0, g_yyyyyyy_0_yzzzzzzz_0, g_yyyyyyy_0_zzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyyy_0_xxxxxxxx_0[i] = 6.0 * g_yyyyy_0_xxxxxxxx_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxxxx_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxxxy_0[i] = 6.0 * g_yyyyy_0_xxxxxxxy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxxxy_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxxy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxxxz_0[i] = 6.0 * g_yyyyy_0_xxxxxxxz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxxxz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxxyy_0[i] = 6.0 * g_yyyyy_0_xxxxxxyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxxyy_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxxyz_0[i] = 6.0 * g_yyyyy_0_xxxxxxyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxxyz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxxzz_0[i] = 6.0 * g_yyyyy_0_xxxxxxzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxxzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxyyy_0[i] = 6.0 * g_yyyyy_0_xxxxxyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxyyy_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxyyz_0[i] = 6.0 * g_yyyyy_0_xxxxxyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxyzz_0[i] = 6.0 * g_yyyyy_0_xxxxxyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxyzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxzzz_0[i] = 6.0 * g_yyyyy_0_xxxxxzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxyyyy_0[i] = 6.0 * g_yyyyy_0_xxxxyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxyyyz_0[i] = 6.0 * g_yyyyy_0_xxxxyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxyyzz_0[i] = 6.0 * g_yyyyy_0_xxxxyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxyzzz_0[i] = 6.0 * g_yyyyy_0_xxxxyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxyzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxzzzz_0[i] = 6.0 * g_yyyyy_0_xxxxzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyyyyy_0[i] = 6.0 * g_yyyyy_0_xxxyyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyyyyy_1[i] * fz_be_0 + 5.0 * g_yyyyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyyyyz_0[i] = 6.0 * g_yyyyy_0_xxxyyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyyyzz_0[i] = 6.0 * g_yyyyy_0_xxxyyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyyzzz_0[i] = 6.0 * g_yyyyy_0_xxxyyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyzzzz_0[i] = 6.0 * g_yyyyy_0_xxxyzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxzzzzz_0[i] = 6.0 * g_yyyyy_0_xxxzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyyyyy_0[i] = 6.0 * g_yyyyy_0_xxyyyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyyyyy_1[i] * fz_be_0 + 6.0 * g_yyyyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyyyyz_0[i] = 6.0 * g_yyyyy_0_xxyyyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyyyzz_0[i] = 6.0 * g_yyyyy_0_xxyyyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyyzzz_0[i] = 6.0 * g_yyyyy_0_xxyyyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyzzzz_0[i] = 6.0 * g_yyyyy_0_xxyyzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyzzzzz_0[i] = 6.0 * g_yyyyy_0_xxyzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxzzzzzz_0[i] = 6.0 * g_yyyyy_0_xxzzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxzzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyyyyy_0[i] = 6.0 * g_yyyyy_0_xyyyyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyyyyy_1[i] * fz_be_0 + 7.0 * g_yyyyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyyyyz_0[i] = 6.0 * g_yyyyy_0_xyyyyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyyyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyyyzz_0[i] = 6.0 * g_yyyyy_0_xyyyyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyyyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyyzzz_0[i] = 6.0 * g_yyyyy_0_xyyyyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyzzzz_0[i] = 6.0 * g_yyyyy_0_xyyyzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyzzzzz_0[i] = 6.0 * g_yyyyy_0_xyyzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyzzzzzz_0[i] = 6.0 * g_yyyyy_0_xyzzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyzzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xzzzzzzz_0[i] = 6.0 * g_yyyyy_0_xzzzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xzzzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyyyyy_0[i] = 6.0 * g_yyyyy_0_yyyyyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyyyyy_1[i] * fz_be_0 + 8.0 * g_yyyyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyyyyz_0[i] = 6.0 * g_yyyyy_0_yyyyyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yyyyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyyyzz_0[i] = 6.0 * g_yyyyy_0_yyyyyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yyyyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyyzzz_0[i] = 6.0 * g_yyyyy_0_yyyyyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yyyyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyzzzz_0[i] = 6.0 * g_yyyyy_0_yyyyzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyzzzzz_0[i] = 6.0 * g_yyyyy_0_yyyzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyzzzzzz_0[i] = 6.0 * g_yyyyy_0_yyzzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yzzzzzzz_0[i] = 6.0 * g_yyyyy_0_yzzzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yzzzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_zzzzzzzz_0[i] = 6.0 * g_yyyyy_0_zzzzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_zzzzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1305-1350 components of targeted buffer : KSL

    auto g_yyyyyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1305);

    auto g_yyyyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1306);

    auto g_yyyyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1307);

    auto g_yyyyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1308);

    auto g_yyyyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1309);

    auto g_yyyyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1310);

    auto g_yyyyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1311);

    auto g_yyyyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1312);

    auto g_yyyyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1313);

    auto g_yyyyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1314);

    auto g_yyyyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1315);

    auto g_yyyyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1316);

    auto g_yyyyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1317);

    auto g_yyyyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1318);

    auto g_yyyyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1319);

    auto g_yyyyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1320);

    auto g_yyyyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1321);

    auto g_yyyyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1322);

    auto g_yyyyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1323);

    auto g_yyyyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1324);

    auto g_yyyyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1325);

    auto g_yyyyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1326);

    auto g_yyyyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1327);

    auto g_yyyyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1328);

    auto g_yyyyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1329);

    auto g_yyyyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1330);

    auto g_yyyyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1331);

    auto g_yyyyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1332);

    auto g_yyyyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1333);

    auto g_yyyyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1334);

    auto g_yyyyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1335);

    auto g_yyyyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1336);

    auto g_yyyyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1337);

    auto g_yyyyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1338);

    auto g_yyyyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1339);

    auto g_yyyyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1340);

    auto g_yyyyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1341);

    auto g_yyyyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1342);

    auto g_yyyyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1343);

    auto g_yyyyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1344);

    auto g_yyyyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1345);

    auto g_yyyyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1346);

    auto g_yyyyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1347);

    auto g_yyyyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1348);

    auto g_yyyyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1349);

    #pragma omp simd aligned(g_yyyyyy_0_xxxxxxx_1, g_yyyyyy_0_xxxxxxxx_1, g_yyyyyy_0_xxxxxxxy_1, g_yyyyyy_0_xxxxxxxz_1, g_yyyyyy_0_xxxxxxy_1, g_yyyyyy_0_xxxxxxyy_1, g_yyyyyy_0_xxxxxxyz_1, g_yyyyyy_0_xxxxxxz_1, g_yyyyyy_0_xxxxxxzz_1, g_yyyyyy_0_xxxxxyy_1, g_yyyyyy_0_xxxxxyyy_1, g_yyyyyy_0_xxxxxyyz_1, g_yyyyyy_0_xxxxxyz_1, g_yyyyyy_0_xxxxxyzz_1, g_yyyyyy_0_xxxxxzz_1, g_yyyyyy_0_xxxxxzzz_1, g_yyyyyy_0_xxxxyyy_1, g_yyyyyy_0_xxxxyyyy_1, g_yyyyyy_0_xxxxyyyz_1, g_yyyyyy_0_xxxxyyz_1, g_yyyyyy_0_xxxxyyzz_1, g_yyyyyy_0_xxxxyzz_1, g_yyyyyy_0_xxxxyzzz_1, g_yyyyyy_0_xxxxzzz_1, g_yyyyyy_0_xxxxzzzz_1, g_yyyyyy_0_xxxyyyy_1, g_yyyyyy_0_xxxyyyyy_1, g_yyyyyy_0_xxxyyyyz_1, g_yyyyyy_0_xxxyyyz_1, g_yyyyyy_0_xxxyyyzz_1, g_yyyyyy_0_xxxyyzz_1, g_yyyyyy_0_xxxyyzzz_1, g_yyyyyy_0_xxxyzzz_1, g_yyyyyy_0_xxxyzzzz_1, g_yyyyyy_0_xxxzzzz_1, g_yyyyyy_0_xxxzzzzz_1, g_yyyyyy_0_xxyyyyy_1, g_yyyyyy_0_xxyyyyyy_1, g_yyyyyy_0_xxyyyyyz_1, g_yyyyyy_0_xxyyyyz_1, g_yyyyyy_0_xxyyyyzz_1, g_yyyyyy_0_xxyyyzz_1, g_yyyyyy_0_xxyyyzzz_1, g_yyyyyy_0_xxyyzzz_1, g_yyyyyy_0_xxyyzzzz_1, g_yyyyyy_0_xxyzzzz_1, g_yyyyyy_0_xxyzzzzz_1, g_yyyyyy_0_xxzzzzz_1, g_yyyyyy_0_xxzzzzzz_1, g_yyyyyy_0_xyyyyyy_1, g_yyyyyy_0_xyyyyyyy_1, g_yyyyyy_0_xyyyyyyz_1, g_yyyyyy_0_xyyyyyz_1, g_yyyyyy_0_xyyyyyzz_1, g_yyyyyy_0_xyyyyzz_1, g_yyyyyy_0_xyyyyzzz_1, g_yyyyyy_0_xyyyzzz_1, g_yyyyyy_0_xyyyzzzz_1, g_yyyyyy_0_xyyzzzz_1, g_yyyyyy_0_xyyzzzzz_1, g_yyyyyy_0_xyzzzzz_1, g_yyyyyy_0_xyzzzzzz_1, g_yyyyyy_0_xzzzzzz_1, g_yyyyyy_0_xzzzzzzz_1, g_yyyyyy_0_yyyyyyy_1, g_yyyyyy_0_yyyyyyyy_1, g_yyyyyy_0_yyyyyyyz_1, g_yyyyyy_0_yyyyyyz_1, g_yyyyyy_0_yyyyyyzz_1, g_yyyyyy_0_yyyyyzz_1, g_yyyyyy_0_yyyyyzzz_1, g_yyyyyy_0_yyyyzzz_1, g_yyyyyy_0_yyyyzzzz_1, g_yyyyyy_0_yyyzzzz_1, g_yyyyyy_0_yyyzzzzz_1, g_yyyyyy_0_yyzzzzz_1, g_yyyyyy_0_yyzzzzzz_1, g_yyyyyy_0_yzzzzzz_1, g_yyyyyy_0_yzzzzzzz_1, g_yyyyyy_0_zzzzzzz_1, g_yyyyyy_0_zzzzzzzz_1, g_yyyyyyz_0_xxxxxxxx_0, g_yyyyyyz_0_xxxxxxxy_0, g_yyyyyyz_0_xxxxxxxz_0, g_yyyyyyz_0_xxxxxxyy_0, g_yyyyyyz_0_xxxxxxyz_0, g_yyyyyyz_0_xxxxxxzz_0, g_yyyyyyz_0_xxxxxyyy_0, g_yyyyyyz_0_xxxxxyyz_0, g_yyyyyyz_0_xxxxxyzz_0, g_yyyyyyz_0_xxxxxzzz_0, g_yyyyyyz_0_xxxxyyyy_0, g_yyyyyyz_0_xxxxyyyz_0, g_yyyyyyz_0_xxxxyyzz_0, g_yyyyyyz_0_xxxxyzzz_0, g_yyyyyyz_0_xxxxzzzz_0, g_yyyyyyz_0_xxxyyyyy_0, g_yyyyyyz_0_xxxyyyyz_0, g_yyyyyyz_0_xxxyyyzz_0, g_yyyyyyz_0_xxxyyzzz_0, g_yyyyyyz_0_xxxyzzzz_0, g_yyyyyyz_0_xxxzzzzz_0, g_yyyyyyz_0_xxyyyyyy_0, g_yyyyyyz_0_xxyyyyyz_0, g_yyyyyyz_0_xxyyyyzz_0, g_yyyyyyz_0_xxyyyzzz_0, g_yyyyyyz_0_xxyyzzzz_0, g_yyyyyyz_0_xxyzzzzz_0, g_yyyyyyz_0_xxzzzzzz_0, g_yyyyyyz_0_xyyyyyyy_0, g_yyyyyyz_0_xyyyyyyz_0, g_yyyyyyz_0_xyyyyyzz_0, g_yyyyyyz_0_xyyyyzzz_0, g_yyyyyyz_0_xyyyzzzz_0, g_yyyyyyz_0_xyyzzzzz_0, g_yyyyyyz_0_xyzzzzzz_0, g_yyyyyyz_0_xzzzzzzz_0, g_yyyyyyz_0_yyyyyyyy_0, g_yyyyyyz_0_yyyyyyyz_0, g_yyyyyyz_0_yyyyyyzz_0, g_yyyyyyz_0_yyyyyzzz_0, g_yyyyyyz_0_yyyyzzzz_0, g_yyyyyyz_0_yyyzzzzz_0, g_yyyyyyz_0_yyzzzzzz_0, g_yyyyyyz_0_yzzzzzzz_0, g_yyyyyyz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyyz_0_xxxxxxxx_0[i] = g_yyyyyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxxxy_0[i] = g_yyyyyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxxxz_0[i] = g_yyyyyy_0_xxxxxxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxxz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxxyy_0[i] = g_yyyyyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxxyz_0[i] = g_yyyyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxxzz_0[i] = 2.0 * g_yyyyyy_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxyyy_0[i] = g_yyyyyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxyyz_0[i] = g_yyyyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxyzz_0[i] = 2.0 * g_yyyyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxzzz_0[i] = 3.0 * g_yyyyyy_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxyyyy_0[i] = g_yyyyyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxyyyz_0[i] = g_yyyyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxyyzz_0[i] = 2.0 * g_yyyyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxyzzz_0[i] = 3.0 * g_yyyyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxzzzz_0[i] = 4.0 * g_yyyyyy_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyyyyy_0[i] = g_yyyyyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyyyyz_0[i] = g_yyyyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyyyzz_0[i] = 2.0 * g_yyyyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyyzzz_0[i] = 3.0 * g_yyyyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyzzzz_0[i] = 4.0 * g_yyyyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxzzzzz_0[i] = 5.0 * g_yyyyyy_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyyyyy_0[i] = g_yyyyyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyyyyz_0[i] = g_yyyyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyyyzz_0[i] = 2.0 * g_yyyyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyyzzz_0[i] = 3.0 * g_yyyyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyzzzz_0[i] = 4.0 * g_yyyyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyzzzzz_0[i] = 5.0 * g_yyyyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxzzzzzz_0[i] = 6.0 * g_yyyyyy_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxzzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyyyyy_0[i] = g_yyyyyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyyyyz_0[i] = g_yyyyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyyyzz_0[i] = 2.0 * g_yyyyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyyzzz_0[i] = 3.0 * g_yyyyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyzzzz_0[i] = 4.0 * g_yyyyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyzzzzz_0[i] = 5.0 * g_yyyyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyzzzzzz_0[i] = 6.0 * g_yyyyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xzzzzzzz_0[i] = 7.0 * g_yyyyyy_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xzzzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyyyyy_0[i] = g_yyyyyy_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyyyyz_0[i] = g_yyyyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyyyzz_0[i] = 2.0 * g_yyyyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyyzzz_0[i] = 3.0 * g_yyyyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyzzzz_0[i] = 4.0 * g_yyyyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyzzzzz_0[i] = 5.0 * g_yyyyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyzzzzzz_0[i] = 6.0 * g_yyyyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyzzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yzzzzzzz_0[i] = 7.0 * g_yyyyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yzzzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_zzzzzzzz_0[i] = 8.0 * g_yyyyyy_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 1350-1395 components of targeted buffer : KSL

    auto g_yyyyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1350);

    auto g_yyyyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1351);

    auto g_yyyyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1352);

    auto g_yyyyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1353);

    auto g_yyyyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1354);

    auto g_yyyyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1355);

    auto g_yyyyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1356);

    auto g_yyyyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1357);

    auto g_yyyyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1358);

    auto g_yyyyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1359);

    auto g_yyyyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1360);

    auto g_yyyyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1361);

    auto g_yyyyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1362);

    auto g_yyyyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1363);

    auto g_yyyyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1364);

    auto g_yyyyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1365);

    auto g_yyyyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1366);

    auto g_yyyyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1367);

    auto g_yyyyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1368);

    auto g_yyyyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1369);

    auto g_yyyyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1370);

    auto g_yyyyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1371);

    auto g_yyyyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1372);

    auto g_yyyyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1373);

    auto g_yyyyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1374);

    auto g_yyyyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1375);

    auto g_yyyyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1376);

    auto g_yyyyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1377);

    auto g_yyyyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1378);

    auto g_yyyyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1379);

    auto g_yyyyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1380);

    auto g_yyyyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1381);

    auto g_yyyyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1382);

    auto g_yyyyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1383);

    auto g_yyyyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1384);

    auto g_yyyyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1385);

    auto g_yyyyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1386);

    auto g_yyyyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1387);

    auto g_yyyyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1388);

    auto g_yyyyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1389);

    auto g_yyyyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1390);

    auto g_yyyyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1391);

    auto g_yyyyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1392);

    auto g_yyyyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1393);

    auto g_yyyyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1394);

    #pragma omp simd aligned(g_yyyyy_0_xxxxxxxy_0, g_yyyyy_0_xxxxxxxy_1, g_yyyyy_0_xxxxxxyy_0, g_yyyyy_0_xxxxxxyy_1, g_yyyyy_0_xxxxxyyy_0, g_yyyyy_0_xxxxxyyy_1, g_yyyyy_0_xxxxyyyy_0, g_yyyyy_0_xxxxyyyy_1, g_yyyyy_0_xxxyyyyy_0, g_yyyyy_0_xxxyyyyy_1, g_yyyyy_0_xxyyyyyy_0, g_yyyyy_0_xxyyyyyy_1, g_yyyyy_0_xyyyyyyy_0, g_yyyyy_0_xyyyyyyy_1, g_yyyyy_0_yyyyyyyy_0, g_yyyyy_0_yyyyyyyy_1, g_yyyyyz_0_xxxxxxxy_1, g_yyyyyz_0_xxxxxxyy_1, g_yyyyyz_0_xxxxxyyy_1, g_yyyyyz_0_xxxxyyyy_1, g_yyyyyz_0_xxxyyyyy_1, g_yyyyyz_0_xxyyyyyy_1, g_yyyyyz_0_xyyyyyyy_1, g_yyyyyz_0_yyyyyyyy_1, g_yyyyyzz_0_xxxxxxxx_0, g_yyyyyzz_0_xxxxxxxy_0, g_yyyyyzz_0_xxxxxxxz_0, g_yyyyyzz_0_xxxxxxyy_0, g_yyyyyzz_0_xxxxxxyz_0, g_yyyyyzz_0_xxxxxxzz_0, g_yyyyyzz_0_xxxxxyyy_0, g_yyyyyzz_0_xxxxxyyz_0, g_yyyyyzz_0_xxxxxyzz_0, g_yyyyyzz_0_xxxxxzzz_0, g_yyyyyzz_0_xxxxyyyy_0, g_yyyyyzz_0_xxxxyyyz_0, g_yyyyyzz_0_xxxxyyzz_0, g_yyyyyzz_0_xxxxyzzz_0, g_yyyyyzz_0_xxxxzzzz_0, g_yyyyyzz_0_xxxyyyyy_0, g_yyyyyzz_0_xxxyyyyz_0, g_yyyyyzz_0_xxxyyyzz_0, g_yyyyyzz_0_xxxyyzzz_0, g_yyyyyzz_0_xxxyzzzz_0, g_yyyyyzz_0_xxxzzzzz_0, g_yyyyyzz_0_xxyyyyyy_0, g_yyyyyzz_0_xxyyyyyz_0, g_yyyyyzz_0_xxyyyyzz_0, g_yyyyyzz_0_xxyyyzzz_0, g_yyyyyzz_0_xxyyzzzz_0, g_yyyyyzz_0_xxyzzzzz_0, g_yyyyyzz_0_xxzzzzzz_0, g_yyyyyzz_0_xyyyyyyy_0, g_yyyyyzz_0_xyyyyyyz_0, g_yyyyyzz_0_xyyyyyzz_0, g_yyyyyzz_0_xyyyyzzz_0, g_yyyyyzz_0_xyyyzzzz_0, g_yyyyyzz_0_xyyzzzzz_0, g_yyyyyzz_0_xyzzzzzz_0, g_yyyyyzz_0_xzzzzzzz_0, g_yyyyyzz_0_yyyyyyyy_0, g_yyyyyzz_0_yyyyyyyz_0, g_yyyyyzz_0_yyyyyyzz_0, g_yyyyyzz_0_yyyyyzzz_0, g_yyyyyzz_0_yyyyzzzz_0, g_yyyyyzz_0_yyyzzzzz_0, g_yyyyyzz_0_yyzzzzzz_0, g_yyyyyzz_0_yzzzzzzz_0, g_yyyyyzz_0_zzzzzzzz_0, g_yyyyzz_0_xxxxxxxx_1, g_yyyyzz_0_xxxxxxxz_1, g_yyyyzz_0_xxxxxxyz_1, g_yyyyzz_0_xxxxxxz_1, g_yyyyzz_0_xxxxxxzz_1, g_yyyyzz_0_xxxxxyyz_1, g_yyyyzz_0_xxxxxyz_1, g_yyyyzz_0_xxxxxyzz_1, g_yyyyzz_0_xxxxxzz_1, g_yyyyzz_0_xxxxxzzz_1, g_yyyyzz_0_xxxxyyyz_1, g_yyyyzz_0_xxxxyyz_1, g_yyyyzz_0_xxxxyyzz_1, g_yyyyzz_0_xxxxyzz_1, g_yyyyzz_0_xxxxyzzz_1, g_yyyyzz_0_xxxxzzz_1, g_yyyyzz_0_xxxxzzzz_1, g_yyyyzz_0_xxxyyyyz_1, g_yyyyzz_0_xxxyyyz_1, g_yyyyzz_0_xxxyyyzz_1, g_yyyyzz_0_xxxyyzz_1, g_yyyyzz_0_xxxyyzzz_1, g_yyyyzz_0_xxxyzzz_1, g_yyyyzz_0_xxxyzzzz_1, g_yyyyzz_0_xxxzzzz_1, g_yyyyzz_0_xxxzzzzz_1, g_yyyyzz_0_xxyyyyyz_1, g_yyyyzz_0_xxyyyyz_1, g_yyyyzz_0_xxyyyyzz_1, g_yyyyzz_0_xxyyyzz_1, g_yyyyzz_0_xxyyyzzz_1, g_yyyyzz_0_xxyyzzz_1, g_yyyyzz_0_xxyyzzzz_1, g_yyyyzz_0_xxyzzzz_1, g_yyyyzz_0_xxyzzzzz_1, g_yyyyzz_0_xxzzzzz_1, g_yyyyzz_0_xxzzzzzz_1, g_yyyyzz_0_xyyyyyyz_1, g_yyyyzz_0_xyyyyyz_1, g_yyyyzz_0_xyyyyyzz_1, g_yyyyzz_0_xyyyyzz_1, g_yyyyzz_0_xyyyyzzz_1, g_yyyyzz_0_xyyyzzz_1, g_yyyyzz_0_xyyyzzzz_1, g_yyyyzz_0_xyyzzzz_1, g_yyyyzz_0_xyyzzzzz_1, g_yyyyzz_0_xyzzzzz_1, g_yyyyzz_0_xyzzzzzz_1, g_yyyyzz_0_xzzzzzz_1, g_yyyyzz_0_xzzzzzzz_1, g_yyyyzz_0_yyyyyyyz_1, g_yyyyzz_0_yyyyyyz_1, g_yyyyzz_0_yyyyyyzz_1, g_yyyyzz_0_yyyyyzz_1, g_yyyyzz_0_yyyyyzzz_1, g_yyyyzz_0_yyyyzzz_1, g_yyyyzz_0_yyyyzzzz_1, g_yyyyzz_0_yyyzzzz_1, g_yyyyzz_0_yyyzzzzz_1, g_yyyyzz_0_yyzzzzz_1, g_yyyyzz_0_yyzzzzzz_1, g_yyyyzz_0_yzzzzzz_1, g_yyyyzz_0_yzzzzzzz_1, g_yyyyzz_0_zzzzzzz_1, g_yyyyzz_0_zzzzzzzz_1, g_yyyzz_0_xxxxxxxx_0, g_yyyzz_0_xxxxxxxx_1, g_yyyzz_0_xxxxxxxz_0, g_yyyzz_0_xxxxxxxz_1, g_yyyzz_0_xxxxxxyz_0, g_yyyzz_0_xxxxxxyz_1, g_yyyzz_0_xxxxxxzz_0, g_yyyzz_0_xxxxxxzz_1, g_yyyzz_0_xxxxxyyz_0, g_yyyzz_0_xxxxxyyz_1, g_yyyzz_0_xxxxxyzz_0, g_yyyzz_0_xxxxxyzz_1, g_yyyzz_0_xxxxxzzz_0, g_yyyzz_0_xxxxxzzz_1, g_yyyzz_0_xxxxyyyz_0, g_yyyzz_0_xxxxyyyz_1, g_yyyzz_0_xxxxyyzz_0, g_yyyzz_0_xxxxyyzz_1, g_yyyzz_0_xxxxyzzz_0, g_yyyzz_0_xxxxyzzz_1, g_yyyzz_0_xxxxzzzz_0, g_yyyzz_0_xxxxzzzz_1, g_yyyzz_0_xxxyyyyz_0, g_yyyzz_0_xxxyyyyz_1, g_yyyzz_0_xxxyyyzz_0, g_yyyzz_0_xxxyyyzz_1, g_yyyzz_0_xxxyyzzz_0, g_yyyzz_0_xxxyyzzz_1, g_yyyzz_0_xxxyzzzz_0, g_yyyzz_0_xxxyzzzz_1, g_yyyzz_0_xxxzzzzz_0, g_yyyzz_0_xxxzzzzz_1, g_yyyzz_0_xxyyyyyz_0, g_yyyzz_0_xxyyyyyz_1, g_yyyzz_0_xxyyyyzz_0, g_yyyzz_0_xxyyyyzz_1, g_yyyzz_0_xxyyyzzz_0, g_yyyzz_0_xxyyyzzz_1, g_yyyzz_0_xxyyzzzz_0, g_yyyzz_0_xxyyzzzz_1, g_yyyzz_0_xxyzzzzz_0, g_yyyzz_0_xxyzzzzz_1, g_yyyzz_0_xxzzzzzz_0, g_yyyzz_0_xxzzzzzz_1, g_yyyzz_0_xyyyyyyz_0, g_yyyzz_0_xyyyyyyz_1, g_yyyzz_0_xyyyyyzz_0, g_yyyzz_0_xyyyyyzz_1, g_yyyzz_0_xyyyyzzz_0, g_yyyzz_0_xyyyyzzz_1, g_yyyzz_0_xyyyzzzz_0, g_yyyzz_0_xyyyzzzz_1, g_yyyzz_0_xyyzzzzz_0, g_yyyzz_0_xyyzzzzz_1, g_yyyzz_0_xyzzzzzz_0, g_yyyzz_0_xyzzzzzz_1, g_yyyzz_0_xzzzzzzz_0, g_yyyzz_0_xzzzzzzz_1, g_yyyzz_0_yyyyyyyz_0, g_yyyzz_0_yyyyyyyz_1, g_yyyzz_0_yyyyyyzz_0, g_yyyzz_0_yyyyyyzz_1, g_yyyzz_0_yyyyyzzz_0, g_yyyzz_0_yyyyyzzz_1, g_yyyzz_0_yyyyzzzz_0, g_yyyzz_0_yyyyzzzz_1, g_yyyzz_0_yyyzzzzz_0, g_yyyzz_0_yyyzzzzz_1, g_yyyzz_0_yyzzzzzz_0, g_yyyzz_0_yyzzzzzz_1, g_yyyzz_0_yzzzzzzz_0, g_yyyzz_0_yzzzzzzz_1, g_yyyzz_0_zzzzzzzz_0, g_yyyzz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyzz_0_xxxxxxxx_0[i] = 4.0 * g_yyyzz_0_xxxxxxxx_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxxxx_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxxxxy_0[i] = g_yyyyy_0_xxxxxxxy_0[i] * fbe_0 - g_yyyyy_0_xxxxxxxy_1[i] * fz_be_0 + g_yyyyyz_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxxxxxz_0[i] = 4.0 * g_yyyzz_0_xxxxxxxz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxxxz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxxxyy_0[i] = g_yyyyy_0_xxxxxxyy_0[i] * fbe_0 - g_yyyyy_0_xxxxxxyy_1[i] * fz_be_0 + g_yyyyyz_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxxxxyz_0[i] = 4.0 * g_yyyzz_0_xxxxxxyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxxyz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxxxzz_0[i] = 4.0 * g_yyyzz_0_xxxxxxzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxxzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxxyyy_0[i] = g_yyyyy_0_xxxxxyyy_0[i] * fbe_0 - g_yyyyy_0_xxxxxyyy_1[i] * fz_be_0 + g_yyyyyz_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxxxyyz_0[i] = 4.0 * g_yyyzz_0_xxxxxyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxxyzz_0[i] = 4.0 * g_yyyzz_0_xxxxxyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxyzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxxzzz_0[i] = 4.0 * g_yyyzz_0_xxxxxzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxyyyy_0[i] = g_yyyyy_0_xxxxyyyy_0[i] * fbe_0 - g_yyyyy_0_xxxxyyyy_1[i] * fz_be_0 + g_yyyyyz_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxxyyyz_0[i] = 4.0 * g_yyyzz_0_xxxxyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxyyzz_0[i] = 4.0 * g_yyyzz_0_xxxxyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxyzzz_0[i] = 4.0 * g_yyyzz_0_xxxxyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxyzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxzzzz_0[i] = 4.0 * g_yyyzz_0_xxxxzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxyyyyy_0[i] = g_yyyyy_0_xxxyyyyy_0[i] * fbe_0 - g_yyyyy_0_xxxyyyyy_1[i] * fz_be_0 + g_yyyyyz_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxyyyyz_0[i] = 4.0 * g_yyyzz_0_xxxyyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxyyyzz_0[i] = 4.0 * g_yyyzz_0_xxxyyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxyyzzz_0[i] = 4.0 * g_yyyzz_0_xxxyyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxyzzzz_0[i] = 4.0 * g_yyyzz_0_xxxyzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxyzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxzzzzz_0[i] = 4.0 * g_yyyzz_0_xxxzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyyyyyy_0[i] = g_yyyyy_0_xxyyyyyy_0[i] * fbe_0 - g_yyyyy_0_xxyyyyyy_1[i] * fz_be_0 + g_yyyyyz_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxyyyyyz_0[i] = 4.0 * g_yyyzz_0_xxyyyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyyzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyyyyzz_0[i] = 4.0 * g_yyyzz_0_xxyyyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyyyzzz_0[i] = 4.0 * g_yyyzz_0_xxyyyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyyzzzz_0[i] = 4.0 * g_yyyzz_0_xxyyzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyzzzzz_0[i] = 4.0 * g_yyyzz_0_xxyzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxzzzzzz_0[i] = 4.0 * g_yyyzz_0_xxzzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxzzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyyyyyy_0[i] = g_yyyyy_0_xyyyyyyy_0[i] * fbe_0 - g_yyyyy_0_xyyyyyyy_1[i] * fz_be_0 + g_yyyyyz_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xyyyyyyz_0[i] = 4.0 * g_yyyzz_0_xyyyyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyyyzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyyyyzz_0[i] = 4.0 * g_yyyzz_0_xyyyyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyyyzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyyyzzz_0[i] = 4.0 * g_yyyzz_0_xyyyyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyyzzzz_0[i] = 4.0 * g_yyyzz_0_xyyyzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyzzzzz_0[i] = 4.0 * g_yyyzz_0_xyyzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyzzzzzz_0[i] = 4.0 * g_yyyzz_0_xyzzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyzzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xzzzzzzz_0[i] = 4.0 * g_yyyzz_0_xzzzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xzzzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyyyyyy_0[i] = g_yyyyy_0_yyyyyyyy_0[i] * fbe_0 - g_yyyyy_0_yyyyyyyy_1[i] * fz_be_0 + g_yyyyyz_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_yyyyyyyz_0[i] = 4.0 * g_yyyzz_0_yyyyyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yyyyzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyyyyzz_0[i] = 4.0 * g_yyyzz_0_yyyyyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yyyyzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyyyzzz_0[i] = 4.0 * g_yyyzz_0_yyyyyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yyyyzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyyzzzz_0[i] = 4.0 * g_yyyzz_0_yyyyzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyzzzzz_0[i] = 4.0 * g_yyyzz_0_yyyzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyzzzzzz_0[i] = 4.0 * g_yyyzz_0_yyzzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yzzzzzzz_0[i] = 4.0 * g_yyyzz_0_yzzzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yzzzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_zzzzzzzz_0[i] = 4.0 * g_yyyzz_0_zzzzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_zzzzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1395-1440 components of targeted buffer : KSL

    auto g_yyyyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1395);

    auto g_yyyyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1396);

    auto g_yyyyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1397);

    auto g_yyyyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1398);

    auto g_yyyyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1399);

    auto g_yyyyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1400);

    auto g_yyyyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1401);

    auto g_yyyyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1402);

    auto g_yyyyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1403);

    auto g_yyyyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1404);

    auto g_yyyyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1405);

    auto g_yyyyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1406);

    auto g_yyyyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1407);

    auto g_yyyyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1408);

    auto g_yyyyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1409);

    auto g_yyyyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1410);

    auto g_yyyyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1411);

    auto g_yyyyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1412);

    auto g_yyyyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1413);

    auto g_yyyyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1414);

    auto g_yyyyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1415);

    auto g_yyyyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1416);

    auto g_yyyyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1417);

    auto g_yyyyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1418);

    auto g_yyyyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1419);

    auto g_yyyyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1420);

    auto g_yyyyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1421);

    auto g_yyyyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1422);

    auto g_yyyyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1423);

    auto g_yyyyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1424);

    auto g_yyyyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1425);

    auto g_yyyyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1426);

    auto g_yyyyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1427);

    auto g_yyyyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1428);

    auto g_yyyyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1429);

    auto g_yyyyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1430);

    auto g_yyyyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1431);

    auto g_yyyyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1432);

    auto g_yyyyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1433);

    auto g_yyyyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1434);

    auto g_yyyyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1435);

    auto g_yyyyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1436);

    auto g_yyyyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1437);

    auto g_yyyyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1438);

    auto g_yyyyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1439);

    #pragma omp simd aligned(g_yyyyz_0_xxxxxxxy_0, g_yyyyz_0_xxxxxxxy_1, g_yyyyz_0_xxxxxxyy_0, g_yyyyz_0_xxxxxxyy_1, g_yyyyz_0_xxxxxyyy_0, g_yyyyz_0_xxxxxyyy_1, g_yyyyz_0_xxxxyyyy_0, g_yyyyz_0_xxxxyyyy_1, g_yyyyz_0_xxxyyyyy_0, g_yyyyz_0_xxxyyyyy_1, g_yyyyz_0_xxyyyyyy_0, g_yyyyz_0_xxyyyyyy_1, g_yyyyz_0_xyyyyyyy_0, g_yyyyz_0_xyyyyyyy_1, g_yyyyz_0_yyyyyyyy_0, g_yyyyz_0_yyyyyyyy_1, g_yyyyzz_0_xxxxxxxy_1, g_yyyyzz_0_xxxxxxyy_1, g_yyyyzz_0_xxxxxyyy_1, g_yyyyzz_0_xxxxyyyy_1, g_yyyyzz_0_xxxyyyyy_1, g_yyyyzz_0_xxyyyyyy_1, g_yyyyzz_0_xyyyyyyy_1, g_yyyyzz_0_yyyyyyyy_1, g_yyyyzzz_0_xxxxxxxx_0, g_yyyyzzz_0_xxxxxxxy_0, g_yyyyzzz_0_xxxxxxxz_0, g_yyyyzzz_0_xxxxxxyy_0, g_yyyyzzz_0_xxxxxxyz_0, g_yyyyzzz_0_xxxxxxzz_0, g_yyyyzzz_0_xxxxxyyy_0, g_yyyyzzz_0_xxxxxyyz_0, g_yyyyzzz_0_xxxxxyzz_0, g_yyyyzzz_0_xxxxxzzz_0, g_yyyyzzz_0_xxxxyyyy_0, g_yyyyzzz_0_xxxxyyyz_0, g_yyyyzzz_0_xxxxyyzz_0, g_yyyyzzz_0_xxxxyzzz_0, g_yyyyzzz_0_xxxxzzzz_0, g_yyyyzzz_0_xxxyyyyy_0, g_yyyyzzz_0_xxxyyyyz_0, g_yyyyzzz_0_xxxyyyzz_0, g_yyyyzzz_0_xxxyyzzz_0, g_yyyyzzz_0_xxxyzzzz_0, g_yyyyzzz_0_xxxzzzzz_0, g_yyyyzzz_0_xxyyyyyy_0, g_yyyyzzz_0_xxyyyyyz_0, g_yyyyzzz_0_xxyyyyzz_0, g_yyyyzzz_0_xxyyyzzz_0, g_yyyyzzz_0_xxyyzzzz_0, g_yyyyzzz_0_xxyzzzzz_0, g_yyyyzzz_0_xxzzzzzz_0, g_yyyyzzz_0_xyyyyyyy_0, g_yyyyzzz_0_xyyyyyyz_0, g_yyyyzzz_0_xyyyyyzz_0, g_yyyyzzz_0_xyyyyzzz_0, g_yyyyzzz_0_xyyyzzzz_0, g_yyyyzzz_0_xyyzzzzz_0, g_yyyyzzz_0_xyzzzzzz_0, g_yyyyzzz_0_xzzzzzzz_0, g_yyyyzzz_0_yyyyyyyy_0, g_yyyyzzz_0_yyyyyyyz_0, g_yyyyzzz_0_yyyyyyzz_0, g_yyyyzzz_0_yyyyyzzz_0, g_yyyyzzz_0_yyyyzzzz_0, g_yyyyzzz_0_yyyzzzzz_0, g_yyyyzzz_0_yyzzzzzz_0, g_yyyyzzz_0_yzzzzzzz_0, g_yyyyzzz_0_zzzzzzzz_0, g_yyyzzz_0_xxxxxxxx_1, g_yyyzzz_0_xxxxxxxz_1, g_yyyzzz_0_xxxxxxyz_1, g_yyyzzz_0_xxxxxxz_1, g_yyyzzz_0_xxxxxxzz_1, g_yyyzzz_0_xxxxxyyz_1, g_yyyzzz_0_xxxxxyz_1, g_yyyzzz_0_xxxxxyzz_1, g_yyyzzz_0_xxxxxzz_1, g_yyyzzz_0_xxxxxzzz_1, g_yyyzzz_0_xxxxyyyz_1, g_yyyzzz_0_xxxxyyz_1, g_yyyzzz_0_xxxxyyzz_1, g_yyyzzz_0_xxxxyzz_1, g_yyyzzz_0_xxxxyzzz_1, g_yyyzzz_0_xxxxzzz_1, g_yyyzzz_0_xxxxzzzz_1, g_yyyzzz_0_xxxyyyyz_1, g_yyyzzz_0_xxxyyyz_1, g_yyyzzz_0_xxxyyyzz_1, g_yyyzzz_0_xxxyyzz_1, g_yyyzzz_0_xxxyyzzz_1, g_yyyzzz_0_xxxyzzz_1, g_yyyzzz_0_xxxyzzzz_1, g_yyyzzz_0_xxxzzzz_1, g_yyyzzz_0_xxxzzzzz_1, g_yyyzzz_0_xxyyyyyz_1, g_yyyzzz_0_xxyyyyz_1, g_yyyzzz_0_xxyyyyzz_1, g_yyyzzz_0_xxyyyzz_1, g_yyyzzz_0_xxyyyzzz_1, g_yyyzzz_0_xxyyzzz_1, g_yyyzzz_0_xxyyzzzz_1, g_yyyzzz_0_xxyzzzz_1, g_yyyzzz_0_xxyzzzzz_1, g_yyyzzz_0_xxzzzzz_1, g_yyyzzz_0_xxzzzzzz_1, g_yyyzzz_0_xyyyyyyz_1, g_yyyzzz_0_xyyyyyz_1, g_yyyzzz_0_xyyyyyzz_1, g_yyyzzz_0_xyyyyzz_1, g_yyyzzz_0_xyyyyzzz_1, g_yyyzzz_0_xyyyzzz_1, g_yyyzzz_0_xyyyzzzz_1, g_yyyzzz_0_xyyzzzz_1, g_yyyzzz_0_xyyzzzzz_1, g_yyyzzz_0_xyzzzzz_1, g_yyyzzz_0_xyzzzzzz_1, g_yyyzzz_0_xzzzzzz_1, g_yyyzzz_0_xzzzzzzz_1, g_yyyzzz_0_yyyyyyyz_1, g_yyyzzz_0_yyyyyyz_1, g_yyyzzz_0_yyyyyyzz_1, g_yyyzzz_0_yyyyyzz_1, g_yyyzzz_0_yyyyyzzz_1, g_yyyzzz_0_yyyyzzz_1, g_yyyzzz_0_yyyyzzzz_1, g_yyyzzz_0_yyyzzzz_1, g_yyyzzz_0_yyyzzzzz_1, g_yyyzzz_0_yyzzzzz_1, g_yyyzzz_0_yyzzzzzz_1, g_yyyzzz_0_yzzzzzz_1, g_yyyzzz_0_yzzzzzzz_1, g_yyyzzz_0_zzzzzzz_1, g_yyyzzz_0_zzzzzzzz_1, g_yyzzz_0_xxxxxxxx_0, g_yyzzz_0_xxxxxxxx_1, g_yyzzz_0_xxxxxxxz_0, g_yyzzz_0_xxxxxxxz_1, g_yyzzz_0_xxxxxxyz_0, g_yyzzz_0_xxxxxxyz_1, g_yyzzz_0_xxxxxxzz_0, g_yyzzz_0_xxxxxxzz_1, g_yyzzz_0_xxxxxyyz_0, g_yyzzz_0_xxxxxyyz_1, g_yyzzz_0_xxxxxyzz_0, g_yyzzz_0_xxxxxyzz_1, g_yyzzz_0_xxxxxzzz_0, g_yyzzz_0_xxxxxzzz_1, g_yyzzz_0_xxxxyyyz_0, g_yyzzz_0_xxxxyyyz_1, g_yyzzz_0_xxxxyyzz_0, g_yyzzz_0_xxxxyyzz_1, g_yyzzz_0_xxxxyzzz_0, g_yyzzz_0_xxxxyzzz_1, g_yyzzz_0_xxxxzzzz_0, g_yyzzz_0_xxxxzzzz_1, g_yyzzz_0_xxxyyyyz_0, g_yyzzz_0_xxxyyyyz_1, g_yyzzz_0_xxxyyyzz_0, g_yyzzz_0_xxxyyyzz_1, g_yyzzz_0_xxxyyzzz_0, g_yyzzz_0_xxxyyzzz_1, g_yyzzz_0_xxxyzzzz_0, g_yyzzz_0_xxxyzzzz_1, g_yyzzz_0_xxxzzzzz_0, g_yyzzz_0_xxxzzzzz_1, g_yyzzz_0_xxyyyyyz_0, g_yyzzz_0_xxyyyyyz_1, g_yyzzz_0_xxyyyyzz_0, g_yyzzz_0_xxyyyyzz_1, g_yyzzz_0_xxyyyzzz_0, g_yyzzz_0_xxyyyzzz_1, g_yyzzz_0_xxyyzzzz_0, g_yyzzz_0_xxyyzzzz_1, g_yyzzz_0_xxyzzzzz_0, g_yyzzz_0_xxyzzzzz_1, g_yyzzz_0_xxzzzzzz_0, g_yyzzz_0_xxzzzzzz_1, g_yyzzz_0_xyyyyyyz_0, g_yyzzz_0_xyyyyyyz_1, g_yyzzz_0_xyyyyyzz_0, g_yyzzz_0_xyyyyyzz_1, g_yyzzz_0_xyyyyzzz_0, g_yyzzz_0_xyyyyzzz_1, g_yyzzz_0_xyyyzzzz_0, g_yyzzz_0_xyyyzzzz_1, g_yyzzz_0_xyyzzzzz_0, g_yyzzz_0_xyyzzzzz_1, g_yyzzz_0_xyzzzzzz_0, g_yyzzz_0_xyzzzzzz_1, g_yyzzz_0_xzzzzzzz_0, g_yyzzz_0_xzzzzzzz_1, g_yyzzz_0_yyyyyyyz_0, g_yyzzz_0_yyyyyyyz_1, g_yyzzz_0_yyyyyyzz_0, g_yyzzz_0_yyyyyyzz_1, g_yyzzz_0_yyyyyzzz_0, g_yyzzz_0_yyyyyzzz_1, g_yyzzz_0_yyyyzzzz_0, g_yyzzz_0_yyyyzzzz_1, g_yyzzz_0_yyyzzzzz_0, g_yyzzz_0_yyyzzzzz_1, g_yyzzz_0_yyzzzzzz_0, g_yyzzz_0_yyzzzzzz_1, g_yyzzz_0_yzzzzzzz_0, g_yyzzz_0_yzzzzzzz_1, g_yyzzz_0_zzzzzzzz_0, g_yyzzz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzzz_0_xxxxxxxx_0[i] = 3.0 * g_yyzzz_0_xxxxxxxx_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxxxx_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxxxxy_0[i] = 2.0 * g_yyyyz_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxxxxxy_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxxxxxz_0[i] = 3.0 * g_yyzzz_0_xxxxxxxz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxxxz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxxxyy_0[i] = 2.0 * g_yyyyz_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxxxxyy_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxxxxyz_0[i] = 3.0 * g_yyzzz_0_xxxxxxyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxxyz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxxxzz_0[i] = 3.0 * g_yyzzz_0_xxxxxxzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxxzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxxyyy_0[i] = 2.0 * g_yyyyz_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxxxyyy_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxxxyyz_0[i] = 3.0 * g_yyzzz_0_xxxxxyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxxyzz_0[i] = 3.0 * g_yyzzz_0_xxxxxyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxyzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxxzzz_0[i] = 3.0 * g_yyzzz_0_xxxxxzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxyyyy_0[i] = 2.0 * g_yyyyz_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxxyyyy_1[i] * fz_be_0 + g_yyyyzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxxyyyz_0[i] = 3.0 * g_yyzzz_0_xxxxyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxyyzz_0[i] = 3.0 * g_yyzzz_0_xxxxyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxyzzz_0[i] = 3.0 * g_yyzzz_0_xxxxyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxyzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxzzzz_0[i] = 3.0 * g_yyzzz_0_xxxxzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxyyyyy_0[i] = 2.0 * g_yyyyz_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxyyyyy_1[i] * fz_be_0 + g_yyyyzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxyyyyz_0[i] = 3.0 * g_yyzzz_0_xxxyyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxyyyzz_0[i] = 3.0 * g_yyzzz_0_xxxyyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxyyzzz_0[i] = 3.0 * g_yyzzz_0_xxxyyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxyzzzz_0[i] = 3.0 * g_yyzzz_0_xxxyzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxyzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxzzzzz_0[i] = 3.0 * g_yyzzz_0_xxxzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyyyyyy_0[i] = 2.0 * g_yyyyz_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxyyyyyy_1[i] * fz_be_0 + g_yyyyzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxyyyyyz_0[i] = 3.0 * g_yyzzz_0_xxyyyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyyyyzz_0[i] = 3.0 * g_yyzzz_0_xxyyyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyyyzzz_0[i] = 3.0 * g_yyzzz_0_xxyyyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyyzzzz_0[i] = 3.0 * g_yyzzz_0_xxyyzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyzzzzz_0[i] = 3.0 * g_yyzzz_0_xxyzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxzzzzzz_0[i] = 3.0 * g_yyzzz_0_xxzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxzzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyyyyyy_0[i] = 2.0 * g_yyyyz_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xyyyyyyy_1[i] * fz_be_0 + g_yyyyzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xyyyyyyz_0[i] = 3.0 * g_yyzzz_0_xyyyyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyyzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyyyyzz_0[i] = 3.0 * g_yyzzz_0_xyyyyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyyzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyyyzzz_0[i] = 3.0 * g_yyzzz_0_xyyyyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyyzzzz_0[i] = 3.0 * g_yyzzz_0_xyyyzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyzzzzz_0[i] = 3.0 * g_yyzzz_0_xyyzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyzzzzzz_0[i] = 3.0 * g_yyzzz_0_xyzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyzzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xzzzzzzz_0[i] = 3.0 * g_yyzzz_0_xzzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xzzzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyyyyyy_0[i] = 2.0 * g_yyyyz_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_yyyyyyyy_1[i] * fz_be_0 + g_yyyyzz_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_yyyyyyyz_0[i] = 3.0 * g_yyzzz_0_yyyyyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yyyzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyyyyzz_0[i] = 3.0 * g_yyzzz_0_yyyyyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yyyzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyyyzzz_0[i] = 3.0 * g_yyzzz_0_yyyyyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yyyzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyyzzzz_0[i] = 3.0 * g_yyzzz_0_yyyyzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyzzzzz_0[i] = 3.0 * g_yyzzz_0_yyyzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyzzzzzz_0[i] = 3.0 * g_yyzzz_0_yyzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yzzzzzzz_0[i] = 3.0 * g_yyzzz_0_yzzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yzzzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_zzzzzzzz_0[i] = 3.0 * g_yyzzz_0_zzzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_zzzzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1440-1485 components of targeted buffer : KSL

    auto g_yyyzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1440);

    auto g_yyyzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1441);

    auto g_yyyzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1442);

    auto g_yyyzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1443);

    auto g_yyyzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1444);

    auto g_yyyzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1445);

    auto g_yyyzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1446);

    auto g_yyyzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1447);

    auto g_yyyzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1448);

    auto g_yyyzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1449);

    auto g_yyyzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1450);

    auto g_yyyzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1451);

    auto g_yyyzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1452);

    auto g_yyyzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1453);

    auto g_yyyzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1454);

    auto g_yyyzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1455);

    auto g_yyyzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1456);

    auto g_yyyzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1457);

    auto g_yyyzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1458);

    auto g_yyyzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1459);

    auto g_yyyzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1460);

    auto g_yyyzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1461);

    auto g_yyyzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1462);

    auto g_yyyzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1463);

    auto g_yyyzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1464);

    auto g_yyyzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1465);

    auto g_yyyzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1466);

    auto g_yyyzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1467);

    auto g_yyyzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1468);

    auto g_yyyzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1469);

    auto g_yyyzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1470);

    auto g_yyyzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1471);

    auto g_yyyzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1472);

    auto g_yyyzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1473);

    auto g_yyyzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1474);

    auto g_yyyzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1475);

    auto g_yyyzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1476);

    auto g_yyyzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1477);

    auto g_yyyzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1478);

    auto g_yyyzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1479);

    auto g_yyyzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1480);

    auto g_yyyzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1481);

    auto g_yyyzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1482);

    auto g_yyyzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1483);

    auto g_yyyzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1484);

    #pragma omp simd aligned(g_yyyzz_0_xxxxxxxy_0, g_yyyzz_0_xxxxxxxy_1, g_yyyzz_0_xxxxxxyy_0, g_yyyzz_0_xxxxxxyy_1, g_yyyzz_0_xxxxxyyy_0, g_yyyzz_0_xxxxxyyy_1, g_yyyzz_0_xxxxyyyy_0, g_yyyzz_0_xxxxyyyy_1, g_yyyzz_0_xxxyyyyy_0, g_yyyzz_0_xxxyyyyy_1, g_yyyzz_0_xxyyyyyy_0, g_yyyzz_0_xxyyyyyy_1, g_yyyzz_0_xyyyyyyy_0, g_yyyzz_0_xyyyyyyy_1, g_yyyzz_0_yyyyyyyy_0, g_yyyzz_0_yyyyyyyy_1, g_yyyzzz_0_xxxxxxxy_1, g_yyyzzz_0_xxxxxxyy_1, g_yyyzzz_0_xxxxxyyy_1, g_yyyzzz_0_xxxxyyyy_1, g_yyyzzz_0_xxxyyyyy_1, g_yyyzzz_0_xxyyyyyy_1, g_yyyzzz_0_xyyyyyyy_1, g_yyyzzz_0_yyyyyyyy_1, g_yyyzzzz_0_xxxxxxxx_0, g_yyyzzzz_0_xxxxxxxy_0, g_yyyzzzz_0_xxxxxxxz_0, g_yyyzzzz_0_xxxxxxyy_0, g_yyyzzzz_0_xxxxxxyz_0, g_yyyzzzz_0_xxxxxxzz_0, g_yyyzzzz_0_xxxxxyyy_0, g_yyyzzzz_0_xxxxxyyz_0, g_yyyzzzz_0_xxxxxyzz_0, g_yyyzzzz_0_xxxxxzzz_0, g_yyyzzzz_0_xxxxyyyy_0, g_yyyzzzz_0_xxxxyyyz_0, g_yyyzzzz_0_xxxxyyzz_0, g_yyyzzzz_0_xxxxyzzz_0, g_yyyzzzz_0_xxxxzzzz_0, g_yyyzzzz_0_xxxyyyyy_0, g_yyyzzzz_0_xxxyyyyz_0, g_yyyzzzz_0_xxxyyyzz_0, g_yyyzzzz_0_xxxyyzzz_0, g_yyyzzzz_0_xxxyzzzz_0, g_yyyzzzz_0_xxxzzzzz_0, g_yyyzzzz_0_xxyyyyyy_0, g_yyyzzzz_0_xxyyyyyz_0, g_yyyzzzz_0_xxyyyyzz_0, g_yyyzzzz_0_xxyyyzzz_0, g_yyyzzzz_0_xxyyzzzz_0, g_yyyzzzz_0_xxyzzzzz_0, g_yyyzzzz_0_xxzzzzzz_0, g_yyyzzzz_0_xyyyyyyy_0, g_yyyzzzz_0_xyyyyyyz_0, g_yyyzzzz_0_xyyyyyzz_0, g_yyyzzzz_0_xyyyyzzz_0, g_yyyzzzz_0_xyyyzzzz_0, g_yyyzzzz_0_xyyzzzzz_0, g_yyyzzzz_0_xyzzzzzz_0, g_yyyzzzz_0_xzzzzzzz_0, g_yyyzzzz_0_yyyyyyyy_0, g_yyyzzzz_0_yyyyyyyz_0, g_yyyzzzz_0_yyyyyyzz_0, g_yyyzzzz_0_yyyyyzzz_0, g_yyyzzzz_0_yyyyzzzz_0, g_yyyzzzz_0_yyyzzzzz_0, g_yyyzzzz_0_yyzzzzzz_0, g_yyyzzzz_0_yzzzzzzz_0, g_yyyzzzz_0_zzzzzzzz_0, g_yyzzzz_0_xxxxxxxx_1, g_yyzzzz_0_xxxxxxxz_1, g_yyzzzz_0_xxxxxxyz_1, g_yyzzzz_0_xxxxxxz_1, g_yyzzzz_0_xxxxxxzz_1, g_yyzzzz_0_xxxxxyyz_1, g_yyzzzz_0_xxxxxyz_1, g_yyzzzz_0_xxxxxyzz_1, g_yyzzzz_0_xxxxxzz_1, g_yyzzzz_0_xxxxxzzz_1, g_yyzzzz_0_xxxxyyyz_1, g_yyzzzz_0_xxxxyyz_1, g_yyzzzz_0_xxxxyyzz_1, g_yyzzzz_0_xxxxyzz_1, g_yyzzzz_0_xxxxyzzz_1, g_yyzzzz_0_xxxxzzz_1, g_yyzzzz_0_xxxxzzzz_1, g_yyzzzz_0_xxxyyyyz_1, g_yyzzzz_0_xxxyyyz_1, g_yyzzzz_0_xxxyyyzz_1, g_yyzzzz_0_xxxyyzz_1, g_yyzzzz_0_xxxyyzzz_1, g_yyzzzz_0_xxxyzzz_1, g_yyzzzz_0_xxxyzzzz_1, g_yyzzzz_0_xxxzzzz_1, g_yyzzzz_0_xxxzzzzz_1, g_yyzzzz_0_xxyyyyyz_1, g_yyzzzz_0_xxyyyyz_1, g_yyzzzz_0_xxyyyyzz_1, g_yyzzzz_0_xxyyyzz_1, g_yyzzzz_0_xxyyyzzz_1, g_yyzzzz_0_xxyyzzz_1, g_yyzzzz_0_xxyyzzzz_1, g_yyzzzz_0_xxyzzzz_1, g_yyzzzz_0_xxyzzzzz_1, g_yyzzzz_0_xxzzzzz_1, g_yyzzzz_0_xxzzzzzz_1, g_yyzzzz_0_xyyyyyyz_1, g_yyzzzz_0_xyyyyyz_1, g_yyzzzz_0_xyyyyyzz_1, g_yyzzzz_0_xyyyyzz_1, g_yyzzzz_0_xyyyyzzz_1, g_yyzzzz_0_xyyyzzz_1, g_yyzzzz_0_xyyyzzzz_1, g_yyzzzz_0_xyyzzzz_1, g_yyzzzz_0_xyyzzzzz_1, g_yyzzzz_0_xyzzzzz_1, g_yyzzzz_0_xyzzzzzz_1, g_yyzzzz_0_xzzzzzz_1, g_yyzzzz_0_xzzzzzzz_1, g_yyzzzz_0_yyyyyyyz_1, g_yyzzzz_0_yyyyyyz_1, g_yyzzzz_0_yyyyyyzz_1, g_yyzzzz_0_yyyyyzz_1, g_yyzzzz_0_yyyyyzzz_1, g_yyzzzz_0_yyyyzzz_1, g_yyzzzz_0_yyyyzzzz_1, g_yyzzzz_0_yyyzzzz_1, g_yyzzzz_0_yyyzzzzz_1, g_yyzzzz_0_yyzzzzz_1, g_yyzzzz_0_yyzzzzzz_1, g_yyzzzz_0_yzzzzzz_1, g_yyzzzz_0_yzzzzzzz_1, g_yyzzzz_0_zzzzzzz_1, g_yyzzzz_0_zzzzzzzz_1, g_yzzzz_0_xxxxxxxx_0, g_yzzzz_0_xxxxxxxx_1, g_yzzzz_0_xxxxxxxz_0, g_yzzzz_0_xxxxxxxz_1, g_yzzzz_0_xxxxxxyz_0, g_yzzzz_0_xxxxxxyz_1, g_yzzzz_0_xxxxxxzz_0, g_yzzzz_0_xxxxxxzz_1, g_yzzzz_0_xxxxxyyz_0, g_yzzzz_0_xxxxxyyz_1, g_yzzzz_0_xxxxxyzz_0, g_yzzzz_0_xxxxxyzz_1, g_yzzzz_0_xxxxxzzz_0, g_yzzzz_0_xxxxxzzz_1, g_yzzzz_0_xxxxyyyz_0, g_yzzzz_0_xxxxyyyz_1, g_yzzzz_0_xxxxyyzz_0, g_yzzzz_0_xxxxyyzz_1, g_yzzzz_0_xxxxyzzz_0, g_yzzzz_0_xxxxyzzz_1, g_yzzzz_0_xxxxzzzz_0, g_yzzzz_0_xxxxzzzz_1, g_yzzzz_0_xxxyyyyz_0, g_yzzzz_0_xxxyyyyz_1, g_yzzzz_0_xxxyyyzz_0, g_yzzzz_0_xxxyyyzz_1, g_yzzzz_0_xxxyyzzz_0, g_yzzzz_0_xxxyyzzz_1, g_yzzzz_0_xxxyzzzz_0, g_yzzzz_0_xxxyzzzz_1, g_yzzzz_0_xxxzzzzz_0, g_yzzzz_0_xxxzzzzz_1, g_yzzzz_0_xxyyyyyz_0, g_yzzzz_0_xxyyyyyz_1, g_yzzzz_0_xxyyyyzz_0, g_yzzzz_0_xxyyyyzz_1, g_yzzzz_0_xxyyyzzz_0, g_yzzzz_0_xxyyyzzz_1, g_yzzzz_0_xxyyzzzz_0, g_yzzzz_0_xxyyzzzz_1, g_yzzzz_0_xxyzzzzz_0, g_yzzzz_0_xxyzzzzz_1, g_yzzzz_0_xxzzzzzz_0, g_yzzzz_0_xxzzzzzz_1, g_yzzzz_0_xyyyyyyz_0, g_yzzzz_0_xyyyyyyz_1, g_yzzzz_0_xyyyyyzz_0, g_yzzzz_0_xyyyyyzz_1, g_yzzzz_0_xyyyyzzz_0, g_yzzzz_0_xyyyyzzz_1, g_yzzzz_0_xyyyzzzz_0, g_yzzzz_0_xyyyzzzz_1, g_yzzzz_0_xyyzzzzz_0, g_yzzzz_0_xyyzzzzz_1, g_yzzzz_0_xyzzzzzz_0, g_yzzzz_0_xyzzzzzz_1, g_yzzzz_0_xzzzzzzz_0, g_yzzzz_0_xzzzzzzz_1, g_yzzzz_0_yyyyyyyz_0, g_yzzzz_0_yyyyyyyz_1, g_yzzzz_0_yyyyyyzz_0, g_yzzzz_0_yyyyyyzz_1, g_yzzzz_0_yyyyyzzz_0, g_yzzzz_0_yyyyyzzz_1, g_yzzzz_0_yyyyzzzz_0, g_yzzzz_0_yyyyzzzz_1, g_yzzzz_0_yyyzzzzz_0, g_yzzzz_0_yyyzzzzz_1, g_yzzzz_0_yyzzzzzz_0, g_yzzzz_0_yyzzzzzz_1, g_yzzzz_0_yzzzzzzz_0, g_yzzzz_0_yzzzzzzz_1, g_yzzzz_0_zzzzzzzz_0, g_yzzzz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzzz_0_xxxxxxxx_0[i] = 2.0 * g_yzzzz_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxxxx_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxxxxy_0[i] = 3.0 * g_yyyzz_0_xxxxxxxy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxxxxxy_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxxxxxz_0[i] = 2.0 * g_yzzzz_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxxxz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxxxyy_0[i] = 3.0 * g_yyyzz_0_xxxxxxyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxxxxyy_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxxxxyz_0[i] = 2.0 * g_yzzzz_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxxyz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxxxzz_0[i] = 2.0 * g_yzzzz_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxxzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxxyyy_0[i] = 3.0 * g_yyyzz_0_xxxxxyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxxxyyy_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxxxyyz_0[i] = 2.0 * g_yzzzz_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxxyzz_0[i] = 2.0 * g_yzzzz_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxyzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxxzzz_0[i] = 2.0 * g_yzzzz_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxyyyy_0[i] = 3.0 * g_yyyzz_0_xxxxyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxxyyyy_1[i] * fz_be_0 + g_yyyzzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxxyyyz_0[i] = 2.0 * g_yzzzz_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxyyzz_0[i] = 2.0 * g_yzzzz_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxyzzz_0[i] = 2.0 * g_yzzzz_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxyzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxzzzz_0[i] = 2.0 * g_yzzzz_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxyyyyy_0[i] = 3.0 * g_yyyzz_0_xxxyyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxyyyyy_1[i] * fz_be_0 + g_yyyzzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxyyyyz_0[i] = 2.0 * g_yzzzz_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxyyyzz_0[i] = 2.0 * g_yzzzz_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxyyzzz_0[i] = 2.0 * g_yzzzz_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxyzzzz_0[i] = 2.0 * g_yzzzz_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxyzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxzzzzz_0[i] = 2.0 * g_yzzzz_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyyyyyy_0[i] = 3.0 * g_yyyzz_0_xxyyyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxyyyyyy_1[i] * fz_be_0 + g_yyyzzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxyyyyyz_0[i] = 2.0 * g_yzzzz_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyyyyzz_0[i] = 2.0 * g_yzzzz_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyyyzzz_0[i] = 2.0 * g_yzzzz_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyyzzzz_0[i] = 2.0 * g_yzzzz_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyzzzzz_0[i] = 2.0 * g_yzzzz_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxzzzzzz_0[i] = 2.0 * g_yzzzz_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxzzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyyyyyy_0[i] = 3.0 * g_yyyzz_0_xyyyyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xyyyyyyy_1[i] * fz_be_0 + g_yyyzzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xyyyyyyz_0[i] = 2.0 * g_yzzzz_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyyyyzz_0[i] = 2.0 * g_yzzzz_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyyyzzz_0[i] = 2.0 * g_yzzzz_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyyzzzz_0[i] = 2.0 * g_yzzzz_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyzzzzz_0[i] = 2.0 * g_yzzzz_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyzzzzzz_0[i] = 2.0 * g_yzzzz_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyzzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xzzzzzzz_0[i] = 2.0 * g_yzzzz_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xzzzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyyyyyy_0[i] = 3.0 * g_yyyzz_0_yyyyyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_yyyyyyyy_1[i] * fz_be_0 + g_yyyzzz_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_yyyyyyyz_0[i] = 2.0 * g_yzzzz_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yyzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyyyyzz_0[i] = 2.0 * g_yzzzz_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yyzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyyyzzz_0[i] = 2.0 * g_yzzzz_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yyzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyyzzzz_0[i] = 2.0 * g_yzzzz_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyzzzzz_0[i] = 2.0 * g_yzzzz_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyzzzzzz_0[i] = 2.0 * g_yzzzz_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yzzzzzzz_0[i] = 2.0 * g_yzzzz_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yzzzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_zzzzzzzz_0[i] = 2.0 * g_yzzzz_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_zzzzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1485-1530 components of targeted buffer : KSL

    auto g_yyzzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1485);

    auto g_yyzzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1486);

    auto g_yyzzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1487);

    auto g_yyzzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1488);

    auto g_yyzzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1489);

    auto g_yyzzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1490);

    auto g_yyzzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1491);

    auto g_yyzzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1492);

    auto g_yyzzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1493);

    auto g_yyzzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1494);

    auto g_yyzzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1495);

    auto g_yyzzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1496);

    auto g_yyzzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1497);

    auto g_yyzzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1498);

    auto g_yyzzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1499);

    auto g_yyzzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1500);

    auto g_yyzzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1501);

    auto g_yyzzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1502);

    auto g_yyzzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1503);

    auto g_yyzzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1504);

    auto g_yyzzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1505);

    auto g_yyzzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1506);

    auto g_yyzzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1507);

    auto g_yyzzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1508);

    auto g_yyzzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1509);

    auto g_yyzzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1510);

    auto g_yyzzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1511);

    auto g_yyzzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1512);

    auto g_yyzzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1513);

    auto g_yyzzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1514);

    auto g_yyzzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1515);

    auto g_yyzzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1516);

    auto g_yyzzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1517);

    auto g_yyzzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1518);

    auto g_yyzzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1519);

    auto g_yyzzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1520);

    auto g_yyzzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1521);

    auto g_yyzzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1522);

    auto g_yyzzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1523);

    auto g_yyzzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1524);

    auto g_yyzzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1525);

    auto g_yyzzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1526);

    auto g_yyzzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1527);

    auto g_yyzzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1528);

    auto g_yyzzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1529);

    #pragma omp simd aligned(g_yyzzz_0_xxxxxxxy_0, g_yyzzz_0_xxxxxxxy_1, g_yyzzz_0_xxxxxxyy_0, g_yyzzz_0_xxxxxxyy_1, g_yyzzz_0_xxxxxyyy_0, g_yyzzz_0_xxxxxyyy_1, g_yyzzz_0_xxxxyyyy_0, g_yyzzz_0_xxxxyyyy_1, g_yyzzz_0_xxxyyyyy_0, g_yyzzz_0_xxxyyyyy_1, g_yyzzz_0_xxyyyyyy_0, g_yyzzz_0_xxyyyyyy_1, g_yyzzz_0_xyyyyyyy_0, g_yyzzz_0_xyyyyyyy_1, g_yyzzz_0_yyyyyyyy_0, g_yyzzz_0_yyyyyyyy_1, g_yyzzzz_0_xxxxxxxy_1, g_yyzzzz_0_xxxxxxyy_1, g_yyzzzz_0_xxxxxyyy_1, g_yyzzzz_0_xxxxyyyy_1, g_yyzzzz_0_xxxyyyyy_1, g_yyzzzz_0_xxyyyyyy_1, g_yyzzzz_0_xyyyyyyy_1, g_yyzzzz_0_yyyyyyyy_1, g_yyzzzzz_0_xxxxxxxx_0, g_yyzzzzz_0_xxxxxxxy_0, g_yyzzzzz_0_xxxxxxxz_0, g_yyzzzzz_0_xxxxxxyy_0, g_yyzzzzz_0_xxxxxxyz_0, g_yyzzzzz_0_xxxxxxzz_0, g_yyzzzzz_0_xxxxxyyy_0, g_yyzzzzz_0_xxxxxyyz_0, g_yyzzzzz_0_xxxxxyzz_0, g_yyzzzzz_0_xxxxxzzz_0, g_yyzzzzz_0_xxxxyyyy_0, g_yyzzzzz_0_xxxxyyyz_0, g_yyzzzzz_0_xxxxyyzz_0, g_yyzzzzz_0_xxxxyzzz_0, g_yyzzzzz_0_xxxxzzzz_0, g_yyzzzzz_0_xxxyyyyy_0, g_yyzzzzz_0_xxxyyyyz_0, g_yyzzzzz_0_xxxyyyzz_0, g_yyzzzzz_0_xxxyyzzz_0, g_yyzzzzz_0_xxxyzzzz_0, g_yyzzzzz_0_xxxzzzzz_0, g_yyzzzzz_0_xxyyyyyy_0, g_yyzzzzz_0_xxyyyyyz_0, g_yyzzzzz_0_xxyyyyzz_0, g_yyzzzzz_0_xxyyyzzz_0, g_yyzzzzz_0_xxyyzzzz_0, g_yyzzzzz_0_xxyzzzzz_0, g_yyzzzzz_0_xxzzzzzz_0, g_yyzzzzz_0_xyyyyyyy_0, g_yyzzzzz_0_xyyyyyyz_0, g_yyzzzzz_0_xyyyyyzz_0, g_yyzzzzz_0_xyyyyzzz_0, g_yyzzzzz_0_xyyyzzzz_0, g_yyzzzzz_0_xyyzzzzz_0, g_yyzzzzz_0_xyzzzzzz_0, g_yyzzzzz_0_xzzzzzzz_0, g_yyzzzzz_0_yyyyyyyy_0, g_yyzzzzz_0_yyyyyyyz_0, g_yyzzzzz_0_yyyyyyzz_0, g_yyzzzzz_0_yyyyyzzz_0, g_yyzzzzz_0_yyyyzzzz_0, g_yyzzzzz_0_yyyzzzzz_0, g_yyzzzzz_0_yyzzzzzz_0, g_yyzzzzz_0_yzzzzzzz_0, g_yyzzzzz_0_zzzzzzzz_0, g_yzzzzz_0_xxxxxxxx_1, g_yzzzzz_0_xxxxxxxz_1, g_yzzzzz_0_xxxxxxyz_1, g_yzzzzz_0_xxxxxxz_1, g_yzzzzz_0_xxxxxxzz_1, g_yzzzzz_0_xxxxxyyz_1, g_yzzzzz_0_xxxxxyz_1, g_yzzzzz_0_xxxxxyzz_1, g_yzzzzz_0_xxxxxzz_1, g_yzzzzz_0_xxxxxzzz_1, g_yzzzzz_0_xxxxyyyz_1, g_yzzzzz_0_xxxxyyz_1, g_yzzzzz_0_xxxxyyzz_1, g_yzzzzz_0_xxxxyzz_1, g_yzzzzz_0_xxxxyzzz_1, g_yzzzzz_0_xxxxzzz_1, g_yzzzzz_0_xxxxzzzz_1, g_yzzzzz_0_xxxyyyyz_1, g_yzzzzz_0_xxxyyyz_1, g_yzzzzz_0_xxxyyyzz_1, g_yzzzzz_0_xxxyyzz_1, g_yzzzzz_0_xxxyyzzz_1, g_yzzzzz_0_xxxyzzz_1, g_yzzzzz_0_xxxyzzzz_1, g_yzzzzz_0_xxxzzzz_1, g_yzzzzz_0_xxxzzzzz_1, g_yzzzzz_0_xxyyyyyz_1, g_yzzzzz_0_xxyyyyz_1, g_yzzzzz_0_xxyyyyzz_1, g_yzzzzz_0_xxyyyzz_1, g_yzzzzz_0_xxyyyzzz_1, g_yzzzzz_0_xxyyzzz_1, g_yzzzzz_0_xxyyzzzz_1, g_yzzzzz_0_xxyzzzz_1, g_yzzzzz_0_xxyzzzzz_1, g_yzzzzz_0_xxzzzzz_1, g_yzzzzz_0_xxzzzzzz_1, g_yzzzzz_0_xyyyyyyz_1, g_yzzzzz_0_xyyyyyz_1, g_yzzzzz_0_xyyyyyzz_1, g_yzzzzz_0_xyyyyzz_1, g_yzzzzz_0_xyyyyzzz_1, g_yzzzzz_0_xyyyzzz_1, g_yzzzzz_0_xyyyzzzz_1, g_yzzzzz_0_xyyzzzz_1, g_yzzzzz_0_xyyzzzzz_1, g_yzzzzz_0_xyzzzzz_1, g_yzzzzz_0_xyzzzzzz_1, g_yzzzzz_0_xzzzzzz_1, g_yzzzzz_0_xzzzzzzz_1, g_yzzzzz_0_yyyyyyyz_1, g_yzzzzz_0_yyyyyyz_1, g_yzzzzz_0_yyyyyyzz_1, g_yzzzzz_0_yyyyyzz_1, g_yzzzzz_0_yyyyyzzz_1, g_yzzzzz_0_yyyyzzz_1, g_yzzzzz_0_yyyyzzzz_1, g_yzzzzz_0_yyyzzzz_1, g_yzzzzz_0_yyyzzzzz_1, g_yzzzzz_0_yyzzzzz_1, g_yzzzzz_0_yyzzzzzz_1, g_yzzzzz_0_yzzzzzz_1, g_yzzzzz_0_yzzzzzzz_1, g_yzzzzz_0_zzzzzzz_1, g_yzzzzz_0_zzzzzzzz_1, g_zzzzz_0_xxxxxxxx_0, g_zzzzz_0_xxxxxxxx_1, g_zzzzz_0_xxxxxxxz_0, g_zzzzz_0_xxxxxxxz_1, g_zzzzz_0_xxxxxxyz_0, g_zzzzz_0_xxxxxxyz_1, g_zzzzz_0_xxxxxxzz_0, g_zzzzz_0_xxxxxxzz_1, g_zzzzz_0_xxxxxyyz_0, g_zzzzz_0_xxxxxyyz_1, g_zzzzz_0_xxxxxyzz_0, g_zzzzz_0_xxxxxyzz_1, g_zzzzz_0_xxxxxzzz_0, g_zzzzz_0_xxxxxzzz_1, g_zzzzz_0_xxxxyyyz_0, g_zzzzz_0_xxxxyyyz_1, g_zzzzz_0_xxxxyyzz_0, g_zzzzz_0_xxxxyyzz_1, g_zzzzz_0_xxxxyzzz_0, g_zzzzz_0_xxxxyzzz_1, g_zzzzz_0_xxxxzzzz_0, g_zzzzz_0_xxxxzzzz_1, g_zzzzz_0_xxxyyyyz_0, g_zzzzz_0_xxxyyyyz_1, g_zzzzz_0_xxxyyyzz_0, g_zzzzz_0_xxxyyyzz_1, g_zzzzz_0_xxxyyzzz_0, g_zzzzz_0_xxxyyzzz_1, g_zzzzz_0_xxxyzzzz_0, g_zzzzz_0_xxxyzzzz_1, g_zzzzz_0_xxxzzzzz_0, g_zzzzz_0_xxxzzzzz_1, g_zzzzz_0_xxyyyyyz_0, g_zzzzz_0_xxyyyyyz_1, g_zzzzz_0_xxyyyyzz_0, g_zzzzz_0_xxyyyyzz_1, g_zzzzz_0_xxyyyzzz_0, g_zzzzz_0_xxyyyzzz_1, g_zzzzz_0_xxyyzzzz_0, g_zzzzz_0_xxyyzzzz_1, g_zzzzz_0_xxyzzzzz_0, g_zzzzz_0_xxyzzzzz_1, g_zzzzz_0_xxzzzzzz_0, g_zzzzz_0_xxzzzzzz_1, g_zzzzz_0_xyyyyyyz_0, g_zzzzz_0_xyyyyyyz_1, g_zzzzz_0_xyyyyyzz_0, g_zzzzz_0_xyyyyyzz_1, g_zzzzz_0_xyyyyzzz_0, g_zzzzz_0_xyyyyzzz_1, g_zzzzz_0_xyyyzzzz_0, g_zzzzz_0_xyyyzzzz_1, g_zzzzz_0_xyyzzzzz_0, g_zzzzz_0_xyyzzzzz_1, g_zzzzz_0_xyzzzzzz_0, g_zzzzz_0_xyzzzzzz_1, g_zzzzz_0_xzzzzzzz_0, g_zzzzz_0_xzzzzzzz_1, g_zzzzz_0_yyyyyyyz_0, g_zzzzz_0_yyyyyyyz_1, g_zzzzz_0_yyyyyyzz_0, g_zzzzz_0_yyyyyyzz_1, g_zzzzz_0_yyyyyzzz_0, g_zzzzz_0_yyyyyzzz_1, g_zzzzz_0_yyyyzzzz_0, g_zzzzz_0_yyyyzzzz_1, g_zzzzz_0_yyyzzzzz_0, g_zzzzz_0_yyyzzzzz_1, g_zzzzz_0_yyzzzzzz_0, g_zzzzz_0_yyzzzzzz_1, g_zzzzz_0_yzzzzzzz_0, g_zzzzz_0_yzzzzzzz_1, g_zzzzz_0_zzzzzzzz_0, g_zzzzz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzzz_0_xxxxxxxx_0[i] = g_zzzzz_0_xxxxxxxx_0[i] * fbe_0 - g_zzzzz_0_xxxxxxxx_1[i] * fz_be_0 + g_yzzzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxxxxy_0[i] = 4.0 * g_yyzzz_0_xxxxxxxy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxxxxxy_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxxxxxz_0[i] = g_zzzzz_0_xxxxxxxz_0[i] * fbe_0 - g_zzzzz_0_xxxxxxxz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxxxyy_0[i] = 4.0 * g_yyzzz_0_xxxxxxyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxxxxyy_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxxxxyz_0[i] = g_zzzzz_0_xxxxxxyz_0[i] * fbe_0 - g_zzzzz_0_xxxxxxyz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxxxzz_0[i] = g_zzzzz_0_xxxxxxzz_0[i] * fbe_0 - g_zzzzz_0_xxxxxxzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxxyyy_0[i] = 4.0 * g_yyzzz_0_xxxxxyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxxxyyy_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxxxyyz_0[i] = g_zzzzz_0_xxxxxyyz_0[i] * fbe_0 - g_zzzzz_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxxyzz_0[i] = g_zzzzz_0_xxxxxyzz_0[i] * fbe_0 - g_zzzzz_0_xxxxxyzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxxzzz_0[i] = g_zzzzz_0_xxxxxzzz_0[i] * fbe_0 - g_zzzzz_0_xxxxxzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxyyyy_0[i] = 4.0 * g_yyzzz_0_xxxxyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxxyyyy_1[i] * fz_be_0 + g_yyzzzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxxyyyz_0[i] = g_zzzzz_0_xxxxyyyz_0[i] * fbe_0 - g_zzzzz_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxyyzz_0[i] = g_zzzzz_0_xxxxyyzz_0[i] * fbe_0 - g_zzzzz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxyzzz_0[i] = g_zzzzz_0_xxxxyzzz_0[i] * fbe_0 - g_zzzzz_0_xxxxyzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxzzzz_0[i] = g_zzzzz_0_xxxxzzzz_0[i] * fbe_0 - g_zzzzz_0_xxxxzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxyyyyy_0[i] = 4.0 * g_yyzzz_0_xxxyyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxyyyyy_1[i] * fz_be_0 + g_yyzzzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxyyyyz_0[i] = g_zzzzz_0_xxxyyyyz_0[i] * fbe_0 - g_zzzzz_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxyyyzz_0[i] = g_zzzzz_0_xxxyyyzz_0[i] * fbe_0 - g_zzzzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxyyzzz_0[i] = g_zzzzz_0_xxxyyzzz_0[i] * fbe_0 - g_zzzzz_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxyzzzz_0[i] = g_zzzzz_0_xxxyzzzz_0[i] * fbe_0 - g_zzzzz_0_xxxyzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxzzzzz_0[i] = g_zzzzz_0_xxxzzzzz_0[i] * fbe_0 - g_zzzzz_0_xxxzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyyyyyy_0[i] = 4.0 * g_yyzzz_0_xxyyyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxyyyyyy_1[i] * fz_be_0 + g_yyzzzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxyyyyyz_0[i] = g_zzzzz_0_xxyyyyyz_0[i] * fbe_0 - g_zzzzz_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yzzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyyyyzz_0[i] = g_zzzzz_0_xxyyyyzz_0[i] * fbe_0 - g_zzzzz_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyyyzzz_0[i] = g_zzzzz_0_xxyyyzzz_0[i] * fbe_0 - g_zzzzz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyyzzzz_0[i] = g_zzzzz_0_xxyyzzzz_0[i] * fbe_0 - g_zzzzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyzzzzz_0[i] = g_zzzzz_0_xxyzzzzz_0[i] * fbe_0 - g_zzzzz_0_xxyzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxzzzzzz_0[i] = g_zzzzz_0_xxzzzzzz_0[i] * fbe_0 - g_zzzzz_0_xxzzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyyyyyy_0[i] = 4.0 * g_yyzzz_0_xyyyyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xyyyyyyy_1[i] * fz_be_0 + g_yyzzzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xyyyyyyz_0[i] = g_zzzzz_0_xyyyyyyz_0[i] * fbe_0 - g_zzzzz_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yzzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyyyyzz_0[i] = g_zzzzz_0_xyyyyyzz_0[i] * fbe_0 - g_zzzzz_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yzzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyyyzzz_0[i] = g_zzzzz_0_xyyyyzzz_0[i] * fbe_0 - g_zzzzz_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyyzzzz_0[i] = g_zzzzz_0_xyyyzzzz_0[i] * fbe_0 - g_zzzzz_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyzzzzz_0[i] = g_zzzzz_0_xyyzzzzz_0[i] * fbe_0 - g_zzzzz_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyzzzzzz_0[i] = g_zzzzz_0_xyzzzzzz_0[i] * fbe_0 - g_zzzzz_0_xyzzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xzzzzzzz_0[i] = g_zzzzz_0_xzzzzzzz_0[i] * fbe_0 - g_zzzzz_0_xzzzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyyyyyy_0[i] = 4.0 * g_yyzzz_0_yyyyyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_yyyyyyyy_1[i] * fz_be_0 + g_yyzzzz_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_yyyyyyyz_0[i] = g_zzzzz_0_yyyyyyyz_0[i] * fbe_0 - g_zzzzz_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yzzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyyyyzz_0[i] = g_zzzzz_0_yyyyyyzz_0[i] * fbe_0 - g_zzzzz_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yzzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyyyzzz_0[i] = g_zzzzz_0_yyyyyzzz_0[i] * fbe_0 - g_zzzzz_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yzzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyyzzzz_0[i] = g_zzzzz_0_yyyyzzzz_0[i] * fbe_0 - g_zzzzz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyzzzzz_0[i] = g_zzzzz_0_yyyzzzzz_0[i] * fbe_0 - g_zzzzz_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyzzzzzz_0[i] = g_zzzzz_0_yyzzzzzz_0[i] * fbe_0 - g_zzzzz_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yzzzzzzz_0[i] = g_zzzzz_0_yzzzzzzz_0[i] * fbe_0 - g_zzzzz_0_yzzzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_zzzzzzzz_0[i] = g_zzzzz_0_zzzzzzzz_0[i] * fbe_0 - g_zzzzz_0_zzzzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1530-1575 components of targeted buffer : KSL

    auto g_yzzzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1530);

    auto g_yzzzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1531);

    auto g_yzzzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1532);

    auto g_yzzzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1533);

    auto g_yzzzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1534);

    auto g_yzzzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1535);

    auto g_yzzzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1536);

    auto g_yzzzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1537);

    auto g_yzzzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1538);

    auto g_yzzzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1539);

    auto g_yzzzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1540);

    auto g_yzzzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1541);

    auto g_yzzzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1542);

    auto g_yzzzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1543);

    auto g_yzzzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1544);

    auto g_yzzzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1545);

    auto g_yzzzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1546);

    auto g_yzzzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1547);

    auto g_yzzzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1548);

    auto g_yzzzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1549);

    auto g_yzzzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1550);

    auto g_yzzzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1551);

    auto g_yzzzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1552);

    auto g_yzzzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1553);

    auto g_yzzzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1554);

    auto g_yzzzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1555);

    auto g_yzzzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1556);

    auto g_yzzzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1557);

    auto g_yzzzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1558);

    auto g_yzzzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1559);

    auto g_yzzzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1560);

    auto g_yzzzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1561);

    auto g_yzzzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1562);

    auto g_yzzzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1563);

    auto g_yzzzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1564);

    auto g_yzzzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1565);

    auto g_yzzzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1566);

    auto g_yzzzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1567);

    auto g_yzzzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1568);

    auto g_yzzzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1569);

    auto g_yzzzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1570);

    auto g_yzzzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1571);

    auto g_yzzzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1572);

    auto g_yzzzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1573);

    auto g_yzzzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1574);

    #pragma omp simd aligned(g_yzzzzzz_0_xxxxxxxx_0, g_yzzzzzz_0_xxxxxxxy_0, g_yzzzzzz_0_xxxxxxxz_0, g_yzzzzzz_0_xxxxxxyy_0, g_yzzzzzz_0_xxxxxxyz_0, g_yzzzzzz_0_xxxxxxzz_0, g_yzzzzzz_0_xxxxxyyy_0, g_yzzzzzz_0_xxxxxyyz_0, g_yzzzzzz_0_xxxxxyzz_0, g_yzzzzzz_0_xxxxxzzz_0, g_yzzzzzz_0_xxxxyyyy_0, g_yzzzzzz_0_xxxxyyyz_0, g_yzzzzzz_0_xxxxyyzz_0, g_yzzzzzz_0_xxxxyzzz_0, g_yzzzzzz_0_xxxxzzzz_0, g_yzzzzzz_0_xxxyyyyy_0, g_yzzzzzz_0_xxxyyyyz_0, g_yzzzzzz_0_xxxyyyzz_0, g_yzzzzzz_0_xxxyyzzz_0, g_yzzzzzz_0_xxxyzzzz_0, g_yzzzzzz_0_xxxzzzzz_0, g_yzzzzzz_0_xxyyyyyy_0, g_yzzzzzz_0_xxyyyyyz_0, g_yzzzzzz_0_xxyyyyzz_0, g_yzzzzzz_0_xxyyyzzz_0, g_yzzzzzz_0_xxyyzzzz_0, g_yzzzzzz_0_xxyzzzzz_0, g_yzzzzzz_0_xxzzzzzz_0, g_yzzzzzz_0_xyyyyyyy_0, g_yzzzzzz_0_xyyyyyyz_0, g_yzzzzzz_0_xyyyyyzz_0, g_yzzzzzz_0_xyyyyzzz_0, g_yzzzzzz_0_xyyyzzzz_0, g_yzzzzzz_0_xyyzzzzz_0, g_yzzzzzz_0_xyzzzzzz_0, g_yzzzzzz_0_xzzzzzzz_0, g_yzzzzzz_0_yyyyyyyy_0, g_yzzzzzz_0_yyyyyyyz_0, g_yzzzzzz_0_yyyyyyzz_0, g_yzzzzzz_0_yyyyyzzz_0, g_yzzzzzz_0_yyyyzzzz_0, g_yzzzzzz_0_yyyzzzzz_0, g_yzzzzzz_0_yyzzzzzz_0, g_yzzzzzz_0_yzzzzzzz_0, g_yzzzzzz_0_zzzzzzzz_0, g_zzzzzz_0_xxxxxxx_1, g_zzzzzz_0_xxxxxxxx_1, g_zzzzzz_0_xxxxxxxy_1, g_zzzzzz_0_xxxxxxxz_1, g_zzzzzz_0_xxxxxxy_1, g_zzzzzz_0_xxxxxxyy_1, g_zzzzzz_0_xxxxxxyz_1, g_zzzzzz_0_xxxxxxz_1, g_zzzzzz_0_xxxxxxzz_1, g_zzzzzz_0_xxxxxyy_1, g_zzzzzz_0_xxxxxyyy_1, g_zzzzzz_0_xxxxxyyz_1, g_zzzzzz_0_xxxxxyz_1, g_zzzzzz_0_xxxxxyzz_1, g_zzzzzz_0_xxxxxzz_1, g_zzzzzz_0_xxxxxzzz_1, g_zzzzzz_0_xxxxyyy_1, g_zzzzzz_0_xxxxyyyy_1, g_zzzzzz_0_xxxxyyyz_1, g_zzzzzz_0_xxxxyyz_1, g_zzzzzz_0_xxxxyyzz_1, g_zzzzzz_0_xxxxyzz_1, g_zzzzzz_0_xxxxyzzz_1, g_zzzzzz_0_xxxxzzz_1, g_zzzzzz_0_xxxxzzzz_1, g_zzzzzz_0_xxxyyyy_1, g_zzzzzz_0_xxxyyyyy_1, g_zzzzzz_0_xxxyyyyz_1, g_zzzzzz_0_xxxyyyz_1, g_zzzzzz_0_xxxyyyzz_1, g_zzzzzz_0_xxxyyzz_1, g_zzzzzz_0_xxxyyzzz_1, g_zzzzzz_0_xxxyzzz_1, g_zzzzzz_0_xxxyzzzz_1, g_zzzzzz_0_xxxzzzz_1, g_zzzzzz_0_xxxzzzzz_1, g_zzzzzz_0_xxyyyyy_1, g_zzzzzz_0_xxyyyyyy_1, g_zzzzzz_0_xxyyyyyz_1, g_zzzzzz_0_xxyyyyz_1, g_zzzzzz_0_xxyyyyzz_1, g_zzzzzz_0_xxyyyzz_1, g_zzzzzz_0_xxyyyzzz_1, g_zzzzzz_0_xxyyzzz_1, g_zzzzzz_0_xxyyzzzz_1, g_zzzzzz_0_xxyzzzz_1, g_zzzzzz_0_xxyzzzzz_1, g_zzzzzz_0_xxzzzzz_1, g_zzzzzz_0_xxzzzzzz_1, g_zzzzzz_0_xyyyyyy_1, g_zzzzzz_0_xyyyyyyy_1, g_zzzzzz_0_xyyyyyyz_1, g_zzzzzz_0_xyyyyyz_1, g_zzzzzz_0_xyyyyyzz_1, g_zzzzzz_0_xyyyyzz_1, g_zzzzzz_0_xyyyyzzz_1, g_zzzzzz_0_xyyyzzz_1, g_zzzzzz_0_xyyyzzzz_1, g_zzzzzz_0_xyyzzzz_1, g_zzzzzz_0_xyyzzzzz_1, g_zzzzzz_0_xyzzzzz_1, g_zzzzzz_0_xyzzzzzz_1, g_zzzzzz_0_xzzzzzz_1, g_zzzzzz_0_xzzzzzzz_1, g_zzzzzz_0_yyyyyyy_1, g_zzzzzz_0_yyyyyyyy_1, g_zzzzzz_0_yyyyyyyz_1, g_zzzzzz_0_yyyyyyz_1, g_zzzzzz_0_yyyyyyzz_1, g_zzzzzz_0_yyyyyzz_1, g_zzzzzz_0_yyyyyzzz_1, g_zzzzzz_0_yyyyzzz_1, g_zzzzzz_0_yyyyzzzz_1, g_zzzzzz_0_yyyzzzz_1, g_zzzzzz_0_yyyzzzzz_1, g_zzzzzz_0_yyzzzzz_1, g_zzzzzz_0_yyzzzzzz_1, g_zzzzzz_0_yzzzzzz_1, g_zzzzzz_0_yzzzzzzz_1, g_zzzzzz_0_zzzzzzz_1, g_zzzzzz_0_zzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzzz_0_xxxxxxxx_0[i] = g_zzzzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxxxy_0[i] = g_zzzzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxxy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxxxz_0[i] = g_zzzzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxxyy_0[i] = 2.0 * g_zzzzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxxyz_0[i] = g_zzzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxxzz_0[i] = g_zzzzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxyyy_0[i] = 3.0 * g_zzzzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxyyz_0[i] = 2.0 * g_zzzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxyzz_0[i] = g_zzzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxzzz_0[i] = g_zzzzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxyyyy_0[i] = 4.0 * g_zzzzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxyyyz_0[i] = 3.0 * g_zzzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxyyzz_0[i] = 2.0 * g_zzzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxyzzz_0[i] = g_zzzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxzzzz_0[i] = g_zzzzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyyyyy_0[i] = 5.0 * g_zzzzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyyyyz_0[i] = 4.0 * g_zzzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyyyzz_0[i] = 3.0 * g_zzzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyyzzz_0[i] = 2.0 * g_zzzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyzzzz_0[i] = g_zzzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxzzzzz_0[i] = g_zzzzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyyyyy_0[i] = 6.0 * g_zzzzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyyyyz_0[i] = 5.0 * g_zzzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyyyzz_0[i] = 4.0 * g_zzzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyyzzz_0[i] = 3.0 * g_zzzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyzzzz_0[i] = 2.0 * g_zzzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyzzzzz_0[i] = g_zzzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxzzzzzz_0[i] = g_zzzzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyyyyy_0[i] = 7.0 * g_zzzzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyyyyz_0[i] = 6.0 * g_zzzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyyyzz_0[i] = 5.0 * g_zzzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyyzzz_0[i] = 4.0 * g_zzzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyzzzz_0[i] = 3.0 * g_zzzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyzzzzz_0[i] = 2.0 * g_zzzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyzzzzzz_0[i] = g_zzzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xzzzzzzz_0[i] = g_zzzzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyyyyy_0[i] = 8.0 * g_zzzzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyyyyz_0[i] = 7.0 * g_zzzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyyyzz_0[i] = 6.0 * g_zzzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyyzzz_0[i] = 5.0 * g_zzzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyzzzz_0[i] = 4.0 * g_zzzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyzzzzz_0[i] = 3.0 * g_zzzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyzzzzzz_0[i] = 2.0 * g_zzzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yzzzzzzz_0[i] = g_zzzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_zzzzzzzz_0[i] = g_zzzzzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1575-1620 components of targeted buffer : KSL

    auto g_zzzzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ksl + 1575);

    auto g_zzzzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_ksl + 1576);

    auto g_zzzzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ksl + 1577);

    auto g_zzzzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ksl + 1578);

    auto g_zzzzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_ksl + 1579);

    auto g_zzzzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ksl + 1580);

    auto g_zzzzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ksl + 1581);

    auto g_zzzzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_ksl + 1582);

    auto g_zzzzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_ksl + 1583);

    auto g_zzzzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ksl + 1584);

    auto g_zzzzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1585);

    auto g_zzzzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1586);

    auto g_zzzzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1587);

    auto g_zzzzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1588);

    auto g_zzzzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1589);

    auto g_zzzzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1590);

    auto g_zzzzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1591);

    auto g_zzzzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1592);

    auto g_zzzzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1593);

    auto g_zzzzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1594);

    auto g_zzzzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1595);

    auto g_zzzzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1596);

    auto g_zzzzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1597);

    auto g_zzzzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1598);

    auto g_zzzzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1599);

    auto g_zzzzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1600);

    auto g_zzzzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1601);

    auto g_zzzzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1602);

    auto g_zzzzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1603);

    auto g_zzzzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1604);

    auto g_zzzzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1605);

    auto g_zzzzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1606);

    auto g_zzzzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1607);

    auto g_zzzzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1608);

    auto g_zzzzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1609);

    auto g_zzzzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1610);

    auto g_zzzzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ksl + 1611);

    auto g_zzzzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ksl + 1612);

    auto g_zzzzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ksl + 1613);

    auto g_zzzzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ksl + 1614);

    auto g_zzzzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1615);

    auto g_zzzzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1616);

    auto g_zzzzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1617);

    auto g_zzzzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1618);

    auto g_zzzzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ksl + 1619);

    #pragma omp simd aligned(g_zzzzz_0_xxxxxxxx_0, g_zzzzz_0_xxxxxxxx_1, g_zzzzz_0_xxxxxxxy_0, g_zzzzz_0_xxxxxxxy_1, g_zzzzz_0_xxxxxxxz_0, g_zzzzz_0_xxxxxxxz_1, g_zzzzz_0_xxxxxxyy_0, g_zzzzz_0_xxxxxxyy_1, g_zzzzz_0_xxxxxxyz_0, g_zzzzz_0_xxxxxxyz_1, g_zzzzz_0_xxxxxxzz_0, g_zzzzz_0_xxxxxxzz_1, g_zzzzz_0_xxxxxyyy_0, g_zzzzz_0_xxxxxyyy_1, g_zzzzz_0_xxxxxyyz_0, g_zzzzz_0_xxxxxyyz_1, g_zzzzz_0_xxxxxyzz_0, g_zzzzz_0_xxxxxyzz_1, g_zzzzz_0_xxxxxzzz_0, g_zzzzz_0_xxxxxzzz_1, g_zzzzz_0_xxxxyyyy_0, g_zzzzz_0_xxxxyyyy_1, g_zzzzz_0_xxxxyyyz_0, g_zzzzz_0_xxxxyyyz_1, g_zzzzz_0_xxxxyyzz_0, g_zzzzz_0_xxxxyyzz_1, g_zzzzz_0_xxxxyzzz_0, g_zzzzz_0_xxxxyzzz_1, g_zzzzz_0_xxxxzzzz_0, g_zzzzz_0_xxxxzzzz_1, g_zzzzz_0_xxxyyyyy_0, g_zzzzz_0_xxxyyyyy_1, g_zzzzz_0_xxxyyyyz_0, g_zzzzz_0_xxxyyyyz_1, g_zzzzz_0_xxxyyyzz_0, g_zzzzz_0_xxxyyyzz_1, g_zzzzz_0_xxxyyzzz_0, g_zzzzz_0_xxxyyzzz_1, g_zzzzz_0_xxxyzzzz_0, g_zzzzz_0_xxxyzzzz_1, g_zzzzz_0_xxxzzzzz_0, g_zzzzz_0_xxxzzzzz_1, g_zzzzz_0_xxyyyyyy_0, g_zzzzz_0_xxyyyyyy_1, g_zzzzz_0_xxyyyyyz_0, g_zzzzz_0_xxyyyyyz_1, g_zzzzz_0_xxyyyyzz_0, g_zzzzz_0_xxyyyyzz_1, g_zzzzz_0_xxyyyzzz_0, g_zzzzz_0_xxyyyzzz_1, g_zzzzz_0_xxyyzzzz_0, g_zzzzz_0_xxyyzzzz_1, g_zzzzz_0_xxyzzzzz_0, g_zzzzz_0_xxyzzzzz_1, g_zzzzz_0_xxzzzzzz_0, g_zzzzz_0_xxzzzzzz_1, g_zzzzz_0_xyyyyyyy_0, g_zzzzz_0_xyyyyyyy_1, g_zzzzz_0_xyyyyyyz_0, g_zzzzz_0_xyyyyyyz_1, g_zzzzz_0_xyyyyyzz_0, g_zzzzz_0_xyyyyyzz_1, g_zzzzz_0_xyyyyzzz_0, g_zzzzz_0_xyyyyzzz_1, g_zzzzz_0_xyyyzzzz_0, g_zzzzz_0_xyyyzzzz_1, g_zzzzz_0_xyyzzzzz_0, g_zzzzz_0_xyyzzzzz_1, g_zzzzz_0_xyzzzzzz_0, g_zzzzz_0_xyzzzzzz_1, g_zzzzz_0_xzzzzzzz_0, g_zzzzz_0_xzzzzzzz_1, g_zzzzz_0_yyyyyyyy_0, g_zzzzz_0_yyyyyyyy_1, g_zzzzz_0_yyyyyyyz_0, g_zzzzz_0_yyyyyyyz_1, g_zzzzz_0_yyyyyyzz_0, g_zzzzz_0_yyyyyyzz_1, g_zzzzz_0_yyyyyzzz_0, g_zzzzz_0_yyyyyzzz_1, g_zzzzz_0_yyyyzzzz_0, g_zzzzz_0_yyyyzzzz_1, g_zzzzz_0_yyyzzzzz_0, g_zzzzz_0_yyyzzzzz_1, g_zzzzz_0_yyzzzzzz_0, g_zzzzz_0_yyzzzzzz_1, g_zzzzz_0_yzzzzzzz_0, g_zzzzz_0_yzzzzzzz_1, g_zzzzz_0_zzzzzzzz_0, g_zzzzz_0_zzzzzzzz_1, g_zzzzzz_0_xxxxxxx_1, g_zzzzzz_0_xxxxxxxx_1, g_zzzzzz_0_xxxxxxxy_1, g_zzzzzz_0_xxxxxxxz_1, g_zzzzzz_0_xxxxxxy_1, g_zzzzzz_0_xxxxxxyy_1, g_zzzzzz_0_xxxxxxyz_1, g_zzzzzz_0_xxxxxxz_1, g_zzzzzz_0_xxxxxxzz_1, g_zzzzzz_0_xxxxxyy_1, g_zzzzzz_0_xxxxxyyy_1, g_zzzzzz_0_xxxxxyyz_1, g_zzzzzz_0_xxxxxyz_1, g_zzzzzz_0_xxxxxyzz_1, g_zzzzzz_0_xxxxxzz_1, g_zzzzzz_0_xxxxxzzz_1, g_zzzzzz_0_xxxxyyy_1, g_zzzzzz_0_xxxxyyyy_1, g_zzzzzz_0_xxxxyyyz_1, g_zzzzzz_0_xxxxyyz_1, g_zzzzzz_0_xxxxyyzz_1, g_zzzzzz_0_xxxxyzz_1, g_zzzzzz_0_xxxxyzzz_1, g_zzzzzz_0_xxxxzzz_1, g_zzzzzz_0_xxxxzzzz_1, g_zzzzzz_0_xxxyyyy_1, g_zzzzzz_0_xxxyyyyy_1, g_zzzzzz_0_xxxyyyyz_1, g_zzzzzz_0_xxxyyyz_1, g_zzzzzz_0_xxxyyyzz_1, g_zzzzzz_0_xxxyyzz_1, g_zzzzzz_0_xxxyyzzz_1, g_zzzzzz_0_xxxyzzz_1, g_zzzzzz_0_xxxyzzzz_1, g_zzzzzz_0_xxxzzzz_1, g_zzzzzz_0_xxxzzzzz_1, g_zzzzzz_0_xxyyyyy_1, g_zzzzzz_0_xxyyyyyy_1, g_zzzzzz_0_xxyyyyyz_1, g_zzzzzz_0_xxyyyyz_1, g_zzzzzz_0_xxyyyyzz_1, g_zzzzzz_0_xxyyyzz_1, g_zzzzzz_0_xxyyyzzz_1, g_zzzzzz_0_xxyyzzz_1, g_zzzzzz_0_xxyyzzzz_1, g_zzzzzz_0_xxyzzzz_1, g_zzzzzz_0_xxyzzzzz_1, g_zzzzzz_0_xxzzzzz_1, g_zzzzzz_0_xxzzzzzz_1, g_zzzzzz_0_xyyyyyy_1, g_zzzzzz_0_xyyyyyyy_1, g_zzzzzz_0_xyyyyyyz_1, g_zzzzzz_0_xyyyyyz_1, g_zzzzzz_0_xyyyyyzz_1, g_zzzzzz_0_xyyyyzz_1, g_zzzzzz_0_xyyyyzzz_1, g_zzzzzz_0_xyyyzzz_1, g_zzzzzz_0_xyyyzzzz_1, g_zzzzzz_0_xyyzzzz_1, g_zzzzzz_0_xyyzzzzz_1, g_zzzzzz_0_xyzzzzz_1, g_zzzzzz_0_xyzzzzzz_1, g_zzzzzz_0_xzzzzzz_1, g_zzzzzz_0_xzzzzzzz_1, g_zzzzzz_0_yyyyyyy_1, g_zzzzzz_0_yyyyyyyy_1, g_zzzzzz_0_yyyyyyyz_1, g_zzzzzz_0_yyyyyyz_1, g_zzzzzz_0_yyyyyyzz_1, g_zzzzzz_0_yyyyyzz_1, g_zzzzzz_0_yyyyyzzz_1, g_zzzzzz_0_yyyyzzz_1, g_zzzzzz_0_yyyyzzzz_1, g_zzzzzz_0_yyyzzzz_1, g_zzzzzz_0_yyyzzzzz_1, g_zzzzzz_0_yyzzzzz_1, g_zzzzzz_0_yyzzzzzz_1, g_zzzzzz_0_yzzzzzz_1, g_zzzzzz_0_yzzzzzzz_1, g_zzzzzz_0_zzzzzzz_1, g_zzzzzz_0_zzzzzzzz_1, g_zzzzzzz_0_xxxxxxxx_0, g_zzzzzzz_0_xxxxxxxy_0, g_zzzzzzz_0_xxxxxxxz_0, g_zzzzzzz_0_xxxxxxyy_0, g_zzzzzzz_0_xxxxxxyz_0, g_zzzzzzz_0_xxxxxxzz_0, g_zzzzzzz_0_xxxxxyyy_0, g_zzzzzzz_0_xxxxxyyz_0, g_zzzzzzz_0_xxxxxyzz_0, g_zzzzzzz_0_xxxxxzzz_0, g_zzzzzzz_0_xxxxyyyy_0, g_zzzzzzz_0_xxxxyyyz_0, g_zzzzzzz_0_xxxxyyzz_0, g_zzzzzzz_0_xxxxyzzz_0, g_zzzzzzz_0_xxxxzzzz_0, g_zzzzzzz_0_xxxyyyyy_0, g_zzzzzzz_0_xxxyyyyz_0, g_zzzzzzz_0_xxxyyyzz_0, g_zzzzzzz_0_xxxyyzzz_0, g_zzzzzzz_0_xxxyzzzz_0, g_zzzzzzz_0_xxxzzzzz_0, g_zzzzzzz_0_xxyyyyyy_0, g_zzzzzzz_0_xxyyyyyz_0, g_zzzzzzz_0_xxyyyyzz_0, g_zzzzzzz_0_xxyyyzzz_0, g_zzzzzzz_0_xxyyzzzz_0, g_zzzzzzz_0_xxyzzzzz_0, g_zzzzzzz_0_xxzzzzzz_0, g_zzzzzzz_0_xyyyyyyy_0, g_zzzzzzz_0_xyyyyyyz_0, g_zzzzzzz_0_xyyyyyzz_0, g_zzzzzzz_0_xyyyyzzz_0, g_zzzzzzz_0_xyyyzzzz_0, g_zzzzzzz_0_xyyzzzzz_0, g_zzzzzzz_0_xyzzzzzz_0, g_zzzzzzz_0_xzzzzzzz_0, g_zzzzzzz_0_yyyyyyyy_0, g_zzzzzzz_0_yyyyyyyz_0, g_zzzzzzz_0_yyyyyyzz_0, g_zzzzzzz_0_yyyyyzzz_0, g_zzzzzzz_0_yyyyzzzz_0, g_zzzzzzz_0_yyyzzzzz_0, g_zzzzzzz_0_yyzzzzzz_0, g_zzzzzzz_0_yzzzzzzz_0, g_zzzzzzz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzzz_0_xxxxxxxx_0[i] = 6.0 * g_zzzzz_0_xxxxxxxx_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxxxx_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxxxx_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxxxy_0[i] = 6.0 * g_zzzzz_0_xxxxxxxy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxxxy_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxxxz_0[i] = 6.0 * g_zzzzz_0_xxxxxxxz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxxxz_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxxz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxxyy_0[i] = 6.0 * g_zzzzz_0_xxxxxxyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxxyy_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxxyz_0[i] = 6.0 * g_zzzzz_0_xxxxxxyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxxyz_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxxzz_0[i] = 6.0 * g_zzzzz_0_xxxxxxzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxxzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxyyy_0[i] = 6.0 * g_zzzzz_0_xxxxxyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxyyy_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxyyz_0[i] = 6.0 * g_zzzzz_0_xxxxxyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxyyz_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxyzz_0[i] = 6.0 * g_zzzzz_0_xxxxxyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxzzz_0[i] = 6.0 * g_zzzzz_0_xxxxxzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxyyyy_0[i] = 6.0 * g_zzzzz_0_xxxxyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxyyyy_1[i] * fz_be_0 + g_zzzzzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxyyyz_0[i] = 6.0 * g_zzzzz_0_xxxxyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxyyyz_1[i] * fz_be_0 + g_zzzzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxyyzz_0[i] = 6.0 * g_zzzzz_0_xxxxyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxyzzz_0[i] = 6.0 * g_zzzzz_0_xxxxyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxzzzz_0[i] = 6.0 * g_zzzzz_0_xxxxzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyyyyy_0[i] = 6.0 * g_zzzzz_0_xxxyyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyyyyy_1[i] * fz_be_0 + g_zzzzzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyyyyz_0[i] = 6.0 * g_zzzzz_0_xxxyyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyyyyz_1[i] * fz_be_0 + g_zzzzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyyyzz_0[i] = 6.0 * g_zzzzz_0_xxxyyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyyzzz_0[i] = 6.0 * g_zzzzz_0_xxxyyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyzzzz_0[i] = 6.0 * g_zzzzz_0_xxxyzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxzzzzz_0[i] = 6.0 * g_zzzzz_0_xxxzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyyyyy_0[i] = 6.0 * g_zzzzz_0_xxyyyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyyyyy_1[i] * fz_be_0 + g_zzzzzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyyyyz_0[i] = 6.0 * g_zzzzz_0_xxyyyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyyyyz_1[i] * fz_be_0 + g_zzzzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyyyzz_0[i] = 6.0 * g_zzzzz_0_xxyyyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyyzzz_0[i] = 6.0 * g_zzzzz_0_xxyyyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyzzzz_0[i] = 6.0 * g_zzzzz_0_xxyyzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyzzzzz_0[i] = 6.0 * g_zzzzz_0_xxyzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxzzzzzz_0[i] = 6.0 * g_zzzzz_0_xxzzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxzzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyyyyy_0[i] = 6.0 * g_zzzzz_0_xyyyyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyyyyy_1[i] * fz_be_0 + g_zzzzzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyyyyz_0[i] = 6.0 * g_zzzzz_0_xyyyyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyyyyz_1[i] * fz_be_0 + g_zzzzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyyyzz_0[i] = 6.0 * g_zzzzz_0_xyyyyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyyzzz_0[i] = 6.0 * g_zzzzz_0_xyyyyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyzzzz_0[i] = 6.0 * g_zzzzz_0_xyyyzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyzzzzz_0[i] = 6.0 * g_zzzzz_0_xyyzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyzzzzzz_0[i] = 6.0 * g_zzzzz_0_xyzzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xzzzzzzz_0[i] = 6.0 * g_zzzzz_0_xzzzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzzzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xzzzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyyyyy_0[i] = 6.0 * g_zzzzz_0_yyyyyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyyyyy_1[i] * fz_be_0 + g_zzzzzz_0_yyyyyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyyyyz_0[i] = 6.0 * g_zzzzz_0_yyyyyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyyyyz_1[i] * fz_be_0 + g_zzzzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyyyzz_0[i] = 6.0 * g_zzzzz_0_yyyyyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyyzzz_0[i] = 6.0 * g_zzzzz_0_yyyyyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyzzzz_0[i] = 6.0 * g_zzzzz_0_yyyyzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyzzzzz_0[i] = 6.0 * g_zzzzz_0_yyyzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyzzzzzz_0[i] = 6.0 * g_zzzzz_0_yyzzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyzzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yzzzzzzz_0[i] = 6.0 * g_zzzzz_0_yzzzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzzzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yzzzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_zzzzzzzz_0[i] = 6.0 * g_zzzzz_0_zzzzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_zzzzzzzz_1[i] * fz_be_0 + 8.0 * g_zzzzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_zzzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

