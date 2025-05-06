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

#include "ThreeCenterElectronRepulsionPrimRecISI.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_isi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_isi,
                                 size_t idx_eri_0_gsi,
                                 size_t idx_eri_1_gsi,
                                 size_t idx_eri_1_hsh,
                                 size_t idx_eri_1_hsi,
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

    /// Set up components of auxilary buffer : GSI

    auto g_xxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi);

    auto g_xxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 1);

    auto g_xxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 2);

    auto g_xxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 3);

    auto g_xxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 4);

    auto g_xxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 5);

    auto g_xxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 6);

    auto g_xxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 7);

    auto g_xxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 8);

    auto g_xxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 9);

    auto g_xxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 10);

    auto g_xxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 11);

    auto g_xxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 12);

    auto g_xxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 13);

    auto g_xxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 14);

    auto g_xxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 15);

    auto g_xxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 16);

    auto g_xxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 17);

    auto g_xxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 18);

    auto g_xxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 19);

    auto g_xxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 20);

    auto g_xxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 21);

    auto g_xxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 22);

    auto g_xxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 23);

    auto g_xxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 24);

    auto g_xxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 25);

    auto g_xxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 26);

    auto g_xxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 27);

    auto g_xxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 28);

    auto g_xxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 30);

    auto g_xxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 33);

    auto g_xxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 37);

    auto g_xxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 42);

    auto g_xxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 48);

    auto g_xxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 56);

    auto g_xxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 57);

    auto g_xxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 59);

    auto g_xxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 62);

    auto g_xxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 66);

    auto g_xxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 71);

    auto g_xxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 84);

    auto g_xxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 85);

    auto g_xxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 86);

    auto g_xxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 87);

    auto g_xxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 88);

    auto g_xxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 89);

    auto g_xxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 90);

    auto g_xxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 91);

    auto g_xxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 92);

    auto g_xxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 93);

    auto g_xxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 94);

    auto g_xxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 95);

    auto g_xxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 96);

    auto g_xxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 97);

    auto g_xxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 98);

    auto g_xxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 99);

    auto g_xxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 100);

    auto g_xxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 101);

    auto g_xxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 102);

    auto g_xxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 103);

    auto g_xxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 104);

    auto g_xxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 105);

    auto g_xxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 106);

    auto g_xxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 107);

    auto g_xxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 108);

    auto g_xxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 109);

    auto g_xxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 110);

    auto g_xxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 111);

    auto g_xxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 140);

    auto g_xxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 141);

    auto g_xxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 142);

    auto g_xxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 143);

    auto g_xxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 144);

    auto g_xxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 145);

    auto g_xxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 146);

    auto g_xxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 147);

    auto g_xxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 148);

    auto g_xxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 149);

    auto g_xxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 150);

    auto g_xxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 151);

    auto g_xxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 152);

    auto g_xxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 153);

    auto g_xxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 154);

    auto g_xxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 155);

    auto g_xxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 156);

    auto g_xxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 157);

    auto g_xxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 158);

    auto g_xxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 159);

    auto g_xxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 160);

    auto g_xxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 161);

    auto g_xxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 162);

    auto g_xxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 163);

    auto g_xxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 164);

    auto g_xxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 165);

    auto g_xxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 166);

    auto g_xxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 167);

    auto g_xyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 169);

    auto g_xyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 171);

    auto g_xyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 172);

    auto g_xyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 174);

    auto g_xyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 175);

    auto g_xyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 176);

    auto g_xyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 178);

    auto g_xyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 179);

    auto g_xyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 180);

    auto g_xyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 181);

    auto g_xyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 183);

    auto g_xyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 184);

    auto g_xyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 185);

    auto g_xyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 186);

    auto g_xyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 187);

    auto g_xyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 189);

    auto g_xyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 190);

    auto g_xyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 191);

    auto g_xyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 192);

    auto g_xyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 193);

    auto g_xyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 194);

    auto g_xyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 195);

    auto g_xzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 254);

    auto g_xzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 256);

    auto g_xzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 257);

    auto g_xzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 259);

    auto g_xzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 260);

    auto g_xzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 261);

    auto g_xzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 263);

    auto g_xzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 264);

    auto g_xzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 265);

    auto g_xzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 266);

    auto g_xzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 268);

    auto g_xzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 269);

    auto g_xzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 270);

    auto g_xzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 271);

    auto g_xzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 272);

    auto g_xzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 273);

    auto g_xzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 274);

    auto g_xzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 275);

    auto g_xzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 276);

    auto g_xzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 277);

    auto g_xzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 278);

    auto g_xzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 279);

    auto g_yyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 280);

    auto g_yyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 281);

    auto g_yyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 282);

    auto g_yyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 283);

    auto g_yyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 284);

    auto g_yyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 285);

    auto g_yyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 286);

    auto g_yyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 287);

    auto g_yyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 288);

    auto g_yyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 289);

    auto g_yyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 290);

    auto g_yyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 291);

    auto g_yyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 292);

    auto g_yyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 293);

    auto g_yyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 294);

    auto g_yyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 295);

    auto g_yyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 296);

    auto g_yyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 297);

    auto g_yyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 298);

    auto g_yyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 299);

    auto g_yyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 300);

    auto g_yyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 301);

    auto g_yyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 302);

    auto g_yyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 303);

    auto g_yyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 304);

    auto g_yyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 305);

    auto g_yyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 306);

    auto g_yyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 307);

    auto g_yyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 309);

    auto g_yyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 311);

    auto g_yyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 314);

    auto g_yyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 318);

    auto g_yyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 323);

    auto g_yyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 329);

    auto g_yyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 336);

    auto g_yyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 337);

    auto g_yyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 338);

    auto g_yyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 339);

    auto g_yyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 340);

    auto g_yyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 341);

    auto g_yyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 342);

    auto g_yyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 343);

    auto g_yyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 344);

    auto g_yyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 345);

    auto g_yyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 346);

    auto g_yyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 347);

    auto g_yyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 348);

    auto g_yyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 349);

    auto g_yyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 350);

    auto g_yyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 351);

    auto g_yyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 352);

    auto g_yyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 353);

    auto g_yyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 354);

    auto g_yyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 355);

    auto g_yyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 356);

    auto g_yyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 357);

    auto g_yyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 358);

    auto g_yyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 359);

    auto g_yyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 360);

    auto g_yyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 361);

    auto g_yyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 362);

    auto g_yyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 363);

    auto g_yzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 364);

    auto g_yzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 366);

    auto g_yzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 368);

    auto g_yzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 369);

    auto g_yzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 371);

    auto g_yzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 372);

    auto g_yzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 373);

    auto g_yzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 375);

    auto g_yzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 376);

    auto g_yzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 377);

    auto g_yzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 378);

    auto g_yzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 380);

    auto g_yzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 381);

    auto g_yzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 382);

    auto g_yzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 383);

    auto g_yzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 384);

    auto g_yzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 386);

    auto g_yzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 387);

    auto g_yzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 388);

    auto g_yzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 389);

    auto g_yzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 390);

    auto g_yzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 391);

    auto g_zzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 392);

    auto g_zzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 393);

    auto g_zzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 394);

    auto g_zzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 395);

    auto g_zzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 396);

    auto g_zzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 397);

    auto g_zzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 398);

    auto g_zzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 399);

    auto g_zzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 400);

    auto g_zzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 401);

    auto g_zzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 402);

    auto g_zzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 403);

    auto g_zzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 404);

    auto g_zzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 405);

    auto g_zzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 406);

    auto g_zzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 407);

    auto g_zzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 408);

    auto g_zzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 409);

    auto g_zzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 410);

    auto g_zzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 411);

    auto g_zzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 412);

    auto g_zzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 413);

    auto g_zzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 414);

    auto g_zzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 415);

    auto g_zzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 416);

    auto g_zzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 417);

    auto g_zzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 418);

    auto g_zzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 419);

    /// Set up components of auxilary buffer : GSI

    auto g_xxxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi);

    auto g_xxxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 1);

    auto g_xxxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 2);

    auto g_xxxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 3);

    auto g_xxxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 4);

    auto g_xxxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 5);

    auto g_xxxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 6);

    auto g_xxxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 7);

    auto g_xxxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 8);

    auto g_xxxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 9);

    auto g_xxxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 10);

    auto g_xxxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 11);

    auto g_xxxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 12);

    auto g_xxxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 13);

    auto g_xxxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 14);

    auto g_xxxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 15);

    auto g_xxxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 16);

    auto g_xxxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 17);

    auto g_xxxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 18);

    auto g_xxxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 19);

    auto g_xxxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 20);

    auto g_xxxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 21);

    auto g_xxxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 22);

    auto g_xxxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 23);

    auto g_xxxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 24);

    auto g_xxxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 25);

    auto g_xxxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 26);

    auto g_xxxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 27);

    auto g_xxxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 28);

    auto g_xxxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 30);

    auto g_xxxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 33);

    auto g_xxxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 37);

    auto g_xxxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 42);

    auto g_xxxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 48);

    auto g_xxxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 56);

    auto g_xxxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 57);

    auto g_xxxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 59);

    auto g_xxxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 62);

    auto g_xxxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 66);

    auto g_xxxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 71);

    auto g_xxyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 84);

    auto g_xxyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 85);

    auto g_xxyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 86);

    auto g_xxyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 87);

    auto g_xxyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 88);

    auto g_xxyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 89);

    auto g_xxyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 90);

    auto g_xxyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 91);

    auto g_xxyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 92);

    auto g_xxyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 93);

    auto g_xxyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 94);

    auto g_xxyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 95);

    auto g_xxyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 96);

    auto g_xxyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 97);

    auto g_xxyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 98);

    auto g_xxyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 99);

    auto g_xxyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 100);

    auto g_xxyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 101);

    auto g_xxyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 102);

    auto g_xxyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 103);

    auto g_xxyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 104);

    auto g_xxyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 105);

    auto g_xxyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 106);

    auto g_xxyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 107);

    auto g_xxyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 108);

    auto g_xxyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 109);

    auto g_xxyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 110);

    auto g_xxyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 111);

    auto g_xxzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 140);

    auto g_xxzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 141);

    auto g_xxzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 142);

    auto g_xxzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 143);

    auto g_xxzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 144);

    auto g_xxzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 145);

    auto g_xxzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 146);

    auto g_xxzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 147);

    auto g_xxzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 148);

    auto g_xxzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 149);

    auto g_xxzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 150);

    auto g_xxzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 151);

    auto g_xxzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 152);

    auto g_xxzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 153);

    auto g_xxzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 154);

    auto g_xxzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 155);

    auto g_xxzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 156);

    auto g_xxzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 157);

    auto g_xxzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 158);

    auto g_xxzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 159);

    auto g_xxzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 160);

    auto g_xxzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 161);

    auto g_xxzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 162);

    auto g_xxzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 163);

    auto g_xxzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 164);

    auto g_xxzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 165);

    auto g_xxzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 166);

    auto g_xxzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 167);

    auto g_xyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 169);

    auto g_xyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 171);

    auto g_xyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 172);

    auto g_xyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 174);

    auto g_xyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 175);

    auto g_xyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 176);

    auto g_xyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 178);

    auto g_xyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 179);

    auto g_xyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 180);

    auto g_xyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 181);

    auto g_xyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 183);

    auto g_xyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 184);

    auto g_xyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 185);

    auto g_xyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 186);

    auto g_xyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 187);

    auto g_xyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 189);

    auto g_xyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 190);

    auto g_xyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 191);

    auto g_xyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 192);

    auto g_xyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 193);

    auto g_xyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 194);

    auto g_xyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 195);

    auto g_xzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 254);

    auto g_xzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 256);

    auto g_xzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 257);

    auto g_xzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 259);

    auto g_xzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 260);

    auto g_xzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 261);

    auto g_xzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 263);

    auto g_xzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 264);

    auto g_xzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 265);

    auto g_xzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 266);

    auto g_xzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 268);

    auto g_xzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 269);

    auto g_xzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 270);

    auto g_xzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 271);

    auto g_xzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 272);

    auto g_xzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 273);

    auto g_xzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 274);

    auto g_xzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 275);

    auto g_xzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 276);

    auto g_xzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 277);

    auto g_xzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 278);

    auto g_xzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 279);

    auto g_yyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 280);

    auto g_yyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 281);

    auto g_yyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 282);

    auto g_yyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 283);

    auto g_yyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 284);

    auto g_yyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 285);

    auto g_yyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 286);

    auto g_yyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 287);

    auto g_yyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 288);

    auto g_yyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 289);

    auto g_yyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 290);

    auto g_yyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 291);

    auto g_yyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 292);

    auto g_yyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 293);

    auto g_yyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 294);

    auto g_yyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 295);

    auto g_yyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 296);

    auto g_yyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 297);

    auto g_yyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 298);

    auto g_yyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 299);

    auto g_yyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 300);

    auto g_yyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 301);

    auto g_yyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 302);

    auto g_yyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 303);

    auto g_yyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 304);

    auto g_yyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 305);

    auto g_yyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 306);

    auto g_yyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 307);

    auto g_yyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 309);

    auto g_yyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 311);

    auto g_yyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 314);

    auto g_yyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 318);

    auto g_yyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 323);

    auto g_yyyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 329);

    auto g_yyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 336);

    auto g_yyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 337);

    auto g_yyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 338);

    auto g_yyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 339);

    auto g_yyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 340);

    auto g_yyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 341);

    auto g_yyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 342);

    auto g_yyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 343);

    auto g_yyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 344);

    auto g_yyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 345);

    auto g_yyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 346);

    auto g_yyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 347);

    auto g_yyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 348);

    auto g_yyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 349);

    auto g_yyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 350);

    auto g_yyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 351);

    auto g_yyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 352);

    auto g_yyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 353);

    auto g_yyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 354);

    auto g_yyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 355);

    auto g_yyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 356);

    auto g_yyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 357);

    auto g_yyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 358);

    auto g_yyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 359);

    auto g_yyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 360);

    auto g_yyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 361);

    auto g_yyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 362);

    auto g_yyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 363);

    auto g_yzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 364);

    auto g_yzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 366);

    auto g_yzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 368);

    auto g_yzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 369);

    auto g_yzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 371);

    auto g_yzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 372);

    auto g_yzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 373);

    auto g_yzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 375);

    auto g_yzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 376);

    auto g_yzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 377);

    auto g_yzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 378);

    auto g_yzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 380);

    auto g_yzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 381);

    auto g_yzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 382);

    auto g_yzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 383);

    auto g_yzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 384);

    auto g_yzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 386);

    auto g_yzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 387);

    auto g_yzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 388);

    auto g_yzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 389);

    auto g_yzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 390);

    auto g_yzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 391);

    auto g_zzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_gsi + 392);

    auto g_zzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_gsi + 393);

    auto g_zzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_gsi + 394);

    auto g_zzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_gsi + 395);

    auto g_zzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_gsi + 396);

    auto g_zzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_gsi + 397);

    auto g_zzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_gsi + 398);

    auto g_zzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_gsi + 399);

    auto g_zzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_gsi + 400);

    auto g_zzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_gsi + 401);

    auto g_zzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_gsi + 402);

    auto g_zzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_gsi + 403);

    auto g_zzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_gsi + 404);

    auto g_zzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_gsi + 405);

    auto g_zzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_gsi + 406);

    auto g_zzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 407);

    auto g_zzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 408);

    auto g_zzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 409);

    auto g_zzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 410);

    auto g_zzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 411);

    auto g_zzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 412);

    auto g_zzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_gsi + 413);

    auto g_zzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_gsi + 414);

    auto g_zzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_gsi + 415);

    auto g_zzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_gsi + 416);

    auto g_zzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_gsi + 417);

    auto g_zzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 418);

    auto g_zzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_gsi + 419);

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

    auto g_xxxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 44);

    auto g_xxxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 46);

    auto g_xxxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 47);

    auto g_xxxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 49);

    auto g_xxxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 50);

    auto g_xxxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 51);

    auto g_xxxxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 53);

    auto g_xxxxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 54);

    auto g_xxxxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 55);

    auto g_xxxxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 56);

    auto g_xxxxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 58);

    auto g_xxxxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 59);

    auto g_xxxxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 60);

    auto g_xxxxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 61);

    auto g_xxxxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 62);

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

    auto g_xyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 256);

    auto g_xyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 259);

    auto g_xyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 260);

    auto g_xyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 263);

    auto g_xyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 264);

    auto g_xyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 265);

    auto g_xyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 268);

    auto g_xyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 269);

    auto g_xyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 270);

    auto g_xyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 271);

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

    auto g_yyyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 338);

    auto g_yyyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 340);

    auto g_yyyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 341);

    auto g_yyyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 343);

    auto g_yyyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 344);

    auto g_yyyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 345);

    auto g_yyyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 347);

    auto g_yyyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 348);

    auto g_yyyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 349);

    auto g_yyyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 350);

    auto g_yyyyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_hsh + 352);

    auto g_yyyyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_hsh + 353);

    auto g_yyyyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_hsh + 354);

    auto g_yyyyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_hsh + 355);

    auto g_yyyyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_hsh + 356);

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

    auto g_yzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_hsh + 400);

    auto g_yzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_hsh + 401);

    auto g_yzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_hsh + 402);

    auto g_yzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_hsh + 403);

    auto g_yzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_hsh + 404);

    auto g_yzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_hsh + 405);

    auto g_yzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_hsh + 406);

    auto g_yzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_hsh + 407);

    auto g_yzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_hsh + 408);

    auto g_yzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_hsh + 409);

    auto g_yzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_hsh + 410);

    auto g_yzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_hsh + 411);

    auto g_yzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_hsh + 412);

    auto g_yzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_hsh + 413);

    auto g_yzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_hsh + 414);

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

    /// Set up components of auxilary buffer : HSI

    auto g_xxxxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi);

    auto g_xxxxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 1);

    auto g_xxxxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 2);

    auto g_xxxxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 3);

    auto g_xxxxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 4);

    auto g_xxxxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 5);

    auto g_xxxxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 6);

    auto g_xxxxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 7);

    auto g_xxxxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 8);

    auto g_xxxxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 9);

    auto g_xxxxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 10);

    auto g_xxxxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 11);

    auto g_xxxxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 12);

    auto g_xxxxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 13);

    auto g_xxxxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 14);

    auto g_xxxxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 15);

    auto g_xxxxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 16);

    auto g_xxxxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 17);

    auto g_xxxxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 18);

    auto g_xxxxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 19);

    auto g_xxxxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 20);

    auto g_xxxxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 21);

    auto g_xxxxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 22);

    auto g_xxxxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 23);

    auto g_xxxxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 24);

    auto g_xxxxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 25);

    auto g_xxxxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 26);

    auto g_xxxxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 27);

    auto g_xxxxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 28);

    auto g_xxxxy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 29);

    auto g_xxxxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 30);

    auto g_xxxxy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 31);

    auto g_xxxxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 33);

    auto g_xxxxy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 34);

    auto g_xxxxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 37);

    auto g_xxxxy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 38);

    auto g_xxxxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 42);

    auto g_xxxxy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 43);

    auto g_xxxxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 48);

    auto g_xxxxy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 49);

    auto g_xxxxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 56);

    auto g_xxxxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 57);

    auto g_xxxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 58);

    auto g_xxxxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 59);

    auto g_xxxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 60);

    auto g_xxxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 61);

    auto g_xxxxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 62);

    auto g_xxxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 63);

    auto g_xxxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 64);

    auto g_xxxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 65);

    auto g_xxxxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 66);

    auto g_xxxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 67);

    auto g_xxxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 68);

    auto g_xxxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 69);

    auto g_xxxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 70);

    auto g_xxxxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 71);

    auto g_xxxxz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 72);

    auto g_xxxxz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 73);

    auto g_xxxxz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 74);

    auto g_xxxxz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 75);

    auto g_xxxxz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 76);

    auto g_xxxxz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 78);

    auto g_xxxxz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 79);

    auto g_xxxxz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 80);

    auto g_xxxxz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 81);

    auto g_xxxxz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 82);

    auto g_xxxxz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 83);

    auto g_xxxyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 84);

    auto g_xxxyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 85);

    auto g_xxxyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 86);

    auto g_xxxyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 87);

    auto g_xxxyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 88);

    auto g_xxxyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 89);

    auto g_xxxyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 90);

    auto g_xxxyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 91);

    auto g_xxxyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 92);

    auto g_xxxyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 93);

    auto g_xxxyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 94);

    auto g_xxxyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 95);

    auto g_xxxyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 96);

    auto g_xxxyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 97);

    auto g_xxxyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 98);

    auto g_xxxyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 99);

    auto g_xxxyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 100);

    auto g_xxxyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 101);

    auto g_xxxyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 102);

    auto g_xxxyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 103);

    auto g_xxxyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 104);

    auto g_xxxyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 105);

    auto g_xxxyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 106);

    auto g_xxxyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 107);

    auto g_xxxyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 108);

    auto g_xxxyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 109);

    auto g_xxxyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 110);

    auto g_xxxyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 111);

    auto g_xxxzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 140);

    auto g_xxxzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 141);

    auto g_xxxzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 142);

    auto g_xxxzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 143);

    auto g_xxxzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 144);

    auto g_xxxzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 145);

    auto g_xxxzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 146);

    auto g_xxxzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 147);

    auto g_xxxzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 148);

    auto g_xxxzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 149);

    auto g_xxxzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 150);

    auto g_xxxzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 151);

    auto g_xxxzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 152);

    auto g_xxxzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 153);

    auto g_xxxzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 154);

    auto g_xxxzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 155);

    auto g_xxxzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 156);

    auto g_xxxzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 157);

    auto g_xxxzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 158);

    auto g_xxxzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 159);

    auto g_xxxzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 160);

    auto g_xxxzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 161);

    auto g_xxxzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 162);

    auto g_xxxzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 163);

    auto g_xxxzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 164);

    auto g_xxxzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 165);

    auto g_xxxzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 166);

    auto g_xxxzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 167);

    auto g_xxyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 168);

    auto g_xxyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 169);

    auto g_xxyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 170);

    auto g_xxyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 171);

    auto g_xxyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 172);

    auto g_xxyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 173);

    auto g_xxyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 174);

    auto g_xxyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 175);

    auto g_xxyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 176);

    auto g_xxyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 177);

    auto g_xxyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 178);

    auto g_xxyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 179);

    auto g_xxyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 180);

    auto g_xxyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 181);

    auto g_xxyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 182);

    auto g_xxyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 183);

    auto g_xxyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 184);

    auto g_xxyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 185);

    auto g_xxyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 186);

    auto g_xxyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 187);

    auto g_xxyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 188);

    auto g_xxyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 189);

    auto g_xxyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 190);

    auto g_xxyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 191);

    auto g_xxyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 192);

    auto g_xxyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 193);

    auto g_xxyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 194);

    auto g_xxyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 195);

    auto g_xxyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 197);

    auto g_xxyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 199);

    auto g_xxyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 202);

    auto g_xxyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 206);

    auto g_xxyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 211);

    auto g_xxyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 224);

    auto g_xxyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 226);

    auto g_xxyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 229);

    auto g_xxyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 233);

    auto g_xxyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 238);

    auto g_xxyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 244);

    auto g_xxzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 252);

    auto g_xxzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 253);

    auto g_xxzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 254);

    auto g_xxzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 255);

    auto g_xxzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 256);

    auto g_xxzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 257);

    auto g_xxzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 258);

    auto g_xxzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 259);

    auto g_xxzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 260);

    auto g_xxzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 261);

    auto g_xxzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 262);

    auto g_xxzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 263);

    auto g_xxzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 264);

    auto g_xxzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 265);

    auto g_xxzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 266);

    auto g_xxzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 267);

    auto g_xxzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 268);

    auto g_xxzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 269);

    auto g_xxzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 270);

    auto g_xxzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 271);

    auto g_xxzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 272);

    auto g_xxzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 273);

    auto g_xxzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 274);

    auto g_xxzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 275);

    auto g_xxzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 276);

    auto g_xxzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 277);

    auto g_xxzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 278);

    auto g_xxzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 279);

    auto g_xyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 280);

    auto g_xyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 281);

    auto g_xyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 283);

    auto g_xyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 284);

    auto g_xyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 286);

    auto g_xyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 287);

    auto g_xyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 288);

    auto g_xyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 290);

    auto g_xyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 291);

    auto g_xyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 292);

    auto g_xyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 293);

    auto g_xyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 295);

    auto g_xyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 296);

    auto g_xyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 297);

    auto g_xyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 298);

    auto g_xyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 299);

    auto g_xyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 301);

    auto g_xyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 302);

    auto g_xyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 303);

    auto g_xyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 304);

    auto g_xyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 305);

    auto g_xyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 306);

    auto g_xyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 307);

    auto g_xyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 340);

    auto g_xyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 343);

    auto g_xyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 344);

    auto g_xyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 347);

    auto g_xyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 348);

    auto g_xyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 349);

    auto g_xyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 352);

    auto g_xyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 353);

    auto g_xyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 354);

    auto g_xyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 355);

    auto g_xyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 357);

    auto g_xyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 358);

    auto g_xyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 359);

    auto g_xyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 360);

    auto g_xyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 361);

    auto g_xyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 362);

    auto g_xyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 363);

    auto g_xzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 392);

    auto g_xzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 394);

    auto g_xzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 396);

    auto g_xzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 397);

    auto g_xzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 399);

    auto g_xzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 400);

    auto g_xzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 401);

    auto g_xzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 403);

    auto g_xzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 404);

    auto g_xzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 405);

    auto g_xzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 406);

    auto g_xzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 408);

    auto g_xzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 409);

    auto g_xzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 410);

    auto g_xzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 411);

    auto g_xzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 412);

    auto g_xzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 413);

    auto g_xzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 414);

    auto g_xzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 415);

    auto g_xzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 416);

    auto g_xzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 417);

    auto g_xzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 418);

    auto g_xzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 419);

    auto g_yyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 420);

    auto g_yyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 421);

    auto g_yyyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 422);

    auto g_yyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 423);

    auto g_yyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 424);

    auto g_yyyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 425);

    auto g_yyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 426);

    auto g_yyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 427);

    auto g_yyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 428);

    auto g_yyyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 429);

    auto g_yyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 430);

    auto g_yyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 431);

    auto g_yyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 432);

    auto g_yyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 433);

    auto g_yyyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 434);

    auto g_yyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 435);

    auto g_yyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 436);

    auto g_yyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 437);

    auto g_yyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 438);

    auto g_yyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 439);

    auto g_yyyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 440);

    auto g_yyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 441);

    auto g_yyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 442);

    auto g_yyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 443);

    auto g_yyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 444);

    auto g_yyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 445);

    auto g_yyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 446);

    auto g_yyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 447);

    auto g_yyyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 449);

    auto g_yyyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 450);

    auto g_yyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 451);

    auto g_yyyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 452);

    auto g_yyyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 453);

    auto g_yyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 454);

    auto g_yyyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 455);

    auto g_yyyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 456);

    auto g_yyyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 457);

    auto g_yyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 458);

    auto g_yyyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 459);

    auto g_yyyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 460);

    auto g_yyyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 461);

    auto g_yyyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 462);

    auto g_yyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 463);

    auto g_yyyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 464);

    auto g_yyyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 465);

    auto g_yyyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 466);

    auto g_yyyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 467);

    auto g_yyyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 468);

    auto g_yyyyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 469);

    auto g_yyyyz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 470);

    auto g_yyyyz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 471);

    auto g_yyyyz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 472);

    auto g_yyyyz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 473);

    auto g_yyyyz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 474);

    auto g_yyyyz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 475);

    auto g_yyyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 476);

    auto g_yyyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 477);

    auto g_yyyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 478);

    auto g_yyyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 479);

    auto g_yyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 480);

    auto g_yyyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 481);

    auto g_yyyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 482);

    auto g_yyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 483);

    auto g_yyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 484);

    auto g_yyyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 485);

    auto g_yyyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 486);

    auto g_yyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 487);

    auto g_yyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 488);

    auto g_yyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 489);

    auto g_yyyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 490);

    auto g_yyyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 491);

    auto g_yyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 492);

    auto g_yyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 493);

    auto g_yyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 494);

    auto g_yyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 495);

    auto g_yyyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 496);

    auto g_yyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 497);

    auto g_yyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 498);

    auto g_yyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 499);

    auto g_yyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 500);

    auto g_yyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 501);

    auto g_yyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 502);

    auto g_yyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 503);

    auto g_yyzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 504);

    auto g_yyzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 505);

    auto g_yyzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 506);

    auto g_yyzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 507);

    auto g_yyzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 508);

    auto g_yyzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 509);

    auto g_yyzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 510);

    auto g_yyzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 511);

    auto g_yyzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 512);

    auto g_yyzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 513);

    auto g_yyzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 514);

    auto g_yyzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 515);

    auto g_yyzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 516);

    auto g_yyzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 517);

    auto g_yyzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 518);

    auto g_yyzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 519);

    auto g_yyzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 520);

    auto g_yyzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 521);

    auto g_yyzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 522);

    auto g_yyzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 523);

    auto g_yyzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 524);

    auto g_yyzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 525);

    auto g_yyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 526);

    auto g_yyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 527);

    auto g_yyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 528);

    auto g_yyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 529);

    auto g_yyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 530);

    auto g_yyzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 531);

    auto g_yzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 532);

    auto g_yzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 533);

    auto g_yzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 534);

    auto g_yzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 535);

    auto g_yzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 536);

    auto g_yzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 537);

    auto g_yzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 538);

    auto g_yzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 539);

    auto g_yzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 540);

    auto g_yzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 541);

    auto g_yzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 542);

    auto g_yzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 543);

    auto g_yzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 544);

    auto g_yzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 545);

    auto g_yzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 546);

    auto g_yzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 547);

    auto g_yzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 548);

    auto g_yzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 549);

    auto g_yzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 550);

    auto g_yzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 551);

    auto g_yzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 552);

    auto g_yzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 553);

    auto g_yzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 554);

    auto g_yzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 555);

    auto g_yzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 556);

    auto g_yzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 557);

    auto g_yzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 558);

    auto g_yzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 559);

    auto g_zzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 560);

    auto g_zzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 561);

    auto g_zzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 562);

    auto g_zzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 563);

    auto g_zzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 564);

    auto g_zzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 565);

    auto g_zzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 566);

    auto g_zzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 567);

    auto g_zzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 568);

    auto g_zzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 569);

    auto g_zzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 570);

    auto g_zzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 571);

    auto g_zzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 572);

    auto g_zzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 573);

    auto g_zzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 574);

    auto g_zzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 575);

    auto g_zzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 576);

    auto g_zzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 577);

    auto g_zzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 578);

    auto g_zzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 579);

    auto g_zzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 580);

    auto g_zzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 581);

    auto g_zzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 582);

    auto g_zzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 583);

    auto g_zzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 584);

    auto g_zzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 585);

    auto g_zzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 586);

    auto g_zzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 587);

    /// Set up 0-28 components of targeted buffer : ISI

    auto g_xxxxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi);

    auto g_xxxxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 1);

    auto g_xxxxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 2);

    auto g_xxxxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 3);

    auto g_xxxxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 4);

    auto g_xxxxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 5);

    auto g_xxxxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 6);

    auto g_xxxxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 7);

    auto g_xxxxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 8);

    auto g_xxxxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 9);

    auto g_xxxxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 10);

    auto g_xxxxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 11);

    auto g_xxxxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 12);

    auto g_xxxxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 13);

    auto g_xxxxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 14);

    auto g_xxxxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 15);

    auto g_xxxxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 16);

    auto g_xxxxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 17);

    auto g_xxxxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 18);

    auto g_xxxxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 19);

    auto g_xxxxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 20);

    auto g_xxxxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 21);

    auto g_xxxxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 22);

    auto g_xxxxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 23);

    auto g_xxxxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 24);

    auto g_xxxxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 25);

    auto g_xxxxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 26);

    auto g_xxxxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 27);

    #pragma omp simd aligned(g_xxxx_0_xxxxxx_0, g_xxxx_0_xxxxxx_1, g_xxxx_0_xxxxxy_0, g_xxxx_0_xxxxxy_1, g_xxxx_0_xxxxxz_0, g_xxxx_0_xxxxxz_1, g_xxxx_0_xxxxyy_0, g_xxxx_0_xxxxyy_1, g_xxxx_0_xxxxyz_0, g_xxxx_0_xxxxyz_1, g_xxxx_0_xxxxzz_0, g_xxxx_0_xxxxzz_1, g_xxxx_0_xxxyyy_0, g_xxxx_0_xxxyyy_1, g_xxxx_0_xxxyyz_0, g_xxxx_0_xxxyyz_1, g_xxxx_0_xxxyzz_0, g_xxxx_0_xxxyzz_1, g_xxxx_0_xxxzzz_0, g_xxxx_0_xxxzzz_1, g_xxxx_0_xxyyyy_0, g_xxxx_0_xxyyyy_1, g_xxxx_0_xxyyyz_0, g_xxxx_0_xxyyyz_1, g_xxxx_0_xxyyzz_0, g_xxxx_0_xxyyzz_1, g_xxxx_0_xxyzzz_0, g_xxxx_0_xxyzzz_1, g_xxxx_0_xxzzzz_0, g_xxxx_0_xxzzzz_1, g_xxxx_0_xyyyyy_0, g_xxxx_0_xyyyyy_1, g_xxxx_0_xyyyyz_0, g_xxxx_0_xyyyyz_1, g_xxxx_0_xyyyzz_0, g_xxxx_0_xyyyzz_1, g_xxxx_0_xyyzzz_0, g_xxxx_0_xyyzzz_1, g_xxxx_0_xyzzzz_0, g_xxxx_0_xyzzzz_1, g_xxxx_0_xzzzzz_0, g_xxxx_0_xzzzzz_1, g_xxxx_0_yyyyyy_0, g_xxxx_0_yyyyyy_1, g_xxxx_0_yyyyyz_0, g_xxxx_0_yyyyyz_1, g_xxxx_0_yyyyzz_0, g_xxxx_0_yyyyzz_1, g_xxxx_0_yyyzzz_0, g_xxxx_0_yyyzzz_1, g_xxxx_0_yyzzzz_0, g_xxxx_0_yyzzzz_1, g_xxxx_0_yzzzzz_0, g_xxxx_0_yzzzzz_1, g_xxxx_0_zzzzzz_0, g_xxxx_0_zzzzzz_1, g_xxxxx_0_xxxxx_1, g_xxxxx_0_xxxxxx_1, g_xxxxx_0_xxxxxy_1, g_xxxxx_0_xxxxxz_1, g_xxxxx_0_xxxxy_1, g_xxxxx_0_xxxxyy_1, g_xxxxx_0_xxxxyz_1, g_xxxxx_0_xxxxz_1, g_xxxxx_0_xxxxzz_1, g_xxxxx_0_xxxyy_1, g_xxxxx_0_xxxyyy_1, g_xxxxx_0_xxxyyz_1, g_xxxxx_0_xxxyz_1, g_xxxxx_0_xxxyzz_1, g_xxxxx_0_xxxzz_1, g_xxxxx_0_xxxzzz_1, g_xxxxx_0_xxyyy_1, g_xxxxx_0_xxyyyy_1, g_xxxxx_0_xxyyyz_1, g_xxxxx_0_xxyyz_1, g_xxxxx_0_xxyyzz_1, g_xxxxx_0_xxyzz_1, g_xxxxx_0_xxyzzz_1, g_xxxxx_0_xxzzz_1, g_xxxxx_0_xxzzzz_1, g_xxxxx_0_xyyyy_1, g_xxxxx_0_xyyyyy_1, g_xxxxx_0_xyyyyz_1, g_xxxxx_0_xyyyz_1, g_xxxxx_0_xyyyzz_1, g_xxxxx_0_xyyzz_1, g_xxxxx_0_xyyzzz_1, g_xxxxx_0_xyzzz_1, g_xxxxx_0_xyzzzz_1, g_xxxxx_0_xzzzz_1, g_xxxxx_0_xzzzzz_1, g_xxxxx_0_yyyyy_1, g_xxxxx_0_yyyyyy_1, g_xxxxx_0_yyyyyz_1, g_xxxxx_0_yyyyz_1, g_xxxxx_0_yyyyzz_1, g_xxxxx_0_yyyzz_1, g_xxxxx_0_yyyzzz_1, g_xxxxx_0_yyzzz_1, g_xxxxx_0_yyzzzz_1, g_xxxxx_0_yzzzz_1, g_xxxxx_0_yzzzzz_1, g_xxxxx_0_zzzzz_1, g_xxxxx_0_zzzzzz_1, g_xxxxxx_0_xxxxxx_0, g_xxxxxx_0_xxxxxy_0, g_xxxxxx_0_xxxxxz_0, g_xxxxxx_0_xxxxyy_0, g_xxxxxx_0_xxxxyz_0, g_xxxxxx_0_xxxxzz_0, g_xxxxxx_0_xxxyyy_0, g_xxxxxx_0_xxxyyz_0, g_xxxxxx_0_xxxyzz_0, g_xxxxxx_0_xxxzzz_0, g_xxxxxx_0_xxyyyy_0, g_xxxxxx_0_xxyyyz_0, g_xxxxxx_0_xxyyzz_0, g_xxxxxx_0_xxyzzz_0, g_xxxxxx_0_xxzzzz_0, g_xxxxxx_0_xyyyyy_0, g_xxxxxx_0_xyyyyz_0, g_xxxxxx_0_xyyyzz_0, g_xxxxxx_0_xyyzzz_0, g_xxxxxx_0_xyzzzz_0, g_xxxxxx_0_xzzzzz_0, g_xxxxxx_0_yyyyyy_0, g_xxxxxx_0_yyyyyz_0, g_xxxxxx_0_yyyyzz_0, g_xxxxxx_0_yyyzzz_0, g_xxxxxx_0_yyzzzz_0, g_xxxxxx_0_yzzzzz_0, g_xxxxxx_0_zzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxx_0_xxxxxx_0[i] = 5.0 * g_xxxx_0_xxxxxx_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxx_1[i] * fz_be_0 + 6.0 * g_xxxxx_0_xxxxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxx_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxy_0[i] = 5.0 * g_xxxx_0_xxxxxy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxxxx_0_xxxxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxz_0[i] = 5.0 * g_xxxx_0_xxxxxz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxxxx_0_xxxxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxyy_0[i] = 5.0 * g_xxxx_0_xxxxyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxyz_0[i] = 5.0 * g_xxxx_0_xxxxyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxzz_0[i] = 5.0 * g_xxxx_0_xxxxzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyyy_0[i] = 5.0 * g_xxxx_0_xxxyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyyz_0[i] = 5.0 * g_xxxx_0_xxxyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyzz_0[i] = 5.0 * g_xxxx_0_xxxyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxzzz_0[i] = 5.0 * g_xxxx_0_xxxzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyyy_0[i] = 5.0 * g_xxxx_0_xxyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyyz_0[i] = 5.0 * g_xxxx_0_xxyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyzz_0[i] = 5.0 * g_xxxx_0_xxyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyzzz_0[i] = 5.0 * g_xxxx_0_xxyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxzzzz_0[i] = 5.0 * g_xxxx_0_xxzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyyy_0[i] = 5.0 * g_xxxx_0_xyyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyyy_1[i] * fz_be_0 + g_xxxxx_0_yyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyyz_0[i] = 5.0 * g_xxxx_0_xyyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyyz_1[i] * fz_be_0 + g_xxxxx_0_yyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyzz_0[i] = 5.0 * g_xxxx_0_xyyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyzz_1[i] * fz_be_0 + g_xxxxx_0_yyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyzzz_0[i] = 5.0 * g_xxxx_0_xyyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyzzz_1[i] * fz_be_0 + g_xxxxx_0_yyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyzzzz_0[i] = 5.0 * g_xxxx_0_xyzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyzzzz_1[i] * fz_be_0 + g_xxxxx_0_yzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xzzzzz_0[i] = 5.0 * g_xxxx_0_xzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xzzzzz_1[i] * fz_be_0 + g_xxxxx_0_zzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyyy_0[i] = 5.0 * g_xxxx_0_yyyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyyy_1[i] * fz_be_0 + g_xxxxx_0_yyyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyyz_0[i] = 5.0 * g_xxxx_0_yyyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyyz_1[i] * fz_be_0 + g_xxxxx_0_yyyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyzz_0[i] = 5.0 * g_xxxx_0_yyyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyzz_1[i] * fz_be_0 + g_xxxxx_0_yyyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyzzz_0[i] = 5.0 * g_xxxx_0_yyyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyzzz_1[i] * fz_be_0 + g_xxxxx_0_yyyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyzzzz_0[i] = 5.0 * g_xxxx_0_yyzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyzzzz_1[i] * fz_be_0 + g_xxxxx_0_yyzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yzzzzz_0[i] = 5.0 * g_xxxx_0_yzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yzzzzz_1[i] * fz_be_0 + g_xxxxx_0_yzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_zzzzzz_0[i] = 5.0 * g_xxxx_0_zzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_zzzzzz_1[i] * fz_be_0 + g_xxxxx_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 28-56 components of targeted buffer : ISI

    auto g_xxxxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 28);

    auto g_xxxxxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 29);

    auto g_xxxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 30);

    auto g_xxxxxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 31);

    auto g_xxxxxy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 32);

    auto g_xxxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 33);

    auto g_xxxxxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 34);

    auto g_xxxxxy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 35);

    auto g_xxxxxy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 36);

    auto g_xxxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 37);

    auto g_xxxxxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 38);

    auto g_xxxxxy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 39);

    auto g_xxxxxy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 40);

    auto g_xxxxxy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 41);

    auto g_xxxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 42);

    auto g_xxxxxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 43);

    auto g_xxxxxy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 44);

    auto g_xxxxxy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 45);

    auto g_xxxxxy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 46);

    auto g_xxxxxy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 47);

    auto g_xxxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 48);

    auto g_xxxxxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 49);

    auto g_xxxxxy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 50);

    auto g_xxxxxy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 51);

    auto g_xxxxxy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 52);

    auto g_xxxxxy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 53);

    auto g_xxxxxy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 54);

    auto g_xxxxxy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 55);

    #pragma omp simd aligned(g_xxxxx_0_xxxxx_1, g_xxxxx_0_xxxxxx_1, g_xxxxx_0_xxxxxy_1, g_xxxxx_0_xxxxxz_1, g_xxxxx_0_xxxxy_1, g_xxxxx_0_xxxxyy_1, g_xxxxx_0_xxxxyz_1, g_xxxxx_0_xxxxz_1, g_xxxxx_0_xxxxzz_1, g_xxxxx_0_xxxyy_1, g_xxxxx_0_xxxyyy_1, g_xxxxx_0_xxxyyz_1, g_xxxxx_0_xxxyz_1, g_xxxxx_0_xxxyzz_1, g_xxxxx_0_xxxzz_1, g_xxxxx_0_xxxzzz_1, g_xxxxx_0_xxyyy_1, g_xxxxx_0_xxyyyy_1, g_xxxxx_0_xxyyyz_1, g_xxxxx_0_xxyyz_1, g_xxxxx_0_xxyyzz_1, g_xxxxx_0_xxyzz_1, g_xxxxx_0_xxyzzz_1, g_xxxxx_0_xxzzz_1, g_xxxxx_0_xxzzzz_1, g_xxxxx_0_xyyyy_1, g_xxxxx_0_xyyyyy_1, g_xxxxx_0_xyyyyz_1, g_xxxxx_0_xyyyz_1, g_xxxxx_0_xyyyzz_1, g_xxxxx_0_xyyzz_1, g_xxxxx_0_xyyzzz_1, g_xxxxx_0_xyzzz_1, g_xxxxx_0_xyzzzz_1, g_xxxxx_0_xzzzz_1, g_xxxxx_0_xzzzzz_1, g_xxxxx_0_yyyyy_1, g_xxxxx_0_yyyyyy_1, g_xxxxx_0_yyyyyz_1, g_xxxxx_0_yyyyz_1, g_xxxxx_0_yyyyzz_1, g_xxxxx_0_yyyzz_1, g_xxxxx_0_yyyzzz_1, g_xxxxx_0_yyzzz_1, g_xxxxx_0_yyzzzz_1, g_xxxxx_0_yzzzz_1, g_xxxxx_0_yzzzzz_1, g_xxxxx_0_zzzzz_1, g_xxxxx_0_zzzzzz_1, g_xxxxxy_0_xxxxxx_0, g_xxxxxy_0_xxxxxy_0, g_xxxxxy_0_xxxxxz_0, g_xxxxxy_0_xxxxyy_0, g_xxxxxy_0_xxxxyz_0, g_xxxxxy_0_xxxxzz_0, g_xxxxxy_0_xxxyyy_0, g_xxxxxy_0_xxxyyz_0, g_xxxxxy_0_xxxyzz_0, g_xxxxxy_0_xxxzzz_0, g_xxxxxy_0_xxyyyy_0, g_xxxxxy_0_xxyyyz_0, g_xxxxxy_0_xxyyzz_0, g_xxxxxy_0_xxyzzz_0, g_xxxxxy_0_xxzzzz_0, g_xxxxxy_0_xyyyyy_0, g_xxxxxy_0_xyyyyz_0, g_xxxxxy_0_xyyyzz_0, g_xxxxxy_0_xyyzzz_0, g_xxxxxy_0_xyzzzz_0, g_xxxxxy_0_xzzzzz_0, g_xxxxxy_0_yyyyyy_0, g_xxxxxy_0_yyyyyz_0, g_xxxxxy_0_yyyyzz_0, g_xxxxxy_0_yyyzzz_0, g_xxxxxy_0_yyzzzz_0, g_xxxxxy_0_yzzzzz_0, g_xxxxxy_0_zzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxy_0_xxxxxx_0[i] = g_xxxxx_0_xxxxxx_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxy_0[i] = g_xxxxx_0_xxxxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxz_0[i] = g_xxxxx_0_xxxxxz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxyy_0[i] = 2.0 * g_xxxxx_0_xxxxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxyz_0[i] = g_xxxxx_0_xxxxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxzz_0[i] = g_xxxxx_0_xxxxzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyyy_0[i] = 3.0 * g_xxxxx_0_xxxyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyyz_0[i] = 2.0 * g_xxxxx_0_xxxyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyzz_0[i] = g_xxxxx_0_xxxzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxzzz_0[i] = g_xxxxx_0_xxxzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyyy_0[i] = 4.0 * g_xxxxx_0_xxyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyyz_0[i] = 3.0 * g_xxxxx_0_xxyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyzz_0[i] = 2.0 * g_xxxxx_0_xxyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyzzz_0[i] = g_xxxxx_0_xxzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxzzzz_0[i] = g_xxxxx_0_xxzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyyy_0[i] = 5.0 * g_xxxxx_0_xyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyyz_0[i] = 4.0 * g_xxxxx_0_xyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyzz_0[i] = 3.0 * g_xxxxx_0_xyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyzzz_0[i] = 2.0 * g_xxxxx_0_xyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyzzzz_0[i] = g_xxxxx_0_xzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xzzzzz_0[i] = g_xxxxx_0_xzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyyy_0[i] = 6.0 * g_xxxxx_0_yyyyy_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyyz_0[i] = 5.0 * g_xxxxx_0_yyyyz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyzz_0[i] = 4.0 * g_xxxxx_0_yyyzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyzzz_0[i] = 3.0 * g_xxxxx_0_yyzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyzzzz_0[i] = 2.0 * g_xxxxx_0_yzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yzzzzz_0[i] = g_xxxxx_0_zzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_zzzzzz_0[i] = g_xxxxx_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 56-84 components of targeted buffer : ISI

    auto g_xxxxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 56);

    auto g_xxxxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 57);

    auto g_xxxxxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 58);

    auto g_xxxxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 59);

    auto g_xxxxxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 60);

    auto g_xxxxxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 61);

    auto g_xxxxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 62);

    auto g_xxxxxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 63);

    auto g_xxxxxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 64);

    auto g_xxxxxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 65);

    auto g_xxxxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 66);

    auto g_xxxxxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 67);

    auto g_xxxxxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 68);

    auto g_xxxxxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 69);

    auto g_xxxxxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 70);

    auto g_xxxxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 71);

    auto g_xxxxxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 72);

    auto g_xxxxxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 73);

    auto g_xxxxxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 74);

    auto g_xxxxxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 75);

    auto g_xxxxxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 76);

    auto g_xxxxxz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 77);

    auto g_xxxxxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 78);

    auto g_xxxxxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 79);

    auto g_xxxxxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 80);

    auto g_xxxxxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 81);

    auto g_xxxxxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 82);

    auto g_xxxxxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 83);

    #pragma omp simd aligned(g_xxxxx_0_xxxxx_1, g_xxxxx_0_xxxxxx_1, g_xxxxx_0_xxxxxy_1, g_xxxxx_0_xxxxxz_1, g_xxxxx_0_xxxxy_1, g_xxxxx_0_xxxxyy_1, g_xxxxx_0_xxxxyz_1, g_xxxxx_0_xxxxz_1, g_xxxxx_0_xxxxzz_1, g_xxxxx_0_xxxyy_1, g_xxxxx_0_xxxyyy_1, g_xxxxx_0_xxxyyz_1, g_xxxxx_0_xxxyz_1, g_xxxxx_0_xxxyzz_1, g_xxxxx_0_xxxzz_1, g_xxxxx_0_xxxzzz_1, g_xxxxx_0_xxyyy_1, g_xxxxx_0_xxyyyy_1, g_xxxxx_0_xxyyyz_1, g_xxxxx_0_xxyyz_1, g_xxxxx_0_xxyyzz_1, g_xxxxx_0_xxyzz_1, g_xxxxx_0_xxyzzz_1, g_xxxxx_0_xxzzz_1, g_xxxxx_0_xxzzzz_1, g_xxxxx_0_xyyyy_1, g_xxxxx_0_xyyyyy_1, g_xxxxx_0_xyyyyz_1, g_xxxxx_0_xyyyz_1, g_xxxxx_0_xyyyzz_1, g_xxxxx_0_xyyzz_1, g_xxxxx_0_xyyzzz_1, g_xxxxx_0_xyzzz_1, g_xxxxx_0_xyzzzz_1, g_xxxxx_0_xzzzz_1, g_xxxxx_0_xzzzzz_1, g_xxxxx_0_yyyyy_1, g_xxxxx_0_yyyyyy_1, g_xxxxx_0_yyyyyz_1, g_xxxxx_0_yyyyz_1, g_xxxxx_0_yyyyzz_1, g_xxxxx_0_yyyzz_1, g_xxxxx_0_yyyzzz_1, g_xxxxx_0_yyzzz_1, g_xxxxx_0_yyzzzz_1, g_xxxxx_0_yzzzz_1, g_xxxxx_0_yzzzzz_1, g_xxxxx_0_zzzzz_1, g_xxxxx_0_zzzzzz_1, g_xxxxxz_0_xxxxxx_0, g_xxxxxz_0_xxxxxy_0, g_xxxxxz_0_xxxxxz_0, g_xxxxxz_0_xxxxyy_0, g_xxxxxz_0_xxxxyz_0, g_xxxxxz_0_xxxxzz_0, g_xxxxxz_0_xxxyyy_0, g_xxxxxz_0_xxxyyz_0, g_xxxxxz_0_xxxyzz_0, g_xxxxxz_0_xxxzzz_0, g_xxxxxz_0_xxyyyy_0, g_xxxxxz_0_xxyyyz_0, g_xxxxxz_0_xxyyzz_0, g_xxxxxz_0_xxyzzz_0, g_xxxxxz_0_xxzzzz_0, g_xxxxxz_0_xyyyyy_0, g_xxxxxz_0_xyyyyz_0, g_xxxxxz_0_xyyyzz_0, g_xxxxxz_0_xyyzzz_0, g_xxxxxz_0_xyzzzz_0, g_xxxxxz_0_xzzzzz_0, g_xxxxxz_0_yyyyyy_0, g_xxxxxz_0_yyyyyz_0, g_xxxxxz_0_yyyyzz_0, g_xxxxxz_0_yyyzzz_0, g_xxxxxz_0_yyzzzz_0, g_xxxxxz_0_yzzzzz_0, g_xxxxxz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxz_0_xxxxxx_0[i] = g_xxxxx_0_xxxxxx_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxy_0[i] = g_xxxxx_0_xxxxxy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxz_0[i] = g_xxxxx_0_xxxxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxyy_0[i] = g_xxxxx_0_xxxxyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxyz_0[i] = g_xxxxx_0_xxxxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxzz_0[i] = 2.0 * g_xxxxx_0_xxxxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyyy_0[i] = g_xxxxx_0_xxxyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyyz_0[i] = g_xxxxx_0_xxxyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyzz_0[i] = 2.0 * g_xxxxx_0_xxxyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxzzz_0[i] = 3.0 * g_xxxxx_0_xxxzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyyy_0[i] = g_xxxxx_0_xxyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyyz_0[i] = g_xxxxx_0_xxyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyzz_0[i] = 2.0 * g_xxxxx_0_xxyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyzzz_0[i] = 3.0 * g_xxxxx_0_xxyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxzzzz_0[i] = 4.0 * g_xxxxx_0_xxzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyyy_0[i] = g_xxxxx_0_xyyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyyz_0[i] = g_xxxxx_0_xyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyzz_0[i] = 2.0 * g_xxxxx_0_xyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyzzz_0[i] = 3.0 * g_xxxxx_0_xyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyzzzz_0[i] = 4.0 * g_xxxxx_0_xyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xzzzzz_0[i] = 5.0 * g_xxxxx_0_xzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyyy_0[i] = g_xxxxx_0_yyyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyyz_0[i] = g_xxxxx_0_yyyyy_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyzz_0[i] = 2.0 * g_xxxxx_0_yyyyz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyzzz_0[i] = 3.0 * g_xxxxx_0_yyyzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyzzzz_0[i] = 4.0 * g_xxxxx_0_yyzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yzzzzz_0[i] = 5.0 * g_xxxxx_0_yzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_zzzzzz_0[i] = 6.0 * g_xxxxx_0_zzzzz_1[i] * fi_acd_0 + g_xxxxx_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 84-112 components of targeted buffer : ISI

    auto g_xxxxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 84);

    auto g_xxxxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 85);

    auto g_xxxxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 86);

    auto g_xxxxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 87);

    auto g_xxxxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 88);

    auto g_xxxxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 89);

    auto g_xxxxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 90);

    auto g_xxxxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 91);

    auto g_xxxxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 92);

    auto g_xxxxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 93);

    auto g_xxxxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 94);

    auto g_xxxxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 95);

    auto g_xxxxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 96);

    auto g_xxxxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 97);

    auto g_xxxxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 98);

    auto g_xxxxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 99);

    auto g_xxxxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 100);

    auto g_xxxxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 101);

    auto g_xxxxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 102);

    auto g_xxxxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 103);

    auto g_xxxxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 104);

    auto g_xxxxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 105);

    auto g_xxxxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 106);

    auto g_xxxxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 107);

    auto g_xxxxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 108);

    auto g_xxxxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 109);

    auto g_xxxxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 110);

    auto g_xxxxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 111);

    #pragma omp simd aligned(g_xxxx_0_xxxxxx_0, g_xxxx_0_xxxxxx_1, g_xxxx_0_xxxxxz_0, g_xxxx_0_xxxxxz_1, g_xxxx_0_xxxxzz_0, g_xxxx_0_xxxxzz_1, g_xxxx_0_xxxzzz_0, g_xxxx_0_xxxzzz_1, g_xxxx_0_xxzzzz_0, g_xxxx_0_xxzzzz_1, g_xxxx_0_xzzzzz_0, g_xxxx_0_xzzzzz_1, g_xxxxy_0_xxxxxx_1, g_xxxxy_0_xxxxxz_1, g_xxxxy_0_xxxxzz_1, g_xxxxy_0_xxxzzz_1, g_xxxxy_0_xxzzzz_1, g_xxxxy_0_xzzzzz_1, g_xxxxyy_0_xxxxxx_0, g_xxxxyy_0_xxxxxy_0, g_xxxxyy_0_xxxxxz_0, g_xxxxyy_0_xxxxyy_0, g_xxxxyy_0_xxxxyz_0, g_xxxxyy_0_xxxxzz_0, g_xxxxyy_0_xxxyyy_0, g_xxxxyy_0_xxxyyz_0, g_xxxxyy_0_xxxyzz_0, g_xxxxyy_0_xxxzzz_0, g_xxxxyy_0_xxyyyy_0, g_xxxxyy_0_xxyyyz_0, g_xxxxyy_0_xxyyzz_0, g_xxxxyy_0_xxyzzz_0, g_xxxxyy_0_xxzzzz_0, g_xxxxyy_0_xyyyyy_0, g_xxxxyy_0_xyyyyz_0, g_xxxxyy_0_xyyyzz_0, g_xxxxyy_0_xyyzzz_0, g_xxxxyy_0_xyzzzz_0, g_xxxxyy_0_xzzzzz_0, g_xxxxyy_0_yyyyyy_0, g_xxxxyy_0_yyyyyz_0, g_xxxxyy_0_yyyyzz_0, g_xxxxyy_0_yyyzzz_0, g_xxxxyy_0_yyzzzz_0, g_xxxxyy_0_yzzzzz_0, g_xxxxyy_0_zzzzzz_0, g_xxxyy_0_xxxxxy_1, g_xxxyy_0_xxxxy_1, g_xxxyy_0_xxxxyy_1, g_xxxyy_0_xxxxyz_1, g_xxxyy_0_xxxyy_1, g_xxxyy_0_xxxyyy_1, g_xxxyy_0_xxxyyz_1, g_xxxyy_0_xxxyz_1, g_xxxyy_0_xxxyzz_1, g_xxxyy_0_xxyyy_1, g_xxxyy_0_xxyyyy_1, g_xxxyy_0_xxyyyz_1, g_xxxyy_0_xxyyz_1, g_xxxyy_0_xxyyzz_1, g_xxxyy_0_xxyzz_1, g_xxxyy_0_xxyzzz_1, g_xxxyy_0_xyyyy_1, g_xxxyy_0_xyyyyy_1, g_xxxyy_0_xyyyyz_1, g_xxxyy_0_xyyyz_1, g_xxxyy_0_xyyyzz_1, g_xxxyy_0_xyyzz_1, g_xxxyy_0_xyyzzz_1, g_xxxyy_0_xyzzz_1, g_xxxyy_0_xyzzzz_1, g_xxxyy_0_yyyyy_1, g_xxxyy_0_yyyyyy_1, g_xxxyy_0_yyyyyz_1, g_xxxyy_0_yyyyz_1, g_xxxyy_0_yyyyzz_1, g_xxxyy_0_yyyzz_1, g_xxxyy_0_yyyzzz_1, g_xxxyy_0_yyzzz_1, g_xxxyy_0_yyzzzz_1, g_xxxyy_0_yzzzz_1, g_xxxyy_0_yzzzzz_1, g_xxxyy_0_zzzzzz_1, g_xxyy_0_xxxxxy_0, g_xxyy_0_xxxxxy_1, g_xxyy_0_xxxxyy_0, g_xxyy_0_xxxxyy_1, g_xxyy_0_xxxxyz_0, g_xxyy_0_xxxxyz_1, g_xxyy_0_xxxyyy_0, g_xxyy_0_xxxyyy_1, g_xxyy_0_xxxyyz_0, g_xxyy_0_xxxyyz_1, g_xxyy_0_xxxyzz_0, g_xxyy_0_xxxyzz_1, g_xxyy_0_xxyyyy_0, g_xxyy_0_xxyyyy_1, g_xxyy_0_xxyyyz_0, g_xxyy_0_xxyyyz_1, g_xxyy_0_xxyyzz_0, g_xxyy_0_xxyyzz_1, g_xxyy_0_xxyzzz_0, g_xxyy_0_xxyzzz_1, g_xxyy_0_xyyyyy_0, g_xxyy_0_xyyyyy_1, g_xxyy_0_xyyyyz_0, g_xxyy_0_xyyyyz_1, g_xxyy_0_xyyyzz_0, g_xxyy_0_xyyyzz_1, g_xxyy_0_xyyzzz_0, g_xxyy_0_xyyzzz_1, g_xxyy_0_xyzzzz_0, g_xxyy_0_xyzzzz_1, g_xxyy_0_yyyyyy_0, g_xxyy_0_yyyyyy_1, g_xxyy_0_yyyyyz_0, g_xxyy_0_yyyyyz_1, g_xxyy_0_yyyyzz_0, g_xxyy_0_yyyyzz_1, g_xxyy_0_yyyzzz_0, g_xxyy_0_yyyzzz_1, g_xxyy_0_yyzzzz_0, g_xxyy_0_yyzzzz_1, g_xxyy_0_yzzzzz_0, g_xxyy_0_yzzzzz_1, g_xxyy_0_zzzzzz_0, g_xxyy_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyy_0_xxxxxx_0[i] = g_xxxx_0_xxxxxx_0[i] * fbe_0 - g_xxxx_0_xxxxxx_1[i] * fz_be_0 + g_xxxxy_0_xxxxxx_1[i] * wa_y[i];

        g_xxxxyy_0_xxxxxy_0[i] = 3.0 * g_xxyy_0_xxxxxy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxxyy_0_xxxxy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxxz_0[i] = g_xxxx_0_xxxxxz_0[i] * fbe_0 - g_xxxx_0_xxxxxz_1[i] * fz_be_0 + g_xxxxy_0_xxxxxz_1[i] * wa_y[i];

        g_xxxxyy_0_xxxxyy_0[i] = 3.0 * g_xxyy_0_xxxxyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxxyy_0_xxxyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxyz_0[i] = 3.0 * g_xxyy_0_xxxxyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxyy_0_xxxyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxzz_0[i] = g_xxxx_0_xxxxzz_0[i] * fbe_0 - g_xxxx_0_xxxxzz_1[i] * fz_be_0 + g_xxxxy_0_xxxxzz_1[i] * wa_y[i];

        g_xxxxyy_0_xxxyyy_0[i] = 3.0 * g_xxyy_0_xxxyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxyyz_0[i] = 3.0 * g_xxyy_0_xxxyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxyzz_0[i] = 3.0 * g_xxyy_0_xxxyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxzzz_0[i] = g_xxxx_0_xxxzzz_0[i] * fbe_0 - g_xxxx_0_xxxzzz_1[i] * fz_be_0 + g_xxxxy_0_xxxzzz_1[i] * wa_y[i];

        g_xxxxyy_0_xxyyyy_0[i] = 3.0 * g_xxyy_0_xxyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxyyyz_0[i] = 3.0 * g_xxyy_0_xxyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxyyzz_0[i] = 3.0 * g_xxyy_0_xxyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxyzzz_0[i] = 3.0 * g_xxyy_0_xxyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxzzzz_0[i] = g_xxxx_0_xxzzzz_0[i] * fbe_0 - g_xxxx_0_xxzzzz_1[i] * fz_be_0 + g_xxxxy_0_xxzzzz_1[i] * wa_y[i];

        g_xxxxyy_0_xyyyyy_0[i] = 3.0 * g_xxyy_0_xyyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyyy_1[i] * fz_be_0 + g_xxxyy_0_yyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xyyyyz_0[i] = 3.0 * g_xxyy_0_xyyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyyz_1[i] * fz_be_0 + g_xxxyy_0_yyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xyyyzz_0[i] = 3.0 * g_xxyy_0_xyyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyzz_1[i] * fz_be_0 + g_xxxyy_0_yyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xyyzzz_0[i] = 3.0 * g_xxyy_0_xyyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyzzz_1[i] * fz_be_0 + g_xxxyy_0_yyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xyzzzz_0[i] = 3.0 * g_xxyy_0_xyzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyzzzz_1[i] * fz_be_0 + g_xxxyy_0_yzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xzzzzz_0[i] = g_xxxx_0_xzzzzz_0[i] * fbe_0 - g_xxxx_0_xzzzzz_1[i] * fz_be_0 + g_xxxxy_0_xzzzzz_1[i] * wa_y[i];

        g_xxxxyy_0_yyyyyy_0[i] = 3.0 * g_xxyy_0_yyyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyyy_1[i] * fz_be_0 + g_xxxyy_0_yyyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_yyyyyz_0[i] = 3.0 * g_xxyy_0_yyyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyyz_1[i] * fz_be_0 + g_xxxyy_0_yyyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_yyyyzz_0[i] = 3.0 * g_xxyy_0_yyyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyzz_1[i] * fz_be_0 + g_xxxyy_0_yyyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_yyyzzz_0[i] = 3.0 * g_xxyy_0_yyyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyzzz_1[i] * fz_be_0 + g_xxxyy_0_yyyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_yyzzzz_0[i] = 3.0 * g_xxyy_0_yyzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyzzzz_1[i] * fz_be_0 + g_xxxyy_0_yyzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_yzzzzz_0[i] = 3.0 * g_xxyy_0_yzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yzzzzz_1[i] * fz_be_0 + g_xxxyy_0_yzzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_zzzzzz_0[i] = 3.0 * g_xxyy_0_zzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_zzzzzz_1[i] * fz_be_0 + g_xxxyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 112-140 components of targeted buffer : ISI

    auto g_xxxxyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 112);

    auto g_xxxxyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 113);

    auto g_xxxxyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 114);

    auto g_xxxxyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 115);

    auto g_xxxxyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 116);

    auto g_xxxxyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 117);

    auto g_xxxxyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 118);

    auto g_xxxxyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 119);

    auto g_xxxxyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 120);

    auto g_xxxxyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 121);

    auto g_xxxxyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 122);

    auto g_xxxxyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 123);

    auto g_xxxxyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 124);

    auto g_xxxxyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 125);

    auto g_xxxxyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 126);

    auto g_xxxxyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 127);

    auto g_xxxxyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 128);

    auto g_xxxxyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 129);

    auto g_xxxxyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 130);

    auto g_xxxxyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 131);

    auto g_xxxxyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 132);

    auto g_xxxxyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 133);

    auto g_xxxxyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 134);

    auto g_xxxxyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 135);

    auto g_xxxxyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 136);

    auto g_xxxxyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 137);

    auto g_xxxxyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 138);

    auto g_xxxxyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 139);

    #pragma omp simd aligned(g_xxxxy_0_xxxxxy_1, g_xxxxy_0_xxxxyy_1, g_xxxxy_0_xxxyyy_1, g_xxxxy_0_xxyyyy_1, g_xxxxy_0_xyyyyy_1, g_xxxxy_0_yyyyyy_1, g_xxxxyz_0_xxxxxx_0, g_xxxxyz_0_xxxxxy_0, g_xxxxyz_0_xxxxxz_0, g_xxxxyz_0_xxxxyy_0, g_xxxxyz_0_xxxxyz_0, g_xxxxyz_0_xxxxzz_0, g_xxxxyz_0_xxxyyy_0, g_xxxxyz_0_xxxyyz_0, g_xxxxyz_0_xxxyzz_0, g_xxxxyz_0_xxxzzz_0, g_xxxxyz_0_xxyyyy_0, g_xxxxyz_0_xxyyyz_0, g_xxxxyz_0_xxyyzz_0, g_xxxxyz_0_xxyzzz_0, g_xxxxyz_0_xxzzzz_0, g_xxxxyz_0_xyyyyy_0, g_xxxxyz_0_xyyyyz_0, g_xxxxyz_0_xyyyzz_0, g_xxxxyz_0_xyyzzz_0, g_xxxxyz_0_xyzzzz_0, g_xxxxyz_0_xzzzzz_0, g_xxxxyz_0_yyyyyy_0, g_xxxxyz_0_yyyyyz_0, g_xxxxyz_0_yyyyzz_0, g_xxxxyz_0_yyyzzz_0, g_xxxxyz_0_yyzzzz_0, g_xxxxyz_0_yzzzzz_0, g_xxxxyz_0_zzzzzz_0, g_xxxxz_0_xxxxxx_1, g_xxxxz_0_xxxxxz_1, g_xxxxz_0_xxxxyz_1, g_xxxxz_0_xxxxz_1, g_xxxxz_0_xxxxzz_1, g_xxxxz_0_xxxyyz_1, g_xxxxz_0_xxxyz_1, g_xxxxz_0_xxxyzz_1, g_xxxxz_0_xxxzz_1, g_xxxxz_0_xxxzzz_1, g_xxxxz_0_xxyyyz_1, g_xxxxz_0_xxyyz_1, g_xxxxz_0_xxyyzz_1, g_xxxxz_0_xxyzz_1, g_xxxxz_0_xxyzzz_1, g_xxxxz_0_xxzzz_1, g_xxxxz_0_xxzzzz_1, g_xxxxz_0_xyyyyz_1, g_xxxxz_0_xyyyz_1, g_xxxxz_0_xyyyzz_1, g_xxxxz_0_xyyzz_1, g_xxxxz_0_xyyzzz_1, g_xxxxz_0_xyzzz_1, g_xxxxz_0_xyzzzz_1, g_xxxxz_0_xzzzz_1, g_xxxxz_0_xzzzzz_1, g_xxxxz_0_yyyyyz_1, g_xxxxz_0_yyyyz_1, g_xxxxz_0_yyyyzz_1, g_xxxxz_0_yyyzz_1, g_xxxxz_0_yyyzzz_1, g_xxxxz_0_yyzzz_1, g_xxxxz_0_yyzzzz_1, g_xxxxz_0_yzzzz_1, g_xxxxz_0_yzzzzz_1, g_xxxxz_0_zzzzz_1, g_xxxxz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyz_0_xxxxxx_0[i] = g_xxxxz_0_xxxxxx_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxxy_0[i] = g_xxxxy_0_xxxxxy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxxxz_0[i] = g_xxxxz_0_xxxxxz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxyy_0[i] = g_xxxxy_0_xxxxyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxxyz_0[i] = g_xxxxz_0_xxxxz_1[i] * fi_acd_0 + g_xxxxz_0_xxxxyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxzz_0[i] = g_xxxxz_0_xxxxzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxyyy_0[i] = g_xxxxy_0_xxxyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxyyz_0[i] = 2.0 * g_xxxxz_0_xxxyz_1[i] * fi_acd_0 + g_xxxxz_0_xxxyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxyzz_0[i] = g_xxxxz_0_xxxzz_1[i] * fi_acd_0 + g_xxxxz_0_xxxyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxzzz_0[i] = g_xxxxz_0_xxxzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyyyy_0[i] = g_xxxxy_0_xxyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxyyyz_0[i] = 3.0 * g_xxxxz_0_xxyyz_1[i] * fi_acd_0 + g_xxxxz_0_xxyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyyzz_0[i] = 2.0 * g_xxxxz_0_xxyzz_1[i] * fi_acd_0 + g_xxxxz_0_xxyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyzzz_0[i] = g_xxxxz_0_xxzzz_1[i] * fi_acd_0 + g_xxxxz_0_xxyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxzzzz_0[i] = g_xxxxz_0_xxzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyyyy_0[i] = g_xxxxy_0_xyyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xyyyyz_0[i] = 4.0 * g_xxxxz_0_xyyyz_1[i] * fi_acd_0 + g_xxxxz_0_xyyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyyzz_0[i] = 3.0 * g_xxxxz_0_xyyzz_1[i] * fi_acd_0 + g_xxxxz_0_xyyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyzzz_0[i] = 2.0 * g_xxxxz_0_xyzzz_1[i] * fi_acd_0 + g_xxxxz_0_xyyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyzzzz_0[i] = g_xxxxz_0_xzzzz_1[i] * fi_acd_0 + g_xxxxz_0_xyzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xzzzzz_0[i] = g_xxxxz_0_xzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyyyy_0[i] = g_xxxxy_0_yyyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_yyyyyz_0[i] = 5.0 * g_xxxxz_0_yyyyz_1[i] * fi_acd_0 + g_xxxxz_0_yyyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyyzz_0[i] = 4.0 * g_xxxxz_0_yyyzz_1[i] * fi_acd_0 + g_xxxxz_0_yyyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyzzz_0[i] = 3.0 * g_xxxxz_0_yyzzz_1[i] * fi_acd_0 + g_xxxxz_0_yyyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyzzzz_0[i] = 2.0 * g_xxxxz_0_yzzzz_1[i] * fi_acd_0 + g_xxxxz_0_yyzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yzzzzz_0[i] = g_xxxxz_0_zzzzz_1[i] * fi_acd_0 + g_xxxxz_0_yzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_zzzzzz_0[i] = g_xxxxz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 140-168 components of targeted buffer : ISI

    auto g_xxxxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 140);

    auto g_xxxxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 141);

    auto g_xxxxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 142);

    auto g_xxxxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 143);

    auto g_xxxxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 144);

    auto g_xxxxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 145);

    auto g_xxxxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 146);

    auto g_xxxxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 147);

    auto g_xxxxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 148);

    auto g_xxxxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 149);

    auto g_xxxxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 150);

    auto g_xxxxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 151);

    auto g_xxxxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 152);

    auto g_xxxxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 153);

    auto g_xxxxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 154);

    auto g_xxxxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 155);

    auto g_xxxxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 156);

    auto g_xxxxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 157);

    auto g_xxxxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 158);

    auto g_xxxxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 159);

    auto g_xxxxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 160);

    auto g_xxxxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 161);

    auto g_xxxxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 162);

    auto g_xxxxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 163);

    auto g_xxxxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 164);

    auto g_xxxxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 165);

    auto g_xxxxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 166);

    auto g_xxxxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 167);

    #pragma omp simd aligned(g_xxxx_0_xxxxxx_0, g_xxxx_0_xxxxxx_1, g_xxxx_0_xxxxxy_0, g_xxxx_0_xxxxxy_1, g_xxxx_0_xxxxyy_0, g_xxxx_0_xxxxyy_1, g_xxxx_0_xxxyyy_0, g_xxxx_0_xxxyyy_1, g_xxxx_0_xxyyyy_0, g_xxxx_0_xxyyyy_1, g_xxxx_0_xyyyyy_0, g_xxxx_0_xyyyyy_1, g_xxxxz_0_xxxxxx_1, g_xxxxz_0_xxxxxy_1, g_xxxxz_0_xxxxyy_1, g_xxxxz_0_xxxyyy_1, g_xxxxz_0_xxyyyy_1, g_xxxxz_0_xyyyyy_1, g_xxxxzz_0_xxxxxx_0, g_xxxxzz_0_xxxxxy_0, g_xxxxzz_0_xxxxxz_0, g_xxxxzz_0_xxxxyy_0, g_xxxxzz_0_xxxxyz_0, g_xxxxzz_0_xxxxzz_0, g_xxxxzz_0_xxxyyy_0, g_xxxxzz_0_xxxyyz_0, g_xxxxzz_0_xxxyzz_0, g_xxxxzz_0_xxxzzz_0, g_xxxxzz_0_xxyyyy_0, g_xxxxzz_0_xxyyyz_0, g_xxxxzz_0_xxyyzz_0, g_xxxxzz_0_xxyzzz_0, g_xxxxzz_0_xxzzzz_0, g_xxxxzz_0_xyyyyy_0, g_xxxxzz_0_xyyyyz_0, g_xxxxzz_0_xyyyzz_0, g_xxxxzz_0_xyyzzz_0, g_xxxxzz_0_xyzzzz_0, g_xxxxzz_0_xzzzzz_0, g_xxxxzz_0_yyyyyy_0, g_xxxxzz_0_yyyyyz_0, g_xxxxzz_0_yyyyzz_0, g_xxxxzz_0_yyyzzz_0, g_xxxxzz_0_yyzzzz_0, g_xxxxzz_0_yzzzzz_0, g_xxxxzz_0_zzzzzz_0, g_xxxzz_0_xxxxxz_1, g_xxxzz_0_xxxxyz_1, g_xxxzz_0_xxxxz_1, g_xxxzz_0_xxxxzz_1, g_xxxzz_0_xxxyyz_1, g_xxxzz_0_xxxyz_1, g_xxxzz_0_xxxyzz_1, g_xxxzz_0_xxxzz_1, g_xxxzz_0_xxxzzz_1, g_xxxzz_0_xxyyyz_1, g_xxxzz_0_xxyyz_1, g_xxxzz_0_xxyyzz_1, g_xxxzz_0_xxyzz_1, g_xxxzz_0_xxyzzz_1, g_xxxzz_0_xxzzz_1, g_xxxzz_0_xxzzzz_1, g_xxxzz_0_xyyyyz_1, g_xxxzz_0_xyyyz_1, g_xxxzz_0_xyyyzz_1, g_xxxzz_0_xyyzz_1, g_xxxzz_0_xyyzzz_1, g_xxxzz_0_xyzzz_1, g_xxxzz_0_xyzzzz_1, g_xxxzz_0_xzzzz_1, g_xxxzz_0_xzzzzz_1, g_xxxzz_0_yyyyyy_1, g_xxxzz_0_yyyyyz_1, g_xxxzz_0_yyyyz_1, g_xxxzz_0_yyyyzz_1, g_xxxzz_0_yyyzz_1, g_xxxzz_0_yyyzzz_1, g_xxxzz_0_yyzzz_1, g_xxxzz_0_yyzzzz_1, g_xxxzz_0_yzzzz_1, g_xxxzz_0_yzzzzz_1, g_xxxzz_0_zzzzz_1, g_xxxzz_0_zzzzzz_1, g_xxzz_0_xxxxxz_0, g_xxzz_0_xxxxxz_1, g_xxzz_0_xxxxyz_0, g_xxzz_0_xxxxyz_1, g_xxzz_0_xxxxzz_0, g_xxzz_0_xxxxzz_1, g_xxzz_0_xxxyyz_0, g_xxzz_0_xxxyyz_1, g_xxzz_0_xxxyzz_0, g_xxzz_0_xxxyzz_1, g_xxzz_0_xxxzzz_0, g_xxzz_0_xxxzzz_1, g_xxzz_0_xxyyyz_0, g_xxzz_0_xxyyyz_1, g_xxzz_0_xxyyzz_0, g_xxzz_0_xxyyzz_1, g_xxzz_0_xxyzzz_0, g_xxzz_0_xxyzzz_1, g_xxzz_0_xxzzzz_0, g_xxzz_0_xxzzzz_1, g_xxzz_0_xyyyyz_0, g_xxzz_0_xyyyyz_1, g_xxzz_0_xyyyzz_0, g_xxzz_0_xyyyzz_1, g_xxzz_0_xyyzzz_0, g_xxzz_0_xyyzzz_1, g_xxzz_0_xyzzzz_0, g_xxzz_0_xyzzzz_1, g_xxzz_0_xzzzzz_0, g_xxzz_0_xzzzzz_1, g_xxzz_0_yyyyyy_0, g_xxzz_0_yyyyyy_1, g_xxzz_0_yyyyyz_0, g_xxzz_0_yyyyyz_1, g_xxzz_0_yyyyzz_0, g_xxzz_0_yyyyzz_1, g_xxzz_0_yyyzzz_0, g_xxzz_0_yyyzzz_1, g_xxzz_0_yyzzzz_0, g_xxzz_0_yyzzzz_1, g_xxzz_0_yzzzzz_0, g_xxzz_0_yzzzzz_1, g_xxzz_0_zzzzzz_0, g_xxzz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzz_0_xxxxxx_0[i] = g_xxxx_0_xxxxxx_0[i] * fbe_0 - g_xxxx_0_xxxxxx_1[i] * fz_be_0 + g_xxxxz_0_xxxxxx_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxxy_0[i] = g_xxxx_0_xxxxxy_0[i] * fbe_0 - g_xxxx_0_xxxxxy_1[i] * fz_be_0 + g_xxxxz_0_xxxxxy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxxz_0[i] = 3.0 * g_xxzz_0_xxxxxz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxxzz_0_xxxxz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxyy_0[i] = g_xxxx_0_xxxxyy_0[i] * fbe_0 - g_xxxx_0_xxxxyy_1[i] * fz_be_0 + g_xxxxz_0_xxxxyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxyz_0[i] = 3.0 * g_xxzz_0_xxxxyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxzz_0_xxxyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxzz_0[i] = 3.0 * g_xxzz_0_xxxxzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxxzz_0_xxxzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxyyy_0[i] = g_xxxx_0_xxxyyy_0[i] * fbe_0 - g_xxxx_0_xxxyyy_1[i] * fz_be_0 + g_xxxxz_0_xxxyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxyyz_0[i] = 3.0 * g_xxzz_0_xxxyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxyzz_0[i] = 3.0 * g_xxzz_0_xxxyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxzzz_0[i] = 3.0 * g_xxzz_0_xxxzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyyyy_0[i] = g_xxxx_0_xxyyyy_0[i] * fbe_0 - g_xxxx_0_xxyyyy_1[i] * fz_be_0 + g_xxxxz_0_xxyyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxyyyz_0[i] = 3.0 * g_xxzz_0_xxyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyyzz_0[i] = 3.0 * g_xxzz_0_xxyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyzzz_0[i] = 3.0 * g_xxzz_0_xxyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxzzzz_0[i] = 3.0 * g_xxzz_0_xxzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyyyy_0[i] = g_xxxx_0_xyyyyy_0[i] * fbe_0 - g_xxxx_0_xyyyyy_1[i] * fz_be_0 + g_xxxxz_0_xyyyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xyyyyz_0[i] = 3.0 * g_xxzz_0_xyyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyyz_1[i] * fz_be_0 + g_xxxzz_0_yyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyyzz_0[i] = 3.0 * g_xxzz_0_xyyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyzz_1[i] * fz_be_0 + g_xxxzz_0_yyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyzzz_0[i] = 3.0 * g_xxzz_0_xyyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyzzz_1[i] * fz_be_0 + g_xxxzz_0_yyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyzzzz_0[i] = 3.0 * g_xxzz_0_xyzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyzzzz_1[i] * fz_be_0 + g_xxxzz_0_yzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xzzzzz_0[i] = 3.0 * g_xxzz_0_xzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xzzzzz_1[i] * fz_be_0 + g_xxxzz_0_zzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyyy_0[i] = 3.0 * g_xxzz_0_yyyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyyy_1[i] * fz_be_0 + g_xxxzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyyz_0[i] = 3.0 * g_xxzz_0_yyyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyyz_1[i] * fz_be_0 + g_xxxzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyzz_0[i] = 3.0 * g_xxzz_0_yyyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyzz_1[i] * fz_be_0 + g_xxxzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyzzz_0[i] = 3.0 * g_xxzz_0_yyyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyzzz_1[i] * fz_be_0 + g_xxxzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyzzzz_0[i] = 3.0 * g_xxzz_0_yyzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyzzzz_1[i] * fz_be_0 + g_xxxzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yzzzzz_0[i] = 3.0 * g_xxzz_0_yzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yzzzzz_1[i] * fz_be_0 + g_xxxzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_zzzzzz_0[i] = 3.0 * g_xxzz_0_zzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_zzzzzz_1[i] * fz_be_0 + g_xxxzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 168-196 components of targeted buffer : ISI

    auto g_xxxyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 168);

    auto g_xxxyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 169);

    auto g_xxxyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 170);

    auto g_xxxyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 171);

    auto g_xxxyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 172);

    auto g_xxxyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 173);

    auto g_xxxyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 174);

    auto g_xxxyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 175);

    auto g_xxxyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 176);

    auto g_xxxyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 177);

    auto g_xxxyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 178);

    auto g_xxxyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 179);

    auto g_xxxyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 180);

    auto g_xxxyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 181);

    auto g_xxxyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 182);

    auto g_xxxyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 183);

    auto g_xxxyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 184);

    auto g_xxxyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 185);

    auto g_xxxyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 186);

    auto g_xxxyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 187);

    auto g_xxxyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 188);

    auto g_xxxyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 189);

    auto g_xxxyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 190);

    auto g_xxxyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 191);

    auto g_xxxyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 192);

    auto g_xxxyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 193);

    auto g_xxxyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 194);

    auto g_xxxyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 195);

    #pragma omp simd aligned(g_xxxy_0_xxxxxx_0, g_xxxy_0_xxxxxx_1, g_xxxy_0_xxxxxz_0, g_xxxy_0_xxxxxz_1, g_xxxy_0_xxxxzz_0, g_xxxy_0_xxxxzz_1, g_xxxy_0_xxxzzz_0, g_xxxy_0_xxxzzz_1, g_xxxy_0_xxzzzz_0, g_xxxy_0_xxzzzz_1, g_xxxy_0_xzzzzz_0, g_xxxy_0_xzzzzz_1, g_xxxyy_0_xxxxxx_1, g_xxxyy_0_xxxxxz_1, g_xxxyy_0_xxxxzz_1, g_xxxyy_0_xxxzzz_1, g_xxxyy_0_xxzzzz_1, g_xxxyy_0_xzzzzz_1, g_xxxyyy_0_xxxxxx_0, g_xxxyyy_0_xxxxxy_0, g_xxxyyy_0_xxxxxz_0, g_xxxyyy_0_xxxxyy_0, g_xxxyyy_0_xxxxyz_0, g_xxxyyy_0_xxxxzz_0, g_xxxyyy_0_xxxyyy_0, g_xxxyyy_0_xxxyyz_0, g_xxxyyy_0_xxxyzz_0, g_xxxyyy_0_xxxzzz_0, g_xxxyyy_0_xxyyyy_0, g_xxxyyy_0_xxyyyz_0, g_xxxyyy_0_xxyyzz_0, g_xxxyyy_0_xxyzzz_0, g_xxxyyy_0_xxzzzz_0, g_xxxyyy_0_xyyyyy_0, g_xxxyyy_0_xyyyyz_0, g_xxxyyy_0_xyyyzz_0, g_xxxyyy_0_xyyzzz_0, g_xxxyyy_0_xyzzzz_0, g_xxxyyy_0_xzzzzz_0, g_xxxyyy_0_yyyyyy_0, g_xxxyyy_0_yyyyyz_0, g_xxxyyy_0_yyyyzz_0, g_xxxyyy_0_yyyzzz_0, g_xxxyyy_0_yyzzzz_0, g_xxxyyy_0_yzzzzz_0, g_xxxyyy_0_zzzzzz_0, g_xxyyy_0_xxxxxy_1, g_xxyyy_0_xxxxy_1, g_xxyyy_0_xxxxyy_1, g_xxyyy_0_xxxxyz_1, g_xxyyy_0_xxxyy_1, g_xxyyy_0_xxxyyy_1, g_xxyyy_0_xxxyyz_1, g_xxyyy_0_xxxyz_1, g_xxyyy_0_xxxyzz_1, g_xxyyy_0_xxyyy_1, g_xxyyy_0_xxyyyy_1, g_xxyyy_0_xxyyyz_1, g_xxyyy_0_xxyyz_1, g_xxyyy_0_xxyyzz_1, g_xxyyy_0_xxyzz_1, g_xxyyy_0_xxyzzz_1, g_xxyyy_0_xyyyy_1, g_xxyyy_0_xyyyyy_1, g_xxyyy_0_xyyyyz_1, g_xxyyy_0_xyyyz_1, g_xxyyy_0_xyyyzz_1, g_xxyyy_0_xyyzz_1, g_xxyyy_0_xyyzzz_1, g_xxyyy_0_xyzzz_1, g_xxyyy_0_xyzzzz_1, g_xxyyy_0_yyyyy_1, g_xxyyy_0_yyyyyy_1, g_xxyyy_0_yyyyyz_1, g_xxyyy_0_yyyyz_1, g_xxyyy_0_yyyyzz_1, g_xxyyy_0_yyyzz_1, g_xxyyy_0_yyyzzz_1, g_xxyyy_0_yyzzz_1, g_xxyyy_0_yyzzzz_1, g_xxyyy_0_yzzzz_1, g_xxyyy_0_yzzzzz_1, g_xxyyy_0_zzzzzz_1, g_xyyy_0_xxxxxy_0, g_xyyy_0_xxxxxy_1, g_xyyy_0_xxxxyy_0, g_xyyy_0_xxxxyy_1, g_xyyy_0_xxxxyz_0, g_xyyy_0_xxxxyz_1, g_xyyy_0_xxxyyy_0, g_xyyy_0_xxxyyy_1, g_xyyy_0_xxxyyz_0, g_xyyy_0_xxxyyz_1, g_xyyy_0_xxxyzz_0, g_xyyy_0_xxxyzz_1, g_xyyy_0_xxyyyy_0, g_xyyy_0_xxyyyy_1, g_xyyy_0_xxyyyz_0, g_xyyy_0_xxyyyz_1, g_xyyy_0_xxyyzz_0, g_xyyy_0_xxyyzz_1, g_xyyy_0_xxyzzz_0, g_xyyy_0_xxyzzz_1, g_xyyy_0_xyyyyy_0, g_xyyy_0_xyyyyy_1, g_xyyy_0_xyyyyz_0, g_xyyy_0_xyyyyz_1, g_xyyy_0_xyyyzz_0, g_xyyy_0_xyyyzz_1, g_xyyy_0_xyyzzz_0, g_xyyy_0_xyyzzz_1, g_xyyy_0_xyzzzz_0, g_xyyy_0_xyzzzz_1, g_xyyy_0_yyyyyy_0, g_xyyy_0_yyyyyy_1, g_xyyy_0_yyyyyz_0, g_xyyy_0_yyyyyz_1, g_xyyy_0_yyyyzz_0, g_xyyy_0_yyyyzz_1, g_xyyy_0_yyyzzz_0, g_xyyy_0_yyyzzz_1, g_xyyy_0_yyzzzz_0, g_xyyy_0_yyzzzz_1, g_xyyy_0_yzzzzz_0, g_xyyy_0_yzzzzz_1, g_xyyy_0_zzzzzz_0, g_xyyy_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyy_0_xxxxxx_0[i] = 2.0 * g_xxxy_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxxx_1[i] * fz_be_0 + g_xxxyy_0_xxxxxx_1[i] * wa_y[i];

        g_xxxyyy_0_xxxxxy_0[i] = 2.0 * g_xyyy_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxyyy_0_xxxxy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxxz_0[i] = 2.0 * g_xxxy_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxxz_1[i] * fz_be_0 + g_xxxyy_0_xxxxxz_1[i] * wa_y[i];

        g_xxxyyy_0_xxxxyy_0[i] = 2.0 * g_xyyy_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxyyy_0_xxxyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxyz_0[i] = 2.0 * g_xyyy_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxyyy_0_xxxyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxzz_0[i] = 2.0 * g_xxxy_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxzz_1[i] * fz_be_0 + g_xxxyy_0_xxxxzz_1[i] * wa_y[i];

        g_xxxyyy_0_xxxyyy_0[i] = 2.0 * g_xyyy_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxyyz_0[i] = 2.0 * g_xyyy_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxyzz_0[i] = 2.0 * g_xyyy_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxzzz_0[i] = 2.0 * g_xxxy_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxzzz_1[i] * fz_be_0 + g_xxxyy_0_xxxzzz_1[i] * wa_y[i];

        g_xxxyyy_0_xxyyyy_0[i] = 2.0 * g_xyyy_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxyyyz_0[i] = 2.0 * g_xyyy_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxyyzz_0[i] = 2.0 * g_xyyy_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxyzzz_0[i] = 2.0 * g_xyyy_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxzzzz_0[i] = 2.0 * g_xxxy_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxzzzz_1[i] * fz_be_0 + g_xxxyy_0_xxzzzz_1[i] * wa_y[i];

        g_xxxyyy_0_xyyyyy_0[i] = 2.0 * g_xyyy_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyyy_1[i] * fz_be_0 + g_xxyyy_0_yyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xyyyyz_0[i] = 2.0 * g_xyyy_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyyz_1[i] * fz_be_0 + g_xxyyy_0_yyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xyyyzz_0[i] = 2.0 * g_xyyy_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyzz_1[i] * fz_be_0 + g_xxyyy_0_yyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xyyzzz_0[i] = 2.0 * g_xyyy_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyzzz_1[i] * fz_be_0 + g_xxyyy_0_yyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xyzzzz_0[i] = 2.0 * g_xyyy_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyzzzz_1[i] * fz_be_0 + g_xxyyy_0_yzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xzzzzz_0[i] = 2.0 * g_xxxy_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xzzzzz_1[i] * fz_be_0 + g_xxxyy_0_xzzzzz_1[i] * wa_y[i];

        g_xxxyyy_0_yyyyyy_0[i] = 2.0 * g_xyyy_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyyy_1[i] * fz_be_0 + g_xxyyy_0_yyyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_yyyyyz_0[i] = 2.0 * g_xyyy_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyyz_1[i] * fz_be_0 + g_xxyyy_0_yyyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_yyyyzz_0[i] = 2.0 * g_xyyy_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyzz_1[i] * fz_be_0 + g_xxyyy_0_yyyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_yyyzzz_0[i] = 2.0 * g_xyyy_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyzzz_1[i] * fz_be_0 + g_xxyyy_0_yyyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_yyzzzz_0[i] = 2.0 * g_xyyy_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyzzzz_1[i] * fz_be_0 + g_xxyyy_0_yyzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_yzzzzz_0[i] = 2.0 * g_xyyy_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yzzzzz_1[i] * fz_be_0 + g_xxyyy_0_yzzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_zzzzzz_0[i] = 2.0 * g_xyyy_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_zzzzzz_1[i] * fz_be_0 + g_xxyyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 196-224 components of targeted buffer : ISI

    auto g_xxxyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 196);

    auto g_xxxyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 197);

    auto g_xxxyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 198);

    auto g_xxxyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 199);

    auto g_xxxyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 200);

    auto g_xxxyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 201);

    auto g_xxxyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 202);

    auto g_xxxyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 203);

    auto g_xxxyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 204);

    auto g_xxxyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 205);

    auto g_xxxyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 206);

    auto g_xxxyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 207);

    auto g_xxxyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 208);

    auto g_xxxyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 209);

    auto g_xxxyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 210);

    auto g_xxxyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 211);

    auto g_xxxyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 212);

    auto g_xxxyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 213);

    auto g_xxxyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 214);

    auto g_xxxyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 215);

    auto g_xxxyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 216);

    auto g_xxxyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 217);

    auto g_xxxyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 218);

    auto g_xxxyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 219);

    auto g_xxxyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 220);

    auto g_xxxyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 221);

    auto g_xxxyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 222);

    auto g_xxxyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 223);

    #pragma omp simd aligned(g_xxxyy_0_xxxxx_1, g_xxxyy_0_xxxxxx_1, g_xxxyy_0_xxxxxy_1, g_xxxyy_0_xxxxxz_1, g_xxxyy_0_xxxxy_1, g_xxxyy_0_xxxxyy_1, g_xxxyy_0_xxxxyz_1, g_xxxyy_0_xxxxz_1, g_xxxyy_0_xxxxzz_1, g_xxxyy_0_xxxyy_1, g_xxxyy_0_xxxyyy_1, g_xxxyy_0_xxxyyz_1, g_xxxyy_0_xxxyz_1, g_xxxyy_0_xxxyzz_1, g_xxxyy_0_xxxzz_1, g_xxxyy_0_xxxzzz_1, g_xxxyy_0_xxyyy_1, g_xxxyy_0_xxyyyy_1, g_xxxyy_0_xxyyyz_1, g_xxxyy_0_xxyyz_1, g_xxxyy_0_xxyyzz_1, g_xxxyy_0_xxyzz_1, g_xxxyy_0_xxyzzz_1, g_xxxyy_0_xxzzz_1, g_xxxyy_0_xxzzzz_1, g_xxxyy_0_xyyyy_1, g_xxxyy_0_xyyyyy_1, g_xxxyy_0_xyyyyz_1, g_xxxyy_0_xyyyz_1, g_xxxyy_0_xyyyzz_1, g_xxxyy_0_xyyzz_1, g_xxxyy_0_xyyzzz_1, g_xxxyy_0_xyzzz_1, g_xxxyy_0_xyzzzz_1, g_xxxyy_0_xzzzz_1, g_xxxyy_0_xzzzzz_1, g_xxxyy_0_yyyyy_1, g_xxxyy_0_yyyyyy_1, g_xxxyy_0_yyyyyz_1, g_xxxyy_0_yyyyz_1, g_xxxyy_0_yyyyzz_1, g_xxxyy_0_yyyzz_1, g_xxxyy_0_yyyzzz_1, g_xxxyy_0_yyzzz_1, g_xxxyy_0_yyzzzz_1, g_xxxyy_0_yzzzz_1, g_xxxyy_0_yzzzzz_1, g_xxxyy_0_zzzzz_1, g_xxxyy_0_zzzzzz_1, g_xxxyyz_0_xxxxxx_0, g_xxxyyz_0_xxxxxy_0, g_xxxyyz_0_xxxxxz_0, g_xxxyyz_0_xxxxyy_0, g_xxxyyz_0_xxxxyz_0, g_xxxyyz_0_xxxxzz_0, g_xxxyyz_0_xxxyyy_0, g_xxxyyz_0_xxxyyz_0, g_xxxyyz_0_xxxyzz_0, g_xxxyyz_0_xxxzzz_0, g_xxxyyz_0_xxyyyy_0, g_xxxyyz_0_xxyyyz_0, g_xxxyyz_0_xxyyzz_0, g_xxxyyz_0_xxyzzz_0, g_xxxyyz_0_xxzzzz_0, g_xxxyyz_0_xyyyyy_0, g_xxxyyz_0_xyyyyz_0, g_xxxyyz_0_xyyyzz_0, g_xxxyyz_0_xyyzzz_0, g_xxxyyz_0_xyzzzz_0, g_xxxyyz_0_xzzzzz_0, g_xxxyyz_0_yyyyyy_0, g_xxxyyz_0_yyyyyz_0, g_xxxyyz_0_yyyyzz_0, g_xxxyyz_0_yyyzzz_0, g_xxxyyz_0_yyzzzz_0, g_xxxyyz_0_yzzzzz_0, g_xxxyyz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyz_0_xxxxxx_0[i] = g_xxxyy_0_xxxxxx_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxy_0[i] = g_xxxyy_0_xxxxxy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxz_0[i] = g_xxxyy_0_xxxxx_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxyy_0[i] = g_xxxyy_0_xxxxyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxyz_0[i] = g_xxxyy_0_xxxxy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxzz_0[i] = 2.0 * g_xxxyy_0_xxxxz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyyy_0[i] = g_xxxyy_0_xxxyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyyz_0[i] = g_xxxyy_0_xxxyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyzz_0[i] = 2.0 * g_xxxyy_0_xxxyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxzzz_0[i] = 3.0 * g_xxxyy_0_xxxzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyyy_0[i] = g_xxxyy_0_xxyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyyz_0[i] = g_xxxyy_0_xxyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyzz_0[i] = 2.0 * g_xxxyy_0_xxyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyzzz_0[i] = 3.0 * g_xxxyy_0_xxyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxzzzz_0[i] = 4.0 * g_xxxyy_0_xxzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyyy_0[i] = g_xxxyy_0_xyyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyyz_0[i] = g_xxxyy_0_xyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyzz_0[i] = 2.0 * g_xxxyy_0_xyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyzzz_0[i] = 3.0 * g_xxxyy_0_xyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyzzzz_0[i] = 4.0 * g_xxxyy_0_xyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xzzzzz_0[i] = 5.0 * g_xxxyy_0_xzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyyy_0[i] = g_xxxyy_0_yyyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyyz_0[i] = g_xxxyy_0_yyyyy_1[i] * fi_acd_0 + g_xxxyy_0_yyyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyzz_0[i] = 2.0 * g_xxxyy_0_yyyyz_1[i] * fi_acd_0 + g_xxxyy_0_yyyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyzzz_0[i] = 3.0 * g_xxxyy_0_yyyzz_1[i] * fi_acd_0 + g_xxxyy_0_yyyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyzzzz_0[i] = 4.0 * g_xxxyy_0_yyzzz_1[i] * fi_acd_0 + g_xxxyy_0_yyzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yzzzzz_0[i] = 5.0 * g_xxxyy_0_yzzzz_1[i] * fi_acd_0 + g_xxxyy_0_yzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_zzzzzz_0[i] = 6.0 * g_xxxyy_0_zzzzz_1[i] * fi_acd_0 + g_xxxyy_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 224-252 components of targeted buffer : ISI

    auto g_xxxyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 224);

    auto g_xxxyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 225);

    auto g_xxxyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 226);

    auto g_xxxyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 227);

    auto g_xxxyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 228);

    auto g_xxxyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 229);

    auto g_xxxyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 230);

    auto g_xxxyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 231);

    auto g_xxxyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 232);

    auto g_xxxyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 233);

    auto g_xxxyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 234);

    auto g_xxxyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 235);

    auto g_xxxyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 236);

    auto g_xxxyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 237);

    auto g_xxxyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 238);

    auto g_xxxyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 239);

    auto g_xxxyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 240);

    auto g_xxxyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 241);

    auto g_xxxyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 242);

    auto g_xxxyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 243);

    auto g_xxxyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 244);

    auto g_xxxyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 245);

    auto g_xxxyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 246);

    auto g_xxxyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 247);

    auto g_xxxyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 248);

    auto g_xxxyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 249);

    auto g_xxxyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 250);

    auto g_xxxyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 251);

    #pragma omp simd aligned(g_xxxyzz_0_xxxxxx_0, g_xxxyzz_0_xxxxxy_0, g_xxxyzz_0_xxxxxz_0, g_xxxyzz_0_xxxxyy_0, g_xxxyzz_0_xxxxyz_0, g_xxxyzz_0_xxxxzz_0, g_xxxyzz_0_xxxyyy_0, g_xxxyzz_0_xxxyyz_0, g_xxxyzz_0_xxxyzz_0, g_xxxyzz_0_xxxzzz_0, g_xxxyzz_0_xxyyyy_0, g_xxxyzz_0_xxyyyz_0, g_xxxyzz_0_xxyyzz_0, g_xxxyzz_0_xxyzzz_0, g_xxxyzz_0_xxzzzz_0, g_xxxyzz_0_xyyyyy_0, g_xxxyzz_0_xyyyyz_0, g_xxxyzz_0_xyyyzz_0, g_xxxyzz_0_xyyzzz_0, g_xxxyzz_0_xyzzzz_0, g_xxxyzz_0_xzzzzz_0, g_xxxyzz_0_yyyyyy_0, g_xxxyzz_0_yyyyyz_0, g_xxxyzz_0_yyyyzz_0, g_xxxyzz_0_yyyzzz_0, g_xxxyzz_0_yyzzzz_0, g_xxxyzz_0_yzzzzz_0, g_xxxyzz_0_zzzzzz_0, g_xxxzz_0_xxxxx_1, g_xxxzz_0_xxxxxx_1, g_xxxzz_0_xxxxxy_1, g_xxxzz_0_xxxxxz_1, g_xxxzz_0_xxxxy_1, g_xxxzz_0_xxxxyy_1, g_xxxzz_0_xxxxyz_1, g_xxxzz_0_xxxxz_1, g_xxxzz_0_xxxxzz_1, g_xxxzz_0_xxxyy_1, g_xxxzz_0_xxxyyy_1, g_xxxzz_0_xxxyyz_1, g_xxxzz_0_xxxyz_1, g_xxxzz_0_xxxyzz_1, g_xxxzz_0_xxxzz_1, g_xxxzz_0_xxxzzz_1, g_xxxzz_0_xxyyy_1, g_xxxzz_0_xxyyyy_1, g_xxxzz_0_xxyyyz_1, g_xxxzz_0_xxyyz_1, g_xxxzz_0_xxyyzz_1, g_xxxzz_0_xxyzz_1, g_xxxzz_0_xxyzzz_1, g_xxxzz_0_xxzzz_1, g_xxxzz_0_xxzzzz_1, g_xxxzz_0_xyyyy_1, g_xxxzz_0_xyyyyy_1, g_xxxzz_0_xyyyyz_1, g_xxxzz_0_xyyyz_1, g_xxxzz_0_xyyyzz_1, g_xxxzz_0_xyyzz_1, g_xxxzz_0_xyyzzz_1, g_xxxzz_0_xyzzz_1, g_xxxzz_0_xyzzzz_1, g_xxxzz_0_xzzzz_1, g_xxxzz_0_xzzzzz_1, g_xxxzz_0_yyyyy_1, g_xxxzz_0_yyyyyy_1, g_xxxzz_0_yyyyyz_1, g_xxxzz_0_yyyyz_1, g_xxxzz_0_yyyyzz_1, g_xxxzz_0_yyyzz_1, g_xxxzz_0_yyyzzz_1, g_xxxzz_0_yyzzz_1, g_xxxzz_0_yyzzzz_1, g_xxxzz_0_yzzzz_1, g_xxxzz_0_yzzzzz_1, g_xxxzz_0_zzzzz_1, g_xxxzz_0_zzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzz_0_xxxxxx_0[i] = g_xxxzz_0_xxxxxx_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxy_0[i] = g_xxxzz_0_xxxxx_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxz_0[i] = g_xxxzz_0_xxxxxz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxyy_0[i] = 2.0 * g_xxxzz_0_xxxxy_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxyz_0[i] = g_xxxzz_0_xxxxz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxzz_0[i] = g_xxxzz_0_xxxxzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyyy_0[i] = 3.0 * g_xxxzz_0_xxxyy_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyyz_0[i] = 2.0 * g_xxxzz_0_xxxyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyzz_0[i] = g_xxxzz_0_xxxzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxzzz_0[i] = g_xxxzz_0_xxxzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyyy_0[i] = 4.0 * g_xxxzz_0_xxyyy_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyyz_0[i] = 3.0 * g_xxxzz_0_xxyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyzz_0[i] = 2.0 * g_xxxzz_0_xxyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyzzz_0[i] = g_xxxzz_0_xxzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxzzzz_0[i] = g_xxxzz_0_xxzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyyy_0[i] = 5.0 * g_xxxzz_0_xyyyy_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyyz_0[i] = 4.0 * g_xxxzz_0_xyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyzz_0[i] = 3.0 * g_xxxzz_0_xyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyzzz_0[i] = 2.0 * g_xxxzz_0_xyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyzzzz_0[i] = g_xxxzz_0_xzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xzzzzz_0[i] = g_xxxzz_0_xzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyyy_0[i] = 6.0 * g_xxxzz_0_yyyyy_1[i] * fi_acd_0 + g_xxxzz_0_yyyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyyz_0[i] = 5.0 * g_xxxzz_0_yyyyz_1[i] * fi_acd_0 + g_xxxzz_0_yyyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyzz_0[i] = 4.0 * g_xxxzz_0_yyyzz_1[i] * fi_acd_0 + g_xxxzz_0_yyyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyzzz_0[i] = 3.0 * g_xxxzz_0_yyzzz_1[i] * fi_acd_0 + g_xxxzz_0_yyyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyzzzz_0[i] = 2.0 * g_xxxzz_0_yzzzz_1[i] * fi_acd_0 + g_xxxzz_0_yyzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yzzzzz_0[i] = g_xxxzz_0_zzzzz_1[i] * fi_acd_0 + g_xxxzz_0_yzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_zzzzzz_0[i] = g_xxxzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 252-280 components of targeted buffer : ISI

    auto g_xxxzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 252);

    auto g_xxxzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 253);

    auto g_xxxzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 254);

    auto g_xxxzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 255);

    auto g_xxxzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 256);

    auto g_xxxzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 257);

    auto g_xxxzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 258);

    auto g_xxxzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 259);

    auto g_xxxzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 260);

    auto g_xxxzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 261);

    auto g_xxxzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 262);

    auto g_xxxzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 263);

    auto g_xxxzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 264);

    auto g_xxxzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 265);

    auto g_xxxzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 266);

    auto g_xxxzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 267);

    auto g_xxxzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 268);

    auto g_xxxzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 269);

    auto g_xxxzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 270);

    auto g_xxxzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 271);

    auto g_xxxzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 272);

    auto g_xxxzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 273);

    auto g_xxxzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 274);

    auto g_xxxzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 275);

    auto g_xxxzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 276);

    auto g_xxxzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 277);

    auto g_xxxzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 278);

    auto g_xxxzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 279);

    #pragma omp simd aligned(g_xxxz_0_xxxxxx_0, g_xxxz_0_xxxxxx_1, g_xxxz_0_xxxxxy_0, g_xxxz_0_xxxxxy_1, g_xxxz_0_xxxxyy_0, g_xxxz_0_xxxxyy_1, g_xxxz_0_xxxyyy_0, g_xxxz_0_xxxyyy_1, g_xxxz_0_xxyyyy_0, g_xxxz_0_xxyyyy_1, g_xxxz_0_xyyyyy_0, g_xxxz_0_xyyyyy_1, g_xxxzz_0_xxxxxx_1, g_xxxzz_0_xxxxxy_1, g_xxxzz_0_xxxxyy_1, g_xxxzz_0_xxxyyy_1, g_xxxzz_0_xxyyyy_1, g_xxxzz_0_xyyyyy_1, g_xxxzzz_0_xxxxxx_0, g_xxxzzz_0_xxxxxy_0, g_xxxzzz_0_xxxxxz_0, g_xxxzzz_0_xxxxyy_0, g_xxxzzz_0_xxxxyz_0, g_xxxzzz_0_xxxxzz_0, g_xxxzzz_0_xxxyyy_0, g_xxxzzz_0_xxxyyz_0, g_xxxzzz_0_xxxyzz_0, g_xxxzzz_0_xxxzzz_0, g_xxxzzz_0_xxyyyy_0, g_xxxzzz_0_xxyyyz_0, g_xxxzzz_0_xxyyzz_0, g_xxxzzz_0_xxyzzz_0, g_xxxzzz_0_xxzzzz_0, g_xxxzzz_0_xyyyyy_0, g_xxxzzz_0_xyyyyz_0, g_xxxzzz_0_xyyyzz_0, g_xxxzzz_0_xyyzzz_0, g_xxxzzz_0_xyzzzz_0, g_xxxzzz_0_xzzzzz_0, g_xxxzzz_0_yyyyyy_0, g_xxxzzz_0_yyyyyz_0, g_xxxzzz_0_yyyyzz_0, g_xxxzzz_0_yyyzzz_0, g_xxxzzz_0_yyzzzz_0, g_xxxzzz_0_yzzzzz_0, g_xxxzzz_0_zzzzzz_0, g_xxzzz_0_xxxxxz_1, g_xxzzz_0_xxxxyz_1, g_xxzzz_0_xxxxz_1, g_xxzzz_0_xxxxzz_1, g_xxzzz_0_xxxyyz_1, g_xxzzz_0_xxxyz_1, g_xxzzz_0_xxxyzz_1, g_xxzzz_0_xxxzz_1, g_xxzzz_0_xxxzzz_1, g_xxzzz_0_xxyyyz_1, g_xxzzz_0_xxyyz_1, g_xxzzz_0_xxyyzz_1, g_xxzzz_0_xxyzz_1, g_xxzzz_0_xxyzzz_1, g_xxzzz_0_xxzzz_1, g_xxzzz_0_xxzzzz_1, g_xxzzz_0_xyyyyz_1, g_xxzzz_0_xyyyz_1, g_xxzzz_0_xyyyzz_1, g_xxzzz_0_xyyzz_1, g_xxzzz_0_xyyzzz_1, g_xxzzz_0_xyzzz_1, g_xxzzz_0_xyzzzz_1, g_xxzzz_0_xzzzz_1, g_xxzzz_0_xzzzzz_1, g_xxzzz_0_yyyyyy_1, g_xxzzz_0_yyyyyz_1, g_xxzzz_0_yyyyz_1, g_xxzzz_0_yyyyzz_1, g_xxzzz_0_yyyzz_1, g_xxzzz_0_yyyzzz_1, g_xxzzz_0_yyzzz_1, g_xxzzz_0_yyzzzz_1, g_xxzzz_0_yzzzz_1, g_xxzzz_0_yzzzzz_1, g_xxzzz_0_zzzzz_1, g_xxzzz_0_zzzzzz_1, g_xzzz_0_xxxxxz_0, g_xzzz_0_xxxxxz_1, g_xzzz_0_xxxxyz_0, g_xzzz_0_xxxxyz_1, g_xzzz_0_xxxxzz_0, g_xzzz_0_xxxxzz_1, g_xzzz_0_xxxyyz_0, g_xzzz_0_xxxyyz_1, g_xzzz_0_xxxyzz_0, g_xzzz_0_xxxyzz_1, g_xzzz_0_xxxzzz_0, g_xzzz_0_xxxzzz_1, g_xzzz_0_xxyyyz_0, g_xzzz_0_xxyyyz_1, g_xzzz_0_xxyyzz_0, g_xzzz_0_xxyyzz_1, g_xzzz_0_xxyzzz_0, g_xzzz_0_xxyzzz_1, g_xzzz_0_xxzzzz_0, g_xzzz_0_xxzzzz_1, g_xzzz_0_xyyyyz_0, g_xzzz_0_xyyyyz_1, g_xzzz_0_xyyyzz_0, g_xzzz_0_xyyyzz_1, g_xzzz_0_xyyzzz_0, g_xzzz_0_xyyzzz_1, g_xzzz_0_xyzzzz_0, g_xzzz_0_xyzzzz_1, g_xzzz_0_xzzzzz_0, g_xzzz_0_xzzzzz_1, g_xzzz_0_yyyyyy_0, g_xzzz_0_yyyyyy_1, g_xzzz_0_yyyyyz_0, g_xzzz_0_yyyyyz_1, g_xzzz_0_yyyyzz_0, g_xzzz_0_yyyyzz_1, g_xzzz_0_yyyzzz_0, g_xzzz_0_yyyzzz_1, g_xzzz_0_yyzzzz_0, g_xzzz_0_yyzzzz_1, g_xzzz_0_yzzzzz_0, g_xzzz_0_yzzzzz_1, g_xzzz_0_zzzzzz_0, g_xzzz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzz_0_xxxxxx_0[i] = 2.0 * g_xxxz_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxxx_1[i] * fz_be_0 + g_xxxzz_0_xxxxxx_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxxy_0[i] = 2.0 * g_xxxz_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxxy_1[i] * fz_be_0 + g_xxxzz_0_xxxxxy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxxz_0[i] = 2.0 * g_xzzz_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxzzz_0_xxxxz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxyy_0[i] = 2.0 * g_xxxz_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxyy_1[i] * fz_be_0 + g_xxxzz_0_xxxxyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxyz_0[i] = 2.0 * g_xzzz_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxzzz_0_xxxyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxzz_0[i] = 2.0 * g_xzzz_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxzzz_0_xxxzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxyyy_0[i] = 2.0 * g_xxxz_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxyyy_1[i] * fz_be_0 + g_xxxzz_0_xxxyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxyyz_0[i] = 2.0 * g_xzzz_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxyzz_0[i] = 2.0 * g_xzzz_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxzzz_0[i] = 2.0 * g_xzzz_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyyyy_0[i] = 2.0 * g_xxxz_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxyyyy_1[i] * fz_be_0 + g_xxxzz_0_xxyyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxyyyz_0[i] = 2.0 * g_xzzz_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyyzz_0[i] = 2.0 * g_xzzz_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyzzz_0[i] = 2.0 * g_xzzz_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxzzzz_0[i] = 2.0 * g_xzzz_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyyyy_0[i] = 2.0 * g_xxxz_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xyyyyy_1[i] * fz_be_0 + g_xxxzz_0_xyyyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xyyyyz_0[i] = 2.0 * g_xzzz_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyyyz_1[i] * fz_be_0 + g_xxzzz_0_yyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyyzz_0[i] = 2.0 * g_xzzz_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyyzz_1[i] * fz_be_0 + g_xxzzz_0_yyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyzzz_0[i] = 2.0 * g_xzzz_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyzzz_1[i] * fz_be_0 + g_xxzzz_0_yyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyzzzz_0[i] = 2.0 * g_xzzz_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyzzzz_1[i] * fz_be_0 + g_xxzzz_0_yzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xzzzzz_0[i] = 2.0 * g_xzzz_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xzzzzz_1[i] * fz_be_0 + g_xxzzz_0_zzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyyy_0[i] = 2.0 * g_xzzz_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyyy_1[i] * fz_be_0 + g_xxzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyyz_0[i] = 2.0 * g_xzzz_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyyz_1[i] * fz_be_0 + g_xxzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyzz_0[i] = 2.0 * g_xzzz_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyzz_1[i] * fz_be_0 + g_xxzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyzzz_0[i] = 2.0 * g_xzzz_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyzzz_1[i] * fz_be_0 + g_xxzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyzzzz_0[i] = 2.0 * g_xzzz_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyzzzz_1[i] * fz_be_0 + g_xxzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yzzzzz_0[i] = 2.0 * g_xzzz_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yzzzzz_1[i] * fz_be_0 + g_xxzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_zzzzzz_0[i] = 2.0 * g_xzzz_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_zzzzzz_1[i] * fz_be_0 + g_xxzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 280-308 components of targeted buffer : ISI

    auto g_xxyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 280);

    auto g_xxyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 281);

    auto g_xxyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 282);

    auto g_xxyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 283);

    auto g_xxyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 284);

    auto g_xxyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 285);

    auto g_xxyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 286);

    auto g_xxyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 287);

    auto g_xxyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 288);

    auto g_xxyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 289);

    auto g_xxyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 290);

    auto g_xxyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 291);

    auto g_xxyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 292);

    auto g_xxyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 293);

    auto g_xxyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 294);

    auto g_xxyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 295);

    auto g_xxyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 296);

    auto g_xxyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 297);

    auto g_xxyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 298);

    auto g_xxyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 299);

    auto g_xxyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 300);

    auto g_xxyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 301);

    auto g_xxyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 302);

    auto g_xxyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 303);

    auto g_xxyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 304);

    auto g_xxyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 305);

    auto g_xxyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 306);

    auto g_xxyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 307);

    #pragma omp simd aligned(g_xxyy_0_xxxxxx_0, g_xxyy_0_xxxxxx_1, g_xxyy_0_xxxxxz_0, g_xxyy_0_xxxxxz_1, g_xxyy_0_xxxxzz_0, g_xxyy_0_xxxxzz_1, g_xxyy_0_xxxzzz_0, g_xxyy_0_xxxzzz_1, g_xxyy_0_xxzzzz_0, g_xxyy_0_xxzzzz_1, g_xxyy_0_xzzzzz_0, g_xxyy_0_xzzzzz_1, g_xxyyy_0_xxxxxx_1, g_xxyyy_0_xxxxxz_1, g_xxyyy_0_xxxxzz_1, g_xxyyy_0_xxxzzz_1, g_xxyyy_0_xxzzzz_1, g_xxyyy_0_xzzzzz_1, g_xxyyyy_0_xxxxxx_0, g_xxyyyy_0_xxxxxy_0, g_xxyyyy_0_xxxxxz_0, g_xxyyyy_0_xxxxyy_0, g_xxyyyy_0_xxxxyz_0, g_xxyyyy_0_xxxxzz_0, g_xxyyyy_0_xxxyyy_0, g_xxyyyy_0_xxxyyz_0, g_xxyyyy_0_xxxyzz_0, g_xxyyyy_0_xxxzzz_0, g_xxyyyy_0_xxyyyy_0, g_xxyyyy_0_xxyyyz_0, g_xxyyyy_0_xxyyzz_0, g_xxyyyy_0_xxyzzz_0, g_xxyyyy_0_xxzzzz_0, g_xxyyyy_0_xyyyyy_0, g_xxyyyy_0_xyyyyz_0, g_xxyyyy_0_xyyyzz_0, g_xxyyyy_0_xyyzzz_0, g_xxyyyy_0_xyzzzz_0, g_xxyyyy_0_xzzzzz_0, g_xxyyyy_0_yyyyyy_0, g_xxyyyy_0_yyyyyz_0, g_xxyyyy_0_yyyyzz_0, g_xxyyyy_0_yyyzzz_0, g_xxyyyy_0_yyzzzz_0, g_xxyyyy_0_yzzzzz_0, g_xxyyyy_0_zzzzzz_0, g_xyyyy_0_xxxxxy_1, g_xyyyy_0_xxxxy_1, g_xyyyy_0_xxxxyy_1, g_xyyyy_0_xxxxyz_1, g_xyyyy_0_xxxyy_1, g_xyyyy_0_xxxyyy_1, g_xyyyy_0_xxxyyz_1, g_xyyyy_0_xxxyz_1, g_xyyyy_0_xxxyzz_1, g_xyyyy_0_xxyyy_1, g_xyyyy_0_xxyyyy_1, g_xyyyy_0_xxyyyz_1, g_xyyyy_0_xxyyz_1, g_xyyyy_0_xxyyzz_1, g_xyyyy_0_xxyzz_1, g_xyyyy_0_xxyzzz_1, g_xyyyy_0_xyyyy_1, g_xyyyy_0_xyyyyy_1, g_xyyyy_0_xyyyyz_1, g_xyyyy_0_xyyyz_1, g_xyyyy_0_xyyyzz_1, g_xyyyy_0_xyyzz_1, g_xyyyy_0_xyyzzz_1, g_xyyyy_0_xyzzz_1, g_xyyyy_0_xyzzzz_1, g_xyyyy_0_yyyyy_1, g_xyyyy_0_yyyyyy_1, g_xyyyy_0_yyyyyz_1, g_xyyyy_0_yyyyz_1, g_xyyyy_0_yyyyzz_1, g_xyyyy_0_yyyzz_1, g_xyyyy_0_yyyzzz_1, g_xyyyy_0_yyzzz_1, g_xyyyy_0_yyzzzz_1, g_xyyyy_0_yzzzz_1, g_xyyyy_0_yzzzzz_1, g_xyyyy_0_zzzzzz_1, g_yyyy_0_xxxxxy_0, g_yyyy_0_xxxxxy_1, g_yyyy_0_xxxxyy_0, g_yyyy_0_xxxxyy_1, g_yyyy_0_xxxxyz_0, g_yyyy_0_xxxxyz_1, g_yyyy_0_xxxyyy_0, g_yyyy_0_xxxyyy_1, g_yyyy_0_xxxyyz_0, g_yyyy_0_xxxyyz_1, g_yyyy_0_xxxyzz_0, g_yyyy_0_xxxyzz_1, g_yyyy_0_xxyyyy_0, g_yyyy_0_xxyyyy_1, g_yyyy_0_xxyyyz_0, g_yyyy_0_xxyyyz_1, g_yyyy_0_xxyyzz_0, g_yyyy_0_xxyyzz_1, g_yyyy_0_xxyzzz_0, g_yyyy_0_xxyzzz_1, g_yyyy_0_xyyyyy_0, g_yyyy_0_xyyyyy_1, g_yyyy_0_xyyyyz_0, g_yyyy_0_xyyyyz_1, g_yyyy_0_xyyyzz_0, g_yyyy_0_xyyyzz_1, g_yyyy_0_xyyzzz_0, g_yyyy_0_xyyzzz_1, g_yyyy_0_xyzzzz_0, g_yyyy_0_xyzzzz_1, g_yyyy_0_yyyyyy_0, g_yyyy_0_yyyyyy_1, g_yyyy_0_yyyyyz_0, g_yyyy_0_yyyyyz_1, g_yyyy_0_yyyyzz_0, g_yyyy_0_yyyyzz_1, g_yyyy_0_yyyzzz_0, g_yyyy_0_yyyzzz_1, g_yyyy_0_yyzzzz_0, g_yyyy_0_yyzzzz_1, g_yyyy_0_yzzzzz_0, g_yyyy_0_yzzzzz_1, g_yyyy_0_zzzzzz_0, g_yyyy_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyy_0_xxxxxx_0[i] = 3.0 * g_xxyy_0_xxxxxx_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxx_1[i] * fz_be_0 + g_xxyyy_0_xxxxxx_1[i] * wa_y[i];

        g_xxyyyy_0_xxxxxy_0[i] = g_yyyy_0_xxxxxy_0[i] * fbe_0 - g_yyyy_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xyyyy_0_xxxxy_1[i] * fi_acd_0 + g_xyyyy_0_xxxxxy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxxz_0[i] = 3.0 * g_xxyy_0_xxxxxz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxz_1[i] * fz_be_0 + g_xxyyy_0_xxxxxz_1[i] * wa_y[i];

        g_xxyyyy_0_xxxxyy_0[i] = g_yyyy_0_xxxxyy_0[i] * fbe_0 - g_yyyy_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xyyyy_0_xxxyy_1[i] * fi_acd_0 + g_xyyyy_0_xxxxyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxyz_0[i] = g_yyyy_0_xxxxyz_0[i] * fbe_0 - g_yyyy_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyyyy_0_xxxyz_1[i] * fi_acd_0 + g_xyyyy_0_xxxxyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxzz_0[i] = 3.0 * g_xxyy_0_xxxxzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxzz_1[i] * fz_be_0 + g_xxyyy_0_xxxxzz_1[i] * wa_y[i];

        g_xxyyyy_0_xxxyyy_0[i] = g_yyyy_0_xxxyyy_0[i] * fbe_0 - g_yyyy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyyy_1[i] * fi_acd_0 + g_xyyyy_0_xxxyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxyyz_0[i] = g_yyyy_0_xxxyyz_0[i] * fbe_0 - g_yyyy_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyyz_1[i] * fi_acd_0 + g_xyyyy_0_xxxyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxyzz_0[i] = g_yyyy_0_xxxyzz_0[i] * fbe_0 - g_yyyy_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyzz_1[i] * fi_acd_0 + g_xyyyy_0_xxxyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxzzz_0[i] = 3.0 * g_xxyy_0_xxxzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxzzz_1[i] * fz_be_0 + g_xxyyy_0_xxxzzz_1[i] * wa_y[i];

        g_xxyyyy_0_xxyyyy_0[i] = g_yyyy_0_xxyyyy_0[i] * fbe_0 - g_yyyy_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyyy_1[i] * fi_acd_0 + g_xyyyy_0_xxyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxyyyz_0[i] = g_yyyy_0_xxyyyz_0[i] * fbe_0 - g_yyyy_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyyz_1[i] * fi_acd_0 + g_xyyyy_0_xxyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxyyzz_0[i] = g_yyyy_0_xxyyzz_0[i] * fbe_0 - g_yyyy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyzz_1[i] * fi_acd_0 + g_xyyyy_0_xxyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxyzzz_0[i] = g_yyyy_0_xxyzzz_0[i] * fbe_0 - g_yyyy_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyzzz_1[i] * fi_acd_0 + g_xyyyy_0_xxyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxzzzz_0[i] = 3.0 * g_xxyy_0_xxzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxzzzz_1[i] * fz_be_0 + g_xxyyy_0_xxzzzz_1[i] * wa_y[i];

        g_xxyyyy_0_xyyyyy_0[i] = g_yyyy_0_xyyyyy_0[i] * fbe_0 - g_yyyy_0_xyyyyy_1[i] * fz_be_0 + g_xyyyy_0_yyyyy_1[i] * fi_acd_0 + g_xyyyy_0_xyyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xyyyyz_0[i] = g_yyyy_0_xyyyyz_0[i] * fbe_0 - g_yyyy_0_xyyyyz_1[i] * fz_be_0 + g_xyyyy_0_yyyyz_1[i] * fi_acd_0 + g_xyyyy_0_xyyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xyyyzz_0[i] = g_yyyy_0_xyyyzz_0[i] * fbe_0 - g_yyyy_0_xyyyzz_1[i] * fz_be_0 + g_xyyyy_0_yyyzz_1[i] * fi_acd_0 + g_xyyyy_0_xyyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xyyzzz_0[i] = g_yyyy_0_xyyzzz_0[i] * fbe_0 - g_yyyy_0_xyyzzz_1[i] * fz_be_0 + g_xyyyy_0_yyzzz_1[i] * fi_acd_0 + g_xyyyy_0_xyyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xyzzzz_0[i] = g_yyyy_0_xyzzzz_0[i] * fbe_0 - g_yyyy_0_xyzzzz_1[i] * fz_be_0 + g_xyyyy_0_yzzzz_1[i] * fi_acd_0 + g_xyyyy_0_xyzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xzzzzz_0[i] = 3.0 * g_xxyy_0_xzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xzzzzz_1[i] * fz_be_0 + g_xxyyy_0_xzzzzz_1[i] * wa_y[i];

        g_xxyyyy_0_yyyyyy_0[i] = g_yyyy_0_yyyyyy_0[i] * fbe_0 - g_yyyy_0_yyyyyy_1[i] * fz_be_0 + g_xyyyy_0_yyyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_yyyyyz_0[i] = g_yyyy_0_yyyyyz_0[i] * fbe_0 - g_yyyy_0_yyyyyz_1[i] * fz_be_0 + g_xyyyy_0_yyyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_yyyyzz_0[i] = g_yyyy_0_yyyyzz_0[i] * fbe_0 - g_yyyy_0_yyyyzz_1[i] * fz_be_0 + g_xyyyy_0_yyyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_yyyzzz_0[i] = g_yyyy_0_yyyzzz_0[i] * fbe_0 - g_yyyy_0_yyyzzz_1[i] * fz_be_0 + g_xyyyy_0_yyyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_yyzzzz_0[i] = g_yyyy_0_yyzzzz_0[i] * fbe_0 - g_yyyy_0_yyzzzz_1[i] * fz_be_0 + g_xyyyy_0_yyzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_yzzzzz_0[i] = g_yyyy_0_yzzzzz_0[i] * fbe_0 - g_yyyy_0_yzzzzz_1[i] * fz_be_0 + g_xyyyy_0_yzzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_zzzzzz_0[i] = g_yyyy_0_zzzzzz_0[i] * fbe_0 - g_yyyy_0_zzzzzz_1[i] * fz_be_0 + g_xyyyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 308-336 components of targeted buffer : ISI

    auto g_xxyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 308);

    auto g_xxyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 309);

    auto g_xxyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 310);

    auto g_xxyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 311);

    auto g_xxyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 312);

    auto g_xxyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 313);

    auto g_xxyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 314);

    auto g_xxyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 315);

    auto g_xxyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 316);

    auto g_xxyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 317);

    auto g_xxyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 318);

    auto g_xxyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 319);

    auto g_xxyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 320);

    auto g_xxyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 321);

    auto g_xxyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 322);

    auto g_xxyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 323);

    auto g_xxyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 324);

    auto g_xxyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 325);

    auto g_xxyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 326);

    auto g_xxyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 327);

    auto g_xxyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 328);

    auto g_xxyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 329);

    auto g_xxyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 330);

    auto g_xxyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 331);

    auto g_xxyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 332);

    auto g_xxyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 333);

    auto g_xxyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 334);

    auto g_xxyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 335);

    #pragma omp simd aligned(g_xxyyy_0_xxxxx_1, g_xxyyy_0_xxxxxx_1, g_xxyyy_0_xxxxxy_1, g_xxyyy_0_xxxxxz_1, g_xxyyy_0_xxxxy_1, g_xxyyy_0_xxxxyy_1, g_xxyyy_0_xxxxyz_1, g_xxyyy_0_xxxxz_1, g_xxyyy_0_xxxxzz_1, g_xxyyy_0_xxxyy_1, g_xxyyy_0_xxxyyy_1, g_xxyyy_0_xxxyyz_1, g_xxyyy_0_xxxyz_1, g_xxyyy_0_xxxyzz_1, g_xxyyy_0_xxxzz_1, g_xxyyy_0_xxxzzz_1, g_xxyyy_0_xxyyy_1, g_xxyyy_0_xxyyyy_1, g_xxyyy_0_xxyyyz_1, g_xxyyy_0_xxyyz_1, g_xxyyy_0_xxyyzz_1, g_xxyyy_0_xxyzz_1, g_xxyyy_0_xxyzzz_1, g_xxyyy_0_xxzzz_1, g_xxyyy_0_xxzzzz_1, g_xxyyy_0_xyyyy_1, g_xxyyy_0_xyyyyy_1, g_xxyyy_0_xyyyyz_1, g_xxyyy_0_xyyyz_1, g_xxyyy_0_xyyyzz_1, g_xxyyy_0_xyyzz_1, g_xxyyy_0_xyyzzz_1, g_xxyyy_0_xyzzz_1, g_xxyyy_0_xyzzzz_1, g_xxyyy_0_xzzzz_1, g_xxyyy_0_xzzzzz_1, g_xxyyy_0_yyyyy_1, g_xxyyy_0_yyyyyy_1, g_xxyyy_0_yyyyyz_1, g_xxyyy_0_yyyyz_1, g_xxyyy_0_yyyyzz_1, g_xxyyy_0_yyyzz_1, g_xxyyy_0_yyyzzz_1, g_xxyyy_0_yyzzz_1, g_xxyyy_0_yyzzzz_1, g_xxyyy_0_yzzzz_1, g_xxyyy_0_yzzzzz_1, g_xxyyy_0_zzzzz_1, g_xxyyy_0_zzzzzz_1, g_xxyyyz_0_xxxxxx_0, g_xxyyyz_0_xxxxxy_0, g_xxyyyz_0_xxxxxz_0, g_xxyyyz_0_xxxxyy_0, g_xxyyyz_0_xxxxyz_0, g_xxyyyz_0_xxxxzz_0, g_xxyyyz_0_xxxyyy_0, g_xxyyyz_0_xxxyyz_0, g_xxyyyz_0_xxxyzz_0, g_xxyyyz_0_xxxzzz_0, g_xxyyyz_0_xxyyyy_0, g_xxyyyz_0_xxyyyz_0, g_xxyyyz_0_xxyyzz_0, g_xxyyyz_0_xxyzzz_0, g_xxyyyz_0_xxzzzz_0, g_xxyyyz_0_xyyyyy_0, g_xxyyyz_0_xyyyyz_0, g_xxyyyz_0_xyyyzz_0, g_xxyyyz_0_xyyzzz_0, g_xxyyyz_0_xyzzzz_0, g_xxyyyz_0_xzzzzz_0, g_xxyyyz_0_yyyyyy_0, g_xxyyyz_0_yyyyyz_0, g_xxyyyz_0_yyyyzz_0, g_xxyyyz_0_yyyzzz_0, g_xxyyyz_0_yyzzzz_0, g_xxyyyz_0_yzzzzz_0, g_xxyyyz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyz_0_xxxxxx_0[i] = g_xxyyy_0_xxxxxx_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxy_0[i] = g_xxyyy_0_xxxxxy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxz_0[i] = g_xxyyy_0_xxxxx_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxyy_0[i] = g_xxyyy_0_xxxxyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxyz_0[i] = g_xxyyy_0_xxxxy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxzz_0[i] = 2.0 * g_xxyyy_0_xxxxz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyyy_0[i] = g_xxyyy_0_xxxyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyyz_0[i] = g_xxyyy_0_xxxyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyzz_0[i] = 2.0 * g_xxyyy_0_xxxyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxzzz_0[i] = 3.0 * g_xxyyy_0_xxxzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyyy_0[i] = g_xxyyy_0_xxyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyyz_0[i] = g_xxyyy_0_xxyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyzz_0[i] = 2.0 * g_xxyyy_0_xxyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyzzz_0[i] = 3.0 * g_xxyyy_0_xxyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxzzzz_0[i] = 4.0 * g_xxyyy_0_xxzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyyy_0[i] = g_xxyyy_0_xyyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyyz_0[i] = g_xxyyy_0_xyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyzz_0[i] = 2.0 * g_xxyyy_0_xyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyzzz_0[i] = 3.0 * g_xxyyy_0_xyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyzzzz_0[i] = 4.0 * g_xxyyy_0_xyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xzzzzz_0[i] = 5.0 * g_xxyyy_0_xzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyyy_0[i] = g_xxyyy_0_yyyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyyz_0[i] = g_xxyyy_0_yyyyy_1[i] * fi_acd_0 + g_xxyyy_0_yyyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyzz_0[i] = 2.0 * g_xxyyy_0_yyyyz_1[i] * fi_acd_0 + g_xxyyy_0_yyyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyzzz_0[i] = 3.0 * g_xxyyy_0_yyyzz_1[i] * fi_acd_0 + g_xxyyy_0_yyyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyzzzz_0[i] = 4.0 * g_xxyyy_0_yyzzz_1[i] * fi_acd_0 + g_xxyyy_0_yyzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yzzzzz_0[i] = 5.0 * g_xxyyy_0_yzzzz_1[i] * fi_acd_0 + g_xxyyy_0_yzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_zzzzzz_0[i] = 6.0 * g_xxyyy_0_zzzzz_1[i] * fi_acd_0 + g_xxyyy_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 336-364 components of targeted buffer : ISI

    auto g_xxyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 336);

    auto g_xxyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 337);

    auto g_xxyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 338);

    auto g_xxyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 339);

    auto g_xxyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 340);

    auto g_xxyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 341);

    auto g_xxyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 342);

    auto g_xxyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 343);

    auto g_xxyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 344);

    auto g_xxyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 345);

    auto g_xxyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 346);

    auto g_xxyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 347);

    auto g_xxyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 348);

    auto g_xxyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 349);

    auto g_xxyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 350);

    auto g_xxyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 351);

    auto g_xxyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 352);

    auto g_xxyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 353);

    auto g_xxyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 354);

    auto g_xxyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 355);

    auto g_xxyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 356);

    auto g_xxyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 357);

    auto g_xxyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 358);

    auto g_xxyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 359);

    auto g_xxyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 360);

    auto g_xxyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 361);

    auto g_xxyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 362);

    auto g_xxyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 363);

    #pragma omp simd aligned(g_xxyy_0_xxxxxy_0, g_xxyy_0_xxxxxy_1, g_xxyy_0_xxxxyy_0, g_xxyy_0_xxxxyy_1, g_xxyy_0_xxxyyy_0, g_xxyy_0_xxxyyy_1, g_xxyy_0_xxyyyy_0, g_xxyy_0_xxyyyy_1, g_xxyy_0_xyyyyy_0, g_xxyy_0_xyyyyy_1, g_xxyyz_0_xxxxxy_1, g_xxyyz_0_xxxxyy_1, g_xxyyz_0_xxxyyy_1, g_xxyyz_0_xxyyyy_1, g_xxyyz_0_xyyyyy_1, g_xxyyzz_0_xxxxxx_0, g_xxyyzz_0_xxxxxy_0, g_xxyyzz_0_xxxxxz_0, g_xxyyzz_0_xxxxyy_0, g_xxyyzz_0_xxxxyz_0, g_xxyyzz_0_xxxxzz_0, g_xxyyzz_0_xxxyyy_0, g_xxyyzz_0_xxxyyz_0, g_xxyyzz_0_xxxyzz_0, g_xxyyzz_0_xxxzzz_0, g_xxyyzz_0_xxyyyy_0, g_xxyyzz_0_xxyyyz_0, g_xxyyzz_0_xxyyzz_0, g_xxyyzz_0_xxyzzz_0, g_xxyyzz_0_xxzzzz_0, g_xxyyzz_0_xyyyyy_0, g_xxyyzz_0_xyyyyz_0, g_xxyyzz_0_xyyyzz_0, g_xxyyzz_0_xyyzzz_0, g_xxyyzz_0_xyzzzz_0, g_xxyyzz_0_xzzzzz_0, g_xxyyzz_0_yyyyyy_0, g_xxyyzz_0_yyyyyz_0, g_xxyyzz_0_yyyyzz_0, g_xxyyzz_0_yyyzzz_0, g_xxyyzz_0_yyzzzz_0, g_xxyyzz_0_yzzzzz_0, g_xxyyzz_0_zzzzzz_0, g_xxyzz_0_xxxxxx_1, g_xxyzz_0_xxxxxz_1, g_xxyzz_0_xxxxzz_1, g_xxyzz_0_xxxzzz_1, g_xxyzz_0_xxzzzz_1, g_xxyzz_0_xzzzzz_1, g_xxzz_0_xxxxxx_0, g_xxzz_0_xxxxxx_1, g_xxzz_0_xxxxxz_0, g_xxzz_0_xxxxxz_1, g_xxzz_0_xxxxzz_0, g_xxzz_0_xxxxzz_1, g_xxzz_0_xxxzzz_0, g_xxzz_0_xxxzzz_1, g_xxzz_0_xxzzzz_0, g_xxzz_0_xxzzzz_1, g_xxzz_0_xzzzzz_0, g_xxzz_0_xzzzzz_1, g_xyyzz_0_xxxxyz_1, g_xyyzz_0_xxxyyz_1, g_xyyzz_0_xxxyz_1, g_xyyzz_0_xxxyzz_1, g_xyyzz_0_xxyyyz_1, g_xyyzz_0_xxyyz_1, g_xyyzz_0_xxyyzz_1, g_xyyzz_0_xxyzz_1, g_xyyzz_0_xxyzzz_1, g_xyyzz_0_xyyyyz_1, g_xyyzz_0_xyyyz_1, g_xyyzz_0_xyyyzz_1, g_xyyzz_0_xyyzz_1, g_xyyzz_0_xyyzzz_1, g_xyyzz_0_xyzzz_1, g_xyyzz_0_xyzzzz_1, g_xyyzz_0_yyyyyy_1, g_xyyzz_0_yyyyyz_1, g_xyyzz_0_yyyyz_1, g_xyyzz_0_yyyyzz_1, g_xyyzz_0_yyyzz_1, g_xyyzz_0_yyyzzz_1, g_xyyzz_0_yyzzz_1, g_xyyzz_0_yyzzzz_1, g_xyyzz_0_yzzzz_1, g_xyyzz_0_yzzzzz_1, g_xyyzz_0_zzzzzz_1, g_yyzz_0_xxxxyz_0, g_yyzz_0_xxxxyz_1, g_yyzz_0_xxxyyz_0, g_yyzz_0_xxxyyz_1, g_yyzz_0_xxxyzz_0, g_yyzz_0_xxxyzz_1, g_yyzz_0_xxyyyz_0, g_yyzz_0_xxyyyz_1, g_yyzz_0_xxyyzz_0, g_yyzz_0_xxyyzz_1, g_yyzz_0_xxyzzz_0, g_yyzz_0_xxyzzz_1, g_yyzz_0_xyyyyz_0, g_yyzz_0_xyyyyz_1, g_yyzz_0_xyyyzz_0, g_yyzz_0_xyyyzz_1, g_yyzz_0_xyyzzz_0, g_yyzz_0_xyyzzz_1, g_yyzz_0_xyzzzz_0, g_yyzz_0_xyzzzz_1, g_yyzz_0_yyyyyy_0, g_yyzz_0_yyyyyy_1, g_yyzz_0_yyyyyz_0, g_yyzz_0_yyyyyz_1, g_yyzz_0_yyyyzz_0, g_yyzz_0_yyyyzz_1, g_yyzz_0_yyyzzz_0, g_yyzz_0_yyyzzz_1, g_yyzz_0_yyzzzz_0, g_yyzz_0_yyzzzz_1, g_yyzz_0_yzzzzz_0, g_yyzz_0_yzzzzz_1, g_yyzz_0_zzzzzz_0, g_yyzz_0_zzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyzz_0_xxxxxx_0[i] = g_xxzz_0_xxxxxx_0[i] * fbe_0 - g_xxzz_0_xxxxxx_1[i] * fz_be_0 + g_xxyzz_0_xxxxxx_1[i] * wa_y[i];

        g_xxyyzz_0_xxxxxy_0[i] = g_xxyy_0_xxxxxy_0[i] * fbe_0 - g_xxyy_0_xxxxxy_1[i] * fz_be_0 + g_xxyyz_0_xxxxxy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxxxz_0[i] = g_xxzz_0_xxxxxz_0[i] * fbe_0 - g_xxzz_0_xxxxxz_1[i] * fz_be_0 + g_xxyzz_0_xxxxxz_1[i] * wa_y[i];

        g_xxyyzz_0_xxxxyy_0[i] = g_xxyy_0_xxxxyy_0[i] * fbe_0 - g_xxyy_0_xxxxyy_1[i] * fz_be_0 + g_xxyyz_0_xxxxyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxxyz_0[i] = g_yyzz_0_xxxxyz_0[i] * fbe_0 - g_yyzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyyzz_0_xxxyz_1[i] * fi_acd_0 + g_xyyzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxxzz_0[i] = g_xxzz_0_xxxxzz_0[i] * fbe_0 - g_xxzz_0_xxxxzz_1[i] * fz_be_0 + g_xxyzz_0_xxxxzz_1[i] * wa_y[i];

        g_xxyyzz_0_xxxyyy_0[i] = g_xxyy_0_xxxyyy_0[i] * fbe_0 - g_xxyy_0_xxxyyy_1[i] * fz_be_0 + g_xxyyz_0_xxxyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxyyz_0[i] = g_yyzz_0_xxxyyz_0[i] * fbe_0 - g_yyzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyyzz_0_xxyyz_1[i] * fi_acd_0 + g_xyyzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxyzz_0[i] = g_yyzz_0_xxxyzz_0[i] * fbe_0 - g_yyzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyyzz_0_xxyzz_1[i] * fi_acd_0 + g_xyyzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxzzz_0[i] = g_xxzz_0_xxxzzz_0[i] * fbe_0 - g_xxzz_0_xxxzzz_1[i] * fz_be_0 + g_xxyzz_0_xxxzzz_1[i] * wa_y[i];

        g_xxyyzz_0_xxyyyy_0[i] = g_xxyy_0_xxyyyy_0[i] * fbe_0 - g_xxyy_0_xxyyyy_1[i] * fz_be_0 + g_xxyyz_0_xxyyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxyyyz_0[i] = g_yyzz_0_xxyyyz_0[i] * fbe_0 - g_yyzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyyyz_1[i] * fi_acd_0 + g_xyyzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxyyzz_0[i] = g_yyzz_0_xxyyzz_0[i] * fbe_0 - g_yyzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyyzz_1[i] * fi_acd_0 + g_xyyzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxyzzz_0[i] = g_yyzz_0_xxyzzz_0[i] * fbe_0 - g_yyzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyzzz_1[i] * fi_acd_0 + g_xyyzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxzzzz_0[i] = g_xxzz_0_xxzzzz_0[i] * fbe_0 - g_xxzz_0_xxzzzz_1[i] * fz_be_0 + g_xxyzz_0_xxzzzz_1[i] * wa_y[i];

        g_xxyyzz_0_xyyyyy_0[i] = g_xxyy_0_xyyyyy_0[i] * fbe_0 - g_xxyy_0_xyyyyy_1[i] * fz_be_0 + g_xxyyz_0_xyyyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xyyyyz_0[i] = g_yyzz_0_xyyyyz_0[i] * fbe_0 - g_yyzz_0_xyyyyz_1[i] * fz_be_0 + g_xyyzz_0_yyyyz_1[i] * fi_acd_0 + g_xyyzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xyyyzz_0[i] = g_yyzz_0_xyyyzz_0[i] * fbe_0 - g_yyzz_0_xyyyzz_1[i] * fz_be_0 + g_xyyzz_0_yyyzz_1[i] * fi_acd_0 + g_xyyzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xyyzzz_0[i] = g_yyzz_0_xyyzzz_0[i] * fbe_0 - g_yyzz_0_xyyzzz_1[i] * fz_be_0 + g_xyyzz_0_yyzzz_1[i] * fi_acd_0 + g_xyyzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xyzzzz_0[i] = g_yyzz_0_xyzzzz_0[i] * fbe_0 - g_yyzz_0_xyzzzz_1[i] * fz_be_0 + g_xyyzz_0_yzzzz_1[i] * fi_acd_0 + g_xyyzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xzzzzz_0[i] = g_xxzz_0_xzzzzz_0[i] * fbe_0 - g_xxzz_0_xzzzzz_1[i] * fz_be_0 + g_xxyzz_0_xzzzzz_1[i] * wa_y[i];

        g_xxyyzz_0_yyyyyy_0[i] = g_yyzz_0_yyyyyy_0[i] * fbe_0 - g_yyzz_0_yyyyyy_1[i] * fz_be_0 + g_xyyzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxyyzz_0_yyyyyz_0[i] = g_yyzz_0_yyyyyz_0[i] * fbe_0 - g_yyzz_0_yyyyyz_1[i] * fz_be_0 + g_xyyzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_yyyyzz_0[i] = g_yyzz_0_yyyyzz_0[i] * fbe_0 - g_yyzz_0_yyyyzz_1[i] * fz_be_0 + g_xyyzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_yyyzzz_0[i] = g_yyzz_0_yyyzzz_0[i] * fbe_0 - g_yyzz_0_yyyzzz_1[i] * fz_be_0 + g_xyyzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_yyzzzz_0[i] = g_yyzz_0_yyzzzz_0[i] * fbe_0 - g_yyzz_0_yyzzzz_1[i] * fz_be_0 + g_xyyzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_yzzzzz_0[i] = g_yyzz_0_yzzzzz_0[i] * fbe_0 - g_yyzz_0_yzzzzz_1[i] * fz_be_0 + g_xyyzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_zzzzzz_0[i] = g_yyzz_0_zzzzzz_0[i] * fbe_0 - g_yyzz_0_zzzzzz_1[i] * fz_be_0 + g_xyyzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 364-392 components of targeted buffer : ISI

    auto g_xxyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 364);

    auto g_xxyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 365);

    auto g_xxyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 366);

    auto g_xxyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 367);

    auto g_xxyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 368);

    auto g_xxyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 369);

    auto g_xxyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 370);

    auto g_xxyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 371);

    auto g_xxyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 372);

    auto g_xxyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 373);

    auto g_xxyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 374);

    auto g_xxyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 375);

    auto g_xxyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 376);

    auto g_xxyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 377);

    auto g_xxyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 378);

    auto g_xxyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 379);

    auto g_xxyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 380);

    auto g_xxyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 381);

    auto g_xxyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 382);

    auto g_xxyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 383);

    auto g_xxyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 384);

    auto g_xxyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 385);

    auto g_xxyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 386);

    auto g_xxyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 387);

    auto g_xxyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 388);

    auto g_xxyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 389);

    auto g_xxyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 390);

    auto g_xxyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 391);

    #pragma omp simd aligned(g_xxyzzz_0_xxxxxx_0, g_xxyzzz_0_xxxxxy_0, g_xxyzzz_0_xxxxxz_0, g_xxyzzz_0_xxxxyy_0, g_xxyzzz_0_xxxxyz_0, g_xxyzzz_0_xxxxzz_0, g_xxyzzz_0_xxxyyy_0, g_xxyzzz_0_xxxyyz_0, g_xxyzzz_0_xxxyzz_0, g_xxyzzz_0_xxxzzz_0, g_xxyzzz_0_xxyyyy_0, g_xxyzzz_0_xxyyyz_0, g_xxyzzz_0_xxyyzz_0, g_xxyzzz_0_xxyzzz_0, g_xxyzzz_0_xxzzzz_0, g_xxyzzz_0_xyyyyy_0, g_xxyzzz_0_xyyyyz_0, g_xxyzzz_0_xyyyzz_0, g_xxyzzz_0_xyyzzz_0, g_xxyzzz_0_xyzzzz_0, g_xxyzzz_0_xzzzzz_0, g_xxyzzz_0_yyyyyy_0, g_xxyzzz_0_yyyyyz_0, g_xxyzzz_0_yyyyzz_0, g_xxyzzz_0_yyyzzz_0, g_xxyzzz_0_yyzzzz_0, g_xxyzzz_0_yzzzzz_0, g_xxyzzz_0_zzzzzz_0, g_xxzzz_0_xxxxx_1, g_xxzzz_0_xxxxxx_1, g_xxzzz_0_xxxxxy_1, g_xxzzz_0_xxxxxz_1, g_xxzzz_0_xxxxy_1, g_xxzzz_0_xxxxyy_1, g_xxzzz_0_xxxxyz_1, g_xxzzz_0_xxxxz_1, g_xxzzz_0_xxxxzz_1, g_xxzzz_0_xxxyy_1, g_xxzzz_0_xxxyyy_1, g_xxzzz_0_xxxyyz_1, g_xxzzz_0_xxxyz_1, g_xxzzz_0_xxxyzz_1, g_xxzzz_0_xxxzz_1, g_xxzzz_0_xxxzzz_1, g_xxzzz_0_xxyyy_1, g_xxzzz_0_xxyyyy_1, g_xxzzz_0_xxyyyz_1, g_xxzzz_0_xxyyz_1, g_xxzzz_0_xxyyzz_1, g_xxzzz_0_xxyzz_1, g_xxzzz_0_xxyzzz_1, g_xxzzz_0_xxzzz_1, g_xxzzz_0_xxzzzz_1, g_xxzzz_0_xyyyy_1, g_xxzzz_0_xyyyyy_1, g_xxzzz_0_xyyyyz_1, g_xxzzz_0_xyyyz_1, g_xxzzz_0_xyyyzz_1, g_xxzzz_0_xyyzz_1, g_xxzzz_0_xyyzzz_1, g_xxzzz_0_xyzzz_1, g_xxzzz_0_xyzzzz_1, g_xxzzz_0_xzzzz_1, g_xxzzz_0_xzzzzz_1, g_xxzzz_0_yyyyy_1, g_xxzzz_0_yyyyyy_1, g_xxzzz_0_yyyyyz_1, g_xxzzz_0_yyyyz_1, g_xxzzz_0_yyyyzz_1, g_xxzzz_0_yyyzz_1, g_xxzzz_0_yyyzzz_1, g_xxzzz_0_yyzzz_1, g_xxzzz_0_yyzzzz_1, g_xxzzz_0_yzzzz_1, g_xxzzz_0_yzzzzz_1, g_xxzzz_0_zzzzz_1, g_xxzzz_0_zzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzz_0_xxxxxx_0[i] = g_xxzzz_0_xxxxxx_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxy_0[i] = g_xxzzz_0_xxxxx_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxz_0[i] = g_xxzzz_0_xxxxxz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxyy_0[i] = 2.0 * g_xxzzz_0_xxxxy_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxyz_0[i] = g_xxzzz_0_xxxxz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxzz_0[i] = g_xxzzz_0_xxxxzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyyy_0[i] = 3.0 * g_xxzzz_0_xxxyy_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyyz_0[i] = 2.0 * g_xxzzz_0_xxxyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyzz_0[i] = g_xxzzz_0_xxxzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxzzz_0[i] = g_xxzzz_0_xxxzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyyy_0[i] = 4.0 * g_xxzzz_0_xxyyy_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyyz_0[i] = 3.0 * g_xxzzz_0_xxyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyzz_0[i] = 2.0 * g_xxzzz_0_xxyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyzzz_0[i] = g_xxzzz_0_xxzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxzzzz_0[i] = g_xxzzz_0_xxzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyyy_0[i] = 5.0 * g_xxzzz_0_xyyyy_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyyz_0[i] = 4.0 * g_xxzzz_0_xyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyzz_0[i] = 3.0 * g_xxzzz_0_xyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyzzz_0[i] = 2.0 * g_xxzzz_0_xyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyzzzz_0[i] = g_xxzzz_0_xzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xzzzzz_0[i] = g_xxzzz_0_xzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyyy_0[i] = 6.0 * g_xxzzz_0_yyyyy_1[i] * fi_acd_0 + g_xxzzz_0_yyyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyyz_0[i] = 5.0 * g_xxzzz_0_yyyyz_1[i] * fi_acd_0 + g_xxzzz_0_yyyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyzz_0[i] = 4.0 * g_xxzzz_0_yyyzz_1[i] * fi_acd_0 + g_xxzzz_0_yyyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyzzz_0[i] = 3.0 * g_xxzzz_0_yyzzz_1[i] * fi_acd_0 + g_xxzzz_0_yyyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyzzzz_0[i] = 2.0 * g_xxzzz_0_yzzzz_1[i] * fi_acd_0 + g_xxzzz_0_yyzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yzzzzz_0[i] = g_xxzzz_0_zzzzz_1[i] * fi_acd_0 + g_xxzzz_0_yzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_zzzzzz_0[i] = g_xxzzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 392-420 components of targeted buffer : ISI

    auto g_xxzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 392);

    auto g_xxzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 393);

    auto g_xxzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 394);

    auto g_xxzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 395);

    auto g_xxzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 396);

    auto g_xxzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 397);

    auto g_xxzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 398);

    auto g_xxzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 399);

    auto g_xxzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 400);

    auto g_xxzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 401);

    auto g_xxzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 402);

    auto g_xxzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 403);

    auto g_xxzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 404);

    auto g_xxzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 405);

    auto g_xxzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 406);

    auto g_xxzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 407);

    auto g_xxzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 408);

    auto g_xxzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 409);

    auto g_xxzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 410);

    auto g_xxzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 411);

    auto g_xxzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 412);

    auto g_xxzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 413);

    auto g_xxzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 414);

    auto g_xxzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 415);

    auto g_xxzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 416);

    auto g_xxzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 417);

    auto g_xxzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 418);

    auto g_xxzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 419);

    #pragma omp simd aligned(g_xxzz_0_xxxxxx_0, g_xxzz_0_xxxxxx_1, g_xxzz_0_xxxxxy_0, g_xxzz_0_xxxxxy_1, g_xxzz_0_xxxxyy_0, g_xxzz_0_xxxxyy_1, g_xxzz_0_xxxyyy_0, g_xxzz_0_xxxyyy_1, g_xxzz_0_xxyyyy_0, g_xxzz_0_xxyyyy_1, g_xxzz_0_xyyyyy_0, g_xxzz_0_xyyyyy_1, g_xxzzz_0_xxxxxx_1, g_xxzzz_0_xxxxxy_1, g_xxzzz_0_xxxxyy_1, g_xxzzz_0_xxxyyy_1, g_xxzzz_0_xxyyyy_1, g_xxzzz_0_xyyyyy_1, g_xxzzzz_0_xxxxxx_0, g_xxzzzz_0_xxxxxy_0, g_xxzzzz_0_xxxxxz_0, g_xxzzzz_0_xxxxyy_0, g_xxzzzz_0_xxxxyz_0, g_xxzzzz_0_xxxxzz_0, g_xxzzzz_0_xxxyyy_0, g_xxzzzz_0_xxxyyz_0, g_xxzzzz_0_xxxyzz_0, g_xxzzzz_0_xxxzzz_0, g_xxzzzz_0_xxyyyy_0, g_xxzzzz_0_xxyyyz_0, g_xxzzzz_0_xxyyzz_0, g_xxzzzz_0_xxyzzz_0, g_xxzzzz_0_xxzzzz_0, g_xxzzzz_0_xyyyyy_0, g_xxzzzz_0_xyyyyz_0, g_xxzzzz_0_xyyyzz_0, g_xxzzzz_0_xyyzzz_0, g_xxzzzz_0_xyzzzz_0, g_xxzzzz_0_xzzzzz_0, g_xxzzzz_0_yyyyyy_0, g_xxzzzz_0_yyyyyz_0, g_xxzzzz_0_yyyyzz_0, g_xxzzzz_0_yyyzzz_0, g_xxzzzz_0_yyzzzz_0, g_xxzzzz_0_yzzzzz_0, g_xxzzzz_0_zzzzzz_0, g_xzzzz_0_xxxxxz_1, g_xzzzz_0_xxxxyz_1, g_xzzzz_0_xxxxz_1, g_xzzzz_0_xxxxzz_1, g_xzzzz_0_xxxyyz_1, g_xzzzz_0_xxxyz_1, g_xzzzz_0_xxxyzz_1, g_xzzzz_0_xxxzz_1, g_xzzzz_0_xxxzzz_1, g_xzzzz_0_xxyyyz_1, g_xzzzz_0_xxyyz_1, g_xzzzz_0_xxyyzz_1, g_xzzzz_0_xxyzz_1, g_xzzzz_0_xxyzzz_1, g_xzzzz_0_xxzzz_1, g_xzzzz_0_xxzzzz_1, g_xzzzz_0_xyyyyz_1, g_xzzzz_0_xyyyz_1, g_xzzzz_0_xyyyzz_1, g_xzzzz_0_xyyzz_1, g_xzzzz_0_xyyzzz_1, g_xzzzz_0_xyzzz_1, g_xzzzz_0_xyzzzz_1, g_xzzzz_0_xzzzz_1, g_xzzzz_0_xzzzzz_1, g_xzzzz_0_yyyyyy_1, g_xzzzz_0_yyyyyz_1, g_xzzzz_0_yyyyz_1, g_xzzzz_0_yyyyzz_1, g_xzzzz_0_yyyzz_1, g_xzzzz_0_yyyzzz_1, g_xzzzz_0_yyzzz_1, g_xzzzz_0_yyzzzz_1, g_xzzzz_0_yzzzz_1, g_xzzzz_0_yzzzzz_1, g_xzzzz_0_zzzzz_1, g_xzzzz_0_zzzzzz_1, g_zzzz_0_xxxxxz_0, g_zzzz_0_xxxxxz_1, g_zzzz_0_xxxxyz_0, g_zzzz_0_xxxxyz_1, g_zzzz_0_xxxxzz_0, g_zzzz_0_xxxxzz_1, g_zzzz_0_xxxyyz_0, g_zzzz_0_xxxyyz_1, g_zzzz_0_xxxyzz_0, g_zzzz_0_xxxyzz_1, g_zzzz_0_xxxzzz_0, g_zzzz_0_xxxzzz_1, g_zzzz_0_xxyyyz_0, g_zzzz_0_xxyyyz_1, g_zzzz_0_xxyyzz_0, g_zzzz_0_xxyyzz_1, g_zzzz_0_xxyzzz_0, g_zzzz_0_xxyzzz_1, g_zzzz_0_xxzzzz_0, g_zzzz_0_xxzzzz_1, g_zzzz_0_xyyyyz_0, g_zzzz_0_xyyyyz_1, g_zzzz_0_xyyyzz_0, g_zzzz_0_xyyyzz_1, g_zzzz_0_xyyzzz_0, g_zzzz_0_xyyzzz_1, g_zzzz_0_xyzzzz_0, g_zzzz_0_xyzzzz_1, g_zzzz_0_xzzzzz_0, g_zzzz_0_xzzzzz_1, g_zzzz_0_yyyyyy_0, g_zzzz_0_yyyyyy_1, g_zzzz_0_yyyyyz_0, g_zzzz_0_yyyyyz_1, g_zzzz_0_yyyyzz_0, g_zzzz_0_yyyyzz_1, g_zzzz_0_yyyzzz_0, g_zzzz_0_yyyzzz_1, g_zzzz_0_yyzzzz_0, g_zzzz_0_yyzzzz_1, g_zzzz_0_yzzzzz_0, g_zzzz_0_yzzzzz_1, g_zzzz_0_zzzzzz_0, g_zzzz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzz_0_xxxxxx_0[i] = 3.0 * g_xxzz_0_xxxxxx_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxx_1[i] * fz_be_0 + g_xxzzz_0_xxxxxx_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxxy_0[i] = 3.0 * g_xxzz_0_xxxxxy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxy_1[i] * fz_be_0 + g_xxzzz_0_xxxxxy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxxz_0[i] = g_zzzz_0_xxxxxz_0[i] * fbe_0 - g_zzzz_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xzzzz_0_xxxxz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxxz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxyy_0[i] = 3.0 * g_xxzz_0_xxxxyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxyy_1[i] * fz_be_0 + g_xxzzz_0_xxxxyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxyz_0[i] = g_zzzz_0_xxxxyz_0[i] * fbe_0 - g_zzzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xzzzz_0_xxxyz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxzz_0[i] = g_zzzz_0_xxxxzz_0[i] * fbe_0 - g_zzzz_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xzzzz_0_xxxzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxyyy_0[i] = 3.0 * g_xxzz_0_xxxyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyyy_1[i] * fz_be_0 + g_xxzzz_0_xxxyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxyyz_0[i] = g_zzzz_0_xxxyyz_0[i] * fbe_0 - g_zzzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxyyz_1[i] * fi_acd_0 + g_xzzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxyzz_0[i] = g_zzzz_0_xxxyzz_0[i] * fbe_0 - g_zzzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxyzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxzzz_0[i] = g_zzzz_0_xxxzzz_0[i] * fbe_0 - g_zzzz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyyyy_0[i] = 3.0 * g_xxzz_0_xxyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyyy_1[i] * fz_be_0 + g_xxzzz_0_xxyyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxyyyz_0[i] = g_zzzz_0_xxyyyz_0[i] * fbe_0 - g_zzzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyyyz_1[i] * fi_acd_0 + g_xzzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyyzz_0[i] = g_zzzz_0_xxyyzz_0[i] * fbe_0 - g_zzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyyzz_1[i] * fi_acd_0 + g_xzzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyzzz_0[i] = g_zzzz_0_xxyzzz_0[i] * fbe_0 - g_zzzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxzzzz_0[i] = g_zzzz_0_xxzzzz_0[i] * fbe_0 - g_zzzz_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyyyy_0[i] = 3.0 * g_xxzz_0_xyyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyyy_1[i] * fz_be_0 + g_xxzzz_0_xyyyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xyyyyz_0[i] = g_zzzz_0_xyyyyz_0[i] * fbe_0 - g_zzzz_0_xyyyyz_1[i] * fz_be_0 + g_xzzzz_0_yyyyz_1[i] * fi_acd_0 + g_xzzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyyzz_0[i] = g_zzzz_0_xyyyzz_0[i] * fbe_0 - g_zzzz_0_xyyyzz_1[i] * fz_be_0 + g_xzzzz_0_yyyzz_1[i] * fi_acd_0 + g_xzzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyzzz_0[i] = g_zzzz_0_xyyzzz_0[i] * fbe_0 - g_zzzz_0_xyyzzz_1[i] * fz_be_0 + g_xzzzz_0_yyzzz_1[i] * fi_acd_0 + g_xzzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyzzzz_0[i] = g_zzzz_0_xyzzzz_0[i] * fbe_0 - g_zzzz_0_xyzzzz_1[i] * fz_be_0 + g_xzzzz_0_yzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xzzzzz_0[i] = g_zzzz_0_xzzzzz_0[i] * fbe_0 - g_zzzz_0_xzzzzz_1[i] * fz_be_0 + g_xzzzz_0_zzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyyy_0[i] = g_zzzz_0_yyyyyy_0[i] * fbe_0 - g_zzzz_0_yyyyyy_1[i] * fz_be_0 + g_xzzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyyz_0[i] = g_zzzz_0_yyyyyz_0[i] * fbe_0 - g_zzzz_0_yyyyyz_1[i] * fz_be_0 + g_xzzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyzz_0[i] = g_zzzz_0_yyyyzz_0[i] * fbe_0 - g_zzzz_0_yyyyzz_1[i] * fz_be_0 + g_xzzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyzzz_0[i] = g_zzzz_0_yyyzzz_0[i] * fbe_0 - g_zzzz_0_yyyzzz_1[i] * fz_be_0 + g_xzzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyzzzz_0[i] = g_zzzz_0_yyzzzz_0[i] * fbe_0 - g_zzzz_0_yyzzzz_1[i] * fz_be_0 + g_xzzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yzzzzz_0[i] = g_zzzz_0_yzzzzz_0[i] * fbe_0 - g_zzzz_0_yzzzzz_1[i] * fz_be_0 + g_xzzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_zzzzzz_0[i] = g_zzzz_0_zzzzzz_0[i] * fbe_0 - g_zzzz_0_zzzzzz_1[i] * fz_be_0 + g_xzzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 420-448 components of targeted buffer : ISI

    auto g_xyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 420);

    auto g_xyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 421);

    auto g_xyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 422);

    auto g_xyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 423);

    auto g_xyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 424);

    auto g_xyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 425);

    auto g_xyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 426);

    auto g_xyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 427);

    auto g_xyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 428);

    auto g_xyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 429);

    auto g_xyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 430);

    auto g_xyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 431);

    auto g_xyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 432);

    auto g_xyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 433);

    auto g_xyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 434);

    auto g_xyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 435);

    auto g_xyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 436);

    auto g_xyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 437);

    auto g_xyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 438);

    auto g_xyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 439);

    auto g_xyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 440);

    auto g_xyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 441);

    auto g_xyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 442);

    auto g_xyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 443);

    auto g_xyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 444);

    auto g_xyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 445);

    auto g_xyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 446);

    auto g_xyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 447);

    #pragma omp simd aligned(g_xyyyyy_0_xxxxxx_0, g_xyyyyy_0_xxxxxy_0, g_xyyyyy_0_xxxxxz_0, g_xyyyyy_0_xxxxyy_0, g_xyyyyy_0_xxxxyz_0, g_xyyyyy_0_xxxxzz_0, g_xyyyyy_0_xxxyyy_0, g_xyyyyy_0_xxxyyz_0, g_xyyyyy_0_xxxyzz_0, g_xyyyyy_0_xxxzzz_0, g_xyyyyy_0_xxyyyy_0, g_xyyyyy_0_xxyyyz_0, g_xyyyyy_0_xxyyzz_0, g_xyyyyy_0_xxyzzz_0, g_xyyyyy_0_xxzzzz_0, g_xyyyyy_0_xyyyyy_0, g_xyyyyy_0_xyyyyz_0, g_xyyyyy_0_xyyyzz_0, g_xyyyyy_0_xyyzzz_0, g_xyyyyy_0_xyzzzz_0, g_xyyyyy_0_xzzzzz_0, g_xyyyyy_0_yyyyyy_0, g_xyyyyy_0_yyyyyz_0, g_xyyyyy_0_yyyyzz_0, g_xyyyyy_0_yyyzzz_0, g_xyyyyy_0_yyzzzz_0, g_xyyyyy_0_yzzzzz_0, g_xyyyyy_0_zzzzzz_0, g_yyyyy_0_xxxxx_1, g_yyyyy_0_xxxxxx_1, g_yyyyy_0_xxxxxy_1, g_yyyyy_0_xxxxxz_1, g_yyyyy_0_xxxxy_1, g_yyyyy_0_xxxxyy_1, g_yyyyy_0_xxxxyz_1, g_yyyyy_0_xxxxz_1, g_yyyyy_0_xxxxzz_1, g_yyyyy_0_xxxyy_1, g_yyyyy_0_xxxyyy_1, g_yyyyy_0_xxxyyz_1, g_yyyyy_0_xxxyz_1, g_yyyyy_0_xxxyzz_1, g_yyyyy_0_xxxzz_1, g_yyyyy_0_xxxzzz_1, g_yyyyy_0_xxyyy_1, g_yyyyy_0_xxyyyy_1, g_yyyyy_0_xxyyyz_1, g_yyyyy_0_xxyyz_1, g_yyyyy_0_xxyyzz_1, g_yyyyy_0_xxyzz_1, g_yyyyy_0_xxyzzz_1, g_yyyyy_0_xxzzz_1, g_yyyyy_0_xxzzzz_1, g_yyyyy_0_xyyyy_1, g_yyyyy_0_xyyyyy_1, g_yyyyy_0_xyyyyz_1, g_yyyyy_0_xyyyz_1, g_yyyyy_0_xyyyzz_1, g_yyyyy_0_xyyzz_1, g_yyyyy_0_xyyzzz_1, g_yyyyy_0_xyzzz_1, g_yyyyy_0_xyzzzz_1, g_yyyyy_0_xzzzz_1, g_yyyyy_0_xzzzzz_1, g_yyyyy_0_yyyyy_1, g_yyyyy_0_yyyyyy_1, g_yyyyy_0_yyyyyz_1, g_yyyyy_0_yyyyz_1, g_yyyyy_0_yyyyzz_1, g_yyyyy_0_yyyzz_1, g_yyyyy_0_yyyzzz_1, g_yyyyy_0_yyzzz_1, g_yyyyy_0_yyzzzz_1, g_yyyyy_0_yzzzz_1, g_yyyyy_0_yzzzzz_1, g_yyyyy_0_zzzzz_1, g_yyyyy_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyy_0_xxxxxx_0[i] = 6.0 * g_yyyyy_0_xxxxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxx_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxy_0[i] = 5.0 * g_yyyyy_0_xxxxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxz_0[i] = 5.0 * g_yyyyy_0_xxxxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxyy_0[i] = 4.0 * g_yyyyy_0_xxxyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxyz_0[i] = 4.0 * g_yyyyy_0_xxxyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxzz_0[i] = 4.0 * g_yyyyy_0_xxxzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyyy_0[i] = 3.0 * g_yyyyy_0_xxyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyyz_0[i] = 3.0 * g_yyyyy_0_xxyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyzz_0[i] = 3.0 * g_yyyyy_0_xxyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxzzz_0[i] = 3.0 * g_yyyyy_0_xxzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyyy_0[i] = 2.0 * g_yyyyy_0_xyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyyz_0[i] = 2.0 * g_yyyyy_0_xyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyzz_0[i] = 2.0 * g_yyyyy_0_xyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyzzz_0[i] = 2.0 * g_yyyyy_0_xyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxzzzz_0[i] = 2.0 * g_yyyyy_0_xzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyyy_0[i] = g_yyyyy_0_yyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyyz_0[i] = g_yyyyy_0_yyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyzz_0[i] = g_yyyyy_0_yyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyzzz_0[i] = g_yyyyy_0_yyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyzzzz_0[i] = g_yyyyy_0_yzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xzzzzz_0[i] = g_yyyyy_0_zzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyyy_0[i] = g_yyyyy_0_yyyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyyz_0[i] = g_yyyyy_0_yyyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyzz_0[i] = g_yyyyy_0_yyyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyzzz_0[i] = g_yyyyy_0_yyyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyzzzz_0[i] = g_yyyyy_0_yyzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yzzzzz_0[i] = g_yyyyy_0_yzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_zzzzzz_0[i] = g_yyyyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 448-476 components of targeted buffer : ISI

    auto g_xyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 448);

    auto g_xyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 449);

    auto g_xyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 450);

    auto g_xyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 451);

    auto g_xyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 452);

    auto g_xyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 453);

    auto g_xyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 454);

    auto g_xyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 455);

    auto g_xyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 456);

    auto g_xyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 457);

    auto g_xyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 458);

    auto g_xyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 459);

    auto g_xyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 460);

    auto g_xyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 461);

    auto g_xyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 462);

    auto g_xyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 463);

    auto g_xyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 464);

    auto g_xyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 465);

    auto g_xyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 466);

    auto g_xyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 467);

    auto g_xyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 468);

    auto g_xyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 469);

    auto g_xyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 470);

    auto g_xyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 471);

    auto g_xyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 472);

    auto g_xyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 473);

    auto g_xyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 474);

    auto g_xyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 475);

    #pragma omp simd aligned(g_xyyyy_0_xxxxxx_1, g_xyyyy_0_xxxxxy_1, g_xyyyy_0_xxxxyy_1, g_xyyyy_0_xxxyyy_1, g_xyyyy_0_xxyyyy_1, g_xyyyy_0_xyyyyy_1, g_xyyyyz_0_xxxxxx_0, g_xyyyyz_0_xxxxxy_0, g_xyyyyz_0_xxxxxz_0, g_xyyyyz_0_xxxxyy_0, g_xyyyyz_0_xxxxyz_0, g_xyyyyz_0_xxxxzz_0, g_xyyyyz_0_xxxyyy_0, g_xyyyyz_0_xxxyyz_0, g_xyyyyz_0_xxxyzz_0, g_xyyyyz_0_xxxzzz_0, g_xyyyyz_0_xxyyyy_0, g_xyyyyz_0_xxyyyz_0, g_xyyyyz_0_xxyyzz_0, g_xyyyyz_0_xxyzzz_0, g_xyyyyz_0_xxzzzz_0, g_xyyyyz_0_xyyyyy_0, g_xyyyyz_0_xyyyyz_0, g_xyyyyz_0_xyyyzz_0, g_xyyyyz_0_xyyzzz_0, g_xyyyyz_0_xyzzzz_0, g_xyyyyz_0_xzzzzz_0, g_xyyyyz_0_yyyyyy_0, g_xyyyyz_0_yyyyyz_0, g_xyyyyz_0_yyyyzz_0, g_xyyyyz_0_yyyzzz_0, g_xyyyyz_0_yyzzzz_0, g_xyyyyz_0_yzzzzz_0, g_xyyyyz_0_zzzzzz_0, g_yyyyz_0_xxxxxz_1, g_yyyyz_0_xxxxyz_1, g_yyyyz_0_xxxxz_1, g_yyyyz_0_xxxxzz_1, g_yyyyz_0_xxxyyz_1, g_yyyyz_0_xxxyz_1, g_yyyyz_0_xxxyzz_1, g_yyyyz_0_xxxzz_1, g_yyyyz_0_xxxzzz_1, g_yyyyz_0_xxyyyz_1, g_yyyyz_0_xxyyz_1, g_yyyyz_0_xxyyzz_1, g_yyyyz_0_xxyzz_1, g_yyyyz_0_xxyzzz_1, g_yyyyz_0_xxzzz_1, g_yyyyz_0_xxzzzz_1, g_yyyyz_0_xyyyyz_1, g_yyyyz_0_xyyyz_1, g_yyyyz_0_xyyyzz_1, g_yyyyz_0_xyyzz_1, g_yyyyz_0_xyyzzz_1, g_yyyyz_0_xyzzz_1, g_yyyyz_0_xyzzzz_1, g_yyyyz_0_xzzzz_1, g_yyyyz_0_xzzzzz_1, g_yyyyz_0_yyyyyy_1, g_yyyyz_0_yyyyyz_1, g_yyyyz_0_yyyyz_1, g_yyyyz_0_yyyyzz_1, g_yyyyz_0_yyyzz_1, g_yyyyz_0_yyyzzz_1, g_yyyyz_0_yyzzz_1, g_yyyyz_0_yyzzzz_1, g_yyyyz_0_yzzzz_1, g_yyyyz_0_yzzzzz_1, g_yyyyz_0_zzzzz_1, g_yyyyz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyz_0_xxxxxx_0[i] = g_xyyyy_0_xxxxxx_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxxy_0[i] = g_xyyyy_0_xxxxxy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxxz_0[i] = 5.0 * g_yyyyz_0_xxxxz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxxz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxyy_0[i] = g_xyyyy_0_xxxxyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxyz_0[i] = 4.0 * g_yyyyz_0_xxxyz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxzz_0[i] = 4.0 * g_yyyyz_0_xxxzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxyyy_0[i] = g_xyyyy_0_xxxyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxyyz_0[i] = 3.0 * g_yyyyz_0_xxyyz_1[i] * fi_acd_0 + g_yyyyz_0_xxxyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxyzz_0[i] = 3.0 * g_yyyyz_0_xxyzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxzzz_0[i] = 3.0 * g_yyyyz_0_xxzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyyyy_0[i] = g_xyyyy_0_xxyyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxyyyz_0[i] = 2.0 * g_yyyyz_0_xyyyz_1[i] * fi_acd_0 + g_yyyyz_0_xxyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyyzz_0[i] = 2.0 * g_yyyyz_0_xyyzz_1[i] * fi_acd_0 + g_yyyyz_0_xxyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyzzz_0[i] = 2.0 * g_yyyyz_0_xyzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxzzzz_0[i] = 2.0 * g_yyyyz_0_xzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyyyy_0[i] = g_xyyyy_0_xyyyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xyyyyz_0[i] = g_yyyyz_0_yyyyz_1[i] * fi_acd_0 + g_yyyyz_0_xyyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyyzz_0[i] = g_yyyyz_0_yyyzz_1[i] * fi_acd_0 + g_yyyyz_0_xyyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyzzz_0[i] = g_yyyyz_0_yyzzz_1[i] * fi_acd_0 + g_yyyyz_0_xyyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyzzzz_0[i] = g_yyyyz_0_yzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xyzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xzzzzz_0[i] = g_yyyyz_0_zzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyyy_0[i] = g_yyyyz_0_yyyyyy_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyyz_0[i] = g_yyyyz_0_yyyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyzz_0[i] = g_yyyyz_0_yyyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyzzz_0[i] = g_yyyyz_0_yyyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyzzzz_0[i] = g_yyyyz_0_yyzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yzzzzz_0[i] = g_yyyyz_0_yzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_zzzzzz_0[i] = g_yyyyz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 476-504 components of targeted buffer : ISI

    auto g_xyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 476);

    auto g_xyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 477);

    auto g_xyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 478);

    auto g_xyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 479);

    auto g_xyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 480);

    auto g_xyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 481);

    auto g_xyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 482);

    auto g_xyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 483);

    auto g_xyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 484);

    auto g_xyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 485);

    auto g_xyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 486);

    auto g_xyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 487);

    auto g_xyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 488);

    auto g_xyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 489);

    auto g_xyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 490);

    auto g_xyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 491);

    auto g_xyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 492);

    auto g_xyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 493);

    auto g_xyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 494);

    auto g_xyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 495);

    auto g_xyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 496);

    auto g_xyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 497);

    auto g_xyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 498);

    auto g_xyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 499);

    auto g_xyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 500);

    auto g_xyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 501);

    auto g_xyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 502);

    auto g_xyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 503);

    #pragma omp simd aligned(g_xyyyzz_0_xxxxxx_0, g_xyyyzz_0_xxxxxy_0, g_xyyyzz_0_xxxxxz_0, g_xyyyzz_0_xxxxyy_0, g_xyyyzz_0_xxxxyz_0, g_xyyyzz_0_xxxxzz_0, g_xyyyzz_0_xxxyyy_0, g_xyyyzz_0_xxxyyz_0, g_xyyyzz_0_xxxyzz_0, g_xyyyzz_0_xxxzzz_0, g_xyyyzz_0_xxyyyy_0, g_xyyyzz_0_xxyyyz_0, g_xyyyzz_0_xxyyzz_0, g_xyyyzz_0_xxyzzz_0, g_xyyyzz_0_xxzzzz_0, g_xyyyzz_0_xyyyyy_0, g_xyyyzz_0_xyyyyz_0, g_xyyyzz_0_xyyyzz_0, g_xyyyzz_0_xyyzzz_0, g_xyyyzz_0_xyzzzz_0, g_xyyyzz_0_xzzzzz_0, g_xyyyzz_0_yyyyyy_0, g_xyyyzz_0_yyyyyz_0, g_xyyyzz_0_yyyyzz_0, g_xyyyzz_0_yyyzzz_0, g_xyyyzz_0_yyzzzz_0, g_xyyyzz_0_yzzzzz_0, g_xyyyzz_0_zzzzzz_0, g_yyyzz_0_xxxxx_1, g_yyyzz_0_xxxxxx_1, g_yyyzz_0_xxxxxy_1, g_yyyzz_0_xxxxxz_1, g_yyyzz_0_xxxxy_1, g_yyyzz_0_xxxxyy_1, g_yyyzz_0_xxxxyz_1, g_yyyzz_0_xxxxz_1, g_yyyzz_0_xxxxzz_1, g_yyyzz_0_xxxyy_1, g_yyyzz_0_xxxyyy_1, g_yyyzz_0_xxxyyz_1, g_yyyzz_0_xxxyz_1, g_yyyzz_0_xxxyzz_1, g_yyyzz_0_xxxzz_1, g_yyyzz_0_xxxzzz_1, g_yyyzz_0_xxyyy_1, g_yyyzz_0_xxyyyy_1, g_yyyzz_0_xxyyyz_1, g_yyyzz_0_xxyyz_1, g_yyyzz_0_xxyyzz_1, g_yyyzz_0_xxyzz_1, g_yyyzz_0_xxyzzz_1, g_yyyzz_0_xxzzz_1, g_yyyzz_0_xxzzzz_1, g_yyyzz_0_xyyyy_1, g_yyyzz_0_xyyyyy_1, g_yyyzz_0_xyyyyz_1, g_yyyzz_0_xyyyz_1, g_yyyzz_0_xyyyzz_1, g_yyyzz_0_xyyzz_1, g_yyyzz_0_xyyzzz_1, g_yyyzz_0_xyzzz_1, g_yyyzz_0_xyzzzz_1, g_yyyzz_0_xzzzz_1, g_yyyzz_0_xzzzzz_1, g_yyyzz_0_yyyyy_1, g_yyyzz_0_yyyyyy_1, g_yyyzz_0_yyyyyz_1, g_yyyzz_0_yyyyz_1, g_yyyzz_0_yyyyzz_1, g_yyyzz_0_yyyzz_1, g_yyyzz_0_yyyzzz_1, g_yyyzz_0_yyzzz_1, g_yyyzz_0_yyzzzz_1, g_yyyzz_0_yzzzz_1, g_yyyzz_0_yzzzzz_1, g_yyyzz_0_zzzzz_1, g_yyyzz_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzz_0_xxxxxx_0[i] = 6.0 * g_yyyzz_0_xxxxx_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxx_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxy_0[i] = 5.0 * g_yyyzz_0_xxxxy_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxz_0[i] = 5.0 * g_yyyzz_0_xxxxz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxyy_0[i] = 4.0 * g_yyyzz_0_xxxyy_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxyz_0[i] = 4.0 * g_yyyzz_0_xxxyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxzz_0[i] = 4.0 * g_yyyzz_0_xxxzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyyy_0[i] = 3.0 * g_yyyzz_0_xxyyy_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyyz_0[i] = 3.0 * g_yyyzz_0_xxyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyzz_0[i] = 3.0 * g_yyyzz_0_xxyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxzzz_0[i] = 3.0 * g_yyyzz_0_xxzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyyy_0[i] = 2.0 * g_yyyzz_0_xyyyy_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyyz_0[i] = 2.0 * g_yyyzz_0_xyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyzz_0[i] = 2.0 * g_yyyzz_0_xyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyzzz_0[i] = 2.0 * g_yyyzz_0_xyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxzzzz_0[i] = 2.0 * g_yyyzz_0_xzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyyy_0[i] = g_yyyzz_0_yyyyy_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyyz_0[i] = g_yyyzz_0_yyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyzz_0[i] = g_yyyzz_0_yyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyzzz_0[i] = g_yyyzz_0_yyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyzzzz_0[i] = g_yyyzz_0_yzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xzzzzz_0[i] = g_yyyzz_0_zzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyyy_0[i] = g_yyyzz_0_yyyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyyz_0[i] = g_yyyzz_0_yyyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyzz_0[i] = g_yyyzz_0_yyyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyzzz_0[i] = g_yyyzz_0_yyyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyzzzz_0[i] = g_yyyzz_0_yyzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yzzzzz_0[i] = g_yyyzz_0_yzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_zzzzzz_0[i] = g_yyyzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 504-532 components of targeted buffer : ISI

    auto g_xyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 504);

    auto g_xyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 505);

    auto g_xyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 506);

    auto g_xyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 507);

    auto g_xyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 508);

    auto g_xyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 509);

    auto g_xyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 510);

    auto g_xyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 511);

    auto g_xyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 512);

    auto g_xyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 513);

    auto g_xyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 514);

    auto g_xyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 515);

    auto g_xyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 516);

    auto g_xyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 517);

    auto g_xyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 518);

    auto g_xyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 519);

    auto g_xyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 520);

    auto g_xyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 521);

    auto g_xyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 522);

    auto g_xyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 523);

    auto g_xyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 524);

    auto g_xyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 525);

    auto g_xyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 526);

    auto g_xyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 527);

    auto g_xyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 528);

    auto g_xyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 529);

    auto g_xyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 530);

    auto g_xyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 531);

    #pragma omp simd aligned(g_xyyzzz_0_xxxxxx_0, g_xyyzzz_0_xxxxxy_0, g_xyyzzz_0_xxxxxz_0, g_xyyzzz_0_xxxxyy_0, g_xyyzzz_0_xxxxyz_0, g_xyyzzz_0_xxxxzz_0, g_xyyzzz_0_xxxyyy_0, g_xyyzzz_0_xxxyyz_0, g_xyyzzz_0_xxxyzz_0, g_xyyzzz_0_xxxzzz_0, g_xyyzzz_0_xxyyyy_0, g_xyyzzz_0_xxyyyz_0, g_xyyzzz_0_xxyyzz_0, g_xyyzzz_0_xxyzzz_0, g_xyyzzz_0_xxzzzz_0, g_xyyzzz_0_xyyyyy_0, g_xyyzzz_0_xyyyyz_0, g_xyyzzz_0_xyyyzz_0, g_xyyzzz_0_xyyzzz_0, g_xyyzzz_0_xyzzzz_0, g_xyyzzz_0_xzzzzz_0, g_xyyzzz_0_yyyyyy_0, g_xyyzzz_0_yyyyyz_0, g_xyyzzz_0_yyyyzz_0, g_xyyzzz_0_yyyzzz_0, g_xyyzzz_0_yyzzzz_0, g_xyyzzz_0_yzzzzz_0, g_xyyzzz_0_zzzzzz_0, g_yyzzz_0_xxxxx_1, g_yyzzz_0_xxxxxx_1, g_yyzzz_0_xxxxxy_1, g_yyzzz_0_xxxxxz_1, g_yyzzz_0_xxxxy_1, g_yyzzz_0_xxxxyy_1, g_yyzzz_0_xxxxyz_1, g_yyzzz_0_xxxxz_1, g_yyzzz_0_xxxxzz_1, g_yyzzz_0_xxxyy_1, g_yyzzz_0_xxxyyy_1, g_yyzzz_0_xxxyyz_1, g_yyzzz_0_xxxyz_1, g_yyzzz_0_xxxyzz_1, g_yyzzz_0_xxxzz_1, g_yyzzz_0_xxxzzz_1, g_yyzzz_0_xxyyy_1, g_yyzzz_0_xxyyyy_1, g_yyzzz_0_xxyyyz_1, g_yyzzz_0_xxyyz_1, g_yyzzz_0_xxyyzz_1, g_yyzzz_0_xxyzz_1, g_yyzzz_0_xxyzzz_1, g_yyzzz_0_xxzzz_1, g_yyzzz_0_xxzzzz_1, g_yyzzz_0_xyyyy_1, g_yyzzz_0_xyyyyy_1, g_yyzzz_0_xyyyyz_1, g_yyzzz_0_xyyyz_1, g_yyzzz_0_xyyyzz_1, g_yyzzz_0_xyyzz_1, g_yyzzz_0_xyyzzz_1, g_yyzzz_0_xyzzz_1, g_yyzzz_0_xyzzzz_1, g_yyzzz_0_xzzzz_1, g_yyzzz_0_xzzzzz_1, g_yyzzz_0_yyyyy_1, g_yyzzz_0_yyyyyy_1, g_yyzzz_0_yyyyyz_1, g_yyzzz_0_yyyyz_1, g_yyzzz_0_yyyyzz_1, g_yyzzz_0_yyyzz_1, g_yyzzz_0_yyyzzz_1, g_yyzzz_0_yyzzz_1, g_yyzzz_0_yyzzzz_1, g_yyzzz_0_yzzzz_1, g_yyzzz_0_yzzzzz_1, g_yyzzz_0_zzzzz_1, g_yyzzz_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzz_0_xxxxxx_0[i] = 6.0 * g_yyzzz_0_xxxxx_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxx_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxy_0[i] = 5.0 * g_yyzzz_0_xxxxy_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxz_0[i] = 5.0 * g_yyzzz_0_xxxxz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxyy_0[i] = 4.0 * g_yyzzz_0_xxxyy_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxyz_0[i] = 4.0 * g_yyzzz_0_xxxyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxzz_0[i] = 4.0 * g_yyzzz_0_xxxzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyyy_0[i] = 3.0 * g_yyzzz_0_xxyyy_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyyz_0[i] = 3.0 * g_yyzzz_0_xxyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyzz_0[i] = 3.0 * g_yyzzz_0_xxyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxzzz_0[i] = 3.0 * g_yyzzz_0_xxzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyyy_0[i] = 2.0 * g_yyzzz_0_xyyyy_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyyz_0[i] = 2.0 * g_yyzzz_0_xyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyzz_0[i] = 2.0 * g_yyzzz_0_xyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyzzz_0[i] = 2.0 * g_yyzzz_0_xyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxzzzz_0[i] = 2.0 * g_yyzzz_0_xzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyyy_0[i] = g_yyzzz_0_yyyyy_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyyz_0[i] = g_yyzzz_0_yyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyzz_0[i] = g_yyzzz_0_yyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyzzz_0[i] = g_yyzzz_0_yyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyzzzz_0[i] = g_yyzzz_0_yzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xzzzzz_0[i] = g_yyzzz_0_zzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyyy_0[i] = g_yyzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyyz_0[i] = g_yyzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyzz_0[i] = g_yyzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyzzz_0[i] = g_yyzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyzzzz_0[i] = g_yyzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yzzzzz_0[i] = g_yyzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_zzzzzz_0[i] = g_yyzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 532-560 components of targeted buffer : ISI

    auto g_xyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 532);

    auto g_xyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 533);

    auto g_xyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 534);

    auto g_xyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 535);

    auto g_xyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 536);

    auto g_xyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 537);

    auto g_xyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 538);

    auto g_xyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 539);

    auto g_xyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 540);

    auto g_xyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 541);

    auto g_xyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 542);

    auto g_xyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 543);

    auto g_xyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 544);

    auto g_xyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 545);

    auto g_xyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 546);

    auto g_xyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 547);

    auto g_xyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 548);

    auto g_xyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 549);

    auto g_xyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 550);

    auto g_xyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 551);

    auto g_xyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 552);

    auto g_xyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 553);

    auto g_xyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 554);

    auto g_xyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 555);

    auto g_xyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 556);

    auto g_xyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 557);

    auto g_xyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 558);

    auto g_xyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 559);

    #pragma omp simd aligned(g_xyzzzz_0_xxxxxx_0, g_xyzzzz_0_xxxxxy_0, g_xyzzzz_0_xxxxxz_0, g_xyzzzz_0_xxxxyy_0, g_xyzzzz_0_xxxxyz_0, g_xyzzzz_0_xxxxzz_0, g_xyzzzz_0_xxxyyy_0, g_xyzzzz_0_xxxyyz_0, g_xyzzzz_0_xxxyzz_0, g_xyzzzz_0_xxxzzz_0, g_xyzzzz_0_xxyyyy_0, g_xyzzzz_0_xxyyyz_0, g_xyzzzz_0_xxyyzz_0, g_xyzzzz_0_xxyzzz_0, g_xyzzzz_0_xxzzzz_0, g_xyzzzz_0_xyyyyy_0, g_xyzzzz_0_xyyyyz_0, g_xyzzzz_0_xyyyzz_0, g_xyzzzz_0_xyyzzz_0, g_xyzzzz_0_xyzzzz_0, g_xyzzzz_0_xzzzzz_0, g_xyzzzz_0_yyyyyy_0, g_xyzzzz_0_yyyyyz_0, g_xyzzzz_0_yyyyzz_0, g_xyzzzz_0_yyyzzz_0, g_xyzzzz_0_yyzzzz_0, g_xyzzzz_0_yzzzzz_0, g_xyzzzz_0_zzzzzz_0, g_xzzzz_0_xxxxxx_1, g_xzzzz_0_xxxxxz_1, g_xzzzz_0_xxxxzz_1, g_xzzzz_0_xxxzzz_1, g_xzzzz_0_xxzzzz_1, g_xzzzz_0_xzzzzz_1, g_yzzzz_0_xxxxxy_1, g_yzzzz_0_xxxxy_1, g_yzzzz_0_xxxxyy_1, g_yzzzz_0_xxxxyz_1, g_yzzzz_0_xxxyy_1, g_yzzzz_0_xxxyyy_1, g_yzzzz_0_xxxyyz_1, g_yzzzz_0_xxxyz_1, g_yzzzz_0_xxxyzz_1, g_yzzzz_0_xxyyy_1, g_yzzzz_0_xxyyyy_1, g_yzzzz_0_xxyyyz_1, g_yzzzz_0_xxyyz_1, g_yzzzz_0_xxyyzz_1, g_yzzzz_0_xxyzz_1, g_yzzzz_0_xxyzzz_1, g_yzzzz_0_xyyyy_1, g_yzzzz_0_xyyyyy_1, g_yzzzz_0_xyyyyz_1, g_yzzzz_0_xyyyz_1, g_yzzzz_0_xyyyzz_1, g_yzzzz_0_xyyzz_1, g_yzzzz_0_xyyzzz_1, g_yzzzz_0_xyzzz_1, g_yzzzz_0_xyzzzz_1, g_yzzzz_0_yyyyy_1, g_yzzzz_0_yyyyyy_1, g_yzzzz_0_yyyyyz_1, g_yzzzz_0_yyyyz_1, g_yzzzz_0_yyyyzz_1, g_yzzzz_0_yyyzz_1, g_yzzzz_0_yyyzzz_1, g_yzzzz_0_yyzzz_1, g_yzzzz_0_yyzzzz_1, g_yzzzz_0_yzzzz_1, g_yzzzz_0_yzzzzz_1, g_yzzzz_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzz_0_xxxxxx_0[i] = g_xzzzz_0_xxxxxx_1[i] * wa_y[i];

        g_xyzzzz_0_xxxxxy_0[i] = 5.0 * g_yzzzz_0_xxxxy_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxxz_0[i] = g_xzzzz_0_xxxxxz_1[i] * wa_y[i];

        g_xyzzzz_0_xxxxyy_0[i] = 4.0 * g_yzzzz_0_xxxyy_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxyz_0[i] = 4.0 * g_yzzzz_0_xxxyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxzz_0[i] = g_xzzzz_0_xxxxzz_1[i] * wa_y[i];

        g_xyzzzz_0_xxxyyy_0[i] = 3.0 * g_yzzzz_0_xxyyy_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxyyz_0[i] = 3.0 * g_yzzzz_0_xxyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxyzz_0[i] = 3.0 * g_yzzzz_0_xxyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxzzz_0[i] = g_xzzzz_0_xxxzzz_1[i] * wa_y[i];

        g_xyzzzz_0_xxyyyy_0[i] = 2.0 * g_yzzzz_0_xyyyy_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxyyyz_0[i] = 2.0 * g_yzzzz_0_xyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxyyzz_0[i] = 2.0 * g_yzzzz_0_xyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxyzzz_0[i] = 2.0 * g_yzzzz_0_xyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxzzzz_0[i] = g_xzzzz_0_xxzzzz_1[i] * wa_y[i];

        g_xyzzzz_0_xyyyyy_0[i] = g_yzzzz_0_yyyyy_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xyyyyz_0[i] = g_yzzzz_0_yyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xyyyzz_0[i] = g_yzzzz_0_yyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xyyzzz_0[i] = g_yzzzz_0_yyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xyzzzz_0[i] = g_yzzzz_0_yzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xzzzzz_0[i] = g_xzzzz_0_xzzzzz_1[i] * wa_y[i];

        g_xyzzzz_0_yyyyyy_0[i] = g_yzzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_yyyyyz_0[i] = g_yzzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_yyyyzz_0[i] = g_yzzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_yyyzzz_0[i] = g_yzzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_yyzzzz_0[i] = g_yzzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_yzzzzz_0[i] = g_yzzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_zzzzzz_0[i] = g_yzzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 560-588 components of targeted buffer : ISI

    auto g_xzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 560);

    auto g_xzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 561);

    auto g_xzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 562);

    auto g_xzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 563);

    auto g_xzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 564);

    auto g_xzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 565);

    auto g_xzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 566);

    auto g_xzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 567);

    auto g_xzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 568);

    auto g_xzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 569);

    auto g_xzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 570);

    auto g_xzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 571);

    auto g_xzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 572);

    auto g_xzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 573);

    auto g_xzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 574);

    auto g_xzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 575);

    auto g_xzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 576);

    auto g_xzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 577);

    auto g_xzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 578);

    auto g_xzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 579);

    auto g_xzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 580);

    auto g_xzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 581);

    auto g_xzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 582);

    auto g_xzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 583);

    auto g_xzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 584);

    auto g_xzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 585);

    auto g_xzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 586);

    auto g_xzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 587);

    #pragma omp simd aligned(g_xzzzzz_0_xxxxxx_0, g_xzzzzz_0_xxxxxy_0, g_xzzzzz_0_xxxxxz_0, g_xzzzzz_0_xxxxyy_0, g_xzzzzz_0_xxxxyz_0, g_xzzzzz_0_xxxxzz_0, g_xzzzzz_0_xxxyyy_0, g_xzzzzz_0_xxxyyz_0, g_xzzzzz_0_xxxyzz_0, g_xzzzzz_0_xxxzzz_0, g_xzzzzz_0_xxyyyy_0, g_xzzzzz_0_xxyyyz_0, g_xzzzzz_0_xxyyzz_0, g_xzzzzz_0_xxyzzz_0, g_xzzzzz_0_xxzzzz_0, g_xzzzzz_0_xyyyyy_0, g_xzzzzz_0_xyyyyz_0, g_xzzzzz_0_xyyyzz_0, g_xzzzzz_0_xyyzzz_0, g_xzzzzz_0_xyzzzz_0, g_xzzzzz_0_xzzzzz_0, g_xzzzzz_0_yyyyyy_0, g_xzzzzz_0_yyyyyz_0, g_xzzzzz_0_yyyyzz_0, g_xzzzzz_0_yyyzzz_0, g_xzzzzz_0_yyzzzz_0, g_xzzzzz_0_yzzzzz_0, g_xzzzzz_0_zzzzzz_0, g_zzzzz_0_xxxxx_1, g_zzzzz_0_xxxxxx_1, g_zzzzz_0_xxxxxy_1, g_zzzzz_0_xxxxxz_1, g_zzzzz_0_xxxxy_1, g_zzzzz_0_xxxxyy_1, g_zzzzz_0_xxxxyz_1, g_zzzzz_0_xxxxz_1, g_zzzzz_0_xxxxzz_1, g_zzzzz_0_xxxyy_1, g_zzzzz_0_xxxyyy_1, g_zzzzz_0_xxxyyz_1, g_zzzzz_0_xxxyz_1, g_zzzzz_0_xxxyzz_1, g_zzzzz_0_xxxzz_1, g_zzzzz_0_xxxzzz_1, g_zzzzz_0_xxyyy_1, g_zzzzz_0_xxyyyy_1, g_zzzzz_0_xxyyyz_1, g_zzzzz_0_xxyyz_1, g_zzzzz_0_xxyyzz_1, g_zzzzz_0_xxyzz_1, g_zzzzz_0_xxyzzz_1, g_zzzzz_0_xxzzz_1, g_zzzzz_0_xxzzzz_1, g_zzzzz_0_xyyyy_1, g_zzzzz_0_xyyyyy_1, g_zzzzz_0_xyyyyz_1, g_zzzzz_0_xyyyz_1, g_zzzzz_0_xyyyzz_1, g_zzzzz_0_xyyzz_1, g_zzzzz_0_xyyzzz_1, g_zzzzz_0_xyzzz_1, g_zzzzz_0_xyzzzz_1, g_zzzzz_0_xzzzz_1, g_zzzzz_0_xzzzzz_1, g_zzzzz_0_yyyyy_1, g_zzzzz_0_yyyyyy_1, g_zzzzz_0_yyyyyz_1, g_zzzzz_0_yyyyz_1, g_zzzzz_0_yyyyzz_1, g_zzzzz_0_yyyzz_1, g_zzzzz_0_yyyzzz_1, g_zzzzz_0_yyzzz_1, g_zzzzz_0_yyzzzz_1, g_zzzzz_0_yzzzz_1, g_zzzzz_0_yzzzzz_1, g_zzzzz_0_zzzzz_1, g_zzzzz_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzz_0_xxxxxx_0[i] = 6.0 * g_zzzzz_0_xxxxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxx_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxy_0[i] = 5.0 * g_zzzzz_0_xxxxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxz_0[i] = 5.0 * g_zzzzz_0_xxxxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxyy_0[i] = 4.0 * g_zzzzz_0_xxxyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxyz_0[i] = 4.0 * g_zzzzz_0_xxxyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxzz_0[i] = 4.0 * g_zzzzz_0_xxxzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyyy_0[i] = 3.0 * g_zzzzz_0_xxyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyyz_0[i] = 3.0 * g_zzzzz_0_xxyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyzz_0[i] = 3.0 * g_zzzzz_0_xxyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxzzz_0[i] = 3.0 * g_zzzzz_0_xxzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyyy_0[i] = 2.0 * g_zzzzz_0_xyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyyz_0[i] = 2.0 * g_zzzzz_0_xyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyzz_0[i] = 2.0 * g_zzzzz_0_xyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyzzz_0[i] = 2.0 * g_zzzzz_0_xyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxzzzz_0[i] = 2.0 * g_zzzzz_0_xzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyyy_0[i] = g_zzzzz_0_yyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyyz_0[i] = g_zzzzz_0_yyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyzz_0[i] = g_zzzzz_0_yyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyzzz_0[i] = g_zzzzz_0_yyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyzzzz_0[i] = g_zzzzz_0_yzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xzzzzz_0[i] = g_zzzzz_0_zzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyyy_0[i] = g_zzzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyyz_0[i] = g_zzzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyzz_0[i] = g_zzzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyzzz_0[i] = g_zzzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyzzzz_0[i] = g_zzzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yzzzzz_0[i] = g_zzzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_zzzzzz_0[i] = g_zzzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 588-616 components of targeted buffer : ISI

    auto g_yyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 588);

    auto g_yyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 589);

    auto g_yyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 590);

    auto g_yyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 591);

    auto g_yyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 592);

    auto g_yyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 593);

    auto g_yyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 594);

    auto g_yyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 595);

    auto g_yyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 596);

    auto g_yyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 597);

    auto g_yyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 598);

    auto g_yyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 599);

    auto g_yyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 600);

    auto g_yyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 601);

    auto g_yyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 602);

    auto g_yyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 603);

    auto g_yyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 604);

    auto g_yyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 605);

    auto g_yyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 606);

    auto g_yyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 607);

    auto g_yyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 608);

    auto g_yyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 609);

    auto g_yyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 610);

    auto g_yyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 611);

    auto g_yyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 612);

    auto g_yyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 613);

    auto g_yyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 614);

    auto g_yyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 615);

    #pragma omp simd aligned(g_yyyy_0_xxxxxx_0, g_yyyy_0_xxxxxx_1, g_yyyy_0_xxxxxy_0, g_yyyy_0_xxxxxy_1, g_yyyy_0_xxxxxz_0, g_yyyy_0_xxxxxz_1, g_yyyy_0_xxxxyy_0, g_yyyy_0_xxxxyy_1, g_yyyy_0_xxxxyz_0, g_yyyy_0_xxxxyz_1, g_yyyy_0_xxxxzz_0, g_yyyy_0_xxxxzz_1, g_yyyy_0_xxxyyy_0, g_yyyy_0_xxxyyy_1, g_yyyy_0_xxxyyz_0, g_yyyy_0_xxxyyz_1, g_yyyy_0_xxxyzz_0, g_yyyy_0_xxxyzz_1, g_yyyy_0_xxxzzz_0, g_yyyy_0_xxxzzz_1, g_yyyy_0_xxyyyy_0, g_yyyy_0_xxyyyy_1, g_yyyy_0_xxyyyz_0, g_yyyy_0_xxyyyz_1, g_yyyy_0_xxyyzz_0, g_yyyy_0_xxyyzz_1, g_yyyy_0_xxyzzz_0, g_yyyy_0_xxyzzz_1, g_yyyy_0_xxzzzz_0, g_yyyy_0_xxzzzz_1, g_yyyy_0_xyyyyy_0, g_yyyy_0_xyyyyy_1, g_yyyy_0_xyyyyz_0, g_yyyy_0_xyyyyz_1, g_yyyy_0_xyyyzz_0, g_yyyy_0_xyyyzz_1, g_yyyy_0_xyyzzz_0, g_yyyy_0_xyyzzz_1, g_yyyy_0_xyzzzz_0, g_yyyy_0_xyzzzz_1, g_yyyy_0_xzzzzz_0, g_yyyy_0_xzzzzz_1, g_yyyy_0_yyyyyy_0, g_yyyy_0_yyyyyy_1, g_yyyy_0_yyyyyz_0, g_yyyy_0_yyyyyz_1, g_yyyy_0_yyyyzz_0, g_yyyy_0_yyyyzz_1, g_yyyy_0_yyyzzz_0, g_yyyy_0_yyyzzz_1, g_yyyy_0_yyzzzz_0, g_yyyy_0_yyzzzz_1, g_yyyy_0_yzzzzz_0, g_yyyy_0_yzzzzz_1, g_yyyy_0_zzzzzz_0, g_yyyy_0_zzzzzz_1, g_yyyyy_0_xxxxx_1, g_yyyyy_0_xxxxxx_1, g_yyyyy_0_xxxxxy_1, g_yyyyy_0_xxxxxz_1, g_yyyyy_0_xxxxy_1, g_yyyyy_0_xxxxyy_1, g_yyyyy_0_xxxxyz_1, g_yyyyy_0_xxxxz_1, g_yyyyy_0_xxxxzz_1, g_yyyyy_0_xxxyy_1, g_yyyyy_0_xxxyyy_1, g_yyyyy_0_xxxyyz_1, g_yyyyy_0_xxxyz_1, g_yyyyy_0_xxxyzz_1, g_yyyyy_0_xxxzz_1, g_yyyyy_0_xxxzzz_1, g_yyyyy_0_xxyyy_1, g_yyyyy_0_xxyyyy_1, g_yyyyy_0_xxyyyz_1, g_yyyyy_0_xxyyz_1, g_yyyyy_0_xxyyzz_1, g_yyyyy_0_xxyzz_1, g_yyyyy_0_xxyzzz_1, g_yyyyy_0_xxzzz_1, g_yyyyy_0_xxzzzz_1, g_yyyyy_0_xyyyy_1, g_yyyyy_0_xyyyyy_1, g_yyyyy_0_xyyyyz_1, g_yyyyy_0_xyyyz_1, g_yyyyy_0_xyyyzz_1, g_yyyyy_0_xyyzz_1, g_yyyyy_0_xyyzzz_1, g_yyyyy_0_xyzzz_1, g_yyyyy_0_xyzzzz_1, g_yyyyy_0_xzzzz_1, g_yyyyy_0_xzzzzz_1, g_yyyyy_0_yyyyy_1, g_yyyyy_0_yyyyyy_1, g_yyyyy_0_yyyyyz_1, g_yyyyy_0_yyyyz_1, g_yyyyy_0_yyyyzz_1, g_yyyyy_0_yyyzz_1, g_yyyyy_0_yyyzzz_1, g_yyyyy_0_yyzzz_1, g_yyyyy_0_yyzzzz_1, g_yyyyy_0_yzzzz_1, g_yyyyy_0_yzzzzz_1, g_yyyyy_0_zzzzz_1, g_yyyyy_0_zzzzzz_1, g_yyyyyy_0_xxxxxx_0, g_yyyyyy_0_xxxxxy_0, g_yyyyyy_0_xxxxxz_0, g_yyyyyy_0_xxxxyy_0, g_yyyyyy_0_xxxxyz_0, g_yyyyyy_0_xxxxzz_0, g_yyyyyy_0_xxxyyy_0, g_yyyyyy_0_xxxyyz_0, g_yyyyyy_0_xxxyzz_0, g_yyyyyy_0_xxxzzz_0, g_yyyyyy_0_xxyyyy_0, g_yyyyyy_0_xxyyyz_0, g_yyyyyy_0_xxyyzz_0, g_yyyyyy_0_xxyzzz_0, g_yyyyyy_0_xxzzzz_0, g_yyyyyy_0_xyyyyy_0, g_yyyyyy_0_xyyyyz_0, g_yyyyyy_0_xyyyzz_0, g_yyyyyy_0_xyyzzz_0, g_yyyyyy_0_xyzzzz_0, g_yyyyyy_0_xzzzzz_0, g_yyyyyy_0_yyyyyy_0, g_yyyyyy_0_yyyyyz_0, g_yyyyyy_0_yyyyzz_0, g_yyyyyy_0_yyyzzz_0, g_yyyyyy_0_yyzzzz_0, g_yyyyyy_0_yzzzzz_0, g_yyyyyy_0_zzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyy_0_xxxxxx_0[i] = 5.0 * g_yyyy_0_xxxxxx_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxx_1[i] * fz_be_0 + g_yyyyy_0_xxxxxx_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxy_0[i] = 5.0 * g_yyyy_0_xxxxxy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxy_1[i] * fz_be_0 + g_yyyyy_0_xxxxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxz_0[i] = 5.0 * g_yyyy_0_xxxxxz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxz_1[i] * fz_be_0 + g_yyyyy_0_xxxxxz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxyy_0[i] = 5.0 * g_yyyy_0_xxxxyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxyy_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxxxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxyz_0[i] = 5.0 * g_yyyy_0_xxxxyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxyz_1[i] * fz_be_0 + g_yyyyy_0_xxxxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxzz_0[i] = 5.0 * g_yyyy_0_xxxxzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxzz_1[i] * fz_be_0 + g_yyyyy_0_xxxxzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyyy_0[i] = 5.0 * g_yyyy_0_xxxyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xxxyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyyz_0[i] = 5.0 * g_yyyy_0_xxxyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxxyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyzz_0[i] = 5.0 * g_yyyy_0_xxxyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyzz_1[i] * fz_be_0 + g_yyyyy_0_xxxzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxzzz_0[i] = 5.0 * g_yyyy_0_xxxzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxzzz_1[i] * fz_be_0 + g_yyyyy_0_xxxzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyyy_0[i] = 5.0 * g_yyyy_0_xxyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyyy_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_xxyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyyz_0[i] = 5.0 * g_yyyy_0_xxyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xxyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyzz_0[i] = 5.0 * g_yyyy_0_xxyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyzzz_0[i] = 5.0 * g_yyyy_0_xxyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyzzz_1[i] * fz_be_0 + g_yyyyy_0_xxzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxzzzz_0[i] = 5.0 * g_yyyy_0_xxzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxzzzz_1[i] * fz_be_0 + g_yyyyy_0_xxzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyyy_0[i] = 5.0 * g_yyyy_0_xyyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyyy_1[i] * fz_be_0 + 5.0 * g_yyyyy_0_xyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyyz_0[i] = 5.0 * g_yyyy_0_xyyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_xyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyzz_0[i] = 5.0 * g_yyyy_0_xyyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyzzz_0[i] = 5.0 * g_yyyy_0_xyyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyzzzz_0[i] = 5.0 * g_yyyy_0_xyzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyzzzz_1[i] * fz_be_0 + g_yyyyy_0_xzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xzzzzz_0[i] = 5.0 * g_yyyy_0_xzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xzzzzz_1[i] * fz_be_0 + g_yyyyy_0_xzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyyy_0[i] = 5.0 * g_yyyy_0_yyyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyyy_1[i] * fz_be_0 + 6.0 * g_yyyyy_0_yyyyy_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyyz_0[i] = 5.0 * g_yyyy_0_yyyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyyy_0_yyyyz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyzz_0[i] = 5.0 * g_yyyy_0_yyyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_yyyzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyzzz_0[i] = 5.0 * g_yyyy_0_yyyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_yyzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyzzzz_0[i] = 5.0 * g_yyyy_0_yyzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_yzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yzzzzz_0[i] = 5.0 * g_yyyy_0_yzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yzzzzz_1[i] * fz_be_0 + g_yyyyy_0_zzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_zzzzzz_0[i] = 5.0 * g_yyyy_0_zzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_zzzzzz_1[i] * fz_be_0 + g_yyyyy_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 616-644 components of targeted buffer : ISI

    auto g_yyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 616);

    auto g_yyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 617);

    auto g_yyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 618);

    auto g_yyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 619);

    auto g_yyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 620);

    auto g_yyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 621);

    auto g_yyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 622);

    auto g_yyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 623);

    auto g_yyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 624);

    auto g_yyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 625);

    auto g_yyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 626);

    auto g_yyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 627);

    auto g_yyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 628);

    auto g_yyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 629);

    auto g_yyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 630);

    auto g_yyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 631);

    auto g_yyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 632);

    auto g_yyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 633);

    auto g_yyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 634);

    auto g_yyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 635);

    auto g_yyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 636);

    auto g_yyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 637);

    auto g_yyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 638);

    auto g_yyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 639);

    auto g_yyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 640);

    auto g_yyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 641);

    auto g_yyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 642);

    auto g_yyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 643);

    #pragma omp simd aligned(g_yyyyy_0_xxxxx_1, g_yyyyy_0_xxxxxx_1, g_yyyyy_0_xxxxxy_1, g_yyyyy_0_xxxxxz_1, g_yyyyy_0_xxxxy_1, g_yyyyy_0_xxxxyy_1, g_yyyyy_0_xxxxyz_1, g_yyyyy_0_xxxxz_1, g_yyyyy_0_xxxxzz_1, g_yyyyy_0_xxxyy_1, g_yyyyy_0_xxxyyy_1, g_yyyyy_0_xxxyyz_1, g_yyyyy_0_xxxyz_1, g_yyyyy_0_xxxyzz_1, g_yyyyy_0_xxxzz_1, g_yyyyy_0_xxxzzz_1, g_yyyyy_0_xxyyy_1, g_yyyyy_0_xxyyyy_1, g_yyyyy_0_xxyyyz_1, g_yyyyy_0_xxyyz_1, g_yyyyy_0_xxyyzz_1, g_yyyyy_0_xxyzz_1, g_yyyyy_0_xxyzzz_1, g_yyyyy_0_xxzzz_1, g_yyyyy_0_xxzzzz_1, g_yyyyy_0_xyyyy_1, g_yyyyy_0_xyyyyy_1, g_yyyyy_0_xyyyyz_1, g_yyyyy_0_xyyyz_1, g_yyyyy_0_xyyyzz_1, g_yyyyy_0_xyyzz_1, g_yyyyy_0_xyyzzz_1, g_yyyyy_0_xyzzz_1, g_yyyyy_0_xyzzzz_1, g_yyyyy_0_xzzzz_1, g_yyyyy_0_xzzzzz_1, g_yyyyy_0_yyyyy_1, g_yyyyy_0_yyyyyy_1, g_yyyyy_0_yyyyyz_1, g_yyyyy_0_yyyyz_1, g_yyyyy_0_yyyyzz_1, g_yyyyy_0_yyyzz_1, g_yyyyy_0_yyyzzz_1, g_yyyyy_0_yyzzz_1, g_yyyyy_0_yyzzzz_1, g_yyyyy_0_yzzzz_1, g_yyyyy_0_yzzzzz_1, g_yyyyy_0_zzzzz_1, g_yyyyy_0_zzzzzz_1, g_yyyyyz_0_xxxxxx_0, g_yyyyyz_0_xxxxxy_0, g_yyyyyz_0_xxxxxz_0, g_yyyyyz_0_xxxxyy_0, g_yyyyyz_0_xxxxyz_0, g_yyyyyz_0_xxxxzz_0, g_yyyyyz_0_xxxyyy_0, g_yyyyyz_0_xxxyyz_0, g_yyyyyz_0_xxxyzz_0, g_yyyyyz_0_xxxzzz_0, g_yyyyyz_0_xxyyyy_0, g_yyyyyz_0_xxyyyz_0, g_yyyyyz_0_xxyyzz_0, g_yyyyyz_0_xxyzzz_0, g_yyyyyz_0_xxzzzz_0, g_yyyyyz_0_xyyyyy_0, g_yyyyyz_0_xyyyyz_0, g_yyyyyz_0_xyyyzz_0, g_yyyyyz_0_xyyzzz_0, g_yyyyyz_0_xyzzzz_0, g_yyyyyz_0_xzzzzz_0, g_yyyyyz_0_yyyyyy_0, g_yyyyyz_0_yyyyyz_0, g_yyyyyz_0_yyyyzz_0, g_yyyyyz_0_yyyzzz_0, g_yyyyyz_0_yyzzzz_0, g_yyyyyz_0_yzzzzz_0, g_yyyyyz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyz_0_xxxxxx_0[i] = g_yyyyy_0_xxxxxx_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxy_0[i] = g_yyyyy_0_xxxxxy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxz_0[i] = g_yyyyy_0_xxxxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxyy_0[i] = g_yyyyy_0_xxxxyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxyz_0[i] = g_yyyyy_0_xxxxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxzz_0[i] = 2.0 * g_yyyyy_0_xxxxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyyy_0[i] = g_yyyyy_0_xxxyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyyz_0[i] = g_yyyyy_0_xxxyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyzz_0[i] = 2.0 * g_yyyyy_0_xxxyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxzzz_0[i] = 3.0 * g_yyyyy_0_xxxzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyyy_0[i] = g_yyyyy_0_xxyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyyz_0[i] = g_yyyyy_0_xxyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyzz_0[i] = 2.0 * g_yyyyy_0_xxyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyzzz_0[i] = 3.0 * g_yyyyy_0_xxyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxzzzz_0[i] = 4.0 * g_yyyyy_0_xxzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyyy_0[i] = g_yyyyy_0_xyyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyyz_0[i] = g_yyyyy_0_xyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyzz_0[i] = 2.0 * g_yyyyy_0_xyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyzzz_0[i] = 3.0 * g_yyyyy_0_xyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyzzzz_0[i] = 4.0 * g_yyyyy_0_xyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xzzzzz_0[i] = 5.0 * g_yyyyy_0_xzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyyy_0[i] = g_yyyyy_0_yyyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyyz_0[i] = g_yyyyy_0_yyyyy_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyzz_0[i] = 2.0 * g_yyyyy_0_yyyyz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyzzz_0[i] = 3.0 * g_yyyyy_0_yyyzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyzzzz_0[i] = 4.0 * g_yyyyy_0_yyzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yzzzzz_0[i] = 5.0 * g_yyyyy_0_yzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_zzzzzz_0[i] = 6.0 * g_yyyyy_0_zzzzz_1[i] * fi_acd_0 + g_yyyyy_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 644-672 components of targeted buffer : ISI

    auto g_yyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 644);

    auto g_yyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 645);

    auto g_yyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 646);

    auto g_yyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 647);

    auto g_yyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 648);

    auto g_yyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 649);

    auto g_yyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 650);

    auto g_yyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 651);

    auto g_yyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 652);

    auto g_yyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 653);

    auto g_yyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 654);

    auto g_yyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 655);

    auto g_yyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 656);

    auto g_yyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 657);

    auto g_yyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 658);

    auto g_yyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 659);

    auto g_yyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 660);

    auto g_yyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 661);

    auto g_yyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 662);

    auto g_yyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 663);

    auto g_yyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 664);

    auto g_yyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 665);

    auto g_yyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 666);

    auto g_yyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 667);

    auto g_yyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 668);

    auto g_yyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 669);

    auto g_yyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 670);

    auto g_yyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 671);

    #pragma omp simd aligned(g_yyyy_0_xxxxxy_0, g_yyyy_0_xxxxxy_1, g_yyyy_0_xxxxyy_0, g_yyyy_0_xxxxyy_1, g_yyyy_0_xxxyyy_0, g_yyyy_0_xxxyyy_1, g_yyyy_0_xxyyyy_0, g_yyyy_0_xxyyyy_1, g_yyyy_0_xyyyyy_0, g_yyyy_0_xyyyyy_1, g_yyyy_0_yyyyyy_0, g_yyyy_0_yyyyyy_1, g_yyyyz_0_xxxxxy_1, g_yyyyz_0_xxxxyy_1, g_yyyyz_0_xxxyyy_1, g_yyyyz_0_xxyyyy_1, g_yyyyz_0_xyyyyy_1, g_yyyyz_0_yyyyyy_1, g_yyyyzz_0_xxxxxx_0, g_yyyyzz_0_xxxxxy_0, g_yyyyzz_0_xxxxxz_0, g_yyyyzz_0_xxxxyy_0, g_yyyyzz_0_xxxxyz_0, g_yyyyzz_0_xxxxzz_0, g_yyyyzz_0_xxxyyy_0, g_yyyyzz_0_xxxyyz_0, g_yyyyzz_0_xxxyzz_0, g_yyyyzz_0_xxxzzz_0, g_yyyyzz_0_xxyyyy_0, g_yyyyzz_0_xxyyyz_0, g_yyyyzz_0_xxyyzz_0, g_yyyyzz_0_xxyzzz_0, g_yyyyzz_0_xxzzzz_0, g_yyyyzz_0_xyyyyy_0, g_yyyyzz_0_xyyyyz_0, g_yyyyzz_0_xyyyzz_0, g_yyyyzz_0_xyyzzz_0, g_yyyyzz_0_xyzzzz_0, g_yyyyzz_0_xzzzzz_0, g_yyyyzz_0_yyyyyy_0, g_yyyyzz_0_yyyyyz_0, g_yyyyzz_0_yyyyzz_0, g_yyyyzz_0_yyyzzz_0, g_yyyyzz_0_yyzzzz_0, g_yyyyzz_0_yzzzzz_0, g_yyyyzz_0_zzzzzz_0, g_yyyzz_0_xxxxxx_1, g_yyyzz_0_xxxxxz_1, g_yyyzz_0_xxxxyz_1, g_yyyzz_0_xxxxz_1, g_yyyzz_0_xxxxzz_1, g_yyyzz_0_xxxyyz_1, g_yyyzz_0_xxxyz_1, g_yyyzz_0_xxxyzz_1, g_yyyzz_0_xxxzz_1, g_yyyzz_0_xxxzzz_1, g_yyyzz_0_xxyyyz_1, g_yyyzz_0_xxyyz_1, g_yyyzz_0_xxyyzz_1, g_yyyzz_0_xxyzz_1, g_yyyzz_0_xxyzzz_1, g_yyyzz_0_xxzzz_1, g_yyyzz_0_xxzzzz_1, g_yyyzz_0_xyyyyz_1, g_yyyzz_0_xyyyz_1, g_yyyzz_0_xyyyzz_1, g_yyyzz_0_xyyzz_1, g_yyyzz_0_xyyzzz_1, g_yyyzz_0_xyzzz_1, g_yyyzz_0_xyzzzz_1, g_yyyzz_0_xzzzz_1, g_yyyzz_0_xzzzzz_1, g_yyyzz_0_yyyyyz_1, g_yyyzz_0_yyyyz_1, g_yyyzz_0_yyyyzz_1, g_yyyzz_0_yyyzz_1, g_yyyzz_0_yyyzzz_1, g_yyyzz_0_yyzzz_1, g_yyyzz_0_yyzzzz_1, g_yyyzz_0_yzzzz_1, g_yyyzz_0_yzzzzz_1, g_yyyzz_0_zzzzz_1, g_yyyzz_0_zzzzzz_1, g_yyzz_0_xxxxxx_0, g_yyzz_0_xxxxxx_1, g_yyzz_0_xxxxxz_0, g_yyzz_0_xxxxxz_1, g_yyzz_0_xxxxyz_0, g_yyzz_0_xxxxyz_1, g_yyzz_0_xxxxzz_0, g_yyzz_0_xxxxzz_1, g_yyzz_0_xxxyyz_0, g_yyzz_0_xxxyyz_1, g_yyzz_0_xxxyzz_0, g_yyzz_0_xxxyzz_1, g_yyzz_0_xxxzzz_0, g_yyzz_0_xxxzzz_1, g_yyzz_0_xxyyyz_0, g_yyzz_0_xxyyyz_1, g_yyzz_0_xxyyzz_0, g_yyzz_0_xxyyzz_1, g_yyzz_0_xxyzzz_0, g_yyzz_0_xxyzzz_1, g_yyzz_0_xxzzzz_0, g_yyzz_0_xxzzzz_1, g_yyzz_0_xyyyyz_0, g_yyzz_0_xyyyyz_1, g_yyzz_0_xyyyzz_0, g_yyzz_0_xyyyzz_1, g_yyzz_0_xyyzzz_0, g_yyzz_0_xyyzzz_1, g_yyzz_0_xyzzzz_0, g_yyzz_0_xyzzzz_1, g_yyzz_0_xzzzzz_0, g_yyzz_0_xzzzzz_1, g_yyzz_0_yyyyyz_0, g_yyzz_0_yyyyyz_1, g_yyzz_0_yyyyzz_0, g_yyzz_0_yyyyzz_1, g_yyzz_0_yyyzzz_0, g_yyzz_0_yyyzzz_1, g_yyzz_0_yyzzzz_0, g_yyzz_0_yyzzzz_1, g_yyzz_0_yzzzzz_0, g_yyzz_0_yzzzzz_1, g_yyzz_0_zzzzzz_0, g_yyzz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzz_0_xxxxxx_0[i] = 3.0 * g_yyzz_0_xxxxxx_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxx_1[i] * fz_be_0 + g_yyyzz_0_xxxxxx_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxxy_0[i] = g_yyyy_0_xxxxxy_0[i] * fbe_0 - g_yyyy_0_xxxxxy_1[i] * fz_be_0 + g_yyyyz_0_xxxxxy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxxxz_0[i] = 3.0 * g_yyzz_0_xxxxxz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxz_1[i] * fz_be_0 + g_yyyzz_0_xxxxxz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxyy_0[i] = g_yyyy_0_xxxxyy_0[i] * fbe_0 - g_yyyy_0_xxxxyy_1[i] * fz_be_0 + g_yyyyz_0_xxxxyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxxyz_0[i] = 3.0 * g_yyzz_0_xxxxyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxyz_1[i] * fz_be_0 + g_yyyzz_0_xxxxz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxzz_0[i] = 3.0 * g_yyzz_0_xxxxzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxzz_1[i] * fz_be_0 + g_yyyzz_0_xxxxzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxyyy_0[i] = g_yyyy_0_xxxyyy_0[i] * fbe_0 - g_yyyy_0_xxxyyy_1[i] * fz_be_0 + g_yyyyz_0_xxxyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxyyz_0[i] = 3.0 * g_yyzz_0_xxxyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xxxyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxyzz_0[i] = 3.0 * g_yyzz_0_xxxyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyzz_1[i] * fz_be_0 + g_yyyzz_0_xxxzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxzzz_0[i] = 3.0 * g_yyzz_0_xxxzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxzzz_1[i] * fz_be_0 + g_yyyzz_0_xxxzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyyyy_0[i] = g_yyyy_0_xxyyyy_0[i] * fbe_0 - g_yyyy_0_xxyyyy_1[i] * fz_be_0 + g_yyyyz_0_xxyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxyyyz_0[i] = 3.0 * g_yyzz_0_xxyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_xxyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyyzz_0[i] = 3.0 * g_yyzz_0_xxyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xxyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyzzz_0[i] = 3.0 * g_yyzz_0_xxyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyzzz_1[i] * fz_be_0 + g_yyyzz_0_xxzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxzzzz_0[i] = 3.0 * g_yyzz_0_xxzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxzzzz_1[i] * fz_be_0 + g_yyyzz_0_xxzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyyyy_0[i] = g_yyyy_0_xyyyyy_0[i] * fbe_0 - g_yyyy_0_xyyyyy_1[i] * fz_be_0 + g_yyyyz_0_xyyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xyyyyz_0[i] = 3.0 * g_yyzz_0_xyyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyzz_0_xyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyyzz_0[i] = 3.0 * g_yyzz_0_xyyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_xyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyzzz_0[i] = 3.0 * g_yyzz_0_xyyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyzzzz_0[i] = 3.0 * g_yyzz_0_xyzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyzzzz_1[i] * fz_be_0 + g_yyyzz_0_xzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xzzzzz_0[i] = 3.0 * g_yyzz_0_xzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xzzzzz_1[i] * fz_be_0 + g_yyyzz_0_xzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyyyy_0[i] = g_yyyy_0_yyyyyy_0[i] * fbe_0 - g_yyyy_0_yyyyyy_1[i] * fz_be_0 + g_yyyyz_0_yyyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_yyyyyz_0[i] = 3.0 * g_yyzz_0_yyyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyzz_0_yyyyz_1[i] * fi_acd_0 + g_yyyzz_0_yyyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyyzz_0[i] = 3.0 * g_yyzz_0_yyyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyzz_0_yyyzz_1[i] * fi_acd_0 + g_yyyzz_0_yyyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyzzz_0[i] = 3.0 * g_yyzz_0_yyyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_yyzzz_1[i] * fi_acd_0 + g_yyyzz_0_yyyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyzzzz_0[i] = 3.0 * g_yyzz_0_yyzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_yzzzz_1[i] * fi_acd_0 + g_yyyzz_0_yyzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yzzzzz_0[i] = 3.0 * g_yyzz_0_yzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yzzzzz_1[i] * fz_be_0 + g_yyyzz_0_zzzzz_1[i] * fi_acd_0 + g_yyyzz_0_yzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_zzzzzz_0[i] = 3.0 * g_yyzz_0_zzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_zzzzzz_1[i] * fz_be_0 + g_yyyzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 672-700 components of targeted buffer : ISI

    auto g_yyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 672);

    auto g_yyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 673);

    auto g_yyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 674);

    auto g_yyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 675);

    auto g_yyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 676);

    auto g_yyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 677);

    auto g_yyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 678);

    auto g_yyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 679);

    auto g_yyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 680);

    auto g_yyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 681);

    auto g_yyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 682);

    auto g_yyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 683);

    auto g_yyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 684);

    auto g_yyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 685);

    auto g_yyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 686);

    auto g_yyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 687);

    auto g_yyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 688);

    auto g_yyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 689);

    auto g_yyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 690);

    auto g_yyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 691);

    auto g_yyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 692);

    auto g_yyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 693);

    auto g_yyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 694);

    auto g_yyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 695);

    auto g_yyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 696);

    auto g_yyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 697);

    auto g_yyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 698);

    auto g_yyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 699);

    #pragma omp simd aligned(g_yyyz_0_xxxxxy_0, g_yyyz_0_xxxxxy_1, g_yyyz_0_xxxxyy_0, g_yyyz_0_xxxxyy_1, g_yyyz_0_xxxyyy_0, g_yyyz_0_xxxyyy_1, g_yyyz_0_xxyyyy_0, g_yyyz_0_xxyyyy_1, g_yyyz_0_xyyyyy_0, g_yyyz_0_xyyyyy_1, g_yyyz_0_yyyyyy_0, g_yyyz_0_yyyyyy_1, g_yyyzz_0_xxxxxy_1, g_yyyzz_0_xxxxyy_1, g_yyyzz_0_xxxyyy_1, g_yyyzz_0_xxyyyy_1, g_yyyzz_0_xyyyyy_1, g_yyyzz_0_yyyyyy_1, g_yyyzzz_0_xxxxxx_0, g_yyyzzz_0_xxxxxy_0, g_yyyzzz_0_xxxxxz_0, g_yyyzzz_0_xxxxyy_0, g_yyyzzz_0_xxxxyz_0, g_yyyzzz_0_xxxxzz_0, g_yyyzzz_0_xxxyyy_0, g_yyyzzz_0_xxxyyz_0, g_yyyzzz_0_xxxyzz_0, g_yyyzzz_0_xxxzzz_0, g_yyyzzz_0_xxyyyy_0, g_yyyzzz_0_xxyyyz_0, g_yyyzzz_0_xxyyzz_0, g_yyyzzz_0_xxyzzz_0, g_yyyzzz_0_xxzzzz_0, g_yyyzzz_0_xyyyyy_0, g_yyyzzz_0_xyyyyz_0, g_yyyzzz_0_xyyyzz_0, g_yyyzzz_0_xyyzzz_0, g_yyyzzz_0_xyzzzz_0, g_yyyzzz_0_xzzzzz_0, g_yyyzzz_0_yyyyyy_0, g_yyyzzz_0_yyyyyz_0, g_yyyzzz_0_yyyyzz_0, g_yyyzzz_0_yyyzzz_0, g_yyyzzz_0_yyzzzz_0, g_yyyzzz_0_yzzzzz_0, g_yyyzzz_0_zzzzzz_0, g_yyzzz_0_xxxxxx_1, g_yyzzz_0_xxxxxz_1, g_yyzzz_0_xxxxyz_1, g_yyzzz_0_xxxxz_1, g_yyzzz_0_xxxxzz_1, g_yyzzz_0_xxxyyz_1, g_yyzzz_0_xxxyz_1, g_yyzzz_0_xxxyzz_1, g_yyzzz_0_xxxzz_1, g_yyzzz_0_xxxzzz_1, g_yyzzz_0_xxyyyz_1, g_yyzzz_0_xxyyz_1, g_yyzzz_0_xxyyzz_1, g_yyzzz_0_xxyzz_1, g_yyzzz_0_xxyzzz_1, g_yyzzz_0_xxzzz_1, g_yyzzz_0_xxzzzz_1, g_yyzzz_0_xyyyyz_1, g_yyzzz_0_xyyyz_1, g_yyzzz_0_xyyyzz_1, g_yyzzz_0_xyyzz_1, g_yyzzz_0_xyyzzz_1, g_yyzzz_0_xyzzz_1, g_yyzzz_0_xyzzzz_1, g_yyzzz_0_xzzzz_1, g_yyzzz_0_xzzzzz_1, g_yyzzz_0_yyyyyz_1, g_yyzzz_0_yyyyz_1, g_yyzzz_0_yyyyzz_1, g_yyzzz_0_yyyzz_1, g_yyzzz_0_yyyzzz_1, g_yyzzz_0_yyzzz_1, g_yyzzz_0_yyzzzz_1, g_yyzzz_0_yzzzz_1, g_yyzzz_0_yzzzzz_1, g_yyzzz_0_zzzzz_1, g_yyzzz_0_zzzzzz_1, g_yzzz_0_xxxxxx_0, g_yzzz_0_xxxxxx_1, g_yzzz_0_xxxxxz_0, g_yzzz_0_xxxxxz_1, g_yzzz_0_xxxxyz_0, g_yzzz_0_xxxxyz_1, g_yzzz_0_xxxxzz_0, g_yzzz_0_xxxxzz_1, g_yzzz_0_xxxyyz_0, g_yzzz_0_xxxyyz_1, g_yzzz_0_xxxyzz_0, g_yzzz_0_xxxyzz_1, g_yzzz_0_xxxzzz_0, g_yzzz_0_xxxzzz_1, g_yzzz_0_xxyyyz_0, g_yzzz_0_xxyyyz_1, g_yzzz_0_xxyyzz_0, g_yzzz_0_xxyyzz_1, g_yzzz_0_xxyzzz_0, g_yzzz_0_xxyzzz_1, g_yzzz_0_xxzzzz_0, g_yzzz_0_xxzzzz_1, g_yzzz_0_xyyyyz_0, g_yzzz_0_xyyyyz_1, g_yzzz_0_xyyyzz_0, g_yzzz_0_xyyyzz_1, g_yzzz_0_xyyzzz_0, g_yzzz_0_xyyzzz_1, g_yzzz_0_xyzzzz_0, g_yzzz_0_xyzzzz_1, g_yzzz_0_xzzzzz_0, g_yzzz_0_xzzzzz_1, g_yzzz_0_yyyyyz_0, g_yzzz_0_yyyyyz_1, g_yzzz_0_yyyyzz_0, g_yzzz_0_yyyyzz_1, g_yzzz_0_yyyzzz_0, g_yzzz_0_yyyzzz_1, g_yzzz_0_yyzzzz_0, g_yzzz_0_yyzzzz_1, g_yzzz_0_yzzzzz_0, g_yzzz_0_yzzzzz_1, g_yzzz_0_zzzzzz_0, g_yzzz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzz_0_xxxxxx_0[i] = 2.0 * g_yzzz_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxx_1[i] * fz_be_0 + g_yyzzz_0_xxxxxx_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxxy_0[i] = 2.0 * g_yyyz_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxxxy_1[i] * fz_be_0 + g_yyyzz_0_xxxxxy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxxxz_0[i] = 2.0 * g_yzzz_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxz_1[i] * fz_be_0 + g_yyzzz_0_xxxxxz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxyy_0[i] = 2.0 * g_yyyz_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxxyy_1[i] * fz_be_0 + g_yyyzz_0_xxxxyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxxyz_0[i] = 2.0 * g_yzzz_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxyz_1[i] * fz_be_0 + g_yyzzz_0_xxxxz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxzz_0[i] = 2.0 * g_yzzz_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxzz_1[i] * fz_be_0 + g_yyzzz_0_xxxxzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxyyy_0[i] = 2.0 * g_yyyz_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxyyy_1[i] * fz_be_0 + g_yyyzz_0_xxxyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxyyz_0[i] = 2.0 * g_yzzz_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xxxyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxyzz_0[i] = 2.0 * g_yzzz_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxyzz_1[i] * fz_be_0 + g_yyzzz_0_xxxzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxzzz_0[i] = 2.0 * g_yzzz_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxzzz_1[i] * fz_be_0 + g_yyzzz_0_xxxzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyyyy_0[i] = 2.0 * g_yyyz_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxyyyy_1[i] * fz_be_0 + g_yyyzz_0_xxyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxyyyz_0[i] = 2.0 * g_yzzz_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_xxyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyyzz_0[i] = 2.0 * g_yzzz_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xxyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyzzz_0[i] = 2.0 * g_yzzz_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyzzz_1[i] * fz_be_0 + g_yyzzz_0_xxzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxzzzz_0[i] = 2.0 * g_yzzz_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxzzzz_1[i] * fz_be_0 + g_yyzzz_0_xxzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyyyy_0[i] = 2.0 * g_yyyz_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xyyyyy_1[i] * fz_be_0 + g_yyyzz_0_xyyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xyyyyz_0[i] = 2.0 * g_yzzz_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyzzz_0_xyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyyzz_0[i] = 2.0 * g_yzzz_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_xyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyzzz_0[i] = 2.0 * g_yzzz_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyzzzz_0[i] = 2.0 * g_yzzz_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyzzzz_1[i] * fz_be_0 + g_yyzzz_0_xzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xzzzzz_0[i] = 2.0 * g_yzzz_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xzzzzz_1[i] * fz_be_0 + g_yyzzz_0_xzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyyyy_0[i] = 2.0 * g_yyyz_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_yyyyyy_1[i] * fz_be_0 + g_yyyzz_0_yyyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_yyyyyz_0[i] = 2.0 * g_yzzz_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyzzz_0_yyyyz_1[i] * fi_acd_0 + g_yyzzz_0_yyyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyyzz_0[i] = 2.0 * g_yzzz_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyzzz_0_yyyzz_1[i] * fi_acd_0 + g_yyzzz_0_yyyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyzzz_0[i] = 2.0 * g_yzzz_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_yyzzz_1[i] * fi_acd_0 + g_yyzzz_0_yyyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyzzzz_0[i] = 2.0 * g_yzzz_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_yzzzz_1[i] * fi_acd_0 + g_yyzzz_0_yyzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yzzzzz_0[i] = 2.0 * g_yzzz_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yzzzzz_1[i] * fz_be_0 + g_yyzzz_0_zzzzz_1[i] * fi_acd_0 + g_yyzzz_0_yzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_zzzzzz_0[i] = 2.0 * g_yzzz_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_zzzzzz_1[i] * fz_be_0 + g_yyzzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 700-728 components of targeted buffer : ISI

    auto g_yyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 700);

    auto g_yyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 701);

    auto g_yyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 702);

    auto g_yyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 703);

    auto g_yyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 704);

    auto g_yyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 705);

    auto g_yyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 706);

    auto g_yyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 707);

    auto g_yyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 708);

    auto g_yyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 709);

    auto g_yyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 710);

    auto g_yyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 711);

    auto g_yyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 712);

    auto g_yyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 713);

    auto g_yyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 714);

    auto g_yyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 715);

    auto g_yyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 716);

    auto g_yyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 717);

    auto g_yyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 718);

    auto g_yyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 719);

    auto g_yyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 720);

    auto g_yyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 721);

    auto g_yyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 722);

    auto g_yyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 723);

    auto g_yyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 724);

    auto g_yyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 725);

    auto g_yyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 726);

    auto g_yyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 727);

    #pragma omp simd aligned(g_yyzz_0_xxxxxy_0, g_yyzz_0_xxxxxy_1, g_yyzz_0_xxxxyy_0, g_yyzz_0_xxxxyy_1, g_yyzz_0_xxxyyy_0, g_yyzz_0_xxxyyy_1, g_yyzz_0_xxyyyy_0, g_yyzz_0_xxyyyy_1, g_yyzz_0_xyyyyy_0, g_yyzz_0_xyyyyy_1, g_yyzz_0_yyyyyy_0, g_yyzz_0_yyyyyy_1, g_yyzzz_0_xxxxxy_1, g_yyzzz_0_xxxxyy_1, g_yyzzz_0_xxxyyy_1, g_yyzzz_0_xxyyyy_1, g_yyzzz_0_xyyyyy_1, g_yyzzz_0_yyyyyy_1, g_yyzzzz_0_xxxxxx_0, g_yyzzzz_0_xxxxxy_0, g_yyzzzz_0_xxxxxz_0, g_yyzzzz_0_xxxxyy_0, g_yyzzzz_0_xxxxyz_0, g_yyzzzz_0_xxxxzz_0, g_yyzzzz_0_xxxyyy_0, g_yyzzzz_0_xxxyyz_0, g_yyzzzz_0_xxxyzz_0, g_yyzzzz_0_xxxzzz_0, g_yyzzzz_0_xxyyyy_0, g_yyzzzz_0_xxyyyz_0, g_yyzzzz_0_xxyyzz_0, g_yyzzzz_0_xxyzzz_0, g_yyzzzz_0_xxzzzz_0, g_yyzzzz_0_xyyyyy_0, g_yyzzzz_0_xyyyyz_0, g_yyzzzz_0_xyyyzz_0, g_yyzzzz_0_xyyzzz_0, g_yyzzzz_0_xyzzzz_0, g_yyzzzz_0_xzzzzz_0, g_yyzzzz_0_yyyyyy_0, g_yyzzzz_0_yyyyyz_0, g_yyzzzz_0_yyyyzz_0, g_yyzzzz_0_yyyzzz_0, g_yyzzzz_0_yyzzzz_0, g_yyzzzz_0_yzzzzz_0, g_yyzzzz_0_zzzzzz_0, g_yzzzz_0_xxxxxx_1, g_yzzzz_0_xxxxxz_1, g_yzzzz_0_xxxxyz_1, g_yzzzz_0_xxxxz_1, g_yzzzz_0_xxxxzz_1, g_yzzzz_0_xxxyyz_1, g_yzzzz_0_xxxyz_1, g_yzzzz_0_xxxyzz_1, g_yzzzz_0_xxxzz_1, g_yzzzz_0_xxxzzz_1, g_yzzzz_0_xxyyyz_1, g_yzzzz_0_xxyyz_1, g_yzzzz_0_xxyyzz_1, g_yzzzz_0_xxyzz_1, g_yzzzz_0_xxyzzz_1, g_yzzzz_0_xxzzz_1, g_yzzzz_0_xxzzzz_1, g_yzzzz_0_xyyyyz_1, g_yzzzz_0_xyyyz_1, g_yzzzz_0_xyyyzz_1, g_yzzzz_0_xyyzz_1, g_yzzzz_0_xyyzzz_1, g_yzzzz_0_xyzzz_1, g_yzzzz_0_xyzzzz_1, g_yzzzz_0_xzzzz_1, g_yzzzz_0_xzzzzz_1, g_yzzzz_0_yyyyyz_1, g_yzzzz_0_yyyyz_1, g_yzzzz_0_yyyyzz_1, g_yzzzz_0_yyyzz_1, g_yzzzz_0_yyyzzz_1, g_yzzzz_0_yyzzz_1, g_yzzzz_0_yyzzzz_1, g_yzzzz_0_yzzzz_1, g_yzzzz_0_yzzzzz_1, g_yzzzz_0_zzzzz_1, g_yzzzz_0_zzzzzz_1, g_zzzz_0_xxxxxx_0, g_zzzz_0_xxxxxx_1, g_zzzz_0_xxxxxz_0, g_zzzz_0_xxxxxz_1, g_zzzz_0_xxxxyz_0, g_zzzz_0_xxxxyz_1, g_zzzz_0_xxxxzz_0, g_zzzz_0_xxxxzz_1, g_zzzz_0_xxxyyz_0, g_zzzz_0_xxxyyz_1, g_zzzz_0_xxxyzz_0, g_zzzz_0_xxxyzz_1, g_zzzz_0_xxxzzz_0, g_zzzz_0_xxxzzz_1, g_zzzz_0_xxyyyz_0, g_zzzz_0_xxyyyz_1, g_zzzz_0_xxyyzz_0, g_zzzz_0_xxyyzz_1, g_zzzz_0_xxyzzz_0, g_zzzz_0_xxyzzz_1, g_zzzz_0_xxzzzz_0, g_zzzz_0_xxzzzz_1, g_zzzz_0_xyyyyz_0, g_zzzz_0_xyyyyz_1, g_zzzz_0_xyyyzz_0, g_zzzz_0_xyyyzz_1, g_zzzz_0_xyyzzz_0, g_zzzz_0_xyyzzz_1, g_zzzz_0_xyzzzz_0, g_zzzz_0_xyzzzz_1, g_zzzz_0_xzzzzz_0, g_zzzz_0_xzzzzz_1, g_zzzz_0_yyyyyz_0, g_zzzz_0_yyyyyz_1, g_zzzz_0_yyyyzz_0, g_zzzz_0_yyyyzz_1, g_zzzz_0_yyyzzz_0, g_zzzz_0_yyyzzz_1, g_zzzz_0_yyzzzz_0, g_zzzz_0_yyzzzz_1, g_zzzz_0_yzzzzz_0, g_zzzz_0_yzzzzz_1, g_zzzz_0_zzzzzz_0, g_zzzz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzz_0_xxxxxx_0[i] = g_zzzz_0_xxxxxx_0[i] * fbe_0 - g_zzzz_0_xxxxxx_1[i] * fz_be_0 + g_yzzzz_0_xxxxxx_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxxy_0[i] = 3.0 * g_yyzz_0_xxxxxy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxy_1[i] * fz_be_0 + g_yyzzz_0_xxxxxy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxxxz_0[i] = g_zzzz_0_xxxxxz_0[i] * fbe_0 - g_zzzz_0_xxxxxz_1[i] * fz_be_0 + g_yzzzz_0_xxxxxz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxyy_0[i] = 3.0 * g_yyzz_0_xxxxyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxyy_1[i] * fz_be_0 + g_yyzzz_0_xxxxyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxxyz_0[i] = g_zzzz_0_xxxxyz_0[i] * fbe_0 - g_zzzz_0_xxxxyz_1[i] * fz_be_0 + g_yzzzz_0_xxxxz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxzz_0[i] = g_zzzz_0_xxxxzz_0[i] * fbe_0 - g_zzzz_0_xxxxzz_1[i] * fz_be_0 + g_yzzzz_0_xxxxzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxyyy_0[i] = 3.0 * g_yyzz_0_xxxyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyyy_1[i] * fz_be_0 + g_yyzzz_0_xxxyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxyyz_0[i] = g_zzzz_0_xxxyyz_0[i] * fbe_0 - g_zzzz_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xxxyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxyzz_0[i] = g_zzzz_0_xxxyzz_0[i] * fbe_0 - g_zzzz_0_xxxyzz_1[i] * fz_be_0 + g_yzzzz_0_xxxzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxzzz_0[i] = g_zzzz_0_xxxzzz_0[i] * fbe_0 - g_zzzz_0_xxxzzz_1[i] * fz_be_0 + g_yzzzz_0_xxxzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyyyy_0[i] = 3.0 * g_yyzz_0_xxyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyyy_1[i] * fz_be_0 + g_yyzzz_0_xxyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxyyyz_0[i] = g_zzzz_0_xxyyyz_0[i] * fbe_0 - g_zzzz_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_xxyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyyzz_0[i] = g_zzzz_0_xxyyzz_0[i] * fbe_0 - g_zzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xxyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyzzz_0[i] = g_zzzz_0_xxyzzz_0[i] * fbe_0 - g_zzzz_0_xxyzzz_1[i] * fz_be_0 + g_yzzzz_0_xxzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxzzzz_0[i] = g_zzzz_0_xxzzzz_0[i] * fbe_0 - g_zzzz_0_xxzzzz_1[i] * fz_be_0 + g_yzzzz_0_xxzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyyyy_0[i] = 3.0 * g_yyzz_0_xyyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyyy_1[i] * fz_be_0 + g_yyzzz_0_xyyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xyyyyz_0[i] = g_zzzz_0_xyyyyz_0[i] * fbe_0 - g_zzzz_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yzzzz_0_xyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyyzz_0[i] = g_zzzz_0_xyyyzz_0[i] * fbe_0 - g_zzzz_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_xyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyzzz_0[i] = g_zzzz_0_xyyzzz_0[i] * fbe_0 - g_zzzz_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyzzzz_0[i] = g_zzzz_0_xyzzzz_0[i] * fbe_0 - g_zzzz_0_xyzzzz_1[i] * fz_be_0 + g_yzzzz_0_xzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xzzzzz_0[i] = g_zzzz_0_xzzzzz_0[i] * fbe_0 - g_zzzz_0_xzzzzz_1[i] * fz_be_0 + g_yzzzz_0_xzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyyyy_0[i] = 3.0 * g_yyzz_0_yyyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyyy_1[i] * fz_be_0 + g_yyzzz_0_yyyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_yyyyyz_0[i] = g_zzzz_0_yyyyyz_0[i] * fbe_0 - g_zzzz_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yzzzz_0_yyyyz_1[i] * fi_acd_0 + g_yzzzz_0_yyyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyyzz_0[i] = g_zzzz_0_yyyyzz_0[i] * fbe_0 - g_zzzz_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yzzzz_0_yyyzz_1[i] * fi_acd_0 + g_yzzzz_0_yyyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyzzz_0[i] = g_zzzz_0_yyyzzz_0[i] * fbe_0 - g_zzzz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_yyzzz_1[i] * fi_acd_0 + g_yzzzz_0_yyyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyzzzz_0[i] = g_zzzz_0_yyzzzz_0[i] * fbe_0 - g_zzzz_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_yzzzz_1[i] * fi_acd_0 + g_yzzzz_0_yyzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yzzzzz_0[i] = g_zzzz_0_yzzzzz_0[i] * fbe_0 - g_zzzz_0_yzzzzz_1[i] * fz_be_0 + g_yzzzz_0_zzzzz_1[i] * fi_acd_0 + g_yzzzz_0_yzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_zzzzzz_0[i] = g_zzzz_0_zzzzzz_0[i] * fbe_0 - g_zzzz_0_zzzzzz_1[i] * fz_be_0 + g_yzzzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 728-756 components of targeted buffer : ISI

    auto g_yzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 728);

    auto g_yzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 729);

    auto g_yzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 730);

    auto g_yzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 731);

    auto g_yzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 732);

    auto g_yzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 733);

    auto g_yzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 734);

    auto g_yzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 735);

    auto g_yzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 736);

    auto g_yzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 737);

    auto g_yzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 738);

    auto g_yzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 739);

    auto g_yzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 740);

    auto g_yzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 741);

    auto g_yzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 742);

    auto g_yzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 743);

    auto g_yzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 744);

    auto g_yzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 745);

    auto g_yzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 746);

    auto g_yzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 747);

    auto g_yzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 748);

    auto g_yzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 749);

    auto g_yzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 750);

    auto g_yzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 751);

    auto g_yzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 752);

    auto g_yzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 753);

    auto g_yzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 754);

    auto g_yzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 755);

    #pragma omp simd aligned(g_yzzzzz_0_xxxxxx_0, g_yzzzzz_0_xxxxxy_0, g_yzzzzz_0_xxxxxz_0, g_yzzzzz_0_xxxxyy_0, g_yzzzzz_0_xxxxyz_0, g_yzzzzz_0_xxxxzz_0, g_yzzzzz_0_xxxyyy_0, g_yzzzzz_0_xxxyyz_0, g_yzzzzz_0_xxxyzz_0, g_yzzzzz_0_xxxzzz_0, g_yzzzzz_0_xxyyyy_0, g_yzzzzz_0_xxyyyz_0, g_yzzzzz_0_xxyyzz_0, g_yzzzzz_0_xxyzzz_0, g_yzzzzz_0_xxzzzz_0, g_yzzzzz_0_xyyyyy_0, g_yzzzzz_0_xyyyyz_0, g_yzzzzz_0_xyyyzz_0, g_yzzzzz_0_xyyzzz_0, g_yzzzzz_0_xyzzzz_0, g_yzzzzz_0_xzzzzz_0, g_yzzzzz_0_yyyyyy_0, g_yzzzzz_0_yyyyyz_0, g_yzzzzz_0_yyyyzz_0, g_yzzzzz_0_yyyzzz_0, g_yzzzzz_0_yyzzzz_0, g_yzzzzz_0_yzzzzz_0, g_yzzzzz_0_zzzzzz_0, g_zzzzz_0_xxxxx_1, g_zzzzz_0_xxxxxx_1, g_zzzzz_0_xxxxxy_1, g_zzzzz_0_xxxxxz_1, g_zzzzz_0_xxxxy_1, g_zzzzz_0_xxxxyy_1, g_zzzzz_0_xxxxyz_1, g_zzzzz_0_xxxxz_1, g_zzzzz_0_xxxxzz_1, g_zzzzz_0_xxxyy_1, g_zzzzz_0_xxxyyy_1, g_zzzzz_0_xxxyyz_1, g_zzzzz_0_xxxyz_1, g_zzzzz_0_xxxyzz_1, g_zzzzz_0_xxxzz_1, g_zzzzz_0_xxxzzz_1, g_zzzzz_0_xxyyy_1, g_zzzzz_0_xxyyyy_1, g_zzzzz_0_xxyyyz_1, g_zzzzz_0_xxyyz_1, g_zzzzz_0_xxyyzz_1, g_zzzzz_0_xxyzz_1, g_zzzzz_0_xxyzzz_1, g_zzzzz_0_xxzzz_1, g_zzzzz_0_xxzzzz_1, g_zzzzz_0_xyyyy_1, g_zzzzz_0_xyyyyy_1, g_zzzzz_0_xyyyyz_1, g_zzzzz_0_xyyyz_1, g_zzzzz_0_xyyyzz_1, g_zzzzz_0_xyyzz_1, g_zzzzz_0_xyyzzz_1, g_zzzzz_0_xyzzz_1, g_zzzzz_0_xyzzzz_1, g_zzzzz_0_xzzzz_1, g_zzzzz_0_xzzzzz_1, g_zzzzz_0_yyyyy_1, g_zzzzz_0_yyyyyy_1, g_zzzzz_0_yyyyyz_1, g_zzzzz_0_yyyyz_1, g_zzzzz_0_yyyyzz_1, g_zzzzz_0_yyyzz_1, g_zzzzz_0_yyyzzz_1, g_zzzzz_0_yyzzz_1, g_zzzzz_0_yyzzzz_1, g_zzzzz_0_yzzzz_1, g_zzzzz_0_yzzzzz_1, g_zzzzz_0_zzzzz_1, g_zzzzz_0_zzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzz_0_xxxxxx_0[i] = g_zzzzz_0_xxxxxx_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxy_0[i] = g_zzzzz_0_xxxxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxz_0[i] = g_zzzzz_0_xxxxxz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxyy_0[i] = 2.0 * g_zzzzz_0_xxxxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxyz_0[i] = g_zzzzz_0_xxxxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxzz_0[i] = g_zzzzz_0_xxxxzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyyy_0[i] = 3.0 * g_zzzzz_0_xxxyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyyz_0[i] = 2.0 * g_zzzzz_0_xxxyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyzz_0[i] = g_zzzzz_0_xxxzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxzzz_0[i] = g_zzzzz_0_xxxzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyyy_0[i] = 4.0 * g_zzzzz_0_xxyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyyz_0[i] = 3.0 * g_zzzzz_0_xxyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyzz_0[i] = 2.0 * g_zzzzz_0_xxyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyzzz_0[i] = g_zzzzz_0_xxzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxzzzz_0[i] = g_zzzzz_0_xxzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyyy_0[i] = 5.0 * g_zzzzz_0_xyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyyz_0[i] = 4.0 * g_zzzzz_0_xyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyzz_0[i] = 3.0 * g_zzzzz_0_xyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyzzz_0[i] = 2.0 * g_zzzzz_0_xyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyzzzz_0[i] = g_zzzzz_0_xzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xzzzzz_0[i] = g_zzzzz_0_xzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyyy_0[i] = 6.0 * g_zzzzz_0_yyyyy_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyyz_0[i] = 5.0 * g_zzzzz_0_yyyyz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyzz_0[i] = 4.0 * g_zzzzz_0_yyyzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyzzz_0[i] = 3.0 * g_zzzzz_0_yyzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyzzzz_0[i] = 2.0 * g_zzzzz_0_yzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yzzzzz_0[i] = g_zzzzz_0_zzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_zzzzzz_0[i] = g_zzzzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 756-784 components of targeted buffer : ISI

    auto g_zzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_isi + 756);

    auto g_zzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_isi + 757);

    auto g_zzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_isi + 758);

    auto g_zzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_isi + 759);

    auto g_zzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_isi + 760);

    auto g_zzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_isi + 761);

    auto g_zzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_isi + 762);

    auto g_zzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_isi + 763);

    auto g_zzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_isi + 764);

    auto g_zzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_isi + 765);

    auto g_zzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_isi + 766);

    auto g_zzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_isi + 767);

    auto g_zzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_isi + 768);

    auto g_zzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_isi + 769);

    auto g_zzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_isi + 770);

    auto g_zzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_isi + 771);

    auto g_zzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_isi + 772);

    auto g_zzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_isi + 773);

    auto g_zzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_isi + 774);

    auto g_zzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_isi + 775);

    auto g_zzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_isi + 776);

    auto g_zzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_isi + 777);

    auto g_zzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_isi + 778);

    auto g_zzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_isi + 779);

    auto g_zzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_isi + 780);

    auto g_zzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_isi + 781);

    auto g_zzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_isi + 782);

    auto g_zzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_isi + 783);

    #pragma omp simd aligned(g_zzzz_0_xxxxxx_0, g_zzzz_0_xxxxxx_1, g_zzzz_0_xxxxxy_0, g_zzzz_0_xxxxxy_1, g_zzzz_0_xxxxxz_0, g_zzzz_0_xxxxxz_1, g_zzzz_0_xxxxyy_0, g_zzzz_0_xxxxyy_1, g_zzzz_0_xxxxyz_0, g_zzzz_0_xxxxyz_1, g_zzzz_0_xxxxzz_0, g_zzzz_0_xxxxzz_1, g_zzzz_0_xxxyyy_0, g_zzzz_0_xxxyyy_1, g_zzzz_0_xxxyyz_0, g_zzzz_0_xxxyyz_1, g_zzzz_0_xxxyzz_0, g_zzzz_0_xxxyzz_1, g_zzzz_0_xxxzzz_0, g_zzzz_0_xxxzzz_1, g_zzzz_0_xxyyyy_0, g_zzzz_0_xxyyyy_1, g_zzzz_0_xxyyyz_0, g_zzzz_0_xxyyyz_1, g_zzzz_0_xxyyzz_0, g_zzzz_0_xxyyzz_1, g_zzzz_0_xxyzzz_0, g_zzzz_0_xxyzzz_1, g_zzzz_0_xxzzzz_0, g_zzzz_0_xxzzzz_1, g_zzzz_0_xyyyyy_0, g_zzzz_0_xyyyyy_1, g_zzzz_0_xyyyyz_0, g_zzzz_0_xyyyyz_1, g_zzzz_0_xyyyzz_0, g_zzzz_0_xyyyzz_1, g_zzzz_0_xyyzzz_0, g_zzzz_0_xyyzzz_1, g_zzzz_0_xyzzzz_0, g_zzzz_0_xyzzzz_1, g_zzzz_0_xzzzzz_0, g_zzzz_0_xzzzzz_1, g_zzzz_0_yyyyyy_0, g_zzzz_0_yyyyyy_1, g_zzzz_0_yyyyyz_0, g_zzzz_0_yyyyyz_1, g_zzzz_0_yyyyzz_0, g_zzzz_0_yyyyzz_1, g_zzzz_0_yyyzzz_0, g_zzzz_0_yyyzzz_1, g_zzzz_0_yyzzzz_0, g_zzzz_0_yyzzzz_1, g_zzzz_0_yzzzzz_0, g_zzzz_0_yzzzzz_1, g_zzzz_0_zzzzzz_0, g_zzzz_0_zzzzzz_1, g_zzzzz_0_xxxxx_1, g_zzzzz_0_xxxxxx_1, g_zzzzz_0_xxxxxy_1, g_zzzzz_0_xxxxxz_1, g_zzzzz_0_xxxxy_1, g_zzzzz_0_xxxxyy_1, g_zzzzz_0_xxxxyz_1, g_zzzzz_0_xxxxz_1, g_zzzzz_0_xxxxzz_1, g_zzzzz_0_xxxyy_1, g_zzzzz_0_xxxyyy_1, g_zzzzz_0_xxxyyz_1, g_zzzzz_0_xxxyz_1, g_zzzzz_0_xxxyzz_1, g_zzzzz_0_xxxzz_1, g_zzzzz_0_xxxzzz_1, g_zzzzz_0_xxyyy_1, g_zzzzz_0_xxyyyy_1, g_zzzzz_0_xxyyyz_1, g_zzzzz_0_xxyyz_1, g_zzzzz_0_xxyyzz_1, g_zzzzz_0_xxyzz_1, g_zzzzz_0_xxyzzz_1, g_zzzzz_0_xxzzz_1, g_zzzzz_0_xxzzzz_1, g_zzzzz_0_xyyyy_1, g_zzzzz_0_xyyyyy_1, g_zzzzz_0_xyyyyz_1, g_zzzzz_0_xyyyz_1, g_zzzzz_0_xyyyzz_1, g_zzzzz_0_xyyzz_1, g_zzzzz_0_xyyzzz_1, g_zzzzz_0_xyzzz_1, g_zzzzz_0_xyzzzz_1, g_zzzzz_0_xzzzz_1, g_zzzzz_0_xzzzzz_1, g_zzzzz_0_yyyyy_1, g_zzzzz_0_yyyyyy_1, g_zzzzz_0_yyyyyz_1, g_zzzzz_0_yyyyz_1, g_zzzzz_0_yyyyzz_1, g_zzzzz_0_yyyzz_1, g_zzzzz_0_yyyzzz_1, g_zzzzz_0_yyzzz_1, g_zzzzz_0_yyzzzz_1, g_zzzzz_0_yzzzz_1, g_zzzzz_0_yzzzzz_1, g_zzzzz_0_zzzzz_1, g_zzzzz_0_zzzzzz_1, g_zzzzzz_0_xxxxxx_0, g_zzzzzz_0_xxxxxy_0, g_zzzzzz_0_xxxxxz_0, g_zzzzzz_0_xxxxyy_0, g_zzzzzz_0_xxxxyz_0, g_zzzzzz_0_xxxxzz_0, g_zzzzzz_0_xxxyyy_0, g_zzzzzz_0_xxxyyz_0, g_zzzzzz_0_xxxyzz_0, g_zzzzzz_0_xxxzzz_0, g_zzzzzz_0_xxyyyy_0, g_zzzzzz_0_xxyyyz_0, g_zzzzzz_0_xxyyzz_0, g_zzzzzz_0_xxyzzz_0, g_zzzzzz_0_xxzzzz_0, g_zzzzzz_0_xyyyyy_0, g_zzzzzz_0_xyyyyz_0, g_zzzzzz_0_xyyyzz_0, g_zzzzzz_0_xyyzzz_0, g_zzzzzz_0_xyzzzz_0, g_zzzzzz_0_xzzzzz_0, g_zzzzzz_0_yyyyyy_0, g_zzzzzz_0_yyyyyz_0, g_zzzzzz_0_yyyyzz_0, g_zzzzzz_0_yyyzzz_0, g_zzzzzz_0_yyzzzz_0, g_zzzzzz_0_yzzzzz_0, g_zzzzzz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzz_0_xxxxxx_0[i] = 5.0 * g_zzzz_0_xxxxxx_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxx_1[i] * fz_be_0 + g_zzzzz_0_xxxxxx_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxy_0[i] = 5.0 * g_zzzz_0_xxxxxy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxy_1[i] * fz_be_0 + g_zzzzz_0_xxxxxy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxz_0[i] = 5.0 * g_zzzz_0_xxxxxz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxz_1[i] * fz_be_0 + g_zzzzz_0_xxxxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxyy_0[i] = 5.0 * g_zzzz_0_xxxxyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxyy_1[i] * fz_be_0 + g_zzzzz_0_xxxxyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxyz_0[i] = 5.0 * g_zzzz_0_xxxxyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxyz_1[i] * fz_be_0 + g_zzzzz_0_xxxxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxzz_0[i] = 5.0 * g_zzzz_0_xxxxzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxxxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyyy_0[i] = 5.0 * g_zzzz_0_xxxyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyyy_1[i] * fz_be_0 + g_zzzzz_0_xxxyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyyz_0[i] = 5.0 * g_zzzz_0_xxxyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyyz_1[i] * fz_be_0 + g_zzzzz_0_xxxyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyzz_0[i] = 5.0 * g_zzzz_0_xxxyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxxyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxzzz_0[i] = 5.0 * g_zzzz_0_xxxzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xxxzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyyy_0[i] = 5.0 * g_zzzz_0_xxyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyyy_1[i] * fz_be_0 + g_zzzzz_0_xxyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyyz_0[i] = 5.0 * g_zzzz_0_xxyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyyz_1[i] * fz_be_0 + g_zzzzz_0_xxyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyzz_0[i] = 5.0 * g_zzzz_0_xxyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyzzz_0[i] = 5.0 * g_zzzz_0_xxyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xxyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxzzzz_0[i] = 5.0 * g_zzzz_0_xxzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_xxzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyyy_0[i] = 5.0 * g_zzzz_0_xyyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyyy_1[i] * fz_be_0 + g_zzzzz_0_xyyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyyz_0[i] = 5.0 * g_zzzz_0_xyyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyyz_1[i] * fz_be_0 + g_zzzzz_0_xyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyzz_0[i] = 5.0 * g_zzzz_0_xyyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyzzz_0[i] = 5.0 * g_zzzz_0_xyyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyzzzz_0[i] = 5.0 * g_zzzz_0_xyzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_xyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xzzzzz_0[i] = 5.0 * g_zzzz_0_xzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_0_xzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyyy_0[i] = 5.0 * g_zzzz_0_yyyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyyy_1[i] * fz_be_0 + g_zzzzz_0_yyyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyyz_0[i] = 5.0 * g_zzzz_0_yyyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyyz_1[i] * fz_be_0 + g_zzzzz_0_yyyyy_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyzz_0[i] = 5.0 * g_zzzz_0_yyyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_yyyyz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyzzz_0[i] = 5.0 * g_zzzz_0_yyyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_yyyzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyzzzz_0[i] = 5.0 * g_zzzz_0_yyzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_yyzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yzzzzz_0[i] = 5.0 * g_zzzz_0_yzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_0_yzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_zzzzzz_0[i] = 5.0 * g_zzzz_0_zzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_zzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzz_0_zzzzz_1[i] * fi_acd_0 + g_zzzzz_0_zzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

