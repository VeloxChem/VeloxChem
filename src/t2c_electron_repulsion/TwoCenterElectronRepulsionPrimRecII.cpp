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

#include "TwoCenterElectronRepulsionPrimRecII.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ii(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ii,
                                const size_t idx_eri_0_gi,
                                const size_t idx_eri_1_gi,
                                const size_t idx_eri_1_hh,
                                const size_t idx_eri_1_hi,
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

    // Set up components of auxiliary buffer : GI

    auto g_xxxx_xxxxxx_0 = pbuffer.data(idx_eri_0_gi);

    auto g_xxxx_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 1);

    auto g_xxxx_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 2);

    auto g_xxxx_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 3);

    auto g_xxxx_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 4);

    auto g_xxxx_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 5);

    auto g_xxxx_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 6);

    auto g_xxxx_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 7);

    auto g_xxxx_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 8);

    auto g_xxxx_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 9);

    auto g_xxxx_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 10);

    auto g_xxxx_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 11);

    auto g_xxxx_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 12);

    auto g_xxxx_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 13);

    auto g_xxxx_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 14);

    auto g_xxxx_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 15);

    auto g_xxxx_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 16);

    auto g_xxxx_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 17);

    auto g_xxxx_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 18);

    auto g_xxxx_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 19);

    auto g_xxxx_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 20);

    auto g_xxxx_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 21);

    auto g_xxxx_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 22);

    auto g_xxxx_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 23);

    auto g_xxxx_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 24);

    auto g_xxxx_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 25);

    auto g_xxxx_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 26);

    auto g_xxxx_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 27);

    auto g_xxxy_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 28);

    auto g_xxxy_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 30);

    auto g_xxxy_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 33);

    auto g_xxxy_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 37);

    auto g_xxxy_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 42);

    auto g_xxxy_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 48);

    auto g_xxxz_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 56);

    auto g_xxxz_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 57);

    auto g_xxxz_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 59);

    auto g_xxxz_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 62);

    auto g_xxxz_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 66);

    auto g_xxxz_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 71);

    auto g_xxyy_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 84);

    auto g_xxyy_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 85);

    auto g_xxyy_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 86);

    auto g_xxyy_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 87);

    auto g_xxyy_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 88);

    auto g_xxyy_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 89);

    auto g_xxyy_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 90);

    auto g_xxyy_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 91);

    auto g_xxyy_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 92);

    auto g_xxyy_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 93);

    auto g_xxyy_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 94);

    auto g_xxyy_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 95);

    auto g_xxyy_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 96);

    auto g_xxyy_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 97);

    auto g_xxyy_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 98);

    auto g_xxyy_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 99);

    auto g_xxyy_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 100);

    auto g_xxyy_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 101);

    auto g_xxyy_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 102);

    auto g_xxyy_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 103);

    auto g_xxyy_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 104);

    auto g_xxyy_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 105);

    auto g_xxyy_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 106);

    auto g_xxyy_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 107);

    auto g_xxyy_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 108);

    auto g_xxyy_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 109);

    auto g_xxyy_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 110);

    auto g_xxyy_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 111);

    auto g_xxzz_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 140);

    auto g_xxzz_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 141);

    auto g_xxzz_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 142);

    auto g_xxzz_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 143);

    auto g_xxzz_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 144);

    auto g_xxzz_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 145);

    auto g_xxzz_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 146);

    auto g_xxzz_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 147);

    auto g_xxzz_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 148);

    auto g_xxzz_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 149);

    auto g_xxzz_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 150);

    auto g_xxzz_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 151);

    auto g_xxzz_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 152);

    auto g_xxzz_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 153);

    auto g_xxzz_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 154);

    auto g_xxzz_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 155);

    auto g_xxzz_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 156);

    auto g_xxzz_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 157);

    auto g_xxzz_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 158);

    auto g_xxzz_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 159);

    auto g_xxzz_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 160);

    auto g_xxzz_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 161);

    auto g_xxzz_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 162);

    auto g_xxzz_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 163);

    auto g_xxzz_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 164);

    auto g_xxzz_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 165);

    auto g_xxzz_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 166);

    auto g_xxzz_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 167);

    auto g_xyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 169);

    auto g_xyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 171);

    auto g_xyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 172);

    auto g_xyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 174);

    auto g_xyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 175);

    auto g_xyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 176);

    auto g_xyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 178);

    auto g_xyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 179);

    auto g_xyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 180);

    auto g_xyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 181);

    auto g_xyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 183);

    auto g_xyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 184);

    auto g_xyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 185);

    auto g_xyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 186);

    auto g_xyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 187);

    auto g_xyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 189);

    auto g_xyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 190);

    auto g_xyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 191);

    auto g_xyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 192);

    auto g_xyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 193);

    auto g_xyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 194);

    auto g_xyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 195);

    auto g_xzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 254);

    auto g_xzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 256);

    auto g_xzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 257);

    auto g_xzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 259);

    auto g_xzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 260);

    auto g_xzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 261);

    auto g_xzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 263);

    auto g_xzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 264);

    auto g_xzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 265);

    auto g_xzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 266);

    auto g_xzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 268);

    auto g_xzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 269);

    auto g_xzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 270);

    auto g_xzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 271);

    auto g_xzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 272);

    auto g_xzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 273);

    auto g_xzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 274);

    auto g_xzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 275);

    auto g_xzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 276);

    auto g_xzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 277);

    auto g_xzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 278);

    auto g_xzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 279);

    auto g_yyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 280);

    auto g_yyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 281);

    auto g_yyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 282);

    auto g_yyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 283);

    auto g_yyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 284);

    auto g_yyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 285);

    auto g_yyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 286);

    auto g_yyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 287);

    auto g_yyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 288);

    auto g_yyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 289);

    auto g_yyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 290);

    auto g_yyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 291);

    auto g_yyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 292);

    auto g_yyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 293);

    auto g_yyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 294);

    auto g_yyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 295);

    auto g_yyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 296);

    auto g_yyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 297);

    auto g_yyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 298);

    auto g_yyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 299);

    auto g_yyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 300);

    auto g_yyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 301);

    auto g_yyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 302);

    auto g_yyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 303);

    auto g_yyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 304);

    auto g_yyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 305);

    auto g_yyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 306);

    auto g_yyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 307);

    auto g_yyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 309);

    auto g_yyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 311);

    auto g_yyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 314);

    auto g_yyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 318);

    auto g_yyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 323);

    auto g_yyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 329);

    auto g_yyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 336);

    auto g_yyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 337);

    auto g_yyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 338);

    auto g_yyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 339);

    auto g_yyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 340);

    auto g_yyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 341);

    auto g_yyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 342);

    auto g_yyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 343);

    auto g_yyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 344);

    auto g_yyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 345);

    auto g_yyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 346);

    auto g_yyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 347);

    auto g_yyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 348);

    auto g_yyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 349);

    auto g_yyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 350);

    auto g_yyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 351);

    auto g_yyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 352);

    auto g_yyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 353);

    auto g_yyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 354);

    auto g_yyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 355);

    auto g_yyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 356);

    auto g_yyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 357);

    auto g_yyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 358);

    auto g_yyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 359);

    auto g_yyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 360);

    auto g_yyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 361);

    auto g_yyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 362);

    auto g_yyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 363);

    auto g_yzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 364);

    auto g_yzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 366);

    auto g_yzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 368);

    auto g_yzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 369);

    auto g_yzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 371);

    auto g_yzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 372);

    auto g_yzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 373);

    auto g_yzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 375);

    auto g_yzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 376);

    auto g_yzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 377);

    auto g_yzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 378);

    auto g_yzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 380);

    auto g_yzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 381);

    auto g_yzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 382);

    auto g_yzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 383);

    auto g_yzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 384);

    auto g_yzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 386);

    auto g_yzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 387);

    auto g_yzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 388);

    auto g_yzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 389);

    auto g_yzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 390);

    auto g_yzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 391);

    auto g_zzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 392);

    auto g_zzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 393);

    auto g_zzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 394);

    auto g_zzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 395);

    auto g_zzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 396);

    auto g_zzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 397);

    auto g_zzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 398);

    auto g_zzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 399);

    auto g_zzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 400);

    auto g_zzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 401);

    auto g_zzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 402);

    auto g_zzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 403);

    auto g_zzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 404);

    auto g_zzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 405);

    auto g_zzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 406);

    auto g_zzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 407);

    auto g_zzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 408);

    auto g_zzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 409);

    auto g_zzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 410);

    auto g_zzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 411);

    auto g_zzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 412);

    auto g_zzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 413);

    auto g_zzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 414);

    auto g_zzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 415);

    auto g_zzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 416);

    auto g_zzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 417);

    auto g_zzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 418);

    auto g_zzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 419);

    // Set up components of auxiliary buffer : GI

    auto g_xxxx_xxxxxx_1 = pbuffer.data(idx_eri_1_gi);

    auto g_xxxx_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 1);

    auto g_xxxx_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 2);

    auto g_xxxx_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 3);

    auto g_xxxx_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 4);

    auto g_xxxx_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 5);

    auto g_xxxx_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 6);

    auto g_xxxx_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 7);

    auto g_xxxx_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 8);

    auto g_xxxx_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 9);

    auto g_xxxx_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 10);

    auto g_xxxx_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 11);

    auto g_xxxx_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 12);

    auto g_xxxx_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 13);

    auto g_xxxx_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 14);

    auto g_xxxx_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 15);

    auto g_xxxx_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 16);

    auto g_xxxx_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 17);

    auto g_xxxx_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 18);

    auto g_xxxx_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 19);

    auto g_xxxx_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 20);

    auto g_xxxx_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 21);

    auto g_xxxx_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 22);

    auto g_xxxx_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 23);

    auto g_xxxx_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 24);

    auto g_xxxx_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 25);

    auto g_xxxx_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 26);

    auto g_xxxx_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 27);

    auto g_xxxy_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 28);

    auto g_xxxy_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 30);

    auto g_xxxy_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 33);

    auto g_xxxy_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 37);

    auto g_xxxy_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 42);

    auto g_xxxy_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 48);

    auto g_xxxz_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 56);

    auto g_xxxz_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 57);

    auto g_xxxz_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 59);

    auto g_xxxz_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 62);

    auto g_xxxz_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 66);

    auto g_xxxz_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 71);

    auto g_xxyy_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 84);

    auto g_xxyy_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 85);

    auto g_xxyy_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 86);

    auto g_xxyy_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 87);

    auto g_xxyy_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 88);

    auto g_xxyy_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 89);

    auto g_xxyy_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 90);

    auto g_xxyy_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 91);

    auto g_xxyy_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 92);

    auto g_xxyy_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 93);

    auto g_xxyy_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 94);

    auto g_xxyy_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 95);

    auto g_xxyy_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 96);

    auto g_xxyy_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 97);

    auto g_xxyy_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 98);

    auto g_xxyy_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 99);

    auto g_xxyy_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 100);

    auto g_xxyy_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 101);

    auto g_xxyy_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 102);

    auto g_xxyy_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 103);

    auto g_xxyy_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 104);

    auto g_xxyy_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 105);

    auto g_xxyy_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 106);

    auto g_xxyy_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 107);

    auto g_xxyy_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 108);

    auto g_xxyy_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 109);

    auto g_xxyy_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 110);

    auto g_xxyy_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 111);

    auto g_xxzz_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 140);

    auto g_xxzz_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 141);

    auto g_xxzz_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 142);

    auto g_xxzz_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 143);

    auto g_xxzz_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 144);

    auto g_xxzz_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 145);

    auto g_xxzz_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 146);

    auto g_xxzz_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 147);

    auto g_xxzz_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 148);

    auto g_xxzz_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 149);

    auto g_xxzz_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 150);

    auto g_xxzz_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 151);

    auto g_xxzz_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 152);

    auto g_xxzz_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 153);

    auto g_xxzz_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 154);

    auto g_xxzz_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 155);

    auto g_xxzz_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 156);

    auto g_xxzz_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 157);

    auto g_xxzz_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 158);

    auto g_xxzz_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 159);

    auto g_xxzz_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 160);

    auto g_xxzz_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 161);

    auto g_xxzz_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 162);

    auto g_xxzz_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 163);

    auto g_xxzz_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 164);

    auto g_xxzz_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 165);

    auto g_xxzz_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 166);

    auto g_xxzz_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 167);

    auto g_xyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 169);

    auto g_xyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 171);

    auto g_xyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 172);

    auto g_xyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 174);

    auto g_xyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 175);

    auto g_xyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 176);

    auto g_xyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 178);

    auto g_xyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 179);

    auto g_xyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 180);

    auto g_xyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 181);

    auto g_xyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 183);

    auto g_xyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 184);

    auto g_xyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 185);

    auto g_xyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 186);

    auto g_xyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 187);

    auto g_xyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 189);

    auto g_xyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 190);

    auto g_xyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 191);

    auto g_xyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 192);

    auto g_xyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 193);

    auto g_xyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 194);

    auto g_xyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 195);

    auto g_xzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 254);

    auto g_xzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 256);

    auto g_xzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 257);

    auto g_xzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 259);

    auto g_xzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 260);

    auto g_xzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 261);

    auto g_xzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 263);

    auto g_xzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 264);

    auto g_xzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 265);

    auto g_xzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 266);

    auto g_xzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 268);

    auto g_xzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 269);

    auto g_xzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 270);

    auto g_xzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 271);

    auto g_xzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 272);

    auto g_xzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 273);

    auto g_xzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 274);

    auto g_xzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 275);

    auto g_xzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 276);

    auto g_xzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 277);

    auto g_xzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 278);

    auto g_xzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 279);

    auto g_yyyy_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 280);

    auto g_yyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 281);

    auto g_yyyy_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 282);

    auto g_yyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 283);

    auto g_yyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 284);

    auto g_yyyy_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 285);

    auto g_yyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 286);

    auto g_yyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 287);

    auto g_yyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 288);

    auto g_yyyy_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 289);

    auto g_yyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 290);

    auto g_yyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 291);

    auto g_yyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 292);

    auto g_yyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 293);

    auto g_yyyy_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 294);

    auto g_yyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 295);

    auto g_yyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 296);

    auto g_yyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 297);

    auto g_yyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 298);

    auto g_yyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 299);

    auto g_yyyy_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 300);

    auto g_yyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 301);

    auto g_yyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 302);

    auto g_yyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 303);

    auto g_yyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 304);

    auto g_yyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 305);

    auto g_yyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 306);

    auto g_yyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 307);

    auto g_yyyz_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 309);

    auto g_yyyz_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 311);

    auto g_yyyz_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 314);

    auto g_yyyz_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 318);

    auto g_yyyz_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 323);

    auto g_yyyz_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 329);

    auto g_yyzz_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 336);

    auto g_yyzz_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 337);

    auto g_yyzz_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 338);

    auto g_yyzz_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 339);

    auto g_yyzz_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 340);

    auto g_yyzz_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 341);

    auto g_yyzz_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 342);

    auto g_yyzz_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 343);

    auto g_yyzz_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 344);

    auto g_yyzz_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 345);

    auto g_yyzz_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 346);

    auto g_yyzz_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 347);

    auto g_yyzz_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 348);

    auto g_yyzz_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 349);

    auto g_yyzz_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 350);

    auto g_yyzz_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 351);

    auto g_yyzz_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 352);

    auto g_yyzz_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 353);

    auto g_yyzz_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 354);

    auto g_yyzz_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 355);

    auto g_yyzz_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 356);

    auto g_yyzz_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 357);

    auto g_yyzz_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 358);

    auto g_yyzz_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 359);

    auto g_yyzz_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 360);

    auto g_yyzz_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 361);

    auto g_yyzz_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 362);

    auto g_yyzz_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 363);

    auto g_yzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 364);

    auto g_yzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 366);

    auto g_yzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 368);

    auto g_yzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 369);

    auto g_yzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 371);

    auto g_yzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 372);

    auto g_yzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 373);

    auto g_yzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 375);

    auto g_yzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 376);

    auto g_yzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 377);

    auto g_yzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 378);

    auto g_yzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 380);

    auto g_yzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 381);

    auto g_yzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 382);

    auto g_yzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 383);

    auto g_yzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 384);

    auto g_yzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 386);

    auto g_yzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 387);

    auto g_yzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 388);

    auto g_yzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 389);

    auto g_yzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 390);

    auto g_yzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 391);

    auto g_zzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 392);

    auto g_zzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 393);

    auto g_zzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 394);

    auto g_zzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 395);

    auto g_zzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 396);

    auto g_zzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 397);

    auto g_zzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 398);

    auto g_zzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 399);

    auto g_zzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 400);

    auto g_zzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 401);

    auto g_zzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 402);

    auto g_zzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 403);

    auto g_zzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 404);

    auto g_zzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 405);

    auto g_zzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 406);

    auto g_zzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 407);

    auto g_zzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 408);

    auto g_zzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 409);

    auto g_zzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 410);

    auto g_zzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 411);

    auto g_zzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 412);

    auto g_zzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 413);

    auto g_zzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 414);

    auto g_zzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 415);

    auto g_zzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 416);

    auto g_zzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 417);

    auto g_zzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 418);

    auto g_zzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 419);

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

    auto g_xxxxz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 44);

    auto g_xxxxz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 46);

    auto g_xxxxz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 47);

    auto g_xxxxz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 49);

    auto g_xxxxz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 50);

    auto g_xxxxz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 51);

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

    auto g_xyyzz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 256);

    auto g_xyyzz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 259);

    auto g_xyyzz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 260);

    auto g_xyyzz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 263);

    auto g_xyyzz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 264);

    auto g_xyyzz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 265);

    auto g_xyyzz_yyyyz_1 = pbuffer.data(idx_eri_1_hh + 268);

    auto g_xyyzz_yyyzz_1 = pbuffer.data(idx_eri_1_hh + 269);

    auto g_xyyzz_yyzzz_1 = pbuffer.data(idx_eri_1_hh + 270);

    auto g_xyyzz_yzzzz_1 = pbuffer.data(idx_eri_1_hh + 271);

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

    auto g_yyyyz_xxxxz_1 = pbuffer.data(idx_eri_1_hh + 338);

    auto g_yyyyz_xxxyz_1 = pbuffer.data(idx_eri_1_hh + 340);

    auto g_yyyyz_xxxzz_1 = pbuffer.data(idx_eri_1_hh + 341);

    auto g_yyyyz_xxyyz_1 = pbuffer.data(idx_eri_1_hh + 343);

    auto g_yyyyz_xxyzz_1 = pbuffer.data(idx_eri_1_hh + 344);

    auto g_yyyyz_xxzzz_1 = pbuffer.data(idx_eri_1_hh + 345);

    auto g_yyyyz_xyyyz_1 = pbuffer.data(idx_eri_1_hh + 347);

    auto g_yyyyz_xyyzz_1 = pbuffer.data(idx_eri_1_hh + 348);

    auto g_yyyyz_xyzzz_1 = pbuffer.data(idx_eri_1_hh + 349);

    auto g_yyyyz_xzzzz_1 = pbuffer.data(idx_eri_1_hh + 350);

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

    // Set up components of auxiliary buffer : HI

    auto g_xxxxx_xxxxxx_1 = pbuffer.data(idx_eri_1_hi);

    auto g_xxxxx_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 1);

    auto g_xxxxx_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 2);

    auto g_xxxxx_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 3);

    auto g_xxxxx_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 4);

    auto g_xxxxx_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 5);

    auto g_xxxxx_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 6);

    auto g_xxxxx_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 7);

    auto g_xxxxx_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 8);

    auto g_xxxxx_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 9);

    auto g_xxxxx_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 10);

    auto g_xxxxx_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 11);

    auto g_xxxxx_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 12);

    auto g_xxxxx_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 13);

    auto g_xxxxx_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 14);

    auto g_xxxxx_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 15);

    auto g_xxxxx_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 16);

    auto g_xxxxx_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 17);

    auto g_xxxxx_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 18);

    auto g_xxxxx_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 19);

    auto g_xxxxx_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 20);

    auto g_xxxxx_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 21);

    auto g_xxxxx_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 22);

    auto g_xxxxx_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 23);

    auto g_xxxxx_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 24);

    auto g_xxxxx_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 25);

    auto g_xxxxx_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 26);

    auto g_xxxxx_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 27);

    auto g_xxxxy_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 28);

    auto g_xxxxy_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 29);

    auto g_xxxxy_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 30);

    auto g_xxxxy_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 31);

    auto g_xxxxy_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 33);

    auto g_xxxxy_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 34);

    auto g_xxxxy_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 37);

    auto g_xxxxy_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 38);

    auto g_xxxxy_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 42);

    auto g_xxxxy_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 43);

    auto g_xxxxy_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 48);

    auto g_xxxxy_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 49);

    auto g_xxxxz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 56);

    auto g_xxxxz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 57);

    auto g_xxxxz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 58);

    auto g_xxxxz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 59);

    auto g_xxxxz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 60);

    auto g_xxxxz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 61);

    auto g_xxxxz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 62);

    auto g_xxxxz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 63);

    auto g_xxxxz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 64);

    auto g_xxxxz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 65);

    auto g_xxxxz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 66);

    auto g_xxxxz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 67);

    auto g_xxxxz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 68);

    auto g_xxxxz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 69);

    auto g_xxxxz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 70);

    auto g_xxxxz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 71);

    auto g_xxxxz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 72);

    auto g_xxxxz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 73);

    auto g_xxxxz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 74);

    auto g_xxxxz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 75);

    auto g_xxxxz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 76);

    auto g_xxxxz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 78);

    auto g_xxxxz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 79);

    auto g_xxxxz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 80);

    auto g_xxxxz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 81);

    auto g_xxxxz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 82);

    auto g_xxxxz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 83);

    auto g_xxxyy_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 84);

    auto g_xxxyy_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 85);

    auto g_xxxyy_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 86);

    auto g_xxxyy_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 87);

    auto g_xxxyy_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 88);

    auto g_xxxyy_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 89);

    auto g_xxxyy_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 90);

    auto g_xxxyy_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 91);

    auto g_xxxyy_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 92);

    auto g_xxxyy_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 93);

    auto g_xxxyy_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 94);

    auto g_xxxyy_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 95);

    auto g_xxxyy_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 96);

    auto g_xxxyy_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 97);

    auto g_xxxyy_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 98);

    auto g_xxxyy_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 99);

    auto g_xxxyy_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 100);

    auto g_xxxyy_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 101);

    auto g_xxxyy_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 102);

    auto g_xxxyy_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 103);

    auto g_xxxyy_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 104);

    auto g_xxxyy_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 105);

    auto g_xxxyy_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 106);

    auto g_xxxyy_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 107);

    auto g_xxxyy_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 108);

    auto g_xxxyy_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 109);

    auto g_xxxyy_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 110);

    auto g_xxxyy_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 111);

    auto g_xxxzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 140);

    auto g_xxxzz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 141);

    auto g_xxxzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 142);

    auto g_xxxzz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 143);

    auto g_xxxzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 144);

    auto g_xxxzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 145);

    auto g_xxxzz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 146);

    auto g_xxxzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 147);

    auto g_xxxzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 148);

    auto g_xxxzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 149);

    auto g_xxxzz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 150);

    auto g_xxxzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 151);

    auto g_xxxzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 152);

    auto g_xxxzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 153);

    auto g_xxxzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 154);

    auto g_xxxzz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 155);

    auto g_xxxzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 156);

    auto g_xxxzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 157);

    auto g_xxxzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 158);

    auto g_xxxzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 159);

    auto g_xxxzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 160);

    auto g_xxxzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 161);

    auto g_xxxzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 162);

    auto g_xxxzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 163);

    auto g_xxxzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 164);

    auto g_xxxzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 165);

    auto g_xxxzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 166);

    auto g_xxxzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 167);

    auto g_xxyyy_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 168);

    auto g_xxyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 169);

    auto g_xxyyy_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 170);

    auto g_xxyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 171);

    auto g_xxyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 172);

    auto g_xxyyy_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 173);

    auto g_xxyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 174);

    auto g_xxyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 175);

    auto g_xxyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 176);

    auto g_xxyyy_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 177);

    auto g_xxyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 178);

    auto g_xxyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 179);

    auto g_xxyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 180);

    auto g_xxyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 181);

    auto g_xxyyy_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 182);

    auto g_xxyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 183);

    auto g_xxyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 184);

    auto g_xxyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 185);

    auto g_xxyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 186);

    auto g_xxyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 187);

    auto g_xxyyy_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 188);

    auto g_xxyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 189);

    auto g_xxyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 190);

    auto g_xxyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 191);

    auto g_xxyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 192);

    auto g_xxyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 193);

    auto g_xxyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 194);

    auto g_xxyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 195);

    auto g_xxyyz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 197);

    auto g_xxyyz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 199);

    auto g_xxyyz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 202);

    auto g_xxyyz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 206);

    auto g_xxyyz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 211);

    auto g_xxyzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 224);

    auto g_xxyzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 226);

    auto g_xxyzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 229);

    auto g_xxyzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 233);

    auto g_xxyzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 238);

    auto g_xxyzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 244);

    auto g_xxzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 252);

    auto g_xxzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 253);

    auto g_xxzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 254);

    auto g_xxzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 255);

    auto g_xxzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 256);

    auto g_xxzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 257);

    auto g_xxzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 258);

    auto g_xxzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 259);

    auto g_xxzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 260);

    auto g_xxzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 261);

    auto g_xxzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 262);

    auto g_xxzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 263);

    auto g_xxzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 264);

    auto g_xxzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 265);

    auto g_xxzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 266);

    auto g_xxzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 267);

    auto g_xxzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 268);

    auto g_xxzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 269);

    auto g_xxzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 270);

    auto g_xxzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 271);

    auto g_xxzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 272);

    auto g_xxzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 273);

    auto g_xxzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 274);

    auto g_xxzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 275);

    auto g_xxzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 276);

    auto g_xxzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 277);

    auto g_xxzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 278);

    auto g_xxzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 279);

    auto g_xyyyy_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 280);

    auto g_xyyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 281);

    auto g_xyyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 283);

    auto g_xyyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 284);

    auto g_xyyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 286);

    auto g_xyyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 287);

    auto g_xyyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 288);

    auto g_xyyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 290);

    auto g_xyyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 291);

    auto g_xyyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 292);

    auto g_xyyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 293);

    auto g_xyyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 295);

    auto g_xyyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 296);

    auto g_xyyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 297);

    auto g_xyyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 298);

    auto g_xyyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 299);

    auto g_xyyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 301);

    auto g_xyyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 302);

    auto g_xyyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 303);

    auto g_xyyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 304);

    auto g_xyyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 305);

    auto g_xyyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 306);

    auto g_xyyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 307);

    auto g_xyyzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 340);

    auto g_xyyzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 343);

    auto g_xyyzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 344);

    auto g_xyyzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 347);

    auto g_xyyzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 348);

    auto g_xyyzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 349);

    auto g_xyyzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 352);

    auto g_xyyzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 353);

    auto g_xyyzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 354);

    auto g_xyyzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 355);

    auto g_xyyzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 357);

    auto g_xyyzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 358);

    auto g_xyyzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 359);

    auto g_xyyzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 360);

    auto g_xyyzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 361);

    auto g_xyyzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 362);

    auto g_xyyzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 363);

    auto g_xzzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 392);

    auto g_xzzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 394);

    auto g_xzzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 396);

    auto g_xzzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 397);

    auto g_xzzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 399);

    auto g_xzzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 400);

    auto g_xzzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 401);

    auto g_xzzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 403);

    auto g_xzzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 404);

    auto g_xzzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 405);

    auto g_xzzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 406);

    auto g_xzzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 408);

    auto g_xzzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 409);

    auto g_xzzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 410);

    auto g_xzzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 411);

    auto g_xzzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 412);

    auto g_xzzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 413);

    auto g_xzzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 414);

    auto g_xzzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 415);

    auto g_xzzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 416);

    auto g_xzzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 417);

    auto g_xzzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 418);

    auto g_xzzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 419);

    auto g_yyyyy_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 420);

    auto g_yyyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 421);

    auto g_yyyyy_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 422);

    auto g_yyyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 423);

    auto g_yyyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 424);

    auto g_yyyyy_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 425);

    auto g_yyyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 426);

    auto g_yyyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 427);

    auto g_yyyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 428);

    auto g_yyyyy_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 429);

    auto g_yyyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 430);

    auto g_yyyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 431);

    auto g_yyyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 432);

    auto g_yyyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 433);

    auto g_yyyyy_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 434);

    auto g_yyyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 435);

    auto g_yyyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 436);

    auto g_yyyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 437);

    auto g_yyyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 438);

    auto g_yyyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 439);

    auto g_yyyyy_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 440);

    auto g_yyyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 441);

    auto g_yyyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 442);

    auto g_yyyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 443);

    auto g_yyyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 444);

    auto g_yyyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 445);

    auto g_yyyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 446);

    auto g_yyyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 447);

    auto g_yyyyz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 449);

    auto g_yyyyz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 450);

    auto g_yyyyz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 451);

    auto g_yyyyz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 452);

    auto g_yyyyz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 453);

    auto g_yyyyz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 454);

    auto g_yyyyz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 455);

    auto g_yyyyz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 456);

    auto g_yyyyz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 457);

    auto g_yyyyz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 458);

    auto g_yyyyz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 459);

    auto g_yyyyz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 460);

    auto g_yyyyz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 461);

    auto g_yyyyz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 462);

    auto g_yyyyz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 463);

    auto g_yyyyz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 464);

    auto g_yyyyz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 465);

    auto g_yyyyz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 466);

    auto g_yyyyz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 467);

    auto g_yyyyz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 468);

    auto g_yyyyz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 469);

    auto g_yyyyz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 470);

    auto g_yyyyz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 471);

    auto g_yyyyz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 472);

    auto g_yyyyz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 473);

    auto g_yyyyz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 474);

    auto g_yyyyz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 475);

    auto g_yyyzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 476);

    auto g_yyyzz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 477);

    auto g_yyyzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 478);

    auto g_yyyzz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 479);

    auto g_yyyzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 480);

    auto g_yyyzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 481);

    auto g_yyyzz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 482);

    auto g_yyyzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 483);

    auto g_yyyzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 484);

    auto g_yyyzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 485);

    auto g_yyyzz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 486);

    auto g_yyyzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 487);

    auto g_yyyzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 488);

    auto g_yyyzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 489);

    auto g_yyyzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 490);

    auto g_yyyzz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 491);

    auto g_yyyzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 492);

    auto g_yyyzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 493);

    auto g_yyyzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 494);

    auto g_yyyzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 495);

    auto g_yyyzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 496);

    auto g_yyyzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 497);

    auto g_yyyzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 498);

    auto g_yyyzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 499);

    auto g_yyyzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 500);

    auto g_yyyzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 501);

    auto g_yyyzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 502);

    auto g_yyyzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 503);

    auto g_yyzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 504);

    auto g_yyzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 505);

    auto g_yyzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 506);

    auto g_yyzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 507);

    auto g_yyzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 508);

    auto g_yyzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 509);

    auto g_yyzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 510);

    auto g_yyzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 511);

    auto g_yyzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 512);

    auto g_yyzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 513);

    auto g_yyzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 514);

    auto g_yyzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 515);

    auto g_yyzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 516);

    auto g_yyzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 517);

    auto g_yyzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 518);

    auto g_yyzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 519);

    auto g_yyzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 520);

    auto g_yyzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 521);

    auto g_yyzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 522);

    auto g_yyzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 523);

    auto g_yyzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 524);

    auto g_yyzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 525);

    auto g_yyzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 526);

    auto g_yyzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 527);

    auto g_yyzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 528);

    auto g_yyzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 529);

    auto g_yyzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 530);

    auto g_yyzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 531);

    auto g_yzzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 532);

    auto g_yzzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 533);

    auto g_yzzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 534);

    auto g_yzzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 535);

    auto g_yzzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 536);

    auto g_yzzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 537);

    auto g_yzzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 538);

    auto g_yzzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 539);

    auto g_yzzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 540);

    auto g_yzzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 541);

    auto g_yzzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 542);

    auto g_yzzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 543);

    auto g_yzzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 544);

    auto g_yzzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 545);

    auto g_yzzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 546);

    auto g_yzzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 547);

    auto g_yzzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 548);

    auto g_yzzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 549);

    auto g_yzzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 550);

    auto g_yzzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 551);

    auto g_yzzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 552);

    auto g_yzzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 553);

    auto g_yzzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 554);

    auto g_yzzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 555);

    auto g_yzzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 556);

    auto g_yzzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 557);

    auto g_yzzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 558);

    auto g_yzzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 559);

    auto g_zzzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_hi + 560);

    auto g_zzzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_hi + 561);

    auto g_zzzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_hi + 562);

    auto g_zzzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_hi + 563);

    auto g_zzzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_hi + 564);

    auto g_zzzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_hi + 565);

    auto g_zzzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_hi + 566);

    auto g_zzzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_hi + 567);

    auto g_zzzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_hi + 568);

    auto g_zzzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_hi + 569);

    auto g_zzzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_hi + 570);

    auto g_zzzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_hi + 571);

    auto g_zzzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_hi + 572);

    auto g_zzzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_hi + 573);

    auto g_zzzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_hi + 574);

    auto g_zzzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_hi + 575);

    auto g_zzzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_hi + 576);

    auto g_zzzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_hi + 577);

    auto g_zzzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_hi + 578);

    auto g_zzzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_hi + 579);

    auto g_zzzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_hi + 580);

    auto g_zzzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_hi + 581);

    auto g_zzzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_hi + 582);

    auto g_zzzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_hi + 583);

    auto g_zzzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_hi + 584);

    auto g_zzzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_hi + 585);

    auto g_zzzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_hi + 586);

    auto g_zzzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_hi + 587);

    // Set up 0-28 components of targeted buffer : II

    auto g_xxxxxx_xxxxxx_0 = pbuffer.data(idx_eri_0_ii);

    auto g_xxxxxx_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 1);

    auto g_xxxxxx_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 2);

    auto g_xxxxxx_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 3);

    auto g_xxxxxx_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 4);

    auto g_xxxxxx_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 5);

    auto g_xxxxxx_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 6);

    auto g_xxxxxx_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 7);

    auto g_xxxxxx_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 8);

    auto g_xxxxxx_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 9);

    auto g_xxxxxx_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 10);

    auto g_xxxxxx_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 11);

    auto g_xxxxxx_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 12);

    auto g_xxxxxx_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 13);

    auto g_xxxxxx_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 14);

    auto g_xxxxxx_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 15);

    auto g_xxxxxx_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 16);

    auto g_xxxxxx_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 17);

    auto g_xxxxxx_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 18);

    auto g_xxxxxx_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 19);

    auto g_xxxxxx_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 20);

    auto g_xxxxxx_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 21);

    auto g_xxxxxx_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 22);

    auto g_xxxxxx_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 23);

    auto g_xxxxxx_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 24);

    auto g_xxxxxx_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 25);

    auto g_xxxxxx_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 26);

    auto g_xxxxxx_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 27);

    #pragma omp simd aligned(g_xxxx_xxxxxx_0, g_xxxx_xxxxxx_1, g_xxxx_xxxxxy_0, g_xxxx_xxxxxy_1, g_xxxx_xxxxxz_0, g_xxxx_xxxxxz_1, g_xxxx_xxxxyy_0, g_xxxx_xxxxyy_1, g_xxxx_xxxxyz_0, g_xxxx_xxxxyz_1, g_xxxx_xxxxzz_0, g_xxxx_xxxxzz_1, g_xxxx_xxxyyy_0, g_xxxx_xxxyyy_1, g_xxxx_xxxyyz_0, g_xxxx_xxxyyz_1, g_xxxx_xxxyzz_0, g_xxxx_xxxyzz_1, g_xxxx_xxxzzz_0, g_xxxx_xxxzzz_1, g_xxxx_xxyyyy_0, g_xxxx_xxyyyy_1, g_xxxx_xxyyyz_0, g_xxxx_xxyyyz_1, g_xxxx_xxyyzz_0, g_xxxx_xxyyzz_1, g_xxxx_xxyzzz_0, g_xxxx_xxyzzz_1, g_xxxx_xxzzzz_0, g_xxxx_xxzzzz_1, g_xxxx_xyyyyy_0, g_xxxx_xyyyyy_1, g_xxxx_xyyyyz_0, g_xxxx_xyyyyz_1, g_xxxx_xyyyzz_0, g_xxxx_xyyyzz_1, g_xxxx_xyyzzz_0, g_xxxx_xyyzzz_1, g_xxxx_xyzzzz_0, g_xxxx_xyzzzz_1, g_xxxx_xzzzzz_0, g_xxxx_xzzzzz_1, g_xxxx_yyyyyy_0, g_xxxx_yyyyyy_1, g_xxxx_yyyyyz_0, g_xxxx_yyyyyz_1, g_xxxx_yyyyzz_0, g_xxxx_yyyyzz_1, g_xxxx_yyyzzz_0, g_xxxx_yyyzzz_1, g_xxxx_yyzzzz_0, g_xxxx_yyzzzz_1, g_xxxx_yzzzzz_0, g_xxxx_yzzzzz_1, g_xxxx_zzzzzz_0, g_xxxx_zzzzzz_1, g_xxxxx_xxxxx_1, g_xxxxx_xxxxxx_1, g_xxxxx_xxxxxy_1, g_xxxxx_xxxxxz_1, g_xxxxx_xxxxy_1, g_xxxxx_xxxxyy_1, g_xxxxx_xxxxyz_1, g_xxxxx_xxxxz_1, g_xxxxx_xxxxzz_1, g_xxxxx_xxxyy_1, g_xxxxx_xxxyyy_1, g_xxxxx_xxxyyz_1, g_xxxxx_xxxyz_1, g_xxxxx_xxxyzz_1, g_xxxxx_xxxzz_1, g_xxxxx_xxxzzz_1, g_xxxxx_xxyyy_1, g_xxxxx_xxyyyy_1, g_xxxxx_xxyyyz_1, g_xxxxx_xxyyz_1, g_xxxxx_xxyyzz_1, g_xxxxx_xxyzz_1, g_xxxxx_xxyzzz_1, g_xxxxx_xxzzz_1, g_xxxxx_xxzzzz_1, g_xxxxx_xyyyy_1, g_xxxxx_xyyyyy_1, g_xxxxx_xyyyyz_1, g_xxxxx_xyyyz_1, g_xxxxx_xyyyzz_1, g_xxxxx_xyyzz_1, g_xxxxx_xyyzzz_1, g_xxxxx_xyzzz_1, g_xxxxx_xyzzzz_1, g_xxxxx_xzzzz_1, g_xxxxx_xzzzzz_1, g_xxxxx_yyyyy_1, g_xxxxx_yyyyyy_1, g_xxxxx_yyyyyz_1, g_xxxxx_yyyyz_1, g_xxxxx_yyyyzz_1, g_xxxxx_yyyzz_1, g_xxxxx_yyyzzz_1, g_xxxxx_yyzzz_1, g_xxxxx_yyzzzz_1, g_xxxxx_yzzzz_1, g_xxxxx_yzzzzz_1, g_xxxxx_zzzzz_1, g_xxxxx_zzzzzz_1, g_xxxxxx_xxxxxx_0, g_xxxxxx_xxxxxy_0, g_xxxxxx_xxxxxz_0, g_xxxxxx_xxxxyy_0, g_xxxxxx_xxxxyz_0, g_xxxxxx_xxxxzz_0, g_xxxxxx_xxxyyy_0, g_xxxxxx_xxxyyz_0, g_xxxxxx_xxxyzz_0, g_xxxxxx_xxxzzz_0, g_xxxxxx_xxyyyy_0, g_xxxxxx_xxyyyz_0, g_xxxxxx_xxyyzz_0, g_xxxxxx_xxyzzz_0, g_xxxxxx_xxzzzz_0, g_xxxxxx_xyyyyy_0, g_xxxxxx_xyyyyz_0, g_xxxxxx_xyyyzz_0, g_xxxxxx_xyyzzz_0, g_xxxxxx_xyzzzz_0, g_xxxxxx_xzzzzz_0, g_xxxxxx_yyyyyy_0, g_xxxxxx_yyyyyz_0, g_xxxxxx_yyyyzz_0, g_xxxxxx_yyyzzz_0, g_xxxxxx_yyzzzz_0, g_xxxxxx_yzzzzz_0, g_xxxxxx_zzzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxx_xxxxxx_0[i] = 5.0 * g_xxxx_xxxxxx_0[i] * fbe_0 - 5.0 * g_xxxx_xxxxxx_1[i] * fz_be_0 + 6.0 * g_xxxxx_xxxxx_1[i] * fe_0 + g_xxxxx_xxxxxx_1[i] * pa_x[i];

        g_xxxxxx_xxxxxy_0[i] = 5.0 * g_xxxx_xxxxxy_0[i] * fbe_0 - 5.0 * g_xxxx_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxxxx_xxxxy_1[i] * fe_0 + g_xxxxx_xxxxxy_1[i] * pa_x[i];

        g_xxxxxx_xxxxxz_0[i] = 5.0 * g_xxxx_xxxxxz_0[i] * fbe_0 - 5.0 * g_xxxx_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxxxx_xxxxz_1[i] * fe_0 + g_xxxxx_xxxxxz_1[i] * pa_x[i];

        g_xxxxxx_xxxxyy_0[i] = 5.0 * g_xxxx_xxxxyy_0[i] * fbe_0 - 5.0 * g_xxxx_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxxxx_xxxyy_1[i] * fe_0 + g_xxxxx_xxxxyy_1[i] * pa_x[i];

        g_xxxxxx_xxxxyz_0[i] = 5.0 * g_xxxx_xxxxyz_0[i] * fbe_0 - 5.0 * g_xxxx_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxxx_xxxyz_1[i] * fe_0 + g_xxxxx_xxxxyz_1[i] * pa_x[i];

        g_xxxxxx_xxxxzz_0[i] = 5.0 * g_xxxx_xxxxzz_0[i] * fbe_0 - 5.0 * g_xxxx_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxxxx_xxxzz_1[i] * fe_0 + g_xxxxx_xxxxzz_1[i] * pa_x[i];

        g_xxxxxx_xxxyyy_0[i] = 5.0 * g_xxxx_xxxyyy_0[i] * fbe_0 - 5.0 * g_xxxx_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxxxx_xxyyy_1[i] * fe_0 + g_xxxxx_xxxyyy_1[i] * pa_x[i];

        g_xxxxxx_xxxyyz_0[i] = 5.0 * g_xxxx_xxxyyz_0[i] * fbe_0 - 5.0 * g_xxxx_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxxx_xxyyz_1[i] * fe_0 + g_xxxxx_xxxyyz_1[i] * pa_x[i];

        g_xxxxxx_xxxyzz_0[i] = 5.0 * g_xxxx_xxxyzz_0[i] * fbe_0 - 5.0 * g_xxxx_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_xxyzz_1[i] * fe_0 + g_xxxxx_xxxyzz_1[i] * pa_x[i];

        g_xxxxxx_xxxzzz_0[i] = 5.0 * g_xxxx_xxxzzz_0[i] * fbe_0 - 5.0 * g_xxxx_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_xxzzz_1[i] * fe_0 + g_xxxxx_xxxzzz_1[i] * pa_x[i];

        g_xxxxxx_xxyyyy_0[i] = 5.0 * g_xxxx_xxyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxxxx_xyyyy_1[i] * fe_0 + g_xxxxx_xxyyyy_1[i] * pa_x[i];

        g_xxxxxx_xxyyyz_0[i] = 5.0 * g_xxxx_xxyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxx_xyyyz_1[i] * fe_0 + g_xxxxx_xxyyyz_1[i] * pa_x[i];

        g_xxxxxx_xxyyzz_0[i] = 5.0 * g_xxxx_xxyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_xyyzz_1[i] * fe_0 + g_xxxxx_xxyyzz_1[i] * pa_x[i];

        g_xxxxxx_xxyzzz_0[i] = 5.0 * g_xxxx_xxyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_xyzzz_1[i] * fe_0 + g_xxxxx_xxyzzz_1[i] * pa_x[i];

        g_xxxxxx_xxzzzz_0[i] = 5.0 * g_xxxx_xxzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_xzzzz_1[i] * fe_0 + g_xxxxx_xxzzzz_1[i] * pa_x[i];

        g_xxxxxx_xyyyyy_0[i] = 5.0 * g_xxxx_xyyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_xyyyyy_1[i] * fz_be_0 + g_xxxxx_yyyyy_1[i] * fe_0 + g_xxxxx_xyyyyy_1[i] * pa_x[i];

        g_xxxxxx_xyyyyz_0[i] = 5.0 * g_xxxx_xyyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_xyyyyz_1[i] * fz_be_0 + g_xxxxx_yyyyz_1[i] * fe_0 + g_xxxxx_xyyyyz_1[i] * pa_x[i];

        g_xxxxxx_xyyyzz_0[i] = 5.0 * g_xxxx_xyyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_xyyyzz_1[i] * fz_be_0 + g_xxxxx_yyyzz_1[i] * fe_0 + g_xxxxx_xyyyzz_1[i] * pa_x[i];

        g_xxxxxx_xyyzzz_0[i] = 5.0 * g_xxxx_xyyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_xyyzzz_1[i] * fz_be_0 + g_xxxxx_yyzzz_1[i] * fe_0 + g_xxxxx_xyyzzz_1[i] * pa_x[i];

        g_xxxxxx_xyzzzz_0[i] = 5.0 * g_xxxx_xyzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_xyzzzz_1[i] * fz_be_0 + g_xxxxx_yzzzz_1[i] * fe_0 + g_xxxxx_xyzzzz_1[i] * pa_x[i];

        g_xxxxxx_xzzzzz_0[i] = 5.0 * g_xxxx_xzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_xzzzzz_1[i] * fz_be_0 + g_xxxxx_zzzzz_1[i] * fe_0 + g_xxxxx_xzzzzz_1[i] * pa_x[i];

        g_xxxxxx_yyyyyy_0[i] = 5.0 * g_xxxx_yyyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_yyyyyy_1[i] * fz_be_0 + g_xxxxx_yyyyyy_1[i] * pa_x[i];

        g_xxxxxx_yyyyyz_0[i] = 5.0 * g_xxxx_yyyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_yyyyyz_1[i] * fz_be_0 + g_xxxxx_yyyyyz_1[i] * pa_x[i];

        g_xxxxxx_yyyyzz_0[i] = 5.0 * g_xxxx_yyyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_yyyyzz_1[i] * fz_be_0 + g_xxxxx_yyyyzz_1[i] * pa_x[i];

        g_xxxxxx_yyyzzz_0[i] = 5.0 * g_xxxx_yyyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_yyyzzz_1[i] * fz_be_0 + g_xxxxx_yyyzzz_1[i] * pa_x[i];

        g_xxxxxx_yyzzzz_0[i] = 5.0 * g_xxxx_yyzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_yyzzzz_1[i] * fz_be_0 + g_xxxxx_yyzzzz_1[i] * pa_x[i];

        g_xxxxxx_yzzzzz_0[i] = 5.0 * g_xxxx_yzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_yzzzzz_1[i] * fz_be_0 + g_xxxxx_yzzzzz_1[i] * pa_x[i];

        g_xxxxxx_zzzzzz_0[i] = 5.0 * g_xxxx_zzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_zzzzzz_1[i] * fz_be_0 + g_xxxxx_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : II

    auto g_xxxxxy_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 28);

    auto g_xxxxxy_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 29);

    auto g_xxxxxy_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 30);

    auto g_xxxxxy_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 31);

    auto g_xxxxxy_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 32);

    auto g_xxxxxy_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 33);

    auto g_xxxxxy_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 34);

    auto g_xxxxxy_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 35);

    auto g_xxxxxy_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 36);

    auto g_xxxxxy_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 37);

    auto g_xxxxxy_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 38);

    auto g_xxxxxy_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 39);

    auto g_xxxxxy_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 40);

    auto g_xxxxxy_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 41);

    auto g_xxxxxy_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 42);

    auto g_xxxxxy_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 43);

    auto g_xxxxxy_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 44);

    auto g_xxxxxy_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 45);

    auto g_xxxxxy_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 46);

    auto g_xxxxxy_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 47);

    auto g_xxxxxy_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 48);

    auto g_xxxxxy_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 49);

    auto g_xxxxxy_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 50);

    auto g_xxxxxy_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 51);

    auto g_xxxxxy_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 52);

    auto g_xxxxxy_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 53);

    auto g_xxxxxy_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 54);

    auto g_xxxxxy_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 55);

    #pragma omp simd aligned(g_xxxxx_xxxxx_1, g_xxxxx_xxxxxx_1, g_xxxxx_xxxxxy_1, g_xxxxx_xxxxxz_1, g_xxxxx_xxxxy_1, g_xxxxx_xxxxyy_1, g_xxxxx_xxxxyz_1, g_xxxxx_xxxxz_1, g_xxxxx_xxxxzz_1, g_xxxxx_xxxyy_1, g_xxxxx_xxxyyy_1, g_xxxxx_xxxyyz_1, g_xxxxx_xxxyz_1, g_xxxxx_xxxyzz_1, g_xxxxx_xxxzz_1, g_xxxxx_xxxzzz_1, g_xxxxx_xxyyy_1, g_xxxxx_xxyyyy_1, g_xxxxx_xxyyyz_1, g_xxxxx_xxyyz_1, g_xxxxx_xxyyzz_1, g_xxxxx_xxyzz_1, g_xxxxx_xxyzzz_1, g_xxxxx_xxzzz_1, g_xxxxx_xxzzzz_1, g_xxxxx_xyyyy_1, g_xxxxx_xyyyyy_1, g_xxxxx_xyyyyz_1, g_xxxxx_xyyyz_1, g_xxxxx_xyyyzz_1, g_xxxxx_xyyzz_1, g_xxxxx_xyyzzz_1, g_xxxxx_xyzzz_1, g_xxxxx_xyzzzz_1, g_xxxxx_xzzzz_1, g_xxxxx_xzzzzz_1, g_xxxxx_yyyyy_1, g_xxxxx_yyyyyy_1, g_xxxxx_yyyyyz_1, g_xxxxx_yyyyz_1, g_xxxxx_yyyyzz_1, g_xxxxx_yyyzz_1, g_xxxxx_yyyzzz_1, g_xxxxx_yyzzz_1, g_xxxxx_yyzzzz_1, g_xxxxx_yzzzz_1, g_xxxxx_yzzzzz_1, g_xxxxx_zzzzz_1, g_xxxxx_zzzzzz_1, g_xxxxxy_xxxxxx_0, g_xxxxxy_xxxxxy_0, g_xxxxxy_xxxxxz_0, g_xxxxxy_xxxxyy_0, g_xxxxxy_xxxxyz_0, g_xxxxxy_xxxxzz_0, g_xxxxxy_xxxyyy_0, g_xxxxxy_xxxyyz_0, g_xxxxxy_xxxyzz_0, g_xxxxxy_xxxzzz_0, g_xxxxxy_xxyyyy_0, g_xxxxxy_xxyyyz_0, g_xxxxxy_xxyyzz_0, g_xxxxxy_xxyzzz_0, g_xxxxxy_xxzzzz_0, g_xxxxxy_xyyyyy_0, g_xxxxxy_xyyyyz_0, g_xxxxxy_xyyyzz_0, g_xxxxxy_xyyzzz_0, g_xxxxxy_xyzzzz_0, g_xxxxxy_xzzzzz_0, g_xxxxxy_yyyyyy_0, g_xxxxxy_yyyyyz_0, g_xxxxxy_yyyyzz_0, g_xxxxxy_yyyzzz_0, g_xxxxxy_yyzzzz_0, g_xxxxxy_yzzzzz_0, g_xxxxxy_zzzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxy_xxxxxx_0[i] = g_xxxxx_xxxxxx_1[i] * pa_y[i];

        g_xxxxxy_xxxxxy_0[i] = g_xxxxx_xxxxx_1[i] * fe_0 + g_xxxxx_xxxxxy_1[i] * pa_y[i];

        g_xxxxxy_xxxxxz_0[i] = g_xxxxx_xxxxxz_1[i] * pa_y[i];

        g_xxxxxy_xxxxyy_0[i] = 2.0 * g_xxxxx_xxxxy_1[i] * fe_0 + g_xxxxx_xxxxyy_1[i] * pa_y[i];

        g_xxxxxy_xxxxyz_0[i] = g_xxxxx_xxxxz_1[i] * fe_0 + g_xxxxx_xxxxyz_1[i] * pa_y[i];

        g_xxxxxy_xxxxzz_0[i] = g_xxxxx_xxxxzz_1[i] * pa_y[i];

        g_xxxxxy_xxxyyy_0[i] = 3.0 * g_xxxxx_xxxyy_1[i] * fe_0 + g_xxxxx_xxxyyy_1[i] * pa_y[i];

        g_xxxxxy_xxxyyz_0[i] = 2.0 * g_xxxxx_xxxyz_1[i] * fe_0 + g_xxxxx_xxxyyz_1[i] * pa_y[i];

        g_xxxxxy_xxxyzz_0[i] = g_xxxxx_xxxzz_1[i] * fe_0 + g_xxxxx_xxxyzz_1[i] * pa_y[i];

        g_xxxxxy_xxxzzz_0[i] = g_xxxxx_xxxzzz_1[i] * pa_y[i];

        g_xxxxxy_xxyyyy_0[i] = 4.0 * g_xxxxx_xxyyy_1[i] * fe_0 + g_xxxxx_xxyyyy_1[i] * pa_y[i];

        g_xxxxxy_xxyyyz_0[i] = 3.0 * g_xxxxx_xxyyz_1[i] * fe_0 + g_xxxxx_xxyyyz_1[i] * pa_y[i];

        g_xxxxxy_xxyyzz_0[i] = 2.0 * g_xxxxx_xxyzz_1[i] * fe_0 + g_xxxxx_xxyyzz_1[i] * pa_y[i];

        g_xxxxxy_xxyzzz_0[i] = g_xxxxx_xxzzz_1[i] * fe_0 + g_xxxxx_xxyzzz_1[i] * pa_y[i];

        g_xxxxxy_xxzzzz_0[i] = g_xxxxx_xxzzzz_1[i] * pa_y[i];

        g_xxxxxy_xyyyyy_0[i] = 5.0 * g_xxxxx_xyyyy_1[i] * fe_0 + g_xxxxx_xyyyyy_1[i] * pa_y[i];

        g_xxxxxy_xyyyyz_0[i] = 4.0 * g_xxxxx_xyyyz_1[i] * fe_0 + g_xxxxx_xyyyyz_1[i] * pa_y[i];

        g_xxxxxy_xyyyzz_0[i] = 3.0 * g_xxxxx_xyyzz_1[i] * fe_0 + g_xxxxx_xyyyzz_1[i] * pa_y[i];

        g_xxxxxy_xyyzzz_0[i] = 2.0 * g_xxxxx_xyzzz_1[i] * fe_0 + g_xxxxx_xyyzzz_1[i] * pa_y[i];

        g_xxxxxy_xyzzzz_0[i] = g_xxxxx_xzzzz_1[i] * fe_0 + g_xxxxx_xyzzzz_1[i] * pa_y[i];

        g_xxxxxy_xzzzzz_0[i] = g_xxxxx_xzzzzz_1[i] * pa_y[i];

        g_xxxxxy_yyyyyy_0[i] = 6.0 * g_xxxxx_yyyyy_1[i] * fe_0 + g_xxxxx_yyyyyy_1[i] * pa_y[i];

        g_xxxxxy_yyyyyz_0[i] = 5.0 * g_xxxxx_yyyyz_1[i] * fe_0 + g_xxxxx_yyyyyz_1[i] * pa_y[i];

        g_xxxxxy_yyyyzz_0[i] = 4.0 * g_xxxxx_yyyzz_1[i] * fe_0 + g_xxxxx_yyyyzz_1[i] * pa_y[i];

        g_xxxxxy_yyyzzz_0[i] = 3.0 * g_xxxxx_yyzzz_1[i] * fe_0 + g_xxxxx_yyyzzz_1[i] * pa_y[i];

        g_xxxxxy_yyzzzz_0[i] = 2.0 * g_xxxxx_yzzzz_1[i] * fe_0 + g_xxxxx_yyzzzz_1[i] * pa_y[i];

        g_xxxxxy_yzzzzz_0[i] = g_xxxxx_zzzzz_1[i] * fe_0 + g_xxxxx_yzzzzz_1[i] * pa_y[i];

        g_xxxxxy_zzzzzz_0[i] = g_xxxxx_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : II

    auto g_xxxxxz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 56);

    auto g_xxxxxz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 57);

    auto g_xxxxxz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 58);

    auto g_xxxxxz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 59);

    auto g_xxxxxz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 60);

    auto g_xxxxxz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 61);

    auto g_xxxxxz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 62);

    auto g_xxxxxz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 63);

    auto g_xxxxxz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 64);

    auto g_xxxxxz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 65);

    auto g_xxxxxz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 66);

    auto g_xxxxxz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 67);

    auto g_xxxxxz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 68);

    auto g_xxxxxz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 69);

    auto g_xxxxxz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 70);

    auto g_xxxxxz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 71);

    auto g_xxxxxz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 72);

    auto g_xxxxxz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 73);

    auto g_xxxxxz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 74);

    auto g_xxxxxz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 75);

    auto g_xxxxxz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 76);

    auto g_xxxxxz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 77);

    auto g_xxxxxz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 78);

    auto g_xxxxxz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 79);

    auto g_xxxxxz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 80);

    auto g_xxxxxz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 81);

    auto g_xxxxxz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 82);

    auto g_xxxxxz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 83);

    #pragma omp simd aligned(g_xxxxx_xxxxx_1, g_xxxxx_xxxxxx_1, g_xxxxx_xxxxxy_1, g_xxxxx_xxxxxz_1, g_xxxxx_xxxxy_1, g_xxxxx_xxxxyy_1, g_xxxxx_xxxxyz_1, g_xxxxx_xxxxz_1, g_xxxxx_xxxxzz_1, g_xxxxx_xxxyy_1, g_xxxxx_xxxyyy_1, g_xxxxx_xxxyyz_1, g_xxxxx_xxxyz_1, g_xxxxx_xxxyzz_1, g_xxxxx_xxxzz_1, g_xxxxx_xxxzzz_1, g_xxxxx_xxyyy_1, g_xxxxx_xxyyyy_1, g_xxxxx_xxyyyz_1, g_xxxxx_xxyyz_1, g_xxxxx_xxyyzz_1, g_xxxxx_xxyzz_1, g_xxxxx_xxyzzz_1, g_xxxxx_xxzzz_1, g_xxxxx_xxzzzz_1, g_xxxxx_xyyyy_1, g_xxxxx_xyyyyy_1, g_xxxxx_xyyyyz_1, g_xxxxx_xyyyz_1, g_xxxxx_xyyyzz_1, g_xxxxx_xyyzz_1, g_xxxxx_xyyzzz_1, g_xxxxx_xyzzz_1, g_xxxxx_xyzzzz_1, g_xxxxx_xzzzz_1, g_xxxxx_xzzzzz_1, g_xxxxx_yyyyy_1, g_xxxxx_yyyyyy_1, g_xxxxx_yyyyyz_1, g_xxxxx_yyyyz_1, g_xxxxx_yyyyzz_1, g_xxxxx_yyyzz_1, g_xxxxx_yyyzzz_1, g_xxxxx_yyzzz_1, g_xxxxx_yyzzzz_1, g_xxxxx_yzzzz_1, g_xxxxx_yzzzzz_1, g_xxxxx_zzzzz_1, g_xxxxx_zzzzzz_1, g_xxxxxz_xxxxxx_0, g_xxxxxz_xxxxxy_0, g_xxxxxz_xxxxxz_0, g_xxxxxz_xxxxyy_0, g_xxxxxz_xxxxyz_0, g_xxxxxz_xxxxzz_0, g_xxxxxz_xxxyyy_0, g_xxxxxz_xxxyyz_0, g_xxxxxz_xxxyzz_0, g_xxxxxz_xxxzzz_0, g_xxxxxz_xxyyyy_0, g_xxxxxz_xxyyyz_0, g_xxxxxz_xxyyzz_0, g_xxxxxz_xxyzzz_0, g_xxxxxz_xxzzzz_0, g_xxxxxz_xyyyyy_0, g_xxxxxz_xyyyyz_0, g_xxxxxz_xyyyzz_0, g_xxxxxz_xyyzzz_0, g_xxxxxz_xyzzzz_0, g_xxxxxz_xzzzzz_0, g_xxxxxz_yyyyyy_0, g_xxxxxz_yyyyyz_0, g_xxxxxz_yyyyzz_0, g_xxxxxz_yyyzzz_0, g_xxxxxz_yyzzzz_0, g_xxxxxz_yzzzzz_0, g_xxxxxz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxz_xxxxxx_0[i] = g_xxxxx_xxxxxx_1[i] * pa_z[i];

        g_xxxxxz_xxxxxy_0[i] = g_xxxxx_xxxxxy_1[i] * pa_z[i];

        g_xxxxxz_xxxxxz_0[i] = g_xxxxx_xxxxx_1[i] * fe_0 + g_xxxxx_xxxxxz_1[i] * pa_z[i];

        g_xxxxxz_xxxxyy_0[i] = g_xxxxx_xxxxyy_1[i] * pa_z[i];

        g_xxxxxz_xxxxyz_0[i] = g_xxxxx_xxxxy_1[i] * fe_0 + g_xxxxx_xxxxyz_1[i] * pa_z[i];

        g_xxxxxz_xxxxzz_0[i] = 2.0 * g_xxxxx_xxxxz_1[i] * fe_0 + g_xxxxx_xxxxzz_1[i] * pa_z[i];

        g_xxxxxz_xxxyyy_0[i] = g_xxxxx_xxxyyy_1[i] * pa_z[i];

        g_xxxxxz_xxxyyz_0[i] = g_xxxxx_xxxyy_1[i] * fe_0 + g_xxxxx_xxxyyz_1[i] * pa_z[i];

        g_xxxxxz_xxxyzz_0[i] = 2.0 * g_xxxxx_xxxyz_1[i] * fe_0 + g_xxxxx_xxxyzz_1[i] * pa_z[i];

        g_xxxxxz_xxxzzz_0[i] = 3.0 * g_xxxxx_xxxzz_1[i] * fe_0 + g_xxxxx_xxxzzz_1[i] * pa_z[i];

        g_xxxxxz_xxyyyy_0[i] = g_xxxxx_xxyyyy_1[i] * pa_z[i];

        g_xxxxxz_xxyyyz_0[i] = g_xxxxx_xxyyy_1[i] * fe_0 + g_xxxxx_xxyyyz_1[i] * pa_z[i];

        g_xxxxxz_xxyyzz_0[i] = 2.0 * g_xxxxx_xxyyz_1[i] * fe_0 + g_xxxxx_xxyyzz_1[i] * pa_z[i];

        g_xxxxxz_xxyzzz_0[i] = 3.0 * g_xxxxx_xxyzz_1[i] * fe_0 + g_xxxxx_xxyzzz_1[i] * pa_z[i];

        g_xxxxxz_xxzzzz_0[i] = 4.0 * g_xxxxx_xxzzz_1[i] * fe_0 + g_xxxxx_xxzzzz_1[i] * pa_z[i];

        g_xxxxxz_xyyyyy_0[i] = g_xxxxx_xyyyyy_1[i] * pa_z[i];

        g_xxxxxz_xyyyyz_0[i] = g_xxxxx_xyyyy_1[i] * fe_0 + g_xxxxx_xyyyyz_1[i] * pa_z[i];

        g_xxxxxz_xyyyzz_0[i] = 2.0 * g_xxxxx_xyyyz_1[i] * fe_0 + g_xxxxx_xyyyzz_1[i] * pa_z[i];

        g_xxxxxz_xyyzzz_0[i] = 3.0 * g_xxxxx_xyyzz_1[i] * fe_0 + g_xxxxx_xyyzzz_1[i] * pa_z[i];

        g_xxxxxz_xyzzzz_0[i] = 4.0 * g_xxxxx_xyzzz_1[i] * fe_0 + g_xxxxx_xyzzzz_1[i] * pa_z[i];

        g_xxxxxz_xzzzzz_0[i] = 5.0 * g_xxxxx_xzzzz_1[i] * fe_0 + g_xxxxx_xzzzzz_1[i] * pa_z[i];

        g_xxxxxz_yyyyyy_0[i] = g_xxxxx_yyyyyy_1[i] * pa_z[i];

        g_xxxxxz_yyyyyz_0[i] = g_xxxxx_yyyyy_1[i] * fe_0 + g_xxxxx_yyyyyz_1[i] * pa_z[i];

        g_xxxxxz_yyyyzz_0[i] = 2.0 * g_xxxxx_yyyyz_1[i] * fe_0 + g_xxxxx_yyyyzz_1[i] * pa_z[i];

        g_xxxxxz_yyyzzz_0[i] = 3.0 * g_xxxxx_yyyzz_1[i] * fe_0 + g_xxxxx_yyyzzz_1[i] * pa_z[i];

        g_xxxxxz_yyzzzz_0[i] = 4.0 * g_xxxxx_yyzzz_1[i] * fe_0 + g_xxxxx_yyzzzz_1[i] * pa_z[i];

        g_xxxxxz_yzzzzz_0[i] = 5.0 * g_xxxxx_yzzzz_1[i] * fe_0 + g_xxxxx_yzzzzz_1[i] * pa_z[i];

        g_xxxxxz_zzzzzz_0[i] = 6.0 * g_xxxxx_zzzzz_1[i] * fe_0 + g_xxxxx_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 84-112 components of targeted buffer : II

    auto g_xxxxyy_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 84);

    auto g_xxxxyy_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 85);

    auto g_xxxxyy_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 86);

    auto g_xxxxyy_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 87);

    auto g_xxxxyy_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 88);

    auto g_xxxxyy_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 89);

    auto g_xxxxyy_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 90);

    auto g_xxxxyy_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 91);

    auto g_xxxxyy_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 92);

    auto g_xxxxyy_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 93);

    auto g_xxxxyy_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 94);

    auto g_xxxxyy_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 95);

    auto g_xxxxyy_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 96);

    auto g_xxxxyy_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 97);

    auto g_xxxxyy_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 98);

    auto g_xxxxyy_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 99);

    auto g_xxxxyy_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 100);

    auto g_xxxxyy_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 101);

    auto g_xxxxyy_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 102);

    auto g_xxxxyy_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 103);

    auto g_xxxxyy_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 104);

    auto g_xxxxyy_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 105);

    auto g_xxxxyy_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 106);

    auto g_xxxxyy_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 107);

    auto g_xxxxyy_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 108);

    auto g_xxxxyy_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 109);

    auto g_xxxxyy_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 110);

    auto g_xxxxyy_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 111);

    #pragma omp simd aligned(g_xxxx_xxxxxx_0, g_xxxx_xxxxxx_1, g_xxxx_xxxxxz_0, g_xxxx_xxxxxz_1, g_xxxx_xxxxzz_0, g_xxxx_xxxxzz_1, g_xxxx_xxxzzz_0, g_xxxx_xxxzzz_1, g_xxxx_xxzzzz_0, g_xxxx_xxzzzz_1, g_xxxx_xzzzzz_0, g_xxxx_xzzzzz_1, g_xxxxy_xxxxxx_1, g_xxxxy_xxxxxz_1, g_xxxxy_xxxxzz_1, g_xxxxy_xxxzzz_1, g_xxxxy_xxzzzz_1, g_xxxxy_xzzzzz_1, g_xxxxyy_xxxxxx_0, g_xxxxyy_xxxxxy_0, g_xxxxyy_xxxxxz_0, g_xxxxyy_xxxxyy_0, g_xxxxyy_xxxxyz_0, g_xxxxyy_xxxxzz_0, g_xxxxyy_xxxyyy_0, g_xxxxyy_xxxyyz_0, g_xxxxyy_xxxyzz_0, g_xxxxyy_xxxzzz_0, g_xxxxyy_xxyyyy_0, g_xxxxyy_xxyyyz_0, g_xxxxyy_xxyyzz_0, g_xxxxyy_xxyzzz_0, g_xxxxyy_xxzzzz_0, g_xxxxyy_xyyyyy_0, g_xxxxyy_xyyyyz_0, g_xxxxyy_xyyyzz_0, g_xxxxyy_xyyzzz_0, g_xxxxyy_xyzzzz_0, g_xxxxyy_xzzzzz_0, g_xxxxyy_yyyyyy_0, g_xxxxyy_yyyyyz_0, g_xxxxyy_yyyyzz_0, g_xxxxyy_yyyzzz_0, g_xxxxyy_yyzzzz_0, g_xxxxyy_yzzzzz_0, g_xxxxyy_zzzzzz_0, g_xxxyy_xxxxxy_1, g_xxxyy_xxxxy_1, g_xxxyy_xxxxyy_1, g_xxxyy_xxxxyz_1, g_xxxyy_xxxyy_1, g_xxxyy_xxxyyy_1, g_xxxyy_xxxyyz_1, g_xxxyy_xxxyz_1, g_xxxyy_xxxyzz_1, g_xxxyy_xxyyy_1, g_xxxyy_xxyyyy_1, g_xxxyy_xxyyyz_1, g_xxxyy_xxyyz_1, g_xxxyy_xxyyzz_1, g_xxxyy_xxyzz_1, g_xxxyy_xxyzzz_1, g_xxxyy_xyyyy_1, g_xxxyy_xyyyyy_1, g_xxxyy_xyyyyz_1, g_xxxyy_xyyyz_1, g_xxxyy_xyyyzz_1, g_xxxyy_xyyzz_1, g_xxxyy_xyyzzz_1, g_xxxyy_xyzzz_1, g_xxxyy_xyzzzz_1, g_xxxyy_yyyyy_1, g_xxxyy_yyyyyy_1, g_xxxyy_yyyyyz_1, g_xxxyy_yyyyz_1, g_xxxyy_yyyyzz_1, g_xxxyy_yyyzz_1, g_xxxyy_yyyzzz_1, g_xxxyy_yyzzz_1, g_xxxyy_yyzzzz_1, g_xxxyy_yzzzz_1, g_xxxyy_yzzzzz_1, g_xxxyy_zzzzzz_1, g_xxyy_xxxxxy_0, g_xxyy_xxxxxy_1, g_xxyy_xxxxyy_0, g_xxyy_xxxxyy_1, g_xxyy_xxxxyz_0, g_xxyy_xxxxyz_1, g_xxyy_xxxyyy_0, g_xxyy_xxxyyy_1, g_xxyy_xxxyyz_0, g_xxyy_xxxyyz_1, g_xxyy_xxxyzz_0, g_xxyy_xxxyzz_1, g_xxyy_xxyyyy_0, g_xxyy_xxyyyy_1, g_xxyy_xxyyyz_0, g_xxyy_xxyyyz_1, g_xxyy_xxyyzz_0, g_xxyy_xxyyzz_1, g_xxyy_xxyzzz_0, g_xxyy_xxyzzz_1, g_xxyy_xyyyyy_0, g_xxyy_xyyyyy_1, g_xxyy_xyyyyz_0, g_xxyy_xyyyyz_1, g_xxyy_xyyyzz_0, g_xxyy_xyyyzz_1, g_xxyy_xyyzzz_0, g_xxyy_xyyzzz_1, g_xxyy_xyzzzz_0, g_xxyy_xyzzzz_1, g_xxyy_yyyyyy_0, g_xxyy_yyyyyy_1, g_xxyy_yyyyyz_0, g_xxyy_yyyyyz_1, g_xxyy_yyyyzz_0, g_xxyy_yyyyzz_1, g_xxyy_yyyzzz_0, g_xxyy_yyyzzz_1, g_xxyy_yyzzzz_0, g_xxyy_yyzzzz_1, g_xxyy_yzzzzz_0, g_xxyy_yzzzzz_1, g_xxyy_zzzzzz_0, g_xxyy_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxyy_xxxxxx_0[i] = g_xxxx_xxxxxx_0[i] * fbe_0 - g_xxxx_xxxxxx_1[i] * fz_be_0 + g_xxxxy_xxxxxx_1[i] * pa_y[i];

        g_xxxxyy_xxxxxy_0[i] = 3.0 * g_xxyy_xxxxxy_0[i] * fbe_0 - 3.0 * g_xxyy_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxxyy_xxxxy_1[i] * fe_0 + g_xxxyy_xxxxxy_1[i] * pa_x[i];

        g_xxxxyy_xxxxxz_0[i] = g_xxxx_xxxxxz_0[i] * fbe_0 - g_xxxx_xxxxxz_1[i] * fz_be_0 + g_xxxxy_xxxxxz_1[i] * pa_y[i];

        g_xxxxyy_xxxxyy_0[i] = 3.0 * g_xxyy_xxxxyy_0[i] * fbe_0 - 3.0 * g_xxyy_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxxyy_xxxyy_1[i] * fe_0 + g_xxxyy_xxxxyy_1[i] * pa_x[i];

        g_xxxxyy_xxxxyz_0[i] = 3.0 * g_xxyy_xxxxyz_0[i] * fbe_0 - 3.0 * g_xxyy_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxyy_xxxyz_1[i] * fe_0 + g_xxxyy_xxxxyz_1[i] * pa_x[i];

        g_xxxxyy_xxxxzz_0[i] = g_xxxx_xxxxzz_0[i] * fbe_0 - g_xxxx_xxxxzz_1[i] * fz_be_0 + g_xxxxy_xxxxzz_1[i] * pa_y[i];

        g_xxxxyy_xxxyyy_0[i] = 3.0 * g_xxyy_xxxyyy_0[i] * fbe_0 - 3.0 * g_xxyy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxxyy_xxyyy_1[i] * fe_0 + g_xxxyy_xxxyyy_1[i] * pa_x[i];

        g_xxxxyy_xxxyyz_0[i] = 3.0 * g_xxyy_xxxyyz_0[i] * fbe_0 - 3.0 * g_xxyy_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxyy_xxyyz_1[i] * fe_0 + g_xxxyy_xxxyyz_1[i] * pa_x[i];

        g_xxxxyy_xxxyzz_0[i] = 3.0 * g_xxyy_xxxyzz_0[i] * fbe_0 - 3.0 * g_xxyy_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxyy_xxyzz_1[i] * fe_0 + g_xxxyy_xxxyzz_1[i] * pa_x[i];

        g_xxxxyy_xxxzzz_0[i] = g_xxxx_xxxzzz_0[i] * fbe_0 - g_xxxx_xxxzzz_1[i] * fz_be_0 + g_xxxxy_xxxzzz_1[i] * pa_y[i];

        g_xxxxyy_xxyyyy_0[i] = 3.0 * g_xxyy_xxyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxxyy_xyyyy_1[i] * fe_0 + g_xxxyy_xxyyyy_1[i] * pa_x[i];

        g_xxxxyy_xxyyyz_0[i] = 3.0 * g_xxyy_xxyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxyy_xyyyz_1[i] * fe_0 + g_xxxyy_xxyyyz_1[i] * pa_x[i];

        g_xxxxyy_xxyyzz_0[i] = 3.0 * g_xxyy_xxyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_xyyzz_1[i] * fe_0 + g_xxxyy_xxyyzz_1[i] * pa_x[i];

        g_xxxxyy_xxyzzz_0[i] = 3.0 * g_xxyy_xxyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_xyzzz_1[i] * fe_0 + g_xxxyy_xxyzzz_1[i] * pa_x[i];

        g_xxxxyy_xxzzzz_0[i] = g_xxxx_xxzzzz_0[i] * fbe_0 - g_xxxx_xxzzzz_1[i] * fz_be_0 + g_xxxxy_xxzzzz_1[i] * pa_y[i];

        g_xxxxyy_xyyyyy_0[i] = 3.0 * g_xxyy_xyyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_xyyyyy_1[i] * fz_be_0 + g_xxxyy_yyyyy_1[i] * fe_0 + g_xxxyy_xyyyyy_1[i] * pa_x[i];

        g_xxxxyy_xyyyyz_0[i] = 3.0 * g_xxyy_xyyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_xyyyyz_1[i] * fz_be_0 + g_xxxyy_yyyyz_1[i] * fe_0 + g_xxxyy_xyyyyz_1[i] * pa_x[i];

        g_xxxxyy_xyyyzz_0[i] = 3.0 * g_xxyy_xyyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_xyyyzz_1[i] * fz_be_0 + g_xxxyy_yyyzz_1[i] * fe_0 + g_xxxyy_xyyyzz_1[i] * pa_x[i];

        g_xxxxyy_xyyzzz_0[i] = 3.0 * g_xxyy_xyyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_xyyzzz_1[i] * fz_be_0 + g_xxxyy_yyzzz_1[i] * fe_0 + g_xxxyy_xyyzzz_1[i] * pa_x[i];

        g_xxxxyy_xyzzzz_0[i] = 3.0 * g_xxyy_xyzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_xyzzzz_1[i] * fz_be_0 + g_xxxyy_yzzzz_1[i] * fe_0 + g_xxxyy_xyzzzz_1[i] * pa_x[i];

        g_xxxxyy_xzzzzz_0[i] = g_xxxx_xzzzzz_0[i] * fbe_0 - g_xxxx_xzzzzz_1[i] * fz_be_0 + g_xxxxy_xzzzzz_1[i] * pa_y[i];

        g_xxxxyy_yyyyyy_0[i] = 3.0 * g_xxyy_yyyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_yyyyyy_1[i] * fz_be_0 + g_xxxyy_yyyyyy_1[i] * pa_x[i];

        g_xxxxyy_yyyyyz_0[i] = 3.0 * g_xxyy_yyyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_yyyyyz_1[i] * fz_be_0 + g_xxxyy_yyyyyz_1[i] * pa_x[i];

        g_xxxxyy_yyyyzz_0[i] = 3.0 * g_xxyy_yyyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_yyyyzz_1[i] * fz_be_0 + g_xxxyy_yyyyzz_1[i] * pa_x[i];

        g_xxxxyy_yyyzzz_0[i] = 3.0 * g_xxyy_yyyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_yyyzzz_1[i] * fz_be_0 + g_xxxyy_yyyzzz_1[i] * pa_x[i];

        g_xxxxyy_yyzzzz_0[i] = 3.0 * g_xxyy_yyzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_yyzzzz_1[i] * fz_be_0 + g_xxxyy_yyzzzz_1[i] * pa_x[i];

        g_xxxxyy_yzzzzz_0[i] = 3.0 * g_xxyy_yzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_yzzzzz_1[i] * fz_be_0 + g_xxxyy_yzzzzz_1[i] * pa_x[i];

        g_xxxxyy_zzzzzz_0[i] = 3.0 * g_xxyy_zzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_zzzzzz_1[i] * fz_be_0 + g_xxxyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 112-140 components of targeted buffer : II

    auto g_xxxxyz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 112);

    auto g_xxxxyz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 113);

    auto g_xxxxyz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 114);

    auto g_xxxxyz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 115);

    auto g_xxxxyz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 116);

    auto g_xxxxyz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 117);

    auto g_xxxxyz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 118);

    auto g_xxxxyz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 119);

    auto g_xxxxyz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 120);

    auto g_xxxxyz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 121);

    auto g_xxxxyz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 122);

    auto g_xxxxyz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 123);

    auto g_xxxxyz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 124);

    auto g_xxxxyz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 125);

    auto g_xxxxyz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 126);

    auto g_xxxxyz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 127);

    auto g_xxxxyz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 128);

    auto g_xxxxyz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 129);

    auto g_xxxxyz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 130);

    auto g_xxxxyz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 131);

    auto g_xxxxyz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 132);

    auto g_xxxxyz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 133);

    auto g_xxxxyz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 134);

    auto g_xxxxyz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 135);

    auto g_xxxxyz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 136);

    auto g_xxxxyz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 137);

    auto g_xxxxyz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 138);

    auto g_xxxxyz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 139);

    #pragma omp simd aligned(g_xxxxy_xxxxxy_1, g_xxxxy_xxxxyy_1, g_xxxxy_xxxyyy_1, g_xxxxy_xxyyyy_1, g_xxxxy_xyyyyy_1, g_xxxxy_yyyyyy_1, g_xxxxyz_xxxxxx_0, g_xxxxyz_xxxxxy_0, g_xxxxyz_xxxxxz_0, g_xxxxyz_xxxxyy_0, g_xxxxyz_xxxxyz_0, g_xxxxyz_xxxxzz_0, g_xxxxyz_xxxyyy_0, g_xxxxyz_xxxyyz_0, g_xxxxyz_xxxyzz_0, g_xxxxyz_xxxzzz_0, g_xxxxyz_xxyyyy_0, g_xxxxyz_xxyyyz_0, g_xxxxyz_xxyyzz_0, g_xxxxyz_xxyzzz_0, g_xxxxyz_xxzzzz_0, g_xxxxyz_xyyyyy_0, g_xxxxyz_xyyyyz_0, g_xxxxyz_xyyyzz_0, g_xxxxyz_xyyzzz_0, g_xxxxyz_xyzzzz_0, g_xxxxyz_xzzzzz_0, g_xxxxyz_yyyyyy_0, g_xxxxyz_yyyyyz_0, g_xxxxyz_yyyyzz_0, g_xxxxyz_yyyzzz_0, g_xxxxyz_yyzzzz_0, g_xxxxyz_yzzzzz_0, g_xxxxyz_zzzzzz_0, g_xxxxz_xxxxxx_1, g_xxxxz_xxxxxz_1, g_xxxxz_xxxxyz_1, g_xxxxz_xxxxz_1, g_xxxxz_xxxxzz_1, g_xxxxz_xxxyyz_1, g_xxxxz_xxxyz_1, g_xxxxz_xxxyzz_1, g_xxxxz_xxxzz_1, g_xxxxz_xxxzzz_1, g_xxxxz_xxyyyz_1, g_xxxxz_xxyyz_1, g_xxxxz_xxyyzz_1, g_xxxxz_xxyzz_1, g_xxxxz_xxyzzz_1, g_xxxxz_xxzzz_1, g_xxxxz_xxzzzz_1, g_xxxxz_xyyyyz_1, g_xxxxz_xyyyz_1, g_xxxxz_xyyyzz_1, g_xxxxz_xyyzz_1, g_xxxxz_xyyzzz_1, g_xxxxz_xyzzz_1, g_xxxxz_xyzzzz_1, g_xxxxz_xzzzz_1, g_xxxxz_xzzzzz_1, g_xxxxz_yyyyyz_1, g_xxxxz_yyyyz_1, g_xxxxz_yyyyzz_1, g_xxxxz_yyyzz_1, g_xxxxz_yyyzzz_1, g_xxxxz_yyzzz_1, g_xxxxz_yyzzzz_1, g_xxxxz_yzzzz_1, g_xxxxz_yzzzzz_1, g_xxxxz_zzzzz_1, g_xxxxz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyz_xxxxxx_0[i] = g_xxxxz_xxxxxx_1[i] * pa_y[i];

        g_xxxxyz_xxxxxy_0[i] = g_xxxxy_xxxxxy_1[i] * pa_z[i];

        g_xxxxyz_xxxxxz_0[i] = g_xxxxz_xxxxxz_1[i] * pa_y[i];

        g_xxxxyz_xxxxyy_0[i] = g_xxxxy_xxxxyy_1[i] * pa_z[i];

        g_xxxxyz_xxxxyz_0[i] = g_xxxxz_xxxxz_1[i] * fe_0 + g_xxxxz_xxxxyz_1[i] * pa_y[i];

        g_xxxxyz_xxxxzz_0[i] = g_xxxxz_xxxxzz_1[i] * pa_y[i];

        g_xxxxyz_xxxyyy_0[i] = g_xxxxy_xxxyyy_1[i] * pa_z[i];

        g_xxxxyz_xxxyyz_0[i] = 2.0 * g_xxxxz_xxxyz_1[i] * fe_0 + g_xxxxz_xxxyyz_1[i] * pa_y[i];

        g_xxxxyz_xxxyzz_0[i] = g_xxxxz_xxxzz_1[i] * fe_0 + g_xxxxz_xxxyzz_1[i] * pa_y[i];

        g_xxxxyz_xxxzzz_0[i] = g_xxxxz_xxxzzz_1[i] * pa_y[i];

        g_xxxxyz_xxyyyy_0[i] = g_xxxxy_xxyyyy_1[i] * pa_z[i];

        g_xxxxyz_xxyyyz_0[i] = 3.0 * g_xxxxz_xxyyz_1[i] * fe_0 + g_xxxxz_xxyyyz_1[i] * pa_y[i];

        g_xxxxyz_xxyyzz_0[i] = 2.0 * g_xxxxz_xxyzz_1[i] * fe_0 + g_xxxxz_xxyyzz_1[i] * pa_y[i];

        g_xxxxyz_xxyzzz_0[i] = g_xxxxz_xxzzz_1[i] * fe_0 + g_xxxxz_xxyzzz_1[i] * pa_y[i];

        g_xxxxyz_xxzzzz_0[i] = g_xxxxz_xxzzzz_1[i] * pa_y[i];

        g_xxxxyz_xyyyyy_0[i] = g_xxxxy_xyyyyy_1[i] * pa_z[i];

        g_xxxxyz_xyyyyz_0[i] = 4.0 * g_xxxxz_xyyyz_1[i] * fe_0 + g_xxxxz_xyyyyz_1[i] * pa_y[i];

        g_xxxxyz_xyyyzz_0[i] = 3.0 * g_xxxxz_xyyzz_1[i] * fe_0 + g_xxxxz_xyyyzz_1[i] * pa_y[i];

        g_xxxxyz_xyyzzz_0[i] = 2.0 * g_xxxxz_xyzzz_1[i] * fe_0 + g_xxxxz_xyyzzz_1[i] * pa_y[i];

        g_xxxxyz_xyzzzz_0[i] = g_xxxxz_xzzzz_1[i] * fe_0 + g_xxxxz_xyzzzz_1[i] * pa_y[i];

        g_xxxxyz_xzzzzz_0[i] = g_xxxxz_xzzzzz_1[i] * pa_y[i];

        g_xxxxyz_yyyyyy_0[i] = g_xxxxy_yyyyyy_1[i] * pa_z[i];

        g_xxxxyz_yyyyyz_0[i] = 5.0 * g_xxxxz_yyyyz_1[i] * fe_0 + g_xxxxz_yyyyyz_1[i] * pa_y[i];

        g_xxxxyz_yyyyzz_0[i] = 4.0 * g_xxxxz_yyyzz_1[i] * fe_0 + g_xxxxz_yyyyzz_1[i] * pa_y[i];

        g_xxxxyz_yyyzzz_0[i] = 3.0 * g_xxxxz_yyzzz_1[i] * fe_0 + g_xxxxz_yyyzzz_1[i] * pa_y[i];

        g_xxxxyz_yyzzzz_0[i] = 2.0 * g_xxxxz_yzzzz_1[i] * fe_0 + g_xxxxz_yyzzzz_1[i] * pa_y[i];

        g_xxxxyz_yzzzzz_0[i] = g_xxxxz_zzzzz_1[i] * fe_0 + g_xxxxz_yzzzzz_1[i] * pa_y[i];

        g_xxxxyz_zzzzzz_0[i] = g_xxxxz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : II

    auto g_xxxxzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 140);

    auto g_xxxxzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 141);

    auto g_xxxxzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 142);

    auto g_xxxxzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 143);

    auto g_xxxxzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 144);

    auto g_xxxxzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 145);

    auto g_xxxxzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 146);

    auto g_xxxxzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 147);

    auto g_xxxxzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 148);

    auto g_xxxxzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 149);

    auto g_xxxxzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 150);

    auto g_xxxxzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 151);

    auto g_xxxxzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 152);

    auto g_xxxxzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 153);

    auto g_xxxxzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 154);

    auto g_xxxxzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 155);

    auto g_xxxxzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 156);

    auto g_xxxxzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 157);

    auto g_xxxxzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 158);

    auto g_xxxxzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 159);

    auto g_xxxxzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 160);

    auto g_xxxxzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 161);

    auto g_xxxxzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 162);

    auto g_xxxxzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 163);

    auto g_xxxxzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 164);

    auto g_xxxxzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 165);

    auto g_xxxxzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 166);

    auto g_xxxxzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 167);

    #pragma omp simd aligned(g_xxxx_xxxxxx_0, g_xxxx_xxxxxx_1, g_xxxx_xxxxxy_0, g_xxxx_xxxxxy_1, g_xxxx_xxxxyy_0, g_xxxx_xxxxyy_1, g_xxxx_xxxyyy_0, g_xxxx_xxxyyy_1, g_xxxx_xxyyyy_0, g_xxxx_xxyyyy_1, g_xxxx_xyyyyy_0, g_xxxx_xyyyyy_1, g_xxxxz_xxxxxx_1, g_xxxxz_xxxxxy_1, g_xxxxz_xxxxyy_1, g_xxxxz_xxxyyy_1, g_xxxxz_xxyyyy_1, g_xxxxz_xyyyyy_1, g_xxxxzz_xxxxxx_0, g_xxxxzz_xxxxxy_0, g_xxxxzz_xxxxxz_0, g_xxxxzz_xxxxyy_0, g_xxxxzz_xxxxyz_0, g_xxxxzz_xxxxzz_0, g_xxxxzz_xxxyyy_0, g_xxxxzz_xxxyyz_0, g_xxxxzz_xxxyzz_0, g_xxxxzz_xxxzzz_0, g_xxxxzz_xxyyyy_0, g_xxxxzz_xxyyyz_0, g_xxxxzz_xxyyzz_0, g_xxxxzz_xxyzzz_0, g_xxxxzz_xxzzzz_0, g_xxxxzz_xyyyyy_0, g_xxxxzz_xyyyyz_0, g_xxxxzz_xyyyzz_0, g_xxxxzz_xyyzzz_0, g_xxxxzz_xyzzzz_0, g_xxxxzz_xzzzzz_0, g_xxxxzz_yyyyyy_0, g_xxxxzz_yyyyyz_0, g_xxxxzz_yyyyzz_0, g_xxxxzz_yyyzzz_0, g_xxxxzz_yyzzzz_0, g_xxxxzz_yzzzzz_0, g_xxxxzz_zzzzzz_0, g_xxxzz_xxxxxz_1, g_xxxzz_xxxxyz_1, g_xxxzz_xxxxz_1, g_xxxzz_xxxxzz_1, g_xxxzz_xxxyyz_1, g_xxxzz_xxxyz_1, g_xxxzz_xxxyzz_1, g_xxxzz_xxxzz_1, g_xxxzz_xxxzzz_1, g_xxxzz_xxyyyz_1, g_xxxzz_xxyyz_1, g_xxxzz_xxyyzz_1, g_xxxzz_xxyzz_1, g_xxxzz_xxyzzz_1, g_xxxzz_xxzzz_1, g_xxxzz_xxzzzz_1, g_xxxzz_xyyyyz_1, g_xxxzz_xyyyz_1, g_xxxzz_xyyyzz_1, g_xxxzz_xyyzz_1, g_xxxzz_xyyzzz_1, g_xxxzz_xyzzz_1, g_xxxzz_xyzzzz_1, g_xxxzz_xzzzz_1, g_xxxzz_xzzzzz_1, g_xxxzz_yyyyyy_1, g_xxxzz_yyyyyz_1, g_xxxzz_yyyyz_1, g_xxxzz_yyyyzz_1, g_xxxzz_yyyzz_1, g_xxxzz_yyyzzz_1, g_xxxzz_yyzzz_1, g_xxxzz_yyzzzz_1, g_xxxzz_yzzzz_1, g_xxxzz_yzzzzz_1, g_xxxzz_zzzzz_1, g_xxxzz_zzzzzz_1, g_xxzz_xxxxxz_0, g_xxzz_xxxxxz_1, g_xxzz_xxxxyz_0, g_xxzz_xxxxyz_1, g_xxzz_xxxxzz_0, g_xxzz_xxxxzz_1, g_xxzz_xxxyyz_0, g_xxzz_xxxyyz_1, g_xxzz_xxxyzz_0, g_xxzz_xxxyzz_1, g_xxzz_xxxzzz_0, g_xxzz_xxxzzz_1, g_xxzz_xxyyyz_0, g_xxzz_xxyyyz_1, g_xxzz_xxyyzz_0, g_xxzz_xxyyzz_1, g_xxzz_xxyzzz_0, g_xxzz_xxyzzz_1, g_xxzz_xxzzzz_0, g_xxzz_xxzzzz_1, g_xxzz_xyyyyz_0, g_xxzz_xyyyyz_1, g_xxzz_xyyyzz_0, g_xxzz_xyyyzz_1, g_xxzz_xyyzzz_0, g_xxzz_xyyzzz_1, g_xxzz_xyzzzz_0, g_xxzz_xyzzzz_1, g_xxzz_xzzzzz_0, g_xxzz_xzzzzz_1, g_xxzz_yyyyyy_0, g_xxzz_yyyyyy_1, g_xxzz_yyyyyz_0, g_xxzz_yyyyyz_1, g_xxzz_yyyyzz_0, g_xxzz_yyyyzz_1, g_xxzz_yyyzzz_0, g_xxzz_yyyzzz_1, g_xxzz_yyzzzz_0, g_xxzz_yyzzzz_1, g_xxzz_yzzzzz_0, g_xxzz_yzzzzz_1, g_xxzz_zzzzzz_0, g_xxzz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxzz_xxxxxx_0[i] = g_xxxx_xxxxxx_0[i] * fbe_0 - g_xxxx_xxxxxx_1[i] * fz_be_0 + g_xxxxz_xxxxxx_1[i] * pa_z[i];

        g_xxxxzz_xxxxxy_0[i] = g_xxxx_xxxxxy_0[i] * fbe_0 - g_xxxx_xxxxxy_1[i] * fz_be_0 + g_xxxxz_xxxxxy_1[i] * pa_z[i];

        g_xxxxzz_xxxxxz_0[i] = 3.0 * g_xxzz_xxxxxz_0[i] * fbe_0 - 3.0 * g_xxzz_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxxzz_xxxxz_1[i] * fe_0 + g_xxxzz_xxxxxz_1[i] * pa_x[i];

        g_xxxxzz_xxxxyy_0[i] = g_xxxx_xxxxyy_0[i] * fbe_0 - g_xxxx_xxxxyy_1[i] * fz_be_0 + g_xxxxz_xxxxyy_1[i] * pa_z[i];

        g_xxxxzz_xxxxyz_0[i] = 3.0 * g_xxzz_xxxxyz_0[i] * fbe_0 - 3.0 * g_xxzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxzz_xxxyz_1[i] * fe_0 + g_xxxzz_xxxxyz_1[i] * pa_x[i];

        g_xxxxzz_xxxxzz_0[i] = 3.0 * g_xxzz_xxxxzz_0[i] * fbe_0 - 3.0 * g_xxzz_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxxzz_xxxzz_1[i] * fe_0 + g_xxxzz_xxxxzz_1[i] * pa_x[i];

        g_xxxxzz_xxxyyy_0[i] = g_xxxx_xxxyyy_0[i] * fbe_0 - g_xxxx_xxxyyy_1[i] * fz_be_0 + g_xxxxz_xxxyyy_1[i] * pa_z[i];

        g_xxxxzz_xxxyyz_0[i] = 3.0 * g_xxzz_xxxyyz_0[i] * fbe_0 - 3.0 * g_xxzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxzz_xxyyz_1[i] * fe_0 + g_xxxzz_xxxyyz_1[i] * pa_x[i];

        g_xxxxzz_xxxyzz_0[i] = 3.0 * g_xxzz_xxxyzz_0[i] * fbe_0 - 3.0 * g_xxzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_xxyzz_1[i] * fe_0 + g_xxxzz_xxxyzz_1[i] * pa_x[i];

        g_xxxxzz_xxxzzz_0[i] = 3.0 * g_xxzz_xxxzzz_0[i] * fbe_0 - 3.0 * g_xxzz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_xxzzz_1[i] * fe_0 + g_xxxzz_xxxzzz_1[i] * pa_x[i];

        g_xxxxzz_xxyyyy_0[i] = g_xxxx_xxyyyy_0[i] * fbe_0 - g_xxxx_xxyyyy_1[i] * fz_be_0 + g_xxxxz_xxyyyy_1[i] * pa_z[i];

        g_xxxxzz_xxyyyz_0[i] = 3.0 * g_xxzz_xxyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxzz_xyyyz_1[i] * fe_0 + g_xxxzz_xxyyyz_1[i] * pa_x[i];

        g_xxxxzz_xxyyzz_0[i] = 3.0 * g_xxzz_xxyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_xyyzz_1[i] * fe_0 + g_xxxzz_xxyyzz_1[i] * pa_x[i];

        g_xxxxzz_xxyzzz_0[i] = 3.0 * g_xxzz_xxyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_xyzzz_1[i] * fe_0 + g_xxxzz_xxyzzz_1[i] * pa_x[i];

        g_xxxxzz_xxzzzz_0[i] = 3.0 * g_xxzz_xxzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_xzzzz_1[i] * fe_0 + g_xxxzz_xxzzzz_1[i] * pa_x[i];

        g_xxxxzz_xyyyyy_0[i] = g_xxxx_xyyyyy_0[i] * fbe_0 - g_xxxx_xyyyyy_1[i] * fz_be_0 + g_xxxxz_xyyyyy_1[i] * pa_z[i];

        g_xxxxzz_xyyyyz_0[i] = 3.0 * g_xxzz_xyyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_xyyyyz_1[i] * fz_be_0 + g_xxxzz_yyyyz_1[i] * fe_0 + g_xxxzz_xyyyyz_1[i] * pa_x[i];

        g_xxxxzz_xyyyzz_0[i] = 3.0 * g_xxzz_xyyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_xyyyzz_1[i] * fz_be_0 + g_xxxzz_yyyzz_1[i] * fe_0 + g_xxxzz_xyyyzz_1[i] * pa_x[i];

        g_xxxxzz_xyyzzz_0[i] = 3.0 * g_xxzz_xyyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_xyyzzz_1[i] * fz_be_0 + g_xxxzz_yyzzz_1[i] * fe_0 + g_xxxzz_xyyzzz_1[i] * pa_x[i];

        g_xxxxzz_xyzzzz_0[i] = 3.0 * g_xxzz_xyzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_xyzzzz_1[i] * fz_be_0 + g_xxxzz_yzzzz_1[i] * fe_0 + g_xxxzz_xyzzzz_1[i] * pa_x[i];

        g_xxxxzz_xzzzzz_0[i] = 3.0 * g_xxzz_xzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_xzzzzz_1[i] * fz_be_0 + g_xxxzz_zzzzz_1[i] * fe_0 + g_xxxzz_xzzzzz_1[i] * pa_x[i];

        g_xxxxzz_yyyyyy_0[i] = 3.0 * g_xxzz_yyyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_yyyyyy_1[i] * fz_be_0 + g_xxxzz_yyyyyy_1[i] * pa_x[i];

        g_xxxxzz_yyyyyz_0[i] = 3.0 * g_xxzz_yyyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_yyyyyz_1[i] * fz_be_0 + g_xxxzz_yyyyyz_1[i] * pa_x[i];

        g_xxxxzz_yyyyzz_0[i] = 3.0 * g_xxzz_yyyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_yyyyzz_1[i] * fz_be_0 + g_xxxzz_yyyyzz_1[i] * pa_x[i];

        g_xxxxzz_yyyzzz_0[i] = 3.0 * g_xxzz_yyyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_yyyzzz_1[i] * fz_be_0 + g_xxxzz_yyyzzz_1[i] * pa_x[i];

        g_xxxxzz_yyzzzz_0[i] = 3.0 * g_xxzz_yyzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_yyzzzz_1[i] * fz_be_0 + g_xxxzz_yyzzzz_1[i] * pa_x[i];

        g_xxxxzz_yzzzzz_0[i] = 3.0 * g_xxzz_yzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_yzzzzz_1[i] * fz_be_0 + g_xxxzz_yzzzzz_1[i] * pa_x[i];

        g_xxxxzz_zzzzzz_0[i] = 3.0 * g_xxzz_zzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_zzzzzz_1[i] * fz_be_0 + g_xxxzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : II

    auto g_xxxyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 168);

    auto g_xxxyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 169);

    auto g_xxxyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 170);

    auto g_xxxyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 171);

    auto g_xxxyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 172);

    auto g_xxxyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 173);

    auto g_xxxyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 174);

    auto g_xxxyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 175);

    auto g_xxxyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 176);

    auto g_xxxyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 177);

    auto g_xxxyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 178);

    auto g_xxxyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 179);

    auto g_xxxyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 180);

    auto g_xxxyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 181);

    auto g_xxxyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 182);

    auto g_xxxyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 183);

    auto g_xxxyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 184);

    auto g_xxxyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 185);

    auto g_xxxyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 186);

    auto g_xxxyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 187);

    auto g_xxxyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 188);

    auto g_xxxyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 189);

    auto g_xxxyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 190);

    auto g_xxxyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 191);

    auto g_xxxyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 192);

    auto g_xxxyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 193);

    auto g_xxxyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 194);

    auto g_xxxyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 195);

    #pragma omp simd aligned(g_xxxy_xxxxxx_0, g_xxxy_xxxxxx_1, g_xxxy_xxxxxz_0, g_xxxy_xxxxxz_1, g_xxxy_xxxxzz_0, g_xxxy_xxxxzz_1, g_xxxy_xxxzzz_0, g_xxxy_xxxzzz_1, g_xxxy_xxzzzz_0, g_xxxy_xxzzzz_1, g_xxxy_xzzzzz_0, g_xxxy_xzzzzz_1, g_xxxyy_xxxxxx_1, g_xxxyy_xxxxxz_1, g_xxxyy_xxxxzz_1, g_xxxyy_xxxzzz_1, g_xxxyy_xxzzzz_1, g_xxxyy_xzzzzz_1, g_xxxyyy_xxxxxx_0, g_xxxyyy_xxxxxy_0, g_xxxyyy_xxxxxz_0, g_xxxyyy_xxxxyy_0, g_xxxyyy_xxxxyz_0, g_xxxyyy_xxxxzz_0, g_xxxyyy_xxxyyy_0, g_xxxyyy_xxxyyz_0, g_xxxyyy_xxxyzz_0, g_xxxyyy_xxxzzz_0, g_xxxyyy_xxyyyy_0, g_xxxyyy_xxyyyz_0, g_xxxyyy_xxyyzz_0, g_xxxyyy_xxyzzz_0, g_xxxyyy_xxzzzz_0, g_xxxyyy_xyyyyy_0, g_xxxyyy_xyyyyz_0, g_xxxyyy_xyyyzz_0, g_xxxyyy_xyyzzz_0, g_xxxyyy_xyzzzz_0, g_xxxyyy_xzzzzz_0, g_xxxyyy_yyyyyy_0, g_xxxyyy_yyyyyz_0, g_xxxyyy_yyyyzz_0, g_xxxyyy_yyyzzz_0, g_xxxyyy_yyzzzz_0, g_xxxyyy_yzzzzz_0, g_xxxyyy_zzzzzz_0, g_xxyyy_xxxxxy_1, g_xxyyy_xxxxy_1, g_xxyyy_xxxxyy_1, g_xxyyy_xxxxyz_1, g_xxyyy_xxxyy_1, g_xxyyy_xxxyyy_1, g_xxyyy_xxxyyz_1, g_xxyyy_xxxyz_1, g_xxyyy_xxxyzz_1, g_xxyyy_xxyyy_1, g_xxyyy_xxyyyy_1, g_xxyyy_xxyyyz_1, g_xxyyy_xxyyz_1, g_xxyyy_xxyyzz_1, g_xxyyy_xxyzz_1, g_xxyyy_xxyzzz_1, g_xxyyy_xyyyy_1, g_xxyyy_xyyyyy_1, g_xxyyy_xyyyyz_1, g_xxyyy_xyyyz_1, g_xxyyy_xyyyzz_1, g_xxyyy_xyyzz_1, g_xxyyy_xyyzzz_1, g_xxyyy_xyzzz_1, g_xxyyy_xyzzzz_1, g_xxyyy_yyyyy_1, g_xxyyy_yyyyyy_1, g_xxyyy_yyyyyz_1, g_xxyyy_yyyyz_1, g_xxyyy_yyyyzz_1, g_xxyyy_yyyzz_1, g_xxyyy_yyyzzz_1, g_xxyyy_yyzzz_1, g_xxyyy_yyzzzz_1, g_xxyyy_yzzzz_1, g_xxyyy_yzzzzz_1, g_xxyyy_zzzzzz_1, g_xyyy_xxxxxy_0, g_xyyy_xxxxxy_1, g_xyyy_xxxxyy_0, g_xyyy_xxxxyy_1, g_xyyy_xxxxyz_0, g_xyyy_xxxxyz_1, g_xyyy_xxxyyy_0, g_xyyy_xxxyyy_1, g_xyyy_xxxyyz_0, g_xyyy_xxxyyz_1, g_xyyy_xxxyzz_0, g_xyyy_xxxyzz_1, g_xyyy_xxyyyy_0, g_xyyy_xxyyyy_1, g_xyyy_xxyyyz_0, g_xyyy_xxyyyz_1, g_xyyy_xxyyzz_0, g_xyyy_xxyyzz_1, g_xyyy_xxyzzz_0, g_xyyy_xxyzzz_1, g_xyyy_xyyyyy_0, g_xyyy_xyyyyy_1, g_xyyy_xyyyyz_0, g_xyyy_xyyyyz_1, g_xyyy_xyyyzz_0, g_xyyy_xyyyzz_1, g_xyyy_xyyzzz_0, g_xyyy_xyyzzz_1, g_xyyy_xyzzzz_0, g_xyyy_xyzzzz_1, g_xyyy_yyyyyy_0, g_xyyy_yyyyyy_1, g_xyyy_yyyyyz_0, g_xyyy_yyyyyz_1, g_xyyy_yyyyzz_0, g_xyyy_yyyyzz_1, g_xyyy_yyyzzz_0, g_xyyy_yyyzzz_1, g_xyyy_yyzzzz_0, g_xyyy_yyzzzz_1, g_xyyy_yzzzzz_0, g_xyyy_yzzzzz_1, g_xyyy_zzzzzz_0, g_xyyy_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyy_xxxxxx_0[i] = 2.0 * g_xxxy_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxxy_xxxxxx_1[i] * fz_be_0 + g_xxxyy_xxxxxx_1[i] * pa_y[i];

        g_xxxyyy_xxxxxy_0[i] = 2.0 * g_xyyy_xxxxxy_0[i] * fbe_0 - 2.0 * g_xyyy_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxyyy_xxxxy_1[i] * fe_0 + g_xxyyy_xxxxxy_1[i] * pa_x[i];

        g_xxxyyy_xxxxxz_0[i] = 2.0 * g_xxxy_xxxxxz_0[i] * fbe_0 - 2.0 * g_xxxy_xxxxxz_1[i] * fz_be_0 + g_xxxyy_xxxxxz_1[i] * pa_y[i];

        g_xxxyyy_xxxxyy_0[i] = 2.0 * g_xyyy_xxxxyy_0[i] * fbe_0 - 2.0 * g_xyyy_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxyyy_xxxyy_1[i] * fe_0 + g_xxyyy_xxxxyy_1[i] * pa_x[i];

        g_xxxyyy_xxxxyz_0[i] = 2.0 * g_xyyy_xxxxyz_0[i] * fbe_0 - 2.0 * g_xyyy_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxyyy_xxxyz_1[i] * fe_0 + g_xxyyy_xxxxyz_1[i] * pa_x[i];

        g_xxxyyy_xxxxzz_0[i] = 2.0 * g_xxxy_xxxxzz_0[i] * fbe_0 - 2.0 * g_xxxy_xxxxzz_1[i] * fz_be_0 + g_xxxyy_xxxxzz_1[i] * pa_y[i];

        g_xxxyyy_xxxyyy_0[i] = 2.0 * g_xyyy_xxxyyy_0[i] * fbe_0 - 2.0 * g_xyyy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxyyy_xxyyy_1[i] * fe_0 + g_xxyyy_xxxyyy_1[i] * pa_x[i];

        g_xxxyyy_xxxyyz_0[i] = 2.0 * g_xyyy_xxxyyz_0[i] * fbe_0 - 2.0 * g_xyyy_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxyyy_xxyyz_1[i] * fe_0 + g_xxyyy_xxxyyz_1[i] * pa_x[i];

        g_xxxyyy_xxxyzz_0[i] = 2.0 * g_xyyy_xxxyzz_0[i] * fbe_0 - 2.0 * g_xyyy_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxyyy_xxyzz_1[i] * fe_0 + g_xxyyy_xxxyzz_1[i] * pa_x[i];

        g_xxxyyy_xxxzzz_0[i] = 2.0 * g_xxxy_xxxzzz_0[i] * fbe_0 - 2.0 * g_xxxy_xxxzzz_1[i] * fz_be_0 + g_xxxyy_xxxzzz_1[i] * pa_y[i];

        g_xxxyyy_xxyyyy_0[i] = 2.0 * g_xyyy_xxyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxyyy_xyyyy_1[i] * fe_0 + g_xxyyy_xxyyyy_1[i] * pa_x[i];

        g_xxxyyy_xxyyyz_0[i] = 2.0 * g_xyyy_xxyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxyyy_xyyyz_1[i] * fe_0 + g_xxyyy_xxyyyz_1[i] * pa_x[i];

        g_xxxyyy_xxyyzz_0[i] = 2.0 * g_xyyy_xxyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_xyyzz_1[i] * fe_0 + g_xxyyy_xxyyzz_1[i] * pa_x[i];

        g_xxxyyy_xxyzzz_0[i] = 2.0 * g_xyyy_xxyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_xyzzz_1[i] * fe_0 + g_xxyyy_xxyzzz_1[i] * pa_x[i];

        g_xxxyyy_xxzzzz_0[i] = 2.0 * g_xxxy_xxzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_xxzzzz_1[i] * fz_be_0 + g_xxxyy_xxzzzz_1[i] * pa_y[i];

        g_xxxyyy_xyyyyy_0[i] = 2.0 * g_xyyy_xyyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_xyyyyy_1[i] * fz_be_0 + g_xxyyy_yyyyy_1[i] * fe_0 + g_xxyyy_xyyyyy_1[i] * pa_x[i];

        g_xxxyyy_xyyyyz_0[i] = 2.0 * g_xyyy_xyyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_xyyyyz_1[i] * fz_be_0 + g_xxyyy_yyyyz_1[i] * fe_0 + g_xxyyy_xyyyyz_1[i] * pa_x[i];

        g_xxxyyy_xyyyzz_0[i] = 2.0 * g_xyyy_xyyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_xyyyzz_1[i] * fz_be_0 + g_xxyyy_yyyzz_1[i] * fe_0 + g_xxyyy_xyyyzz_1[i] * pa_x[i];

        g_xxxyyy_xyyzzz_0[i] = 2.0 * g_xyyy_xyyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_xyyzzz_1[i] * fz_be_0 + g_xxyyy_yyzzz_1[i] * fe_0 + g_xxyyy_xyyzzz_1[i] * pa_x[i];

        g_xxxyyy_xyzzzz_0[i] = 2.0 * g_xyyy_xyzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_xyzzzz_1[i] * fz_be_0 + g_xxyyy_yzzzz_1[i] * fe_0 + g_xxyyy_xyzzzz_1[i] * pa_x[i];

        g_xxxyyy_xzzzzz_0[i] = 2.0 * g_xxxy_xzzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_xzzzzz_1[i] * fz_be_0 + g_xxxyy_xzzzzz_1[i] * pa_y[i];

        g_xxxyyy_yyyyyy_0[i] = 2.0 * g_xyyy_yyyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_yyyyyy_1[i] * fz_be_0 + g_xxyyy_yyyyyy_1[i] * pa_x[i];

        g_xxxyyy_yyyyyz_0[i] = 2.0 * g_xyyy_yyyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_yyyyyz_1[i] * fz_be_0 + g_xxyyy_yyyyyz_1[i] * pa_x[i];

        g_xxxyyy_yyyyzz_0[i] = 2.0 * g_xyyy_yyyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_yyyyzz_1[i] * fz_be_0 + g_xxyyy_yyyyzz_1[i] * pa_x[i];

        g_xxxyyy_yyyzzz_0[i] = 2.0 * g_xyyy_yyyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_yyyzzz_1[i] * fz_be_0 + g_xxyyy_yyyzzz_1[i] * pa_x[i];

        g_xxxyyy_yyzzzz_0[i] = 2.0 * g_xyyy_yyzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_yyzzzz_1[i] * fz_be_0 + g_xxyyy_yyzzzz_1[i] * pa_x[i];

        g_xxxyyy_yzzzzz_0[i] = 2.0 * g_xyyy_yzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_yzzzzz_1[i] * fz_be_0 + g_xxyyy_yzzzzz_1[i] * pa_x[i];

        g_xxxyyy_zzzzzz_0[i] = 2.0 * g_xyyy_zzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_zzzzzz_1[i] * fz_be_0 + g_xxyyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 196-224 components of targeted buffer : II

    auto g_xxxyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 196);

    auto g_xxxyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 197);

    auto g_xxxyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 198);

    auto g_xxxyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 199);

    auto g_xxxyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 200);

    auto g_xxxyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 201);

    auto g_xxxyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 202);

    auto g_xxxyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 203);

    auto g_xxxyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 204);

    auto g_xxxyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 205);

    auto g_xxxyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 206);

    auto g_xxxyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 207);

    auto g_xxxyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 208);

    auto g_xxxyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 209);

    auto g_xxxyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 210);

    auto g_xxxyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 211);

    auto g_xxxyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 212);

    auto g_xxxyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 213);

    auto g_xxxyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 214);

    auto g_xxxyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 215);

    auto g_xxxyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 216);

    auto g_xxxyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 217);

    auto g_xxxyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 218);

    auto g_xxxyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 219);

    auto g_xxxyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 220);

    auto g_xxxyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 221);

    auto g_xxxyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 222);

    auto g_xxxyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 223);

    #pragma omp simd aligned(g_xxxyy_xxxxx_1, g_xxxyy_xxxxxx_1, g_xxxyy_xxxxxy_1, g_xxxyy_xxxxxz_1, g_xxxyy_xxxxy_1, g_xxxyy_xxxxyy_1, g_xxxyy_xxxxyz_1, g_xxxyy_xxxxz_1, g_xxxyy_xxxxzz_1, g_xxxyy_xxxyy_1, g_xxxyy_xxxyyy_1, g_xxxyy_xxxyyz_1, g_xxxyy_xxxyz_1, g_xxxyy_xxxyzz_1, g_xxxyy_xxxzz_1, g_xxxyy_xxxzzz_1, g_xxxyy_xxyyy_1, g_xxxyy_xxyyyy_1, g_xxxyy_xxyyyz_1, g_xxxyy_xxyyz_1, g_xxxyy_xxyyzz_1, g_xxxyy_xxyzz_1, g_xxxyy_xxyzzz_1, g_xxxyy_xxzzz_1, g_xxxyy_xxzzzz_1, g_xxxyy_xyyyy_1, g_xxxyy_xyyyyy_1, g_xxxyy_xyyyyz_1, g_xxxyy_xyyyz_1, g_xxxyy_xyyyzz_1, g_xxxyy_xyyzz_1, g_xxxyy_xyyzzz_1, g_xxxyy_xyzzz_1, g_xxxyy_xyzzzz_1, g_xxxyy_xzzzz_1, g_xxxyy_xzzzzz_1, g_xxxyy_yyyyy_1, g_xxxyy_yyyyyy_1, g_xxxyy_yyyyyz_1, g_xxxyy_yyyyz_1, g_xxxyy_yyyyzz_1, g_xxxyy_yyyzz_1, g_xxxyy_yyyzzz_1, g_xxxyy_yyzzz_1, g_xxxyy_yyzzzz_1, g_xxxyy_yzzzz_1, g_xxxyy_yzzzzz_1, g_xxxyy_zzzzz_1, g_xxxyy_zzzzzz_1, g_xxxyyz_xxxxxx_0, g_xxxyyz_xxxxxy_0, g_xxxyyz_xxxxxz_0, g_xxxyyz_xxxxyy_0, g_xxxyyz_xxxxyz_0, g_xxxyyz_xxxxzz_0, g_xxxyyz_xxxyyy_0, g_xxxyyz_xxxyyz_0, g_xxxyyz_xxxyzz_0, g_xxxyyz_xxxzzz_0, g_xxxyyz_xxyyyy_0, g_xxxyyz_xxyyyz_0, g_xxxyyz_xxyyzz_0, g_xxxyyz_xxyzzz_0, g_xxxyyz_xxzzzz_0, g_xxxyyz_xyyyyy_0, g_xxxyyz_xyyyyz_0, g_xxxyyz_xyyyzz_0, g_xxxyyz_xyyzzz_0, g_xxxyyz_xyzzzz_0, g_xxxyyz_xzzzzz_0, g_xxxyyz_yyyyyy_0, g_xxxyyz_yyyyyz_0, g_xxxyyz_yyyyzz_0, g_xxxyyz_yyyzzz_0, g_xxxyyz_yyzzzz_0, g_xxxyyz_yzzzzz_0, g_xxxyyz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyyz_xxxxxx_0[i] = g_xxxyy_xxxxxx_1[i] * pa_z[i];

        g_xxxyyz_xxxxxy_0[i] = g_xxxyy_xxxxxy_1[i] * pa_z[i];

        g_xxxyyz_xxxxxz_0[i] = g_xxxyy_xxxxx_1[i] * fe_0 + g_xxxyy_xxxxxz_1[i] * pa_z[i];

        g_xxxyyz_xxxxyy_0[i] = g_xxxyy_xxxxyy_1[i] * pa_z[i];

        g_xxxyyz_xxxxyz_0[i] = g_xxxyy_xxxxy_1[i] * fe_0 + g_xxxyy_xxxxyz_1[i] * pa_z[i];

        g_xxxyyz_xxxxzz_0[i] = 2.0 * g_xxxyy_xxxxz_1[i] * fe_0 + g_xxxyy_xxxxzz_1[i] * pa_z[i];

        g_xxxyyz_xxxyyy_0[i] = g_xxxyy_xxxyyy_1[i] * pa_z[i];

        g_xxxyyz_xxxyyz_0[i] = g_xxxyy_xxxyy_1[i] * fe_0 + g_xxxyy_xxxyyz_1[i] * pa_z[i];

        g_xxxyyz_xxxyzz_0[i] = 2.0 * g_xxxyy_xxxyz_1[i] * fe_0 + g_xxxyy_xxxyzz_1[i] * pa_z[i];

        g_xxxyyz_xxxzzz_0[i] = 3.0 * g_xxxyy_xxxzz_1[i] * fe_0 + g_xxxyy_xxxzzz_1[i] * pa_z[i];

        g_xxxyyz_xxyyyy_0[i] = g_xxxyy_xxyyyy_1[i] * pa_z[i];

        g_xxxyyz_xxyyyz_0[i] = g_xxxyy_xxyyy_1[i] * fe_0 + g_xxxyy_xxyyyz_1[i] * pa_z[i];

        g_xxxyyz_xxyyzz_0[i] = 2.0 * g_xxxyy_xxyyz_1[i] * fe_0 + g_xxxyy_xxyyzz_1[i] * pa_z[i];

        g_xxxyyz_xxyzzz_0[i] = 3.0 * g_xxxyy_xxyzz_1[i] * fe_0 + g_xxxyy_xxyzzz_1[i] * pa_z[i];

        g_xxxyyz_xxzzzz_0[i] = 4.0 * g_xxxyy_xxzzz_1[i] * fe_0 + g_xxxyy_xxzzzz_1[i] * pa_z[i];

        g_xxxyyz_xyyyyy_0[i] = g_xxxyy_xyyyyy_1[i] * pa_z[i];

        g_xxxyyz_xyyyyz_0[i] = g_xxxyy_xyyyy_1[i] * fe_0 + g_xxxyy_xyyyyz_1[i] * pa_z[i];

        g_xxxyyz_xyyyzz_0[i] = 2.0 * g_xxxyy_xyyyz_1[i] * fe_0 + g_xxxyy_xyyyzz_1[i] * pa_z[i];

        g_xxxyyz_xyyzzz_0[i] = 3.0 * g_xxxyy_xyyzz_1[i] * fe_0 + g_xxxyy_xyyzzz_1[i] * pa_z[i];

        g_xxxyyz_xyzzzz_0[i] = 4.0 * g_xxxyy_xyzzz_1[i] * fe_0 + g_xxxyy_xyzzzz_1[i] * pa_z[i];

        g_xxxyyz_xzzzzz_0[i] = 5.0 * g_xxxyy_xzzzz_1[i] * fe_0 + g_xxxyy_xzzzzz_1[i] * pa_z[i];

        g_xxxyyz_yyyyyy_0[i] = g_xxxyy_yyyyyy_1[i] * pa_z[i];

        g_xxxyyz_yyyyyz_0[i] = g_xxxyy_yyyyy_1[i] * fe_0 + g_xxxyy_yyyyyz_1[i] * pa_z[i];

        g_xxxyyz_yyyyzz_0[i] = 2.0 * g_xxxyy_yyyyz_1[i] * fe_0 + g_xxxyy_yyyyzz_1[i] * pa_z[i];

        g_xxxyyz_yyyzzz_0[i] = 3.0 * g_xxxyy_yyyzz_1[i] * fe_0 + g_xxxyy_yyyzzz_1[i] * pa_z[i];

        g_xxxyyz_yyzzzz_0[i] = 4.0 * g_xxxyy_yyzzz_1[i] * fe_0 + g_xxxyy_yyzzzz_1[i] * pa_z[i];

        g_xxxyyz_yzzzzz_0[i] = 5.0 * g_xxxyy_yzzzz_1[i] * fe_0 + g_xxxyy_yzzzzz_1[i] * pa_z[i];

        g_xxxyyz_zzzzzz_0[i] = 6.0 * g_xxxyy_zzzzz_1[i] * fe_0 + g_xxxyy_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 224-252 components of targeted buffer : II

    auto g_xxxyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 224);

    auto g_xxxyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 225);

    auto g_xxxyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 226);

    auto g_xxxyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 227);

    auto g_xxxyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 228);

    auto g_xxxyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 229);

    auto g_xxxyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 230);

    auto g_xxxyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 231);

    auto g_xxxyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 232);

    auto g_xxxyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 233);

    auto g_xxxyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 234);

    auto g_xxxyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 235);

    auto g_xxxyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 236);

    auto g_xxxyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 237);

    auto g_xxxyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 238);

    auto g_xxxyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 239);

    auto g_xxxyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 240);

    auto g_xxxyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 241);

    auto g_xxxyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 242);

    auto g_xxxyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 243);

    auto g_xxxyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 244);

    auto g_xxxyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 245);

    auto g_xxxyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 246);

    auto g_xxxyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 247);

    auto g_xxxyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 248);

    auto g_xxxyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 249);

    auto g_xxxyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 250);

    auto g_xxxyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 251);

    #pragma omp simd aligned(g_xxxyzz_xxxxxx_0, g_xxxyzz_xxxxxy_0, g_xxxyzz_xxxxxz_0, g_xxxyzz_xxxxyy_0, g_xxxyzz_xxxxyz_0, g_xxxyzz_xxxxzz_0, g_xxxyzz_xxxyyy_0, g_xxxyzz_xxxyyz_0, g_xxxyzz_xxxyzz_0, g_xxxyzz_xxxzzz_0, g_xxxyzz_xxyyyy_0, g_xxxyzz_xxyyyz_0, g_xxxyzz_xxyyzz_0, g_xxxyzz_xxyzzz_0, g_xxxyzz_xxzzzz_0, g_xxxyzz_xyyyyy_0, g_xxxyzz_xyyyyz_0, g_xxxyzz_xyyyzz_0, g_xxxyzz_xyyzzz_0, g_xxxyzz_xyzzzz_0, g_xxxyzz_xzzzzz_0, g_xxxyzz_yyyyyy_0, g_xxxyzz_yyyyyz_0, g_xxxyzz_yyyyzz_0, g_xxxyzz_yyyzzz_0, g_xxxyzz_yyzzzz_0, g_xxxyzz_yzzzzz_0, g_xxxyzz_zzzzzz_0, g_xxxzz_xxxxx_1, g_xxxzz_xxxxxx_1, g_xxxzz_xxxxxy_1, g_xxxzz_xxxxxz_1, g_xxxzz_xxxxy_1, g_xxxzz_xxxxyy_1, g_xxxzz_xxxxyz_1, g_xxxzz_xxxxz_1, g_xxxzz_xxxxzz_1, g_xxxzz_xxxyy_1, g_xxxzz_xxxyyy_1, g_xxxzz_xxxyyz_1, g_xxxzz_xxxyz_1, g_xxxzz_xxxyzz_1, g_xxxzz_xxxzz_1, g_xxxzz_xxxzzz_1, g_xxxzz_xxyyy_1, g_xxxzz_xxyyyy_1, g_xxxzz_xxyyyz_1, g_xxxzz_xxyyz_1, g_xxxzz_xxyyzz_1, g_xxxzz_xxyzz_1, g_xxxzz_xxyzzz_1, g_xxxzz_xxzzz_1, g_xxxzz_xxzzzz_1, g_xxxzz_xyyyy_1, g_xxxzz_xyyyyy_1, g_xxxzz_xyyyyz_1, g_xxxzz_xyyyz_1, g_xxxzz_xyyyzz_1, g_xxxzz_xyyzz_1, g_xxxzz_xyyzzz_1, g_xxxzz_xyzzz_1, g_xxxzz_xyzzzz_1, g_xxxzz_xzzzz_1, g_xxxzz_xzzzzz_1, g_xxxzz_yyyyy_1, g_xxxzz_yyyyyy_1, g_xxxzz_yyyyyz_1, g_xxxzz_yyyyz_1, g_xxxzz_yyyyzz_1, g_xxxzz_yyyzz_1, g_xxxzz_yyyzzz_1, g_xxxzz_yyzzz_1, g_xxxzz_yyzzzz_1, g_xxxzz_yzzzz_1, g_xxxzz_yzzzzz_1, g_xxxzz_zzzzz_1, g_xxxzz_zzzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyzz_xxxxxx_0[i] = g_xxxzz_xxxxxx_1[i] * pa_y[i];

        g_xxxyzz_xxxxxy_0[i] = g_xxxzz_xxxxx_1[i] * fe_0 + g_xxxzz_xxxxxy_1[i] * pa_y[i];

        g_xxxyzz_xxxxxz_0[i] = g_xxxzz_xxxxxz_1[i] * pa_y[i];

        g_xxxyzz_xxxxyy_0[i] = 2.0 * g_xxxzz_xxxxy_1[i] * fe_0 + g_xxxzz_xxxxyy_1[i] * pa_y[i];

        g_xxxyzz_xxxxyz_0[i] = g_xxxzz_xxxxz_1[i] * fe_0 + g_xxxzz_xxxxyz_1[i] * pa_y[i];

        g_xxxyzz_xxxxzz_0[i] = g_xxxzz_xxxxzz_1[i] * pa_y[i];

        g_xxxyzz_xxxyyy_0[i] = 3.0 * g_xxxzz_xxxyy_1[i] * fe_0 + g_xxxzz_xxxyyy_1[i] * pa_y[i];

        g_xxxyzz_xxxyyz_0[i] = 2.0 * g_xxxzz_xxxyz_1[i] * fe_0 + g_xxxzz_xxxyyz_1[i] * pa_y[i];

        g_xxxyzz_xxxyzz_0[i] = g_xxxzz_xxxzz_1[i] * fe_0 + g_xxxzz_xxxyzz_1[i] * pa_y[i];

        g_xxxyzz_xxxzzz_0[i] = g_xxxzz_xxxzzz_1[i] * pa_y[i];

        g_xxxyzz_xxyyyy_0[i] = 4.0 * g_xxxzz_xxyyy_1[i] * fe_0 + g_xxxzz_xxyyyy_1[i] * pa_y[i];

        g_xxxyzz_xxyyyz_0[i] = 3.0 * g_xxxzz_xxyyz_1[i] * fe_0 + g_xxxzz_xxyyyz_1[i] * pa_y[i];

        g_xxxyzz_xxyyzz_0[i] = 2.0 * g_xxxzz_xxyzz_1[i] * fe_0 + g_xxxzz_xxyyzz_1[i] * pa_y[i];

        g_xxxyzz_xxyzzz_0[i] = g_xxxzz_xxzzz_1[i] * fe_0 + g_xxxzz_xxyzzz_1[i] * pa_y[i];

        g_xxxyzz_xxzzzz_0[i] = g_xxxzz_xxzzzz_1[i] * pa_y[i];

        g_xxxyzz_xyyyyy_0[i] = 5.0 * g_xxxzz_xyyyy_1[i] * fe_0 + g_xxxzz_xyyyyy_1[i] * pa_y[i];

        g_xxxyzz_xyyyyz_0[i] = 4.0 * g_xxxzz_xyyyz_1[i] * fe_0 + g_xxxzz_xyyyyz_1[i] * pa_y[i];

        g_xxxyzz_xyyyzz_0[i] = 3.0 * g_xxxzz_xyyzz_1[i] * fe_0 + g_xxxzz_xyyyzz_1[i] * pa_y[i];

        g_xxxyzz_xyyzzz_0[i] = 2.0 * g_xxxzz_xyzzz_1[i] * fe_0 + g_xxxzz_xyyzzz_1[i] * pa_y[i];

        g_xxxyzz_xyzzzz_0[i] = g_xxxzz_xzzzz_1[i] * fe_0 + g_xxxzz_xyzzzz_1[i] * pa_y[i];

        g_xxxyzz_xzzzzz_0[i] = g_xxxzz_xzzzzz_1[i] * pa_y[i];

        g_xxxyzz_yyyyyy_0[i] = 6.0 * g_xxxzz_yyyyy_1[i] * fe_0 + g_xxxzz_yyyyyy_1[i] * pa_y[i];

        g_xxxyzz_yyyyyz_0[i] = 5.0 * g_xxxzz_yyyyz_1[i] * fe_0 + g_xxxzz_yyyyyz_1[i] * pa_y[i];

        g_xxxyzz_yyyyzz_0[i] = 4.0 * g_xxxzz_yyyzz_1[i] * fe_0 + g_xxxzz_yyyyzz_1[i] * pa_y[i];

        g_xxxyzz_yyyzzz_0[i] = 3.0 * g_xxxzz_yyzzz_1[i] * fe_0 + g_xxxzz_yyyzzz_1[i] * pa_y[i];

        g_xxxyzz_yyzzzz_0[i] = 2.0 * g_xxxzz_yzzzz_1[i] * fe_0 + g_xxxzz_yyzzzz_1[i] * pa_y[i];

        g_xxxyzz_yzzzzz_0[i] = g_xxxzz_zzzzz_1[i] * fe_0 + g_xxxzz_yzzzzz_1[i] * pa_y[i];

        g_xxxyzz_zzzzzz_0[i] = g_xxxzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 252-280 components of targeted buffer : II

    auto g_xxxzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 252);

    auto g_xxxzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 253);

    auto g_xxxzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 254);

    auto g_xxxzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 255);

    auto g_xxxzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 256);

    auto g_xxxzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 257);

    auto g_xxxzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 258);

    auto g_xxxzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 259);

    auto g_xxxzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 260);

    auto g_xxxzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 261);

    auto g_xxxzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 262);

    auto g_xxxzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 263);

    auto g_xxxzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 264);

    auto g_xxxzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 265);

    auto g_xxxzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 266);

    auto g_xxxzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 267);

    auto g_xxxzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 268);

    auto g_xxxzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 269);

    auto g_xxxzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 270);

    auto g_xxxzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 271);

    auto g_xxxzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 272);

    auto g_xxxzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 273);

    auto g_xxxzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 274);

    auto g_xxxzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 275);

    auto g_xxxzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 276);

    auto g_xxxzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 277);

    auto g_xxxzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 278);

    auto g_xxxzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 279);

    #pragma omp simd aligned(g_xxxz_xxxxxx_0, g_xxxz_xxxxxx_1, g_xxxz_xxxxxy_0, g_xxxz_xxxxxy_1, g_xxxz_xxxxyy_0, g_xxxz_xxxxyy_1, g_xxxz_xxxyyy_0, g_xxxz_xxxyyy_1, g_xxxz_xxyyyy_0, g_xxxz_xxyyyy_1, g_xxxz_xyyyyy_0, g_xxxz_xyyyyy_1, g_xxxzz_xxxxxx_1, g_xxxzz_xxxxxy_1, g_xxxzz_xxxxyy_1, g_xxxzz_xxxyyy_1, g_xxxzz_xxyyyy_1, g_xxxzz_xyyyyy_1, g_xxxzzz_xxxxxx_0, g_xxxzzz_xxxxxy_0, g_xxxzzz_xxxxxz_0, g_xxxzzz_xxxxyy_0, g_xxxzzz_xxxxyz_0, g_xxxzzz_xxxxzz_0, g_xxxzzz_xxxyyy_0, g_xxxzzz_xxxyyz_0, g_xxxzzz_xxxyzz_0, g_xxxzzz_xxxzzz_0, g_xxxzzz_xxyyyy_0, g_xxxzzz_xxyyyz_0, g_xxxzzz_xxyyzz_0, g_xxxzzz_xxyzzz_0, g_xxxzzz_xxzzzz_0, g_xxxzzz_xyyyyy_0, g_xxxzzz_xyyyyz_0, g_xxxzzz_xyyyzz_0, g_xxxzzz_xyyzzz_0, g_xxxzzz_xyzzzz_0, g_xxxzzz_xzzzzz_0, g_xxxzzz_yyyyyy_0, g_xxxzzz_yyyyyz_0, g_xxxzzz_yyyyzz_0, g_xxxzzz_yyyzzz_0, g_xxxzzz_yyzzzz_0, g_xxxzzz_yzzzzz_0, g_xxxzzz_zzzzzz_0, g_xxzzz_xxxxxz_1, g_xxzzz_xxxxyz_1, g_xxzzz_xxxxz_1, g_xxzzz_xxxxzz_1, g_xxzzz_xxxyyz_1, g_xxzzz_xxxyz_1, g_xxzzz_xxxyzz_1, g_xxzzz_xxxzz_1, g_xxzzz_xxxzzz_1, g_xxzzz_xxyyyz_1, g_xxzzz_xxyyz_1, g_xxzzz_xxyyzz_1, g_xxzzz_xxyzz_1, g_xxzzz_xxyzzz_1, g_xxzzz_xxzzz_1, g_xxzzz_xxzzzz_1, g_xxzzz_xyyyyz_1, g_xxzzz_xyyyz_1, g_xxzzz_xyyyzz_1, g_xxzzz_xyyzz_1, g_xxzzz_xyyzzz_1, g_xxzzz_xyzzz_1, g_xxzzz_xyzzzz_1, g_xxzzz_xzzzz_1, g_xxzzz_xzzzzz_1, g_xxzzz_yyyyyy_1, g_xxzzz_yyyyyz_1, g_xxzzz_yyyyz_1, g_xxzzz_yyyyzz_1, g_xxzzz_yyyzz_1, g_xxzzz_yyyzzz_1, g_xxzzz_yyzzz_1, g_xxzzz_yyzzzz_1, g_xxzzz_yzzzz_1, g_xxzzz_yzzzzz_1, g_xxzzz_zzzzz_1, g_xxzzz_zzzzzz_1, g_xzzz_xxxxxz_0, g_xzzz_xxxxxz_1, g_xzzz_xxxxyz_0, g_xzzz_xxxxyz_1, g_xzzz_xxxxzz_0, g_xzzz_xxxxzz_1, g_xzzz_xxxyyz_0, g_xzzz_xxxyyz_1, g_xzzz_xxxyzz_0, g_xzzz_xxxyzz_1, g_xzzz_xxxzzz_0, g_xzzz_xxxzzz_1, g_xzzz_xxyyyz_0, g_xzzz_xxyyyz_1, g_xzzz_xxyyzz_0, g_xzzz_xxyyzz_1, g_xzzz_xxyzzz_0, g_xzzz_xxyzzz_1, g_xzzz_xxzzzz_0, g_xzzz_xxzzzz_1, g_xzzz_xyyyyz_0, g_xzzz_xyyyyz_1, g_xzzz_xyyyzz_0, g_xzzz_xyyyzz_1, g_xzzz_xyyzzz_0, g_xzzz_xyyzzz_1, g_xzzz_xyzzzz_0, g_xzzz_xyzzzz_1, g_xzzz_xzzzzz_0, g_xzzz_xzzzzz_1, g_xzzz_yyyyyy_0, g_xzzz_yyyyyy_1, g_xzzz_yyyyyz_0, g_xzzz_yyyyyz_1, g_xzzz_yyyyzz_0, g_xzzz_yyyyzz_1, g_xzzz_yyyzzz_0, g_xzzz_yyyzzz_1, g_xzzz_yyzzzz_0, g_xzzz_yyzzzz_1, g_xzzz_yzzzzz_0, g_xzzz_yzzzzz_1, g_xzzz_zzzzzz_0, g_xzzz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzzz_xxxxxx_0[i] = 2.0 * g_xxxz_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxxz_xxxxxx_1[i] * fz_be_0 + g_xxxzz_xxxxxx_1[i] * pa_z[i];

        g_xxxzzz_xxxxxy_0[i] = 2.0 * g_xxxz_xxxxxy_0[i] * fbe_0 - 2.0 * g_xxxz_xxxxxy_1[i] * fz_be_0 + g_xxxzz_xxxxxy_1[i] * pa_z[i];

        g_xxxzzz_xxxxxz_0[i] = 2.0 * g_xzzz_xxxxxz_0[i] * fbe_0 - 2.0 * g_xzzz_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxzzz_xxxxz_1[i] * fe_0 + g_xxzzz_xxxxxz_1[i] * pa_x[i];

        g_xxxzzz_xxxxyy_0[i] = 2.0 * g_xxxz_xxxxyy_0[i] * fbe_0 - 2.0 * g_xxxz_xxxxyy_1[i] * fz_be_0 + g_xxxzz_xxxxyy_1[i] * pa_z[i];

        g_xxxzzz_xxxxyz_0[i] = 2.0 * g_xzzz_xxxxyz_0[i] * fbe_0 - 2.0 * g_xzzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxzzz_xxxyz_1[i] * fe_0 + g_xxzzz_xxxxyz_1[i] * pa_x[i];

        g_xxxzzz_xxxxzz_0[i] = 2.0 * g_xzzz_xxxxzz_0[i] * fbe_0 - 2.0 * g_xzzz_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxzzz_xxxzz_1[i] * fe_0 + g_xxzzz_xxxxzz_1[i] * pa_x[i];

        g_xxxzzz_xxxyyy_0[i] = 2.0 * g_xxxz_xxxyyy_0[i] * fbe_0 - 2.0 * g_xxxz_xxxyyy_1[i] * fz_be_0 + g_xxxzz_xxxyyy_1[i] * pa_z[i];

        g_xxxzzz_xxxyyz_0[i] = 2.0 * g_xzzz_xxxyyz_0[i] * fbe_0 - 2.0 * g_xzzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxzzz_xxyyz_1[i] * fe_0 + g_xxzzz_xxxyyz_1[i] * pa_x[i];

        g_xxxzzz_xxxyzz_0[i] = 2.0 * g_xzzz_xxxyzz_0[i] * fbe_0 - 2.0 * g_xzzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_xxyzz_1[i] * fe_0 + g_xxzzz_xxxyzz_1[i] * pa_x[i];

        g_xxxzzz_xxxzzz_0[i] = 2.0 * g_xzzz_xxxzzz_0[i] * fbe_0 - 2.0 * g_xzzz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_xxzzz_1[i] * fe_0 + g_xxzzz_xxxzzz_1[i] * pa_x[i];

        g_xxxzzz_xxyyyy_0[i] = 2.0 * g_xxxz_xxyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_xxyyyy_1[i] * fz_be_0 + g_xxxzz_xxyyyy_1[i] * pa_z[i];

        g_xxxzzz_xxyyyz_0[i] = 2.0 * g_xzzz_xxyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxzzz_xyyyz_1[i] * fe_0 + g_xxzzz_xxyyyz_1[i] * pa_x[i];

        g_xxxzzz_xxyyzz_0[i] = 2.0 * g_xzzz_xxyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_xyyzz_1[i] * fe_0 + g_xxzzz_xxyyzz_1[i] * pa_x[i];

        g_xxxzzz_xxyzzz_0[i] = 2.0 * g_xzzz_xxyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_xyzzz_1[i] * fe_0 + g_xxzzz_xxyzzz_1[i] * pa_x[i];

        g_xxxzzz_xxzzzz_0[i] = 2.0 * g_xzzz_xxzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_xzzzz_1[i] * fe_0 + g_xxzzz_xxzzzz_1[i] * pa_x[i];

        g_xxxzzz_xyyyyy_0[i] = 2.0 * g_xxxz_xyyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_xyyyyy_1[i] * fz_be_0 + g_xxxzz_xyyyyy_1[i] * pa_z[i];

        g_xxxzzz_xyyyyz_0[i] = 2.0 * g_xzzz_xyyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_xyyyyz_1[i] * fz_be_0 + g_xxzzz_yyyyz_1[i] * fe_0 + g_xxzzz_xyyyyz_1[i] * pa_x[i];

        g_xxxzzz_xyyyzz_0[i] = 2.0 * g_xzzz_xyyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_xyyyzz_1[i] * fz_be_0 + g_xxzzz_yyyzz_1[i] * fe_0 + g_xxzzz_xyyyzz_1[i] * pa_x[i];

        g_xxxzzz_xyyzzz_0[i] = 2.0 * g_xzzz_xyyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_xyyzzz_1[i] * fz_be_0 + g_xxzzz_yyzzz_1[i] * fe_0 + g_xxzzz_xyyzzz_1[i] * pa_x[i];

        g_xxxzzz_xyzzzz_0[i] = 2.0 * g_xzzz_xyzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_xyzzzz_1[i] * fz_be_0 + g_xxzzz_yzzzz_1[i] * fe_0 + g_xxzzz_xyzzzz_1[i] * pa_x[i];

        g_xxxzzz_xzzzzz_0[i] = 2.0 * g_xzzz_xzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_xzzzzz_1[i] * fz_be_0 + g_xxzzz_zzzzz_1[i] * fe_0 + g_xxzzz_xzzzzz_1[i] * pa_x[i];

        g_xxxzzz_yyyyyy_0[i] = 2.0 * g_xzzz_yyyyyy_0[i] * fbe_0 - 2.0 * g_xzzz_yyyyyy_1[i] * fz_be_0 + g_xxzzz_yyyyyy_1[i] * pa_x[i];

        g_xxxzzz_yyyyyz_0[i] = 2.0 * g_xzzz_yyyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_yyyyyz_1[i] * fz_be_0 + g_xxzzz_yyyyyz_1[i] * pa_x[i];

        g_xxxzzz_yyyyzz_0[i] = 2.0 * g_xzzz_yyyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_yyyyzz_1[i] * fz_be_0 + g_xxzzz_yyyyzz_1[i] * pa_x[i];

        g_xxxzzz_yyyzzz_0[i] = 2.0 * g_xzzz_yyyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_yyyzzz_1[i] * fz_be_0 + g_xxzzz_yyyzzz_1[i] * pa_x[i];

        g_xxxzzz_yyzzzz_0[i] = 2.0 * g_xzzz_yyzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_yyzzzz_1[i] * fz_be_0 + g_xxzzz_yyzzzz_1[i] * pa_x[i];

        g_xxxzzz_yzzzzz_0[i] = 2.0 * g_xzzz_yzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_yzzzzz_1[i] * fz_be_0 + g_xxzzz_yzzzzz_1[i] * pa_x[i];

        g_xxxzzz_zzzzzz_0[i] = 2.0 * g_xzzz_zzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_zzzzzz_1[i] * fz_be_0 + g_xxzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 280-308 components of targeted buffer : II

    auto g_xxyyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 280);

    auto g_xxyyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 281);

    auto g_xxyyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 282);

    auto g_xxyyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 283);

    auto g_xxyyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 284);

    auto g_xxyyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 285);

    auto g_xxyyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 286);

    auto g_xxyyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 287);

    auto g_xxyyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 288);

    auto g_xxyyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 289);

    auto g_xxyyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 290);

    auto g_xxyyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 291);

    auto g_xxyyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 292);

    auto g_xxyyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 293);

    auto g_xxyyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 294);

    auto g_xxyyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 295);

    auto g_xxyyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 296);

    auto g_xxyyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 297);

    auto g_xxyyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 298);

    auto g_xxyyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 299);

    auto g_xxyyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 300);

    auto g_xxyyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 301);

    auto g_xxyyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 302);

    auto g_xxyyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 303);

    auto g_xxyyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 304);

    auto g_xxyyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 305);

    auto g_xxyyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 306);

    auto g_xxyyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 307);

    #pragma omp simd aligned(g_xxyy_xxxxxx_0, g_xxyy_xxxxxx_1, g_xxyy_xxxxxz_0, g_xxyy_xxxxxz_1, g_xxyy_xxxxzz_0, g_xxyy_xxxxzz_1, g_xxyy_xxxzzz_0, g_xxyy_xxxzzz_1, g_xxyy_xxzzzz_0, g_xxyy_xxzzzz_1, g_xxyy_xzzzzz_0, g_xxyy_xzzzzz_1, g_xxyyy_xxxxxx_1, g_xxyyy_xxxxxz_1, g_xxyyy_xxxxzz_1, g_xxyyy_xxxzzz_1, g_xxyyy_xxzzzz_1, g_xxyyy_xzzzzz_1, g_xxyyyy_xxxxxx_0, g_xxyyyy_xxxxxy_0, g_xxyyyy_xxxxxz_0, g_xxyyyy_xxxxyy_0, g_xxyyyy_xxxxyz_0, g_xxyyyy_xxxxzz_0, g_xxyyyy_xxxyyy_0, g_xxyyyy_xxxyyz_0, g_xxyyyy_xxxyzz_0, g_xxyyyy_xxxzzz_0, g_xxyyyy_xxyyyy_0, g_xxyyyy_xxyyyz_0, g_xxyyyy_xxyyzz_0, g_xxyyyy_xxyzzz_0, g_xxyyyy_xxzzzz_0, g_xxyyyy_xyyyyy_0, g_xxyyyy_xyyyyz_0, g_xxyyyy_xyyyzz_0, g_xxyyyy_xyyzzz_0, g_xxyyyy_xyzzzz_0, g_xxyyyy_xzzzzz_0, g_xxyyyy_yyyyyy_0, g_xxyyyy_yyyyyz_0, g_xxyyyy_yyyyzz_0, g_xxyyyy_yyyzzz_0, g_xxyyyy_yyzzzz_0, g_xxyyyy_yzzzzz_0, g_xxyyyy_zzzzzz_0, g_xyyyy_xxxxxy_1, g_xyyyy_xxxxy_1, g_xyyyy_xxxxyy_1, g_xyyyy_xxxxyz_1, g_xyyyy_xxxyy_1, g_xyyyy_xxxyyy_1, g_xyyyy_xxxyyz_1, g_xyyyy_xxxyz_1, g_xyyyy_xxxyzz_1, g_xyyyy_xxyyy_1, g_xyyyy_xxyyyy_1, g_xyyyy_xxyyyz_1, g_xyyyy_xxyyz_1, g_xyyyy_xxyyzz_1, g_xyyyy_xxyzz_1, g_xyyyy_xxyzzz_1, g_xyyyy_xyyyy_1, g_xyyyy_xyyyyy_1, g_xyyyy_xyyyyz_1, g_xyyyy_xyyyz_1, g_xyyyy_xyyyzz_1, g_xyyyy_xyyzz_1, g_xyyyy_xyyzzz_1, g_xyyyy_xyzzz_1, g_xyyyy_xyzzzz_1, g_xyyyy_yyyyy_1, g_xyyyy_yyyyyy_1, g_xyyyy_yyyyyz_1, g_xyyyy_yyyyz_1, g_xyyyy_yyyyzz_1, g_xyyyy_yyyzz_1, g_xyyyy_yyyzzz_1, g_xyyyy_yyzzz_1, g_xyyyy_yyzzzz_1, g_xyyyy_yzzzz_1, g_xyyyy_yzzzzz_1, g_xyyyy_zzzzzz_1, g_yyyy_xxxxxy_0, g_yyyy_xxxxxy_1, g_yyyy_xxxxyy_0, g_yyyy_xxxxyy_1, g_yyyy_xxxxyz_0, g_yyyy_xxxxyz_1, g_yyyy_xxxyyy_0, g_yyyy_xxxyyy_1, g_yyyy_xxxyyz_0, g_yyyy_xxxyyz_1, g_yyyy_xxxyzz_0, g_yyyy_xxxyzz_1, g_yyyy_xxyyyy_0, g_yyyy_xxyyyy_1, g_yyyy_xxyyyz_0, g_yyyy_xxyyyz_1, g_yyyy_xxyyzz_0, g_yyyy_xxyyzz_1, g_yyyy_xxyzzz_0, g_yyyy_xxyzzz_1, g_yyyy_xyyyyy_0, g_yyyy_xyyyyy_1, g_yyyy_xyyyyz_0, g_yyyy_xyyyyz_1, g_yyyy_xyyyzz_0, g_yyyy_xyyyzz_1, g_yyyy_xyyzzz_0, g_yyyy_xyyzzz_1, g_yyyy_xyzzzz_0, g_yyyy_xyzzzz_1, g_yyyy_yyyyyy_0, g_yyyy_yyyyyy_1, g_yyyy_yyyyyz_0, g_yyyy_yyyyyz_1, g_yyyy_yyyyzz_0, g_yyyy_yyyyzz_1, g_yyyy_yyyzzz_0, g_yyyy_yyyzzz_1, g_yyyy_yyzzzz_0, g_yyyy_yyzzzz_1, g_yyyy_yzzzzz_0, g_yyyy_yzzzzz_1, g_yyyy_zzzzzz_0, g_yyyy_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyy_xxxxxx_0[i] = 3.0 * g_xxyy_xxxxxx_0[i] * fbe_0 - 3.0 * g_xxyy_xxxxxx_1[i] * fz_be_0 + g_xxyyy_xxxxxx_1[i] * pa_y[i];

        g_xxyyyy_xxxxxy_0[i] = g_yyyy_xxxxxy_0[i] * fbe_0 - g_yyyy_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xyyyy_xxxxy_1[i] * fe_0 + g_xyyyy_xxxxxy_1[i] * pa_x[i];

        g_xxyyyy_xxxxxz_0[i] = 3.0 * g_xxyy_xxxxxz_0[i] * fbe_0 - 3.0 * g_xxyy_xxxxxz_1[i] * fz_be_0 + g_xxyyy_xxxxxz_1[i] * pa_y[i];

        g_xxyyyy_xxxxyy_0[i] = g_yyyy_xxxxyy_0[i] * fbe_0 - g_yyyy_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xyyyy_xxxyy_1[i] * fe_0 + g_xyyyy_xxxxyy_1[i] * pa_x[i];

        g_xxyyyy_xxxxyz_0[i] = g_yyyy_xxxxyz_0[i] * fbe_0 - g_yyyy_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyyyy_xxxyz_1[i] * fe_0 + g_xyyyy_xxxxyz_1[i] * pa_x[i];

        g_xxyyyy_xxxxzz_0[i] = 3.0 * g_xxyy_xxxxzz_0[i] * fbe_0 - 3.0 * g_xxyy_xxxxzz_1[i] * fz_be_0 + g_xxyyy_xxxxzz_1[i] * pa_y[i];

        g_xxyyyy_xxxyyy_0[i] = g_yyyy_xxxyyy_0[i] * fbe_0 - g_yyyy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xyyyy_xxyyy_1[i] * fe_0 + g_xyyyy_xxxyyy_1[i] * pa_x[i];

        g_xxyyyy_xxxyyz_0[i] = g_yyyy_xxxyyz_0[i] * fbe_0 - g_yyyy_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyyyy_xxyyz_1[i] * fe_0 + g_xyyyy_xxxyyz_1[i] * pa_x[i];

        g_xxyyyy_xxxyzz_0[i] = g_yyyy_xxxyzz_0[i] * fbe_0 - g_yyyy_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyyyy_xxyzz_1[i] * fe_0 + g_xyyyy_xxxyzz_1[i] * pa_x[i];

        g_xxyyyy_xxxzzz_0[i] = 3.0 * g_xxyy_xxxzzz_0[i] * fbe_0 - 3.0 * g_xxyy_xxxzzz_1[i] * fz_be_0 + g_xxyyy_xxxzzz_1[i] * pa_y[i];

        g_xxyyyy_xxyyyy_0[i] = g_yyyy_xxyyyy_0[i] * fbe_0 - g_yyyy_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xyyyy_xyyyy_1[i] * fe_0 + g_xyyyy_xxyyyy_1[i] * pa_x[i];

        g_xxyyyy_xxyyyz_0[i] = g_yyyy_xxyyyz_0[i] * fbe_0 - g_yyyy_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyyyy_xyyyz_1[i] * fe_0 + g_xyyyy_xxyyyz_1[i] * pa_x[i];

        g_xxyyyy_xxyyzz_0[i] = g_yyyy_xxyyzz_0[i] * fbe_0 - g_yyyy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_xyyzz_1[i] * fe_0 + g_xyyyy_xxyyzz_1[i] * pa_x[i];

        g_xxyyyy_xxyzzz_0[i] = g_yyyy_xxyzzz_0[i] * fbe_0 - g_yyyy_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_xyzzz_1[i] * fe_0 + g_xyyyy_xxyzzz_1[i] * pa_x[i];

        g_xxyyyy_xxzzzz_0[i] = 3.0 * g_xxyy_xxzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_xxzzzz_1[i] * fz_be_0 + g_xxyyy_xxzzzz_1[i] * pa_y[i];

        g_xxyyyy_xyyyyy_0[i] = g_yyyy_xyyyyy_0[i] * fbe_0 - g_yyyy_xyyyyy_1[i] * fz_be_0 + g_xyyyy_yyyyy_1[i] * fe_0 + g_xyyyy_xyyyyy_1[i] * pa_x[i];

        g_xxyyyy_xyyyyz_0[i] = g_yyyy_xyyyyz_0[i] * fbe_0 - g_yyyy_xyyyyz_1[i] * fz_be_0 + g_xyyyy_yyyyz_1[i] * fe_0 + g_xyyyy_xyyyyz_1[i] * pa_x[i];

        g_xxyyyy_xyyyzz_0[i] = g_yyyy_xyyyzz_0[i] * fbe_0 - g_yyyy_xyyyzz_1[i] * fz_be_0 + g_xyyyy_yyyzz_1[i] * fe_0 + g_xyyyy_xyyyzz_1[i] * pa_x[i];

        g_xxyyyy_xyyzzz_0[i] = g_yyyy_xyyzzz_0[i] * fbe_0 - g_yyyy_xyyzzz_1[i] * fz_be_0 + g_xyyyy_yyzzz_1[i] * fe_0 + g_xyyyy_xyyzzz_1[i] * pa_x[i];

        g_xxyyyy_xyzzzz_0[i] = g_yyyy_xyzzzz_0[i] * fbe_0 - g_yyyy_xyzzzz_1[i] * fz_be_0 + g_xyyyy_yzzzz_1[i] * fe_0 + g_xyyyy_xyzzzz_1[i] * pa_x[i];

        g_xxyyyy_xzzzzz_0[i] = 3.0 * g_xxyy_xzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_xzzzzz_1[i] * fz_be_0 + g_xxyyy_xzzzzz_1[i] * pa_y[i];

        g_xxyyyy_yyyyyy_0[i] = g_yyyy_yyyyyy_0[i] * fbe_0 - g_yyyy_yyyyyy_1[i] * fz_be_0 + g_xyyyy_yyyyyy_1[i] * pa_x[i];

        g_xxyyyy_yyyyyz_0[i] = g_yyyy_yyyyyz_0[i] * fbe_0 - g_yyyy_yyyyyz_1[i] * fz_be_0 + g_xyyyy_yyyyyz_1[i] * pa_x[i];

        g_xxyyyy_yyyyzz_0[i] = g_yyyy_yyyyzz_0[i] * fbe_0 - g_yyyy_yyyyzz_1[i] * fz_be_0 + g_xyyyy_yyyyzz_1[i] * pa_x[i];

        g_xxyyyy_yyyzzz_0[i] = g_yyyy_yyyzzz_0[i] * fbe_0 - g_yyyy_yyyzzz_1[i] * fz_be_0 + g_xyyyy_yyyzzz_1[i] * pa_x[i];

        g_xxyyyy_yyzzzz_0[i] = g_yyyy_yyzzzz_0[i] * fbe_0 - g_yyyy_yyzzzz_1[i] * fz_be_0 + g_xyyyy_yyzzzz_1[i] * pa_x[i];

        g_xxyyyy_yzzzzz_0[i] = g_yyyy_yzzzzz_0[i] * fbe_0 - g_yyyy_yzzzzz_1[i] * fz_be_0 + g_xyyyy_yzzzzz_1[i] * pa_x[i];

        g_xxyyyy_zzzzzz_0[i] = g_yyyy_zzzzzz_0[i] * fbe_0 - g_yyyy_zzzzzz_1[i] * fz_be_0 + g_xyyyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 308-336 components of targeted buffer : II

    auto g_xxyyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 308);

    auto g_xxyyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 309);

    auto g_xxyyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 310);

    auto g_xxyyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 311);

    auto g_xxyyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 312);

    auto g_xxyyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 313);

    auto g_xxyyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 314);

    auto g_xxyyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 315);

    auto g_xxyyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 316);

    auto g_xxyyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 317);

    auto g_xxyyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 318);

    auto g_xxyyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 319);

    auto g_xxyyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 320);

    auto g_xxyyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 321);

    auto g_xxyyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 322);

    auto g_xxyyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 323);

    auto g_xxyyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 324);

    auto g_xxyyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 325);

    auto g_xxyyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 326);

    auto g_xxyyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 327);

    auto g_xxyyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 328);

    auto g_xxyyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 329);

    auto g_xxyyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 330);

    auto g_xxyyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 331);

    auto g_xxyyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 332);

    auto g_xxyyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 333);

    auto g_xxyyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 334);

    auto g_xxyyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 335);

    #pragma omp simd aligned(g_xxyyy_xxxxx_1, g_xxyyy_xxxxxx_1, g_xxyyy_xxxxxy_1, g_xxyyy_xxxxxz_1, g_xxyyy_xxxxy_1, g_xxyyy_xxxxyy_1, g_xxyyy_xxxxyz_1, g_xxyyy_xxxxz_1, g_xxyyy_xxxxzz_1, g_xxyyy_xxxyy_1, g_xxyyy_xxxyyy_1, g_xxyyy_xxxyyz_1, g_xxyyy_xxxyz_1, g_xxyyy_xxxyzz_1, g_xxyyy_xxxzz_1, g_xxyyy_xxxzzz_1, g_xxyyy_xxyyy_1, g_xxyyy_xxyyyy_1, g_xxyyy_xxyyyz_1, g_xxyyy_xxyyz_1, g_xxyyy_xxyyzz_1, g_xxyyy_xxyzz_1, g_xxyyy_xxyzzz_1, g_xxyyy_xxzzz_1, g_xxyyy_xxzzzz_1, g_xxyyy_xyyyy_1, g_xxyyy_xyyyyy_1, g_xxyyy_xyyyyz_1, g_xxyyy_xyyyz_1, g_xxyyy_xyyyzz_1, g_xxyyy_xyyzz_1, g_xxyyy_xyyzzz_1, g_xxyyy_xyzzz_1, g_xxyyy_xyzzzz_1, g_xxyyy_xzzzz_1, g_xxyyy_xzzzzz_1, g_xxyyy_yyyyy_1, g_xxyyy_yyyyyy_1, g_xxyyy_yyyyyz_1, g_xxyyy_yyyyz_1, g_xxyyy_yyyyzz_1, g_xxyyy_yyyzz_1, g_xxyyy_yyyzzz_1, g_xxyyy_yyzzz_1, g_xxyyy_yyzzzz_1, g_xxyyy_yzzzz_1, g_xxyyy_yzzzzz_1, g_xxyyy_zzzzz_1, g_xxyyy_zzzzzz_1, g_xxyyyz_xxxxxx_0, g_xxyyyz_xxxxxy_0, g_xxyyyz_xxxxxz_0, g_xxyyyz_xxxxyy_0, g_xxyyyz_xxxxyz_0, g_xxyyyz_xxxxzz_0, g_xxyyyz_xxxyyy_0, g_xxyyyz_xxxyyz_0, g_xxyyyz_xxxyzz_0, g_xxyyyz_xxxzzz_0, g_xxyyyz_xxyyyy_0, g_xxyyyz_xxyyyz_0, g_xxyyyz_xxyyzz_0, g_xxyyyz_xxyzzz_0, g_xxyyyz_xxzzzz_0, g_xxyyyz_xyyyyy_0, g_xxyyyz_xyyyyz_0, g_xxyyyz_xyyyzz_0, g_xxyyyz_xyyzzz_0, g_xxyyyz_xyzzzz_0, g_xxyyyz_xzzzzz_0, g_xxyyyz_yyyyyy_0, g_xxyyyz_yyyyyz_0, g_xxyyyz_yyyyzz_0, g_xxyyyz_yyyzzz_0, g_xxyyyz_yyzzzz_0, g_xxyyyz_yzzzzz_0, g_xxyyyz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyyz_xxxxxx_0[i] = g_xxyyy_xxxxxx_1[i] * pa_z[i];

        g_xxyyyz_xxxxxy_0[i] = g_xxyyy_xxxxxy_1[i] * pa_z[i];

        g_xxyyyz_xxxxxz_0[i] = g_xxyyy_xxxxx_1[i] * fe_0 + g_xxyyy_xxxxxz_1[i] * pa_z[i];

        g_xxyyyz_xxxxyy_0[i] = g_xxyyy_xxxxyy_1[i] * pa_z[i];

        g_xxyyyz_xxxxyz_0[i] = g_xxyyy_xxxxy_1[i] * fe_0 + g_xxyyy_xxxxyz_1[i] * pa_z[i];

        g_xxyyyz_xxxxzz_0[i] = 2.0 * g_xxyyy_xxxxz_1[i] * fe_0 + g_xxyyy_xxxxzz_1[i] * pa_z[i];

        g_xxyyyz_xxxyyy_0[i] = g_xxyyy_xxxyyy_1[i] * pa_z[i];

        g_xxyyyz_xxxyyz_0[i] = g_xxyyy_xxxyy_1[i] * fe_0 + g_xxyyy_xxxyyz_1[i] * pa_z[i];

        g_xxyyyz_xxxyzz_0[i] = 2.0 * g_xxyyy_xxxyz_1[i] * fe_0 + g_xxyyy_xxxyzz_1[i] * pa_z[i];

        g_xxyyyz_xxxzzz_0[i] = 3.0 * g_xxyyy_xxxzz_1[i] * fe_0 + g_xxyyy_xxxzzz_1[i] * pa_z[i];

        g_xxyyyz_xxyyyy_0[i] = g_xxyyy_xxyyyy_1[i] * pa_z[i];

        g_xxyyyz_xxyyyz_0[i] = g_xxyyy_xxyyy_1[i] * fe_0 + g_xxyyy_xxyyyz_1[i] * pa_z[i];

        g_xxyyyz_xxyyzz_0[i] = 2.0 * g_xxyyy_xxyyz_1[i] * fe_0 + g_xxyyy_xxyyzz_1[i] * pa_z[i];

        g_xxyyyz_xxyzzz_0[i] = 3.0 * g_xxyyy_xxyzz_1[i] * fe_0 + g_xxyyy_xxyzzz_1[i] * pa_z[i];

        g_xxyyyz_xxzzzz_0[i] = 4.0 * g_xxyyy_xxzzz_1[i] * fe_0 + g_xxyyy_xxzzzz_1[i] * pa_z[i];

        g_xxyyyz_xyyyyy_0[i] = g_xxyyy_xyyyyy_1[i] * pa_z[i];

        g_xxyyyz_xyyyyz_0[i] = g_xxyyy_xyyyy_1[i] * fe_0 + g_xxyyy_xyyyyz_1[i] * pa_z[i];

        g_xxyyyz_xyyyzz_0[i] = 2.0 * g_xxyyy_xyyyz_1[i] * fe_0 + g_xxyyy_xyyyzz_1[i] * pa_z[i];

        g_xxyyyz_xyyzzz_0[i] = 3.0 * g_xxyyy_xyyzz_1[i] * fe_0 + g_xxyyy_xyyzzz_1[i] * pa_z[i];

        g_xxyyyz_xyzzzz_0[i] = 4.0 * g_xxyyy_xyzzz_1[i] * fe_0 + g_xxyyy_xyzzzz_1[i] * pa_z[i];

        g_xxyyyz_xzzzzz_0[i] = 5.0 * g_xxyyy_xzzzz_1[i] * fe_0 + g_xxyyy_xzzzzz_1[i] * pa_z[i];

        g_xxyyyz_yyyyyy_0[i] = g_xxyyy_yyyyyy_1[i] * pa_z[i];

        g_xxyyyz_yyyyyz_0[i] = g_xxyyy_yyyyy_1[i] * fe_0 + g_xxyyy_yyyyyz_1[i] * pa_z[i];

        g_xxyyyz_yyyyzz_0[i] = 2.0 * g_xxyyy_yyyyz_1[i] * fe_0 + g_xxyyy_yyyyzz_1[i] * pa_z[i];

        g_xxyyyz_yyyzzz_0[i] = 3.0 * g_xxyyy_yyyzz_1[i] * fe_0 + g_xxyyy_yyyzzz_1[i] * pa_z[i];

        g_xxyyyz_yyzzzz_0[i] = 4.0 * g_xxyyy_yyzzz_1[i] * fe_0 + g_xxyyy_yyzzzz_1[i] * pa_z[i];

        g_xxyyyz_yzzzzz_0[i] = 5.0 * g_xxyyy_yzzzz_1[i] * fe_0 + g_xxyyy_yzzzzz_1[i] * pa_z[i];

        g_xxyyyz_zzzzzz_0[i] = 6.0 * g_xxyyy_zzzzz_1[i] * fe_0 + g_xxyyy_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 336-364 components of targeted buffer : II

    auto g_xxyyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 336);

    auto g_xxyyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 337);

    auto g_xxyyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 338);

    auto g_xxyyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 339);

    auto g_xxyyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 340);

    auto g_xxyyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 341);

    auto g_xxyyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 342);

    auto g_xxyyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 343);

    auto g_xxyyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 344);

    auto g_xxyyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 345);

    auto g_xxyyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 346);

    auto g_xxyyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 347);

    auto g_xxyyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 348);

    auto g_xxyyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 349);

    auto g_xxyyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 350);

    auto g_xxyyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 351);

    auto g_xxyyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 352);

    auto g_xxyyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 353);

    auto g_xxyyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 354);

    auto g_xxyyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 355);

    auto g_xxyyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 356);

    auto g_xxyyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 357);

    auto g_xxyyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 358);

    auto g_xxyyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 359);

    auto g_xxyyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 360);

    auto g_xxyyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 361);

    auto g_xxyyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 362);

    auto g_xxyyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 363);

    #pragma omp simd aligned(g_xxyy_xxxxxy_0, g_xxyy_xxxxxy_1, g_xxyy_xxxxyy_0, g_xxyy_xxxxyy_1, g_xxyy_xxxyyy_0, g_xxyy_xxxyyy_1, g_xxyy_xxyyyy_0, g_xxyy_xxyyyy_1, g_xxyy_xyyyyy_0, g_xxyy_xyyyyy_1, g_xxyyz_xxxxxy_1, g_xxyyz_xxxxyy_1, g_xxyyz_xxxyyy_1, g_xxyyz_xxyyyy_1, g_xxyyz_xyyyyy_1, g_xxyyzz_xxxxxx_0, g_xxyyzz_xxxxxy_0, g_xxyyzz_xxxxxz_0, g_xxyyzz_xxxxyy_0, g_xxyyzz_xxxxyz_0, g_xxyyzz_xxxxzz_0, g_xxyyzz_xxxyyy_0, g_xxyyzz_xxxyyz_0, g_xxyyzz_xxxyzz_0, g_xxyyzz_xxxzzz_0, g_xxyyzz_xxyyyy_0, g_xxyyzz_xxyyyz_0, g_xxyyzz_xxyyzz_0, g_xxyyzz_xxyzzz_0, g_xxyyzz_xxzzzz_0, g_xxyyzz_xyyyyy_0, g_xxyyzz_xyyyyz_0, g_xxyyzz_xyyyzz_0, g_xxyyzz_xyyzzz_0, g_xxyyzz_xyzzzz_0, g_xxyyzz_xzzzzz_0, g_xxyyzz_yyyyyy_0, g_xxyyzz_yyyyyz_0, g_xxyyzz_yyyyzz_0, g_xxyyzz_yyyzzz_0, g_xxyyzz_yyzzzz_0, g_xxyyzz_yzzzzz_0, g_xxyyzz_zzzzzz_0, g_xxyzz_xxxxxx_1, g_xxyzz_xxxxxz_1, g_xxyzz_xxxxzz_1, g_xxyzz_xxxzzz_1, g_xxyzz_xxzzzz_1, g_xxyzz_xzzzzz_1, g_xxzz_xxxxxx_0, g_xxzz_xxxxxx_1, g_xxzz_xxxxxz_0, g_xxzz_xxxxxz_1, g_xxzz_xxxxzz_0, g_xxzz_xxxxzz_1, g_xxzz_xxxzzz_0, g_xxzz_xxxzzz_1, g_xxzz_xxzzzz_0, g_xxzz_xxzzzz_1, g_xxzz_xzzzzz_0, g_xxzz_xzzzzz_1, g_xyyzz_xxxxyz_1, g_xyyzz_xxxyyz_1, g_xyyzz_xxxyz_1, g_xyyzz_xxxyzz_1, g_xyyzz_xxyyyz_1, g_xyyzz_xxyyz_1, g_xyyzz_xxyyzz_1, g_xyyzz_xxyzz_1, g_xyyzz_xxyzzz_1, g_xyyzz_xyyyyz_1, g_xyyzz_xyyyz_1, g_xyyzz_xyyyzz_1, g_xyyzz_xyyzz_1, g_xyyzz_xyyzzz_1, g_xyyzz_xyzzz_1, g_xyyzz_xyzzzz_1, g_xyyzz_yyyyyy_1, g_xyyzz_yyyyyz_1, g_xyyzz_yyyyz_1, g_xyyzz_yyyyzz_1, g_xyyzz_yyyzz_1, g_xyyzz_yyyzzz_1, g_xyyzz_yyzzz_1, g_xyyzz_yyzzzz_1, g_xyyzz_yzzzz_1, g_xyyzz_yzzzzz_1, g_xyyzz_zzzzzz_1, g_yyzz_xxxxyz_0, g_yyzz_xxxxyz_1, g_yyzz_xxxyyz_0, g_yyzz_xxxyyz_1, g_yyzz_xxxyzz_0, g_yyzz_xxxyzz_1, g_yyzz_xxyyyz_0, g_yyzz_xxyyyz_1, g_yyzz_xxyyzz_0, g_yyzz_xxyyzz_1, g_yyzz_xxyzzz_0, g_yyzz_xxyzzz_1, g_yyzz_xyyyyz_0, g_yyzz_xyyyyz_1, g_yyzz_xyyyzz_0, g_yyzz_xyyyzz_1, g_yyzz_xyyzzz_0, g_yyzz_xyyzzz_1, g_yyzz_xyzzzz_0, g_yyzz_xyzzzz_1, g_yyzz_yyyyyy_0, g_yyzz_yyyyyy_1, g_yyzz_yyyyyz_0, g_yyzz_yyyyyz_1, g_yyzz_yyyyzz_0, g_yyzz_yyyyzz_1, g_yyzz_yyyzzz_0, g_yyzz_yyyzzz_1, g_yyzz_yyzzzz_0, g_yyzz_yyzzzz_1, g_yyzz_yzzzzz_0, g_yyzz_yzzzzz_1, g_yyzz_zzzzzz_0, g_yyzz_zzzzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyzz_xxxxxx_0[i] = g_xxzz_xxxxxx_0[i] * fbe_0 - g_xxzz_xxxxxx_1[i] * fz_be_0 + g_xxyzz_xxxxxx_1[i] * pa_y[i];

        g_xxyyzz_xxxxxy_0[i] = g_xxyy_xxxxxy_0[i] * fbe_0 - g_xxyy_xxxxxy_1[i] * fz_be_0 + g_xxyyz_xxxxxy_1[i] * pa_z[i];

        g_xxyyzz_xxxxxz_0[i] = g_xxzz_xxxxxz_0[i] * fbe_0 - g_xxzz_xxxxxz_1[i] * fz_be_0 + g_xxyzz_xxxxxz_1[i] * pa_y[i];

        g_xxyyzz_xxxxyy_0[i] = g_xxyy_xxxxyy_0[i] * fbe_0 - g_xxyy_xxxxyy_1[i] * fz_be_0 + g_xxyyz_xxxxyy_1[i] * pa_z[i];

        g_xxyyzz_xxxxyz_0[i] = g_yyzz_xxxxyz_0[i] * fbe_0 - g_yyzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyyzz_xxxyz_1[i] * fe_0 + g_xyyzz_xxxxyz_1[i] * pa_x[i];

        g_xxyyzz_xxxxzz_0[i] = g_xxzz_xxxxzz_0[i] * fbe_0 - g_xxzz_xxxxzz_1[i] * fz_be_0 + g_xxyzz_xxxxzz_1[i] * pa_y[i];

        g_xxyyzz_xxxyyy_0[i] = g_xxyy_xxxyyy_0[i] * fbe_0 - g_xxyy_xxxyyy_1[i] * fz_be_0 + g_xxyyz_xxxyyy_1[i] * pa_z[i];

        g_xxyyzz_xxxyyz_0[i] = g_yyzz_xxxyyz_0[i] * fbe_0 - g_yyzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyyzz_xxyyz_1[i] * fe_0 + g_xyyzz_xxxyyz_1[i] * pa_x[i];

        g_xxyyzz_xxxyzz_0[i] = g_yyzz_xxxyzz_0[i] * fbe_0 - g_yyzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyyzz_xxyzz_1[i] * fe_0 + g_xyyzz_xxxyzz_1[i] * pa_x[i];

        g_xxyyzz_xxxzzz_0[i] = g_xxzz_xxxzzz_0[i] * fbe_0 - g_xxzz_xxxzzz_1[i] * fz_be_0 + g_xxyzz_xxxzzz_1[i] * pa_y[i];

        g_xxyyzz_xxyyyy_0[i] = g_xxyy_xxyyyy_0[i] * fbe_0 - g_xxyy_xxyyyy_1[i] * fz_be_0 + g_xxyyz_xxyyyy_1[i] * pa_z[i];

        g_xxyyzz_xxyyyz_0[i] = g_yyzz_xxyyyz_0[i] * fbe_0 - g_yyzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyyzz_xyyyz_1[i] * fe_0 + g_xyyzz_xxyyyz_1[i] * pa_x[i];

        g_xxyyzz_xxyyzz_0[i] = g_yyzz_xxyyzz_0[i] * fbe_0 - g_yyzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_xyyzz_1[i] * fe_0 + g_xyyzz_xxyyzz_1[i] * pa_x[i];

        g_xxyyzz_xxyzzz_0[i] = g_yyzz_xxyzzz_0[i] * fbe_0 - g_yyzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_xyzzz_1[i] * fe_0 + g_xyyzz_xxyzzz_1[i] * pa_x[i];

        g_xxyyzz_xxzzzz_0[i] = g_xxzz_xxzzzz_0[i] * fbe_0 - g_xxzz_xxzzzz_1[i] * fz_be_0 + g_xxyzz_xxzzzz_1[i] * pa_y[i];

        g_xxyyzz_xyyyyy_0[i] = g_xxyy_xyyyyy_0[i] * fbe_0 - g_xxyy_xyyyyy_1[i] * fz_be_0 + g_xxyyz_xyyyyy_1[i] * pa_z[i];

        g_xxyyzz_xyyyyz_0[i] = g_yyzz_xyyyyz_0[i] * fbe_0 - g_yyzz_xyyyyz_1[i] * fz_be_0 + g_xyyzz_yyyyz_1[i] * fe_0 + g_xyyzz_xyyyyz_1[i] * pa_x[i];

        g_xxyyzz_xyyyzz_0[i] = g_yyzz_xyyyzz_0[i] * fbe_0 - g_yyzz_xyyyzz_1[i] * fz_be_0 + g_xyyzz_yyyzz_1[i] * fe_0 + g_xyyzz_xyyyzz_1[i] * pa_x[i];

        g_xxyyzz_xyyzzz_0[i] = g_yyzz_xyyzzz_0[i] * fbe_0 - g_yyzz_xyyzzz_1[i] * fz_be_0 + g_xyyzz_yyzzz_1[i] * fe_0 + g_xyyzz_xyyzzz_1[i] * pa_x[i];

        g_xxyyzz_xyzzzz_0[i] = g_yyzz_xyzzzz_0[i] * fbe_0 - g_yyzz_xyzzzz_1[i] * fz_be_0 + g_xyyzz_yzzzz_1[i] * fe_0 + g_xyyzz_xyzzzz_1[i] * pa_x[i];

        g_xxyyzz_xzzzzz_0[i] = g_xxzz_xzzzzz_0[i] * fbe_0 - g_xxzz_xzzzzz_1[i] * fz_be_0 + g_xxyzz_xzzzzz_1[i] * pa_y[i];

        g_xxyyzz_yyyyyy_0[i] = g_yyzz_yyyyyy_0[i] * fbe_0 - g_yyzz_yyyyyy_1[i] * fz_be_0 + g_xyyzz_yyyyyy_1[i] * pa_x[i];

        g_xxyyzz_yyyyyz_0[i] = g_yyzz_yyyyyz_0[i] * fbe_0 - g_yyzz_yyyyyz_1[i] * fz_be_0 + g_xyyzz_yyyyyz_1[i] * pa_x[i];

        g_xxyyzz_yyyyzz_0[i] = g_yyzz_yyyyzz_0[i] * fbe_0 - g_yyzz_yyyyzz_1[i] * fz_be_0 + g_xyyzz_yyyyzz_1[i] * pa_x[i];

        g_xxyyzz_yyyzzz_0[i] = g_yyzz_yyyzzz_0[i] * fbe_0 - g_yyzz_yyyzzz_1[i] * fz_be_0 + g_xyyzz_yyyzzz_1[i] * pa_x[i];

        g_xxyyzz_yyzzzz_0[i] = g_yyzz_yyzzzz_0[i] * fbe_0 - g_yyzz_yyzzzz_1[i] * fz_be_0 + g_xyyzz_yyzzzz_1[i] * pa_x[i];

        g_xxyyzz_yzzzzz_0[i] = g_yyzz_yzzzzz_0[i] * fbe_0 - g_yyzz_yzzzzz_1[i] * fz_be_0 + g_xyyzz_yzzzzz_1[i] * pa_x[i];

        g_xxyyzz_zzzzzz_0[i] = g_yyzz_zzzzzz_0[i] * fbe_0 - g_yyzz_zzzzzz_1[i] * fz_be_0 + g_xyyzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 364-392 components of targeted buffer : II

    auto g_xxyzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 364);

    auto g_xxyzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 365);

    auto g_xxyzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 366);

    auto g_xxyzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 367);

    auto g_xxyzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 368);

    auto g_xxyzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 369);

    auto g_xxyzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 370);

    auto g_xxyzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 371);

    auto g_xxyzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 372);

    auto g_xxyzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 373);

    auto g_xxyzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 374);

    auto g_xxyzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 375);

    auto g_xxyzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 376);

    auto g_xxyzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 377);

    auto g_xxyzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 378);

    auto g_xxyzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 379);

    auto g_xxyzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 380);

    auto g_xxyzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 381);

    auto g_xxyzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 382);

    auto g_xxyzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 383);

    auto g_xxyzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 384);

    auto g_xxyzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 385);

    auto g_xxyzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 386);

    auto g_xxyzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 387);

    auto g_xxyzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 388);

    auto g_xxyzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 389);

    auto g_xxyzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 390);

    auto g_xxyzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 391);

    #pragma omp simd aligned(g_xxyzzz_xxxxxx_0, g_xxyzzz_xxxxxy_0, g_xxyzzz_xxxxxz_0, g_xxyzzz_xxxxyy_0, g_xxyzzz_xxxxyz_0, g_xxyzzz_xxxxzz_0, g_xxyzzz_xxxyyy_0, g_xxyzzz_xxxyyz_0, g_xxyzzz_xxxyzz_0, g_xxyzzz_xxxzzz_0, g_xxyzzz_xxyyyy_0, g_xxyzzz_xxyyyz_0, g_xxyzzz_xxyyzz_0, g_xxyzzz_xxyzzz_0, g_xxyzzz_xxzzzz_0, g_xxyzzz_xyyyyy_0, g_xxyzzz_xyyyyz_0, g_xxyzzz_xyyyzz_0, g_xxyzzz_xyyzzz_0, g_xxyzzz_xyzzzz_0, g_xxyzzz_xzzzzz_0, g_xxyzzz_yyyyyy_0, g_xxyzzz_yyyyyz_0, g_xxyzzz_yyyyzz_0, g_xxyzzz_yyyzzz_0, g_xxyzzz_yyzzzz_0, g_xxyzzz_yzzzzz_0, g_xxyzzz_zzzzzz_0, g_xxzzz_xxxxx_1, g_xxzzz_xxxxxx_1, g_xxzzz_xxxxxy_1, g_xxzzz_xxxxxz_1, g_xxzzz_xxxxy_1, g_xxzzz_xxxxyy_1, g_xxzzz_xxxxyz_1, g_xxzzz_xxxxz_1, g_xxzzz_xxxxzz_1, g_xxzzz_xxxyy_1, g_xxzzz_xxxyyy_1, g_xxzzz_xxxyyz_1, g_xxzzz_xxxyz_1, g_xxzzz_xxxyzz_1, g_xxzzz_xxxzz_1, g_xxzzz_xxxzzz_1, g_xxzzz_xxyyy_1, g_xxzzz_xxyyyy_1, g_xxzzz_xxyyyz_1, g_xxzzz_xxyyz_1, g_xxzzz_xxyyzz_1, g_xxzzz_xxyzz_1, g_xxzzz_xxyzzz_1, g_xxzzz_xxzzz_1, g_xxzzz_xxzzzz_1, g_xxzzz_xyyyy_1, g_xxzzz_xyyyyy_1, g_xxzzz_xyyyyz_1, g_xxzzz_xyyyz_1, g_xxzzz_xyyyzz_1, g_xxzzz_xyyzz_1, g_xxzzz_xyyzzz_1, g_xxzzz_xyzzz_1, g_xxzzz_xyzzzz_1, g_xxzzz_xzzzz_1, g_xxzzz_xzzzzz_1, g_xxzzz_yyyyy_1, g_xxzzz_yyyyyy_1, g_xxzzz_yyyyyz_1, g_xxzzz_yyyyz_1, g_xxzzz_yyyyzz_1, g_xxzzz_yyyzz_1, g_xxzzz_yyyzzz_1, g_xxzzz_yyzzz_1, g_xxzzz_yyzzzz_1, g_xxzzz_yzzzz_1, g_xxzzz_yzzzzz_1, g_xxzzz_zzzzz_1, g_xxzzz_zzzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzzz_xxxxxx_0[i] = g_xxzzz_xxxxxx_1[i] * pa_y[i];

        g_xxyzzz_xxxxxy_0[i] = g_xxzzz_xxxxx_1[i] * fe_0 + g_xxzzz_xxxxxy_1[i] * pa_y[i];

        g_xxyzzz_xxxxxz_0[i] = g_xxzzz_xxxxxz_1[i] * pa_y[i];

        g_xxyzzz_xxxxyy_0[i] = 2.0 * g_xxzzz_xxxxy_1[i] * fe_0 + g_xxzzz_xxxxyy_1[i] * pa_y[i];

        g_xxyzzz_xxxxyz_0[i] = g_xxzzz_xxxxz_1[i] * fe_0 + g_xxzzz_xxxxyz_1[i] * pa_y[i];

        g_xxyzzz_xxxxzz_0[i] = g_xxzzz_xxxxzz_1[i] * pa_y[i];

        g_xxyzzz_xxxyyy_0[i] = 3.0 * g_xxzzz_xxxyy_1[i] * fe_0 + g_xxzzz_xxxyyy_1[i] * pa_y[i];

        g_xxyzzz_xxxyyz_0[i] = 2.0 * g_xxzzz_xxxyz_1[i] * fe_0 + g_xxzzz_xxxyyz_1[i] * pa_y[i];

        g_xxyzzz_xxxyzz_0[i] = g_xxzzz_xxxzz_1[i] * fe_0 + g_xxzzz_xxxyzz_1[i] * pa_y[i];

        g_xxyzzz_xxxzzz_0[i] = g_xxzzz_xxxzzz_1[i] * pa_y[i];

        g_xxyzzz_xxyyyy_0[i] = 4.0 * g_xxzzz_xxyyy_1[i] * fe_0 + g_xxzzz_xxyyyy_1[i] * pa_y[i];

        g_xxyzzz_xxyyyz_0[i] = 3.0 * g_xxzzz_xxyyz_1[i] * fe_0 + g_xxzzz_xxyyyz_1[i] * pa_y[i];

        g_xxyzzz_xxyyzz_0[i] = 2.0 * g_xxzzz_xxyzz_1[i] * fe_0 + g_xxzzz_xxyyzz_1[i] * pa_y[i];

        g_xxyzzz_xxyzzz_0[i] = g_xxzzz_xxzzz_1[i] * fe_0 + g_xxzzz_xxyzzz_1[i] * pa_y[i];

        g_xxyzzz_xxzzzz_0[i] = g_xxzzz_xxzzzz_1[i] * pa_y[i];

        g_xxyzzz_xyyyyy_0[i] = 5.0 * g_xxzzz_xyyyy_1[i] * fe_0 + g_xxzzz_xyyyyy_1[i] * pa_y[i];

        g_xxyzzz_xyyyyz_0[i] = 4.0 * g_xxzzz_xyyyz_1[i] * fe_0 + g_xxzzz_xyyyyz_1[i] * pa_y[i];

        g_xxyzzz_xyyyzz_0[i] = 3.0 * g_xxzzz_xyyzz_1[i] * fe_0 + g_xxzzz_xyyyzz_1[i] * pa_y[i];

        g_xxyzzz_xyyzzz_0[i] = 2.0 * g_xxzzz_xyzzz_1[i] * fe_0 + g_xxzzz_xyyzzz_1[i] * pa_y[i];

        g_xxyzzz_xyzzzz_0[i] = g_xxzzz_xzzzz_1[i] * fe_0 + g_xxzzz_xyzzzz_1[i] * pa_y[i];

        g_xxyzzz_xzzzzz_0[i] = g_xxzzz_xzzzzz_1[i] * pa_y[i];

        g_xxyzzz_yyyyyy_0[i] = 6.0 * g_xxzzz_yyyyy_1[i] * fe_0 + g_xxzzz_yyyyyy_1[i] * pa_y[i];

        g_xxyzzz_yyyyyz_0[i] = 5.0 * g_xxzzz_yyyyz_1[i] * fe_0 + g_xxzzz_yyyyyz_1[i] * pa_y[i];

        g_xxyzzz_yyyyzz_0[i] = 4.0 * g_xxzzz_yyyzz_1[i] * fe_0 + g_xxzzz_yyyyzz_1[i] * pa_y[i];

        g_xxyzzz_yyyzzz_0[i] = 3.0 * g_xxzzz_yyzzz_1[i] * fe_0 + g_xxzzz_yyyzzz_1[i] * pa_y[i];

        g_xxyzzz_yyzzzz_0[i] = 2.0 * g_xxzzz_yzzzz_1[i] * fe_0 + g_xxzzz_yyzzzz_1[i] * pa_y[i];

        g_xxyzzz_yzzzzz_0[i] = g_xxzzz_zzzzz_1[i] * fe_0 + g_xxzzz_yzzzzz_1[i] * pa_y[i];

        g_xxyzzz_zzzzzz_0[i] = g_xxzzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 392-420 components of targeted buffer : II

    auto g_xxzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 392);

    auto g_xxzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 393);

    auto g_xxzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 394);

    auto g_xxzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 395);

    auto g_xxzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 396);

    auto g_xxzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 397);

    auto g_xxzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 398);

    auto g_xxzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 399);

    auto g_xxzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 400);

    auto g_xxzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 401);

    auto g_xxzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 402);

    auto g_xxzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 403);

    auto g_xxzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 404);

    auto g_xxzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 405);

    auto g_xxzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 406);

    auto g_xxzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 407);

    auto g_xxzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 408);

    auto g_xxzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 409);

    auto g_xxzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 410);

    auto g_xxzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 411);

    auto g_xxzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 412);

    auto g_xxzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 413);

    auto g_xxzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 414);

    auto g_xxzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 415);

    auto g_xxzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 416);

    auto g_xxzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 417);

    auto g_xxzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 418);

    auto g_xxzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 419);

    #pragma omp simd aligned(g_xxzz_xxxxxx_0, g_xxzz_xxxxxx_1, g_xxzz_xxxxxy_0, g_xxzz_xxxxxy_1, g_xxzz_xxxxyy_0, g_xxzz_xxxxyy_1, g_xxzz_xxxyyy_0, g_xxzz_xxxyyy_1, g_xxzz_xxyyyy_0, g_xxzz_xxyyyy_1, g_xxzz_xyyyyy_0, g_xxzz_xyyyyy_1, g_xxzzz_xxxxxx_1, g_xxzzz_xxxxxy_1, g_xxzzz_xxxxyy_1, g_xxzzz_xxxyyy_1, g_xxzzz_xxyyyy_1, g_xxzzz_xyyyyy_1, g_xxzzzz_xxxxxx_0, g_xxzzzz_xxxxxy_0, g_xxzzzz_xxxxxz_0, g_xxzzzz_xxxxyy_0, g_xxzzzz_xxxxyz_0, g_xxzzzz_xxxxzz_0, g_xxzzzz_xxxyyy_0, g_xxzzzz_xxxyyz_0, g_xxzzzz_xxxyzz_0, g_xxzzzz_xxxzzz_0, g_xxzzzz_xxyyyy_0, g_xxzzzz_xxyyyz_0, g_xxzzzz_xxyyzz_0, g_xxzzzz_xxyzzz_0, g_xxzzzz_xxzzzz_0, g_xxzzzz_xyyyyy_0, g_xxzzzz_xyyyyz_0, g_xxzzzz_xyyyzz_0, g_xxzzzz_xyyzzz_0, g_xxzzzz_xyzzzz_0, g_xxzzzz_xzzzzz_0, g_xxzzzz_yyyyyy_0, g_xxzzzz_yyyyyz_0, g_xxzzzz_yyyyzz_0, g_xxzzzz_yyyzzz_0, g_xxzzzz_yyzzzz_0, g_xxzzzz_yzzzzz_0, g_xxzzzz_zzzzzz_0, g_xzzzz_xxxxxz_1, g_xzzzz_xxxxyz_1, g_xzzzz_xxxxz_1, g_xzzzz_xxxxzz_1, g_xzzzz_xxxyyz_1, g_xzzzz_xxxyz_1, g_xzzzz_xxxyzz_1, g_xzzzz_xxxzz_1, g_xzzzz_xxxzzz_1, g_xzzzz_xxyyyz_1, g_xzzzz_xxyyz_1, g_xzzzz_xxyyzz_1, g_xzzzz_xxyzz_1, g_xzzzz_xxyzzz_1, g_xzzzz_xxzzz_1, g_xzzzz_xxzzzz_1, g_xzzzz_xyyyyz_1, g_xzzzz_xyyyz_1, g_xzzzz_xyyyzz_1, g_xzzzz_xyyzz_1, g_xzzzz_xyyzzz_1, g_xzzzz_xyzzz_1, g_xzzzz_xyzzzz_1, g_xzzzz_xzzzz_1, g_xzzzz_xzzzzz_1, g_xzzzz_yyyyyy_1, g_xzzzz_yyyyyz_1, g_xzzzz_yyyyz_1, g_xzzzz_yyyyzz_1, g_xzzzz_yyyzz_1, g_xzzzz_yyyzzz_1, g_xzzzz_yyzzz_1, g_xzzzz_yyzzzz_1, g_xzzzz_yzzzz_1, g_xzzzz_yzzzzz_1, g_xzzzz_zzzzz_1, g_xzzzz_zzzzzz_1, g_zzzz_xxxxxz_0, g_zzzz_xxxxxz_1, g_zzzz_xxxxyz_0, g_zzzz_xxxxyz_1, g_zzzz_xxxxzz_0, g_zzzz_xxxxzz_1, g_zzzz_xxxyyz_0, g_zzzz_xxxyyz_1, g_zzzz_xxxyzz_0, g_zzzz_xxxyzz_1, g_zzzz_xxxzzz_0, g_zzzz_xxxzzz_1, g_zzzz_xxyyyz_0, g_zzzz_xxyyyz_1, g_zzzz_xxyyzz_0, g_zzzz_xxyyzz_1, g_zzzz_xxyzzz_0, g_zzzz_xxyzzz_1, g_zzzz_xxzzzz_0, g_zzzz_xxzzzz_1, g_zzzz_xyyyyz_0, g_zzzz_xyyyyz_1, g_zzzz_xyyyzz_0, g_zzzz_xyyyzz_1, g_zzzz_xyyzzz_0, g_zzzz_xyyzzz_1, g_zzzz_xyzzzz_0, g_zzzz_xyzzzz_1, g_zzzz_xzzzzz_0, g_zzzz_xzzzzz_1, g_zzzz_yyyyyy_0, g_zzzz_yyyyyy_1, g_zzzz_yyyyyz_0, g_zzzz_yyyyyz_1, g_zzzz_yyyyzz_0, g_zzzz_yyyyzz_1, g_zzzz_yyyzzz_0, g_zzzz_yyyzzz_1, g_zzzz_yyzzzz_0, g_zzzz_yyzzzz_1, g_zzzz_yzzzzz_0, g_zzzz_yzzzzz_1, g_zzzz_zzzzzz_0, g_zzzz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzzz_xxxxxx_0[i] = 3.0 * g_xxzz_xxxxxx_0[i] * fbe_0 - 3.0 * g_xxzz_xxxxxx_1[i] * fz_be_0 + g_xxzzz_xxxxxx_1[i] * pa_z[i];

        g_xxzzzz_xxxxxy_0[i] = 3.0 * g_xxzz_xxxxxy_0[i] * fbe_0 - 3.0 * g_xxzz_xxxxxy_1[i] * fz_be_0 + g_xxzzz_xxxxxy_1[i] * pa_z[i];

        g_xxzzzz_xxxxxz_0[i] = g_zzzz_xxxxxz_0[i] * fbe_0 - g_zzzz_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xzzzz_xxxxz_1[i] * fe_0 + g_xzzzz_xxxxxz_1[i] * pa_x[i];

        g_xxzzzz_xxxxyy_0[i] = 3.0 * g_xxzz_xxxxyy_0[i] * fbe_0 - 3.0 * g_xxzz_xxxxyy_1[i] * fz_be_0 + g_xxzzz_xxxxyy_1[i] * pa_z[i];

        g_xxzzzz_xxxxyz_0[i] = g_zzzz_xxxxyz_0[i] * fbe_0 - g_zzzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xzzzz_xxxyz_1[i] * fe_0 + g_xzzzz_xxxxyz_1[i] * pa_x[i];

        g_xxzzzz_xxxxzz_0[i] = g_zzzz_xxxxzz_0[i] * fbe_0 - g_zzzz_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xzzzz_xxxzz_1[i] * fe_0 + g_xzzzz_xxxxzz_1[i] * pa_x[i];

        g_xxzzzz_xxxyyy_0[i] = 3.0 * g_xxzz_xxxyyy_0[i] * fbe_0 - 3.0 * g_xxzz_xxxyyy_1[i] * fz_be_0 + g_xxzzz_xxxyyy_1[i] * pa_z[i];

        g_xxzzzz_xxxyyz_0[i] = g_zzzz_xxxyyz_0[i] * fbe_0 - g_zzzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xzzzz_xxyyz_1[i] * fe_0 + g_xzzzz_xxxyyz_1[i] * pa_x[i];

        g_xxzzzz_xxxyzz_0[i] = g_zzzz_xxxyzz_0[i] * fbe_0 - g_zzzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_xxyzz_1[i] * fe_0 + g_xzzzz_xxxyzz_1[i] * pa_x[i];

        g_xxzzzz_xxxzzz_0[i] = g_zzzz_xxxzzz_0[i] * fbe_0 - g_zzzz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_xxzzz_1[i] * fe_0 + g_xzzzz_xxxzzz_1[i] * pa_x[i];

        g_xxzzzz_xxyyyy_0[i] = 3.0 * g_xxzz_xxyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_xxyyyy_1[i] * fz_be_0 + g_xxzzz_xxyyyy_1[i] * pa_z[i];

        g_xxzzzz_xxyyyz_0[i] = g_zzzz_xxyyyz_0[i] * fbe_0 - g_zzzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xzzzz_xyyyz_1[i] * fe_0 + g_xzzzz_xxyyyz_1[i] * pa_x[i];

        g_xxzzzz_xxyyzz_0[i] = g_zzzz_xxyyzz_0[i] * fbe_0 - g_zzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_xyyzz_1[i] * fe_0 + g_xzzzz_xxyyzz_1[i] * pa_x[i];

        g_xxzzzz_xxyzzz_0[i] = g_zzzz_xxyzzz_0[i] * fbe_0 - g_zzzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_xyzzz_1[i] * fe_0 + g_xzzzz_xxyzzz_1[i] * pa_x[i];

        g_xxzzzz_xxzzzz_0[i] = g_zzzz_xxzzzz_0[i] * fbe_0 - g_zzzz_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_xzzzz_1[i] * fe_0 + g_xzzzz_xxzzzz_1[i] * pa_x[i];

        g_xxzzzz_xyyyyy_0[i] = 3.0 * g_xxzz_xyyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_xyyyyy_1[i] * fz_be_0 + g_xxzzz_xyyyyy_1[i] * pa_z[i];

        g_xxzzzz_xyyyyz_0[i] = g_zzzz_xyyyyz_0[i] * fbe_0 - g_zzzz_xyyyyz_1[i] * fz_be_0 + g_xzzzz_yyyyz_1[i] * fe_0 + g_xzzzz_xyyyyz_1[i] * pa_x[i];

        g_xxzzzz_xyyyzz_0[i] = g_zzzz_xyyyzz_0[i] * fbe_0 - g_zzzz_xyyyzz_1[i] * fz_be_0 + g_xzzzz_yyyzz_1[i] * fe_0 + g_xzzzz_xyyyzz_1[i] * pa_x[i];

        g_xxzzzz_xyyzzz_0[i] = g_zzzz_xyyzzz_0[i] * fbe_0 - g_zzzz_xyyzzz_1[i] * fz_be_0 + g_xzzzz_yyzzz_1[i] * fe_0 + g_xzzzz_xyyzzz_1[i] * pa_x[i];

        g_xxzzzz_xyzzzz_0[i] = g_zzzz_xyzzzz_0[i] * fbe_0 - g_zzzz_xyzzzz_1[i] * fz_be_0 + g_xzzzz_yzzzz_1[i] * fe_0 + g_xzzzz_xyzzzz_1[i] * pa_x[i];

        g_xxzzzz_xzzzzz_0[i] = g_zzzz_xzzzzz_0[i] * fbe_0 - g_zzzz_xzzzzz_1[i] * fz_be_0 + g_xzzzz_zzzzz_1[i] * fe_0 + g_xzzzz_xzzzzz_1[i] * pa_x[i];

        g_xxzzzz_yyyyyy_0[i] = g_zzzz_yyyyyy_0[i] * fbe_0 - g_zzzz_yyyyyy_1[i] * fz_be_0 + g_xzzzz_yyyyyy_1[i] * pa_x[i];

        g_xxzzzz_yyyyyz_0[i] = g_zzzz_yyyyyz_0[i] * fbe_0 - g_zzzz_yyyyyz_1[i] * fz_be_0 + g_xzzzz_yyyyyz_1[i] * pa_x[i];

        g_xxzzzz_yyyyzz_0[i] = g_zzzz_yyyyzz_0[i] * fbe_0 - g_zzzz_yyyyzz_1[i] * fz_be_0 + g_xzzzz_yyyyzz_1[i] * pa_x[i];

        g_xxzzzz_yyyzzz_0[i] = g_zzzz_yyyzzz_0[i] * fbe_0 - g_zzzz_yyyzzz_1[i] * fz_be_0 + g_xzzzz_yyyzzz_1[i] * pa_x[i];

        g_xxzzzz_yyzzzz_0[i] = g_zzzz_yyzzzz_0[i] * fbe_0 - g_zzzz_yyzzzz_1[i] * fz_be_0 + g_xzzzz_yyzzzz_1[i] * pa_x[i];

        g_xxzzzz_yzzzzz_0[i] = g_zzzz_yzzzzz_0[i] * fbe_0 - g_zzzz_yzzzzz_1[i] * fz_be_0 + g_xzzzz_yzzzzz_1[i] * pa_x[i];

        g_xxzzzz_zzzzzz_0[i] = g_zzzz_zzzzzz_0[i] * fbe_0 - g_zzzz_zzzzzz_1[i] * fz_be_0 + g_xzzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 420-448 components of targeted buffer : II

    auto g_xyyyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 420);

    auto g_xyyyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 421);

    auto g_xyyyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 422);

    auto g_xyyyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 423);

    auto g_xyyyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 424);

    auto g_xyyyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 425);

    auto g_xyyyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 426);

    auto g_xyyyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 427);

    auto g_xyyyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 428);

    auto g_xyyyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 429);

    auto g_xyyyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 430);

    auto g_xyyyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 431);

    auto g_xyyyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 432);

    auto g_xyyyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 433);

    auto g_xyyyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 434);

    auto g_xyyyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 435);

    auto g_xyyyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 436);

    auto g_xyyyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 437);

    auto g_xyyyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 438);

    auto g_xyyyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 439);

    auto g_xyyyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 440);

    auto g_xyyyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 441);

    auto g_xyyyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 442);

    auto g_xyyyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 443);

    auto g_xyyyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 444);

    auto g_xyyyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 445);

    auto g_xyyyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 446);

    auto g_xyyyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 447);

    #pragma omp simd aligned(g_xyyyyy_xxxxxx_0, g_xyyyyy_xxxxxy_0, g_xyyyyy_xxxxxz_0, g_xyyyyy_xxxxyy_0, g_xyyyyy_xxxxyz_0, g_xyyyyy_xxxxzz_0, g_xyyyyy_xxxyyy_0, g_xyyyyy_xxxyyz_0, g_xyyyyy_xxxyzz_0, g_xyyyyy_xxxzzz_0, g_xyyyyy_xxyyyy_0, g_xyyyyy_xxyyyz_0, g_xyyyyy_xxyyzz_0, g_xyyyyy_xxyzzz_0, g_xyyyyy_xxzzzz_0, g_xyyyyy_xyyyyy_0, g_xyyyyy_xyyyyz_0, g_xyyyyy_xyyyzz_0, g_xyyyyy_xyyzzz_0, g_xyyyyy_xyzzzz_0, g_xyyyyy_xzzzzz_0, g_xyyyyy_yyyyyy_0, g_xyyyyy_yyyyyz_0, g_xyyyyy_yyyyzz_0, g_xyyyyy_yyyzzz_0, g_xyyyyy_yyzzzz_0, g_xyyyyy_yzzzzz_0, g_xyyyyy_zzzzzz_0, g_yyyyy_xxxxx_1, g_yyyyy_xxxxxx_1, g_yyyyy_xxxxxy_1, g_yyyyy_xxxxxz_1, g_yyyyy_xxxxy_1, g_yyyyy_xxxxyy_1, g_yyyyy_xxxxyz_1, g_yyyyy_xxxxz_1, g_yyyyy_xxxxzz_1, g_yyyyy_xxxyy_1, g_yyyyy_xxxyyy_1, g_yyyyy_xxxyyz_1, g_yyyyy_xxxyz_1, g_yyyyy_xxxyzz_1, g_yyyyy_xxxzz_1, g_yyyyy_xxxzzz_1, g_yyyyy_xxyyy_1, g_yyyyy_xxyyyy_1, g_yyyyy_xxyyyz_1, g_yyyyy_xxyyz_1, g_yyyyy_xxyyzz_1, g_yyyyy_xxyzz_1, g_yyyyy_xxyzzz_1, g_yyyyy_xxzzz_1, g_yyyyy_xxzzzz_1, g_yyyyy_xyyyy_1, g_yyyyy_xyyyyy_1, g_yyyyy_xyyyyz_1, g_yyyyy_xyyyz_1, g_yyyyy_xyyyzz_1, g_yyyyy_xyyzz_1, g_yyyyy_xyyzzz_1, g_yyyyy_xyzzz_1, g_yyyyy_xyzzzz_1, g_yyyyy_xzzzz_1, g_yyyyy_xzzzzz_1, g_yyyyy_yyyyy_1, g_yyyyy_yyyyyy_1, g_yyyyy_yyyyyz_1, g_yyyyy_yyyyz_1, g_yyyyy_yyyyzz_1, g_yyyyy_yyyzz_1, g_yyyyy_yyyzzz_1, g_yyyyy_yyzzz_1, g_yyyyy_yyzzzz_1, g_yyyyy_yzzzz_1, g_yyyyy_yzzzzz_1, g_yyyyy_zzzzz_1, g_yyyyy_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyy_xxxxxx_0[i] = 6.0 * g_yyyyy_xxxxx_1[i] * fe_0 + g_yyyyy_xxxxxx_1[i] * pa_x[i];

        g_xyyyyy_xxxxxy_0[i] = 5.0 * g_yyyyy_xxxxy_1[i] * fe_0 + g_yyyyy_xxxxxy_1[i] * pa_x[i];

        g_xyyyyy_xxxxxz_0[i] = 5.0 * g_yyyyy_xxxxz_1[i] * fe_0 + g_yyyyy_xxxxxz_1[i] * pa_x[i];

        g_xyyyyy_xxxxyy_0[i] = 4.0 * g_yyyyy_xxxyy_1[i] * fe_0 + g_yyyyy_xxxxyy_1[i] * pa_x[i];

        g_xyyyyy_xxxxyz_0[i] = 4.0 * g_yyyyy_xxxyz_1[i] * fe_0 + g_yyyyy_xxxxyz_1[i] * pa_x[i];

        g_xyyyyy_xxxxzz_0[i] = 4.0 * g_yyyyy_xxxzz_1[i] * fe_0 + g_yyyyy_xxxxzz_1[i] * pa_x[i];

        g_xyyyyy_xxxyyy_0[i] = 3.0 * g_yyyyy_xxyyy_1[i] * fe_0 + g_yyyyy_xxxyyy_1[i] * pa_x[i];

        g_xyyyyy_xxxyyz_0[i] = 3.0 * g_yyyyy_xxyyz_1[i] * fe_0 + g_yyyyy_xxxyyz_1[i] * pa_x[i];

        g_xyyyyy_xxxyzz_0[i] = 3.0 * g_yyyyy_xxyzz_1[i] * fe_0 + g_yyyyy_xxxyzz_1[i] * pa_x[i];

        g_xyyyyy_xxxzzz_0[i] = 3.0 * g_yyyyy_xxzzz_1[i] * fe_0 + g_yyyyy_xxxzzz_1[i] * pa_x[i];

        g_xyyyyy_xxyyyy_0[i] = 2.0 * g_yyyyy_xyyyy_1[i] * fe_0 + g_yyyyy_xxyyyy_1[i] * pa_x[i];

        g_xyyyyy_xxyyyz_0[i] = 2.0 * g_yyyyy_xyyyz_1[i] * fe_0 + g_yyyyy_xxyyyz_1[i] * pa_x[i];

        g_xyyyyy_xxyyzz_0[i] = 2.0 * g_yyyyy_xyyzz_1[i] * fe_0 + g_yyyyy_xxyyzz_1[i] * pa_x[i];

        g_xyyyyy_xxyzzz_0[i] = 2.0 * g_yyyyy_xyzzz_1[i] * fe_0 + g_yyyyy_xxyzzz_1[i] * pa_x[i];

        g_xyyyyy_xxzzzz_0[i] = 2.0 * g_yyyyy_xzzzz_1[i] * fe_0 + g_yyyyy_xxzzzz_1[i] * pa_x[i];

        g_xyyyyy_xyyyyy_0[i] = g_yyyyy_yyyyy_1[i] * fe_0 + g_yyyyy_xyyyyy_1[i] * pa_x[i];

        g_xyyyyy_xyyyyz_0[i] = g_yyyyy_yyyyz_1[i] * fe_0 + g_yyyyy_xyyyyz_1[i] * pa_x[i];

        g_xyyyyy_xyyyzz_0[i] = g_yyyyy_yyyzz_1[i] * fe_0 + g_yyyyy_xyyyzz_1[i] * pa_x[i];

        g_xyyyyy_xyyzzz_0[i] = g_yyyyy_yyzzz_1[i] * fe_0 + g_yyyyy_xyyzzz_1[i] * pa_x[i];

        g_xyyyyy_xyzzzz_0[i] = g_yyyyy_yzzzz_1[i] * fe_0 + g_yyyyy_xyzzzz_1[i] * pa_x[i];

        g_xyyyyy_xzzzzz_0[i] = g_yyyyy_zzzzz_1[i] * fe_0 + g_yyyyy_xzzzzz_1[i] * pa_x[i];

        g_xyyyyy_yyyyyy_0[i] = g_yyyyy_yyyyyy_1[i] * pa_x[i];

        g_xyyyyy_yyyyyz_0[i] = g_yyyyy_yyyyyz_1[i] * pa_x[i];

        g_xyyyyy_yyyyzz_0[i] = g_yyyyy_yyyyzz_1[i] * pa_x[i];

        g_xyyyyy_yyyzzz_0[i] = g_yyyyy_yyyzzz_1[i] * pa_x[i];

        g_xyyyyy_yyzzzz_0[i] = g_yyyyy_yyzzzz_1[i] * pa_x[i];

        g_xyyyyy_yzzzzz_0[i] = g_yyyyy_yzzzzz_1[i] * pa_x[i];

        g_xyyyyy_zzzzzz_0[i] = g_yyyyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 448-476 components of targeted buffer : II

    auto g_xyyyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 448);

    auto g_xyyyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 449);

    auto g_xyyyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 450);

    auto g_xyyyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 451);

    auto g_xyyyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 452);

    auto g_xyyyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 453);

    auto g_xyyyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 454);

    auto g_xyyyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 455);

    auto g_xyyyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 456);

    auto g_xyyyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 457);

    auto g_xyyyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 458);

    auto g_xyyyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 459);

    auto g_xyyyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 460);

    auto g_xyyyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 461);

    auto g_xyyyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 462);

    auto g_xyyyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 463);

    auto g_xyyyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 464);

    auto g_xyyyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 465);

    auto g_xyyyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 466);

    auto g_xyyyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 467);

    auto g_xyyyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 468);

    auto g_xyyyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 469);

    auto g_xyyyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 470);

    auto g_xyyyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 471);

    auto g_xyyyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 472);

    auto g_xyyyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 473);

    auto g_xyyyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 474);

    auto g_xyyyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 475);

    #pragma omp simd aligned(g_xyyyy_xxxxxx_1, g_xyyyy_xxxxxy_1, g_xyyyy_xxxxyy_1, g_xyyyy_xxxyyy_1, g_xyyyy_xxyyyy_1, g_xyyyy_xyyyyy_1, g_xyyyyz_xxxxxx_0, g_xyyyyz_xxxxxy_0, g_xyyyyz_xxxxxz_0, g_xyyyyz_xxxxyy_0, g_xyyyyz_xxxxyz_0, g_xyyyyz_xxxxzz_0, g_xyyyyz_xxxyyy_0, g_xyyyyz_xxxyyz_0, g_xyyyyz_xxxyzz_0, g_xyyyyz_xxxzzz_0, g_xyyyyz_xxyyyy_0, g_xyyyyz_xxyyyz_0, g_xyyyyz_xxyyzz_0, g_xyyyyz_xxyzzz_0, g_xyyyyz_xxzzzz_0, g_xyyyyz_xyyyyy_0, g_xyyyyz_xyyyyz_0, g_xyyyyz_xyyyzz_0, g_xyyyyz_xyyzzz_0, g_xyyyyz_xyzzzz_0, g_xyyyyz_xzzzzz_0, g_xyyyyz_yyyyyy_0, g_xyyyyz_yyyyyz_0, g_xyyyyz_yyyyzz_0, g_xyyyyz_yyyzzz_0, g_xyyyyz_yyzzzz_0, g_xyyyyz_yzzzzz_0, g_xyyyyz_zzzzzz_0, g_yyyyz_xxxxxz_1, g_yyyyz_xxxxyz_1, g_yyyyz_xxxxz_1, g_yyyyz_xxxxzz_1, g_yyyyz_xxxyyz_1, g_yyyyz_xxxyz_1, g_yyyyz_xxxyzz_1, g_yyyyz_xxxzz_1, g_yyyyz_xxxzzz_1, g_yyyyz_xxyyyz_1, g_yyyyz_xxyyz_1, g_yyyyz_xxyyzz_1, g_yyyyz_xxyzz_1, g_yyyyz_xxyzzz_1, g_yyyyz_xxzzz_1, g_yyyyz_xxzzzz_1, g_yyyyz_xyyyyz_1, g_yyyyz_xyyyz_1, g_yyyyz_xyyyzz_1, g_yyyyz_xyyzz_1, g_yyyyz_xyyzzz_1, g_yyyyz_xyzzz_1, g_yyyyz_xyzzzz_1, g_yyyyz_xzzzz_1, g_yyyyz_xzzzzz_1, g_yyyyz_yyyyyy_1, g_yyyyz_yyyyyz_1, g_yyyyz_yyyyz_1, g_yyyyz_yyyyzz_1, g_yyyyz_yyyzz_1, g_yyyyz_yyyzzz_1, g_yyyyz_yyzzz_1, g_yyyyz_yyzzzz_1, g_yyyyz_yzzzz_1, g_yyyyz_yzzzzz_1, g_yyyyz_zzzzz_1, g_yyyyz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyz_xxxxxx_0[i] = g_xyyyy_xxxxxx_1[i] * pa_z[i];

        g_xyyyyz_xxxxxy_0[i] = g_xyyyy_xxxxxy_1[i] * pa_z[i];

        g_xyyyyz_xxxxxz_0[i] = 5.0 * g_yyyyz_xxxxz_1[i] * fe_0 + g_yyyyz_xxxxxz_1[i] * pa_x[i];

        g_xyyyyz_xxxxyy_0[i] = g_xyyyy_xxxxyy_1[i] * pa_z[i];

        g_xyyyyz_xxxxyz_0[i] = 4.0 * g_yyyyz_xxxyz_1[i] * fe_0 + g_yyyyz_xxxxyz_1[i] * pa_x[i];

        g_xyyyyz_xxxxzz_0[i] = 4.0 * g_yyyyz_xxxzz_1[i] * fe_0 + g_yyyyz_xxxxzz_1[i] * pa_x[i];

        g_xyyyyz_xxxyyy_0[i] = g_xyyyy_xxxyyy_1[i] * pa_z[i];

        g_xyyyyz_xxxyyz_0[i] = 3.0 * g_yyyyz_xxyyz_1[i] * fe_0 + g_yyyyz_xxxyyz_1[i] * pa_x[i];

        g_xyyyyz_xxxyzz_0[i] = 3.0 * g_yyyyz_xxyzz_1[i] * fe_0 + g_yyyyz_xxxyzz_1[i] * pa_x[i];

        g_xyyyyz_xxxzzz_0[i] = 3.0 * g_yyyyz_xxzzz_1[i] * fe_0 + g_yyyyz_xxxzzz_1[i] * pa_x[i];

        g_xyyyyz_xxyyyy_0[i] = g_xyyyy_xxyyyy_1[i] * pa_z[i];

        g_xyyyyz_xxyyyz_0[i] = 2.0 * g_yyyyz_xyyyz_1[i] * fe_0 + g_yyyyz_xxyyyz_1[i] * pa_x[i];

        g_xyyyyz_xxyyzz_0[i] = 2.0 * g_yyyyz_xyyzz_1[i] * fe_0 + g_yyyyz_xxyyzz_1[i] * pa_x[i];

        g_xyyyyz_xxyzzz_0[i] = 2.0 * g_yyyyz_xyzzz_1[i] * fe_0 + g_yyyyz_xxyzzz_1[i] * pa_x[i];

        g_xyyyyz_xxzzzz_0[i] = 2.0 * g_yyyyz_xzzzz_1[i] * fe_0 + g_yyyyz_xxzzzz_1[i] * pa_x[i];

        g_xyyyyz_xyyyyy_0[i] = g_xyyyy_xyyyyy_1[i] * pa_z[i];

        g_xyyyyz_xyyyyz_0[i] = g_yyyyz_yyyyz_1[i] * fe_0 + g_yyyyz_xyyyyz_1[i] * pa_x[i];

        g_xyyyyz_xyyyzz_0[i] = g_yyyyz_yyyzz_1[i] * fe_0 + g_yyyyz_xyyyzz_1[i] * pa_x[i];

        g_xyyyyz_xyyzzz_0[i] = g_yyyyz_yyzzz_1[i] * fe_0 + g_yyyyz_xyyzzz_1[i] * pa_x[i];

        g_xyyyyz_xyzzzz_0[i] = g_yyyyz_yzzzz_1[i] * fe_0 + g_yyyyz_xyzzzz_1[i] * pa_x[i];

        g_xyyyyz_xzzzzz_0[i] = g_yyyyz_zzzzz_1[i] * fe_0 + g_yyyyz_xzzzzz_1[i] * pa_x[i];

        g_xyyyyz_yyyyyy_0[i] = g_yyyyz_yyyyyy_1[i] * pa_x[i];

        g_xyyyyz_yyyyyz_0[i] = g_yyyyz_yyyyyz_1[i] * pa_x[i];

        g_xyyyyz_yyyyzz_0[i] = g_yyyyz_yyyyzz_1[i] * pa_x[i];

        g_xyyyyz_yyyzzz_0[i] = g_yyyyz_yyyzzz_1[i] * pa_x[i];

        g_xyyyyz_yyzzzz_0[i] = g_yyyyz_yyzzzz_1[i] * pa_x[i];

        g_xyyyyz_yzzzzz_0[i] = g_yyyyz_yzzzzz_1[i] * pa_x[i];

        g_xyyyyz_zzzzzz_0[i] = g_yyyyz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 476-504 components of targeted buffer : II

    auto g_xyyyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 476);

    auto g_xyyyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 477);

    auto g_xyyyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 478);

    auto g_xyyyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 479);

    auto g_xyyyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 480);

    auto g_xyyyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 481);

    auto g_xyyyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 482);

    auto g_xyyyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 483);

    auto g_xyyyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 484);

    auto g_xyyyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 485);

    auto g_xyyyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 486);

    auto g_xyyyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 487);

    auto g_xyyyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 488);

    auto g_xyyyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 489);

    auto g_xyyyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 490);

    auto g_xyyyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 491);

    auto g_xyyyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 492);

    auto g_xyyyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 493);

    auto g_xyyyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 494);

    auto g_xyyyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 495);

    auto g_xyyyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 496);

    auto g_xyyyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 497);

    auto g_xyyyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 498);

    auto g_xyyyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 499);

    auto g_xyyyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 500);

    auto g_xyyyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 501);

    auto g_xyyyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 502);

    auto g_xyyyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 503);

    #pragma omp simd aligned(g_xyyyzz_xxxxxx_0, g_xyyyzz_xxxxxy_0, g_xyyyzz_xxxxxz_0, g_xyyyzz_xxxxyy_0, g_xyyyzz_xxxxyz_0, g_xyyyzz_xxxxzz_0, g_xyyyzz_xxxyyy_0, g_xyyyzz_xxxyyz_0, g_xyyyzz_xxxyzz_0, g_xyyyzz_xxxzzz_0, g_xyyyzz_xxyyyy_0, g_xyyyzz_xxyyyz_0, g_xyyyzz_xxyyzz_0, g_xyyyzz_xxyzzz_0, g_xyyyzz_xxzzzz_0, g_xyyyzz_xyyyyy_0, g_xyyyzz_xyyyyz_0, g_xyyyzz_xyyyzz_0, g_xyyyzz_xyyzzz_0, g_xyyyzz_xyzzzz_0, g_xyyyzz_xzzzzz_0, g_xyyyzz_yyyyyy_0, g_xyyyzz_yyyyyz_0, g_xyyyzz_yyyyzz_0, g_xyyyzz_yyyzzz_0, g_xyyyzz_yyzzzz_0, g_xyyyzz_yzzzzz_0, g_xyyyzz_zzzzzz_0, g_yyyzz_xxxxx_1, g_yyyzz_xxxxxx_1, g_yyyzz_xxxxxy_1, g_yyyzz_xxxxxz_1, g_yyyzz_xxxxy_1, g_yyyzz_xxxxyy_1, g_yyyzz_xxxxyz_1, g_yyyzz_xxxxz_1, g_yyyzz_xxxxzz_1, g_yyyzz_xxxyy_1, g_yyyzz_xxxyyy_1, g_yyyzz_xxxyyz_1, g_yyyzz_xxxyz_1, g_yyyzz_xxxyzz_1, g_yyyzz_xxxzz_1, g_yyyzz_xxxzzz_1, g_yyyzz_xxyyy_1, g_yyyzz_xxyyyy_1, g_yyyzz_xxyyyz_1, g_yyyzz_xxyyz_1, g_yyyzz_xxyyzz_1, g_yyyzz_xxyzz_1, g_yyyzz_xxyzzz_1, g_yyyzz_xxzzz_1, g_yyyzz_xxzzzz_1, g_yyyzz_xyyyy_1, g_yyyzz_xyyyyy_1, g_yyyzz_xyyyyz_1, g_yyyzz_xyyyz_1, g_yyyzz_xyyyzz_1, g_yyyzz_xyyzz_1, g_yyyzz_xyyzzz_1, g_yyyzz_xyzzz_1, g_yyyzz_xyzzzz_1, g_yyyzz_xzzzz_1, g_yyyzz_xzzzzz_1, g_yyyzz_yyyyy_1, g_yyyzz_yyyyyy_1, g_yyyzz_yyyyyz_1, g_yyyzz_yyyyz_1, g_yyyzz_yyyyzz_1, g_yyyzz_yyyzz_1, g_yyyzz_yyyzzz_1, g_yyyzz_yyzzz_1, g_yyyzz_yyzzzz_1, g_yyyzz_yzzzz_1, g_yyyzz_yzzzzz_1, g_yyyzz_zzzzz_1, g_yyyzz_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyzz_xxxxxx_0[i] = 6.0 * g_yyyzz_xxxxx_1[i] * fe_0 + g_yyyzz_xxxxxx_1[i] * pa_x[i];

        g_xyyyzz_xxxxxy_0[i] = 5.0 * g_yyyzz_xxxxy_1[i] * fe_0 + g_yyyzz_xxxxxy_1[i] * pa_x[i];

        g_xyyyzz_xxxxxz_0[i] = 5.0 * g_yyyzz_xxxxz_1[i] * fe_0 + g_yyyzz_xxxxxz_1[i] * pa_x[i];

        g_xyyyzz_xxxxyy_0[i] = 4.0 * g_yyyzz_xxxyy_1[i] * fe_0 + g_yyyzz_xxxxyy_1[i] * pa_x[i];

        g_xyyyzz_xxxxyz_0[i] = 4.0 * g_yyyzz_xxxyz_1[i] * fe_0 + g_yyyzz_xxxxyz_1[i] * pa_x[i];

        g_xyyyzz_xxxxzz_0[i] = 4.0 * g_yyyzz_xxxzz_1[i] * fe_0 + g_yyyzz_xxxxzz_1[i] * pa_x[i];

        g_xyyyzz_xxxyyy_0[i] = 3.0 * g_yyyzz_xxyyy_1[i] * fe_0 + g_yyyzz_xxxyyy_1[i] * pa_x[i];

        g_xyyyzz_xxxyyz_0[i] = 3.0 * g_yyyzz_xxyyz_1[i] * fe_0 + g_yyyzz_xxxyyz_1[i] * pa_x[i];

        g_xyyyzz_xxxyzz_0[i] = 3.0 * g_yyyzz_xxyzz_1[i] * fe_0 + g_yyyzz_xxxyzz_1[i] * pa_x[i];

        g_xyyyzz_xxxzzz_0[i] = 3.0 * g_yyyzz_xxzzz_1[i] * fe_0 + g_yyyzz_xxxzzz_1[i] * pa_x[i];

        g_xyyyzz_xxyyyy_0[i] = 2.0 * g_yyyzz_xyyyy_1[i] * fe_0 + g_yyyzz_xxyyyy_1[i] * pa_x[i];

        g_xyyyzz_xxyyyz_0[i] = 2.0 * g_yyyzz_xyyyz_1[i] * fe_0 + g_yyyzz_xxyyyz_1[i] * pa_x[i];

        g_xyyyzz_xxyyzz_0[i] = 2.0 * g_yyyzz_xyyzz_1[i] * fe_0 + g_yyyzz_xxyyzz_1[i] * pa_x[i];

        g_xyyyzz_xxyzzz_0[i] = 2.0 * g_yyyzz_xyzzz_1[i] * fe_0 + g_yyyzz_xxyzzz_1[i] * pa_x[i];

        g_xyyyzz_xxzzzz_0[i] = 2.0 * g_yyyzz_xzzzz_1[i] * fe_0 + g_yyyzz_xxzzzz_1[i] * pa_x[i];

        g_xyyyzz_xyyyyy_0[i] = g_yyyzz_yyyyy_1[i] * fe_0 + g_yyyzz_xyyyyy_1[i] * pa_x[i];

        g_xyyyzz_xyyyyz_0[i] = g_yyyzz_yyyyz_1[i] * fe_0 + g_yyyzz_xyyyyz_1[i] * pa_x[i];

        g_xyyyzz_xyyyzz_0[i] = g_yyyzz_yyyzz_1[i] * fe_0 + g_yyyzz_xyyyzz_1[i] * pa_x[i];

        g_xyyyzz_xyyzzz_0[i] = g_yyyzz_yyzzz_1[i] * fe_0 + g_yyyzz_xyyzzz_1[i] * pa_x[i];

        g_xyyyzz_xyzzzz_0[i] = g_yyyzz_yzzzz_1[i] * fe_0 + g_yyyzz_xyzzzz_1[i] * pa_x[i];

        g_xyyyzz_xzzzzz_0[i] = g_yyyzz_zzzzz_1[i] * fe_0 + g_yyyzz_xzzzzz_1[i] * pa_x[i];

        g_xyyyzz_yyyyyy_0[i] = g_yyyzz_yyyyyy_1[i] * pa_x[i];

        g_xyyyzz_yyyyyz_0[i] = g_yyyzz_yyyyyz_1[i] * pa_x[i];

        g_xyyyzz_yyyyzz_0[i] = g_yyyzz_yyyyzz_1[i] * pa_x[i];

        g_xyyyzz_yyyzzz_0[i] = g_yyyzz_yyyzzz_1[i] * pa_x[i];

        g_xyyyzz_yyzzzz_0[i] = g_yyyzz_yyzzzz_1[i] * pa_x[i];

        g_xyyyzz_yzzzzz_0[i] = g_yyyzz_yzzzzz_1[i] * pa_x[i];

        g_xyyyzz_zzzzzz_0[i] = g_yyyzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 504-532 components of targeted buffer : II

    auto g_xyyzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 504);

    auto g_xyyzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 505);

    auto g_xyyzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 506);

    auto g_xyyzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 507);

    auto g_xyyzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 508);

    auto g_xyyzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 509);

    auto g_xyyzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 510);

    auto g_xyyzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 511);

    auto g_xyyzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 512);

    auto g_xyyzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 513);

    auto g_xyyzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 514);

    auto g_xyyzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 515);

    auto g_xyyzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 516);

    auto g_xyyzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 517);

    auto g_xyyzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 518);

    auto g_xyyzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 519);

    auto g_xyyzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 520);

    auto g_xyyzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 521);

    auto g_xyyzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 522);

    auto g_xyyzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 523);

    auto g_xyyzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 524);

    auto g_xyyzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 525);

    auto g_xyyzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 526);

    auto g_xyyzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 527);

    auto g_xyyzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 528);

    auto g_xyyzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 529);

    auto g_xyyzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 530);

    auto g_xyyzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 531);

    #pragma omp simd aligned(g_xyyzzz_xxxxxx_0, g_xyyzzz_xxxxxy_0, g_xyyzzz_xxxxxz_0, g_xyyzzz_xxxxyy_0, g_xyyzzz_xxxxyz_0, g_xyyzzz_xxxxzz_0, g_xyyzzz_xxxyyy_0, g_xyyzzz_xxxyyz_0, g_xyyzzz_xxxyzz_0, g_xyyzzz_xxxzzz_0, g_xyyzzz_xxyyyy_0, g_xyyzzz_xxyyyz_0, g_xyyzzz_xxyyzz_0, g_xyyzzz_xxyzzz_0, g_xyyzzz_xxzzzz_0, g_xyyzzz_xyyyyy_0, g_xyyzzz_xyyyyz_0, g_xyyzzz_xyyyzz_0, g_xyyzzz_xyyzzz_0, g_xyyzzz_xyzzzz_0, g_xyyzzz_xzzzzz_0, g_xyyzzz_yyyyyy_0, g_xyyzzz_yyyyyz_0, g_xyyzzz_yyyyzz_0, g_xyyzzz_yyyzzz_0, g_xyyzzz_yyzzzz_0, g_xyyzzz_yzzzzz_0, g_xyyzzz_zzzzzz_0, g_yyzzz_xxxxx_1, g_yyzzz_xxxxxx_1, g_yyzzz_xxxxxy_1, g_yyzzz_xxxxxz_1, g_yyzzz_xxxxy_1, g_yyzzz_xxxxyy_1, g_yyzzz_xxxxyz_1, g_yyzzz_xxxxz_1, g_yyzzz_xxxxzz_1, g_yyzzz_xxxyy_1, g_yyzzz_xxxyyy_1, g_yyzzz_xxxyyz_1, g_yyzzz_xxxyz_1, g_yyzzz_xxxyzz_1, g_yyzzz_xxxzz_1, g_yyzzz_xxxzzz_1, g_yyzzz_xxyyy_1, g_yyzzz_xxyyyy_1, g_yyzzz_xxyyyz_1, g_yyzzz_xxyyz_1, g_yyzzz_xxyyzz_1, g_yyzzz_xxyzz_1, g_yyzzz_xxyzzz_1, g_yyzzz_xxzzz_1, g_yyzzz_xxzzzz_1, g_yyzzz_xyyyy_1, g_yyzzz_xyyyyy_1, g_yyzzz_xyyyyz_1, g_yyzzz_xyyyz_1, g_yyzzz_xyyyzz_1, g_yyzzz_xyyzz_1, g_yyzzz_xyyzzz_1, g_yyzzz_xyzzz_1, g_yyzzz_xyzzzz_1, g_yyzzz_xzzzz_1, g_yyzzz_xzzzzz_1, g_yyzzz_yyyyy_1, g_yyzzz_yyyyyy_1, g_yyzzz_yyyyyz_1, g_yyzzz_yyyyz_1, g_yyzzz_yyyyzz_1, g_yyzzz_yyyzz_1, g_yyzzz_yyyzzz_1, g_yyzzz_yyzzz_1, g_yyzzz_yyzzzz_1, g_yyzzz_yzzzz_1, g_yyzzz_yzzzzz_1, g_yyzzz_zzzzz_1, g_yyzzz_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzzz_xxxxxx_0[i] = 6.0 * g_yyzzz_xxxxx_1[i] * fe_0 + g_yyzzz_xxxxxx_1[i] * pa_x[i];

        g_xyyzzz_xxxxxy_0[i] = 5.0 * g_yyzzz_xxxxy_1[i] * fe_0 + g_yyzzz_xxxxxy_1[i] * pa_x[i];

        g_xyyzzz_xxxxxz_0[i] = 5.0 * g_yyzzz_xxxxz_1[i] * fe_0 + g_yyzzz_xxxxxz_1[i] * pa_x[i];

        g_xyyzzz_xxxxyy_0[i] = 4.0 * g_yyzzz_xxxyy_1[i] * fe_0 + g_yyzzz_xxxxyy_1[i] * pa_x[i];

        g_xyyzzz_xxxxyz_0[i] = 4.0 * g_yyzzz_xxxyz_1[i] * fe_0 + g_yyzzz_xxxxyz_1[i] * pa_x[i];

        g_xyyzzz_xxxxzz_0[i] = 4.0 * g_yyzzz_xxxzz_1[i] * fe_0 + g_yyzzz_xxxxzz_1[i] * pa_x[i];

        g_xyyzzz_xxxyyy_0[i] = 3.0 * g_yyzzz_xxyyy_1[i] * fe_0 + g_yyzzz_xxxyyy_1[i] * pa_x[i];

        g_xyyzzz_xxxyyz_0[i] = 3.0 * g_yyzzz_xxyyz_1[i] * fe_0 + g_yyzzz_xxxyyz_1[i] * pa_x[i];

        g_xyyzzz_xxxyzz_0[i] = 3.0 * g_yyzzz_xxyzz_1[i] * fe_0 + g_yyzzz_xxxyzz_1[i] * pa_x[i];

        g_xyyzzz_xxxzzz_0[i] = 3.0 * g_yyzzz_xxzzz_1[i] * fe_0 + g_yyzzz_xxxzzz_1[i] * pa_x[i];

        g_xyyzzz_xxyyyy_0[i] = 2.0 * g_yyzzz_xyyyy_1[i] * fe_0 + g_yyzzz_xxyyyy_1[i] * pa_x[i];

        g_xyyzzz_xxyyyz_0[i] = 2.0 * g_yyzzz_xyyyz_1[i] * fe_0 + g_yyzzz_xxyyyz_1[i] * pa_x[i];

        g_xyyzzz_xxyyzz_0[i] = 2.0 * g_yyzzz_xyyzz_1[i] * fe_0 + g_yyzzz_xxyyzz_1[i] * pa_x[i];

        g_xyyzzz_xxyzzz_0[i] = 2.0 * g_yyzzz_xyzzz_1[i] * fe_0 + g_yyzzz_xxyzzz_1[i] * pa_x[i];

        g_xyyzzz_xxzzzz_0[i] = 2.0 * g_yyzzz_xzzzz_1[i] * fe_0 + g_yyzzz_xxzzzz_1[i] * pa_x[i];

        g_xyyzzz_xyyyyy_0[i] = g_yyzzz_yyyyy_1[i] * fe_0 + g_yyzzz_xyyyyy_1[i] * pa_x[i];

        g_xyyzzz_xyyyyz_0[i] = g_yyzzz_yyyyz_1[i] * fe_0 + g_yyzzz_xyyyyz_1[i] * pa_x[i];

        g_xyyzzz_xyyyzz_0[i] = g_yyzzz_yyyzz_1[i] * fe_0 + g_yyzzz_xyyyzz_1[i] * pa_x[i];

        g_xyyzzz_xyyzzz_0[i] = g_yyzzz_yyzzz_1[i] * fe_0 + g_yyzzz_xyyzzz_1[i] * pa_x[i];

        g_xyyzzz_xyzzzz_0[i] = g_yyzzz_yzzzz_1[i] * fe_0 + g_yyzzz_xyzzzz_1[i] * pa_x[i];

        g_xyyzzz_xzzzzz_0[i] = g_yyzzz_zzzzz_1[i] * fe_0 + g_yyzzz_xzzzzz_1[i] * pa_x[i];

        g_xyyzzz_yyyyyy_0[i] = g_yyzzz_yyyyyy_1[i] * pa_x[i];

        g_xyyzzz_yyyyyz_0[i] = g_yyzzz_yyyyyz_1[i] * pa_x[i];

        g_xyyzzz_yyyyzz_0[i] = g_yyzzz_yyyyzz_1[i] * pa_x[i];

        g_xyyzzz_yyyzzz_0[i] = g_yyzzz_yyyzzz_1[i] * pa_x[i];

        g_xyyzzz_yyzzzz_0[i] = g_yyzzz_yyzzzz_1[i] * pa_x[i];

        g_xyyzzz_yzzzzz_0[i] = g_yyzzz_yzzzzz_1[i] * pa_x[i];

        g_xyyzzz_zzzzzz_0[i] = g_yyzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 532-560 components of targeted buffer : II

    auto g_xyzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 532);

    auto g_xyzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 533);

    auto g_xyzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 534);

    auto g_xyzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 535);

    auto g_xyzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 536);

    auto g_xyzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 537);

    auto g_xyzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 538);

    auto g_xyzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 539);

    auto g_xyzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 540);

    auto g_xyzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 541);

    auto g_xyzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 542);

    auto g_xyzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 543);

    auto g_xyzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 544);

    auto g_xyzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 545);

    auto g_xyzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 546);

    auto g_xyzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 547);

    auto g_xyzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 548);

    auto g_xyzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 549);

    auto g_xyzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 550);

    auto g_xyzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 551);

    auto g_xyzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 552);

    auto g_xyzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 553);

    auto g_xyzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 554);

    auto g_xyzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 555);

    auto g_xyzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 556);

    auto g_xyzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 557);

    auto g_xyzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 558);

    auto g_xyzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 559);

    #pragma omp simd aligned(g_xyzzzz_xxxxxx_0, g_xyzzzz_xxxxxy_0, g_xyzzzz_xxxxxz_0, g_xyzzzz_xxxxyy_0, g_xyzzzz_xxxxyz_0, g_xyzzzz_xxxxzz_0, g_xyzzzz_xxxyyy_0, g_xyzzzz_xxxyyz_0, g_xyzzzz_xxxyzz_0, g_xyzzzz_xxxzzz_0, g_xyzzzz_xxyyyy_0, g_xyzzzz_xxyyyz_0, g_xyzzzz_xxyyzz_0, g_xyzzzz_xxyzzz_0, g_xyzzzz_xxzzzz_0, g_xyzzzz_xyyyyy_0, g_xyzzzz_xyyyyz_0, g_xyzzzz_xyyyzz_0, g_xyzzzz_xyyzzz_0, g_xyzzzz_xyzzzz_0, g_xyzzzz_xzzzzz_0, g_xyzzzz_yyyyyy_0, g_xyzzzz_yyyyyz_0, g_xyzzzz_yyyyzz_0, g_xyzzzz_yyyzzz_0, g_xyzzzz_yyzzzz_0, g_xyzzzz_yzzzzz_0, g_xyzzzz_zzzzzz_0, g_xzzzz_xxxxxx_1, g_xzzzz_xxxxxz_1, g_xzzzz_xxxxzz_1, g_xzzzz_xxxzzz_1, g_xzzzz_xxzzzz_1, g_xzzzz_xzzzzz_1, g_yzzzz_xxxxxy_1, g_yzzzz_xxxxy_1, g_yzzzz_xxxxyy_1, g_yzzzz_xxxxyz_1, g_yzzzz_xxxyy_1, g_yzzzz_xxxyyy_1, g_yzzzz_xxxyyz_1, g_yzzzz_xxxyz_1, g_yzzzz_xxxyzz_1, g_yzzzz_xxyyy_1, g_yzzzz_xxyyyy_1, g_yzzzz_xxyyyz_1, g_yzzzz_xxyyz_1, g_yzzzz_xxyyzz_1, g_yzzzz_xxyzz_1, g_yzzzz_xxyzzz_1, g_yzzzz_xyyyy_1, g_yzzzz_xyyyyy_1, g_yzzzz_xyyyyz_1, g_yzzzz_xyyyz_1, g_yzzzz_xyyyzz_1, g_yzzzz_xyyzz_1, g_yzzzz_xyyzzz_1, g_yzzzz_xyzzz_1, g_yzzzz_xyzzzz_1, g_yzzzz_yyyyy_1, g_yzzzz_yyyyyy_1, g_yzzzz_yyyyyz_1, g_yzzzz_yyyyz_1, g_yzzzz_yyyyzz_1, g_yzzzz_yyyzz_1, g_yzzzz_yyyzzz_1, g_yzzzz_yyzzz_1, g_yzzzz_yyzzzz_1, g_yzzzz_yzzzz_1, g_yzzzz_yzzzzz_1, g_yzzzz_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzzz_xxxxxx_0[i] = g_xzzzz_xxxxxx_1[i] * pa_y[i];

        g_xyzzzz_xxxxxy_0[i] = 5.0 * g_yzzzz_xxxxy_1[i] * fe_0 + g_yzzzz_xxxxxy_1[i] * pa_x[i];

        g_xyzzzz_xxxxxz_0[i] = g_xzzzz_xxxxxz_1[i] * pa_y[i];

        g_xyzzzz_xxxxyy_0[i] = 4.0 * g_yzzzz_xxxyy_1[i] * fe_0 + g_yzzzz_xxxxyy_1[i] * pa_x[i];

        g_xyzzzz_xxxxyz_0[i] = 4.0 * g_yzzzz_xxxyz_1[i] * fe_0 + g_yzzzz_xxxxyz_1[i] * pa_x[i];

        g_xyzzzz_xxxxzz_0[i] = g_xzzzz_xxxxzz_1[i] * pa_y[i];

        g_xyzzzz_xxxyyy_0[i] = 3.0 * g_yzzzz_xxyyy_1[i] * fe_0 + g_yzzzz_xxxyyy_1[i] * pa_x[i];

        g_xyzzzz_xxxyyz_0[i] = 3.0 * g_yzzzz_xxyyz_1[i] * fe_0 + g_yzzzz_xxxyyz_1[i] * pa_x[i];

        g_xyzzzz_xxxyzz_0[i] = 3.0 * g_yzzzz_xxyzz_1[i] * fe_0 + g_yzzzz_xxxyzz_1[i] * pa_x[i];

        g_xyzzzz_xxxzzz_0[i] = g_xzzzz_xxxzzz_1[i] * pa_y[i];

        g_xyzzzz_xxyyyy_0[i] = 2.0 * g_yzzzz_xyyyy_1[i] * fe_0 + g_yzzzz_xxyyyy_1[i] * pa_x[i];

        g_xyzzzz_xxyyyz_0[i] = 2.0 * g_yzzzz_xyyyz_1[i] * fe_0 + g_yzzzz_xxyyyz_1[i] * pa_x[i];

        g_xyzzzz_xxyyzz_0[i] = 2.0 * g_yzzzz_xyyzz_1[i] * fe_0 + g_yzzzz_xxyyzz_1[i] * pa_x[i];

        g_xyzzzz_xxyzzz_0[i] = 2.0 * g_yzzzz_xyzzz_1[i] * fe_0 + g_yzzzz_xxyzzz_1[i] * pa_x[i];

        g_xyzzzz_xxzzzz_0[i] = g_xzzzz_xxzzzz_1[i] * pa_y[i];

        g_xyzzzz_xyyyyy_0[i] = g_yzzzz_yyyyy_1[i] * fe_0 + g_yzzzz_xyyyyy_1[i] * pa_x[i];

        g_xyzzzz_xyyyyz_0[i] = g_yzzzz_yyyyz_1[i] * fe_0 + g_yzzzz_xyyyyz_1[i] * pa_x[i];

        g_xyzzzz_xyyyzz_0[i] = g_yzzzz_yyyzz_1[i] * fe_0 + g_yzzzz_xyyyzz_1[i] * pa_x[i];

        g_xyzzzz_xyyzzz_0[i] = g_yzzzz_yyzzz_1[i] * fe_0 + g_yzzzz_xyyzzz_1[i] * pa_x[i];

        g_xyzzzz_xyzzzz_0[i] = g_yzzzz_yzzzz_1[i] * fe_0 + g_yzzzz_xyzzzz_1[i] * pa_x[i];

        g_xyzzzz_xzzzzz_0[i] = g_xzzzz_xzzzzz_1[i] * pa_y[i];

        g_xyzzzz_yyyyyy_0[i] = g_yzzzz_yyyyyy_1[i] * pa_x[i];

        g_xyzzzz_yyyyyz_0[i] = g_yzzzz_yyyyyz_1[i] * pa_x[i];

        g_xyzzzz_yyyyzz_0[i] = g_yzzzz_yyyyzz_1[i] * pa_x[i];

        g_xyzzzz_yyyzzz_0[i] = g_yzzzz_yyyzzz_1[i] * pa_x[i];

        g_xyzzzz_yyzzzz_0[i] = g_yzzzz_yyzzzz_1[i] * pa_x[i];

        g_xyzzzz_yzzzzz_0[i] = g_yzzzz_yzzzzz_1[i] * pa_x[i];

        g_xyzzzz_zzzzzz_0[i] = g_yzzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 560-588 components of targeted buffer : II

    auto g_xzzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 560);

    auto g_xzzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 561);

    auto g_xzzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 562);

    auto g_xzzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 563);

    auto g_xzzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 564);

    auto g_xzzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 565);

    auto g_xzzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 566);

    auto g_xzzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 567);

    auto g_xzzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 568);

    auto g_xzzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 569);

    auto g_xzzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 570);

    auto g_xzzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 571);

    auto g_xzzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 572);

    auto g_xzzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 573);

    auto g_xzzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 574);

    auto g_xzzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 575);

    auto g_xzzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 576);

    auto g_xzzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 577);

    auto g_xzzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 578);

    auto g_xzzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 579);

    auto g_xzzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 580);

    auto g_xzzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 581);

    auto g_xzzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 582);

    auto g_xzzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 583);

    auto g_xzzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 584);

    auto g_xzzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 585);

    auto g_xzzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 586);

    auto g_xzzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 587);

    #pragma omp simd aligned(g_xzzzzz_xxxxxx_0, g_xzzzzz_xxxxxy_0, g_xzzzzz_xxxxxz_0, g_xzzzzz_xxxxyy_0, g_xzzzzz_xxxxyz_0, g_xzzzzz_xxxxzz_0, g_xzzzzz_xxxyyy_0, g_xzzzzz_xxxyyz_0, g_xzzzzz_xxxyzz_0, g_xzzzzz_xxxzzz_0, g_xzzzzz_xxyyyy_0, g_xzzzzz_xxyyyz_0, g_xzzzzz_xxyyzz_0, g_xzzzzz_xxyzzz_0, g_xzzzzz_xxzzzz_0, g_xzzzzz_xyyyyy_0, g_xzzzzz_xyyyyz_0, g_xzzzzz_xyyyzz_0, g_xzzzzz_xyyzzz_0, g_xzzzzz_xyzzzz_0, g_xzzzzz_xzzzzz_0, g_xzzzzz_yyyyyy_0, g_xzzzzz_yyyyyz_0, g_xzzzzz_yyyyzz_0, g_xzzzzz_yyyzzz_0, g_xzzzzz_yyzzzz_0, g_xzzzzz_yzzzzz_0, g_xzzzzz_zzzzzz_0, g_zzzzz_xxxxx_1, g_zzzzz_xxxxxx_1, g_zzzzz_xxxxxy_1, g_zzzzz_xxxxxz_1, g_zzzzz_xxxxy_1, g_zzzzz_xxxxyy_1, g_zzzzz_xxxxyz_1, g_zzzzz_xxxxz_1, g_zzzzz_xxxxzz_1, g_zzzzz_xxxyy_1, g_zzzzz_xxxyyy_1, g_zzzzz_xxxyyz_1, g_zzzzz_xxxyz_1, g_zzzzz_xxxyzz_1, g_zzzzz_xxxzz_1, g_zzzzz_xxxzzz_1, g_zzzzz_xxyyy_1, g_zzzzz_xxyyyy_1, g_zzzzz_xxyyyz_1, g_zzzzz_xxyyz_1, g_zzzzz_xxyyzz_1, g_zzzzz_xxyzz_1, g_zzzzz_xxyzzz_1, g_zzzzz_xxzzz_1, g_zzzzz_xxzzzz_1, g_zzzzz_xyyyy_1, g_zzzzz_xyyyyy_1, g_zzzzz_xyyyyz_1, g_zzzzz_xyyyz_1, g_zzzzz_xyyyzz_1, g_zzzzz_xyyzz_1, g_zzzzz_xyyzzz_1, g_zzzzz_xyzzz_1, g_zzzzz_xyzzzz_1, g_zzzzz_xzzzz_1, g_zzzzz_xzzzzz_1, g_zzzzz_yyyyy_1, g_zzzzz_yyyyyy_1, g_zzzzz_yyyyyz_1, g_zzzzz_yyyyz_1, g_zzzzz_yyyyzz_1, g_zzzzz_yyyzz_1, g_zzzzz_yyyzzz_1, g_zzzzz_yyzzz_1, g_zzzzz_yyzzzz_1, g_zzzzz_yzzzz_1, g_zzzzz_yzzzzz_1, g_zzzzz_zzzzz_1, g_zzzzz_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzzz_xxxxxx_0[i] = 6.0 * g_zzzzz_xxxxx_1[i] * fe_0 + g_zzzzz_xxxxxx_1[i] * pa_x[i];

        g_xzzzzz_xxxxxy_0[i] = 5.0 * g_zzzzz_xxxxy_1[i] * fe_0 + g_zzzzz_xxxxxy_1[i] * pa_x[i];

        g_xzzzzz_xxxxxz_0[i] = 5.0 * g_zzzzz_xxxxz_1[i] * fe_0 + g_zzzzz_xxxxxz_1[i] * pa_x[i];

        g_xzzzzz_xxxxyy_0[i] = 4.0 * g_zzzzz_xxxyy_1[i] * fe_0 + g_zzzzz_xxxxyy_1[i] * pa_x[i];

        g_xzzzzz_xxxxyz_0[i] = 4.0 * g_zzzzz_xxxyz_1[i] * fe_0 + g_zzzzz_xxxxyz_1[i] * pa_x[i];

        g_xzzzzz_xxxxzz_0[i] = 4.0 * g_zzzzz_xxxzz_1[i] * fe_0 + g_zzzzz_xxxxzz_1[i] * pa_x[i];

        g_xzzzzz_xxxyyy_0[i] = 3.0 * g_zzzzz_xxyyy_1[i] * fe_0 + g_zzzzz_xxxyyy_1[i] * pa_x[i];

        g_xzzzzz_xxxyyz_0[i] = 3.0 * g_zzzzz_xxyyz_1[i] * fe_0 + g_zzzzz_xxxyyz_1[i] * pa_x[i];

        g_xzzzzz_xxxyzz_0[i] = 3.0 * g_zzzzz_xxyzz_1[i] * fe_0 + g_zzzzz_xxxyzz_1[i] * pa_x[i];

        g_xzzzzz_xxxzzz_0[i] = 3.0 * g_zzzzz_xxzzz_1[i] * fe_0 + g_zzzzz_xxxzzz_1[i] * pa_x[i];

        g_xzzzzz_xxyyyy_0[i] = 2.0 * g_zzzzz_xyyyy_1[i] * fe_0 + g_zzzzz_xxyyyy_1[i] * pa_x[i];

        g_xzzzzz_xxyyyz_0[i] = 2.0 * g_zzzzz_xyyyz_1[i] * fe_0 + g_zzzzz_xxyyyz_1[i] * pa_x[i];

        g_xzzzzz_xxyyzz_0[i] = 2.0 * g_zzzzz_xyyzz_1[i] * fe_0 + g_zzzzz_xxyyzz_1[i] * pa_x[i];

        g_xzzzzz_xxyzzz_0[i] = 2.0 * g_zzzzz_xyzzz_1[i] * fe_0 + g_zzzzz_xxyzzz_1[i] * pa_x[i];

        g_xzzzzz_xxzzzz_0[i] = 2.0 * g_zzzzz_xzzzz_1[i] * fe_0 + g_zzzzz_xxzzzz_1[i] * pa_x[i];

        g_xzzzzz_xyyyyy_0[i] = g_zzzzz_yyyyy_1[i] * fe_0 + g_zzzzz_xyyyyy_1[i] * pa_x[i];

        g_xzzzzz_xyyyyz_0[i] = g_zzzzz_yyyyz_1[i] * fe_0 + g_zzzzz_xyyyyz_1[i] * pa_x[i];

        g_xzzzzz_xyyyzz_0[i] = g_zzzzz_yyyzz_1[i] * fe_0 + g_zzzzz_xyyyzz_1[i] * pa_x[i];

        g_xzzzzz_xyyzzz_0[i] = g_zzzzz_yyzzz_1[i] * fe_0 + g_zzzzz_xyyzzz_1[i] * pa_x[i];

        g_xzzzzz_xyzzzz_0[i] = g_zzzzz_yzzzz_1[i] * fe_0 + g_zzzzz_xyzzzz_1[i] * pa_x[i];

        g_xzzzzz_xzzzzz_0[i] = g_zzzzz_zzzzz_1[i] * fe_0 + g_zzzzz_xzzzzz_1[i] * pa_x[i];

        g_xzzzzz_yyyyyy_0[i] = g_zzzzz_yyyyyy_1[i] * pa_x[i];

        g_xzzzzz_yyyyyz_0[i] = g_zzzzz_yyyyyz_1[i] * pa_x[i];

        g_xzzzzz_yyyyzz_0[i] = g_zzzzz_yyyyzz_1[i] * pa_x[i];

        g_xzzzzz_yyyzzz_0[i] = g_zzzzz_yyyzzz_1[i] * pa_x[i];

        g_xzzzzz_yyzzzz_0[i] = g_zzzzz_yyzzzz_1[i] * pa_x[i];

        g_xzzzzz_yzzzzz_0[i] = g_zzzzz_yzzzzz_1[i] * pa_x[i];

        g_xzzzzz_zzzzzz_0[i] = g_zzzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 588-616 components of targeted buffer : II

    auto g_yyyyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 588);

    auto g_yyyyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 589);

    auto g_yyyyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 590);

    auto g_yyyyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 591);

    auto g_yyyyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 592);

    auto g_yyyyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 593);

    auto g_yyyyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 594);

    auto g_yyyyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 595);

    auto g_yyyyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 596);

    auto g_yyyyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 597);

    auto g_yyyyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 598);

    auto g_yyyyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 599);

    auto g_yyyyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 600);

    auto g_yyyyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 601);

    auto g_yyyyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 602);

    auto g_yyyyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 603);

    auto g_yyyyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 604);

    auto g_yyyyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 605);

    auto g_yyyyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 606);

    auto g_yyyyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 607);

    auto g_yyyyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 608);

    auto g_yyyyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 609);

    auto g_yyyyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 610);

    auto g_yyyyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 611);

    auto g_yyyyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 612);

    auto g_yyyyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 613);

    auto g_yyyyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 614);

    auto g_yyyyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 615);

    #pragma omp simd aligned(g_yyyy_xxxxxx_0, g_yyyy_xxxxxx_1, g_yyyy_xxxxxy_0, g_yyyy_xxxxxy_1, g_yyyy_xxxxxz_0, g_yyyy_xxxxxz_1, g_yyyy_xxxxyy_0, g_yyyy_xxxxyy_1, g_yyyy_xxxxyz_0, g_yyyy_xxxxyz_1, g_yyyy_xxxxzz_0, g_yyyy_xxxxzz_1, g_yyyy_xxxyyy_0, g_yyyy_xxxyyy_1, g_yyyy_xxxyyz_0, g_yyyy_xxxyyz_1, g_yyyy_xxxyzz_0, g_yyyy_xxxyzz_1, g_yyyy_xxxzzz_0, g_yyyy_xxxzzz_1, g_yyyy_xxyyyy_0, g_yyyy_xxyyyy_1, g_yyyy_xxyyyz_0, g_yyyy_xxyyyz_1, g_yyyy_xxyyzz_0, g_yyyy_xxyyzz_1, g_yyyy_xxyzzz_0, g_yyyy_xxyzzz_1, g_yyyy_xxzzzz_0, g_yyyy_xxzzzz_1, g_yyyy_xyyyyy_0, g_yyyy_xyyyyy_1, g_yyyy_xyyyyz_0, g_yyyy_xyyyyz_1, g_yyyy_xyyyzz_0, g_yyyy_xyyyzz_1, g_yyyy_xyyzzz_0, g_yyyy_xyyzzz_1, g_yyyy_xyzzzz_0, g_yyyy_xyzzzz_1, g_yyyy_xzzzzz_0, g_yyyy_xzzzzz_1, g_yyyy_yyyyyy_0, g_yyyy_yyyyyy_1, g_yyyy_yyyyyz_0, g_yyyy_yyyyyz_1, g_yyyy_yyyyzz_0, g_yyyy_yyyyzz_1, g_yyyy_yyyzzz_0, g_yyyy_yyyzzz_1, g_yyyy_yyzzzz_0, g_yyyy_yyzzzz_1, g_yyyy_yzzzzz_0, g_yyyy_yzzzzz_1, g_yyyy_zzzzzz_0, g_yyyy_zzzzzz_1, g_yyyyy_xxxxx_1, g_yyyyy_xxxxxx_1, g_yyyyy_xxxxxy_1, g_yyyyy_xxxxxz_1, g_yyyyy_xxxxy_1, g_yyyyy_xxxxyy_1, g_yyyyy_xxxxyz_1, g_yyyyy_xxxxz_1, g_yyyyy_xxxxzz_1, g_yyyyy_xxxyy_1, g_yyyyy_xxxyyy_1, g_yyyyy_xxxyyz_1, g_yyyyy_xxxyz_1, g_yyyyy_xxxyzz_1, g_yyyyy_xxxzz_1, g_yyyyy_xxxzzz_1, g_yyyyy_xxyyy_1, g_yyyyy_xxyyyy_1, g_yyyyy_xxyyyz_1, g_yyyyy_xxyyz_1, g_yyyyy_xxyyzz_1, g_yyyyy_xxyzz_1, g_yyyyy_xxyzzz_1, g_yyyyy_xxzzz_1, g_yyyyy_xxzzzz_1, g_yyyyy_xyyyy_1, g_yyyyy_xyyyyy_1, g_yyyyy_xyyyyz_1, g_yyyyy_xyyyz_1, g_yyyyy_xyyyzz_1, g_yyyyy_xyyzz_1, g_yyyyy_xyyzzz_1, g_yyyyy_xyzzz_1, g_yyyyy_xyzzzz_1, g_yyyyy_xzzzz_1, g_yyyyy_xzzzzz_1, g_yyyyy_yyyyy_1, g_yyyyy_yyyyyy_1, g_yyyyy_yyyyyz_1, g_yyyyy_yyyyz_1, g_yyyyy_yyyyzz_1, g_yyyyy_yyyzz_1, g_yyyyy_yyyzzz_1, g_yyyyy_yyzzz_1, g_yyyyy_yyzzzz_1, g_yyyyy_yzzzz_1, g_yyyyy_yzzzzz_1, g_yyyyy_zzzzz_1, g_yyyyy_zzzzzz_1, g_yyyyyy_xxxxxx_0, g_yyyyyy_xxxxxy_0, g_yyyyyy_xxxxxz_0, g_yyyyyy_xxxxyy_0, g_yyyyyy_xxxxyz_0, g_yyyyyy_xxxxzz_0, g_yyyyyy_xxxyyy_0, g_yyyyyy_xxxyyz_0, g_yyyyyy_xxxyzz_0, g_yyyyyy_xxxzzz_0, g_yyyyyy_xxyyyy_0, g_yyyyyy_xxyyyz_0, g_yyyyyy_xxyyzz_0, g_yyyyyy_xxyzzz_0, g_yyyyyy_xxzzzz_0, g_yyyyyy_xyyyyy_0, g_yyyyyy_xyyyyz_0, g_yyyyyy_xyyyzz_0, g_yyyyyy_xyyzzz_0, g_yyyyyy_xyzzzz_0, g_yyyyyy_xzzzzz_0, g_yyyyyy_yyyyyy_0, g_yyyyyy_yyyyyz_0, g_yyyyyy_yyyyzz_0, g_yyyyyy_yyyzzz_0, g_yyyyyy_yyzzzz_0, g_yyyyyy_yzzzzz_0, g_yyyyyy_zzzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyy_xxxxxx_0[i] = 5.0 * g_yyyy_xxxxxx_0[i] * fbe_0 - 5.0 * g_yyyy_xxxxxx_1[i] * fz_be_0 + g_yyyyy_xxxxxx_1[i] * pa_y[i];

        g_yyyyyy_xxxxxy_0[i] = 5.0 * g_yyyy_xxxxxy_0[i] * fbe_0 - 5.0 * g_yyyy_xxxxxy_1[i] * fz_be_0 + g_yyyyy_xxxxx_1[i] * fe_0 + g_yyyyy_xxxxxy_1[i] * pa_y[i];

        g_yyyyyy_xxxxxz_0[i] = 5.0 * g_yyyy_xxxxxz_0[i] * fbe_0 - 5.0 * g_yyyy_xxxxxz_1[i] * fz_be_0 + g_yyyyy_xxxxxz_1[i] * pa_y[i];

        g_yyyyyy_xxxxyy_0[i] = 5.0 * g_yyyy_xxxxyy_0[i] * fbe_0 - 5.0 * g_yyyy_xxxxyy_1[i] * fz_be_0 + 2.0 * g_yyyyy_xxxxy_1[i] * fe_0 + g_yyyyy_xxxxyy_1[i] * pa_y[i];

        g_yyyyyy_xxxxyz_0[i] = 5.0 * g_yyyy_xxxxyz_0[i] * fbe_0 - 5.0 * g_yyyy_xxxxyz_1[i] * fz_be_0 + g_yyyyy_xxxxz_1[i] * fe_0 + g_yyyyy_xxxxyz_1[i] * pa_y[i];

        g_yyyyyy_xxxxzz_0[i] = 5.0 * g_yyyy_xxxxzz_0[i] * fbe_0 - 5.0 * g_yyyy_xxxxzz_1[i] * fz_be_0 + g_yyyyy_xxxxzz_1[i] * pa_y[i];

        g_yyyyyy_xxxyyy_0[i] = 5.0 * g_yyyy_xxxyyy_0[i] * fbe_0 - 5.0 * g_yyyy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_yyyyy_xxxyy_1[i] * fe_0 + g_yyyyy_xxxyyy_1[i] * pa_y[i];

        g_yyyyyy_xxxyyz_0[i] = 5.0 * g_yyyy_xxxyyz_0[i] * fbe_0 - 5.0 * g_yyyy_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyy_xxxyz_1[i] * fe_0 + g_yyyyy_xxxyyz_1[i] * pa_y[i];

        g_yyyyyy_xxxyzz_0[i] = 5.0 * g_yyyy_xxxyzz_0[i] * fbe_0 - 5.0 * g_yyyy_xxxyzz_1[i] * fz_be_0 + g_yyyyy_xxxzz_1[i] * fe_0 + g_yyyyy_xxxyzz_1[i] * pa_y[i];

        g_yyyyyy_xxxzzz_0[i] = 5.0 * g_yyyy_xxxzzz_0[i] * fbe_0 - 5.0 * g_yyyy_xxxzzz_1[i] * fz_be_0 + g_yyyyy_xxxzzz_1[i] * pa_y[i];

        g_yyyyyy_xxyyyy_0[i] = 5.0 * g_yyyy_xxyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_xxyyyy_1[i] * fz_be_0 + 4.0 * g_yyyyy_xxyyy_1[i] * fe_0 + g_yyyyy_xxyyyy_1[i] * pa_y[i];

        g_yyyyyy_xxyyyz_0[i] = 5.0 * g_yyyy_xxyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyy_xxyyz_1[i] * fe_0 + g_yyyyy_xxyyyz_1[i] * pa_y[i];

        g_yyyyyy_xxyyzz_0[i] = 5.0 * g_yyyy_xxyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_xxyzz_1[i] * fe_0 + g_yyyyy_xxyyzz_1[i] * pa_y[i];

        g_yyyyyy_xxyzzz_0[i] = 5.0 * g_yyyy_xxyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_xxyzzz_1[i] * fz_be_0 + g_yyyyy_xxzzz_1[i] * fe_0 + g_yyyyy_xxyzzz_1[i] * pa_y[i];

        g_yyyyyy_xxzzzz_0[i] = 5.0 * g_yyyy_xxzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_xxzzzz_1[i] * fz_be_0 + g_yyyyy_xxzzzz_1[i] * pa_y[i];

        g_yyyyyy_xyyyyy_0[i] = 5.0 * g_yyyy_xyyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_xyyyyy_1[i] * fz_be_0 + 5.0 * g_yyyyy_xyyyy_1[i] * fe_0 + g_yyyyy_xyyyyy_1[i] * pa_y[i];

        g_yyyyyy_xyyyyz_0[i] = 5.0 * g_yyyy_xyyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyy_xyyyz_1[i] * fe_0 + g_yyyyy_xyyyyz_1[i] * pa_y[i];

        g_yyyyyy_xyyyzz_0[i] = 5.0 * g_yyyy_xyyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_xyyzz_1[i] * fe_0 + g_yyyyy_xyyyzz_1[i] * pa_y[i];

        g_yyyyyy_xyyzzz_0[i] = 5.0 * g_yyyy_xyyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_xyzzz_1[i] * fe_0 + g_yyyyy_xyyzzz_1[i] * pa_y[i];

        g_yyyyyy_xyzzzz_0[i] = 5.0 * g_yyyy_xyzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_xyzzzz_1[i] * fz_be_0 + g_yyyyy_xzzzz_1[i] * fe_0 + g_yyyyy_xyzzzz_1[i] * pa_y[i];

        g_yyyyyy_xzzzzz_0[i] = 5.0 * g_yyyy_xzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_xzzzzz_1[i] * fz_be_0 + g_yyyyy_xzzzzz_1[i] * pa_y[i];

        g_yyyyyy_yyyyyy_0[i] = 5.0 * g_yyyy_yyyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_yyyyyy_1[i] * fz_be_0 + 6.0 * g_yyyyy_yyyyy_1[i] * fe_0 + g_yyyyy_yyyyyy_1[i] * pa_y[i];

        g_yyyyyy_yyyyyz_0[i] = 5.0 * g_yyyy_yyyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyyy_yyyyz_1[i] * fe_0 + g_yyyyy_yyyyyz_1[i] * pa_y[i];

        g_yyyyyy_yyyyzz_0[i] = 5.0 * g_yyyy_yyyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyyy_yyyzz_1[i] * fe_0 + g_yyyyy_yyyyzz_1[i] * pa_y[i];

        g_yyyyyy_yyyzzz_0[i] = 5.0 * g_yyyy_yyyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_yyzzz_1[i] * fe_0 + g_yyyyy_yyyzzz_1[i] * pa_y[i];

        g_yyyyyy_yyzzzz_0[i] = 5.0 * g_yyyy_yyzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_yzzzz_1[i] * fe_0 + g_yyyyy_yyzzzz_1[i] * pa_y[i];

        g_yyyyyy_yzzzzz_0[i] = 5.0 * g_yyyy_yzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_yzzzzz_1[i] * fz_be_0 + g_yyyyy_zzzzz_1[i] * fe_0 + g_yyyyy_yzzzzz_1[i] * pa_y[i];

        g_yyyyyy_zzzzzz_0[i] = 5.0 * g_yyyy_zzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_zzzzzz_1[i] * fz_be_0 + g_yyyyy_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 616-644 components of targeted buffer : II

    auto g_yyyyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 616);

    auto g_yyyyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 617);

    auto g_yyyyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 618);

    auto g_yyyyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 619);

    auto g_yyyyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 620);

    auto g_yyyyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 621);

    auto g_yyyyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 622);

    auto g_yyyyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 623);

    auto g_yyyyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 624);

    auto g_yyyyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 625);

    auto g_yyyyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 626);

    auto g_yyyyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 627);

    auto g_yyyyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 628);

    auto g_yyyyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 629);

    auto g_yyyyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 630);

    auto g_yyyyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 631);

    auto g_yyyyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 632);

    auto g_yyyyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 633);

    auto g_yyyyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 634);

    auto g_yyyyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 635);

    auto g_yyyyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 636);

    auto g_yyyyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 637);

    auto g_yyyyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 638);

    auto g_yyyyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 639);

    auto g_yyyyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 640);

    auto g_yyyyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 641);

    auto g_yyyyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 642);

    auto g_yyyyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 643);

    #pragma omp simd aligned(g_yyyyy_xxxxx_1, g_yyyyy_xxxxxx_1, g_yyyyy_xxxxxy_1, g_yyyyy_xxxxxz_1, g_yyyyy_xxxxy_1, g_yyyyy_xxxxyy_1, g_yyyyy_xxxxyz_1, g_yyyyy_xxxxz_1, g_yyyyy_xxxxzz_1, g_yyyyy_xxxyy_1, g_yyyyy_xxxyyy_1, g_yyyyy_xxxyyz_1, g_yyyyy_xxxyz_1, g_yyyyy_xxxyzz_1, g_yyyyy_xxxzz_1, g_yyyyy_xxxzzz_1, g_yyyyy_xxyyy_1, g_yyyyy_xxyyyy_1, g_yyyyy_xxyyyz_1, g_yyyyy_xxyyz_1, g_yyyyy_xxyyzz_1, g_yyyyy_xxyzz_1, g_yyyyy_xxyzzz_1, g_yyyyy_xxzzz_1, g_yyyyy_xxzzzz_1, g_yyyyy_xyyyy_1, g_yyyyy_xyyyyy_1, g_yyyyy_xyyyyz_1, g_yyyyy_xyyyz_1, g_yyyyy_xyyyzz_1, g_yyyyy_xyyzz_1, g_yyyyy_xyyzzz_1, g_yyyyy_xyzzz_1, g_yyyyy_xyzzzz_1, g_yyyyy_xzzzz_1, g_yyyyy_xzzzzz_1, g_yyyyy_yyyyy_1, g_yyyyy_yyyyyy_1, g_yyyyy_yyyyyz_1, g_yyyyy_yyyyz_1, g_yyyyy_yyyyzz_1, g_yyyyy_yyyzz_1, g_yyyyy_yyyzzz_1, g_yyyyy_yyzzz_1, g_yyyyy_yyzzzz_1, g_yyyyy_yzzzz_1, g_yyyyy_yzzzzz_1, g_yyyyy_zzzzz_1, g_yyyyy_zzzzzz_1, g_yyyyyz_xxxxxx_0, g_yyyyyz_xxxxxy_0, g_yyyyyz_xxxxxz_0, g_yyyyyz_xxxxyy_0, g_yyyyyz_xxxxyz_0, g_yyyyyz_xxxxzz_0, g_yyyyyz_xxxyyy_0, g_yyyyyz_xxxyyz_0, g_yyyyyz_xxxyzz_0, g_yyyyyz_xxxzzz_0, g_yyyyyz_xxyyyy_0, g_yyyyyz_xxyyyz_0, g_yyyyyz_xxyyzz_0, g_yyyyyz_xxyzzz_0, g_yyyyyz_xxzzzz_0, g_yyyyyz_xyyyyy_0, g_yyyyyz_xyyyyz_0, g_yyyyyz_xyyyzz_0, g_yyyyyz_xyyzzz_0, g_yyyyyz_xyzzzz_0, g_yyyyyz_xzzzzz_0, g_yyyyyz_yyyyyy_0, g_yyyyyz_yyyyyz_0, g_yyyyyz_yyyyzz_0, g_yyyyyz_yyyzzz_0, g_yyyyyz_yyzzzz_0, g_yyyyyz_yzzzzz_0, g_yyyyyz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyyz_xxxxxx_0[i] = g_yyyyy_xxxxxx_1[i] * pa_z[i];

        g_yyyyyz_xxxxxy_0[i] = g_yyyyy_xxxxxy_1[i] * pa_z[i];

        g_yyyyyz_xxxxxz_0[i] = g_yyyyy_xxxxx_1[i] * fe_0 + g_yyyyy_xxxxxz_1[i] * pa_z[i];

        g_yyyyyz_xxxxyy_0[i] = g_yyyyy_xxxxyy_1[i] * pa_z[i];

        g_yyyyyz_xxxxyz_0[i] = g_yyyyy_xxxxy_1[i] * fe_0 + g_yyyyy_xxxxyz_1[i] * pa_z[i];

        g_yyyyyz_xxxxzz_0[i] = 2.0 * g_yyyyy_xxxxz_1[i] * fe_0 + g_yyyyy_xxxxzz_1[i] * pa_z[i];

        g_yyyyyz_xxxyyy_0[i] = g_yyyyy_xxxyyy_1[i] * pa_z[i];

        g_yyyyyz_xxxyyz_0[i] = g_yyyyy_xxxyy_1[i] * fe_0 + g_yyyyy_xxxyyz_1[i] * pa_z[i];

        g_yyyyyz_xxxyzz_0[i] = 2.0 * g_yyyyy_xxxyz_1[i] * fe_0 + g_yyyyy_xxxyzz_1[i] * pa_z[i];

        g_yyyyyz_xxxzzz_0[i] = 3.0 * g_yyyyy_xxxzz_1[i] * fe_0 + g_yyyyy_xxxzzz_1[i] * pa_z[i];

        g_yyyyyz_xxyyyy_0[i] = g_yyyyy_xxyyyy_1[i] * pa_z[i];

        g_yyyyyz_xxyyyz_0[i] = g_yyyyy_xxyyy_1[i] * fe_0 + g_yyyyy_xxyyyz_1[i] * pa_z[i];

        g_yyyyyz_xxyyzz_0[i] = 2.0 * g_yyyyy_xxyyz_1[i] * fe_0 + g_yyyyy_xxyyzz_1[i] * pa_z[i];

        g_yyyyyz_xxyzzz_0[i] = 3.0 * g_yyyyy_xxyzz_1[i] * fe_0 + g_yyyyy_xxyzzz_1[i] * pa_z[i];

        g_yyyyyz_xxzzzz_0[i] = 4.0 * g_yyyyy_xxzzz_1[i] * fe_0 + g_yyyyy_xxzzzz_1[i] * pa_z[i];

        g_yyyyyz_xyyyyy_0[i] = g_yyyyy_xyyyyy_1[i] * pa_z[i];

        g_yyyyyz_xyyyyz_0[i] = g_yyyyy_xyyyy_1[i] * fe_0 + g_yyyyy_xyyyyz_1[i] * pa_z[i];

        g_yyyyyz_xyyyzz_0[i] = 2.0 * g_yyyyy_xyyyz_1[i] * fe_0 + g_yyyyy_xyyyzz_1[i] * pa_z[i];

        g_yyyyyz_xyyzzz_0[i] = 3.0 * g_yyyyy_xyyzz_1[i] * fe_0 + g_yyyyy_xyyzzz_1[i] * pa_z[i];

        g_yyyyyz_xyzzzz_0[i] = 4.0 * g_yyyyy_xyzzz_1[i] * fe_0 + g_yyyyy_xyzzzz_1[i] * pa_z[i];

        g_yyyyyz_xzzzzz_0[i] = 5.0 * g_yyyyy_xzzzz_1[i] * fe_0 + g_yyyyy_xzzzzz_1[i] * pa_z[i];

        g_yyyyyz_yyyyyy_0[i] = g_yyyyy_yyyyyy_1[i] * pa_z[i];

        g_yyyyyz_yyyyyz_0[i] = g_yyyyy_yyyyy_1[i] * fe_0 + g_yyyyy_yyyyyz_1[i] * pa_z[i];

        g_yyyyyz_yyyyzz_0[i] = 2.0 * g_yyyyy_yyyyz_1[i] * fe_0 + g_yyyyy_yyyyzz_1[i] * pa_z[i];

        g_yyyyyz_yyyzzz_0[i] = 3.0 * g_yyyyy_yyyzz_1[i] * fe_0 + g_yyyyy_yyyzzz_1[i] * pa_z[i];

        g_yyyyyz_yyzzzz_0[i] = 4.0 * g_yyyyy_yyzzz_1[i] * fe_0 + g_yyyyy_yyzzzz_1[i] * pa_z[i];

        g_yyyyyz_yzzzzz_0[i] = 5.0 * g_yyyyy_yzzzz_1[i] * fe_0 + g_yyyyy_yzzzzz_1[i] * pa_z[i];

        g_yyyyyz_zzzzzz_0[i] = 6.0 * g_yyyyy_zzzzz_1[i] * fe_0 + g_yyyyy_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 644-672 components of targeted buffer : II

    auto g_yyyyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 644);

    auto g_yyyyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 645);

    auto g_yyyyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 646);

    auto g_yyyyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 647);

    auto g_yyyyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 648);

    auto g_yyyyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 649);

    auto g_yyyyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 650);

    auto g_yyyyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 651);

    auto g_yyyyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 652);

    auto g_yyyyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 653);

    auto g_yyyyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 654);

    auto g_yyyyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 655);

    auto g_yyyyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 656);

    auto g_yyyyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 657);

    auto g_yyyyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 658);

    auto g_yyyyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 659);

    auto g_yyyyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 660);

    auto g_yyyyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 661);

    auto g_yyyyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 662);

    auto g_yyyyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 663);

    auto g_yyyyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 664);

    auto g_yyyyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 665);

    auto g_yyyyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 666);

    auto g_yyyyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 667);

    auto g_yyyyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 668);

    auto g_yyyyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 669);

    auto g_yyyyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 670);

    auto g_yyyyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 671);

    #pragma omp simd aligned(g_yyyy_xxxxxy_0, g_yyyy_xxxxxy_1, g_yyyy_xxxxyy_0, g_yyyy_xxxxyy_1, g_yyyy_xxxyyy_0, g_yyyy_xxxyyy_1, g_yyyy_xxyyyy_0, g_yyyy_xxyyyy_1, g_yyyy_xyyyyy_0, g_yyyy_xyyyyy_1, g_yyyy_yyyyyy_0, g_yyyy_yyyyyy_1, g_yyyyz_xxxxxy_1, g_yyyyz_xxxxyy_1, g_yyyyz_xxxyyy_1, g_yyyyz_xxyyyy_1, g_yyyyz_xyyyyy_1, g_yyyyz_yyyyyy_1, g_yyyyzz_xxxxxx_0, g_yyyyzz_xxxxxy_0, g_yyyyzz_xxxxxz_0, g_yyyyzz_xxxxyy_0, g_yyyyzz_xxxxyz_0, g_yyyyzz_xxxxzz_0, g_yyyyzz_xxxyyy_0, g_yyyyzz_xxxyyz_0, g_yyyyzz_xxxyzz_0, g_yyyyzz_xxxzzz_0, g_yyyyzz_xxyyyy_0, g_yyyyzz_xxyyyz_0, g_yyyyzz_xxyyzz_0, g_yyyyzz_xxyzzz_0, g_yyyyzz_xxzzzz_0, g_yyyyzz_xyyyyy_0, g_yyyyzz_xyyyyz_0, g_yyyyzz_xyyyzz_0, g_yyyyzz_xyyzzz_0, g_yyyyzz_xyzzzz_0, g_yyyyzz_xzzzzz_0, g_yyyyzz_yyyyyy_0, g_yyyyzz_yyyyyz_0, g_yyyyzz_yyyyzz_0, g_yyyyzz_yyyzzz_0, g_yyyyzz_yyzzzz_0, g_yyyyzz_yzzzzz_0, g_yyyyzz_zzzzzz_0, g_yyyzz_xxxxxx_1, g_yyyzz_xxxxxz_1, g_yyyzz_xxxxyz_1, g_yyyzz_xxxxz_1, g_yyyzz_xxxxzz_1, g_yyyzz_xxxyyz_1, g_yyyzz_xxxyz_1, g_yyyzz_xxxyzz_1, g_yyyzz_xxxzz_1, g_yyyzz_xxxzzz_1, g_yyyzz_xxyyyz_1, g_yyyzz_xxyyz_1, g_yyyzz_xxyyzz_1, g_yyyzz_xxyzz_1, g_yyyzz_xxyzzz_1, g_yyyzz_xxzzz_1, g_yyyzz_xxzzzz_1, g_yyyzz_xyyyyz_1, g_yyyzz_xyyyz_1, g_yyyzz_xyyyzz_1, g_yyyzz_xyyzz_1, g_yyyzz_xyyzzz_1, g_yyyzz_xyzzz_1, g_yyyzz_xyzzzz_1, g_yyyzz_xzzzz_1, g_yyyzz_xzzzzz_1, g_yyyzz_yyyyyz_1, g_yyyzz_yyyyz_1, g_yyyzz_yyyyzz_1, g_yyyzz_yyyzz_1, g_yyyzz_yyyzzz_1, g_yyyzz_yyzzz_1, g_yyyzz_yyzzzz_1, g_yyyzz_yzzzz_1, g_yyyzz_yzzzzz_1, g_yyyzz_zzzzz_1, g_yyyzz_zzzzzz_1, g_yyzz_xxxxxx_0, g_yyzz_xxxxxx_1, g_yyzz_xxxxxz_0, g_yyzz_xxxxxz_1, g_yyzz_xxxxyz_0, g_yyzz_xxxxyz_1, g_yyzz_xxxxzz_0, g_yyzz_xxxxzz_1, g_yyzz_xxxyyz_0, g_yyzz_xxxyyz_1, g_yyzz_xxxyzz_0, g_yyzz_xxxyzz_1, g_yyzz_xxxzzz_0, g_yyzz_xxxzzz_1, g_yyzz_xxyyyz_0, g_yyzz_xxyyyz_1, g_yyzz_xxyyzz_0, g_yyzz_xxyyzz_1, g_yyzz_xxyzzz_0, g_yyzz_xxyzzz_1, g_yyzz_xxzzzz_0, g_yyzz_xxzzzz_1, g_yyzz_xyyyyz_0, g_yyzz_xyyyyz_1, g_yyzz_xyyyzz_0, g_yyzz_xyyyzz_1, g_yyzz_xyyzzz_0, g_yyzz_xyyzzz_1, g_yyzz_xyzzzz_0, g_yyzz_xyzzzz_1, g_yyzz_xzzzzz_0, g_yyzz_xzzzzz_1, g_yyzz_yyyyyz_0, g_yyzz_yyyyyz_1, g_yyzz_yyyyzz_0, g_yyzz_yyyyzz_1, g_yyzz_yyyzzz_0, g_yyzz_yyyzzz_1, g_yyzz_yyzzzz_0, g_yyzz_yyzzzz_1, g_yyzz_yzzzzz_0, g_yyzz_yzzzzz_1, g_yyzz_zzzzzz_0, g_yyzz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyzz_xxxxxx_0[i] = 3.0 * g_yyzz_xxxxxx_0[i] * fbe_0 - 3.0 * g_yyzz_xxxxxx_1[i] * fz_be_0 + g_yyyzz_xxxxxx_1[i] * pa_y[i];

        g_yyyyzz_xxxxxy_0[i] = g_yyyy_xxxxxy_0[i] * fbe_0 - g_yyyy_xxxxxy_1[i] * fz_be_0 + g_yyyyz_xxxxxy_1[i] * pa_z[i];

        g_yyyyzz_xxxxxz_0[i] = 3.0 * g_yyzz_xxxxxz_0[i] * fbe_0 - 3.0 * g_yyzz_xxxxxz_1[i] * fz_be_0 + g_yyyzz_xxxxxz_1[i] * pa_y[i];

        g_yyyyzz_xxxxyy_0[i] = g_yyyy_xxxxyy_0[i] * fbe_0 - g_yyyy_xxxxyy_1[i] * fz_be_0 + g_yyyyz_xxxxyy_1[i] * pa_z[i];

        g_yyyyzz_xxxxyz_0[i] = 3.0 * g_yyzz_xxxxyz_0[i] * fbe_0 - 3.0 * g_yyzz_xxxxyz_1[i] * fz_be_0 + g_yyyzz_xxxxz_1[i] * fe_0 + g_yyyzz_xxxxyz_1[i] * pa_y[i];

        g_yyyyzz_xxxxzz_0[i] = 3.0 * g_yyzz_xxxxzz_0[i] * fbe_0 - 3.0 * g_yyzz_xxxxzz_1[i] * fz_be_0 + g_yyyzz_xxxxzz_1[i] * pa_y[i];

        g_yyyyzz_xxxyyy_0[i] = g_yyyy_xxxyyy_0[i] * fbe_0 - g_yyyy_xxxyyy_1[i] * fz_be_0 + g_yyyyz_xxxyyy_1[i] * pa_z[i];

        g_yyyyzz_xxxyyz_0[i] = 3.0 * g_yyzz_xxxyyz_0[i] * fbe_0 - 3.0 * g_yyzz_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyzz_xxxyz_1[i] * fe_0 + g_yyyzz_xxxyyz_1[i] * pa_y[i];

        g_yyyyzz_xxxyzz_0[i] = 3.0 * g_yyzz_xxxyzz_0[i] * fbe_0 - 3.0 * g_yyzz_xxxyzz_1[i] * fz_be_0 + g_yyyzz_xxxzz_1[i] * fe_0 + g_yyyzz_xxxyzz_1[i] * pa_y[i];

        g_yyyyzz_xxxzzz_0[i] = 3.0 * g_yyzz_xxxzzz_0[i] * fbe_0 - 3.0 * g_yyzz_xxxzzz_1[i] * fz_be_0 + g_yyyzz_xxxzzz_1[i] * pa_y[i];

        g_yyyyzz_xxyyyy_0[i] = g_yyyy_xxyyyy_0[i] * fbe_0 - g_yyyy_xxyyyy_1[i] * fz_be_0 + g_yyyyz_xxyyyy_1[i] * pa_z[i];

        g_yyyyzz_xxyyyz_0[i] = 3.0 * g_yyzz_xxyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyzz_xxyyz_1[i] * fe_0 + g_yyyzz_xxyyyz_1[i] * pa_y[i];

        g_yyyyzz_xxyyzz_0[i] = 3.0 * g_yyzz_xxyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_xxyzz_1[i] * fe_0 + g_yyyzz_xxyyzz_1[i] * pa_y[i];

        g_yyyyzz_xxyzzz_0[i] = 3.0 * g_yyzz_xxyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_xxyzzz_1[i] * fz_be_0 + g_yyyzz_xxzzz_1[i] * fe_0 + g_yyyzz_xxyzzz_1[i] * pa_y[i];

        g_yyyyzz_xxzzzz_0[i] = 3.0 * g_yyzz_xxzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_xxzzzz_1[i] * fz_be_0 + g_yyyzz_xxzzzz_1[i] * pa_y[i];

        g_yyyyzz_xyyyyy_0[i] = g_yyyy_xyyyyy_0[i] * fbe_0 - g_yyyy_xyyyyy_1[i] * fz_be_0 + g_yyyyz_xyyyyy_1[i] * pa_z[i];

        g_yyyyzz_xyyyyz_0[i] = 3.0 * g_yyzz_xyyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyzz_xyyyz_1[i] * fe_0 + g_yyyzz_xyyyyz_1[i] * pa_y[i];

        g_yyyyzz_xyyyzz_0[i] = 3.0 * g_yyzz_xyyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_xyyzz_1[i] * fe_0 + g_yyyzz_xyyyzz_1[i] * pa_y[i];

        g_yyyyzz_xyyzzz_0[i] = 3.0 * g_yyzz_xyyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_xyzzz_1[i] * fe_0 + g_yyyzz_xyyzzz_1[i] * pa_y[i];

        g_yyyyzz_xyzzzz_0[i] = 3.0 * g_yyzz_xyzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_xyzzzz_1[i] * fz_be_0 + g_yyyzz_xzzzz_1[i] * fe_0 + g_yyyzz_xyzzzz_1[i] * pa_y[i];

        g_yyyyzz_xzzzzz_0[i] = 3.0 * g_yyzz_xzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_xzzzzz_1[i] * fz_be_0 + g_yyyzz_xzzzzz_1[i] * pa_y[i];

        g_yyyyzz_yyyyyy_0[i] = g_yyyy_yyyyyy_0[i] * fbe_0 - g_yyyy_yyyyyy_1[i] * fz_be_0 + g_yyyyz_yyyyyy_1[i] * pa_z[i];

        g_yyyyzz_yyyyyz_0[i] = 3.0 * g_yyzz_yyyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyzz_yyyyz_1[i] * fe_0 + g_yyyzz_yyyyyz_1[i] * pa_y[i];

        g_yyyyzz_yyyyzz_0[i] = 3.0 * g_yyzz_yyyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyzz_yyyzz_1[i] * fe_0 + g_yyyzz_yyyyzz_1[i] * pa_y[i];

        g_yyyyzz_yyyzzz_0[i] = 3.0 * g_yyzz_yyyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_yyzzz_1[i] * fe_0 + g_yyyzz_yyyzzz_1[i] * pa_y[i];

        g_yyyyzz_yyzzzz_0[i] = 3.0 * g_yyzz_yyzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_yzzzz_1[i] * fe_0 + g_yyyzz_yyzzzz_1[i] * pa_y[i];

        g_yyyyzz_yzzzzz_0[i] = 3.0 * g_yyzz_yzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_yzzzzz_1[i] * fz_be_0 + g_yyyzz_zzzzz_1[i] * fe_0 + g_yyyzz_yzzzzz_1[i] * pa_y[i];

        g_yyyyzz_zzzzzz_0[i] = 3.0 * g_yyzz_zzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_zzzzzz_1[i] * fz_be_0 + g_yyyzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 672-700 components of targeted buffer : II

    auto g_yyyzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 672);

    auto g_yyyzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 673);

    auto g_yyyzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 674);

    auto g_yyyzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 675);

    auto g_yyyzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 676);

    auto g_yyyzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 677);

    auto g_yyyzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 678);

    auto g_yyyzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 679);

    auto g_yyyzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 680);

    auto g_yyyzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 681);

    auto g_yyyzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 682);

    auto g_yyyzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 683);

    auto g_yyyzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 684);

    auto g_yyyzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 685);

    auto g_yyyzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 686);

    auto g_yyyzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 687);

    auto g_yyyzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 688);

    auto g_yyyzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 689);

    auto g_yyyzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 690);

    auto g_yyyzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 691);

    auto g_yyyzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 692);

    auto g_yyyzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 693);

    auto g_yyyzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 694);

    auto g_yyyzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 695);

    auto g_yyyzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 696);

    auto g_yyyzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 697);

    auto g_yyyzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 698);

    auto g_yyyzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 699);

    #pragma omp simd aligned(g_yyyz_xxxxxy_0, g_yyyz_xxxxxy_1, g_yyyz_xxxxyy_0, g_yyyz_xxxxyy_1, g_yyyz_xxxyyy_0, g_yyyz_xxxyyy_1, g_yyyz_xxyyyy_0, g_yyyz_xxyyyy_1, g_yyyz_xyyyyy_0, g_yyyz_xyyyyy_1, g_yyyz_yyyyyy_0, g_yyyz_yyyyyy_1, g_yyyzz_xxxxxy_1, g_yyyzz_xxxxyy_1, g_yyyzz_xxxyyy_1, g_yyyzz_xxyyyy_1, g_yyyzz_xyyyyy_1, g_yyyzz_yyyyyy_1, g_yyyzzz_xxxxxx_0, g_yyyzzz_xxxxxy_0, g_yyyzzz_xxxxxz_0, g_yyyzzz_xxxxyy_0, g_yyyzzz_xxxxyz_0, g_yyyzzz_xxxxzz_0, g_yyyzzz_xxxyyy_0, g_yyyzzz_xxxyyz_0, g_yyyzzz_xxxyzz_0, g_yyyzzz_xxxzzz_0, g_yyyzzz_xxyyyy_0, g_yyyzzz_xxyyyz_0, g_yyyzzz_xxyyzz_0, g_yyyzzz_xxyzzz_0, g_yyyzzz_xxzzzz_0, g_yyyzzz_xyyyyy_0, g_yyyzzz_xyyyyz_0, g_yyyzzz_xyyyzz_0, g_yyyzzz_xyyzzz_0, g_yyyzzz_xyzzzz_0, g_yyyzzz_xzzzzz_0, g_yyyzzz_yyyyyy_0, g_yyyzzz_yyyyyz_0, g_yyyzzz_yyyyzz_0, g_yyyzzz_yyyzzz_0, g_yyyzzz_yyzzzz_0, g_yyyzzz_yzzzzz_0, g_yyyzzz_zzzzzz_0, g_yyzzz_xxxxxx_1, g_yyzzz_xxxxxz_1, g_yyzzz_xxxxyz_1, g_yyzzz_xxxxz_1, g_yyzzz_xxxxzz_1, g_yyzzz_xxxyyz_1, g_yyzzz_xxxyz_1, g_yyzzz_xxxyzz_1, g_yyzzz_xxxzz_1, g_yyzzz_xxxzzz_1, g_yyzzz_xxyyyz_1, g_yyzzz_xxyyz_1, g_yyzzz_xxyyzz_1, g_yyzzz_xxyzz_1, g_yyzzz_xxyzzz_1, g_yyzzz_xxzzz_1, g_yyzzz_xxzzzz_1, g_yyzzz_xyyyyz_1, g_yyzzz_xyyyz_1, g_yyzzz_xyyyzz_1, g_yyzzz_xyyzz_1, g_yyzzz_xyyzzz_1, g_yyzzz_xyzzz_1, g_yyzzz_xyzzzz_1, g_yyzzz_xzzzz_1, g_yyzzz_xzzzzz_1, g_yyzzz_yyyyyz_1, g_yyzzz_yyyyz_1, g_yyzzz_yyyyzz_1, g_yyzzz_yyyzz_1, g_yyzzz_yyyzzz_1, g_yyzzz_yyzzz_1, g_yyzzz_yyzzzz_1, g_yyzzz_yzzzz_1, g_yyzzz_yzzzzz_1, g_yyzzz_zzzzz_1, g_yyzzz_zzzzzz_1, g_yzzz_xxxxxx_0, g_yzzz_xxxxxx_1, g_yzzz_xxxxxz_0, g_yzzz_xxxxxz_1, g_yzzz_xxxxyz_0, g_yzzz_xxxxyz_1, g_yzzz_xxxxzz_0, g_yzzz_xxxxzz_1, g_yzzz_xxxyyz_0, g_yzzz_xxxyyz_1, g_yzzz_xxxyzz_0, g_yzzz_xxxyzz_1, g_yzzz_xxxzzz_0, g_yzzz_xxxzzz_1, g_yzzz_xxyyyz_0, g_yzzz_xxyyyz_1, g_yzzz_xxyyzz_0, g_yzzz_xxyyzz_1, g_yzzz_xxyzzz_0, g_yzzz_xxyzzz_1, g_yzzz_xxzzzz_0, g_yzzz_xxzzzz_1, g_yzzz_xyyyyz_0, g_yzzz_xyyyyz_1, g_yzzz_xyyyzz_0, g_yzzz_xyyyzz_1, g_yzzz_xyyzzz_0, g_yzzz_xyyzzz_1, g_yzzz_xyzzzz_0, g_yzzz_xyzzzz_1, g_yzzz_xzzzzz_0, g_yzzz_xzzzzz_1, g_yzzz_yyyyyz_0, g_yzzz_yyyyyz_1, g_yzzz_yyyyzz_0, g_yzzz_yyyyzz_1, g_yzzz_yyyzzz_0, g_yzzz_yyyzzz_1, g_yzzz_yyzzzz_0, g_yzzz_yyzzzz_1, g_yzzz_yzzzzz_0, g_yzzz_yzzzzz_1, g_yzzz_zzzzzz_0, g_yzzz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzzz_xxxxxx_0[i] = 2.0 * g_yzzz_xxxxxx_0[i] * fbe_0 - 2.0 * g_yzzz_xxxxxx_1[i] * fz_be_0 + g_yyzzz_xxxxxx_1[i] * pa_y[i];

        g_yyyzzz_xxxxxy_0[i] = 2.0 * g_yyyz_xxxxxy_0[i] * fbe_0 - 2.0 * g_yyyz_xxxxxy_1[i] * fz_be_0 + g_yyyzz_xxxxxy_1[i] * pa_z[i];

        g_yyyzzz_xxxxxz_0[i] = 2.0 * g_yzzz_xxxxxz_0[i] * fbe_0 - 2.0 * g_yzzz_xxxxxz_1[i] * fz_be_0 + g_yyzzz_xxxxxz_1[i] * pa_y[i];

        g_yyyzzz_xxxxyy_0[i] = 2.0 * g_yyyz_xxxxyy_0[i] * fbe_0 - 2.0 * g_yyyz_xxxxyy_1[i] * fz_be_0 + g_yyyzz_xxxxyy_1[i] * pa_z[i];

        g_yyyzzz_xxxxyz_0[i] = 2.0 * g_yzzz_xxxxyz_0[i] * fbe_0 - 2.0 * g_yzzz_xxxxyz_1[i] * fz_be_0 + g_yyzzz_xxxxz_1[i] * fe_0 + g_yyzzz_xxxxyz_1[i] * pa_y[i];

        g_yyyzzz_xxxxzz_0[i] = 2.0 * g_yzzz_xxxxzz_0[i] * fbe_0 - 2.0 * g_yzzz_xxxxzz_1[i] * fz_be_0 + g_yyzzz_xxxxzz_1[i] * pa_y[i];

        g_yyyzzz_xxxyyy_0[i] = 2.0 * g_yyyz_xxxyyy_0[i] * fbe_0 - 2.0 * g_yyyz_xxxyyy_1[i] * fz_be_0 + g_yyyzz_xxxyyy_1[i] * pa_z[i];

        g_yyyzzz_xxxyyz_0[i] = 2.0 * g_yzzz_xxxyyz_0[i] * fbe_0 - 2.0 * g_yzzz_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyzzz_xxxyz_1[i] * fe_0 + g_yyzzz_xxxyyz_1[i] * pa_y[i];

        g_yyyzzz_xxxyzz_0[i] = 2.0 * g_yzzz_xxxyzz_0[i] * fbe_0 - 2.0 * g_yzzz_xxxyzz_1[i] * fz_be_0 + g_yyzzz_xxxzz_1[i] * fe_0 + g_yyzzz_xxxyzz_1[i] * pa_y[i];

        g_yyyzzz_xxxzzz_0[i] = 2.0 * g_yzzz_xxxzzz_0[i] * fbe_0 - 2.0 * g_yzzz_xxxzzz_1[i] * fz_be_0 + g_yyzzz_xxxzzz_1[i] * pa_y[i];

        g_yyyzzz_xxyyyy_0[i] = 2.0 * g_yyyz_xxyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_xxyyyy_1[i] * fz_be_0 + g_yyyzz_xxyyyy_1[i] * pa_z[i];

        g_yyyzzz_xxyyyz_0[i] = 2.0 * g_yzzz_xxyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyzzz_xxyyz_1[i] * fe_0 + g_yyzzz_xxyyyz_1[i] * pa_y[i];

        g_yyyzzz_xxyyzz_0[i] = 2.0 * g_yzzz_xxyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_xxyzz_1[i] * fe_0 + g_yyzzz_xxyyzz_1[i] * pa_y[i];

        g_yyyzzz_xxyzzz_0[i] = 2.0 * g_yzzz_xxyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_xxyzzz_1[i] * fz_be_0 + g_yyzzz_xxzzz_1[i] * fe_0 + g_yyzzz_xxyzzz_1[i] * pa_y[i];

        g_yyyzzz_xxzzzz_0[i] = 2.0 * g_yzzz_xxzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_xxzzzz_1[i] * fz_be_0 + g_yyzzz_xxzzzz_1[i] * pa_y[i];

        g_yyyzzz_xyyyyy_0[i] = 2.0 * g_yyyz_xyyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_xyyyyy_1[i] * fz_be_0 + g_yyyzz_xyyyyy_1[i] * pa_z[i];

        g_yyyzzz_xyyyyz_0[i] = 2.0 * g_yzzz_xyyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyzzz_xyyyz_1[i] * fe_0 + g_yyzzz_xyyyyz_1[i] * pa_y[i];

        g_yyyzzz_xyyyzz_0[i] = 2.0 * g_yzzz_xyyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_xyyzz_1[i] * fe_0 + g_yyzzz_xyyyzz_1[i] * pa_y[i];

        g_yyyzzz_xyyzzz_0[i] = 2.0 * g_yzzz_xyyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_xyzzz_1[i] * fe_0 + g_yyzzz_xyyzzz_1[i] * pa_y[i];

        g_yyyzzz_xyzzzz_0[i] = 2.0 * g_yzzz_xyzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_xyzzzz_1[i] * fz_be_0 + g_yyzzz_xzzzz_1[i] * fe_0 + g_yyzzz_xyzzzz_1[i] * pa_y[i];

        g_yyyzzz_xzzzzz_0[i] = 2.0 * g_yzzz_xzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_xzzzzz_1[i] * fz_be_0 + g_yyzzz_xzzzzz_1[i] * pa_y[i];

        g_yyyzzz_yyyyyy_0[i] = 2.0 * g_yyyz_yyyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_yyyyyy_1[i] * fz_be_0 + g_yyyzz_yyyyyy_1[i] * pa_z[i];

        g_yyyzzz_yyyyyz_0[i] = 2.0 * g_yzzz_yyyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyzzz_yyyyz_1[i] * fe_0 + g_yyzzz_yyyyyz_1[i] * pa_y[i];

        g_yyyzzz_yyyyzz_0[i] = 2.0 * g_yzzz_yyyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyzzz_yyyzz_1[i] * fe_0 + g_yyzzz_yyyyzz_1[i] * pa_y[i];

        g_yyyzzz_yyyzzz_0[i] = 2.0 * g_yzzz_yyyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_yyzzz_1[i] * fe_0 + g_yyzzz_yyyzzz_1[i] * pa_y[i];

        g_yyyzzz_yyzzzz_0[i] = 2.0 * g_yzzz_yyzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_yzzzz_1[i] * fe_0 + g_yyzzz_yyzzzz_1[i] * pa_y[i];

        g_yyyzzz_yzzzzz_0[i] = 2.0 * g_yzzz_yzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_yzzzzz_1[i] * fz_be_0 + g_yyzzz_zzzzz_1[i] * fe_0 + g_yyzzz_yzzzzz_1[i] * pa_y[i];

        g_yyyzzz_zzzzzz_0[i] = 2.0 * g_yzzz_zzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_zzzzzz_1[i] * fz_be_0 + g_yyzzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 700-728 components of targeted buffer : II

    auto g_yyzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 700);

    auto g_yyzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 701);

    auto g_yyzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 702);

    auto g_yyzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 703);

    auto g_yyzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 704);

    auto g_yyzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 705);

    auto g_yyzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 706);

    auto g_yyzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 707);

    auto g_yyzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 708);

    auto g_yyzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 709);

    auto g_yyzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 710);

    auto g_yyzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 711);

    auto g_yyzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 712);

    auto g_yyzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 713);

    auto g_yyzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 714);

    auto g_yyzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 715);

    auto g_yyzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 716);

    auto g_yyzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 717);

    auto g_yyzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 718);

    auto g_yyzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 719);

    auto g_yyzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 720);

    auto g_yyzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 721);

    auto g_yyzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 722);

    auto g_yyzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 723);

    auto g_yyzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 724);

    auto g_yyzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 725);

    auto g_yyzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 726);

    auto g_yyzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 727);

    #pragma omp simd aligned(g_yyzz_xxxxxy_0, g_yyzz_xxxxxy_1, g_yyzz_xxxxyy_0, g_yyzz_xxxxyy_1, g_yyzz_xxxyyy_0, g_yyzz_xxxyyy_1, g_yyzz_xxyyyy_0, g_yyzz_xxyyyy_1, g_yyzz_xyyyyy_0, g_yyzz_xyyyyy_1, g_yyzz_yyyyyy_0, g_yyzz_yyyyyy_1, g_yyzzz_xxxxxy_1, g_yyzzz_xxxxyy_1, g_yyzzz_xxxyyy_1, g_yyzzz_xxyyyy_1, g_yyzzz_xyyyyy_1, g_yyzzz_yyyyyy_1, g_yyzzzz_xxxxxx_0, g_yyzzzz_xxxxxy_0, g_yyzzzz_xxxxxz_0, g_yyzzzz_xxxxyy_0, g_yyzzzz_xxxxyz_0, g_yyzzzz_xxxxzz_0, g_yyzzzz_xxxyyy_0, g_yyzzzz_xxxyyz_0, g_yyzzzz_xxxyzz_0, g_yyzzzz_xxxzzz_0, g_yyzzzz_xxyyyy_0, g_yyzzzz_xxyyyz_0, g_yyzzzz_xxyyzz_0, g_yyzzzz_xxyzzz_0, g_yyzzzz_xxzzzz_0, g_yyzzzz_xyyyyy_0, g_yyzzzz_xyyyyz_0, g_yyzzzz_xyyyzz_0, g_yyzzzz_xyyzzz_0, g_yyzzzz_xyzzzz_0, g_yyzzzz_xzzzzz_0, g_yyzzzz_yyyyyy_0, g_yyzzzz_yyyyyz_0, g_yyzzzz_yyyyzz_0, g_yyzzzz_yyyzzz_0, g_yyzzzz_yyzzzz_0, g_yyzzzz_yzzzzz_0, g_yyzzzz_zzzzzz_0, g_yzzzz_xxxxxx_1, g_yzzzz_xxxxxz_1, g_yzzzz_xxxxyz_1, g_yzzzz_xxxxz_1, g_yzzzz_xxxxzz_1, g_yzzzz_xxxyyz_1, g_yzzzz_xxxyz_1, g_yzzzz_xxxyzz_1, g_yzzzz_xxxzz_1, g_yzzzz_xxxzzz_1, g_yzzzz_xxyyyz_1, g_yzzzz_xxyyz_1, g_yzzzz_xxyyzz_1, g_yzzzz_xxyzz_1, g_yzzzz_xxyzzz_1, g_yzzzz_xxzzz_1, g_yzzzz_xxzzzz_1, g_yzzzz_xyyyyz_1, g_yzzzz_xyyyz_1, g_yzzzz_xyyyzz_1, g_yzzzz_xyyzz_1, g_yzzzz_xyyzzz_1, g_yzzzz_xyzzz_1, g_yzzzz_xyzzzz_1, g_yzzzz_xzzzz_1, g_yzzzz_xzzzzz_1, g_yzzzz_yyyyyz_1, g_yzzzz_yyyyz_1, g_yzzzz_yyyyzz_1, g_yzzzz_yyyzz_1, g_yzzzz_yyyzzz_1, g_yzzzz_yyzzz_1, g_yzzzz_yyzzzz_1, g_yzzzz_yzzzz_1, g_yzzzz_yzzzzz_1, g_yzzzz_zzzzz_1, g_yzzzz_zzzzzz_1, g_zzzz_xxxxxx_0, g_zzzz_xxxxxx_1, g_zzzz_xxxxxz_0, g_zzzz_xxxxxz_1, g_zzzz_xxxxyz_0, g_zzzz_xxxxyz_1, g_zzzz_xxxxzz_0, g_zzzz_xxxxzz_1, g_zzzz_xxxyyz_0, g_zzzz_xxxyyz_1, g_zzzz_xxxyzz_0, g_zzzz_xxxyzz_1, g_zzzz_xxxzzz_0, g_zzzz_xxxzzz_1, g_zzzz_xxyyyz_0, g_zzzz_xxyyyz_1, g_zzzz_xxyyzz_0, g_zzzz_xxyyzz_1, g_zzzz_xxyzzz_0, g_zzzz_xxyzzz_1, g_zzzz_xxzzzz_0, g_zzzz_xxzzzz_1, g_zzzz_xyyyyz_0, g_zzzz_xyyyyz_1, g_zzzz_xyyyzz_0, g_zzzz_xyyyzz_1, g_zzzz_xyyzzz_0, g_zzzz_xyyzzz_1, g_zzzz_xyzzzz_0, g_zzzz_xyzzzz_1, g_zzzz_xzzzzz_0, g_zzzz_xzzzzz_1, g_zzzz_yyyyyz_0, g_zzzz_yyyyyz_1, g_zzzz_yyyyzz_0, g_zzzz_yyyyzz_1, g_zzzz_yyyzzz_0, g_zzzz_yyyzzz_1, g_zzzz_yyzzzz_0, g_zzzz_yyzzzz_1, g_zzzz_yzzzzz_0, g_zzzz_yzzzzz_1, g_zzzz_zzzzzz_0, g_zzzz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzzz_xxxxxx_0[i] = g_zzzz_xxxxxx_0[i] * fbe_0 - g_zzzz_xxxxxx_1[i] * fz_be_0 + g_yzzzz_xxxxxx_1[i] * pa_y[i];

        g_yyzzzz_xxxxxy_0[i] = 3.0 * g_yyzz_xxxxxy_0[i] * fbe_0 - 3.0 * g_yyzz_xxxxxy_1[i] * fz_be_0 + g_yyzzz_xxxxxy_1[i] * pa_z[i];

        g_yyzzzz_xxxxxz_0[i] = g_zzzz_xxxxxz_0[i] * fbe_0 - g_zzzz_xxxxxz_1[i] * fz_be_0 + g_yzzzz_xxxxxz_1[i] * pa_y[i];

        g_yyzzzz_xxxxyy_0[i] = 3.0 * g_yyzz_xxxxyy_0[i] * fbe_0 - 3.0 * g_yyzz_xxxxyy_1[i] * fz_be_0 + g_yyzzz_xxxxyy_1[i] * pa_z[i];

        g_yyzzzz_xxxxyz_0[i] = g_zzzz_xxxxyz_0[i] * fbe_0 - g_zzzz_xxxxyz_1[i] * fz_be_0 + g_yzzzz_xxxxz_1[i] * fe_0 + g_yzzzz_xxxxyz_1[i] * pa_y[i];

        g_yyzzzz_xxxxzz_0[i] = g_zzzz_xxxxzz_0[i] * fbe_0 - g_zzzz_xxxxzz_1[i] * fz_be_0 + g_yzzzz_xxxxzz_1[i] * pa_y[i];

        g_yyzzzz_xxxyyy_0[i] = 3.0 * g_yyzz_xxxyyy_0[i] * fbe_0 - 3.0 * g_yyzz_xxxyyy_1[i] * fz_be_0 + g_yyzzz_xxxyyy_1[i] * pa_z[i];

        g_yyzzzz_xxxyyz_0[i] = g_zzzz_xxxyyz_0[i] * fbe_0 - g_zzzz_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yzzzz_xxxyz_1[i] * fe_0 + g_yzzzz_xxxyyz_1[i] * pa_y[i];

        g_yyzzzz_xxxyzz_0[i] = g_zzzz_xxxyzz_0[i] * fbe_0 - g_zzzz_xxxyzz_1[i] * fz_be_0 + g_yzzzz_xxxzz_1[i] * fe_0 + g_yzzzz_xxxyzz_1[i] * pa_y[i];

        g_yyzzzz_xxxzzz_0[i] = g_zzzz_xxxzzz_0[i] * fbe_0 - g_zzzz_xxxzzz_1[i] * fz_be_0 + g_yzzzz_xxxzzz_1[i] * pa_y[i];

        g_yyzzzz_xxyyyy_0[i] = 3.0 * g_yyzz_xxyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_xxyyyy_1[i] * fz_be_0 + g_yyzzz_xxyyyy_1[i] * pa_z[i];

        g_yyzzzz_xxyyyz_0[i] = g_zzzz_xxyyyz_0[i] * fbe_0 - g_zzzz_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yzzzz_xxyyz_1[i] * fe_0 + g_yzzzz_xxyyyz_1[i] * pa_y[i];

        g_yyzzzz_xxyyzz_0[i] = g_zzzz_xxyyzz_0[i] * fbe_0 - g_zzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_xxyzz_1[i] * fe_0 + g_yzzzz_xxyyzz_1[i] * pa_y[i];

        g_yyzzzz_xxyzzz_0[i] = g_zzzz_xxyzzz_0[i] * fbe_0 - g_zzzz_xxyzzz_1[i] * fz_be_0 + g_yzzzz_xxzzz_1[i] * fe_0 + g_yzzzz_xxyzzz_1[i] * pa_y[i];

        g_yyzzzz_xxzzzz_0[i] = g_zzzz_xxzzzz_0[i] * fbe_0 - g_zzzz_xxzzzz_1[i] * fz_be_0 + g_yzzzz_xxzzzz_1[i] * pa_y[i];

        g_yyzzzz_xyyyyy_0[i] = 3.0 * g_yyzz_xyyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_xyyyyy_1[i] * fz_be_0 + g_yyzzz_xyyyyy_1[i] * pa_z[i];

        g_yyzzzz_xyyyyz_0[i] = g_zzzz_xyyyyz_0[i] * fbe_0 - g_zzzz_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yzzzz_xyyyz_1[i] * fe_0 + g_yzzzz_xyyyyz_1[i] * pa_y[i];

        g_yyzzzz_xyyyzz_0[i] = g_zzzz_xyyyzz_0[i] * fbe_0 - g_zzzz_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_xyyzz_1[i] * fe_0 + g_yzzzz_xyyyzz_1[i] * pa_y[i];

        g_yyzzzz_xyyzzz_0[i] = g_zzzz_xyyzzz_0[i] * fbe_0 - g_zzzz_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_xyzzz_1[i] * fe_0 + g_yzzzz_xyyzzz_1[i] * pa_y[i];

        g_yyzzzz_xyzzzz_0[i] = g_zzzz_xyzzzz_0[i] * fbe_0 - g_zzzz_xyzzzz_1[i] * fz_be_0 + g_yzzzz_xzzzz_1[i] * fe_0 + g_yzzzz_xyzzzz_1[i] * pa_y[i];

        g_yyzzzz_xzzzzz_0[i] = g_zzzz_xzzzzz_0[i] * fbe_0 - g_zzzz_xzzzzz_1[i] * fz_be_0 + g_yzzzz_xzzzzz_1[i] * pa_y[i];

        g_yyzzzz_yyyyyy_0[i] = 3.0 * g_yyzz_yyyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_yyyyyy_1[i] * fz_be_0 + g_yyzzz_yyyyyy_1[i] * pa_z[i];

        g_yyzzzz_yyyyyz_0[i] = g_zzzz_yyyyyz_0[i] * fbe_0 - g_zzzz_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yzzzz_yyyyz_1[i] * fe_0 + g_yzzzz_yyyyyz_1[i] * pa_y[i];

        g_yyzzzz_yyyyzz_0[i] = g_zzzz_yyyyzz_0[i] * fbe_0 - g_zzzz_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yzzzz_yyyzz_1[i] * fe_0 + g_yzzzz_yyyyzz_1[i] * pa_y[i];

        g_yyzzzz_yyyzzz_0[i] = g_zzzz_yyyzzz_0[i] * fbe_0 - g_zzzz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_yyzzz_1[i] * fe_0 + g_yzzzz_yyyzzz_1[i] * pa_y[i];

        g_yyzzzz_yyzzzz_0[i] = g_zzzz_yyzzzz_0[i] * fbe_0 - g_zzzz_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_yzzzz_1[i] * fe_0 + g_yzzzz_yyzzzz_1[i] * pa_y[i];

        g_yyzzzz_yzzzzz_0[i] = g_zzzz_yzzzzz_0[i] * fbe_0 - g_zzzz_yzzzzz_1[i] * fz_be_0 + g_yzzzz_zzzzz_1[i] * fe_0 + g_yzzzz_yzzzzz_1[i] * pa_y[i];

        g_yyzzzz_zzzzzz_0[i] = g_zzzz_zzzzzz_0[i] * fbe_0 - g_zzzz_zzzzzz_1[i] * fz_be_0 + g_yzzzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 728-756 components of targeted buffer : II

    auto g_yzzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 728);

    auto g_yzzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 729);

    auto g_yzzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 730);

    auto g_yzzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 731);

    auto g_yzzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 732);

    auto g_yzzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 733);

    auto g_yzzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 734);

    auto g_yzzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 735);

    auto g_yzzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 736);

    auto g_yzzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 737);

    auto g_yzzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 738);

    auto g_yzzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 739);

    auto g_yzzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 740);

    auto g_yzzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 741);

    auto g_yzzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 742);

    auto g_yzzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 743);

    auto g_yzzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 744);

    auto g_yzzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 745);

    auto g_yzzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 746);

    auto g_yzzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 747);

    auto g_yzzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 748);

    auto g_yzzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 749);

    auto g_yzzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 750);

    auto g_yzzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 751);

    auto g_yzzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 752);

    auto g_yzzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 753);

    auto g_yzzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 754);

    auto g_yzzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 755);

    #pragma omp simd aligned(g_yzzzzz_xxxxxx_0, g_yzzzzz_xxxxxy_0, g_yzzzzz_xxxxxz_0, g_yzzzzz_xxxxyy_0, g_yzzzzz_xxxxyz_0, g_yzzzzz_xxxxzz_0, g_yzzzzz_xxxyyy_0, g_yzzzzz_xxxyyz_0, g_yzzzzz_xxxyzz_0, g_yzzzzz_xxxzzz_0, g_yzzzzz_xxyyyy_0, g_yzzzzz_xxyyyz_0, g_yzzzzz_xxyyzz_0, g_yzzzzz_xxyzzz_0, g_yzzzzz_xxzzzz_0, g_yzzzzz_xyyyyy_0, g_yzzzzz_xyyyyz_0, g_yzzzzz_xyyyzz_0, g_yzzzzz_xyyzzz_0, g_yzzzzz_xyzzzz_0, g_yzzzzz_xzzzzz_0, g_yzzzzz_yyyyyy_0, g_yzzzzz_yyyyyz_0, g_yzzzzz_yyyyzz_0, g_yzzzzz_yyyzzz_0, g_yzzzzz_yyzzzz_0, g_yzzzzz_yzzzzz_0, g_yzzzzz_zzzzzz_0, g_zzzzz_xxxxx_1, g_zzzzz_xxxxxx_1, g_zzzzz_xxxxxy_1, g_zzzzz_xxxxxz_1, g_zzzzz_xxxxy_1, g_zzzzz_xxxxyy_1, g_zzzzz_xxxxyz_1, g_zzzzz_xxxxz_1, g_zzzzz_xxxxzz_1, g_zzzzz_xxxyy_1, g_zzzzz_xxxyyy_1, g_zzzzz_xxxyyz_1, g_zzzzz_xxxyz_1, g_zzzzz_xxxyzz_1, g_zzzzz_xxxzz_1, g_zzzzz_xxxzzz_1, g_zzzzz_xxyyy_1, g_zzzzz_xxyyyy_1, g_zzzzz_xxyyyz_1, g_zzzzz_xxyyz_1, g_zzzzz_xxyyzz_1, g_zzzzz_xxyzz_1, g_zzzzz_xxyzzz_1, g_zzzzz_xxzzz_1, g_zzzzz_xxzzzz_1, g_zzzzz_xyyyy_1, g_zzzzz_xyyyyy_1, g_zzzzz_xyyyyz_1, g_zzzzz_xyyyz_1, g_zzzzz_xyyyzz_1, g_zzzzz_xyyzz_1, g_zzzzz_xyyzzz_1, g_zzzzz_xyzzz_1, g_zzzzz_xyzzzz_1, g_zzzzz_xzzzz_1, g_zzzzz_xzzzzz_1, g_zzzzz_yyyyy_1, g_zzzzz_yyyyyy_1, g_zzzzz_yyyyyz_1, g_zzzzz_yyyyz_1, g_zzzzz_yyyyzz_1, g_zzzzz_yyyzz_1, g_zzzzz_yyyzzz_1, g_zzzzz_yyzzz_1, g_zzzzz_yyzzzz_1, g_zzzzz_yzzzz_1, g_zzzzz_yzzzzz_1, g_zzzzz_zzzzz_1, g_zzzzz_zzzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzzz_xxxxxx_0[i] = g_zzzzz_xxxxxx_1[i] * pa_y[i];

        g_yzzzzz_xxxxxy_0[i] = g_zzzzz_xxxxx_1[i] * fe_0 + g_zzzzz_xxxxxy_1[i] * pa_y[i];

        g_yzzzzz_xxxxxz_0[i] = g_zzzzz_xxxxxz_1[i] * pa_y[i];

        g_yzzzzz_xxxxyy_0[i] = 2.0 * g_zzzzz_xxxxy_1[i] * fe_0 + g_zzzzz_xxxxyy_1[i] * pa_y[i];

        g_yzzzzz_xxxxyz_0[i] = g_zzzzz_xxxxz_1[i] * fe_0 + g_zzzzz_xxxxyz_1[i] * pa_y[i];

        g_yzzzzz_xxxxzz_0[i] = g_zzzzz_xxxxzz_1[i] * pa_y[i];

        g_yzzzzz_xxxyyy_0[i] = 3.0 * g_zzzzz_xxxyy_1[i] * fe_0 + g_zzzzz_xxxyyy_1[i] * pa_y[i];

        g_yzzzzz_xxxyyz_0[i] = 2.0 * g_zzzzz_xxxyz_1[i] * fe_0 + g_zzzzz_xxxyyz_1[i] * pa_y[i];

        g_yzzzzz_xxxyzz_0[i] = g_zzzzz_xxxzz_1[i] * fe_0 + g_zzzzz_xxxyzz_1[i] * pa_y[i];

        g_yzzzzz_xxxzzz_0[i] = g_zzzzz_xxxzzz_1[i] * pa_y[i];

        g_yzzzzz_xxyyyy_0[i] = 4.0 * g_zzzzz_xxyyy_1[i] * fe_0 + g_zzzzz_xxyyyy_1[i] * pa_y[i];

        g_yzzzzz_xxyyyz_0[i] = 3.0 * g_zzzzz_xxyyz_1[i] * fe_0 + g_zzzzz_xxyyyz_1[i] * pa_y[i];

        g_yzzzzz_xxyyzz_0[i] = 2.0 * g_zzzzz_xxyzz_1[i] * fe_0 + g_zzzzz_xxyyzz_1[i] * pa_y[i];

        g_yzzzzz_xxyzzz_0[i] = g_zzzzz_xxzzz_1[i] * fe_0 + g_zzzzz_xxyzzz_1[i] * pa_y[i];

        g_yzzzzz_xxzzzz_0[i] = g_zzzzz_xxzzzz_1[i] * pa_y[i];

        g_yzzzzz_xyyyyy_0[i] = 5.0 * g_zzzzz_xyyyy_1[i] * fe_0 + g_zzzzz_xyyyyy_1[i] * pa_y[i];

        g_yzzzzz_xyyyyz_0[i] = 4.0 * g_zzzzz_xyyyz_1[i] * fe_0 + g_zzzzz_xyyyyz_1[i] * pa_y[i];

        g_yzzzzz_xyyyzz_0[i] = 3.0 * g_zzzzz_xyyzz_1[i] * fe_0 + g_zzzzz_xyyyzz_1[i] * pa_y[i];

        g_yzzzzz_xyyzzz_0[i] = 2.0 * g_zzzzz_xyzzz_1[i] * fe_0 + g_zzzzz_xyyzzz_1[i] * pa_y[i];

        g_yzzzzz_xyzzzz_0[i] = g_zzzzz_xzzzz_1[i] * fe_0 + g_zzzzz_xyzzzz_1[i] * pa_y[i];

        g_yzzzzz_xzzzzz_0[i] = g_zzzzz_xzzzzz_1[i] * pa_y[i];

        g_yzzzzz_yyyyyy_0[i] = 6.0 * g_zzzzz_yyyyy_1[i] * fe_0 + g_zzzzz_yyyyyy_1[i] * pa_y[i];

        g_yzzzzz_yyyyyz_0[i] = 5.0 * g_zzzzz_yyyyz_1[i] * fe_0 + g_zzzzz_yyyyyz_1[i] * pa_y[i];

        g_yzzzzz_yyyyzz_0[i] = 4.0 * g_zzzzz_yyyzz_1[i] * fe_0 + g_zzzzz_yyyyzz_1[i] * pa_y[i];

        g_yzzzzz_yyyzzz_0[i] = 3.0 * g_zzzzz_yyzzz_1[i] * fe_0 + g_zzzzz_yyyzzz_1[i] * pa_y[i];

        g_yzzzzz_yyzzzz_0[i] = 2.0 * g_zzzzz_yzzzz_1[i] * fe_0 + g_zzzzz_yyzzzz_1[i] * pa_y[i];

        g_yzzzzz_yzzzzz_0[i] = g_zzzzz_zzzzz_1[i] * fe_0 + g_zzzzz_yzzzzz_1[i] * pa_y[i];

        g_yzzzzz_zzzzzz_0[i] = g_zzzzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 756-784 components of targeted buffer : II

    auto g_zzzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_ii + 756);

    auto g_zzzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_ii + 757);

    auto g_zzzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_ii + 758);

    auto g_zzzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_ii + 759);

    auto g_zzzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_ii + 760);

    auto g_zzzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_ii + 761);

    auto g_zzzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_ii + 762);

    auto g_zzzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_ii + 763);

    auto g_zzzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_ii + 764);

    auto g_zzzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_ii + 765);

    auto g_zzzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_ii + 766);

    auto g_zzzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_ii + 767);

    auto g_zzzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_ii + 768);

    auto g_zzzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_ii + 769);

    auto g_zzzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_ii + 770);

    auto g_zzzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_ii + 771);

    auto g_zzzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_ii + 772);

    auto g_zzzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_ii + 773);

    auto g_zzzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_ii + 774);

    auto g_zzzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_ii + 775);

    auto g_zzzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_ii + 776);

    auto g_zzzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_ii + 777);

    auto g_zzzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_ii + 778);

    auto g_zzzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_ii + 779);

    auto g_zzzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_ii + 780);

    auto g_zzzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_ii + 781);

    auto g_zzzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_ii + 782);

    auto g_zzzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_ii + 783);

    #pragma omp simd aligned(g_zzzz_xxxxxx_0, g_zzzz_xxxxxx_1, g_zzzz_xxxxxy_0, g_zzzz_xxxxxy_1, g_zzzz_xxxxxz_0, g_zzzz_xxxxxz_1, g_zzzz_xxxxyy_0, g_zzzz_xxxxyy_1, g_zzzz_xxxxyz_0, g_zzzz_xxxxyz_1, g_zzzz_xxxxzz_0, g_zzzz_xxxxzz_1, g_zzzz_xxxyyy_0, g_zzzz_xxxyyy_1, g_zzzz_xxxyyz_0, g_zzzz_xxxyyz_1, g_zzzz_xxxyzz_0, g_zzzz_xxxyzz_1, g_zzzz_xxxzzz_0, g_zzzz_xxxzzz_1, g_zzzz_xxyyyy_0, g_zzzz_xxyyyy_1, g_zzzz_xxyyyz_0, g_zzzz_xxyyyz_1, g_zzzz_xxyyzz_0, g_zzzz_xxyyzz_1, g_zzzz_xxyzzz_0, g_zzzz_xxyzzz_1, g_zzzz_xxzzzz_0, g_zzzz_xxzzzz_1, g_zzzz_xyyyyy_0, g_zzzz_xyyyyy_1, g_zzzz_xyyyyz_0, g_zzzz_xyyyyz_1, g_zzzz_xyyyzz_0, g_zzzz_xyyyzz_1, g_zzzz_xyyzzz_0, g_zzzz_xyyzzz_1, g_zzzz_xyzzzz_0, g_zzzz_xyzzzz_1, g_zzzz_xzzzzz_0, g_zzzz_xzzzzz_1, g_zzzz_yyyyyy_0, g_zzzz_yyyyyy_1, g_zzzz_yyyyyz_0, g_zzzz_yyyyyz_1, g_zzzz_yyyyzz_0, g_zzzz_yyyyzz_1, g_zzzz_yyyzzz_0, g_zzzz_yyyzzz_1, g_zzzz_yyzzzz_0, g_zzzz_yyzzzz_1, g_zzzz_yzzzzz_0, g_zzzz_yzzzzz_1, g_zzzz_zzzzzz_0, g_zzzz_zzzzzz_1, g_zzzzz_xxxxx_1, g_zzzzz_xxxxxx_1, g_zzzzz_xxxxxy_1, g_zzzzz_xxxxxz_1, g_zzzzz_xxxxy_1, g_zzzzz_xxxxyy_1, g_zzzzz_xxxxyz_1, g_zzzzz_xxxxz_1, g_zzzzz_xxxxzz_1, g_zzzzz_xxxyy_1, g_zzzzz_xxxyyy_1, g_zzzzz_xxxyyz_1, g_zzzzz_xxxyz_1, g_zzzzz_xxxyzz_1, g_zzzzz_xxxzz_1, g_zzzzz_xxxzzz_1, g_zzzzz_xxyyy_1, g_zzzzz_xxyyyy_1, g_zzzzz_xxyyyz_1, g_zzzzz_xxyyz_1, g_zzzzz_xxyyzz_1, g_zzzzz_xxyzz_1, g_zzzzz_xxyzzz_1, g_zzzzz_xxzzz_1, g_zzzzz_xxzzzz_1, g_zzzzz_xyyyy_1, g_zzzzz_xyyyyy_1, g_zzzzz_xyyyyz_1, g_zzzzz_xyyyz_1, g_zzzzz_xyyyzz_1, g_zzzzz_xyyzz_1, g_zzzzz_xyyzzz_1, g_zzzzz_xyzzz_1, g_zzzzz_xyzzzz_1, g_zzzzz_xzzzz_1, g_zzzzz_xzzzzz_1, g_zzzzz_yyyyy_1, g_zzzzz_yyyyyy_1, g_zzzzz_yyyyyz_1, g_zzzzz_yyyyz_1, g_zzzzz_yyyyzz_1, g_zzzzz_yyyzz_1, g_zzzzz_yyyzzz_1, g_zzzzz_yyzzz_1, g_zzzzz_yyzzzz_1, g_zzzzz_yzzzz_1, g_zzzzz_yzzzzz_1, g_zzzzz_zzzzz_1, g_zzzzz_zzzzzz_1, g_zzzzzz_xxxxxx_0, g_zzzzzz_xxxxxy_0, g_zzzzzz_xxxxxz_0, g_zzzzzz_xxxxyy_0, g_zzzzzz_xxxxyz_0, g_zzzzzz_xxxxzz_0, g_zzzzzz_xxxyyy_0, g_zzzzzz_xxxyyz_0, g_zzzzzz_xxxyzz_0, g_zzzzzz_xxxzzz_0, g_zzzzzz_xxyyyy_0, g_zzzzzz_xxyyyz_0, g_zzzzzz_xxyyzz_0, g_zzzzzz_xxyzzz_0, g_zzzzzz_xxzzzz_0, g_zzzzzz_xyyyyy_0, g_zzzzzz_xyyyyz_0, g_zzzzzz_xyyyzz_0, g_zzzzzz_xyyzzz_0, g_zzzzzz_xyzzzz_0, g_zzzzzz_xzzzzz_0, g_zzzzzz_yyyyyy_0, g_zzzzzz_yyyyyz_0, g_zzzzzz_yyyyzz_0, g_zzzzzz_yyyzzz_0, g_zzzzzz_yyzzzz_0, g_zzzzzz_yzzzzz_0, g_zzzzzz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzzz_xxxxxx_0[i] = 5.0 * g_zzzz_xxxxxx_0[i] * fbe_0 - 5.0 * g_zzzz_xxxxxx_1[i] * fz_be_0 + g_zzzzz_xxxxxx_1[i] * pa_z[i];

        g_zzzzzz_xxxxxy_0[i] = 5.0 * g_zzzz_xxxxxy_0[i] * fbe_0 - 5.0 * g_zzzz_xxxxxy_1[i] * fz_be_0 + g_zzzzz_xxxxxy_1[i] * pa_z[i];

        g_zzzzzz_xxxxxz_0[i] = 5.0 * g_zzzz_xxxxxz_0[i] * fbe_0 - 5.0 * g_zzzz_xxxxxz_1[i] * fz_be_0 + g_zzzzz_xxxxx_1[i] * fe_0 + g_zzzzz_xxxxxz_1[i] * pa_z[i];

        g_zzzzzz_xxxxyy_0[i] = 5.0 * g_zzzz_xxxxyy_0[i] * fbe_0 - 5.0 * g_zzzz_xxxxyy_1[i] * fz_be_0 + g_zzzzz_xxxxyy_1[i] * pa_z[i];

        g_zzzzzz_xxxxyz_0[i] = 5.0 * g_zzzz_xxxxyz_0[i] * fbe_0 - 5.0 * g_zzzz_xxxxyz_1[i] * fz_be_0 + g_zzzzz_xxxxy_1[i] * fe_0 + g_zzzzz_xxxxyz_1[i] * pa_z[i];

        g_zzzzzz_xxxxzz_0[i] = 5.0 * g_zzzz_xxxxzz_0[i] * fbe_0 - 5.0 * g_zzzz_xxxxzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_xxxxz_1[i] * fe_0 + g_zzzzz_xxxxzz_1[i] * pa_z[i];

        g_zzzzzz_xxxyyy_0[i] = 5.0 * g_zzzz_xxxyyy_0[i] * fbe_0 - 5.0 * g_zzzz_xxxyyy_1[i] * fz_be_0 + g_zzzzz_xxxyyy_1[i] * pa_z[i];

        g_zzzzzz_xxxyyz_0[i] = 5.0 * g_zzzz_xxxyyz_0[i] * fbe_0 - 5.0 * g_zzzz_xxxyyz_1[i] * fz_be_0 + g_zzzzz_xxxyy_1[i] * fe_0 + g_zzzzz_xxxyyz_1[i] * pa_z[i];

        g_zzzzzz_xxxyzz_0[i] = 5.0 * g_zzzz_xxxyzz_0[i] * fbe_0 - 5.0 * g_zzzz_xxxyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_xxxyz_1[i] * fe_0 + g_zzzzz_xxxyzz_1[i] * pa_z[i];

        g_zzzzzz_xxxzzz_0[i] = 5.0 * g_zzzz_xxxzzz_0[i] * fbe_0 - 5.0 * g_zzzz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_xxxzz_1[i] * fe_0 + g_zzzzz_xxxzzz_1[i] * pa_z[i];

        g_zzzzzz_xxyyyy_0[i] = 5.0 * g_zzzz_xxyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_xxyyyy_1[i] * fz_be_0 + g_zzzzz_xxyyyy_1[i] * pa_z[i];

        g_zzzzzz_xxyyyz_0[i] = 5.0 * g_zzzz_xxyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_xxyyyz_1[i] * fz_be_0 + g_zzzzz_xxyyy_1[i] * fe_0 + g_zzzzz_xxyyyz_1[i] * pa_z[i];

        g_zzzzzz_xxyyzz_0[i] = 5.0 * g_zzzz_xxyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_xxyyz_1[i] * fe_0 + g_zzzzz_xxyyzz_1[i] * pa_z[i];

        g_zzzzzz_xxyzzz_0[i] = 5.0 * g_zzzz_xxyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_xxyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_xxyzz_1[i] * fe_0 + g_zzzzz_xxyzzz_1[i] * pa_z[i];

        g_zzzzzz_xxzzzz_0[i] = 5.0 * g_zzzz_xxzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_xxzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_xxzzz_1[i] * fe_0 + g_zzzzz_xxzzzz_1[i] * pa_z[i];

        g_zzzzzz_xyyyyy_0[i] = 5.0 * g_zzzz_xyyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_xyyyyy_1[i] * fz_be_0 + g_zzzzz_xyyyyy_1[i] * pa_z[i];

        g_zzzzzz_xyyyyz_0[i] = 5.0 * g_zzzz_xyyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_xyyyyz_1[i] * fz_be_0 + g_zzzzz_xyyyy_1[i] * fe_0 + g_zzzzz_xyyyyz_1[i] * pa_z[i];

        g_zzzzzz_xyyyzz_0[i] = 5.0 * g_zzzz_xyyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_xyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_xyyyz_1[i] * fe_0 + g_zzzzz_xyyyzz_1[i] * pa_z[i];

        g_zzzzzz_xyyzzz_0[i] = 5.0 * g_zzzz_xyyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_xyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_xyyzz_1[i] * fe_0 + g_zzzzz_xyyzzz_1[i] * pa_z[i];

        g_zzzzzz_xyzzzz_0[i] = 5.0 * g_zzzz_xyzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_xyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_xyzzz_1[i] * fe_0 + g_zzzzz_xyzzzz_1[i] * pa_z[i];

        g_zzzzzz_xzzzzz_0[i] = 5.0 * g_zzzz_xzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_xzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_xzzzz_1[i] * fe_0 + g_zzzzz_xzzzzz_1[i] * pa_z[i];

        g_zzzzzz_yyyyyy_0[i] = 5.0 * g_zzzz_yyyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_yyyyyy_1[i] * fz_be_0 + g_zzzzz_yyyyyy_1[i] * pa_z[i];

        g_zzzzzz_yyyyyz_0[i] = 5.0 * g_zzzz_yyyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_yyyyyz_1[i] * fz_be_0 + g_zzzzz_yyyyy_1[i] * fe_0 + g_zzzzz_yyyyyz_1[i] * pa_z[i];

        g_zzzzzz_yyyyzz_0[i] = 5.0 * g_zzzz_yyyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_yyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_yyyyz_1[i] * fe_0 + g_zzzzz_yyyyzz_1[i] * pa_z[i];

        g_zzzzzz_yyyzzz_0[i] = 5.0 * g_zzzz_yyyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_yyyzz_1[i] * fe_0 + g_zzzzz_yyyzzz_1[i] * pa_z[i];

        g_zzzzzz_yyzzzz_0[i] = 5.0 * g_zzzz_yyzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_yyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_yyzzz_1[i] * fe_0 + g_zzzzz_yyzzzz_1[i] * pa_z[i];

        g_zzzzzz_yzzzzz_0[i] = 5.0 * g_zzzz_yzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_yzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_yzzzz_1[i] * fe_0 + g_zzzzz_yzzzzz_1[i] * pa_z[i];

        g_zzzzzz_zzzzzz_0[i] = 5.0 * g_zzzz_zzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_zzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzz_zzzzz_1[i] * fe_0 + g_zzzzz_zzzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

