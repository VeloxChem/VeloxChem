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

#include "ElectronRepulsionPrimRecSISI.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sisi(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sisi,
                                  size_t                idx_eri_0_sgsi,
                                  size_t                idx_eri_1_sgsi,
                                  size_t                idx_eri_1_shsh,
                                  size_t                idx_eri_0_shsi,
                                  size_t                idx_eri_1_shsi,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SGSI

    auto g_0_xxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi);

    auto g_0_xxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 1);

    auto g_0_xxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 2);

    auto g_0_xxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 3);

    auto g_0_xxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 4);

    auto g_0_xxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 5);

    auto g_0_xxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 6);

    auto g_0_xxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 7);

    auto g_0_xxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 8);

    auto g_0_xxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 9);

    auto g_0_xxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 10);

    auto g_0_xxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 11);

    auto g_0_xxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 12);

    auto g_0_xxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 13);

    auto g_0_xxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 14);

    auto g_0_xxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 15);

    auto g_0_xxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 16);

    auto g_0_xxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 17);

    auto g_0_xxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 18);

    auto g_0_xxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 19);

    auto g_0_xxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 20);

    auto g_0_xxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 21);

    auto g_0_xxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 22);

    auto g_0_xxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 23);

    auto g_0_xxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 24);

    auto g_0_xxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 25);

    auto g_0_xxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 26);

    auto g_0_xxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 27);

    auto g_0_xxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 28);

    auto g_0_xxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 30);

    auto g_0_xxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 33);

    auto g_0_xxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 37);

    auto g_0_xxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 42);

    auto g_0_xxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 48);

    auto g_0_xxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 56);

    auto g_0_xxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 57);

    auto g_0_xxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 59);

    auto g_0_xxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 62);

    auto g_0_xxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 66);

    auto g_0_xxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 71);

    auto g_0_xxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 84);

    auto g_0_xxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 85);

    auto g_0_xxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 86);

    auto g_0_xxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 87);

    auto g_0_xxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 88);

    auto g_0_xxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 89);

    auto g_0_xxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 90);

    auto g_0_xxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 91);

    auto g_0_xxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 92);

    auto g_0_xxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 93);

    auto g_0_xxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 94);

    auto g_0_xxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 95);

    auto g_0_xxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 96);

    auto g_0_xxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 97);

    auto g_0_xxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 98);

    auto g_0_xxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 99);

    auto g_0_xxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 100);

    auto g_0_xxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 101);

    auto g_0_xxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 102);

    auto g_0_xxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 103);

    auto g_0_xxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 104);

    auto g_0_xxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 105);

    auto g_0_xxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 106);

    auto g_0_xxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 107);

    auto g_0_xxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 108);

    auto g_0_xxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 109);

    auto g_0_xxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 110);

    auto g_0_xxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 111);

    auto g_0_xxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 140);

    auto g_0_xxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 141);

    auto g_0_xxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 142);

    auto g_0_xxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 143);

    auto g_0_xxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 144);

    auto g_0_xxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 145);

    auto g_0_xxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 146);

    auto g_0_xxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 147);

    auto g_0_xxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 148);

    auto g_0_xxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 149);

    auto g_0_xxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 150);

    auto g_0_xxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 151);

    auto g_0_xxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 152);

    auto g_0_xxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 153);

    auto g_0_xxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 154);

    auto g_0_xxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 155);

    auto g_0_xxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 156);

    auto g_0_xxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 157);

    auto g_0_xxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 158);

    auto g_0_xxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 159);

    auto g_0_xxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 160);

    auto g_0_xxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 161);

    auto g_0_xxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 162);

    auto g_0_xxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 163);

    auto g_0_xxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 164);

    auto g_0_xxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 165);

    auto g_0_xxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 166);

    auto g_0_xxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 167);

    auto g_0_xyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 169);

    auto g_0_xyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 171);

    auto g_0_xyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 172);

    auto g_0_xyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 174);

    auto g_0_xyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 175);

    auto g_0_xyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 176);

    auto g_0_xyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 178);

    auto g_0_xyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 179);

    auto g_0_xyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 180);

    auto g_0_xyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 181);

    auto g_0_xyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 183);

    auto g_0_xyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 184);

    auto g_0_xyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 185);

    auto g_0_xyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 186);

    auto g_0_xyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 187);

    auto g_0_xyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 189);

    auto g_0_xyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 190);

    auto g_0_xyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 191);

    auto g_0_xyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 192);

    auto g_0_xyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 193);

    auto g_0_xyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 194);

    auto g_0_xyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 195);

    auto g_0_xzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 254);

    auto g_0_xzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 256);

    auto g_0_xzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 257);

    auto g_0_xzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 259);

    auto g_0_xzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 260);

    auto g_0_xzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 261);

    auto g_0_xzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 263);

    auto g_0_xzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 264);

    auto g_0_xzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 265);

    auto g_0_xzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 266);

    auto g_0_xzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 268);

    auto g_0_xzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 269);

    auto g_0_xzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 270);

    auto g_0_xzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 271);

    auto g_0_xzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 272);

    auto g_0_xzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 273);

    auto g_0_xzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 274);

    auto g_0_xzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 275);

    auto g_0_xzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 276);

    auto g_0_xzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 277);

    auto g_0_xzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 278);

    auto g_0_xzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 279);

    auto g_0_yyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 280);

    auto g_0_yyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 281);

    auto g_0_yyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 282);

    auto g_0_yyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 283);

    auto g_0_yyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 284);

    auto g_0_yyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 285);

    auto g_0_yyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 286);

    auto g_0_yyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 287);

    auto g_0_yyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 288);

    auto g_0_yyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 289);

    auto g_0_yyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 290);

    auto g_0_yyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 291);

    auto g_0_yyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 292);

    auto g_0_yyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 293);

    auto g_0_yyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 294);

    auto g_0_yyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 295);

    auto g_0_yyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 296);

    auto g_0_yyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 297);

    auto g_0_yyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 298);

    auto g_0_yyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 299);

    auto g_0_yyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 300);

    auto g_0_yyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 301);

    auto g_0_yyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 302);

    auto g_0_yyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 303);

    auto g_0_yyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 304);

    auto g_0_yyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 305);

    auto g_0_yyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 306);

    auto g_0_yyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 307);

    auto g_0_yyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 309);

    auto g_0_yyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 311);

    auto g_0_yyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 314);

    auto g_0_yyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 318);

    auto g_0_yyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 323);

    auto g_0_yyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 329);

    auto g_0_yyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 336);

    auto g_0_yyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 337);

    auto g_0_yyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 338);

    auto g_0_yyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 339);

    auto g_0_yyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 340);

    auto g_0_yyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 341);

    auto g_0_yyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 342);

    auto g_0_yyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 343);

    auto g_0_yyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 344);

    auto g_0_yyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 345);

    auto g_0_yyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 346);

    auto g_0_yyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 347);

    auto g_0_yyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 348);

    auto g_0_yyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 349);

    auto g_0_yyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 350);

    auto g_0_yyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 351);

    auto g_0_yyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 352);

    auto g_0_yyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 353);

    auto g_0_yyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 354);

    auto g_0_yyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 355);

    auto g_0_yyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 356);

    auto g_0_yyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 357);

    auto g_0_yyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 358);

    auto g_0_yyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 359);

    auto g_0_yyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 360);

    auto g_0_yyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 361);

    auto g_0_yyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 362);

    auto g_0_yyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 363);

    auto g_0_yzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 364);

    auto g_0_yzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 366);

    auto g_0_yzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 368);

    auto g_0_yzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 369);

    auto g_0_yzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 371);

    auto g_0_yzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 372);

    auto g_0_yzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 373);

    auto g_0_yzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 375);

    auto g_0_yzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 376);

    auto g_0_yzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 377);

    auto g_0_yzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 378);

    auto g_0_yzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 380);

    auto g_0_yzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 381);

    auto g_0_yzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 382);

    auto g_0_yzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 383);

    auto g_0_yzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 384);

    auto g_0_yzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 386);

    auto g_0_yzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 387);

    auto g_0_yzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 388);

    auto g_0_yzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 389);

    auto g_0_yzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 390);

    auto g_0_yzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 391);

    auto g_0_zzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sgsi + 392);

    auto g_0_zzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sgsi + 393);

    auto g_0_zzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sgsi + 394);

    auto g_0_zzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sgsi + 395);

    auto g_0_zzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sgsi + 396);

    auto g_0_zzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sgsi + 397);

    auto g_0_zzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sgsi + 398);

    auto g_0_zzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sgsi + 399);

    auto g_0_zzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sgsi + 400);

    auto g_0_zzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sgsi + 401);

    auto g_0_zzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 402);

    auto g_0_zzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 403);

    auto g_0_zzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 404);

    auto g_0_zzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 405);

    auto g_0_zzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 406);

    auto g_0_zzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 407);

    auto g_0_zzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 408);

    auto g_0_zzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 409);

    auto g_0_zzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 410);

    auto g_0_zzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 411);

    auto g_0_zzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 412);

    auto g_0_zzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sgsi + 413);

    auto g_0_zzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sgsi + 414);

    auto g_0_zzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sgsi + 415);

    auto g_0_zzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sgsi + 416);

    auto g_0_zzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 417);

    auto g_0_zzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 418);

    auto g_0_zzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sgsi + 419);

    /// Set up components of auxilary buffer : SGSI

    auto g_0_xxxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi);

    auto g_0_xxxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 1);

    auto g_0_xxxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 2);

    auto g_0_xxxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 3);

    auto g_0_xxxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 4);

    auto g_0_xxxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 5);

    auto g_0_xxxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 6);

    auto g_0_xxxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 7);

    auto g_0_xxxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 8);

    auto g_0_xxxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 9);

    auto g_0_xxxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 10);

    auto g_0_xxxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 11);

    auto g_0_xxxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 12);

    auto g_0_xxxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 13);

    auto g_0_xxxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 14);

    auto g_0_xxxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 15);

    auto g_0_xxxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 16);

    auto g_0_xxxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 17);

    auto g_0_xxxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 18);

    auto g_0_xxxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 19);

    auto g_0_xxxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 20);

    auto g_0_xxxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 21);

    auto g_0_xxxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 22);

    auto g_0_xxxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 23);

    auto g_0_xxxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 24);

    auto g_0_xxxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 25);

    auto g_0_xxxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 26);

    auto g_0_xxxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 27);

    auto g_0_xxxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 28);

    auto g_0_xxxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 30);

    auto g_0_xxxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 33);

    auto g_0_xxxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 37);

    auto g_0_xxxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 42);

    auto g_0_xxxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 48);

    auto g_0_xxxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 56);

    auto g_0_xxxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 57);

    auto g_0_xxxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 59);

    auto g_0_xxxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 62);

    auto g_0_xxxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 66);

    auto g_0_xxxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 71);

    auto g_0_xxyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 84);

    auto g_0_xxyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 85);

    auto g_0_xxyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 86);

    auto g_0_xxyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 87);

    auto g_0_xxyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 88);

    auto g_0_xxyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 89);

    auto g_0_xxyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 90);

    auto g_0_xxyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 91);

    auto g_0_xxyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 92);

    auto g_0_xxyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 93);

    auto g_0_xxyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 94);

    auto g_0_xxyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 95);

    auto g_0_xxyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 96);

    auto g_0_xxyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 97);

    auto g_0_xxyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 98);

    auto g_0_xxyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 99);

    auto g_0_xxyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 100);

    auto g_0_xxyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 101);

    auto g_0_xxyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 102);

    auto g_0_xxyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 103);

    auto g_0_xxyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 104);

    auto g_0_xxyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 105);

    auto g_0_xxyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 106);

    auto g_0_xxyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 107);

    auto g_0_xxyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 108);

    auto g_0_xxyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 109);

    auto g_0_xxyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 110);

    auto g_0_xxyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 111);

    auto g_0_xxzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 140);

    auto g_0_xxzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 141);

    auto g_0_xxzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 142);

    auto g_0_xxzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 143);

    auto g_0_xxzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 144);

    auto g_0_xxzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 145);

    auto g_0_xxzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 146);

    auto g_0_xxzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 147);

    auto g_0_xxzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 148);

    auto g_0_xxzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 149);

    auto g_0_xxzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 150);

    auto g_0_xxzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 151);

    auto g_0_xxzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 152);

    auto g_0_xxzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 153);

    auto g_0_xxzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 154);

    auto g_0_xxzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 155);

    auto g_0_xxzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 156);

    auto g_0_xxzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 157);

    auto g_0_xxzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 158);

    auto g_0_xxzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 159);

    auto g_0_xxzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 160);

    auto g_0_xxzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 161);

    auto g_0_xxzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 162);

    auto g_0_xxzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 163);

    auto g_0_xxzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 164);

    auto g_0_xxzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 165);

    auto g_0_xxzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 166);

    auto g_0_xxzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 167);

    auto g_0_xyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 169);

    auto g_0_xyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 171);

    auto g_0_xyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 172);

    auto g_0_xyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 174);

    auto g_0_xyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 175);

    auto g_0_xyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 176);

    auto g_0_xyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 178);

    auto g_0_xyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 179);

    auto g_0_xyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 180);

    auto g_0_xyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 181);

    auto g_0_xyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 183);

    auto g_0_xyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 184);

    auto g_0_xyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 185);

    auto g_0_xyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 186);

    auto g_0_xyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 187);

    auto g_0_xyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 189);

    auto g_0_xyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 190);

    auto g_0_xyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 191);

    auto g_0_xyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 192);

    auto g_0_xyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 193);

    auto g_0_xyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 194);

    auto g_0_xyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 195);

    auto g_0_xzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 254);

    auto g_0_xzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 256);

    auto g_0_xzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 257);

    auto g_0_xzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 259);

    auto g_0_xzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 260);

    auto g_0_xzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 261);

    auto g_0_xzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 263);

    auto g_0_xzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 264);

    auto g_0_xzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 265);

    auto g_0_xzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 266);

    auto g_0_xzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 268);

    auto g_0_xzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 269);

    auto g_0_xzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 270);

    auto g_0_xzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 271);

    auto g_0_xzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 272);

    auto g_0_xzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 273);

    auto g_0_xzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 274);

    auto g_0_xzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 275);

    auto g_0_xzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 276);

    auto g_0_xzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 277);

    auto g_0_xzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 278);

    auto g_0_xzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 279);

    auto g_0_yyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 280);

    auto g_0_yyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 281);

    auto g_0_yyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 282);

    auto g_0_yyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 283);

    auto g_0_yyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 284);

    auto g_0_yyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 285);

    auto g_0_yyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 286);

    auto g_0_yyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 287);

    auto g_0_yyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 288);

    auto g_0_yyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 289);

    auto g_0_yyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 290);

    auto g_0_yyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 291);

    auto g_0_yyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 292);

    auto g_0_yyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 293);

    auto g_0_yyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 294);

    auto g_0_yyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 295);

    auto g_0_yyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 296);

    auto g_0_yyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 297);

    auto g_0_yyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 298);

    auto g_0_yyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 299);

    auto g_0_yyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 300);

    auto g_0_yyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 301);

    auto g_0_yyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 302);

    auto g_0_yyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 303);

    auto g_0_yyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 304);

    auto g_0_yyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 305);

    auto g_0_yyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 306);

    auto g_0_yyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 307);

    auto g_0_yyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 309);

    auto g_0_yyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 311);

    auto g_0_yyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 314);

    auto g_0_yyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 318);

    auto g_0_yyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 323);

    auto g_0_yyyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 329);

    auto g_0_yyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 336);

    auto g_0_yyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 337);

    auto g_0_yyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 338);

    auto g_0_yyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 339);

    auto g_0_yyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 340);

    auto g_0_yyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 341);

    auto g_0_yyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 342);

    auto g_0_yyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 343);

    auto g_0_yyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 344);

    auto g_0_yyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 345);

    auto g_0_yyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 346);

    auto g_0_yyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 347);

    auto g_0_yyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 348);

    auto g_0_yyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 349);

    auto g_0_yyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 350);

    auto g_0_yyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 351);

    auto g_0_yyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 352);

    auto g_0_yyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 353);

    auto g_0_yyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 354);

    auto g_0_yyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 355);

    auto g_0_yyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 356);

    auto g_0_yyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 357);

    auto g_0_yyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 358);

    auto g_0_yyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 359);

    auto g_0_yyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 360);

    auto g_0_yyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 361);

    auto g_0_yyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 362);

    auto g_0_yyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 363);

    auto g_0_yzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 364);

    auto g_0_yzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 366);

    auto g_0_yzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 368);

    auto g_0_yzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 369);

    auto g_0_yzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 371);

    auto g_0_yzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 372);

    auto g_0_yzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 373);

    auto g_0_yzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 375);

    auto g_0_yzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 376);

    auto g_0_yzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 377);

    auto g_0_yzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 378);

    auto g_0_yzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 380);

    auto g_0_yzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 381);

    auto g_0_yzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 382);

    auto g_0_yzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 383);

    auto g_0_yzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 384);

    auto g_0_yzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 386);

    auto g_0_yzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 387);

    auto g_0_yzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 388);

    auto g_0_yzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 389);

    auto g_0_yzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 390);

    auto g_0_yzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 391);

    auto g_0_zzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sgsi + 392);

    auto g_0_zzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sgsi + 393);

    auto g_0_zzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sgsi + 394);

    auto g_0_zzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sgsi + 395);

    auto g_0_zzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sgsi + 396);

    auto g_0_zzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sgsi + 397);

    auto g_0_zzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sgsi + 398);

    auto g_0_zzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sgsi + 399);

    auto g_0_zzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sgsi + 400);

    auto g_0_zzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sgsi + 401);

    auto g_0_zzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 402);

    auto g_0_zzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 403);

    auto g_0_zzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 404);

    auto g_0_zzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 405);

    auto g_0_zzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 406);

    auto g_0_zzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 407);

    auto g_0_zzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 408);

    auto g_0_zzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 409);

    auto g_0_zzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 410);

    auto g_0_zzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 411);

    auto g_0_zzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 412);

    auto g_0_zzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sgsi + 413);

    auto g_0_zzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sgsi + 414);

    auto g_0_zzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sgsi + 415);

    auto g_0_zzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sgsi + 416);

    auto g_0_zzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 417);

    auto g_0_zzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 418);

    auto g_0_zzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sgsi + 419);

    /// Set up components of auxilary buffer : SHSH

    auto g_0_xxxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh);

    auto g_0_xxxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 1);

    auto g_0_xxxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 2);

    auto g_0_xxxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 3);

    auto g_0_xxxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 4);

    auto g_0_xxxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 5);

    auto g_0_xxxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 6);

    auto g_0_xxxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 7);

    auto g_0_xxxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 8);

    auto g_0_xxxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 9);

    auto g_0_xxxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 10);

    auto g_0_xxxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 11);

    auto g_0_xxxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 12);

    auto g_0_xxxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 13);

    auto g_0_xxxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 14);

    auto g_0_xxxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 15);

    auto g_0_xxxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 16);

    auto g_0_xxxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 17);

    auto g_0_xxxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 18);

    auto g_0_xxxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 19);

    auto g_0_xxxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 20);

    auto g_0_xxxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 44);

    auto g_0_xxxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 46);

    auto g_0_xxxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 47);

    auto g_0_xxxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 49);

    auto g_0_xxxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 50);

    auto g_0_xxxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 51);

    auto g_0_xxxxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 53);

    auto g_0_xxxxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 54);

    auto g_0_xxxxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 55);

    auto g_0_xxxxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 56);

    auto g_0_xxxxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 58);

    auto g_0_xxxxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 59);

    auto g_0_xxxxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 60);

    auto g_0_xxxxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 61);

    auto g_0_xxxxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 62);

    auto g_0_xxxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 63);

    auto g_0_xxxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 64);

    auto g_0_xxxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 65);

    auto g_0_xxxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 66);

    auto g_0_xxxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 67);

    auto g_0_xxxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 68);

    auto g_0_xxxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 69);

    auto g_0_xxxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 70);

    auto g_0_xxxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 71);

    auto g_0_xxxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 72);

    auto g_0_xxxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 73);

    auto g_0_xxxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 74);

    auto g_0_xxxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 75);

    auto g_0_xxxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 76);

    auto g_0_xxxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 77);

    auto g_0_xxxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 78);

    auto g_0_xxxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 79);

    auto g_0_xxxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 80);

    auto g_0_xxxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 81);

    auto g_0_xxxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 82);

    auto g_0_xxxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 83);

    auto g_0_xxxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 105);

    auto g_0_xxxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 106);

    auto g_0_xxxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 107);

    auto g_0_xxxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 108);

    auto g_0_xxxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 109);

    auto g_0_xxxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 110);

    auto g_0_xxxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 111);

    auto g_0_xxxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 112);

    auto g_0_xxxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 113);

    auto g_0_xxxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 114);

    auto g_0_xxxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 115);

    auto g_0_xxxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 116);

    auto g_0_xxxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 117);

    auto g_0_xxxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 118);

    auto g_0_xxxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 119);

    auto g_0_xxxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 120);

    auto g_0_xxxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 121);

    auto g_0_xxxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 122);

    auto g_0_xxxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 123);

    auto g_0_xxxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 124);

    auto g_0_xxxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 125);

    auto g_0_xxyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 126);

    auto g_0_xxyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 127);

    auto g_0_xxyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 128);

    auto g_0_xxyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 129);

    auto g_0_xxyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 130);

    auto g_0_xxyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 131);

    auto g_0_xxyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 132);

    auto g_0_xxyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 133);

    auto g_0_xxyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 134);

    auto g_0_xxyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 135);

    auto g_0_xxyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 136);

    auto g_0_xxyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 137);

    auto g_0_xxyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 138);

    auto g_0_xxyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 139);

    auto g_0_xxyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 140);

    auto g_0_xxyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 141);

    auto g_0_xxyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 142);

    auto g_0_xxyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 143);

    auto g_0_xxyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 144);

    auto g_0_xxyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 145);

    auto g_0_xxyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 146);

    auto g_0_xxzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 189);

    auto g_0_xxzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 190);

    auto g_0_xxzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 191);

    auto g_0_xxzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 192);

    auto g_0_xxzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 193);

    auto g_0_xxzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 194);

    auto g_0_xxzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 195);

    auto g_0_xxzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 196);

    auto g_0_xxzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 197);

    auto g_0_xxzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 198);

    auto g_0_xxzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 199);

    auto g_0_xxzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 200);

    auto g_0_xxzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 201);

    auto g_0_xxzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 202);

    auto g_0_xxzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 203);

    auto g_0_xxzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 204);

    auto g_0_xxzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 205);

    auto g_0_xxzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 206);

    auto g_0_xxzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 207);

    auto g_0_xxzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 208);

    auto g_0_xxzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 209);

    auto g_0_xyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 211);

    auto g_0_xyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 213);

    auto g_0_xyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 214);

    auto g_0_xyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 216);

    auto g_0_xyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 217);

    auto g_0_xyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 218);

    auto g_0_xyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 220);

    auto g_0_xyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 221);

    auto g_0_xyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 222);

    auto g_0_xyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 223);

    auto g_0_xyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 225);

    auto g_0_xyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 226);

    auto g_0_xyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 227);

    auto g_0_xyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 228);

    auto g_0_xyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 229);

    auto g_0_xyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 256);

    auto g_0_xyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 259);

    auto g_0_xyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 260);

    auto g_0_xyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 263);

    auto g_0_xyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 264);

    auto g_0_xyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 265);

    auto g_0_xyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 268);

    auto g_0_xyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 269);

    auto g_0_xyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 270);

    auto g_0_xyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 271);

    auto g_0_xzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 296);

    auto g_0_xzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 298);

    auto g_0_xzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 299);

    auto g_0_xzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 301);

    auto g_0_xzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 302);

    auto g_0_xzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 303);

    auto g_0_xzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 305);

    auto g_0_xzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 306);

    auto g_0_xzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 307);

    auto g_0_xzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 308);

    auto g_0_xzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 310);

    auto g_0_xzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 311);

    auto g_0_xzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 312);

    auto g_0_xzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 313);

    auto g_0_xzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 314);

    auto g_0_yyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 315);

    auto g_0_yyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 316);

    auto g_0_yyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 317);

    auto g_0_yyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 318);

    auto g_0_yyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 319);

    auto g_0_yyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 320);

    auto g_0_yyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 321);

    auto g_0_yyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 322);

    auto g_0_yyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 323);

    auto g_0_yyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 324);

    auto g_0_yyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 325);

    auto g_0_yyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 326);

    auto g_0_yyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 327);

    auto g_0_yyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 328);

    auto g_0_yyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 329);

    auto g_0_yyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 330);

    auto g_0_yyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 331);

    auto g_0_yyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 332);

    auto g_0_yyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 333);

    auto g_0_yyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 334);

    auto g_0_yyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 335);

    auto g_0_yyyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 338);

    auto g_0_yyyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 340);

    auto g_0_yyyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 341);

    auto g_0_yyyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 343);

    auto g_0_yyyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 344);

    auto g_0_yyyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 345);

    auto g_0_yyyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 347);

    auto g_0_yyyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 348);

    auto g_0_yyyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 349);

    auto g_0_yyyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 350);

    auto g_0_yyyyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 352);

    auto g_0_yyyyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 353);

    auto g_0_yyyyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 354);

    auto g_0_yyyyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 355);

    auto g_0_yyyyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 356);

    auto g_0_yyyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 357);

    auto g_0_yyyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 358);

    auto g_0_yyyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 359);

    auto g_0_yyyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 360);

    auto g_0_yyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 361);

    auto g_0_yyyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 362);

    auto g_0_yyyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 363);

    auto g_0_yyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 364);

    auto g_0_yyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 365);

    auto g_0_yyyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 366);

    auto g_0_yyyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 367);

    auto g_0_yyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 368);

    auto g_0_yyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 369);

    auto g_0_yyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 370);

    auto g_0_yyyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 371);

    auto g_0_yyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 372);

    auto g_0_yyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 373);

    auto g_0_yyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 374);

    auto g_0_yyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 375);

    auto g_0_yyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 376);

    auto g_0_yyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 377);

    auto g_0_yyzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 378);

    auto g_0_yyzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 379);

    auto g_0_yyzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 380);

    auto g_0_yyzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 381);

    auto g_0_yyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 382);

    auto g_0_yyzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 383);

    auto g_0_yyzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 384);

    auto g_0_yyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 385);

    auto g_0_yyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 386);

    auto g_0_yyzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 387);

    auto g_0_yyzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 388);

    auto g_0_yyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 389);

    auto g_0_yyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 390);

    auto g_0_yyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 391);

    auto g_0_yyzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 392);

    auto g_0_yyzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 393);

    auto g_0_yyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 394);

    auto g_0_yyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 395);

    auto g_0_yyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 396);

    auto g_0_yyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 397);

    auto g_0_yyzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 398);

    auto g_0_yzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 400);

    auto g_0_yzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 401);

    auto g_0_yzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 402);

    auto g_0_yzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 403);

    auto g_0_yzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 404);

    auto g_0_yzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 405);

    auto g_0_yzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 406);

    auto g_0_yzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 407);

    auto g_0_yzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 408);

    auto g_0_yzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 409);

    auto g_0_yzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 410);

    auto g_0_yzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 411);

    auto g_0_yzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 412);

    auto g_0_yzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 413);

    auto g_0_yzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 414);

    auto g_0_yzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 415);

    auto g_0_yzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 416);

    auto g_0_yzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 417);

    auto g_0_yzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 418);

    auto g_0_yzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 419);

    auto g_0_zzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_shsh + 420);

    auto g_0_zzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_shsh + 421);

    auto g_0_zzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_shsh + 422);

    auto g_0_zzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_shsh + 423);

    auto g_0_zzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_shsh + 424);

    auto g_0_zzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_shsh + 425);

    auto g_0_zzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_shsh + 426);

    auto g_0_zzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_shsh + 427);

    auto g_0_zzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_shsh + 428);

    auto g_0_zzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_shsh + 429);

    auto g_0_zzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_shsh + 430);

    auto g_0_zzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_shsh + 431);

    auto g_0_zzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_shsh + 432);

    auto g_0_zzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_shsh + 433);

    auto g_0_zzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_shsh + 434);

    auto g_0_zzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_shsh + 435);

    auto g_0_zzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_shsh + 436);

    auto g_0_zzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_shsh + 437);

    auto g_0_zzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_shsh + 438);

    auto g_0_zzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_shsh + 439);

    auto g_0_zzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_shsh + 440);

    /// Set up components of auxilary buffer : SHSI

    auto g_0_xxxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi);

    auto g_0_xxxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 1);

    auto g_0_xxxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 2);

    auto g_0_xxxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 3);

    auto g_0_xxxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 4);

    auto g_0_xxxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 5);

    auto g_0_xxxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 6);

    auto g_0_xxxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 7);

    auto g_0_xxxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 8);

    auto g_0_xxxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 9);

    auto g_0_xxxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 10);

    auto g_0_xxxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 11);

    auto g_0_xxxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 12);

    auto g_0_xxxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 13);

    auto g_0_xxxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 14);

    auto g_0_xxxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 15);

    auto g_0_xxxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 16);

    auto g_0_xxxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 17);

    auto g_0_xxxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 18);

    auto g_0_xxxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 19);

    auto g_0_xxxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 20);

    auto g_0_xxxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 21);

    auto g_0_xxxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 22);

    auto g_0_xxxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 23);

    auto g_0_xxxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 24);

    auto g_0_xxxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 25);

    auto g_0_xxxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 26);

    auto g_0_xxxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 27);

    auto g_0_xxxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 28);

    auto g_0_xxxxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 29);

    auto g_0_xxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 30);

    auto g_0_xxxxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 31);

    auto g_0_xxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 33);

    auto g_0_xxxxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 34);

    auto g_0_xxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 37);

    auto g_0_xxxxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 38);

    auto g_0_xxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 42);

    auto g_0_xxxxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 43);

    auto g_0_xxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 48);

    auto g_0_xxxxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 49);

    auto g_0_xxxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 56);

    auto g_0_xxxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 57);

    auto g_0_xxxxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 58);

    auto g_0_xxxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 59);

    auto g_0_xxxxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 60);

    auto g_0_xxxxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 61);

    auto g_0_xxxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 62);

    auto g_0_xxxxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 63);

    auto g_0_xxxxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 64);

    auto g_0_xxxxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 65);

    auto g_0_xxxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 66);

    auto g_0_xxxxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 67);

    auto g_0_xxxxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 68);

    auto g_0_xxxxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 69);

    auto g_0_xxxxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 70);

    auto g_0_xxxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 71);

    auto g_0_xxxxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 72);

    auto g_0_xxxxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 73);

    auto g_0_xxxxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 74);

    auto g_0_xxxxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 75);

    auto g_0_xxxxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 76);

    auto g_0_xxxxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 78);

    auto g_0_xxxxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 79);

    auto g_0_xxxxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 80);

    auto g_0_xxxxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 81);

    auto g_0_xxxxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 82);

    auto g_0_xxxxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 83);

    auto g_0_xxxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 84);

    auto g_0_xxxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 85);

    auto g_0_xxxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 86);

    auto g_0_xxxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 87);

    auto g_0_xxxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 88);

    auto g_0_xxxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 89);

    auto g_0_xxxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 90);

    auto g_0_xxxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 91);

    auto g_0_xxxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 92);

    auto g_0_xxxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 93);

    auto g_0_xxxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 94);

    auto g_0_xxxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 95);

    auto g_0_xxxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 96);

    auto g_0_xxxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 97);

    auto g_0_xxxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 98);

    auto g_0_xxxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 99);

    auto g_0_xxxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 100);

    auto g_0_xxxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 101);

    auto g_0_xxxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 102);

    auto g_0_xxxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 103);

    auto g_0_xxxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 104);

    auto g_0_xxxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 105);

    auto g_0_xxxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 106);

    auto g_0_xxxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 107);

    auto g_0_xxxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 108);

    auto g_0_xxxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 109);

    auto g_0_xxxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 110);

    auto g_0_xxxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 111);

    auto g_0_xxxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 140);

    auto g_0_xxxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 141);

    auto g_0_xxxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 142);

    auto g_0_xxxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 143);

    auto g_0_xxxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 144);

    auto g_0_xxxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 145);

    auto g_0_xxxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 146);

    auto g_0_xxxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 147);

    auto g_0_xxxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 148);

    auto g_0_xxxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 149);

    auto g_0_xxxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 150);

    auto g_0_xxxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 151);

    auto g_0_xxxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 152);

    auto g_0_xxxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 153);

    auto g_0_xxxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 154);

    auto g_0_xxxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 155);

    auto g_0_xxxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 156);

    auto g_0_xxxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 157);

    auto g_0_xxxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 158);

    auto g_0_xxxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 159);

    auto g_0_xxxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 160);

    auto g_0_xxxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 161);

    auto g_0_xxxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 162);

    auto g_0_xxxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 163);

    auto g_0_xxxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 164);

    auto g_0_xxxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 165);

    auto g_0_xxxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 166);

    auto g_0_xxxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 167);

    auto g_0_xxyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 168);

    auto g_0_xxyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 169);

    auto g_0_xxyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 170);

    auto g_0_xxyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 171);

    auto g_0_xxyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 172);

    auto g_0_xxyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 173);

    auto g_0_xxyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 174);

    auto g_0_xxyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 175);

    auto g_0_xxyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 176);

    auto g_0_xxyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 177);

    auto g_0_xxyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 178);

    auto g_0_xxyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 179);

    auto g_0_xxyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 180);

    auto g_0_xxyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 181);

    auto g_0_xxyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 182);

    auto g_0_xxyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 183);

    auto g_0_xxyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 184);

    auto g_0_xxyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 185);

    auto g_0_xxyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 186);

    auto g_0_xxyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 187);

    auto g_0_xxyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 188);

    auto g_0_xxyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 189);

    auto g_0_xxyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 190);

    auto g_0_xxyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 191);

    auto g_0_xxyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 192);

    auto g_0_xxyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 193);

    auto g_0_xxyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 194);

    auto g_0_xxyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 195);

    auto g_0_xxyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 197);

    auto g_0_xxyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 199);

    auto g_0_xxyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 202);

    auto g_0_xxyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 206);

    auto g_0_xxyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 211);

    auto g_0_xxyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 224);

    auto g_0_xxyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 226);

    auto g_0_xxyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 229);

    auto g_0_xxyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 233);

    auto g_0_xxyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 238);

    auto g_0_xxyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 244);

    auto g_0_xxzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 252);

    auto g_0_xxzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 253);

    auto g_0_xxzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 254);

    auto g_0_xxzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 255);

    auto g_0_xxzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 256);

    auto g_0_xxzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 257);

    auto g_0_xxzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 258);

    auto g_0_xxzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 259);

    auto g_0_xxzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 260);

    auto g_0_xxzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 261);

    auto g_0_xxzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 262);

    auto g_0_xxzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 263);

    auto g_0_xxzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 264);

    auto g_0_xxzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 265);

    auto g_0_xxzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 266);

    auto g_0_xxzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 267);

    auto g_0_xxzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 268);

    auto g_0_xxzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 269);

    auto g_0_xxzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 270);

    auto g_0_xxzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 271);

    auto g_0_xxzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 272);

    auto g_0_xxzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 273);

    auto g_0_xxzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 274);

    auto g_0_xxzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 275);

    auto g_0_xxzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 276);

    auto g_0_xxzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 277);

    auto g_0_xxzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 278);

    auto g_0_xxzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 279);

    auto g_0_xyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 280);

    auto g_0_xyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 281);

    auto g_0_xyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 283);

    auto g_0_xyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 284);

    auto g_0_xyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 286);

    auto g_0_xyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 287);

    auto g_0_xyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 288);

    auto g_0_xyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 290);

    auto g_0_xyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 291);

    auto g_0_xyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 292);

    auto g_0_xyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 293);

    auto g_0_xyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 295);

    auto g_0_xyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 296);

    auto g_0_xyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 297);

    auto g_0_xyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 298);

    auto g_0_xyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 299);

    auto g_0_xyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 301);

    auto g_0_xyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 302);

    auto g_0_xyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 303);

    auto g_0_xyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 304);

    auto g_0_xyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 305);

    auto g_0_xyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 306);

    auto g_0_xyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 307);

    auto g_0_xyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 340);

    auto g_0_xyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 343);

    auto g_0_xyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 344);

    auto g_0_xyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 347);

    auto g_0_xyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 348);

    auto g_0_xyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 349);

    auto g_0_xyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 352);

    auto g_0_xyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 353);

    auto g_0_xyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 354);

    auto g_0_xyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 355);

    auto g_0_xyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 357);

    auto g_0_xyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 358);

    auto g_0_xyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 359);

    auto g_0_xyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 360);

    auto g_0_xyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 361);

    auto g_0_xyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 362);

    auto g_0_xyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 363);

    auto g_0_xzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 392);

    auto g_0_xzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 394);

    auto g_0_xzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 396);

    auto g_0_xzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 397);

    auto g_0_xzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 399);

    auto g_0_xzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 400);

    auto g_0_xzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 401);

    auto g_0_xzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 403);

    auto g_0_xzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 404);

    auto g_0_xzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 405);

    auto g_0_xzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 406);

    auto g_0_xzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 408);

    auto g_0_xzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 409);

    auto g_0_xzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 410);

    auto g_0_xzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 411);

    auto g_0_xzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 412);

    auto g_0_xzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 413);

    auto g_0_xzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 414);

    auto g_0_xzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 415);

    auto g_0_xzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 416);

    auto g_0_xzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 417);

    auto g_0_xzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 418);

    auto g_0_xzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 419);

    auto g_0_yyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 420);

    auto g_0_yyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 421);

    auto g_0_yyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 422);

    auto g_0_yyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 423);

    auto g_0_yyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 424);

    auto g_0_yyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 425);

    auto g_0_yyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 426);

    auto g_0_yyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 427);

    auto g_0_yyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 428);

    auto g_0_yyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 429);

    auto g_0_yyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 430);

    auto g_0_yyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 431);

    auto g_0_yyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 432);

    auto g_0_yyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 433);

    auto g_0_yyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 434);

    auto g_0_yyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 435);

    auto g_0_yyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 436);

    auto g_0_yyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 437);

    auto g_0_yyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 438);

    auto g_0_yyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 439);

    auto g_0_yyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 440);

    auto g_0_yyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 441);

    auto g_0_yyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 442);

    auto g_0_yyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 443);

    auto g_0_yyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 444);

    auto g_0_yyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 445);

    auto g_0_yyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 446);

    auto g_0_yyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 447);

    auto g_0_yyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 449);

    auto g_0_yyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 450);

    auto g_0_yyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 451);

    auto g_0_yyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 452);

    auto g_0_yyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 453);

    auto g_0_yyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 454);

    auto g_0_yyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 455);

    auto g_0_yyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 456);

    auto g_0_yyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 457);

    auto g_0_yyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 458);

    auto g_0_yyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 459);

    auto g_0_yyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 460);

    auto g_0_yyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 461);

    auto g_0_yyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 462);

    auto g_0_yyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 463);

    auto g_0_yyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 464);

    auto g_0_yyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 465);

    auto g_0_yyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 466);

    auto g_0_yyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 467);

    auto g_0_yyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 468);

    auto g_0_yyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 469);

    auto g_0_yyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 470);

    auto g_0_yyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 471);

    auto g_0_yyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 472);

    auto g_0_yyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 473);

    auto g_0_yyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 474);

    auto g_0_yyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 475);

    auto g_0_yyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 476);

    auto g_0_yyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 477);

    auto g_0_yyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 478);

    auto g_0_yyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 479);

    auto g_0_yyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 480);

    auto g_0_yyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 481);

    auto g_0_yyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 482);

    auto g_0_yyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 483);

    auto g_0_yyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 484);

    auto g_0_yyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 485);

    auto g_0_yyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 486);

    auto g_0_yyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 487);

    auto g_0_yyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 488);

    auto g_0_yyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 489);

    auto g_0_yyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 490);

    auto g_0_yyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 491);

    auto g_0_yyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 492);

    auto g_0_yyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 493);

    auto g_0_yyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 494);

    auto g_0_yyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 495);

    auto g_0_yyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 496);

    auto g_0_yyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 497);

    auto g_0_yyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 498);

    auto g_0_yyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 499);

    auto g_0_yyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 500);

    auto g_0_yyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 501);

    auto g_0_yyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 502);

    auto g_0_yyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 503);

    auto g_0_yyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 504);

    auto g_0_yyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 505);

    auto g_0_yyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 506);

    auto g_0_yyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 507);

    auto g_0_yyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 508);

    auto g_0_yyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 509);

    auto g_0_yyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 510);

    auto g_0_yyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 511);

    auto g_0_yyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 512);

    auto g_0_yyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 513);

    auto g_0_yyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 514);

    auto g_0_yyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 515);

    auto g_0_yyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 516);

    auto g_0_yyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 517);

    auto g_0_yyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 518);

    auto g_0_yyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 519);

    auto g_0_yyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 520);

    auto g_0_yyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 521);

    auto g_0_yyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 522);

    auto g_0_yyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 523);

    auto g_0_yyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 524);

    auto g_0_yyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 525);

    auto g_0_yyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 526);

    auto g_0_yyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 527);

    auto g_0_yyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 528);

    auto g_0_yyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 529);

    auto g_0_yyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 530);

    auto g_0_yyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 531);

    auto g_0_yzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 532);

    auto g_0_yzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 533);

    auto g_0_yzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 534);

    auto g_0_yzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 535);

    auto g_0_yzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 536);

    auto g_0_yzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 537);

    auto g_0_yzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 538);

    auto g_0_yzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 539);

    auto g_0_yzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 540);

    auto g_0_yzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 541);

    auto g_0_yzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 542);

    auto g_0_yzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 543);

    auto g_0_yzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 544);

    auto g_0_yzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 545);

    auto g_0_yzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 546);

    auto g_0_yzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 547);

    auto g_0_yzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 548);

    auto g_0_yzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 549);

    auto g_0_yzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 550);

    auto g_0_yzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 551);

    auto g_0_yzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 552);

    auto g_0_yzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 553);

    auto g_0_yzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 554);

    auto g_0_yzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 555);

    auto g_0_yzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 556);

    auto g_0_yzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 557);

    auto g_0_yzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 558);

    auto g_0_yzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 559);

    auto g_0_zzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 560);

    auto g_0_zzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 561);

    auto g_0_zzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 562);

    auto g_0_zzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 563);

    auto g_0_zzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 564);

    auto g_0_zzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 565);

    auto g_0_zzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 566);

    auto g_0_zzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 567);

    auto g_0_zzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 568);

    auto g_0_zzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 569);

    auto g_0_zzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 570);

    auto g_0_zzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 571);

    auto g_0_zzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 572);

    auto g_0_zzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 573);

    auto g_0_zzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 574);

    auto g_0_zzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 575);

    auto g_0_zzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 576);

    auto g_0_zzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 577);

    auto g_0_zzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 578);

    auto g_0_zzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 579);

    auto g_0_zzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 580);

    auto g_0_zzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 581);

    auto g_0_zzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 582);

    auto g_0_zzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 583);

    auto g_0_zzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 584);

    auto g_0_zzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 585);

    auto g_0_zzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 586);

    auto g_0_zzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 587);

    /// Set up components of auxilary buffer : SHSI

    auto g_0_xxxxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi);

    auto g_0_xxxxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 1);

    auto g_0_xxxxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 2);

    auto g_0_xxxxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 3);

    auto g_0_xxxxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 4);

    auto g_0_xxxxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 5);

    auto g_0_xxxxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 6);

    auto g_0_xxxxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 7);

    auto g_0_xxxxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 8);

    auto g_0_xxxxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 9);

    auto g_0_xxxxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 10);

    auto g_0_xxxxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 11);

    auto g_0_xxxxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 12);

    auto g_0_xxxxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 13);

    auto g_0_xxxxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 14);

    auto g_0_xxxxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 15);

    auto g_0_xxxxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 16);

    auto g_0_xxxxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 17);

    auto g_0_xxxxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 18);

    auto g_0_xxxxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 19);

    auto g_0_xxxxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 20);

    auto g_0_xxxxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 21);

    auto g_0_xxxxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 22);

    auto g_0_xxxxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 23);

    auto g_0_xxxxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 24);

    auto g_0_xxxxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 25);

    auto g_0_xxxxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 26);

    auto g_0_xxxxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 27);

    auto g_0_xxxxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 28);

    auto g_0_xxxxy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 29);

    auto g_0_xxxxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 30);

    auto g_0_xxxxy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 31);

    auto g_0_xxxxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 33);

    auto g_0_xxxxy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 34);

    auto g_0_xxxxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 37);

    auto g_0_xxxxy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 38);

    auto g_0_xxxxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 42);

    auto g_0_xxxxy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 43);

    auto g_0_xxxxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 48);

    auto g_0_xxxxy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 49);

    auto g_0_xxxxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 56);

    auto g_0_xxxxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 57);

    auto g_0_xxxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 58);

    auto g_0_xxxxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 59);

    auto g_0_xxxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 60);

    auto g_0_xxxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 61);

    auto g_0_xxxxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 62);

    auto g_0_xxxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 63);

    auto g_0_xxxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 64);

    auto g_0_xxxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 65);

    auto g_0_xxxxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 66);

    auto g_0_xxxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 67);

    auto g_0_xxxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 68);

    auto g_0_xxxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 69);

    auto g_0_xxxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 70);

    auto g_0_xxxxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 71);

    auto g_0_xxxxz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 72);

    auto g_0_xxxxz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 73);

    auto g_0_xxxxz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 74);

    auto g_0_xxxxz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 75);

    auto g_0_xxxxz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 76);

    auto g_0_xxxxz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 78);

    auto g_0_xxxxz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 79);

    auto g_0_xxxxz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 80);

    auto g_0_xxxxz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 81);

    auto g_0_xxxxz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 82);

    auto g_0_xxxxz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 83);

    auto g_0_xxxyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 84);

    auto g_0_xxxyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 85);

    auto g_0_xxxyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 86);

    auto g_0_xxxyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 87);

    auto g_0_xxxyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 88);

    auto g_0_xxxyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 89);

    auto g_0_xxxyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 90);

    auto g_0_xxxyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 91);

    auto g_0_xxxyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 92);

    auto g_0_xxxyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 93);

    auto g_0_xxxyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 94);

    auto g_0_xxxyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 95);

    auto g_0_xxxyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 96);

    auto g_0_xxxyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 97);

    auto g_0_xxxyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 98);

    auto g_0_xxxyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 99);

    auto g_0_xxxyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 100);

    auto g_0_xxxyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 101);

    auto g_0_xxxyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 102);

    auto g_0_xxxyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 103);

    auto g_0_xxxyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 104);

    auto g_0_xxxyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 105);

    auto g_0_xxxyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 106);

    auto g_0_xxxyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 107);

    auto g_0_xxxyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 108);

    auto g_0_xxxyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 109);

    auto g_0_xxxyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 110);

    auto g_0_xxxyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 111);

    auto g_0_xxxzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 140);

    auto g_0_xxxzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 141);

    auto g_0_xxxzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 142);

    auto g_0_xxxzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 143);

    auto g_0_xxxzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 144);

    auto g_0_xxxzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 145);

    auto g_0_xxxzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 146);

    auto g_0_xxxzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 147);

    auto g_0_xxxzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 148);

    auto g_0_xxxzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 149);

    auto g_0_xxxzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 150);

    auto g_0_xxxzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 151);

    auto g_0_xxxzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 152);

    auto g_0_xxxzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 153);

    auto g_0_xxxzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 154);

    auto g_0_xxxzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 155);

    auto g_0_xxxzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 156);

    auto g_0_xxxzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 157);

    auto g_0_xxxzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 158);

    auto g_0_xxxzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 159);

    auto g_0_xxxzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 160);

    auto g_0_xxxzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 161);

    auto g_0_xxxzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 162);

    auto g_0_xxxzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 163);

    auto g_0_xxxzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 164);

    auto g_0_xxxzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 165);

    auto g_0_xxxzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 166);

    auto g_0_xxxzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 167);

    auto g_0_xxyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 168);

    auto g_0_xxyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 169);

    auto g_0_xxyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 170);

    auto g_0_xxyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 171);

    auto g_0_xxyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 172);

    auto g_0_xxyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 173);

    auto g_0_xxyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 174);

    auto g_0_xxyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 175);

    auto g_0_xxyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 176);

    auto g_0_xxyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 177);

    auto g_0_xxyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 178);

    auto g_0_xxyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 179);

    auto g_0_xxyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 180);

    auto g_0_xxyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 181);

    auto g_0_xxyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 182);

    auto g_0_xxyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 183);

    auto g_0_xxyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 184);

    auto g_0_xxyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 185);

    auto g_0_xxyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 186);

    auto g_0_xxyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 187);

    auto g_0_xxyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 188);

    auto g_0_xxyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 189);

    auto g_0_xxyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 190);

    auto g_0_xxyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 191);

    auto g_0_xxyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 192);

    auto g_0_xxyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 193);

    auto g_0_xxyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 194);

    auto g_0_xxyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 195);

    auto g_0_xxyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 197);

    auto g_0_xxyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 199);

    auto g_0_xxyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 202);

    auto g_0_xxyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 206);

    auto g_0_xxyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 211);

    auto g_0_xxyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 224);

    auto g_0_xxyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 226);

    auto g_0_xxyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 229);

    auto g_0_xxyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 233);

    auto g_0_xxyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 238);

    auto g_0_xxyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 244);

    auto g_0_xxzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 252);

    auto g_0_xxzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 253);

    auto g_0_xxzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 254);

    auto g_0_xxzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 255);

    auto g_0_xxzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 256);

    auto g_0_xxzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 257);

    auto g_0_xxzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 258);

    auto g_0_xxzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 259);

    auto g_0_xxzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 260);

    auto g_0_xxzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 261);

    auto g_0_xxzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 262);

    auto g_0_xxzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 263);

    auto g_0_xxzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 264);

    auto g_0_xxzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 265);

    auto g_0_xxzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 266);

    auto g_0_xxzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 267);

    auto g_0_xxzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 268);

    auto g_0_xxzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 269);

    auto g_0_xxzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 270);

    auto g_0_xxzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 271);

    auto g_0_xxzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 272);

    auto g_0_xxzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 273);

    auto g_0_xxzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 274);

    auto g_0_xxzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 275);

    auto g_0_xxzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 276);

    auto g_0_xxzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 277);

    auto g_0_xxzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 278);

    auto g_0_xxzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 279);

    auto g_0_xyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 280);

    auto g_0_xyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 281);

    auto g_0_xyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 283);

    auto g_0_xyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 284);

    auto g_0_xyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 286);

    auto g_0_xyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 287);

    auto g_0_xyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 288);

    auto g_0_xyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 290);

    auto g_0_xyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 291);

    auto g_0_xyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 292);

    auto g_0_xyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 293);

    auto g_0_xyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 295);

    auto g_0_xyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 296);

    auto g_0_xyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 297);

    auto g_0_xyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 298);

    auto g_0_xyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 299);

    auto g_0_xyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 301);

    auto g_0_xyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 302);

    auto g_0_xyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 303);

    auto g_0_xyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 304);

    auto g_0_xyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 305);

    auto g_0_xyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 306);

    auto g_0_xyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 307);

    auto g_0_xyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 340);

    auto g_0_xyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 343);

    auto g_0_xyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 344);

    auto g_0_xyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 347);

    auto g_0_xyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 348);

    auto g_0_xyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 349);

    auto g_0_xyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 352);

    auto g_0_xyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 353);

    auto g_0_xyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 354);

    auto g_0_xyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 355);

    auto g_0_xyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 357);

    auto g_0_xyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 358);

    auto g_0_xyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 359);

    auto g_0_xyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 360);

    auto g_0_xyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 361);

    auto g_0_xyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 362);

    auto g_0_xyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 363);

    auto g_0_xzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 392);

    auto g_0_xzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 394);

    auto g_0_xzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 396);

    auto g_0_xzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 397);

    auto g_0_xzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 399);

    auto g_0_xzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 400);

    auto g_0_xzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 401);

    auto g_0_xzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 403);

    auto g_0_xzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 404);

    auto g_0_xzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 405);

    auto g_0_xzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 406);

    auto g_0_xzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 408);

    auto g_0_xzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 409);

    auto g_0_xzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 410);

    auto g_0_xzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 411);

    auto g_0_xzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 412);

    auto g_0_xzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 413);

    auto g_0_xzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 414);

    auto g_0_xzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 415);

    auto g_0_xzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 416);

    auto g_0_xzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 417);

    auto g_0_xzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 418);

    auto g_0_xzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 419);

    auto g_0_yyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 420);

    auto g_0_yyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 421);

    auto g_0_yyyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 422);

    auto g_0_yyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 423);

    auto g_0_yyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 424);

    auto g_0_yyyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 425);

    auto g_0_yyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 426);

    auto g_0_yyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 427);

    auto g_0_yyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 428);

    auto g_0_yyyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 429);

    auto g_0_yyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 430);

    auto g_0_yyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 431);

    auto g_0_yyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 432);

    auto g_0_yyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 433);

    auto g_0_yyyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 434);

    auto g_0_yyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 435);

    auto g_0_yyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 436);

    auto g_0_yyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 437);

    auto g_0_yyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 438);

    auto g_0_yyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 439);

    auto g_0_yyyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 440);

    auto g_0_yyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 441);

    auto g_0_yyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 442);

    auto g_0_yyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 443);

    auto g_0_yyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 444);

    auto g_0_yyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 445);

    auto g_0_yyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 446);

    auto g_0_yyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 447);

    auto g_0_yyyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 449);

    auto g_0_yyyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 450);

    auto g_0_yyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 451);

    auto g_0_yyyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 452);

    auto g_0_yyyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 453);

    auto g_0_yyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 454);

    auto g_0_yyyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 455);

    auto g_0_yyyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 456);

    auto g_0_yyyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 457);

    auto g_0_yyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 458);

    auto g_0_yyyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 459);

    auto g_0_yyyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 460);

    auto g_0_yyyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 461);

    auto g_0_yyyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 462);

    auto g_0_yyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 463);

    auto g_0_yyyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 464);

    auto g_0_yyyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 465);

    auto g_0_yyyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 466);

    auto g_0_yyyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 467);

    auto g_0_yyyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 468);

    auto g_0_yyyyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 469);

    auto g_0_yyyyz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 470);

    auto g_0_yyyyz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 471);

    auto g_0_yyyyz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 472);

    auto g_0_yyyyz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 473);

    auto g_0_yyyyz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 474);

    auto g_0_yyyyz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 475);

    auto g_0_yyyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 476);

    auto g_0_yyyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 477);

    auto g_0_yyyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 478);

    auto g_0_yyyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 479);

    auto g_0_yyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 480);

    auto g_0_yyyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 481);

    auto g_0_yyyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 482);

    auto g_0_yyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 483);

    auto g_0_yyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 484);

    auto g_0_yyyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 485);

    auto g_0_yyyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 486);

    auto g_0_yyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 487);

    auto g_0_yyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 488);

    auto g_0_yyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 489);

    auto g_0_yyyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 490);

    auto g_0_yyyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 491);

    auto g_0_yyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 492);

    auto g_0_yyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 493);

    auto g_0_yyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 494);

    auto g_0_yyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 495);

    auto g_0_yyyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 496);

    auto g_0_yyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 497);

    auto g_0_yyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 498);

    auto g_0_yyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 499);

    auto g_0_yyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 500);

    auto g_0_yyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 501);

    auto g_0_yyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 502);

    auto g_0_yyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 503);

    auto g_0_yyzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 504);

    auto g_0_yyzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 505);

    auto g_0_yyzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 506);

    auto g_0_yyzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 507);

    auto g_0_yyzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 508);

    auto g_0_yyzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 509);

    auto g_0_yyzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 510);

    auto g_0_yyzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 511);

    auto g_0_yyzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 512);

    auto g_0_yyzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 513);

    auto g_0_yyzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 514);

    auto g_0_yyzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 515);

    auto g_0_yyzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 516);

    auto g_0_yyzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 517);

    auto g_0_yyzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 518);

    auto g_0_yyzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 519);

    auto g_0_yyzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 520);

    auto g_0_yyzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 521);

    auto g_0_yyzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 522);

    auto g_0_yyzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 523);

    auto g_0_yyzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 524);

    auto g_0_yyzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 525);

    auto g_0_yyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 526);

    auto g_0_yyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 527);

    auto g_0_yyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 528);

    auto g_0_yyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 529);

    auto g_0_yyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 530);

    auto g_0_yyzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 531);

    auto g_0_yzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 532);

    auto g_0_yzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 533);

    auto g_0_yzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 534);

    auto g_0_yzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 535);

    auto g_0_yzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 536);

    auto g_0_yzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 537);

    auto g_0_yzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 538);

    auto g_0_yzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 539);

    auto g_0_yzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 540);

    auto g_0_yzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 541);

    auto g_0_yzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 542);

    auto g_0_yzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 543);

    auto g_0_yzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 544);

    auto g_0_yzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 545);

    auto g_0_yzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 546);

    auto g_0_yzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 547);

    auto g_0_yzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 548);

    auto g_0_yzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 549);

    auto g_0_yzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 550);

    auto g_0_yzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 551);

    auto g_0_yzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 552);

    auto g_0_yzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 553);

    auto g_0_yzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 554);

    auto g_0_yzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 555);

    auto g_0_yzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 556);

    auto g_0_yzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 557);

    auto g_0_yzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 558);

    auto g_0_yzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 559);

    auto g_0_zzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 560);

    auto g_0_zzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 561);

    auto g_0_zzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 562);

    auto g_0_zzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 563);

    auto g_0_zzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 564);

    auto g_0_zzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 565);

    auto g_0_zzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 566);

    auto g_0_zzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 567);

    auto g_0_zzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 568);

    auto g_0_zzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 569);

    auto g_0_zzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 570);

    auto g_0_zzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 571);

    auto g_0_zzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 572);

    auto g_0_zzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 573);

    auto g_0_zzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 574);

    auto g_0_zzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 575);

    auto g_0_zzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 576);

    auto g_0_zzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 577);

    auto g_0_zzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 578);

    auto g_0_zzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 579);

    auto g_0_zzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 580);

    auto g_0_zzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 581);

    auto g_0_zzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 582);

    auto g_0_zzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 583);

    auto g_0_zzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 584);

    auto g_0_zzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 585);

    auto g_0_zzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 586);

    auto g_0_zzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 587);

    /// Set up 0-28 components of targeted buffer : SISI

    auto g_0_xxxxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi);

    auto g_0_xxxxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 1);

    auto g_0_xxxxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 2);

    auto g_0_xxxxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 3);

    auto g_0_xxxxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 4);

    auto g_0_xxxxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 5);

    auto g_0_xxxxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 6);

    auto g_0_xxxxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 7);

    auto g_0_xxxxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 8);

    auto g_0_xxxxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 9);

    auto g_0_xxxxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 10);

    auto g_0_xxxxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 11);

    auto g_0_xxxxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 12);

    auto g_0_xxxxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 13);

    auto g_0_xxxxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 14);

    auto g_0_xxxxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 15);

    auto g_0_xxxxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 16);

    auto g_0_xxxxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 17);

    auto g_0_xxxxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 18);

    auto g_0_xxxxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 19);

    auto g_0_xxxxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 20);

    auto g_0_xxxxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 21);

    auto g_0_xxxxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 22);

    auto g_0_xxxxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 23);

    auto g_0_xxxxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 24);

    auto g_0_xxxxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 25);

    auto g_0_xxxxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 26);

    auto g_0_xxxxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 27);

#pragma omp simd aligned(g_0_xxxx_0_xxxxxx_0,       \
                             g_0_xxxx_0_xxxxxx_1,   \
                             g_0_xxxx_0_xxxxxy_0,   \
                             g_0_xxxx_0_xxxxxy_1,   \
                             g_0_xxxx_0_xxxxxz_0,   \
                             g_0_xxxx_0_xxxxxz_1,   \
                             g_0_xxxx_0_xxxxyy_0,   \
                             g_0_xxxx_0_xxxxyy_1,   \
                             g_0_xxxx_0_xxxxyz_0,   \
                             g_0_xxxx_0_xxxxyz_1,   \
                             g_0_xxxx_0_xxxxzz_0,   \
                             g_0_xxxx_0_xxxxzz_1,   \
                             g_0_xxxx_0_xxxyyy_0,   \
                             g_0_xxxx_0_xxxyyy_1,   \
                             g_0_xxxx_0_xxxyyz_0,   \
                             g_0_xxxx_0_xxxyyz_1,   \
                             g_0_xxxx_0_xxxyzz_0,   \
                             g_0_xxxx_0_xxxyzz_1,   \
                             g_0_xxxx_0_xxxzzz_0,   \
                             g_0_xxxx_0_xxxzzz_1,   \
                             g_0_xxxx_0_xxyyyy_0,   \
                             g_0_xxxx_0_xxyyyy_1,   \
                             g_0_xxxx_0_xxyyyz_0,   \
                             g_0_xxxx_0_xxyyyz_1,   \
                             g_0_xxxx_0_xxyyzz_0,   \
                             g_0_xxxx_0_xxyyzz_1,   \
                             g_0_xxxx_0_xxyzzz_0,   \
                             g_0_xxxx_0_xxyzzz_1,   \
                             g_0_xxxx_0_xxzzzz_0,   \
                             g_0_xxxx_0_xxzzzz_1,   \
                             g_0_xxxx_0_xyyyyy_0,   \
                             g_0_xxxx_0_xyyyyy_1,   \
                             g_0_xxxx_0_xyyyyz_0,   \
                             g_0_xxxx_0_xyyyyz_1,   \
                             g_0_xxxx_0_xyyyzz_0,   \
                             g_0_xxxx_0_xyyyzz_1,   \
                             g_0_xxxx_0_xyyzzz_0,   \
                             g_0_xxxx_0_xyyzzz_1,   \
                             g_0_xxxx_0_xyzzzz_0,   \
                             g_0_xxxx_0_xyzzzz_1,   \
                             g_0_xxxx_0_xzzzzz_0,   \
                             g_0_xxxx_0_xzzzzz_1,   \
                             g_0_xxxx_0_yyyyyy_0,   \
                             g_0_xxxx_0_yyyyyy_1,   \
                             g_0_xxxx_0_yyyyyz_0,   \
                             g_0_xxxx_0_yyyyyz_1,   \
                             g_0_xxxx_0_yyyyzz_0,   \
                             g_0_xxxx_0_yyyyzz_1,   \
                             g_0_xxxx_0_yyyzzz_0,   \
                             g_0_xxxx_0_yyyzzz_1,   \
                             g_0_xxxx_0_yyzzzz_0,   \
                             g_0_xxxx_0_yyzzzz_1,   \
                             g_0_xxxx_0_yzzzzz_0,   \
                             g_0_xxxx_0_yzzzzz_1,   \
                             g_0_xxxx_0_zzzzzz_0,   \
                             g_0_xxxx_0_zzzzzz_1,   \
                             g_0_xxxxx_0_xxxxx_1,   \
                             g_0_xxxxx_0_xxxxxx_0,  \
                             g_0_xxxxx_0_xxxxxx_1,  \
                             g_0_xxxxx_0_xxxxxy_0,  \
                             g_0_xxxxx_0_xxxxxy_1,  \
                             g_0_xxxxx_0_xxxxxz_0,  \
                             g_0_xxxxx_0_xxxxxz_1,  \
                             g_0_xxxxx_0_xxxxy_1,   \
                             g_0_xxxxx_0_xxxxyy_0,  \
                             g_0_xxxxx_0_xxxxyy_1,  \
                             g_0_xxxxx_0_xxxxyz_0,  \
                             g_0_xxxxx_0_xxxxyz_1,  \
                             g_0_xxxxx_0_xxxxz_1,   \
                             g_0_xxxxx_0_xxxxzz_0,  \
                             g_0_xxxxx_0_xxxxzz_1,  \
                             g_0_xxxxx_0_xxxyy_1,   \
                             g_0_xxxxx_0_xxxyyy_0,  \
                             g_0_xxxxx_0_xxxyyy_1,  \
                             g_0_xxxxx_0_xxxyyz_0,  \
                             g_0_xxxxx_0_xxxyyz_1,  \
                             g_0_xxxxx_0_xxxyz_1,   \
                             g_0_xxxxx_0_xxxyzz_0,  \
                             g_0_xxxxx_0_xxxyzz_1,  \
                             g_0_xxxxx_0_xxxzz_1,   \
                             g_0_xxxxx_0_xxxzzz_0,  \
                             g_0_xxxxx_0_xxxzzz_1,  \
                             g_0_xxxxx_0_xxyyy_1,   \
                             g_0_xxxxx_0_xxyyyy_0,  \
                             g_0_xxxxx_0_xxyyyy_1,  \
                             g_0_xxxxx_0_xxyyyz_0,  \
                             g_0_xxxxx_0_xxyyyz_1,  \
                             g_0_xxxxx_0_xxyyz_1,   \
                             g_0_xxxxx_0_xxyyzz_0,  \
                             g_0_xxxxx_0_xxyyzz_1,  \
                             g_0_xxxxx_0_xxyzz_1,   \
                             g_0_xxxxx_0_xxyzzz_0,  \
                             g_0_xxxxx_0_xxyzzz_1,  \
                             g_0_xxxxx_0_xxzzz_1,   \
                             g_0_xxxxx_0_xxzzzz_0,  \
                             g_0_xxxxx_0_xxzzzz_1,  \
                             g_0_xxxxx_0_xyyyy_1,   \
                             g_0_xxxxx_0_xyyyyy_0,  \
                             g_0_xxxxx_0_xyyyyy_1,  \
                             g_0_xxxxx_0_xyyyyz_0,  \
                             g_0_xxxxx_0_xyyyyz_1,  \
                             g_0_xxxxx_0_xyyyz_1,   \
                             g_0_xxxxx_0_xyyyzz_0,  \
                             g_0_xxxxx_0_xyyyzz_1,  \
                             g_0_xxxxx_0_xyyzz_1,   \
                             g_0_xxxxx_0_xyyzzz_0,  \
                             g_0_xxxxx_0_xyyzzz_1,  \
                             g_0_xxxxx_0_xyzzz_1,   \
                             g_0_xxxxx_0_xyzzzz_0,  \
                             g_0_xxxxx_0_xyzzzz_1,  \
                             g_0_xxxxx_0_xzzzz_1,   \
                             g_0_xxxxx_0_xzzzzz_0,  \
                             g_0_xxxxx_0_xzzzzz_1,  \
                             g_0_xxxxx_0_yyyyy_1,   \
                             g_0_xxxxx_0_yyyyyy_0,  \
                             g_0_xxxxx_0_yyyyyy_1,  \
                             g_0_xxxxx_0_yyyyyz_0,  \
                             g_0_xxxxx_0_yyyyyz_1,  \
                             g_0_xxxxx_0_yyyyz_1,   \
                             g_0_xxxxx_0_yyyyzz_0,  \
                             g_0_xxxxx_0_yyyyzz_1,  \
                             g_0_xxxxx_0_yyyzz_1,   \
                             g_0_xxxxx_0_yyyzzz_0,  \
                             g_0_xxxxx_0_yyyzzz_1,  \
                             g_0_xxxxx_0_yyzzz_1,   \
                             g_0_xxxxx_0_yyzzzz_0,  \
                             g_0_xxxxx_0_yyzzzz_1,  \
                             g_0_xxxxx_0_yzzzz_1,   \
                             g_0_xxxxx_0_yzzzzz_0,  \
                             g_0_xxxxx_0_yzzzzz_1,  \
                             g_0_xxxxx_0_zzzzz_1,   \
                             g_0_xxxxx_0_zzzzzz_0,  \
                             g_0_xxxxx_0_zzzzzz_1,  \
                             g_0_xxxxxx_0_xxxxxx_0, \
                             g_0_xxxxxx_0_xxxxxy_0, \
                             g_0_xxxxxx_0_xxxxxz_0, \
                             g_0_xxxxxx_0_xxxxyy_0, \
                             g_0_xxxxxx_0_xxxxyz_0, \
                             g_0_xxxxxx_0_xxxxzz_0, \
                             g_0_xxxxxx_0_xxxyyy_0, \
                             g_0_xxxxxx_0_xxxyyz_0, \
                             g_0_xxxxxx_0_xxxyzz_0, \
                             g_0_xxxxxx_0_xxxzzz_0, \
                             g_0_xxxxxx_0_xxyyyy_0, \
                             g_0_xxxxxx_0_xxyyyz_0, \
                             g_0_xxxxxx_0_xxyyzz_0, \
                             g_0_xxxxxx_0_xxyzzz_0, \
                             g_0_xxxxxx_0_xxzzzz_0, \
                             g_0_xxxxxx_0_xyyyyy_0, \
                             g_0_xxxxxx_0_xyyyyz_0, \
                             g_0_xxxxxx_0_xyyyzz_0, \
                             g_0_xxxxxx_0_xyyzzz_0, \
                             g_0_xxxxxx_0_xyzzzz_0, \
                             g_0_xxxxxx_0_xzzzzz_0, \
                             g_0_xxxxxx_0_yyyyyy_0, \
                             g_0_xxxxxx_0_yyyyyz_0, \
                             g_0_xxxxxx_0_yyyyzz_0, \
                             g_0_xxxxxx_0_yyyzzz_0, \
                             g_0_xxxxxx_0_yyzzzz_0, \
                             g_0_xxxxxx_0_yzzzzz_0, \
                             g_0_xxxxxx_0_zzzzzz_0, \
                             wp_x,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxx_0_xxxxxx_0[i] = 5.0 * g_0_xxxx_0_xxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxx_1[i] * fti_ab_0 +
                                   6.0 * g_0_xxxxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxx_0[i] * pb_x + g_0_xxxxx_0_xxxxxx_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxy_0[i] = 5.0 * g_0_xxxx_0_xxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxy_1[i] * fti_ab_0 +
                                   5.0 * g_0_xxxxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxy_0[i] * pb_x + g_0_xxxxx_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxz_0[i] = 5.0 * g_0_xxxx_0_xxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxz_1[i] * fti_ab_0 +
                                   5.0 * g_0_xxxxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxz_0[i] * pb_x + g_0_xxxxx_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxyy_0[i] = 5.0 * g_0_xxxx_0_xxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxyy_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxxxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyy_0[i] * pb_x + g_0_xxxxx_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxyz_0[i] = 5.0 * g_0_xxxx_0_xxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxyz_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxxxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyz_0[i] * pb_x + g_0_xxxxx_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxzz_0[i] = 5.0 * g_0_xxxx_0_xxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxxxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxzz_0[i] * pb_x + g_0_xxxxx_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyyy_0[i] = 5.0 * g_0_xxxx_0_xxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyyy_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxxxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyy_0[i] * pb_x + g_0_xxxxx_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyyz_0[i] = 5.0 * g_0_xxxx_0_xxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyyz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxxxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyz_0[i] * pb_x + g_0_xxxxx_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyzz_0[i] = 5.0 * g_0_xxxx_0_xxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxxxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyzz_0[i] * pb_x + g_0_xxxxx_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxzzz_0[i] = 5.0 * g_0_xxxx_0_xxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxxxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxzzz_0[i] * pb_x + g_0_xxxxx_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyyy_0[i] = 5.0 * g_0_xxxx_0_xxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyyy_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyy_0[i] * pb_x + g_0_xxxxx_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyyz_0[i] = 5.0 * g_0_xxxx_0_xxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyyz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyz_0[i] * pb_x + g_0_xxxxx_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyzz_0[i] = 5.0 * g_0_xxxx_0_xxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyzz_0[i] * pb_x + g_0_xxxxx_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyzzz_0[i] = 5.0 * g_0_xxxx_0_xxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyzzz_0[i] * pb_x + g_0_xxxxx_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxzzzz_0[i] = 5.0 * g_0_xxxx_0_xxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxzzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxzzzz_0[i] * pb_x + g_0_xxxxx_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyyy_0[i] = 5.0 * g_0_xxxx_0_xyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyyy_1[i] * fti_ab_0 +
                                   g_0_xxxxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyy_0[i] * pb_x + g_0_xxxxx_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyyz_0[i] = 5.0 * g_0_xxxx_0_xyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyyz_1[i] * fti_ab_0 +
                                   g_0_xxxxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyz_0[i] * pb_x + g_0_xxxxx_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyzz_0[i] = 5.0 * g_0_xxxx_0_xyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyzz_1[i] * fti_ab_0 +
                                   g_0_xxxxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyzz_0[i] * pb_x + g_0_xxxxx_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyzzz_0[i] = 5.0 * g_0_xxxx_0_xyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyzzz_1[i] * fti_ab_0 +
                                   g_0_xxxxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyzzz_0[i] * pb_x + g_0_xxxxx_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyzzzz_0[i] = 5.0 * g_0_xxxx_0_xyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyzzzz_1[i] * fti_ab_0 +
                                   g_0_xxxxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzzzz_0[i] * pb_x + g_0_xxxxx_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xzzzzz_0[i] = 5.0 * g_0_xxxx_0_xzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xzzzzz_1[i] * fti_ab_0 +
                                   g_0_xxxxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xzzzzz_0[i] * pb_x + g_0_xxxxx_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyyy_0[i] = 5.0 * g_0_xxxx_0_yyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyyy_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyyyy_0[i] * pb_x +
                                   g_0_xxxxx_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyyz_0[i] = 5.0 * g_0_xxxx_0_yyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyyz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyyyz_0[i] * pb_x +
                                   g_0_xxxxx_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyzz_0[i] = 5.0 * g_0_xxxx_0_yyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyyzz_0[i] * pb_x +
                                   g_0_xxxxx_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyzzz_0[i] = 5.0 * g_0_xxxx_0_yyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyzzz_0[i] * pb_x +
                                   g_0_xxxxx_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyzzzz_0[i] = 5.0 * g_0_xxxx_0_yyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyzzzz_0[i] * pb_x +
                                   g_0_xxxxx_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yzzzzz_0[i] = 5.0 * g_0_xxxx_0_yzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yzzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yzzzzz_0[i] * pb_x +
                                   g_0_xxxxx_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_zzzzzz_0[i] = 5.0 * g_0_xxxx_0_zzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_zzzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_zzzzzz_0[i] * pb_x +
                                   g_0_xxxxx_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 28-56 components of targeted buffer : SISI

    auto g_0_xxxxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 28);

    auto g_0_xxxxxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 29);

    auto g_0_xxxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 30);

    auto g_0_xxxxxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 31);

    auto g_0_xxxxxy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 32);

    auto g_0_xxxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 33);

    auto g_0_xxxxxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 34);

    auto g_0_xxxxxy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 35);

    auto g_0_xxxxxy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 36);

    auto g_0_xxxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 37);

    auto g_0_xxxxxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 38);

    auto g_0_xxxxxy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 39);

    auto g_0_xxxxxy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 40);

    auto g_0_xxxxxy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 41);

    auto g_0_xxxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 42);

    auto g_0_xxxxxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 43);

    auto g_0_xxxxxy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 44);

    auto g_0_xxxxxy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 45);

    auto g_0_xxxxxy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 46);

    auto g_0_xxxxxy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 47);

    auto g_0_xxxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 48);

    auto g_0_xxxxxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 49);

    auto g_0_xxxxxy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 50);

    auto g_0_xxxxxy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 51);

    auto g_0_xxxxxy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 52);

    auto g_0_xxxxxy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 53);

    auto g_0_xxxxxy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 54);

    auto g_0_xxxxxy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 55);

#pragma omp simd aligned(g_0_xxxxx_0_xxxxx_1,       \
                             g_0_xxxxx_0_xxxxxx_0,  \
                             g_0_xxxxx_0_xxxxxx_1,  \
                             g_0_xxxxx_0_xxxxxy_0,  \
                             g_0_xxxxx_0_xxxxxy_1,  \
                             g_0_xxxxx_0_xxxxxz_0,  \
                             g_0_xxxxx_0_xxxxxz_1,  \
                             g_0_xxxxx_0_xxxxy_1,   \
                             g_0_xxxxx_0_xxxxyy_0,  \
                             g_0_xxxxx_0_xxxxyy_1,  \
                             g_0_xxxxx_0_xxxxyz_0,  \
                             g_0_xxxxx_0_xxxxyz_1,  \
                             g_0_xxxxx_0_xxxxz_1,   \
                             g_0_xxxxx_0_xxxxzz_0,  \
                             g_0_xxxxx_0_xxxxzz_1,  \
                             g_0_xxxxx_0_xxxyy_1,   \
                             g_0_xxxxx_0_xxxyyy_0,  \
                             g_0_xxxxx_0_xxxyyy_1,  \
                             g_0_xxxxx_0_xxxyyz_0,  \
                             g_0_xxxxx_0_xxxyyz_1,  \
                             g_0_xxxxx_0_xxxyz_1,   \
                             g_0_xxxxx_0_xxxyzz_0,  \
                             g_0_xxxxx_0_xxxyzz_1,  \
                             g_0_xxxxx_0_xxxzz_1,   \
                             g_0_xxxxx_0_xxxzzz_0,  \
                             g_0_xxxxx_0_xxxzzz_1,  \
                             g_0_xxxxx_0_xxyyy_1,   \
                             g_0_xxxxx_0_xxyyyy_0,  \
                             g_0_xxxxx_0_xxyyyy_1,  \
                             g_0_xxxxx_0_xxyyyz_0,  \
                             g_0_xxxxx_0_xxyyyz_1,  \
                             g_0_xxxxx_0_xxyyz_1,   \
                             g_0_xxxxx_0_xxyyzz_0,  \
                             g_0_xxxxx_0_xxyyzz_1,  \
                             g_0_xxxxx_0_xxyzz_1,   \
                             g_0_xxxxx_0_xxyzzz_0,  \
                             g_0_xxxxx_0_xxyzzz_1,  \
                             g_0_xxxxx_0_xxzzz_1,   \
                             g_0_xxxxx_0_xxzzzz_0,  \
                             g_0_xxxxx_0_xxzzzz_1,  \
                             g_0_xxxxx_0_xyyyy_1,   \
                             g_0_xxxxx_0_xyyyyy_0,  \
                             g_0_xxxxx_0_xyyyyy_1,  \
                             g_0_xxxxx_0_xyyyyz_0,  \
                             g_0_xxxxx_0_xyyyyz_1,  \
                             g_0_xxxxx_0_xyyyz_1,   \
                             g_0_xxxxx_0_xyyyzz_0,  \
                             g_0_xxxxx_0_xyyyzz_1,  \
                             g_0_xxxxx_0_xyyzz_1,   \
                             g_0_xxxxx_0_xyyzzz_0,  \
                             g_0_xxxxx_0_xyyzzz_1,  \
                             g_0_xxxxx_0_xyzzz_1,   \
                             g_0_xxxxx_0_xyzzzz_0,  \
                             g_0_xxxxx_0_xyzzzz_1,  \
                             g_0_xxxxx_0_xzzzz_1,   \
                             g_0_xxxxx_0_xzzzzz_0,  \
                             g_0_xxxxx_0_xzzzzz_1,  \
                             g_0_xxxxx_0_yyyyy_1,   \
                             g_0_xxxxx_0_yyyyyy_0,  \
                             g_0_xxxxx_0_yyyyyy_1,  \
                             g_0_xxxxx_0_yyyyyz_0,  \
                             g_0_xxxxx_0_yyyyyz_1,  \
                             g_0_xxxxx_0_yyyyz_1,   \
                             g_0_xxxxx_0_yyyyzz_0,  \
                             g_0_xxxxx_0_yyyyzz_1,  \
                             g_0_xxxxx_0_yyyzz_1,   \
                             g_0_xxxxx_0_yyyzzz_0,  \
                             g_0_xxxxx_0_yyyzzz_1,  \
                             g_0_xxxxx_0_yyzzz_1,   \
                             g_0_xxxxx_0_yyzzzz_0,  \
                             g_0_xxxxx_0_yyzzzz_1,  \
                             g_0_xxxxx_0_yzzzz_1,   \
                             g_0_xxxxx_0_yzzzzz_0,  \
                             g_0_xxxxx_0_yzzzzz_1,  \
                             g_0_xxxxx_0_zzzzz_1,   \
                             g_0_xxxxx_0_zzzzzz_0,  \
                             g_0_xxxxx_0_zzzzzz_1,  \
                             g_0_xxxxxy_0_xxxxxx_0, \
                             g_0_xxxxxy_0_xxxxxy_0, \
                             g_0_xxxxxy_0_xxxxxz_0, \
                             g_0_xxxxxy_0_xxxxyy_0, \
                             g_0_xxxxxy_0_xxxxyz_0, \
                             g_0_xxxxxy_0_xxxxzz_0, \
                             g_0_xxxxxy_0_xxxyyy_0, \
                             g_0_xxxxxy_0_xxxyyz_0, \
                             g_0_xxxxxy_0_xxxyzz_0, \
                             g_0_xxxxxy_0_xxxzzz_0, \
                             g_0_xxxxxy_0_xxyyyy_0, \
                             g_0_xxxxxy_0_xxyyyz_0, \
                             g_0_xxxxxy_0_xxyyzz_0, \
                             g_0_xxxxxy_0_xxyzzz_0, \
                             g_0_xxxxxy_0_xxzzzz_0, \
                             g_0_xxxxxy_0_xyyyyy_0, \
                             g_0_xxxxxy_0_xyyyyz_0, \
                             g_0_xxxxxy_0_xyyyzz_0, \
                             g_0_xxxxxy_0_xyyzzz_0, \
                             g_0_xxxxxy_0_xyzzzz_0, \
                             g_0_xxxxxy_0_xzzzzz_0, \
                             g_0_xxxxxy_0_yyyyyy_0, \
                             g_0_xxxxxy_0_yyyyyz_0, \
                             g_0_xxxxxy_0_yyyyzz_0, \
                             g_0_xxxxxy_0_yyyzzz_0, \
                             g_0_xxxxxy_0_yyzzzz_0, \
                             g_0_xxxxxy_0_yzzzzz_0, \
                             g_0_xxxxxy_0_zzzzzz_0, \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxy_0_xxxxxx_0[i] = g_0_xxxxx_0_xxxxxx_0[i] * pb_y + g_0_xxxxx_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxy_0[i] = g_0_xxxxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxy_0[i] * pb_y + g_0_xxxxx_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxz_0[i] = g_0_xxxxx_0_xxxxxz_0[i] * pb_y + g_0_xxxxx_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxyy_0[i] = 2.0 * g_0_xxxxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyy_0[i] * pb_y + g_0_xxxxx_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxyz_0[i] = g_0_xxxxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyz_0[i] * pb_y + g_0_xxxxx_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxzz_0[i] = g_0_xxxxx_0_xxxxzz_0[i] * pb_y + g_0_xxxxx_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyyy_0[i] = 3.0 * g_0_xxxxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyy_0[i] * pb_y + g_0_xxxxx_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyyz_0[i] = 2.0 * g_0_xxxxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyz_0[i] * pb_y + g_0_xxxxx_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyzz_0[i] = g_0_xxxxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyzz_0[i] * pb_y + g_0_xxxxx_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxzzz_0[i] = g_0_xxxxx_0_xxxzzz_0[i] * pb_y + g_0_xxxxx_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyyy_0[i] = 4.0 * g_0_xxxxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyy_0[i] * pb_y + g_0_xxxxx_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyyz_0[i] = 3.0 * g_0_xxxxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyz_0[i] * pb_y + g_0_xxxxx_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyzz_0[i] = 2.0 * g_0_xxxxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyzz_0[i] * pb_y + g_0_xxxxx_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyzzz_0[i] = g_0_xxxxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyzzz_0[i] * pb_y + g_0_xxxxx_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxzzzz_0[i] = g_0_xxxxx_0_xxzzzz_0[i] * pb_y + g_0_xxxxx_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyyy_0[i] = 5.0 * g_0_xxxxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyy_0[i] * pb_y + g_0_xxxxx_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyyz_0[i] = 4.0 * g_0_xxxxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyz_0[i] * pb_y + g_0_xxxxx_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyzz_0[i] = 3.0 * g_0_xxxxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyzz_0[i] * pb_y + g_0_xxxxx_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyzzz_0[i] = 2.0 * g_0_xxxxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyzzz_0[i] * pb_y + g_0_xxxxx_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyzzzz_0[i] = g_0_xxxxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzzzz_0[i] * pb_y + g_0_xxxxx_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xzzzzz_0[i] = g_0_xxxxx_0_xzzzzz_0[i] * pb_y + g_0_xxxxx_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyyy_0[i] = 6.0 * g_0_xxxxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyy_0[i] * pb_y + g_0_xxxxx_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyyz_0[i] = 5.0 * g_0_xxxxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyz_0[i] * pb_y + g_0_xxxxx_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyzz_0[i] = 4.0 * g_0_xxxxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyzz_0[i] * pb_y + g_0_xxxxx_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyzzz_0[i] = 3.0 * g_0_xxxxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyzzz_0[i] * pb_y + g_0_xxxxx_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyzzzz_0[i] = 2.0 * g_0_xxxxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyzzzz_0[i] * pb_y + g_0_xxxxx_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yzzzzz_0[i] = g_0_xxxxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yzzzzz_0[i] * pb_y + g_0_xxxxx_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_zzzzzz_0[i] = g_0_xxxxx_0_zzzzzz_0[i] * pb_y + g_0_xxxxx_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 56-84 components of targeted buffer : SISI

    auto g_0_xxxxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 56);

    auto g_0_xxxxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 57);

    auto g_0_xxxxxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 58);

    auto g_0_xxxxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 59);

    auto g_0_xxxxxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 60);

    auto g_0_xxxxxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 61);

    auto g_0_xxxxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 62);

    auto g_0_xxxxxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 63);

    auto g_0_xxxxxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 64);

    auto g_0_xxxxxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 65);

    auto g_0_xxxxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 66);

    auto g_0_xxxxxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 67);

    auto g_0_xxxxxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 68);

    auto g_0_xxxxxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 69);

    auto g_0_xxxxxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 70);

    auto g_0_xxxxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 71);

    auto g_0_xxxxxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 72);

    auto g_0_xxxxxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 73);

    auto g_0_xxxxxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 74);

    auto g_0_xxxxxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 75);

    auto g_0_xxxxxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 76);

    auto g_0_xxxxxz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 77);

    auto g_0_xxxxxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 78);

    auto g_0_xxxxxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 79);

    auto g_0_xxxxxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 80);

    auto g_0_xxxxxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 81);

    auto g_0_xxxxxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 82);

    auto g_0_xxxxxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 83);

#pragma omp simd aligned(g_0_xxxxx_0_xxxxx_1,       \
                             g_0_xxxxx_0_xxxxxx_0,  \
                             g_0_xxxxx_0_xxxxxx_1,  \
                             g_0_xxxxx_0_xxxxxy_0,  \
                             g_0_xxxxx_0_xxxxxy_1,  \
                             g_0_xxxxx_0_xxxxxz_0,  \
                             g_0_xxxxx_0_xxxxxz_1,  \
                             g_0_xxxxx_0_xxxxy_1,   \
                             g_0_xxxxx_0_xxxxyy_0,  \
                             g_0_xxxxx_0_xxxxyy_1,  \
                             g_0_xxxxx_0_xxxxyz_0,  \
                             g_0_xxxxx_0_xxxxyz_1,  \
                             g_0_xxxxx_0_xxxxz_1,   \
                             g_0_xxxxx_0_xxxxzz_0,  \
                             g_0_xxxxx_0_xxxxzz_1,  \
                             g_0_xxxxx_0_xxxyy_1,   \
                             g_0_xxxxx_0_xxxyyy_0,  \
                             g_0_xxxxx_0_xxxyyy_1,  \
                             g_0_xxxxx_0_xxxyyz_0,  \
                             g_0_xxxxx_0_xxxyyz_1,  \
                             g_0_xxxxx_0_xxxyz_1,   \
                             g_0_xxxxx_0_xxxyzz_0,  \
                             g_0_xxxxx_0_xxxyzz_1,  \
                             g_0_xxxxx_0_xxxzz_1,   \
                             g_0_xxxxx_0_xxxzzz_0,  \
                             g_0_xxxxx_0_xxxzzz_1,  \
                             g_0_xxxxx_0_xxyyy_1,   \
                             g_0_xxxxx_0_xxyyyy_0,  \
                             g_0_xxxxx_0_xxyyyy_1,  \
                             g_0_xxxxx_0_xxyyyz_0,  \
                             g_0_xxxxx_0_xxyyyz_1,  \
                             g_0_xxxxx_0_xxyyz_1,   \
                             g_0_xxxxx_0_xxyyzz_0,  \
                             g_0_xxxxx_0_xxyyzz_1,  \
                             g_0_xxxxx_0_xxyzz_1,   \
                             g_0_xxxxx_0_xxyzzz_0,  \
                             g_0_xxxxx_0_xxyzzz_1,  \
                             g_0_xxxxx_0_xxzzz_1,   \
                             g_0_xxxxx_0_xxzzzz_0,  \
                             g_0_xxxxx_0_xxzzzz_1,  \
                             g_0_xxxxx_0_xyyyy_1,   \
                             g_0_xxxxx_0_xyyyyy_0,  \
                             g_0_xxxxx_0_xyyyyy_1,  \
                             g_0_xxxxx_0_xyyyyz_0,  \
                             g_0_xxxxx_0_xyyyyz_1,  \
                             g_0_xxxxx_0_xyyyz_1,   \
                             g_0_xxxxx_0_xyyyzz_0,  \
                             g_0_xxxxx_0_xyyyzz_1,  \
                             g_0_xxxxx_0_xyyzz_1,   \
                             g_0_xxxxx_0_xyyzzz_0,  \
                             g_0_xxxxx_0_xyyzzz_1,  \
                             g_0_xxxxx_0_xyzzz_1,   \
                             g_0_xxxxx_0_xyzzzz_0,  \
                             g_0_xxxxx_0_xyzzzz_1,  \
                             g_0_xxxxx_0_xzzzz_1,   \
                             g_0_xxxxx_0_xzzzzz_0,  \
                             g_0_xxxxx_0_xzzzzz_1,  \
                             g_0_xxxxx_0_yyyyy_1,   \
                             g_0_xxxxx_0_yyyyyy_0,  \
                             g_0_xxxxx_0_yyyyyy_1,  \
                             g_0_xxxxx_0_yyyyyz_0,  \
                             g_0_xxxxx_0_yyyyyz_1,  \
                             g_0_xxxxx_0_yyyyz_1,   \
                             g_0_xxxxx_0_yyyyzz_0,  \
                             g_0_xxxxx_0_yyyyzz_1,  \
                             g_0_xxxxx_0_yyyzz_1,   \
                             g_0_xxxxx_0_yyyzzz_0,  \
                             g_0_xxxxx_0_yyyzzz_1,  \
                             g_0_xxxxx_0_yyzzz_1,   \
                             g_0_xxxxx_0_yyzzzz_0,  \
                             g_0_xxxxx_0_yyzzzz_1,  \
                             g_0_xxxxx_0_yzzzz_1,   \
                             g_0_xxxxx_0_yzzzzz_0,  \
                             g_0_xxxxx_0_yzzzzz_1,  \
                             g_0_xxxxx_0_zzzzz_1,   \
                             g_0_xxxxx_0_zzzzzz_0,  \
                             g_0_xxxxx_0_zzzzzz_1,  \
                             g_0_xxxxxz_0_xxxxxx_0, \
                             g_0_xxxxxz_0_xxxxxy_0, \
                             g_0_xxxxxz_0_xxxxxz_0, \
                             g_0_xxxxxz_0_xxxxyy_0, \
                             g_0_xxxxxz_0_xxxxyz_0, \
                             g_0_xxxxxz_0_xxxxzz_0, \
                             g_0_xxxxxz_0_xxxyyy_0, \
                             g_0_xxxxxz_0_xxxyyz_0, \
                             g_0_xxxxxz_0_xxxyzz_0, \
                             g_0_xxxxxz_0_xxxzzz_0, \
                             g_0_xxxxxz_0_xxyyyy_0, \
                             g_0_xxxxxz_0_xxyyyz_0, \
                             g_0_xxxxxz_0_xxyyzz_0, \
                             g_0_xxxxxz_0_xxyzzz_0, \
                             g_0_xxxxxz_0_xxzzzz_0, \
                             g_0_xxxxxz_0_xyyyyy_0, \
                             g_0_xxxxxz_0_xyyyyz_0, \
                             g_0_xxxxxz_0_xyyyzz_0, \
                             g_0_xxxxxz_0_xyyzzz_0, \
                             g_0_xxxxxz_0_xyzzzz_0, \
                             g_0_xxxxxz_0_xzzzzz_0, \
                             g_0_xxxxxz_0_yyyyyy_0, \
                             g_0_xxxxxz_0_yyyyyz_0, \
                             g_0_xxxxxz_0_yyyyzz_0, \
                             g_0_xxxxxz_0_yyyzzz_0, \
                             g_0_xxxxxz_0_yyzzzz_0, \
                             g_0_xxxxxz_0_yzzzzz_0, \
                             g_0_xxxxxz_0_zzzzzz_0, \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxz_0_xxxxxx_0[i] = g_0_xxxxx_0_xxxxxx_0[i] * pb_z + g_0_xxxxx_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxy_0[i] = g_0_xxxxx_0_xxxxxy_0[i] * pb_z + g_0_xxxxx_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxz_0[i] = g_0_xxxxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxz_0[i] * pb_z + g_0_xxxxx_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxyy_0[i] = g_0_xxxxx_0_xxxxyy_0[i] * pb_z + g_0_xxxxx_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxyz_0[i] = g_0_xxxxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyz_0[i] * pb_z + g_0_xxxxx_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxzz_0[i] = 2.0 * g_0_xxxxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxzz_0[i] * pb_z + g_0_xxxxx_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyyy_0[i] = g_0_xxxxx_0_xxxyyy_0[i] * pb_z + g_0_xxxxx_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyyz_0[i] = g_0_xxxxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyz_0[i] * pb_z + g_0_xxxxx_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyzz_0[i] = 2.0 * g_0_xxxxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyzz_0[i] * pb_z + g_0_xxxxx_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxzzz_0[i] = 3.0 * g_0_xxxxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxzzz_0[i] * pb_z + g_0_xxxxx_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyyy_0[i] = g_0_xxxxx_0_xxyyyy_0[i] * pb_z + g_0_xxxxx_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyyz_0[i] = g_0_xxxxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyz_0[i] * pb_z + g_0_xxxxx_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyzz_0[i] = 2.0 * g_0_xxxxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyzz_0[i] * pb_z + g_0_xxxxx_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyzzz_0[i] = 3.0 * g_0_xxxxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyzzz_0[i] * pb_z + g_0_xxxxx_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxzzzz_0[i] = 4.0 * g_0_xxxxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxzzzz_0[i] * pb_z + g_0_xxxxx_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyyy_0[i] = g_0_xxxxx_0_xyyyyy_0[i] * pb_z + g_0_xxxxx_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyyz_0[i] = g_0_xxxxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyz_0[i] * pb_z + g_0_xxxxx_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyzz_0[i] = 2.0 * g_0_xxxxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyzz_0[i] * pb_z + g_0_xxxxx_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyzzz_0[i] = 3.0 * g_0_xxxxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyzzz_0[i] * pb_z + g_0_xxxxx_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyzzzz_0[i] = 4.0 * g_0_xxxxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzzzz_0[i] * pb_z + g_0_xxxxx_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xzzzzz_0[i] = 5.0 * g_0_xxxxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xzzzzz_0[i] * pb_z + g_0_xxxxx_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyyy_0[i] = g_0_xxxxx_0_yyyyyy_0[i] * pb_z + g_0_xxxxx_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyyz_0[i] = g_0_xxxxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyz_0[i] * pb_z + g_0_xxxxx_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyzz_0[i] = 2.0 * g_0_xxxxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyzz_0[i] * pb_z + g_0_xxxxx_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyzzz_0[i] = 3.0 * g_0_xxxxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyzzz_0[i] * pb_z + g_0_xxxxx_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyzzzz_0[i] = 4.0 * g_0_xxxxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyzzzz_0[i] * pb_z + g_0_xxxxx_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yzzzzz_0[i] = 5.0 * g_0_xxxxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yzzzzz_0[i] * pb_z + g_0_xxxxx_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_zzzzzz_0[i] = 6.0 * g_0_xxxxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_zzzzzz_0[i] * pb_z + g_0_xxxxx_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 84-112 components of targeted buffer : SISI

    auto g_0_xxxxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 84);

    auto g_0_xxxxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 85);

    auto g_0_xxxxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 86);

    auto g_0_xxxxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 87);

    auto g_0_xxxxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 88);

    auto g_0_xxxxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 89);

    auto g_0_xxxxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 90);

    auto g_0_xxxxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 91);

    auto g_0_xxxxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 92);

    auto g_0_xxxxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 93);

    auto g_0_xxxxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 94);

    auto g_0_xxxxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 95);

    auto g_0_xxxxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 96);

    auto g_0_xxxxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 97);

    auto g_0_xxxxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 98);

    auto g_0_xxxxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 99);

    auto g_0_xxxxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 100);

    auto g_0_xxxxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 101);

    auto g_0_xxxxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 102);

    auto g_0_xxxxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 103);

    auto g_0_xxxxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 104);

    auto g_0_xxxxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 105);

    auto g_0_xxxxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 106);

    auto g_0_xxxxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 107);

    auto g_0_xxxxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 108);

    auto g_0_xxxxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 109);

    auto g_0_xxxxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 110);

    auto g_0_xxxxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 111);

#pragma omp simd aligned(g_0_xxxx_0_xxxxxx_0,       \
                             g_0_xxxx_0_xxxxxx_1,   \
                             g_0_xxxx_0_xxxxxz_0,   \
                             g_0_xxxx_0_xxxxxz_1,   \
                             g_0_xxxx_0_xxxxzz_0,   \
                             g_0_xxxx_0_xxxxzz_1,   \
                             g_0_xxxx_0_xxxzzz_0,   \
                             g_0_xxxx_0_xxxzzz_1,   \
                             g_0_xxxx_0_xxzzzz_0,   \
                             g_0_xxxx_0_xxzzzz_1,   \
                             g_0_xxxx_0_xzzzzz_0,   \
                             g_0_xxxx_0_xzzzzz_1,   \
                             g_0_xxxxy_0_xxxxxx_0,  \
                             g_0_xxxxy_0_xxxxxx_1,  \
                             g_0_xxxxy_0_xxxxxz_0,  \
                             g_0_xxxxy_0_xxxxxz_1,  \
                             g_0_xxxxy_0_xxxxzz_0,  \
                             g_0_xxxxy_0_xxxxzz_1,  \
                             g_0_xxxxy_0_xxxzzz_0,  \
                             g_0_xxxxy_0_xxxzzz_1,  \
                             g_0_xxxxy_0_xxzzzz_0,  \
                             g_0_xxxxy_0_xxzzzz_1,  \
                             g_0_xxxxy_0_xzzzzz_0,  \
                             g_0_xxxxy_0_xzzzzz_1,  \
                             g_0_xxxxyy_0_xxxxxx_0, \
                             g_0_xxxxyy_0_xxxxxy_0, \
                             g_0_xxxxyy_0_xxxxxz_0, \
                             g_0_xxxxyy_0_xxxxyy_0, \
                             g_0_xxxxyy_0_xxxxyz_0, \
                             g_0_xxxxyy_0_xxxxzz_0, \
                             g_0_xxxxyy_0_xxxyyy_0, \
                             g_0_xxxxyy_0_xxxyyz_0, \
                             g_0_xxxxyy_0_xxxyzz_0, \
                             g_0_xxxxyy_0_xxxzzz_0, \
                             g_0_xxxxyy_0_xxyyyy_0, \
                             g_0_xxxxyy_0_xxyyyz_0, \
                             g_0_xxxxyy_0_xxyyzz_0, \
                             g_0_xxxxyy_0_xxyzzz_0, \
                             g_0_xxxxyy_0_xxzzzz_0, \
                             g_0_xxxxyy_0_xyyyyy_0, \
                             g_0_xxxxyy_0_xyyyyz_0, \
                             g_0_xxxxyy_0_xyyyzz_0, \
                             g_0_xxxxyy_0_xyyzzz_0, \
                             g_0_xxxxyy_0_xyzzzz_0, \
                             g_0_xxxxyy_0_xzzzzz_0, \
                             g_0_xxxxyy_0_yyyyyy_0, \
                             g_0_xxxxyy_0_yyyyyz_0, \
                             g_0_xxxxyy_0_yyyyzz_0, \
                             g_0_xxxxyy_0_yyyzzz_0, \
                             g_0_xxxxyy_0_yyzzzz_0, \
                             g_0_xxxxyy_0_yzzzzz_0, \
                             g_0_xxxxyy_0_zzzzzz_0, \
                             g_0_xxxyy_0_xxxxxy_0,  \
                             g_0_xxxyy_0_xxxxxy_1,  \
                             g_0_xxxyy_0_xxxxy_1,   \
                             g_0_xxxyy_0_xxxxyy_0,  \
                             g_0_xxxyy_0_xxxxyy_1,  \
                             g_0_xxxyy_0_xxxxyz_0,  \
                             g_0_xxxyy_0_xxxxyz_1,  \
                             g_0_xxxyy_0_xxxyy_1,   \
                             g_0_xxxyy_0_xxxyyy_0,  \
                             g_0_xxxyy_0_xxxyyy_1,  \
                             g_0_xxxyy_0_xxxyyz_0,  \
                             g_0_xxxyy_0_xxxyyz_1,  \
                             g_0_xxxyy_0_xxxyz_1,   \
                             g_0_xxxyy_0_xxxyzz_0,  \
                             g_0_xxxyy_0_xxxyzz_1,  \
                             g_0_xxxyy_0_xxyyy_1,   \
                             g_0_xxxyy_0_xxyyyy_0,  \
                             g_0_xxxyy_0_xxyyyy_1,  \
                             g_0_xxxyy_0_xxyyyz_0,  \
                             g_0_xxxyy_0_xxyyyz_1,  \
                             g_0_xxxyy_0_xxyyz_1,   \
                             g_0_xxxyy_0_xxyyzz_0,  \
                             g_0_xxxyy_0_xxyyzz_1,  \
                             g_0_xxxyy_0_xxyzz_1,   \
                             g_0_xxxyy_0_xxyzzz_0,  \
                             g_0_xxxyy_0_xxyzzz_1,  \
                             g_0_xxxyy_0_xyyyy_1,   \
                             g_0_xxxyy_0_xyyyyy_0,  \
                             g_0_xxxyy_0_xyyyyy_1,  \
                             g_0_xxxyy_0_xyyyyz_0,  \
                             g_0_xxxyy_0_xyyyyz_1,  \
                             g_0_xxxyy_0_xyyyz_1,   \
                             g_0_xxxyy_0_xyyyzz_0,  \
                             g_0_xxxyy_0_xyyyzz_1,  \
                             g_0_xxxyy_0_xyyzz_1,   \
                             g_0_xxxyy_0_xyyzzz_0,  \
                             g_0_xxxyy_0_xyyzzz_1,  \
                             g_0_xxxyy_0_xyzzz_1,   \
                             g_0_xxxyy_0_xyzzzz_0,  \
                             g_0_xxxyy_0_xyzzzz_1,  \
                             g_0_xxxyy_0_yyyyy_1,   \
                             g_0_xxxyy_0_yyyyyy_0,  \
                             g_0_xxxyy_0_yyyyyy_1,  \
                             g_0_xxxyy_0_yyyyyz_0,  \
                             g_0_xxxyy_0_yyyyyz_1,  \
                             g_0_xxxyy_0_yyyyz_1,   \
                             g_0_xxxyy_0_yyyyzz_0,  \
                             g_0_xxxyy_0_yyyyzz_1,  \
                             g_0_xxxyy_0_yyyzz_1,   \
                             g_0_xxxyy_0_yyyzzz_0,  \
                             g_0_xxxyy_0_yyyzzz_1,  \
                             g_0_xxxyy_0_yyzzz_1,   \
                             g_0_xxxyy_0_yyzzzz_0,  \
                             g_0_xxxyy_0_yyzzzz_1,  \
                             g_0_xxxyy_0_yzzzz_1,   \
                             g_0_xxxyy_0_yzzzzz_0,  \
                             g_0_xxxyy_0_yzzzzz_1,  \
                             g_0_xxxyy_0_zzzzzz_0,  \
                             g_0_xxxyy_0_zzzzzz_1,  \
                             g_0_xxyy_0_xxxxxy_0,   \
                             g_0_xxyy_0_xxxxxy_1,   \
                             g_0_xxyy_0_xxxxyy_0,   \
                             g_0_xxyy_0_xxxxyy_1,   \
                             g_0_xxyy_0_xxxxyz_0,   \
                             g_0_xxyy_0_xxxxyz_1,   \
                             g_0_xxyy_0_xxxyyy_0,   \
                             g_0_xxyy_0_xxxyyy_1,   \
                             g_0_xxyy_0_xxxyyz_0,   \
                             g_0_xxyy_0_xxxyyz_1,   \
                             g_0_xxyy_0_xxxyzz_0,   \
                             g_0_xxyy_0_xxxyzz_1,   \
                             g_0_xxyy_0_xxyyyy_0,   \
                             g_0_xxyy_0_xxyyyy_1,   \
                             g_0_xxyy_0_xxyyyz_0,   \
                             g_0_xxyy_0_xxyyyz_1,   \
                             g_0_xxyy_0_xxyyzz_0,   \
                             g_0_xxyy_0_xxyyzz_1,   \
                             g_0_xxyy_0_xxyzzz_0,   \
                             g_0_xxyy_0_xxyzzz_1,   \
                             g_0_xxyy_0_xyyyyy_0,   \
                             g_0_xxyy_0_xyyyyy_1,   \
                             g_0_xxyy_0_xyyyyz_0,   \
                             g_0_xxyy_0_xyyyyz_1,   \
                             g_0_xxyy_0_xyyyzz_0,   \
                             g_0_xxyy_0_xyyyzz_1,   \
                             g_0_xxyy_0_xyyzzz_0,   \
                             g_0_xxyy_0_xyyzzz_1,   \
                             g_0_xxyy_0_xyzzzz_0,   \
                             g_0_xxyy_0_xyzzzz_1,   \
                             g_0_xxyy_0_yyyyyy_0,   \
                             g_0_xxyy_0_yyyyyy_1,   \
                             g_0_xxyy_0_yyyyyz_0,   \
                             g_0_xxyy_0_yyyyyz_1,   \
                             g_0_xxyy_0_yyyyzz_0,   \
                             g_0_xxyy_0_yyyyzz_1,   \
                             g_0_xxyy_0_yyyzzz_0,   \
                             g_0_xxyy_0_yyyzzz_1,   \
                             g_0_xxyy_0_yyzzzz_0,   \
                             g_0_xxyy_0_yyzzzz_1,   \
                             g_0_xxyy_0_yzzzzz_0,   \
                             g_0_xxyy_0_yzzzzz_1,   \
                             g_0_xxyy_0_zzzzzz_0,   \
                             g_0_xxyy_0_zzzzzz_1,   \
                             wp_x,                  \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyy_0_xxxxxx_0[i] =
            g_0_xxxx_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxxx_0[i] * pb_y + g_0_xxxxy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxxxy_0[i] = 3.0 * g_0_xxyy_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxy_1[i] * fti_ab_0 +
                                   5.0 * g_0_xxxyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxy_0[i] * pb_x + g_0_xxxyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxxz_0[i] =
            g_0_xxxx_0_xxxxxz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxxz_0[i] * pb_y + g_0_xxxxy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxxyy_0[i] = 3.0 * g_0_xxyy_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxyy_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxxyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyy_0[i] * pb_x + g_0_xxxyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxyz_0[i] = 3.0 * g_0_xxyy_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxyz_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxxyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyz_0[i] * pb_x + g_0_xxxyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxzz_0[i] =
            g_0_xxxx_0_xxxxzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxzz_0[i] * pb_y + g_0_xxxxy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxyyy_0[i] = 3.0 * g_0_xxyy_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyyy_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxxyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyy_0[i] * pb_x + g_0_xxxyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxyyz_0[i] = 3.0 * g_0_xxyy_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyyz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxxyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyz_0[i] * pb_x + g_0_xxxyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxyzz_0[i] = 3.0 * g_0_xxyy_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxxyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyzz_0[i] * pb_x + g_0_xxxyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxzzz_0[i] =
            g_0_xxxx_0_xxxzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxzzz_0[i] * pb_y + g_0_xxxxy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxyyyy_0[i] = 3.0 * g_0_xxyy_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyyy_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyy_0[i] * pb_x + g_0_xxxyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyyyz_0[i] = 3.0 * g_0_xxyy_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyyz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyz_0[i] * pb_x + g_0_xxxyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyyzz_0[i] = 3.0 * g_0_xxyy_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyzz_0[i] * pb_x + g_0_xxxyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyzzz_0[i] = 3.0 * g_0_xxyy_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyzzz_0[i] * pb_x + g_0_xxxyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxzzzz_0[i] =
            g_0_xxxx_0_xxzzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxzzzz_0[i] * pb_y + g_0_xxxxy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xyyyyy_0[i] = 3.0 * g_0_xxyy_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyyy_1[i] * fti_ab_0 +
                                   g_0_xxxyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyy_0[i] * pb_x + g_0_xxxyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyyyz_0[i] = 3.0 * g_0_xxyy_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyyz_1[i] * fti_ab_0 +
                                   g_0_xxxyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyz_0[i] * pb_x + g_0_xxxyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyyzz_0[i] = 3.0 * g_0_xxyy_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyzz_1[i] * fti_ab_0 +
                                   g_0_xxxyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyzz_0[i] * pb_x + g_0_xxxyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyzzz_0[i] = 3.0 * g_0_xxyy_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyzzz_1[i] * fti_ab_0 +
                                   g_0_xxxyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyzzz_0[i] * pb_x + g_0_xxxyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyzzzz_0[i] = 3.0 * g_0_xxyy_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyzzzz_1[i] * fti_ab_0 +
                                   g_0_xxxyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyzzzz_0[i] * pb_x + g_0_xxxyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xzzzzz_0[i] =
            g_0_xxxx_0_xzzzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xzzzzz_0[i] * pb_y + g_0_xxxxy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_yyyyyy_0[i] = 3.0 * g_0_xxyy_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyyy_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyyyy_0[i] * pb_x +
                                   g_0_xxxyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyyyz_0[i] = 3.0 * g_0_xxyy_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyyz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyyyz_0[i] * pb_x +
                                   g_0_xxxyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyyzz_0[i] = 3.0 * g_0_xxyy_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyyzz_0[i] * pb_x +
                                   g_0_xxxyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyzzz_0[i] = 3.0 * g_0_xxyy_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyzzz_0[i] * pb_x +
                                   g_0_xxxyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyzzzz_0[i] = 3.0 * g_0_xxyy_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyzzzz_0[i] * pb_x +
                                   g_0_xxxyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yzzzzz_0[i] = 3.0 * g_0_xxyy_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yzzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yzzzzz_0[i] * pb_x +
                                   g_0_xxxyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_zzzzzz_0[i] = 3.0 * g_0_xxyy_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_zzzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_zzzzzz_0[i] * pb_x +
                                   g_0_xxxyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 112-140 components of targeted buffer : SISI

    auto g_0_xxxxyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 112);

    auto g_0_xxxxyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 113);

    auto g_0_xxxxyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 114);

    auto g_0_xxxxyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 115);

    auto g_0_xxxxyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 116);

    auto g_0_xxxxyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 117);

    auto g_0_xxxxyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 118);

    auto g_0_xxxxyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 119);

    auto g_0_xxxxyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 120);

    auto g_0_xxxxyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 121);

    auto g_0_xxxxyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 122);

    auto g_0_xxxxyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 123);

    auto g_0_xxxxyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 124);

    auto g_0_xxxxyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 125);

    auto g_0_xxxxyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 126);

    auto g_0_xxxxyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 127);

    auto g_0_xxxxyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 128);

    auto g_0_xxxxyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 129);

    auto g_0_xxxxyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 130);

    auto g_0_xxxxyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 131);

    auto g_0_xxxxyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 132);

    auto g_0_xxxxyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 133);

    auto g_0_xxxxyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 134);

    auto g_0_xxxxyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 135);

    auto g_0_xxxxyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 136);

    auto g_0_xxxxyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 137);

    auto g_0_xxxxyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 138);

    auto g_0_xxxxyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 139);

#pragma omp simd aligned(g_0_xxxxy_0_xxxxxy_0,      \
                             g_0_xxxxy_0_xxxxxy_1,  \
                             g_0_xxxxy_0_xxxxyy_0,  \
                             g_0_xxxxy_0_xxxxyy_1,  \
                             g_0_xxxxy_0_xxxyyy_0,  \
                             g_0_xxxxy_0_xxxyyy_1,  \
                             g_0_xxxxy_0_xxyyyy_0,  \
                             g_0_xxxxy_0_xxyyyy_1,  \
                             g_0_xxxxy_0_xyyyyy_0,  \
                             g_0_xxxxy_0_xyyyyy_1,  \
                             g_0_xxxxy_0_yyyyyy_0,  \
                             g_0_xxxxy_0_yyyyyy_1,  \
                             g_0_xxxxyz_0_xxxxxx_0, \
                             g_0_xxxxyz_0_xxxxxy_0, \
                             g_0_xxxxyz_0_xxxxxz_0, \
                             g_0_xxxxyz_0_xxxxyy_0, \
                             g_0_xxxxyz_0_xxxxyz_0, \
                             g_0_xxxxyz_0_xxxxzz_0, \
                             g_0_xxxxyz_0_xxxyyy_0, \
                             g_0_xxxxyz_0_xxxyyz_0, \
                             g_0_xxxxyz_0_xxxyzz_0, \
                             g_0_xxxxyz_0_xxxzzz_0, \
                             g_0_xxxxyz_0_xxyyyy_0, \
                             g_0_xxxxyz_0_xxyyyz_0, \
                             g_0_xxxxyz_0_xxyyzz_0, \
                             g_0_xxxxyz_0_xxyzzz_0, \
                             g_0_xxxxyz_0_xxzzzz_0, \
                             g_0_xxxxyz_0_xyyyyy_0, \
                             g_0_xxxxyz_0_xyyyyz_0, \
                             g_0_xxxxyz_0_xyyyzz_0, \
                             g_0_xxxxyz_0_xyyzzz_0, \
                             g_0_xxxxyz_0_xyzzzz_0, \
                             g_0_xxxxyz_0_xzzzzz_0, \
                             g_0_xxxxyz_0_yyyyyy_0, \
                             g_0_xxxxyz_0_yyyyyz_0, \
                             g_0_xxxxyz_0_yyyyzz_0, \
                             g_0_xxxxyz_0_yyyzzz_0, \
                             g_0_xxxxyz_0_yyzzzz_0, \
                             g_0_xxxxyz_0_yzzzzz_0, \
                             g_0_xxxxyz_0_zzzzzz_0, \
                             g_0_xxxxz_0_xxxxxx_0,  \
                             g_0_xxxxz_0_xxxxxx_1,  \
                             g_0_xxxxz_0_xxxxxz_0,  \
                             g_0_xxxxz_0_xxxxxz_1,  \
                             g_0_xxxxz_0_xxxxyz_0,  \
                             g_0_xxxxz_0_xxxxyz_1,  \
                             g_0_xxxxz_0_xxxxz_1,   \
                             g_0_xxxxz_0_xxxxzz_0,  \
                             g_0_xxxxz_0_xxxxzz_1,  \
                             g_0_xxxxz_0_xxxyyz_0,  \
                             g_0_xxxxz_0_xxxyyz_1,  \
                             g_0_xxxxz_0_xxxyz_1,   \
                             g_0_xxxxz_0_xxxyzz_0,  \
                             g_0_xxxxz_0_xxxyzz_1,  \
                             g_0_xxxxz_0_xxxzz_1,   \
                             g_0_xxxxz_0_xxxzzz_0,  \
                             g_0_xxxxz_0_xxxzzz_1,  \
                             g_0_xxxxz_0_xxyyyz_0,  \
                             g_0_xxxxz_0_xxyyyz_1,  \
                             g_0_xxxxz_0_xxyyz_1,   \
                             g_0_xxxxz_0_xxyyzz_0,  \
                             g_0_xxxxz_0_xxyyzz_1,  \
                             g_0_xxxxz_0_xxyzz_1,   \
                             g_0_xxxxz_0_xxyzzz_0,  \
                             g_0_xxxxz_0_xxyzzz_1,  \
                             g_0_xxxxz_0_xxzzz_1,   \
                             g_0_xxxxz_0_xxzzzz_0,  \
                             g_0_xxxxz_0_xxzzzz_1,  \
                             g_0_xxxxz_0_xyyyyz_0,  \
                             g_0_xxxxz_0_xyyyyz_1,  \
                             g_0_xxxxz_0_xyyyz_1,   \
                             g_0_xxxxz_0_xyyyzz_0,  \
                             g_0_xxxxz_0_xyyyzz_1,  \
                             g_0_xxxxz_0_xyyzz_1,   \
                             g_0_xxxxz_0_xyyzzz_0,  \
                             g_0_xxxxz_0_xyyzzz_1,  \
                             g_0_xxxxz_0_xyzzz_1,   \
                             g_0_xxxxz_0_xyzzzz_0,  \
                             g_0_xxxxz_0_xyzzzz_1,  \
                             g_0_xxxxz_0_xzzzz_1,   \
                             g_0_xxxxz_0_xzzzzz_0,  \
                             g_0_xxxxz_0_xzzzzz_1,  \
                             g_0_xxxxz_0_yyyyyz_0,  \
                             g_0_xxxxz_0_yyyyyz_1,  \
                             g_0_xxxxz_0_yyyyz_1,   \
                             g_0_xxxxz_0_yyyyzz_0,  \
                             g_0_xxxxz_0_yyyyzz_1,  \
                             g_0_xxxxz_0_yyyzz_1,   \
                             g_0_xxxxz_0_yyyzzz_0,  \
                             g_0_xxxxz_0_yyyzzz_1,  \
                             g_0_xxxxz_0_yyzzz_1,   \
                             g_0_xxxxz_0_yyzzzz_0,  \
                             g_0_xxxxz_0_yyzzzz_1,  \
                             g_0_xxxxz_0_yzzzz_1,   \
                             g_0_xxxxz_0_yzzzzz_0,  \
                             g_0_xxxxz_0_yzzzzz_1,  \
                             g_0_xxxxz_0_zzzzz_1,   \
                             g_0_xxxxz_0_zzzzzz_0,  \
                             g_0_xxxxz_0_zzzzzz_1,  \
                             wp_y,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyz_0_xxxxxx_0[i] = g_0_xxxxz_0_xxxxxx_0[i] * pb_y + g_0_xxxxz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxxy_0[i] = g_0_xxxxy_0_xxxxxy_0[i] * pb_z + g_0_xxxxy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxxxz_0[i] = g_0_xxxxz_0_xxxxxz_0[i] * pb_y + g_0_xxxxz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxyy_0[i] = g_0_xxxxy_0_xxxxyy_0[i] * pb_z + g_0_xxxxy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxxyz_0[i] = g_0_xxxxz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxxyz_0[i] * pb_y + g_0_xxxxz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxzz_0[i] = g_0_xxxxz_0_xxxxzz_0[i] * pb_y + g_0_xxxxz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxyyy_0[i] = g_0_xxxxy_0_xxxyyy_0[i] * pb_z + g_0_xxxxy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxyyz_0[i] = 2.0 * g_0_xxxxz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxyyz_0[i] * pb_y + g_0_xxxxz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxyzz_0[i] = g_0_xxxxz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxyzz_0[i] * pb_y + g_0_xxxxz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxzzz_0[i] = g_0_xxxxz_0_xxxzzz_0[i] * pb_y + g_0_xxxxz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyyyy_0[i] = g_0_xxxxy_0_xxyyyy_0[i] * pb_z + g_0_xxxxy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxyyyz_0[i] = 3.0 * g_0_xxxxz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyyyz_0[i] * pb_y + g_0_xxxxz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyyzz_0[i] = 2.0 * g_0_xxxxz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyyzz_0[i] * pb_y + g_0_xxxxz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyzzz_0[i] = g_0_xxxxz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyzzz_0[i] * pb_y + g_0_xxxxz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxzzzz_0[i] = g_0_xxxxz_0_xxzzzz_0[i] * pb_y + g_0_xxxxz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyyyy_0[i] = g_0_xxxxy_0_xyyyyy_0[i] * pb_z + g_0_xxxxy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xyyyyz_0[i] = 4.0 * g_0_xxxxz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyyyz_0[i] * pb_y + g_0_xxxxz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyyzz_0[i] = 3.0 * g_0_xxxxz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyyzz_0[i] * pb_y + g_0_xxxxz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyzzz_0[i] = 2.0 * g_0_xxxxz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyzzz_0[i] * pb_y + g_0_xxxxz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyzzzz_0[i] = g_0_xxxxz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyzzzz_0[i] * pb_y + g_0_xxxxz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xzzzzz_0[i] = g_0_xxxxz_0_xzzzzz_0[i] * pb_y + g_0_xxxxz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyyyy_0[i] = g_0_xxxxy_0_yyyyyy_0[i] * pb_z + g_0_xxxxy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_yyyyyz_0[i] = 5.0 * g_0_xxxxz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyyyz_0[i] * pb_y + g_0_xxxxz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyyzz_0[i] = 4.0 * g_0_xxxxz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyyzz_0[i] * pb_y + g_0_xxxxz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyzzz_0[i] = 3.0 * g_0_xxxxz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyzzz_0[i] * pb_y + g_0_xxxxz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyzzzz_0[i] = 2.0 * g_0_xxxxz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyzzzz_0[i] * pb_y + g_0_xxxxz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yzzzzz_0[i] = g_0_xxxxz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yzzzzz_0[i] * pb_y + g_0_xxxxz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_zzzzzz_0[i] = g_0_xxxxz_0_zzzzzz_0[i] * pb_y + g_0_xxxxz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 140-168 components of targeted buffer : SISI

    auto g_0_xxxxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 140);

    auto g_0_xxxxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 141);

    auto g_0_xxxxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 142);

    auto g_0_xxxxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 143);

    auto g_0_xxxxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 144);

    auto g_0_xxxxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 145);

    auto g_0_xxxxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 146);

    auto g_0_xxxxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 147);

    auto g_0_xxxxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 148);

    auto g_0_xxxxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 149);

    auto g_0_xxxxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 150);

    auto g_0_xxxxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 151);

    auto g_0_xxxxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 152);

    auto g_0_xxxxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 153);

    auto g_0_xxxxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 154);

    auto g_0_xxxxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 155);

    auto g_0_xxxxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 156);

    auto g_0_xxxxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 157);

    auto g_0_xxxxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 158);

    auto g_0_xxxxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 159);

    auto g_0_xxxxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 160);

    auto g_0_xxxxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 161);

    auto g_0_xxxxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 162);

    auto g_0_xxxxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 163);

    auto g_0_xxxxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 164);

    auto g_0_xxxxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 165);

    auto g_0_xxxxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 166);

    auto g_0_xxxxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 167);

#pragma omp simd aligned(g_0_xxxx_0_xxxxxx_0,       \
                             g_0_xxxx_0_xxxxxx_1,   \
                             g_0_xxxx_0_xxxxxy_0,   \
                             g_0_xxxx_0_xxxxxy_1,   \
                             g_0_xxxx_0_xxxxyy_0,   \
                             g_0_xxxx_0_xxxxyy_1,   \
                             g_0_xxxx_0_xxxyyy_0,   \
                             g_0_xxxx_0_xxxyyy_1,   \
                             g_0_xxxx_0_xxyyyy_0,   \
                             g_0_xxxx_0_xxyyyy_1,   \
                             g_0_xxxx_0_xyyyyy_0,   \
                             g_0_xxxx_0_xyyyyy_1,   \
                             g_0_xxxxz_0_xxxxxx_0,  \
                             g_0_xxxxz_0_xxxxxx_1,  \
                             g_0_xxxxz_0_xxxxxy_0,  \
                             g_0_xxxxz_0_xxxxxy_1,  \
                             g_0_xxxxz_0_xxxxyy_0,  \
                             g_0_xxxxz_0_xxxxyy_1,  \
                             g_0_xxxxz_0_xxxyyy_0,  \
                             g_0_xxxxz_0_xxxyyy_1,  \
                             g_0_xxxxz_0_xxyyyy_0,  \
                             g_0_xxxxz_0_xxyyyy_1,  \
                             g_0_xxxxz_0_xyyyyy_0,  \
                             g_0_xxxxz_0_xyyyyy_1,  \
                             g_0_xxxxzz_0_xxxxxx_0, \
                             g_0_xxxxzz_0_xxxxxy_0, \
                             g_0_xxxxzz_0_xxxxxz_0, \
                             g_0_xxxxzz_0_xxxxyy_0, \
                             g_0_xxxxzz_0_xxxxyz_0, \
                             g_0_xxxxzz_0_xxxxzz_0, \
                             g_0_xxxxzz_0_xxxyyy_0, \
                             g_0_xxxxzz_0_xxxyyz_0, \
                             g_0_xxxxzz_0_xxxyzz_0, \
                             g_0_xxxxzz_0_xxxzzz_0, \
                             g_0_xxxxzz_0_xxyyyy_0, \
                             g_0_xxxxzz_0_xxyyyz_0, \
                             g_0_xxxxzz_0_xxyyzz_0, \
                             g_0_xxxxzz_0_xxyzzz_0, \
                             g_0_xxxxzz_0_xxzzzz_0, \
                             g_0_xxxxzz_0_xyyyyy_0, \
                             g_0_xxxxzz_0_xyyyyz_0, \
                             g_0_xxxxzz_0_xyyyzz_0, \
                             g_0_xxxxzz_0_xyyzzz_0, \
                             g_0_xxxxzz_0_xyzzzz_0, \
                             g_0_xxxxzz_0_xzzzzz_0, \
                             g_0_xxxxzz_0_yyyyyy_0, \
                             g_0_xxxxzz_0_yyyyyz_0, \
                             g_0_xxxxzz_0_yyyyzz_0, \
                             g_0_xxxxzz_0_yyyzzz_0, \
                             g_0_xxxxzz_0_yyzzzz_0, \
                             g_0_xxxxzz_0_yzzzzz_0, \
                             g_0_xxxxzz_0_zzzzzz_0, \
                             g_0_xxxzz_0_xxxxxz_0,  \
                             g_0_xxxzz_0_xxxxxz_1,  \
                             g_0_xxxzz_0_xxxxyz_0,  \
                             g_0_xxxzz_0_xxxxyz_1,  \
                             g_0_xxxzz_0_xxxxz_1,   \
                             g_0_xxxzz_0_xxxxzz_0,  \
                             g_0_xxxzz_0_xxxxzz_1,  \
                             g_0_xxxzz_0_xxxyyz_0,  \
                             g_0_xxxzz_0_xxxyyz_1,  \
                             g_0_xxxzz_0_xxxyz_1,   \
                             g_0_xxxzz_0_xxxyzz_0,  \
                             g_0_xxxzz_0_xxxyzz_1,  \
                             g_0_xxxzz_0_xxxzz_1,   \
                             g_0_xxxzz_0_xxxzzz_0,  \
                             g_0_xxxzz_0_xxxzzz_1,  \
                             g_0_xxxzz_0_xxyyyz_0,  \
                             g_0_xxxzz_0_xxyyyz_1,  \
                             g_0_xxxzz_0_xxyyz_1,   \
                             g_0_xxxzz_0_xxyyzz_0,  \
                             g_0_xxxzz_0_xxyyzz_1,  \
                             g_0_xxxzz_0_xxyzz_1,   \
                             g_0_xxxzz_0_xxyzzz_0,  \
                             g_0_xxxzz_0_xxyzzz_1,  \
                             g_0_xxxzz_0_xxzzz_1,   \
                             g_0_xxxzz_0_xxzzzz_0,  \
                             g_0_xxxzz_0_xxzzzz_1,  \
                             g_0_xxxzz_0_xyyyyz_0,  \
                             g_0_xxxzz_0_xyyyyz_1,  \
                             g_0_xxxzz_0_xyyyz_1,   \
                             g_0_xxxzz_0_xyyyzz_0,  \
                             g_0_xxxzz_0_xyyyzz_1,  \
                             g_0_xxxzz_0_xyyzz_1,   \
                             g_0_xxxzz_0_xyyzzz_0,  \
                             g_0_xxxzz_0_xyyzzz_1,  \
                             g_0_xxxzz_0_xyzzz_1,   \
                             g_0_xxxzz_0_xyzzzz_0,  \
                             g_0_xxxzz_0_xyzzzz_1,  \
                             g_0_xxxzz_0_xzzzz_1,   \
                             g_0_xxxzz_0_xzzzzz_0,  \
                             g_0_xxxzz_0_xzzzzz_1,  \
                             g_0_xxxzz_0_yyyyyy_0,  \
                             g_0_xxxzz_0_yyyyyy_1,  \
                             g_0_xxxzz_0_yyyyyz_0,  \
                             g_0_xxxzz_0_yyyyyz_1,  \
                             g_0_xxxzz_0_yyyyz_1,   \
                             g_0_xxxzz_0_yyyyzz_0,  \
                             g_0_xxxzz_0_yyyyzz_1,  \
                             g_0_xxxzz_0_yyyzz_1,   \
                             g_0_xxxzz_0_yyyzzz_0,  \
                             g_0_xxxzz_0_yyyzzz_1,  \
                             g_0_xxxzz_0_yyzzz_1,   \
                             g_0_xxxzz_0_yyzzzz_0,  \
                             g_0_xxxzz_0_yyzzzz_1,  \
                             g_0_xxxzz_0_yzzzz_1,   \
                             g_0_xxxzz_0_yzzzzz_0,  \
                             g_0_xxxzz_0_yzzzzz_1,  \
                             g_0_xxxzz_0_zzzzz_1,   \
                             g_0_xxxzz_0_zzzzzz_0,  \
                             g_0_xxxzz_0_zzzzzz_1,  \
                             g_0_xxzz_0_xxxxxz_0,   \
                             g_0_xxzz_0_xxxxxz_1,   \
                             g_0_xxzz_0_xxxxyz_0,   \
                             g_0_xxzz_0_xxxxyz_1,   \
                             g_0_xxzz_0_xxxxzz_0,   \
                             g_0_xxzz_0_xxxxzz_1,   \
                             g_0_xxzz_0_xxxyyz_0,   \
                             g_0_xxzz_0_xxxyyz_1,   \
                             g_0_xxzz_0_xxxyzz_0,   \
                             g_0_xxzz_0_xxxyzz_1,   \
                             g_0_xxzz_0_xxxzzz_0,   \
                             g_0_xxzz_0_xxxzzz_1,   \
                             g_0_xxzz_0_xxyyyz_0,   \
                             g_0_xxzz_0_xxyyyz_1,   \
                             g_0_xxzz_0_xxyyzz_0,   \
                             g_0_xxzz_0_xxyyzz_1,   \
                             g_0_xxzz_0_xxyzzz_0,   \
                             g_0_xxzz_0_xxyzzz_1,   \
                             g_0_xxzz_0_xxzzzz_0,   \
                             g_0_xxzz_0_xxzzzz_1,   \
                             g_0_xxzz_0_xyyyyz_0,   \
                             g_0_xxzz_0_xyyyyz_1,   \
                             g_0_xxzz_0_xyyyzz_0,   \
                             g_0_xxzz_0_xyyyzz_1,   \
                             g_0_xxzz_0_xyyzzz_0,   \
                             g_0_xxzz_0_xyyzzz_1,   \
                             g_0_xxzz_0_xyzzzz_0,   \
                             g_0_xxzz_0_xyzzzz_1,   \
                             g_0_xxzz_0_xzzzzz_0,   \
                             g_0_xxzz_0_xzzzzz_1,   \
                             g_0_xxzz_0_yyyyyy_0,   \
                             g_0_xxzz_0_yyyyyy_1,   \
                             g_0_xxzz_0_yyyyyz_0,   \
                             g_0_xxzz_0_yyyyyz_1,   \
                             g_0_xxzz_0_yyyyzz_0,   \
                             g_0_xxzz_0_yyyyzz_1,   \
                             g_0_xxzz_0_yyyzzz_0,   \
                             g_0_xxzz_0_yyyzzz_1,   \
                             g_0_xxzz_0_yyzzzz_0,   \
                             g_0_xxzz_0_yyzzzz_1,   \
                             g_0_xxzz_0_yzzzzz_0,   \
                             g_0_xxzz_0_yzzzzz_1,   \
                             g_0_xxzz_0_zzzzzz_0,   \
                             g_0_xxzz_0_zzzzzz_1,   \
                             wp_x,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzz_0_xxxxxx_0[i] =
            g_0_xxxx_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxxx_0[i] * pb_z + g_0_xxxxz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxxy_0[i] =
            g_0_xxxx_0_xxxxxy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxxy_0[i] * pb_z + g_0_xxxxz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxxz_0[i] = 3.0 * g_0_xxzz_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxz_1[i] * fti_ab_0 +
                                   5.0 * g_0_xxxzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxz_0[i] * pb_x + g_0_xxxzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxyy_0[i] =
            g_0_xxxx_0_xxxxyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxyy_0[i] * pb_z + g_0_xxxxz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxyz_0[i] = 3.0 * g_0_xxzz_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxyz_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxxzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyz_0[i] * pb_x + g_0_xxxzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxzz_0[i] = 3.0 * g_0_xxzz_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxxzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxzz_0[i] * pb_x + g_0_xxxzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxyyy_0[i] =
            g_0_xxxx_0_xxxyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxyyy_0[i] * pb_z + g_0_xxxxz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxyyz_0[i] = 3.0 * g_0_xxzz_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyyz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxxzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyz_0[i] * pb_x + g_0_xxxzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxyzz_0[i] = 3.0 * g_0_xxzz_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxxzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyzz_0[i] * pb_x + g_0_xxxzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxzzz_0[i] = 3.0 * g_0_xxzz_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxxzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxzzz_0[i] * pb_x + g_0_xxxzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyyyy_0[i] =
            g_0_xxxx_0_xxyyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxyyyy_0[i] * pb_z + g_0_xxxxz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxyyyz_0[i] = 3.0 * g_0_xxzz_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyyz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyz_0[i] * pb_x + g_0_xxxzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyyzz_0[i] = 3.0 * g_0_xxzz_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyzz_0[i] * pb_x + g_0_xxxzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyzzz_0[i] = 3.0 * g_0_xxzz_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyzzz_0[i] * pb_x + g_0_xxxzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxzzzz_0[i] = 3.0 * g_0_xxzz_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxzzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxxzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxzzzz_0[i] * pb_x + g_0_xxxzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyyyy_0[i] =
            g_0_xxxx_0_xyyyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xyyyyy_0[i] * pb_z + g_0_xxxxz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xyyyyz_0[i] = 3.0 * g_0_xxzz_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyyz_1[i] * fti_ab_0 +
                                   g_0_xxxzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyz_0[i] * pb_x + g_0_xxxzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyyzz_0[i] = 3.0 * g_0_xxzz_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyzz_1[i] * fti_ab_0 +
                                   g_0_xxxzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyzz_0[i] * pb_x + g_0_xxxzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyzzz_0[i] = 3.0 * g_0_xxzz_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyzzz_1[i] * fti_ab_0 +
                                   g_0_xxxzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyzzz_0[i] * pb_x + g_0_xxxzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyzzzz_0[i] = 3.0 * g_0_xxzz_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyzzzz_1[i] * fti_ab_0 +
                                   g_0_xxxzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyzzzz_0[i] * pb_x + g_0_xxxzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xzzzzz_0[i] = 3.0 * g_0_xxzz_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xzzzzz_1[i] * fti_ab_0 +
                                   g_0_xxxzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xzzzzz_0[i] * pb_x + g_0_xxxzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyyy_0[i] = 3.0 * g_0_xxzz_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyyyy_0[i] * pb_x +
                                   g_0_xxxzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyyz_0[i] = 3.0 * g_0_xxzz_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyyyz_0[i] * pb_x +
                                   g_0_xxxzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyzz_0[i] = 3.0 * g_0_xxzz_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyyzz_0[i] * pb_x +
                                   g_0_xxxzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyzzz_0[i] = 3.0 * g_0_xxzz_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyzzz_0[i] * pb_x +
                                   g_0_xxxzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyzzzz_0[i] = 3.0 * g_0_xxzz_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyzzzz_0[i] * pb_x +
                                   g_0_xxxzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yzzzzz_0[i] = 3.0 * g_0_xxzz_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yzzzzz_0[i] * pb_x +
                                   g_0_xxxzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_zzzzzz_0[i] = 3.0 * g_0_xxzz_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_zzzzzz_0[i] * pb_x +
                                   g_0_xxxzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 168-196 components of targeted buffer : SISI

    auto g_0_xxxyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 168);

    auto g_0_xxxyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 169);

    auto g_0_xxxyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 170);

    auto g_0_xxxyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 171);

    auto g_0_xxxyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 172);

    auto g_0_xxxyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 173);

    auto g_0_xxxyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 174);

    auto g_0_xxxyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 175);

    auto g_0_xxxyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 176);

    auto g_0_xxxyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 177);

    auto g_0_xxxyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 178);

    auto g_0_xxxyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 179);

    auto g_0_xxxyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 180);

    auto g_0_xxxyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 181);

    auto g_0_xxxyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 182);

    auto g_0_xxxyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 183);

    auto g_0_xxxyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 184);

    auto g_0_xxxyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 185);

    auto g_0_xxxyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 186);

    auto g_0_xxxyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 187);

    auto g_0_xxxyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 188);

    auto g_0_xxxyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 189);

    auto g_0_xxxyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 190);

    auto g_0_xxxyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 191);

    auto g_0_xxxyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 192);

    auto g_0_xxxyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 193);

    auto g_0_xxxyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 194);

    auto g_0_xxxyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 195);

#pragma omp simd aligned(g_0_xxxy_0_xxxxxx_0,       \
                             g_0_xxxy_0_xxxxxx_1,   \
                             g_0_xxxy_0_xxxxxz_0,   \
                             g_0_xxxy_0_xxxxxz_1,   \
                             g_0_xxxy_0_xxxxzz_0,   \
                             g_0_xxxy_0_xxxxzz_1,   \
                             g_0_xxxy_0_xxxzzz_0,   \
                             g_0_xxxy_0_xxxzzz_1,   \
                             g_0_xxxy_0_xxzzzz_0,   \
                             g_0_xxxy_0_xxzzzz_1,   \
                             g_0_xxxy_0_xzzzzz_0,   \
                             g_0_xxxy_0_xzzzzz_1,   \
                             g_0_xxxyy_0_xxxxxx_0,  \
                             g_0_xxxyy_0_xxxxxx_1,  \
                             g_0_xxxyy_0_xxxxxz_0,  \
                             g_0_xxxyy_0_xxxxxz_1,  \
                             g_0_xxxyy_0_xxxxzz_0,  \
                             g_0_xxxyy_0_xxxxzz_1,  \
                             g_0_xxxyy_0_xxxzzz_0,  \
                             g_0_xxxyy_0_xxxzzz_1,  \
                             g_0_xxxyy_0_xxzzzz_0,  \
                             g_0_xxxyy_0_xxzzzz_1,  \
                             g_0_xxxyy_0_xzzzzz_0,  \
                             g_0_xxxyy_0_xzzzzz_1,  \
                             g_0_xxxyyy_0_xxxxxx_0, \
                             g_0_xxxyyy_0_xxxxxy_0, \
                             g_0_xxxyyy_0_xxxxxz_0, \
                             g_0_xxxyyy_0_xxxxyy_0, \
                             g_0_xxxyyy_0_xxxxyz_0, \
                             g_0_xxxyyy_0_xxxxzz_0, \
                             g_0_xxxyyy_0_xxxyyy_0, \
                             g_0_xxxyyy_0_xxxyyz_0, \
                             g_0_xxxyyy_0_xxxyzz_0, \
                             g_0_xxxyyy_0_xxxzzz_0, \
                             g_0_xxxyyy_0_xxyyyy_0, \
                             g_0_xxxyyy_0_xxyyyz_0, \
                             g_0_xxxyyy_0_xxyyzz_0, \
                             g_0_xxxyyy_0_xxyzzz_0, \
                             g_0_xxxyyy_0_xxzzzz_0, \
                             g_0_xxxyyy_0_xyyyyy_0, \
                             g_0_xxxyyy_0_xyyyyz_0, \
                             g_0_xxxyyy_0_xyyyzz_0, \
                             g_0_xxxyyy_0_xyyzzz_0, \
                             g_0_xxxyyy_0_xyzzzz_0, \
                             g_0_xxxyyy_0_xzzzzz_0, \
                             g_0_xxxyyy_0_yyyyyy_0, \
                             g_0_xxxyyy_0_yyyyyz_0, \
                             g_0_xxxyyy_0_yyyyzz_0, \
                             g_0_xxxyyy_0_yyyzzz_0, \
                             g_0_xxxyyy_0_yyzzzz_0, \
                             g_0_xxxyyy_0_yzzzzz_0, \
                             g_0_xxxyyy_0_zzzzzz_0, \
                             g_0_xxyyy_0_xxxxxy_0,  \
                             g_0_xxyyy_0_xxxxxy_1,  \
                             g_0_xxyyy_0_xxxxy_1,   \
                             g_0_xxyyy_0_xxxxyy_0,  \
                             g_0_xxyyy_0_xxxxyy_1,  \
                             g_0_xxyyy_0_xxxxyz_0,  \
                             g_0_xxyyy_0_xxxxyz_1,  \
                             g_0_xxyyy_0_xxxyy_1,   \
                             g_0_xxyyy_0_xxxyyy_0,  \
                             g_0_xxyyy_0_xxxyyy_1,  \
                             g_0_xxyyy_0_xxxyyz_0,  \
                             g_0_xxyyy_0_xxxyyz_1,  \
                             g_0_xxyyy_0_xxxyz_1,   \
                             g_0_xxyyy_0_xxxyzz_0,  \
                             g_0_xxyyy_0_xxxyzz_1,  \
                             g_0_xxyyy_0_xxyyy_1,   \
                             g_0_xxyyy_0_xxyyyy_0,  \
                             g_0_xxyyy_0_xxyyyy_1,  \
                             g_0_xxyyy_0_xxyyyz_0,  \
                             g_0_xxyyy_0_xxyyyz_1,  \
                             g_0_xxyyy_0_xxyyz_1,   \
                             g_0_xxyyy_0_xxyyzz_0,  \
                             g_0_xxyyy_0_xxyyzz_1,  \
                             g_0_xxyyy_0_xxyzz_1,   \
                             g_0_xxyyy_0_xxyzzz_0,  \
                             g_0_xxyyy_0_xxyzzz_1,  \
                             g_0_xxyyy_0_xyyyy_1,   \
                             g_0_xxyyy_0_xyyyyy_0,  \
                             g_0_xxyyy_0_xyyyyy_1,  \
                             g_0_xxyyy_0_xyyyyz_0,  \
                             g_0_xxyyy_0_xyyyyz_1,  \
                             g_0_xxyyy_0_xyyyz_1,   \
                             g_0_xxyyy_0_xyyyzz_0,  \
                             g_0_xxyyy_0_xyyyzz_1,  \
                             g_0_xxyyy_0_xyyzz_1,   \
                             g_0_xxyyy_0_xyyzzz_0,  \
                             g_0_xxyyy_0_xyyzzz_1,  \
                             g_0_xxyyy_0_xyzzz_1,   \
                             g_0_xxyyy_0_xyzzzz_0,  \
                             g_0_xxyyy_0_xyzzzz_1,  \
                             g_0_xxyyy_0_yyyyy_1,   \
                             g_0_xxyyy_0_yyyyyy_0,  \
                             g_0_xxyyy_0_yyyyyy_1,  \
                             g_0_xxyyy_0_yyyyyz_0,  \
                             g_0_xxyyy_0_yyyyyz_1,  \
                             g_0_xxyyy_0_yyyyz_1,   \
                             g_0_xxyyy_0_yyyyzz_0,  \
                             g_0_xxyyy_0_yyyyzz_1,  \
                             g_0_xxyyy_0_yyyzz_1,   \
                             g_0_xxyyy_0_yyyzzz_0,  \
                             g_0_xxyyy_0_yyyzzz_1,  \
                             g_0_xxyyy_0_yyzzz_1,   \
                             g_0_xxyyy_0_yyzzzz_0,  \
                             g_0_xxyyy_0_yyzzzz_1,  \
                             g_0_xxyyy_0_yzzzz_1,   \
                             g_0_xxyyy_0_yzzzzz_0,  \
                             g_0_xxyyy_0_yzzzzz_1,  \
                             g_0_xxyyy_0_zzzzzz_0,  \
                             g_0_xxyyy_0_zzzzzz_1,  \
                             g_0_xyyy_0_xxxxxy_0,   \
                             g_0_xyyy_0_xxxxxy_1,   \
                             g_0_xyyy_0_xxxxyy_0,   \
                             g_0_xyyy_0_xxxxyy_1,   \
                             g_0_xyyy_0_xxxxyz_0,   \
                             g_0_xyyy_0_xxxxyz_1,   \
                             g_0_xyyy_0_xxxyyy_0,   \
                             g_0_xyyy_0_xxxyyy_1,   \
                             g_0_xyyy_0_xxxyyz_0,   \
                             g_0_xyyy_0_xxxyyz_1,   \
                             g_0_xyyy_0_xxxyzz_0,   \
                             g_0_xyyy_0_xxxyzz_1,   \
                             g_0_xyyy_0_xxyyyy_0,   \
                             g_0_xyyy_0_xxyyyy_1,   \
                             g_0_xyyy_0_xxyyyz_0,   \
                             g_0_xyyy_0_xxyyyz_1,   \
                             g_0_xyyy_0_xxyyzz_0,   \
                             g_0_xyyy_0_xxyyzz_1,   \
                             g_0_xyyy_0_xxyzzz_0,   \
                             g_0_xyyy_0_xxyzzz_1,   \
                             g_0_xyyy_0_xyyyyy_0,   \
                             g_0_xyyy_0_xyyyyy_1,   \
                             g_0_xyyy_0_xyyyyz_0,   \
                             g_0_xyyy_0_xyyyyz_1,   \
                             g_0_xyyy_0_xyyyzz_0,   \
                             g_0_xyyy_0_xyyyzz_1,   \
                             g_0_xyyy_0_xyyzzz_0,   \
                             g_0_xyyy_0_xyyzzz_1,   \
                             g_0_xyyy_0_xyzzzz_0,   \
                             g_0_xyyy_0_xyzzzz_1,   \
                             g_0_xyyy_0_yyyyyy_0,   \
                             g_0_xyyy_0_yyyyyy_1,   \
                             g_0_xyyy_0_yyyyyz_0,   \
                             g_0_xyyy_0_yyyyyz_1,   \
                             g_0_xyyy_0_yyyyzz_0,   \
                             g_0_xyyy_0_yyyyzz_1,   \
                             g_0_xyyy_0_yyyzzz_0,   \
                             g_0_xyyy_0_yyyzzz_1,   \
                             g_0_xyyy_0_yyzzzz_0,   \
                             g_0_xyyy_0_yyzzzz_1,   \
                             g_0_xyyy_0_yzzzzz_0,   \
                             g_0_xyyy_0_yzzzzz_1,   \
                             g_0_xyyy_0_zzzzzz_0,   \
                             g_0_xyyy_0_zzzzzz_1,   \
                             wp_x,                  \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyy_0_xxxxxx_0[i] = 2.0 * g_0_xxxy_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxxxx_0[i] * pb_y +
                                   g_0_xxxyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxxxy_0[i] = 2.0 * g_0_xyyy_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxxy_1[i] * fti_ab_0 +
                                   5.0 * g_0_xxyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxy_0[i] * pb_x + g_0_xxyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxxz_0[i] = 2.0 * g_0_xxxy_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxxxz_0[i] * pb_y +
                                   g_0_xxxyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxxyy_0[i] = 2.0 * g_0_xyyy_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxyy_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyy_0[i] * pb_x + g_0_xxyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxyz_0[i] = 2.0 * g_0_xyyy_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxyz_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyz_0[i] * pb_x + g_0_xxyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxzz_0[i] = 2.0 * g_0_xxxy_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxxzz_0[i] * pb_y +
                                   g_0_xxxyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxyyy_0[i] = 2.0 * g_0_xyyy_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyyy_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyy_0[i] * pb_x + g_0_xxyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxyyz_0[i] = 2.0 * g_0_xyyy_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyyz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyz_0[i] * pb_x + g_0_xxyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxyzz_0[i] = 2.0 * g_0_xyyy_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyzz_0[i] * pb_x + g_0_xxyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxzzz_0[i] = 2.0 * g_0_xxxy_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxzzz_0[i] * pb_y +
                                   g_0_xxxyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxyyyy_0[i] = 2.0 * g_0_xyyy_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyyy_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyy_0[i] * pb_x + g_0_xxyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyyyz_0[i] = 2.0 * g_0_xyyy_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyyz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyz_0[i] * pb_x + g_0_xxyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyyzz_0[i] = 2.0 * g_0_xyyy_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyzz_0[i] * pb_x + g_0_xxyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyzzz_0[i] = 2.0 * g_0_xyyy_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyzzz_0[i] * pb_x + g_0_xxyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxzzzz_0[i] = 2.0 * g_0_xxxy_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxzzzz_0[i] * pb_y +
                                   g_0_xxxyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xyyyyy_0[i] = 2.0 * g_0_xyyy_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyyy_1[i] * fti_ab_0 +
                                   g_0_xxyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyy_0[i] * pb_x + g_0_xxyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyyyz_0[i] = 2.0 * g_0_xyyy_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyyz_1[i] * fti_ab_0 +
                                   g_0_xxyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyz_0[i] * pb_x + g_0_xxyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyyzz_0[i] = 2.0 * g_0_xyyy_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyzz_1[i] * fti_ab_0 +
                                   g_0_xxyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyzz_0[i] * pb_x + g_0_xxyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyzzz_0[i] = 2.0 * g_0_xyyy_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyzzz_1[i] * fti_ab_0 +
                                   g_0_xxyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyzzz_0[i] * pb_x + g_0_xxyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyzzzz_0[i] = 2.0 * g_0_xyyy_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyzzzz_1[i] * fti_ab_0 +
                                   g_0_xxyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyzzzz_0[i] * pb_x + g_0_xxyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xzzzzz_0[i] = 2.0 * g_0_xxxy_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xzzzzz_0[i] * pb_y +
                                   g_0_xxxyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_yyyyyy_0[i] = 2.0 * g_0_xyyy_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyyy_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyyyy_0[i] * pb_x +
                                   g_0_xxyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyyyz_0[i] = 2.0 * g_0_xyyy_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyyz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyyyz_0[i] * pb_x +
                                   g_0_xxyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyyzz_0[i] = 2.0 * g_0_xyyy_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyyzz_0[i] * pb_x +
                                   g_0_xxyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyzzz_0[i] = 2.0 * g_0_xyyy_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyzzz_0[i] * pb_x +
                                   g_0_xxyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyzzzz_0[i] = 2.0 * g_0_xyyy_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyzzzz_0[i] * pb_x +
                                   g_0_xxyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yzzzzz_0[i] = 2.0 * g_0_xyyy_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yzzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yzzzzz_0[i] * pb_x +
                                   g_0_xxyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_zzzzzz_0[i] = 2.0 * g_0_xyyy_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_zzzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_zzzzzz_0[i] * pb_x +
                                   g_0_xxyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 196-224 components of targeted buffer : SISI

    auto g_0_xxxyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 196);

    auto g_0_xxxyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 197);

    auto g_0_xxxyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 198);

    auto g_0_xxxyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 199);

    auto g_0_xxxyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 200);

    auto g_0_xxxyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 201);

    auto g_0_xxxyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 202);

    auto g_0_xxxyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 203);

    auto g_0_xxxyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 204);

    auto g_0_xxxyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 205);

    auto g_0_xxxyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 206);

    auto g_0_xxxyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 207);

    auto g_0_xxxyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 208);

    auto g_0_xxxyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 209);

    auto g_0_xxxyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 210);

    auto g_0_xxxyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 211);

    auto g_0_xxxyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 212);

    auto g_0_xxxyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 213);

    auto g_0_xxxyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 214);

    auto g_0_xxxyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 215);

    auto g_0_xxxyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 216);

    auto g_0_xxxyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 217);

    auto g_0_xxxyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 218);

    auto g_0_xxxyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 219);

    auto g_0_xxxyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 220);

    auto g_0_xxxyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 221);

    auto g_0_xxxyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 222);

    auto g_0_xxxyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 223);

#pragma omp simd aligned(g_0_xxxyy_0_xxxxx_1,       \
                             g_0_xxxyy_0_xxxxxx_0,  \
                             g_0_xxxyy_0_xxxxxx_1,  \
                             g_0_xxxyy_0_xxxxxy_0,  \
                             g_0_xxxyy_0_xxxxxy_1,  \
                             g_0_xxxyy_0_xxxxxz_0,  \
                             g_0_xxxyy_0_xxxxxz_1,  \
                             g_0_xxxyy_0_xxxxy_1,   \
                             g_0_xxxyy_0_xxxxyy_0,  \
                             g_0_xxxyy_0_xxxxyy_1,  \
                             g_0_xxxyy_0_xxxxyz_0,  \
                             g_0_xxxyy_0_xxxxyz_1,  \
                             g_0_xxxyy_0_xxxxz_1,   \
                             g_0_xxxyy_0_xxxxzz_0,  \
                             g_0_xxxyy_0_xxxxzz_1,  \
                             g_0_xxxyy_0_xxxyy_1,   \
                             g_0_xxxyy_0_xxxyyy_0,  \
                             g_0_xxxyy_0_xxxyyy_1,  \
                             g_0_xxxyy_0_xxxyyz_0,  \
                             g_0_xxxyy_0_xxxyyz_1,  \
                             g_0_xxxyy_0_xxxyz_1,   \
                             g_0_xxxyy_0_xxxyzz_0,  \
                             g_0_xxxyy_0_xxxyzz_1,  \
                             g_0_xxxyy_0_xxxzz_1,   \
                             g_0_xxxyy_0_xxxzzz_0,  \
                             g_0_xxxyy_0_xxxzzz_1,  \
                             g_0_xxxyy_0_xxyyy_1,   \
                             g_0_xxxyy_0_xxyyyy_0,  \
                             g_0_xxxyy_0_xxyyyy_1,  \
                             g_0_xxxyy_0_xxyyyz_0,  \
                             g_0_xxxyy_0_xxyyyz_1,  \
                             g_0_xxxyy_0_xxyyz_1,   \
                             g_0_xxxyy_0_xxyyzz_0,  \
                             g_0_xxxyy_0_xxyyzz_1,  \
                             g_0_xxxyy_0_xxyzz_1,   \
                             g_0_xxxyy_0_xxyzzz_0,  \
                             g_0_xxxyy_0_xxyzzz_1,  \
                             g_0_xxxyy_0_xxzzz_1,   \
                             g_0_xxxyy_0_xxzzzz_0,  \
                             g_0_xxxyy_0_xxzzzz_1,  \
                             g_0_xxxyy_0_xyyyy_1,   \
                             g_0_xxxyy_0_xyyyyy_0,  \
                             g_0_xxxyy_0_xyyyyy_1,  \
                             g_0_xxxyy_0_xyyyyz_0,  \
                             g_0_xxxyy_0_xyyyyz_1,  \
                             g_0_xxxyy_0_xyyyz_1,   \
                             g_0_xxxyy_0_xyyyzz_0,  \
                             g_0_xxxyy_0_xyyyzz_1,  \
                             g_0_xxxyy_0_xyyzz_1,   \
                             g_0_xxxyy_0_xyyzzz_0,  \
                             g_0_xxxyy_0_xyyzzz_1,  \
                             g_0_xxxyy_0_xyzzz_1,   \
                             g_0_xxxyy_0_xyzzzz_0,  \
                             g_0_xxxyy_0_xyzzzz_1,  \
                             g_0_xxxyy_0_xzzzz_1,   \
                             g_0_xxxyy_0_xzzzzz_0,  \
                             g_0_xxxyy_0_xzzzzz_1,  \
                             g_0_xxxyy_0_yyyyy_1,   \
                             g_0_xxxyy_0_yyyyyy_0,  \
                             g_0_xxxyy_0_yyyyyy_1,  \
                             g_0_xxxyy_0_yyyyyz_0,  \
                             g_0_xxxyy_0_yyyyyz_1,  \
                             g_0_xxxyy_0_yyyyz_1,   \
                             g_0_xxxyy_0_yyyyzz_0,  \
                             g_0_xxxyy_0_yyyyzz_1,  \
                             g_0_xxxyy_0_yyyzz_1,   \
                             g_0_xxxyy_0_yyyzzz_0,  \
                             g_0_xxxyy_0_yyyzzz_1,  \
                             g_0_xxxyy_0_yyzzz_1,   \
                             g_0_xxxyy_0_yyzzzz_0,  \
                             g_0_xxxyy_0_yyzzzz_1,  \
                             g_0_xxxyy_0_yzzzz_1,   \
                             g_0_xxxyy_0_yzzzzz_0,  \
                             g_0_xxxyy_0_yzzzzz_1,  \
                             g_0_xxxyy_0_zzzzz_1,   \
                             g_0_xxxyy_0_zzzzzz_0,  \
                             g_0_xxxyy_0_zzzzzz_1,  \
                             g_0_xxxyyz_0_xxxxxx_0, \
                             g_0_xxxyyz_0_xxxxxy_0, \
                             g_0_xxxyyz_0_xxxxxz_0, \
                             g_0_xxxyyz_0_xxxxyy_0, \
                             g_0_xxxyyz_0_xxxxyz_0, \
                             g_0_xxxyyz_0_xxxxzz_0, \
                             g_0_xxxyyz_0_xxxyyy_0, \
                             g_0_xxxyyz_0_xxxyyz_0, \
                             g_0_xxxyyz_0_xxxyzz_0, \
                             g_0_xxxyyz_0_xxxzzz_0, \
                             g_0_xxxyyz_0_xxyyyy_0, \
                             g_0_xxxyyz_0_xxyyyz_0, \
                             g_0_xxxyyz_0_xxyyzz_0, \
                             g_0_xxxyyz_0_xxyzzz_0, \
                             g_0_xxxyyz_0_xxzzzz_0, \
                             g_0_xxxyyz_0_xyyyyy_0, \
                             g_0_xxxyyz_0_xyyyyz_0, \
                             g_0_xxxyyz_0_xyyyzz_0, \
                             g_0_xxxyyz_0_xyyzzz_0, \
                             g_0_xxxyyz_0_xyzzzz_0, \
                             g_0_xxxyyz_0_xzzzzz_0, \
                             g_0_xxxyyz_0_yyyyyy_0, \
                             g_0_xxxyyz_0_yyyyyz_0, \
                             g_0_xxxyyz_0_yyyyzz_0, \
                             g_0_xxxyyz_0_yyyzzz_0, \
                             g_0_xxxyyz_0_yyzzzz_0, \
                             g_0_xxxyyz_0_yzzzzz_0, \
                             g_0_xxxyyz_0_zzzzzz_0, \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyz_0_xxxxxx_0[i] = g_0_xxxyy_0_xxxxxx_0[i] * pb_z + g_0_xxxyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxy_0[i] = g_0_xxxyy_0_xxxxxy_0[i] * pb_z + g_0_xxxyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxz_0[i] = g_0_xxxyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxz_0[i] * pb_z + g_0_xxxyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxyy_0[i] = g_0_xxxyy_0_xxxxyy_0[i] * pb_z + g_0_xxxyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxyz_0[i] = g_0_xxxyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyz_0[i] * pb_z + g_0_xxxyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxzz_0[i] = 2.0 * g_0_xxxyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxzz_0[i] * pb_z + g_0_xxxyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyyy_0[i] = g_0_xxxyy_0_xxxyyy_0[i] * pb_z + g_0_xxxyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyyz_0[i] = g_0_xxxyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyz_0[i] * pb_z + g_0_xxxyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyzz_0[i] = 2.0 * g_0_xxxyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyzz_0[i] * pb_z + g_0_xxxyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxzzz_0[i] = 3.0 * g_0_xxxyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxzzz_0[i] * pb_z + g_0_xxxyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyyy_0[i] = g_0_xxxyy_0_xxyyyy_0[i] * pb_z + g_0_xxxyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyyz_0[i] = g_0_xxxyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyz_0[i] * pb_z + g_0_xxxyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyzz_0[i] = 2.0 * g_0_xxxyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyzz_0[i] * pb_z + g_0_xxxyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyzzz_0[i] = 3.0 * g_0_xxxyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyzzz_0[i] * pb_z + g_0_xxxyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxzzzz_0[i] = 4.0 * g_0_xxxyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxzzzz_0[i] * pb_z + g_0_xxxyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyyy_0[i] = g_0_xxxyy_0_xyyyyy_0[i] * pb_z + g_0_xxxyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyyz_0[i] = g_0_xxxyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyz_0[i] * pb_z + g_0_xxxyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyzz_0[i] = 2.0 * g_0_xxxyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyzz_0[i] * pb_z + g_0_xxxyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyzzz_0[i] = 3.0 * g_0_xxxyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyzzz_0[i] * pb_z + g_0_xxxyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyzzzz_0[i] = 4.0 * g_0_xxxyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyzzzz_0[i] * pb_z + g_0_xxxyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xzzzzz_0[i] = 5.0 * g_0_xxxyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xzzzzz_0[i] * pb_z + g_0_xxxyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyyy_0[i] = g_0_xxxyy_0_yyyyyy_0[i] * pb_z + g_0_xxxyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyyz_0[i] = g_0_xxxyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyyyz_0[i] * pb_z + g_0_xxxyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyzz_0[i] = 2.0 * g_0_xxxyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyyzz_0[i] * pb_z + g_0_xxxyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyzzz_0[i] = 3.0 * g_0_xxxyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyzzz_0[i] * pb_z + g_0_xxxyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyzzzz_0[i] = 4.0 * g_0_xxxyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyzzzz_0[i] * pb_z + g_0_xxxyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yzzzzz_0[i] = 5.0 * g_0_xxxyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yzzzzz_0[i] * pb_z + g_0_xxxyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_zzzzzz_0[i] = 6.0 * g_0_xxxyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_zzzzzz_0[i] * pb_z + g_0_xxxyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 224-252 components of targeted buffer : SISI

    auto g_0_xxxyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 224);

    auto g_0_xxxyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 225);

    auto g_0_xxxyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 226);

    auto g_0_xxxyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 227);

    auto g_0_xxxyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 228);

    auto g_0_xxxyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 229);

    auto g_0_xxxyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 230);

    auto g_0_xxxyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 231);

    auto g_0_xxxyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 232);

    auto g_0_xxxyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 233);

    auto g_0_xxxyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 234);

    auto g_0_xxxyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 235);

    auto g_0_xxxyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 236);

    auto g_0_xxxyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 237);

    auto g_0_xxxyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 238);

    auto g_0_xxxyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 239);

    auto g_0_xxxyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 240);

    auto g_0_xxxyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 241);

    auto g_0_xxxyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 242);

    auto g_0_xxxyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 243);

    auto g_0_xxxyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 244);

    auto g_0_xxxyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 245);

    auto g_0_xxxyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 246);

    auto g_0_xxxyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 247);

    auto g_0_xxxyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 248);

    auto g_0_xxxyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 249);

    auto g_0_xxxyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 250);

    auto g_0_xxxyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 251);

#pragma omp simd aligned(g_0_xxxyzz_0_xxxxxx_0,     \
                             g_0_xxxyzz_0_xxxxxy_0, \
                             g_0_xxxyzz_0_xxxxxz_0, \
                             g_0_xxxyzz_0_xxxxyy_0, \
                             g_0_xxxyzz_0_xxxxyz_0, \
                             g_0_xxxyzz_0_xxxxzz_0, \
                             g_0_xxxyzz_0_xxxyyy_0, \
                             g_0_xxxyzz_0_xxxyyz_0, \
                             g_0_xxxyzz_0_xxxyzz_0, \
                             g_0_xxxyzz_0_xxxzzz_0, \
                             g_0_xxxyzz_0_xxyyyy_0, \
                             g_0_xxxyzz_0_xxyyyz_0, \
                             g_0_xxxyzz_0_xxyyzz_0, \
                             g_0_xxxyzz_0_xxyzzz_0, \
                             g_0_xxxyzz_0_xxzzzz_0, \
                             g_0_xxxyzz_0_xyyyyy_0, \
                             g_0_xxxyzz_0_xyyyyz_0, \
                             g_0_xxxyzz_0_xyyyzz_0, \
                             g_0_xxxyzz_0_xyyzzz_0, \
                             g_0_xxxyzz_0_xyzzzz_0, \
                             g_0_xxxyzz_0_xzzzzz_0, \
                             g_0_xxxyzz_0_yyyyyy_0, \
                             g_0_xxxyzz_0_yyyyyz_0, \
                             g_0_xxxyzz_0_yyyyzz_0, \
                             g_0_xxxyzz_0_yyyzzz_0, \
                             g_0_xxxyzz_0_yyzzzz_0, \
                             g_0_xxxyzz_0_yzzzzz_0, \
                             g_0_xxxyzz_0_zzzzzz_0, \
                             g_0_xxxzz_0_xxxxx_1,   \
                             g_0_xxxzz_0_xxxxxx_0,  \
                             g_0_xxxzz_0_xxxxxx_1,  \
                             g_0_xxxzz_0_xxxxxy_0,  \
                             g_0_xxxzz_0_xxxxxy_1,  \
                             g_0_xxxzz_0_xxxxxz_0,  \
                             g_0_xxxzz_0_xxxxxz_1,  \
                             g_0_xxxzz_0_xxxxy_1,   \
                             g_0_xxxzz_0_xxxxyy_0,  \
                             g_0_xxxzz_0_xxxxyy_1,  \
                             g_0_xxxzz_0_xxxxyz_0,  \
                             g_0_xxxzz_0_xxxxyz_1,  \
                             g_0_xxxzz_0_xxxxz_1,   \
                             g_0_xxxzz_0_xxxxzz_0,  \
                             g_0_xxxzz_0_xxxxzz_1,  \
                             g_0_xxxzz_0_xxxyy_1,   \
                             g_0_xxxzz_0_xxxyyy_0,  \
                             g_0_xxxzz_0_xxxyyy_1,  \
                             g_0_xxxzz_0_xxxyyz_0,  \
                             g_0_xxxzz_0_xxxyyz_1,  \
                             g_0_xxxzz_0_xxxyz_1,   \
                             g_0_xxxzz_0_xxxyzz_0,  \
                             g_0_xxxzz_0_xxxyzz_1,  \
                             g_0_xxxzz_0_xxxzz_1,   \
                             g_0_xxxzz_0_xxxzzz_0,  \
                             g_0_xxxzz_0_xxxzzz_1,  \
                             g_0_xxxzz_0_xxyyy_1,   \
                             g_0_xxxzz_0_xxyyyy_0,  \
                             g_0_xxxzz_0_xxyyyy_1,  \
                             g_0_xxxzz_0_xxyyyz_0,  \
                             g_0_xxxzz_0_xxyyyz_1,  \
                             g_0_xxxzz_0_xxyyz_1,   \
                             g_0_xxxzz_0_xxyyzz_0,  \
                             g_0_xxxzz_0_xxyyzz_1,  \
                             g_0_xxxzz_0_xxyzz_1,   \
                             g_0_xxxzz_0_xxyzzz_0,  \
                             g_0_xxxzz_0_xxyzzz_1,  \
                             g_0_xxxzz_0_xxzzz_1,   \
                             g_0_xxxzz_0_xxzzzz_0,  \
                             g_0_xxxzz_0_xxzzzz_1,  \
                             g_0_xxxzz_0_xyyyy_1,   \
                             g_0_xxxzz_0_xyyyyy_0,  \
                             g_0_xxxzz_0_xyyyyy_1,  \
                             g_0_xxxzz_0_xyyyyz_0,  \
                             g_0_xxxzz_0_xyyyyz_1,  \
                             g_0_xxxzz_0_xyyyz_1,   \
                             g_0_xxxzz_0_xyyyzz_0,  \
                             g_0_xxxzz_0_xyyyzz_1,  \
                             g_0_xxxzz_0_xyyzz_1,   \
                             g_0_xxxzz_0_xyyzzz_0,  \
                             g_0_xxxzz_0_xyyzzz_1,  \
                             g_0_xxxzz_0_xyzzz_1,   \
                             g_0_xxxzz_0_xyzzzz_0,  \
                             g_0_xxxzz_0_xyzzzz_1,  \
                             g_0_xxxzz_0_xzzzz_1,   \
                             g_0_xxxzz_0_xzzzzz_0,  \
                             g_0_xxxzz_0_xzzzzz_1,  \
                             g_0_xxxzz_0_yyyyy_1,   \
                             g_0_xxxzz_0_yyyyyy_0,  \
                             g_0_xxxzz_0_yyyyyy_1,  \
                             g_0_xxxzz_0_yyyyyz_0,  \
                             g_0_xxxzz_0_yyyyyz_1,  \
                             g_0_xxxzz_0_yyyyz_1,   \
                             g_0_xxxzz_0_yyyyzz_0,  \
                             g_0_xxxzz_0_yyyyzz_1,  \
                             g_0_xxxzz_0_yyyzz_1,   \
                             g_0_xxxzz_0_yyyzzz_0,  \
                             g_0_xxxzz_0_yyyzzz_1,  \
                             g_0_xxxzz_0_yyzzz_1,   \
                             g_0_xxxzz_0_yyzzzz_0,  \
                             g_0_xxxzz_0_yyzzzz_1,  \
                             g_0_xxxzz_0_yzzzz_1,   \
                             g_0_xxxzz_0_yzzzzz_0,  \
                             g_0_xxxzz_0_yzzzzz_1,  \
                             g_0_xxxzz_0_zzzzz_1,   \
                             g_0_xxxzz_0_zzzzzz_0,  \
                             g_0_xxxzz_0_zzzzzz_1,  \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzz_0_xxxxxx_0[i] = g_0_xxxzz_0_xxxxxx_0[i] * pb_y + g_0_xxxzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxy_0[i] = g_0_xxxzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxy_0[i] * pb_y + g_0_xxxzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxz_0[i] = g_0_xxxzz_0_xxxxxz_0[i] * pb_y + g_0_xxxzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxyy_0[i] = 2.0 * g_0_xxxzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyy_0[i] * pb_y + g_0_xxxzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxyz_0[i] = g_0_xxxzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyz_0[i] * pb_y + g_0_xxxzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxzz_0[i] = g_0_xxxzz_0_xxxxzz_0[i] * pb_y + g_0_xxxzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyyy_0[i] = 3.0 * g_0_xxxzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyy_0[i] * pb_y + g_0_xxxzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyyz_0[i] = 2.0 * g_0_xxxzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyz_0[i] * pb_y + g_0_xxxzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyzz_0[i] = g_0_xxxzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyzz_0[i] * pb_y + g_0_xxxzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxzzz_0[i] = g_0_xxxzz_0_xxxzzz_0[i] * pb_y + g_0_xxxzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyyy_0[i] = 4.0 * g_0_xxxzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyy_0[i] * pb_y + g_0_xxxzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyyz_0[i] = 3.0 * g_0_xxxzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyz_0[i] * pb_y + g_0_xxxzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyzz_0[i] = 2.0 * g_0_xxxzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyzz_0[i] * pb_y + g_0_xxxzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyzzz_0[i] = g_0_xxxzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyzzz_0[i] * pb_y + g_0_xxxzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxzzzz_0[i] = g_0_xxxzz_0_xxzzzz_0[i] * pb_y + g_0_xxxzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyyy_0[i] = 5.0 * g_0_xxxzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyy_0[i] * pb_y + g_0_xxxzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyyz_0[i] = 4.0 * g_0_xxxzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyz_0[i] * pb_y + g_0_xxxzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyzz_0[i] = 3.0 * g_0_xxxzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyzz_0[i] * pb_y + g_0_xxxzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyzzz_0[i] = 2.0 * g_0_xxxzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyzzz_0[i] * pb_y + g_0_xxxzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyzzzz_0[i] = g_0_xxxzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyzzzz_0[i] * pb_y + g_0_xxxzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xzzzzz_0[i] = g_0_xxxzz_0_xzzzzz_0[i] * pb_y + g_0_xxxzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyyy_0[i] = 6.0 * g_0_xxxzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyyy_0[i] * pb_y + g_0_xxxzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyyz_0[i] = 5.0 * g_0_xxxzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyyz_0[i] * pb_y + g_0_xxxzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyzz_0[i] = 4.0 * g_0_xxxzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyzz_0[i] * pb_y + g_0_xxxzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyzzz_0[i] = 3.0 * g_0_xxxzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyzzz_0[i] * pb_y + g_0_xxxzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyzzzz_0[i] = 2.0 * g_0_xxxzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyzzzz_0[i] * pb_y + g_0_xxxzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yzzzzz_0[i] = g_0_xxxzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yzzzzz_0[i] * pb_y + g_0_xxxzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_zzzzzz_0[i] = g_0_xxxzz_0_zzzzzz_0[i] * pb_y + g_0_xxxzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 252-280 components of targeted buffer : SISI

    auto g_0_xxxzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 252);

    auto g_0_xxxzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 253);

    auto g_0_xxxzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 254);

    auto g_0_xxxzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 255);

    auto g_0_xxxzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 256);

    auto g_0_xxxzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 257);

    auto g_0_xxxzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 258);

    auto g_0_xxxzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 259);

    auto g_0_xxxzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 260);

    auto g_0_xxxzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 261);

    auto g_0_xxxzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 262);

    auto g_0_xxxzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 263);

    auto g_0_xxxzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 264);

    auto g_0_xxxzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 265);

    auto g_0_xxxzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 266);

    auto g_0_xxxzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 267);

    auto g_0_xxxzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 268);

    auto g_0_xxxzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 269);

    auto g_0_xxxzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 270);

    auto g_0_xxxzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 271);

    auto g_0_xxxzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 272);

    auto g_0_xxxzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 273);

    auto g_0_xxxzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 274);

    auto g_0_xxxzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 275);

    auto g_0_xxxzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 276);

    auto g_0_xxxzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 277);

    auto g_0_xxxzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 278);

    auto g_0_xxxzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 279);

#pragma omp simd aligned(g_0_xxxz_0_xxxxxx_0,       \
                             g_0_xxxz_0_xxxxxx_1,   \
                             g_0_xxxz_0_xxxxxy_0,   \
                             g_0_xxxz_0_xxxxxy_1,   \
                             g_0_xxxz_0_xxxxyy_0,   \
                             g_0_xxxz_0_xxxxyy_1,   \
                             g_0_xxxz_0_xxxyyy_0,   \
                             g_0_xxxz_0_xxxyyy_1,   \
                             g_0_xxxz_0_xxyyyy_0,   \
                             g_0_xxxz_0_xxyyyy_1,   \
                             g_0_xxxz_0_xyyyyy_0,   \
                             g_0_xxxz_0_xyyyyy_1,   \
                             g_0_xxxzz_0_xxxxxx_0,  \
                             g_0_xxxzz_0_xxxxxx_1,  \
                             g_0_xxxzz_0_xxxxxy_0,  \
                             g_0_xxxzz_0_xxxxxy_1,  \
                             g_0_xxxzz_0_xxxxyy_0,  \
                             g_0_xxxzz_0_xxxxyy_1,  \
                             g_0_xxxzz_0_xxxyyy_0,  \
                             g_0_xxxzz_0_xxxyyy_1,  \
                             g_0_xxxzz_0_xxyyyy_0,  \
                             g_0_xxxzz_0_xxyyyy_1,  \
                             g_0_xxxzz_0_xyyyyy_0,  \
                             g_0_xxxzz_0_xyyyyy_1,  \
                             g_0_xxxzzz_0_xxxxxx_0, \
                             g_0_xxxzzz_0_xxxxxy_0, \
                             g_0_xxxzzz_0_xxxxxz_0, \
                             g_0_xxxzzz_0_xxxxyy_0, \
                             g_0_xxxzzz_0_xxxxyz_0, \
                             g_0_xxxzzz_0_xxxxzz_0, \
                             g_0_xxxzzz_0_xxxyyy_0, \
                             g_0_xxxzzz_0_xxxyyz_0, \
                             g_0_xxxzzz_0_xxxyzz_0, \
                             g_0_xxxzzz_0_xxxzzz_0, \
                             g_0_xxxzzz_0_xxyyyy_0, \
                             g_0_xxxzzz_0_xxyyyz_0, \
                             g_0_xxxzzz_0_xxyyzz_0, \
                             g_0_xxxzzz_0_xxyzzz_0, \
                             g_0_xxxzzz_0_xxzzzz_0, \
                             g_0_xxxzzz_0_xyyyyy_0, \
                             g_0_xxxzzz_0_xyyyyz_0, \
                             g_0_xxxzzz_0_xyyyzz_0, \
                             g_0_xxxzzz_0_xyyzzz_0, \
                             g_0_xxxzzz_0_xyzzzz_0, \
                             g_0_xxxzzz_0_xzzzzz_0, \
                             g_0_xxxzzz_0_yyyyyy_0, \
                             g_0_xxxzzz_0_yyyyyz_0, \
                             g_0_xxxzzz_0_yyyyzz_0, \
                             g_0_xxxzzz_0_yyyzzz_0, \
                             g_0_xxxzzz_0_yyzzzz_0, \
                             g_0_xxxzzz_0_yzzzzz_0, \
                             g_0_xxxzzz_0_zzzzzz_0, \
                             g_0_xxzzz_0_xxxxxz_0,  \
                             g_0_xxzzz_0_xxxxxz_1,  \
                             g_0_xxzzz_0_xxxxyz_0,  \
                             g_0_xxzzz_0_xxxxyz_1,  \
                             g_0_xxzzz_0_xxxxz_1,   \
                             g_0_xxzzz_0_xxxxzz_0,  \
                             g_0_xxzzz_0_xxxxzz_1,  \
                             g_0_xxzzz_0_xxxyyz_0,  \
                             g_0_xxzzz_0_xxxyyz_1,  \
                             g_0_xxzzz_0_xxxyz_1,   \
                             g_0_xxzzz_0_xxxyzz_0,  \
                             g_0_xxzzz_0_xxxyzz_1,  \
                             g_0_xxzzz_0_xxxzz_1,   \
                             g_0_xxzzz_0_xxxzzz_0,  \
                             g_0_xxzzz_0_xxxzzz_1,  \
                             g_0_xxzzz_0_xxyyyz_0,  \
                             g_0_xxzzz_0_xxyyyz_1,  \
                             g_0_xxzzz_0_xxyyz_1,   \
                             g_0_xxzzz_0_xxyyzz_0,  \
                             g_0_xxzzz_0_xxyyzz_1,  \
                             g_0_xxzzz_0_xxyzz_1,   \
                             g_0_xxzzz_0_xxyzzz_0,  \
                             g_0_xxzzz_0_xxyzzz_1,  \
                             g_0_xxzzz_0_xxzzz_1,   \
                             g_0_xxzzz_0_xxzzzz_0,  \
                             g_0_xxzzz_0_xxzzzz_1,  \
                             g_0_xxzzz_0_xyyyyz_0,  \
                             g_0_xxzzz_0_xyyyyz_1,  \
                             g_0_xxzzz_0_xyyyz_1,   \
                             g_0_xxzzz_0_xyyyzz_0,  \
                             g_0_xxzzz_0_xyyyzz_1,  \
                             g_0_xxzzz_0_xyyzz_1,   \
                             g_0_xxzzz_0_xyyzzz_0,  \
                             g_0_xxzzz_0_xyyzzz_1,  \
                             g_0_xxzzz_0_xyzzz_1,   \
                             g_0_xxzzz_0_xyzzzz_0,  \
                             g_0_xxzzz_0_xyzzzz_1,  \
                             g_0_xxzzz_0_xzzzz_1,   \
                             g_0_xxzzz_0_xzzzzz_0,  \
                             g_0_xxzzz_0_xzzzzz_1,  \
                             g_0_xxzzz_0_yyyyyy_0,  \
                             g_0_xxzzz_0_yyyyyy_1,  \
                             g_0_xxzzz_0_yyyyyz_0,  \
                             g_0_xxzzz_0_yyyyyz_1,  \
                             g_0_xxzzz_0_yyyyz_1,   \
                             g_0_xxzzz_0_yyyyzz_0,  \
                             g_0_xxzzz_0_yyyyzz_1,  \
                             g_0_xxzzz_0_yyyzz_1,   \
                             g_0_xxzzz_0_yyyzzz_0,  \
                             g_0_xxzzz_0_yyyzzz_1,  \
                             g_0_xxzzz_0_yyzzz_1,   \
                             g_0_xxzzz_0_yyzzzz_0,  \
                             g_0_xxzzz_0_yyzzzz_1,  \
                             g_0_xxzzz_0_yzzzz_1,   \
                             g_0_xxzzz_0_yzzzzz_0,  \
                             g_0_xxzzz_0_yzzzzz_1,  \
                             g_0_xxzzz_0_zzzzz_1,   \
                             g_0_xxzzz_0_zzzzzz_0,  \
                             g_0_xxzzz_0_zzzzzz_1,  \
                             g_0_xzzz_0_xxxxxz_0,   \
                             g_0_xzzz_0_xxxxxz_1,   \
                             g_0_xzzz_0_xxxxyz_0,   \
                             g_0_xzzz_0_xxxxyz_1,   \
                             g_0_xzzz_0_xxxxzz_0,   \
                             g_0_xzzz_0_xxxxzz_1,   \
                             g_0_xzzz_0_xxxyyz_0,   \
                             g_0_xzzz_0_xxxyyz_1,   \
                             g_0_xzzz_0_xxxyzz_0,   \
                             g_0_xzzz_0_xxxyzz_1,   \
                             g_0_xzzz_0_xxxzzz_0,   \
                             g_0_xzzz_0_xxxzzz_1,   \
                             g_0_xzzz_0_xxyyyz_0,   \
                             g_0_xzzz_0_xxyyyz_1,   \
                             g_0_xzzz_0_xxyyzz_0,   \
                             g_0_xzzz_0_xxyyzz_1,   \
                             g_0_xzzz_0_xxyzzz_0,   \
                             g_0_xzzz_0_xxyzzz_1,   \
                             g_0_xzzz_0_xxzzzz_0,   \
                             g_0_xzzz_0_xxzzzz_1,   \
                             g_0_xzzz_0_xyyyyz_0,   \
                             g_0_xzzz_0_xyyyyz_1,   \
                             g_0_xzzz_0_xyyyzz_0,   \
                             g_0_xzzz_0_xyyyzz_1,   \
                             g_0_xzzz_0_xyyzzz_0,   \
                             g_0_xzzz_0_xyyzzz_1,   \
                             g_0_xzzz_0_xyzzzz_0,   \
                             g_0_xzzz_0_xyzzzz_1,   \
                             g_0_xzzz_0_xzzzzz_0,   \
                             g_0_xzzz_0_xzzzzz_1,   \
                             g_0_xzzz_0_yyyyyy_0,   \
                             g_0_xzzz_0_yyyyyy_1,   \
                             g_0_xzzz_0_yyyyyz_0,   \
                             g_0_xzzz_0_yyyyyz_1,   \
                             g_0_xzzz_0_yyyyzz_0,   \
                             g_0_xzzz_0_yyyyzz_1,   \
                             g_0_xzzz_0_yyyzzz_0,   \
                             g_0_xzzz_0_yyyzzz_1,   \
                             g_0_xzzz_0_yyzzzz_0,   \
                             g_0_xzzz_0_yyzzzz_1,   \
                             g_0_xzzz_0_yzzzzz_0,   \
                             g_0_xzzz_0_yzzzzz_1,   \
                             g_0_xzzz_0_zzzzzz_0,   \
                             g_0_xzzz_0_zzzzzz_1,   \
                             wp_x,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzz_0_xxxxxx_0[i] = 2.0 * g_0_xxxz_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxxxx_0[i] * pb_z +
                                   g_0_xxxzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxxy_0[i] = 2.0 * g_0_xxxz_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxxxy_0[i] * pb_z +
                                   g_0_xxxzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxxz_0[i] = 2.0 * g_0_xzzz_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                   5.0 * g_0_xxzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxz_0[i] * pb_x + g_0_xxzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxyy_0[i] = 2.0 * g_0_xxxz_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxxyy_0[i] * pb_z +
                                   g_0_xxxzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxyz_0[i] = 2.0 * g_0_xzzz_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyz_0[i] * pb_x + g_0_xxzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxzz_0[i] = 2.0 * g_0_xzzz_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_xxzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxzz_0[i] * pb_x + g_0_xxzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxyyy_0[i] = 2.0 * g_0_xxxz_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxyyy_0[i] * pb_z +
                                   g_0_xxxzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxyyz_0[i] = 2.0 * g_0_xzzz_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyz_0[i] * pb_x + g_0_xxzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxyzz_0[i] = 2.0 * g_0_xzzz_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyzz_0[i] * pb_x + g_0_xxzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxzzz_0[i] = 2.0 * g_0_xzzz_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_xxzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxzzz_0[i] * pb_x + g_0_xxzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyyyy_0[i] = 2.0 * g_0_xxxz_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxyyyy_0[i] * pb_z +
                                   g_0_xxxzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxyyyz_0[i] = 2.0 * g_0_xzzz_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyz_0[i] * pb_x + g_0_xxzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyyzz_0[i] = 2.0 * g_0_xzzz_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyzz_0[i] * pb_x + g_0_xxzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyzzz_0[i] = 2.0 * g_0_xzzz_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyzzz_0[i] * pb_x + g_0_xxzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxzzzz_0[i] = 2.0 * g_0_xzzz_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_xxzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxzzzz_0[i] * pb_x + g_0_xxzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyyyy_0[i] = 2.0 * g_0_xxxz_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xyyyyy_0[i] * pb_z +
                                   g_0_xxxzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xyyyyz_0[i] = 2.0 * g_0_xzzz_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                   g_0_xxzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyz_0[i] * pb_x + g_0_xxzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyyzz_0[i] = 2.0 * g_0_xzzz_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                   g_0_xxzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyzz_0[i] * pb_x + g_0_xxzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyzzz_0[i] = 2.0 * g_0_xzzz_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                   g_0_xxzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyzzz_0[i] * pb_x + g_0_xxzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyzzzz_0[i] = 2.0 * g_0_xzzz_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                   g_0_xxzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyzzzz_0[i] * pb_x + g_0_xxzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xzzzzz_0[i] = 2.0 * g_0_xzzz_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                   g_0_xxzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xzzzzz_0[i] * pb_x + g_0_xxzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyyy_0[i] = 2.0 * g_0_xzzz_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyyyy_0[i] * pb_x +
                                   g_0_xxzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyyz_0[i] = 2.0 * g_0_xzzz_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyyyz_0[i] * pb_x +
                                   g_0_xxzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyzz_0[i] = 2.0 * g_0_xzzz_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyyzz_0[i] * pb_x +
                                   g_0_xxzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyzzz_0[i] = 2.0 * g_0_xzzz_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyzzz_0[i] * pb_x +
                                   g_0_xxzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyzzzz_0[i] = 2.0 * g_0_xzzz_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyzzzz_0[i] * pb_x +
                                   g_0_xxzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yzzzzz_0[i] = 2.0 * g_0_xzzz_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yzzzzz_0[i] * pb_x +
                                   g_0_xxzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_zzzzzz_0[i] = 2.0 * g_0_xzzz_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_zzzzzz_0[i] * pb_x +
                                   g_0_xxzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 280-308 components of targeted buffer : SISI

    auto g_0_xxyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 280);

    auto g_0_xxyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 281);

    auto g_0_xxyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 282);

    auto g_0_xxyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 283);

    auto g_0_xxyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 284);

    auto g_0_xxyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 285);

    auto g_0_xxyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 286);

    auto g_0_xxyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 287);

    auto g_0_xxyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 288);

    auto g_0_xxyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 289);

    auto g_0_xxyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 290);

    auto g_0_xxyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 291);

    auto g_0_xxyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 292);

    auto g_0_xxyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 293);

    auto g_0_xxyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 294);

    auto g_0_xxyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 295);

    auto g_0_xxyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 296);

    auto g_0_xxyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 297);

    auto g_0_xxyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 298);

    auto g_0_xxyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 299);

    auto g_0_xxyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 300);

    auto g_0_xxyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 301);

    auto g_0_xxyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 302);

    auto g_0_xxyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 303);

    auto g_0_xxyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 304);

    auto g_0_xxyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 305);

    auto g_0_xxyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 306);

    auto g_0_xxyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 307);

#pragma omp simd aligned(g_0_xxyy_0_xxxxxx_0,       \
                             g_0_xxyy_0_xxxxxx_1,   \
                             g_0_xxyy_0_xxxxxz_0,   \
                             g_0_xxyy_0_xxxxxz_1,   \
                             g_0_xxyy_0_xxxxzz_0,   \
                             g_0_xxyy_0_xxxxzz_1,   \
                             g_0_xxyy_0_xxxzzz_0,   \
                             g_0_xxyy_0_xxxzzz_1,   \
                             g_0_xxyy_0_xxzzzz_0,   \
                             g_0_xxyy_0_xxzzzz_1,   \
                             g_0_xxyy_0_xzzzzz_0,   \
                             g_0_xxyy_0_xzzzzz_1,   \
                             g_0_xxyyy_0_xxxxxx_0,  \
                             g_0_xxyyy_0_xxxxxx_1,  \
                             g_0_xxyyy_0_xxxxxz_0,  \
                             g_0_xxyyy_0_xxxxxz_1,  \
                             g_0_xxyyy_0_xxxxzz_0,  \
                             g_0_xxyyy_0_xxxxzz_1,  \
                             g_0_xxyyy_0_xxxzzz_0,  \
                             g_0_xxyyy_0_xxxzzz_1,  \
                             g_0_xxyyy_0_xxzzzz_0,  \
                             g_0_xxyyy_0_xxzzzz_1,  \
                             g_0_xxyyy_0_xzzzzz_0,  \
                             g_0_xxyyy_0_xzzzzz_1,  \
                             g_0_xxyyyy_0_xxxxxx_0, \
                             g_0_xxyyyy_0_xxxxxy_0, \
                             g_0_xxyyyy_0_xxxxxz_0, \
                             g_0_xxyyyy_0_xxxxyy_0, \
                             g_0_xxyyyy_0_xxxxyz_0, \
                             g_0_xxyyyy_0_xxxxzz_0, \
                             g_0_xxyyyy_0_xxxyyy_0, \
                             g_0_xxyyyy_0_xxxyyz_0, \
                             g_0_xxyyyy_0_xxxyzz_0, \
                             g_0_xxyyyy_0_xxxzzz_0, \
                             g_0_xxyyyy_0_xxyyyy_0, \
                             g_0_xxyyyy_0_xxyyyz_0, \
                             g_0_xxyyyy_0_xxyyzz_0, \
                             g_0_xxyyyy_0_xxyzzz_0, \
                             g_0_xxyyyy_0_xxzzzz_0, \
                             g_0_xxyyyy_0_xyyyyy_0, \
                             g_0_xxyyyy_0_xyyyyz_0, \
                             g_0_xxyyyy_0_xyyyzz_0, \
                             g_0_xxyyyy_0_xyyzzz_0, \
                             g_0_xxyyyy_0_xyzzzz_0, \
                             g_0_xxyyyy_0_xzzzzz_0, \
                             g_0_xxyyyy_0_yyyyyy_0, \
                             g_0_xxyyyy_0_yyyyyz_0, \
                             g_0_xxyyyy_0_yyyyzz_0, \
                             g_0_xxyyyy_0_yyyzzz_0, \
                             g_0_xxyyyy_0_yyzzzz_0, \
                             g_0_xxyyyy_0_yzzzzz_0, \
                             g_0_xxyyyy_0_zzzzzz_0, \
                             g_0_xyyyy_0_xxxxxy_0,  \
                             g_0_xyyyy_0_xxxxxy_1,  \
                             g_0_xyyyy_0_xxxxy_1,   \
                             g_0_xyyyy_0_xxxxyy_0,  \
                             g_0_xyyyy_0_xxxxyy_1,  \
                             g_0_xyyyy_0_xxxxyz_0,  \
                             g_0_xyyyy_0_xxxxyz_1,  \
                             g_0_xyyyy_0_xxxyy_1,   \
                             g_0_xyyyy_0_xxxyyy_0,  \
                             g_0_xyyyy_0_xxxyyy_1,  \
                             g_0_xyyyy_0_xxxyyz_0,  \
                             g_0_xyyyy_0_xxxyyz_1,  \
                             g_0_xyyyy_0_xxxyz_1,   \
                             g_0_xyyyy_0_xxxyzz_0,  \
                             g_0_xyyyy_0_xxxyzz_1,  \
                             g_0_xyyyy_0_xxyyy_1,   \
                             g_0_xyyyy_0_xxyyyy_0,  \
                             g_0_xyyyy_0_xxyyyy_1,  \
                             g_0_xyyyy_0_xxyyyz_0,  \
                             g_0_xyyyy_0_xxyyyz_1,  \
                             g_0_xyyyy_0_xxyyz_1,   \
                             g_0_xyyyy_0_xxyyzz_0,  \
                             g_0_xyyyy_0_xxyyzz_1,  \
                             g_0_xyyyy_0_xxyzz_1,   \
                             g_0_xyyyy_0_xxyzzz_0,  \
                             g_0_xyyyy_0_xxyzzz_1,  \
                             g_0_xyyyy_0_xyyyy_1,   \
                             g_0_xyyyy_0_xyyyyy_0,  \
                             g_0_xyyyy_0_xyyyyy_1,  \
                             g_0_xyyyy_0_xyyyyz_0,  \
                             g_0_xyyyy_0_xyyyyz_1,  \
                             g_0_xyyyy_0_xyyyz_1,   \
                             g_0_xyyyy_0_xyyyzz_0,  \
                             g_0_xyyyy_0_xyyyzz_1,  \
                             g_0_xyyyy_0_xyyzz_1,   \
                             g_0_xyyyy_0_xyyzzz_0,  \
                             g_0_xyyyy_0_xyyzzz_1,  \
                             g_0_xyyyy_0_xyzzz_1,   \
                             g_0_xyyyy_0_xyzzzz_0,  \
                             g_0_xyyyy_0_xyzzzz_1,  \
                             g_0_xyyyy_0_yyyyy_1,   \
                             g_0_xyyyy_0_yyyyyy_0,  \
                             g_0_xyyyy_0_yyyyyy_1,  \
                             g_0_xyyyy_0_yyyyyz_0,  \
                             g_0_xyyyy_0_yyyyyz_1,  \
                             g_0_xyyyy_0_yyyyz_1,   \
                             g_0_xyyyy_0_yyyyzz_0,  \
                             g_0_xyyyy_0_yyyyzz_1,  \
                             g_0_xyyyy_0_yyyzz_1,   \
                             g_0_xyyyy_0_yyyzzz_0,  \
                             g_0_xyyyy_0_yyyzzz_1,  \
                             g_0_xyyyy_0_yyzzz_1,   \
                             g_0_xyyyy_0_yyzzzz_0,  \
                             g_0_xyyyy_0_yyzzzz_1,  \
                             g_0_xyyyy_0_yzzzz_1,   \
                             g_0_xyyyy_0_yzzzzz_0,  \
                             g_0_xyyyy_0_yzzzzz_1,  \
                             g_0_xyyyy_0_zzzzzz_0,  \
                             g_0_xyyyy_0_zzzzzz_1,  \
                             g_0_yyyy_0_xxxxxy_0,   \
                             g_0_yyyy_0_xxxxxy_1,   \
                             g_0_yyyy_0_xxxxyy_0,   \
                             g_0_yyyy_0_xxxxyy_1,   \
                             g_0_yyyy_0_xxxxyz_0,   \
                             g_0_yyyy_0_xxxxyz_1,   \
                             g_0_yyyy_0_xxxyyy_0,   \
                             g_0_yyyy_0_xxxyyy_1,   \
                             g_0_yyyy_0_xxxyyz_0,   \
                             g_0_yyyy_0_xxxyyz_1,   \
                             g_0_yyyy_0_xxxyzz_0,   \
                             g_0_yyyy_0_xxxyzz_1,   \
                             g_0_yyyy_0_xxyyyy_0,   \
                             g_0_yyyy_0_xxyyyy_1,   \
                             g_0_yyyy_0_xxyyyz_0,   \
                             g_0_yyyy_0_xxyyyz_1,   \
                             g_0_yyyy_0_xxyyzz_0,   \
                             g_0_yyyy_0_xxyyzz_1,   \
                             g_0_yyyy_0_xxyzzz_0,   \
                             g_0_yyyy_0_xxyzzz_1,   \
                             g_0_yyyy_0_xyyyyy_0,   \
                             g_0_yyyy_0_xyyyyy_1,   \
                             g_0_yyyy_0_xyyyyz_0,   \
                             g_0_yyyy_0_xyyyyz_1,   \
                             g_0_yyyy_0_xyyyzz_0,   \
                             g_0_yyyy_0_xyyyzz_1,   \
                             g_0_yyyy_0_xyyzzz_0,   \
                             g_0_yyyy_0_xyyzzz_1,   \
                             g_0_yyyy_0_xyzzzz_0,   \
                             g_0_yyyy_0_xyzzzz_1,   \
                             g_0_yyyy_0_yyyyyy_0,   \
                             g_0_yyyy_0_yyyyyy_1,   \
                             g_0_yyyy_0_yyyyyz_0,   \
                             g_0_yyyy_0_yyyyyz_1,   \
                             g_0_yyyy_0_yyyyzz_0,   \
                             g_0_yyyy_0_yyyyzz_1,   \
                             g_0_yyyy_0_yyyzzz_0,   \
                             g_0_yyyy_0_yyyzzz_1,   \
                             g_0_yyyy_0_yyzzzz_0,   \
                             g_0_yyyy_0_yyzzzz_1,   \
                             g_0_yyyy_0_yzzzzz_0,   \
                             g_0_yyyy_0_yzzzzz_1,   \
                             g_0_yyyy_0_zzzzzz_0,   \
                             g_0_yyyy_0_zzzzzz_1,   \
                             wp_x,                  \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyy_0_xxxxxx_0[i] = 3.0 * g_0_xxyy_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxxxx_0[i] * pb_y +
                                   g_0_xxyyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxxxy_0[i] = g_0_yyyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxy_1[i] * fti_ab_0 + 5.0 * g_0_xyyyy_0_xxxxy_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xxxxxy_0[i] * pb_x + g_0_xyyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxxz_0[i] = 3.0 * g_0_xxyy_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxxxz_0[i] * pb_y +
                                   g_0_xxyyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxxyy_0[i] = g_0_yyyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxyy_1[i] * fti_ab_0 + 4.0 * g_0_xyyyy_0_xxxyy_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xxxxyy_0[i] * pb_x + g_0_xyyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxyz_0[i] = g_0_yyyy_0_xxxxyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyy_0_xxxyz_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xxxxyz_0[i] * pb_x + g_0_xyyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxzz_0[i] = 3.0 * g_0_xxyy_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxxzz_0[i] * pb_y +
                                   g_0_xxyyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxyyy_0[i] = g_0_yyyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_xyyyy_0_xxyyy_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xxxyyy_0[i] * pb_x + g_0_xyyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxyyz_0[i] = g_0_yyyy_0_xxxyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyy_0_xxyyz_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xxxyyz_0[i] * pb_x + g_0_xyyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxyzz_0[i] = g_0_yyyy_0_xxxyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyy_0_xxyzz_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xxxyzz_0[i] * pb_x + g_0_xyyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxzzz_0[i] = 3.0 * g_0_xxyy_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxzzz_0[i] * pb_y +
                                   g_0_xxyyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxyyyy_0[i] = g_0_yyyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyyyy_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xxyyyy_0[i] * pb_x + g_0_xyyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyyyz_0[i] = g_0_yyyy_0_xxyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyyyz_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xxyyyz_0[i] * pb_x + g_0_xyyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyyzz_0[i] = g_0_yyyy_0_xxyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyyzz_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xxyyzz_0[i] * pb_x + g_0_xyyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyzzz_0[i] = g_0_yyyy_0_xxyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyzzz_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xxyzzz_0[i] * pb_x + g_0_xyyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxzzzz_0[i] = 3.0 * g_0_xxyy_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxzzzz_0[i] * pb_y +
                                   g_0_xxyyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xyyyyy_0[i] = g_0_yyyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyy_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xyyyyy_0[i] * pb_x + g_0_xyyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyyyz_0[i] = g_0_yyyy_0_xyyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyyz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyz_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xyyyyz_0[i] * pb_x + g_0_xyyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyyzz_0[i] = g_0_yyyy_0_xyyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyzz_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xyyyzz_0[i] * pb_x + g_0_xyyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyzzz_0[i] = g_0_yyyy_0_xyyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyzzz_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xyyzzz_0[i] * pb_x + g_0_xyyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyzzzz_0[i] = g_0_yyyy_0_xyzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yzzzz_1[i] * fi_abcd_0 +
                                   g_0_xyyyy_0_xyzzzz_0[i] * pb_x + g_0_xyyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xzzzzz_0[i] = 3.0 * g_0_xxyy_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xzzzzz_0[i] * pb_y +
                                   g_0_xxyyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_yyyyyy_0[i] =
            g_0_yyyy_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyyy_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyy_0[i] * pb_x + g_0_xyyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyyyz_0[i] =
            g_0_yyyy_0_yyyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyyz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyz_0[i] * pb_x + g_0_xyyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyyzz_0[i] =
            g_0_yyyy_0_yyyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyzz_0[i] * pb_x + g_0_xyyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyzzz_0[i] =
            g_0_yyyy_0_yyyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyzzz_0[i] * pb_x + g_0_xyyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyzzzz_0[i] =
            g_0_yyyy_0_yyzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyzzzz_0[i] * pb_x + g_0_xyyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yzzzzz_0[i] =
            g_0_yyyy_0_yzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yzzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yzzzzz_0[i] * pb_x + g_0_xyyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_zzzzzz_0[i] =
            g_0_yyyy_0_zzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_zzzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_zzzzzz_0[i] * pb_x + g_0_xyyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 308-336 components of targeted buffer : SISI

    auto g_0_xxyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 308);

    auto g_0_xxyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 309);

    auto g_0_xxyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 310);

    auto g_0_xxyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 311);

    auto g_0_xxyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 312);

    auto g_0_xxyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 313);

    auto g_0_xxyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 314);

    auto g_0_xxyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 315);

    auto g_0_xxyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 316);

    auto g_0_xxyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 317);

    auto g_0_xxyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 318);

    auto g_0_xxyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 319);

    auto g_0_xxyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 320);

    auto g_0_xxyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 321);

    auto g_0_xxyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 322);

    auto g_0_xxyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 323);

    auto g_0_xxyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 324);

    auto g_0_xxyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 325);

    auto g_0_xxyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 326);

    auto g_0_xxyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 327);

    auto g_0_xxyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 328);

    auto g_0_xxyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 329);

    auto g_0_xxyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 330);

    auto g_0_xxyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 331);

    auto g_0_xxyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 332);

    auto g_0_xxyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 333);

    auto g_0_xxyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 334);

    auto g_0_xxyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 335);

#pragma omp simd aligned(g_0_xxyyy_0_xxxxx_1,       \
                             g_0_xxyyy_0_xxxxxx_0,  \
                             g_0_xxyyy_0_xxxxxx_1,  \
                             g_0_xxyyy_0_xxxxxy_0,  \
                             g_0_xxyyy_0_xxxxxy_1,  \
                             g_0_xxyyy_0_xxxxxz_0,  \
                             g_0_xxyyy_0_xxxxxz_1,  \
                             g_0_xxyyy_0_xxxxy_1,   \
                             g_0_xxyyy_0_xxxxyy_0,  \
                             g_0_xxyyy_0_xxxxyy_1,  \
                             g_0_xxyyy_0_xxxxyz_0,  \
                             g_0_xxyyy_0_xxxxyz_1,  \
                             g_0_xxyyy_0_xxxxz_1,   \
                             g_0_xxyyy_0_xxxxzz_0,  \
                             g_0_xxyyy_0_xxxxzz_1,  \
                             g_0_xxyyy_0_xxxyy_1,   \
                             g_0_xxyyy_0_xxxyyy_0,  \
                             g_0_xxyyy_0_xxxyyy_1,  \
                             g_0_xxyyy_0_xxxyyz_0,  \
                             g_0_xxyyy_0_xxxyyz_1,  \
                             g_0_xxyyy_0_xxxyz_1,   \
                             g_0_xxyyy_0_xxxyzz_0,  \
                             g_0_xxyyy_0_xxxyzz_1,  \
                             g_0_xxyyy_0_xxxzz_1,   \
                             g_0_xxyyy_0_xxxzzz_0,  \
                             g_0_xxyyy_0_xxxzzz_1,  \
                             g_0_xxyyy_0_xxyyy_1,   \
                             g_0_xxyyy_0_xxyyyy_0,  \
                             g_0_xxyyy_0_xxyyyy_1,  \
                             g_0_xxyyy_0_xxyyyz_0,  \
                             g_0_xxyyy_0_xxyyyz_1,  \
                             g_0_xxyyy_0_xxyyz_1,   \
                             g_0_xxyyy_0_xxyyzz_0,  \
                             g_0_xxyyy_0_xxyyzz_1,  \
                             g_0_xxyyy_0_xxyzz_1,   \
                             g_0_xxyyy_0_xxyzzz_0,  \
                             g_0_xxyyy_0_xxyzzz_1,  \
                             g_0_xxyyy_0_xxzzz_1,   \
                             g_0_xxyyy_0_xxzzzz_0,  \
                             g_0_xxyyy_0_xxzzzz_1,  \
                             g_0_xxyyy_0_xyyyy_1,   \
                             g_0_xxyyy_0_xyyyyy_0,  \
                             g_0_xxyyy_0_xyyyyy_1,  \
                             g_0_xxyyy_0_xyyyyz_0,  \
                             g_0_xxyyy_0_xyyyyz_1,  \
                             g_0_xxyyy_0_xyyyz_1,   \
                             g_0_xxyyy_0_xyyyzz_0,  \
                             g_0_xxyyy_0_xyyyzz_1,  \
                             g_0_xxyyy_0_xyyzz_1,   \
                             g_0_xxyyy_0_xyyzzz_0,  \
                             g_0_xxyyy_0_xyyzzz_1,  \
                             g_0_xxyyy_0_xyzzz_1,   \
                             g_0_xxyyy_0_xyzzzz_0,  \
                             g_0_xxyyy_0_xyzzzz_1,  \
                             g_0_xxyyy_0_xzzzz_1,   \
                             g_0_xxyyy_0_xzzzzz_0,  \
                             g_0_xxyyy_0_xzzzzz_1,  \
                             g_0_xxyyy_0_yyyyy_1,   \
                             g_0_xxyyy_0_yyyyyy_0,  \
                             g_0_xxyyy_0_yyyyyy_1,  \
                             g_0_xxyyy_0_yyyyyz_0,  \
                             g_0_xxyyy_0_yyyyyz_1,  \
                             g_0_xxyyy_0_yyyyz_1,   \
                             g_0_xxyyy_0_yyyyzz_0,  \
                             g_0_xxyyy_0_yyyyzz_1,  \
                             g_0_xxyyy_0_yyyzz_1,   \
                             g_0_xxyyy_0_yyyzzz_0,  \
                             g_0_xxyyy_0_yyyzzz_1,  \
                             g_0_xxyyy_0_yyzzz_1,   \
                             g_0_xxyyy_0_yyzzzz_0,  \
                             g_0_xxyyy_0_yyzzzz_1,  \
                             g_0_xxyyy_0_yzzzz_1,   \
                             g_0_xxyyy_0_yzzzzz_0,  \
                             g_0_xxyyy_0_yzzzzz_1,  \
                             g_0_xxyyy_0_zzzzz_1,   \
                             g_0_xxyyy_0_zzzzzz_0,  \
                             g_0_xxyyy_0_zzzzzz_1,  \
                             g_0_xxyyyz_0_xxxxxx_0, \
                             g_0_xxyyyz_0_xxxxxy_0, \
                             g_0_xxyyyz_0_xxxxxz_0, \
                             g_0_xxyyyz_0_xxxxyy_0, \
                             g_0_xxyyyz_0_xxxxyz_0, \
                             g_0_xxyyyz_0_xxxxzz_0, \
                             g_0_xxyyyz_0_xxxyyy_0, \
                             g_0_xxyyyz_0_xxxyyz_0, \
                             g_0_xxyyyz_0_xxxyzz_0, \
                             g_0_xxyyyz_0_xxxzzz_0, \
                             g_0_xxyyyz_0_xxyyyy_0, \
                             g_0_xxyyyz_0_xxyyyz_0, \
                             g_0_xxyyyz_0_xxyyzz_0, \
                             g_0_xxyyyz_0_xxyzzz_0, \
                             g_0_xxyyyz_0_xxzzzz_0, \
                             g_0_xxyyyz_0_xyyyyy_0, \
                             g_0_xxyyyz_0_xyyyyz_0, \
                             g_0_xxyyyz_0_xyyyzz_0, \
                             g_0_xxyyyz_0_xyyzzz_0, \
                             g_0_xxyyyz_0_xyzzzz_0, \
                             g_0_xxyyyz_0_xzzzzz_0, \
                             g_0_xxyyyz_0_yyyyyy_0, \
                             g_0_xxyyyz_0_yyyyyz_0, \
                             g_0_xxyyyz_0_yyyyzz_0, \
                             g_0_xxyyyz_0_yyyzzz_0, \
                             g_0_xxyyyz_0_yyzzzz_0, \
                             g_0_xxyyyz_0_yzzzzz_0, \
                             g_0_xxyyyz_0_zzzzzz_0, \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyz_0_xxxxxx_0[i] = g_0_xxyyy_0_xxxxxx_0[i] * pb_z + g_0_xxyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxy_0[i] = g_0_xxyyy_0_xxxxxy_0[i] * pb_z + g_0_xxyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxz_0[i] = g_0_xxyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxz_0[i] * pb_z + g_0_xxyyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxyy_0[i] = g_0_xxyyy_0_xxxxyy_0[i] * pb_z + g_0_xxyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxyz_0[i] = g_0_xxyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyz_0[i] * pb_z + g_0_xxyyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxzz_0[i] = 2.0 * g_0_xxyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxzz_0[i] * pb_z + g_0_xxyyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyyy_0[i] = g_0_xxyyy_0_xxxyyy_0[i] * pb_z + g_0_xxyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyyz_0[i] = g_0_xxyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyz_0[i] * pb_z + g_0_xxyyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyzz_0[i] = 2.0 * g_0_xxyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyzz_0[i] * pb_z + g_0_xxyyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxzzz_0[i] = 3.0 * g_0_xxyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxzzz_0[i] * pb_z + g_0_xxyyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyyy_0[i] = g_0_xxyyy_0_xxyyyy_0[i] * pb_z + g_0_xxyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyyz_0[i] = g_0_xxyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyz_0[i] * pb_z + g_0_xxyyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyzz_0[i] = 2.0 * g_0_xxyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyzz_0[i] * pb_z + g_0_xxyyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyzzz_0[i] = 3.0 * g_0_xxyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyzzz_0[i] * pb_z + g_0_xxyyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxzzzz_0[i] = 4.0 * g_0_xxyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxzzzz_0[i] * pb_z + g_0_xxyyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyyy_0[i] = g_0_xxyyy_0_xyyyyy_0[i] * pb_z + g_0_xxyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyyz_0[i] = g_0_xxyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyz_0[i] * pb_z + g_0_xxyyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyzz_0[i] = 2.0 * g_0_xxyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyzz_0[i] * pb_z + g_0_xxyyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyzzz_0[i] = 3.0 * g_0_xxyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyzzz_0[i] * pb_z + g_0_xxyyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyzzzz_0[i] = 4.0 * g_0_xxyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyzzzz_0[i] * pb_z + g_0_xxyyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xzzzzz_0[i] = 5.0 * g_0_xxyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xzzzzz_0[i] * pb_z + g_0_xxyyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyyy_0[i] = g_0_xxyyy_0_yyyyyy_0[i] * pb_z + g_0_xxyyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyyz_0[i] = g_0_xxyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyyyz_0[i] * pb_z + g_0_xxyyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyzz_0[i] = 2.0 * g_0_xxyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyyzz_0[i] * pb_z + g_0_xxyyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyzzz_0[i] = 3.0 * g_0_xxyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyzzz_0[i] * pb_z + g_0_xxyyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyzzzz_0[i] = 4.0 * g_0_xxyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyzzzz_0[i] * pb_z + g_0_xxyyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yzzzzz_0[i] = 5.0 * g_0_xxyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yzzzzz_0[i] * pb_z + g_0_xxyyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_zzzzzz_0[i] = 6.0 * g_0_xxyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_zzzzzz_0[i] * pb_z + g_0_xxyyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 336-364 components of targeted buffer : SISI

    auto g_0_xxyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 336);

    auto g_0_xxyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 337);

    auto g_0_xxyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 338);

    auto g_0_xxyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 339);

    auto g_0_xxyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 340);

    auto g_0_xxyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 341);

    auto g_0_xxyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 342);

    auto g_0_xxyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 343);

    auto g_0_xxyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 344);

    auto g_0_xxyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 345);

    auto g_0_xxyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 346);

    auto g_0_xxyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 347);

    auto g_0_xxyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 348);

    auto g_0_xxyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 349);

    auto g_0_xxyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 350);

    auto g_0_xxyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 351);

    auto g_0_xxyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 352);

    auto g_0_xxyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 353);

    auto g_0_xxyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 354);

    auto g_0_xxyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 355);

    auto g_0_xxyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 356);

    auto g_0_xxyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 357);

    auto g_0_xxyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 358);

    auto g_0_xxyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 359);

    auto g_0_xxyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 360);

    auto g_0_xxyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 361);

    auto g_0_xxyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 362);

    auto g_0_xxyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 363);

#pragma omp simd aligned(g_0_xxyy_0_xxxxxy_0,       \
                             g_0_xxyy_0_xxxxxy_1,   \
                             g_0_xxyy_0_xxxxyy_0,   \
                             g_0_xxyy_0_xxxxyy_1,   \
                             g_0_xxyy_0_xxxyyy_0,   \
                             g_0_xxyy_0_xxxyyy_1,   \
                             g_0_xxyy_0_xxyyyy_0,   \
                             g_0_xxyy_0_xxyyyy_1,   \
                             g_0_xxyy_0_xyyyyy_0,   \
                             g_0_xxyy_0_xyyyyy_1,   \
                             g_0_xxyyz_0_xxxxxy_0,  \
                             g_0_xxyyz_0_xxxxxy_1,  \
                             g_0_xxyyz_0_xxxxyy_0,  \
                             g_0_xxyyz_0_xxxxyy_1,  \
                             g_0_xxyyz_0_xxxyyy_0,  \
                             g_0_xxyyz_0_xxxyyy_1,  \
                             g_0_xxyyz_0_xxyyyy_0,  \
                             g_0_xxyyz_0_xxyyyy_1,  \
                             g_0_xxyyz_0_xyyyyy_0,  \
                             g_0_xxyyz_0_xyyyyy_1,  \
                             g_0_xxyyzz_0_xxxxxx_0, \
                             g_0_xxyyzz_0_xxxxxy_0, \
                             g_0_xxyyzz_0_xxxxxz_0, \
                             g_0_xxyyzz_0_xxxxyy_0, \
                             g_0_xxyyzz_0_xxxxyz_0, \
                             g_0_xxyyzz_0_xxxxzz_0, \
                             g_0_xxyyzz_0_xxxyyy_0, \
                             g_0_xxyyzz_0_xxxyyz_0, \
                             g_0_xxyyzz_0_xxxyzz_0, \
                             g_0_xxyyzz_0_xxxzzz_0, \
                             g_0_xxyyzz_0_xxyyyy_0, \
                             g_0_xxyyzz_0_xxyyyz_0, \
                             g_0_xxyyzz_0_xxyyzz_0, \
                             g_0_xxyyzz_0_xxyzzz_0, \
                             g_0_xxyyzz_0_xxzzzz_0, \
                             g_0_xxyyzz_0_xyyyyy_0, \
                             g_0_xxyyzz_0_xyyyyz_0, \
                             g_0_xxyyzz_0_xyyyzz_0, \
                             g_0_xxyyzz_0_xyyzzz_0, \
                             g_0_xxyyzz_0_xyzzzz_0, \
                             g_0_xxyyzz_0_xzzzzz_0, \
                             g_0_xxyyzz_0_yyyyyy_0, \
                             g_0_xxyyzz_0_yyyyyz_0, \
                             g_0_xxyyzz_0_yyyyzz_0, \
                             g_0_xxyyzz_0_yyyzzz_0, \
                             g_0_xxyyzz_0_yyzzzz_0, \
                             g_0_xxyyzz_0_yzzzzz_0, \
                             g_0_xxyyzz_0_zzzzzz_0, \
                             g_0_xxyzz_0_xxxxxx_0,  \
                             g_0_xxyzz_0_xxxxxx_1,  \
                             g_0_xxyzz_0_xxxxxz_0,  \
                             g_0_xxyzz_0_xxxxxz_1,  \
                             g_0_xxyzz_0_xxxxzz_0,  \
                             g_0_xxyzz_0_xxxxzz_1,  \
                             g_0_xxyzz_0_xxxzzz_0,  \
                             g_0_xxyzz_0_xxxzzz_1,  \
                             g_0_xxyzz_0_xxzzzz_0,  \
                             g_0_xxyzz_0_xxzzzz_1,  \
                             g_0_xxyzz_0_xzzzzz_0,  \
                             g_0_xxyzz_0_xzzzzz_1,  \
                             g_0_xxzz_0_xxxxxx_0,   \
                             g_0_xxzz_0_xxxxxx_1,   \
                             g_0_xxzz_0_xxxxxz_0,   \
                             g_0_xxzz_0_xxxxxz_1,   \
                             g_0_xxzz_0_xxxxzz_0,   \
                             g_0_xxzz_0_xxxxzz_1,   \
                             g_0_xxzz_0_xxxzzz_0,   \
                             g_0_xxzz_0_xxxzzz_1,   \
                             g_0_xxzz_0_xxzzzz_0,   \
                             g_0_xxzz_0_xxzzzz_1,   \
                             g_0_xxzz_0_xzzzzz_0,   \
                             g_0_xxzz_0_xzzzzz_1,   \
                             g_0_xyyzz_0_xxxxyz_0,  \
                             g_0_xyyzz_0_xxxxyz_1,  \
                             g_0_xyyzz_0_xxxyyz_0,  \
                             g_0_xyyzz_0_xxxyyz_1,  \
                             g_0_xyyzz_0_xxxyz_1,   \
                             g_0_xyyzz_0_xxxyzz_0,  \
                             g_0_xyyzz_0_xxxyzz_1,  \
                             g_0_xyyzz_0_xxyyyz_0,  \
                             g_0_xyyzz_0_xxyyyz_1,  \
                             g_0_xyyzz_0_xxyyz_1,   \
                             g_0_xyyzz_0_xxyyzz_0,  \
                             g_0_xyyzz_0_xxyyzz_1,  \
                             g_0_xyyzz_0_xxyzz_1,   \
                             g_0_xyyzz_0_xxyzzz_0,  \
                             g_0_xyyzz_0_xxyzzz_1,  \
                             g_0_xyyzz_0_xyyyyz_0,  \
                             g_0_xyyzz_0_xyyyyz_1,  \
                             g_0_xyyzz_0_xyyyz_1,   \
                             g_0_xyyzz_0_xyyyzz_0,  \
                             g_0_xyyzz_0_xyyyzz_1,  \
                             g_0_xyyzz_0_xyyzz_1,   \
                             g_0_xyyzz_0_xyyzzz_0,  \
                             g_0_xyyzz_0_xyyzzz_1,  \
                             g_0_xyyzz_0_xyzzz_1,   \
                             g_0_xyyzz_0_xyzzzz_0,  \
                             g_0_xyyzz_0_xyzzzz_1,  \
                             g_0_xyyzz_0_yyyyyy_0,  \
                             g_0_xyyzz_0_yyyyyy_1,  \
                             g_0_xyyzz_0_yyyyyz_0,  \
                             g_0_xyyzz_0_yyyyyz_1,  \
                             g_0_xyyzz_0_yyyyz_1,   \
                             g_0_xyyzz_0_yyyyzz_0,  \
                             g_0_xyyzz_0_yyyyzz_1,  \
                             g_0_xyyzz_0_yyyzz_1,   \
                             g_0_xyyzz_0_yyyzzz_0,  \
                             g_0_xyyzz_0_yyyzzz_1,  \
                             g_0_xyyzz_0_yyzzz_1,   \
                             g_0_xyyzz_0_yyzzzz_0,  \
                             g_0_xyyzz_0_yyzzzz_1,  \
                             g_0_xyyzz_0_yzzzz_1,   \
                             g_0_xyyzz_0_yzzzzz_0,  \
                             g_0_xyyzz_0_yzzzzz_1,  \
                             g_0_xyyzz_0_zzzzzz_0,  \
                             g_0_xyyzz_0_zzzzzz_1,  \
                             g_0_yyzz_0_xxxxyz_0,   \
                             g_0_yyzz_0_xxxxyz_1,   \
                             g_0_yyzz_0_xxxyyz_0,   \
                             g_0_yyzz_0_xxxyyz_1,   \
                             g_0_yyzz_0_xxxyzz_0,   \
                             g_0_yyzz_0_xxxyzz_1,   \
                             g_0_yyzz_0_xxyyyz_0,   \
                             g_0_yyzz_0_xxyyyz_1,   \
                             g_0_yyzz_0_xxyyzz_0,   \
                             g_0_yyzz_0_xxyyzz_1,   \
                             g_0_yyzz_0_xxyzzz_0,   \
                             g_0_yyzz_0_xxyzzz_1,   \
                             g_0_yyzz_0_xyyyyz_0,   \
                             g_0_yyzz_0_xyyyyz_1,   \
                             g_0_yyzz_0_xyyyzz_0,   \
                             g_0_yyzz_0_xyyyzz_1,   \
                             g_0_yyzz_0_xyyzzz_0,   \
                             g_0_yyzz_0_xyyzzz_1,   \
                             g_0_yyzz_0_xyzzzz_0,   \
                             g_0_yyzz_0_xyzzzz_1,   \
                             g_0_yyzz_0_yyyyyy_0,   \
                             g_0_yyzz_0_yyyyyy_1,   \
                             g_0_yyzz_0_yyyyyz_0,   \
                             g_0_yyzz_0_yyyyyz_1,   \
                             g_0_yyzz_0_yyyyzz_0,   \
                             g_0_yyzz_0_yyyyzz_1,   \
                             g_0_yyzz_0_yyyzzz_0,   \
                             g_0_yyzz_0_yyyzzz_1,   \
                             g_0_yyzz_0_yyzzzz_0,   \
                             g_0_yyzz_0_yyzzzz_1,   \
                             g_0_yyzz_0_yzzzzz_0,   \
                             g_0_yyzz_0_yzzzzz_1,   \
                             g_0_yyzz_0_zzzzzz_0,   \
                             g_0_yyzz_0_zzzzzz_1,   \
                             wp_x,                  \
                             wp_y,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzz_0_xxxxxx_0[i] =
            g_0_xxzz_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxxx_0[i] * pb_y + g_0_xxyzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxxxy_0[i] =
            g_0_xxyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxxxy_0[i] * pb_z + g_0_xxyyz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxxxz_0[i] =
            g_0_xxzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxxz_0[i] * pb_y + g_0_xxyzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxxyy_0[i] =
            g_0_xxyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxxyy_0[i] * pb_z + g_0_xxyyz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxxyz_0[i] = g_0_yyzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxxyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyzz_0_xxxyz_1[i] * fi_abcd_0 +
                                   g_0_xyyzz_0_xxxxyz_0[i] * pb_x + g_0_xyyzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxxzz_0[i] =
            g_0_xxzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxzz_0[i] * pb_y + g_0_xxyzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxyyy_0[i] =
            g_0_xxyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxyyy_0[i] * pb_z + g_0_xxyyz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxyyz_0[i] = g_0_yyzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzz_0_xxyyz_1[i] * fi_abcd_0 +
                                   g_0_xyyzz_0_xxxyyz_0[i] * pb_x + g_0_xyyzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxyzz_0[i] = g_0_yyzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzz_0_xxyzz_1[i] * fi_abcd_0 +
                                   g_0_xyyzz_0_xxxyzz_0[i] * pb_x + g_0_xyyzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxzzz_0[i] =
            g_0_xxzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxzzz_0[i] * pb_y + g_0_xxyzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxyyyy_0[i] =
            g_0_xxyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxyyyy_0[i] * pb_z + g_0_xxyyz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxyyyz_0[i] = g_0_yyzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzz_0_xyyyz_1[i] * fi_abcd_0 +
                                   g_0_xyyzz_0_xxyyyz_0[i] * pb_x + g_0_xyyzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxyyzz_0[i] = g_0_yyzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzz_0_xyyzz_1[i] * fi_abcd_0 +
                                   g_0_xyyzz_0_xxyyzz_0[i] * pb_x + g_0_xyyzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxyzzz_0[i] = g_0_yyzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzz_0_xyzzz_1[i] * fi_abcd_0 +
                                   g_0_xyyzz_0_xxyzzz_0[i] * pb_x + g_0_xyyzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxzzzz_0[i] =
            g_0_xxzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxzzzz_0[i] * pb_y + g_0_xxyzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xyyyyy_0[i] =
            g_0_xxyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xyyyyy_0[i] * pb_z + g_0_xxyyz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xyyyyz_0[i] = g_0_yyzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyyyz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyz_1[i] * fi_abcd_0 +
                                   g_0_xyyzz_0_xyyyyz_0[i] * pb_x + g_0_xyyzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyyyzz_0[i] = g_0_yyzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyyzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyzz_1[i] * fi_abcd_0 +
                                   g_0_xyyzz_0_xyyyzz_0[i] * pb_x + g_0_xyyzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyyzzz_0[i] = g_0_yyzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyzzz_1[i] * fi_abcd_0 +
                                   g_0_xyyzz_0_xyyzzz_0[i] * pb_x + g_0_xyyzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyzzzz_0[i] = g_0_yyzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yzzzz_1[i] * fi_abcd_0 +
                                   g_0_xyyzz_0_xyzzzz_0[i] * pb_x + g_0_xyyzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xzzzzz_0[i] =
            g_0_xxzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xzzzzz_0[i] * pb_y + g_0_xxyzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_yyyyyy_0[i] =
            g_0_yyzz_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyyy_0[i] * pb_x + g_0_xyyzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyyyz_0[i] =
            g_0_yyzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyyz_0[i] * pb_x + g_0_xyyzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyyzz_0[i] =
            g_0_yyzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyzz_0[i] * pb_x + g_0_xyyzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyzzz_0[i] =
            g_0_yyzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyzzz_0[i] * pb_x + g_0_xyyzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyzzzz_0[i] =
            g_0_yyzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyzzzz_0[i] * pb_x + g_0_xyyzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yzzzzz_0[i] =
            g_0_yyzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yzzzzz_0[i] * pb_x + g_0_xyyzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_zzzzzz_0[i] =
            g_0_yyzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_zzzzzz_0[i] * pb_x + g_0_xyyzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 364-392 components of targeted buffer : SISI

    auto g_0_xxyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 364);

    auto g_0_xxyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 365);

    auto g_0_xxyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 366);

    auto g_0_xxyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 367);

    auto g_0_xxyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 368);

    auto g_0_xxyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 369);

    auto g_0_xxyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 370);

    auto g_0_xxyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 371);

    auto g_0_xxyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 372);

    auto g_0_xxyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 373);

    auto g_0_xxyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 374);

    auto g_0_xxyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 375);

    auto g_0_xxyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 376);

    auto g_0_xxyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 377);

    auto g_0_xxyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 378);

    auto g_0_xxyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 379);

    auto g_0_xxyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 380);

    auto g_0_xxyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 381);

    auto g_0_xxyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 382);

    auto g_0_xxyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 383);

    auto g_0_xxyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 384);

    auto g_0_xxyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 385);

    auto g_0_xxyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 386);

    auto g_0_xxyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 387);

    auto g_0_xxyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 388);

    auto g_0_xxyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 389);

    auto g_0_xxyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 390);

    auto g_0_xxyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 391);

#pragma omp simd aligned(g_0_xxyzzz_0_xxxxxx_0,     \
                             g_0_xxyzzz_0_xxxxxy_0, \
                             g_0_xxyzzz_0_xxxxxz_0, \
                             g_0_xxyzzz_0_xxxxyy_0, \
                             g_0_xxyzzz_0_xxxxyz_0, \
                             g_0_xxyzzz_0_xxxxzz_0, \
                             g_0_xxyzzz_0_xxxyyy_0, \
                             g_0_xxyzzz_0_xxxyyz_0, \
                             g_0_xxyzzz_0_xxxyzz_0, \
                             g_0_xxyzzz_0_xxxzzz_0, \
                             g_0_xxyzzz_0_xxyyyy_0, \
                             g_0_xxyzzz_0_xxyyyz_0, \
                             g_0_xxyzzz_0_xxyyzz_0, \
                             g_0_xxyzzz_0_xxyzzz_0, \
                             g_0_xxyzzz_0_xxzzzz_0, \
                             g_0_xxyzzz_0_xyyyyy_0, \
                             g_0_xxyzzz_0_xyyyyz_0, \
                             g_0_xxyzzz_0_xyyyzz_0, \
                             g_0_xxyzzz_0_xyyzzz_0, \
                             g_0_xxyzzz_0_xyzzzz_0, \
                             g_0_xxyzzz_0_xzzzzz_0, \
                             g_0_xxyzzz_0_yyyyyy_0, \
                             g_0_xxyzzz_0_yyyyyz_0, \
                             g_0_xxyzzz_0_yyyyzz_0, \
                             g_0_xxyzzz_0_yyyzzz_0, \
                             g_0_xxyzzz_0_yyzzzz_0, \
                             g_0_xxyzzz_0_yzzzzz_0, \
                             g_0_xxyzzz_0_zzzzzz_0, \
                             g_0_xxzzz_0_xxxxx_1,   \
                             g_0_xxzzz_0_xxxxxx_0,  \
                             g_0_xxzzz_0_xxxxxx_1,  \
                             g_0_xxzzz_0_xxxxxy_0,  \
                             g_0_xxzzz_0_xxxxxy_1,  \
                             g_0_xxzzz_0_xxxxxz_0,  \
                             g_0_xxzzz_0_xxxxxz_1,  \
                             g_0_xxzzz_0_xxxxy_1,   \
                             g_0_xxzzz_0_xxxxyy_0,  \
                             g_0_xxzzz_0_xxxxyy_1,  \
                             g_0_xxzzz_0_xxxxyz_0,  \
                             g_0_xxzzz_0_xxxxyz_1,  \
                             g_0_xxzzz_0_xxxxz_1,   \
                             g_0_xxzzz_0_xxxxzz_0,  \
                             g_0_xxzzz_0_xxxxzz_1,  \
                             g_0_xxzzz_0_xxxyy_1,   \
                             g_0_xxzzz_0_xxxyyy_0,  \
                             g_0_xxzzz_0_xxxyyy_1,  \
                             g_0_xxzzz_0_xxxyyz_0,  \
                             g_0_xxzzz_0_xxxyyz_1,  \
                             g_0_xxzzz_0_xxxyz_1,   \
                             g_0_xxzzz_0_xxxyzz_0,  \
                             g_0_xxzzz_0_xxxyzz_1,  \
                             g_0_xxzzz_0_xxxzz_1,   \
                             g_0_xxzzz_0_xxxzzz_0,  \
                             g_0_xxzzz_0_xxxzzz_1,  \
                             g_0_xxzzz_0_xxyyy_1,   \
                             g_0_xxzzz_0_xxyyyy_0,  \
                             g_0_xxzzz_0_xxyyyy_1,  \
                             g_0_xxzzz_0_xxyyyz_0,  \
                             g_0_xxzzz_0_xxyyyz_1,  \
                             g_0_xxzzz_0_xxyyz_1,   \
                             g_0_xxzzz_0_xxyyzz_0,  \
                             g_0_xxzzz_0_xxyyzz_1,  \
                             g_0_xxzzz_0_xxyzz_1,   \
                             g_0_xxzzz_0_xxyzzz_0,  \
                             g_0_xxzzz_0_xxyzzz_1,  \
                             g_0_xxzzz_0_xxzzz_1,   \
                             g_0_xxzzz_0_xxzzzz_0,  \
                             g_0_xxzzz_0_xxzzzz_1,  \
                             g_0_xxzzz_0_xyyyy_1,   \
                             g_0_xxzzz_0_xyyyyy_0,  \
                             g_0_xxzzz_0_xyyyyy_1,  \
                             g_0_xxzzz_0_xyyyyz_0,  \
                             g_0_xxzzz_0_xyyyyz_1,  \
                             g_0_xxzzz_0_xyyyz_1,   \
                             g_0_xxzzz_0_xyyyzz_0,  \
                             g_0_xxzzz_0_xyyyzz_1,  \
                             g_0_xxzzz_0_xyyzz_1,   \
                             g_0_xxzzz_0_xyyzzz_0,  \
                             g_0_xxzzz_0_xyyzzz_1,  \
                             g_0_xxzzz_0_xyzzz_1,   \
                             g_0_xxzzz_0_xyzzzz_0,  \
                             g_0_xxzzz_0_xyzzzz_1,  \
                             g_0_xxzzz_0_xzzzz_1,   \
                             g_0_xxzzz_0_xzzzzz_0,  \
                             g_0_xxzzz_0_xzzzzz_1,  \
                             g_0_xxzzz_0_yyyyy_1,   \
                             g_0_xxzzz_0_yyyyyy_0,  \
                             g_0_xxzzz_0_yyyyyy_1,  \
                             g_0_xxzzz_0_yyyyyz_0,  \
                             g_0_xxzzz_0_yyyyyz_1,  \
                             g_0_xxzzz_0_yyyyz_1,   \
                             g_0_xxzzz_0_yyyyzz_0,  \
                             g_0_xxzzz_0_yyyyzz_1,  \
                             g_0_xxzzz_0_yyyzz_1,   \
                             g_0_xxzzz_0_yyyzzz_0,  \
                             g_0_xxzzz_0_yyyzzz_1,  \
                             g_0_xxzzz_0_yyzzz_1,   \
                             g_0_xxzzz_0_yyzzzz_0,  \
                             g_0_xxzzz_0_yyzzzz_1,  \
                             g_0_xxzzz_0_yzzzz_1,   \
                             g_0_xxzzz_0_yzzzzz_0,  \
                             g_0_xxzzz_0_yzzzzz_1,  \
                             g_0_xxzzz_0_zzzzz_1,   \
                             g_0_xxzzz_0_zzzzzz_0,  \
                             g_0_xxzzz_0_zzzzzz_1,  \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzz_0_xxxxxx_0[i] = g_0_xxzzz_0_xxxxxx_0[i] * pb_y + g_0_xxzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxy_0[i] = g_0_xxzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxy_0[i] * pb_y + g_0_xxzzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxz_0[i] = g_0_xxzzz_0_xxxxxz_0[i] * pb_y + g_0_xxzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxyy_0[i] = 2.0 * g_0_xxzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyy_0[i] * pb_y + g_0_xxzzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxyz_0[i] = g_0_xxzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyz_0[i] * pb_y + g_0_xxzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxzz_0[i] = g_0_xxzzz_0_xxxxzz_0[i] * pb_y + g_0_xxzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyyy_0[i] = 3.0 * g_0_xxzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyy_0[i] * pb_y + g_0_xxzzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyyz_0[i] = 2.0 * g_0_xxzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyz_0[i] * pb_y + g_0_xxzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyzz_0[i] = g_0_xxzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyzz_0[i] * pb_y + g_0_xxzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxzzz_0[i] = g_0_xxzzz_0_xxxzzz_0[i] * pb_y + g_0_xxzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyyy_0[i] = 4.0 * g_0_xxzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyy_0[i] * pb_y + g_0_xxzzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyyz_0[i] = 3.0 * g_0_xxzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyz_0[i] * pb_y + g_0_xxzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyzz_0[i] = 2.0 * g_0_xxzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyzz_0[i] * pb_y + g_0_xxzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyzzz_0[i] = g_0_xxzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyzzz_0[i] * pb_y + g_0_xxzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxzzzz_0[i] = g_0_xxzzz_0_xxzzzz_0[i] * pb_y + g_0_xxzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyyy_0[i] = 5.0 * g_0_xxzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyy_0[i] * pb_y + g_0_xxzzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyyz_0[i] = 4.0 * g_0_xxzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyz_0[i] * pb_y + g_0_xxzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyzz_0[i] = 3.0 * g_0_xxzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyzz_0[i] * pb_y + g_0_xxzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyzzz_0[i] = 2.0 * g_0_xxzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyzzz_0[i] * pb_y + g_0_xxzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyzzzz_0[i] = g_0_xxzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyzzzz_0[i] * pb_y + g_0_xxzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xzzzzz_0[i] = g_0_xxzzz_0_xzzzzz_0[i] * pb_y + g_0_xxzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyyy_0[i] = 6.0 * g_0_xxzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyyy_0[i] * pb_y + g_0_xxzzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyyz_0[i] = 5.0 * g_0_xxzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyyz_0[i] * pb_y + g_0_xxzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyzz_0[i] = 4.0 * g_0_xxzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyzz_0[i] * pb_y + g_0_xxzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyzzz_0[i] = 3.0 * g_0_xxzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyzzz_0[i] * pb_y + g_0_xxzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyzzzz_0[i] = 2.0 * g_0_xxzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyzzzz_0[i] * pb_y + g_0_xxzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yzzzzz_0[i] = g_0_xxzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yzzzzz_0[i] * pb_y + g_0_xxzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_zzzzzz_0[i] = g_0_xxzzz_0_zzzzzz_0[i] * pb_y + g_0_xxzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 392-420 components of targeted buffer : SISI

    auto g_0_xxzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 392);

    auto g_0_xxzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 393);

    auto g_0_xxzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 394);

    auto g_0_xxzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 395);

    auto g_0_xxzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 396);

    auto g_0_xxzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 397);

    auto g_0_xxzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 398);

    auto g_0_xxzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 399);

    auto g_0_xxzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 400);

    auto g_0_xxzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 401);

    auto g_0_xxzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 402);

    auto g_0_xxzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 403);

    auto g_0_xxzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 404);

    auto g_0_xxzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 405);

    auto g_0_xxzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 406);

    auto g_0_xxzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 407);

    auto g_0_xxzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 408);

    auto g_0_xxzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 409);

    auto g_0_xxzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 410);

    auto g_0_xxzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 411);

    auto g_0_xxzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 412);

    auto g_0_xxzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 413);

    auto g_0_xxzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 414);

    auto g_0_xxzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 415);

    auto g_0_xxzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 416);

    auto g_0_xxzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 417);

    auto g_0_xxzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 418);

    auto g_0_xxzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 419);

#pragma omp simd aligned(g_0_xxzz_0_xxxxxx_0,       \
                             g_0_xxzz_0_xxxxxx_1,   \
                             g_0_xxzz_0_xxxxxy_0,   \
                             g_0_xxzz_0_xxxxxy_1,   \
                             g_0_xxzz_0_xxxxyy_0,   \
                             g_0_xxzz_0_xxxxyy_1,   \
                             g_0_xxzz_0_xxxyyy_0,   \
                             g_0_xxzz_0_xxxyyy_1,   \
                             g_0_xxzz_0_xxyyyy_0,   \
                             g_0_xxzz_0_xxyyyy_1,   \
                             g_0_xxzz_0_xyyyyy_0,   \
                             g_0_xxzz_0_xyyyyy_1,   \
                             g_0_xxzzz_0_xxxxxx_0,  \
                             g_0_xxzzz_0_xxxxxx_1,  \
                             g_0_xxzzz_0_xxxxxy_0,  \
                             g_0_xxzzz_0_xxxxxy_1,  \
                             g_0_xxzzz_0_xxxxyy_0,  \
                             g_0_xxzzz_0_xxxxyy_1,  \
                             g_0_xxzzz_0_xxxyyy_0,  \
                             g_0_xxzzz_0_xxxyyy_1,  \
                             g_0_xxzzz_0_xxyyyy_0,  \
                             g_0_xxzzz_0_xxyyyy_1,  \
                             g_0_xxzzz_0_xyyyyy_0,  \
                             g_0_xxzzz_0_xyyyyy_1,  \
                             g_0_xxzzzz_0_xxxxxx_0, \
                             g_0_xxzzzz_0_xxxxxy_0, \
                             g_0_xxzzzz_0_xxxxxz_0, \
                             g_0_xxzzzz_0_xxxxyy_0, \
                             g_0_xxzzzz_0_xxxxyz_0, \
                             g_0_xxzzzz_0_xxxxzz_0, \
                             g_0_xxzzzz_0_xxxyyy_0, \
                             g_0_xxzzzz_0_xxxyyz_0, \
                             g_0_xxzzzz_0_xxxyzz_0, \
                             g_0_xxzzzz_0_xxxzzz_0, \
                             g_0_xxzzzz_0_xxyyyy_0, \
                             g_0_xxzzzz_0_xxyyyz_0, \
                             g_0_xxzzzz_0_xxyyzz_0, \
                             g_0_xxzzzz_0_xxyzzz_0, \
                             g_0_xxzzzz_0_xxzzzz_0, \
                             g_0_xxzzzz_0_xyyyyy_0, \
                             g_0_xxzzzz_0_xyyyyz_0, \
                             g_0_xxzzzz_0_xyyyzz_0, \
                             g_0_xxzzzz_0_xyyzzz_0, \
                             g_0_xxzzzz_0_xyzzzz_0, \
                             g_0_xxzzzz_0_xzzzzz_0, \
                             g_0_xxzzzz_0_yyyyyy_0, \
                             g_0_xxzzzz_0_yyyyyz_0, \
                             g_0_xxzzzz_0_yyyyzz_0, \
                             g_0_xxzzzz_0_yyyzzz_0, \
                             g_0_xxzzzz_0_yyzzzz_0, \
                             g_0_xxzzzz_0_yzzzzz_0, \
                             g_0_xxzzzz_0_zzzzzz_0, \
                             g_0_xzzzz_0_xxxxxz_0,  \
                             g_0_xzzzz_0_xxxxxz_1,  \
                             g_0_xzzzz_0_xxxxyz_0,  \
                             g_0_xzzzz_0_xxxxyz_1,  \
                             g_0_xzzzz_0_xxxxz_1,   \
                             g_0_xzzzz_0_xxxxzz_0,  \
                             g_0_xzzzz_0_xxxxzz_1,  \
                             g_0_xzzzz_0_xxxyyz_0,  \
                             g_0_xzzzz_0_xxxyyz_1,  \
                             g_0_xzzzz_0_xxxyz_1,   \
                             g_0_xzzzz_0_xxxyzz_0,  \
                             g_0_xzzzz_0_xxxyzz_1,  \
                             g_0_xzzzz_0_xxxzz_1,   \
                             g_0_xzzzz_0_xxxzzz_0,  \
                             g_0_xzzzz_0_xxxzzz_1,  \
                             g_0_xzzzz_0_xxyyyz_0,  \
                             g_0_xzzzz_0_xxyyyz_1,  \
                             g_0_xzzzz_0_xxyyz_1,   \
                             g_0_xzzzz_0_xxyyzz_0,  \
                             g_0_xzzzz_0_xxyyzz_1,  \
                             g_0_xzzzz_0_xxyzz_1,   \
                             g_0_xzzzz_0_xxyzzz_0,  \
                             g_0_xzzzz_0_xxyzzz_1,  \
                             g_0_xzzzz_0_xxzzz_1,   \
                             g_0_xzzzz_0_xxzzzz_0,  \
                             g_0_xzzzz_0_xxzzzz_1,  \
                             g_0_xzzzz_0_xyyyyz_0,  \
                             g_0_xzzzz_0_xyyyyz_1,  \
                             g_0_xzzzz_0_xyyyz_1,   \
                             g_0_xzzzz_0_xyyyzz_0,  \
                             g_0_xzzzz_0_xyyyzz_1,  \
                             g_0_xzzzz_0_xyyzz_1,   \
                             g_0_xzzzz_0_xyyzzz_0,  \
                             g_0_xzzzz_0_xyyzzz_1,  \
                             g_0_xzzzz_0_xyzzz_1,   \
                             g_0_xzzzz_0_xyzzzz_0,  \
                             g_0_xzzzz_0_xyzzzz_1,  \
                             g_0_xzzzz_0_xzzzz_1,   \
                             g_0_xzzzz_0_xzzzzz_0,  \
                             g_0_xzzzz_0_xzzzzz_1,  \
                             g_0_xzzzz_0_yyyyyy_0,  \
                             g_0_xzzzz_0_yyyyyy_1,  \
                             g_0_xzzzz_0_yyyyyz_0,  \
                             g_0_xzzzz_0_yyyyyz_1,  \
                             g_0_xzzzz_0_yyyyz_1,   \
                             g_0_xzzzz_0_yyyyzz_0,  \
                             g_0_xzzzz_0_yyyyzz_1,  \
                             g_0_xzzzz_0_yyyzz_1,   \
                             g_0_xzzzz_0_yyyzzz_0,  \
                             g_0_xzzzz_0_yyyzzz_1,  \
                             g_0_xzzzz_0_yyzzz_1,   \
                             g_0_xzzzz_0_yyzzzz_0,  \
                             g_0_xzzzz_0_yyzzzz_1,  \
                             g_0_xzzzz_0_yzzzz_1,   \
                             g_0_xzzzz_0_yzzzzz_0,  \
                             g_0_xzzzz_0_yzzzzz_1,  \
                             g_0_xzzzz_0_zzzzz_1,   \
                             g_0_xzzzz_0_zzzzzz_0,  \
                             g_0_xzzzz_0_zzzzzz_1,  \
                             g_0_zzzz_0_xxxxxz_0,   \
                             g_0_zzzz_0_xxxxxz_1,   \
                             g_0_zzzz_0_xxxxyz_0,   \
                             g_0_zzzz_0_xxxxyz_1,   \
                             g_0_zzzz_0_xxxxzz_0,   \
                             g_0_zzzz_0_xxxxzz_1,   \
                             g_0_zzzz_0_xxxyyz_0,   \
                             g_0_zzzz_0_xxxyyz_1,   \
                             g_0_zzzz_0_xxxyzz_0,   \
                             g_0_zzzz_0_xxxyzz_1,   \
                             g_0_zzzz_0_xxxzzz_0,   \
                             g_0_zzzz_0_xxxzzz_1,   \
                             g_0_zzzz_0_xxyyyz_0,   \
                             g_0_zzzz_0_xxyyyz_1,   \
                             g_0_zzzz_0_xxyyzz_0,   \
                             g_0_zzzz_0_xxyyzz_1,   \
                             g_0_zzzz_0_xxyzzz_0,   \
                             g_0_zzzz_0_xxyzzz_1,   \
                             g_0_zzzz_0_xxzzzz_0,   \
                             g_0_zzzz_0_xxzzzz_1,   \
                             g_0_zzzz_0_xyyyyz_0,   \
                             g_0_zzzz_0_xyyyyz_1,   \
                             g_0_zzzz_0_xyyyzz_0,   \
                             g_0_zzzz_0_xyyyzz_1,   \
                             g_0_zzzz_0_xyyzzz_0,   \
                             g_0_zzzz_0_xyyzzz_1,   \
                             g_0_zzzz_0_xyzzzz_0,   \
                             g_0_zzzz_0_xyzzzz_1,   \
                             g_0_zzzz_0_xzzzzz_0,   \
                             g_0_zzzz_0_xzzzzz_1,   \
                             g_0_zzzz_0_yyyyyy_0,   \
                             g_0_zzzz_0_yyyyyy_1,   \
                             g_0_zzzz_0_yyyyyz_0,   \
                             g_0_zzzz_0_yyyyyz_1,   \
                             g_0_zzzz_0_yyyyzz_0,   \
                             g_0_zzzz_0_yyyyzz_1,   \
                             g_0_zzzz_0_yyyzzz_0,   \
                             g_0_zzzz_0_yyyzzz_1,   \
                             g_0_zzzz_0_yyzzzz_0,   \
                             g_0_zzzz_0_yyzzzz_1,   \
                             g_0_zzzz_0_yzzzzz_0,   \
                             g_0_zzzz_0_yzzzzz_1,   \
                             g_0_zzzz_0_zzzzzz_0,   \
                             g_0_zzzz_0_zzzzzz_1,   \
                             wp_x,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzz_0_xxxxxx_0[i] = 3.0 * g_0_xxzz_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxxxx_0[i] * pb_z +
                                   g_0_xxzzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxxy_0[i] = 3.0 * g_0_xxzz_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxxxy_0[i] * pb_z +
                                   g_0_xxzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxxz_0[i] = g_0_zzzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxz_1[i] * fti_ab_0 + 5.0 * g_0_xzzzz_0_xxxxz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xxxxxz_0[i] * pb_x + g_0_xzzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxyy_0[i] = 3.0 * g_0_xxzz_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxxyy_0[i] * pb_z +
                                   g_0_xxzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxyz_0[i] = g_0_zzzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxyz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzz_0_xxxyz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xxxxyz_0[i] * pb_x + g_0_xzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxzz_0[i] = g_0_zzzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzz_0_xxxzz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xxxxzz_0[i] * pb_x + g_0_xzzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxyyy_0[i] = 3.0 * g_0_xxzz_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxyyy_0[i] * pb_z +
                                   g_0_xxzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxyyz_0[i] = g_0_zzzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyyz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzz_0_xxyyz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xxxyyz_0[i] * pb_x + g_0_xzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxyzz_0[i] = g_0_zzzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzz_0_xxyzz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xxxyzz_0[i] * pb_x + g_0_xzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxzzz_0[i] = g_0_zzzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzz_0_xxzzz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xxxzzz_0[i] * pb_x + g_0_xzzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyyyy_0[i] = 3.0 * g_0_xxzz_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxyyyy_0[i] * pb_z +
                                   g_0_xxzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxyyyz_0[i] = g_0_zzzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xyyyz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xxyyyz_0[i] * pb_x + g_0_xzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyyzz_0[i] = g_0_zzzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xyyzz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xxyyzz_0[i] * pb_x + g_0_xzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyzzz_0[i] = g_0_zzzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xyzzz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xxyzzz_0[i] * pb_x + g_0_xzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxzzzz_0[i] = g_0_zzzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xxzzzz_0[i] * pb_x + g_0_xzzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyyyy_0[i] = 3.0 * g_0_xxzz_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xyyyyy_0[i] * pb_z +
                                   g_0_xxzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xyyyyz_0[i] = g_0_zzzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyyz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xyyyyz_0[i] * pb_x + g_0_xzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyyzz_0[i] = g_0_zzzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyzz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xyyyzz_0[i] * pb_x + g_0_xzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyzzz_0[i] = g_0_zzzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyzzz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xyyzzz_0[i] * pb_x + g_0_xzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyzzzz_0[i] = g_0_zzzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xyzzzz_0[i] * pb_x + g_0_xzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xzzzzz_0[i] = g_0_zzzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_zzzzz_1[i] * fi_abcd_0 +
                                   g_0_xzzzz_0_xzzzzz_0[i] * pb_x + g_0_xzzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyyy_0[i] =
            g_0_zzzz_0_yyyyyy_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyyy_0[i] * pb_x + g_0_xzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyyz_0[i] =
            g_0_zzzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyyz_0[i] * pb_x + g_0_xzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyzz_0[i] =
            g_0_zzzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyzz_0[i] * pb_x + g_0_xzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyzzz_0[i] =
            g_0_zzzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyzzz_0[i] * pb_x + g_0_xzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyzzzz_0[i] =
            g_0_zzzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyzzzz_0[i] * pb_x + g_0_xzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yzzzzz_0[i] =
            g_0_zzzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yzzzzz_0[i] * pb_x + g_0_xzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_zzzzzz_0[i] =
            g_0_zzzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_zzzzzz_0[i] * pb_x + g_0_xzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 420-448 components of targeted buffer : SISI

    auto g_0_xyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 420);

    auto g_0_xyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 421);

    auto g_0_xyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 422);

    auto g_0_xyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 423);

    auto g_0_xyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 424);

    auto g_0_xyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 425);

    auto g_0_xyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 426);

    auto g_0_xyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 427);

    auto g_0_xyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 428);

    auto g_0_xyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 429);

    auto g_0_xyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 430);

    auto g_0_xyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 431);

    auto g_0_xyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 432);

    auto g_0_xyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 433);

    auto g_0_xyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 434);

    auto g_0_xyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 435);

    auto g_0_xyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 436);

    auto g_0_xyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 437);

    auto g_0_xyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 438);

    auto g_0_xyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 439);

    auto g_0_xyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 440);

    auto g_0_xyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 441);

    auto g_0_xyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 442);

    auto g_0_xyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 443);

    auto g_0_xyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 444);

    auto g_0_xyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 445);

    auto g_0_xyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 446);

    auto g_0_xyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 447);

#pragma omp simd aligned(g_0_xyyyyy_0_xxxxxx_0,     \
                             g_0_xyyyyy_0_xxxxxy_0, \
                             g_0_xyyyyy_0_xxxxxz_0, \
                             g_0_xyyyyy_0_xxxxyy_0, \
                             g_0_xyyyyy_0_xxxxyz_0, \
                             g_0_xyyyyy_0_xxxxzz_0, \
                             g_0_xyyyyy_0_xxxyyy_0, \
                             g_0_xyyyyy_0_xxxyyz_0, \
                             g_0_xyyyyy_0_xxxyzz_0, \
                             g_0_xyyyyy_0_xxxzzz_0, \
                             g_0_xyyyyy_0_xxyyyy_0, \
                             g_0_xyyyyy_0_xxyyyz_0, \
                             g_0_xyyyyy_0_xxyyzz_0, \
                             g_0_xyyyyy_0_xxyzzz_0, \
                             g_0_xyyyyy_0_xxzzzz_0, \
                             g_0_xyyyyy_0_xyyyyy_0, \
                             g_0_xyyyyy_0_xyyyyz_0, \
                             g_0_xyyyyy_0_xyyyzz_0, \
                             g_0_xyyyyy_0_xyyzzz_0, \
                             g_0_xyyyyy_0_xyzzzz_0, \
                             g_0_xyyyyy_0_xzzzzz_0, \
                             g_0_xyyyyy_0_yyyyyy_0, \
                             g_0_xyyyyy_0_yyyyyz_0, \
                             g_0_xyyyyy_0_yyyyzz_0, \
                             g_0_xyyyyy_0_yyyzzz_0, \
                             g_0_xyyyyy_0_yyzzzz_0, \
                             g_0_xyyyyy_0_yzzzzz_0, \
                             g_0_xyyyyy_0_zzzzzz_0, \
                             g_0_yyyyy_0_xxxxx_1,   \
                             g_0_yyyyy_0_xxxxxx_0,  \
                             g_0_yyyyy_0_xxxxxx_1,  \
                             g_0_yyyyy_0_xxxxxy_0,  \
                             g_0_yyyyy_0_xxxxxy_1,  \
                             g_0_yyyyy_0_xxxxxz_0,  \
                             g_0_yyyyy_0_xxxxxz_1,  \
                             g_0_yyyyy_0_xxxxy_1,   \
                             g_0_yyyyy_0_xxxxyy_0,  \
                             g_0_yyyyy_0_xxxxyy_1,  \
                             g_0_yyyyy_0_xxxxyz_0,  \
                             g_0_yyyyy_0_xxxxyz_1,  \
                             g_0_yyyyy_0_xxxxz_1,   \
                             g_0_yyyyy_0_xxxxzz_0,  \
                             g_0_yyyyy_0_xxxxzz_1,  \
                             g_0_yyyyy_0_xxxyy_1,   \
                             g_0_yyyyy_0_xxxyyy_0,  \
                             g_0_yyyyy_0_xxxyyy_1,  \
                             g_0_yyyyy_0_xxxyyz_0,  \
                             g_0_yyyyy_0_xxxyyz_1,  \
                             g_0_yyyyy_0_xxxyz_1,   \
                             g_0_yyyyy_0_xxxyzz_0,  \
                             g_0_yyyyy_0_xxxyzz_1,  \
                             g_0_yyyyy_0_xxxzz_1,   \
                             g_0_yyyyy_0_xxxzzz_0,  \
                             g_0_yyyyy_0_xxxzzz_1,  \
                             g_0_yyyyy_0_xxyyy_1,   \
                             g_0_yyyyy_0_xxyyyy_0,  \
                             g_0_yyyyy_0_xxyyyy_1,  \
                             g_0_yyyyy_0_xxyyyz_0,  \
                             g_0_yyyyy_0_xxyyyz_1,  \
                             g_0_yyyyy_0_xxyyz_1,   \
                             g_0_yyyyy_0_xxyyzz_0,  \
                             g_0_yyyyy_0_xxyyzz_1,  \
                             g_0_yyyyy_0_xxyzz_1,   \
                             g_0_yyyyy_0_xxyzzz_0,  \
                             g_0_yyyyy_0_xxyzzz_1,  \
                             g_0_yyyyy_0_xxzzz_1,   \
                             g_0_yyyyy_0_xxzzzz_0,  \
                             g_0_yyyyy_0_xxzzzz_1,  \
                             g_0_yyyyy_0_xyyyy_1,   \
                             g_0_yyyyy_0_xyyyyy_0,  \
                             g_0_yyyyy_0_xyyyyy_1,  \
                             g_0_yyyyy_0_xyyyyz_0,  \
                             g_0_yyyyy_0_xyyyyz_1,  \
                             g_0_yyyyy_0_xyyyz_1,   \
                             g_0_yyyyy_0_xyyyzz_0,  \
                             g_0_yyyyy_0_xyyyzz_1,  \
                             g_0_yyyyy_0_xyyzz_1,   \
                             g_0_yyyyy_0_xyyzzz_0,  \
                             g_0_yyyyy_0_xyyzzz_1,  \
                             g_0_yyyyy_0_xyzzz_1,   \
                             g_0_yyyyy_0_xyzzzz_0,  \
                             g_0_yyyyy_0_xyzzzz_1,  \
                             g_0_yyyyy_0_xzzzz_1,   \
                             g_0_yyyyy_0_xzzzzz_0,  \
                             g_0_yyyyy_0_xzzzzz_1,  \
                             g_0_yyyyy_0_yyyyy_1,   \
                             g_0_yyyyy_0_yyyyyy_0,  \
                             g_0_yyyyy_0_yyyyyy_1,  \
                             g_0_yyyyy_0_yyyyyz_0,  \
                             g_0_yyyyy_0_yyyyyz_1,  \
                             g_0_yyyyy_0_yyyyz_1,   \
                             g_0_yyyyy_0_yyyyzz_0,  \
                             g_0_yyyyy_0_yyyyzz_1,  \
                             g_0_yyyyy_0_yyyzz_1,   \
                             g_0_yyyyy_0_yyyzzz_0,  \
                             g_0_yyyyy_0_yyyzzz_1,  \
                             g_0_yyyyy_0_yyzzz_1,   \
                             g_0_yyyyy_0_yyzzzz_0,  \
                             g_0_yyyyy_0_yyzzzz_1,  \
                             g_0_yyyyy_0_yzzzz_1,   \
                             g_0_yyyyy_0_yzzzzz_0,  \
                             g_0_yyyyy_0_yzzzzz_1,  \
                             g_0_yyyyy_0_zzzzz_1,   \
                             g_0_yyyyy_0_zzzzzz_0,  \
                             g_0_yyyyy_0_zzzzzz_1,  \
                             wp_x,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyy_0_xxxxxx_0[i] = 6.0 * g_0_yyyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxx_0[i] * pb_x + g_0_yyyyy_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxy_0[i] = 5.0 * g_0_yyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxy_0[i] * pb_x + g_0_yyyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxz_0[i] = 5.0 * g_0_yyyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxz_0[i] * pb_x + g_0_yyyyy_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxyy_0[i] = 4.0 * g_0_yyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyy_0[i] * pb_x + g_0_yyyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxyz_0[i] = 4.0 * g_0_yyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyz_0[i] * pb_x + g_0_yyyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxzz_0[i] = 4.0 * g_0_yyyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxzz_0[i] * pb_x + g_0_yyyyy_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyyy_0[i] = 3.0 * g_0_yyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyy_0[i] * pb_x + g_0_yyyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyyz_0[i] = 3.0 * g_0_yyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyz_0[i] * pb_x + g_0_yyyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyzz_0[i] = 3.0 * g_0_yyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyzz_0[i] * pb_x + g_0_yyyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxzzz_0[i] = 3.0 * g_0_yyyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxzzz_0[i] * pb_x + g_0_yyyyy_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyyy_0[i] = 2.0 * g_0_yyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyy_0[i] * pb_x + g_0_yyyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyyz_0[i] = 2.0 * g_0_yyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyz_0[i] * pb_x + g_0_yyyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyzz_0[i] = 2.0 * g_0_yyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyzz_0[i] * pb_x + g_0_yyyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyzzz_0[i] = 2.0 * g_0_yyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyzzz_0[i] * pb_x + g_0_yyyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxzzzz_0[i] = 2.0 * g_0_yyyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxzzzz_0[i] * pb_x + g_0_yyyyy_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyyy_0[i] = g_0_yyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyy_0[i] * pb_x + g_0_yyyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyyz_0[i] = g_0_yyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyz_0[i] * pb_x + g_0_yyyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyzz_0[i] = g_0_yyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyzz_0[i] * pb_x + g_0_yyyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyzzz_0[i] = g_0_yyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyzzz_0[i] * pb_x + g_0_yyyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyzzzz_0[i] = g_0_yyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzzzz_0[i] * pb_x + g_0_yyyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xzzzzz_0[i] = g_0_yyyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xzzzzz_0[i] * pb_x + g_0_yyyyy_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyyy_0[i] = g_0_yyyyy_0_yyyyyy_0[i] * pb_x + g_0_yyyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyyz_0[i] = g_0_yyyyy_0_yyyyyz_0[i] * pb_x + g_0_yyyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyzz_0[i] = g_0_yyyyy_0_yyyyzz_0[i] * pb_x + g_0_yyyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyzzz_0[i] = g_0_yyyyy_0_yyyzzz_0[i] * pb_x + g_0_yyyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyzzzz_0[i] = g_0_yyyyy_0_yyzzzz_0[i] * pb_x + g_0_yyyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yzzzzz_0[i] = g_0_yyyyy_0_yzzzzz_0[i] * pb_x + g_0_yyyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_zzzzzz_0[i] = g_0_yyyyy_0_zzzzzz_0[i] * pb_x + g_0_yyyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 448-476 components of targeted buffer : SISI

    auto g_0_xyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 448);

    auto g_0_xyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 449);

    auto g_0_xyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 450);

    auto g_0_xyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 451);

    auto g_0_xyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 452);

    auto g_0_xyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 453);

    auto g_0_xyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 454);

    auto g_0_xyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 455);

    auto g_0_xyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 456);

    auto g_0_xyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 457);

    auto g_0_xyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 458);

    auto g_0_xyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 459);

    auto g_0_xyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 460);

    auto g_0_xyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 461);

    auto g_0_xyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 462);

    auto g_0_xyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 463);

    auto g_0_xyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 464);

    auto g_0_xyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 465);

    auto g_0_xyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 466);

    auto g_0_xyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 467);

    auto g_0_xyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 468);

    auto g_0_xyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 469);

    auto g_0_xyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 470);

    auto g_0_xyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 471);

    auto g_0_xyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 472);

    auto g_0_xyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 473);

    auto g_0_xyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 474);

    auto g_0_xyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 475);

#pragma omp simd aligned(g_0_xyyyy_0_xxxxxx_0,      \
                             g_0_xyyyy_0_xxxxxx_1,  \
                             g_0_xyyyy_0_xxxxxy_0,  \
                             g_0_xyyyy_0_xxxxxy_1,  \
                             g_0_xyyyy_0_xxxxyy_0,  \
                             g_0_xyyyy_0_xxxxyy_1,  \
                             g_0_xyyyy_0_xxxyyy_0,  \
                             g_0_xyyyy_0_xxxyyy_1,  \
                             g_0_xyyyy_0_xxyyyy_0,  \
                             g_0_xyyyy_0_xxyyyy_1,  \
                             g_0_xyyyy_0_xyyyyy_0,  \
                             g_0_xyyyy_0_xyyyyy_1,  \
                             g_0_xyyyyz_0_xxxxxx_0, \
                             g_0_xyyyyz_0_xxxxxy_0, \
                             g_0_xyyyyz_0_xxxxxz_0, \
                             g_0_xyyyyz_0_xxxxyy_0, \
                             g_0_xyyyyz_0_xxxxyz_0, \
                             g_0_xyyyyz_0_xxxxzz_0, \
                             g_0_xyyyyz_0_xxxyyy_0, \
                             g_0_xyyyyz_0_xxxyyz_0, \
                             g_0_xyyyyz_0_xxxyzz_0, \
                             g_0_xyyyyz_0_xxxzzz_0, \
                             g_0_xyyyyz_0_xxyyyy_0, \
                             g_0_xyyyyz_0_xxyyyz_0, \
                             g_0_xyyyyz_0_xxyyzz_0, \
                             g_0_xyyyyz_0_xxyzzz_0, \
                             g_0_xyyyyz_0_xxzzzz_0, \
                             g_0_xyyyyz_0_xyyyyy_0, \
                             g_0_xyyyyz_0_xyyyyz_0, \
                             g_0_xyyyyz_0_xyyyzz_0, \
                             g_0_xyyyyz_0_xyyzzz_0, \
                             g_0_xyyyyz_0_xyzzzz_0, \
                             g_0_xyyyyz_0_xzzzzz_0, \
                             g_0_xyyyyz_0_yyyyyy_0, \
                             g_0_xyyyyz_0_yyyyyz_0, \
                             g_0_xyyyyz_0_yyyyzz_0, \
                             g_0_xyyyyz_0_yyyzzz_0, \
                             g_0_xyyyyz_0_yyzzzz_0, \
                             g_0_xyyyyz_0_yzzzzz_0, \
                             g_0_xyyyyz_0_zzzzzz_0, \
                             g_0_yyyyz_0_xxxxxz_0,  \
                             g_0_yyyyz_0_xxxxxz_1,  \
                             g_0_yyyyz_0_xxxxyz_0,  \
                             g_0_yyyyz_0_xxxxyz_1,  \
                             g_0_yyyyz_0_xxxxz_1,   \
                             g_0_yyyyz_0_xxxxzz_0,  \
                             g_0_yyyyz_0_xxxxzz_1,  \
                             g_0_yyyyz_0_xxxyyz_0,  \
                             g_0_yyyyz_0_xxxyyz_1,  \
                             g_0_yyyyz_0_xxxyz_1,   \
                             g_0_yyyyz_0_xxxyzz_0,  \
                             g_0_yyyyz_0_xxxyzz_1,  \
                             g_0_yyyyz_0_xxxzz_1,   \
                             g_0_yyyyz_0_xxxzzz_0,  \
                             g_0_yyyyz_0_xxxzzz_1,  \
                             g_0_yyyyz_0_xxyyyz_0,  \
                             g_0_yyyyz_0_xxyyyz_1,  \
                             g_0_yyyyz_0_xxyyz_1,   \
                             g_0_yyyyz_0_xxyyzz_0,  \
                             g_0_yyyyz_0_xxyyzz_1,  \
                             g_0_yyyyz_0_xxyzz_1,   \
                             g_0_yyyyz_0_xxyzzz_0,  \
                             g_0_yyyyz_0_xxyzzz_1,  \
                             g_0_yyyyz_0_xxzzz_1,   \
                             g_0_yyyyz_0_xxzzzz_0,  \
                             g_0_yyyyz_0_xxzzzz_1,  \
                             g_0_yyyyz_0_xyyyyz_0,  \
                             g_0_yyyyz_0_xyyyyz_1,  \
                             g_0_yyyyz_0_xyyyz_1,   \
                             g_0_yyyyz_0_xyyyzz_0,  \
                             g_0_yyyyz_0_xyyyzz_1,  \
                             g_0_yyyyz_0_xyyzz_1,   \
                             g_0_yyyyz_0_xyyzzz_0,  \
                             g_0_yyyyz_0_xyyzzz_1,  \
                             g_0_yyyyz_0_xyzzz_1,   \
                             g_0_yyyyz_0_xyzzzz_0,  \
                             g_0_yyyyz_0_xyzzzz_1,  \
                             g_0_yyyyz_0_xzzzz_1,   \
                             g_0_yyyyz_0_xzzzzz_0,  \
                             g_0_yyyyz_0_xzzzzz_1,  \
                             g_0_yyyyz_0_yyyyyy_0,  \
                             g_0_yyyyz_0_yyyyyy_1,  \
                             g_0_yyyyz_0_yyyyyz_0,  \
                             g_0_yyyyz_0_yyyyyz_1,  \
                             g_0_yyyyz_0_yyyyz_1,   \
                             g_0_yyyyz_0_yyyyzz_0,  \
                             g_0_yyyyz_0_yyyyzz_1,  \
                             g_0_yyyyz_0_yyyzz_1,   \
                             g_0_yyyyz_0_yyyzzz_0,  \
                             g_0_yyyyz_0_yyyzzz_1,  \
                             g_0_yyyyz_0_yyzzz_1,   \
                             g_0_yyyyz_0_yyzzzz_0,  \
                             g_0_yyyyz_0_yyzzzz_1,  \
                             g_0_yyyyz_0_yzzzz_1,   \
                             g_0_yyyyz_0_yzzzzz_0,  \
                             g_0_yyyyz_0_yzzzzz_1,  \
                             g_0_yyyyz_0_zzzzz_1,   \
                             g_0_yyyyz_0_zzzzzz_0,  \
                             g_0_yyyyz_0_zzzzzz_1,  \
                             wp_x,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyz_0_xxxxxx_0[i] = g_0_xyyyy_0_xxxxxx_0[i] * pb_z + g_0_xyyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxxy_0[i] = g_0_xyyyy_0_xxxxxy_0[i] * pb_z + g_0_xyyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxxz_0[i] = 5.0 * g_0_yyyyz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxxz_0[i] * pb_x + g_0_yyyyz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxyy_0[i] = g_0_xyyyy_0_xxxxyy_0[i] * pb_z + g_0_xyyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxyz_0[i] = 4.0 * g_0_yyyyz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxyz_0[i] * pb_x + g_0_yyyyz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxzz_0[i] = 4.0 * g_0_yyyyz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxzz_0[i] * pb_x + g_0_yyyyz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxyyy_0[i] = g_0_xyyyy_0_xxxyyy_0[i] * pb_z + g_0_xyyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxyyz_0[i] = 3.0 * g_0_yyyyz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxyyz_0[i] * pb_x + g_0_yyyyz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxyzz_0[i] = 3.0 * g_0_yyyyz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxyzz_0[i] * pb_x + g_0_yyyyz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxzzz_0[i] = 3.0 * g_0_yyyyz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxzzz_0[i] * pb_x + g_0_yyyyz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyyyy_0[i] = g_0_xyyyy_0_xxyyyy_0[i] * pb_z + g_0_xyyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxyyyz_0[i] = 2.0 * g_0_yyyyz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyyyz_0[i] * pb_x + g_0_yyyyz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyyzz_0[i] = 2.0 * g_0_yyyyz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyyzz_0[i] * pb_x + g_0_yyyyz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyzzz_0[i] = 2.0 * g_0_yyyyz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyzzz_0[i] * pb_x + g_0_yyyyz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxzzzz_0[i] = 2.0 * g_0_yyyyz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxzzzz_0[i] * pb_x + g_0_yyyyz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyyyy_0[i] = g_0_xyyyy_0_xyyyyy_0[i] * pb_z + g_0_xyyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xyyyyz_0[i] = g_0_yyyyz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyyyz_0[i] * pb_x + g_0_yyyyz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyyzz_0[i] = g_0_yyyyz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyyzz_0[i] * pb_x + g_0_yyyyz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyzzz_0[i] = g_0_yyyyz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyzzz_0[i] * pb_x + g_0_yyyyz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyzzzz_0[i] = g_0_yyyyz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyzzzz_0[i] * pb_x + g_0_yyyyz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xzzzzz_0[i] = g_0_yyyyz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xzzzzz_0[i] * pb_x + g_0_yyyyz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyyy_0[i] = g_0_yyyyz_0_yyyyyy_0[i] * pb_x + g_0_yyyyz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyyz_0[i] = g_0_yyyyz_0_yyyyyz_0[i] * pb_x + g_0_yyyyz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyzz_0[i] = g_0_yyyyz_0_yyyyzz_0[i] * pb_x + g_0_yyyyz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyzzz_0[i] = g_0_yyyyz_0_yyyzzz_0[i] * pb_x + g_0_yyyyz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyzzzz_0[i] = g_0_yyyyz_0_yyzzzz_0[i] * pb_x + g_0_yyyyz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yzzzzz_0[i] = g_0_yyyyz_0_yzzzzz_0[i] * pb_x + g_0_yyyyz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_zzzzzz_0[i] = g_0_yyyyz_0_zzzzzz_0[i] * pb_x + g_0_yyyyz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 476-504 components of targeted buffer : SISI

    auto g_0_xyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 476);

    auto g_0_xyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 477);

    auto g_0_xyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 478);

    auto g_0_xyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 479);

    auto g_0_xyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 480);

    auto g_0_xyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 481);

    auto g_0_xyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 482);

    auto g_0_xyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 483);

    auto g_0_xyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 484);

    auto g_0_xyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 485);

    auto g_0_xyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 486);

    auto g_0_xyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 487);

    auto g_0_xyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 488);

    auto g_0_xyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 489);

    auto g_0_xyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 490);

    auto g_0_xyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 491);

    auto g_0_xyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 492);

    auto g_0_xyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 493);

    auto g_0_xyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 494);

    auto g_0_xyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 495);

    auto g_0_xyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 496);

    auto g_0_xyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 497);

    auto g_0_xyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 498);

    auto g_0_xyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 499);

    auto g_0_xyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 500);

    auto g_0_xyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 501);

    auto g_0_xyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 502);

    auto g_0_xyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 503);

#pragma omp simd aligned(g_0_xyyyzz_0_xxxxxx_0,     \
                             g_0_xyyyzz_0_xxxxxy_0, \
                             g_0_xyyyzz_0_xxxxxz_0, \
                             g_0_xyyyzz_0_xxxxyy_0, \
                             g_0_xyyyzz_0_xxxxyz_0, \
                             g_0_xyyyzz_0_xxxxzz_0, \
                             g_0_xyyyzz_0_xxxyyy_0, \
                             g_0_xyyyzz_0_xxxyyz_0, \
                             g_0_xyyyzz_0_xxxyzz_0, \
                             g_0_xyyyzz_0_xxxzzz_0, \
                             g_0_xyyyzz_0_xxyyyy_0, \
                             g_0_xyyyzz_0_xxyyyz_0, \
                             g_0_xyyyzz_0_xxyyzz_0, \
                             g_0_xyyyzz_0_xxyzzz_0, \
                             g_0_xyyyzz_0_xxzzzz_0, \
                             g_0_xyyyzz_0_xyyyyy_0, \
                             g_0_xyyyzz_0_xyyyyz_0, \
                             g_0_xyyyzz_0_xyyyzz_0, \
                             g_0_xyyyzz_0_xyyzzz_0, \
                             g_0_xyyyzz_0_xyzzzz_0, \
                             g_0_xyyyzz_0_xzzzzz_0, \
                             g_0_xyyyzz_0_yyyyyy_0, \
                             g_0_xyyyzz_0_yyyyyz_0, \
                             g_0_xyyyzz_0_yyyyzz_0, \
                             g_0_xyyyzz_0_yyyzzz_0, \
                             g_0_xyyyzz_0_yyzzzz_0, \
                             g_0_xyyyzz_0_yzzzzz_0, \
                             g_0_xyyyzz_0_zzzzzz_0, \
                             g_0_yyyzz_0_xxxxx_1,   \
                             g_0_yyyzz_0_xxxxxx_0,  \
                             g_0_yyyzz_0_xxxxxx_1,  \
                             g_0_yyyzz_0_xxxxxy_0,  \
                             g_0_yyyzz_0_xxxxxy_1,  \
                             g_0_yyyzz_0_xxxxxz_0,  \
                             g_0_yyyzz_0_xxxxxz_1,  \
                             g_0_yyyzz_0_xxxxy_1,   \
                             g_0_yyyzz_0_xxxxyy_0,  \
                             g_0_yyyzz_0_xxxxyy_1,  \
                             g_0_yyyzz_0_xxxxyz_0,  \
                             g_0_yyyzz_0_xxxxyz_1,  \
                             g_0_yyyzz_0_xxxxz_1,   \
                             g_0_yyyzz_0_xxxxzz_0,  \
                             g_0_yyyzz_0_xxxxzz_1,  \
                             g_0_yyyzz_0_xxxyy_1,   \
                             g_0_yyyzz_0_xxxyyy_0,  \
                             g_0_yyyzz_0_xxxyyy_1,  \
                             g_0_yyyzz_0_xxxyyz_0,  \
                             g_0_yyyzz_0_xxxyyz_1,  \
                             g_0_yyyzz_0_xxxyz_1,   \
                             g_0_yyyzz_0_xxxyzz_0,  \
                             g_0_yyyzz_0_xxxyzz_1,  \
                             g_0_yyyzz_0_xxxzz_1,   \
                             g_0_yyyzz_0_xxxzzz_0,  \
                             g_0_yyyzz_0_xxxzzz_1,  \
                             g_0_yyyzz_0_xxyyy_1,   \
                             g_0_yyyzz_0_xxyyyy_0,  \
                             g_0_yyyzz_0_xxyyyy_1,  \
                             g_0_yyyzz_0_xxyyyz_0,  \
                             g_0_yyyzz_0_xxyyyz_1,  \
                             g_0_yyyzz_0_xxyyz_1,   \
                             g_0_yyyzz_0_xxyyzz_0,  \
                             g_0_yyyzz_0_xxyyzz_1,  \
                             g_0_yyyzz_0_xxyzz_1,   \
                             g_0_yyyzz_0_xxyzzz_0,  \
                             g_0_yyyzz_0_xxyzzz_1,  \
                             g_0_yyyzz_0_xxzzz_1,   \
                             g_0_yyyzz_0_xxzzzz_0,  \
                             g_0_yyyzz_0_xxzzzz_1,  \
                             g_0_yyyzz_0_xyyyy_1,   \
                             g_0_yyyzz_0_xyyyyy_0,  \
                             g_0_yyyzz_0_xyyyyy_1,  \
                             g_0_yyyzz_0_xyyyyz_0,  \
                             g_0_yyyzz_0_xyyyyz_1,  \
                             g_0_yyyzz_0_xyyyz_1,   \
                             g_0_yyyzz_0_xyyyzz_0,  \
                             g_0_yyyzz_0_xyyyzz_1,  \
                             g_0_yyyzz_0_xyyzz_1,   \
                             g_0_yyyzz_0_xyyzzz_0,  \
                             g_0_yyyzz_0_xyyzzz_1,  \
                             g_0_yyyzz_0_xyzzz_1,   \
                             g_0_yyyzz_0_xyzzzz_0,  \
                             g_0_yyyzz_0_xyzzzz_1,  \
                             g_0_yyyzz_0_xzzzz_1,   \
                             g_0_yyyzz_0_xzzzzz_0,  \
                             g_0_yyyzz_0_xzzzzz_1,  \
                             g_0_yyyzz_0_yyyyy_1,   \
                             g_0_yyyzz_0_yyyyyy_0,  \
                             g_0_yyyzz_0_yyyyyy_1,  \
                             g_0_yyyzz_0_yyyyyz_0,  \
                             g_0_yyyzz_0_yyyyyz_1,  \
                             g_0_yyyzz_0_yyyyz_1,   \
                             g_0_yyyzz_0_yyyyzz_0,  \
                             g_0_yyyzz_0_yyyyzz_1,  \
                             g_0_yyyzz_0_yyyzz_1,   \
                             g_0_yyyzz_0_yyyzzz_0,  \
                             g_0_yyyzz_0_yyyzzz_1,  \
                             g_0_yyyzz_0_yyzzz_1,   \
                             g_0_yyyzz_0_yyzzzz_0,  \
                             g_0_yyyzz_0_yyzzzz_1,  \
                             g_0_yyyzz_0_yzzzz_1,   \
                             g_0_yyyzz_0_yzzzzz_0,  \
                             g_0_yyyzz_0_yzzzzz_1,  \
                             g_0_yyyzz_0_zzzzz_1,   \
                             g_0_yyyzz_0_zzzzzz_0,  \
                             g_0_yyyzz_0_zzzzzz_1,  \
                             wp_x,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzz_0_xxxxxx_0[i] = 6.0 * g_0_yyyzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxx_0[i] * pb_x + g_0_yyyzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxy_0[i] = 5.0 * g_0_yyyzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxy_0[i] * pb_x + g_0_yyyzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxz_0[i] = 5.0 * g_0_yyyzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxz_0[i] * pb_x + g_0_yyyzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxyy_0[i] = 4.0 * g_0_yyyzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyy_0[i] * pb_x + g_0_yyyzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxyz_0[i] = 4.0 * g_0_yyyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyz_0[i] * pb_x + g_0_yyyzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxzz_0[i] = 4.0 * g_0_yyyzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxzz_0[i] * pb_x + g_0_yyyzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyyy_0[i] = 3.0 * g_0_yyyzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyy_0[i] * pb_x + g_0_yyyzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyyz_0[i] = 3.0 * g_0_yyyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyz_0[i] * pb_x + g_0_yyyzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyzz_0[i] = 3.0 * g_0_yyyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyzz_0[i] * pb_x + g_0_yyyzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxzzz_0[i] = 3.0 * g_0_yyyzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxzzz_0[i] * pb_x + g_0_yyyzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyyy_0[i] = 2.0 * g_0_yyyzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyy_0[i] * pb_x + g_0_yyyzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyyz_0[i] = 2.0 * g_0_yyyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyz_0[i] * pb_x + g_0_yyyzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyzz_0[i] = 2.0 * g_0_yyyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyzz_0[i] * pb_x + g_0_yyyzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyzzz_0[i] = 2.0 * g_0_yyyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyzzz_0[i] * pb_x + g_0_yyyzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxzzzz_0[i] = 2.0 * g_0_yyyzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxzzzz_0[i] * pb_x + g_0_yyyzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyyy_0[i] = g_0_yyyzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyy_0[i] * pb_x + g_0_yyyzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyyz_0[i] = g_0_yyyzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyz_0[i] * pb_x + g_0_yyyzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyzz_0[i] = g_0_yyyzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyzz_0[i] * pb_x + g_0_yyyzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyzzz_0[i] = g_0_yyyzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyzzz_0[i] * pb_x + g_0_yyyzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyzzzz_0[i] = g_0_yyyzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyzzzz_0[i] * pb_x + g_0_yyyzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xzzzzz_0[i] = g_0_yyyzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xzzzzz_0[i] * pb_x + g_0_yyyzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyyy_0[i] = g_0_yyyzz_0_yyyyyy_0[i] * pb_x + g_0_yyyzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyyz_0[i] = g_0_yyyzz_0_yyyyyz_0[i] * pb_x + g_0_yyyzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyzz_0[i] = g_0_yyyzz_0_yyyyzz_0[i] * pb_x + g_0_yyyzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyzzz_0[i] = g_0_yyyzz_0_yyyzzz_0[i] * pb_x + g_0_yyyzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyzzzz_0[i] = g_0_yyyzz_0_yyzzzz_0[i] * pb_x + g_0_yyyzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yzzzzz_0[i] = g_0_yyyzz_0_yzzzzz_0[i] * pb_x + g_0_yyyzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_zzzzzz_0[i] = g_0_yyyzz_0_zzzzzz_0[i] * pb_x + g_0_yyyzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 504-532 components of targeted buffer : SISI

    auto g_0_xyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 504);

    auto g_0_xyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 505);

    auto g_0_xyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 506);

    auto g_0_xyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 507);

    auto g_0_xyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 508);

    auto g_0_xyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 509);

    auto g_0_xyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 510);

    auto g_0_xyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 511);

    auto g_0_xyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 512);

    auto g_0_xyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 513);

    auto g_0_xyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 514);

    auto g_0_xyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 515);

    auto g_0_xyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 516);

    auto g_0_xyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 517);

    auto g_0_xyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 518);

    auto g_0_xyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 519);

    auto g_0_xyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 520);

    auto g_0_xyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 521);

    auto g_0_xyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 522);

    auto g_0_xyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 523);

    auto g_0_xyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 524);

    auto g_0_xyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 525);

    auto g_0_xyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 526);

    auto g_0_xyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 527);

    auto g_0_xyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 528);

    auto g_0_xyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 529);

    auto g_0_xyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 530);

    auto g_0_xyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 531);

#pragma omp simd aligned(g_0_xyyzzz_0_xxxxxx_0,     \
                             g_0_xyyzzz_0_xxxxxy_0, \
                             g_0_xyyzzz_0_xxxxxz_0, \
                             g_0_xyyzzz_0_xxxxyy_0, \
                             g_0_xyyzzz_0_xxxxyz_0, \
                             g_0_xyyzzz_0_xxxxzz_0, \
                             g_0_xyyzzz_0_xxxyyy_0, \
                             g_0_xyyzzz_0_xxxyyz_0, \
                             g_0_xyyzzz_0_xxxyzz_0, \
                             g_0_xyyzzz_0_xxxzzz_0, \
                             g_0_xyyzzz_0_xxyyyy_0, \
                             g_0_xyyzzz_0_xxyyyz_0, \
                             g_0_xyyzzz_0_xxyyzz_0, \
                             g_0_xyyzzz_0_xxyzzz_0, \
                             g_0_xyyzzz_0_xxzzzz_0, \
                             g_0_xyyzzz_0_xyyyyy_0, \
                             g_0_xyyzzz_0_xyyyyz_0, \
                             g_0_xyyzzz_0_xyyyzz_0, \
                             g_0_xyyzzz_0_xyyzzz_0, \
                             g_0_xyyzzz_0_xyzzzz_0, \
                             g_0_xyyzzz_0_xzzzzz_0, \
                             g_0_xyyzzz_0_yyyyyy_0, \
                             g_0_xyyzzz_0_yyyyyz_0, \
                             g_0_xyyzzz_0_yyyyzz_0, \
                             g_0_xyyzzz_0_yyyzzz_0, \
                             g_0_xyyzzz_0_yyzzzz_0, \
                             g_0_xyyzzz_0_yzzzzz_0, \
                             g_0_xyyzzz_0_zzzzzz_0, \
                             g_0_yyzzz_0_xxxxx_1,   \
                             g_0_yyzzz_0_xxxxxx_0,  \
                             g_0_yyzzz_0_xxxxxx_1,  \
                             g_0_yyzzz_0_xxxxxy_0,  \
                             g_0_yyzzz_0_xxxxxy_1,  \
                             g_0_yyzzz_0_xxxxxz_0,  \
                             g_0_yyzzz_0_xxxxxz_1,  \
                             g_0_yyzzz_0_xxxxy_1,   \
                             g_0_yyzzz_0_xxxxyy_0,  \
                             g_0_yyzzz_0_xxxxyy_1,  \
                             g_0_yyzzz_0_xxxxyz_0,  \
                             g_0_yyzzz_0_xxxxyz_1,  \
                             g_0_yyzzz_0_xxxxz_1,   \
                             g_0_yyzzz_0_xxxxzz_0,  \
                             g_0_yyzzz_0_xxxxzz_1,  \
                             g_0_yyzzz_0_xxxyy_1,   \
                             g_0_yyzzz_0_xxxyyy_0,  \
                             g_0_yyzzz_0_xxxyyy_1,  \
                             g_0_yyzzz_0_xxxyyz_0,  \
                             g_0_yyzzz_0_xxxyyz_1,  \
                             g_0_yyzzz_0_xxxyz_1,   \
                             g_0_yyzzz_0_xxxyzz_0,  \
                             g_0_yyzzz_0_xxxyzz_1,  \
                             g_0_yyzzz_0_xxxzz_1,   \
                             g_0_yyzzz_0_xxxzzz_0,  \
                             g_0_yyzzz_0_xxxzzz_1,  \
                             g_0_yyzzz_0_xxyyy_1,   \
                             g_0_yyzzz_0_xxyyyy_0,  \
                             g_0_yyzzz_0_xxyyyy_1,  \
                             g_0_yyzzz_0_xxyyyz_0,  \
                             g_0_yyzzz_0_xxyyyz_1,  \
                             g_0_yyzzz_0_xxyyz_1,   \
                             g_0_yyzzz_0_xxyyzz_0,  \
                             g_0_yyzzz_0_xxyyzz_1,  \
                             g_0_yyzzz_0_xxyzz_1,   \
                             g_0_yyzzz_0_xxyzzz_0,  \
                             g_0_yyzzz_0_xxyzzz_1,  \
                             g_0_yyzzz_0_xxzzz_1,   \
                             g_0_yyzzz_0_xxzzzz_0,  \
                             g_0_yyzzz_0_xxzzzz_1,  \
                             g_0_yyzzz_0_xyyyy_1,   \
                             g_0_yyzzz_0_xyyyyy_0,  \
                             g_0_yyzzz_0_xyyyyy_1,  \
                             g_0_yyzzz_0_xyyyyz_0,  \
                             g_0_yyzzz_0_xyyyyz_1,  \
                             g_0_yyzzz_0_xyyyz_1,   \
                             g_0_yyzzz_0_xyyyzz_0,  \
                             g_0_yyzzz_0_xyyyzz_1,  \
                             g_0_yyzzz_0_xyyzz_1,   \
                             g_0_yyzzz_0_xyyzzz_0,  \
                             g_0_yyzzz_0_xyyzzz_1,  \
                             g_0_yyzzz_0_xyzzz_1,   \
                             g_0_yyzzz_0_xyzzzz_0,  \
                             g_0_yyzzz_0_xyzzzz_1,  \
                             g_0_yyzzz_0_xzzzz_1,   \
                             g_0_yyzzz_0_xzzzzz_0,  \
                             g_0_yyzzz_0_xzzzzz_1,  \
                             g_0_yyzzz_0_yyyyy_1,   \
                             g_0_yyzzz_0_yyyyyy_0,  \
                             g_0_yyzzz_0_yyyyyy_1,  \
                             g_0_yyzzz_0_yyyyyz_0,  \
                             g_0_yyzzz_0_yyyyyz_1,  \
                             g_0_yyzzz_0_yyyyz_1,   \
                             g_0_yyzzz_0_yyyyzz_0,  \
                             g_0_yyzzz_0_yyyyzz_1,  \
                             g_0_yyzzz_0_yyyzz_1,   \
                             g_0_yyzzz_0_yyyzzz_0,  \
                             g_0_yyzzz_0_yyyzzz_1,  \
                             g_0_yyzzz_0_yyzzz_1,   \
                             g_0_yyzzz_0_yyzzzz_0,  \
                             g_0_yyzzz_0_yyzzzz_1,  \
                             g_0_yyzzz_0_yzzzz_1,   \
                             g_0_yyzzz_0_yzzzzz_0,  \
                             g_0_yyzzz_0_yzzzzz_1,  \
                             g_0_yyzzz_0_zzzzz_1,   \
                             g_0_yyzzz_0_zzzzzz_0,  \
                             g_0_yyzzz_0_zzzzzz_1,  \
                             wp_x,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzz_0_xxxxxx_0[i] = 6.0 * g_0_yyzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxx_0[i] * pb_x + g_0_yyzzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxy_0[i] = 5.0 * g_0_yyzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxy_0[i] * pb_x + g_0_yyzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxz_0[i] = 5.0 * g_0_yyzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxz_0[i] * pb_x + g_0_yyzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxyy_0[i] = 4.0 * g_0_yyzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyy_0[i] * pb_x + g_0_yyzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxyz_0[i] = 4.0 * g_0_yyzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyz_0[i] * pb_x + g_0_yyzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxzz_0[i] = 4.0 * g_0_yyzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxzz_0[i] * pb_x + g_0_yyzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyyy_0[i] = 3.0 * g_0_yyzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyy_0[i] * pb_x + g_0_yyzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyyz_0[i] = 3.0 * g_0_yyzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyz_0[i] * pb_x + g_0_yyzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyzz_0[i] = 3.0 * g_0_yyzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyzz_0[i] * pb_x + g_0_yyzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxzzz_0[i] = 3.0 * g_0_yyzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxzzz_0[i] * pb_x + g_0_yyzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyyy_0[i] = 2.0 * g_0_yyzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyy_0[i] * pb_x + g_0_yyzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyyz_0[i] = 2.0 * g_0_yyzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyz_0[i] * pb_x + g_0_yyzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyzz_0[i] = 2.0 * g_0_yyzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyzz_0[i] * pb_x + g_0_yyzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyzzz_0[i] = 2.0 * g_0_yyzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyzzz_0[i] * pb_x + g_0_yyzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxzzzz_0[i] = 2.0 * g_0_yyzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxzzzz_0[i] * pb_x + g_0_yyzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyyy_0[i] = g_0_yyzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyy_0[i] * pb_x + g_0_yyzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyyz_0[i] = g_0_yyzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyz_0[i] * pb_x + g_0_yyzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyzz_0[i] = g_0_yyzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyzz_0[i] * pb_x + g_0_yyzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyzzz_0[i] = g_0_yyzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyzzz_0[i] * pb_x + g_0_yyzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyzzzz_0[i] = g_0_yyzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyzzzz_0[i] * pb_x + g_0_yyzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xzzzzz_0[i] = g_0_yyzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xzzzzz_0[i] * pb_x + g_0_yyzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyyy_0[i] = g_0_yyzzz_0_yyyyyy_0[i] * pb_x + g_0_yyzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyyz_0[i] = g_0_yyzzz_0_yyyyyz_0[i] * pb_x + g_0_yyzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyzz_0[i] = g_0_yyzzz_0_yyyyzz_0[i] * pb_x + g_0_yyzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyzzz_0[i] = g_0_yyzzz_0_yyyzzz_0[i] * pb_x + g_0_yyzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyzzzz_0[i] = g_0_yyzzz_0_yyzzzz_0[i] * pb_x + g_0_yyzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yzzzzz_0[i] = g_0_yyzzz_0_yzzzzz_0[i] * pb_x + g_0_yyzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_zzzzzz_0[i] = g_0_yyzzz_0_zzzzzz_0[i] * pb_x + g_0_yyzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 532-560 components of targeted buffer : SISI

    auto g_0_xyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 532);

    auto g_0_xyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 533);

    auto g_0_xyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 534);

    auto g_0_xyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 535);

    auto g_0_xyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 536);

    auto g_0_xyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 537);

    auto g_0_xyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 538);

    auto g_0_xyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 539);

    auto g_0_xyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 540);

    auto g_0_xyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 541);

    auto g_0_xyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 542);

    auto g_0_xyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 543);

    auto g_0_xyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 544);

    auto g_0_xyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 545);

    auto g_0_xyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 546);

    auto g_0_xyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 547);

    auto g_0_xyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 548);

    auto g_0_xyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 549);

    auto g_0_xyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 550);

    auto g_0_xyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 551);

    auto g_0_xyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 552);

    auto g_0_xyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 553);

    auto g_0_xyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 554);

    auto g_0_xyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 555);

    auto g_0_xyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 556);

    auto g_0_xyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 557);

    auto g_0_xyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 558);

    auto g_0_xyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 559);

#pragma omp simd aligned(g_0_xyzzzz_0_xxxxxx_0,     \
                             g_0_xyzzzz_0_xxxxxy_0, \
                             g_0_xyzzzz_0_xxxxxz_0, \
                             g_0_xyzzzz_0_xxxxyy_0, \
                             g_0_xyzzzz_0_xxxxyz_0, \
                             g_0_xyzzzz_0_xxxxzz_0, \
                             g_0_xyzzzz_0_xxxyyy_0, \
                             g_0_xyzzzz_0_xxxyyz_0, \
                             g_0_xyzzzz_0_xxxyzz_0, \
                             g_0_xyzzzz_0_xxxzzz_0, \
                             g_0_xyzzzz_0_xxyyyy_0, \
                             g_0_xyzzzz_0_xxyyyz_0, \
                             g_0_xyzzzz_0_xxyyzz_0, \
                             g_0_xyzzzz_0_xxyzzz_0, \
                             g_0_xyzzzz_0_xxzzzz_0, \
                             g_0_xyzzzz_0_xyyyyy_0, \
                             g_0_xyzzzz_0_xyyyyz_0, \
                             g_0_xyzzzz_0_xyyyzz_0, \
                             g_0_xyzzzz_0_xyyzzz_0, \
                             g_0_xyzzzz_0_xyzzzz_0, \
                             g_0_xyzzzz_0_xzzzzz_0, \
                             g_0_xyzzzz_0_yyyyyy_0, \
                             g_0_xyzzzz_0_yyyyyz_0, \
                             g_0_xyzzzz_0_yyyyzz_0, \
                             g_0_xyzzzz_0_yyyzzz_0, \
                             g_0_xyzzzz_0_yyzzzz_0, \
                             g_0_xyzzzz_0_yzzzzz_0, \
                             g_0_xyzzzz_0_zzzzzz_0, \
                             g_0_xzzzz_0_xxxxxx_0,  \
                             g_0_xzzzz_0_xxxxxx_1,  \
                             g_0_xzzzz_0_xxxxxz_0,  \
                             g_0_xzzzz_0_xxxxxz_1,  \
                             g_0_xzzzz_0_xxxxzz_0,  \
                             g_0_xzzzz_0_xxxxzz_1,  \
                             g_0_xzzzz_0_xxxzzz_0,  \
                             g_0_xzzzz_0_xxxzzz_1,  \
                             g_0_xzzzz_0_xxzzzz_0,  \
                             g_0_xzzzz_0_xxzzzz_1,  \
                             g_0_xzzzz_0_xzzzzz_0,  \
                             g_0_xzzzz_0_xzzzzz_1,  \
                             g_0_yzzzz_0_xxxxxy_0,  \
                             g_0_yzzzz_0_xxxxxy_1,  \
                             g_0_yzzzz_0_xxxxy_1,   \
                             g_0_yzzzz_0_xxxxyy_0,  \
                             g_0_yzzzz_0_xxxxyy_1,  \
                             g_0_yzzzz_0_xxxxyz_0,  \
                             g_0_yzzzz_0_xxxxyz_1,  \
                             g_0_yzzzz_0_xxxyy_1,   \
                             g_0_yzzzz_0_xxxyyy_0,  \
                             g_0_yzzzz_0_xxxyyy_1,  \
                             g_0_yzzzz_0_xxxyyz_0,  \
                             g_0_yzzzz_0_xxxyyz_1,  \
                             g_0_yzzzz_0_xxxyz_1,   \
                             g_0_yzzzz_0_xxxyzz_0,  \
                             g_0_yzzzz_0_xxxyzz_1,  \
                             g_0_yzzzz_0_xxyyy_1,   \
                             g_0_yzzzz_0_xxyyyy_0,  \
                             g_0_yzzzz_0_xxyyyy_1,  \
                             g_0_yzzzz_0_xxyyyz_0,  \
                             g_0_yzzzz_0_xxyyyz_1,  \
                             g_0_yzzzz_0_xxyyz_1,   \
                             g_0_yzzzz_0_xxyyzz_0,  \
                             g_0_yzzzz_0_xxyyzz_1,  \
                             g_0_yzzzz_0_xxyzz_1,   \
                             g_0_yzzzz_0_xxyzzz_0,  \
                             g_0_yzzzz_0_xxyzzz_1,  \
                             g_0_yzzzz_0_xyyyy_1,   \
                             g_0_yzzzz_0_xyyyyy_0,  \
                             g_0_yzzzz_0_xyyyyy_1,  \
                             g_0_yzzzz_0_xyyyyz_0,  \
                             g_0_yzzzz_0_xyyyyz_1,  \
                             g_0_yzzzz_0_xyyyz_1,   \
                             g_0_yzzzz_0_xyyyzz_0,  \
                             g_0_yzzzz_0_xyyyzz_1,  \
                             g_0_yzzzz_0_xyyzz_1,   \
                             g_0_yzzzz_0_xyyzzz_0,  \
                             g_0_yzzzz_0_xyyzzz_1,  \
                             g_0_yzzzz_0_xyzzz_1,   \
                             g_0_yzzzz_0_xyzzzz_0,  \
                             g_0_yzzzz_0_xyzzzz_1,  \
                             g_0_yzzzz_0_yyyyy_1,   \
                             g_0_yzzzz_0_yyyyyy_0,  \
                             g_0_yzzzz_0_yyyyyy_1,  \
                             g_0_yzzzz_0_yyyyyz_0,  \
                             g_0_yzzzz_0_yyyyyz_1,  \
                             g_0_yzzzz_0_yyyyz_1,   \
                             g_0_yzzzz_0_yyyyzz_0,  \
                             g_0_yzzzz_0_yyyyzz_1,  \
                             g_0_yzzzz_0_yyyzz_1,   \
                             g_0_yzzzz_0_yyyzzz_0,  \
                             g_0_yzzzz_0_yyyzzz_1,  \
                             g_0_yzzzz_0_yyzzz_1,   \
                             g_0_yzzzz_0_yyzzzz_0,  \
                             g_0_yzzzz_0_yyzzzz_1,  \
                             g_0_yzzzz_0_yzzzz_1,   \
                             g_0_yzzzz_0_yzzzzz_0,  \
                             g_0_yzzzz_0_yzzzzz_1,  \
                             g_0_yzzzz_0_zzzzzz_0,  \
                             g_0_yzzzz_0_zzzzzz_1,  \
                             wp_x,                  \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzz_0_xxxxxx_0[i] = g_0_xzzzz_0_xxxxxx_0[i] * pb_y + g_0_xzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxxxy_0[i] = 5.0 * g_0_yzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxxy_0[i] * pb_x + g_0_yzzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxxz_0[i] = g_0_xzzzz_0_xxxxxz_0[i] * pb_y + g_0_xzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxxyy_0[i] = 4.0 * g_0_yzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyy_0[i] * pb_x + g_0_yzzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxyz_0[i] = 4.0 * g_0_yzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyz_0[i] * pb_x + g_0_yzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxzz_0[i] = g_0_xzzzz_0_xxxxzz_0[i] * pb_y + g_0_xzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxyyy_0[i] = 3.0 * g_0_yzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyy_0[i] * pb_x + g_0_yzzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxyyz_0[i] = 3.0 * g_0_yzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyz_0[i] * pb_x + g_0_yzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxyzz_0[i] = 3.0 * g_0_yzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyzz_0[i] * pb_x + g_0_yzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxzzz_0[i] = g_0_xzzzz_0_xxxzzz_0[i] * pb_y + g_0_xzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxyyyy_0[i] = 2.0 * g_0_yzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyy_0[i] * pb_x + g_0_yzzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyyyz_0[i] = 2.0 * g_0_yzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyz_0[i] * pb_x + g_0_yzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyyzz_0[i] = 2.0 * g_0_yzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyzz_0[i] * pb_x + g_0_yzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyzzz_0[i] = 2.0 * g_0_yzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyzzz_0[i] * pb_x + g_0_yzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxzzzz_0[i] = g_0_xzzzz_0_xxzzzz_0[i] * pb_y + g_0_xzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xyyyyy_0[i] = g_0_yzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyy_0[i] * pb_x + g_0_yzzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyyyz_0[i] = g_0_yzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyz_0[i] * pb_x + g_0_yzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyyzz_0[i] = g_0_yzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyzz_0[i] * pb_x + g_0_yzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyzzz_0[i] = g_0_yzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyzzz_0[i] * pb_x + g_0_yzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyzzzz_0[i] = g_0_yzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyzzzz_0[i] * pb_x + g_0_yzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xzzzzz_0[i] = g_0_xzzzz_0_xzzzzz_0[i] * pb_y + g_0_xzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_yyyyyy_0[i] = g_0_yzzzz_0_yyyyyy_0[i] * pb_x + g_0_yzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyyyz_0[i] = g_0_yzzzz_0_yyyyyz_0[i] * pb_x + g_0_yzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyyzz_0[i] = g_0_yzzzz_0_yyyyzz_0[i] * pb_x + g_0_yzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyzzz_0[i] = g_0_yzzzz_0_yyyzzz_0[i] * pb_x + g_0_yzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyzzzz_0[i] = g_0_yzzzz_0_yyzzzz_0[i] * pb_x + g_0_yzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yzzzzz_0[i] = g_0_yzzzz_0_yzzzzz_0[i] * pb_x + g_0_yzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_zzzzzz_0[i] = g_0_yzzzz_0_zzzzzz_0[i] * pb_x + g_0_yzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 560-588 components of targeted buffer : SISI

    auto g_0_xzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 560);

    auto g_0_xzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 561);

    auto g_0_xzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 562);

    auto g_0_xzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 563);

    auto g_0_xzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 564);

    auto g_0_xzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 565);

    auto g_0_xzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 566);

    auto g_0_xzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 567);

    auto g_0_xzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 568);

    auto g_0_xzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 569);

    auto g_0_xzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 570);

    auto g_0_xzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 571);

    auto g_0_xzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 572);

    auto g_0_xzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 573);

    auto g_0_xzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 574);

    auto g_0_xzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 575);

    auto g_0_xzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 576);

    auto g_0_xzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 577);

    auto g_0_xzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 578);

    auto g_0_xzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 579);

    auto g_0_xzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 580);

    auto g_0_xzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 581);

    auto g_0_xzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 582);

    auto g_0_xzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 583);

    auto g_0_xzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 584);

    auto g_0_xzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 585);

    auto g_0_xzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 586);

    auto g_0_xzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 587);

#pragma omp simd aligned(g_0_xzzzzz_0_xxxxxx_0,     \
                             g_0_xzzzzz_0_xxxxxy_0, \
                             g_0_xzzzzz_0_xxxxxz_0, \
                             g_0_xzzzzz_0_xxxxyy_0, \
                             g_0_xzzzzz_0_xxxxyz_0, \
                             g_0_xzzzzz_0_xxxxzz_0, \
                             g_0_xzzzzz_0_xxxyyy_0, \
                             g_0_xzzzzz_0_xxxyyz_0, \
                             g_0_xzzzzz_0_xxxyzz_0, \
                             g_0_xzzzzz_0_xxxzzz_0, \
                             g_0_xzzzzz_0_xxyyyy_0, \
                             g_0_xzzzzz_0_xxyyyz_0, \
                             g_0_xzzzzz_0_xxyyzz_0, \
                             g_0_xzzzzz_0_xxyzzz_0, \
                             g_0_xzzzzz_0_xxzzzz_0, \
                             g_0_xzzzzz_0_xyyyyy_0, \
                             g_0_xzzzzz_0_xyyyyz_0, \
                             g_0_xzzzzz_0_xyyyzz_0, \
                             g_0_xzzzzz_0_xyyzzz_0, \
                             g_0_xzzzzz_0_xyzzzz_0, \
                             g_0_xzzzzz_0_xzzzzz_0, \
                             g_0_xzzzzz_0_yyyyyy_0, \
                             g_0_xzzzzz_0_yyyyyz_0, \
                             g_0_xzzzzz_0_yyyyzz_0, \
                             g_0_xzzzzz_0_yyyzzz_0, \
                             g_0_xzzzzz_0_yyzzzz_0, \
                             g_0_xzzzzz_0_yzzzzz_0, \
                             g_0_xzzzzz_0_zzzzzz_0, \
                             g_0_zzzzz_0_xxxxx_1,   \
                             g_0_zzzzz_0_xxxxxx_0,  \
                             g_0_zzzzz_0_xxxxxx_1,  \
                             g_0_zzzzz_0_xxxxxy_0,  \
                             g_0_zzzzz_0_xxxxxy_1,  \
                             g_0_zzzzz_0_xxxxxz_0,  \
                             g_0_zzzzz_0_xxxxxz_1,  \
                             g_0_zzzzz_0_xxxxy_1,   \
                             g_0_zzzzz_0_xxxxyy_0,  \
                             g_0_zzzzz_0_xxxxyy_1,  \
                             g_0_zzzzz_0_xxxxyz_0,  \
                             g_0_zzzzz_0_xxxxyz_1,  \
                             g_0_zzzzz_0_xxxxz_1,   \
                             g_0_zzzzz_0_xxxxzz_0,  \
                             g_0_zzzzz_0_xxxxzz_1,  \
                             g_0_zzzzz_0_xxxyy_1,   \
                             g_0_zzzzz_0_xxxyyy_0,  \
                             g_0_zzzzz_0_xxxyyy_1,  \
                             g_0_zzzzz_0_xxxyyz_0,  \
                             g_0_zzzzz_0_xxxyyz_1,  \
                             g_0_zzzzz_0_xxxyz_1,   \
                             g_0_zzzzz_0_xxxyzz_0,  \
                             g_0_zzzzz_0_xxxyzz_1,  \
                             g_0_zzzzz_0_xxxzz_1,   \
                             g_0_zzzzz_0_xxxzzz_0,  \
                             g_0_zzzzz_0_xxxzzz_1,  \
                             g_0_zzzzz_0_xxyyy_1,   \
                             g_0_zzzzz_0_xxyyyy_0,  \
                             g_0_zzzzz_0_xxyyyy_1,  \
                             g_0_zzzzz_0_xxyyyz_0,  \
                             g_0_zzzzz_0_xxyyyz_1,  \
                             g_0_zzzzz_0_xxyyz_1,   \
                             g_0_zzzzz_0_xxyyzz_0,  \
                             g_0_zzzzz_0_xxyyzz_1,  \
                             g_0_zzzzz_0_xxyzz_1,   \
                             g_0_zzzzz_0_xxyzzz_0,  \
                             g_0_zzzzz_0_xxyzzz_1,  \
                             g_0_zzzzz_0_xxzzz_1,   \
                             g_0_zzzzz_0_xxzzzz_0,  \
                             g_0_zzzzz_0_xxzzzz_1,  \
                             g_0_zzzzz_0_xyyyy_1,   \
                             g_0_zzzzz_0_xyyyyy_0,  \
                             g_0_zzzzz_0_xyyyyy_1,  \
                             g_0_zzzzz_0_xyyyyz_0,  \
                             g_0_zzzzz_0_xyyyyz_1,  \
                             g_0_zzzzz_0_xyyyz_1,   \
                             g_0_zzzzz_0_xyyyzz_0,  \
                             g_0_zzzzz_0_xyyyzz_1,  \
                             g_0_zzzzz_0_xyyzz_1,   \
                             g_0_zzzzz_0_xyyzzz_0,  \
                             g_0_zzzzz_0_xyyzzz_1,  \
                             g_0_zzzzz_0_xyzzz_1,   \
                             g_0_zzzzz_0_xyzzzz_0,  \
                             g_0_zzzzz_0_xyzzzz_1,  \
                             g_0_zzzzz_0_xzzzz_1,   \
                             g_0_zzzzz_0_xzzzzz_0,  \
                             g_0_zzzzz_0_xzzzzz_1,  \
                             g_0_zzzzz_0_yyyyy_1,   \
                             g_0_zzzzz_0_yyyyyy_0,  \
                             g_0_zzzzz_0_yyyyyy_1,  \
                             g_0_zzzzz_0_yyyyyz_0,  \
                             g_0_zzzzz_0_yyyyyz_1,  \
                             g_0_zzzzz_0_yyyyz_1,   \
                             g_0_zzzzz_0_yyyyzz_0,  \
                             g_0_zzzzz_0_yyyyzz_1,  \
                             g_0_zzzzz_0_yyyzz_1,   \
                             g_0_zzzzz_0_yyyzzz_0,  \
                             g_0_zzzzz_0_yyyzzz_1,  \
                             g_0_zzzzz_0_yyzzz_1,   \
                             g_0_zzzzz_0_yyzzzz_0,  \
                             g_0_zzzzz_0_yyzzzz_1,  \
                             g_0_zzzzz_0_yzzzz_1,   \
                             g_0_zzzzz_0_yzzzzz_0,  \
                             g_0_zzzzz_0_yzzzzz_1,  \
                             g_0_zzzzz_0_zzzzz_1,   \
                             g_0_zzzzz_0_zzzzzz_0,  \
                             g_0_zzzzz_0_zzzzzz_1,  \
                             wp_x,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzz_0_xxxxxx_0[i] = 6.0 * g_0_zzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxx_0[i] * pb_x + g_0_zzzzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxy_0[i] = 5.0 * g_0_zzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxy_0[i] * pb_x + g_0_zzzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxz_0[i] = 5.0 * g_0_zzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxz_0[i] * pb_x + g_0_zzzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxyy_0[i] = 4.0 * g_0_zzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyy_0[i] * pb_x + g_0_zzzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxyz_0[i] = 4.0 * g_0_zzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyz_0[i] * pb_x + g_0_zzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxzz_0[i] = 4.0 * g_0_zzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxzz_0[i] * pb_x + g_0_zzzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyyy_0[i] = 3.0 * g_0_zzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyy_0[i] * pb_x + g_0_zzzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyyz_0[i] = 3.0 * g_0_zzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyz_0[i] * pb_x + g_0_zzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyzz_0[i] = 3.0 * g_0_zzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyzz_0[i] * pb_x + g_0_zzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxzzz_0[i] = 3.0 * g_0_zzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxzzz_0[i] * pb_x + g_0_zzzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyyy_0[i] = 2.0 * g_0_zzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyy_0[i] * pb_x + g_0_zzzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyyz_0[i] = 2.0 * g_0_zzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyz_0[i] * pb_x + g_0_zzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyzz_0[i] = 2.0 * g_0_zzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyzz_0[i] * pb_x + g_0_zzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyzzz_0[i] = 2.0 * g_0_zzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyzzz_0[i] * pb_x + g_0_zzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxzzzz_0[i] = 2.0 * g_0_zzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxzzzz_0[i] * pb_x + g_0_zzzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyyy_0[i] = g_0_zzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyy_0[i] * pb_x + g_0_zzzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyyz_0[i] = g_0_zzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyz_0[i] * pb_x + g_0_zzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyzz_0[i] = g_0_zzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyzz_0[i] * pb_x + g_0_zzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyzzz_0[i] = g_0_zzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyzzz_0[i] * pb_x + g_0_zzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyzzzz_0[i] = g_0_zzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzzzz_0[i] * pb_x + g_0_zzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xzzzzz_0[i] = g_0_zzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xzzzzz_0[i] * pb_x + g_0_zzzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyyy_0[i] = g_0_zzzzz_0_yyyyyy_0[i] * pb_x + g_0_zzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyyz_0[i] = g_0_zzzzz_0_yyyyyz_0[i] * pb_x + g_0_zzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyzz_0[i] = g_0_zzzzz_0_yyyyzz_0[i] * pb_x + g_0_zzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyzzz_0[i] = g_0_zzzzz_0_yyyzzz_0[i] * pb_x + g_0_zzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyzzzz_0[i] = g_0_zzzzz_0_yyzzzz_0[i] * pb_x + g_0_zzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yzzzzz_0[i] = g_0_zzzzz_0_yzzzzz_0[i] * pb_x + g_0_zzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_zzzzzz_0[i] = g_0_zzzzz_0_zzzzzz_0[i] * pb_x + g_0_zzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 588-616 components of targeted buffer : SISI

    auto g_0_yyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 588);

    auto g_0_yyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 589);

    auto g_0_yyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 590);

    auto g_0_yyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 591);

    auto g_0_yyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 592);

    auto g_0_yyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 593);

    auto g_0_yyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 594);

    auto g_0_yyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 595);

    auto g_0_yyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 596);

    auto g_0_yyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 597);

    auto g_0_yyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 598);

    auto g_0_yyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 599);

    auto g_0_yyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 600);

    auto g_0_yyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 601);

    auto g_0_yyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 602);

    auto g_0_yyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 603);

    auto g_0_yyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 604);

    auto g_0_yyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 605);

    auto g_0_yyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 606);

    auto g_0_yyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 607);

    auto g_0_yyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 608);

    auto g_0_yyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 609);

    auto g_0_yyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 610);

    auto g_0_yyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 611);

    auto g_0_yyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 612);

    auto g_0_yyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 613);

    auto g_0_yyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 614);

    auto g_0_yyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 615);

#pragma omp simd aligned(g_0_yyyy_0_xxxxxx_0,       \
                             g_0_yyyy_0_xxxxxx_1,   \
                             g_0_yyyy_0_xxxxxy_0,   \
                             g_0_yyyy_0_xxxxxy_1,   \
                             g_0_yyyy_0_xxxxxz_0,   \
                             g_0_yyyy_0_xxxxxz_1,   \
                             g_0_yyyy_0_xxxxyy_0,   \
                             g_0_yyyy_0_xxxxyy_1,   \
                             g_0_yyyy_0_xxxxyz_0,   \
                             g_0_yyyy_0_xxxxyz_1,   \
                             g_0_yyyy_0_xxxxzz_0,   \
                             g_0_yyyy_0_xxxxzz_1,   \
                             g_0_yyyy_0_xxxyyy_0,   \
                             g_0_yyyy_0_xxxyyy_1,   \
                             g_0_yyyy_0_xxxyyz_0,   \
                             g_0_yyyy_0_xxxyyz_1,   \
                             g_0_yyyy_0_xxxyzz_0,   \
                             g_0_yyyy_0_xxxyzz_1,   \
                             g_0_yyyy_0_xxxzzz_0,   \
                             g_0_yyyy_0_xxxzzz_1,   \
                             g_0_yyyy_0_xxyyyy_0,   \
                             g_0_yyyy_0_xxyyyy_1,   \
                             g_0_yyyy_0_xxyyyz_0,   \
                             g_0_yyyy_0_xxyyyz_1,   \
                             g_0_yyyy_0_xxyyzz_0,   \
                             g_0_yyyy_0_xxyyzz_1,   \
                             g_0_yyyy_0_xxyzzz_0,   \
                             g_0_yyyy_0_xxyzzz_1,   \
                             g_0_yyyy_0_xxzzzz_0,   \
                             g_0_yyyy_0_xxzzzz_1,   \
                             g_0_yyyy_0_xyyyyy_0,   \
                             g_0_yyyy_0_xyyyyy_1,   \
                             g_0_yyyy_0_xyyyyz_0,   \
                             g_0_yyyy_0_xyyyyz_1,   \
                             g_0_yyyy_0_xyyyzz_0,   \
                             g_0_yyyy_0_xyyyzz_1,   \
                             g_0_yyyy_0_xyyzzz_0,   \
                             g_0_yyyy_0_xyyzzz_1,   \
                             g_0_yyyy_0_xyzzzz_0,   \
                             g_0_yyyy_0_xyzzzz_1,   \
                             g_0_yyyy_0_xzzzzz_0,   \
                             g_0_yyyy_0_xzzzzz_1,   \
                             g_0_yyyy_0_yyyyyy_0,   \
                             g_0_yyyy_0_yyyyyy_1,   \
                             g_0_yyyy_0_yyyyyz_0,   \
                             g_0_yyyy_0_yyyyyz_1,   \
                             g_0_yyyy_0_yyyyzz_0,   \
                             g_0_yyyy_0_yyyyzz_1,   \
                             g_0_yyyy_0_yyyzzz_0,   \
                             g_0_yyyy_0_yyyzzz_1,   \
                             g_0_yyyy_0_yyzzzz_0,   \
                             g_0_yyyy_0_yyzzzz_1,   \
                             g_0_yyyy_0_yzzzzz_0,   \
                             g_0_yyyy_0_yzzzzz_1,   \
                             g_0_yyyy_0_zzzzzz_0,   \
                             g_0_yyyy_0_zzzzzz_1,   \
                             g_0_yyyyy_0_xxxxx_1,   \
                             g_0_yyyyy_0_xxxxxx_0,  \
                             g_0_yyyyy_0_xxxxxx_1,  \
                             g_0_yyyyy_0_xxxxxy_0,  \
                             g_0_yyyyy_0_xxxxxy_1,  \
                             g_0_yyyyy_0_xxxxxz_0,  \
                             g_0_yyyyy_0_xxxxxz_1,  \
                             g_0_yyyyy_0_xxxxy_1,   \
                             g_0_yyyyy_0_xxxxyy_0,  \
                             g_0_yyyyy_0_xxxxyy_1,  \
                             g_0_yyyyy_0_xxxxyz_0,  \
                             g_0_yyyyy_0_xxxxyz_1,  \
                             g_0_yyyyy_0_xxxxz_1,   \
                             g_0_yyyyy_0_xxxxzz_0,  \
                             g_0_yyyyy_0_xxxxzz_1,  \
                             g_0_yyyyy_0_xxxyy_1,   \
                             g_0_yyyyy_0_xxxyyy_0,  \
                             g_0_yyyyy_0_xxxyyy_1,  \
                             g_0_yyyyy_0_xxxyyz_0,  \
                             g_0_yyyyy_0_xxxyyz_1,  \
                             g_0_yyyyy_0_xxxyz_1,   \
                             g_0_yyyyy_0_xxxyzz_0,  \
                             g_0_yyyyy_0_xxxyzz_1,  \
                             g_0_yyyyy_0_xxxzz_1,   \
                             g_0_yyyyy_0_xxxzzz_0,  \
                             g_0_yyyyy_0_xxxzzz_1,  \
                             g_0_yyyyy_0_xxyyy_1,   \
                             g_0_yyyyy_0_xxyyyy_0,  \
                             g_0_yyyyy_0_xxyyyy_1,  \
                             g_0_yyyyy_0_xxyyyz_0,  \
                             g_0_yyyyy_0_xxyyyz_1,  \
                             g_0_yyyyy_0_xxyyz_1,   \
                             g_0_yyyyy_0_xxyyzz_0,  \
                             g_0_yyyyy_0_xxyyzz_1,  \
                             g_0_yyyyy_0_xxyzz_1,   \
                             g_0_yyyyy_0_xxyzzz_0,  \
                             g_0_yyyyy_0_xxyzzz_1,  \
                             g_0_yyyyy_0_xxzzz_1,   \
                             g_0_yyyyy_0_xxzzzz_0,  \
                             g_0_yyyyy_0_xxzzzz_1,  \
                             g_0_yyyyy_0_xyyyy_1,   \
                             g_0_yyyyy_0_xyyyyy_0,  \
                             g_0_yyyyy_0_xyyyyy_1,  \
                             g_0_yyyyy_0_xyyyyz_0,  \
                             g_0_yyyyy_0_xyyyyz_1,  \
                             g_0_yyyyy_0_xyyyz_1,   \
                             g_0_yyyyy_0_xyyyzz_0,  \
                             g_0_yyyyy_0_xyyyzz_1,  \
                             g_0_yyyyy_0_xyyzz_1,   \
                             g_0_yyyyy_0_xyyzzz_0,  \
                             g_0_yyyyy_0_xyyzzz_1,  \
                             g_0_yyyyy_0_xyzzz_1,   \
                             g_0_yyyyy_0_xyzzzz_0,  \
                             g_0_yyyyy_0_xyzzzz_1,  \
                             g_0_yyyyy_0_xzzzz_1,   \
                             g_0_yyyyy_0_xzzzzz_0,  \
                             g_0_yyyyy_0_xzzzzz_1,  \
                             g_0_yyyyy_0_yyyyy_1,   \
                             g_0_yyyyy_0_yyyyyy_0,  \
                             g_0_yyyyy_0_yyyyyy_1,  \
                             g_0_yyyyy_0_yyyyyz_0,  \
                             g_0_yyyyy_0_yyyyyz_1,  \
                             g_0_yyyyy_0_yyyyz_1,   \
                             g_0_yyyyy_0_yyyyzz_0,  \
                             g_0_yyyyy_0_yyyyzz_1,  \
                             g_0_yyyyy_0_yyyzz_1,   \
                             g_0_yyyyy_0_yyyzzz_0,  \
                             g_0_yyyyy_0_yyyzzz_1,  \
                             g_0_yyyyy_0_yyzzz_1,   \
                             g_0_yyyyy_0_yyzzzz_0,  \
                             g_0_yyyyy_0_yyzzzz_1,  \
                             g_0_yyyyy_0_yzzzz_1,   \
                             g_0_yyyyy_0_yzzzzz_0,  \
                             g_0_yyyyy_0_yzzzzz_1,  \
                             g_0_yyyyy_0_zzzzz_1,   \
                             g_0_yyyyy_0_zzzzzz_0,  \
                             g_0_yyyyy_0_zzzzzz_1,  \
                             g_0_yyyyyy_0_xxxxxx_0, \
                             g_0_yyyyyy_0_xxxxxy_0, \
                             g_0_yyyyyy_0_xxxxxz_0, \
                             g_0_yyyyyy_0_xxxxyy_0, \
                             g_0_yyyyyy_0_xxxxyz_0, \
                             g_0_yyyyyy_0_xxxxzz_0, \
                             g_0_yyyyyy_0_xxxyyy_0, \
                             g_0_yyyyyy_0_xxxyyz_0, \
                             g_0_yyyyyy_0_xxxyzz_0, \
                             g_0_yyyyyy_0_xxxzzz_0, \
                             g_0_yyyyyy_0_xxyyyy_0, \
                             g_0_yyyyyy_0_xxyyyz_0, \
                             g_0_yyyyyy_0_xxyyzz_0, \
                             g_0_yyyyyy_0_xxyzzz_0, \
                             g_0_yyyyyy_0_xxzzzz_0, \
                             g_0_yyyyyy_0_xyyyyy_0, \
                             g_0_yyyyyy_0_xyyyyz_0, \
                             g_0_yyyyyy_0_xyyyzz_0, \
                             g_0_yyyyyy_0_xyyzzz_0, \
                             g_0_yyyyyy_0_xyzzzz_0, \
                             g_0_yyyyyy_0_xzzzzz_0, \
                             g_0_yyyyyy_0_yyyyyy_0, \
                             g_0_yyyyyy_0_yyyyyz_0, \
                             g_0_yyyyyy_0_yyyyzz_0, \
                             g_0_yyyyyy_0_yyyzzz_0, \
                             g_0_yyyyyy_0_yyzzzz_0, \
                             g_0_yyyyyy_0_yzzzzz_0, \
                             g_0_yyyyyy_0_zzzzzz_0, \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyy_0_xxxxxx_0[i] = 5.0 * g_0_yyyy_0_xxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxx_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxxxx_0[i] * pb_y +
                                   g_0_yyyyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxy_0[i] = 5.0 * g_0_yyyy_0_xxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxy_1[i] * fti_ab_0 +
                                   g_0_yyyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxy_0[i] * pb_y + g_0_yyyyy_0_xxxxxy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxz_0[i] = 5.0 * g_0_yyyy_0_xxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxxxz_0[i] * pb_y +
                                   g_0_yyyyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxyy_0[i] = 5.0 * g_0_yyyy_0_xxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxyy_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyy_0[i] * pb_y + g_0_yyyyy_0_xxxxyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxyz_0[i] = 5.0 * g_0_yyyy_0_xxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxyz_1[i] * fti_ab_0 +
                                   g_0_yyyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyz_0[i] * pb_y + g_0_yyyyy_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxzz_0[i] = 5.0 * g_0_yyyy_0_xxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxxzz_0[i] * pb_y +
                                   g_0_yyyyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyyy_0[i] = 5.0 * g_0_yyyy_0_xxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyyy_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyy_0[i] * pb_y + g_0_yyyyy_0_xxxyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyyz_0[i] = 5.0 * g_0_yyyy_0_xxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyyz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyz_0[i] * pb_y + g_0_yyyyy_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyzz_0[i] = 5.0 * g_0_yyyy_0_xxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyzz_1[i] * fti_ab_0 +
                                   g_0_yyyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyzz_0[i] * pb_y + g_0_yyyyy_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxzzz_0[i] = 5.0 * g_0_yyyy_0_xxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxzzz_0[i] * pb_y +
                                   g_0_yyyyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyyy_0[i] = 5.0 * g_0_yyyy_0_xxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyyy_1[i] * fti_ab_0 +
                                   4.0 * g_0_yyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyy_0[i] * pb_y + g_0_yyyyy_0_xxyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyyz_0[i] = 5.0 * g_0_yyyy_0_xxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyyz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyz_0[i] * pb_y + g_0_yyyyy_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyzz_0[i] = 5.0 * g_0_yyyy_0_xxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyzz_0[i] * pb_y + g_0_yyyyy_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyzzz_0[i] = 5.0 * g_0_yyyy_0_xxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyzzz_1[i] * fti_ab_0 +
                                   g_0_yyyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyzzz_0[i] * pb_y + g_0_yyyyy_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxzzzz_0[i] = 5.0 * g_0_yyyy_0_xxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxzzzz_0[i] * pb_y +
                                   g_0_yyyyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyyy_0[i] = 5.0 * g_0_yyyy_0_xyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyyy_1[i] * fti_ab_0 +
                                   5.0 * g_0_yyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyy_0[i] * pb_y + g_0_yyyyy_0_xyyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyyz_0[i] = 5.0 * g_0_yyyy_0_xyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyyz_1[i] * fti_ab_0 +
                                   4.0 * g_0_yyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyz_0[i] * pb_y + g_0_yyyyy_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyzz_0[i] = 5.0 * g_0_yyyy_0_xyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyzz_0[i] * pb_y + g_0_yyyyy_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyzzz_0[i] = 5.0 * g_0_yyyy_0_xyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyzzz_0[i] * pb_y + g_0_yyyyy_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyzzzz_0[i] = 5.0 * g_0_yyyy_0_xyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyzzzz_1[i] * fti_ab_0 +
                                   g_0_yyyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzzzz_0[i] * pb_y + g_0_yyyyy_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xzzzzz_0[i] = 5.0 * g_0_yyyy_0_xzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xzzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xzzzzz_0[i] * pb_y +
                                   g_0_yyyyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyyy_0[i] = 5.0 * g_0_yyyy_0_yyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyyy_1[i] * fti_ab_0 +
                                   6.0 * g_0_yyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyy_0[i] * pb_y + g_0_yyyyy_0_yyyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyyz_0[i] = 5.0 * g_0_yyyy_0_yyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyyz_1[i] * fti_ab_0 +
                                   5.0 * g_0_yyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyz_0[i] * pb_y + g_0_yyyyy_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyzz_0[i] = 5.0 * g_0_yyyy_0_yyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_yyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyzz_0[i] * pb_y + g_0_yyyyy_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyzzz_0[i] = 5.0 * g_0_yyyy_0_yyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyzzz_0[i] * pb_y + g_0_yyyyy_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyzzzz_0[i] = 5.0 * g_0_yyyy_0_yyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyzzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyzzzz_0[i] * pb_y + g_0_yyyyy_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yzzzzz_0[i] = 5.0 * g_0_yyyy_0_yzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yzzzzz_1[i] * fti_ab_0 +
                                   g_0_yyyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yzzzzz_0[i] * pb_y + g_0_yyyyy_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_zzzzzz_0[i] = 5.0 * g_0_yyyy_0_zzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_zzzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_zzzzzz_0[i] * pb_y +
                                   g_0_yyyyy_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 616-644 components of targeted buffer : SISI

    auto g_0_yyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 616);

    auto g_0_yyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 617);

    auto g_0_yyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 618);

    auto g_0_yyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 619);

    auto g_0_yyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 620);

    auto g_0_yyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 621);

    auto g_0_yyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 622);

    auto g_0_yyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 623);

    auto g_0_yyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 624);

    auto g_0_yyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 625);

    auto g_0_yyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 626);

    auto g_0_yyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 627);

    auto g_0_yyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 628);

    auto g_0_yyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 629);

    auto g_0_yyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 630);

    auto g_0_yyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 631);

    auto g_0_yyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 632);

    auto g_0_yyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 633);

    auto g_0_yyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 634);

    auto g_0_yyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 635);

    auto g_0_yyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 636);

    auto g_0_yyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 637);

    auto g_0_yyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 638);

    auto g_0_yyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 639);

    auto g_0_yyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 640);

    auto g_0_yyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 641);

    auto g_0_yyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 642);

    auto g_0_yyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 643);

#pragma omp simd aligned(g_0_yyyyy_0_xxxxx_1,       \
                             g_0_yyyyy_0_xxxxxx_0,  \
                             g_0_yyyyy_0_xxxxxx_1,  \
                             g_0_yyyyy_0_xxxxxy_0,  \
                             g_0_yyyyy_0_xxxxxy_1,  \
                             g_0_yyyyy_0_xxxxxz_0,  \
                             g_0_yyyyy_0_xxxxxz_1,  \
                             g_0_yyyyy_0_xxxxy_1,   \
                             g_0_yyyyy_0_xxxxyy_0,  \
                             g_0_yyyyy_0_xxxxyy_1,  \
                             g_0_yyyyy_0_xxxxyz_0,  \
                             g_0_yyyyy_0_xxxxyz_1,  \
                             g_0_yyyyy_0_xxxxz_1,   \
                             g_0_yyyyy_0_xxxxzz_0,  \
                             g_0_yyyyy_0_xxxxzz_1,  \
                             g_0_yyyyy_0_xxxyy_1,   \
                             g_0_yyyyy_0_xxxyyy_0,  \
                             g_0_yyyyy_0_xxxyyy_1,  \
                             g_0_yyyyy_0_xxxyyz_0,  \
                             g_0_yyyyy_0_xxxyyz_1,  \
                             g_0_yyyyy_0_xxxyz_1,   \
                             g_0_yyyyy_0_xxxyzz_0,  \
                             g_0_yyyyy_0_xxxyzz_1,  \
                             g_0_yyyyy_0_xxxzz_1,   \
                             g_0_yyyyy_0_xxxzzz_0,  \
                             g_0_yyyyy_0_xxxzzz_1,  \
                             g_0_yyyyy_0_xxyyy_1,   \
                             g_0_yyyyy_0_xxyyyy_0,  \
                             g_0_yyyyy_0_xxyyyy_1,  \
                             g_0_yyyyy_0_xxyyyz_0,  \
                             g_0_yyyyy_0_xxyyyz_1,  \
                             g_0_yyyyy_0_xxyyz_1,   \
                             g_0_yyyyy_0_xxyyzz_0,  \
                             g_0_yyyyy_0_xxyyzz_1,  \
                             g_0_yyyyy_0_xxyzz_1,   \
                             g_0_yyyyy_0_xxyzzz_0,  \
                             g_0_yyyyy_0_xxyzzz_1,  \
                             g_0_yyyyy_0_xxzzz_1,   \
                             g_0_yyyyy_0_xxzzzz_0,  \
                             g_0_yyyyy_0_xxzzzz_1,  \
                             g_0_yyyyy_0_xyyyy_1,   \
                             g_0_yyyyy_0_xyyyyy_0,  \
                             g_0_yyyyy_0_xyyyyy_1,  \
                             g_0_yyyyy_0_xyyyyz_0,  \
                             g_0_yyyyy_0_xyyyyz_1,  \
                             g_0_yyyyy_0_xyyyz_1,   \
                             g_0_yyyyy_0_xyyyzz_0,  \
                             g_0_yyyyy_0_xyyyzz_1,  \
                             g_0_yyyyy_0_xyyzz_1,   \
                             g_0_yyyyy_0_xyyzzz_0,  \
                             g_0_yyyyy_0_xyyzzz_1,  \
                             g_0_yyyyy_0_xyzzz_1,   \
                             g_0_yyyyy_0_xyzzzz_0,  \
                             g_0_yyyyy_0_xyzzzz_1,  \
                             g_0_yyyyy_0_xzzzz_1,   \
                             g_0_yyyyy_0_xzzzzz_0,  \
                             g_0_yyyyy_0_xzzzzz_1,  \
                             g_0_yyyyy_0_yyyyy_1,   \
                             g_0_yyyyy_0_yyyyyy_0,  \
                             g_0_yyyyy_0_yyyyyy_1,  \
                             g_0_yyyyy_0_yyyyyz_0,  \
                             g_0_yyyyy_0_yyyyyz_1,  \
                             g_0_yyyyy_0_yyyyz_1,   \
                             g_0_yyyyy_0_yyyyzz_0,  \
                             g_0_yyyyy_0_yyyyzz_1,  \
                             g_0_yyyyy_0_yyyzz_1,   \
                             g_0_yyyyy_0_yyyzzz_0,  \
                             g_0_yyyyy_0_yyyzzz_1,  \
                             g_0_yyyyy_0_yyzzz_1,   \
                             g_0_yyyyy_0_yyzzzz_0,  \
                             g_0_yyyyy_0_yyzzzz_1,  \
                             g_0_yyyyy_0_yzzzz_1,   \
                             g_0_yyyyy_0_yzzzzz_0,  \
                             g_0_yyyyy_0_yzzzzz_1,  \
                             g_0_yyyyy_0_zzzzz_1,   \
                             g_0_yyyyy_0_zzzzzz_0,  \
                             g_0_yyyyy_0_zzzzzz_1,  \
                             g_0_yyyyyz_0_xxxxxx_0, \
                             g_0_yyyyyz_0_xxxxxy_0, \
                             g_0_yyyyyz_0_xxxxxz_0, \
                             g_0_yyyyyz_0_xxxxyy_0, \
                             g_0_yyyyyz_0_xxxxyz_0, \
                             g_0_yyyyyz_0_xxxxzz_0, \
                             g_0_yyyyyz_0_xxxyyy_0, \
                             g_0_yyyyyz_0_xxxyyz_0, \
                             g_0_yyyyyz_0_xxxyzz_0, \
                             g_0_yyyyyz_0_xxxzzz_0, \
                             g_0_yyyyyz_0_xxyyyy_0, \
                             g_0_yyyyyz_0_xxyyyz_0, \
                             g_0_yyyyyz_0_xxyyzz_0, \
                             g_0_yyyyyz_0_xxyzzz_0, \
                             g_0_yyyyyz_0_xxzzzz_0, \
                             g_0_yyyyyz_0_xyyyyy_0, \
                             g_0_yyyyyz_0_xyyyyz_0, \
                             g_0_yyyyyz_0_xyyyzz_0, \
                             g_0_yyyyyz_0_xyyzzz_0, \
                             g_0_yyyyyz_0_xyzzzz_0, \
                             g_0_yyyyyz_0_xzzzzz_0, \
                             g_0_yyyyyz_0_yyyyyy_0, \
                             g_0_yyyyyz_0_yyyyyz_0, \
                             g_0_yyyyyz_0_yyyyzz_0, \
                             g_0_yyyyyz_0_yyyzzz_0, \
                             g_0_yyyyyz_0_yyzzzz_0, \
                             g_0_yyyyyz_0_yzzzzz_0, \
                             g_0_yyyyyz_0_zzzzzz_0, \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyz_0_xxxxxx_0[i] = g_0_yyyyy_0_xxxxxx_0[i] * pb_z + g_0_yyyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxy_0[i] = g_0_yyyyy_0_xxxxxy_0[i] * pb_z + g_0_yyyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxz_0[i] = g_0_yyyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxz_0[i] * pb_z + g_0_yyyyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxyy_0[i] = g_0_yyyyy_0_xxxxyy_0[i] * pb_z + g_0_yyyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxyz_0[i] = g_0_yyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyz_0[i] * pb_z + g_0_yyyyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxzz_0[i] = 2.0 * g_0_yyyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxzz_0[i] * pb_z + g_0_yyyyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyyy_0[i] = g_0_yyyyy_0_xxxyyy_0[i] * pb_z + g_0_yyyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyyz_0[i] = g_0_yyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyz_0[i] * pb_z + g_0_yyyyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyzz_0[i] = 2.0 * g_0_yyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyzz_0[i] * pb_z + g_0_yyyyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxzzz_0[i] = 3.0 * g_0_yyyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxzzz_0[i] * pb_z + g_0_yyyyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyyy_0[i] = g_0_yyyyy_0_xxyyyy_0[i] * pb_z + g_0_yyyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyyz_0[i] = g_0_yyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyz_0[i] * pb_z + g_0_yyyyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyzz_0[i] = 2.0 * g_0_yyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyzz_0[i] * pb_z + g_0_yyyyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyzzz_0[i] = 3.0 * g_0_yyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyzzz_0[i] * pb_z + g_0_yyyyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxzzzz_0[i] = 4.0 * g_0_yyyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxzzzz_0[i] * pb_z + g_0_yyyyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyyy_0[i] = g_0_yyyyy_0_xyyyyy_0[i] * pb_z + g_0_yyyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyyz_0[i] = g_0_yyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyz_0[i] * pb_z + g_0_yyyyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyzz_0[i] = 2.0 * g_0_yyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyzz_0[i] * pb_z + g_0_yyyyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyzzz_0[i] = 3.0 * g_0_yyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyzzz_0[i] * pb_z + g_0_yyyyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyzzzz_0[i] = 4.0 * g_0_yyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzzzz_0[i] * pb_z + g_0_yyyyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xzzzzz_0[i] = 5.0 * g_0_yyyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xzzzzz_0[i] * pb_z + g_0_yyyyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyyy_0[i] = g_0_yyyyy_0_yyyyyy_0[i] * pb_z + g_0_yyyyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyyz_0[i] = g_0_yyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyz_0[i] * pb_z + g_0_yyyyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyzz_0[i] = 2.0 * g_0_yyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyzz_0[i] * pb_z + g_0_yyyyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyzzz_0[i] = 3.0 * g_0_yyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyzzz_0[i] * pb_z + g_0_yyyyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyzzzz_0[i] = 4.0 * g_0_yyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyzzzz_0[i] * pb_z + g_0_yyyyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yzzzzz_0[i] = 5.0 * g_0_yyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yzzzzz_0[i] * pb_z + g_0_yyyyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_zzzzzz_0[i] = 6.0 * g_0_yyyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_zzzzzz_0[i] * pb_z + g_0_yyyyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 644-672 components of targeted buffer : SISI

    auto g_0_yyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 644);

    auto g_0_yyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 645);

    auto g_0_yyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 646);

    auto g_0_yyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 647);

    auto g_0_yyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 648);

    auto g_0_yyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 649);

    auto g_0_yyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 650);

    auto g_0_yyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 651);

    auto g_0_yyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 652);

    auto g_0_yyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 653);

    auto g_0_yyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 654);

    auto g_0_yyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 655);

    auto g_0_yyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 656);

    auto g_0_yyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 657);

    auto g_0_yyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 658);

    auto g_0_yyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 659);

    auto g_0_yyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 660);

    auto g_0_yyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 661);

    auto g_0_yyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 662);

    auto g_0_yyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 663);

    auto g_0_yyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 664);

    auto g_0_yyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 665);

    auto g_0_yyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 666);

    auto g_0_yyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 667);

    auto g_0_yyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 668);

    auto g_0_yyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 669);

    auto g_0_yyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 670);

    auto g_0_yyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 671);

#pragma omp simd aligned(g_0_yyyy_0_xxxxxy_0,       \
                             g_0_yyyy_0_xxxxxy_1,   \
                             g_0_yyyy_0_xxxxyy_0,   \
                             g_0_yyyy_0_xxxxyy_1,   \
                             g_0_yyyy_0_xxxyyy_0,   \
                             g_0_yyyy_0_xxxyyy_1,   \
                             g_0_yyyy_0_xxyyyy_0,   \
                             g_0_yyyy_0_xxyyyy_1,   \
                             g_0_yyyy_0_xyyyyy_0,   \
                             g_0_yyyy_0_xyyyyy_1,   \
                             g_0_yyyy_0_yyyyyy_0,   \
                             g_0_yyyy_0_yyyyyy_1,   \
                             g_0_yyyyz_0_xxxxxy_0,  \
                             g_0_yyyyz_0_xxxxxy_1,  \
                             g_0_yyyyz_0_xxxxyy_0,  \
                             g_0_yyyyz_0_xxxxyy_1,  \
                             g_0_yyyyz_0_xxxyyy_0,  \
                             g_0_yyyyz_0_xxxyyy_1,  \
                             g_0_yyyyz_0_xxyyyy_0,  \
                             g_0_yyyyz_0_xxyyyy_1,  \
                             g_0_yyyyz_0_xyyyyy_0,  \
                             g_0_yyyyz_0_xyyyyy_1,  \
                             g_0_yyyyz_0_yyyyyy_0,  \
                             g_0_yyyyz_0_yyyyyy_1,  \
                             g_0_yyyyzz_0_xxxxxx_0, \
                             g_0_yyyyzz_0_xxxxxy_0, \
                             g_0_yyyyzz_0_xxxxxz_0, \
                             g_0_yyyyzz_0_xxxxyy_0, \
                             g_0_yyyyzz_0_xxxxyz_0, \
                             g_0_yyyyzz_0_xxxxzz_0, \
                             g_0_yyyyzz_0_xxxyyy_0, \
                             g_0_yyyyzz_0_xxxyyz_0, \
                             g_0_yyyyzz_0_xxxyzz_0, \
                             g_0_yyyyzz_0_xxxzzz_0, \
                             g_0_yyyyzz_0_xxyyyy_0, \
                             g_0_yyyyzz_0_xxyyyz_0, \
                             g_0_yyyyzz_0_xxyyzz_0, \
                             g_0_yyyyzz_0_xxyzzz_0, \
                             g_0_yyyyzz_0_xxzzzz_0, \
                             g_0_yyyyzz_0_xyyyyy_0, \
                             g_0_yyyyzz_0_xyyyyz_0, \
                             g_0_yyyyzz_0_xyyyzz_0, \
                             g_0_yyyyzz_0_xyyzzz_0, \
                             g_0_yyyyzz_0_xyzzzz_0, \
                             g_0_yyyyzz_0_xzzzzz_0, \
                             g_0_yyyyzz_0_yyyyyy_0, \
                             g_0_yyyyzz_0_yyyyyz_0, \
                             g_0_yyyyzz_0_yyyyzz_0, \
                             g_0_yyyyzz_0_yyyzzz_0, \
                             g_0_yyyyzz_0_yyzzzz_0, \
                             g_0_yyyyzz_0_yzzzzz_0, \
                             g_0_yyyyzz_0_zzzzzz_0, \
                             g_0_yyyzz_0_xxxxxx_0,  \
                             g_0_yyyzz_0_xxxxxx_1,  \
                             g_0_yyyzz_0_xxxxxz_0,  \
                             g_0_yyyzz_0_xxxxxz_1,  \
                             g_0_yyyzz_0_xxxxyz_0,  \
                             g_0_yyyzz_0_xxxxyz_1,  \
                             g_0_yyyzz_0_xxxxz_1,   \
                             g_0_yyyzz_0_xxxxzz_0,  \
                             g_0_yyyzz_0_xxxxzz_1,  \
                             g_0_yyyzz_0_xxxyyz_0,  \
                             g_0_yyyzz_0_xxxyyz_1,  \
                             g_0_yyyzz_0_xxxyz_1,   \
                             g_0_yyyzz_0_xxxyzz_0,  \
                             g_0_yyyzz_0_xxxyzz_1,  \
                             g_0_yyyzz_0_xxxzz_1,   \
                             g_0_yyyzz_0_xxxzzz_0,  \
                             g_0_yyyzz_0_xxxzzz_1,  \
                             g_0_yyyzz_0_xxyyyz_0,  \
                             g_0_yyyzz_0_xxyyyz_1,  \
                             g_0_yyyzz_0_xxyyz_1,   \
                             g_0_yyyzz_0_xxyyzz_0,  \
                             g_0_yyyzz_0_xxyyzz_1,  \
                             g_0_yyyzz_0_xxyzz_1,   \
                             g_0_yyyzz_0_xxyzzz_0,  \
                             g_0_yyyzz_0_xxyzzz_1,  \
                             g_0_yyyzz_0_xxzzz_1,   \
                             g_0_yyyzz_0_xxzzzz_0,  \
                             g_0_yyyzz_0_xxzzzz_1,  \
                             g_0_yyyzz_0_xyyyyz_0,  \
                             g_0_yyyzz_0_xyyyyz_1,  \
                             g_0_yyyzz_0_xyyyz_1,   \
                             g_0_yyyzz_0_xyyyzz_0,  \
                             g_0_yyyzz_0_xyyyzz_1,  \
                             g_0_yyyzz_0_xyyzz_1,   \
                             g_0_yyyzz_0_xyyzzz_0,  \
                             g_0_yyyzz_0_xyyzzz_1,  \
                             g_0_yyyzz_0_xyzzz_1,   \
                             g_0_yyyzz_0_xyzzzz_0,  \
                             g_0_yyyzz_0_xyzzzz_1,  \
                             g_0_yyyzz_0_xzzzz_1,   \
                             g_0_yyyzz_0_xzzzzz_0,  \
                             g_0_yyyzz_0_xzzzzz_1,  \
                             g_0_yyyzz_0_yyyyyz_0,  \
                             g_0_yyyzz_0_yyyyyz_1,  \
                             g_0_yyyzz_0_yyyyz_1,   \
                             g_0_yyyzz_0_yyyyzz_0,  \
                             g_0_yyyzz_0_yyyyzz_1,  \
                             g_0_yyyzz_0_yyyzz_1,   \
                             g_0_yyyzz_0_yyyzzz_0,  \
                             g_0_yyyzz_0_yyyzzz_1,  \
                             g_0_yyyzz_0_yyzzz_1,   \
                             g_0_yyyzz_0_yyzzzz_0,  \
                             g_0_yyyzz_0_yyzzzz_1,  \
                             g_0_yyyzz_0_yzzzz_1,   \
                             g_0_yyyzz_0_yzzzzz_0,  \
                             g_0_yyyzz_0_yzzzzz_1,  \
                             g_0_yyyzz_0_zzzzz_1,   \
                             g_0_yyyzz_0_zzzzzz_0,  \
                             g_0_yyyzz_0_zzzzzz_1,  \
                             g_0_yyzz_0_xxxxxx_0,   \
                             g_0_yyzz_0_xxxxxx_1,   \
                             g_0_yyzz_0_xxxxxz_0,   \
                             g_0_yyzz_0_xxxxxz_1,   \
                             g_0_yyzz_0_xxxxyz_0,   \
                             g_0_yyzz_0_xxxxyz_1,   \
                             g_0_yyzz_0_xxxxzz_0,   \
                             g_0_yyzz_0_xxxxzz_1,   \
                             g_0_yyzz_0_xxxyyz_0,   \
                             g_0_yyzz_0_xxxyyz_1,   \
                             g_0_yyzz_0_xxxyzz_0,   \
                             g_0_yyzz_0_xxxyzz_1,   \
                             g_0_yyzz_0_xxxzzz_0,   \
                             g_0_yyzz_0_xxxzzz_1,   \
                             g_0_yyzz_0_xxyyyz_0,   \
                             g_0_yyzz_0_xxyyyz_1,   \
                             g_0_yyzz_0_xxyyzz_0,   \
                             g_0_yyzz_0_xxyyzz_1,   \
                             g_0_yyzz_0_xxyzzz_0,   \
                             g_0_yyzz_0_xxyzzz_1,   \
                             g_0_yyzz_0_xxzzzz_0,   \
                             g_0_yyzz_0_xxzzzz_1,   \
                             g_0_yyzz_0_xyyyyz_0,   \
                             g_0_yyzz_0_xyyyyz_1,   \
                             g_0_yyzz_0_xyyyzz_0,   \
                             g_0_yyzz_0_xyyyzz_1,   \
                             g_0_yyzz_0_xyyzzz_0,   \
                             g_0_yyzz_0_xyyzzz_1,   \
                             g_0_yyzz_0_xyzzzz_0,   \
                             g_0_yyzz_0_xyzzzz_1,   \
                             g_0_yyzz_0_xzzzzz_0,   \
                             g_0_yyzz_0_xzzzzz_1,   \
                             g_0_yyzz_0_yyyyyz_0,   \
                             g_0_yyzz_0_yyyyyz_1,   \
                             g_0_yyzz_0_yyyyzz_0,   \
                             g_0_yyzz_0_yyyyzz_1,   \
                             g_0_yyzz_0_yyyzzz_0,   \
                             g_0_yyzz_0_yyyzzz_1,   \
                             g_0_yyzz_0_yyzzzz_0,   \
                             g_0_yyzz_0_yyzzzz_1,   \
                             g_0_yyzz_0_yzzzzz_0,   \
                             g_0_yyzz_0_yzzzzz_1,   \
                             g_0_yyzz_0_zzzzzz_0,   \
                             g_0_yyzz_0_zzzzzz_1,   \
                             wp_y,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzz_0_xxxxxx_0[i] = 3.0 * g_0_yyzz_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxxx_0[i] * pb_y +
                                   g_0_yyyzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxxy_0[i] =
            g_0_yyyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxxxy_0[i] * pb_z + g_0_yyyyz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxxxz_0[i] = 3.0 * g_0_yyzz_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxxz_0[i] * pb_y +
                                   g_0_yyyzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxyy_0[i] =
            g_0_yyyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxxyy_0[i] * pb_z + g_0_yyyyz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxxyz_0[i] = 3.0 * g_0_yyzz_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxyz_1[i] * fti_ab_0 +
                                   g_0_yyyzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyz_0[i] * pb_y + g_0_yyyzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxzz_0[i] = 3.0 * g_0_yyzz_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxzz_0[i] * pb_y +
                                   g_0_yyyzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxyyy_0[i] =
            g_0_yyyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxyyy_0[i] * pb_z + g_0_yyyyz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxyyz_0[i] = 3.0 * g_0_yyzz_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyyz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyz_0[i] * pb_y + g_0_yyyzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxyzz_0[i] = 3.0 * g_0_yyzz_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyzz_1[i] * fti_ab_0 +
                                   g_0_yyyzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyzz_0[i] * pb_y + g_0_yyyzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxzzz_0[i] = 3.0 * g_0_yyzz_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxzzz_0[i] * pb_y +
                                   g_0_yyyzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyyyy_0[i] =
            g_0_yyyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxyyyy_0[i] * pb_z + g_0_yyyyz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxyyyz_0[i] = 3.0 * g_0_yyzz_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyyz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyz_0[i] * pb_y + g_0_yyyzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyyzz_0[i] = 3.0 * g_0_yyzz_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyzz_0[i] * pb_y + g_0_yyyzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyzzz_0[i] = 3.0 * g_0_yyzz_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyzzz_1[i] * fti_ab_0 +
                                   g_0_yyyzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyzzz_0[i] * pb_y + g_0_yyyzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxzzzz_0[i] = 3.0 * g_0_yyzz_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxzzzz_0[i] * pb_y +
                                   g_0_yyyzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyyyy_0[i] =
            g_0_yyyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xyyyyy_0[i] * pb_z + g_0_yyyyz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xyyyyz_0[i] = 3.0 * g_0_yyzz_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyyz_1[i] * fti_ab_0 +
                                   4.0 * g_0_yyyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyz_0[i] * pb_y + g_0_yyyzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyyzz_0[i] = 3.0 * g_0_yyzz_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyzz_0[i] * pb_y + g_0_yyyzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyzzz_0[i] = 3.0 * g_0_yyzz_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyzzz_0[i] * pb_y + g_0_yyyzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyzzzz_0[i] = 3.0 * g_0_yyzz_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyzzzz_1[i] * fti_ab_0 +
                                   g_0_yyyzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyzzzz_0[i] * pb_y + g_0_yyyzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xzzzzz_0[i] = 3.0 * g_0_yyzz_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xzzzzz_0[i] * pb_y +
                                   g_0_yyyzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyyyy_0[i] =
            g_0_yyyy_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_yyyyyy_0[i] * pb_z + g_0_yyyyz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_yyyyyz_0[i] = 3.0 * g_0_yyzz_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyyz_1[i] * fti_ab_0 +
                                   5.0 * g_0_yyyzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyyyz_0[i] * pb_y + g_0_yyyzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyyzz_0[i] = 3.0 * g_0_yyzz_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_yyyzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyyzz_0[i] * pb_y + g_0_yyyzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyzzz_0[i] = 3.0 * g_0_yyzz_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyyzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyzzz_0[i] * pb_y + g_0_yyyzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyzzzz_0[i] = 3.0 * g_0_yyzz_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyzzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyyzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyzzzz_0[i] * pb_y + g_0_yyyzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yzzzzz_0[i] = 3.0 * g_0_yyzz_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yzzzzz_1[i] * fti_ab_0 +
                                   g_0_yyyzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yzzzzz_0[i] * pb_y + g_0_yyyzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_zzzzzz_0[i] = 3.0 * g_0_yyzz_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_zzzzzz_0[i] * pb_y +
                                   g_0_yyyzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 672-700 components of targeted buffer : SISI

    auto g_0_yyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 672);

    auto g_0_yyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 673);

    auto g_0_yyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 674);

    auto g_0_yyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 675);

    auto g_0_yyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 676);

    auto g_0_yyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 677);

    auto g_0_yyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 678);

    auto g_0_yyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 679);

    auto g_0_yyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 680);

    auto g_0_yyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 681);

    auto g_0_yyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 682);

    auto g_0_yyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 683);

    auto g_0_yyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 684);

    auto g_0_yyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 685);

    auto g_0_yyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 686);

    auto g_0_yyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 687);

    auto g_0_yyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 688);

    auto g_0_yyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 689);

    auto g_0_yyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 690);

    auto g_0_yyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 691);

    auto g_0_yyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 692);

    auto g_0_yyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 693);

    auto g_0_yyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 694);

    auto g_0_yyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 695);

    auto g_0_yyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 696);

    auto g_0_yyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 697);

    auto g_0_yyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 698);

    auto g_0_yyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 699);

#pragma omp simd aligned(g_0_yyyz_0_xxxxxy_0,       \
                             g_0_yyyz_0_xxxxxy_1,   \
                             g_0_yyyz_0_xxxxyy_0,   \
                             g_0_yyyz_0_xxxxyy_1,   \
                             g_0_yyyz_0_xxxyyy_0,   \
                             g_0_yyyz_0_xxxyyy_1,   \
                             g_0_yyyz_0_xxyyyy_0,   \
                             g_0_yyyz_0_xxyyyy_1,   \
                             g_0_yyyz_0_xyyyyy_0,   \
                             g_0_yyyz_0_xyyyyy_1,   \
                             g_0_yyyz_0_yyyyyy_0,   \
                             g_0_yyyz_0_yyyyyy_1,   \
                             g_0_yyyzz_0_xxxxxy_0,  \
                             g_0_yyyzz_0_xxxxxy_1,  \
                             g_0_yyyzz_0_xxxxyy_0,  \
                             g_0_yyyzz_0_xxxxyy_1,  \
                             g_0_yyyzz_0_xxxyyy_0,  \
                             g_0_yyyzz_0_xxxyyy_1,  \
                             g_0_yyyzz_0_xxyyyy_0,  \
                             g_0_yyyzz_0_xxyyyy_1,  \
                             g_0_yyyzz_0_xyyyyy_0,  \
                             g_0_yyyzz_0_xyyyyy_1,  \
                             g_0_yyyzz_0_yyyyyy_0,  \
                             g_0_yyyzz_0_yyyyyy_1,  \
                             g_0_yyyzzz_0_xxxxxx_0, \
                             g_0_yyyzzz_0_xxxxxy_0, \
                             g_0_yyyzzz_0_xxxxxz_0, \
                             g_0_yyyzzz_0_xxxxyy_0, \
                             g_0_yyyzzz_0_xxxxyz_0, \
                             g_0_yyyzzz_0_xxxxzz_0, \
                             g_0_yyyzzz_0_xxxyyy_0, \
                             g_0_yyyzzz_0_xxxyyz_0, \
                             g_0_yyyzzz_0_xxxyzz_0, \
                             g_0_yyyzzz_0_xxxzzz_0, \
                             g_0_yyyzzz_0_xxyyyy_0, \
                             g_0_yyyzzz_0_xxyyyz_0, \
                             g_0_yyyzzz_0_xxyyzz_0, \
                             g_0_yyyzzz_0_xxyzzz_0, \
                             g_0_yyyzzz_0_xxzzzz_0, \
                             g_0_yyyzzz_0_xyyyyy_0, \
                             g_0_yyyzzz_0_xyyyyz_0, \
                             g_0_yyyzzz_0_xyyyzz_0, \
                             g_0_yyyzzz_0_xyyzzz_0, \
                             g_0_yyyzzz_0_xyzzzz_0, \
                             g_0_yyyzzz_0_xzzzzz_0, \
                             g_0_yyyzzz_0_yyyyyy_0, \
                             g_0_yyyzzz_0_yyyyyz_0, \
                             g_0_yyyzzz_0_yyyyzz_0, \
                             g_0_yyyzzz_0_yyyzzz_0, \
                             g_0_yyyzzz_0_yyzzzz_0, \
                             g_0_yyyzzz_0_yzzzzz_0, \
                             g_0_yyyzzz_0_zzzzzz_0, \
                             g_0_yyzzz_0_xxxxxx_0,  \
                             g_0_yyzzz_0_xxxxxx_1,  \
                             g_0_yyzzz_0_xxxxxz_0,  \
                             g_0_yyzzz_0_xxxxxz_1,  \
                             g_0_yyzzz_0_xxxxyz_0,  \
                             g_0_yyzzz_0_xxxxyz_1,  \
                             g_0_yyzzz_0_xxxxz_1,   \
                             g_0_yyzzz_0_xxxxzz_0,  \
                             g_0_yyzzz_0_xxxxzz_1,  \
                             g_0_yyzzz_0_xxxyyz_0,  \
                             g_0_yyzzz_0_xxxyyz_1,  \
                             g_0_yyzzz_0_xxxyz_1,   \
                             g_0_yyzzz_0_xxxyzz_0,  \
                             g_0_yyzzz_0_xxxyzz_1,  \
                             g_0_yyzzz_0_xxxzz_1,   \
                             g_0_yyzzz_0_xxxzzz_0,  \
                             g_0_yyzzz_0_xxxzzz_1,  \
                             g_0_yyzzz_0_xxyyyz_0,  \
                             g_0_yyzzz_0_xxyyyz_1,  \
                             g_0_yyzzz_0_xxyyz_1,   \
                             g_0_yyzzz_0_xxyyzz_0,  \
                             g_0_yyzzz_0_xxyyzz_1,  \
                             g_0_yyzzz_0_xxyzz_1,   \
                             g_0_yyzzz_0_xxyzzz_0,  \
                             g_0_yyzzz_0_xxyzzz_1,  \
                             g_0_yyzzz_0_xxzzz_1,   \
                             g_0_yyzzz_0_xxzzzz_0,  \
                             g_0_yyzzz_0_xxzzzz_1,  \
                             g_0_yyzzz_0_xyyyyz_0,  \
                             g_0_yyzzz_0_xyyyyz_1,  \
                             g_0_yyzzz_0_xyyyz_1,   \
                             g_0_yyzzz_0_xyyyzz_0,  \
                             g_0_yyzzz_0_xyyyzz_1,  \
                             g_0_yyzzz_0_xyyzz_1,   \
                             g_0_yyzzz_0_xyyzzz_0,  \
                             g_0_yyzzz_0_xyyzzz_1,  \
                             g_0_yyzzz_0_xyzzz_1,   \
                             g_0_yyzzz_0_xyzzzz_0,  \
                             g_0_yyzzz_0_xyzzzz_1,  \
                             g_0_yyzzz_0_xzzzz_1,   \
                             g_0_yyzzz_0_xzzzzz_0,  \
                             g_0_yyzzz_0_xzzzzz_1,  \
                             g_0_yyzzz_0_yyyyyz_0,  \
                             g_0_yyzzz_0_yyyyyz_1,  \
                             g_0_yyzzz_0_yyyyz_1,   \
                             g_0_yyzzz_0_yyyyzz_0,  \
                             g_0_yyzzz_0_yyyyzz_1,  \
                             g_0_yyzzz_0_yyyzz_1,   \
                             g_0_yyzzz_0_yyyzzz_0,  \
                             g_0_yyzzz_0_yyyzzz_1,  \
                             g_0_yyzzz_0_yyzzz_1,   \
                             g_0_yyzzz_0_yyzzzz_0,  \
                             g_0_yyzzz_0_yyzzzz_1,  \
                             g_0_yyzzz_0_yzzzz_1,   \
                             g_0_yyzzz_0_yzzzzz_0,  \
                             g_0_yyzzz_0_yzzzzz_1,  \
                             g_0_yyzzz_0_zzzzz_1,   \
                             g_0_yyzzz_0_zzzzzz_0,  \
                             g_0_yyzzz_0_zzzzzz_1,  \
                             g_0_yzzz_0_xxxxxx_0,   \
                             g_0_yzzz_0_xxxxxx_1,   \
                             g_0_yzzz_0_xxxxxz_0,   \
                             g_0_yzzz_0_xxxxxz_1,   \
                             g_0_yzzz_0_xxxxyz_0,   \
                             g_0_yzzz_0_xxxxyz_1,   \
                             g_0_yzzz_0_xxxxzz_0,   \
                             g_0_yzzz_0_xxxxzz_1,   \
                             g_0_yzzz_0_xxxyyz_0,   \
                             g_0_yzzz_0_xxxyyz_1,   \
                             g_0_yzzz_0_xxxyzz_0,   \
                             g_0_yzzz_0_xxxyzz_1,   \
                             g_0_yzzz_0_xxxzzz_0,   \
                             g_0_yzzz_0_xxxzzz_1,   \
                             g_0_yzzz_0_xxyyyz_0,   \
                             g_0_yzzz_0_xxyyyz_1,   \
                             g_0_yzzz_0_xxyyzz_0,   \
                             g_0_yzzz_0_xxyyzz_1,   \
                             g_0_yzzz_0_xxyzzz_0,   \
                             g_0_yzzz_0_xxyzzz_1,   \
                             g_0_yzzz_0_xxzzzz_0,   \
                             g_0_yzzz_0_xxzzzz_1,   \
                             g_0_yzzz_0_xyyyyz_0,   \
                             g_0_yzzz_0_xyyyyz_1,   \
                             g_0_yzzz_0_xyyyzz_0,   \
                             g_0_yzzz_0_xyyyzz_1,   \
                             g_0_yzzz_0_xyyzzz_0,   \
                             g_0_yzzz_0_xyyzzz_1,   \
                             g_0_yzzz_0_xyzzzz_0,   \
                             g_0_yzzz_0_xyzzzz_1,   \
                             g_0_yzzz_0_xzzzzz_0,   \
                             g_0_yzzz_0_xzzzzz_1,   \
                             g_0_yzzz_0_yyyyyz_0,   \
                             g_0_yzzz_0_yyyyyz_1,   \
                             g_0_yzzz_0_yyyyzz_0,   \
                             g_0_yzzz_0_yyyyzz_1,   \
                             g_0_yzzz_0_yyyzzz_0,   \
                             g_0_yzzz_0_yyyzzz_1,   \
                             g_0_yzzz_0_yyzzzz_0,   \
                             g_0_yzzz_0_yyzzzz_1,   \
                             g_0_yzzz_0_yzzzzz_0,   \
                             g_0_yzzz_0_yzzzzz_1,   \
                             g_0_yzzz_0_zzzzzz_0,   \
                             g_0_yzzz_0_zzzzzz_1,   \
                             wp_y,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzz_0_xxxxxx_0[i] = 2.0 * g_0_yzzz_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxxx_0[i] * pb_y +
                                   g_0_yyzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxxy_0[i] = 2.0 * g_0_yyyz_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxxxy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxxy_0[i] * pb_z +
                                   g_0_yyyzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxxxz_0[i] = 2.0 * g_0_yzzz_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxxz_0[i] * pb_y +
                                   g_0_yyzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxyy_0[i] = 2.0 * g_0_yyyz_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxxyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxyy_0[i] * pb_z +
                                   g_0_yyyzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxxyz_0[i] = 2.0 * g_0_yzzz_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                   g_0_yyzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyz_0[i] * pb_y + g_0_yyzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxzz_0[i] = 2.0 * g_0_yzzz_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxzz_0[i] * pb_y +
                                   g_0_yyzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxyyy_0[i] = 2.0 * g_0_yyyz_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxyyy_0[i] * pb_z +
                                   g_0_yyyzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxyyz_0[i] = 2.0 * g_0_yzzz_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyz_0[i] * pb_y + g_0_yyzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxyzz_0[i] = 2.0 * g_0_yzzz_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                   g_0_yyzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyzz_0[i] * pb_y + g_0_yyzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxzzz_0[i] = 2.0 * g_0_yzzz_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxzzz_0[i] * pb_y +
                                   g_0_yyzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyyyy_0[i] = 2.0 * g_0_yyyz_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxyyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxyyyy_0[i] * pb_z +
                                   g_0_yyyzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxyyyz_0[i] = 2.0 * g_0_yzzz_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyz_0[i] * pb_y + g_0_yyzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyyzz_0[i] = 2.0 * g_0_yzzz_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyzz_0[i] * pb_y + g_0_yyzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyzzz_0[i] = 2.0 * g_0_yzzz_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                   g_0_yyzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyzzz_0[i] * pb_y + g_0_yyzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxzzzz_0[i] = 2.0 * g_0_yzzz_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxzzzz_0[i] * pb_y +
                                   g_0_yyzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyyyy_0[i] = 2.0 * g_0_yyyz_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xyyyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xyyyyy_0[i] * pb_z +
                                   g_0_yyyzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xyyyyz_0[i] = 2.0 * g_0_yzzz_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                   4.0 * g_0_yyzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyz_0[i] * pb_y + g_0_yyzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyyzz_0[i] = 2.0 * g_0_yzzz_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyzz_0[i] * pb_y + g_0_yyzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyzzz_0[i] = 2.0 * g_0_yzzz_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyzzz_0[i] * pb_y + g_0_yyzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyzzzz_0[i] = 2.0 * g_0_yzzz_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                   g_0_yyzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyzzzz_0[i] * pb_y + g_0_yyzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xzzzzz_0[i] = 2.0 * g_0_yzzz_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xzzzzz_0[i] * pb_y +
                                   g_0_yyzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyyyy_0[i] = 2.0 * g_0_yyyz_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_yyyyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_yyyyyy_0[i] * pb_z +
                                   g_0_yyyzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_yyyyyz_0[i] = 2.0 * g_0_yzzz_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                   5.0 * g_0_yyzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyyyz_0[i] * pb_y + g_0_yyzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyyzz_0[i] = 2.0 * g_0_yzzz_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_yyzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyyzz_0[i] * pb_y + g_0_yyzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyzzz_0[i] = 2.0 * g_0_yzzz_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_yyzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyzzz_0[i] * pb_y + g_0_yyzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyzzzz_0[i] = 2.0 * g_0_yzzz_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_yyzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyzzzz_0[i] * pb_y + g_0_yyzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yzzzzz_0[i] = 2.0 * g_0_yzzz_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                   g_0_yyzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yzzzzz_0[i] * pb_y + g_0_yyzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_zzzzzz_0[i] = 2.0 * g_0_yzzz_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_zzzzzz_0[i] * pb_y +
                                   g_0_yyzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 700-728 components of targeted buffer : SISI

    auto g_0_yyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 700);

    auto g_0_yyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 701);

    auto g_0_yyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 702);

    auto g_0_yyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 703);

    auto g_0_yyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 704);

    auto g_0_yyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 705);

    auto g_0_yyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 706);

    auto g_0_yyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 707);

    auto g_0_yyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 708);

    auto g_0_yyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 709);

    auto g_0_yyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 710);

    auto g_0_yyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 711);

    auto g_0_yyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 712);

    auto g_0_yyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 713);

    auto g_0_yyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 714);

    auto g_0_yyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 715);

    auto g_0_yyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 716);

    auto g_0_yyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 717);

    auto g_0_yyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 718);

    auto g_0_yyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 719);

    auto g_0_yyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 720);

    auto g_0_yyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 721);

    auto g_0_yyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 722);

    auto g_0_yyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 723);

    auto g_0_yyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 724);

    auto g_0_yyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 725);

    auto g_0_yyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 726);

    auto g_0_yyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 727);

#pragma omp simd aligned(g_0_yyzz_0_xxxxxy_0,       \
                             g_0_yyzz_0_xxxxxy_1,   \
                             g_0_yyzz_0_xxxxyy_0,   \
                             g_0_yyzz_0_xxxxyy_1,   \
                             g_0_yyzz_0_xxxyyy_0,   \
                             g_0_yyzz_0_xxxyyy_1,   \
                             g_0_yyzz_0_xxyyyy_0,   \
                             g_0_yyzz_0_xxyyyy_1,   \
                             g_0_yyzz_0_xyyyyy_0,   \
                             g_0_yyzz_0_xyyyyy_1,   \
                             g_0_yyzz_0_yyyyyy_0,   \
                             g_0_yyzz_0_yyyyyy_1,   \
                             g_0_yyzzz_0_xxxxxy_0,  \
                             g_0_yyzzz_0_xxxxxy_1,  \
                             g_0_yyzzz_0_xxxxyy_0,  \
                             g_0_yyzzz_0_xxxxyy_1,  \
                             g_0_yyzzz_0_xxxyyy_0,  \
                             g_0_yyzzz_0_xxxyyy_1,  \
                             g_0_yyzzz_0_xxyyyy_0,  \
                             g_0_yyzzz_0_xxyyyy_1,  \
                             g_0_yyzzz_0_xyyyyy_0,  \
                             g_0_yyzzz_0_xyyyyy_1,  \
                             g_0_yyzzz_0_yyyyyy_0,  \
                             g_0_yyzzz_0_yyyyyy_1,  \
                             g_0_yyzzzz_0_xxxxxx_0, \
                             g_0_yyzzzz_0_xxxxxy_0, \
                             g_0_yyzzzz_0_xxxxxz_0, \
                             g_0_yyzzzz_0_xxxxyy_0, \
                             g_0_yyzzzz_0_xxxxyz_0, \
                             g_0_yyzzzz_0_xxxxzz_0, \
                             g_0_yyzzzz_0_xxxyyy_0, \
                             g_0_yyzzzz_0_xxxyyz_0, \
                             g_0_yyzzzz_0_xxxyzz_0, \
                             g_0_yyzzzz_0_xxxzzz_0, \
                             g_0_yyzzzz_0_xxyyyy_0, \
                             g_0_yyzzzz_0_xxyyyz_0, \
                             g_0_yyzzzz_0_xxyyzz_0, \
                             g_0_yyzzzz_0_xxyzzz_0, \
                             g_0_yyzzzz_0_xxzzzz_0, \
                             g_0_yyzzzz_0_xyyyyy_0, \
                             g_0_yyzzzz_0_xyyyyz_0, \
                             g_0_yyzzzz_0_xyyyzz_0, \
                             g_0_yyzzzz_0_xyyzzz_0, \
                             g_0_yyzzzz_0_xyzzzz_0, \
                             g_0_yyzzzz_0_xzzzzz_0, \
                             g_0_yyzzzz_0_yyyyyy_0, \
                             g_0_yyzzzz_0_yyyyyz_0, \
                             g_0_yyzzzz_0_yyyyzz_0, \
                             g_0_yyzzzz_0_yyyzzz_0, \
                             g_0_yyzzzz_0_yyzzzz_0, \
                             g_0_yyzzzz_0_yzzzzz_0, \
                             g_0_yyzzzz_0_zzzzzz_0, \
                             g_0_yzzzz_0_xxxxxx_0,  \
                             g_0_yzzzz_0_xxxxxx_1,  \
                             g_0_yzzzz_0_xxxxxz_0,  \
                             g_0_yzzzz_0_xxxxxz_1,  \
                             g_0_yzzzz_0_xxxxyz_0,  \
                             g_0_yzzzz_0_xxxxyz_1,  \
                             g_0_yzzzz_0_xxxxz_1,   \
                             g_0_yzzzz_0_xxxxzz_0,  \
                             g_0_yzzzz_0_xxxxzz_1,  \
                             g_0_yzzzz_0_xxxyyz_0,  \
                             g_0_yzzzz_0_xxxyyz_1,  \
                             g_0_yzzzz_0_xxxyz_1,   \
                             g_0_yzzzz_0_xxxyzz_0,  \
                             g_0_yzzzz_0_xxxyzz_1,  \
                             g_0_yzzzz_0_xxxzz_1,   \
                             g_0_yzzzz_0_xxxzzz_0,  \
                             g_0_yzzzz_0_xxxzzz_1,  \
                             g_0_yzzzz_0_xxyyyz_0,  \
                             g_0_yzzzz_0_xxyyyz_1,  \
                             g_0_yzzzz_0_xxyyz_1,   \
                             g_0_yzzzz_0_xxyyzz_0,  \
                             g_0_yzzzz_0_xxyyzz_1,  \
                             g_0_yzzzz_0_xxyzz_1,   \
                             g_0_yzzzz_0_xxyzzz_0,  \
                             g_0_yzzzz_0_xxyzzz_1,  \
                             g_0_yzzzz_0_xxzzz_1,   \
                             g_0_yzzzz_0_xxzzzz_0,  \
                             g_0_yzzzz_0_xxzzzz_1,  \
                             g_0_yzzzz_0_xyyyyz_0,  \
                             g_0_yzzzz_0_xyyyyz_1,  \
                             g_0_yzzzz_0_xyyyz_1,   \
                             g_0_yzzzz_0_xyyyzz_0,  \
                             g_0_yzzzz_0_xyyyzz_1,  \
                             g_0_yzzzz_0_xyyzz_1,   \
                             g_0_yzzzz_0_xyyzzz_0,  \
                             g_0_yzzzz_0_xyyzzz_1,  \
                             g_0_yzzzz_0_xyzzz_1,   \
                             g_0_yzzzz_0_xyzzzz_0,  \
                             g_0_yzzzz_0_xyzzzz_1,  \
                             g_0_yzzzz_0_xzzzz_1,   \
                             g_0_yzzzz_0_xzzzzz_0,  \
                             g_0_yzzzz_0_xzzzzz_1,  \
                             g_0_yzzzz_0_yyyyyz_0,  \
                             g_0_yzzzz_0_yyyyyz_1,  \
                             g_0_yzzzz_0_yyyyz_1,   \
                             g_0_yzzzz_0_yyyyzz_0,  \
                             g_0_yzzzz_0_yyyyzz_1,  \
                             g_0_yzzzz_0_yyyzz_1,   \
                             g_0_yzzzz_0_yyyzzz_0,  \
                             g_0_yzzzz_0_yyyzzz_1,  \
                             g_0_yzzzz_0_yyzzz_1,   \
                             g_0_yzzzz_0_yyzzzz_0,  \
                             g_0_yzzzz_0_yyzzzz_1,  \
                             g_0_yzzzz_0_yzzzz_1,   \
                             g_0_yzzzz_0_yzzzzz_0,  \
                             g_0_yzzzz_0_yzzzzz_1,  \
                             g_0_yzzzz_0_zzzzz_1,   \
                             g_0_yzzzz_0_zzzzzz_0,  \
                             g_0_yzzzz_0_zzzzzz_1,  \
                             g_0_zzzz_0_xxxxxx_0,   \
                             g_0_zzzz_0_xxxxxx_1,   \
                             g_0_zzzz_0_xxxxxz_0,   \
                             g_0_zzzz_0_xxxxxz_1,   \
                             g_0_zzzz_0_xxxxyz_0,   \
                             g_0_zzzz_0_xxxxyz_1,   \
                             g_0_zzzz_0_xxxxzz_0,   \
                             g_0_zzzz_0_xxxxzz_1,   \
                             g_0_zzzz_0_xxxyyz_0,   \
                             g_0_zzzz_0_xxxyyz_1,   \
                             g_0_zzzz_0_xxxyzz_0,   \
                             g_0_zzzz_0_xxxyzz_1,   \
                             g_0_zzzz_0_xxxzzz_0,   \
                             g_0_zzzz_0_xxxzzz_1,   \
                             g_0_zzzz_0_xxyyyz_0,   \
                             g_0_zzzz_0_xxyyyz_1,   \
                             g_0_zzzz_0_xxyyzz_0,   \
                             g_0_zzzz_0_xxyyzz_1,   \
                             g_0_zzzz_0_xxyzzz_0,   \
                             g_0_zzzz_0_xxyzzz_1,   \
                             g_0_zzzz_0_xxzzzz_0,   \
                             g_0_zzzz_0_xxzzzz_1,   \
                             g_0_zzzz_0_xyyyyz_0,   \
                             g_0_zzzz_0_xyyyyz_1,   \
                             g_0_zzzz_0_xyyyzz_0,   \
                             g_0_zzzz_0_xyyyzz_1,   \
                             g_0_zzzz_0_xyyzzz_0,   \
                             g_0_zzzz_0_xyyzzz_1,   \
                             g_0_zzzz_0_xyzzzz_0,   \
                             g_0_zzzz_0_xyzzzz_1,   \
                             g_0_zzzz_0_xzzzzz_0,   \
                             g_0_zzzz_0_xzzzzz_1,   \
                             g_0_zzzz_0_yyyyyz_0,   \
                             g_0_zzzz_0_yyyyyz_1,   \
                             g_0_zzzz_0_yyyyzz_0,   \
                             g_0_zzzz_0_yyyyzz_1,   \
                             g_0_zzzz_0_yyyzzz_0,   \
                             g_0_zzzz_0_yyyzzz_1,   \
                             g_0_zzzz_0_yyzzzz_0,   \
                             g_0_zzzz_0_yyzzzz_1,   \
                             g_0_zzzz_0_yzzzzz_0,   \
                             g_0_zzzz_0_yzzzzz_1,   \
                             g_0_zzzz_0_zzzzzz_0,   \
                             g_0_zzzz_0_zzzzzz_1,   \
                             wp_y,                  \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzz_0_xxxxxx_0[i] =
            g_0_zzzz_0_xxxxxx_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxxx_0[i] * pb_y + g_0_yzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxxy_0[i] = 3.0 * g_0_yyzz_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxxy_0[i] * pb_z +
                                   g_0_yyzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxxxz_0[i] =
            g_0_zzzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxxz_0[i] * pb_y + g_0_yzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxyy_0[i] = 3.0 * g_0_yyzz_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxyy_0[i] * pb_z +
                                   g_0_yyzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxxyz_0[i] = g_0_zzzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxyz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_xxxxyz_0[i] * pb_y + g_0_yzzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxzz_0[i] =
            g_0_zzzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxzz_0[i] * pb_y + g_0_yzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxyyy_0[i] = 3.0 * g_0_yyzz_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxyyy_0[i] * pb_z +
                                   g_0_yyzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxyyz_0[i] = g_0_zzzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_xxxyz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_xxxyyz_0[i] * pb_y + g_0_yzzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxyzz_0[i] = g_0_zzzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxzz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_xxxyzz_0[i] * pb_y + g_0_yzzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxzzz_0[i] =
            g_0_zzzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxzzz_0[i] * pb_y + g_0_yzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyyyy_0[i] = 3.0 * g_0_yyzz_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxyyyy_0[i] * pb_z +
                                   g_0_yyzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxyyyz_0[i] = g_0_zzzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzz_0_xxyyz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_xxyyyz_0[i] * pb_y + g_0_yzzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyyzz_0[i] = g_0_zzzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_xxyzz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_xxyyzz_0[i] * pb_y + g_0_yzzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyzzz_0[i] = g_0_zzzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxzzz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_xxyzzz_0[i] * pb_y + g_0_yzzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxzzzz_0[i] =
            g_0_zzzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxzzzz_0[i] * pb_y + g_0_yzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyyyy_0[i] = 3.0 * g_0_yyzz_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xyyyyy_0[i] * pb_z +
                                   g_0_yyzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xyyyyz_0[i] = g_0_zzzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzz_0_xyyyz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_xyyyyz_0[i] * pb_y + g_0_yzzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyyzz_0[i] = g_0_zzzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzz_0_xyyzz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_xyyyzz_0[i] * pb_y + g_0_yzzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyzzz_0[i] = g_0_zzzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_xyzzz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_xyyzzz_0[i] * pb_y + g_0_yzzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyzzzz_0[i] = g_0_zzzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_xyzzzz_0[i] * pb_y + g_0_yzzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xzzzzz_0[i] =
            g_0_zzzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xzzzzz_0[i] * pb_y + g_0_yzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyyyy_0[i] = 3.0 * g_0_yyzz_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_yyyyyy_0[i] * pb_z +
                                   g_0_yyzzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_yyyyyz_0[i] = g_0_zzzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yzzzz_0_yyyyz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_yyyyyz_0[i] * pb_y + g_0_yzzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyyzz_0[i] = g_0_zzzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzz_0_yyyzz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_yyyyzz_0[i] * pb_y + g_0_yzzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyzzz_0[i] = g_0_zzzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzz_0_yyzzz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_yyyzzz_0[i] * pb_y + g_0_yzzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyzzzz_0[i] = g_0_zzzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_yzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_yyzzzz_0[i] * pb_y + g_0_yzzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yzzzzz_0[i] = g_0_zzzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_zzzzz_1[i] * fi_abcd_0 +
                                   g_0_yzzzz_0_yzzzzz_0[i] * pb_y + g_0_yzzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_zzzzzz_0[i] =
            g_0_zzzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_zzzzzz_0[i] * pb_y + g_0_yzzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 728-756 components of targeted buffer : SISI

    auto g_0_yzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 728);

    auto g_0_yzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 729);

    auto g_0_yzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 730);

    auto g_0_yzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 731);

    auto g_0_yzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 732);

    auto g_0_yzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 733);

    auto g_0_yzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 734);

    auto g_0_yzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 735);

    auto g_0_yzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 736);

    auto g_0_yzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 737);

    auto g_0_yzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 738);

    auto g_0_yzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 739);

    auto g_0_yzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 740);

    auto g_0_yzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 741);

    auto g_0_yzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 742);

    auto g_0_yzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 743);

    auto g_0_yzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 744);

    auto g_0_yzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 745);

    auto g_0_yzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 746);

    auto g_0_yzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 747);

    auto g_0_yzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 748);

    auto g_0_yzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 749);

    auto g_0_yzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 750);

    auto g_0_yzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 751);

    auto g_0_yzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 752);

    auto g_0_yzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 753);

    auto g_0_yzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 754);

    auto g_0_yzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 755);

#pragma omp simd aligned(g_0_yzzzzz_0_xxxxxx_0,     \
                             g_0_yzzzzz_0_xxxxxy_0, \
                             g_0_yzzzzz_0_xxxxxz_0, \
                             g_0_yzzzzz_0_xxxxyy_0, \
                             g_0_yzzzzz_0_xxxxyz_0, \
                             g_0_yzzzzz_0_xxxxzz_0, \
                             g_0_yzzzzz_0_xxxyyy_0, \
                             g_0_yzzzzz_0_xxxyyz_0, \
                             g_0_yzzzzz_0_xxxyzz_0, \
                             g_0_yzzzzz_0_xxxzzz_0, \
                             g_0_yzzzzz_0_xxyyyy_0, \
                             g_0_yzzzzz_0_xxyyyz_0, \
                             g_0_yzzzzz_0_xxyyzz_0, \
                             g_0_yzzzzz_0_xxyzzz_0, \
                             g_0_yzzzzz_0_xxzzzz_0, \
                             g_0_yzzzzz_0_xyyyyy_0, \
                             g_0_yzzzzz_0_xyyyyz_0, \
                             g_0_yzzzzz_0_xyyyzz_0, \
                             g_0_yzzzzz_0_xyyzzz_0, \
                             g_0_yzzzzz_0_xyzzzz_0, \
                             g_0_yzzzzz_0_xzzzzz_0, \
                             g_0_yzzzzz_0_yyyyyy_0, \
                             g_0_yzzzzz_0_yyyyyz_0, \
                             g_0_yzzzzz_0_yyyyzz_0, \
                             g_0_yzzzzz_0_yyyzzz_0, \
                             g_0_yzzzzz_0_yyzzzz_0, \
                             g_0_yzzzzz_0_yzzzzz_0, \
                             g_0_yzzzzz_0_zzzzzz_0, \
                             g_0_zzzzz_0_xxxxx_1,   \
                             g_0_zzzzz_0_xxxxxx_0,  \
                             g_0_zzzzz_0_xxxxxx_1,  \
                             g_0_zzzzz_0_xxxxxy_0,  \
                             g_0_zzzzz_0_xxxxxy_1,  \
                             g_0_zzzzz_0_xxxxxz_0,  \
                             g_0_zzzzz_0_xxxxxz_1,  \
                             g_0_zzzzz_0_xxxxy_1,   \
                             g_0_zzzzz_0_xxxxyy_0,  \
                             g_0_zzzzz_0_xxxxyy_1,  \
                             g_0_zzzzz_0_xxxxyz_0,  \
                             g_0_zzzzz_0_xxxxyz_1,  \
                             g_0_zzzzz_0_xxxxz_1,   \
                             g_0_zzzzz_0_xxxxzz_0,  \
                             g_0_zzzzz_0_xxxxzz_1,  \
                             g_0_zzzzz_0_xxxyy_1,   \
                             g_0_zzzzz_0_xxxyyy_0,  \
                             g_0_zzzzz_0_xxxyyy_1,  \
                             g_0_zzzzz_0_xxxyyz_0,  \
                             g_0_zzzzz_0_xxxyyz_1,  \
                             g_0_zzzzz_0_xxxyz_1,   \
                             g_0_zzzzz_0_xxxyzz_0,  \
                             g_0_zzzzz_0_xxxyzz_1,  \
                             g_0_zzzzz_0_xxxzz_1,   \
                             g_0_zzzzz_0_xxxzzz_0,  \
                             g_0_zzzzz_0_xxxzzz_1,  \
                             g_0_zzzzz_0_xxyyy_1,   \
                             g_0_zzzzz_0_xxyyyy_0,  \
                             g_0_zzzzz_0_xxyyyy_1,  \
                             g_0_zzzzz_0_xxyyyz_0,  \
                             g_0_zzzzz_0_xxyyyz_1,  \
                             g_0_zzzzz_0_xxyyz_1,   \
                             g_0_zzzzz_0_xxyyzz_0,  \
                             g_0_zzzzz_0_xxyyzz_1,  \
                             g_0_zzzzz_0_xxyzz_1,   \
                             g_0_zzzzz_0_xxyzzz_0,  \
                             g_0_zzzzz_0_xxyzzz_1,  \
                             g_0_zzzzz_0_xxzzz_1,   \
                             g_0_zzzzz_0_xxzzzz_0,  \
                             g_0_zzzzz_0_xxzzzz_1,  \
                             g_0_zzzzz_0_xyyyy_1,   \
                             g_0_zzzzz_0_xyyyyy_0,  \
                             g_0_zzzzz_0_xyyyyy_1,  \
                             g_0_zzzzz_0_xyyyyz_0,  \
                             g_0_zzzzz_0_xyyyyz_1,  \
                             g_0_zzzzz_0_xyyyz_1,   \
                             g_0_zzzzz_0_xyyyzz_0,  \
                             g_0_zzzzz_0_xyyyzz_1,  \
                             g_0_zzzzz_0_xyyzz_1,   \
                             g_0_zzzzz_0_xyyzzz_0,  \
                             g_0_zzzzz_0_xyyzzz_1,  \
                             g_0_zzzzz_0_xyzzz_1,   \
                             g_0_zzzzz_0_xyzzzz_0,  \
                             g_0_zzzzz_0_xyzzzz_1,  \
                             g_0_zzzzz_0_xzzzz_1,   \
                             g_0_zzzzz_0_xzzzzz_0,  \
                             g_0_zzzzz_0_xzzzzz_1,  \
                             g_0_zzzzz_0_yyyyy_1,   \
                             g_0_zzzzz_0_yyyyyy_0,  \
                             g_0_zzzzz_0_yyyyyy_1,  \
                             g_0_zzzzz_0_yyyyyz_0,  \
                             g_0_zzzzz_0_yyyyyz_1,  \
                             g_0_zzzzz_0_yyyyz_1,   \
                             g_0_zzzzz_0_yyyyzz_0,  \
                             g_0_zzzzz_0_yyyyzz_1,  \
                             g_0_zzzzz_0_yyyzz_1,   \
                             g_0_zzzzz_0_yyyzzz_0,  \
                             g_0_zzzzz_0_yyyzzz_1,  \
                             g_0_zzzzz_0_yyzzz_1,   \
                             g_0_zzzzz_0_yyzzzz_0,  \
                             g_0_zzzzz_0_yyzzzz_1,  \
                             g_0_zzzzz_0_yzzzz_1,   \
                             g_0_zzzzz_0_yzzzzz_0,  \
                             g_0_zzzzz_0_yzzzzz_1,  \
                             g_0_zzzzz_0_zzzzz_1,   \
                             g_0_zzzzz_0_zzzzzz_0,  \
                             g_0_zzzzz_0_zzzzzz_1,  \
                             wp_y,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzz_0_xxxxxx_0[i] = g_0_zzzzz_0_xxxxxx_0[i] * pb_y + g_0_zzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxy_0[i] = g_0_zzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxy_0[i] * pb_y + g_0_zzzzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxz_0[i] = g_0_zzzzz_0_xxxxxz_0[i] * pb_y + g_0_zzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxyy_0[i] = 2.0 * g_0_zzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyy_0[i] * pb_y + g_0_zzzzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxyz_0[i] = g_0_zzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyz_0[i] * pb_y + g_0_zzzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxzz_0[i] = g_0_zzzzz_0_xxxxzz_0[i] * pb_y + g_0_zzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyyy_0[i] = 3.0 * g_0_zzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyy_0[i] * pb_y + g_0_zzzzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyyz_0[i] = 2.0 * g_0_zzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyz_0[i] * pb_y + g_0_zzzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyzz_0[i] = g_0_zzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyzz_0[i] * pb_y + g_0_zzzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxzzz_0[i] = g_0_zzzzz_0_xxxzzz_0[i] * pb_y + g_0_zzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyyy_0[i] = 4.0 * g_0_zzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyy_0[i] * pb_y + g_0_zzzzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyyz_0[i] = 3.0 * g_0_zzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyz_0[i] * pb_y + g_0_zzzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyzz_0[i] = 2.0 * g_0_zzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyzz_0[i] * pb_y + g_0_zzzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyzzz_0[i] = g_0_zzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyzzz_0[i] * pb_y + g_0_zzzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxzzzz_0[i] = g_0_zzzzz_0_xxzzzz_0[i] * pb_y + g_0_zzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyyy_0[i] = 5.0 * g_0_zzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyy_0[i] * pb_y + g_0_zzzzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyyz_0[i] = 4.0 * g_0_zzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyz_0[i] * pb_y + g_0_zzzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyzz_0[i] = 3.0 * g_0_zzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyzz_0[i] * pb_y + g_0_zzzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyzzz_0[i] = 2.0 * g_0_zzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyzzz_0[i] * pb_y + g_0_zzzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyzzzz_0[i] = g_0_zzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzzzz_0[i] * pb_y + g_0_zzzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xzzzzz_0[i] = g_0_zzzzz_0_xzzzzz_0[i] * pb_y + g_0_zzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyyy_0[i] = 6.0 * g_0_zzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyy_0[i] * pb_y + g_0_zzzzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyyz_0[i] = 5.0 * g_0_zzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyz_0[i] * pb_y + g_0_zzzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyzz_0[i] = 4.0 * g_0_zzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyzz_0[i] * pb_y + g_0_zzzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyzzz_0[i] = 3.0 * g_0_zzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyzzz_0[i] * pb_y + g_0_zzzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyzzzz_0[i] = 2.0 * g_0_zzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyzzzz_0[i] * pb_y + g_0_zzzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yzzzzz_0[i] = g_0_zzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yzzzzz_0[i] * pb_y + g_0_zzzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_zzzzzz_0[i] = g_0_zzzzz_0_zzzzzz_0[i] * pb_y + g_0_zzzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 756-784 components of targeted buffer : SISI

    auto g_0_zzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 756);

    auto g_0_zzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 757);

    auto g_0_zzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 758);

    auto g_0_zzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 759);

    auto g_0_zzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 760);

    auto g_0_zzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 761);

    auto g_0_zzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 762);

    auto g_0_zzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 763);

    auto g_0_zzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 764);

    auto g_0_zzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 765);

    auto g_0_zzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 766);

    auto g_0_zzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 767);

    auto g_0_zzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 768);

    auto g_0_zzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 769);

    auto g_0_zzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 770);

    auto g_0_zzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 771);

    auto g_0_zzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 772);

    auto g_0_zzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 773);

    auto g_0_zzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 774);

    auto g_0_zzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 775);

    auto g_0_zzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 776);

    auto g_0_zzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 777);

    auto g_0_zzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 778);

    auto g_0_zzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 779);

    auto g_0_zzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 780);

    auto g_0_zzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 781);

    auto g_0_zzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 782);

    auto g_0_zzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 783);

#pragma omp simd aligned(g_0_zzzz_0_xxxxxx_0,       \
                             g_0_zzzz_0_xxxxxx_1,   \
                             g_0_zzzz_0_xxxxxy_0,   \
                             g_0_zzzz_0_xxxxxy_1,   \
                             g_0_zzzz_0_xxxxxz_0,   \
                             g_0_zzzz_0_xxxxxz_1,   \
                             g_0_zzzz_0_xxxxyy_0,   \
                             g_0_zzzz_0_xxxxyy_1,   \
                             g_0_zzzz_0_xxxxyz_0,   \
                             g_0_zzzz_0_xxxxyz_1,   \
                             g_0_zzzz_0_xxxxzz_0,   \
                             g_0_zzzz_0_xxxxzz_1,   \
                             g_0_zzzz_0_xxxyyy_0,   \
                             g_0_zzzz_0_xxxyyy_1,   \
                             g_0_zzzz_0_xxxyyz_0,   \
                             g_0_zzzz_0_xxxyyz_1,   \
                             g_0_zzzz_0_xxxyzz_0,   \
                             g_0_zzzz_0_xxxyzz_1,   \
                             g_0_zzzz_0_xxxzzz_0,   \
                             g_0_zzzz_0_xxxzzz_1,   \
                             g_0_zzzz_0_xxyyyy_0,   \
                             g_0_zzzz_0_xxyyyy_1,   \
                             g_0_zzzz_0_xxyyyz_0,   \
                             g_0_zzzz_0_xxyyyz_1,   \
                             g_0_zzzz_0_xxyyzz_0,   \
                             g_0_zzzz_0_xxyyzz_1,   \
                             g_0_zzzz_0_xxyzzz_0,   \
                             g_0_zzzz_0_xxyzzz_1,   \
                             g_0_zzzz_0_xxzzzz_0,   \
                             g_0_zzzz_0_xxzzzz_1,   \
                             g_0_zzzz_0_xyyyyy_0,   \
                             g_0_zzzz_0_xyyyyy_1,   \
                             g_0_zzzz_0_xyyyyz_0,   \
                             g_0_zzzz_0_xyyyyz_1,   \
                             g_0_zzzz_0_xyyyzz_0,   \
                             g_0_zzzz_0_xyyyzz_1,   \
                             g_0_zzzz_0_xyyzzz_0,   \
                             g_0_zzzz_0_xyyzzz_1,   \
                             g_0_zzzz_0_xyzzzz_0,   \
                             g_0_zzzz_0_xyzzzz_1,   \
                             g_0_zzzz_0_xzzzzz_0,   \
                             g_0_zzzz_0_xzzzzz_1,   \
                             g_0_zzzz_0_yyyyyy_0,   \
                             g_0_zzzz_0_yyyyyy_1,   \
                             g_0_zzzz_0_yyyyyz_0,   \
                             g_0_zzzz_0_yyyyyz_1,   \
                             g_0_zzzz_0_yyyyzz_0,   \
                             g_0_zzzz_0_yyyyzz_1,   \
                             g_0_zzzz_0_yyyzzz_0,   \
                             g_0_zzzz_0_yyyzzz_1,   \
                             g_0_zzzz_0_yyzzzz_0,   \
                             g_0_zzzz_0_yyzzzz_1,   \
                             g_0_zzzz_0_yzzzzz_0,   \
                             g_0_zzzz_0_yzzzzz_1,   \
                             g_0_zzzz_0_zzzzzz_0,   \
                             g_0_zzzz_0_zzzzzz_1,   \
                             g_0_zzzzz_0_xxxxx_1,   \
                             g_0_zzzzz_0_xxxxxx_0,  \
                             g_0_zzzzz_0_xxxxxx_1,  \
                             g_0_zzzzz_0_xxxxxy_0,  \
                             g_0_zzzzz_0_xxxxxy_1,  \
                             g_0_zzzzz_0_xxxxxz_0,  \
                             g_0_zzzzz_0_xxxxxz_1,  \
                             g_0_zzzzz_0_xxxxy_1,   \
                             g_0_zzzzz_0_xxxxyy_0,  \
                             g_0_zzzzz_0_xxxxyy_1,  \
                             g_0_zzzzz_0_xxxxyz_0,  \
                             g_0_zzzzz_0_xxxxyz_1,  \
                             g_0_zzzzz_0_xxxxz_1,   \
                             g_0_zzzzz_0_xxxxzz_0,  \
                             g_0_zzzzz_0_xxxxzz_1,  \
                             g_0_zzzzz_0_xxxyy_1,   \
                             g_0_zzzzz_0_xxxyyy_0,  \
                             g_0_zzzzz_0_xxxyyy_1,  \
                             g_0_zzzzz_0_xxxyyz_0,  \
                             g_0_zzzzz_0_xxxyyz_1,  \
                             g_0_zzzzz_0_xxxyz_1,   \
                             g_0_zzzzz_0_xxxyzz_0,  \
                             g_0_zzzzz_0_xxxyzz_1,  \
                             g_0_zzzzz_0_xxxzz_1,   \
                             g_0_zzzzz_0_xxxzzz_0,  \
                             g_0_zzzzz_0_xxxzzz_1,  \
                             g_0_zzzzz_0_xxyyy_1,   \
                             g_0_zzzzz_0_xxyyyy_0,  \
                             g_0_zzzzz_0_xxyyyy_1,  \
                             g_0_zzzzz_0_xxyyyz_0,  \
                             g_0_zzzzz_0_xxyyyz_1,  \
                             g_0_zzzzz_0_xxyyz_1,   \
                             g_0_zzzzz_0_xxyyzz_0,  \
                             g_0_zzzzz_0_xxyyzz_1,  \
                             g_0_zzzzz_0_xxyzz_1,   \
                             g_0_zzzzz_0_xxyzzz_0,  \
                             g_0_zzzzz_0_xxyzzz_1,  \
                             g_0_zzzzz_0_xxzzz_1,   \
                             g_0_zzzzz_0_xxzzzz_0,  \
                             g_0_zzzzz_0_xxzzzz_1,  \
                             g_0_zzzzz_0_xyyyy_1,   \
                             g_0_zzzzz_0_xyyyyy_0,  \
                             g_0_zzzzz_0_xyyyyy_1,  \
                             g_0_zzzzz_0_xyyyyz_0,  \
                             g_0_zzzzz_0_xyyyyz_1,  \
                             g_0_zzzzz_0_xyyyz_1,   \
                             g_0_zzzzz_0_xyyyzz_0,  \
                             g_0_zzzzz_0_xyyyzz_1,  \
                             g_0_zzzzz_0_xyyzz_1,   \
                             g_0_zzzzz_0_xyyzzz_0,  \
                             g_0_zzzzz_0_xyyzzz_1,  \
                             g_0_zzzzz_0_xyzzz_1,   \
                             g_0_zzzzz_0_xyzzzz_0,  \
                             g_0_zzzzz_0_xyzzzz_1,  \
                             g_0_zzzzz_0_xzzzz_1,   \
                             g_0_zzzzz_0_xzzzzz_0,  \
                             g_0_zzzzz_0_xzzzzz_1,  \
                             g_0_zzzzz_0_yyyyy_1,   \
                             g_0_zzzzz_0_yyyyyy_0,  \
                             g_0_zzzzz_0_yyyyyy_1,  \
                             g_0_zzzzz_0_yyyyyz_0,  \
                             g_0_zzzzz_0_yyyyyz_1,  \
                             g_0_zzzzz_0_yyyyz_1,   \
                             g_0_zzzzz_0_yyyyzz_0,  \
                             g_0_zzzzz_0_yyyyzz_1,  \
                             g_0_zzzzz_0_yyyzz_1,   \
                             g_0_zzzzz_0_yyyzzz_0,  \
                             g_0_zzzzz_0_yyyzzz_1,  \
                             g_0_zzzzz_0_yyzzz_1,   \
                             g_0_zzzzz_0_yyzzzz_0,  \
                             g_0_zzzzz_0_yyzzzz_1,  \
                             g_0_zzzzz_0_yzzzz_1,   \
                             g_0_zzzzz_0_yzzzzz_0,  \
                             g_0_zzzzz_0_yzzzzz_1,  \
                             g_0_zzzzz_0_zzzzz_1,   \
                             g_0_zzzzz_0_zzzzzz_0,  \
                             g_0_zzzzz_0_zzzzzz_1,  \
                             g_0_zzzzzz_0_xxxxxx_0, \
                             g_0_zzzzzz_0_xxxxxy_0, \
                             g_0_zzzzzz_0_xxxxxz_0, \
                             g_0_zzzzzz_0_xxxxyy_0, \
                             g_0_zzzzzz_0_xxxxyz_0, \
                             g_0_zzzzzz_0_xxxxzz_0, \
                             g_0_zzzzzz_0_xxxyyy_0, \
                             g_0_zzzzzz_0_xxxyyz_0, \
                             g_0_zzzzzz_0_xxxyzz_0, \
                             g_0_zzzzzz_0_xxxzzz_0, \
                             g_0_zzzzzz_0_xxyyyy_0, \
                             g_0_zzzzzz_0_xxyyyz_0, \
                             g_0_zzzzzz_0_xxyyzz_0, \
                             g_0_zzzzzz_0_xxyzzz_0, \
                             g_0_zzzzzz_0_xxzzzz_0, \
                             g_0_zzzzzz_0_xyyyyy_0, \
                             g_0_zzzzzz_0_xyyyyz_0, \
                             g_0_zzzzzz_0_xyyyzz_0, \
                             g_0_zzzzzz_0_xyyzzz_0, \
                             g_0_zzzzzz_0_xyzzzz_0, \
                             g_0_zzzzzz_0_xzzzzz_0, \
                             g_0_zzzzzz_0_yyyyyy_0, \
                             g_0_zzzzzz_0_yyyyyz_0, \
                             g_0_zzzzzz_0_yyyyzz_0, \
                             g_0_zzzzzz_0_yyyzzz_0, \
                             g_0_zzzzzz_0_yyzzzz_0, \
                             g_0_zzzzzz_0_yzzzzz_0, \
                             g_0_zzzzzz_0_zzzzzz_0, \
                             wp_z,                  \
                             c_exps,                \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzz_0_xxxxxx_0[i] = 5.0 * g_0_zzzz_0_xxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxxxx_0[i] * pb_z +
                                   g_0_zzzzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxy_0[i] = 5.0 * g_0_zzzz_0_xxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxxxy_0[i] * pb_z +
                                   g_0_zzzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxz_0[i] = 5.0 * g_0_zzzz_0_xxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                   g_0_zzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxz_0[i] * pb_z + g_0_zzzzz_0_xxxxxz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxyy_0[i] = 5.0 * g_0_zzzz_0_xxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxxyy_0[i] * pb_z +
                                   g_0_zzzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxyz_0[i] = 5.0 * g_0_zzzz_0_xxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                   g_0_zzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyz_0[i] * pb_z + g_0_zzzzz_0_xxxxyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxzz_0[i] = 5.0 * g_0_zzzz_0_xxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_zzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxzz_0[i] * pb_z + g_0_zzzzz_0_xxxxzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyyy_0[i] = 5.0 * g_0_zzzz_0_xxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxyyy_0[i] * pb_z +
                                   g_0_zzzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyyz_0[i] = 5.0 * g_0_zzzz_0_xxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                   g_0_zzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyz_0[i] * pb_z + g_0_zzzzz_0_xxxyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyzz_0[i] = 5.0 * g_0_zzzz_0_xxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_zzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyzz_0[i] * pb_z + g_0_zzzzz_0_xxxyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxzzz_0[i] = 5.0 * g_0_zzzz_0_xxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_zzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxzzz_0[i] * pb_z + g_0_zzzzz_0_xxxzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyyy_0[i] = 5.0 * g_0_zzzz_0_xxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxyyyy_0[i] * pb_z +
                                   g_0_zzzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyyz_0[i] = 5.0 * g_0_zzzz_0_xxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                   g_0_zzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyz_0[i] * pb_z + g_0_zzzzz_0_xxyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyzz_0[i] = 5.0 * g_0_zzzz_0_xxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_zzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyzz_0[i] * pb_z + g_0_zzzzz_0_xxyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyzzz_0[i] = 5.0 * g_0_zzzz_0_xxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_zzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyzzz_0[i] * pb_z + g_0_zzzzz_0_xxyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxzzzz_0[i] = 5.0 * g_0_zzzz_0_xxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_zzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxzzzz_0[i] * pb_z + g_0_zzzzz_0_xxzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyyy_0[i] = 5.0 * g_0_zzzz_0_xyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xyyyyy_0[i] * pb_z +
                                   g_0_zzzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyyz_0[i] = 5.0 * g_0_zzzz_0_xyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                   g_0_zzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyz_0[i] * pb_z + g_0_zzzzz_0_xyyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyzz_0[i] = 5.0 * g_0_zzzz_0_xyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_zzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyzz_0[i] * pb_z + g_0_zzzzz_0_xyyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyzzz_0[i] = 5.0 * g_0_zzzz_0_xyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_zzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyzzz_0[i] * pb_z + g_0_zzzzz_0_xyyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyzzzz_0[i] = 5.0 * g_0_zzzz_0_xyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_zzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzzzz_0[i] * pb_z + g_0_zzzzz_0_xyzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xzzzzz_0[i] = 5.0 * g_0_zzzz_0_xzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                   5.0 * g_0_zzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xzzzzz_0[i] * pb_z + g_0_zzzzz_0_xzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyyy_0[i] = 5.0 * g_0_zzzz_0_yyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_yyyyyy_0[i] * pb_z +
                                   g_0_zzzzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyyz_0[i] = 5.0 * g_0_zzzz_0_yyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                   g_0_zzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyz_0[i] * pb_z + g_0_zzzzz_0_yyyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyzz_0[i] = 5.0 * g_0_zzzz_0_yyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                   2.0 * g_0_zzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyzz_0[i] * pb_z + g_0_zzzzz_0_yyyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyzzz_0[i] = 5.0 * g_0_zzzz_0_yyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                   3.0 * g_0_zzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyzzz_0[i] * pb_z + g_0_zzzzz_0_yyyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyzzzz_0[i] = 5.0 * g_0_zzzz_0_yyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                   4.0 * g_0_zzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyzzzz_0[i] * pb_z + g_0_zzzzz_0_yyzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yzzzzz_0[i] = 5.0 * g_0_zzzz_0_yzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                   5.0 * g_0_zzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yzzzzz_0[i] * pb_z + g_0_zzzzz_0_yzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_zzzzzz_0[i] = 5.0 * g_0_zzzz_0_zzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                   6.0 * g_0_zzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_zzzzzz_0[i] * pb_z + g_0_zzzzz_0_zzzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
