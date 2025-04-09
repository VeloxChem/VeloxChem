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

#include "ThreeCenterElectronRepulsionPrimRecKSI.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ksi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksi,
                                 size_t idx_eri_0_hsi,
                                 size_t idx_eri_1_hsi,
                                 size_t idx_eri_1_ish,
                                 size_t idx_eri_1_isi,
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

    /// Set up components of auxilary buffer : HSI

    auto g_xxxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi);

    auto g_xxxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 1);

    auto g_xxxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 2);

    auto g_xxxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 3);

    auto g_xxxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 4);

    auto g_xxxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 5);

    auto g_xxxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 6);

    auto g_xxxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 7);

    auto g_xxxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 8);

    auto g_xxxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 9);

    auto g_xxxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 10);

    auto g_xxxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 11);

    auto g_xxxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 12);

    auto g_xxxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 13);

    auto g_xxxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 14);

    auto g_xxxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 15);

    auto g_xxxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 16);

    auto g_xxxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 17);

    auto g_xxxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 18);

    auto g_xxxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 19);

    auto g_xxxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 20);

    auto g_xxxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 21);

    auto g_xxxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 22);

    auto g_xxxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 23);

    auto g_xxxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 24);

    auto g_xxxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 25);

    auto g_xxxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 26);

    auto g_xxxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 27);

    auto g_xxxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 28);

    auto g_xxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 30);

    auto g_xxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 33);

    auto g_xxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 37);

    auto g_xxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 42);

    auto g_xxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 48);

    auto g_xxxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 56);

    auto g_xxxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 57);

    auto g_xxxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 59);

    auto g_xxxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 62);

    auto g_xxxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 66);

    auto g_xxxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 71);

    auto g_xxxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 84);

    auto g_xxxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 85);

    auto g_xxxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 86);

    auto g_xxxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 87);

    auto g_xxxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 88);

    auto g_xxxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 89);

    auto g_xxxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 90);

    auto g_xxxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 91);

    auto g_xxxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 92);

    auto g_xxxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 93);

    auto g_xxxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 94);

    auto g_xxxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 95);

    auto g_xxxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 96);

    auto g_xxxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 97);

    auto g_xxxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 98);

    auto g_xxxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 99);

    auto g_xxxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 100);

    auto g_xxxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 101);

    auto g_xxxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 102);

    auto g_xxxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 103);

    auto g_xxxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 104);

    auto g_xxxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 105);

    auto g_xxxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 106);

    auto g_xxxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 107);

    auto g_xxxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 108);

    auto g_xxxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 109);

    auto g_xxxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 110);

    auto g_xxxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 111);

    auto g_xxxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 140);

    auto g_xxxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 141);

    auto g_xxxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 142);

    auto g_xxxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 143);

    auto g_xxxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 144);

    auto g_xxxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 145);

    auto g_xxxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 146);

    auto g_xxxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 147);

    auto g_xxxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 148);

    auto g_xxxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 149);

    auto g_xxxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 150);

    auto g_xxxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 151);

    auto g_xxxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 152);

    auto g_xxxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 153);

    auto g_xxxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 154);

    auto g_xxxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 155);

    auto g_xxxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 156);

    auto g_xxxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 157);

    auto g_xxxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 158);

    auto g_xxxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 159);

    auto g_xxxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 160);

    auto g_xxxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 161);

    auto g_xxxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 162);

    auto g_xxxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 163);

    auto g_xxxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 164);

    auto g_xxxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 165);

    auto g_xxxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 166);

    auto g_xxxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 167);

    auto g_xxyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 168);

    auto g_xxyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 169);

    auto g_xxyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 170);

    auto g_xxyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 171);

    auto g_xxyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 172);

    auto g_xxyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 173);

    auto g_xxyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 174);

    auto g_xxyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 175);

    auto g_xxyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 176);

    auto g_xxyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 177);

    auto g_xxyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 178);

    auto g_xxyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 179);

    auto g_xxyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 180);

    auto g_xxyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 181);

    auto g_xxyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 182);

    auto g_xxyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 183);

    auto g_xxyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 184);

    auto g_xxyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 185);

    auto g_xxyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 186);

    auto g_xxyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 187);

    auto g_xxyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 188);

    auto g_xxyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 189);

    auto g_xxyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 190);

    auto g_xxyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 191);

    auto g_xxyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 192);

    auto g_xxyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 193);

    auto g_xxyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 194);

    auto g_xxyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 195);

    auto g_xxyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 197);

    auto g_xxyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 199);

    auto g_xxyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 202);

    auto g_xxyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 206);

    auto g_xxyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 211);

    auto g_xxyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 224);

    auto g_xxyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 226);

    auto g_xxyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 229);

    auto g_xxyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 233);

    auto g_xxyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 238);

    auto g_xxyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 244);

    auto g_xxzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 252);

    auto g_xxzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 253);

    auto g_xxzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 254);

    auto g_xxzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 255);

    auto g_xxzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 256);

    auto g_xxzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 257);

    auto g_xxzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 258);

    auto g_xxzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 259);

    auto g_xxzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 260);

    auto g_xxzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 261);

    auto g_xxzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 262);

    auto g_xxzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 263);

    auto g_xxzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 264);

    auto g_xxzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 265);

    auto g_xxzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 266);

    auto g_xxzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 267);

    auto g_xxzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 268);

    auto g_xxzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 269);

    auto g_xxzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 270);

    auto g_xxzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 271);

    auto g_xxzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 272);

    auto g_xxzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 273);

    auto g_xxzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 274);

    auto g_xxzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 275);

    auto g_xxzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 276);

    auto g_xxzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 277);

    auto g_xxzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 278);

    auto g_xxzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 279);

    auto g_xyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 281);

    auto g_xyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 283);

    auto g_xyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 284);

    auto g_xyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 286);

    auto g_xyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 287);

    auto g_xyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 288);

    auto g_xyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 290);

    auto g_xyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 291);

    auto g_xyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 292);

    auto g_xyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 293);

    auto g_xyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 295);

    auto g_xyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 296);

    auto g_xyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 297);

    auto g_xyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 298);

    auto g_xyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 299);

    auto g_xyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 301);

    auto g_xyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 302);

    auto g_xyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 303);

    auto g_xyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 304);

    auto g_xyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 305);

    auto g_xyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 306);

    auto g_xyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 307);

    auto g_xyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 340);

    auto g_xyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 343);

    auto g_xyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 344);

    auto g_xyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 347);

    auto g_xyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 348);

    auto g_xyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 349);

    auto g_xyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 352);

    auto g_xyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 353);

    auto g_xyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 354);

    auto g_xyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 355);

    auto g_xyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 357);

    auto g_xyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 358);

    auto g_xyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 359);

    auto g_xyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 360);

    auto g_xyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 361);

    auto g_xyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 362);

    auto g_xyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 363);

    auto g_xzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 394);

    auto g_xzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 396);

    auto g_xzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 397);

    auto g_xzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 399);

    auto g_xzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 400);

    auto g_xzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 401);

    auto g_xzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 403);

    auto g_xzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 404);

    auto g_xzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 405);

    auto g_xzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 406);

    auto g_xzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 408);

    auto g_xzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 409);

    auto g_xzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 410);

    auto g_xzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 411);

    auto g_xzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 412);

    auto g_xzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 413);

    auto g_xzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 414);

    auto g_xzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 415);

    auto g_xzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 416);

    auto g_xzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 417);

    auto g_xzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 418);

    auto g_xzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 419);

    auto g_yyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 420);

    auto g_yyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 421);

    auto g_yyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 422);

    auto g_yyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 423);

    auto g_yyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 424);

    auto g_yyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 425);

    auto g_yyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 426);

    auto g_yyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 427);

    auto g_yyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 428);

    auto g_yyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 429);

    auto g_yyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 430);

    auto g_yyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 431);

    auto g_yyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 432);

    auto g_yyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 433);

    auto g_yyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 434);

    auto g_yyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 435);

    auto g_yyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 436);

    auto g_yyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 437);

    auto g_yyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 438);

    auto g_yyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 439);

    auto g_yyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 440);

    auto g_yyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 441);

    auto g_yyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 442);

    auto g_yyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 443);

    auto g_yyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 444);

    auto g_yyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 445);

    auto g_yyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 446);

    auto g_yyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 447);

    auto g_yyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 449);

    auto g_yyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 451);

    auto g_yyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 454);

    auto g_yyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 458);

    auto g_yyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 463);

    auto g_yyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 469);

    auto g_yyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 476);

    auto g_yyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 477);

    auto g_yyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 478);

    auto g_yyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 479);

    auto g_yyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 480);

    auto g_yyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 481);

    auto g_yyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 482);

    auto g_yyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 483);

    auto g_yyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 484);

    auto g_yyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 485);

    auto g_yyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 486);

    auto g_yyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 487);

    auto g_yyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 488);

    auto g_yyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 489);

    auto g_yyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 490);

    auto g_yyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 491);

    auto g_yyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 492);

    auto g_yyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 493);

    auto g_yyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 494);

    auto g_yyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 495);

    auto g_yyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 496);

    auto g_yyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 497);

    auto g_yyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 498);

    auto g_yyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 499);

    auto g_yyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 500);

    auto g_yyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 501);

    auto g_yyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 502);

    auto g_yyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 503);

    auto g_yyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 504);

    auto g_yyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 505);

    auto g_yyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 506);

    auto g_yyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 507);

    auto g_yyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 508);

    auto g_yyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 509);

    auto g_yyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 510);

    auto g_yyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 511);

    auto g_yyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 512);

    auto g_yyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 513);

    auto g_yyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 514);

    auto g_yyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 515);

    auto g_yyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 516);

    auto g_yyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 517);

    auto g_yyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 518);

    auto g_yyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 519);

    auto g_yyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 520);

    auto g_yyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 521);

    auto g_yyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 522);

    auto g_yyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 523);

    auto g_yyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 524);

    auto g_yyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 525);

    auto g_yyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 526);

    auto g_yyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 527);

    auto g_yyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 528);

    auto g_yyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 529);

    auto g_yyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 530);

    auto g_yyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 531);

    auto g_yzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 532);

    auto g_yzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 534);

    auto g_yzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 536);

    auto g_yzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 537);

    auto g_yzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 539);

    auto g_yzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 540);

    auto g_yzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 541);

    auto g_yzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 543);

    auto g_yzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 544);

    auto g_yzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 545);

    auto g_yzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 546);

    auto g_yzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 548);

    auto g_yzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 549);

    auto g_yzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 550);

    auto g_yzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 551);

    auto g_yzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 552);

    auto g_yzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 554);

    auto g_yzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 555);

    auto g_yzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 556);

    auto g_yzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 557);

    auto g_yzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 558);

    auto g_yzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 559);

    auto g_zzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_hsi + 560);

    auto g_zzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_hsi + 561);

    auto g_zzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_hsi + 562);

    auto g_zzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_hsi + 563);

    auto g_zzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_hsi + 564);

    auto g_zzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_hsi + 565);

    auto g_zzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_hsi + 566);

    auto g_zzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_hsi + 567);

    auto g_zzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_hsi + 568);

    auto g_zzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_hsi + 569);

    auto g_zzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_hsi + 570);

    auto g_zzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_hsi + 571);

    auto g_zzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_hsi + 572);

    auto g_zzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_hsi + 573);

    auto g_zzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_hsi + 574);

    auto g_zzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 575);

    auto g_zzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 576);

    auto g_zzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 577);

    auto g_zzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 578);

    auto g_zzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 579);

    auto g_zzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 580);

    auto g_zzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_hsi + 581);

    auto g_zzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_hsi + 582);

    auto g_zzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_hsi + 583);

    auto g_zzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_hsi + 584);

    auto g_zzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_hsi + 585);

    auto g_zzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 586);

    auto g_zzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_hsi + 587);

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

    auto g_xxxxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 30);

    auto g_xxxxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 33);

    auto g_xxxxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 37);

    auto g_xxxxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 42);

    auto g_xxxxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 48);

    auto g_xxxxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_hsi + 56);

    auto g_xxxxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_hsi + 57);

    auto g_xxxxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 59);

    auto g_xxxxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 62);

    auto g_xxxxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 66);

    auto g_xxxxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 71);

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

    auto g_yyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_hsi + 451);

    auto g_yyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_hsi + 454);

    auto g_yyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_hsi + 458);

    auto g_yyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 463);

    auto g_yyyyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_hsi + 469);

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

    auto g_yzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 534);

    auto g_yzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 536);

    auto g_yzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 537);

    auto g_yzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 539);

    auto g_yzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 540);

    auto g_yzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 541);

    auto g_yzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 543);

    auto g_yzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 544);

    auto g_yzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 545);

    auto g_yzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 546);

    auto g_yzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 548);

    auto g_yzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 549);

    auto g_yzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 550);

    auto g_yzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 551);

    auto g_yzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 552);

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

    auto g_xxxxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 44);

    auto g_xxxxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 46);

    auto g_xxxxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 47);

    auto g_xxxxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 49);

    auto g_xxxxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 50);

    auto g_xxxxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 51);

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

    auto g_xxyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 256);

    auto g_xxyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 259);

    auto g_xxyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 260);

    auto g_xxyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 263);

    auto g_xxyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 264);

    auto g_xxyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 265);

    auto g_xxyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 268);

    auto g_xxyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 269);

    auto g_xxyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 270);

    auto g_xxyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 271);

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

    auto g_xyyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 361);

    auto g_xyyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 364);

    auto g_xyyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 365);

    auto g_xyyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 368);

    auto g_xyyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 369);

    auto g_xyyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 370);

    auto g_xyyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 373);

    auto g_xyyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 374);

    auto g_xyyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 375);

    auto g_xyyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 376);

    auto g_xyyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 382);

    auto g_xyyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 385);

    auto g_xyyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 386);

    auto g_xyyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 389);

    auto g_xyyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 390);

    auto g_xyyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 391);

    auto g_xyyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_ish + 394);

    auto g_xyyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_ish + 395);

    auto g_xyyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_ish + 396);

    auto g_xyyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_ish + 397);

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

    auto g_yyyyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_ish + 464);

    auto g_yyyyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_ish + 466);

    auto g_yyyyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_ish + 467);

    auto g_yyyyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_ish + 469);

    auto g_yyyyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_ish + 470);

    auto g_yyyyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_ish + 471);

    auto g_yyyyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_ish + 473);

    auto g_yyyyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_ish + 474);

    auto g_yyyyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_ish + 475);

    auto g_yyyyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_ish + 476);

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

    /// Set up components of auxilary buffer : ISI

    auto g_xxxxxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi);

    auto g_xxxxxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 1);

    auto g_xxxxxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 2);

    auto g_xxxxxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 3);

    auto g_xxxxxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 4);

    auto g_xxxxxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 5);

    auto g_xxxxxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 6);

    auto g_xxxxxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 7);

    auto g_xxxxxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 8);

    auto g_xxxxxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 9);

    auto g_xxxxxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 10);

    auto g_xxxxxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 11);

    auto g_xxxxxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 12);

    auto g_xxxxxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 13);

    auto g_xxxxxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 14);

    auto g_xxxxxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 15);

    auto g_xxxxxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 16);

    auto g_xxxxxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 17);

    auto g_xxxxxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 18);

    auto g_xxxxxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 19);

    auto g_xxxxxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 20);

    auto g_xxxxxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 21);

    auto g_xxxxxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 22);

    auto g_xxxxxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 23);

    auto g_xxxxxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 24);

    auto g_xxxxxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 25);

    auto g_xxxxxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 26);

    auto g_xxxxxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 27);

    auto g_xxxxxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 28);

    auto g_xxxxxy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 29);

    auto g_xxxxxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 30);

    auto g_xxxxxy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 31);

    auto g_xxxxxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 33);

    auto g_xxxxxy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 34);

    auto g_xxxxxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 37);

    auto g_xxxxxy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 38);

    auto g_xxxxxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 42);

    auto g_xxxxxy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 43);

    auto g_xxxxxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 48);

    auto g_xxxxxy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 49);

    auto g_xxxxxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 56);

    auto g_xxxxxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 57);

    auto g_xxxxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 58);

    auto g_xxxxxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 59);

    auto g_xxxxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 60);

    auto g_xxxxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 61);

    auto g_xxxxxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 62);

    auto g_xxxxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 63);

    auto g_xxxxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 64);

    auto g_xxxxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 65);

    auto g_xxxxxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 66);

    auto g_xxxxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 67);

    auto g_xxxxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 68);

    auto g_xxxxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 69);

    auto g_xxxxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 70);

    auto g_xxxxxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 71);

    auto g_xxxxxz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 72);

    auto g_xxxxxz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 73);

    auto g_xxxxxz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 74);

    auto g_xxxxxz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 75);

    auto g_xxxxxz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 76);

    auto g_xxxxxz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 78);

    auto g_xxxxxz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 79);

    auto g_xxxxxz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 80);

    auto g_xxxxxz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 81);

    auto g_xxxxxz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 82);

    auto g_xxxxxz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 83);

    auto g_xxxxyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 84);

    auto g_xxxxyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 85);

    auto g_xxxxyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 86);

    auto g_xxxxyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 87);

    auto g_xxxxyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 88);

    auto g_xxxxyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 89);

    auto g_xxxxyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 90);

    auto g_xxxxyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 91);

    auto g_xxxxyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 92);

    auto g_xxxxyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 93);

    auto g_xxxxyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 94);

    auto g_xxxxyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 95);

    auto g_xxxxyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 96);

    auto g_xxxxyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 97);

    auto g_xxxxyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 98);

    auto g_xxxxyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 99);

    auto g_xxxxyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 100);

    auto g_xxxxyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 101);

    auto g_xxxxyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 102);

    auto g_xxxxyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 103);

    auto g_xxxxyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 104);

    auto g_xxxxyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 105);

    auto g_xxxxyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 106);

    auto g_xxxxyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 107);

    auto g_xxxxyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 108);

    auto g_xxxxyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 109);

    auto g_xxxxyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 110);

    auto g_xxxxyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 111);

    auto g_xxxxzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 140);

    auto g_xxxxzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 141);

    auto g_xxxxzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 142);

    auto g_xxxxzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 143);

    auto g_xxxxzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 144);

    auto g_xxxxzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 145);

    auto g_xxxxzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 146);

    auto g_xxxxzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 147);

    auto g_xxxxzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 148);

    auto g_xxxxzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 149);

    auto g_xxxxzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 150);

    auto g_xxxxzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 151);

    auto g_xxxxzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 152);

    auto g_xxxxzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 153);

    auto g_xxxxzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 154);

    auto g_xxxxzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 155);

    auto g_xxxxzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 156);

    auto g_xxxxzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 157);

    auto g_xxxxzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 158);

    auto g_xxxxzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 159);

    auto g_xxxxzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 160);

    auto g_xxxxzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 161);

    auto g_xxxxzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 162);

    auto g_xxxxzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 163);

    auto g_xxxxzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 164);

    auto g_xxxxzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 165);

    auto g_xxxxzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 166);

    auto g_xxxxzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 167);

    auto g_xxxyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 168);

    auto g_xxxyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 169);

    auto g_xxxyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 170);

    auto g_xxxyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 171);

    auto g_xxxyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 172);

    auto g_xxxyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 173);

    auto g_xxxyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 174);

    auto g_xxxyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 175);

    auto g_xxxyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 176);

    auto g_xxxyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 177);

    auto g_xxxyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 178);

    auto g_xxxyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 179);

    auto g_xxxyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 180);

    auto g_xxxyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 181);

    auto g_xxxyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 182);

    auto g_xxxyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 183);

    auto g_xxxyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 184);

    auto g_xxxyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 185);

    auto g_xxxyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 186);

    auto g_xxxyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 187);

    auto g_xxxyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 188);

    auto g_xxxyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 189);

    auto g_xxxyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 190);

    auto g_xxxyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 191);

    auto g_xxxyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 192);

    auto g_xxxyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 193);

    auto g_xxxyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 194);

    auto g_xxxyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 195);

    auto g_xxxyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 197);

    auto g_xxxyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 199);

    auto g_xxxyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 202);

    auto g_xxxyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 206);

    auto g_xxxyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 211);

    auto g_xxxyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 224);

    auto g_xxxyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 226);

    auto g_xxxyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 229);

    auto g_xxxyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 233);

    auto g_xxxyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 238);

    auto g_xxxyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 244);

    auto g_xxxzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 252);

    auto g_xxxzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 253);

    auto g_xxxzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 254);

    auto g_xxxzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 255);

    auto g_xxxzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 256);

    auto g_xxxzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 257);

    auto g_xxxzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 258);

    auto g_xxxzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 259);

    auto g_xxxzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 260);

    auto g_xxxzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 261);

    auto g_xxxzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 262);

    auto g_xxxzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 263);

    auto g_xxxzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 264);

    auto g_xxxzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 265);

    auto g_xxxzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 266);

    auto g_xxxzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 267);

    auto g_xxxzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 268);

    auto g_xxxzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 269);

    auto g_xxxzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 270);

    auto g_xxxzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 271);

    auto g_xxxzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 272);

    auto g_xxxzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 273);

    auto g_xxxzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 274);

    auto g_xxxzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 275);

    auto g_xxxzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 276);

    auto g_xxxzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 277);

    auto g_xxxzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 278);

    auto g_xxxzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 279);

    auto g_xxyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 280);

    auto g_xxyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 281);

    auto g_xxyyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 282);

    auto g_xxyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 283);

    auto g_xxyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 284);

    auto g_xxyyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 285);

    auto g_xxyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 286);

    auto g_xxyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 287);

    auto g_xxyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 288);

    auto g_xxyyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 289);

    auto g_xxyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 290);

    auto g_xxyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 291);

    auto g_xxyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 292);

    auto g_xxyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 293);

    auto g_xxyyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 294);

    auto g_xxyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 295);

    auto g_xxyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 296);

    auto g_xxyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 297);

    auto g_xxyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 298);

    auto g_xxyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 299);

    auto g_xxyyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 300);

    auto g_xxyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 301);

    auto g_xxyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 302);

    auto g_xxyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 303);

    auto g_xxyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 304);

    auto g_xxyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 305);

    auto g_xxyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 306);

    auto g_xxyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 307);

    auto g_xxyyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 309);

    auto g_xxyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 311);

    auto g_xxyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 314);

    auto g_xxyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 318);

    auto g_xxyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 323);

    auto g_xxyyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 336);

    auto g_xxyyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 337);

    auto g_xxyyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 338);

    auto g_xxyyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 339);

    auto g_xxyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 340);

    auto g_xxyyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 341);

    auto g_xxyyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 342);

    auto g_xxyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 343);

    auto g_xxyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 344);

    auto g_xxyyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 345);

    auto g_xxyyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 346);

    auto g_xxyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 347);

    auto g_xxyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 348);

    auto g_xxyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 349);

    auto g_xxyyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 350);

    auto g_xxyyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 351);

    auto g_xxyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 352);

    auto g_xxyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 353);

    auto g_xxyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 354);

    auto g_xxyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 355);

    auto g_xxyyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 356);

    auto g_xxyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 357);

    auto g_xxyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 358);

    auto g_xxyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 359);

    auto g_xxyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 360);

    auto g_xxyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 361);

    auto g_xxyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 362);

    auto g_xxyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 363);

    auto g_xxyzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 364);

    auto g_xxyzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 366);

    auto g_xxyzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 369);

    auto g_xxyzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 373);

    auto g_xxyzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 378);

    auto g_xxyzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 384);

    auto g_xxzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 392);

    auto g_xxzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 393);

    auto g_xxzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 394);

    auto g_xxzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 395);

    auto g_xxzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 396);

    auto g_xxzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 397);

    auto g_xxzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 398);

    auto g_xxzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 399);

    auto g_xxzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 400);

    auto g_xxzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 401);

    auto g_xxzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 402);

    auto g_xxzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 403);

    auto g_xxzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 404);

    auto g_xxzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 405);

    auto g_xxzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 406);

    auto g_xxzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 407);

    auto g_xxzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 408);

    auto g_xxzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 409);

    auto g_xxzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 410);

    auto g_xxzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 411);

    auto g_xxzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 412);

    auto g_xxzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 413);

    auto g_xxzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 414);

    auto g_xxzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 415);

    auto g_xxzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 416);

    auto g_xxzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 417);

    auto g_xxzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 418);

    auto g_xxzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 419);

    auto g_xyyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 420);

    auto g_xyyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 421);

    auto g_xyyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 423);

    auto g_xyyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 424);

    auto g_xyyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 426);

    auto g_xyyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 427);

    auto g_xyyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 428);

    auto g_xyyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 430);

    auto g_xyyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 431);

    auto g_xyyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 432);

    auto g_xyyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 433);

    auto g_xyyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 435);

    auto g_xyyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 436);

    auto g_xyyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 437);

    auto g_xyyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 438);

    auto g_xyyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 439);

    auto g_xyyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 441);

    auto g_xyyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 442);

    auto g_xyyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 443);

    auto g_xyyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 444);

    auto g_xyyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 445);

    auto g_xyyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 446);

    auto g_xyyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 447);

    auto g_xyyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 480);

    auto g_xyyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 483);

    auto g_xyyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 484);

    auto g_xyyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 487);

    auto g_xyyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 488);

    auto g_xyyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 489);

    auto g_xyyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 492);

    auto g_xyyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 493);

    auto g_xyyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 494);

    auto g_xyyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 495);

    auto g_xyyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 497);

    auto g_xyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 498);

    auto g_xyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 499);

    auto g_xyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 500);

    auto g_xyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 501);

    auto g_xyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 502);

    auto g_xyyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 503);

    auto g_xyyzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 508);

    auto g_xyyzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 511);

    auto g_xyyzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 512);

    auto g_xyyzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 515);

    auto g_xyyzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 516);

    auto g_xyyzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 517);

    auto g_xyyzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 520);

    auto g_xyyzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 521);

    auto g_xyyzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 522);

    auto g_xyyzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 523);

    auto g_xyyzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 525);

    auto g_xyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 526);

    auto g_xyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 527);

    auto g_xyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 528);

    auto g_xyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 529);

    auto g_xyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 530);

    auto g_xyyzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 531);

    auto g_xzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 560);

    auto g_xzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 562);

    auto g_xzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 564);

    auto g_xzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 565);

    auto g_xzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 567);

    auto g_xzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 568);

    auto g_xzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 569);

    auto g_xzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 571);

    auto g_xzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 572);

    auto g_xzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 573);

    auto g_xzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 574);

    auto g_xzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 576);

    auto g_xzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 577);

    auto g_xzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 578);

    auto g_xzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 579);

    auto g_xzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 580);

    auto g_xzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 581);

    auto g_xzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 582);

    auto g_xzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 583);

    auto g_xzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 584);

    auto g_xzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 585);

    auto g_xzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 586);

    auto g_xzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 587);

    auto g_yyyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 588);

    auto g_yyyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 589);

    auto g_yyyyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 590);

    auto g_yyyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 591);

    auto g_yyyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 592);

    auto g_yyyyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 593);

    auto g_yyyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 594);

    auto g_yyyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 595);

    auto g_yyyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 596);

    auto g_yyyyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 597);

    auto g_yyyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 598);

    auto g_yyyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 599);

    auto g_yyyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 600);

    auto g_yyyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 601);

    auto g_yyyyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 602);

    auto g_yyyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 603);

    auto g_yyyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 604);

    auto g_yyyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 605);

    auto g_yyyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 606);

    auto g_yyyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 607);

    auto g_yyyyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 608);

    auto g_yyyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 609);

    auto g_yyyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 610);

    auto g_yyyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 611);

    auto g_yyyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 612);

    auto g_yyyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 613);

    auto g_yyyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 614);

    auto g_yyyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 615);

    auto g_yyyyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 617);

    auto g_yyyyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 618);

    auto g_yyyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 619);

    auto g_yyyyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 620);

    auto g_yyyyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 621);

    auto g_yyyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 622);

    auto g_yyyyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 623);

    auto g_yyyyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 624);

    auto g_yyyyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 625);

    auto g_yyyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 626);

    auto g_yyyyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 627);

    auto g_yyyyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 628);

    auto g_yyyyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 629);

    auto g_yyyyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 630);

    auto g_yyyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 631);

    auto g_yyyyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 632);

    auto g_yyyyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 633);

    auto g_yyyyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 634);

    auto g_yyyyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 635);

    auto g_yyyyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 636);

    auto g_yyyyyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 637);

    auto g_yyyyyz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 638);

    auto g_yyyyyz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 639);

    auto g_yyyyyz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 640);

    auto g_yyyyyz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 641);

    auto g_yyyyyz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 642);

    auto g_yyyyyz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 643);

    auto g_yyyyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 644);

    auto g_yyyyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 645);

    auto g_yyyyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 646);

    auto g_yyyyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 647);

    auto g_yyyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 648);

    auto g_yyyyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 649);

    auto g_yyyyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 650);

    auto g_yyyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 651);

    auto g_yyyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 652);

    auto g_yyyyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 653);

    auto g_yyyyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 654);

    auto g_yyyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 655);

    auto g_yyyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 656);

    auto g_yyyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 657);

    auto g_yyyyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 658);

    auto g_yyyyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 659);

    auto g_yyyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 660);

    auto g_yyyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 661);

    auto g_yyyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 662);

    auto g_yyyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 663);

    auto g_yyyyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 664);

    auto g_yyyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 665);

    auto g_yyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 666);

    auto g_yyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 667);

    auto g_yyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 668);

    auto g_yyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 669);

    auto g_yyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 670);

    auto g_yyyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 671);

    auto g_yyyzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 672);

    auto g_yyyzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 673);

    auto g_yyyzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 674);

    auto g_yyyzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 675);

    auto g_yyyzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 676);

    auto g_yyyzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 677);

    auto g_yyyzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 678);

    auto g_yyyzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 679);

    auto g_yyyzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 680);

    auto g_yyyzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 681);

    auto g_yyyzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 682);

    auto g_yyyzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 683);

    auto g_yyyzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 684);

    auto g_yyyzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 685);

    auto g_yyyzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 686);

    auto g_yyyzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 687);

    auto g_yyyzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 688);

    auto g_yyyzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 689);

    auto g_yyyzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 690);

    auto g_yyyzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 691);

    auto g_yyyzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 692);

    auto g_yyyzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 693);

    auto g_yyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 694);

    auto g_yyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 695);

    auto g_yyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 696);

    auto g_yyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 697);

    auto g_yyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 698);

    auto g_yyyzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 699);

    auto g_yyzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 700);

    auto g_yyzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 701);

    auto g_yyzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 702);

    auto g_yyzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 703);

    auto g_yyzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 704);

    auto g_yyzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 705);

    auto g_yyzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 706);

    auto g_yyzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 707);

    auto g_yyzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 708);

    auto g_yyzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 709);

    auto g_yyzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 710);

    auto g_yyzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 711);

    auto g_yyzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 712);

    auto g_yyzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 713);

    auto g_yyzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 714);

    auto g_yyzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 715);

    auto g_yyzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 716);

    auto g_yyzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 717);

    auto g_yyzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 718);

    auto g_yyzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 719);

    auto g_yyzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 720);

    auto g_yyzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 721);

    auto g_yyzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 722);

    auto g_yyzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 723);

    auto g_yyzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 724);

    auto g_yyzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 725);

    auto g_yyzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 726);

    auto g_yyzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 727);

    auto g_yzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 728);

    auto g_yzzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 729);

    auto g_yzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 730);

    auto g_yzzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 731);

    auto g_yzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 732);

    auto g_yzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 733);

    auto g_yzzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 734);

    auto g_yzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 735);

    auto g_yzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 736);

    auto g_yzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 737);

    auto g_yzzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 738);

    auto g_yzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 739);

    auto g_yzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 740);

    auto g_yzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 741);

    auto g_yzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 742);

    auto g_yzzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 743);

    auto g_yzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 744);

    auto g_yzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 745);

    auto g_yzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 746);

    auto g_yzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 747);

    auto g_yzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 748);

    auto g_yzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 749);

    auto g_yzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 750);

    auto g_yzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 751);

    auto g_yzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 752);

    auto g_yzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 753);

    auto g_yzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 754);

    auto g_yzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 755);

    auto g_zzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_isi + 756);

    auto g_zzzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_isi + 757);

    auto g_zzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 758);

    auto g_zzzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_isi + 759);

    auto g_zzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 760);

    auto g_zzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 761);

    auto g_zzzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_isi + 762);

    auto g_zzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 763);

    auto g_zzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 764);

    auto g_zzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 765);

    auto g_zzzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_isi + 766);

    auto g_zzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 767);

    auto g_zzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 768);

    auto g_zzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 769);

    auto g_zzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 770);

    auto g_zzzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_isi + 771);

    auto g_zzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 772);

    auto g_zzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 773);

    auto g_zzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 774);

    auto g_zzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 775);

    auto g_zzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 776);

    auto g_zzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_isi + 777);

    auto g_zzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 778);

    auto g_zzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 779);

    auto g_zzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 780);

    auto g_zzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 781);

    auto g_zzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 782);

    auto g_zzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_isi + 783);

    /// Set up 0-28 components of targeted buffer : KSI

    auto g_xxxxxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi);

    auto g_xxxxxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 1);

    auto g_xxxxxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 2);

    auto g_xxxxxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 3);

    auto g_xxxxxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 4);

    auto g_xxxxxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 5);

    auto g_xxxxxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 6);

    auto g_xxxxxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 7);

    auto g_xxxxxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 8);

    auto g_xxxxxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 9);

    auto g_xxxxxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 10);

    auto g_xxxxxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 11);

    auto g_xxxxxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 12);

    auto g_xxxxxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 13);

    auto g_xxxxxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 14);

    auto g_xxxxxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 15);

    auto g_xxxxxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 16);

    auto g_xxxxxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 17);

    auto g_xxxxxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 18);

    auto g_xxxxxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 19);

    auto g_xxxxxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 20);

    auto g_xxxxxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 21);

    auto g_xxxxxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 22);

    auto g_xxxxxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 23);

    auto g_xxxxxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 24);

    auto g_xxxxxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 25);

    auto g_xxxxxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 26);

    auto g_xxxxxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 27);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxx_0, g_xxxxx_0_xxxxxx_1, g_xxxxx_0_xxxxxy_0, g_xxxxx_0_xxxxxy_1, g_xxxxx_0_xxxxxz_0, g_xxxxx_0_xxxxxz_1, g_xxxxx_0_xxxxyy_0, g_xxxxx_0_xxxxyy_1, g_xxxxx_0_xxxxyz_0, g_xxxxx_0_xxxxyz_1, g_xxxxx_0_xxxxzz_0, g_xxxxx_0_xxxxzz_1, g_xxxxx_0_xxxyyy_0, g_xxxxx_0_xxxyyy_1, g_xxxxx_0_xxxyyz_0, g_xxxxx_0_xxxyyz_1, g_xxxxx_0_xxxyzz_0, g_xxxxx_0_xxxyzz_1, g_xxxxx_0_xxxzzz_0, g_xxxxx_0_xxxzzz_1, g_xxxxx_0_xxyyyy_0, g_xxxxx_0_xxyyyy_1, g_xxxxx_0_xxyyyz_0, g_xxxxx_0_xxyyyz_1, g_xxxxx_0_xxyyzz_0, g_xxxxx_0_xxyyzz_1, g_xxxxx_0_xxyzzz_0, g_xxxxx_0_xxyzzz_1, g_xxxxx_0_xxzzzz_0, g_xxxxx_0_xxzzzz_1, g_xxxxx_0_xyyyyy_0, g_xxxxx_0_xyyyyy_1, g_xxxxx_0_xyyyyz_0, g_xxxxx_0_xyyyyz_1, g_xxxxx_0_xyyyzz_0, g_xxxxx_0_xyyyzz_1, g_xxxxx_0_xyyzzz_0, g_xxxxx_0_xyyzzz_1, g_xxxxx_0_xyzzzz_0, g_xxxxx_0_xyzzzz_1, g_xxxxx_0_xzzzzz_0, g_xxxxx_0_xzzzzz_1, g_xxxxx_0_yyyyyy_0, g_xxxxx_0_yyyyyy_1, g_xxxxx_0_yyyyyz_0, g_xxxxx_0_yyyyyz_1, g_xxxxx_0_yyyyzz_0, g_xxxxx_0_yyyyzz_1, g_xxxxx_0_yyyzzz_0, g_xxxxx_0_yyyzzz_1, g_xxxxx_0_yyzzzz_0, g_xxxxx_0_yyzzzz_1, g_xxxxx_0_yzzzzz_0, g_xxxxx_0_yzzzzz_1, g_xxxxx_0_zzzzzz_0, g_xxxxx_0_zzzzzz_1, g_xxxxxx_0_xxxxx_1, g_xxxxxx_0_xxxxxx_1, g_xxxxxx_0_xxxxxy_1, g_xxxxxx_0_xxxxxz_1, g_xxxxxx_0_xxxxy_1, g_xxxxxx_0_xxxxyy_1, g_xxxxxx_0_xxxxyz_1, g_xxxxxx_0_xxxxz_1, g_xxxxxx_0_xxxxzz_1, g_xxxxxx_0_xxxyy_1, g_xxxxxx_0_xxxyyy_1, g_xxxxxx_0_xxxyyz_1, g_xxxxxx_0_xxxyz_1, g_xxxxxx_0_xxxyzz_1, g_xxxxxx_0_xxxzz_1, g_xxxxxx_0_xxxzzz_1, g_xxxxxx_0_xxyyy_1, g_xxxxxx_0_xxyyyy_1, g_xxxxxx_0_xxyyyz_1, g_xxxxxx_0_xxyyz_1, g_xxxxxx_0_xxyyzz_1, g_xxxxxx_0_xxyzz_1, g_xxxxxx_0_xxyzzz_1, g_xxxxxx_0_xxzzz_1, g_xxxxxx_0_xxzzzz_1, g_xxxxxx_0_xyyyy_1, g_xxxxxx_0_xyyyyy_1, g_xxxxxx_0_xyyyyz_1, g_xxxxxx_0_xyyyz_1, g_xxxxxx_0_xyyyzz_1, g_xxxxxx_0_xyyzz_1, g_xxxxxx_0_xyyzzz_1, g_xxxxxx_0_xyzzz_1, g_xxxxxx_0_xyzzzz_1, g_xxxxxx_0_xzzzz_1, g_xxxxxx_0_xzzzzz_1, g_xxxxxx_0_yyyyy_1, g_xxxxxx_0_yyyyyy_1, g_xxxxxx_0_yyyyyz_1, g_xxxxxx_0_yyyyz_1, g_xxxxxx_0_yyyyzz_1, g_xxxxxx_0_yyyzz_1, g_xxxxxx_0_yyyzzz_1, g_xxxxxx_0_yyzzz_1, g_xxxxxx_0_yyzzzz_1, g_xxxxxx_0_yzzzz_1, g_xxxxxx_0_yzzzzz_1, g_xxxxxx_0_zzzzz_1, g_xxxxxx_0_zzzzzz_1, g_xxxxxxx_0_xxxxxx_0, g_xxxxxxx_0_xxxxxy_0, g_xxxxxxx_0_xxxxxz_0, g_xxxxxxx_0_xxxxyy_0, g_xxxxxxx_0_xxxxyz_0, g_xxxxxxx_0_xxxxzz_0, g_xxxxxxx_0_xxxyyy_0, g_xxxxxxx_0_xxxyyz_0, g_xxxxxxx_0_xxxyzz_0, g_xxxxxxx_0_xxxzzz_0, g_xxxxxxx_0_xxyyyy_0, g_xxxxxxx_0_xxyyyz_0, g_xxxxxxx_0_xxyyzz_0, g_xxxxxxx_0_xxyzzz_0, g_xxxxxxx_0_xxzzzz_0, g_xxxxxxx_0_xyyyyy_0, g_xxxxxxx_0_xyyyyz_0, g_xxxxxxx_0_xyyyzz_0, g_xxxxxxx_0_xyyzzz_0, g_xxxxxxx_0_xyzzzz_0, g_xxxxxxx_0_xzzzzz_0, g_xxxxxxx_0_yyyyyy_0, g_xxxxxxx_0_yyyyyz_0, g_xxxxxxx_0_yyyyzz_0, g_xxxxxxx_0_yyyzzz_0, g_xxxxxxx_0_yyzzzz_0, g_xxxxxxx_0_yzzzzz_0, g_xxxxxxx_0_zzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxxx_0_xxxxxx_0[i] = 6.0 * g_xxxxx_0_xxxxxx_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxx_1[i] * fz_be_0 + 6.0 * g_xxxxxx_0_xxxxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxx_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxy_0[i] = 6.0 * g_xxxxx_0_xxxxxy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxxxxx_0_xxxxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxz_0[i] = 6.0 * g_xxxxx_0_xxxxxz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxxxxx_0_xxxxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxyy_0[i] = 6.0 * g_xxxxx_0_xxxxyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxyz_0[i] = 6.0 * g_xxxxx_0_xxxxyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxzz_0[i] = 6.0 * g_xxxxx_0_xxxxzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyyy_0[i] = 6.0 * g_xxxxx_0_xxxyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyyz_0[i] = 6.0 * g_xxxxx_0_xxxyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyzz_0[i] = 6.0 * g_xxxxx_0_xxxyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxzzz_0[i] = 6.0 * g_xxxxx_0_xxxzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyyy_0[i] = 6.0 * g_xxxxx_0_xxyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyyz_0[i] = 6.0 * g_xxxxx_0_xxyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyzz_0[i] = 6.0 * g_xxxxx_0_xxyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyzzz_0[i] = 6.0 * g_xxxxx_0_xxyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxzzzz_0[i] = 6.0 * g_xxxxx_0_xxzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyyy_0[i] = 6.0 * g_xxxxx_0_xyyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyyy_1[i] * fz_be_0 + g_xxxxxx_0_yyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyyz_0[i] = 6.0 * g_xxxxx_0_xyyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyyz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyzz_0[i] = 6.0 * g_xxxxx_0_xyyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyzzz_0[i] = 6.0 * g_xxxxx_0_xyyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyzzzz_0[i] = 6.0 * g_xxxxx_0_xyzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xzzzzz_0[i] = 6.0 * g_xxxxx_0_xzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_zzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyyy_0[i] = 6.0 * g_xxxxx_0_yyyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyyy_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyyz_0[i] = 6.0 * g_xxxxx_0_yyyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyyz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyzz_0[i] = 6.0 * g_xxxxx_0_yyyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyzzz_0[i] = 6.0 * g_xxxxx_0_yyyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyzzzz_0[i] = 6.0 * g_xxxxx_0_yyzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yzzzzz_0[i] = 6.0 * g_xxxxx_0_yzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_zzzzzz_0[i] = 6.0 * g_xxxxx_0_zzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_zzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 28-56 components of targeted buffer : KSI

    auto g_xxxxxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 28);

    auto g_xxxxxxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 29);

    auto g_xxxxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 30);

    auto g_xxxxxxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 31);

    auto g_xxxxxxy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 32);

    auto g_xxxxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 33);

    auto g_xxxxxxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 34);

    auto g_xxxxxxy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 35);

    auto g_xxxxxxy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 36);

    auto g_xxxxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 37);

    auto g_xxxxxxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 38);

    auto g_xxxxxxy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 39);

    auto g_xxxxxxy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 40);

    auto g_xxxxxxy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 41);

    auto g_xxxxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 42);

    auto g_xxxxxxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 43);

    auto g_xxxxxxy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 44);

    auto g_xxxxxxy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 45);

    auto g_xxxxxxy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 46);

    auto g_xxxxxxy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 47);

    auto g_xxxxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 48);

    auto g_xxxxxxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 49);

    auto g_xxxxxxy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 50);

    auto g_xxxxxxy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 51);

    auto g_xxxxxxy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 52);

    auto g_xxxxxxy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 53);

    auto g_xxxxxxy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 54);

    auto g_xxxxxxy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 55);

    #pragma omp simd aligned(g_xxxxxx_0_xxxxx_1, g_xxxxxx_0_xxxxxx_1, g_xxxxxx_0_xxxxxy_1, g_xxxxxx_0_xxxxxz_1, g_xxxxxx_0_xxxxy_1, g_xxxxxx_0_xxxxyy_1, g_xxxxxx_0_xxxxyz_1, g_xxxxxx_0_xxxxz_1, g_xxxxxx_0_xxxxzz_1, g_xxxxxx_0_xxxyy_1, g_xxxxxx_0_xxxyyy_1, g_xxxxxx_0_xxxyyz_1, g_xxxxxx_0_xxxyz_1, g_xxxxxx_0_xxxyzz_1, g_xxxxxx_0_xxxzz_1, g_xxxxxx_0_xxxzzz_1, g_xxxxxx_0_xxyyy_1, g_xxxxxx_0_xxyyyy_1, g_xxxxxx_0_xxyyyz_1, g_xxxxxx_0_xxyyz_1, g_xxxxxx_0_xxyyzz_1, g_xxxxxx_0_xxyzz_1, g_xxxxxx_0_xxyzzz_1, g_xxxxxx_0_xxzzz_1, g_xxxxxx_0_xxzzzz_1, g_xxxxxx_0_xyyyy_1, g_xxxxxx_0_xyyyyy_1, g_xxxxxx_0_xyyyyz_1, g_xxxxxx_0_xyyyz_1, g_xxxxxx_0_xyyyzz_1, g_xxxxxx_0_xyyzz_1, g_xxxxxx_0_xyyzzz_1, g_xxxxxx_0_xyzzz_1, g_xxxxxx_0_xyzzzz_1, g_xxxxxx_0_xzzzz_1, g_xxxxxx_0_xzzzzz_1, g_xxxxxx_0_yyyyy_1, g_xxxxxx_0_yyyyyy_1, g_xxxxxx_0_yyyyyz_1, g_xxxxxx_0_yyyyz_1, g_xxxxxx_0_yyyyzz_1, g_xxxxxx_0_yyyzz_1, g_xxxxxx_0_yyyzzz_1, g_xxxxxx_0_yyzzz_1, g_xxxxxx_0_yyzzzz_1, g_xxxxxx_0_yzzzz_1, g_xxxxxx_0_yzzzzz_1, g_xxxxxx_0_zzzzz_1, g_xxxxxx_0_zzzzzz_1, g_xxxxxxy_0_xxxxxx_0, g_xxxxxxy_0_xxxxxy_0, g_xxxxxxy_0_xxxxxz_0, g_xxxxxxy_0_xxxxyy_0, g_xxxxxxy_0_xxxxyz_0, g_xxxxxxy_0_xxxxzz_0, g_xxxxxxy_0_xxxyyy_0, g_xxxxxxy_0_xxxyyz_0, g_xxxxxxy_0_xxxyzz_0, g_xxxxxxy_0_xxxzzz_0, g_xxxxxxy_0_xxyyyy_0, g_xxxxxxy_0_xxyyyz_0, g_xxxxxxy_0_xxyyzz_0, g_xxxxxxy_0_xxyzzz_0, g_xxxxxxy_0_xxzzzz_0, g_xxxxxxy_0_xyyyyy_0, g_xxxxxxy_0_xyyyyz_0, g_xxxxxxy_0_xyyyzz_0, g_xxxxxxy_0_xyyzzz_0, g_xxxxxxy_0_xyzzzz_0, g_xxxxxxy_0_xzzzzz_0, g_xxxxxxy_0_yyyyyy_0, g_xxxxxxy_0_yyyyyz_0, g_xxxxxxy_0_yyyyzz_0, g_xxxxxxy_0_yyyzzz_0, g_xxxxxxy_0_yyzzzz_0, g_xxxxxxy_0_yzzzzz_0, g_xxxxxxy_0_zzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxy_0_xxxxxx_0[i] = g_xxxxxx_0_xxxxxx_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxy_0[i] = g_xxxxxx_0_xxxxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxz_0[i] = g_xxxxxx_0_xxxxxz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxyy_0[i] = 2.0 * g_xxxxxx_0_xxxxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxyz_0[i] = g_xxxxxx_0_xxxxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxzz_0[i] = g_xxxxxx_0_xxxxzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyyy_0[i] = 3.0 * g_xxxxxx_0_xxxyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyyz_0[i] = 2.0 * g_xxxxxx_0_xxxyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyzz_0[i] = g_xxxxxx_0_xxxzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxzzz_0[i] = g_xxxxxx_0_xxxzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyyy_0[i] = 4.0 * g_xxxxxx_0_xxyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyyz_0[i] = 3.0 * g_xxxxxx_0_xxyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyzz_0[i] = 2.0 * g_xxxxxx_0_xxyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyzzz_0[i] = g_xxxxxx_0_xxzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxzzzz_0[i] = g_xxxxxx_0_xxzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyyy_0[i] = 5.0 * g_xxxxxx_0_xyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyyz_0[i] = 4.0 * g_xxxxxx_0_xyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyzz_0[i] = 3.0 * g_xxxxxx_0_xyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyzzz_0[i] = 2.0 * g_xxxxxx_0_xyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyzzzz_0[i] = g_xxxxxx_0_xzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xzzzzz_0[i] = g_xxxxxx_0_xzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyyy_0[i] = 6.0 * g_xxxxxx_0_yyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyyz_0[i] = 5.0 * g_xxxxxx_0_yyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyzz_0[i] = 4.0 * g_xxxxxx_0_yyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyzzz_0[i] = 3.0 * g_xxxxxx_0_yyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyzzzz_0[i] = 2.0 * g_xxxxxx_0_yzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yzzzzz_0[i] = g_xxxxxx_0_zzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_zzzzzz_0[i] = g_xxxxxx_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 56-84 components of targeted buffer : KSI

    auto g_xxxxxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 56);

    auto g_xxxxxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 57);

    auto g_xxxxxxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 58);

    auto g_xxxxxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 59);

    auto g_xxxxxxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 60);

    auto g_xxxxxxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 61);

    auto g_xxxxxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 62);

    auto g_xxxxxxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 63);

    auto g_xxxxxxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 64);

    auto g_xxxxxxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 65);

    auto g_xxxxxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 66);

    auto g_xxxxxxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 67);

    auto g_xxxxxxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 68);

    auto g_xxxxxxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 69);

    auto g_xxxxxxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 70);

    auto g_xxxxxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 71);

    auto g_xxxxxxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 72);

    auto g_xxxxxxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 73);

    auto g_xxxxxxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 74);

    auto g_xxxxxxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 75);

    auto g_xxxxxxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 76);

    auto g_xxxxxxz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 77);

    auto g_xxxxxxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 78);

    auto g_xxxxxxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 79);

    auto g_xxxxxxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 80);

    auto g_xxxxxxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 81);

    auto g_xxxxxxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 82);

    auto g_xxxxxxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 83);

    #pragma omp simd aligned(g_xxxxxx_0_xxxxx_1, g_xxxxxx_0_xxxxxx_1, g_xxxxxx_0_xxxxxy_1, g_xxxxxx_0_xxxxxz_1, g_xxxxxx_0_xxxxy_1, g_xxxxxx_0_xxxxyy_1, g_xxxxxx_0_xxxxyz_1, g_xxxxxx_0_xxxxz_1, g_xxxxxx_0_xxxxzz_1, g_xxxxxx_0_xxxyy_1, g_xxxxxx_0_xxxyyy_1, g_xxxxxx_0_xxxyyz_1, g_xxxxxx_0_xxxyz_1, g_xxxxxx_0_xxxyzz_1, g_xxxxxx_0_xxxzz_1, g_xxxxxx_0_xxxzzz_1, g_xxxxxx_0_xxyyy_1, g_xxxxxx_0_xxyyyy_1, g_xxxxxx_0_xxyyyz_1, g_xxxxxx_0_xxyyz_1, g_xxxxxx_0_xxyyzz_1, g_xxxxxx_0_xxyzz_1, g_xxxxxx_0_xxyzzz_1, g_xxxxxx_0_xxzzz_1, g_xxxxxx_0_xxzzzz_1, g_xxxxxx_0_xyyyy_1, g_xxxxxx_0_xyyyyy_1, g_xxxxxx_0_xyyyyz_1, g_xxxxxx_0_xyyyz_1, g_xxxxxx_0_xyyyzz_1, g_xxxxxx_0_xyyzz_1, g_xxxxxx_0_xyyzzz_1, g_xxxxxx_0_xyzzz_1, g_xxxxxx_0_xyzzzz_1, g_xxxxxx_0_xzzzz_1, g_xxxxxx_0_xzzzzz_1, g_xxxxxx_0_yyyyy_1, g_xxxxxx_0_yyyyyy_1, g_xxxxxx_0_yyyyyz_1, g_xxxxxx_0_yyyyz_1, g_xxxxxx_0_yyyyzz_1, g_xxxxxx_0_yyyzz_1, g_xxxxxx_0_yyyzzz_1, g_xxxxxx_0_yyzzz_1, g_xxxxxx_0_yyzzzz_1, g_xxxxxx_0_yzzzz_1, g_xxxxxx_0_yzzzzz_1, g_xxxxxx_0_zzzzz_1, g_xxxxxx_0_zzzzzz_1, g_xxxxxxz_0_xxxxxx_0, g_xxxxxxz_0_xxxxxy_0, g_xxxxxxz_0_xxxxxz_0, g_xxxxxxz_0_xxxxyy_0, g_xxxxxxz_0_xxxxyz_0, g_xxxxxxz_0_xxxxzz_0, g_xxxxxxz_0_xxxyyy_0, g_xxxxxxz_0_xxxyyz_0, g_xxxxxxz_0_xxxyzz_0, g_xxxxxxz_0_xxxzzz_0, g_xxxxxxz_0_xxyyyy_0, g_xxxxxxz_0_xxyyyz_0, g_xxxxxxz_0_xxyyzz_0, g_xxxxxxz_0_xxyzzz_0, g_xxxxxxz_0_xxzzzz_0, g_xxxxxxz_0_xyyyyy_0, g_xxxxxxz_0_xyyyyz_0, g_xxxxxxz_0_xyyyzz_0, g_xxxxxxz_0_xyyzzz_0, g_xxxxxxz_0_xyzzzz_0, g_xxxxxxz_0_xzzzzz_0, g_xxxxxxz_0_yyyyyy_0, g_xxxxxxz_0_yyyyyz_0, g_xxxxxxz_0_yyyyzz_0, g_xxxxxxz_0_yyyzzz_0, g_xxxxxxz_0_yyzzzz_0, g_xxxxxxz_0_yzzzzz_0, g_xxxxxxz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxz_0_xxxxxx_0[i] = g_xxxxxx_0_xxxxxx_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxy_0[i] = g_xxxxxx_0_xxxxxy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxz_0[i] = g_xxxxxx_0_xxxxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxyy_0[i] = g_xxxxxx_0_xxxxyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxyz_0[i] = g_xxxxxx_0_xxxxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxzz_0[i] = 2.0 * g_xxxxxx_0_xxxxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyyy_0[i] = g_xxxxxx_0_xxxyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyyz_0[i] = g_xxxxxx_0_xxxyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyzz_0[i] = 2.0 * g_xxxxxx_0_xxxyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxzzz_0[i] = 3.0 * g_xxxxxx_0_xxxzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyyy_0[i] = g_xxxxxx_0_xxyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyyz_0[i] = g_xxxxxx_0_xxyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyzz_0[i] = 2.0 * g_xxxxxx_0_xxyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyzzz_0[i] = 3.0 * g_xxxxxx_0_xxyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxzzzz_0[i] = 4.0 * g_xxxxxx_0_xxzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyyy_0[i] = g_xxxxxx_0_xyyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyyz_0[i] = g_xxxxxx_0_xyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyzz_0[i] = 2.0 * g_xxxxxx_0_xyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyzzz_0[i] = 3.0 * g_xxxxxx_0_xyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyzzzz_0[i] = 4.0 * g_xxxxxx_0_xyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xzzzzz_0[i] = 5.0 * g_xxxxxx_0_xzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyyy_0[i] = g_xxxxxx_0_yyyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyyz_0[i] = g_xxxxxx_0_yyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyzz_0[i] = 2.0 * g_xxxxxx_0_yyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyzzz_0[i] = 3.0 * g_xxxxxx_0_yyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyzzzz_0[i] = 4.0 * g_xxxxxx_0_yyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yzzzzz_0[i] = 5.0 * g_xxxxxx_0_yzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_zzzzzz_0[i] = 6.0 * g_xxxxxx_0_zzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 84-112 components of targeted buffer : KSI

    auto g_xxxxxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 84);

    auto g_xxxxxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 85);

    auto g_xxxxxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 86);

    auto g_xxxxxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 87);

    auto g_xxxxxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 88);

    auto g_xxxxxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 89);

    auto g_xxxxxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 90);

    auto g_xxxxxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 91);

    auto g_xxxxxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 92);

    auto g_xxxxxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 93);

    auto g_xxxxxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 94);

    auto g_xxxxxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 95);

    auto g_xxxxxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 96);

    auto g_xxxxxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 97);

    auto g_xxxxxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 98);

    auto g_xxxxxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 99);

    auto g_xxxxxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 100);

    auto g_xxxxxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 101);

    auto g_xxxxxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 102);

    auto g_xxxxxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 103);

    auto g_xxxxxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 104);

    auto g_xxxxxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 105);

    auto g_xxxxxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 106);

    auto g_xxxxxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 107);

    auto g_xxxxxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 108);

    auto g_xxxxxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 109);

    auto g_xxxxxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 110);

    auto g_xxxxxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 111);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxx_0, g_xxxxx_0_xxxxxx_1, g_xxxxx_0_xxxxxz_0, g_xxxxx_0_xxxxxz_1, g_xxxxx_0_xxxxzz_0, g_xxxxx_0_xxxxzz_1, g_xxxxx_0_xxxzzz_0, g_xxxxx_0_xxxzzz_1, g_xxxxx_0_xxzzzz_0, g_xxxxx_0_xxzzzz_1, g_xxxxx_0_xzzzzz_0, g_xxxxx_0_xzzzzz_1, g_xxxxxy_0_xxxxxx_1, g_xxxxxy_0_xxxxxz_1, g_xxxxxy_0_xxxxzz_1, g_xxxxxy_0_xxxzzz_1, g_xxxxxy_0_xxzzzz_1, g_xxxxxy_0_xzzzzz_1, g_xxxxxyy_0_xxxxxx_0, g_xxxxxyy_0_xxxxxy_0, g_xxxxxyy_0_xxxxxz_0, g_xxxxxyy_0_xxxxyy_0, g_xxxxxyy_0_xxxxyz_0, g_xxxxxyy_0_xxxxzz_0, g_xxxxxyy_0_xxxyyy_0, g_xxxxxyy_0_xxxyyz_0, g_xxxxxyy_0_xxxyzz_0, g_xxxxxyy_0_xxxzzz_0, g_xxxxxyy_0_xxyyyy_0, g_xxxxxyy_0_xxyyyz_0, g_xxxxxyy_0_xxyyzz_0, g_xxxxxyy_0_xxyzzz_0, g_xxxxxyy_0_xxzzzz_0, g_xxxxxyy_0_xyyyyy_0, g_xxxxxyy_0_xyyyyz_0, g_xxxxxyy_0_xyyyzz_0, g_xxxxxyy_0_xyyzzz_0, g_xxxxxyy_0_xyzzzz_0, g_xxxxxyy_0_xzzzzz_0, g_xxxxxyy_0_yyyyyy_0, g_xxxxxyy_0_yyyyyz_0, g_xxxxxyy_0_yyyyzz_0, g_xxxxxyy_0_yyyzzz_0, g_xxxxxyy_0_yyzzzz_0, g_xxxxxyy_0_yzzzzz_0, g_xxxxxyy_0_zzzzzz_0, g_xxxxyy_0_xxxxxy_1, g_xxxxyy_0_xxxxy_1, g_xxxxyy_0_xxxxyy_1, g_xxxxyy_0_xxxxyz_1, g_xxxxyy_0_xxxyy_1, g_xxxxyy_0_xxxyyy_1, g_xxxxyy_0_xxxyyz_1, g_xxxxyy_0_xxxyz_1, g_xxxxyy_0_xxxyzz_1, g_xxxxyy_0_xxyyy_1, g_xxxxyy_0_xxyyyy_1, g_xxxxyy_0_xxyyyz_1, g_xxxxyy_0_xxyyz_1, g_xxxxyy_0_xxyyzz_1, g_xxxxyy_0_xxyzz_1, g_xxxxyy_0_xxyzzz_1, g_xxxxyy_0_xyyyy_1, g_xxxxyy_0_xyyyyy_1, g_xxxxyy_0_xyyyyz_1, g_xxxxyy_0_xyyyz_1, g_xxxxyy_0_xyyyzz_1, g_xxxxyy_0_xyyzz_1, g_xxxxyy_0_xyyzzz_1, g_xxxxyy_0_xyzzz_1, g_xxxxyy_0_xyzzzz_1, g_xxxxyy_0_yyyyy_1, g_xxxxyy_0_yyyyyy_1, g_xxxxyy_0_yyyyyz_1, g_xxxxyy_0_yyyyz_1, g_xxxxyy_0_yyyyzz_1, g_xxxxyy_0_yyyzz_1, g_xxxxyy_0_yyyzzz_1, g_xxxxyy_0_yyzzz_1, g_xxxxyy_0_yyzzzz_1, g_xxxxyy_0_yzzzz_1, g_xxxxyy_0_yzzzzz_1, g_xxxxyy_0_zzzzzz_1, g_xxxyy_0_xxxxxy_0, g_xxxyy_0_xxxxxy_1, g_xxxyy_0_xxxxyy_0, g_xxxyy_0_xxxxyy_1, g_xxxyy_0_xxxxyz_0, g_xxxyy_0_xxxxyz_1, g_xxxyy_0_xxxyyy_0, g_xxxyy_0_xxxyyy_1, g_xxxyy_0_xxxyyz_0, g_xxxyy_0_xxxyyz_1, g_xxxyy_0_xxxyzz_0, g_xxxyy_0_xxxyzz_1, g_xxxyy_0_xxyyyy_0, g_xxxyy_0_xxyyyy_1, g_xxxyy_0_xxyyyz_0, g_xxxyy_0_xxyyyz_1, g_xxxyy_0_xxyyzz_0, g_xxxyy_0_xxyyzz_1, g_xxxyy_0_xxyzzz_0, g_xxxyy_0_xxyzzz_1, g_xxxyy_0_xyyyyy_0, g_xxxyy_0_xyyyyy_1, g_xxxyy_0_xyyyyz_0, g_xxxyy_0_xyyyyz_1, g_xxxyy_0_xyyyzz_0, g_xxxyy_0_xyyyzz_1, g_xxxyy_0_xyyzzz_0, g_xxxyy_0_xyyzzz_1, g_xxxyy_0_xyzzzz_0, g_xxxyy_0_xyzzzz_1, g_xxxyy_0_yyyyyy_0, g_xxxyy_0_yyyyyy_1, g_xxxyy_0_yyyyyz_0, g_xxxyy_0_yyyyyz_1, g_xxxyy_0_yyyyzz_0, g_xxxyy_0_yyyyzz_1, g_xxxyy_0_yyyzzz_0, g_xxxyy_0_yyyzzz_1, g_xxxyy_0_yyzzzz_0, g_xxxyy_0_yyzzzz_1, g_xxxyy_0_yzzzzz_0, g_xxxyy_0_yzzzzz_1, g_xxxyy_0_zzzzzz_0, g_xxxyy_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxyy_0_xxxxxx_0[i] = g_xxxxx_0_xxxxxx_0[i] * fbe_0 - g_xxxxx_0_xxxxxx_1[i] * fz_be_0 + g_xxxxxy_0_xxxxxx_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxxxy_0[i] = 4.0 * g_xxxyy_0_xxxxxy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxxxyy_0_xxxxy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxxz_0[i] = g_xxxxx_0_xxxxxz_0[i] * fbe_0 - g_xxxxx_0_xxxxxz_1[i] * fz_be_0 + g_xxxxxy_0_xxxxxz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxxyy_0[i] = 4.0 * g_xxxyy_0_xxxxyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxxxyy_0_xxxyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxyz_0[i] = 4.0 * g_xxxyy_0_xxxxyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxxyy_0_xxxyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxzz_0[i] = g_xxxxx_0_xxxxzz_0[i] * fbe_0 - g_xxxxx_0_xxxxzz_1[i] * fz_be_0 + g_xxxxxy_0_xxxxzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxyyy_0[i] = 4.0 * g_xxxyy_0_xxxyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxyyz_0[i] = 4.0 * g_xxxyy_0_xxxyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxyzz_0[i] = 4.0 * g_xxxyy_0_xxxyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxzzz_0[i] = g_xxxxx_0_xxxzzz_0[i] * fbe_0 - g_xxxxx_0_xxxzzz_1[i] * fz_be_0 + g_xxxxxy_0_xxxzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxyyyy_0[i] = 4.0 * g_xxxyy_0_xxyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyyyz_0[i] = 4.0 * g_xxxyy_0_xxyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyyzz_0[i] = 4.0 * g_xxxyy_0_xxyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyzzz_0[i] = 4.0 * g_xxxyy_0_xxyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxzzzz_0[i] = g_xxxxx_0_xxzzzz_0[i] * fbe_0 - g_xxxxx_0_xxzzzz_1[i] * fz_be_0 + g_xxxxxy_0_xxzzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xyyyyy_0[i] = 4.0 * g_xxxyy_0_xyyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyyy_1[i] * fz_be_0 + g_xxxxyy_0_yyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyyyz_0[i] = 4.0 * g_xxxyy_0_xyyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyyz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyyzz_0[i] = 4.0 * g_xxxyy_0_xyyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyzzz_0[i] = 4.0 * g_xxxyy_0_xyyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyzzzz_0[i] = 4.0 * g_xxxyy_0_xyzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xzzzzz_0[i] = g_xxxxx_0_xzzzzz_0[i] * fbe_0 - g_xxxxx_0_xzzzzz_1[i] * fz_be_0 + g_xxxxxy_0_xzzzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_yyyyyy_0[i] = 4.0 * g_xxxyy_0_yyyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyyy_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyyyz_0[i] = 4.0 * g_xxxyy_0_yyyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyyz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyyzz_0[i] = 4.0 * g_xxxyy_0_yyyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyzzz_0[i] = 4.0 * g_xxxyy_0_yyyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyzzzz_0[i] = 4.0 * g_xxxyy_0_yyzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yzzzzz_0[i] = 4.0 * g_xxxyy_0_yzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yzzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_zzzzzz_0[i] = 4.0 * g_xxxyy_0_zzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_zzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 112-140 components of targeted buffer : KSI

    auto g_xxxxxyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 112);

    auto g_xxxxxyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 113);

    auto g_xxxxxyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 114);

    auto g_xxxxxyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 115);

    auto g_xxxxxyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 116);

    auto g_xxxxxyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 117);

    auto g_xxxxxyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 118);

    auto g_xxxxxyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 119);

    auto g_xxxxxyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 120);

    auto g_xxxxxyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 121);

    auto g_xxxxxyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 122);

    auto g_xxxxxyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 123);

    auto g_xxxxxyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 124);

    auto g_xxxxxyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 125);

    auto g_xxxxxyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 126);

    auto g_xxxxxyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 127);

    auto g_xxxxxyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 128);

    auto g_xxxxxyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 129);

    auto g_xxxxxyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 130);

    auto g_xxxxxyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 131);

    auto g_xxxxxyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 132);

    auto g_xxxxxyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 133);

    auto g_xxxxxyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 134);

    auto g_xxxxxyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 135);

    auto g_xxxxxyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 136);

    auto g_xxxxxyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 137);

    auto g_xxxxxyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 138);

    auto g_xxxxxyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 139);

    #pragma omp simd aligned(g_xxxxxy_0_xxxxxy_1, g_xxxxxy_0_xxxxyy_1, g_xxxxxy_0_xxxyyy_1, g_xxxxxy_0_xxyyyy_1, g_xxxxxy_0_xyyyyy_1, g_xxxxxy_0_yyyyyy_1, g_xxxxxyz_0_xxxxxx_0, g_xxxxxyz_0_xxxxxy_0, g_xxxxxyz_0_xxxxxz_0, g_xxxxxyz_0_xxxxyy_0, g_xxxxxyz_0_xxxxyz_0, g_xxxxxyz_0_xxxxzz_0, g_xxxxxyz_0_xxxyyy_0, g_xxxxxyz_0_xxxyyz_0, g_xxxxxyz_0_xxxyzz_0, g_xxxxxyz_0_xxxzzz_0, g_xxxxxyz_0_xxyyyy_0, g_xxxxxyz_0_xxyyyz_0, g_xxxxxyz_0_xxyyzz_0, g_xxxxxyz_0_xxyzzz_0, g_xxxxxyz_0_xxzzzz_0, g_xxxxxyz_0_xyyyyy_0, g_xxxxxyz_0_xyyyyz_0, g_xxxxxyz_0_xyyyzz_0, g_xxxxxyz_0_xyyzzz_0, g_xxxxxyz_0_xyzzzz_0, g_xxxxxyz_0_xzzzzz_0, g_xxxxxyz_0_yyyyyy_0, g_xxxxxyz_0_yyyyyz_0, g_xxxxxyz_0_yyyyzz_0, g_xxxxxyz_0_yyyzzz_0, g_xxxxxyz_0_yyzzzz_0, g_xxxxxyz_0_yzzzzz_0, g_xxxxxyz_0_zzzzzz_0, g_xxxxxz_0_xxxxxx_1, g_xxxxxz_0_xxxxxz_1, g_xxxxxz_0_xxxxyz_1, g_xxxxxz_0_xxxxz_1, g_xxxxxz_0_xxxxzz_1, g_xxxxxz_0_xxxyyz_1, g_xxxxxz_0_xxxyz_1, g_xxxxxz_0_xxxyzz_1, g_xxxxxz_0_xxxzz_1, g_xxxxxz_0_xxxzzz_1, g_xxxxxz_0_xxyyyz_1, g_xxxxxz_0_xxyyz_1, g_xxxxxz_0_xxyyzz_1, g_xxxxxz_0_xxyzz_1, g_xxxxxz_0_xxyzzz_1, g_xxxxxz_0_xxzzz_1, g_xxxxxz_0_xxzzzz_1, g_xxxxxz_0_xyyyyz_1, g_xxxxxz_0_xyyyz_1, g_xxxxxz_0_xyyyzz_1, g_xxxxxz_0_xyyzz_1, g_xxxxxz_0_xyyzzz_1, g_xxxxxz_0_xyzzz_1, g_xxxxxz_0_xyzzzz_1, g_xxxxxz_0_xzzzz_1, g_xxxxxz_0_xzzzzz_1, g_xxxxxz_0_yyyyyz_1, g_xxxxxz_0_yyyyz_1, g_xxxxxz_0_yyyyzz_1, g_xxxxxz_0_yyyzz_1, g_xxxxxz_0_yyyzzz_1, g_xxxxxz_0_yyzzz_1, g_xxxxxz_0_yyzzzz_1, g_xxxxxz_0_yzzzz_1, g_xxxxxz_0_yzzzzz_1, g_xxxxxz_0_zzzzz_1, g_xxxxxz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxyz_0_xxxxxx_0[i] = g_xxxxxz_0_xxxxxx_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxxy_0[i] = g_xxxxxy_0_xxxxxy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxxxz_0[i] = g_xxxxxz_0_xxxxxz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxyy_0[i] = g_xxxxxy_0_xxxxyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxxyz_0[i] = g_xxxxxz_0_xxxxz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxxyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxzz_0[i] = g_xxxxxz_0_xxxxzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxyyy_0[i] = g_xxxxxy_0_xxxyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxyyz_0[i] = 2.0 * g_xxxxxz_0_xxxyz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxyzz_0[i] = g_xxxxxz_0_xxxzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxzzz_0[i] = g_xxxxxz_0_xxxzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyyyy_0[i] = g_xxxxxy_0_xxyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxyyyz_0[i] = 3.0 * g_xxxxxz_0_xxyyz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyyzz_0[i] = 2.0 * g_xxxxxz_0_xxyzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyzzz_0[i] = g_xxxxxz_0_xxzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxzzzz_0[i] = g_xxxxxz_0_xxzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyyyy_0[i] = g_xxxxxy_0_xyyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xyyyyz_0[i] = 4.0 * g_xxxxxz_0_xyyyz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyyzz_0[i] = 3.0 * g_xxxxxz_0_xyyzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyzzz_0[i] = 2.0 * g_xxxxxz_0_xyzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyzzzz_0[i] = g_xxxxxz_0_xzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xzzzzz_0[i] = g_xxxxxz_0_xzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyyyy_0[i] = g_xxxxxy_0_yyyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_yyyyyz_0[i] = 5.0 * g_xxxxxz_0_yyyyz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyyzz_0[i] = 4.0 * g_xxxxxz_0_yyyzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyzzz_0[i] = 3.0 * g_xxxxxz_0_yyzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyzzzz_0[i] = 2.0 * g_xxxxxz_0_yzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yzzzzz_0[i] = g_xxxxxz_0_zzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_zzzzzz_0[i] = g_xxxxxz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 140-168 components of targeted buffer : KSI

    auto g_xxxxxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 140);

    auto g_xxxxxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 141);

    auto g_xxxxxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 142);

    auto g_xxxxxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 143);

    auto g_xxxxxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 144);

    auto g_xxxxxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 145);

    auto g_xxxxxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 146);

    auto g_xxxxxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 147);

    auto g_xxxxxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 148);

    auto g_xxxxxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 149);

    auto g_xxxxxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 150);

    auto g_xxxxxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 151);

    auto g_xxxxxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 152);

    auto g_xxxxxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 153);

    auto g_xxxxxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 154);

    auto g_xxxxxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 155);

    auto g_xxxxxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 156);

    auto g_xxxxxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 157);

    auto g_xxxxxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 158);

    auto g_xxxxxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 159);

    auto g_xxxxxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 160);

    auto g_xxxxxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 161);

    auto g_xxxxxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 162);

    auto g_xxxxxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 163);

    auto g_xxxxxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 164);

    auto g_xxxxxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 165);

    auto g_xxxxxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 166);

    auto g_xxxxxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 167);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxx_0, g_xxxxx_0_xxxxxx_1, g_xxxxx_0_xxxxxy_0, g_xxxxx_0_xxxxxy_1, g_xxxxx_0_xxxxyy_0, g_xxxxx_0_xxxxyy_1, g_xxxxx_0_xxxyyy_0, g_xxxxx_0_xxxyyy_1, g_xxxxx_0_xxyyyy_0, g_xxxxx_0_xxyyyy_1, g_xxxxx_0_xyyyyy_0, g_xxxxx_0_xyyyyy_1, g_xxxxxz_0_xxxxxx_1, g_xxxxxz_0_xxxxxy_1, g_xxxxxz_0_xxxxyy_1, g_xxxxxz_0_xxxyyy_1, g_xxxxxz_0_xxyyyy_1, g_xxxxxz_0_xyyyyy_1, g_xxxxxzz_0_xxxxxx_0, g_xxxxxzz_0_xxxxxy_0, g_xxxxxzz_0_xxxxxz_0, g_xxxxxzz_0_xxxxyy_0, g_xxxxxzz_0_xxxxyz_0, g_xxxxxzz_0_xxxxzz_0, g_xxxxxzz_0_xxxyyy_0, g_xxxxxzz_0_xxxyyz_0, g_xxxxxzz_0_xxxyzz_0, g_xxxxxzz_0_xxxzzz_0, g_xxxxxzz_0_xxyyyy_0, g_xxxxxzz_0_xxyyyz_0, g_xxxxxzz_0_xxyyzz_0, g_xxxxxzz_0_xxyzzz_0, g_xxxxxzz_0_xxzzzz_0, g_xxxxxzz_0_xyyyyy_0, g_xxxxxzz_0_xyyyyz_0, g_xxxxxzz_0_xyyyzz_0, g_xxxxxzz_0_xyyzzz_0, g_xxxxxzz_0_xyzzzz_0, g_xxxxxzz_0_xzzzzz_0, g_xxxxxzz_0_yyyyyy_0, g_xxxxxzz_0_yyyyyz_0, g_xxxxxzz_0_yyyyzz_0, g_xxxxxzz_0_yyyzzz_0, g_xxxxxzz_0_yyzzzz_0, g_xxxxxzz_0_yzzzzz_0, g_xxxxxzz_0_zzzzzz_0, g_xxxxzz_0_xxxxxz_1, g_xxxxzz_0_xxxxyz_1, g_xxxxzz_0_xxxxz_1, g_xxxxzz_0_xxxxzz_1, g_xxxxzz_0_xxxyyz_1, g_xxxxzz_0_xxxyz_1, g_xxxxzz_0_xxxyzz_1, g_xxxxzz_0_xxxzz_1, g_xxxxzz_0_xxxzzz_1, g_xxxxzz_0_xxyyyz_1, g_xxxxzz_0_xxyyz_1, g_xxxxzz_0_xxyyzz_1, g_xxxxzz_0_xxyzz_1, g_xxxxzz_0_xxyzzz_1, g_xxxxzz_0_xxzzz_1, g_xxxxzz_0_xxzzzz_1, g_xxxxzz_0_xyyyyz_1, g_xxxxzz_0_xyyyz_1, g_xxxxzz_0_xyyyzz_1, g_xxxxzz_0_xyyzz_1, g_xxxxzz_0_xyyzzz_1, g_xxxxzz_0_xyzzz_1, g_xxxxzz_0_xyzzzz_1, g_xxxxzz_0_xzzzz_1, g_xxxxzz_0_xzzzzz_1, g_xxxxzz_0_yyyyyy_1, g_xxxxzz_0_yyyyyz_1, g_xxxxzz_0_yyyyz_1, g_xxxxzz_0_yyyyzz_1, g_xxxxzz_0_yyyzz_1, g_xxxxzz_0_yyyzzz_1, g_xxxxzz_0_yyzzz_1, g_xxxxzz_0_yyzzzz_1, g_xxxxzz_0_yzzzz_1, g_xxxxzz_0_yzzzzz_1, g_xxxxzz_0_zzzzz_1, g_xxxxzz_0_zzzzzz_1, g_xxxzz_0_xxxxxz_0, g_xxxzz_0_xxxxxz_1, g_xxxzz_0_xxxxyz_0, g_xxxzz_0_xxxxyz_1, g_xxxzz_0_xxxxzz_0, g_xxxzz_0_xxxxzz_1, g_xxxzz_0_xxxyyz_0, g_xxxzz_0_xxxyyz_1, g_xxxzz_0_xxxyzz_0, g_xxxzz_0_xxxyzz_1, g_xxxzz_0_xxxzzz_0, g_xxxzz_0_xxxzzz_1, g_xxxzz_0_xxyyyz_0, g_xxxzz_0_xxyyyz_1, g_xxxzz_0_xxyyzz_0, g_xxxzz_0_xxyyzz_1, g_xxxzz_0_xxyzzz_0, g_xxxzz_0_xxyzzz_1, g_xxxzz_0_xxzzzz_0, g_xxxzz_0_xxzzzz_1, g_xxxzz_0_xyyyyz_0, g_xxxzz_0_xyyyyz_1, g_xxxzz_0_xyyyzz_0, g_xxxzz_0_xyyyzz_1, g_xxxzz_0_xyyzzz_0, g_xxxzz_0_xyyzzz_1, g_xxxzz_0_xyzzzz_0, g_xxxzz_0_xyzzzz_1, g_xxxzz_0_xzzzzz_0, g_xxxzz_0_xzzzzz_1, g_xxxzz_0_yyyyyy_0, g_xxxzz_0_yyyyyy_1, g_xxxzz_0_yyyyyz_0, g_xxxzz_0_yyyyyz_1, g_xxxzz_0_yyyyzz_0, g_xxxzz_0_yyyyzz_1, g_xxxzz_0_yyyzzz_0, g_xxxzz_0_yyyzzz_1, g_xxxzz_0_yyzzzz_0, g_xxxzz_0_yyzzzz_1, g_xxxzz_0_yzzzzz_0, g_xxxzz_0_yzzzzz_1, g_xxxzz_0_zzzzzz_0, g_xxxzz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxzz_0_xxxxxx_0[i] = g_xxxxx_0_xxxxxx_0[i] * fbe_0 - g_xxxxx_0_xxxxxx_1[i] * fz_be_0 + g_xxxxxz_0_xxxxxx_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxxy_0[i] = g_xxxxx_0_xxxxxy_0[i] * fbe_0 - g_xxxxx_0_xxxxxy_1[i] * fz_be_0 + g_xxxxxz_0_xxxxxy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxxz_0[i] = 4.0 * g_xxxzz_0_xxxxxz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxxxzz_0_xxxxz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxyy_0[i] = g_xxxxx_0_xxxxyy_0[i] * fbe_0 - g_xxxxx_0_xxxxyy_1[i] * fz_be_0 + g_xxxxxz_0_xxxxyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxyz_0[i] = 4.0 * g_xxxzz_0_xxxxyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_0_xxxyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxzz_0[i] = 4.0 * g_xxxzz_0_xxxxzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_0_xxxzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxyyy_0[i] = g_xxxxx_0_xxxyyy_0[i] * fbe_0 - g_xxxxx_0_xxxyyy_1[i] * fz_be_0 + g_xxxxxz_0_xxxyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxyyz_0[i] = 4.0 * g_xxxzz_0_xxxyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxyzz_0[i] = 4.0 * g_xxxzz_0_xxxyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxzzz_0[i] = 4.0 * g_xxxzz_0_xxxzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyyyy_0[i] = g_xxxxx_0_xxyyyy_0[i] * fbe_0 - g_xxxxx_0_xxyyyy_1[i] * fz_be_0 + g_xxxxxz_0_xxyyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxyyyz_0[i] = 4.0 * g_xxxzz_0_xxyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyyzz_0[i] = 4.0 * g_xxxzz_0_xxyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyzzz_0[i] = 4.0 * g_xxxzz_0_xxyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxzzzz_0[i] = 4.0 * g_xxxzz_0_xxzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyyyy_0[i] = g_xxxxx_0_xyyyyy_0[i] * fbe_0 - g_xxxxx_0_xyyyyy_1[i] * fz_be_0 + g_xxxxxz_0_xyyyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xyyyyz_0[i] = 4.0 * g_xxxzz_0_xyyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyyyz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyyzz_0[i] = 4.0 * g_xxxzz_0_xyyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyyzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyzzz_0[i] = 4.0 * g_xxxzz_0_xyyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyzzzz_0[i] = 4.0 * g_xxxzz_0_xyzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xzzzzz_0[i] = 4.0 * g_xxxzz_0_xzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_zzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyyy_0[i] = 4.0 * g_xxxzz_0_yyyyyy_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyyy_1[i] * fz_be_0 + g_xxxxzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyyz_0[i] = 4.0 * g_xxxzz_0_yyyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyyz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyzz_0[i] = 4.0 * g_xxxzz_0_yyyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyzzz_0[i] = 4.0 * g_xxxzz_0_yyyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyzzzz_0[i] = 4.0 * g_xxxzz_0_yyzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yzzzzz_0[i] = 4.0 * g_xxxzz_0_yzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_zzzzzz_0[i] = 4.0 * g_xxxzz_0_zzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_zzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 168-196 components of targeted buffer : KSI

    auto g_xxxxyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 168);

    auto g_xxxxyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 169);

    auto g_xxxxyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 170);

    auto g_xxxxyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 171);

    auto g_xxxxyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 172);

    auto g_xxxxyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 173);

    auto g_xxxxyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 174);

    auto g_xxxxyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 175);

    auto g_xxxxyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 176);

    auto g_xxxxyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 177);

    auto g_xxxxyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 178);

    auto g_xxxxyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 179);

    auto g_xxxxyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 180);

    auto g_xxxxyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 181);

    auto g_xxxxyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 182);

    auto g_xxxxyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 183);

    auto g_xxxxyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 184);

    auto g_xxxxyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 185);

    auto g_xxxxyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 186);

    auto g_xxxxyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 187);

    auto g_xxxxyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 188);

    auto g_xxxxyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 189);

    auto g_xxxxyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 190);

    auto g_xxxxyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 191);

    auto g_xxxxyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 192);

    auto g_xxxxyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 193);

    auto g_xxxxyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 194);

    auto g_xxxxyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 195);

    #pragma omp simd aligned(g_xxxxy_0_xxxxxx_0, g_xxxxy_0_xxxxxx_1, g_xxxxy_0_xxxxxz_0, g_xxxxy_0_xxxxxz_1, g_xxxxy_0_xxxxzz_0, g_xxxxy_0_xxxxzz_1, g_xxxxy_0_xxxzzz_0, g_xxxxy_0_xxxzzz_1, g_xxxxy_0_xxzzzz_0, g_xxxxy_0_xxzzzz_1, g_xxxxy_0_xzzzzz_0, g_xxxxy_0_xzzzzz_1, g_xxxxyy_0_xxxxxx_1, g_xxxxyy_0_xxxxxz_1, g_xxxxyy_0_xxxxzz_1, g_xxxxyy_0_xxxzzz_1, g_xxxxyy_0_xxzzzz_1, g_xxxxyy_0_xzzzzz_1, g_xxxxyyy_0_xxxxxx_0, g_xxxxyyy_0_xxxxxy_0, g_xxxxyyy_0_xxxxxz_0, g_xxxxyyy_0_xxxxyy_0, g_xxxxyyy_0_xxxxyz_0, g_xxxxyyy_0_xxxxzz_0, g_xxxxyyy_0_xxxyyy_0, g_xxxxyyy_0_xxxyyz_0, g_xxxxyyy_0_xxxyzz_0, g_xxxxyyy_0_xxxzzz_0, g_xxxxyyy_0_xxyyyy_0, g_xxxxyyy_0_xxyyyz_0, g_xxxxyyy_0_xxyyzz_0, g_xxxxyyy_0_xxyzzz_0, g_xxxxyyy_0_xxzzzz_0, g_xxxxyyy_0_xyyyyy_0, g_xxxxyyy_0_xyyyyz_0, g_xxxxyyy_0_xyyyzz_0, g_xxxxyyy_0_xyyzzz_0, g_xxxxyyy_0_xyzzzz_0, g_xxxxyyy_0_xzzzzz_0, g_xxxxyyy_0_yyyyyy_0, g_xxxxyyy_0_yyyyyz_0, g_xxxxyyy_0_yyyyzz_0, g_xxxxyyy_0_yyyzzz_0, g_xxxxyyy_0_yyzzzz_0, g_xxxxyyy_0_yzzzzz_0, g_xxxxyyy_0_zzzzzz_0, g_xxxyyy_0_xxxxxy_1, g_xxxyyy_0_xxxxy_1, g_xxxyyy_0_xxxxyy_1, g_xxxyyy_0_xxxxyz_1, g_xxxyyy_0_xxxyy_1, g_xxxyyy_0_xxxyyy_1, g_xxxyyy_0_xxxyyz_1, g_xxxyyy_0_xxxyz_1, g_xxxyyy_0_xxxyzz_1, g_xxxyyy_0_xxyyy_1, g_xxxyyy_0_xxyyyy_1, g_xxxyyy_0_xxyyyz_1, g_xxxyyy_0_xxyyz_1, g_xxxyyy_0_xxyyzz_1, g_xxxyyy_0_xxyzz_1, g_xxxyyy_0_xxyzzz_1, g_xxxyyy_0_xyyyy_1, g_xxxyyy_0_xyyyyy_1, g_xxxyyy_0_xyyyyz_1, g_xxxyyy_0_xyyyz_1, g_xxxyyy_0_xyyyzz_1, g_xxxyyy_0_xyyzz_1, g_xxxyyy_0_xyyzzz_1, g_xxxyyy_0_xyzzz_1, g_xxxyyy_0_xyzzzz_1, g_xxxyyy_0_yyyyy_1, g_xxxyyy_0_yyyyyy_1, g_xxxyyy_0_yyyyyz_1, g_xxxyyy_0_yyyyz_1, g_xxxyyy_0_yyyyzz_1, g_xxxyyy_0_yyyzz_1, g_xxxyyy_0_yyyzzz_1, g_xxxyyy_0_yyzzz_1, g_xxxyyy_0_yyzzzz_1, g_xxxyyy_0_yzzzz_1, g_xxxyyy_0_yzzzzz_1, g_xxxyyy_0_zzzzzz_1, g_xxyyy_0_xxxxxy_0, g_xxyyy_0_xxxxxy_1, g_xxyyy_0_xxxxyy_0, g_xxyyy_0_xxxxyy_1, g_xxyyy_0_xxxxyz_0, g_xxyyy_0_xxxxyz_1, g_xxyyy_0_xxxyyy_0, g_xxyyy_0_xxxyyy_1, g_xxyyy_0_xxxyyz_0, g_xxyyy_0_xxxyyz_1, g_xxyyy_0_xxxyzz_0, g_xxyyy_0_xxxyzz_1, g_xxyyy_0_xxyyyy_0, g_xxyyy_0_xxyyyy_1, g_xxyyy_0_xxyyyz_0, g_xxyyy_0_xxyyyz_1, g_xxyyy_0_xxyyzz_0, g_xxyyy_0_xxyyzz_1, g_xxyyy_0_xxyzzz_0, g_xxyyy_0_xxyzzz_1, g_xxyyy_0_xyyyyy_0, g_xxyyy_0_xyyyyy_1, g_xxyyy_0_xyyyyz_0, g_xxyyy_0_xyyyyz_1, g_xxyyy_0_xyyyzz_0, g_xxyyy_0_xyyyzz_1, g_xxyyy_0_xyyzzz_0, g_xxyyy_0_xyyzzz_1, g_xxyyy_0_xyzzzz_0, g_xxyyy_0_xyzzzz_1, g_xxyyy_0_yyyyyy_0, g_xxyyy_0_yyyyyy_1, g_xxyyy_0_yyyyyz_0, g_xxyyy_0_yyyyyz_1, g_xxyyy_0_yyyyzz_0, g_xxyyy_0_yyyyzz_1, g_xxyyy_0_yyyzzz_0, g_xxyyy_0_yyyzzz_1, g_xxyyy_0_yyzzzz_0, g_xxyyy_0_yyzzzz_1, g_xxyyy_0_yzzzzz_0, g_xxyyy_0_yzzzzz_1, g_xxyyy_0_zzzzzz_0, g_xxyyy_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyyy_0_xxxxxx_0[i] = 2.0 * g_xxxxy_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxxx_1[i] * fz_be_0 + g_xxxxyy_0_xxxxxx_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxxxy_0[i] = 3.0 * g_xxyyy_0_xxxxxy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxxyyy_0_xxxxy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxxz_0[i] = 2.0 * g_xxxxy_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxxz_1[i] * fz_be_0 + g_xxxxyy_0_xxxxxz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxxyy_0[i] = 3.0 * g_xxyyy_0_xxxxyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxxyyy_0_xxxyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxyz_0[i] = 3.0 * g_xxyyy_0_xxxxyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxyyy_0_xxxyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxzz_0[i] = 2.0 * g_xxxxy_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxzz_1[i] * fz_be_0 + g_xxxxyy_0_xxxxzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxyyy_0[i] = 3.0 * g_xxyyy_0_xxxyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxyyz_0[i] = 3.0 * g_xxyyy_0_xxxyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxyzz_0[i] = 3.0 * g_xxyyy_0_xxxyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxzzz_0[i] = 2.0 * g_xxxxy_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxzzz_1[i] * fz_be_0 + g_xxxxyy_0_xxxzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxyyyy_0[i] = 3.0 * g_xxyyy_0_xxyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyyyz_0[i] = 3.0 * g_xxyyy_0_xxyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyyzz_0[i] = 3.0 * g_xxyyy_0_xxyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyzzz_0[i] = 3.0 * g_xxyyy_0_xxyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxzzzz_0[i] = 2.0 * g_xxxxy_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxzzzz_1[i] * fz_be_0 + g_xxxxyy_0_xxzzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xyyyyy_0[i] = 3.0 * g_xxyyy_0_xyyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyyy_1[i] * fz_be_0 + g_xxxyyy_0_yyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyyyz_0[i] = 3.0 * g_xxyyy_0_xyyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyyz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyyzz_0[i] = 3.0 * g_xxyyy_0_xyyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyzzz_0[i] = 3.0 * g_xxyyy_0_xyyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyzzzz_0[i] = 3.0 * g_xxyyy_0_xyzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xzzzzz_0[i] = 2.0 * g_xxxxy_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_xzzzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_yyyyyy_0[i] = 3.0 * g_xxyyy_0_yyyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyyy_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyyyz_0[i] = 3.0 * g_xxyyy_0_yyyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyyz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyyzz_0[i] = 3.0 * g_xxyyy_0_yyyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyzzz_0[i] = 3.0 * g_xxyyy_0_yyyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyzzzz_0[i] = 3.0 * g_xxyyy_0_yyzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yzzzzz_0[i] = 3.0 * g_xxyyy_0_yzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yzzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_zzzzzz_0[i] = 3.0 * g_xxyyy_0_zzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_zzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 196-224 components of targeted buffer : KSI

    auto g_xxxxyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 196);

    auto g_xxxxyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 197);

    auto g_xxxxyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 198);

    auto g_xxxxyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 199);

    auto g_xxxxyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 200);

    auto g_xxxxyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 201);

    auto g_xxxxyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 202);

    auto g_xxxxyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 203);

    auto g_xxxxyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 204);

    auto g_xxxxyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 205);

    auto g_xxxxyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 206);

    auto g_xxxxyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 207);

    auto g_xxxxyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 208);

    auto g_xxxxyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 209);

    auto g_xxxxyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 210);

    auto g_xxxxyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 211);

    auto g_xxxxyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 212);

    auto g_xxxxyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 213);

    auto g_xxxxyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 214);

    auto g_xxxxyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 215);

    auto g_xxxxyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 216);

    auto g_xxxxyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 217);

    auto g_xxxxyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 218);

    auto g_xxxxyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 219);

    auto g_xxxxyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 220);

    auto g_xxxxyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 221);

    auto g_xxxxyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 222);

    auto g_xxxxyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 223);

    #pragma omp simd aligned(g_xxxxyy_0_xxxxx_1, g_xxxxyy_0_xxxxxx_1, g_xxxxyy_0_xxxxxy_1, g_xxxxyy_0_xxxxxz_1, g_xxxxyy_0_xxxxy_1, g_xxxxyy_0_xxxxyy_1, g_xxxxyy_0_xxxxyz_1, g_xxxxyy_0_xxxxz_1, g_xxxxyy_0_xxxxzz_1, g_xxxxyy_0_xxxyy_1, g_xxxxyy_0_xxxyyy_1, g_xxxxyy_0_xxxyyz_1, g_xxxxyy_0_xxxyz_1, g_xxxxyy_0_xxxyzz_1, g_xxxxyy_0_xxxzz_1, g_xxxxyy_0_xxxzzz_1, g_xxxxyy_0_xxyyy_1, g_xxxxyy_0_xxyyyy_1, g_xxxxyy_0_xxyyyz_1, g_xxxxyy_0_xxyyz_1, g_xxxxyy_0_xxyyzz_1, g_xxxxyy_0_xxyzz_1, g_xxxxyy_0_xxyzzz_1, g_xxxxyy_0_xxzzz_1, g_xxxxyy_0_xxzzzz_1, g_xxxxyy_0_xyyyy_1, g_xxxxyy_0_xyyyyy_1, g_xxxxyy_0_xyyyyz_1, g_xxxxyy_0_xyyyz_1, g_xxxxyy_0_xyyyzz_1, g_xxxxyy_0_xyyzz_1, g_xxxxyy_0_xyyzzz_1, g_xxxxyy_0_xyzzz_1, g_xxxxyy_0_xyzzzz_1, g_xxxxyy_0_xzzzz_1, g_xxxxyy_0_xzzzzz_1, g_xxxxyy_0_yyyyy_1, g_xxxxyy_0_yyyyyy_1, g_xxxxyy_0_yyyyyz_1, g_xxxxyy_0_yyyyz_1, g_xxxxyy_0_yyyyzz_1, g_xxxxyy_0_yyyzz_1, g_xxxxyy_0_yyyzzz_1, g_xxxxyy_0_yyzzz_1, g_xxxxyy_0_yyzzzz_1, g_xxxxyy_0_yzzzz_1, g_xxxxyy_0_yzzzzz_1, g_xxxxyy_0_zzzzz_1, g_xxxxyy_0_zzzzzz_1, g_xxxxyyz_0_xxxxxx_0, g_xxxxyyz_0_xxxxxy_0, g_xxxxyyz_0_xxxxxz_0, g_xxxxyyz_0_xxxxyy_0, g_xxxxyyz_0_xxxxyz_0, g_xxxxyyz_0_xxxxzz_0, g_xxxxyyz_0_xxxyyy_0, g_xxxxyyz_0_xxxyyz_0, g_xxxxyyz_0_xxxyzz_0, g_xxxxyyz_0_xxxzzz_0, g_xxxxyyz_0_xxyyyy_0, g_xxxxyyz_0_xxyyyz_0, g_xxxxyyz_0_xxyyzz_0, g_xxxxyyz_0_xxyzzz_0, g_xxxxyyz_0_xxzzzz_0, g_xxxxyyz_0_xyyyyy_0, g_xxxxyyz_0_xyyyyz_0, g_xxxxyyz_0_xyyyzz_0, g_xxxxyyz_0_xyyzzz_0, g_xxxxyyz_0_xyzzzz_0, g_xxxxyyz_0_xzzzzz_0, g_xxxxyyz_0_yyyyyy_0, g_xxxxyyz_0_yyyyyz_0, g_xxxxyyz_0_yyyyzz_0, g_xxxxyyz_0_yyyzzz_0, g_xxxxyyz_0_yyzzzz_0, g_xxxxyyz_0_yzzzzz_0, g_xxxxyyz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyyz_0_xxxxxx_0[i] = g_xxxxyy_0_xxxxxx_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxy_0[i] = g_xxxxyy_0_xxxxxy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxz_0[i] = g_xxxxyy_0_xxxxx_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxyy_0[i] = g_xxxxyy_0_xxxxyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxyz_0[i] = g_xxxxyy_0_xxxxy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxzz_0[i] = 2.0 * g_xxxxyy_0_xxxxz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyyy_0[i] = g_xxxxyy_0_xxxyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyyz_0[i] = g_xxxxyy_0_xxxyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyzz_0[i] = 2.0 * g_xxxxyy_0_xxxyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxzzz_0[i] = 3.0 * g_xxxxyy_0_xxxzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyyy_0[i] = g_xxxxyy_0_xxyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyyz_0[i] = g_xxxxyy_0_xxyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyzz_0[i] = 2.0 * g_xxxxyy_0_xxyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyzzz_0[i] = 3.0 * g_xxxxyy_0_xxyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxzzzz_0[i] = 4.0 * g_xxxxyy_0_xxzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyyy_0[i] = g_xxxxyy_0_xyyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyyz_0[i] = g_xxxxyy_0_xyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyzz_0[i] = 2.0 * g_xxxxyy_0_xyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyzzz_0[i] = 3.0 * g_xxxxyy_0_xyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyzzzz_0[i] = 4.0 * g_xxxxyy_0_xyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xzzzzz_0[i] = 5.0 * g_xxxxyy_0_xzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyyy_0[i] = g_xxxxyy_0_yyyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyyz_0[i] = g_xxxxyy_0_yyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_yyyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyzz_0[i] = 2.0 * g_xxxxyy_0_yyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_yyyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyzzz_0[i] = 3.0 * g_xxxxyy_0_yyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_yyyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyzzzz_0[i] = 4.0 * g_xxxxyy_0_yyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_yyzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yzzzzz_0[i] = 5.0 * g_xxxxyy_0_yzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_yzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_zzzzzz_0[i] = 6.0 * g_xxxxyy_0_zzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 224-252 components of targeted buffer : KSI

    auto g_xxxxyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 224);

    auto g_xxxxyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 225);

    auto g_xxxxyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 226);

    auto g_xxxxyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 227);

    auto g_xxxxyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 228);

    auto g_xxxxyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 229);

    auto g_xxxxyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 230);

    auto g_xxxxyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 231);

    auto g_xxxxyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 232);

    auto g_xxxxyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 233);

    auto g_xxxxyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 234);

    auto g_xxxxyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 235);

    auto g_xxxxyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 236);

    auto g_xxxxyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 237);

    auto g_xxxxyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 238);

    auto g_xxxxyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 239);

    auto g_xxxxyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 240);

    auto g_xxxxyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 241);

    auto g_xxxxyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 242);

    auto g_xxxxyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 243);

    auto g_xxxxyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 244);

    auto g_xxxxyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 245);

    auto g_xxxxyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 246);

    auto g_xxxxyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 247);

    auto g_xxxxyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 248);

    auto g_xxxxyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 249);

    auto g_xxxxyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 250);

    auto g_xxxxyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 251);

    #pragma omp simd aligned(g_xxxxyzz_0_xxxxxx_0, g_xxxxyzz_0_xxxxxy_0, g_xxxxyzz_0_xxxxxz_0, g_xxxxyzz_0_xxxxyy_0, g_xxxxyzz_0_xxxxyz_0, g_xxxxyzz_0_xxxxzz_0, g_xxxxyzz_0_xxxyyy_0, g_xxxxyzz_0_xxxyyz_0, g_xxxxyzz_0_xxxyzz_0, g_xxxxyzz_0_xxxzzz_0, g_xxxxyzz_0_xxyyyy_0, g_xxxxyzz_0_xxyyyz_0, g_xxxxyzz_0_xxyyzz_0, g_xxxxyzz_0_xxyzzz_0, g_xxxxyzz_0_xxzzzz_0, g_xxxxyzz_0_xyyyyy_0, g_xxxxyzz_0_xyyyyz_0, g_xxxxyzz_0_xyyyzz_0, g_xxxxyzz_0_xyyzzz_0, g_xxxxyzz_0_xyzzzz_0, g_xxxxyzz_0_xzzzzz_0, g_xxxxyzz_0_yyyyyy_0, g_xxxxyzz_0_yyyyyz_0, g_xxxxyzz_0_yyyyzz_0, g_xxxxyzz_0_yyyzzz_0, g_xxxxyzz_0_yyzzzz_0, g_xxxxyzz_0_yzzzzz_0, g_xxxxyzz_0_zzzzzz_0, g_xxxxzz_0_xxxxx_1, g_xxxxzz_0_xxxxxx_1, g_xxxxzz_0_xxxxxy_1, g_xxxxzz_0_xxxxxz_1, g_xxxxzz_0_xxxxy_1, g_xxxxzz_0_xxxxyy_1, g_xxxxzz_0_xxxxyz_1, g_xxxxzz_0_xxxxz_1, g_xxxxzz_0_xxxxzz_1, g_xxxxzz_0_xxxyy_1, g_xxxxzz_0_xxxyyy_1, g_xxxxzz_0_xxxyyz_1, g_xxxxzz_0_xxxyz_1, g_xxxxzz_0_xxxyzz_1, g_xxxxzz_0_xxxzz_1, g_xxxxzz_0_xxxzzz_1, g_xxxxzz_0_xxyyy_1, g_xxxxzz_0_xxyyyy_1, g_xxxxzz_0_xxyyyz_1, g_xxxxzz_0_xxyyz_1, g_xxxxzz_0_xxyyzz_1, g_xxxxzz_0_xxyzz_1, g_xxxxzz_0_xxyzzz_1, g_xxxxzz_0_xxzzz_1, g_xxxxzz_0_xxzzzz_1, g_xxxxzz_0_xyyyy_1, g_xxxxzz_0_xyyyyy_1, g_xxxxzz_0_xyyyyz_1, g_xxxxzz_0_xyyyz_1, g_xxxxzz_0_xyyyzz_1, g_xxxxzz_0_xyyzz_1, g_xxxxzz_0_xyyzzz_1, g_xxxxzz_0_xyzzz_1, g_xxxxzz_0_xyzzzz_1, g_xxxxzz_0_xzzzz_1, g_xxxxzz_0_xzzzzz_1, g_xxxxzz_0_yyyyy_1, g_xxxxzz_0_yyyyyy_1, g_xxxxzz_0_yyyyyz_1, g_xxxxzz_0_yyyyz_1, g_xxxxzz_0_yyyyzz_1, g_xxxxzz_0_yyyzz_1, g_xxxxzz_0_yyyzzz_1, g_xxxxzz_0_yyzzz_1, g_xxxxzz_0_yyzzzz_1, g_xxxxzz_0_yzzzz_1, g_xxxxzz_0_yzzzzz_1, g_xxxxzz_0_zzzzz_1, g_xxxxzz_0_zzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyzz_0_xxxxxx_0[i] = g_xxxxzz_0_xxxxxx_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxy_0[i] = g_xxxxzz_0_xxxxx_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxz_0[i] = g_xxxxzz_0_xxxxxz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxyy_0[i] = 2.0 * g_xxxxzz_0_xxxxy_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxyz_0[i] = g_xxxxzz_0_xxxxz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxzz_0[i] = g_xxxxzz_0_xxxxzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyyy_0[i] = 3.0 * g_xxxxzz_0_xxxyy_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyyz_0[i] = 2.0 * g_xxxxzz_0_xxxyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyzz_0[i] = g_xxxxzz_0_xxxzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxzzz_0[i] = g_xxxxzz_0_xxxzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyyy_0[i] = 4.0 * g_xxxxzz_0_xxyyy_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyyz_0[i] = 3.0 * g_xxxxzz_0_xxyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyzz_0[i] = 2.0 * g_xxxxzz_0_xxyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyzzz_0[i] = g_xxxxzz_0_xxzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxzzzz_0[i] = g_xxxxzz_0_xxzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyyy_0[i] = 5.0 * g_xxxxzz_0_xyyyy_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyyz_0[i] = 4.0 * g_xxxxzz_0_xyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyzz_0[i] = 3.0 * g_xxxxzz_0_xyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyzzz_0[i] = 2.0 * g_xxxxzz_0_xyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyzzzz_0[i] = g_xxxxzz_0_xzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xzzzzz_0[i] = g_xxxxzz_0_xzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyyy_0[i] = 6.0 * g_xxxxzz_0_yyyyy_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyyz_0[i] = 5.0 * g_xxxxzz_0_yyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyzz_0[i] = 4.0 * g_xxxxzz_0_yyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyzzz_0[i] = 3.0 * g_xxxxzz_0_yyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyzzzz_0[i] = 2.0 * g_xxxxzz_0_yzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yzzzzz_0[i] = g_xxxxzz_0_zzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_zzzzzz_0[i] = g_xxxxzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 252-280 components of targeted buffer : KSI

    auto g_xxxxzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 252);

    auto g_xxxxzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 253);

    auto g_xxxxzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 254);

    auto g_xxxxzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 255);

    auto g_xxxxzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 256);

    auto g_xxxxzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 257);

    auto g_xxxxzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 258);

    auto g_xxxxzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 259);

    auto g_xxxxzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 260);

    auto g_xxxxzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 261);

    auto g_xxxxzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 262);

    auto g_xxxxzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 263);

    auto g_xxxxzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 264);

    auto g_xxxxzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 265);

    auto g_xxxxzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 266);

    auto g_xxxxzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 267);

    auto g_xxxxzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 268);

    auto g_xxxxzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 269);

    auto g_xxxxzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 270);

    auto g_xxxxzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 271);

    auto g_xxxxzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 272);

    auto g_xxxxzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 273);

    auto g_xxxxzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 274);

    auto g_xxxxzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 275);

    auto g_xxxxzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 276);

    auto g_xxxxzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 277);

    auto g_xxxxzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 278);

    auto g_xxxxzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 279);

    #pragma omp simd aligned(g_xxxxz_0_xxxxxx_0, g_xxxxz_0_xxxxxx_1, g_xxxxz_0_xxxxxy_0, g_xxxxz_0_xxxxxy_1, g_xxxxz_0_xxxxyy_0, g_xxxxz_0_xxxxyy_1, g_xxxxz_0_xxxyyy_0, g_xxxxz_0_xxxyyy_1, g_xxxxz_0_xxyyyy_0, g_xxxxz_0_xxyyyy_1, g_xxxxz_0_xyyyyy_0, g_xxxxz_0_xyyyyy_1, g_xxxxzz_0_xxxxxx_1, g_xxxxzz_0_xxxxxy_1, g_xxxxzz_0_xxxxyy_1, g_xxxxzz_0_xxxyyy_1, g_xxxxzz_0_xxyyyy_1, g_xxxxzz_0_xyyyyy_1, g_xxxxzzz_0_xxxxxx_0, g_xxxxzzz_0_xxxxxy_0, g_xxxxzzz_0_xxxxxz_0, g_xxxxzzz_0_xxxxyy_0, g_xxxxzzz_0_xxxxyz_0, g_xxxxzzz_0_xxxxzz_0, g_xxxxzzz_0_xxxyyy_0, g_xxxxzzz_0_xxxyyz_0, g_xxxxzzz_0_xxxyzz_0, g_xxxxzzz_0_xxxzzz_0, g_xxxxzzz_0_xxyyyy_0, g_xxxxzzz_0_xxyyyz_0, g_xxxxzzz_0_xxyyzz_0, g_xxxxzzz_0_xxyzzz_0, g_xxxxzzz_0_xxzzzz_0, g_xxxxzzz_0_xyyyyy_0, g_xxxxzzz_0_xyyyyz_0, g_xxxxzzz_0_xyyyzz_0, g_xxxxzzz_0_xyyzzz_0, g_xxxxzzz_0_xyzzzz_0, g_xxxxzzz_0_xzzzzz_0, g_xxxxzzz_0_yyyyyy_0, g_xxxxzzz_0_yyyyyz_0, g_xxxxzzz_0_yyyyzz_0, g_xxxxzzz_0_yyyzzz_0, g_xxxxzzz_0_yyzzzz_0, g_xxxxzzz_0_yzzzzz_0, g_xxxxzzz_0_zzzzzz_0, g_xxxzzz_0_xxxxxz_1, g_xxxzzz_0_xxxxyz_1, g_xxxzzz_0_xxxxz_1, g_xxxzzz_0_xxxxzz_1, g_xxxzzz_0_xxxyyz_1, g_xxxzzz_0_xxxyz_1, g_xxxzzz_0_xxxyzz_1, g_xxxzzz_0_xxxzz_1, g_xxxzzz_0_xxxzzz_1, g_xxxzzz_0_xxyyyz_1, g_xxxzzz_0_xxyyz_1, g_xxxzzz_0_xxyyzz_1, g_xxxzzz_0_xxyzz_1, g_xxxzzz_0_xxyzzz_1, g_xxxzzz_0_xxzzz_1, g_xxxzzz_0_xxzzzz_1, g_xxxzzz_0_xyyyyz_1, g_xxxzzz_0_xyyyz_1, g_xxxzzz_0_xyyyzz_1, g_xxxzzz_0_xyyzz_1, g_xxxzzz_0_xyyzzz_1, g_xxxzzz_0_xyzzz_1, g_xxxzzz_0_xyzzzz_1, g_xxxzzz_0_xzzzz_1, g_xxxzzz_0_xzzzzz_1, g_xxxzzz_0_yyyyyy_1, g_xxxzzz_0_yyyyyz_1, g_xxxzzz_0_yyyyz_1, g_xxxzzz_0_yyyyzz_1, g_xxxzzz_0_yyyzz_1, g_xxxzzz_0_yyyzzz_1, g_xxxzzz_0_yyzzz_1, g_xxxzzz_0_yyzzzz_1, g_xxxzzz_0_yzzzz_1, g_xxxzzz_0_yzzzzz_1, g_xxxzzz_0_zzzzz_1, g_xxxzzz_0_zzzzzz_1, g_xxzzz_0_xxxxxz_0, g_xxzzz_0_xxxxxz_1, g_xxzzz_0_xxxxyz_0, g_xxzzz_0_xxxxyz_1, g_xxzzz_0_xxxxzz_0, g_xxzzz_0_xxxxzz_1, g_xxzzz_0_xxxyyz_0, g_xxzzz_0_xxxyyz_1, g_xxzzz_0_xxxyzz_0, g_xxzzz_0_xxxyzz_1, g_xxzzz_0_xxxzzz_0, g_xxzzz_0_xxxzzz_1, g_xxzzz_0_xxyyyz_0, g_xxzzz_0_xxyyyz_1, g_xxzzz_0_xxyyzz_0, g_xxzzz_0_xxyyzz_1, g_xxzzz_0_xxyzzz_0, g_xxzzz_0_xxyzzz_1, g_xxzzz_0_xxzzzz_0, g_xxzzz_0_xxzzzz_1, g_xxzzz_0_xyyyyz_0, g_xxzzz_0_xyyyyz_1, g_xxzzz_0_xyyyzz_0, g_xxzzz_0_xyyyzz_1, g_xxzzz_0_xyyzzz_0, g_xxzzz_0_xyyzzz_1, g_xxzzz_0_xyzzzz_0, g_xxzzz_0_xyzzzz_1, g_xxzzz_0_xzzzzz_0, g_xxzzz_0_xzzzzz_1, g_xxzzz_0_yyyyyy_0, g_xxzzz_0_yyyyyy_1, g_xxzzz_0_yyyyyz_0, g_xxzzz_0_yyyyyz_1, g_xxzzz_0_yyyyzz_0, g_xxzzz_0_yyyyzz_1, g_xxzzz_0_yyyzzz_0, g_xxzzz_0_yyyzzz_1, g_xxzzz_0_yyzzzz_0, g_xxzzz_0_yyzzzz_1, g_xxzzz_0_yzzzzz_0, g_xxzzz_0_yzzzzz_1, g_xxzzz_0_zzzzzz_0, g_xxzzz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzzz_0_xxxxxx_0[i] = 2.0 * g_xxxxz_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxxx_1[i] * fz_be_0 + g_xxxxzz_0_xxxxxx_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxxy_0[i] = 2.0 * g_xxxxz_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxxy_1[i] * fz_be_0 + g_xxxxzz_0_xxxxxy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxxz_0[i] = 3.0 * g_xxzzz_0_xxxxxz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxxzzz_0_xxxxz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxyy_0[i] = 2.0 * g_xxxxz_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxyy_1[i] * fz_be_0 + g_xxxxzz_0_xxxxyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxyz_0[i] = 3.0 * g_xxzzz_0_xxxxyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_0_xxxyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxzz_0[i] = 3.0 * g_xxzzz_0_xxxxzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_0_xxxzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxyyy_0[i] = 2.0 * g_xxxxz_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxyyy_1[i] * fz_be_0 + g_xxxxzz_0_xxxyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxyyz_0[i] = 3.0 * g_xxzzz_0_xxxyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxyzz_0[i] = 3.0 * g_xxzzz_0_xxxyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxzzz_0[i] = 3.0 * g_xxzzz_0_xxxzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyyyy_0[i] = 2.0 * g_xxxxz_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxyyyy_1[i] * fz_be_0 + g_xxxxzz_0_xxyyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxyyyz_0[i] = 3.0 * g_xxzzz_0_xxyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyyzz_0[i] = 3.0 * g_xxzzz_0_xxyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyzzz_0[i] = 3.0 * g_xxzzz_0_xxyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxzzzz_0[i] = 3.0 * g_xxzzz_0_xxzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyyyy_0[i] = 2.0 * g_xxxxz_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xyyyyy_1[i] * fz_be_0 + g_xxxxzz_0_xyyyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xyyyyz_0[i] = 3.0 * g_xxzzz_0_xyyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyyyz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyyzz_0[i] = 3.0 * g_xxzzz_0_xyyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyyzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyzzz_0[i] = 3.0 * g_xxzzz_0_xyyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyzzzz_0[i] = 3.0 * g_xxzzz_0_xyzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xzzzzz_0[i] = 3.0 * g_xxzzz_0_xzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_zzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyyy_0[i] = 3.0 * g_xxzzz_0_yyyyyy_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyyy_1[i] * fz_be_0 + g_xxxzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyyz_0[i] = 3.0 * g_xxzzz_0_yyyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyyz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyzz_0[i] = 3.0 * g_xxzzz_0_yyyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyzzz_0[i] = 3.0 * g_xxzzz_0_yyyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyzzzz_0[i] = 3.0 * g_xxzzz_0_yyzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yzzzzz_0[i] = 3.0 * g_xxzzz_0_yzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_zzzzzz_0[i] = 3.0 * g_xxzzz_0_zzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_zzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 280-308 components of targeted buffer : KSI

    auto g_xxxyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 280);

    auto g_xxxyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 281);

    auto g_xxxyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 282);

    auto g_xxxyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 283);

    auto g_xxxyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 284);

    auto g_xxxyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 285);

    auto g_xxxyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 286);

    auto g_xxxyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 287);

    auto g_xxxyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 288);

    auto g_xxxyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 289);

    auto g_xxxyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 290);

    auto g_xxxyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 291);

    auto g_xxxyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 292);

    auto g_xxxyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 293);

    auto g_xxxyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 294);

    auto g_xxxyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 295);

    auto g_xxxyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 296);

    auto g_xxxyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 297);

    auto g_xxxyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 298);

    auto g_xxxyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 299);

    auto g_xxxyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 300);

    auto g_xxxyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 301);

    auto g_xxxyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 302);

    auto g_xxxyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 303);

    auto g_xxxyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 304);

    auto g_xxxyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 305);

    auto g_xxxyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 306);

    auto g_xxxyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 307);

    #pragma omp simd aligned(g_xxxyy_0_xxxxxx_0, g_xxxyy_0_xxxxxx_1, g_xxxyy_0_xxxxxz_0, g_xxxyy_0_xxxxxz_1, g_xxxyy_0_xxxxzz_0, g_xxxyy_0_xxxxzz_1, g_xxxyy_0_xxxzzz_0, g_xxxyy_0_xxxzzz_1, g_xxxyy_0_xxzzzz_0, g_xxxyy_0_xxzzzz_1, g_xxxyy_0_xzzzzz_0, g_xxxyy_0_xzzzzz_1, g_xxxyyy_0_xxxxxx_1, g_xxxyyy_0_xxxxxz_1, g_xxxyyy_0_xxxxzz_1, g_xxxyyy_0_xxxzzz_1, g_xxxyyy_0_xxzzzz_1, g_xxxyyy_0_xzzzzz_1, g_xxxyyyy_0_xxxxxx_0, g_xxxyyyy_0_xxxxxy_0, g_xxxyyyy_0_xxxxxz_0, g_xxxyyyy_0_xxxxyy_0, g_xxxyyyy_0_xxxxyz_0, g_xxxyyyy_0_xxxxzz_0, g_xxxyyyy_0_xxxyyy_0, g_xxxyyyy_0_xxxyyz_0, g_xxxyyyy_0_xxxyzz_0, g_xxxyyyy_0_xxxzzz_0, g_xxxyyyy_0_xxyyyy_0, g_xxxyyyy_0_xxyyyz_0, g_xxxyyyy_0_xxyyzz_0, g_xxxyyyy_0_xxyzzz_0, g_xxxyyyy_0_xxzzzz_0, g_xxxyyyy_0_xyyyyy_0, g_xxxyyyy_0_xyyyyz_0, g_xxxyyyy_0_xyyyzz_0, g_xxxyyyy_0_xyyzzz_0, g_xxxyyyy_0_xyzzzz_0, g_xxxyyyy_0_xzzzzz_0, g_xxxyyyy_0_yyyyyy_0, g_xxxyyyy_0_yyyyyz_0, g_xxxyyyy_0_yyyyzz_0, g_xxxyyyy_0_yyyzzz_0, g_xxxyyyy_0_yyzzzz_0, g_xxxyyyy_0_yzzzzz_0, g_xxxyyyy_0_zzzzzz_0, g_xxyyyy_0_xxxxxy_1, g_xxyyyy_0_xxxxy_1, g_xxyyyy_0_xxxxyy_1, g_xxyyyy_0_xxxxyz_1, g_xxyyyy_0_xxxyy_1, g_xxyyyy_0_xxxyyy_1, g_xxyyyy_0_xxxyyz_1, g_xxyyyy_0_xxxyz_1, g_xxyyyy_0_xxxyzz_1, g_xxyyyy_0_xxyyy_1, g_xxyyyy_0_xxyyyy_1, g_xxyyyy_0_xxyyyz_1, g_xxyyyy_0_xxyyz_1, g_xxyyyy_0_xxyyzz_1, g_xxyyyy_0_xxyzz_1, g_xxyyyy_0_xxyzzz_1, g_xxyyyy_0_xyyyy_1, g_xxyyyy_0_xyyyyy_1, g_xxyyyy_0_xyyyyz_1, g_xxyyyy_0_xyyyz_1, g_xxyyyy_0_xyyyzz_1, g_xxyyyy_0_xyyzz_1, g_xxyyyy_0_xyyzzz_1, g_xxyyyy_0_xyzzz_1, g_xxyyyy_0_xyzzzz_1, g_xxyyyy_0_yyyyy_1, g_xxyyyy_0_yyyyyy_1, g_xxyyyy_0_yyyyyz_1, g_xxyyyy_0_yyyyz_1, g_xxyyyy_0_yyyyzz_1, g_xxyyyy_0_yyyzz_1, g_xxyyyy_0_yyyzzz_1, g_xxyyyy_0_yyzzz_1, g_xxyyyy_0_yyzzzz_1, g_xxyyyy_0_yzzzz_1, g_xxyyyy_0_yzzzzz_1, g_xxyyyy_0_zzzzzz_1, g_xyyyy_0_xxxxxy_0, g_xyyyy_0_xxxxxy_1, g_xyyyy_0_xxxxyy_0, g_xyyyy_0_xxxxyy_1, g_xyyyy_0_xxxxyz_0, g_xyyyy_0_xxxxyz_1, g_xyyyy_0_xxxyyy_0, g_xyyyy_0_xxxyyy_1, g_xyyyy_0_xxxyyz_0, g_xyyyy_0_xxxyyz_1, g_xyyyy_0_xxxyzz_0, g_xyyyy_0_xxxyzz_1, g_xyyyy_0_xxyyyy_0, g_xyyyy_0_xxyyyy_1, g_xyyyy_0_xxyyyz_0, g_xyyyy_0_xxyyyz_1, g_xyyyy_0_xxyyzz_0, g_xyyyy_0_xxyyzz_1, g_xyyyy_0_xxyzzz_0, g_xyyyy_0_xxyzzz_1, g_xyyyy_0_xyyyyy_0, g_xyyyy_0_xyyyyy_1, g_xyyyy_0_xyyyyz_0, g_xyyyy_0_xyyyyz_1, g_xyyyy_0_xyyyzz_0, g_xyyyy_0_xyyyzz_1, g_xyyyy_0_xyyzzz_0, g_xyyyy_0_xyyzzz_1, g_xyyyy_0_xyzzzz_0, g_xyyyy_0_xyzzzz_1, g_xyyyy_0_yyyyyy_0, g_xyyyy_0_yyyyyy_1, g_xyyyy_0_yyyyyz_0, g_xyyyy_0_yyyyyz_1, g_xyyyy_0_yyyyzz_0, g_xyyyy_0_yyyyzz_1, g_xyyyy_0_yyyzzz_0, g_xyyyy_0_yyyzzz_1, g_xyyyy_0_yyzzzz_0, g_xyyyy_0_yyzzzz_1, g_xyyyy_0_yzzzzz_0, g_xyyyy_0_yzzzzz_1, g_xyyyy_0_zzzzzz_0, g_xyyyy_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyyy_0_xxxxxx_0[i] = 3.0 * g_xxxyy_0_xxxxxx_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxxx_1[i] * fz_be_0 + g_xxxyyy_0_xxxxxx_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxxxy_0[i] = 2.0 * g_xyyyy_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxyyyy_0_xxxxy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxxz_0[i] = 3.0 * g_xxxyy_0_xxxxxz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxxz_1[i] * fz_be_0 + g_xxxyyy_0_xxxxxz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxxyy_0[i] = 2.0 * g_xyyyy_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxyyyy_0_xxxyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxyz_0[i] = 2.0 * g_xyyyy_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxyyyy_0_xxxyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxzz_0[i] = 3.0 * g_xxxyy_0_xxxxzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxzz_1[i] * fz_be_0 + g_xxxyyy_0_xxxxzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxyyy_0[i] = 2.0 * g_xyyyy_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxyyz_0[i] = 2.0 * g_xyyyy_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxyzz_0[i] = 2.0 * g_xyyyy_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxzzz_0[i] = 3.0 * g_xxxyy_0_xxxzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxzzz_1[i] * fz_be_0 + g_xxxyyy_0_xxxzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxyyyy_0[i] = 2.0 * g_xyyyy_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyyyz_0[i] = 2.0 * g_xyyyy_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyyzz_0[i] = 2.0 * g_xyyyy_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyzzz_0[i] = 2.0 * g_xyyyy_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxzzzz_0[i] = 3.0 * g_xxxyy_0_xxzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxzzzz_1[i] * fz_be_0 + g_xxxyyy_0_xxzzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xyyyyy_0[i] = 2.0 * g_xyyyy_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyyy_1[i] * fz_be_0 + g_xxyyyy_0_yyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyyyz_0[i] = 2.0 * g_xyyyy_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyyz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyyzz_0[i] = 2.0 * g_xyyyy_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyzzz_0[i] = 2.0 * g_xyyyy_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyzzzz_0[i] = 2.0 * g_xyyyy_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xzzzzz_0[i] = 3.0 * g_xxxyy_0_xzzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_xzzzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_yyyyyy_0[i] = 2.0 * g_xyyyy_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyyy_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyyyz_0[i] = 2.0 * g_xyyyy_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyyz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyyzz_0[i] = 2.0 * g_xyyyy_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyzzz_0[i] = 2.0 * g_xyyyy_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyzzzz_0[i] = 2.0 * g_xyyyy_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yzzzzz_0[i] = 2.0 * g_xyyyy_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yzzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_zzzzzz_0[i] = 2.0 * g_xyyyy_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_zzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 308-336 components of targeted buffer : KSI

    auto g_xxxyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 308);

    auto g_xxxyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 309);

    auto g_xxxyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 310);

    auto g_xxxyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 311);

    auto g_xxxyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 312);

    auto g_xxxyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 313);

    auto g_xxxyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 314);

    auto g_xxxyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 315);

    auto g_xxxyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 316);

    auto g_xxxyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 317);

    auto g_xxxyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 318);

    auto g_xxxyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 319);

    auto g_xxxyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 320);

    auto g_xxxyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 321);

    auto g_xxxyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 322);

    auto g_xxxyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 323);

    auto g_xxxyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 324);

    auto g_xxxyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 325);

    auto g_xxxyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 326);

    auto g_xxxyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 327);

    auto g_xxxyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 328);

    auto g_xxxyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 329);

    auto g_xxxyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 330);

    auto g_xxxyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 331);

    auto g_xxxyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 332);

    auto g_xxxyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 333);

    auto g_xxxyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 334);

    auto g_xxxyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 335);

    #pragma omp simd aligned(g_xxxyyy_0_xxxxx_1, g_xxxyyy_0_xxxxxx_1, g_xxxyyy_0_xxxxxy_1, g_xxxyyy_0_xxxxxz_1, g_xxxyyy_0_xxxxy_1, g_xxxyyy_0_xxxxyy_1, g_xxxyyy_0_xxxxyz_1, g_xxxyyy_0_xxxxz_1, g_xxxyyy_0_xxxxzz_1, g_xxxyyy_0_xxxyy_1, g_xxxyyy_0_xxxyyy_1, g_xxxyyy_0_xxxyyz_1, g_xxxyyy_0_xxxyz_1, g_xxxyyy_0_xxxyzz_1, g_xxxyyy_0_xxxzz_1, g_xxxyyy_0_xxxzzz_1, g_xxxyyy_0_xxyyy_1, g_xxxyyy_0_xxyyyy_1, g_xxxyyy_0_xxyyyz_1, g_xxxyyy_0_xxyyz_1, g_xxxyyy_0_xxyyzz_1, g_xxxyyy_0_xxyzz_1, g_xxxyyy_0_xxyzzz_1, g_xxxyyy_0_xxzzz_1, g_xxxyyy_0_xxzzzz_1, g_xxxyyy_0_xyyyy_1, g_xxxyyy_0_xyyyyy_1, g_xxxyyy_0_xyyyyz_1, g_xxxyyy_0_xyyyz_1, g_xxxyyy_0_xyyyzz_1, g_xxxyyy_0_xyyzz_1, g_xxxyyy_0_xyyzzz_1, g_xxxyyy_0_xyzzz_1, g_xxxyyy_0_xyzzzz_1, g_xxxyyy_0_xzzzz_1, g_xxxyyy_0_xzzzzz_1, g_xxxyyy_0_yyyyy_1, g_xxxyyy_0_yyyyyy_1, g_xxxyyy_0_yyyyyz_1, g_xxxyyy_0_yyyyz_1, g_xxxyyy_0_yyyyzz_1, g_xxxyyy_0_yyyzz_1, g_xxxyyy_0_yyyzzz_1, g_xxxyyy_0_yyzzz_1, g_xxxyyy_0_yyzzzz_1, g_xxxyyy_0_yzzzz_1, g_xxxyyy_0_yzzzzz_1, g_xxxyyy_0_zzzzz_1, g_xxxyyy_0_zzzzzz_1, g_xxxyyyz_0_xxxxxx_0, g_xxxyyyz_0_xxxxxy_0, g_xxxyyyz_0_xxxxxz_0, g_xxxyyyz_0_xxxxyy_0, g_xxxyyyz_0_xxxxyz_0, g_xxxyyyz_0_xxxxzz_0, g_xxxyyyz_0_xxxyyy_0, g_xxxyyyz_0_xxxyyz_0, g_xxxyyyz_0_xxxyzz_0, g_xxxyyyz_0_xxxzzz_0, g_xxxyyyz_0_xxyyyy_0, g_xxxyyyz_0_xxyyyz_0, g_xxxyyyz_0_xxyyzz_0, g_xxxyyyz_0_xxyzzz_0, g_xxxyyyz_0_xxzzzz_0, g_xxxyyyz_0_xyyyyy_0, g_xxxyyyz_0_xyyyyz_0, g_xxxyyyz_0_xyyyzz_0, g_xxxyyyz_0_xyyzzz_0, g_xxxyyyz_0_xyzzzz_0, g_xxxyyyz_0_xzzzzz_0, g_xxxyyyz_0_yyyyyy_0, g_xxxyyyz_0_yyyyyz_0, g_xxxyyyz_0_yyyyzz_0, g_xxxyyyz_0_yyyzzz_0, g_xxxyyyz_0_yyzzzz_0, g_xxxyyyz_0_yzzzzz_0, g_xxxyyyz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyyz_0_xxxxxx_0[i] = g_xxxyyy_0_xxxxxx_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxy_0[i] = g_xxxyyy_0_xxxxxy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxz_0[i] = g_xxxyyy_0_xxxxx_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxyy_0[i] = g_xxxyyy_0_xxxxyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxyz_0[i] = g_xxxyyy_0_xxxxy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxzz_0[i] = 2.0 * g_xxxyyy_0_xxxxz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyyy_0[i] = g_xxxyyy_0_xxxyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyyz_0[i] = g_xxxyyy_0_xxxyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyzz_0[i] = 2.0 * g_xxxyyy_0_xxxyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxzzz_0[i] = 3.0 * g_xxxyyy_0_xxxzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyyy_0[i] = g_xxxyyy_0_xxyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyyz_0[i] = g_xxxyyy_0_xxyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyzz_0[i] = 2.0 * g_xxxyyy_0_xxyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyzzz_0[i] = 3.0 * g_xxxyyy_0_xxyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxzzzz_0[i] = 4.0 * g_xxxyyy_0_xxzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyyy_0[i] = g_xxxyyy_0_xyyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyyz_0[i] = g_xxxyyy_0_xyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyzz_0[i] = 2.0 * g_xxxyyy_0_xyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyzzz_0[i] = 3.0 * g_xxxyyy_0_xyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyzzzz_0[i] = 4.0 * g_xxxyyy_0_xyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xzzzzz_0[i] = 5.0 * g_xxxyyy_0_xzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyyy_0[i] = g_xxxyyy_0_yyyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyyz_0[i] = g_xxxyyy_0_yyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_yyyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyzz_0[i] = 2.0 * g_xxxyyy_0_yyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_yyyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyzzz_0[i] = 3.0 * g_xxxyyy_0_yyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_yyyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyzzzz_0[i] = 4.0 * g_xxxyyy_0_yyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_yyzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yzzzzz_0[i] = 5.0 * g_xxxyyy_0_yzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_yzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_zzzzzz_0[i] = 6.0 * g_xxxyyy_0_zzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 336-364 components of targeted buffer : KSI

    auto g_xxxyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 336);

    auto g_xxxyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 337);

    auto g_xxxyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 338);

    auto g_xxxyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 339);

    auto g_xxxyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 340);

    auto g_xxxyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 341);

    auto g_xxxyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 342);

    auto g_xxxyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 343);

    auto g_xxxyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 344);

    auto g_xxxyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 345);

    auto g_xxxyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 346);

    auto g_xxxyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 347);

    auto g_xxxyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 348);

    auto g_xxxyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 349);

    auto g_xxxyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 350);

    auto g_xxxyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 351);

    auto g_xxxyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 352);

    auto g_xxxyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 353);

    auto g_xxxyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 354);

    auto g_xxxyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 355);

    auto g_xxxyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 356);

    auto g_xxxyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 357);

    auto g_xxxyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 358);

    auto g_xxxyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 359);

    auto g_xxxyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 360);

    auto g_xxxyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 361);

    auto g_xxxyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 362);

    auto g_xxxyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 363);

    #pragma omp simd aligned(g_xxxyy_0_xxxxxy_0, g_xxxyy_0_xxxxxy_1, g_xxxyy_0_xxxxyy_0, g_xxxyy_0_xxxxyy_1, g_xxxyy_0_xxxyyy_0, g_xxxyy_0_xxxyyy_1, g_xxxyy_0_xxyyyy_0, g_xxxyy_0_xxyyyy_1, g_xxxyy_0_xyyyyy_0, g_xxxyy_0_xyyyyy_1, g_xxxyyz_0_xxxxxy_1, g_xxxyyz_0_xxxxyy_1, g_xxxyyz_0_xxxyyy_1, g_xxxyyz_0_xxyyyy_1, g_xxxyyz_0_xyyyyy_1, g_xxxyyzz_0_xxxxxx_0, g_xxxyyzz_0_xxxxxy_0, g_xxxyyzz_0_xxxxxz_0, g_xxxyyzz_0_xxxxyy_0, g_xxxyyzz_0_xxxxyz_0, g_xxxyyzz_0_xxxxzz_0, g_xxxyyzz_0_xxxyyy_0, g_xxxyyzz_0_xxxyyz_0, g_xxxyyzz_0_xxxyzz_0, g_xxxyyzz_0_xxxzzz_0, g_xxxyyzz_0_xxyyyy_0, g_xxxyyzz_0_xxyyyz_0, g_xxxyyzz_0_xxyyzz_0, g_xxxyyzz_0_xxyzzz_0, g_xxxyyzz_0_xxzzzz_0, g_xxxyyzz_0_xyyyyy_0, g_xxxyyzz_0_xyyyyz_0, g_xxxyyzz_0_xyyyzz_0, g_xxxyyzz_0_xyyzzz_0, g_xxxyyzz_0_xyzzzz_0, g_xxxyyzz_0_xzzzzz_0, g_xxxyyzz_0_yyyyyy_0, g_xxxyyzz_0_yyyyyz_0, g_xxxyyzz_0_yyyyzz_0, g_xxxyyzz_0_yyyzzz_0, g_xxxyyzz_0_yyzzzz_0, g_xxxyyzz_0_yzzzzz_0, g_xxxyyzz_0_zzzzzz_0, g_xxxyzz_0_xxxxxx_1, g_xxxyzz_0_xxxxxz_1, g_xxxyzz_0_xxxxzz_1, g_xxxyzz_0_xxxzzz_1, g_xxxyzz_0_xxzzzz_1, g_xxxyzz_0_xzzzzz_1, g_xxxzz_0_xxxxxx_0, g_xxxzz_0_xxxxxx_1, g_xxxzz_0_xxxxxz_0, g_xxxzz_0_xxxxxz_1, g_xxxzz_0_xxxxzz_0, g_xxxzz_0_xxxxzz_1, g_xxxzz_0_xxxzzz_0, g_xxxzz_0_xxxzzz_1, g_xxxzz_0_xxzzzz_0, g_xxxzz_0_xxzzzz_1, g_xxxzz_0_xzzzzz_0, g_xxxzz_0_xzzzzz_1, g_xxyyzz_0_xxxxyz_1, g_xxyyzz_0_xxxyyz_1, g_xxyyzz_0_xxxyz_1, g_xxyyzz_0_xxxyzz_1, g_xxyyzz_0_xxyyyz_1, g_xxyyzz_0_xxyyz_1, g_xxyyzz_0_xxyyzz_1, g_xxyyzz_0_xxyzz_1, g_xxyyzz_0_xxyzzz_1, g_xxyyzz_0_xyyyyz_1, g_xxyyzz_0_xyyyz_1, g_xxyyzz_0_xyyyzz_1, g_xxyyzz_0_xyyzz_1, g_xxyyzz_0_xyyzzz_1, g_xxyyzz_0_xyzzz_1, g_xxyyzz_0_xyzzzz_1, g_xxyyzz_0_yyyyyy_1, g_xxyyzz_0_yyyyyz_1, g_xxyyzz_0_yyyyz_1, g_xxyyzz_0_yyyyzz_1, g_xxyyzz_0_yyyzz_1, g_xxyyzz_0_yyyzzz_1, g_xxyyzz_0_yyzzz_1, g_xxyyzz_0_yyzzzz_1, g_xxyyzz_0_yzzzz_1, g_xxyyzz_0_yzzzzz_1, g_xxyyzz_0_zzzzzz_1, g_xyyzz_0_xxxxyz_0, g_xyyzz_0_xxxxyz_1, g_xyyzz_0_xxxyyz_0, g_xyyzz_0_xxxyyz_1, g_xyyzz_0_xxxyzz_0, g_xyyzz_0_xxxyzz_1, g_xyyzz_0_xxyyyz_0, g_xyyzz_0_xxyyyz_1, g_xyyzz_0_xxyyzz_0, g_xyyzz_0_xxyyzz_1, g_xyyzz_0_xxyzzz_0, g_xyyzz_0_xxyzzz_1, g_xyyzz_0_xyyyyz_0, g_xyyzz_0_xyyyyz_1, g_xyyzz_0_xyyyzz_0, g_xyyzz_0_xyyyzz_1, g_xyyzz_0_xyyzzz_0, g_xyyzz_0_xyyzzz_1, g_xyyzz_0_xyzzzz_0, g_xyyzz_0_xyzzzz_1, g_xyyzz_0_yyyyyy_0, g_xyyzz_0_yyyyyy_1, g_xyyzz_0_yyyyyz_0, g_xyyzz_0_yyyyyz_1, g_xyyzz_0_yyyyzz_0, g_xyyzz_0_yyyyzz_1, g_xyyzz_0_yyyzzz_0, g_xyyzz_0_yyyzzz_1, g_xyyzz_0_yyzzzz_0, g_xyyzz_0_yyzzzz_1, g_xyyzz_0_yzzzzz_0, g_xyyzz_0_yzzzzz_1, g_xyyzz_0_zzzzzz_0, g_xyyzz_0_zzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyzz_0_xxxxxx_0[i] = g_xxxzz_0_xxxxxx_0[i] * fbe_0 - g_xxxzz_0_xxxxxx_1[i] * fz_be_0 + g_xxxyzz_0_xxxxxx_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxxxy_0[i] = g_xxxyy_0_xxxxxy_0[i] * fbe_0 - g_xxxyy_0_xxxxxy_1[i] * fz_be_0 + g_xxxyyz_0_xxxxxy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxxxz_0[i] = g_xxxzz_0_xxxxxz_0[i] * fbe_0 - g_xxxzz_0_xxxxxz_1[i] * fz_be_0 + g_xxxyzz_0_xxxxxz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxxyy_0[i] = g_xxxyy_0_xxxxyy_0[i] * fbe_0 - g_xxxyy_0_xxxxyy_1[i] * fz_be_0 + g_xxxyyz_0_xxxxyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxxyz_0[i] = 2.0 * g_xyyzz_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxyyzz_0_xxxyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxxzz_0[i] = g_xxxzz_0_xxxxzz_0[i] * fbe_0 - g_xxxzz_0_xxxxzz_1[i] * fz_be_0 + g_xxxyzz_0_xxxxzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxyyy_0[i] = g_xxxyy_0_xxxyyy_0[i] * fbe_0 - g_xxxyy_0_xxxyyy_1[i] * fz_be_0 + g_xxxyyz_0_xxxyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxyyz_0[i] = 2.0 * g_xyyzz_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_0_xxyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxyzz_0[i] = 2.0 * g_xyyzz_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_0_xxyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxzzz_0[i] = g_xxxzz_0_xxxzzz_0[i] * fbe_0 - g_xxxzz_0_xxxzzz_1[i] * fz_be_0 + g_xxxyzz_0_xxxzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxyyyy_0[i] = g_xxxyy_0_xxyyyy_0[i] * fbe_0 - g_xxxyy_0_xxyyyy_1[i] * fz_be_0 + g_xxxyyz_0_xxyyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxyyyz_0[i] = 2.0 * g_xyyzz_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxyyzz_0[i] = 2.0 * g_xyyzz_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxyzzz_0[i] = 2.0 * g_xyyzz_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxzzzz_0[i] = g_xxxzz_0_xxzzzz_0[i] * fbe_0 - g_xxxzz_0_xxzzzz_1[i] * fz_be_0 + g_xxxyzz_0_xxzzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xyyyyy_0[i] = g_xxxyy_0_xyyyyy_0[i] * fbe_0 - g_xxxyy_0_xyyyyy_1[i] * fz_be_0 + g_xxxyyz_0_xyyyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xyyyyz_0[i] = 2.0 * g_xyyzz_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyyyz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyyyzz_0[i] = 2.0 * g_xyyzz_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyyzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyyzzz_0[i] = 2.0 * g_xyyzz_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyzzzz_0[i] = 2.0 * g_xyyzz_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yzzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xzzzzz_0[i] = g_xxxzz_0_xzzzzz_0[i] * fbe_0 - g_xxxzz_0_xzzzzz_1[i] * fz_be_0 + g_xxxyzz_0_xzzzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_yyyyyy_0[i] = 2.0 * g_xyyzz_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyyy_1[i] * fz_be_0 + g_xxyyzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyyyz_0[i] = 2.0 * g_xyyzz_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyyz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyyzz_0[i] = 2.0 * g_xyyzz_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyzzz_0[i] = 2.0 * g_xyyzz_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyzzzz_0[i] = 2.0 * g_xyyzz_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yzzzzz_0[i] = 2.0 * g_xyyzz_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_zzzzzz_0[i] = 2.0 * g_xyyzz_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_zzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 364-392 components of targeted buffer : KSI

    auto g_xxxyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 364);

    auto g_xxxyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 365);

    auto g_xxxyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 366);

    auto g_xxxyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 367);

    auto g_xxxyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 368);

    auto g_xxxyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 369);

    auto g_xxxyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 370);

    auto g_xxxyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 371);

    auto g_xxxyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 372);

    auto g_xxxyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 373);

    auto g_xxxyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 374);

    auto g_xxxyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 375);

    auto g_xxxyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 376);

    auto g_xxxyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 377);

    auto g_xxxyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 378);

    auto g_xxxyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 379);

    auto g_xxxyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 380);

    auto g_xxxyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 381);

    auto g_xxxyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 382);

    auto g_xxxyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 383);

    auto g_xxxyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 384);

    auto g_xxxyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 385);

    auto g_xxxyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 386);

    auto g_xxxyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 387);

    auto g_xxxyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 388);

    auto g_xxxyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 389);

    auto g_xxxyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 390);

    auto g_xxxyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 391);

    #pragma omp simd aligned(g_xxxyzzz_0_xxxxxx_0, g_xxxyzzz_0_xxxxxy_0, g_xxxyzzz_0_xxxxxz_0, g_xxxyzzz_0_xxxxyy_0, g_xxxyzzz_0_xxxxyz_0, g_xxxyzzz_0_xxxxzz_0, g_xxxyzzz_0_xxxyyy_0, g_xxxyzzz_0_xxxyyz_0, g_xxxyzzz_0_xxxyzz_0, g_xxxyzzz_0_xxxzzz_0, g_xxxyzzz_0_xxyyyy_0, g_xxxyzzz_0_xxyyyz_0, g_xxxyzzz_0_xxyyzz_0, g_xxxyzzz_0_xxyzzz_0, g_xxxyzzz_0_xxzzzz_0, g_xxxyzzz_0_xyyyyy_0, g_xxxyzzz_0_xyyyyz_0, g_xxxyzzz_0_xyyyzz_0, g_xxxyzzz_0_xyyzzz_0, g_xxxyzzz_0_xyzzzz_0, g_xxxyzzz_0_xzzzzz_0, g_xxxyzzz_0_yyyyyy_0, g_xxxyzzz_0_yyyyyz_0, g_xxxyzzz_0_yyyyzz_0, g_xxxyzzz_0_yyyzzz_0, g_xxxyzzz_0_yyzzzz_0, g_xxxyzzz_0_yzzzzz_0, g_xxxyzzz_0_zzzzzz_0, g_xxxzzz_0_xxxxx_1, g_xxxzzz_0_xxxxxx_1, g_xxxzzz_0_xxxxxy_1, g_xxxzzz_0_xxxxxz_1, g_xxxzzz_0_xxxxy_1, g_xxxzzz_0_xxxxyy_1, g_xxxzzz_0_xxxxyz_1, g_xxxzzz_0_xxxxz_1, g_xxxzzz_0_xxxxzz_1, g_xxxzzz_0_xxxyy_1, g_xxxzzz_0_xxxyyy_1, g_xxxzzz_0_xxxyyz_1, g_xxxzzz_0_xxxyz_1, g_xxxzzz_0_xxxyzz_1, g_xxxzzz_0_xxxzz_1, g_xxxzzz_0_xxxzzz_1, g_xxxzzz_0_xxyyy_1, g_xxxzzz_0_xxyyyy_1, g_xxxzzz_0_xxyyyz_1, g_xxxzzz_0_xxyyz_1, g_xxxzzz_0_xxyyzz_1, g_xxxzzz_0_xxyzz_1, g_xxxzzz_0_xxyzzz_1, g_xxxzzz_0_xxzzz_1, g_xxxzzz_0_xxzzzz_1, g_xxxzzz_0_xyyyy_1, g_xxxzzz_0_xyyyyy_1, g_xxxzzz_0_xyyyyz_1, g_xxxzzz_0_xyyyz_1, g_xxxzzz_0_xyyyzz_1, g_xxxzzz_0_xyyzz_1, g_xxxzzz_0_xyyzzz_1, g_xxxzzz_0_xyzzz_1, g_xxxzzz_0_xyzzzz_1, g_xxxzzz_0_xzzzz_1, g_xxxzzz_0_xzzzzz_1, g_xxxzzz_0_yyyyy_1, g_xxxzzz_0_yyyyyy_1, g_xxxzzz_0_yyyyyz_1, g_xxxzzz_0_yyyyz_1, g_xxxzzz_0_yyyyzz_1, g_xxxzzz_0_yyyzz_1, g_xxxzzz_0_yyyzzz_1, g_xxxzzz_0_yyzzz_1, g_xxxzzz_0_yyzzzz_1, g_xxxzzz_0_yzzzz_1, g_xxxzzz_0_yzzzzz_1, g_xxxzzz_0_zzzzz_1, g_xxxzzz_0_zzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzzz_0_xxxxxx_0[i] = g_xxxzzz_0_xxxxxx_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxy_0[i] = g_xxxzzz_0_xxxxx_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxz_0[i] = g_xxxzzz_0_xxxxxz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxyy_0[i] = 2.0 * g_xxxzzz_0_xxxxy_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxyz_0[i] = g_xxxzzz_0_xxxxz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxzz_0[i] = g_xxxzzz_0_xxxxzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyyy_0[i] = 3.0 * g_xxxzzz_0_xxxyy_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyyz_0[i] = 2.0 * g_xxxzzz_0_xxxyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyzz_0[i] = g_xxxzzz_0_xxxzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxzzz_0[i] = g_xxxzzz_0_xxxzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyyy_0[i] = 4.0 * g_xxxzzz_0_xxyyy_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyyz_0[i] = 3.0 * g_xxxzzz_0_xxyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyzz_0[i] = 2.0 * g_xxxzzz_0_xxyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyzzz_0[i] = g_xxxzzz_0_xxzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxzzzz_0[i] = g_xxxzzz_0_xxzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyyy_0[i] = 5.0 * g_xxxzzz_0_xyyyy_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyyz_0[i] = 4.0 * g_xxxzzz_0_xyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyzz_0[i] = 3.0 * g_xxxzzz_0_xyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyzzz_0[i] = 2.0 * g_xxxzzz_0_xyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyzzzz_0[i] = g_xxxzzz_0_xzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xzzzzz_0[i] = g_xxxzzz_0_xzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyyy_0[i] = 6.0 * g_xxxzzz_0_yyyyy_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyyz_0[i] = 5.0 * g_xxxzzz_0_yyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyzz_0[i] = 4.0 * g_xxxzzz_0_yyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyzzz_0[i] = 3.0 * g_xxxzzz_0_yyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyzzzz_0[i] = 2.0 * g_xxxzzz_0_yzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yzzzzz_0[i] = g_xxxzzz_0_zzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_zzzzzz_0[i] = g_xxxzzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 392-420 components of targeted buffer : KSI

    auto g_xxxzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 392);

    auto g_xxxzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 393);

    auto g_xxxzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 394);

    auto g_xxxzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 395);

    auto g_xxxzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 396);

    auto g_xxxzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 397);

    auto g_xxxzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 398);

    auto g_xxxzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 399);

    auto g_xxxzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 400);

    auto g_xxxzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 401);

    auto g_xxxzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 402);

    auto g_xxxzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 403);

    auto g_xxxzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 404);

    auto g_xxxzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 405);

    auto g_xxxzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 406);

    auto g_xxxzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 407);

    auto g_xxxzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 408);

    auto g_xxxzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 409);

    auto g_xxxzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 410);

    auto g_xxxzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 411);

    auto g_xxxzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 412);

    auto g_xxxzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 413);

    auto g_xxxzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 414);

    auto g_xxxzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 415);

    auto g_xxxzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 416);

    auto g_xxxzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 417);

    auto g_xxxzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 418);

    auto g_xxxzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 419);

    #pragma omp simd aligned(g_xxxzz_0_xxxxxx_0, g_xxxzz_0_xxxxxx_1, g_xxxzz_0_xxxxxy_0, g_xxxzz_0_xxxxxy_1, g_xxxzz_0_xxxxyy_0, g_xxxzz_0_xxxxyy_1, g_xxxzz_0_xxxyyy_0, g_xxxzz_0_xxxyyy_1, g_xxxzz_0_xxyyyy_0, g_xxxzz_0_xxyyyy_1, g_xxxzz_0_xyyyyy_0, g_xxxzz_0_xyyyyy_1, g_xxxzzz_0_xxxxxx_1, g_xxxzzz_0_xxxxxy_1, g_xxxzzz_0_xxxxyy_1, g_xxxzzz_0_xxxyyy_1, g_xxxzzz_0_xxyyyy_1, g_xxxzzz_0_xyyyyy_1, g_xxxzzzz_0_xxxxxx_0, g_xxxzzzz_0_xxxxxy_0, g_xxxzzzz_0_xxxxxz_0, g_xxxzzzz_0_xxxxyy_0, g_xxxzzzz_0_xxxxyz_0, g_xxxzzzz_0_xxxxzz_0, g_xxxzzzz_0_xxxyyy_0, g_xxxzzzz_0_xxxyyz_0, g_xxxzzzz_0_xxxyzz_0, g_xxxzzzz_0_xxxzzz_0, g_xxxzzzz_0_xxyyyy_0, g_xxxzzzz_0_xxyyyz_0, g_xxxzzzz_0_xxyyzz_0, g_xxxzzzz_0_xxyzzz_0, g_xxxzzzz_0_xxzzzz_0, g_xxxzzzz_0_xyyyyy_0, g_xxxzzzz_0_xyyyyz_0, g_xxxzzzz_0_xyyyzz_0, g_xxxzzzz_0_xyyzzz_0, g_xxxzzzz_0_xyzzzz_0, g_xxxzzzz_0_xzzzzz_0, g_xxxzzzz_0_yyyyyy_0, g_xxxzzzz_0_yyyyyz_0, g_xxxzzzz_0_yyyyzz_0, g_xxxzzzz_0_yyyzzz_0, g_xxxzzzz_0_yyzzzz_0, g_xxxzzzz_0_yzzzzz_0, g_xxxzzzz_0_zzzzzz_0, g_xxzzzz_0_xxxxxz_1, g_xxzzzz_0_xxxxyz_1, g_xxzzzz_0_xxxxz_1, g_xxzzzz_0_xxxxzz_1, g_xxzzzz_0_xxxyyz_1, g_xxzzzz_0_xxxyz_1, g_xxzzzz_0_xxxyzz_1, g_xxzzzz_0_xxxzz_1, g_xxzzzz_0_xxxzzz_1, g_xxzzzz_0_xxyyyz_1, g_xxzzzz_0_xxyyz_1, g_xxzzzz_0_xxyyzz_1, g_xxzzzz_0_xxyzz_1, g_xxzzzz_0_xxyzzz_1, g_xxzzzz_0_xxzzz_1, g_xxzzzz_0_xxzzzz_1, g_xxzzzz_0_xyyyyz_1, g_xxzzzz_0_xyyyz_1, g_xxzzzz_0_xyyyzz_1, g_xxzzzz_0_xyyzz_1, g_xxzzzz_0_xyyzzz_1, g_xxzzzz_0_xyzzz_1, g_xxzzzz_0_xyzzzz_1, g_xxzzzz_0_xzzzz_1, g_xxzzzz_0_xzzzzz_1, g_xxzzzz_0_yyyyyy_1, g_xxzzzz_0_yyyyyz_1, g_xxzzzz_0_yyyyz_1, g_xxzzzz_0_yyyyzz_1, g_xxzzzz_0_yyyzz_1, g_xxzzzz_0_yyyzzz_1, g_xxzzzz_0_yyzzz_1, g_xxzzzz_0_yyzzzz_1, g_xxzzzz_0_yzzzz_1, g_xxzzzz_0_yzzzzz_1, g_xxzzzz_0_zzzzz_1, g_xxzzzz_0_zzzzzz_1, g_xzzzz_0_xxxxxz_0, g_xzzzz_0_xxxxxz_1, g_xzzzz_0_xxxxyz_0, g_xzzzz_0_xxxxyz_1, g_xzzzz_0_xxxxzz_0, g_xzzzz_0_xxxxzz_1, g_xzzzz_0_xxxyyz_0, g_xzzzz_0_xxxyyz_1, g_xzzzz_0_xxxyzz_0, g_xzzzz_0_xxxyzz_1, g_xzzzz_0_xxxzzz_0, g_xzzzz_0_xxxzzz_1, g_xzzzz_0_xxyyyz_0, g_xzzzz_0_xxyyyz_1, g_xzzzz_0_xxyyzz_0, g_xzzzz_0_xxyyzz_1, g_xzzzz_0_xxyzzz_0, g_xzzzz_0_xxyzzz_1, g_xzzzz_0_xxzzzz_0, g_xzzzz_0_xxzzzz_1, g_xzzzz_0_xyyyyz_0, g_xzzzz_0_xyyyyz_1, g_xzzzz_0_xyyyzz_0, g_xzzzz_0_xyyyzz_1, g_xzzzz_0_xyyzzz_0, g_xzzzz_0_xyyzzz_1, g_xzzzz_0_xyzzzz_0, g_xzzzz_0_xyzzzz_1, g_xzzzz_0_xzzzzz_0, g_xzzzz_0_xzzzzz_1, g_xzzzz_0_yyyyyy_0, g_xzzzz_0_yyyyyy_1, g_xzzzz_0_yyyyyz_0, g_xzzzz_0_yyyyyz_1, g_xzzzz_0_yyyyzz_0, g_xzzzz_0_yyyyzz_1, g_xzzzz_0_yyyzzz_0, g_xzzzz_0_yyyzzz_1, g_xzzzz_0_yyzzzz_0, g_xzzzz_0_yyzzzz_1, g_xzzzz_0_yzzzzz_0, g_xzzzz_0_yzzzzz_1, g_xzzzz_0_zzzzzz_0, g_xzzzz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzzz_0_xxxxxx_0[i] = 3.0 * g_xxxzz_0_xxxxxx_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxxx_1[i] * fz_be_0 + g_xxxzzz_0_xxxxxx_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxxy_0[i] = 3.0 * g_xxxzz_0_xxxxxy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxxy_1[i] * fz_be_0 + g_xxxzzz_0_xxxxxy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxxz_0[i] = 2.0 * g_xzzzz_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxzzzz_0_xxxxz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxyy_0[i] = 3.0 * g_xxxzz_0_xxxxyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxyy_1[i] * fz_be_0 + g_xxxzzz_0_xxxxyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxyz_0[i] = 2.0 * g_xzzzz_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_0_xxxyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxzz_0[i] = 2.0 * g_xzzzz_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_0_xxxzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxyyy_0[i] = 3.0 * g_xxxzz_0_xxxyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxyyy_1[i] * fz_be_0 + g_xxxzzz_0_xxxyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxyyz_0[i] = 2.0 * g_xzzzz_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxyzz_0[i] = 2.0 * g_xzzzz_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxzzz_0[i] = 2.0 * g_xzzzz_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyyyy_0[i] = 3.0 * g_xxxzz_0_xxyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxyyyy_1[i] * fz_be_0 + g_xxxzzz_0_xxyyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxyyyz_0[i] = 2.0 * g_xzzzz_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyyzz_0[i] = 2.0 * g_xzzzz_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyzzz_0[i] = 2.0 * g_xzzzz_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxzzzz_0[i] = 2.0 * g_xzzzz_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyyyy_0[i] = 3.0 * g_xxxzz_0_xyyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xyyyyy_1[i] * fz_be_0 + g_xxxzzz_0_xyyyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xyyyyz_0[i] = 2.0 * g_xzzzz_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyyyz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyyzz_0[i] = 2.0 * g_xzzzz_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyyzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyzzz_0[i] = 2.0 * g_xzzzz_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyzzzz_0[i] = 2.0 * g_xzzzz_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xzzzzz_0[i] = 2.0 * g_xzzzz_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_zzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyyy_0[i] = 2.0 * g_xzzzz_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyyy_1[i] * fz_be_0 + g_xxzzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyyz_0[i] = 2.0 * g_xzzzz_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyyz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyzz_0[i] = 2.0 * g_xzzzz_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyzzz_0[i] = 2.0 * g_xzzzz_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyzzzz_0[i] = 2.0 * g_xzzzz_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yzzzzz_0[i] = 2.0 * g_xzzzz_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_zzzzzz_0[i] = 2.0 * g_xzzzz_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_zzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 420-448 components of targeted buffer : KSI

    auto g_xxyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 420);

    auto g_xxyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 421);

    auto g_xxyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 422);

    auto g_xxyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 423);

    auto g_xxyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 424);

    auto g_xxyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 425);

    auto g_xxyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 426);

    auto g_xxyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 427);

    auto g_xxyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 428);

    auto g_xxyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 429);

    auto g_xxyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 430);

    auto g_xxyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 431);

    auto g_xxyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 432);

    auto g_xxyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 433);

    auto g_xxyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 434);

    auto g_xxyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 435);

    auto g_xxyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 436);

    auto g_xxyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 437);

    auto g_xxyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 438);

    auto g_xxyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 439);

    auto g_xxyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 440);

    auto g_xxyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 441);

    auto g_xxyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 442);

    auto g_xxyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 443);

    auto g_xxyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 444);

    auto g_xxyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 445);

    auto g_xxyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 446);

    auto g_xxyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 447);

    #pragma omp simd aligned(g_xxyyy_0_xxxxxx_0, g_xxyyy_0_xxxxxx_1, g_xxyyy_0_xxxxxz_0, g_xxyyy_0_xxxxxz_1, g_xxyyy_0_xxxxzz_0, g_xxyyy_0_xxxxzz_1, g_xxyyy_0_xxxzzz_0, g_xxyyy_0_xxxzzz_1, g_xxyyy_0_xxzzzz_0, g_xxyyy_0_xxzzzz_1, g_xxyyy_0_xzzzzz_0, g_xxyyy_0_xzzzzz_1, g_xxyyyy_0_xxxxxx_1, g_xxyyyy_0_xxxxxz_1, g_xxyyyy_0_xxxxzz_1, g_xxyyyy_0_xxxzzz_1, g_xxyyyy_0_xxzzzz_1, g_xxyyyy_0_xzzzzz_1, g_xxyyyyy_0_xxxxxx_0, g_xxyyyyy_0_xxxxxy_0, g_xxyyyyy_0_xxxxxz_0, g_xxyyyyy_0_xxxxyy_0, g_xxyyyyy_0_xxxxyz_0, g_xxyyyyy_0_xxxxzz_0, g_xxyyyyy_0_xxxyyy_0, g_xxyyyyy_0_xxxyyz_0, g_xxyyyyy_0_xxxyzz_0, g_xxyyyyy_0_xxxzzz_0, g_xxyyyyy_0_xxyyyy_0, g_xxyyyyy_0_xxyyyz_0, g_xxyyyyy_0_xxyyzz_0, g_xxyyyyy_0_xxyzzz_0, g_xxyyyyy_0_xxzzzz_0, g_xxyyyyy_0_xyyyyy_0, g_xxyyyyy_0_xyyyyz_0, g_xxyyyyy_0_xyyyzz_0, g_xxyyyyy_0_xyyzzz_0, g_xxyyyyy_0_xyzzzz_0, g_xxyyyyy_0_xzzzzz_0, g_xxyyyyy_0_yyyyyy_0, g_xxyyyyy_0_yyyyyz_0, g_xxyyyyy_0_yyyyzz_0, g_xxyyyyy_0_yyyzzz_0, g_xxyyyyy_0_yyzzzz_0, g_xxyyyyy_0_yzzzzz_0, g_xxyyyyy_0_zzzzzz_0, g_xyyyyy_0_xxxxxy_1, g_xyyyyy_0_xxxxy_1, g_xyyyyy_0_xxxxyy_1, g_xyyyyy_0_xxxxyz_1, g_xyyyyy_0_xxxyy_1, g_xyyyyy_0_xxxyyy_1, g_xyyyyy_0_xxxyyz_1, g_xyyyyy_0_xxxyz_1, g_xyyyyy_0_xxxyzz_1, g_xyyyyy_0_xxyyy_1, g_xyyyyy_0_xxyyyy_1, g_xyyyyy_0_xxyyyz_1, g_xyyyyy_0_xxyyz_1, g_xyyyyy_0_xxyyzz_1, g_xyyyyy_0_xxyzz_1, g_xyyyyy_0_xxyzzz_1, g_xyyyyy_0_xyyyy_1, g_xyyyyy_0_xyyyyy_1, g_xyyyyy_0_xyyyyz_1, g_xyyyyy_0_xyyyz_1, g_xyyyyy_0_xyyyzz_1, g_xyyyyy_0_xyyzz_1, g_xyyyyy_0_xyyzzz_1, g_xyyyyy_0_xyzzz_1, g_xyyyyy_0_xyzzzz_1, g_xyyyyy_0_yyyyy_1, g_xyyyyy_0_yyyyyy_1, g_xyyyyy_0_yyyyyz_1, g_xyyyyy_0_yyyyz_1, g_xyyyyy_0_yyyyzz_1, g_xyyyyy_0_yyyzz_1, g_xyyyyy_0_yyyzzz_1, g_xyyyyy_0_yyzzz_1, g_xyyyyy_0_yyzzzz_1, g_xyyyyy_0_yzzzz_1, g_xyyyyy_0_yzzzzz_1, g_xyyyyy_0_zzzzzz_1, g_yyyyy_0_xxxxxy_0, g_yyyyy_0_xxxxxy_1, g_yyyyy_0_xxxxyy_0, g_yyyyy_0_xxxxyy_1, g_yyyyy_0_xxxxyz_0, g_yyyyy_0_xxxxyz_1, g_yyyyy_0_xxxyyy_0, g_yyyyy_0_xxxyyy_1, g_yyyyy_0_xxxyyz_0, g_yyyyy_0_xxxyyz_1, g_yyyyy_0_xxxyzz_0, g_yyyyy_0_xxxyzz_1, g_yyyyy_0_xxyyyy_0, g_yyyyy_0_xxyyyy_1, g_yyyyy_0_xxyyyz_0, g_yyyyy_0_xxyyyz_1, g_yyyyy_0_xxyyzz_0, g_yyyyy_0_xxyyzz_1, g_yyyyy_0_xxyzzz_0, g_yyyyy_0_xxyzzz_1, g_yyyyy_0_xyyyyy_0, g_yyyyy_0_xyyyyy_1, g_yyyyy_0_xyyyyz_0, g_yyyyy_0_xyyyyz_1, g_yyyyy_0_xyyyzz_0, g_yyyyy_0_xyyyzz_1, g_yyyyy_0_xyyzzz_0, g_yyyyy_0_xyyzzz_1, g_yyyyy_0_xyzzzz_0, g_yyyyy_0_xyzzzz_1, g_yyyyy_0_yyyyyy_0, g_yyyyy_0_yyyyyy_1, g_yyyyy_0_yyyyyz_0, g_yyyyy_0_yyyyyz_1, g_yyyyy_0_yyyyzz_0, g_yyyyy_0_yyyyzz_1, g_yyyyy_0_yyyzzz_0, g_yyyyy_0_yyyzzz_1, g_yyyyy_0_yyzzzz_0, g_yyyyy_0_yyzzzz_1, g_yyyyy_0_yzzzzz_0, g_yyyyy_0_yzzzzz_1, g_yyyyy_0_zzzzzz_0, g_yyyyy_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyyy_0_xxxxxx_0[i] = 4.0 * g_xxyyy_0_xxxxxx_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxxx_1[i] * fz_be_0 + g_xxyyyy_0_xxxxxx_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxxxy_0[i] = g_yyyyy_0_xxxxxy_0[i] * fbe_0 - g_yyyyy_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xyyyyy_0_xxxxy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxxy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxxz_0[i] = 4.0 * g_xxyyy_0_xxxxxz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxxz_1[i] * fz_be_0 + g_xxyyyy_0_xxxxxz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxxyy_0[i] = g_yyyyy_0_xxxxyy_0[i] * fbe_0 - g_yyyyy_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xyyyyy_0_xxxyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxyz_0[i] = g_yyyyy_0_xxxxyz_0[i] * fbe_0 - g_yyyyy_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyyyyy_0_xxxyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxzz_0[i] = 4.0 * g_xxyyy_0_xxxxzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxzz_1[i] * fz_be_0 + g_xxyyyy_0_xxxxzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxyyy_0[i] = g_yyyyy_0_xxxyyy_0[i] * fbe_0 - g_yyyyy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxyyz_0[i] = g_yyyyy_0_xxxyyz_0[i] * fbe_0 - g_yyyyy_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxyzz_0[i] = g_yyyyy_0_xxxyzz_0[i] * fbe_0 - g_yyyyy_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxzzz_0[i] = 4.0 * g_xxyyy_0_xxxzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxzzz_1[i] * fz_be_0 + g_xxyyyy_0_xxxzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxyyyy_0[i] = g_yyyyy_0_xxyyyy_0[i] * fbe_0 - g_yyyyy_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyyyz_0[i] = g_yyyyy_0_xxyyyz_0[i] * fbe_0 - g_yyyyy_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyyzz_0[i] = g_yyyyy_0_xxyyzz_0[i] * fbe_0 - g_yyyyy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyzzz_0[i] = g_yyyyy_0_xxyzzz_0[i] * fbe_0 - g_yyyyy_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxzzzz_0[i] = 4.0 * g_xxyyy_0_xxzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxzzzz_1[i] * fz_be_0 + g_xxyyyy_0_xxzzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xyyyyy_0[i] = g_yyyyy_0_xyyyyy_0[i] * fbe_0 - g_yyyyy_0_xyyyyy_1[i] * fz_be_0 + g_xyyyyy_0_yyyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyyyz_0[i] = g_yyyyy_0_xyyyyz_0[i] * fbe_0 - g_yyyyy_0_xyyyyz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyyzz_0[i] = g_yyyyy_0_xyyyzz_0[i] * fbe_0 - g_yyyyy_0_xyyyzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyzzz_0[i] = g_yyyyy_0_xyyzzz_0[i] * fbe_0 - g_yyyyy_0_xyyzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyzzzz_0[i] = g_yyyyy_0_xyzzzz_0[i] * fbe_0 - g_yyyyy_0_xyzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yzzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xzzzzz_0[i] = 4.0 * g_xxyyy_0_xzzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_xzzzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_yyyyyy_0[i] = g_yyyyy_0_yyyyyy_0[i] * fbe_0 - g_yyyyy_0_yyyyyy_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyyyz_0[i] = g_yyyyy_0_yyyyyz_0[i] * fbe_0 - g_yyyyy_0_yyyyyz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyyzz_0[i] = g_yyyyy_0_yyyyzz_0[i] * fbe_0 - g_yyyyy_0_yyyyzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyzzz_0[i] = g_yyyyy_0_yyyzzz_0[i] * fbe_0 - g_yyyyy_0_yyyzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyzzzz_0[i] = g_yyyyy_0_yyzzzz_0[i] * fbe_0 - g_yyyyy_0_yyzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yzzzzz_0[i] = g_yyyyy_0_yzzzzz_0[i] * fbe_0 - g_yyyyy_0_yzzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yzzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_zzzzzz_0[i] = g_yyyyy_0_zzzzzz_0[i] * fbe_0 - g_yyyyy_0_zzzzzz_1[i] * fz_be_0 + g_xyyyyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 448-476 components of targeted buffer : KSI

    auto g_xxyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 448);

    auto g_xxyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 449);

    auto g_xxyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 450);

    auto g_xxyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 451);

    auto g_xxyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 452);

    auto g_xxyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 453);

    auto g_xxyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 454);

    auto g_xxyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 455);

    auto g_xxyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 456);

    auto g_xxyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 457);

    auto g_xxyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 458);

    auto g_xxyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 459);

    auto g_xxyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 460);

    auto g_xxyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 461);

    auto g_xxyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 462);

    auto g_xxyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 463);

    auto g_xxyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 464);

    auto g_xxyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 465);

    auto g_xxyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 466);

    auto g_xxyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 467);

    auto g_xxyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 468);

    auto g_xxyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 469);

    auto g_xxyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 470);

    auto g_xxyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 471);

    auto g_xxyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 472);

    auto g_xxyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 473);

    auto g_xxyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 474);

    auto g_xxyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 475);

    #pragma omp simd aligned(g_xxyyyy_0_xxxxx_1, g_xxyyyy_0_xxxxxx_1, g_xxyyyy_0_xxxxxy_1, g_xxyyyy_0_xxxxxz_1, g_xxyyyy_0_xxxxy_1, g_xxyyyy_0_xxxxyy_1, g_xxyyyy_0_xxxxyz_1, g_xxyyyy_0_xxxxz_1, g_xxyyyy_0_xxxxzz_1, g_xxyyyy_0_xxxyy_1, g_xxyyyy_0_xxxyyy_1, g_xxyyyy_0_xxxyyz_1, g_xxyyyy_0_xxxyz_1, g_xxyyyy_0_xxxyzz_1, g_xxyyyy_0_xxxzz_1, g_xxyyyy_0_xxxzzz_1, g_xxyyyy_0_xxyyy_1, g_xxyyyy_0_xxyyyy_1, g_xxyyyy_0_xxyyyz_1, g_xxyyyy_0_xxyyz_1, g_xxyyyy_0_xxyyzz_1, g_xxyyyy_0_xxyzz_1, g_xxyyyy_0_xxyzzz_1, g_xxyyyy_0_xxzzz_1, g_xxyyyy_0_xxzzzz_1, g_xxyyyy_0_xyyyy_1, g_xxyyyy_0_xyyyyy_1, g_xxyyyy_0_xyyyyz_1, g_xxyyyy_0_xyyyz_1, g_xxyyyy_0_xyyyzz_1, g_xxyyyy_0_xyyzz_1, g_xxyyyy_0_xyyzzz_1, g_xxyyyy_0_xyzzz_1, g_xxyyyy_0_xyzzzz_1, g_xxyyyy_0_xzzzz_1, g_xxyyyy_0_xzzzzz_1, g_xxyyyy_0_yyyyy_1, g_xxyyyy_0_yyyyyy_1, g_xxyyyy_0_yyyyyz_1, g_xxyyyy_0_yyyyz_1, g_xxyyyy_0_yyyyzz_1, g_xxyyyy_0_yyyzz_1, g_xxyyyy_0_yyyzzz_1, g_xxyyyy_0_yyzzz_1, g_xxyyyy_0_yyzzzz_1, g_xxyyyy_0_yzzzz_1, g_xxyyyy_0_yzzzzz_1, g_xxyyyy_0_zzzzz_1, g_xxyyyy_0_zzzzzz_1, g_xxyyyyz_0_xxxxxx_0, g_xxyyyyz_0_xxxxxy_0, g_xxyyyyz_0_xxxxxz_0, g_xxyyyyz_0_xxxxyy_0, g_xxyyyyz_0_xxxxyz_0, g_xxyyyyz_0_xxxxzz_0, g_xxyyyyz_0_xxxyyy_0, g_xxyyyyz_0_xxxyyz_0, g_xxyyyyz_0_xxxyzz_0, g_xxyyyyz_0_xxxzzz_0, g_xxyyyyz_0_xxyyyy_0, g_xxyyyyz_0_xxyyyz_0, g_xxyyyyz_0_xxyyzz_0, g_xxyyyyz_0_xxyzzz_0, g_xxyyyyz_0_xxzzzz_0, g_xxyyyyz_0_xyyyyy_0, g_xxyyyyz_0_xyyyyz_0, g_xxyyyyz_0_xyyyzz_0, g_xxyyyyz_0_xyyzzz_0, g_xxyyyyz_0_xyzzzz_0, g_xxyyyyz_0_xzzzzz_0, g_xxyyyyz_0_yyyyyy_0, g_xxyyyyz_0_yyyyyz_0, g_xxyyyyz_0_yyyyzz_0, g_xxyyyyz_0_yyyzzz_0, g_xxyyyyz_0_yyzzzz_0, g_xxyyyyz_0_yzzzzz_0, g_xxyyyyz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyyz_0_xxxxxx_0[i] = g_xxyyyy_0_xxxxxx_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxy_0[i] = g_xxyyyy_0_xxxxxy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxz_0[i] = g_xxyyyy_0_xxxxx_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxyy_0[i] = g_xxyyyy_0_xxxxyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxyz_0[i] = g_xxyyyy_0_xxxxy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxzz_0[i] = 2.0 * g_xxyyyy_0_xxxxz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyyy_0[i] = g_xxyyyy_0_xxxyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyyz_0[i] = g_xxyyyy_0_xxxyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyzz_0[i] = 2.0 * g_xxyyyy_0_xxxyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxzzz_0[i] = 3.0 * g_xxyyyy_0_xxxzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyyy_0[i] = g_xxyyyy_0_xxyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyyz_0[i] = g_xxyyyy_0_xxyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyzz_0[i] = 2.0 * g_xxyyyy_0_xxyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyzzz_0[i] = 3.0 * g_xxyyyy_0_xxyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxzzzz_0[i] = 4.0 * g_xxyyyy_0_xxzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyyy_0[i] = g_xxyyyy_0_xyyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyyz_0[i] = g_xxyyyy_0_xyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyzz_0[i] = 2.0 * g_xxyyyy_0_xyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyzzz_0[i] = 3.0 * g_xxyyyy_0_xyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyzzzz_0[i] = 4.0 * g_xxyyyy_0_xyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xzzzzz_0[i] = 5.0 * g_xxyyyy_0_xzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyyy_0[i] = g_xxyyyy_0_yyyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyyz_0[i] = g_xxyyyy_0_yyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_yyyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyzz_0[i] = 2.0 * g_xxyyyy_0_yyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_yyyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyzzz_0[i] = 3.0 * g_xxyyyy_0_yyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_yyyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyzzzz_0[i] = 4.0 * g_xxyyyy_0_yyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_yyzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yzzzzz_0[i] = 5.0 * g_xxyyyy_0_yzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_yzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_zzzzzz_0[i] = 6.0 * g_xxyyyy_0_zzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 476-504 components of targeted buffer : KSI

    auto g_xxyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 476);

    auto g_xxyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 477);

    auto g_xxyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 478);

    auto g_xxyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 479);

    auto g_xxyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 480);

    auto g_xxyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 481);

    auto g_xxyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 482);

    auto g_xxyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 483);

    auto g_xxyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 484);

    auto g_xxyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 485);

    auto g_xxyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 486);

    auto g_xxyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 487);

    auto g_xxyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 488);

    auto g_xxyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 489);

    auto g_xxyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 490);

    auto g_xxyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 491);

    auto g_xxyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 492);

    auto g_xxyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 493);

    auto g_xxyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 494);

    auto g_xxyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 495);

    auto g_xxyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 496);

    auto g_xxyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 497);

    auto g_xxyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 498);

    auto g_xxyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 499);

    auto g_xxyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 500);

    auto g_xxyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 501);

    auto g_xxyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 502);

    auto g_xxyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 503);

    #pragma omp simd aligned(g_xxyyy_0_xxxxxy_0, g_xxyyy_0_xxxxxy_1, g_xxyyy_0_xxxxyy_0, g_xxyyy_0_xxxxyy_1, g_xxyyy_0_xxxyyy_0, g_xxyyy_0_xxxyyy_1, g_xxyyy_0_xxyyyy_0, g_xxyyy_0_xxyyyy_1, g_xxyyy_0_xyyyyy_0, g_xxyyy_0_xyyyyy_1, g_xxyyyz_0_xxxxxy_1, g_xxyyyz_0_xxxxyy_1, g_xxyyyz_0_xxxyyy_1, g_xxyyyz_0_xxyyyy_1, g_xxyyyz_0_xyyyyy_1, g_xxyyyzz_0_xxxxxx_0, g_xxyyyzz_0_xxxxxy_0, g_xxyyyzz_0_xxxxxz_0, g_xxyyyzz_0_xxxxyy_0, g_xxyyyzz_0_xxxxyz_0, g_xxyyyzz_0_xxxxzz_0, g_xxyyyzz_0_xxxyyy_0, g_xxyyyzz_0_xxxyyz_0, g_xxyyyzz_0_xxxyzz_0, g_xxyyyzz_0_xxxzzz_0, g_xxyyyzz_0_xxyyyy_0, g_xxyyyzz_0_xxyyyz_0, g_xxyyyzz_0_xxyyzz_0, g_xxyyyzz_0_xxyzzz_0, g_xxyyyzz_0_xxzzzz_0, g_xxyyyzz_0_xyyyyy_0, g_xxyyyzz_0_xyyyyz_0, g_xxyyyzz_0_xyyyzz_0, g_xxyyyzz_0_xyyzzz_0, g_xxyyyzz_0_xyzzzz_0, g_xxyyyzz_0_xzzzzz_0, g_xxyyyzz_0_yyyyyy_0, g_xxyyyzz_0_yyyyyz_0, g_xxyyyzz_0_yyyyzz_0, g_xxyyyzz_0_yyyzzz_0, g_xxyyyzz_0_yyzzzz_0, g_xxyyyzz_0_yzzzzz_0, g_xxyyyzz_0_zzzzzz_0, g_xxyyzz_0_xxxxxx_1, g_xxyyzz_0_xxxxxz_1, g_xxyyzz_0_xxxxzz_1, g_xxyyzz_0_xxxzzz_1, g_xxyyzz_0_xxzzzz_1, g_xxyyzz_0_xzzzzz_1, g_xxyzz_0_xxxxxx_0, g_xxyzz_0_xxxxxx_1, g_xxyzz_0_xxxxxz_0, g_xxyzz_0_xxxxxz_1, g_xxyzz_0_xxxxzz_0, g_xxyzz_0_xxxxzz_1, g_xxyzz_0_xxxzzz_0, g_xxyzz_0_xxxzzz_1, g_xxyzz_0_xxzzzz_0, g_xxyzz_0_xxzzzz_1, g_xxyzz_0_xzzzzz_0, g_xxyzz_0_xzzzzz_1, g_xyyyzz_0_xxxxyz_1, g_xyyyzz_0_xxxyyz_1, g_xyyyzz_0_xxxyz_1, g_xyyyzz_0_xxxyzz_1, g_xyyyzz_0_xxyyyz_1, g_xyyyzz_0_xxyyz_1, g_xyyyzz_0_xxyyzz_1, g_xyyyzz_0_xxyzz_1, g_xyyyzz_0_xxyzzz_1, g_xyyyzz_0_xyyyyz_1, g_xyyyzz_0_xyyyz_1, g_xyyyzz_0_xyyyzz_1, g_xyyyzz_0_xyyzz_1, g_xyyyzz_0_xyyzzz_1, g_xyyyzz_0_xyzzz_1, g_xyyyzz_0_xyzzzz_1, g_xyyyzz_0_yyyyyy_1, g_xyyyzz_0_yyyyyz_1, g_xyyyzz_0_yyyyz_1, g_xyyyzz_0_yyyyzz_1, g_xyyyzz_0_yyyzz_1, g_xyyyzz_0_yyyzzz_1, g_xyyyzz_0_yyzzz_1, g_xyyyzz_0_yyzzzz_1, g_xyyyzz_0_yzzzz_1, g_xyyyzz_0_yzzzzz_1, g_xyyyzz_0_zzzzzz_1, g_yyyzz_0_xxxxyz_0, g_yyyzz_0_xxxxyz_1, g_yyyzz_0_xxxyyz_0, g_yyyzz_0_xxxyyz_1, g_yyyzz_0_xxxyzz_0, g_yyyzz_0_xxxyzz_1, g_yyyzz_0_xxyyyz_0, g_yyyzz_0_xxyyyz_1, g_yyyzz_0_xxyyzz_0, g_yyyzz_0_xxyyzz_1, g_yyyzz_0_xxyzzz_0, g_yyyzz_0_xxyzzz_1, g_yyyzz_0_xyyyyz_0, g_yyyzz_0_xyyyyz_1, g_yyyzz_0_xyyyzz_0, g_yyyzz_0_xyyyzz_1, g_yyyzz_0_xyyzzz_0, g_yyyzz_0_xyyzzz_1, g_yyyzz_0_xyzzzz_0, g_yyyzz_0_xyzzzz_1, g_yyyzz_0_yyyyyy_0, g_yyyzz_0_yyyyyy_1, g_yyyzz_0_yyyyyz_0, g_yyyzz_0_yyyyyz_1, g_yyyzz_0_yyyyzz_0, g_yyyzz_0_yyyyzz_1, g_yyyzz_0_yyyzzz_0, g_yyyzz_0_yyyzzz_1, g_yyyzz_0_yyzzzz_0, g_yyyzz_0_yyzzzz_1, g_yyyzz_0_yzzzzz_0, g_yyyzz_0_yzzzzz_1, g_yyyzz_0_zzzzzz_0, g_yyyzz_0_zzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyzz_0_xxxxxx_0[i] = 2.0 * g_xxyzz_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxxx_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxx_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxxxy_0[i] = g_xxyyy_0_xxxxxy_0[i] * fbe_0 - g_xxyyy_0_xxxxxy_1[i] * fz_be_0 + g_xxyyyz_0_xxxxxy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxxxz_0[i] = 2.0 * g_xxyzz_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxxz_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxxyy_0[i] = g_xxyyy_0_xxxxyy_0[i] * fbe_0 - g_xxyyy_0_xxxxyy_1[i] * fz_be_0 + g_xxyyyz_0_xxxxyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxxyz_0[i] = g_yyyzz_0_xxxxyz_0[i] * fbe_0 - g_yyyzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyyyzz_0_xxxyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxxzz_0[i] = 2.0 * g_xxyzz_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxzz_1[i] * fz_be_0 + g_xxyyzz_0_xxxxzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxyyy_0[i] = g_xxyyy_0_xxxyyy_0[i] * fbe_0 - g_xxyyy_0_xxxyyy_1[i] * fz_be_0 + g_xxyyyz_0_xxxyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxyyz_0[i] = g_yyyzz_0_xxxyyz_0[i] * fbe_0 - g_yyyzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_0_xxyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxyzz_0[i] = g_yyyzz_0_xxxyzz_0[i] * fbe_0 - g_yyyzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_0_xxyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxzzz_0[i] = 2.0 * g_xxyzz_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxzzz_1[i] * fz_be_0 + g_xxyyzz_0_xxxzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxyyyy_0[i] = g_xxyyy_0_xxyyyy_0[i] * fbe_0 - g_xxyyy_0_xxyyyy_1[i] * fz_be_0 + g_xxyyyz_0_xxyyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxyyyz_0[i] = g_yyyzz_0_xxyyyz_0[i] * fbe_0 - g_yyyzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxyyzz_0[i] = g_yyyzz_0_xxyyzz_0[i] * fbe_0 - g_yyyzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxyzzz_0[i] = g_yyyzz_0_xxyzzz_0[i] * fbe_0 - g_yyyzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxzzzz_0[i] = 2.0 * g_xxyzz_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxzzzz_1[i] * fz_be_0 + g_xxyyzz_0_xxzzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xyyyyy_0[i] = g_xxyyy_0_xyyyyy_0[i] * fbe_0 - g_xxyyy_0_xyyyyy_1[i] * fz_be_0 + g_xxyyyz_0_xyyyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xyyyyz_0[i] = g_yyyzz_0_xyyyyz_0[i] * fbe_0 - g_yyyzz_0_xyyyyz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyyyzz_0[i] = g_yyyzz_0_xyyyzz_0[i] * fbe_0 - g_yyyzz_0_xyyyzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyyzzz_0[i] = g_yyyzz_0_xyyzzz_0[i] * fbe_0 - g_yyyzz_0_xyyzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyzzzz_0[i] = g_yyyzz_0_xyzzzz_0[i] * fbe_0 - g_yyyzz_0_xyzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yzzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xzzzzz_0[i] = 2.0 * g_xxyzz_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_xzzzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_yyyyyy_0[i] = g_yyyzz_0_yyyyyy_0[i] * fbe_0 - g_yyyzz_0_yyyyyy_1[i] * fz_be_0 + g_xyyyzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyyyz_0[i] = g_yyyzz_0_yyyyyz_0[i] * fbe_0 - g_yyyzz_0_yyyyyz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyyzz_0[i] = g_yyyzz_0_yyyyzz_0[i] * fbe_0 - g_yyyzz_0_yyyyzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyzzz_0[i] = g_yyyzz_0_yyyzzz_0[i] * fbe_0 - g_yyyzz_0_yyyzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyzzzz_0[i] = g_yyyzz_0_yyzzzz_0[i] * fbe_0 - g_yyyzz_0_yyzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yzzzzz_0[i] = g_yyyzz_0_yzzzzz_0[i] * fbe_0 - g_yyyzz_0_yzzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_zzzzzz_0[i] = g_yyyzz_0_zzzzzz_0[i] * fbe_0 - g_yyyzz_0_zzzzzz_1[i] * fz_be_0 + g_xyyyzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 504-532 components of targeted buffer : KSI

    auto g_xxyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 504);

    auto g_xxyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 505);

    auto g_xxyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 506);

    auto g_xxyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 507);

    auto g_xxyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 508);

    auto g_xxyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 509);

    auto g_xxyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 510);

    auto g_xxyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 511);

    auto g_xxyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 512);

    auto g_xxyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 513);

    auto g_xxyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 514);

    auto g_xxyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 515);

    auto g_xxyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 516);

    auto g_xxyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 517);

    auto g_xxyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 518);

    auto g_xxyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 519);

    auto g_xxyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 520);

    auto g_xxyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 521);

    auto g_xxyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 522);

    auto g_xxyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 523);

    auto g_xxyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 524);

    auto g_xxyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 525);

    auto g_xxyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 526);

    auto g_xxyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 527);

    auto g_xxyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 528);

    auto g_xxyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 529);

    auto g_xxyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 530);

    auto g_xxyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 531);

    #pragma omp simd aligned(g_xxyyz_0_xxxxxy_0, g_xxyyz_0_xxxxxy_1, g_xxyyz_0_xxxxyy_0, g_xxyyz_0_xxxxyy_1, g_xxyyz_0_xxxyyy_0, g_xxyyz_0_xxxyyy_1, g_xxyyz_0_xxyyyy_0, g_xxyyz_0_xxyyyy_1, g_xxyyz_0_xyyyyy_0, g_xxyyz_0_xyyyyy_1, g_xxyyzz_0_xxxxxy_1, g_xxyyzz_0_xxxxyy_1, g_xxyyzz_0_xxxyyy_1, g_xxyyzz_0_xxyyyy_1, g_xxyyzz_0_xyyyyy_1, g_xxyyzzz_0_xxxxxx_0, g_xxyyzzz_0_xxxxxy_0, g_xxyyzzz_0_xxxxxz_0, g_xxyyzzz_0_xxxxyy_0, g_xxyyzzz_0_xxxxyz_0, g_xxyyzzz_0_xxxxzz_0, g_xxyyzzz_0_xxxyyy_0, g_xxyyzzz_0_xxxyyz_0, g_xxyyzzz_0_xxxyzz_0, g_xxyyzzz_0_xxxzzz_0, g_xxyyzzz_0_xxyyyy_0, g_xxyyzzz_0_xxyyyz_0, g_xxyyzzz_0_xxyyzz_0, g_xxyyzzz_0_xxyzzz_0, g_xxyyzzz_0_xxzzzz_0, g_xxyyzzz_0_xyyyyy_0, g_xxyyzzz_0_xyyyyz_0, g_xxyyzzz_0_xyyyzz_0, g_xxyyzzz_0_xyyzzz_0, g_xxyyzzz_0_xyzzzz_0, g_xxyyzzz_0_xzzzzz_0, g_xxyyzzz_0_yyyyyy_0, g_xxyyzzz_0_yyyyyz_0, g_xxyyzzz_0_yyyyzz_0, g_xxyyzzz_0_yyyzzz_0, g_xxyyzzz_0_yyzzzz_0, g_xxyyzzz_0_yzzzzz_0, g_xxyyzzz_0_zzzzzz_0, g_xxyzzz_0_xxxxxx_1, g_xxyzzz_0_xxxxxz_1, g_xxyzzz_0_xxxxzz_1, g_xxyzzz_0_xxxzzz_1, g_xxyzzz_0_xxzzzz_1, g_xxyzzz_0_xzzzzz_1, g_xxzzz_0_xxxxxx_0, g_xxzzz_0_xxxxxx_1, g_xxzzz_0_xxxxxz_0, g_xxzzz_0_xxxxxz_1, g_xxzzz_0_xxxxzz_0, g_xxzzz_0_xxxxzz_1, g_xxzzz_0_xxxzzz_0, g_xxzzz_0_xxxzzz_1, g_xxzzz_0_xxzzzz_0, g_xxzzz_0_xxzzzz_1, g_xxzzz_0_xzzzzz_0, g_xxzzz_0_xzzzzz_1, g_xyyzzz_0_xxxxyz_1, g_xyyzzz_0_xxxyyz_1, g_xyyzzz_0_xxxyz_1, g_xyyzzz_0_xxxyzz_1, g_xyyzzz_0_xxyyyz_1, g_xyyzzz_0_xxyyz_1, g_xyyzzz_0_xxyyzz_1, g_xyyzzz_0_xxyzz_1, g_xyyzzz_0_xxyzzz_1, g_xyyzzz_0_xyyyyz_1, g_xyyzzz_0_xyyyz_1, g_xyyzzz_0_xyyyzz_1, g_xyyzzz_0_xyyzz_1, g_xyyzzz_0_xyyzzz_1, g_xyyzzz_0_xyzzz_1, g_xyyzzz_0_xyzzzz_1, g_xyyzzz_0_yyyyyy_1, g_xyyzzz_0_yyyyyz_1, g_xyyzzz_0_yyyyz_1, g_xyyzzz_0_yyyyzz_1, g_xyyzzz_0_yyyzz_1, g_xyyzzz_0_yyyzzz_1, g_xyyzzz_0_yyzzz_1, g_xyyzzz_0_yyzzzz_1, g_xyyzzz_0_yzzzz_1, g_xyyzzz_0_yzzzzz_1, g_xyyzzz_0_zzzzzz_1, g_yyzzz_0_xxxxyz_0, g_yyzzz_0_xxxxyz_1, g_yyzzz_0_xxxyyz_0, g_yyzzz_0_xxxyyz_1, g_yyzzz_0_xxxyzz_0, g_yyzzz_0_xxxyzz_1, g_yyzzz_0_xxyyyz_0, g_yyzzz_0_xxyyyz_1, g_yyzzz_0_xxyyzz_0, g_yyzzz_0_xxyyzz_1, g_yyzzz_0_xxyzzz_0, g_yyzzz_0_xxyzzz_1, g_yyzzz_0_xyyyyz_0, g_yyzzz_0_xyyyyz_1, g_yyzzz_0_xyyyzz_0, g_yyzzz_0_xyyyzz_1, g_yyzzz_0_xyyzzz_0, g_yyzzz_0_xyyzzz_1, g_yyzzz_0_xyzzzz_0, g_yyzzz_0_xyzzzz_1, g_yyzzz_0_yyyyyy_0, g_yyzzz_0_yyyyyy_1, g_yyzzz_0_yyyyyz_0, g_yyzzz_0_yyyyyz_1, g_yyzzz_0_yyyyzz_0, g_yyzzz_0_yyyyzz_1, g_yyzzz_0_yyyzzz_0, g_yyzzz_0_yyyzzz_1, g_yyzzz_0_yyzzzz_0, g_yyzzz_0_yyzzzz_1, g_yyzzz_0_yzzzzz_0, g_yyzzz_0_yzzzzz_1, g_yyzzz_0_zzzzzz_0, g_yyzzz_0_zzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyzzz_0_xxxxxx_0[i] = g_xxzzz_0_xxxxxx_0[i] * fbe_0 - g_xxzzz_0_xxxxxx_1[i] * fz_be_0 + g_xxyzzz_0_xxxxxx_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxxxy_0[i] = 2.0 * g_xxyyz_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxxxy_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxxxz_0[i] = g_xxzzz_0_xxxxxz_0[i] * fbe_0 - g_xxzzz_0_xxxxxz_1[i] * fz_be_0 + g_xxyzzz_0_xxxxxz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxxyy_0[i] = 2.0 * g_xxyyz_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxxyy_1[i] * fz_be_0 + g_xxyyzz_0_xxxxyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxxyz_0[i] = g_yyzzz_0_xxxxyz_0[i] * fbe_0 - g_yyzzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyyzzz_0_xxxyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxxzz_0[i] = g_xxzzz_0_xxxxzz_0[i] * fbe_0 - g_xxzzz_0_xxxxzz_1[i] * fz_be_0 + g_xxyzzz_0_xxxxzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxyyy_0[i] = 2.0 * g_xxyyz_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxyyy_1[i] * fz_be_0 + g_xxyyzz_0_xxxyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxyyz_0[i] = g_yyzzz_0_xxxyyz_0[i] * fbe_0 - g_yyzzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_0_xxyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxyzz_0[i] = g_yyzzz_0_xxxyzz_0[i] * fbe_0 - g_yyzzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_0_xxyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxzzz_0[i] = g_xxzzz_0_xxxzzz_0[i] * fbe_0 - g_xxzzz_0_xxxzzz_1[i] * fz_be_0 + g_xxyzzz_0_xxxzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxyyyy_0[i] = 2.0 * g_xxyyz_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxyyyy_1[i] * fz_be_0 + g_xxyyzz_0_xxyyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxyyyz_0[i] = g_yyzzz_0_xxyyyz_0[i] * fbe_0 - g_yyzzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxyyzz_0[i] = g_yyzzz_0_xxyyzz_0[i] * fbe_0 - g_yyzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxyzzz_0[i] = g_yyzzz_0_xxyzzz_0[i] * fbe_0 - g_yyzzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxzzzz_0[i] = g_xxzzz_0_xxzzzz_0[i] * fbe_0 - g_xxzzz_0_xxzzzz_1[i] * fz_be_0 + g_xxyzzz_0_xxzzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xyyyyy_0[i] = 2.0 * g_xxyyz_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xyyyyy_1[i] * fz_be_0 + g_xxyyzz_0_xyyyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xyyyyz_0[i] = g_yyzzz_0_xyyyyz_0[i] * fbe_0 - g_yyzzz_0_xyyyyz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyyyzz_0[i] = g_yyzzz_0_xyyyzz_0[i] * fbe_0 - g_yyzzz_0_xyyyzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyyzzz_0[i] = g_yyzzz_0_xyyzzz_0[i] * fbe_0 - g_yyzzz_0_xyyzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyzzzz_0[i] = g_yyzzz_0_xyzzzz_0[i] * fbe_0 - g_yyzzz_0_xyzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yzzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xzzzzz_0[i] = g_xxzzz_0_xzzzzz_0[i] * fbe_0 - g_xxzzz_0_xzzzzz_1[i] * fz_be_0 + g_xxyzzz_0_xzzzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_yyyyyy_0[i] = g_yyzzz_0_yyyyyy_0[i] * fbe_0 - g_yyzzz_0_yyyyyy_1[i] * fz_be_0 + g_xyyzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyyyz_0[i] = g_yyzzz_0_yyyyyz_0[i] * fbe_0 - g_yyzzz_0_yyyyyz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyyzz_0[i] = g_yyzzz_0_yyyyzz_0[i] * fbe_0 - g_yyzzz_0_yyyyzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyzzz_0[i] = g_yyzzz_0_yyyzzz_0[i] * fbe_0 - g_yyzzz_0_yyyzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyzzzz_0[i] = g_yyzzz_0_yyzzzz_0[i] * fbe_0 - g_yyzzz_0_yyzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yzzzzz_0[i] = g_yyzzz_0_yzzzzz_0[i] * fbe_0 - g_yyzzz_0_yzzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_zzzzzz_0[i] = g_yyzzz_0_zzzzzz_0[i] * fbe_0 - g_yyzzz_0_zzzzzz_1[i] * fz_be_0 + g_xyyzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 532-560 components of targeted buffer : KSI

    auto g_xxyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 532);

    auto g_xxyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 533);

    auto g_xxyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 534);

    auto g_xxyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 535);

    auto g_xxyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 536);

    auto g_xxyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 537);

    auto g_xxyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 538);

    auto g_xxyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 539);

    auto g_xxyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 540);

    auto g_xxyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 541);

    auto g_xxyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 542);

    auto g_xxyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 543);

    auto g_xxyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 544);

    auto g_xxyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 545);

    auto g_xxyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 546);

    auto g_xxyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 547);

    auto g_xxyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 548);

    auto g_xxyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 549);

    auto g_xxyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 550);

    auto g_xxyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 551);

    auto g_xxyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 552);

    auto g_xxyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 553);

    auto g_xxyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 554);

    auto g_xxyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 555);

    auto g_xxyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 556);

    auto g_xxyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 557);

    auto g_xxyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 558);

    auto g_xxyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 559);

    #pragma omp simd aligned(g_xxyzzzz_0_xxxxxx_0, g_xxyzzzz_0_xxxxxy_0, g_xxyzzzz_0_xxxxxz_0, g_xxyzzzz_0_xxxxyy_0, g_xxyzzzz_0_xxxxyz_0, g_xxyzzzz_0_xxxxzz_0, g_xxyzzzz_0_xxxyyy_0, g_xxyzzzz_0_xxxyyz_0, g_xxyzzzz_0_xxxyzz_0, g_xxyzzzz_0_xxxzzz_0, g_xxyzzzz_0_xxyyyy_0, g_xxyzzzz_0_xxyyyz_0, g_xxyzzzz_0_xxyyzz_0, g_xxyzzzz_0_xxyzzz_0, g_xxyzzzz_0_xxzzzz_0, g_xxyzzzz_0_xyyyyy_0, g_xxyzzzz_0_xyyyyz_0, g_xxyzzzz_0_xyyyzz_0, g_xxyzzzz_0_xyyzzz_0, g_xxyzzzz_0_xyzzzz_0, g_xxyzzzz_0_xzzzzz_0, g_xxyzzzz_0_yyyyyy_0, g_xxyzzzz_0_yyyyyz_0, g_xxyzzzz_0_yyyyzz_0, g_xxyzzzz_0_yyyzzz_0, g_xxyzzzz_0_yyzzzz_0, g_xxyzzzz_0_yzzzzz_0, g_xxyzzzz_0_zzzzzz_0, g_xxzzzz_0_xxxxx_1, g_xxzzzz_0_xxxxxx_1, g_xxzzzz_0_xxxxxy_1, g_xxzzzz_0_xxxxxz_1, g_xxzzzz_0_xxxxy_1, g_xxzzzz_0_xxxxyy_1, g_xxzzzz_0_xxxxyz_1, g_xxzzzz_0_xxxxz_1, g_xxzzzz_0_xxxxzz_1, g_xxzzzz_0_xxxyy_1, g_xxzzzz_0_xxxyyy_1, g_xxzzzz_0_xxxyyz_1, g_xxzzzz_0_xxxyz_1, g_xxzzzz_0_xxxyzz_1, g_xxzzzz_0_xxxzz_1, g_xxzzzz_0_xxxzzz_1, g_xxzzzz_0_xxyyy_1, g_xxzzzz_0_xxyyyy_1, g_xxzzzz_0_xxyyyz_1, g_xxzzzz_0_xxyyz_1, g_xxzzzz_0_xxyyzz_1, g_xxzzzz_0_xxyzz_1, g_xxzzzz_0_xxyzzz_1, g_xxzzzz_0_xxzzz_1, g_xxzzzz_0_xxzzzz_1, g_xxzzzz_0_xyyyy_1, g_xxzzzz_0_xyyyyy_1, g_xxzzzz_0_xyyyyz_1, g_xxzzzz_0_xyyyz_1, g_xxzzzz_0_xyyyzz_1, g_xxzzzz_0_xyyzz_1, g_xxzzzz_0_xyyzzz_1, g_xxzzzz_0_xyzzz_1, g_xxzzzz_0_xyzzzz_1, g_xxzzzz_0_xzzzz_1, g_xxzzzz_0_xzzzzz_1, g_xxzzzz_0_yyyyy_1, g_xxzzzz_0_yyyyyy_1, g_xxzzzz_0_yyyyyz_1, g_xxzzzz_0_yyyyz_1, g_xxzzzz_0_yyyyzz_1, g_xxzzzz_0_yyyzz_1, g_xxzzzz_0_yyyzzz_1, g_xxzzzz_0_yyzzz_1, g_xxzzzz_0_yyzzzz_1, g_xxzzzz_0_yzzzz_1, g_xxzzzz_0_yzzzzz_1, g_xxzzzz_0_zzzzz_1, g_xxzzzz_0_zzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzzz_0_xxxxxx_0[i] = g_xxzzzz_0_xxxxxx_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxy_0[i] = g_xxzzzz_0_xxxxx_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxz_0[i] = g_xxzzzz_0_xxxxxz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxyy_0[i] = 2.0 * g_xxzzzz_0_xxxxy_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxyz_0[i] = g_xxzzzz_0_xxxxz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxzz_0[i] = g_xxzzzz_0_xxxxzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyyy_0[i] = 3.0 * g_xxzzzz_0_xxxyy_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyyz_0[i] = 2.0 * g_xxzzzz_0_xxxyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyzz_0[i] = g_xxzzzz_0_xxxzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxzzz_0[i] = g_xxzzzz_0_xxxzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyyy_0[i] = 4.0 * g_xxzzzz_0_xxyyy_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyyz_0[i] = 3.0 * g_xxzzzz_0_xxyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyzz_0[i] = 2.0 * g_xxzzzz_0_xxyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyzzz_0[i] = g_xxzzzz_0_xxzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxzzzz_0[i] = g_xxzzzz_0_xxzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyyy_0[i] = 5.0 * g_xxzzzz_0_xyyyy_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyyz_0[i] = 4.0 * g_xxzzzz_0_xyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyzz_0[i] = 3.0 * g_xxzzzz_0_xyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyzzz_0[i] = 2.0 * g_xxzzzz_0_xyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyzzzz_0[i] = g_xxzzzz_0_xzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xzzzzz_0[i] = g_xxzzzz_0_xzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyyy_0[i] = 6.0 * g_xxzzzz_0_yyyyy_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyyz_0[i] = 5.0 * g_xxzzzz_0_yyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyzz_0[i] = 4.0 * g_xxzzzz_0_yyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyzzz_0[i] = 3.0 * g_xxzzzz_0_yyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyzzzz_0[i] = 2.0 * g_xxzzzz_0_yzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yzzzzz_0[i] = g_xxzzzz_0_zzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_zzzzzz_0[i] = g_xxzzzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 560-588 components of targeted buffer : KSI

    auto g_xxzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 560);

    auto g_xxzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 561);

    auto g_xxzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 562);

    auto g_xxzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 563);

    auto g_xxzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 564);

    auto g_xxzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 565);

    auto g_xxzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 566);

    auto g_xxzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 567);

    auto g_xxzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 568);

    auto g_xxzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 569);

    auto g_xxzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 570);

    auto g_xxzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 571);

    auto g_xxzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 572);

    auto g_xxzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 573);

    auto g_xxzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 574);

    auto g_xxzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 575);

    auto g_xxzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 576);

    auto g_xxzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 577);

    auto g_xxzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 578);

    auto g_xxzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 579);

    auto g_xxzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 580);

    auto g_xxzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 581);

    auto g_xxzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 582);

    auto g_xxzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 583);

    auto g_xxzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 584);

    auto g_xxzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 585);

    auto g_xxzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 586);

    auto g_xxzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 587);

    #pragma omp simd aligned(g_xxzzz_0_xxxxxx_0, g_xxzzz_0_xxxxxx_1, g_xxzzz_0_xxxxxy_0, g_xxzzz_0_xxxxxy_1, g_xxzzz_0_xxxxyy_0, g_xxzzz_0_xxxxyy_1, g_xxzzz_0_xxxyyy_0, g_xxzzz_0_xxxyyy_1, g_xxzzz_0_xxyyyy_0, g_xxzzz_0_xxyyyy_1, g_xxzzz_0_xyyyyy_0, g_xxzzz_0_xyyyyy_1, g_xxzzzz_0_xxxxxx_1, g_xxzzzz_0_xxxxxy_1, g_xxzzzz_0_xxxxyy_1, g_xxzzzz_0_xxxyyy_1, g_xxzzzz_0_xxyyyy_1, g_xxzzzz_0_xyyyyy_1, g_xxzzzzz_0_xxxxxx_0, g_xxzzzzz_0_xxxxxy_0, g_xxzzzzz_0_xxxxxz_0, g_xxzzzzz_0_xxxxyy_0, g_xxzzzzz_0_xxxxyz_0, g_xxzzzzz_0_xxxxzz_0, g_xxzzzzz_0_xxxyyy_0, g_xxzzzzz_0_xxxyyz_0, g_xxzzzzz_0_xxxyzz_0, g_xxzzzzz_0_xxxzzz_0, g_xxzzzzz_0_xxyyyy_0, g_xxzzzzz_0_xxyyyz_0, g_xxzzzzz_0_xxyyzz_0, g_xxzzzzz_0_xxyzzz_0, g_xxzzzzz_0_xxzzzz_0, g_xxzzzzz_0_xyyyyy_0, g_xxzzzzz_0_xyyyyz_0, g_xxzzzzz_0_xyyyzz_0, g_xxzzzzz_0_xyyzzz_0, g_xxzzzzz_0_xyzzzz_0, g_xxzzzzz_0_xzzzzz_0, g_xxzzzzz_0_yyyyyy_0, g_xxzzzzz_0_yyyyyz_0, g_xxzzzzz_0_yyyyzz_0, g_xxzzzzz_0_yyyzzz_0, g_xxzzzzz_0_yyzzzz_0, g_xxzzzzz_0_yzzzzz_0, g_xxzzzzz_0_zzzzzz_0, g_xzzzzz_0_xxxxxz_1, g_xzzzzz_0_xxxxyz_1, g_xzzzzz_0_xxxxz_1, g_xzzzzz_0_xxxxzz_1, g_xzzzzz_0_xxxyyz_1, g_xzzzzz_0_xxxyz_1, g_xzzzzz_0_xxxyzz_1, g_xzzzzz_0_xxxzz_1, g_xzzzzz_0_xxxzzz_1, g_xzzzzz_0_xxyyyz_1, g_xzzzzz_0_xxyyz_1, g_xzzzzz_0_xxyyzz_1, g_xzzzzz_0_xxyzz_1, g_xzzzzz_0_xxyzzz_1, g_xzzzzz_0_xxzzz_1, g_xzzzzz_0_xxzzzz_1, g_xzzzzz_0_xyyyyz_1, g_xzzzzz_0_xyyyz_1, g_xzzzzz_0_xyyyzz_1, g_xzzzzz_0_xyyzz_1, g_xzzzzz_0_xyyzzz_1, g_xzzzzz_0_xyzzz_1, g_xzzzzz_0_xyzzzz_1, g_xzzzzz_0_xzzzz_1, g_xzzzzz_0_xzzzzz_1, g_xzzzzz_0_yyyyyy_1, g_xzzzzz_0_yyyyyz_1, g_xzzzzz_0_yyyyz_1, g_xzzzzz_0_yyyyzz_1, g_xzzzzz_0_yyyzz_1, g_xzzzzz_0_yyyzzz_1, g_xzzzzz_0_yyzzz_1, g_xzzzzz_0_yyzzzz_1, g_xzzzzz_0_yzzzz_1, g_xzzzzz_0_yzzzzz_1, g_xzzzzz_0_zzzzz_1, g_xzzzzz_0_zzzzzz_1, g_zzzzz_0_xxxxxz_0, g_zzzzz_0_xxxxxz_1, g_zzzzz_0_xxxxyz_0, g_zzzzz_0_xxxxyz_1, g_zzzzz_0_xxxxzz_0, g_zzzzz_0_xxxxzz_1, g_zzzzz_0_xxxyyz_0, g_zzzzz_0_xxxyyz_1, g_zzzzz_0_xxxyzz_0, g_zzzzz_0_xxxyzz_1, g_zzzzz_0_xxxzzz_0, g_zzzzz_0_xxxzzz_1, g_zzzzz_0_xxyyyz_0, g_zzzzz_0_xxyyyz_1, g_zzzzz_0_xxyyzz_0, g_zzzzz_0_xxyyzz_1, g_zzzzz_0_xxyzzz_0, g_zzzzz_0_xxyzzz_1, g_zzzzz_0_xxzzzz_0, g_zzzzz_0_xxzzzz_1, g_zzzzz_0_xyyyyz_0, g_zzzzz_0_xyyyyz_1, g_zzzzz_0_xyyyzz_0, g_zzzzz_0_xyyyzz_1, g_zzzzz_0_xyyzzz_0, g_zzzzz_0_xyyzzz_1, g_zzzzz_0_xyzzzz_0, g_zzzzz_0_xyzzzz_1, g_zzzzz_0_xzzzzz_0, g_zzzzz_0_xzzzzz_1, g_zzzzz_0_yyyyyy_0, g_zzzzz_0_yyyyyy_1, g_zzzzz_0_yyyyyz_0, g_zzzzz_0_yyyyyz_1, g_zzzzz_0_yyyyzz_0, g_zzzzz_0_yyyyzz_1, g_zzzzz_0_yyyzzz_0, g_zzzzz_0_yyyzzz_1, g_zzzzz_0_yyzzzz_0, g_zzzzz_0_yyzzzz_1, g_zzzzz_0_yzzzzz_0, g_zzzzz_0_yzzzzz_1, g_zzzzz_0_zzzzzz_0, g_zzzzz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzzz_0_xxxxxx_0[i] = 4.0 * g_xxzzz_0_xxxxxx_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxxx_1[i] * fz_be_0 + g_xxzzzz_0_xxxxxx_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxxy_0[i] = 4.0 * g_xxzzz_0_xxxxxy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxxy_1[i] * fz_be_0 + g_xxzzzz_0_xxxxxy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxxz_0[i] = g_zzzzz_0_xxxxxz_0[i] * fbe_0 - g_zzzzz_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xzzzzz_0_xxxxz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxxz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxyy_0[i] = 4.0 * g_xxzzz_0_xxxxyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxyy_1[i] * fz_be_0 + g_xxzzzz_0_xxxxyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxyz_0[i] = g_zzzzz_0_xxxxyz_0[i] * fbe_0 - g_zzzzz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_0_xxxyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxzz_0[i] = g_zzzzz_0_xxxxzz_0[i] * fbe_0 - g_zzzzz_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_0_xxxzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxyyy_0[i] = 4.0 * g_xxzzz_0_xxxyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxyyy_1[i] * fz_be_0 + g_xxzzzz_0_xxxyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxyyz_0[i] = g_zzzzz_0_xxxyyz_0[i] * fbe_0 - g_zzzzz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxyzz_0[i] = g_zzzzz_0_xxxyzz_0[i] * fbe_0 - g_zzzzz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxzzz_0[i] = g_zzzzz_0_xxxzzz_0[i] * fbe_0 - g_zzzzz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyyyy_0[i] = 4.0 * g_xxzzz_0_xxyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxyyyy_1[i] * fz_be_0 + g_xxzzzz_0_xxyyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxyyyz_0[i] = g_zzzzz_0_xxyyyz_0[i] * fbe_0 - g_zzzzz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyyzz_0[i] = g_zzzzz_0_xxyyzz_0[i] * fbe_0 - g_zzzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyzzz_0[i] = g_zzzzz_0_xxyzzz_0[i] * fbe_0 - g_zzzzz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxzzzz_0[i] = g_zzzzz_0_xxzzzz_0[i] * fbe_0 - g_zzzzz_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyyyy_0[i] = 4.0 * g_xxzzz_0_xyyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xyyyyy_1[i] * fz_be_0 + g_xxzzzz_0_xyyyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xyyyyz_0[i] = g_zzzzz_0_xyyyyz_0[i] * fbe_0 - g_zzzzz_0_xyyyyz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyyzz_0[i] = g_zzzzz_0_xyyyzz_0[i] * fbe_0 - g_zzzzz_0_xyyyzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyzzz_0[i] = g_zzzzz_0_xyyzzz_0[i] * fbe_0 - g_zzzzz_0_xyyzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyzzzz_0[i] = g_zzzzz_0_xyzzzz_0[i] * fbe_0 - g_zzzzz_0_xyzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xzzzzz_0[i] = g_zzzzz_0_xzzzzz_0[i] * fbe_0 - g_zzzzz_0_xzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_zzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyyy_0[i] = g_zzzzz_0_yyyyyy_0[i] * fbe_0 - g_zzzzz_0_yyyyyy_1[i] * fz_be_0 + g_xzzzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyyz_0[i] = g_zzzzz_0_yyyyyz_0[i] * fbe_0 - g_zzzzz_0_yyyyyz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyzz_0[i] = g_zzzzz_0_yyyyzz_0[i] * fbe_0 - g_zzzzz_0_yyyyzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyzzz_0[i] = g_zzzzz_0_yyyzzz_0[i] * fbe_0 - g_zzzzz_0_yyyzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyzzzz_0[i] = g_zzzzz_0_yyzzzz_0[i] * fbe_0 - g_zzzzz_0_yyzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yzzzzz_0[i] = g_zzzzz_0_yzzzzz_0[i] * fbe_0 - g_zzzzz_0_yzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_zzzzzz_0[i] = g_zzzzz_0_zzzzzz_0[i] * fbe_0 - g_zzzzz_0_zzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 588-616 components of targeted buffer : KSI

    auto g_xyyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 588);

    auto g_xyyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 589);

    auto g_xyyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 590);

    auto g_xyyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 591);

    auto g_xyyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 592);

    auto g_xyyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 593);

    auto g_xyyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 594);

    auto g_xyyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 595);

    auto g_xyyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 596);

    auto g_xyyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 597);

    auto g_xyyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 598);

    auto g_xyyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 599);

    auto g_xyyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 600);

    auto g_xyyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 601);

    auto g_xyyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 602);

    auto g_xyyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 603);

    auto g_xyyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 604);

    auto g_xyyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 605);

    auto g_xyyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 606);

    auto g_xyyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 607);

    auto g_xyyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 608);

    auto g_xyyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 609);

    auto g_xyyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 610);

    auto g_xyyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 611);

    auto g_xyyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 612);

    auto g_xyyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 613);

    auto g_xyyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 614);

    auto g_xyyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 615);

    #pragma omp simd aligned(g_xyyyyyy_0_xxxxxx_0, g_xyyyyyy_0_xxxxxy_0, g_xyyyyyy_0_xxxxxz_0, g_xyyyyyy_0_xxxxyy_0, g_xyyyyyy_0_xxxxyz_0, g_xyyyyyy_0_xxxxzz_0, g_xyyyyyy_0_xxxyyy_0, g_xyyyyyy_0_xxxyyz_0, g_xyyyyyy_0_xxxyzz_0, g_xyyyyyy_0_xxxzzz_0, g_xyyyyyy_0_xxyyyy_0, g_xyyyyyy_0_xxyyyz_0, g_xyyyyyy_0_xxyyzz_0, g_xyyyyyy_0_xxyzzz_0, g_xyyyyyy_0_xxzzzz_0, g_xyyyyyy_0_xyyyyy_0, g_xyyyyyy_0_xyyyyz_0, g_xyyyyyy_0_xyyyzz_0, g_xyyyyyy_0_xyyzzz_0, g_xyyyyyy_0_xyzzzz_0, g_xyyyyyy_0_xzzzzz_0, g_xyyyyyy_0_yyyyyy_0, g_xyyyyyy_0_yyyyyz_0, g_xyyyyyy_0_yyyyzz_0, g_xyyyyyy_0_yyyzzz_0, g_xyyyyyy_0_yyzzzz_0, g_xyyyyyy_0_yzzzzz_0, g_xyyyyyy_0_zzzzzz_0, g_yyyyyy_0_xxxxx_1, g_yyyyyy_0_xxxxxx_1, g_yyyyyy_0_xxxxxy_1, g_yyyyyy_0_xxxxxz_1, g_yyyyyy_0_xxxxy_1, g_yyyyyy_0_xxxxyy_1, g_yyyyyy_0_xxxxyz_1, g_yyyyyy_0_xxxxz_1, g_yyyyyy_0_xxxxzz_1, g_yyyyyy_0_xxxyy_1, g_yyyyyy_0_xxxyyy_1, g_yyyyyy_0_xxxyyz_1, g_yyyyyy_0_xxxyz_1, g_yyyyyy_0_xxxyzz_1, g_yyyyyy_0_xxxzz_1, g_yyyyyy_0_xxxzzz_1, g_yyyyyy_0_xxyyy_1, g_yyyyyy_0_xxyyyy_1, g_yyyyyy_0_xxyyyz_1, g_yyyyyy_0_xxyyz_1, g_yyyyyy_0_xxyyzz_1, g_yyyyyy_0_xxyzz_1, g_yyyyyy_0_xxyzzz_1, g_yyyyyy_0_xxzzz_1, g_yyyyyy_0_xxzzzz_1, g_yyyyyy_0_xyyyy_1, g_yyyyyy_0_xyyyyy_1, g_yyyyyy_0_xyyyyz_1, g_yyyyyy_0_xyyyz_1, g_yyyyyy_0_xyyyzz_1, g_yyyyyy_0_xyyzz_1, g_yyyyyy_0_xyyzzz_1, g_yyyyyy_0_xyzzz_1, g_yyyyyy_0_xyzzzz_1, g_yyyyyy_0_xzzzz_1, g_yyyyyy_0_xzzzzz_1, g_yyyyyy_0_yyyyy_1, g_yyyyyy_0_yyyyyy_1, g_yyyyyy_0_yyyyyz_1, g_yyyyyy_0_yyyyz_1, g_yyyyyy_0_yyyyzz_1, g_yyyyyy_0_yyyzz_1, g_yyyyyy_0_yyyzzz_1, g_yyyyyy_0_yyzzz_1, g_yyyyyy_0_yyzzzz_1, g_yyyyyy_0_yzzzz_1, g_yyyyyy_0_yzzzzz_1, g_yyyyyy_0_zzzzz_1, g_yyyyyy_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyy_0_xxxxxx_0[i] = 6.0 * g_yyyyyy_0_xxxxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxx_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxy_0[i] = 5.0 * g_yyyyyy_0_xxxxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxz_0[i] = 5.0 * g_yyyyyy_0_xxxxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxyy_0[i] = 4.0 * g_yyyyyy_0_xxxyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxyz_0[i] = 4.0 * g_yyyyyy_0_xxxyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxzz_0[i] = 4.0 * g_yyyyyy_0_xxxzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyyy_0[i] = 3.0 * g_yyyyyy_0_xxyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyyz_0[i] = 3.0 * g_yyyyyy_0_xxyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyzz_0[i] = 3.0 * g_yyyyyy_0_xxyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxzzz_0[i] = 3.0 * g_yyyyyy_0_xxzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyyy_0[i] = 2.0 * g_yyyyyy_0_xyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyyz_0[i] = 2.0 * g_yyyyyy_0_xyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyzz_0[i] = 2.0 * g_yyyyyy_0_xyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyzzz_0[i] = 2.0 * g_yyyyyy_0_xyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxzzzz_0[i] = 2.0 * g_yyyyyy_0_xzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyyy_0[i] = g_yyyyyy_0_yyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyyz_0[i] = g_yyyyyy_0_yyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyzz_0[i] = g_yyyyyy_0_yyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyzzz_0[i] = g_yyyyyy_0_yyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyzzzz_0[i] = g_yyyyyy_0_yzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xzzzzz_0[i] = g_yyyyyy_0_zzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyyy_0[i] = g_yyyyyy_0_yyyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyyz_0[i] = g_yyyyyy_0_yyyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyzz_0[i] = g_yyyyyy_0_yyyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyzzz_0[i] = g_yyyyyy_0_yyyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyzzzz_0[i] = g_yyyyyy_0_yyzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yzzzzz_0[i] = g_yyyyyy_0_yzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_zzzzzz_0[i] = g_yyyyyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 616-644 components of targeted buffer : KSI

    auto g_xyyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 616);

    auto g_xyyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 617);

    auto g_xyyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 618);

    auto g_xyyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 619);

    auto g_xyyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 620);

    auto g_xyyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 621);

    auto g_xyyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 622);

    auto g_xyyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 623);

    auto g_xyyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 624);

    auto g_xyyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 625);

    auto g_xyyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 626);

    auto g_xyyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 627);

    auto g_xyyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 628);

    auto g_xyyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 629);

    auto g_xyyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 630);

    auto g_xyyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 631);

    auto g_xyyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 632);

    auto g_xyyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 633);

    auto g_xyyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 634);

    auto g_xyyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 635);

    auto g_xyyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 636);

    auto g_xyyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 637);

    auto g_xyyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 638);

    auto g_xyyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 639);

    auto g_xyyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 640);

    auto g_xyyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 641);

    auto g_xyyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 642);

    auto g_xyyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 643);

    #pragma omp simd aligned(g_xyyyyy_0_xxxxxx_1, g_xyyyyy_0_xxxxxy_1, g_xyyyyy_0_xxxxyy_1, g_xyyyyy_0_xxxyyy_1, g_xyyyyy_0_xxyyyy_1, g_xyyyyy_0_xyyyyy_1, g_xyyyyyz_0_xxxxxx_0, g_xyyyyyz_0_xxxxxy_0, g_xyyyyyz_0_xxxxxz_0, g_xyyyyyz_0_xxxxyy_0, g_xyyyyyz_0_xxxxyz_0, g_xyyyyyz_0_xxxxzz_0, g_xyyyyyz_0_xxxyyy_0, g_xyyyyyz_0_xxxyyz_0, g_xyyyyyz_0_xxxyzz_0, g_xyyyyyz_0_xxxzzz_0, g_xyyyyyz_0_xxyyyy_0, g_xyyyyyz_0_xxyyyz_0, g_xyyyyyz_0_xxyyzz_0, g_xyyyyyz_0_xxyzzz_0, g_xyyyyyz_0_xxzzzz_0, g_xyyyyyz_0_xyyyyy_0, g_xyyyyyz_0_xyyyyz_0, g_xyyyyyz_0_xyyyzz_0, g_xyyyyyz_0_xyyzzz_0, g_xyyyyyz_0_xyzzzz_0, g_xyyyyyz_0_xzzzzz_0, g_xyyyyyz_0_yyyyyy_0, g_xyyyyyz_0_yyyyyz_0, g_xyyyyyz_0_yyyyzz_0, g_xyyyyyz_0_yyyzzz_0, g_xyyyyyz_0_yyzzzz_0, g_xyyyyyz_0_yzzzzz_0, g_xyyyyyz_0_zzzzzz_0, g_yyyyyz_0_xxxxxz_1, g_yyyyyz_0_xxxxyz_1, g_yyyyyz_0_xxxxz_1, g_yyyyyz_0_xxxxzz_1, g_yyyyyz_0_xxxyyz_1, g_yyyyyz_0_xxxyz_1, g_yyyyyz_0_xxxyzz_1, g_yyyyyz_0_xxxzz_1, g_yyyyyz_0_xxxzzz_1, g_yyyyyz_0_xxyyyz_1, g_yyyyyz_0_xxyyz_1, g_yyyyyz_0_xxyyzz_1, g_yyyyyz_0_xxyzz_1, g_yyyyyz_0_xxyzzz_1, g_yyyyyz_0_xxzzz_1, g_yyyyyz_0_xxzzzz_1, g_yyyyyz_0_xyyyyz_1, g_yyyyyz_0_xyyyz_1, g_yyyyyz_0_xyyyzz_1, g_yyyyyz_0_xyyzz_1, g_yyyyyz_0_xyyzzz_1, g_yyyyyz_0_xyzzz_1, g_yyyyyz_0_xyzzzz_1, g_yyyyyz_0_xzzzz_1, g_yyyyyz_0_xzzzzz_1, g_yyyyyz_0_yyyyyy_1, g_yyyyyz_0_yyyyyz_1, g_yyyyyz_0_yyyyz_1, g_yyyyyz_0_yyyyzz_1, g_yyyyyz_0_yyyzz_1, g_yyyyyz_0_yyyzzz_1, g_yyyyyz_0_yyzzz_1, g_yyyyyz_0_yyzzzz_1, g_yyyyyz_0_yzzzz_1, g_yyyyyz_0_yzzzzz_1, g_yyyyyz_0_zzzzz_1, g_yyyyyz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyz_0_xxxxxx_0[i] = g_xyyyyy_0_xxxxxx_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxxy_0[i] = g_xyyyyy_0_xxxxxy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxxz_0[i] = 5.0 * g_yyyyyz_0_xxxxz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxxz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxyy_0[i] = g_xyyyyy_0_xxxxyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxyz_0[i] = 4.0 * g_yyyyyz_0_xxxyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxzz_0[i] = 4.0 * g_yyyyyz_0_xxxzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxyyy_0[i] = g_xyyyyy_0_xxxyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxyyz_0[i] = 3.0 * g_yyyyyz_0_xxyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxyzz_0[i] = 3.0 * g_yyyyyz_0_xxyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxzzz_0[i] = 3.0 * g_yyyyyz_0_xxzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyyyy_0[i] = g_xyyyyy_0_xxyyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxyyyz_0[i] = 2.0 * g_yyyyyz_0_xyyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyyzz_0[i] = 2.0 * g_yyyyyz_0_xyyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyzzz_0[i] = 2.0 * g_yyyyyz_0_xyzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxzzzz_0[i] = 2.0 * g_yyyyyz_0_xzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyyyy_0[i] = g_xyyyyy_0_xyyyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xyyyyz_0[i] = g_yyyyyz_0_yyyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyyzz_0[i] = g_yyyyyz_0_yyyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyzzz_0[i] = g_yyyyyz_0_yyzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyzzzz_0[i] = g_yyyyyz_0_yzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xzzzzz_0[i] = g_yyyyyz_0_zzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyyy_0[i] = g_yyyyyz_0_yyyyyy_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyyz_0[i] = g_yyyyyz_0_yyyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyzz_0[i] = g_yyyyyz_0_yyyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyzzz_0[i] = g_yyyyyz_0_yyyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyzzzz_0[i] = g_yyyyyz_0_yyzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yzzzzz_0[i] = g_yyyyyz_0_yzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_zzzzzz_0[i] = g_yyyyyz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 644-672 components of targeted buffer : KSI

    auto g_xyyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 644);

    auto g_xyyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 645);

    auto g_xyyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 646);

    auto g_xyyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 647);

    auto g_xyyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 648);

    auto g_xyyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 649);

    auto g_xyyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 650);

    auto g_xyyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 651);

    auto g_xyyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 652);

    auto g_xyyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 653);

    auto g_xyyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 654);

    auto g_xyyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 655);

    auto g_xyyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 656);

    auto g_xyyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 657);

    auto g_xyyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 658);

    auto g_xyyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 659);

    auto g_xyyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 660);

    auto g_xyyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 661);

    auto g_xyyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 662);

    auto g_xyyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 663);

    auto g_xyyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 664);

    auto g_xyyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 665);

    auto g_xyyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 666);

    auto g_xyyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 667);

    auto g_xyyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 668);

    auto g_xyyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 669);

    auto g_xyyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 670);

    auto g_xyyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 671);

    #pragma omp simd aligned(g_xyyyyzz_0_xxxxxx_0, g_xyyyyzz_0_xxxxxy_0, g_xyyyyzz_0_xxxxxz_0, g_xyyyyzz_0_xxxxyy_0, g_xyyyyzz_0_xxxxyz_0, g_xyyyyzz_0_xxxxzz_0, g_xyyyyzz_0_xxxyyy_0, g_xyyyyzz_0_xxxyyz_0, g_xyyyyzz_0_xxxyzz_0, g_xyyyyzz_0_xxxzzz_0, g_xyyyyzz_0_xxyyyy_0, g_xyyyyzz_0_xxyyyz_0, g_xyyyyzz_0_xxyyzz_0, g_xyyyyzz_0_xxyzzz_0, g_xyyyyzz_0_xxzzzz_0, g_xyyyyzz_0_xyyyyy_0, g_xyyyyzz_0_xyyyyz_0, g_xyyyyzz_0_xyyyzz_0, g_xyyyyzz_0_xyyzzz_0, g_xyyyyzz_0_xyzzzz_0, g_xyyyyzz_0_xzzzzz_0, g_xyyyyzz_0_yyyyyy_0, g_xyyyyzz_0_yyyyyz_0, g_xyyyyzz_0_yyyyzz_0, g_xyyyyzz_0_yyyzzz_0, g_xyyyyzz_0_yyzzzz_0, g_xyyyyzz_0_yzzzzz_0, g_xyyyyzz_0_zzzzzz_0, g_yyyyzz_0_xxxxx_1, g_yyyyzz_0_xxxxxx_1, g_yyyyzz_0_xxxxxy_1, g_yyyyzz_0_xxxxxz_1, g_yyyyzz_0_xxxxy_1, g_yyyyzz_0_xxxxyy_1, g_yyyyzz_0_xxxxyz_1, g_yyyyzz_0_xxxxz_1, g_yyyyzz_0_xxxxzz_1, g_yyyyzz_0_xxxyy_1, g_yyyyzz_0_xxxyyy_1, g_yyyyzz_0_xxxyyz_1, g_yyyyzz_0_xxxyz_1, g_yyyyzz_0_xxxyzz_1, g_yyyyzz_0_xxxzz_1, g_yyyyzz_0_xxxzzz_1, g_yyyyzz_0_xxyyy_1, g_yyyyzz_0_xxyyyy_1, g_yyyyzz_0_xxyyyz_1, g_yyyyzz_0_xxyyz_1, g_yyyyzz_0_xxyyzz_1, g_yyyyzz_0_xxyzz_1, g_yyyyzz_0_xxyzzz_1, g_yyyyzz_0_xxzzz_1, g_yyyyzz_0_xxzzzz_1, g_yyyyzz_0_xyyyy_1, g_yyyyzz_0_xyyyyy_1, g_yyyyzz_0_xyyyyz_1, g_yyyyzz_0_xyyyz_1, g_yyyyzz_0_xyyyzz_1, g_yyyyzz_0_xyyzz_1, g_yyyyzz_0_xyyzzz_1, g_yyyyzz_0_xyzzz_1, g_yyyyzz_0_xyzzzz_1, g_yyyyzz_0_xzzzz_1, g_yyyyzz_0_xzzzzz_1, g_yyyyzz_0_yyyyy_1, g_yyyyzz_0_yyyyyy_1, g_yyyyzz_0_yyyyyz_1, g_yyyyzz_0_yyyyz_1, g_yyyyzz_0_yyyyzz_1, g_yyyyzz_0_yyyzz_1, g_yyyyzz_0_yyyzzz_1, g_yyyyzz_0_yyzzz_1, g_yyyyzz_0_yyzzzz_1, g_yyyyzz_0_yzzzz_1, g_yyyyzz_0_yzzzzz_1, g_yyyyzz_0_zzzzz_1, g_yyyyzz_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyzz_0_xxxxxx_0[i] = 6.0 * g_yyyyzz_0_xxxxx_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxx_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxy_0[i] = 5.0 * g_yyyyzz_0_xxxxy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxz_0[i] = 5.0 * g_yyyyzz_0_xxxxz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxyy_0[i] = 4.0 * g_yyyyzz_0_xxxyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxyz_0[i] = 4.0 * g_yyyyzz_0_xxxyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxzz_0[i] = 4.0 * g_yyyyzz_0_xxxzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyyy_0[i] = 3.0 * g_yyyyzz_0_xxyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyyz_0[i] = 3.0 * g_yyyyzz_0_xxyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyzz_0[i] = 3.0 * g_yyyyzz_0_xxyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxzzz_0[i] = 3.0 * g_yyyyzz_0_xxzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyyy_0[i] = 2.0 * g_yyyyzz_0_xyyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyyz_0[i] = 2.0 * g_yyyyzz_0_xyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyzz_0[i] = 2.0 * g_yyyyzz_0_xyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyzzz_0[i] = 2.0 * g_yyyyzz_0_xyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxzzzz_0[i] = 2.0 * g_yyyyzz_0_xzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyyy_0[i] = g_yyyyzz_0_yyyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyyz_0[i] = g_yyyyzz_0_yyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyzz_0[i] = g_yyyyzz_0_yyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyzzz_0[i] = g_yyyyzz_0_yyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyzzzz_0[i] = g_yyyyzz_0_yzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xzzzzz_0[i] = g_yyyyzz_0_zzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyyy_0[i] = g_yyyyzz_0_yyyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyyz_0[i] = g_yyyyzz_0_yyyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyzz_0[i] = g_yyyyzz_0_yyyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyzzz_0[i] = g_yyyyzz_0_yyyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyzzzz_0[i] = g_yyyyzz_0_yyzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yzzzzz_0[i] = g_yyyyzz_0_yzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_zzzzzz_0[i] = g_yyyyzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 672-700 components of targeted buffer : KSI

    auto g_xyyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 672);

    auto g_xyyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 673);

    auto g_xyyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 674);

    auto g_xyyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 675);

    auto g_xyyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 676);

    auto g_xyyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 677);

    auto g_xyyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 678);

    auto g_xyyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 679);

    auto g_xyyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 680);

    auto g_xyyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 681);

    auto g_xyyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 682);

    auto g_xyyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 683);

    auto g_xyyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 684);

    auto g_xyyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 685);

    auto g_xyyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 686);

    auto g_xyyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 687);

    auto g_xyyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 688);

    auto g_xyyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 689);

    auto g_xyyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 690);

    auto g_xyyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 691);

    auto g_xyyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 692);

    auto g_xyyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 693);

    auto g_xyyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 694);

    auto g_xyyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 695);

    auto g_xyyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 696);

    auto g_xyyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 697);

    auto g_xyyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 698);

    auto g_xyyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 699);

    #pragma omp simd aligned(g_xyyyzzz_0_xxxxxx_0, g_xyyyzzz_0_xxxxxy_0, g_xyyyzzz_0_xxxxxz_0, g_xyyyzzz_0_xxxxyy_0, g_xyyyzzz_0_xxxxyz_0, g_xyyyzzz_0_xxxxzz_0, g_xyyyzzz_0_xxxyyy_0, g_xyyyzzz_0_xxxyyz_0, g_xyyyzzz_0_xxxyzz_0, g_xyyyzzz_0_xxxzzz_0, g_xyyyzzz_0_xxyyyy_0, g_xyyyzzz_0_xxyyyz_0, g_xyyyzzz_0_xxyyzz_0, g_xyyyzzz_0_xxyzzz_0, g_xyyyzzz_0_xxzzzz_0, g_xyyyzzz_0_xyyyyy_0, g_xyyyzzz_0_xyyyyz_0, g_xyyyzzz_0_xyyyzz_0, g_xyyyzzz_0_xyyzzz_0, g_xyyyzzz_0_xyzzzz_0, g_xyyyzzz_0_xzzzzz_0, g_xyyyzzz_0_yyyyyy_0, g_xyyyzzz_0_yyyyyz_0, g_xyyyzzz_0_yyyyzz_0, g_xyyyzzz_0_yyyzzz_0, g_xyyyzzz_0_yyzzzz_0, g_xyyyzzz_0_yzzzzz_0, g_xyyyzzz_0_zzzzzz_0, g_yyyzzz_0_xxxxx_1, g_yyyzzz_0_xxxxxx_1, g_yyyzzz_0_xxxxxy_1, g_yyyzzz_0_xxxxxz_1, g_yyyzzz_0_xxxxy_1, g_yyyzzz_0_xxxxyy_1, g_yyyzzz_0_xxxxyz_1, g_yyyzzz_0_xxxxz_1, g_yyyzzz_0_xxxxzz_1, g_yyyzzz_0_xxxyy_1, g_yyyzzz_0_xxxyyy_1, g_yyyzzz_0_xxxyyz_1, g_yyyzzz_0_xxxyz_1, g_yyyzzz_0_xxxyzz_1, g_yyyzzz_0_xxxzz_1, g_yyyzzz_0_xxxzzz_1, g_yyyzzz_0_xxyyy_1, g_yyyzzz_0_xxyyyy_1, g_yyyzzz_0_xxyyyz_1, g_yyyzzz_0_xxyyz_1, g_yyyzzz_0_xxyyzz_1, g_yyyzzz_0_xxyzz_1, g_yyyzzz_0_xxyzzz_1, g_yyyzzz_0_xxzzz_1, g_yyyzzz_0_xxzzzz_1, g_yyyzzz_0_xyyyy_1, g_yyyzzz_0_xyyyyy_1, g_yyyzzz_0_xyyyyz_1, g_yyyzzz_0_xyyyz_1, g_yyyzzz_0_xyyyzz_1, g_yyyzzz_0_xyyzz_1, g_yyyzzz_0_xyyzzz_1, g_yyyzzz_0_xyzzz_1, g_yyyzzz_0_xyzzzz_1, g_yyyzzz_0_xzzzz_1, g_yyyzzz_0_xzzzzz_1, g_yyyzzz_0_yyyyy_1, g_yyyzzz_0_yyyyyy_1, g_yyyzzz_0_yyyyyz_1, g_yyyzzz_0_yyyyz_1, g_yyyzzz_0_yyyyzz_1, g_yyyzzz_0_yyyzz_1, g_yyyzzz_0_yyyzzz_1, g_yyyzzz_0_yyzzz_1, g_yyyzzz_0_yyzzzz_1, g_yyyzzz_0_yzzzz_1, g_yyyzzz_0_yzzzzz_1, g_yyyzzz_0_zzzzz_1, g_yyyzzz_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzzz_0_xxxxxx_0[i] = 6.0 * g_yyyzzz_0_xxxxx_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxx_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxy_0[i] = 5.0 * g_yyyzzz_0_xxxxy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxz_0[i] = 5.0 * g_yyyzzz_0_xxxxz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxyy_0[i] = 4.0 * g_yyyzzz_0_xxxyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxyz_0[i] = 4.0 * g_yyyzzz_0_xxxyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxzz_0[i] = 4.0 * g_yyyzzz_0_xxxzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyyy_0[i] = 3.0 * g_yyyzzz_0_xxyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyyz_0[i] = 3.0 * g_yyyzzz_0_xxyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyzz_0[i] = 3.0 * g_yyyzzz_0_xxyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxzzz_0[i] = 3.0 * g_yyyzzz_0_xxzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyyy_0[i] = 2.0 * g_yyyzzz_0_xyyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyyz_0[i] = 2.0 * g_yyyzzz_0_xyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyzz_0[i] = 2.0 * g_yyyzzz_0_xyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyzzz_0[i] = 2.0 * g_yyyzzz_0_xyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxzzzz_0[i] = 2.0 * g_yyyzzz_0_xzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyyy_0[i] = g_yyyzzz_0_yyyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyyz_0[i] = g_yyyzzz_0_yyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyzz_0[i] = g_yyyzzz_0_yyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyzzz_0[i] = g_yyyzzz_0_yyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyzzzz_0[i] = g_yyyzzz_0_yzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xzzzzz_0[i] = g_yyyzzz_0_zzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyyy_0[i] = g_yyyzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyyz_0[i] = g_yyyzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyzz_0[i] = g_yyyzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyzzz_0[i] = g_yyyzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyzzzz_0[i] = g_yyyzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yzzzzz_0[i] = g_yyyzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_zzzzzz_0[i] = g_yyyzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 700-728 components of targeted buffer : KSI

    auto g_xyyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 700);

    auto g_xyyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 701);

    auto g_xyyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 702);

    auto g_xyyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 703);

    auto g_xyyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 704);

    auto g_xyyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 705);

    auto g_xyyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 706);

    auto g_xyyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 707);

    auto g_xyyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 708);

    auto g_xyyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 709);

    auto g_xyyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 710);

    auto g_xyyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 711);

    auto g_xyyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 712);

    auto g_xyyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 713);

    auto g_xyyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 714);

    auto g_xyyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 715);

    auto g_xyyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 716);

    auto g_xyyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 717);

    auto g_xyyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 718);

    auto g_xyyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 719);

    auto g_xyyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 720);

    auto g_xyyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 721);

    auto g_xyyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 722);

    auto g_xyyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 723);

    auto g_xyyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 724);

    auto g_xyyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 725);

    auto g_xyyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 726);

    auto g_xyyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 727);

    #pragma omp simd aligned(g_xyyzzzz_0_xxxxxx_0, g_xyyzzzz_0_xxxxxy_0, g_xyyzzzz_0_xxxxxz_0, g_xyyzzzz_0_xxxxyy_0, g_xyyzzzz_0_xxxxyz_0, g_xyyzzzz_0_xxxxzz_0, g_xyyzzzz_0_xxxyyy_0, g_xyyzzzz_0_xxxyyz_0, g_xyyzzzz_0_xxxyzz_0, g_xyyzzzz_0_xxxzzz_0, g_xyyzzzz_0_xxyyyy_0, g_xyyzzzz_0_xxyyyz_0, g_xyyzzzz_0_xxyyzz_0, g_xyyzzzz_0_xxyzzz_0, g_xyyzzzz_0_xxzzzz_0, g_xyyzzzz_0_xyyyyy_0, g_xyyzzzz_0_xyyyyz_0, g_xyyzzzz_0_xyyyzz_0, g_xyyzzzz_0_xyyzzz_0, g_xyyzzzz_0_xyzzzz_0, g_xyyzzzz_0_xzzzzz_0, g_xyyzzzz_0_yyyyyy_0, g_xyyzzzz_0_yyyyyz_0, g_xyyzzzz_0_yyyyzz_0, g_xyyzzzz_0_yyyzzz_0, g_xyyzzzz_0_yyzzzz_0, g_xyyzzzz_0_yzzzzz_0, g_xyyzzzz_0_zzzzzz_0, g_yyzzzz_0_xxxxx_1, g_yyzzzz_0_xxxxxx_1, g_yyzzzz_0_xxxxxy_1, g_yyzzzz_0_xxxxxz_1, g_yyzzzz_0_xxxxy_1, g_yyzzzz_0_xxxxyy_1, g_yyzzzz_0_xxxxyz_1, g_yyzzzz_0_xxxxz_1, g_yyzzzz_0_xxxxzz_1, g_yyzzzz_0_xxxyy_1, g_yyzzzz_0_xxxyyy_1, g_yyzzzz_0_xxxyyz_1, g_yyzzzz_0_xxxyz_1, g_yyzzzz_0_xxxyzz_1, g_yyzzzz_0_xxxzz_1, g_yyzzzz_0_xxxzzz_1, g_yyzzzz_0_xxyyy_1, g_yyzzzz_0_xxyyyy_1, g_yyzzzz_0_xxyyyz_1, g_yyzzzz_0_xxyyz_1, g_yyzzzz_0_xxyyzz_1, g_yyzzzz_0_xxyzz_1, g_yyzzzz_0_xxyzzz_1, g_yyzzzz_0_xxzzz_1, g_yyzzzz_0_xxzzzz_1, g_yyzzzz_0_xyyyy_1, g_yyzzzz_0_xyyyyy_1, g_yyzzzz_0_xyyyyz_1, g_yyzzzz_0_xyyyz_1, g_yyzzzz_0_xyyyzz_1, g_yyzzzz_0_xyyzz_1, g_yyzzzz_0_xyyzzz_1, g_yyzzzz_0_xyzzz_1, g_yyzzzz_0_xyzzzz_1, g_yyzzzz_0_xzzzz_1, g_yyzzzz_0_xzzzzz_1, g_yyzzzz_0_yyyyy_1, g_yyzzzz_0_yyyyyy_1, g_yyzzzz_0_yyyyyz_1, g_yyzzzz_0_yyyyz_1, g_yyzzzz_0_yyyyzz_1, g_yyzzzz_0_yyyzz_1, g_yyzzzz_0_yyyzzz_1, g_yyzzzz_0_yyzzz_1, g_yyzzzz_0_yyzzzz_1, g_yyzzzz_0_yzzzz_1, g_yyzzzz_0_yzzzzz_1, g_yyzzzz_0_zzzzz_1, g_yyzzzz_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzzz_0_xxxxxx_0[i] = 6.0 * g_yyzzzz_0_xxxxx_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxx_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxy_0[i] = 5.0 * g_yyzzzz_0_xxxxy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxz_0[i] = 5.0 * g_yyzzzz_0_xxxxz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxyy_0[i] = 4.0 * g_yyzzzz_0_xxxyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxyz_0[i] = 4.0 * g_yyzzzz_0_xxxyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxzz_0[i] = 4.0 * g_yyzzzz_0_xxxzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyyy_0[i] = 3.0 * g_yyzzzz_0_xxyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyyz_0[i] = 3.0 * g_yyzzzz_0_xxyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyzz_0[i] = 3.0 * g_yyzzzz_0_xxyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxzzz_0[i] = 3.0 * g_yyzzzz_0_xxzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyyy_0[i] = 2.0 * g_yyzzzz_0_xyyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyyz_0[i] = 2.0 * g_yyzzzz_0_xyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyzz_0[i] = 2.0 * g_yyzzzz_0_xyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyzzz_0[i] = 2.0 * g_yyzzzz_0_xyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxzzzz_0[i] = 2.0 * g_yyzzzz_0_xzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyyy_0[i] = g_yyzzzz_0_yyyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyyz_0[i] = g_yyzzzz_0_yyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyzz_0[i] = g_yyzzzz_0_yyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyzzz_0[i] = g_yyzzzz_0_yyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyzzzz_0[i] = g_yyzzzz_0_yzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xzzzzz_0[i] = g_yyzzzz_0_zzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyyy_0[i] = g_yyzzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyyz_0[i] = g_yyzzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyzz_0[i] = g_yyzzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyzzz_0[i] = g_yyzzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyzzzz_0[i] = g_yyzzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yzzzzz_0[i] = g_yyzzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_zzzzzz_0[i] = g_yyzzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 728-756 components of targeted buffer : KSI

    auto g_xyzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 728);

    auto g_xyzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 729);

    auto g_xyzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 730);

    auto g_xyzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 731);

    auto g_xyzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 732);

    auto g_xyzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 733);

    auto g_xyzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 734);

    auto g_xyzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 735);

    auto g_xyzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 736);

    auto g_xyzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 737);

    auto g_xyzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 738);

    auto g_xyzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 739);

    auto g_xyzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 740);

    auto g_xyzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 741);

    auto g_xyzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 742);

    auto g_xyzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 743);

    auto g_xyzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 744);

    auto g_xyzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 745);

    auto g_xyzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 746);

    auto g_xyzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 747);

    auto g_xyzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 748);

    auto g_xyzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 749);

    auto g_xyzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 750);

    auto g_xyzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 751);

    auto g_xyzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 752);

    auto g_xyzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 753);

    auto g_xyzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 754);

    auto g_xyzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 755);

    #pragma omp simd aligned(g_xyzzzzz_0_xxxxxx_0, g_xyzzzzz_0_xxxxxy_0, g_xyzzzzz_0_xxxxxz_0, g_xyzzzzz_0_xxxxyy_0, g_xyzzzzz_0_xxxxyz_0, g_xyzzzzz_0_xxxxzz_0, g_xyzzzzz_0_xxxyyy_0, g_xyzzzzz_0_xxxyyz_0, g_xyzzzzz_0_xxxyzz_0, g_xyzzzzz_0_xxxzzz_0, g_xyzzzzz_0_xxyyyy_0, g_xyzzzzz_0_xxyyyz_0, g_xyzzzzz_0_xxyyzz_0, g_xyzzzzz_0_xxyzzz_0, g_xyzzzzz_0_xxzzzz_0, g_xyzzzzz_0_xyyyyy_0, g_xyzzzzz_0_xyyyyz_0, g_xyzzzzz_0_xyyyzz_0, g_xyzzzzz_0_xyyzzz_0, g_xyzzzzz_0_xyzzzz_0, g_xyzzzzz_0_xzzzzz_0, g_xyzzzzz_0_yyyyyy_0, g_xyzzzzz_0_yyyyyz_0, g_xyzzzzz_0_yyyyzz_0, g_xyzzzzz_0_yyyzzz_0, g_xyzzzzz_0_yyzzzz_0, g_xyzzzzz_0_yzzzzz_0, g_xyzzzzz_0_zzzzzz_0, g_xzzzzz_0_xxxxxx_1, g_xzzzzz_0_xxxxxz_1, g_xzzzzz_0_xxxxzz_1, g_xzzzzz_0_xxxzzz_1, g_xzzzzz_0_xxzzzz_1, g_xzzzzz_0_xzzzzz_1, g_yzzzzz_0_xxxxxy_1, g_yzzzzz_0_xxxxy_1, g_yzzzzz_0_xxxxyy_1, g_yzzzzz_0_xxxxyz_1, g_yzzzzz_0_xxxyy_1, g_yzzzzz_0_xxxyyy_1, g_yzzzzz_0_xxxyyz_1, g_yzzzzz_0_xxxyz_1, g_yzzzzz_0_xxxyzz_1, g_yzzzzz_0_xxyyy_1, g_yzzzzz_0_xxyyyy_1, g_yzzzzz_0_xxyyyz_1, g_yzzzzz_0_xxyyz_1, g_yzzzzz_0_xxyyzz_1, g_yzzzzz_0_xxyzz_1, g_yzzzzz_0_xxyzzz_1, g_yzzzzz_0_xyyyy_1, g_yzzzzz_0_xyyyyy_1, g_yzzzzz_0_xyyyyz_1, g_yzzzzz_0_xyyyz_1, g_yzzzzz_0_xyyyzz_1, g_yzzzzz_0_xyyzz_1, g_yzzzzz_0_xyyzzz_1, g_yzzzzz_0_xyzzz_1, g_yzzzzz_0_xyzzzz_1, g_yzzzzz_0_yyyyy_1, g_yzzzzz_0_yyyyyy_1, g_yzzzzz_0_yyyyyz_1, g_yzzzzz_0_yyyyz_1, g_yzzzzz_0_yyyyzz_1, g_yzzzzz_0_yyyzz_1, g_yzzzzz_0_yyyzzz_1, g_yzzzzz_0_yyzzz_1, g_yzzzzz_0_yyzzzz_1, g_yzzzzz_0_yzzzz_1, g_yzzzzz_0_yzzzzz_1, g_yzzzzz_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzzz_0_xxxxxx_0[i] = g_xzzzzz_0_xxxxxx_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxxxy_0[i] = 5.0 * g_yzzzzz_0_xxxxy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxxz_0[i] = g_xzzzzz_0_xxxxxz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxxyy_0[i] = 4.0 * g_yzzzzz_0_xxxyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxyz_0[i] = 4.0 * g_yzzzzz_0_xxxyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxzz_0[i] = g_xzzzzz_0_xxxxzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxyyy_0[i] = 3.0 * g_yzzzzz_0_xxyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxyyz_0[i] = 3.0 * g_yzzzzz_0_xxyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxyzz_0[i] = 3.0 * g_yzzzzz_0_xxyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxzzz_0[i] = g_xzzzzz_0_xxxzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxyyyy_0[i] = 2.0 * g_yzzzzz_0_xyyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyyyz_0[i] = 2.0 * g_yzzzzz_0_xyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyyzz_0[i] = 2.0 * g_yzzzzz_0_xyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyzzz_0[i] = 2.0 * g_yzzzzz_0_xyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxzzzz_0[i] = g_xzzzzz_0_xxzzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xyyyyy_0[i] = g_yzzzzz_0_yyyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyyyz_0[i] = g_yzzzzz_0_yyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyyzz_0[i] = g_yzzzzz_0_yyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyzzz_0[i] = g_yzzzzz_0_yyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyzzzz_0[i] = g_yzzzzz_0_yzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xzzzzz_0[i] = g_xzzzzz_0_xzzzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_yyyyyy_0[i] = g_yzzzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyyyz_0[i] = g_yzzzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyyzz_0[i] = g_yzzzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyzzz_0[i] = g_yzzzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyzzzz_0[i] = g_yzzzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yzzzzz_0[i] = g_yzzzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_zzzzzz_0[i] = g_yzzzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 756-784 components of targeted buffer : KSI

    auto g_xzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 756);

    auto g_xzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 757);

    auto g_xzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 758);

    auto g_xzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 759);

    auto g_xzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 760);

    auto g_xzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 761);

    auto g_xzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 762);

    auto g_xzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 763);

    auto g_xzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 764);

    auto g_xzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 765);

    auto g_xzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 766);

    auto g_xzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 767);

    auto g_xzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 768);

    auto g_xzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 769);

    auto g_xzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 770);

    auto g_xzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 771);

    auto g_xzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 772);

    auto g_xzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 773);

    auto g_xzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 774);

    auto g_xzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 775);

    auto g_xzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 776);

    auto g_xzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 777);

    auto g_xzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 778);

    auto g_xzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 779);

    auto g_xzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 780);

    auto g_xzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 781);

    auto g_xzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 782);

    auto g_xzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 783);

    #pragma omp simd aligned(g_xzzzzzz_0_xxxxxx_0, g_xzzzzzz_0_xxxxxy_0, g_xzzzzzz_0_xxxxxz_0, g_xzzzzzz_0_xxxxyy_0, g_xzzzzzz_0_xxxxyz_0, g_xzzzzzz_0_xxxxzz_0, g_xzzzzzz_0_xxxyyy_0, g_xzzzzzz_0_xxxyyz_0, g_xzzzzzz_0_xxxyzz_0, g_xzzzzzz_0_xxxzzz_0, g_xzzzzzz_0_xxyyyy_0, g_xzzzzzz_0_xxyyyz_0, g_xzzzzzz_0_xxyyzz_0, g_xzzzzzz_0_xxyzzz_0, g_xzzzzzz_0_xxzzzz_0, g_xzzzzzz_0_xyyyyy_0, g_xzzzzzz_0_xyyyyz_0, g_xzzzzzz_0_xyyyzz_0, g_xzzzzzz_0_xyyzzz_0, g_xzzzzzz_0_xyzzzz_0, g_xzzzzzz_0_xzzzzz_0, g_xzzzzzz_0_yyyyyy_0, g_xzzzzzz_0_yyyyyz_0, g_xzzzzzz_0_yyyyzz_0, g_xzzzzzz_0_yyyzzz_0, g_xzzzzzz_0_yyzzzz_0, g_xzzzzzz_0_yzzzzz_0, g_xzzzzzz_0_zzzzzz_0, g_zzzzzz_0_xxxxx_1, g_zzzzzz_0_xxxxxx_1, g_zzzzzz_0_xxxxxy_1, g_zzzzzz_0_xxxxxz_1, g_zzzzzz_0_xxxxy_1, g_zzzzzz_0_xxxxyy_1, g_zzzzzz_0_xxxxyz_1, g_zzzzzz_0_xxxxz_1, g_zzzzzz_0_xxxxzz_1, g_zzzzzz_0_xxxyy_1, g_zzzzzz_0_xxxyyy_1, g_zzzzzz_0_xxxyyz_1, g_zzzzzz_0_xxxyz_1, g_zzzzzz_0_xxxyzz_1, g_zzzzzz_0_xxxzz_1, g_zzzzzz_0_xxxzzz_1, g_zzzzzz_0_xxyyy_1, g_zzzzzz_0_xxyyyy_1, g_zzzzzz_0_xxyyyz_1, g_zzzzzz_0_xxyyz_1, g_zzzzzz_0_xxyyzz_1, g_zzzzzz_0_xxyzz_1, g_zzzzzz_0_xxyzzz_1, g_zzzzzz_0_xxzzz_1, g_zzzzzz_0_xxzzzz_1, g_zzzzzz_0_xyyyy_1, g_zzzzzz_0_xyyyyy_1, g_zzzzzz_0_xyyyyz_1, g_zzzzzz_0_xyyyz_1, g_zzzzzz_0_xyyyzz_1, g_zzzzzz_0_xyyzz_1, g_zzzzzz_0_xyyzzz_1, g_zzzzzz_0_xyzzz_1, g_zzzzzz_0_xyzzzz_1, g_zzzzzz_0_xzzzz_1, g_zzzzzz_0_xzzzzz_1, g_zzzzzz_0_yyyyy_1, g_zzzzzz_0_yyyyyy_1, g_zzzzzz_0_yyyyyz_1, g_zzzzzz_0_yyyyz_1, g_zzzzzz_0_yyyyzz_1, g_zzzzzz_0_yyyzz_1, g_zzzzzz_0_yyyzzz_1, g_zzzzzz_0_yyzzz_1, g_zzzzzz_0_yyzzzz_1, g_zzzzzz_0_yzzzz_1, g_zzzzzz_0_yzzzzz_1, g_zzzzzz_0_zzzzz_1, g_zzzzzz_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzzz_0_xxxxxx_0[i] = 6.0 * g_zzzzzz_0_xxxxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxx_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxy_0[i] = 5.0 * g_zzzzzz_0_xxxxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxz_0[i] = 5.0 * g_zzzzzz_0_xxxxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxyy_0[i] = 4.0 * g_zzzzzz_0_xxxyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxyz_0[i] = 4.0 * g_zzzzzz_0_xxxyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxzz_0[i] = 4.0 * g_zzzzzz_0_xxxzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyyy_0[i] = 3.0 * g_zzzzzz_0_xxyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyyz_0[i] = 3.0 * g_zzzzzz_0_xxyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyzz_0[i] = 3.0 * g_zzzzzz_0_xxyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxzzz_0[i] = 3.0 * g_zzzzzz_0_xxzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyyy_0[i] = 2.0 * g_zzzzzz_0_xyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyyz_0[i] = 2.0 * g_zzzzzz_0_xyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyzz_0[i] = 2.0 * g_zzzzzz_0_xyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyzzz_0[i] = 2.0 * g_zzzzzz_0_xyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxzzzz_0[i] = 2.0 * g_zzzzzz_0_xzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyyy_0[i] = g_zzzzzz_0_yyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyyz_0[i] = g_zzzzzz_0_yyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyzz_0[i] = g_zzzzzz_0_yyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyzzz_0[i] = g_zzzzzz_0_yyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyzzzz_0[i] = g_zzzzzz_0_yzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xzzzzz_0[i] = g_zzzzzz_0_zzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyyy_0[i] = g_zzzzzz_0_yyyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyyz_0[i] = g_zzzzzz_0_yyyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyzz_0[i] = g_zzzzzz_0_yyyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyzzz_0[i] = g_zzzzzz_0_yyyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyzzzz_0[i] = g_zzzzzz_0_yyzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yzzzzz_0[i] = g_zzzzzz_0_yzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_zzzzzz_0[i] = g_zzzzzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 784-812 components of targeted buffer : KSI

    auto g_yyyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 784);

    auto g_yyyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 785);

    auto g_yyyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 786);

    auto g_yyyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 787);

    auto g_yyyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 788);

    auto g_yyyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 789);

    auto g_yyyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 790);

    auto g_yyyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 791);

    auto g_yyyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 792);

    auto g_yyyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 793);

    auto g_yyyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 794);

    auto g_yyyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 795);

    auto g_yyyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 796);

    auto g_yyyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 797);

    auto g_yyyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 798);

    auto g_yyyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 799);

    auto g_yyyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 800);

    auto g_yyyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 801);

    auto g_yyyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 802);

    auto g_yyyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 803);

    auto g_yyyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 804);

    auto g_yyyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 805);

    auto g_yyyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 806);

    auto g_yyyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 807);

    auto g_yyyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 808);

    auto g_yyyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 809);

    auto g_yyyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 810);

    auto g_yyyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 811);

    #pragma omp simd aligned(g_yyyyy_0_xxxxxx_0, g_yyyyy_0_xxxxxx_1, g_yyyyy_0_xxxxxy_0, g_yyyyy_0_xxxxxy_1, g_yyyyy_0_xxxxxz_0, g_yyyyy_0_xxxxxz_1, g_yyyyy_0_xxxxyy_0, g_yyyyy_0_xxxxyy_1, g_yyyyy_0_xxxxyz_0, g_yyyyy_0_xxxxyz_1, g_yyyyy_0_xxxxzz_0, g_yyyyy_0_xxxxzz_1, g_yyyyy_0_xxxyyy_0, g_yyyyy_0_xxxyyy_1, g_yyyyy_0_xxxyyz_0, g_yyyyy_0_xxxyyz_1, g_yyyyy_0_xxxyzz_0, g_yyyyy_0_xxxyzz_1, g_yyyyy_0_xxxzzz_0, g_yyyyy_0_xxxzzz_1, g_yyyyy_0_xxyyyy_0, g_yyyyy_0_xxyyyy_1, g_yyyyy_0_xxyyyz_0, g_yyyyy_0_xxyyyz_1, g_yyyyy_0_xxyyzz_0, g_yyyyy_0_xxyyzz_1, g_yyyyy_0_xxyzzz_0, g_yyyyy_0_xxyzzz_1, g_yyyyy_0_xxzzzz_0, g_yyyyy_0_xxzzzz_1, g_yyyyy_0_xyyyyy_0, g_yyyyy_0_xyyyyy_1, g_yyyyy_0_xyyyyz_0, g_yyyyy_0_xyyyyz_1, g_yyyyy_0_xyyyzz_0, g_yyyyy_0_xyyyzz_1, g_yyyyy_0_xyyzzz_0, g_yyyyy_0_xyyzzz_1, g_yyyyy_0_xyzzzz_0, g_yyyyy_0_xyzzzz_1, g_yyyyy_0_xzzzzz_0, g_yyyyy_0_xzzzzz_1, g_yyyyy_0_yyyyyy_0, g_yyyyy_0_yyyyyy_1, g_yyyyy_0_yyyyyz_0, g_yyyyy_0_yyyyyz_1, g_yyyyy_0_yyyyzz_0, g_yyyyy_0_yyyyzz_1, g_yyyyy_0_yyyzzz_0, g_yyyyy_0_yyyzzz_1, g_yyyyy_0_yyzzzz_0, g_yyyyy_0_yyzzzz_1, g_yyyyy_0_yzzzzz_0, g_yyyyy_0_yzzzzz_1, g_yyyyy_0_zzzzzz_0, g_yyyyy_0_zzzzzz_1, g_yyyyyy_0_xxxxx_1, g_yyyyyy_0_xxxxxx_1, g_yyyyyy_0_xxxxxy_1, g_yyyyyy_0_xxxxxz_1, g_yyyyyy_0_xxxxy_1, g_yyyyyy_0_xxxxyy_1, g_yyyyyy_0_xxxxyz_1, g_yyyyyy_0_xxxxz_1, g_yyyyyy_0_xxxxzz_1, g_yyyyyy_0_xxxyy_1, g_yyyyyy_0_xxxyyy_1, g_yyyyyy_0_xxxyyz_1, g_yyyyyy_0_xxxyz_1, g_yyyyyy_0_xxxyzz_1, g_yyyyyy_0_xxxzz_1, g_yyyyyy_0_xxxzzz_1, g_yyyyyy_0_xxyyy_1, g_yyyyyy_0_xxyyyy_1, g_yyyyyy_0_xxyyyz_1, g_yyyyyy_0_xxyyz_1, g_yyyyyy_0_xxyyzz_1, g_yyyyyy_0_xxyzz_1, g_yyyyyy_0_xxyzzz_1, g_yyyyyy_0_xxzzz_1, g_yyyyyy_0_xxzzzz_1, g_yyyyyy_0_xyyyy_1, g_yyyyyy_0_xyyyyy_1, g_yyyyyy_0_xyyyyz_1, g_yyyyyy_0_xyyyz_1, g_yyyyyy_0_xyyyzz_1, g_yyyyyy_0_xyyzz_1, g_yyyyyy_0_xyyzzz_1, g_yyyyyy_0_xyzzz_1, g_yyyyyy_0_xyzzzz_1, g_yyyyyy_0_xzzzz_1, g_yyyyyy_0_xzzzzz_1, g_yyyyyy_0_yyyyy_1, g_yyyyyy_0_yyyyyy_1, g_yyyyyy_0_yyyyyz_1, g_yyyyyy_0_yyyyz_1, g_yyyyyy_0_yyyyzz_1, g_yyyyyy_0_yyyzz_1, g_yyyyyy_0_yyyzzz_1, g_yyyyyy_0_yyzzz_1, g_yyyyyy_0_yyzzzz_1, g_yyyyyy_0_yzzzz_1, g_yyyyyy_0_yzzzzz_1, g_yyyyyy_0_zzzzz_1, g_yyyyyy_0_zzzzzz_1, g_yyyyyyy_0_xxxxxx_0, g_yyyyyyy_0_xxxxxy_0, g_yyyyyyy_0_xxxxxz_0, g_yyyyyyy_0_xxxxyy_0, g_yyyyyyy_0_xxxxyz_0, g_yyyyyyy_0_xxxxzz_0, g_yyyyyyy_0_xxxyyy_0, g_yyyyyyy_0_xxxyyz_0, g_yyyyyyy_0_xxxyzz_0, g_yyyyyyy_0_xxxzzz_0, g_yyyyyyy_0_xxyyyy_0, g_yyyyyyy_0_xxyyyz_0, g_yyyyyyy_0_xxyyzz_0, g_yyyyyyy_0_xxyzzz_0, g_yyyyyyy_0_xxzzzz_0, g_yyyyyyy_0_xyyyyy_0, g_yyyyyyy_0_xyyyyz_0, g_yyyyyyy_0_xyyyzz_0, g_yyyyyyy_0_xyyzzz_0, g_yyyyyyy_0_xyzzzz_0, g_yyyyyyy_0_xzzzzz_0, g_yyyyyyy_0_yyyyyy_0, g_yyyyyyy_0_yyyyyz_0, g_yyyyyyy_0_yyyyzz_0, g_yyyyyyy_0_yyyzzz_0, g_yyyyyyy_0_yyzzzz_0, g_yyyyyyy_0_yzzzzz_0, g_yyyyyyy_0_zzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyyy_0_xxxxxx_0[i] = 6.0 * g_yyyyy_0_xxxxxx_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxx_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxx_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxy_0[i] = 6.0 * g_yyyyy_0_xxxxxy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxy_1[i] * fz_be_0 + g_yyyyyy_0_xxxxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxz_0[i] = 6.0 * g_yyyyy_0_xxxxxz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxyy_0[i] = 6.0 * g_yyyyy_0_xxxxyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxyy_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxxxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxyz_0[i] = 6.0 * g_yyyyy_0_xxxxyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxyz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxzz_0[i] = 6.0 * g_yyyyy_0_xxxxzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyyy_0[i] = 6.0 * g_yyyyy_0_xxxyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xxxyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyyz_0[i] = 6.0 * g_yyyyy_0_xxxyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxxyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyzz_0[i] = 6.0 * g_yyyyy_0_xxxyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxzzz_0[i] = 6.0 * g_yyyyy_0_xxxzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyyy_0[i] = 6.0 * g_yyyyy_0_xxyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyyy_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_xxyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyyz_0[i] = 6.0 * g_yyyyy_0_xxyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xxyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyzz_0[i] = 6.0 * g_yyyyy_0_xxyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyzzz_0[i] = 6.0 * g_yyyyy_0_xxyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxzzzz_0[i] = 6.0 * g_yyyyy_0_xxzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyyy_0[i] = 6.0 * g_yyyyy_0_xyyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyyy_1[i] * fz_be_0 + 5.0 * g_yyyyyy_0_xyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyyz_0[i] = 6.0 * g_yyyyy_0_xyyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_xyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyzz_0[i] = 6.0 * g_yyyyy_0_xyyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyzzz_0[i] = 6.0 * g_yyyyy_0_xyyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyzzzz_0[i] = 6.0 * g_yyyyy_0_xyzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xzzzzz_0[i] = 6.0 * g_yyyyy_0_xzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyyy_0[i] = 6.0 * g_yyyyy_0_yyyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyyy_1[i] * fz_be_0 + 6.0 * g_yyyyyy_0_yyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyyz_0[i] = 6.0 * g_yyyyy_0_yyyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyyyy_0_yyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyzz_0[i] = 6.0 * g_yyyyy_0_yyyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_yyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyzzz_0[i] = 6.0 * g_yyyyy_0_yyyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_yyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyzzzz_0[i] = 6.0 * g_yyyyy_0_yyzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_yzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yzzzzz_0[i] = 6.0 * g_yyyyy_0_yzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_zzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_zzzzzz_0[i] = 6.0 * g_yyyyy_0_zzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_zzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 812-840 components of targeted buffer : KSI

    auto g_yyyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 812);

    auto g_yyyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 813);

    auto g_yyyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 814);

    auto g_yyyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 815);

    auto g_yyyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 816);

    auto g_yyyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 817);

    auto g_yyyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 818);

    auto g_yyyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 819);

    auto g_yyyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 820);

    auto g_yyyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 821);

    auto g_yyyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 822);

    auto g_yyyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 823);

    auto g_yyyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 824);

    auto g_yyyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 825);

    auto g_yyyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 826);

    auto g_yyyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 827);

    auto g_yyyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 828);

    auto g_yyyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 829);

    auto g_yyyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 830);

    auto g_yyyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 831);

    auto g_yyyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 832);

    auto g_yyyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 833);

    auto g_yyyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 834);

    auto g_yyyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 835);

    auto g_yyyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 836);

    auto g_yyyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 837);

    auto g_yyyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 838);

    auto g_yyyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 839);

    #pragma omp simd aligned(g_yyyyyy_0_xxxxx_1, g_yyyyyy_0_xxxxxx_1, g_yyyyyy_0_xxxxxy_1, g_yyyyyy_0_xxxxxz_1, g_yyyyyy_0_xxxxy_1, g_yyyyyy_0_xxxxyy_1, g_yyyyyy_0_xxxxyz_1, g_yyyyyy_0_xxxxz_1, g_yyyyyy_0_xxxxzz_1, g_yyyyyy_0_xxxyy_1, g_yyyyyy_0_xxxyyy_1, g_yyyyyy_0_xxxyyz_1, g_yyyyyy_0_xxxyz_1, g_yyyyyy_0_xxxyzz_1, g_yyyyyy_0_xxxzz_1, g_yyyyyy_0_xxxzzz_1, g_yyyyyy_0_xxyyy_1, g_yyyyyy_0_xxyyyy_1, g_yyyyyy_0_xxyyyz_1, g_yyyyyy_0_xxyyz_1, g_yyyyyy_0_xxyyzz_1, g_yyyyyy_0_xxyzz_1, g_yyyyyy_0_xxyzzz_1, g_yyyyyy_0_xxzzz_1, g_yyyyyy_0_xxzzzz_1, g_yyyyyy_0_xyyyy_1, g_yyyyyy_0_xyyyyy_1, g_yyyyyy_0_xyyyyz_1, g_yyyyyy_0_xyyyz_1, g_yyyyyy_0_xyyyzz_1, g_yyyyyy_0_xyyzz_1, g_yyyyyy_0_xyyzzz_1, g_yyyyyy_0_xyzzz_1, g_yyyyyy_0_xyzzzz_1, g_yyyyyy_0_xzzzz_1, g_yyyyyy_0_xzzzzz_1, g_yyyyyy_0_yyyyy_1, g_yyyyyy_0_yyyyyy_1, g_yyyyyy_0_yyyyyz_1, g_yyyyyy_0_yyyyz_1, g_yyyyyy_0_yyyyzz_1, g_yyyyyy_0_yyyzz_1, g_yyyyyy_0_yyyzzz_1, g_yyyyyy_0_yyzzz_1, g_yyyyyy_0_yyzzzz_1, g_yyyyyy_0_yzzzz_1, g_yyyyyy_0_yzzzzz_1, g_yyyyyy_0_zzzzz_1, g_yyyyyy_0_zzzzzz_1, g_yyyyyyz_0_xxxxxx_0, g_yyyyyyz_0_xxxxxy_0, g_yyyyyyz_0_xxxxxz_0, g_yyyyyyz_0_xxxxyy_0, g_yyyyyyz_0_xxxxyz_0, g_yyyyyyz_0_xxxxzz_0, g_yyyyyyz_0_xxxyyy_0, g_yyyyyyz_0_xxxyyz_0, g_yyyyyyz_0_xxxyzz_0, g_yyyyyyz_0_xxxzzz_0, g_yyyyyyz_0_xxyyyy_0, g_yyyyyyz_0_xxyyyz_0, g_yyyyyyz_0_xxyyzz_0, g_yyyyyyz_0_xxyzzz_0, g_yyyyyyz_0_xxzzzz_0, g_yyyyyyz_0_xyyyyy_0, g_yyyyyyz_0_xyyyyz_0, g_yyyyyyz_0_xyyyzz_0, g_yyyyyyz_0_xyyzzz_0, g_yyyyyyz_0_xyzzzz_0, g_yyyyyyz_0_xzzzzz_0, g_yyyyyyz_0_yyyyyy_0, g_yyyyyyz_0_yyyyyz_0, g_yyyyyyz_0_yyyyzz_0, g_yyyyyyz_0_yyyzzz_0, g_yyyyyyz_0_yyzzzz_0, g_yyyyyyz_0_yzzzzz_0, g_yyyyyyz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyyz_0_xxxxxx_0[i] = g_yyyyyy_0_xxxxxx_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxy_0[i] = g_yyyyyy_0_xxxxxy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxz_0[i] = g_yyyyyy_0_xxxxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxyy_0[i] = g_yyyyyy_0_xxxxyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxyz_0[i] = g_yyyyyy_0_xxxxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxzz_0[i] = 2.0 * g_yyyyyy_0_xxxxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyyy_0[i] = g_yyyyyy_0_xxxyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyyz_0[i] = g_yyyyyy_0_xxxyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyzz_0[i] = 2.0 * g_yyyyyy_0_xxxyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxzzz_0[i] = 3.0 * g_yyyyyy_0_xxxzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyyy_0[i] = g_yyyyyy_0_xxyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyyz_0[i] = g_yyyyyy_0_xxyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyzz_0[i] = 2.0 * g_yyyyyy_0_xxyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyzzz_0[i] = 3.0 * g_yyyyyy_0_xxyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxzzzz_0[i] = 4.0 * g_yyyyyy_0_xxzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyyy_0[i] = g_yyyyyy_0_xyyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyyz_0[i] = g_yyyyyy_0_xyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyzz_0[i] = 2.0 * g_yyyyyy_0_xyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyzzz_0[i] = 3.0 * g_yyyyyy_0_xyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyzzzz_0[i] = 4.0 * g_yyyyyy_0_xyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xzzzzz_0[i] = 5.0 * g_yyyyyy_0_xzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyyy_0[i] = g_yyyyyy_0_yyyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyyz_0[i] = g_yyyyyy_0_yyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyzz_0[i] = 2.0 * g_yyyyyy_0_yyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyzzz_0[i] = 3.0 * g_yyyyyy_0_yyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyzzzz_0[i] = 4.0 * g_yyyyyy_0_yyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yzzzzz_0[i] = 5.0 * g_yyyyyy_0_yzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_zzzzzz_0[i] = 6.0 * g_yyyyyy_0_zzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 840-868 components of targeted buffer : KSI

    auto g_yyyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 840);

    auto g_yyyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 841);

    auto g_yyyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 842);

    auto g_yyyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 843);

    auto g_yyyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 844);

    auto g_yyyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 845);

    auto g_yyyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 846);

    auto g_yyyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 847);

    auto g_yyyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 848);

    auto g_yyyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 849);

    auto g_yyyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 850);

    auto g_yyyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 851);

    auto g_yyyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 852);

    auto g_yyyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 853);

    auto g_yyyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 854);

    auto g_yyyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 855);

    auto g_yyyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 856);

    auto g_yyyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 857);

    auto g_yyyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 858);

    auto g_yyyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 859);

    auto g_yyyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 860);

    auto g_yyyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 861);

    auto g_yyyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 862);

    auto g_yyyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 863);

    auto g_yyyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 864);

    auto g_yyyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 865);

    auto g_yyyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 866);

    auto g_yyyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 867);

    #pragma omp simd aligned(g_yyyyy_0_xxxxxy_0, g_yyyyy_0_xxxxxy_1, g_yyyyy_0_xxxxyy_0, g_yyyyy_0_xxxxyy_1, g_yyyyy_0_xxxyyy_0, g_yyyyy_0_xxxyyy_1, g_yyyyy_0_xxyyyy_0, g_yyyyy_0_xxyyyy_1, g_yyyyy_0_xyyyyy_0, g_yyyyy_0_xyyyyy_1, g_yyyyy_0_yyyyyy_0, g_yyyyy_0_yyyyyy_1, g_yyyyyz_0_xxxxxy_1, g_yyyyyz_0_xxxxyy_1, g_yyyyyz_0_xxxyyy_1, g_yyyyyz_0_xxyyyy_1, g_yyyyyz_0_xyyyyy_1, g_yyyyyz_0_yyyyyy_1, g_yyyyyzz_0_xxxxxx_0, g_yyyyyzz_0_xxxxxy_0, g_yyyyyzz_0_xxxxxz_0, g_yyyyyzz_0_xxxxyy_0, g_yyyyyzz_0_xxxxyz_0, g_yyyyyzz_0_xxxxzz_0, g_yyyyyzz_0_xxxyyy_0, g_yyyyyzz_0_xxxyyz_0, g_yyyyyzz_0_xxxyzz_0, g_yyyyyzz_0_xxxzzz_0, g_yyyyyzz_0_xxyyyy_0, g_yyyyyzz_0_xxyyyz_0, g_yyyyyzz_0_xxyyzz_0, g_yyyyyzz_0_xxyzzz_0, g_yyyyyzz_0_xxzzzz_0, g_yyyyyzz_0_xyyyyy_0, g_yyyyyzz_0_xyyyyz_0, g_yyyyyzz_0_xyyyzz_0, g_yyyyyzz_0_xyyzzz_0, g_yyyyyzz_0_xyzzzz_0, g_yyyyyzz_0_xzzzzz_0, g_yyyyyzz_0_yyyyyy_0, g_yyyyyzz_0_yyyyyz_0, g_yyyyyzz_0_yyyyzz_0, g_yyyyyzz_0_yyyzzz_0, g_yyyyyzz_0_yyzzzz_0, g_yyyyyzz_0_yzzzzz_0, g_yyyyyzz_0_zzzzzz_0, g_yyyyzz_0_xxxxxx_1, g_yyyyzz_0_xxxxxz_1, g_yyyyzz_0_xxxxyz_1, g_yyyyzz_0_xxxxz_1, g_yyyyzz_0_xxxxzz_1, g_yyyyzz_0_xxxyyz_1, g_yyyyzz_0_xxxyz_1, g_yyyyzz_0_xxxyzz_1, g_yyyyzz_0_xxxzz_1, g_yyyyzz_0_xxxzzz_1, g_yyyyzz_0_xxyyyz_1, g_yyyyzz_0_xxyyz_1, g_yyyyzz_0_xxyyzz_1, g_yyyyzz_0_xxyzz_1, g_yyyyzz_0_xxyzzz_1, g_yyyyzz_0_xxzzz_1, g_yyyyzz_0_xxzzzz_1, g_yyyyzz_0_xyyyyz_1, g_yyyyzz_0_xyyyz_1, g_yyyyzz_0_xyyyzz_1, g_yyyyzz_0_xyyzz_1, g_yyyyzz_0_xyyzzz_1, g_yyyyzz_0_xyzzz_1, g_yyyyzz_0_xyzzzz_1, g_yyyyzz_0_xzzzz_1, g_yyyyzz_0_xzzzzz_1, g_yyyyzz_0_yyyyyz_1, g_yyyyzz_0_yyyyz_1, g_yyyyzz_0_yyyyzz_1, g_yyyyzz_0_yyyzz_1, g_yyyyzz_0_yyyzzz_1, g_yyyyzz_0_yyzzz_1, g_yyyyzz_0_yyzzzz_1, g_yyyyzz_0_yzzzz_1, g_yyyyzz_0_yzzzzz_1, g_yyyyzz_0_zzzzz_1, g_yyyyzz_0_zzzzzz_1, g_yyyzz_0_xxxxxx_0, g_yyyzz_0_xxxxxx_1, g_yyyzz_0_xxxxxz_0, g_yyyzz_0_xxxxxz_1, g_yyyzz_0_xxxxyz_0, g_yyyzz_0_xxxxyz_1, g_yyyzz_0_xxxxzz_0, g_yyyzz_0_xxxxzz_1, g_yyyzz_0_xxxyyz_0, g_yyyzz_0_xxxyyz_1, g_yyyzz_0_xxxyzz_0, g_yyyzz_0_xxxyzz_1, g_yyyzz_0_xxxzzz_0, g_yyyzz_0_xxxzzz_1, g_yyyzz_0_xxyyyz_0, g_yyyzz_0_xxyyyz_1, g_yyyzz_0_xxyyzz_0, g_yyyzz_0_xxyyzz_1, g_yyyzz_0_xxyzzz_0, g_yyyzz_0_xxyzzz_1, g_yyyzz_0_xxzzzz_0, g_yyyzz_0_xxzzzz_1, g_yyyzz_0_xyyyyz_0, g_yyyzz_0_xyyyyz_1, g_yyyzz_0_xyyyzz_0, g_yyyzz_0_xyyyzz_1, g_yyyzz_0_xyyzzz_0, g_yyyzz_0_xyyzzz_1, g_yyyzz_0_xyzzzz_0, g_yyyzz_0_xyzzzz_1, g_yyyzz_0_xzzzzz_0, g_yyyzz_0_xzzzzz_1, g_yyyzz_0_yyyyyz_0, g_yyyzz_0_yyyyyz_1, g_yyyzz_0_yyyyzz_0, g_yyyzz_0_yyyyzz_1, g_yyyzz_0_yyyzzz_0, g_yyyzz_0_yyyzzz_1, g_yyyzz_0_yyzzzz_0, g_yyyzz_0_yyzzzz_1, g_yyyzz_0_yzzzzz_0, g_yyyzz_0_yzzzzz_1, g_yyyzz_0_zzzzzz_0, g_yyyzz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyzz_0_xxxxxx_0[i] = 4.0 * g_yyyzz_0_xxxxxx_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxx_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxx_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxxy_0[i] = g_yyyyy_0_xxxxxy_0[i] * fbe_0 - g_yyyyy_0_xxxxxy_1[i] * fz_be_0 + g_yyyyyz_0_xxxxxy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxxxz_0[i] = 4.0 * g_yyyzz_0_xxxxxz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxyy_0[i] = g_yyyyy_0_xxxxyy_0[i] * fbe_0 - g_yyyyy_0_xxxxyy_1[i] * fz_be_0 + g_yyyyyz_0_xxxxyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxxyz_0[i] = 4.0 * g_yyyzz_0_xxxxyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxyz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxzz_0[i] = 4.0 * g_yyyzz_0_xxxxzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxyyy_0[i] = g_yyyyy_0_xxxyyy_0[i] * fbe_0 - g_yyyyy_0_xxxyyy_1[i] * fz_be_0 + g_yyyyyz_0_xxxyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxyyz_0[i] = 4.0 * g_yyyzz_0_xxxyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xxxyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxyzz_0[i] = 4.0 * g_yyyzz_0_xxxyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxyzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxzzz_0[i] = 4.0 * g_yyyzz_0_xxxzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyyyy_0[i] = g_yyyyy_0_xxyyyy_0[i] * fbe_0 - g_yyyyy_0_xxyyyy_1[i] * fz_be_0 + g_yyyyyz_0_xxyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxyyyz_0[i] = 4.0 * g_yyyzz_0_xxyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_xxyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyyzz_0[i] = 4.0 * g_yyyzz_0_xxyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xxyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyzzz_0[i] = 4.0 * g_yyyzz_0_xxyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxzzzz_0[i] = 4.0 * g_yyyzz_0_xxzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyyyy_0[i] = g_yyyyy_0_xyyyyy_0[i] * fbe_0 - g_yyyyy_0_xyyyyy_1[i] * fz_be_0 + g_yyyyyz_0_xyyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xyyyyz_0[i] = 4.0 * g_yyyzz_0_xyyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_0_xyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyyzz_0[i] = 4.0 * g_yyyzz_0_xyyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_xyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyzzz_0[i] = 4.0 * g_yyyzz_0_xyyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyzzzz_0[i] = 4.0 * g_yyyzz_0_xyzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xzzzzz_0[i] = 4.0 * g_yyyzz_0_xzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyyyy_0[i] = g_yyyyy_0_yyyyyy_0[i] * fbe_0 - g_yyyyy_0_yyyyyy_1[i] * fz_be_0 + g_yyyyyz_0_yyyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_yyyyyz_0[i] = 4.0 * g_yyyzz_0_yyyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyyzz_0_yyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyyzz_0[i] = 4.0 * g_yyyzz_0_yyyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_0_yyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyzzz_0[i] = 4.0 * g_yyyzz_0_yyyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_yyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyzzzz_0[i] = 4.0 * g_yyyzz_0_yyzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_yzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yzzzzz_0[i] = 4.0 * g_yyyzz_0_yzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_zzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_zzzzzz_0[i] = 4.0 * g_yyyzz_0_zzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_zzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 868-896 components of targeted buffer : KSI

    auto g_yyyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 868);

    auto g_yyyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 869);

    auto g_yyyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 870);

    auto g_yyyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 871);

    auto g_yyyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 872);

    auto g_yyyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 873);

    auto g_yyyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 874);

    auto g_yyyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 875);

    auto g_yyyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 876);

    auto g_yyyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 877);

    auto g_yyyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 878);

    auto g_yyyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 879);

    auto g_yyyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 880);

    auto g_yyyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 881);

    auto g_yyyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 882);

    auto g_yyyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 883);

    auto g_yyyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 884);

    auto g_yyyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 885);

    auto g_yyyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 886);

    auto g_yyyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 887);

    auto g_yyyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 888);

    auto g_yyyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 889);

    auto g_yyyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 890);

    auto g_yyyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 891);

    auto g_yyyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 892);

    auto g_yyyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 893);

    auto g_yyyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 894);

    auto g_yyyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 895);

    #pragma omp simd aligned(g_yyyyz_0_xxxxxy_0, g_yyyyz_0_xxxxxy_1, g_yyyyz_0_xxxxyy_0, g_yyyyz_0_xxxxyy_1, g_yyyyz_0_xxxyyy_0, g_yyyyz_0_xxxyyy_1, g_yyyyz_0_xxyyyy_0, g_yyyyz_0_xxyyyy_1, g_yyyyz_0_xyyyyy_0, g_yyyyz_0_xyyyyy_1, g_yyyyz_0_yyyyyy_0, g_yyyyz_0_yyyyyy_1, g_yyyyzz_0_xxxxxy_1, g_yyyyzz_0_xxxxyy_1, g_yyyyzz_0_xxxyyy_1, g_yyyyzz_0_xxyyyy_1, g_yyyyzz_0_xyyyyy_1, g_yyyyzz_0_yyyyyy_1, g_yyyyzzz_0_xxxxxx_0, g_yyyyzzz_0_xxxxxy_0, g_yyyyzzz_0_xxxxxz_0, g_yyyyzzz_0_xxxxyy_0, g_yyyyzzz_0_xxxxyz_0, g_yyyyzzz_0_xxxxzz_0, g_yyyyzzz_0_xxxyyy_0, g_yyyyzzz_0_xxxyyz_0, g_yyyyzzz_0_xxxyzz_0, g_yyyyzzz_0_xxxzzz_0, g_yyyyzzz_0_xxyyyy_0, g_yyyyzzz_0_xxyyyz_0, g_yyyyzzz_0_xxyyzz_0, g_yyyyzzz_0_xxyzzz_0, g_yyyyzzz_0_xxzzzz_0, g_yyyyzzz_0_xyyyyy_0, g_yyyyzzz_0_xyyyyz_0, g_yyyyzzz_0_xyyyzz_0, g_yyyyzzz_0_xyyzzz_0, g_yyyyzzz_0_xyzzzz_0, g_yyyyzzz_0_xzzzzz_0, g_yyyyzzz_0_yyyyyy_0, g_yyyyzzz_0_yyyyyz_0, g_yyyyzzz_0_yyyyzz_0, g_yyyyzzz_0_yyyzzz_0, g_yyyyzzz_0_yyzzzz_0, g_yyyyzzz_0_yzzzzz_0, g_yyyyzzz_0_zzzzzz_0, g_yyyzzz_0_xxxxxx_1, g_yyyzzz_0_xxxxxz_1, g_yyyzzz_0_xxxxyz_1, g_yyyzzz_0_xxxxz_1, g_yyyzzz_0_xxxxzz_1, g_yyyzzz_0_xxxyyz_1, g_yyyzzz_0_xxxyz_1, g_yyyzzz_0_xxxyzz_1, g_yyyzzz_0_xxxzz_1, g_yyyzzz_0_xxxzzz_1, g_yyyzzz_0_xxyyyz_1, g_yyyzzz_0_xxyyz_1, g_yyyzzz_0_xxyyzz_1, g_yyyzzz_0_xxyzz_1, g_yyyzzz_0_xxyzzz_1, g_yyyzzz_0_xxzzz_1, g_yyyzzz_0_xxzzzz_1, g_yyyzzz_0_xyyyyz_1, g_yyyzzz_0_xyyyz_1, g_yyyzzz_0_xyyyzz_1, g_yyyzzz_0_xyyzz_1, g_yyyzzz_0_xyyzzz_1, g_yyyzzz_0_xyzzz_1, g_yyyzzz_0_xyzzzz_1, g_yyyzzz_0_xzzzz_1, g_yyyzzz_0_xzzzzz_1, g_yyyzzz_0_yyyyyz_1, g_yyyzzz_0_yyyyz_1, g_yyyzzz_0_yyyyzz_1, g_yyyzzz_0_yyyzz_1, g_yyyzzz_0_yyyzzz_1, g_yyyzzz_0_yyzzz_1, g_yyyzzz_0_yyzzzz_1, g_yyyzzz_0_yzzzz_1, g_yyyzzz_0_yzzzzz_1, g_yyyzzz_0_zzzzz_1, g_yyyzzz_0_zzzzzz_1, g_yyzzz_0_xxxxxx_0, g_yyzzz_0_xxxxxx_1, g_yyzzz_0_xxxxxz_0, g_yyzzz_0_xxxxxz_1, g_yyzzz_0_xxxxyz_0, g_yyzzz_0_xxxxyz_1, g_yyzzz_0_xxxxzz_0, g_yyzzz_0_xxxxzz_1, g_yyzzz_0_xxxyyz_0, g_yyzzz_0_xxxyyz_1, g_yyzzz_0_xxxyzz_0, g_yyzzz_0_xxxyzz_1, g_yyzzz_0_xxxzzz_0, g_yyzzz_0_xxxzzz_1, g_yyzzz_0_xxyyyz_0, g_yyzzz_0_xxyyyz_1, g_yyzzz_0_xxyyzz_0, g_yyzzz_0_xxyyzz_1, g_yyzzz_0_xxyzzz_0, g_yyzzz_0_xxyzzz_1, g_yyzzz_0_xxzzzz_0, g_yyzzz_0_xxzzzz_1, g_yyzzz_0_xyyyyz_0, g_yyzzz_0_xyyyyz_1, g_yyzzz_0_xyyyzz_0, g_yyzzz_0_xyyyzz_1, g_yyzzz_0_xyyzzz_0, g_yyzzz_0_xyyzzz_1, g_yyzzz_0_xyzzzz_0, g_yyzzz_0_xyzzzz_1, g_yyzzz_0_xzzzzz_0, g_yyzzz_0_xzzzzz_1, g_yyzzz_0_yyyyyz_0, g_yyzzz_0_yyyyyz_1, g_yyzzz_0_yyyyzz_0, g_yyzzz_0_yyyyzz_1, g_yyzzz_0_yyyzzz_0, g_yyzzz_0_yyyzzz_1, g_yyzzz_0_yyzzzz_0, g_yyzzz_0_yyzzzz_1, g_yyzzz_0_yzzzzz_0, g_yyzzz_0_yzzzzz_1, g_yyzzz_0_zzzzzz_0, g_yyzzz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzzz_0_xxxxxx_0[i] = 3.0 * g_yyzzz_0_xxxxxx_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxx_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxx_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxxy_0[i] = 2.0 * g_yyyyz_0_xxxxxy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxxxy_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxxxz_0[i] = 3.0 * g_yyzzz_0_xxxxxz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxyy_0[i] = 2.0 * g_yyyyz_0_xxxxyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxxyy_1[i] * fz_be_0 + g_yyyyzz_0_xxxxyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxxyz_0[i] = 3.0 * g_yyzzz_0_xxxxyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxyz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxzz_0[i] = 3.0 * g_yyzzz_0_xxxxzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxyyy_0[i] = 2.0 * g_yyyyz_0_xxxyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxyyy_1[i] * fz_be_0 + g_yyyyzz_0_xxxyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxyyz_0[i] = 3.0 * g_yyzzz_0_xxxyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xxxyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxyzz_0[i] = 3.0 * g_yyzzz_0_xxxyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxyzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxzzz_0[i] = 3.0 * g_yyzzz_0_xxxzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyyyy_0[i] = 2.0 * g_yyyyz_0_xxyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxyyyy_1[i] * fz_be_0 + g_yyyyzz_0_xxyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxyyyz_0[i] = 3.0 * g_yyzzz_0_xxyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_xxyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyyzz_0[i] = 3.0 * g_yyzzz_0_xxyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xxyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyzzz_0[i] = 3.0 * g_yyzzz_0_xxyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxzzzz_0[i] = 3.0 * g_yyzzz_0_xxzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyyyy_0[i] = 2.0 * g_yyyyz_0_xyyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xyyyyy_1[i] * fz_be_0 + g_yyyyzz_0_xyyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xyyyyz_0[i] = 3.0 * g_yyzzz_0_xyyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_0_xyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyyzz_0[i] = 3.0 * g_yyzzz_0_xyyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_xyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyzzz_0[i] = 3.0 * g_yyzzz_0_xyyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyzzzz_0[i] = 3.0 * g_yyzzz_0_xyzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xzzzzz_0[i] = 3.0 * g_yyzzz_0_xzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyyyy_0[i] = 2.0 * g_yyyyz_0_yyyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_yyyyyy_1[i] * fz_be_0 + g_yyyyzz_0_yyyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_yyyyyz_0[i] = 3.0 * g_yyzzz_0_yyyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyzzz_0_yyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyyzz_0[i] = 3.0 * g_yyzzz_0_yyyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_0_yyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyzzz_0[i] = 3.0 * g_yyzzz_0_yyyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_yyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyzzzz_0[i] = 3.0 * g_yyzzz_0_yyzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_yzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yzzzzz_0[i] = 3.0 * g_yyzzz_0_yzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_zzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_zzzzzz_0[i] = 3.0 * g_yyzzz_0_zzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_zzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 896-924 components of targeted buffer : KSI

    auto g_yyyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 896);

    auto g_yyyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 897);

    auto g_yyyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 898);

    auto g_yyyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 899);

    auto g_yyyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 900);

    auto g_yyyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 901);

    auto g_yyyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 902);

    auto g_yyyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 903);

    auto g_yyyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 904);

    auto g_yyyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 905);

    auto g_yyyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 906);

    auto g_yyyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 907);

    auto g_yyyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 908);

    auto g_yyyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 909);

    auto g_yyyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 910);

    auto g_yyyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 911);

    auto g_yyyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 912);

    auto g_yyyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 913);

    auto g_yyyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 914);

    auto g_yyyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 915);

    auto g_yyyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 916);

    auto g_yyyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 917);

    auto g_yyyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 918);

    auto g_yyyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 919);

    auto g_yyyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 920);

    auto g_yyyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 921);

    auto g_yyyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 922);

    auto g_yyyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 923);

    #pragma omp simd aligned(g_yyyzz_0_xxxxxy_0, g_yyyzz_0_xxxxxy_1, g_yyyzz_0_xxxxyy_0, g_yyyzz_0_xxxxyy_1, g_yyyzz_0_xxxyyy_0, g_yyyzz_0_xxxyyy_1, g_yyyzz_0_xxyyyy_0, g_yyyzz_0_xxyyyy_1, g_yyyzz_0_xyyyyy_0, g_yyyzz_0_xyyyyy_1, g_yyyzz_0_yyyyyy_0, g_yyyzz_0_yyyyyy_1, g_yyyzzz_0_xxxxxy_1, g_yyyzzz_0_xxxxyy_1, g_yyyzzz_0_xxxyyy_1, g_yyyzzz_0_xxyyyy_1, g_yyyzzz_0_xyyyyy_1, g_yyyzzz_0_yyyyyy_1, g_yyyzzzz_0_xxxxxx_0, g_yyyzzzz_0_xxxxxy_0, g_yyyzzzz_0_xxxxxz_0, g_yyyzzzz_0_xxxxyy_0, g_yyyzzzz_0_xxxxyz_0, g_yyyzzzz_0_xxxxzz_0, g_yyyzzzz_0_xxxyyy_0, g_yyyzzzz_0_xxxyyz_0, g_yyyzzzz_0_xxxyzz_0, g_yyyzzzz_0_xxxzzz_0, g_yyyzzzz_0_xxyyyy_0, g_yyyzzzz_0_xxyyyz_0, g_yyyzzzz_0_xxyyzz_0, g_yyyzzzz_0_xxyzzz_0, g_yyyzzzz_0_xxzzzz_0, g_yyyzzzz_0_xyyyyy_0, g_yyyzzzz_0_xyyyyz_0, g_yyyzzzz_0_xyyyzz_0, g_yyyzzzz_0_xyyzzz_0, g_yyyzzzz_0_xyzzzz_0, g_yyyzzzz_0_xzzzzz_0, g_yyyzzzz_0_yyyyyy_0, g_yyyzzzz_0_yyyyyz_0, g_yyyzzzz_0_yyyyzz_0, g_yyyzzzz_0_yyyzzz_0, g_yyyzzzz_0_yyzzzz_0, g_yyyzzzz_0_yzzzzz_0, g_yyyzzzz_0_zzzzzz_0, g_yyzzzz_0_xxxxxx_1, g_yyzzzz_0_xxxxxz_1, g_yyzzzz_0_xxxxyz_1, g_yyzzzz_0_xxxxz_1, g_yyzzzz_0_xxxxzz_1, g_yyzzzz_0_xxxyyz_1, g_yyzzzz_0_xxxyz_1, g_yyzzzz_0_xxxyzz_1, g_yyzzzz_0_xxxzz_1, g_yyzzzz_0_xxxzzz_1, g_yyzzzz_0_xxyyyz_1, g_yyzzzz_0_xxyyz_1, g_yyzzzz_0_xxyyzz_1, g_yyzzzz_0_xxyzz_1, g_yyzzzz_0_xxyzzz_1, g_yyzzzz_0_xxzzz_1, g_yyzzzz_0_xxzzzz_1, g_yyzzzz_0_xyyyyz_1, g_yyzzzz_0_xyyyz_1, g_yyzzzz_0_xyyyzz_1, g_yyzzzz_0_xyyzz_1, g_yyzzzz_0_xyyzzz_1, g_yyzzzz_0_xyzzz_1, g_yyzzzz_0_xyzzzz_1, g_yyzzzz_0_xzzzz_1, g_yyzzzz_0_xzzzzz_1, g_yyzzzz_0_yyyyyz_1, g_yyzzzz_0_yyyyz_1, g_yyzzzz_0_yyyyzz_1, g_yyzzzz_0_yyyzz_1, g_yyzzzz_0_yyyzzz_1, g_yyzzzz_0_yyzzz_1, g_yyzzzz_0_yyzzzz_1, g_yyzzzz_0_yzzzz_1, g_yyzzzz_0_yzzzzz_1, g_yyzzzz_0_zzzzz_1, g_yyzzzz_0_zzzzzz_1, g_yzzzz_0_xxxxxx_0, g_yzzzz_0_xxxxxx_1, g_yzzzz_0_xxxxxz_0, g_yzzzz_0_xxxxxz_1, g_yzzzz_0_xxxxyz_0, g_yzzzz_0_xxxxyz_1, g_yzzzz_0_xxxxzz_0, g_yzzzz_0_xxxxzz_1, g_yzzzz_0_xxxyyz_0, g_yzzzz_0_xxxyyz_1, g_yzzzz_0_xxxyzz_0, g_yzzzz_0_xxxyzz_1, g_yzzzz_0_xxxzzz_0, g_yzzzz_0_xxxzzz_1, g_yzzzz_0_xxyyyz_0, g_yzzzz_0_xxyyyz_1, g_yzzzz_0_xxyyzz_0, g_yzzzz_0_xxyyzz_1, g_yzzzz_0_xxyzzz_0, g_yzzzz_0_xxyzzz_1, g_yzzzz_0_xxzzzz_0, g_yzzzz_0_xxzzzz_1, g_yzzzz_0_xyyyyz_0, g_yzzzz_0_xyyyyz_1, g_yzzzz_0_xyyyzz_0, g_yzzzz_0_xyyyzz_1, g_yzzzz_0_xyyzzz_0, g_yzzzz_0_xyyzzz_1, g_yzzzz_0_xyzzzz_0, g_yzzzz_0_xyzzzz_1, g_yzzzz_0_xzzzzz_0, g_yzzzz_0_xzzzzz_1, g_yzzzz_0_yyyyyz_0, g_yzzzz_0_yyyyyz_1, g_yzzzz_0_yyyyzz_0, g_yzzzz_0_yyyyzz_1, g_yzzzz_0_yyyzzz_0, g_yzzzz_0_yyyzzz_1, g_yzzzz_0_yyzzzz_0, g_yzzzz_0_yyzzzz_1, g_yzzzz_0_yzzzzz_0, g_yzzzz_0_yzzzzz_1, g_yzzzz_0_zzzzzz_0, g_yzzzz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzzz_0_xxxxxx_0[i] = 2.0 * g_yzzzz_0_xxxxxx_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxx_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxx_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxxy_0[i] = 3.0 * g_yyyzz_0_xxxxxy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxxxy_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxxxz_0[i] = 2.0 * g_yzzzz_0_xxxxxz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxyy_0[i] = 3.0 * g_yyyzz_0_xxxxyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxxyy_1[i] * fz_be_0 + g_yyyzzz_0_xxxxyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxxyz_0[i] = 2.0 * g_yzzzz_0_xxxxyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxyz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxzz_0[i] = 2.0 * g_yzzzz_0_xxxxzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxyyy_0[i] = 3.0 * g_yyyzz_0_xxxyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxyyy_1[i] * fz_be_0 + g_yyyzzz_0_xxxyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxyyz_0[i] = 2.0 * g_yzzzz_0_xxxyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xxxyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxyzz_0[i] = 2.0 * g_yzzzz_0_xxxyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxyzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxzzz_0[i] = 2.0 * g_yzzzz_0_xxxzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyyyy_0[i] = 3.0 * g_yyyzz_0_xxyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxyyyy_1[i] * fz_be_0 + g_yyyzzz_0_xxyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxyyyz_0[i] = 2.0 * g_yzzzz_0_xxyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_xxyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyyzz_0[i] = 2.0 * g_yzzzz_0_xxyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xxyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyzzz_0[i] = 2.0 * g_yzzzz_0_xxyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxzzzz_0[i] = 2.0 * g_yzzzz_0_xxzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyyyy_0[i] = 3.0 * g_yyyzz_0_xyyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xyyyyy_1[i] * fz_be_0 + g_yyyzzz_0_xyyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xyyyyz_0[i] = 2.0 * g_yzzzz_0_xyyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_0_xyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyyzz_0[i] = 2.0 * g_yzzzz_0_xyyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_xyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyzzz_0[i] = 2.0 * g_yzzzz_0_xyyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyzzzz_0[i] = 2.0 * g_yzzzz_0_xyzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xzzzzz_0[i] = 2.0 * g_yzzzz_0_xzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyyyy_0[i] = 3.0 * g_yyyzz_0_yyyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_yyyyyy_1[i] * fz_be_0 + g_yyyzzz_0_yyyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_yyyyyz_0[i] = 2.0 * g_yzzzz_0_yyyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyzzzz_0_yyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyyzz_0[i] = 2.0 * g_yzzzz_0_yyyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_0_yyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyzzz_0[i] = 2.0 * g_yzzzz_0_yyyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_yyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyzzzz_0[i] = 2.0 * g_yzzzz_0_yyzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_yzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yzzzzz_0[i] = 2.0 * g_yzzzz_0_yzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_zzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_zzzzzz_0[i] = 2.0 * g_yzzzz_0_zzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_zzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 924-952 components of targeted buffer : KSI

    auto g_yyzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 924);

    auto g_yyzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 925);

    auto g_yyzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 926);

    auto g_yyzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 927);

    auto g_yyzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 928);

    auto g_yyzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 929);

    auto g_yyzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 930);

    auto g_yyzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 931);

    auto g_yyzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 932);

    auto g_yyzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 933);

    auto g_yyzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 934);

    auto g_yyzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 935);

    auto g_yyzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 936);

    auto g_yyzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 937);

    auto g_yyzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 938);

    auto g_yyzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 939);

    auto g_yyzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 940);

    auto g_yyzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 941);

    auto g_yyzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 942);

    auto g_yyzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 943);

    auto g_yyzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 944);

    auto g_yyzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 945);

    auto g_yyzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 946);

    auto g_yyzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 947);

    auto g_yyzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 948);

    auto g_yyzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 949);

    auto g_yyzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 950);

    auto g_yyzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 951);

    #pragma omp simd aligned(g_yyzzz_0_xxxxxy_0, g_yyzzz_0_xxxxxy_1, g_yyzzz_0_xxxxyy_0, g_yyzzz_0_xxxxyy_1, g_yyzzz_0_xxxyyy_0, g_yyzzz_0_xxxyyy_1, g_yyzzz_0_xxyyyy_0, g_yyzzz_0_xxyyyy_1, g_yyzzz_0_xyyyyy_0, g_yyzzz_0_xyyyyy_1, g_yyzzz_0_yyyyyy_0, g_yyzzz_0_yyyyyy_1, g_yyzzzz_0_xxxxxy_1, g_yyzzzz_0_xxxxyy_1, g_yyzzzz_0_xxxyyy_1, g_yyzzzz_0_xxyyyy_1, g_yyzzzz_0_xyyyyy_1, g_yyzzzz_0_yyyyyy_1, g_yyzzzzz_0_xxxxxx_0, g_yyzzzzz_0_xxxxxy_0, g_yyzzzzz_0_xxxxxz_0, g_yyzzzzz_0_xxxxyy_0, g_yyzzzzz_0_xxxxyz_0, g_yyzzzzz_0_xxxxzz_0, g_yyzzzzz_0_xxxyyy_0, g_yyzzzzz_0_xxxyyz_0, g_yyzzzzz_0_xxxyzz_0, g_yyzzzzz_0_xxxzzz_0, g_yyzzzzz_0_xxyyyy_0, g_yyzzzzz_0_xxyyyz_0, g_yyzzzzz_0_xxyyzz_0, g_yyzzzzz_0_xxyzzz_0, g_yyzzzzz_0_xxzzzz_0, g_yyzzzzz_0_xyyyyy_0, g_yyzzzzz_0_xyyyyz_0, g_yyzzzzz_0_xyyyzz_0, g_yyzzzzz_0_xyyzzz_0, g_yyzzzzz_0_xyzzzz_0, g_yyzzzzz_0_xzzzzz_0, g_yyzzzzz_0_yyyyyy_0, g_yyzzzzz_0_yyyyyz_0, g_yyzzzzz_0_yyyyzz_0, g_yyzzzzz_0_yyyzzz_0, g_yyzzzzz_0_yyzzzz_0, g_yyzzzzz_0_yzzzzz_0, g_yyzzzzz_0_zzzzzz_0, g_yzzzzz_0_xxxxxx_1, g_yzzzzz_0_xxxxxz_1, g_yzzzzz_0_xxxxyz_1, g_yzzzzz_0_xxxxz_1, g_yzzzzz_0_xxxxzz_1, g_yzzzzz_0_xxxyyz_1, g_yzzzzz_0_xxxyz_1, g_yzzzzz_0_xxxyzz_1, g_yzzzzz_0_xxxzz_1, g_yzzzzz_0_xxxzzz_1, g_yzzzzz_0_xxyyyz_1, g_yzzzzz_0_xxyyz_1, g_yzzzzz_0_xxyyzz_1, g_yzzzzz_0_xxyzz_1, g_yzzzzz_0_xxyzzz_1, g_yzzzzz_0_xxzzz_1, g_yzzzzz_0_xxzzzz_1, g_yzzzzz_0_xyyyyz_1, g_yzzzzz_0_xyyyz_1, g_yzzzzz_0_xyyyzz_1, g_yzzzzz_0_xyyzz_1, g_yzzzzz_0_xyyzzz_1, g_yzzzzz_0_xyzzz_1, g_yzzzzz_0_xyzzzz_1, g_yzzzzz_0_xzzzz_1, g_yzzzzz_0_xzzzzz_1, g_yzzzzz_0_yyyyyz_1, g_yzzzzz_0_yyyyz_1, g_yzzzzz_0_yyyyzz_1, g_yzzzzz_0_yyyzz_1, g_yzzzzz_0_yyyzzz_1, g_yzzzzz_0_yyzzz_1, g_yzzzzz_0_yyzzzz_1, g_yzzzzz_0_yzzzz_1, g_yzzzzz_0_yzzzzz_1, g_yzzzzz_0_zzzzz_1, g_yzzzzz_0_zzzzzz_1, g_zzzzz_0_xxxxxx_0, g_zzzzz_0_xxxxxx_1, g_zzzzz_0_xxxxxz_0, g_zzzzz_0_xxxxxz_1, g_zzzzz_0_xxxxyz_0, g_zzzzz_0_xxxxyz_1, g_zzzzz_0_xxxxzz_0, g_zzzzz_0_xxxxzz_1, g_zzzzz_0_xxxyyz_0, g_zzzzz_0_xxxyyz_1, g_zzzzz_0_xxxyzz_0, g_zzzzz_0_xxxyzz_1, g_zzzzz_0_xxxzzz_0, g_zzzzz_0_xxxzzz_1, g_zzzzz_0_xxyyyz_0, g_zzzzz_0_xxyyyz_1, g_zzzzz_0_xxyyzz_0, g_zzzzz_0_xxyyzz_1, g_zzzzz_0_xxyzzz_0, g_zzzzz_0_xxyzzz_1, g_zzzzz_0_xxzzzz_0, g_zzzzz_0_xxzzzz_1, g_zzzzz_0_xyyyyz_0, g_zzzzz_0_xyyyyz_1, g_zzzzz_0_xyyyzz_0, g_zzzzz_0_xyyyzz_1, g_zzzzz_0_xyyzzz_0, g_zzzzz_0_xyyzzz_1, g_zzzzz_0_xyzzzz_0, g_zzzzz_0_xyzzzz_1, g_zzzzz_0_xzzzzz_0, g_zzzzz_0_xzzzzz_1, g_zzzzz_0_yyyyyz_0, g_zzzzz_0_yyyyyz_1, g_zzzzz_0_yyyyzz_0, g_zzzzz_0_yyyyzz_1, g_zzzzz_0_yyyzzz_0, g_zzzzz_0_yyyzzz_1, g_zzzzz_0_yyzzzz_0, g_zzzzz_0_yyzzzz_1, g_zzzzz_0_yzzzzz_0, g_zzzzz_0_yzzzzz_1, g_zzzzz_0_zzzzzz_0, g_zzzzz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzzz_0_xxxxxx_0[i] = g_zzzzz_0_xxxxxx_0[i] * fbe_0 - g_zzzzz_0_xxxxxx_1[i] * fz_be_0 + g_yzzzzz_0_xxxxxx_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxxy_0[i] = 4.0 * g_yyzzz_0_xxxxxy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxxxy_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxxxz_0[i] = g_zzzzz_0_xxxxxz_0[i] * fbe_0 - g_zzzzz_0_xxxxxz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxxz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxyy_0[i] = 4.0 * g_yyzzz_0_xxxxyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxxyy_1[i] * fz_be_0 + g_yyzzzz_0_xxxxyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxxyz_0[i] = g_zzzzz_0_xxxxyz_0[i] * fbe_0 - g_zzzzz_0_xxxxyz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxzz_0[i] = g_zzzzz_0_xxxxzz_0[i] * fbe_0 - g_zzzzz_0_xxxxzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxyyy_0[i] = 4.0 * g_yyzzz_0_xxxyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxyyy_1[i] * fz_be_0 + g_yyzzzz_0_xxxyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxyyz_0[i] = g_zzzzz_0_xxxyyz_0[i] * fbe_0 - g_zzzzz_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xxxyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxyzz_0[i] = g_zzzzz_0_xxxyzz_0[i] * fbe_0 - g_zzzzz_0_xxxyzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxzzz_0[i] = g_zzzzz_0_xxxzzz_0[i] * fbe_0 - g_zzzzz_0_xxxzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyyyy_0[i] = 4.0 * g_yyzzz_0_xxyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxyyyy_1[i] * fz_be_0 + g_yyzzzz_0_xxyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxyyyz_0[i] = g_zzzzz_0_xxyyyz_0[i] * fbe_0 - g_zzzzz_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_xxyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyyzz_0[i] = g_zzzzz_0_xxyyzz_0[i] * fbe_0 - g_zzzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xxyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyzzz_0[i] = g_zzzzz_0_xxyzzz_0[i] * fbe_0 - g_zzzzz_0_xxyzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxzzzz_0[i] = g_zzzzz_0_xxzzzz_0[i] * fbe_0 - g_zzzzz_0_xxzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyyyy_0[i] = 4.0 * g_yyzzz_0_xyyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xyyyyy_1[i] * fz_be_0 + g_yyzzzz_0_xyyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xyyyyz_0[i] = g_zzzzz_0_xyyyyz_0[i] * fbe_0 - g_zzzzz_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_0_xyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyyzz_0[i] = g_zzzzz_0_xyyyzz_0[i] * fbe_0 - g_zzzzz_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_xyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyzzz_0[i] = g_zzzzz_0_xyyzzz_0[i] * fbe_0 - g_zzzzz_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyzzzz_0[i] = g_zzzzz_0_xyzzzz_0[i] * fbe_0 - g_zzzzz_0_xyzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xzzzzz_0[i] = g_zzzzz_0_xzzzzz_0[i] * fbe_0 - g_zzzzz_0_xzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyyyy_0[i] = 4.0 * g_yyzzz_0_yyyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_yyyyyy_1[i] * fz_be_0 + g_yyzzzz_0_yyyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_yyyyyz_0[i] = g_zzzzz_0_yyyyyz_0[i] * fbe_0 - g_zzzzz_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yzzzzz_0_yyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyyzz_0[i] = g_zzzzz_0_yyyyzz_0[i] * fbe_0 - g_zzzzz_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_0_yyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyzzz_0[i] = g_zzzzz_0_yyyzzz_0[i] * fbe_0 - g_zzzzz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_yyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyzzzz_0[i] = g_zzzzz_0_yyzzzz_0[i] * fbe_0 - g_zzzzz_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_yzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yzzzzz_0[i] = g_zzzzz_0_yzzzzz_0[i] * fbe_0 - g_zzzzz_0_yzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_zzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_zzzzzz_0[i] = g_zzzzz_0_zzzzzz_0[i] * fbe_0 - g_zzzzz_0_zzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 952-980 components of targeted buffer : KSI

    auto g_yzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 952);

    auto g_yzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 953);

    auto g_yzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 954);

    auto g_yzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 955);

    auto g_yzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 956);

    auto g_yzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 957);

    auto g_yzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 958);

    auto g_yzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 959);

    auto g_yzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 960);

    auto g_yzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 961);

    auto g_yzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 962);

    auto g_yzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 963);

    auto g_yzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 964);

    auto g_yzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 965);

    auto g_yzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 966);

    auto g_yzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 967);

    auto g_yzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 968);

    auto g_yzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 969);

    auto g_yzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 970);

    auto g_yzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 971);

    auto g_yzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 972);

    auto g_yzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 973);

    auto g_yzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 974);

    auto g_yzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 975);

    auto g_yzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 976);

    auto g_yzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 977);

    auto g_yzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 978);

    auto g_yzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 979);

    #pragma omp simd aligned(g_yzzzzzz_0_xxxxxx_0, g_yzzzzzz_0_xxxxxy_0, g_yzzzzzz_0_xxxxxz_0, g_yzzzzzz_0_xxxxyy_0, g_yzzzzzz_0_xxxxyz_0, g_yzzzzzz_0_xxxxzz_0, g_yzzzzzz_0_xxxyyy_0, g_yzzzzzz_0_xxxyyz_0, g_yzzzzzz_0_xxxyzz_0, g_yzzzzzz_0_xxxzzz_0, g_yzzzzzz_0_xxyyyy_0, g_yzzzzzz_0_xxyyyz_0, g_yzzzzzz_0_xxyyzz_0, g_yzzzzzz_0_xxyzzz_0, g_yzzzzzz_0_xxzzzz_0, g_yzzzzzz_0_xyyyyy_0, g_yzzzzzz_0_xyyyyz_0, g_yzzzzzz_0_xyyyzz_0, g_yzzzzzz_0_xyyzzz_0, g_yzzzzzz_0_xyzzzz_0, g_yzzzzzz_0_xzzzzz_0, g_yzzzzzz_0_yyyyyy_0, g_yzzzzzz_0_yyyyyz_0, g_yzzzzzz_0_yyyyzz_0, g_yzzzzzz_0_yyyzzz_0, g_yzzzzzz_0_yyzzzz_0, g_yzzzzzz_0_yzzzzz_0, g_yzzzzzz_0_zzzzzz_0, g_zzzzzz_0_xxxxx_1, g_zzzzzz_0_xxxxxx_1, g_zzzzzz_0_xxxxxy_1, g_zzzzzz_0_xxxxxz_1, g_zzzzzz_0_xxxxy_1, g_zzzzzz_0_xxxxyy_1, g_zzzzzz_0_xxxxyz_1, g_zzzzzz_0_xxxxz_1, g_zzzzzz_0_xxxxzz_1, g_zzzzzz_0_xxxyy_1, g_zzzzzz_0_xxxyyy_1, g_zzzzzz_0_xxxyyz_1, g_zzzzzz_0_xxxyz_1, g_zzzzzz_0_xxxyzz_1, g_zzzzzz_0_xxxzz_1, g_zzzzzz_0_xxxzzz_1, g_zzzzzz_0_xxyyy_1, g_zzzzzz_0_xxyyyy_1, g_zzzzzz_0_xxyyyz_1, g_zzzzzz_0_xxyyz_1, g_zzzzzz_0_xxyyzz_1, g_zzzzzz_0_xxyzz_1, g_zzzzzz_0_xxyzzz_1, g_zzzzzz_0_xxzzz_1, g_zzzzzz_0_xxzzzz_1, g_zzzzzz_0_xyyyy_1, g_zzzzzz_0_xyyyyy_1, g_zzzzzz_0_xyyyyz_1, g_zzzzzz_0_xyyyz_1, g_zzzzzz_0_xyyyzz_1, g_zzzzzz_0_xyyzz_1, g_zzzzzz_0_xyyzzz_1, g_zzzzzz_0_xyzzz_1, g_zzzzzz_0_xyzzzz_1, g_zzzzzz_0_xzzzz_1, g_zzzzzz_0_xzzzzz_1, g_zzzzzz_0_yyyyy_1, g_zzzzzz_0_yyyyyy_1, g_zzzzzz_0_yyyyyz_1, g_zzzzzz_0_yyyyz_1, g_zzzzzz_0_yyyyzz_1, g_zzzzzz_0_yyyzz_1, g_zzzzzz_0_yyyzzz_1, g_zzzzzz_0_yyzzz_1, g_zzzzzz_0_yyzzzz_1, g_zzzzzz_0_yzzzz_1, g_zzzzzz_0_yzzzzz_1, g_zzzzzz_0_zzzzz_1, g_zzzzzz_0_zzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzzz_0_xxxxxx_0[i] = g_zzzzzz_0_xxxxxx_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxy_0[i] = g_zzzzzz_0_xxxxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxz_0[i] = g_zzzzzz_0_xxxxxz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxyy_0[i] = 2.0 * g_zzzzzz_0_xxxxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxyz_0[i] = g_zzzzzz_0_xxxxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxzz_0[i] = g_zzzzzz_0_xxxxzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyyy_0[i] = 3.0 * g_zzzzzz_0_xxxyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyyz_0[i] = 2.0 * g_zzzzzz_0_xxxyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyzz_0[i] = g_zzzzzz_0_xxxzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxzzz_0[i] = g_zzzzzz_0_xxxzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyyy_0[i] = 4.0 * g_zzzzzz_0_xxyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyyz_0[i] = 3.0 * g_zzzzzz_0_xxyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyzz_0[i] = 2.0 * g_zzzzzz_0_xxyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyzzz_0[i] = g_zzzzzz_0_xxzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxzzzz_0[i] = g_zzzzzz_0_xxzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyyy_0[i] = 5.0 * g_zzzzzz_0_xyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyyz_0[i] = 4.0 * g_zzzzzz_0_xyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyzz_0[i] = 3.0 * g_zzzzzz_0_xyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyzzz_0[i] = 2.0 * g_zzzzzz_0_xyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyzzzz_0[i] = g_zzzzzz_0_xzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xzzzzz_0[i] = g_zzzzzz_0_xzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyyy_0[i] = 6.0 * g_zzzzzz_0_yyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyyz_0[i] = 5.0 * g_zzzzzz_0_yyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyzz_0[i] = 4.0 * g_zzzzzz_0_yyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyzzz_0[i] = 3.0 * g_zzzzzz_0_yyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyzzzz_0[i] = 2.0 * g_zzzzzz_0_yzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yzzzzz_0[i] = g_zzzzzz_0_zzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_zzzzzz_0[i] = g_zzzzzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 980-1008 components of targeted buffer : KSI

    auto g_zzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ksi + 980);

    auto g_zzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ksi + 981);

    auto g_zzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ksi + 982);

    auto g_zzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ksi + 983);

    auto g_zzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ksi + 984);

    auto g_zzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ksi + 985);

    auto g_zzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ksi + 986);

    auto g_zzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ksi + 987);

    auto g_zzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ksi + 988);

    auto g_zzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ksi + 989);

    auto g_zzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ksi + 990);

    auto g_zzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ksi + 991);

    auto g_zzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ksi + 992);

    auto g_zzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ksi + 993);

    auto g_zzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ksi + 994);

    auto g_zzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 995);

    auto g_zzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 996);

    auto g_zzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 997);

    auto g_zzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 998);

    auto g_zzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 999);

    auto g_zzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 1000);

    auto g_zzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ksi + 1001);

    auto g_zzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ksi + 1002);

    auto g_zzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ksi + 1003);

    auto g_zzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ksi + 1004);

    auto g_zzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ksi + 1005);

    auto g_zzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 1006);

    auto g_zzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ksi + 1007);

    #pragma omp simd aligned(g_zzzzz_0_xxxxxx_0, g_zzzzz_0_xxxxxx_1, g_zzzzz_0_xxxxxy_0, g_zzzzz_0_xxxxxy_1, g_zzzzz_0_xxxxxz_0, g_zzzzz_0_xxxxxz_1, g_zzzzz_0_xxxxyy_0, g_zzzzz_0_xxxxyy_1, g_zzzzz_0_xxxxyz_0, g_zzzzz_0_xxxxyz_1, g_zzzzz_0_xxxxzz_0, g_zzzzz_0_xxxxzz_1, g_zzzzz_0_xxxyyy_0, g_zzzzz_0_xxxyyy_1, g_zzzzz_0_xxxyyz_0, g_zzzzz_0_xxxyyz_1, g_zzzzz_0_xxxyzz_0, g_zzzzz_0_xxxyzz_1, g_zzzzz_0_xxxzzz_0, g_zzzzz_0_xxxzzz_1, g_zzzzz_0_xxyyyy_0, g_zzzzz_0_xxyyyy_1, g_zzzzz_0_xxyyyz_0, g_zzzzz_0_xxyyyz_1, g_zzzzz_0_xxyyzz_0, g_zzzzz_0_xxyyzz_1, g_zzzzz_0_xxyzzz_0, g_zzzzz_0_xxyzzz_1, g_zzzzz_0_xxzzzz_0, g_zzzzz_0_xxzzzz_1, g_zzzzz_0_xyyyyy_0, g_zzzzz_0_xyyyyy_1, g_zzzzz_0_xyyyyz_0, g_zzzzz_0_xyyyyz_1, g_zzzzz_0_xyyyzz_0, g_zzzzz_0_xyyyzz_1, g_zzzzz_0_xyyzzz_0, g_zzzzz_0_xyyzzz_1, g_zzzzz_0_xyzzzz_0, g_zzzzz_0_xyzzzz_1, g_zzzzz_0_xzzzzz_0, g_zzzzz_0_xzzzzz_1, g_zzzzz_0_yyyyyy_0, g_zzzzz_0_yyyyyy_1, g_zzzzz_0_yyyyyz_0, g_zzzzz_0_yyyyyz_1, g_zzzzz_0_yyyyzz_0, g_zzzzz_0_yyyyzz_1, g_zzzzz_0_yyyzzz_0, g_zzzzz_0_yyyzzz_1, g_zzzzz_0_yyzzzz_0, g_zzzzz_0_yyzzzz_1, g_zzzzz_0_yzzzzz_0, g_zzzzz_0_yzzzzz_1, g_zzzzz_0_zzzzzz_0, g_zzzzz_0_zzzzzz_1, g_zzzzzz_0_xxxxx_1, g_zzzzzz_0_xxxxxx_1, g_zzzzzz_0_xxxxxy_1, g_zzzzzz_0_xxxxxz_1, g_zzzzzz_0_xxxxy_1, g_zzzzzz_0_xxxxyy_1, g_zzzzzz_0_xxxxyz_1, g_zzzzzz_0_xxxxz_1, g_zzzzzz_0_xxxxzz_1, g_zzzzzz_0_xxxyy_1, g_zzzzzz_0_xxxyyy_1, g_zzzzzz_0_xxxyyz_1, g_zzzzzz_0_xxxyz_1, g_zzzzzz_0_xxxyzz_1, g_zzzzzz_0_xxxzz_1, g_zzzzzz_0_xxxzzz_1, g_zzzzzz_0_xxyyy_1, g_zzzzzz_0_xxyyyy_1, g_zzzzzz_0_xxyyyz_1, g_zzzzzz_0_xxyyz_1, g_zzzzzz_0_xxyyzz_1, g_zzzzzz_0_xxyzz_1, g_zzzzzz_0_xxyzzz_1, g_zzzzzz_0_xxzzz_1, g_zzzzzz_0_xxzzzz_1, g_zzzzzz_0_xyyyy_1, g_zzzzzz_0_xyyyyy_1, g_zzzzzz_0_xyyyyz_1, g_zzzzzz_0_xyyyz_1, g_zzzzzz_0_xyyyzz_1, g_zzzzzz_0_xyyzz_1, g_zzzzzz_0_xyyzzz_1, g_zzzzzz_0_xyzzz_1, g_zzzzzz_0_xyzzzz_1, g_zzzzzz_0_xzzzz_1, g_zzzzzz_0_xzzzzz_1, g_zzzzzz_0_yyyyy_1, g_zzzzzz_0_yyyyyy_1, g_zzzzzz_0_yyyyyz_1, g_zzzzzz_0_yyyyz_1, g_zzzzzz_0_yyyyzz_1, g_zzzzzz_0_yyyzz_1, g_zzzzzz_0_yyyzzz_1, g_zzzzzz_0_yyzzz_1, g_zzzzzz_0_yyzzzz_1, g_zzzzzz_0_yzzzz_1, g_zzzzzz_0_yzzzzz_1, g_zzzzzz_0_zzzzz_1, g_zzzzzz_0_zzzzzz_1, g_zzzzzzz_0_xxxxxx_0, g_zzzzzzz_0_xxxxxy_0, g_zzzzzzz_0_xxxxxz_0, g_zzzzzzz_0_xxxxyy_0, g_zzzzzzz_0_xxxxyz_0, g_zzzzzzz_0_xxxxzz_0, g_zzzzzzz_0_xxxyyy_0, g_zzzzzzz_0_xxxyyz_0, g_zzzzzzz_0_xxxyzz_0, g_zzzzzzz_0_xxxzzz_0, g_zzzzzzz_0_xxyyyy_0, g_zzzzzzz_0_xxyyyz_0, g_zzzzzzz_0_xxyyzz_0, g_zzzzzzz_0_xxyzzz_0, g_zzzzzzz_0_xxzzzz_0, g_zzzzzzz_0_xyyyyy_0, g_zzzzzzz_0_xyyyyz_0, g_zzzzzzz_0_xyyyzz_0, g_zzzzzzz_0_xyyzzz_0, g_zzzzzzz_0_xyzzzz_0, g_zzzzzzz_0_xzzzzz_0, g_zzzzzzz_0_yyyyyy_0, g_zzzzzzz_0_yyyyyz_0, g_zzzzzzz_0_yyyyzz_0, g_zzzzzzz_0_yyyzzz_0, g_zzzzzzz_0_yyzzzz_0, g_zzzzzzz_0_yzzzzz_0, g_zzzzzzz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzzz_0_xxxxxx_0[i] = 6.0 * g_zzzzz_0_xxxxxx_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxx_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxx_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxy_0[i] = 6.0 * g_zzzzz_0_xxxxxy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxy_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxz_0[i] = 6.0 * g_zzzzz_0_xxxxxz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxz_1[i] * fz_be_0 + g_zzzzzz_0_xxxxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxyy_0[i] = 6.0 * g_zzzzz_0_xxxxyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxyy_1[i] * fz_be_0 + g_zzzzzz_0_xxxxyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxyz_0[i] = 6.0 * g_zzzzz_0_xxxxyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxyz_1[i] * fz_be_0 + g_zzzzzz_0_xxxxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxzz_0[i] = 6.0 * g_zzzzz_0_xxxxzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxxxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyyy_0[i] = 6.0 * g_zzzzz_0_xxxyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyyy_1[i] * fz_be_0 + g_zzzzzz_0_xxxyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyyz_0[i] = 6.0 * g_zzzzz_0_xxxyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyyz_1[i] * fz_be_0 + g_zzzzzz_0_xxxyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyzz_0[i] = 6.0 * g_zzzzz_0_xxxyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxxyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxzzz_0[i] = 6.0 * g_zzzzz_0_xxxzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xxxzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyyy_0[i] = 6.0 * g_zzzzz_0_xxyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyyy_1[i] * fz_be_0 + g_zzzzzz_0_xxyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyyz_0[i] = 6.0 * g_zzzzz_0_xxyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyyz_1[i] * fz_be_0 + g_zzzzzz_0_xxyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyzz_0[i] = 6.0 * g_zzzzz_0_xxyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyzzz_0[i] = 6.0 * g_zzzzz_0_xxyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xxyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxzzzz_0[i] = 6.0 * g_zzzzz_0_xxzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_xxzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyyy_0[i] = 6.0 * g_zzzzz_0_xyyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyyy_1[i] * fz_be_0 + g_zzzzzz_0_xyyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyyz_0[i] = 6.0 * g_zzzzz_0_xyyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyyz_1[i] * fz_be_0 + g_zzzzzz_0_xyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyzz_0[i] = 6.0 * g_zzzzz_0_xyyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyzzz_0[i] = 6.0 * g_zzzzz_0_xyyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyzzzz_0[i] = 6.0 * g_zzzzz_0_xyzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_xyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xzzzzz_0[i] = 6.0 * g_zzzzz_0_xzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_0_xzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyyy_0[i] = 6.0 * g_zzzzz_0_yyyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyyy_1[i] * fz_be_0 + g_zzzzzz_0_yyyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyyz_0[i] = 6.0 * g_zzzzz_0_yyyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyyz_1[i] * fz_be_0 + g_zzzzzz_0_yyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyzz_0[i] = 6.0 * g_zzzzz_0_yyyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_yyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyzzz_0[i] = 6.0 * g_zzzzz_0_yyyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_yyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyzzzz_0[i] = 6.0 * g_zzzzz_0_yyzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_yyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yzzzzz_0[i] = 6.0 * g_zzzzz_0_yzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_0_yzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_zzzzzz_0[i] = 6.0 * g_zzzzz_0_zzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_zzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzzz_0_zzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_zzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

