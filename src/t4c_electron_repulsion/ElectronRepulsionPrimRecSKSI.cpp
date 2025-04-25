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

#include "ElectronRepulsionPrimRecSKSI.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sksi(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sksi,
                                  size_t                idx_eri_0_shsi,
                                  size_t                idx_eri_1_shsi,
                                  size_t                idx_eri_1_sish,
                                  size_t                idx_eri_0_sisi,
                                  size_t                idx_eri_1_sisi,
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

    auto g_0_xxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 30);

    auto g_0_xxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 33);

    auto g_0_xxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 37);

    auto g_0_xxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 42);

    auto g_0_xxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 48);

    auto g_0_xxxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_shsi + 56);

    auto g_0_xxxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_shsi + 57);

    auto g_0_xxxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 59);

    auto g_0_xxxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 62);

    auto g_0_xxxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 66);

    auto g_0_xxxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 71);

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

    auto g_0_yyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_shsi + 451);

    auto g_0_yyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_shsi + 454);

    auto g_0_yyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_shsi + 458);

    auto g_0_yyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 463);

    auto g_0_yyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_shsi + 469);

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

    auto g_0_yzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_shsi + 534);

    auto g_0_yzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_shsi + 536);

    auto g_0_yzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_shsi + 537);

    auto g_0_yzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_shsi + 539);

    auto g_0_yzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_shsi + 540);

    auto g_0_yzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_shsi + 541);

    auto g_0_yzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_shsi + 543);

    auto g_0_yzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_shsi + 544);

    auto g_0_yzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_shsi + 545);

    auto g_0_yzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_shsi + 546);

    auto g_0_yzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_shsi + 548);

    auto g_0_yzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_shsi + 549);

    auto g_0_yzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_shsi + 550);

    auto g_0_yzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_shsi + 551);

    auto g_0_yzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_shsi + 552);

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

    auto g_0_xxxxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 30);

    auto g_0_xxxxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 33);

    auto g_0_xxxxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 37);

    auto g_0_xxxxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 42);

    auto g_0_xxxxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 48);

    auto g_0_xxxxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_shsi + 56);

    auto g_0_xxxxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_shsi + 57);

    auto g_0_xxxxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 59);

    auto g_0_xxxxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 62);

    auto g_0_xxxxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 66);

    auto g_0_xxxxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 71);

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

    auto g_0_yyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_shsi + 451);

    auto g_0_yyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_shsi + 454);

    auto g_0_yyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_shsi + 458);

    auto g_0_yyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 463);

    auto g_0_yyyyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_shsi + 469);

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

    auto g_0_yzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_shsi + 534);

    auto g_0_yzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_shsi + 536);

    auto g_0_yzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_shsi + 537);

    auto g_0_yzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_shsi + 539);

    auto g_0_yzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_shsi + 540);

    auto g_0_yzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_shsi + 541);

    auto g_0_yzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_shsi + 543);

    auto g_0_yzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_shsi + 544);

    auto g_0_yzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_shsi + 545);

    auto g_0_yzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_shsi + 546);

    auto g_0_yzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_shsi + 548);

    auto g_0_yzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_shsi + 549);

    auto g_0_yzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_shsi + 550);

    auto g_0_yzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_shsi + 551);

    auto g_0_yzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_shsi + 552);

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

    /// Set up components of auxilary buffer : SISH

    auto g_0_xxxxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish);

    auto g_0_xxxxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 1);

    auto g_0_xxxxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 2);

    auto g_0_xxxxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 3);

    auto g_0_xxxxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 4);

    auto g_0_xxxxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 5);

    auto g_0_xxxxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 6);

    auto g_0_xxxxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 7);

    auto g_0_xxxxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 8);

    auto g_0_xxxxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 9);

    auto g_0_xxxxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 10);

    auto g_0_xxxxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 11);

    auto g_0_xxxxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 12);

    auto g_0_xxxxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 13);

    auto g_0_xxxxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 14);

    auto g_0_xxxxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 15);

    auto g_0_xxxxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 16);

    auto g_0_xxxxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 17);

    auto g_0_xxxxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 18);

    auto g_0_xxxxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 19);

    auto g_0_xxxxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 20);

    auto g_0_xxxxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 44);

    auto g_0_xxxxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 46);

    auto g_0_xxxxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 47);

    auto g_0_xxxxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 49);

    auto g_0_xxxxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 50);

    auto g_0_xxxxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 51);

    auto g_0_xxxxxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 53);

    auto g_0_xxxxxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 54);

    auto g_0_xxxxxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 55);

    auto g_0_xxxxxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 56);

    auto g_0_xxxxxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 58);

    auto g_0_xxxxxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 59);

    auto g_0_xxxxxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 60);

    auto g_0_xxxxxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 61);

    auto g_0_xxxxxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 62);

    auto g_0_xxxxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 63);

    auto g_0_xxxxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 64);

    auto g_0_xxxxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 65);

    auto g_0_xxxxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 66);

    auto g_0_xxxxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 67);

    auto g_0_xxxxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 68);

    auto g_0_xxxxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 69);

    auto g_0_xxxxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 70);

    auto g_0_xxxxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 71);

    auto g_0_xxxxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 72);

    auto g_0_xxxxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 73);

    auto g_0_xxxxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 74);

    auto g_0_xxxxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 75);

    auto g_0_xxxxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 76);

    auto g_0_xxxxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 77);

    auto g_0_xxxxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 78);

    auto g_0_xxxxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 79);

    auto g_0_xxxxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 80);

    auto g_0_xxxxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 81);

    auto g_0_xxxxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 82);

    auto g_0_xxxxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 83);

    auto g_0_xxxxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 105);

    auto g_0_xxxxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 106);

    auto g_0_xxxxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 107);

    auto g_0_xxxxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 108);

    auto g_0_xxxxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 109);

    auto g_0_xxxxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 110);

    auto g_0_xxxxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 111);

    auto g_0_xxxxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 112);

    auto g_0_xxxxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 113);

    auto g_0_xxxxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 114);

    auto g_0_xxxxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 115);

    auto g_0_xxxxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 116);

    auto g_0_xxxxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 117);

    auto g_0_xxxxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 118);

    auto g_0_xxxxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 119);

    auto g_0_xxxxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 120);

    auto g_0_xxxxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 121);

    auto g_0_xxxxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 122);

    auto g_0_xxxxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 123);

    auto g_0_xxxxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 124);

    auto g_0_xxxxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 125);

    auto g_0_xxxyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 126);

    auto g_0_xxxyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 127);

    auto g_0_xxxyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 128);

    auto g_0_xxxyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 129);

    auto g_0_xxxyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 130);

    auto g_0_xxxyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 131);

    auto g_0_xxxyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 132);

    auto g_0_xxxyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 133);

    auto g_0_xxxyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 134);

    auto g_0_xxxyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 135);

    auto g_0_xxxyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 136);

    auto g_0_xxxyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 137);

    auto g_0_xxxyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 138);

    auto g_0_xxxyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 139);

    auto g_0_xxxyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 140);

    auto g_0_xxxyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 141);

    auto g_0_xxxyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 142);

    auto g_0_xxxyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 143);

    auto g_0_xxxyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 144);

    auto g_0_xxxyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 145);

    auto g_0_xxxyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 146);

    auto g_0_xxxzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 189);

    auto g_0_xxxzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 190);

    auto g_0_xxxzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 191);

    auto g_0_xxxzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 192);

    auto g_0_xxxzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 193);

    auto g_0_xxxzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 194);

    auto g_0_xxxzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 195);

    auto g_0_xxxzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 196);

    auto g_0_xxxzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 197);

    auto g_0_xxxzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 198);

    auto g_0_xxxzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 199);

    auto g_0_xxxzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 200);

    auto g_0_xxxzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 201);

    auto g_0_xxxzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 202);

    auto g_0_xxxzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 203);

    auto g_0_xxxzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 204);

    auto g_0_xxxzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 205);

    auto g_0_xxxzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 206);

    auto g_0_xxxzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 207);

    auto g_0_xxxzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 208);

    auto g_0_xxxzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 209);

    auto g_0_xxyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 210);

    auto g_0_xxyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 211);

    auto g_0_xxyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 212);

    auto g_0_xxyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 213);

    auto g_0_xxyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 214);

    auto g_0_xxyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 215);

    auto g_0_xxyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 216);

    auto g_0_xxyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 217);

    auto g_0_xxyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 218);

    auto g_0_xxyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 219);

    auto g_0_xxyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 220);

    auto g_0_xxyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 221);

    auto g_0_xxyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 222);

    auto g_0_xxyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 223);

    auto g_0_xxyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 224);

    auto g_0_xxyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 225);

    auto g_0_xxyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 226);

    auto g_0_xxyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 227);

    auto g_0_xxyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 228);

    auto g_0_xxyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 229);

    auto g_0_xxyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 230);

    auto g_0_xxyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 256);

    auto g_0_xxyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 259);

    auto g_0_xxyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 260);

    auto g_0_xxyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 263);

    auto g_0_xxyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 264);

    auto g_0_xxyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 265);

    auto g_0_xxyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 268);

    auto g_0_xxyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 269);

    auto g_0_xxyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 270);

    auto g_0_xxyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 271);

    auto g_0_xxzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 294);

    auto g_0_xxzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 295);

    auto g_0_xxzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 296);

    auto g_0_xxzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 297);

    auto g_0_xxzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 298);

    auto g_0_xxzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 299);

    auto g_0_xxzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 300);

    auto g_0_xxzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 301);

    auto g_0_xxzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 302);

    auto g_0_xxzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 303);

    auto g_0_xxzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 304);

    auto g_0_xxzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 305);

    auto g_0_xxzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 306);

    auto g_0_xxzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 307);

    auto g_0_xxzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 308);

    auto g_0_xxzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 309);

    auto g_0_xxzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 310);

    auto g_0_xxzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 311);

    auto g_0_xxzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 312);

    auto g_0_xxzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 313);

    auto g_0_xxzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 314);

    auto g_0_xyyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 316);

    auto g_0_xyyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 318);

    auto g_0_xyyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 319);

    auto g_0_xyyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 321);

    auto g_0_xyyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 322);

    auto g_0_xyyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 323);

    auto g_0_xyyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 325);

    auto g_0_xyyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 326);

    auto g_0_xyyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 327);

    auto g_0_xyyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 328);

    auto g_0_xyyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 330);

    auto g_0_xyyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 331);

    auto g_0_xyyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 332);

    auto g_0_xyyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 333);

    auto g_0_xyyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 334);

    auto g_0_xyyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 361);

    auto g_0_xyyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 364);

    auto g_0_xyyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 365);

    auto g_0_xyyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 368);

    auto g_0_xyyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 369);

    auto g_0_xyyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 370);

    auto g_0_xyyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 373);

    auto g_0_xyyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 374);

    auto g_0_xyyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 375);

    auto g_0_xyyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 376);

    auto g_0_xyyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 382);

    auto g_0_xyyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 385);

    auto g_0_xyyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 386);

    auto g_0_xyyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 389);

    auto g_0_xyyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 390);

    auto g_0_xyyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 391);

    auto g_0_xyyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 394);

    auto g_0_xyyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 395);

    auto g_0_xyyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 396);

    auto g_0_xyyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 397);

    auto g_0_xzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 422);

    auto g_0_xzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 424);

    auto g_0_xzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 425);

    auto g_0_xzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 427);

    auto g_0_xzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 428);

    auto g_0_xzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 429);

    auto g_0_xzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 431);

    auto g_0_xzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 432);

    auto g_0_xzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 433);

    auto g_0_xzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 434);

    auto g_0_xzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 436);

    auto g_0_xzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 437);

    auto g_0_xzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 438);

    auto g_0_xzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 439);

    auto g_0_xzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 440);

    auto g_0_yyyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 441);

    auto g_0_yyyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 442);

    auto g_0_yyyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 443);

    auto g_0_yyyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 444);

    auto g_0_yyyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 445);

    auto g_0_yyyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 446);

    auto g_0_yyyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 447);

    auto g_0_yyyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 448);

    auto g_0_yyyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 449);

    auto g_0_yyyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 450);

    auto g_0_yyyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 451);

    auto g_0_yyyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 452);

    auto g_0_yyyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 453);

    auto g_0_yyyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 454);

    auto g_0_yyyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 455);

    auto g_0_yyyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 456);

    auto g_0_yyyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 457);

    auto g_0_yyyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 458);

    auto g_0_yyyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 459);

    auto g_0_yyyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 460);

    auto g_0_yyyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 461);

    auto g_0_yyyyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 464);

    auto g_0_yyyyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 466);

    auto g_0_yyyyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 467);

    auto g_0_yyyyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 469);

    auto g_0_yyyyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 470);

    auto g_0_yyyyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 471);

    auto g_0_yyyyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 473);

    auto g_0_yyyyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 474);

    auto g_0_yyyyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 475);

    auto g_0_yyyyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 476);

    auto g_0_yyyyyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 478);

    auto g_0_yyyyyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 479);

    auto g_0_yyyyyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 480);

    auto g_0_yyyyyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 481);

    auto g_0_yyyyyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 482);

    auto g_0_yyyyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 483);

    auto g_0_yyyyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 484);

    auto g_0_yyyyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 485);

    auto g_0_yyyyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 486);

    auto g_0_yyyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 487);

    auto g_0_yyyyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 488);

    auto g_0_yyyyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 489);

    auto g_0_yyyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 490);

    auto g_0_yyyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 491);

    auto g_0_yyyyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 492);

    auto g_0_yyyyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 493);

    auto g_0_yyyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 494);

    auto g_0_yyyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 495);

    auto g_0_yyyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 496);

    auto g_0_yyyyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 497);

    auto g_0_yyyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 498);

    auto g_0_yyyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 499);

    auto g_0_yyyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 500);

    auto g_0_yyyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 501);

    auto g_0_yyyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 502);

    auto g_0_yyyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 503);

    auto g_0_yyyzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 504);

    auto g_0_yyyzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 505);

    auto g_0_yyyzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 506);

    auto g_0_yyyzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 507);

    auto g_0_yyyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 508);

    auto g_0_yyyzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 509);

    auto g_0_yyyzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 510);

    auto g_0_yyyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 511);

    auto g_0_yyyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 512);

    auto g_0_yyyzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 513);

    auto g_0_yyyzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 514);

    auto g_0_yyyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 515);

    auto g_0_yyyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 516);

    auto g_0_yyyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 517);

    auto g_0_yyyzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 518);

    auto g_0_yyyzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 519);

    auto g_0_yyyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 520);

    auto g_0_yyyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 521);

    auto g_0_yyyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 522);

    auto g_0_yyyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 523);

    auto g_0_yyyzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 524);

    auto g_0_yyzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 525);

    auto g_0_yyzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 526);

    auto g_0_yyzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 527);

    auto g_0_yyzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 528);

    auto g_0_yyzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 529);

    auto g_0_yyzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 530);

    auto g_0_yyzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 531);

    auto g_0_yyzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 532);

    auto g_0_yyzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 533);

    auto g_0_yyzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 534);

    auto g_0_yyzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 535);

    auto g_0_yyzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 536);

    auto g_0_yyzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 537);

    auto g_0_yyzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 538);

    auto g_0_yyzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 539);

    auto g_0_yyzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 540);

    auto g_0_yyzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 541);

    auto g_0_yyzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 542);

    auto g_0_yyzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 543);

    auto g_0_yyzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 544);

    auto g_0_yyzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 545);

    auto g_0_yzzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 547);

    auto g_0_yzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 548);

    auto g_0_yzzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 549);

    auto g_0_yzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 550);

    auto g_0_yzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 551);

    auto g_0_yzzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 552);

    auto g_0_yzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 553);

    auto g_0_yzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 554);

    auto g_0_yzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 555);

    auto g_0_yzzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 556);

    auto g_0_yzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 557);

    auto g_0_yzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 558);

    auto g_0_yzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 559);

    auto g_0_yzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 560);

    auto g_0_yzzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 561);

    auto g_0_yzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 562);

    auto g_0_yzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 563);

    auto g_0_yzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 564);

    auto g_0_yzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 565);

    auto g_0_yzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 566);

    auto g_0_zzzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sish + 567);

    auto g_0_zzzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sish + 568);

    auto g_0_zzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sish + 569);

    auto g_0_zzzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sish + 570);

    auto g_0_zzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sish + 571);

    auto g_0_zzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sish + 572);

    auto g_0_zzzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sish + 573);

    auto g_0_zzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sish + 574);

    auto g_0_zzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sish + 575);

    auto g_0_zzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sish + 576);

    auto g_0_zzzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sish + 577);

    auto g_0_zzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sish + 578);

    auto g_0_zzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sish + 579);

    auto g_0_zzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sish + 580);

    auto g_0_zzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sish + 581);

    auto g_0_zzzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sish + 582);

    auto g_0_zzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sish + 583);

    auto g_0_zzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sish + 584);

    auto g_0_zzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sish + 585);

    auto g_0_zzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sish + 586);

    auto g_0_zzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sish + 587);

    /// Set up components of auxilary buffer : SISI

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

    auto g_0_xxxxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 28);

    auto g_0_xxxxxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 29);

    auto g_0_xxxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 30);

    auto g_0_xxxxxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 31);

    auto g_0_xxxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 33);

    auto g_0_xxxxxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 34);

    auto g_0_xxxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 37);

    auto g_0_xxxxxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 38);

    auto g_0_xxxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 42);

    auto g_0_xxxxxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 43);

    auto g_0_xxxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 48);

    auto g_0_xxxxxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 49);

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

    auto g_0_xxxxxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 78);

    auto g_0_xxxxxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 79);

    auto g_0_xxxxxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 80);

    auto g_0_xxxxxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 81);

    auto g_0_xxxxxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 82);

    auto g_0_xxxxxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 83);

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

    auto g_0_xxxyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 197);

    auto g_0_xxxyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 199);

    auto g_0_xxxyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 202);

    auto g_0_xxxyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 206);

    auto g_0_xxxyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 211);

    auto g_0_xxxyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 224);

    auto g_0_xxxyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 226);

    auto g_0_xxxyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 229);

    auto g_0_xxxyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 233);

    auto g_0_xxxyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 238);

    auto g_0_xxxyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 244);

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

    auto g_0_xxyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 309);

    auto g_0_xxyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 311);

    auto g_0_xxyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 314);

    auto g_0_xxyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 318);

    auto g_0_xxyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 323);

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

    auto g_0_xxyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 364);

    auto g_0_xxyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 366);

    auto g_0_xxyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 369);

    auto g_0_xxyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 373);

    auto g_0_xxyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 378);

    auto g_0_xxyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 384);

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

    auto g_0_xyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 420);

    auto g_0_xyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 421);

    auto g_0_xyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 423);

    auto g_0_xyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 424);

    auto g_0_xyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 426);

    auto g_0_xyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 427);

    auto g_0_xyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 428);

    auto g_0_xyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 430);

    auto g_0_xyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 431);

    auto g_0_xyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 432);

    auto g_0_xyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 433);

    auto g_0_xyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 435);

    auto g_0_xyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 436);

    auto g_0_xyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 437);

    auto g_0_xyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 438);

    auto g_0_xyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 439);

    auto g_0_xyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 441);

    auto g_0_xyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 442);

    auto g_0_xyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 443);

    auto g_0_xyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 444);

    auto g_0_xyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 445);

    auto g_0_xyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 446);

    auto g_0_xyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 447);

    auto g_0_xyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 480);

    auto g_0_xyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 483);

    auto g_0_xyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 484);

    auto g_0_xyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 487);

    auto g_0_xyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 488);

    auto g_0_xyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 489);

    auto g_0_xyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 492);

    auto g_0_xyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 493);

    auto g_0_xyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 494);

    auto g_0_xyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 495);

    auto g_0_xyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 497);

    auto g_0_xyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 498);

    auto g_0_xyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 499);

    auto g_0_xyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 500);

    auto g_0_xyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 501);

    auto g_0_xyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 502);

    auto g_0_xyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 503);

    auto g_0_xyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 508);

    auto g_0_xyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 511);

    auto g_0_xyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 512);

    auto g_0_xyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 515);

    auto g_0_xyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 516);

    auto g_0_xyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 517);

    auto g_0_xyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 520);

    auto g_0_xyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 521);

    auto g_0_xyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 522);

    auto g_0_xyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 523);

    auto g_0_xyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 525);

    auto g_0_xyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 526);

    auto g_0_xyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 527);

    auto g_0_xyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 528);

    auto g_0_xyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 529);

    auto g_0_xyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 530);

    auto g_0_xyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 531);

    auto g_0_xzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 560);

    auto g_0_xzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 562);

    auto g_0_xzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 564);

    auto g_0_xzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 565);

    auto g_0_xzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 567);

    auto g_0_xzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 568);

    auto g_0_xzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 569);

    auto g_0_xzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 571);

    auto g_0_xzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 572);

    auto g_0_xzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 573);

    auto g_0_xzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 574);

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

    /// Set up components of auxilary buffer : SISI

    auto g_0_xxxxxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi);

    auto g_0_xxxxxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 1);

    auto g_0_xxxxxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 2);

    auto g_0_xxxxxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 3);

    auto g_0_xxxxxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 4);

    auto g_0_xxxxxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 5);

    auto g_0_xxxxxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 6);

    auto g_0_xxxxxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 7);

    auto g_0_xxxxxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 8);

    auto g_0_xxxxxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 9);

    auto g_0_xxxxxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 10);

    auto g_0_xxxxxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 11);

    auto g_0_xxxxxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 12);

    auto g_0_xxxxxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 13);

    auto g_0_xxxxxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 14);

    auto g_0_xxxxxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 15);

    auto g_0_xxxxxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 16);

    auto g_0_xxxxxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 17);

    auto g_0_xxxxxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 18);

    auto g_0_xxxxxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 19);

    auto g_0_xxxxxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 20);

    auto g_0_xxxxxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 21);

    auto g_0_xxxxxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 22);

    auto g_0_xxxxxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 23);

    auto g_0_xxxxxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 24);

    auto g_0_xxxxxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 25);

    auto g_0_xxxxxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 26);

    auto g_0_xxxxxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 27);

    auto g_0_xxxxxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 28);

    auto g_0_xxxxxy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 29);

    auto g_0_xxxxxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 30);

    auto g_0_xxxxxy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 31);

    auto g_0_xxxxxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 33);

    auto g_0_xxxxxy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 34);

    auto g_0_xxxxxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 37);

    auto g_0_xxxxxy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 38);

    auto g_0_xxxxxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 42);

    auto g_0_xxxxxy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 43);

    auto g_0_xxxxxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 48);

    auto g_0_xxxxxy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 49);

    auto g_0_xxxxxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 56);

    auto g_0_xxxxxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 57);

    auto g_0_xxxxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 58);

    auto g_0_xxxxxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 59);

    auto g_0_xxxxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 60);

    auto g_0_xxxxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 61);

    auto g_0_xxxxxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 62);

    auto g_0_xxxxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 63);

    auto g_0_xxxxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 64);

    auto g_0_xxxxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 65);

    auto g_0_xxxxxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 66);

    auto g_0_xxxxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 67);

    auto g_0_xxxxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 68);

    auto g_0_xxxxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 69);

    auto g_0_xxxxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 70);

    auto g_0_xxxxxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 71);

    auto g_0_xxxxxz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 72);

    auto g_0_xxxxxz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 73);

    auto g_0_xxxxxz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 74);

    auto g_0_xxxxxz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 75);

    auto g_0_xxxxxz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 76);

    auto g_0_xxxxxz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 78);

    auto g_0_xxxxxz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 79);

    auto g_0_xxxxxz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 80);

    auto g_0_xxxxxz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 81);

    auto g_0_xxxxxz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 82);

    auto g_0_xxxxxz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 83);

    auto g_0_xxxxyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 84);

    auto g_0_xxxxyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 85);

    auto g_0_xxxxyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 86);

    auto g_0_xxxxyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 87);

    auto g_0_xxxxyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 88);

    auto g_0_xxxxyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 89);

    auto g_0_xxxxyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 90);

    auto g_0_xxxxyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 91);

    auto g_0_xxxxyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 92);

    auto g_0_xxxxyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 93);

    auto g_0_xxxxyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 94);

    auto g_0_xxxxyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 95);

    auto g_0_xxxxyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 96);

    auto g_0_xxxxyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 97);

    auto g_0_xxxxyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 98);

    auto g_0_xxxxyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 99);

    auto g_0_xxxxyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 100);

    auto g_0_xxxxyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 101);

    auto g_0_xxxxyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 102);

    auto g_0_xxxxyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 103);

    auto g_0_xxxxyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 104);

    auto g_0_xxxxyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 105);

    auto g_0_xxxxyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 106);

    auto g_0_xxxxyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 107);

    auto g_0_xxxxyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 108);

    auto g_0_xxxxyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 109);

    auto g_0_xxxxyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 110);

    auto g_0_xxxxyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 111);

    auto g_0_xxxxzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 140);

    auto g_0_xxxxzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 141);

    auto g_0_xxxxzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 142);

    auto g_0_xxxxzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 143);

    auto g_0_xxxxzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 144);

    auto g_0_xxxxzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 145);

    auto g_0_xxxxzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 146);

    auto g_0_xxxxzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 147);

    auto g_0_xxxxzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 148);

    auto g_0_xxxxzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 149);

    auto g_0_xxxxzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 150);

    auto g_0_xxxxzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 151);

    auto g_0_xxxxzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 152);

    auto g_0_xxxxzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 153);

    auto g_0_xxxxzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 154);

    auto g_0_xxxxzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 155);

    auto g_0_xxxxzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 156);

    auto g_0_xxxxzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 157);

    auto g_0_xxxxzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 158);

    auto g_0_xxxxzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 159);

    auto g_0_xxxxzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 160);

    auto g_0_xxxxzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 161);

    auto g_0_xxxxzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 162);

    auto g_0_xxxxzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 163);

    auto g_0_xxxxzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 164);

    auto g_0_xxxxzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 165);

    auto g_0_xxxxzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 166);

    auto g_0_xxxxzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 167);

    auto g_0_xxxyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 168);

    auto g_0_xxxyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 169);

    auto g_0_xxxyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 170);

    auto g_0_xxxyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 171);

    auto g_0_xxxyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 172);

    auto g_0_xxxyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 173);

    auto g_0_xxxyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 174);

    auto g_0_xxxyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 175);

    auto g_0_xxxyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 176);

    auto g_0_xxxyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 177);

    auto g_0_xxxyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 178);

    auto g_0_xxxyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 179);

    auto g_0_xxxyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 180);

    auto g_0_xxxyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 181);

    auto g_0_xxxyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 182);

    auto g_0_xxxyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 183);

    auto g_0_xxxyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 184);

    auto g_0_xxxyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 185);

    auto g_0_xxxyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 186);

    auto g_0_xxxyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 187);

    auto g_0_xxxyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 188);

    auto g_0_xxxyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 189);

    auto g_0_xxxyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 190);

    auto g_0_xxxyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 191);

    auto g_0_xxxyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 192);

    auto g_0_xxxyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 193);

    auto g_0_xxxyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 194);

    auto g_0_xxxyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 195);

    auto g_0_xxxyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 197);

    auto g_0_xxxyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 199);

    auto g_0_xxxyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 202);

    auto g_0_xxxyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 206);

    auto g_0_xxxyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 211);

    auto g_0_xxxyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 224);

    auto g_0_xxxyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 226);

    auto g_0_xxxyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 229);

    auto g_0_xxxyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 233);

    auto g_0_xxxyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 238);

    auto g_0_xxxyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 244);

    auto g_0_xxxzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 252);

    auto g_0_xxxzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 253);

    auto g_0_xxxzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 254);

    auto g_0_xxxzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 255);

    auto g_0_xxxzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 256);

    auto g_0_xxxzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 257);

    auto g_0_xxxzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 258);

    auto g_0_xxxzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 259);

    auto g_0_xxxzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 260);

    auto g_0_xxxzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 261);

    auto g_0_xxxzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 262);

    auto g_0_xxxzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 263);

    auto g_0_xxxzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 264);

    auto g_0_xxxzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 265);

    auto g_0_xxxzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 266);

    auto g_0_xxxzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 267);

    auto g_0_xxxzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 268);

    auto g_0_xxxzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 269);

    auto g_0_xxxzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 270);

    auto g_0_xxxzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 271);

    auto g_0_xxxzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 272);

    auto g_0_xxxzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 273);

    auto g_0_xxxzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 274);

    auto g_0_xxxzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 275);

    auto g_0_xxxzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 276);

    auto g_0_xxxzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 277);

    auto g_0_xxxzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 278);

    auto g_0_xxxzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 279);

    auto g_0_xxyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 280);

    auto g_0_xxyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 281);

    auto g_0_xxyyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 282);

    auto g_0_xxyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 283);

    auto g_0_xxyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 284);

    auto g_0_xxyyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 285);

    auto g_0_xxyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 286);

    auto g_0_xxyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 287);

    auto g_0_xxyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 288);

    auto g_0_xxyyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 289);

    auto g_0_xxyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 290);

    auto g_0_xxyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 291);

    auto g_0_xxyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 292);

    auto g_0_xxyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 293);

    auto g_0_xxyyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 294);

    auto g_0_xxyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 295);

    auto g_0_xxyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 296);

    auto g_0_xxyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 297);

    auto g_0_xxyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 298);

    auto g_0_xxyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 299);

    auto g_0_xxyyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 300);

    auto g_0_xxyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 301);

    auto g_0_xxyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 302);

    auto g_0_xxyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 303);

    auto g_0_xxyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 304);

    auto g_0_xxyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 305);

    auto g_0_xxyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 306);

    auto g_0_xxyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 307);

    auto g_0_xxyyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 309);

    auto g_0_xxyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 311);

    auto g_0_xxyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 314);

    auto g_0_xxyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 318);

    auto g_0_xxyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 323);

    auto g_0_xxyyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 336);

    auto g_0_xxyyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 337);

    auto g_0_xxyyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 338);

    auto g_0_xxyyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 339);

    auto g_0_xxyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 340);

    auto g_0_xxyyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 341);

    auto g_0_xxyyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 342);

    auto g_0_xxyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 343);

    auto g_0_xxyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 344);

    auto g_0_xxyyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 345);

    auto g_0_xxyyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 346);

    auto g_0_xxyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 347);

    auto g_0_xxyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 348);

    auto g_0_xxyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 349);

    auto g_0_xxyyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 350);

    auto g_0_xxyyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 351);

    auto g_0_xxyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 352);

    auto g_0_xxyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 353);

    auto g_0_xxyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 354);

    auto g_0_xxyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 355);

    auto g_0_xxyyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 356);

    auto g_0_xxyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 357);

    auto g_0_xxyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 358);

    auto g_0_xxyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 359);

    auto g_0_xxyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 360);

    auto g_0_xxyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 361);

    auto g_0_xxyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 362);

    auto g_0_xxyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 363);

    auto g_0_xxyzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 364);

    auto g_0_xxyzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 366);

    auto g_0_xxyzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 369);

    auto g_0_xxyzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 373);

    auto g_0_xxyzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 378);

    auto g_0_xxyzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 384);

    auto g_0_xxzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 392);

    auto g_0_xxzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 393);

    auto g_0_xxzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 394);

    auto g_0_xxzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 395);

    auto g_0_xxzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 396);

    auto g_0_xxzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 397);

    auto g_0_xxzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 398);

    auto g_0_xxzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 399);

    auto g_0_xxzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 400);

    auto g_0_xxzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 401);

    auto g_0_xxzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 402);

    auto g_0_xxzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 403);

    auto g_0_xxzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 404);

    auto g_0_xxzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 405);

    auto g_0_xxzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 406);

    auto g_0_xxzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 407);

    auto g_0_xxzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 408);

    auto g_0_xxzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 409);

    auto g_0_xxzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 410);

    auto g_0_xxzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 411);

    auto g_0_xxzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 412);

    auto g_0_xxzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 413);

    auto g_0_xxzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 414);

    auto g_0_xxzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 415);

    auto g_0_xxzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 416);

    auto g_0_xxzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 417);

    auto g_0_xxzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 418);

    auto g_0_xxzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 419);

    auto g_0_xyyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 420);

    auto g_0_xyyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 421);

    auto g_0_xyyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 423);

    auto g_0_xyyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 424);

    auto g_0_xyyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 426);

    auto g_0_xyyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 427);

    auto g_0_xyyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 428);

    auto g_0_xyyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 430);

    auto g_0_xyyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 431);

    auto g_0_xyyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 432);

    auto g_0_xyyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 433);

    auto g_0_xyyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 435);

    auto g_0_xyyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 436);

    auto g_0_xyyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 437);

    auto g_0_xyyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 438);

    auto g_0_xyyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 439);

    auto g_0_xyyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 441);

    auto g_0_xyyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 442);

    auto g_0_xyyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 443);

    auto g_0_xyyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 444);

    auto g_0_xyyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 445);

    auto g_0_xyyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 446);

    auto g_0_xyyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 447);

    auto g_0_xyyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 480);

    auto g_0_xyyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 483);

    auto g_0_xyyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 484);

    auto g_0_xyyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 487);

    auto g_0_xyyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 488);

    auto g_0_xyyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 489);

    auto g_0_xyyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 492);

    auto g_0_xyyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 493);

    auto g_0_xyyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 494);

    auto g_0_xyyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 495);

    auto g_0_xyyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 497);

    auto g_0_xyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 498);

    auto g_0_xyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 499);

    auto g_0_xyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 500);

    auto g_0_xyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 501);

    auto g_0_xyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 502);

    auto g_0_xyyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 503);

    auto g_0_xyyzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 508);

    auto g_0_xyyzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 511);

    auto g_0_xyyzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 512);

    auto g_0_xyyzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 515);

    auto g_0_xyyzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 516);

    auto g_0_xyyzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 517);

    auto g_0_xyyzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 520);

    auto g_0_xyyzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 521);

    auto g_0_xyyzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 522);

    auto g_0_xyyzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 523);

    auto g_0_xyyzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 525);

    auto g_0_xyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 526);

    auto g_0_xyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 527);

    auto g_0_xyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 528);

    auto g_0_xyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 529);

    auto g_0_xyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 530);

    auto g_0_xyyzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 531);

    auto g_0_xzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 560);

    auto g_0_xzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 562);

    auto g_0_xzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 564);

    auto g_0_xzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 565);

    auto g_0_xzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 567);

    auto g_0_xzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 568);

    auto g_0_xzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 569);

    auto g_0_xzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 571);

    auto g_0_xzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 572);

    auto g_0_xzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 573);

    auto g_0_xzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 574);

    auto g_0_xzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 576);

    auto g_0_xzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 577);

    auto g_0_xzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 578);

    auto g_0_xzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 579);

    auto g_0_xzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 580);

    auto g_0_xzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 581);

    auto g_0_xzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 582);

    auto g_0_xzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 583);

    auto g_0_xzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 584);

    auto g_0_xzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 585);

    auto g_0_xzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 586);

    auto g_0_xzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 587);

    auto g_0_yyyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 588);

    auto g_0_yyyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 589);

    auto g_0_yyyyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 590);

    auto g_0_yyyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 591);

    auto g_0_yyyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 592);

    auto g_0_yyyyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 593);

    auto g_0_yyyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 594);

    auto g_0_yyyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 595);

    auto g_0_yyyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 596);

    auto g_0_yyyyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 597);

    auto g_0_yyyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 598);

    auto g_0_yyyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 599);

    auto g_0_yyyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 600);

    auto g_0_yyyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 601);

    auto g_0_yyyyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 602);

    auto g_0_yyyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 603);

    auto g_0_yyyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 604);

    auto g_0_yyyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 605);

    auto g_0_yyyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 606);

    auto g_0_yyyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 607);

    auto g_0_yyyyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 608);

    auto g_0_yyyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 609);

    auto g_0_yyyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 610);

    auto g_0_yyyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 611);

    auto g_0_yyyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 612);

    auto g_0_yyyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 613);

    auto g_0_yyyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 614);

    auto g_0_yyyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 615);

    auto g_0_yyyyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 617);

    auto g_0_yyyyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 618);

    auto g_0_yyyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 619);

    auto g_0_yyyyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 620);

    auto g_0_yyyyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 621);

    auto g_0_yyyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 622);

    auto g_0_yyyyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 623);

    auto g_0_yyyyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 624);

    auto g_0_yyyyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 625);

    auto g_0_yyyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 626);

    auto g_0_yyyyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 627);

    auto g_0_yyyyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 628);

    auto g_0_yyyyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 629);

    auto g_0_yyyyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 630);

    auto g_0_yyyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 631);

    auto g_0_yyyyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 632);

    auto g_0_yyyyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 633);

    auto g_0_yyyyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 634);

    auto g_0_yyyyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 635);

    auto g_0_yyyyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 636);

    auto g_0_yyyyyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 637);

    auto g_0_yyyyyz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 638);

    auto g_0_yyyyyz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 639);

    auto g_0_yyyyyz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 640);

    auto g_0_yyyyyz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 641);

    auto g_0_yyyyyz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 642);

    auto g_0_yyyyyz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 643);

    auto g_0_yyyyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 644);

    auto g_0_yyyyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 645);

    auto g_0_yyyyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 646);

    auto g_0_yyyyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 647);

    auto g_0_yyyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 648);

    auto g_0_yyyyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 649);

    auto g_0_yyyyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 650);

    auto g_0_yyyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 651);

    auto g_0_yyyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 652);

    auto g_0_yyyyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 653);

    auto g_0_yyyyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 654);

    auto g_0_yyyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 655);

    auto g_0_yyyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 656);

    auto g_0_yyyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 657);

    auto g_0_yyyyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 658);

    auto g_0_yyyyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 659);

    auto g_0_yyyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 660);

    auto g_0_yyyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 661);

    auto g_0_yyyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 662);

    auto g_0_yyyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 663);

    auto g_0_yyyyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 664);

    auto g_0_yyyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 665);

    auto g_0_yyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 666);

    auto g_0_yyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 667);

    auto g_0_yyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 668);

    auto g_0_yyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 669);

    auto g_0_yyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 670);

    auto g_0_yyyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 671);

    auto g_0_yyyzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 672);

    auto g_0_yyyzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 673);

    auto g_0_yyyzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 674);

    auto g_0_yyyzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 675);

    auto g_0_yyyzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 676);

    auto g_0_yyyzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 677);

    auto g_0_yyyzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 678);

    auto g_0_yyyzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 679);

    auto g_0_yyyzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 680);

    auto g_0_yyyzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 681);

    auto g_0_yyyzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 682);

    auto g_0_yyyzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 683);

    auto g_0_yyyzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 684);

    auto g_0_yyyzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 685);

    auto g_0_yyyzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 686);

    auto g_0_yyyzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 687);

    auto g_0_yyyzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 688);

    auto g_0_yyyzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 689);

    auto g_0_yyyzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 690);

    auto g_0_yyyzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 691);

    auto g_0_yyyzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 692);

    auto g_0_yyyzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 693);

    auto g_0_yyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 694);

    auto g_0_yyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 695);

    auto g_0_yyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 696);

    auto g_0_yyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 697);

    auto g_0_yyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 698);

    auto g_0_yyyzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 699);

    auto g_0_yyzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 700);

    auto g_0_yyzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 701);

    auto g_0_yyzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 702);

    auto g_0_yyzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 703);

    auto g_0_yyzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 704);

    auto g_0_yyzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 705);

    auto g_0_yyzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 706);

    auto g_0_yyzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 707);

    auto g_0_yyzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 708);

    auto g_0_yyzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 709);

    auto g_0_yyzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 710);

    auto g_0_yyzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 711);

    auto g_0_yyzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 712);

    auto g_0_yyzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 713);

    auto g_0_yyzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 714);

    auto g_0_yyzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 715);

    auto g_0_yyzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 716);

    auto g_0_yyzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 717);

    auto g_0_yyzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 718);

    auto g_0_yyzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 719);

    auto g_0_yyzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 720);

    auto g_0_yyzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 721);

    auto g_0_yyzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 722);

    auto g_0_yyzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 723);

    auto g_0_yyzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 724);

    auto g_0_yyzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 725);

    auto g_0_yyzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 726);

    auto g_0_yyzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 727);

    auto g_0_yzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 728);

    auto g_0_yzzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 729);

    auto g_0_yzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 730);

    auto g_0_yzzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 731);

    auto g_0_yzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 732);

    auto g_0_yzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 733);

    auto g_0_yzzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 734);

    auto g_0_yzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 735);

    auto g_0_yzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 736);

    auto g_0_yzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 737);

    auto g_0_yzzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 738);

    auto g_0_yzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 739);

    auto g_0_yzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 740);

    auto g_0_yzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 741);

    auto g_0_yzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 742);

    auto g_0_yzzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 743);

    auto g_0_yzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 744);

    auto g_0_yzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 745);

    auto g_0_yzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 746);

    auto g_0_yzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 747);

    auto g_0_yzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 748);

    auto g_0_yzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 749);

    auto g_0_yzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 750);

    auto g_0_yzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 751);

    auto g_0_yzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 752);

    auto g_0_yzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 753);

    auto g_0_yzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 754);

    auto g_0_yzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 755);

    auto g_0_zzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 756);

    auto g_0_zzzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 757);

    auto g_0_zzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 758);

    auto g_0_zzzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 759);

    auto g_0_zzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 760);

    auto g_0_zzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 761);

    auto g_0_zzzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 762);

    auto g_0_zzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 763);

    auto g_0_zzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 764);

    auto g_0_zzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 765);

    auto g_0_zzzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 766);

    auto g_0_zzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 767);

    auto g_0_zzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 768);

    auto g_0_zzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 769);

    auto g_0_zzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 770);

    auto g_0_zzzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 771);

    auto g_0_zzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 772);

    auto g_0_zzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 773);

    auto g_0_zzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 774);

    auto g_0_zzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 775);

    auto g_0_zzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 776);

    auto g_0_zzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 777);

    auto g_0_zzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 778);

    auto g_0_zzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 779);

    auto g_0_zzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 780);

    auto g_0_zzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 781);

    auto g_0_zzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 782);

    auto g_0_zzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 783);

    /// Set up 0-28 components of targeted buffer : SKSI

    auto g_0_xxxxxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi);

    auto g_0_xxxxxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 1);

    auto g_0_xxxxxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 2);

    auto g_0_xxxxxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 3);

    auto g_0_xxxxxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 4);

    auto g_0_xxxxxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 5);

    auto g_0_xxxxxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 6);

    auto g_0_xxxxxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 7);

    auto g_0_xxxxxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 8);

    auto g_0_xxxxxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 9);

    auto g_0_xxxxxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 10);

    auto g_0_xxxxxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 11);

    auto g_0_xxxxxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 12);

    auto g_0_xxxxxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 13);

    auto g_0_xxxxxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 14);

    auto g_0_xxxxxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 15);

    auto g_0_xxxxxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 16);

    auto g_0_xxxxxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 17);

    auto g_0_xxxxxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 18);

    auto g_0_xxxxxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 19);

    auto g_0_xxxxxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 20);

    auto g_0_xxxxxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 21);

    auto g_0_xxxxxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 22);

    auto g_0_xxxxxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 23);

    auto g_0_xxxxxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 24);

    auto g_0_xxxxxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 25);

    auto g_0_xxxxxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 26);

    auto g_0_xxxxxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 27);

#pragma omp simd aligned(g_0_xxxxx_0_xxxxxx_0,       \
                             g_0_xxxxx_0_xxxxxx_1,   \
                             g_0_xxxxx_0_xxxxxy_0,   \
                             g_0_xxxxx_0_xxxxxy_1,   \
                             g_0_xxxxx_0_xxxxxz_0,   \
                             g_0_xxxxx_0_xxxxxz_1,   \
                             g_0_xxxxx_0_xxxxyy_0,   \
                             g_0_xxxxx_0_xxxxyy_1,   \
                             g_0_xxxxx_0_xxxxyz_0,   \
                             g_0_xxxxx_0_xxxxyz_1,   \
                             g_0_xxxxx_0_xxxxzz_0,   \
                             g_0_xxxxx_0_xxxxzz_1,   \
                             g_0_xxxxx_0_xxxyyy_0,   \
                             g_0_xxxxx_0_xxxyyy_1,   \
                             g_0_xxxxx_0_xxxyyz_0,   \
                             g_0_xxxxx_0_xxxyyz_1,   \
                             g_0_xxxxx_0_xxxyzz_0,   \
                             g_0_xxxxx_0_xxxyzz_1,   \
                             g_0_xxxxx_0_xxxzzz_0,   \
                             g_0_xxxxx_0_xxxzzz_1,   \
                             g_0_xxxxx_0_xxyyyy_0,   \
                             g_0_xxxxx_0_xxyyyy_1,   \
                             g_0_xxxxx_0_xxyyyz_0,   \
                             g_0_xxxxx_0_xxyyyz_1,   \
                             g_0_xxxxx_0_xxyyzz_0,   \
                             g_0_xxxxx_0_xxyyzz_1,   \
                             g_0_xxxxx_0_xxyzzz_0,   \
                             g_0_xxxxx_0_xxyzzz_1,   \
                             g_0_xxxxx_0_xxzzzz_0,   \
                             g_0_xxxxx_0_xxzzzz_1,   \
                             g_0_xxxxx_0_xyyyyy_0,   \
                             g_0_xxxxx_0_xyyyyy_1,   \
                             g_0_xxxxx_0_xyyyyz_0,   \
                             g_0_xxxxx_0_xyyyyz_1,   \
                             g_0_xxxxx_0_xyyyzz_0,   \
                             g_0_xxxxx_0_xyyyzz_1,   \
                             g_0_xxxxx_0_xyyzzz_0,   \
                             g_0_xxxxx_0_xyyzzz_1,   \
                             g_0_xxxxx_0_xyzzzz_0,   \
                             g_0_xxxxx_0_xyzzzz_1,   \
                             g_0_xxxxx_0_xzzzzz_0,   \
                             g_0_xxxxx_0_xzzzzz_1,   \
                             g_0_xxxxx_0_yyyyyy_0,   \
                             g_0_xxxxx_0_yyyyyy_1,   \
                             g_0_xxxxx_0_yyyyyz_0,   \
                             g_0_xxxxx_0_yyyyyz_1,   \
                             g_0_xxxxx_0_yyyyzz_0,   \
                             g_0_xxxxx_0_yyyyzz_1,   \
                             g_0_xxxxx_0_yyyzzz_0,   \
                             g_0_xxxxx_0_yyyzzz_1,   \
                             g_0_xxxxx_0_yyzzzz_0,   \
                             g_0_xxxxx_0_yyzzzz_1,   \
                             g_0_xxxxx_0_yzzzzz_0,   \
                             g_0_xxxxx_0_yzzzzz_1,   \
                             g_0_xxxxx_0_zzzzzz_0,   \
                             g_0_xxxxx_0_zzzzzz_1,   \
                             g_0_xxxxxx_0_xxxxx_1,   \
                             g_0_xxxxxx_0_xxxxxx_0,  \
                             g_0_xxxxxx_0_xxxxxx_1,  \
                             g_0_xxxxxx_0_xxxxxy_0,  \
                             g_0_xxxxxx_0_xxxxxy_1,  \
                             g_0_xxxxxx_0_xxxxxz_0,  \
                             g_0_xxxxxx_0_xxxxxz_1,  \
                             g_0_xxxxxx_0_xxxxy_1,   \
                             g_0_xxxxxx_0_xxxxyy_0,  \
                             g_0_xxxxxx_0_xxxxyy_1,  \
                             g_0_xxxxxx_0_xxxxyz_0,  \
                             g_0_xxxxxx_0_xxxxyz_1,  \
                             g_0_xxxxxx_0_xxxxz_1,   \
                             g_0_xxxxxx_0_xxxxzz_0,  \
                             g_0_xxxxxx_0_xxxxzz_1,  \
                             g_0_xxxxxx_0_xxxyy_1,   \
                             g_0_xxxxxx_0_xxxyyy_0,  \
                             g_0_xxxxxx_0_xxxyyy_1,  \
                             g_0_xxxxxx_0_xxxyyz_0,  \
                             g_0_xxxxxx_0_xxxyyz_1,  \
                             g_0_xxxxxx_0_xxxyz_1,   \
                             g_0_xxxxxx_0_xxxyzz_0,  \
                             g_0_xxxxxx_0_xxxyzz_1,  \
                             g_0_xxxxxx_0_xxxzz_1,   \
                             g_0_xxxxxx_0_xxxzzz_0,  \
                             g_0_xxxxxx_0_xxxzzz_1,  \
                             g_0_xxxxxx_0_xxyyy_1,   \
                             g_0_xxxxxx_0_xxyyyy_0,  \
                             g_0_xxxxxx_0_xxyyyy_1,  \
                             g_0_xxxxxx_0_xxyyyz_0,  \
                             g_0_xxxxxx_0_xxyyyz_1,  \
                             g_0_xxxxxx_0_xxyyz_1,   \
                             g_0_xxxxxx_0_xxyyzz_0,  \
                             g_0_xxxxxx_0_xxyyzz_1,  \
                             g_0_xxxxxx_0_xxyzz_1,   \
                             g_0_xxxxxx_0_xxyzzz_0,  \
                             g_0_xxxxxx_0_xxyzzz_1,  \
                             g_0_xxxxxx_0_xxzzz_1,   \
                             g_0_xxxxxx_0_xxzzzz_0,  \
                             g_0_xxxxxx_0_xxzzzz_1,  \
                             g_0_xxxxxx_0_xyyyy_1,   \
                             g_0_xxxxxx_0_xyyyyy_0,  \
                             g_0_xxxxxx_0_xyyyyy_1,  \
                             g_0_xxxxxx_0_xyyyyz_0,  \
                             g_0_xxxxxx_0_xyyyyz_1,  \
                             g_0_xxxxxx_0_xyyyz_1,   \
                             g_0_xxxxxx_0_xyyyzz_0,  \
                             g_0_xxxxxx_0_xyyyzz_1,  \
                             g_0_xxxxxx_0_xyyzz_1,   \
                             g_0_xxxxxx_0_xyyzzz_0,  \
                             g_0_xxxxxx_0_xyyzzz_1,  \
                             g_0_xxxxxx_0_xyzzz_1,   \
                             g_0_xxxxxx_0_xyzzzz_0,  \
                             g_0_xxxxxx_0_xyzzzz_1,  \
                             g_0_xxxxxx_0_xzzzz_1,   \
                             g_0_xxxxxx_0_xzzzzz_0,  \
                             g_0_xxxxxx_0_xzzzzz_1,  \
                             g_0_xxxxxx_0_yyyyy_1,   \
                             g_0_xxxxxx_0_yyyyyy_0,  \
                             g_0_xxxxxx_0_yyyyyy_1,  \
                             g_0_xxxxxx_0_yyyyyz_0,  \
                             g_0_xxxxxx_0_yyyyyz_1,  \
                             g_0_xxxxxx_0_yyyyz_1,   \
                             g_0_xxxxxx_0_yyyyzz_0,  \
                             g_0_xxxxxx_0_yyyyzz_1,  \
                             g_0_xxxxxx_0_yyyzz_1,   \
                             g_0_xxxxxx_0_yyyzzz_0,  \
                             g_0_xxxxxx_0_yyyzzz_1,  \
                             g_0_xxxxxx_0_yyzzz_1,   \
                             g_0_xxxxxx_0_yyzzzz_0,  \
                             g_0_xxxxxx_0_yyzzzz_1,  \
                             g_0_xxxxxx_0_yzzzz_1,   \
                             g_0_xxxxxx_0_yzzzzz_0,  \
                             g_0_xxxxxx_0_yzzzzz_1,  \
                             g_0_xxxxxx_0_zzzzz_1,   \
                             g_0_xxxxxx_0_zzzzzz_0,  \
                             g_0_xxxxxx_0_zzzzzz_1,  \
                             g_0_xxxxxxx_0_xxxxxx_0, \
                             g_0_xxxxxxx_0_xxxxxy_0, \
                             g_0_xxxxxxx_0_xxxxxz_0, \
                             g_0_xxxxxxx_0_xxxxyy_0, \
                             g_0_xxxxxxx_0_xxxxyz_0, \
                             g_0_xxxxxxx_0_xxxxzz_0, \
                             g_0_xxxxxxx_0_xxxyyy_0, \
                             g_0_xxxxxxx_0_xxxyyz_0, \
                             g_0_xxxxxxx_0_xxxyzz_0, \
                             g_0_xxxxxxx_0_xxxzzz_0, \
                             g_0_xxxxxxx_0_xxyyyy_0, \
                             g_0_xxxxxxx_0_xxyyyz_0, \
                             g_0_xxxxxxx_0_xxyyzz_0, \
                             g_0_xxxxxxx_0_xxyzzz_0, \
                             g_0_xxxxxxx_0_xxzzzz_0, \
                             g_0_xxxxxxx_0_xyyyyy_0, \
                             g_0_xxxxxxx_0_xyyyyz_0, \
                             g_0_xxxxxxx_0_xyyyzz_0, \
                             g_0_xxxxxxx_0_xyyzzz_0, \
                             g_0_xxxxxxx_0_xyzzzz_0, \
                             g_0_xxxxxxx_0_xzzzzz_0, \
                             g_0_xxxxxxx_0_yyyyyy_0, \
                             g_0_xxxxxxx_0_yyyyyz_0, \
                             g_0_xxxxxxx_0_yyyyzz_0, \
                             g_0_xxxxxxx_0_yyyzzz_0, \
                             g_0_xxxxxxx_0_yyzzzz_0, \
                             g_0_xxxxxxx_0_yzzzzz_0, \
                             g_0_xxxxxxx_0_zzzzzz_0, \
                             wp_x,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxx_0_xxxxxx_0[i] = 6.0 * g_0_xxxxx_0_xxxxxx_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxxx_1[i] * fti_ab_0 +
                                    6.0 * g_0_xxxxxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxx_0[i] * pb_x + g_0_xxxxxx_0_xxxxxx_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxxy_0[i] = 6.0 * g_0_xxxxx_0_xxxxxy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxxy_1[i] * fti_ab_0 +
                                    5.0 * g_0_xxxxxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxy_0[i] * pb_x + g_0_xxxxxx_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxxz_0[i] = 6.0 * g_0_xxxxx_0_xxxxxz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxxz_1[i] * fti_ab_0 +
                                    5.0 * g_0_xxxxxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxz_0[i] * pb_x + g_0_xxxxxx_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxyy_0[i] = 6.0 * g_0_xxxxx_0_xxxxyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxyy_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxxxxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyy_0[i] * pb_x + g_0_xxxxxx_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxyz_0[i] = 6.0 * g_0_xxxxx_0_xxxxyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxxxxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyz_0[i] * pb_x + g_0_xxxxxx_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxxzz_0[i] = 6.0 * g_0_xxxxx_0_xxxxzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxxzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxxxxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxzz_0[i] * pb_x + g_0_xxxxxx_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxyyy_0[i] = 6.0 * g_0_xxxxx_0_xxxyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxyyy_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxxxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyy_0[i] * pb_x + g_0_xxxxxx_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxyyz_0[i] = 6.0 * g_0_xxxxx_0_xxxyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxxxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyz_0[i] * pb_x + g_0_xxxxxx_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxyzz_0[i] = 6.0 * g_0_xxxxx_0_xxxyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxxxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyzz_0[i] * pb_x + g_0_xxxxxx_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxxzzz_0[i] = 6.0 * g_0_xxxxx_0_xxxzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxxzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxxxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxzzz_0[i] * pb_x + g_0_xxxxxx_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyyyy_0[i] = 6.0 * g_0_xxxxx_0_xxyyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyyyy_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyy_0[i] * pb_x + g_0_xxxxxx_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyyyz_0[i] = 6.0 * g_0_xxxxx_0_xxyyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyz_0[i] * pb_x + g_0_xxxxxx_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyyzz_0[i] = 6.0 * g_0_xxxxx_0_xxyyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyzz_0[i] * pb_x + g_0_xxxxxx_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxyzzz_0[i] = 6.0 * g_0_xxxxx_0_xxyzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyzzz_0[i] * pb_x + g_0_xxxxxx_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xxzzzz_0[i] = 6.0 * g_0_xxxxx_0_xxzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xxzzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxzzzz_0[i] * pb_x + g_0_xxxxxx_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyyyy_0[i] = 6.0 * g_0_xxxxx_0_xyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyy_0[i] * pb_x + g_0_xxxxxx_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyyyz_0[i] = 6.0 * g_0_xxxxx_0_xyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyz_0[i] * pb_x + g_0_xxxxxx_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyyzz_0[i] = 6.0 * g_0_xxxxx_0_xyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyzz_0[i] * pb_x + g_0_xxxxxx_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyyzzz_0[i] = 6.0 * g_0_xxxxx_0_xyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyzzz_0[i] * pb_x + g_0_xxxxxx_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xyzzzz_0[i] = 6.0 * g_0_xxxxx_0_xyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyzzzz_0[i] * pb_x + g_0_xxxxxx_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_xzzzzz_0[i] = 6.0 * g_0_xxxxx_0_xzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_xzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xzzzzz_0[i] * pb_x + g_0_xxxxxx_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyyyy_0[i] = 6.0 * g_0_xxxxx_0_yyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_yyyyyy_0[i] * pb_x + g_0_xxxxxx_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyyyz_0[i] = 6.0 * g_0_xxxxx_0_yyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_yyyyyz_0[i] * pb_x + g_0_xxxxxx_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyyzz_0[i] = 6.0 * g_0_xxxxx_0_yyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_yyyyzz_0[i] * pb_x + g_0_xxxxxx_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyyzzz_0[i] = 6.0 * g_0_xxxxx_0_yyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_yyyzzz_0[i] * pb_x + g_0_xxxxxx_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yyzzzz_0[i] = 6.0 * g_0_xxxxx_0_yyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_yyzzzz_0[i] * pb_x + g_0_xxxxxx_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_yzzzzz_0[i] = 6.0 * g_0_xxxxx_0_yzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_yzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_yzzzzz_0[i] * pb_x + g_0_xxxxxx_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxx_0_zzzzzz_0[i] = 6.0 * g_0_xxxxx_0_zzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_zzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxxx_0_zzzzzz_0[i] * pb_x + g_0_xxxxxx_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 28-56 components of targeted buffer : SKSI

    auto g_0_xxxxxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 28);

    auto g_0_xxxxxxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 29);

    auto g_0_xxxxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 30);

    auto g_0_xxxxxxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 31);

    auto g_0_xxxxxxy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 32);

    auto g_0_xxxxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 33);

    auto g_0_xxxxxxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 34);

    auto g_0_xxxxxxy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 35);

    auto g_0_xxxxxxy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 36);

    auto g_0_xxxxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 37);

    auto g_0_xxxxxxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 38);

    auto g_0_xxxxxxy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 39);

    auto g_0_xxxxxxy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 40);

    auto g_0_xxxxxxy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 41);

    auto g_0_xxxxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 42);

    auto g_0_xxxxxxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 43);

    auto g_0_xxxxxxy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 44);

    auto g_0_xxxxxxy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 45);

    auto g_0_xxxxxxy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 46);

    auto g_0_xxxxxxy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 47);

    auto g_0_xxxxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 48);

    auto g_0_xxxxxxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 49);

    auto g_0_xxxxxxy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 50);

    auto g_0_xxxxxxy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 51);

    auto g_0_xxxxxxy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 52);

    auto g_0_xxxxxxy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 53);

    auto g_0_xxxxxxy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 54);

    auto g_0_xxxxxxy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 55);

#pragma omp simd aligned(g_0_xxxxxx_0_xxxxx_1,       \
                             g_0_xxxxxx_0_xxxxxx_0,  \
                             g_0_xxxxxx_0_xxxxxx_1,  \
                             g_0_xxxxxx_0_xxxxxy_0,  \
                             g_0_xxxxxx_0_xxxxxy_1,  \
                             g_0_xxxxxx_0_xxxxxz_0,  \
                             g_0_xxxxxx_0_xxxxxz_1,  \
                             g_0_xxxxxx_0_xxxxy_1,   \
                             g_0_xxxxxx_0_xxxxyy_0,  \
                             g_0_xxxxxx_0_xxxxyy_1,  \
                             g_0_xxxxxx_0_xxxxyz_0,  \
                             g_0_xxxxxx_0_xxxxyz_1,  \
                             g_0_xxxxxx_0_xxxxz_1,   \
                             g_0_xxxxxx_0_xxxxzz_0,  \
                             g_0_xxxxxx_0_xxxxzz_1,  \
                             g_0_xxxxxx_0_xxxyy_1,   \
                             g_0_xxxxxx_0_xxxyyy_0,  \
                             g_0_xxxxxx_0_xxxyyy_1,  \
                             g_0_xxxxxx_0_xxxyyz_0,  \
                             g_0_xxxxxx_0_xxxyyz_1,  \
                             g_0_xxxxxx_0_xxxyz_1,   \
                             g_0_xxxxxx_0_xxxyzz_0,  \
                             g_0_xxxxxx_0_xxxyzz_1,  \
                             g_0_xxxxxx_0_xxxzz_1,   \
                             g_0_xxxxxx_0_xxxzzz_0,  \
                             g_0_xxxxxx_0_xxxzzz_1,  \
                             g_0_xxxxxx_0_xxyyy_1,   \
                             g_0_xxxxxx_0_xxyyyy_0,  \
                             g_0_xxxxxx_0_xxyyyy_1,  \
                             g_0_xxxxxx_0_xxyyyz_0,  \
                             g_0_xxxxxx_0_xxyyyz_1,  \
                             g_0_xxxxxx_0_xxyyz_1,   \
                             g_0_xxxxxx_0_xxyyzz_0,  \
                             g_0_xxxxxx_0_xxyyzz_1,  \
                             g_0_xxxxxx_0_xxyzz_1,   \
                             g_0_xxxxxx_0_xxyzzz_0,  \
                             g_0_xxxxxx_0_xxyzzz_1,  \
                             g_0_xxxxxx_0_xxzzz_1,   \
                             g_0_xxxxxx_0_xxzzzz_0,  \
                             g_0_xxxxxx_0_xxzzzz_1,  \
                             g_0_xxxxxx_0_xyyyy_1,   \
                             g_0_xxxxxx_0_xyyyyy_0,  \
                             g_0_xxxxxx_0_xyyyyy_1,  \
                             g_0_xxxxxx_0_xyyyyz_0,  \
                             g_0_xxxxxx_0_xyyyyz_1,  \
                             g_0_xxxxxx_0_xyyyz_1,   \
                             g_0_xxxxxx_0_xyyyzz_0,  \
                             g_0_xxxxxx_0_xyyyzz_1,  \
                             g_0_xxxxxx_0_xyyzz_1,   \
                             g_0_xxxxxx_0_xyyzzz_0,  \
                             g_0_xxxxxx_0_xyyzzz_1,  \
                             g_0_xxxxxx_0_xyzzz_1,   \
                             g_0_xxxxxx_0_xyzzzz_0,  \
                             g_0_xxxxxx_0_xyzzzz_1,  \
                             g_0_xxxxxx_0_xzzzz_1,   \
                             g_0_xxxxxx_0_xzzzzz_0,  \
                             g_0_xxxxxx_0_xzzzzz_1,  \
                             g_0_xxxxxx_0_yyyyy_1,   \
                             g_0_xxxxxx_0_yyyyyy_0,  \
                             g_0_xxxxxx_0_yyyyyy_1,  \
                             g_0_xxxxxx_0_yyyyyz_0,  \
                             g_0_xxxxxx_0_yyyyyz_1,  \
                             g_0_xxxxxx_0_yyyyz_1,   \
                             g_0_xxxxxx_0_yyyyzz_0,  \
                             g_0_xxxxxx_0_yyyyzz_1,  \
                             g_0_xxxxxx_0_yyyzz_1,   \
                             g_0_xxxxxx_0_yyyzzz_0,  \
                             g_0_xxxxxx_0_yyyzzz_1,  \
                             g_0_xxxxxx_0_yyzzz_1,   \
                             g_0_xxxxxx_0_yyzzzz_0,  \
                             g_0_xxxxxx_0_yyzzzz_1,  \
                             g_0_xxxxxx_0_yzzzz_1,   \
                             g_0_xxxxxx_0_yzzzzz_0,  \
                             g_0_xxxxxx_0_yzzzzz_1,  \
                             g_0_xxxxxx_0_zzzzz_1,   \
                             g_0_xxxxxx_0_zzzzzz_0,  \
                             g_0_xxxxxx_0_zzzzzz_1,  \
                             g_0_xxxxxxy_0_xxxxxx_0, \
                             g_0_xxxxxxy_0_xxxxxy_0, \
                             g_0_xxxxxxy_0_xxxxxz_0, \
                             g_0_xxxxxxy_0_xxxxyy_0, \
                             g_0_xxxxxxy_0_xxxxyz_0, \
                             g_0_xxxxxxy_0_xxxxzz_0, \
                             g_0_xxxxxxy_0_xxxyyy_0, \
                             g_0_xxxxxxy_0_xxxyyz_0, \
                             g_0_xxxxxxy_0_xxxyzz_0, \
                             g_0_xxxxxxy_0_xxxzzz_0, \
                             g_0_xxxxxxy_0_xxyyyy_0, \
                             g_0_xxxxxxy_0_xxyyyz_0, \
                             g_0_xxxxxxy_0_xxyyzz_0, \
                             g_0_xxxxxxy_0_xxyzzz_0, \
                             g_0_xxxxxxy_0_xxzzzz_0, \
                             g_0_xxxxxxy_0_xyyyyy_0, \
                             g_0_xxxxxxy_0_xyyyyz_0, \
                             g_0_xxxxxxy_0_xyyyzz_0, \
                             g_0_xxxxxxy_0_xyyzzz_0, \
                             g_0_xxxxxxy_0_xyzzzz_0, \
                             g_0_xxxxxxy_0_xzzzzz_0, \
                             g_0_xxxxxxy_0_yyyyyy_0, \
                             g_0_xxxxxxy_0_yyyyyz_0, \
                             g_0_xxxxxxy_0_yyyyzz_0, \
                             g_0_xxxxxxy_0_yyyzzz_0, \
                             g_0_xxxxxxy_0_yyzzzz_0, \
                             g_0_xxxxxxy_0_yzzzzz_0, \
                             g_0_xxxxxxy_0_zzzzzz_0, \
                             wp_y,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxy_0_xxxxxx_0[i] = g_0_xxxxxx_0_xxxxxx_0[i] * pb_y + g_0_xxxxxx_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxxy_0[i] = g_0_xxxxxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxy_0[i] * pb_y + g_0_xxxxxx_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxxz_0[i] = g_0_xxxxxx_0_xxxxxz_0[i] * pb_y + g_0_xxxxxx_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxyy_0[i] = 2.0 * g_0_xxxxxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyy_0[i] * pb_y + g_0_xxxxxx_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxyz_0[i] = g_0_xxxxxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyz_0[i] * pb_y + g_0_xxxxxx_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxxzz_0[i] = g_0_xxxxxx_0_xxxxzz_0[i] * pb_y + g_0_xxxxxx_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxyyy_0[i] = 3.0 * g_0_xxxxxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyy_0[i] * pb_y + g_0_xxxxxx_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxyyz_0[i] = 2.0 * g_0_xxxxxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyz_0[i] * pb_y + g_0_xxxxxx_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxyzz_0[i] = g_0_xxxxxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyzz_0[i] * pb_y + g_0_xxxxxx_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxxzzz_0[i] = g_0_xxxxxx_0_xxxzzz_0[i] * pb_y + g_0_xxxxxx_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyyyy_0[i] = 4.0 * g_0_xxxxxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyy_0[i] * pb_y + g_0_xxxxxx_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyyyz_0[i] = 3.0 * g_0_xxxxxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyz_0[i] * pb_y + g_0_xxxxxx_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyyzz_0[i] = 2.0 * g_0_xxxxxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyzz_0[i] * pb_y + g_0_xxxxxx_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxyzzz_0[i] = g_0_xxxxxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyzzz_0[i] * pb_y + g_0_xxxxxx_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xxzzzz_0[i] = g_0_xxxxxx_0_xxzzzz_0[i] * pb_y + g_0_xxxxxx_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyyyy_0[i] = 5.0 * g_0_xxxxxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyy_0[i] * pb_y + g_0_xxxxxx_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyyyz_0[i] = 4.0 * g_0_xxxxxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyz_0[i] * pb_y + g_0_xxxxxx_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyyzz_0[i] = 3.0 * g_0_xxxxxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyzz_0[i] * pb_y + g_0_xxxxxx_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyyzzz_0[i] = 2.0 * g_0_xxxxxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyzzz_0[i] * pb_y + g_0_xxxxxx_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xyzzzz_0[i] = g_0_xxxxxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyzzzz_0[i] * pb_y + g_0_xxxxxx_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_xzzzzz_0[i] = g_0_xxxxxx_0_xzzzzz_0[i] * pb_y + g_0_xxxxxx_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyyyy_0[i] = 6.0 * g_0_xxxxxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyyy_0[i] * pb_y + g_0_xxxxxx_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyyyz_0[i] = 5.0 * g_0_xxxxxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyyz_0[i] * pb_y + g_0_xxxxxx_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyyzz_0[i] = 4.0 * g_0_xxxxxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyzz_0[i] * pb_y + g_0_xxxxxx_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyyzzz_0[i] = 3.0 * g_0_xxxxxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyzzz_0[i] * pb_y + g_0_xxxxxx_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yyzzzz_0[i] = 2.0 * g_0_xxxxxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyzzzz_0[i] * pb_y + g_0_xxxxxx_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_yzzzzz_0[i] = g_0_xxxxxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yzzzzz_0[i] * pb_y + g_0_xxxxxx_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxy_0_zzzzzz_0[i] = g_0_xxxxxx_0_zzzzzz_0[i] * pb_y + g_0_xxxxxx_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 56-84 components of targeted buffer : SKSI

    auto g_0_xxxxxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 56);

    auto g_0_xxxxxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 57);

    auto g_0_xxxxxxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 58);

    auto g_0_xxxxxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 59);

    auto g_0_xxxxxxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 60);

    auto g_0_xxxxxxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 61);

    auto g_0_xxxxxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 62);

    auto g_0_xxxxxxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 63);

    auto g_0_xxxxxxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 64);

    auto g_0_xxxxxxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 65);

    auto g_0_xxxxxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 66);

    auto g_0_xxxxxxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 67);

    auto g_0_xxxxxxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 68);

    auto g_0_xxxxxxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 69);

    auto g_0_xxxxxxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 70);

    auto g_0_xxxxxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 71);

    auto g_0_xxxxxxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 72);

    auto g_0_xxxxxxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 73);

    auto g_0_xxxxxxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 74);

    auto g_0_xxxxxxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 75);

    auto g_0_xxxxxxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 76);

    auto g_0_xxxxxxz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 77);

    auto g_0_xxxxxxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 78);

    auto g_0_xxxxxxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 79);

    auto g_0_xxxxxxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 80);

    auto g_0_xxxxxxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 81);

    auto g_0_xxxxxxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 82);

    auto g_0_xxxxxxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 83);

#pragma omp simd aligned(g_0_xxxxxx_0_xxxxx_1,       \
                             g_0_xxxxxx_0_xxxxxx_0,  \
                             g_0_xxxxxx_0_xxxxxx_1,  \
                             g_0_xxxxxx_0_xxxxxy_0,  \
                             g_0_xxxxxx_0_xxxxxy_1,  \
                             g_0_xxxxxx_0_xxxxxz_0,  \
                             g_0_xxxxxx_0_xxxxxz_1,  \
                             g_0_xxxxxx_0_xxxxy_1,   \
                             g_0_xxxxxx_0_xxxxyy_0,  \
                             g_0_xxxxxx_0_xxxxyy_1,  \
                             g_0_xxxxxx_0_xxxxyz_0,  \
                             g_0_xxxxxx_0_xxxxyz_1,  \
                             g_0_xxxxxx_0_xxxxz_1,   \
                             g_0_xxxxxx_0_xxxxzz_0,  \
                             g_0_xxxxxx_0_xxxxzz_1,  \
                             g_0_xxxxxx_0_xxxyy_1,   \
                             g_0_xxxxxx_0_xxxyyy_0,  \
                             g_0_xxxxxx_0_xxxyyy_1,  \
                             g_0_xxxxxx_0_xxxyyz_0,  \
                             g_0_xxxxxx_0_xxxyyz_1,  \
                             g_0_xxxxxx_0_xxxyz_1,   \
                             g_0_xxxxxx_0_xxxyzz_0,  \
                             g_0_xxxxxx_0_xxxyzz_1,  \
                             g_0_xxxxxx_0_xxxzz_1,   \
                             g_0_xxxxxx_0_xxxzzz_0,  \
                             g_0_xxxxxx_0_xxxzzz_1,  \
                             g_0_xxxxxx_0_xxyyy_1,   \
                             g_0_xxxxxx_0_xxyyyy_0,  \
                             g_0_xxxxxx_0_xxyyyy_1,  \
                             g_0_xxxxxx_0_xxyyyz_0,  \
                             g_0_xxxxxx_0_xxyyyz_1,  \
                             g_0_xxxxxx_0_xxyyz_1,   \
                             g_0_xxxxxx_0_xxyyzz_0,  \
                             g_0_xxxxxx_0_xxyyzz_1,  \
                             g_0_xxxxxx_0_xxyzz_1,   \
                             g_0_xxxxxx_0_xxyzzz_0,  \
                             g_0_xxxxxx_0_xxyzzz_1,  \
                             g_0_xxxxxx_0_xxzzz_1,   \
                             g_0_xxxxxx_0_xxzzzz_0,  \
                             g_0_xxxxxx_0_xxzzzz_1,  \
                             g_0_xxxxxx_0_xyyyy_1,   \
                             g_0_xxxxxx_0_xyyyyy_0,  \
                             g_0_xxxxxx_0_xyyyyy_1,  \
                             g_0_xxxxxx_0_xyyyyz_0,  \
                             g_0_xxxxxx_0_xyyyyz_1,  \
                             g_0_xxxxxx_0_xyyyz_1,   \
                             g_0_xxxxxx_0_xyyyzz_0,  \
                             g_0_xxxxxx_0_xyyyzz_1,  \
                             g_0_xxxxxx_0_xyyzz_1,   \
                             g_0_xxxxxx_0_xyyzzz_0,  \
                             g_0_xxxxxx_0_xyyzzz_1,  \
                             g_0_xxxxxx_0_xyzzz_1,   \
                             g_0_xxxxxx_0_xyzzzz_0,  \
                             g_0_xxxxxx_0_xyzzzz_1,  \
                             g_0_xxxxxx_0_xzzzz_1,   \
                             g_0_xxxxxx_0_xzzzzz_0,  \
                             g_0_xxxxxx_0_xzzzzz_1,  \
                             g_0_xxxxxx_0_yyyyy_1,   \
                             g_0_xxxxxx_0_yyyyyy_0,  \
                             g_0_xxxxxx_0_yyyyyy_1,  \
                             g_0_xxxxxx_0_yyyyyz_0,  \
                             g_0_xxxxxx_0_yyyyyz_1,  \
                             g_0_xxxxxx_0_yyyyz_1,   \
                             g_0_xxxxxx_0_yyyyzz_0,  \
                             g_0_xxxxxx_0_yyyyzz_1,  \
                             g_0_xxxxxx_0_yyyzz_1,   \
                             g_0_xxxxxx_0_yyyzzz_0,  \
                             g_0_xxxxxx_0_yyyzzz_1,  \
                             g_0_xxxxxx_0_yyzzz_1,   \
                             g_0_xxxxxx_0_yyzzzz_0,  \
                             g_0_xxxxxx_0_yyzzzz_1,  \
                             g_0_xxxxxx_0_yzzzz_1,   \
                             g_0_xxxxxx_0_yzzzzz_0,  \
                             g_0_xxxxxx_0_yzzzzz_1,  \
                             g_0_xxxxxx_0_zzzzz_1,   \
                             g_0_xxxxxx_0_zzzzzz_0,  \
                             g_0_xxxxxx_0_zzzzzz_1,  \
                             g_0_xxxxxxz_0_xxxxxx_0, \
                             g_0_xxxxxxz_0_xxxxxy_0, \
                             g_0_xxxxxxz_0_xxxxxz_0, \
                             g_0_xxxxxxz_0_xxxxyy_0, \
                             g_0_xxxxxxz_0_xxxxyz_0, \
                             g_0_xxxxxxz_0_xxxxzz_0, \
                             g_0_xxxxxxz_0_xxxyyy_0, \
                             g_0_xxxxxxz_0_xxxyyz_0, \
                             g_0_xxxxxxz_0_xxxyzz_0, \
                             g_0_xxxxxxz_0_xxxzzz_0, \
                             g_0_xxxxxxz_0_xxyyyy_0, \
                             g_0_xxxxxxz_0_xxyyyz_0, \
                             g_0_xxxxxxz_0_xxyyzz_0, \
                             g_0_xxxxxxz_0_xxyzzz_0, \
                             g_0_xxxxxxz_0_xxzzzz_0, \
                             g_0_xxxxxxz_0_xyyyyy_0, \
                             g_0_xxxxxxz_0_xyyyyz_0, \
                             g_0_xxxxxxz_0_xyyyzz_0, \
                             g_0_xxxxxxz_0_xyyzzz_0, \
                             g_0_xxxxxxz_0_xyzzzz_0, \
                             g_0_xxxxxxz_0_xzzzzz_0, \
                             g_0_xxxxxxz_0_yyyyyy_0, \
                             g_0_xxxxxxz_0_yyyyyz_0, \
                             g_0_xxxxxxz_0_yyyyzz_0, \
                             g_0_xxxxxxz_0_yyyzzz_0, \
                             g_0_xxxxxxz_0_yyzzzz_0, \
                             g_0_xxxxxxz_0_yzzzzz_0, \
                             g_0_xxxxxxz_0_zzzzzz_0, \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxz_0_xxxxxx_0[i] = g_0_xxxxxx_0_xxxxxx_0[i] * pb_z + g_0_xxxxxx_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxxy_0[i] = g_0_xxxxxx_0_xxxxxy_0[i] * pb_z + g_0_xxxxxx_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxxz_0[i] = g_0_xxxxxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxxz_0[i] * pb_z + g_0_xxxxxx_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxyy_0[i] = g_0_xxxxxx_0_xxxxyy_0[i] * pb_z + g_0_xxxxxx_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxyz_0[i] = g_0_xxxxxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxyz_0[i] * pb_z + g_0_xxxxxx_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxxzz_0[i] = 2.0 * g_0_xxxxxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxxzz_0[i] * pb_z + g_0_xxxxxx_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxyyy_0[i] = g_0_xxxxxx_0_xxxyyy_0[i] * pb_z + g_0_xxxxxx_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxyyz_0[i] = g_0_xxxxxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyyz_0[i] * pb_z + g_0_xxxxxx_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxyzz_0[i] = 2.0 * g_0_xxxxxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxyzz_0[i] * pb_z + g_0_xxxxxx_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxxzzz_0[i] = 3.0 * g_0_xxxxxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxxzzz_0[i] * pb_z + g_0_xxxxxx_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyyyy_0[i] = g_0_xxxxxx_0_xxyyyy_0[i] * pb_z + g_0_xxxxxx_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyyyz_0[i] = g_0_xxxxxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyyz_0[i] * pb_z + g_0_xxxxxx_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyyzz_0[i] = 2.0 * g_0_xxxxxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyyzz_0[i] * pb_z + g_0_xxxxxx_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxyzzz_0[i] = 3.0 * g_0_xxxxxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxyzzz_0[i] * pb_z + g_0_xxxxxx_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xxzzzz_0[i] = 4.0 * g_0_xxxxxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xxzzzz_0[i] * pb_z + g_0_xxxxxx_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyyyy_0[i] = g_0_xxxxxx_0_xyyyyy_0[i] * pb_z + g_0_xxxxxx_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyyyz_0[i] = g_0_xxxxxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyyz_0[i] * pb_z + g_0_xxxxxx_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyyzz_0[i] = 2.0 * g_0_xxxxxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyyzz_0[i] * pb_z + g_0_xxxxxx_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyyzzz_0[i] = 3.0 * g_0_xxxxxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyyzzz_0[i] * pb_z + g_0_xxxxxx_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xyzzzz_0[i] = 4.0 * g_0_xxxxxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xyzzzz_0[i] * pb_z + g_0_xxxxxx_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_xzzzzz_0[i] = 5.0 * g_0_xxxxxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_xzzzzz_0[i] * pb_z + g_0_xxxxxx_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyyyy_0[i] = g_0_xxxxxx_0_yyyyyy_0[i] * pb_z + g_0_xxxxxx_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyyyz_0[i] = g_0_xxxxxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyyz_0[i] * pb_z + g_0_xxxxxx_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyyzz_0[i] = 2.0 * g_0_xxxxxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyyzz_0[i] * pb_z + g_0_xxxxxx_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyyzzz_0[i] = 3.0 * g_0_xxxxxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyyzzz_0[i] * pb_z + g_0_xxxxxx_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yyzzzz_0[i] = 4.0 * g_0_xxxxxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yyzzzz_0[i] * pb_z + g_0_xxxxxx_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_yzzzzz_0[i] = 5.0 * g_0_xxxxxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_yzzzzz_0[i] * pb_z + g_0_xxxxxx_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxz_0_zzzzzz_0[i] = 6.0 * g_0_xxxxxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxxx_0_zzzzzz_0[i] * pb_z + g_0_xxxxxx_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 84-112 components of targeted buffer : SKSI

    auto g_0_xxxxxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 84);

    auto g_0_xxxxxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 85);

    auto g_0_xxxxxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 86);

    auto g_0_xxxxxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 87);

    auto g_0_xxxxxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 88);

    auto g_0_xxxxxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 89);

    auto g_0_xxxxxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 90);

    auto g_0_xxxxxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 91);

    auto g_0_xxxxxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 92);

    auto g_0_xxxxxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 93);

    auto g_0_xxxxxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 94);

    auto g_0_xxxxxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 95);

    auto g_0_xxxxxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 96);

    auto g_0_xxxxxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 97);

    auto g_0_xxxxxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 98);

    auto g_0_xxxxxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 99);

    auto g_0_xxxxxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 100);

    auto g_0_xxxxxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 101);

    auto g_0_xxxxxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 102);

    auto g_0_xxxxxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 103);

    auto g_0_xxxxxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 104);

    auto g_0_xxxxxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 105);

    auto g_0_xxxxxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 106);

    auto g_0_xxxxxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 107);

    auto g_0_xxxxxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 108);

    auto g_0_xxxxxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 109);

    auto g_0_xxxxxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 110);

    auto g_0_xxxxxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 111);

#pragma omp simd aligned(g_0_xxxxx_0_xxxxxx_0,       \
                             g_0_xxxxx_0_xxxxxx_1,   \
                             g_0_xxxxx_0_xxxxxz_0,   \
                             g_0_xxxxx_0_xxxxxz_1,   \
                             g_0_xxxxx_0_xxxxzz_0,   \
                             g_0_xxxxx_0_xxxxzz_1,   \
                             g_0_xxxxx_0_xxxzzz_0,   \
                             g_0_xxxxx_0_xxxzzz_1,   \
                             g_0_xxxxx_0_xxzzzz_0,   \
                             g_0_xxxxx_0_xxzzzz_1,   \
                             g_0_xxxxx_0_xzzzzz_0,   \
                             g_0_xxxxx_0_xzzzzz_1,   \
                             g_0_xxxxxy_0_xxxxxx_0,  \
                             g_0_xxxxxy_0_xxxxxx_1,  \
                             g_0_xxxxxy_0_xxxxxz_0,  \
                             g_0_xxxxxy_0_xxxxxz_1,  \
                             g_0_xxxxxy_0_xxxxzz_0,  \
                             g_0_xxxxxy_0_xxxxzz_1,  \
                             g_0_xxxxxy_0_xxxzzz_0,  \
                             g_0_xxxxxy_0_xxxzzz_1,  \
                             g_0_xxxxxy_0_xxzzzz_0,  \
                             g_0_xxxxxy_0_xxzzzz_1,  \
                             g_0_xxxxxy_0_xzzzzz_0,  \
                             g_0_xxxxxy_0_xzzzzz_1,  \
                             g_0_xxxxxyy_0_xxxxxx_0, \
                             g_0_xxxxxyy_0_xxxxxy_0, \
                             g_0_xxxxxyy_0_xxxxxz_0, \
                             g_0_xxxxxyy_0_xxxxyy_0, \
                             g_0_xxxxxyy_0_xxxxyz_0, \
                             g_0_xxxxxyy_0_xxxxzz_0, \
                             g_0_xxxxxyy_0_xxxyyy_0, \
                             g_0_xxxxxyy_0_xxxyyz_0, \
                             g_0_xxxxxyy_0_xxxyzz_0, \
                             g_0_xxxxxyy_0_xxxzzz_0, \
                             g_0_xxxxxyy_0_xxyyyy_0, \
                             g_0_xxxxxyy_0_xxyyyz_0, \
                             g_0_xxxxxyy_0_xxyyzz_0, \
                             g_0_xxxxxyy_0_xxyzzz_0, \
                             g_0_xxxxxyy_0_xxzzzz_0, \
                             g_0_xxxxxyy_0_xyyyyy_0, \
                             g_0_xxxxxyy_0_xyyyyz_0, \
                             g_0_xxxxxyy_0_xyyyzz_0, \
                             g_0_xxxxxyy_0_xyyzzz_0, \
                             g_0_xxxxxyy_0_xyzzzz_0, \
                             g_0_xxxxxyy_0_xzzzzz_0, \
                             g_0_xxxxxyy_0_yyyyyy_0, \
                             g_0_xxxxxyy_0_yyyyyz_0, \
                             g_0_xxxxxyy_0_yyyyzz_0, \
                             g_0_xxxxxyy_0_yyyzzz_0, \
                             g_0_xxxxxyy_0_yyzzzz_0, \
                             g_0_xxxxxyy_0_yzzzzz_0, \
                             g_0_xxxxxyy_0_zzzzzz_0, \
                             g_0_xxxxyy_0_xxxxxy_0,  \
                             g_0_xxxxyy_0_xxxxxy_1,  \
                             g_0_xxxxyy_0_xxxxy_1,   \
                             g_0_xxxxyy_0_xxxxyy_0,  \
                             g_0_xxxxyy_0_xxxxyy_1,  \
                             g_0_xxxxyy_0_xxxxyz_0,  \
                             g_0_xxxxyy_0_xxxxyz_1,  \
                             g_0_xxxxyy_0_xxxyy_1,   \
                             g_0_xxxxyy_0_xxxyyy_0,  \
                             g_0_xxxxyy_0_xxxyyy_1,  \
                             g_0_xxxxyy_0_xxxyyz_0,  \
                             g_0_xxxxyy_0_xxxyyz_1,  \
                             g_0_xxxxyy_0_xxxyz_1,   \
                             g_0_xxxxyy_0_xxxyzz_0,  \
                             g_0_xxxxyy_0_xxxyzz_1,  \
                             g_0_xxxxyy_0_xxyyy_1,   \
                             g_0_xxxxyy_0_xxyyyy_0,  \
                             g_0_xxxxyy_0_xxyyyy_1,  \
                             g_0_xxxxyy_0_xxyyyz_0,  \
                             g_0_xxxxyy_0_xxyyyz_1,  \
                             g_0_xxxxyy_0_xxyyz_1,   \
                             g_0_xxxxyy_0_xxyyzz_0,  \
                             g_0_xxxxyy_0_xxyyzz_1,  \
                             g_0_xxxxyy_0_xxyzz_1,   \
                             g_0_xxxxyy_0_xxyzzz_0,  \
                             g_0_xxxxyy_0_xxyzzz_1,  \
                             g_0_xxxxyy_0_xyyyy_1,   \
                             g_0_xxxxyy_0_xyyyyy_0,  \
                             g_0_xxxxyy_0_xyyyyy_1,  \
                             g_0_xxxxyy_0_xyyyyz_0,  \
                             g_0_xxxxyy_0_xyyyyz_1,  \
                             g_0_xxxxyy_0_xyyyz_1,   \
                             g_0_xxxxyy_0_xyyyzz_0,  \
                             g_0_xxxxyy_0_xyyyzz_1,  \
                             g_0_xxxxyy_0_xyyzz_1,   \
                             g_0_xxxxyy_0_xyyzzz_0,  \
                             g_0_xxxxyy_0_xyyzzz_1,  \
                             g_0_xxxxyy_0_xyzzz_1,   \
                             g_0_xxxxyy_0_xyzzzz_0,  \
                             g_0_xxxxyy_0_xyzzzz_1,  \
                             g_0_xxxxyy_0_yyyyy_1,   \
                             g_0_xxxxyy_0_yyyyyy_0,  \
                             g_0_xxxxyy_0_yyyyyy_1,  \
                             g_0_xxxxyy_0_yyyyyz_0,  \
                             g_0_xxxxyy_0_yyyyyz_1,  \
                             g_0_xxxxyy_0_yyyyz_1,   \
                             g_0_xxxxyy_0_yyyyzz_0,  \
                             g_0_xxxxyy_0_yyyyzz_1,  \
                             g_0_xxxxyy_0_yyyzz_1,   \
                             g_0_xxxxyy_0_yyyzzz_0,  \
                             g_0_xxxxyy_0_yyyzzz_1,  \
                             g_0_xxxxyy_0_yyzzz_1,   \
                             g_0_xxxxyy_0_yyzzzz_0,  \
                             g_0_xxxxyy_0_yyzzzz_1,  \
                             g_0_xxxxyy_0_yzzzz_1,   \
                             g_0_xxxxyy_0_yzzzzz_0,  \
                             g_0_xxxxyy_0_yzzzzz_1,  \
                             g_0_xxxxyy_0_zzzzzz_0,  \
                             g_0_xxxxyy_0_zzzzzz_1,  \
                             g_0_xxxyy_0_xxxxxy_0,   \
                             g_0_xxxyy_0_xxxxxy_1,   \
                             g_0_xxxyy_0_xxxxyy_0,   \
                             g_0_xxxyy_0_xxxxyy_1,   \
                             g_0_xxxyy_0_xxxxyz_0,   \
                             g_0_xxxyy_0_xxxxyz_1,   \
                             g_0_xxxyy_0_xxxyyy_0,   \
                             g_0_xxxyy_0_xxxyyy_1,   \
                             g_0_xxxyy_0_xxxyyz_0,   \
                             g_0_xxxyy_0_xxxyyz_1,   \
                             g_0_xxxyy_0_xxxyzz_0,   \
                             g_0_xxxyy_0_xxxyzz_1,   \
                             g_0_xxxyy_0_xxyyyy_0,   \
                             g_0_xxxyy_0_xxyyyy_1,   \
                             g_0_xxxyy_0_xxyyyz_0,   \
                             g_0_xxxyy_0_xxyyyz_1,   \
                             g_0_xxxyy_0_xxyyzz_0,   \
                             g_0_xxxyy_0_xxyyzz_1,   \
                             g_0_xxxyy_0_xxyzzz_0,   \
                             g_0_xxxyy_0_xxyzzz_1,   \
                             g_0_xxxyy_0_xyyyyy_0,   \
                             g_0_xxxyy_0_xyyyyy_1,   \
                             g_0_xxxyy_0_xyyyyz_0,   \
                             g_0_xxxyy_0_xyyyyz_1,   \
                             g_0_xxxyy_0_xyyyzz_0,   \
                             g_0_xxxyy_0_xyyyzz_1,   \
                             g_0_xxxyy_0_xyyzzz_0,   \
                             g_0_xxxyy_0_xyyzzz_1,   \
                             g_0_xxxyy_0_xyzzzz_0,   \
                             g_0_xxxyy_0_xyzzzz_1,   \
                             g_0_xxxyy_0_yyyyyy_0,   \
                             g_0_xxxyy_0_yyyyyy_1,   \
                             g_0_xxxyy_0_yyyyyz_0,   \
                             g_0_xxxyy_0_yyyyyz_1,   \
                             g_0_xxxyy_0_yyyyzz_0,   \
                             g_0_xxxyy_0_yyyyzz_1,   \
                             g_0_xxxyy_0_yyyzzz_0,   \
                             g_0_xxxyy_0_yyyzzz_1,   \
                             g_0_xxxyy_0_yyzzzz_0,   \
                             g_0_xxxyy_0_yyzzzz_1,   \
                             g_0_xxxyy_0_yzzzzz_0,   \
                             g_0_xxxyy_0_yzzzzz_1,   \
                             g_0_xxxyy_0_zzzzzz_0,   \
                             g_0_xxxyy_0_zzzzzz_1,   \
                             wp_x,                   \
                             wp_y,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxyy_0_xxxxxx_0[i] = g_0_xxxxx_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxxxx_0[i] * pb_y +
                                    g_0_xxxxxy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxxxxy_0[i] = 4.0 * g_0_xxxyy_0_xxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxxxy_1[i] * fti_ab_0 +
                                    5.0 * g_0_xxxxyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxxy_0[i] * pb_x + g_0_xxxxyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxxxz_0[i] = g_0_xxxxx_0_xxxxxz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxxxz_0[i] * pb_y +
                                    g_0_xxxxxy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxxxyy_0[i] = 4.0 * g_0_xxxyy_0_xxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxxyy_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxxxyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxyy_0[i] * pb_x + g_0_xxxxyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxxyz_0[i] = 4.0 * g_0_xxxyy_0_xxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxxyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxxxyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxyz_0[i] * pb_x + g_0_xxxxyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxxzz_0[i] = g_0_xxxxx_0_xxxxzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxxzz_0[i] * pb_y +
                                    g_0_xxxxxy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxxyyy_0[i] = 4.0 * g_0_xxxyy_0_xxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxyyy_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxxyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyyy_0[i] * pb_x + g_0_xxxxyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxyyz_0[i] = 4.0 * g_0_xxxyy_0_xxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxxyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyyz_0[i] * pb_x + g_0_xxxxyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxyzz_0[i] = 4.0 * g_0_xxxyy_0_xxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxxyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxxyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyzz_0[i] * pb_x + g_0_xxxxyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxxzzz_0[i] = g_0_xxxxx_0_xxxzzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxxzzz_0[i] * pb_y +
                                    g_0_xxxxxy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xxyyyy_0[i] = 4.0 * g_0_xxxyy_0_xxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyyyy_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyyy_0[i] * pb_x + g_0_xxxxyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxyyyz_0[i] = 4.0 * g_0_xxxyy_0_xxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyyz_0[i] * pb_x + g_0_xxxxyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxyyzz_0[i] = 4.0 * g_0_xxxyy_0_xxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyzz_0[i] * pb_x + g_0_xxxxyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxyzzz_0[i] = 4.0 * g_0_xxxyy_0_xxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xxyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyzzz_0[i] * pb_x + g_0_xxxxyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xxzzzz_0[i] = g_0_xxxxx_0_xxzzzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xxzzzz_0[i] * pb_y +
                                    g_0_xxxxxy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_xyyyyy_0[i] = 4.0 * g_0_xxxyy_0_xyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyyy_0[i] * pb_x + g_0_xxxxyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyyyyz_0[i] = 4.0 * g_0_xxxyy_0_xyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyyz_0[i] * pb_x + g_0_xxxxyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyyyzz_0[i] = 4.0 * g_0_xxxyy_0_xyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyzz_0[i] * pb_x + g_0_xxxxyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyyzzz_0[i] = 4.0 * g_0_xxxyy_0_xyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyzzz_0[i] * pb_x + g_0_xxxxyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xyzzzz_0[i] = 4.0 * g_0_xxxyy_0_xyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_xyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyzzzz_0[i] * pb_x + g_0_xxxxyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_xzzzzz_0[i] = g_0_xxxxx_0_xzzzzz_0[i] * fi_ab_0 - g_0_xxxxx_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxxxxy_0_xzzzzz_0[i] * pb_y +
                                    g_0_xxxxxy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyy_0_yyyyyy_0[i] = 4.0 * g_0_xxxyy_0_yyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_yyyyyy_0[i] * pb_x + g_0_xxxxyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyyyyz_0[i] = 4.0 * g_0_xxxyy_0_yyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_yyyyyz_0[i] * pb_x + g_0_xxxxyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyyyzz_0[i] = 4.0 * g_0_xxxyy_0_yyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_yyyyzz_0[i] * pb_x + g_0_xxxxyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyyzzz_0[i] = 4.0 * g_0_xxxyy_0_yyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_yyyzzz_0[i] * pb_x + g_0_xxxxyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yyzzzz_0[i] = 4.0 * g_0_xxxyy_0_yyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_yyzzzz_0[i] * pb_x + g_0_xxxxyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_yzzzzz_0[i] = 4.0 * g_0_xxxyy_0_yzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_yzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_yzzzzz_0[i] * pb_x + g_0_xxxxyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxxyy_0_zzzzzz_0[i] = 4.0 * g_0_xxxyy_0_zzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_zzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_zzzzzz_0[i] * pb_x + g_0_xxxxyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 112-140 components of targeted buffer : SKSI

    auto g_0_xxxxxyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 112);

    auto g_0_xxxxxyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 113);

    auto g_0_xxxxxyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 114);

    auto g_0_xxxxxyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 115);

    auto g_0_xxxxxyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 116);

    auto g_0_xxxxxyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 117);

    auto g_0_xxxxxyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 118);

    auto g_0_xxxxxyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 119);

    auto g_0_xxxxxyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 120);

    auto g_0_xxxxxyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 121);

    auto g_0_xxxxxyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 122);

    auto g_0_xxxxxyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 123);

    auto g_0_xxxxxyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 124);

    auto g_0_xxxxxyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 125);

    auto g_0_xxxxxyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 126);

    auto g_0_xxxxxyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 127);

    auto g_0_xxxxxyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 128);

    auto g_0_xxxxxyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 129);

    auto g_0_xxxxxyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 130);

    auto g_0_xxxxxyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 131);

    auto g_0_xxxxxyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 132);

    auto g_0_xxxxxyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 133);

    auto g_0_xxxxxyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 134);

    auto g_0_xxxxxyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 135);

    auto g_0_xxxxxyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 136);

    auto g_0_xxxxxyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 137);

    auto g_0_xxxxxyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 138);

    auto g_0_xxxxxyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 139);

#pragma omp simd aligned(g_0_xxxxxy_0_xxxxxy_0,      \
                             g_0_xxxxxy_0_xxxxxy_1,  \
                             g_0_xxxxxy_0_xxxxyy_0,  \
                             g_0_xxxxxy_0_xxxxyy_1,  \
                             g_0_xxxxxy_0_xxxyyy_0,  \
                             g_0_xxxxxy_0_xxxyyy_1,  \
                             g_0_xxxxxy_0_xxyyyy_0,  \
                             g_0_xxxxxy_0_xxyyyy_1,  \
                             g_0_xxxxxy_0_xyyyyy_0,  \
                             g_0_xxxxxy_0_xyyyyy_1,  \
                             g_0_xxxxxy_0_yyyyyy_0,  \
                             g_0_xxxxxy_0_yyyyyy_1,  \
                             g_0_xxxxxyz_0_xxxxxx_0, \
                             g_0_xxxxxyz_0_xxxxxy_0, \
                             g_0_xxxxxyz_0_xxxxxz_0, \
                             g_0_xxxxxyz_0_xxxxyy_0, \
                             g_0_xxxxxyz_0_xxxxyz_0, \
                             g_0_xxxxxyz_0_xxxxzz_0, \
                             g_0_xxxxxyz_0_xxxyyy_0, \
                             g_0_xxxxxyz_0_xxxyyz_0, \
                             g_0_xxxxxyz_0_xxxyzz_0, \
                             g_0_xxxxxyz_0_xxxzzz_0, \
                             g_0_xxxxxyz_0_xxyyyy_0, \
                             g_0_xxxxxyz_0_xxyyyz_0, \
                             g_0_xxxxxyz_0_xxyyzz_0, \
                             g_0_xxxxxyz_0_xxyzzz_0, \
                             g_0_xxxxxyz_0_xxzzzz_0, \
                             g_0_xxxxxyz_0_xyyyyy_0, \
                             g_0_xxxxxyz_0_xyyyyz_0, \
                             g_0_xxxxxyz_0_xyyyzz_0, \
                             g_0_xxxxxyz_0_xyyzzz_0, \
                             g_0_xxxxxyz_0_xyzzzz_0, \
                             g_0_xxxxxyz_0_xzzzzz_0, \
                             g_0_xxxxxyz_0_yyyyyy_0, \
                             g_0_xxxxxyz_0_yyyyyz_0, \
                             g_0_xxxxxyz_0_yyyyzz_0, \
                             g_0_xxxxxyz_0_yyyzzz_0, \
                             g_0_xxxxxyz_0_yyzzzz_0, \
                             g_0_xxxxxyz_0_yzzzzz_0, \
                             g_0_xxxxxyz_0_zzzzzz_0, \
                             g_0_xxxxxz_0_xxxxxx_0,  \
                             g_0_xxxxxz_0_xxxxxx_1,  \
                             g_0_xxxxxz_0_xxxxxz_0,  \
                             g_0_xxxxxz_0_xxxxxz_1,  \
                             g_0_xxxxxz_0_xxxxyz_0,  \
                             g_0_xxxxxz_0_xxxxyz_1,  \
                             g_0_xxxxxz_0_xxxxz_1,   \
                             g_0_xxxxxz_0_xxxxzz_0,  \
                             g_0_xxxxxz_0_xxxxzz_1,  \
                             g_0_xxxxxz_0_xxxyyz_0,  \
                             g_0_xxxxxz_0_xxxyyz_1,  \
                             g_0_xxxxxz_0_xxxyz_1,   \
                             g_0_xxxxxz_0_xxxyzz_0,  \
                             g_0_xxxxxz_0_xxxyzz_1,  \
                             g_0_xxxxxz_0_xxxzz_1,   \
                             g_0_xxxxxz_0_xxxzzz_0,  \
                             g_0_xxxxxz_0_xxxzzz_1,  \
                             g_0_xxxxxz_0_xxyyyz_0,  \
                             g_0_xxxxxz_0_xxyyyz_1,  \
                             g_0_xxxxxz_0_xxyyz_1,   \
                             g_0_xxxxxz_0_xxyyzz_0,  \
                             g_0_xxxxxz_0_xxyyzz_1,  \
                             g_0_xxxxxz_0_xxyzz_1,   \
                             g_0_xxxxxz_0_xxyzzz_0,  \
                             g_0_xxxxxz_0_xxyzzz_1,  \
                             g_0_xxxxxz_0_xxzzz_1,   \
                             g_0_xxxxxz_0_xxzzzz_0,  \
                             g_0_xxxxxz_0_xxzzzz_1,  \
                             g_0_xxxxxz_0_xyyyyz_0,  \
                             g_0_xxxxxz_0_xyyyyz_1,  \
                             g_0_xxxxxz_0_xyyyz_1,   \
                             g_0_xxxxxz_0_xyyyzz_0,  \
                             g_0_xxxxxz_0_xyyyzz_1,  \
                             g_0_xxxxxz_0_xyyzz_1,   \
                             g_0_xxxxxz_0_xyyzzz_0,  \
                             g_0_xxxxxz_0_xyyzzz_1,  \
                             g_0_xxxxxz_0_xyzzz_1,   \
                             g_0_xxxxxz_0_xyzzzz_0,  \
                             g_0_xxxxxz_0_xyzzzz_1,  \
                             g_0_xxxxxz_0_xzzzz_1,   \
                             g_0_xxxxxz_0_xzzzzz_0,  \
                             g_0_xxxxxz_0_xzzzzz_1,  \
                             g_0_xxxxxz_0_yyyyyz_0,  \
                             g_0_xxxxxz_0_yyyyyz_1,  \
                             g_0_xxxxxz_0_yyyyz_1,   \
                             g_0_xxxxxz_0_yyyyzz_0,  \
                             g_0_xxxxxz_0_yyyyzz_1,  \
                             g_0_xxxxxz_0_yyyzz_1,   \
                             g_0_xxxxxz_0_yyyzzz_0,  \
                             g_0_xxxxxz_0_yyyzzz_1,  \
                             g_0_xxxxxz_0_yyzzz_1,   \
                             g_0_xxxxxz_0_yyzzzz_0,  \
                             g_0_xxxxxz_0_yyzzzz_1,  \
                             g_0_xxxxxz_0_yzzzz_1,   \
                             g_0_xxxxxz_0_yzzzzz_0,  \
                             g_0_xxxxxz_0_yzzzzz_1,  \
                             g_0_xxxxxz_0_zzzzz_1,   \
                             g_0_xxxxxz_0_zzzzzz_0,  \
                             g_0_xxxxxz_0_zzzzzz_1,  \
                             wp_y,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyz_0_xxxxxx_0[i] = g_0_xxxxxz_0_xxxxxx_0[i] * pb_y + g_0_xxxxxz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxxxy_0[i] = g_0_xxxxxy_0_xxxxxy_0[i] * pb_z + g_0_xxxxxy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxxxxz_0[i] = g_0_xxxxxz_0_xxxxxz_0[i] * pb_y + g_0_xxxxxz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxxyy_0[i] = g_0_xxxxxy_0_xxxxyy_0[i] * pb_z + g_0_xxxxxy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxxxyz_0[i] = g_0_xxxxxz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxxxyz_0[i] * pb_y + g_0_xxxxxz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxxzz_0[i] = g_0_xxxxxz_0_xxxxzz_0[i] * pb_y + g_0_xxxxxz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxyyy_0[i] = g_0_xxxxxy_0_xxxyyy_0[i] * pb_z + g_0_xxxxxy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxxyyz_0[i] = 2.0 * g_0_xxxxxz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxxyyz_0[i] * pb_y + g_0_xxxxxz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxyzz_0[i] = g_0_xxxxxz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxxyzz_0[i] * pb_y + g_0_xxxxxz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxxzzz_0[i] = g_0_xxxxxz_0_xxxzzz_0[i] * pb_y + g_0_xxxxxz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxyyyy_0[i] = g_0_xxxxxy_0_xxyyyy_0[i] * pb_z + g_0_xxxxxy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xxyyyz_0[i] = 3.0 * g_0_xxxxxz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxyyyz_0[i] * pb_y + g_0_xxxxxz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxyyzz_0[i] = 2.0 * g_0_xxxxxz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxyyzz_0[i] * pb_y + g_0_xxxxxz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxyzzz_0[i] = g_0_xxxxxz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xxyzzz_0[i] * pb_y + g_0_xxxxxz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xxzzzz_0[i] = g_0_xxxxxz_0_xxzzzz_0[i] * pb_y + g_0_xxxxxz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyyyyy_0[i] = g_0_xxxxxy_0_xyyyyy_0[i] * pb_z + g_0_xxxxxy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_xyyyyz_0[i] = 4.0 * g_0_xxxxxz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyyyyz_0[i] * pb_y + g_0_xxxxxz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyyyzz_0[i] = 3.0 * g_0_xxxxxz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyyyzz_0[i] * pb_y + g_0_xxxxxz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyyzzz_0[i] = 2.0 * g_0_xxxxxz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyyzzz_0[i] * pb_y + g_0_xxxxxz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xyzzzz_0[i] = g_0_xxxxxz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_xyzzzz_0[i] * pb_y + g_0_xxxxxz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_xzzzzz_0[i] = g_0_xxxxxz_0_xzzzzz_0[i] * pb_y + g_0_xxxxxz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyyyyy_0[i] = g_0_xxxxxy_0_yyyyyy_0[i] * pb_z + g_0_xxxxxy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxxxyz_0_yyyyyz_0[i] = 5.0 * g_0_xxxxxz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyyyyz_0[i] * pb_y + g_0_xxxxxz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyyyzz_0[i] = 4.0 * g_0_xxxxxz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyyyzz_0[i] * pb_y + g_0_xxxxxz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyyzzz_0[i] = 3.0 * g_0_xxxxxz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyyzzz_0[i] * pb_y + g_0_xxxxxz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yyzzzz_0[i] = 2.0 * g_0_xxxxxz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yyzzzz_0[i] * pb_y + g_0_xxxxxz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_yzzzzz_0[i] = g_0_xxxxxz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxxz_0_yzzzzz_0[i] * pb_y + g_0_xxxxxz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyz_0_zzzzzz_0[i] = g_0_xxxxxz_0_zzzzzz_0[i] * pb_y + g_0_xxxxxz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 140-168 components of targeted buffer : SKSI

    auto g_0_xxxxxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 140);

    auto g_0_xxxxxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 141);

    auto g_0_xxxxxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 142);

    auto g_0_xxxxxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 143);

    auto g_0_xxxxxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 144);

    auto g_0_xxxxxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 145);

    auto g_0_xxxxxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 146);

    auto g_0_xxxxxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 147);

    auto g_0_xxxxxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 148);

    auto g_0_xxxxxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 149);

    auto g_0_xxxxxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 150);

    auto g_0_xxxxxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 151);

    auto g_0_xxxxxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 152);

    auto g_0_xxxxxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 153);

    auto g_0_xxxxxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 154);

    auto g_0_xxxxxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 155);

    auto g_0_xxxxxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 156);

    auto g_0_xxxxxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 157);

    auto g_0_xxxxxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 158);

    auto g_0_xxxxxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 159);

    auto g_0_xxxxxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 160);

    auto g_0_xxxxxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 161);

    auto g_0_xxxxxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 162);

    auto g_0_xxxxxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 163);

    auto g_0_xxxxxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 164);

    auto g_0_xxxxxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 165);

    auto g_0_xxxxxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 166);

    auto g_0_xxxxxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 167);

#pragma omp simd aligned(g_0_xxxxx_0_xxxxxx_0,       \
                             g_0_xxxxx_0_xxxxxx_1,   \
                             g_0_xxxxx_0_xxxxxy_0,   \
                             g_0_xxxxx_0_xxxxxy_1,   \
                             g_0_xxxxx_0_xxxxyy_0,   \
                             g_0_xxxxx_0_xxxxyy_1,   \
                             g_0_xxxxx_0_xxxyyy_0,   \
                             g_0_xxxxx_0_xxxyyy_1,   \
                             g_0_xxxxx_0_xxyyyy_0,   \
                             g_0_xxxxx_0_xxyyyy_1,   \
                             g_0_xxxxx_0_xyyyyy_0,   \
                             g_0_xxxxx_0_xyyyyy_1,   \
                             g_0_xxxxxz_0_xxxxxx_0,  \
                             g_0_xxxxxz_0_xxxxxx_1,  \
                             g_0_xxxxxz_0_xxxxxy_0,  \
                             g_0_xxxxxz_0_xxxxxy_1,  \
                             g_0_xxxxxz_0_xxxxyy_0,  \
                             g_0_xxxxxz_0_xxxxyy_1,  \
                             g_0_xxxxxz_0_xxxyyy_0,  \
                             g_0_xxxxxz_0_xxxyyy_1,  \
                             g_0_xxxxxz_0_xxyyyy_0,  \
                             g_0_xxxxxz_0_xxyyyy_1,  \
                             g_0_xxxxxz_0_xyyyyy_0,  \
                             g_0_xxxxxz_0_xyyyyy_1,  \
                             g_0_xxxxxzz_0_xxxxxx_0, \
                             g_0_xxxxxzz_0_xxxxxy_0, \
                             g_0_xxxxxzz_0_xxxxxz_0, \
                             g_0_xxxxxzz_0_xxxxyy_0, \
                             g_0_xxxxxzz_0_xxxxyz_0, \
                             g_0_xxxxxzz_0_xxxxzz_0, \
                             g_0_xxxxxzz_0_xxxyyy_0, \
                             g_0_xxxxxzz_0_xxxyyz_0, \
                             g_0_xxxxxzz_0_xxxyzz_0, \
                             g_0_xxxxxzz_0_xxxzzz_0, \
                             g_0_xxxxxzz_0_xxyyyy_0, \
                             g_0_xxxxxzz_0_xxyyyz_0, \
                             g_0_xxxxxzz_0_xxyyzz_0, \
                             g_0_xxxxxzz_0_xxyzzz_0, \
                             g_0_xxxxxzz_0_xxzzzz_0, \
                             g_0_xxxxxzz_0_xyyyyy_0, \
                             g_0_xxxxxzz_0_xyyyyz_0, \
                             g_0_xxxxxzz_0_xyyyzz_0, \
                             g_0_xxxxxzz_0_xyyzzz_0, \
                             g_0_xxxxxzz_0_xyzzzz_0, \
                             g_0_xxxxxzz_0_xzzzzz_0, \
                             g_0_xxxxxzz_0_yyyyyy_0, \
                             g_0_xxxxxzz_0_yyyyyz_0, \
                             g_0_xxxxxzz_0_yyyyzz_0, \
                             g_0_xxxxxzz_0_yyyzzz_0, \
                             g_0_xxxxxzz_0_yyzzzz_0, \
                             g_0_xxxxxzz_0_yzzzzz_0, \
                             g_0_xxxxxzz_0_zzzzzz_0, \
                             g_0_xxxxzz_0_xxxxxz_0,  \
                             g_0_xxxxzz_0_xxxxxz_1,  \
                             g_0_xxxxzz_0_xxxxyz_0,  \
                             g_0_xxxxzz_0_xxxxyz_1,  \
                             g_0_xxxxzz_0_xxxxz_1,   \
                             g_0_xxxxzz_0_xxxxzz_0,  \
                             g_0_xxxxzz_0_xxxxzz_1,  \
                             g_0_xxxxzz_0_xxxyyz_0,  \
                             g_0_xxxxzz_0_xxxyyz_1,  \
                             g_0_xxxxzz_0_xxxyz_1,   \
                             g_0_xxxxzz_0_xxxyzz_0,  \
                             g_0_xxxxzz_0_xxxyzz_1,  \
                             g_0_xxxxzz_0_xxxzz_1,   \
                             g_0_xxxxzz_0_xxxzzz_0,  \
                             g_0_xxxxzz_0_xxxzzz_1,  \
                             g_0_xxxxzz_0_xxyyyz_0,  \
                             g_0_xxxxzz_0_xxyyyz_1,  \
                             g_0_xxxxzz_0_xxyyz_1,   \
                             g_0_xxxxzz_0_xxyyzz_0,  \
                             g_0_xxxxzz_0_xxyyzz_1,  \
                             g_0_xxxxzz_0_xxyzz_1,   \
                             g_0_xxxxzz_0_xxyzzz_0,  \
                             g_0_xxxxzz_0_xxyzzz_1,  \
                             g_0_xxxxzz_0_xxzzz_1,   \
                             g_0_xxxxzz_0_xxzzzz_0,  \
                             g_0_xxxxzz_0_xxzzzz_1,  \
                             g_0_xxxxzz_0_xyyyyz_0,  \
                             g_0_xxxxzz_0_xyyyyz_1,  \
                             g_0_xxxxzz_0_xyyyz_1,   \
                             g_0_xxxxzz_0_xyyyzz_0,  \
                             g_0_xxxxzz_0_xyyyzz_1,  \
                             g_0_xxxxzz_0_xyyzz_1,   \
                             g_0_xxxxzz_0_xyyzzz_0,  \
                             g_0_xxxxzz_0_xyyzzz_1,  \
                             g_0_xxxxzz_0_xyzzz_1,   \
                             g_0_xxxxzz_0_xyzzzz_0,  \
                             g_0_xxxxzz_0_xyzzzz_1,  \
                             g_0_xxxxzz_0_xzzzz_1,   \
                             g_0_xxxxzz_0_xzzzzz_0,  \
                             g_0_xxxxzz_0_xzzzzz_1,  \
                             g_0_xxxxzz_0_yyyyyy_0,  \
                             g_0_xxxxzz_0_yyyyyy_1,  \
                             g_0_xxxxzz_0_yyyyyz_0,  \
                             g_0_xxxxzz_0_yyyyyz_1,  \
                             g_0_xxxxzz_0_yyyyz_1,   \
                             g_0_xxxxzz_0_yyyyzz_0,  \
                             g_0_xxxxzz_0_yyyyzz_1,  \
                             g_0_xxxxzz_0_yyyzz_1,   \
                             g_0_xxxxzz_0_yyyzzz_0,  \
                             g_0_xxxxzz_0_yyyzzz_1,  \
                             g_0_xxxxzz_0_yyzzz_1,   \
                             g_0_xxxxzz_0_yyzzzz_0,  \
                             g_0_xxxxzz_0_yyzzzz_1,  \
                             g_0_xxxxzz_0_yzzzz_1,   \
                             g_0_xxxxzz_0_yzzzzz_0,  \
                             g_0_xxxxzz_0_yzzzzz_1,  \
                             g_0_xxxxzz_0_zzzzz_1,   \
                             g_0_xxxxzz_0_zzzzzz_0,  \
                             g_0_xxxxzz_0_zzzzzz_1,  \
                             g_0_xxxzz_0_xxxxxz_0,   \
                             g_0_xxxzz_0_xxxxxz_1,   \
                             g_0_xxxzz_0_xxxxyz_0,   \
                             g_0_xxxzz_0_xxxxyz_1,   \
                             g_0_xxxzz_0_xxxxzz_0,   \
                             g_0_xxxzz_0_xxxxzz_1,   \
                             g_0_xxxzz_0_xxxyyz_0,   \
                             g_0_xxxzz_0_xxxyyz_1,   \
                             g_0_xxxzz_0_xxxyzz_0,   \
                             g_0_xxxzz_0_xxxyzz_1,   \
                             g_0_xxxzz_0_xxxzzz_0,   \
                             g_0_xxxzz_0_xxxzzz_1,   \
                             g_0_xxxzz_0_xxyyyz_0,   \
                             g_0_xxxzz_0_xxyyyz_1,   \
                             g_0_xxxzz_0_xxyyzz_0,   \
                             g_0_xxxzz_0_xxyyzz_1,   \
                             g_0_xxxzz_0_xxyzzz_0,   \
                             g_0_xxxzz_0_xxyzzz_1,   \
                             g_0_xxxzz_0_xxzzzz_0,   \
                             g_0_xxxzz_0_xxzzzz_1,   \
                             g_0_xxxzz_0_xyyyyz_0,   \
                             g_0_xxxzz_0_xyyyyz_1,   \
                             g_0_xxxzz_0_xyyyzz_0,   \
                             g_0_xxxzz_0_xyyyzz_1,   \
                             g_0_xxxzz_0_xyyzzz_0,   \
                             g_0_xxxzz_0_xyyzzz_1,   \
                             g_0_xxxzz_0_xyzzzz_0,   \
                             g_0_xxxzz_0_xyzzzz_1,   \
                             g_0_xxxzz_0_xzzzzz_0,   \
                             g_0_xxxzz_0_xzzzzz_1,   \
                             g_0_xxxzz_0_yyyyyy_0,   \
                             g_0_xxxzz_0_yyyyyy_1,   \
                             g_0_xxxzz_0_yyyyyz_0,   \
                             g_0_xxxzz_0_yyyyyz_1,   \
                             g_0_xxxzz_0_yyyyzz_0,   \
                             g_0_xxxzz_0_yyyyzz_1,   \
                             g_0_xxxzz_0_yyyzzz_0,   \
                             g_0_xxxzz_0_yyyzzz_1,   \
                             g_0_xxxzz_0_yyzzzz_0,   \
                             g_0_xxxzz_0_yyzzzz_1,   \
                             g_0_xxxzz_0_yzzzzz_0,   \
                             g_0_xxxzz_0_yzzzzz_1,   \
                             g_0_xxxzz_0_zzzzzz_0,   \
                             g_0_xxxzz_0_zzzzzz_1,   \
                             wp_x,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxzz_0_xxxxxx_0[i] = g_0_xxxxx_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxxxx_0[i] * pb_z +
                                    g_0_xxxxxz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxxxy_0[i] = g_0_xxxxx_0_xxxxxy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxxxy_0[i] * pb_z +
                                    g_0_xxxxxz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxxxz_0[i] = 4.0 * g_0_xxxzz_0_xxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxxxz_1[i] * fti_ab_0 +
                                    5.0 * g_0_xxxxzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxxz_0[i] * pb_x + g_0_xxxxzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxxyy_0[i] = g_0_xxxxx_0_xxxxyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxxyy_0[i] * pb_z +
                                    g_0_xxxxxz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxxyz_0[i] = 4.0 * g_0_xxxzz_0_xxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxxyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxxxzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxyz_0[i] * pb_x + g_0_xxxxzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxxzz_0[i] = 4.0 * g_0_xxxzz_0_xxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxxzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxxxzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxzz_0[i] * pb_x + g_0_xxxxzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxyyy_0[i] = g_0_xxxxx_0_xxxyyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxxyyy_0[i] * pb_z +
                                    g_0_xxxxxz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxxyyz_0[i] = 4.0 * g_0_xxxzz_0_xxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxxzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyyz_0[i] * pb_x + g_0_xxxxzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxyzz_0[i] = 4.0 * g_0_xxxzz_0_xxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxxzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyzz_0[i] * pb_x + g_0_xxxxzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxxzzz_0[i] = 4.0 * g_0_xxxzz_0_xxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxxzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxxzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxzzz_0[i] * pb_x + g_0_xxxxzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxyyyy_0[i] = g_0_xxxxx_0_xxyyyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xxyyyy_0[i] * pb_z +
                                    g_0_xxxxxz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xxyyyz_0[i] = 4.0 * g_0_xxxzz_0_xxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxyyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyyz_0[i] * pb_x + g_0_xxxxzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxyyzz_0[i] = 4.0 * g_0_xxxzz_0_xxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyzz_0[i] * pb_x + g_0_xxxxzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxyzzz_0[i] = 4.0 * g_0_xxxzz_0_xxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyzzz_0[i] * pb_x + g_0_xxxxzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xxzzzz_0[i] = 4.0 * g_0_xxxzz_0_xxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xxzzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxxzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxzzzz_0[i] * pb_x + g_0_xxxxzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyyyyy_0[i] = g_0_xxxxx_0_xyyyyy_0[i] * fi_ab_0 - g_0_xxxxx_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxxxxz_0_xyyyyy_0[i] * pb_z +
                                    g_0_xxxxxz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxxzz_0_xyyyyz_0[i] = 4.0 * g_0_xxxzz_0_xyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyyz_0[i] * pb_x + g_0_xxxxzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyyyzz_0[i] = 4.0 * g_0_xxxzz_0_xyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyzz_0[i] * pb_x + g_0_xxxxzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyyzzz_0[i] = 4.0 * g_0_xxxzz_0_xyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyzzz_0[i] * pb_x + g_0_xxxxzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xyzzzz_0[i] = 4.0 * g_0_xxxzz_0_xyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyzzzz_0[i] * pb_x + g_0_xxxxzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_xzzzzz_0[i] = 4.0 * g_0_xxxzz_0_xzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_xzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xzzzzz_0[i] * pb_x + g_0_xxxxzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyyyy_0[i] = 4.0 * g_0_xxxzz_0_yyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_yyyyyy_0[i] * pb_x + g_0_xxxxzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyyyz_0[i] = 4.0 * g_0_xxxzz_0_yyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_yyyyyz_0[i] * pb_x + g_0_xxxxzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyyzz_0[i] = 4.0 * g_0_xxxzz_0_yyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_yyyyzz_0[i] * pb_x + g_0_xxxxzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyyzzz_0[i] = 4.0 * g_0_xxxzz_0_yyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_yyyzzz_0[i] * pb_x + g_0_xxxxzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yyzzzz_0[i] = 4.0 * g_0_xxxzz_0_yyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_yyzzzz_0[i] * pb_x + g_0_xxxxzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_yzzzzz_0[i] = 4.0 * g_0_xxxzz_0_yzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_yzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_yzzzzz_0[i] * pb_x + g_0_xxxxzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzz_0_zzzzzz_0[i] = 4.0 * g_0_xxxzz_0_zzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_zzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_zzzzzz_0[i] * pb_x + g_0_xxxxzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 168-196 components of targeted buffer : SKSI

    auto g_0_xxxxyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 168);

    auto g_0_xxxxyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 169);

    auto g_0_xxxxyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 170);

    auto g_0_xxxxyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 171);

    auto g_0_xxxxyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 172);

    auto g_0_xxxxyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 173);

    auto g_0_xxxxyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 174);

    auto g_0_xxxxyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 175);

    auto g_0_xxxxyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 176);

    auto g_0_xxxxyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 177);

    auto g_0_xxxxyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 178);

    auto g_0_xxxxyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 179);

    auto g_0_xxxxyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 180);

    auto g_0_xxxxyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 181);

    auto g_0_xxxxyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 182);

    auto g_0_xxxxyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 183);

    auto g_0_xxxxyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 184);

    auto g_0_xxxxyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 185);

    auto g_0_xxxxyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 186);

    auto g_0_xxxxyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 187);

    auto g_0_xxxxyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 188);

    auto g_0_xxxxyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 189);

    auto g_0_xxxxyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 190);

    auto g_0_xxxxyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 191);

    auto g_0_xxxxyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 192);

    auto g_0_xxxxyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 193);

    auto g_0_xxxxyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 194);

    auto g_0_xxxxyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 195);

#pragma omp simd aligned(g_0_xxxxy_0_xxxxxx_0,       \
                             g_0_xxxxy_0_xxxxxx_1,   \
                             g_0_xxxxy_0_xxxxxz_0,   \
                             g_0_xxxxy_0_xxxxxz_1,   \
                             g_0_xxxxy_0_xxxxzz_0,   \
                             g_0_xxxxy_0_xxxxzz_1,   \
                             g_0_xxxxy_0_xxxzzz_0,   \
                             g_0_xxxxy_0_xxxzzz_1,   \
                             g_0_xxxxy_0_xxzzzz_0,   \
                             g_0_xxxxy_0_xxzzzz_1,   \
                             g_0_xxxxy_0_xzzzzz_0,   \
                             g_0_xxxxy_0_xzzzzz_1,   \
                             g_0_xxxxyy_0_xxxxxx_0,  \
                             g_0_xxxxyy_0_xxxxxx_1,  \
                             g_0_xxxxyy_0_xxxxxz_0,  \
                             g_0_xxxxyy_0_xxxxxz_1,  \
                             g_0_xxxxyy_0_xxxxzz_0,  \
                             g_0_xxxxyy_0_xxxxzz_1,  \
                             g_0_xxxxyy_0_xxxzzz_0,  \
                             g_0_xxxxyy_0_xxxzzz_1,  \
                             g_0_xxxxyy_0_xxzzzz_0,  \
                             g_0_xxxxyy_0_xxzzzz_1,  \
                             g_0_xxxxyy_0_xzzzzz_0,  \
                             g_0_xxxxyy_0_xzzzzz_1,  \
                             g_0_xxxxyyy_0_xxxxxx_0, \
                             g_0_xxxxyyy_0_xxxxxy_0, \
                             g_0_xxxxyyy_0_xxxxxz_0, \
                             g_0_xxxxyyy_0_xxxxyy_0, \
                             g_0_xxxxyyy_0_xxxxyz_0, \
                             g_0_xxxxyyy_0_xxxxzz_0, \
                             g_0_xxxxyyy_0_xxxyyy_0, \
                             g_0_xxxxyyy_0_xxxyyz_0, \
                             g_0_xxxxyyy_0_xxxyzz_0, \
                             g_0_xxxxyyy_0_xxxzzz_0, \
                             g_0_xxxxyyy_0_xxyyyy_0, \
                             g_0_xxxxyyy_0_xxyyyz_0, \
                             g_0_xxxxyyy_0_xxyyzz_0, \
                             g_0_xxxxyyy_0_xxyzzz_0, \
                             g_0_xxxxyyy_0_xxzzzz_0, \
                             g_0_xxxxyyy_0_xyyyyy_0, \
                             g_0_xxxxyyy_0_xyyyyz_0, \
                             g_0_xxxxyyy_0_xyyyzz_0, \
                             g_0_xxxxyyy_0_xyyzzz_0, \
                             g_0_xxxxyyy_0_xyzzzz_0, \
                             g_0_xxxxyyy_0_xzzzzz_0, \
                             g_0_xxxxyyy_0_yyyyyy_0, \
                             g_0_xxxxyyy_0_yyyyyz_0, \
                             g_0_xxxxyyy_0_yyyyzz_0, \
                             g_0_xxxxyyy_0_yyyzzz_0, \
                             g_0_xxxxyyy_0_yyzzzz_0, \
                             g_0_xxxxyyy_0_yzzzzz_0, \
                             g_0_xxxxyyy_0_zzzzzz_0, \
                             g_0_xxxyyy_0_xxxxxy_0,  \
                             g_0_xxxyyy_0_xxxxxy_1,  \
                             g_0_xxxyyy_0_xxxxy_1,   \
                             g_0_xxxyyy_0_xxxxyy_0,  \
                             g_0_xxxyyy_0_xxxxyy_1,  \
                             g_0_xxxyyy_0_xxxxyz_0,  \
                             g_0_xxxyyy_0_xxxxyz_1,  \
                             g_0_xxxyyy_0_xxxyy_1,   \
                             g_0_xxxyyy_0_xxxyyy_0,  \
                             g_0_xxxyyy_0_xxxyyy_1,  \
                             g_0_xxxyyy_0_xxxyyz_0,  \
                             g_0_xxxyyy_0_xxxyyz_1,  \
                             g_0_xxxyyy_0_xxxyz_1,   \
                             g_0_xxxyyy_0_xxxyzz_0,  \
                             g_0_xxxyyy_0_xxxyzz_1,  \
                             g_0_xxxyyy_0_xxyyy_1,   \
                             g_0_xxxyyy_0_xxyyyy_0,  \
                             g_0_xxxyyy_0_xxyyyy_1,  \
                             g_0_xxxyyy_0_xxyyyz_0,  \
                             g_0_xxxyyy_0_xxyyyz_1,  \
                             g_0_xxxyyy_0_xxyyz_1,   \
                             g_0_xxxyyy_0_xxyyzz_0,  \
                             g_0_xxxyyy_0_xxyyzz_1,  \
                             g_0_xxxyyy_0_xxyzz_1,   \
                             g_0_xxxyyy_0_xxyzzz_0,  \
                             g_0_xxxyyy_0_xxyzzz_1,  \
                             g_0_xxxyyy_0_xyyyy_1,   \
                             g_0_xxxyyy_0_xyyyyy_0,  \
                             g_0_xxxyyy_0_xyyyyy_1,  \
                             g_0_xxxyyy_0_xyyyyz_0,  \
                             g_0_xxxyyy_0_xyyyyz_1,  \
                             g_0_xxxyyy_0_xyyyz_1,   \
                             g_0_xxxyyy_0_xyyyzz_0,  \
                             g_0_xxxyyy_0_xyyyzz_1,  \
                             g_0_xxxyyy_0_xyyzz_1,   \
                             g_0_xxxyyy_0_xyyzzz_0,  \
                             g_0_xxxyyy_0_xyyzzz_1,  \
                             g_0_xxxyyy_0_xyzzz_1,   \
                             g_0_xxxyyy_0_xyzzzz_0,  \
                             g_0_xxxyyy_0_xyzzzz_1,  \
                             g_0_xxxyyy_0_yyyyy_1,   \
                             g_0_xxxyyy_0_yyyyyy_0,  \
                             g_0_xxxyyy_0_yyyyyy_1,  \
                             g_0_xxxyyy_0_yyyyyz_0,  \
                             g_0_xxxyyy_0_yyyyyz_1,  \
                             g_0_xxxyyy_0_yyyyz_1,   \
                             g_0_xxxyyy_0_yyyyzz_0,  \
                             g_0_xxxyyy_0_yyyyzz_1,  \
                             g_0_xxxyyy_0_yyyzz_1,   \
                             g_0_xxxyyy_0_yyyzzz_0,  \
                             g_0_xxxyyy_0_yyyzzz_1,  \
                             g_0_xxxyyy_0_yyzzz_1,   \
                             g_0_xxxyyy_0_yyzzzz_0,  \
                             g_0_xxxyyy_0_yyzzzz_1,  \
                             g_0_xxxyyy_0_yzzzz_1,   \
                             g_0_xxxyyy_0_yzzzzz_0,  \
                             g_0_xxxyyy_0_yzzzzz_1,  \
                             g_0_xxxyyy_0_zzzzzz_0,  \
                             g_0_xxxyyy_0_zzzzzz_1,  \
                             g_0_xxyyy_0_xxxxxy_0,   \
                             g_0_xxyyy_0_xxxxxy_1,   \
                             g_0_xxyyy_0_xxxxyy_0,   \
                             g_0_xxyyy_0_xxxxyy_1,   \
                             g_0_xxyyy_0_xxxxyz_0,   \
                             g_0_xxyyy_0_xxxxyz_1,   \
                             g_0_xxyyy_0_xxxyyy_0,   \
                             g_0_xxyyy_0_xxxyyy_1,   \
                             g_0_xxyyy_0_xxxyyz_0,   \
                             g_0_xxyyy_0_xxxyyz_1,   \
                             g_0_xxyyy_0_xxxyzz_0,   \
                             g_0_xxyyy_0_xxxyzz_1,   \
                             g_0_xxyyy_0_xxyyyy_0,   \
                             g_0_xxyyy_0_xxyyyy_1,   \
                             g_0_xxyyy_0_xxyyyz_0,   \
                             g_0_xxyyy_0_xxyyyz_1,   \
                             g_0_xxyyy_0_xxyyzz_0,   \
                             g_0_xxyyy_0_xxyyzz_1,   \
                             g_0_xxyyy_0_xxyzzz_0,   \
                             g_0_xxyyy_0_xxyzzz_1,   \
                             g_0_xxyyy_0_xyyyyy_0,   \
                             g_0_xxyyy_0_xyyyyy_1,   \
                             g_0_xxyyy_0_xyyyyz_0,   \
                             g_0_xxyyy_0_xyyyyz_1,   \
                             g_0_xxyyy_0_xyyyzz_0,   \
                             g_0_xxyyy_0_xyyyzz_1,   \
                             g_0_xxyyy_0_xyyzzz_0,   \
                             g_0_xxyyy_0_xyyzzz_1,   \
                             g_0_xxyyy_0_xyzzzz_0,   \
                             g_0_xxyyy_0_xyzzzz_1,   \
                             g_0_xxyyy_0_yyyyyy_0,   \
                             g_0_xxyyy_0_yyyyyy_1,   \
                             g_0_xxyyy_0_yyyyyz_0,   \
                             g_0_xxyyy_0_yyyyyz_1,   \
                             g_0_xxyyy_0_yyyyzz_0,   \
                             g_0_xxyyy_0_yyyyzz_1,   \
                             g_0_xxyyy_0_yyyzzz_0,   \
                             g_0_xxyyy_0_yyyzzz_1,   \
                             g_0_xxyyy_0_yyzzzz_0,   \
                             g_0_xxyyy_0_yyzzzz_1,   \
                             g_0_xxyyy_0_yzzzzz_0,   \
                             g_0_xxyyy_0_yzzzzz_1,   \
                             g_0_xxyyy_0_zzzzzz_0,   \
                             g_0_xxyyy_0_zzzzzz_1,   \
                             wp_x,                   \
                             wp_y,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyy_0_xxxxxx_0[i] = 2.0 * g_0_xxxxy_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxxxx_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_xxxxxx_0[i] * pb_y + g_0_xxxxyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxxxxy_0[i] = 3.0 * g_0_xxyyy_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxxxy_1[i] * fti_ab_0 +
                                    5.0 * g_0_xxxyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxxy_0[i] * pb_x + g_0_xxxyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxxxz_0[i] = 2.0 * g_0_xxxxy_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxxxz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_xxxxxz_0[i] * pb_y + g_0_xxxxyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxxxyy_0[i] = 3.0 * g_0_xxyyy_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxxyy_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxxyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxyy_0[i] * pb_x + g_0_xxxyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxxyz_0[i] = 3.0 * g_0_xxyyy_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxxyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxxyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxyz_0[i] * pb_x + g_0_xxxyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxxzz_0[i] = 2.0 * g_0_xxxxy_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxxzz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_xxxxzz_0[i] * pb_y + g_0_xxxxyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxxyyy_0[i] = 3.0 * g_0_xxyyy_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxyyy_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyyy_0[i] * pb_x + g_0_xxxyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxyyz_0[i] = 3.0 * g_0_xxyyy_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyyz_0[i] * pb_x + g_0_xxxyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxyzz_0[i] = 3.0 * g_0_xxyyy_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxxyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyzz_0[i] * pb_x + g_0_xxxyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxxzzz_0[i] = 2.0 * g_0_xxxxy_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxxzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_xxxzzz_0[i] * pb_y + g_0_xxxxyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xxyyyy_0[i] = 3.0 * g_0_xxyyy_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyyyy_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyyy_0[i] * pb_x + g_0_xxxyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxyyyz_0[i] = 3.0 * g_0_xxyyy_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyyz_0[i] * pb_x + g_0_xxxyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxyyzz_0[i] = 3.0 * g_0_xxyyy_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyzz_0[i] * pb_x + g_0_xxxyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxyzzz_0[i] = 3.0 * g_0_xxyyy_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xxyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyzzz_0[i] * pb_x + g_0_xxxyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xxzzzz_0[i] = 2.0 * g_0_xxxxy_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xxzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_xxzzzz_0[i] * pb_y + g_0_xxxxyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_xyyyyy_0[i] = 3.0 * g_0_xxyyy_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyyy_0[i] * pb_x + g_0_xxxyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyyyyz_0[i] = 3.0 * g_0_xxyyy_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyyz_0[i] * pb_x + g_0_xxxyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyyyzz_0[i] = 3.0 * g_0_xxyyy_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyzz_0[i] * pb_x + g_0_xxxyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyyzzz_0[i] = 3.0 * g_0_xxyyy_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyzzz_0[i] * pb_x + g_0_xxxyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xyzzzz_0[i] = 3.0 * g_0_xxyyy_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_xyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyzzzz_0[i] * pb_x + g_0_xxxyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_xzzzzz_0[i] = 2.0 * g_0_xxxxy_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxy_0_xzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxxyy_0_xzzzzz_0[i] * pb_y + g_0_xxxxyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxyyy_0_yyyyyy_0[i] = 3.0 * g_0_xxyyy_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_yyyyyy_0[i] * pb_x + g_0_xxxyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyyyyz_0[i] = 3.0 * g_0_xxyyy_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_yyyyyz_0[i] * pb_x + g_0_xxxyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyyyzz_0[i] = 3.0 * g_0_xxyyy_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_yyyyzz_0[i] * pb_x + g_0_xxxyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyyzzz_0[i] = 3.0 * g_0_xxyyy_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_yyyzzz_0[i] * pb_x + g_0_xxxyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yyzzzz_0[i] = 3.0 * g_0_xxyyy_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_yyzzzz_0[i] * pb_x + g_0_xxxyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_yzzzzz_0[i] = 3.0 * g_0_xxyyy_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_yzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_yzzzzz_0[i] * pb_x + g_0_xxxyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxyyy_0_zzzzzz_0[i] = 3.0 * g_0_xxyyy_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_zzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_zzzzzz_0[i] * pb_x + g_0_xxxyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 196-224 components of targeted buffer : SKSI

    auto g_0_xxxxyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 196);

    auto g_0_xxxxyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 197);

    auto g_0_xxxxyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 198);

    auto g_0_xxxxyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 199);

    auto g_0_xxxxyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 200);

    auto g_0_xxxxyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 201);

    auto g_0_xxxxyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 202);

    auto g_0_xxxxyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 203);

    auto g_0_xxxxyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 204);

    auto g_0_xxxxyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 205);

    auto g_0_xxxxyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 206);

    auto g_0_xxxxyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 207);

    auto g_0_xxxxyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 208);

    auto g_0_xxxxyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 209);

    auto g_0_xxxxyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 210);

    auto g_0_xxxxyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 211);

    auto g_0_xxxxyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 212);

    auto g_0_xxxxyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 213);

    auto g_0_xxxxyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 214);

    auto g_0_xxxxyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 215);

    auto g_0_xxxxyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 216);

    auto g_0_xxxxyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 217);

    auto g_0_xxxxyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 218);

    auto g_0_xxxxyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 219);

    auto g_0_xxxxyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 220);

    auto g_0_xxxxyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 221);

    auto g_0_xxxxyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 222);

    auto g_0_xxxxyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 223);

#pragma omp simd aligned(g_0_xxxxyy_0_xxxxx_1,       \
                             g_0_xxxxyy_0_xxxxxx_0,  \
                             g_0_xxxxyy_0_xxxxxx_1,  \
                             g_0_xxxxyy_0_xxxxxy_0,  \
                             g_0_xxxxyy_0_xxxxxy_1,  \
                             g_0_xxxxyy_0_xxxxxz_0,  \
                             g_0_xxxxyy_0_xxxxxz_1,  \
                             g_0_xxxxyy_0_xxxxy_1,   \
                             g_0_xxxxyy_0_xxxxyy_0,  \
                             g_0_xxxxyy_0_xxxxyy_1,  \
                             g_0_xxxxyy_0_xxxxyz_0,  \
                             g_0_xxxxyy_0_xxxxyz_1,  \
                             g_0_xxxxyy_0_xxxxz_1,   \
                             g_0_xxxxyy_0_xxxxzz_0,  \
                             g_0_xxxxyy_0_xxxxzz_1,  \
                             g_0_xxxxyy_0_xxxyy_1,   \
                             g_0_xxxxyy_0_xxxyyy_0,  \
                             g_0_xxxxyy_0_xxxyyy_1,  \
                             g_0_xxxxyy_0_xxxyyz_0,  \
                             g_0_xxxxyy_0_xxxyyz_1,  \
                             g_0_xxxxyy_0_xxxyz_1,   \
                             g_0_xxxxyy_0_xxxyzz_0,  \
                             g_0_xxxxyy_0_xxxyzz_1,  \
                             g_0_xxxxyy_0_xxxzz_1,   \
                             g_0_xxxxyy_0_xxxzzz_0,  \
                             g_0_xxxxyy_0_xxxzzz_1,  \
                             g_0_xxxxyy_0_xxyyy_1,   \
                             g_0_xxxxyy_0_xxyyyy_0,  \
                             g_0_xxxxyy_0_xxyyyy_1,  \
                             g_0_xxxxyy_0_xxyyyz_0,  \
                             g_0_xxxxyy_0_xxyyyz_1,  \
                             g_0_xxxxyy_0_xxyyz_1,   \
                             g_0_xxxxyy_0_xxyyzz_0,  \
                             g_0_xxxxyy_0_xxyyzz_1,  \
                             g_0_xxxxyy_0_xxyzz_1,   \
                             g_0_xxxxyy_0_xxyzzz_0,  \
                             g_0_xxxxyy_0_xxyzzz_1,  \
                             g_0_xxxxyy_0_xxzzz_1,   \
                             g_0_xxxxyy_0_xxzzzz_0,  \
                             g_0_xxxxyy_0_xxzzzz_1,  \
                             g_0_xxxxyy_0_xyyyy_1,   \
                             g_0_xxxxyy_0_xyyyyy_0,  \
                             g_0_xxxxyy_0_xyyyyy_1,  \
                             g_0_xxxxyy_0_xyyyyz_0,  \
                             g_0_xxxxyy_0_xyyyyz_1,  \
                             g_0_xxxxyy_0_xyyyz_1,   \
                             g_0_xxxxyy_0_xyyyzz_0,  \
                             g_0_xxxxyy_0_xyyyzz_1,  \
                             g_0_xxxxyy_0_xyyzz_1,   \
                             g_0_xxxxyy_0_xyyzzz_0,  \
                             g_0_xxxxyy_0_xyyzzz_1,  \
                             g_0_xxxxyy_0_xyzzz_1,   \
                             g_0_xxxxyy_0_xyzzzz_0,  \
                             g_0_xxxxyy_0_xyzzzz_1,  \
                             g_0_xxxxyy_0_xzzzz_1,   \
                             g_0_xxxxyy_0_xzzzzz_0,  \
                             g_0_xxxxyy_0_xzzzzz_1,  \
                             g_0_xxxxyy_0_yyyyy_1,   \
                             g_0_xxxxyy_0_yyyyyy_0,  \
                             g_0_xxxxyy_0_yyyyyy_1,  \
                             g_0_xxxxyy_0_yyyyyz_0,  \
                             g_0_xxxxyy_0_yyyyyz_1,  \
                             g_0_xxxxyy_0_yyyyz_1,   \
                             g_0_xxxxyy_0_yyyyzz_0,  \
                             g_0_xxxxyy_0_yyyyzz_1,  \
                             g_0_xxxxyy_0_yyyzz_1,   \
                             g_0_xxxxyy_0_yyyzzz_0,  \
                             g_0_xxxxyy_0_yyyzzz_1,  \
                             g_0_xxxxyy_0_yyzzz_1,   \
                             g_0_xxxxyy_0_yyzzzz_0,  \
                             g_0_xxxxyy_0_yyzzzz_1,  \
                             g_0_xxxxyy_0_yzzzz_1,   \
                             g_0_xxxxyy_0_yzzzzz_0,  \
                             g_0_xxxxyy_0_yzzzzz_1,  \
                             g_0_xxxxyy_0_zzzzz_1,   \
                             g_0_xxxxyy_0_zzzzzz_0,  \
                             g_0_xxxxyy_0_zzzzzz_1,  \
                             g_0_xxxxyyz_0_xxxxxx_0, \
                             g_0_xxxxyyz_0_xxxxxy_0, \
                             g_0_xxxxyyz_0_xxxxxz_0, \
                             g_0_xxxxyyz_0_xxxxyy_0, \
                             g_0_xxxxyyz_0_xxxxyz_0, \
                             g_0_xxxxyyz_0_xxxxzz_0, \
                             g_0_xxxxyyz_0_xxxyyy_0, \
                             g_0_xxxxyyz_0_xxxyyz_0, \
                             g_0_xxxxyyz_0_xxxyzz_0, \
                             g_0_xxxxyyz_0_xxxzzz_0, \
                             g_0_xxxxyyz_0_xxyyyy_0, \
                             g_0_xxxxyyz_0_xxyyyz_0, \
                             g_0_xxxxyyz_0_xxyyzz_0, \
                             g_0_xxxxyyz_0_xxyzzz_0, \
                             g_0_xxxxyyz_0_xxzzzz_0, \
                             g_0_xxxxyyz_0_xyyyyy_0, \
                             g_0_xxxxyyz_0_xyyyyz_0, \
                             g_0_xxxxyyz_0_xyyyzz_0, \
                             g_0_xxxxyyz_0_xyyzzz_0, \
                             g_0_xxxxyyz_0_xyzzzz_0, \
                             g_0_xxxxyyz_0_xzzzzz_0, \
                             g_0_xxxxyyz_0_yyyyyy_0, \
                             g_0_xxxxyyz_0_yyyyyz_0, \
                             g_0_xxxxyyz_0_yyyyzz_0, \
                             g_0_xxxxyyz_0_yyyzzz_0, \
                             g_0_xxxxyyz_0_yyzzzz_0, \
                             g_0_xxxxyyz_0_yzzzzz_0, \
                             g_0_xxxxyyz_0_zzzzzz_0, \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyz_0_xxxxxx_0[i] = g_0_xxxxyy_0_xxxxxx_0[i] * pb_z + g_0_xxxxyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxxy_0[i] = g_0_xxxxyy_0_xxxxxy_0[i] * pb_z + g_0_xxxxyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxxz_0[i] = g_0_xxxxyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxxz_0[i] * pb_z + g_0_xxxxyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxyy_0[i] = g_0_xxxxyy_0_xxxxyy_0[i] * pb_z + g_0_xxxxyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxyz_0[i] = g_0_xxxxyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxyz_0[i] * pb_z + g_0_xxxxyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxxzz_0[i] = 2.0 * g_0_xxxxyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxxzz_0[i] * pb_z + g_0_xxxxyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxyyy_0[i] = g_0_xxxxyy_0_xxxyyy_0[i] * pb_z + g_0_xxxxyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxyyz_0[i] = g_0_xxxxyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyyz_0[i] * pb_z + g_0_xxxxyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxyzz_0[i] = 2.0 * g_0_xxxxyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxyzz_0[i] * pb_z + g_0_xxxxyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxxzzz_0[i] = 3.0 * g_0_xxxxyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxxzzz_0[i] * pb_z + g_0_xxxxyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyyyy_0[i] = g_0_xxxxyy_0_xxyyyy_0[i] * pb_z + g_0_xxxxyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyyyz_0[i] = g_0_xxxxyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyyz_0[i] * pb_z + g_0_xxxxyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyyzz_0[i] = 2.0 * g_0_xxxxyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyyzz_0[i] * pb_z + g_0_xxxxyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxyzzz_0[i] = 3.0 * g_0_xxxxyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxyzzz_0[i] * pb_z + g_0_xxxxyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xxzzzz_0[i] = 4.0 * g_0_xxxxyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xxzzzz_0[i] * pb_z + g_0_xxxxyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyyyy_0[i] = g_0_xxxxyy_0_xyyyyy_0[i] * pb_z + g_0_xxxxyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyyyz_0[i] = g_0_xxxxyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyyz_0[i] * pb_z + g_0_xxxxyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyyzz_0[i] = 2.0 * g_0_xxxxyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyyzz_0[i] * pb_z + g_0_xxxxyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyyzzz_0[i] = 3.0 * g_0_xxxxyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyyzzz_0[i] * pb_z + g_0_xxxxyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xyzzzz_0[i] = 4.0 * g_0_xxxxyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xyzzzz_0[i] * pb_z + g_0_xxxxyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_xzzzzz_0[i] = 5.0 * g_0_xxxxyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_xzzzzz_0[i] * pb_z + g_0_xxxxyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyyyy_0[i] = g_0_xxxxyy_0_yyyyyy_0[i] * pb_z + g_0_xxxxyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyyyz_0[i] = g_0_xxxxyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyyyyz_0[i] * pb_z + g_0_xxxxyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyyzz_0[i] = 2.0 * g_0_xxxxyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyyyzz_0[i] * pb_z + g_0_xxxxyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyyzzz_0[i] = 3.0 * g_0_xxxxyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyyzzz_0[i] * pb_z + g_0_xxxxyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yyzzzz_0[i] = 4.0 * g_0_xxxxyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yyzzzz_0[i] * pb_z + g_0_xxxxyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_yzzzzz_0[i] = 5.0 * g_0_xxxxyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_yzzzzz_0[i] * pb_z + g_0_xxxxyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyz_0_zzzzzz_0[i] = 6.0 * g_0_xxxxyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxyy_0_zzzzzz_0[i] * pb_z + g_0_xxxxyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 224-252 components of targeted buffer : SKSI

    auto g_0_xxxxyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 224);

    auto g_0_xxxxyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 225);

    auto g_0_xxxxyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 226);

    auto g_0_xxxxyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 227);

    auto g_0_xxxxyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 228);

    auto g_0_xxxxyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 229);

    auto g_0_xxxxyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 230);

    auto g_0_xxxxyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 231);

    auto g_0_xxxxyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 232);

    auto g_0_xxxxyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 233);

    auto g_0_xxxxyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 234);

    auto g_0_xxxxyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 235);

    auto g_0_xxxxyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 236);

    auto g_0_xxxxyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 237);

    auto g_0_xxxxyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 238);

    auto g_0_xxxxyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 239);

    auto g_0_xxxxyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 240);

    auto g_0_xxxxyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 241);

    auto g_0_xxxxyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 242);

    auto g_0_xxxxyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 243);

    auto g_0_xxxxyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 244);

    auto g_0_xxxxyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 245);

    auto g_0_xxxxyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 246);

    auto g_0_xxxxyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 247);

    auto g_0_xxxxyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 248);

    auto g_0_xxxxyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 249);

    auto g_0_xxxxyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 250);

    auto g_0_xxxxyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 251);

#pragma omp simd aligned(g_0_xxxxyzz_0_xxxxxx_0,     \
                             g_0_xxxxyzz_0_xxxxxy_0, \
                             g_0_xxxxyzz_0_xxxxxz_0, \
                             g_0_xxxxyzz_0_xxxxyy_0, \
                             g_0_xxxxyzz_0_xxxxyz_0, \
                             g_0_xxxxyzz_0_xxxxzz_0, \
                             g_0_xxxxyzz_0_xxxyyy_0, \
                             g_0_xxxxyzz_0_xxxyyz_0, \
                             g_0_xxxxyzz_0_xxxyzz_0, \
                             g_0_xxxxyzz_0_xxxzzz_0, \
                             g_0_xxxxyzz_0_xxyyyy_0, \
                             g_0_xxxxyzz_0_xxyyyz_0, \
                             g_0_xxxxyzz_0_xxyyzz_0, \
                             g_0_xxxxyzz_0_xxyzzz_0, \
                             g_0_xxxxyzz_0_xxzzzz_0, \
                             g_0_xxxxyzz_0_xyyyyy_0, \
                             g_0_xxxxyzz_0_xyyyyz_0, \
                             g_0_xxxxyzz_0_xyyyzz_0, \
                             g_0_xxxxyzz_0_xyyzzz_0, \
                             g_0_xxxxyzz_0_xyzzzz_0, \
                             g_0_xxxxyzz_0_xzzzzz_0, \
                             g_0_xxxxyzz_0_yyyyyy_0, \
                             g_0_xxxxyzz_0_yyyyyz_0, \
                             g_0_xxxxyzz_0_yyyyzz_0, \
                             g_0_xxxxyzz_0_yyyzzz_0, \
                             g_0_xxxxyzz_0_yyzzzz_0, \
                             g_0_xxxxyzz_0_yzzzzz_0, \
                             g_0_xxxxyzz_0_zzzzzz_0, \
                             g_0_xxxxzz_0_xxxxx_1,   \
                             g_0_xxxxzz_0_xxxxxx_0,  \
                             g_0_xxxxzz_0_xxxxxx_1,  \
                             g_0_xxxxzz_0_xxxxxy_0,  \
                             g_0_xxxxzz_0_xxxxxy_1,  \
                             g_0_xxxxzz_0_xxxxxz_0,  \
                             g_0_xxxxzz_0_xxxxxz_1,  \
                             g_0_xxxxzz_0_xxxxy_1,   \
                             g_0_xxxxzz_0_xxxxyy_0,  \
                             g_0_xxxxzz_0_xxxxyy_1,  \
                             g_0_xxxxzz_0_xxxxyz_0,  \
                             g_0_xxxxzz_0_xxxxyz_1,  \
                             g_0_xxxxzz_0_xxxxz_1,   \
                             g_0_xxxxzz_0_xxxxzz_0,  \
                             g_0_xxxxzz_0_xxxxzz_1,  \
                             g_0_xxxxzz_0_xxxyy_1,   \
                             g_0_xxxxzz_0_xxxyyy_0,  \
                             g_0_xxxxzz_0_xxxyyy_1,  \
                             g_0_xxxxzz_0_xxxyyz_0,  \
                             g_0_xxxxzz_0_xxxyyz_1,  \
                             g_0_xxxxzz_0_xxxyz_1,   \
                             g_0_xxxxzz_0_xxxyzz_0,  \
                             g_0_xxxxzz_0_xxxyzz_1,  \
                             g_0_xxxxzz_0_xxxzz_1,   \
                             g_0_xxxxzz_0_xxxzzz_0,  \
                             g_0_xxxxzz_0_xxxzzz_1,  \
                             g_0_xxxxzz_0_xxyyy_1,   \
                             g_0_xxxxzz_0_xxyyyy_0,  \
                             g_0_xxxxzz_0_xxyyyy_1,  \
                             g_0_xxxxzz_0_xxyyyz_0,  \
                             g_0_xxxxzz_0_xxyyyz_1,  \
                             g_0_xxxxzz_0_xxyyz_1,   \
                             g_0_xxxxzz_0_xxyyzz_0,  \
                             g_0_xxxxzz_0_xxyyzz_1,  \
                             g_0_xxxxzz_0_xxyzz_1,   \
                             g_0_xxxxzz_0_xxyzzz_0,  \
                             g_0_xxxxzz_0_xxyzzz_1,  \
                             g_0_xxxxzz_0_xxzzz_1,   \
                             g_0_xxxxzz_0_xxzzzz_0,  \
                             g_0_xxxxzz_0_xxzzzz_1,  \
                             g_0_xxxxzz_0_xyyyy_1,   \
                             g_0_xxxxzz_0_xyyyyy_0,  \
                             g_0_xxxxzz_0_xyyyyy_1,  \
                             g_0_xxxxzz_0_xyyyyz_0,  \
                             g_0_xxxxzz_0_xyyyyz_1,  \
                             g_0_xxxxzz_0_xyyyz_1,   \
                             g_0_xxxxzz_0_xyyyzz_0,  \
                             g_0_xxxxzz_0_xyyyzz_1,  \
                             g_0_xxxxzz_0_xyyzz_1,   \
                             g_0_xxxxzz_0_xyyzzz_0,  \
                             g_0_xxxxzz_0_xyyzzz_1,  \
                             g_0_xxxxzz_0_xyzzz_1,   \
                             g_0_xxxxzz_0_xyzzzz_0,  \
                             g_0_xxxxzz_0_xyzzzz_1,  \
                             g_0_xxxxzz_0_xzzzz_1,   \
                             g_0_xxxxzz_0_xzzzzz_0,  \
                             g_0_xxxxzz_0_xzzzzz_1,  \
                             g_0_xxxxzz_0_yyyyy_1,   \
                             g_0_xxxxzz_0_yyyyyy_0,  \
                             g_0_xxxxzz_0_yyyyyy_1,  \
                             g_0_xxxxzz_0_yyyyyz_0,  \
                             g_0_xxxxzz_0_yyyyyz_1,  \
                             g_0_xxxxzz_0_yyyyz_1,   \
                             g_0_xxxxzz_0_yyyyzz_0,  \
                             g_0_xxxxzz_0_yyyyzz_1,  \
                             g_0_xxxxzz_0_yyyzz_1,   \
                             g_0_xxxxzz_0_yyyzzz_0,  \
                             g_0_xxxxzz_0_yyyzzz_1,  \
                             g_0_xxxxzz_0_yyzzz_1,   \
                             g_0_xxxxzz_0_yyzzzz_0,  \
                             g_0_xxxxzz_0_yyzzzz_1,  \
                             g_0_xxxxzz_0_yzzzz_1,   \
                             g_0_xxxxzz_0_yzzzzz_0,  \
                             g_0_xxxxzz_0_yzzzzz_1,  \
                             g_0_xxxxzz_0_zzzzz_1,   \
                             g_0_xxxxzz_0_zzzzzz_0,  \
                             g_0_xxxxzz_0_zzzzzz_1,  \
                             wp_y,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzz_0_xxxxxx_0[i] = g_0_xxxxzz_0_xxxxxx_0[i] * pb_y + g_0_xxxxzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxxy_0[i] = g_0_xxxxzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxxy_0[i] * pb_y + g_0_xxxxzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxxz_0[i] = g_0_xxxxzz_0_xxxxxz_0[i] * pb_y + g_0_xxxxzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxyy_0[i] = 2.0 * g_0_xxxxzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxyy_0[i] * pb_y + g_0_xxxxzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxyz_0[i] = g_0_xxxxzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxxyz_0[i] * pb_y + g_0_xxxxzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxxzz_0[i] = g_0_xxxxzz_0_xxxxzz_0[i] * pb_y + g_0_xxxxzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxyyy_0[i] = 3.0 * g_0_xxxxzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyyy_0[i] * pb_y + g_0_xxxxzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxyyz_0[i] = 2.0 * g_0_xxxxzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyyz_0[i] * pb_y + g_0_xxxxzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxyzz_0[i] = g_0_xxxxzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxxyzz_0[i] * pb_y + g_0_xxxxzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxxzzz_0[i] = g_0_xxxxzz_0_xxxzzz_0[i] * pb_y + g_0_xxxxzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyyyy_0[i] = 4.0 * g_0_xxxxzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyyy_0[i] * pb_y + g_0_xxxxzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyyyz_0[i] = 3.0 * g_0_xxxxzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyyz_0[i] * pb_y + g_0_xxxxzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyyzz_0[i] = 2.0 * g_0_xxxxzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyyzz_0[i] * pb_y + g_0_xxxxzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxyzzz_0[i] = g_0_xxxxzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xxyzzz_0[i] * pb_y + g_0_xxxxzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xxzzzz_0[i] = g_0_xxxxzz_0_xxzzzz_0[i] * pb_y + g_0_xxxxzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyyyy_0[i] = 5.0 * g_0_xxxxzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyyy_0[i] * pb_y + g_0_xxxxzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyyyz_0[i] = 4.0 * g_0_xxxxzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyyz_0[i] * pb_y + g_0_xxxxzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyyzz_0[i] = 3.0 * g_0_xxxxzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyyzz_0[i] * pb_y + g_0_xxxxzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyyzzz_0[i] = 2.0 * g_0_xxxxzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyyzzz_0[i] * pb_y + g_0_xxxxzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xyzzzz_0[i] = g_0_xxxxzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_xyzzzz_0[i] * pb_y + g_0_xxxxzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_xzzzzz_0[i] = g_0_xxxxzz_0_xzzzzz_0[i] * pb_y + g_0_xxxxzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyyyy_0[i] = 6.0 * g_0_xxxxzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyyyy_0[i] * pb_y + g_0_xxxxzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyyyz_0[i] = 5.0 * g_0_xxxxzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyyyz_0[i] * pb_y + g_0_xxxxzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyyzz_0[i] = 4.0 * g_0_xxxxzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyyzz_0[i] * pb_y + g_0_xxxxzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyyzzz_0[i] = 3.0 * g_0_xxxxzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyyzzz_0[i] * pb_y + g_0_xxxxzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yyzzzz_0[i] = 2.0 * g_0_xxxxzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yyzzzz_0[i] * pb_y + g_0_xxxxzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_yzzzzz_0[i] = g_0_xxxxzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxzz_0_yzzzzz_0[i] * pb_y + g_0_xxxxzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzz_0_zzzzzz_0[i] = g_0_xxxxzz_0_zzzzzz_0[i] * pb_y + g_0_xxxxzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 252-280 components of targeted buffer : SKSI

    auto g_0_xxxxzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 252);

    auto g_0_xxxxzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 253);

    auto g_0_xxxxzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 254);

    auto g_0_xxxxzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 255);

    auto g_0_xxxxzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 256);

    auto g_0_xxxxzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 257);

    auto g_0_xxxxzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 258);

    auto g_0_xxxxzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 259);

    auto g_0_xxxxzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 260);

    auto g_0_xxxxzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 261);

    auto g_0_xxxxzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 262);

    auto g_0_xxxxzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 263);

    auto g_0_xxxxzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 264);

    auto g_0_xxxxzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 265);

    auto g_0_xxxxzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 266);

    auto g_0_xxxxzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 267);

    auto g_0_xxxxzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 268);

    auto g_0_xxxxzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 269);

    auto g_0_xxxxzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 270);

    auto g_0_xxxxzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 271);

    auto g_0_xxxxzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 272);

    auto g_0_xxxxzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 273);

    auto g_0_xxxxzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 274);

    auto g_0_xxxxzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 275);

    auto g_0_xxxxzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 276);

    auto g_0_xxxxzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 277);

    auto g_0_xxxxzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 278);

    auto g_0_xxxxzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 279);

#pragma omp simd aligned(g_0_xxxxz_0_xxxxxx_0,       \
                             g_0_xxxxz_0_xxxxxx_1,   \
                             g_0_xxxxz_0_xxxxxy_0,   \
                             g_0_xxxxz_0_xxxxxy_1,   \
                             g_0_xxxxz_0_xxxxyy_0,   \
                             g_0_xxxxz_0_xxxxyy_1,   \
                             g_0_xxxxz_0_xxxyyy_0,   \
                             g_0_xxxxz_0_xxxyyy_1,   \
                             g_0_xxxxz_0_xxyyyy_0,   \
                             g_0_xxxxz_0_xxyyyy_1,   \
                             g_0_xxxxz_0_xyyyyy_0,   \
                             g_0_xxxxz_0_xyyyyy_1,   \
                             g_0_xxxxzz_0_xxxxxx_0,  \
                             g_0_xxxxzz_0_xxxxxx_1,  \
                             g_0_xxxxzz_0_xxxxxy_0,  \
                             g_0_xxxxzz_0_xxxxxy_1,  \
                             g_0_xxxxzz_0_xxxxyy_0,  \
                             g_0_xxxxzz_0_xxxxyy_1,  \
                             g_0_xxxxzz_0_xxxyyy_0,  \
                             g_0_xxxxzz_0_xxxyyy_1,  \
                             g_0_xxxxzz_0_xxyyyy_0,  \
                             g_0_xxxxzz_0_xxyyyy_1,  \
                             g_0_xxxxzz_0_xyyyyy_0,  \
                             g_0_xxxxzz_0_xyyyyy_1,  \
                             g_0_xxxxzzz_0_xxxxxx_0, \
                             g_0_xxxxzzz_0_xxxxxy_0, \
                             g_0_xxxxzzz_0_xxxxxz_0, \
                             g_0_xxxxzzz_0_xxxxyy_0, \
                             g_0_xxxxzzz_0_xxxxyz_0, \
                             g_0_xxxxzzz_0_xxxxzz_0, \
                             g_0_xxxxzzz_0_xxxyyy_0, \
                             g_0_xxxxzzz_0_xxxyyz_0, \
                             g_0_xxxxzzz_0_xxxyzz_0, \
                             g_0_xxxxzzz_0_xxxzzz_0, \
                             g_0_xxxxzzz_0_xxyyyy_0, \
                             g_0_xxxxzzz_0_xxyyyz_0, \
                             g_0_xxxxzzz_0_xxyyzz_0, \
                             g_0_xxxxzzz_0_xxyzzz_0, \
                             g_0_xxxxzzz_0_xxzzzz_0, \
                             g_0_xxxxzzz_0_xyyyyy_0, \
                             g_0_xxxxzzz_0_xyyyyz_0, \
                             g_0_xxxxzzz_0_xyyyzz_0, \
                             g_0_xxxxzzz_0_xyyzzz_0, \
                             g_0_xxxxzzz_0_xyzzzz_0, \
                             g_0_xxxxzzz_0_xzzzzz_0, \
                             g_0_xxxxzzz_0_yyyyyy_0, \
                             g_0_xxxxzzz_0_yyyyyz_0, \
                             g_0_xxxxzzz_0_yyyyzz_0, \
                             g_0_xxxxzzz_0_yyyzzz_0, \
                             g_0_xxxxzzz_0_yyzzzz_0, \
                             g_0_xxxxzzz_0_yzzzzz_0, \
                             g_0_xxxxzzz_0_zzzzzz_0, \
                             g_0_xxxzzz_0_xxxxxz_0,  \
                             g_0_xxxzzz_0_xxxxxz_1,  \
                             g_0_xxxzzz_0_xxxxyz_0,  \
                             g_0_xxxzzz_0_xxxxyz_1,  \
                             g_0_xxxzzz_0_xxxxz_1,   \
                             g_0_xxxzzz_0_xxxxzz_0,  \
                             g_0_xxxzzz_0_xxxxzz_1,  \
                             g_0_xxxzzz_0_xxxyyz_0,  \
                             g_0_xxxzzz_0_xxxyyz_1,  \
                             g_0_xxxzzz_0_xxxyz_1,   \
                             g_0_xxxzzz_0_xxxyzz_0,  \
                             g_0_xxxzzz_0_xxxyzz_1,  \
                             g_0_xxxzzz_0_xxxzz_1,   \
                             g_0_xxxzzz_0_xxxzzz_0,  \
                             g_0_xxxzzz_0_xxxzzz_1,  \
                             g_0_xxxzzz_0_xxyyyz_0,  \
                             g_0_xxxzzz_0_xxyyyz_1,  \
                             g_0_xxxzzz_0_xxyyz_1,   \
                             g_0_xxxzzz_0_xxyyzz_0,  \
                             g_0_xxxzzz_0_xxyyzz_1,  \
                             g_0_xxxzzz_0_xxyzz_1,   \
                             g_0_xxxzzz_0_xxyzzz_0,  \
                             g_0_xxxzzz_0_xxyzzz_1,  \
                             g_0_xxxzzz_0_xxzzz_1,   \
                             g_0_xxxzzz_0_xxzzzz_0,  \
                             g_0_xxxzzz_0_xxzzzz_1,  \
                             g_0_xxxzzz_0_xyyyyz_0,  \
                             g_0_xxxzzz_0_xyyyyz_1,  \
                             g_0_xxxzzz_0_xyyyz_1,   \
                             g_0_xxxzzz_0_xyyyzz_0,  \
                             g_0_xxxzzz_0_xyyyzz_1,  \
                             g_0_xxxzzz_0_xyyzz_1,   \
                             g_0_xxxzzz_0_xyyzzz_0,  \
                             g_0_xxxzzz_0_xyyzzz_1,  \
                             g_0_xxxzzz_0_xyzzz_1,   \
                             g_0_xxxzzz_0_xyzzzz_0,  \
                             g_0_xxxzzz_0_xyzzzz_1,  \
                             g_0_xxxzzz_0_xzzzz_1,   \
                             g_0_xxxzzz_0_xzzzzz_0,  \
                             g_0_xxxzzz_0_xzzzzz_1,  \
                             g_0_xxxzzz_0_yyyyyy_0,  \
                             g_0_xxxzzz_0_yyyyyy_1,  \
                             g_0_xxxzzz_0_yyyyyz_0,  \
                             g_0_xxxzzz_0_yyyyyz_1,  \
                             g_0_xxxzzz_0_yyyyz_1,   \
                             g_0_xxxzzz_0_yyyyzz_0,  \
                             g_0_xxxzzz_0_yyyyzz_1,  \
                             g_0_xxxzzz_0_yyyzz_1,   \
                             g_0_xxxzzz_0_yyyzzz_0,  \
                             g_0_xxxzzz_0_yyyzzz_1,  \
                             g_0_xxxzzz_0_yyzzz_1,   \
                             g_0_xxxzzz_0_yyzzzz_0,  \
                             g_0_xxxzzz_0_yyzzzz_1,  \
                             g_0_xxxzzz_0_yzzzz_1,   \
                             g_0_xxxzzz_0_yzzzzz_0,  \
                             g_0_xxxzzz_0_yzzzzz_1,  \
                             g_0_xxxzzz_0_zzzzz_1,   \
                             g_0_xxxzzz_0_zzzzzz_0,  \
                             g_0_xxxzzz_0_zzzzzz_1,  \
                             g_0_xxzzz_0_xxxxxz_0,   \
                             g_0_xxzzz_0_xxxxxz_1,   \
                             g_0_xxzzz_0_xxxxyz_0,   \
                             g_0_xxzzz_0_xxxxyz_1,   \
                             g_0_xxzzz_0_xxxxzz_0,   \
                             g_0_xxzzz_0_xxxxzz_1,   \
                             g_0_xxzzz_0_xxxyyz_0,   \
                             g_0_xxzzz_0_xxxyyz_1,   \
                             g_0_xxzzz_0_xxxyzz_0,   \
                             g_0_xxzzz_0_xxxyzz_1,   \
                             g_0_xxzzz_0_xxxzzz_0,   \
                             g_0_xxzzz_0_xxxzzz_1,   \
                             g_0_xxzzz_0_xxyyyz_0,   \
                             g_0_xxzzz_0_xxyyyz_1,   \
                             g_0_xxzzz_0_xxyyzz_0,   \
                             g_0_xxzzz_0_xxyyzz_1,   \
                             g_0_xxzzz_0_xxyzzz_0,   \
                             g_0_xxzzz_0_xxyzzz_1,   \
                             g_0_xxzzz_0_xxzzzz_0,   \
                             g_0_xxzzz_0_xxzzzz_1,   \
                             g_0_xxzzz_0_xyyyyz_0,   \
                             g_0_xxzzz_0_xyyyyz_1,   \
                             g_0_xxzzz_0_xyyyzz_0,   \
                             g_0_xxzzz_0_xyyyzz_1,   \
                             g_0_xxzzz_0_xyyzzz_0,   \
                             g_0_xxzzz_0_xyyzzz_1,   \
                             g_0_xxzzz_0_xyzzzz_0,   \
                             g_0_xxzzz_0_xyzzzz_1,   \
                             g_0_xxzzz_0_xzzzzz_0,   \
                             g_0_xxzzz_0_xzzzzz_1,   \
                             g_0_xxzzz_0_yyyyyy_0,   \
                             g_0_xxzzz_0_yyyyyy_1,   \
                             g_0_xxzzz_0_yyyyyz_0,   \
                             g_0_xxzzz_0_yyyyyz_1,   \
                             g_0_xxzzz_0_yyyyzz_0,   \
                             g_0_xxzzz_0_yyyyzz_1,   \
                             g_0_xxzzz_0_yyyzzz_0,   \
                             g_0_xxzzz_0_yyyzzz_1,   \
                             g_0_xxzzz_0_yyzzzz_0,   \
                             g_0_xxzzz_0_yyzzzz_1,   \
                             g_0_xxzzz_0_yzzzzz_0,   \
                             g_0_xxzzz_0_yzzzzz_1,   \
                             g_0_xxzzz_0_zzzzzz_0,   \
                             g_0_xxzzz_0_zzzzzz_1,   \
                             wp_x,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzzz_0_xxxxxx_0[i] = 2.0 * g_0_xxxxz_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxxxx_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_xxxxxx_0[i] * pb_z + g_0_xxxxzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxxxy_0[i] = 2.0 * g_0_xxxxz_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxxxy_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_xxxxxy_0[i] * pb_z + g_0_xxxxzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxxxz_0[i] = 3.0 * g_0_xxzzz_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                    5.0 * g_0_xxxzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxxz_0[i] * pb_x + g_0_xxxzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxxyy_0[i] = 2.0 * g_0_xxxxz_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxxyy_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_xxxxyy_0[i] * pb_z + g_0_xxxxzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxxyz_0[i] = 3.0 * g_0_xxzzz_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxxzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxyz_0[i] * pb_x + g_0_xxxzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxxzz_0[i] = 3.0 * g_0_xxzzz_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxxzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxzz_0[i] * pb_x + g_0_xxxzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxyyy_0[i] = 2.0 * g_0_xxxxz_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxxyyy_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_xxxyyy_0[i] * pb_z + g_0_xxxxzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxxyyz_0[i] = 3.0 * g_0_xxzzz_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyyz_0[i] * pb_x + g_0_xxxzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxyzz_0[i] = 3.0 * g_0_xxzzz_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyzz_0[i] * pb_x + g_0_xxxzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxxzzz_0[i] = 3.0 * g_0_xxzzz_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxxzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxzzz_0[i] * pb_x + g_0_xxxzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxyyyy_0[i] = 2.0 * g_0_xxxxz_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xxyyyy_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_xxyyyy_0[i] * pb_z + g_0_xxxxzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xxyyyz_0[i] = 3.0 * g_0_xxzzz_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyyz_0[i] * pb_x + g_0_xxxzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxyyzz_0[i] = 3.0 * g_0_xxzzz_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyzz_0[i] * pb_x + g_0_xxxzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxyzzz_0[i] = 3.0 * g_0_xxzzz_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyzzz_0[i] * pb_x + g_0_xxxzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xxzzzz_0[i] = 3.0 * g_0_xxzzz_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxxzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxzzzz_0[i] * pb_x + g_0_xxxzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyyyyy_0[i] = 2.0 * g_0_xxxxz_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxz_0_xyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxxxzz_0_xyyyyy_0[i] * pb_z + g_0_xxxxzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxzzz_0_xyyyyz_0[i] = 3.0 * g_0_xxzzz_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyyz_0[i] * pb_x + g_0_xxxzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyyyzz_0[i] = 3.0 * g_0_xxzzz_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyzz_0[i] * pb_x + g_0_xxxzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyyzzz_0[i] = 3.0 * g_0_xxzzz_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyzzz_0[i] * pb_x + g_0_xxxzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xyzzzz_0[i] = 3.0 * g_0_xxzzz_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyzzzz_0[i] * pb_x + g_0_xxxzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_xzzzzz_0[i] = 3.0 * g_0_xxzzz_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xzzzzz_0[i] * pb_x + g_0_xxxzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyyyy_0[i] = 3.0 * g_0_xxzzz_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_yyyyyy_0[i] * pb_x + g_0_xxxzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyyyz_0[i] = 3.0 * g_0_xxzzz_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_yyyyyz_0[i] * pb_x + g_0_xxxzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyyzz_0[i] = 3.0 * g_0_xxzzz_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_yyyyzz_0[i] * pb_x + g_0_xxxzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyyzzz_0[i] = 3.0 * g_0_xxzzz_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_yyyzzz_0[i] * pb_x + g_0_xxxzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yyzzzz_0[i] = 3.0 * g_0_xxzzz_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_yyzzzz_0[i] * pb_x + g_0_xxxzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_yzzzzz_0[i] = 3.0 * g_0_xxzzz_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_yzzzzz_0[i] * pb_x + g_0_xxxzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzz_0_zzzzzz_0[i] = 3.0 * g_0_xxzzz_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_zzzzzz_0[i] * pb_x + g_0_xxxzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 280-308 components of targeted buffer : SKSI

    auto g_0_xxxyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 280);

    auto g_0_xxxyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 281);

    auto g_0_xxxyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 282);

    auto g_0_xxxyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 283);

    auto g_0_xxxyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 284);

    auto g_0_xxxyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 285);

    auto g_0_xxxyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 286);

    auto g_0_xxxyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 287);

    auto g_0_xxxyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 288);

    auto g_0_xxxyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 289);

    auto g_0_xxxyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 290);

    auto g_0_xxxyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 291);

    auto g_0_xxxyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 292);

    auto g_0_xxxyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 293);

    auto g_0_xxxyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 294);

    auto g_0_xxxyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 295);

    auto g_0_xxxyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 296);

    auto g_0_xxxyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 297);

    auto g_0_xxxyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 298);

    auto g_0_xxxyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 299);

    auto g_0_xxxyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 300);

    auto g_0_xxxyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 301);

    auto g_0_xxxyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 302);

    auto g_0_xxxyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 303);

    auto g_0_xxxyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 304);

    auto g_0_xxxyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 305);

    auto g_0_xxxyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 306);

    auto g_0_xxxyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 307);

#pragma omp simd aligned(g_0_xxxyy_0_xxxxxx_0,       \
                             g_0_xxxyy_0_xxxxxx_1,   \
                             g_0_xxxyy_0_xxxxxz_0,   \
                             g_0_xxxyy_0_xxxxxz_1,   \
                             g_0_xxxyy_0_xxxxzz_0,   \
                             g_0_xxxyy_0_xxxxzz_1,   \
                             g_0_xxxyy_0_xxxzzz_0,   \
                             g_0_xxxyy_0_xxxzzz_1,   \
                             g_0_xxxyy_0_xxzzzz_0,   \
                             g_0_xxxyy_0_xxzzzz_1,   \
                             g_0_xxxyy_0_xzzzzz_0,   \
                             g_0_xxxyy_0_xzzzzz_1,   \
                             g_0_xxxyyy_0_xxxxxx_0,  \
                             g_0_xxxyyy_0_xxxxxx_1,  \
                             g_0_xxxyyy_0_xxxxxz_0,  \
                             g_0_xxxyyy_0_xxxxxz_1,  \
                             g_0_xxxyyy_0_xxxxzz_0,  \
                             g_0_xxxyyy_0_xxxxzz_1,  \
                             g_0_xxxyyy_0_xxxzzz_0,  \
                             g_0_xxxyyy_0_xxxzzz_1,  \
                             g_0_xxxyyy_0_xxzzzz_0,  \
                             g_0_xxxyyy_0_xxzzzz_1,  \
                             g_0_xxxyyy_0_xzzzzz_0,  \
                             g_0_xxxyyy_0_xzzzzz_1,  \
                             g_0_xxxyyyy_0_xxxxxx_0, \
                             g_0_xxxyyyy_0_xxxxxy_0, \
                             g_0_xxxyyyy_0_xxxxxz_0, \
                             g_0_xxxyyyy_0_xxxxyy_0, \
                             g_0_xxxyyyy_0_xxxxyz_0, \
                             g_0_xxxyyyy_0_xxxxzz_0, \
                             g_0_xxxyyyy_0_xxxyyy_0, \
                             g_0_xxxyyyy_0_xxxyyz_0, \
                             g_0_xxxyyyy_0_xxxyzz_0, \
                             g_0_xxxyyyy_0_xxxzzz_0, \
                             g_0_xxxyyyy_0_xxyyyy_0, \
                             g_0_xxxyyyy_0_xxyyyz_0, \
                             g_0_xxxyyyy_0_xxyyzz_0, \
                             g_0_xxxyyyy_0_xxyzzz_0, \
                             g_0_xxxyyyy_0_xxzzzz_0, \
                             g_0_xxxyyyy_0_xyyyyy_0, \
                             g_0_xxxyyyy_0_xyyyyz_0, \
                             g_0_xxxyyyy_0_xyyyzz_0, \
                             g_0_xxxyyyy_0_xyyzzz_0, \
                             g_0_xxxyyyy_0_xyzzzz_0, \
                             g_0_xxxyyyy_0_xzzzzz_0, \
                             g_0_xxxyyyy_0_yyyyyy_0, \
                             g_0_xxxyyyy_0_yyyyyz_0, \
                             g_0_xxxyyyy_0_yyyyzz_0, \
                             g_0_xxxyyyy_0_yyyzzz_0, \
                             g_0_xxxyyyy_0_yyzzzz_0, \
                             g_0_xxxyyyy_0_yzzzzz_0, \
                             g_0_xxxyyyy_0_zzzzzz_0, \
                             g_0_xxyyyy_0_xxxxxy_0,  \
                             g_0_xxyyyy_0_xxxxxy_1,  \
                             g_0_xxyyyy_0_xxxxy_1,   \
                             g_0_xxyyyy_0_xxxxyy_0,  \
                             g_0_xxyyyy_0_xxxxyy_1,  \
                             g_0_xxyyyy_0_xxxxyz_0,  \
                             g_0_xxyyyy_0_xxxxyz_1,  \
                             g_0_xxyyyy_0_xxxyy_1,   \
                             g_0_xxyyyy_0_xxxyyy_0,  \
                             g_0_xxyyyy_0_xxxyyy_1,  \
                             g_0_xxyyyy_0_xxxyyz_0,  \
                             g_0_xxyyyy_0_xxxyyz_1,  \
                             g_0_xxyyyy_0_xxxyz_1,   \
                             g_0_xxyyyy_0_xxxyzz_0,  \
                             g_0_xxyyyy_0_xxxyzz_1,  \
                             g_0_xxyyyy_0_xxyyy_1,   \
                             g_0_xxyyyy_0_xxyyyy_0,  \
                             g_0_xxyyyy_0_xxyyyy_1,  \
                             g_0_xxyyyy_0_xxyyyz_0,  \
                             g_0_xxyyyy_0_xxyyyz_1,  \
                             g_0_xxyyyy_0_xxyyz_1,   \
                             g_0_xxyyyy_0_xxyyzz_0,  \
                             g_0_xxyyyy_0_xxyyzz_1,  \
                             g_0_xxyyyy_0_xxyzz_1,   \
                             g_0_xxyyyy_0_xxyzzz_0,  \
                             g_0_xxyyyy_0_xxyzzz_1,  \
                             g_0_xxyyyy_0_xyyyy_1,   \
                             g_0_xxyyyy_0_xyyyyy_0,  \
                             g_0_xxyyyy_0_xyyyyy_1,  \
                             g_0_xxyyyy_0_xyyyyz_0,  \
                             g_0_xxyyyy_0_xyyyyz_1,  \
                             g_0_xxyyyy_0_xyyyz_1,   \
                             g_0_xxyyyy_0_xyyyzz_0,  \
                             g_0_xxyyyy_0_xyyyzz_1,  \
                             g_0_xxyyyy_0_xyyzz_1,   \
                             g_0_xxyyyy_0_xyyzzz_0,  \
                             g_0_xxyyyy_0_xyyzzz_1,  \
                             g_0_xxyyyy_0_xyzzz_1,   \
                             g_0_xxyyyy_0_xyzzzz_0,  \
                             g_0_xxyyyy_0_xyzzzz_1,  \
                             g_0_xxyyyy_0_yyyyy_1,   \
                             g_0_xxyyyy_0_yyyyyy_0,  \
                             g_0_xxyyyy_0_yyyyyy_1,  \
                             g_0_xxyyyy_0_yyyyyz_0,  \
                             g_0_xxyyyy_0_yyyyyz_1,  \
                             g_0_xxyyyy_0_yyyyz_1,   \
                             g_0_xxyyyy_0_yyyyzz_0,  \
                             g_0_xxyyyy_0_yyyyzz_1,  \
                             g_0_xxyyyy_0_yyyzz_1,   \
                             g_0_xxyyyy_0_yyyzzz_0,  \
                             g_0_xxyyyy_0_yyyzzz_1,  \
                             g_0_xxyyyy_0_yyzzz_1,   \
                             g_0_xxyyyy_0_yyzzzz_0,  \
                             g_0_xxyyyy_0_yyzzzz_1,  \
                             g_0_xxyyyy_0_yzzzz_1,   \
                             g_0_xxyyyy_0_yzzzzz_0,  \
                             g_0_xxyyyy_0_yzzzzz_1,  \
                             g_0_xxyyyy_0_zzzzzz_0,  \
                             g_0_xxyyyy_0_zzzzzz_1,  \
                             g_0_xyyyy_0_xxxxxy_0,   \
                             g_0_xyyyy_0_xxxxxy_1,   \
                             g_0_xyyyy_0_xxxxyy_0,   \
                             g_0_xyyyy_0_xxxxyy_1,   \
                             g_0_xyyyy_0_xxxxyz_0,   \
                             g_0_xyyyy_0_xxxxyz_1,   \
                             g_0_xyyyy_0_xxxyyy_0,   \
                             g_0_xyyyy_0_xxxyyy_1,   \
                             g_0_xyyyy_0_xxxyyz_0,   \
                             g_0_xyyyy_0_xxxyyz_1,   \
                             g_0_xyyyy_0_xxxyzz_0,   \
                             g_0_xyyyy_0_xxxyzz_1,   \
                             g_0_xyyyy_0_xxyyyy_0,   \
                             g_0_xyyyy_0_xxyyyy_1,   \
                             g_0_xyyyy_0_xxyyyz_0,   \
                             g_0_xyyyy_0_xxyyyz_1,   \
                             g_0_xyyyy_0_xxyyzz_0,   \
                             g_0_xyyyy_0_xxyyzz_1,   \
                             g_0_xyyyy_0_xxyzzz_0,   \
                             g_0_xyyyy_0_xxyzzz_1,   \
                             g_0_xyyyy_0_xyyyyy_0,   \
                             g_0_xyyyy_0_xyyyyy_1,   \
                             g_0_xyyyy_0_xyyyyz_0,   \
                             g_0_xyyyy_0_xyyyyz_1,   \
                             g_0_xyyyy_0_xyyyzz_0,   \
                             g_0_xyyyy_0_xyyyzz_1,   \
                             g_0_xyyyy_0_xyyzzz_0,   \
                             g_0_xyyyy_0_xyyzzz_1,   \
                             g_0_xyyyy_0_xyzzzz_0,   \
                             g_0_xyyyy_0_xyzzzz_1,   \
                             g_0_xyyyy_0_yyyyyy_0,   \
                             g_0_xyyyy_0_yyyyyy_1,   \
                             g_0_xyyyy_0_yyyyyz_0,   \
                             g_0_xyyyy_0_yyyyyz_1,   \
                             g_0_xyyyy_0_yyyyzz_0,   \
                             g_0_xyyyy_0_yyyyzz_1,   \
                             g_0_xyyyy_0_yyyzzz_0,   \
                             g_0_xyyyy_0_yyyzzz_1,   \
                             g_0_xyyyy_0_yyzzzz_0,   \
                             g_0_xyyyy_0_yyzzzz_1,   \
                             g_0_xyyyy_0_yzzzzz_0,   \
                             g_0_xyyyy_0_yzzzzz_1,   \
                             g_0_xyyyy_0_zzzzzz_0,   \
                             g_0_xyyyy_0_zzzzzz_1,   \
                             wp_x,                   \
                             wp_y,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyy_0_xxxxxx_0[i] = 3.0 * g_0_xxxyy_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxxxx_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_xxxxxx_0[i] * pb_y + g_0_xxxyyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxxxxy_0[i] = 2.0 * g_0_xyyyy_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxxxy_1[i] * fti_ab_0 +
                                    5.0 * g_0_xxyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxxy_0[i] * pb_x + g_0_xxyyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxxxz_0[i] = 3.0 * g_0_xxxyy_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxxxz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_xxxxxz_0[i] * pb_y + g_0_xxxyyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxxxyy_0[i] = 2.0 * g_0_xyyyy_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxxyy_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxyy_0[i] * pb_x + g_0_xxyyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxxyz_0[i] = 2.0 * g_0_xyyyy_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxxyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxyz_0[i] * pb_x + g_0_xxyyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxxzz_0[i] = 3.0 * g_0_xxxyy_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxxzz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_xxxxzz_0[i] * pb_y + g_0_xxxyyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxxyyy_0[i] = 2.0 * g_0_xyyyy_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxyyy_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyyy_0[i] * pb_x + g_0_xxyyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxyyz_0[i] = 2.0 * g_0_xyyyy_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyyz_0[i] * pb_x + g_0_xxyyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxyzz_0[i] = 2.0 * g_0_xyyyy_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxxyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyzz_0[i] * pb_x + g_0_xxyyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxxzzz_0[i] = 3.0 * g_0_xxxyy_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxxzzz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_xxxzzz_0[i] * pb_y + g_0_xxxyyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xxyyyy_0[i] = 2.0 * g_0_xyyyy_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyyyy_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyyy_0[i] * pb_x + g_0_xxyyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxyyyz_0[i] = 2.0 * g_0_xyyyy_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyyz_0[i] * pb_x + g_0_xxyyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxyyzz_0[i] = 2.0 * g_0_xyyyy_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyzz_0[i] * pb_x + g_0_xxyyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxyzzz_0[i] = 2.0 * g_0_xyyyy_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xxyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyzzz_0[i] * pb_x + g_0_xxyyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xxzzzz_0[i] = 3.0 * g_0_xxxyy_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xxzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_xxzzzz_0[i] * pb_y + g_0_xxxyyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_xyyyyy_0[i] = 2.0 * g_0_xyyyy_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyyy_0[i] * pb_x + g_0_xxyyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyyyyz_0[i] = 2.0 * g_0_xyyyy_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyyz_0[i] * pb_x + g_0_xxyyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyyyzz_0[i] = 2.0 * g_0_xyyyy_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyzz_0[i] * pb_x + g_0_xxyyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyyzzz_0[i] = 2.0 * g_0_xyyyy_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyzzz_0[i] * pb_x + g_0_xxyyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xyzzzz_0[i] = 2.0 * g_0_xyyyy_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_xyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyzzzz_0[i] * pb_x + g_0_xxyyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_xzzzzz_0[i] = 3.0 * g_0_xxxyy_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxyy_0_xzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxxyyy_0_xzzzzz_0[i] * pb_y + g_0_xxxyyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxyyyy_0_yyyyyy_0[i] = 2.0 * g_0_xyyyy_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_yyyyyy_0[i] * pb_x + g_0_xxyyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyyyyz_0[i] = 2.0 * g_0_xyyyy_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_yyyyyz_0[i] * pb_x + g_0_xxyyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyyyzz_0[i] = 2.0 * g_0_xyyyy_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_yyyyzz_0[i] * pb_x + g_0_xxyyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyyzzz_0[i] = 2.0 * g_0_xyyyy_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_yyyzzz_0[i] * pb_x + g_0_xxyyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yyzzzz_0[i] = 2.0 * g_0_xyyyy_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_yyzzzz_0[i] * pb_x + g_0_xxyyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_yzzzzz_0[i] = 2.0 * g_0_xyyyy_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_yzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_yzzzzz_0[i] * pb_x + g_0_xxyyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxyyyy_0_zzzzzz_0[i] = 2.0 * g_0_xyyyy_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_zzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_zzzzzz_0[i] * pb_x + g_0_xxyyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 308-336 components of targeted buffer : SKSI

    auto g_0_xxxyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 308);

    auto g_0_xxxyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 309);

    auto g_0_xxxyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 310);

    auto g_0_xxxyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 311);

    auto g_0_xxxyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 312);

    auto g_0_xxxyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 313);

    auto g_0_xxxyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 314);

    auto g_0_xxxyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 315);

    auto g_0_xxxyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 316);

    auto g_0_xxxyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 317);

    auto g_0_xxxyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 318);

    auto g_0_xxxyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 319);

    auto g_0_xxxyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 320);

    auto g_0_xxxyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 321);

    auto g_0_xxxyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 322);

    auto g_0_xxxyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 323);

    auto g_0_xxxyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 324);

    auto g_0_xxxyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 325);

    auto g_0_xxxyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 326);

    auto g_0_xxxyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 327);

    auto g_0_xxxyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 328);

    auto g_0_xxxyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 329);

    auto g_0_xxxyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 330);

    auto g_0_xxxyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 331);

    auto g_0_xxxyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 332);

    auto g_0_xxxyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 333);

    auto g_0_xxxyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 334);

    auto g_0_xxxyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 335);

#pragma omp simd aligned(g_0_xxxyyy_0_xxxxx_1,       \
                             g_0_xxxyyy_0_xxxxxx_0,  \
                             g_0_xxxyyy_0_xxxxxx_1,  \
                             g_0_xxxyyy_0_xxxxxy_0,  \
                             g_0_xxxyyy_0_xxxxxy_1,  \
                             g_0_xxxyyy_0_xxxxxz_0,  \
                             g_0_xxxyyy_0_xxxxxz_1,  \
                             g_0_xxxyyy_0_xxxxy_1,   \
                             g_0_xxxyyy_0_xxxxyy_0,  \
                             g_0_xxxyyy_0_xxxxyy_1,  \
                             g_0_xxxyyy_0_xxxxyz_0,  \
                             g_0_xxxyyy_0_xxxxyz_1,  \
                             g_0_xxxyyy_0_xxxxz_1,   \
                             g_0_xxxyyy_0_xxxxzz_0,  \
                             g_0_xxxyyy_0_xxxxzz_1,  \
                             g_0_xxxyyy_0_xxxyy_1,   \
                             g_0_xxxyyy_0_xxxyyy_0,  \
                             g_0_xxxyyy_0_xxxyyy_1,  \
                             g_0_xxxyyy_0_xxxyyz_0,  \
                             g_0_xxxyyy_0_xxxyyz_1,  \
                             g_0_xxxyyy_0_xxxyz_1,   \
                             g_0_xxxyyy_0_xxxyzz_0,  \
                             g_0_xxxyyy_0_xxxyzz_1,  \
                             g_0_xxxyyy_0_xxxzz_1,   \
                             g_0_xxxyyy_0_xxxzzz_0,  \
                             g_0_xxxyyy_0_xxxzzz_1,  \
                             g_0_xxxyyy_0_xxyyy_1,   \
                             g_0_xxxyyy_0_xxyyyy_0,  \
                             g_0_xxxyyy_0_xxyyyy_1,  \
                             g_0_xxxyyy_0_xxyyyz_0,  \
                             g_0_xxxyyy_0_xxyyyz_1,  \
                             g_0_xxxyyy_0_xxyyz_1,   \
                             g_0_xxxyyy_0_xxyyzz_0,  \
                             g_0_xxxyyy_0_xxyyzz_1,  \
                             g_0_xxxyyy_0_xxyzz_1,   \
                             g_0_xxxyyy_0_xxyzzz_0,  \
                             g_0_xxxyyy_0_xxyzzz_1,  \
                             g_0_xxxyyy_0_xxzzz_1,   \
                             g_0_xxxyyy_0_xxzzzz_0,  \
                             g_0_xxxyyy_0_xxzzzz_1,  \
                             g_0_xxxyyy_0_xyyyy_1,   \
                             g_0_xxxyyy_0_xyyyyy_0,  \
                             g_0_xxxyyy_0_xyyyyy_1,  \
                             g_0_xxxyyy_0_xyyyyz_0,  \
                             g_0_xxxyyy_0_xyyyyz_1,  \
                             g_0_xxxyyy_0_xyyyz_1,   \
                             g_0_xxxyyy_0_xyyyzz_0,  \
                             g_0_xxxyyy_0_xyyyzz_1,  \
                             g_0_xxxyyy_0_xyyzz_1,   \
                             g_0_xxxyyy_0_xyyzzz_0,  \
                             g_0_xxxyyy_0_xyyzzz_1,  \
                             g_0_xxxyyy_0_xyzzz_1,   \
                             g_0_xxxyyy_0_xyzzzz_0,  \
                             g_0_xxxyyy_0_xyzzzz_1,  \
                             g_0_xxxyyy_0_xzzzz_1,   \
                             g_0_xxxyyy_0_xzzzzz_0,  \
                             g_0_xxxyyy_0_xzzzzz_1,  \
                             g_0_xxxyyy_0_yyyyy_1,   \
                             g_0_xxxyyy_0_yyyyyy_0,  \
                             g_0_xxxyyy_0_yyyyyy_1,  \
                             g_0_xxxyyy_0_yyyyyz_0,  \
                             g_0_xxxyyy_0_yyyyyz_1,  \
                             g_0_xxxyyy_0_yyyyz_1,   \
                             g_0_xxxyyy_0_yyyyzz_0,  \
                             g_0_xxxyyy_0_yyyyzz_1,  \
                             g_0_xxxyyy_0_yyyzz_1,   \
                             g_0_xxxyyy_0_yyyzzz_0,  \
                             g_0_xxxyyy_0_yyyzzz_1,  \
                             g_0_xxxyyy_0_yyzzz_1,   \
                             g_0_xxxyyy_0_yyzzzz_0,  \
                             g_0_xxxyyy_0_yyzzzz_1,  \
                             g_0_xxxyyy_0_yzzzz_1,   \
                             g_0_xxxyyy_0_yzzzzz_0,  \
                             g_0_xxxyyy_0_yzzzzz_1,  \
                             g_0_xxxyyy_0_zzzzz_1,   \
                             g_0_xxxyyy_0_zzzzzz_0,  \
                             g_0_xxxyyy_0_zzzzzz_1,  \
                             g_0_xxxyyyz_0_xxxxxx_0, \
                             g_0_xxxyyyz_0_xxxxxy_0, \
                             g_0_xxxyyyz_0_xxxxxz_0, \
                             g_0_xxxyyyz_0_xxxxyy_0, \
                             g_0_xxxyyyz_0_xxxxyz_0, \
                             g_0_xxxyyyz_0_xxxxzz_0, \
                             g_0_xxxyyyz_0_xxxyyy_0, \
                             g_0_xxxyyyz_0_xxxyyz_0, \
                             g_0_xxxyyyz_0_xxxyzz_0, \
                             g_0_xxxyyyz_0_xxxzzz_0, \
                             g_0_xxxyyyz_0_xxyyyy_0, \
                             g_0_xxxyyyz_0_xxyyyz_0, \
                             g_0_xxxyyyz_0_xxyyzz_0, \
                             g_0_xxxyyyz_0_xxyzzz_0, \
                             g_0_xxxyyyz_0_xxzzzz_0, \
                             g_0_xxxyyyz_0_xyyyyy_0, \
                             g_0_xxxyyyz_0_xyyyyz_0, \
                             g_0_xxxyyyz_0_xyyyzz_0, \
                             g_0_xxxyyyz_0_xyyzzz_0, \
                             g_0_xxxyyyz_0_xyzzzz_0, \
                             g_0_xxxyyyz_0_xzzzzz_0, \
                             g_0_xxxyyyz_0_yyyyyy_0, \
                             g_0_xxxyyyz_0_yyyyyz_0, \
                             g_0_xxxyyyz_0_yyyyzz_0, \
                             g_0_xxxyyyz_0_yyyzzz_0, \
                             g_0_xxxyyyz_0_yyzzzz_0, \
                             g_0_xxxyyyz_0_yzzzzz_0, \
                             g_0_xxxyyyz_0_zzzzzz_0, \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyz_0_xxxxxx_0[i] = g_0_xxxyyy_0_xxxxxx_0[i] * pb_z + g_0_xxxyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxxy_0[i] = g_0_xxxyyy_0_xxxxxy_0[i] * pb_z + g_0_xxxyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxxz_0[i] = g_0_xxxyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxxz_0[i] * pb_z + g_0_xxxyyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxyy_0[i] = g_0_xxxyyy_0_xxxxyy_0[i] * pb_z + g_0_xxxyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxyz_0[i] = g_0_xxxyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxyz_0[i] * pb_z + g_0_xxxyyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxxzz_0[i] = 2.0 * g_0_xxxyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxxzz_0[i] * pb_z + g_0_xxxyyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxyyy_0[i] = g_0_xxxyyy_0_xxxyyy_0[i] * pb_z + g_0_xxxyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxyyz_0[i] = g_0_xxxyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyyz_0[i] * pb_z + g_0_xxxyyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxyzz_0[i] = 2.0 * g_0_xxxyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxyzz_0[i] * pb_z + g_0_xxxyyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxxzzz_0[i] = 3.0 * g_0_xxxyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxxzzz_0[i] * pb_z + g_0_xxxyyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyyyy_0[i] = g_0_xxxyyy_0_xxyyyy_0[i] * pb_z + g_0_xxxyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyyyz_0[i] = g_0_xxxyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyyz_0[i] * pb_z + g_0_xxxyyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyyzz_0[i] = 2.0 * g_0_xxxyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyyzz_0[i] * pb_z + g_0_xxxyyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxyzzz_0[i] = 3.0 * g_0_xxxyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxyzzz_0[i] * pb_z + g_0_xxxyyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xxzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xxzzzz_0[i] * pb_z + g_0_xxxyyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyyyy_0[i] = g_0_xxxyyy_0_xyyyyy_0[i] * pb_z + g_0_xxxyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyyyz_0[i] = g_0_xxxyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyyz_0[i] * pb_z + g_0_xxxyyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyyzz_0[i] = 2.0 * g_0_xxxyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyyzz_0[i] * pb_z + g_0_xxxyyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyyzzz_0[i] = 3.0 * g_0_xxxyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyyzzz_0[i] * pb_z + g_0_xxxyyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xyzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xyzzzz_0[i] * pb_z + g_0_xxxyyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_xzzzzz_0[i] = 5.0 * g_0_xxxyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_xzzzzz_0[i] * pb_z + g_0_xxxyyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyyyy_0[i] = g_0_xxxyyy_0_yyyyyy_0[i] * pb_z + g_0_xxxyyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyyyz_0[i] = g_0_xxxyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyyyyz_0[i] * pb_z + g_0_xxxyyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyyzz_0[i] = 2.0 * g_0_xxxyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyyyzz_0[i] * pb_z + g_0_xxxyyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyyzzz_0[i] = 3.0 * g_0_xxxyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyyzzz_0[i] * pb_z + g_0_xxxyyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yyzzzz_0[i] = 4.0 * g_0_xxxyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yyzzzz_0[i] * pb_z + g_0_xxxyyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_yzzzzz_0[i] = 5.0 * g_0_xxxyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_yzzzzz_0[i] * pb_z + g_0_xxxyyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyz_0_zzzzzz_0[i] = 6.0 * g_0_xxxyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxyyy_0_zzzzzz_0[i] * pb_z + g_0_xxxyyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 336-364 components of targeted buffer : SKSI

    auto g_0_xxxyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 336);

    auto g_0_xxxyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 337);

    auto g_0_xxxyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 338);

    auto g_0_xxxyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 339);

    auto g_0_xxxyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 340);

    auto g_0_xxxyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 341);

    auto g_0_xxxyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 342);

    auto g_0_xxxyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 343);

    auto g_0_xxxyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 344);

    auto g_0_xxxyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 345);

    auto g_0_xxxyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 346);

    auto g_0_xxxyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 347);

    auto g_0_xxxyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 348);

    auto g_0_xxxyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 349);

    auto g_0_xxxyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 350);

    auto g_0_xxxyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 351);

    auto g_0_xxxyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 352);

    auto g_0_xxxyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 353);

    auto g_0_xxxyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 354);

    auto g_0_xxxyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 355);

    auto g_0_xxxyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 356);

    auto g_0_xxxyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 357);

    auto g_0_xxxyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 358);

    auto g_0_xxxyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 359);

    auto g_0_xxxyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 360);

    auto g_0_xxxyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 361);

    auto g_0_xxxyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 362);

    auto g_0_xxxyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 363);

#pragma omp simd aligned(g_0_xxxyy_0_xxxxxy_0,       \
                             g_0_xxxyy_0_xxxxxy_1,   \
                             g_0_xxxyy_0_xxxxyy_0,   \
                             g_0_xxxyy_0_xxxxyy_1,   \
                             g_0_xxxyy_0_xxxyyy_0,   \
                             g_0_xxxyy_0_xxxyyy_1,   \
                             g_0_xxxyy_0_xxyyyy_0,   \
                             g_0_xxxyy_0_xxyyyy_1,   \
                             g_0_xxxyy_0_xyyyyy_0,   \
                             g_0_xxxyy_0_xyyyyy_1,   \
                             g_0_xxxyyz_0_xxxxxy_0,  \
                             g_0_xxxyyz_0_xxxxxy_1,  \
                             g_0_xxxyyz_0_xxxxyy_0,  \
                             g_0_xxxyyz_0_xxxxyy_1,  \
                             g_0_xxxyyz_0_xxxyyy_0,  \
                             g_0_xxxyyz_0_xxxyyy_1,  \
                             g_0_xxxyyz_0_xxyyyy_0,  \
                             g_0_xxxyyz_0_xxyyyy_1,  \
                             g_0_xxxyyz_0_xyyyyy_0,  \
                             g_0_xxxyyz_0_xyyyyy_1,  \
                             g_0_xxxyyzz_0_xxxxxx_0, \
                             g_0_xxxyyzz_0_xxxxxy_0, \
                             g_0_xxxyyzz_0_xxxxxz_0, \
                             g_0_xxxyyzz_0_xxxxyy_0, \
                             g_0_xxxyyzz_0_xxxxyz_0, \
                             g_0_xxxyyzz_0_xxxxzz_0, \
                             g_0_xxxyyzz_0_xxxyyy_0, \
                             g_0_xxxyyzz_0_xxxyyz_0, \
                             g_0_xxxyyzz_0_xxxyzz_0, \
                             g_0_xxxyyzz_0_xxxzzz_0, \
                             g_0_xxxyyzz_0_xxyyyy_0, \
                             g_0_xxxyyzz_0_xxyyyz_0, \
                             g_0_xxxyyzz_0_xxyyzz_0, \
                             g_0_xxxyyzz_0_xxyzzz_0, \
                             g_0_xxxyyzz_0_xxzzzz_0, \
                             g_0_xxxyyzz_0_xyyyyy_0, \
                             g_0_xxxyyzz_0_xyyyyz_0, \
                             g_0_xxxyyzz_0_xyyyzz_0, \
                             g_0_xxxyyzz_0_xyyzzz_0, \
                             g_0_xxxyyzz_0_xyzzzz_0, \
                             g_0_xxxyyzz_0_xzzzzz_0, \
                             g_0_xxxyyzz_0_yyyyyy_0, \
                             g_0_xxxyyzz_0_yyyyyz_0, \
                             g_0_xxxyyzz_0_yyyyzz_0, \
                             g_0_xxxyyzz_0_yyyzzz_0, \
                             g_0_xxxyyzz_0_yyzzzz_0, \
                             g_0_xxxyyzz_0_yzzzzz_0, \
                             g_0_xxxyyzz_0_zzzzzz_0, \
                             g_0_xxxyzz_0_xxxxxx_0,  \
                             g_0_xxxyzz_0_xxxxxx_1,  \
                             g_0_xxxyzz_0_xxxxxz_0,  \
                             g_0_xxxyzz_0_xxxxxz_1,  \
                             g_0_xxxyzz_0_xxxxzz_0,  \
                             g_0_xxxyzz_0_xxxxzz_1,  \
                             g_0_xxxyzz_0_xxxzzz_0,  \
                             g_0_xxxyzz_0_xxxzzz_1,  \
                             g_0_xxxyzz_0_xxzzzz_0,  \
                             g_0_xxxyzz_0_xxzzzz_1,  \
                             g_0_xxxyzz_0_xzzzzz_0,  \
                             g_0_xxxyzz_0_xzzzzz_1,  \
                             g_0_xxxzz_0_xxxxxx_0,   \
                             g_0_xxxzz_0_xxxxxx_1,   \
                             g_0_xxxzz_0_xxxxxz_0,   \
                             g_0_xxxzz_0_xxxxxz_1,   \
                             g_0_xxxzz_0_xxxxzz_0,   \
                             g_0_xxxzz_0_xxxxzz_1,   \
                             g_0_xxxzz_0_xxxzzz_0,   \
                             g_0_xxxzz_0_xxxzzz_1,   \
                             g_0_xxxzz_0_xxzzzz_0,   \
                             g_0_xxxzz_0_xxzzzz_1,   \
                             g_0_xxxzz_0_xzzzzz_0,   \
                             g_0_xxxzz_0_xzzzzz_1,   \
                             g_0_xxyyzz_0_xxxxyz_0,  \
                             g_0_xxyyzz_0_xxxxyz_1,  \
                             g_0_xxyyzz_0_xxxyyz_0,  \
                             g_0_xxyyzz_0_xxxyyz_1,  \
                             g_0_xxyyzz_0_xxxyz_1,   \
                             g_0_xxyyzz_0_xxxyzz_0,  \
                             g_0_xxyyzz_0_xxxyzz_1,  \
                             g_0_xxyyzz_0_xxyyyz_0,  \
                             g_0_xxyyzz_0_xxyyyz_1,  \
                             g_0_xxyyzz_0_xxyyz_1,   \
                             g_0_xxyyzz_0_xxyyzz_0,  \
                             g_0_xxyyzz_0_xxyyzz_1,  \
                             g_0_xxyyzz_0_xxyzz_1,   \
                             g_0_xxyyzz_0_xxyzzz_0,  \
                             g_0_xxyyzz_0_xxyzzz_1,  \
                             g_0_xxyyzz_0_xyyyyz_0,  \
                             g_0_xxyyzz_0_xyyyyz_1,  \
                             g_0_xxyyzz_0_xyyyz_1,   \
                             g_0_xxyyzz_0_xyyyzz_0,  \
                             g_0_xxyyzz_0_xyyyzz_1,  \
                             g_0_xxyyzz_0_xyyzz_1,   \
                             g_0_xxyyzz_0_xyyzzz_0,  \
                             g_0_xxyyzz_0_xyyzzz_1,  \
                             g_0_xxyyzz_0_xyzzz_1,   \
                             g_0_xxyyzz_0_xyzzzz_0,  \
                             g_0_xxyyzz_0_xyzzzz_1,  \
                             g_0_xxyyzz_0_yyyyyy_0,  \
                             g_0_xxyyzz_0_yyyyyy_1,  \
                             g_0_xxyyzz_0_yyyyyz_0,  \
                             g_0_xxyyzz_0_yyyyyz_1,  \
                             g_0_xxyyzz_0_yyyyz_1,   \
                             g_0_xxyyzz_0_yyyyzz_0,  \
                             g_0_xxyyzz_0_yyyyzz_1,  \
                             g_0_xxyyzz_0_yyyzz_1,   \
                             g_0_xxyyzz_0_yyyzzz_0,  \
                             g_0_xxyyzz_0_yyyzzz_1,  \
                             g_0_xxyyzz_0_yyzzz_1,   \
                             g_0_xxyyzz_0_yyzzzz_0,  \
                             g_0_xxyyzz_0_yyzzzz_1,  \
                             g_0_xxyyzz_0_yzzzz_1,   \
                             g_0_xxyyzz_0_yzzzzz_0,  \
                             g_0_xxyyzz_0_yzzzzz_1,  \
                             g_0_xxyyzz_0_zzzzzz_0,  \
                             g_0_xxyyzz_0_zzzzzz_1,  \
                             g_0_xyyzz_0_xxxxyz_0,   \
                             g_0_xyyzz_0_xxxxyz_1,   \
                             g_0_xyyzz_0_xxxyyz_0,   \
                             g_0_xyyzz_0_xxxyyz_1,   \
                             g_0_xyyzz_0_xxxyzz_0,   \
                             g_0_xyyzz_0_xxxyzz_1,   \
                             g_0_xyyzz_0_xxyyyz_0,   \
                             g_0_xyyzz_0_xxyyyz_1,   \
                             g_0_xyyzz_0_xxyyzz_0,   \
                             g_0_xyyzz_0_xxyyzz_1,   \
                             g_0_xyyzz_0_xxyzzz_0,   \
                             g_0_xyyzz_0_xxyzzz_1,   \
                             g_0_xyyzz_0_xyyyyz_0,   \
                             g_0_xyyzz_0_xyyyyz_1,   \
                             g_0_xyyzz_0_xyyyzz_0,   \
                             g_0_xyyzz_0_xyyyzz_1,   \
                             g_0_xyyzz_0_xyyzzz_0,   \
                             g_0_xyyzz_0_xyyzzz_1,   \
                             g_0_xyyzz_0_xyzzzz_0,   \
                             g_0_xyyzz_0_xyzzzz_1,   \
                             g_0_xyyzz_0_yyyyyy_0,   \
                             g_0_xyyzz_0_yyyyyy_1,   \
                             g_0_xyyzz_0_yyyyyz_0,   \
                             g_0_xyyzz_0_yyyyyz_1,   \
                             g_0_xyyzz_0_yyyyzz_0,   \
                             g_0_xyyzz_0_yyyyzz_1,   \
                             g_0_xyyzz_0_yyyzzz_0,   \
                             g_0_xyyzz_0_yyyzzz_1,   \
                             g_0_xyyzz_0_yyzzzz_0,   \
                             g_0_xyyzz_0_yyzzzz_1,   \
                             g_0_xyyzz_0_yzzzzz_0,   \
                             g_0_xyyzz_0_yzzzzz_1,   \
                             g_0_xyyzz_0_zzzzzz_0,   \
                             g_0_xyyzz_0_zzzzzz_1,   \
                             wp_x,                   \
                             wp_y,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyzz_0_xxxxxx_0[i] = g_0_xxxzz_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxxxx_0[i] * pb_y +
                                    g_0_xxxyzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxxxxy_0[i] = g_0_xxxyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxxxxy_0[i] * pb_z +
                                    g_0_xxxyyz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxxxxz_0[i] = g_0_xxxzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxxxz_0[i] * pb_y +
                                    g_0_xxxyzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxxxyy_0[i] = g_0_xxxyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxxxyy_0[i] * pb_z +
                                    g_0_xxxyyz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxxxyz_0[i] = 2.0 * g_0_xyyzz_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxxxyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxyyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxxxyz_0[i] * pb_x + g_0_xxyyzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxxxzz_0[i] = g_0_xxxzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxxzz_0[i] * pb_y +
                                    g_0_xxxyzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxxyyy_0[i] = g_0_xxxyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxxyyy_0[i] * pb_z +
                                    g_0_xxxyyz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxxyyz_0[i] = 2.0 * g_0_xyyzz_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxxyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxyyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxxyyz_0[i] * pb_x + g_0_xxyyzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxxyzz_0[i] = 2.0 * g_0_xyyzz_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxxyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxyyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxxyzz_0[i] * pb_x + g_0_xxyyzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxxzzz_0[i] = g_0_xxxzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxxzzz_0[i] * pb_y +
                                    g_0_xxxyzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xxyyyy_0[i] = g_0_xxxyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xxyyyy_0[i] * pb_z +
                                    g_0_xxxyyz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xxyyyz_0[i] = 2.0 * g_0_xyyzz_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxyyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxyyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxyyyz_0[i] * pb_x + g_0_xxyyzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxyyzz_0[i] = 2.0 * g_0_xyyzz_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxyyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxyyzz_0[i] * pb_x + g_0_xxyyzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxyzzz_0[i] = 2.0 * g_0_xyyzz_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xxyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxyyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xxyzzz_0[i] * pb_x + g_0_xxyyzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xxzzzz_0[i] = g_0_xxxzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xxzzzz_0[i] * pb_y +
                                    g_0_xxxyzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_xyyyyy_0[i] = g_0_xxxyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_xxxyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxxyyz_0_xyyyyy_0[i] * pb_z +
                                    g_0_xxxyyz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxyyzz_0_xyyyyz_0[i] = 2.0 * g_0_xyyzz_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xyyyyz_0[i] * pb_x + g_0_xxyyzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xyyyzz_0[i] = 2.0 * g_0_xyyzz_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xyyyzz_0[i] * pb_x + g_0_xxyyzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xyyzzz_0[i] = 2.0 * g_0_xyyzz_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xyyzzz_0[i] * pb_x + g_0_xxyyzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xyzzzz_0[i] = 2.0 * g_0_xyyzz_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_xyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxyyzz_0_xyzzzz_0[i] * pb_x + g_0_xxyyzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_xzzzzz_0[i] = g_0_xxxzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_xxxzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxxyzz_0_xzzzzz_0[i] * pb_y +
                                    g_0_xxxyzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxyyzz_0_yyyyyy_0[i] = 2.0 * g_0_xyyzz_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_yyyyyy_0[i] * pb_x + g_0_xxyyzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyyyyz_0[i] = 2.0 * g_0_xyyzz_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_yyyyyz_0[i] * pb_x + g_0_xxyyzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyyyzz_0[i] = 2.0 * g_0_xyyzz_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_yyyyzz_0[i] * pb_x + g_0_xxyyzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyyzzz_0[i] = 2.0 * g_0_xyyzz_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_yyyzzz_0[i] * pb_x + g_0_xxyyzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yyzzzz_0[i] = 2.0 * g_0_xyyzz_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_yyzzzz_0[i] * pb_x + g_0_xxyyzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_yzzzzz_0[i] = 2.0 * g_0_xyyzz_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_yzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_yzzzzz_0[i] * pb_x + g_0_xxyyzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxyyzz_0_zzzzzz_0[i] = 2.0 * g_0_xyyzz_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_zzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_zzzzzz_0[i] * pb_x + g_0_xxyyzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 364-392 components of targeted buffer : SKSI

    auto g_0_xxxyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 364);

    auto g_0_xxxyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 365);

    auto g_0_xxxyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 366);

    auto g_0_xxxyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 367);

    auto g_0_xxxyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 368);

    auto g_0_xxxyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 369);

    auto g_0_xxxyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 370);

    auto g_0_xxxyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 371);

    auto g_0_xxxyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 372);

    auto g_0_xxxyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 373);

    auto g_0_xxxyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 374);

    auto g_0_xxxyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 375);

    auto g_0_xxxyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 376);

    auto g_0_xxxyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 377);

    auto g_0_xxxyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 378);

    auto g_0_xxxyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 379);

    auto g_0_xxxyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 380);

    auto g_0_xxxyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 381);

    auto g_0_xxxyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 382);

    auto g_0_xxxyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 383);

    auto g_0_xxxyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 384);

    auto g_0_xxxyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 385);

    auto g_0_xxxyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 386);

    auto g_0_xxxyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 387);

    auto g_0_xxxyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 388);

    auto g_0_xxxyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 389);

    auto g_0_xxxyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 390);

    auto g_0_xxxyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 391);

#pragma omp simd aligned(g_0_xxxyzzz_0_xxxxxx_0,     \
                             g_0_xxxyzzz_0_xxxxxy_0, \
                             g_0_xxxyzzz_0_xxxxxz_0, \
                             g_0_xxxyzzz_0_xxxxyy_0, \
                             g_0_xxxyzzz_0_xxxxyz_0, \
                             g_0_xxxyzzz_0_xxxxzz_0, \
                             g_0_xxxyzzz_0_xxxyyy_0, \
                             g_0_xxxyzzz_0_xxxyyz_0, \
                             g_0_xxxyzzz_0_xxxyzz_0, \
                             g_0_xxxyzzz_0_xxxzzz_0, \
                             g_0_xxxyzzz_0_xxyyyy_0, \
                             g_0_xxxyzzz_0_xxyyyz_0, \
                             g_0_xxxyzzz_0_xxyyzz_0, \
                             g_0_xxxyzzz_0_xxyzzz_0, \
                             g_0_xxxyzzz_0_xxzzzz_0, \
                             g_0_xxxyzzz_0_xyyyyy_0, \
                             g_0_xxxyzzz_0_xyyyyz_0, \
                             g_0_xxxyzzz_0_xyyyzz_0, \
                             g_0_xxxyzzz_0_xyyzzz_0, \
                             g_0_xxxyzzz_0_xyzzzz_0, \
                             g_0_xxxyzzz_0_xzzzzz_0, \
                             g_0_xxxyzzz_0_yyyyyy_0, \
                             g_0_xxxyzzz_0_yyyyyz_0, \
                             g_0_xxxyzzz_0_yyyyzz_0, \
                             g_0_xxxyzzz_0_yyyzzz_0, \
                             g_0_xxxyzzz_0_yyzzzz_0, \
                             g_0_xxxyzzz_0_yzzzzz_0, \
                             g_0_xxxyzzz_0_zzzzzz_0, \
                             g_0_xxxzzz_0_xxxxx_1,   \
                             g_0_xxxzzz_0_xxxxxx_0,  \
                             g_0_xxxzzz_0_xxxxxx_1,  \
                             g_0_xxxzzz_0_xxxxxy_0,  \
                             g_0_xxxzzz_0_xxxxxy_1,  \
                             g_0_xxxzzz_0_xxxxxz_0,  \
                             g_0_xxxzzz_0_xxxxxz_1,  \
                             g_0_xxxzzz_0_xxxxy_1,   \
                             g_0_xxxzzz_0_xxxxyy_0,  \
                             g_0_xxxzzz_0_xxxxyy_1,  \
                             g_0_xxxzzz_0_xxxxyz_0,  \
                             g_0_xxxzzz_0_xxxxyz_1,  \
                             g_0_xxxzzz_0_xxxxz_1,   \
                             g_0_xxxzzz_0_xxxxzz_0,  \
                             g_0_xxxzzz_0_xxxxzz_1,  \
                             g_0_xxxzzz_0_xxxyy_1,   \
                             g_0_xxxzzz_0_xxxyyy_0,  \
                             g_0_xxxzzz_0_xxxyyy_1,  \
                             g_0_xxxzzz_0_xxxyyz_0,  \
                             g_0_xxxzzz_0_xxxyyz_1,  \
                             g_0_xxxzzz_0_xxxyz_1,   \
                             g_0_xxxzzz_0_xxxyzz_0,  \
                             g_0_xxxzzz_0_xxxyzz_1,  \
                             g_0_xxxzzz_0_xxxzz_1,   \
                             g_0_xxxzzz_0_xxxzzz_0,  \
                             g_0_xxxzzz_0_xxxzzz_1,  \
                             g_0_xxxzzz_0_xxyyy_1,   \
                             g_0_xxxzzz_0_xxyyyy_0,  \
                             g_0_xxxzzz_0_xxyyyy_1,  \
                             g_0_xxxzzz_0_xxyyyz_0,  \
                             g_0_xxxzzz_0_xxyyyz_1,  \
                             g_0_xxxzzz_0_xxyyz_1,   \
                             g_0_xxxzzz_0_xxyyzz_0,  \
                             g_0_xxxzzz_0_xxyyzz_1,  \
                             g_0_xxxzzz_0_xxyzz_1,   \
                             g_0_xxxzzz_0_xxyzzz_0,  \
                             g_0_xxxzzz_0_xxyzzz_1,  \
                             g_0_xxxzzz_0_xxzzz_1,   \
                             g_0_xxxzzz_0_xxzzzz_0,  \
                             g_0_xxxzzz_0_xxzzzz_1,  \
                             g_0_xxxzzz_0_xyyyy_1,   \
                             g_0_xxxzzz_0_xyyyyy_0,  \
                             g_0_xxxzzz_0_xyyyyy_1,  \
                             g_0_xxxzzz_0_xyyyyz_0,  \
                             g_0_xxxzzz_0_xyyyyz_1,  \
                             g_0_xxxzzz_0_xyyyz_1,   \
                             g_0_xxxzzz_0_xyyyzz_0,  \
                             g_0_xxxzzz_0_xyyyzz_1,  \
                             g_0_xxxzzz_0_xyyzz_1,   \
                             g_0_xxxzzz_0_xyyzzz_0,  \
                             g_0_xxxzzz_0_xyyzzz_1,  \
                             g_0_xxxzzz_0_xyzzz_1,   \
                             g_0_xxxzzz_0_xyzzzz_0,  \
                             g_0_xxxzzz_0_xyzzzz_1,  \
                             g_0_xxxzzz_0_xzzzz_1,   \
                             g_0_xxxzzz_0_xzzzzz_0,  \
                             g_0_xxxzzz_0_xzzzzz_1,  \
                             g_0_xxxzzz_0_yyyyy_1,   \
                             g_0_xxxzzz_0_yyyyyy_0,  \
                             g_0_xxxzzz_0_yyyyyy_1,  \
                             g_0_xxxzzz_0_yyyyyz_0,  \
                             g_0_xxxzzz_0_yyyyyz_1,  \
                             g_0_xxxzzz_0_yyyyz_1,   \
                             g_0_xxxzzz_0_yyyyzz_0,  \
                             g_0_xxxzzz_0_yyyyzz_1,  \
                             g_0_xxxzzz_0_yyyzz_1,   \
                             g_0_xxxzzz_0_yyyzzz_0,  \
                             g_0_xxxzzz_0_yyyzzz_1,  \
                             g_0_xxxzzz_0_yyzzz_1,   \
                             g_0_xxxzzz_0_yyzzzz_0,  \
                             g_0_xxxzzz_0_yyzzzz_1,  \
                             g_0_xxxzzz_0_yzzzz_1,   \
                             g_0_xxxzzz_0_yzzzzz_0,  \
                             g_0_xxxzzz_0_yzzzzz_1,  \
                             g_0_xxxzzz_0_zzzzz_1,   \
                             g_0_xxxzzz_0_zzzzzz_0,  \
                             g_0_xxxzzz_0_zzzzzz_1,  \
                             wp_y,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzz_0_xxxxxx_0[i] = g_0_xxxzzz_0_xxxxxx_0[i] * pb_y + g_0_xxxzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxxy_0[i] = g_0_xxxzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxxy_0[i] * pb_y + g_0_xxxzzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxxz_0[i] = g_0_xxxzzz_0_xxxxxz_0[i] * pb_y + g_0_xxxzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxyy_0[i] = 2.0 * g_0_xxxzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxyy_0[i] * pb_y + g_0_xxxzzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxyz_0[i] = g_0_xxxzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxxyz_0[i] * pb_y + g_0_xxxzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxxzz_0[i] = g_0_xxxzzz_0_xxxxzz_0[i] * pb_y + g_0_xxxzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxyyy_0[i] = 3.0 * g_0_xxxzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyyy_0[i] * pb_y + g_0_xxxzzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxyyz_0[i] = 2.0 * g_0_xxxzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyyz_0[i] * pb_y + g_0_xxxzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxyzz_0[i] = g_0_xxxzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxxyzz_0[i] * pb_y + g_0_xxxzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxxzzz_0[i] = g_0_xxxzzz_0_xxxzzz_0[i] * pb_y + g_0_xxxzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyyyy_0[i] = 4.0 * g_0_xxxzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyyy_0[i] * pb_y + g_0_xxxzzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyyyz_0[i] = 3.0 * g_0_xxxzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyyz_0[i] * pb_y + g_0_xxxzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyyzz_0[i] = 2.0 * g_0_xxxzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyyzz_0[i] * pb_y + g_0_xxxzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxyzzz_0[i] = g_0_xxxzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xxyzzz_0[i] * pb_y + g_0_xxxzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xxzzzz_0[i] = g_0_xxxzzz_0_xxzzzz_0[i] * pb_y + g_0_xxxzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyyyy_0[i] = 5.0 * g_0_xxxzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyyy_0[i] * pb_y + g_0_xxxzzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyyyz_0[i] = 4.0 * g_0_xxxzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyyz_0[i] * pb_y + g_0_xxxzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyyzz_0[i] = 3.0 * g_0_xxxzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyyzz_0[i] * pb_y + g_0_xxxzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyyzzz_0[i] = 2.0 * g_0_xxxzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyyzzz_0[i] * pb_y + g_0_xxxzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xyzzzz_0[i] = g_0_xxxzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_xyzzzz_0[i] * pb_y + g_0_xxxzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_xzzzzz_0[i] = g_0_xxxzzz_0_xzzzzz_0[i] * pb_y + g_0_xxxzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyyyy_0[i] = 6.0 * g_0_xxxzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyyyy_0[i] * pb_y + g_0_xxxzzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyyyz_0[i] = 5.0 * g_0_xxxzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyyyz_0[i] * pb_y + g_0_xxxzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyyzz_0[i] = 4.0 * g_0_xxxzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyyzz_0[i] * pb_y + g_0_xxxzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyyzzz_0[i] = 3.0 * g_0_xxxzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyyzzz_0[i] * pb_y + g_0_xxxzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yyzzzz_0[i] = 2.0 * g_0_xxxzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yyzzzz_0[i] * pb_y + g_0_xxxzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_yzzzzz_0[i] = g_0_xxxzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxzzz_0_yzzzzz_0[i] * pb_y + g_0_xxxzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzz_0_zzzzzz_0[i] = g_0_xxxzzz_0_zzzzzz_0[i] * pb_y + g_0_xxxzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 392-420 components of targeted buffer : SKSI

    auto g_0_xxxzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 392);

    auto g_0_xxxzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 393);

    auto g_0_xxxzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 394);

    auto g_0_xxxzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 395);

    auto g_0_xxxzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 396);

    auto g_0_xxxzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 397);

    auto g_0_xxxzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 398);

    auto g_0_xxxzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 399);

    auto g_0_xxxzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 400);

    auto g_0_xxxzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 401);

    auto g_0_xxxzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 402);

    auto g_0_xxxzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 403);

    auto g_0_xxxzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 404);

    auto g_0_xxxzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 405);

    auto g_0_xxxzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 406);

    auto g_0_xxxzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 407);

    auto g_0_xxxzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 408);

    auto g_0_xxxzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 409);

    auto g_0_xxxzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 410);

    auto g_0_xxxzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 411);

    auto g_0_xxxzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 412);

    auto g_0_xxxzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 413);

    auto g_0_xxxzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 414);

    auto g_0_xxxzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 415);

    auto g_0_xxxzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 416);

    auto g_0_xxxzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 417);

    auto g_0_xxxzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 418);

    auto g_0_xxxzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 419);

#pragma omp simd aligned(g_0_xxxzz_0_xxxxxx_0,       \
                             g_0_xxxzz_0_xxxxxx_1,   \
                             g_0_xxxzz_0_xxxxxy_0,   \
                             g_0_xxxzz_0_xxxxxy_1,   \
                             g_0_xxxzz_0_xxxxyy_0,   \
                             g_0_xxxzz_0_xxxxyy_1,   \
                             g_0_xxxzz_0_xxxyyy_0,   \
                             g_0_xxxzz_0_xxxyyy_1,   \
                             g_0_xxxzz_0_xxyyyy_0,   \
                             g_0_xxxzz_0_xxyyyy_1,   \
                             g_0_xxxzz_0_xyyyyy_0,   \
                             g_0_xxxzz_0_xyyyyy_1,   \
                             g_0_xxxzzz_0_xxxxxx_0,  \
                             g_0_xxxzzz_0_xxxxxx_1,  \
                             g_0_xxxzzz_0_xxxxxy_0,  \
                             g_0_xxxzzz_0_xxxxxy_1,  \
                             g_0_xxxzzz_0_xxxxyy_0,  \
                             g_0_xxxzzz_0_xxxxyy_1,  \
                             g_0_xxxzzz_0_xxxyyy_0,  \
                             g_0_xxxzzz_0_xxxyyy_1,  \
                             g_0_xxxzzz_0_xxyyyy_0,  \
                             g_0_xxxzzz_0_xxyyyy_1,  \
                             g_0_xxxzzz_0_xyyyyy_0,  \
                             g_0_xxxzzz_0_xyyyyy_1,  \
                             g_0_xxxzzzz_0_xxxxxx_0, \
                             g_0_xxxzzzz_0_xxxxxy_0, \
                             g_0_xxxzzzz_0_xxxxxz_0, \
                             g_0_xxxzzzz_0_xxxxyy_0, \
                             g_0_xxxzzzz_0_xxxxyz_0, \
                             g_0_xxxzzzz_0_xxxxzz_0, \
                             g_0_xxxzzzz_0_xxxyyy_0, \
                             g_0_xxxzzzz_0_xxxyyz_0, \
                             g_0_xxxzzzz_0_xxxyzz_0, \
                             g_0_xxxzzzz_0_xxxzzz_0, \
                             g_0_xxxzzzz_0_xxyyyy_0, \
                             g_0_xxxzzzz_0_xxyyyz_0, \
                             g_0_xxxzzzz_0_xxyyzz_0, \
                             g_0_xxxzzzz_0_xxyzzz_0, \
                             g_0_xxxzzzz_0_xxzzzz_0, \
                             g_0_xxxzzzz_0_xyyyyy_0, \
                             g_0_xxxzzzz_0_xyyyyz_0, \
                             g_0_xxxzzzz_0_xyyyzz_0, \
                             g_0_xxxzzzz_0_xyyzzz_0, \
                             g_0_xxxzzzz_0_xyzzzz_0, \
                             g_0_xxxzzzz_0_xzzzzz_0, \
                             g_0_xxxzzzz_0_yyyyyy_0, \
                             g_0_xxxzzzz_0_yyyyyz_0, \
                             g_0_xxxzzzz_0_yyyyzz_0, \
                             g_0_xxxzzzz_0_yyyzzz_0, \
                             g_0_xxxzzzz_0_yyzzzz_0, \
                             g_0_xxxzzzz_0_yzzzzz_0, \
                             g_0_xxxzzzz_0_zzzzzz_0, \
                             g_0_xxzzzz_0_xxxxxz_0,  \
                             g_0_xxzzzz_0_xxxxxz_1,  \
                             g_0_xxzzzz_0_xxxxyz_0,  \
                             g_0_xxzzzz_0_xxxxyz_1,  \
                             g_0_xxzzzz_0_xxxxz_1,   \
                             g_0_xxzzzz_0_xxxxzz_0,  \
                             g_0_xxzzzz_0_xxxxzz_1,  \
                             g_0_xxzzzz_0_xxxyyz_0,  \
                             g_0_xxzzzz_0_xxxyyz_1,  \
                             g_0_xxzzzz_0_xxxyz_1,   \
                             g_0_xxzzzz_0_xxxyzz_0,  \
                             g_0_xxzzzz_0_xxxyzz_1,  \
                             g_0_xxzzzz_0_xxxzz_1,   \
                             g_0_xxzzzz_0_xxxzzz_0,  \
                             g_0_xxzzzz_0_xxxzzz_1,  \
                             g_0_xxzzzz_0_xxyyyz_0,  \
                             g_0_xxzzzz_0_xxyyyz_1,  \
                             g_0_xxzzzz_0_xxyyz_1,   \
                             g_0_xxzzzz_0_xxyyzz_0,  \
                             g_0_xxzzzz_0_xxyyzz_1,  \
                             g_0_xxzzzz_0_xxyzz_1,   \
                             g_0_xxzzzz_0_xxyzzz_0,  \
                             g_0_xxzzzz_0_xxyzzz_1,  \
                             g_0_xxzzzz_0_xxzzz_1,   \
                             g_0_xxzzzz_0_xxzzzz_0,  \
                             g_0_xxzzzz_0_xxzzzz_1,  \
                             g_0_xxzzzz_0_xyyyyz_0,  \
                             g_0_xxzzzz_0_xyyyyz_1,  \
                             g_0_xxzzzz_0_xyyyz_1,   \
                             g_0_xxzzzz_0_xyyyzz_0,  \
                             g_0_xxzzzz_0_xyyyzz_1,  \
                             g_0_xxzzzz_0_xyyzz_1,   \
                             g_0_xxzzzz_0_xyyzzz_0,  \
                             g_0_xxzzzz_0_xyyzzz_1,  \
                             g_0_xxzzzz_0_xyzzz_1,   \
                             g_0_xxzzzz_0_xyzzzz_0,  \
                             g_0_xxzzzz_0_xyzzzz_1,  \
                             g_0_xxzzzz_0_xzzzz_1,   \
                             g_0_xxzzzz_0_xzzzzz_0,  \
                             g_0_xxzzzz_0_xzzzzz_1,  \
                             g_0_xxzzzz_0_yyyyyy_0,  \
                             g_0_xxzzzz_0_yyyyyy_1,  \
                             g_0_xxzzzz_0_yyyyyz_0,  \
                             g_0_xxzzzz_0_yyyyyz_1,  \
                             g_0_xxzzzz_0_yyyyz_1,   \
                             g_0_xxzzzz_0_yyyyzz_0,  \
                             g_0_xxzzzz_0_yyyyzz_1,  \
                             g_0_xxzzzz_0_yyyzz_1,   \
                             g_0_xxzzzz_0_yyyzzz_0,  \
                             g_0_xxzzzz_0_yyyzzz_1,  \
                             g_0_xxzzzz_0_yyzzz_1,   \
                             g_0_xxzzzz_0_yyzzzz_0,  \
                             g_0_xxzzzz_0_yyzzzz_1,  \
                             g_0_xxzzzz_0_yzzzz_1,   \
                             g_0_xxzzzz_0_yzzzzz_0,  \
                             g_0_xxzzzz_0_yzzzzz_1,  \
                             g_0_xxzzzz_0_zzzzz_1,   \
                             g_0_xxzzzz_0_zzzzzz_0,  \
                             g_0_xxzzzz_0_zzzzzz_1,  \
                             g_0_xzzzz_0_xxxxxz_0,   \
                             g_0_xzzzz_0_xxxxxz_1,   \
                             g_0_xzzzz_0_xxxxyz_0,   \
                             g_0_xzzzz_0_xxxxyz_1,   \
                             g_0_xzzzz_0_xxxxzz_0,   \
                             g_0_xzzzz_0_xxxxzz_1,   \
                             g_0_xzzzz_0_xxxyyz_0,   \
                             g_0_xzzzz_0_xxxyyz_1,   \
                             g_0_xzzzz_0_xxxyzz_0,   \
                             g_0_xzzzz_0_xxxyzz_1,   \
                             g_0_xzzzz_0_xxxzzz_0,   \
                             g_0_xzzzz_0_xxxzzz_1,   \
                             g_0_xzzzz_0_xxyyyz_0,   \
                             g_0_xzzzz_0_xxyyyz_1,   \
                             g_0_xzzzz_0_xxyyzz_0,   \
                             g_0_xzzzz_0_xxyyzz_1,   \
                             g_0_xzzzz_0_xxyzzz_0,   \
                             g_0_xzzzz_0_xxyzzz_1,   \
                             g_0_xzzzz_0_xxzzzz_0,   \
                             g_0_xzzzz_0_xxzzzz_1,   \
                             g_0_xzzzz_0_xyyyyz_0,   \
                             g_0_xzzzz_0_xyyyyz_1,   \
                             g_0_xzzzz_0_xyyyzz_0,   \
                             g_0_xzzzz_0_xyyyzz_1,   \
                             g_0_xzzzz_0_xyyzzz_0,   \
                             g_0_xzzzz_0_xyyzzz_1,   \
                             g_0_xzzzz_0_xyzzzz_0,   \
                             g_0_xzzzz_0_xyzzzz_1,   \
                             g_0_xzzzz_0_xzzzzz_0,   \
                             g_0_xzzzz_0_xzzzzz_1,   \
                             g_0_xzzzz_0_yyyyyy_0,   \
                             g_0_xzzzz_0_yyyyyy_1,   \
                             g_0_xzzzz_0_yyyyyz_0,   \
                             g_0_xzzzz_0_yyyyyz_1,   \
                             g_0_xzzzz_0_yyyyzz_0,   \
                             g_0_xzzzz_0_yyyyzz_1,   \
                             g_0_xzzzz_0_yyyzzz_0,   \
                             g_0_xzzzz_0_yyyzzz_1,   \
                             g_0_xzzzz_0_yyzzzz_0,   \
                             g_0_xzzzz_0_yyzzzz_1,   \
                             g_0_xzzzz_0_yzzzzz_0,   \
                             g_0_xzzzz_0_yzzzzz_1,   \
                             g_0_xzzzz_0_zzzzzz_0,   \
                             g_0_xzzzz_0_zzzzzz_1,   \
                             wp_x,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzzz_0_xxxxxx_0[i] = 3.0 * g_0_xxxzz_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxxxx_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_xxxxxx_0[i] * pb_z + g_0_xxxzzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxxxy_0[i] = 3.0 * g_0_xxxzz_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxxxy_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_xxxxxy_0[i] * pb_z + g_0_xxxzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxxxz_0[i] = 2.0 * g_0_xzzzz_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                    5.0 * g_0_xxzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxxz_0[i] * pb_x + g_0_xxzzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxxyy_0[i] = 3.0 * g_0_xxxzz_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxxyy_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_xxxxyy_0[i] * pb_z + g_0_xxxzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxxyz_0[i] = 2.0 * g_0_xzzzz_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxyz_0[i] * pb_x + g_0_xxzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxxzz_0[i] = 2.0 * g_0_xzzzz_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xxzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxzz_0[i] * pb_x + g_0_xxzzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxyyy_0[i] = 3.0 * g_0_xxxzz_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxxyyy_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_xxxyyy_0[i] * pb_z + g_0_xxxzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxxyyz_0[i] = 2.0 * g_0_xzzzz_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyyz_0[i] * pb_x + g_0_xxzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxyzz_0[i] = 2.0 * g_0_xzzzz_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyzz_0[i] * pb_x + g_0_xxzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxxzzz_0[i] = 2.0 * g_0_xzzzz_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xxzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxzzz_0[i] * pb_x + g_0_xxzzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxyyyy_0[i] = 3.0 * g_0_xxxzz_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xxyyyy_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_xxyyyy_0[i] * pb_z + g_0_xxxzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xxyyyz_0[i] = 2.0 * g_0_xzzzz_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyyz_0[i] * pb_x + g_0_xxzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxyyzz_0[i] = 2.0 * g_0_xzzzz_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyzz_0[i] * pb_x + g_0_xxzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxyzzz_0[i] = 2.0 * g_0_xzzzz_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyzzz_0[i] * pb_x + g_0_xxzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xxzzzz_0[i] = 2.0 * g_0_xzzzz_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xxzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxzzzz_0[i] * pb_x + g_0_xxzzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyyyyy_0[i] = 3.0 * g_0_xxxzz_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxzz_0_xyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxxzzz_0_xyyyyy_0[i] * pb_z + g_0_xxxzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxzzzz_0_xyyyyz_0[i] = 2.0 * g_0_xzzzz_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyyz_0[i] * pb_x + g_0_xxzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyyyzz_0[i] = 2.0 * g_0_xzzzz_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyzz_0[i] * pb_x + g_0_xxzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyyzzz_0[i] = 2.0 * g_0_xzzzz_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyzzz_0[i] * pb_x + g_0_xxzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xyzzzz_0[i] = 2.0 * g_0_xzzzz_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyzzzz_0[i] * pb_x + g_0_xxzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_xzzzzz_0[i] = 2.0 * g_0_xzzzz_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xzzzzz_0[i] * pb_x + g_0_xxzzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyyyy_0[i] = 2.0 * g_0_xzzzz_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_yyyyyy_0[i] * pb_x + g_0_xxzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyyyz_0[i] = 2.0 * g_0_xzzzz_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_yyyyyz_0[i] * pb_x + g_0_xxzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyyzz_0[i] = 2.0 * g_0_xzzzz_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_yyyyzz_0[i] * pb_x + g_0_xxzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyyzzz_0[i] = 2.0 * g_0_xzzzz_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_yyyzzz_0[i] * pb_x + g_0_xxzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yyzzzz_0[i] = 2.0 * g_0_xzzzz_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_yyzzzz_0[i] * pb_x + g_0_xxzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_yzzzzz_0[i] = 2.0 * g_0_xzzzz_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_yzzzzz_0[i] * pb_x + g_0_xxzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzz_0_zzzzzz_0[i] = 2.0 * g_0_xzzzz_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_zzzzzz_0[i] * pb_x + g_0_xxzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 420-448 components of targeted buffer : SKSI

    auto g_0_xxyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 420);

    auto g_0_xxyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 421);

    auto g_0_xxyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 422);

    auto g_0_xxyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 423);

    auto g_0_xxyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 424);

    auto g_0_xxyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 425);

    auto g_0_xxyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 426);

    auto g_0_xxyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 427);

    auto g_0_xxyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 428);

    auto g_0_xxyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 429);

    auto g_0_xxyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 430);

    auto g_0_xxyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 431);

    auto g_0_xxyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 432);

    auto g_0_xxyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 433);

    auto g_0_xxyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 434);

    auto g_0_xxyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 435);

    auto g_0_xxyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 436);

    auto g_0_xxyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 437);

    auto g_0_xxyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 438);

    auto g_0_xxyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 439);

    auto g_0_xxyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 440);

    auto g_0_xxyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 441);

    auto g_0_xxyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 442);

    auto g_0_xxyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 443);

    auto g_0_xxyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 444);

    auto g_0_xxyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 445);

    auto g_0_xxyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 446);

    auto g_0_xxyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 447);

#pragma omp simd aligned(g_0_xxyyy_0_xxxxxx_0,       \
                             g_0_xxyyy_0_xxxxxx_1,   \
                             g_0_xxyyy_0_xxxxxz_0,   \
                             g_0_xxyyy_0_xxxxxz_1,   \
                             g_0_xxyyy_0_xxxxzz_0,   \
                             g_0_xxyyy_0_xxxxzz_1,   \
                             g_0_xxyyy_0_xxxzzz_0,   \
                             g_0_xxyyy_0_xxxzzz_1,   \
                             g_0_xxyyy_0_xxzzzz_0,   \
                             g_0_xxyyy_0_xxzzzz_1,   \
                             g_0_xxyyy_0_xzzzzz_0,   \
                             g_0_xxyyy_0_xzzzzz_1,   \
                             g_0_xxyyyy_0_xxxxxx_0,  \
                             g_0_xxyyyy_0_xxxxxx_1,  \
                             g_0_xxyyyy_0_xxxxxz_0,  \
                             g_0_xxyyyy_0_xxxxxz_1,  \
                             g_0_xxyyyy_0_xxxxzz_0,  \
                             g_0_xxyyyy_0_xxxxzz_1,  \
                             g_0_xxyyyy_0_xxxzzz_0,  \
                             g_0_xxyyyy_0_xxxzzz_1,  \
                             g_0_xxyyyy_0_xxzzzz_0,  \
                             g_0_xxyyyy_0_xxzzzz_1,  \
                             g_0_xxyyyy_0_xzzzzz_0,  \
                             g_0_xxyyyy_0_xzzzzz_1,  \
                             g_0_xxyyyyy_0_xxxxxx_0, \
                             g_0_xxyyyyy_0_xxxxxy_0, \
                             g_0_xxyyyyy_0_xxxxxz_0, \
                             g_0_xxyyyyy_0_xxxxyy_0, \
                             g_0_xxyyyyy_0_xxxxyz_0, \
                             g_0_xxyyyyy_0_xxxxzz_0, \
                             g_0_xxyyyyy_0_xxxyyy_0, \
                             g_0_xxyyyyy_0_xxxyyz_0, \
                             g_0_xxyyyyy_0_xxxyzz_0, \
                             g_0_xxyyyyy_0_xxxzzz_0, \
                             g_0_xxyyyyy_0_xxyyyy_0, \
                             g_0_xxyyyyy_0_xxyyyz_0, \
                             g_0_xxyyyyy_0_xxyyzz_0, \
                             g_0_xxyyyyy_0_xxyzzz_0, \
                             g_0_xxyyyyy_0_xxzzzz_0, \
                             g_0_xxyyyyy_0_xyyyyy_0, \
                             g_0_xxyyyyy_0_xyyyyz_0, \
                             g_0_xxyyyyy_0_xyyyzz_0, \
                             g_0_xxyyyyy_0_xyyzzz_0, \
                             g_0_xxyyyyy_0_xyzzzz_0, \
                             g_0_xxyyyyy_0_xzzzzz_0, \
                             g_0_xxyyyyy_0_yyyyyy_0, \
                             g_0_xxyyyyy_0_yyyyyz_0, \
                             g_0_xxyyyyy_0_yyyyzz_0, \
                             g_0_xxyyyyy_0_yyyzzz_0, \
                             g_0_xxyyyyy_0_yyzzzz_0, \
                             g_0_xxyyyyy_0_yzzzzz_0, \
                             g_0_xxyyyyy_0_zzzzzz_0, \
                             g_0_xyyyyy_0_xxxxxy_0,  \
                             g_0_xyyyyy_0_xxxxxy_1,  \
                             g_0_xyyyyy_0_xxxxy_1,   \
                             g_0_xyyyyy_0_xxxxyy_0,  \
                             g_0_xyyyyy_0_xxxxyy_1,  \
                             g_0_xyyyyy_0_xxxxyz_0,  \
                             g_0_xyyyyy_0_xxxxyz_1,  \
                             g_0_xyyyyy_0_xxxyy_1,   \
                             g_0_xyyyyy_0_xxxyyy_0,  \
                             g_0_xyyyyy_0_xxxyyy_1,  \
                             g_0_xyyyyy_0_xxxyyz_0,  \
                             g_0_xyyyyy_0_xxxyyz_1,  \
                             g_0_xyyyyy_0_xxxyz_1,   \
                             g_0_xyyyyy_0_xxxyzz_0,  \
                             g_0_xyyyyy_0_xxxyzz_1,  \
                             g_0_xyyyyy_0_xxyyy_1,   \
                             g_0_xyyyyy_0_xxyyyy_0,  \
                             g_0_xyyyyy_0_xxyyyy_1,  \
                             g_0_xyyyyy_0_xxyyyz_0,  \
                             g_0_xyyyyy_0_xxyyyz_1,  \
                             g_0_xyyyyy_0_xxyyz_1,   \
                             g_0_xyyyyy_0_xxyyzz_0,  \
                             g_0_xyyyyy_0_xxyyzz_1,  \
                             g_0_xyyyyy_0_xxyzz_1,   \
                             g_0_xyyyyy_0_xxyzzz_0,  \
                             g_0_xyyyyy_0_xxyzzz_1,  \
                             g_0_xyyyyy_0_xyyyy_1,   \
                             g_0_xyyyyy_0_xyyyyy_0,  \
                             g_0_xyyyyy_0_xyyyyy_1,  \
                             g_0_xyyyyy_0_xyyyyz_0,  \
                             g_0_xyyyyy_0_xyyyyz_1,  \
                             g_0_xyyyyy_0_xyyyz_1,   \
                             g_0_xyyyyy_0_xyyyzz_0,  \
                             g_0_xyyyyy_0_xyyyzz_1,  \
                             g_0_xyyyyy_0_xyyzz_1,   \
                             g_0_xyyyyy_0_xyyzzz_0,  \
                             g_0_xyyyyy_0_xyyzzz_1,  \
                             g_0_xyyyyy_0_xyzzz_1,   \
                             g_0_xyyyyy_0_xyzzzz_0,  \
                             g_0_xyyyyy_0_xyzzzz_1,  \
                             g_0_xyyyyy_0_yyyyy_1,   \
                             g_0_xyyyyy_0_yyyyyy_0,  \
                             g_0_xyyyyy_0_yyyyyy_1,  \
                             g_0_xyyyyy_0_yyyyyz_0,  \
                             g_0_xyyyyy_0_yyyyyz_1,  \
                             g_0_xyyyyy_0_yyyyz_1,   \
                             g_0_xyyyyy_0_yyyyzz_0,  \
                             g_0_xyyyyy_0_yyyyzz_1,  \
                             g_0_xyyyyy_0_yyyzz_1,   \
                             g_0_xyyyyy_0_yyyzzz_0,  \
                             g_0_xyyyyy_0_yyyzzz_1,  \
                             g_0_xyyyyy_0_yyzzz_1,   \
                             g_0_xyyyyy_0_yyzzzz_0,  \
                             g_0_xyyyyy_0_yyzzzz_1,  \
                             g_0_xyyyyy_0_yzzzz_1,   \
                             g_0_xyyyyy_0_yzzzzz_0,  \
                             g_0_xyyyyy_0_yzzzzz_1,  \
                             g_0_xyyyyy_0_zzzzzz_0,  \
                             g_0_xyyyyy_0_zzzzzz_1,  \
                             g_0_yyyyy_0_xxxxxy_0,   \
                             g_0_yyyyy_0_xxxxxy_1,   \
                             g_0_yyyyy_0_xxxxyy_0,   \
                             g_0_yyyyy_0_xxxxyy_1,   \
                             g_0_yyyyy_0_xxxxyz_0,   \
                             g_0_yyyyy_0_xxxxyz_1,   \
                             g_0_yyyyy_0_xxxyyy_0,   \
                             g_0_yyyyy_0_xxxyyy_1,   \
                             g_0_yyyyy_0_xxxyyz_0,   \
                             g_0_yyyyy_0_xxxyyz_1,   \
                             g_0_yyyyy_0_xxxyzz_0,   \
                             g_0_yyyyy_0_xxxyzz_1,   \
                             g_0_yyyyy_0_xxyyyy_0,   \
                             g_0_yyyyy_0_xxyyyy_1,   \
                             g_0_yyyyy_0_xxyyyz_0,   \
                             g_0_yyyyy_0_xxyyyz_1,   \
                             g_0_yyyyy_0_xxyyzz_0,   \
                             g_0_yyyyy_0_xxyyzz_1,   \
                             g_0_yyyyy_0_xxyzzz_0,   \
                             g_0_yyyyy_0_xxyzzz_1,   \
                             g_0_yyyyy_0_xyyyyy_0,   \
                             g_0_yyyyy_0_xyyyyy_1,   \
                             g_0_yyyyy_0_xyyyyz_0,   \
                             g_0_yyyyy_0_xyyyyz_1,   \
                             g_0_yyyyy_0_xyyyzz_0,   \
                             g_0_yyyyy_0_xyyyzz_1,   \
                             g_0_yyyyy_0_xyyzzz_0,   \
                             g_0_yyyyy_0_xyyzzz_1,   \
                             g_0_yyyyy_0_xyzzzz_0,   \
                             g_0_yyyyy_0_xyzzzz_1,   \
                             g_0_yyyyy_0_yyyyyy_0,   \
                             g_0_yyyyy_0_yyyyyy_1,   \
                             g_0_yyyyy_0_yyyyyz_0,   \
                             g_0_yyyyy_0_yyyyyz_1,   \
                             g_0_yyyyy_0_yyyyzz_0,   \
                             g_0_yyyyy_0_yyyyzz_1,   \
                             g_0_yyyyy_0_yyyzzz_0,   \
                             g_0_yyyyy_0_yyyzzz_1,   \
                             g_0_yyyyy_0_yyzzzz_0,   \
                             g_0_yyyyy_0_yyzzzz_1,   \
                             g_0_yyyyy_0_yzzzzz_0,   \
                             g_0_yyyyy_0_yzzzzz_1,   \
                             g_0_yyyyy_0_zzzzzz_0,   \
                             g_0_yyyyy_0_zzzzzz_1,   \
                             wp_x,                   \
                             wp_y,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyy_0_xxxxxx_0[i] = 4.0 * g_0_xxyyy_0_xxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxxxx_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_xxxxxx_0[i] * pb_y + g_0_xxyyyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxxxxy_0[i] = g_0_yyyyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxxy_1[i] * fti_ab_0 +
                                    5.0 * g_0_xyyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxxxy_0[i] * pb_x + g_0_xyyyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxxxz_0[i] = 4.0 * g_0_xxyyy_0_xxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxxxz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_xxxxxz_0[i] * pb_y + g_0_xxyyyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxxxyy_0[i] = g_0_yyyyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxyy_1[i] * fti_ab_0 +
                                    4.0 * g_0_xyyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxxyy_0[i] * pb_x + g_0_xyyyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxxyz_0[i] = g_0_yyyyy_0_xxxxyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xyyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxxyz_0[i] * pb_x + g_0_xyyyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxxzz_0[i] = 4.0 * g_0_xxyyy_0_xxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxxzz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_xxxxzz_0[i] * pb_y + g_0_xxyyyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxxyyy_0[i] = g_0_yyyyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxyyy_1[i] * fti_ab_0 +
                                    3.0 * g_0_xyyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxyyy_0[i] * pb_x + g_0_xyyyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxyyz_0[i] = g_0_yyyyy_0_xxxyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xyyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxyyz_0[i] * pb_x + g_0_xyyyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxyzz_0[i] = g_0_yyyyy_0_xxxyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xyyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxxyzz_0[i] * pb_x + g_0_xyyyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxxzzz_0[i] = 4.0 * g_0_xxyyy_0_xxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxxzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_xxxzzz_0[i] * pb_y + g_0_xxyyyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xxyyyy_0[i] = g_0_yyyyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyyyy_1[i] * fti_ab_0 +
                                    2.0 * g_0_xyyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxyyyy_0[i] * pb_x + g_0_xyyyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxyyyz_0[i] = g_0_yyyyy_0_xxyyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xyyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxyyyz_0[i] * pb_x + g_0_xyyyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxyyzz_0[i] = g_0_yyyyy_0_xxyyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xyyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxyyzz_0[i] * pb_x + g_0_xyyyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxyzzz_0[i] = g_0_yyyyy_0_xxyzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xyyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xyyyyy_0_xxyzzz_0[i] * pb_x + g_0_xyyyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xxzzzz_0[i] = 4.0 * g_0_xxyyy_0_xxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xxzzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_xxzzzz_0[i] * pb_y + g_0_xxyyyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_xyyyyy_0[i] = g_0_yyyyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyy_1[i] * fi_abcd_0 +
                                    g_0_xyyyyy_0_xyyyyy_0[i] * pb_x + g_0_xyyyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyyyyz_0[i] = g_0_yyyyy_0_xyyyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyyyz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyz_1[i] * fi_abcd_0 +
                                    g_0_xyyyyy_0_xyyyyz_0[i] * pb_x + g_0_xyyyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyyyzz_0[i] = g_0_yyyyy_0_xyyyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyyzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyzz_1[i] * fi_abcd_0 +
                                    g_0_xyyyyy_0_xyyyzz_0[i] * pb_x + g_0_xyyyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyyzzz_0[i] = g_0_yyyyy_0_xyyzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyzzz_1[i] * fi_abcd_0 +
                                    g_0_xyyyyy_0_xyyzzz_0[i] * pb_x + g_0_xyyyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xyzzzz_0[i] = g_0_yyyyy_0_xyzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_xyzzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yzzzz_1[i] * fi_abcd_0 +
                                    g_0_xyyyyy_0_xyzzzz_0[i] * pb_x + g_0_xyyyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_xzzzzz_0[i] = 4.0 * g_0_xxyyy_0_xzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxyyy_0_xzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyyy_0_xzzzzz_0[i] * pb_y + g_0_xxyyyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyyyyy_0_yyyyyy_0[i] = g_0_yyyyy_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyyyy_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyyy_0[i] * pb_x +
                                    g_0_xyyyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyyyyz_0[i] = g_0_yyyyy_0_yyyyyz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyyyz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyyz_0[i] * pb_x +
                                    g_0_xyyyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyyyzz_0[i] = g_0_yyyyy_0_yyyyzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyyzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyyzz_0[i] * pb_x +
                                    g_0_xyyyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyyzzz_0[i] = g_0_yyyyy_0_yyyzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyyzzz_0[i] * pb_x +
                                    g_0_xyyyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yyzzzz_0[i] = g_0_yyyyy_0_yyzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yyzzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yyzzzz_0[i] * pb_x +
                                    g_0_xyyyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_yzzzzz_0[i] = g_0_yyyyy_0_yzzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_yzzzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_yzzzzz_0[i] * pb_x +
                                    g_0_xyyyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxyyyyy_0_zzzzzz_0[i] = g_0_yyyyy_0_zzzzzz_0[i] * fi_ab_0 - g_0_yyyyy_0_zzzzzz_1[i] * fti_ab_0 + g_0_xyyyyy_0_zzzzzz_0[i] * pb_x +
                                    g_0_xyyyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 448-476 components of targeted buffer : SKSI

    auto g_0_xxyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 448);

    auto g_0_xxyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 449);

    auto g_0_xxyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 450);

    auto g_0_xxyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 451);

    auto g_0_xxyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 452);

    auto g_0_xxyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 453);

    auto g_0_xxyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 454);

    auto g_0_xxyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 455);

    auto g_0_xxyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 456);

    auto g_0_xxyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 457);

    auto g_0_xxyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 458);

    auto g_0_xxyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 459);

    auto g_0_xxyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 460);

    auto g_0_xxyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 461);

    auto g_0_xxyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 462);

    auto g_0_xxyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 463);

    auto g_0_xxyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 464);

    auto g_0_xxyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 465);

    auto g_0_xxyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 466);

    auto g_0_xxyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 467);

    auto g_0_xxyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 468);

    auto g_0_xxyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 469);

    auto g_0_xxyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 470);

    auto g_0_xxyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 471);

    auto g_0_xxyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 472);

    auto g_0_xxyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 473);

    auto g_0_xxyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 474);

    auto g_0_xxyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 475);

#pragma omp simd aligned(g_0_xxyyyy_0_xxxxx_1,       \
                             g_0_xxyyyy_0_xxxxxx_0,  \
                             g_0_xxyyyy_0_xxxxxx_1,  \
                             g_0_xxyyyy_0_xxxxxy_0,  \
                             g_0_xxyyyy_0_xxxxxy_1,  \
                             g_0_xxyyyy_0_xxxxxz_0,  \
                             g_0_xxyyyy_0_xxxxxz_1,  \
                             g_0_xxyyyy_0_xxxxy_1,   \
                             g_0_xxyyyy_0_xxxxyy_0,  \
                             g_0_xxyyyy_0_xxxxyy_1,  \
                             g_0_xxyyyy_0_xxxxyz_0,  \
                             g_0_xxyyyy_0_xxxxyz_1,  \
                             g_0_xxyyyy_0_xxxxz_1,   \
                             g_0_xxyyyy_0_xxxxzz_0,  \
                             g_0_xxyyyy_0_xxxxzz_1,  \
                             g_0_xxyyyy_0_xxxyy_1,   \
                             g_0_xxyyyy_0_xxxyyy_0,  \
                             g_0_xxyyyy_0_xxxyyy_1,  \
                             g_0_xxyyyy_0_xxxyyz_0,  \
                             g_0_xxyyyy_0_xxxyyz_1,  \
                             g_0_xxyyyy_0_xxxyz_1,   \
                             g_0_xxyyyy_0_xxxyzz_0,  \
                             g_0_xxyyyy_0_xxxyzz_1,  \
                             g_0_xxyyyy_0_xxxzz_1,   \
                             g_0_xxyyyy_0_xxxzzz_0,  \
                             g_0_xxyyyy_0_xxxzzz_1,  \
                             g_0_xxyyyy_0_xxyyy_1,   \
                             g_0_xxyyyy_0_xxyyyy_0,  \
                             g_0_xxyyyy_0_xxyyyy_1,  \
                             g_0_xxyyyy_0_xxyyyz_0,  \
                             g_0_xxyyyy_0_xxyyyz_1,  \
                             g_0_xxyyyy_0_xxyyz_1,   \
                             g_0_xxyyyy_0_xxyyzz_0,  \
                             g_0_xxyyyy_0_xxyyzz_1,  \
                             g_0_xxyyyy_0_xxyzz_1,   \
                             g_0_xxyyyy_0_xxyzzz_0,  \
                             g_0_xxyyyy_0_xxyzzz_1,  \
                             g_0_xxyyyy_0_xxzzz_1,   \
                             g_0_xxyyyy_0_xxzzzz_0,  \
                             g_0_xxyyyy_0_xxzzzz_1,  \
                             g_0_xxyyyy_0_xyyyy_1,   \
                             g_0_xxyyyy_0_xyyyyy_0,  \
                             g_0_xxyyyy_0_xyyyyy_1,  \
                             g_0_xxyyyy_0_xyyyyz_0,  \
                             g_0_xxyyyy_0_xyyyyz_1,  \
                             g_0_xxyyyy_0_xyyyz_1,   \
                             g_0_xxyyyy_0_xyyyzz_0,  \
                             g_0_xxyyyy_0_xyyyzz_1,  \
                             g_0_xxyyyy_0_xyyzz_1,   \
                             g_0_xxyyyy_0_xyyzzz_0,  \
                             g_0_xxyyyy_0_xyyzzz_1,  \
                             g_0_xxyyyy_0_xyzzz_1,   \
                             g_0_xxyyyy_0_xyzzzz_0,  \
                             g_0_xxyyyy_0_xyzzzz_1,  \
                             g_0_xxyyyy_0_xzzzz_1,   \
                             g_0_xxyyyy_0_xzzzzz_0,  \
                             g_0_xxyyyy_0_xzzzzz_1,  \
                             g_0_xxyyyy_0_yyyyy_1,   \
                             g_0_xxyyyy_0_yyyyyy_0,  \
                             g_0_xxyyyy_0_yyyyyy_1,  \
                             g_0_xxyyyy_0_yyyyyz_0,  \
                             g_0_xxyyyy_0_yyyyyz_1,  \
                             g_0_xxyyyy_0_yyyyz_1,   \
                             g_0_xxyyyy_0_yyyyzz_0,  \
                             g_0_xxyyyy_0_yyyyzz_1,  \
                             g_0_xxyyyy_0_yyyzz_1,   \
                             g_0_xxyyyy_0_yyyzzz_0,  \
                             g_0_xxyyyy_0_yyyzzz_1,  \
                             g_0_xxyyyy_0_yyzzz_1,   \
                             g_0_xxyyyy_0_yyzzzz_0,  \
                             g_0_xxyyyy_0_yyzzzz_1,  \
                             g_0_xxyyyy_0_yzzzz_1,   \
                             g_0_xxyyyy_0_yzzzzz_0,  \
                             g_0_xxyyyy_0_yzzzzz_1,  \
                             g_0_xxyyyy_0_zzzzz_1,   \
                             g_0_xxyyyy_0_zzzzzz_0,  \
                             g_0_xxyyyy_0_zzzzzz_1,  \
                             g_0_xxyyyyz_0_xxxxxx_0, \
                             g_0_xxyyyyz_0_xxxxxy_0, \
                             g_0_xxyyyyz_0_xxxxxz_0, \
                             g_0_xxyyyyz_0_xxxxyy_0, \
                             g_0_xxyyyyz_0_xxxxyz_0, \
                             g_0_xxyyyyz_0_xxxxzz_0, \
                             g_0_xxyyyyz_0_xxxyyy_0, \
                             g_0_xxyyyyz_0_xxxyyz_0, \
                             g_0_xxyyyyz_0_xxxyzz_0, \
                             g_0_xxyyyyz_0_xxxzzz_0, \
                             g_0_xxyyyyz_0_xxyyyy_0, \
                             g_0_xxyyyyz_0_xxyyyz_0, \
                             g_0_xxyyyyz_0_xxyyzz_0, \
                             g_0_xxyyyyz_0_xxyzzz_0, \
                             g_0_xxyyyyz_0_xxzzzz_0, \
                             g_0_xxyyyyz_0_xyyyyy_0, \
                             g_0_xxyyyyz_0_xyyyyz_0, \
                             g_0_xxyyyyz_0_xyyyzz_0, \
                             g_0_xxyyyyz_0_xyyzzz_0, \
                             g_0_xxyyyyz_0_xyzzzz_0, \
                             g_0_xxyyyyz_0_xzzzzz_0, \
                             g_0_xxyyyyz_0_yyyyyy_0, \
                             g_0_xxyyyyz_0_yyyyyz_0, \
                             g_0_xxyyyyz_0_yyyyzz_0, \
                             g_0_xxyyyyz_0_yyyzzz_0, \
                             g_0_xxyyyyz_0_yyzzzz_0, \
                             g_0_xxyyyyz_0_yzzzzz_0, \
                             g_0_xxyyyyz_0_zzzzzz_0, \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyz_0_xxxxxx_0[i] = g_0_xxyyyy_0_xxxxxx_0[i] * pb_z + g_0_xxyyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxxy_0[i] = g_0_xxyyyy_0_xxxxxy_0[i] * pb_z + g_0_xxyyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxxz_0[i] = g_0_xxyyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxxz_0[i] * pb_z + g_0_xxyyyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxyy_0[i] = g_0_xxyyyy_0_xxxxyy_0[i] * pb_z + g_0_xxyyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxyz_0[i] = g_0_xxyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxyz_0[i] * pb_z + g_0_xxyyyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxxzz_0[i] = 2.0 * g_0_xxyyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxxzz_0[i] * pb_z + g_0_xxyyyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxyyy_0[i] = g_0_xxyyyy_0_xxxyyy_0[i] * pb_z + g_0_xxyyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxyyz_0[i] = g_0_xxyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyyz_0[i] * pb_z + g_0_xxyyyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxyzz_0[i] = 2.0 * g_0_xxyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxyzz_0[i] * pb_z + g_0_xxyyyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxxzzz_0[i] = 3.0 * g_0_xxyyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxxzzz_0[i] * pb_z + g_0_xxyyyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyyyy_0[i] = g_0_xxyyyy_0_xxyyyy_0[i] * pb_z + g_0_xxyyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyyyz_0[i] = g_0_xxyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyyz_0[i] * pb_z + g_0_xxyyyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyyzz_0[i] = 2.0 * g_0_xxyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyyzz_0[i] * pb_z + g_0_xxyyyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxyzzz_0[i] = 3.0 * g_0_xxyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxyzzz_0[i] * pb_z + g_0_xxyyyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xxzzzz_0[i] = 4.0 * g_0_xxyyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xxzzzz_0[i] * pb_z + g_0_xxyyyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyyyy_0[i] = g_0_xxyyyy_0_xyyyyy_0[i] * pb_z + g_0_xxyyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyyyz_0[i] = g_0_xxyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyyz_0[i] * pb_z + g_0_xxyyyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyyzz_0[i] = 2.0 * g_0_xxyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyyzz_0[i] * pb_z + g_0_xxyyyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyyzzz_0[i] = 3.0 * g_0_xxyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyyzzz_0[i] * pb_z + g_0_xxyyyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xyzzzz_0[i] = 4.0 * g_0_xxyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xyzzzz_0[i] * pb_z + g_0_xxyyyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_xzzzzz_0[i] = 5.0 * g_0_xxyyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_xzzzzz_0[i] * pb_z + g_0_xxyyyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyyyy_0[i] = g_0_xxyyyy_0_yyyyyy_0[i] * pb_z + g_0_xxyyyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyyyz_0[i] = g_0_xxyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyyyyz_0[i] * pb_z + g_0_xxyyyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyyzz_0[i] = 2.0 * g_0_xxyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyyyzz_0[i] * pb_z + g_0_xxyyyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyyzzz_0[i] = 3.0 * g_0_xxyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyyzzz_0[i] * pb_z + g_0_xxyyyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yyzzzz_0[i] = 4.0 * g_0_xxyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yyzzzz_0[i] * pb_z + g_0_xxyyyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_yzzzzz_0[i] = 5.0 * g_0_xxyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_yzzzzz_0[i] * pb_z + g_0_xxyyyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyz_0_zzzzzz_0[i] = 6.0 * g_0_xxyyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxyyyy_0_zzzzzz_0[i] * pb_z + g_0_xxyyyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 476-504 components of targeted buffer : SKSI

    auto g_0_xxyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 476);

    auto g_0_xxyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 477);

    auto g_0_xxyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 478);

    auto g_0_xxyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 479);

    auto g_0_xxyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 480);

    auto g_0_xxyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 481);

    auto g_0_xxyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 482);

    auto g_0_xxyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 483);

    auto g_0_xxyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 484);

    auto g_0_xxyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 485);

    auto g_0_xxyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 486);

    auto g_0_xxyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 487);

    auto g_0_xxyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 488);

    auto g_0_xxyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 489);

    auto g_0_xxyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 490);

    auto g_0_xxyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 491);

    auto g_0_xxyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 492);

    auto g_0_xxyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 493);

    auto g_0_xxyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 494);

    auto g_0_xxyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 495);

    auto g_0_xxyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 496);

    auto g_0_xxyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 497);

    auto g_0_xxyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 498);

    auto g_0_xxyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 499);

    auto g_0_xxyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 500);

    auto g_0_xxyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 501);

    auto g_0_xxyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 502);

    auto g_0_xxyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 503);

#pragma omp simd aligned(g_0_xxyyy_0_xxxxxy_0,       \
                             g_0_xxyyy_0_xxxxxy_1,   \
                             g_0_xxyyy_0_xxxxyy_0,   \
                             g_0_xxyyy_0_xxxxyy_1,   \
                             g_0_xxyyy_0_xxxyyy_0,   \
                             g_0_xxyyy_0_xxxyyy_1,   \
                             g_0_xxyyy_0_xxyyyy_0,   \
                             g_0_xxyyy_0_xxyyyy_1,   \
                             g_0_xxyyy_0_xyyyyy_0,   \
                             g_0_xxyyy_0_xyyyyy_1,   \
                             g_0_xxyyyz_0_xxxxxy_0,  \
                             g_0_xxyyyz_0_xxxxxy_1,  \
                             g_0_xxyyyz_0_xxxxyy_0,  \
                             g_0_xxyyyz_0_xxxxyy_1,  \
                             g_0_xxyyyz_0_xxxyyy_0,  \
                             g_0_xxyyyz_0_xxxyyy_1,  \
                             g_0_xxyyyz_0_xxyyyy_0,  \
                             g_0_xxyyyz_0_xxyyyy_1,  \
                             g_0_xxyyyz_0_xyyyyy_0,  \
                             g_0_xxyyyz_0_xyyyyy_1,  \
                             g_0_xxyyyzz_0_xxxxxx_0, \
                             g_0_xxyyyzz_0_xxxxxy_0, \
                             g_0_xxyyyzz_0_xxxxxz_0, \
                             g_0_xxyyyzz_0_xxxxyy_0, \
                             g_0_xxyyyzz_0_xxxxyz_0, \
                             g_0_xxyyyzz_0_xxxxzz_0, \
                             g_0_xxyyyzz_0_xxxyyy_0, \
                             g_0_xxyyyzz_0_xxxyyz_0, \
                             g_0_xxyyyzz_0_xxxyzz_0, \
                             g_0_xxyyyzz_0_xxxzzz_0, \
                             g_0_xxyyyzz_0_xxyyyy_0, \
                             g_0_xxyyyzz_0_xxyyyz_0, \
                             g_0_xxyyyzz_0_xxyyzz_0, \
                             g_0_xxyyyzz_0_xxyzzz_0, \
                             g_0_xxyyyzz_0_xxzzzz_0, \
                             g_0_xxyyyzz_0_xyyyyy_0, \
                             g_0_xxyyyzz_0_xyyyyz_0, \
                             g_0_xxyyyzz_0_xyyyzz_0, \
                             g_0_xxyyyzz_0_xyyzzz_0, \
                             g_0_xxyyyzz_0_xyzzzz_0, \
                             g_0_xxyyyzz_0_xzzzzz_0, \
                             g_0_xxyyyzz_0_yyyyyy_0, \
                             g_0_xxyyyzz_0_yyyyyz_0, \
                             g_0_xxyyyzz_0_yyyyzz_0, \
                             g_0_xxyyyzz_0_yyyzzz_0, \
                             g_0_xxyyyzz_0_yyzzzz_0, \
                             g_0_xxyyyzz_0_yzzzzz_0, \
                             g_0_xxyyyzz_0_zzzzzz_0, \
                             g_0_xxyyzz_0_xxxxxx_0,  \
                             g_0_xxyyzz_0_xxxxxx_1,  \
                             g_0_xxyyzz_0_xxxxxz_0,  \
                             g_0_xxyyzz_0_xxxxxz_1,  \
                             g_0_xxyyzz_0_xxxxzz_0,  \
                             g_0_xxyyzz_0_xxxxzz_1,  \
                             g_0_xxyyzz_0_xxxzzz_0,  \
                             g_0_xxyyzz_0_xxxzzz_1,  \
                             g_0_xxyyzz_0_xxzzzz_0,  \
                             g_0_xxyyzz_0_xxzzzz_1,  \
                             g_0_xxyyzz_0_xzzzzz_0,  \
                             g_0_xxyyzz_0_xzzzzz_1,  \
                             g_0_xxyzz_0_xxxxxx_0,   \
                             g_0_xxyzz_0_xxxxxx_1,   \
                             g_0_xxyzz_0_xxxxxz_0,   \
                             g_0_xxyzz_0_xxxxxz_1,   \
                             g_0_xxyzz_0_xxxxzz_0,   \
                             g_0_xxyzz_0_xxxxzz_1,   \
                             g_0_xxyzz_0_xxxzzz_0,   \
                             g_0_xxyzz_0_xxxzzz_1,   \
                             g_0_xxyzz_0_xxzzzz_0,   \
                             g_0_xxyzz_0_xxzzzz_1,   \
                             g_0_xxyzz_0_xzzzzz_0,   \
                             g_0_xxyzz_0_xzzzzz_1,   \
                             g_0_xyyyzz_0_xxxxyz_0,  \
                             g_0_xyyyzz_0_xxxxyz_1,  \
                             g_0_xyyyzz_0_xxxyyz_0,  \
                             g_0_xyyyzz_0_xxxyyz_1,  \
                             g_0_xyyyzz_0_xxxyz_1,   \
                             g_0_xyyyzz_0_xxxyzz_0,  \
                             g_0_xyyyzz_0_xxxyzz_1,  \
                             g_0_xyyyzz_0_xxyyyz_0,  \
                             g_0_xyyyzz_0_xxyyyz_1,  \
                             g_0_xyyyzz_0_xxyyz_1,   \
                             g_0_xyyyzz_0_xxyyzz_0,  \
                             g_0_xyyyzz_0_xxyyzz_1,  \
                             g_0_xyyyzz_0_xxyzz_1,   \
                             g_0_xyyyzz_0_xxyzzz_0,  \
                             g_0_xyyyzz_0_xxyzzz_1,  \
                             g_0_xyyyzz_0_xyyyyz_0,  \
                             g_0_xyyyzz_0_xyyyyz_1,  \
                             g_0_xyyyzz_0_xyyyz_1,   \
                             g_0_xyyyzz_0_xyyyzz_0,  \
                             g_0_xyyyzz_0_xyyyzz_1,  \
                             g_0_xyyyzz_0_xyyzz_1,   \
                             g_0_xyyyzz_0_xyyzzz_0,  \
                             g_0_xyyyzz_0_xyyzzz_1,  \
                             g_0_xyyyzz_0_xyzzz_1,   \
                             g_0_xyyyzz_0_xyzzzz_0,  \
                             g_0_xyyyzz_0_xyzzzz_1,  \
                             g_0_xyyyzz_0_yyyyyy_0,  \
                             g_0_xyyyzz_0_yyyyyy_1,  \
                             g_0_xyyyzz_0_yyyyyz_0,  \
                             g_0_xyyyzz_0_yyyyyz_1,  \
                             g_0_xyyyzz_0_yyyyz_1,   \
                             g_0_xyyyzz_0_yyyyzz_0,  \
                             g_0_xyyyzz_0_yyyyzz_1,  \
                             g_0_xyyyzz_0_yyyzz_1,   \
                             g_0_xyyyzz_0_yyyzzz_0,  \
                             g_0_xyyyzz_0_yyyzzz_1,  \
                             g_0_xyyyzz_0_yyzzz_1,   \
                             g_0_xyyyzz_0_yyzzzz_0,  \
                             g_0_xyyyzz_0_yyzzzz_1,  \
                             g_0_xyyyzz_0_yzzzz_1,   \
                             g_0_xyyyzz_0_yzzzzz_0,  \
                             g_0_xyyyzz_0_yzzzzz_1,  \
                             g_0_xyyyzz_0_zzzzzz_0,  \
                             g_0_xyyyzz_0_zzzzzz_1,  \
                             g_0_yyyzz_0_xxxxyz_0,   \
                             g_0_yyyzz_0_xxxxyz_1,   \
                             g_0_yyyzz_0_xxxyyz_0,   \
                             g_0_yyyzz_0_xxxyyz_1,   \
                             g_0_yyyzz_0_xxxyzz_0,   \
                             g_0_yyyzz_0_xxxyzz_1,   \
                             g_0_yyyzz_0_xxyyyz_0,   \
                             g_0_yyyzz_0_xxyyyz_1,   \
                             g_0_yyyzz_0_xxyyzz_0,   \
                             g_0_yyyzz_0_xxyyzz_1,   \
                             g_0_yyyzz_0_xxyzzz_0,   \
                             g_0_yyyzz_0_xxyzzz_1,   \
                             g_0_yyyzz_0_xyyyyz_0,   \
                             g_0_yyyzz_0_xyyyyz_1,   \
                             g_0_yyyzz_0_xyyyzz_0,   \
                             g_0_yyyzz_0_xyyyzz_1,   \
                             g_0_yyyzz_0_xyyzzz_0,   \
                             g_0_yyyzz_0_xyyzzz_1,   \
                             g_0_yyyzz_0_xyzzzz_0,   \
                             g_0_yyyzz_0_xyzzzz_1,   \
                             g_0_yyyzz_0_yyyyyy_0,   \
                             g_0_yyyzz_0_yyyyyy_1,   \
                             g_0_yyyzz_0_yyyyyz_0,   \
                             g_0_yyyzz_0_yyyyyz_1,   \
                             g_0_yyyzz_0_yyyyzz_0,   \
                             g_0_yyyzz_0_yyyyzz_1,   \
                             g_0_yyyzz_0_yyyzzz_0,   \
                             g_0_yyyzz_0_yyyzzz_1,   \
                             g_0_yyyzz_0_yyzzzz_0,   \
                             g_0_yyyzz_0_yyzzzz_1,   \
                             g_0_yyyzz_0_yzzzzz_0,   \
                             g_0_yyyzz_0_yzzzzz_1,   \
                             g_0_yyyzz_0_zzzzzz_0,   \
                             g_0_yyyzz_0_zzzzzz_1,   \
                             wp_x,                   \
                             wp_y,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyzz_0_xxxxxx_0[i] = 2.0 * g_0_xxyzz_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxxxx_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_xxxxxx_0[i] * pb_y + g_0_xxyyzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxxxxy_0[i] = g_0_xxyyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxxxxy_0[i] * pb_z +
                                    g_0_xxyyyz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxxxxz_0[i] = 2.0 * g_0_xxyzz_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxxxz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_xxxxxz_0[i] * pb_y + g_0_xxyyzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxxxyy_0[i] = g_0_xxyyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxxxyy_0[i] * pb_z +
                                    g_0_xxyyyz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxxxyz_0[i] = g_0_yyyzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxxxyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xyyyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxxxyz_0[i] * pb_x + g_0_xyyyzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxxxzz_0[i] = 2.0 * g_0_xxyzz_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxxzz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_xxxxzz_0[i] * pb_y + g_0_xxyyzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxxyyy_0[i] = g_0_xxyyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxxyyy_0[i] * pb_z +
                                    g_0_xxyyyz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxxyyz_0[i] = g_0_yyyzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxxyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xyyyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxxyyz_0[i] * pb_x + g_0_xyyyzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxxyzz_0[i] = g_0_yyyzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxxyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xyyyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxxyzz_0[i] * pb_x + g_0_xyyyzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxxzzz_0[i] = 2.0 * g_0_xxyzz_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxxzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_xxxzzz_0[i] * pb_y + g_0_xxyyzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xxyyyy_0[i] = g_0_xxyyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xxyyyy_0[i] * pb_z +
                                    g_0_xxyyyz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xxyyyz_0[i] = g_0_yyyzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxyyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xyyyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxyyyz_0[i] * pb_x + g_0_xyyyzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxyyzz_0[i] = g_0_yyyzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xyyyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxyyzz_0[i] * pb_x + g_0_xyyyzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxyzzz_0[i] = g_0_yyyzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xxyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xyyyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xyyyzz_0_xxyzzz_0[i] * pb_x + g_0_xyyyzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xxzzzz_0[i] = 2.0 * g_0_xxyzz_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xxzzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_xxzzzz_0[i] * pb_y + g_0_xxyyzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_xyyyyy_0[i] = g_0_xxyyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_xxyyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxyyyz_0_xyyyyy_0[i] * pb_z +
                                    g_0_xxyyyz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxyyyzz_0_xyyyyz_0[i] = g_0_yyyzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyyyyz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyyz_1[i] * fi_abcd_0 +
                                    g_0_xyyyzz_0_xyyyyz_0[i] * pb_x + g_0_xyyyzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xyyyzz_0[i] = g_0_yyyzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyyyzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyzz_1[i] * fi_abcd_0 +
                                    g_0_xyyyzz_0_xyyyzz_0[i] * pb_x + g_0_xyyyzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xyyzzz_0[i] = g_0_yyyzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyyzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyzzz_1[i] * fi_abcd_0 +
                                    g_0_xyyyzz_0_xyyzzz_0[i] * pb_x + g_0_xyyyzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xyzzzz_0[i] = g_0_yyyzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yzzzz_1[i] * fi_abcd_0 +
                                    g_0_xyyyzz_0_xyzzzz_0[i] * pb_x + g_0_xyyyzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_xzzzzz_0[i] = 2.0 * g_0_xxyzz_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzz_0_xzzzzz_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_xzzzzz_0[i] * pb_y + g_0_xxyyzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyyyzz_0_yyyyyy_0[i] = g_0_yyyzz_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyyyy_0[i] * pb_x +
                                    g_0_xyyyzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyyyyz_0[i] = g_0_yyyzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyyyz_0[i] * pb_x +
                                    g_0_xyyyzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyyyzz_0[i] = g_0_yyyzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyyzz_0[i] * pb_x +
                                    g_0_xyyyzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyyzzz_0[i] = g_0_yyyzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyyzzz_0[i] * pb_x +
                                    g_0_xyyyzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yyzzzz_0[i] = g_0_yyyzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yyzzzz_0[i] * pb_x +
                                    g_0_xyyyzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_yzzzzz_0[i] = g_0_yyyzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_yzzzzz_0[i] * pb_x +
                                    g_0_xyyyzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxyyyzz_0_zzzzzz_0[i] = g_0_yyyzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_yyyzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xyyyzz_0_zzzzzz_0[i] * pb_x +
                                    g_0_xyyyzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 504-532 components of targeted buffer : SKSI

    auto g_0_xxyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 504);

    auto g_0_xxyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 505);

    auto g_0_xxyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 506);

    auto g_0_xxyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 507);

    auto g_0_xxyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 508);

    auto g_0_xxyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 509);

    auto g_0_xxyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 510);

    auto g_0_xxyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 511);

    auto g_0_xxyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 512);

    auto g_0_xxyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 513);

    auto g_0_xxyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 514);

    auto g_0_xxyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 515);

    auto g_0_xxyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 516);

    auto g_0_xxyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 517);

    auto g_0_xxyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 518);

    auto g_0_xxyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 519);

    auto g_0_xxyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 520);

    auto g_0_xxyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 521);

    auto g_0_xxyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 522);

    auto g_0_xxyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 523);

    auto g_0_xxyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 524);

    auto g_0_xxyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 525);

    auto g_0_xxyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 526);

    auto g_0_xxyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 527);

    auto g_0_xxyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 528);

    auto g_0_xxyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 529);

    auto g_0_xxyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 530);

    auto g_0_xxyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 531);

#pragma omp simd aligned(g_0_xxyyz_0_xxxxxy_0,       \
                             g_0_xxyyz_0_xxxxxy_1,   \
                             g_0_xxyyz_0_xxxxyy_0,   \
                             g_0_xxyyz_0_xxxxyy_1,   \
                             g_0_xxyyz_0_xxxyyy_0,   \
                             g_0_xxyyz_0_xxxyyy_1,   \
                             g_0_xxyyz_0_xxyyyy_0,   \
                             g_0_xxyyz_0_xxyyyy_1,   \
                             g_0_xxyyz_0_xyyyyy_0,   \
                             g_0_xxyyz_0_xyyyyy_1,   \
                             g_0_xxyyzz_0_xxxxxy_0,  \
                             g_0_xxyyzz_0_xxxxxy_1,  \
                             g_0_xxyyzz_0_xxxxyy_0,  \
                             g_0_xxyyzz_0_xxxxyy_1,  \
                             g_0_xxyyzz_0_xxxyyy_0,  \
                             g_0_xxyyzz_0_xxxyyy_1,  \
                             g_0_xxyyzz_0_xxyyyy_0,  \
                             g_0_xxyyzz_0_xxyyyy_1,  \
                             g_0_xxyyzz_0_xyyyyy_0,  \
                             g_0_xxyyzz_0_xyyyyy_1,  \
                             g_0_xxyyzzz_0_xxxxxx_0, \
                             g_0_xxyyzzz_0_xxxxxy_0, \
                             g_0_xxyyzzz_0_xxxxxz_0, \
                             g_0_xxyyzzz_0_xxxxyy_0, \
                             g_0_xxyyzzz_0_xxxxyz_0, \
                             g_0_xxyyzzz_0_xxxxzz_0, \
                             g_0_xxyyzzz_0_xxxyyy_0, \
                             g_0_xxyyzzz_0_xxxyyz_0, \
                             g_0_xxyyzzz_0_xxxyzz_0, \
                             g_0_xxyyzzz_0_xxxzzz_0, \
                             g_0_xxyyzzz_0_xxyyyy_0, \
                             g_0_xxyyzzz_0_xxyyyz_0, \
                             g_0_xxyyzzz_0_xxyyzz_0, \
                             g_0_xxyyzzz_0_xxyzzz_0, \
                             g_0_xxyyzzz_0_xxzzzz_0, \
                             g_0_xxyyzzz_0_xyyyyy_0, \
                             g_0_xxyyzzz_0_xyyyyz_0, \
                             g_0_xxyyzzz_0_xyyyzz_0, \
                             g_0_xxyyzzz_0_xyyzzz_0, \
                             g_0_xxyyzzz_0_xyzzzz_0, \
                             g_0_xxyyzzz_0_xzzzzz_0, \
                             g_0_xxyyzzz_0_yyyyyy_0, \
                             g_0_xxyyzzz_0_yyyyyz_0, \
                             g_0_xxyyzzz_0_yyyyzz_0, \
                             g_0_xxyyzzz_0_yyyzzz_0, \
                             g_0_xxyyzzz_0_yyzzzz_0, \
                             g_0_xxyyzzz_0_yzzzzz_0, \
                             g_0_xxyyzzz_0_zzzzzz_0, \
                             g_0_xxyzzz_0_xxxxxx_0,  \
                             g_0_xxyzzz_0_xxxxxx_1,  \
                             g_0_xxyzzz_0_xxxxxz_0,  \
                             g_0_xxyzzz_0_xxxxxz_1,  \
                             g_0_xxyzzz_0_xxxxzz_0,  \
                             g_0_xxyzzz_0_xxxxzz_1,  \
                             g_0_xxyzzz_0_xxxzzz_0,  \
                             g_0_xxyzzz_0_xxxzzz_1,  \
                             g_0_xxyzzz_0_xxzzzz_0,  \
                             g_0_xxyzzz_0_xxzzzz_1,  \
                             g_0_xxyzzz_0_xzzzzz_0,  \
                             g_0_xxyzzz_0_xzzzzz_1,  \
                             g_0_xxzzz_0_xxxxxx_0,   \
                             g_0_xxzzz_0_xxxxxx_1,   \
                             g_0_xxzzz_0_xxxxxz_0,   \
                             g_0_xxzzz_0_xxxxxz_1,   \
                             g_0_xxzzz_0_xxxxzz_0,   \
                             g_0_xxzzz_0_xxxxzz_1,   \
                             g_0_xxzzz_0_xxxzzz_0,   \
                             g_0_xxzzz_0_xxxzzz_1,   \
                             g_0_xxzzz_0_xxzzzz_0,   \
                             g_0_xxzzz_0_xxzzzz_1,   \
                             g_0_xxzzz_0_xzzzzz_0,   \
                             g_0_xxzzz_0_xzzzzz_1,   \
                             g_0_xyyzzz_0_xxxxyz_0,  \
                             g_0_xyyzzz_0_xxxxyz_1,  \
                             g_0_xyyzzz_0_xxxyyz_0,  \
                             g_0_xyyzzz_0_xxxyyz_1,  \
                             g_0_xyyzzz_0_xxxyz_1,   \
                             g_0_xyyzzz_0_xxxyzz_0,  \
                             g_0_xyyzzz_0_xxxyzz_1,  \
                             g_0_xyyzzz_0_xxyyyz_0,  \
                             g_0_xyyzzz_0_xxyyyz_1,  \
                             g_0_xyyzzz_0_xxyyz_1,   \
                             g_0_xyyzzz_0_xxyyzz_0,  \
                             g_0_xyyzzz_0_xxyyzz_1,  \
                             g_0_xyyzzz_0_xxyzz_1,   \
                             g_0_xyyzzz_0_xxyzzz_0,  \
                             g_0_xyyzzz_0_xxyzzz_1,  \
                             g_0_xyyzzz_0_xyyyyz_0,  \
                             g_0_xyyzzz_0_xyyyyz_1,  \
                             g_0_xyyzzz_0_xyyyz_1,   \
                             g_0_xyyzzz_0_xyyyzz_0,  \
                             g_0_xyyzzz_0_xyyyzz_1,  \
                             g_0_xyyzzz_0_xyyzz_1,   \
                             g_0_xyyzzz_0_xyyzzz_0,  \
                             g_0_xyyzzz_0_xyyzzz_1,  \
                             g_0_xyyzzz_0_xyzzz_1,   \
                             g_0_xyyzzz_0_xyzzzz_0,  \
                             g_0_xyyzzz_0_xyzzzz_1,  \
                             g_0_xyyzzz_0_yyyyyy_0,  \
                             g_0_xyyzzz_0_yyyyyy_1,  \
                             g_0_xyyzzz_0_yyyyyz_0,  \
                             g_0_xyyzzz_0_yyyyyz_1,  \
                             g_0_xyyzzz_0_yyyyz_1,   \
                             g_0_xyyzzz_0_yyyyzz_0,  \
                             g_0_xyyzzz_0_yyyyzz_1,  \
                             g_0_xyyzzz_0_yyyzz_1,   \
                             g_0_xyyzzz_0_yyyzzz_0,  \
                             g_0_xyyzzz_0_yyyzzz_1,  \
                             g_0_xyyzzz_0_yyzzz_1,   \
                             g_0_xyyzzz_0_yyzzzz_0,  \
                             g_0_xyyzzz_0_yyzzzz_1,  \
                             g_0_xyyzzz_0_yzzzz_1,   \
                             g_0_xyyzzz_0_yzzzzz_0,  \
                             g_0_xyyzzz_0_yzzzzz_1,  \
                             g_0_xyyzzz_0_zzzzzz_0,  \
                             g_0_xyyzzz_0_zzzzzz_1,  \
                             g_0_yyzzz_0_xxxxyz_0,   \
                             g_0_yyzzz_0_xxxxyz_1,   \
                             g_0_yyzzz_0_xxxyyz_0,   \
                             g_0_yyzzz_0_xxxyyz_1,   \
                             g_0_yyzzz_0_xxxyzz_0,   \
                             g_0_yyzzz_0_xxxyzz_1,   \
                             g_0_yyzzz_0_xxyyyz_0,   \
                             g_0_yyzzz_0_xxyyyz_1,   \
                             g_0_yyzzz_0_xxyyzz_0,   \
                             g_0_yyzzz_0_xxyyzz_1,   \
                             g_0_yyzzz_0_xxyzzz_0,   \
                             g_0_yyzzz_0_xxyzzz_1,   \
                             g_0_yyzzz_0_xyyyyz_0,   \
                             g_0_yyzzz_0_xyyyyz_1,   \
                             g_0_yyzzz_0_xyyyzz_0,   \
                             g_0_yyzzz_0_xyyyzz_1,   \
                             g_0_yyzzz_0_xyyzzz_0,   \
                             g_0_yyzzz_0_xyyzzz_1,   \
                             g_0_yyzzz_0_xyzzzz_0,   \
                             g_0_yyzzz_0_xyzzzz_1,   \
                             g_0_yyzzz_0_yyyyyy_0,   \
                             g_0_yyzzz_0_yyyyyy_1,   \
                             g_0_yyzzz_0_yyyyyz_0,   \
                             g_0_yyzzz_0_yyyyyz_1,   \
                             g_0_yyzzz_0_yyyyzz_0,   \
                             g_0_yyzzz_0_yyyyzz_1,   \
                             g_0_yyzzz_0_yyyzzz_0,   \
                             g_0_yyzzz_0_yyyzzz_1,   \
                             g_0_yyzzz_0_yyzzzz_0,   \
                             g_0_yyzzz_0_yyzzzz_1,   \
                             g_0_yyzzz_0_yzzzzz_0,   \
                             g_0_yyzzz_0_yzzzzz_1,   \
                             g_0_yyzzz_0_zzzzzz_0,   \
                             g_0_yyzzz_0_zzzzzz_1,   \
                             wp_x,                   \
                             wp_y,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzzz_0_xxxxxx_0[i] = g_0_xxzzz_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxxxx_0[i] * pb_y +
                                    g_0_xxyzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxxxxy_0[i] = 2.0 * g_0_xxyyz_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxxxxy_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_xxxxxy_0[i] * pb_z + g_0_xxyyzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxxxxz_0[i] = g_0_xxzzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxxxz_0[i] * pb_y +
                                    g_0_xxyzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxxxyy_0[i] = 2.0 * g_0_xxyyz_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxxxyy_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_xxxxyy_0[i] * pb_z + g_0_xxyyzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxxxyz_0[i] = g_0_yyzzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xyyzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxxxyz_0[i] * pb_x + g_0_xyyzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxxxzz_0[i] = g_0_xxzzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxxzz_0[i] * pb_y +
                                    g_0_xxyzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxxyyy_0[i] = 2.0 * g_0_xxyyz_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxxyyy_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_xxxyyy_0[i] * pb_z + g_0_xxyyzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxxyyz_0[i] = g_0_yyzzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xyyzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxxyyz_0[i] * pb_x + g_0_xyyzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxxyzz_0[i] = g_0_yyzzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xyyzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxxyzz_0[i] * pb_x + g_0_xyyzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxxzzz_0[i] = g_0_xxzzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxxzzz_0[i] * pb_y +
                                    g_0_xxyzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xxyyyy_0[i] = 2.0 * g_0_xxyyz_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xxyyyy_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_xxyyyy_0[i] * pb_z + g_0_xxyyzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xxyyyz_0[i] = g_0_yyzzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xyyzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxyyyz_0[i] * pb_x + g_0_xyyzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxyyzz_0[i] = g_0_yyzzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xyyzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxyyzz_0[i] * pb_x + g_0_xyyzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxyzzz_0[i] = g_0_yyzzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xyyzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xyyzzz_0_xxyzzz_0[i] * pb_x + g_0_xyyzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xxzzzz_0[i] = g_0_xxzzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xxzzzz_0[i] * pb_y +
                                    g_0_xxyzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_xyyyyy_0[i] = 2.0 * g_0_xxyyz_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyz_0_xyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxyyzz_0_xyyyyy_0[i] * pb_z + g_0_xxyyzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxyyzzz_0_xyyyyz_0[i] = g_0_yyzzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyyyyz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyyz_1[i] * fi_abcd_0 +
                                    g_0_xyyzzz_0_xyyyyz_0[i] * pb_x + g_0_xyyzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xyyyzz_0[i] = g_0_yyzzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyyyzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyzz_1[i] * fi_abcd_0 +
                                    g_0_xyyzzz_0_xyyyzz_0[i] * pb_x + g_0_xyyzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xyyzzz_0[i] = g_0_yyzzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyyzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyzzz_1[i] * fi_abcd_0 +
                                    g_0_xyyzzz_0_xyyzzz_0[i] * pb_x + g_0_xyyzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xyzzzz_0[i] = g_0_yyzzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yzzzz_1[i] * fi_abcd_0 +
                                    g_0_xyyzzz_0_xyzzzz_0[i] * pb_x + g_0_xyyzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_xzzzzz_0[i] = g_0_xxzzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_xxzzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxyzzz_0_xzzzzz_0[i] * pb_y +
                                    g_0_xxyzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyyzzz_0_yyyyyy_0[i] = g_0_yyzzz_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyyyy_0[i] * pb_x +
                                    g_0_xyyzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyyyyz_0[i] = g_0_yyzzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyyyz_0[i] * pb_x +
                                    g_0_xyyzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyyyzz_0[i] = g_0_yyzzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyyzz_0[i] * pb_x +
                                    g_0_xyyzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyyzzz_0[i] = g_0_yyzzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyyzzz_0[i] * pb_x +
                                    g_0_xyyzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yyzzzz_0[i] = g_0_yyzzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yyzzzz_0[i] * pb_x +
                                    g_0_xyyzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_yzzzzz_0[i] = g_0_yyzzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_yzzzzz_0[i] * pb_x +
                                    g_0_xyyzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxyyzzz_0_zzzzzz_0[i] = g_0_yyzzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_yyzzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xyyzzz_0_zzzzzz_0[i] * pb_x +
                                    g_0_xyyzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 532-560 components of targeted buffer : SKSI

    auto g_0_xxyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 532);

    auto g_0_xxyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 533);

    auto g_0_xxyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 534);

    auto g_0_xxyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 535);

    auto g_0_xxyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 536);

    auto g_0_xxyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 537);

    auto g_0_xxyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 538);

    auto g_0_xxyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 539);

    auto g_0_xxyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 540);

    auto g_0_xxyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 541);

    auto g_0_xxyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 542);

    auto g_0_xxyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 543);

    auto g_0_xxyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 544);

    auto g_0_xxyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 545);

    auto g_0_xxyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 546);

    auto g_0_xxyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 547);

    auto g_0_xxyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 548);

    auto g_0_xxyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 549);

    auto g_0_xxyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 550);

    auto g_0_xxyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 551);

    auto g_0_xxyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 552);

    auto g_0_xxyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 553);

    auto g_0_xxyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 554);

    auto g_0_xxyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 555);

    auto g_0_xxyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 556);

    auto g_0_xxyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 557);

    auto g_0_xxyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 558);

    auto g_0_xxyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 559);

#pragma omp simd aligned(g_0_xxyzzzz_0_xxxxxx_0,     \
                             g_0_xxyzzzz_0_xxxxxy_0, \
                             g_0_xxyzzzz_0_xxxxxz_0, \
                             g_0_xxyzzzz_0_xxxxyy_0, \
                             g_0_xxyzzzz_0_xxxxyz_0, \
                             g_0_xxyzzzz_0_xxxxzz_0, \
                             g_0_xxyzzzz_0_xxxyyy_0, \
                             g_0_xxyzzzz_0_xxxyyz_0, \
                             g_0_xxyzzzz_0_xxxyzz_0, \
                             g_0_xxyzzzz_0_xxxzzz_0, \
                             g_0_xxyzzzz_0_xxyyyy_0, \
                             g_0_xxyzzzz_0_xxyyyz_0, \
                             g_0_xxyzzzz_0_xxyyzz_0, \
                             g_0_xxyzzzz_0_xxyzzz_0, \
                             g_0_xxyzzzz_0_xxzzzz_0, \
                             g_0_xxyzzzz_0_xyyyyy_0, \
                             g_0_xxyzzzz_0_xyyyyz_0, \
                             g_0_xxyzzzz_0_xyyyzz_0, \
                             g_0_xxyzzzz_0_xyyzzz_0, \
                             g_0_xxyzzzz_0_xyzzzz_0, \
                             g_0_xxyzzzz_0_xzzzzz_0, \
                             g_0_xxyzzzz_0_yyyyyy_0, \
                             g_0_xxyzzzz_0_yyyyyz_0, \
                             g_0_xxyzzzz_0_yyyyzz_0, \
                             g_0_xxyzzzz_0_yyyzzz_0, \
                             g_0_xxyzzzz_0_yyzzzz_0, \
                             g_0_xxyzzzz_0_yzzzzz_0, \
                             g_0_xxyzzzz_0_zzzzzz_0, \
                             g_0_xxzzzz_0_xxxxx_1,   \
                             g_0_xxzzzz_0_xxxxxx_0,  \
                             g_0_xxzzzz_0_xxxxxx_1,  \
                             g_0_xxzzzz_0_xxxxxy_0,  \
                             g_0_xxzzzz_0_xxxxxy_1,  \
                             g_0_xxzzzz_0_xxxxxz_0,  \
                             g_0_xxzzzz_0_xxxxxz_1,  \
                             g_0_xxzzzz_0_xxxxy_1,   \
                             g_0_xxzzzz_0_xxxxyy_0,  \
                             g_0_xxzzzz_0_xxxxyy_1,  \
                             g_0_xxzzzz_0_xxxxyz_0,  \
                             g_0_xxzzzz_0_xxxxyz_1,  \
                             g_0_xxzzzz_0_xxxxz_1,   \
                             g_0_xxzzzz_0_xxxxzz_0,  \
                             g_0_xxzzzz_0_xxxxzz_1,  \
                             g_0_xxzzzz_0_xxxyy_1,   \
                             g_0_xxzzzz_0_xxxyyy_0,  \
                             g_0_xxzzzz_0_xxxyyy_1,  \
                             g_0_xxzzzz_0_xxxyyz_0,  \
                             g_0_xxzzzz_0_xxxyyz_1,  \
                             g_0_xxzzzz_0_xxxyz_1,   \
                             g_0_xxzzzz_0_xxxyzz_0,  \
                             g_0_xxzzzz_0_xxxyzz_1,  \
                             g_0_xxzzzz_0_xxxzz_1,   \
                             g_0_xxzzzz_0_xxxzzz_0,  \
                             g_0_xxzzzz_0_xxxzzz_1,  \
                             g_0_xxzzzz_0_xxyyy_1,   \
                             g_0_xxzzzz_0_xxyyyy_0,  \
                             g_0_xxzzzz_0_xxyyyy_1,  \
                             g_0_xxzzzz_0_xxyyyz_0,  \
                             g_0_xxzzzz_0_xxyyyz_1,  \
                             g_0_xxzzzz_0_xxyyz_1,   \
                             g_0_xxzzzz_0_xxyyzz_0,  \
                             g_0_xxzzzz_0_xxyyzz_1,  \
                             g_0_xxzzzz_0_xxyzz_1,   \
                             g_0_xxzzzz_0_xxyzzz_0,  \
                             g_0_xxzzzz_0_xxyzzz_1,  \
                             g_0_xxzzzz_0_xxzzz_1,   \
                             g_0_xxzzzz_0_xxzzzz_0,  \
                             g_0_xxzzzz_0_xxzzzz_1,  \
                             g_0_xxzzzz_0_xyyyy_1,   \
                             g_0_xxzzzz_0_xyyyyy_0,  \
                             g_0_xxzzzz_0_xyyyyy_1,  \
                             g_0_xxzzzz_0_xyyyyz_0,  \
                             g_0_xxzzzz_0_xyyyyz_1,  \
                             g_0_xxzzzz_0_xyyyz_1,   \
                             g_0_xxzzzz_0_xyyyzz_0,  \
                             g_0_xxzzzz_0_xyyyzz_1,  \
                             g_0_xxzzzz_0_xyyzz_1,   \
                             g_0_xxzzzz_0_xyyzzz_0,  \
                             g_0_xxzzzz_0_xyyzzz_1,  \
                             g_0_xxzzzz_0_xyzzz_1,   \
                             g_0_xxzzzz_0_xyzzzz_0,  \
                             g_0_xxzzzz_0_xyzzzz_1,  \
                             g_0_xxzzzz_0_xzzzz_1,   \
                             g_0_xxzzzz_0_xzzzzz_0,  \
                             g_0_xxzzzz_0_xzzzzz_1,  \
                             g_0_xxzzzz_0_yyyyy_1,   \
                             g_0_xxzzzz_0_yyyyyy_0,  \
                             g_0_xxzzzz_0_yyyyyy_1,  \
                             g_0_xxzzzz_0_yyyyyz_0,  \
                             g_0_xxzzzz_0_yyyyyz_1,  \
                             g_0_xxzzzz_0_yyyyz_1,   \
                             g_0_xxzzzz_0_yyyyzz_0,  \
                             g_0_xxzzzz_0_yyyyzz_1,  \
                             g_0_xxzzzz_0_yyyzz_1,   \
                             g_0_xxzzzz_0_yyyzzz_0,  \
                             g_0_xxzzzz_0_yyyzzz_1,  \
                             g_0_xxzzzz_0_yyzzz_1,   \
                             g_0_xxzzzz_0_yyzzzz_0,  \
                             g_0_xxzzzz_0_yyzzzz_1,  \
                             g_0_xxzzzz_0_yzzzz_1,   \
                             g_0_xxzzzz_0_yzzzzz_0,  \
                             g_0_xxzzzz_0_yzzzzz_1,  \
                             g_0_xxzzzz_0_zzzzz_1,   \
                             g_0_xxzzzz_0_zzzzzz_0,  \
                             g_0_xxzzzz_0_zzzzzz_1,  \
                             wp_y,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzz_0_xxxxxx_0[i] = g_0_xxzzzz_0_xxxxxx_0[i] * pb_y + g_0_xxzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxxy_0[i] = g_0_xxzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxxy_0[i] * pb_y + g_0_xxzzzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxxz_0[i] = g_0_xxzzzz_0_xxxxxz_0[i] * pb_y + g_0_xxzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxyy_0[i] = 2.0 * g_0_xxzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxyy_0[i] * pb_y + g_0_xxzzzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxyz_0[i] = g_0_xxzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxxyz_0[i] * pb_y + g_0_xxzzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxxzz_0[i] = g_0_xxzzzz_0_xxxxzz_0[i] * pb_y + g_0_xxzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxyyy_0[i] = 3.0 * g_0_xxzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyyy_0[i] * pb_y + g_0_xxzzzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxyyz_0[i] = 2.0 * g_0_xxzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyyz_0[i] * pb_y + g_0_xxzzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxyzz_0[i] = g_0_xxzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxxyzz_0[i] * pb_y + g_0_xxzzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxxzzz_0[i] = g_0_xxzzzz_0_xxxzzz_0[i] * pb_y + g_0_xxzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyyyy_0[i] = 4.0 * g_0_xxzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyyy_0[i] * pb_y + g_0_xxzzzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyyyz_0[i] = 3.0 * g_0_xxzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyyz_0[i] * pb_y + g_0_xxzzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyyzz_0[i] = 2.0 * g_0_xxzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyyzz_0[i] * pb_y + g_0_xxzzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxyzzz_0[i] = g_0_xxzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xxyzzz_0[i] * pb_y + g_0_xxzzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xxzzzz_0[i] = g_0_xxzzzz_0_xxzzzz_0[i] * pb_y + g_0_xxzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyyyy_0[i] = 5.0 * g_0_xxzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyyy_0[i] * pb_y + g_0_xxzzzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyyyz_0[i] = 4.0 * g_0_xxzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyyz_0[i] * pb_y + g_0_xxzzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyyzz_0[i] = 3.0 * g_0_xxzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyyzz_0[i] * pb_y + g_0_xxzzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyyzzz_0[i] = 2.0 * g_0_xxzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyyzzz_0[i] * pb_y + g_0_xxzzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xyzzzz_0[i] = g_0_xxzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_xyzzzz_0[i] * pb_y + g_0_xxzzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_xzzzzz_0[i] = g_0_xxzzzz_0_xzzzzz_0[i] * pb_y + g_0_xxzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyyyy_0[i] = 6.0 * g_0_xxzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyyyy_0[i] * pb_y + g_0_xxzzzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyyyz_0[i] = 5.0 * g_0_xxzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyyyz_0[i] * pb_y + g_0_xxzzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyyzz_0[i] = 4.0 * g_0_xxzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyyzz_0[i] * pb_y + g_0_xxzzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyyzzz_0[i] = 3.0 * g_0_xxzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyyzzz_0[i] * pb_y + g_0_xxzzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yyzzzz_0[i] = 2.0 * g_0_xxzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yyzzzz_0[i] * pb_y + g_0_xxzzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_yzzzzz_0[i] = g_0_xxzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxzzzz_0_yzzzzz_0[i] * pb_y + g_0_xxzzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzz_0_zzzzzz_0[i] = g_0_xxzzzz_0_zzzzzz_0[i] * pb_y + g_0_xxzzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 560-588 components of targeted buffer : SKSI

    auto g_0_xxzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 560);

    auto g_0_xxzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 561);

    auto g_0_xxzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 562);

    auto g_0_xxzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 563);

    auto g_0_xxzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 564);

    auto g_0_xxzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 565);

    auto g_0_xxzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 566);

    auto g_0_xxzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 567);

    auto g_0_xxzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 568);

    auto g_0_xxzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 569);

    auto g_0_xxzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 570);

    auto g_0_xxzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 571);

    auto g_0_xxzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 572);

    auto g_0_xxzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 573);

    auto g_0_xxzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 574);

    auto g_0_xxzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 575);

    auto g_0_xxzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 576);

    auto g_0_xxzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 577);

    auto g_0_xxzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 578);

    auto g_0_xxzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 579);

    auto g_0_xxzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 580);

    auto g_0_xxzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 581);

    auto g_0_xxzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 582);

    auto g_0_xxzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 583);

    auto g_0_xxzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 584);

    auto g_0_xxzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 585);

    auto g_0_xxzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 586);

    auto g_0_xxzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 587);

#pragma omp simd aligned(g_0_xxzzz_0_xxxxxx_0,       \
                             g_0_xxzzz_0_xxxxxx_1,   \
                             g_0_xxzzz_0_xxxxxy_0,   \
                             g_0_xxzzz_0_xxxxxy_1,   \
                             g_0_xxzzz_0_xxxxyy_0,   \
                             g_0_xxzzz_0_xxxxyy_1,   \
                             g_0_xxzzz_0_xxxyyy_0,   \
                             g_0_xxzzz_0_xxxyyy_1,   \
                             g_0_xxzzz_0_xxyyyy_0,   \
                             g_0_xxzzz_0_xxyyyy_1,   \
                             g_0_xxzzz_0_xyyyyy_0,   \
                             g_0_xxzzz_0_xyyyyy_1,   \
                             g_0_xxzzzz_0_xxxxxx_0,  \
                             g_0_xxzzzz_0_xxxxxx_1,  \
                             g_0_xxzzzz_0_xxxxxy_0,  \
                             g_0_xxzzzz_0_xxxxxy_1,  \
                             g_0_xxzzzz_0_xxxxyy_0,  \
                             g_0_xxzzzz_0_xxxxyy_1,  \
                             g_0_xxzzzz_0_xxxyyy_0,  \
                             g_0_xxzzzz_0_xxxyyy_1,  \
                             g_0_xxzzzz_0_xxyyyy_0,  \
                             g_0_xxzzzz_0_xxyyyy_1,  \
                             g_0_xxzzzz_0_xyyyyy_0,  \
                             g_0_xxzzzz_0_xyyyyy_1,  \
                             g_0_xxzzzzz_0_xxxxxx_0, \
                             g_0_xxzzzzz_0_xxxxxy_0, \
                             g_0_xxzzzzz_0_xxxxxz_0, \
                             g_0_xxzzzzz_0_xxxxyy_0, \
                             g_0_xxzzzzz_0_xxxxyz_0, \
                             g_0_xxzzzzz_0_xxxxzz_0, \
                             g_0_xxzzzzz_0_xxxyyy_0, \
                             g_0_xxzzzzz_0_xxxyyz_0, \
                             g_0_xxzzzzz_0_xxxyzz_0, \
                             g_0_xxzzzzz_0_xxxzzz_0, \
                             g_0_xxzzzzz_0_xxyyyy_0, \
                             g_0_xxzzzzz_0_xxyyyz_0, \
                             g_0_xxzzzzz_0_xxyyzz_0, \
                             g_0_xxzzzzz_0_xxyzzz_0, \
                             g_0_xxzzzzz_0_xxzzzz_0, \
                             g_0_xxzzzzz_0_xyyyyy_0, \
                             g_0_xxzzzzz_0_xyyyyz_0, \
                             g_0_xxzzzzz_0_xyyyzz_0, \
                             g_0_xxzzzzz_0_xyyzzz_0, \
                             g_0_xxzzzzz_0_xyzzzz_0, \
                             g_0_xxzzzzz_0_xzzzzz_0, \
                             g_0_xxzzzzz_0_yyyyyy_0, \
                             g_0_xxzzzzz_0_yyyyyz_0, \
                             g_0_xxzzzzz_0_yyyyzz_0, \
                             g_0_xxzzzzz_0_yyyzzz_0, \
                             g_0_xxzzzzz_0_yyzzzz_0, \
                             g_0_xxzzzzz_0_yzzzzz_0, \
                             g_0_xxzzzzz_0_zzzzzz_0, \
                             g_0_xzzzzz_0_xxxxxz_0,  \
                             g_0_xzzzzz_0_xxxxxz_1,  \
                             g_0_xzzzzz_0_xxxxyz_0,  \
                             g_0_xzzzzz_0_xxxxyz_1,  \
                             g_0_xzzzzz_0_xxxxz_1,   \
                             g_0_xzzzzz_0_xxxxzz_0,  \
                             g_0_xzzzzz_0_xxxxzz_1,  \
                             g_0_xzzzzz_0_xxxyyz_0,  \
                             g_0_xzzzzz_0_xxxyyz_1,  \
                             g_0_xzzzzz_0_xxxyz_1,   \
                             g_0_xzzzzz_0_xxxyzz_0,  \
                             g_0_xzzzzz_0_xxxyzz_1,  \
                             g_0_xzzzzz_0_xxxzz_1,   \
                             g_0_xzzzzz_0_xxxzzz_0,  \
                             g_0_xzzzzz_0_xxxzzz_1,  \
                             g_0_xzzzzz_0_xxyyyz_0,  \
                             g_0_xzzzzz_0_xxyyyz_1,  \
                             g_0_xzzzzz_0_xxyyz_1,   \
                             g_0_xzzzzz_0_xxyyzz_0,  \
                             g_0_xzzzzz_0_xxyyzz_1,  \
                             g_0_xzzzzz_0_xxyzz_1,   \
                             g_0_xzzzzz_0_xxyzzz_0,  \
                             g_0_xzzzzz_0_xxyzzz_1,  \
                             g_0_xzzzzz_0_xxzzz_1,   \
                             g_0_xzzzzz_0_xxzzzz_0,  \
                             g_0_xzzzzz_0_xxzzzz_1,  \
                             g_0_xzzzzz_0_xyyyyz_0,  \
                             g_0_xzzzzz_0_xyyyyz_1,  \
                             g_0_xzzzzz_0_xyyyz_1,   \
                             g_0_xzzzzz_0_xyyyzz_0,  \
                             g_0_xzzzzz_0_xyyyzz_1,  \
                             g_0_xzzzzz_0_xyyzz_1,   \
                             g_0_xzzzzz_0_xyyzzz_0,  \
                             g_0_xzzzzz_0_xyyzzz_1,  \
                             g_0_xzzzzz_0_xyzzz_1,   \
                             g_0_xzzzzz_0_xyzzzz_0,  \
                             g_0_xzzzzz_0_xyzzzz_1,  \
                             g_0_xzzzzz_0_xzzzz_1,   \
                             g_0_xzzzzz_0_xzzzzz_0,  \
                             g_0_xzzzzz_0_xzzzzz_1,  \
                             g_0_xzzzzz_0_yyyyyy_0,  \
                             g_0_xzzzzz_0_yyyyyy_1,  \
                             g_0_xzzzzz_0_yyyyyz_0,  \
                             g_0_xzzzzz_0_yyyyyz_1,  \
                             g_0_xzzzzz_0_yyyyz_1,   \
                             g_0_xzzzzz_0_yyyyzz_0,  \
                             g_0_xzzzzz_0_yyyyzz_1,  \
                             g_0_xzzzzz_0_yyyzz_1,   \
                             g_0_xzzzzz_0_yyyzzz_0,  \
                             g_0_xzzzzz_0_yyyzzz_1,  \
                             g_0_xzzzzz_0_yyzzz_1,   \
                             g_0_xzzzzz_0_yyzzzz_0,  \
                             g_0_xzzzzz_0_yyzzzz_1,  \
                             g_0_xzzzzz_0_yzzzz_1,   \
                             g_0_xzzzzz_0_yzzzzz_0,  \
                             g_0_xzzzzz_0_yzzzzz_1,  \
                             g_0_xzzzzz_0_zzzzz_1,   \
                             g_0_xzzzzz_0_zzzzzz_0,  \
                             g_0_xzzzzz_0_zzzzzz_1,  \
                             g_0_zzzzz_0_xxxxxz_0,   \
                             g_0_zzzzz_0_xxxxxz_1,   \
                             g_0_zzzzz_0_xxxxyz_0,   \
                             g_0_zzzzz_0_xxxxyz_1,   \
                             g_0_zzzzz_0_xxxxzz_0,   \
                             g_0_zzzzz_0_xxxxzz_1,   \
                             g_0_zzzzz_0_xxxyyz_0,   \
                             g_0_zzzzz_0_xxxyyz_1,   \
                             g_0_zzzzz_0_xxxyzz_0,   \
                             g_0_zzzzz_0_xxxyzz_1,   \
                             g_0_zzzzz_0_xxxzzz_0,   \
                             g_0_zzzzz_0_xxxzzz_1,   \
                             g_0_zzzzz_0_xxyyyz_0,   \
                             g_0_zzzzz_0_xxyyyz_1,   \
                             g_0_zzzzz_0_xxyyzz_0,   \
                             g_0_zzzzz_0_xxyyzz_1,   \
                             g_0_zzzzz_0_xxyzzz_0,   \
                             g_0_zzzzz_0_xxyzzz_1,   \
                             g_0_zzzzz_0_xxzzzz_0,   \
                             g_0_zzzzz_0_xxzzzz_1,   \
                             g_0_zzzzz_0_xyyyyz_0,   \
                             g_0_zzzzz_0_xyyyyz_1,   \
                             g_0_zzzzz_0_xyyyzz_0,   \
                             g_0_zzzzz_0_xyyyzz_1,   \
                             g_0_zzzzz_0_xyyzzz_0,   \
                             g_0_zzzzz_0_xyyzzz_1,   \
                             g_0_zzzzz_0_xyzzzz_0,   \
                             g_0_zzzzz_0_xyzzzz_1,   \
                             g_0_zzzzz_0_xzzzzz_0,   \
                             g_0_zzzzz_0_xzzzzz_1,   \
                             g_0_zzzzz_0_yyyyyy_0,   \
                             g_0_zzzzz_0_yyyyyy_1,   \
                             g_0_zzzzz_0_yyyyyz_0,   \
                             g_0_zzzzz_0_yyyyyz_1,   \
                             g_0_zzzzz_0_yyyyzz_0,   \
                             g_0_zzzzz_0_yyyyzz_1,   \
                             g_0_zzzzz_0_yyyzzz_0,   \
                             g_0_zzzzz_0_yyyzzz_1,   \
                             g_0_zzzzz_0_yyzzzz_0,   \
                             g_0_zzzzz_0_yyzzzz_1,   \
                             g_0_zzzzz_0_yzzzzz_0,   \
                             g_0_zzzzz_0_yzzzzz_1,   \
                             g_0_zzzzz_0_zzzzzz_0,   \
                             g_0_zzzzz_0_zzzzzz_1,   \
                             wp_x,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzzz_0_xxxxxx_0[i] = 4.0 * g_0_xxzzz_0_xxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxxxx_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_xxxxxx_0[i] * pb_z + g_0_xxzzzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxxxy_0[i] = 4.0 * g_0_xxzzz_0_xxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxxxy_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_xxxxxy_0[i] * pb_z + g_0_xxzzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxxxz_0[i] = g_0_zzzzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                    5.0 * g_0_xzzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxxxz_0[i] * pb_x + g_0_xzzzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxxyy_0[i] = 4.0 * g_0_xxzzz_0_xxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxxyy_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_xxxxyy_0[i] * pb_z + g_0_xxzzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxxyz_0[i] = g_0_zzzzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxxyz_0[i] * pb_x + g_0_xzzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxxzz_0[i] = g_0_zzzzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_xzzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxxzz_0[i] * pb_x + g_0_xzzzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxyyy_0[i] = 4.0 * g_0_xxzzz_0_xxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxxyyy_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_xxxyyy_0[i] * pb_z + g_0_xxzzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxxyyz_0[i] = g_0_zzzzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxyyz_0[i] * pb_x + g_0_xzzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxyzz_0[i] = g_0_zzzzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxyzz_0[i] * pb_x + g_0_xzzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxxzzz_0[i] = g_0_zzzzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_xzzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxxzzz_0[i] * pb_x + g_0_xzzzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxyyyy_0[i] = 4.0 * g_0_xxzzz_0_xxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xxyyyy_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_xxyyyy_0[i] * pb_z + g_0_xxzzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xxyyyz_0[i] = g_0_zzzzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxyyyz_0[i] * pb_x + g_0_xzzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxyyzz_0[i] = g_0_zzzzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxyyzz_0[i] * pb_x + g_0_xzzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxyzzz_0[i] = g_0_zzzzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxyzzz_0[i] * pb_x + g_0_xzzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xxzzzz_0[i] = g_0_zzzzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_xzzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xzzzzz_0_xxzzzz_0[i] * pb_x + g_0_xzzzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyyyyy_0[i] = 4.0 * g_0_xxzzz_0_xyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxzzz_0_xyyyyy_1[i] * fti_ab_0 +
                                    g_0_xxzzzz_0_xyyyyy_0[i] * pb_z + g_0_xxzzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxzzzzz_0_xyyyyz_0[i] = g_0_zzzzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyyyz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyyz_1[i] * fi_abcd_0 +
                                    g_0_xzzzzz_0_xyyyyz_0[i] * pb_x + g_0_xzzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyyyzz_0[i] = g_0_zzzzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyyzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyzz_1[i] * fi_abcd_0 +
                                    g_0_xzzzzz_0_xyyyzz_0[i] * pb_x + g_0_xzzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyyzzz_0[i] = g_0_zzzzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyzzz_1[i] * fi_abcd_0 +
                                    g_0_xzzzzz_0_xyyzzz_0[i] * pb_x + g_0_xzzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xyzzzz_0[i] = g_0_zzzzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yzzzz_1[i] * fi_abcd_0 +
                                    g_0_xzzzzz_0_xyzzzz_0[i] * pb_x + g_0_xzzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_xzzzzz_0[i] = g_0_zzzzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_zzzzz_1[i] * fi_abcd_0 +
                                    g_0_xzzzzz_0_xzzzzz_0[i] * pb_x + g_0_xzzzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyyyy_0[i] = g_0_zzzzz_0_yyyyyy_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyyyy_0[i] * pb_x +
                                    g_0_xzzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyyyz_0[i] = g_0_zzzzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyyyz_0[i] * pb_x +
                                    g_0_xzzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyyzz_0[i] = g_0_zzzzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyyzz_0[i] * pb_x +
                                    g_0_xzzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyyzzz_0[i] = g_0_zzzzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyyzzz_0[i] * pb_x +
                                    g_0_xzzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yyzzzz_0[i] = g_0_zzzzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yyzzzz_0[i] * pb_x +
                                    g_0_xzzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_yzzzzz_0[i] = g_0_zzzzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_yzzzzz_0[i] * pb_x +
                                    g_0_xzzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzz_0_zzzzzz_0[i] = g_0_zzzzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xzzzzz_0_zzzzzz_0[i] * pb_x +
                                    g_0_xzzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 588-616 components of targeted buffer : SKSI

    auto g_0_xyyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 588);

    auto g_0_xyyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 589);

    auto g_0_xyyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 590);

    auto g_0_xyyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 591);

    auto g_0_xyyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 592);

    auto g_0_xyyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 593);

    auto g_0_xyyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 594);

    auto g_0_xyyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 595);

    auto g_0_xyyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 596);

    auto g_0_xyyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 597);

    auto g_0_xyyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 598);

    auto g_0_xyyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 599);

    auto g_0_xyyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 600);

    auto g_0_xyyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 601);

    auto g_0_xyyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 602);

    auto g_0_xyyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 603);

    auto g_0_xyyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 604);

    auto g_0_xyyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 605);

    auto g_0_xyyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 606);

    auto g_0_xyyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 607);

    auto g_0_xyyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 608);

    auto g_0_xyyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 609);

    auto g_0_xyyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 610);

    auto g_0_xyyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 611);

    auto g_0_xyyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 612);

    auto g_0_xyyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 613);

    auto g_0_xyyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 614);

    auto g_0_xyyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 615);

#pragma omp simd aligned(g_0_xyyyyyy_0_xxxxxx_0,     \
                             g_0_xyyyyyy_0_xxxxxy_0, \
                             g_0_xyyyyyy_0_xxxxxz_0, \
                             g_0_xyyyyyy_0_xxxxyy_0, \
                             g_0_xyyyyyy_0_xxxxyz_0, \
                             g_0_xyyyyyy_0_xxxxzz_0, \
                             g_0_xyyyyyy_0_xxxyyy_0, \
                             g_0_xyyyyyy_0_xxxyyz_0, \
                             g_0_xyyyyyy_0_xxxyzz_0, \
                             g_0_xyyyyyy_0_xxxzzz_0, \
                             g_0_xyyyyyy_0_xxyyyy_0, \
                             g_0_xyyyyyy_0_xxyyyz_0, \
                             g_0_xyyyyyy_0_xxyyzz_0, \
                             g_0_xyyyyyy_0_xxyzzz_0, \
                             g_0_xyyyyyy_0_xxzzzz_0, \
                             g_0_xyyyyyy_0_xyyyyy_0, \
                             g_0_xyyyyyy_0_xyyyyz_0, \
                             g_0_xyyyyyy_0_xyyyzz_0, \
                             g_0_xyyyyyy_0_xyyzzz_0, \
                             g_0_xyyyyyy_0_xyzzzz_0, \
                             g_0_xyyyyyy_0_xzzzzz_0, \
                             g_0_xyyyyyy_0_yyyyyy_0, \
                             g_0_xyyyyyy_0_yyyyyz_0, \
                             g_0_xyyyyyy_0_yyyyzz_0, \
                             g_0_xyyyyyy_0_yyyzzz_0, \
                             g_0_xyyyyyy_0_yyzzzz_0, \
                             g_0_xyyyyyy_0_yzzzzz_0, \
                             g_0_xyyyyyy_0_zzzzzz_0, \
                             g_0_yyyyyy_0_xxxxx_1,   \
                             g_0_yyyyyy_0_xxxxxx_0,  \
                             g_0_yyyyyy_0_xxxxxx_1,  \
                             g_0_yyyyyy_0_xxxxxy_0,  \
                             g_0_yyyyyy_0_xxxxxy_1,  \
                             g_0_yyyyyy_0_xxxxxz_0,  \
                             g_0_yyyyyy_0_xxxxxz_1,  \
                             g_0_yyyyyy_0_xxxxy_1,   \
                             g_0_yyyyyy_0_xxxxyy_0,  \
                             g_0_yyyyyy_0_xxxxyy_1,  \
                             g_0_yyyyyy_0_xxxxyz_0,  \
                             g_0_yyyyyy_0_xxxxyz_1,  \
                             g_0_yyyyyy_0_xxxxz_1,   \
                             g_0_yyyyyy_0_xxxxzz_0,  \
                             g_0_yyyyyy_0_xxxxzz_1,  \
                             g_0_yyyyyy_0_xxxyy_1,   \
                             g_0_yyyyyy_0_xxxyyy_0,  \
                             g_0_yyyyyy_0_xxxyyy_1,  \
                             g_0_yyyyyy_0_xxxyyz_0,  \
                             g_0_yyyyyy_0_xxxyyz_1,  \
                             g_0_yyyyyy_0_xxxyz_1,   \
                             g_0_yyyyyy_0_xxxyzz_0,  \
                             g_0_yyyyyy_0_xxxyzz_1,  \
                             g_0_yyyyyy_0_xxxzz_1,   \
                             g_0_yyyyyy_0_xxxzzz_0,  \
                             g_0_yyyyyy_0_xxxzzz_1,  \
                             g_0_yyyyyy_0_xxyyy_1,   \
                             g_0_yyyyyy_0_xxyyyy_0,  \
                             g_0_yyyyyy_0_xxyyyy_1,  \
                             g_0_yyyyyy_0_xxyyyz_0,  \
                             g_0_yyyyyy_0_xxyyyz_1,  \
                             g_0_yyyyyy_0_xxyyz_1,   \
                             g_0_yyyyyy_0_xxyyzz_0,  \
                             g_0_yyyyyy_0_xxyyzz_1,  \
                             g_0_yyyyyy_0_xxyzz_1,   \
                             g_0_yyyyyy_0_xxyzzz_0,  \
                             g_0_yyyyyy_0_xxyzzz_1,  \
                             g_0_yyyyyy_0_xxzzz_1,   \
                             g_0_yyyyyy_0_xxzzzz_0,  \
                             g_0_yyyyyy_0_xxzzzz_1,  \
                             g_0_yyyyyy_0_xyyyy_1,   \
                             g_0_yyyyyy_0_xyyyyy_0,  \
                             g_0_yyyyyy_0_xyyyyy_1,  \
                             g_0_yyyyyy_0_xyyyyz_0,  \
                             g_0_yyyyyy_0_xyyyyz_1,  \
                             g_0_yyyyyy_0_xyyyz_1,   \
                             g_0_yyyyyy_0_xyyyzz_0,  \
                             g_0_yyyyyy_0_xyyyzz_1,  \
                             g_0_yyyyyy_0_xyyzz_1,   \
                             g_0_yyyyyy_0_xyyzzz_0,  \
                             g_0_yyyyyy_0_xyyzzz_1,  \
                             g_0_yyyyyy_0_xyzzz_1,   \
                             g_0_yyyyyy_0_xyzzzz_0,  \
                             g_0_yyyyyy_0_xyzzzz_1,  \
                             g_0_yyyyyy_0_xzzzz_1,   \
                             g_0_yyyyyy_0_xzzzzz_0,  \
                             g_0_yyyyyy_0_xzzzzz_1,  \
                             g_0_yyyyyy_0_yyyyy_1,   \
                             g_0_yyyyyy_0_yyyyyy_0,  \
                             g_0_yyyyyy_0_yyyyyy_1,  \
                             g_0_yyyyyy_0_yyyyyz_0,  \
                             g_0_yyyyyy_0_yyyyyz_1,  \
                             g_0_yyyyyy_0_yyyyz_1,   \
                             g_0_yyyyyy_0_yyyyzz_0,  \
                             g_0_yyyyyy_0_yyyyzz_1,  \
                             g_0_yyyyyy_0_yyyzz_1,   \
                             g_0_yyyyyy_0_yyyzzz_0,  \
                             g_0_yyyyyy_0_yyyzzz_1,  \
                             g_0_yyyyyy_0_yyzzz_1,   \
                             g_0_yyyyyy_0_yyzzzz_0,  \
                             g_0_yyyyyy_0_yyzzzz_1,  \
                             g_0_yyyyyy_0_yzzzz_1,   \
                             g_0_yyyyyy_0_yzzzzz_0,  \
                             g_0_yyyyyy_0_yzzzzz_1,  \
                             g_0_yyyyyy_0_zzzzz_1,   \
                             g_0_yyyyyy_0_zzzzzz_0,  \
                             g_0_yyyyyy_0_zzzzzz_1,  \
                             wp_x,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyy_0_xxxxxx_0[i] = 6.0 * g_0_yyyyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxx_0[i] * pb_x + g_0_yyyyyy_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxxy_0[i] = 5.0 * g_0_yyyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxy_0[i] * pb_x + g_0_yyyyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxxz_0[i] = 5.0 * g_0_yyyyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxz_0[i] * pb_x + g_0_yyyyyy_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxyy_0[i] = 4.0 * g_0_yyyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyy_0[i] * pb_x + g_0_yyyyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxyz_0[i] = 4.0 * g_0_yyyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyz_0[i] * pb_x + g_0_yyyyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxxzz_0[i] = 4.0 * g_0_yyyyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxzz_0[i] * pb_x + g_0_yyyyyy_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxyyy_0[i] = 3.0 * g_0_yyyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyy_0[i] * pb_x + g_0_yyyyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxyyz_0[i] = 3.0 * g_0_yyyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyz_0[i] * pb_x + g_0_yyyyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxyzz_0[i] = 3.0 * g_0_yyyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyzz_0[i] * pb_x + g_0_yyyyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxxzzz_0[i] = 3.0 * g_0_yyyyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxzzz_0[i] * pb_x + g_0_yyyyyy_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyyyy_0[i] = 2.0 * g_0_yyyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyy_0[i] * pb_x + g_0_yyyyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyyyz_0[i] = 2.0 * g_0_yyyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyz_0[i] * pb_x + g_0_yyyyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyyzz_0[i] = 2.0 * g_0_yyyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyzz_0[i] * pb_x + g_0_yyyyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxyzzz_0[i] = 2.0 * g_0_yyyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyzzz_0[i] * pb_x + g_0_yyyyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xxzzzz_0[i] = 2.0 * g_0_yyyyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxzzzz_0[i] * pb_x + g_0_yyyyyy_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyyyy_0[i] = g_0_yyyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyy_0[i] * pb_x + g_0_yyyyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyyyz_0[i] = g_0_yyyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyz_0[i] * pb_x + g_0_yyyyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyyzz_0[i] = g_0_yyyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyzz_0[i] * pb_x + g_0_yyyyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyyzzz_0[i] = g_0_yyyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyzzz_0[i] * pb_x + g_0_yyyyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xyzzzz_0[i] = g_0_yyyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyzzzz_0[i] * pb_x + g_0_yyyyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_xzzzzz_0[i] = g_0_yyyyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xzzzzz_0[i] * pb_x + g_0_yyyyyy_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyyyy_0[i] = g_0_yyyyyy_0_yyyyyy_0[i] * pb_x + g_0_yyyyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyyyz_0[i] = g_0_yyyyyy_0_yyyyyz_0[i] * pb_x + g_0_yyyyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyyzz_0[i] = g_0_yyyyyy_0_yyyyzz_0[i] * pb_x + g_0_yyyyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyyzzz_0[i] = g_0_yyyyyy_0_yyyzzz_0[i] * pb_x + g_0_yyyyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yyzzzz_0[i] = g_0_yyyyyy_0_yyzzzz_0[i] * pb_x + g_0_yyyyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_yzzzzz_0[i] = g_0_yyyyyy_0_yzzzzz_0[i] * pb_x + g_0_yyyyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyy_0_zzzzzz_0[i] = g_0_yyyyyy_0_zzzzzz_0[i] * pb_x + g_0_yyyyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 616-644 components of targeted buffer : SKSI

    auto g_0_xyyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 616);

    auto g_0_xyyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 617);

    auto g_0_xyyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 618);

    auto g_0_xyyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 619);

    auto g_0_xyyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 620);

    auto g_0_xyyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 621);

    auto g_0_xyyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 622);

    auto g_0_xyyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 623);

    auto g_0_xyyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 624);

    auto g_0_xyyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 625);

    auto g_0_xyyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 626);

    auto g_0_xyyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 627);

    auto g_0_xyyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 628);

    auto g_0_xyyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 629);

    auto g_0_xyyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 630);

    auto g_0_xyyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 631);

    auto g_0_xyyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 632);

    auto g_0_xyyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 633);

    auto g_0_xyyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 634);

    auto g_0_xyyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 635);

    auto g_0_xyyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 636);

    auto g_0_xyyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 637);

    auto g_0_xyyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 638);

    auto g_0_xyyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 639);

    auto g_0_xyyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 640);

    auto g_0_xyyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 641);

    auto g_0_xyyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 642);

    auto g_0_xyyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 643);

#pragma omp simd aligned(g_0_xyyyyy_0_xxxxxx_0,      \
                             g_0_xyyyyy_0_xxxxxx_1,  \
                             g_0_xyyyyy_0_xxxxxy_0,  \
                             g_0_xyyyyy_0_xxxxxy_1,  \
                             g_0_xyyyyy_0_xxxxyy_0,  \
                             g_0_xyyyyy_0_xxxxyy_1,  \
                             g_0_xyyyyy_0_xxxyyy_0,  \
                             g_0_xyyyyy_0_xxxyyy_1,  \
                             g_0_xyyyyy_0_xxyyyy_0,  \
                             g_0_xyyyyy_0_xxyyyy_1,  \
                             g_0_xyyyyy_0_xyyyyy_0,  \
                             g_0_xyyyyy_0_xyyyyy_1,  \
                             g_0_xyyyyyz_0_xxxxxx_0, \
                             g_0_xyyyyyz_0_xxxxxy_0, \
                             g_0_xyyyyyz_0_xxxxxz_0, \
                             g_0_xyyyyyz_0_xxxxyy_0, \
                             g_0_xyyyyyz_0_xxxxyz_0, \
                             g_0_xyyyyyz_0_xxxxzz_0, \
                             g_0_xyyyyyz_0_xxxyyy_0, \
                             g_0_xyyyyyz_0_xxxyyz_0, \
                             g_0_xyyyyyz_0_xxxyzz_0, \
                             g_0_xyyyyyz_0_xxxzzz_0, \
                             g_0_xyyyyyz_0_xxyyyy_0, \
                             g_0_xyyyyyz_0_xxyyyz_0, \
                             g_0_xyyyyyz_0_xxyyzz_0, \
                             g_0_xyyyyyz_0_xxyzzz_0, \
                             g_0_xyyyyyz_0_xxzzzz_0, \
                             g_0_xyyyyyz_0_xyyyyy_0, \
                             g_0_xyyyyyz_0_xyyyyz_0, \
                             g_0_xyyyyyz_0_xyyyzz_0, \
                             g_0_xyyyyyz_0_xyyzzz_0, \
                             g_0_xyyyyyz_0_xyzzzz_0, \
                             g_0_xyyyyyz_0_xzzzzz_0, \
                             g_0_xyyyyyz_0_yyyyyy_0, \
                             g_0_xyyyyyz_0_yyyyyz_0, \
                             g_0_xyyyyyz_0_yyyyzz_0, \
                             g_0_xyyyyyz_0_yyyzzz_0, \
                             g_0_xyyyyyz_0_yyzzzz_0, \
                             g_0_xyyyyyz_0_yzzzzz_0, \
                             g_0_xyyyyyz_0_zzzzzz_0, \
                             g_0_yyyyyz_0_xxxxxz_0,  \
                             g_0_yyyyyz_0_xxxxxz_1,  \
                             g_0_yyyyyz_0_xxxxyz_0,  \
                             g_0_yyyyyz_0_xxxxyz_1,  \
                             g_0_yyyyyz_0_xxxxz_1,   \
                             g_0_yyyyyz_0_xxxxzz_0,  \
                             g_0_yyyyyz_0_xxxxzz_1,  \
                             g_0_yyyyyz_0_xxxyyz_0,  \
                             g_0_yyyyyz_0_xxxyyz_1,  \
                             g_0_yyyyyz_0_xxxyz_1,   \
                             g_0_yyyyyz_0_xxxyzz_0,  \
                             g_0_yyyyyz_0_xxxyzz_1,  \
                             g_0_yyyyyz_0_xxxzz_1,   \
                             g_0_yyyyyz_0_xxxzzz_0,  \
                             g_0_yyyyyz_0_xxxzzz_1,  \
                             g_0_yyyyyz_0_xxyyyz_0,  \
                             g_0_yyyyyz_0_xxyyyz_1,  \
                             g_0_yyyyyz_0_xxyyz_1,   \
                             g_0_yyyyyz_0_xxyyzz_0,  \
                             g_0_yyyyyz_0_xxyyzz_1,  \
                             g_0_yyyyyz_0_xxyzz_1,   \
                             g_0_yyyyyz_0_xxyzzz_0,  \
                             g_0_yyyyyz_0_xxyzzz_1,  \
                             g_0_yyyyyz_0_xxzzz_1,   \
                             g_0_yyyyyz_0_xxzzzz_0,  \
                             g_0_yyyyyz_0_xxzzzz_1,  \
                             g_0_yyyyyz_0_xyyyyz_0,  \
                             g_0_yyyyyz_0_xyyyyz_1,  \
                             g_0_yyyyyz_0_xyyyz_1,   \
                             g_0_yyyyyz_0_xyyyzz_0,  \
                             g_0_yyyyyz_0_xyyyzz_1,  \
                             g_0_yyyyyz_0_xyyzz_1,   \
                             g_0_yyyyyz_0_xyyzzz_0,  \
                             g_0_yyyyyz_0_xyyzzz_1,  \
                             g_0_yyyyyz_0_xyzzz_1,   \
                             g_0_yyyyyz_0_xyzzzz_0,  \
                             g_0_yyyyyz_0_xyzzzz_1,  \
                             g_0_yyyyyz_0_xzzzz_1,   \
                             g_0_yyyyyz_0_xzzzzz_0,  \
                             g_0_yyyyyz_0_xzzzzz_1,  \
                             g_0_yyyyyz_0_yyyyyy_0,  \
                             g_0_yyyyyz_0_yyyyyy_1,  \
                             g_0_yyyyyz_0_yyyyyz_0,  \
                             g_0_yyyyyz_0_yyyyyz_1,  \
                             g_0_yyyyyz_0_yyyyz_1,   \
                             g_0_yyyyyz_0_yyyyzz_0,  \
                             g_0_yyyyyz_0_yyyyzz_1,  \
                             g_0_yyyyyz_0_yyyzz_1,   \
                             g_0_yyyyyz_0_yyyzzz_0,  \
                             g_0_yyyyyz_0_yyyzzz_1,  \
                             g_0_yyyyyz_0_yyzzz_1,   \
                             g_0_yyyyyz_0_yyzzzz_0,  \
                             g_0_yyyyyz_0_yyzzzz_1,  \
                             g_0_yyyyyz_0_yzzzz_1,   \
                             g_0_yyyyyz_0_yzzzzz_0,  \
                             g_0_yyyyyz_0_yzzzzz_1,  \
                             g_0_yyyyyz_0_zzzzz_1,   \
                             g_0_yyyyyz_0_zzzzzz_0,  \
                             g_0_yyyyyz_0_zzzzzz_1,  \
                             wp_x,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyz_0_xxxxxx_0[i] = g_0_xyyyyy_0_xxxxxx_0[i] * pb_z + g_0_xyyyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxxxy_0[i] = g_0_xyyyyy_0_xxxxxy_0[i] * pb_z + g_0_xyyyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxxxz_0[i] = 5.0 * g_0_yyyyyz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxxxz_0[i] * pb_x + g_0_yyyyyz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxxyy_0[i] = g_0_xyyyyy_0_xxxxyy_0[i] * pb_z + g_0_xyyyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxxyz_0[i] = 4.0 * g_0_yyyyyz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxxyz_0[i] * pb_x + g_0_yyyyyz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxxzz_0[i] = 4.0 * g_0_yyyyyz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxxzz_0[i] * pb_x + g_0_yyyyyz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxyyy_0[i] = g_0_xyyyyy_0_xxxyyy_0[i] * pb_z + g_0_xyyyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxxyyz_0[i] = 3.0 * g_0_yyyyyz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxyyz_0[i] * pb_x + g_0_yyyyyz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxyzz_0[i] = 3.0 * g_0_yyyyyz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxyzz_0[i] * pb_x + g_0_yyyyyz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxxzzz_0[i] = 3.0 * g_0_yyyyyz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxxzzz_0[i] * pb_x + g_0_yyyyyz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxyyyy_0[i] = g_0_xyyyyy_0_xxyyyy_0[i] * pb_z + g_0_xyyyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xxyyyz_0[i] = 2.0 * g_0_yyyyyz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxyyyz_0[i] * pb_x + g_0_yyyyyz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxyyzz_0[i] = 2.0 * g_0_yyyyyz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxyyzz_0[i] * pb_x + g_0_yyyyyz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxyzzz_0[i] = 2.0 * g_0_yyyyyz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxyzzz_0[i] * pb_x + g_0_yyyyyz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xxzzzz_0[i] = 2.0 * g_0_yyyyyz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xxzzzz_0[i] * pb_x + g_0_yyyyyz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyyyyy_0[i] = g_0_xyyyyy_0_xyyyyy_0[i] * pb_z + g_0_xyyyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xyyyyyz_0_xyyyyz_0[i] = g_0_yyyyyz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyyyyz_0[i] * pb_x + g_0_yyyyyz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyyyzz_0[i] = g_0_yyyyyz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyyyzz_0[i] * pb_x + g_0_yyyyyz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyyzzz_0[i] = g_0_yyyyyz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyyzzz_0[i] * pb_x + g_0_yyyyyz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xyzzzz_0[i] = g_0_yyyyyz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xyzzzz_0[i] * pb_x + g_0_yyyyyz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_xzzzzz_0[i] = g_0_yyyyyz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyyz_0_xzzzzz_0[i] * pb_x + g_0_yyyyyz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyyyy_0[i] = g_0_yyyyyz_0_yyyyyy_0[i] * pb_x + g_0_yyyyyz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyyyz_0[i] = g_0_yyyyyz_0_yyyyyz_0[i] * pb_x + g_0_yyyyyz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyyzz_0[i] = g_0_yyyyyz_0_yyyyzz_0[i] * pb_x + g_0_yyyyyz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyyzzz_0[i] = g_0_yyyyyz_0_yyyzzz_0[i] * pb_x + g_0_yyyyyz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yyzzzz_0[i] = g_0_yyyyyz_0_yyzzzz_0[i] * pb_x + g_0_yyyyyz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_yzzzzz_0[i] = g_0_yyyyyz_0_yzzzzz_0[i] * pb_x + g_0_yyyyyz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyz_0_zzzzzz_0[i] = g_0_yyyyyz_0_zzzzzz_0[i] * pb_x + g_0_yyyyyz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 644-672 components of targeted buffer : SKSI

    auto g_0_xyyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 644);

    auto g_0_xyyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 645);

    auto g_0_xyyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 646);

    auto g_0_xyyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 647);

    auto g_0_xyyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 648);

    auto g_0_xyyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 649);

    auto g_0_xyyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 650);

    auto g_0_xyyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 651);

    auto g_0_xyyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 652);

    auto g_0_xyyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 653);

    auto g_0_xyyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 654);

    auto g_0_xyyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 655);

    auto g_0_xyyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 656);

    auto g_0_xyyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 657);

    auto g_0_xyyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 658);

    auto g_0_xyyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 659);

    auto g_0_xyyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 660);

    auto g_0_xyyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 661);

    auto g_0_xyyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 662);

    auto g_0_xyyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 663);

    auto g_0_xyyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 664);

    auto g_0_xyyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 665);

    auto g_0_xyyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 666);

    auto g_0_xyyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 667);

    auto g_0_xyyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 668);

    auto g_0_xyyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 669);

    auto g_0_xyyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 670);

    auto g_0_xyyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 671);

#pragma omp simd aligned(g_0_xyyyyzz_0_xxxxxx_0,     \
                             g_0_xyyyyzz_0_xxxxxy_0, \
                             g_0_xyyyyzz_0_xxxxxz_0, \
                             g_0_xyyyyzz_0_xxxxyy_0, \
                             g_0_xyyyyzz_0_xxxxyz_0, \
                             g_0_xyyyyzz_0_xxxxzz_0, \
                             g_0_xyyyyzz_0_xxxyyy_0, \
                             g_0_xyyyyzz_0_xxxyyz_0, \
                             g_0_xyyyyzz_0_xxxyzz_0, \
                             g_0_xyyyyzz_0_xxxzzz_0, \
                             g_0_xyyyyzz_0_xxyyyy_0, \
                             g_0_xyyyyzz_0_xxyyyz_0, \
                             g_0_xyyyyzz_0_xxyyzz_0, \
                             g_0_xyyyyzz_0_xxyzzz_0, \
                             g_0_xyyyyzz_0_xxzzzz_0, \
                             g_0_xyyyyzz_0_xyyyyy_0, \
                             g_0_xyyyyzz_0_xyyyyz_0, \
                             g_0_xyyyyzz_0_xyyyzz_0, \
                             g_0_xyyyyzz_0_xyyzzz_0, \
                             g_0_xyyyyzz_0_xyzzzz_0, \
                             g_0_xyyyyzz_0_xzzzzz_0, \
                             g_0_xyyyyzz_0_yyyyyy_0, \
                             g_0_xyyyyzz_0_yyyyyz_0, \
                             g_0_xyyyyzz_0_yyyyzz_0, \
                             g_0_xyyyyzz_0_yyyzzz_0, \
                             g_0_xyyyyzz_0_yyzzzz_0, \
                             g_0_xyyyyzz_0_yzzzzz_0, \
                             g_0_xyyyyzz_0_zzzzzz_0, \
                             g_0_yyyyzz_0_xxxxx_1,   \
                             g_0_yyyyzz_0_xxxxxx_0,  \
                             g_0_yyyyzz_0_xxxxxx_1,  \
                             g_0_yyyyzz_0_xxxxxy_0,  \
                             g_0_yyyyzz_0_xxxxxy_1,  \
                             g_0_yyyyzz_0_xxxxxz_0,  \
                             g_0_yyyyzz_0_xxxxxz_1,  \
                             g_0_yyyyzz_0_xxxxy_1,   \
                             g_0_yyyyzz_0_xxxxyy_0,  \
                             g_0_yyyyzz_0_xxxxyy_1,  \
                             g_0_yyyyzz_0_xxxxyz_0,  \
                             g_0_yyyyzz_0_xxxxyz_1,  \
                             g_0_yyyyzz_0_xxxxz_1,   \
                             g_0_yyyyzz_0_xxxxzz_0,  \
                             g_0_yyyyzz_0_xxxxzz_1,  \
                             g_0_yyyyzz_0_xxxyy_1,   \
                             g_0_yyyyzz_0_xxxyyy_0,  \
                             g_0_yyyyzz_0_xxxyyy_1,  \
                             g_0_yyyyzz_0_xxxyyz_0,  \
                             g_0_yyyyzz_0_xxxyyz_1,  \
                             g_0_yyyyzz_0_xxxyz_1,   \
                             g_0_yyyyzz_0_xxxyzz_0,  \
                             g_0_yyyyzz_0_xxxyzz_1,  \
                             g_0_yyyyzz_0_xxxzz_1,   \
                             g_0_yyyyzz_0_xxxzzz_0,  \
                             g_0_yyyyzz_0_xxxzzz_1,  \
                             g_0_yyyyzz_0_xxyyy_1,   \
                             g_0_yyyyzz_0_xxyyyy_0,  \
                             g_0_yyyyzz_0_xxyyyy_1,  \
                             g_0_yyyyzz_0_xxyyyz_0,  \
                             g_0_yyyyzz_0_xxyyyz_1,  \
                             g_0_yyyyzz_0_xxyyz_1,   \
                             g_0_yyyyzz_0_xxyyzz_0,  \
                             g_0_yyyyzz_0_xxyyzz_1,  \
                             g_0_yyyyzz_0_xxyzz_1,   \
                             g_0_yyyyzz_0_xxyzzz_0,  \
                             g_0_yyyyzz_0_xxyzzz_1,  \
                             g_0_yyyyzz_0_xxzzz_1,   \
                             g_0_yyyyzz_0_xxzzzz_0,  \
                             g_0_yyyyzz_0_xxzzzz_1,  \
                             g_0_yyyyzz_0_xyyyy_1,   \
                             g_0_yyyyzz_0_xyyyyy_0,  \
                             g_0_yyyyzz_0_xyyyyy_1,  \
                             g_0_yyyyzz_0_xyyyyz_0,  \
                             g_0_yyyyzz_0_xyyyyz_1,  \
                             g_0_yyyyzz_0_xyyyz_1,   \
                             g_0_yyyyzz_0_xyyyzz_0,  \
                             g_0_yyyyzz_0_xyyyzz_1,  \
                             g_0_yyyyzz_0_xyyzz_1,   \
                             g_0_yyyyzz_0_xyyzzz_0,  \
                             g_0_yyyyzz_0_xyyzzz_1,  \
                             g_0_yyyyzz_0_xyzzz_1,   \
                             g_0_yyyyzz_0_xyzzzz_0,  \
                             g_0_yyyyzz_0_xyzzzz_1,  \
                             g_0_yyyyzz_0_xzzzz_1,   \
                             g_0_yyyyzz_0_xzzzzz_0,  \
                             g_0_yyyyzz_0_xzzzzz_1,  \
                             g_0_yyyyzz_0_yyyyy_1,   \
                             g_0_yyyyzz_0_yyyyyy_0,  \
                             g_0_yyyyzz_0_yyyyyy_1,  \
                             g_0_yyyyzz_0_yyyyyz_0,  \
                             g_0_yyyyzz_0_yyyyyz_1,  \
                             g_0_yyyyzz_0_yyyyz_1,   \
                             g_0_yyyyzz_0_yyyyzz_0,  \
                             g_0_yyyyzz_0_yyyyzz_1,  \
                             g_0_yyyyzz_0_yyyzz_1,   \
                             g_0_yyyyzz_0_yyyzzz_0,  \
                             g_0_yyyyzz_0_yyyzzz_1,  \
                             g_0_yyyyzz_0_yyzzz_1,   \
                             g_0_yyyyzz_0_yyzzzz_0,  \
                             g_0_yyyyzz_0_yyzzzz_1,  \
                             g_0_yyyyzz_0_yzzzz_1,   \
                             g_0_yyyyzz_0_yzzzzz_0,  \
                             g_0_yyyyzz_0_yzzzzz_1,  \
                             g_0_yyyyzz_0_zzzzz_1,   \
                             g_0_yyyyzz_0_zzzzzz_0,  \
                             g_0_yyyyzz_0_zzzzzz_1,  \
                             wp_x,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzz_0_xxxxxx_0[i] = 6.0 * g_0_yyyyzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxxx_0[i] * pb_x + g_0_yyyyzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxxy_0[i] = 5.0 * g_0_yyyyzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxxy_0[i] * pb_x + g_0_yyyyzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxxz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxxz_0[i] * pb_x + g_0_yyyyzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxyy_0[i] = 4.0 * g_0_yyyyzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxyy_0[i] * pb_x + g_0_yyyyzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxyz_0[i] = 4.0 * g_0_yyyyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxyz_0[i] * pb_x + g_0_yyyyzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxxzz_0[i] = 4.0 * g_0_yyyyzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxzz_0[i] * pb_x + g_0_yyyyzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxyyy_0[i] = 3.0 * g_0_yyyyzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyyy_0[i] * pb_x + g_0_yyyyzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxyyz_0[i] = 3.0 * g_0_yyyyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyyz_0[i] * pb_x + g_0_yyyyzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxyzz_0[i] = 3.0 * g_0_yyyyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyzz_0[i] * pb_x + g_0_yyyyzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxxzzz_0[i] = 3.0 * g_0_yyyyzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxzzz_0[i] * pb_x + g_0_yyyyzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyyyy_0[i] = 2.0 * g_0_yyyyzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyyy_0[i] * pb_x + g_0_yyyyzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyyyz_0[i] = 2.0 * g_0_yyyyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyyz_0[i] * pb_x + g_0_yyyyzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyyzz_0[i] = 2.0 * g_0_yyyyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyzz_0[i] * pb_x + g_0_yyyyzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxyzzz_0[i] = 2.0 * g_0_yyyyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyzzz_0[i] * pb_x + g_0_yyyyzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xxzzzz_0[i] = 2.0 * g_0_yyyyzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxzzzz_0[i] * pb_x + g_0_yyyyzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyyyy_0[i] = g_0_yyyyzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyyy_0[i] * pb_x + g_0_yyyyzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyyyz_0[i] = g_0_yyyyzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyyz_0[i] * pb_x + g_0_yyyyzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyyzz_0[i] = g_0_yyyyzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyzz_0[i] * pb_x + g_0_yyyyzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyyzzz_0[i] = g_0_yyyyzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyzzz_0[i] * pb_x + g_0_yyyyzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xyzzzz_0[i] = g_0_yyyyzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyzzzz_0[i] * pb_x + g_0_yyyyzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_xzzzzz_0[i] = g_0_yyyyzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xzzzzz_0[i] * pb_x + g_0_yyyyzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyyyy_0[i] = g_0_yyyyzz_0_yyyyyy_0[i] * pb_x + g_0_yyyyzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyyyz_0[i] = g_0_yyyyzz_0_yyyyyz_0[i] * pb_x + g_0_yyyyzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyyzz_0[i] = g_0_yyyyzz_0_yyyyzz_0[i] * pb_x + g_0_yyyyzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyyzzz_0[i] = g_0_yyyyzz_0_yyyzzz_0[i] * pb_x + g_0_yyyyzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yyzzzz_0[i] = g_0_yyyyzz_0_yyzzzz_0[i] * pb_x + g_0_yyyyzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_yzzzzz_0[i] = g_0_yyyyzz_0_yzzzzz_0[i] * pb_x + g_0_yyyyzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzz_0_zzzzzz_0[i] = g_0_yyyyzz_0_zzzzzz_0[i] * pb_x + g_0_yyyyzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 672-700 components of targeted buffer : SKSI

    auto g_0_xyyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 672);

    auto g_0_xyyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 673);

    auto g_0_xyyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 674);

    auto g_0_xyyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 675);

    auto g_0_xyyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 676);

    auto g_0_xyyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 677);

    auto g_0_xyyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 678);

    auto g_0_xyyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 679);

    auto g_0_xyyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 680);

    auto g_0_xyyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 681);

    auto g_0_xyyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 682);

    auto g_0_xyyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 683);

    auto g_0_xyyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 684);

    auto g_0_xyyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 685);

    auto g_0_xyyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 686);

    auto g_0_xyyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 687);

    auto g_0_xyyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 688);

    auto g_0_xyyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 689);

    auto g_0_xyyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 690);

    auto g_0_xyyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 691);

    auto g_0_xyyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 692);

    auto g_0_xyyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 693);

    auto g_0_xyyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 694);

    auto g_0_xyyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 695);

    auto g_0_xyyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 696);

    auto g_0_xyyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 697);

    auto g_0_xyyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 698);

    auto g_0_xyyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 699);

#pragma omp simd aligned(g_0_xyyyzzz_0_xxxxxx_0,     \
                             g_0_xyyyzzz_0_xxxxxy_0, \
                             g_0_xyyyzzz_0_xxxxxz_0, \
                             g_0_xyyyzzz_0_xxxxyy_0, \
                             g_0_xyyyzzz_0_xxxxyz_0, \
                             g_0_xyyyzzz_0_xxxxzz_0, \
                             g_0_xyyyzzz_0_xxxyyy_0, \
                             g_0_xyyyzzz_0_xxxyyz_0, \
                             g_0_xyyyzzz_0_xxxyzz_0, \
                             g_0_xyyyzzz_0_xxxzzz_0, \
                             g_0_xyyyzzz_0_xxyyyy_0, \
                             g_0_xyyyzzz_0_xxyyyz_0, \
                             g_0_xyyyzzz_0_xxyyzz_0, \
                             g_0_xyyyzzz_0_xxyzzz_0, \
                             g_0_xyyyzzz_0_xxzzzz_0, \
                             g_0_xyyyzzz_0_xyyyyy_0, \
                             g_0_xyyyzzz_0_xyyyyz_0, \
                             g_0_xyyyzzz_0_xyyyzz_0, \
                             g_0_xyyyzzz_0_xyyzzz_0, \
                             g_0_xyyyzzz_0_xyzzzz_0, \
                             g_0_xyyyzzz_0_xzzzzz_0, \
                             g_0_xyyyzzz_0_yyyyyy_0, \
                             g_0_xyyyzzz_0_yyyyyz_0, \
                             g_0_xyyyzzz_0_yyyyzz_0, \
                             g_0_xyyyzzz_0_yyyzzz_0, \
                             g_0_xyyyzzz_0_yyzzzz_0, \
                             g_0_xyyyzzz_0_yzzzzz_0, \
                             g_0_xyyyzzz_0_zzzzzz_0, \
                             g_0_yyyzzz_0_xxxxx_1,   \
                             g_0_yyyzzz_0_xxxxxx_0,  \
                             g_0_yyyzzz_0_xxxxxx_1,  \
                             g_0_yyyzzz_0_xxxxxy_0,  \
                             g_0_yyyzzz_0_xxxxxy_1,  \
                             g_0_yyyzzz_0_xxxxxz_0,  \
                             g_0_yyyzzz_0_xxxxxz_1,  \
                             g_0_yyyzzz_0_xxxxy_1,   \
                             g_0_yyyzzz_0_xxxxyy_0,  \
                             g_0_yyyzzz_0_xxxxyy_1,  \
                             g_0_yyyzzz_0_xxxxyz_0,  \
                             g_0_yyyzzz_0_xxxxyz_1,  \
                             g_0_yyyzzz_0_xxxxz_1,   \
                             g_0_yyyzzz_0_xxxxzz_0,  \
                             g_0_yyyzzz_0_xxxxzz_1,  \
                             g_0_yyyzzz_0_xxxyy_1,   \
                             g_0_yyyzzz_0_xxxyyy_0,  \
                             g_0_yyyzzz_0_xxxyyy_1,  \
                             g_0_yyyzzz_0_xxxyyz_0,  \
                             g_0_yyyzzz_0_xxxyyz_1,  \
                             g_0_yyyzzz_0_xxxyz_1,   \
                             g_0_yyyzzz_0_xxxyzz_0,  \
                             g_0_yyyzzz_0_xxxyzz_1,  \
                             g_0_yyyzzz_0_xxxzz_1,   \
                             g_0_yyyzzz_0_xxxzzz_0,  \
                             g_0_yyyzzz_0_xxxzzz_1,  \
                             g_0_yyyzzz_0_xxyyy_1,   \
                             g_0_yyyzzz_0_xxyyyy_0,  \
                             g_0_yyyzzz_0_xxyyyy_1,  \
                             g_0_yyyzzz_0_xxyyyz_0,  \
                             g_0_yyyzzz_0_xxyyyz_1,  \
                             g_0_yyyzzz_0_xxyyz_1,   \
                             g_0_yyyzzz_0_xxyyzz_0,  \
                             g_0_yyyzzz_0_xxyyzz_1,  \
                             g_0_yyyzzz_0_xxyzz_1,   \
                             g_0_yyyzzz_0_xxyzzz_0,  \
                             g_0_yyyzzz_0_xxyzzz_1,  \
                             g_0_yyyzzz_0_xxzzz_1,   \
                             g_0_yyyzzz_0_xxzzzz_0,  \
                             g_0_yyyzzz_0_xxzzzz_1,  \
                             g_0_yyyzzz_0_xyyyy_1,   \
                             g_0_yyyzzz_0_xyyyyy_0,  \
                             g_0_yyyzzz_0_xyyyyy_1,  \
                             g_0_yyyzzz_0_xyyyyz_0,  \
                             g_0_yyyzzz_0_xyyyyz_1,  \
                             g_0_yyyzzz_0_xyyyz_1,   \
                             g_0_yyyzzz_0_xyyyzz_0,  \
                             g_0_yyyzzz_0_xyyyzz_1,  \
                             g_0_yyyzzz_0_xyyzz_1,   \
                             g_0_yyyzzz_0_xyyzzz_0,  \
                             g_0_yyyzzz_0_xyyzzz_1,  \
                             g_0_yyyzzz_0_xyzzz_1,   \
                             g_0_yyyzzz_0_xyzzzz_0,  \
                             g_0_yyyzzz_0_xyzzzz_1,  \
                             g_0_yyyzzz_0_xzzzz_1,   \
                             g_0_yyyzzz_0_xzzzzz_0,  \
                             g_0_yyyzzz_0_xzzzzz_1,  \
                             g_0_yyyzzz_0_yyyyy_1,   \
                             g_0_yyyzzz_0_yyyyyy_0,  \
                             g_0_yyyzzz_0_yyyyyy_1,  \
                             g_0_yyyzzz_0_yyyyyz_0,  \
                             g_0_yyyzzz_0_yyyyyz_1,  \
                             g_0_yyyzzz_0_yyyyz_1,   \
                             g_0_yyyzzz_0_yyyyzz_0,  \
                             g_0_yyyzzz_0_yyyyzz_1,  \
                             g_0_yyyzzz_0_yyyzz_1,   \
                             g_0_yyyzzz_0_yyyzzz_0,  \
                             g_0_yyyzzz_0_yyyzzz_1,  \
                             g_0_yyyzzz_0_yyzzz_1,   \
                             g_0_yyyzzz_0_yyzzzz_0,  \
                             g_0_yyyzzz_0_yyzzzz_1,  \
                             g_0_yyyzzz_0_yzzzz_1,   \
                             g_0_yyyzzz_0_yzzzzz_0,  \
                             g_0_yyyzzz_0_yzzzzz_1,  \
                             g_0_yyyzzz_0_zzzzz_1,   \
                             g_0_yyyzzz_0_zzzzzz_0,  \
                             g_0_yyyzzz_0_zzzzzz_1,  \
                             wp_x,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzz_0_xxxxxx_0[i] = 6.0 * g_0_yyyzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxxx_0[i] * pb_x + g_0_yyyzzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxxy_0[i] = 5.0 * g_0_yyyzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxxy_0[i] * pb_x + g_0_yyyzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxxz_0[i] = 5.0 * g_0_yyyzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxxz_0[i] * pb_x + g_0_yyyzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxyy_0[i] = 4.0 * g_0_yyyzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxyy_0[i] * pb_x + g_0_yyyzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxyz_0[i] = 4.0 * g_0_yyyzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxyz_0[i] * pb_x + g_0_yyyzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxxzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxzz_0[i] * pb_x + g_0_yyyzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxyyy_0[i] = 3.0 * g_0_yyyzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyyy_0[i] * pb_x + g_0_yyyzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxyyz_0[i] = 3.0 * g_0_yyyzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyyz_0[i] * pb_x + g_0_yyyzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxyzz_0[i] = 3.0 * g_0_yyyzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyzz_0[i] * pb_x + g_0_yyyzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxxzzz_0[i] = 3.0 * g_0_yyyzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxzzz_0[i] * pb_x + g_0_yyyzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyyyy_0[i] = 2.0 * g_0_yyyzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyyy_0[i] * pb_x + g_0_yyyzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyyyz_0[i] = 2.0 * g_0_yyyzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyyz_0[i] * pb_x + g_0_yyyzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyyzz_0[i] = 2.0 * g_0_yyyzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyzz_0[i] * pb_x + g_0_yyyzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxyzzz_0[i] = 2.0 * g_0_yyyzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyzzz_0[i] * pb_x + g_0_yyyzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xxzzzz_0[i] = 2.0 * g_0_yyyzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxzzzz_0[i] * pb_x + g_0_yyyzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyyyy_0[i] = g_0_yyyzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyyy_0[i] * pb_x + g_0_yyyzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyyyz_0[i] = g_0_yyyzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyyz_0[i] * pb_x + g_0_yyyzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyyzz_0[i] = g_0_yyyzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyzz_0[i] * pb_x + g_0_yyyzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyyzzz_0[i] = g_0_yyyzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyzzz_0[i] * pb_x + g_0_yyyzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xyzzzz_0[i] = g_0_yyyzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyzzzz_0[i] * pb_x + g_0_yyyzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_xzzzzz_0[i] = g_0_yyyzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xzzzzz_0[i] * pb_x + g_0_yyyzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyyyy_0[i] = g_0_yyyzzz_0_yyyyyy_0[i] * pb_x + g_0_yyyzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyyyz_0[i] = g_0_yyyzzz_0_yyyyyz_0[i] * pb_x + g_0_yyyzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyyzz_0[i] = g_0_yyyzzz_0_yyyyzz_0[i] * pb_x + g_0_yyyzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyyzzz_0[i] = g_0_yyyzzz_0_yyyzzz_0[i] * pb_x + g_0_yyyzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yyzzzz_0[i] = g_0_yyyzzz_0_yyzzzz_0[i] * pb_x + g_0_yyyzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_yzzzzz_0[i] = g_0_yyyzzz_0_yzzzzz_0[i] * pb_x + g_0_yyyzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzz_0_zzzzzz_0[i] = g_0_yyyzzz_0_zzzzzz_0[i] * pb_x + g_0_yyyzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 700-728 components of targeted buffer : SKSI

    auto g_0_xyyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 700);

    auto g_0_xyyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 701);

    auto g_0_xyyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 702);

    auto g_0_xyyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 703);

    auto g_0_xyyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 704);

    auto g_0_xyyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 705);

    auto g_0_xyyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 706);

    auto g_0_xyyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 707);

    auto g_0_xyyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 708);

    auto g_0_xyyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 709);

    auto g_0_xyyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 710);

    auto g_0_xyyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 711);

    auto g_0_xyyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 712);

    auto g_0_xyyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 713);

    auto g_0_xyyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 714);

    auto g_0_xyyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 715);

    auto g_0_xyyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 716);

    auto g_0_xyyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 717);

    auto g_0_xyyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 718);

    auto g_0_xyyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 719);

    auto g_0_xyyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 720);

    auto g_0_xyyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 721);

    auto g_0_xyyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 722);

    auto g_0_xyyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 723);

    auto g_0_xyyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 724);

    auto g_0_xyyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 725);

    auto g_0_xyyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 726);

    auto g_0_xyyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 727);

#pragma omp simd aligned(g_0_xyyzzzz_0_xxxxxx_0,     \
                             g_0_xyyzzzz_0_xxxxxy_0, \
                             g_0_xyyzzzz_0_xxxxxz_0, \
                             g_0_xyyzzzz_0_xxxxyy_0, \
                             g_0_xyyzzzz_0_xxxxyz_0, \
                             g_0_xyyzzzz_0_xxxxzz_0, \
                             g_0_xyyzzzz_0_xxxyyy_0, \
                             g_0_xyyzzzz_0_xxxyyz_0, \
                             g_0_xyyzzzz_0_xxxyzz_0, \
                             g_0_xyyzzzz_0_xxxzzz_0, \
                             g_0_xyyzzzz_0_xxyyyy_0, \
                             g_0_xyyzzzz_0_xxyyyz_0, \
                             g_0_xyyzzzz_0_xxyyzz_0, \
                             g_0_xyyzzzz_0_xxyzzz_0, \
                             g_0_xyyzzzz_0_xxzzzz_0, \
                             g_0_xyyzzzz_0_xyyyyy_0, \
                             g_0_xyyzzzz_0_xyyyyz_0, \
                             g_0_xyyzzzz_0_xyyyzz_0, \
                             g_0_xyyzzzz_0_xyyzzz_0, \
                             g_0_xyyzzzz_0_xyzzzz_0, \
                             g_0_xyyzzzz_0_xzzzzz_0, \
                             g_0_xyyzzzz_0_yyyyyy_0, \
                             g_0_xyyzzzz_0_yyyyyz_0, \
                             g_0_xyyzzzz_0_yyyyzz_0, \
                             g_0_xyyzzzz_0_yyyzzz_0, \
                             g_0_xyyzzzz_0_yyzzzz_0, \
                             g_0_xyyzzzz_0_yzzzzz_0, \
                             g_0_xyyzzzz_0_zzzzzz_0, \
                             g_0_yyzzzz_0_xxxxx_1,   \
                             g_0_yyzzzz_0_xxxxxx_0,  \
                             g_0_yyzzzz_0_xxxxxx_1,  \
                             g_0_yyzzzz_0_xxxxxy_0,  \
                             g_0_yyzzzz_0_xxxxxy_1,  \
                             g_0_yyzzzz_0_xxxxxz_0,  \
                             g_0_yyzzzz_0_xxxxxz_1,  \
                             g_0_yyzzzz_0_xxxxy_1,   \
                             g_0_yyzzzz_0_xxxxyy_0,  \
                             g_0_yyzzzz_0_xxxxyy_1,  \
                             g_0_yyzzzz_0_xxxxyz_0,  \
                             g_0_yyzzzz_0_xxxxyz_1,  \
                             g_0_yyzzzz_0_xxxxz_1,   \
                             g_0_yyzzzz_0_xxxxzz_0,  \
                             g_0_yyzzzz_0_xxxxzz_1,  \
                             g_0_yyzzzz_0_xxxyy_1,   \
                             g_0_yyzzzz_0_xxxyyy_0,  \
                             g_0_yyzzzz_0_xxxyyy_1,  \
                             g_0_yyzzzz_0_xxxyyz_0,  \
                             g_0_yyzzzz_0_xxxyyz_1,  \
                             g_0_yyzzzz_0_xxxyz_1,   \
                             g_0_yyzzzz_0_xxxyzz_0,  \
                             g_0_yyzzzz_0_xxxyzz_1,  \
                             g_0_yyzzzz_0_xxxzz_1,   \
                             g_0_yyzzzz_0_xxxzzz_0,  \
                             g_0_yyzzzz_0_xxxzzz_1,  \
                             g_0_yyzzzz_0_xxyyy_1,   \
                             g_0_yyzzzz_0_xxyyyy_0,  \
                             g_0_yyzzzz_0_xxyyyy_1,  \
                             g_0_yyzzzz_0_xxyyyz_0,  \
                             g_0_yyzzzz_0_xxyyyz_1,  \
                             g_0_yyzzzz_0_xxyyz_1,   \
                             g_0_yyzzzz_0_xxyyzz_0,  \
                             g_0_yyzzzz_0_xxyyzz_1,  \
                             g_0_yyzzzz_0_xxyzz_1,   \
                             g_0_yyzzzz_0_xxyzzz_0,  \
                             g_0_yyzzzz_0_xxyzzz_1,  \
                             g_0_yyzzzz_0_xxzzz_1,   \
                             g_0_yyzzzz_0_xxzzzz_0,  \
                             g_0_yyzzzz_0_xxzzzz_1,  \
                             g_0_yyzzzz_0_xyyyy_1,   \
                             g_0_yyzzzz_0_xyyyyy_0,  \
                             g_0_yyzzzz_0_xyyyyy_1,  \
                             g_0_yyzzzz_0_xyyyyz_0,  \
                             g_0_yyzzzz_0_xyyyyz_1,  \
                             g_0_yyzzzz_0_xyyyz_1,   \
                             g_0_yyzzzz_0_xyyyzz_0,  \
                             g_0_yyzzzz_0_xyyyzz_1,  \
                             g_0_yyzzzz_0_xyyzz_1,   \
                             g_0_yyzzzz_0_xyyzzz_0,  \
                             g_0_yyzzzz_0_xyyzzz_1,  \
                             g_0_yyzzzz_0_xyzzz_1,   \
                             g_0_yyzzzz_0_xyzzzz_0,  \
                             g_0_yyzzzz_0_xyzzzz_1,  \
                             g_0_yyzzzz_0_xzzzz_1,   \
                             g_0_yyzzzz_0_xzzzzz_0,  \
                             g_0_yyzzzz_0_xzzzzz_1,  \
                             g_0_yyzzzz_0_yyyyy_1,   \
                             g_0_yyzzzz_0_yyyyyy_0,  \
                             g_0_yyzzzz_0_yyyyyy_1,  \
                             g_0_yyzzzz_0_yyyyyz_0,  \
                             g_0_yyzzzz_0_yyyyyz_1,  \
                             g_0_yyzzzz_0_yyyyz_1,   \
                             g_0_yyzzzz_0_yyyyzz_0,  \
                             g_0_yyzzzz_0_yyyyzz_1,  \
                             g_0_yyzzzz_0_yyyzz_1,   \
                             g_0_yyzzzz_0_yyyzzz_0,  \
                             g_0_yyzzzz_0_yyyzzz_1,  \
                             g_0_yyzzzz_0_yyzzz_1,   \
                             g_0_yyzzzz_0_yyzzzz_0,  \
                             g_0_yyzzzz_0_yyzzzz_1,  \
                             g_0_yyzzzz_0_yzzzz_1,   \
                             g_0_yyzzzz_0_yzzzzz_0,  \
                             g_0_yyzzzz_0_yzzzzz_1,  \
                             g_0_yyzzzz_0_zzzzz_1,   \
                             g_0_yyzzzz_0_zzzzzz_0,  \
                             g_0_yyzzzz_0_zzzzzz_1,  \
                             wp_x,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzz_0_xxxxxx_0[i] = 6.0 * g_0_yyzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxxx_0[i] * pb_x + g_0_yyzzzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxxy_0[i] = 5.0 * g_0_yyzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxxy_0[i] * pb_x + g_0_yyzzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxxz_0[i] = 5.0 * g_0_yyzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxxz_0[i] * pb_x + g_0_yyzzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxyy_0[i] = 4.0 * g_0_yyzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxyy_0[i] * pb_x + g_0_yyzzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxyz_0[i] = 4.0 * g_0_yyzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxyz_0[i] * pb_x + g_0_yyzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxxzz_0[i] = 4.0 * g_0_yyzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxzz_0[i] * pb_x + g_0_yyzzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxyyy_0[i] = 3.0 * g_0_yyzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyyy_0[i] * pb_x + g_0_yyzzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxyyz_0[i] = 3.0 * g_0_yyzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyyz_0[i] * pb_x + g_0_yyzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxyzz_0[i] = 3.0 * g_0_yyzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyzz_0[i] * pb_x + g_0_yyzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxxzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxzzz_0[i] * pb_x + g_0_yyzzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyyyy_0[i] = 2.0 * g_0_yyzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyyy_0[i] * pb_x + g_0_yyzzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyyyz_0[i] = 2.0 * g_0_yyzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyyz_0[i] * pb_x + g_0_yyzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyyzz_0[i] = 2.0 * g_0_yyzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyzz_0[i] * pb_x + g_0_yyzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxyzzz_0[i] = 2.0 * g_0_yyzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyzzz_0[i] * pb_x + g_0_yyzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xxzzzz_0[i] = 2.0 * g_0_yyzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxzzzz_0[i] * pb_x + g_0_yyzzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyyyy_0[i] = g_0_yyzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyyy_0[i] * pb_x + g_0_yyzzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyyyz_0[i] = g_0_yyzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyyz_0[i] * pb_x + g_0_yyzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyyzz_0[i] = g_0_yyzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyzz_0[i] * pb_x + g_0_yyzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyyzzz_0[i] = g_0_yyzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyzzz_0[i] * pb_x + g_0_yyzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xyzzzz_0[i] = g_0_yyzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyzzzz_0[i] * pb_x + g_0_yyzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_xzzzzz_0[i] = g_0_yyzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xzzzzz_0[i] * pb_x + g_0_yyzzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyyyy_0[i] = g_0_yyzzzz_0_yyyyyy_0[i] * pb_x + g_0_yyzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyyyz_0[i] = g_0_yyzzzz_0_yyyyyz_0[i] * pb_x + g_0_yyzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyyzz_0[i] = g_0_yyzzzz_0_yyyyzz_0[i] * pb_x + g_0_yyzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyyzzz_0[i] = g_0_yyzzzz_0_yyyzzz_0[i] * pb_x + g_0_yyzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yyzzzz_0[i] = g_0_yyzzzz_0_yyzzzz_0[i] * pb_x + g_0_yyzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_yzzzzz_0[i] = g_0_yyzzzz_0_yzzzzz_0[i] * pb_x + g_0_yyzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzz_0_zzzzzz_0[i] = g_0_yyzzzz_0_zzzzzz_0[i] * pb_x + g_0_yyzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 728-756 components of targeted buffer : SKSI

    auto g_0_xyzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 728);

    auto g_0_xyzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 729);

    auto g_0_xyzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 730);

    auto g_0_xyzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 731);

    auto g_0_xyzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 732);

    auto g_0_xyzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 733);

    auto g_0_xyzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 734);

    auto g_0_xyzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 735);

    auto g_0_xyzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 736);

    auto g_0_xyzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 737);

    auto g_0_xyzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 738);

    auto g_0_xyzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 739);

    auto g_0_xyzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 740);

    auto g_0_xyzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 741);

    auto g_0_xyzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 742);

    auto g_0_xyzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 743);

    auto g_0_xyzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 744);

    auto g_0_xyzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 745);

    auto g_0_xyzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 746);

    auto g_0_xyzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 747);

    auto g_0_xyzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 748);

    auto g_0_xyzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 749);

    auto g_0_xyzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 750);

    auto g_0_xyzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 751);

    auto g_0_xyzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 752);

    auto g_0_xyzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 753);

    auto g_0_xyzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 754);

    auto g_0_xyzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 755);

#pragma omp simd aligned(g_0_xyzzzzz_0_xxxxxx_0,     \
                             g_0_xyzzzzz_0_xxxxxy_0, \
                             g_0_xyzzzzz_0_xxxxxz_0, \
                             g_0_xyzzzzz_0_xxxxyy_0, \
                             g_0_xyzzzzz_0_xxxxyz_0, \
                             g_0_xyzzzzz_0_xxxxzz_0, \
                             g_0_xyzzzzz_0_xxxyyy_0, \
                             g_0_xyzzzzz_0_xxxyyz_0, \
                             g_0_xyzzzzz_0_xxxyzz_0, \
                             g_0_xyzzzzz_0_xxxzzz_0, \
                             g_0_xyzzzzz_0_xxyyyy_0, \
                             g_0_xyzzzzz_0_xxyyyz_0, \
                             g_0_xyzzzzz_0_xxyyzz_0, \
                             g_0_xyzzzzz_0_xxyzzz_0, \
                             g_0_xyzzzzz_0_xxzzzz_0, \
                             g_0_xyzzzzz_0_xyyyyy_0, \
                             g_0_xyzzzzz_0_xyyyyz_0, \
                             g_0_xyzzzzz_0_xyyyzz_0, \
                             g_0_xyzzzzz_0_xyyzzz_0, \
                             g_0_xyzzzzz_0_xyzzzz_0, \
                             g_0_xyzzzzz_0_xzzzzz_0, \
                             g_0_xyzzzzz_0_yyyyyy_0, \
                             g_0_xyzzzzz_0_yyyyyz_0, \
                             g_0_xyzzzzz_0_yyyyzz_0, \
                             g_0_xyzzzzz_0_yyyzzz_0, \
                             g_0_xyzzzzz_0_yyzzzz_0, \
                             g_0_xyzzzzz_0_yzzzzz_0, \
                             g_0_xyzzzzz_0_zzzzzz_0, \
                             g_0_xzzzzz_0_xxxxxx_0,  \
                             g_0_xzzzzz_0_xxxxxx_1,  \
                             g_0_xzzzzz_0_xxxxxz_0,  \
                             g_0_xzzzzz_0_xxxxxz_1,  \
                             g_0_xzzzzz_0_xxxxzz_0,  \
                             g_0_xzzzzz_0_xxxxzz_1,  \
                             g_0_xzzzzz_0_xxxzzz_0,  \
                             g_0_xzzzzz_0_xxxzzz_1,  \
                             g_0_xzzzzz_0_xxzzzz_0,  \
                             g_0_xzzzzz_0_xxzzzz_1,  \
                             g_0_xzzzzz_0_xzzzzz_0,  \
                             g_0_xzzzzz_0_xzzzzz_1,  \
                             g_0_yzzzzz_0_xxxxxy_0,  \
                             g_0_yzzzzz_0_xxxxxy_1,  \
                             g_0_yzzzzz_0_xxxxy_1,   \
                             g_0_yzzzzz_0_xxxxyy_0,  \
                             g_0_yzzzzz_0_xxxxyy_1,  \
                             g_0_yzzzzz_0_xxxxyz_0,  \
                             g_0_yzzzzz_0_xxxxyz_1,  \
                             g_0_yzzzzz_0_xxxyy_1,   \
                             g_0_yzzzzz_0_xxxyyy_0,  \
                             g_0_yzzzzz_0_xxxyyy_1,  \
                             g_0_yzzzzz_0_xxxyyz_0,  \
                             g_0_yzzzzz_0_xxxyyz_1,  \
                             g_0_yzzzzz_0_xxxyz_1,   \
                             g_0_yzzzzz_0_xxxyzz_0,  \
                             g_0_yzzzzz_0_xxxyzz_1,  \
                             g_0_yzzzzz_0_xxyyy_1,   \
                             g_0_yzzzzz_0_xxyyyy_0,  \
                             g_0_yzzzzz_0_xxyyyy_1,  \
                             g_0_yzzzzz_0_xxyyyz_0,  \
                             g_0_yzzzzz_0_xxyyyz_1,  \
                             g_0_yzzzzz_0_xxyyz_1,   \
                             g_0_yzzzzz_0_xxyyzz_0,  \
                             g_0_yzzzzz_0_xxyyzz_1,  \
                             g_0_yzzzzz_0_xxyzz_1,   \
                             g_0_yzzzzz_0_xxyzzz_0,  \
                             g_0_yzzzzz_0_xxyzzz_1,  \
                             g_0_yzzzzz_0_xyyyy_1,   \
                             g_0_yzzzzz_0_xyyyyy_0,  \
                             g_0_yzzzzz_0_xyyyyy_1,  \
                             g_0_yzzzzz_0_xyyyyz_0,  \
                             g_0_yzzzzz_0_xyyyyz_1,  \
                             g_0_yzzzzz_0_xyyyz_1,   \
                             g_0_yzzzzz_0_xyyyzz_0,  \
                             g_0_yzzzzz_0_xyyyzz_1,  \
                             g_0_yzzzzz_0_xyyzz_1,   \
                             g_0_yzzzzz_0_xyyzzz_0,  \
                             g_0_yzzzzz_0_xyyzzz_1,  \
                             g_0_yzzzzz_0_xyzzz_1,   \
                             g_0_yzzzzz_0_xyzzzz_0,  \
                             g_0_yzzzzz_0_xyzzzz_1,  \
                             g_0_yzzzzz_0_yyyyy_1,   \
                             g_0_yzzzzz_0_yyyyyy_0,  \
                             g_0_yzzzzz_0_yyyyyy_1,  \
                             g_0_yzzzzz_0_yyyyyz_0,  \
                             g_0_yzzzzz_0_yyyyyz_1,  \
                             g_0_yzzzzz_0_yyyyz_1,   \
                             g_0_yzzzzz_0_yyyyzz_0,  \
                             g_0_yzzzzz_0_yyyyzz_1,  \
                             g_0_yzzzzz_0_yyyzz_1,   \
                             g_0_yzzzzz_0_yyyzzz_0,  \
                             g_0_yzzzzz_0_yyyzzz_1,  \
                             g_0_yzzzzz_0_yyzzz_1,   \
                             g_0_yzzzzz_0_yyzzzz_0,  \
                             g_0_yzzzzz_0_yyzzzz_1,  \
                             g_0_yzzzzz_0_yzzzz_1,   \
                             g_0_yzzzzz_0_yzzzzz_0,  \
                             g_0_yzzzzz_0_yzzzzz_1,  \
                             g_0_yzzzzz_0_zzzzzz_0,  \
                             g_0_yzzzzz_0_zzzzzz_1,  \
                             wp_x,                   \
                             wp_y,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzzz_0_xxxxxx_0[i] = g_0_xzzzzz_0_xxxxxx_0[i] * pb_y + g_0_xzzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxxxxy_0[i] = 5.0 * g_0_yzzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxxy_0[i] * pb_x + g_0_yzzzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxxxz_0[i] = g_0_xzzzzz_0_xxxxxz_0[i] * pb_y + g_0_xzzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxxxyy_0[i] = 4.0 * g_0_yzzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxyy_0[i] * pb_x + g_0_yzzzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxxyz_0[i] = 4.0 * g_0_yzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxxyz_0[i] * pb_x + g_0_yzzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxxzz_0[i] = g_0_xzzzzz_0_xxxxzz_0[i] * pb_y + g_0_xzzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxxyyy_0[i] = 3.0 * g_0_yzzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyyy_0[i] * pb_x + g_0_yzzzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxyyz_0[i] = 3.0 * g_0_yzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyyz_0[i] * pb_x + g_0_yzzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxyzz_0[i] = 3.0 * g_0_yzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyzz_0[i] * pb_x + g_0_yzzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxxzzz_0[i] = g_0_xzzzzz_0_xxxzzz_0[i] * pb_y + g_0_xzzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xxyyyy_0[i] = 2.0 * g_0_yzzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyyy_0[i] * pb_x + g_0_yzzzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxyyyz_0[i] = 2.0 * g_0_yzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyyz_0[i] * pb_x + g_0_yzzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxyyzz_0[i] = 2.0 * g_0_yzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyzz_0[i] * pb_x + g_0_yzzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxyzzz_0[i] = 2.0 * g_0_yzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyzzz_0[i] * pb_x + g_0_yzzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xxzzzz_0[i] = g_0_xzzzzz_0_xxzzzz_0[i] * pb_y + g_0_xzzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_xyyyyy_0[i] = g_0_yzzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyyy_0[i] * pb_x + g_0_yzzzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyyyyz_0[i] = g_0_yzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyyz_0[i] * pb_x + g_0_yzzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyyyzz_0[i] = g_0_yzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyzz_0[i] * pb_x + g_0_yzzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyyzzz_0[i] = g_0_yzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyzzz_0[i] * pb_x + g_0_yzzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xyzzzz_0[i] = g_0_yzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyzzzz_0[i] * pb_x + g_0_yzzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_xzzzzz_0[i] = g_0_xzzzzz_0_xzzzzz_0[i] * pb_y + g_0_xzzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xyzzzzz_0_yyyyyy_0[i] = g_0_yzzzzz_0_yyyyyy_0[i] * pb_x + g_0_yzzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyyyyz_0[i] = g_0_yzzzzz_0_yyyyyz_0[i] * pb_x + g_0_yzzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyyyzz_0[i] = g_0_yzzzzz_0_yyyyzz_0[i] * pb_x + g_0_yzzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyyzzz_0[i] = g_0_yzzzzz_0_yyyzzz_0[i] * pb_x + g_0_yzzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yyzzzz_0[i] = g_0_yzzzzz_0_yyzzzz_0[i] * pb_x + g_0_yzzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_yzzzzz_0[i] = g_0_yzzzzz_0_yzzzzz_0[i] * pb_x + g_0_yzzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyzzzzz_0_zzzzzz_0[i] = g_0_yzzzzz_0_zzzzzz_0[i] * pb_x + g_0_yzzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 756-784 components of targeted buffer : SKSI

    auto g_0_xzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 756);

    auto g_0_xzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 757);

    auto g_0_xzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 758);

    auto g_0_xzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 759);

    auto g_0_xzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 760);

    auto g_0_xzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 761);

    auto g_0_xzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 762);

    auto g_0_xzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 763);

    auto g_0_xzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 764);

    auto g_0_xzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 765);

    auto g_0_xzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 766);

    auto g_0_xzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 767);

    auto g_0_xzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 768);

    auto g_0_xzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 769);

    auto g_0_xzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 770);

    auto g_0_xzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 771);

    auto g_0_xzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 772);

    auto g_0_xzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 773);

    auto g_0_xzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 774);

    auto g_0_xzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 775);

    auto g_0_xzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 776);

    auto g_0_xzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 777);

    auto g_0_xzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 778);

    auto g_0_xzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 779);

    auto g_0_xzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 780);

    auto g_0_xzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 781);

    auto g_0_xzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 782);

    auto g_0_xzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 783);

#pragma omp simd aligned(g_0_xzzzzzz_0_xxxxxx_0,     \
                             g_0_xzzzzzz_0_xxxxxy_0, \
                             g_0_xzzzzzz_0_xxxxxz_0, \
                             g_0_xzzzzzz_0_xxxxyy_0, \
                             g_0_xzzzzzz_0_xxxxyz_0, \
                             g_0_xzzzzzz_0_xxxxzz_0, \
                             g_0_xzzzzzz_0_xxxyyy_0, \
                             g_0_xzzzzzz_0_xxxyyz_0, \
                             g_0_xzzzzzz_0_xxxyzz_0, \
                             g_0_xzzzzzz_0_xxxzzz_0, \
                             g_0_xzzzzzz_0_xxyyyy_0, \
                             g_0_xzzzzzz_0_xxyyyz_0, \
                             g_0_xzzzzzz_0_xxyyzz_0, \
                             g_0_xzzzzzz_0_xxyzzz_0, \
                             g_0_xzzzzzz_0_xxzzzz_0, \
                             g_0_xzzzzzz_0_xyyyyy_0, \
                             g_0_xzzzzzz_0_xyyyyz_0, \
                             g_0_xzzzzzz_0_xyyyzz_0, \
                             g_0_xzzzzzz_0_xyyzzz_0, \
                             g_0_xzzzzzz_0_xyzzzz_0, \
                             g_0_xzzzzzz_0_xzzzzz_0, \
                             g_0_xzzzzzz_0_yyyyyy_0, \
                             g_0_xzzzzzz_0_yyyyyz_0, \
                             g_0_xzzzzzz_0_yyyyzz_0, \
                             g_0_xzzzzzz_0_yyyzzz_0, \
                             g_0_xzzzzzz_0_yyzzzz_0, \
                             g_0_xzzzzzz_0_yzzzzz_0, \
                             g_0_xzzzzzz_0_zzzzzz_0, \
                             g_0_zzzzzz_0_xxxxx_1,   \
                             g_0_zzzzzz_0_xxxxxx_0,  \
                             g_0_zzzzzz_0_xxxxxx_1,  \
                             g_0_zzzzzz_0_xxxxxy_0,  \
                             g_0_zzzzzz_0_xxxxxy_1,  \
                             g_0_zzzzzz_0_xxxxxz_0,  \
                             g_0_zzzzzz_0_xxxxxz_1,  \
                             g_0_zzzzzz_0_xxxxy_1,   \
                             g_0_zzzzzz_0_xxxxyy_0,  \
                             g_0_zzzzzz_0_xxxxyy_1,  \
                             g_0_zzzzzz_0_xxxxyz_0,  \
                             g_0_zzzzzz_0_xxxxyz_1,  \
                             g_0_zzzzzz_0_xxxxz_1,   \
                             g_0_zzzzzz_0_xxxxzz_0,  \
                             g_0_zzzzzz_0_xxxxzz_1,  \
                             g_0_zzzzzz_0_xxxyy_1,   \
                             g_0_zzzzzz_0_xxxyyy_0,  \
                             g_0_zzzzzz_0_xxxyyy_1,  \
                             g_0_zzzzzz_0_xxxyyz_0,  \
                             g_0_zzzzzz_0_xxxyyz_1,  \
                             g_0_zzzzzz_0_xxxyz_1,   \
                             g_0_zzzzzz_0_xxxyzz_0,  \
                             g_0_zzzzzz_0_xxxyzz_1,  \
                             g_0_zzzzzz_0_xxxzz_1,   \
                             g_0_zzzzzz_0_xxxzzz_0,  \
                             g_0_zzzzzz_0_xxxzzz_1,  \
                             g_0_zzzzzz_0_xxyyy_1,   \
                             g_0_zzzzzz_0_xxyyyy_0,  \
                             g_0_zzzzzz_0_xxyyyy_1,  \
                             g_0_zzzzzz_0_xxyyyz_0,  \
                             g_0_zzzzzz_0_xxyyyz_1,  \
                             g_0_zzzzzz_0_xxyyz_1,   \
                             g_0_zzzzzz_0_xxyyzz_0,  \
                             g_0_zzzzzz_0_xxyyzz_1,  \
                             g_0_zzzzzz_0_xxyzz_1,   \
                             g_0_zzzzzz_0_xxyzzz_0,  \
                             g_0_zzzzzz_0_xxyzzz_1,  \
                             g_0_zzzzzz_0_xxzzz_1,   \
                             g_0_zzzzzz_0_xxzzzz_0,  \
                             g_0_zzzzzz_0_xxzzzz_1,  \
                             g_0_zzzzzz_0_xyyyy_1,   \
                             g_0_zzzzzz_0_xyyyyy_0,  \
                             g_0_zzzzzz_0_xyyyyy_1,  \
                             g_0_zzzzzz_0_xyyyyz_0,  \
                             g_0_zzzzzz_0_xyyyyz_1,  \
                             g_0_zzzzzz_0_xyyyz_1,   \
                             g_0_zzzzzz_0_xyyyzz_0,  \
                             g_0_zzzzzz_0_xyyyzz_1,  \
                             g_0_zzzzzz_0_xyyzz_1,   \
                             g_0_zzzzzz_0_xyyzzz_0,  \
                             g_0_zzzzzz_0_xyyzzz_1,  \
                             g_0_zzzzzz_0_xyzzz_1,   \
                             g_0_zzzzzz_0_xyzzzz_0,  \
                             g_0_zzzzzz_0_xyzzzz_1,  \
                             g_0_zzzzzz_0_xzzzz_1,   \
                             g_0_zzzzzz_0_xzzzzz_0,  \
                             g_0_zzzzzz_0_xzzzzz_1,  \
                             g_0_zzzzzz_0_yyyyy_1,   \
                             g_0_zzzzzz_0_yyyyyy_0,  \
                             g_0_zzzzzz_0_yyyyyy_1,  \
                             g_0_zzzzzz_0_yyyyyz_0,  \
                             g_0_zzzzzz_0_yyyyyz_1,  \
                             g_0_zzzzzz_0_yyyyz_1,   \
                             g_0_zzzzzz_0_yyyyzz_0,  \
                             g_0_zzzzzz_0_yyyyzz_1,  \
                             g_0_zzzzzz_0_yyyzz_1,   \
                             g_0_zzzzzz_0_yyyzzz_0,  \
                             g_0_zzzzzz_0_yyyzzz_1,  \
                             g_0_zzzzzz_0_yyzzz_1,   \
                             g_0_zzzzzz_0_yyzzzz_0,  \
                             g_0_zzzzzz_0_yyzzzz_1,  \
                             g_0_zzzzzz_0_yzzzz_1,   \
                             g_0_zzzzzz_0_yzzzzz_0,  \
                             g_0_zzzzzz_0_yzzzzz_1,  \
                             g_0_zzzzzz_0_zzzzz_1,   \
                             g_0_zzzzzz_0_zzzzzz_0,  \
                             g_0_zzzzzz_0_zzzzzz_1,  \
                             wp_x,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzz_0_xxxxxx_0[i] = 6.0 * g_0_zzzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxx_0[i] * pb_x + g_0_zzzzzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxxy_0[i] = 5.0 * g_0_zzzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxy_0[i] * pb_x + g_0_zzzzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxxz_0[i] = 5.0 * g_0_zzzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxz_0[i] * pb_x + g_0_zzzzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxyy_0[i] = 4.0 * g_0_zzzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyy_0[i] * pb_x + g_0_zzzzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxyz_0[i] = 4.0 * g_0_zzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyz_0[i] * pb_x + g_0_zzzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxxzz_0[i] = 4.0 * g_0_zzzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxzz_0[i] * pb_x + g_0_zzzzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxyyy_0[i] = 3.0 * g_0_zzzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyy_0[i] * pb_x + g_0_zzzzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxyyz_0[i] = 3.0 * g_0_zzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyz_0[i] * pb_x + g_0_zzzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxyzz_0[i] = 3.0 * g_0_zzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyzz_0[i] * pb_x + g_0_zzzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxxzzz_0[i] = 3.0 * g_0_zzzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxzzz_0[i] * pb_x + g_0_zzzzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyyyy_0[i] = 2.0 * g_0_zzzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyy_0[i] * pb_x + g_0_zzzzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyyyz_0[i] = 2.0 * g_0_zzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyz_0[i] * pb_x + g_0_zzzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyyzz_0[i] = 2.0 * g_0_zzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyzz_0[i] * pb_x + g_0_zzzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxyzzz_0[i] = 2.0 * g_0_zzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyzzz_0[i] * pb_x + g_0_zzzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xxzzzz_0[i] = 2.0 * g_0_zzzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxzzzz_0[i] * pb_x + g_0_zzzzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyyyy_0[i] = g_0_zzzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyy_0[i] * pb_x + g_0_zzzzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyyyz_0[i] = g_0_zzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyz_0[i] * pb_x + g_0_zzzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyyzz_0[i] = g_0_zzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyzz_0[i] * pb_x + g_0_zzzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyyzzz_0[i] = g_0_zzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyzzz_0[i] * pb_x + g_0_zzzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xyzzzz_0[i] = g_0_zzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyzzzz_0[i] * pb_x + g_0_zzzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_xzzzzz_0[i] = g_0_zzzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xzzzzz_0[i] * pb_x + g_0_zzzzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyyyy_0[i] = g_0_zzzzzz_0_yyyyyy_0[i] * pb_x + g_0_zzzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyyyz_0[i] = g_0_zzzzzz_0_yyyyyz_0[i] * pb_x + g_0_zzzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyyzz_0[i] = g_0_zzzzzz_0_yyyyzz_0[i] * pb_x + g_0_zzzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyyzzz_0[i] = g_0_zzzzzz_0_yyyzzz_0[i] * pb_x + g_0_zzzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yyzzzz_0[i] = g_0_zzzzzz_0_yyzzzz_0[i] * pb_x + g_0_zzzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_yzzzzz_0[i] = g_0_zzzzzz_0_yzzzzz_0[i] * pb_x + g_0_zzzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzz_0_zzzzzz_0[i] = g_0_zzzzzz_0_zzzzzz_0[i] * pb_x + g_0_zzzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 784-812 components of targeted buffer : SKSI

    auto g_0_yyyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 784);

    auto g_0_yyyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 785);

    auto g_0_yyyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 786);

    auto g_0_yyyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 787);

    auto g_0_yyyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 788);

    auto g_0_yyyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 789);

    auto g_0_yyyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 790);

    auto g_0_yyyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 791);

    auto g_0_yyyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 792);

    auto g_0_yyyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 793);

    auto g_0_yyyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 794);

    auto g_0_yyyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 795);

    auto g_0_yyyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 796);

    auto g_0_yyyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 797);

    auto g_0_yyyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 798);

    auto g_0_yyyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 799);

    auto g_0_yyyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 800);

    auto g_0_yyyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 801);

    auto g_0_yyyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 802);

    auto g_0_yyyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 803);

    auto g_0_yyyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 804);

    auto g_0_yyyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 805);

    auto g_0_yyyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 806);

    auto g_0_yyyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 807);

    auto g_0_yyyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 808);

    auto g_0_yyyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 809);

    auto g_0_yyyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 810);

    auto g_0_yyyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 811);

#pragma omp simd aligned(g_0_yyyyy_0_xxxxxx_0,       \
                             g_0_yyyyy_0_xxxxxx_1,   \
                             g_0_yyyyy_0_xxxxxy_0,   \
                             g_0_yyyyy_0_xxxxxy_1,   \
                             g_0_yyyyy_0_xxxxxz_0,   \
                             g_0_yyyyy_0_xxxxxz_1,   \
                             g_0_yyyyy_0_xxxxyy_0,   \
                             g_0_yyyyy_0_xxxxyy_1,   \
                             g_0_yyyyy_0_xxxxyz_0,   \
                             g_0_yyyyy_0_xxxxyz_1,   \
                             g_0_yyyyy_0_xxxxzz_0,   \
                             g_0_yyyyy_0_xxxxzz_1,   \
                             g_0_yyyyy_0_xxxyyy_0,   \
                             g_0_yyyyy_0_xxxyyy_1,   \
                             g_0_yyyyy_0_xxxyyz_0,   \
                             g_0_yyyyy_0_xxxyyz_1,   \
                             g_0_yyyyy_0_xxxyzz_0,   \
                             g_0_yyyyy_0_xxxyzz_1,   \
                             g_0_yyyyy_0_xxxzzz_0,   \
                             g_0_yyyyy_0_xxxzzz_1,   \
                             g_0_yyyyy_0_xxyyyy_0,   \
                             g_0_yyyyy_0_xxyyyy_1,   \
                             g_0_yyyyy_0_xxyyyz_0,   \
                             g_0_yyyyy_0_xxyyyz_1,   \
                             g_0_yyyyy_0_xxyyzz_0,   \
                             g_0_yyyyy_0_xxyyzz_1,   \
                             g_0_yyyyy_0_xxyzzz_0,   \
                             g_0_yyyyy_0_xxyzzz_1,   \
                             g_0_yyyyy_0_xxzzzz_0,   \
                             g_0_yyyyy_0_xxzzzz_1,   \
                             g_0_yyyyy_0_xyyyyy_0,   \
                             g_0_yyyyy_0_xyyyyy_1,   \
                             g_0_yyyyy_0_xyyyyz_0,   \
                             g_0_yyyyy_0_xyyyyz_1,   \
                             g_0_yyyyy_0_xyyyzz_0,   \
                             g_0_yyyyy_0_xyyyzz_1,   \
                             g_0_yyyyy_0_xyyzzz_0,   \
                             g_0_yyyyy_0_xyyzzz_1,   \
                             g_0_yyyyy_0_xyzzzz_0,   \
                             g_0_yyyyy_0_xyzzzz_1,   \
                             g_0_yyyyy_0_xzzzzz_0,   \
                             g_0_yyyyy_0_xzzzzz_1,   \
                             g_0_yyyyy_0_yyyyyy_0,   \
                             g_0_yyyyy_0_yyyyyy_1,   \
                             g_0_yyyyy_0_yyyyyz_0,   \
                             g_0_yyyyy_0_yyyyyz_1,   \
                             g_0_yyyyy_0_yyyyzz_0,   \
                             g_0_yyyyy_0_yyyyzz_1,   \
                             g_0_yyyyy_0_yyyzzz_0,   \
                             g_0_yyyyy_0_yyyzzz_1,   \
                             g_0_yyyyy_0_yyzzzz_0,   \
                             g_0_yyyyy_0_yyzzzz_1,   \
                             g_0_yyyyy_0_yzzzzz_0,   \
                             g_0_yyyyy_0_yzzzzz_1,   \
                             g_0_yyyyy_0_zzzzzz_0,   \
                             g_0_yyyyy_0_zzzzzz_1,   \
                             g_0_yyyyyy_0_xxxxx_1,   \
                             g_0_yyyyyy_0_xxxxxx_0,  \
                             g_0_yyyyyy_0_xxxxxx_1,  \
                             g_0_yyyyyy_0_xxxxxy_0,  \
                             g_0_yyyyyy_0_xxxxxy_1,  \
                             g_0_yyyyyy_0_xxxxxz_0,  \
                             g_0_yyyyyy_0_xxxxxz_1,  \
                             g_0_yyyyyy_0_xxxxy_1,   \
                             g_0_yyyyyy_0_xxxxyy_0,  \
                             g_0_yyyyyy_0_xxxxyy_1,  \
                             g_0_yyyyyy_0_xxxxyz_0,  \
                             g_0_yyyyyy_0_xxxxyz_1,  \
                             g_0_yyyyyy_0_xxxxz_1,   \
                             g_0_yyyyyy_0_xxxxzz_0,  \
                             g_0_yyyyyy_0_xxxxzz_1,  \
                             g_0_yyyyyy_0_xxxyy_1,   \
                             g_0_yyyyyy_0_xxxyyy_0,  \
                             g_0_yyyyyy_0_xxxyyy_1,  \
                             g_0_yyyyyy_0_xxxyyz_0,  \
                             g_0_yyyyyy_0_xxxyyz_1,  \
                             g_0_yyyyyy_0_xxxyz_1,   \
                             g_0_yyyyyy_0_xxxyzz_0,  \
                             g_0_yyyyyy_0_xxxyzz_1,  \
                             g_0_yyyyyy_0_xxxzz_1,   \
                             g_0_yyyyyy_0_xxxzzz_0,  \
                             g_0_yyyyyy_0_xxxzzz_1,  \
                             g_0_yyyyyy_0_xxyyy_1,   \
                             g_0_yyyyyy_0_xxyyyy_0,  \
                             g_0_yyyyyy_0_xxyyyy_1,  \
                             g_0_yyyyyy_0_xxyyyz_0,  \
                             g_0_yyyyyy_0_xxyyyz_1,  \
                             g_0_yyyyyy_0_xxyyz_1,   \
                             g_0_yyyyyy_0_xxyyzz_0,  \
                             g_0_yyyyyy_0_xxyyzz_1,  \
                             g_0_yyyyyy_0_xxyzz_1,   \
                             g_0_yyyyyy_0_xxyzzz_0,  \
                             g_0_yyyyyy_0_xxyzzz_1,  \
                             g_0_yyyyyy_0_xxzzz_1,   \
                             g_0_yyyyyy_0_xxzzzz_0,  \
                             g_0_yyyyyy_0_xxzzzz_1,  \
                             g_0_yyyyyy_0_xyyyy_1,   \
                             g_0_yyyyyy_0_xyyyyy_0,  \
                             g_0_yyyyyy_0_xyyyyy_1,  \
                             g_0_yyyyyy_0_xyyyyz_0,  \
                             g_0_yyyyyy_0_xyyyyz_1,  \
                             g_0_yyyyyy_0_xyyyz_1,   \
                             g_0_yyyyyy_0_xyyyzz_0,  \
                             g_0_yyyyyy_0_xyyyzz_1,  \
                             g_0_yyyyyy_0_xyyzz_1,   \
                             g_0_yyyyyy_0_xyyzzz_0,  \
                             g_0_yyyyyy_0_xyyzzz_1,  \
                             g_0_yyyyyy_0_xyzzz_1,   \
                             g_0_yyyyyy_0_xyzzzz_0,  \
                             g_0_yyyyyy_0_xyzzzz_1,  \
                             g_0_yyyyyy_0_xzzzz_1,   \
                             g_0_yyyyyy_0_xzzzzz_0,  \
                             g_0_yyyyyy_0_xzzzzz_1,  \
                             g_0_yyyyyy_0_yyyyy_1,   \
                             g_0_yyyyyy_0_yyyyyy_0,  \
                             g_0_yyyyyy_0_yyyyyy_1,  \
                             g_0_yyyyyy_0_yyyyyz_0,  \
                             g_0_yyyyyy_0_yyyyyz_1,  \
                             g_0_yyyyyy_0_yyyyz_1,   \
                             g_0_yyyyyy_0_yyyyzz_0,  \
                             g_0_yyyyyy_0_yyyyzz_1,  \
                             g_0_yyyyyy_0_yyyzz_1,   \
                             g_0_yyyyyy_0_yyyzzz_0,  \
                             g_0_yyyyyy_0_yyyzzz_1,  \
                             g_0_yyyyyy_0_yyzzz_1,   \
                             g_0_yyyyyy_0_yyzzzz_0,  \
                             g_0_yyyyyy_0_yyzzzz_1,  \
                             g_0_yyyyyy_0_yzzzz_1,   \
                             g_0_yyyyyy_0_yzzzzz_0,  \
                             g_0_yyyyyy_0_yzzzzz_1,  \
                             g_0_yyyyyy_0_zzzzz_1,   \
                             g_0_yyyyyy_0_zzzzzz_0,  \
                             g_0_yyyyyy_0_zzzzzz_1,  \
                             g_0_yyyyyyy_0_xxxxxx_0, \
                             g_0_yyyyyyy_0_xxxxxy_0, \
                             g_0_yyyyyyy_0_xxxxxz_0, \
                             g_0_yyyyyyy_0_xxxxyy_0, \
                             g_0_yyyyyyy_0_xxxxyz_0, \
                             g_0_yyyyyyy_0_xxxxzz_0, \
                             g_0_yyyyyyy_0_xxxyyy_0, \
                             g_0_yyyyyyy_0_xxxyyz_0, \
                             g_0_yyyyyyy_0_xxxyzz_0, \
                             g_0_yyyyyyy_0_xxxzzz_0, \
                             g_0_yyyyyyy_0_xxyyyy_0, \
                             g_0_yyyyyyy_0_xxyyyz_0, \
                             g_0_yyyyyyy_0_xxyyzz_0, \
                             g_0_yyyyyyy_0_xxyzzz_0, \
                             g_0_yyyyyyy_0_xxzzzz_0, \
                             g_0_yyyyyyy_0_xyyyyy_0, \
                             g_0_yyyyyyy_0_xyyyyz_0, \
                             g_0_yyyyyyy_0_xyyyzz_0, \
                             g_0_yyyyyyy_0_xyyzzz_0, \
                             g_0_yyyyyyy_0_xyzzzz_0, \
                             g_0_yyyyyyy_0_xzzzzz_0, \
                             g_0_yyyyyyy_0_yyyyyy_0, \
                             g_0_yyyyyyy_0_yyyyyz_0, \
                             g_0_yyyyyyy_0_yyyyzz_0, \
                             g_0_yyyyyyy_0_yyyzzz_0, \
                             g_0_yyyyyyy_0_yyzzzz_0, \
                             g_0_yyyyyyy_0_yzzzzz_0, \
                             g_0_yyyyyyy_0_zzzzzz_0, \
                             wp_y,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyy_0_xxxxxx_0[i] = 6.0 * g_0_yyyyy_0_xxxxxx_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxxx_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_xxxxxx_0[i] * pb_y + g_0_yyyyyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxxy_0[i] = 6.0 * g_0_yyyyy_0_xxxxxy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxxy_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxy_0[i] * pb_y + g_0_yyyyyy_0_xxxxxy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxxz_0[i] = 6.0 * g_0_yyyyy_0_xxxxxz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxxz_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_xxxxxz_0[i] * pb_y + g_0_yyyyyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxyy_0[i] = 6.0 * g_0_yyyyy_0_xxxxyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxyy_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyy_0[i] * pb_y + g_0_yyyyyy_0_xxxxyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxyz_0[i] = 6.0 * g_0_yyyyy_0_xxxxyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxyz_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyz_0[i] * pb_y + g_0_yyyyyy_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxxzz_0[i] = 6.0 * g_0_yyyyy_0_xxxxzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxxzz_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_xxxxzz_0[i] * pb_y + g_0_yyyyyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxyyy_0[i] = 6.0 * g_0_yyyyy_0_xxxyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxyyy_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyy_0[i] * pb_y + g_0_yyyyyy_0_xxxyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxyyz_0[i] = 6.0 * g_0_yyyyy_0_xxxyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyz_0[i] * pb_y + g_0_yyyyyy_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxyzz_0[i] = 6.0 * g_0_yyyyy_0_xxxyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxyzz_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyzz_0[i] * pb_y + g_0_yyyyyy_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxxzzz_0[i] = 6.0 * g_0_yyyyy_0_xxxzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxxzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_xxxzzz_0[i] * pb_y + g_0_yyyyyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyyyy_0[i] = 6.0 * g_0_yyyyy_0_xxyyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyyyy_1[i] * fti_ab_0 +
                                    4.0 * g_0_yyyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyy_0[i] * pb_y + g_0_yyyyyy_0_xxyyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyyyz_0[i] = 6.0 * g_0_yyyyy_0_xxyyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyz_0[i] * pb_y + g_0_yyyyyy_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyyzz_0[i] = 6.0 * g_0_yyyyy_0_xxyyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyzz_0[i] * pb_y + g_0_yyyyyy_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxyzzz_0[i] = 6.0 * g_0_yyyyy_0_xxyzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxyzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyzzz_0[i] * pb_y + g_0_yyyyyy_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xxzzzz_0[i] = 6.0 * g_0_yyyyy_0_xxzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xxzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_xxzzzz_0[i] * pb_y + g_0_yyyyyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyyyy_0[i] = 6.0 * g_0_yyyyy_0_xyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyyyy_1[i] * fti_ab_0 +
                                    5.0 * g_0_yyyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyy_0[i] * pb_y + g_0_yyyyyy_0_xyyyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyyyz_0[i] = 6.0 * g_0_yyyyy_0_xyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyyyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_yyyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyz_0[i] * pb_y + g_0_yyyyyy_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyyzz_0[i] = 6.0 * g_0_yyyyy_0_xyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyzz_0[i] * pb_y + g_0_yyyyyy_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyyzzz_0[i] = 6.0 * g_0_yyyyy_0_xyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyzzz_0[i] * pb_y + g_0_yyyyyy_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xyzzzz_0[i] = 6.0 * g_0_yyyyy_0_xyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xyzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyzzzz_0[i] * pb_y + g_0_yyyyyy_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_xzzzzz_0[i] = 6.0 * g_0_yyyyy_0_xzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_xzzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_xzzzzz_0[i] * pb_y + g_0_yyyyyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyyyy_0[i] = 6.0 * g_0_yyyyy_0_yyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyyyy_1[i] * fti_ab_0 +
                                    6.0 * g_0_yyyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyyy_0[i] * pb_y + g_0_yyyyyy_0_yyyyyy_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyyyz_0[i] = 6.0 * g_0_yyyyy_0_yyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyyyz_1[i] * fti_ab_0 +
                                    5.0 * g_0_yyyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyyz_0[i] * pb_y + g_0_yyyyyy_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyyzz_0[i] = 6.0 * g_0_yyyyy_0_yyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyyzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_yyyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyzz_0[i] * pb_y + g_0_yyyyyy_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyyzzz_0[i] = 6.0 * g_0_yyyyy_0_yyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyyzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyzzz_0[i] * pb_y + g_0_yyyyyy_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yyzzzz_0[i] = 6.0 * g_0_yyyyy_0_yyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yyzzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyzzzz_0[i] * pb_y + g_0_yyyyyy_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_yzzzzz_0[i] = 6.0 * g_0_yyyyy_0_yzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_yzzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yzzzzz_0[i] * pb_y + g_0_yyyyyy_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyy_0_zzzzzz_0[i] = 6.0 * g_0_yyyyy_0_zzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_zzzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyyy_0_zzzzzz_0[i] * pb_y + g_0_yyyyyy_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 812-840 components of targeted buffer : SKSI

    auto g_0_yyyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 812);

    auto g_0_yyyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 813);

    auto g_0_yyyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 814);

    auto g_0_yyyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 815);

    auto g_0_yyyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 816);

    auto g_0_yyyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 817);

    auto g_0_yyyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 818);

    auto g_0_yyyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 819);

    auto g_0_yyyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 820);

    auto g_0_yyyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 821);

    auto g_0_yyyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 822);

    auto g_0_yyyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 823);

    auto g_0_yyyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 824);

    auto g_0_yyyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 825);

    auto g_0_yyyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 826);

    auto g_0_yyyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 827);

    auto g_0_yyyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 828);

    auto g_0_yyyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 829);

    auto g_0_yyyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 830);

    auto g_0_yyyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 831);

    auto g_0_yyyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 832);

    auto g_0_yyyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 833);

    auto g_0_yyyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 834);

    auto g_0_yyyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 835);

    auto g_0_yyyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 836);

    auto g_0_yyyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 837);

    auto g_0_yyyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 838);

    auto g_0_yyyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 839);

#pragma omp simd aligned(g_0_yyyyyy_0_xxxxx_1,       \
                             g_0_yyyyyy_0_xxxxxx_0,  \
                             g_0_yyyyyy_0_xxxxxx_1,  \
                             g_0_yyyyyy_0_xxxxxy_0,  \
                             g_0_yyyyyy_0_xxxxxy_1,  \
                             g_0_yyyyyy_0_xxxxxz_0,  \
                             g_0_yyyyyy_0_xxxxxz_1,  \
                             g_0_yyyyyy_0_xxxxy_1,   \
                             g_0_yyyyyy_0_xxxxyy_0,  \
                             g_0_yyyyyy_0_xxxxyy_1,  \
                             g_0_yyyyyy_0_xxxxyz_0,  \
                             g_0_yyyyyy_0_xxxxyz_1,  \
                             g_0_yyyyyy_0_xxxxz_1,   \
                             g_0_yyyyyy_0_xxxxzz_0,  \
                             g_0_yyyyyy_0_xxxxzz_1,  \
                             g_0_yyyyyy_0_xxxyy_1,   \
                             g_0_yyyyyy_0_xxxyyy_0,  \
                             g_0_yyyyyy_0_xxxyyy_1,  \
                             g_0_yyyyyy_0_xxxyyz_0,  \
                             g_0_yyyyyy_0_xxxyyz_1,  \
                             g_0_yyyyyy_0_xxxyz_1,   \
                             g_0_yyyyyy_0_xxxyzz_0,  \
                             g_0_yyyyyy_0_xxxyzz_1,  \
                             g_0_yyyyyy_0_xxxzz_1,   \
                             g_0_yyyyyy_0_xxxzzz_0,  \
                             g_0_yyyyyy_0_xxxzzz_1,  \
                             g_0_yyyyyy_0_xxyyy_1,   \
                             g_0_yyyyyy_0_xxyyyy_0,  \
                             g_0_yyyyyy_0_xxyyyy_1,  \
                             g_0_yyyyyy_0_xxyyyz_0,  \
                             g_0_yyyyyy_0_xxyyyz_1,  \
                             g_0_yyyyyy_0_xxyyz_1,   \
                             g_0_yyyyyy_0_xxyyzz_0,  \
                             g_0_yyyyyy_0_xxyyzz_1,  \
                             g_0_yyyyyy_0_xxyzz_1,   \
                             g_0_yyyyyy_0_xxyzzz_0,  \
                             g_0_yyyyyy_0_xxyzzz_1,  \
                             g_0_yyyyyy_0_xxzzz_1,   \
                             g_0_yyyyyy_0_xxzzzz_0,  \
                             g_0_yyyyyy_0_xxzzzz_1,  \
                             g_0_yyyyyy_0_xyyyy_1,   \
                             g_0_yyyyyy_0_xyyyyy_0,  \
                             g_0_yyyyyy_0_xyyyyy_1,  \
                             g_0_yyyyyy_0_xyyyyz_0,  \
                             g_0_yyyyyy_0_xyyyyz_1,  \
                             g_0_yyyyyy_0_xyyyz_1,   \
                             g_0_yyyyyy_0_xyyyzz_0,  \
                             g_0_yyyyyy_0_xyyyzz_1,  \
                             g_0_yyyyyy_0_xyyzz_1,   \
                             g_0_yyyyyy_0_xyyzzz_0,  \
                             g_0_yyyyyy_0_xyyzzz_1,  \
                             g_0_yyyyyy_0_xyzzz_1,   \
                             g_0_yyyyyy_0_xyzzzz_0,  \
                             g_0_yyyyyy_0_xyzzzz_1,  \
                             g_0_yyyyyy_0_xzzzz_1,   \
                             g_0_yyyyyy_0_xzzzzz_0,  \
                             g_0_yyyyyy_0_xzzzzz_1,  \
                             g_0_yyyyyy_0_yyyyy_1,   \
                             g_0_yyyyyy_0_yyyyyy_0,  \
                             g_0_yyyyyy_0_yyyyyy_1,  \
                             g_0_yyyyyy_0_yyyyyz_0,  \
                             g_0_yyyyyy_0_yyyyyz_1,  \
                             g_0_yyyyyy_0_yyyyz_1,   \
                             g_0_yyyyyy_0_yyyyzz_0,  \
                             g_0_yyyyyy_0_yyyyzz_1,  \
                             g_0_yyyyyy_0_yyyzz_1,   \
                             g_0_yyyyyy_0_yyyzzz_0,  \
                             g_0_yyyyyy_0_yyyzzz_1,  \
                             g_0_yyyyyy_0_yyzzz_1,   \
                             g_0_yyyyyy_0_yyzzzz_0,  \
                             g_0_yyyyyy_0_yyzzzz_1,  \
                             g_0_yyyyyy_0_yzzzz_1,   \
                             g_0_yyyyyy_0_yzzzzz_0,  \
                             g_0_yyyyyy_0_yzzzzz_1,  \
                             g_0_yyyyyy_0_zzzzz_1,   \
                             g_0_yyyyyy_0_zzzzzz_0,  \
                             g_0_yyyyyy_0_zzzzzz_1,  \
                             g_0_yyyyyyz_0_xxxxxx_0, \
                             g_0_yyyyyyz_0_xxxxxy_0, \
                             g_0_yyyyyyz_0_xxxxxz_0, \
                             g_0_yyyyyyz_0_xxxxyy_0, \
                             g_0_yyyyyyz_0_xxxxyz_0, \
                             g_0_yyyyyyz_0_xxxxzz_0, \
                             g_0_yyyyyyz_0_xxxyyy_0, \
                             g_0_yyyyyyz_0_xxxyyz_0, \
                             g_0_yyyyyyz_0_xxxyzz_0, \
                             g_0_yyyyyyz_0_xxxzzz_0, \
                             g_0_yyyyyyz_0_xxyyyy_0, \
                             g_0_yyyyyyz_0_xxyyyz_0, \
                             g_0_yyyyyyz_0_xxyyzz_0, \
                             g_0_yyyyyyz_0_xxyzzz_0, \
                             g_0_yyyyyyz_0_xxzzzz_0, \
                             g_0_yyyyyyz_0_xyyyyy_0, \
                             g_0_yyyyyyz_0_xyyyyz_0, \
                             g_0_yyyyyyz_0_xyyyzz_0, \
                             g_0_yyyyyyz_0_xyyzzz_0, \
                             g_0_yyyyyyz_0_xyzzzz_0, \
                             g_0_yyyyyyz_0_xzzzzz_0, \
                             g_0_yyyyyyz_0_yyyyyy_0, \
                             g_0_yyyyyyz_0_yyyyyz_0, \
                             g_0_yyyyyyz_0_yyyyzz_0, \
                             g_0_yyyyyyz_0_yyyzzz_0, \
                             g_0_yyyyyyz_0_yyzzzz_0, \
                             g_0_yyyyyyz_0_yzzzzz_0, \
                             g_0_yyyyyyz_0_zzzzzz_0, \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyz_0_xxxxxx_0[i] = g_0_yyyyyy_0_xxxxxx_0[i] * pb_z + g_0_yyyyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxxy_0[i] = g_0_yyyyyy_0_xxxxxy_0[i] * pb_z + g_0_yyyyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxxz_0[i] = g_0_yyyyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxxz_0[i] * pb_z + g_0_yyyyyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxyy_0[i] = g_0_yyyyyy_0_xxxxyy_0[i] * pb_z + g_0_yyyyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxyz_0[i] = g_0_yyyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxyz_0[i] * pb_z + g_0_yyyyyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxxzz_0[i] = 2.0 * g_0_yyyyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxxzz_0[i] * pb_z + g_0_yyyyyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxyyy_0[i] = g_0_yyyyyy_0_xxxyyy_0[i] * pb_z + g_0_yyyyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxyyz_0[i] = g_0_yyyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyyz_0[i] * pb_z + g_0_yyyyyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxyzz_0[i] = 2.0 * g_0_yyyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxyzz_0[i] * pb_z + g_0_yyyyyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxxzzz_0[i] = 3.0 * g_0_yyyyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxxzzz_0[i] * pb_z + g_0_yyyyyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyyyy_0[i] = g_0_yyyyyy_0_xxyyyy_0[i] * pb_z + g_0_yyyyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyyyz_0[i] = g_0_yyyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyyz_0[i] * pb_z + g_0_yyyyyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyyzz_0[i] = 2.0 * g_0_yyyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyyzz_0[i] * pb_z + g_0_yyyyyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxyzzz_0[i] = 3.0 * g_0_yyyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxyzzz_0[i] * pb_z + g_0_yyyyyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xxzzzz_0[i] = 4.0 * g_0_yyyyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xxzzzz_0[i] * pb_z + g_0_yyyyyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyyyy_0[i] = g_0_yyyyyy_0_xyyyyy_0[i] * pb_z + g_0_yyyyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyyyz_0[i] = g_0_yyyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyyz_0[i] * pb_z + g_0_yyyyyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyyzz_0[i] = 2.0 * g_0_yyyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyyzz_0[i] * pb_z + g_0_yyyyyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyyzzz_0[i] = 3.0 * g_0_yyyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyyzzz_0[i] * pb_z + g_0_yyyyyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xyzzzz_0[i] = 4.0 * g_0_yyyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xyzzzz_0[i] * pb_z + g_0_yyyyyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_xzzzzz_0[i] = 5.0 * g_0_yyyyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_xzzzzz_0[i] * pb_z + g_0_yyyyyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyyyy_0[i] = g_0_yyyyyy_0_yyyyyy_0[i] * pb_z + g_0_yyyyyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyyyz_0[i] = g_0_yyyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyyz_0[i] * pb_z + g_0_yyyyyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyyzz_0[i] = 2.0 * g_0_yyyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyyzz_0[i] * pb_z + g_0_yyyyyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyyzzz_0[i] = 3.0 * g_0_yyyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyyzzz_0[i] * pb_z + g_0_yyyyyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yyzzzz_0[i] = 4.0 * g_0_yyyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yyzzzz_0[i] * pb_z + g_0_yyyyyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_yzzzzz_0[i] = 5.0 * g_0_yyyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_yzzzzz_0[i] * pb_z + g_0_yyyyyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyz_0_zzzzzz_0[i] = 6.0 * g_0_yyyyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyyy_0_zzzzzz_0[i] * pb_z + g_0_yyyyyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 840-868 components of targeted buffer : SKSI

    auto g_0_yyyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 840);

    auto g_0_yyyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 841);

    auto g_0_yyyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 842);

    auto g_0_yyyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 843);

    auto g_0_yyyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 844);

    auto g_0_yyyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 845);

    auto g_0_yyyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 846);

    auto g_0_yyyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 847);

    auto g_0_yyyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 848);

    auto g_0_yyyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 849);

    auto g_0_yyyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 850);

    auto g_0_yyyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 851);

    auto g_0_yyyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 852);

    auto g_0_yyyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 853);

    auto g_0_yyyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 854);

    auto g_0_yyyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 855);

    auto g_0_yyyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 856);

    auto g_0_yyyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 857);

    auto g_0_yyyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 858);

    auto g_0_yyyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 859);

    auto g_0_yyyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 860);

    auto g_0_yyyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 861);

    auto g_0_yyyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 862);

    auto g_0_yyyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 863);

    auto g_0_yyyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 864);

    auto g_0_yyyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 865);

    auto g_0_yyyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 866);

    auto g_0_yyyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 867);

#pragma omp simd aligned(g_0_yyyyy_0_xxxxxy_0,       \
                             g_0_yyyyy_0_xxxxxy_1,   \
                             g_0_yyyyy_0_xxxxyy_0,   \
                             g_0_yyyyy_0_xxxxyy_1,   \
                             g_0_yyyyy_0_xxxyyy_0,   \
                             g_0_yyyyy_0_xxxyyy_1,   \
                             g_0_yyyyy_0_xxyyyy_0,   \
                             g_0_yyyyy_0_xxyyyy_1,   \
                             g_0_yyyyy_0_xyyyyy_0,   \
                             g_0_yyyyy_0_xyyyyy_1,   \
                             g_0_yyyyy_0_yyyyyy_0,   \
                             g_0_yyyyy_0_yyyyyy_1,   \
                             g_0_yyyyyz_0_xxxxxy_0,  \
                             g_0_yyyyyz_0_xxxxxy_1,  \
                             g_0_yyyyyz_0_xxxxyy_0,  \
                             g_0_yyyyyz_0_xxxxyy_1,  \
                             g_0_yyyyyz_0_xxxyyy_0,  \
                             g_0_yyyyyz_0_xxxyyy_1,  \
                             g_0_yyyyyz_0_xxyyyy_0,  \
                             g_0_yyyyyz_0_xxyyyy_1,  \
                             g_0_yyyyyz_0_xyyyyy_0,  \
                             g_0_yyyyyz_0_xyyyyy_1,  \
                             g_0_yyyyyz_0_yyyyyy_0,  \
                             g_0_yyyyyz_0_yyyyyy_1,  \
                             g_0_yyyyyzz_0_xxxxxx_0, \
                             g_0_yyyyyzz_0_xxxxxy_0, \
                             g_0_yyyyyzz_0_xxxxxz_0, \
                             g_0_yyyyyzz_0_xxxxyy_0, \
                             g_0_yyyyyzz_0_xxxxyz_0, \
                             g_0_yyyyyzz_0_xxxxzz_0, \
                             g_0_yyyyyzz_0_xxxyyy_0, \
                             g_0_yyyyyzz_0_xxxyyz_0, \
                             g_0_yyyyyzz_0_xxxyzz_0, \
                             g_0_yyyyyzz_0_xxxzzz_0, \
                             g_0_yyyyyzz_0_xxyyyy_0, \
                             g_0_yyyyyzz_0_xxyyyz_0, \
                             g_0_yyyyyzz_0_xxyyzz_0, \
                             g_0_yyyyyzz_0_xxyzzz_0, \
                             g_0_yyyyyzz_0_xxzzzz_0, \
                             g_0_yyyyyzz_0_xyyyyy_0, \
                             g_0_yyyyyzz_0_xyyyyz_0, \
                             g_0_yyyyyzz_0_xyyyzz_0, \
                             g_0_yyyyyzz_0_xyyzzz_0, \
                             g_0_yyyyyzz_0_xyzzzz_0, \
                             g_0_yyyyyzz_0_xzzzzz_0, \
                             g_0_yyyyyzz_0_yyyyyy_0, \
                             g_0_yyyyyzz_0_yyyyyz_0, \
                             g_0_yyyyyzz_0_yyyyzz_0, \
                             g_0_yyyyyzz_0_yyyzzz_0, \
                             g_0_yyyyyzz_0_yyzzzz_0, \
                             g_0_yyyyyzz_0_yzzzzz_0, \
                             g_0_yyyyyzz_0_zzzzzz_0, \
                             g_0_yyyyzz_0_xxxxxx_0,  \
                             g_0_yyyyzz_0_xxxxxx_1,  \
                             g_0_yyyyzz_0_xxxxxz_0,  \
                             g_0_yyyyzz_0_xxxxxz_1,  \
                             g_0_yyyyzz_0_xxxxyz_0,  \
                             g_0_yyyyzz_0_xxxxyz_1,  \
                             g_0_yyyyzz_0_xxxxz_1,   \
                             g_0_yyyyzz_0_xxxxzz_0,  \
                             g_0_yyyyzz_0_xxxxzz_1,  \
                             g_0_yyyyzz_0_xxxyyz_0,  \
                             g_0_yyyyzz_0_xxxyyz_1,  \
                             g_0_yyyyzz_0_xxxyz_1,   \
                             g_0_yyyyzz_0_xxxyzz_0,  \
                             g_0_yyyyzz_0_xxxyzz_1,  \
                             g_0_yyyyzz_0_xxxzz_1,   \
                             g_0_yyyyzz_0_xxxzzz_0,  \
                             g_0_yyyyzz_0_xxxzzz_1,  \
                             g_0_yyyyzz_0_xxyyyz_0,  \
                             g_0_yyyyzz_0_xxyyyz_1,  \
                             g_0_yyyyzz_0_xxyyz_1,   \
                             g_0_yyyyzz_0_xxyyzz_0,  \
                             g_0_yyyyzz_0_xxyyzz_1,  \
                             g_0_yyyyzz_0_xxyzz_1,   \
                             g_0_yyyyzz_0_xxyzzz_0,  \
                             g_0_yyyyzz_0_xxyzzz_1,  \
                             g_0_yyyyzz_0_xxzzz_1,   \
                             g_0_yyyyzz_0_xxzzzz_0,  \
                             g_0_yyyyzz_0_xxzzzz_1,  \
                             g_0_yyyyzz_0_xyyyyz_0,  \
                             g_0_yyyyzz_0_xyyyyz_1,  \
                             g_0_yyyyzz_0_xyyyz_1,   \
                             g_0_yyyyzz_0_xyyyzz_0,  \
                             g_0_yyyyzz_0_xyyyzz_1,  \
                             g_0_yyyyzz_0_xyyzz_1,   \
                             g_0_yyyyzz_0_xyyzzz_0,  \
                             g_0_yyyyzz_0_xyyzzz_1,  \
                             g_0_yyyyzz_0_xyzzz_1,   \
                             g_0_yyyyzz_0_xyzzzz_0,  \
                             g_0_yyyyzz_0_xyzzzz_1,  \
                             g_0_yyyyzz_0_xzzzz_1,   \
                             g_0_yyyyzz_0_xzzzzz_0,  \
                             g_0_yyyyzz_0_xzzzzz_1,  \
                             g_0_yyyyzz_0_yyyyyz_0,  \
                             g_0_yyyyzz_0_yyyyyz_1,  \
                             g_0_yyyyzz_0_yyyyz_1,   \
                             g_0_yyyyzz_0_yyyyzz_0,  \
                             g_0_yyyyzz_0_yyyyzz_1,  \
                             g_0_yyyyzz_0_yyyzz_1,   \
                             g_0_yyyyzz_0_yyyzzz_0,  \
                             g_0_yyyyzz_0_yyyzzz_1,  \
                             g_0_yyyyzz_0_yyzzz_1,   \
                             g_0_yyyyzz_0_yyzzzz_0,  \
                             g_0_yyyyzz_0_yyzzzz_1,  \
                             g_0_yyyyzz_0_yzzzz_1,   \
                             g_0_yyyyzz_0_yzzzzz_0,  \
                             g_0_yyyyzz_0_yzzzzz_1,  \
                             g_0_yyyyzz_0_zzzzz_1,   \
                             g_0_yyyyzz_0_zzzzzz_0,  \
                             g_0_yyyyzz_0_zzzzzz_1,  \
                             g_0_yyyzz_0_xxxxxx_0,   \
                             g_0_yyyzz_0_xxxxxx_1,   \
                             g_0_yyyzz_0_xxxxxz_0,   \
                             g_0_yyyzz_0_xxxxxz_1,   \
                             g_0_yyyzz_0_xxxxyz_0,   \
                             g_0_yyyzz_0_xxxxyz_1,   \
                             g_0_yyyzz_0_xxxxzz_0,   \
                             g_0_yyyzz_0_xxxxzz_1,   \
                             g_0_yyyzz_0_xxxyyz_0,   \
                             g_0_yyyzz_0_xxxyyz_1,   \
                             g_0_yyyzz_0_xxxyzz_0,   \
                             g_0_yyyzz_0_xxxyzz_1,   \
                             g_0_yyyzz_0_xxxzzz_0,   \
                             g_0_yyyzz_0_xxxzzz_1,   \
                             g_0_yyyzz_0_xxyyyz_0,   \
                             g_0_yyyzz_0_xxyyyz_1,   \
                             g_0_yyyzz_0_xxyyzz_0,   \
                             g_0_yyyzz_0_xxyyzz_1,   \
                             g_0_yyyzz_0_xxyzzz_0,   \
                             g_0_yyyzz_0_xxyzzz_1,   \
                             g_0_yyyzz_0_xxzzzz_0,   \
                             g_0_yyyzz_0_xxzzzz_1,   \
                             g_0_yyyzz_0_xyyyyz_0,   \
                             g_0_yyyzz_0_xyyyyz_1,   \
                             g_0_yyyzz_0_xyyyzz_0,   \
                             g_0_yyyzz_0_xyyyzz_1,   \
                             g_0_yyyzz_0_xyyzzz_0,   \
                             g_0_yyyzz_0_xyyzzz_1,   \
                             g_0_yyyzz_0_xyzzzz_0,   \
                             g_0_yyyzz_0_xyzzzz_1,   \
                             g_0_yyyzz_0_xzzzzz_0,   \
                             g_0_yyyzz_0_xzzzzz_1,   \
                             g_0_yyyzz_0_yyyyyz_0,   \
                             g_0_yyyzz_0_yyyyyz_1,   \
                             g_0_yyyzz_0_yyyyzz_0,   \
                             g_0_yyyzz_0_yyyyzz_1,   \
                             g_0_yyyzz_0_yyyzzz_0,   \
                             g_0_yyyzz_0_yyyzzz_1,   \
                             g_0_yyyzz_0_yyzzzz_0,   \
                             g_0_yyyzz_0_yyzzzz_1,   \
                             g_0_yyyzz_0_yzzzzz_0,   \
                             g_0_yyyzz_0_yzzzzz_1,   \
                             g_0_yyyzz_0_zzzzzz_0,   \
                             g_0_yyyzz_0_zzzzzz_1,   \
                             wp_y,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyzz_0_xxxxxx_0[i] = 4.0 * g_0_yyyzz_0_xxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxxx_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xxxxxx_0[i] * pb_y + g_0_yyyyzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxxxy_0[i] = g_0_yyyyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxxy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxxxxy_0[i] * pb_z +
                                    g_0_yyyyyz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxxxxz_0[i] = 4.0 * g_0_yyyzz_0_xxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxxz_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xxxxxz_0[i] * pb_y + g_0_yyyyzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxxyy_0[i] = g_0_yyyyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxxyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxxxyy_0[i] * pb_z +
                                    g_0_yyyyyz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxxxyz_0[i] = 4.0 * g_0_yyyzz_0_xxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxyz_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxxyz_0[i] * pb_y + g_0_yyyyzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxxzz_0[i] = 4.0 * g_0_yyyzz_0_xxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxxzz_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xxxxzz_0[i] * pb_y + g_0_yyyyzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxyyy_0[i] = g_0_yyyyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxxyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxxyyy_0[i] * pb_z +
                                    g_0_yyyyyz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxxyyz_0[i] = 4.0 * g_0_yyyzz_0_xxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyyz_0[i] * pb_y + g_0_yyyyzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxyzz_0[i] = 4.0 * g_0_yyyzz_0_xxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxyzz_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxxyzz_0[i] * pb_y + g_0_yyyyzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxxzzz_0[i] = 4.0 * g_0_yyyzz_0_xxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxxzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xxxzzz_0[i] * pb_y + g_0_yyyyzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxyyyy_0[i] = g_0_yyyyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xxyyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xxyyyy_0[i] * pb_z +
                                    g_0_yyyyyz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xxyyyz_0[i] = 4.0 * g_0_yyyzz_0_xxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxyyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyyyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyyz_0[i] * pb_y + g_0_yyyyzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxyyzz_0[i] = 4.0 * g_0_yyyzz_0_xxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyyzz_0[i] * pb_y + g_0_yyyyzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxyzzz_0[i] = 4.0 * g_0_yyyzz_0_xxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxyzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xxyzzz_0[i] * pb_y + g_0_yyyyzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xxzzzz_0[i] = 4.0 * g_0_yyyzz_0_xxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xxzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xxzzzz_0[i] * pb_y + g_0_yyyyzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyyyyy_0[i] = g_0_yyyyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_xyyyyy_0[i] * pb_z +
                                    g_0_yyyyyz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_xyyyyz_0[i] = 4.0 * g_0_yyyzz_0_xyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyyyyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_yyyyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyyz_0[i] * pb_y + g_0_yyyyzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyyyzz_0[i] = 4.0 * g_0_yyyzz_0_xyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyyyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyyyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyyzz_0[i] * pb_y + g_0_yyyyzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyyzzz_0[i] = 4.0 * g_0_yyyzz_0_xyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyyzzz_0[i] * pb_y + g_0_yyyyzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xyzzzz_0[i] = 4.0 * g_0_yyyzz_0_xyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xyzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_xyzzzz_0[i] * pb_y + g_0_yyyyzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_xzzzzz_0[i] = 4.0 * g_0_yyyzz_0_xzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_xzzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xzzzzz_0[i] * pb_y + g_0_yyyyzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyyyyy_0[i] = g_0_yyyyy_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyyyy_0_yyyyyy_1[i] * fti_ab_0 + g_0_yyyyyz_0_yyyyyy_0[i] * pb_z +
                                    g_0_yyyyyz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyyyzz_0_yyyyyz_0[i] = 4.0 * g_0_yyyzz_0_yyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyyyyz_1[i] * fti_ab_0 +
                                    5.0 * g_0_yyyyzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyyyyz_0[i] * pb_y + g_0_yyyyzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyyyzz_0[i] = 4.0 * g_0_yyyzz_0_yyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyyyzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_yyyyzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyyyzz_0[i] * pb_y + g_0_yyyyzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyyzzz_0[i] = 4.0 * g_0_yyyzz_0_yyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyyzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyyyzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyyzzz_0[i] * pb_y + g_0_yyyyzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yyzzzz_0[i] = 4.0 * g_0_yyyzz_0_yyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yyzzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyyzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yyzzzz_0[i] * pb_y + g_0_yyyyzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_yzzzzz_0[i] = 4.0 * g_0_yyyzz_0_yzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_yzzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyzz_0_yzzzzz_0[i] * pb_y + g_0_yyyyzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzz_0_zzzzzz_0[i] = 4.0 * g_0_yyyzz_0_zzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_zzzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_zzzzzz_0[i] * pb_y + g_0_yyyyzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 868-896 components of targeted buffer : SKSI

    auto g_0_yyyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 868);

    auto g_0_yyyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 869);

    auto g_0_yyyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 870);

    auto g_0_yyyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 871);

    auto g_0_yyyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 872);

    auto g_0_yyyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 873);

    auto g_0_yyyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 874);

    auto g_0_yyyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 875);

    auto g_0_yyyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 876);

    auto g_0_yyyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 877);

    auto g_0_yyyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 878);

    auto g_0_yyyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 879);

    auto g_0_yyyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 880);

    auto g_0_yyyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 881);

    auto g_0_yyyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 882);

    auto g_0_yyyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 883);

    auto g_0_yyyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 884);

    auto g_0_yyyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 885);

    auto g_0_yyyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 886);

    auto g_0_yyyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 887);

    auto g_0_yyyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 888);

    auto g_0_yyyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 889);

    auto g_0_yyyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 890);

    auto g_0_yyyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 891);

    auto g_0_yyyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 892);

    auto g_0_yyyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 893);

    auto g_0_yyyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 894);

    auto g_0_yyyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 895);

#pragma omp simd aligned(g_0_yyyyz_0_xxxxxy_0,       \
                             g_0_yyyyz_0_xxxxxy_1,   \
                             g_0_yyyyz_0_xxxxyy_0,   \
                             g_0_yyyyz_0_xxxxyy_1,   \
                             g_0_yyyyz_0_xxxyyy_0,   \
                             g_0_yyyyz_0_xxxyyy_1,   \
                             g_0_yyyyz_0_xxyyyy_0,   \
                             g_0_yyyyz_0_xxyyyy_1,   \
                             g_0_yyyyz_0_xyyyyy_0,   \
                             g_0_yyyyz_0_xyyyyy_1,   \
                             g_0_yyyyz_0_yyyyyy_0,   \
                             g_0_yyyyz_0_yyyyyy_1,   \
                             g_0_yyyyzz_0_xxxxxy_0,  \
                             g_0_yyyyzz_0_xxxxxy_1,  \
                             g_0_yyyyzz_0_xxxxyy_0,  \
                             g_0_yyyyzz_0_xxxxyy_1,  \
                             g_0_yyyyzz_0_xxxyyy_0,  \
                             g_0_yyyyzz_0_xxxyyy_1,  \
                             g_0_yyyyzz_0_xxyyyy_0,  \
                             g_0_yyyyzz_0_xxyyyy_1,  \
                             g_0_yyyyzz_0_xyyyyy_0,  \
                             g_0_yyyyzz_0_xyyyyy_1,  \
                             g_0_yyyyzz_0_yyyyyy_0,  \
                             g_0_yyyyzz_0_yyyyyy_1,  \
                             g_0_yyyyzzz_0_xxxxxx_0, \
                             g_0_yyyyzzz_0_xxxxxy_0, \
                             g_0_yyyyzzz_0_xxxxxz_0, \
                             g_0_yyyyzzz_0_xxxxyy_0, \
                             g_0_yyyyzzz_0_xxxxyz_0, \
                             g_0_yyyyzzz_0_xxxxzz_0, \
                             g_0_yyyyzzz_0_xxxyyy_0, \
                             g_0_yyyyzzz_0_xxxyyz_0, \
                             g_0_yyyyzzz_0_xxxyzz_0, \
                             g_0_yyyyzzz_0_xxxzzz_0, \
                             g_0_yyyyzzz_0_xxyyyy_0, \
                             g_0_yyyyzzz_0_xxyyyz_0, \
                             g_0_yyyyzzz_0_xxyyzz_0, \
                             g_0_yyyyzzz_0_xxyzzz_0, \
                             g_0_yyyyzzz_0_xxzzzz_0, \
                             g_0_yyyyzzz_0_xyyyyy_0, \
                             g_0_yyyyzzz_0_xyyyyz_0, \
                             g_0_yyyyzzz_0_xyyyzz_0, \
                             g_0_yyyyzzz_0_xyyzzz_0, \
                             g_0_yyyyzzz_0_xyzzzz_0, \
                             g_0_yyyyzzz_0_xzzzzz_0, \
                             g_0_yyyyzzz_0_yyyyyy_0, \
                             g_0_yyyyzzz_0_yyyyyz_0, \
                             g_0_yyyyzzz_0_yyyyzz_0, \
                             g_0_yyyyzzz_0_yyyzzz_0, \
                             g_0_yyyyzzz_0_yyzzzz_0, \
                             g_0_yyyyzzz_0_yzzzzz_0, \
                             g_0_yyyyzzz_0_zzzzzz_0, \
                             g_0_yyyzzz_0_xxxxxx_0,  \
                             g_0_yyyzzz_0_xxxxxx_1,  \
                             g_0_yyyzzz_0_xxxxxz_0,  \
                             g_0_yyyzzz_0_xxxxxz_1,  \
                             g_0_yyyzzz_0_xxxxyz_0,  \
                             g_0_yyyzzz_0_xxxxyz_1,  \
                             g_0_yyyzzz_0_xxxxz_1,   \
                             g_0_yyyzzz_0_xxxxzz_0,  \
                             g_0_yyyzzz_0_xxxxzz_1,  \
                             g_0_yyyzzz_0_xxxyyz_0,  \
                             g_0_yyyzzz_0_xxxyyz_1,  \
                             g_0_yyyzzz_0_xxxyz_1,   \
                             g_0_yyyzzz_0_xxxyzz_0,  \
                             g_0_yyyzzz_0_xxxyzz_1,  \
                             g_0_yyyzzz_0_xxxzz_1,   \
                             g_0_yyyzzz_0_xxxzzz_0,  \
                             g_0_yyyzzz_0_xxxzzz_1,  \
                             g_0_yyyzzz_0_xxyyyz_0,  \
                             g_0_yyyzzz_0_xxyyyz_1,  \
                             g_0_yyyzzz_0_xxyyz_1,   \
                             g_0_yyyzzz_0_xxyyzz_0,  \
                             g_0_yyyzzz_0_xxyyzz_1,  \
                             g_0_yyyzzz_0_xxyzz_1,   \
                             g_0_yyyzzz_0_xxyzzz_0,  \
                             g_0_yyyzzz_0_xxyzzz_1,  \
                             g_0_yyyzzz_0_xxzzz_1,   \
                             g_0_yyyzzz_0_xxzzzz_0,  \
                             g_0_yyyzzz_0_xxzzzz_1,  \
                             g_0_yyyzzz_0_xyyyyz_0,  \
                             g_0_yyyzzz_0_xyyyyz_1,  \
                             g_0_yyyzzz_0_xyyyz_1,   \
                             g_0_yyyzzz_0_xyyyzz_0,  \
                             g_0_yyyzzz_0_xyyyzz_1,  \
                             g_0_yyyzzz_0_xyyzz_1,   \
                             g_0_yyyzzz_0_xyyzzz_0,  \
                             g_0_yyyzzz_0_xyyzzz_1,  \
                             g_0_yyyzzz_0_xyzzz_1,   \
                             g_0_yyyzzz_0_xyzzzz_0,  \
                             g_0_yyyzzz_0_xyzzzz_1,  \
                             g_0_yyyzzz_0_xzzzz_1,   \
                             g_0_yyyzzz_0_xzzzzz_0,  \
                             g_0_yyyzzz_0_xzzzzz_1,  \
                             g_0_yyyzzz_0_yyyyyz_0,  \
                             g_0_yyyzzz_0_yyyyyz_1,  \
                             g_0_yyyzzz_0_yyyyz_1,   \
                             g_0_yyyzzz_0_yyyyzz_0,  \
                             g_0_yyyzzz_0_yyyyzz_1,  \
                             g_0_yyyzzz_0_yyyzz_1,   \
                             g_0_yyyzzz_0_yyyzzz_0,  \
                             g_0_yyyzzz_0_yyyzzz_1,  \
                             g_0_yyyzzz_0_yyzzz_1,   \
                             g_0_yyyzzz_0_yyzzzz_0,  \
                             g_0_yyyzzz_0_yyzzzz_1,  \
                             g_0_yyyzzz_0_yzzzz_1,   \
                             g_0_yyyzzz_0_yzzzzz_0,  \
                             g_0_yyyzzz_0_yzzzzz_1,  \
                             g_0_yyyzzz_0_zzzzz_1,   \
                             g_0_yyyzzz_0_zzzzzz_0,  \
                             g_0_yyyzzz_0_zzzzzz_1,  \
                             g_0_yyzzz_0_xxxxxx_0,   \
                             g_0_yyzzz_0_xxxxxx_1,   \
                             g_0_yyzzz_0_xxxxxz_0,   \
                             g_0_yyzzz_0_xxxxxz_1,   \
                             g_0_yyzzz_0_xxxxyz_0,   \
                             g_0_yyzzz_0_xxxxyz_1,   \
                             g_0_yyzzz_0_xxxxzz_0,   \
                             g_0_yyzzz_0_xxxxzz_1,   \
                             g_0_yyzzz_0_xxxyyz_0,   \
                             g_0_yyzzz_0_xxxyyz_1,   \
                             g_0_yyzzz_0_xxxyzz_0,   \
                             g_0_yyzzz_0_xxxyzz_1,   \
                             g_0_yyzzz_0_xxxzzz_0,   \
                             g_0_yyzzz_0_xxxzzz_1,   \
                             g_0_yyzzz_0_xxyyyz_0,   \
                             g_0_yyzzz_0_xxyyyz_1,   \
                             g_0_yyzzz_0_xxyyzz_0,   \
                             g_0_yyzzz_0_xxyyzz_1,   \
                             g_0_yyzzz_0_xxyzzz_0,   \
                             g_0_yyzzz_0_xxyzzz_1,   \
                             g_0_yyzzz_0_xxzzzz_0,   \
                             g_0_yyzzz_0_xxzzzz_1,   \
                             g_0_yyzzz_0_xyyyyz_0,   \
                             g_0_yyzzz_0_xyyyyz_1,   \
                             g_0_yyzzz_0_xyyyzz_0,   \
                             g_0_yyzzz_0_xyyyzz_1,   \
                             g_0_yyzzz_0_xyyzzz_0,   \
                             g_0_yyzzz_0_xyyzzz_1,   \
                             g_0_yyzzz_0_xyzzzz_0,   \
                             g_0_yyzzz_0_xyzzzz_1,   \
                             g_0_yyzzz_0_xzzzzz_0,   \
                             g_0_yyzzz_0_xzzzzz_1,   \
                             g_0_yyzzz_0_yyyyyz_0,   \
                             g_0_yyzzz_0_yyyyyz_1,   \
                             g_0_yyzzz_0_yyyyzz_0,   \
                             g_0_yyzzz_0_yyyyzz_1,   \
                             g_0_yyzzz_0_yyyzzz_0,   \
                             g_0_yyzzz_0_yyyzzz_1,   \
                             g_0_yyzzz_0_yyzzzz_0,   \
                             g_0_yyzzz_0_yyzzzz_1,   \
                             g_0_yyzzz_0_yzzzzz_0,   \
                             g_0_yyzzz_0_yzzzzz_1,   \
                             g_0_yyzzz_0_zzzzzz_0,   \
                             g_0_yyzzz_0_zzzzzz_1,   \
                             wp_y,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzzz_0_xxxxxx_0[i] = 3.0 * g_0_yyzzz_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxxx_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xxxxxx_0[i] * pb_y + g_0_yyyzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxxxy_0[i] = 2.0 * g_0_yyyyz_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxxxxy_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xxxxxy_0[i] * pb_z + g_0_yyyyzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxxxxz_0[i] = 3.0 * g_0_yyzzz_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xxxxxz_0[i] * pb_y + g_0_yyyzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxxyy_0[i] = 2.0 * g_0_yyyyz_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxxxyy_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xxxxyy_0[i] * pb_z + g_0_yyyyzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxxxyz_0[i] = 3.0 * g_0_yyzzz_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxxyz_0[i] * pb_y + g_0_yyyzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxxzz_0[i] = 3.0 * g_0_yyzzz_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xxxxzz_0[i] * pb_y + g_0_yyyzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxyyy_0[i] = 2.0 * g_0_yyyyz_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxxyyy_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xxxyyy_0[i] * pb_z + g_0_yyyyzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxxyyz_0[i] = 3.0 * g_0_yyzzz_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyyz_0[i] * pb_y + g_0_yyyzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxyzz_0[i] = 3.0 * g_0_yyzzz_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxxyzz_0[i] * pb_y + g_0_yyyzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxxzzz_0[i] = 3.0 * g_0_yyzzz_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xxxzzz_0[i] * pb_y + g_0_yyyzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxyyyy_0[i] = 2.0 * g_0_yyyyz_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xxyyyy_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xxyyyy_0[i] * pb_z + g_0_yyyyzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xxyyyz_0[i] = 3.0 * g_0_yyzzz_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyyzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyyz_0[i] * pb_y + g_0_yyyzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxyyzz_0[i] = 3.0 * g_0_yyzzz_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyyzz_0[i] * pb_y + g_0_yyyzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxyzzz_0[i] = 3.0 * g_0_yyzzz_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xxyzzz_0[i] * pb_y + g_0_yyyzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xxzzzz_0[i] = 3.0 * g_0_yyzzz_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xxzzzz_0[i] * pb_y + g_0_yyyzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyyyyy_0[i] = 2.0 * g_0_yyyyz_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_xyyyyy_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_xyyyyy_0[i] * pb_z + g_0_yyyyzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_xyyyyz_0[i] = 3.0 * g_0_yyzzz_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_yyyzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyyz_0[i] * pb_y + g_0_yyyzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyyyzz_0[i] = 3.0 * g_0_yyzzz_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyyzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyyzz_0[i] * pb_y + g_0_yyyzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyyzzz_0[i] = 3.0 * g_0_yyzzz_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyyzzz_0[i] * pb_y + g_0_yyyzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xyzzzz_0[i] = 3.0 * g_0_yyzzz_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_xyzzzz_0[i] * pb_y + g_0_yyyzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_xzzzzz_0[i] = 3.0 * g_0_yyzzz_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xzzzzz_0[i] * pb_y + g_0_yyyzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyyyyy_0[i] = 2.0 * g_0_yyyyz_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyz_0_yyyyyy_1[i] * fti_ab_0 +
                                    g_0_yyyyzz_0_yyyyyy_0[i] * pb_z + g_0_yyyyzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyyzzz_0_yyyyyz_0[i] = 3.0 * g_0_yyzzz_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                    5.0 * g_0_yyyzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyyyyz_0[i] * pb_y + g_0_yyyzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyyyzz_0[i] = 3.0 * g_0_yyzzz_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_yyyzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyyyzz_0[i] * pb_y + g_0_yyyzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyyzzz_0[i] = 3.0 * g_0_yyzzz_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyyzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyyzzz_0[i] * pb_y + g_0_yyyzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yyzzzz_0[i] = 3.0 * g_0_yyzzz_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyyzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yyzzzz_0[i] * pb_y + g_0_yyyzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_yzzzzz_0[i] = 3.0 * g_0_yyzzz_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyzzz_0_yzzzzz_0[i] * pb_y + g_0_yyyzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzz_0_zzzzzz_0[i] = 3.0 * g_0_yyzzz_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_zzzzzz_0[i] * pb_y + g_0_yyyzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 896-924 components of targeted buffer : SKSI

    auto g_0_yyyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 896);

    auto g_0_yyyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 897);

    auto g_0_yyyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 898);

    auto g_0_yyyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 899);

    auto g_0_yyyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 900);

    auto g_0_yyyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 901);

    auto g_0_yyyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 902);

    auto g_0_yyyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 903);

    auto g_0_yyyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 904);

    auto g_0_yyyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 905);

    auto g_0_yyyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 906);

    auto g_0_yyyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 907);

    auto g_0_yyyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 908);

    auto g_0_yyyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 909);

    auto g_0_yyyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 910);

    auto g_0_yyyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 911);

    auto g_0_yyyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 912);

    auto g_0_yyyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 913);

    auto g_0_yyyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 914);

    auto g_0_yyyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 915);

    auto g_0_yyyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 916);

    auto g_0_yyyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 917);

    auto g_0_yyyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 918);

    auto g_0_yyyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 919);

    auto g_0_yyyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 920);

    auto g_0_yyyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 921);

    auto g_0_yyyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 922);

    auto g_0_yyyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 923);

#pragma omp simd aligned(g_0_yyyzz_0_xxxxxy_0,       \
                             g_0_yyyzz_0_xxxxxy_1,   \
                             g_0_yyyzz_0_xxxxyy_0,   \
                             g_0_yyyzz_0_xxxxyy_1,   \
                             g_0_yyyzz_0_xxxyyy_0,   \
                             g_0_yyyzz_0_xxxyyy_1,   \
                             g_0_yyyzz_0_xxyyyy_0,   \
                             g_0_yyyzz_0_xxyyyy_1,   \
                             g_0_yyyzz_0_xyyyyy_0,   \
                             g_0_yyyzz_0_xyyyyy_1,   \
                             g_0_yyyzz_0_yyyyyy_0,   \
                             g_0_yyyzz_0_yyyyyy_1,   \
                             g_0_yyyzzz_0_xxxxxy_0,  \
                             g_0_yyyzzz_0_xxxxxy_1,  \
                             g_0_yyyzzz_0_xxxxyy_0,  \
                             g_0_yyyzzz_0_xxxxyy_1,  \
                             g_0_yyyzzz_0_xxxyyy_0,  \
                             g_0_yyyzzz_0_xxxyyy_1,  \
                             g_0_yyyzzz_0_xxyyyy_0,  \
                             g_0_yyyzzz_0_xxyyyy_1,  \
                             g_0_yyyzzz_0_xyyyyy_0,  \
                             g_0_yyyzzz_0_xyyyyy_1,  \
                             g_0_yyyzzz_0_yyyyyy_0,  \
                             g_0_yyyzzz_0_yyyyyy_1,  \
                             g_0_yyyzzzz_0_xxxxxx_0, \
                             g_0_yyyzzzz_0_xxxxxy_0, \
                             g_0_yyyzzzz_0_xxxxxz_0, \
                             g_0_yyyzzzz_0_xxxxyy_0, \
                             g_0_yyyzzzz_0_xxxxyz_0, \
                             g_0_yyyzzzz_0_xxxxzz_0, \
                             g_0_yyyzzzz_0_xxxyyy_0, \
                             g_0_yyyzzzz_0_xxxyyz_0, \
                             g_0_yyyzzzz_0_xxxyzz_0, \
                             g_0_yyyzzzz_0_xxxzzz_0, \
                             g_0_yyyzzzz_0_xxyyyy_0, \
                             g_0_yyyzzzz_0_xxyyyz_0, \
                             g_0_yyyzzzz_0_xxyyzz_0, \
                             g_0_yyyzzzz_0_xxyzzz_0, \
                             g_0_yyyzzzz_0_xxzzzz_0, \
                             g_0_yyyzzzz_0_xyyyyy_0, \
                             g_0_yyyzzzz_0_xyyyyz_0, \
                             g_0_yyyzzzz_0_xyyyzz_0, \
                             g_0_yyyzzzz_0_xyyzzz_0, \
                             g_0_yyyzzzz_0_xyzzzz_0, \
                             g_0_yyyzzzz_0_xzzzzz_0, \
                             g_0_yyyzzzz_0_yyyyyy_0, \
                             g_0_yyyzzzz_0_yyyyyz_0, \
                             g_0_yyyzzzz_0_yyyyzz_0, \
                             g_0_yyyzzzz_0_yyyzzz_0, \
                             g_0_yyyzzzz_0_yyzzzz_0, \
                             g_0_yyyzzzz_0_yzzzzz_0, \
                             g_0_yyyzzzz_0_zzzzzz_0, \
                             g_0_yyzzzz_0_xxxxxx_0,  \
                             g_0_yyzzzz_0_xxxxxx_1,  \
                             g_0_yyzzzz_0_xxxxxz_0,  \
                             g_0_yyzzzz_0_xxxxxz_1,  \
                             g_0_yyzzzz_0_xxxxyz_0,  \
                             g_0_yyzzzz_0_xxxxyz_1,  \
                             g_0_yyzzzz_0_xxxxz_1,   \
                             g_0_yyzzzz_0_xxxxzz_0,  \
                             g_0_yyzzzz_0_xxxxzz_1,  \
                             g_0_yyzzzz_0_xxxyyz_0,  \
                             g_0_yyzzzz_0_xxxyyz_1,  \
                             g_0_yyzzzz_0_xxxyz_1,   \
                             g_0_yyzzzz_0_xxxyzz_0,  \
                             g_0_yyzzzz_0_xxxyzz_1,  \
                             g_0_yyzzzz_0_xxxzz_1,   \
                             g_0_yyzzzz_0_xxxzzz_0,  \
                             g_0_yyzzzz_0_xxxzzz_1,  \
                             g_0_yyzzzz_0_xxyyyz_0,  \
                             g_0_yyzzzz_0_xxyyyz_1,  \
                             g_0_yyzzzz_0_xxyyz_1,   \
                             g_0_yyzzzz_0_xxyyzz_0,  \
                             g_0_yyzzzz_0_xxyyzz_1,  \
                             g_0_yyzzzz_0_xxyzz_1,   \
                             g_0_yyzzzz_0_xxyzzz_0,  \
                             g_0_yyzzzz_0_xxyzzz_1,  \
                             g_0_yyzzzz_0_xxzzz_1,   \
                             g_0_yyzzzz_0_xxzzzz_0,  \
                             g_0_yyzzzz_0_xxzzzz_1,  \
                             g_0_yyzzzz_0_xyyyyz_0,  \
                             g_0_yyzzzz_0_xyyyyz_1,  \
                             g_0_yyzzzz_0_xyyyz_1,   \
                             g_0_yyzzzz_0_xyyyzz_0,  \
                             g_0_yyzzzz_0_xyyyzz_1,  \
                             g_0_yyzzzz_0_xyyzz_1,   \
                             g_0_yyzzzz_0_xyyzzz_0,  \
                             g_0_yyzzzz_0_xyyzzz_1,  \
                             g_0_yyzzzz_0_xyzzz_1,   \
                             g_0_yyzzzz_0_xyzzzz_0,  \
                             g_0_yyzzzz_0_xyzzzz_1,  \
                             g_0_yyzzzz_0_xzzzz_1,   \
                             g_0_yyzzzz_0_xzzzzz_0,  \
                             g_0_yyzzzz_0_xzzzzz_1,  \
                             g_0_yyzzzz_0_yyyyyz_0,  \
                             g_0_yyzzzz_0_yyyyyz_1,  \
                             g_0_yyzzzz_0_yyyyz_1,   \
                             g_0_yyzzzz_0_yyyyzz_0,  \
                             g_0_yyzzzz_0_yyyyzz_1,  \
                             g_0_yyzzzz_0_yyyzz_1,   \
                             g_0_yyzzzz_0_yyyzzz_0,  \
                             g_0_yyzzzz_0_yyyzzz_1,  \
                             g_0_yyzzzz_0_yyzzz_1,   \
                             g_0_yyzzzz_0_yyzzzz_0,  \
                             g_0_yyzzzz_0_yyzzzz_1,  \
                             g_0_yyzzzz_0_yzzzz_1,   \
                             g_0_yyzzzz_0_yzzzzz_0,  \
                             g_0_yyzzzz_0_yzzzzz_1,  \
                             g_0_yyzzzz_0_zzzzz_1,   \
                             g_0_yyzzzz_0_zzzzzz_0,  \
                             g_0_yyzzzz_0_zzzzzz_1,  \
                             g_0_yzzzz_0_xxxxxx_0,   \
                             g_0_yzzzz_0_xxxxxx_1,   \
                             g_0_yzzzz_0_xxxxxz_0,   \
                             g_0_yzzzz_0_xxxxxz_1,   \
                             g_0_yzzzz_0_xxxxyz_0,   \
                             g_0_yzzzz_0_xxxxyz_1,   \
                             g_0_yzzzz_0_xxxxzz_0,   \
                             g_0_yzzzz_0_xxxxzz_1,   \
                             g_0_yzzzz_0_xxxyyz_0,   \
                             g_0_yzzzz_0_xxxyyz_1,   \
                             g_0_yzzzz_0_xxxyzz_0,   \
                             g_0_yzzzz_0_xxxyzz_1,   \
                             g_0_yzzzz_0_xxxzzz_0,   \
                             g_0_yzzzz_0_xxxzzz_1,   \
                             g_0_yzzzz_0_xxyyyz_0,   \
                             g_0_yzzzz_0_xxyyyz_1,   \
                             g_0_yzzzz_0_xxyyzz_0,   \
                             g_0_yzzzz_0_xxyyzz_1,   \
                             g_0_yzzzz_0_xxyzzz_0,   \
                             g_0_yzzzz_0_xxyzzz_1,   \
                             g_0_yzzzz_0_xxzzzz_0,   \
                             g_0_yzzzz_0_xxzzzz_1,   \
                             g_0_yzzzz_0_xyyyyz_0,   \
                             g_0_yzzzz_0_xyyyyz_1,   \
                             g_0_yzzzz_0_xyyyzz_0,   \
                             g_0_yzzzz_0_xyyyzz_1,   \
                             g_0_yzzzz_0_xyyzzz_0,   \
                             g_0_yzzzz_0_xyyzzz_1,   \
                             g_0_yzzzz_0_xyzzzz_0,   \
                             g_0_yzzzz_0_xyzzzz_1,   \
                             g_0_yzzzz_0_xzzzzz_0,   \
                             g_0_yzzzz_0_xzzzzz_1,   \
                             g_0_yzzzz_0_yyyyyz_0,   \
                             g_0_yzzzz_0_yyyyyz_1,   \
                             g_0_yzzzz_0_yyyyzz_0,   \
                             g_0_yzzzz_0_yyyyzz_1,   \
                             g_0_yzzzz_0_yyyzzz_0,   \
                             g_0_yzzzz_0_yyyzzz_1,   \
                             g_0_yzzzz_0_yyzzzz_0,   \
                             g_0_yzzzz_0_yyzzzz_1,   \
                             g_0_yzzzz_0_yzzzzz_0,   \
                             g_0_yzzzz_0_yzzzzz_1,   \
                             g_0_yzzzz_0_zzzzzz_0,   \
                             g_0_yzzzz_0_zzzzzz_1,   \
                             wp_y,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzzz_0_xxxxxx_0[i] = 2.0 * g_0_yzzzz_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxxx_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xxxxxx_0[i] * pb_y + g_0_yyzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxxxy_0[i] = 3.0 * g_0_yyyzz_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxxxxy_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xxxxxy_0[i] * pb_z + g_0_yyyzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxxxxz_0[i] = 2.0 * g_0_yzzzz_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xxxxxz_0[i] * pb_y + g_0_yyzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxxyy_0[i] = 3.0 * g_0_yyyzz_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxxxyy_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xxxxyy_0[i] * pb_z + g_0_yyyzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxxxyz_0[i] = 2.0 * g_0_yzzzz_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxxyz_0[i] * pb_y + g_0_yyzzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxxzz_0[i] = 2.0 * g_0_yzzzz_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xxxxzz_0[i] * pb_y + g_0_yyzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxyyy_0[i] = 3.0 * g_0_yyyzz_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxxyyy_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xxxyyy_0[i] * pb_z + g_0_yyyzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxxyyz_0[i] = 2.0 * g_0_yzzzz_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyyz_0[i] * pb_y + g_0_yyzzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxyzz_0[i] = 2.0 * g_0_yzzzz_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxxyzz_0[i] * pb_y + g_0_yyzzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxxzzz_0[i] = 2.0 * g_0_yzzzz_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xxxzzz_0[i] * pb_y + g_0_yyzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxyyyy_0[i] = 3.0 * g_0_yyyzz_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xxyyyy_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xxyyyy_0[i] * pb_z + g_0_yyyzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xxyyyz_0[i] = 2.0 * g_0_yzzzz_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyyz_0[i] * pb_y + g_0_yyzzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxyyzz_0[i] = 2.0 * g_0_yzzzz_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyyzz_0[i] * pb_y + g_0_yyzzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxyzzz_0[i] = 2.0 * g_0_yzzzz_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xxyzzz_0[i] * pb_y + g_0_yyzzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xxzzzz_0[i] = 2.0 * g_0_yzzzz_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xxzzzz_0[i] * pb_y + g_0_yyzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyyyyy_0[i] = 3.0 * g_0_yyyzz_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_xyyyyy_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_xyyyyy_0[i] * pb_z + g_0_yyyzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_xyyyyz_0[i] = 2.0 * g_0_yzzzz_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_yyzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyyz_0[i] * pb_y + g_0_yyzzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyyyzz_0[i] = 2.0 * g_0_yzzzz_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyyzz_0[i] * pb_y + g_0_yyzzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyyzzz_0[i] = 2.0 * g_0_yzzzz_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyyzzz_0[i] * pb_y + g_0_yyzzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xyzzzz_0[i] = 2.0 * g_0_yzzzz_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_xyzzzz_0[i] * pb_y + g_0_yyzzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_xzzzzz_0[i] = 2.0 * g_0_yzzzz_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xzzzzz_0[i] * pb_y + g_0_yyzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyyyyy_0[i] = 3.0 * g_0_yyyzz_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyzz_0_yyyyyy_1[i] * fti_ab_0 +
                                    g_0_yyyzzz_0_yyyyyy_0[i] * pb_z + g_0_yyyzzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyzzzz_0_yyyyyz_0[i] = 2.0 * g_0_yzzzz_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                    5.0 * g_0_yyzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyyyyz_0[i] * pb_y + g_0_yyzzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyyyzz_0[i] = 2.0 * g_0_yzzzz_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_yyzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyyyzz_0[i] * pb_y + g_0_yyzzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyyzzz_0[i] = 2.0 * g_0_yzzzz_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yyzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyyzzz_0[i] * pb_y + g_0_yyzzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yyzzzz_0[i] = 2.0 * g_0_yzzzz_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yyzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yyzzzz_0[i] * pb_y + g_0_yyzzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_yzzzzz_0[i] = 2.0 * g_0_yzzzz_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyzzzz_0_yzzzzz_0[i] * pb_y + g_0_yyzzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzz_0_zzzzzz_0[i] = 2.0 * g_0_yzzzz_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_zzzzzz_0[i] * pb_y + g_0_yyzzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 924-952 components of targeted buffer : SKSI

    auto g_0_yyzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 924);

    auto g_0_yyzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 925);

    auto g_0_yyzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 926);

    auto g_0_yyzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 927);

    auto g_0_yyzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 928);

    auto g_0_yyzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 929);

    auto g_0_yyzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 930);

    auto g_0_yyzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 931);

    auto g_0_yyzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 932);

    auto g_0_yyzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 933);

    auto g_0_yyzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 934);

    auto g_0_yyzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 935);

    auto g_0_yyzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 936);

    auto g_0_yyzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 937);

    auto g_0_yyzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 938);

    auto g_0_yyzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 939);

    auto g_0_yyzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 940);

    auto g_0_yyzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 941);

    auto g_0_yyzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 942);

    auto g_0_yyzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 943);

    auto g_0_yyzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 944);

    auto g_0_yyzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 945);

    auto g_0_yyzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 946);

    auto g_0_yyzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 947);

    auto g_0_yyzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 948);

    auto g_0_yyzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 949);

    auto g_0_yyzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 950);

    auto g_0_yyzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 951);

#pragma omp simd aligned(g_0_yyzzz_0_xxxxxy_0,       \
                             g_0_yyzzz_0_xxxxxy_1,   \
                             g_0_yyzzz_0_xxxxyy_0,   \
                             g_0_yyzzz_0_xxxxyy_1,   \
                             g_0_yyzzz_0_xxxyyy_0,   \
                             g_0_yyzzz_0_xxxyyy_1,   \
                             g_0_yyzzz_0_xxyyyy_0,   \
                             g_0_yyzzz_0_xxyyyy_1,   \
                             g_0_yyzzz_0_xyyyyy_0,   \
                             g_0_yyzzz_0_xyyyyy_1,   \
                             g_0_yyzzz_0_yyyyyy_0,   \
                             g_0_yyzzz_0_yyyyyy_1,   \
                             g_0_yyzzzz_0_xxxxxy_0,  \
                             g_0_yyzzzz_0_xxxxxy_1,  \
                             g_0_yyzzzz_0_xxxxyy_0,  \
                             g_0_yyzzzz_0_xxxxyy_1,  \
                             g_0_yyzzzz_0_xxxyyy_0,  \
                             g_0_yyzzzz_0_xxxyyy_1,  \
                             g_0_yyzzzz_0_xxyyyy_0,  \
                             g_0_yyzzzz_0_xxyyyy_1,  \
                             g_0_yyzzzz_0_xyyyyy_0,  \
                             g_0_yyzzzz_0_xyyyyy_1,  \
                             g_0_yyzzzz_0_yyyyyy_0,  \
                             g_0_yyzzzz_0_yyyyyy_1,  \
                             g_0_yyzzzzz_0_xxxxxx_0, \
                             g_0_yyzzzzz_0_xxxxxy_0, \
                             g_0_yyzzzzz_0_xxxxxz_0, \
                             g_0_yyzzzzz_0_xxxxyy_0, \
                             g_0_yyzzzzz_0_xxxxyz_0, \
                             g_0_yyzzzzz_0_xxxxzz_0, \
                             g_0_yyzzzzz_0_xxxyyy_0, \
                             g_0_yyzzzzz_0_xxxyyz_0, \
                             g_0_yyzzzzz_0_xxxyzz_0, \
                             g_0_yyzzzzz_0_xxxzzz_0, \
                             g_0_yyzzzzz_0_xxyyyy_0, \
                             g_0_yyzzzzz_0_xxyyyz_0, \
                             g_0_yyzzzzz_0_xxyyzz_0, \
                             g_0_yyzzzzz_0_xxyzzz_0, \
                             g_0_yyzzzzz_0_xxzzzz_0, \
                             g_0_yyzzzzz_0_xyyyyy_0, \
                             g_0_yyzzzzz_0_xyyyyz_0, \
                             g_0_yyzzzzz_0_xyyyzz_0, \
                             g_0_yyzzzzz_0_xyyzzz_0, \
                             g_0_yyzzzzz_0_xyzzzz_0, \
                             g_0_yyzzzzz_0_xzzzzz_0, \
                             g_0_yyzzzzz_0_yyyyyy_0, \
                             g_0_yyzzzzz_0_yyyyyz_0, \
                             g_0_yyzzzzz_0_yyyyzz_0, \
                             g_0_yyzzzzz_0_yyyzzz_0, \
                             g_0_yyzzzzz_0_yyzzzz_0, \
                             g_0_yyzzzzz_0_yzzzzz_0, \
                             g_0_yyzzzzz_0_zzzzzz_0, \
                             g_0_yzzzzz_0_xxxxxx_0,  \
                             g_0_yzzzzz_0_xxxxxx_1,  \
                             g_0_yzzzzz_0_xxxxxz_0,  \
                             g_0_yzzzzz_0_xxxxxz_1,  \
                             g_0_yzzzzz_0_xxxxyz_0,  \
                             g_0_yzzzzz_0_xxxxyz_1,  \
                             g_0_yzzzzz_0_xxxxz_1,   \
                             g_0_yzzzzz_0_xxxxzz_0,  \
                             g_0_yzzzzz_0_xxxxzz_1,  \
                             g_0_yzzzzz_0_xxxyyz_0,  \
                             g_0_yzzzzz_0_xxxyyz_1,  \
                             g_0_yzzzzz_0_xxxyz_1,   \
                             g_0_yzzzzz_0_xxxyzz_0,  \
                             g_0_yzzzzz_0_xxxyzz_1,  \
                             g_0_yzzzzz_0_xxxzz_1,   \
                             g_0_yzzzzz_0_xxxzzz_0,  \
                             g_0_yzzzzz_0_xxxzzz_1,  \
                             g_0_yzzzzz_0_xxyyyz_0,  \
                             g_0_yzzzzz_0_xxyyyz_1,  \
                             g_0_yzzzzz_0_xxyyz_1,   \
                             g_0_yzzzzz_0_xxyyzz_0,  \
                             g_0_yzzzzz_0_xxyyzz_1,  \
                             g_0_yzzzzz_0_xxyzz_1,   \
                             g_0_yzzzzz_0_xxyzzz_0,  \
                             g_0_yzzzzz_0_xxyzzz_1,  \
                             g_0_yzzzzz_0_xxzzz_1,   \
                             g_0_yzzzzz_0_xxzzzz_0,  \
                             g_0_yzzzzz_0_xxzzzz_1,  \
                             g_0_yzzzzz_0_xyyyyz_0,  \
                             g_0_yzzzzz_0_xyyyyz_1,  \
                             g_0_yzzzzz_0_xyyyz_1,   \
                             g_0_yzzzzz_0_xyyyzz_0,  \
                             g_0_yzzzzz_0_xyyyzz_1,  \
                             g_0_yzzzzz_0_xyyzz_1,   \
                             g_0_yzzzzz_0_xyyzzz_0,  \
                             g_0_yzzzzz_0_xyyzzz_1,  \
                             g_0_yzzzzz_0_xyzzz_1,   \
                             g_0_yzzzzz_0_xyzzzz_0,  \
                             g_0_yzzzzz_0_xyzzzz_1,  \
                             g_0_yzzzzz_0_xzzzz_1,   \
                             g_0_yzzzzz_0_xzzzzz_0,  \
                             g_0_yzzzzz_0_xzzzzz_1,  \
                             g_0_yzzzzz_0_yyyyyz_0,  \
                             g_0_yzzzzz_0_yyyyyz_1,  \
                             g_0_yzzzzz_0_yyyyz_1,   \
                             g_0_yzzzzz_0_yyyyzz_0,  \
                             g_0_yzzzzz_0_yyyyzz_1,  \
                             g_0_yzzzzz_0_yyyzz_1,   \
                             g_0_yzzzzz_0_yyyzzz_0,  \
                             g_0_yzzzzz_0_yyyzzz_1,  \
                             g_0_yzzzzz_0_yyzzz_1,   \
                             g_0_yzzzzz_0_yyzzzz_0,  \
                             g_0_yzzzzz_0_yyzzzz_1,  \
                             g_0_yzzzzz_0_yzzzz_1,   \
                             g_0_yzzzzz_0_yzzzzz_0,  \
                             g_0_yzzzzz_0_yzzzzz_1,  \
                             g_0_yzzzzz_0_zzzzz_1,   \
                             g_0_yzzzzz_0_zzzzzz_0,  \
                             g_0_yzzzzz_0_zzzzzz_1,  \
                             g_0_zzzzz_0_xxxxxx_0,   \
                             g_0_zzzzz_0_xxxxxx_1,   \
                             g_0_zzzzz_0_xxxxxz_0,   \
                             g_0_zzzzz_0_xxxxxz_1,   \
                             g_0_zzzzz_0_xxxxyz_0,   \
                             g_0_zzzzz_0_xxxxyz_1,   \
                             g_0_zzzzz_0_xxxxzz_0,   \
                             g_0_zzzzz_0_xxxxzz_1,   \
                             g_0_zzzzz_0_xxxyyz_0,   \
                             g_0_zzzzz_0_xxxyyz_1,   \
                             g_0_zzzzz_0_xxxyzz_0,   \
                             g_0_zzzzz_0_xxxyzz_1,   \
                             g_0_zzzzz_0_xxxzzz_0,   \
                             g_0_zzzzz_0_xxxzzz_1,   \
                             g_0_zzzzz_0_xxyyyz_0,   \
                             g_0_zzzzz_0_xxyyyz_1,   \
                             g_0_zzzzz_0_xxyyzz_0,   \
                             g_0_zzzzz_0_xxyyzz_1,   \
                             g_0_zzzzz_0_xxyzzz_0,   \
                             g_0_zzzzz_0_xxyzzz_1,   \
                             g_0_zzzzz_0_xxzzzz_0,   \
                             g_0_zzzzz_0_xxzzzz_1,   \
                             g_0_zzzzz_0_xyyyyz_0,   \
                             g_0_zzzzz_0_xyyyyz_1,   \
                             g_0_zzzzz_0_xyyyzz_0,   \
                             g_0_zzzzz_0_xyyyzz_1,   \
                             g_0_zzzzz_0_xyyzzz_0,   \
                             g_0_zzzzz_0_xyyzzz_1,   \
                             g_0_zzzzz_0_xyzzzz_0,   \
                             g_0_zzzzz_0_xyzzzz_1,   \
                             g_0_zzzzz_0_xzzzzz_0,   \
                             g_0_zzzzz_0_xzzzzz_1,   \
                             g_0_zzzzz_0_yyyyyz_0,   \
                             g_0_zzzzz_0_yyyyyz_1,   \
                             g_0_zzzzz_0_yyyyzz_0,   \
                             g_0_zzzzz_0_yyyyzz_1,   \
                             g_0_zzzzz_0_yyyzzz_0,   \
                             g_0_zzzzz_0_yyyzzz_1,   \
                             g_0_zzzzz_0_yyzzzz_0,   \
                             g_0_zzzzz_0_yyzzzz_1,   \
                             g_0_zzzzz_0_yzzzzz_0,   \
                             g_0_zzzzz_0_yzzzzz_1,   \
                             g_0_zzzzz_0_zzzzzz_0,   \
                             g_0_zzzzz_0_zzzzzz_1,   \
                             wp_y,                   \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzzz_0_xxxxxx_0[i] = g_0_zzzzz_0_xxxxxx_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxxxx_0[i] * pb_y +
                                    g_0_yzzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxxxy_0[i] = 4.0 * g_0_yyzzz_0_xxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxxxxy_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xxxxxy_0[i] * pb_z + g_0_yyzzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxxxxz_0[i] = g_0_zzzzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxxxz_0[i] * pb_y +
                                    g_0_yzzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxxyy_0[i] = 4.0 * g_0_yyzzz_0_xxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxxxyy_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xxxxyy_0[i] * pb_z + g_0_yyzzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxxxyz_0[i] = g_0_zzzzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxyz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxxz_1[i] * fi_abcd_0 +
                                    g_0_yzzzzz_0_xxxxyz_0[i] * pb_y + g_0_yzzzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxxzz_0[i] = g_0_zzzzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxxzz_0[i] * pb_y +
                                    g_0_yzzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxyyy_0[i] = 4.0 * g_0_yyzzz_0_xxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxxyyy_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xxxyyy_0[i] * pb_z + g_0_yyzzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxxyyz_0[i] = g_0_zzzzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxxyyz_0[i] * pb_y + g_0_yzzzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxyzz_0[i] = g_0_zzzzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxyzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxzz_1[i] * fi_abcd_0 +
                                    g_0_yzzzzz_0_xxxyzz_0[i] * pb_y + g_0_yzzzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxxzzz_0[i] = g_0_zzzzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxxzzz_0[i] * pb_y +
                                    g_0_yzzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxyyyy_0[i] = 4.0 * g_0_yyzzz_0_xxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xxyyyy_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xxyyyy_0[i] * pb_z + g_0_yyzzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xxyyyz_0[i] = g_0_zzzzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyyz_0[i] * pb_y + g_0_yzzzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxyyzz_0[i] = g_0_zzzzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xxyyzz_0[i] * pb_y + g_0_yzzzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxyzzz_0[i] = g_0_zzzzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxyzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxzzz_1[i] * fi_abcd_0 +
                                    g_0_yzzzzz_0_xxyzzz_0[i] * pb_y + g_0_yzzzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xxzzzz_0[i] = g_0_zzzzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xxzzzz_0[i] * pb_y +
                                    g_0_yzzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyyyyy_0[i] = 4.0 * g_0_yyzzz_0_xyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_xyyyyy_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_xyyyyy_0[i] * pb_z + g_0_yyzzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_xyyyyz_0[i] = g_0_zzzzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                    4.0 * g_0_yzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyyz_0[i] * pb_y + g_0_yzzzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyyyzz_0[i] = g_0_zzzzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyyzz_0[i] * pb_y + g_0_yzzzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyyzzz_0[i] = g_0_zzzzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_xyyzzz_0[i] * pb_y + g_0_yzzzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xyzzzz_0[i] = g_0_zzzzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xzzzz_1[i] * fi_abcd_0 +
                                    g_0_yzzzzz_0_xyzzzz_0[i] * pb_y + g_0_yzzzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_xzzzzz_0[i] = g_0_zzzzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_xzzzzz_0[i] * pb_y +
                                    g_0_yzzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyyyyy_0[i] = 4.0 * g_0_yyzzz_0_yyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyzzz_0_yyyyyy_1[i] * fti_ab_0 +
                                    g_0_yyzzzz_0_yyyyyy_0[i] * pb_z + g_0_yyzzzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyzzzzz_0_yyyyyz_0[i] = g_0_zzzzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                    5.0 * g_0_yzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yyyyyz_0[i] * pb_y + g_0_yzzzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyyyzz_0[i] = g_0_zzzzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_yzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yyyyzz_0[i] * pb_y + g_0_yzzzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyyzzz_0[i] = g_0_zzzzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_yzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yyyzzz_0[i] * pb_y + g_0_yzzzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yyzzzz_0[i] = g_0_zzzzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_yzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yzzzzz_0_yyzzzz_0[i] * pb_y + g_0_yzzzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_yzzzzz_0[i] = g_0_zzzzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_zzzzz_1[i] * fi_abcd_0 +
                                    g_0_yzzzzz_0_yzzzzz_0[i] * pb_y + g_0_yzzzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzz_0_zzzzzz_0[i] = g_0_zzzzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_zzzzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_yzzzzz_0_zzzzzz_0[i] * pb_y +
                                    g_0_yzzzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 952-980 components of targeted buffer : SKSI

    auto g_0_yzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 952);

    auto g_0_yzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 953);

    auto g_0_yzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 954);

    auto g_0_yzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 955);

    auto g_0_yzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 956);

    auto g_0_yzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 957);

    auto g_0_yzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 958);

    auto g_0_yzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 959);

    auto g_0_yzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 960);

    auto g_0_yzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 961);

    auto g_0_yzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 962);

    auto g_0_yzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 963);

    auto g_0_yzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 964);

    auto g_0_yzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 965);

    auto g_0_yzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 966);

    auto g_0_yzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 967);

    auto g_0_yzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 968);

    auto g_0_yzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 969);

    auto g_0_yzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 970);

    auto g_0_yzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 971);

    auto g_0_yzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 972);

    auto g_0_yzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 973);

    auto g_0_yzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 974);

    auto g_0_yzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 975);

    auto g_0_yzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 976);

    auto g_0_yzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 977);

    auto g_0_yzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 978);

    auto g_0_yzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 979);

#pragma omp simd aligned(g_0_yzzzzzz_0_xxxxxx_0,     \
                             g_0_yzzzzzz_0_xxxxxy_0, \
                             g_0_yzzzzzz_0_xxxxxz_0, \
                             g_0_yzzzzzz_0_xxxxyy_0, \
                             g_0_yzzzzzz_0_xxxxyz_0, \
                             g_0_yzzzzzz_0_xxxxzz_0, \
                             g_0_yzzzzzz_0_xxxyyy_0, \
                             g_0_yzzzzzz_0_xxxyyz_0, \
                             g_0_yzzzzzz_0_xxxyzz_0, \
                             g_0_yzzzzzz_0_xxxzzz_0, \
                             g_0_yzzzzzz_0_xxyyyy_0, \
                             g_0_yzzzzzz_0_xxyyyz_0, \
                             g_0_yzzzzzz_0_xxyyzz_0, \
                             g_0_yzzzzzz_0_xxyzzz_0, \
                             g_0_yzzzzzz_0_xxzzzz_0, \
                             g_0_yzzzzzz_0_xyyyyy_0, \
                             g_0_yzzzzzz_0_xyyyyz_0, \
                             g_0_yzzzzzz_0_xyyyzz_0, \
                             g_0_yzzzzzz_0_xyyzzz_0, \
                             g_0_yzzzzzz_0_xyzzzz_0, \
                             g_0_yzzzzzz_0_xzzzzz_0, \
                             g_0_yzzzzzz_0_yyyyyy_0, \
                             g_0_yzzzzzz_0_yyyyyz_0, \
                             g_0_yzzzzzz_0_yyyyzz_0, \
                             g_0_yzzzzzz_0_yyyzzz_0, \
                             g_0_yzzzzzz_0_yyzzzz_0, \
                             g_0_yzzzzzz_0_yzzzzz_0, \
                             g_0_yzzzzzz_0_zzzzzz_0, \
                             g_0_zzzzzz_0_xxxxx_1,   \
                             g_0_zzzzzz_0_xxxxxx_0,  \
                             g_0_zzzzzz_0_xxxxxx_1,  \
                             g_0_zzzzzz_0_xxxxxy_0,  \
                             g_0_zzzzzz_0_xxxxxy_1,  \
                             g_0_zzzzzz_0_xxxxxz_0,  \
                             g_0_zzzzzz_0_xxxxxz_1,  \
                             g_0_zzzzzz_0_xxxxy_1,   \
                             g_0_zzzzzz_0_xxxxyy_0,  \
                             g_0_zzzzzz_0_xxxxyy_1,  \
                             g_0_zzzzzz_0_xxxxyz_0,  \
                             g_0_zzzzzz_0_xxxxyz_1,  \
                             g_0_zzzzzz_0_xxxxz_1,   \
                             g_0_zzzzzz_0_xxxxzz_0,  \
                             g_0_zzzzzz_0_xxxxzz_1,  \
                             g_0_zzzzzz_0_xxxyy_1,   \
                             g_0_zzzzzz_0_xxxyyy_0,  \
                             g_0_zzzzzz_0_xxxyyy_1,  \
                             g_0_zzzzzz_0_xxxyyz_0,  \
                             g_0_zzzzzz_0_xxxyyz_1,  \
                             g_0_zzzzzz_0_xxxyz_1,   \
                             g_0_zzzzzz_0_xxxyzz_0,  \
                             g_0_zzzzzz_0_xxxyzz_1,  \
                             g_0_zzzzzz_0_xxxzz_1,   \
                             g_0_zzzzzz_0_xxxzzz_0,  \
                             g_0_zzzzzz_0_xxxzzz_1,  \
                             g_0_zzzzzz_0_xxyyy_1,   \
                             g_0_zzzzzz_0_xxyyyy_0,  \
                             g_0_zzzzzz_0_xxyyyy_1,  \
                             g_0_zzzzzz_0_xxyyyz_0,  \
                             g_0_zzzzzz_0_xxyyyz_1,  \
                             g_0_zzzzzz_0_xxyyz_1,   \
                             g_0_zzzzzz_0_xxyyzz_0,  \
                             g_0_zzzzzz_0_xxyyzz_1,  \
                             g_0_zzzzzz_0_xxyzz_1,   \
                             g_0_zzzzzz_0_xxyzzz_0,  \
                             g_0_zzzzzz_0_xxyzzz_1,  \
                             g_0_zzzzzz_0_xxzzz_1,   \
                             g_0_zzzzzz_0_xxzzzz_0,  \
                             g_0_zzzzzz_0_xxzzzz_1,  \
                             g_0_zzzzzz_0_xyyyy_1,   \
                             g_0_zzzzzz_0_xyyyyy_0,  \
                             g_0_zzzzzz_0_xyyyyy_1,  \
                             g_0_zzzzzz_0_xyyyyz_0,  \
                             g_0_zzzzzz_0_xyyyyz_1,  \
                             g_0_zzzzzz_0_xyyyz_1,   \
                             g_0_zzzzzz_0_xyyyzz_0,  \
                             g_0_zzzzzz_0_xyyyzz_1,  \
                             g_0_zzzzzz_0_xyyzz_1,   \
                             g_0_zzzzzz_0_xyyzzz_0,  \
                             g_0_zzzzzz_0_xyyzzz_1,  \
                             g_0_zzzzzz_0_xyzzz_1,   \
                             g_0_zzzzzz_0_xyzzzz_0,  \
                             g_0_zzzzzz_0_xyzzzz_1,  \
                             g_0_zzzzzz_0_xzzzz_1,   \
                             g_0_zzzzzz_0_xzzzzz_0,  \
                             g_0_zzzzzz_0_xzzzzz_1,  \
                             g_0_zzzzzz_0_yyyyy_1,   \
                             g_0_zzzzzz_0_yyyyyy_0,  \
                             g_0_zzzzzz_0_yyyyyy_1,  \
                             g_0_zzzzzz_0_yyyyyz_0,  \
                             g_0_zzzzzz_0_yyyyyz_1,  \
                             g_0_zzzzzz_0_yyyyz_1,   \
                             g_0_zzzzzz_0_yyyyzz_0,  \
                             g_0_zzzzzz_0_yyyyzz_1,  \
                             g_0_zzzzzz_0_yyyzz_1,   \
                             g_0_zzzzzz_0_yyyzzz_0,  \
                             g_0_zzzzzz_0_yyyzzz_1,  \
                             g_0_zzzzzz_0_yyzzz_1,   \
                             g_0_zzzzzz_0_yyzzzz_0,  \
                             g_0_zzzzzz_0_yyzzzz_1,  \
                             g_0_zzzzzz_0_yzzzz_1,   \
                             g_0_zzzzzz_0_yzzzzz_0,  \
                             g_0_zzzzzz_0_yzzzzz_1,  \
                             g_0_zzzzzz_0_zzzzz_1,   \
                             g_0_zzzzzz_0_zzzzzz_0,  \
                             g_0_zzzzzz_0_zzzzzz_1,  \
                             wp_y,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzz_0_xxxxxx_0[i] = g_0_zzzzzz_0_xxxxxx_0[i] * pb_y + g_0_zzzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxxy_0[i] = g_0_zzzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxy_0[i] * pb_y + g_0_zzzzzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxxz_0[i] = g_0_zzzzzz_0_xxxxxz_0[i] * pb_y + g_0_zzzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxyy_0[i] = 2.0 * g_0_zzzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyy_0[i] * pb_y + g_0_zzzzzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxyz_0[i] = g_0_zzzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyz_0[i] * pb_y + g_0_zzzzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxxzz_0[i] = g_0_zzzzzz_0_xxxxzz_0[i] * pb_y + g_0_zzzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxyyy_0[i] = 3.0 * g_0_zzzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyy_0[i] * pb_y + g_0_zzzzzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxyyz_0[i] = 2.0 * g_0_zzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyz_0[i] * pb_y + g_0_zzzzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxyzz_0[i] = g_0_zzzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyzz_0[i] * pb_y + g_0_zzzzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxxzzz_0[i] = g_0_zzzzzz_0_xxxzzz_0[i] * pb_y + g_0_zzzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyyyy_0[i] = 4.0 * g_0_zzzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyy_0[i] * pb_y + g_0_zzzzzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyyyz_0[i] = 3.0 * g_0_zzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyz_0[i] * pb_y + g_0_zzzzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyyzz_0[i] = 2.0 * g_0_zzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyzz_0[i] * pb_y + g_0_zzzzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxyzzz_0[i] = g_0_zzzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyzzz_0[i] * pb_y + g_0_zzzzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xxzzzz_0[i] = g_0_zzzzzz_0_xxzzzz_0[i] * pb_y + g_0_zzzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyyyy_0[i] = 5.0 * g_0_zzzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyy_0[i] * pb_y + g_0_zzzzzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyyyz_0[i] = 4.0 * g_0_zzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyz_0[i] * pb_y + g_0_zzzzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyyzz_0[i] = 3.0 * g_0_zzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyzz_0[i] * pb_y + g_0_zzzzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyyzzz_0[i] = 2.0 * g_0_zzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyzzz_0[i] * pb_y + g_0_zzzzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xyzzzz_0[i] = g_0_zzzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyzzzz_0[i] * pb_y + g_0_zzzzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_xzzzzz_0[i] = g_0_zzzzzz_0_xzzzzz_0[i] * pb_y + g_0_zzzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyyyy_0[i] = 6.0 * g_0_zzzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyyy_0[i] * pb_y + g_0_zzzzzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyyyz_0[i] = 5.0 * g_0_zzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyyz_0[i] * pb_y + g_0_zzzzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyyzz_0[i] = 4.0 * g_0_zzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyzz_0[i] * pb_y + g_0_zzzzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyyzzz_0[i] = 3.0 * g_0_zzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyzzz_0[i] * pb_y + g_0_zzzzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yyzzzz_0[i] = 2.0 * g_0_zzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyzzzz_0[i] * pb_y + g_0_zzzzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_yzzzzz_0[i] = g_0_zzzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yzzzzz_0[i] * pb_y + g_0_zzzzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzz_0_zzzzzz_0[i] = g_0_zzzzzz_0_zzzzzz_0[i] * pb_y + g_0_zzzzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 980-1008 components of targeted buffer : SKSI

    auto g_0_zzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 980);

    auto g_0_zzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 981);

    auto g_0_zzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 982);

    auto g_0_zzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 983);

    auto g_0_zzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 984);

    auto g_0_zzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 985);

    auto g_0_zzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 986);

    auto g_0_zzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 987);

    auto g_0_zzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 988);

    auto g_0_zzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 989);

    auto g_0_zzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 990);

    auto g_0_zzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 991);

    auto g_0_zzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 992);

    auto g_0_zzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 993);

    auto g_0_zzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 994);

    auto g_0_zzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 995);

    auto g_0_zzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 996);

    auto g_0_zzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 997);

    auto g_0_zzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 998);

    auto g_0_zzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 999);

    auto g_0_zzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 1000);

    auto g_0_zzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 1001);

    auto g_0_zzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 1002);

    auto g_0_zzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 1003);

    auto g_0_zzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 1004);

    auto g_0_zzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 1005);

    auto g_0_zzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 1006);

    auto g_0_zzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 1007);

#pragma omp simd aligned(g_0_zzzzz_0_xxxxxx_0,       \
                             g_0_zzzzz_0_xxxxxx_1,   \
                             g_0_zzzzz_0_xxxxxy_0,   \
                             g_0_zzzzz_0_xxxxxy_1,   \
                             g_0_zzzzz_0_xxxxxz_0,   \
                             g_0_zzzzz_0_xxxxxz_1,   \
                             g_0_zzzzz_0_xxxxyy_0,   \
                             g_0_zzzzz_0_xxxxyy_1,   \
                             g_0_zzzzz_0_xxxxyz_0,   \
                             g_0_zzzzz_0_xxxxyz_1,   \
                             g_0_zzzzz_0_xxxxzz_0,   \
                             g_0_zzzzz_0_xxxxzz_1,   \
                             g_0_zzzzz_0_xxxyyy_0,   \
                             g_0_zzzzz_0_xxxyyy_1,   \
                             g_0_zzzzz_0_xxxyyz_0,   \
                             g_0_zzzzz_0_xxxyyz_1,   \
                             g_0_zzzzz_0_xxxyzz_0,   \
                             g_0_zzzzz_0_xxxyzz_1,   \
                             g_0_zzzzz_0_xxxzzz_0,   \
                             g_0_zzzzz_0_xxxzzz_1,   \
                             g_0_zzzzz_0_xxyyyy_0,   \
                             g_0_zzzzz_0_xxyyyy_1,   \
                             g_0_zzzzz_0_xxyyyz_0,   \
                             g_0_zzzzz_0_xxyyyz_1,   \
                             g_0_zzzzz_0_xxyyzz_0,   \
                             g_0_zzzzz_0_xxyyzz_1,   \
                             g_0_zzzzz_0_xxyzzz_0,   \
                             g_0_zzzzz_0_xxyzzz_1,   \
                             g_0_zzzzz_0_xxzzzz_0,   \
                             g_0_zzzzz_0_xxzzzz_1,   \
                             g_0_zzzzz_0_xyyyyy_0,   \
                             g_0_zzzzz_0_xyyyyy_1,   \
                             g_0_zzzzz_0_xyyyyz_0,   \
                             g_0_zzzzz_0_xyyyyz_1,   \
                             g_0_zzzzz_0_xyyyzz_0,   \
                             g_0_zzzzz_0_xyyyzz_1,   \
                             g_0_zzzzz_0_xyyzzz_0,   \
                             g_0_zzzzz_0_xyyzzz_1,   \
                             g_0_zzzzz_0_xyzzzz_0,   \
                             g_0_zzzzz_0_xyzzzz_1,   \
                             g_0_zzzzz_0_xzzzzz_0,   \
                             g_0_zzzzz_0_xzzzzz_1,   \
                             g_0_zzzzz_0_yyyyyy_0,   \
                             g_0_zzzzz_0_yyyyyy_1,   \
                             g_0_zzzzz_0_yyyyyz_0,   \
                             g_0_zzzzz_0_yyyyyz_1,   \
                             g_0_zzzzz_0_yyyyzz_0,   \
                             g_0_zzzzz_0_yyyyzz_1,   \
                             g_0_zzzzz_0_yyyzzz_0,   \
                             g_0_zzzzz_0_yyyzzz_1,   \
                             g_0_zzzzz_0_yyzzzz_0,   \
                             g_0_zzzzz_0_yyzzzz_1,   \
                             g_0_zzzzz_0_yzzzzz_0,   \
                             g_0_zzzzz_0_yzzzzz_1,   \
                             g_0_zzzzz_0_zzzzzz_0,   \
                             g_0_zzzzz_0_zzzzzz_1,   \
                             g_0_zzzzzz_0_xxxxx_1,   \
                             g_0_zzzzzz_0_xxxxxx_0,  \
                             g_0_zzzzzz_0_xxxxxx_1,  \
                             g_0_zzzzzz_0_xxxxxy_0,  \
                             g_0_zzzzzz_0_xxxxxy_1,  \
                             g_0_zzzzzz_0_xxxxxz_0,  \
                             g_0_zzzzzz_0_xxxxxz_1,  \
                             g_0_zzzzzz_0_xxxxy_1,   \
                             g_0_zzzzzz_0_xxxxyy_0,  \
                             g_0_zzzzzz_0_xxxxyy_1,  \
                             g_0_zzzzzz_0_xxxxyz_0,  \
                             g_0_zzzzzz_0_xxxxyz_1,  \
                             g_0_zzzzzz_0_xxxxz_1,   \
                             g_0_zzzzzz_0_xxxxzz_0,  \
                             g_0_zzzzzz_0_xxxxzz_1,  \
                             g_0_zzzzzz_0_xxxyy_1,   \
                             g_0_zzzzzz_0_xxxyyy_0,  \
                             g_0_zzzzzz_0_xxxyyy_1,  \
                             g_0_zzzzzz_0_xxxyyz_0,  \
                             g_0_zzzzzz_0_xxxyyz_1,  \
                             g_0_zzzzzz_0_xxxyz_1,   \
                             g_0_zzzzzz_0_xxxyzz_0,  \
                             g_0_zzzzzz_0_xxxyzz_1,  \
                             g_0_zzzzzz_0_xxxzz_1,   \
                             g_0_zzzzzz_0_xxxzzz_0,  \
                             g_0_zzzzzz_0_xxxzzz_1,  \
                             g_0_zzzzzz_0_xxyyy_1,   \
                             g_0_zzzzzz_0_xxyyyy_0,  \
                             g_0_zzzzzz_0_xxyyyy_1,  \
                             g_0_zzzzzz_0_xxyyyz_0,  \
                             g_0_zzzzzz_0_xxyyyz_1,  \
                             g_0_zzzzzz_0_xxyyz_1,   \
                             g_0_zzzzzz_0_xxyyzz_0,  \
                             g_0_zzzzzz_0_xxyyzz_1,  \
                             g_0_zzzzzz_0_xxyzz_1,   \
                             g_0_zzzzzz_0_xxyzzz_0,  \
                             g_0_zzzzzz_0_xxyzzz_1,  \
                             g_0_zzzzzz_0_xxzzz_1,   \
                             g_0_zzzzzz_0_xxzzzz_0,  \
                             g_0_zzzzzz_0_xxzzzz_1,  \
                             g_0_zzzzzz_0_xyyyy_1,   \
                             g_0_zzzzzz_0_xyyyyy_0,  \
                             g_0_zzzzzz_0_xyyyyy_1,  \
                             g_0_zzzzzz_0_xyyyyz_0,  \
                             g_0_zzzzzz_0_xyyyyz_1,  \
                             g_0_zzzzzz_0_xyyyz_1,   \
                             g_0_zzzzzz_0_xyyyzz_0,  \
                             g_0_zzzzzz_0_xyyyzz_1,  \
                             g_0_zzzzzz_0_xyyzz_1,   \
                             g_0_zzzzzz_0_xyyzzz_0,  \
                             g_0_zzzzzz_0_xyyzzz_1,  \
                             g_0_zzzzzz_0_xyzzz_1,   \
                             g_0_zzzzzz_0_xyzzzz_0,  \
                             g_0_zzzzzz_0_xyzzzz_1,  \
                             g_0_zzzzzz_0_xzzzz_1,   \
                             g_0_zzzzzz_0_xzzzzz_0,  \
                             g_0_zzzzzz_0_xzzzzz_1,  \
                             g_0_zzzzzz_0_yyyyy_1,   \
                             g_0_zzzzzz_0_yyyyyy_0,  \
                             g_0_zzzzzz_0_yyyyyy_1,  \
                             g_0_zzzzzz_0_yyyyyz_0,  \
                             g_0_zzzzzz_0_yyyyyz_1,  \
                             g_0_zzzzzz_0_yyyyz_1,   \
                             g_0_zzzzzz_0_yyyyzz_0,  \
                             g_0_zzzzzz_0_yyyyzz_1,  \
                             g_0_zzzzzz_0_yyyzz_1,   \
                             g_0_zzzzzz_0_yyyzzz_0,  \
                             g_0_zzzzzz_0_yyyzzz_1,  \
                             g_0_zzzzzz_0_yyzzz_1,   \
                             g_0_zzzzzz_0_yyzzzz_0,  \
                             g_0_zzzzzz_0_yyzzzz_1,  \
                             g_0_zzzzzz_0_yzzzz_1,   \
                             g_0_zzzzzz_0_yzzzzz_0,  \
                             g_0_zzzzzz_0_yzzzzz_1,  \
                             g_0_zzzzzz_0_zzzzz_1,   \
                             g_0_zzzzzz_0_zzzzzz_0,  \
                             g_0_zzzzzz_0_zzzzzz_1,  \
                             g_0_zzzzzzz_0_xxxxxx_0, \
                             g_0_zzzzzzz_0_xxxxxy_0, \
                             g_0_zzzzzzz_0_xxxxxz_0, \
                             g_0_zzzzzzz_0_xxxxyy_0, \
                             g_0_zzzzzzz_0_xxxxyz_0, \
                             g_0_zzzzzzz_0_xxxxzz_0, \
                             g_0_zzzzzzz_0_xxxyyy_0, \
                             g_0_zzzzzzz_0_xxxyyz_0, \
                             g_0_zzzzzzz_0_xxxyzz_0, \
                             g_0_zzzzzzz_0_xxxzzz_0, \
                             g_0_zzzzzzz_0_xxyyyy_0, \
                             g_0_zzzzzzz_0_xxyyyz_0, \
                             g_0_zzzzzzz_0_xxyyzz_0, \
                             g_0_zzzzzzz_0_xxyzzz_0, \
                             g_0_zzzzzzz_0_xxzzzz_0, \
                             g_0_zzzzzzz_0_xyyyyy_0, \
                             g_0_zzzzzzz_0_xyyyyz_0, \
                             g_0_zzzzzzz_0_xyyyzz_0, \
                             g_0_zzzzzzz_0_xyyzzz_0, \
                             g_0_zzzzzzz_0_xyzzzz_0, \
                             g_0_zzzzzzz_0_xzzzzz_0, \
                             g_0_zzzzzzz_0_yyyyyy_0, \
                             g_0_zzzzzzz_0_yyyyyz_0, \
                             g_0_zzzzzzz_0_yyyyzz_0, \
                             g_0_zzzzzzz_0_yyyzzz_0, \
                             g_0_zzzzzzz_0_yyzzzz_0, \
                             g_0_zzzzzzz_0_yzzzzz_0, \
                             g_0_zzzzzzz_0_zzzzzz_0, \
                             wp_z,                   \
                             c_exps,                 \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzz_0_xxxxxx_0[i] = 6.0 * g_0_zzzzz_0_xxxxxx_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxxx_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_xxxxxx_0[i] * pb_z + g_0_zzzzzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxxy_0[i] = 6.0 * g_0_zzzzz_0_xxxxxy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxxy_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_xxxxxy_0[i] * pb_z + g_0_zzzzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxxz_0[i] = 6.0 * g_0_zzzzz_0_xxxxxz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxxz_0[i] * pb_z + g_0_zzzzzz_0_xxxxxz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxyy_0[i] = 6.0 * g_0_zzzzz_0_xxxxyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxyy_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_xxxxyy_0[i] * pb_z + g_0_zzzzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxyz_0[i] = 6.0 * g_0_zzzzz_0_xxxxyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxyz_0[i] * pb_z + g_0_zzzzzz_0_xxxxyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxxzz_0[i] = 6.0 * g_0_zzzzz_0_xxxxzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_zzzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxxzz_0[i] * pb_z + g_0_zzzzzz_0_xxxxzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxyyy_0[i] = 6.0 * g_0_zzzzz_0_xxxyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxyyy_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_xxxyyy_0[i] * pb_z + g_0_zzzzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxyyz_0[i] = 6.0 * g_0_zzzzz_0_xxxyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyyz_0[i] * pb_z + g_0_zzzzzz_0_xxxyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxyzz_0[i] = 6.0 * g_0_zzzzz_0_xxxyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_zzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxyzz_0[i] * pb_z + g_0_zzzzzz_0_xxxyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxxzzz_0[i] = 6.0 * g_0_zzzzz_0_xxxzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_zzzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxxzzz_0[i] * pb_z + g_0_zzzzzz_0_xxxzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyyyy_0[i] = 6.0 * g_0_zzzzz_0_xxyyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyyyy_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_xxyyyy_0[i] * pb_z + g_0_zzzzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyyyz_0[i] = 6.0 * g_0_zzzzz_0_xxyyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyyz_0[i] * pb_z + g_0_zzzzzz_0_xxyyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyyzz_0[i] = 6.0 * g_0_zzzzz_0_xxyyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_zzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyyzz_0[i] * pb_z + g_0_zzzzzz_0_xxyyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxyzzz_0[i] = 6.0 * g_0_zzzzz_0_xxyzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_zzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxyzzz_0[i] * pb_z + g_0_zzzzzz_0_xxyzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xxzzzz_0[i] = 6.0 * g_0_zzzzz_0_xxzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_zzzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xxzzzz_0[i] * pb_z + g_0_zzzzzz_0_xxzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyyyy_0[i] = 6.0 * g_0_zzzzz_0_xyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyyyy_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_xyyyyy_0[i] * pb_z + g_0_zzzzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyyyz_0[i] = 6.0 * g_0_zzzzz_0_xyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyyz_0[i] * pb_z + g_0_zzzzzz_0_xyyyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyyzz_0[i] = 6.0 * g_0_zzzzz_0_xyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_zzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyyzz_0[i] * pb_z + g_0_zzzzzz_0_xyyyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyyzzz_0[i] = 6.0 * g_0_zzzzz_0_xyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_zzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyyzzz_0[i] * pb_z + g_0_zzzzzz_0_xyyzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xyzzzz_0[i] = 6.0 * g_0_zzzzz_0_xyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_zzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xyzzzz_0[i] * pb_z + g_0_zzzzzz_0_xyzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_xzzzzz_0[i] = 6.0 * g_0_zzzzz_0_xzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                    5.0 * g_0_zzzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_xzzzzz_0[i] * pb_z + g_0_zzzzzz_0_xzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyyyy_0[i] = 6.0 * g_0_zzzzz_0_yyyyyy_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyyyy_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_yyyyyy_0[i] * pb_z + g_0_zzzzzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyyyz_0[i] = 6.0 * g_0_zzzzz_0_yyyyyz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                    g_0_zzzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyyz_0[i] * pb_z + g_0_zzzzzz_0_yyyyyz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyyzz_0[i] = 6.0 * g_0_zzzzz_0_yyyyzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                    2.0 * g_0_zzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyyzz_0[i] * pb_z + g_0_zzzzzz_0_yyyyzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyyzzz_0[i] = 6.0 * g_0_zzzzz_0_yyyzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                    3.0 * g_0_zzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyyzzz_0[i] * pb_z + g_0_zzzzzz_0_yyyzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yyzzzz_0[i] = 6.0 * g_0_zzzzz_0_yyzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                    4.0 * g_0_zzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yyzzzz_0[i] * pb_z + g_0_zzzzzz_0_yyzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_yzzzzz_0[i] = 6.0 * g_0_zzzzz_0_yzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                    5.0 * g_0_zzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_yzzzzz_0[i] * pb_z + g_0_zzzzzz_0_yzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzz_0_zzzzzz_0[i] = 6.0 * g_0_zzzzz_0_zzzzzz_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                    6.0 * g_0_zzzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzzzzz_0_zzzzzz_0[i] * pb_z + g_0_zzzzzz_0_zzzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
