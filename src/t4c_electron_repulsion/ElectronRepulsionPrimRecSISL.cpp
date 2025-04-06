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

#include "ElectronRepulsionPrimRecSISL.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sisl(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sisl,
                                  size_t                idx_eri_0_sgsl,
                                  size_t                idx_eri_1_sgsl,
                                  size_t                idx_eri_1_shsk,
                                  size_t                idx_eri_0_shsl,
                                  size_t                idx_eri_1_shsl,
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

    /// Set up components of auxilary buffer : SGSL

    auto g_0_xxxx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl);

    auto g_0_xxxx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 1);

    auto g_0_xxxx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 2);

    auto g_0_xxxx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 3);

    auto g_0_xxxx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 4);

    auto g_0_xxxx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 5);

    auto g_0_xxxx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 6);

    auto g_0_xxxx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 7);

    auto g_0_xxxx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 8);

    auto g_0_xxxx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 9);

    auto g_0_xxxx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 10);

    auto g_0_xxxx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 11);

    auto g_0_xxxx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 12);

    auto g_0_xxxx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 13);

    auto g_0_xxxx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 14);

    auto g_0_xxxx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 15);

    auto g_0_xxxx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 16);

    auto g_0_xxxx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 17);

    auto g_0_xxxx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 18);

    auto g_0_xxxx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 19);

    auto g_0_xxxx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 20);

    auto g_0_xxxx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 21);

    auto g_0_xxxx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 22);

    auto g_0_xxxx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 23);

    auto g_0_xxxx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 24);

    auto g_0_xxxx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 25);

    auto g_0_xxxx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 26);

    auto g_0_xxxx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 27);

    auto g_0_xxxx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 28);

    auto g_0_xxxx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 29);

    auto g_0_xxxx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 30);

    auto g_0_xxxx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 31);

    auto g_0_xxxx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 32);

    auto g_0_xxxx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 33);

    auto g_0_xxxx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 34);

    auto g_0_xxxx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 35);

    auto g_0_xxxx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 36);

    auto g_0_xxxx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 37);

    auto g_0_xxxx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 38);

    auto g_0_xxxx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 39);

    auto g_0_xxxx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 40);

    auto g_0_xxxx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 41);

    auto g_0_xxxx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 42);

    auto g_0_xxxx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 43);

    auto g_0_xxxx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 44);

    auto g_0_xxxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 45);

    auto g_0_xxxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 47);

    auto g_0_xxxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 50);

    auto g_0_xxxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 54);

    auto g_0_xxxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 59);

    auto g_0_xxxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 65);

    auto g_0_xxxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 72);

    auto g_0_xxxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 80);

    auto g_0_xxxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 90);

    auto g_0_xxxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 91);

    auto g_0_xxxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 93);

    auto g_0_xxxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 96);

    auto g_0_xxxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 100);

    auto g_0_xxxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 105);

    auto g_0_xxxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 111);

    auto g_0_xxxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 118);

    auto g_0_xxyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 135);

    auto g_0_xxyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 136);

    auto g_0_xxyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 137);

    auto g_0_xxyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 138);

    auto g_0_xxyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 139);

    auto g_0_xxyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 140);

    auto g_0_xxyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 141);

    auto g_0_xxyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 142);

    auto g_0_xxyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 143);

    auto g_0_xxyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 144);

    auto g_0_xxyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 145);

    auto g_0_xxyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 146);

    auto g_0_xxyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 147);

    auto g_0_xxyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 148);

    auto g_0_xxyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 149);

    auto g_0_xxyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 150);

    auto g_0_xxyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 151);

    auto g_0_xxyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 152);

    auto g_0_xxyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 153);

    auto g_0_xxyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 154);

    auto g_0_xxyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 155);

    auto g_0_xxyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 156);

    auto g_0_xxyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 157);

    auto g_0_xxyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 158);

    auto g_0_xxyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 159);

    auto g_0_xxyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 160);

    auto g_0_xxyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 161);

    auto g_0_xxyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 162);

    auto g_0_xxyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 163);

    auto g_0_xxyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 164);

    auto g_0_xxyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 165);

    auto g_0_xxyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 166);

    auto g_0_xxyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 167);

    auto g_0_xxyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 168);

    auto g_0_xxyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 169);

    auto g_0_xxyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 170);

    auto g_0_xxyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 171);

    auto g_0_xxyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 172);

    auto g_0_xxyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 173);

    auto g_0_xxyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 174);

    auto g_0_xxyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 175);

    auto g_0_xxyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 176);

    auto g_0_xxyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 177);

    auto g_0_xxyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 178);

    auto g_0_xxyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 179);

    auto g_0_xxzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 225);

    auto g_0_xxzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 226);

    auto g_0_xxzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 227);

    auto g_0_xxzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 228);

    auto g_0_xxzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 229);

    auto g_0_xxzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 230);

    auto g_0_xxzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 231);

    auto g_0_xxzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 232);

    auto g_0_xxzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 233);

    auto g_0_xxzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 234);

    auto g_0_xxzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 235);

    auto g_0_xxzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 236);

    auto g_0_xxzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 237);

    auto g_0_xxzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 238);

    auto g_0_xxzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 239);

    auto g_0_xxzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 240);

    auto g_0_xxzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 241);

    auto g_0_xxzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 242);

    auto g_0_xxzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 243);

    auto g_0_xxzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 244);

    auto g_0_xxzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 245);

    auto g_0_xxzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 246);

    auto g_0_xxzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 247);

    auto g_0_xxzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 248);

    auto g_0_xxzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 249);

    auto g_0_xxzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 250);

    auto g_0_xxzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 251);

    auto g_0_xxzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 252);

    auto g_0_xxzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 253);

    auto g_0_xxzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 254);

    auto g_0_xxzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 255);

    auto g_0_xxzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 256);

    auto g_0_xxzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 257);

    auto g_0_xxzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 258);

    auto g_0_xxzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 259);

    auto g_0_xxzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 260);

    auto g_0_xxzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 261);

    auto g_0_xxzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 262);

    auto g_0_xxzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 263);

    auto g_0_xxzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 264);

    auto g_0_xxzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 265);

    auto g_0_xxzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 266);

    auto g_0_xxzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 267);

    auto g_0_xxzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 268);

    auto g_0_xxzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 269);

    auto g_0_xyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 271);

    auto g_0_xyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 273);

    auto g_0_xyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 274);

    auto g_0_xyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 276);

    auto g_0_xyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 277);

    auto g_0_xyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 278);

    auto g_0_xyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 280);

    auto g_0_xyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 281);

    auto g_0_xyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 282);

    auto g_0_xyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 283);

    auto g_0_xyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 285);

    auto g_0_xyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 286);

    auto g_0_xyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 287);

    auto g_0_xyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 288);

    auto g_0_xyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 289);

    auto g_0_xyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 291);

    auto g_0_xyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 292);

    auto g_0_xyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 293);

    auto g_0_xyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 294);

    auto g_0_xyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 295);

    auto g_0_xyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 296);

    auto g_0_xyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 298);

    auto g_0_xyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 299);

    auto g_0_xyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 300);

    auto g_0_xyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 301);

    auto g_0_xyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 302);

    auto g_0_xyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 303);

    auto g_0_xyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 304);

    auto g_0_xyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 306);

    auto g_0_xyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 307);

    auto g_0_xyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 308);

    auto g_0_xyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 309);

    auto g_0_xyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 310);

    auto g_0_xyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 311);

    auto g_0_xyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 312);

    auto g_0_xyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 313);

    auto g_0_xyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 314);

    auto g_0_xzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 407);

    auto g_0_xzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 409);

    auto g_0_xzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 410);

    auto g_0_xzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 412);

    auto g_0_xzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 413);

    auto g_0_xzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 414);

    auto g_0_xzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 416);

    auto g_0_xzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 417);

    auto g_0_xzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 418);

    auto g_0_xzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 419);

    auto g_0_xzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 421);

    auto g_0_xzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 422);

    auto g_0_xzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 423);

    auto g_0_xzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 424);

    auto g_0_xzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 425);

    auto g_0_xzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 427);

    auto g_0_xzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 428);

    auto g_0_xzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 429);

    auto g_0_xzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 430);

    auto g_0_xzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 431);

    auto g_0_xzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 432);

    auto g_0_xzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 434);

    auto g_0_xzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 435);

    auto g_0_xzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 436);

    auto g_0_xzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 437);

    auto g_0_xzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 438);

    auto g_0_xzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 439);

    auto g_0_xzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 440);

    auto g_0_xzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 441);

    auto g_0_xzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 442);

    auto g_0_xzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 443);

    auto g_0_xzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 444);

    auto g_0_xzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 445);

    auto g_0_xzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 446);

    auto g_0_xzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 447);

    auto g_0_xzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 448);

    auto g_0_xzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 449);

    auto g_0_yyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 450);

    auto g_0_yyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 451);

    auto g_0_yyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 452);

    auto g_0_yyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 453);

    auto g_0_yyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 454);

    auto g_0_yyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 455);

    auto g_0_yyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 456);

    auto g_0_yyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 457);

    auto g_0_yyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 458);

    auto g_0_yyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 459);

    auto g_0_yyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 460);

    auto g_0_yyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 461);

    auto g_0_yyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 462);

    auto g_0_yyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 463);

    auto g_0_yyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 464);

    auto g_0_yyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 465);

    auto g_0_yyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 466);

    auto g_0_yyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 467);

    auto g_0_yyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 468);

    auto g_0_yyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 469);

    auto g_0_yyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 470);

    auto g_0_yyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 471);

    auto g_0_yyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 472);

    auto g_0_yyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 473);

    auto g_0_yyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 474);

    auto g_0_yyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 475);

    auto g_0_yyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 476);

    auto g_0_yyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 477);

    auto g_0_yyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 478);

    auto g_0_yyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 479);

    auto g_0_yyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 480);

    auto g_0_yyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 481);

    auto g_0_yyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 482);

    auto g_0_yyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 483);

    auto g_0_yyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 484);

    auto g_0_yyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 485);

    auto g_0_yyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 486);

    auto g_0_yyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 487);

    auto g_0_yyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 488);

    auto g_0_yyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 489);

    auto g_0_yyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 490);

    auto g_0_yyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 491);

    auto g_0_yyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 492);

    auto g_0_yyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 493);

    auto g_0_yyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 494);

    auto g_0_yyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 496);

    auto g_0_yyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 498);

    auto g_0_yyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 501);

    auto g_0_yyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 505);

    auto g_0_yyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 510);

    auto g_0_yyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 516);

    auto g_0_yyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 523);

    auto g_0_yyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 531);

    auto g_0_yyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 540);

    auto g_0_yyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 541);

    auto g_0_yyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 542);

    auto g_0_yyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 543);

    auto g_0_yyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 544);

    auto g_0_yyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 545);

    auto g_0_yyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 546);

    auto g_0_yyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 547);

    auto g_0_yyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 548);

    auto g_0_yyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 549);

    auto g_0_yyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 550);

    auto g_0_yyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 551);

    auto g_0_yyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 552);

    auto g_0_yyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 553);

    auto g_0_yyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 554);

    auto g_0_yyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 555);

    auto g_0_yyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 556);

    auto g_0_yyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 557);

    auto g_0_yyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 558);

    auto g_0_yyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 559);

    auto g_0_yyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 560);

    auto g_0_yyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 561);

    auto g_0_yyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 562);

    auto g_0_yyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 563);

    auto g_0_yyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 564);

    auto g_0_yyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 565);

    auto g_0_yyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 566);

    auto g_0_yyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 567);

    auto g_0_yyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 568);

    auto g_0_yyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 569);

    auto g_0_yyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 570);

    auto g_0_yyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 571);

    auto g_0_yyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 572);

    auto g_0_yyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 573);

    auto g_0_yyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 574);

    auto g_0_yyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 575);

    auto g_0_yyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 576);

    auto g_0_yyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 577);

    auto g_0_yyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 578);

    auto g_0_yyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 579);

    auto g_0_yyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 580);

    auto g_0_yyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 581);

    auto g_0_yyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 582);

    auto g_0_yyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 583);

    auto g_0_yyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 584);

    auto g_0_yzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 585);

    auto g_0_yzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 587);

    auto g_0_yzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 589);

    auto g_0_yzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 590);

    auto g_0_yzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 592);

    auto g_0_yzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 593);

    auto g_0_yzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 594);

    auto g_0_yzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 596);

    auto g_0_yzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 597);

    auto g_0_yzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 598);

    auto g_0_yzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 599);

    auto g_0_yzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 601);

    auto g_0_yzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 602);

    auto g_0_yzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 603);

    auto g_0_yzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 604);

    auto g_0_yzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 605);

    auto g_0_yzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 607);

    auto g_0_yzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 608);

    auto g_0_yzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 609);

    auto g_0_yzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 610);

    auto g_0_yzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 611);

    auto g_0_yzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 612);

    auto g_0_yzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 614);

    auto g_0_yzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 615);

    auto g_0_yzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 616);

    auto g_0_yzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 617);

    auto g_0_yzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 618);

    auto g_0_yzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 619);

    auto g_0_yzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 620);

    auto g_0_yzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 622);

    auto g_0_yzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 623);

    auto g_0_yzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 624);

    auto g_0_yzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 625);

    auto g_0_yzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 626);

    auto g_0_yzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 627);

    auto g_0_yzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 628);

    auto g_0_yzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 629);

    auto g_0_zzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sgsl + 630);

    auto g_0_zzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sgsl + 631);

    auto g_0_zzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sgsl + 632);

    auto g_0_zzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sgsl + 633);

    auto g_0_zzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sgsl + 634);

    auto g_0_zzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sgsl + 635);

    auto g_0_zzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sgsl + 636);

    auto g_0_zzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sgsl + 637);

    auto g_0_zzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sgsl + 638);

    auto g_0_zzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sgsl + 639);

    auto g_0_zzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 640);

    auto g_0_zzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 641);

    auto g_0_zzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 642);

    auto g_0_zzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 643);

    auto g_0_zzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 644);

    auto g_0_zzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 645);

    auto g_0_zzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 646);

    auto g_0_zzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 647);

    auto g_0_zzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 648);

    auto g_0_zzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 649);

    auto g_0_zzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 650);

    auto g_0_zzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 651);

    auto g_0_zzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 652);

    auto g_0_zzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 653);

    auto g_0_zzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 654);

    auto g_0_zzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 655);

    auto g_0_zzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 656);

    auto g_0_zzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 657);

    auto g_0_zzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 658);

    auto g_0_zzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 659);

    auto g_0_zzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 660);

    auto g_0_zzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 661);

    auto g_0_zzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 662);

    auto g_0_zzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 663);

    auto g_0_zzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 664);

    auto g_0_zzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 665);

    auto g_0_zzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sgsl + 666);

    auto g_0_zzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sgsl + 667);

    auto g_0_zzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sgsl + 668);

    auto g_0_zzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sgsl + 669);

    auto g_0_zzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 670);

    auto g_0_zzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 671);

    auto g_0_zzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 672);

    auto g_0_zzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 673);

    auto g_0_zzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sgsl + 674);

    /// Set up components of auxilary buffer : SGSL

    auto g_0_xxxx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sgsl);

    auto g_0_xxxx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sgsl + 1);

    auto g_0_xxxx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sgsl + 2);

    auto g_0_xxxx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sgsl + 3);

    auto g_0_xxxx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sgsl + 4);

    auto g_0_xxxx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sgsl + 5);

    auto g_0_xxxx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sgsl + 6);

    auto g_0_xxxx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sgsl + 7);

    auto g_0_xxxx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sgsl + 8);

    auto g_0_xxxx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sgsl + 9);

    auto g_0_xxxx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 10);

    auto g_0_xxxx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 11);

    auto g_0_xxxx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 12);

    auto g_0_xxxx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 13);

    auto g_0_xxxx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 14);

    auto g_0_xxxx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 15);

    auto g_0_xxxx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 16);

    auto g_0_xxxx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 17);

    auto g_0_xxxx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 18);

    auto g_0_xxxx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 19);

    auto g_0_xxxx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 20);

    auto g_0_xxxx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 21);

    auto g_0_xxxx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 22);

    auto g_0_xxxx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 23);

    auto g_0_xxxx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 24);

    auto g_0_xxxx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 25);

    auto g_0_xxxx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 26);

    auto g_0_xxxx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 27);

    auto g_0_xxxx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 28);

    auto g_0_xxxx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 29);

    auto g_0_xxxx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 30);

    auto g_0_xxxx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 31);

    auto g_0_xxxx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 32);

    auto g_0_xxxx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 33);

    auto g_0_xxxx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 34);

    auto g_0_xxxx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 35);

    auto g_0_xxxx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 36);

    auto g_0_xxxx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 37);

    auto g_0_xxxx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 38);

    auto g_0_xxxx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 39);

    auto g_0_xxxx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 40);

    auto g_0_xxxx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 41);

    auto g_0_xxxx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 42);

    auto g_0_xxxx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 43);

    auto g_0_xxxx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 44);

    auto g_0_xxxy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sgsl + 45);

    auto g_0_xxxy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sgsl + 47);

    auto g_0_xxxy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sgsl + 50);

    auto g_0_xxxy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sgsl + 54);

    auto g_0_xxxy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 59);

    auto g_0_xxxy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 65);

    auto g_0_xxxy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 72);

    auto g_0_xxxy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 80);

    auto g_0_xxxz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sgsl + 90);

    auto g_0_xxxz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sgsl + 91);

    auto g_0_xxxz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sgsl + 93);

    auto g_0_xxxz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sgsl + 96);

    auto g_0_xxxz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 100);

    auto g_0_xxxz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 105);

    auto g_0_xxxz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 111);

    auto g_0_xxxz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 118);

    auto g_0_xxyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sgsl + 135);

    auto g_0_xxyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sgsl + 136);

    auto g_0_xxyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sgsl + 137);

    auto g_0_xxyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sgsl + 138);

    auto g_0_xxyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sgsl + 139);

    auto g_0_xxyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sgsl + 140);

    auto g_0_xxyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sgsl + 141);

    auto g_0_xxyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sgsl + 142);

    auto g_0_xxyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sgsl + 143);

    auto g_0_xxyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sgsl + 144);

    auto g_0_xxyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 145);

    auto g_0_xxyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 146);

    auto g_0_xxyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 147);

    auto g_0_xxyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 148);

    auto g_0_xxyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 149);

    auto g_0_xxyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 150);

    auto g_0_xxyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 151);

    auto g_0_xxyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 152);

    auto g_0_xxyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 153);

    auto g_0_xxyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 154);

    auto g_0_xxyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 155);

    auto g_0_xxyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 156);

    auto g_0_xxyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 157);

    auto g_0_xxyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 158);

    auto g_0_xxyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 159);

    auto g_0_xxyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 160);

    auto g_0_xxyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 161);

    auto g_0_xxyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 162);

    auto g_0_xxyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 163);

    auto g_0_xxyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 164);

    auto g_0_xxyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 165);

    auto g_0_xxyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 166);

    auto g_0_xxyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 167);

    auto g_0_xxyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 168);

    auto g_0_xxyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 169);

    auto g_0_xxyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 170);

    auto g_0_xxyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 171);

    auto g_0_xxyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 172);

    auto g_0_xxyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 173);

    auto g_0_xxyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 174);

    auto g_0_xxyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 175);

    auto g_0_xxyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 176);

    auto g_0_xxyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 177);

    auto g_0_xxyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 178);

    auto g_0_xxyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 179);

    auto g_0_xxzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sgsl + 225);

    auto g_0_xxzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sgsl + 226);

    auto g_0_xxzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sgsl + 227);

    auto g_0_xxzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sgsl + 228);

    auto g_0_xxzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sgsl + 229);

    auto g_0_xxzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sgsl + 230);

    auto g_0_xxzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sgsl + 231);

    auto g_0_xxzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sgsl + 232);

    auto g_0_xxzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sgsl + 233);

    auto g_0_xxzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sgsl + 234);

    auto g_0_xxzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 235);

    auto g_0_xxzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 236);

    auto g_0_xxzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 237);

    auto g_0_xxzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 238);

    auto g_0_xxzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 239);

    auto g_0_xxzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 240);

    auto g_0_xxzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 241);

    auto g_0_xxzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 242);

    auto g_0_xxzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 243);

    auto g_0_xxzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 244);

    auto g_0_xxzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 245);

    auto g_0_xxzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 246);

    auto g_0_xxzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 247);

    auto g_0_xxzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 248);

    auto g_0_xxzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 249);

    auto g_0_xxzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 250);

    auto g_0_xxzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 251);

    auto g_0_xxzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 252);

    auto g_0_xxzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 253);

    auto g_0_xxzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 254);

    auto g_0_xxzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 255);

    auto g_0_xxzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 256);

    auto g_0_xxzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 257);

    auto g_0_xxzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 258);

    auto g_0_xxzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 259);

    auto g_0_xxzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 260);

    auto g_0_xxzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 261);

    auto g_0_xxzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 262);

    auto g_0_xxzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 263);

    auto g_0_xxzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 264);

    auto g_0_xxzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 265);

    auto g_0_xxzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 266);

    auto g_0_xxzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 267);

    auto g_0_xxzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 268);

    auto g_0_xxzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 269);

    auto g_0_xyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sgsl + 271);

    auto g_0_xyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sgsl + 273);

    auto g_0_xyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sgsl + 274);

    auto g_0_xyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sgsl + 276);

    auto g_0_xyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sgsl + 277);

    auto g_0_xyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sgsl + 278);

    auto g_0_xyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 280);

    auto g_0_xyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 281);

    auto g_0_xyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 282);

    auto g_0_xyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 283);

    auto g_0_xyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 285);

    auto g_0_xyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 286);

    auto g_0_xyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 287);

    auto g_0_xyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 288);

    auto g_0_xyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 289);

    auto g_0_xyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 291);

    auto g_0_xyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 292);

    auto g_0_xyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 293);

    auto g_0_xyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 294);

    auto g_0_xyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 295);

    auto g_0_xyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 296);

    auto g_0_xyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 298);

    auto g_0_xyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 299);

    auto g_0_xyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 300);

    auto g_0_xyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 301);

    auto g_0_xyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 302);

    auto g_0_xyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 303);

    auto g_0_xyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 304);

    auto g_0_xyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 306);

    auto g_0_xyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 307);

    auto g_0_xyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 308);

    auto g_0_xyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 309);

    auto g_0_xyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 310);

    auto g_0_xyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 311);

    auto g_0_xyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 312);

    auto g_0_xyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 313);

    auto g_0_xyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 314);

    auto g_0_xzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sgsl + 407);

    auto g_0_xzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sgsl + 409);

    auto g_0_xzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sgsl + 410);

    auto g_0_xzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sgsl + 412);

    auto g_0_xzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sgsl + 413);

    auto g_0_xzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sgsl + 414);

    auto g_0_xzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 416);

    auto g_0_xzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 417);

    auto g_0_xzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 418);

    auto g_0_xzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 419);

    auto g_0_xzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 421);

    auto g_0_xzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 422);

    auto g_0_xzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 423);

    auto g_0_xzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 424);

    auto g_0_xzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 425);

    auto g_0_xzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 427);

    auto g_0_xzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 428);

    auto g_0_xzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 429);

    auto g_0_xzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 430);

    auto g_0_xzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 431);

    auto g_0_xzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 432);

    auto g_0_xzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 434);

    auto g_0_xzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 435);

    auto g_0_xzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 436);

    auto g_0_xzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 437);

    auto g_0_xzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 438);

    auto g_0_xzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 439);

    auto g_0_xzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 440);

    auto g_0_xzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 441);

    auto g_0_xzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 442);

    auto g_0_xzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 443);

    auto g_0_xzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 444);

    auto g_0_xzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 445);

    auto g_0_xzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 446);

    auto g_0_xzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 447);

    auto g_0_xzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 448);

    auto g_0_xzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 449);

    auto g_0_yyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sgsl + 450);

    auto g_0_yyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sgsl + 451);

    auto g_0_yyyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sgsl + 452);

    auto g_0_yyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sgsl + 453);

    auto g_0_yyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sgsl + 454);

    auto g_0_yyyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sgsl + 455);

    auto g_0_yyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sgsl + 456);

    auto g_0_yyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sgsl + 457);

    auto g_0_yyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sgsl + 458);

    auto g_0_yyyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sgsl + 459);

    auto g_0_yyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 460);

    auto g_0_yyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 461);

    auto g_0_yyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 462);

    auto g_0_yyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 463);

    auto g_0_yyyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 464);

    auto g_0_yyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 465);

    auto g_0_yyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 466);

    auto g_0_yyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 467);

    auto g_0_yyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 468);

    auto g_0_yyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 469);

    auto g_0_yyyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 470);

    auto g_0_yyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 471);

    auto g_0_yyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 472);

    auto g_0_yyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 473);

    auto g_0_yyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 474);

    auto g_0_yyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 475);

    auto g_0_yyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 476);

    auto g_0_yyyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 477);

    auto g_0_yyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 478);

    auto g_0_yyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 479);

    auto g_0_yyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 480);

    auto g_0_yyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 481);

    auto g_0_yyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 482);

    auto g_0_yyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 483);

    auto g_0_yyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 484);

    auto g_0_yyyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 485);

    auto g_0_yyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 486);

    auto g_0_yyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 487);

    auto g_0_yyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 488);

    auto g_0_yyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 489);

    auto g_0_yyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 490);

    auto g_0_yyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 491);

    auto g_0_yyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 492);

    auto g_0_yyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 493);

    auto g_0_yyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 494);

    auto g_0_yyyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sgsl + 496);

    auto g_0_yyyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sgsl + 498);

    auto g_0_yyyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sgsl + 501);

    auto g_0_yyyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 505);

    auto g_0_yyyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 510);

    auto g_0_yyyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 516);

    auto g_0_yyyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 523);

    auto g_0_yyyz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 531);

    auto g_0_yyzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sgsl + 540);

    auto g_0_yyzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sgsl + 541);

    auto g_0_yyzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sgsl + 542);

    auto g_0_yyzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sgsl + 543);

    auto g_0_yyzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sgsl + 544);

    auto g_0_yyzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sgsl + 545);

    auto g_0_yyzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sgsl + 546);

    auto g_0_yyzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sgsl + 547);

    auto g_0_yyzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sgsl + 548);

    auto g_0_yyzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sgsl + 549);

    auto g_0_yyzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 550);

    auto g_0_yyzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 551);

    auto g_0_yyzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 552);

    auto g_0_yyzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 553);

    auto g_0_yyzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 554);

    auto g_0_yyzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 555);

    auto g_0_yyzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 556);

    auto g_0_yyzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 557);

    auto g_0_yyzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 558);

    auto g_0_yyzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 559);

    auto g_0_yyzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 560);

    auto g_0_yyzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 561);

    auto g_0_yyzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 562);

    auto g_0_yyzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 563);

    auto g_0_yyzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 564);

    auto g_0_yyzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 565);

    auto g_0_yyzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 566);

    auto g_0_yyzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 567);

    auto g_0_yyzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 568);

    auto g_0_yyzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 569);

    auto g_0_yyzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 570);

    auto g_0_yyzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 571);

    auto g_0_yyzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 572);

    auto g_0_yyzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 573);

    auto g_0_yyzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 574);

    auto g_0_yyzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 575);

    auto g_0_yyzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 576);

    auto g_0_yyzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 577);

    auto g_0_yyzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 578);

    auto g_0_yyzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 579);

    auto g_0_yyzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 580);

    auto g_0_yyzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 581);

    auto g_0_yyzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 582);

    auto g_0_yyzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 583);

    auto g_0_yyzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 584);

    auto g_0_yzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sgsl + 585);

    auto g_0_yzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sgsl + 587);

    auto g_0_yzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sgsl + 589);

    auto g_0_yzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sgsl + 590);

    auto g_0_yzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sgsl + 592);

    auto g_0_yzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sgsl + 593);

    auto g_0_yzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sgsl + 594);

    auto g_0_yzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 596);

    auto g_0_yzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 597);

    auto g_0_yzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 598);

    auto g_0_yzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 599);

    auto g_0_yzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 601);

    auto g_0_yzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 602);

    auto g_0_yzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 603);

    auto g_0_yzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 604);

    auto g_0_yzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 605);

    auto g_0_yzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 607);

    auto g_0_yzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 608);

    auto g_0_yzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 609);

    auto g_0_yzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 610);

    auto g_0_yzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 611);

    auto g_0_yzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 612);

    auto g_0_yzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 614);

    auto g_0_yzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 615);

    auto g_0_yzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 616);

    auto g_0_yzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 617);

    auto g_0_yzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 618);

    auto g_0_yzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 619);

    auto g_0_yzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 620);

    auto g_0_yzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 622);

    auto g_0_yzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 623);

    auto g_0_yzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 624);

    auto g_0_yzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 625);

    auto g_0_yzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 626);

    auto g_0_yzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 627);

    auto g_0_yzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 628);

    auto g_0_yzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 629);

    auto g_0_zzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sgsl + 630);

    auto g_0_zzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sgsl + 631);

    auto g_0_zzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sgsl + 632);

    auto g_0_zzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sgsl + 633);

    auto g_0_zzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sgsl + 634);

    auto g_0_zzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sgsl + 635);

    auto g_0_zzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sgsl + 636);

    auto g_0_zzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sgsl + 637);

    auto g_0_zzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sgsl + 638);

    auto g_0_zzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sgsl + 639);

    auto g_0_zzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 640);

    auto g_0_zzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 641);

    auto g_0_zzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 642);

    auto g_0_zzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 643);

    auto g_0_zzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 644);

    auto g_0_zzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 645);

    auto g_0_zzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 646);

    auto g_0_zzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 647);

    auto g_0_zzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 648);

    auto g_0_zzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 649);

    auto g_0_zzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 650);

    auto g_0_zzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 651);

    auto g_0_zzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 652);

    auto g_0_zzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 653);

    auto g_0_zzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 654);

    auto g_0_zzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 655);

    auto g_0_zzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 656);

    auto g_0_zzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 657);

    auto g_0_zzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 658);

    auto g_0_zzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 659);

    auto g_0_zzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 660);

    auto g_0_zzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 661);

    auto g_0_zzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 662);

    auto g_0_zzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 663);

    auto g_0_zzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 664);

    auto g_0_zzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 665);

    auto g_0_zzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sgsl + 666);

    auto g_0_zzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sgsl + 667);

    auto g_0_zzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sgsl + 668);

    auto g_0_zzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sgsl + 669);

    auto g_0_zzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 670);

    auto g_0_zzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 671);

    auto g_0_zzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 672);

    auto g_0_zzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 673);

    auto g_0_zzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sgsl + 674);

    /// Set up components of auxilary buffer : SHSK

    auto g_0_xxxxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk);

    auto g_0_xxxxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 1);

    auto g_0_xxxxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 2);

    auto g_0_xxxxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 3);

    auto g_0_xxxxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 4);

    auto g_0_xxxxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 5);

    auto g_0_xxxxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 6);

    auto g_0_xxxxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 7);

    auto g_0_xxxxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 8);

    auto g_0_xxxxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 9);

    auto g_0_xxxxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 10);

    auto g_0_xxxxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 11);

    auto g_0_xxxxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 12);

    auto g_0_xxxxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 13);

    auto g_0_xxxxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 14);

    auto g_0_xxxxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 15);

    auto g_0_xxxxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 16);

    auto g_0_xxxxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 17);

    auto g_0_xxxxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 18);

    auto g_0_xxxxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 19);

    auto g_0_xxxxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 20);

    auto g_0_xxxxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 21);

    auto g_0_xxxxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 22);

    auto g_0_xxxxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 23);

    auto g_0_xxxxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 24);

    auto g_0_xxxxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 25);

    auto g_0_xxxxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 26);

    auto g_0_xxxxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 27);

    auto g_0_xxxxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 28);

    auto g_0_xxxxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 29);

    auto g_0_xxxxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 30);

    auto g_0_xxxxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 31);

    auto g_0_xxxxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 32);

    auto g_0_xxxxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 33);

    auto g_0_xxxxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 34);

    auto g_0_xxxxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 35);

    auto g_0_xxxxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 74);

    auto g_0_xxxxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 76);

    auto g_0_xxxxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 77);

    auto g_0_xxxxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 79);

    auto g_0_xxxxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 80);

    auto g_0_xxxxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 81);

    auto g_0_xxxxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 83);

    auto g_0_xxxxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 84);

    auto g_0_xxxxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 85);

    auto g_0_xxxxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 86);

    auto g_0_xxxxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 88);

    auto g_0_xxxxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 89);

    auto g_0_xxxxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 90);

    auto g_0_xxxxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 91);

    auto g_0_xxxxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 92);

    auto g_0_xxxxz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 94);

    auto g_0_xxxxz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 95);

    auto g_0_xxxxz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 96);

    auto g_0_xxxxz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 97);

    auto g_0_xxxxz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 98);

    auto g_0_xxxxz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 99);

    auto g_0_xxxxz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 101);

    auto g_0_xxxxz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 102);

    auto g_0_xxxxz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 103);

    auto g_0_xxxxz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 104);

    auto g_0_xxxxz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 105);

    auto g_0_xxxxz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 106);

    auto g_0_xxxxz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 107);

    auto g_0_xxxyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 108);

    auto g_0_xxxyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 109);

    auto g_0_xxxyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 110);

    auto g_0_xxxyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 111);

    auto g_0_xxxyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 112);

    auto g_0_xxxyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 113);

    auto g_0_xxxyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 114);

    auto g_0_xxxyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 115);

    auto g_0_xxxyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 116);

    auto g_0_xxxyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 117);

    auto g_0_xxxyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 118);

    auto g_0_xxxyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 119);

    auto g_0_xxxyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 120);

    auto g_0_xxxyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 121);

    auto g_0_xxxyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 122);

    auto g_0_xxxyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 123);

    auto g_0_xxxyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 124);

    auto g_0_xxxyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 125);

    auto g_0_xxxyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 126);

    auto g_0_xxxyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 127);

    auto g_0_xxxyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 128);

    auto g_0_xxxyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 129);

    auto g_0_xxxyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 130);

    auto g_0_xxxyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 131);

    auto g_0_xxxyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 132);

    auto g_0_xxxyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 133);

    auto g_0_xxxyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 134);

    auto g_0_xxxyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 135);

    auto g_0_xxxyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 136);

    auto g_0_xxxyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 137);

    auto g_0_xxxyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 138);

    auto g_0_xxxyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 139);

    auto g_0_xxxyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 140);

    auto g_0_xxxyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 141);

    auto g_0_xxxyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 142);

    auto g_0_xxxyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 143);

    auto g_0_xxxzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 180);

    auto g_0_xxxzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 181);

    auto g_0_xxxzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 182);

    auto g_0_xxxzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 183);

    auto g_0_xxxzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 184);

    auto g_0_xxxzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 185);

    auto g_0_xxxzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 186);

    auto g_0_xxxzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 187);

    auto g_0_xxxzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 188);

    auto g_0_xxxzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 189);

    auto g_0_xxxzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 190);

    auto g_0_xxxzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 191);

    auto g_0_xxxzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 192);

    auto g_0_xxxzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 193);

    auto g_0_xxxzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 194);

    auto g_0_xxxzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 195);

    auto g_0_xxxzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 196);

    auto g_0_xxxzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 197);

    auto g_0_xxxzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 198);

    auto g_0_xxxzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 199);

    auto g_0_xxxzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 200);

    auto g_0_xxxzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 201);

    auto g_0_xxxzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 202);

    auto g_0_xxxzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 203);

    auto g_0_xxxzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 204);

    auto g_0_xxxzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 205);

    auto g_0_xxxzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 206);

    auto g_0_xxxzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 207);

    auto g_0_xxxzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 208);

    auto g_0_xxxzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 209);

    auto g_0_xxxzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 210);

    auto g_0_xxxzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 211);

    auto g_0_xxxzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 212);

    auto g_0_xxxzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 213);

    auto g_0_xxxzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 214);

    auto g_0_xxxzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 215);

    auto g_0_xxyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 216);

    auto g_0_xxyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 217);

    auto g_0_xxyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 218);

    auto g_0_xxyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 219);

    auto g_0_xxyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 220);

    auto g_0_xxyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 221);

    auto g_0_xxyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 222);

    auto g_0_xxyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 223);

    auto g_0_xxyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 224);

    auto g_0_xxyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 225);

    auto g_0_xxyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 226);

    auto g_0_xxyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 227);

    auto g_0_xxyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 228);

    auto g_0_xxyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 229);

    auto g_0_xxyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 230);

    auto g_0_xxyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 231);

    auto g_0_xxyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 232);

    auto g_0_xxyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 233);

    auto g_0_xxyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 234);

    auto g_0_xxyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 235);

    auto g_0_xxyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 236);

    auto g_0_xxyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 237);

    auto g_0_xxyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 238);

    auto g_0_xxyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 239);

    auto g_0_xxyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 240);

    auto g_0_xxyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 241);

    auto g_0_xxyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 242);

    auto g_0_xxyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 243);

    auto g_0_xxyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 244);

    auto g_0_xxyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 245);

    auto g_0_xxyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 246);

    auto g_0_xxyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 247);

    auto g_0_xxyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 248);

    auto g_0_xxyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 249);

    auto g_0_xxyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 250);

    auto g_0_xxyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 251);

    auto g_0_xxzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 324);

    auto g_0_xxzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 325);

    auto g_0_xxzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 326);

    auto g_0_xxzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 327);

    auto g_0_xxzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 328);

    auto g_0_xxzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 329);

    auto g_0_xxzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 330);

    auto g_0_xxzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 331);

    auto g_0_xxzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 332);

    auto g_0_xxzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 333);

    auto g_0_xxzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 334);

    auto g_0_xxzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 335);

    auto g_0_xxzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 336);

    auto g_0_xxzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 337);

    auto g_0_xxzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 338);

    auto g_0_xxzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 339);

    auto g_0_xxzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 340);

    auto g_0_xxzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 341);

    auto g_0_xxzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 342);

    auto g_0_xxzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 343);

    auto g_0_xxzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 344);

    auto g_0_xxzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 345);

    auto g_0_xxzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 346);

    auto g_0_xxzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 347);

    auto g_0_xxzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 348);

    auto g_0_xxzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 349);

    auto g_0_xxzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 350);

    auto g_0_xxzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 351);

    auto g_0_xxzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 352);

    auto g_0_xxzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 353);

    auto g_0_xxzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 354);

    auto g_0_xxzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 355);

    auto g_0_xxzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 356);

    auto g_0_xxzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 357);

    auto g_0_xxzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 358);

    auto g_0_xxzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 359);

    auto g_0_xyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 361);

    auto g_0_xyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 363);

    auto g_0_xyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 364);

    auto g_0_xyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 366);

    auto g_0_xyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 367);

    auto g_0_xyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 368);

    auto g_0_xyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 370);

    auto g_0_xyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 371);

    auto g_0_xyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 372);

    auto g_0_xyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 373);

    auto g_0_xyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 375);

    auto g_0_xyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 376);

    auto g_0_xyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 377);

    auto g_0_xyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 378);

    auto g_0_xyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 379);

    auto g_0_xyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 381);

    auto g_0_xyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 382);

    auto g_0_xyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 383);

    auto g_0_xyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 384);

    auto g_0_xyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 385);

    auto g_0_xyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 386);

    auto g_0_xyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 388);

    auto g_0_xyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 389);

    auto g_0_xyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 390);

    auto g_0_xyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 391);

    auto g_0_xyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 392);

    auto g_0_xyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 393);

    auto g_0_xyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 394);

    auto g_0_xyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 436);

    auto g_0_xyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 439);

    auto g_0_xyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 440);

    auto g_0_xyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 443);

    auto g_0_xyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 444);

    auto g_0_xyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 445);

    auto g_0_xyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 448);

    auto g_0_xyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 449);

    auto g_0_xyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 450);

    auto g_0_xyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 451);

    auto g_0_xyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 454);

    auto g_0_xyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 455);

    auto g_0_xyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 456);

    auto g_0_xyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 457);

    auto g_0_xyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 458);

    auto g_0_xyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 461);

    auto g_0_xyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 462);

    auto g_0_xyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 463);

    auto g_0_xyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 464);

    auto g_0_xyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 465);

    auto g_0_xyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 466);

    auto g_0_xzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 506);

    auto g_0_xzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 508);

    auto g_0_xzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 509);

    auto g_0_xzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 511);

    auto g_0_xzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 512);

    auto g_0_xzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 513);

    auto g_0_xzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 515);

    auto g_0_xzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 516);

    auto g_0_xzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 517);

    auto g_0_xzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 518);

    auto g_0_xzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 520);

    auto g_0_xzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 521);

    auto g_0_xzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 522);

    auto g_0_xzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 523);

    auto g_0_xzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 524);

    auto g_0_xzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 526);

    auto g_0_xzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 527);

    auto g_0_xzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 528);

    auto g_0_xzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 529);

    auto g_0_xzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 530);

    auto g_0_xzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 531);

    auto g_0_xzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 533);

    auto g_0_xzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 534);

    auto g_0_xzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 535);

    auto g_0_xzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 536);

    auto g_0_xzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 537);

    auto g_0_xzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 538);

    auto g_0_xzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 539);

    auto g_0_yyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 540);

    auto g_0_yyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 541);

    auto g_0_yyyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 542);

    auto g_0_yyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 543);

    auto g_0_yyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 544);

    auto g_0_yyyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 545);

    auto g_0_yyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 546);

    auto g_0_yyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 547);

    auto g_0_yyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 548);

    auto g_0_yyyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 549);

    auto g_0_yyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 550);

    auto g_0_yyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 551);

    auto g_0_yyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 552);

    auto g_0_yyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 553);

    auto g_0_yyyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 554);

    auto g_0_yyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 555);

    auto g_0_yyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 556);

    auto g_0_yyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 557);

    auto g_0_yyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 558);

    auto g_0_yyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 559);

    auto g_0_yyyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 560);

    auto g_0_yyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 561);

    auto g_0_yyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 562);

    auto g_0_yyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 563);

    auto g_0_yyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 564);

    auto g_0_yyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 565);

    auto g_0_yyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 566);

    auto g_0_yyyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 567);

    auto g_0_yyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 568);

    auto g_0_yyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 569);

    auto g_0_yyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 570);

    auto g_0_yyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 571);

    auto g_0_yyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 572);

    auto g_0_yyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 573);

    auto g_0_yyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 574);

    auto g_0_yyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 575);

    auto g_0_yyyyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 578);

    auto g_0_yyyyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 580);

    auto g_0_yyyyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 581);

    auto g_0_yyyyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 583);

    auto g_0_yyyyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 584);

    auto g_0_yyyyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 585);

    auto g_0_yyyyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 587);

    auto g_0_yyyyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 588);

    auto g_0_yyyyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 589);

    auto g_0_yyyyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 590);

    auto g_0_yyyyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 592);

    auto g_0_yyyyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 593);

    auto g_0_yyyyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 594);

    auto g_0_yyyyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 595);

    auto g_0_yyyyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 596);

    auto g_0_yyyyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 598);

    auto g_0_yyyyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 599);

    auto g_0_yyyyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 600);

    auto g_0_yyyyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 601);

    auto g_0_yyyyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 602);

    auto g_0_yyyyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 603);

    auto g_0_yyyyz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 605);

    auto g_0_yyyyz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 606);

    auto g_0_yyyyz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 607);

    auto g_0_yyyyz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 608);

    auto g_0_yyyyz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 609);

    auto g_0_yyyyz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 610);

    auto g_0_yyyyz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 611);

    auto g_0_yyyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 612);

    auto g_0_yyyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 613);

    auto g_0_yyyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 614);

    auto g_0_yyyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 615);

    auto g_0_yyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 616);

    auto g_0_yyyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 617);

    auto g_0_yyyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 618);

    auto g_0_yyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 619);

    auto g_0_yyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 620);

    auto g_0_yyyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 621);

    auto g_0_yyyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 622);

    auto g_0_yyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 623);

    auto g_0_yyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 624);

    auto g_0_yyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 625);

    auto g_0_yyyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 626);

    auto g_0_yyyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 627);

    auto g_0_yyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 628);

    auto g_0_yyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 629);

    auto g_0_yyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 630);

    auto g_0_yyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 631);

    auto g_0_yyyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 632);

    auto g_0_yyyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 633);

    auto g_0_yyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 634);

    auto g_0_yyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 635);

    auto g_0_yyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 636);

    auto g_0_yyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 637);

    auto g_0_yyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 638);

    auto g_0_yyyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 639);

    auto g_0_yyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 640);

    auto g_0_yyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 641);

    auto g_0_yyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 642);

    auto g_0_yyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 643);

    auto g_0_yyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 644);

    auto g_0_yyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 645);

    auto g_0_yyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 646);

    auto g_0_yyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 647);

    auto g_0_yyzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 648);

    auto g_0_yyzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 649);

    auto g_0_yyzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 650);

    auto g_0_yyzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 651);

    auto g_0_yyzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 652);

    auto g_0_yyzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 653);

    auto g_0_yyzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 654);

    auto g_0_yyzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 655);

    auto g_0_yyzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 656);

    auto g_0_yyzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 657);

    auto g_0_yyzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 658);

    auto g_0_yyzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 659);

    auto g_0_yyzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 660);

    auto g_0_yyzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 661);

    auto g_0_yyzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 662);

    auto g_0_yyzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 663);

    auto g_0_yyzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 664);

    auto g_0_yyzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 665);

    auto g_0_yyzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 666);

    auto g_0_yyzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 667);

    auto g_0_yyzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 668);

    auto g_0_yyzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 669);

    auto g_0_yyzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 670);

    auto g_0_yyzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 671);

    auto g_0_yyzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 672);

    auto g_0_yyzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 673);

    auto g_0_yyzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 674);

    auto g_0_yyzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 675);

    auto g_0_yyzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 676);

    auto g_0_yyzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 677);

    auto g_0_yyzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 678);

    auto g_0_yyzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 679);

    auto g_0_yyzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 680);

    auto g_0_yyzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 681);

    auto g_0_yyzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 682);

    auto g_0_yyzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 683);

    auto g_0_yzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 685);

    auto g_0_yzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 686);

    auto g_0_yzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 687);

    auto g_0_yzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 688);

    auto g_0_yzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 689);

    auto g_0_yzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 690);

    auto g_0_yzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 691);

    auto g_0_yzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 692);

    auto g_0_yzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 693);

    auto g_0_yzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 694);

    auto g_0_yzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 695);

    auto g_0_yzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 696);

    auto g_0_yzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 697);

    auto g_0_yzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 698);

    auto g_0_yzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 699);

    auto g_0_yzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 700);

    auto g_0_yzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 701);

    auto g_0_yzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 702);

    auto g_0_yzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 703);

    auto g_0_yzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 704);

    auto g_0_yzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 705);

    auto g_0_yzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 706);

    auto g_0_yzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 707);

    auto g_0_yzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 708);

    auto g_0_yzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 709);

    auto g_0_yzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 710);

    auto g_0_yzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 711);

    auto g_0_yzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 712);

    auto g_0_yzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 713);

    auto g_0_yzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 714);

    auto g_0_yzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 715);

    auto g_0_yzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 716);

    auto g_0_yzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 717);

    auto g_0_yzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 718);

    auto g_0_yzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 719);

    auto g_0_zzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_shsk + 720);

    auto g_0_zzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_shsk + 721);

    auto g_0_zzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_shsk + 722);

    auto g_0_zzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_shsk + 723);

    auto g_0_zzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_shsk + 724);

    auto g_0_zzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_shsk + 725);

    auto g_0_zzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_shsk + 726);

    auto g_0_zzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_shsk + 727);

    auto g_0_zzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_shsk + 728);

    auto g_0_zzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_shsk + 729);

    auto g_0_zzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_shsk + 730);

    auto g_0_zzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_shsk + 731);

    auto g_0_zzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_shsk + 732);

    auto g_0_zzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_shsk + 733);

    auto g_0_zzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_shsk + 734);

    auto g_0_zzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 735);

    auto g_0_zzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 736);

    auto g_0_zzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 737);

    auto g_0_zzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 738);

    auto g_0_zzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 739);

    auto g_0_zzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 740);

    auto g_0_zzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 741);

    auto g_0_zzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 742);

    auto g_0_zzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 743);

    auto g_0_zzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 744);

    auto g_0_zzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 745);

    auto g_0_zzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 746);

    auto g_0_zzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 747);

    auto g_0_zzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_shsk + 748);

    auto g_0_zzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_shsk + 749);

    auto g_0_zzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_shsk + 750);

    auto g_0_zzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_shsk + 751);

    auto g_0_zzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_shsk + 752);

    auto g_0_zzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 753);

    auto g_0_zzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 754);

    auto g_0_zzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_shsk + 755);

    /// Set up components of auxilary buffer : SHSL

    auto g_0_xxxxx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl);

    auto g_0_xxxxx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 1);

    auto g_0_xxxxx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 2);

    auto g_0_xxxxx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 3);

    auto g_0_xxxxx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 4);

    auto g_0_xxxxx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 5);

    auto g_0_xxxxx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 6);

    auto g_0_xxxxx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 7);

    auto g_0_xxxxx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 8);

    auto g_0_xxxxx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 9);

    auto g_0_xxxxx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 10);

    auto g_0_xxxxx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 11);

    auto g_0_xxxxx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 12);

    auto g_0_xxxxx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 13);

    auto g_0_xxxxx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 14);

    auto g_0_xxxxx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 15);

    auto g_0_xxxxx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 16);

    auto g_0_xxxxx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 17);

    auto g_0_xxxxx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 18);

    auto g_0_xxxxx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 19);

    auto g_0_xxxxx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 20);

    auto g_0_xxxxx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 21);

    auto g_0_xxxxx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 22);

    auto g_0_xxxxx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 23);

    auto g_0_xxxxx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 24);

    auto g_0_xxxxx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 25);

    auto g_0_xxxxx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 26);

    auto g_0_xxxxx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 27);

    auto g_0_xxxxx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 28);

    auto g_0_xxxxx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 29);

    auto g_0_xxxxx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 30);

    auto g_0_xxxxx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 31);

    auto g_0_xxxxx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 32);

    auto g_0_xxxxx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 33);

    auto g_0_xxxxx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 34);

    auto g_0_xxxxx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 35);

    auto g_0_xxxxx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 36);

    auto g_0_xxxxx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 37);

    auto g_0_xxxxx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 38);

    auto g_0_xxxxx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 39);

    auto g_0_xxxxx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 40);

    auto g_0_xxxxx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 41);

    auto g_0_xxxxx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 42);

    auto g_0_xxxxx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 43);

    auto g_0_xxxxx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 44);

    auto g_0_xxxxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 45);

    auto g_0_xxxxy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 46);

    auto g_0_xxxxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 47);

    auto g_0_xxxxy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 48);

    auto g_0_xxxxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 50);

    auto g_0_xxxxy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 51);

    auto g_0_xxxxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 54);

    auto g_0_xxxxy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 55);

    auto g_0_xxxxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 59);

    auto g_0_xxxxy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 60);

    auto g_0_xxxxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 65);

    auto g_0_xxxxy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 66);

    auto g_0_xxxxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 72);

    auto g_0_xxxxy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 73);

    auto g_0_xxxxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 80);

    auto g_0_xxxxy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 81);

    auto g_0_xxxxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 90);

    auto g_0_xxxxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 91);

    auto g_0_xxxxz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 92);

    auto g_0_xxxxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 93);

    auto g_0_xxxxz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 94);

    auto g_0_xxxxz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 95);

    auto g_0_xxxxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 96);

    auto g_0_xxxxz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 97);

    auto g_0_xxxxz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 98);

    auto g_0_xxxxz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 99);

    auto g_0_xxxxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 100);

    auto g_0_xxxxz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 101);

    auto g_0_xxxxz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 102);

    auto g_0_xxxxz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 103);

    auto g_0_xxxxz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 104);

    auto g_0_xxxxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 105);

    auto g_0_xxxxz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 106);

    auto g_0_xxxxz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 107);

    auto g_0_xxxxz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 108);

    auto g_0_xxxxz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 109);

    auto g_0_xxxxz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 110);

    auto g_0_xxxxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 111);

    auto g_0_xxxxz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 112);

    auto g_0_xxxxz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 113);

    auto g_0_xxxxz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 114);

    auto g_0_xxxxz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 115);

    auto g_0_xxxxz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 116);

    auto g_0_xxxxz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 117);

    auto g_0_xxxxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 118);

    auto g_0_xxxxz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 119);

    auto g_0_xxxxz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 120);

    auto g_0_xxxxz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 121);

    auto g_0_xxxxz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 122);

    auto g_0_xxxxz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 123);

    auto g_0_xxxxz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 124);

    auto g_0_xxxxz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 125);

    auto g_0_xxxxz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 127);

    auto g_0_xxxxz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 128);

    auto g_0_xxxxz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 129);

    auto g_0_xxxxz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 130);

    auto g_0_xxxxz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 131);

    auto g_0_xxxxz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 132);

    auto g_0_xxxxz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 133);

    auto g_0_xxxxz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 134);

    auto g_0_xxxyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 135);

    auto g_0_xxxyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 136);

    auto g_0_xxxyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 137);

    auto g_0_xxxyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 138);

    auto g_0_xxxyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 139);

    auto g_0_xxxyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 140);

    auto g_0_xxxyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 141);

    auto g_0_xxxyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 142);

    auto g_0_xxxyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 143);

    auto g_0_xxxyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 144);

    auto g_0_xxxyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 145);

    auto g_0_xxxyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 146);

    auto g_0_xxxyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 147);

    auto g_0_xxxyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 148);

    auto g_0_xxxyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 149);

    auto g_0_xxxyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 150);

    auto g_0_xxxyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 151);

    auto g_0_xxxyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 152);

    auto g_0_xxxyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 153);

    auto g_0_xxxyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 154);

    auto g_0_xxxyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 155);

    auto g_0_xxxyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 156);

    auto g_0_xxxyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 157);

    auto g_0_xxxyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 158);

    auto g_0_xxxyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 159);

    auto g_0_xxxyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 160);

    auto g_0_xxxyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 161);

    auto g_0_xxxyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 162);

    auto g_0_xxxyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 163);

    auto g_0_xxxyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 164);

    auto g_0_xxxyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 165);

    auto g_0_xxxyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 166);

    auto g_0_xxxyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 167);

    auto g_0_xxxyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 168);

    auto g_0_xxxyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 169);

    auto g_0_xxxyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 170);

    auto g_0_xxxyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 171);

    auto g_0_xxxyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 172);

    auto g_0_xxxyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 173);

    auto g_0_xxxyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 174);

    auto g_0_xxxyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 175);

    auto g_0_xxxyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 176);

    auto g_0_xxxyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 177);

    auto g_0_xxxyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 178);

    auto g_0_xxxyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 179);

    auto g_0_xxxzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 225);

    auto g_0_xxxzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 226);

    auto g_0_xxxzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 227);

    auto g_0_xxxzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 228);

    auto g_0_xxxzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 229);

    auto g_0_xxxzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 230);

    auto g_0_xxxzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 231);

    auto g_0_xxxzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 232);

    auto g_0_xxxzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 233);

    auto g_0_xxxzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 234);

    auto g_0_xxxzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 235);

    auto g_0_xxxzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 236);

    auto g_0_xxxzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 237);

    auto g_0_xxxzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 238);

    auto g_0_xxxzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 239);

    auto g_0_xxxzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 240);

    auto g_0_xxxzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 241);

    auto g_0_xxxzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 242);

    auto g_0_xxxzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 243);

    auto g_0_xxxzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 244);

    auto g_0_xxxzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 245);

    auto g_0_xxxzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 246);

    auto g_0_xxxzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 247);

    auto g_0_xxxzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 248);

    auto g_0_xxxzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 249);

    auto g_0_xxxzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 250);

    auto g_0_xxxzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 251);

    auto g_0_xxxzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 252);

    auto g_0_xxxzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 253);

    auto g_0_xxxzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 254);

    auto g_0_xxxzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 255);

    auto g_0_xxxzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 256);

    auto g_0_xxxzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 257);

    auto g_0_xxxzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 258);

    auto g_0_xxxzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 259);

    auto g_0_xxxzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 260);

    auto g_0_xxxzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 261);

    auto g_0_xxxzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 262);

    auto g_0_xxxzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 263);

    auto g_0_xxxzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 264);

    auto g_0_xxxzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 265);

    auto g_0_xxxzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 266);

    auto g_0_xxxzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 267);

    auto g_0_xxxzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 268);

    auto g_0_xxxzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 269);

    auto g_0_xxyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 270);

    auto g_0_xxyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 271);

    auto g_0_xxyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 272);

    auto g_0_xxyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 273);

    auto g_0_xxyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 274);

    auto g_0_xxyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 275);

    auto g_0_xxyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 276);

    auto g_0_xxyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 277);

    auto g_0_xxyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 278);

    auto g_0_xxyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 279);

    auto g_0_xxyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 280);

    auto g_0_xxyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 281);

    auto g_0_xxyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 282);

    auto g_0_xxyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 283);

    auto g_0_xxyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 284);

    auto g_0_xxyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 285);

    auto g_0_xxyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 286);

    auto g_0_xxyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 287);

    auto g_0_xxyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 288);

    auto g_0_xxyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 289);

    auto g_0_xxyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 290);

    auto g_0_xxyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 291);

    auto g_0_xxyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 292);

    auto g_0_xxyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 293);

    auto g_0_xxyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 294);

    auto g_0_xxyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 295);

    auto g_0_xxyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 296);

    auto g_0_xxyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 297);

    auto g_0_xxyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 298);

    auto g_0_xxyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 299);

    auto g_0_xxyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 300);

    auto g_0_xxyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 301);

    auto g_0_xxyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 302);

    auto g_0_xxyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 303);

    auto g_0_xxyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 304);

    auto g_0_xxyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 305);

    auto g_0_xxyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 306);

    auto g_0_xxyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 307);

    auto g_0_xxyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 308);

    auto g_0_xxyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 309);

    auto g_0_xxyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 310);

    auto g_0_xxyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 311);

    auto g_0_xxyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 312);

    auto g_0_xxyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 313);

    auto g_0_xxyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 314);

    auto g_0_xxyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 316);

    auto g_0_xxyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 318);

    auto g_0_xxyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 321);

    auto g_0_xxyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 325);

    auto g_0_xxyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 330);

    auto g_0_xxyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 336);

    auto g_0_xxyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 343);

    auto g_0_xxyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 360);

    auto g_0_xxyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 362);

    auto g_0_xxyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 365);

    auto g_0_xxyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 369);

    auto g_0_xxyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 374);

    auto g_0_xxyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 380);

    auto g_0_xxyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 387);

    auto g_0_xxyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 395);

    auto g_0_xxzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 405);

    auto g_0_xxzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 406);

    auto g_0_xxzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 407);

    auto g_0_xxzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 408);

    auto g_0_xxzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 409);

    auto g_0_xxzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 410);

    auto g_0_xxzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 411);

    auto g_0_xxzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 412);

    auto g_0_xxzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 413);

    auto g_0_xxzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 414);

    auto g_0_xxzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 415);

    auto g_0_xxzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 416);

    auto g_0_xxzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 417);

    auto g_0_xxzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 418);

    auto g_0_xxzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 419);

    auto g_0_xxzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 420);

    auto g_0_xxzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 421);

    auto g_0_xxzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 422);

    auto g_0_xxzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 423);

    auto g_0_xxzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 424);

    auto g_0_xxzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 425);

    auto g_0_xxzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 426);

    auto g_0_xxzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 427);

    auto g_0_xxzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 428);

    auto g_0_xxzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 429);

    auto g_0_xxzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 430);

    auto g_0_xxzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 431);

    auto g_0_xxzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 432);

    auto g_0_xxzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 433);

    auto g_0_xxzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 434);

    auto g_0_xxzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 435);

    auto g_0_xxzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 436);

    auto g_0_xxzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 437);

    auto g_0_xxzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 438);

    auto g_0_xxzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 439);

    auto g_0_xxzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 440);

    auto g_0_xxzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 441);

    auto g_0_xxzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 442);

    auto g_0_xxzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 443);

    auto g_0_xxzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 444);

    auto g_0_xxzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 445);

    auto g_0_xxzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 446);

    auto g_0_xxzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 447);

    auto g_0_xxzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 448);

    auto g_0_xxzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 449);

    auto g_0_xyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 450);

    auto g_0_xyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 451);

    auto g_0_xyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 453);

    auto g_0_xyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 454);

    auto g_0_xyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 456);

    auto g_0_xyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 457);

    auto g_0_xyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 458);

    auto g_0_xyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 460);

    auto g_0_xyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 461);

    auto g_0_xyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 462);

    auto g_0_xyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 463);

    auto g_0_xyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 465);

    auto g_0_xyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 466);

    auto g_0_xyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 467);

    auto g_0_xyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 468);

    auto g_0_xyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 469);

    auto g_0_xyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 471);

    auto g_0_xyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 472);

    auto g_0_xyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 473);

    auto g_0_xyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 474);

    auto g_0_xyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 475);

    auto g_0_xyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 476);

    auto g_0_xyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 478);

    auto g_0_xyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 479);

    auto g_0_xyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 480);

    auto g_0_xyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 481);

    auto g_0_xyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 482);

    auto g_0_xyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 483);

    auto g_0_xyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 484);

    auto g_0_xyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 486);

    auto g_0_xyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 487);

    auto g_0_xyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 488);

    auto g_0_xyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 489);

    auto g_0_xyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 490);

    auto g_0_xyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 491);

    auto g_0_xyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 492);

    auto g_0_xyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 493);

    auto g_0_xyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 494);

    auto g_0_xyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 544);

    auto g_0_xyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 547);

    auto g_0_xyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 548);

    auto g_0_xyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 551);

    auto g_0_xyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 552);

    auto g_0_xyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 553);

    auto g_0_xyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 556);

    auto g_0_xyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 557);

    auto g_0_xyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 558);

    auto g_0_xyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 559);

    auto g_0_xyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 562);

    auto g_0_xyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 563);

    auto g_0_xyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 564);

    auto g_0_xyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 565);

    auto g_0_xyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 566);

    auto g_0_xyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 569);

    auto g_0_xyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 570);

    auto g_0_xyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 571);

    auto g_0_xyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 572);

    auto g_0_xyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 573);

    auto g_0_xyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 574);

    auto g_0_xyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 576);

    auto g_0_xyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 577);

    auto g_0_xyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 578);

    auto g_0_xyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 579);

    auto g_0_xyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 580);

    auto g_0_xyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 581);

    auto g_0_xyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 582);

    auto g_0_xyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 583);

    auto g_0_xyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 584);

    auto g_0_xzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 630);

    auto g_0_xzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 632);

    auto g_0_xzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 634);

    auto g_0_xzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 635);

    auto g_0_xzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 637);

    auto g_0_xzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 638);

    auto g_0_xzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 639);

    auto g_0_xzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 641);

    auto g_0_xzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 642);

    auto g_0_xzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 643);

    auto g_0_xzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 644);

    auto g_0_xzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 646);

    auto g_0_xzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 647);

    auto g_0_xzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 648);

    auto g_0_xzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 649);

    auto g_0_xzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 650);

    auto g_0_xzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 652);

    auto g_0_xzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 653);

    auto g_0_xzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 654);

    auto g_0_xzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 655);

    auto g_0_xzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 656);

    auto g_0_xzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 657);

    auto g_0_xzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 659);

    auto g_0_xzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 660);

    auto g_0_xzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 661);

    auto g_0_xzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 662);

    auto g_0_xzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 663);

    auto g_0_xzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 664);

    auto g_0_xzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 665);

    auto g_0_xzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 666);

    auto g_0_xzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 667);

    auto g_0_xzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 668);

    auto g_0_xzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 669);

    auto g_0_xzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 670);

    auto g_0_xzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 671);

    auto g_0_xzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 672);

    auto g_0_xzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 673);

    auto g_0_xzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 674);

    auto g_0_yyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 675);

    auto g_0_yyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 676);

    auto g_0_yyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 677);

    auto g_0_yyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 678);

    auto g_0_yyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 679);

    auto g_0_yyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 680);

    auto g_0_yyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 681);

    auto g_0_yyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 682);

    auto g_0_yyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 683);

    auto g_0_yyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 684);

    auto g_0_yyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 685);

    auto g_0_yyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 686);

    auto g_0_yyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 687);

    auto g_0_yyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 688);

    auto g_0_yyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 689);

    auto g_0_yyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 690);

    auto g_0_yyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 691);

    auto g_0_yyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 692);

    auto g_0_yyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 693);

    auto g_0_yyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 694);

    auto g_0_yyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 695);

    auto g_0_yyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 696);

    auto g_0_yyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 697);

    auto g_0_yyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 698);

    auto g_0_yyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 699);

    auto g_0_yyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 700);

    auto g_0_yyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 701);

    auto g_0_yyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 702);

    auto g_0_yyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 703);

    auto g_0_yyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 704);

    auto g_0_yyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 705);

    auto g_0_yyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 706);

    auto g_0_yyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 707);

    auto g_0_yyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 708);

    auto g_0_yyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 709);

    auto g_0_yyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 710);

    auto g_0_yyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 711);

    auto g_0_yyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 712);

    auto g_0_yyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 713);

    auto g_0_yyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 714);

    auto g_0_yyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 715);

    auto g_0_yyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 716);

    auto g_0_yyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 717);

    auto g_0_yyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 718);

    auto g_0_yyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 719);

    auto g_0_yyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 721);

    auto g_0_yyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 722);

    auto g_0_yyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 723);

    auto g_0_yyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 724);

    auto g_0_yyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 725);

    auto g_0_yyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 726);

    auto g_0_yyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 727);

    auto g_0_yyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 728);

    auto g_0_yyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 729);

    auto g_0_yyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 730);

    auto g_0_yyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 731);

    auto g_0_yyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 732);

    auto g_0_yyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 733);

    auto g_0_yyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 734);

    auto g_0_yyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 735);

    auto g_0_yyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 736);

    auto g_0_yyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 737);

    auto g_0_yyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 738);

    auto g_0_yyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 739);

    auto g_0_yyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 740);

    auto g_0_yyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 741);

    auto g_0_yyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 742);

    auto g_0_yyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 743);

    auto g_0_yyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 744);

    auto g_0_yyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 745);

    auto g_0_yyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 746);

    auto g_0_yyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 747);

    auto g_0_yyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 748);

    auto g_0_yyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 749);

    auto g_0_yyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 750);

    auto g_0_yyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 751);

    auto g_0_yyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 752);

    auto g_0_yyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 753);

    auto g_0_yyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 754);

    auto g_0_yyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 755);

    auto g_0_yyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 756);

    auto g_0_yyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 757);

    auto g_0_yyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 758);

    auto g_0_yyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 759);

    auto g_0_yyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 760);

    auto g_0_yyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 761);

    auto g_0_yyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 762);

    auto g_0_yyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 763);

    auto g_0_yyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 764);

    auto g_0_yyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 765);

    auto g_0_yyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 766);

    auto g_0_yyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 767);

    auto g_0_yyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 768);

    auto g_0_yyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 769);

    auto g_0_yyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 770);

    auto g_0_yyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 771);

    auto g_0_yyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 772);

    auto g_0_yyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 773);

    auto g_0_yyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 774);

    auto g_0_yyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 775);

    auto g_0_yyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 776);

    auto g_0_yyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 777);

    auto g_0_yyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 778);

    auto g_0_yyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 779);

    auto g_0_yyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 780);

    auto g_0_yyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 781);

    auto g_0_yyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 782);

    auto g_0_yyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 783);

    auto g_0_yyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 784);

    auto g_0_yyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 785);

    auto g_0_yyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 786);

    auto g_0_yyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 787);

    auto g_0_yyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 788);

    auto g_0_yyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 789);

    auto g_0_yyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 790);

    auto g_0_yyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 791);

    auto g_0_yyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 792);

    auto g_0_yyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 793);

    auto g_0_yyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 794);

    auto g_0_yyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 795);

    auto g_0_yyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 796);

    auto g_0_yyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 797);

    auto g_0_yyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 798);

    auto g_0_yyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 799);

    auto g_0_yyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 800);

    auto g_0_yyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 801);

    auto g_0_yyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 802);

    auto g_0_yyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 803);

    auto g_0_yyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 804);

    auto g_0_yyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 805);

    auto g_0_yyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 806);

    auto g_0_yyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 807);

    auto g_0_yyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 808);

    auto g_0_yyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 809);

    auto g_0_yyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 810);

    auto g_0_yyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 811);

    auto g_0_yyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 812);

    auto g_0_yyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 813);

    auto g_0_yyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 814);

    auto g_0_yyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 815);

    auto g_0_yyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 816);

    auto g_0_yyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 817);

    auto g_0_yyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 818);

    auto g_0_yyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 819);

    auto g_0_yyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 820);

    auto g_0_yyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 821);

    auto g_0_yyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 822);

    auto g_0_yyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 823);

    auto g_0_yyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 824);

    auto g_0_yyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 825);

    auto g_0_yyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 826);

    auto g_0_yyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 827);

    auto g_0_yyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 828);

    auto g_0_yyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 829);

    auto g_0_yyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 830);

    auto g_0_yyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 831);

    auto g_0_yyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 832);

    auto g_0_yyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 833);

    auto g_0_yyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 834);

    auto g_0_yyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 835);

    auto g_0_yyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 836);

    auto g_0_yyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 837);

    auto g_0_yyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 838);

    auto g_0_yyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 839);

    auto g_0_yyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 840);

    auto g_0_yyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 841);

    auto g_0_yyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 842);

    auto g_0_yyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 843);

    auto g_0_yyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 844);

    auto g_0_yyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 845);

    auto g_0_yyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 846);

    auto g_0_yyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 847);

    auto g_0_yyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 848);

    auto g_0_yyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 849);

    auto g_0_yyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 850);

    auto g_0_yyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 851);

    auto g_0_yyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 852);

    auto g_0_yyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 853);

    auto g_0_yyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 854);

    auto g_0_yzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 855);

    auto g_0_yzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 856);

    auto g_0_yzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 857);

    auto g_0_yzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 858);

    auto g_0_yzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 859);

    auto g_0_yzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 860);

    auto g_0_yzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 861);

    auto g_0_yzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 862);

    auto g_0_yzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 863);

    auto g_0_yzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 864);

    auto g_0_yzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 865);

    auto g_0_yzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 866);

    auto g_0_yzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 867);

    auto g_0_yzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 868);

    auto g_0_yzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 869);

    auto g_0_yzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 870);

    auto g_0_yzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 871);

    auto g_0_yzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 872);

    auto g_0_yzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 873);

    auto g_0_yzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 874);

    auto g_0_yzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 875);

    auto g_0_yzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 876);

    auto g_0_yzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 877);

    auto g_0_yzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 878);

    auto g_0_yzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 879);

    auto g_0_yzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 880);

    auto g_0_yzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 881);

    auto g_0_yzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 882);

    auto g_0_yzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 883);

    auto g_0_yzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 884);

    auto g_0_yzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 885);

    auto g_0_yzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 886);

    auto g_0_yzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 887);

    auto g_0_yzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 888);

    auto g_0_yzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 889);

    auto g_0_yzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 890);

    auto g_0_yzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 891);

    auto g_0_yzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 892);

    auto g_0_yzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 893);

    auto g_0_yzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 894);

    auto g_0_yzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 895);

    auto g_0_yzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 896);

    auto g_0_yzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 897);

    auto g_0_yzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 898);

    auto g_0_yzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 899);

    auto g_0_zzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_shsl + 900);

    auto g_0_zzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_shsl + 901);

    auto g_0_zzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_shsl + 902);

    auto g_0_zzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_shsl + 903);

    auto g_0_zzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_shsl + 904);

    auto g_0_zzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_shsl + 905);

    auto g_0_zzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_shsl + 906);

    auto g_0_zzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_shsl + 907);

    auto g_0_zzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_shsl + 908);

    auto g_0_zzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_shsl + 909);

    auto g_0_zzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_shsl + 910);

    auto g_0_zzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_shsl + 911);

    auto g_0_zzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_shsl + 912);

    auto g_0_zzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_shsl + 913);

    auto g_0_zzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_shsl + 914);

    auto g_0_zzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 915);

    auto g_0_zzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 916);

    auto g_0_zzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 917);

    auto g_0_zzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 918);

    auto g_0_zzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 919);

    auto g_0_zzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 920);

    auto g_0_zzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 921);

    auto g_0_zzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 922);

    auto g_0_zzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 923);

    auto g_0_zzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 924);

    auto g_0_zzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 925);

    auto g_0_zzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 926);

    auto g_0_zzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 927);

    auto g_0_zzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 928);

    auto g_0_zzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 929);

    auto g_0_zzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 930);

    auto g_0_zzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 931);

    auto g_0_zzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 932);

    auto g_0_zzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 933);

    auto g_0_zzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 934);

    auto g_0_zzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 935);

    auto g_0_zzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_shsl + 936);

    auto g_0_zzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_shsl + 937);

    auto g_0_zzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_shsl + 938);

    auto g_0_zzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_shsl + 939);

    auto g_0_zzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_shsl + 940);

    auto g_0_zzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 941);

    auto g_0_zzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 942);

    auto g_0_zzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 943);

    auto g_0_zzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_shsl + 944);

    /// Set up components of auxilary buffer : SHSL

    auto g_0_xxxxx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl);

    auto g_0_xxxxx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 1);

    auto g_0_xxxxx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 2);

    auto g_0_xxxxx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 3);

    auto g_0_xxxxx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 4);

    auto g_0_xxxxx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 5);

    auto g_0_xxxxx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 6);

    auto g_0_xxxxx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 7);

    auto g_0_xxxxx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 8);

    auto g_0_xxxxx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 9);

    auto g_0_xxxxx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 10);

    auto g_0_xxxxx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 11);

    auto g_0_xxxxx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 12);

    auto g_0_xxxxx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 13);

    auto g_0_xxxxx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 14);

    auto g_0_xxxxx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 15);

    auto g_0_xxxxx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 16);

    auto g_0_xxxxx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 17);

    auto g_0_xxxxx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 18);

    auto g_0_xxxxx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 19);

    auto g_0_xxxxx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 20);

    auto g_0_xxxxx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 21);

    auto g_0_xxxxx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 22);

    auto g_0_xxxxx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 23);

    auto g_0_xxxxx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 24);

    auto g_0_xxxxx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 25);

    auto g_0_xxxxx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 26);

    auto g_0_xxxxx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 27);

    auto g_0_xxxxx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 28);

    auto g_0_xxxxx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 29);

    auto g_0_xxxxx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 30);

    auto g_0_xxxxx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 31);

    auto g_0_xxxxx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 32);

    auto g_0_xxxxx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 33);

    auto g_0_xxxxx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 34);

    auto g_0_xxxxx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 35);

    auto g_0_xxxxx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 36);

    auto g_0_xxxxx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 37);

    auto g_0_xxxxx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 38);

    auto g_0_xxxxx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 39);

    auto g_0_xxxxx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 40);

    auto g_0_xxxxx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 41);

    auto g_0_xxxxx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 42);

    auto g_0_xxxxx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 43);

    auto g_0_xxxxx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 44);

    auto g_0_xxxxy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 45);

    auto g_0_xxxxy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 46);

    auto g_0_xxxxy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 47);

    auto g_0_xxxxy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 48);

    auto g_0_xxxxy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 50);

    auto g_0_xxxxy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 51);

    auto g_0_xxxxy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 54);

    auto g_0_xxxxy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 55);

    auto g_0_xxxxy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 59);

    auto g_0_xxxxy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 60);

    auto g_0_xxxxy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 65);

    auto g_0_xxxxy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 66);

    auto g_0_xxxxy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 72);

    auto g_0_xxxxy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 73);

    auto g_0_xxxxy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 80);

    auto g_0_xxxxy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 81);

    auto g_0_xxxxz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 90);

    auto g_0_xxxxz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 91);

    auto g_0_xxxxz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 92);

    auto g_0_xxxxz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 93);

    auto g_0_xxxxz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 94);

    auto g_0_xxxxz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 95);

    auto g_0_xxxxz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 96);

    auto g_0_xxxxz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 97);

    auto g_0_xxxxz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 98);

    auto g_0_xxxxz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 99);

    auto g_0_xxxxz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 100);

    auto g_0_xxxxz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 101);

    auto g_0_xxxxz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 102);

    auto g_0_xxxxz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 103);

    auto g_0_xxxxz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 104);

    auto g_0_xxxxz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 105);

    auto g_0_xxxxz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 106);

    auto g_0_xxxxz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 107);

    auto g_0_xxxxz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 108);

    auto g_0_xxxxz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 109);

    auto g_0_xxxxz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 110);

    auto g_0_xxxxz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 111);

    auto g_0_xxxxz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 112);

    auto g_0_xxxxz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 113);

    auto g_0_xxxxz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 114);

    auto g_0_xxxxz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 115);

    auto g_0_xxxxz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 116);

    auto g_0_xxxxz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 117);

    auto g_0_xxxxz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 118);

    auto g_0_xxxxz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 119);

    auto g_0_xxxxz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 120);

    auto g_0_xxxxz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 121);

    auto g_0_xxxxz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 122);

    auto g_0_xxxxz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 123);

    auto g_0_xxxxz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 124);

    auto g_0_xxxxz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 125);

    auto g_0_xxxxz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 127);

    auto g_0_xxxxz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 128);

    auto g_0_xxxxz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 129);

    auto g_0_xxxxz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 130);

    auto g_0_xxxxz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 131);

    auto g_0_xxxxz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 132);

    auto g_0_xxxxz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 133);

    auto g_0_xxxxz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 134);

    auto g_0_xxxyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 135);

    auto g_0_xxxyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 136);

    auto g_0_xxxyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 137);

    auto g_0_xxxyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 138);

    auto g_0_xxxyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 139);

    auto g_0_xxxyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 140);

    auto g_0_xxxyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 141);

    auto g_0_xxxyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 142);

    auto g_0_xxxyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 143);

    auto g_0_xxxyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 144);

    auto g_0_xxxyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 145);

    auto g_0_xxxyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 146);

    auto g_0_xxxyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 147);

    auto g_0_xxxyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 148);

    auto g_0_xxxyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 149);

    auto g_0_xxxyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 150);

    auto g_0_xxxyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 151);

    auto g_0_xxxyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 152);

    auto g_0_xxxyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 153);

    auto g_0_xxxyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 154);

    auto g_0_xxxyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 155);

    auto g_0_xxxyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 156);

    auto g_0_xxxyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 157);

    auto g_0_xxxyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 158);

    auto g_0_xxxyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 159);

    auto g_0_xxxyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 160);

    auto g_0_xxxyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 161);

    auto g_0_xxxyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 162);

    auto g_0_xxxyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 163);

    auto g_0_xxxyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 164);

    auto g_0_xxxyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 165);

    auto g_0_xxxyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 166);

    auto g_0_xxxyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 167);

    auto g_0_xxxyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 168);

    auto g_0_xxxyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 169);

    auto g_0_xxxyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 170);

    auto g_0_xxxyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 171);

    auto g_0_xxxyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 172);

    auto g_0_xxxyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 173);

    auto g_0_xxxyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 174);

    auto g_0_xxxyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 175);

    auto g_0_xxxyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 176);

    auto g_0_xxxyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 177);

    auto g_0_xxxyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 178);

    auto g_0_xxxyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 179);

    auto g_0_xxxzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 225);

    auto g_0_xxxzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 226);

    auto g_0_xxxzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 227);

    auto g_0_xxxzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 228);

    auto g_0_xxxzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 229);

    auto g_0_xxxzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 230);

    auto g_0_xxxzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 231);

    auto g_0_xxxzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 232);

    auto g_0_xxxzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 233);

    auto g_0_xxxzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 234);

    auto g_0_xxxzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 235);

    auto g_0_xxxzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 236);

    auto g_0_xxxzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 237);

    auto g_0_xxxzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 238);

    auto g_0_xxxzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 239);

    auto g_0_xxxzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 240);

    auto g_0_xxxzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 241);

    auto g_0_xxxzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 242);

    auto g_0_xxxzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 243);

    auto g_0_xxxzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 244);

    auto g_0_xxxzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 245);

    auto g_0_xxxzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 246);

    auto g_0_xxxzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 247);

    auto g_0_xxxzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 248);

    auto g_0_xxxzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 249);

    auto g_0_xxxzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 250);

    auto g_0_xxxzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 251);

    auto g_0_xxxzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 252);

    auto g_0_xxxzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 253);

    auto g_0_xxxzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 254);

    auto g_0_xxxzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 255);

    auto g_0_xxxzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 256);

    auto g_0_xxxzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 257);

    auto g_0_xxxzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 258);

    auto g_0_xxxzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 259);

    auto g_0_xxxzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 260);

    auto g_0_xxxzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 261);

    auto g_0_xxxzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 262);

    auto g_0_xxxzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 263);

    auto g_0_xxxzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 264);

    auto g_0_xxxzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 265);

    auto g_0_xxxzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 266);

    auto g_0_xxxzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 267);

    auto g_0_xxxzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 268);

    auto g_0_xxxzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 269);

    auto g_0_xxyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 270);

    auto g_0_xxyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 271);

    auto g_0_xxyyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 272);

    auto g_0_xxyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 273);

    auto g_0_xxyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 274);

    auto g_0_xxyyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 275);

    auto g_0_xxyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 276);

    auto g_0_xxyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 277);

    auto g_0_xxyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 278);

    auto g_0_xxyyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 279);

    auto g_0_xxyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 280);

    auto g_0_xxyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 281);

    auto g_0_xxyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 282);

    auto g_0_xxyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 283);

    auto g_0_xxyyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 284);

    auto g_0_xxyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 285);

    auto g_0_xxyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 286);

    auto g_0_xxyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 287);

    auto g_0_xxyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 288);

    auto g_0_xxyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 289);

    auto g_0_xxyyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 290);

    auto g_0_xxyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 291);

    auto g_0_xxyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 292);

    auto g_0_xxyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 293);

    auto g_0_xxyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 294);

    auto g_0_xxyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 295);

    auto g_0_xxyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 296);

    auto g_0_xxyyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 297);

    auto g_0_xxyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 298);

    auto g_0_xxyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 299);

    auto g_0_xxyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 300);

    auto g_0_xxyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 301);

    auto g_0_xxyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 302);

    auto g_0_xxyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 303);

    auto g_0_xxyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 304);

    auto g_0_xxyyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 305);

    auto g_0_xxyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 306);

    auto g_0_xxyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 307);

    auto g_0_xxyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 308);

    auto g_0_xxyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 309);

    auto g_0_xxyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 310);

    auto g_0_xxyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 311);

    auto g_0_xxyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 312);

    auto g_0_xxyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 313);

    auto g_0_xxyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 314);

    auto g_0_xxyyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 316);

    auto g_0_xxyyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 318);

    auto g_0_xxyyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 321);

    auto g_0_xxyyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 325);

    auto g_0_xxyyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 330);

    auto g_0_xxyyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 336);

    auto g_0_xxyyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 343);

    auto g_0_xxyzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 360);

    auto g_0_xxyzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 362);

    auto g_0_xxyzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 365);

    auto g_0_xxyzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 369);

    auto g_0_xxyzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 374);

    auto g_0_xxyzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 380);

    auto g_0_xxyzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 387);

    auto g_0_xxyzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 395);

    auto g_0_xxzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 405);

    auto g_0_xxzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 406);

    auto g_0_xxzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 407);

    auto g_0_xxzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 408);

    auto g_0_xxzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 409);

    auto g_0_xxzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 410);

    auto g_0_xxzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 411);

    auto g_0_xxzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 412);

    auto g_0_xxzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 413);

    auto g_0_xxzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 414);

    auto g_0_xxzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 415);

    auto g_0_xxzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 416);

    auto g_0_xxzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 417);

    auto g_0_xxzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 418);

    auto g_0_xxzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 419);

    auto g_0_xxzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 420);

    auto g_0_xxzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 421);

    auto g_0_xxzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 422);

    auto g_0_xxzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 423);

    auto g_0_xxzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 424);

    auto g_0_xxzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 425);

    auto g_0_xxzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 426);

    auto g_0_xxzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 427);

    auto g_0_xxzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 428);

    auto g_0_xxzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 429);

    auto g_0_xxzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 430);

    auto g_0_xxzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 431);

    auto g_0_xxzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 432);

    auto g_0_xxzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 433);

    auto g_0_xxzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 434);

    auto g_0_xxzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 435);

    auto g_0_xxzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 436);

    auto g_0_xxzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 437);

    auto g_0_xxzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 438);

    auto g_0_xxzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 439);

    auto g_0_xxzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 440);

    auto g_0_xxzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 441);

    auto g_0_xxzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 442);

    auto g_0_xxzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 443);

    auto g_0_xxzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 444);

    auto g_0_xxzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 445);

    auto g_0_xxzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 446);

    auto g_0_xxzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 447);

    auto g_0_xxzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 448);

    auto g_0_xxzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 449);

    auto g_0_xyyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 450);

    auto g_0_xyyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 451);

    auto g_0_xyyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 453);

    auto g_0_xyyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 454);

    auto g_0_xyyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 456);

    auto g_0_xyyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 457);

    auto g_0_xyyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 458);

    auto g_0_xyyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 460);

    auto g_0_xyyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 461);

    auto g_0_xyyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 462);

    auto g_0_xyyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 463);

    auto g_0_xyyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 465);

    auto g_0_xyyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 466);

    auto g_0_xyyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 467);

    auto g_0_xyyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 468);

    auto g_0_xyyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 469);

    auto g_0_xyyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 471);

    auto g_0_xyyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 472);

    auto g_0_xyyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 473);

    auto g_0_xyyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 474);

    auto g_0_xyyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 475);

    auto g_0_xyyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 476);

    auto g_0_xyyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 478);

    auto g_0_xyyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 479);

    auto g_0_xyyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 480);

    auto g_0_xyyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 481);

    auto g_0_xyyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 482);

    auto g_0_xyyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 483);

    auto g_0_xyyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 484);

    auto g_0_xyyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 486);

    auto g_0_xyyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 487);

    auto g_0_xyyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 488);

    auto g_0_xyyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 489);

    auto g_0_xyyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 490);

    auto g_0_xyyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 491);

    auto g_0_xyyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 492);

    auto g_0_xyyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 493);

    auto g_0_xyyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 494);

    auto g_0_xyyzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 544);

    auto g_0_xyyzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 547);

    auto g_0_xyyzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 548);

    auto g_0_xyyzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 551);

    auto g_0_xyyzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 552);

    auto g_0_xyyzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 553);

    auto g_0_xyyzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 556);

    auto g_0_xyyzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 557);

    auto g_0_xyyzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 558);

    auto g_0_xyyzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 559);

    auto g_0_xyyzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 562);

    auto g_0_xyyzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 563);

    auto g_0_xyyzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 564);

    auto g_0_xyyzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 565);

    auto g_0_xyyzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 566);

    auto g_0_xyyzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 569);

    auto g_0_xyyzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 570);

    auto g_0_xyyzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 571);

    auto g_0_xyyzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 572);

    auto g_0_xyyzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 573);

    auto g_0_xyyzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 574);

    auto g_0_xyyzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 576);

    auto g_0_xyyzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 577);

    auto g_0_xyyzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 578);

    auto g_0_xyyzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 579);

    auto g_0_xyyzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 580);

    auto g_0_xyyzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 581);

    auto g_0_xyyzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 582);

    auto g_0_xyyzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 583);

    auto g_0_xyyzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 584);

    auto g_0_xzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 630);

    auto g_0_xzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 632);

    auto g_0_xzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 634);

    auto g_0_xzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 635);

    auto g_0_xzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 637);

    auto g_0_xzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 638);

    auto g_0_xzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 639);

    auto g_0_xzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 641);

    auto g_0_xzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 642);

    auto g_0_xzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 643);

    auto g_0_xzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 644);

    auto g_0_xzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 646);

    auto g_0_xzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 647);

    auto g_0_xzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 648);

    auto g_0_xzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 649);

    auto g_0_xzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 650);

    auto g_0_xzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 652);

    auto g_0_xzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 653);

    auto g_0_xzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 654);

    auto g_0_xzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 655);

    auto g_0_xzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 656);

    auto g_0_xzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 657);

    auto g_0_xzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 659);

    auto g_0_xzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 660);

    auto g_0_xzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 661);

    auto g_0_xzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 662);

    auto g_0_xzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 663);

    auto g_0_xzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 664);

    auto g_0_xzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 665);

    auto g_0_xzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 666);

    auto g_0_xzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 667);

    auto g_0_xzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 668);

    auto g_0_xzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 669);

    auto g_0_xzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 670);

    auto g_0_xzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 671);

    auto g_0_xzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 672);

    auto g_0_xzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 673);

    auto g_0_xzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 674);

    auto g_0_yyyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 675);

    auto g_0_yyyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 676);

    auto g_0_yyyyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 677);

    auto g_0_yyyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 678);

    auto g_0_yyyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 679);

    auto g_0_yyyyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 680);

    auto g_0_yyyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 681);

    auto g_0_yyyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 682);

    auto g_0_yyyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 683);

    auto g_0_yyyyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 684);

    auto g_0_yyyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 685);

    auto g_0_yyyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 686);

    auto g_0_yyyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 687);

    auto g_0_yyyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 688);

    auto g_0_yyyyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 689);

    auto g_0_yyyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 690);

    auto g_0_yyyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 691);

    auto g_0_yyyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 692);

    auto g_0_yyyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 693);

    auto g_0_yyyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 694);

    auto g_0_yyyyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 695);

    auto g_0_yyyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 696);

    auto g_0_yyyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 697);

    auto g_0_yyyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 698);

    auto g_0_yyyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 699);

    auto g_0_yyyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 700);

    auto g_0_yyyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 701);

    auto g_0_yyyyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 702);

    auto g_0_yyyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 703);

    auto g_0_yyyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 704);

    auto g_0_yyyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 705);

    auto g_0_yyyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 706);

    auto g_0_yyyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 707);

    auto g_0_yyyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 708);

    auto g_0_yyyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 709);

    auto g_0_yyyyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 710);

    auto g_0_yyyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 711);

    auto g_0_yyyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 712);

    auto g_0_yyyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 713);

    auto g_0_yyyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 714);

    auto g_0_yyyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 715);

    auto g_0_yyyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 716);

    auto g_0_yyyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 717);

    auto g_0_yyyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 718);

    auto g_0_yyyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 719);

    auto g_0_yyyyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 721);

    auto g_0_yyyyz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 722);

    auto g_0_yyyyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 723);

    auto g_0_yyyyz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 724);

    auto g_0_yyyyz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 725);

    auto g_0_yyyyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 726);

    auto g_0_yyyyz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 727);

    auto g_0_yyyyz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 728);

    auto g_0_yyyyz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 729);

    auto g_0_yyyyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 730);

    auto g_0_yyyyz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 731);

    auto g_0_yyyyz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 732);

    auto g_0_yyyyz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 733);

    auto g_0_yyyyz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 734);

    auto g_0_yyyyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 735);

    auto g_0_yyyyz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 736);

    auto g_0_yyyyz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 737);

    auto g_0_yyyyz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 738);

    auto g_0_yyyyz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 739);

    auto g_0_yyyyz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 740);

    auto g_0_yyyyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 741);

    auto g_0_yyyyz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 742);

    auto g_0_yyyyz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 743);

    auto g_0_yyyyz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 744);

    auto g_0_yyyyz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 745);

    auto g_0_yyyyz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 746);

    auto g_0_yyyyz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 747);

    auto g_0_yyyyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 748);

    auto g_0_yyyyz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 749);

    auto g_0_yyyyz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 750);

    auto g_0_yyyyz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 751);

    auto g_0_yyyyz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 752);

    auto g_0_yyyyz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 753);

    auto g_0_yyyyz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 754);

    auto g_0_yyyyz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 755);

    auto g_0_yyyyz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 756);

    auto g_0_yyyyz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 757);

    auto g_0_yyyyz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 758);

    auto g_0_yyyyz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 759);

    auto g_0_yyyyz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 760);

    auto g_0_yyyyz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 761);

    auto g_0_yyyyz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 762);

    auto g_0_yyyyz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 763);

    auto g_0_yyyyz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 764);

    auto g_0_yyyzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 765);

    auto g_0_yyyzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 766);

    auto g_0_yyyzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 767);

    auto g_0_yyyzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 768);

    auto g_0_yyyzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 769);

    auto g_0_yyyzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 770);

    auto g_0_yyyzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 771);

    auto g_0_yyyzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 772);

    auto g_0_yyyzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 773);

    auto g_0_yyyzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 774);

    auto g_0_yyyzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 775);

    auto g_0_yyyzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 776);

    auto g_0_yyyzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 777);

    auto g_0_yyyzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 778);

    auto g_0_yyyzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 779);

    auto g_0_yyyzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 780);

    auto g_0_yyyzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 781);

    auto g_0_yyyzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 782);

    auto g_0_yyyzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 783);

    auto g_0_yyyzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 784);

    auto g_0_yyyzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 785);

    auto g_0_yyyzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 786);

    auto g_0_yyyzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 787);

    auto g_0_yyyzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 788);

    auto g_0_yyyzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 789);

    auto g_0_yyyzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 790);

    auto g_0_yyyzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 791);

    auto g_0_yyyzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 792);

    auto g_0_yyyzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 793);

    auto g_0_yyyzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 794);

    auto g_0_yyyzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 795);

    auto g_0_yyyzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 796);

    auto g_0_yyyzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 797);

    auto g_0_yyyzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 798);

    auto g_0_yyyzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 799);

    auto g_0_yyyzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 800);

    auto g_0_yyyzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 801);

    auto g_0_yyyzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 802);

    auto g_0_yyyzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 803);

    auto g_0_yyyzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 804);

    auto g_0_yyyzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 805);

    auto g_0_yyyzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 806);

    auto g_0_yyyzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 807);

    auto g_0_yyyzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 808);

    auto g_0_yyyzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 809);

    auto g_0_yyzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 810);

    auto g_0_yyzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 811);

    auto g_0_yyzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 812);

    auto g_0_yyzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 813);

    auto g_0_yyzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 814);

    auto g_0_yyzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 815);

    auto g_0_yyzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 816);

    auto g_0_yyzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 817);

    auto g_0_yyzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 818);

    auto g_0_yyzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 819);

    auto g_0_yyzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 820);

    auto g_0_yyzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 821);

    auto g_0_yyzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 822);

    auto g_0_yyzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 823);

    auto g_0_yyzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 824);

    auto g_0_yyzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 825);

    auto g_0_yyzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 826);

    auto g_0_yyzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 827);

    auto g_0_yyzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 828);

    auto g_0_yyzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 829);

    auto g_0_yyzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 830);

    auto g_0_yyzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 831);

    auto g_0_yyzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 832);

    auto g_0_yyzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 833);

    auto g_0_yyzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 834);

    auto g_0_yyzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 835);

    auto g_0_yyzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 836);

    auto g_0_yyzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 837);

    auto g_0_yyzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 838);

    auto g_0_yyzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 839);

    auto g_0_yyzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 840);

    auto g_0_yyzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 841);

    auto g_0_yyzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 842);

    auto g_0_yyzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 843);

    auto g_0_yyzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 844);

    auto g_0_yyzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 845);

    auto g_0_yyzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 846);

    auto g_0_yyzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 847);

    auto g_0_yyzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 848);

    auto g_0_yyzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 849);

    auto g_0_yyzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 850);

    auto g_0_yyzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 851);

    auto g_0_yyzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 852);

    auto g_0_yyzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 853);

    auto g_0_yyzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 854);

    auto g_0_yzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 855);

    auto g_0_yzzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 856);

    auto g_0_yzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 857);

    auto g_0_yzzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 858);

    auto g_0_yzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 859);

    auto g_0_yzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 860);

    auto g_0_yzzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 861);

    auto g_0_yzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 862);

    auto g_0_yzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 863);

    auto g_0_yzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 864);

    auto g_0_yzzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 865);

    auto g_0_yzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 866);

    auto g_0_yzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 867);

    auto g_0_yzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 868);

    auto g_0_yzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 869);

    auto g_0_yzzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 870);

    auto g_0_yzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 871);

    auto g_0_yzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 872);

    auto g_0_yzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 873);

    auto g_0_yzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 874);

    auto g_0_yzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 875);

    auto g_0_yzzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 876);

    auto g_0_yzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 877);

    auto g_0_yzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 878);

    auto g_0_yzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 879);

    auto g_0_yzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 880);

    auto g_0_yzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 881);

    auto g_0_yzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 882);

    auto g_0_yzzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 883);

    auto g_0_yzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 884);

    auto g_0_yzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 885);

    auto g_0_yzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 886);

    auto g_0_yzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 887);

    auto g_0_yzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 888);

    auto g_0_yzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 889);

    auto g_0_yzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 890);

    auto g_0_yzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 891);

    auto g_0_yzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 892);

    auto g_0_yzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 893);

    auto g_0_yzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 894);

    auto g_0_yzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 895);

    auto g_0_yzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 896);

    auto g_0_yzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 897);

    auto g_0_yzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 898);

    auto g_0_yzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 899);

    auto g_0_zzzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_shsl + 900);

    auto g_0_zzzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_shsl + 901);

    auto g_0_zzzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_shsl + 902);

    auto g_0_zzzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_shsl + 903);

    auto g_0_zzzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_shsl + 904);

    auto g_0_zzzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_shsl + 905);

    auto g_0_zzzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_shsl + 906);

    auto g_0_zzzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_shsl + 907);

    auto g_0_zzzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_shsl + 908);

    auto g_0_zzzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_shsl + 909);

    auto g_0_zzzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_shsl + 910);

    auto g_0_zzzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_shsl + 911);

    auto g_0_zzzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_shsl + 912);

    auto g_0_zzzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_shsl + 913);

    auto g_0_zzzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_shsl + 914);

    auto g_0_zzzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 915);

    auto g_0_zzzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 916);

    auto g_0_zzzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 917);

    auto g_0_zzzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 918);

    auto g_0_zzzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 919);

    auto g_0_zzzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 920);

    auto g_0_zzzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 921);

    auto g_0_zzzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 922);

    auto g_0_zzzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 923);

    auto g_0_zzzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 924);

    auto g_0_zzzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 925);

    auto g_0_zzzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 926);

    auto g_0_zzzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 927);

    auto g_0_zzzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 928);

    auto g_0_zzzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 929);

    auto g_0_zzzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 930);

    auto g_0_zzzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 931);

    auto g_0_zzzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 932);

    auto g_0_zzzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 933);

    auto g_0_zzzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 934);

    auto g_0_zzzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 935);

    auto g_0_zzzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_shsl + 936);

    auto g_0_zzzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_shsl + 937);

    auto g_0_zzzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_shsl + 938);

    auto g_0_zzzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_shsl + 939);

    auto g_0_zzzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_shsl + 940);

    auto g_0_zzzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 941);

    auto g_0_zzzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 942);

    auto g_0_zzzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 943);

    auto g_0_zzzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_shsl + 944);

    /// Set up 0-45 components of targeted buffer : SISL

    auto g_0_xxxxxx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl);

    auto g_0_xxxxxx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 1);

    auto g_0_xxxxxx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 2);

    auto g_0_xxxxxx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 3);

    auto g_0_xxxxxx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 4);

    auto g_0_xxxxxx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 5);

    auto g_0_xxxxxx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 6);

    auto g_0_xxxxxx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 7);

    auto g_0_xxxxxx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 8);

    auto g_0_xxxxxx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 9);

    auto g_0_xxxxxx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 10);

    auto g_0_xxxxxx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 11);

    auto g_0_xxxxxx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 12);

    auto g_0_xxxxxx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 13);

    auto g_0_xxxxxx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 14);

    auto g_0_xxxxxx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 15);

    auto g_0_xxxxxx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 16);

    auto g_0_xxxxxx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 17);

    auto g_0_xxxxxx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 18);

    auto g_0_xxxxxx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 19);

    auto g_0_xxxxxx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 20);

    auto g_0_xxxxxx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 21);

    auto g_0_xxxxxx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 22);

    auto g_0_xxxxxx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 23);

    auto g_0_xxxxxx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 24);

    auto g_0_xxxxxx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 25);

    auto g_0_xxxxxx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 26);

    auto g_0_xxxxxx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 27);

    auto g_0_xxxxxx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 28);

    auto g_0_xxxxxx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 29);

    auto g_0_xxxxxx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 30);

    auto g_0_xxxxxx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 31);

    auto g_0_xxxxxx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 32);

    auto g_0_xxxxxx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 33);

    auto g_0_xxxxxx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 34);

    auto g_0_xxxxxx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 35);

    auto g_0_xxxxxx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 36);

    auto g_0_xxxxxx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 37);

    auto g_0_xxxxxx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 38);

    auto g_0_xxxxxx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 39);

    auto g_0_xxxxxx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 40);

    auto g_0_xxxxxx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 41);

    auto g_0_xxxxxx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 42);

    auto g_0_xxxxxx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 43);

    auto g_0_xxxxxx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 44);

#pragma omp simd aligned(g_0_xxxx_0_xxxxxxxx_0,       \
                             g_0_xxxx_0_xxxxxxxx_1,   \
                             g_0_xxxx_0_xxxxxxxy_0,   \
                             g_0_xxxx_0_xxxxxxxy_1,   \
                             g_0_xxxx_0_xxxxxxxz_0,   \
                             g_0_xxxx_0_xxxxxxxz_1,   \
                             g_0_xxxx_0_xxxxxxyy_0,   \
                             g_0_xxxx_0_xxxxxxyy_1,   \
                             g_0_xxxx_0_xxxxxxyz_0,   \
                             g_0_xxxx_0_xxxxxxyz_1,   \
                             g_0_xxxx_0_xxxxxxzz_0,   \
                             g_0_xxxx_0_xxxxxxzz_1,   \
                             g_0_xxxx_0_xxxxxyyy_0,   \
                             g_0_xxxx_0_xxxxxyyy_1,   \
                             g_0_xxxx_0_xxxxxyyz_0,   \
                             g_0_xxxx_0_xxxxxyyz_1,   \
                             g_0_xxxx_0_xxxxxyzz_0,   \
                             g_0_xxxx_0_xxxxxyzz_1,   \
                             g_0_xxxx_0_xxxxxzzz_0,   \
                             g_0_xxxx_0_xxxxxzzz_1,   \
                             g_0_xxxx_0_xxxxyyyy_0,   \
                             g_0_xxxx_0_xxxxyyyy_1,   \
                             g_0_xxxx_0_xxxxyyyz_0,   \
                             g_0_xxxx_0_xxxxyyyz_1,   \
                             g_0_xxxx_0_xxxxyyzz_0,   \
                             g_0_xxxx_0_xxxxyyzz_1,   \
                             g_0_xxxx_0_xxxxyzzz_0,   \
                             g_0_xxxx_0_xxxxyzzz_1,   \
                             g_0_xxxx_0_xxxxzzzz_0,   \
                             g_0_xxxx_0_xxxxzzzz_1,   \
                             g_0_xxxx_0_xxxyyyyy_0,   \
                             g_0_xxxx_0_xxxyyyyy_1,   \
                             g_0_xxxx_0_xxxyyyyz_0,   \
                             g_0_xxxx_0_xxxyyyyz_1,   \
                             g_0_xxxx_0_xxxyyyzz_0,   \
                             g_0_xxxx_0_xxxyyyzz_1,   \
                             g_0_xxxx_0_xxxyyzzz_0,   \
                             g_0_xxxx_0_xxxyyzzz_1,   \
                             g_0_xxxx_0_xxxyzzzz_0,   \
                             g_0_xxxx_0_xxxyzzzz_1,   \
                             g_0_xxxx_0_xxxzzzzz_0,   \
                             g_0_xxxx_0_xxxzzzzz_1,   \
                             g_0_xxxx_0_xxyyyyyy_0,   \
                             g_0_xxxx_0_xxyyyyyy_1,   \
                             g_0_xxxx_0_xxyyyyyz_0,   \
                             g_0_xxxx_0_xxyyyyyz_1,   \
                             g_0_xxxx_0_xxyyyyzz_0,   \
                             g_0_xxxx_0_xxyyyyzz_1,   \
                             g_0_xxxx_0_xxyyyzzz_0,   \
                             g_0_xxxx_0_xxyyyzzz_1,   \
                             g_0_xxxx_0_xxyyzzzz_0,   \
                             g_0_xxxx_0_xxyyzzzz_1,   \
                             g_0_xxxx_0_xxyzzzzz_0,   \
                             g_0_xxxx_0_xxyzzzzz_1,   \
                             g_0_xxxx_0_xxzzzzzz_0,   \
                             g_0_xxxx_0_xxzzzzzz_1,   \
                             g_0_xxxx_0_xyyyyyyy_0,   \
                             g_0_xxxx_0_xyyyyyyy_1,   \
                             g_0_xxxx_0_xyyyyyyz_0,   \
                             g_0_xxxx_0_xyyyyyyz_1,   \
                             g_0_xxxx_0_xyyyyyzz_0,   \
                             g_0_xxxx_0_xyyyyyzz_1,   \
                             g_0_xxxx_0_xyyyyzzz_0,   \
                             g_0_xxxx_0_xyyyyzzz_1,   \
                             g_0_xxxx_0_xyyyzzzz_0,   \
                             g_0_xxxx_0_xyyyzzzz_1,   \
                             g_0_xxxx_0_xyyzzzzz_0,   \
                             g_0_xxxx_0_xyyzzzzz_1,   \
                             g_0_xxxx_0_xyzzzzzz_0,   \
                             g_0_xxxx_0_xyzzzzzz_1,   \
                             g_0_xxxx_0_xzzzzzzz_0,   \
                             g_0_xxxx_0_xzzzzzzz_1,   \
                             g_0_xxxx_0_yyyyyyyy_0,   \
                             g_0_xxxx_0_yyyyyyyy_1,   \
                             g_0_xxxx_0_yyyyyyyz_0,   \
                             g_0_xxxx_0_yyyyyyyz_1,   \
                             g_0_xxxx_0_yyyyyyzz_0,   \
                             g_0_xxxx_0_yyyyyyzz_1,   \
                             g_0_xxxx_0_yyyyyzzz_0,   \
                             g_0_xxxx_0_yyyyyzzz_1,   \
                             g_0_xxxx_0_yyyyzzzz_0,   \
                             g_0_xxxx_0_yyyyzzzz_1,   \
                             g_0_xxxx_0_yyyzzzzz_0,   \
                             g_0_xxxx_0_yyyzzzzz_1,   \
                             g_0_xxxx_0_yyzzzzzz_0,   \
                             g_0_xxxx_0_yyzzzzzz_1,   \
                             g_0_xxxx_0_yzzzzzzz_0,   \
                             g_0_xxxx_0_yzzzzzzz_1,   \
                             g_0_xxxx_0_zzzzzzzz_0,   \
                             g_0_xxxx_0_zzzzzzzz_1,   \
                             g_0_xxxxx_0_xxxxxxx_1,   \
                             g_0_xxxxx_0_xxxxxxxx_0,  \
                             g_0_xxxxx_0_xxxxxxxx_1,  \
                             g_0_xxxxx_0_xxxxxxxy_0,  \
                             g_0_xxxxx_0_xxxxxxxy_1,  \
                             g_0_xxxxx_0_xxxxxxxz_0,  \
                             g_0_xxxxx_0_xxxxxxxz_1,  \
                             g_0_xxxxx_0_xxxxxxy_1,   \
                             g_0_xxxxx_0_xxxxxxyy_0,  \
                             g_0_xxxxx_0_xxxxxxyy_1,  \
                             g_0_xxxxx_0_xxxxxxyz_0,  \
                             g_0_xxxxx_0_xxxxxxyz_1,  \
                             g_0_xxxxx_0_xxxxxxz_1,   \
                             g_0_xxxxx_0_xxxxxxzz_0,  \
                             g_0_xxxxx_0_xxxxxxzz_1,  \
                             g_0_xxxxx_0_xxxxxyy_1,   \
                             g_0_xxxxx_0_xxxxxyyy_0,  \
                             g_0_xxxxx_0_xxxxxyyy_1,  \
                             g_0_xxxxx_0_xxxxxyyz_0,  \
                             g_0_xxxxx_0_xxxxxyyz_1,  \
                             g_0_xxxxx_0_xxxxxyz_1,   \
                             g_0_xxxxx_0_xxxxxyzz_0,  \
                             g_0_xxxxx_0_xxxxxyzz_1,  \
                             g_0_xxxxx_0_xxxxxzz_1,   \
                             g_0_xxxxx_0_xxxxxzzz_0,  \
                             g_0_xxxxx_0_xxxxxzzz_1,  \
                             g_0_xxxxx_0_xxxxyyy_1,   \
                             g_0_xxxxx_0_xxxxyyyy_0,  \
                             g_0_xxxxx_0_xxxxyyyy_1,  \
                             g_0_xxxxx_0_xxxxyyyz_0,  \
                             g_0_xxxxx_0_xxxxyyyz_1,  \
                             g_0_xxxxx_0_xxxxyyz_1,   \
                             g_0_xxxxx_0_xxxxyyzz_0,  \
                             g_0_xxxxx_0_xxxxyyzz_1,  \
                             g_0_xxxxx_0_xxxxyzz_1,   \
                             g_0_xxxxx_0_xxxxyzzz_0,  \
                             g_0_xxxxx_0_xxxxyzzz_1,  \
                             g_0_xxxxx_0_xxxxzzz_1,   \
                             g_0_xxxxx_0_xxxxzzzz_0,  \
                             g_0_xxxxx_0_xxxxzzzz_1,  \
                             g_0_xxxxx_0_xxxyyyy_1,   \
                             g_0_xxxxx_0_xxxyyyyy_0,  \
                             g_0_xxxxx_0_xxxyyyyy_1,  \
                             g_0_xxxxx_0_xxxyyyyz_0,  \
                             g_0_xxxxx_0_xxxyyyyz_1,  \
                             g_0_xxxxx_0_xxxyyyz_1,   \
                             g_0_xxxxx_0_xxxyyyzz_0,  \
                             g_0_xxxxx_0_xxxyyyzz_1,  \
                             g_0_xxxxx_0_xxxyyzz_1,   \
                             g_0_xxxxx_0_xxxyyzzz_0,  \
                             g_0_xxxxx_0_xxxyyzzz_1,  \
                             g_0_xxxxx_0_xxxyzzz_1,   \
                             g_0_xxxxx_0_xxxyzzzz_0,  \
                             g_0_xxxxx_0_xxxyzzzz_1,  \
                             g_0_xxxxx_0_xxxzzzz_1,   \
                             g_0_xxxxx_0_xxxzzzzz_0,  \
                             g_0_xxxxx_0_xxxzzzzz_1,  \
                             g_0_xxxxx_0_xxyyyyy_1,   \
                             g_0_xxxxx_0_xxyyyyyy_0,  \
                             g_0_xxxxx_0_xxyyyyyy_1,  \
                             g_0_xxxxx_0_xxyyyyyz_0,  \
                             g_0_xxxxx_0_xxyyyyyz_1,  \
                             g_0_xxxxx_0_xxyyyyz_1,   \
                             g_0_xxxxx_0_xxyyyyzz_0,  \
                             g_0_xxxxx_0_xxyyyyzz_1,  \
                             g_0_xxxxx_0_xxyyyzz_1,   \
                             g_0_xxxxx_0_xxyyyzzz_0,  \
                             g_0_xxxxx_0_xxyyyzzz_1,  \
                             g_0_xxxxx_0_xxyyzzz_1,   \
                             g_0_xxxxx_0_xxyyzzzz_0,  \
                             g_0_xxxxx_0_xxyyzzzz_1,  \
                             g_0_xxxxx_0_xxyzzzz_1,   \
                             g_0_xxxxx_0_xxyzzzzz_0,  \
                             g_0_xxxxx_0_xxyzzzzz_1,  \
                             g_0_xxxxx_0_xxzzzzz_1,   \
                             g_0_xxxxx_0_xxzzzzzz_0,  \
                             g_0_xxxxx_0_xxzzzzzz_1,  \
                             g_0_xxxxx_0_xyyyyyy_1,   \
                             g_0_xxxxx_0_xyyyyyyy_0,  \
                             g_0_xxxxx_0_xyyyyyyy_1,  \
                             g_0_xxxxx_0_xyyyyyyz_0,  \
                             g_0_xxxxx_0_xyyyyyyz_1,  \
                             g_0_xxxxx_0_xyyyyyz_1,   \
                             g_0_xxxxx_0_xyyyyyzz_0,  \
                             g_0_xxxxx_0_xyyyyyzz_1,  \
                             g_0_xxxxx_0_xyyyyzz_1,   \
                             g_0_xxxxx_0_xyyyyzzz_0,  \
                             g_0_xxxxx_0_xyyyyzzz_1,  \
                             g_0_xxxxx_0_xyyyzzz_1,   \
                             g_0_xxxxx_0_xyyyzzzz_0,  \
                             g_0_xxxxx_0_xyyyzzzz_1,  \
                             g_0_xxxxx_0_xyyzzzz_1,   \
                             g_0_xxxxx_0_xyyzzzzz_0,  \
                             g_0_xxxxx_0_xyyzzzzz_1,  \
                             g_0_xxxxx_0_xyzzzzz_1,   \
                             g_0_xxxxx_0_xyzzzzzz_0,  \
                             g_0_xxxxx_0_xyzzzzzz_1,  \
                             g_0_xxxxx_0_xzzzzzz_1,   \
                             g_0_xxxxx_0_xzzzzzzz_0,  \
                             g_0_xxxxx_0_xzzzzzzz_1,  \
                             g_0_xxxxx_0_yyyyyyy_1,   \
                             g_0_xxxxx_0_yyyyyyyy_0,  \
                             g_0_xxxxx_0_yyyyyyyy_1,  \
                             g_0_xxxxx_0_yyyyyyyz_0,  \
                             g_0_xxxxx_0_yyyyyyyz_1,  \
                             g_0_xxxxx_0_yyyyyyz_1,   \
                             g_0_xxxxx_0_yyyyyyzz_0,  \
                             g_0_xxxxx_0_yyyyyyzz_1,  \
                             g_0_xxxxx_0_yyyyyzz_1,   \
                             g_0_xxxxx_0_yyyyyzzz_0,  \
                             g_0_xxxxx_0_yyyyyzzz_1,  \
                             g_0_xxxxx_0_yyyyzzz_1,   \
                             g_0_xxxxx_0_yyyyzzzz_0,  \
                             g_0_xxxxx_0_yyyyzzzz_1,  \
                             g_0_xxxxx_0_yyyzzzz_1,   \
                             g_0_xxxxx_0_yyyzzzzz_0,  \
                             g_0_xxxxx_0_yyyzzzzz_1,  \
                             g_0_xxxxx_0_yyzzzzz_1,   \
                             g_0_xxxxx_0_yyzzzzzz_0,  \
                             g_0_xxxxx_0_yyzzzzzz_1,  \
                             g_0_xxxxx_0_yzzzzzz_1,   \
                             g_0_xxxxx_0_yzzzzzzz_0,  \
                             g_0_xxxxx_0_yzzzzzzz_1,  \
                             g_0_xxxxx_0_zzzzzzz_1,   \
                             g_0_xxxxx_0_zzzzzzzz_0,  \
                             g_0_xxxxx_0_zzzzzzzz_1,  \
                             g_0_xxxxxx_0_xxxxxxxx_0, \
                             g_0_xxxxxx_0_xxxxxxxy_0, \
                             g_0_xxxxxx_0_xxxxxxxz_0, \
                             g_0_xxxxxx_0_xxxxxxyy_0, \
                             g_0_xxxxxx_0_xxxxxxyz_0, \
                             g_0_xxxxxx_0_xxxxxxzz_0, \
                             g_0_xxxxxx_0_xxxxxyyy_0, \
                             g_0_xxxxxx_0_xxxxxyyz_0, \
                             g_0_xxxxxx_0_xxxxxyzz_0, \
                             g_0_xxxxxx_0_xxxxxzzz_0, \
                             g_0_xxxxxx_0_xxxxyyyy_0, \
                             g_0_xxxxxx_0_xxxxyyyz_0, \
                             g_0_xxxxxx_0_xxxxyyzz_0, \
                             g_0_xxxxxx_0_xxxxyzzz_0, \
                             g_0_xxxxxx_0_xxxxzzzz_0, \
                             g_0_xxxxxx_0_xxxyyyyy_0, \
                             g_0_xxxxxx_0_xxxyyyyz_0, \
                             g_0_xxxxxx_0_xxxyyyzz_0, \
                             g_0_xxxxxx_0_xxxyyzzz_0, \
                             g_0_xxxxxx_0_xxxyzzzz_0, \
                             g_0_xxxxxx_0_xxxzzzzz_0, \
                             g_0_xxxxxx_0_xxyyyyyy_0, \
                             g_0_xxxxxx_0_xxyyyyyz_0, \
                             g_0_xxxxxx_0_xxyyyyzz_0, \
                             g_0_xxxxxx_0_xxyyyzzz_0, \
                             g_0_xxxxxx_0_xxyyzzzz_0, \
                             g_0_xxxxxx_0_xxyzzzzz_0, \
                             g_0_xxxxxx_0_xxzzzzzz_0, \
                             g_0_xxxxxx_0_xyyyyyyy_0, \
                             g_0_xxxxxx_0_xyyyyyyz_0, \
                             g_0_xxxxxx_0_xyyyyyzz_0, \
                             g_0_xxxxxx_0_xyyyyzzz_0, \
                             g_0_xxxxxx_0_xyyyzzzz_0, \
                             g_0_xxxxxx_0_xyyzzzzz_0, \
                             g_0_xxxxxx_0_xyzzzzzz_0, \
                             g_0_xxxxxx_0_xzzzzzzz_0, \
                             g_0_xxxxxx_0_yyyyyyyy_0, \
                             g_0_xxxxxx_0_yyyyyyyz_0, \
                             g_0_xxxxxx_0_yyyyyyzz_0, \
                             g_0_xxxxxx_0_yyyyyzzz_0, \
                             g_0_xxxxxx_0_yyyyzzzz_0, \
                             g_0_xxxxxx_0_yyyzzzzz_0, \
                             g_0_xxxxxx_0_yyzzzzzz_0, \
                             g_0_xxxxxx_0_yzzzzzzz_0, \
                             g_0_xxxxxx_0_zzzzzzzz_0, \
                             wp_x,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxx_0_xxxxxxxx_0[i] = 5.0 * g_0_xxxx_0_xxxxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxxxx_1[i] * fti_ab_0 +
                                     8.0 * g_0_xxxxx_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxxx_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxxxy_0[i] = 5.0 * g_0_xxxx_0_xxxxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxxxy_1[i] * fti_ab_0 +
                                     7.0 * g_0_xxxxx_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxxy_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxxxz_0[i] = 5.0 * g_0_xxxx_0_xxxxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxxxz_1[i] * fti_ab_0 +
                                     7.0 * g_0_xxxxx_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxxz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxxyy_0[i] = 5.0 * g_0_xxxx_0_xxxxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxxyy_1[i] * fti_ab_0 +
                                     6.0 * g_0_xxxxx_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxyy_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxxyz_0[i] = 5.0 * g_0_xxxx_0_xxxxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxxyz_1[i] * fti_ab_0 +
                                     6.0 * g_0_xxxxx_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxyz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxxzz_0[i] = 5.0 * g_0_xxxx_0_xxxxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxxzz_1[i] * fti_ab_0 +
                                     6.0 * g_0_xxxxx_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxyyy_0[i] = 5.0 * g_0_xxxx_0_xxxxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxyyy_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxxx_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyyy_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxyyz_0[i] = 5.0 * g_0_xxxx_0_xxxxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxxx_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyyz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxyzz_0[i] = 5.0 * g_0_xxxx_0_xxxxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxyzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxxx_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxzzz_0[i] = 5.0 * g_0_xxxx_0_xxxxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxxx_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxzzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxyyyy_0[i] = 5.0 * g_0_xxxx_0_xxxxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxyyyy_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxx_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyyy_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxyyyz_0[i] = 5.0 * g_0_xxxx_0_xxxxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxx_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyyz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxyyzz_0[i] = 5.0 * g_0_xxxx_0_xxxxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxx_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxyzzz_0[i] = 5.0 * g_0_xxxx_0_xxxxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxyzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxx_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyzzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxzzzz_0[i] = 5.0 * g_0_xxxx_0_xxxxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxx_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxzzzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyyyyy_0[i] = 5.0 * g_0_xxxx_0_xxxyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyyyyy_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxx_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyyy_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyyyyz_0[i] = 5.0 * g_0_xxxx_0_xxxyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxx_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyyz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyyyzz_0[i] = 5.0 * g_0_xxxx_0_xxxyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxx_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyyzzz_0[i] = 5.0 * g_0_xxxx_0_xxxyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxx_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyzzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyzzzz_0[i] = 5.0 * g_0_xxxx_0_xxxyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxx_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyzzzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxzzzzz_0[i] = 5.0 * g_0_xxxx_0_xxxzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxzzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxx_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxzzzzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyyyyy_0[i] = 5.0 * g_0_xxxx_0_xxyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyyyyy_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxx_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyyy_0[i] * pb_x +
                                     g_0_xxxxx_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyyyyz_0[i] = 5.0 * g_0_xxxx_0_xxyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxx_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyyz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyyyzz_0[i] = 5.0 * g_0_xxxx_0_xxyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxx_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyyzzz_0[i] = 5.0 * g_0_xxxx_0_xxyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxx_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyzzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyzzzz_0[i] = 5.0 * g_0_xxxx_0_xxyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxx_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyzzzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyzzzzz_0[i] = 5.0 * g_0_xxxx_0_xxyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxx_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyzzzzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxzzzzzz_0[i] = 5.0 * g_0_xxxx_0_xxzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxzzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxx_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxzzzzzz_0[i] * pb_x +
                                     g_0_xxxxx_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyyyyy_0[i] = 5.0 * g_0_xxxx_0_xyyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyyy_0[i] * pb_x + g_0_xxxxx_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyyyyz_0[i] = 5.0 * g_0_xxxx_0_xyyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyyz_0[i] * pb_x + g_0_xxxxx_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyyyzz_0[i] = 5.0 * g_0_xxxx_0_xyyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyzz_0[i] * pb_x + g_0_xxxxx_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyyzzz_0[i] = 5.0 * g_0_xxxx_0_xyyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyzzz_0[i] * pb_x + g_0_xxxxx_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyzzzz_0[i] = 5.0 * g_0_xxxx_0_xyyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyzzzz_0[i] * pb_x + g_0_xxxxx_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyzzzzz_0[i] = 5.0 * g_0_xxxx_0_xyyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyzzzzz_0[i] * pb_x + g_0_xxxxx_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyzzzzzz_0[i] = 5.0 * g_0_xxxx_0_xyzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzzzzzz_0[i] * pb_x + g_0_xxxxx_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xzzzzzzz_0[i] = 5.0 * g_0_xxxx_0_xzzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xzzzzzzz_0[i] * pb_x + g_0_xxxxx_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyyyyy_0[i] = 5.0 * g_0_xxxx_0_yyyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyyyyyyy_0[i] * pb_x + g_0_xxxxx_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyyyyz_0[i] = 5.0 * g_0_xxxx_0_yyyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyyyyyyz_0[i] * pb_x + g_0_xxxxx_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyyyzz_0[i] = 5.0 * g_0_xxxx_0_yyyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyyyyyzz_0[i] * pb_x + g_0_xxxxx_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyyzzz_0[i] = 5.0 * g_0_xxxx_0_yyyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyyyyzzz_0[i] * pb_x + g_0_xxxxx_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyzzzz_0[i] = 5.0 * g_0_xxxx_0_yyyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyyyzzzz_0[i] * pb_x + g_0_xxxxx_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyzzzzz_0[i] = 5.0 * g_0_xxxx_0_yyyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyyzzzzz_0[i] * pb_x + g_0_xxxxx_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyzzzzzz_0[i] = 5.0 * g_0_xxxx_0_yyzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yyzzzzzz_0[i] * pb_x + g_0_xxxxx_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yzzzzzzz_0[i] = 5.0 * g_0_xxxx_0_yzzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_yzzzzzzz_0[i] * pb_x + g_0_xxxxx_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_zzzzzzzz_0[i] = 5.0 * g_0_xxxx_0_zzzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_zzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxx_0_zzzzzzzz_0[i] * pb_x + g_0_xxxxx_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 45-90 components of targeted buffer : SISL

    auto g_0_xxxxxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 45);

    auto g_0_xxxxxy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 46);

    auto g_0_xxxxxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 47);

    auto g_0_xxxxxy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 48);

    auto g_0_xxxxxy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 49);

    auto g_0_xxxxxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 50);

    auto g_0_xxxxxy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 51);

    auto g_0_xxxxxy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 52);

    auto g_0_xxxxxy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 53);

    auto g_0_xxxxxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 54);

    auto g_0_xxxxxy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 55);

    auto g_0_xxxxxy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 56);

    auto g_0_xxxxxy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 57);

    auto g_0_xxxxxy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 58);

    auto g_0_xxxxxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 59);

    auto g_0_xxxxxy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 60);

    auto g_0_xxxxxy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 61);

    auto g_0_xxxxxy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 62);

    auto g_0_xxxxxy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 63);

    auto g_0_xxxxxy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 64);

    auto g_0_xxxxxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 65);

    auto g_0_xxxxxy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 66);

    auto g_0_xxxxxy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 67);

    auto g_0_xxxxxy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 68);

    auto g_0_xxxxxy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 69);

    auto g_0_xxxxxy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 70);

    auto g_0_xxxxxy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 71);

    auto g_0_xxxxxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 72);

    auto g_0_xxxxxy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 73);

    auto g_0_xxxxxy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 74);

    auto g_0_xxxxxy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 75);

    auto g_0_xxxxxy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 76);

    auto g_0_xxxxxy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 77);

    auto g_0_xxxxxy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 78);

    auto g_0_xxxxxy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 79);

    auto g_0_xxxxxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 80);

    auto g_0_xxxxxy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 81);

    auto g_0_xxxxxy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 82);

    auto g_0_xxxxxy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 83);

    auto g_0_xxxxxy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 84);

    auto g_0_xxxxxy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 85);

    auto g_0_xxxxxy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 86);

    auto g_0_xxxxxy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 87);

    auto g_0_xxxxxy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 88);

    auto g_0_xxxxxy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 89);

#pragma omp simd aligned(g_0_xxxxx_0_xxxxxxx_1,       \
                             g_0_xxxxx_0_xxxxxxxx_0,  \
                             g_0_xxxxx_0_xxxxxxxx_1,  \
                             g_0_xxxxx_0_xxxxxxxy_0,  \
                             g_0_xxxxx_0_xxxxxxxy_1,  \
                             g_0_xxxxx_0_xxxxxxxz_0,  \
                             g_0_xxxxx_0_xxxxxxxz_1,  \
                             g_0_xxxxx_0_xxxxxxy_1,   \
                             g_0_xxxxx_0_xxxxxxyy_0,  \
                             g_0_xxxxx_0_xxxxxxyy_1,  \
                             g_0_xxxxx_0_xxxxxxyz_0,  \
                             g_0_xxxxx_0_xxxxxxyz_1,  \
                             g_0_xxxxx_0_xxxxxxz_1,   \
                             g_0_xxxxx_0_xxxxxxzz_0,  \
                             g_0_xxxxx_0_xxxxxxzz_1,  \
                             g_0_xxxxx_0_xxxxxyy_1,   \
                             g_0_xxxxx_0_xxxxxyyy_0,  \
                             g_0_xxxxx_0_xxxxxyyy_1,  \
                             g_0_xxxxx_0_xxxxxyyz_0,  \
                             g_0_xxxxx_0_xxxxxyyz_1,  \
                             g_0_xxxxx_0_xxxxxyz_1,   \
                             g_0_xxxxx_0_xxxxxyzz_0,  \
                             g_0_xxxxx_0_xxxxxyzz_1,  \
                             g_0_xxxxx_0_xxxxxzz_1,   \
                             g_0_xxxxx_0_xxxxxzzz_0,  \
                             g_0_xxxxx_0_xxxxxzzz_1,  \
                             g_0_xxxxx_0_xxxxyyy_1,   \
                             g_0_xxxxx_0_xxxxyyyy_0,  \
                             g_0_xxxxx_0_xxxxyyyy_1,  \
                             g_0_xxxxx_0_xxxxyyyz_0,  \
                             g_0_xxxxx_0_xxxxyyyz_1,  \
                             g_0_xxxxx_0_xxxxyyz_1,   \
                             g_0_xxxxx_0_xxxxyyzz_0,  \
                             g_0_xxxxx_0_xxxxyyzz_1,  \
                             g_0_xxxxx_0_xxxxyzz_1,   \
                             g_0_xxxxx_0_xxxxyzzz_0,  \
                             g_0_xxxxx_0_xxxxyzzz_1,  \
                             g_0_xxxxx_0_xxxxzzz_1,   \
                             g_0_xxxxx_0_xxxxzzzz_0,  \
                             g_0_xxxxx_0_xxxxzzzz_1,  \
                             g_0_xxxxx_0_xxxyyyy_1,   \
                             g_0_xxxxx_0_xxxyyyyy_0,  \
                             g_0_xxxxx_0_xxxyyyyy_1,  \
                             g_0_xxxxx_0_xxxyyyyz_0,  \
                             g_0_xxxxx_0_xxxyyyyz_1,  \
                             g_0_xxxxx_0_xxxyyyz_1,   \
                             g_0_xxxxx_0_xxxyyyzz_0,  \
                             g_0_xxxxx_0_xxxyyyzz_1,  \
                             g_0_xxxxx_0_xxxyyzz_1,   \
                             g_0_xxxxx_0_xxxyyzzz_0,  \
                             g_0_xxxxx_0_xxxyyzzz_1,  \
                             g_0_xxxxx_0_xxxyzzz_1,   \
                             g_0_xxxxx_0_xxxyzzzz_0,  \
                             g_0_xxxxx_0_xxxyzzzz_1,  \
                             g_0_xxxxx_0_xxxzzzz_1,   \
                             g_0_xxxxx_0_xxxzzzzz_0,  \
                             g_0_xxxxx_0_xxxzzzzz_1,  \
                             g_0_xxxxx_0_xxyyyyy_1,   \
                             g_0_xxxxx_0_xxyyyyyy_0,  \
                             g_0_xxxxx_0_xxyyyyyy_1,  \
                             g_0_xxxxx_0_xxyyyyyz_0,  \
                             g_0_xxxxx_0_xxyyyyyz_1,  \
                             g_0_xxxxx_0_xxyyyyz_1,   \
                             g_0_xxxxx_0_xxyyyyzz_0,  \
                             g_0_xxxxx_0_xxyyyyzz_1,  \
                             g_0_xxxxx_0_xxyyyzz_1,   \
                             g_0_xxxxx_0_xxyyyzzz_0,  \
                             g_0_xxxxx_0_xxyyyzzz_1,  \
                             g_0_xxxxx_0_xxyyzzz_1,   \
                             g_0_xxxxx_0_xxyyzzzz_0,  \
                             g_0_xxxxx_0_xxyyzzzz_1,  \
                             g_0_xxxxx_0_xxyzzzz_1,   \
                             g_0_xxxxx_0_xxyzzzzz_0,  \
                             g_0_xxxxx_0_xxyzzzzz_1,  \
                             g_0_xxxxx_0_xxzzzzz_1,   \
                             g_0_xxxxx_0_xxzzzzzz_0,  \
                             g_0_xxxxx_0_xxzzzzzz_1,  \
                             g_0_xxxxx_0_xyyyyyy_1,   \
                             g_0_xxxxx_0_xyyyyyyy_0,  \
                             g_0_xxxxx_0_xyyyyyyy_1,  \
                             g_0_xxxxx_0_xyyyyyyz_0,  \
                             g_0_xxxxx_0_xyyyyyyz_1,  \
                             g_0_xxxxx_0_xyyyyyz_1,   \
                             g_0_xxxxx_0_xyyyyyzz_0,  \
                             g_0_xxxxx_0_xyyyyyzz_1,  \
                             g_0_xxxxx_0_xyyyyzz_1,   \
                             g_0_xxxxx_0_xyyyyzzz_0,  \
                             g_0_xxxxx_0_xyyyyzzz_1,  \
                             g_0_xxxxx_0_xyyyzzz_1,   \
                             g_0_xxxxx_0_xyyyzzzz_0,  \
                             g_0_xxxxx_0_xyyyzzzz_1,  \
                             g_0_xxxxx_0_xyyzzzz_1,   \
                             g_0_xxxxx_0_xyyzzzzz_0,  \
                             g_0_xxxxx_0_xyyzzzzz_1,  \
                             g_0_xxxxx_0_xyzzzzz_1,   \
                             g_0_xxxxx_0_xyzzzzzz_0,  \
                             g_0_xxxxx_0_xyzzzzzz_1,  \
                             g_0_xxxxx_0_xzzzzzz_1,   \
                             g_0_xxxxx_0_xzzzzzzz_0,  \
                             g_0_xxxxx_0_xzzzzzzz_1,  \
                             g_0_xxxxx_0_yyyyyyy_1,   \
                             g_0_xxxxx_0_yyyyyyyy_0,  \
                             g_0_xxxxx_0_yyyyyyyy_1,  \
                             g_0_xxxxx_0_yyyyyyyz_0,  \
                             g_0_xxxxx_0_yyyyyyyz_1,  \
                             g_0_xxxxx_0_yyyyyyz_1,   \
                             g_0_xxxxx_0_yyyyyyzz_0,  \
                             g_0_xxxxx_0_yyyyyyzz_1,  \
                             g_0_xxxxx_0_yyyyyzz_1,   \
                             g_0_xxxxx_0_yyyyyzzz_0,  \
                             g_0_xxxxx_0_yyyyyzzz_1,  \
                             g_0_xxxxx_0_yyyyzzz_1,   \
                             g_0_xxxxx_0_yyyyzzzz_0,  \
                             g_0_xxxxx_0_yyyyzzzz_1,  \
                             g_0_xxxxx_0_yyyzzzz_1,   \
                             g_0_xxxxx_0_yyyzzzzz_0,  \
                             g_0_xxxxx_0_yyyzzzzz_1,  \
                             g_0_xxxxx_0_yyzzzzz_1,   \
                             g_0_xxxxx_0_yyzzzzzz_0,  \
                             g_0_xxxxx_0_yyzzzzzz_1,  \
                             g_0_xxxxx_0_yzzzzzz_1,   \
                             g_0_xxxxx_0_yzzzzzzz_0,  \
                             g_0_xxxxx_0_yzzzzzzz_1,  \
                             g_0_xxxxx_0_zzzzzzz_1,   \
                             g_0_xxxxx_0_zzzzzzzz_0,  \
                             g_0_xxxxx_0_zzzzzzzz_1,  \
                             g_0_xxxxxy_0_xxxxxxxx_0, \
                             g_0_xxxxxy_0_xxxxxxxy_0, \
                             g_0_xxxxxy_0_xxxxxxxz_0, \
                             g_0_xxxxxy_0_xxxxxxyy_0, \
                             g_0_xxxxxy_0_xxxxxxyz_0, \
                             g_0_xxxxxy_0_xxxxxxzz_0, \
                             g_0_xxxxxy_0_xxxxxyyy_0, \
                             g_0_xxxxxy_0_xxxxxyyz_0, \
                             g_0_xxxxxy_0_xxxxxyzz_0, \
                             g_0_xxxxxy_0_xxxxxzzz_0, \
                             g_0_xxxxxy_0_xxxxyyyy_0, \
                             g_0_xxxxxy_0_xxxxyyyz_0, \
                             g_0_xxxxxy_0_xxxxyyzz_0, \
                             g_0_xxxxxy_0_xxxxyzzz_0, \
                             g_0_xxxxxy_0_xxxxzzzz_0, \
                             g_0_xxxxxy_0_xxxyyyyy_0, \
                             g_0_xxxxxy_0_xxxyyyyz_0, \
                             g_0_xxxxxy_0_xxxyyyzz_0, \
                             g_0_xxxxxy_0_xxxyyzzz_0, \
                             g_0_xxxxxy_0_xxxyzzzz_0, \
                             g_0_xxxxxy_0_xxxzzzzz_0, \
                             g_0_xxxxxy_0_xxyyyyyy_0, \
                             g_0_xxxxxy_0_xxyyyyyz_0, \
                             g_0_xxxxxy_0_xxyyyyzz_0, \
                             g_0_xxxxxy_0_xxyyyzzz_0, \
                             g_0_xxxxxy_0_xxyyzzzz_0, \
                             g_0_xxxxxy_0_xxyzzzzz_0, \
                             g_0_xxxxxy_0_xxzzzzzz_0, \
                             g_0_xxxxxy_0_xyyyyyyy_0, \
                             g_0_xxxxxy_0_xyyyyyyz_0, \
                             g_0_xxxxxy_0_xyyyyyzz_0, \
                             g_0_xxxxxy_0_xyyyyzzz_0, \
                             g_0_xxxxxy_0_xyyyzzzz_0, \
                             g_0_xxxxxy_0_xyyzzzzz_0, \
                             g_0_xxxxxy_0_xyzzzzzz_0, \
                             g_0_xxxxxy_0_xzzzzzzz_0, \
                             g_0_xxxxxy_0_yyyyyyyy_0, \
                             g_0_xxxxxy_0_yyyyyyyz_0, \
                             g_0_xxxxxy_0_yyyyyyzz_0, \
                             g_0_xxxxxy_0_yyyyyzzz_0, \
                             g_0_xxxxxy_0_yyyyzzzz_0, \
                             g_0_xxxxxy_0_yyyzzzzz_0, \
                             g_0_xxxxxy_0_yyzzzzzz_0, \
                             g_0_xxxxxy_0_yzzzzzzz_0, \
                             g_0_xxxxxy_0_zzzzzzzz_0, \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxy_0_xxxxxxxx_0[i] = g_0_xxxxx_0_xxxxxxxx_0[i] * pb_y + g_0_xxxxx_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxxxy_0[i] = g_0_xxxxx_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxxy_0[i] * pb_y + g_0_xxxxx_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxxxz_0[i] = g_0_xxxxx_0_xxxxxxxz_0[i] * pb_y + g_0_xxxxx_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxxyy_0[i] =
            2.0 * g_0_xxxxx_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxyy_0[i] * pb_y + g_0_xxxxx_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxxyz_0[i] = g_0_xxxxx_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxyz_0[i] * pb_y + g_0_xxxxx_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxxzz_0[i] = g_0_xxxxx_0_xxxxxxzz_0[i] * pb_y + g_0_xxxxx_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxyyy_0[i] =
            3.0 * g_0_xxxxx_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyyy_0[i] * pb_y + g_0_xxxxx_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxyyz_0[i] =
            2.0 * g_0_xxxxx_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyyz_0[i] * pb_y + g_0_xxxxx_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxyzz_0[i] = g_0_xxxxx_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyzz_0[i] * pb_y + g_0_xxxxx_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxzzz_0[i] = g_0_xxxxx_0_xxxxxzzz_0[i] * pb_y + g_0_xxxxx_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxyyyy_0[i] =
            4.0 * g_0_xxxxx_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyyy_0[i] * pb_y + g_0_xxxxx_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxyyyz_0[i] =
            3.0 * g_0_xxxxx_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyyz_0[i] * pb_y + g_0_xxxxx_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxyyzz_0[i] =
            2.0 * g_0_xxxxx_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyzz_0[i] * pb_y + g_0_xxxxx_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxyzzz_0[i] = g_0_xxxxx_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyzzz_0[i] * pb_y + g_0_xxxxx_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxzzzz_0[i] = g_0_xxxxx_0_xxxxzzzz_0[i] * pb_y + g_0_xxxxx_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyyyyy_0[i] =
            5.0 * g_0_xxxxx_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyyy_0[i] * pb_y + g_0_xxxxx_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyyyyz_0[i] =
            4.0 * g_0_xxxxx_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyyz_0[i] * pb_y + g_0_xxxxx_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyyyzz_0[i] =
            3.0 * g_0_xxxxx_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyzz_0[i] * pb_y + g_0_xxxxx_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyyzzz_0[i] =
            2.0 * g_0_xxxxx_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyzzz_0[i] * pb_y + g_0_xxxxx_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyzzzz_0[i] = g_0_xxxxx_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyzzzz_0[i] * pb_y + g_0_xxxxx_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxzzzzz_0[i] = g_0_xxxxx_0_xxxzzzzz_0[i] * pb_y + g_0_xxxxx_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyyyyy_0[i] =
            6.0 * g_0_xxxxx_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyyy_0[i] * pb_y + g_0_xxxxx_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyyyyz_0[i] =
            5.0 * g_0_xxxxx_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyyz_0[i] * pb_y + g_0_xxxxx_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyyyzz_0[i] =
            4.0 * g_0_xxxxx_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyzz_0[i] * pb_y + g_0_xxxxx_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyyzzz_0[i] =
            3.0 * g_0_xxxxx_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyzzz_0[i] * pb_y + g_0_xxxxx_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyzzzz_0[i] =
            2.0 * g_0_xxxxx_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyzzzz_0[i] * pb_y + g_0_xxxxx_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyzzzzz_0[i] = g_0_xxxxx_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyzzzzz_0[i] * pb_y + g_0_xxxxx_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxzzzzzz_0[i] = g_0_xxxxx_0_xxzzzzzz_0[i] * pb_y + g_0_xxxxx_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyyyyy_0[i] =
            7.0 * g_0_xxxxx_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyyy_0[i] * pb_y + g_0_xxxxx_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyyyyz_0[i] =
            6.0 * g_0_xxxxx_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyyz_0[i] * pb_y + g_0_xxxxx_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyyyzz_0[i] =
            5.0 * g_0_xxxxx_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyzz_0[i] * pb_y + g_0_xxxxx_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyyzzz_0[i] =
            4.0 * g_0_xxxxx_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyzzz_0[i] * pb_y + g_0_xxxxx_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyzzzz_0[i] =
            3.0 * g_0_xxxxx_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyzzzz_0[i] * pb_y + g_0_xxxxx_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyzzzzz_0[i] =
            2.0 * g_0_xxxxx_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyzzzzz_0[i] * pb_y + g_0_xxxxx_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyzzzzzz_0[i] = g_0_xxxxx_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzzzzzz_0[i] * pb_y + g_0_xxxxx_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xzzzzzzz_0[i] = g_0_xxxxx_0_xzzzzzzz_0[i] * pb_y + g_0_xxxxx_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyyyyy_0[i] =
            8.0 * g_0_xxxxx_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyyyy_0[i] * pb_y + g_0_xxxxx_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyyyyz_0[i] =
            7.0 * g_0_xxxxx_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyyyz_0[i] * pb_y + g_0_xxxxx_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyyyzz_0[i] =
            6.0 * g_0_xxxxx_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyyzz_0[i] * pb_y + g_0_xxxxx_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyyzzz_0[i] =
            5.0 * g_0_xxxxx_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyzzz_0[i] * pb_y + g_0_xxxxx_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyzzzz_0[i] =
            4.0 * g_0_xxxxx_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyzzzz_0[i] * pb_y + g_0_xxxxx_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyzzzzz_0[i] =
            3.0 * g_0_xxxxx_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyzzzzz_0[i] * pb_y + g_0_xxxxx_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyzzzzzz_0[i] =
            2.0 * g_0_xxxxx_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyzzzzzz_0[i] * pb_y + g_0_xxxxx_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yzzzzzzz_0[i] = g_0_xxxxx_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yzzzzzzz_0[i] * pb_y + g_0_xxxxx_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_zzzzzzzz_0[i] = g_0_xxxxx_0_zzzzzzzz_0[i] * pb_y + g_0_xxxxx_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 90-135 components of targeted buffer : SISL

    auto g_0_xxxxxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 90);

    auto g_0_xxxxxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 91);

    auto g_0_xxxxxz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 92);

    auto g_0_xxxxxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 93);

    auto g_0_xxxxxz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 94);

    auto g_0_xxxxxz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 95);

    auto g_0_xxxxxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 96);

    auto g_0_xxxxxz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 97);

    auto g_0_xxxxxz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 98);

    auto g_0_xxxxxz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 99);

    auto g_0_xxxxxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 100);

    auto g_0_xxxxxz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 101);

    auto g_0_xxxxxz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 102);

    auto g_0_xxxxxz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 103);

    auto g_0_xxxxxz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 104);

    auto g_0_xxxxxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 105);

    auto g_0_xxxxxz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 106);

    auto g_0_xxxxxz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 107);

    auto g_0_xxxxxz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 108);

    auto g_0_xxxxxz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 109);

    auto g_0_xxxxxz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 110);

    auto g_0_xxxxxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 111);

    auto g_0_xxxxxz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 112);

    auto g_0_xxxxxz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 113);

    auto g_0_xxxxxz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 114);

    auto g_0_xxxxxz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 115);

    auto g_0_xxxxxz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 116);

    auto g_0_xxxxxz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 117);

    auto g_0_xxxxxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 118);

    auto g_0_xxxxxz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 119);

    auto g_0_xxxxxz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 120);

    auto g_0_xxxxxz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 121);

    auto g_0_xxxxxz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 122);

    auto g_0_xxxxxz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 123);

    auto g_0_xxxxxz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 124);

    auto g_0_xxxxxz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 125);

    auto g_0_xxxxxz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 126);

    auto g_0_xxxxxz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 127);

    auto g_0_xxxxxz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 128);

    auto g_0_xxxxxz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 129);

    auto g_0_xxxxxz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 130);

    auto g_0_xxxxxz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 131);

    auto g_0_xxxxxz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 132);

    auto g_0_xxxxxz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 133);

    auto g_0_xxxxxz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 134);

#pragma omp simd aligned(g_0_xxxxx_0_xxxxxxx_1,       \
                             g_0_xxxxx_0_xxxxxxxx_0,  \
                             g_0_xxxxx_0_xxxxxxxx_1,  \
                             g_0_xxxxx_0_xxxxxxxy_0,  \
                             g_0_xxxxx_0_xxxxxxxy_1,  \
                             g_0_xxxxx_0_xxxxxxxz_0,  \
                             g_0_xxxxx_0_xxxxxxxz_1,  \
                             g_0_xxxxx_0_xxxxxxy_1,   \
                             g_0_xxxxx_0_xxxxxxyy_0,  \
                             g_0_xxxxx_0_xxxxxxyy_1,  \
                             g_0_xxxxx_0_xxxxxxyz_0,  \
                             g_0_xxxxx_0_xxxxxxyz_1,  \
                             g_0_xxxxx_0_xxxxxxz_1,   \
                             g_0_xxxxx_0_xxxxxxzz_0,  \
                             g_0_xxxxx_0_xxxxxxzz_1,  \
                             g_0_xxxxx_0_xxxxxyy_1,   \
                             g_0_xxxxx_0_xxxxxyyy_0,  \
                             g_0_xxxxx_0_xxxxxyyy_1,  \
                             g_0_xxxxx_0_xxxxxyyz_0,  \
                             g_0_xxxxx_0_xxxxxyyz_1,  \
                             g_0_xxxxx_0_xxxxxyz_1,   \
                             g_0_xxxxx_0_xxxxxyzz_0,  \
                             g_0_xxxxx_0_xxxxxyzz_1,  \
                             g_0_xxxxx_0_xxxxxzz_1,   \
                             g_0_xxxxx_0_xxxxxzzz_0,  \
                             g_0_xxxxx_0_xxxxxzzz_1,  \
                             g_0_xxxxx_0_xxxxyyy_1,   \
                             g_0_xxxxx_0_xxxxyyyy_0,  \
                             g_0_xxxxx_0_xxxxyyyy_1,  \
                             g_0_xxxxx_0_xxxxyyyz_0,  \
                             g_0_xxxxx_0_xxxxyyyz_1,  \
                             g_0_xxxxx_0_xxxxyyz_1,   \
                             g_0_xxxxx_0_xxxxyyzz_0,  \
                             g_0_xxxxx_0_xxxxyyzz_1,  \
                             g_0_xxxxx_0_xxxxyzz_1,   \
                             g_0_xxxxx_0_xxxxyzzz_0,  \
                             g_0_xxxxx_0_xxxxyzzz_1,  \
                             g_0_xxxxx_0_xxxxzzz_1,   \
                             g_0_xxxxx_0_xxxxzzzz_0,  \
                             g_0_xxxxx_0_xxxxzzzz_1,  \
                             g_0_xxxxx_0_xxxyyyy_1,   \
                             g_0_xxxxx_0_xxxyyyyy_0,  \
                             g_0_xxxxx_0_xxxyyyyy_1,  \
                             g_0_xxxxx_0_xxxyyyyz_0,  \
                             g_0_xxxxx_0_xxxyyyyz_1,  \
                             g_0_xxxxx_0_xxxyyyz_1,   \
                             g_0_xxxxx_0_xxxyyyzz_0,  \
                             g_0_xxxxx_0_xxxyyyzz_1,  \
                             g_0_xxxxx_0_xxxyyzz_1,   \
                             g_0_xxxxx_0_xxxyyzzz_0,  \
                             g_0_xxxxx_0_xxxyyzzz_1,  \
                             g_0_xxxxx_0_xxxyzzz_1,   \
                             g_0_xxxxx_0_xxxyzzzz_0,  \
                             g_0_xxxxx_0_xxxyzzzz_1,  \
                             g_0_xxxxx_0_xxxzzzz_1,   \
                             g_0_xxxxx_0_xxxzzzzz_0,  \
                             g_0_xxxxx_0_xxxzzzzz_1,  \
                             g_0_xxxxx_0_xxyyyyy_1,   \
                             g_0_xxxxx_0_xxyyyyyy_0,  \
                             g_0_xxxxx_0_xxyyyyyy_1,  \
                             g_0_xxxxx_0_xxyyyyyz_0,  \
                             g_0_xxxxx_0_xxyyyyyz_1,  \
                             g_0_xxxxx_0_xxyyyyz_1,   \
                             g_0_xxxxx_0_xxyyyyzz_0,  \
                             g_0_xxxxx_0_xxyyyyzz_1,  \
                             g_0_xxxxx_0_xxyyyzz_1,   \
                             g_0_xxxxx_0_xxyyyzzz_0,  \
                             g_0_xxxxx_0_xxyyyzzz_1,  \
                             g_0_xxxxx_0_xxyyzzz_1,   \
                             g_0_xxxxx_0_xxyyzzzz_0,  \
                             g_0_xxxxx_0_xxyyzzzz_1,  \
                             g_0_xxxxx_0_xxyzzzz_1,   \
                             g_0_xxxxx_0_xxyzzzzz_0,  \
                             g_0_xxxxx_0_xxyzzzzz_1,  \
                             g_0_xxxxx_0_xxzzzzz_1,   \
                             g_0_xxxxx_0_xxzzzzzz_0,  \
                             g_0_xxxxx_0_xxzzzzzz_1,  \
                             g_0_xxxxx_0_xyyyyyy_1,   \
                             g_0_xxxxx_0_xyyyyyyy_0,  \
                             g_0_xxxxx_0_xyyyyyyy_1,  \
                             g_0_xxxxx_0_xyyyyyyz_0,  \
                             g_0_xxxxx_0_xyyyyyyz_1,  \
                             g_0_xxxxx_0_xyyyyyz_1,   \
                             g_0_xxxxx_0_xyyyyyzz_0,  \
                             g_0_xxxxx_0_xyyyyyzz_1,  \
                             g_0_xxxxx_0_xyyyyzz_1,   \
                             g_0_xxxxx_0_xyyyyzzz_0,  \
                             g_0_xxxxx_0_xyyyyzzz_1,  \
                             g_0_xxxxx_0_xyyyzzz_1,   \
                             g_0_xxxxx_0_xyyyzzzz_0,  \
                             g_0_xxxxx_0_xyyyzzzz_1,  \
                             g_0_xxxxx_0_xyyzzzz_1,   \
                             g_0_xxxxx_0_xyyzzzzz_0,  \
                             g_0_xxxxx_0_xyyzzzzz_1,  \
                             g_0_xxxxx_0_xyzzzzz_1,   \
                             g_0_xxxxx_0_xyzzzzzz_0,  \
                             g_0_xxxxx_0_xyzzzzzz_1,  \
                             g_0_xxxxx_0_xzzzzzz_1,   \
                             g_0_xxxxx_0_xzzzzzzz_0,  \
                             g_0_xxxxx_0_xzzzzzzz_1,  \
                             g_0_xxxxx_0_yyyyyyy_1,   \
                             g_0_xxxxx_0_yyyyyyyy_0,  \
                             g_0_xxxxx_0_yyyyyyyy_1,  \
                             g_0_xxxxx_0_yyyyyyyz_0,  \
                             g_0_xxxxx_0_yyyyyyyz_1,  \
                             g_0_xxxxx_0_yyyyyyz_1,   \
                             g_0_xxxxx_0_yyyyyyzz_0,  \
                             g_0_xxxxx_0_yyyyyyzz_1,  \
                             g_0_xxxxx_0_yyyyyzz_1,   \
                             g_0_xxxxx_0_yyyyyzzz_0,  \
                             g_0_xxxxx_0_yyyyyzzz_1,  \
                             g_0_xxxxx_0_yyyyzzz_1,   \
                             g_0_xxxxx_0_yyyyzzzz_0,  \
                             g_0_xxxxx_0_yyyyzzzz_1,  \
                             g_0_xxxxx_0_yyyzzzz_1,   \
                             g_0_xxxxx_0_yyyzzzzz_0,  \
                             g_0_xxxxx_0_yyyzzzzz_1,  \
                             g_0_xxxxx_0_yyzzzzz_1,   \
                             g_0_xxxxx_0_yyzzzzzz_0,  \
                             g_0_xxxxx_0_yyzzzzzz_1,  \
                             g_0_xxxxx_0_yzzzzzz_1,   \
                             g_0_xxxxx_0_yzzzzzzz_0,  \
                             g_0_xxxxx_0_yzzzzzzz_1,  \
                             g_0_xxxxx_0_zzzzzzz_1,   \
                             g_0_xxxxx_0_zzzzzzzz_0,  \
                             g_0_xxxxx_0_zzzzzzzz_1,  \
                             g_0_xxxxxz_0_xxxxxxxx_0, \
                             g_0_xxxxxz_0_xxxxxxxy_0, \
                             g_0_xxxxxz_0_xxxxxxxz_0, \
                             g_0_xxxxxz_0_xxxxxxyy_0, \
                             g_0_xxxxxz_0_xxxxxxyz_0, \
                             g_0_xxxxxz_0_xxxxxxzz_0, \
                             g_0_xxxxxz_0_xxxxxyyy_0, \
                             g_0_xxxxxz_0_xxxxxyyz_0, \
                             g_0_xxxxxz_0_xxxxxyzz_0, \
                             g_0_xxxxxz_0_xxxxxzzz_0, \
                             g_0_xxxxxz_0_xxxxyyyy_0, \
                             g_0_xxxxxz_0_xxxxyyyz_0, \
                             g_0_xxxxxz_0_xxxxyyzz_0, \
                             g_0_xxxxxz_0_xxxxyzzz_0, \
                             g_0_xxxxxz_0_xxxxzzzz_0, \
                             g_0_xxxxxz_0_xxxyyyyy_0, \
                             g_0_xxxxxz_0_xxxyyyyz_0, \
                             g_0_xxxxxz_0_xxxyyyzz_0, \
                             g_0_xxxxxz_0_xxxyyzzz_0, \
                             g_0_xxxxxz_0_xxxyzzzz_0, \
                             g_0_xxxxxz_0_xxxzzzzz_0, \
                             g_0_xxxxxz_0_xxyyyyyy_0, \
                             g_0_xxxxxz_0_xxyyyyyz_0, \
                             g_0_xxxxxz_0_xxyyyyzz_0, \
                             g_0_xxxxxz_0_xxyyyzzz_0, \
                             g_0_xxxxxz_0_xxyyzzzz_0, \
                             g_0_xxxxxz_0_xxyzzzzz_0, \
                             g_0_xxxxxz_0_xxzzzzzz_0, \
                             g_0_xxxxxz_0_xyyyyyyy_0, \
                             g_0_xxxxxz_0_xyyyyyyz_0, \
                             g_0_xxxxxz_0_xyyyyyzz_0, \
                             g_0_xxxxxz_0_xyyyyzzz_0, \
                             g_0_xxxxxz_0_xyyyzzzz_0, \
                             g_0_xxxxxz_0_xyyzzzzz_0, \
                             g_0_xxxxxz_0_xyzzzzzz_0, \
                             g_0_xxxxxz_0_xzzzzzzz_0, \
                             g_0_xxxxxz_0_yyyyyyyy_0, \
                             g_0_xxxxxz_0_yyyyyyyz_0, \
                             g_0_xxxxxz_0_yyyyyyzz_0, \
                             g_0_xxxxxz_0_yyyyyzzz_0, \
                             g_0_xxxxxz_0_yyyyzzzz_0, \
                             g_0_xxxxxz_0_yyyzzzzz_0, \
                             g_0_xxxxxz_0_yyzzzzzz_0, \
                             g_0_xxxxxz_0_yzzzzzzz_0, \
                             g_0_xxxxxz_0_zzzzzzzz_0, \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxz_0_xxxxxxxx_0[i] = g_0_xxxxx_0_xxxxxxxx_0[i] * pb_z + g_0_xxxxx_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxxxy_0[i] = g_0_xxxxx_0_xxxxxxxy_0[i] * pb_z + g_0_xxxxx_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxxxz_0[i] = g_0_xxxxx_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxxz_0[i] * pb_z + g_0_xxxxx_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxxyy_0[i] = g_0_xxxxx_0_xxxxxxyy_0[i] * pb_z + g_0_xxxxx_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxxyz_0[i] = g_0_xxxxx_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxyz_0[i] * pb_z + g_0_xxxxx_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxxzz_0[i] =
            2.0 * g_0_xxxxx_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxzz_0[i] * pb_z + g_0_xxxxx_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxyyy_0[i] = g_0_xxxxx_0_xxxxxyyy_0[i] * pb_z + g_0_xxxxx_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxyyz_0[i] = g_0_xxxxx_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyyz_0[i] * pb_z + g_0_xxxxx_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxyzz_0[i] =
            2.0 * g_0_xxxxx_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyzz_0[i] * pb_z + g_0_xxxxx_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxzzz_0[i] =
            3.0 * g_0_xxxxx_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxzzz_0[i] * pb_z + g_0_xxxxx_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxyyyy_0[i] = g_0_xxxxx_0_xxxxyyyy_0[i] * pb_z + g_0_xxxxx_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxyyyz_0[i] = g_0_xxxxx_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyyz_0[i] * pb_z + g_0_xxxxx_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxyyzz_0[i] =
            2.0 * g_0_xxxxx_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyzz_0[i] * pb_z + g_0_xxxxx_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxyzzz_0[i] =
            3.0 * g_0_xxxxx_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyzzz_0[i] * pb_z + g_0_xxxxx_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxzzzz_0[i] =
            4.0 * g_0_xxxxx_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxzzzz_0[i] * pb_z + g_0_xxxxx_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyyyyy_0[i] = g_0_xxxxx_0_xxxyyyyy_0[i] * pb_z + g_0_xxxxx_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyyyyz_0[i] = g_0_xxxxx_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyyz_0[i] * pb_z + g_0_xxxxx_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyyyzz_0[i] =
            2.0 * g_0_xxxxx_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyzz_0[i] * pb_z + g_0_xxxxx_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyyzzz_0[i] =
            3.0 * g_0_xxxxx_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyzzz_0[i] * pb_z + g_0_xxxxx_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyzzzz_0[i] =
            4.0 * g_0_xxxxx_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyzzzz_0[i] * pb_z + g_0_xxxxx_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxzzzzz_0[i] =
            5.0 * g_0_xxxxx_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxzzzzz_0[i] * pb_z + g_0_xxxxx_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyyyyy_0[i] = g_0_xxxxx_0_xxyyyyyy_0[i] * pb_z + g_0_xxxxx_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyyyyz_0[i] = g_0_xxxxx_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyyz_0[i] * pb_z + g_0_xxxxx_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyyyzz_0[i] =
            2.0 * g_0_xxxxx_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyzz_0[i] * pb_z + g_0_xxxxx_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyyzzz_0[i] =
            3.0 * g_0_xxxxx_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyzzz_0[i] * pb_z + g_0_xxxxx_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyzzzz_0[i] =
            4.0 * g_0_xxxxx_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyzzzz_0[i] * pb_z + g_0_xxxxx_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyzzzzz_0[i] =
            5.0 * g_0_xxxxx_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyzzzzz_0[i] * pb_z + g_0_xxxxx_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxzzzzzz_0[i] =
            6.0 * g_0_xxxxx_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxzzzzzz_0[i] * pb_z + g_0_xxxxx_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyyyyy_0[i] = g_0_xxxxx_0_xyyyyyyy_0[i] * pb_z + g_0_xxxxx_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyyyyz_0[i] = g_0_xxxxx_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyyz_0[i] * pb_z + g_0_xxxxx_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyyyzz_0[i] =
            2.0 * g_0_xxxxx_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyzz_0[i] * pb_z + g_0_xxxxx_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyyzzz_0[i] =
            3.0 * g_0_xxxxx_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyzzz_0[i] * pb_z + g_0_xxxxx_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyzzzz_0[i] =
            4.0 * g_0_xxxxx_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyzzzz_0[i] * pb_z + g_0_xxxxx_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyzzzzz_0[i] =
            5.0 * g_0_xxxxx_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyzzzzz_0[i] * pb_z + g_0_xxxxx_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyzzzzzz_0[i] =
            6.0 * g_0_xxxxx_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzzzzzz_0[i] * pb_z + g_0_xxxxx_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xzzzzzzz_0[i] =
            7.0 * g_0_xxxxx_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xzzzzzzz_0[i] * pb_z + g_0_xxxxx_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyyyyy_0[i] = g_0_xxxxx_0_yyyyyyyy_0[i] * pb_z + g_0_xxxxx_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyyyyz_0[i] = g_0_xxxxx_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyyyz_0[i] * pb_z + g_0_xxxxx_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyyyzz_0[i] =
            2.0 * g_0_xxxxx_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyyzz_0[i] * pb_z + g_0_xxxxx_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyyzzz_0[i] =
            3.0 * g_0_xxxxx_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyzzz_0[i] * pb_z + g_0_xxxxx_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyzzzz_0[i] =
            4.0 * g_0_xxxxx_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyzzzz_0[i] * pb_z + g_0_xxxxx_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyzzzzz_0[i] =
            5.0 * g_0_xxxxx_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyzzzzz_0[i] * pb_z + g_0_xxxxx_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyzzzzzz_0[i] =
            6.0 * g_0_xxxxx_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyzzzzzz_0[i] * pb_z + g_0_xxxxx_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yzzzzzzz_0[i] =
            7.0 * g_0_xxxxx_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yzzzzzzz_0[i] * pb_z + g_0_xxxxx_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_zzzzzzzz_0[i] =
            8.0 * g_0_xxxxx_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_zzzzzzzz_0[i] * pb_z + g_0_xxxxx_0_zzzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 135-180 components of targeted buffer : SISL

    auto g_0_xxxxyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 135);

    auto g_0_xxxxyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 136);

    auto g_0_xxxxyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 137);

    auto g_0_xxxxyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 138);

    auto g_0_xxxxyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 139);

    auto g_0_xxxxyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 140);

    auto g_0_xxxxyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 141);

    auto g_0_xxxxyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 142);

    auto g_0_xxxxyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 143);

    auto g_0_xxxxyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 144);

    auto g_0_xxxxyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 145);

    auto g_0_xxxxyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 146);

    auto g_0_xxxxyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 147);

    auto g_0_xxxxyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 148);

    auto g_0_xxxxyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 149);

    auto g_0_xxxxyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 150);

    auto g_0_xxxxyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 151);

    auto g_0_xxxxyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 152);

    auto g_0_xxxxyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 153);

    auto g_0_xxxxyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 154);

    auto g_0_xxxxyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 155);

    auto g_0_xxxxyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 156);

    auto g_0_xxxxyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 157);

    auto g_0_xxxxyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 158);

    auto g_0_xxxxyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 159);

    auto g_0_xxxxyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 160);

    auto g_0_xxxxyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 161);

    auto g_0_xxxxyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 162);

    auto g_0_xxxxyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 163);

    auto g_0_xxxxyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 164);

    auto g_0_xxxxyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 165);

    auto g_0_xxxxyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 166);

    auto g_0_xxxxyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 167);

    auto g_0_xxxxyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 168);

    auto g_0_xxxxyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 169);

    auto g_0_xxxxyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 170);

    auto g_0_xxxxyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 171);

    auto g_0_xxxxyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 172);

    auto g_0_xxxxyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 173);

    auto g_0_xxxxyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 174);

    auto g_0_xxxxyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 175);

    auto g_0_xxxxyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 176);

    auto g_0_xxxxyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 177);

    auto g_0_xxxxyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 178);

    auto g_0_xxxxyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 179);

#pragma omp simd aligned(g_0_xxxx_0_xxxxxxxx_0,       \
                             g_0_xxxx_0_xxxxxxxx_1,   \
                             g_0_xxxx_0_xxxxxxxz_0,   \
                             g_0_xxxx_0_xxxxxxxz_1,   \
                             g_0_xxxx_0_xxxxxxzz_0,   \
                             g_0_xxxx_0_xxxxxxzz_1,   \
                             g_0_xxxx_0_xxxxxzzz_0,   \
                             g_0_xxxx_0_xxxxxzzz_1,   \
                             g_0_xxxx_0_xxxxzzzz_0,   \
                             g_0_xxxx_0_xxxxzzzz_1,   \
                             g_0_xxxx_0_xxxzzzzz_0,   \
                             g_0_xxxx_0_xxxzzzzz_1,   \
                             g_0_xxxx_0_xxzzzzzz_0,   \
                             g_0_xxxx_0_xxzzzzzz_1,   \
                             g_0_xxxx_0_xzzzzzzz_0,   \
                             g_0_xxxx_0_xzzzzzzz_1,   \
                             g_0_xxxxy_0_xxxxxxxx_0,  \
                             g_0_xxxxy_0_xxxxxxxx_1,  \
                             g_0_xxxxy_0_xxxxxxxz_0,  \
                             g_0_xxxxy_0_xxxxxxxz_1,  \
                             g_0_xxxxy_0_xxxxxxzz_0,  \
                             g_0_xxxxy_0_xxxxxxzz_1,  \
                             g_0_xxxxy_0_xxxxxzzz_0,  \
                             g_0_xxxxy_0_xxxxxzzz_1,  \
                             g_0_xxxxy_0_xxxxzzzz_0,  \
                             g_0_xxxxy_0_xxxxzzzz_1,  \
                             g_0_xxxxy_0_xxxzzzzz_0,  \
                             g_0_xxxxy_0_xxxzzzzz_1,  \
                             g_0_xxxxy_0_xxzzzzzz_0,  \
                             g_0_xxxxy_0_xxzzzzzz_1,  \
                             g_0_xxxxy_0_xzzzzzzz_0,  \
                             g_0_xxxxy_0_xzzzzzzz_1,  \
                             g_0_xxxxyy_0_xxxxxxxx_0, \
                             g_0_xxxxyy_0_xxxxxxxy_0, \
                             g_0_xxxxyy_0_xxxxxxxz_0, \
                             g_0_xxxxyy_0_xxxxxxyy_0, \
                             g_0_xxxxyy_0_xxxxxxyz_0, \
                             g_0_xxxxyy_0_xxxxxxzz_0, \
                             g_0_xxxxyy_0_xxxxxyyy_0, \
                             g_0_xxxxyy_0_xxxxxyyz_0, \
                             g_0_xxxxyy_0_xxxxxyzz_0, \
                             g_0_xxxxyy_0_xxxxxzzz_0, \
                             g_0_xxxxyy_0_xxxxyyyy_0, \
                             g_0_xxxxyy_0_xxxxyyyz_0, \
                             g_0_xxxxyy_0_xxxxyyzz_0, \
                             g_0_xxxxyy_0_xxxxyzzz_0, \
                             g_0_xxxxyy_0_xxxxzzzz_0, \
                             g_0_xxxxyy_0_xxxyyyyy_0, \
                             g_0_xxxxyy_0_xxxyyyyz_0, \
                             g_0_xxxxyy_0_xxxyyyzz_0, \
                             g_0_xxxxyy_0_xxxyyzzz_0, \
                             g_0_xxxxyy_0_xxxyzzzz_0, \
                             g_0_xxxxyy_0_xxxzzzzz_0, \
                             g_0_xxxxyy_0_xxyyyyyy_0, \
                             g_0_xxxxyy_0_xxyyyyyz_0, \
                             g_0_xxxxyy_0_xxyyyyzz_0, \
                             g_0_xxxxyy_0_xxyyyzzz_0, \
                             g_0_xxxxyy_0_xxyyzzzz_0, \
                             g_0_xxxxyy_0_xxyzzzzz_0, \
                             g_0_xxxxyy_0_xxzzzzzz_0, \
                             g_0_xxxxyy_0_xyyyyyyy_0, \
                             g_0_xxxxyy_0_xyyyyyyz_0, \
                             g_0_xxxxyy_0_xyyyyyzz_0, \
                             g_0_xxxxyy_0_xyyyyzzz_0, \
                             g_0_xxxxyy_0_xyyyzzzz_0, \
                             g_0_xxxxyy_0_xyyzzzzz_0, \
                             g_0_xxxxyy_0_xyzzzzzz_0, \
                             g_0_xxxxyy_0_xzzzzzzz_0, \
                             g_0_xxxxyy_0_yyyyyyyy_0, \
                             g_0_xxxxyy_0_yyyyyyyz_0, \
                             g_0_xxxxyy_0_yyyyyyzz_0, \
                             g_0_xxxxyy_0_yyyyyzzz_0, \
                             g_0_xxxxyy_0_yyyyzzzz_0, \
                             g_0_xxxxyy_0_yyyzzzzz_0, \
                             g_0_xxxxyy_0_yyzzzzzz_0, \
                             g_0_xxxxyy_0_yzzzzzzz_0, \
                             g_0_xxxxyy_0_zzzzzzzz_0, \
                             g_0_xxxyy_0_xxxxxxxy_0,  \
                             g_0_xxxyy_0_xxxxxxxy_1,  \
                             g_0_xxxyy_0_xxxxxxy_1,   \
                             g_0_xxxyy_0_xxxxxxyy_0,  \
                             g_0_xxxyy_0_xxxxxxyy_1,  \
                             g_0_xxxyy_0_xxxxxxyz_0,  \
                             g_0_xxxyy_0_xxxxxxyz_1,  \
                             g_0_xxxyy_0_xxxxxyy_1,   \
                             g_0_xxxyy_0_xxxxxyyy_0,  \
                             g_0_xxxyy_0_xxxxxyyy_1,  \
                             g_0_xxxyy_0_xxxxxyyz_0,  \
                             g_0_xxxyy_0_xxxxxyyz_1,  \
                             g_0_xxxyy_0_xxxxxyz_1,   \
                             g_0_xxxyy_0_xxxxxyzz_0,  \
                             g_0_xxxyy_0_xxxxxyzz_1,  \
                             g_0_xxxyy_0_xxxxyyy_1,   \
                             g_0_xxxyy_0_xxxxyyyy_0,  \
                             g_0_xxxyy_0_xxxxyyyy_1,  \
                             g_0_xxxyy_0_xxxxyyyz_0,  \
                             g_0_xxxyy_0_xxxxyyyz_1,  \
                             g_0_xxxyy_0_xxxxyyz_1,   \
                             g_0_xxxyy_0_xxxxyyzz_0,  \
                             g_0_xxxyy_0_xxxxyyzz_1,  \
                             g_0_xxxyy_0_xxxxyzz_1,   \
                             g_0_xxxyy_0_xxxxyzzz_0,  \
                             g_0_xxxyy_0_xxxxyzzz_1,  \
                             g_0_xxxyy_0_xxxyyyy_1,   \
                             g_0_xxxyy_0_xxxyyyyy_0,  \
                             g_0_xxxyy_0_xxxyyyyy_1,  \
                             g_0_xxxyy_0_xxxyyyyz_0,  \
                             g_0_xxxyy_0_xxxyyyyz_1,  \
                             g_0_xxxyy_0_xxxyyyz_1,   \
                             g_0_xxxyy_0_xxxyyyzz_0,  \
                             g_0_xxxyy_0_xxxyyyzz_1,  \
                             g_0_xxxyy_0_xxxyyzz_1,   \
                             g_0_xxxyy_0_xxxyyzzz_0,  \
                             g_0_xxxyy_0_xxxyyzzz_1,  \
                             g_0_xxxyy_0_xxxyzzz_1,   \
                             g_0_xxxyy_0_xxxyzzzz_0,  \
                             g_0_xxxyy_0_xxxyzzzz_1,  \
                             g_0_xxxyy_0_xxyyyyy_1,   \
                             g_0_xxxyy_0_xxyyyyyy_0,  \
                             g_0_xxxyy_0_xxyyyyyy_1,  \
                             g_0_xxxyy_0_xxyyyyyz_0,  \
                             g_0_xxxyy_0_xxyyyyyz_1,  \
                             g_0_xxxyy_0_xxyyyyz_1,   \
                             g_0_xxxyy_0_xxyyyyzz_0,  \
                             g_0_xxxyy_0_xxyyyyzz_1,  \
                             g_0_xxxyy_0_xxyyyzz_1,   \
                             g_0_xxxyy_0_xxyyyzzz_0,  \
                             g_0_xxxyy_0_xxyyyzzz_1,  \
                             g_0_xxxyy_0_xxyyzzz_1,   \
                             g_0_xxxyy_0_xxyyzzzz_0,  \
                             g_0_xxxyy_0_xxyyzzzz_1,  \
                             g_0_xxxyy_0_xxyzzzz_1,   \
                             g_0_xxxyy_0_xxyzzzzz_0,  \
                             g_0_xxxyy_0_xxyzzzzz_1,  \
                             g_0_xxxyy_0_xyyyyyy_1,   \
                             g_0_xxxyy_0_xyyyyyyy_0,  \
                             g_0_xxxyy_0_xyyyyyyy_1,  \
                             g_0_xxxyy_0_xyyyyyyz_0,  \
                             g_0_xxxyy_0_xyyyyyyz_1,  \
                             g_0_xxxyy_0_xyyyyyz_1,   \
                             g_0_xxxyy_0_xyyyyyzz_0,  \
                             g_0_xxxyy_0_xyyyyyzz_1,  \
                             g_0_xxxyy_0_xyyyyzz_1,   \
                             g_0_xxxyy_0_xyyyyzzz_0,  \
                             g_0_xxxyy_0_xyyyyzzz_1,  \
                             g_0_xxxyy_0_xyyyzzz_1,   \
                             g_0_xxxyy_0_xyyyzzzz_0,  \
                             g_0_xxxyy_0_xyyyzzzz_1,  \
                             g_0_xxxyy_0_xyyzzzz_1,   \
                             g_0_xxxyy_0_xyyzzzzz_0,  \
                             g_0_xxxyy_0_xyyzzzzz_1,  \
                             g_0_xxxyy_0_xyzzzzz_1,   \
                             g_0_xxxyy_0_xyzzzzzz_0,  \
                             g_0_xxxyy_0_xyzzzzzz_1,  \
                             g_0_xxxyy_0_yyyyyyy_1,   \
                             g_0_xxxyy_0_yyyyyyyy_0,  \
                             g_0_xxxyy_0_yyyyyyyy_1,  \
                             g_0_xxxyy_0_yyyyyyyz_0,  \
                             g_0_xxxyy_0_yyyyyyyz_1,  \
                             g_0_xxxyy_0_yyyyyyz_1,   \
                             g_0_xxxyy_0_yyyyyyzz_0,  \
                             g_0_xxxyy_0_yyyyyyzz_1,  \
                             g_0_xxxyy_0_yyyyyzz_1,   \
                             g_0_xxxyy_0_yyyyyzzz_0,  \
                             g_0_xxxyy_0_yyyyyzzz_1,  \
                             g_0_xxxyy_0_yyyyzzz_1,   \
                             g_0_xxxyy_0_yyyyzzzz_0,  \
                             g_0_xxxyy_0_yyyyzzzz_1,  \
                             g_0_xxxyy_0_yyyzzzz_1,   \
                             g_0_xxxyy_0_yyyzzzzz_0,  \
                             g_0_xxxyy_0_yyyzzzzz_1,  \
                             g_0_xxxyy_0_yyzzzzz_1,   \
                             g_0_xxxyy_0_yyzzzzzz_0,  \
                             g_0_xxxyy_0_yyzzzzzz_1,  \
                             g_0_xxxyy_0_yzzzzzz_1,   \
                             g_0_xxxyy_0_yzzzzzzz_0,  \
                             g_0_xxxyy_0_yzzzzzzz_1,  \
                             g_0_xxxyy_0_zzzzzzzz_0,  \
                             g_0_xxxyy_0_zzzzzzzz_1,  \
                             g_0_xxyy_0_xxxxxxxy_0,   \
                             g_0_xxyy_0_xxxxxxxy_1,   \
                             g_0_xxyy_0_xxxxxxyy_0,   \
                             g_0_xxyy_0_xxxxxxyy_1,   \
                             g_0_xxyy_0_xxxxxxyz_0,   \
                             g_0_xxyy_0_xxxxxxyz_1,   \
                             g_0_xxyy_0_xxxxxyyy_0,   \
                             g_0_xxyy_0_xxxxxyyy_1,   \
                             g_0_xxyy_0_xxxxxyyz_0,   \
                             g_0_xxyy_0_xxxxxyyz_1,   \
                             g_0_xxyy_0_xxxxxyzz_0,   \
                             g_0_xxyy_0_xxxxxyzz_1,   \
                             g_0_xxyy_0_xxxxyyyy_0,   \
                             g_0_xxyy_0_xxxxyyyy_1,   \
                             g_0_xxyy_0_xxxxyyyz_0,   \
                             g_0_xxyy_0_xxxxyyyz_1,   \
                             g_0_xxyy_0_xxxxyyzz_0,   \
                             g_0_xxyy_0_xxxxyyzz_1,   \
                             g_0_xxyy_0_xxxxyzzz_0,   \
                             g_0_xxyy_0_xxxxyzzz_1,   \
                             g_0_xxyy_0_xxxyyyyy_0,   \
                             g_0_xxyy_0_xxxyyyyy_1,   \
                             g_0_xxyy_0_xxxyyyyz_0,   \
                             g_0_xxyy_0_xxxyyyyz_1,   \
                             g_0_xxyy_0_xxxyyyzz_0,   \
                             g_0_xxyy_0_xxxyyyzz_1,   \
                             g_0_xxyy_0_xxxyyzzz_0,   \
                             g_0_xxyy_0_xxxyyzzz_1,   \
                             g_0_xxyy_0_xxxyzzzz_0,   \
                             g_0_xxyy_0_xxxyzzzz_1,   \
                             g_0_xxyy_0_xxyyyyyy_0,   \
                             g_0_xxyy_0_xxyyyyyy_1,   \
                             g_0_xxyy_0_xxyyyyyz_0,   \
                             g_0_xxyy_0_xxyyyyyz_1,   \
                             g_0_xxyy_0_xxyyyyzz_0,   \
                             g_0_xxyy_0_xxyyyyzz_1,   \
                             g_0_xxyy_0_xxyyyzzz_0,   \
                             g_0_xxyy_0_xxyyyzzz_1,   \
                             g_0_xxyy_0_xxyyzzzz_0,   \
                             g_0_xxyy_0_xxyyzzzz_1,   \
                             g_0_xxyy_0_xxyzzzzz_0,   \
                             g_0_xxyy_0_xxyzzzzz_1,   \
                             g_0_xxyy_0_xyyyyyyy_0,   \
                             g_0_xxyy_0_xyyyyyyy_1,   \
                             g_0_xxyy_0_xyyyyyyz_0,   \
                             g_0_xxyy_0_xyyyyyyz_1,   \
                             g_0_xxyy_0_xyyyyyzz_0,   \
                             g_0_xxyy_0_xyyyyyzz_1,   \
                             g_0_xxyy_0_xyyyyzzz_0,   \
                             g_0_xxyy_0_xyyyyzzz_1,   \
                             g_0_xxyy_0_xyyyzzzz_0,   \
                             g_0_xxyy_0_xyyyzzzz_1,   \
                             g_0_xxyy_0_xyyzzzzz_0,   \
                             g_0_xxyy_0_xyyzzzzz_1,   \
                             g_0_xxyy_0_xyzzzzzz_0,   \
                             g_0_xxyy_0_xyzzzzzz_1,   \
                             g_0_xxyy_0_yyyyyyyy_0,   \
                             g_0_xxyy_0_yyyyyyyy_1,   \
                             g_0_xxyy_0_yyyyyyyz_0,   \
                             g_0_xxyy_0_yyyyyyyz_1,   \
                             g_0_xxyy_0_yyyyyyzz_0,   \
                             g_0_xxyy_0_yyyyyyzz_1,   \
                             g_0_xxyy_0_yyyyyzzz_0,   \
                             g_0_xxyy_0_yyyyyzzz_1,   \
                             g_0_xxyy_0_yyyyzzzz_0,   \
                             g_0_xxyy_0_yyyyzzzz_1,   \
                             g_0_xxyy_0_yyyzzzzz_0,   \
                             g_0_xxyy_0_yyyzzzzz_1,   \
                             g_0_xxyy_0_yyzzzzzz_0,   \
                             g_0_xxyy_0_yyzzzzzz_1,   \
                             g_0_xxyy_0_yzzzzzzz_0,   \
                             g_0_xxyy_0_yzzzzzzz_1,   \
                             g_0_xxyy_0_zzzzzzzz_0,   \
                             g_0_xxyy_0_zzzzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyy_0_xxxxxxxx_0[i] = g_0_xxxx_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxxxxx_0[i] * pb_y +
                                     g_0_xxxxy_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxxxxxy_0[i] = 3.0 * g_0_xxyy_0_xxxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxxxy_1[i] * fti_ab_0 +
                                     7.0 * g_0_xxxyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxxxy_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxxxxz_0[i] = g_0_xxxx_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxxxxz_0[i] * pb_y +
                                     g_0_xxxxy_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxxxxyy_0[i] = 3.0 * g_0_xxyy_0_xxxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxxyy_1[i] * fti_ab_0 +
                                     6.0 * g_0_xxxyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxxyy_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxxxyz_0[i] = 3.0 * g_0_xxyy_0_xxxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxxyz_1[i] * fti_ab_0 +
                                     6.0 * g_0_xxxyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxxyz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxxxzz_0[i] = g_0_xxxx_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxxxzz_0[i] * pb_y +
                                     g_0_xxxxy_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxxxyyy_0[i] = 3.0 * g_0_xxyy_0_xxxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxyyy_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxyyy_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxxyyz_0[i] = 3.0 * g_0_xxyy_0_xxxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxyyz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxxyzz_0[i] = 3.0 * g_0_xxyy_0_xxxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxyzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxyzz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxxzzz_0[i] = g_0_xxxx_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxxzzz_0[i] * pb_y +
                                     g_0_xxxxy_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxxyyyy_0[i] = 3.0 * g_0_xxyy_0_xxxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxyyyy_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyyyy_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxyyyz_0[i] = 3.0 * g_0_xxyy_0_xxxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyyyz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxyyzz_0[i] = 3.0 * g_0_xxyy_0_xxxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyyzz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxyzzz_0[i] = 3.0 * g_0_xxyy_0_xxxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxyzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyzzz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxzzzz_0[i] = g_0_xxxx_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxzzzz_0[i] * pb_y +
                                     g_0_xxxxy_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxyyyyy_0[i] = 3.0 * g_0_xxyy_0_xxxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyyyyy_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyyyy_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxyyyyz_0[i] = 3.0 * g_0_xxyy_0_xxxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyyyz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxyyyzz_0[i] = 3.0 * g_0_xxyy_0_xxxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyyzz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxyyzzz_0[i] = 3.0 * g_0_xxyy_0_xxxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyzzz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxyzzzz_0[i] = 3.0 * g_0_xxyy_0_xxxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyzzzz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxzzzzz_0[i] = g_0_xxxx_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxzzzzz_0[i] * pb_y +
                                     g_0_xxxxy_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxyyyyyy_0[i] = 3.0 * g_0_xxyy_0_xxyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyyyyy_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyyyy_0[i] * pb_x +
                                     g_0_xxxyy_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyyyyyz_0[i] = 3.0 * g_0_xxyy_0_xxyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyyyz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyyyyzz_0[i] = 3.0 * g_0_xxyy_0_xxyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyyzz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyyyzzz_0[i] = 3.0 * g_0_xxyy_0_xxyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyzzz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyyzzzz_0[i] = 3.0 * g_0_xxyy_0_xxyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyzzzz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyzzzzz_0[i] = 3.0 * g_0_xxyy_0_xxyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyzzzzz_0[i] * pb_x +
                                     g_0_xxxyy_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxzzzzzz_0[i] = g_0_xxxx_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxzzzzzz_0[i] * pb_y +
                                     g_0_xxxxy_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xyyyyyyy_0[i] = 3.0 * g_0_xxyy_0_xyyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyyyy_0[i] * pb_x + g_0_xxxyy_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyyyyyz_0[i] = 3.0 * g_0_xxyy_0_xyyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyyyz_0[i] * pb_x + g_0_xxxyy_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyyyyzz_0[i] = 3.0 * g_0_xxyy_0_xyyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyyzz_0[i] * pb_x + g_0_xxxyy_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyyyzzz_0[i] = 3.0 * g_0_xxyy_0_xyyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyzzz_0[i] * pb_x + g_0_xxxyy_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyyzzzz_0[i] = 3.0 * g_0_xxyy_0_xyyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyzzzz_0[i] * pb_x + g_0_xxxyy_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyzzzzz_0[i] = 3.0 * g_0_xxyy_0_xyyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyzzzzz_0[i] * pb_x + g_0_xxxyy_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyzzzzzz_0[i] = 3.0 * g_0_xxyy_0_xyzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyzzzzzz_0[i] * pb_x + g_0_xxxyy_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xzzzzzzz_0[i] = g_0_xxxx_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xzzzzzzz_0[i] * pb_y +
                                     g_0_xxxxy_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_yyyyyyyy_0[i] = 3.0 * g_0_xxyy_0_yyyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyyyyyyy_0[i] * pb_x + g_0_xxxyy_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyyyyyz_0[i] = 3.0 * g_0_xxyy_0_yyyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyyyyyyz_0[i] * pb_x + g_0_xxxyy_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyyyyzz_0[i] = 3.0 * g_0_xxyy_0_yyyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyyyyyzz_0[i] * pb_x + g_0_xxxyy_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyyyzzz_0[i] = 3.0 * g_0_xxyy_0_yyyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyyyyzzz_0[i] * pb_x + g_0_xxxyy_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyyzzzz_0[i] = 3.0 * g_0_xxyy_0_yyyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyyyzzzz_0[i] * pb_x + g_0_xxxyy_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyzzzzz_0[i] = 3.0 * g_0_xxyy_0_yyyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyyzzzzz_0[i] * pb_x + g_0_xxxyy_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyzzzzzz_0[i] = 3.0 * g_0_xxyy_0_yyzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yyzzzzzz_0[i] * pb_x + g_0_xxxyy_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yzzzzzzz_0[i] = 3.0 * g_0_xxyy_0_yzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_yzzzzzzz_0[i] * pb_x + g_0_xxxyy_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_zzzzzzzz_0[i] = 3.0 * g_0_xxyy_0_zzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_zzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_zzzzzzzz_0[i] * pb_x + g_0_xxxyy_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 180-225 components of targeted buffer : SISL

    auto g_0_xxxxyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 180);

    auto g_0_xxxxyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 181);

    auto g_0_xxxxyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 182);

    auto g_0_xxxxyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 183);

    auto g_0_xxxxyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 184);

    auto g_0_xxxxyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 185);

    auto g_0_xxxxyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 186);

    auto g_0_xxxxyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 187);

    auto g_0_xxxxyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 188);

    auto g_0_xxxxyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 189);

    auto g_0_xxxxyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 190);

    auto g_0_xxxxyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 191);

    auto g_0_xxxxyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 192);

    auto g_0_xxxxyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 193);

    auto g_0_xxxxyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 194);

    auto g_0_xxxxyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 195);

    auto g_0_xxxxyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 196);

    auto g_0_xxxxyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 197);

    auto g_0_xxxxyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 198);

    auto g_0_xxxxyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 199);

    auto g_0_xxxxyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 200);

    auto g_0_xxxxyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 201);

    auto g_0_xxxxyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 202);

    auto g_0_xxxxyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 203);

    auto g_0_xxxxyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 204);

    auto g_0_xxxxyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 205);

    auto g_0_xxxxyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 206);

    auto g_0_xxxxyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 207);

    auto g_0_xxxxyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 208);

    auto g_0_xxxxyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 209);

    auto g_0_xxxxyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 210);

    auto g_0_xxxxyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 211);

    auto g_0_xxxxyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 212);

    auto g_0_xxxxyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 213);

    auto g_0_xxxxyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 214);

    auto g_0_xxxxyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 215);

    auto g_0_xxxxyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 216);

    auto g_0_xxxxyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 217);

    auto g_0_xxxxyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 218);

    auto g_0_xxxxyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 219);

    auto g_0_xxxxyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 220);

    auto g_0_xxxxyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 221);

    auto g_0_xxxxyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 222);

    auto g_0_xxxxyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 223);

    auto g_0_xxxxyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 224);

#pragma omp simd aligned(g_0_xxxxy_0_xxxxxxxy_0,      \
                             g_0_xxxxy_0_xxxxxxxy_1,  \
                             g_0_xxxxy_0_xxxxxxyy_0,  \
                             g_0_xxxxy_0_xxxxxxyy_1,  \
                             g_0_xxxxy_0_xxxxxyyy_0,  \
                             g_0_xxxxy_0_xxxxxyyy_1,  \
                             g_0_xxxxy_0_xxxxyyyy_0,  \
                             g_0_xxxxy_0_xxxxyyyy_1,  \
                             g_0_xxxxy_0_xxxyyyyy_0,  \
                             g_0_xxxxy_0_xxxyyyyy_1,  \
                             g_0_xxxxy_0_xxyyyyyy_0,  \
                             g_0_xxxxy_0_xxyyyyyy_1,  \
                             g_0_xxxxy_0_xyyyyyyy_0,  \
                             g_0_xxxxy_0_xyyyyyyy_1,  \
                             g_0_xxxxy_0_yyyyyyyy_0,  \
                             g_0_xxxxy_0_yyyyyyyy_1,  \
                             g_0_xxxxyz_0_xxxxxxxx_0, \
                             g_0_xxxxyz_0_xxxxxxxy_0, \
                             g_0_xxxxyz_0_xxxxxxxz_0, \
                             g_0_xxxxyz_0_xxxxxxyy_0, \
                             g_0_xxxxyz_0_xxxxxxyz_0, \
                             g_0_xxxxyz_0_xxxxxxzz_0, \
                             g_0_xxxxyz_0_xxxxxyyy_0, \
                             g_0_xxxxyz_0_xxxxxyyz_0, \
                             g_0_xxxxyz_0_xxxxxyzz_0, \
                             g_0_xxxxyz_0_xxxxxzzz_0, \
                             g_0_xxxxyz_0_xxxxyyyy_0, \
                             g_0_xxxxyz_0_xxxxyyyz_0, \
                             g_0_xxxxyz_0_xxxxyyzz_0, \
                             g_0_xxxxyz_0_xxxxyzzz_0, \
                             g_0_xxxxyz_0_xxxxzzzz_0, \
                             g_0_xxxxyz_0_xxxyyyyy_0, \
                             g_0_xxxxyz_0_xxxyyyyz_0, \
                             g_0_xxxxyz_0_xxxyyyzz_0, \
                             g_0_xxxxyz_0_xxxyyzzz_0, \
                             g_0_xxxxyz_0_xxxyzzzz_0, \
                             g_0_xxxxyz_0_xxxzzzzz_0, \
                             g_0_xxxxyz_0_xxyyyyyy_0, \
                             g_0_xxxxyz_0_xxyyyyyz_0, \
                             g_0_xxxxyz_0_xxyyyyzz_0, \
                             g_0_xxxxyz_0_xxyyyzzz_0, \
                             g_0_xxxxyz_0_xxyyzzzz_0, \
                             g_0_xxxxyz_0_xxyzzzzz_0, \
                             g_0_xxxxyz_0_xxzzzzzz_0, \
                             g_0_xxxxyz_0_xyyyyyyy_0, \
                             g_0_xxxxyz_0_xyyyyyyz_0, \
                             g_0_xxxxyz_0_xyyyyyzz_0, \
                             g_0_xxxxyz_0_xyyyyzzz_0, \
                             g_0_xxxxyz_0_xyyyzzzz_0, \
                             g_0_xxxxyz_0_xyyzzzzz_0, \
                             g_0_xxxxyz_0_xyzzzzzz_0, \
                             g_0_xxxxyz_0_xzzzzzzz_0, \
                             g_0_xxxxyz_0_yyyyyyyy_0, \
                             g_0_xxxxyz_0_yyyyyyyz_0, \
                             g_0_xxxxyz_0_yyyyyyzz_0, \
                             g_0_xxxxyz_0_yyyyyzzz_0, \
                             g_0_xxxxyz_0_yyyyzzzz_0, \
                             g_0_xxxxyz_0_yyyzzzzz_0, \
                             g_0_xxxxyz_0_yyzzzzzz_0, \
                             g_0_xxxxyz_0_yzzzzzzz_0, \
                             g_0_xxxxyz_0_zzzzzzzz_0, \
                             g_0_xxxxz_0_xxxxxxxx_0,  \
                             g_0_xxxxz_0_xxxxxxxx_1,  \
                             g_0_xxxxz_0_xxxxxxxz_0,  \
                             g_0_xxxxz_0_xxxxxxxz_1,  \
                             g_0_xxxxz_0_xxxxxxyz_0,  \
                             g_0_xxxxz_0_xxxxxxyz_1,  \
                             g_0_xxxxz_0_xxxxxxz_1,   \
                             g_0_xxxxz_0_xxxxxxzz_0,  \
                             g_0_xxxxz_0_xxxxxxzz_1,  \
                             g_0_xxxxz_0_xxxxxyyz_0,  \
                             g_0_xxxxz_0_xxxxxyyz_1,  \
                             g_0_xxxxz_0_xxxxxyz_1,   \
                             g_0_xxxxz_0_xxxxxyzz_0,  \
                             g_0_xxxxz_0_xxxxxyzz_1,  \
                             g_0_xxxxz_0_xxxxxzz_1,   \
                             g_0_xxxxz_0_xxxxxzzz_0,  \
                             g_0_xxxxz_0_xxxxxzzz_1,  \
                             g_0_xxxxz_0_xxxxyyyz_0,  \
                             g_0_xxxxz_0_xxxxyyyz_1,  \
                             g_0_xxxxz_0_xxxxyyz_1,   \
                             g_0_xxxxz_0_xxxxyyzz_0,  \
                             g_0_xxxxz_0_xxxxyyzz_1,  \
                             g_0_xxxxz_0_xxxxyzz_1,   \
                             g_0_xxxxz_0_xxxxyzzz_0,  \
                             g_0_xxxxz_0_xxxxyzzz_1,  \
                             g_0_xxxxz_0_xxxxzzz_1,   \
                             g_0_xxxxz_0_xxxxzzzz_0,  \
                             g_0_xxxxz_0_xxxxzzzz_1,  \
                             g_0_xxxxz_0_xxxyyyyz_0,  \
                             g_0_xxxxz_0_xxxyyyyz_1,  \
                             g_0_xxxxz_0_xxxyyyz_1,   \
                             g_0_xxxxz_0_xxxyyyzz_0,  \
                             g_0_xxxxz_0_xxxyyyzz_1,  \
                             g_0_xxxxz_0_xxxyyzz_1,   \
                             g_0_xxxxz_0_xxxyyzzz_0,  \
                             g_0_xxxxz_0_xxxyyzzz_1,  \
                             g_0_xxxxz_0_xxxyzzz_1,   \
                             g_0_xxxxz_0_xxxyzzzz_0,  \
                             g_0_xxxxz_0_xxxyzzzz_1,  \
                             g_0_xxxxz_0_xxxzzzz_1,   \
                             g_0_xxxxz_0_xxxzzzzz_0,  \
                             g_0_xxxxz_0_xxxzzzzz_1,  \
                             g_0_xxxxz_0_xxyyyyyz_0,  \
                             g_0_xxxxz_0_xxyyyyyz_1,  \
                             g_0_xxxxz_0_xxyyyyz_1,   \
                             g_0_xxxxz_0_xxyyyyzz_0,  \
                             g_0_xxxxz_0_xxyyyyzz_1,  \
                             g_0_xxxxz_0_xxyyyzz_1,   \
                             g_0_xxxxz_0_xxyyyzzz_0,  \
                             g_0_xxxxz_0_xxyyyzzz_1,  \
                             g_0_xxxxz_0_xxyyzzz_1,   \
                             g_0_xxxxz_0_xxyyzzzz_0,  \
                             g_0_xxxxz_0_xxyyzzzz_1,  \
                             g_0_xxxxz_0_xxyzzzz_1,   \
                             g_0_xxxxz_0_xxyzzzzz_0,  \
                             g_0_xxxxz_0_xxyzzzzz_1,  \
                             g_0_xxxxz_0_xxzzzzz_1,   \
                             g_0_xxxxz_0_xxzzzzzz_0,  \
                             g_0_xxxxz_0_xxzzzzzz_1,  \
                             g_0_xxxxz_0_xyyyyyyz_0,  \
                             g_0_xxxxz_0_xyyyyyyz_1,  \
                             g_0_xxxxz_0_xyyyyyz_1,   \
                             g_0_xxxxz_0_xyyyyyzz_0,  \
                             g_0_xxxxz_0_xyyyyyzz_1,  \
                             g_0_xxxxz_0_xyyyyzz_1,   \
                             g_0_xxxxz_0_xyyyyzzz_0,  \
                             g_0_xxxxz_0_xyyyyzzz_1,  \
                             g_0_xxxxz_0_xyyyzzz_1,   \
                             g_0_xxxxz_0_xyyyzzzz_0,  \
                             g_0_xxxxz_0_xyyyzzzz_1,  \
                             g_0_xxxxz_0_xyyzzzz_1,   \
                             g_0_xxxxz_0_xyyzzzzz_0,  \
                             g_0_xxxxz_0_xyyzzzzz_1,  \
                             g_0_xxxxz_0_xyzzzzz_1,   \
                             g_0_xxxxz_0_xyzzzzzz_0,  \
                             g_0_xxxxz_0_xyzzzzzz_1,  \
                             g_0_xxxxz_0_xzzzzzz_1,   \
                             g_0_xxxxz_0_xzzzzzzz_0,  \
                             g_0_xxxxz_0_xzzzzzzz_1,  \
                             g_0_xxxxz_0_yyyyyyyz_0,  \
                             g_0_xxxxz_0_yyyyyyyz_1,  \
                             g_0_xxxxz_0_yyyyyyz_1,   \
                             g_0_xxxxz_0_yyyyyyzz_0,  \
                             g_0_xxxxz_0_yyyyyyzz_1,  \
                             g_0_xxxxz_0_yyyyyzz_1,   \
                             g_0_xxxxz_0_yyyyyzzz_0,  \
                             g_0_xxxxz_0_yyyyyzzz_1,  \
                             g_0_xxxxz_0_yyyyzzz_1,   \
                             g_0_xxxxz_0_yyyyzzzz_0,  \
                             g_0_xxxxz_0_yyyyzzzz_1,  \
                             g_0_xxxxz_0_yyyzzzz_1,   \
                             g_0_xxxxz_0_yyyzzzzz_0,  \
                             g_0_xxxxz_0_yyyzzzzz_1,  \
                             g_0_xxxxz_0_yyzzzzz_1,   \
                             g_0_xxxxz_0_yyzzzzzz_0,  \
                             g_0_xxxxz_0_yyzzzzzz_1,  \
                             g_0_xxxxz_0_yzzzzzz_1,   \
                             g_0_xxxxz_0_yzzzzzzz_0,  \
                             g_0_xxxxz_0_yzzzzzzz_1,  \
                             g_0_xxxxz_0_zzzzzzz_1,   \
                             g_0_xxxxz_0_zzzzzzzz_0,  \
                             g_0_xxxxz_0_zzzzzzzz_1,  \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyz_0_xxxxxxxx_0[i] = g_0_xxxxz_0_xxxxxxxx_0[i] * pb_y + g_0_xxxxz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxxxxy_0[i] = g_0_xxxxy_0_xxxxxxxy_0[i] * pb_z + g_0_xxxxy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxxxxxz_0[i] = g_0_xxxxz_0_xxxxxxxz_0[i] * pb_y + g_0_xxxxz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxxxyy_0[i] = g_0_xxxxy_0_xxxxxxyy_0[i] * pb_z + g_0_xxxxy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxxxxyz_0[i] = g_0_xxxxz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxxxxyz_0[i] * pb_y + g_0_xxxxz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxxxzz_0[i] = g_0_xxxxz_0_xxxxxxzz_0[i] * pb_y + g_0_xxxxz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxxyyy_0[i] = g_0_xxxxy_0_xxxxxyyy_0[i] * pb_z + g_0_xxxxy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxxxyyz_0[i] =
            2.0 * g_0_xxxxz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxxxyyz_0[i] * pb_y + g_0_xxxxz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxxyzz_0[i] = g_0_xxxxz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxxxyzz_0[i] * pb_y + g_0_xxxxz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxxzzz_0[i] = g_0_xxxxz_0_xxxxxzzz_0[i] * pb_y + g_0_xxxxz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxyyyy_0[i] = g_0_xxxxy_0_xxxxyyyy_0[i] * pb_z + g_0_xxxxy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxxyyyz_0[i] =
            3.0 * g_0_xxxxz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxxyyyz_0[i] * pb_y + g_0_xxxxz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxyyzz_0[i] =
            2.0 * g_0_xxxxz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxxyyzz_0[i] * pb_y + g_0_xxxxz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxyzzz_0[i] = g_0_xxxxz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxxyzzz_0[i] * pb_y + g_0_xxxxz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxzzzz_0[i] = g_0_xxxxz_0_xxxxzzzz_0[i] * pb_y + g_0_xxxxz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxyyyyy_0[i] = g_0_xxxxy_0_xxxyyyyy_0[i] * pb_z + g_0_xxxxy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxyyyyz_0[i] =
            4.0 * g_0_xxxxz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxyyyyz_0[i] * pb_y + g_0_xxxxz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxyyyzz_0[i] =
            3.0 * g_0_xxxxz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxyyyzz_0[i] * pb_y + g_0_xxxxz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxyyzzz_0[i] =
            2.0 * g_0_xxxxz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxyyzzz_0[i] * pb_y + g_0_xxxxz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxyzzzz_0[i] = g_0_xxxxz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxyzzzz_0[i] * pb_y + g_0_xxxxz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxzzzzz_0[i] = g_0_xxxxz_0_xxxzzzzz_0[i] * pb_y + g_0_xxxxz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyyyyyy_0[i] = g_0_xxxxy_0_xxyyyyyy_0[i] * pb_z + g_0_xxxxy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxyyyyyz_0[i] =
            5.0 * g_0_xxxxz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyyyyyz_0[i] * pb_y + g_0_xxxxz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyyyyzz_0[i] =
            4.0 * g_0_xxxxz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyyyyzz_0[i] * pb_y + g_0_xxxxz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyyyzzz_0[i] =
            3.0 * g_0_xxxxz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyyyzzz_0[i] * pb_y + g_0_xxxxz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyyzzzz_0[i] =
            2.0 * g_0_xxxxz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyyzzzz_0[i] * pb_y + g_0_xxxxz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyzzzzz_0[i] = g_0_xxxxz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyzzzzz_0[i] * pb_y + g_0_xxxxz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxzzzzzz_0[i] = g_0_xxxxz_0_xxzzzzzz_0[i] * pb_y + g_0_xxxxz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyyyyyy_0[i] = g_0_xxxxy_0_xyyyyyyy_0[i] * pb_z + g_0_xxxxy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xyyyyyyz_0[i] =
            6.0 * g_0_xxxxz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyyyyyz_0[i] * pb_y + g_0_xxxxz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyyyyzz_0[i] =
            5.0 * g_0_xxxxz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyyyyzz_0[i] * pb_y + g_0_xxxxz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyyyzzz_0[i] =
            4.0 * g_0_xxxxz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyyyzzz_0[i] * pb_y + g_0_xxxxz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyyzzzz_0[i] =
            3.0 * g_0_xxxxz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyyzzzz_0[i] * pb_y + g_0_xxxxz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyzzzzz_0[i] =
            2.0 * g_0_xxxxz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyzzzzz_0[i] * pb_y + g_0_xxxxz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyzzzzzz_0[i] = g_0_xxxxz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyzzzzzz_0[i] * pb_y + g_0_xxxxz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xzzzzzzz_0[i] = g_0_xxxxz_0_xzzzzzzz_0[i] * pb_y + g_0_xxxxz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyyyyyy_0[i] = g_0_xxxxy_0_yyyyyyyy_0[i] * pb_z + g_0_xxxxy_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_yyyyyyyz_0[i] =
            7.0 * g_0_xxxxz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyyyyyz_0[i] * pb_y + g_0_xxxxz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyyyyzz_0[i] =
            6.0 * g_0_xxxxz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyyyyzz_0[i] * pb_y + g_0_xxxxz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyyyzzz_0[i] =
            5.0 * g_0_xxxxz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyyyzzz_0[i] * pb_y + g_0_xxxxz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyyzzzz_0[i] =
            4.0 * g_0_xxxxz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyyzzzz_0[i] * pb_y + g_0_xxxxz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyzzzzz_0[i] =
            3.0 * g_0_xxxxz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyzzzzz_0[i] * pb_y + g_0_xxxxz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyzzzzzz_0[i] =
            2.0 * g_0_xxxxz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyzzzzzz_0[i] * pb_y + g_0_xxxxz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yzzzzzzz_0[i] = g_0_xxxxz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yzzzzzzz_0[i] * pb_y + g_0_xxxxz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_zzzzzzzz_0[i] = g_0_xxxxz_0_zzzzzzzz_0[i] * pb_y + g_0_xxxxz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 225-270 components of targeted buffer : SISL

    auto g_0_xxxxzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 225);

    auto g_0_xxxxzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 226);

    auto g_0_xxxxzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 227);

    auto g_0_xxxxzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 228);

    auto g_0_xxxxzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 229);

    auto g_0_xxxxzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 230);

    auto g_0_xxxxzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 231);

    auto g_0_xxxxzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 232);

    auto g_0_xxxxzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 233);

    auto g_0_xxxxzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 234);

    auto g_0_xxxxzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 235);

    auto g_0_xxxxzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 236);

    auto g_0_xxxxzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 237);

    auto g_0_xxxxzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 238);

    auto g_0_xxxxzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 239);

    auto g_0_xxxxzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 240);

    auto g_0_xxxxzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 241);

    auto g_0_xxxxzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 242);

    auto g_0_xxxxzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 243);

    auto g_0_xxxxzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 244);

    auto g_0_xxxxzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 245);

    auto g_0_xxxxzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 246);

    auto g_0_xxxxzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 247);

    auto g_0_xxxxzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 248);

    auto g_0_xxxxzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 249);

    auto g_0_xxxxzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 250);

    auto g_0_xxxxzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 251);

    auto g_0_xxxxzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 252);

    auto g_0_xxxxzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 253);

    auto g_0_xxxxzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 254);

    auto g_0_xxxxzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 255);

    auto g_0_xxxxzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 256);

    auto g_0_xxxxzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 257);

    auto g_0_xxxxzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 258);

    auto g_0_xxxxzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 259);

    auto g_0_xxxxzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 260);

    auto g_0_xxxxzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 261);

    auto g_0_xxxxzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 262);

    auto g_0_xxxxzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 263);

    auto g_0_xxxxzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 264);

    auto g_0_xxxxzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 265);

    auto g_0_xxxxzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 266);

    auto g_0_xxxxzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 267);

    auto g_0_xxxxzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 268);

    auto g_0_xxxxzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 269);

#pragma omp simd aligned(g_0_xxxx_0_xxxxxxxx_0,       \
                             g_0_xxxx_0_xxxxxxxx_1,   \
                             g_0_xxxx_0_xxxxxxxy_0,   \
                             g_0_xxxx_0_xxxxxxxy_1,   \
                             g_0_xxxx_0_xxxxxxyy_0,   \
                             g_0_xxxx_0_xxxxxxyy_1,   \
                             g_0_xxxx_0_xxxxxyyy_0,   \
                             g_0_xxxx_0_xxxxxyyy_1,   \
                             g_0_xxxx_0_xxxxyyyy_0,   \
                             g_0_xxxx_0_xxxxyyyy_1,   \
                             g_0_xxxx_0_xxxyyyyy_0,   \
                             g_0_xxxx_0_xxxyyyyy_1,   \
                             g_0_xxxx_0_xxyyyyyy_0,   \
                             g_0_xxxx_0_xxyyyyyy_1,   \
                             g_0_xxxx_0_xyyyyyyy_0,   \
                             g_0_xxxx_0_xyyyyyyy_1,   \
                             g_0_xxxxz_0_xxxxxxxx_0,  \
                             g_0_xxxxz_0_xxxxxxxx_1,  \
                             g_0_xxxxz_0_xxxxxxxy_0,  \
                             g_0_xxxxz_0_xxxxxxxy_1,  \
                             g_0_xxxxz_0_xxxxxxyy_0,  \
                             g_0_xxxxz_0_xxxxxxyy_1,  \
                             g_0_xxxxz_0_xxxxxyyy_0,  \
                             g_0_xxxxz_0_xxxxxyyy_1,  \
                             g_0_xxxxz_0_xxxxyyyy_0,  \
                             g_0_xxxxz_0_xxxxyyyy_1,  \
                             g_0_xxxxz_0_xxxyyyyy_0,  \
                             g_0_xxxxz_0_xxxyyyyy_1,  \
                             g_0_xxxxz_0_xxyyyyyy_0,  \
                             g_0_xxxxz_0_xxyyyyyy_1,  \
                             g_0_xxxxz_0_xyyyyyyy_0,  \
                             g_0_xxxxz_0_xyyyyyyy_1,  \
                             g_0_xxxxzz_0_xxxxxxxx_0, \
                             g_0_xxxxzz_0_xxxxxxxy_0, \
                             g_0_xxxxzz_0_xxxxxxxz_0, \
                             g_0_xxxxzz_0_xxxxxxyy_0, \
                             g_0_xxxxzz_0_xxxxxxyz_0, \
                             g_0_xxxxzz_0_xxxxxxzz_0, \
                             g_0_xxxxzz_0_xxxxxyyy_0, \
                             g_0_xxxxzz_0_xxxxxyyz_0, \
                             g_0_xxxxzz_0_xxxxxyzz_0, \
                             g_0_xxxxzz_0_xxxxxzzz_0, \
                             g_0_xxxxzz_0_xxxxyyyy_0, \
                             g_0_xxxxzz_0_xxxxyyyz_0, \
                             g_0_xxxxzz_0_xxxxyyzz_0, \
                             g_0_xxxxzz_0_xxxxyzzz_0, \
                             g_0_xxxxzz_0_xxxxzzzz_0, \
                             g_0_xxxxzz_0_xxxyyyyy_0, \
                             g_0_xxxxzz_0_xxxyyyyz_0, \
                             g_0_xxxxzz_0_xxxyyyzz_0, \
                             g_0_xxxxzz_0_xxxyyzzz_0, \
                             g_0_xxxxzz_0_xxxyzzzz_0, \
                             g_0_xxxxzz_0_xxxzzzzz_0, \
                             g_0_xxxxzz_0_xxyyyyyy_0, \
                             g_0_xxxxzz_0_xxyyyyyz_0, \
                             g_0_xxxxzz_0_xxyyyyzz_0, \
                             g_0_xxxxzz_0_xxyyyzzz_0, \
                             g_0_xxxxzz_0_xxyyzzzz_0, \
                             g_0_xxxxzz_0_xxyzzzzz_0, \
                             g_0_xxxxzz_0_xxzzzzzz_0, \
                             g_0_xxxxzz_0_xyyyyyyy_0, \
                             g_0_xxxxzz_0_xyyyyyyz_0, \
                             g_0_xxxxzz_0_xyyyyyzz_0, \
                             g_0_xxxxzz_0_xyyyyzzz_0, \
                             g_0_xxxxzz_0_xyyyzzzz_0, \
                             g_0_xxxxzz_0_xyyzzzzz_0, \
                             g_0_xxxxzz_0_xyzzzzzz_0, \
                             g_0_xxxxzz_0_xzzzzzzz_0, \
                             g_0_xxxxzz_0_yyyyyyyy_0, \
                             g_0_xxxxzz_0_yyyyyyyz_0, \
                             g_0_xxxxzz_0_yyyyyyzz_0, \
                             g_0_xxxxzz_0_yyyyyzzz_0, \
                             g_0_xxxxzz_0_yyyyzzzz_0, \
                             g_0_xxxxzz_0_yyyzzzzz_0, \
                             g_0_xxxxzz_0_yyzzzzzz_0, \
                             g_0_xxxxzz_0_yzzzzzzz_0, \
                             g_0_xxxxzz_0_zzzzzzzz_0, \
                             g_0_xxxzz_0_xxxxxxxz_0,  \
                             g_0_xxxzz_0_xxxxxxxz_1,  \
                             g_0_xxxzz_0_xxxxxxyz_0,  \
                             g_0_xxxzz_0_xxxxxxyz_1,  \
                             g_0_xxxzz_0_xxxxxxz_1,   \
                             g_0_xxxzz_0_xxxxxxzz_0,  \
                             g_0_xxxzz_0_xxxxxxzz_1,  \
                             g_0_xxxzz_0_xxxxxyyz_0,  \
                             g_0_xxxzz_0_xxxxxyyz_1,  \
                             g_0_xxxzz_0_xxxxxyz_1,   \
                             g_0_xxxzz_0_xxxxxyzz_0,  \
                             g_0_xxxzz_0_xxxxxyzz_1,  \
                             g_0_xxxzz_0_xxxxxzz_1,   \
                             g_0_xxxzz_0_xxxxxzzz_0,  \
                             g_0_xxxzz_0_xxxxxzzz_1,  \
                             g_0_xxxzz_0_xxxxyyyz_0,  \
                             g_0_xxxzz_0_xxxxyyyz_1,  \
                             g_0_xxxzz_0_xxxxyyz_1,   \
                             g_0_xxxzz_0_xxxxyyzz_0,  \
                             g_0_xxxzz_0_xxxxyyzz_1,  \
                             g_0_xxxzz_0_xxxxyzz_1,   \
                             g_0_xxxzz_0_xxxxyzzz_0,  \
                             g_0_xxxzz_0_xxxxyzzz_1,  \
                             g_0_xxxzz_0_xxxxzzz_1,   \
                             g_0_xxxzz_0_xxxxzzzz_0,  \
                             g_0_xxxzz_0_xxxxzzzz_1,  \
                             g_0_xxxzz_0_xxxyyyyz_0,  \
                             g_0_xxxzz_0_xxxyyyyz_1,  \
                             g_0_xxxzz_0_xxxyyyz_1,   \
                             g_0_xxxzz_0_xxxyyyzz_0,  \
                             g_0_xxxzz_0_xxxyyyzz_1,  \
                             g_0_xxxzz_0_xxxyyzz_1,   \
                             g_0_xxxzz_0_xxxyyzzz_0,  \
                             g_0_xxxzz_0_xxxyyzzz_1,  \
                             g_0_xxxzz_0_xxxyzzz_1,   \
                             g_0_xxxzz_0_xxxyzzzz_0,  \
                             g_0_xxxzz_0_xxxyzzzz_1,  \
                             g_0_xxxzz_0_xxxzzzz_1,   \
                             g_0_xxxzz_0_xxxzzzzz_0,  \
                             g_0_xxxzz_0_xxxzzzzz_1,  \
                             g_0_xxxzz_0_xxyyyyyz_0,  \
                             g_0_xxxzz_0_xxyyyyyz_1,  \
                             g_0_xxxzz_0_xxyyyyz_1,   \
                             g_0_xxxzz_0_xxyyyyzz_0,  \
                             g_0_xxxzz_0_xxyyyyzz_1,  \
                             g_0_xxxzz_0_xxyyyzz_1,   \
                             g_0_xxxzz_0_xxyyyzzz_0,  \
                             g_0_xxxzz_0_xxyyyzzz_1,  \
                             g_0_xxxzz_0_xxyyzzz_1,   \
                             g_0_xxxzz_0_xxyyzzzz_0,  \
                             g_0_xxxzz_0_xxyyzzzz_1,  \
                             g_0_xxxzz_0_xxyzzzz_1,   \
                             g_0_xxxzz_0_xxyzzzzz_0,  \
                             g_0_xxxzz_0_xxyzzzzz_1,  \
                             g_0_xxxzz_0_xxzzzzz_1,   \
                             g_0_xxxzz_0_xxzzzzzz_0,  \
                             g_0_xxxzz_0_xxzzzzzz_1,  \
                             g_0_xxxzz_0_xyyyyyyz_0,  \
                             g_0_xxxzz_0_xyyyyyyz_1,  \
                             g_0_xxxzz_0_xyyyyyz_1,   \
                             g_0_xxxzz_0_xyyyyyzz_0,  \
                             g_0_xxxzz_0_xyyyyyzz_1,  \
                             g_0_xxxzz_0_xyyyyzz_1,   \
                             g_0_xxxzz_0_xyyyyzzz_0,  \
                             g_0_xxxzz_0_xyyyyzzz_1,  \
                             g_0_xxxzz_0_xyyyzzz_1,   \
                             g_0_xxxzz_0_xyyyzzzz_0,  \
                             g_0_xxxzz_0_xyyyzzzz_1,  \
                             g_0_xxxzz_0_xyyzzzz_1,   \
                             g_0_xxxzz_0_xyyzzzzz_0,  \
                             g_0_xxxzz_0_xyyzzzzz_1,  \
                             g_0_xxxzz_0_xyzzzzz_1,   \
                             g_0_xxxzz_0_xyzzzzzz_0,  \
                             g_0_xxxzz_0_xyzzzzzz_1,  \
                             g_0_xxxzz_0_xzzzzzz_1,   \
                             g_0_xxxzz_0_xzzzzzzz_0,  \
                             g_0_xxxzz_0_xzzzzzzz_1,  \
                             g_0_xxxzz_0_yyyyyyyy_0,  \
                             g_0_xxxzz_0_yyyyyyyy_1,  \
                             g_0_xxxzz_0_yyyyyyyz_0,  \
                             g_0_xxxzz_0_yyyyyyyz_1,  \
                             g_0_xxxzz_0_yyyyyyz_1,   \
                             g_0_xxxzz_0_yyyyyyzz_0,  \
                             g_0_xxxzz_0_yyyyyyzz_1,  \
                             g_0_xxxzz_0_yyyyyzz_1,   \
                             g_0_xxxzz_0_yyyyyzzz_0,  \
                             g_0_xxxzz_0_yyyyyzzz_1,  \
                             g_0_xxxzz_0_yyyyzzz_1,   \
                             g_0_xxxzz_0_yyyyzzzz_0,  \
                             g_0_xxxzz_0_yyyyzzzz_1,  \
                             g_0_xxxzz_0_yyyzzzz_1,   \
                             g_0_xxxzz_0_yyyzzzzz_0,  \
                             g_0_xxxzz_0_yyyzzzzz_1,  \
                             g_0_xxxzz_0_yyzzzzz_1,   \
                             g_0_xxxzz_0_yyzzzzzz_0,  \
                             g_0_xxxzz_0_yyzzzzzz_1,  \
                             g_0_xxxzz_0_yzzzzzz_1,   \
                             g_0_xxxzz_0_yzzzzzzz_0,  \
                             g_0_xxxzz_0_yzzzzzzz_1,  \
                             g_0_xxxzz_0_zzzzzzz_1,   \
                             g_0_xxxzz_0_zzzzzzzz_0,  \
                             g_0_xxxzz_0_zzzzzzzz_1,  \
                             g_0_xxzz_0_xxxxxxxz_0,   \
                             g_0_xxzz_0_xxxxxxxz_1,   \
                             g_0_xxzz_0_xxxxxxyz_0,   \
                             g_0_xxzz_0_xxxxxxyz_1,   \
                             g_0_xxzz_0_xxxxxxzz_0,   \
                             g_0_xxzz_0_xxxxxxzz_1,   \
                             g_0_xxzz_0_xxxxxyyz_0,   \
                             g_0_xxzz_0_xxxxxyyz_1,   \
                             g_0_xxzz_0_xxxxxyzz_0,   \
                             g_0_xxzz_0_xxxxxyzz_1,   \
                             g_0_xxzz_0_xxxxxzzz_0,   \
                             g_0_xxzz_0_xxxxxzzz_1,   \
                             g_0_xxzz_0_xxxxyyyz_0,   \
                             g_0_xxzz_0_xxxxyyyz_1,   \
                             g_0_xxzz_0_xxxxyyzz_0,   \
                             g_0_xxzz_0_xxxxyyzz_1,   \
                             g_0_xxzz_0_xxxxyzzz_0,   \
                             g_0_xxzz_0_xxxxyzzz_1,   \
                             g_0_xxzz_0_xxxxzzzz_0,   \
                             g_0_xxzz_0_xxxxzzzz_1,   \
                             g_0_xxzz_0_xxxyyyyz_0,   \
                             g_0_xxzz_0_xxxyyyyz_1,   \
                             g_0_xxzz_0_xxxyyyzz_0,   \
                             g_0_xxzz_0_xxxyyyzz_1,   \
                             g_0_xxzz_0_xxxyyzzz_0,   \
                             g_0_xxzz_0_xxxyyzzz_1,   \
                             g_0_xxzz_0_xxxyzzzz_0,   \
                             g_0_xxzz_0_xxxyzzzz_1,   \
                             g_0_xxzz_0_xxxzzzzz_0,   \
                             g_0_xxzz_0_xxxzzzzz_1,   \
                             g_0_xxzz_0_xxyyyyyz_0,   \
                             g_0_xxzz_0_xxyyyyyz_1,   \
                             g_0_xxzz_0_xxyyyyzz_0,   \
                             g_0_xxzz_0_xxyyyyzz_1,   \
                             g_0_xxzz_0_xxyyyzzz_0,   \
                             g_0_xxzz_0_xxyyyzzz_1,   \
                             g_0_xxzz_0_xxyyzzzz_0,   \
                             g_0_xxzz_0_xxyyzzzz_1,   \
                             g_0_xxzz_0_xxyzzzzz_0,   \
                             g_0_xxzz_0_xxyzzzzz_1,   \
                             g_0_xxzz_0_xxzzzzzz_0,   \
                             g_0_xxzz_0_xxzzzzzz_1,   \
                             g_0_xxzz_0_xyyyyyyz_0,   \
                             g_0_xxzz_0_xyyyyyyz_1,   \
                             g_0_xxzz_0_xyyyyyzz_0,   \
                             g_0_xxzz_0_xyyyyyzz_1,   \
                             g_0_xxzz_0_xyyyyzzz_0,   \
                             g_0_xxzz_0_xyyyyzzz_1,   \
                             g_0_xxzz_0_xyyyzzzz_0,   \
                             g_0_xxzz_0_xyyyzzzz_1,   \
                             g_0_xxzz_0_xyyzzzzz_0,   \
                             g_0_xxzz_0_xyyzzzzz_1,   \
                             g_0_xxzz_0_xyzzzzzz_0,   \
                             g_0_xxzz_0_xyzzzzzz_1,   \
                             g_0_xxzz_0_xzzzzzzz_0,   \
                             g_0_xxzz_0_xzzzzzzz_1,   \
                             g_0_xxzz_0_yyyyyyyy_0,   \
                             g_0_xxzz_0_yyyyyyyy_1,   \
                             g_0_xxzz_0_yyyyyyyz_0,   \
                             g_0_xxzz_0_yyyyyyyz_1,   \
                             g_0_xxzz_0_yyyyyyzz_0,   \
                             g_0_xxzz_0_yyyyyyzz_1,   \
                             g_0_xxzz_0_yyyyyzzz_0,   \
                             g_0_xxzz_0_yyyyyzzz_1,   \
                             g_0_xxzz_0_yyyyzzzz_0,   \
                             g_0_xxzz_0_yyyyzzzz_1,   \
                             g_0_xxzz_0_yyyzzzzz_0,   \
                             g_0_xxzz_0_yyyzzzzz_1,   \
                             g_0_xxzz_0_yyzzzzzz_0,   \
                             g_0_xxzz_0_yyzzzzzz_1,   \
                             g_0_xxzz_0_yzzzzzzz_0,   \
                             g_0_xxzz_0_yzzzzzzz_1,   \
                             g_0_xxzz_0_zzzzzzzz_0,   \
                             g_0_xxzz_0_zzzzzzzz_1,   \
                             wp_x,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzz_0_xxxxxxxx_0[i] = g_0_xxxx_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxxxxx_0[i] * pb_z +
                                     g_0_xxxxz_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxxxxy_0[i] = g_0_xxxx_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxxxxy_0[i] * pb_z +
                                     g_0_xxxxz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxxxxz_0[i] = 3.0 * g_0_xxzz_0_xxxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxxxz_1[i] * fti_ab_0 +
                                     7.0 * g_0_xxxzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxxxz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxxxyy_0[i] = g_0_xxxx_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxxxyy_0[i] * pb_z +
                                     g_0_xxxxz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxxxyz_0[i] = 3.0 * g_0_xxzz_0_xxxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxxyz_1[i] * fti_ab_0 +
                                     6.0 * g_0_xxxzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxxyz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxxxzz_0[i] = 3.0 * g_0_xxzz_0_xxxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxxzz_1[i] * fti_ab_0 +
                                     6.0 * g_0_xxxzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxxzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxxyyy_0[i] = g_0_xxxx_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxxyyy_0[i] * pb_z +
                                     g_0_xxxxz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxxyyz_0[i] = 3.0 * g_0_xxzz_0_xxxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxyyz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxxyzz_0[i] = 3.0 * g_0_xxzz_0_xxxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxyzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxyzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxxzzz_0[i] = 3.0 * g_0_xxzz_0_xxxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxzzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxyyyy_0[i] = g_0_xxxx_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxyyyy_0[i] * pb_z +
                                     g_0_xxxxz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxyyyz_0[i] = 3.0 * g_0_xxzz_0_xxxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyyyz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxyyzz_0[i] = 3.0 * g_0_xxzz_0_xxxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyyzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxyzzz_0[i] = 3.0 * g_0_xxzz_0_xxxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxyzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyzzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxzzzz_0[i] = 3.0 * g_0_xxzz_0_xxxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxzzzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxyyyyy_0[i] = g_0_xxxx_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxyyyyy_0[i] * pb_z +
                                     g_0_xxxxz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxyyyyz_0[i] = 3.0 * g_0_xxzz_0_xxxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyyyz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxyyyzz_0[i] = 3.0 * g_0_xxzz_0_xxxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyyzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxyyzzz_0[i] = 3.0 * g_0_xxzz_0_xxxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyzzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxyzzzz_0[i] = 3.0 * g_0_xxzz_0_xxxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyzzzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxzzzzz_0[i] = 3.0 * g_0_xxzz_0_xxxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxzzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxzzzzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyyyyyy_0[i] = g_0_xxxx_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxyyyyyy_0[i] * pb_z +
                                     g_0_xxxxz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxyyyyyz_0[i] = 3.0 * g_0_xxzz_0_xxyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyyyz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyyyyzz_0[i] = 3.0 * g_0_xxzz_0_xxyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyyzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyyyzzz_0[i] = 3.0 * g_0_xxzz_0_xxyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyzzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyyzzzz_0[i] = 3.0 * g_0_xxzz_0_xxyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyzzzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyzzzzz_0[i] = 3.0 * g_0_xxzz_0_xxyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyzzzzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxzzzzzz_0[i] = 3.0 * g_0_xxzz_0_xxzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxzzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxzzzzzz_0[i] * pb_x +
                                     g_0_xxxzz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyyyyyy_0[i] = g_0_xxxx_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xyyyyyyy_0[i] * pb_z +
                                     g_0_xxxxz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xyyyyyyz_0[i] = 3.0 * g_0_xxzz_0_xyyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyyyz_0[i] * pb_x + g_0_xxxzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyyyyzz_0[i] = 3.0 * g_0_xxzz_0_xyyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyyzz_0[i] * pb_x + g_0_xxxzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyyyzzz_0[i] = 3.0 * g_0_xxzz_0_xyyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyzzz_0[i] * pb_x + g_0_xxxzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyyzzzz_0[i] = 3.0 * g_0_xxzz_0_xyyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyzzzz_0[i] * pb_x + g_0_xxxzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyzzzzz_0[i] = 3.0 * g_0_xxzz_0_xyyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyzzzzz_0[i] * pb_x + g_0_xxxzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyzzzzzz_0[i] = 3.0 * g_0_xxzz_0_xyzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyzzzzzz_0[i] * pb_x + g_0_xxxzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xzzzzzzz_0[i] = 3.0 * g_0_xxzz_0_xzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xzzzzzzz_0[i] * pb_x + g_0_xxxzz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyyyyy_0[i] = 3.0 * g_0_xxzz_0_yyyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yyyyyyyy_0[i] * pb_x + g_0_xxxzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyyyyz_0[i] = 3.0 * g_0_xxzz_0_yyyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yyyyyyyz_0[i] * pb_x + g_0_xxxzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyyyzz_0[i] = 3.0 * g_0_xxzz_0_yyyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yyyyyyzz_0[i] * pb_x + g_0_xxxzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyyzzz_0[i] = 3.0 * g_0_xxzz_0_yyyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yyyyyzzz_0[i] * pb_x + g_0_xxxzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyzzzz_0[i] = 3.0 * g_0_xxzz_0_yyyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yyyyzzzz_0[i] * pb_x + g_0_xxxzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyzzzzz_0[i] = 3.0 * g_0_xxzz_0_yyyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yyyzzzzz_0[i] * pb_x + g_0_xxxzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyzzzzzz_0[i] = 3.0 * g_0_xxzz_0_yyzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yyzzzzzz_0[i] * pb_x + g_0_xxxzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yzzzzzzz_0[i] = 3.0 * g_0_xxzz_0_yzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_yzzzzzzz_0[i] * pb_x + g_0_xxxzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_zzzzzzzz_0[i] = 3.0 * g_0_xxzz_0_zzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_zzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_zzzzzzzz_0[i] * pb_x + g_0_xxxzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 270-315 components of targeted buffer : SISL

    auto g_0_xxxyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 270);

    auto g_0_xxxyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 271);

    auto g_0_xxxyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 272);

    auto g_0_xxxyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 273);

    auto g_0_xxxyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 274);

    auto g_0_xxxyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 275);

    auto g_0_xxxyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 276);

    auto g_0_xxxyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 277);

    auto g_0_xxxyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 278);

    auto g_0_xxxyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 279);

    auto g_0_xxxyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 280);

    auto g_0_xxxyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 281);

    auto g_0_xxxyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 282);

    auto g_0_xxxyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 283);

    auto g_0_xxxyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 284);

    auto g_0_xxxyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 285);

    auto g_0_xxxyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 286);

    auto g_0_xxxyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 287);

    auto g_0_xxxyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 288);

    auto g_0_xxxyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 289);

    auto g_0_xxxyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 290);

    auto g_0_xxxyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 291);

    auto g_0_xxxyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 292);

    auto g_0_xxxyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 293);

    auto g_0_xxxyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 294);

    auto g_0_xxxyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 295);

    auto g_0_xxxyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 296);

    auto g_0_xxxyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 297);

    auto g_0_xxxyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 298);

    auto g_0_xxxyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 299);

    auto g_0_xxxyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 300);

    auto g_0_xxxyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 301);

    auto g_0_xxxyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 302);

    auto g_0_xxxyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 303);

    auto g_0_xxxyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 304);

    auto g_0_xxxyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 305);

    auto g_0_xxxyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 306);

    auto g_0_xxxyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 307);

    auto g_0_xxxyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 308);

    auto g_0_xxxyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 309);

    auto g_0_xxxyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 310);

    auto g_0_xxxyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 311);

    auto g_0_xxxyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 312);

    auto g_0_xxxyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 313);

    auto g_0_xxxyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 314);

#pragma omp simd aligned(g_0_xxxy_0_xxxxxxxx_0,       \
                             g_0_xxxy_0_xxxxxxxx_1,   \
                             g_0_xxxy_0_xxxxxxxz_0,   \
                             g_0_xxxy_0_xxxxxxxz_1,   \
                             g_0_xxxy_0_xxxxxxzz_0,   \
                             g_0_xxxy_0_xxxxxxzz_1,   \
                             g_0_xxxy_0_xxxxxzzz_0,   \
                             g_0_xxxy_0_xxxxxzzz_1,   \
                             g_0_xxxy_0_xxxxzzzz_0,   \
                             g_0_xxxy_0_xxxxzzzz_1,   \
                             g_0_xxxy_0_xxxzzzzz_0,   \
                             g_0_xxxy_0_xxxzzzzz_1,   \
                             g_0_xxxy_0_xxzzzzzz_0,   \
                             g_0_xxxy_0_xxzzzzzz_1,   \
                             g_0_xxxy_0_xzzzzzzz_0,   \
                             g_0_xxxy_0_xzzzzzzz_1,   \
                             g_0_xxxyy_0_xxxxxxxx_0,  \
                             g_0_xxxyy_0_xxxxxxxx_1,  \
                             g_0_xxxyy_0_xxxxxxxz_0,  \
                             g_0_xxxyy_0_xxxxxxxz_1,  \
                             g_0_xxxyy_0_xxxxxxzz_0,  \
                             g_0_xxxyy_0_xxxxxxzz_1,  \
                             g_0_xxxyy_0_xxxxxzzz_0,  \
                             g_0_xxxyy_0_xxxxxzzz_1,  \
                             g_0_xxxyy_0_xxxxzzzz_0,  \
                             g_0_xxxyy_0_xxxxzzzz_1,  \
                             g_0_xxxyy_0_xxxzzzzz_0,  \
                             g_0_xxxyy_0_xxxzzzzz_1,  \
                             g_0_xxxyy_0_xxzzzzzz_0,  \
                             g_0_xxxyy_0_xxzzzzzz_1,  \
                             g_0_xxxyy_0_xzzzzzzz_0,  \
                             g_0_xxxyy_0_xzzzzzzz_1,  \
                             g_0_xxxyyy_0_xxxxxxxx_0, \
                             g_0_xxxyyy_0_xxxxxxxy_0, \
                             g_0_xxxyyy_0_xxxxxxxz_0, \
                             g_0_xxxyyy_0_xxxxxxyy_0, \
                             g_0_xxxyyy_0_xxxxxxyz_0, \
                             g_0_xxxyyy_0_xxxxxxzz_0, \
                             g_0_xxxyyy_0_xxxxxyyy_0, \
                             g_0_xxxyyy_0_xxxxxyyz_0, \
                             g_0_xxxyyy_0_xxxxxyzz_0, \
                             g_0_xxxyyy_0_xxxxxzzz_0, \
                             g_0_xxxyyy_0_xxxxyyyy_0, \
                             g_0_xxxyyy_0_xxxxyyyz_0, \
                             g_0_xxxyyy_0_xxxxyyzz_0, \
                             g_0_xxxyyy_0_xxxxyzzz_0, \
                             g_0_xxxyyy_0_xxxxzzzz_0, \
                             g_0_xxxyyy_0_xxxyyyyy_0, \
                             g_0_xxxyyy_0_xxxyyyyz_0, \
                             g_0_xxxyyy_0_xxxyyyzz_0, \
                             g_0_xxxyyy_0_xxxyyzzz_0, \
                             g_0_xxxyyy_0_xxxyzzzz_0, \
                             g_0_xxxyyy_0_xxxzzzzz_0, \
                             g_0_xxxyyy_0_xxyyyyyy_0, \
                             g_0_xxxyyy_0_xxyyyyyz_0, \
                             g_0_xxxyyy_0_xxyyyyzz_0, \
                             g_0_xxxyyy_0_xxyyyzzz_0, \
                             g_0_xxxyyy_0_xxyyzzzz_0, \
                             g_0_xxxyyy_0_xxyzzzzz_0, \
                             g_0_xxxyyy_0_xxzzzzzz_0, \
                             g_0_xxxyyy_0_xyyyyyyy_0, \
                             g_0_xxxyyy_0_xyyyyyyz_0, \
                             g_0_xxxyyy_0_xyyyyyzz_0, \
                             g_0_xxxyyy_0_xyyyyzzz_0, \
                             g_0_xxxyyy_0_xyyyzzzz_0, \
                             g_0_xxxyyy_0_xyyzzzzz_0, \
                             g_0_xxxyyy_0_xyzzzzzz_0, \
                             g_0_xxxyyy_0_xzzzzzzz_0, \
                             g_0_xxxyyy_0_yyyyyyyy_0, \
                             g_0_xxxyyy_0_yyyyyyyz_0, \
                             g_0_xxxyyy_0_yyyyyyzz_0, \
                             g_0_xxxyyy_0_yyyyyzzz_0, \
                             g_0_xxxyyy_0_yyyyzzzz_0, \
                             g_0_xxxyyy_0_yyyzzzzz_0, \
                             g_0_xxxyyy_0_yyzzzzzz_0, \
                             g_0_xxxyyy_0_yzzzzzzz_0, \
                             g_0_xxxyyy_0_zzzzzzzz_0, \
                             g_0_xxyyy_0_xxxxxxxy_0,  \
                             g_0_xxyyy_0_xxxxxxxy_1,  \
                             g_0_xxyyy_0_xxxxxxy_1,   \
                             g_0_xxyyy_0_xxxxxxyy_0,  \
                             g_0_xxyyy_0_xxxxxxyy_1,  \
                             g_0_xxyyy_0_xxxxxxyz_0,  \
                             g_0_xxyyy_0_xxxxxxyz_1,  \
                             g_0_xxyyy_0_xxxxxyy_1,   \
                             g_0_xxyyy_0_xxxxxyyy_0,  \
                             g_0_xxyyy_0_xxxxxyyy_1,  \
                             g_0_xxyyy_0_xxxxxyyz_0,  \
                             g_0_xxyyy_0_xxxxxyyz_1,  \
                             g_0_xxyyy_0_xxxxxyz_1,   \
                             g_0_xxyyy_0_xxxxxyzz_0,  \
                             g_0_xxyyy_0_xxxxxyzz_1,  \
                             g_0_xxyyy_0_xxxxyyy_1,   \
                             g_0_xxyyy_0_xxxxyyyy_0,  \
                             g_0_xxyyy_0_xxxxyyyy_1,  \
                             g_0_xxyyy_0_xxxxyyyz_0,  \
                             g_0_xxyyy_0_xxxxyyyz_1,  \
                             g_0_xxyyy_0_xxxxyyz_1,   \
                             g_0_xxyyy_0_xxxxyyzz_0,  \
                             g_0_xxyyy_0_xxxxyyzz_1,  \
                             g_0_xxyyy_0_xxxxyzz_1,   \
                             g_0_xxyyy_0_xxxxyzzz_0,  \
                             g_0_xxyyy_0_xxxxyzzz_1,  \
                             g_0_xxyyy_0_xxxyyyy_1,   \
                             g_0_xxyyy_0_xxxyyyyy_0,  \
                             g_0_xxyyy_0_xxxyyyyy_1,  \
                             g_0_xxyyy_0_xxxyyyyz_0,  \
                             g_0_xxyyy_0_xxxyyyyz_1,  \
                             g_0_xxyyy_0_xxxyyyz_1,   \
                             g_0_xxyyy_0_xxxyyyzz_0,  \
                             g_0_xxyyy_0_xxxyyyzz_1,  \
                             g_0_xxyyy_0_xxxyyzz_1,   \
                             g_0_xxyyy_0_xxxyyzzz_0,  \
                             g_0_xxyyy_0_xxxyyzzz_1,  \
                             g_0_xxyyy_0_xxxyzzz_1,   \
                             g_0_xxyyy_0_xxxyzzzz_0,  \
                             g_0_xxyyy_0_xxxyzzzz_1,  \
                             g_0_xxyyy_0_xxyyyyy_1,   \
                             g_0_xxyyy_0_xxyyyyyy_0,  \
                             g_0_xxyyy_0_xxyyyyyy_1,  \
                             g_0_xxyyy_0_xxyyyyyz_0,  \
                             g_0_xxyyy_0_xxyyyyyz_1,  \
                             g_0_xxyyy_0_xxyyyyz_1,   \
                             g_0_xxyyy_0_xxyyyyzz_0,  \
                             g_0_xxyyy_0_xxyyyyzz_1,  \
                             g_0_xxyyy_0_xxyyyzz_1,   \
                             g_0_xxyyy_0_xxyyyzzz_0,  \
                             g_0_xxyyy_0_xxyyyzzz_1,  \
                             g_0_xxyyy_0_xxyyzzz_1,   \
                             g_0_xxyyy_0_xxyyzzzz_0,  \
                             g_0_xxyyy_0_xxyyzzzz_1,  \
                             g_0_xxyyy_0_xxyzzzz_1,   \
                             g_0_xxyyy_0_xxyzzzzz_0,  \
                             g_0_xxyyy_0_xxyzzzzz_1,  \
                             g_0_xxyyy_0_xyyyyyy_1,   \
                             g_0_xxyyy_0_xyyyyyyy_0,  \
                             g_0_xxyyy_0_xyyyyyyy_1,  \
                             g_0_xxyyy_0_xyyyyyyz_0,  \
                             g_0_xxyyy_0_xyyyyyyz_1,  \
                             g_0_xxyyy_0_xyyyyyz_1,   \
                             g_0_xxyyy_0_xyyyyyzz_0,  \
                             g_0_xxyyy_0_xyyyyyzz_1,  \
                             g_0_xxyyy_0_xyyyyzz_1,   \
                             g_0_xxyyy_0_xyyyyzzz_0,  \
                             g_0_xxyyy_0_xyyyyzzz_1,  \
                             g_0_xxyyy_0_xyyyzzz_1,   \
                             g_0_xxyyy_0_xyyyzzzz_0,  \
                             g_0_xxyyy_0_xyyyzzzz_1,  \
                             g_0_xxyyy_0_xyyzzzz_1,   \
                             g_0_xxyyy_0_xyyzzzzz_0,  \
                             g_0_xxyyy_0_xyyzzzzz_1,  \
                             g_0_xxyyy_0_xyzzzzz_1,   \
                             g_0_xxyyy_0_xyzzzzzz_0,  \
                             g_0_xxyyy_0_xyzzzzzz_1,  \
                             g_0_xxyyy_0_yyyyyyy_1,   \
                             g_0_xxyyy_0_yyyyyyyy_0,  \
                             g_0_xxyyy_0_yyyyyyyy_1,  \
                             g_0_xxyyy_0_yyyyyyyz_0,  \
                             g_0_xxyyy_0_yyyyyyyz_1,  \
                             g_0_xxyyy_0_yyyyyyz_1,   \
                             g_0_xxyyy_0_yyyyyyzz_0,  \
                             g_0_xxyyy_0_yyyyyyzz_1,  \
                             g_0_xxyyy_0_yyyyyzz_1,   \
                             g_0_xxyyy_0_yyyyyzzz_0,  \
                             g_0_xxyyy_0_yyyyyzzz_1,  \
                             g_0_xxyyy_0_yyyyzzz_1,   \
                             g_0_xxyyy_0_yyyyzzzz_0,  \
                             g_0_xxyyy_0_yyyyzzzz_1,  \
                             g_0_xxyyy_0_yyyzzzz_1,   \
                             g_0_xxyyy_0_yyyzzzzz_0,  \
                             g_0_xxyyy_0_yyyzzzzz_1,  \
                             g_0_xxyyy_0_yyzzzzz_1,   \
                             g_0_xxyyy_0_yyzzzzzz_0,  \
                             g_0_xxyyy_0_yyzzzzzz_1,  \
                             g_0_xxyyy_0_yzzzzzz_1,   \
                             g_0_xxyyy_0_yzzzzzzz_0,  \
                             g_0_xxyyy_0_yzzzzzzz_1,  \
                             g_0_xxyyy_0_zzzzzzzz_0,  \
                             g_0_xxyyy_0_zzzzzzzz_1,  \
                             g_0_xyyy_0_xxxxxxxy_0,   \
                             g_0_xyyy_0_xxxxxxxy_1,   \
                             g_0_xyyy_0_xxxxxxyy_0,   \
                             g_0_xyyy_0_xxxxxxyy_1,   \
                             g_0_xyyy_0_xxxxxxyz_0,   \
                             g_0_xyyy_0_xxxxxxyz_1,   \
                             g_0_xyyy_0_xxxxxyyy_0,   \
                             g_0_xyyy_0_xxxxxyyy_1,   \
                             g_0_xyyy_0_xxxxxyyz_0,   \
                             g_0_xyyy_0_xxxxxyyz_1,   \
                             g_0_xyyy_0_xxxxxyzz_0,   \
                             g_0_xyyy_0_xxxxxyzz_1,   \
                             g_0_xyyy_0_xxxxyyyy_0,   \
                             g_0_xyyy_0_xxxxyyyy_1,   \
                             g_0_xyyy_0_xxxxyyyz_0,   \
                             g_0_xyyy_0_xxxxyyyz_1,   \
                             g_0_xyyy_0_xxxxyyzz_0,   \
                             g_0_xyyy_0_xxxxyyzz_1,   \
                             g_0_xyyy_0_xxxxyzzz_0,   \
                             g_0_xyyy_0_xxxxyzzz_1,   \
                             g_0_xyyy_0_xxxyyyyy_0,   \
                             g_0_xyyy_0_xxxyyyyy_1,   \
                             g_0_xyyy_0_xxxyyyyz_0,   \
                             g_0_xyyy_0_xxxyyyyz_1,   \
                             g_0_xyyy_0_xxxyyyzz_0,   \
                             g_0_xyyy_0_xxxyyyzz_1,   \
                             g_0_xyyy_0_xxxyyzzz_0,   \
                             g_0_xyyy_0_xxxyyzzz_1,   \
                             g_0_xyyy_0_xxxyzzzz_0,   \
                             g_0_xyyy_0_xxxyzzzz_1,   \
                             g_0_xyyy_0_xxyyyyyy_0,   \
                             g_0_xyyy_0_xxyyyyyy_1,   \
                             g_0_xyyy_0_xxyyyyyz_0,   \
                             g_0_xyyy_0_xxyyyyyz_1,   \
                             g_0_xyyy_0_xxyyyyzz_0,   \
                             g_0_xyyy_0_xxyyyyzz_1,   \
                             g_0_xyyy_0_xxyyyzzz_0,   \
                             g_0_xyyy_0_xxyyyzzz_1,   \
                             g_0_xyyy_0_xxyyzzzz_0,   \
                             g_0_xyyy_0_xxyyzzzz_1,   \
                             g_0_xyyy_0_xxyzzzzz_0,   \
                             g_0_xyyy_0_xxyzzzzz_1,   \
                             g_0_xyyy_0_xyyyyyyy_0,   \
                             g_0_xyyy_0_xyyyyyyy_1,   \
                             g_0_xyyy_0_xyyyyyyz_0,   \
                             g_0_xyyy_0_xyyyyyyz_1,   \
                             g_0_xyyy_0_xyyyyyzz_0,   \
                             g_0_xyyy_0_xyyyyyzz_1,   \
                             g_0_xyyy_0_xyyyyzzz_0,   \
                             g_0_xyyy_0_xyyyyzzz_1,   \
                             g_0_xyyy_0_xyyyzzzz_0,   \
                             g_0_xyyy_0_xyyyzzzz_1,   \
                             g_0_xyyy_0_xyyzzzzz_0,   \
                             g_0_xyyy_0_xyyzzzzz_1,   \
                             g_0_xyyy_0_xyzzzzzz_0,   \
                             g_0_xyyy_0_xyzzzzzz_1,   \
                             g_0_xyyy_0_yyyyyyyy_0,   \
                             g_0_xyyy_0_yyyyyyyy_1,   \
                             g_0_xyyy_0_yyyyyyyz_0,   \
                             g_0_xyyy_0_yyyyyyyz_1,   \
                             g_0_xyyy_0_yyyyyyzz_0,   \
                             g_0_xyyy_0_yyyyyyzz_1,   \
                             g_0_xyyy_0_yyyyyzzz_0,   \
                             g_0_xyyy_0_yyyyyzzz_1,   \
                             g_0_xyyy_0_yyyyzzzz_0,   \
                             g_0_xyyy_0_yyyyzzzz_1,   \
                             g_0_xyyy_0_yyyzzzzz_0,   \
                             g_0_xyyy_0_yyyzzzzz_1,   \
                             g_0_xyyy_0_yyzzzzzz_0,   \
                             g_0_xyyy_0_yyzzzzzz_1,   \
                             g_0_xyyy_0_yzzzzzzz_0,   \
                             g_0_xyyy_0_yzzzzzzz_1,   \
                             g_0_xyyy_0_zzzzzzzz_0,   \
                             g_0_xyyy_0_zzzzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyy_0_xxxxxxxx_0[i] = 2.0 * g_0_xxxy_0_xxxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_xxxxxxxx_0[i] * pb_y + g_0_xxxyy_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxxxxxy_0[i] = 2.0 * g_0_xyyy_0_xxxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxxxxy_1[i] * fti_ab_0 +
                                     7.0 * g_0_xxyyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxxxy_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxxxxz_0[i] = 2.0 * g_0_xxxy_0_xxxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxxxxz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_xxxxxxxz_0[i] * pb_y + g_0_xxxyy_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxxxxyy_0[i] = 2.0 * g_0_xyyy_0_xxxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxxxyy_1[i] * fti_ab_0 +
                                     6.0 * g_0_xxyyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxxyy_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxxxyz_0[i] = 2.0 * g_0_xyyy_0_xxxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxxxyz_1[i] * fti_ab_0 +
                                     6.0 * g_0_xxyyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxxyz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxxxzz_0[i] = 2.0 * g_0_xxxy_0_xxxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxxxzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_xxxxxxzz_0[i] * pb_y + g_0_xxxyy_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxxxyyy_0[i] = 2.0 * g_0_xyyy_0_xxxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxxyyy_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxyyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxyyy_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxxyyz_0[i] = 2.0 * g_0_xyyy_0_xxxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxyyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxyyz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxxyzz_0[i] = 2.0 * g_0_xyyy_0_xxxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxxyzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxyyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxyzz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxxzzz_0[i] = 2.0 * g_0_xxxy_0_xxxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxxzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_xxxxxzzz_0[i] * pb_y + g_0_xxxyy_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxxyyyy_0[i] = 2.0 * g_0_xyyy_0_xxxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxyyyy_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxyyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyyyy_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxyyyz_0[i] = 2.0 * g_0_xyyy_0_xxxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxyyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyyyz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxyyzz_0[i] = 2.0 * g_0_xyyy_0_xxxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxyyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyyzz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxyzzz_0[i] = 2.0 * g_0_xyyy_0_xxxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxyzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxyyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyzzz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxzzzz_0[i] = 2.0 * g_0_xxxy_0_xxxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_xxxxzzzz_0[i] * pb_y + g_0_xxxyy_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxyyyyy_0[i] = 2.0 * g_0_xyyy_0_xxxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyyyyy_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxyyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyyyy_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxyyyyz_0[i] = 2.0 * g_0_xyyy_0_xxxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxyyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyyyz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxyyyzz_0[i] = 2.0 * g_0_xyyy_0_xxxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxyyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyyzz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxyyzzz_0[i] = 2.0 * g_0_xyyy_0_xxxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxyyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyzzz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxyzzzz_0[i] = 2.0 * g_0_xyyy_0_xxxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxyyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyzzzz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxzzzzz_0[i] = 2.0 * g_0_xxxy_0_xxxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_xxxzzzzz_0[i] * pb_y + g_0_xxxyy_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxyyyyyy_0[i] = 2.0 * g_0_xyyy_0_xxyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyyyyy_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyyyy_0[i] * pb_x +
                                     g_0_xxyyy_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyyyyyz_0[i] = 2.0 * g_0_xyyy_0_xxyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyyyz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyyyyzz_0[i] = 2.0 * g_0_xyyy_0_xxyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyyzz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyyyzzz_0[i] = 2.0 * g_0_xyyy_0_xxyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyzzz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyyzzzz_0[i] = 2.0 * g_0_xyyy_0_xxyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyzzzz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyzzzzz_0[i] = 2.0 * g_0_xyyy_0_xxyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyzzzzz_0[i] * pb_x +
                                     g_0_xxyyy_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxzzzzzz_0[i] = 2.0 * g_0_xxxy_0_xxzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_xxzzzzzz_0[i] * pb_y + g_0_xxxyy_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xyyyyyyy_0[i] = 2.0 * g_0_xyyy_0_xyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyyyy_0[i] * pb_x + g_0_xxyyy_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyyyyyz_0[i] = 2.0 * g_0_xyyy_0_xyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyyyz_0[i] * pb_x + g_0_xxyyy_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyyyyzz_0[i] = 2.0 * g_0_xyyy_0_xyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyyzz_0[i] * pb_x + g_0_xxyyy_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyyyzzz_0[i] = 2.0 * g_0_xyyy_0_xyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyzzz_0[i] * pb_x + g_0_xxyyy_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyyzzzz_0[i] = 2.0 * g_0_xyyy_0_xyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyzzzz_0[i] * pb_x + g_0_xxyyy_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyzzzzz_0[i] = 2.0 * g_0_xyyy_0_xyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyzzzzz_0[i] * pb_x + g_0_xxyyy_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyzzzzzz_0[i] = 2.0 * g_0_xyyy_0_xyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyzzzzzz_0[i] * pb_x + g_0_xxyyy_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xzzzzzzz_0[i] = 2.0 * g_0_xxxy_0_xzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyy_0_xzzzzzzz_0[i] * pb_y + g_0_xxxyy_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_yyyyyyyy_0[i] = 2.0 * g_0_xyyy_0_yyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyyyyyyy_0[i] * pb_x + g_0_xxyyy_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyyyyyz_0[i] = 2.0 * g_0_xyyy_0_yyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyyyyyyz_0[i] * pb_x + g_0_xxyyy_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyyyyzz_0[i] = 2.0 * g_0_xyyy_0_yyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyyyyyzz_0[i] * pb_x + g_0_xxyyy_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyyyzzz_0[i] = 2.0 * g_0_xyyy_0_yyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyyyyzzz_0[i] * pb_x + g_0_xxyyy_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyyzzzz_0[i] = 2.0 * g_0_xyyy_0_yyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyyyzzzz_0[i] * pb_x + g_0_xxyyy_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyzzzzz_0[i] = 2.0 * g_0_xyyy_0_yyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyyzzzzz_0[i] * pb_x + g_0_xxyyy_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyzzzzzz_0[i] = 2.0 * g_0_xyyy_0_yyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yyzzzzzz_0[i] * pb_x + g_0_xxyyy_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yzzzzzzz_0[i] = 2.0 * g_0_xyyy_0_yzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_yzzzzzzz_0[i] * pb_x + g_0_xxyyy_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_zzzzzzzz_0[i] = 2.0 * g_0_xyyy_0_zzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_zzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_zzzzzzzz_0[i] * pb_x + g_0_xxyyy_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 315-360 components of targeted buffer : SISL

    auto g_0_xxxyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 315);

    auto g_0_xxxyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 316);

    auto g_0_xxxyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 317);

    auto g_0_xxxyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 318);

    auto g_0_xxxyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 319);

    auto g_0_xxxyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 320);

    auto g_0_xxxyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 321);

    auto g_0_xxxyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 322);

    auto g_0_xxxyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 323);

    auto g_0_xxxyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 324);

    auto g_0_xxxyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 325);

    auto g_0_xxxyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 326);

    auto g_0_xxxyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 327);

    auto g_0_xxxyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 328);

    auto g_0_xxxyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 329);

    auto g_0_xxxyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 330);

    auto g_0_xxxyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 331);

    auto g_0_xxxyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 332);

    auto g_0_xxxyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 333);

    auto g_0_xxxyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 334);

    auto g_0_xxxyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 335);

    auto g_0_xxxyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 336);

    auto g_0_xxxyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 337);

    auto g_0_xxxyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 338);

    auto g_0_xxxyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 339);

    auto g_0_xxxyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 340);

    auto g_0_xxxyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 341);

    auto g_0_xxxyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 342);

    auto g_0_xxxyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 343);

    auto g_0_xxxyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 344);

    auto g_0_xxxyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 345);

    auto g_0_xxxyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 346);

    auto g_0_xxxyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 347);

    auto g_0_xxxyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 348);

    auto g_0_xxxyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 349);

    auto g_0_xxxyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 350);

    auto g_0_xxxyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 351);

    auto g_0_xxxyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 352);

    auto g_0_xxxyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 353);

    auto g_0_xxxyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 354);

    auto g_0_xxxyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 355);

    auto g_0_xxxyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 356);

    auto g_0_xxxyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 357);

    auto g_0_xxxyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 358);

    auto g_0_xxxyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 359);

#pragma omp simd aligned(g_0_xxxyy_0_xxxxxxx_1,       \
                             g_0_xxxyy_0_xxxxxxxx_0,  \
                             g_0_xxxyy_0_xxxxxxxx_1,  \
                             g_0_xxxyy_0_xxxxxxxy_0,  \
                             g_0_xxxyy_0_xxxxxxxy_1,  \
                             g_0_xxxyy_0_xxxxxxxz_0,  \
                             g_0_xxxyy_0_xxxxxxxz_1,  \
                             g_0_xxxyy_0_xxxxxxy_1,   \
                             g_0_xxxyy_0_xxxxxxyy_0,  \
                             g_0_xxxyy_0_xxxxxxyy_1,  \
                             g_0_xxxyy_0_xxxxxxyz_0,  \
                             g_0_xxxyy_0_xxxxxxyz_1,  \
                             g_0_xxxyy_0_xxxxxxz_1,   \
                             g_0_xxxyy_0_xxxxxxzz_0,  \
                             g_0_xxxyy_0_xxxxxxzz_1,  \
                             g_0_xxxyy_0_xxxxxyy_1,   \
                             g_0_xxxyy_0_xxxxxyyy_0,  \
                             g_0_xxxyy_0_xxxxxyyy_1,  \
                             g_0_xxxyy_0_xxxxxyyz_0,  \
                             g_0_xxxyy_0_xxxxxyyz_1,  \
                             g_0_xxxyy_0_xxxxxyz_1,   \
                             g_0_xxxyy_0_xxxxxyzz_0,  \
                             g_0_xxxyy_0_xxxxxyzz_1,  \
                             g_0_xxxyy_0_xxxxxzz_1,   \
                             g_0_xxxyy_0_xxxxxzzz_0,  \
                             g_0_xxxyy_0_xxxxxzzz_1,  \
                             g_0_xxxyy_0_xxxxyyy_1,   \
                             g_0_xxxyy_0_xxxxyyyy_0,  \
                             g_0_xxxyy_0_xxxxyyyy_1,  \
                             g_0_xxxyy_0_xxxxyyyz_0,  \
                             g_0_xxxyy_0_xxxxyyyz_1,  \
                             g_0_xxxyy_0_xxxxyyz_1,   \
                             g_0_xxxyy_0_xxxxyyzz_0,  \
                             g_0_xxxyy_0_xxxxyyzz_1,  \
                             g_0_xxxyy_0_xxxxyzz_1,   \
                             g_0_xxxyy_0_xxxxyzzz_0,  \
                             g_0_xxxyy_0_xxxxyzzz_1,  \
                             g_0_xxxyy_0_xxxxzzz_1,   \
                             g_0_xxxyy_0_xxxxzzzz_0,  \
                             g_0_xxxyy_0_xxxxzzzz_1,  \
                             g_0_xxxyy_0_xxxyyyy_1,   \
                             g_0_xxxyy_0_xxxyyyyy_0,  \
                             g_0_xxxyy_0_xxxyyyyy_1,  \
                             g_0_xxxyy_0_xxxyyyyz_0,  \
                             g_0_xxxyy_0_xxxyyyyz_1,  \
                             g_0_xxxyy_0_xxxyyyz_1,   \
                             g_0_xxxyy_0_xxxyyyzz_0,  \
                             g_0_xxxyy_0_xxxyyyzz_1,  \
                             g_0_xxxyy_0_xxxyyzz_1,   \
                             g_0_xxxyy_0_xxxyyzzz_0,  \
                             g_0_xxxyy_0_xxxyyzzz_1,  \
                             g_0_xxxyy_0_xxxyzzz_1,   \
                             g_0_xxxyy_0_xxxyzzzz_0,  \
                             g_0_xxxyy_0_xxxyzzzz_1,  \
                             g_0_xxxyy_0_xxxzzzz_1,   \
                             g_0_xxxyy_0_xxxzzzzz_0,  \
                             g_0_xxxyy_0_xxxzzzzz_1,  \
                             g_0_xxxyy_0_xxyyyyy_1,   \
                             g_0_xxxyy_0_xxyyyyyy_0,  \
                             g_0_xxxyy_0_xxyyyyyy_1,  \
                             g_0_xxxyy_0_xxyyyyyz_0,  \
                             g_0_xxxyy_0_xxyyyyyz_1,  \
                             g_0_xxxyy_0_xxyyyyz_1,   \
                             g_0_xxxyy_0_xxyyyyzz_0,  \
                             g_0_xxxyy_0_xxyyyyzz_1,  \
                             g_0_xxxyy_0_xxyyyzz_1,   \
                             g_0_xxxyy_0_xxyyyzzz_0,  \
                             g_0_xxxyy_0_xxyyyzzz_1,  \
                             g_0_xxxyy_0_xxyyzzz_1,   \
                             g_0_xxxyy_0_xxyyzzzz_0,  \
                             g_0_xxxyy_0_xxyyzzzz_1,  \
                             g_0_xxxyy_0_xxyzzzz_1,   \
                             g_0_xxxyy_0_xxyzzzzz_0,  \
                             g_0_xxxyy_0_xxyzzzzz_1,  \
                             g_0_xxxyy_0_xxzzzzz_1,   \
                             g_0_xxxyy_0_xxzzzzzz_0,  \
                             g_0_xxxyy_0_xxzzzzzz_1,  \
                             g_0_xxxyy_0_xyyyyyy_1,   \
                             g_0_xxxyy_0_xyyyyyyy_0,  \
                             g_0_xxxyy_0_xyyyyyyy_1,  \
                             g_0_xxxyy_0_xyyyyyyz_0,  \
                             g_0_xxxyy_0_xyyyyyyz_1,  \
                             g_0_xxxyy_0_xyyyyyz_1,   \
                             g_0_xxxyy_0_xyyyyyzz_0,  \
                             g_0_xxxyy_0_xyyyyyzz_1,  \
                             g_0_xxxyy_0_xyyyyzz_1,   \
                             g_0_xxxyy_0_xyyyyzzz_0,  \
                             g_0_xxxyy_0_xyyyyzzz_1,  \
                             g_0_xxxyy_0_xyyyzzz_1,   \
                             g_0_xxxyy_0_xyyyzzzz_0,  \
                             g_0_xxxyy_0_xyyyzzzz_1,  \
                             g_0_xxxyy_0_xyyzzzz_1,   \
                             g_0_xxxyy_0_xyyzzzzz_0,  \
                             g_0_xxxyy_0_xyyzzzzz_1,  \
                             g_0_xxxyy_0_xyzzzzz_1,   \
                             g_0_xxxyy_0_xyzzzzzz_0,  \
                             g_0_xxxyy_0_xyzzzzzz_1,  \
                             g_0_xxxyy_0_xzzzzzz_1,   \
                             g_0_xxxyy_0_xzzzzzzz_0,  \
                             g_0_xxxyy_0_xzzzzzzz_1,  \
                             g_0_xxxyy_0_yyyyyyy_1,   \
                             g_0_xxxyy_0_yyyyyyyy_0,  \
                             g_0_xxxyy_0_yyyyyyyy_1,  \
                             g_0_xxxyy_0_yyyyyyyz_0,  \
                             g_0_xxxyy_0_yyyyyyyz_1,  \
                             g_0_xxxyy_0_yyyyyyz_1,   \
                             g_0_xxxyy_0_yyyyyyzz_0,  \
                             g_0_xxxyy_0_yyyyyyzz_1,  \
                             g_0_xxxyy_0_yyyyyzz_1,   \
                             g_0_xxxyy_0_yyyyyzzz_0,  \
                             g_0_xxxyy_0_yyyyyzzz_1,  \
                             g_0_xxxyy_0_yyyyzzz_1,   \
                             g_0_xxxyy_0_yyyyzzzz_0,  \
                             g_0_xxxyy_0_yyyyzzzz_1,  \
                             g_0_xxxyy_0_yyyzzzz_1,   \
                             g_0_xxxyy_0_yyyzzzzz_0,  \
                             g_0_xxxyy_0_yyyzzzzz_1,  \
                             g_0_xxxyy_0_yyzzzzz_1,   \
                             g_0_xxxyy_0_yyzzzzzz_0,  \
                             g_0_xxxyy_0_yyzzzzzz_1,  \
                             g_0_xxxyy_0_yzzzzzz_1,   \
                             g_0_xxxyy_0_yzzzzzzz_0,  \
                             g_0_xxxyy_0_yzzzzzzz_1,  \
                             g_0_xxxyy_0_zzzzzzz_1,   \
                             g_0_xxxyy_0_zzzzzzzz_0,  \
                             g_0_xxxyy_0_zzzzzzzz_1,  \
                             g_0_xxxyyz_0_xxxxxxxx_0, \
                             g_0_xxxyyz_0_xxxxxxxy_0, \
                             g_0_xxxyyz_0_xxxxxxxz_0, \
                             g_0_xxxyyz_0_xxxxxxyy_0, \
                             g_0_xxxyyz_0_xxxxxxyz_0, \
                             g_0_xxxyyz_0_xxxxxxzz_0, \
                             g_0_xxxyyz_0_xxxxxyyy_0, \
                             g_0_xxxyyz_0_xxxxxyyz_0, \
                             g_0_xxxyyz_0_xxxxxyzz_0, \
                             g_0_xxxyyz_0_xxxxxzzz_0, \
                             g_0_xxxyyz_0_xxxxyyyy_0, \
                             g_0_xxxyyz_0_xxxxyyyz_0, \
                             g_0_xxxyyz_0_xxxxyyzz_0, \
                             g_0_xxxyyz_0_xxxxyzzz_0, \
                             g_0_xxxyyz_0_xxxxzzzz_0, \
                             g_0_xxxyyz_0_xxxyyyyy_0, \
                             g_0_xxxyyz_0_xxxyyyyz_0, \
                             g_0_xxxyyz_0_xxxyyyzz_0, \
                             g_0_xxxyyz_0_xxxyyzzz_0, \
                             g_0_xxxyyz_0_xxxyzzzz_0, \
                             g_0_xxxyyz_0_xxxzzzzz_0, \
                             g_0_xxxyyz_0_xxyyyyyy_0, \
                             g_0_xxxyyz_0_xxyyyyyz_0, \
                             g_0_xxxyyz_0_xxyyyyzz_0, \
                             g_0_xxxyyz_0_xxyyyzzz_0, \
                             g_0_xxxyyz_0_xxyyzzzz_0, \
                             g_0_xxxyyz_0_xxyzzzzz_0, \
                             g_0_xxxyyz_0_xxzzzzzz_0, \
                             g_0_xxxyyz_0_xyyyyyyy_0, \
                             g_0_xxxyyz_0_xyyyyyyz_0, \
                             g_0_xxxyyz_0_xyyyyyzz_0, \
                             g_0_xxxyyz_0_xyyyyzzz_0, \
                             g_0_xxxyyz_0_xyyyzzzz_0, \
                             g_0_xxxyyz_0_xyyzzzzz_0, \
                             g_0_xxxyyz_0_xyzzzzzz_0, \
                             g_0_xxxyyz_0_xzzzzzzz_0, \
                             g_0_xxxyyz_0_yyyyyyyy_0, \
                             g_0_xxxyyz_0_yyyyyyyz_0, \
                             g_0_xxxyyz_0_yyyyyyzz_0, \
                             g_0_xxxyyz_0_yyyyyzzz_0, \
                             g_0_xxxyyz_0_yyyyzzzz_0, \
                             g_0_xxxyyz_0_yyyzzzzz_0, \
                             g_0_xxxyyz_0_yyzzzzzz_0, \
                             g_0_xxxyyz_0_yzzzzzzz_0, \
                             g_0_xxxyyz_0_zzzzzzzz_0, \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyz_0_xxxxxxxx_0[i] = g_0_xxxyy_0_xxxxxxxx_0[i] * pb_z + g_0_xxxyy_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxxxy_0[i] = g_0_xxxyy_0_xxxxxxxy_0[i] * pb_z + g_0_xxxyy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxxxz_0[i] = g_0_xxxyy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxxxz_0[i] * pb_z + g_0_xxxyy_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxxyy_0[i] = g_0_xxxyy_0_xxxxxxyy_0[i] * pb_z + g_0_xxxyy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxxyz_0[i] = g_0_xxxyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxxyz_0[i] * pb_z + g_0_xxxyy_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxxzz_0[i] =
            2.0 * g_0_xxxyy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxxzz_0[i] * pb_z + g_0_xxxyy_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxyyy_0[i] = g_0_xxxyy_0_xxxxxyyy_0[i] * pb_z + g_0_xxxyy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxyyz_0[i] = g_0_xxxyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxyyz_0[i] * pb_z + g_0_xxxyy_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxyzz_0[i] =
            2.0 * g_0_xxxyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxyzz_0[i] * pb_z + g_0_xxxyy_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxzzz_0[i] =
            3.0 * g_0_xxxyy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxzzz_0[i] * pb_z + g_0_xxxyy_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxyyyy_0[i] = g_0_xxxyy_0_xxxxyyyy_0[i] * pb_z + g_0_xxxyy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxyyyz_0[i] = g_0_xxxyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyyyz_0[i] * pb_z + g_0_xxxyy_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxyyzz_0[i] =
            2.0 * g_0_xxxyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyyzz_0[i] * pb_z + g_0_xxxyy_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxyzzz_0[i] =
            3.0 * g_0_xxxyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyzzz_0[i] * pb_z + g_0_xxxyy_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxzzzz_0[i] =
            4.0 * g_0_xxxyy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxzzzz_0[i] * pb_z + g_0_xxxyy_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyyyyy_0[i] = g_0_xxxyy_0_xxxyyyyy_0[i] * pb_z + g_0_xxxyy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyyyyz_0[i] = g_0_xxxyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyyyz_0[i] * pb_z + g_0_xxxyy_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyyyzz_0[i] =
            2.0 * g_0_xxxyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyyzz_0[i] * pb_z + g_0_xxxyy_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyyzzz_0[i] =
            3.0 * g_0_xxxyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyzzz_0[i] * pb_z + g_0_xxxyy_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyzzzz_0[i] =
            4.0 * g_0_xxxyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyzzzz_0[i] * pb_z + g_0_xxxyy_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxzzzzz_0[i] =
            5.0 * g_0_xxxyy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxzzzzz_0[i] * pb_z + g_0_xxxyy_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyyyyy_0[i] = g_0_xxxyy_0_xxyyyyyy_0[i] * pb_z + g_0_xxxyy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyyyyz_0[i] = g_0_xxxyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyyyz_0[i] * pb_z + g_0_xxxyy_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyyyzz_0[i] =
            2.0 * g_0_xxxyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyyzz_0[i] * pb_z + g_0_xxxyy_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyyzzz_0[i] =
            3.0 * g_0_xxxyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyzzz_0[i] * pb_z + g_0_xxxyy_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyzzzz_0[i] =
            4.0 * g_0_xxxyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyzzzz_0[i] * pb_z + g_0_xxxyy_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyzzzzz_0[i] =
            5.0 * g_0_xxxyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyzzzzz_0[i] * pb_z + g_0_xxxyy_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxzzzzzz_0[i] =
            6.0 * g_0_xxxyy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxzzzzzz_0[i] * pb_z + g_0_xxxyy_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyyyyy_0[i] = g_0_xxxyy_0_xyyyyyyy_0[i] * pb_z + g_0_xxxyy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyyyyz_0[i] = g_0_xxxyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyyyz_0[i] * pb_z + g_0_xxxyy_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyyyzz_0[i] =
            2.0 * g_0_xxxyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyyzz_0[i] * pb_z + g_0_xxxyy_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyyzzz_0[i] =
            3.0 * g_0_xxxyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyzzz_0[i] * pb_z + g_0_xxxyy_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyzzzz_0[i] =
            4.0 * g_0_xxxyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyzzzz_0[i] * pb_z + g_0_xxxyy_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyzzzzz_0[i] =
            5.0 * g_0_xxxyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyzzzzz_0[i] * pb_z + g_0_xxxyy_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyzzzzzz_0[i] =
            6.0 * g_0_xxxyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyzzzzzz_0[i] * pb_z + g_0_xxxyy_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xzzzzzzz_0[i] =
            7.0 * g_0_xxxyy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xzzzzzzz_0[i] * pb_z + g_0_xxxyy_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyyyyy_0[i] = g_0_xxxyy_0_yyyyyyyy_0[i] * pb_z + g_0_xxxyy_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyyyyz_0[i] = g_0_xxxyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyyyyyz_0[i] * pb_z + g_0_xxxyy_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyyyzz_0[i] =
            2.0 * g_0_xxxyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyyyyzz_0[i] * pb_z + g_0_xxxyy_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyyzzz_0[i] =
            3.0 * g_0_xxxyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyyyzzz_0[i] * pb_z + g_0_xxxyy_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyzzzz_0[i] =
            4.0 * g_0_xxxyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyyzzzz_0[i] * pb_z + g_0_xxxyy_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyzzzzz_0[i] =
            5.0 * g_0_xxxyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyzzzzz_0[i] * pb_z + g_0_xxxyy_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyzzzzzz_0[i] =
            6.0 * g_0_xxxyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyzzzzzz_0[i] * pb_z + g_0_xxxyy_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yzzzzzzz_0[i] =
            7.0 * g_0_xxxyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yzzzzzzz_0[i] * pb_z + g_0_xxxyy_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_zzzzzzzz_0[i] =
            8.0 * g_0_xxxyy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_zzzzzzzz_0[i] * pb_z + g_0_xxxyy_0_zzzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 360-405 components of targeted buffer : SISL

    auto g_0_xxxyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 360);

    auto g_0_xxxyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 361);

    auto g_0_xxxyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 362);

    auto g_0_xxxyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 363);

    auto g_0_xxxyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 364);

    auto g_0_xxxyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 365);

    auto g_0_xxxyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 366);

    auto g_0_xxxyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 367);

    auto g_0_xxxyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 368);

    auto g_0_xxxyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 369);

    auto g_0_xxxyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 370);

    auto g_0_xxxyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 371);

    auto g_0_xxxyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 372);

    auto g_0_xxxyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 373);

    auto g_0_xxxyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 374);

    auto g_0_xxxyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 375);

    auto g_0_xxxyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 376);

    auto g_0_xxxyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 377);

    auto g_0_xxxyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 378);

    auto g_0_xxxyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 379);

    auto g_0_xxxyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 380);

    auto g_0_xxxyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 381);

    auto g_0_xxxyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 382);

    auto g_0_xxxyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 383);

    auto g_0_xxxyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 384);

    auto g_0_xxxyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 385);

    auto g_0_xxxyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 386);

    auto g_0_xxxyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 387);

    auto g_0_xxxyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 388);

    auto g_0_xxxyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 389);

    auto g_0_xxxyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 390);

    auto g_0_xxxyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 391);

    auto g_0_xxxyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 392);

    auto g_0_xxxyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 393);

    auto g_0_xxxyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 394);

    auto g_0_xxxyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 395);

    auto g_0_xxxyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 396);

    auto g_0_xxxyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 397);

    auto g_0_xxxyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 398);

    auto g_0_xxxyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 399);

    auto g_0_xxxyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 400);

    auto g_0_xxxyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 401);

    auto g_0_xxxyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 402);

    auto g_0_xxxyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 403);

    auto g_0_xxxyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 404);

#pragma omp simd aligned(g_0_xxxyzz_0_xxxxxxxx_0,     \
                             g_0_xxxyzz_0_xxxxxxxy_0, \
                             g_0_xxxyzz_0_xxxxxxxz_0, \
                             g_0_xxxyzz_0_xxxxxxyy_0, \
                             g_0_xxxyzz_0_xxxxxxyz_0, \
                             g_0_xxxyzz_0_xxxxxxzz_0, \
                             g_0_xxxyzz_0_xxxxxyyy_0, \
                             g_0_xxxyzz_0_xxxxxyyz_0, \
                             g_0_xxxyzz_0_xxxxxyzz_0, \
                             g_0_xxxyzz_0_xxxxxzzz_0, \
                             g_0_xxxyzz_0_xxxxyyyy_0, \
                             g_0_xxxyzz_0_xxxxyyyz_0, \
                             g_0_xxxyzz_0_xxxxyyzz_0, \
                             g_0_xxxyzz_0_xxxxyzzz_0, \
                             g_0_xxxyzz_0_xxxxzzzz_0, \
                             g_0_xxxyzz_0_xxxyyyyy_0, \
                             g_0_xxxyzz_0_xxxyyyyz_0, \
                             g_0_xxxyzz_0_xxxyyyzz_0, \
                             g_0_xxxyzz_0_xxxyyzzz_0, \
                             g_0_xxxyzz_0_xxxyzzzz_0, \
                             g_0_xxxyzz_0_xxxzzzzz_0, \
                             g_0_xxxyzz_0_xxyyyyyy_0, \
                             g_0_xxxyzz_0_xxyyyyyz_0, \
                             g_0_xxxyzz_0_xxyyyyzz_0, \
                             g_0_xxxyzz_0_xxyyyzzz_0, \
                             g_0_xxxyzz_0_xxyyzzzz_0, \
                             g_0_xxxyzz_0_xxyzzzzz_0, \
                             g_0_xxxyzz_0_xxzzzzzz_0, \
                             g_0_xxxyzz_0_xyyyyyyy_0, \
                             g_0_xxxyzz_0_xyyyyyyz_0, \
                             g_0_xxxyzz_0_xyyyyyzz_0, \
                             g_0_xxxyzz_0_xyyyyzzz_0, \
                             g_0_xxxyzz_0_xyyyzzzz_0, \
                             g_0_xxxyzz_0_xyyzzzzz_0, \
                             g_0_xxxyzz_0_xyzzzzzz_0, \
                             g_0_xxxyzz_0_xzzzzzzz_0, \
                             g_0_xxxyzz_0_yyyyyyyy_0, \
                             g_0_xxxyzz_0_yyyyyyyz_0, \
                             g_0_xxxyzz_0_yyyyyyzz_0, \
                             g_0_xxxyzz_0_yyyyyzzz_0, \
                             g_0_xxxyzz_0_yyyyzzzz_0, \
                             g_0_xxxyzz_0_yyyzzzzz_0, \
                             g_0_xxxyzz_0_yyzzzzzz_0, \
                             g_0_xxxyzz_0_yzzzzzzz_0, \
                             g_0_xxxyzz_0_zzzzzzzz_0, \
                             g_0_xxxzz_0_xxxxxxx_1,   \
                             g_0_xxxzz_0_xxxxxxxx_0,  \
                             g_0_xxxzz_0_xxxxxxxx_1,  \
                             g_0_xxxzz_0_xxxxxxxy_0,  \
                             g_0_xxxzz_0_xxxxxxxy_1,  \
                             g_0_xxxzz_0_xxxxxxxz_0,  \
                             g_0_xxxzz_0_xxxxxxxz_1,  \
                             g_0_xxxzz_0_xxxxxxy_1,   \
                             g_0_xxxzz_0_xxxxxxyy_0,  \
                             g_0_xxxzz_0_xxxxxxyy_1,  \
                             g_0_xxxzz_0_xxxxxxyz_0,  \
                             g_0_xxxzz_0_xxxxxxyz_1,  \
                             g_0_xxxzz_0_xxxxxxz_1,   \
                             g_0_xxxzz_0_xxxxxxzz_0,  \
                             g_0_xxxzz_0_xxxxxxzz_1,  \
                             g_0_xxxzz_0_xxxxxyy_1,   \
                             g_0_xxxzz_0_xxxxxyyy_0,  \
                             g_0_xxxzz_0_xxxxxyyy_1,  \
                             g_0_xxxzz_0_xxxxxyyz_0,  \
                             g_0_xxxzz_0_xxxxxyyz_1,  \
                             g_0_xxxzz_0_xxxxxyz_1,   \
                             g_0_xxxzz_0_xxxxxyzz_0,  \
                             g_0_xxxzz_0_xxxxxyzz_1,  \
                             g_0_xxxzz_0_xxxxxzz_1,   \
                             g_0_xxxzz_0_xxxxxzzz_0,  \
                             g_0_xxxzz_0_xxxxxzzz_1,  \
                             g_0_xxxzz_0_xxxxyyy_1,   \
                             g_0_xxxzz_0_xxxxyyyy_0,  \
                             g_0_xxxzz_0_xxxxyyyy_1,  \
                             g_0_xxxzz_0_xxxxyyyz_0,  \
                             g_0_xxxzz_0_xxxxyyyz_1,  \
                             g_0_xxxzz_0_xxxxyyz_1,   \
                             g_0_xxxzz_0_xxxxyyzz_0,  \
                             g_0_xxxzz_0_xxxxyyzz_1,  \
                             g_0_xxxzz_0_xxxxyzz_1,   \
                             g_0_xxxzz_0_xxxxyzzz_0,  \
                             g_0_xxxzz_0_xxxxyzzz_1,  \
                             g_0_xxxzz_0_xxxxzzz_1,   \
                             g_0_xxxzz_0_xxxxzzzz_0,  \
                             g_0_xxxzz_0_xxxxzzzz_1,  \
                             g_0_xxxzz_0_xxxyyyy_1,   \
                             g_0_xxxzz_0_xxxyyyyy_0,  \
                             g_0_xxxzz_0_xxxyyyyy_1,  \
                             g_0_xxxzz_0_xxxyyyyz_0,  \
                             g_0_xxxzz_0_xxxyyyyz_1,  \
                             g_0_xxxzz_0_xxxyyyz_1,   \
                             g_0_xxxzz_0_xxxyyyzz_0,  \
                             g_0_xxxzz_0_xxxyyyzz_1,  \
                             g_0_xxxzz_0_xxxyyzz_1,   \
                             g_0_xxxzz_0_xxxyyzzz_0,  \
                             g_0_xxxzz_0_xxxyyzzz_1,  \
                             g_0_xxxzz_0_xxxyzzz_1,   \
                             g_0_xxxzz_0_xxxyzzzz_0,  \
                             g_0_xxxzz_0_xxxyzzzz_1,  \
                             g_0_xxxzz_0_xxxzzzz_1,   \
                             g_0_xxxzz_0_xxxzzzzz_0,  \
                             g_0_xxxzz_0_xxxzzzzz_1,  \
                             g_0_xxxzz_0_xxyyyyy_1,   \
                             g_0_xxxzz_0_xxyyyyyy_0,  \
                             g_0_xxxzz_0_xxyyyyyy_1,  \
                             g_0_xxxzz_0_xxyyyyyz_0,  \
                             g_0_xxxzz_0_xxyyyyyz_1,  \
                             g_0_xxxzz_0_xxyyyyz_1,   \
                             g_0_xxxzz_0_xxyyyyzz_0,  \
                             g_0_xxxzz_0_xxyyyyzz_1,  \
                             g_0_xxxzz_0_xxyyyzz_1,   \
                             g_0_xxxzz_0_xxyyyzzz_0,  \
                             g_0_xxxzz_0_xxyyyzzz_1,  \
                             g_0_xxxzz_0_xxyyzzz_1,   \
                             g_0_xxxzz_0_xxyyzzzz_0,  \
                             g_0_xxxzz_0_xxyyzzzz_1,  \
                             g_0_xxxzz_0_xxyzzzz_1,   \
                             g_0_xxxzz_0_xxyzzzzz_0,  \
                             g_0_xxxzz_0_xxyzzzzz_1,  \
                             g_0_xxxzz_0_xxzzzzz_1,   \
                             g_0_xxxzz_0_xxzzzzzz_0,  \
                             g_0_xxxzz_0_xxzzzzzz_1,  \
                             g_0_xxxzz_0_xyyyyyy_1,   \
                             g_0_xxxzz_0_xyyyyyyy_0,  \
                             g_0_xxxzz_0_xyyyyyyy_1,  \
                             g_0_xxxzz_0_xyyyyyyz_0,  \
                             g_0_xxxzz_0_xyyyyyyz_1,  \
                             g_0_xxxzz_0_xyyyyyz_1,   \
                             g_0_xxxzz_0_xyyyyyzz_0,  \
                             g_0_xxxzz_0_xyyyyyzz_1,  \
                             g_0_xxxzz_0_xyyyyzz_1,   \
                             g_0_xxxzz_0_xyyyyzzz_0,  \
                             g_0_xxxzz_0_xyyyyzzz_1,  \
                             g_0_xxxzz_0_xyyyzzz_1,   \
                             g_0_xxxzz_0_xyyyzzzz_0,  \
                             g_0_xxxzz_0_xyyyzzzz_1,  \
                             g_0_xxxzz_0_xyyzzzz_1,   \
                             g_0_xxxzz_0_xyyzzzzz_0,  \
                             g_0_xxxzz_0_xyyzzzzz_1,  \
                             g_0_xxxzz_0_xyzzzzz_1,   \
                             g_0_xxxzz_0_xyzzzzzz_0,  \
                             g_0_xxxzz_0_xyzzzzzz_1,  \
                             g_0_xxxzz_0_xzzzzzz_1,   \
                             g_0_xxxzz_0_xzzzzzzz_0,  \
                             g_0_xxxzz_0_xzzzzzzz_1,  \
                             g_0_xxxzz_0_yyyyyyy_1,   \
                             g_0_xxxzz_0_yyyyyyyy_0,  \
                             g_0_xxxzz_0_yyyyyyyy_1,  \
                             g_0_xxxzz_0_yyyyyyyz_0,  \
                             g_0_xxxzz_0_yyyyyyyz_1,  \
                             g_0_xxxzz_0_yyyyyyz_1,   \
                             g_0_xxxzz_0_yyyyyyzz_0,  \
                             g_0_xxxzz_0_yyyyyyzz_1,  \
                             g_0_xxxzz_0_yyyyyzz_1,   \
                             g_0_xxxzz_0_yyyyyzzz_0,  \
                             g_0_xxxzz_0_yyyyyzzz_1,  \
                             g_0_xxxzz_0_yyyyzzz_1,   \
                             g_0_xxxzz_0_yyyyzzzz_0,  \
                             g_0_xxxzz_0_yyyyzzzz_1,  \
                             g_0_xxxzz_0_yyyzzzz_1,   \
                             g_0_xxxzz_0_yyyzzzzz_0,  \
                             g_0_xxxzz_0_yyyzzzzz_1,  \
                             g_0_xxxzz_0_yyzzzzz_1,   \
                             g_0_xxxzz_0_yyzzzzzz_0,  \
                             g_0_xxxzz_0_yyzzzzzz_1,  \
                             g_0_xxxzz_0_yzzzzzz_1,   \
                             g_0_xxxzz_0_yzzzzzzz_0,  \
                             g_0_xxxzz_0_yzzzzzzz_1,  \
                             g_0_xxxzz_0_zzzzzzz_1,   \
                             g_0_xxxzz_0_zzzzzzzz_0,  \
                             g_0_xxxzz_0_zzzzzzzz_1,  \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzz_0_xxxxxxxx_0[i] = g_0_xxxzz_0_xxxxxxxx_0[i] * pb_y + g_0_xxxzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxxxy_0[i] = g_0_xxxzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxxxy_0[i] * pb_y + g_0_xxxzz_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxxxz_0[i] = g_0_xxxzz_0_xxxxxxxz_0[i] * pb_y + g_0_xxxzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxxyy_0[i] =
            2.0 * g_0_xxxzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxxyy_0[i] * pb_y + g_0_xxxzz_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxxyz_0[i] = g_0_xxxzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxxyz_0[i] * pb_y + g_0_xxxzz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxxzz_0[i] = g_0_xxxzz_0_xxxxxxzz_0[i] * pb_y + g_0_xxxzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxyyy_0[i] =
            3.0 * g_0_xxxzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxyyy_0[i] * pb_y + g_0_xxxzz_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxyyz_0[i] =
            2.0 * g_0_xxxzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxyyz_0[i] * pb_y + g_0_xxxzz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxyzz_0[i] = g_0_xxxzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxyzz_0[i] * pb_y + g_0_xxxzz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxzzz_0[i] = g_0_xxxzz_0_xxxxxzzz_0[i] * pb_y + g_0_xxxzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxyyyy_0[i] =
            4.0 * g_0_xxxzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyyyy_0[i] * pb_y + g_0_xxxzz_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxyyyz_0[i] =
            3.0 * g_0_xxxzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyyyz_0[i] * pb_y + g_0_xxxzz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxyyzz_0[i] =
            2.0 * g_0_xxxzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyyzz_0[i] * pb_y + g_0_xxxzz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxyzzz_0[i] = g_0_xxxzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyzzz_0[i] * pb_y + g_0_xxxzz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxzzzz_0[i] = g_0_xxxzz_0_xxxxzzzz_0[i] * pb_y + g_0_xxxzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyyyyy_0[i] =
            5.0 * g_0_xxxzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyyyy_0[i] * pb_y + g_0_xxxzz_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyyyyz_0[i] =
            4.0 * g_0_xxxzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyyyz_0[i] * pb_y + g_0_xxxzz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyyyzz_0[i] =
            3.0 * g_0_xxxzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyyzz_0[i] * pb_y + g_0_xxxzz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyyzzz_0[i] =
            2.0 * g_0_xxxzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyzzz_0[i] * pb_y + g_0_xxxzz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyzzzz_0[i] = g_0_xxxzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyzzzz_0[i] * pb_y + g_0_xxxzz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxzzzzz_0[i] = g_0_xxxzz_0_xxxzzzzz_0[i] * pb_y + g_0_xxxzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyyyyy_0[i] =
            6.0 * g_0_xxxzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyyyy_0[i] * pb_y + g_0_xxxzz_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyyyyz_0[i] =
            5.0 * g_0_xxxzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyyyz_0[i] * pb_y + g_0_xxxzz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyyyzz_0[i] =
            4.0 * g_0_xxxzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyyzz_0[i] * pb_y + g_0_xxxzz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyyzzz_0[i] =
            3.0 * g_0_xxxzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyzzz_0[i] * pb_y + g_0_xxxzz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyzzzz_0[i] =
            2.0 * g_0_xxxzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyzzzz_0[i] * pb_y + g_0_xxxzz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyzzzzz_0[i] = g_0_xxxzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyzzzzz_0[i] * pb_y + g_0_xxxzz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxzzzzzz_0[i] = g_0_xxxzz_0_xxzzzzzz_0[i] * pb_y + g_0_xxxzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyyyyy_0[i] =
            7.0 * g_0_xxxzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyyyy_0[i] * pb_y + g_0_xxxzz_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyyyyz_0[i] =
            6.0 * g_0_xxxzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyyyz_0[i] * pb_y + g_0_xxxzz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyyyzz_0[i] =
            5.0 * g_0_xxxzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyyzz_0[i] * pb_y + g_0_xxxzz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyyzzz_0[i] =
            4.0 * g_0_xxxzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyzzz_0[i] * pb_y + g_0_xxxzz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyzzzz_0[i] =
            3.0 * g_0_xxxzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyzzzz_0[i] * pb_y + g_0_xxxzz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyzzzzz_0[i] =
            2.0 * g_0_xxxzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyzzzzz_0[i] * pb_y + g_0_xxxzz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyzzzzzz_0[i] = g_0_xxxzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyzzzzzz_0[i] * pb_y + g_0_xxxzz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xzzzzzzz_0[i] = g_0_xxxzz_0_xzzzzzzz_0[i] * pb_y + g_0_xxxzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyyyyy_0[i] =
            8.0 * g_0_xxxzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyyyyy_0[i] * pb_y + g_0_xxxzz_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyyyyz_0[i] =
            7.0 * g_0_xxxzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyyyyz_0[i] * pb_y + g_0_xxxzz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyyyzz_0[i] =
            6.0 * g_0_xxxzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyyyzz_0[i] * pb_y + g_0_xxxzz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyyzzz_0[i] =
            5.0 * g_0_xxxzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyyzzz_0[i] * pb_y + g_0_xxxzz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyzzzz_0[i] =
            4.0 * g_0_xxxzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyzzzz_0[i] * pb_y + g_0_xxxzz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyzzzzz_0[i] =
            3.0 * g_0_xxxzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyzzzzz_0[i] * pb_y + g_0_xxxzz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyzzzzzz_0[i] =
            2.0 * g_0_xxxzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyzzzzzz_0[i] * pb_y + g_0_xxxzz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yzzzzzzz_0[i] = g_0_xxxzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yzzzzzzz_0[i] * pb_y + g_0_xxxzz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_zzzzzzzz_0[i] = g_0_xxxzz_0_zzzzzzzz_0[i] * pb_y + g_0_xxxzz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 405-450 components of targeted buffer : SISL

    auto g_0_xxxzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 405);

    auto g_0_xxxzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 406);

    auto g_0_xxxzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 407);

    auto g_0_xxxzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 408);

    auto g_0_xxxzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 409);

    auto g_0_xxxzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 410);

    auto g_0_xxxzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 411);

    auto g_0_xxxzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 412);

    auto g_0_xxxzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 413);

    auto g_0_xxxzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 414);

    auto g_0_xxxzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 415);

    auto g_0_xxxzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 416);

    auto g_0_xxxzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 417);

    auto g_0_xxxzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 418);

    auto g_0_xxxzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 419);

    auto g_0_xxxzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 420);

    auto g_0_xxxzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 421);

    auto g_0_xxxzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 422);

    auto g_0_xxxzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 423);

    auto g_0_xxxzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 424);

    auto g_0_xxxzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 425);

    auto g_0_xxxzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 426);

    auto g_0_xxxzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 427);

    auto g_0_xxxzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 428);

    auto g_0_xxxzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 429);

    auto g_0_xxxzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 430);

    auto g_0_xxxzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 431);

    auto g_0_xxxzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 432);

    auto g_0_xxxzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 433);

    auto g_0_xxxzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 434);

    auto g_0_xxxzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 435);

    auto g_0_xxxzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 436);

    auto g_0_xxxzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 437);

    auto g_0_xxxzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 438);

    auto g_0_xxxzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 439);

    auto g_0_xxxzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 440);

    auto g_0_xxxzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 441);

    auto g_0_xxxzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 442);

    auto g_0_xxxzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 443);

    auto g_0_xxxzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 444);

    auto g_0_xxxzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 445);

    auto g_0_xxxzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 446);

    auto g_0_xxxzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 447);

    auto g_0_xxxzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 448);

    auto g_0_xxxzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 449);

#pragma omp simd aligned(g_0_xxxz_0_xxxxxxxx_0,       \
                             g_0_xxxz_0_xxxxxxxx_1,   \
                             g_0_xxxz_0_xxxxxxxy_0,   \
                             g_0_xxxz_0_xxxxxxxy_1,   \
                             g_0_xxxz_0_xxxxxxyy_0,   \
                             g_0_xxxz_0_xxxxxxyy_1,   \
                             g_0_xxxz_0_xxxxxyyy_0,   \
                             g_0_xxxz_0_xxxxxyyy_1,   \
                             g_0_xxxz_0_xxxxyyyy_0,   \
                             g_0_xxxz_0_xxxxyyyy_1,   \
                             g_0_xxxz_0_xxxyyyyy_0,   \
                             g_0_xxxz_0_xxxyyyyy_1,   \
                             g_0_xxxz_0_xxyyyyyy_0,   \
                             g_0_xxxz_0_xxyyyyyy_1,   \
                             g_0_xxxz_0_xyyyyyyy_0,   \
                             g_0_xxxz_0_xyyyyyyy_1,   \
                             g_0_xxxzz_0_xxxxxxxx_0,  \
                             g_0_xxxzz_0_xxxxxxxx_1,  \
                             g_0_xxxzz_0_xxxxxxxy_0,  \
                             g_0_xxxzz_0_xxxxxxxy_1,  \
                             g_0_xxxzz_0_xxxxxxyy_0,  \
                             g_0_xxxzz_0_xxxxxxyy_1,  \
                             g_0_xxxzz_0_xxxxxyyy_0,  \
                             g_0_xxxzz_0_xxxxxyyy_1,  \
                             g_0_xxxzz_0_xxxxyyyy_0,  \
                             g_0_xxxzz_0_xxxxyyyy_1,  \
                             g_0_xxxzz_0_xxxyyyyy_0,  \
                             g_0_xxxzz_0_xxxyyyyy_1,  \
                             g_0_xxxzz_0_xxyyyyyy_0,  \
                             g_0_xxxzz_0_xxyyyyyy_1,  \
                             g_0_xxxzz_0_xyyyyyyy_0,  \
                             g_0_xxxzz_0_xyyyyyyy_1,  \
                             g_0_xxxzzz_0_xxxxxxxx_0, \
                             g_0_xxxzzz_0_xxxxxxxy_0, \
                             g_0_xxxzzz_0_xxxxxxxz_0, \
                             g_0_xxxzzz_0_xxxxxxyy_0, \
                             g_0_xxxzzz_0_xxxxxxyz_0, \
                             g_0_xxxzzz_0_xxxxxxzz_0, \
                             g_0_xxxzzz_0_xxxxxyyy_0, \
                             g_0_xxxzzz_0_xxxxxyyz_0, \
                             g_0_xxxzzz_0_xxxxxyzz_0, \
                             g_0_xxxzzz_0_xxxxxzzz_0, \
                             g_0_xxxzzz_0_xxxxyyyy_0, \
                             g_0_xxxzzz_0_xxxxyyyz_0, \
                             g_0_xxxzzz_0_xxxxyyzz_0, \
                             g_0_xxxzzz_0_xxxxyzzz_0, \
                             g_0_xxxzzz_0_xxxxzzzz_0, \
                             g_0_xxxzzz_0_xxxyyyyy_0, \
                             g_0_xxxzzz_0_xxxyyyyz_0, \
                             g_0_xxxzzz_0_xxxyyyzz_0, \
                             g_0_xxxzzz_0_xxxyyzzz_0, \
                             g_0_xxxzzz_0_xxxyzzzz_0, \
                             g_0_xxxzzz_0_xxxzzzzz_0, \
                             g_0_xxxzzz_0_xxyyyyyy_0, \
                             g_0_xxxzzz_0_xxyyyyyz_0, \
                             g_0_xxxzzz_0_xxyyyyzz_0, \
                             g_0_xxxzzz_0_xxyyyzzz_0, \
                             g_0_xxxzzz_0_xxyyzzzz_0, \
                             g_0_xxxzzz_0_xxyzzzzz_0, \
                             g_0_xxxzzz_0_xxzzzzzz_0, \
                             g_0_xxxzzz_0_xyyyyyyy_0, \
                             g_0_xxxzzz_0_xyyyyyyz_0, \
                             g_0_xxxzzz_0_xyyyyyzz_0, \
                             g_0_xxxzzz_0_xyyyyzzz_0, \
                             g_0_xxxzzz_0_xyyyzzzz_0, \
                             g_0_xxxzzz_0_xyyzzzzz_0, \
                             g_0_xxxzzz_0_xyzzzzzz_0, \
                             g_0_xxxzzz_0_xzzzzzzz_0, \
                             g_0_xxxzzz_0_yyyyyyyy_0, \
                             g_0_xxxzzz_0_yyyyyyyz_0, \
                             g_0_xxxzzz_0_yyyyyyzz_0, \
                             g_0_xxxzzz_0_yyyyyzzz_0, \
                             g_0_xxxzzz_0_yyyyzzzz_0, \
                             g_0_xxxzzz_0_yyyzzzzz_0, \
                             g_0_xxxzzz_0_yyzzzzzz_0, \
                             g_0_xxxzzz_0_yzzzzzzz_0, \
                             g_0_xxxzzz_0_zzzzzzzz_0, \
                             g_0_xxzzz_0_xxxxxxxz_0,  \
                             g_0_xxzzz_0_xxxxxxxz_1,  \
                             g_0_xxzzz_0_xxxxxxyz_0,  \
                             g_0_xxzzz_0_xxxxxxyz_1,  \
                             g_0_xxzzz_0_xxxxxxz_1,   \
                             g_0_xxzzz_0_xxxxxxzz_0,  \
                             g_0_xxzzz_0_xxxxxxzz_1,  \
                             g_0_xxzzz_0_xxxxxyyz_0,  \
                             g_0_xxzzz_0_xxxxxyyz_1,  \
                             g_0_xxzzz_0_xxxxxyz_1,   \
                             g_0_xxzzz_0_xxxxxyzz_0,  \
                             g_0_xxzzz_0_xxxxxyzz_1,  \
                             g_0_xxzzz_0_xxxxxzz_1,   \
                             g_0_xxzzz_0_xxxxxzzz_0,  \
                             g_0_xxzzz_0_xxxxxzzz_1,  \
                             g_0_xxzzz_0_xxxxyyyz_0,  \
                             g_0_xxzzz_0_xxxxyyyz_1,  \
                             g_0_xxzzz_0_xxxxyyz_1,   \
                             g_0_xxzzz_0_xxxxyyzz_0,  \
                             g_0_xxzzz_0_xxxxyyzz_1,  \
                             g_0_xxzzz_0_xxxxyzz_1,   \
                             g_0_xxzzz_0_xxxxyzzz_0,  \
                             g_0_xxzzz_0_xxxxyzzz_1,  \
                             g_0_xxzzz_0_xxxxzzz_1,   \
                             g_0_xxzzz_0_xxxxzzzz_0,  \
                             g_0_xxzzz_0_xxxxzzzz_1,  \
                             g_0_xxzzz_0_xxxyyyyz_0,  \
                             g_0_xxzzz_0_xxxyyyyz_1,  \
                             g_0_xxzzz_0_xxxyyyz_1,   \
                             g_0_xxzzz_0_xxxyyyzz_0,  \
                             g_0_xxzzz_0_xxxyyyzz_1,  \
                             g_0_xxzzz_0_xxxyyzz_1,   \
                             g_0_xxzzz_0_xxxyyzzz_0,  \
                             g_0_xxzzz_0_xxxyyzzz_1,  \
                             g_0_xxzzz_0_xxxyzzz_1,   \
                             g_0_xxzzz_0_xxxyzzzz_0,  \
                             g_0_xxzzz_0_xxxyzzzz_1,  \
                             g_0_xxzzz_0_xxxzzzz_1,   \
                             g_0_xxzzz_0_xxxzzzzz_0,  \
                             g_0_xxzzz_0_xxxzzzzz_1,  \
                             g_0_xxzzz_0_xxyyyyyz_0,  \
                             g_0_xxzzz_0_xxyyyyyz_1,  \
                             g_0_xxzzz_0_xxyyyyz_1,   \
                             g_0_xxzzz_0_xxyyyyzz_0,  \
                             g_0_xxzzz_0_xxyyyyzz_1,  \
                             g_0_xxzzz_0_xxyyyzz_1,   \
                             g_0_xxzzz_0_xxyyyzzz_0,  \
                             g_0_xxzzz_0_xxyyyzzz_1,  \
                             g_0_xxzzz_0_xxyyzzz_1,   \
                             g_0_xxzzz_0_xxyyzzzz_0,  \
                             g_0_xxzzz_0_xxyyzzzz_1,  \
                             g_0_xxzzz_0_xxyzzzz_1,   \
                             g_0_xxzzz_0_xxyzzzzz_0,  \
                             g_0_xxzzz_0_xxyzzzzz_1,  \
                             g_0_xxzzz_0_xxzzzzz_1,   \
                             g_0_xxzzz_0_xxzzzzzz_0,  \
                             g_0_xxzzz_0_xxzzzzzz_1,  \
                             g_0_xxzzz_0_xyyyyyyz_0,  \
                             g_0_xxzzz_0_xyyyyyyz_1,  \
                             g_0_xxzzz_0_xyyyyyz_1,   \
                             g_0_xxzzz_0_xyyyyyzz_0,  \
                             g_0_xxzzz_0_xyyyyyzz_1,  \
                             g_0_xxzzz_0_xyyyyzz_1,   \
                             g_0_xxzzz_0_xyyyyzzz_0,  \
                             g_0_xxzzz_0_xyyyyzzz_1,  \
                             g_0_xxzzz_0_xyyyzzz_1,   \
                             g_0_xxzzz_0_xyyyzzzz_0,  \
                             g_0_xxzzz_0_xyyyzzzz_1,  \
                             g_0_xxzzz_0_xyyzzzz_1,   \
                             g_0_xxzzz_0_xyyzzzzz_0,  \
                             g_0_xxzzz_0_xyyzzzzz_1,  \
                             g_0_xxzzz_0_xyzzzzz_1,   \
                             g_0_xxzzz_0_xyzzzzzz_0,  \
                             g_0_xxzzz_0_xyzzzzzz_1,  \
                             g_0_xxzzz_0_xzzzzzz_1,   \
                             g_0_xxzzz_0_xzzzzzzz_0,  \
                             g_0_xxzzz_0_xzzzzzzz_1,  \
                             g_0_xxzzz_0_yyyyyyyy_0,  \
                             g_0_xxzzz_0_yyyyyyyy_1,  \
                             g_0_xxzzz_0_yyyyyyyz_0,  \
                             g_0_xxzzz_0_yyyyyyyz_1,  \
                             g_0_xxzzz_0_yyyyyyz_1,   \
                             g_0_xxzzz_0_yyyyyyzz_0,  \
                             g_0_xxzzz_0_yyyyyyzz_1,  \
                             g_0_xxzzz_0_yyyyyzz_1,   \
                             g_0_xxzzz_0_yyyyyzzz_0,  \
                             g_0_xxzzz_0_yyyyyzzz_1,  \
                             g_0_xxzzz_0_yyyyzzz_1,   \
                             g_0_xxzzz_0_yyyyzzzz_0,  \
                             g_0_xxzzz_0_yyyyzzzz_1,  \
                             g_0_xxzzz_0_yyyzzzz_1,   \
                             g_0_xxzzz_0_yyyzzzzz_0,  \
                             g_0_xxzzz_0_yyyzzzzz_1,  \
                             g_0_xxzzz_0_yyzzzzz_1,   \
                             g_0_xxzzz_0_yyzzzzzz_0,  \
                             g_0_xxzzz_0_yyzzzzzz_1,  \
                             g_0_xxzzz_0_yzzzzzz_1,   \
                             g_0_xxzzz_0_yzzzzzzz_0,  \
                             g_0_xxzzz_0_yzzzzzzz_1,  \
                             g_0_xxzzz_0_zzzzzzz_1,   \
                             g_0_xxzzz_0_zzzzzzzz_0,  \
                             g_0_xxzzz_0_zzzzzzzz_1,  \
                             g_0_xzzz_0_xxxxxxxz_0,   \
                             g_0_xzzz_0_xxxxxxxz_1,   \
                             g_0_xzzz_0_xxxxxxyz_0,   \
                             g_0_xzzz_0_xxxxxxyz_1,   \
                             g_0_xzzz_0_xxxxxxzz_0,   \
                             g_0_xzzz_0_xxxxxxzz_1,   \
                             g_0_xzzz_0_xxxxxyyz_0,   \
                             g_0_xzzz_0_xxxxxyyz_1,   \
                             g_0_xzzz_0_xxxxxyzz_0,   \
                             g_0_xzzz_0_xxxxxyzz_1,   \
                             g_0_xzzz_0_xxxxxzzz_0,   \
                             g_0_xzzz_0_xxxxxzzz_1,   \
                             g_0_xzzz_0_xxxxyyyz_0,   \
                             g_0_xzzz_0_xxxxyyyz_1,   \
                             g_0_xzzz_0_xxxxyyzz_0,   \
                             g_0_xzzz_0_xxxxyyzz_1,   \
                             g_0_xzzz_0_xxxxyzzz_0,   \
                             g_0_xzzz_0_xxxxyzzz_1,   \
                             g_0_xzzz_0_xxxxzzzz_0,   \
                             g_0_xzzz_0_xxxxzzzz_1,   \
                             g_0_xzzz_0_xxxyyyyz_0,   \
                             g_0_xzzz_0_xxxyyyyz_1,   \
                             g_0_xzzz_0_xxxyyyzz_0,   \
                             g_0_xzzz_0_xxxyyyzz_1,   \
                             g_0_xzzz_0_xxxyyzzz_0,   \
                             g_0_xzzz_0_xxxyyzzz_1,   \
                             g_0_xzzz_0_xxxyzzzz_0,   \
                             g_0_xzzz_0_xxxyzzzz_1,   \
                             g_0_xzzz_0_xxxzzzzz_0,   \
                             g_0_xzzz_0_xxxzzzzz_1,   \
                             g_0_xzzz_0_xxyyyyyz_0,   \
                             g_0_xzzz_0_xxyyyyyz_1,   \
                             g_0_xzzz_0_xxyyyyzz_0,   \
                             g_0_xzzz_0_xxyyyyzz_1,   \
                             g_0_xzzz_0_xxyyyzzz_0,   \
                             g_0_xzzz_0_xxyyyzzz_1,   \
                             g_0_xzzz_0_xxyyzzzz_0,   \
                             g_0_xzzz_0_xxyyzzzz_1,   \
                             g_0_xzzz_0_xxyzzzzz_0,   \
                             g_0_xzzz_0_xxyzzzzz_1,   \
                             g_0_xzzz_0_xxzzzzzz_0,   \
                             g_0_xzzz_0_xxzzzzzz_1,   \
                             g_0_xzzz_0_xyyyyyyz_0,   \
                             g_0_xzzz_0_xyyyyyyz_1,   \
                             g_0_xzzz_0_xyyyyyzz_0,   \
                             g_0_xzzz_0_xyyyyyzz_1,   \
                             g_0_xzzz_0_xyyyyzzz_0,   \
                             g_0_xzzz_0_xyyyyzzz_1,   \
                             g_0_xzzz_0_xyyyzzzz_0,   \
                             g_0_xzzz_0_xyyyzzzz_1,   \
                             g_0_xzzz_0_xyyzzzzz_0,   \
                             g_0_xzzz_0_xyyzzzzz_1,   \
                             g_0_xzzz_0_xyzzzzzz_0,   \
                             g_0_xzzz_0_xyzzzzzz_1,   \
                             g_0_xzzz_0_xzzzzzzz_0,   \
                             g_0_xzzz_0_xzzzzzzz_1,   \
                             g_0_xzzz_0_yyyyyyyy_0,   \
                             g_0_xzzz_0_yyyyyyyy_1,   \
                             g_0_xzzz_0_yyyyyyyz_0,   \
                             g_0_xzzz_0_yyyyyyyz_1,   \
                             g_0_xzzz_0_yyyyyyzz_0,   \
                             g_0_xzzz_0_yyyyyyzz_1,   \
                             g_0_xzzz_0_yyyyyzzz_0,   \
                             g_0_xzzz_0_yyyyyzzz_1,   \
                             g_0_xzzz_0_yyyyzzzz_0,   \
                             g_0_xzzz_0_yyyyzzzz_1,   \
                             g_0_xzzz_0_yyyzzzzz_0,   \
                             g_0_xzzz_0_yyyzzzzz_1,   \
                             g_0_xzzz_0_yyzzzzzz_0,   \
                             g_0_xzzz_0_yyzzzzzz_1,   \
                             g_0_xzzz_0_yzzzzzzz_0,   \
                             g_0_xzzz_0_yzzzzzzz_1,   \
                             g_0_xzzz_0_zzzzzzzz_0,   \
                             g_0_xzzz_0_zzzzzzzz_1,   \
                             wp_x,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzz_0_xxxxxxxx_0[i] = 2.0 * g_0_xxxz_0_xxxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_xxxxxxxx_0[i] * pb_z + g_0_xxxzz_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxxxxy_0[i] = 2.0 * g_0_xxxz_0_xxxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxxxxy_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_xxxxxxxy_0[i] * pb_z + g_0_xxxzz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxxxxz_0[i] = 2.0 * g_0_xzzz_0_xxxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxxxxz_1[i] * fti_ab_0 +
                                     7.0 * g_0_xxzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxxxz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxxxyy_0[i] = 2.0 * g_0_xxxz_0_xxxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxxxyy_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_xxxxxxyy_0[i] * pb_z + g_0_xxxzz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxxxyz_0[i] = 2.0 * g_0_xzzz_0_xxxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxxxyz_1[i] * fti_ab_0 +
                                     6.0 * g_0_xxzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxxyz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxxxzz_0[i] = 2.0 * g_0_xzzz_0_xxxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxxxzz_1[i] * fti_ab_0 +
                                     6.0 * g_0_xxzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxxzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxxyyy_0[i] = 2.0 * g_0_xxxz_0_xxxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxxyyy_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_xxxxxyyy_0[i] * pb_z + g_0_xxxzz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxxyyz_0[i] = 2.0 * g_0_xzzz_0_xxxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxyyz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxxyzz_0[i] = 2.0 * g_0_xzzz_0_xxxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxxyzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxyzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxxzzz_0[i] = 2.0 * g_0_xzzz_0_xxxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxxzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxzzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxyyyy_0[i] = 2.0 * g_0_xxxz_0_xxxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_xxxxyyyy_0[i] * pb_z + g_0_xxxzz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxyyyz_0[i] = 2.0 * g_0_xzzz_0_xxxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyyyz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxyyzz_0[i] = 2.0 * g_0_xzzz_0_xxxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyyzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxyzzz_0[i] = 2.0 * g_0_xzzz_0_xxxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxyzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyzzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxzzzz_0[i] = 2.0 * g_0_xzzz_0_xxxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxzzzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxyyyyy_0[i] = 2.0 * g_0_xxxz_0_xxxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_xxxyyyyy_0[i] * pb_z + g_0_xxxzz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxyyyyz_0[i] = 2.0 * g_0_xzzz_0_xxxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyyyz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxyyyzz_0[i] = 2.0 * g_0_xzzz_0_xxxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyyzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxyyzzz_0[i] = 2.0 * g_0_xzzz_0_xxxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyzzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxyzzzz_0[i] = 2.0 * g_0_xzzz_0_xxxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxyzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyzzzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxzzzzz_0[i] = 2.0 * g_0_xzzz_0_xxxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxzzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxzzzzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyyyyyy_0[i] = 2.0 * g_0_xxxz_0_xxyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxyyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_xxyyyyyy_0[i] * pb_z + g_0_xxxzz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxyyyyyz_0[i] = 2.0 * g_0_xzzz_0_xxyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyyyz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyyyyzz_0[i] = 2.0 * g_0_xzzz_0_xxyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyyzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyyyzzz_0[i] = 2.0 * g_0_xzzz_0_xxyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyzzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyyzzzz_0[i] = 2.0 * g_0_xzzz_0_xxyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyzzzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyzzzzz_0[i] = 2.0 * g_0_xzzz_0_xxyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyzzzzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxzzzzzz_0[i] = 2.0 * g_0_xzzz_0_xxzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxzzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxzzzzzz_0[i] * pb_x +
                                     g_0_xxzzz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyyyyyy_0[i] = 2.0 * g_0_xxxz_0_xyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxzz_0_xyyyyyyy_0[i] * pb_z + g_0_xxxzz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xyyyyyyz_0[i] = 2.0 * g_0_xzzz_0_xyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyyyz_0[i] * pb_x + g_0_xxzzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyyyyzz_0[i] = 2.0 * g_0_xzzz_0_xyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyyzz_0[i] * pb_x + g_0_xxzzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyyyzzz_0[i] = 2.0 * g_0_xzzz_0_xyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyzzz_0[i] * pb_x + g_0_xxzzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyyzzzz_0[i] = 2.0 * g_0_xzzz_0_xyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyzzzz_0[i] * pb_x + g_0_xxzzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyzzzzz_0[i] = 2.0 * g_0_xzzz_0_xyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyzzzzz_0[i] * pb_x + g_0_xxzzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyzzzzzz_0[i] = 2.0 * g_0_xzzz_0_xyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyzzzzzz_0[i] * pb_x + g_0_xxzzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xzzzzzzz_0[i] = 2.0 * g_0_xzzz_0_xzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xzzzzzzz_0[i] * pb_x + g_0_xxzzz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyyyyy_0[i] = 2.0 * g_0_xzzz_0_yyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yyyyyyyy_0[i] * pb_x + g_0_xxzzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyyyyz_0[i] = 2.0 * g_0_xzzz_0_yyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yyyyyyyz_0[i] * pb_x + g_0_xxzzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyyyzz_0[i] = 2.0 * g_0_xzzz_0_yyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yyyyyyzz_0[i] * pb_x + g_0_xxzzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyyzzz_0[i] = 2.0 * g_0_xzzz_0_yyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yyyyyzzz_0[i] * pb_x + g_0_xxzzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyzzzz_0[i] = 2.0 * g_0_xzzz_0_yyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yyyyzzzz_0[i] * pb_x + g_0_xxzzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyzzzzz_0[i] = 2.0 * g_0_xzzz_0_yyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yyyzzzzz_0[i] * pb_x + g_0_xxzzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyzzzzzz_0[i] = 2.0 * g_0_xzzz_0_yyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yyzzzzzz_0[i] * pb_x + g_0_xxzzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yzzzzzzz_0[i] = 2.0 * g_0_xzzz_0_yzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_yzzzzzzz_0[i] * pb_x + g_0_xxzzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_zzzzzzzz_0[i] = 2.0 * g_0_xzzz_0_zzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_zzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_zzzzzzzz_0[i] * pb_x + g_0_xxzzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 450-495 components of targeted buffer : SISL

    auto g_0_xxyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 450);

    auto g_0_xxyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 451);

    auto g_0_xxyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 452);

    auto g_0_xxyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 453);

    auto g_0_xxyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 454);

    auto g_0_xxyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 455);

    auto g_0_xxyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 456);

    auto g_0_xxyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 457);

    auto g_0_xxyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 458);

    auto g_0_xxyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 459);

    auto g_0_xxyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 460);

    auto g_0_xxyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 461);

    auto g_0_xxyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 462);

    auto g_0_xxyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 463);

    auto g_0_xxyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 464);

    auto g_0_xxyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 465);

    auto g_0_xxyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 466);

    auto g_0_xxyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 467);

    auto g_0_xxyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 468);

    auto g_0_xxyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 469);

    auto g_0_xxyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 470);

    auto g_0_xxyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 471);

    auto g_0_xxyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 472);

    auto g_0_xxyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 473);

    auto g_0_xxyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 474);

    auto g_0_xxyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 475);

    auto g_0_xxyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 476);

    auto g_0_xxyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 477);

    auto g_0_xxyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 478);

    auto g_0_xxyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 479);

    auto g_0_xxyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 480);

    auto g_0_xxyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 481);

    auto g_0_xxyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 482);

    auto g_0_xxyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 483);

    auto g_0_xxyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 484);

    auto g_0_xxyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 485);

    auto g_0_xxyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 486);

    auto g_0_xxyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 487);

    auto g_0_xxyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 488);

    auto g_0_xxyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 489);

    auto g_0_xxyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 490);

    auto g_0_xxyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 491);

    auto g_0_xxyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 492);

    auto g_0_xxyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 493);

    auto g_0_xxyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 494);

#pragma omp simd aligned(g_0_xxyy_0_xxxxxxxx_0,       \
                             g_0_xxyy_0_xxxxxxxx_1,   \
                             g_0_xxyy_0_xxxxxxxz_0,   \
                             g_0_xxyy_0_xxxxxxxz_1,   \
                             g_0_xxyy_0_xxxxxxzz_0,   \
                             g_0_xxyy_0_xxxxxxzz_1,   \
                             g_0_xxyy_0_xxxxxzzz_0,   \
                             g_0_xxyy_0_xxxxxzzz_1,   \
                             g_0_xxyy_0_xxxxzzzz_0,   \
                             g_0_xxyy_0_xxxxzzzz_1,   \
                             g_0_xxyy_0_xxxzzzzz_0,   \
                             g_0_xxyy_0_xxxzzzzz_1,   \
                             g_0_xxyy_0_xxzzzzzz_0,   \
                             g_0_xxyy_0_xxzzzzzz_1,   \
                             g_0_xxyy_0_xzzzzzzz_0,   \
                             g_0_xxyy_0_xzzzzzzz_1,   \
                             g_0_xxyyy_0_xxxxxxxx_0,  \
                             g_0_xxyyy_0_xxxxxxxx_1,  \
                             g_0_xxyyy_0_xxxxxxxz_0,  \
                             g_0_xxyyy_0_xxxxxxxz_1,  \
                             g_0_xxyyy_0_xxxxxxzz_0,  \
                             g_0_xxyyy_0_xxxxxxzz_1,  \
                             g_0_xxyyy_0_xxxxxzzz_0,  \
                             g_0_xxyyy_0_xxxxxzzz_1,  \
                             g_0_xxyyy_0_xxxxzzzz_0,  \
                             g_0_xxyyy_0_xxxxzzzz_1,  \
                             g_0_xxyyy_0_xxxzzzzz_0,  \
                             g_0_xxyyy_0_xxxzzzzz_1,  \
                             g_0_xxyyy_0_xxzzzzzz_0,  \
                             g_0_xxyyy_0_xxzzzzzz_1,  \
                             g_0_xxyyy_0_xzzzzzzz_0,  \
                             g_0_xxyyy_0_xzzzzzzz_1,  \
                             g_0_xxyyyy_0_xxxxxxxx_0, \
                             g_0_xxyyyy_0_xxxxxxxy_0, \
                             g_0_xxyyyy_0_xxxxxxxz_0, \
                             g_0_xxyyyy_0_xxxxxxyy_0, \
                             g_0_xxyyyy_0_xxxxxxyz_0, \
                             g_0_xxyyyy_0_xxxxxxzz_0, \
                             g_0_xxyyyy_0_xxxxxyyy_0, \
                             g_0_xxyyyy_0_xxxxxyyz_0, \
                             g_0_xxyyyy_0_xxxxxyzz_0, \
                             g_0_xxyyyy_0_xxxxxzzz_0, \
                             g_0_xxyyyy_0_xxxxyyyy_0, \
                             g_0_xxyyyy_0_xxxxyyyz_0, \
                             g_0_xxyyyy_0_xxxxyyzz_0, \
                             g_0_xxyyyy_0_xxxxyzzz_0, \
                             g_0_xxyyyy_0_xxxxzzzz_0, \
                             g_0_xxyyyy_0_xxxyyyyy_0, \
                             g_0_xxyyyy_0_xxxyyyyz_0, \
                             g_0_xxyyyy_0_xxxyyyzz_0, \
                             g_0_xxyyyy_0_xxxyyzzz_0, \
                             g_0_xxyyyy_0_xxxyzzzz_0, \
                             g_0_xxyyyy_0_xxxzzzzz_0, \
                             g_0_xxyyyy_0_xxyyyyyy_0, \
                             g_0_xxyyyy_0_xxyyyyyz_0, \
                             g_0_xxyyyy_0_xxyyyyzz_0, \
                             g_0_xxyyyy_0_xxyyyzzz_0, \
                             g_0_xxyyyy_0_xxyyzzzz_0, \
                             g_0_xxyyyy_0_xxyzzzzz_0, \
                             g_0_xxyyyy_0_xxzzzzzz_0, \
                             g_0_xxyyyy_0_xyyyyyyy_0, \
                             g_0_xxyyyy_0_xyyyyyyz_0, \
                             g_0_xxyyyy_0_xyyyyyzz_0, \
                             g_0_xxyyyy_0_xyyyyzzz_0, \
                             g_0_xxyyyy_0_xyyyzzzz_0, \
                             g_0_xxyyyy_0_xyyzzzzz_0, \
                             g_0_xxyyyy_0_xyzzzzzz_0, \
                             g_0_xxyyyy_0_xzzzzzzz_0, \
                             g_0_xxyyyy_0_yyyyyyyy_0, \
                             g_0_xxyyyy_0_yyyyyyyz_0, \
                             g_0_xxyyyy_0_yyyyyyzz_0, \
                             g_0_xxyyyy_0_yyyyyzzz_0, \
                             g_0_xxyyyy_0_yyyyzzzz_0, \
                             g_0_xxyyyy_0_yyyzzzzz_0, \
                             g_0_xxyyyy_0_yyzzzzzz_0, \
                             g_0_xxyyyy_0_yzzzzzzz_0, \
                             g_0_xxyyyy_0_zzzzzzzz_0, \
                             g_0_xyyyy_0_xxxxxxxy_0,  \
                             g_0_xyyyy_0_xxxxxxxy_1,  \
                             g_0_xyyyy_0_xxxxxxy_1,   \
                             g_0_xyyyy_0_xxxxxxyy_0,  \
                             g_0_xyyyy_0_xxxxxxyy_1,  \
                             g_0_xyyyy_0_xxxxxxyz_0,  \
                             g_0_xyyyy_0_xxxxxxyz_1,  \
                             g_0_xyyyy_0_xxxxxyy_1,   \
                             g_0_xyyyy_0_xxxxxyyy_0,  \
                             g_0_xyyyy_0_xxxxxyyy_1,  \
                             g_0_xyyyy_0_xxxxxyyz_0,  \
                             g_0_xyyyy_0_xxxxxyyz_1,  \
                             g_0_xyyyy_0_xxxxxyz_1,   \
                             g_0_xyyyy_0_xxxxxyzz_0,  \
                             g_0_xyyyy_0_xxxxxyzz_1,  \
                             g_0_xyyyy_0_xxxxyyy_1,   \
                             g_0_xyyyy_0_xxxxyyyy_0,  \
                             g_0_xyyyy_0_xxxxyyyy_1,  \
                             g_0_xyyyy_0_xxxxyyyz_0,  \
                             g_0_xyyyy_0_xxxxyyyz_1,  \
                             g_0_xyyyy_0_xxxxyyz_1,   \
                             g_0_xyyyy_0_xxxxyyzz_0,  \
                             g_0_xyyyy_0_xxxxyyzz_1,  \
                             g_0_xyyyy_0_xxxxyzz_1,   \
                             g_0_xyyyy_0_xxxxyzzz_0,  \
                             g_0_xyyyy_0_xxxxyzzz_1,  \
                             g_0_xyyyy_0_xxxyyyy_1,   \
                             g_0_xyyyy_0_xxxyyyyy_0,  \
                             g_0_xyyyy_0_xxxyyyyy_1,  \
                             g_0_xyyyy_0_xxxyyyyz_0,  \
                             g_0_xyyyy_0_xxxyyyyz_1,  \
                             g_0_xyyyy_0_xxxyyyz_1,   \
                             g_0_xyyyy_0_xxxyyyzz_0,  \
                             g_0_xyyyy_0_xxxyyyzz_1,  \
                             g_0_xyyyy_0_xxxyyzz_1,   \
                             g_0_xyyyy_0_xxxyyzzz_0,  \
                             g_0_xyyyy_0_xxxyyzzz_1,  \
                             g_0_xyyyy_0_xxxyzzz_1,   \
                             g_0_xyyyy_0_xxxyzzzz_0,  \
                             g_0_xyyyy_0_xxxyzzzz_1,  \
                             g_0_xyyyy_0_xxyyyyy_1,   \
                             g_0_xyyyy_0_xxyyyyyy_0,  \
                             g_0_xyyyy_0_xxyyyyyy_1,  \
                             g_0_xyyyy_0_xxyyyyyz_0,  \
                             g_0_xyyyy_0_xxyyyyyz_1,  \
                             g_0_xyyyy_0_xxyyyyz_1,   \
                             g_0_xyyyy_0_xxyyyyzz_0,  \
                             g_0_xyyyy_0_xxyyyyzz_1,  \
                             g_0_xyyyy_0_xxyyyzz_1,   \
                             g_0_xyyyy_0_xxyyyzzz_0,  \
                             g_0_xyyyy_0_xxyyyzzz_1,  \
                             g_0_xyyyy_0_xxyyzzz_1,   \
                             g_0_xyyyy_0_xxyyzzzz_0,  \
                             g_0_xyyyy_0_xxyyzzzz_1,  \
                             g_0_xyyyy_0_xxyzzzz_1,   \
                             g_0_xyyyy_0_xxyzzzzz_0,  \
                             g_0_xyyyy_0_xxyzzzzz_1,  \
                             g_0_xyyyy_0_xyyyyyy_1,   \
                             g_0_xyyyy_0_xyyyyyyy_0,  \
                             g_0_xyyyy_0_xyyyyyyy_1,  \
                             g_0_xyyyy_0_xyyyyyyz_0,  \
                             g_0_xyyyy_0_xyyyyyyz_1,  \
                             g_0_xyyyy_0_xyyyyyz_1,   \
                             g_0_xyyyy_0_xyyyyyzz_0,  \
                             g_0_xyyyy_0_xyyyyyzz_1,  \
                             g_0_xyyyy_0_xyyyyzz_1,   \
                             g_0_xyyyy_0_xyyyyzzz_0,  \
                             g_0_xyyyy_0_xyyyyzzz_1,  \
                             g_0_xyyyy_0_xyyyzzz_1,   \
                             g_0_xyyyy_0_xyyyzzzz_0,  \
                             g_0_xyyyy_0_xyyyzzzz_1,  \
                             g_0_xyyyy_0_xyyzzzz_1,   \
                             g_0_xyyyy_0_xyyzzzzz_0,  \
                             g_0_xyyyy_0_xyyzzzzz_1,  \
                             g_0_xyyyy_0_xyzzzzz_1,   \
                             g_0_xyyyy_0_xyzzzzzz_0,  \
                             g_0_xyyyy_0_xyzzzzzz_1,  \
                             g_0_xyyyy_0_yyyyyyy_1,   \
                             g_0_xyyyy_0_yyyyyyyy_0,  \
                             g_0_xyyyy_0_yyyyyyyy_1,  \
                             g_0_xyyyy_0_yyyyyyyz_0,  \
                             g_0_xyyyy_0_yyyyyyyz_1,  \
                             g_0_xyyyy_0_yyyyyyz_1,   \
                             g_0_xyyyy_0_yyyyyyzz_0,  \
                             g_0_xyyyy_0_yyyyyyzz_1,  \
                             g_0_xyyyy_0_yyyyyzz_1,   \
                             g_0_xyyyy_0_yyyyyzzz_0,  \
                             g_0_xyyyy_0_yyyyyzzz_1,  \
                             g_0_xyyyy_0_yyyyzzz_1,   \
                             g_0_xyyyy_0_yyyyzzzz_0,  \
                             g_0_xyyyy_0_yyyyzzzz_1,  \
                             g_0_xyyyy_0_yyyzzzz_1,   \
                             g_0_xyyyy_0_yyyzzzzz_0,  \
                             g_0_xyyyy_0_yyyzzzzz_1,  \
                             g_0_xyyyy_0_yyzzzzz_1,   \
                             g_0_xyyyy_0_yyzzzzzz_0,  \
                             g_0_xyyyy_0_yyzzzzzz_1,  \
                             g_0_xyyyy_0_yzzzzzz_1,   \
                             g_0_xyyyy_0_yzzzzzzz_0,  \
                             g_0_xyyyy_0_yzzzzzzz_1,  \
                             g_0_xyyyy_0_zzzzzzzz_0,  \
                             g_0_xyyyy_0_zzzzzzzz_1,  \
                             g_0_yyyy_0_xxxxxxxy_0,   \
                             g_0_yyyy_0_xxxxxxxy_1,   \
                             g_0_yyyy_0_xxxxxxyy_0,   \
                             g_0_yyyy_0_xxxxxxyy_1,   \
                             g_0_yyyy_0_xxxxxxyz_0,   \
                             g_0_yyyy_0_xxxxxxyz_1,   \
                             g_0_yyyy_0_xxxxxyyy_0,   \
                             g_0_yyyy_0_xxxxxyyy_1,   \
                             g_0_yyyy_0_xxxxxyyz_0,   \
                             g_0_yyyy_0_xxxxxyyz_1,   \
                             g_0_yyyy_0_xxxxxyzz_0,   \
                             g_0_yyyy_0_xxxxxyzz_1,   \
                             g_0_yyyy_0_xxxxyyyy_0,   \
                             g_0_yyyy_0_xxxxyyyy_1,   \
                             g_0_yyyy_0_xxxxyyyz_0,   \
                             g_0_yyyy_0_xxxxyyyz_1,   \
                             g_0_yyyy_0_xxxxyyzz_0,   \
                             g_0_yyyy_0_xxxxyyzz_1,   \
                             g_0_yyyy_0_xxxxyzzz_0,   \
                             g_0_yyyy_0_xxxxyzzz_1,   \
                             g_0_yyyy_0_xxxyyyyy_0,   \
                             g_0_yyyy_0_xxxyyyyy_1,   \
                             g_0_yyyy_0_xxxyyyyz_0,   \
                             g_0_yyyy_0_xxxyyyyz_1,   \
                             g_0_yyyy_0_xxxyyyzz_0,   \
                             g_0_yyyy_0_xxxyyyzz_1,   \
                             g_0_yyyy_0_xxxyyzzz_0,   \
                             g_0_yyyy_0_xxxyyzzz_1,   \
                             g_0_yyyy_0_xxxyzzzz_0,   \
                             g_0_yyyy_0_xxxyzzzz_1,   \
                             g_0_yyyy_0_xxyyyyyy_0,   \
                             g_0_yyyy_0_xxyyyyyy_1,   \
                             g_0_yyyy_0_xxyyyyyz_0,   \
                             g_0_yyyy_0_xxyyyyyz_1,   \
                             g_0_yyyy_0_xxyyyyzz_0,   \
                             g_0_yyyy_0_xxyyyyzz_1,   \
                             g_0_yyyy_0_xxyyyzzz_0,   \
                             g_0_yyyy_0_xxyyyzzz_1,   \
                             g_0_yyyy_0_xxyyzzzz_0,   \
                             g_0_yyyy_0_xxyyzzzz_1,   \
                             g_0_yyyy_0_xxyzzzzz_0,   \
                             g_0_yyyy_0_xxyzzzzz_1,   \
                             g_0_yyyy_0_xyyyyyyy_0,   \
                             g_0_yyyy_0_xyyyyyyy_1,   \
                             g_0_yyyy_0_xyyyyyyz_0,   \
                             g_0_yyyy_0_xyyyyyyz_1,   \
                             g_0_yyyy_0_xyyyyyzz_0,   \
                             g_0_yyyy_0_xyyyyyzz_1,   \
                             g_0_yyyy_0_xyyyyzzz_0,   \
                             g_0_yyyy_0_xyyyyzzz_1,   \
                             g_0_yyyy_0_xyyyzzzz_0,   \
                             g_0_yyyy_0_xyyyzzzz_1,   \
                             g_0_yyyy_0_xyyzzzzz_0,   \
                             g_0_yyyy_0_xyyzzzzz_1,   \
                             g_0_yyyy_0_xyzzzzzz_0,   \
                             g_0_yyyy_0_xyzzzzzz_1,   \
                             g_0_yyyy_0_yyyyyyyy_0,   \
                             g_0_yyyy_0_yyyyyyyy_1,   \
                             g_0_yyyy_0_yyyyyyyz_0,   \
                             g_0_yyyy_0_yyyyyyyz_1,   \
                             g_0_yyyy_0_yyyyyyzz_0,   \
                             g_0_yyyy_0_yyyyyyzz_1,   \
                             g_0_yyyy_0_yyyyyzzz_0,   \
                             g_0_yyyy_0_yyyyyzzz_1,   \
                             g_0_yyyy_0_yyyyzzzz_0,   \
                             g_0_yyyy_0_yyyyzzzz_1,   \
                             g_0_yyyy_0_yyyzzzzz_0,   \
                             g_0_yyyy_0_yyyzzzzz_1,   \
                             g_0_yyyy_0_yyzzzzzz_0,   \
                             g_0_yyyy_0_yyzzzzzz_1,   \
                             g_0_yyyy_0_yzzzzzzz_0,   \
                             g_0_yyyy_0_yzzzzzzz_1,   \
                             g_0_yyyy_0_zzzzzzzz_0,   \
                             g_0_yyyy_0_zzzzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyy_0_xxxxxxxx_0[i] = 3.0 * g_0_xxyy_0_xxxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_xxxxxxxx_0[i] * pb_y + g_0_xxyyy_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxxxxxy_0[i] = g_0_yyyy_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxxxy_1[i] * fti_ab_0 +
                                     7.0 * g_0_xyyyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxxxxy_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxxxxz_0[i] = 3.0 * g_0_xxyy_0_xxxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxxxz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_xxxxxxxz_0[i] * pb_y + g_0_xxyyy_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxxxxyy_0[i] = g_0_yyyy_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxxyy_1[i] * fti_ab_0 +
                                     6.0 * g_0_xyyyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxxxyy_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxxxyz_0[i] = g_0_yyyy_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxxyz_1[i] * fti_ab_0 +
                                     6.0 * g_0_xyyyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxxxyz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxxxzz_0[i] = 3.0 * g_0_xxyy_0_xxxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxxzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_xxxxxxzz_0[i] * pb_y + g_0_xxyyy_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxxxyyy_0[i] = g_0_yyyy_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxyyy_1[i] * fti_ab_0 +
                                     5.0 * g_0_xyyyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxxyyy_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxxyyz_0[i] = g_0_yyyy_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xyyyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxxyyz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxxyzz_0[i] = g_0_yyyy_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxyzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xyyyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxxyzz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxxzzz_0[i] = 3.0 * g_0_xxyy_0_xxxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_xxxxxzzz_0[i] * pb_y + g_0_xxyyy_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxxyyyy_0[i] = g_0_yyyy_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxyyyy_1[i] * fti_ab_0 +
                                     4.0 * g_0_xyyyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxyyyy_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxyyyz_0[i] = g_0_yyyy_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xyyyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxyyyz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxyyzz_0[i] = g_0_yyyy_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xyyyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxyyzz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxyzzz_0[i] = g_0_yyyy_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxyzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xyyyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxyzzz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxzzzz_0[i] = 3.0 * g_0_xxyy_0_xxxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_xxxxzzzz_0[i] * pb_y + g_0_xxyyy_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxyyyyy_0[i] = g_0_yyyy_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyyyyy_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxyyyyy_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxyyyyz_0[i] = g_0_yyyy_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxyyyyz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxyyyzz_0[i] = g_0_yyyy_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxyyyzz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxyyzzz_0[i] = g_0_yyyy_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxyyzzz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxyzzzz_0[i] = g_0_yyyy_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxyzzzz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxzzzzz_0[i] = 3.0 * g_0_xxyy_0_xxxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_xxxzzzzz_0[i] * pb_y + g_0_xxyyy_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxyyyyyy_0[i] = g_0_yyyy_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyyyyy_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxyyyyyy_0[i] * pb_x +
                                     g_0_xyyyy_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyyyyyz_0[i] = g_0_yyyy_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxyyyyyz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyyyyzz_0[i] = g_0_yyyy_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxyyyyzz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyyyzzz_0[i] = g_0_yyyy_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxyyyzzz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyyzzzz_0[i] = g_0_yyyy_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxyyzzzz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyzzzzz_0[i] = g_0_yyyy_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxyzzzzz_0[i] * pb_x +
                                     g_0_xyyyy_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxzzzzzz_0[i] = 3.0 * g_0_xxyy_0_xxzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_xxzzzzzz_0[i] * pb_y + g_0_xxyyy_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xyyyyyyy_0[i] = g_0_yyyy_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyyy_1[i] * fi_abcd_0 +
                                     g_0_xyyyy_0_xyyyyyyy_0[i] * pb_x + g_0_xyyyy_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyyyyyz_0[i] = g_0_yyyy_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyyz_1[i] * fi_abcd_0 +
                                     g_0_xyyyy_0_xyyyyyyz_0[i] * pb_x + g_0_xyyyy_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyyyyzz_0[i] = g_0_yyyy_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyyyzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyy_0_xyyyyyzz_0[i] * pb_x + g_0_xyyyy_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyyyzzz_0[i] = g_0_yyyy_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyyzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyy_0_xyyyyzzz_0[i] * pb_x + g_0_xyyyy_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyyzzzz_0[i] = g_0_yyyy_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyzzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyy_0_xyyyzzzz_0[i] * pb_x + g_0_xyyyy_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyzzzzz_0[i] = g_0_yyyy_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyzzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyzzzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyy_0_xyyzzzzz_0[i] * pb_x + g_0_xyyyy_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyzzzzzz_0[i] = g_0_yyyy_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yzzzzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyy_0_xyzzzzzz_0[i] * pb_x + g_0_xyyyy_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xzzzzzzz_0[i] = 3.0 * g_0_xxyy_0_xzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyy_0_xzzzzzzz_0[i] * pb_y + g_0_xxyyy_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_yyyyyyyy_0[i] = g_0_yyyy_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyyyy_0[i] * pb_x +
                                     g_0_xyyyy_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyyyyyz_0[i] = g_0_yyyy_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyyyz_0[i] * pb_x +
                                     g_0_xyyyy_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyyyyzz_0[i] = g_0_yyyy_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyyzz_0[i] * pb_x +
                                     g_0_xyyyy_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyyyzzz_0[i] = g_0_yyyy_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyzzz_0[i] * pb_x +
                                     g_0_xyyyy_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyyzzzz_0[i] = g_0_yyyy_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyzzzz_0[i] * pb_x +
                                     g_0_xyyyy_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyzzzzz_0[i] = g_0_yyyy_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyzzzzz_0[i] * pb_x +
                                     g_0_xyyyy_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyzzzzzz_0[i] = g_0_yyyy_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyzzzzzz_0[i] * pb_x +
                                     g_0_xyyyy_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yzzzzzzz_0[i] = g_0_yyyy_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yzzzzzzz_0[i] * pb_x +
                                     g_0_xyyyy_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_zzzzzzzz_0[i] = g_0_yyyy_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_zzzzzzzz_0[i] * pb_x +
                                     g_0_xyyyy_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 495-540 components of targeted buffer : SISL

    auto g_0_xxyyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 495);

    auto g_0_xxyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 496);

    auto g_0_xxyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 497);

    auto g_0_xxyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 498);

    auto g_0_xxyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 499);

    auto g_0_xxyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 500);

    auto g_0_xxyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 501);

    auto g_0_xxyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 502);

    auto g_0_xxyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 503);

    auto g_0_xxyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 504);

    auto g_0_xxyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 505);

    auto g_0_xxyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 506);

    auto g_0_xxyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 507);

    auto g_0_xxyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 508);

    auto g_0_xxyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 509);

    auto g_0_xxyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 510);

    auto g_0_xxyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 511);

    auto g_0_xxyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 512);

    auto g_0_xxyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 513);

    auto g_0_xxyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 514);

    auto g_0_xxyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 515);

    auto g_0_xxyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 516);

    auto g_0_xxyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 517);

    auto g_0_xxyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 518);

    auto g_0_xxyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 519);

    auto g_0_xxyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 520);

    auto g_0_xxyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 521);

    auto g_0_xxyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 522);

    auto g_0_xxyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 523);

    auto g_0_xxyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 524);

    auto g_0_xxyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 525);

    auto g_0_xxyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 526);

    auto g_0_xxyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 527);

    auto g_0_xxyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 528);

    auto g_0_xxyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 529);

    auto g_0_xxyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 530);

    auto g_0_xxyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 531);

    auto g_0_xxyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 532);

    auto g_0_xxyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 533);

    auto g_0_xxyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 534);

    auto g_0_xxyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 535);

    auto g_0_xxyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 536);

    auto g_0_xxyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 537);

    auto g_0_xxyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 538);

    auto g_0_xxyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 539);

#pragma omp simd aligned(g_0_xxyyy_0_xxxxxxx_1,       \
                             g_0_xxyyy_0_xxxxxxxx_0,  \
                             g_0_xxyyy_0_xxxxxxxx_1,  \
                             g_0_xxyyy_0_xxxxxxxy_0,  \
                             g_0_xxyyy_0_xxxxxxxy_1,  \
                             g_0_xxyyy_0_xxxxxxxz_0,  \
                             g_0_xxyyy_0_xxxxxxxz_1,  \
                             g_0_xxyyy_0_xxxxxxy_1,   \
                             g_0_xxyyy_0_xxxxxxyy_0,  \
                             g_0_xxyyy_0_xxxxxxyy_1,  \
                             g_0_xxyyy_0_xxxxxxyz_0,  \
                             g_0_xxyyy_0_xxxxxxyz_1,  \
                             g_0_xxyyy_0_xxxxxxz_1,   \
                             g_0_xxyyy_0_xxxxxxzz_0,  \
                             g_0_xxyyy_0_xxxxxxzz_1,  \
                             g_0_xxyyy_0_xxxxxyy_1,   \
                             g_0_xxyyy_0_xxxxxyyy_0,  \
                             g_0_xxyyy_0_xxxxxyyy_1,  \
                             g_0_xxyyy_0_xxxxxyyz_0,  \
                             g_0_xxyyy_0_xxxxxyyz_1,  \
                             g_0_xxyyy_0_xxxxxyz_1,   \
                             g_0_xxyyy_0_xxxxxyzz_0,  \
                             g_0_xxyyy_0_xxxxxyzz_1,  \
                             g_0_xxyyy_0_xxxxxzz_1,   \
                             g_0_xxyyy_0_xxxxxzzz_0,  \
                             g_0_xxyyy_0_xxxxxzzz_1,  \
                             g_0_xxyyy_0_xxxxyyy_1,   \
                             g_0_xxyyy_0_xxxxyyyy_0,  \
                             g_0_xxyyy_0_xxxxyyyy_1,  \
                             g_0_xxyyy_0_xxxxyyyz_0,  \
                             g_0_xxyyy_0_xxxxyyyz_1,  \
                             g_0_xxyyy_0_xxxxyyz_1,   \
                             g_0_xxyyy_0_xxxxyyzz_0,  \
                             g_0_xxyyy_0_xxxxyyzz_1,  \
                             g_0_xxyyy_0_xxxxyzz_1,   \
                             g_0_xxyyy_0_xxxxyzzz_0,  \
                             g_0_xxyyy_0_xxxxyzzz_1,  \
                             g_0_xxyyy_0_xxxxzzz_1,   \
                             g_0_xxyyy_0_xxxxzzzz_0,  \
                             g_0_xxyyy_0_xxxxzzzz_1,  \
                             g_0_xxyyy_0_xxxyyyy_1,   \
                             g_0_xxyyy_0_xxxyyyyy_0,  \
                             g_0_xxyyy_0_xxxyyyyy_1,  \
                             g_0_xxyyy_0_xxxyyyyz_0,  \
                             g_0_xxyyy_0_xxxyyyyz_1,  \
                             g_0_xxyyy_0_xxxyyyz_1,   \
                             g_0_xxyyy_0_xxxyyyzz_0,  \
                             g_0_xxyyy_0_xxxyyyzz_1,  \
                             g_0_xxyyy_0_xxxyyzz_1,   \
                             g_0_xxyyy_0_xxxyyzzz_0,  \
                             g_0_xxyyy_0_xxxyyzzz_1,  \
                             g_0_xxyyy_0_xxxyzzz_1,   \
                             g_0_xxyyy_0_xxxyzzzz_0,  \
                             g_0_xxyyy_0_xxxyzzzz_1,  \
                             g_0_xxyyy_0_xxxzzzz_1,   \
                             g_0_xxyyy_0_xxxzzzzz_0,  \
                             g_0_xxyyy_0_xxxzzzzz_1,  \
                             g_0_xxyyy_0_xxyyyyy_1,   \
                             g_0_xxyyy_0_xxyyyyyy_0,  \
                             g_0_xxyyy_0_xxyyyyyy_1,  \
                             g_0_xxyyy_0_xxyyyyyz_0,  \
                             g_0_xxyyy_0_xxyyyyyz_1,  \
                             g_0_xxyyy_0_xxyyyyz_1,   \
                             g_0_xxyyy_0_xxyyyyzz_0,  \
                             g_0_xxyyy_0_xxyyyyzz_1,  \
                             g_0_xxyyy_0_xxyyyzz_1,   \
                             g_0_xxyyy_0_xxyyyzzz_0,  \
                             g_0_xxyyy_0_xxyyyzzz_1,  \
                             g_0_xxyyy_0_xxyyzzz_1,   \
                             g_0_xxyyy_0_xxyyzzzz_0,  \
                             g_0_xxyyy_0_xxyyzzzz_1,  \
                             g_0_xxyyy_0_xxyzzzz_1,   \
                             g_0_xxyyy_0_xxyzzzzz_0,  \
                             g_0_xxyyy_0_xxyzzzzz_1,  \
                             g_0_xxyyy_0_xxzzzzz_1,   \
                             g_0_xxyyy_0_xxzzzzzz_0,  \
                             g_0_xxyyy_0_xxzzzzzz_1,  \
                             g_0_xxyyy_0_xyyyyyy_1,   \
                             g_0_xxyyy_0_xyyyyyyy_0,  \
                             g_0_xxyyy_0_xyyyyyyy_1,  \
                             g_0_xxyyy_0_xyyyyyyz_0,  \
                             g_0_xxyyy_0_xyyyyyyz_1,  \
                             g_0_xxyyy_0_xyyyyyz_1,   \
                             g_0_xxyyy_0_xyyyyyzz_0,  \
                             g_0_xxyyy_0_xyyyyyzz_1,  \
                             g_0_xxyyy_0_xyyyyzz_1,   \
                             g_0_xxyyy_0_xyyyyzzz_0,  \
                             g_0_xxyyy_0_xyyyyzzz_1,  \
                             g_0_xxyyy_0_xyyyzzz_1,   \
                             g_0_xxyyy_0_xyyyzzzz_0,  \
                             g_0_xxyyy_0_xyyyzzzz_1,  \
                             g_0_xxyyy_0_xyyzzzz_1,   \
                             g_0_xxyyy_0_xyyzzzzz_0,  \
                             g_0_xxyyy_0_xyyzzzzz_1,  \
                             g_0_xxyyy_0_xyzzzzz_1,   \
                             g_0_xxyyy_0_xyzzzzzz_0,  \
                             g_0_xxyyy_0_xyzzzzzz_1,  \
                             g_0_xxyyy_0_xzzzzzz_1,   \
                             g_0_xxyyy_0_xzzzzzzz_0,  \
                             g_0_xxyyy_0_xzzzzzzz_1,  \
                             g_0_xxyyy_0_yyyyyyy_1,   \
                             g_0_xxyyy_0_yyyyyyyy_0,  \
                             g_0_xxyyy_0_yyyyyyyy_1,  \
                             g_0_xxyyy_0_yyyyyyyz_0,  \
                             g_0_xxyyy_0_yyyyyyyz_1,  \
                             g_0_xxyyy_0_yyyyyyz_1,   \
                             g_0_xxyyy_0_yyyyyyzz_0,  \
                             g_0_xxyyy_0_yyyyyyzz_1,  \
                             g_0_xxyyy_0_yyyyyzz_1,   \
                             g_0_xxyyy_0_yyyyyzzz_0,  \
                             g_0_xxyyy_0_yyyyyzzz_1,  \
                             g_0_xxyyy_0_yyyyzzz_1,   \
                             g_0_xxyyy_0_yyyyzzzz_0,  \
                             g_0_xxyyy_0_yyyyzzzz_1,  \
                             g_0_xxyyy_0_yyyzzzz_1,   \
                             g_0_xxyyy_0_yyyzzzzz_0,  \
                             g_0_xxyyy_0_yyyzzzzz_1,  \
                             g_0_xxyyy_0_yyzzzzz_1,   \
                             g_0_xxyyy_0_yyzzzzzz_0,  \
                             g_0_xxyyy_0_yyzzzzzz_1,  \
                             g_0_xxyyy_0_yzzzzzz_1,   \
                             g_0_xxyyy_0_yzzzzzzz_0,  \
                             g_0_xxyyy_0_yzzzzzzz_1,  \
                             g_0_xxyyy_0_zzzzzzz_1,   \
                             g_0_xxyyy_0_zzzzzzzz_0,  \
                             g_0_xxyyy_0_zzzzzzzz_1,  \
                             g_0_xxyyyz_0_xxxxxxxx_0, \
                             g_0_xxyyyz_0_xxxxxxxy_0, \
                             g_0_xxyyyz_0_xxxxxxxz_0, \
                             g_0_xxyyyz_0_xxxxxxyy_0, \
                             g_0_xxyyyz_0_xxxxxxyz_0, \
                             g_0_xxyyyz_0_xxxxxxzz_0, \
                             g_0_xxyyyz_0_xxxxxyyy_0, \
                             g_0_xxyyyz_0_xxxxxyyz_0, \
                             g_0_xxyyyz_0_xxxxxyzz_0, \
                             g_0_xxyyyz_0_xxxxxzzz_0, \
                             g_0_xxyyyz_0_xxxxyyyy_0, \
                             g_0_xxyyyz_0_xxxxyyyz_0, \
                             g_0_xxyyyz_0_xxxxyyzz_0, \
                             g_0_xxyyyz_0_xxxxyzzz_0, \
                             g_0_xxyyyz_0_xxxxzzzz_0, \
                             g_0_xxyyyz_0_xxxyyyyy_0, \
                             g_0_xxyyyz_0_xxxyyyyz_0, \
                             g_0_xxyyyz_0_xxxyyyzz_0, \
                             g_0_xxyyyz_0_xxxyyzzz_0, \
                             g_0_xxyyyz_0_xxxyzzzz_0, \
                             g_0_xxyyyz_0_xxxzzzzz_0, \
                             g_0_xxyyyz_0_xxyyyyyy_0, \
                             g_0_xxyyyz_0_xxyyyyyz_0, \
                             g_0_xxyyyz_0_xxyyyyzz_0, \
                             g_0_xxyyyz_0_xxyyyzzz_0, \
                             g_0_xxyyyz_0_xxyyzzzz_0, \
                             g_0_xxyyyz_0_xxyzzzzz_0, \
                             g_0_xxyyyz_0_xxzzzzzz_0, \
                             g_0_xxyyyz_0_xyyyyyyy_0, \
                             g_0_xxyyyz_0_xyyyyyyz_0, \
                             g_0_xxyyyz_0_xyyyyyzz_0, \
                             g_0_xxyyyz_0_xyyyyzzz_0, \
                             g_0_xxyyyz_0_xyyyzzzz_0, \
                             g_0_xxyyyz_0_xyyzzzzz_0, \
                             g_0_xxyyyz_0_xyzzzzzz_0, \
                             g_0_xxyyyz_0_xzzzzzzz_0, \
                             g_0_xxyyyz_0_yyyyyyyy_0, \
                             g_0_xxyyyz_0_yyyyyyyz_0, \
                             g_0_xxyyyz_0_yyyyyyzz_0, \
                             g_0_xxyyyz_0_yyyyyzzz_0, \
                             g_0_xxyyyz_0_yyyyzzzz_0, \
                             g_0_xxyyyz_0_yyyzzzzz_0, \
                             g_0_xxyyyz_0_yyzzzzzz_0, \
                             g_0_xxyyyz_0_yzzzzzzz_0, \
                             g_0_xxyyyz_0_zzzzzzzz_0, \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyz_0_xxxxxxxx_0[i] = g_0_xxyyy_0_xxxxxxxx_0[i] * pb_z + g_0_xxyyy_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxxxy_0[i] = g_0_xxyyy_0_xxxxxxxy_0[i] * pb_z + g_0_xxyyy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxxxz_0[i] = g_0_xxyyy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxxxz_0[i] * pb_z + g_0_xxyyy_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxxyy_0[i] = g_0_xxyyy_0_xxxxxxyy_0[i] * pb_z + g_0_xxyyy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxxyz_0[i] = g_0_xxyyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxxyz_0[i] * pb_z + g_0_xxyyy_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxxzz_0[i] =
            2.0 * g_0_xxyyy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxxzz_0[i] * pb_z + g_0_xxyyy_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxyyy_0[i] = g_0_xxyyy_0_xxxxxyyy_0[i] * pb_z + g_0_xxyyy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxyyz_0[i] = g_0_xxyyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxyyz_0[i] * pb_z + g_0_xxyyy_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxyzz_0[i] =
            2.0 * g_0_xxyyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxyzz_0[i] * pb_z + g_0_xxyyy_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxzzz_0[i] =
            3.0 * g_0_xxyyy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxzzz_0[i] * pb_z + g_0_xxyyy_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxyyyy_0[i] = g_0_xxyyy_0_xxxxyyyy_0[i] * pb_z + g_0_xxyyy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxyyyz_0[i] = g_0_xxyyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyyyz_0[i] * pb_z + g_0_xxyyy_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxyyzz_0[i] =
            2.0 * g_0_xxyyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyyzz_0[i] * pb_z + g_0_xxyyy_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxyzzz_0[i] =
            3.0 * g_0_xxyyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyzzz_0[i] * pb_z + g_0_xxyyy_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxzzzz_0[i] =
            4.0 * g_0_xxyyy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxzzzz_0[i] * pb_z + g_0_xxyyy_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyyyyy_0[i] = g_0_xxyyy_0_xxxyyyyy_0[i] * pb_z + g_0_xxyyy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyyyyz_0[i] = g_0_xxyyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyyyz_0[i] * pb_z + g_0_xxyyy_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyyyzz_0[i] =
            2.0 * g_0_xxyyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyyzz_0[i] * pb_z + g_0_xxyyy_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyyzzz_0[i] =
            3.0 * g_0_xxyyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyzzz_0[i] * pb_z + g_0_xxyyy_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyzzzz_0[i] =
            4.0 * g_0_xxyyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyzzzz_0[i] * pb_z + g_0_xxyyy_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxzzzzz_0[i] =
            5.0 * g_0_xxyyy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxzzzzz_0[i] * pb_z + g_0_xxyyy_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyyyyy_0[i] = g_0_xxyyy_0_xxyyyyyy_0[i] * pb_z + g_0_xxyyy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyyyyz_0[i] = g_0_xxyyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyyyz_0[i] * pb_z + g_0_xxyyy_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyyyzz_0[i] =
            2.0 * g_0_xxyyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyyzz_0[i] * pb_z + g_0_xxyyy_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyyzzz_0[i] =
            3.0 * g_0_xxyyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyzzz_0[i] * pb_z + g_0_xxyyy_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyzzzz_0[i] =
            4.0 * g_0_xxyyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyzzzz_0[i] * pb_z + g_0_xxyyy_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyzzzzz_0[i] =
            5.0 * g_0_xxyyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyzzzzz_0[i] * pb_z + g_0_xxyyy_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxzzzzzz_0[i] =
            6.0 * g_0_xxyyy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxzzzzzz_0[i] * pb_z + g_0_xxyyy_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyyyyy_0[i] = g_0_xxyyy_0_xyyyyyyy_0[i] * pb_z + g_0_xxyyy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyyyyz_0[i] = g_0_xxyyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyyyz_0[i] * pb_z + g_0_xxyyy_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyyyzz_0[i] =
            2.0 * g_0_xxyyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyyzz_0[i] * pb_z + g_0_xxyyy_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyyzzz_0[i] =
            3.0 * g_0_xxyyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyzzz_0[i] * pb_z + g_0_xxyyy_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyzzzz_0[i] =
            4.0 * g_0_xxyyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyzzzz_0[i] * pb_z + g_0_xxyyy_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyzzzzz_0[i] =
            5.0 * g_0_xxyyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyzzzzz_0[i] * pb_z + g_0_xxyyy_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyzzzzzz_0[i] =
            6.0 * g_0_xxyyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyzzzzzz_0[i] * pb_z + g_0_xxyyy_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xzzzzzzz_0[i] =
            7.0 * g_0_xxyyy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xzzzzzzz_0[i] * pb_z + g_0_xxyyy_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyyyyy_0[i] = g_0_xxyyy_0_yyyyyyyy_0[i] * pb_z + g_0_xxyyy_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyyyyz_0[i] = g_0_xxyyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyyyyyz_0[i] * pb_z + g_0_xxyyy_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyyyzz_0[i] =
            2.0 * g_0_xxyyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyyyyzz_0[i] * pb_z + g_0_xxyyy_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyyzzz_0[i] =
            3.0 * g_0_xxyyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyyyzzz_0[i] * pb_z + g_0_xxyyy_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyzzzz_0[i] =
            4.0 * g_0_xxyyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyyzzzz_0[i] * pb_z + g_0_xxyyy_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyzzzzz_0[i] =
            5.0 * g_0_xxyyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyzzzzz_0[i] * pb_z + g_0_xxyyy_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyzzzzzz_0[i] =
            6.0 * g_0_xxyyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyzzzzzz_0[i] * pb_z + g_0_xxyyy_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yzzzzzzz_0[i] =
            7.0 * g_0_xxyyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yzzzzzzz_0[i] * pb_z + g_0_xxyyy_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_zzzzzzzz_0[i] =
            8.0 * g_0_xxyyy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_zzzzzzzz_0[i] * pb_z + g_0_xxyyy_0_zzzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 540-585 components of targeted buffer : SISL

    auto g_0_xxyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 540);

    auto g_0_xxyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 541);

    auto g_0_xxyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 542);

    auto g_0_xxyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 543);

    auto g_0_xxyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 544);

    auto g_0_xxyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 545);

    auto g_0_xxyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 546);

    auto g_0_xxyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 547);

    auto g_0_xxyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 548);

    auto g_0_xxyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 549);

    auto g_0_xxyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 550);

    auto g_0_xxyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 551);

    auto g_0_xxyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 552);

    auto g_0_xxyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 553);

    auto g_0_xxyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 554);

    auto g_0_xxyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 555);

    auto g_0_xxyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 556);

    auto g_0_xxyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 557);

    auto g_0_xxyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 558);

    auto g_0_xxyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 559);

    auto g_0_xxyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 560);

    auto g_0_xxyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 561);

    auto g_0_xxyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 562);

    auto g_0_xxyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 563);

    auto g_0_xxyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 564);

    auto g_0_xxyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 565);

    auto g_0_xxyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 566);

    auto g_0_xxyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 567);

    auto g_0_xxyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 568);

    auto g_0_xxyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 569);

    auto g_0_xxyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 570);

    auto g_0_xxyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 571);

    auto g_0_xxyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 572);

    auto g_0_xxyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 573);

    auto g_0_xxyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 574);

    auto g_0_xxyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 575);

    auto g_0_xxyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 576);

    auto g_0_xxyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 577);

    auto g_0_xxyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 578);

    auto g_0_xxyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 579);

    auto g_0_xxyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 580);

    auto g_0_xxyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 581);

    auto g_0_xxyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 582);

    auto g_0_xxyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 583);

    auto g_0_xxyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 584);

#pragma omp simd aligned(g_0_xxyy_0_xxxxxxxy_0,       \
                             g_0_xxyy_0_xxxxxxxy_1,   \
                             g_0_xxyy_0_xxxxxxyy_0,   \
                             g_0_xxyy_0_xxxxxxyy_1,   \
                             g_0_xxyy_0_xxxxxyyy_0,   \
                             g_0_xxyy_0_xxxxxyyy_1,   \
                             g_0_xxyy_0_xxxxyyyy_0,   \
                             g_0_xxyy_0_xxxxyyyy_1,   \
                             g_0_xxyy_0_xxxyyyyy_0,   \
                             g_0_xxyy_0_xxxyyyyy_1,   \
                             g_0_xxyy_0_xxyyyyyy_0,   \
                             g_0_xxyy_0_xxyyyyyy_1,   \
                             g_0_xxyy_0_xyyyyyyy_0,   \
                             g_0_xxyy_0_xyyyyyyy_1,   \
                             g_0_xxyyz_0_xxxxxxxy_0,  \
                             g_0_xxyyz_0_xxxxxxxy_1,  \
                             g_0_xxyyz_0_xxxxxxyy_0,  \
                             g_0_xxyyz_0_xxxxxxyy_1,  \
                             g_0_xxyyz_0_xxxxxyyy_0,  \
                             g_0_xxyyz_0_xxxxxyyy_1,  \
                             g_0_xxyyz_0_xxxxyyyy_0,  \
                             g_0_xxyyz_0_xxxxyyyy_1,  \
                             g_0_xxyyz_0_xxxyyyyy_0,  \
                             g_0_xxyyz_0_xxxyyyyy_1,  \
                             g_0_xxyyz_0_xxyyyyyy_0,  \
                             g_0_xxyyz_0_xxyyyyyy_1,  \
                             g_0_xxyyz_0_xyyyyyyy_0,  \
                             g_0_xxyyz_0_xyyyyyyy_1,  \
                             g_0_xxyyzz_0_xxxxxxxx_0, \
                             g_0_xxyyzz_0_xxxxxxxy_0, \
                             g_0_xxyyzz_0_xxxxxxxz_0, \
                             g_0_xxyyzz_0_xxxxxxyy_0, \
                             g_0_xxyyzz_0_xxxxxxyz_0, \
                             g_0_xxyyzz_0_xxxxxxzz_0, \
                             g_0_xxyyzz_0_xxxxxyyy_0, \
                             g_0_xxyyzz_0_xxxxxyyz_0, \
                             g_0_xxyyzz_0_xxxxxyzz_0, \
                             g_0_xxyyzz_0_xxxxxzzz_0, \
                             g_0_xxyyzz_0_xxxxyyyy_0, \
                             g_0_xxyyzz_0_xxxxyyyz_0, \
                             g_0_xxyyzz_0_xxxxyyzz_0, \
                             g_0_xxyyzz_0_xxxxyzzz_0, \
                             g_0_xxyyzz_0_xxxxzzzz_0, \
                             g_0_xxyyzz_0_xxxyyyyy_0, \
                             g_0_xxyyzz_0_xxxyyyyz_0, \
                             g_0_xxyyzz_0_xxxyyyzz_0, \
                             g_0_xxyyzz_0_xxxyyzzz_0, \
                             g_0_xxyyzz_0_xxxyzzzz_0, \
                             g_0_xxyyzz_0_xxxzzzzz_0, \
                             g_0_xxyyzz_0_xxyyyyyy_0, \
                             g_0_xxyyzz_0_xxyyyyyz_0, \
                             g_0_xxyyzz_0_xxyyyyzz_0, \
                             g_0_xxyyzz_0_xxyyyzzz_0, \
                             g_0_xxyyzz_0_xxyyzzzz_0, \
                             g_0_xxyyzz_0_xxyzzzzz_0, \
                             g_0_xxyyzz_0_xxzzzzzz_0, \
                             g_0_xxyyzz_0_xyyyyyyy_0, \
                             g_0_xxyyzz_0_xyyyyyyz_0, \
                             g_0_xxyyzz_0_xyyyyyzz_0, \
                             g_0_xxyyzz_0_xyyyyzzz_0, \
                             g_0_xxyyzz_0_xyyyzzzz_0, \
                             g_0_xxyyzz_0_xyyzzzzz_0, \
                             g_0_xxyyzz_0_xyzzzzzz_0, \
                             g_0_xxyyzz_0_xzzzzzzz_0, \
                             g_0_xxyyzz_0_yyyyyyyy_0, \
                             g_0_xxyyzz_0_yyyyyyyz_0, \
                             g_0_xxyyzz_0_yyyyyyzz_0, \
                             g_0_xxyyzz_0_yyyyyzzz_0, \
                             g_0_xxyyzz_0_yyyyzzzz_0, \
                             g_0_xxyyzz_0_yyyzzzzz_0, \
                             g_0_xxyyzz_0_yyzzzzzz_0, \
                             g_0_xxyyzz_0_yzzzzzzz_0, \
                             g_0_xxyyzz_0_zzzzzzzz_0, \
                             g_0_xxyzz_0_xxxxxxxx_0,  \
                             g_0_xxyzz_0_xxxxxxxx_1,  \
                             g_0_xxyzz_0_xxxxxxxz_0,  \
                             g_0_xxyzz_0_xxxxxxxz_1,  \
                             g_0_xxyzz_0_xxxxxxzz_0,  \
                             g_0_xxyzz_0_xxxxxxzz_1,  \
                             g_0_xxyzz_0_xxxxxzzz_0,  \
                             g_0_xxyzz_0_xxxxxzzz_1,  \
                             g_0_xxyzz_0_xxxxzzzz_0,  \
                             g_0_xxyzz_0_xxxxzzzz_1,  \
                             g_0_xxyzz_0_xxxzzzzz_0,  \
                             g_0_xxyzz_0_xxxzzzzz_1,  \
                             g_0_xxyzz_0_xxzzzzzz_0,  \
                             g_0_xxyzz_0_xxzzzzzz_1,  \
                             g_0_xxyzz_0_xzzzzzzz_0,  \
                             g_0_xxyzz_0_xzzzzzzz_1,  \
                             g_0_xxzz_0_xxxxxxxx_0,   \
                             g_0_xxzz_0_xxxxxxxx_1,   \
                             g_0_xxzz_0_xxxxxxxz_0,   \
                             g_0_xxzz_0_xxxxxxxz_1,   \
                             g_0_xxzz_0_xxxxxxzz_0,   \
                             g_0_xxzz_0_xxxxxxzz_1,   \
                             g_0_xxzz_0_xxxxxzzz_0,   \
                             g_0_xxzz_0_xxxxxzzz_1,   \
                             g_0_xxzz_0_xxxxzzzz_0,   \
                             g_0_xxzz_0_xxxxzzzz_1,   \
                             g_0_xxzz_0_xxxzzzzz_0,   \
                             g_0_xxzz_0_xxxzzzzz_1,   \
                             g_0_xxzz_0_xxzzzzzz_0,   \
                             g_0_xxzz_0_xxzzzzzz_1,   \
                             g_0_xxzz_0_xzzzzzzz_0,   \
                             g_0_xxzz_0_xzzzzzzz_1,   \
                             g_0_xyyzz_0_xxxxxxyz_0,  \
                             g_0_xyyzz_0_xxxxxxyz_1,  \
                             g_0_xyyzz_0_xxxxxyyz_0,  \
                             g_0_xyyzz_0_xxxxxyyz_1,  \
                             g_0_xyyzz_0_xxxxxyz_1,   \
                             g_0_xyyzz_0_xxxxxyzz_0,  \
                             g_0_xyyzz_0_xxxxxyzz_1,  \
                             g_0_xyyzz_0_xxxxyyyz_0,  \
                             g_0_xyyzz_0_xxxxyyyz_1,  \
                             g_0_xyyzz_0_xxxxyyz_1,   \
                             g_0_xyyzz_0_xxxxyyzz_0,  \
                             g_0_xyyzz_0_xxxxyyzz_1,  \
                             g_0_xyyzz_0_xxxxyzz_1,   \
                             g_0_xyyzz_0_xxxxyzzz_0,  \
                             g_0_xyyzz_0_xxxxyzzz_1,  \
                             g_0_xyyzz_0_xxxyyyyz_0,  \
                             g_0_xyyzz_0_xxxyyyyz_1,  \
                             g_0_xyyzz_0_xxxyyyz_1,   \
                             g_0_xyyzz_0_xxxyyyzz_0,  \
                             g_0_xyyzz_0_xxxyyyzz_1,  \
                             g_0_xyyzz_0_xxxyyzz_1,   \
                             g_0_xyyzz_0_xxxyyzzz_0,  \
                             g_0_xyyzz_0_xxxyyzzz_1,  \
                             g_0_xyyzz_0_xxxyzzz_1,   \
                             g_0_xyyzz_0_xxxyzzzz_0,  \
                             g_0_xyyzz_0_xxxyzzzz_1,  \
                             g_0_xyyzz_0_xxyyyyyz_0,  \
                             g_0_xyyzz_0_xxyyyyyz_1,  \
                             g_0_xyyzz_0_xxyyyyz_1,   \
                             g_0_xyyzz_0_xxyyyyzz_0,  \
                             g_0_xyyzz_0_xxyyyyzz_1,  \
                             g_0_xyyzz_0_xxyyyzz_1,   \
                             g_0_xyyzz_0_xxyyyzzz_0,  \
                             g_0_xyyzz_0_xxyyyzzz_1,  \
                             g_0_xyyzz_0_xxyyzzz_1,   \
                             g_0_xyyzz_0_xxyyzzzz_0,  \
                             g_0_xyyzz_0_xxyyzzzz_1,  \
                             g_0_xyyzz_0_xxyzzzz_1,   \
                             g_0_xyyzz_0_xxyzzzzz_0,  \
                             g_0_xyyzz_0_xxyzzzzz_1,  \
                             g_0_xyyzz_0_xyyyyyyz_0,  \
                             g_0_xyyzz_0_xyyyyyyz_1,  \
                             g_0_xyyzz_0_xyyyyyz_1,   \
                             g_0_xyyzz_0_xyyyyyzz_0,  \
                             g_0_xyyzz_0_xyyyyyzz_1,  \
                             g_0_xyyzz_0_xyyyyzz_1,   \
                             g_0_xyyzz_0_xyyyyzzz_0,  \
                             g_0_xyyzz_0_xyyyyzzz_1,  \
                             g_0_xyyzz_0_xyyyzzz_1,   \
                             g_0_xyyzz_0_xyyyzzzz_0,  \
                             g_0_xyyzz_0_xyyyzzzz_1,  \
                             g_0_xyyzz_0_xyyzzzz_1,   \
                             g_0_xyyzz_0_xyyzzzzz_0,  \
                             g_0_xyyzz_0_xyyzzzzz_1,  \
                             g_0_xyyzz_0_xyzzzzz_1,   \
                             g_0_xyyzz_0_xyzzzzzz_0,  \
                             g_0_xyyzz_0_xyzzzzzz_1,  \
                             g_0_xyyzz_0_yyyyyyyy_0,  \
                             g_0_xyyzz_0_yyyyyyyy_1,  \
                             g_0_xyyzz_0_yyyyyyyz_0,  \
                             g_0_xyyzz_0_yyyyyyyz_1,  \
                             g_0_xyyzz_0_yyyyyyz_1,   \
                             g_0_xyyzz_0_yyyyyyzz_0,  \
                             g_0_xyyzz_0_yyyyyyzz_1,  \
                             g_0_xyyzz_0_yyyyyzz_1,   \
                             g_0_xyyzz_0_yyyyyzzz_0,  \
                             g_0_xyyzz_0_yyyyyzzz_1,  \
                             g_0_xyyzz_0_yyyyzzz_1,   \
                             g_0_xyyzz_0_yyyyzzzz_0,  \
                             g_0_xyyzz_0_yyyyzzzz_1,  \
                             g_0_xyyzz_0_yyyzzzz_1,   \
                             g_0_xyyzz_0_yyyzzzzz_0,  \
                             g_0_xyyzz_0_yyyzzzzz_1,  \
                             g_0_xyyzz_0_yyzzzzz_1,   \
                             g_0_xyyzz_0_yyzzzzzz_0,  \
                             g_0_xyyzz_0_yyzzzzzz_1,  \
                             g_0_xyyzz_0_yzzzzzz_1,   \
                             g_0_xyyzz_0_yzzzzzzz_0,  \
                             g_0_xyyzz_0_yzzzzzzz_1,  \
                             g_0_xyyzz_0_zzzzzzzz_0,  \
                             g_0_xyyzz_0_zzzzzzzz_1,  \
                             g_0_yyzz_0_xxxxxxyz_0,   \
                             g_0_yyzz_0_xxxxxxyz_1,   \
                             g_0_yyzz_0_xxxxxyyz_0,   \
                             g_0_yyzz_0_xxxxxyyz_1,   \
                             g_0_yyzz_0_xxxxxyzz_0,   \
                             g_0_yyzz_0_xxxxxyzz_1,   \
                             g_0_yyzz_0_xxxxyyyz_0,   \
                             g_0_yyzz_0_xxxxyyyz_1,   \
                             g_0_yyzz_0_xxxxyyzz_0,   \
                             g_0_yyzz_0_xxxxyyzz_1,   \
                             g_0_yyzz_0_xxxxyzzz_0,   \
                             g_0_yyzz_0_xxxxyzzz_1,   \
                             g_0_yyzz_0_xxxyyyyz_0,   \
                             g_0_yyzz_0_xxxyyyyz_1,   \
                             g_0_yyzz_0_xxxyyyzz_0,   \
                             g_0_yyzz_0_xxxyyyzz_1,   \
                             g_0_yyzz_0_xxxyyzzz_0,   \
                             g_0_yyzz_0_xxxyyzzz_1,   \
                             g_0_yyzz_0_xxxyzzzz_0,   \
                             g_0_yyzz_0_xxxyzzzz_1,   \
                             g_0_yyzz_0_xxyyyyyz_0,   \
                             g_0_yyzz_0_xxyyyyyz_1,   \
                             g_0_yyzz_0_xxyyyyzz_0,   \
                             g_0_yyzz_0_xxyyyyzz_1,   \
                             g_0_yyzz_0_xxyyyzzz_0,   \
                             g_0_yyzz_0_xxyyyzzz_1,   \
                             g_0_yyzz_0_xxyyzzzz_0,   \
                             g_0_yyzz_0_xxyyzzzz_1,   \
                             g_0_yyzz_0_xxyzzzzz_0,   \
                             g_0_yyzz_0_xxyzzzzz_1,   \
                             g_0_yyzz_0_xyyyyyyz_0,   \
                             g_0_yyzz_0_xyyyyyyz_1,   \
                             g_0_yyzz_0_xyyyyyzz_0,   \
                             g_0_yyzz_0_xyyyyyzz_1,   \
                             g_0_yyzz_0_xyyyyzzz_0,   \
                             g_0_yyzz_0_xyyyyzzz_1,   \
                             g_0_yyzz_0_xyyyzzzz_0,   \
                             g_0_yyzz_0_xyyyzzzz_1,   \
                             g_0_yyzz_0_xyyzzzzz_0,   \
                             g_0_yyzz_0_xyyzzzzz_1,   \
                             g_0_yyzz_0_xyzzzzzz_0,   \
                             g_0_yyzz_0_xyzzzzzz_1,   \
                             g_0_yyzz_0_yyyyyyyy_0,   \
                             g_0_yyzz_0_yyyyyyyy_1,   \
                             g_0_yyzz_0_yyyyyyyz_0,   \
                             g_0_yyzz_0_yyyyyyyz_1,   \
                             g_0_yyzz_0_yyyyyyzz_0,   \
                             g_0_yyzz_0_yyyyyyzz_1,   \
                             g_0_yyzz_0_yyyyyzzz_0,   \
                             g_0_yyzz_0_yyyyyzzz_1,   \
                             g_0_yyzz_0_yyyyzzzz_0,   \
                             g_0_yyzz_0_yyyyzzzz_1,   \
                             g_0_yyzz_0_yyyzzzzz_0,   \
                             g_0_yyzz_0_yyyzzzzz_1,   \
                             g_0_yyzz_0_yyzzzzzz_0,   \
                             g_0_yyzz_0_yyzzzzzz_1,   \
                             g_0_yyzz_0_yzzzzzzz_0,   \
                             g_0_yyzz_0_yzzzzzzz_1,   \
                             g_0_yyzz_0_zzzzzzzz_0,   \
                             g_0_yyzz_0_zzzzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzz_0_xxxxxxxx_0[i] = g_0_xxzz_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxxxxx_0[i] * pb_y +
                                     g_0_xxyzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxxxxxy_0[i] = g_0_xxyy_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxxxxxy_0[i] * pb_z +
                                     g_0_xxyyz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxxxxxz_0[i] = g_0_xxzz_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxxxxz_0[i] * pb_y +
                                     g_0_xxyzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxxxxyy_0[i] = g_0_xxyy_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxxxxyy_0[i] * pb_z +
                                     g_0_xxyyz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxxxxyz_0[i] = g_0_yyzz_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxxxxyz_1[i] * fti_ab_0 +
                                     6.0 * g_0_xyyzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxxxxyz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxxxxzz_0[i] = g_0_xxzz_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxxxzz_0[i] * pb_y +
                                     g_0_xxyzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxxxyyy_0[i] = g_0_xxyy_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxxxyyy_0[i] * pb_z +
                                     g_0_xxyyz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxxxyyz_0[i] = g_0_yyzz_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xyyzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxxxyyz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxxxyzz_0[i] = g_0_yyzz_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxxxyzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xyyzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxxxyzz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxxxzzz_0[i] = g_0_xxzz_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxxzzz_0[i] * pb_y +
                                     g_0_xxyzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxxyyyy_0[i] = g_0_xxyy_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxxyyyy_0[i] * pb_z +
                                     g_0_xxyyz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxxyyyz_0[i] = g_0_yyzz_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xyyzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxxyyyz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxxyyzz_0[i] = g_0_yyzz_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xyyzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxxyyzz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxxyzzz_0[i] = g_0_yyzz_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxxyzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xyyzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxxyzzz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxxzzzz_0[i] = g_0_xxzz_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxzzzz_0[i] * pb_y +
                                     g_0_xxyzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxyyyyy_0[i] = g_0_xxyy_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxyyyyy_0[i] * pb_z +
                                     g_0_xxyyz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxyyyyz_0[i] = g_0_yyzz_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxyyyyz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxyyyzz_0[i] = g_0_yyzz_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxyyyzz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxyyzzz_0[i] = g_0_yyzz_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxyyzzz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxyzzzz_0[i] = g_0_yyzz_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxyzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxyzzzz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxzzzzz_0[i] = g_0_xxzz_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxzzzzz_0[i] * pb_y +
                                     g_0_xxyzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxyyyyyy_0[i] = g_0_xxyy_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxyyyyyy_0[i] * pb_z +
                                     g_0_xxyyz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxyyyyyz_0[i] = g_0_yyzz_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxyyyyyz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxyyyyzz_0[i] = g_0_yyzz_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxyyyyzz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxyyyzzz_0[i] = g_0_yyzz_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxyyyzzz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxyyzzzz_0[i] = g_0_yyzz_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxyyzzzz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxyzzzzz_0[i] = g_0_yyzz_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxyzzzzz_0[i] * pb_x +
                                     g_0_xyyzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxzzzzzz_0[i] = g_0_xxzz_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxzzzzzz_0[i] * pb_y +
                                     g_0_xxyzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xyyyyyyy_0[i] = g_0_xxyy_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xyyyyyyy_0[i] * pb_z +
                                     g_0_xxyyz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xyyyyyyz_0[i] = g_0_yyzz_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyyyz_1[i] * fi_abcd_0 +
                                     g_0_xyyzz_0_xyyyyyyz_0[i] * pb_x + g_0_xyyzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyyyyyzz_0[i] = g_0_yyzz_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyyyyzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyyzz_1[i] * fi_abcd_0 +
                                     g_0_xyyzz_0_xyyyyyzz_0[i] * pb_x + g_0_xyyzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyyyyzzz_0[i] = g_0_yyzz_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyyyzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyzz_0_xyyyyzzz_0[i] * pb_x + g_0_xyyzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyyyzzzz_0[i] = g_0_yyzz_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyyzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyzzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyzz_0_xyyyzzzz_0[i] * pb_x + g_0_xyyzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyyzzzzz_0[i] = g_0_yyzz_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyzzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyzzzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyzz_0_xyyzzzzz_0[i] * pb_x + g_0_xyyzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyzzzzzz_0[i] = g_0_yyzz_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yzzzzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyzz_0_xyzzzzzz_0[i] * pb_x + g_0_xyyzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xzzzzzzz_0[i] = g_0_xxzz_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xzzzzzzz_0[i] * pb_y +
                                     g_0_xxyzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_yyyyyyyy_0[i] = g_0_yyzz_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyyyyy_0[i] * pb_x +
                                     g_0_xyyzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyyyyyz_0[i] = g_0_yyzz_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyyyyz_0[i] * pb_x +
                                     g_0_xyyzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyyyyzz_0[i] = g_0_yyzz_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyyyzz_0[i] * pb_x +
                                     g_0_xyyzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyyyzzz_0[i] = g_0_yyzz_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyyzzz_0[i] * pb_x +
                                     g_0_xyyzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyyzzzz_0[i] = g_0_yyzz_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyzzzz_0[i] * pb_x +
                                     g_0_xyyzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyzzzzz_0[i] = g_0_yyzz_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyzzzzz_0[i] * pb_x +
                                     g_0_xyyzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyzzzzzz_0[i] = g_0_yyzz_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyzzzzzz_0[i] * pb_x +
                                     g_0_xyyzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yzzzzzzz_0[i] = g_0_yyzz_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yzzzzzzz_0[i] * pb_x +
                                     g_0_xyyzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_zzzzzzzz_0[i] = g_0_yyzz_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_zzzzzzzz_0[i] * pb_x +
                                     g_0_xyyzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 585-630 components of targeted buffer : SISL

    auto g_0_xxyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 585);

    auto g_0_xxyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 586);

    auto g_0_xxyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 587);

    auto g_0_xxyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 588);

    auto g_0_xxyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 589);

    auto g_0_xxyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 590);

    auto g_0_xxyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 591);

    auto g_0_xxyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 592);

    auto g_0_xxyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 593);

    auto g_0_xxyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 594);

    auto g_0_xxyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 595);

    auto g_0_xxyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 596);

    auto g_0_xxyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 597);

    auto g_0_xxyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 598);

    auto g_0_xxyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 599);

    auto g_0_xxyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 600);

    auto g_0_xxyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 601);

    auto g_0_xxyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 602);

    auto g_0_xxyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 603);

    auto g_0_xxyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 604);

    auto g_0_xxyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 605);

    auto g_0_xxyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 606);

    auto g_0_xxyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 607);

    auto g_0_xxyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 608);

    auto g_0_xxyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 609);

    auto g_0_xxyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 610);

    auto g_0_xxyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 611);

    auto g_0_xxyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 612);

    auto g_0_xxyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 613);

    auto g_0_xxyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 614);

    auto g_0_xxyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 615);

    auto g_0_xxyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 616);

    auto g_0_xxyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 617);

    auto g_0_xxyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 618);

    auto g_0_xxyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 619);

    auto g_0_xxyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 620);

    auto g_0_xxyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 621);

    auto g_0_xxyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 622);

    auto g_0_xxyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 623);

    auto g_0_xxyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 624);

    auto g_0_xxyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 625);

    auto g_0_xxyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 626);

    auto g_0_xxyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 627);

    auto g_0_xxyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 628);

    auto g_0_xxyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 629);

#pragma omp simd aligned(g_0_xxyzzz_0_xxxxxxxx_0,     \
                             g_0_xxyzzz_0_xxxxxxxy_0, \
                             g_0_xxyzzz_0_xxxxxxxz_0, \
                             g_0_xxyzzz_0_xxxxxxyy_0, \
                             g_0_xxyzzz_0_xxxxxxyz_0, \
                             g_0_xxyzzz_0_xxxxxxzz_0, \
                             g_0_xxyzzz_0_xxxxxyyy_0, \
                             g_0_xxyzzz_0_xxxxxyyz_0, \
                             g_0_xxyzzz_0_xxxxxyzz_0, \
                             g_0_xxyzzz_0_xxxxxzzz_0, \
                             g_0_xxyzzz_0_xxxxyyyy_0, \
                             g_0_xxyzzz_0_xxxxyyyz_0, \
                             g_0_xxyzzz_0_xxxxyyzz_0, \
                             g_0_xxyzzz_0_xxxxyzzz_0, \
                             g_0_xxyzzz_0_xxxxzzzz_0, \
                             g_0_xxyzzz_0_xxxyyyyy_0, \
                             g_0_xxyzzz_0_xxxyyyyz_0, \
                             g_0_xxyzzz_0_xxxyyyzz_0, \
                             g_0_xxyzzz_0_xxxyyzzz_0, \
                             g_0_xxyzzz_0_xxxyzzzz_0, \
                             g_0_xxyzzz_0_xxxzzzzz_0, \
                             g_0_xxyzzz_0_xxyyyyyy_0, \
                             g_0_xxyzzz_0_xxyyyyyz_0, \
                             g_0_xxyzzz_0_xxyyyyzz_0, \
                             g_0_xxyzzz_0_xxyyyzzz_0, \
                             g_0_xxyzzz_0_xxyyzzzz_0, \
                             g_0_xxyzzz_0_xxyzzzzz_0, \
                             g_0_xxyzzz_0_xxzzzzzz_0, \
                             g_0_xxyzzz_0_xyyyyyyy_0, \
                             g_0_xxyzzz_0_xyyyyyyz_0, \
                             g_0_xxyzzz_0_xyyyyyzz_0, \
                             g_0_xxyzzz_0_xyyyyzzz_0, \
                             g_0_xxyzzz_0_xyyyzzzz_0, \
                             g_0_xxyzzz_0_xyyzzzzz_0, \
                             g_0_xxyzzz_0_xyzzzzzz_0, \
                             g_0_xxyzzz_0_xzzzzzzz_0, \
                             g_0_xxyzzz_0_yyyyyyyy_0, \
                             g_0_xxyzzz_0_yyyyyyyz_0, \
                             g_0_xxyzzz_0_yyyyyyzz_0, \
                             g_0_xxyzzz_0_yyyyyzzz_0, \
                             g_0_xxyzzz_0_yyyyzzzz_0, \
                             g_0_xxyzzz_0_yyyzzzzz_0, \
                             g_0_xxyzzz_0_yyzzzzzz_0, \
                             g_0_xxyzzz_0_yzzzzzzz_0, \
                             g_0_xxyzzz_0_zzzzzzzz_0, \
                             g_0_xxzzz_0_xxxxxxx_1,   \
                             g_0_xxzzz_0_xxxxxxxx_0,  \
                             g_0_xxzzz_0_xxxxxxxx_1,  \
                             g_0_xxzzz_0_xxxxxxxy_0,  \
                             g_0_xxzzz_0_xxxxxxxy_1,  \
                             g_0_xxzzz_0_xxxxxxxz_0,  \
                             g_0_xxzzz_0_xxxxxxxz_1,  \
                             g_0_xxzzz_0_xxxxxxy_1,   \
                             g_0_xxzzz_0_xxxxxxyy_0,  \
                             g_0_xxzzz_0_xxxxxxyy_1,  \
                             g_0_xxzzz_0_xxxxxxyz_0,  \
                             g_0_xxzzz_0_xxxxxxyz_1,  \
                             g_0_xxzzz_0_xxxxxxz_1,   \
                             g_0_xxzzz_0_xxxxxxzz_0,  \
                             g_0_xxzzz_0_xxxxxxzz_1,  \
                             g_0_xxzzz_0_xxxxxyy_1,   \
                             g_0_xxzzz_0_xxxxxyyy_0,  \
                             g_0_xxzzz_0_xxxxxyyy_1,  \
                             g_0_xxzzz_0_xxxxxyyz_0,  \
                             g_0_xxzzz_0_xxxxxyyz_1,  \
                             g_0_xxzzz_0_xxxxxyz_1,   \
                             g_0_xxzzz_0_xxxxxyzz_0,  \
                             g_0_xxzzz_0_xxxxxyzz_1,  \
                             g_0_xxzzz_0_xxxxxzz_1,   \
                             g_0_xxzzz_0_xxxxxzzz_0,  \
                             g_0_xxzzz_0_xxxxxzzz_1,  \
                             g_0_xxzzz_0_xxxxyyy_1,   \
                             g_0_xxzzz_0_xxxxyyyy_0,  \
                             g_0_xxzzz_0_xxxxyyyy_1,  \
                             g_0_xxzzz_0_xxxxyyyz_0,  \
                             g_0_xxzzz_0_xxxxyyyz_1,  \
                             g_0_xxzzz_0_xxxxyyz_1,   \
                             g_0_xxzzz_0_xxxxyyzz_0,  \
                             g_0_xxzzz_0_xxxxyyzz_1,  \
                             g_0_xxzzz_0_xxxxyzz_1,   \
                             g_0_xxzzz_0_xxxxyzzz_0,  \
                             g_0_xxzzz_0_xxxxyzzz_1,  \
                             g_0_xxzzz_0_xxxxzzz_1,   \
                             g_0_xxzzz_0_xxxxzzzz_0,  \
                             g_0_xxzzz_0_xxxxzzzz_1,  \
                             g_0_xxzzz_0_xxxyyyy_1,   \
                             g_0_xxzzz_0_xxxyyyyy_0,  \
                             g_0_xxzzz_0_xxxyyyyy_1,  \
                             g_0_xxzzz_0_xxxyyyyz_0,  \
                             g_0_xxzzz_0_xxxyyyyz_1,  \
                             g_0_xxzzz_0_xxxyyyz_1,   \
                             g_0_xxzzz_0_xxxyyyzz_0,  \
                             g_0_xxzzz_0_xxxyyyzz_1,  \
                             g_0_xxzzz_0_xxxyyzz_1,   \
                             g_0_xxzzz_0_xxxyyzzz_0,  \
                             g_0_xxzzz_0_xxxyyzzz_1,  \
                             g_0_xxzzz_0_xxxyzzz_1,   \
                             g_0_xxzzz_0_xxxyzzzz_0,  \
                             g_0_xxzzz_0_xxxyzzzz_1,  \
                             g_0_xxzzz_0_xxxzzzz_1,   \
                             g_0_xxzzz_0_xxxzzzzz_0,  \
                             g_0_xxzzz_0_xxxzzzzz_1,  \
                             g_0_xxzzz_0_xxyyyyy_1,   \
                             g_0_xxzzz_0_xxyyyyyy_0,  \
                             g_0_xxzzz_0_xxyyyyyy_1,  \
                             g_0_xxzzz_0_xxyyyyyz_0,  \
                             g_0_xxzzz_0_xxyyyyyz_1,  \
                             g_0_xxzzz_0_xxyyyyz_1,   \
                             g_0_xxzzz_0_xxyyyyzz_0,  \
                             g_0_xxzzz_0_xxyyyyzz_1,  \
                             g_0_xxzzz_0_xxyyyzz_1,   \
                             g_0_xxzzz_0_xxyyyzzz_0,  \
                             g_0_xxzzz_0_xxyyyzzz_1,  \
                             g_0_xxzzz_0_xxyyzzz_1,   \
                             g_0_xxzzz_0_xxyyzzzz_0,  \
                             g_0_xxzzz_0_xxyyzzzz_1,  \
                             g_0_xxzzz_0_xxyzzzz_1,   \
                             g_0_xxzzz_0_xxyzzzzz_0,  \
                             g_0_xxzzz_0_xxyzzzzz_1,  \
                             g_0_xxzzz_0_xxzzzzz_1,   \
                             g_0_xxzzz_0_xxzzzzzz_0,  \
                             g_0_xxzzz_0_xxzzzzzz_1,  \
                             g_0_xxzzz_0_xyyyyyy_1,   \
                             g_0_xxzzz_0_xyyyyyyy_0,  \
                             g_0_xxzzz_0_xyyyyyyy_1,  \
                             g_0_xxzzz_0_xyyyyyyz_0,  \
                             g_0_xxzzz_0_xyyyyyyz_1,  \
                             g_0_xxzzz_0_xyyyyyz_1,   \
                             g_0_xxzzz_0_xyyyyyzz_0,  \
                             g_0_xxzzz_0_xyyyyyzz_1,  \
                             g_0_xxzzz_0_xyyyyzz_1,   \
                             g_0_xxzzz_0_xyyyyzzz_0,  \
                             g_0_xxzzz_0_xyyyyzzz_1,  \
                             g_0_xxzzz_0_xyyyzzz_1,   \
                             g_0_xxzzz_0_xyyyzzzz_0,  \
                             g_0_xxzzz_0_xyyyzzzz_1,  \
                             g_0_xxzzz_0_xyyzzzz_1,   \
                             g_0_xxzzz_0_xyyzzzzz_0,  \
                             g_0_xxzzz_0_xyyzzzzz_1,  \
                             g_0_xxzzz_0_xyzzzzz_1,   \
                             g_0_xxzzz_0_xyzzzzzz_0,  \
                             g_0_xxzzz_0_xyzzzzzz_1,  \
                             g_0_xxzzz_0_xzzzzzz_1,   \
                             g_0_xxzzz_0_xzzzzzzz_0,  \
                             g_0_xxzzz_0_xzzzzzzz_1,  \
                             g_0_xxzzz_0_yyyyyyy_1,   \
                             g_0_xxzzz_0_yyyyyyyy_0,  \
                             g_0_xxzzz_0_yyyyyyyy_1,  \
                             g_0_xxzzz_0_yyyyyyyz_0,  \
                             g_0_xxzzz_0_yyyyyyyz_1,  \
                             g_0_xxzzz_0_yyyyyyz_1,   \
                             g_0_xxzzz_0_yyyyyyzz_0,  \
                             g_0_xxzzz_0_yyyyyyzz_1,  \
                             g_0_xxzzz_0_yyyyyzz_1,   \
                             g_0_xxzzz_0_yyyyyzzz_0,  \
                             g_0_xxzzz_0_yyyyyzzz_1,  \
                             g_0_xxzzz_0_yyyyzzz_1,   \
                             g_0_xxzzz_0_yyyyzzzz_0,  \
                             g_0_xxzzz_0_yyyyzzzz_1,  \
                             g_0_xxzzz_0_yyyzzzz_1,   \
                             g_0_xxzzz_0_yyyzzzzz_0,  \
                             g_0_xxzzz_0_yyyzzzzz_1,  \
                             g_0_xxzzz_0_yyzzzzz_1,   \
                             g_0_xxzzz_0_yyzzzzzz_0,  \
                             g_0_xxzzz_0_yyzzzzzz_1,  \
                             g_0_xxzzz_0_yzzzzzz_1,   \
                             g_0_xxzzz_0_yzzzzzzz_0,  \
                             g_0_xxzzz_0_yzzzzzzz_1,  \
                             g_0_xxzzz_0_zzzzzzz_1,   \
                             g_0_xxzzz_0_zzzzzzzz_0,  \
                             g_0_xxzzz_0_zzzzzzzz_1,  \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzz_0_xxxxxxxx_0[i] = g_0_xxzzz_0_xxxxxxxx_0[i] * pb_y + g_0_xxzzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxxxy_0[i] = g_0_xxzzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxxxy_0[i] * pb_y + g_0_xxzzz_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxxxz_0[i] = g_0_xxzzz_0_xxxxxxxz_0[i] * pb_y + g_0_xxzzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxxyy_0[i] =
            2.0 * g_0_xxzzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxxyy_0[i] * pb_y + g_0_xxzzz_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxxyz_0[i] = g_0_xxzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxxyz_0[i] * pb_y + g_0_xxzzz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxxzz_0[i] = g_0_xxzzz_0_xxxxxxzz_0[i] * pb_y + g_0_xxzzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxyyy_0[i] =
            3.0 * g_0_xxzzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxyyy_0[i] * pb_y + g_0_xxzzz_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxyyz_0[i] =
            2.0 * g_0_xxzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxyyz_0[i] * pb_y + g_0_xxzzz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxyzz_0[i] = g_0_xxzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxyzz_0[i] * pb_y + g_0_xxzzz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxzzz_0[i] = g_0_xxzzz_0_xxxxxzzz_0[i] * pb_y + g_0_xxzzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxyyyy_0[i] =
            4.0 * g_0_xxzzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyyyy_0[i] * pb_y + g_0_xxzzz_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxyyyz_0[i] =
            3.0 * g_0_xxzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyyyz_0[i] * pb_y + g_0_xxzzz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxyyzz_0[i] =
            2.0 * g_0_xxzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyyzz_0[i] * pb_y + g_0_xxzzz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxyzzz_0[i] = g_0_xxzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyzzz_0[i] * pb_y + g_0_xxzzz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxzzzz_0[i] = g_0_xxzzz_0_xxxxzzzz_0[i] * pb_y + g_0_xxzzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyyyyy_0[i] =
            5.0 * g_0_xxzzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyyyy_0[i] * pb_y + g_0_xxzzz_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyyyyz_0[i] =
            4.0 * g_0_xxzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyyyz_0[i] * pb_y + g_0_xxzzz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyyyzz_0[i] =
            3.0 * g_0_xxzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyyzz_0[i] * pb_y + g_0_xxzzz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyyzzz_0[i] =
            2.0 * g_0_xxzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyzzz_0[i] * pb_y + g_0_xxzzz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyzzzz_0[i] = g_0_xxzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyzzzz_0[i] * pb_y + g_0_xxzzz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxzzzzz_0[i] = g_0_xxzzz_0_xxxzzzzz_0[i] * pb_y + g_0_xxzzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyyyyy_0[i] =
            6.0 * g_0_xxzzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyyyy_0[i] * pb_y + g_0_xxzzz_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyyyyz_0[i] =
            5.0 * g_0_xxzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyyyz_0[i] * pb_y + g_0_xxzzz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyyyzz_0[i] =
            4.0 * g_0_xxzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyyzz_0[i] * pb_y + g_0_xxzzz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyyzzz_0[i] =
            3.0 * g_0_xxzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyzzz_0[i] * pb_y + g_0_xxzzz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyzzzz_0[i] =
            2.0 * g_0_xxzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyzzzz_0[i] * pb_y + g_0_xxzzz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyzzzzz_0[i] = g_0_xxzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyzzzzz_0[i] * pb_y + g_0_xxzzz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxzzzzzz_0[i] = g_0_xxzzz_0_xxzzzzzz_0[i] * pb_y + g_0_xxzzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyyyyy_0[i] =
            7.0 * g_0_xxzzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyyyy_0[i] * pb_y + g_0_xxzzz_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyyyyz_0[i] =
            6.0 * g_0_xxzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyyyz_0[i] * pb_y + g_0_xxzzz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyyyzz_0[i] =
            5.0 * g_0_xxzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyyzz_0[i] * pb_y + g_0_xxzzz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyyzzz_0[i] =
            4.0 * g_0_xxzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyzzz_0[i] * pb_y + g_0_xxzzz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyzzzz_0[i] =
            3.0 * g_0_xxzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyzzzz_0[i] * pb_y + g_0_xxzzz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyzzzzz_0[i] =
            2.0 * g_0_xxzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyzzzzz_0[i] * pb_y + g_0_xxzzz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyzzzzzz_0[i] = g_0_xxzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyzzzzzz_0[i] * pb_y + g_0_xxzzz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xzzzzzzz_0[i] = g_0_xxzzz_0_xzzzzzzz_0[i] * pb_y + g_0_xxzzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyyyyy_0[i] =
            8.0 * g_0_xxzzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyyyyy_0[i] * pb_y + g_0_xxzzz_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyyyyz_0[i] =
            7.0 * g_0_xxzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyyyyz_0[i] * pb_y + g_0_xxzzz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyyyzz_0[i] =
            6.0 * g_0_xxzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyyyzz_0[i] * pb_y + g_0_xxzzz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyyzzz_0[i] =
            5.0 * g_0_xxzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyyzzz_0[i] * pb_y + g_0_xxzzz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyzzzz_0[i] =
            4.0 * g_0_xxzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyzzzz_0[i] * pb_y + g_0_xxzzz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyzzzzz_0[i] =
            3.0 * g_0_xxzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyzzzzz_0[i] * pb_y + g_0_xxzzz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyzzzzzz_0[i] =
            2.0 * g_0_xxzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyzzzzzz_0[i] * pb_y + g_0_xxzzz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yzzzzzzz_0[i] = g_0_xxzzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yzzzzzzz_0[i] * pb_y + g_0_xxzzz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_zzzzzzzz_0[i] = g_0_xxzzz_0_zzzzzzzz_0[i] * pb_y + g_0_xxzzz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 630-675 components of targeted buffer : SISL

    auto g_0_xxzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 630);

    auto g_0_xxzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 631);

    auto g_0_xxzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 632);

    auto g_0_xxzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 633);

    auto g_0_xxzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 634);

    auto g_0_xxzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 635);

    auto g_0_xxzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 636);

    auto g_0_xxzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 637);

    auto g_0_xxzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 638);

    auto g_0_xxzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 639);

    auto g_0_xxzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 640);

    auto g_0_xxzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 641);

    auto g_0_xxzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 642);

    auto g_0_xxzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 643);

    auto g_0_xxzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 644);

    auto g_0_xxzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 645);

    auto g_0_xxzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 646);

    auto g_0_xxzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 647);

    auto g_0_xxzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 648);

    auto g_0_xxzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 649);

    auto g_0_xxzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 650);

    auto g_0_xxzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 651);

    auto g_0_xxzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 652);

    auto g_0_xxzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 653);

    auto g_0_xxzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 654);

    auto g_0_xxzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 655);

    auto g_0_xxzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 656);

    auto g_0_xxzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 657);

    auto g_0_xxzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 658);

    auto g_0_xxzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 659);

    auto g_0_xxzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 660);

    auto g_0_xxzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 661);

    auto g_0_xxzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 662);

    auto g_0_xxzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 663);

    auto g_0_xxzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 664);

    auto g_0_xxzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 665);

    auto g_0_xxzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 666);

    auto g_0_xxzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 667);

    auto g_0_xxzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 668);

    auto g_0_xxzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 669);

    auto g_0_xxzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 670);

    auto g_0_xxzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 671);

    auto g_0_xxzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 672);

    auto g_0_xxzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 673);

    auto g_0_xxzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 674);

#pragma omp simd aligned(g_0_xxzz_0_xxxxxxxx_0,       \
                             g_0_xxzz_0_xxxxxxxx_1,   \
                             g_0_xxzz_0_xxxxxxxy_0,   \
                             g_0_xxzz_0_xxxxxxxy_1,   \
                             g_0_xxzz_0_xxxxxxyy_0,   \
                             g_0_xxzz_0_xxxxxxyy_1,   \
                             g_0_xxzz_0_xxxxxyyy_0,   \
                             g_0_xxzz_0_xxxxxyyy_1,   \
                             g_0_xxzz_0_xxxxyyyy_0,   \
                             g_0_xxzz_0_xxxxyyyy_1,   \
                             g_0_xxzz_0_xxxyyyyy_0,   \
                             g_0_xxzz_0_xxxyyyyy_1,   \
                             g_0_xxzz_0_xxyyyyyy_0,   \
                             g_0_xxzz_0_xxyyyyyy_1,   \
                             g_0_xxzz_0_xyyyyyyy_0,   \
                             g_0_xxzz_0_xyyyyyyy_1,   \
                             g_0_xxzzz_0_xxxxxxxx_0,  \
                             g_0_xxzzz_0_xxxxxxxx_1,  \
                             g_0_xxzzz_0_xxxxxxxy_0,  \
                             g_0_xxzzz_0_xxxxxxxy_1,  \
                             g_0_xxzzz_0_xxxxxxyy_0,  \
                             g_0_xxzzz_0_xxxxxxyy_1,  \
                             g_0_xxzzz_0_xxxxxyyy_0,  \
                             g_0_xxzzz_0_xxxxxyyy_1,  \
                             g_0_xxzzz_0_xxxxyyyy_0,  \
                             g_0_xxzzz_0_xxxxyyyy_1,  \
                             g_0_xxzzz_0_xxxyyyyy_0,  \
                             g_0_xxzzz_0_xxxyyyyy_1,  \
                             g_0_xxzzz_0_xxyyyyyy_0,  \
                             g_0_xxzzz_0_xxyyyyyy_1,  \
                             g_0_xxzzz_0_xyyyyyyy_0,  \
                             g_0_xxzzz_0_xyyyyyyy_1,  \
                             g_0_xxzzzz_0_xxxxxxxx_0, \
                             g_0_xxzzzz_0_xxxxxxxy_0, \
                             g_0_xxzzzz_0_xxxxxxxz_0, \
                             g_0_xxzzzz_0_xxxxxxyy_0, \
                             g_0_xxzzzz_0_xxxxxxyz_0, \
                             g_0_xxzzzz_0_xxxxxxzz_0, \
                             g_0_xxzzzz_0_xxxxxyyy_0, \
                             g_0_xxzzzz_0_xxxxxyyz_0, \
                             g_0_xxzzzz_0_xxxxxyzz_0, \
                             g_0_xxzzzz_0_xxxxxzzz_0, \
                             g_0_xxzzzz_0_xxxxyyyy_0, \
                             g_0_xxzzzz_0_xxxxyyyz_0, \
                             g_0_xxzzzz_0_xxxxyyzz_0, \
                             g_0_xxzzzz_0_xxxxyzzz_0, \
                             g_0_xxzzzz_0_xxxxzzzz_0, \
                             g_0_xxzzzz_0_xxxyyyyy_0, \
                             g_0_xxzzzz_0_xxxyyyyz_0, \
                             g_0_xxzzzz_0_xxxyyyzz_0, \
                             g_0_xxzzzz_0_xxxyyzzz_0, \
                             g_0_xxzzzz_0_xxxyzzzz_0, \
                             g_0_xxzzzz_0_xxxzzzzz_0, \
                             g_0_xxzzzz_0_xxyyyyyy_0, \
                             g_0_xxzzzz_0_xxyyyyyz_0, \
                             g_0_xxzzzz_0_xxyyyyzz_0, \
                             g_0_xxzzzz_0_xxyyyzzz_0, \
                             g_0_xxzzzz_0_xxyyzzzz_0, \
                             g_0_xxzzzz_0_xxyzzzzz_0, \
                             g_0_xxzzzz_0_xxzzzzzz_0, \
                             g_0_xxzzzz_0_xyyyyyyy_0, \
                             g_0_xxzzzz_0_xyyyyyyz_0, \
                             g_0_xxzzzz_0_xyyyyyzz_0, \
                             g_0_xxzzzz_0_xyyyyzzz_0, \
                             g_0_xxzzzz_0_xyyyzzzz_0, \
                             g_0_xxzzzz_0_xyyzzzzz_0, \
                             g_0_xxzzzz_0_xyzzzzzz_0, \
                             g_0_xxzzzz_0_xzzzzzzz_0, \
                             g_0_xxzzzz_0_yyyyyyyy_0, \
                             g_0_xxzzzz_0_yyyyyyyz_0, \
                             g_0_xxzzzz_0_yyyyyyzz_0, \
                             g_0_xxzzzz_0_yyyyyzzz_0, \
                             g_0_xxzzzz_0_yyyyzzzz_0, \
                             g_0_xxzzzz_0_yyyzzzzz_0, \
                             g_0_xxzzzz_0_yyzzzzzz_0, \
                             g_0_xxzzzz_0_yzzzzzzz_0, \
                             g_0_xxzzzz_0_zzzzzzzz_0, \
                             g_0_xzzzz_0_xxxxxxxz_0,  \
                             g_0_xzzzz_0_xxxxxxxz_1,  \
                             g_0_xzzzz_0_xxxxxxyz_0,  \
                             g_0_xzzzz_0_xxxxxxyz_1,  \
                             g_0_xzzzz_0_xxxxxxz_1,   \
                             g_0_xzzzz_0_xxxxxxzz_0,  \
                             g_0_xzzzz_0_xxxxxxzz_1,  \
                             g_0_xzzzz_0_xxxxxyyz_0,  \
                             g_0_xzzzz_0_xxxxxyyz_1,  \
                             g_0_xzzzz_0_xxxxxyz_1,   \
                             g_0_xzzzz_0_xxxxxyzz_0,  \
                             g_0_xzzzz_0_xxxxxyzz_1,  \
                             g_0_xzzzz_0_xxxxxzz_1,   \
                             g_0_xzzzz_0_xxxxxzzz_0,  \
                             g_0_xzzzz_0_xxxxxzzz_1,  \
                             g_0_xzzzz_0_xxxxyyyz_0,  \
                             g_0_xzzzz_0_xxxxyyyz_1,  \
                             g_0_xzzzz_0_xxxxyyz_1,   \
                             g_0_xzzzz_0_xxxxyyzz_0,  \
                             g_0_xzzzz_0_xxxxyyzz_1,  \
                             g_0_xzzzz_0_xxxxyzz_1,   \
                             g_0_xzzzz_0_xxxxyzzz_0,  \
                             g_0_xzzzz_0_xxxxyzzz_1,  \
                             g_0_xzzzz_0_xxxxzzz_1,   \
                             g_0_xzzzz_0_xxxxzzzz_0,  \
                             g_0_xzzzz_0_xxxxzzzz_1,  \
                             g_0_xzzzz_0_xxxyyyyz_0,  \
                             g_0_xzzzz_0_xxxyyyyz_1,  \
                             g_0_xzzzz_0_xxxyyyz_1,   \
                             g_0_xzzzz_0_xxxyyyzz_0,  \
                             g_0_xzzzz_0_xxxyyyzz_1,  \
                             g_0_xzzzz_0_xxxyyzz_1,   \
                             g_0_xzzzz_0_xxxyyzzz_0,  \
                             g_0_xzzzz_0_xxxyyzzz_1,  \
                             g_0_xzzzz_0_xxxyzzz_1,   \
                             g_0_xzzzz_0_xxxyzzzz_0,  \
                             g_0_xzzzz_0_xxxyzzzz_1,  \
                             g_0_xzzzz_0_xxxzzzz_1,   \
                             g_0_xzzzz_0_xxxzzzzz_0,  \
                             g_0_xzzzz_0_xxxzzzzz_1,  \
                             g_0_xzzzz_0_xxyyyyyz_0,  \
                             g_0_xzzzz_0_xxyyyyyz_1,  \
                             g_0_xzzzz_0_xxyyyyz_1,   \
                             g_0_xzzzz_0_xxyyyyzz_0,  \
                             g_0_xzzzz_0_xxyyyyzz_1,  \
                             g_0_xzzzz_0_xxyyyzz_1,   \
                             g_0_xzzzz_0_xxyyyzzz_0,  \
                             g_0_xzzzz_0_xxyyyzzz_1,  \
                             g_0_xzzzz_0_xxyyzzz_1,   \
                             g_0_xzzzz_0_xxyyzzzz_0,  \
                             g_0_xzzzz_0_xxyyzzzz_1,  \
                             g_0_xzzzz_0_xxyzzzz_1,   \
                             g_0_xzzzz_0_xxyzzzzz_0,  \
                             g_0_xzzzz_0_xxyzzzzz_1,  \
                             g_0_xzzzz_0_xxzzzzz_1,   \
                             g_0_xzzzz_0_xxzzzzzz_0,  \
                             g_0_xzzzz_0_xxzzzzzz_1,  \
                             g_0_xzzzz_0_xyyyyyyz_0,  \
                             g_0_xzzzz_0_xyyyyyyz_1,  \
                             g_0_xzzzz_0_xyyyyyz_1,   \
                             g_0_xzzzz_0_xyyyyyzz_0,  \
                             g_0_xzzzz_0_xyyyyyzz_1,  \
                             g_0_xzzzz_0_xyyyyzz_1,   \
                             g_0_xzzzz_0_xyyyyzzz_0,  \
                             g_0_xzzzz_0_xyyyyzzz_1,  \
                             g_0_xzzzz_0_xyyyzzz_1,   \
                             g_0_xzzzz_0_xyyyzzzz_0,  \
                             g_0_xzzzz_0_xyyyzzzz_1,  \
                             g_0_xzzzz_0_xyyzzzz_1,   \
                             g_0_xzzzz_0_xyyzzzzz_0,  \
                             g_0_xzzzz_0_xyyzzzzz_1,  \
                             g_0_xzzzz_0_xyzzzzz_1,   \
                             g_0_xzzzz_0_xyzzzzzz_0,  \
                             g_0_xzzzz_0_xyzzzzzz_1,  \
                             g_0_xzzzz_0_xzzzzzz_1,   \
                             g_0_xzzzz_0_xzzzzzzz_0,  \
                             g_0_xzzzz_0_xzzzzzzz_1,  \
                             g_0_xzzzz_0_yyyyyyyy_0,  \
                             g_0_xzzzz_0_yyyyyyyy_1,  \
                             g_0_xzzzz_0_yyyyyyyz_0,  \
                             g_0_xzzzz_0_yyyyyyyz_1,  \
                             g_0_xzzzz_0_yyyyyyz_1,   \
                             g_0_xzzzz_0_yyyyyyzz_0,  \
                             g_0_xzzzz_0_yyyyyyzz_1,  \
                             g_0_xzzzz_0_yyyyyzz_1,   \
                             g_0_xzzzz_0_yyyyyzzz_0,  \
                             g_0_xzzzz_0_yyyyyzzz_1,  \
                             g_0_xzzzz_0_yyyyzzz_1,   \
                             g_0_xzzzz_0_yyyyzzzz_0,  \
                             g_0_xzzzz_0_yyyyzzzz_1,  \
                             g_0_xzzzz_0_yyyzzzz_1,   \
                             g_0_xzzzz_0_yyyzzzzz_0,  \
                             g_0_xzzzz_0_yyyzzzzz_1,  \
                             g_0_xzzzz_0_yyzzzzz_1,   \
                             g_0_xzzzz_0_yyzzzzzz_0,  \
                             g_0_xzzzz_0_yyzzzzzz_1,  \
                             g_0_xzzzz_0_yzzzzzz_1,   \
                             g_0_xzzzz_0_yzzzzzzz_0,  \
                             g_0_xzzzz_0_yzzzzzzz_1,  \
                             g_0_xzzzz_0_zzzzzzz_1,   \
                             g_0_xzzzz_0_zzzzzzzz_0,  \
                             g_0_xzzzz_0_zzzzzzzz_1,  \
                             g_0_zzzz_0_xxxxxxxz_0,   \
                             g_0_zzzz_0_xxxxxxxz_1,   \
                             g_0_zzzz_0_xxxxxxyz_0,   \
                             g_0_zzzz_0_xxxxxxyz_1,   \
                             g_0_zzzz_0_xxxxxxzz_0,   \
                             g_0_zzzz_0_xxxxxxzz_1,   \
                             g_0_zzzz_0_xxxxxyyz_0,   \
                             g_0_zzzz_0_xxxxxyyz_1,   \
                             g_0_zzzz_0_xxxxxyzz_0,   \
                             g_0_zzzz_0_xxxxxyzz_1,   \
                             g_0_zzzz_0_xxxxxzzz_0,   \
                             g_0_zzzz_0_xxxxxzzz_1,   \
                             g_0_zzzz_0_xxxxyyyz_0,   \
                             g_0_zzzz_0_xxxxyyyz_1,   \
                             g_0_zzzz_0_xxxxyyzz_0,   \
                             g_0_zzzz_0_xxxxyyzz_1,   \
                             g_0_zzzz_0_xxxxyzzz_0,   \
                             g_0_zzzz_0_xxxxyzzz_1,   \
                             g_0_zzzz_0_xxxxzzzz_0,   \
                             g_0_zzzz_0_xxxxzzzz_1,   \
                             g_0_zzzz_0_xxxyyyyz_0,   \
                             g_0_zzzz_0_xxxyyyyz_1,   \
                             g_0_zzzz_0_xxxyyyzz_0,   \
                             g_0_zzzz_0_xxxyyyzz_1,   \
                             g_0_zzzz_0_xxxyyzzz_0,   \
                             g_0_zzzz_0_xxxyyzzz_1,   \
                             g_0_zzzz_0_xxxyzzzz_0,   \
                             g_0_zzzz_0_xxxyzzzz_1,   \
                             g_0_zzzz_0_xxxzzzzz_0,   \
                             g_0_zzzz_0_xxxzzzzz_1,   \
                             g_0_zzzz_0_xxyyyyyz_0,   \
                             g_0_zzzz_0_xxyyyyyz_1,   \
                             g_0_zzzz_0_xxyyyyzz_0,   \
                             g_0_zzzz_0_xxyyyyzz_1,   \
                             g_0_zzzz_0_xxyyyzzz_0,   \
                             g_0_zzzz_0_xxyyyzzz_1,   \
                             g_0_zzzz_0_xxyyzzzz_0,   \
                             g_0_zzzz_0_xxyyzzzz_1,   \
                             g_0_zzzz_0_xxyzzzzz_0,   \
                             g_0_zzzz_0_xxyzzzzz_1,   \
                             g_0_zzzz_0_xxzzzzzz_0,   \
                             g_0_zzzz_0_xxzzzzzz_1,   \
                             g_0_zzzz_0_xyyyyyyz_0,   \
                             g_0_zzzz_0_xyyyyyyz_1,   \
                             g_0_zzzz_0_xyyyyyzz_0,   \
                             g_0_zzzz_0_xyyyyyzz_1,   \
                             g_0_zzzz_0_xyyyyzzz_0,   \
                             g_0_zzzz_0_xyyyyzzz_1,   \
                             g_0_zzzz_0_xyyyzzzz_0,   \
                             g_0_zzzz_0_xyyyzzzz_1,   \
                             g_0_zzzz_0_xyyzzzzz_0,   \
                             g_0_zzzz_0_xyyzzzzz_1,   \
                             g_0_zzzz_0_xyzzzzzz_0,   \
                             g_0_zzzz_0_xyzzzzzz_1,   \
                             g_0_zzzz_0_xzzzzzzz_0,   \
                             g_0_zzzz_0_xzzzzzzz_1,   \
                             g_0_zzzz_0_yyyyyyyy_0,   \
                             g_0_zzzz_0_yyyyyyyy_1,   \
                             g_0_zzzz_0_yyyyyyyz_0,   \
                             g_0_zzzz_0_yyyyyyyz_1,   \
                             g_0_zzzz_0_yyyyyyzz_0,   \
                             g_0_zzzz_0_yyyyyyzz_1,   \
                             g_0_zzzz_0_yyyyyzzz_0,   \
                             g_0_zzzz_0_yyyyyzzz_1,   \
                             g_0_zzzz_0_yyyyzzzz_0,   \
                             g_0_zzzz_0_yyyyzzzz_1,   \
                             g_0_zzzz_0_yyyzzzzz_0,   \
                             g_0_zzzz_0_yyyzzzzz_1,   \
                             g_0_zzzz_0_yyzzzzzz_0,   \
                             g_0_zzzz_0_yyzzzzzz_1,   \
                             g_0_zzzz_0_yzzzzzzz_0,   \
                             g_0_zzzz_0_yzzzzzzz_1,   \
                             g_0_zzzz_0_zzzzzzzz_0,   \
                             g_0_zzzz_0_zzzzzzzz_1,   \
                             wp_x,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzz_0_xxxxxxxx_0[i] = 3.0 * g_0_xxzz_0_xxxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_xxxxxxxx_0[i] * pb_z + g_0_xxzzz_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxxxxy_0[i] = 3.0 * g_0_xxzz_0_xxxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxxxy_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_xxxxxxxy_0[i] * pb_z + g_0_xxzzz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxxxxz_0[i] = g_0_zzzz_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxxxz_1[i] * fti_ab_0 +
                                     7.0 * g_0_xzzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxxxxz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxxxyy_0[i] = 3.0 * g_0_xxzz_0_xxxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxxyy_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_xxxxxxyy_0[i] * pb_z + g_0_xxzzz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxxxyz_0[i] = g_0_zzzz_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxxyz_1[i] * fti_ab_0 +
                                     6.0 * g_0_xzzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxxxyz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxxxzz_0[i] = g_0_zzzz_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxxzz_1[i] * fti_ab_0 +
                                     6.0 * g_0_xzzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxxxzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxxyyy_0[i] = 3.0 * g_0_xxzz_0_xxxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxyyy_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_xxxxxyyy_0[i] * pb_z + g_0_xxzzz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxxyyz_0[i] = g_0_zzzz_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xzzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxxyyz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxxyzz_0[i] = g_0_zzzz_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxyzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xzzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxxyzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxxzzz_0[i] = g_0_zzzz_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xzzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxxzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxyyyy_0[i] = 3.0 * g_0_xxzz_0_xxxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxyyyy_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_xxxxyyyy_0[i] * pb_z + g_0_xxzzz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxyyyz_0[i] = g_0_zzzz_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xzzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxyyyz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxyyzz_0[i] = g_0_zzzz_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xzzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxyyzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxyzzz_0[i] = g_0_zzzz_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxyzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xzzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxyzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxzzzz_0[i] = g_0_zzzz_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xzzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxzzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxyyyyy_0[i] = 3.0 * g_0_xxzz_0_xxxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_xxxyyyyy_0[i] * pb_z + g_0_xxzzz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxyyyyz_0[i] = g_0_zzzz_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xzzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxyyyyz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxyyyzz_0[i] = g_0_zzzz_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xzzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxyyyzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxyyzzz_0[i] = g_0_zzzz_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xzzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxyyzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxyzzzz_0[i] = g_0_zzzz_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xzzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxyzzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxzzzzz_0[i] = g_0_zzzz_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxzzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xzzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxzzzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyyyyyy_0[i] = 3.0 * g_0_xxzz_0_xxyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_xxyyyyyy_0[i] * pb_z + g_0_xxzzz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxyyyyyz_0[i] = g_0_zzzz_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xzzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxyyyyyz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyyyyzz_0[i] = g_0_zzzz_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xzzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxyyyyzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyyyzzz_0[i] = g_0_zzzz_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xzzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxyyyzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyyzzzz_0[i] = g_0_zzzz_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xzzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxyyzzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyzzzzz_0[i] = g_0_zzzz_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xzzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxyzzzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxzzzzzz_0[i] = g_0_zzzz_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxzzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xzzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxzzzzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyyyyyy_0[i] = 3.0 * g_0_xxzz_0_xyyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxzzz_0_xyyyyyyy_0[i] * pb_z + g_0_xxzzz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xyyyyyyz_0[i] = g_0_zzzz_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyyyz_1[i] * fi_abcd_0 +
                                     g_0_xzzzz_0_xyyyyyyz_0[i] * pb_x + g_0_xzzzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyyyyzz_0[i] = g_0_zzzz_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyyyzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyyzz_1[i] * fi_abcd_0 +
                                     g_0_xzzzz_0_xyyyyyzz_0[i] * pb_x + g_0_xzzzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyyyzzz_0[i] = g_0_zzzz_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyyzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyzzz_1[i] * fi_abcd_0 +
                                     g_0_xzzzz_0_xyyyyzzz_0[i] * pb_x + g_0_xzzzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyyzzzz_0[i] = g_0_zzzz_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyzzzz_1[i] * fi_abcd_0 +
                                     g_0_xzzzz_0_xyyyzzzz_0[i] * pb_x + g_0_xzzzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyzzzzz_0[i] = g_0_zzzz_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyzzzzz_1[i] * fi_abcd_0 +
                                     g_0_xzzzz_0_xyyzzzzz_0[i] * pb_x + g_0_xzzzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyzzzzzz_0[i] = g_0_zzzz_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yzzzzzz_1[i] * fi_abcd_0 +
                                     g_0_xzzzz_0_xyzzzzzz_0[i] * pb_x + g_0_xzzzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xzzzzzzz_0[i] = g_0_zzzz_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_zzzzzzz_1[i] * fi_abcd_0 +
                                     g_0_xzzzz_0_xzzzzzzz_0[i] * pb_x + g_0_xzzzz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyyyyy_0[i] = g_0_zzzz_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyyyyy_0[i] * pb_x +
                                     g_0_xzzzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyyyyz_0[i] = g_0_zzzz_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyyyyz_0[i] * pb_x +
                                     g_0_xzzzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyyyzz_0[i] = g_0_zzzz_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyyyzz_0[i] * pb_x +
                                     g_0_xzzzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyyzzz_0[i] = g_0_zzzz_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyyzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyzzzz_0[i] = g_0_zzzz_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyzzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyzzzzz_0[i] = g_0_zzzz_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyzzzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyzzzzzz_0[i] = g_0_zzzz_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyzzzzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yzzzzzzz_0[i] = g_0_zzzz_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yzzzzzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_zzzzzzzz_0[i] = g_0_zzzz_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_zzzzzzzz_0[i] * pb_x +
                                     g_0_xzzzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 675-720 components of targeted buffer : SISL

    auto g_0_xyyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 675);

    auto g_0_xyyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 676);

    auto g_0_xyyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 677);

    auto g_0_xyyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 678);

    auto g_0_xyyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 679);

    auto g_0_xyyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 680);

    auto g_0_xyyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 681);

    auto g_0_xyyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 682);

    auto g_0_xyyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 683);

    auto g_0_xyyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 684);

    auto g_0_xyyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 685);

    auto g_0_xyyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 686);

    auto g_0_xyyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 687);

    auto g_0_xyyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 688);

    auto g_0_xyyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 689);

    auto g_0_xyyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 690);

    auto g_0_xyyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 691);

    auto g_0_xyyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 692);

    auto g_0_xyyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 693);

    auto g_0_xyyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 694);

    auto g_0_xyyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 695);

    auto g_0_xyyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 696);

    auto g_0_xyyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 697);

    auto g_0_xyyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 698);

    auto g_0_xyyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 699);

    auto g_0_xyyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 700);

    auto g_0_xyyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 701);

    auto g_0_xyyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 702);

    auto g_0_xyyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 703);

    auto g_0_xyyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 704);

    auto g_0_xyyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 705);

    auto g_0_xyyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 706);

    auto g_0_xyyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 707);

    auto g_0_xyyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 708);

    auto g_0_xyyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 709);

    auto g_0_xyyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 710);

    auto g_0_xyyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 711);

    auto g_0_xyyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 712);

    auto g_0_xyyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 713);

    auto g_0_xyyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 714);

    auto g_0_xyyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 715);

    auto g_0_xyyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 716);

    auto g_0_xyyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 717);

    auto g_0_xyyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 718);

    auto g_0_xyyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 719);

#pragma omp simd aligned(g_0_xyyyyy_0_xxxxxxxx_0,     \
                             g_0_xyyyyy_0_xxxxxxxy_0, \
                             g_0_xyyyyy_0_xxxxxxxz_0, \
                             g_0_xyyyyy_0_xxxxxxyy_0, \
                             g_0_xyyyyy_0_xxxxxxyz_0, \
                             g_0_xyyyyy_0_xxxxxxzz_0, \
                             g_0_xyyyyy_0_xxxxxyyy_0, \
                             g_0_xyyyyy_0_xxxxxyyz_0, \
                             g_0_xyyyyy_0_xxxxxyzz_0, \
                             g_0_xyyyyy_0_xxxxxzzz_0, \
                             g_0_xyyyyy_0_xxxxyyyy_0, \
                             g_0_xyyyyy_0_xxxxyyyz_0, \
                             g_0_xyyyyy_0_xxxxyyzz_0, \
                             g_0_xyyyyy_0_xxxxyzzz_0, \
                             g_0_xyyyyy_0_xxxxzzzz_0, \
                             g_0_xyyyyy_0_xxxyyyyy_0, \
                             g_0_xyyyyy_0_xxxyyyyz_0, \
                             g_0_xyyyyy_0_xxxyyyzz_0, \
                             g_0_xyyyyy_0_xxxyyzzz_0, \
                             g_0_xyyyyy_0_xxxyzzzz_0, \
                             g_0_xyyyyy_0_xxxzzzzz_0, \
                             g_0_xyyyyy_0_xxyyyyyy_0, \
                             g_0_xyyyyy_0_xxyyyyyz_0, \
                             g_0_xyyyyy_0_xxyyyyzz_0, \
                             g_0_xyyyyy_0_xxyyyzzz_0, \
                             g_0_xyyyyy_0_xxyyzzzz_0, \
                             g_0_xyyyyy_0_xxyzzzzz_0, \
                             g_0_xyyyyy_0_xxzzzzzz_0, \
                             g_0_xyyyyy_0_xyyyyyyy_0, \
                             g_0_xyyyyy_0_xyyyyyyz_0, \
                             g_0_xyyyyy_0_xyyyyyzz_0, \
                             g_0_xyyyyy_0_xyyyyzzz_0, \
                             g_0_xyyyyy_0_xyyyzzzz_0, \
                             g_0_xyyyyy_0_xyyzzzzz_0, \
                             g_0_xyyyyy_0_xyzzzzzz_0, \
                             g_0_xyyyyy_0_xzzzzzzz_0, \
                             g_0_xyyyyy_0_yyyyyyyy_0, \
                             g_0_xyyyyy_0_yyyyyyyz_0, \
                             g_0_xyyyyy_0_yyyyyyzz_0, \
                             g_0_xyyyyy_0_yyyyyzzz_0, \
                             g_0_xyyyyy_0_yyyyzzzz_0, \
                             g_0_xyyyyy_0_yyyzzzzz_0, \
                             g_0_xyyyyy_0_yyzzzzzz_0, \
                             g_0_xyyyyy_0_yzzzzzzz_0, \
                             g_0_xyyyyy_0_zzzzzzzz_0, \
                             g_0_yyyyy_0_xxxxxxx_1,   \
                             g_0_yyyyy_0_xxxxxxxx_0,  \
                             g_0_yyyyy_0_xxxxxxxx_1,  \
                             g_0_yyyyy_0_xxxxxxxy_0,  \
                             g_0_yyyyy_0_xxxxxxxy_1,  \
                             g_0_yyyyy_0_xxxxxxxz_0,  \
                             g_0_yyyyy_0_xxxxxxxz_1,  \
                             g_0_yyyyy_0_xxxxxxy_1,   \
                             g_0_yyyyy_0_xxxxxxyy_0,  \
                             g_0_yyyyy_0_xxxxxxyy_1,  \
                             g_0_yyyyy_0_xxxxxxyz_0,  \
                             g_0_yyyyy_0_xxxxxxyz_1,  \
                             g_0_yyyyy_0_xxxxxxz_1,   \
                             g_0_yyyyy_0_xxxxxxzz_0,  \
                             g_0_yyyyy_0_xxxxxxzz_1,  \
                             g_0_yyyyy_0_xxxxxyy_1,   \
                             g_0_yyyyy_0_xxxxxyyy_0,  \
                             g_0_yyyyy_0_xxxxxyyy_1,  \
                             g_0_yyyyy_0_xxxxxyyz_0,  \
                             g_0_yyyyy_0_xxxxxyyz_1,  \
                             g_0_yyyyy_0_xxxxxyz_1,   \
                             g_0_yyyyy_0_xxxxxyzz_0,  \
                             g_0_yyyyy_0_xxxxxyzz_1,  \
                             g_0_yyyyy_0_xxxxxzz_1,   \
                             g_0_yyyyy_0_xxxxxzzz_0,  \
                             g_0_yyyyy_0_xxxxxzzz_1,  \
                             g_0_yyyyy_0_xxxxyyy_1,   \
                             g_0_yyyyy_0_xxxxyyyy_0,  \
                             g_0_yyyyy_0_xxxxyyyy_1,  \
                             g_0_yyyyy_0_xxxxyyyz_0,  \
                             g_0_yyyyy_0_xxxxyyyz_1,  \
                             g_0_yyyyy_0_xxxxyyz_1,   \
                             g_0_yyyyy_0_xxxxyyzz_0,  \
                             g_0_yyyyy_0_xxxxyyzz_1,  \
                             g_0_yyyyy_0_xxxxyzz_1,   \
                             g_0_yyyyy_0_xxxxyzzz_0,  \
                             g_0_yyyyy_0_xxxxyzzz_1,  \
                             g_0_yyyyy_0_xxxxzzz_1,   \
                             g_0_yyyyy_0_xxxxzzzz_0,  \
                             g_0_yyyyy_0_xxxxzzzz_1,  \
                             g_0_yyyyy_0_xxxyyyy_1,   \
                             g_0_yyyyy_0_xxxyyyyy_0,  \
                             g_0_yyyyy_0_xxxyyyyy_1,  \
                             g_0_yyyyy_0_xxxyyyyz_0,  \
                             g_0_yyyyy_0_xxxyyyyz_1,  \
                             g_0_yyyyy_0_xxxyyyz_1,   \
                             g_0_yyyyy_0_xxxyyyzz_0,  \
                             g_0_yyyyy_0_xxxyyyzz_1,  \
                             g_0_yyyyy_0_xxxyyzz_1,   \
                             g_0_yyyyy_0_xxxyyzzz_0,  \
                             g_0_yyyyy_0_xxxyyzzz_1,  \
                             g_0_yyyyy_0_xxxyzzz_1,   \
                             g_0_yyyyy_0_xxxyzzzz_0,  \
                             g_0_yyyyy_0_xxxyzzzz_1,  \
                             g_0_yyyyy_0_xxxzzzz_1,   \
                             g_0_yyyyy_0_xxxzzzzz_0,  \
                             g_0_yyyyy_0_xxxzzzzz_1,  \
                             g_0_yyyyy_0_xxyyyyy_1,   \
                             g_0_yyyyy_0_xxyyyyyy_0,  \
                             g_0_yyyyy_0_xxyyyyyy_1,  \
                             g_0_yyyyy_0_xxyyyyyz_0,  \
                             g_0_yyyyy_0_xxyyyyyz_1,  \
                             g_0_yyyyy_0_xxyyyyz_1,   \
                             g_0_yyyyy_0_xxyyyyzz_0,  \
                             g_0_yyyyy_0_xxyyyyzz_1,  \
                             g_0_yyyyy_0_xxyyyzz_1,   \
                             g_0_yyyyy_0_xxyyyzzz_0,  \
                             g_0_yyyyy_0_xxyyyzzz_1,  \
                             g_0_yyyyy_0_xxyyzzz_1,   \
                             g_0_yyyyy_0_xxyyzzzz_0,  \
                             g_0_yyyyy_0_xxyyzzzz_1,  \
                             g_0_yyyyy_0_xxyzzzz_1,   \
                             g_0_yyyyy_0_xxyzzzzz_0,  \
                             g_0_yyyyy_0_xxyzzzzz_1,  \
                             g_0_yyyyy_0_xxzzzzz_1,   \
                             g_0_yyyyy_0_xxzzzzzz_0,  \
                             g_0_yyyyy_0_xxzzzzzz_1,  \
                             g_0_yyyyy_0_xyyyyyy_1,   \
                             g_0_yyyyy_0_xyyyyyyy_0,  \
                             g_0_yyyyy_0_xyyyyyyy_1,  \
                             g_0_yyyyy_0_xyyyyyyz_0,  \
                             g_0_yyyyy_0_xyyyyyyz_1,  \
                             g_0_yyyyy_0_xyyyyyz_1,   \
                             g_0_yyyyy_0_xyyyyyzz_0,  \
                             g_0_yyyyy_0_xyyyyyzz_1,  \
                             g_0_yyyyy_0_xyyyyzz_1,   \
                             g_0_yyyyy_0_xyyyyzzz_0,  \
                             g_0_yyyyy_0_xyyyyzzz_1,  \
                             g_0_yyyyy_0_xyyyzzz_1,   \
                             g_0_yyyyy_0_xyyyzzzz_0,  \
                             g_0_yyyyy_0_xyyyzzzz_1,  \
                             g_0_yyyyy_0_xyyzzzz_1,   \
                             g_0_yyyyy_0_xyyzzzzz_0,  \
                             g_0_yyyyy_0_xyyzzzzz_1,  \
                             g_0_yyyyy_0_xyzzzzz_1,   \
                             g_0_yyyyy_0_xyzzzzzz_0,  \
                             g_0_yyyyy_0_xyzzzzzz_1,  \
                             g_0_yyyyy_0_xzzzzzz_1,   \
                             g_0_yyyyy_0_xzzzzzzz_0,  \
                             g_0_yyyyy_0_xzzzzzzz_1,  \
                             g_0_yyyyy_0_yyyyyyy_1,   \
                             g_0_yyyyy_0_yyyyyyyy_0,  \
                             g_0_yyyyy_0_yyyyyyyy_1,  \
                             g_0_yyyyy_0_yyyyyyyz_0,  \
                             g_0_yyyyy_0_yyyyyyyz_1,  \
                             g_0_yyyyy_0_yyyyyyz_1,   \
                             g_0_yyyyy_0_yyyyyyzz_0,  \
                             g_0_yyyyy_0_yyyyyyzz_1,  \
                             g_0_yyyyy_0_yyyyyzz_1,   \
                             g_0_yyyyy_0_yyyyyzzz_0,  \
                             g_0_yyyyy_0_yyyyyzzz_1,  \
                             g_0_yyyyy_0_yyyyzzz_1,   \
                             g_0_yyyyy_0_yyyyzzzz_0,  \
                             g_0_yyyyy_0_yyyyzzzz_1,  \
                             g_0_yyyyy_0_yyyzzzz_1,   \
                             g_0_yyyyy_0_yyyzzzzz_0,  \
                             g_0_yyyyy_0_yyyzzzzz_1,  \
                             g_0_yyyyy_0_yyzzzzz_1,   \
                             g_0_yyyyy_0_yyzzzzzz_0,  \
                             g_0_yyyyy_0_yyzzzzzz_1,  \
                             g_0_yyyyy_0_yzzzzzz_1,   \
                             g_0_yyyyy_0_yzzzzzzz_0,  \
                             g_0_yyyyy_0_yzzzzzzz_1,  \
                             g_0_yyyyy_0_zzzzzzz_1,   \
                             g_0_yyyyy_0_zzzzzzzz_0,  \
                             g_0_yyyyy_0_zzzzzzzz_1,  \
                             wp_x,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyy_0_xxxxxxxx_0[i] =
            8.0 * g_0_yyyyy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxxx_0[i] * pb_x + g_0_yyyyy_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxxxy_0[i] =
            7.0 * g_0_yyyyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxxy_0[i] * pb_x + g_0_yyyyy_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxxxz_0[i] =
            7.0 * g_0_yyyyy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxxz_0[i] * pb_x + g_0_yyyyy_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxxyy_0[i] =
            6.0 * g_0_yyyyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxyy_0[i] * pb_x + g_0_yyyyy_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxxyz_0[i] =
            6.0 * g_0_yyyyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxyz_0[i] * pb_x + g_0_yyyyy_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxxzz_0[i] =
            6.0 * g_0_yyyyy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxzz_0[i] * pb_x + g_0_yyyyy_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxyyy_0[i] =
            5.0 * g_0_yyyyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyyy_0[i] * pb_x + g_0_yyyyy_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxyyz_0[i] =
            5.0 * g_0_yyyyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyyz_0[i] * pb_x + g_0_yyyyy_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxyzz_0[i] =
            5.0 * g_0_yyyyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyzz_0[i] * pb_x + g_0_yyyyy_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxzzz_0[i] =
            5.0 * g_0_yyyyy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxzzz_0[i] * pb_x + g_0_yyyyy_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxyyyy_0[i] =
            4.0 * g_0_yyyyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyyy_0[i] * pb_x + g_0_yyyyy_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxyyyz_0[i] =
            4.0 * g_0_yyyyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyyz_0[i] * pb_x + g_0_yyyyy_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxyyzz_0[i] =
            4.0 * g_0_yyyyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyzz_0[i] * pb_x + g_0_yyyyy_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxyzzz_0[i] =
            4.0 * g_0_yyyyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyzzz_0[i] * pb_x + g_0_yyyyy_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxzzzz_0[i] =
            4.0 * g_0_yyyyy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxzzzz_0[i] * pb_x + g_0_yyyyy_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyyyyy_0[i] =
            3.0 * g_0_yyyyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyyy_0[i] * pb_x + g_0_yyyyy_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyyyyz_0[i] =
            3.0 * g_0_yyyyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyyz_0[i] * pb_x + g_0_yyyyy_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyyyzz_0[i] =
            3.0 * g_0_yyyyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyzz_0[i] * pb_x + g_0_yyyyy_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyyzzz_0[i] =
            3.0 * g_0_yyyyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyzzz_0[i] * pb_x + g_0_yyyyy_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyzzzz_0[i] =
            3.0 * g_0_yyyyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyzzzz_0[i] * pb_x + g_0_yyyyy_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxzzzzz_0[i] =
            3.0 * g_0_yyyyy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxzzzzz_0[i] * pb_x + g_0_yyyyy_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyyyyy_0[i] =
            2.0 * g_0_yyyyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyyy_0[i] * pb_x + g_0_yyyyy_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyyyyz_0[i] =
            2.0 * g_0_yyyyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyyz_0[i] * pb_x + g_0_yyyyy_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyyyzz_0[i] =
            2.0 * g_0_yyyyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyzz_0[i] * pb_x + g_0_yyyyy_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyyzzz_0[i] =
            2.0 * g_0_yyyyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyzzz_0[i] * pb_x + g_0_yyyyy_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyzzzz_0[i] =
            2.0 * g_0_yyyyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyzzzz_0[i] * pb_x + g_0_yyyyy_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyzzzzz_0[i] =
            2.0 * g_0_yyyyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyzzzzz_0[i] * pb_x + g_0_yyyyy_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxzzzzzz_0[i] =
            2.0 * g_0_yyyyy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxzzzzzz_0[i] * pb_x + g_0_yyyyy_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyyyyy_0[i] = g_0_yyyyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyyy_0[i] * pb_x + g_0_yyyyy_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyyyyz_0[i] = g_0_yyyyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyyz_0[i] * pb_x + g_0_yyyyy_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyyyzz_0[i] = g_0_yyyyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyzz_0[i] * pb_x + g_0_yyyyy_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyyzzz_0[i] = g_0_yyyyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyzzz_0[i] * pb_x + g_0_yyyyy_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyzzzz_0[i] = g_0_yyyyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyzzzz_0[i] * pb_x + g_0_yyyyy_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyzzzzz_0[i] = g_0_yyyyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyzzzzz_0[i] * pb_x + g_0_yyyyy_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyzzzzzz_0[i] = g_0_yyyyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzzzzzz_0[i] * pb_x + g_0_yyyyy_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xzzzzzzz_0[i] = g_0_yyyyy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xzzzzzzz_0[i] * pb_x + g_0_yyyyy_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyyyyy_0[i] = g_0_yyyyy_0_yyyyyyyy_0[i] * pb_x + g_0_yyyyy_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyyyyz_0[i] = g_0_yyyyy_0_yyyyyyyz_0[i] * pb_x + g_0_yyyyy_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyyyzz_0[i] = g_0_yyyyy_0_yyyyyyzz_0[i] * pb_x + g_0_yyyyy_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyyzzz_0[i] = g_0_yyyyy_0_yyyyyzzz_0[i] * pb_x + g_0_yyyyy_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyzzzz_0[i] = g_0_yyyyy_0_yyyyzzzz_0[i] * pb_x + g_0_yyyyy_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyzzzzz_0[i] = g_0_yyyyy_0_yyyzzzzz_0[i] * pb_x + g_0_yyyyy_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyzzzzzz_0[i] = g_0_yyyyy_0_yyzzzzzz_0[i] * pb_x + g_0_yyyyy_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yzzzzzzz_0[i] = g_0_yyyyy_0_yzzzzzzz_0[i] * pb_x + g_0_yyyyy_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_zzzzzzzz_0[i] = g_0_yyyyy_0_zzzzzzzz_0[i] * pb_x + g_0_yyyyy_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 720-765 components of targeted buffer : SISL

    auto g_0_xyyyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 720);

    auto g_0_xyyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 721);

    auto g_0_xyyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 722);

    auto g_0_xyyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 723);

    auto g_0_xyyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 724);

    auto g_0_xyyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 725);

    auto g_0_xyyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 726);

    auto g_0_xyyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 727);

    auto g_0_xyyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 728);

    auto g_0_xyyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 729);

    auto g_0_xyyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 730);

    auto g_0_xyyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 731);

    auto g_0_xyyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 732);

    auto g_0_xyyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 733);

    auto g_0_xyyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 734);

    auto g_0_xyyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 735);

    auto g_0_xyyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 736);

    auto g_0_xyyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 737);

    auto g_0_xyyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 738);

    auto g_0_xyyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 739);

    auto g_0_xyyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 740);

    auto g_0_xyyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 741);

    auto g_0_xyyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 742);

    auto g_0_xyyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 743);

    auto g_0_xyyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 744);

    auto g_0_xyyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 745);

    auto g_0_xyyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 746);

    auto g_0_xyyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 747);

    auto g_0_xyyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 748);

    auto g_0_xyyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 749);

    auto g_0_xyyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 750);

    auto g_0_xyyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 751);

    auto g_0_xyyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 752);

    auto g_0_xyyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 753);

    auto g_0_xyyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 754);

    auto g_0_xyyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 755);

    auto g_0_xyyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 756);

    auto g_0_xyyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 757);

    auto g_0_xyyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 758);

    auto g_0_xyyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 759);

    auto g_0_xyyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 760);

    auto g_0_xyyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 761);

    auto g_0_xyyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 762);

    auto g_0_xyyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 763);

    auto g_0_xyyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 764);

#pragma omp simd aligned(g_0_xyyyy_0_xxxxxxxx_0,      \
                             g_0_xyyyy_0_xxxxxxxx_1,  \
                             g_0_xyyyy_0_xxxxxxxy_0,  \
                             g_0_xyyyy_0_xxxxxxxy_1,  \
                             g_0_xyyyy_0_xxxxxxyy_0,  \
                             g_0_xyyyy_0_xxxxxxyy_1,  \
                             g_0_xyyyy_0_xxxxxyyy_0,  \
                             g_0_xyyyy_0_xxxxxyyy_1,  \
                             g_0_xyyyy_0_xxxxyyyy_0,  \
                             g_0_xyyyy_0_xxxxyyyy_1,  \
                             g_0_xyyyy_0_xxxyyyyy_0,  \
                             g_0_xyyyy_0_xxxyyyyy_1,  \
                             g_0_xyyyy_0_xxyyyyyy_0,  \
                             g_0_xyyyy_0_xxyyyyyy_1,  \
                             g_0_xyyyy_0_xyyyyyyy_0,  \
                             g_0_xyyyy_0_xyyyyyyy_1,  \
                             g_0_xyyyyz_0_xxxxxxxx_0, \
                             g_0_xyyyyz_0_xxxxxxxy_0, \
                             g_0_xyyyyz_0_xxxxxxxz_0, \
                             g_0_xyyyyz_0_xxxxxxyy_0, \
                             g_0_xyyyyz_0_xxxxxxyz_0, \
                             g_0_xyyyyz_0_xxxxxxzz_0, \
                             g_0_xyyyyz_0_xxxxxyyy_0, \
                             g_0_xyyyyz_0_xxxxxyyz_0, \
                             g_0_xyyyyz_0_xxxxxyzz_0, \
                             g_0_xyyyyz_0_xxxxxzzz_0, \
                             g_0_xyyyyz_0_xxxxyyyy_0, \
                             g_0_xyyyyz_0_xxxxyyyz_0, \
                             g_0_xyyyyz_0_xxxxyyzz_0, \
                             g_0_xyyyyz_0_xxxxyzzz_0, \
                             g_0_xyyyyz_0_xxxxzzzz_0, \
                             g_0_xyyyyz_0_xxxyyyyy_0, \
                             g_0_xyyyyz_0_xxxyyyyz_0, \
                             g_0_xyyyyz_0_xxxyyyzz_0, \
                             g_0_xyyyyz_0_xxxyyzzz_0, \
                             g_0_xyyyyz_0_xxxyzzzz_0, \
                             g_0_xyyyyz_0_xxxzzzzz_0, \
                             g_0_xyyyyz_0_xxyyyyyy_0, \
                             g_0_xyyyyz_0_xxyyyyyz_0, \
                             g_0_xyyyyz_0_xxyyyyzz_0, \
                             g_0_xyyyyz_0_xxyyyzzz_0, \
                             g_0_xyyyyz_0_xxyyzzzz_0, \
                             g_0_xyyyyz_0_xxyzzzzz_0, \
                             g_0_xyyyyz_0_xxzzzzzz_0, \
                             g_0_xyyyyz_0_xyyyyyyy_0, \
                             g_0_xyyyyz_0_xyyyyyyz_0, \
                             g_0_xyyyyz_0_xyyyyyzz_0, \
                             g_0_xyyyyz_0_xyyyyzzz_0, \
                             g_0_xyyyyz_0_xyyyzzzz_0, \
                             g_0_xyyyyz_0_xyyzzzzz_0, \
                             g_0_xyyyyz_0_xyzzzzzz_0, \
                             g_0_xyyyyz_0_xzzzzzzz_0, \
                             g_0_xyyyyz_0_yyyyyyyy_0, \
                             g_0_xyyyyz_0_yyyyyyyz_0, \
                             g_0_xyyyyz_0_yyyyyyzz_0, \
                             g_0_xyyyyz_0_yyyyyzzz_0, \
                             g_0_xyyyyz_0_yyyyzzzz_0, \
                             g_0_xyyyyz_0_yyyzzzzz_0, \
                             g_0_xyyyyz_0_yyzzzzzz_0, \
                             g_0_xyyyyz_0_yzzzzzzz_0, \
                             g_0_xyyyyz_0_zzzzzzzz_0, \
                             g_0_yyyyz_0_xxxxxxxz_0,  \
                             g_0_yyyyz_0_xxxxxxxz_1,  \
                             g_0_yyyyz_0_xxxxxxyz_0,  \
                             g_0_yyyyz_0_xxxxxxyz_1,  \
                             g_0_yyyyz_0_xxxxxxz_1,   \
                             g_0_yyyyz_0_xxxxxxzz_0,  \
                             g_0_yyyyz_0_xxxxxxzz_1,  \
                             g_0_yyyyz_0_xxxxxyyz_0,  \
                             g_0_yyyyz_0_xxxxxyyz_1,  \
                             g_0_yyyyz_0_xxxxxyz_1,   \
                             g_0_yyyyz_0_xxxxxyzz_0,  \
                             g_0_yyyyz_0_xxxxxyzz_1,  \
                             g_0_yyyyz_0_xxxxxzz_1,   \
                             g_0_yyyyz_0_xxxxxzzz_0,  \
                             g_0_yyyyz_0_xxxxxzzz_1,  \
                             g_0_yyyyz_0_xxxxyyyz_0,  \
                             g_0_yyyyz_0_xxxxyyyz_1,  \
                             g_0_yyyyz_0_xxxxyyz_1,   \
                             g_0_yyyyz_0_xxxxyyzz_0,  \
                             g_0_yyyyz_0_xxxxyyzz_1,  \
                             g_0_yyyyz_0_xxxxyzz_1,   \
                             g_0_yyyyz_0_xxxxyzzz_0,  \
                             g_0_yyyyz_0_xxxxyzzz_1,  \
                             g_0_yyyyz_0_xxxxzzz_1,   \
                             g_0_yyyyz_0_xxxxzzzz_0,  \
                             g_0_yyyyz_0_xxxxzzzz_1,  \
                             g_0_yyyyz_0_xxxyyyyz_0,  \
                             g_0_yyyyz_0_xxxyyyyz_1,  \
                             g_0_yyyyz_0_xxxyyyz_1,   \
                             g_0_yyyyz_0_xxxyyyzz_0,  \
                             g_0_yyyyz_0_xxxyyyzz_1,  \
                             g_0_yyyyz_0_xxxyyzz_1,   \
                             g_0_yyyyz_0_xxxyyzzz_0,  \
                             g_0_yyyyz_0_xxxyyzzz_1,  \
                             g_0_yyyyz_0_xxxyzzz_1,   \
                             g_0_yyyyz_0_xxxyzzzz_0,  \
                             g_0_yyyyz_0_xxxyzzzz_1,  \
                             g_0_yyyyz_0_xxxzzzz_1,   \
                             g_0_yyyyz_0_xxxzzzzz_0,  \
                             g_0_yyyyz_0_xxxzzzzz_1,  \
                             g_0_yyyyz_0_xxyyyyyz_0,  \
                             g_0_yyyyz_0_xxyyyyyz_1,  \
                             g_0_yyyyz_0_xxyyyyz_1,   \
                             g_0_yyyyz_0_xxyyyyzz_0,  \
                             g_0_yyyyz_0_xxyyyyzz_1,  \
                             g_0_yyyyz_0_xxyyyzz_1,   \
                             g_0_yyyyz_0_xxyyyzzz_0,  \
                             g_0_yyyyz_0_xxyyyzzz_1,  \
                             g_0_yyyyz_0_xxyyzzz_1,   \
                             g_0_yyyyz_0_xxyyzzzz_0,  \
                             g_0_yyyyz_0_xxyyzzzz_1,  \
                             g_0_yyyyz_0_xxyzzzz_1,   \
                             g_0_yyyyz_0_xxyzzzzz_0,  \
                             g_0_yyyyz_0_xxyzzzzz_1,  \
                             g_0_yyyyz_0_xxzzzzz_1,   \
                             g_0_yyyyz_0_xxzzzzzz_0,  \
                             g_0_yyyyz_0_xxzzzzzz_1,  \
                             g_0_yyyyz_0_xyyyyyyz_0,  \
                             g_0_yyyyz_0_xyyyyyyz_1,  \
                             g_0_yyyyz_0_xyyyyyz_1,   \
                             g_0_yyyyz_0_xyyyyyzz_0,  \
                             g_0_yyyyz_0_xyyyyyzz_1,  \
                             g_0_yyyyz_0_xyyyyzz_1,   \
                             g_0_yyyyz_0_xyyyyzzz_0,  \
                             g_0_yyyyz_0_xyyyyzzz_1,  \
                             g_0_yyyyz_0_xyyyzzz_1,   \
                             g_0_yyyyz_0_xyyyzzzz_0,  \
                             g_0_yyyyz_0_xyyyzzzz_1,  \
                             g_0_yyyyz_0_xyyzzzz_1,   \
                             g_0_yyyyz_0_xyyzzzzz_0,  \
                             g_0_yyyyz_0_xyyzzzzz_1,  \
                             g_0_yyyyz_0_xyzzzzz_1,   \
                             g_0_yyyyz_0_xyzzzzzz_0,  \
                             g_0_yyyyz_0_xyzzzzzz_1,  \
                             g_0_yyyyz_0_xzzzzzz_1,   \
                             g_0_yyyyz_0_xzzzzzzz_0,  \
                             g_0_yyyyz_0_xzzzzzzz_1,  \
                             g_0_yyyyz_0_yyyyyyyy_0,  \
                             g_0_yyyyz_0_yyyyyyyy_1,  \
                             g_0_yyyyz_0_yyyyyyyz_0,  \
                             g_0_yyyyz_0_yyyyyyyz_1,  \
                             g_0_yyyyz_0_yyyyyyz_1,   \
                             g_0_yyyyz_0_yyyyyyzz_0,  \
                             g_0_yyyyz_0_yyyyyyzz_1,  \
                             g_0_yyyyz_0_yyyyyzz_1,   \
                             g_0_yyyyz_0_yyyyyzzz_0,  \
                             g_0_yyyyz_0_yyyyyzzz_1,  \
                             g_0_yyyyz_0_yyyyzzz_1,   \
                             g_0_yyyyz_0_yyyyzzzz_0,  \
                             g_0_yyyyz_0_yyyyzzzz_1,  \
                             g_0_yyyyz_0_yyyzzzz_1,   \
                             g_0_yyyyz_0_yyyzzzzz_0,  \
                             g_0_yyyyz_0_yyyzzzzz_1,  \
                             g_0_yyyyz_0_yyzzzzz_1,   \
                             g_0_yyyyz_0_yyzzzzzz_0,  \
                             g_0_yyyyz_0_yyzzzzzz_1,  \
                             g_0_yyyyz_0_yzzzzzz_1,   \
                             g_0_yyyyz_0_yzzzzzzz_0,  \
                             g_0_yyyyz_0_yzzzzzzz_1,  \
                             g_0_yyyyz_0_zzzzzzz_1,   \
                             g_0_yyyyz_0_zzzzzzzz_0,  \
                             g_0_yyyyz_0_zzzzzzzz_1,  \
                             wp_x,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyz_0_xxxxxxxx_0[i] = g_0_xyyyy_0_xxxxxxxx_0[i] * pb_z + g_0_xyyyy_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxxxxy_0[i] = g_0_xyyyy_0_xxxxxxxy_0[i] * pb_z + g_0_xyyyy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxxxxz_0[i] =
            7.0 * g_0_yyyyz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxxxxz_0[i] * pb_x + g_0_yyyyz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxxxyy_0[i] = g_0_xyyyy_0_xxxxxxyy_0[i] * pb_z + g_0_xyyyy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxxxyz_0[i] =
            6.0 * g_0_yyyyz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxxxyz_0[i] * pb_x + g_0_yyyyz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxxxzz_0[i] =
            6.0 * g_0_yyyyz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxxxzz_0[i] * pb_x + g_0_yyyyz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxxyyy_0[i] = g_0_xyyyy_0_xxxxxyyy_0[i] * pb_z + g_0_xyyyy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxxyyz_0[i] =
            5.0 * g_0_yyyyz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxxyyz_0[i] * pb_x + g_0_yyyyz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxxyzz_0[i] =
            5.0 * g_0_yyyyz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxxyzz_0[i] * pb_x + g_0_yyyyz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxxzzz_0[i] =
            5.0 * g_0_yyyyz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxxzzz_0[i] * pb_x + g_0_yyyyz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxyyyy_0[i] = g_0_xyyyy_0_xxxxyyyy_0[i] * pb_z + g_0_xyyyy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxyyyz_0[i] =
            4.0 * g_0_yyyyz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxyyyz_0[i] * pb_x + g_0_yyyyz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxyyzz_0[i] =
            4.0 * g_0_yyyyz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxyyzz_0[i] * pb_x + g_0_yyyyz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxyzzz_0[i] =
            4.0 * g_0_yyyyz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxyzzz_0[i] * pb_x + g_0_yyyyz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxzzzz_0[i] =
            4.0 * g_0_yyyyz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxzzzz_0[i] * pb_x + g_0_yyyyz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxyyyyy_0[i] = g_0_xyyyy_0_xxxyyyyy_0[i] * pb_z + g_0_xyyyy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxyyyyz_0[i] =
            3.0 * g_0_yyyyz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxyyyyz_0[i] * pb_x + g_0_yyyyz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxyyyzz_0[i] =
            3.0 * g_0_yyyyz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxyyyzz_0[i] * pb_x + g_0_yyyyz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxyyzzz_0[i] =
            3.0 * g_0_yyyyz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxyyzzz_0[i] * pb_x + g_0_yyyyz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxyzzzz_0[i] =
            3.0 * g_0_yyyyz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxyzzzz_0[i] * pb_x + g_0_yyyyz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxzzzzz_0[i] =
            3.0 * g_0_yyyyz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxzzzzz_0[i] * pb_x + g_0_yyyyz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyyyyyy_0[i] = g_0_xyyyy_0_xxyyyyyy_0[i] * pb_z + g_0_xyyyy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxyyyyyz_0[i] =
            2.0 * g_0_yyyyz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyyyyyz_0[i] * pb_x + g_0_yyyyz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyyyyzz_0[i] =
            2.0 * g_0_yyyyz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyyyyzz_0[i] * pb_x + g_0_yyyyz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyyyzzz_0[i] =
            2.0 * g_0_yyyyz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyyyzzz_0[i] * pb_x + g_0_yyyyz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyyzzzz_0[i] =
            2.0 * g_0_yyyyz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyyzzzz_0[i] * pb_x + g_0_yyyyz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyzzzzz_0[i] =
            2.0 * g_0_yyyyz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyzzzzz_0[i] * pb_x + g_0_yyyyz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxzzzzzz_0[i] =
            2.0 * g_0_yyyyz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxzzzzzz_0[i] * pb_x + g_0_yyyyz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyyyyyy_0[i] = g_0_xyyyy_0_xyyyyyyy_0[i] * pb_z + g_0_xyyyy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xyyyyyyz_0[i] = g_0_yyyyz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyyyyyz_0[i] * pb_x + g_0_yyyyz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyyyyzz_0[i] = g_0_yyyyz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyyyyzz_0[i] * pb_x + g_0_yyyyz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyyyzzz_0[i] = g_0_yyyyz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyyyzzz_0[i] * pb_x + g_0_yyyyz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyyzzzz_0[i] = g_0_yyyyz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyyzzzz_0[i] * pb_x + g_0_yyyyz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyzzzzz_0[i] = g_0_yyyyz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyzzzzz_0[i] * pb_x + g_0_yyyyz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyzzzzzz_0[i] = g_0_yyyyz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyzzzzzz_0[i] * pb_x + g_0_yyyyz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xzzzzzzz_0[i] = g_0_yyyyz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xzzzzzzz_0[i] * pb_x + g_0_yyyyz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyyyyy_0[i] = g_0_yyyyz_0_yyyyyyyy_0[i] * pb_x + g_0_yyyyz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyyyyz_0[i] = g_0_yyyyz_0_yyyyyyyz_0[i] * pb_x + g_0_yyyyz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyyyzz_0[i] = g_0_yyyyz_0_yyyyyyzz_0[i] * pb_x + g_0_yyyyz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyyzzz_0[i] = g_0_yyyyz_0_yyyyyzzz_0[i] * pb_x + g_0_yyyyz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyzzzz_0[i] = g_0_yyyyz_0_yyyyzzzz_0[i] * pb_x + g_0_yyyyz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyzzzzz_0[i] = g_0_yyyyz_0_yyyzzzzz_0[i] * pb_x + g_0_yyyyz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyzzzzzz_0[i] = g_0_yyyyz_0_yyzzzzzz_0[i] * pb_x + g_0_yyyyz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yzzzzzzz_0[i] = g_0_yyyyz_0_yzzzzzzz_0[i] * pb_x + g_0_yyyyz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_zzzzzzzz_0[i] = g_0_yyyyz_0_zzzzzzzz_0[i] * pb_x + g_0_yyyyz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 765-810 components of targeted buffer : SISL

    auto g_0_xyyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 765);

    auto g_0_xyyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 766);

    auto g_0_xyyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 767);

    auto g_0_xyyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 768);

    auto g_0_xyyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 769);

    auto g_0_xyyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 770);

    auto g_0_xyyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 771);

    auto g_0_xyyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 772);

    auto g_0_xyyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 773);

    auto g_0_xyyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 774);

    auto g_0_xyyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 775);

    auto g_0_xyyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 776);

    auto g_0_xyyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 777);

    auto g_0_xyyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 778);

    auto g_0_xyyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 779);

    auto g_0_xyyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 780);

    auto g_0_xyyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 781);

    auto g_0_xyyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 782);

    auto g_0_xyyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 783);

    auto g_0_xyyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 784);

    auto g_0_xyyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 785);

    auto g_0_xyyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 786);

    auto g_0_xyyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 787);

    auto g_0_xyyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 788);

    auto g_0_xyyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 789);

    auto g_0_xyyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 790);

    auto g_0_xyyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 791);

    auto g_0_xyyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 792);

    auto g_0_xyyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 793);

    auto g_0_xyyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 794);

    auto g_0_xyyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 795);

    auto g_0_xyyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 796);

    auto g_0_xyyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 797);

    auto g_0_xyyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 798);

    auto g_0_xyyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 799);

    auto g_0_xyyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 800);

    auto g_0_xyyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 801);

    auto g_0_xyyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 802);

    auto g_0_xyyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 803);

    auto g_0_xyyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 804);

    auto g_0_xyyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 805);

    auto g_0_xyyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 806);

    auto g_0_xyyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 807);

    auto g_0_xyyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 808);

    auto g_0_xyyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 809);

#pragma omp simd aligned(g_0_xyyyzz_0_xxxxxxxx_0,     \
                             g_0_xyyyzz_0_xxxxxxxy_0, \
                             g_0_xyyyzz_0_xxxxxxxz_0, \
                             g_0_xyyyzz_0_xxxxxxyy_0, \
                             g_0_xyyyzz_0_xxxxxxyz_0, \
                             g_0_xyyyzz_0_xxxxxxzz_0, \
                             g_0_xyyyzz_0_xxxxxyyy_0, \
                             g_0_xyyyzz_0_xxxxxyyz_0, \
                             g_0_xyyyzz_0_xxxxxyzz_0, \
                             g_0_xyyyzz_0_xxxxxzzz_0, \
                             g_0_xyyyzz_0_xxxxyyyy_0, \
                             g_0_xyyyzz_0_xxxxyyyz_0, \
                             g_0_xyyyzz_0_xxxxyyzz_0, \
                             g_0_xyyyzz_0_xxxxyzzz_0, \
                             g_0_xyyyzz_0_xxxxzzzz_0, \
                             g_0_xyyyzz_0_xxxyyyyy_0, \
                             g_0_xyyyzz_0_xxxyyyyz_0, \
                             g_0_xyyyzz_0_xxxyyyzz_0, \
                             g_0_xyyyzz_0_xxxyyzzz_0, \
                             g_0_xyyyzz_0_xxxyzzzz_0, \
                             g_0_xyyyzz_0_xxxzzzzz_0, \
                             g_0_xyyyzz_0_xxyyyyyy_0, \
                             g_0_xyyyzz_0_xxyyyyyz_0, \
                             g_0_xyyyzz_0_xxyyyyzz_0, \
                             g_0_xyyyzz_0_xxyyyzzz_0, \
                             g_0_xyyyzz_0_xxyyzzzz_0, \
                             g_0_xyyyzz_0_xxyzzzzz_0, \
                             g_0_xyyyzz_0_xxzzzzzz_0, \
                             g_0_xyyyzz_0_xyyyyyyy_0, \
                             g_0_xyyyzz_0_xyyyyyyz_0, \
                             g_0_xyyyzz_0_xyyyyyzz_0, \
                             g_0_xyyyzz_0_xyyyyzzz_0, \
                             g_0_xyyyzz_0_xyyyzzzz_0, \
                             g_0_xyyyzz_0_xyyzzzzz_0, \
                             g_0_xyyyzz_0_xyzzzzzz_0, \
                             g_0_xyyyzz_0_xzzzzzzz_0, \
                             g_0_xyyyzz_0_yyyyyyyy_0, \
                             g_0_xyyyzz_0_yyyyyyyz_0, \
                             g_0_xyyyzz_0_yyyyyyzz_0, \
                             g_0_xyyyzz_0_yyyyyzzz_0, \
                             g_0_xyyyzz_0_yyyyzzzz_0, \
                             g_0_xyyyzz_0_yyyzzzzz_0, \
                             g_0_xyyyzz_0_yyzzzzzz_0, \
                             g_0_xyyyzz_0_yzzzzzzz_0, \
                             g_0_xyyyzz_0_zzzzzzzz_0, \
                             g_0_yyyzz_0_xxxxxxx_1,   \
                             g_0_yyyzz_0_xxxxxxxx_0,  \
                             g_0_yyyzz_0_xxxxxxxx_1,  \
                             g_0_yyyzz_0_xxxxxxxy_0,  \
                             g_0_yyyzz_0_xxxxxxxy_1,  \
                             g_0_yyyzz_0_xxxxxxxz_0,  \
                             g_0_yyyzz_0_xxxxxxxz_1,  \
                             g_0_yyyzz_0_xxxxxxy_1,   \
                             g_0_yyyzz_0_xxxxxxyy_0,  \
                             g_0_yyyzz_0_xxxxxxyy_1,  \
                             g_0_yyyzz_0_xxxxxxyz_0,  \
                             g_0_yyyzz_0_xxxxxxyz_1,  \
                             g_0_yyyzz_0_xxxxxxz_1,   \
                             g_0_yyyzz_0_xxxxxxzz_0,  \
                             g_0_yyyzz_0_xxxxxxzz_1,  \
                             g_0_yyyzz_0_xxxxxyy_1,   \
                             g_0_yyyzz_0_xxxxxyyy_0,  \
                             g_0_yyyzz_0_xxxxxyyy_1,  \
                             g_0_yyyzz_0_xxxxxyyz_0,  \
                             g_0_yyyzz_0_xxxxxyyz_1,  \
                             g_0_yyyzz_0_xxxxxyz_1,   \
                             g_0_yyyzz_0_xxxxxyzz_0,  \
                             g_0_yyyzz_0_xxxxxyzz_1,  \
                             g_0_yyyzz_0_xxxxxzz_1,   \
                             g_0_yyyzz_0_xxxxxzzz_0,  \
                             g_0_yyyzz_0_xxxxxzzz_1,  \
                             g_0_yyyzz_0_xxxxyyy_1,   \
                             g_0_yyyzz_0_xxxxyyyy_0,  \
                             g_0_yyyzz_0_xxxxyyyy_1,  \
                             g_0_yyyzz_0_xxxxyyyz_0,  \
                             g_0_yyyzz_0_xxxxyyyz_1,  \
                             g_0_yyyzz_0_xxxxyyz_1,   \
                             g_0_yyyzz_0_xxxxyyzz_0,  \
                             g_0_yyyzz_0_xxxxyyzz_1,  \
                             g_0_yyyzz_0_xxxxyzz_1,   \
                             g_0_yyyzz_0_xxxxyzzz_0,  \
                             g_0_yyyzz_0_xxxxyzzz_1,  \
                             g_0_yyyzz_0_xxxxzzz_1,   \
                             g_0_yyyzz_0_xxxxzzzz_0,  \
                             g_0_yyyzz_0_xxxxzzzz_1,  \
                             g_0_yyyzz_0_xxxyyyy_1,   \
                             g_0_yyyzz_0_xxxyyyyy_0,  \
                             g_0_yyyzz_0_xxxyyyyy_1,  \
                             g_0_yyyzz_0_xxxyyyyz_0,  \
                             g_0_yyyzz_0_xxxyyyyz_1,  \
                             g_0_yyyzz_0_xxxyyyz_1,   \
                             g_0_yyyzz_0_xxxyyyzz_0,  \
                             g_0_yyyzz_0_xxxyyyzz_1,  \
                             g_0_yyyzz_0_xxxyyzz_1,   \
                             g_0_yyyzz_0_xxxyyzzz_0,  \
                             g_0_yyyzz_0_xxxyyzzz_1,  \
                             g_0_yyyzz_0_xxxyzzz_1,   \
                             g_0_yyyzz_0_xxxyzzzz_0,  \
                             g_0_yyyzz_0_xxxyzzzz_1,  \
                             g_0_yyyzz_0_xxxzzzz_1,   \
                             g_0_yyyzz_0_xxxzzzzz_0,  \
                             g_0_yyyzz_0_xxxzzzzz_1,  \
                             g_0_yyyzz_0_xxyyyyy_1,   \
                             g_0_yyyzz_0_xxyyyyyy_0,  \
                             g_0_yyyzz_0_xxyyyyyy_1,  \
                             g_0_yyyzz_0_xxyyyyyz_0,  \
                             g_0_yyyzz_0_xxyyyyyz_1,  \
                             g_0_yyyzz_0_xxyyyyz_1,   \
                             g_0_yyyzz_0_xxyyyyzz_0,  \
                             g_0_yyyzz_0_xxyyyyzz_1,  \
                             g_0_yyyzz_0_xxyyyzz_1,   \
                             g_0_yyyzz_0_xxyyyzzz_0,  \
                             g_0_yyyzz_0_xxyyyzzz_1,  \
                             g_0_yyyzz_0_xxyyzzz_1,   \
                             g_0_yyyzz_0_xxyyzzzz_0,  \
                             g_0_yyyzz_0_xxyyzzzz_1,  \
                             g_0_yyyzz_0_xxyzzzz_1,   \
                             g_0_yyyzz_0_xxyzzzzz_0,  \
                             g_0_yyyzz_0_xxyzzzzz_1,  \
                             g_0_yyyzz_0_xxzzzzz_1,   \
                             g_0_yyyzz_0_xxzzzzzz_0,  \
                             g_0_yyyzz_0_xxzzzzzz_1,  \
                             g_0_yyyzz_0_xyyyyyy_1,   \
                             g_0_yyyzz_0_xyyyyyyy_0,  \
                             g_0_yyyzz_0_xyyyyyyy_1,  \
                             g_0_yyyzz_0_xyyyyyyz_0,  \
                             g_0_yyyzz_0_xyyyyyyz_1,  \
                             g_0_yyyzz_0_xyyyyyz_1,   \
                             g_0_yyyzz_0_xyyyyyzz_0,  \
                             g_0_yyyzz_0_xyyyyyzz_1,  \
                             g_0_yyyzz_0_xyyyyzz_1,   \
                             g_0_yyyzz_0_xyyyyzzz_0,  \
                             g_0_yyyzz_0_xyyyyzzz_1,  \
                             g_0_yyyzz_0_xyyyzzz_1,   \
                             g_0_yyyzz_0_xyyyzzzz_0,  \
                             g_0_yyyzz_0_xyyyzzzz_1,  \
                             g_0_yyyzz_0_xyyzzzz_1,   \
                             g_0_yyyzz_0_xyyzzzzz_0,  \
                             g_0_yyyzz_0_xyyzzzzz_1,  \
                             g_0_yyyzz_0_xyzzzzz_1,   \
                             g_0_yyyzz_0_xyzzzzzz_0,  \
                             g_0_yyyzz_0_xyzzzzzz_1,  \
                             g_0_yyyzz_0_xzzzzzz_1,   \
                             g_0_yyyzz_0_xzzzzzzz_0,  \
                             g_0_yyyzz_0_xzzzzzzz_1,  \
                             g_0_yyyzz_0_yyyyyyy_1,   \
                             g_0_yyyzz_0_yyyyyyyy_0,  \
                             g_0_yyyzz_0_yyyyyyyy_1,  \
                             g_0_yyyzz_0_yyyyyyyz_0,  \
                             g_0_yyyzz_0_yyyyyyyz_1,  \
                             g_0_yyyzz_0_yyyyyyz_1,   \
                             g_0_yyyzz_0_yyyyyyzz_0,  \
                             g_0_yyyzz_0_yyyyyyzz_1,  \
                             g_0_yyyzz_0_yyyyyzz_1,   \
                             g_0_yyyzz_0_yyyyyzzz_0,  \
                             g_0_yyyzz_0_yyyyyzzz_1,  \
                             g_0_yyyzz_0_yyyyzzz_1,   \
                             g_0_yyyzz_0_yyyyzzzz_0,  \
                             g_0_yyyzz_0_yyyyzzzz_1,  \
                             g_0_yyyzz_0_yyyzzzz_1,   \
                             g_0_yyyzz_0_yyyzzzzz_0,  \
                             g_0_yyyzz_0_yyyzzzzz_1,  \
                             g_0_yyyzz_0_yyzzzzz_1,   \
                             g_0_yyyzz_0_yyzzzzzz_0,  \
                             g_0_yyyzz_0_yyzzzzzz_1,  \
                             g_0_yyyzz_0_yzzzzzz_1,   \
                             g_0_yyyzz_0_yzzzzzzz_0,  \
                             g_0_yyyzz_0_yzzzzzzz_1,  \
                             g_0_yyyzz_0_zzzzzzz_1,   \
                             g_0_yyyzz_0_zzzzzzzz_0,  \
                             g_0_yyyzz_0_zzzzzzzz_1,  \
                             wp_x,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzz_0_xxxxxxxx_0[i] =
            8.0 * g_0_yyyzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxxxx_0[i] * pb_x + g_0_yyyzz_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxxxy_0[i] =
            7.0 * g_0_yyyzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxxxy_0[i] * pb_x + g_0_yyyzz_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxxxz_0[i] =
            7.0 * g_0_yyyzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxxxz_0[i] * pb_x + g_0_yyyzz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxxyy_0[i] =
            6.0 * g_0_yyyzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxxyy_0[i] * pb_x + g_0_yyyzz_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxxyz_0[i] =
            6.0 * g_0_yyyzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxxyz_0[i] * pb_x + g_0_yyyzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxxzz_0[i] =
            6.0 * g_0_yyyzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxxzz_0[i] * pb_x + g_0_yyyzz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxyyy_0[i] =
            5.0 * g_0_yyyzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxyyy_0[i] * pb_x + g_0_yyyzz_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxyyz_0[i] =
            5.0 * g_0_yyyzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxyyz_0[i] * pb_x + g_0_yyyzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxyzz_0[i] =
            5.0 * g_0_yyyzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxyzz_0[i] * pb_x + g_0_yyyzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxzzz_0[i] =
            5.0 * g_0_yyyzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxzzz_0[i] * pb_x + g_0_yyyzz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxyyyy_0[i] =
            4.0 * g_0_yyyzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyyyy_0[i] * pb_x + g_0_yyyzz_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxyyyz_0[i] =
            4.0 * g_0_yyyzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyyyz_0[i] * pb_x + g_0_yyyzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxyyzz_0[i] =
            4.0 * g_0_yyyzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyyzz_0[i] * pb_x + g_0_yyyzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxyzzz_0[i] =
            4.0 * g_0_yyyzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyzzz_0[i] * pb_x + g_0_yyyzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxzzzz_0[i] =
            4.0 * g_0_yyyzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxzzzz_0[i] * pb_x + g_0_yyyzz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyyyyy_0[i] =
            3.0 * g_0_yyyzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyyyy_0[i] * pb_x + g_0_yyyzz_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyyyyz_0[i] =
            3.0 * g_0_yyyzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyyyz_0[i] * pb_x + g_0_yyyzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyyyzz_0[i] =
            3.0 * g_0_yyyzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyyzz_0[i] * pb_x + g_0_yyyzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyyzzz_0[i] =
            3.0 * g_0_yyyzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyzzz_0[i] * pb_x + g_0_yyyzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyzzzz_0[i] =
            3.0 * g_0_yyyzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyzzzz_0[i] * pb_x + g_0_yyyzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxzzzzz_0[i] =
            3.0 * g_0_yyyzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxzzzzz_0[i] * pb_x + g_0_yyyzz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyyyyy_0[i] =
            2.0 * g_0_yyyzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyyyy_0[i] * pb_x + g_0_yyyzz_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyyyyz_0[i] =
            2.0 * g_0_yyyzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyyyz_0[i] * pb_x + g_0_yyyzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyyyzz_0[i] =
            2.0 * g_0_yyyzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyyzz_0[i] * pb_x + g_0_yyyzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyyzzz_0[i] =
            2.0 * g_0_yyyzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyzzz_0[i] * pb_x + g_0_yyyzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyzzzz_0[i] =
            2.0 * g_0_yyyzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyzzzz_0[i] * pb_x + g_0_yyyzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyzzzzz_0[i] =
            2.0 * g_0_yyyzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyzzzzz_0[i] * pb_x + g_0_yyyzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxzzzzzz_0[i] =
            2.0 * g_0_yyyzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxzzzzzz_0[i] * pb_x + g_0_yyyzz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyyyyy_0[i] = g_0_yyyzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyyyy_0[i] * pb_x + g_0_yyyzz_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyyyyz_0[i] = g_0_yyyzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyyyz_0[i] * pb_x + g_0_yyyzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyyyzz_0[i] = g_0_yyyzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyyzz_0[i] * pb_x + g_0_yyyzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyyzzz_0[i] = g_0_yyyzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyzzz_0[i] * pb_x + g_0_yyyzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyzzzz_0[i] = g_0_yyyzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyzzzz_0[i] * pb_x + g_0_yyyzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyzzzzz_0[i] = g_0_yyyzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyzzzzz_0[i] * pb_x + g_0_yyyzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyzzzzzz_0[i] = g_0_yyyzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyzzzzzz_0[i] * pb_x + g_0_yyyzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xzzzzzzz_0[i] = g_0_yyyzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xzzzzzzz_0[i] * pb_x + g_0_yyyzz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyyyyy_0[i] = g_0_yyyzz_0_yyyyyyyy_0[i] * pb_x + g_0_yyyzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyyyyz_0[i] = g_0_yyyzz_0_yyyyyyyz_0[i] * pb_x + g_0_yyyzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyyyzz_0[i] = g_0_yyyzz_0_yyyyyyzz_0[i] * pb_x + g_0_yyyzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyyzzz_0[i] = g_0_yyyzz_0_yyyyyzzz_0[i] * pb_x + g_0_yyyzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyzzzz_0[i] = g_0_yyyzz_0_yyyyzzzz_0[i] * pb_x + g_0_yyyzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyzzzzz_0[i] = g_0_yyyzz_0_yyyzzzzz_0[i] * pb_x + g_0_yyyzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyzzzzzz_0[i] = g_0_yyyzz_0_yyzzzzzz_0[i] * pb_x + g_0_yyyzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yzzzzzzz_0[i] = g_0_yyyzz_0_yzzzzzzz_0[i] * pb_x + g_0_yyyzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_zzzzzzzz_0[i] = g_0_yyyzz_0_zzzzzzzz_0[i] * pb_x + g_0_yyyzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 810-855 components of targeted buffer : SISL

    auto g_0_xyyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 810);

    auto g_0_xyyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 811);

    auto g_0_xyyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 812);

    auto g_0_xyyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 813);

    auto g_0_xyyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 814);

    auto g_0_xyyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 815);

    auto g_0_xyyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 816);

    auto g_0_xyyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 817);

    auto g_0_xyyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 818);

    auto g_0_xyyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 819);

    auto g_0_xyyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 820);

    auto g_0_xyyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 821);

    auto g_0_xyyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 822);

    auto g_0_xyyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 823);

    auto g_0_xyyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 824);

    auto g_0_xyyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 825);

    auto g_0_xyyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 826);

    auto g_0_xyyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 827);

    auto g_0_xyyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 828);

    auto g_0_xyyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 829);

    auto g_0_xyyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 830);

    auto g_0_xyyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 831);

    auto g_0_xyyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 832);

    auto g_0_xyyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 833);

    auto g_0_xyyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 834);

    auto g_0_xyyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 835);

    auto g_0_xyyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 836);

    auto g_0_xyyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 837);

    auto g_0_xyyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 838);

    auto g_0_xyyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 839);

    auto g_0_xyyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 840);

    auto g_0_xyyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 841);

    auto g_0_xyyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 842);

    auto g_0_xyyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 843);

    auto g_0_xyyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 844);

    auto g_0_xyyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 845);

    auto g_0_xyyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 846);

    auto g_0_xyyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 847);

    auto g_0_xyyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 848);

    auto g_0_xyyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 849);

    auto g_0_xyyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 850);

    auto g_0_xyyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 851);

    auto g_0_xyyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 852);

    auto g_0_xyyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 853);

    auto g_0_xyyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 854);

#pragma omp simd aligned(g_0_xyyzzz_0_xxxxxxxx_0,     \
                             g_0_xyyzzz_0_xxxxxxxy_0, \
                             g_0_xyyzzz_0_xxxxxxxz_0, \
                             g_0_xyyzzz_0_xxxxxxyy_0, \
                             g_0_xyyzzz_0_xxxxxxyz_0, \
                             g_0_xyyzzz_0_xxxxxxzz_0, \
                             g_0_xyyzzz_0_xxxxxyyy_0, \
                             g_0_xyyzzz_0_xxxxxyyz_0, \
                             g_0_xyyzzz_0_xxxxxyzz_0, \
                             g_0_xyyzzz_0_xxxxxzzz_0, \
                             g_0_xyyzzz_0_xxxxyyyy_0, \
                             g_0_xyyzzz_0_xxxxyyyz_0, \
                             g_0_xyyzzz_0_xxxxyyzz_0, \
                             g_0_xyyzzz_0_xxxxyzzz_0, \
                             g_0_xyyzzz_0_xxxxzzzz_0, \
                             g_0_xyyzzz_0_xxxyyyyy_0, \
                             g_0_xyyzzz_0_xxxyyyyz_0, \
                             g_0_xyyzzz_0_xxxyyyzz_0, \
                             g_0_xyyzzz_0_xxxyyzzz_0, \
                             g_0_xyyzzz_0_xxxyzzzz_0, \
                             g_0_xyyzzz_0_xxxzzzzz_0, \
                             g_0_xyyzzz_0_xxyyyyyy_0, \
                             g_0_xyyzzz_0_xxyyyyyz_0, \
                             g_0_xyyzzz_0_xxyyyyzz_0, \
                             g_0_xyyzzz_0_xxyyyzzz_0, \
                             g_0_xyyzzz_0_xxyyzzzz_0, \
                             g_0_xyyzzz_0_xxyzzzzz_0, \
                             g_0_xyyzzz_0_xxzzzzzz_0, \
                             g_0_xyyzzz_0_xyyyyyyy_0, \
                             g_0_xyyzzz_0_xyyyyyyz_0, \
                             g_0_xyyzzz_0_xyyyyyzz_0, \
                             g_0_xyyzzz_0_xyyyyzzz_0, \
                             g_0_xyyzzz_0_xyyyzzzz_0, \
                             g_0_xyyzzz_0_xyyzzzzz_0, \
                             g_0_xyyzzz_0_xyzzzzzz_0, \
                             g_0_xyyzzz_0_xzzzzzzz_0, \
                             g_0_xyyzzz_0_yyyyyyyy_0, \
                             g_0_xyyzzz_0_yyyyyyyz_0, \
                             g_0_xyyzzz_0_yyyyyyzz_0, \
                             g_0_xyyzzz_0_yyyyyzzz_0, \
                             g_0_xyyzzz_0_yyyyzzzz_0, \
                             g_0_xyyzzz_0_yyyzzzzz_0, \
                             g_0_xyyzzz_0_yyzzzzzz_0, \
                             g_0_xyyzzz_0_yzzzzzzz_0, \
                             g_0_xyyzzz_0_zzzzzzzz_0, \
                             g_0_yyzzz_0_xxxxxxx_1,   \
                             g_0_yyzzz_0_xxxxxxxx_0,  \
                             g_0_yyzzz_0_xxxxxxxx_1,  \
                             g_0_yyzzz_0_xxxxxxxy_0,  \
                             g_0_yyzzz_0_xxxxxxxy_1,  \
                             g_0_yyzzz_0_xxxxxxxz_0,  \
                             g_0_yyzzz_0_xxxxxxxz_1,  \
                             g_0_yyzzz_0_xxxxxxy_1,   \
                             g_0_yyzzz_0_xxxxxxyy_0,  \
                             g_0_yyzzz_0_xxxxxxyy_1,  \
                             g_0_yyzzz_0_xxxxxxyz_0,  \
                             g_0_yyzzz_0_xxxxxxyz_1,  \
                             g_0_yyzzz_0_xxxxxxz_1,   \
                             g_0_yyzzz_0_xxxxxxzz_0,  \
                             g_0_yyzzz_0_xxxxxxzz_1,  \
                             g_0_yyzzz_0_xxxxxyy_1,   \
                             g_0_yyzzz_0_xxxxxyyy_0,  \
                             g_0_yyzzz_0_xxxxxyyy_1,  \
                             g_0_yyzzz_0_xxxxxyyz_0,  \
                             g_0_yyzzz_0_xxxxxyyz_1,  \
                             g_0_yyzzz_0_xxxxxyz_1,   \
                             g_0_yyzzz_0_xxxxxyzz_0,  \
                             g_0_yyzzz_0_xxxxxyzz_1,  \
                             g_0_yyzzz_0_xxxxxzz_1,   \
                             g_0_yyzzz_0_xxxxxzzz_0,  \
                             g_0_yyzzz_0_xxxxxzzz_1,  \
                             g_0_yyzzz_0_xxxxyyy_1,   \
                             g_0_yyzzz_0_xxxxyyyy_0,  \
                             g_0_yyzzz_0_xxxxyyyy_1,  \
                             g_0_yyzzz_0_xxxxyyyz_0,  \
                             g_0_yyzzz_0_xxxxyyyz_1,  \
                             g_0_yyzzz_0_xxxxyyz_1,   \
                             g_0_yyzzz_0_xxxxyyzz_0,  \
                             g_0_yyzzz_0_xxxxyyzz_1,  \
                             g_0_yyzzz_0_xxxxyzz_1,   \
                             g_0_yyzzz_0_xxxxyzzz_0,  \
                             g_0_yyzzz_0_xxxxyzzz_1,  \
                             g_0_yyzzz_0_xxxxzzz_1,   \
                             g_0_yyzzz_0_xxxxzzzz_0,  \
                             g_0_yyzzz_0_xxxxzzzz_1,  \
                             g_0_yyzzz_0_xxxyyyy_1,   \
                             g_0_yyzzz_0_xxxyyyyy_0,  \
                             g_0_yyzzz_0_xxxyyyyy_1,  \
                             g_0_yyzzz_0_xxxyyyyz_0,  \
                             g_0_yyzzz_0_xxxyyyyz_1,  \
                             g_0_yyzzz_0_xxxyyyz_1,   \
                             g_0_yyzzz_0_xxxyyyzz_0,  \
                             g_0_yyzzz_0_xxxyyyzz_1,  \
                             g_0_yyzzz_0_xxxyyzz_1,   \
                             g_0_yyzzz_0_xxxyyzzz_0,  \
                             g_0_yyzzz_0_xxxyyzzz_1,  \
                             g_0_yyzzz_0_xxxyzzz_1,   \
                             g_0_yyzzz_0_xxxyzzzz_0,  \
                             g_0_yyzzz_0_xxxyzzzz_1,  \
                             g_0_yyzzz_0_xxxzzzz_1,   \
                             g_0_yyzzz_0_xxxzzzzz_0,  \
                             g_0_yyzzz_0_xxxzzzzz_1,  \
                             g_0_yyzzz_0_xxyyyyy_1,   \
                             g_0_yyzzz_0_xxyyyyyy_0,  \
                             g_0_yyzzz_0_xxyyyyyy_1,  \
                             g_0_yyzzz_0_xxyyyyyz_0,  \
                             g_0_yyzzz_0_xxyyyyyz_1,  \
                             g_0_yyzzz_0_xxyyyyz_1,   \
                             g_0_yyzzz_0_xxyyyyzz_0,  \
                             g_0_yyzzz_0_xxyyyyzz_1,  \
                             g_0_yyzzz_0_xxyyyzz_1,   \
                             g_0_yyzzz_0_xxyyyzzz_0,  \
                             g_0_yyzzz_0_xxyyyzzz_1,  \
                             g_0_yyzzz_0_xxyyzzz_1,   \
                             g_0_yyzzz_0_xxyyzzzz_0,  \
                             g_0_yyzzz_0_xxyyzzzz_1,  \
                             g_0_yyzzz_0_xxyzzzz_1,   \
                             g_0_yyzzz_0_xxyzzzzz_0,  \
                             g_0_yyzzz_0_xxyzzzzz_1,  \
                             g_0_yyzzz_0_xxzzzzz_1,   \
                             g_0_yyzzz_0_xxzzzzzz_0,  \
                             g_0_yyzzz_0_xxzzzzzz_1,  \
                             g_0_yyzzz_0_xyyyyyy_1,   \
                             g_0_yyzzz_0_xyyyyyyy_0,  \
                             g_0_yyzzz_0_xyyyyyyy_1,  \
                             g_0_yyzzz_0_xyyyyyyz_0,  \
                             g_0_yyzzz_0_xyyyyyyz_1,  \
                             g_0_yyzzz_0_xyyyyyz_1,   \
                             g_0_yyzzz_0_xyyyyyzz_0,  \
                             g_0_yyzzz_0_xyyyyyzz_1,  \
                             g_0_yyzzz_0_xyyyyzz_1,   \
                             g_0_yyzzz_0_xyyyyzzz_0,  \
                             g_0_yyzzz_0_xyyyyzzz_1,  \
                             g_0_yyzzz_0_xyyyzzz_1,   \
                             g_0_yyzzz_0_xyyyzzzz_0,  \
                             g_0_yyzzz_0_xyyyzzzz_1,  \
                             g_0_yyzzz_0_xyyzzzz_1,   \
                             g_0_yyzzz_0_xyyzzzzz_0,  \
                             g_0_yyzzz_0_xyyzzzzz_1,  \
                             g_0_yyzzz_0_xyzzzzz_1,   \
                             g_0_yyzzz_0_xyzzzzzz_0,  \
                             g_0_yyzzz_0_xyzzzzzz_1,  \
                             g_0_yyzzz_0_xzzzzzz_1,   \
                             g_0_yyzzz_0_xzzzzzzz_0,  \
                             g_0_yyzzz_0_xzzzzzzz_1,  \
                             g_0_yyzzz_0_yyyyyyy_1,   \
                             g_0_yyzzz_0_yyyyyyyy_0,  \
                             g_0_yyzzz_0_yyyyyyyy_1,  \
                             g_0_yyzzz_0_yyyyyyyz_0,  \
                             g_0_yyzzz_0_yyyyyyyz_1,  \
                             g_0_yyzzz_0_yyyyyyz_1,   \
                             g_0_yyzzz_0_yyyyyyzz_0,  \
                             g_0_yyzzz_0_yyyyyyzz_1,  \
                             g_0_yyzzz_0_yyyyyzz_1,   \
                             g_0_yyzzz_0_yyyyyzzz_0,  \
                             g_0_yyzzz_0_yyyyyzzz_1,  \
                             g_0_yyzzz_0_yyyyzzz_1,   \
                             g_0_yyzzz_0_yyyyzzzz_0,  \
                             g_0_yyzzz_0_yyyyzzzz_1,  \
                             g_0_yyzzz_0_yyyzzzz_1,   \
                             g_0_yyzzz_0_yyyzzzzz_0,  \
                             g_0_yyzzz_0_yyyzzzzz_1,  \
                             g_0_yyzzz_0_yyzzzzz_1,   \
                             g_0_yyzzz_0_yyzzzzzz_0,  \
                             g_0_yyzzz_0_yyzzzzzz_1,  \
                             g_0_yyzzz_0_yzzzzzz_1,   \
                             g_0_yyzzz_0_yzzzzzzz_0,  \
                             g_0_yyzzz_0_yzzzzzzz_1,  \
                             g_0_yyzzz_0_zzzzzzz_1,   \
                             g_0_yyzzz_0_zzzzzzzz_0,  \
                             g_0_yyzzz_0_zzzzzzzz_1,  \
                             wp_x,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzz_0_xxxxxxxx_0[i] =
            8.0 * g_0_yyzzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxxxx_0[i] * pb_x + g_0_yyzzz_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxxxy_0[i] =
            7.0 * g_0_yyzzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxxxy_0[i] * pb_x + g_0_yyzzz_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxxxz_0[i] =
            7.0 * g_0_yyzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxxxz_0[i] * pb_x + g_0_yyzzz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxxyy_0[i] =
            6.0 * g_0_yyzzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxxyy_0[i] * pb_x + g_0_yyzzz_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxxyz_0[i] =
            6.0 * g_0_yyzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxxyz_0[i] * pb_x + g_0_yyzzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxxzz_0[i] =
            6.0 * g_0_yyzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxxzz_0[i] * pb_x + g_0_yyzzz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxyyy_0[i] =
            5.0 * g_0_yyzzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxyyy_0[i] * pb_x + g_0_yyzzz_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxyyz_0[i] =
            5.0 * g_0_yyzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxyyz_0[i] * pb_x + g_0_yyzzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxyzz_0[i] =
            5.0 * g_0_yyzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxyzz_0[i] * pb_x + g_0_yyzzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxzzz_0[i] =
            5.0 * g_0_yyzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxzzz_0[i] * pb_x + g_0_yyzzz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxyyyy_0[i] =
            4.0 * g_0_yyzzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyyyy_0[i] * pb_x + g_0_yyzzz_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxyyyz_0[i] =
            4.0 * g_0_yyzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyyyz_0[i] * pb_x + g_0_yyzzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxyyzz_0[i] =
            4.0 * g_0_yyzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyyzz_0[i] * pb_x + g_0_yyzzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxyzzz_0[i] =
            4.0 * g_0_yyzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyzzz_0[i] * pb_x + g_0_yyzzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxzzzz_0[i] =
            4.0 * g_0_yyzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxzzzz_0[i] * pb_x + g_0_yyzzz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyyyyy_0[i] =
            3.0 * g_0_yyzzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyyyy_0[i] * pb_x + g_0_yyzzz_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyyyyz_0[i] =
            3.0 * g_0_yyzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyyyz_0[i] * pb_x + g_0_yyzzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyyyzz_0[i] =
            3.0 * g_0_yyzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyyzz_0[i] * pb_x + g_0_yyzzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyyzzz_0[i] =
            3.0 * g_0_yyzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyzzz_0[i] * pb_x + g_0_yyzzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyzzzz_0[i] =
            3.0 * g_0_yyzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyzzzz_0[i] * pb_x + g_0_yyzzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxzzzzz_0[i] =
            3.0 * g_0_yyzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxzzzzz_0[i] * pb_x + g_0_yyzzz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyyyyy_0[i] =
            2.0 * g_0_yyzzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyyyy_0[i] * pb_x + g_0_yyzzz_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyyyyz_0[i] =
            2.0 * g_0_yyzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyyyz_0[i] * pb_x + g_0_yyzzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyyyzz_0[i] =
            2.0 * g_0_yyzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyyzz_0[i] * pb_x + g_0_yyzzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyyzzz_0[i] =
            2.0 * g_0_yyzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyzzz_0[i] * pb_x + g_0_yyzzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyzzzz_0[i] =
            2.0 * g_0_yyzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyzzzz_0[i] * pb_x + g_0_yyzzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyzzzzz_0[i] =
            2.0 * g_0_yyzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyzzzzz_0[i] * pb_x + g_0_yyzzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxzzzzzz_0[i] =
            2.0 * g_0_yyzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxzzzzzz_0[i] * pb_x + g_0_yyzzz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyyyyy_0[i] = g_0_yyzzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyyyy_0[i] * pb_x + g_0_yyzzz_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyyyyz_0[i] = g_0_yyzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyyyz_0[i] * pb_x + g_0_yyzzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyyyzz_0[i] = g_0_yyzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyyzz_0[i] * pb_x + g_0_yyzzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyyzzz_0[i] = g_0_yyzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyzzz_0[i] * pb_x + g_0_yyzzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyzzzz_0[i] = g_0_yyzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyzzzz_0[i] * pb_x + g_0_yyzzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyzzzzz_0[i] = g_0_yyzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyzzzzz_0[i] * pb_x + g_0_yyzzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyzzzzzz_0[i] = g_0_yyzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyzzzzzz_0[i] * pb_x + g_0_yyzzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xzzzzzzz_0[i] = g_0_yyzzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xzzzzzzz_0[i] * pb_x + g_0_yyzzz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyyyyy_0[i] = g_0_yyzzz_0_yyyyyyyy_0[i] * pb_x + g_0_yyzzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyyyyz_0[i] = g_0_yyzzz_0_yyyyyyyz_0[i] * pb_x + g_0_yyzzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyyyzz_0[i] = g_0_yyzzz_0_yyyyyyzz_0[i] * pb_x + g_0_yyzzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyyzzz_0[i] = g_0_yyzzz_0_yyyyyzzz_0[i] * pb_x + g_0_yyzzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyzzzz_0[i] = g_0_yyzzz_0_yyyyzzzz_0[i] * pb_x + g_0_yyzzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyzzzzz_0[i] = g_0_yyzzz_0_yyyzzzzz_0[i] * pb_x + g_0_yyzzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyzzzzzz_0[i] = g_0_yyzzz_0_yyzzzzzz_0[i] * pb_x + g_0_yyzzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yzzzzzzz_0[i] = g_0_yyzzz_0_yzzzzzzz_0[i] * pb_x + g_0_yyzzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_zzzzzzzz_0[i] = g_0_yyzzz_0_zzzzzzzz_0[i] * pb_x + g_0_yyzzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 855-900 components of targeted buffer : SISL

    auto g_0_xyzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 855);

    auto g_0_xyzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 856);

    auto g_0_xyzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 857);

    auto g_0_xyzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 858);

    auto g_0_xyzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 859);

    auto g_0_xyzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 860);

    auto g_0_xyzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 861);

    auto g_0_xyzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 862);

    auto g_0_xyzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 863);

    auto g_0_xyzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 864);

    auto g_0_xyzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 865);

    auto g_0_xyzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 866);

    auto g_0_xyzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 867);

    auto g_0_xyzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 868);

    auto g_0_xyzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 869);

    auto g_0_xyzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 870);

    auto g_0_xyzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 871);

    auto g_0_xyzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 872);

    auto g_0_xyzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 873);

    auto g_0_xyzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 874);

    auto g_0_xyzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 875);

    auto g_0_xyzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 876);

    auto g_0_xyzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 877);

    auto g_0_xyzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 878);

    auto g_0_xyzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 879);

    auto g_0_xyzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 880);

    auto g_0_xyzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 881);

    auto g_0_xyzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 882);

    auto g_0_xyzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 883);

    auto g_0_xyzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 884);

    auto g_0_xyzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 885);

    auto g_0_xyzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 886);

    auto g_0_xyzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 887);

    auto g_0_xyzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 888);

    auto g_0_xyzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 889);

    auto g_0_xyzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 890);

    auto g_0_xyzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 891);

    auto g_0_xyzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 892);

    auto g_0_xyzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 893);

    auto g_0_xyzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 894);

    auto g_0_xyzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 895);

    auto g_0_xyzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 896);

    auto g_0_xyzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 897);

    auto g_0_xyzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 898);

    auto g_0_xyzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 899);

#pragma omp simd aligned(g_0_xyzzzz_0_xxxxxxxx_0,     \
                             g_0_xyzzzz_0_xxxxxxxy_0, \
                             g_0_xyzzzz_0_xxxxxxxz_0, \
                             g_0_xyzzzz_0_xxxxxxyy_0, \
                             g_0_xyzzzz_0_xxxxxxyz_0, \
                             g_0_xyzzzz_0_xxxxxxzz_0, \
                             g_0_xyzzzz_0_xxxxxyyy_0, \
                             g_0_xyzzzz_0_xxxxxyyz_0, \
                             g_0_xyzzzz_0_xxxxxyzz_0, \
                             g_0_xyzzzz_0_xxxxxzzz_0, \
                             g_0_xyzzzz_0_xxxxyyyy_0, \
                             g_0_xyzzzz_0_xxxxyyyz_0, \
                             g_0_xyzzzz_0_xxxxyyzz_0, \
                             g_0_xyzzzz_0_xxxxyzzz_0, \
                             g_0_xyzzzz_0_xxxxzzzz_0, \
                             g_0_xyzzzz_0_xxxyyyyy_0, \
                             g_0_xyzzzz_0_xxxyyyyz_0, \
                             g_0_xyzzzz_0_xxxyyyzz_0, \
                             g_0_xyzzzz_0_xxxyyzzz_0, \
                             g_0_xyzzzz_0_xxxyzzzz_0, \
                             g_0_xyzzzz_0_xxxzzzzz_0, \
                             g_0_xyzzzz_0_xxyyyyyy_0, \
                             g_0_xyzzzz_0_xxyyyyyz_0, \
                             g_0_xyzzzz_0_xxyyyyzz_0, \
                             g_0_xyzzzz_0_xxyyyzzz_0, \
                             g_0_xyzzzz_0_xxyyzzzz_0, \
                             g_0_xyzzzz_0_xxyzzzzz_0, \
                             g_0_xyzzzz_0_xxzzzzzz_0, \
                             g_0_xyzzzz_0_xyyyyyyy_0, \
                             g_0_xyzzzz_0_xyyyyyyz_0, \
                             g_0_xyzzzz_0_xyyyyyzz_0, \
                             g_0_xyzzzz_0_xyyyyzzz_0, \
                             g_0_xyzzzz_0_xyyyzzzz_0, \
                             g_0_xyzzzz_0_xyyzzzzz_0, \
                             g_0_xyzzzz_0_xyzzzzzz_0, \
                             g_0_xyzzzz_0_xzzzzzzz_0, \
                             g_0_xyzzzz_0_yyyyyyyy_0, \
                             g_0_xyzzzz_0_yyyyyyyz_0, \
                             g_0_xyzzzz_0_yyyyyyzz_0, \
                             g_0_xyzzzz_0_yyyyyzzz_0, \
                             g_0_xyzzzz_0_yyyyzzzz_0, \
                             g_0_xyzzzz_0_yyyzzzzz_0, \
                             g_0_xyzzzz_0_yyzzzzzz_0, \
                             g_0_xyzzzz_0_yzzzzzzz_0, \
                             g_0_xyzzzz_0_zzzzzzzz_0, \
                             g_0_xzzzz_0_xxxxxxxx_0,  \
                             g_0_xzzzz_0_xxxxxxxx_1,  \
                             g_0_xzzzz_0_xxxxxxxz_0,  \
                             g_0_xzzzz_0_xxxxxxxz_1,  \
                             g_0_xzzzz_0_xxxxxxzz_0,  \
                             g_0_xzzzz_0_xxxxxxzz_1,  \
                             g_0_xzzzz_0_xxxxxzzz_0,  \
                             g_0_xzzzz_0_xxxxxzzz_1,  \
                             g_0_xzzzz_0_xxxxzzzz_0,  \
                             g_0_xzzzz_0_xxxxzzzz_1,  \
                             g_0_xzzzz_0_xxxzzzzz_0,  \
                             g_0_xzzzz_0_xxxzzzzz_1,  \
                             g_0_xzzzz_0_xxzzzzzz_0,  \
                             g_0_xzzzz_0_xxzzzzzz_1,  \
                             g_0_xzzzz_0_xzzzzzzz_0,  \
                             g_0_xzzzz_0_xzzzzzzz_1,  \
                             g_0_yzzzz_0_xxxxxxxy_0,  \
                             g_0_yzzzz_0_xxxxxxxy_1,  \
                             g_0_yzzzz_0_xxxxxxy_1,   \
                             g_0_yzzzz_0_xxxxxxyy_0,  \
                             g_0_yzzzz_0_xxxxxxyy_1,  \
                             g_0_yzzzz_0_xxxxxxyz_0,  \
                             g_0_yzzzz_0_xxxxxxyz_1,  \
                             g_0_yzzzz_0_xxxxxyy_1,   \
                             g_0_yzzzz_0_xxxxxyyy_0,  \
                             g_0_yzzzz_0_xxxxxyyy_1,  \
                             g_0_yzzzz_0_xxxxxyyz_0,  \
                             g_0_yzzzz_0_xxxxxyyz_1,  \
                             g_0_yzzzz_0_xxxxxyz_1,   \
                             g_0_yzzzz_0_xxxxxyzz_0,  \
                             g_0_yzzzz_0_xxxxxyzz_1,  \
                             g_0_yzzzz_0_xxxxyyy_1,   \
                             g_0_yzzzz_0_xxxxyyyy_0,  \
                             g_0_yzzzz_0_xxxxyyyy_1,  \
                             g_0_yzzzz_0_xxxxyyyz_0,  \
                             g_0_yzzzz_0_xxxxyyyz_1,  \
                             g_0_yzzzz_0_xxxxyyz_1,   \
                             g_0_yzzzz_0_xxxxyyzz_0,  \
                             g_0_yzzzz_0_xxxxyyzz_1,  \
                             g_0_yzzzz_0_xxxxyzz_1,   \
                             g_0_yzzzz_0_xxxxyzzz_0,  \
                             g_0_yzzzz_0_xxxxyzzz_1,  \
                             g_0_yzzzz_0_xxxyyyy_1,   \
                             g_0_yzzzz_0_xxxyyyyy_0,  \
                             g_0_yzzzz_0_xxxyyyyy_1,  \
                             g_0_yzzzz_0_xxxyyyyz_0,  \
                             g_0_yzzzz_0_xxxyyyyz_1,  \
                             g_0_yzzzz_0_xxxyyyz_1,   \
                             g_0_yzzzz_0_xxxyyyzz_0,  \
                             g_0_yzzzz_0_xxxyyyzz_1,  \
                             g_0_yzzzz_0_xxxyyzz_1,   \
                             g_0_yzzzz_0_xxxyyzzz_0,  \
                             g_0_yzzzz_0_xxxyyzzz_1,  \
                             g_0_yzzzz_0_xxxyzzz_1,   \
                             g_0_yzzzz_0_xxxyzzzz_0,  \
                             g_0_yzzzz_0_xxxyzzzz_1,  \
                             g_0_yzzzz_0_xxyyyyy_1,   \
                             g_0_yzzzz_0_xxyyyyyy_0,  \
                             g_0_yzzzz_0_xxyyyyyy_1,  \
                             g_0_yzzzz_0_xxyyyyyz_0,  \
                             g_0_yzzzz_0_xxyyyyyz_1,  \
                             g_0_yzzzz_0_xxyyyyz_1,   \
                             g_0_yzzzz_0_xxyyyyzz_0,  \
                             g_0_yzzzz_0_xxyyyyzz_1,  \
                             g_0_yzzzz_0_xxyyyzz_1,   \
                             g_0_yzzzz_0_xxyyyzzz_0,  \
                             g_0_yzzzz_0_xxyyyzzz_1,  \
                             g_0_yzzzz_0_xxyyzzz_1,   \
                             g_0_yzzzz_0_xxyyzzzz_0,  \
                             g_0_yzzzz_0_xxyyzzzz_1,  \
                             g_0_yzzzz_0_xxyzzzz_1,   \
                             g_0_yzzzz_0_xxyzzzzz_0,  \
                             g_0_yzzzz_0_xxyzzzzz_1,  \
                             g_0_yzzzz_0_xyyyyyy_1,   \
                             g_0_yzzzz_0_xyyyyyyy_0,  \
                             g_0_yzzzz_0_xyyyyyyy_1,  \
                             g_0_yzzzz_0_xyyyyyyz_0,  \
                             g_0_yzzzz_0_xyyyyyyz_1,  \
                             g_0_yzzzz_0_xyyyyyz_1,   \
                             g_0_yzzzz_0_xyyyyyzz_0,  \
                             g_0_yzzzz_0_xyyyyyzz_1,  \
                             g_0_yzzzz_0_xyyyyzz_1,   \
                             g_0_yzzzz_0_xyyyyzzz_0,  \
                             g_0_yzzzz_0_xyyyyzzz_1,  \
                             g_0_yzzzz_0_xyyyzzz_1,   \
                             g_0_yzzzz_0_xyyyzzzz_0,  \
                             g_0_yzzzz_0_xyyyzzzz_1,  \
                             g_0_yzzzz_0_xyyzzzz_1,   \
                             g_0_yzzzz_0_xyyzzzzz_0,  \
                             g_0_yzzzz_0_xyyzzzzz_1,  \
                             g_0_yzzzz_0_xyzzzzz_1,   \
                             g_0_yzzzz_0_xyzzzzzz_0,  \
                             g_0_yzzzz_0_xyzzzzzz_1,  \
                             g_0_yzzzz_0_yyyyyyy_1,   \
                             g_0_yzzzz_0_yyyyyyyy_0,  \
                             g_0_yzzzz_0_yyyyyyyy_1,  \
                             g_0_yzzzz_0_yyyyyyyz_0,  \
                             g_0_yzzzz_0_yyyyyyyz_1,  \
                             g_0_yzzzz_0_yyyyyyz_1,   \
                             g_0_yzzzz_0_yyyyyyzz_0,  \
                             g_0_yzzzz_0_yyyyyyzz_1,  \
                             g_0_yzzzz_0_yyyyyzz_1,   \
                             g_0_yzzzz_0_yyyyyzzz_0,  \
                             g_0_yzzzz_0_yyyyyzzz_1,  \
                             g_0_yzzzz_0_yyyyzzz_1,   \
                             g_0_yzzzz_0_yyyyzzzz_0,  \
                             g_0_yzzzz_0_yyyyzzzz_1,  \
                             g_0_yzzzz_0_yyyzzzz_1,   \
                             g_0_yzzzz_0_yyyzzzzz_0,  \
                             g_0_yzzzz_0_yyyzzzzz_1,  \
                             g_0_yzzzz_0_yyzzzzz_1,   \
                             g_0_yzzzz_0_yyzzzzzz_0,  \
                             g_0_yzzzz_0_yyzzzzzz_1,  \
                             g_0_yzzzz_0_yzzzzzz_1,   \
                             g_0_yzzzz_0_yzzzzzzz_0,  \
                             g_0_yzzzz_0_yzzzzzzz_1,  \
                             g_0_yzzzz_0_zzzzzzzz_0,  \
                             g_0_yzzzz_0_zzzzzzzz_1,  \
                             wp_x,                    \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzz_0_xxxxxxxx_0[i] = g_0_xzzzz_0_xxxxxxxx_0[i] * pb_y + g_0_xzzzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxxxxxy_0[i] =
            7.0 * g_0_yzzzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxxxxy_0[i] * pb_x + g_0_yzzzz_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxxxxz_0[i] = g_0_xzzzz_0_xxxxxxxz_0[i] * pb_y + g_0_xzzzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxxxxyy_0[i] =
            6.0 * g_0_yzzzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxxxyy_0[i] * pb_x + g_0_yzzzz_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxxxyz_0[i] =
            6.0 * g_0_yzzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxxxyz_0[i] * pb_x + g_0_yzzzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxxxzz_0[i] = g_0_xzzzz_0_xxxxxxzz_0[i] * pb_y + g_0_xzzzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxxxyyy_0[i] =
            5.0 * g_0_yzzzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxxyyy_0[i] * pb_x + g_0_yzzzz_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxxyyz_0[i] =
            5.0 * g_0_yzzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxxyyz_0[i] * pb_x + g_0_yzzzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxxyzz_0[i] =
            5.0 * g_0_yzzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxxyzz_0[i] * pb_x + g_0_yzzzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxxzzz_0[i] = g_0_xzzzz_0_xxxxxzzz_0[i] * pb_y + g_0_xzzzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxxyyyy_0[i] =
            4.0 * g_0_yzzzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyyyy_0[i] * pb_x + g_0_yzzzz_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxyyyz_0[i] =
            4.0 * g_0_yzzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyyyz_0[i] * pb_x + g_0_yzzzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxyyzz_0[i] =
            4.0 * g_0_yzzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyyzz_0[i] * pb_x + g_0_yzzzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxyzzz_0[i] =
            4.0 * g_0_yzzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyzzz_0[i] * pb_x + g_0_yzzzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxzzzz_0[i] = g_0_xzzzz_0_xxxxzzzz_0[i] * pb_y + g_0_xzzzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxyyyyy_0[i] =
            3.0 * g_0_yzzzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyyyy_0[i] * pb_x + g_0_yzzzz_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxyyyyz_0[i] =
            3.0 * g_0_yzzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyyyz_0[i] * pb_x + g_0_yzzzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxyyyzz_0[i] =
            3.0 * g_0_yzzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyyzz_0[i] * pb_x + g_0_yzzzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxyyzzz_0[i] =
            3.0 * g_0_yzzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyzzz_0[i] * pb_x + g_0_yzzzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxyzzzz_0[i] =
            3.0 * g_0_yzzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyzzzz_0[i] * pb_x + g_0_yzzzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxzzzzz_0[i] = g_0_xzzzz_0_xxxzzzzz_0[i] * pb_y + g_0_xzzzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxyyyyyy_0[i] =
            2.0 * g_0_yzzzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyyyy_0[i] * pb_x + g_0_yzzzz_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyyyyyz_0[i] =
            2.0 * g_0_yzzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyyyz_0[i] * pb_x + g_0_yzzzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyyyyzz_0[i] =
            2.0 * g_0_yzzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyyzz_0[i] * pb_x + g_0_yzzzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyyyzzz_0[i] =
            2.0 * g_0_yzzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyzzz_0[i] * pb_x + g_0_yzzzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyyzzzz_0[i] =
            2.0 * g_0_yzzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyzzzz_0[i] * pb_x + g_0_yzzzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyzzzzz_0[i] =
            2.0 * g_0_yzzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyzzzzz_0[i] * pb_x + g_0_yzzzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxzzzzzz_0[i] = g_0_xzzzz_0_xxzzzzzz_0[i] * pb_y + g_0_xzzzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xyyyyyyy_0[i] = g_0_yzzzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyyyy_0[i] * pb_x + g_0_yzzzz_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyyyyyz_0[i] = g_0_yzzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyyyz_0[i] * pb_x + g_0_yzzzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyyyyzz_0[i] = g_0_yzzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyyzz_0[i] * pb_x + g_0_yzzzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyyyzzz_0[i] = g_0_yzzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyzzz_0[i] * pb_x + g_0_yzzzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyyzzzz_0[i] = g_0_yzzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyzzzz_0[i] * pb_x + g_0_yzzzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyzzzzz_0[i] = g_0_yzzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyzzzzz_0[i] * pb_x + g_0_yzzzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyzzzzzz_0[i] = g_0_yzzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyzzzzzz_0[i] * pb_x + g_0_yzzzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xzzzzzzz_0[i] = g_0_xzzzz_0_xzzzzzzz_0[i] * pb_y + g_0_xzzzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_yyyyyyyy_0[i] = g_0_yzzzz_0_yyyyyyyy_0[i] * pb_x + g_0_yzzzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyyyyyz_0[i] = g_0_yzzzz_0_yyyyyyyz_0[i] * pb_x + g_0_yzzzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyyyyzz_0[i] = g_0_yzzzz_0_yyyyyyzz_0[i] * pb_x + g_0_yzzzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyyyzzz_0[i] = g_0_yzzzz_0_yyyyyzzz_0[i] * pb_x + g_0_yzzzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyyzzzz_0[i] = g_0_yzzzz_0_yyyyzzzz_0[i] * pb_x + g_0_yzzzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyzzzzz_0[i] = g_0_yzzzz_0_yyyzzzzz_0[i] * pb_x + g_0_yzzzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyzzzzzz_0[i] = g_0_yzzzz_0_yyzzzzzz_0[i] * pb_x + g_0_yzzzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yzzzzzzz_0[i] = g_0_yzzzz_0_yzzzzzzz_0[i] * pb_x + g_0_yzzzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_zzzzzzzz_0[i] = g_0_yzzzz_0_zzzzzzzz_0[i] * pb_x + g_0_yzzzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 900-945 components of targeted buffer : SISL

    auto g_0_xzzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 900);

    auto g_0_xzzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 901);

    auto g_0_xzzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 902);

    auto g_0_xzzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 903);

    auto g_0_xzzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 904);

    auto g_0_xzzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 905);

    auto g_0_xzzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 906);

    auto g_0_xzzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 907);

    auto g_0_xzzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 908);

    auto g_0_xzzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 909);

    auto g_0_xzzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 910);

    auto g_0_xzzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 911);

    auto g_0_xzzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 912);

    auto g_0_xzzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 913);

    auto g_0_xzzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 914);

    auto g_0_xzzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 915);

    auto g_0_xzzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 916);

    auto g_0_xzzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 917);

    auto g_0_xzzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 918);

    auto g_0_xzzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 919);

    auto g_0_xzzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 920);

    auto g_0_xzzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 921);

    auto g_0_xzzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 922);

    auto g_0_xzzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 923);

    auto g_0_xzzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 924);

    auto g_0_xzzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 925);

    auto g_0_xzzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 926);

    auto g_0_xzzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 927);

    auto g_0_xzzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 928);

    auto g_0_xzzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 929);

    auto g_0_xzzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 930);

    auto g_0_xzzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 931);

    auto g_0_xzzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 932);

    auto g_0_xzzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 933);

    auto g_0_xzzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 934);

    auto g_0_xzzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 935);

    auto g_0_xzzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 936);

    auto g_0_xzzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 937);

    auto g_0_xzzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 938);

    auto g_0_xzzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 939);

    auto g_0_xzzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 940);

    auto g_0_xzzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 941);

    auto g_0_xzzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 942);

    auto g_0_xzzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 943);

    auto g_0_xzzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 944);

#pragma omp simd aligned(g_0_xzzzzz_0_xxxxxxxx_0,     \
                             g_0_xzzzzz_0_xxxxxxxy_0, \
                             g_0_xzzzzz_0_xxxxxxxz_0, \
                             g_0_xzzzzz_0_xxxxxxyy_0, \
                             g_0_xzzzzz_0_xxxxxxyz_0, \
                             g_0_xzzzzz_0_xxxxxxzz_0, \
                             g_0_xzzzzz_0_xxxxxyyy_0, \
                             g_0_xzzzzz_0_xxxxxyyz_0, \
                             g_0_xzzzzz_0_xxxxxyzz_0, \
                             g_0_xzzzzz_0_xxxxxzzz_0, \
                             g_0_xzzzzz_0_xxxxyyyy_0, \
                             g_0_xzzzzz_0_xxxxyyyz_0, \
                             g_0_xzzzzz_0_xxxxyyzz_0, \
                             g_0_xzzzzz_0_xxxxyzzz_0, \
                             g_0_xzzzzz_0_xxxxzzzz_0, \
                             g_0_xzzzzz_0_xxxyyyyy_0, \
                             g_0_xzzzzz_0_xxxyyyyz_0, \
                             g_0_xzzzzz_0_xxxyyyzz_0, \
                             g_0_xzzzzz_0_xxxyyzzz_0, \
                             g_0_xzzzzz_0_xxxyzzzz_0, \
                             g_0_xzzzzz_0_xxxzzzzz_0, \
                             g_0_xzzzzz_0_xxyyyyyy_0, \
                             g_0_xzzzzz_0_xxyyyyyz_0, \
                             g_0_xzzzzz_0_xxyyyyzz_0, \
                             g_0_xzzzzz_0_xxyyyzzz_0, \
                             g_0_xzzzzz_0_xxyyzzzz_0, \
                             g_0_xzzzzz_0_xxyzzzzz_0, \
                             g_0_xzzzzz_0_xxzzzzzz_0, \
                             g_0_xzzzzz_0_xyyyyyyy_0, \
                             g_0_xzzzzz_0_xyyyyyyz_0, \
                             g_0_xzzzzz_0_xyyyyyzz_0, \
                             g_0_xzzzzz_0_xyyyyzzz_0, \
                             g_0_xzzzzz_0_xyyyzzzz_0, \
                             g_0_xzzzzz_0_xyyzzzzz_0, \
                             g_0_xzzzzz_0_xyzzzzzz_0, \
                             g_0_xzzzzz_0_xzzzzzzz_0, \
                             g_0_xzzzzz_0_yyyyyyyy_0, \
                             g_0_xzzzzz_0_yyyyyyyz_0, \
                             g_0_xzzzzz_0_yyyyyyzz_0, \
                             g_0_xzzzzz_0_yyyyyzzz_0, \
                             g_0_xzzzzz_0_yyyyzzzz_0, \
                             g_0_xzzzzz_0_yyyzzzzz_0, \
                             g_0_xzzzzz_0_yyzzzzzz_0, \
                             g_0_xzzzzz_0_yzzzzzzz_0, \
                             g_0_xzzzzz_0_zzzzzzzz_0, \
                             g_0_zzzzz_0_xxxxxxx_1,   \
                             g_0_zzzzz_0_xxxxxxxx_0,  \
                             g_0_zzzzz_0_xxxxxxxx_1,  \
                             g_0_zzzzz_0_xxxxxxxy_0,  \
                             g_0_zzzzz_0_xxxxxxxy_1,  \
                             g_0_zzzzz_0_xxxxxxxz_0,  \
                             g_0_zzzzz_0_xxxxxxxz_1,  \
                             g_0_zzzzz_0_xxxxxxy_1,   \
                             g_0_zzzzz_0_xxxxxxyy_0,  \
                             g_0_zzzzz_0_xxxxxxyy_1,  \
                             g_0_zzzzz_0_xxxxxxyz_0,  \
                             g_0_zzzzz_0_xxxxxxyz_1,  \
                             g_0_zzzzz_0_xxxxxxz_1,   \
                             g_0_zzzzz_0_xxxxxxzz_0,  \
                             g_0_zzzzz_0_xxxxxxzz_1,  \
                             g_0_zzzzz_0_xxxxxyy_1,   \
                             g_0_zzzzz_0_xxxxxyyy_0,  \
                             g_0_zzzzz_0_xxxxxyyy_1,  \
                             g_0_zzzzz_0_xxxxxyyz_0,  \
                             g_0_zzzzz_0_xxxxxyyz_1,  \
                             g_0_zzzzz_0_xxxxxyz_1,   \
                             g_0_zzzzz_0_xxxxxyzz_0,  \
                             g_0_zzzzz_0_xxxxxyzz_1,  \
                             g_0_zzzzz_0_xxxxxzz_1,   \
                             g_0_zzzzz_0_xxxxxzzz_0,  \
                             g_0_zzzzz_0_xxxxxzzz_1,  \
                             g_0_zzzzz_0_xxxxyyy_1,   \
                             g_0_zzzzz_0_xxxxyyyy_0,  \
                             g_0_zzzzz_0_xxxxyyyy_1,  \
                             g_0_zzzzz_0_xxxxyyyz_0,  \
                             g_0_zzzzz_0_xxxxyyyz_1,  \
                             g_0_zzzzz_0_xxxxyyz_1,   \
                             g_0_zzzzz_0_xxxxyyzz_0,  \
                             g_0_zzzzz_0_xxxxyyzz_1,  \
                             g_0_zzzzz_0_xxxxyzz_1,   \
                             g_0_zzzzz_0_xxxxyzzz_0,  \
                             g_0_zzzzz_0_xxxxyzzz_1,  \
                             g_0_zzzzz_0_xxxxzzz_1,   \
                             g_0_zzzzz_0_xxxxzzzz_0,  \
                             g_0_zzzzz_0_xxxxzzzz_1,  \
                             g_0_zzzzz_0_xxxyyyy_1,   \
                             g_0_zzzzz_0_xxxyyyyy_0,  \
                             g_0_zzzzz_0_xxxyyyyy_1,  \
                             g_0_zzzzz_0_xxxyyyyz_0,  \
                             g_0_zzzzz_0_xxxyyyyz_1,  \
                             g_0_zzzzz_0_xxxyyyz_1,   \
                             g_0_zzzzz_0_xxxyyyzz_0,  \
                             g_0_zzzzz_0_xxxyyyzz_1,  \
                             g_0_zzzzz_0_xxxyyzz_1,   \
                             g_0_zzzzz_0_xxxyyzzz_0,  \
                             g_0_zzzzz_0_xxxyyzzz_1,  \
                             g_0_zzzzz_0_xxxyzzz_1,   \
                             g_0_zzzzz_0_xxxyzzzz_0,  \
                             g_0_zzzzz_0_xxxyzzzz_1,  \
                             g_0_zzzzz_0_xxxzzzz_1,   \
                             g_0_zzzzz_0_xxxzzzzz_0,  \
                             g_0_zzzzz_0_xxxzzzzz_1,  \
                             g_0_zzzzz_0_xxyyyyy_1,   \
                             g_0_zzzzz_0_xxyyyyyy_0,  \
                             g_0_zzzzz_0_xxyyyyyy_1,  \
                             g_0_zzzzz_0_xxyyyyyz_0,  \
                             g_0_zzzzz_0_xxyyyyyz_1,  \
                             g_0_zzzzz_0_xxyyyyz_1,   \
                             g_0_zzzzz_0_xxyyyyzz_0,  \
                             g_0_zzzzz_0_xxyyyyzz_1,  \
                             g_0_zzzzz_0_xxyyyzz_1,   \
                             g_0_zzzzz_0_xxyyyzzz_0,  \
                             g_0_zzzzz_0_xxyyyzzz_1,  \
                             g_0_zzzzz_0_xxyyzzz_1,   \
                             g_0_zzzzz_0_xxyyzzzz_0,  \
                             g_0_zzzzz_0_xxyyzzzz_1,  \
                             g_0_zzzzz_0_xxyzzzz_1,   \
                             g_0_zzzzz_0_xxyzzzzz_0,  \
                             g_0_zzzzz_0_xxyzzzzz_1,  \
                             g_0_zzzzz_0_xxzzzzz_1,   \
                             g_0_zzzzz_0_xxzzzzzz_0,  \
                             g_0_zzzzz_0_xxzzzzzz_1,  \
                             g_0_zzzzz_0_xyyyyyy_1,   \
                             g_0_zzzzz_0_xyyyyyyy_0,  \
                             g_0_zzzzz_0_xyyyyyyy_1,  \
                             g_0_zzzzz_0_xyyyyyyz_0,  \
                             g_0_zzzzz_0_xyyyyyyz_1,  \
                             g_0_zzzzz_0_xyyyyyz_1,   \
                             g_0_zzzzz_0_xyyyyyzz_0,  \
                             g_0_zzzzz_0_xyyyyyzz_1,  \
                             g_0_zzzzz_0_xyyyyzz_1,   \
                             g_0_zzzzz_0_xyyyyzzz_0,  \
                             g_0_zzzzz_0_xyyyyzzz_1,  \
                             g_0_zzzzz_0_xyyyzzz_1,   \
                             g_0_zzzzz_0_xyyyzzzz_0,  \
                             g_0_zzzzz_0_xyyyzzzz_1,  \
                             g_0_zzzzz_0_xyyzzzz_1,   \
                             g_0_zzzzz_0_xyyzzzzz_0,  \
                             g_0_zzzzz_0_xyyzzzzz_1,  \
                             g_0_zzzzz_0_xyzzzzz_1,   \
                             g_0_zzzzz_0_xyzzzzzz_0,  \
                             g_0_zzzzz_0_xyzzzzzz_1,  \
                             g_0_zzzzz_0_xzzzzzz_1,   \
                             g_0_zzzzz_0_xzzzzzzz_0,  \
                             g_0_zzzzz_0_xzzzzzzz_1,  \
                             g_0_zzzzz_0_yyyyyyy_1,   \
                             g_0_zzzzz_0_yyyyyyyy_0,  \
                             g_0_zzzzz_0_yyyyyyyy_1,  \
                             g_0_zzzzz_0_yyyyyyyz_0,  \
                             g_0_zzzzz_0_yyyyyyyz_1,  \
                             g_0_zzzzz_0_yyyyyyz_1,   \
                             g_0_zzzzz_0_yyyyyyzz_0,  \
                             g_0_zzzzz_0_yyyyyyzz_1,  \
                             g_0_zzzzz_0_yyyyyzz_1,   \
                             g_0_zzzzz_0_yyyyyzzz_0,  \
                             g_0_zzzzz_0_yyyyyzzz_1,  \
                             g_0_zzzzz_0_yyyyzzz_1,   \
                             g_0_zzzzz_0_yyyyzzzz_0,  \
                             g_0_zzzzz_0_yyyyzzzz_1,  \
                             g_0_zzzzz_0_yyyzzzz_1,   \
                             g_0_zzzzz_0_yyyzzzzz_0,  \
                             g_0_zzzzz_0_yyyzzzzz_1,  \
                             g_0_zzzzz_0_yyzzzzz_1,   \
                             g_0_zzzzz_0_yyzzzzzz_0,  \
                             g_0_zzzzz_0_yyzzzzzz_1,  \
                             g_0_zzzzz_0_yzzzzzz_1,   \
                             g_0_zzzzz_0_yzzzzzzz_0,  \
                             g_0_zzzzz_0_yzzzzzzz_1,  \
                             g_0_zzzzz_0_zzzzzzz_1,   \
                             g_0_zzzzz_0_zzzzzzzz_0,  \
                             g_0_zzzzz_0_zzzzzzzz_1,  \
                             wp_x,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzz_0_xxxxxxxx_0[i] =
            8.0 * g_0_zzzzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxxx_0[i] * pb_x + g_0_zzzzz_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxxxy_0[i] =
            7.0 * g_0_zzzzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxxy_0[i] * pb_x + g_0_zzzzz_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxxxz_0[i] =
            7.0 * g_0_zzzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxxz_0[i] * pb_x + g_0_zzzzz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxxyy_0[i] =
            6.0 * g_0_zzzzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxyy_0[i] * pb_x + g_0_zzzzz_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxxyz_0[i] =
            6.0 * g_0_zzzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxyz_0[i] * pb_x + g_0_zzzzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxxzz_0[i] =
            6.0 * g_0_zzzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxzz_0[i] * pb_x + g_0_zzzzz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxyyy_0[i] =
            5.0 * g_0_zzzzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyyy_0[i] * pb_x + g_0_zzzzz_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxyyz_0[i] =
            5.0 * g_0_zzzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyyz_0[i] * pb_x + g_0_zzzzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxyzz_0[i] =
            5.0 * g_0_zzzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyzz_0[i] * pb_x + g_0_zzzzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxzzz_0[i] =
            5.0 * g_0_zzzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxzzz_0[i] * pb_x + g_0_zzzzz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxyyyy_0[i] =
            4.0 * g_0_zzzzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyyy_0[i] * pb_x + g_0_zzzzz_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxyyyz_0[i] =
            4.0 * g_0_zzzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyyz_0[i] * pb_x + g_0_zzzzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxyyzz_0[i] =
            4.0 * g_0_zzzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyzz_0[i] * pb_x + g_0_zzzzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxyzzz_0[i] =
            4.0 * g_0_zzzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyzzz_0[i] * pb_x + g_0_zzzzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxzzzz_0[i] =
            4.0 * g_0_zzzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxzzzz_0[i] * pb_x + g_0_zzzzz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyyyyy_0[i] =
            3.0 * g_0_zzzzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyyy_0[i] * pb_x + g_0_zzzzz_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyyyyz_0[i] =
            3.0 * g_0_zzzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyyz_0[i] * pb_x + g_0_zzzzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyyyzz_0[i] =
            3.0 * g_0_zzzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyzz_0[i] * pb_x + g_0_zzzzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyyzzz_0[i] =
            3.0 * g_0_zzzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyzzz_0[i] * pb_x + g_0_zzzzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyzzzz_0[i] =
            3.0 * g_0_zzzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyzzzz_0[i] * pb_x + g_0_zzzzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxzzzzz_0[i] =
            3.0 * g_0_zzzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxzzzzz_0[i] * pb_x + g_0_zzzzz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyyyyy_0[i] =
            2.0 * g_0_zzzzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyyy_0[i] * pb_x + g_0_zzzzz_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyyyyz_0[i] =
            2.0 * g_0_zzzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyyz_0[i] * pb_x + g_0_zzzzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyyyzz_0[i] =
            2.0 * g_0_zzzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyzz_0[i] * pb_x + g_0_zzzzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyyzzz_0[i] =
            2.0 * g_0_zzzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyzzz_0[i] * pb_x + g_0_zzzzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyzzzz_0[i] =
            2.0 * g_0_zzzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyzzzz_0[i] * pb_x + g_0_zzzzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyzzzzz_0[i] =
            2.0 * g_0_zzzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyzzzzz_0[i] * pb_x + g_0_zzzzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxzzzzzz_0[i] =
            2.0 * g_0_zzzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxzzzzzz_0[i] * pb_x + g_0_zzzzz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyyyyy_0[i] = g_0_zzzzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyyy_0[i] * pb_x + g_0_zzzzz_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyyyyz_0[i] = g_0_zzzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyyz_0[i] * pb_x + g_0_zzzzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyyyzz_0[i] = g_0_zzzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyzz_0[i] * pb_x + g_0_zzzzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyyzzz_0[i] = g_0_zzzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyzzz_0[i] * pb_x + g_0_zzzzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyzzzz_0[i] = g_0_zzzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyzzzz_0[i] * pb_x + g_0_zzzzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyzzzzz_0[i] = g_0_zzzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyzzzzz_0[i] * pb_x + g_0_zzzzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyzzzzzz_0[i] = g_0_zzzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzzzzzz_0[i] * pb_x + g_0_zzzzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xzzzzzzz_0[i] = g_0_zzzzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xzzzzzzz_0[i] * pb_x + g_0_zzzzz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyyyyy_0[i] = g_0_zzzzz_0_yyyyyyyy_0[i] * pb_x + g_0_zzzzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyyyyz_0[i] = g_0_zzzzz_0_yyyyyyyz_0[i] * pb_x + g_0_zzzzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyyyzz_0[i] = g_0_zzzzz_0_yyyyyyzz_0[i] * pb_x + g_0_zzzzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyyzzz_0[i] = g_0_zzzzz_0_yyyyyzzz_0[i] * pb_x + g_0_zzzzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyzzzz_0[i] = g_0_zzzzz_0_yyyyzzzz_0[i] * pb_x + g_0_zzzzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyzzzzz_0[i] = g_0_zzzzz_0_yyyzzzzz_0[i] * pb_x + g_0_zzzzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyzzzzzz_0[i] = g_0_zzzzz_0_yyzzzzzz_0[i] * pb_x + g_0_zzzzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yzzzzzzz_0[i] = g_0_zzzzz_0_yzzzzzzz_0[i] * pb_x + g_0_zzzzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_zzzzzzzz_0[i] = g_0_zzzzz_0_zzzzzzzz_0[i] * pb_x + g_0_zzzzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 945-990 components of targeted buffer : SISL

    auto g_0_yyyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 945);

    auto g_0_yyyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 946);

    auto g_0_yyyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 947);

    auto g_0_yyyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 948);

    auto g_0_yyyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 949);

    auto g_0_yyyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 950);

    auto g_0_yyyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 951);

    auto g_0_yyyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 952);

    auto g_0_yyyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 953);

    auto g_0_yyyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 954);

    auto g_0_yyyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 955);

    auto g_0_yyyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 956);

    auto g_0_yyyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 957);

    auto g_0_yyyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 958);

    auto g_0_yyyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 959);

    auto g_0_yyyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 960);

    auto g_0_yyyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 961);

    auto g_0_yyyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 962);

    auto g_0_yyyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 963);

    auto g_0_yyyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 964);

    auto g_0_yyyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 965);

    auto g_0_yyyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 966);

    auto g_0_yyyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 967);

    auto g_0_yyyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 968);

    auto g_0_yyyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 969);

    auto g_0_yyyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 970);

    auto g_0_yyyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 971);

    auto g_0_yyyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 972);

    auto g_0_yyyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 973);

    auto g_0_yyyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 974);

    auto g_0_yyyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 975);

    auto g_0_yyyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 976);

    auto g_0_yyyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 977);

    auto g_0_yyyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 978);

    auto g_0_yyyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 979);

    auto g_0_yyyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 980);

    auto g_0_yyyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 981);

    auto g_0_yyyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 982);

    auto g_0_yyyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 983);

    auto g_0_yyyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 984);

    auto g_0_yyyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 985);

    auto g_0_yyyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 986);

    auto g_0_yyyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 987);

    auto g_0_yyyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 988);

    auto g_0_yyyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 989);

#pragma omp simd aligned(g_0_yyyy_0_xxxxxxxx_0,       \
                             g_0_yyyy_0_xxxxxxxx_1,   \
                             g_0_yyyy_0_xxxxxxxy_0,   \
                             g_0_yyyy_0_xxxxxxxy_1,   \
                             g_0_yyyy_0_xxxxxxxz_0,   \
                             g_0_yyyy_0_xxxxxxxz_1,   \
                             g_0_yyyy_0_xxxxxxyy_0,   \
                             g_0_yyyy_0_xxxxxxyy_1,   \
                             g_0_yyyy_0_xxxxxxyz_0,   \
                             g_0_yyyy_0_xxxxxxyz_1,   \
                             g_0_yyyy_0_xxxxxxzz_0,   \
                             g_0_yyyy_0_xxxxxxzz_1,   \
                             g_0_yyyy_0_xxxxxyyy_0,   \
                             g_0_yyyy_0_xxxxxyyy_1,   \
                             g_0_yyyy_0_xxxxxyyz_0,   \
                             g_0_yyyy_0_xxxxxyyz_1,   \
                             g_0_yyyy_0_xxxxxyzz_0,   \
                             g_0_yyyy_0_xxxxxyzz_1,   \
                             g_0_yyyy_0_xxxxxzzz_0,   \
                             g_0_yyyy_0_xxxxxzzz_1,   \
                             g_0_yyyy_0_xxxxyyyy_0,   \
                             g_0_yyyy_0_xxxxyyyy_1,   \
                             g_0_yyyy_0_xxxxyyyz_0,   \
                             g_0_yyyy_0_xxxxyyyz_1,   \
                             g_0_yyyy_0_xxxxyyzz_0,   \
                             g_0_yyyy_0_xxxxyyzz_1,   \
                             g_0_yyyy_0_xxxxyzzz_0,   \
                             g_0_yyyy_0_xxxxyzzz_1,   \
                             g_0_yyyy_0_xxxxzzzz_0,   \
                             g_0_yyyy_0_xxxxzzzz_1,   \
                             g_0_yyyy_0_xxxyyyyy_0,   \
                             g_0_yyyy_0_xxxyyyyy_1,   \
                             g_0_yyyy_0_xxxyyyyz_0,   \
                             g_0_yyyy_0_xxxyyyyz_1,   \
                             g_0_yyyy_0_xxxyyyzz_0,   \
                             g_0_yyyy_0_xxxyyyzz_1,   \
                             g_0_yyyy_0_xxxyyzzz_0,   \
                             g_0_yyyy_0_xxxyyzzz_1,   \
                             g_0_yyyy_0_xxxyzzzz_0,   \
                             g_0_yyyy_0_xxxyzzzz_1,   \
                             g_0_yyyy_0_xxxzzzzz_0,   \
                             g_0_yyyy_0_xxxzzzzz_1,   \
                             g_0_yyyy_0_xxyyyyyy_0,   \
                             g_0_yyyy_0_xxyyyyyy_1,   \
                             g_0_yyyy_0_xxyyyyyz_0,   \
                             g_0_yyyy_0_xxyyyyyz_1,   \
                             g_0_yyyy_0_xxyyyyzz_0,   \
                             g_0_yyyy_0_xxyyyyzz_1,   \
                             g_0_yyyy_0_xxyyyzzz_0,   \
                             g_0_yyyy_0_xxyyyzzz_1,   \
                             g_0_yyyy_0_xxyyzzzz_0,   \
                             g_0_yyyy_0_xxyyzzzz_1,   \
                             g_0_yyyy_0_xxyzzzzz_0,   \
                             g_0_yyyy_0_xxyzzzzz_1,   \
                             g_0_yyyy_0_xxzzzzzz_0,   \
                             g_0_yyyy_0_xxzzzzzz_1,   \
                             g_0_yyyy_0_xyyyyyyy_0,   \
                             g_0_yyyy_0_xyyyyyyy_1,   \
                             g_0_yyyy_0_xyyyyyyz_0,   \
                             g_0_yyyy_0_xyyyyyyz_1,   \
                             g_0_yyyy_0_xyyyyyzz_0,   \
                             g_0_yyyy_0_xyyyyyzz_1,   \
                             g_0_yyyy_0_xyyyyzzz_0,   \
                             g_0_yyyy_0_xyyyyzzz_1,   \
                             g_0_yyyy_0_xyyyzzzz_0,   \
                             g_0_yyyy_0_xyyyzzzz_1,   \
                             g_0_yyyy_0_xyyzzzzz_0,   \
                             g_0_yyyy_0_xyyzzzzz_1,   \
                             g_0_yyyy_0_xyzzzzzz_0,   \
                             g_0_yyyy_0_xyzzzzzz_1,   \
                             g_0_yyyy_0_xzzzzzzz_0,   \
                             g_0_yyyy_0_xzzzzzzz_1,   \
                             g_0_yyyy_0_yyyyyyyy_0,   \
                             g_0_yyyy_0_yyyyyyyy_1,   \
                             g_0_yyyy_0_yyyyyyyz_0,   \
                             g_0_yyyy_0_yyyyyyyz_1,   \
                             g_0_yyyy_0_yyyyyyzz_0,   \
                             g_0_yyyy_0_yyyyyyzz_1,   \
                             g_0_yyyy_0_yyyyyzzz_0,   \
                             g_0_yyyy_0_yyyyyzzz_1,   \
                             g_0_yyyy_0_yyyyzzzz_0,   \
                             g_0_yyyy_0_yyyyzzzz_1,   \
                             g_0_yyyy_0_yyyzzzzz_0,   \
                             g_0_yyyy_0_yyyzzzzz_1,   \
                             g_0_yyyy_0_yyzzzzzz_0,   \
                             g_0_yyyy_0_yyzzzzzz_1,   \
                             g_0_yyyy_0_yzzzzzzz_0,   \
                             g_0_yyyy_0_yzzzzzzz_1,   \
                             g_0_yyyy_0_zzzzzzzz_0,   \
                             g_0_yyyy_0_zzzzzzzz_1,   \
                             g_0_yyyyy_0_xxxxxxx_1,   \
                             g_0_yyyyy_0_xxxxxxxx_0,  \
                             g_0_yyyyy_0_xxxxxxxx_1,  \
                             g_0_yyyyy_0_xxxxxxxy_0,  \
                             g_0_yyyyy_0_xxxxxxxy_1,  \
                             g_0_yyyyy_0_xxxxxxxz_0,  \
                             g_0_yyyyy_0_xxxxxxxz_1,  \
                             g_0_yyyyy_0_xxxxxxy_1,   \
                             g_0_yyyyy_0_xxxxxxyy_0,  \
                             g_0_yyyyy_0_xxxxxxyy_1,  \
                             g_0_yyyyy_0_xxxxxxyz_0,  \
                             g_0_yyyyy_0_xxxxxxyz_1,  \
                             g_0_yyyyy_0_xxxxxxz_1,   \
                             g_0_yyyyy_0_xxxxxxzz_0,  \
                             g_0_yyyyy_0_xxxxxxzz_1,  \
                             g_0_yyyyy_0_xxxxxyy_1,   \
                             g_0_yyyyy_0_xxxxxyyy_0,  \
                             g_0_yyyyy_0_xxxxxyyy_1,  \
                             g_0_yyyyy_0_xxxxxyyz_0,  \
                             g_0_yyyyy_0_xxxxxyyz_1,  \
                             g_0_yyyyy_0_xxxxxyz_1,   \
                             g_0_yyyyy_0_xxxxxyzz_0,  \
                             g_0_yyyyy_0_xxxxxyzz_1,  \
                             g_0_yyyyy_0_xxxxxzz_1,   \
                             g_0_yyyyy_0_xxxxxzzz_0,  \
                             g_0_yyyyy_0_xxxxxzzz_1,  \
                             g_0_yyyyy_0_xxxxyyy_1,   \
                             g_0_yyyyy_0_xxxxyyyy_0,  \
                             g_0_yyyyy_0_xxxxyyyy_1,  \
                             g_0_yyyyy_0_xxxxyyyz_0,  \
                             g_0_yyyyy_0_xxxxyyyz_1,  \
                             g_0_yyyyy_0_xxxxyyz_1,   \
                             g_0_yyyyy_0_xxxxyyzz_0,  \
                             g_0_yyyyy_0_xxxxyyzz_1,  \
                             g_0_yyyyy_0_xxxxyzz_1,   \
                             g_0_yyyyy_0_xxxxyzzz_0,  \
                             g_0_yyyyy_0_xxxxyzzz_1,  \
                             g_0_yyyyy_0_xxxxzzz_1,   \
                             g_0_yyyyy_0_xxxxzzzz_0,  \
                             g_0_yyyyy_0_xxxxzzzz_1,  \
                             g_0_yyyyy_0_xxxyyyy_1,   \
                             g_0_yyyyy_0_xxxyyyyy_0,  \
                             g_0_yyyyy_0_xxxyyyyy_1,  \
                             g_0_yyyyy_0_xxxyyyyz_0,  \
                             g_0_yyyyy_0_xxxyyyyz_1,  \
                             g_0_yyyyy_0_xxxyyyz_1,   \
                             g_0_yyyyy_0_xxxyyyzz_0,  \
                             g_0_yyyyy_0_xxxyyyzz_1,  \
                             g_0_yyyyy_0_xxxyyzz_1,   \
                             g_0_yyyyy_0_xxxyyzzz_0,  \
                             g_0_yyyyy_0_xxxyyzzz_1,  \
                             g_0_yyyyy_0_xxxyzzz_1,   \
                             g_0_yyyyy_0_xxxyzzzz_0,  \
                             g_0_yyyyy_0_xxxyzzzz_1,  \
                             g_0_yyyyy_0_xxxzzzz_1,   \
                             g_0_yyyyy_0_xxxzzzzz_0,  \
                             g_0_yyyyy_0_xxxzzzzz_1,  \
                             g_0_yyyyy_0_xxyyyyy_1,   \
                             g_0_yyyyy_0_xxyyyyyy_0,  \
                             g_0_yyyyy_0_xxyyyyyy_1,  \
                             g_0_yyyyy_0_xxyyyyyz_0,  \
                             g_0_yyyyy_0_xxyyyyyz_1,  \
                             g_0_yyyyy_0_xxyyyyz_1,   \
                             g_0_yyyyy_0_xxyyyyzz_0,  \
                             g_0_yyyyy_0_xxyyyyzz_1,  \
                             g_0_yyyyy_0_xxyyyzz_1,   \
                             g_0_yyyyy_0_xxyyyzzz_0,  \
                             g_0_yyyyy_0_xxyyyzzz_1,  \
                             g_0_yyyyy_0_xxyyzzz_1,   \
                             g_0_yyyyy_0_xxyyzzzz_0,  \
                             g_0_yyyyy_0_xxyyzzzz_1,  \
                             g_0_yyyyy_0_xxyzzzz_1,   \
                             g_0_yyyyy_0_xxyzzzzz_0,  \
                             g_0_yyyyy_0_xxyzzzzz_1,  \
                             g_0_yyyyy_0_xxzzzzz_1,   \
                             g_0_yyyyy_0_xxzzzzzz_0,  \
                             g_0_yyyyy_0_xxzzzzzz_1,  \
                             g_0_yyyyy_0_xyyyyyy_1,   \
                             g_0_yyyyy_0_xyyyyyyy_0,  \
                             g_0_yyyyy_0_xyyyyyyy_1,  \
                             g_0_yyyyy_0_xyyyyyyz_0,  \
                             g_0_yyyyy_0_xyyyyyyz_1,  \
                             g_0_yyyyy_0_xyyyyyz_1,   \
                             g_0_yyyyy_0_xyyyyyzz_0,  \
                             g_0_yyyyy_0_xyyyyyzz_1,  \
                             g_0_yyyyy_0_xyyyyzz_1,   \
                             g_0_yyyyy_0_xyyyyzzz_0,  \
                             g_0_yyyyy_0_xyyyyzzz_1,  \
                             g_0_yyyyy_0_xyyyzzz_1,   \
                             g_0_yyyyy_0_xyyyzzzz_0,  \
                             g_0_yyyyy_0_xyyyzzzz_1,  \
                             g_0_yyyyy_0_xyyzzzz_1,   \
                             g_0_yyyyy_0_xyyzzzzz_0,  \
                             g_0_yyyyy_0_xyyzzzzz_1,  \
                             g_0_yyyyy_0_xyzzzzz_1,   \
                             g_0_yyyyy_0_xyzzzzzz_0,  \
                             g_0_yyyyy_0_xyzzzzzz_1,  \
                             g_0_yyyyy_0_xzzzzzz_1,   \
                             g_0_yyyyy_0_xzzzzzzz_0,  \
                             g_0_yyyyy_0_xzzzzzzz_1,  \
                             g_0_yyyyy_0_yyyyyyy_1,   \
                             g_0_yyyyy_0_yyyyyyyy_0,  \
                             g_0_yyyyy_0_yyyyyyyy_1,  \
                             g_0_yyyyy_0_yyyyyyyz_0,  \
                             g_0_yyyyy_0_yyyyyyyz_1,  \
                             g_0_yyyyy_0_yyyyyyz_1,   \
                             g_0_yyyyy_0_yyyyyyzz_0,  \
                             g_0_yyyyy_0_yyyyyyzz_1,  \
                             g_0_yyyyy_0_yyyyyzz_1,   \
                             g_0_yyyyy_0_yyyyyzzz_0,  \
                             g_0_yyyyy_0_yyyyyzzz_1,  \
                             g_0_yyyyy_0_yyyyzzz_1,   \
                             g_0_yyyyy_0_yyyyzzzz_0,  \
                             g_0_yyyyy_0_yyyyzzzz_1,  \
                             g_0_yyyyy_0_yyyzzzz_1,   \
                             g_0_yyyyy_0_yyyzzzzz_0,  \
                             g_0_yyyyy_0_yyyzzzzz_1,  \
                             g_0_yyyyy_0_yyzzzzz_1,   \
                             g_0_yyyyy_0_yyzzzzzz_0,  \
                             g_0_yyyyy_0_yyzzzzzz_1,  \
                             g_0_yyyyy_0_yzzzzzz_1,   \
                             g_0_yyyyy_0_yzzzzzzz_0,  \
                             g_0_yyyyy_0_yzzzzzzz_1,  \
                             g_0_yyyyy_0_zzzzzzz_1,   \
                             g_0_yyyyy_0_zzzzzzzz_0,  \
                             g_0_yyyyy_0_zzzzzzzz_1,  \
                             g_0_yyyyyy_0_xxxxxxxx_0, \
                             g_0_yyyyyy_0_xxxxxxxy_0, \
                             g_0_yyyyyy_0_xxxxxxxz_0, \
                             g_0_yyyyyy_0_xxxxxxyy_0, \
                             g_0_yyyyyy_0_xxxxxxyz_0, \
                             g_0_yyyyyy_0_xxxxxxzz_0, \
                             g_0_yyyyyy_0_xxxxxyyy_0, \
                             g_0_yyyyyy_0_xxxxxyyz_0, \
                             g_0_yyyyyy_0_xxxxxyzz_0, \
                             g_0_yyyyyy_0_xxxxxzzz_0, \
                             g_0_yyyyyy_0_xxxxyyyy_0, \
                             g_0_yyyyyy_0_xxxxyyyz_0, \
                             g_0_yyyyyy_0_xxxxyyzz_0, \
                             g_0_yyyyyy_0_xxxxyzzz_0, \
                             g_0_yyyyyy_0_xxxxzzzz_0, \
                             g_0_yyyyyy_0_xxxyyyyy_0, \
                             g_0_yyyyyy_0_xxxyyyyz_0, \
                             g_0_yyyyyy_0_xxxyyyzz_0, \
                             g_0_yyyyyy_0_xxxyyzzz_0, \
                             g_0_yyyyyy_0_xxxyzzzz_0, \
                             g_0_yyyyyy_0_xxxzzzzz_0, \
                             g_0_yyyyyy_0_xxyyyyyy_0, \
                             g_0_yyyyyy_0_xxyyyyyz_0, \
                             g_0_yyyyyy_0_xxyyyyzz_0, \
                             g_0_yyyyyy_0_xxyyyzzz_0, \
                             g_0_yyyyyy_0_xxyyzzzz_0, \
                             g_0_yyyyyy_0_xxyzzzzz_0, \
                             g_0_yyyyyy_0_xxzzzzzz_0, \
                             g_0_yyyyyy_0_xyyyyyyy_0, \
                             g_0_yyyyyy_0_xyyyyyyz_0, \
                             g_0_yyyyyy_0_xyyyyyzz_0, \
                             g_0_yyyyyy_0_xyyyyzzz_0, \
                             g_0_yyyyyy_0_xyyyzzzz_0, \
                             g_0_yyyyyy_0_xyyzzzzz_0, \
                             g_0_yyyyyy_0_xyzzzzzz_0, \
                             g_0_yyyyyy_0_xzzzzzzz_0, \
                             g_0_yyyyyy_0_yyyyyyyy_0, \
                             g_0_yyyyyy_0_yyyyyyyz_0, \
                             g_0_yyyyyy_0_yyyyyyzz_0, \
                             g_0_yyyyyy_0_yyyyyzzz_0, \
                             g_0_yyyyyy_0_yyyyzzzz_0, \
                             g_0_yyyyyy_0_yyyzzzzz_0, \
                             g_0_yyyyyy_0_yyzzzzzz_0, \
                             g_0_yyyyyy_0_yzzzzzzz_0, \
                             g_0_yyyyyy_0_zzzzzzzz_0, \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyy_0_xxxxxxxx_0[i] = 5.0 * g_0_yyyy_0_xxxxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxxxx_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxxxxxxx_0[i] * pb_y + g_0_yyyyy_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxxxy_0[i] = 5.0 * g_0_yyyy_0_xxxxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxxxy_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxxy_0[i] * pb_y + g_0_yyyyy_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxxxz_0[i] = 5.0 * g_0_yyyy_0_xxxxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxxxz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxxxxxxz_0[i] * pb_y + g_0_yyyyy_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxxyy_0[i] = 5.0 * g_0_yyyy_0_xxxxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxxyy_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxyy_0[i] * pb_y +
                                     g_0_yyyyy_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxxyz_0[i] = 5.0 * g_0_yyyy_0_xxxxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxxyz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxyz_0[i] * pb_y + g_0_yyyyy_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxxzz_0[i] = 5.0 * g_0_yyyy_0_xxxxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxxzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxxxxxzz_0[i] * pb_y + g_0_yyyyy_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxyyy_0[i] = 5.0 * g_0_yyyy_0_xxxxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxyyy_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyyy_0[i] * pb_y +
                                     g_0_yyyyy_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxyyz_0[i] = 5.0 * g_0_yyyy_0_xxxxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyyz_0[i] * pb_y +
                                     g_0_yyyyy_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxyzz_0[i] = 5.0 * g_0_yyyy_0_xxxxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxyzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyzz_0[i] * pb_y + g_0_yyyyy_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxzzz_0[i] = 5.0 * g_0_yyyy_0_xxxxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxxxxzzz_0[i] * pb_y + g_0_yyyyy_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxyyyy_0[i] = 5.0 * g_0_yyyy_0_xxxxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxyyyy_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyyy_0[i] * pb_y +
                                     g_0_yyyyy_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxyyyz_0[i] = 5.0 * g_0_yyyy_0_xxxxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyyz_0[i] * pb_y +
                                     g_0_yyyyy_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxyyzz_0[i] = 5.0 * g_0_yyyy_0_xxxxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyzz_0[i] * pb_y +
                                     g_0_yyyyy_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxyzzz_0[i] = 5.0 * g_0_yyyy_0_xxxxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxyzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyzzz_0[i] * pb_y + g_0_yyyyy_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxzzzz_0[i] = 5.0 * g_0_yyyy_0_xxxxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxxxzzzz_0[i] * pb_y + g_0_yyyyy_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyyyyy_0[i] = 5.0 * g_0_yyyy_0_xxxyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyyyyy_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyyyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyyy_0[i] * pb_y +
                                     g_0_yyyyy_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyyyyz_0[i] = 5.0 * g_0_yyyy_0_xxxyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyyz_0[i] * pb_y +
                                     g_0_yyyyy_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyyyzz_0[i] = 5.0 * g_0_yyyy_0_xxxyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyzz_0[i] * pb_y +
                                     g_0_yyyyy_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyyzzz_0[i] = 5.0 * g_0_yyyy_0_xxxyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyzzz_0[i] * pb_y +
                                     g_0_yyyyy_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyzzzz_0[i] = 5.0 * g_0_yyyy_0_xxxyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyzzzz_0[i] * pb_y + g_0_yyyyy_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxzzzzz_0[i] = 5.0 * g_0_yyyy_0_xxxzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxxzzzzz_0[i] * pb_y + g_0_yyyyy_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyyyyy_0[i] = 5.0 * g_0_yyyy_0_xxyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyyyyy_1[i] * fti_ab_0 +
                                     6.0 * g_0_yyyyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyyy_0[i] * pb_y +
                                     g_0_yyyyy_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyyyyz_0[i] = 5.0 * g_0_yyyy_0_xxyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyyyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyyz_0[i] * pb_y +
                                     g_0_yyyyy_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyyyzz_0[i] = 5.0 * g_0_yyyy_0_xxyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyzz_0[i] * pb_y +
                                     g_0_yyyyy_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyyzzz_0[i] = 5.0 * g_0_yyyy_0_xxyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyzzz_0[i] * pb_y +
                                     g_0_yyyyy_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyzzzz_0[i] = 5.0 * g_0_yyyy_0_xxyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyzzzz_0[i] * pb_y +
                                     g_0_yyyyy_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyzzzzz_0[i] = 5.0 * g_0_yyyy_0_xxyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyzzzzz_0[i] * pb_y + g_0_yyyyy_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxzzzzzz_0[i] = 5.0 * g_0_yyyy_0_xxzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xxzzzzzz_0[i] * pb_y + g_0_yyyyy_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyyyyy_0[i] = 5.0 * g_0_yyyy_0_xyyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyyyyy_1[i] * fti_ab_0 +
                                     7.0 * g_0_yyyyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyyy_0[i] * pb_y +
                                     g_0_yyyyy_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyyyyz_0[i] = 5.0 * g_0_yyyy_0_xyyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyyyyz_1[i] * fti_ab_0 +
                                     6.0 * g_0_yyyyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyyz_0[i] * pb_y +
                                     g_0_yyyyy_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyyyzz_0[i] = 5.0 * g_0_yyyy_0_xyyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyyyzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyyyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyzz_0[i] * pb_y +
                                     g_0_yyyyy_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyyzzz_0[i] = 5.0 * g_0_yyyy_0_xyyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyyzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyzzz_0[i] * pb_y +
                                     g_0_yyyyy_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyzzzz_0[i] = 5.0 * g_0_yyyy_0_xyyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyzzzz_0[i] * pb_y +
                                     g_0_yyyyy_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyzzzzz_0[i] = 5.0 * g_0_yyyy_0_xyyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyzzzzz_0[i] * pb_y +
                                     g_0_yyyyy_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyzzzzzz_0[i] = 5.0 * g_0_yyyy_0_xyzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzzzzzz_0[i] * pb_y + g_0_yyyyy_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xzzzzzzz_0[i] = 5.0 * g_0_yyyy_0_xzzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_xzzzzzzz_0[i] * pb_y + g_0_yyyyy_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyyyyy_0[i] = 5.0 * g_0_yyyy_0_yyyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyyyyy_1[i] * fti_ab_0 +
                                     8.0 * g_0_yyyyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyyyy_0[i] * pb_y +
                                     g_0_yyyyy_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyyyyz_0[i] = 5.0 * g_0_yyyy_0_yyyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyyyyz_1[i] * fti_ab_0 +
                                     7.0 * g_0_yyyyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyyyz_0[i] * pb_y +
                                     g_0_yyyyy_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyyyzz_0[i] = 5.0 * g_0_yyyy_0_yyyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyyyzz_1[i] * fti_ab_0 +
                                     6.0 * g_0_yyyyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyyzz_0[i] * pb_y +
                                     g_0_yyyyy_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyyzzz_0[i] = 5.0 * g_0_yyyy_0_yyyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyyzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyyyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyzzz_0[i] * pb_y +
                                     g_0_yyyyy_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyzzzz_0[i] = 5.0 * g_0_yyyy_0_yyyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyzzzz_0[i] * pb_y +
                                     g_0_yyyyy_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyzzzzz_0[i] = 5.0 * g_0_yyyy_0_yyyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyzzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyzzzzz_0[i] * pb_y +
                                     g_0_yyyyy_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyzzzzzz_0[i] = 5.0 * g_0_yyyy_0_yyzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyzzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyzzzzzz_0[i] * pb_y +
                                     g_0_yyyyy_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yzzzzzzz_0[i] = 5.0 * g_0_yyyy_0_yzzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yzzzzzzz_0[i] * pb_y + g_0_yyyyy_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_zzzzzzzz_0[i] = 5.0 * g_0_yyyy_0_zzzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_zzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyy_0_zzzzzzzz_0[i] * pb_y + g_0_yyyyy_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 990-1035 components of targeted buffer : SISL

    auto g_0_yyyyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 990);

    auto g_0_yyyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 991);

    auto g_0_yyyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 992);

    auto g_0_yyyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 993);

    auto g_0_yyyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 994);

    auto g_0_yyyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 995);

    auto g_0_yyyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 996);

    auto g_0_yyyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 997);

    auto g_0_yyyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 998);

    auto g_0_yyyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 999);

    auto g_0_yyyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1000);

    auto g_0_yyyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1001);

    auto g_0_yyyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1002);

    auto g_0_yyyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1003);

    auto g_0_yyyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1004);

    auto g_0_yyyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1005);

    auto g_0_yyyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1006);

    auto g_0_yyyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1007);

    auto g_0_yyyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1008);

    auto g_0_yyyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1009);

    auto g_0_yyyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1010);

    auto g_0_yyyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1011);

    auto g_0_yyyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1012);

    auto g_0_yyyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1013);

    auto g_0_yyyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1014);

    auto g_0_yyyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1015);

    auto g_0_yyyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1016);

    auto g_0_yyyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1017);

    auto g_0_yyyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1018);

    auto g_0_yyyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1019);

    auto g_0_yyyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1020);

    auto g_0_yyyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1021);

    auto g_0_yyyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1022);

    auto g_0_yyyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1023);

    auto g_0_yyyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1024);

    auto g_0_yyyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1025);

    auto g_0_yyyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1026);

    auto g_0_yyyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1027);

    auto g_0_yyyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1028);

    auto g_0_yyyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1029);

    auto g_0_yyyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1030);

    auto g_0_yyyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1031);

    auto g_0_yyyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1032);

    auto g_0_yyyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1033);

    auto g_0_yyyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1034);

#pragma omp simd aligned(g_0_yyyyy_0_xxxxxxx_1,       \
                             g_0_yyyyy_0_xxxxxxxx_0,  \
                             g_0_yyyyy_0_xxxxxxxx_1,  \
                             g_0_yyyyy_0_xxxxxxxy_0,  \
                             g_0_yyyyy_0_xxxxxxxy_1,  \
                             g_0_yyyyy_0_xxxxxxxz_0,  \
                             g_0_yyyyy_0_xxxxxxxz_1,  \
                             g_0_yyyyy_0_xxxxxxy_1,   \
                             g_0_yyyyy_0_xxxxxxyy_0,  \
                             g_0_yyyyy_0_xxxxxxyy_1,  \
                             g_0_yyyyy_0_xxxxxxyz_0,  \
                             g_0_yyyyy_0_xxxxxxyz_1,  \
                             g_0_yyyyy_0_xxxxxxz_1,   \
                             g_0_yyyyy_0_xxxxxxzz_0,  \
                             g_0_yyyyy_0_xxxxxxzz_1,  \
                             g_0_yyyyy_0_xxxxxyy_1,   \
                             g_0_yyyyy_0_xxxxxyyy_0,  \
                             g_0_yyyyy_0_xxxxxyyy_1,  \
                             g_0_yyyyy_0_xxxxxyyz_0,  \
                             g_0_yyyyy_0_xxxxxyyz_1,  \
                             g_0_yyyyy_0_xxxxxyz_1,   \
                             g_0_yyyyy_0_xxxxxyzz_0,  \
                             g_0_yyyyy_0_xxxxxyzz_1,  \
                             g_0_yyyyy_0_xxxxxzz_1,   \
                             g_0_yyyyy_0_xxxxxzzz_0,  \
                             g_0_yyyyy_0_xxxxxzzz_1,  \
                             g_0_yyyyy_0_xxxxyyy_1,   \
                             g_0_yyyyy_0_xxxxyyyy_0,  \
                             g_0_yyyyy_0_xxxxyyyy_1,  \
                             g_0_yyyyy_0_xxxxyyyz_0,  \
                             g_0_yyyyy_0_xxxxyyyz_1,  \
                             g_0_yyyyy_0_xxxxyyz_1,   \
                             g_0_yyyyy_0_xxxxyyzz_0,  \
                             g_0_yyyyy_0_xxxxyyzz_1,  \
                             g_0_yyyyy_0_xxxxyzz_1,   \
                             g_0_yyyyy_0_xxxxyzzz_0,  \
                             g_0_yyyyy_0_xxxxyzzz_1,  \
                             g_0_yyyyy_0_xxxxzzz_1,   \
                             g_0_yyyyy_0_xxxxzzzz_0,  \
                             g_0_yyyyy_0_xxxxzzzz_1,  \
                             g_0_yyyyy_0_xxxyyyy_1,   \
                             g_0_yyyyy_0_xxxyyyyy_0,  \
                             g_0_yyyyy_0_xxxyyyyy_1,  \
                             g_0_yyyyy_0_xxxyyyyz_0,  \
                             g_0_yyyyy_0_xxxyyyyz_1,  \
                             g_0_yyyyy_0_xxxyyyz_1,   \
                             g_0_yyyyy_0_xxxyyyzz_0,  \
                             g_0_yyyyy_0_xxxyyyzz_1,  \
                             g_0_yyyyy_0_xxxyyzz_1,   \
                             g_0_yyyyy_0_xxxyyzzz_0,  \
                             g_0_yyyyy_0_xxxyyzzz_1,  \
                             g_0_yyyyy_0_xxxyzzz_1,   \
                             g_0_yyyyy_0_xxxyzzzz_0,  \
                             g_0_yyyyy_0_xxxyzzzz_1,  \
                             g_0_yyyyy_0_xxxzzzz_1,   \
                             g_0_yyyyy_0_xxxzzzzz_0,  \
                             g_0_yyyyy_0_xxxzzzzz_1,  \
                             g_0_yyyyy_0_xxyyyyy_1,   \
                             g_0_yyyyy_0_xxyyyyyy_0,  \
                             g_0_yyyyy_0_xxyyyyyy_1,  \
                             g_0_yyyyy_0_xxyyyyyz_0,  \
                             g_0_yyyyy_0_xxyyyyyz_1,  \
                             g_0_yyyyy_0_xxyyyyz_1,   \
                             g_0_yyyyy_0_xxyyyyzz_0,  \
                             g_0_yyyyy_0_xxyyyyzz_1,  \
                             g_0_yyyyy_0_xxyyyzz_1,   \
                             g_0_yyyyy_0_xxyyyzzz_0,  \
                             g_0_yyyyy_0_xxyyyzzz_1,  \
                             g_0_yyyyy_0_xxyyzzz_1,   \
                             g_0_yyyyy_0_xxyyzzzz_0,  \
                             g_0_yyyyy_0_xxyyzzzz_1,  \
                             g_0_yyyyy_0_xxyzzzz_1,   \
                             g_0_yyyyy_0_xxyzzzzz_0,  \
                             g_0_yyyyy_0_xxyzzzzz_1,  \
                             g_0_yyyyy_0_xxzzzzz_1,   \
                             g_0_yyyyy_0_xxzzzzzz_0,  \
                             g_0_yyyyy_0_xxzzzzzz_1,  \
                             g_0_yyyyy_0_xyyyyyy_1,   \
                             g_0_yyyyy_0_xyyyyyyy_0,  \
                             g_0_yyyyy_0_xyyyyyyy_1,  \
                             g_0_yyyyy_0_xyyyyyyz_0,  \
                             g_0_yyyyy_0_xyyyyyyz_1,  \
                             g_0_yyyyy_0_xyyyyyz_1,   \
                             g_0_yyyyy_0_xyyyyyzz_0,  \
                             g_0_yyyyy_0_xyyyyyzz_1,  \
                             g_0_yyyyy_0_xyyyyzz_1,   \
                             g_0_yyyyy_0_xyyyyzzz_0,  \
                             g_0_yyyyy_0_xyyyyzzz_1,  \
                             g_0_yyyyy_0_xyyyzzz_1,   \
                             g_0_yyyyy_0_xyyyzzzz_0,  \
                             g_0_yyyyy_0_xyyyzzzz_1,  \
                             g_0_yyyyy_0_xyyzzzz_1,   \
                             g_0_yyyyy_0_xyyzzzzz_0,  \
                             g_0_yyyyy_0_xyyzzzzz_1,  \
                             g_0_yyyyy_0_xyzzzzz_1,   \
                             g_0_yyyyy_0_xyzzzzzz_0,  \
                             g_0_yyyyy_0_xyzzzzzz_1,  \
                             g_0_yyyyy_0_xzzzzzz_1,   \
                             g_0_yyyyy_0_xzzzzzzz_0,  \
                             g_0_yyyyy_0_xzzzzzzz_1,  \
                             g_0_yyyyy_0_yyyyyyy_1,   \
                             g_0_yyyyy_0_yyyyyyyy_0,  \
                             g_0_yyyyy_0_yyyyyyyy_1,  \
                             g_0_yyyyy_0_yyyyyyyz_0,  \
                             g_0_yyyyy_0_yyyyyyyz_1,  \
                             g_0_yyyyy_0_yyyyyyz_1,   \
                             g_0_yyyyy_0_yyyyyyzz_0,  \
                             g_0_yyyyy_0_yyyyyyzz_1,  \
                             g_0_yyyyy_0_yyyyyzz_1,   \
                             g_0_yyyyy_0_yyyyyzzz_0,  \
                             g_0_yyyyy_0_yyyyyzzz_1,  \
                             g_0_yyyyy_0_yyyyzzz_1,   \
                             g_0_yyyyy_0_yyyyzzzz_0,  \
                             g_0_yyyyy_0_yyyyzzzz_1,  \
                             g_0_yyyyy_0_yyyzzzz_1,   \
                             g_0_yyyyy_0_yyyzzzzz_0,  \
                             g_0_yyyyy_0_yyyzzzzz_1,  \
                             g_0_yyyyy_0_yyzzzzz_1,   \
                             g_0_yyyyy_0_yyzzzzzz_0,  \
                             g_0_yyyyy_0_yyzzzzzz_1,  \
                             g_0_yyyyy_0_yzzzzzz_1,   \
                             g_0_yyyyy_0_yzzzzzzz_0,  \
                             g_0_yyyyy_0_yzzzzzzz_1,  \
                             g_0_yyyyy_0_zzzzzzz_1,   \
                             g_0_yyyyy_0_zzzzzzzz_0,  \
                             g_0_yyyyy_0_zzzzzzzz_1,  \
                             g_0_yyyyyz_0_xxxxxxxx_0, \
                             g_0_yyyyyz_0_xxxxxxxy_0, \
                             g_0_yyyyyz_0_xxxxxxxz_0, \
                             g_0_yyyyyz_0_xxxxxxyy_0, \
                             g_0_yyyyyz_0_xxxxxxyz_0, \
                             g_0_yyyyyz_0_xxxxxxzz_0, \
                             g_0_yyyyyz_0_xxxxxyyy_0, \
                             g_0_yyyyyz_0_xxxxxyyz_0, \
                             g_0_yyyyyz_0_xxxxxyzz_0, \
                             g_0_yyyyyz_0_xxxxxzzz_0, \
                             g_0_yyyyyz_0_xxxxyyyy_0, \
                             g_0_yyyyyz_0_xxxxyyyz_0, \
                             g_0_yyyyyz_0_xxxxyyzz_0, \
                             g_0_yyyyyz_0_xxxxyzzz_0, \
                             g_0_yyyyyz_0_xxxxzzzz_0, \
                             g_0_yyyyyz_0_xxxyyyyy_0, \
                             g_0_yyyyyz_0_xxxyyyyz_0, \
                             g_0_yyyyyz_0_xxxyyyzz_0, \
                             g_0_yyyyyz_0_xxxyyzzz_0, \
                             g_0_yyyyyz_0_xxxyzzzz_0, \
                             g_0_yyyyyz_0_xxxzzzzz_0, \
                             g_0_yyyyyz_0_xxyyyyyy_0, \
                             g_0_yyyyyz_0_xxyyyyyz_0, \
                             g_0_yyyyyz_0_xxyyyyzz_0, \
                             g_0_yyyyyz_0_xxyyyzzz_0, \
                             g_0_yyyyyz_0_xxyyzzzz_0, \
                             g_0_yyyyyz_0_xxyzzzzz_0, \
                             g_0_yyyyyz_0_xxzzzzzz_0, \
                             g_0_yyyyyz_0_xyyyyyyy_0, \
                             g_0_yyyyyz_0_xyyyyyyz_0, \
                             g_0_yyyyyz_0_xyyyyyzz_0, \
                             g_0_yyyyyz_0_xyyyyzzz_0, \
                             g_0_yyyyyz_0_xyyyzzzz_0, \
                             g_0_yyyyyz_0_xyyzzzzz_0, \
                             g_0_yyyyyz_0_xyzzzzzz_0, \
                             g_0_yyyyyz_0_xzzzzzzz_0, \
                             g_0_yyyyyz_0_yyyyyyyy_0, \
                             g_0_yyyyyz_0_yyyyyyyz_0, \
                             g_0_yyyyyz_0_yyyyyyzz_0, \
                             g_0_yyyyyz_0_yyyyyzzz_0, \
                             g_0_yyyyyz_0_yyyyzzzz_0, \
                             g_0_yyyyyz_0_yyyzzzzz_0, \
                             g_0_yyyyyz_0_yyzzzzzz_0, \
                             g_0_yyyyyz_0_yzzzzzzz_0, \
                             g_0_yyyyyz_0_zzzzzzzz_0, \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyz_0_xxxxxxxx_0[i] = g_0_yyyyy_0_xxxxxxxx_0[i] * pb_z + g_0_yyyyy_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxxxy_0[i] = g_0_yyyyy_0_xxxxxxxy_0[i] * pb_z + g_0_yyyyy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxxxz_0[i] = g_0_yyyyy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxxz_0[i] * pb_z + g_0_yyyyy_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxxyy_0[i] = g_0_yyyyy_0_xxxxxxyy_0[i] * pb_z + g_0_yyyyy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxxyz_0[i] = g_0_yyyyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxyz_0[i] * pb_z + g_0_yyyyy_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxxzz_0[i] =
            2.0 * g_0_yyyyy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxzz_0[i] * pb_z + g_0_yyyyy_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxyyy_0[i] = g_0_yyyyy_0_xxxxxyyy_0[i] * pb_z + g_0_yyyyy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxyyz_0[i] = g_0_yyyyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyyz_0[i] * pb_z + g_0_yyyyy_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxyzz_0[i] =
            2.0 * g_0_yyyyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyzz_0[i] * pb_z + g_0_yyyyy_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxzzz_0[i] =
            3.0 * g_0_yyyyy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxzzz_0[i] * pb_z + g_0_yyyyy_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxyyyy_0[i] = g_0_yyyyy_0_xxxxyyyy_0[i] * pb_z + g_0_yyyyy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxyyyz_0[i] = g_0_yyyyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyyz_0[i] * pb_z + g_0_yyyyy_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxyyzz_0[i] =
            2.0 * g_0_yyyyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyzz_0[i] * pb_z + g_0_yyyyy_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxyzzz_0[i] =
            3.0 * g_0_yyyyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyzzz_0[i] * pb_z + g_0_yyyyy_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxzzzz_0[i] =
            4.0 * g_0_yyyyy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxzzzz_0[i] * pb_z + g_0_yyyyy_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyyyyy_0[i] = g_0_yyyyy_0_xxxyyyyy_0[i] * pb_z + g_0_yyyyy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyyyyz_0[i] = g_0_yyyyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyyz_0[i] * pb_z + g_0_yyyyy_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyyyzz_0[i] =
            2.0 * g_0_yyyyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyzz_0[i] * pb_z + g_0_yyyyy_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyyzzz_0[i] =
            3.0 * g_0_yyyyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyzzz_0[i] * pb_z + g_0_yyyyy_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyzzzz_0[i] =
            4.0 * g_0_yyyyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyzzzz_0[i] * pb_z + g_0_yyyyy_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxzzzzz_0[i] =
            5.0 * g_0_yyyyy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxzzzzz_0[i] * pb_z + g_0_yyyyy_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyyyyy_0[i] = g_0_yyyyy_0_xxyyyyyy_0[i] * pb_z + g_0_yyyyy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyyyyz_0[i] = g_0_yyyyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyyz_0[i] * pb_z + g_0_yyyyy_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyyyzz_0[i] =
            2.0 * g_0_yyyyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyzz_0[i] * pb_z + g_0_yyyyy_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyyzzz_0[i] =
            3.0 * g_0_yyyyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyzzz_0[i] * pb_z + g_0_yyyyy_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyzzzz_0[i] =
            4.0 * g_0_yyyyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyzzzz_0[i] * pb_z + g_0_yyyyy_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyzzzzz_0[i] =
            5.0 * g_0_yyyyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyzzzzz_0[i] * pb_z + g_0_yyyyy_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxzzzzzz_0[i] =
            6.0 * g_0_yyyyy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxzzzzzz_0[i] * pb_z + g_0_yyyyy_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyyyyy_0[i] = g_0_yyyyy_0_xyyyyyyy_0[i] * pb_z + g_0_yyyyy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyyyyz_0[i] = g_0_yyyyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyyz_0[i] * pb_z + g_0_yyyyy_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyyyzz_0[i] =
            2.0 * g_0_yyyyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyzz_0[i] * pb_z + g_0_yyyyy_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyyzzz_0[i] =
            3.0 * g_0_yyyyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyzzz_0[i] * pb_z + g_0_yyyyy_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyzzzz_0[i] =
            4.0 * g_0_yyyyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyzzzz_0[i] * pb_z + g_0_yyyyy_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyzzzzz_0[i] =
            5.0 * g_0_yyyyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyzzzzz_0[i] * pb_z + g_0_yyyyy_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyzzzzzz_0[i] =
            6.0 * g_0_yyyyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzzzzzz_0[i] * pb_z + g_0_yyyyy_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xzzzzzzz_0[i] =
            7.0 * g_0_yyyyy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xzzzzzzz_0[i] * pb_z + g_0_yyyyy_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyyyyy_0[i] = g_0_yyyyy_0_yyyyyyyy_0[i] * pb_z + g_0_yyyyy_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyyyyz_0[i] = g_0_yyyyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyyyz_0[i] * pb_z + g_0_yyyyy_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyyyzz_0[i] =
            2.0 * g_0_yyyyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyyzz_0[i] * pb_z + g_0_yyyyy_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyyzzz_0[i] =
            3.0 * g_0_yyyyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyzzz_0[i] * pb_z + g_0_yyyyy_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyzzzz_0[i] =
            4.0 * g_0_yyyyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyzzzz_0[i] * pb_z + g_0_yyyyy_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyzzzzz_0[i] =
            5.0 * g_0_yyyyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyzzzzz_0[i] * pb_z + g_0_yyyyy_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyzzzzzz_0[i] =
            6.0 * g_0_yyyyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyzzzzzz_0[i] * pb_z + g_0_yyyyy_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yzzzzzzz_0[i] =
            7.0 * g_0_yyyyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yzzzzzzz_0[i] * pb_z + g_0_yyyyy_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_zzzzzzzz_0[i] =
            8.0 * g_0_yyyyy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_zzzzzzzz_0[i] * pb_z + g_0_yyyyy_0_zzzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 1035-1080 components of targeted buffer : SISL

    auto g_0_yyyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 1035);

    auto g_0_yyyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 1036);

    auto g_0_yyyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 1037);

    auto g_0_yyyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 1038);

    auto g_0_yyyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 1039);

    auto g_0_yyyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 1040);

    auto g_0_yyyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 1041);

    auto g_0_yyyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 1042);

    auto g_0_yyyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 1043);

    auto g_0_yyyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 1044);

    auto g_0_yyyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1045);

    auto g_0_yyyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1046);

    auto g_0_yyyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1047);

    auto g_0_yyyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1048);

    auto g_0_yyyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1049);

    auto g_0_yyyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1050);

    auto g_0_yyyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1051);

    auto g_0_yyyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1052);

    auto g_0_yyyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1053);

    auto g_0_yyyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1054);

    auto g_0_yyyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1055);

    auto g_0_yyyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1056);

    auto g_0_yyyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1057);

    auto g_0_yyyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1058);

    auto g_0_yyyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1059);

    auto g_0_yyyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1060);

    auto g_0_yyyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1061);

    auto g_0_yyyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1062);

    auto g_0_yyyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1063);

    auto g_0_yyyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1064);

    auto g_0_yyyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1065);

    auto g_0_yyyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1066);

    auto g_0_yyyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1067);

    auto g_0_yyyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1068);

    auto g_0_yyyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1069);

    auto g_0_yyyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1070);

    auto g_0_yyyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1071);

    auto g_0_yyyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1072);

    auto g_0_yyyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1073);

    auto g_0_yyyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1074);

    auto g_0_yyyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1075);

    auto g_0_yyyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1076);

    auto g_0_yyyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1077);

    auto g_0_yyyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1078);

    auto g_0_yyyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1079);

#pragma omp simd aligned(g_0_yyyy_0_xxxxxxxy_0,       \
                             g_0_yyyy_0_xxxxxxxy_1,   \
                             g_0_yyyy_0_xxxxxxyy_0,   \
                             g_0_yyyy_0_xxxxxxyy_1,   \
                             g_0_yyyy_0_xxxxxyyy_0,   \
                             g_0_yyyy_0_xxxxxyyy_1,   \
                             g_0_yyyy_0_xxxxyyyy_0,   \
                             g_0_yyyy_0_xxxxyyyy_1,   \
                             g_0_yyyy_0_xxxyyyyy_0,   \
                             g_0_yyyy_0_xxxyyyyy_1,   \
                             g_0_yyyy_0_xxyyyyyy_0,   \
                             g_0_yyyy_0_xxyyyyyy_1,   \
                             g_0_yyyy_0_xyyyyyyy_0,   \
                             g_0_yyyy_0_xyyyyyyy_1,   \
                             g_0_yyyy_0_yyyyyyyy_0,   \
                             g_0_yyyy_0_yyyyyyyy_1,   \
                             g_0_yyyyz_0_xxxxxxxy_0,  \
                             g_0_yyyyz_0_xxxxxxxy_1,  \
                             g_0_yyyyz_0_xxxxxxyy_0,  \
                             g_0_yyyyz_0_xxxxxxyy_1,  \
                             g_0_yyyyz_0_xxxxxyyy_0,  \
                             g_0_yyyyz_0_xxxxxyyy_1,  \
                             g_0_yyyyz_0_xxxxyyyy_0,  \
                             g_0_yyyyz_0_xxxxyyyy_1,  \
                             g_0_yyyyz_0_xxxyyyyy_0,  \
                             g_0_yyyyz_0_xxxyyyyy_1,  \
                             g_0_yyyyz_0_xxyyyyyy_0,  \
                             g_0_yyyyz_0_xxyyyyyy_1,  \
                             g_0_yyyyz_0_xyyyyyyy_0,  \
                             g_0_yyyyz_0_xyyyyyyy_1,  \
                             g_0_yyyyz_0_yyyyyyyy_0,  \
                             g_0_yyyyz_0_yyyyyyyy_1,  \
                             g_0_yyyyzz_0_xxxxxxxx_0, \
                             g_0_yyyyzz_0_xxxxxxxy_0, \
                             g_0_yyyyzz_0_xxxxxxxz_0, \
                             g_0_yyyyzz_0_xxxxxxyy_0, \
                             g_0_yyyyzz_0_xxxxxxyz_0, \
                             g_0_yyyyzz_0_xxxxxxzz_0, \
                             g_0_yyyyzz_0_xxxxxyyy_0, \
                             g_0_yyyyzz_0_xxxxxyyz_0, \
                             g_0_yyyyzz_0_xxxxxyzz_0, \
                             g_0_yyyyzz_0_xxxxxzzz_0, \
                             g_0_yyyyzz_0_xxxxyyyy_0, \
                             g_0_yyyyzz_0_xxxxyyyz_0, \
                             g_0_yyyyzz_0_xxxxyyzz_0, \
                             g_0_yyyyzz_0_xxxxyzzz_0, \
                             g_0_yyyyzz_0_xxxxzzzz_0, \
                             g_0_yyyyzz_0_xxxyyyyy_0, \
                             g_0_yyyyzz_0_xxxyyyyz_0, \
                             g_0_yyyyzz_0_xxxyyyzz_0, \
                             g_0_yyyyzz_0_xxxyyzzz_0, \
                             g_0_yyyyzz_0_xxxyzzzz_0, \
                             g_0_yyyyzz_0_xxxzzzzz_0, \
                             g_0_yyyyzz_0_xxyyyyyy_0, \
                             g_0_yyyyzz_0_xxyyyyyz_0, \
                             g_0_yyyyzz_0_xxyyyyzz_0, \
                             g_0_yyyyzz_0_xxyyyzzz_0, \
                             g_0_yyyyzz_0_xxyyzzzz_0, \
                             g_0_yyyyzz_0_xxyzzzzz_0, \
                             g_0_yyyyzz_0_xxzzzzzz_0, \
                             g_0_yyyyzz_0_xyyyyyyy_0, \
                             g_0_yyyyzz_0_xyyyyyyz_0, \
                             g_0_yyyyzz_0_xyyyyyzz_0, \
                             g_0_yyyyzz_0_xyyyyzzz_0, \
                             g_0_yyyyzz_0_xyyyzzzz_0, \
                             g_0_yyyyzz_0_xyyzzzzz_0, \
                             g_0_yyyyzz_0_xyzzzzzz_0, \
                             g_0_yyyyzz_0_xzzzzzzz_0, \
                             g_0_yyyyzz_0_yyyyyyyy_0, \
                             g_0_yyyyzz_0_yyyyyyyz_0, \
                             g_0_yyyyzz_0_yyyyyyzz_0, \
                             g_0_yyyyzz_0_yyyyyzzz_0, \
                             g_0_yyyyzz_0_yyyyzzzz_0, \
                             g_0_yyyyzz_0_yyyzzzzz_0, \
                             g_0_yyyyzz_0_yyzzzzzz_0, \
                             g_0_yyyyzz_0_yzzzzzzz_0, \
                             g_0_yyyyzz_0_zzzzzzzz_0, \
                             g_0_yyyzz_0_xxxxxxxx_0,  \
                             g_0_yyyzz_0_xxxxxxxx_1,  \
                             g_0_yyyzz_0_xxxxxxxz_0,  \
                             g_0_yyyzz_0_xxxxxxxz_1,  \
                             g_0_yyyzz_0_xxxxxxyz_0,  \
                             g_0_yyyzz_0_xxxxxxyz_1,  \
                             g_0_yyyzz_0_xxxxxxz_1,   \
                             g_0_yyyzz_0_xxxxxxzz_0,  \
                             g_0_yyyzz_0_xxxxxxzz_1,  \
                             g_0_yyyzz_0_xxxxxyyz_0,  \
                             g_0_yyyzz_0_xxxxxyyz_1,  \
                             g_0_yyyzz_0_xxxxxyz_1,   \
                             g_0_yyyzz_0_xxxxxyzz_0,  \
                             g_0_yyyzz_0_xxxxxyzz_1,  \
                             g_0_yyyzz_0_xxxxxzz_1,   \
                             g_0_yyyzz_0_xxxxxzzz_0,  \
                             g_0_yyyzz_0_xxxxxzzz_1,  \
                             g_0_yyyzz_0_xxxxyyyz_0,  \
                             g_0_yyyzz_0_xxxxyyyz_1,  \
                             g_0_yyyzz_0_xxxxyyz_1,   \
                             g_0_yyyzz_0_xxxxyyzz_0,  \
                             g_0_yyyzz_0_xxxxyyzz_1,  \
                             g_0_yyyzz_0_xxxxyzz_1,   \
                             g_0_yyyzz_0_xxxxyzzz_0,  \
                             g_0_yyyzz_0_xxxxyzzz_1,  \
                             g_0_yyyzz_0_xxxxzzz_1,   \
                             g_0_yyyzz_0_xxxxzzzz_0,  \
                             g_0_yyyzz_0_xxxxzzzz_1,  \
                             g_0_yyyzz_0_xxxyyyyz_0,  \
                             g_0_yyyzz_0_xxxyyyyz_1,  \
                             g_0_yyyzz_0_xxxyyyz_1,   \
                             g_0_yyyzz_0_xxxyyyzz_0,  \
                             g_0_yyyzz_0_xxxyyyzz_1,  \
                             g_0_yyyzz_0_xxxyyzz_1,   \
                             g_0_yyyzz_0_xxxyyzzz_0,  \
                             g_0_yyyzz_0_xxxyyzzz_1,  \
                             g_0_yyyzz_0_xxxyzzz_1,   \
                             g_0_yyyzz_0_xxxyzzzz_0,  \
                             g_0_yyyzz_0_xxxyzzzz_1,  \
                             g_0_yyyzz_0_xxxzzzz_1,   \
                             g_0_yyyzz_0_xxxzzzzz_0,  \
                             g_0_yyyzz_0_xxxzzzzz_1,  \
                             g_0_yyyzz_0_xxyyyyyz_0,  \
                             g_0_yyyzz_0_xxyyyyyz_1,  \
                             g_0_yyyzz_0_xxyyyyz_1,   \
                             g_0_yyyzz_0_xxyyyyzz_0,  \
                             g_0_yyyzz_0_xxyyyyzz_1,  \
                             g_0_yyyzz_0_xxyyyzz_1,   \
                             g_0_yyyzz_0_xxyyyzzz_0,  \
                             g_0_yyyzz_0_xxyyyzzz_1,  \
                             g_0_yyyzz_0_xxyyzzz_1,   \
                             g_0_yyyzz_0_xxyyzzzz_0,  \
                             g_0_yyyzz_0_xxyyzzzz_1,  \
                             g_0_yyyzz_0_xxyzzzz_1,   \
                             g_0_yyyzz_0_xxyzzzzz_0,  \
                             g_0_yyyzz_0_xxyzzzzz_1,  \
                             g_0_yyyzz_0_xxzzzzz_1,   \
                             g_0_yyyzz_0_xxzzzzzz_0,  \
                             g_0_yyyzz_0_xxzzzzzz_1,  \
                             g_0_yyyzz_0_xyyyyyyz_0,  \
                             g_0_yyyzz_0_xyyyyyyz_1,  \
                             g_0_yyyzz_0_xyyyyyz_1,   \
                             g_0_yyyzz_0_xyyyyyzz_0,  \
                             g_0_yyyzz_0_xyyyyyzz_1,  \
                             g_0_yyyzz_0_xyyyyzz_1,   \
                             g_0_yyyzz_0_xyyyyzzz_0,  \
                             g_0_yyyzz_0_xyyyyzzz_1,  \
                             g_0_yyyzz_0_xyyyzzz_1,   \
                             g_0_yyyzz_0_xyyyzzzz_0,  \
                             g_0_yyyzz_0_xyyyzzzz_1,  \
                             g_0_yyyzz_0_xyyzzzz_1,   \
                             g_0_yyyzz_0_xyyzzzzz_0,  \
                             g_0_yyyzz_0_xyyzzzzz_1,  \
                             g_0_yyyzz_0_xyzzzzz_1,   \
                             g_0_yyyzz_0_xyzzzzzz_0,  \
                             g_0_yyyzz_0_xyzzzzzz_1,  \
                             g_0_yyyzz_0_xzzzzzz_1,   \
                             g_0_yyyzz_0_xzzzzzzz_0,  \
                             g_0_yyyzz_0_xzzzzzzz_1,  \
                             g_0_yyyzz_0_yyyyyyyz_0,  \
                             g_0_yyyzz_0_yyyyyyyz_1,  \
                             g_0_yyyzz_0_yyyyyyz_1,   \
                             g_0_yyyzz_0_yyyyyyzz_0,  \
                             g_0_yyyzz_0_yyyyyyzz_1,  \
                             g_0_yyyzz_0_yyyyyzz_1,   \
                             g_0_yyyzz_0_yyyyyzzz_0,  \
                             g_0_yyyzz_0_yyyyyzzz_1,  \
                             g_0_yyyzz_0_yyyyzzz_1,   \
                             g_0_yyyzz_0_yyyyzzzz_0,  \
                             g_0_yyyzz_0_yyyyzzzz_1,  \
                             g_0_yyyzz_0_yyyzzzz_1,   \
                             g_0_yyyzz_0_yyyzzzzz_0,  \
                             g_0_yyyzz_0_yyyzzzzz_1,  \
                             g_0_yyyzz_0_yyzzzzz_1,   \
                             g_0_yyyzz_0_yyzzzzzz_0,  \
                             g_0_yyyzz_0_yyzzzzzz_1,  \
                             g_0_yyyzz_0_yzzzzzz_1,   \
                             g_0_yyyzz_0_yzzzzzzz_0,  \
                             g_0_yyyzz_0_yzzzzzzz_1,  \
                             g_0_yyyzz_0_zzzzzzz_1,   \
                             g_0_yyyzz_0_zzzzzzzz_0,  \
                             g_0_yyyzz_0_zzzzzzzz_1,  \
                             g_0_yyzz_0_xxxxxxxx_0,   \
                             g_0_yyzz_0_xxxxxxxx_1,   \
                             g_0_yyzz_0_xxxxxxxz_0,   \
                             g_0_yyzz_0_xxxxxxxz_1,   \
                             g_0_yyzz_0_xxxxxxyz_0,   \
                             g_0_yyzz_0_xxxxxxyz_1,   \
                             g_0_yyzz_0_xxxxxxzz_0,   \
                             g_0_yyzz_0_xxxxxxzz_1,   \
                             g_0_yyzz_0_xxxxxyyz_0,   \
                             g_0_yyzz_0_xxxxxyyz_1,   \
                             g_0_yyzz_0_xxxxxyzz_0,   \
                             g_0_yyzz_0_xxxxxyzz_1,   \
                             g_0_yyzz_0_xxxxxzzz_0,   \
                             g_0_yyzz_0_xxxxxzzz_1,   \
                             g_0_yyzz_0_xxxxyyyz_0,   \
                             g_0_yyzz_0_xxxxyyyz_1,   \
                             g_0_yyzz_0_xxxxyyzz_0,   \
                             g_0_yyzz_0_xxxxyyzz_1,   \
                             g_0_yyzz_0_xxxxyzzz_0,   \
                             g_0_yyzz_0_xxxxyzzz_1,   \
                             g_0_yyzz_0_xxxxzzzz_0,   \
                             g_0_yyzz_0_xxxxzzzz_1,   \
                             g_0_yyzz_0_xxxyyyyz_0,   \
                             g_0_yyzz_0_xxxyyyyz_1,   \
                             g_0_yyzz_0_xxxyyyzz_0,   \
                             g_0_yyzz_0_xxxyyyzz_1,   \
                             g_0_yyzz_0_xxxyyzzz_0,   \
                             g_0_yyzz_0_xxxyyzzz_1,   \
                             g_0_yyzz_0_xxxyzzzz_0,   \
                             g_0_yyzz_0_xxxyzzzz_1,   \
                             g_0_yyzz_0_xxxzzzzz_0,   \
                             g_0_yyzz_0_xxxzzzzz_1,   \
                             g_0_yyzz_0_xxyyyyyz_0,   \
                             g_0_yyzz_0_xxyyyyyz_1,   \
                             g_0_yyzz_0_xxyyyyzz_0,   \
                             g_0_yyzz_0_xxyyyyzz_1,   \
                             g_0_yyzz_0_xxyyyzzz_0,   \
                             g_0_yyzz_0_xxyyyzzz_1,   \
                             g_0_yyzz_0_xxyyzzzz_0,   \
                             g_0_yyzz_0_xxyyzzzz_1,   \
                             g_0_yyzz_0_xxyzzzzz_0,   \
                             g_0_yyzz_0_xxyzzzzz_1,   \
                             g_0_yyzz_0_xxzzzzzz_0,   \
                             g_0_yyzz_0_xxzzzzzz_1,   \
                             g_0_yyzz_0_xyyyyyyz_0,   \
                             g_0_yyzz_0_xyyyyyyz_1,   \
                             g_0_yyzz_0_xyyyyyzz_0,   \
                             g_0_yyzz_0_xyyyyyzz_1,   \
                             g_0_yyzz_0_xyyyyzzz_0,   \
                             g_0_yyzz_0_xyyyyzzz_1,   \
                             g_0_yyzz_0_xyyyzzzz_0,   \
                             g_0_yyzz_0_xyyyzzzz_1,   \
                             g_0_yyzz_0_xyyzzzzz_0,   \
                             g_0_yyzz_0_xyyzzzzz_1,   \
                             g_0_yyzz_0_xyzzzzzz_0,   \
                             g_0_yyzz_0_xyzzzzzz_1,   \
                             g_0_yyzz_0_xzzzzzzz_0,   \
                             g_0_yyzz_0_xzzzzzzz_1,   \
                             g_0_yyzz_0_yyyyyyyz_0,   \
                             g_0_yyzz_0_yyyyyyyz_1,   \
                             g_0_yyzz_0_yyyyyyzz_0,   \
                             g_0_yyzz_0_yyyyyyzz_1,   \
                             g_0_yyzz_0_yyyyyzzz_0,   \
                             g_0_yyzz_0_yyyyyzzz_1,   \
                             g_0_yyzz_0_yyyyzzzz_0,   \
                             g_0_yyzz_0_yyyyzzzz_1,   \
                             g_0_yyzz_0_yyyzzzzz_0,   \
                             g_0_yyzz_0_yyyzzzzz_1,   \
                             g_0_yyzz_0_yyzzzzzz_0,   \
                             g_0_yyzz_0_yyzzzzzz_1,   \
                             g_0_yyzz_0_yzzzzzzz_0,   \
                             g_0_yyzz_0_yzzzzzzz_1,   \
                             g_0_yyzz_0_zzzzzzzz_0,   \
                             g_0_yyzz_0_zzzzzzzz_1,   \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzz_0_xxxxxxxx_0[i] = 3.0 * g_0_yyzz_0_xxxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxxxx_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxxxxxx_0[i] * pb_y + g_0_yyyzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxxxxy_0[i] = g_0_yyyy_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxxxxxy_0[i] * pb_z +
                                     g_0_yyyyz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxxxxxz_0[i] = 3.0 * g_0_yyzz_0_xxxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxxxz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxxxxxz_0[i] * pb_y + g_0_yyyzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxxxyy_0[i] = g_0_yyyy_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxxxxyy_0[i] * pb_z +
                                     g_0_yyyyz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxxxxyz_0[i] = 3.0 * g_0_yyzz_0_xxxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxxyz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxxyz_0[i] * pb_y + g_0_yyyzz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxxxzz_0[i] = 3.0 * g_0_yyzz_0_xxxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxxzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxxxxzz_0[i] * pb_y + g_0_yyyzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxxyyy_0[i] = g_0_yyyy_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxxxyyy_0[i] * pb_z +
                                     g_0_yyyyz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxxxyyz_0[i] = 3.0 * g_0_yyzz_0_xxxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxyyz_0[i] * pb_y +
                                     g_0_yyyzz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxxyzz_0[i] = 3.0 * g_0_yyzz_0_xxxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxyzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxyzz_0[i] * pb_y + g_0_yyyzz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxxzzz_0[i] = 3.0 * g_0_yyzz_0_xxxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxxxzzz_0[i] * pb_y + g_0_yyyzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxyyyy_0[i] = g_0_yyyy_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxxyyyy_0[i] * pb_z +
                                     g_0_yyyyz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxxyyyz_0[i] = 3.0 * g_0_yyzz_0_xxxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyyyz_0[i] * pb_y +
                                     g_0_yyyzz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxyyzz_0[i] = 3.0 * g_0_yyzz_0_xxxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyyzz_0[i] * pb_y +
                                     g_0_yyyzz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxyzzz_0[i] = 3.0 * g_0_yyzz_0_xxxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxyzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyzzz_0[i] * pb_y + g_0_yyyzz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxzzzz_0[i] = 3.0 * g_0_yyzz_0_xxxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxxzzzz_0[i] * pb_y + g_0_yyyzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxyyyyy_0[i] = g_0_yyyy_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxyyyyy_0[i] * pb_z +
                                     g_0_yyyyz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxyyyyz_0[i] = 3.0 * g_0_yyzz_0_xxxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyyyz_0[i] * pb_y +
                                     g_0_yyyzz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxyyyzz_0[i] = 3.0 * g_0_yyzz_0_xxxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyyzz_0[i] * pb_y +
                                     g_0_yyyzz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxyyzzz_0[i] = 3.0 * g_0_yyzz_0_xxxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyzzz_0[i] * pb_y +
                                     g_0_yyyzz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxyzzzz_0[i] = 3.0 * g_0_yyzz_0_xxxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyzzzz_0[i] * pb_y + g_0_yyyzz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxzzzzz_0[i] = 3.0 * g_0_yyzz_0_xxxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxzzzzz_0[i] * pb_y + g_0_yyyzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyyyyyy_0[i] = g_0_yyyy_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxyyyyyy_0[i] * pb_z +
                                     g_0_yyyyz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxyyyyyz_0[i] = 3.0 * g_0_yyzz_0_xxyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyyzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyyyz_0[i] * pb_y +
                                     g_0_yyyzz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyyyyzz_0[i] = 3.0 * g_0_yyzz_0_xxyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyyzz_0[i] * pb_y +
                                     g_0_yyyzz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyyyzzz_0[i] = 3.0 * g_0_yyzz_0_xxyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyzzz_0[i] * pb_y +
                                     g_0_yyyzz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyyzzzz_0[i] = 3.0 * g_0_yyzz_0_xxyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyzzzz_0[i] * pb_y +
                                     g_0_yyyzz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyzzzzz_0[i] = 3.0 * g_0_yyzz_0_xxyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyzzzzz_0[i] * pb_y + g_0_yyyzz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxzzzzzz_0[i] = 3.0 * g_0_yyzz_0_xxzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxzzzzzz_0[i] * pb_y + g_0_yyyzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyyyyyy_0[i] = g_0_yyyy_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xyyyyyyy_0[i] * pb_z +
                                     g_0_yyyyz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xyyyyyyz_0[i] = 3.0 * g_0_yyzz_0_xyyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyyyyz_1[i] * fti_ab_0 +
                                     6.0 * g_0_yyyzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyyyz_0[i] * pb_y +
                                     g_0_yyyzz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyyyyzz_0[i] = 3.0 * g_0_yyzz_0_xyyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyyyzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyyzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyyzz_0[i] * pb_y +
                                     g_0_yyyzz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyyyzzz_0[i] = 3.0 * g_0_yyzz_0_xyyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyyzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyzzz_0[i] * pb_y +
                                     g_0_yyyzz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyyzzzz_0[i] = 3.0 * g_0_yyzz_0_xyyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyzzzz_0[i] * pb_y +
                                     g_0_yyyzz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyzzzzz_0[i] = 3.0 * g_0_yyzz_0_xyyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyzzzzz_0[i] * pb_y +
                                     g_0_yyyzz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyzzzzzz_0[i] = 3.0 * g_0_yyzz_0_xyzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyzzzzzz_0[i] * pb_y + g_0_yyyzz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xzzzzzzz_0[i] = 3.0 * g_0_yyzz_0_xzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xzzzzzzz_0[i] * pb_y + g_0_yyyzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyyyyyy_0[i] = g_0_yyyy_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_yyyyyyyy_0[i] * pb_z +
                                     g_0_yyyyz_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_yyyyyyyz_0[i] = 3.0 * g_0_yyzz_0_yyyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyyyyz_1[i] * fti_ab_0 +
                                     7.0 * g_0_yyyzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyyyyyz_0[i] * pb_y +
                                     g_0_yyyzz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyyyyzz_0[i] = 3.0 * g_0_yyzz_0_yyyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyyyzz_1[i] * fti_ab_0 +
                                     6.0 * g_0_yyyzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyyyyzz_0[i] * pb_y +
                                     g_0_yyyzz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyyyzzz_0[i] = 3.0 * g_0_yyzz_0_yyyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyyzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyyzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyyyzzz_0[i] * pb_y +
                                     g_0_yyyzz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyyzzzz_0[i] = 3.0 * g_0_yyzz_0_yyyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyyzzzz_0[i] * pb_y +
                                     g_0_yyyzz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyzzzzz_0[i] = 3.0 * g_0_yyzz_0_yyyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyzzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyzzzzz_0[i] * pb_y +
                                     g_0_yyyzz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyzzzzzz_0[i] = 3.0 * g_0_yyzz_0_yyzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyzzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyzzzzzz_0[i] * pb_y +
                                     g_0_yyyzz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yzzzzzzz_0[i] = 3.0 * g_0_yyzz_0_yzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yzzzzzzz_0[i] * pb_y + g_0_yyyzz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_zzzzzzzz_0[i] = 3.0 * g_0_yyzz_0_zzzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_zzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_zzzzzzzz_0[i] * pb_y + g_0_yyyzz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1080-1125 components of targeted buffer : SISL

    auto g_0_yyyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 1080);

    auto g_0_yyyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 1081);

    auto g_0_yyyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 1082);

    auto g_0_yyyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 1083);

    auto g_0_yyyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 1084);

    auto g_0_yyyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 1085);

    auto g_0_yyyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 1086);

    auto g_0_yyyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 1087);

    auto g_0_yyyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 1088);

    auto g_0_yyyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 1089);

    auto g_0_yyyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1090);

    auto g_0_yyyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1091);

    auto g_0_yyyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1092);

    auto g_0_yyyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1093);

    auto g_0_yyyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1094);

    auto g_0_yyyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1095);

    auto g_0_yyyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1096);

    auto g_0_yyyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1097);

    auto g_0_yyyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1098);

    auto g_0_yyyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1099);

    auto g_0_yyyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1100);

    auto g_0_yyyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1101);

    auto g_0_yyyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1102);

    auto g_0_yyyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1103);

    auto g_0_yyyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1104);

    auto g_0_yyyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1105);

    auto g_0_yyyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1106);

    auto g_0_yyyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1107);

    auto g_0_yyyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1108);

    auto g_0_yyyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1109);

    auto g_0_yyyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1110);

    auto g_0_yyyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1111);

    auto g_0_yyyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1112);

    auto g_0_yyyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1113);

    auto g_0_yyyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1114);

    auto g_0_yyyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1115);

    auto g_0_yyyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1116);

    auto g_0_yyyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1117);

    auto g_0_yyyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1118);

    auto g_0_yyyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1119);

    auto g_0_yyyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1120);

    auto g_0_yyyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1121);

    auto g_0_yyyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1122);

    auto g_0_yyyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1123);

    auto g_0_yyyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1124);

#pragma omp simd aligned(g_0_yyyz_0_xxxxxxxy_0,       \
                             g_0_yyyz_0_xxxxxxxy_1,   \
                             g_0_yyyz_0_xxxxxxyy_0,   \
                             g_0_yyyz_0_xxxxxxyy_1,   \
                             g_0_yyyz_0_xxxxxyyy_0,   \
                             g_0_yyyz_0_xxxxxyyy_1,   \
                             g_0_yyyz_0_xxxxyyyy_0,   \
                             g_0_yyyz_0_xxxxyyyy_1,   \
                             g_0_yyyz_0_xxxyyyyy_0,   \
                             g_0_yyyz_0_xxxyyyyy_1,   \
                             g_0_yyyz_0_xxyyyyyy_0,   \
                             g_0_yyyz_0_xxyyyyyy_1,   \
                             g_0_yyyz_0_xyyyyyyy_0,   \
                             g_0_yyyz_0_xyyyyyyy_1,   \
                             g_0_yyyz_0_yyyyyyyy_0,   \
                             g_0_yyyz_0_yyyyyyyy_1,   \
                             g_0_yyyzz_0_xxxxxxxy_0,  \
                             g_0_yyyzz_0_xxxxxxxy_1,  \
                             g_0_yyyzz_0_xxxxxxyy_0,  \
                             g_0_yyyzz_0_xxxxxxyy_1,  \
                             g_0_yyyzz_0_xxxxxyyy_0,  \
                             g_0_yyyzz_0_xxxxxyyy_1,  \
                             g_0_yyyzz_0_xxxxyyyy_0,  \
                             g_0_yyyzz_0_xxxxyyyy_1,  \
                             g_0_yyyzz_0_xxxyyyyy_0,  \
                             g_0_yyyzz_0_xxxyyyyy_1,  \
                             g_0_yyyzz_0_xxyyyyyy_0,  \
                             g_0_yyyzz_0_xxyyyyyy_1,  \
                             g_0_yyyzz_0_xyyyyyyy_0,  \
                             g_0_yyyzz_0_xyyyyyyy_1,  \
                             g_0_yyyzz_0_yyyyyyyy_0,  \
                             g_0_yyyzz_0_yyyyyyyy_1,  \
                             g_0_yyyzzz_0_xxxxxxxx_0, \
                             g_0_yyyzzz_0_xxxxxxxy_0, \
                             g_0_yyyzzz_0_xxxxxxxz_0, \
                             g_0_yyyzzz_0_xxxxxxyy_0, \
                             g_0_yyyzzz_0_xxxxxxyz_0, \
                             g_0_yyyzzz_0_xxxxxxzz_0, \
                             g_0_yyyzzz_0_xxxxxyyy_0, \
                             g_0_yyyzzz_0_xxxxxyyz_0, \
                             g_0_yyyzzz_0_xxxxxyzz_0, \
                             g_0_yyyzzz_0_xxxxxzzz_0, \
                             g_0_yyyzzz_0_xxxxyyyy_0, \
                             g_0_yyyzzz_0_xxxxyyyz_0, \
                             g_0_yyyzzz_0_xxxxyyzz_0, \
                             g_0_yyyzzz_0_xxxxyzzz_0, \
                             g_0_yyyzzz_0_xxxxzzzz_0, \
                             g_0_yyyzzz_0_xxxyyyyy_0, \
                             g_0_yyyzzz_0_xxxyyyyz_0, \
                             g_0_yyyzzz_0_xxxyyyzz_0, \
                             g_0_yyyzzz_0_xxxyyzzz_0, \
                             g_0_yyyzzz_0_xxxyzzzz_0, \
                             g_0_yyyzzz_0_xxxzzzzz_0, \
                             g_0_yyyzzz_0_xxyyyyyy_0, \
                             g_0_yyyzzz_0_xxyyyyyz_0, \
                             g_0_yyyzzz_0_xxyyyyzz_0, \
                             g_0_yyyzzz_0_xxyyyzzz_0, \
                             g_0_yyyzzz_0_xxyyzzzz_0, \
                             g_0_yyyzzz_0_xxyzzzzz_0, \
                             g_0_yyyzzz_0_xxzzzzzz_0, \
                             g_0_yyyzzz_0_xyyyyyyy_0, \
                             g_0_yyyzzz_0_xyyyyyyz_0, \
                             g_0_yyyzzz_0_xyyyyyzz_0, \
                             g_0_yyyzzz_0_xyyyyzzz_0, \
                             g_0_yyyzzz_0_xyyyzzzz_0, \
                             g_0_yyyzzz_0_xyyzzzzz_0, \
                             g_0_yyyzzz_0_xyzzzzzz_0, \
                             g_0_yyyzzz_0_xzzzzzzz_0, \
                             g_0_yyyzzz_0_yyyyyyyy_0, \
                             g_0_yyyzzz_0_yyyyyyyz_0, \
                             g_0_yyyzzz_0_yyyyyyzz_0, \
                             g_0_yyyzzz_0_yyyyyzzz_0, \
                             g_0_yyyzzz_0_yyyyzzzz_0, \
                             g_0_yyyzzz_0_yyyzzzzz_0, \
                             g_0_yyyzzz_0_yyzzzzzz_0, \
                             g_0_yyyzzz_0_yzzzzzzz_0, \
                             g_0_yyyzzz_0_zzzzzzzz_0, \
                             g_0_yyzzz_0_xxxxxxxx_0,  \
                             g_0_yyzzz_0_xxxxxxxx_1,  \
                             g_0_yyzzz_0_xxxxxxxz_0,  \
                             g_0_yyzzz_0_xxxxxxxz_1,  \
                             g_0_yyzzz_0_xxxxxxyz_0,  \
                             g_0_yyzzz_0_xxxxxxyz_1,  \
                             g_0_yyzzz_0_xxxxxxz_1,   \
                             g_0_yyzzz_0_xxxxxxzz_0,  \
                             g_0_yyzzz_0_xxxxxxzz_1,  \
                             g_0_yyzzz_0_xxxxxyyz_0,  \
                             g_0_yyzzz_0_xxxxxyyz_1,  \
                             g_0_yyzzz_0_xxxxxyz_1,   \
                             g_0_yyzzz_0_xxxxxyzz_0,  \
                             g_0_yyzzz_0_xxxxxyzz_1,  \
                             g_0_yyzzz_0_xxxxxzz_1,   \
                             g_0_yyzzz_0_xxxxxzzz_0,  \
                             g_0_yyzzz_0_xxxxxzzz_1,  \
                             g_0_yyzzz_0_xxxxyyyz_0,  \
                             g_0_yyzzz_0_xxxxyyyz_1,  \
                             g_0_yyzzz_0_xxxxyyz_1,   \
                             g_0_yyzzz_0_xxxxyyzz_0,  \
                             g_0_yyzzz_0_xxxxyyzz_1,  \
                             g_0_yyzzz_0_xxxxyzz_1,   \
                             g_0_yyzzz_0_xxxxyzzz_0,  \
                             g_0_yyzzz_0_xxxxyzzz_1,  \
                             g_0_yyzzz_0_xxxxzzz_1,   \
                             g_0_yyzzz_0_xxxxzzzz_0,  \
                             g_0_yyzzz_0_xxxxzzzz_1,  \
                             g_0_yyzzz_0_xxxyyyyz_0,  \
                             g_0_yyzzz_0_xxxyyyyz_1,  \
                             g_0_yyzzz_0_xxxyyyz_1,   \
                             g_0_yyzzz_0_xxxyyyzz_0,  \
                             g_0_yyzzz_0_xxxyyyzz_1,  \
                             g_0_yyzzz_0_xxxyyzz_1,   \
                             g_0_yyzzz_0_xxxyyzzz_0,  \
                             g_0_yyzzz_0_xxxyyzzz_1,  \
                             g_0_yyzzz_0_xxxyzzz_1,   \
                             g_0_yyzzz_0_xxxyzzzz_0,  \
                             g_0_yyzzz_0_xxxyzzzz_1,  \
                             g_0_yyzzz_0_xxxzzzz_1,   \
                             g_0_yyzzz_0_xxxzzzzz_0,  \
                             g_0_yyzzz_0_xxxzzzzz_1,  \
                             g_0_yyzzz_0_xxyyyyyz_0,  \
                             g_0_yyzzz_0_xxyyyyyz_1,  \
                             g_0_yyzzz_0_xxyyyyz_1,   \
                             g_0_yyzzz_0_xxyyyyzz_0,  \
                             g_0_yyzzz_0_xxyyyyzz_1,  \
                             g_0_yyzzz_0_xxyyyzz_1,   \
                             g_0_yyzzz_0_xxyyyzzz_0,  \
                             g_0_yyzzz_0_xxyyyzzz_1,  \
                             g_0_yyzzz_0_xxyyzzz_1,   \
                             g_0_yyzzz_0_xxyyzzzz_0,  \
                             g_0_yyzzz_0_xxyyzzzz_1,  \
                             g_0_yyzzz_0_xxyzzzz_1,   \
                             g_0_yyzzz_0_xxyzzzzz_0,  \
                             g_0_yyzzz_0_xxyzzzzz_1,  \
                             g_0_yyzzz_0_xxzzzzz_1,   \
                             g_0_yyzzz_0_xxzzzzzz_0,  \
                             g_0_yyzzz_0_xxzzzzzz_1,  \
                             g_0_yyzzz_0_xyyyyyyz_0,  \
                             g_0_yyzzz_0_xyyyyyyz_1,  \
                             g_0_yyzzz_0_xyyyyyz_1,   \
                             g_0_yyzzz_0_xyyyyyzz_0,  \
                             g_0_yyzzz_0_xyyyyyzz_1,  \
                             g_0_yyzzz_0_xyyyyzz_1,   \
                             g_0_yyzzz_0_xyyyyzzz_0,  \
                             g_0_yyzzz_0_xyyyyzzz_1,  \
                             g_0_yyzzz_0_xyyyzzz_1,   \
                             g_0_yyzzz_0_xyyyzzzz_0,  \
                             g_0_yyzzz_0_xyyyzzzz_1,  \
                             g_0_yyzzz_0_xyyzzzz_1,   \
                             g_0_yyzzz_0_xyyzzzzz_0,  \
                             g_0_yyzzz_0_xyyzzzzz_1,  \
                             g_0_yyzzz_0_xyzzzzz_1,   \
                             g_0_yyzzz_0_xyzzzzzz_0,  \
                             g_0_yyzzz_0_xyzzzzzz_1,  \
                             g_0_yyzzz_0_xzzzzzz_1,   \
                             g_0_yyzzz_0_xzzzzzzz_0,  \
                             g_0_yyzzz_0_xzzzzzzz_1,  \
                             g_0_yyzzz_0_yyyyyyyz_0,  \
                             g_0_yyzzz_0_yyyyyyyz_1,  \
                             g_0_yyzzz_0_yyyyyyz_1,   \
                             g_0_yyzzz_0_yyyyyyzz_0,  \
                             g_0_yyzzz_0_yyyyyyzz_1,  \
                             g_0_yyzzz_0_yyyyyzz_1,   \
                             g_0_yyzzz_0_yyyyyzzz_0,  \
                             g_0_yyzzz_0_yyyyyzzz_1,  \
                             g_0_yyzzz_0_yyyyzzz_1,   \
                             g_0_yyzzz_0_yyyyzzzz_0,  \
                             g_0_yyzzz_0_yyyyzzzz_1,  \
                             g_0_yyzzz_0_yyyzzzz_1,   \
                             g_0_yyzzz_0_yyyzzzzz_0,  \
                             g_0_yyzzz_0_yyyzzzzz_1,  \
                             g_0_yyzzz_0_yyzzzzz_1,   \
                             g_0_yyzzz_0_yyzzzzzz_0,  \
                             g_0_yyzzz_0_yyzzzzzz_1,  \
                             g_0_yyzzz_0_yzzzzzz_1,   \
                             g_0_yyzzz_0_yzzzzzzz_0,  \
                             g_0_yyzzz_0_yzzzzzzz_1,  \
                             g_0_yyzzz_0_zzzzzzz_1,   \
                             g_0_yyzzz_0_zzzzzzzz_0,  \
                             g_0_yyzzz_0_zzzzzzzz_1,  \
                             g_0_yzzz_0_xxxxxxxx_0,   \
                             g_0_yzzz_0_xxxxxxxx_1,   \
                             g_0_yzzz_0_xxxxxxxz_0,   \
                             g_0_yzzz_0_xxxxxxxz_1,   \
                             g_0_yzzz_0_xxxxxxyz_0,   \
                             g_0_yzzz_0_xxxxxxyz_1,   \
                             g_0_yzzz_0_xxxxxxzz_0,   \
                             g_0_yzzz_0_xxxxxxzz_1,   \
                             g_0_yzzz_0_xxxxxyyz_0,   \
                             g_0_yzzz_0_xxxxxyyz_1,   \
                             g_0_yzzz_0_xxxxxyzz_0,   \
                             g_0_yzzz_0_xxxxxyzz_1,   \
                             g_0_yzzz_0_xxxxxzzz_0,   \
                             g_0_yzzz_0_xxxxxzzz_1,   \
                             g_0_yzzz_0_xxxxyyyz_0,   \
                             g_0_yzzz_0_xxxxyyyz_1,   \
                             g_0_yzzz_0_xxxxyyzz_0,   \
                             g_0_yzzz_0_xxxxyyzz_1,   \
                             g_0_yzzz_0_xxxxyzzz_0,   \
                             g_0_yzzz_0_xxxxyzzz_1,   \
                             g_0_yzzz_0_xxxxzzzz_0,   \
                             g_0_yzzz_0_xxxxzzzz_1,   \
                             g_0_yzzz_0_xxxyyyyz_0,   \
                             g_0_yzzz_0_xxxyyyyz_1,   \
                             g_0_yzzz_0_xxxyyyzz_0,   \
                             g_0_yzzz_0_xxxyyyzz_1,   \
                             g_0_yzzz_0_xxxyyzzz_0,   \
                             g_0_yzzz_0_xxxyyzzz_1,   \
                             g_0_yzzz_0_xxxyzzzz_0,   \
                             g_0_yzzz_0_xxxyzzzz_1,   \
                             g_0_yzzz_0_xxxzzzzz_0,   \
                             g_0_yzzz_0_xxxzzzzz_1,   \
                             g_0_yzzz_0_xxyyyyyz_0,   \
                             g_0_yzzz_0_xxyyyyyz_1,   \
                             g_0_yzzz_0_xxyyyyzz_0,   \
                             g_0_yzzz_0_xxyyyyzz_1,   \
                             g_0_yzzz_0_xxyyyzzz_0,   \
                             g_0_yzzz_0_xxyyyzzz_1,   \
                             g_0_yzzz_0_xxyyzzzz_0,   \
                             g_0_yzzz_0_xxyyzzzz_1,   \
                             g_0_yzzz_0_xxyzzzzz_0,   \
                             g_0_yzzz_0_xxyzzzzz_1,   \
                             g_0_yzzz_0_xxzzzzzz_0,   \
                             g_0_yzzz_0_xxzzzzzz_1,   \
                             g_0_yzzz_0_xyyyyyyz_0,   \
                             g_0_yzzz_0_xyyyyyyz_1,   \
                             g_0_yzzz_0_xyyyyyzz_0,   \
                             g_0_yzzz_0_xyyyyyzz_1,   \
                             g_0_yzzz_0_xyyyyzzz_0,   \
                             g_0_yzzz_0_xyyyyzzz_1,   \
                             g_0_yzzz_0_xyyyzzzz_0,   \
                             g_0_yzzz_0_xyyyzzzz_1,   \
                             g_0_yzzz_0_xyyzzzzz_0,   \
                             g_0_yzzz_0_xyyzzzzz_1,   \
                             g_0_yzzz_0_xyzzzzzz_0,   \
                             g_0_yzzz_0_xyzzzzzz_1,   \
                             g_0_yzzz_0_xzzzzzzz_0,   \
                             g_0_yzzz_0_xzzzzzzz_1,   \
                             g_0_yzzz_0_yyyyyyyz_0,   \
                             g_0_yzzz_0_yyyyyyyz_1,   \
                             g_0_yzzz_0_yyyyyyzz_0,   \
                             g_0_yzzz_0_yyyyyyzz_1,   \
                             g_0_yzzz_0_yyyyyzzz_0,   \
                             g_0_yzzz_0_yyyyyzzz_1,   \
                             g_0_yzzz_0_yyyyzzzz_0,   \
                             g_0_yzzz_0_yyyyzzzz_1,   \
                             g_0_yzzz_0_yyyzzzzz_0,   \
                             g_0_yzzz_0_yyyzzzzz_1,   \
                             g_0_yzzz_0_yyzzzzzz_0,   \
                             g_0_yzzz_0_yyzzzzzz_1,   \
                             g_0_yzzz_0_yzzzzzzz_0,   \
                             g_0_yzzz_0_yzzzzzzz_1,   \
                             g_0_yzzz_0_zzzzzzzz_0,   \
                             g_0_yzzz_0_zzzzzzzz_1,   \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzz_0_xxxxxxxx_0[i] = 2.0 * g_0_yzzz_0_xxxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxxxx_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxxxxxx_0[i] * pb_y + g_0_yyzzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxxxxy_0[i] = 2.0 * g_0_yyyz_0_xxxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxxxxxy_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxxxxxy_0[i] * pb_z + g_0_yyyzz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxxxxxz_0[i] = 2.0 * g_0_yzzz_0_xxxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxxxz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxxxxxz_0[i] * pb_y + g_0_yyzzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxxxyy_0[i] = 2.0 * g_0_yyyz_0_xxxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxxxxyy_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxxxxyy_0[i] * pb_z + g_0_yyyzz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxxxxyz_0[i] = 2.0 * g_0_yzzz_0_xxxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxxyz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxxyz_0[i] * pb_y + g_0_yyzzz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxxxzz_0[i] = 2.0 * g_0_yzzz_0_xxxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxxzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxxxxzz_0[i] * pb_y + g_0_yyzzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxxyyy_0[i] = 2.0 * g_0_yyyz_0_xxxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxxxyyy_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxxxyyy_0[i] * pb_z + g_0_yyyzz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxxxyyz_0[i] = 2.0 * g_0_yzzz_0_xxxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxyyz_0[i] * pb_y +
                                     g_0_yyzzz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxxyzz_0[i] = 2.0 * g_0_yzzz_0_xxxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxyzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxyzz_0[i] * pb_y + g_0_yyzzz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxxzzz_0[i] = 2.0 * g_0_yzzz_0_xxxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxxxzzz_0[i] * pb_y + g_0_yyzzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxyyyy_0[i] = 2.0 * g_0_yyyz_0_xxxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxxyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxxyyyy_0[i] * pb_z + g_0_yyyzz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxxyyyz_0[i] = 2.0 * g_0_yzzz_0_xxxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyyyz_0[i] * pb_y +
                                     g_0_yyzzz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxyyzz_0[i] = 2.0 * g_0_yzzz_0_xxxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyyzz_0[i] * pb_y +
                                     g_0_yyzzz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxyzzz_0[i] = 2.0 * g_0_yzzz_0_xxxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxyzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyzzz_0[i] * pb_y + g_0_yyzzz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxzzzz_0[i] = 2.0 * g_0_yzzz_0_xxxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxxzzzz_0[i] * pb_y + g_0_yyzzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxyyyyy_0[i] = 2.0 * g_0_yyyz_0_xxxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxxyyyyy_0[i] * pb_z + g_0_yyyzz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxyyyyz_0[i] = 2.0 * g_0_yzzz_0_xxxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyyyz_0[i] * pb_y +
                                     g_0_yyzzz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxyyyzz_0[i] = 2.0 * g_0_yzzz_0_xxxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyyzz_0[i] * pb_y +
                                     g_0_yyzzz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxyyzzz_0[i] = 2.0 * g_0_yzzz_0_xxxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyzzz_0[i] * pb_y +
                                     g_0_yyzzz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxyzzzz_0[i] = 2.0 * g_0_yzzz_0_xxxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxyzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyzzzz_0[i] * pb_y + g_0_yyzzz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxzzzzz_0[i] = 2.0 * g_0_yzzz_0_xxxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxzzzzz_0[i] * pb_y + g_0_yyzzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyyyyyy_0[i] = 2.0 * g_0_yyyz_0_xxyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxyyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xxyyyyyy_0[i] * pb_z + g_0_yyyzz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxyyyyyz_0[i] = 2.0 * g_0_yzzz_0_xxyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyyyz_0[i] * pb_y +
                                     g_0_yyzzz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyyyyzz_0[i] = 2.0 * g_0_yzzz_0_xxyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyyzz_0[i] * pb_y +
                                     g_0_yyzzz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyyyzzz_0[i] = 2.0 * g_0_yzzz_0_xxyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyzzz_0[i] * pb_y +
                                     g_0_yyzzz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyyzzzz_0[i] = 2.0 * g_0_yzzz_0_xxyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyzzzz_0[i] * pb_y +
                                     g_0_yyzzz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyzzzzz_0[i] = 2.0 * g_0_yzzz_0_xxyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyzzzzz_0[i] * pb_y + g_0_yyzzz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxzzzzzz_0[i] = 2.0 * g_0_yzzz_0_xxzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxzzzzzz_0[i] * pb_y + g_0_yyzzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyyyyyy_0[i] = 2.0 * g_0_yyyz_0_xyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_xyyyyyyy_0[i] * pb_z + g_0_yyyzz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xyyyyyyz_0[i] = 2.0 * g_0_yzzz_0_xyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyyyyyz_1[i] * fti_ab_0 +
                                     6.0 * g_0_yyzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyyyz_0[i] * pb_y +
                                     g_0_yyzzz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyyyyzz_0[i] = 2.0 * g_0_yzzz_0_xyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyyyyzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyyzz_0[i] * pb_y +
                                     g_0_yyzzz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyyyzzz_0[i] = 2.0 * g_0_yzzz_0_xyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyyyzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyzzz_0[i] * pb_y +
                                     g_0_yyzzz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyyzzzz_0[i] = 2.0 * g_0_yzzz_0_xyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyyzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyzzzz_0[i] * pb_y +
                                     g_0_yyzzz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyzzzzz_0[i] = 2.0 * g_0_yzzz_0_xyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyzzzzz_0[i] * pb_y +
                                     g_0_yyzzz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyzzzzzz_0[i] = 2.0 * g_0_yzzz_0_xyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyzzzzzz_0[i] * pb_y + g_0_yyzzz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xzzzzzzz_0[i] = 2.0 * g_0_yzzz_0_xzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xzzzzzzz_0[i] * pb_y + g_0_yyzzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyyyyyy_0[i] = 2.0 * g_0_yyyz_0_yyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_yyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyzz_0_yyyyyyyy_0[i] * pb_z + g_0_yyyzz_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_yyyyyyyz_0[i] = 2.0 * g_0_yzzz_0_yyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyyyyyz_1[i] * fti_ab_0 +
                                     7.0 * g_0_yyzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyyyyyz_0[i] * pb_y +
                                     g_0_yyzzz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyyyyzz_0[i] = 2.0 * g_0_yzzz_0_yyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyyyyzz_1[i] * fti_ab_0 +
                                     6.0 * g_0_yyzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyyyyzz_0[i] * pb_y +
                                     g_0_yyzzz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyyyzzz_0[i] = 2.0 * g_0_yzzz_0_yyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyyyzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyyyzzz_0[i] * pb_y +
                                     g_0_yyzzz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyyzzzz_0[i] = 2.0 * g_0_yzzz_0_yyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyyzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyyzzzz_0[i] * pb_y +
                                     g_0_yyzzz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyzzzzz_0[i] = 2.0 * g_0_yzzz_0_yyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyzzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyzzzzz_0[i] * pb_y +
                                     g_0_yyzzz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyzzzzzz_0[i] = 2.0 * g_0_yzzz_0_yyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyzzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyzzzzzz_0[i] * pb_y +
                                     g_0_yyzzz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yzzzzzzz_0[i] = 2.0 * g_0_yzzz_0_yzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yzzzzzzz_0[i] * pb_y + g_0_yyzzz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_zzzzzzzz_0[i] = 2.0 * g_0_yzzz_0_zzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_zzzzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_zzzzzzzz_0[i] * pb_y + g_0_yyzzz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1125-1170 components of targeted buffer : SISL

    auto g_0_yyzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 1125);

    auto g_0_yyzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 1126);

    auto g_0_yyzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 1127);

    auto g_0_yyzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 1128);

    auto g_0_yyzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 1129);

    auto g_0_yyzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 1130);

    auto g_0_yyzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 1131);

    auto g_0_yyzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 1132);

    auto g_0_yyzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 1133);

    auto g_0_yyzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 1134);

    auto g_0_yyzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1135);

    auto g_0_yyzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1136);

    auto g_0_yyzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1137);

    auto g_0_yyzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1138);

    auto g_0_yyzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1139);

    auto g_0_yyzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1140);

    auto g_0_yyzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1141);

    auto g_0_yyzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1142);

    auto g_0_yyzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1143);

    auto g_0_yyzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1144);

    auto g_0_yyzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1145);

    auto g_0_yyzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1146);

    auto g_0_yyzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1147);

    auto g_0_yyzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1148);

    auto g_0_yyzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1149);

    auto g_0_yyzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1150);

    auto g_0_yyzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1151);

    auto g_0_yyzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1152);

    auto g_0_yyzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1153);

    auto g_0_yyzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1154);

    auto g_0_yyzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1155);

    auto g_0_yyzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1156);

    auto g_0_yyzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1157);

    auto g_0_yyzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1158);

    auto g_0_yyzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1159);

    auto g_0_yyzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1160);

    auto g_0_yyzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1161);

    auto g_0_yyzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1162);

    auto g_0_yyzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1163);

    auto g_0_yyzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1164);

    auto g_0_yyzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1165);

    auto g_0_yyzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1166);

    auto g_0_yyzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1167);

    auto g_0_yyzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1168);

    auto g_0_yyzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1169);

#pragma omp simd aligned(g_0_yyzz_0_xxxxxxxy_0,       \
                             g_0_yyzz_0_xxxxxxxy_1,   \
                             g_0_yyzz_0_xxxxxxyy_0,   \
                             g_0_yyzz_0_xxxxxxyy_1,   \
                             g_0_yyzz_0_xxxxxyyy_0,   \
                             g_0_yyzz_0_xxxxxyyy_1,   \
                             g_0_yyzz_0_xxxxyyyy_0,   \
                             g_0_yyzz_0_xxxxyyyy_1,   \
                             g_0_yyzz_0_xxxyyyyy_0,   \
                             g_0_yyzz_0_xxxyyyyy_1,   \
                             g_0_yyzz_0_xxyyyyyy_0,   \
                             g_0_yyzz_0_xxyyyyyy_1,   \
                             g_0_yyzz_0_xyyyyyyy_0,   \
                             g_0_yyzz_0_xyyyyyyy_1,   \
                             g_0_yyzz_0_yyyyyyyy_0,   \
                             g_0_yyzz_0_yyyyyyyy_1,   \
                             g_0_yyzzz_0_xxxxxxxy_0,  \
                             g_0_yyzzz_0_xxxxxxxy_1,  \
                             g_0_yyzzz_0_xxxxxxyy_0,  \
                             g_0_yyzzz_0_xxxxxxyy_1,  \
                             g_0_yyzzz_0_xxxxxyyy_0,  \
                             g_0_yyzzz_0_xxxxxyyy_1,  \
                             g_0_yyzzz_0_xxxxyyyy_0,  \
                             g_0_yyzzz_0_xxxxyyyy_1,  \
                             g_0_yyzzz_0_xxxyyyyy_0,  \
                             g_0_yyzzz_0_xxxyyyyy_1,  \
                             g_0_yyzzz_0_xxyyyyyy_0,  \
                             g_0_yyzzz_0_xxyyyyyy_1,  \
                             g_0_yyzzz_0_xyyyyyyy_0,  \
                             g_0_yyzzz_0_xyyyyyyy_1,  \
                             g_0_yyzzz_0_yyyyyyyy_0,  \
                             g_0_yyzzz_0_yyyyyyyy_1,  \
                             g_0_yyzzzz_0_xxxxxxxx_0, \
                             g_0_yyzzzz_0_xxxxxxxy_0, \
                             g_0_yyzzzz_0_xxxxxxxz_0, \
                             g_0_yyzzzz_0_xxxxxxyy_0, \
                             g_0_yyzzzz_0_xxxxxxyz_0, \
                             g_0_yyzzzz_0_xxxxxxzz_0, \
                             g_0_yyzzzz_0_xxxxxyyy_0, \
                             g_0_yyzzzz_0_xxxxxyyz_0, \
                             g_0_yyzzzz_0_xxxxxyzz_0, \
                             g_0_yyzzzz_0_xxxxxzzz_0, \
                             g_0_yyzzzz_0_xxxxyyyy_0, \
                             g_0_yyzzzz_0_xxxxyyyz_0, \
                             g_0_yyzzzz_0_xxxxyyzz_0, \
                             g_0_yyzzzz_0_xxxxyzzz_0, \
                             g_0_yyzzzz_0_xxxxzzzz_0, \
                             g_0_yyzzzz_0_xxxyyyyy_0, \
                             g_0_yyzzzz_0_xxxyyyyz_0, \
                             g_0_yyzzzz_0_xxxyyyzz_0, \
                             g_0_yyzzzz_0_xxxyyzzz_0, \
                             g_0_yyzzzz_0_xxxyzzzz_0, \
                             g_0_yyzzzz_0_xxxzzzzz_0, \
                             g_0_yyzzzz_0_xxyyyyyy_0, \
                             g_0_yyzzzz_0_xxyyyyyz_0, \
                             g_0_yyzzzz_0_xxyyyyzz_0, \
                             g_0_yyzzzz_0_xxyyyzzz_0, \
                             g_0_yyzzzz_0_xxyyzzzz_0, \
                             g_0_yyzzzz_0_xxyzzzzz_0, \
                             g_0_yyzzzz_0_xxzzzzzz_0, \
                             g_0_yyzzzz_0_xyyyyyyy_0, \
                             g_0_yyzzzz_0_xyyyyyyz_0, \
                             g_0_yyzzzz_0_xyyyyyzz_0, \
                             g_0_yyzzzz_0_xyyyyzzz_0, \
                             g_0_yyzzzz_0_xyyyzzzz_0, \
                             g_0_yyzzzz_0_xyyzzzzz_0, \
                             g_0_yyzzzz_0_xyzzzzzz_0, \
                             g_0_yyzzzz_0_xzzzzzzz_0, \
                             g_0_yyzzzz_0_yyyyyyyy_0, \
                             g_0_yyzzzz_0_yyyyyyyz_0, \
                             g_0_yyzzzz_0_yyyyyyzz_0, \
                             g_0_yyzzzz_0_yyyyyzzz_0, \
                             g_0_yyzzzz_0_yyyyzzzz_0, \
                             g_0_yyzzzz_0_yyyzzzzz_0, \
                             g_0_yyzzzz_0_yyzzzzzz_0, \
                             g_0_yyzzzz_0_yzzzzzzz_0, \
                             g_0_yyzzzz_0_zzzzzzzz_0, \
                             g_0_yzzzz_0_xxxxxxxx_0,  \
                             g_0_yzzzz_0_xxxxxxxx_1,  \
                             g_0_yzzzz_0_xxxxxxxz_0,  \
                             g_0_yzzzz_0_xxxxxxxz_1,  \
                             g_0_yzzzz_0_xxxxxxyz_0,  \
                             g_0_yzzzz_0_xxxxxxyz_1,  \
                             g_0_yzzzz_0_xxxxxxz_1,   \
                             g_0_yzzzz_0_xxxxxxzz_0,  \
                             g_0_yzzzz_0_xxxxxxzz_1,  \
                             g_0_yzzzz_0_xxxxxyyz_0,  \
                             g_0_yzzzz_0_xxxxxyyz_1,  \
                             g_0_yzzzz_0_xxxxxyz_1,   \
                             g_0_yzzzz_0_xxxxxyzz_0,  \
                             g_0_yzzzz_0_xxxxxyzz_1,  \
                             g_0_yzzzz_0_xxxxxzz_1,   \
                             g_0_yzzzz_0_xxxxxzzz_0,  \
                             g_0_yzzzz_0_xxxxxzzz_1,  \
                             g_0_yzzzz_0_xxxxyyyz_0,  \
                             g_0_yzzzz_0_xxxxyyyz_1,  \
                             g_0_yzzzz_0_xxxxyyz_1,   \
                             g_0_yzzzz_0_xxxxyyzz_0,  \
                             g_0_yzzzz_0_xxxxyyzz_1,  \
                             g_0_yzzzz_0_xxxxyzz_1,   \
                             g_0_yzzzz_0_xxxxyzzz_0,  \
                             g_0_yzzzz_0_xxxxyzzz_1,  \
                             g_0_yzzzz_0_xxxxzzz_1,   \
                             g_0_yzzzz_0_xxxxzzzz_0,  \
                             g_0_yzzzz_0_xxxxzzzz_1,  \
                             g_0_yzzzz_0_xxxyyyyz_0,  \
                             g_0_yzzzz_0_xxxyyyyz_1,  \
                             g_0_yzzzz_0_xxxyyyz_1,   \
                             g_0_yzzzz_0_xxxyyyzz_0,  \
                             g_0_yzzzz_0_xxxyyyzz_1,  \
                             g_0_yzzzz_0_xxxyyzz_1,   \
                             g_0_yzzzz_0_xxxyyzzz_0,  \
                             g_0_yzzzz_0_xxxyyzzz_1,  \
                             g_0_yzzzz_0_xxxyzzz_1,   \
                             g_0_yzzzz_0_xxxyzzzz_0,  \
                             g_0_yzzzz_0_xxxyzzzz_1,  \
                             g_0_yzzzz_0_xxxzzzz_1,   \
                             g_0_yzzzz_0_xxxzzzzz_0,  \
                             g_0_yzzzz_0_xxxzzzzz_1,  \
                             g_0_yzzzz_0_xxyyyyyz_0,  \
                             g_0_yzzzz_0_xxyyyyyz_1,  \
                             g_0_yzzzz_0_xxyyyyz_1,   \
                             g_0_yzzzz_0_xxyyyyzz_0,  \
                             g_0_yzzzz_0_xxyyyyzz_1,  \
                             g_0_yzzzz_0_xxyyyzz_1,   \
                             g_0_yzzzz_0_xxyyyzzz_0,  \
                             g_0_yzzzz_0_xxyyyzzz_1,  \
                             g_0_yzzzz_0_xxyyzzz_1,   \
                             g_0_yzzzz_0_xxyyzzzz_0,  \
                             g_0_yzzzz_0_xxyyzzzz_1,  \
                             g_0_yzzzz_0_xxyzzzz_1,   \
                             g_0_yzzzz_0_xxyzzzzz_0,  \
                             g_0_yzzzz_0_xxyzzzzz_1,  \
                             g_0_yzzzz_0_xxzzzzz_1,   \
                             g_0_yzzzz_0_xxzzzzzz_0,  \
                             g_0_yzzzz_0_xxzzzzzz_1,  \
                             g_0_yzzzz_0_xyyyyyyz_0,  \
                             g_0_yzzzz_0_xyyyyyyz_1,  \
                             g_0_yzzzz_0_xyyyyyz_1,   \
                             g_0_yzzzz_0_xyyyyyzz_0,  \
                             g_0_yzzzz_0_xyyyyyzz_1,  \
                             g_0_yzzzz_0_xyyyyzz_1,   \
                             g_0_yzzzz_0_xyyyyzzz_0,  \
                             g_0_yzzzz_0_xyyyyzzz_1,  \
                             g_0_yzzzz_0_xyyyzzz_1,   \
                             g_0_yzzzz_0_xyyyzzzz_0,  \
                             g_0_yzzzz_0_xyyyzzzz_1,  \
                             g_0_yzzzz_0_xyyzzzz_1,   \
                             g_0_yzzzz_0_xyyzzzzz_0,  \
                             g_0_yzzzz_0_xyyzzzzz_1,  \
                             g_0_yzzzz_0_xyzzzzz_1,   \
                             g_0_yzzzz_0_xyzzzzzz_0,  \
                             g_0_yzzzz_0_xyzzzzzz_1,  \
                             g_0_yzzzz_0_xzzzzzz_1,   \
                             g_0_yzzzz_0_xzzzzzzz_0,  \
                             g_0_yzzzz_0_xzzzzzzz_1,  \
                             g_0_yzzzz_0_yyyyyyyz_0,  \
                             g_0_yzzzz_0_yyyyyyyz_1,  \
                             g_0_yzzzz_0_yyyyyyz_1,   \
                             g_0_yzzzz_0_yyyyyyzz_0,  \
                             g_0_yzzzz_0_yyyyyyzz_1,  \
                             g_0_yzzzz_0_yyyyyzz_1,   \
                             g_0_yzzzz_0_yyyyyzzz_0,  \
                             g_0_yzzzz_0_yyyyyzzz_1,  \
                             g_0_yzzzz_0_yyyyzzz_1,   \
                             g_0_yzzzz_0_yyyyzzzz_0,  \
                             g_0_yzzzz_0_yyyyzzzz_1,  \
                             g_0_yzzzz_0_yyyzzzz_1,   \
                             g_0_yzzzz_0_yyyzzzzz_0,  \
                             g_0_yzzzz_0_yyyzzzzz_1,  \
                             g_0_yzzzz_0_yyzzzzz_1,   \
                             g_0_yzzzz_0_yyzzzzzz_0,  \
                             g_0_yzzzz_0_yyzzzzzz_1,  \
                             g_0_yzzzz_0_yzzzzzz_1,   \
                             g_0_yzzzz_0_yzzzzzzz_0,  \
                             g_0_yzzzz_0_yzzzzzzz_1,  \
                             g_0_yzzzz_0_zzzzzzz_1,   \
                             g_0_yzzzz_0_zzzzzzzz_0,  \
                             g_0_yzzzz_0_zzzzzzzz_1,  \
                             g_0_zzzz_0_xxxxxxxx_0,   \
                             g_0_zzzz_0_xxxxxxxx_1,   \
                             g_0_zzzz_0_xxxxxxxz_0,   \
                             g_0_zzzz_0_xxxxxxxz_1,   \
                             g_0_zzzz_0_xxxxxxyz_0,   \
                             g_0_zzzz_0_xxxxxxyz_1,   \
                             g_0_zzzz_0_xxxxxxzz_0,   \
                             g_0_zzzz_0_xxxxxxzz_1,   \
                             g_0_zzzz_0_xxxxxyyz_0,   \
                             g_0_zzzz_0_xxxxxyyz_1,   \
                             g_0_zzzz_0_xxxxxyzz_0,   \
                             g_0_zzzz_0_xxxxxyzz_1,   \
                             g_0_zzzz_0_xxxxxzzz_0,   \
                             g_0_zzzz_0_xxxxxzzz_1,   \
                             g_0_zzzz_0_xxxxyyyz_0,   \
                             g_0_zzzz_0_xxxxyyyz_1,   \
                             g_0_zzzz_0_xxxxyyzz_0,   \
                             g_0_zzzz_0_xxxxyyzz_1,   \
                             g_0_zzzz_0_xxxxyzzz_0,   \
                             g_0_zzzz_0_xxxxyzzz_1,   \
                             g_0_zzzz_0_xxxxzzzz_0,   \
                             g_0_zzzz_0_xxxxzzzz_1,   \
                             g_0_zzzz_0_xxxyyyyz_0,   \
                             g_0_zzzz_0_xxxyyyyz_1,   \
                             g_0_zzzz_0_xxxyyyzz_0,   \
                             g_0_zzzz_0_xxxyyyzz_1,   \
                             g_0_zzzz_0_xxxyyzzz_0,   \
                             g_0_zzzz_0_xxxyyzzz_1,   \
                             g_0_zzzz_0_xxxyzzzz_0,   \
                             g_0_zzzz_0_xxxyzzzz_1,   \
                             g_0_zzzz_0_xxxzzzzz_0,   \
                             g_0_zzzz_0_xxxzzzzz_1,   \
                             g_0_zzzz_0_xxyyyyyz_0,   \
                             g_0_zzzz_0_xxyyyyyz_1,   \
                             g_0_zzzz_0_xxyyyyzz_0,   \
                             g_0_zzzz_0_xxyyyyzz_1,   \
                             g_0_zzzz_0_xxyyyzzz_0,   \
                             g_0_zzzz_0_xxyyyzzz_1,   \
                             g_0_zzzz_0_xxyyzzzz_0,   \
                             g_0_zzzz_0_xxyyzzzz_1,   \
                             g_0_zzzz_0_xxyzzzzz_0,   \
                             g_0_zzzz_0_xxyzzzzz_1,   \
                             g_0_zzzz_0_xxzzzzzz_0,   \
                             g_0_zzzz_0_xxzzzzzz_1,   \
                             g_0_zzzz_0_xyyyyyyz_0,   \
                             g_0_zzzz_0_xyyyyyyz_1,   \
                             g_0_zzzz_0_xyyyyyzz_0,   \
                             g_0_zzzz_0_xyyyyyzz_1,   \
                             g_0_zzzz_0_xyyyyzzz_0,   \
                             g_0_zzzz_0_xyyyyzzz_1,   \
                             g_0_zzzz_0_xyyyzzzz_0,   \
                             g_0_zzzz_0_xyyyzzzz_1,   \
                             g_0_zzzz_0_xyyzzzzz_0,   \
                             g_0_zzzz_0_xyyzzzzz_1,   \
                             g_0_zzzz_0_xyzzzzzz_0,   \
                             g_0_zzzz_0_xyzzzzzz_1,   \
                             g_0_zzzz_0_xzzzzzzz_0,   \
                             g_0_zzzz_0_xzzzzzzz_1,   \
                             g_0_zzzz_0_yyyyyyyz_0,   \
                             g_0_zzzz_0_yyyyyyyz_1,   \
                             g_0_zzzz_0_yyyyyyzz_0,   \
                             g_0_zzzz_0_yyyyyyzz_1,   \
                             g_0_zzzz_0_yyyyyzzz_0,   \
                             g_0_zzzz_0_yyyyyzzz_1,   \
                             g_0_zzzz_0_yyyyzzzz_0,   \
                             g_0_zzzz_0_yyyyzzzz_1,   \
                             g_0_zzzz_0_yyyzzzzz_0,   \
                             g_0_zzzz_0_yyyzzzzz_1,   \
                             g_0_zzzz_0_yyzzzzzz_0,   \
                             g_0_zzzz_0_yyzzzzzz_1,   \
                             g_0_zzzz_0_yzzzzzzz_0,   \
                             g_0_zzzz_0_yzzzzzzz_1,   \
                             g_0_zzzz_0_zzzzzzzz_0,   \
                             g_0_zzzz_0_zzzzzzzz_1,   \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzz_0_xxxxxxxx_0[i] = g_0_zzzz_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxxxxx_0[i] * pb_y +
                                     g_0_yzzzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxxxxy_0[i] = 3.0 * g_0_yyzz_0_xxxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxxxy_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxxxxxy_0[i] * pb_z + g_0_yyzzz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxxxxxz_0[i] = g_0_zzzz_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxxxxz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxxxyy_0[i] = 3.0 * g_0_yyzz_0_xxxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxxyy_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxxxxyy_0[i] * pb_z + g_0_yyzzz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxxxxyz_0[i] = g_0_zzzz_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxxyz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxxxz_1[i] * fi_abcd_0 +
                                     g_0_yzzzz_0_xxxxxxyz_0[i] * pb_y + g_0_yzzzz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxxxzz_0[i] = g_0_zzzz_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxxxzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxxyyy_0[i] = 3.0 * g_0_yyzz_0_xxxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxyyy_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxxxyyy_0[i] * pb_z + g_0_yyzzz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxxxyyz_0[i] = g_0_zzzz_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yzzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxxyyz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxxyzz_0[i] = g_0_zzzz_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxyzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxxzz_1[i] * fi_abcd_0 +
                                     g_0_yzzzz_0_xxxxxyzz_0[i] * pb_y + g_0_yzzzz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxxzzz_0[i] = g_0_zzzz_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxxzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxyyyy_0[i] = 3.0 * g_0_yyzz_0_xxxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxyyyy_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxxyyyy_0[i] * pb_z + g_0_yyzzz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxxyyyz_0[i] = g_0_zzzz_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yzzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyyyz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxyyzz_0[i] = g_0_zzzz_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yzzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyyzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxyzzz_0[i] = g_0_zzzz_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxyzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxzzz_1[i] * fi_abcd_0 +
                                     g_0_yzzzz_0_xxxxyzzz_0[i] * pb_y + g_0_yzzzz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxzzzz_0[i] = g_0_zzzz_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxzzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxyyyyy_0[i] = 3.0 * g_0_yyzz_0_xxxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxxyyyyy_0[i] * pb_z + g_0_yyzzz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxyyyyz_0[i] = g_0_zzzz_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yzzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyyyz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxyyyzz_0[i] = g_0_zzzz_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yzzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyyzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxyyzzz_0[i] = g_0_zzzz_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yzzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxyzzzz_0[i] = g_0_zzzz_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxzzzz_1[i] * fi_abcd_0 +
                                     g_0_yzzzz_0_xxxyzzzz_0[i] * pb_y + g_0_yzzzz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxzzzzz_0[i] = g_0_zzzz_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxzzzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyyyyyy_0[i] = 3.0 * g_0_yyzz_0_xxyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xxyyyyyy_0[i] * pb_z + g_0_yyzzz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxyyyyyz_0[i] = g_0_zzzz_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yzzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyyyz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyyyyzz_0[i] = g_0_zzzz_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yzzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyyzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyyyzzz_0[i] = g_0_zzzz_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yzzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyyzzzz_0[i] = g_0_zzzz_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yzzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyzzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyzzzzz_0[i] = g_0_zzzz_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxzzzzz_1[i] * fi_abcd_0 +
                                     g_0_yzzzz_0_xxyzzzzz_0[i] * pb_y + g_0_yzzzz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxzzzzzz_0[i] = g_0_zzzz_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxzzzzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyyyyyy_0[i] = 3.0 * g_0_yyzz_0_xyyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_xyyyyyyy_0[i] * pb_z + g_0_yyzzz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xyyyyyyz_0[i] = g_0_zzzz_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyyyyz_1[i] * fti_ab_0 +
                                     6.0 * g_0_yzzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyyyz_0[i] * pb_y +
                                     g_0_yzzzz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyyyyzz_0[i] = g_0_zzzz_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyyyzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yzzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyyzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyyyzzz_0[i] = g_0_zzzz_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyyzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yzzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyyzzzz_0[i] = g_0_zzzz_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yzzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyzzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyzzzzz_0[i] = g_0_zzzz_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yzzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyzzzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyzzzzzz_0[i] = g_0_zzzz_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xzzzzzz_1[i] * fi_abcd_0 +
                                     g_0_yzzzz_0_xyzzzzzz_0[i] * pb_y + g_0_yzzzz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xzzzzzzz_0[i] = g_0_zzzz_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xzzzzzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyyyyyy_0[i] = 3.0 * g_0_yyzz_0_yyyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyzzz_0_yyyyyyyy_0[i] * pb_z + g_0_yyzzz_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_yyyyyyyz_0[i] = g_0_zzzz_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyyyz_1[i] * fti_ab_0 +
                                     7.0 * g_0_yzzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yyyyyyyz_0[i] * pb_y +
                                     g_0_yzzzz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyyyyzz_0[i] = g_0_zzzz_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyyzz_1[i] * fti_ab_0 +
                                     6.0 * g_0_yzzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yyyyyyzz_0[i] * pb_y +
                                     g_0_yzzzz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyyyzzz_0[i] = g_0_zzzz_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yzzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yyyyyzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyyzzzz_0[i] = g_0_zzzz_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yzzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yyyyzzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyzzzzz_0[i] = g_0_zzzz_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyzzzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yzzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yyyzzzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyzzzzzz_0[i] = g_0_zzzz_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyzzzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yzzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yyzzzzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yzzzzzzz_0[i] = g_0_zzzz_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_zzzzzzz_1[i] * fi_abcd_0 +
                                     g_0_yzzzz_0_yzzzzzzz_0[i] * pb_y + g_0_yzzzz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_zzzzzzzz_0[i] = g_0_zzzz_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_zzzzzzzz_0[i] * pb_y +
                                     g_0_yzzzz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1170-1215 components of targeted buffer : SISL

    auto g_0_yzzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 1170);

    auto g_0_yzzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 1171);

    auto g_0_yzzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 1172);

    auto g_0_yzzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 1173);

    auto g_0_yzzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 1174);

    auto g_0_yzzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 1175);

    auto g_0_yzzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 1176);

    auto g_0_yzzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 1177);

    auto g_0_yzzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 1178);

    auto g_0_yzzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 1179);

    auto g_0_yzzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1180);

    auto g_0_yzzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1181);

    auto g_0_yzzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1182);

    auto g_0_yzzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1183);

    auto g_0_yzzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1184);

    auto g_0_yzzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1185);

    auto g_0_yzzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1186);

    auto g_0_yzzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1187);

    auto g_0_yzzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1188);

    auto g_0_yzzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1189);

    auto g_0_yzzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1190);

    auto g_0_yzzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1191);

    auto g_0_yzzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1192);

    auto g_0_yzzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1193);

    auto g_0_yzzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1194);

    auto g_0_yzzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1195);

    auto g_0_yzzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1196);

    auto g_0_yzzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1197);

    auto g_0_yzzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1198);

    auto g_0_yzzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1199);

    auto g_0_yzzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1200);

    auto g_0_yzzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1201);

    auto g_0_yzzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1202);

    auto g_0_yzzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1203);

    auto g_0_yzzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1204);

    auto g_0_yzzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1205);

    auto g_0_yzzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1206);

    auto g_0_yzzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1207);

    auto g_0_yzzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1208);

    auto g_0_yzzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1209);

    auto g_0_yzzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1210);

    auto g_0_yzzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1211);

    auto g_0_yzzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1212);

    auto g_0_yzzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1213);

    auto g_0_yzzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1214);

#pragma omp simd aligned(g_0_yzzzzz_0_xxxxxxxx_0,     \
                             g_0_yzzzzz_0_xxxxxxxy_0, \
                             g_0_yzzzzz_0_xxxxxxxz_0, \
                             g_0_yzzzzz_0_xxxxxxyy_0, \
                             g_0_yzzzzz_0_xxxxxxyz_0, \
                             g_0_yzzzzz_0_xxxxxxzz_0, \
                             g_0_yzzzzz_0_xxxxxyyy_0, \
                             g_0_yzzzzz_0_xxxxxyyz_0, \
                             g_0_yzzzzz_0_xxxxxyzz_0, \
                             g_0_yzzzzz_0_xxxxxzzz_0, \
                             g_0_yzzzzz_0_xxxxyyyy_0, \
                             g_0_yzzzzz_0_xxxxyyyz_0, \
                             g_0_yzzzzz_0_xxxxyyzz_0, \
                             g_0_yzzzzz_0_xxxxyzzz_0, \
                             g_0_yzzzzz_0_xxxxzzzz_0, \
                             g_0_yzzzzz_0_xxxyyyyy_0, \
                             g_0_yzzzzz_0_xxxyyyyz_0, \
                             g_0_yzzzzz_0_xxxyyyzz_0, \
                             g_0_yzzzzz_0_xxxyyzzz_0, \
                             g_0_yzzzzz_0_xxxyzzzz_0, \
                             g_0_yzzzzz_0_xxxzzzzz_0, \
                             g_0_yzzzzz_0_xxyyyyyy_0, \
                             g_0_yzzzzz_0_xxyyyyyz_0, \
                             g_0_yzzzzz_0_xxyyyyzz_0, \
                             g_0_yzzzzz_0_xxyyyzzz_0, \
                             g_0_yzzzzz_0_xxyyzzzz_0, \
                             g_0_yzzzzz_0_xxyzzzzz_0, \
                             g_0_yzzzzz_0_xxzzzzzz_0, \
                             g_0_yzzzzz_0_xyyyyyyy_0, \
                             g_0_yzzzzz_0_xyyyyyyz_0, \
                             g_0_yzzzzz_0_xyyyyyzz_0, \
                             g_0_yzzzzz_0_xyyyyzzz_0, \
                             g_0_yzzzzz_0_xyyyzzzz_0, \
                             g_0_yzzzzz_0_xyyzzzzz_0, \
                             g_0_yzzzzz_0_xyzzzzzz_0, \
                             g_0_yzzzzz_0_xzzzzzzz_0, \
                             g_0_yzzzzz_0_yyyyyyyy_0, \
                             g_0_yzzzzz_0_yyyyyyyz_0, \
                             g_0_yzzzzz_0_yyyyyyzz_0, \
                             g_0_yzzzzz_0_yyyyyzzz_0, \
                             g_0_yzzzzz_0_yyyyzzzz_0, \
                             g_0_yzzzzz_0_yyyzzzzz_0, \
                             g_0_yzzzzz_0_yyzzzzzz_0, \
                             g_0_yzzzzz_0_yzzzzzzz_0, \
                             g_0_yzzzzz_0_zzzzzzzz_0, \
                             g_0_zzzzz_0_xxxxxxx_1,   \
                             g_0_zzzzz_0_xxxxxxxx_0,  \
                             g_0_zzzzz_0_xxxxxxxx_1,  \
                             g_0_zzzzz_0_xxxxxxxy_0,  \
                             g_0_zzzzz_0_xxxxxxxy_1,  \
                             g_0_zzzzz_0_xxxxxxxz_0,  \
                             g_0_zzzzz_0_xxxxxxxz_1,  \
                             g_0_zzzzz_0_xxxxxxy_1,   \
                             g_0_zzzzz_0_xxxxxxyy_0,  \
                             g_0_zzzzz_0_xxxxxxyy_1,  \
                             g_0_zzzzz_0_xxxxxxyz_0,  \
                             g_0_zzzzz_0_xxxxxxyz_1,  \
                             g_0_zzzzz_0_xxxxxxz_1,   \
                             g_0_zzzzz_0_xxxxxxzz_0,  \
                             g_0_zzzzz_0_xxxxxxzz_1,  \
                             g_0_zzzzz_0_xxxxxyy_1,   \
                             g_0_zzzzz_0_xxxxxyyy_0,  \
                             g_0_zzzzz_0_xxxxxyyy_1,  \
                             g_0_zzzzz_0_xxxxxyyz_0,  \
                             g_0_zzzzz_0_xxxxxyyz_1,  \
                             g_0_zzzzz_0_xxxxxyz_1,   \
                             g_0_zzzzz_0_xxxxxyzz_0,  \
                             g_0_zzzzz_0_xxxxxyzz_1,  \
                             g_0_zzzzz_0_xxxxxzz_1,   \
                             g_0_zzzzz_0_xxxxxzzz_0,  \
                             g_0_zzzzz_0_xxxxxzzz_1,  \
                             g_0_zzzzz_0_xxxxyyy_1,   \
                             g_0_zzzzz_0_xxxxyyyy_0,  \
                             g_0_zzzzz_0_xxxxyyyy_1,  \
                             g_0_zzzzz_0_xxxxyyyz_0,  \
                             g_0_zzzzz_0_xxxxyyyz_1,  \
                             g_0_zzzzz_0_xxxxyyz_1,   \
                             g_0_zzzzz_0_xxxxyyzz_0,  \
                             g_0_zzzzz_0_xxxxyyzz_1,  \
                             g_0_zzzzz_0_xxxxyzz_1,   \
                             g_0_zzzzz_0_xxxxyzzz_0,  \
                             g_0_zzzzz_0_xxxxyzzz_1,  \
                             g_0_zzzzz_0_xxxxzzz_1,   \
                             g_0_zzzzz_0_xxxxzzzz_0,  \
                             g_0_zzzzz_0_xxxxzzzz_1,  \
                             g_0_zzzzz_0_xxxyyyy_1,   \
                             g_0_zzzzz_0_xxxyyyyy_0,  \
                             g_0_zzzzz_0_xxxyyyyy_1,  \
                             g_0_zzzzz_0_xxxyyyyz_0,  \
                             g_0_zzzzz_0_xxxyyyyz_1,  \
                             g_0_zzzzz_0_xxxyyyz_1,   \
                             g_0_zzzzz_0_xxxyyyzz_0,  \
                             g_0_zzzzz_0_xxxyyyzz_1,  \
                             g_0_zzzzz_0_xxxyyzz_1,   \
                             g_0_zzzzz_0_xxxyyzzz_0,  \
                             g_0_zzzzz_0_xxxyyzzz_1,  \
                             g_0_zzzzz_0_xxxyzzz_1,   \
                             g_0_zzzzz_0_xxxyzzzz_0,  \
                             g_0_zzzzz_0_xxxyzzzz_1,  \
                             g_0_zzzzz_0_xxxzzzz_1,   \
                             g_0_zzzzz_0_xxxzzzzz_0,  \
                             g_0_zzzzz_0_xxxzzzzz_1,  \
                             g_0_zzzzz_0_xxyyyyy_1,   \
                             g_0_zzzzz_0_xxyyyyyy_0,  \
                             g_0_zzzzz_0_xxyyyyyy_1,  \
                             g_0_zzzzz_0_xxyyyyyz_0,  \
                             g_0_zzzzz_0_xxyyyyyz_1,  \
                             g_0_zzzzz_0_xxyyyyz_1,   \
                             g_0_zzzzz_0_xxyyyyzz_0,  \
                             g_0_zzzzz_0_xxyyyyzz_1,  \
                             g_0_zzzzz_0_xxyyyzz_1,   \
                             g_0_zzzzz_0_xxyyyzzz_0,  \
                             g_0_zzzzz_0_xxyyyzzz_1,  \
                             g_0_zzzzz_0_xxyyzzz_1,   \
                             g_0_zzzzz_0_xxyyzzzz_0,  \
                             g_0_zzzzz_0_xxyyzzzz_1,  \
                             g_0_zzzzz_0_xxyzzzz_1,   \
                             g_0_zzzzz_0_xxyzzzzz_0,  \
                             g_0_zzzzz_0_xxyzzzzz_1,  \
                             g_0_zzzzz_0_xxzzzzz_1,   \
                             g_0_zzzzz_0_xxzzzzzz_0,  \
                             g_0_zzzzz_0_xxzzzzzz_1,  \
                             g_0_zzzzz_0_xyyyyyy_1,   \
                             g_0_zzzzz_0_xyyyyyyy_0,  \
                             g_0_zzzzz_0_xyyyyyyy_1,  \
                             g_0_zzzzz_0_xyyyyyyz_0,  \
                             g_0_zzzzz_0_xyyyyyyz_1,  \
                             g_0_zzzzz_0_xyyyyyz_1,   \
                             g_0_zzzzz_0_xyyyyyzz_0,  \
                             g_0_zzzzz_0_xyyyyyzz_1,  \
                             g_0_zzzzz_0_xyyyyzz_1,   \
                             g_0_zzzzz_0_xyyyyzzz_0,  \
                             g_0_zzzzz_0_xyyyyzzz_1,  \
                             g_0_zzzzz_0_xyyyzzz_1,   \
                             g_0_zzzzz_0_xyyyzzzz_0,  \
                             g_0_zzzzz_0_xyyyzzzz_1,  \
                             g_0_zzzzz_0_xyyzzzz_1,   \
                             g_0_zzzzz_0_xyyzzzzz_0,  \
                             g_0_zzzzz_0_xyyzzzzz_1,  \
                             g_0_zzzzz_0_xyzzzzz_1,   \
                             g_0_zzzzz_0_xyzzzzzz_0,  \
                             g_0_zzzzz_0_xyzzzzzz_1,  \
                             g_0_zzzzz_0_xzzzzzz_1,   \
                             g_0_zzzzz_0_xzzzzzzz_0,  \
                             g_0_zzzzz_0_xzzzzzzz_1,  \
                             g_0_zzzzz_0_yyyyyyy_1,   \
                             g_0_zzzzz_0_yyyyyyyy_0,  \
                             g_0_zzzzz_0_yyyyyyyy_1,  \
                             g_0_zzzzz_0_yyyyyyyz_0,  \
                             g_0_zzzzz_0_yyyyyyyz_1,  \
                             g_0_zzzzz_0_yyyyyyz_1,   \
                             g_0_zzzzz_0_yyyyyyzz_0,  \
                             g_0_zzzzz_0_yyyyyyzz_1,  \
                             g_0_zzzzz_0_yyyyyzz_1,   \
                             g_0_zzzzz_0_yyyyyzzz_0,  \
                             g_0_zzzzz_0_yyyyyzzz_1,  \
                             g_0_zzzzz_0_yyyyzzz_1,   \
                             g_0_zzzzz_0_yyyyzzzz_0,  \
                             g_0_zzzzz_0_yyyyzzzz_1,  \
                             g_0_zzzzz_0_yyyzzzz_1,   \
                             g_0_zzzzz_0_yyyzzzzz_0,  \
                             g_0_zzzzz_0_yyyzzzzz_1,  \
                             g_0_zzzzz_0_yyzzzzz_1,   \
                             g_0_zzzzz_0_yyzzzzzz_0,  \
                             g_0_zzzzz_0_yyzzzzzz_1,  \
                             g_0_zzzzz_0_yzzzzzz_1,   \
                             g_0_zzzzz_0_yzzzzzzz_0,  \
                             g_0_zzzzz_0_yzzzzzzz_1,  \
                             g_0_zzzzz_0_zzzzzzz_1,   \
                             g_0_zzzzz_0_zzzzzzzz_0,  \
                             g_0_zzzzz_0_zzzzzzzz_1,  \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzz_0_xxxxxxxx_0[i] = g_0_zzzzz_0_xxxxxxxx_0[i] * pb_y + g_0_zzzzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxxxy_0[i] = g_0_zzzzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxxy_0[i] * pb_y + g_0_zzzzz_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxxxz_0[i] = g_0_zzzzz_0_xxxxxxxz_0[i] * pb_y + g_0_zzzzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxxyy_0[i] =
            2.0 * g_0_zzzzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxyy_0[i] * pb_y + g_0_zzzzz_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxxyz_0[i] = g_0_zzzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxyz_0[i] * pb_y + g_0_zzzzz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxxzz_0[i] = g_0_zzzzz_0_xxxxxxzz_0[i] * pb_y + g_0_zzzzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxyyy_0[i] =
            3.0 * g_0_zzzzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyyy_0[i] * pb_y + g_0_zzzzz_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxyyz_0[i] =
            2.0 * g_0_zzzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyyz_0[i] * pb_y + g_0_zzzzz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxyzz_0[i] = g_0_zzzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyzz_0[i] * pb_y + g_0_zzzzz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxzzz_0[i] = g_0_zzzzz_0_xxxxxzzz_0[i] * pb_y + g_0_zzzzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxyyyy_0[i] =
            4.0 * g_0_zzzzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyyy_0[i] * pb_y + g_0_zzzzz_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxyyyz_0[i] =
            3.0 * g_0_zzzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyyz_0[i] * pb_y + g_0_zzzzz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxyyzz_0[i] =
            2.0 * g_0_zzzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyzz_0[i] * pb_y + g_0_zzzzz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxyzzz_0[i] = g_0_zzzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyzzz_0[i] * pb_y + g_0_zzzzz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxzzzz_0[i] = g_0_zzzzz_0_xxxxzzzz_0[i] * pb_y + g_0_zzzzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyyyyy_0[i] =
            5.0 * g_0_zzzzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyyy_0[i] * pb_y + g_0_zzzzz_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyyyyz_0[i] =
            4.0 * g_0_zzzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyyz_0[i] * pb_y + g_0_zzzzz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyyyzz_0[i] =
            3.0 * g_0_zzzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyzz_0[i] * pb_y + g_0_zzzzz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyyzzz_0[i] =
            2.0 * g_0_zzzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyzzz_0[i] * pb_y + g_0_zzzzz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyzzzz_0[i] = g_0_zzzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyzzzz_0[i] * pb_y + g_0_zzzzz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxzzzzz_0[i] = g_0_zzzzz_0_xxxzzzzz_0[i] * pb_y + g_0_zzzzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyyyyy_0[i] =
            6.0 * g_0_zzzzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyyy_0[i] * pb_y + g_0_zzzzz_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyyyyz_0[i] =
            5.0 * g_0_zzzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyyz_0[i] * pb_y + g_0_zzzzz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyyyzz_0[i] =
            4.0 * g_0_zzzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyzz_0[i] * pb_y + g_0_zzzzz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyyzzz_0[i] =
            3.0 * g_0_zzzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyzzz_0[i] * pb_y + g_0_zzzzz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyzzzz_0[i] =
            2.0 * g_0_zzzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyzzzz_0[i] * pb_y + g_0_zzzzz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyzzzzz_0[i] = g_0_zzzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyzzzzz_0[i] * pb_y + g_0_zzzzz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxzzzzzz_0[i] = g_0_zzzzz_0_xxzzzzzz_0[i] * pb_y + g_0_zzzzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyyyyy_0[i] =
            7.0 * g_0_zzzzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyyy_0[i] * pb_y + g_0_zzzzz_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyyyyz_0[i] =
            6.0 * g_0_zzzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyyz_0[i] * pb_y + g_0_zzzzz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyyyzz_0[i] =
            5.0 * g_0_zzzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyzz_0[i] * pb_y + g_0_zzzzz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyyzzz_0[i] =
            4.0 * g_0_zzzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyzzz_0[i] * pb_y + g_0_zzzzz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyzzzz_0[i] =
            3.0 * g_0_zzzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyzzzz_0[i] * pb_y + g_0_zzzzz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyzzzzz_0[i] =
            2.0 * g_0_zzzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyzzzzz_0[i] * pb_y + g_0_zzzzz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyzzzzzz_0[i] = g_0_zzzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzzzzzz_0[i] * pb_y + g_0_zzzzz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xzzzzzzz_0[i] = g_0_zzzzz_0_xzzzzzzz_0[i] * pb_y + g_0_zzzzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyyyyy_0[i] =
            8.0 * g_0_zzzzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyyyy_0[i] * pb_y + g_0_zzzzz_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyyyyz_0[i] =
            7.0 * g_0_zzzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyyyz_0[i] * pb_y + g_0_zzzzz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyyyzz_0[i] =
            6.0 * g_0_zzzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyyzz_0[i] * pb_y + g_0_zzzzz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyyzzz_0[i] =
            5.0 * g_0_zzzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyzzz_0[i] * pb_y + g_0_zzzzz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyzzzz_0[i] =
            4.0 * g_0_zzzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyzzzz_0[i] * pb_y + g_0_zzzzz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyzzzzz_0[i] =
            3.0 * g_0_zzzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyzzzzz_0[i] * pb_y + g_0_zzzzz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyzzzzzz_0[i] =
            2.0 * g_0_zzzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyzzzzzz_0[i] * pb_y + g_0_zzzzz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yzzzzzzz_0[i] = g_0_zzzzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yzzzzzzz_0[i] * pb_y + g_0_zzzzz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_zzzzzzzz_0[i] = g_0_zzzzz_0_zzzzzzzz_0[i] * pb_y + g_0_zzzzz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1215-1260 components of targeted buffer : SISL

    auto g_0_zzzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sisl + 1215);

    auto g_0_zzzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sisl + 1216);

    auto g_0_zzzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sisl + 1217);

    auto g_0_zzzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sisl + 1218);

    auto g_0_zzzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sisl + 1219);

    auto g_0_zzzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sisl + 1220);

    auto g_0_zzzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sisl + 1221);

    auto g_0_zzzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sisl + 1222);

    auto g_0_zzzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sisl + 1223);

    auto g_0_zzzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sisl + 1224);

    auto g_0_zzzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1225);

    auto g_0_zzzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1226);

    auto g_0_zzzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1227);

    auto g_0_zzzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1228);

    auto g_0_zzzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1229);

    auto g_0_zzzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1230);

    auto g_0_zzzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1231);

    auto g_0_zzzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1232);

    auto g_0_zzzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1233);

    auto g_0_zzzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1234);

    auto g_0_zzzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1235);

    auto g_0_zzzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1236);

    auto g_0_zzzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1237);

    auto g_0_zzzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1238);

    auto g_0_zzzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1239);

    auto g_0_zzzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1240);

    auto g_0_zzzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1241);

    auto g_0_zzzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1242);

    auto g_0_zzzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1243);

    auto g_0_zzzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1244);

    auto g_0_zzzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1245);

    auto g_0_zzzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1246);

    auto g_0_zzzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1247);

    auto g_0_zzzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1248);

    auto g_0_zzzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1249);

    auto g_0_zzzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1250);

    auto g_0_zzzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sisl + 1251);

    auto g_0_zzzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sisl + 1252);

    auto g_0_zzzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sisl + 1253);

    auto g_0_zzzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sisl + 1254);

    auto g_0_zzzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1255);

    auto g_0_zzzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1256);

    auto g_0_zzzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1257);

    auto g_0_zzzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1258);

    auto g_0_zzzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sisl + 1259);

#pragma omp simd aligned(g_0_zzzz_0_xxxxxxxx_0,       \
                             g_0_zzzz_0_xxxxxxxx_1,   \
                             g_0_zzzz_0_xxxxxxxy_0,   \
                             g_0_zzzz_0_xxxxxxxy_1,   \
                             g_0_zzzz_0_xxxxxxxz_0,   \
                             g_0_zzzz_0_xxxxxxxz_1,   \
                             g_0_zzzz_0_xxxxxxyy_0,   \
                             g_0_zzzz_0_xxxxxxyy_1,   \
                             g_0_zzzz_0_xxxxxxyz_0,   \
                             g_0_zzzz_0_xxxxxxyz_1,   \
                             g_0_zzzz_0_xxxxxxzz_0,   \
                             g_0_zzzz_0_xxxxxxzz_1,   \
                             g_0_zzzz_0_xxxxxyyy_0,   \
                             g_0_zzzz_0_xxxxxyyy_1,   \
                             g_0_zzzz_0_xxxxxyyz_0,   \
                             g_0_zzzz_0_xxxxxyyz_1,   \
                             g_0_zzzz_0_xxxxxyzz_0,   \
                             g_0_zzzz_0_xxxxxyzz_1,   \
                             g_0_zzzz_0_xxxxxzzz_0,   \
                             g_0_zzzz_0_xxxxxzzz_1,   \
                             g_0_zzzz_0_xxxxyyyy_0,   \
                             g_0_zzzz_0_xxxxyyyy_1,   \
                             g_0_zzzz_0_xxxxyyyz_0,   \
                             g_0_zzzz_0_xxxxyyyz_1,   \
                             g_0_zzzz_0_xxxxyyzz_0,   \
                             g_0_zzzz_0_xxxxyyzz_1,   \
                             g_0_zzzz_0_xxxxyzzz_0,   \
                             g_0_zzzz_0_xxxxyzzz_1,   \
                             g_0_zzzz_0_xxxxzzzz_0,   \
                             g_0_zzzz_0_xxxxzzzz_1,   \
                             g_0_zzzz_0_xxxyyyyy_0,   \
                             g_0_zzzz_0_xxxyyyyy_1,   \
                             g_0_zzzz_0_xxxyyyyz_0,   \
                             g_0_zzzz_0_xxxyyyyz_1,   \
                             g_0_zzzz_0_xxxyyyzz_0,   \
                             g_0_zzzz_0_xxxyyyzz_1,   \
                             g_0_zzzz_0_xxxyyzzz_0,   \
                             g_0_zzzz_0_xxxyyzzz_1,   \
                             g_0_zzzz_0_xxxyzzzz_0,   \
                             g_0_zzzz_0_xxxyzzzz_1,   \
                             g_0_zzzz_0_xxxzzzzz_0,   \
                             g_0_zzzz_0_xxxzzzzz_1,   \
                             g_0_zzzz_0_xxyyyyyy_0,   \
                             g_0_zzzz_0_xxyyyyyy_1,   \
                             g_0_zzzz_0_xxyyyyyz_0,   \
                             g_0_zzzz_0_xxyyyyyz_1,   \
                             g_0_zzzz_0_xxyyyyzz_0,   \
                             g_0_zzzz_0_xxyyyyzz_1,   \
                             g_0_zzzz_0_xxyyyzzz_0,   \
                             g_0_zzzz_0_xxyyyzzz_1,   \
                             g_0_zzzz_0_xxyyzzzz_0,   \
                             g_0_zzzz_0_xxyyzzzz_1,   \
                             g_0_zzzz_0_xxyzzzzz_0,   \
                             g_0_zzzz_0_xxyzzzzz_1,   \
                             g_0_zzzz_0_xxzzzzzz_0,   \
                             g_0_zzzz_0_xxzzzzzz_1,   \
                             g_0_zzzz_0_xyyyyyyy_0,   \
                             g_0_zzzz_0_xyyyyyyy_1,   \
                             g_0_zzzz_0_xyyyyyyz_0,   \
                             g_0_zzzz_0_xyyyyyyz_1,   \
                             g_0_zzzz_0_xyyyyyzz_0,   \
                             g_0_zzzz_0_xyyyyyzz_1,   \
                             g_0_zzzz_0_xyyyyzzz_0,   \
                             g_0_zzzz_0_xyyyyzzz_1,   \
                             g_0_zzzz_0_xyyyzzzz_0,   \
                             g_0_zzzz_0_xyyyzzzz_1,   \
                             g_0_zzzz_0_xyyzzzzz_0,   \
                             g_0_zzzz_0_xyyzzzzz_1,   \
                             g_0_zzzz_0_xyzzzzzz_0,   \
                             g_0_zzzz_0_xyzzzzzz_1,   \
                             g_0_zzzz_0_xzzzzzzz_0,   \
                             g_0_zzzz_0_xzzzzzzz_1,   \
                             g_0_zzzz_0_yyyyyyyy_0,   \
                             g_0_zzzz_0_yyyyyyyy_1,   \
                             g_0_zzzz_0_yyyyyyyz_0,   \
                             g_0_zzzz_0_yyyyyyyz_1,   \
                             g_0_zzzz_0_yyyyyyzz_0,   \
                             g_0_zzzz_0_yyyyyyzz_1,   \
                             g_0_zzzz_0_yyyyyzzz_0,   \
                             g_0_zzzz_0_yyyyyzzz_1,   \
                             g_0_zzzz_0_yyyyzzzz_0,   \
                             g_0_zzzz_0_yyyyzzzz_1,   \
                             g_0_zzzz_0_yyyzzzzz_0,   \
                             g_0_zzzz_0_yyyzzzzz_1,   \
                             g_0_zzzz_0_yyzzzzzz_0,   \
                             g_0_zzzz_0_yyzzzzzz_1,   \
                             g_0_zzzz_0_yzzzzzzz_0,   \
                             g_0_zzzz_0_yzzzzzzz_1,   \
                             g_0_zzzz_0_zzzzzzzz_0,   \
                             g_0_zzzz_0_zzzzzzzz_1,   \
                             g_0_zzzzz_0_xxxxxxx_1,   \
                             g_0_zzzzz_0_xxxxxxxx_0,  \
                             g_0_zzzzz_0_xxxxxxxx_1,  \
                             g_0_zzzzz_0_xxxxxxxy_0,  \
                             g_0_zzzzz_0_xxxxxxxy_1,  \
                             g_0_zzzzz_0_xxxxxxxz_0,  \
                             g_0_zzzzz_0_xxxxxxxz_1,  \
                             g_0_zzzzz_0_xxxxxxy_1,   \
                             g_0_zzzzz_0_xxxxxxyy_0,  \
                             g_0_zzzzz_0_xxxxxxyy_1,  \
                             g_0_zzzzz_0_xxxxxxyz_0,  \
                             g_0_zzzzz_0_xxxxxxyz_1,  \
                             g_0_zzzzz_0_xxxxxxz_1,   \
                             g_0_zzzzz_0_xxxxxxzz_0,  \
                             g_0_zzzzz_0_xxxxxxzz_1,  \
                             g_0_zzzzz_0_xxxxxyy_1,   \
                             g_0_zzzzz_0_xxxxxyyy_0,  \
                             g_0_zzzzz_0_xxxxxyyy_1,  \
                             g_0_zzzzz_0_xxxxxyyz_0,  \
                             g_0_zzzzz_0_xxxxxyyz_1,  \
                             g_0_zzzzz_0_xxxxxyz_1,   \
                             g_0_zzzzz_0_xxxxxyzz_0,  \
                             g_0_zzzzz_0_xxxxxyzz_1,  \
                             g_0_zzzzz_0_xxxxxzz_1,   \
                             g_0_zzzzz_0_xxxxxzzz_0,  \
                             g_0_zzzzz_0_xxxxxzzz_1,  \
                             g_0_zzzzz_0_xxxxyyy_1,   \
                             g_0_zzzzz_0_xxxxyyyy_0,  \
                             g_0_zzzzz_0_xxxxyyyy_1,  \
                             g_0_zzzzz_0_xxxxyyyz_0,  \
                             g_0_zzzzz_0_xxxxyyyz_1,  \
                             g_0_zzzzz_0_xxxxyyz_1,   \
                             g_0_zzzzz_0_xxxxyyzz_0,  \
                             g_0_zzzzz_0_xxxxyyzz_1,  \
                             g_0_zzzzz_0_xxxxyzz_1,   \
                             g_0_zzzzz_0_xxxxyzzz_0,  \
                             g_0_zzzzz_0_xxxxyzzz_1,  \
                             g_0_zzzzz_0_xxxxzzz_1,   \
                             g_0_zzzzz_0_xxxxzzzz_0,  \
                             g_0_zzzzz_0_xxxxzzzz_1,  \
                             g_0_zzzzz_0_xxxyyyy_1,   \
                             g_0_zzzzz_0_xxxyyyyy_0,  \
                             g_0_zzzzz_0_xxxyyyyy_1,  \
                             g_0_zzzzz_0_xxxyyyyz_0,  \
                             g_0_zzzzz_0_xxxyyyyz_1,  \
                             g_0_zzzzz_0_xxxyyyz_1,   \
                             g_0_zzzzz_0_xxxyyyzz_0,  \
                             g_0_zzzzz_0_xxxyyyzz_1,  \
                             g_0_zzzzz_0_xxxyyzz_1,   \
                             g_0_zzzzz_0_xxxyyzzz_0,  \
                             g_0_zzzzz_0_xxxyyzzz_1,  \
                             g_0_zzzzz_0_xxxyzzz_1,   \
                             g_0_zzzzz_0_xxxyzzzz_0,  \
                             g_0_zzzzz_0_xxxyzzzz_1,  \
                             g_0_zzzzz_0_xxxzzzz_1,   \
                             g_0_zzzzz_0_xxxzzzzz_0,  \
                             g_0_zzzzz_0_xxxzzzzz_1,  \
                             g_0_zzzzz_0_xxyyyyy_1,   \
                             g_0_zzzzz_0_xxyyyyyy_0,  \
                             g_0_zzzzz_0_xxyyyyyy_1,  \
                             g_0_zzzzz_0_xxyyyyyz_0,  \
                             g_0_zzzzz_0_xxyyyyyz_1,  \
                             g_0_zzzzz_0_xxyyyyz_1,   \
                             g_0_zzzzz_0_xxyyyyzz_0,  \
                             g_0_zzzzz_0_xxyyyyzz_1,  \
                             g_0_zzzzz_0_xxyyyzz_1,   \
                             g_0_zzzzz_0_xxyyyzzz_0,  \
                             g_0_zzzzz_0_xxyyyzzz_1,  \
                             g_0_zzzzz_0_xxyyzzz_1,   \
                             g_0_zzzzz_0_xxyyzzzz_0,  \
                             g_0_zzzzz_0_xxyyzzzz_1,  \
                             g_0_zzzzz_0_xxyzzzz_1,   \
                             g_0_zzzzz_0_xxyzzzzz_0,  \
                             g_0_zzzzz_0_xxyzzzzz_1,  \
                             g_0_zzzzz_0_xxzzzzz_1,   \
                             g_0_zzzzz_0_xxzzzzzz_0,  \
                             g_0_zzzzz_0_xxzzzzzz_1,  \
                             g_0_zzzzz_0_xyyyyyy_1,   \
                             g_0_zzzzz_0_xyyyyyyy_0,  \
                             g_0_zzzzz_0_xyyyyyyy_1,  \
                             g_0_zzzzz_0_xyyyyyyz_0,  \
                             g_0_zzzzz_0_xyyyyyyz_1,  \
                             g_0_zzzzz_0_xyyyyyz_1,   \
                             g_0_zzzzz_0_xyyyyyzz_0,  \
                             g_0_zzzzz_0_xyyyyyzz_1,  \
                             g_0_zzzzz_0_xyyyyzz_1,   \
                             g_0_zzzzz_0_xyyyyzzz_0,  \
                             g_0_zzzzz_0_xyyyyzzz_1,  \
                             g_0_zzzzz_0_xyyyzzz_1,   \
                             g_0_zzzzz_0_xyyyzzzz_0,  \
                             g_0_zzzzz_0_xyyyzzzz_1,  \
                             g_0_zzzzz_0_xyyzzzz_1,   \
                             g_0_zzzzz_0_xyyzzzzz_0,  \
                             g_0_zzzzz_0_xyyzzzzz_1,  \
                             g_0_zzzzz_0_xyzzzzz_1,   \
                             g_0_zzzzz_0_xyzzzzzz_0,  \
                             g_0_zzzzz_0_xyzzzzzz_1,  \
                             g_0_zzzzz_0_xzzzzzz_1,   \
                             g_0_zzzzz_0_xzzzzzzz_0,  \
                             g_0_zzzzz_0_xzzzzzzz_1,  \
                             g_0_zzzzz_0_yyyyyyy_1,   \
                             g_0_zzzzz_0_yyyyyyyy_0,  \
                             g_0_zzzzz_0_yyyyyyyy_1,  \
                             g_0_zzzzz_0_yyyyyyyz_0,  \
                             g_0_zzzzz_0_yyyyyyyz_1,  \
                             g_0_zzzzz_0_yyyyyyz_1,   \
                             g_0_zzzzz_0_yyyyyyzz_0,  \
                             g_0_zzzzz_0_yyyyyyzz_1,  \
                             g_0_zzzzz_0_yyyyyzz_1,   \
                             g_0_zzzzz_0_yyyyyzzz_0,  \
                             g_0_zzzzz_0_yyyyyzzz_1,  \
                             g_0_zzzzz_0_yyyyzzz_1,   \
                             g_0_zzzzz_0_yyyyzzzz_0,  \
                             g_0_zzzzz_0_yyyyzzzz_1,  \
                             g_0_zzzzz_0_yyyzzzz_1,   \
                             g_0_zzzzz_0_yyyzzzzz_0,  \
                             g_0_zzzzz_0_yyyzzzzz_1,  \
                             g_0_zzzzz_0_yyzzzzz_1,   \
                             g_0_zzzzz_0_yyzzzzzz_0,  \
                             g_0_zzzzz_0_yyzzzzzz_1,  \
                             g_0_zzzzz_0_yzzzzzz_1,   \
                             g_0_zzzzz_0_yzzzzzzz_0,  \
                             g_0_zzzzz_0_yzzzzzzz_1,  \
                             g_0_zzzzz_0_zzzzzzz_1,   \
                             g_0_zzzzz_0_zzzzzzzz_0,  \
                             g_0_zzzzz_0_zzzzzzzz_1,  \
                             g_0_zzzzzz_0_xxxxxxxx_0, \
                             g_0_zzzzzz_0_xxxxxxxy_0, \
                             g_0_zzzzzz_0_xxxxxxxz_0, \
                             g_0_zzzzzz_0_xxxxxxyy_0, \
                             g_0_zzzzzz_0_xxxxxxyz_0, \
                             g_0_zzzzzz_0_xxxxxxzz_0, \
                             g_0_zzzzzz_0_xxxxxyyy_0, \
                             g_0_zzzzzz_0_xxxxxyyz_0, \
                             g_0_zzzzzz_0_xxxxxyzz_0, \
                             g_0_zzzzzz_0_xxxxxzzz_0, \
                             g_0_zzzzzz_0_xxxxyyyy_0, \
                             g_0_zzzzzz_0_xxxxyyyz_0, \
                             g_0_zzzzzz_0_xxxxyyzz_0, \
                             g_0_zzzzzz_0_xxxxyzzz_0, \
                             g_0_zzzzzz_0_xxxxzzzz_0, \
                             g_0_zzzzzz_0_xxxyyyyy_0, \
                             g_0_zzzzzz_0_xxxyyyyz_0, \
                             g_0_zzzzzz_0_xxxyyyzz_0, \
                             g_0_zzzzzz_0_xxxyyzzz_0, \
                             g_0_zzzzzz_0_xxxyzzzz_0, \
                             g_0_zzzzzz_0_xxxzzzzz_0, \
                             g_0_zzzzzz_0_xxyyyyyy_0, \
                             g_0_zzzzzz_0_xxyyyyyz_0, \
                             g_0_zzzzzz_0_xxyyyyzz_0, \
                             g_0_zzzzzz_0_xxyyyzzz_0, \
                             g_0_zzzzzz_0_xxyyzzzz_0, \
                             g_0_zzzzzz_0_xxyzzzzz_0, \
                             g_0_zzzzzz_0_xxzzzzzz_0, \
                             g_0_zzzzzz_0_xyyyyyyy_0, \
                             g_0_zzzzzz_0_xyyyyyyz_0, \
                             g_0_zzzzzz_0_xyyyyyzz_0, \
                             g_0_zzzzzz_0_xyyyyzzz_0, \
                             g_0_zzzzzz_0_xyyyzzzz_0, \
                             g_0_zzzzzz_0_xyyzzzzz_0, \
                             g_0_zzzzzz_0_xyzzzzzz_0, \
                             g_0_zzzzzz_0_xzzzzzzz_0, \
                             g_0_zzzzzz_0_yyyyyyyy_0, \
                             g_0_zzzzzz_0_yyyyyyyz_0, \
                             g_0_zzzzzz_0_yyyyyyzz_0, \
                             g_0_zzzzzz_0_yyyyyzzz_0, \
                             g_0_zzzzzz_0_yyyyzzzz_0, \
                             g_0_zzzzzz_0_yyyzzzzz_0, \
                             g_0_zzzzzz_0_yyzzzzzz_0, \
                             g_0_zzzzzz_0_yzzzzzzz_0, \
                             g_0_zzzzzz_0_zzzzzzzz_0, \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzz_0_xxxxxxxx_0[i] = 5.0 * g_0_zzzz_0_xxxxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxxxx_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxxxxxxx_0[i] * pb_z + g_0_zzzzz_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxxxy_0[i] = 5.0 * g_0_zzzz_0_xxxxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxxxy_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxxxxxxy_0[i] * pb_z + g_0_zzzzz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxxxz_0[i] = 5.0 * g_0_zzzz_0_xxxxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxxxz_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxxz_0[i] * pb_z + g_0_zzzzz_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxxyy_0[i] = 5.0 * g_0_zzzz_0_xxxxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxxyy_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxxxxxyy_0[i] * pb_z + g_0_zzzzz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxxyz_0[i] = 5.0 * g_0_zzzz_0_xxxxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxxyz_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxyz_0[i] * pb_z + g_0_zzzzz_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxxzz_0[i] = 5.0 * g_0_zzzz_0_xxxxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxxzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_zzzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxyyy_0[i] = 5.0 * g_0_zzzz_0_xxxxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxyyy_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxxxxyyy_0[i] * pb_z + g_0_zzzzz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxyyz_0[i] = 5.0 * g_0_zzzz_0_xxxxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxyyz_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyyz_0[i] * pb_z + g_0_zzzzz_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxyzz_0[i] = 5.0 * g_0_zzzz_0_xxxxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_zzzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxzzz_0[i] = 5.0 * g_0_zzzz_0_xxxxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_zzzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxyyyy_0[i] = 5.0 * g_0_zzzz_0_xxxxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxyyyy_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxxxyyyy_0[i] * pb_z + g_0_zzzzz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxyyyz_0[i] = 5.0 * g_0_zzzz_0_xxxxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxyyyz_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyyz_0[i] * pb_z + g_0_zzzzz_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxyyzz_0[i] = 5.0 * g_0_zzzz_0_xxxxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_zzzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxyzzz_0[i] = 5.0 * g_0_zzzz_0_xxxxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_zzzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxzzzz_0[i] = 5.0 * g_0_zzzz_0_xxxxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_zzzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyyyyy_0[i] = 5.0 * g_0_zzzz_0_xxxyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyyyyy_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxxyyyyy_0[i] * pb_z + g_0_zzzzz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyyyyz_0[i] = 5.0 * g_0_zzzz_0_xxxyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyyyyz_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyyz_0[i] * pb_z + g_0_zzzzz_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyyyzz_0[i] = 5.0 * g_0_zzzz_0_xxxyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_zzzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyyzzz_0[i] = 5.0 * g_0_zzzz_0_xxxyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_zzzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyzzzz_0[i] = 5.0 * g_0_zzzz_0_xxxyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_zzzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxzzzzz_0[i] = 5.0 * g_0_zzzz_0_xxxzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxzzzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_zzzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxzzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyyyyy_0[i] = 5.0 * g_0_zzzz_0_xxyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyyyyy_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxyyyyyy_0[i] * pb_z + g_0_zzzzz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyyyyz_0[i] = 5.0 * g_0_zzzz_0_xxyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyyyyz_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyyz_0[i] * pb_z + g_0_zzzzz_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyyyzz_0[i] = 5.0 * g_0_zzzz_0_xxyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_zzzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyyzzz_0[i] = 5.0 * g_0_zzzz_0_xxyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_zzzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyzzzz_0[i] = 5.0 * g_0_zzzz_0_xxyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_zzzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyzzzzz_0[i] = 5.0 * g_0_zzzz_0_xxyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyzzzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_zzzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyzzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxzzzzzz_0[i] = 5.0 * g_0_zzzz_0_xxzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxzzzzzz_1[i] * fti_ab_0 +
                                     6.0 * g_0_zzzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxzzzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyyyyy_0[i] = 5.0 * g_0_zzzz_0_xyyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xyyyyyyy_0[i] * pb_z + g_0_zzzzz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyyyyz_0[i] = 5.0 * g_0_zzzz_0_xyyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyyyyz_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyyz_0[i] * pb_z + g_0_zzzzz_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyyyzz_0[i] = 5.0 * g_0_zzzz_0_xyyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_zzzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyyzzz_0[i] = 5.0 * g_0_zzzz_0_xyyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_zzzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyzzzz_0[i] = 5.0 * g_0_zzzz_0_xyyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_zzzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyzzzzz_0[i] = 5.0 * g_0_zzzz_0_xyyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyzzzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_zzzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyzzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyzzzzzz_0[i] = 5.0 * g_0_zzzz_0_xyzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyzzzzzz_1[i] * fti_ab_0 +
                                     6.0 * g_0_zzzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzzzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xzzzzzzz_0[i] = 5.0 * g_0_zzzz_0_xzzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xzzzzzzz_1[i] * fti_ab_0 +
                                     7.0 * g_0_zzzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xzzzzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyyyyy_0[i] = 5.0 * g_0_zzzz_0_yyyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyyyyy_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_yyyyyyyy_0[i] * pb_z + g_0_zzzzz_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyyyyz_0[i] = 5.0 * g_0_zzzz_0_yyyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyyyyz_1[i] * fti_ab_0 +
                                     g_0_zzzzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyyyz_0[i] * pb_z + g_0_zzzzz_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyyyzz_0[i] = 5.0 * g_0_zzzz_0_yyyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_zzzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyyzz_0[i] * pb_z +
                                     g_0_zzzzz_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyyzzz_0[i] = 5.0 * g_0_zzzz_0_yyyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_zzzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyzzzz_0[i] = 5.0 * g_0_zzzz_0_yyyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_zzzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyzzzzz_0[i] = 5.0 * g_0_zzzz_0_yyyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyzzzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_zzzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyzzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyzzzzzz_0[i] = 5.0 * g_0_zzzz_0_yyzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyzzzzzz_1[i] * fti_ab_0 +
                                     6.0 * g_0_zzzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyzzzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yzzzzzzz_0[i] = 5.0 * g_0_zzzz_0_yzzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yzzzzzzz_1[i] * fti_ab_0 +
                                     7.0 * g_0_zzzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yzzzzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_zzzzzzzz_0[i] = 5.0 * g_0_zzzz_0_zzzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_zzzzzzzz_1[i] * fti_ab_0 +
                                     8.0 * g_0_zzzzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_zzzzzzzz_0[i] * pb_z +
                                     g_0_zzzzz_0_zzzzzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
