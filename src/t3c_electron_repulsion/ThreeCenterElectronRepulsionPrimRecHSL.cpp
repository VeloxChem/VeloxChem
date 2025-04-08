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

#include "ThreeCenterElectronRepulsionPrimRecHSL.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_hsl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsl,
                                 size_t idx_eri_0_fsl,
                                 size_t idx_eri_1_fsl,
                                 size_t idx_eri_1_gsk,
                                 size_t idx_eri_1_gsl,
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

    /// Set up components of auxilary buffer : FSL

    auto g_xxx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl);

    auto g_xxx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 1);

    auto g_xxx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 2);

    auto g_xxx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 3);

    auto g_xxx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 4);

    auto g_xxx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 5);

    auto g_xxx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 6);

    auto g_xxx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 7);

    auto g_xxx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 8);

    auto g_xxx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 9);

    auto g_xxx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 10);

    auto g_xxx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 11);

    auto g_xxx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 12);

    auto g_xxx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 13);

    auto g_xxx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 14);

    auto g_xxx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 15);

    auto g_xxx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 16);

    auto g_xxx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 17);

    auto g_xxx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 18);

    auto g_xxx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 19);

    auto g_xxx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 20);

    auto g_xxx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 21);

    auto g_xxx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 22);

    auto g_xxx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 23);

    auto g_xxx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 24);

    auto g_xxx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 25);

    auto g_xxx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 26);

    auto g_xxx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 27);

    auto g_xxx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 28);

    auto g_xxx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 29);

    auto g_xxx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 30);

    auto g_xxx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 31);

    auto g_xxx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 32);

    auto g_xxx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 33);

    auto g_xxx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 34);

    auto g_xxx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 35);

    auto g_xxx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 36);

    auto g_xxx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 37);

    auto g_xxx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 38);

    auto g_xxx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 39);

    auto g_xxx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 40);

    auto g_xxx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 41);

    auto g_xxx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 42);

    auto g_xxx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 43);

    auto g_xxx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 44);

    auto g_xxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 45);

    auto g_xxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 47);

    auto g_xxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 50);

    auto g_xxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 54);

    auto g_xxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 59);

    auto g_xxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 65);

    auto g_xxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 72);

    auto g_xxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 80);

    auto g_xxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 90);

    auto g_xxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 91);

    auto g_xxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 93);

    auto g_xxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 96);

    auto g_xxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 100);

    auto g_xxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 105);

    auto g_xxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 111);

    auto g_xxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 118);

    auto g_xyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 136);

    auto g_xyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 138);

    auto g_xyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 139);

    auto g_xyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 141);

    auto g_xyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 142);

    auto g_xyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 143);

    auto g_xyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 145);

    auto g_xyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 146);

    auto g_xyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 147);

    auto g_xyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 148);

    auto g_xyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 150);

    auto g_xyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 151);

    auto g_xyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 152);

    auto g_xyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 153);

    auto g_xyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 154);

    auto g_xyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 156);

    auto g_xyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 157);

    auto g_xyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 158);

    auto g_xyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 159);

    auto g_xyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 160);

    auto g_xyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 161);

    auto g_xyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 163);

    auto g_xyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 164);

    auto g_xyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 165);

    auto g_xyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 166);

    auto g_xyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 167);

    auto g_xyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 168);

    auto g_xyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 169);

    auto g_xyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 171);

    auto g_xyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 172);

    auto g_xyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 173);

    auto g_xyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 174);

    auto g_xyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 175);

    auto g_xyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 176);

    auto g_xyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 177);

    auto g_xyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 178);

    auto g_xyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 179);

    auto g_xzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 227);

    auto g_xzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 229);

    auto g_xzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 230);

    auto g_xzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 232);

    auto g_xzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 233);

    auto g_xzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 234);

    auto g_xzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 236);

    auto g_xzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 237);

    auto g_xzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 238);

    auto g_xzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 239);

    auto g_xzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 241);

    auto g_xzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 242);

    auto g_xzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 243);

    auto g_xzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 244);

    auto g_xzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 245);

    auto g_xzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 247);

    auto g_xzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 248);

    auto g_xzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 249);

    auto g_xzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 250);

    auto g_xzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 251);

    auto g_xzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 252);

    auto g_xzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 254);

    auto g_xzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 255);

    auto g_xzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 256);

    auto g_xzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 257);

    auto g_xzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 258);

    auto g_xzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 259);

    auto g_xzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 260);

    auto g_xzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 261);

    auto g_xzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 262);

    auto g_xzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 263);

    auto g_xzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 264);

    auto g_xzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 265);

    auto g_xzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 266);

    auto g_xzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 267);

    auto g_xzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 268);

    auto g_xzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 269);

    auto g_yyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 270);

    auto g_yyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 271);

    auto g_yyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 272);

    auto g_yyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 273);

    auto g_yyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 274);

    auto g_yyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 275);

    auto g_yyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 276);

    auto g_yyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 277);

    auto g_yyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 278);

    auto g_yyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 279);

    auto g_yyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 280);

    auto g_yyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 281);

    auto g_yyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 282);

    auto g_yyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 283);

    auto g_yyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 284);

    auto g_yyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 285);

    auto g_yyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 286);

    auto g_yyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 287);

    auto g_yyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 288);

    auto g_yyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 289);

    auto g_yyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 290);

    auto g_yyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 291);

    auto g_yyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 292);

    auto g_yyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 293);

    auto g_yyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 294);

    auto g_yyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 295);

    auto g_yyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 296);

    auto g_yyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 297);

    auto g_yyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 298);

    auto g_yyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 299);

    auto g_yyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 300);

    auto g_yyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 301);

    auto g_yyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 302);

    auto g_yyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 303);

    auto g_yyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 304);

    auto g_yyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 305);

    auto g_yyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 306);

    auto g_yyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 307);

    auto g_yyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 308);

    auto g_yyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 309);

    auto g_yyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 310);

    auto g_yyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 311);

    auto g_yyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 312);

    auto g_yyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 313);

    auto g_yyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 314);

    auto g_yyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 316);

    auto g_yyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 318);

    auto g_yyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 321);

    auto g_yyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 325);

    auto g_yyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 330);

    auto g_yyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 336);

    auto g_yyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 343);

    auto g_yyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 351);

    auto g_yzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 360);

    auto g_yzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 362);

    auto g_yzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 364);

    auto g_yzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 365);

    auto g_yzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 367);

    auto g_yzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 368);

    auto g_yzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 369);

    auto g_yzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 371);

    auto g_yzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 372);

    auto g_yzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 373);

    auto g_yzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 374);

    auto g_yzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 376);

    auto g_yzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 377);

    auto g_yzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 378);

    auto g_yzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 379);

    auto g_yzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 380);

    auto g_yzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 382);

    auto g_yzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 383);

    auto g_yzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 384);

    auto g_yzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 385);

    auto g_yzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 386);

    auto g_yzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 387);

    auto g_yzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 389);

    auto g_yzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 390);

    auto g_yzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 391);

    auto g_yzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 392);

    auto g_yzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 393);

    auto g_yzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 394);

    auto g_yzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 395);

    auto g_yzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 397);

    auto g_yzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 398);

    auto g_yzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 399);

    auto g_yzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 400);

    auto g_yzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 401);

    auto g_yzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 402);

    auto g_yzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 403);

    auto g_yzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 404);

    auto g_zzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_fsl + 405);

    auto g_zzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_fsl + 406);

    auto g_zzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_fsl + 407);

    auto g_zzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_fsl + 408);

    auto g_zzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_fsl + 409);

    auto g_zzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_fsl + 410);

    auto g_zzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_fsl + 411);

    auto g_zzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_fsl + 412);

    auto g_zzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_fsl + 413);

    auto g_zzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_fsl + 414);

    auto g_zzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_fsl + 415);

    auto g_zzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_fsl + 416);

    auto g_zzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_fsl + 417);

    auto g_zzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_fsl + 418);

    auto g_zzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_fsl + 419);

    auto g_zzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 420);

    auto g_zzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 421);

    auto g_zzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 422);

    auto g_zzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 423);

    auto g_zzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 424);

    auto g_zzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 425);

    auto g_zzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 426);

    auto g_zzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 427);

    auto g_zzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 428);

    auto g_zzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 429);

    auto g_zzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 430);

    auto g_zzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 431);

    auto g_zzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 432);

    auto g_zzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 433);

    auto g_zzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 434);

    auto g_zzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 435);

    auto g_zzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 436);

    auto g_zzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 437);

    auto g_zzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 438);

    auto g_zzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 439);

    auto g_zzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 440);

    auto g_zzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_fsl + 441);

    auto g_zzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_fsl + 442);

    auto g_zzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_fsl + 443);

    auto g_zzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_fsl + 444);

    auto g_zzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_fsl + 445);

    auto g_zzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 446);

    auto g_zzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 447);

    auto g_zzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 448);

    auto g_zzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_fsl + 449);

    /// Set up components of auxilary buffer : FSL

    auto g_xxx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl);

    auto g_xxx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 1);

    auto g_xxx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 2);

    auto g_xxx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 3);

    auto g_xxx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 4);

    auto g_xxx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 5);

    auto g_xxx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 6);

    auto g_xxx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 7);

    auto g_xxx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 8);

    auto g_xxx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 9);

    auto g_xxx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 10);

    auto g_xxx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 11);

    auto g_xxx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 12);

    auto g_xxx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 13);

    auto g_xxx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 14);

    auto g_xxx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 15);

    auto g_xxx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 16);

    auto g_xxx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 17);

    auto g_xxx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 18);

    auto g_xxx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 19);

    auto g_xxx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 20);

    auto g_xxx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 21);

    auto g_xxx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 22);

    auto g_xxx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 23);

    auto g_xxx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 24);

    auto g_xxx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 25);

    auto g_xxx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 26);

    auto g_xxx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 27);

    auto g_xxx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 28);

    auto g_xxx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 29);

    auto g_xxx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 30);

    auto g_xxx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 31);

    auto g_xxx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 32);

    auto g_xxx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 33);

    auto g_xxx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 34);

    auto g_xxx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 35);

    auto g_xxx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 36);

    auto g_xxx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 37);

    auto g_xxx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 38);

    auto g_xxx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 39);

    auto g_xxx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 40);

    auto g_xxx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 41);

    auto g_xxx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 42);

    auto g_xxx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 43);

    auto g_xxx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 44);

    auto g_xxy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl + 45);

    auto g_xxy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 47);

    auto g_xxy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 50);

    auto g_xxy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 54);

    auto g_xxy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 59);

    auto g_xxy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 65);

    auto g_xxy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 72);

    auto g_xxy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 80);

    auto g_xxz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl + 90);

    auto g_xxz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 91);

    auto g_xxz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 93);

    auto g_xxz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 96);

    auto g_xxz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 100);

    auto g_xxz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 105);

    auto g_xxz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 111);

    auto g_xxz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 118);

    auto g_xyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 136);

    auto g_xyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 138);

    auto g_xyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 139);

    auto g_xyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 141);

    auto g_xyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 142);

    auto g_xyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 143);

    auto g_xyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 145);

    auto g_xyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 146);

    auto g_xyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 147);

    auto g_xyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 148);

    auto g_xyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 150);

    auto g_xyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 151);

    auto g_xyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 152);

    auto g_xyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 153);

    auto g_xyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 154);

    auto g_xyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 156);

    auto g_xyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 157);

    auto g_xyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 158);

    auto g_xyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 159);

    auto g_xyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 160);

    auto g_xyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 161);

    auto g_xyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 163);

    auto g_xyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 164);

    auto g_xyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 165);

    auto g_xyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 166);

    auto g_xyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 167);

    auto g_xyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 168);

    auto g_xyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 169);

    auto g_xyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 171);

    auto g_xyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 172);

    auto g_xyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 173);

    auto g_xyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 174);

    auto g_xyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 175);

    auto g_xyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 176);

    auto g_xyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 177);

    auto g_xyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 178);

    auto g_xyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 179);

    auto g_xzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 227);

    auto g_xzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 229);

    auto g_xzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 230);

    auto g_xzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 232);

    auto g_xzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 233);

    auto g_xzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 234);

    auto g_xzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 236);

    auto g_xzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 237);

    auto g_xzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 238);

    auto g_xzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 239);

    auto g_xzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 241);

    auto g_xzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 242);

    auto g_xzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 243);

    auto g_xzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 244);

    auto g_xzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 245);

    auto g_xzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 247);

    auto g_xzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 248);

    auto g_xzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 249);

    auto g_xzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 250);

    auto g_xzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 251);

    auto g_xzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 252);

    auto g_xzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 254);

    auto g_xzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 255);

    auto g_xzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 256);

    auto g_xzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 257);

    auto g_xzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 258);

    auto g_xzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 259);

    auto g_xzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 260);

    auto g_xzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 261);

    auto g_xzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 262);

    auto g_xzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 263);

    auto g_xzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 264);

    auto g_xzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 265);

    auto g_xzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 266);

    auto g_xzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 267);

    auto g_xzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 268);

    auto g_xzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 269);

    auto g_yyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl + 270);

    auto g_yyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 271);

    auto g_yyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 272);

    auto g_yyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 273);

    auto g_yyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 274);

    auto g_yyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 275);

    auto g_yyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 276);

    auto g_yyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 277);

    auto g_yyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 278);

    auto g_yyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 279);

    auto g_yyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 280);

    auto g_yyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 281);

    auto g_yyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 282);

    auto g_yyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 283);

    auto g_yyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 284);

    auto g_yyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 285);

    auto g_yyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 286);

    auto g_yyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 287);

    auto g_yyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 288);

    auto g_yyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 289);

    auto g_yyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 290);

    auto g_yyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 291);

    auto g_yyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 292);

    auto g_yyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 293);

    auto g_yyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 294);

    auto g_yyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 295);

    auto g_yyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 296);

    auto g_yyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 297);

    auto g_yyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 298);

    auto g_yyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 299);

    auto g_yyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 300);

    auto g_yyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 301);

    auto g_yyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 302);

    auto g_yyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 303);

    auto g_yyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 304);

    auto g_yyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 305);

    auto g_yyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 306);

    auto g_yyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 307);

    auto g_yyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 308);

    auto g_yyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 309);

    auto g_yyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 310);

    auto g_yyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 311);

    auto g_yyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 312);

    auto g_yyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 313);

    auto g_yyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 314);

    auto g_yyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 316);

    auto g_yyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 318);

    auto g_yyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 321);

    auto g_yyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 325);

    auto g_yyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 330);

    auto g_yyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 336);

    auto g_yyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 343);

    auto g_yyz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 351);

    auto g_yzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl + 360);

    auto g_yzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 362);

    auto g_yzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 364);

    auto g_yzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 365);

    auto g_yzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 367);

    auto g_yzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 368);

    auto g_yzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 369);

    auto g_yzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 371);

    auto g_yzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 372);

    auto g_yzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 373);

    auto g_yzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 374);

    auto g_yzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 376);

    auto g_yzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 377);

    auto g_yzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 378);

    auto g_yzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 379);

    auto g_yzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 380);

    auto g_yzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 382);

    auto g_yzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 383);

    auto g_yzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 384);

    auto g_yzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 385);

    auto g_yzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 386);

    auto g_yzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 387);

    auto g_yzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 389);

    auto g_yzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 390);

    auto g_yzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 391);

    auto g_yzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 392);

    auto g_yzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 393);

    auto g_yzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 394);

    auto g_yzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 395);

    auto g_yzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 397);

    auto g_yzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 398);

    auto g_yzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 399);

    auto g_yzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 400);

    auto g_yzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 401);

    auto g_yzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 402);

    auto g_yzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 403);

    auto g_yzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 404);

    auto g_zzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl + 405);

    auto g_zzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 406);

    auto g_zzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 407);

    auto g_zzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 408);

    auto g_zzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 409);

    auto g_zzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 410);

    auto g_zzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 411);

    auto g_zzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 412);

    auto g_zzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 413);

    auto g_zzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 414);

    auto g_zzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 415);

    auto g_zzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 416);

    auto g_zzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 417);

    auto g_zzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 418);

    auto g_zzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 419);

    auto g_zzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 420);

    auto g_zzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 421);

    auto g_zzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 422);

    auto g_zzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 423);

    auto g_zzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 424);

    auto g_zzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 425);

    auto g_zzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 426);

    auto g_zzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 427);

    auto g_zzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 428);

    auto g_zzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 429);

    auto g_zzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 430);

    auto g_zzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 431);

    auto g_zzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 432);

    auto g_zzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 433);

    auto g_zzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 434);

    auto g_zzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 435);

    auto g_zzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 436);

    auto g_zzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 437);

    auto g_zzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 438);

    auto g_zzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 439);

    auto g_zzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 440);

    auto g_zzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 441);

    auto g_zzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 442);

    auto g_zzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 443);

    auto g_zzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 444);

    auto g_zzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 445);

    auto g_zzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 446);

    auto g_zzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 447);

    auto g_zzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 448);

    auto g_zzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 449);

    /// Set up components of auxilary buffer : GSK

    auto g_xxxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk);

    auto g_xxxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 1);

    auto g_xxxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 2);

    auto g_xxxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 3);

    auto g_xxxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 4);

    auto g_xxxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 5);

    auto g_xxxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 6);

    auto g_xxxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 7);

    auto g_xxxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 8);

    auto g_xxxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 9);

    auto g_xxxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 10);

    auto g_xxxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 11);

    auto g_xxxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 12);

    auto g_xxxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 13);

    auto g_xxxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 14);

    auto g_xxxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 15);

    auto g_xxxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 16);

    auto g_xxxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 17);

    auto g_xxxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 18);

    auto g_xxxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 19);

    auto g_xxxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 20);

    auto g_xxxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 21);

    auto g_xxxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 22);

    auto g_xxxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 23);

    auto g_xxxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 24);

    auto g_xxxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 25);

    auto g_xxxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 26);

    auto g_xxxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 27);

    auto g_xxxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 28);

    auto g_xxxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 29);

    auto g_xxxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 30);

    auto g_xxxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 31);

    auto g_xxxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 32);

    auto g_xxxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 33);

    auto g_xxxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 34);

    auto g_xxxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 35);

    auto g_xxxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 74);

    auto g_xxxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 76);

    auto g_xxxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 77);

    auto g_xxxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 79);

    auto g_xxxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 80);

    auto g_xxxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 81);

    auto g_xxxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 83);

    auto g_xxxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 84);

    auto g_xxxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 85);

    auto g_xxxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 86);

    auto g_xxxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 88);

    auto g_xxxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 89);

    auto g_xxxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 90);

    auto g_xxxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 91);

    auto g_xxxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 92);

    auto g_xxxz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 94);

    auto g_xxxz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 95);

    auto g_xxxz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 96);

    auto g_xxxz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 97);

    auto g_xxxz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 98);

    auto g_xxxz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 99);

    auto g_xxxz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 101);

    auto g_xxxz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 102);

    auto g_xxxz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 103);

    auto g_xxxz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 104);

    auto g_xxxz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 105);

    auto g_xxxz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 106);

    auto g_xxxz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 107);

    auto g_xxyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 108);

    auto g_xxyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 109);

    auto g_xxyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 110);

    auto g_xxyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 111);

    auto g_xxyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 112);

    auto g_xxyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 113);

    auto g_xxyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 114);

    auto g_xxyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 115);

    auto g_xxyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 116);

    auto g_xxyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 117);

    auto g_xxyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 118);

    auto g_xxyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 119);

    auto g_xxyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 120);

    auto g_xxyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 121);

    auto g_xxyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 122);

    auto g_xxyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 123);

    auto g_xxyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 124);

    auto g_xxyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 125);

    auto g_xxyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 126);

    auto g_xxyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 127);

    auto g_xxyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 128);

    auto g_xxyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 129);

    auto g_xxyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 130);

    auto g_xxyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 131);

    auto g_xxyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 132);

    auto g_xxyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 133);

    auto g_xxyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 134);

    auto g_xxyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 135);

    auto g_xxyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 136);

    auto g_xxyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 137);

    auto g_xxyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 138);

    auto g_xxyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 139);

    auto g_xxyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 140);

    auto g_xxyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 141);

    auto g_xxyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 142);

    auto g_xxyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 143);

    auto g_xxzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 180);

    auto g_xxzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 181);

    auto g_xxzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 182);

    auto g_xxzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 183);

    auto g_xxzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 184);

    auto g_xxzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 185);

    auto g_xxzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 186);

    auto g_xxzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 187);

    auto g_xxzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 188);

    auto g_xxzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 189);

    auto g_xxzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 190);

    auto g_xxzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 191);

    auto g_xxzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 192);

    auto g_xxzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 193);

    auto g_xxzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 194);

    auto g_xxzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 195);

    auto g_xxzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 196);

    auto g_xxzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 197);

    auto g_xxzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 198);

    auto g_xxzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 199);

    auto g_xxzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 200);

    auto g_xxzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 201);

    auto g_xxzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 202);

    auto g_xxzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 203);

    auto g_xxzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 204);

    auto g_xxzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 205);

    auto g_xxzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 206);

    auto g_xxzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 207);

    auto g_xxzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 208);

    auto g_xxzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 209);

    auto g_xxzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 210);

    auto g_xxzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 211);

    auto g_xxzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 212);

    auto g_xxzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 213);

    auto g_xxzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 214);

    auto g_xxzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 215);

    auto g_xyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 217);

    auto g_xyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 219);

    auto g_xyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 220);

    auto g_xyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 222);

    auto g_xyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 223);

    auto g_xyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 224);

    auto g_xyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 226);

    auto g_xyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 227);

    auto g_xyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 228);

    auto g_xyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 229);

    auto g_xyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 231);

    auto g_xyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 232);

    auto g_xyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 233);

    auto g_xyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 234);

    auto g_xyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 235);

    auto g_xyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 237);

    auto g_xyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 238);

    auto g_xyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 239);

    auto g_xyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 240);

    auto g_xyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 241);

    auto g_xyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 242);

    auto g_xyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 244);

    auto g_xyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 245);

    auto g_xyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 246);

    auto g_xyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 247);

    auto g_xyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 248);

    auto g_xyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 249);

    auto g_xyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 250);

    auto g_xzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 326);

    auto g_xzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 328);

    auto g_xzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 329);

    auto g_xzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 331);

    auto g_xzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 332);

    auto g_xzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 333);

    auto g_xzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 335);

    auto g_xzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 336);

    auto g_xzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 337);

    auto g_xzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 338);

    auto g_xzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 340);

    auto g_xzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 341);

    auto g_xzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 342);

    auto g_xzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 343);

    auto g_xzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 344);

    auto g_xzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 346);

    auto g_xzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 347);

    auto g_xzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 348);

    auto g_xzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 349);

    auto g_xzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 350);

    auto g_xzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 351);

    auto g_xzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 353);

    auto g_xzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 354);

    auto g_xzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 355);

    auto g_xzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 356);

    auto g_xzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 357);

    auto g_xzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 358);

    auto g_xzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 359);

    auto g_yyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 360);

    auto g_yyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 361);

    auto g_yyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 362);

    auto g_yyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 363);

    auto g_yyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 364);

    auto g_yyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 365);

    auto g_yyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 366);

    auto g_yyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 367);

    auto g_yyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 368);

    auto g_yyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 369);

    auto g_yyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 370);

    auto g_yyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 371);

    auto g_yyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 372);

    auto g_yyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 373);

    auto g_yyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 374);

    auto g_yyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 375);

    auto g_yyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 376);

    auto g_yyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 377);

    auto g_yyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 378);

    auto g_yyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 379);

    auto g_yyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 380);

    auto g_yyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 381);

    auto g_yyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 382);

    auto g_yyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 383);

    auto g_yyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 384);

    auto g_yyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 385);

    auto g_yyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 386);

    auto g_yyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 387);

    auto g_yyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 388);

    auto g_yyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 389);

    auto g_yyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 390);

    auto g_yyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 391);

    auto g_yyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 392);

    auto g_yyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 393);

    auto g_yyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 394);

    auto g_yyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 395);

    auto g_yyyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 398);

    auto g_yyyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 400);

    auto g_yyyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 401);

    auto g_yyyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 403);

    auto g_yyyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 404);

    auto g_yyyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 405);

    auto g_yyyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 407);

    auto g_yyyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 408);

    auto g_yyyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 409);

    auto g_yyyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 410);

    auto g_yyyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 412);

    auto g_yyyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 413);

    auto g_yyyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 414);

    auto g_yyyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 415);

    auto g_yyyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 416);

    auto g_yyyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 418);

    auto g_yyyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 419);

    auto g_yyyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 420);

    auto g_yyyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 421);

    auto g_yyyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 422);

    auto g_yyyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 423);

    auto g_yyyz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 425);

    auto g_yyyz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 426);

    auto g_yyyz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 427);

    auto g_yyyz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 428);

    auto g_yyyz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 429);

    auto g_yyyz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 430);

    auto g_yyyz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 431);

    auto g_yyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 432);

    auto g_yyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 433);

    auto g_yyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 434);

    auto g_yyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 435);

    auto g_yyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 436);

    auto g_yyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 437);

    auto g_yyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 438);

    auto g_yyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 439);

    auto g_yyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 440);

    auto g_yyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 441);

    auto g_yyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 442);

    auto g_yyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 443);

    auto g_yyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 444);

    auto g_yyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 445);

    auto g_yyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 446);

    auto g_yyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 447);

    auto g_yyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 448);

    auto g_yyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 449);

    auto g_yyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 450);

    auto g_yyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 451);

    auto g_yyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 452);

    auto g_yyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 453);

    auto g_yyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 454);

    auto g_yyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 455);

    auto g_yyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 456);

    auto g_yyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 457);

    auto g_yyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 458);

    auto g_yyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 459);

    auto g_yyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 460);

    auto g_yyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 461);

    auto g_yyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 462);

    auto g_yyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 463);

    auto g_yyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 464);

    auto g_yyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 465);

    auto g_yyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 466);

    auto g_yyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 467);

    auto g_yzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 469);

    auto g_yzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 470);

    auto g_yzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 471);

    auto g_yzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 472);

    auto g_yzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 473);

    auto g_yzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 474);

    auto g_yzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 475);

    auto g_yzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 476);

    auto g_yzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 477);

    auto g_yzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 478);

    auto g_yzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 479);

    auto g_yzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 480);

    auto g_yzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 481);

    auto g_yzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 482);

    auto g_yzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 483);

    auto g_yzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 484);

    auto g_yzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 485);

    auto g_yzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 486);

    auto g_yzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 487);

    auto g_yzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 488);

    auto g_yzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 489);

    auto g_yzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 490);

    auto g_yzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 491);

    auto g_yzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 492);

    auto g_yzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 493);

    auto g_yzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 494);

    auto g_yzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 495);

    auto g_yzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 496);

    auto g_yzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 497);

    auto g_yzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 498);

    auto g_yzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 499);

    auto g_yzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 500);

    auto g_yzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 501);

    auto g_yzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 502);

    auto g_yzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 503);

    auto g_zzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 504);

    auto g_zzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 505);

    auto g_zzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 506);

    auto g_zzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 507);

    auto g_zzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 508);

    auto g_zzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 509);

    auto g_zzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 510);

    auto g_zzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 511);

    auto g_zzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 512);

    auto g_zzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 513);

    auto g_zzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 514);

    auto g_zzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 515);

    auto g_zzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 516);

    auto g_zzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 517);

    auto g_zzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 518);

    auto g_zzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 519);

    auto g_zzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 520);

    auto g_zzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 521);

    auto g_zzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 522);

    auto g_zzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 523);

    auto g_zzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 524);

    auto g_zzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 525);

    auto g_zzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 526);

    auto g_zzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 527);

    auto g_zzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 528);

    auto g_zzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 529);

    auto g_zzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 530);

    auto g_zzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 531);

    auto g_zzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 532);

    auto g_zzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 533);

    auto g_zzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 534);

    auto g_zzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 535);

    auto g_zzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 536);

    auto g_zzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 537);

    auto g_zzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 538);

    auto g_zzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 539);

    /// Set up components of auxilary buffer : GSL

    auto g_xxxx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl);

    auto g_xxxx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 1);

    auto g_xxxx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 2);

    auto g_xxxx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 3);

    auto g_xxxx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 4);

    auto g_xxxx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 5);

    auto g_xxxx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 6);

    auto g_xxxx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 7);

    auto g_xxxx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 8);

    auto g_xxxx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 9);

    auto g_xxxx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 10);

    auto g_xxxx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 11);

    auto g_xxxx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 12);

    auto g_xxxx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 13);

    auto g_xxxx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 14);

    auto g_xxxx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 15);

    auto g_xxxx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 16);

    auto g_xxxx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 17);

    auto g_xxxx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 18);

    auto g_xxxx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 19);

    auto g_xxxx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 20);

    auto g_xxxx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 21);

    auto g_xxxx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 22);

    auto g_xxxx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 23);

    auto g_xxxx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 24);

    auto g_xxxx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 25);

    auto g_xxxx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 26);

    auto g_xxxx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 27);

    auto g_xxxx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 28);

    auto g_xxxx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 29);

    auto g_xxxx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 30);

    auto g_xxxx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 31);

    auto g_xxxx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 32);

    auto g_xxxx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 33);

    auto g_xxxx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 34);

    auto g_xxxx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 35);

    auto g_xxxx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 36);

    auto g_xxxx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 37);

    auto g_xxxx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 38);

    auto g_xxxx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 39);

    auto g_xxxx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 40);

    auto g_xxxx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 41);

    auto g_xxxx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 42);

    auto g_xxxx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 43);

    auto g_xxxx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 44);

    auto g_xxxy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 45);

    auto g_xxxy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 46);

    auto g_xxxy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 47);

    auto g_xxxy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 48);

    auto g_xxxy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 50);

    auto g_xxxy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 51);

    auto g_xxxy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 54);

    auto g_xxxy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 55);

    auto g_xxxy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 59);

    auto g_xxxy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 60);

    auto g_xxxy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 65);

    auto g_xxxy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 66);

    auto g_xxxy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 72);

    auto g_xxxy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 73);

    auto g_xxxy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 80);

    auto g_xxxy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 81);

    auto g_xxxz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 90);

    auto g_xxxz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 91);

    auto g_xxxz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 92);

    auto g_xxxz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 93);

    auto g_xxxz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 94);

    auto g_xxxz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 95);

    auto g_xxxz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 96);

    auto g_xxxz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 97);

    auto g_xxxz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 98);

    auto g_xxxz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 99);

    auto g_xxxz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 100);

    auto g_xxxz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 101);

    auto g_xxxz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 102);

    auto g_xxxz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 103);

    auto g_xxxz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 104);

    auto g_xxxz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 105);

    auto g_xxxz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 106);

    auto g_xxxz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 107);

    auto g_xxxz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 108);

    auto g_xxxz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 109);

    auto g_xxxz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 110);

    auto g_xxxz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 111);

    auto g_xxxz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 112);

    auto g_xxxz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 113);

    auto g_xxxz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 114);

    auto g_xxxz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 115);

    auto g_xxxz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 116);

    auto g_xxxz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 117);

    auto g_xxxz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 118);

    auto g_xxxz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 119);

    auto g_xxxz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 120);

    auto g_xxxz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 121);

    auto g_xxxz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 122);

    auto g_xxxz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 123);

    auto g_xxxz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 124);

    auto g_xxxz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 125);

    auto g_xxxz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 127);

    auto g_xxxz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 128);

    auto g_xxxz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 129);

    auto g_xxxz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 130);

    auto g_xxxz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 131);

    auto g_xxxz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 132);

    auto g_xxxz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 133);

    auto g_xxxz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 134);

    auto g_xxyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 135);

    auto g_xxyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 136);

    auto g_xxyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 137);

    auto g_xxyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 138);

    auto g_xxyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 139);

    auto g_xxyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 140);

    auto g_xxyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 141);

    auto g_xxyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 142);

    auto g_xxyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 143);

    auto g_xxyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 144);

    auto g_xxyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 145);

    auto g_xxyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 146);

    auto g_xxyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 147);

    auto g_xxyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 148);

    auto g_xxyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 149);

    auto g_xxyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 150);

    auto g_xxyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 151);

    auto g_xxyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 152);

    auto g_xxyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 153);

    auto g_xxyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 154);

    auto g_xxyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 155);

    auto g_xxyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 156);

    auto g_xxyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 157);

    auto g_xxyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 158);

    auto g_xxyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 159);

    auto g_xxyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 160);

    auto g_xxyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 161);

    auto g_xxyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 162);

    auto g_xxyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 163);

    auto g_xxyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 164);

    auto g_xxyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 165);

    auto g_xxyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 166);

    auto g_xxyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 167);

    auto g_xxyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 168);

    auto g_xxyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 169);

    auto g_xxyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 170);

    auto g_xxyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 171);

    auto g_xxyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 172);

    auto g_xxyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 173);

    auto g_xxyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 174);

    auto g_xxyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 175);

    auto g_xxyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 176);

    auto g_xxyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 177);

    auto g_xxyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 178);

    auto g_xxyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 179);

    auto g_xxzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 225);

    auto g_xxzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 226);

    auto g_xxzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 227);

    auto g_xxzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 228);

    auto g_xxzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 229);

    auto g_xxzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 230);

    auto g_xxzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 231);

    auto g_xxzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 232);

    auto g_xxzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 233);

    auto g_xxzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 234);

    auto g_xxzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 235);

    auto g_xxzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 236);

    auto g_xxzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 237);

    auto g_xxzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 238);

    auto g_xxzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 239);

    auto g_xxzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 240);

    auto g_xxzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 241);

    auto g_xxzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 242);

    auto g_xxzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 243);

    auto g_xxzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 244);

    auto g_xxzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 245);

    auto g_xxzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 246);

    auto g_xxzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 247);

    auto g_xxzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 248);

    auto g_xxzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 249);

    auto g_xxzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 250);

    auto g_xxzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 251);

    auto g_xxzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 252);

    auto g_xxzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 253);

    auto g_xxzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 254);

    auto g_xxzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 255);

    auto g_xxzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 256);

    auto g_xxzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 257);

    auto g_xxzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 258);

    auto g_xxzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 259);

    auto g_xxzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 260);

    auto g_xxzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 261);

    auto g_xxzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 262);

    auto g_xxzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 263);

    auto g_xxzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 264);

    auto g_xxzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 265);

    auto g_xxzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 266);

    auto g_xxzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 267);

    auto g_xxzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 268);

    auto g_xxzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 269);

    auto g_xyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 270);

    auto g_xyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 271);

    auto g_xyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 273);

    auto g_xyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 274);

    auto g_xyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 276);

    auto g_xyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 277);

    auto g_xyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 278);

    auto g_xyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 280);

    auto g_xyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 281);

    auto g_xyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 282);

    auto g_xyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 283);

    auto g_xyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 285);

    auto g_xyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 286);

    auto g_xyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 287);

    auto g_xyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 288);

    auto g_xyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 289);

    auto g_xyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 291);

    auto g_xyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 292);

    auto g_xyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 293);

    auto g_xyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 294);

    auto g_xyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 295);

    auto g_xyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 296);

    auto g_xyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 298);

    auto g_xyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 299);

    auto g_xyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 300);

    auto g_xyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 301);

    auto g_xyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 302);

    auto g_xyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 303);

    auto g_xyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 304);

    auto g_xyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 306);

    auto g_xyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 307);

    auto g_xyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 308);

    auto g_xyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 309);

    auto g_xyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 310);

    auto g_xyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 311);

    auto g_xyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 312);

    auto g_xyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 313);

    auto g_xyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 314);

    auto g_xzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 405);

    auto g_xzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 407);

    auto g_xzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 409);

    auto g_xzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 410);

    auto g_xzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 412);

    auto g_xzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 413);

    auto g_xzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 414);

    auto g_xzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 416);

    auto g_xzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 417);

    auto g_xzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 418);

    auto g_xzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 419);

    auto g_xzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 421);

    auto g_xzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 422);

    auto g_xzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 423);

    auto g_xzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 424);

    auto g_xzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 425);

    auto g_xzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 427);

    auto g_xzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 428);

    auto g_xzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 429);

    auto g_xzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 430);

    auto g_xzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 431);

    auto g_xzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 432);

    auto g_xzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 434);

    auto g_xzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 435);

    auto g_xzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 436);

    auto g_xzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 437);

    auto g_xzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 438);

    auto g_xzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 439);

    auto g_xzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 440);

    auto g_xzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 441);

    auto g_xzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 442);

    auto g_xzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 443);

    auto g_xzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 444);

    auto g_xzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 445);

    auto g_xzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 446);

    auto g_xzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 447);

    auto g_xzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 448);

    auto g_xzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 449);

    auto g_yyyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 450);

    auto g_yyyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 451);

    auto g_yyyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 452);

    auto g_yyyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 453);

    auto g_yyyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 454);

    auto g_yyyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 455);

    auto g_yyyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 456);

    auto g_yyyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 457);

    auto g_yyyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 458);

    auto g_yyyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 459);

    auto g_yyyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 460);

    auto g_yyyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 461);

    auto g_yyyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 462);

    auto g_yyyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 463);

    auto g_yyyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 464);

    auto g_yyyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 465);

    auto g_yyyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 466);

    auto g_yyyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 467);

    auto g_yyyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 468);

    auto g_yyyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 469);

    auto g_yyyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 470);

    auto g_yyyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 471);

    auto g_yyyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 472);

    auto g_yyyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 473);

    auto g_yyyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 474);

    auto g_yyyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 475);

    auto g_yyyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 476);

    auto g_yyyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 477);

    auto g_yyyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 478);

    auto g_yyyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 479);

    auto g_yyyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 480);

    auto g_yyyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 481);

    auto g_yyyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 482);

    auto g_yyyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 483);

    auto g_yyyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 484);

    auto g_yyyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 485);

    auto g_yyyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 486);

    auto g_yyyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 487);

    auto g_yyyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 488);

    auto g_yyyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 489);

    auto g_yyyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 490);

    auto g_yyyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 491);

    auto g_yyyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 492);

    auto g_yyyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 493);

    auto g_yyyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 494);

    auto g_yyyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 496);

    auto g_yyyz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 497);

    auto g_yyyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 498);

    auto g_yyyz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 499);

    auto g_yyyz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 500);

    auto g_yyyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 501);

    auto g_yyyz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 502);

    auto g_yyyz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 503);

    auto g_yyyz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 504);

    auto g_yyyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 505);

    auto g_yyyz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 506);

    auto g_yyyz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 507);

    auto g_yyyz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 508);

    auto g_yyyz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 509);

    auto g_yyyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 510);

    auto g_yyyz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 511);

    auto g_yyyz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 512);

    auto g_yyyz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 513);

    auto g_yyyz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 514);

    auto g_yyyz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 515);

    auto g_yyyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 516);

    auto g_yyyz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 517);

    auto g_yyyz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 518);

    auto g_yyyz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 519);

    auto g_yyyz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 520);

    auto g_yyyz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 521);

    auto g_yyyz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 522);

    auto g_yyyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 523);

    auto g_yyyz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 524);

    auto g_yyyz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 525);

    auto g_yyyz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 526);

    auto g_yyyz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 527);

    auto g_yyyz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 528);

    auto g_yyyz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 529);

    auto g_yyyz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 530);

    auto g_yyyz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 531);

    auto g_yyyz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 532);

    auto g_yyyz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 533);

    auto g_yyyz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 534);

    auto g_yyyz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 535);

    auto g_yyyz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 536);

    auto g_yyyz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 537);

    auto g_yyyz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 538);

    auto g_yyyz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 539);

    auto g_yyzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 540);

    auto g_yyzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 541);

    auto g_yyzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 542);

    auto g_yyzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 543);

    auto g_yyzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 544);

    auto g_yyzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 545);

    auto g_yyzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 546);

    auto g_yyzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 547);

    auto g_yyzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 548);

    auto g_yyzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 549);

    auto g_yyzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 550);

    auto g_yyzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 551);

    auto g_yyzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 552);

    auto g_yyzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 553);

    auto g_yyzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 554);

    auto g_yyzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 555);

    auto g_yyzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 556);

    auto g_yyzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 557);

    auto g_yyzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 558);

    auto g_yyzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 559);

    auto g_yyzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 560);

    auto g_yyzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 561);

    auto g_yyzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 562);

    auto g_yyzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 563);

    auto g_yyzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 564);

    auto g_yyzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 565);

    auto g_yyzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 566);

    auto g_yyzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 567);

    auto g_yyzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 568);

    auto g_yyzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 569);

    auto g_yyzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 570);

    auto g_yyzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 571);

    auto g_yyzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 572);

    auto g_yyzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 573);

    auto g_yyzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 574);

    auto g_yyzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 575);

    auto g_yyzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 576);

    auto g_yyzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 577);

    auto g_yyzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 578);

    auto g_yyzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 579);

    auto g_yyzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 580);

    auto g_yyzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 581);

    auto g_yyzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 582);

    auto g_yyzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 583);

    auto g_yyzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 584);

    auto g_yzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 585);

    auto g_yzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 586);

    auto g_yzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 587);

    auto g_yzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 588);

    auto g_yzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 589);

    auto g_yzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 590);

    auto g_yzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 591);

    auto g_yzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 592);

    auto g_yzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 593);

    auto g_yzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 594);

    auto g_yzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 595);

    auto g_yzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 596);

    auto g_yzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 597);

    auto g_yzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 598);

    auto g_yzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 599);

    auto g_yzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 600);

    auto g_yzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 601);

    auto g_yzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 602);

    auto g_yzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 603);

    auto g_yzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 604);

    auto g_yzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 605);

    auto g_yzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 606);

    auto g_yzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 607);

    auto g_yzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 608);

    auto g_yzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 609);

    auto g_yzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 610);

    auto g_yzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 611);

    auto g_yzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 612);

    auto g_yzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 613);

    auto g_yzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 614);

    auto g_yzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 615);

    auto g_yzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 616);

    auto g_yzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 617);

    auto g_yzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 618);

    auto g_yzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 619);

    auto g_yzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 620);

    auto g_yzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 621);

    auto g_yzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 622);

    auto g_yzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 623);

    auto g_yzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 624);

    auto g_yzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 625);

    auto g_yzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 626);

    auto g_yzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 627);

    auto g_yzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 628);

    auto g_yzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 629);

    auto g_zzzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_gsl + 630);

    auto g_zzzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_gsl + 631);

    auto g_zzzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_gsl + 632);

    auto g_zzzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_gsl + 633);

    auto g_zzzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_gsl + 634);

    auto g_zzzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_gsl + 635);

    auto g_zzzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_gsl + 636);

    auto g_zzzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_gsl + 637);

    auto g_zzzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_gsl + 638);

    auto g_zzzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_gsl + 639);

    auto g_zzzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_gsl + 640);

    auto g_zzzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_gsl + 641);

    auto g_zzzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_gsl + 642);

    auto g_zzzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_gsl + 643);

    auto g_zzzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_gsl + 644);

    auto g_zzzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 645);

    auto g_zzzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 646);

    auto g_zzzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 647);

    auto g_zzzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 648);

    auto g_zzzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 649);

    auto g_zzzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 650);

    auto g_zzzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 651);

    auto g_zzzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 652);

    auto g_zzzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 653);

    auto g_zzzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 654);

    auto g_zzzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 655);

    auto g_zzzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 656);

    auto g_zzzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 657);

    auto g_zzzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 658);

    auto g_zzzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 659);

    auto g_zzzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 660);

    auto g_zzzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 661);

    auto g_zzzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 662);

    auto g_zzzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 663);

    auto g_zzzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 664);

    auto g_zzzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 665);

    auto g_zzzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_gsl + 666);

    auto g_zzzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_gsl + 667);

    auto g_zzzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_gsl + 668);

    auto g_zzzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_gsl + 669);

    auto g_zzzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_gsl + 670);

    auto g_zzzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 671);

    auto g_zzzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 672);

    auto g_zzzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 673);

    auto g_zzzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_gsl + 674);

    /// Set up 0-45 components of targeted buffer : HSL

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

    #pragma omp simd aligned(g_xxx_0_xxxxxxxx_0, g_xxx_0_xxxxxxxx_1, g_xxx_0_xxxxxxxy_0, g_xxx_0_xxxxxxxy_1, g_xxx_0_xxxxxxxz_0, g_xxx_0_xxxxxxxz_1, g_xxx_0_xxxxxxyy_0, g_xxx_0_xxxxxxyy_1, g_xxx_0_xxxxxxyz_0, g_xxx_0_xxxxxxyz_1, g_xxx_0_xxxxxxzz_0, g_xxx_0_xxxxxxzz_1, g_xxx_0_xxxxxyyy_0, g_xxx_0_xxxxxyyy_1, g_xxx_0_xxxxxyyz_0, g_xxx_0_xxxxxyyz_1, g_xxx_0_xxxxxyzz_0, g_xxx_0_xxxxxyzz_1, g_xxx_0_xxxxxzzz_0, g_xxx_0_xxxxxzzz_1, g_xxx_0_xxxxyyyy_0, g_xxx_0_xxxxyyyy_1, g_xxx_0_xxxxyyyz_0, g_xxx_0_xxxxyyyz_1, g_xxx_0_xxxxyyzz_0, g_xxx_0_xxxxyyzz_1, g_xxx_0_xxxxyzzz_0, g_xxx_0_xxxxyzzz_1, g_xxx_0_xxxxzzzz_0, g_xxx_0_xxxxzzzz_1, g_xxx_0_xxxyyyyy_0, g_xxx_0_xxxyyyyy_1, g_xxx_0_xxxyyyyz_0, g_xxx_0_xxxyyyyz_1, g_xxx_0_xxxyyyzz_0, g_xxx_0_xxxyyyzz_1, g_xxx_0_xxxyyzzz_0, g_xxx_0_xxxyyzzz_1, g_xxx_0_xxxyzzzz_0, g_xxx_0_xxxyzzzz_1, g_xxx_0_xxxzzzzz_0, g_xxx_0_xxxzzzzz_1, g_xxx_0_xxyyyyyy_0, g_xxx_0_xxyyyyyy_1, g_xxx_0_xxyyyyyz_0, g_xxx_0_xxyyyyyz_1, g_xxx_0_xxyyyyzz_0, g_xxx_0_xxyyyyzz_1, g_xxx_0_xxyyyzzz_0, g_xxx_0_xxyyyzzz_1, g_xxx_0_xxyyzzzz_0, g_xxx_0_xxyyzzzz_1, g_xxx_0_xxyzzzzz_0, g_xxx_0_xxyzzzzz_1, g_xxx_0_xxzzzzzz_0, g_xxx_0_xxzzzzzz_1, g_xxx_0_xyyyyyyy_0, g_xxx_0_xyyyyyyy_1, g_xxx_0_xyyyyyyz_0, g_xxx_0_xyyyyyyz_1, g_xxx_0_xyyyyyzz_0, g_xxx_0_xyyyyyzz_1, g_xxx_0_xyyyyzzz_0, g_xxx_0_xyyyyzzz_1, g_xxx_0_xyyyzzzz_0, g_xxx_0_xyyyzzzz_1, g_xxx_0_xyyzzzzz_0, g_xxx_0_xyyzzzzz_1, g_xxx_0_xyzzzzzz_0, g_xxx_0_xyzzzzzz_1, g_xxx_0_xzzzzzzz_0, g_xxx_0_xzzzzzzz_1, g_xxx_0_yyyyyyyy_0, g_xxx_0_yyyyyyyy_1, g_xxx_0_yyyyyyyz_0, g_xxx_0_yyyyyyyz_1, g_xxx_0_yyyyyyzz_0, g_xxx_0_yyyyyyzz_1, g_xxx_0_yyyyyzzz_0, g_xxx_0_yyyyyzzz_1, g_xxx_0_yyyyzzzz_0, g_xxx_0_yyyyzzzz_1, g_xxx_0_yyyzzzzz_0, g_xxx_0_yyyzzzzz_1, g_xxx_0_yyzzzzzz_0, g_xxx_0_yyzzzzzz_1, g_xxx_0_yzzzzzzz_0, g_xxx_0_yzzzzzzz_1, g_xxx_0_zzzzzzzz_0, g_xxx_0_zzzzzzzz_1, g_xxxx_0_xxxxxxx_1, g_xxxx_0_xxxxxxxx_1, g_xxxx_0_xxxxxxxy_1, g_xxxx_0_xxxxxxxz_1, g_xxxx_0_xxxxxxy_1, g_xxxx_0_xxxxxxyy_1, g_xxxx_0_xxxxxxyz_1, g_xxxx_0_xxxxxxz_1, g_xxxx_0_xxxxxxzz_1, g_xxxx_0_xxxxxyy_1, g_xxxx_0_xxxxxyyy_1, g_xxxx_0_xxxxxyyz_1, g_xxxx_0_xxxxxyz_1, g_xxxx_0_xxxxxyzz_1, g_xxxx_0_xxxxxzz_1, g_xxxx_0_xxxxxzzz_1, g_xxxx_0_xxxxyyy_1, g_xxxx_0_xxxxyyyy_1, g_xxxx_0_xxxxyyyz_1, g_xxxx_0_xxxxyyz_1, g_xxxx_0_xxxxyyzz_1, g_xxxx_0_xxxxyzz_1, g_xxxx_0_xxxxyzzz_1, g_xxxx_0_xxxxzzz_1, g_xxxx_0_xxxxzzzz_1, g_xxxx_0_xxxyyyy_1, g_xxxx_0_xxxyyyyy_1, g_xxxx_0_xxxyyyyz_1, g_xxxx_0_xxxyyyz_1, g_xxxx_0_xxxyyyzz_1, g_xxxx_0_xxxyyzz_1, g_xxxx_0_xxxyyzzz_1, g_xxxx_0_xxxyzzz_1, g_xxxx_0_xxxyzzzz_1, g_xxxx_0_xxxzzzz_1, g_xxxx_0_xxxzzzzz_1, g_xxxx_0_xxyyyyy_1, g_xxxx_0_xxyyyyyy_1, g_xxxx_0_xxyyyyyz_1, g_xxxx_0_xxyyyyz_1, g_xxxx_0_xxyyyyzz_1, g_xxxx_0_xxyyyzz_1, g_xxxx_0_xxyyyzzz_1, g_xxxx_0_xxyyzzz_1, g_xxxx_0_xxyyzzzz_1, g_xxxx_0_xxyzzzz_1, g_xxxx_0_xxyzzzzz_1, g_xxxx_0_xxzzzzz_1, g_xxxx_0_xxzzzzzz_1, g_xxxx_0_xyyyyyy_1, g_xxxx_0_xyyyyyyy_1, g_xxxx_0_xyyyyyyz_1, g_xxxx_0_xyyyyyz_1, g_xxxx_0_xyyyyyzz_1, g_xxxx_0_xyyyyzz_1, g_xxxx_0_xyyyyzzz_1, g_xxxx_0_xyyyzzz_1, g_xxxx_0_xyyyzzzz_1, g_xxxx_0_xyyzzzz_1, g_xxxx_0_xyyzzzzz_1, g_xxxx_0_xyzzzzz_1, g_xxxx_0_xyzzzzzz_1, g_xxxx_0_xzzzzzz_1, g_xxxx_0_xzzzzzzz_1, g_xxxx_0_yyyyyyy_1, g_xxxx_0_yyyyyyyy_1, g_xxxx_0_yyyyyyyz_1, g_xxxx_0_yyyyyyz_1, g_xxxx_0_yyyyyyzz_1, g_xxxx_0_yyyyyzz_1, g_xxxx_0_yyyyyzzz_1, g_xxxx_0_yyyyzzz_1, g_xxxx_0_yyyyzzzz_1, g_xxxx_0_yyyzzzz_1, g_xxxx_0_yyyzzzzz_1, g_xxxx_0_yyzzzzz_1, g_xxxx_0_yyzzzzzz_1, g_xxxx_0_yzzzzzz_1, g_xxxx_0_yzzzzzzz_1, g_xxxx_0_zzzzzzz_1, g_xxxx_0_zzzzzzzz_1, g_xxxxx_0_xxxxxxxx_0, g_xxxxx_0_xxxxxxxy_0, g_xxxxx_0_xxxxxxxz_0, g_xxxxx_0_xxxxxxyy_0, g_xxxxx_0_xxxxxxyz_0, g_xxxxx_0_xxxxxxzz_0, g_xxxxx_0_xxxxxyyy_0, g_xxxxx_0_xxxxxyyz_0, g_xxxxx_0_xxxxxyzz_0, g_xxxxx_0_xxxxxzzz_0, g_xxxxx_0_xxxxyyyy_0, g_xxxxx_0_xxxxyyyz_0, g_xxxxx_0_xxxxyyzz_0, g_xxxxx_0_xxxxyzzz_0, g_xxxxx_0_xxxxzzzz_0, g_xxxxx_0_xxxyyyyy_0, g_xxxxx_0_xxxyyyyz_0, g_xxxxx_0_xxxyyyzz_0, g_xxxxx_0_xxxyyzzz_0, g_xxxxx_0_xxxyzzzz_0, g_xxxxx_0_xxxzzzzz_0, g_xxxxx_0_xxyyyyyy_0, g_xxxxx_0_xxyyyyyz_0, g_xxxxx_0_xxyyyyzz_0, g_xxxxx_0_xxyyyzzz_0, g_xxxxx_0_xxyyzzzz_0, g_xxxxx_0_xxyzzzzz_0, g_xxxxx_0_xxzzzzzz_0, g_xxxxx_0_xyyyyyyy_0, g_xxxxx_0_xyyyyyyz_0, g_xxxxx_0_xyyyyyzz_0, g_xxxxx_0_xyyyyzzz_0, g_xxxxx_0_xyyyzzzz_0, g_xxxxx_0_xyyzzzzz_0, g_xxxxx_0_xyzzzzzz_0, g_xxxxx_0_xzzzzzzz_0, g_xxxxx_0_yyyyyyyy_0, g_xxxxx_0_yyyyyyyz_0, g_xxxxx_0_yyyyyyzz_0, g_xxxxx_0_yyyyyzzz_0, g_xxxxx_0_yyyyzzzz_0, g_xxxxx_0_yyyzzzzz_0, g_xxxxx_0_yyzzzzzz_0, g_xxxxx_0_yzzzzzzz_0, g_xxxxx_0_zzzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxx_0_xxxxxxxx_0[i] = 4.0 * g_xxx_0_xxxxxxxx_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxxxx_1[i] * fz_be_0 + 8.0 * g_xxxx_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxxx_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxxxy_0[i] = 4.0 * g_xxx_0_xxxxxxxy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xxxx_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxxxz_0[i] = 4.0 * g_xxx_0_xxxxxxxz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xxxx_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxxyy_0[i] = 4.0 * g_xxx_0_xxxxxxyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xxxx_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxxyz_0[i] = 4.0 * g_xxx_0_xxxxxxyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxxx_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxxzz_0[i] = 4.0 * g_xxx_0_xxxxxxzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xxxx_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxyyy_0[i] = 4.0 * g_xxx_0_xxxxxyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xxxx_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxyyz_0[i] = 4.0 * g_xxx_0_xxxxxyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxxx_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxyzz_0[i] = 4.0 * g_xxx_0_xxxxxyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxxx_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxxzzz_0[i] = 4.0 * g_xxx_0_xxxxxzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xxxx_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxyyyy_0[i] = 4.0 * g_xxx_0_xxxxyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxxxx_0_xxxxyyyz_0[i] = 4.0 * g_xxx_0_xxxxyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxyyzz_0[i] = 4.0 * g_xxx_0_xxxxyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxyzzz_0[i] = 4.0 * g_xxx_0_xxxxyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxxzzzz_0[i] = 4.0 * g_xxx_0_xxxxzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xxxx_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxyyyyy_0[i] = 4.0 * g_xxx_0_xxxyyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxxxx_0_xxxyyyyz_0[i] = 4.0 * g_xxx_0_xxxyyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxxx_0_xxxyyyzz_0[i] = 4.0 * g_xxx_0_xxxyyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxyyzzz_0[i] = 4.0 * g_xxx_0_xxxyyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxyzzzz_0[i] = 4.0 * g_xxx_0_xxxyzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxxzzzzz_0[i] = 4.0 * g_xxx_0_xxxzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xxxx_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxyyyyyy_0[i] = 4.0 * g_xxx_0_xxyyyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxxxx_0_xxyyyyyz_0[i] = 4.0 * g_xxx_0_xxyyyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxxx_0_xxyyyyzz_0[i] = 4.0 * g_xxx_0_xxyyyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxxx_0_xxyyyzzz_0[i] = 4.0 * g_xxx_0_xxyyyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxyyzzzz_0[i] = 4.0 * g_xxx_0_xxyyzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxyzzzzz_0[i] = 4.0 * g_xxx_0_xxyzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xxzzzzzz_0[i] = 4.0 * g_xxx_0_xxzzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xyyyyyyy_0[i] = 4.0 * g_xxx_0_xyyyyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxx_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxxxx_0_xyyyyyyz_0[i] = 4.0 * g_xxx_0_xyyyyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyyyyz_1[i] * fz_be_0 + g_xxxx_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxxx_0_xyyyyyzz_0[i] = 4.0 * g_xxx_0_xyyyyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyyyzz_1[i] * fz_be_0 + g_xxxx_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxxx_0_xyyyyzzz_0[i] = 4.0 * g_xxx_0_xyyyyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyyzzz_1[i] * fz_be_0 + g_xxxx_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxxx_0_xyyyzzzz_0[i] = 4.0 * g_xxx_0_xyyyzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyyzzzz_1[i] * fz_be_0 + g_xxxx_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xyyzzzzz_0[i] = 4.0 * g_xxx_0_xyyzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyyzzzzz_1[i] * fz_be_0 + g_xxxx_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xyzzzzzz_0[i] = 4.0 * g_xxx_0_xyzzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xyzzzzzz_1[i] * fz_be_0 + g_xxxx_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_xzzzzzzz_0[i] = 4.0 * g_xxx_0_xzzzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxx_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_yyyyyyyy_0[i] = 4.0 * g_xxx_0_yyyyyyyy_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyyyyy_1[i] * fz_be_0 + g_xxxx_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxxx_0_yyyyyyyz_0[i] = 4.0 * g_xxx_0_yyyyyyyz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyyyyz_1[i] * fz_be_0 + g_xxxx_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxxx_0_yyyyyyzz_0[i] = 4.0 * g_xxx_0_yyyyyyzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyyyzz_1[i] * fz_be_0 + g_xxxx_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxxx_0_yyyyyzzz_0[i] = 4.0 * g_xxx_0_yyyyyzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyyzzz_1[i] * fz_be_0 + g_xxxx_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxxx_0_yyyyzzzz_0[i] = 4.0 * g_xxx_0_yyyyzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyyzzzz_1[i] * fz_be_0 + g_xxxx_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxxx_0_yyyzzzzz_0[i] = 4.0 * g_xxx_0_yyyzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyyzzzzz_1[i] * fz_be_0 + g_xxxx_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_yyzzzzzz_0[i] = 4.0 * g_xxx_0_yyzzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yyzzzzzz_1[i] * fz_be_0 + g_xxxx_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_yzzzzzzz_0[i] = 4.0 * g_xxx_0_yzzzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_yzzzzzzz_1[i] * fz_be_0 + g_xxxx_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxxx_0_zzzzzzzz_0[i] = 4.0 * g_xxx_0_zzzzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_0_zzzzzzzz_1[i] * fz_be_0 + g_xxxx_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 45-90 components of targeted buffer : HSL

    auto g_xxxxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 45);

    auto g_xxxxy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 46);

    auto g_xxxxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 47);

    auto g_xxxxy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 48);

    auto g_xxxxy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 49);

    auto g_xxxxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 50);

    auto g_xxxxy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 51);

    auto g_xxxxy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 52);

    auto g_xxxxy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 53);

    auto g_xxxxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 54);

    auto g_xxxxy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 55);

    auto g_xxxxy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 56);

    auto g_xxxxy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 57);

    auto g_xxxxy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 58);

    auto g_xxxxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 59);

    auto g_xxxxy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 60);

    auto g_xxxxy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 61);

    auto g_xxxxy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 62);

    auto g_xxxxy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 63);

    auto g_xxxxy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 64);

    auto g_xxxxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 65);

    auto g_xxxxy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 66);

    auto g_xxxxy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 67);

    auto g_xxxxy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 68);

    auto g_xxxxy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 69);

    auto g_xxxxy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 70);

    auto g_xxxxy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 71);

    auto g_xxxxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 72);

    auto g_xxxxy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 73);

    auto g_xxxxy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 74);

    auto g_xxxxy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 75);

    auto g_xxxxy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 76);

    auto g_xxxxy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 77);

    auto g_xxxxy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 78);

    auto g_xxxxy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 79);

    auto g_xxxxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 80);

    auto g_xxxxy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 81);

    auto g_xxxxy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 82);

    auto g_xxxxy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 83);

    auto g_xxxxy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 84);

    auto g_xxxxy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 85);

    auto g_xxxxy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 86);

    auto g_xxxxy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 87);

    auto g_xxxxy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 88);

    auto g_xxxxy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 89);

    #pragma omp simd aligned(g_xxxx_0_xxxxxxx_1, g_xxxx_0_xxxxxxxx_1, g_xxxx_0_xxxxxxxy_1, g_xxxx_0_xxxxxxxz_1, g_xxxx_0_xxxxxxy_1, g_xxxx_0_xxxxxxyy_1, g_xxxx_0_xxxxxxyz_1, g_xxxx_0_xxxxxxz_1, g_xxxx_0_xxxxxxzz_1, g_xxxx_0_xxxxxyy_1, g_xxxx_0_xxxxxyyy_1, g_xxxx_0_xxxxxyyz_1, g_xxxx_0_xxxxxyz_1, g_xxxx_0_xxxxxyzz_1, g_xxxx_0_xxxxxzz_1, g_xxxx_0_xxxxxzzz_1, g_xxxx_0_xxxxyyy_1, g_xxxx_0_xxxxyyyy_1, g_xxxx_0_xxxxyyyz_1, g_xxxx_0_xxxxyyz_1, g_xxxx_0_xxxxyyzz_1, g_xxxx_0_xxxxyzz_1, g_xxxx_0_xxxxyzzz_1, g_xxxx_0_xxxxzzz_1, g_xxxx_0_xxxxzzzz_1, g_xxxx_0_xxxyyyy_1, g_xxxx_0_xxxyyyyy_1, g_xxxx_0_xxxyyyyz_1, g_xxxx_0_xxxyyyz_1, g_xxxx_0_xxxyyyzz_1, g_xxxx_0_xxxyyzz_1, g_xxxx_0_xxxyyzzz_1, g_xxxx_0_xxxyzzz_1, g_xxxx_0_xxxyzzzz_1, g_xxxx_0_xxxzzzz_1, g_xxxx_0_xxxzzzzz_1, g_xxxx_0_xxyyyyy_1, g_xxxx_0_xxyyyyyy_1, g_xxxx_0_xxyyyyyz_1, g_xxxx_0_xxyyyyz_1, g_xxxx_0_xxyyyyzz_1, g_xxxx_0_xxyyyzz_1, g_xxxx_0_xxyyyzzz_1, g_xxxx_0_xxyyzzz_1, g_xxxx_0_xxyyzzzz_1, g_xxxx_0_xxyzzzz_1, g_xxxx_0_xxyzzzzz_1, g_xxxx_0_xxzzzzz_1, g_xxxx_0_xxzzzzzz_1, g_xxxx_0_xyyyyyy_1, g_xxxx_0_xyyyyyyy_1, g_xxxx_0_xyyyyyyz_1, g_xxxx_0_xyyyyyz_1, g_xxxx_0_xyyyyyzz_1, g_xxxx_0_xyyyyzz_1, g_xxxx_0_xyyyyzzz_1, g_xxxx_0_xyyyzzz_1, g_xxxx_0_xyyyzzzz_1, g_xxxx_0_xyyzzzz_1, g_xxxx_0_xyyzzzzz_1, g_xxxx_0_xyzzzzz_1, g_xxxx_0_xyzzzzzz_1, g_xxxx_0_xzzzzzz_1, g_xxxx_0_xzzzzzzz_1, g_xxxx_0_yyyyyyy_1, g_xxxx_0_yyyyyyyy_1, g_xxxx_0_yyyyyyyz_1, g_xxxx_0_yyyyyyz_1, g_xxxx_0_yyyyyyzz_1, g_xxxx_0_yyyyyzz_1, g_xxxx_0_yyyyyzzz_1, g_xxxx_0_yyyyzzz_1, g_xxxx_0_yyyyzzzz_1, g_xxxx_0_yyyzzzz_1, g_xxxx_0_yyyzzzzz_1, g_xxxx_0_yyzzzzz_1, g_xxxx_0_yyzzzzzz_1, g_xxxx_0_yzzzzzz_1, g_xxxx_0_yzzzzzzz_1, g_xxxx_0_zzzzzzz_1, g_xxxx_0_zzzzzzzz_1, g_xxxxy_0_xxxxxxxx_0, g_xxxxy_0_xxxxxxxy_0, g_xxxxy_0_xxxxxxxz_0, g_xxxxy_0_xxxxxxyy_0, g_xxxxy_0_xxxxxxyz_0, g_xxxxy_0_xxxxxxzz_0, g_xxxxy_0_xxxxxyyy_0, g_xxxxy_0_xxxxxyyz_0, g_xxxxy_0_xxxxxyzz_0, g_xxxxy_0_xxxxxzzz_0, g_xxxxy_0_xxxxyyyy_0, g_xxxxy_0_xxxxyyyz_0, g_xxxxy_0_xxxxyyzz_0, g_xxxxy_0_xxxxyzzz_0, g_xxxxy_0_xxxxzzzz_0, g_xxxxy_0_xxxyyyyy_0, g_xxxxy_0_xxxyyyyz_0, g_xxxxy_0_xxxyyyzz_0, g_xxxxy_0_xxxyyzzz_0, g_xxxxy_0_xxxyzzzz_0, g_xxxxy_0_xxxzzzzz_0, g_xxxxy_0_xxyyyyyy_0, g_xxxxy_0_xxyyyyyz_0, g_xxxxy_0_xxyyyyzz_0, g_xxxxy_0_xxyyyzzz_0, g_xxxxy_0_xxyyzzzz_0, g_xxxxy_0_xxyzzzzz_0, g_xxxxy_0_xxzzzzzz_0, g_xxxxy_0_xyyyyyyy_0, g_xxxxy_0_xyyyyyyz_0, g_xxxxy_0_xyyyyyzz_0, g_xxxxy_0_xyyyyzzz_0, g_xxxxy_0_xyyyzzzz_0, g_xxxxy_0_xyyzzzzz_0, g_xxxxy_0_xyzzzzzz_0, g_xxxxy_0_xzzzzzzz_0, g_xxxxy_0_yyyyyyyy_0, g_xxxxy_0_yyyyyyyz_0, g_xxxxy_0_yyyyyyzz_0, g_xxxxy_0_yyyyyzzz_0, g_xxxxy_0_yyyyzzzz_0, g_xxxxy_0_yyyzzzzz_0, g_xxxxy_0_yyzzzzzz_0, g_xxxxy_0_yzzzzzzz_0, g_xxxxy_0_zzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxy_0_xxxxxxxx_0[i] = g_xxxx_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxxxy_0[i] = g_xxxx_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxxy_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxxxz_0[i] = g_xxxx_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxxyy_0[i] = 2.0 * g_xxxx_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxyy_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxxyz_0[i] = g_xxxx_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxxzz_0[i] = g_xxxx_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxyyy_0[i] = 3.0 * g_xxxx_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyyy_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxyyz_0[i] = 2.0 * g_xxxx_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxyzz_0[i] = g_xxxx_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxxzzz_0[i] = g_xxxx_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxyyyy_0[i] = 4.0 * g_xxxx_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyyy_1[i] * wa_y[i];

        g_xxxxy_0_xxxxyyyz_0[i] = 3.0 * g_xxxx_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxyyzz_0[i] = 2.0 * g_xxxx_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxyzzz_0[i] = g_xxxx_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxxzzzz_0[i] = g_xxxx_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxyyyyy_0[i] = 5.0 * g_xxxx_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyyy_1[i] * wa_y[i];

        g_xxxxy_0_xxxyyyyz_0[i] = 4.0 * g_xxxx_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxxxy_0_xxxyyyzz_0[i] = 3.0 * g_xxxx_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxyyzzz_0[i] = 2.0 * g_xxxx_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxyzzzz_0[i] = g_xxxx_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxxzzzzz_0[i] = g_xxxx_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxyyyyyy_0[i] = 6.0 * g_xxxx_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyyy_1[i] * wa_y[i];

        g_xxxxy_0_xxyyyyyz_0[i] = 5.0 * g_xxxx_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxxxy_0_xxyyyyzz_0[i] = 4.0 * g_xxxx_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxxxy_0_xxyyyzzz_0[i] = 3.0 * g_xxxx_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxyyzzzz_0[i] = 2.0 * g_xxxx_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxyzzzzz_0[i] = g_xxxx_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xxzzzzzz_0[i] = g_xxxx_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xyyyyyyy_0[i] = 7.0 * g_xxxx_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyyy_1[i] * wa_y[i];

        g_xxxxy_0_xyyyyyyz_0[i] = 6.0 * g_xxxx_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxxxy_0_xyyyyyzz_0[i] = 5.0 * g_xxxx_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxxxy_0_xyyyyzzz_0[i] = 4.0 * g_xxxx_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxxxy_0_xyyyzzzz_0[i] = 3.0 * g_xxxx_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xyyzzzzz_0[i] = 2.0 * g_xxxx_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xyzzzzzz_0[i] = g_xxxx_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_xzzzzzzz_0[i] = g_xxxx_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_yyyyyyyy_0[i] = 8.0 * g_xxxx_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxx_0_yyyyyyyy_1[i] * wa_y[i];

        g_xxxxy_0_yyyyyyyz_0[i] = 7.0 * g_xxxx_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxx_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxxxy_0_yyyyyyzz_0[i] = 6.0 * g_xxxx_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxx_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxxxy_0_yyyyyzzz_0[i] = 5.0 * g_xxxx_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxx_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxxxy_0_yyyyzzzz_0[i] = 4.0 * g_xxxx_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxx_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxxxy_0_yyyzzzzz_0[i] = 3.0 * g_xxxx_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxx_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_yyzzzzzz_0[i] = 2.0 * g_xxxx_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxx_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_yzzzzzzz_0[i] = g_xxxx_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxx_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxxxy_0_zzzzzzzz_0[i] = g_xxxx_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 90-135 components of targeted buffer : HSL

    auto g_xxxxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 90);

    auto g_xxxxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 91);

    auto g_xxxxz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 92);

    auto g_xxxxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 93);

    auto g_xxxxz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 94);

    auto g_xxxxz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 95);

    auto g_xxxxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 96);

    auto g_xxxxz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 97);

    auto g_xxxxz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 98);

    auto g_xxxxz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 99);

    auto g_xxxxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 100);

    auto g_xxxxz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 101);

    auto g_xxxxz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 102);

    auto g_xxxxz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 103);

    auto g_xxxxz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 104);

    auto g_xxxxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 105);

    auto g_xxxxz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 106);

    auto g_xxxxz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 107);

    auto g_xxxxz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 108);

    auto g_xxxxz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 109);

    auto g_xxxxz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 110);

    auto g_xxxxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 111);

    auto g_xxxxz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 112);

    auto g_xxxxz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 113);

    auto g_xxxxz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 114);

    auto g_xxxxz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 115);

    auto g_xxxxz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 116);

    auto g_xxxxz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 117);

    auto g_xxxxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 118);

    auto g_xxxxz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 119);

    auto g_xxxxz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 120);

    auto g_xxxxz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 121);

    auto g_xxxxz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 122);

    auto g_xxxxz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 123);

    auto g_xxxxz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 124);

    auto g_xxxxz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 125);

    auto g_xxxxz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 126);

    auto g_xxxxz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 127);

    auto g_xxxxz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 128);

    auto g_xxxxz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 129);

    auto g_xxxxz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 130);

    auto g_xxxxz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 131);

    auto g_xxxxz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 132);

    auto g_xxxxz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 133);

    auto g_xxxxz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 134);

    #pragma omp simd aligned(g_xxxx_0_xxxxxxx_1, g_xxxx_0_xxxxxxxx_1, g_xxxx_0_xxxxxxxy_1, g_xxxx_0_xxxxxxxz_1, g_xxxx_0_xxxxxxy_1, g_xxxx_0_xxxxxxyy_1, g_xxxx_0_xxxxxxyz_1, g_xxxx_0_xxxxxxz_1, g_xxxx_0_xxxxxxzz_1, g_xxxx_0_xxxxxyy_1, g_xxxx_0_xxxxxyyy_1, g_xxxx_0_xxxxxyyz_1, g_xxxx_0_xxxxxyz_1, g_xxxx_0_xxxxxyzz_1, g_xxxx_0_xxxxxzz_1, g_xxxx_0_xxxxxzzz_1, g_xxxx_0_xxxxyyy_1, g_xxxx_0_xxxxyyyy_1, g_xxxx_0_xxxxyyyz_1, g_xxxx_0_xxxxyyz_1, g_xxxx_0_xxxxyyzz_1, g_xxxx_0_xxxxyzz_1, g_xxxx_0_xxxxyzzz_1, g_xxxx_0_xxxxzzz_1, g_xxxx_0_xxxxzzzz_1, g_xxxx_0_xxxyyyy_1, g_xxxx_0_xxxyyyyy_1, g_xxxx_0_xxxyyyyz_1, g_xxxx_0_xxxyyyz_1, g_xxxx_0_xxxyyyzz_1, g_xxxx_0_xxxyyzz_1, g_xxxx_0_xxxyyzzz_1, g_xxxx_0_xxxyzzz_1, g_xxxx_0_xxxyzzzz_1, g_xxxx_0_xxxzzzz_1, g_xxxx_0_xxxzzzzz_1, g_xxxx_0_xxyyyyy_1, g_xxxx_0_xxyyyyyy_1, g_xxxx_0_xxyyyyyz_1, g_xxxx_0_xxyyyyz_1, g_xxxx_0_xxyyyyzz_1, g_xxxx_0_xxyyyzz_1, g_xxxx_0_xxyyyzzz_1, g_xxxx_0_xxyyzzz_1, g_xxxx_0_xxyyzzzz_1, g_xxxx_0_xxyzzzz_1, g_xxxx_0_xxyzzzzz_1, g_xxxx_0_xxzzzzz_1, g_xxxx_0_xxzzzzzz_1, g_xxxx_0_xyyyyyy_1, g_xxxx_0_xyyyyyyy_1, g_xxxx_0_xyyyyyyz_1, g_xxxx_0_xyyyyyz_1, g_xxxx_0_xyyyyyzz_1, g_xxxx_0_xyyyyzz_1, g_xxxx_0_xyyyyzzz_1, g_xxxx_0_xyyyzzz_1, g_xxxx_0_xyyyzzzz_1, g_xxxx_0_xyyzzzz_1, g_xxxx_0_xyyzzzzz_1, g_xxxx_0_xyzzzzz_1, g_xxxx_0_xyzzzzzz_1, g_xxxx_0_xzzzzzz_1, g_xxxx_0_xzzzzzzz_1, g_xxxx_0_yyyyyyy_1, g_xxxx_0_yyyyyyyy_1, g_xxxx_0_yyyyyyyz_1, g_xxxx_0_yyyyyyz_1, g_xxxx_0_yyyyyyzz_1, g_xxxx_0_yyyyyzz_1, g_xxxx_0_yyyyyzzz_1, g_xxxx_0_yyyyzzz_1, g_xxxx_0_yyyyzzzz_1, g_xxxx_0_yyyzzzz_1, g_xxxx_0_yyyzzzzz_1, g_xxxx_0_yyzzzzz_1, g_xxxx_0_yyzzzzzz_1, g_xxxx_0_yzzzzzz_1, g_xxxx_0_yzzzzzzz_1, g_xxxx_0_zzzzzzz_1, g_xxxx_0_zzzzzzzz_1, g_xxxxz_0_xxxxxxxx_0, g_xxxxz_0_xxxxxxxy_0, g_xxxxz_0_xxxxxxxz_0, g_xxxxz_0_xxxxxxyy_0, g_xxxxz_0_xxxxxxyz_0, g_xxxxz_0_xxxxxxzz_0, g_xxxxz_0_xxxxxyyy_0, g_xxxxz_0_xxxxxyyz_0, g_xxxxz_0_xxxxxyzz_0, g_xxxxz_0_xxxxxzzz_0, g_xxxxz_0_xxxxyyyy_0, g_xxxxz_0_xxxxyyyz_0, g_xxxxz_0_xxxxyyzz_0, g_xxxxz_0_xxxxyzzz_0, g_xxxxz_0_xxxxzzzz_0, g_xxxxz_0_xxxyyyyy_0, g_xxxxz_0_xxxyyyyz_0, g_xxxxz_0_xxxyyyzz_0, g_xxxxz_0_xxxyyzzz_0, g_xxxxz_0_xxxyzzzz_0, g_xxxxz_0_xxxzzzzz_0, g_xxxxz_0_xxyyyyyy_0, g_xxxxz_0_xxyyyyyz_0, g_xxxxz_0_xxyyyyzz_0, g_xxxxz_0_xxyyyzzz_0, g_xxxxz_0_xxyyzzzz_0, g_xxxxz_0_xxyzzzzz_0, g_xxxxz_0_xxzzzzzz_0, g_xxxxz_0_xyyyyyyy_0, g_xxxxz_0_xyyyyyyz_0, g_xxxxz_0_xyyyyyzz_0, g_xxxxz_0_xyyyyzzz_0, g_xxxxz_0_xyyyzzzz_0, g_xxxxz_0_xyyzzzzz_0, g_xxxxz_0_xyzzzzzz_0, g_xxxxz_0_xzzzzzzz_0, g_xxxxz_0_yyyyyyyy_0, g_xxxxz_0_yyyyyyyz_0, g_xxxxz_0_yyyyyyzz_0, g_xxxxz_0_yyyyyzzz_0, g_xxxxz_0_yyyyzzzz_0, g_xxxxz_0_yyyzzzzz_0, g_xxxxz_0_yyzzzzzz_0, g_xxxxz_0_yzzzzzzz_0, g_xxxxz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxz_0_xxxxxxxx_0[i] = g_xxxx_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxxxy_0[i] = g_xxxx_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxxxz_0[i] = g_xxxx_0_xxxxxxx_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxxz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxxyy_0[i] = g_xxxx_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxxyz_0[i] = g_xxxx_0_xxxxxxy_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxyz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxxzz_0[i] = 2.0 * g_xxxx_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxxzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxyyy_0[i] = g_xxxx_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxyyz_0[i] = g_xxxx_0_xxxxxyy_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyyz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxyzz_0[i] = 2.0 * g_xxxx_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxyzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxxzzz_0[i] = 3.0 * g_xxxx_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxxzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxyyyy_0[i] = g_xxxx_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxxz_0_xxxxyyyz_0[i] = g_xxxx_0_xxxxyyy_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyyz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxyyzz_0[i] = 2.0 * g_xxxx_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyyzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxyzzz_0[i] = 3.0 * g_xxxx_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxyzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxxzzzz_0[i] = 4.0 * g_xxxx_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxxzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxyyyyy_0[i] = g_xxxx_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxxz_0_xxxyyyyz_0[i] = g_xxxx_0_xxxyyyy_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyyz_1[i] * wa_z[i];

        g_xxxxz_0_xxxyyyzz_0[i] = 2.0 * g_xxxx_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyyzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxyyzzz_0[i] = 3.0 * g_xxxx_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyyzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxyzzzz_0[i] = 4.0 * g_xxxx_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxyzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxxzzzzz_0[i] = 5.0 * g_xxxx_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxxzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxyyyyyy_0[i] = g_xxxx_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxxz_0_xxyyyyyz_0[i] = g_xxxx_0_xxyyyyy_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyyz_1[i] * wa_z[i];

        g_xxxxz_0_xxyyyyzz_0[i] = 2.0 * g_xxxx_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyyzz_1[i] * wa_z[i];

        g_xxxxz_0_xxyyyzzz_0[i] = 3.0 * g_xxxx_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyyzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxyyzzzz_0[i] = 4.0 * g_xxxx_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyyzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxyzzzzz_0[i] = 5.0 * g_xxxx_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxyzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xxzzzzzz_0[i] = 6.0 * g_xxxx_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xxzzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xyyyyyyy_0[i] = g_xxxx_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxxz_0_xyyyyyyz_0[i] = g_xxxx_0_xyyyyyy_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyyz_1[i] * wa_z[i];

        g_xxxxz_0_xyyyyyzz_0[i] = 2.0 * g_xxxx_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyyzz_1[i] * wa_z[i];

        g_xxxxz_0_xyyyyzzz_0[i] = 3.0 * g_xxxx_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyyzzz_1[i] * wa_z[i];

        g_xxxxz_0_xyyyzzzz_0[i] = 4.0 * g_xxxx_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyyzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xyyzzzzz_0[i] = 5.0 * g_xxxx_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyyzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xyzzzzzz_0[i] = 6.0 * g_xxxx_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xyzzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_xzzzzzzz_0[i] = 7.0 * g_xxxx_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxx_0_xzzzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_yyyyyyyy_0[i] = g_xxxx_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxxxz_0_yyyyyyyz_0[i] = g_xxxx_0_yyyyyyy_1[i] * fi_acd_0 + g_xxxx_0_yyyyyyyz_1[i] * wa_z[i];

        g_xxxxz_0_yyyyyyzz_0[i] = 2.0 * g_xxxx_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxx_0_yyyyyyzz_1[i] * wa_z[i];

        g_xxxxz_0_yyyyyzzz_0[i] = 3.0 * g_xxxx_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxx_0_yyyyyzzz_1[i] * wa_z[i];

        g_xxxxz_0_yyyyzzzz_0[i] = 4.0 * g_xxxx_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxx_0_yyyyzzzz_1[i] * wa_z[i];

        g_xxxxz_0_yyyzzzzz_0[i] = 5.0 * g_xxxx_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxx_0_yyyzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_yyzzzzzz_0[i] = 6.0 * g_xxxx_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxx_0_yyzzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_yzzzzzzz_0[i] = 7.0 * g_xxxx_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxx_0_yzzzzzzz_1[i] * wa_z[i];

        g_xxxxz_0_zzzzzzzz_0[i] = 8.0 * g_xxxx_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxx_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 135-180 components of targeted buffer : HSL

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

    #pragma omp simd aligned(g_xxx_0_xxxxxxxx_0, g_xxx_0_xxxxxxxx_1, g_xxx_0_xxxxxxxz_0, g_xxx_0_xxxxxxxz_1, g_xxx_0_xxxxxxzz_0, g_xxx_0_xxxxxxzz_1, g_xxx_0_xxxxxzzz_0, g_xxx_0_xxxxxzzz_1, g_xxx_0_xxxxzzzz_0, g_xxx_0_xxxxzzzz_1, g_xxx_0_xxxzzzzz_0, g_xxx_0_xxxzzzzz_1, g_xxx_0_xxzzzzzz_0, g_xxx_0_xxzzzzzz_1, g_xxx_0_xzzzzzzz_0, g_xxx_0_xzzzzzzz_1, g_xxxy_0_xxxxxxxx_1, g_xxxy_0_xxxxxxxz_1, g_xxxy_0_xxxxxxzz_1, g_xxxy_0_xxxxxzzz_1, g_xxxy_0_xxxxzzzz_1, g_xxxy_0_xxxzzzzz_1, g_xxxy_0_xxzzzzzz_1, g_xxxy_0_xzzzzzzz_1, g_xxxyy_0_xxxxxxxx_0, g_xxxyy_0_xxxxxxxy_0, g_xxxyy_0_xxxxxxxz_0, g_xxxyy_0_xxxxxxyy_0, g_xxxyy_0_xxxxxxyz_0, g_xxxyy_0_xxxxxxzz_0, g_xxxyy_0_xxxxxyyy_0, g_xxxyy_0_xxxxxyyz_0, g_xxxyy_0_xxxxxyzz_0, g_xxxyy_0_xxxxxzzz_0, g_xxxyy_0_xxxxyyyy_0, g_xxxyy_0_xxxxyyyz_0, g_xxxyy_0_xxxxyyzz_0, g_xxxyy_0_xxxxyzzz_0, g_xxxyy_0_xxxxzzzz_0, g_xxxyy_0_xxxyyyyy_0, g_xxxyy_0_xxxyyyyz_0, g_xxxyy_0_xxxyyyzz_0, g_xxxyy_0_xxxyyzzz_0, g_xxxyy_0_xxxyzzzz_0, g_xxxyy_0_xxxzzzzz_0, g_xxxyy_0_xxyyyyyy_0, g_xxxyy_0_xxyyyyyz_0, g_xxxyy_0_xxyyyyzz_0, g_xxxyy_0_xxyyyzzz_0, g_xxxyy_0_xxyyzzzz_0, g_xxxyy_0_xxyzzzzz_0, g_xxxyy_0_xxzzzzzz_0, g_xxxyy_0_xyyyyyyy_0, g_xxxyy_0_xyyyyyyz_0, g_xxxyy_0_xyyyyyzz_0, g_xxxyy_0_xyyyyzzz_0, g_xxxyy_0_xyyyzzzz_0, g_xxxyy_0_xyyzzzzz_0, g_xxxyy_0_xyzzzzzz_0, g_xxxyy_0_xzzzzzzz_0, g_xxxyy_0_yyyyyyyy_0, g_xxxyy_0_yyyyyyyz_0, g_xxxyy_0_yyyyyyzz_0, g_xxxyy_0_yyyyyzzz_0, g_xxxyy_0_yyyyzzzz_0, g_xxxyy_0_yyyzzzzz_0, g_xxxyy_0_yyzzzzzz_0, g_xxxyy_0_yzzzzzzz_0, g_xxxyy_0_zzzzzzzz_0, g_xxyy_0_xxxxxxxy_1, g_xxyy_0_xxxxxxy_1, g_xxyy_0_xxxxxxyy_1, g_xxyy_0_xxxxxxyz_1, g_xxyy_0_xxxxxyy_1, g_xxyy_0_xxxxxyyy_1, g_xxyy_0_xxxxxyyz_1, g_xxyy_0_xxxxxyz_1, g_xxyy_0_xxxxxyzz_1, g_xxyy_0_xxxxyyy_1, g_xxyy_0_xxxxyyyy_1, g_xxyy_0_xxxxyyyz_1, g_xxyy_0_xxxxyyz_1, g_xxyy_0_xxxxyyzz_1, g_xxyy_0_xxxxyzz_1, g_xxyy_0_xxxxyzzz_1, g_xxyy_0_xxxyyyy_1, g_xxyy_0_xxxyyyyy_1, g_xxyy_0_xxxyyyyz_1, g_xxyy_0_xxxyyyz_1, g_xxyy_0_xxxyyyzz_1, g_xxyy_0_xxxyyzz_1, g_xxyy_0_xxxyyzzz_1, g_xxyy_0_xxxyzzz_1, g_xxyy_0_xxxyzzzz_1, g_xxyy_0_xxyyyyy_1, g_xxyy_0_xxyyyyyy_1, g_xxyy_0_xxyyyyyz_1, g_xxyy_0_xxyyyyz_1, g_xxyy_0_xxyyyyzz_1, g_xxyy_0_xxyyyzz_1, g_xxyy_0_xxyyyzzz_1, g_xxyy_0_xxyyzzz_1, g_xxyy_0_xxyyzzzz_1, g_xxyy_0_xxyzzzz_1, g_xxyy_0_xxyzzzzz_1, g_xxyy_0_xyyyyyy_1, g_xxyy_0_xyyyyyyy_1, g_xxyy_0_xyyyyyyz_1, g_xxyy_0_xyyyyyz_1, g_xxyy_0_xyyyyyzz_1, g_xxyy_0_xyyyyzz_1, g_xxyy_0_xyyyyzzz_1, g_xxyy_0_xyyyzzz_1, g_xxyy_0_xyyyzzzz_1, g_xxyy_0_xyyzzzz_1, g_xxyy_0_xyyzzzzz_1, g_xxyy_0_xyzzzzz_1, g_xxyy_0_xyzzzzzz_1, g_xxyy_0_yyyyyyy_1, g_xxyy_0_yyyyyyyy_1, g_xxyy_0_yyyyyyyz_1, g_xxyy_0_yyyyyyz_1, g_xxyy_0_yyyyyyzz_1, g_xxyy_0_yyyyyzz_1, g_xxyy_0_yyyyyzzz_1, g_xxyy_0_yyyyzzz_1, g_xxyy_0_yyyyzzzz_1, g_xxyy_0_yyyzzzz_1, g_xxyy_0_yyyzzzzz_1, g_xxyy_0_yyzzzzz_1, g_xxyy_0_yyzzzzzz_1, g_xxyy_0_yzzzzzz_1, g_xxyy_0_yzzzzzzz_1, g_xxyy_0_zzzzzzzz_1, g_xyy_0_xxxxxxxy_0, g_xyy_0_xxxxxxxy_1, g_xyy_0_xxxxxxyy_0, g_xyy_0_xxxxxxyy_1, g_xyy_0_xxxxxxyz_0, g_xyy_0_xxxxxxyz_1, g_xyy_0_xxxxxyyy_0, g_xyy_0_xxxxxyyy_1, g_xyy_0_xxxxxyyz_0, g_xyy_0_xxxxxyyz_1, g_xyy_0_xxxxxyzz_0, g_xyy_0_xxxxxyzz_1, g_xyy_0_xxxxyyyy_0, g_xyy_0_xxxxyyyy_1, g_xyy_0_xxxxyyyz_0, g_xyy_0_xxxxyyyz_1, g_xyy_0_xxxxyyzz_0, g_xyy_0_xxxxyyzz_1, g_xyy_0_xxxxyzzz_0, g_xyy_0_xxxxyzzz_1, g_xyy_0_xxxyyyyy_0, g_xyy_0_xxxyyyyy_1, g_xyy_0_xxxyyyyz_0, g_xyy_0_xxxyyyyz_1, g_xyy_0_xxxyyyzz_0, g_xyy_0_xxxyyyzz_1, g_xyy_0_xxxyyzzz_0, g_xyy_0_xxxyyzzz_1, g_xyy_0_xxxyzzzz_0, g_xyy_0_xxxyzzzz_1, g_xyy_0_xxyyyyyy_0, g_xyy_0_xxyyyyyy_1, g_xyy_0_xxyyyyyz_0, g_xyy_0_xxyyyyyz_1, g_xyy_0_xxyyyyzz_0, g_xyy_0_xxyyyyzz_1, g_xyy_0_xxyyyzzz_0, g_xyy_0_xxyyyzzz_1, g_xyy_0_xxyyzzzz_0, g_xyy_0_xxyyzzzz_1, g_xyy_0_xxyzzzzz_0, g_xyy_0_xxyzzzzz_1, g_xyy_0_xyyyyyyy_0, g_xyy_0_xyyyyyyy_1, g_xyy_0_xyyyyyyz_0, g_xyy_0_xyyyyyyz_1, g_xyy_0_xyyyyyzz_0, g_xyy_0_xyyyyyzz_1, g_xyy_0_xyyyyzzz_0, g_xyy_0_xyyyyzzz_1, g_xyy_0_xyyyzzzz_0, g_xyy_0_xyyyzzzz_1, g_xyy_0_xyyzzzzz_0, g_xyy_0_xyyzzzzz_1, g_xyy_0_xyzzzzzz_0, g_xyy_0_xyzzzzzz_1, g_xyy_0_yyyyyyyy_0, g_xyy_0_yyyyyyyy_1, g_xyy_0_yyyyyyyz_0, g_xyy_0_yyyyyyyz_1, g_xyy_0_yyyyyyzz_0, g_xyy_0_yyyyyyzz_1, g_xyy_0_yyyyyzzz_0, g_xyy_0_yyyyyzzz_1, g_xyy_0_yyyyzzzz_0, g_xyy_0_yyyyzzzz_1, g_xyy_0_yyyzzzzz_0, g_xyy_0_yyyzzzzz_1, g_xyy_0_yyzzzzzz_0, g_xyy_0_yyzzzzzz_1, g_xyy_0_yzzzzzzz_0, g_xyy_0_yzzzzzzz_1, g_xyy_0_zzzzzzzz_0, g_xyy_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyy_0_xxxxxxxx_0[i] = g_xxx_0_xxxxxxxx_0[i] * fbe_0 - g_xxx_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxy_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxyy_0_xxxxxxxy_0[i] = 2.0 * g_xyy_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xxyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xxyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxxyy_0_xxxxxxxz_0[i] = g_xxx_0_xxxxxxxz_0[i] * fbe_0 - g_xxx_0_xxxxxxxz_1[i] * fz_be_0 + g_xxxy_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxyy_0_xxxxxxyy_0[i] = 2.0 * g_xyy_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xxyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xxyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxxyy_0_xxxxxxyz_0[i] = 2.0 * g_xyy_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xxyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxyy_0_xxxxxxzz_0[i] = g_xxx_0_xxxxxxzz_0[i] * fbe_0 - g_xxx_0_xxxxxxzz_1[i] * fz_be_0 + g_xxxy_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxyy_0_xxxxxyyy_0[i] = 2.0 * g_xyy_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xxyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xxyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxxyy_0_xxxxxyyz_0[i] = 2.0 * g_xyy_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xxyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxyy_0_xxxxxyzz_0[i] = 2.0 * g_xyy_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xxyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxyy_0_xxxxxzzz_0[i] = g_xxx_0_xxxxxzzz_0[i] * fbe_0 - g_xxx_0_xxxxxzzz_1[i] * fz_be_0 + g_xxxy_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxyy_0_xxxxyyyy_0[i] = 2.0 * g_xyy_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xxyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xxyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxxyy_0_xxxxyyyz_0[i] = 2.0 * g_xyy_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xxyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxyy_0_xxxxyyzz_0[i] = 2.0 * g_xyy_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xxyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxyy_0_xxxxyzzz_0[i] = 2.0 * g_xyy_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xxyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxyy_0_xxxxzzzz_0[i] = g_xxx_0_xxxxzzzz_0[i] * fbe_0 - g_xxx_0_xxxxzzzz_1[i] * fz_be_0 + g_xxxy_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxyy_0_xxxyyyyy_0[i] = 2.0 * g_xyy_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xxyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxxyy_0_xxxyyyyz_0[i] = 2.0 * g_xyy_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xxyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxyy_0_xxxyyyzz_0[i] = 2.0 * g_xyy_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xxyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxyy_0_xxxyyzzz_0[i] = 2.0 * g_xyy_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xxyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxyy_0_xxxyzzzz_0[i] = 2.0 * g_xyy_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xxyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxyy_0_xxxzzzzz_0[i] = g_xxx_0_xxxzzzzz_0[i] * fbe_0 - g_xxx_0_xxxzzzzz_1[i] * fz_be_0 + g_xxxy_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxyy_0_xxyyyyyy_0[i] = 2.0 * g_xyy_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xxyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxxyy_0_xxyyyyyz_0[i] = 2.0 * g_xyy_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xxyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxyy_0_xxyyyyzz_0[i] = 2.0 * g_xyy_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xxyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxyy_0_xxyyyzzz_0[i] = 2.0 * g_xyy_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xxyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxyy_0_xxyyzzzz_0[i] = 2.0 * g_xyy_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xxyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxyy_0_xxyzzzzz_0[i] = 2.0 * g_xyy_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xxyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxyy_0_xxzzzzzz_0[i] = g_xxx_0_xxzzzzzz_0[i] * fbe_0 - g_xxx_0_xxzzzzzz_1[i] * fz_be_0 + g_xxxy_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxyy_0_xyyyyyyy_0[i] = 2.0 * g_xyy_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyyyyy_1[i] * fz_be_0 + g_xxyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xxyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxxyy_0_xyyyyyyz_0[i] = 2.0 * g_xyy_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyyyyz_1[i] * fz_be_0 + g_xxyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xxyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxyy_0_xyyyyyzz_0[i] = 2.0 * g_xyy_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyyyzz_1[i] * fz_be_0 + g_xxyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xxyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxyy_0_xyyyyzzz_0[i] = 2.0 * g_xyy_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyyzzz_1[i] * fz_be_0 + g_xxyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xxyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxyy_0_xyyyzzzz_0[i] = 2.0 * g_xyy_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyyzzzz_1[i] * fz_be_0 + g_xxyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xxyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxyy_0_xyyzzzzz_0[i] = 2.0 * g_xyy_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyyzzzzz_1[i] * fz_be_0 + g_xxyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xxyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxyy_0_xyzzzzzz_0[i] = 2.0 * g_xyy_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_xyzzzzzz_1[i] * fz_be_0 + g_xxyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xxyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxyy_0_xzzzzzzz_0[i] = g_xxx_0_xzzzzzzz_0[i] * fbe_0 - g_xxx_0_xzzzzzzz_1[i] * fz_be_0 + g_xxxy_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxyy_0_yyyyyyyy_0[i] = 2.0 * g_xyy_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyyyyy_1[i] * fz_be_0 + g_xxyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxyy_0_yyyyyyyz_0[i] = 2.0 * g_xyy_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyyyyz_1[i] * fz_be_0 + g_xxyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxyy_0_yyyyyyzz_0[i] = 2.0 * g_xyy_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyyyzz_1[i] * fz_be_0 + g_xxyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxyy_0_yyyyyzzz_0[i] = 2.0 * g_xyy_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyyzzz_1[i] * fz_be_0 + g_xxyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxyy_0_yyyyzzzz_0[i] = 2.0 * g_xyy_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyyzzzz_1[i] * fz_be_0 + g_xxyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxyy_0_yyyzzzzz_0[i] = 2.0 * g_xyy_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyyzzzzz_1[i] * fz_be_0 + g_xxyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxyy_0_yyzzzzzz_0[i] = 2.0 * g_xyy_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yyzzzzzz_1[i] * fz_be_0 + g_xxyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxyy_0_yzzzzzzz_0[i] = 2.0 * g_xyy_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_yzzzzzzz_1[i] * fz_be_0 + g_xxyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxyy_0_zzzzzzzz_0[i] = 2.0 * g_xyy_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_0_zzzzzzzz_1[i] * fz_be_0 + g_xxyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 180-225 components of targeted buffer : HSL

    auto g_xxxyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 180);

    auto g_xxxyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 181);

    auto g_xxxyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 182);

    auto g_xxxyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 183);

    auto g_xxxyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 184);

    auto g_xxxyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 185);

    auto g_xxxyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 186);

    auto g_xxxyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 187);

    auto g_xxxyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 188);

    auto g_xxxyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 189);

    auto g_xxxyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 190);

    auto g_xxxyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 191);

    auto g_xxxyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 192);

    auto g_xxxyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 193);

    auto g_xxxyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 194);

    auto g_xxxyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 195);

    auto g_xxxyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 196);

    auto g_xxxyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 197);

    auto g_xxxyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 198);

    auto g_xxxyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 199);

    auto g_xxxyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 200);

    auto g_xxxyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 201);

    auto g_xxxyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 202);

    auto g_xxxyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 203);

    auto g_xxxyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 204);

    auto g_xxxyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 205);

    auto g_xxxyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 206);

    auto g_xxxyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 207);

    auto g_xxxyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 208);

    auto g_xxxyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 209);

    auto g_xxxyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 210);

    auto g_xxxyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 211);

    auto g_xxxyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 212);

    auto g_xxxyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 213);

    auto g_xxxyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 214);

    auto g_xxxyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 215);

    auto g_xxxyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 216);

    auto g_xxxyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 217);

    auto g_xxxyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 218);

    auto g_xxxyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 219);

    auto g_xxxyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 220);

    auto g_xxxyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 221);

    auto g_xxxyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 222);

    auto g_xxxyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 223);

    auto g_xxxyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 224);

    #pragma omp simd aligned(g_xxxy_0_xxxxxxxy_1, g_xxxy_0_xxxxxxyy_1, g_xxxy_0_xxxxxyyy_1, g_xxxy_0_xxxxyyyy_1, g_xxxy_0_xxxyyyyy_1, g_xxxy_0_xxyyyyyy_1, g_xxxy_0_xyyyyyyy_1, g_xxxy_0_yyyyyyyy_1, g_xxxyz_0_xxxxxxxx_0, g_xxxyz_0_xxxxxxxy_0, g_xxxyz_0_xxxxxxxz_0, g_xxxyz_0_xxxxxxyy_0, g_xxxyz_0_xxxxxxyz_0, g_xxxyz_0_xxxxxxzz_0, g_xxxyz_0_xxxxxyyy_0, g_xxxyz_0_xxxxxyyz_0, g_xxxyz_0_xxxxxyzz_0, g_xxxyz_0_xxxxxzzz_0, g_xxxyz_0_xxxxyyyy_0, g_xxxyz_0_xxxxyyyz_0, g_xxxyz_0_xxxxyyzz_0, g_xxxyz_0_xxxxyzzz_0, g_xxxyz_0_xxxxzzzz_0, g_xxxyz_0_xxxyyyyy_0, g_xxxyz_0_xxxyyyyz_0, g_xxxyz_0_xxxyyyzz_0, g_xxxyz_0_xxxyyzzz_0, g_xxxyz_0_xxxyzzzz_0, g_xxxyz_0_xxxzzzzz_0, g_xxxyz_0_xxyyyyyy_0, g_xxxyz_0_xxyyyyyz_0, g_xxxyz_0_xxyyyyzz_0, g_xxxyz_0_xxyyyzzz_0, g_xxxyz_0_xxyyzzzz_0, g_xxxyz_0_xxyzzzzz_0, g_xxxyz_0_xxzzzzzz_0, g_xxxyz_0_xyyyyyyy_0, g_xxxyz_0_xyyyyyyz_0, g_xxxyz_0_xyyyyyzz_0, g_xxxyz_0_xyyyyzzz_0, g_xxxyz_0_xyyyzzzz_0, g_xxxyz_0_xyyzzzzz_0, g_xxxyz_0_xyzzzzzz_0, g_xxxyz_0_xzzzzzzz_0, g_xxxyz_0_yyyyyyyy_0, g_xxxyz_0_yyyyyyyz_0, g_xxxyz_0_yyyyyyzz_0, g_xxxyz_0_yyyyyzzz_0, g_xxxyz_0_yyyyzzzz_0, g_xxxyz_0_yyyzzzzz_0, g_xxxyz_0_yyzzzzzz_0, g_xxxyz_0_yzzzzzzz_0, g_xxxyz_0_zzzzzzzz_0, g_xxxz_0_xxxxxxxx_1, g_xxxz_0_xxxxxxxz_1, g_xxxz_0_xxxxxxyz_1, g_xxxz_0_xxxxxxz_1, g_xxxz_0_xxxxxxzz_1, g_xxxz_0_xxxxxyyz_1, g_xxxz_0_xxxxxyz_1, g_xxxz_0_xxxxxyzz_1, g_xxxz_0_xxxxxzz_1, g_xxxz_0_xxxxxzzz_1, g_xxxz_0_xxxxyyyz_1, g_xxxz_0_xxxxyyz_1, g_xxxz_0_xxxxyyzz_1, g_xxxz_0_xxxxyzz_1, g_xxxz_0_xxxxyzzz_1, g_xxxz_0_xxxxzzz_1, g_xxxz_0_xxxxzzzz_1, g_xxxz_0_xxxyyyyz_1, g_xxxz_0_xxxyyyz_1, g_xxxz_0_xxxyyyzz_1, g_xxxz_0_xxxyyzz_1, g_xxxz_0_xxxyyzzz_1, g_xxxz_0_xxxyzzz_1, g_xxxz_0_xxxyzzzz_1, g_xxxz_0_xxxzzzz_1, g_xxxz_0_xxxzzzzz_1, g_xxxz_0_xxyyyyyz_1, g_xxxz_0_xxyyyyz_1, g_xxxz_0_xxyyyyzz_1, g_xxxz_0_xxyyyzz_1, g_xxxz_0_xxyyyzzz_1, g_xxxz_0_xxyyzzz_1, g_xxxz_0_xxyyzzzz_1, g_xxxz_0_xxyzzzz_1, g_xxxz_0_xxyzzzzz_1, g_xxxz_0_xxzzzzz_1, g_xxxz_0_xxzzzzzz_1, g_xxxz_0_xyyyyyyz_1, g_xxxz_0_xyyyyyz_1, g_xxxz_0_xyyyyyzz_1, g_xxxz_0_xyyyyzz_1, g_xxxz_0_xyyyyzzz_1, g_xxxz_0_xyyyzzz_1, g_xxxz_0_xyyyzzzz_1, g_xxxz_0_xyyzzzz_1, g_xxxz_0_xyyzzzzz_1, g_xxxz_0_xyzzzzz_1, g_xxxz_0_xyzzzzzz_1, g_xxxz_0_xzzzzzz_1, g_xxxz_0_xzzzzzzz_1, g_xxxz_0_yyyyyyyz_1, g_xxxz_0_yyyyyyz_1, g_xxxz_0_yyyyyyzz_1, g_xxxz_0_yyyyyzz_1, g_xxxz_0_yyyyyzzz_1, g_xxxz_0_yyyyzzz_1, g_xxxz_0_yyyyzzzz_1, g_xxxz_0_yyyzzzz_1, g_xxxz_0_yyyzzzzz_1, g_xxxz_0_yyzzzzz_1, g_xxxz_0_yyzzzzzz_1, g_xxxz_0_yzzzzzz_1, g_xxxz_0_yzzzzzzz_1, g_xxxz_0_zzzzzzz_1, g_xxxz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyz_0_xxxxxxxx_0[i] = g_xxxz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxyz_0_xxxxxxxy_0[i] = g_xxxy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxyz_0_xxxxxxxz_0[i] = g_xxxz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxxxyy_0[i] = g_xxxy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxyz_0_xxxxxxyz_0[i] = g_xxxz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxxz_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxxxzz_0[i] = g_xxxz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxxyyy_0[i] = g_xxxy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxyz_0_xxxxxyyz_0[i] = 2.0 * g_xxxz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxxz_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxxyzz_0[i] = g_xxxz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxxz_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxxzzz_0[i] = g_xxxz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxyyyy_0[i] = g_xxxy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxyz_0_xxxxyyyz_0[i] = 3.0 * g_xxxz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxxz_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxyyzz_0[i] = 2.0 * g_xxxz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxxz_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxyzzz_0[i] = g_xxxz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxxz_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxxzzzz_0[i] = g_xxxz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxyyyyy_0[i] = g_xxxy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxyz_0_xxxyyyyz_0[i] = 4.0 * g_xxxz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxxz_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxxyz_0_xxxyyyzz_0[i] = 3.0 * g_xxxz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxxz_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxyyzzz_0[i] = 2.0 * g_xxxz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxxz_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxyzzzz_0[i] = g_xxxz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxxz_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxxzzzzz_0[i] = g_xxxz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxyyyyyy_0[i] = g_xxxy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxyz_0_xxyyyyyz_0[i] = 5.0 * g_xxxz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxxz_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxxyz_0_xxyyyyzz_0[i] = 4.0 * g_xxxz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxxz_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxxyz_0_xxyyyzzz_0[i] = 3.0 * g_xxxz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxxz_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxyyzzzz_0[i] = 2.0 * g_xxxz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxxz_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxyzzzzz_0[i] = g_xxxz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxxz_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xxzzzzzz_0[i] = g_xxxz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xyyyyyyy_0[i] = g_xxxy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxyz_0_xyyyyyyz_0[i] = 6.0 * g_xxxz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxxz_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxxyz_0_xyyyyyzz_0[i] = 5.0 * g_xxxz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxxz_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxxyz_0_xyyyyzzz_0[i] = 4.0 * g_xxxz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxxz_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxxyz_0_xyyyzzzz_0[i] = 3.0 * g_xxxz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxxz_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xyyzzzzz_0[i] = 2.0 * g_xxxz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxxz_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xyzzzzzz_0[i] = g_xxxz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxxz_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_xzzzzzzz_0[i] = g_xxxz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_yyyyyyyy_0[i] = g_xxxy_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxxyz_0_yyyyyyyz_0[i] = 7.0 * g_xxxz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxxz_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxxyz_0_yyyyyyzz_0[i] = 6.0 * g_xxxz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxxz_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxxyz_0_yyyyyzzz_0[i] = 5.0 * g_xxxz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxxz_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxxyz_0_yyyyzzzz_0[i] = 4.0 * g_xxxz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxxz_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxxyz_0_yyyzzzzz_0[i] = 3.0 * g_xxxz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxxz_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_yyzzzzzz_0[i] = 2.0 * g_xxxz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxxz_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_yzzzzzzz_0[i] = g_xxxz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxxz_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxxyz_0_zzzzzzzz_0[i] = g_xxxz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 225-270 components of targeted buffer : HSL

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

    #pragma omp simd aligned(g_xxx_0_xxxxxxxx_0, g_xxx_0_xxxxxxxx_1, g_xxx_0_xxxxxxxy_0, g_xxx_0_xxxxxxxy_1, g_xxx_0_xxxxxxyy_0, g_xxx_0_xxxxxxyy_1, g_xxx_0_xxxxxyyy_0, g_xxx_0_xxxxxyyy_1, g_xxx_0_xxxxyyyy_0, g_xxx_0_xxxxyyyy_1, g_xxx_0_xxxyyyyy_0, g_xxx_0_xxxyyyyy_1, g_xxx_0_xxyyyyyy_0, g_xxx_0_xxyyyyyy_1, g_xxx_0_xyyyyyyy_0, g_xxx_0_xyyyyyyy_1, g_xxxz_0_xxxxxxxx_1, g_xxxz_0_xxxxxxxy_1, g_xxxz_0_xxxxxxyy_1, g_xxxz_0_xxxxxyyy_1, g_xxxz_0_xxxxyyyy_1, g_xxxz_0_xxxyyyyy_1, g_xxxz_0_xxyyyyyy_1, g_xxxz_0_xyyyyyyy_1, g_xxxzz_0_xxxxxxxx_0, g_xxxzz_0_xxxxxxxy_0, g_xxxzz_0_xxxxxxxz_0, g_xxxzz_0_xxxxxxyy_0, g_xxxzz_0_xxxxxxyz_0, g_xxxzz_0_xxxxxxzz_0, g_xxxzz_0_xxxxxyyy_0, g_xxxzz_0_xxxxxyyz_0, g_xxxzz_0_xxxxxyzz_0, g_xxxzz_0_xxxxxzzz_0, g_xxxzz_0_xxxxyyyy_0, g_xxxzz_0_xxxxyyyz_0, g_xxxzz_0_xxxxyyzz_0, g_xxxzz_0_xxxxyzzz_0, g_xxxzz_0_xxxxzzzz_0, g_xxxzz_0_xxxyyyyy_0, g_xxxzz_0_xxxyyyyz_0, g_xxxzz_0_xxxyyyzz_0, g_xxxzz_0_xxxyyzzz_0, g_xxxzz_0_xxxyzzzz_0, g_xxxzz_0_xxxzzzzz_0, g_xxxzz_0_xxyyyyyy_0, g_xxxzz_0_xxyyyyyz_0, g_xxxzz_0_xxyyyyzz_0, g_xxxzz_0_xxyyyzzz_0, g_xxxzz_0_xxyyzzzz_0, g_xxxzz_0_xxyzzzzz_0, g_xxxzz_0_xxzzzzzz_0, g_xxxzz_0_xyyyyyyy_0, g_xxxzz_0_xyyyyyyz_0, g_xxxzz_0_xyyyyyzz_0, g_xxxzz_0_xyyyyzzz_0, g_xxxzz_0_xyyyzzzz_0, g_xxxzz_0_xyyzzzzz_0, g_xxxzz_0_xyzzzzzz_0, g_xxxzz_0_xzzzzzzz_0, g_xxxzz_0_yyyyyyyy_0, g_xxxzz_0_yyyyyyyz_0, g_xxxzz_0_yyyyyyzz_0, g_xxxzz_0_yyyyyzzz_0, g_xxxzz_0_yyyyzzzz_0, g_xxxzz_0_yyyzzzzz_0, g_xxxzz_0_yyzzzzzz_0, g_xxxzz_0_yzzzzzzz_0, g_xxxzz_0_zzzzzzzz_0, g_xxzz_0_xxxxxxxz_1, g_xxzz_0_xxxxxxyz_1, g_xxzz_0_xxxxxxz_1, g_xxzz_0_xxxxxxzz_1, g_xxzz_0_xxxxxyyz_1, g_xxzz_0_xxxxxyz_1, g_xxzz_0_xxxxxyzz_1, g_xxzz_0_xxxxxzz_1, g_xxzz_0_xxxxxzzz_1, g_xxzz_0_xxxxyyyz_1, g_xxzz_0_xxxxyyz_1, g_xxzz_0_xxxxyyzz_1, g_xxzz_0_xxxxyzz_1, g_xxzz_0_xxxxyzzz_1, g_xxzz_0_xxxxzzz_1, g_xxzz_0_xxxxzzzz_1, g_xxzz_0_xxxyyyyz_1, g_xxzz_0_xxxyyyz_1, g_xxzz_0_xxxyyyzz_1, g_xxzz_0_xxxyyzz_1, g_xxzz_0_xxxyyzzz_1, g_xxzz_0_xxxyzzz_1, g_xxzz_0_xxxyzzzz_1, g_xxzz_0_xxxzzzz_1, g_xxzz_0_xxxzzzzz_1, g_xxzz_0_xxyyyyyz_1, g_xxzz_0_xxyyyyz_1, g_xxzz_0_xxyyyyzz_1, g_xxzz_0_xxyyyzz_1, g_xxzz_0_xxyyyzzz_1, g_xxzz_0_xxyyzzz_1, g_xxzz_0_xxyyzzzz_1, g_xxzz_0_xxyzzzz_1, g_xxzz_0_xxyzzzzz_1, g_xxzz_0_xxzzzzz_1, g_xxzz_0_xxzzzzzz_1, g_xxzz_0_xyyyyyyz_1, g_xxzz_0_xyyyyyz_1, g_xxzz_0_xyyyyyzz_1, g_xxzz_0_xyyyyzz_1, g_xxzz_0_xyyyyzzz_1, g_xxzz_0_xyyyzzz_1, g_xxzz_0_xyyyzzzz_1, g_xxzz_0_xyyzzzz_1, g_xxzz_0_xyyzzzzz_1, g_xxzz_0_xyzzzzz_1, g_xxzz_0_xyzzzzzz_1, g_xxzz_0_xzzzzzz_1, g_xxzz_0_xzzzzzzz_1, g_xxzz_0_yyyyyyyy_1, g_xxzz_0_yyyyyyyz_1, g_xxzz_0_yyyyyyz_1, g_xxzz_0_yyyyyyzz_1, g_xxzz_0_yyyyyzz_1, g_xxzz_0_yyyyyzzz_1, g_xxzz_0_yyyyzzz_1, g_xxzz_0_yyyyzzzz_1, g_xxzz_0_yyyzzzz_1, g_xxzz_0_yyyzzzzz_1, g_xxzz_0_yyzzzzz_1, g_xxzz_0_yyzzzzzz_1, g_xxzz_0_yzzzzzz_1, g_xxzz_0_yzzzzzzz_1, g_xxzz_0_zzzzzzz_1, g_xxzz_0_zzzzzzzz_1, g_xzz_0_xxxxxxxz_0, g_xzz_0_xxxxxxxz_1, g_xzz_0_xxxxxxyz_0, g_xzz_0_xxxxxxyz_1, g_xzz_0_xxxxxxzz_0, g_xzz_0_xxxxxxzz_1, g_xzz_0_xxxxxyyz_0, g_xzz_0_xxxxxyyz_1, g_xzz_0_xxxxxyzz_0, g_xzz_0_xxxxxyzz_1, g_xzz_0_xxxxxzzz_0, g_xzz_0_xxxxxzzz_1, g_xzz_0_xxxxyyyz_0, g_xzz_0_xxxxyyyz_1, g_xzz_0_xxxxyyzz_0, g_xzz_0_xxxxyyzz_1, g_xzz_0_xxxxyzzz_0, g_xzz_0_xxxxyzzz_1, g_xzz_0_xxxxzzzz_0, g_xzz_0_xxxxzzzz_1, g_xzz_0_xxxyyyyz_0, g_xzz_0_xxxyyyyz_1, g_xzz_0_xxxyyyzz_0, g_xzz_0_xxxyyyzz_1, g_xzz_0_xxxyyzzz_0, g_xzz_0_xxxyyzzz_1, g_xzz_0_xxxyzzzz_0, g_xzz_0_xxxyzzzz_1, g_xzz_0_xxxzzzzz_0, g_xzz_0_xxxzzzzz_1, g_xzz_0_xxyyyyyz_0, g_xzz_0_xxyyyyyz_1, g_xzz_0_xxyyyyzz_0, g_xzz_0_xxyyyyzz_1, g_xzz_0_xxyyyzzz_0, g_xzz_0_xxyyyzzz_1, g_xzz_0_xxyyzzzz_0, g_xzz_0_xxyyzzzz_1, g_xzz_0_xxyzzzzz_0, g_xzz_0_xxyzzzzz_1, g_xzz_0_xxzzzzzz_0, g_xzz_0_xxzzzzzz_1, g_xzz_0_xyyyyyyz_0, g_xzz_0_xyyyyyyz_1, g_xzz_0_xyyyyyzz_0, g_xzz_0_xyyyyyzz_1, g_xzz_0_xyyyyzzz_0, g_xzz_0_xyyyyzzz_1, g_xzz_0_xyyyzzzz_0, g_xzz_0_xyyyzzzz_1, g_xzz_0_xyyzzzzz_0, g_xzz_0_xyyzzzzz_1, g_xzz_0_xyzzzzzz_0, g_xzz_0_xyzzzzzz_1, g_xzz_0_xzzzzzzz_0, g_xzz_0_xzzzzzzz_1, g_xzz_0_yyyyyyyy_0, g_xzz_0_yyyyyyyy_1, g_xzz_0_yyyyyyyz_0, g_xzz_0_yyyyyyyz_1, g_xzz_0_yyyyyyzz_0, g_xzz_0_yyyyyyzz_1, g_xzz_0_yyyyyzzz_0, g_xzz_0_yyyyyzzz_1, g_xzz_0_yyyyzzzz_0, g_xzz_0_yyyyzzzz_1, g_xzz_0_yyyzzzzz_0, g_xzz_0_yyyzzzzz_1, g_xzz_0_yyzzzzzz_0, g_xzz_0_yyzzzzzz_1, g_xzz_0_yzzzzzzz_0, g_xzz_0_yzzzzzzz_1, g_xzz_0_zzzzzzzz_0, g_xzz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzz_0_xxxxxxxx_0[i] = g_xxx_0_xxxxxxxx_0[i] * fbe_0 - g_xxx_0_xxxxxxxx_1[i] * fz_be_0 + g_xxxz_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxzz_0_xxxxxxxy_0[i] = g_xxx_0_xxxxxxxy_0[i] * fbe_0 - g_xxx_0_xxxxxxxy_1[i] * fz_be_0 + g_xxxz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxzz_0_xxxxxxxz_0[i] = 2.0 * g_xzz_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xxzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxxxyy_0[i] = g_xxx_0_xxxxxxyy_0[i] * fbe_0 - g_xxx_0_xxxxxxyy_1[i] * fz_be_0 + g_xxxz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxzz_0_xxxxxxyz_0[i] = 2.0 * g_xzz_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxxxzz_0[i] = 2.0 * g_xzz_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xxzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxxyyy_0[i] = g_xxx_0_xxxxxyyy_0[i] * fbe_0 - g_xxx_0_xxxxxyyy_1[i] * fz_be_0 + g_xxxz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxzz_0_xxxxxyyz_0[i] = 2.0 * g_xzz_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxxyzz_0[i] = 2.0 * g_xzz_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxxzzz_0[i] = 2.0 * g_xzz_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xxzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxyyyy_0[i] = g_xxx_0_xxxxyyyy_0[i] * fbe_0 - g_xxx_0_xxxxyyyy_1[i] * fz_be_0 + g_xxxz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxzz_0_xxxxyyyz_0[i] = 2.0 * g_xzz_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxyyzz_0[i] = 2.0 * g_xzz_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxyzzz_0[i] = 2.0 * g_xzz_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxxzzzz_0[i] = 2.0 * g_xzz_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xxzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxyyyyy_0[i] = g_xxx_0_xxxyyyyy_0[i] * fbe_0 - g_xxx_0_xxxyyyyy_1[i] * fz_be_0 + g_xxxz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxzz_0_xxxyyyyz_0[i] = 2.0 * g_xzz_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxzz_0_xxxyyyzz_0[i] = 2.0 * g_xzz_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxyyzzz_0[i] = 2.0 * g_xzz_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxyzzzz_0[i] = 2.0 * g_xzz_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxxzzzzz_0[i] = 2.0 * g_xzz_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xxzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxyyyyyy_0[i] = g_xxx_0_xxyyyyyy_0[i] * fbe_0 - g_xxx_0_xxyyyyyy_1[i] * fz_be_0 + g_xxxz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxzz_0_xxyyyyyz_0[i] = 2.0 * g_xzz_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxzz_0_xxyyyyzz_0[i] = 2.0 * g_xzz_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxzz_0_xxyyyzzz_0[i] = 2.0 * g_xzz_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxyyzzzz_0[i] = 2.0 * g_xzz_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxyzzzzz_0[i] = 2.0 * g_xzz_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xxzzzzzz_0[i] = 2.0 * g_xzz_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xyyyyyyy_0[i] = g_xxx_0_xyyyyyyy_0[i] * fbe_0 - g_xxx_0_xyyyyyyy_1[i] * fz_be_0 + g_xxxz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxzz_0_xyyyyyyz_0[i] = 2.0 * g_xzz_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xxzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxzz_0_xyyyyyzz_0[i] = 2.0 * g_xzz_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xxzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxzz_0_xyyyyzzz_0[i] = 2.0 * g_xzz_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xxzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxzz_0_xyyyzzzz_0[i] = 2.0 * g_xzz_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xxzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xyyzzzzz_0[i] = 2.0 * g_xzz_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xxzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xyzzzzzz_0[i] = 2.0 * g_xzz_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xxzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_xzzzzzzz_0[i] = 2.0 * g_xzz_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xxzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_yyyyyyyy_0[i] = 2.0 * g_xzz_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xxzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxzz_0_yyyyyyyz_0[i] = 2.0 * g_xzz_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xxzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxzz_0_yyyyyyzz_0[i] = 2.0 * g_xzz_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xxzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxzz_0_yyyyyzzz_0[i] = 2.0 * g_xzz_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xxzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxzz_0_yyyyzzzz_0[i] = 2.0 * g_xzz_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xxzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxzz_0_yyyzzzzz_0[i] = 2.0 * g_xzz_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xxzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_yyzzzzzz_0[i] = 2.0 * g_xzz_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xxzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_yzzzzzzz_0[i] = 2.0 * g_xzz_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xxzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxzz_0_zzzzzzzz_0[i] = 2.0 * g_xzz_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xxzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 270-315 components of targeted buffer : HSL

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

    #pragma omp simd aligned(g_xxy_0_xxxxxxxx_0, g_xxy_0_xxxxxxxx_1, g_xxy_0_xxxxxxxz_0, g_xxy_0_xxxxxxxz_1, g_xxy_0_xxxxxxzz_0, g_xxy_0_xxxxxxzz_1, g_xxy_0_xxxxxzzz_0, g_xxy_0_xxxxxzzz_1, g_xxy_0_xxxxzzzz_0, g_xxy_0_xxxxzzzz_1, g_xxy_0_xxxzzzzz_0, g_xxy_0_xxxzzzzz_1, g_xxy_0_xxzzzzzz_0, g_xxy_0_xxzzzzzz_1, g_xxy_0_xzzzzzzz_0, g_xxy_0_xzzzzzzz_1, g_xxyy_0_xxxxxxxx_1, g_xxyy_0_xxxxxxxz_1, g_xxyy_0_xxxxxxzz_1, g_xxyy_0_xxxxxzzz_1, g_xxyy_0_xxxxzzzz_1, g_xxyy_0_xxxzzzzz_1, g_xxyy_0_xxzzzzzz_1, g_xxyy_0_xzzzzzzz_1, g_xxyyy_0_xxxxxxxx_0, g_xxyyy_0_xxxxxxxy_0, g_xxyyy_0_xxxxxxxz_0, g_xxyyy_0_xxxxxxyy_0, g_xxyyy_0_xxxxxxyz_0, g_xxyyy_0_xxxxxxzz_0, g_xxyyy_0_xxxxxyyy_0, g_xxyyy_0_xxxxxyyz_0, g_xxyyy_0_xxxxxyzz_0, g_xxyyy_0_xxxxxzzz_0, g_xxyyy_0_xxxxyyyy_0, g_xxyyy_0_xxxxyyyz_0, g_xxyyy_0_xxxxyyzz_0, g_xxyyy_0_xxxxyzzz_0, g_xxyyy_0_xxxxzzzz_0, g_xxyyy_0_xxxyyyyy_0, g_xxyyy_0_xxxyyyyz_0, g_xxyyy_0_xxxyyyzz_0, g_xxyyy_0_xxxyyzzz_0, g_xxyyy_0_xxxyzzzz_0, g_xxyyy_0_xxxzzzzz_0, g_xxyyy_0_xxyyyyyy_0, g_xxyyy_0_xxyyyyyz_0, g_xxyyy_0_xxyyyyzz_0, g_xxyyy_0_xxyyyzzz_0, g_xxyyy_0_xxyyzzzz_0, g_xxyyy_0_xxyzzzzz_0, g_xxyyy_0_xxzzzzzz_0, g_xxyyy_0_xyyyyyyy_0, g_xxyyy_0_xyyyyyyz_0, g_xxyyy_0_xyyyyyzz_0, g_xxyyy_0_xyyyyzzz_0, g_xxyyy_0_xyyyzzzz_0, g_xxyyy_0_xyyzzzzz_0, g_xxyyy_0_xyzzzzzz_0, g_xxyyy_0_xzzzzzzz_0, g_xxyyy_0_yyyyyyyy_0, g_xxyyy_0_yyyyyyyz_0, g_xxyyy_0_yyyyyyzz_0, g_xxyyy_0_yyyyyzzz_0, g_xxyyy_0_yyyyzzzz_0, g_xxyyy_0_yyyzzzzz_0, g_xxyyy_0_yyzzzzzz_0, g_xxyyy_0_yzzzzzzz_0, g_xxyyy_0_zzzzzzzz_0, g_xyyy_0_xxxxxxxy_1, g_xyyy_0_xxxxxxy_1, g_xyyy_0_xxxxxxyy_1, g_xyyy_0_xxxxxxyz_1, g_xyyy_0_xxxxxyy_1, g_xyyy_0_xxxxxyyy_1, g_xyyy_0_xxxxxyyz_1, g_xyyy_0_xxxxxyz_1, g_xyyy_0_xxxxxyzz_1, g_xyyy_0_xxxxyyy_1, g_xyyy_0_xxxxyyyy_1, g_xyyy_0_xxxxyyyz_1, g_xyyy_0_xxxxyyz_1, g_xyyy_0_xxxxyyzz_1, g_xyyy_0_xxxxyzz_1, g_xyyy_0_xxxxyzzz_1, g_xyyy_0_xxxyyyy_1, g_xyyy_0_xxxyyyyy_1, g_xyyy_0_xxxyyyyz_1, g_xyyy_0_xxxyyyz_1, g_xyyy_0_xxxyyyzz_1, g_xyyy_0_xxxyyzz_1, g_xyyy_0_xxxyyzzz_1, g_xyyy_0_xxxyzzz_1, g_xyyy_0_xxxyzzzz_1, g_xyyy_0_xxyyyyy_1, g_xyyy_0_xxyyyyyy_1, g_xyyy_0_xxyyyyyz_1, g_xyyy_0_xxyyyyz_1, g_xyyy_0_xxyyyyzz_1, g_xyyy_0_xxyyyzz_1, g_xyyy_0_xxyyyzzz_1, g_xyyy_0_xxyyzzz_1, g_xyyy_0_xxyyzzzz_1, g_xyyy_0_xxyzzzz_1, g_xyyy_0_xxyzzzzz_1, g_xyyy_0_xyyyyyy_1, g_xyyy_0_xyyyyyyy_1, g_xyyy_0_xyyyyyyz_1, g_xyyy_0_xyyyyyz_1, g_xyyy_0_xyyyyyzz_1, g_xyyy_0_xyyyyzz_1, g_xyyy_0_xyyyyzzz_1, g_xyyy_0_xyyyzzz_1, g_xyyy_0_xyyyzzzz_1, g_xyyy_0_xyyzzzz_1, g_xyyy_0_xyyzzzzz_1, g_xyyy_0_xyzzzzz_1, g_xyyy_0_xyzzzzzz_1, g_xyyy_0_yyyyyyy_1, g_xyyy_0_yyyyyyyy_1, g_xyyy_0_yyyyyyyz_1, g_xyyy_0_yyyyyyz_1, g_xyyy_0_yyyyyyzz_1, g_xyyy_0_yyyyyzz_1, g_xyyy_0_yyyyyzzz_1, g_xyyy_0_yyyyzzz_1, g_xyyy_0_yyyyzzzz_1, g_xyyy_0_yyyzzzz_1, g_xyyy_0_yyyzzzzz_1, g_xyyy_0_yyzzzzz_1, g_xyyy_0_yyzzzzzz_1, g_xyyy_0_yzzzzzz_1, g_xyyy_0_yzzzzzzz_1, g_xyyy_0_zzzzzzzz_1, g_yyy_0_xxxxxxxy_0, g_yyy_0_xxxxxxxy_1, g_yyy_0_xxxxxxyy_0, g_yyy_0_xxxxxxyy_1, g_yyy_0_xxxxxxyz_0, g_yyy_0_xxxxxxyz_1, g_yyy_0_xxxxxyyy_0, g_yyy_0_xxxxxyyy_1, g_yyy_0_xxxxxyyz_0, g_yyy_0_xxxxxyyz_1, g_yyy_0_xxxxxyzz_0, g_yyy_0_xxxxxyzz_1, g_yyy_0_xxxxyyyy_0, g_yyy_0_xxxxyyyy_1, g_yyy_0_xxxxyyyz_0, g_yyy_0_xxxxyyyz_1, g_yyy_0_xxxxyyzz_0, g_yyy_0_xxxxyyzz_1, g_yyy_0_xxxxyzzz_0, g_yyy_0_xxxxyzzz_1, g_yyy_0_xxxyyyyy_0, g_yyy_0_xxxyyyyy_1, g_yyy_0_xxxyyyyz_0, g_yyy_0_xxxyyyyz_1, g_yyy_0_xxxyyyzz_0, g_yyy_0_xxxyyyzz_1, g_yyy_0_xxxyyzzz_0, g_yyy_0_xxxyyzzz_1, g_yyy_0_xxxyzzzz_0, g_yyy_0_xxxyzzzz_1, g_yyy_0_xxyyyyyy_0, g_yyy_0_xxyyyyyy_1, g_yyy_0_xxyyyyyz_0, g_yyy_0_xxyyyyyz_1, g_yyy_0_xxyyyyzz_0, g_yyy_0_xxyyyyzz_1, g_yyy_0_xxyyyzzz_0, g_yyy_0_xxyyyzzz_1, g_yyy_0_xxyyzzzz_0, g_yyy_0_xxyyzzzz_1, g_yyy_0_xxyzzzzz_0, g_yyy_0_xxyzzzzz_1, g_yyy_0_xyyyyyyy_0, g_yyy_0_xyyyyyyy_1, g_yyy_0_xyyyyyyz_0, g_yyy_0_xyyyyyyz_1, g_yyy_0_xyyyyyzz_0, g_yyy_0_xyyyyyzz_1, g_yyy_0_xyyyyzzz_0, g_yyy_0_xyyyyzzz_1, g_yyy_0_xyyyzzzz_0, g_yyy_0_xyyyzzzz_1, g_yyy_0_xyyzzzzz_0, g_yyy_0_xyyzzzzz_1, g_yyy_0_xyzzzzzz_0, g_yyy_0_xyzzzzzz_1, g_yyy_0_yyyyyyyy_0, g_yyy_0_yyyyyyyy_1, g_yyy_0_yyyyyyyz_0, g_yyy_0_yyyyyyyz_1, g_yyy_0_yyyyyyzz_0, g_yyy_0_yyyyyyzz_1, g_yyy_0_yyyyyzzz_0, g_yyy_0_yyyyyzzz_1, g_yyy_0_yyyyzzzz_0, g_yyy_0_yyyyzzzz_1, g_yyy_0_yyyzzzzz_0, g_yyy_0_yyyzzzzz_1, g_yyy_0_yyzzzzzz_0, g_yyy_0_yyzzzzzz_1, g_yyy_0_yzzzzzzz_0, g_yyy_0_yzzzzzzz_1, g_yyy_0_zzzzzzzz_0, g_yyy_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyy_0_xxxxxxxx_0[i] = 2.0 * g_xxy_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxxxxx_1[i] * fz_be_0 + g_xxyy_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxyyy_0_xxxxxxxy_0[i] = g_yyy_0_xxxxxxxy_0[i] * fbe_0 - g_yyy_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xyyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxyyy_0_xxxxxxxz_0[i] = 2.0 * g_xxy_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxxxxz_1[i] * fz_be_0 + g_xxyy_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxyyy_0_xxxxxxyy_0[i] = g_yyy_0_xxxxxxyy_0[i] * fbe_0 - g_yyy_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xyyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxyyy_0_xxxxxxyz_0[i] = g_yyy_0_xxxxxxyz_0[i] * fbe_0 - g_yyy_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xyyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxyyy_0_xxxxxxzz_0[i] = 2.0 * g_xxy_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxxxzz_1[i] * fz_be_0 + g_xxyy_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxyyy_0_xxxxxyyy_0[i] = g_yyy_0_xxxxxyyy_0[i] * fbe_0 - g_yyy_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xyyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxyyy_0_xxxxxyyz_0[i] = g_yyy_0_xxxxxyyz_0[i] * fbe_0 - g_yyy_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xyyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxyyy_0_xxxxxyzz_0[i] = g_yyy_0_xxxxxyzz_0[i] * fbe_0 - g_yyy_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xyyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxyyy_0_xxxxxzzz_0[i] = 2.0 * g_xxy_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxxzzz_1[i] * fz_be_0 + g_xxyy_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxyyy_0_xxxxyyyy_0[i] = g_yyy_0_xxxxyyyy_0[i] * fbe_0 - g_yyy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xyyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxyyy_0_xxxxyyyz_0[i] = g_yyy_0_xxxxyyyz_0[i] * fbe_0 - g_yyy_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xyyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxyyy_0_xxxxyyzz_0[i] = g_yyy_0_xxxxyyzz_0[i] * fbe_0 - g_yyy_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xyyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxyyy_0_xxxxyzzz_0[i] = g_yyy_0_xxxxyzzz_0[i] * fbe_0 - g_yyy_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xyyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxyyy_0_xxxxzzzz_0[i] = 2.0 * g_xxy_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxxzzzz_1[i] * fz_be_0 + g_xxyy_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxyyy_0_xxxyyyyy_0[i] = g_yyy_0_xxxyyyyy_0[i] * fbe_0 - g_yyy_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xyyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxyyy_0_xxxyyyyz_0[i] = g_yyy_0_xxxyyyyz_0[i] * fbe_0 - g_yyy_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xyyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxyyy_0_xxxyyyzz_0[i] = g_yyy_0_xxxyyyzz_0[i] * fbe_0 - g_yyy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xyyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxyyy_0_xxxyyzzz_0[i] = g_yyy_0_xxxyyzzz_0[i] * fbe_0 - g_yyy_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xyyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxyyy_0_xxxyzzzz_0[i] = g_yyy_0_xxxyzzzz_0[i] * fbe_0 - g_yyy_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xyyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxyyy_0_xxxzzzzz_0[i] = 2.0 * g_xxy_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxxzzzzz_1[i] * fz_be_0 + g_xxyy_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxyyy_0_xxyyyyyy_0[i] = g_yyy_0_xxyyyyyy_0[i] * fbe_0 - g_yyy_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xyyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxyyy_0_xxyyyyyz_0[i] = g_yyy_0_xxyyyyyz_0[i] * fbe_0 - g_yyy_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xyyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxyyy_0_xxyyyyzz_0[i] = g_yyy_0_xxyyyyzz_0[i] * fbe_0 - g_yyy_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xyyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxyyy_0_xxyyyzzz_0[i] = g_yyy_0_xxyyyzzz_0[i] * fbe_0 - g_yyy_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xyyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxyyy_0_xxyyzzzz_0[i] = g_yyy_0_xxyyzzzz_0[i] * fbe_0 - g_yyy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xyyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxyyy_0_xxyzzzzz_0[i] = g_yyy_0_xxyzzzzz_0[i] * fbe_0 - g_yyy_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xyyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxyyy_0_xxzzzzzz_0[i] = 2.0 * g_xxy_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xxzzzzzz_1[i] * fz_be_0 + g_xxyy_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxyyy_0_xyyyyyyy_0[i] = g_yyy_0_xyyyyyyy_0[i] * fbe_0 - g_yyy_0_xyyyyyyy_1[i] * fz_be_0 + g_xyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xyyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxyyy_0_xyyyyyyz_0[i] = g_yyy_0_xyyyyyyz_0[i] * fbe_0 - g_yyy_0_xyyyyyyz_1[i] * fz_be_0 + g_xyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xyyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxyyy_0_xyyyyyzz_0[i] = g_yyy_0_xyyyyyzz_0[i] * fbe_0 - g_yyy_0_xyyyyyzz_1[i] * fz_be_0 + g_xyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xyyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxyyy_0_xyyyyzzz_0[i] = g_yyy_0_xyyyyzzz_0[i] * fbe_0 - g_yyy_0_xyyyyzzz_1[i] * fz_be_0 + g_xyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xyyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxyyy_0_xyyyzzzz_0[i] = g_yyy_0_xyyyzzzz_0[i] * fbe_0 - g_yyy_0_xyyyzzzz_1[i] * fz_be_0 + g_xyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xyyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxyyy_0_xyyzzzzz_0[i] = g_yyy_0_xyyzzzzz_0[i] * fbe_0 - g_yyy_0_xyyzzzzz_1[i] * fz_be_0 + g_xyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xyyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxyyy_0_xyzzzzzz_0[i] = g_yyy_0_xyzzzzzz_0[i] * fbe_0 - g_yyy_0_xyzzzzzz_1[i] * fz_be_0 + g_xyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xyyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxyyy_0_xzzzzzzz_0[i] = 2.0 * g_xxy_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_xxy_0_xzzzzzzz_1[i] * fz_be_0 + g_xxyy_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxyyy_0_yyyyyyyy_0[i] = g_yyy_0_yyyyyyyy_0[i] * fbe_0 - g_yyy_0_yyyyyyyy_1[i] * fz_be_0 + g_xyyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxyyy_0_yyyyyyyz_0[i] = g_yyy_0_yyyyyyyz_0[i] * fbe_0 - g_yyy_0_yyyyyyyz_1[i] * fz_be_0 + g_xyyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxyyy_0_yyyyyyzz_0[i] = g_yyy_0_yyyyyyzz_0[i] * fbe_0 - g_yyy_0_yyyyyyzz_1[i] * fz_be_0 + g_xyyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxyyy_0_yyyyyzzz_0[i] = g_yyy_0_yyyyyzzz_0[i] * fbe_0 - g_yyy_0_yyyyyzzz_1[i] * fz_be_0 + g_xyyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxyyy_0_yyyyzzzz_0[i] = g_yyy_0_yyyyzzzz_0[i] * fbe_0 - g_yyy_0_yyyyzzzz_1[i] * fz_be_0 + g_xyyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxyyy_0_yyyzzzzz_0[i] = g_yyy_0_yyyzzzzz_0[i] * fbe_0 - g_yyy_0_yyyzzzzz_1[i] * fz_be_0 + g_xyyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxyyy_0_yyzzzzzz_0[i] = g_yyy_0_yyzzzzzz_0[i] * fbe_0 - g_yyy_0_yyzzzzzz_1[i] * fz_be_0 + g_xyyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxyyy_0_yzzzzzzz_0[i] = g_yyy_0_yzzzzzzz_0[i] * fbe_0 - g_yyy_0_yzzzzzzz_1[i] * fz_be_0 + g_xyyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxyyy_0_zzzzzzzz_0[i] = g_yyy_0_zzzzzzzz_0[i] * fbe_0 - g_yyy_0_zzzzzzzz_1[i] * fz_be_0 + g_xyyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 315-360 components of targeted buffer : HSL

    auto g_xxyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 315);

    auto g_xxyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 316);

    auto g_xxyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 317);

    auto g_xxyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 318);

    auto g_xxyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 319);

    auto g_xxyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 320);

    auto g_xxyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 321);

    auto g_xxyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 322);

    auto g_xxyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 323);

    auto g_xxyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 324);

    auto g_xxyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 325);

    auto g_xxyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 326);

    auto g_xxyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 327);

    auto g_xxyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 328);

    auto g_xxyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 329);

    auto g_xxyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 330);

    auto g_xxyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 331);

    auto g_xxyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 332);

    auto g_xxyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 333);

    auto g_xxyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 334);

    auto g_xxyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 335);

    auto g_xxyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 336);

    auto g_xxyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 337);

    auto g_xxyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 338);

    auto g_xxyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 339);

    auto g_xxyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 340);

    auto g_xxyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 341);

    auto g_xxyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 342);

    auto g_xxyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 343);

    auto g_xxyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 344);

    auto g_xxyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 345);

    auto g_xxyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 346);

    auto g_xxyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 347);

    auto g_xxyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 348);

    auto g_xxyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 349);

    auto g_xxyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 350);

    auto g_xxyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 351);

    auto g_xxyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 352);

    auto g_xxyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 353);

    auto g_xxyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 354);

    auto g_xxyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 355);

    auto g_xxyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 356);

    auto g_xxyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 357);

    auto g_xxyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 358);

    auto g_xxyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 359);

    #pragma omp simd aligned(g_xxyy_0_xxxxxxx_1, g_xxyy_0_xxxxxxxx_1, g_xxyy_0_xxxxxxxy_1, g_xxyy_0_xxxxxxxz_1, g_xxyy_0_xxxxxxy_1, g_xxyy_0_xxxxxxyy_1, g_xxyy_0_xxxxxxyz_1, g_xxyy_0_xxxxxxz_1, g_xxyy_0_xxxxxxzz_1, g_xxyy_0_xxxxxyy_1, g_xxyy_0_xxxxxyyy_1, g_xxyy_0_xxxxxyyz_1, g_xxyy_0_xxxxxyz_1, g_xxyy_0_xxxxxyzz_1, g_xxyy_0_xxxxxzz_1, g_xxyy_0_xxxxxzzz_1, g_xxyy_0_xxxxyyy_1, g_xxyy_0_xxxxyyyy_1, g_xxyy_0_xxxxyyyz_1, g_xxyy_0_xxxxyyz_1, g_xxyy_0_xxxxyyzz_1, g_xxyy_0_xxxxyzz_1, g_xxyy_0_xxxxyzzz_1, g_xxyy_0_xxxxzzz_1, g_xxyy_0_xxxxzzzz_1, g_xxyy_0_xxxyyyy_1, g_xxyy_0_xxxyyyyy_1, g_xxyy_0_xxxyyyyz_1, g_xxyy_0_xxxyyyz_1, g_xxyy_0_xxxyyyzz_1, g_xxyy_0_xxxyyzz_1, g_xxyy_0_xxxyyzzz_1, g_xxyy_0_xxxyzzz_1, g_xxyy_0_xxxyzzzz_1, g_xxyy_0_xxxzzzz_1, g_xxyy_0_xxxzzzzz_1, g_xxyy_0_xxyyyyy_1, g_xxyy_0_xxyyyyyy_1, g_xxyy_0_xxyyyyyz_1, g_xxyy_0_xxyyyyz_1, g_xxyy_0_xxyyyyzz_1, g_xxyy_0_xxyyyzz_1, g_xxyy_0_xxyyyzzz_1, g_xxyy_0_xxyyzzz_1, g_xxyy_0_xxyyzzzz_1, g_xxyy_0_xxyzzzz_1, g_xxyy_0_xxyzzzzz_1, g_xxyy_0_xxzzzzz_1, g_xxyy_0_xxzzzzzz_1, g_xxyy_0_xyyyyyy_1, g_xxyy_0_xyyyyyyy_1, g_xxyy_0_xyyyyyyz_1, g_xxyy_0_xyyyyyz_1, g_xxyy_0_xyyyyyzz_1, g_xxyy_0_xyyyyzz_1, g_xxyy_0_xyyyyzzz_1, g_xxyy_0_xyyyzzz_1, g_xxyy_0_xyyyzzzz_1, g_xxyy_0_xyyzzzz_1, g_xxyy_0_xyyzzzzz_1, g_xxyy_0_xyzzzzz_1, g_xxyy_0_xyzzzzzz_1, g_xxyy_0_xzzzzzz_1, g_xxyy_0_xzzzzzzz_1, g_xxyy_0_yyyyyyy_1, g_xxyy_0_yyyyyyyy_1, g_xxyy_0_yyyyyyyz_1, g_xxyy_0_yyyyyyz_1, g_xxyy_0_yyyyyyzz_1, g_xxyy_0_yyyyyzz_1, g_xxyy_0_yyyyyzzz_1, g_xxyy_0_yyyyzzz_1, g_xxyy_0_yyyyzzzz_1, g_xxyy_0_yyyzzzz_1, g_xxyy_0_yyyzzzzz_1, g_xxyy_0_yyzzzzz_1, g_xxyy_0_yyzzzzzz_1, g_xxyy_0_yzzzzzz_1, g_xxyy_0_yzzzzzzz_1, g_xxyy_0_zzzzzzz_1, g_xxyy_0_zzzzzzzz_1, g_xxyyz_0_xxxxxxxx_0, g_xxyyz_0_xxxxxxxy_0, g_xxyyz_0_xxxxxxxz_0, g_xxyyz_0_xxxxxxyy_0, g_xxyyz_0_xxxxxxyz_0, g_xxyyz_0_xxxxxxzz_0, g_xxyyz_0_xxxxxyyy_0, g_xxyyz_0_xxxxxyyz_0, g_xxyyz_0_xxxxxyzz_0, g_xxyyz_0_xxxxxzzz_0, g_xxyyz_0_xxxxyyyy_0, g_xxyyz_0_xxxxyyyz_0, g_xxyyz_0_xxxxyyzz_0, g_xxyyz_0_xxxxyzzz_0, g_xxyyz_0_xxxxzzzz_0, g_xxyyz_0_xxxyyyyy_0, g_xxyyz_0_xxxyyyyz_0, g_xxyyz_0_xxxyyyzz_0, g_xxyyz_0_xxxyyzzz_0, g_xxyyz_0_xxxyzzzz_0, g_xxyyz_0_xxxzzzzz_0, g_xxyyz_0_xxyyyyyy_0, g_xxyyz_0_xxyyyyyz_0, g_xxyyz_0_xxyyyyzz_0, g_xxyyz_0_xxyyyzzz_0, g_xxyyz_0_xxyyzzzz_0, g_xxyyz_0_xxyzzzzz_0, g_xxyyz_0_xxzzzzzz_0, g_xxyyz_0_xyyyyyyy_0, g_xxyyz_0_xyyyyyyz_0, g_xxyyz_0_xyyyyyzz_0, g_xxyyz_0_xyyyyzzz_0, g_xxyyz_0_xyyyzzzz_0, g_xxyyz_0_xyyzzzzz_0, g_xxyyz_0_xyzzzzzz_0, g_xxyyz_0_xzzzzzzz_0, g_xxyyz_0_yyyyyyyy_0, g_xxyyz_0_yyyyyyyz_0, g_xxyyz_0_yyyyyyzz_0, g_xxyyz_0_yyyyyzzz_0, g_xxyyz_0_yyyyzzzz_0, g_xxyyz_0_yyyzzzzz_0, g_xxyyz_0_yyzzzzzz_0, g_xxyyz_0_yzzzzzzz_0, g_xxyyz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyz_0_xxxxxxxx_0[i] = g_xxyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxxxy_0[i] = g_xxyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxxxz_0[i] = g_xxyy_0_xxxxxxx_1[i] * fi_acd_0 + g_xxyy_0_xxxxxxxz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxxyy_0[i] = g_xxyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxxyz_0[i] = g_xxyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xxyy_0_xxxxxxyz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxxzz_0[i] = 2.0 * g_xxyy_0_xxxxxxz_1[i] * fi_acd_0 + g_xxyy_0_xxxxxxzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxyyy_0[i] = g_xxyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxyyz_0[i] = g_xxyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xxyy_0_xxxxxyyz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxyzz_0[i] = 2.0 * g_xxyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xxyy_0_xxxxxyzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxxzzz_0[i] = 3.0 * g_xxyy_0_xxxxxzz_1[i] * fi_acd_0 + g_xxyy_0_xxxxxzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxyyyy_0[i] = g_xxyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxyyz_0_xxxxyyyz_0[i] = g_xxyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xxyy_0_xxxxyyyz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxyyzz_0[i] = 2.0 * g_xxyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xxyy_0_xxxxyyzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxyzzz_0[i] = 3.0 * g_xxyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xxyy_0_xxxxyzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxxzzzz_0[i] = 4.0 * g_xxyy_0_xxxxzzz_1[i] * fi_acd_0 + g_xxyy_0_xxxxzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxyyyyy_0[i] = g_xxyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxyyz_0_xxxyyyyz_0[i] = g_xxyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xxyy_0_xxxyyyyz_1[i] * wa_z[i];

        g_xxyyz_0_xxxyyyzz_0[i] = 2.0 * g_xxyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xxyy_0_xxxyyyzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxyyzzz_0[i] = 3.0 * g_xxyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xxyy_0_xxxyyzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxyzzzz_0[i] = 4.0 * g_xxyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xxyy_0_xxxyzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxxzzzzz_0[i] = 5.0 * g_xxyy_0_xxxzzzz_1[i] * fi_acd_0 + g_xxyy_0_xxxzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxyyyyyy_0[i] = g_xxyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxyyz_0_xxyyyyyz_0[i] = g_xxyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xxyy_0_xxyyyyyz_1[i] * wa_z[i];

        g_xxyyz_0_xxyyyyzz_0[i] = 2.0 * g_xxyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xxyy_0_xxyyyyzz_1[i] * wa_z[i];

        g_xxyyz_0_xxyyyzzz_0[i] = 3.0 * g_xxyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xxyy_0_xxyyyzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxyyzzzz_0[i] = 4.0 * g_xxyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xxyy_0_xxyyzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxyzzzzz_0[i] = 5.0 * g_xxyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xxyy_0_xxyzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xxzzzzzz_0[i] = 6.0 * g_xxyy_0_xxzzzzz_1[i] * fi_acd_0 + g_xxyy_0_xxzzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xyyyyyyy_0[i] = g_xxyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxyyz_0_xyyyyyyz_0[i] = g_xxyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xxyy_0_xyyyyyyz_1[i] * wa_z[i];

        g_xxyyz_0_xyyyyyzz_0[i] = 2.0 * g_xxyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xxyy_0_xyyyyyzz_1[i] * wa_z[i];

        g_xxyyz_0_xyyyyzzz_0[i] = 3.0 * g_xxyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xxyy_0_xyyyyzzz_1[i] * wa_z[i];

        g_xxyyz_0_xyyyzzzz_0[i] = 4.0 * g_xxyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xxyy_0_xyyyzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xyyzzzzz_0[i] = 5.0 * g_xxyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xxyy_0_xyyzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xyzzzzzz_0[i] = 6.0 * g_xxyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xxyy_0_xyzzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_xzzzzzzz_0[i] = 7.0 * g_xxyy_0_xzzzzzz_1[i] * fi_acd_0 + g_xxyy_0_xzzzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_yyyyyyyy_0[i] = g_xxyy_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxyyz_0_yyyyyyyz_0[i] = g_xxyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xxyy_0_yyyyyyyz_1[i] * wa_z[i];

        g_xxyyz_0_yyyyyyzz_0[i] = 2.0 * g_xxyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xxyy_0_yyyyyyzz_1[i] * wa_z[i];

        g_xxyyz_0_yyyyyzzz_0[i] = 3.0 * g_xxyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xxyy_0_yyyyyzzz_1[i] * wa_z[i];

        g_xxyyz_0_yyyyzzzz_0[i] = 4.0 * g_xxyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xxyy_0_yyyyzzzz_1[i] * wa_z[i];

        g_xxyyz_0_yyyzzzzz_0[i] = 5.0 * g_xxyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xxyy_0_yyyzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_yyzzzzzz_0[i] = 6.0 * g_xxyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xxyy_0_yyzzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_yzzzzzzz_0[i] = 7.0 * g_xxyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xxyy_0_yzzzzzzz_1[i] * wa_z[i];

        g_xxyyz_0_zzzzzzzz_0[i] = 8.0 * g_xxyy_0_zzzzzzz_1[i] * fi_acd_0 + g_xxyy_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 360-405 components of targeted buffer : HSL

    auto g_xxyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 360);

    auto g_xxyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 361);

    auto g_xxyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 362);

    auto g_xxyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 363);

    auto g_xxyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 364);

    auto g_xxyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 365);

    auto g_xxyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 366);

    auto g_xxyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 367);

    auto g_xxyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 368);

    auto g_xxyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 369);

    auto g_xxyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 370);

    auto g_xxyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 371);

    auto g_xxyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 372);

    auto g_xxyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 373);

    auto g_xxyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 374);

    auto g_xxyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 375);

    auto g_xxyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 376);

    auto g_xxyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 377);

    auto g_xxyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 378);

    auto g_xxyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 379);

    auto g_xxyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 380);

    auto g_xxyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 381);

    auto g_xxyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 382);

    auto g_xxyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 383);

    auto g_xxyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 384);

    auto g_xxyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 385);

    auto g_xxyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 386);

    auto g_xxyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 387);

    auto g_xxyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 388);

    auto g_xxyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 389);

    auto g_xxyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 390);

    auto g_xxyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 391);

    auto g_xxyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 392);

    auto g_xxyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 393);

    auto g_xxyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 394);

    auto g_xxyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 395);

    auto g_xxyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 396);

    auto g_xxyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 397);

    auto g_xxyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 398);

    auto g_xxyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 399);

    auto g_xxyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 400);

    auto g_xxyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 401);

    auto g_xxyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 402);

    auto g_xxyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 403);

    auto g_xxyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 404);

    #pragma omp simd aligned(g_xxyzz_0_xxxxxxxx_0, g_xxyzz_0_xxxxxxxy_0, g_xxyzz_0_xxxxxxxz_0, g_xxyzz_0_xxxxxxyy_0, g_xxyzz_0_xxxxxxyz_0, g_xxyzz_0_xxxxxxzz_0, g_xxyzz_0_xxxxxyyy_0, g_xxyzz_0_xxxxxyyz_0, g_xxyzz_0_xxxxxyzz_0, g_xxyzz_0_xxxxxzzz_0, g_xxyzz_0_xxxxyyyy_0, g_xxyzz_0_xxxxyyyz_0, g_xxyzz_0_xxxxyyzz_0, g_xxyzz_0_xxxxyzzz_0, g_xxyzz_0_xxxxzzzz_0, g_xxyzz_0_xxxyyyyy_0, g_xxyzz_0_xxxyyyyz_0, g_xxyzz_0_xxxyyyzz_0, g_xxyzz_0_xxxyyzzz_0, g_xxyzz_0_xxxyzzzz_0, g_xxyzz_0_xxxzzzzz_0, g_xxyzz_0_xxyyyyyy_0, g_xxyzz_0_xxyyyyyz_0, g_xxyzz_0_xxyyyyzz_0, g_xxyzz_0_xxyyyzzz_0, g_xxyzz_0_xxyyzzzz_0, g_xxyzz_0_xxyzzzzz_0, g_xxyzz_0_xxzzzzzz_0, g_xxyzz_0_xyyyyyyy_0, g_xxyzz_0_xyyyyyyz_0, g_xxyzz_0_xyyyyyzz_0, g_xxyzz_0_xyyyyzzz_0, g_xxyzz_0_xyyyzzzz_0, g_xxyzz_0_xyyzzzzz_0, g_xxyzz_0_xyzzzzzz_0, g_xxyzz_0_xzzzzzzz_0, g_xxyzz_0_yyyyyyyy_0, g_xxyzz_0_yyyyyyyz_0, g_xxyzz_0_yyyyyyzz_0, g_xxyzz_0_yyyyyzzz_0, g_xxyzz_0_yyyyzzzz_0, g_xxyzz_0_yyyzzzzz_0, g_xxyzz_0_yyzzzzzz_0, g_xxyzz_0_yzzzzzzz_0, g_xxyzz_0_zzzzzzzz_0, g_xxzz_0_xxxxxxx_1, g_xxzz_0_xxxxxxxx_1, g_xxzz_0_xxxxxxxy_1, g_xxzz_0_xxxxxxxz_1, g_xxzz_0_xxxxxxy_1, g_xxzz_0_xxxxxxyy_1, g_xxzz_0_xxxxxxyz_1, g_xxzz_0_xxxxxxz_1, g_xxzz_0_xxxxxxzz_1, g_xxzz_0_xxxxxyy_1, g_xxzz_0_xxxxxyyy_1, g_xxzz_0_xxxxxyyz_1, g_xxzz_0_xxxxxyz_1, g_xxzz_0_xxxxxyzz_1, g_xxzz_0_xxxxxzz_1, g_xxzz_0_xxxxxzzz_1, g_xxzz_0_xxxxyyy_1, g_xxzz_0_xxxxyyyy_1, g_xxzz_0_xxxxyyyz_1, g_xxzz_0_xxxxyyz_1, g_xxzz_0_xxxxyyzz_1, g_xxzz_0_xxxxyzz_1, g_xxzz_0_xxxxyzzz_1, g_xxzz_0_xxxxzzz_1, g_xxzz_0_xxxxzzzz_1, g_xxzz_0_xxxyyyy_1, g_xxzz_0_xxxyyyyy_1, g_xxzz_0_xxxyyyyz_1, g_xxzz_0_xxxyyyz_1, g_xxzz_0_xxxyyyzz_1, g_xxzz_0_xxxyyzz_1, g_xxzz_0_xxxyyzzz_1, g_xxzz_0_xxxyzzz_1, g_xxzz_0_xxxyzzzz_1, g_xxzz_0_xxxzzzz_1, g_xxzz_0_xxxzzzzz_1, g_xxzz_0_xxyyyyy_1, g_xxzz_0_xxyyyyyy_1, g_xxzz_0_xxyyyyyz_1, g_xxzz_0_xxyyyyz_1, g_xxzz_0_xxyyyyzz_1, g_xxzz_0_xxyyyzz_1, g_xxzz_0_xxyyyzzz_1, g_xxzz_0_xxyyzzz_1, g_xxzz_0_xxyyzzzz_1, g_xxzz_0_xxyzzzz_1, g_xxzz_0_xxyzzzzz_1, g_xxzz_0_xxzzzzz_1, g_xxzz_0_xxzzzzzz_1, g_xxzz_0_xyyyyyy_1, g_xxzz_0_xyyyyyyy_1, g_xxzz_0_xyyyyyyz_1, g_xxzz_0_xyyyyyz_1, g_xxzz_0_xyyyyyzz_1, g_xxzz_0_xyyyyzz_1, g_xxzz_0_xyyyyzzz_1, g_xxzz_0_xyyyzzz_1, g_xxzz_0_xyyyzzzz_1, g_xxzz_0_xyyzzzz_1, g_xxzz_0_xyyzzzzz_1, g_xxzz_0_xyzzzzz_1, g_xxzz_0_xyzzzzzz_1, g_xxzz_0_xzzzzzz_1, g_xxzz_0_xzzzzzzz_1, g_xxzz_0_yyyyyyy_1, g_xxzz_0_yyyyyyyy_1, g_xxzz_0_yyyyyyyz_1, g_xxzz_0_yyyyyyz_1, g_xxzz_0_yyyyyyzz_1, g_xxzz_0_yyyyyzz_1, g_xxzz_0_yyyyyzzz_1, g_xxzz_0_yyyyzzz_1, g_xxzz_0_yyyyzzzz_1, g_xxzz_0_yyyzzzz_1, g_xxzz_0_yyyzzzzz_1, g_xxzz_0_yyzzzzz_1, g_xxzz_0_yyzzzzzz_1, g_xxzz_0_yzzzzzz_1, g_xxzz_0_yzzzzzzz_1, g_xxzz_0_zzzzzzz_1, g_xxzz_0_zzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzz_0_xxxxxxxx_0[i] = g_xxzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxxxy_0[i] = g_xxzz_0_xxxxxxx_1[i] * fi_acd_0 + g_xxzz_0_xxxxxxxy_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxxxz_0[i] = g_xxzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxxyy_0[i] = 2.0 * g_xxzz_0_xxxxxxy_1[i] * fi_acd_0 + g_xxzz_0_xxxxxxyy_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxxyz_0[i] = g_xxzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxxzz_0[i] = g_xxzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxyyy_0[i] = 3.0 * g_xxzz_0_xxxxxyy_1[i] * fi_acd_0 + g_xxzz_0_xxxxxyyy_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxyyz_0[i] = 2.0 * g_xxzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxyzz_0[i] = g_xxzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxxzzz_0[i] = g_xxzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxyyyy_0[i] = 4.0 * g_xxzz_0_xxxxyyy_1[i] * fi_acd_0 + g_xxzz_0_xxxxyyyy_1[i] * wa_y[i];

        g_xxyzz_0_xxxxyyyz_0[i] = 3.0 * g_xxzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxyyzz_0[i] = 2.0 * g_xxzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxyzzz_0[i] = g_xxzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxxzzzz_0[i] = g_xxzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxyyyyy_0[i] = 5.0 * g_xxzz_0_xxxyyyy_1[i] * fi_acd_0 + g_xxzz_0_xxxyyyyy_1[i] * wa_y[i];

        g_xxyzz_0_xxxyyyyz_0[i] = 4.0 * g_xxzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxyzz_0_xxxyyyzz_0[i] = 3.0 * g_xxzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxyyzzz_0[i] = 2.0 * g_xxzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxyzzzz_0[i] = g_xxzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxxzzzzz_0[i] = g_xxzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxyyyyyy_0[i] = 6.0 * g_xxzz_0_xxyyyyy_1[i] * fi_acd_0 + g_xxzz_0_xxyyyyyy_1[i] * wa_y[i];

        g_xxyzz_0_xxyyyyyz_0[i] = 5.0 * g_xxzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxyzz_0_xxyyyyzz_0[i] = 4.0 * g_xxzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxyzz_0_xxyyyzzz_0[i] = 3.0 * g_xxzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxyyzzzz_0[i] = 2.0 * g_xxzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxyzzzzz_0[i] = g_xxzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xxzzzzzz_0[i] = g_xxzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xyyyyyyy_0[i] = 7.0 * g_xxzz_0_xyyyyyy_1[i] * fi_acd_0 + g_xxzz_0_xyyyyyyy_1[i] * wa_y[i];

        g_xxyzz_0_xyyyyyyz_0[i] = 6.0 * g_xxzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxyzz_0_xyyyyyzz_0[i] = 5.0 * g_xxzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxyzz_0_xyyyyzzz_0[i] = 4.0 * g_xxzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxyzz_0_xyyyzzzz_0[i] = 3.0 * g_xxzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xyyzzzzz_0[i] = 2.0 * g_xxzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xyzzzzzz_0[i] = g_xxzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_xzzzzzzz_0[i] = g_xxzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_yyyyyyyy_0[i] = 8.0 * g_xxzz_0_yyyyyyy_1[i] * fi_acd_0 + g_xxzz_0_yyyyyyyy_1[i] * wa_y[i];

        g_xxyzz_0_yyyyyyyz_0[i] = 7.0 * g_xxzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxyzz_0_yyyyyyzz_0[i] = 6.0 * g_xxzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxyzz_0_yyyyyzzz_0[i] = 5.0 * g_xxzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxyzz_0_yyyyzzzz_0[i] = 4.0 * g_xxzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxyzz_0_yyyzzzzz_0[i] = 3.0 * g_xxzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_yyzzzzzz_0[i] = 2.0 * g_xxzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_yzzzzzzz_0[i] = g_xxzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxyzz_0_zzzzzzzz_0[i] = g_xxzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 405-450 components of targeted buffer : HSL

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

    #pragma omp simd aligned(g_xxz_0_xxxxxxxx_0, g_xxz_0_xxxxxxxx_1, g_xxz_0_xxxxxxxy_0, g_xxz_0_xxxxxxxy_1, g_xxz_0_xxxxxxyy_0, g_xxz_0_xxxxxxyy_1, g_xxz_0_xxxxxyyy_0, g_xxz_0_xxxxxyyy_1, g_xxz_0_xxxxyyyy_0, g_xxz_0_xxxxyyyy_1, g_xxz_0_xxxyyyyy_0, g_xxz_0_xxxyyyyy_1, g_xxz_0_xxyyyyyy_0, g_xxz_0_xxyyyyyy_1, g_xxz_0_xyyyyyyy_0, g_xxz_0_xyyyyyyy_1, g_xxzz_0_xxxxxxxx_1, g_xxzz_0_xxxxxxxy_1, g_xxzz_0_xxxxxxyy_1, g_xxzz_0_xxxxxyyy_1, g_xxzz_0_xxxxyyyy_1, g_xxzz_0_xxxyyyyy_1, g_xxzz_0_xxyyyyyy_1, g_xxzz_0_xyyyyyyy_1, g_xxzzz_0_xxxxxxxx_0, g_xxzzz_0_xxxxxxxy_0, g_xxzzz_0_xxxxxxxz_0, g_xxzzz_0_xxxxxxyy_0, g_xxzzz_0_xxxxxxyz_0, g_xxzzz_0_xxxxxxzz_0, g_xxzzz_0_xxxxxyyy_0, g_xxzzz_0_xxxxxyyz_0, g_xxzzz_0_xxxxxyzz_0, g_xxzzz_0_xxxxxzzz_0, g_xxzzz_0_xxxxyyyy_0, g_xxzzz_0_xxxxyyyz_0, g_xxzzz_0_xxxxyyzz_0, g_xxzzz_0_xxxxyzzz_0, g_xxzzz_0_xxxxzzzz_0, g_xxzzz_0_xxxyyyyy_0, g_xxzzz_0_xxxyyyyz_0, g_xxzzz_0_xxxyyyzz_0, g_xxzzz_0_xxxyyzzz_0, g_xxzzz_0_xxxyzzzz_0, g_xxzzz_0_xxxzzzzz_0, g_xxzzz_0_xxyyyyyy_0, g_xxzzz_0_xxyyyyyz_0, g_xxzzz_0_xxyyyyzz_0, g_xxzzz_0_xxyyyzzz_0, g_xxzzz_0_xxyyzzzz_0, g_xxzzz_0_xxyzzzzz_0, g_xxzzz_0_xxzzzzzz_0, g_xxzzz_0_xyyyyyyy_0, g_xxzzz_0_xyyyyyyz_0, g_xxzzz_0_xyyyyyzz_0, g_xxzzz_0_xyyyyzzz_0, g_xxzzz_0_xyyyzzzz_0, g_xxzzz_0_xyyzzzzz_0, g_xxzzz_0_xyzzzzzz_0, g_xxzzz_0_xzzzzzzz_0, g_xxzzz_0_yyyyyyyy_0, g_xxzzz_0_yyyyyyyz_0, g_xxzzz_0_yyyyyyzz_0, g_xxzzz_0_yyyyyzzz_0, g_xxzzz_0_yyyyzzzz_0, g_xxzzz_0_yyyzzzzz_0, g_xxzzz_0_yyzzzzzz_0, g_xxzzz_0_yzzzzzzz_0, g_xxzzz_0_zzzzzzzz_0, g_xzzz_0_xxxxxxxz_1, g_xzzz_0_xxxxxxyz_1, g_xzzz_0_xxxxxxz_1, g_xzzz_0_xxxxxxzz_1, g_xzzz_0_xxxxxyyz_1, g_xzzz_0_xxxxxyz_1, g_xzzz_0_xxxxxyzz_1, g_xzzz_0_xxxxxzz_1, g_xzzz_0_xxxxxzzz_1, g_xzzz_0_xxxxyyyz_1, g_xzzz_0_xxxxyyz_1, g_xzzz_0_xxxxyyzz_1, g_xzzz_0_xxxxyzz_1, g_xzzz_0_xxxxyzzz_1, g_xzzz_0_xxxxzzz_1, g_xzzz_0_xxxxzzzz_1, g_xzzz_0_xxxyyyyz_1, g_xzzz_0_xxxyyyz_1, g_xzzz_0_xxxyyyzz_1, g_xzzz_0_xxxyyzz_1, g_xzzz_0_xxxyyzzz_1, g_xzzz_0_xxxyzzz_1, g_xzzz_0_xxxyzzzz_1, g_xzzz_0_xxxzzzz_1, g_xzzz_0_xxxzzzzz_1, g_xzzz_0_xxyyyyyz_1, g_xzzz_0_xxyyyyz_1, g_xzzz_0_xxyyyyzz_1, g_xzzz_0_xxyyyzz_1, g_xzzz_0_xxyyyzzz_1, g_xzzz_0_xxyyzzz_1, g_xzzz_0_xxyyzzzz_1, g_xzzz_0_xxyzzzz_1, g_xzzz_0_xxyzzzzz_1, g_xzzz_0_xxzzzzz_1, g_xzzz_0_xxzzzzzz_1, g_xzzz_0_xyyyyyyz_1, g_xzzz_0_xyyyyyz_1, g_xzzz_0_xyyyyyzz_1, g_xzzz_0_xyyyyzz_1, g_xzzz_0_xyyyyzzz_1, g_xzzz_0_xyyyzzz_1, g_xzzz_0_xyyyzzzz_1, g_xzzz_0_xyyzzzz_1, g_xzzz_0_xyyzzzzz_1, g_xzzz_0_xyzzzzz_1, g_xzzz_0_xyzzzzzz_1, g_xzzz_0_xzzzzzz_1, g_xzzz_0_xzzzzzzz_1, g_xzzz_0_yyyyyyyy_1, g_xzzz_0_yyyyyyyz_1, g_xzzz_0_yyyyyyz_1, g_xzzz_0_yyyyyyzz_1, g_xzzz_0_yyyyyzz_1, g_xzzz_0_yyyyyzzz_1, g_xzzz_0_yyyyzzz_1, g_xzzz_0_yyyyzzzz_1, g_xzzz_0_yyyzzzz_1, g_xzzz_0_yyyzzzzz_1, g_xzzz_0_yyzzzzz_1, g_xzzz_0_yyzzzzzz_1, g_xzzz_0_yzzzzzz_1, g_xzzz_0_yzzzzzzz_1, g_xzzz_0_zzzzzzz_1, g_xzzz_0_zzzzzzzz_1, g_zzz_0_xxxxxxxz_0, g_zzz_0_xxxxxxxz_1, g_zzz_0_xxxxxxyz_0, g_zzz_0_xxxxxxyz_1, g_zzz_0_xxxxxxzz_0, g_zzz_0_xxxxxxzz_1, g_zzz_0_xxxxxyyz_0, g_zzz_0_xxxxxyyz_1, g_zzz_0_xxxxxyzz_0, g_zzz_0_xxxxxyzz_1, g_zzz_0_xxxxxzzz_0, g_zzz_0_xxxxxzzz_1, g_zzz_0_xxxxyyyz_0, g_zzz_0_xxxxyyyz_1, g_zzz_0_xxxxyyzz_0, g_zzz_0_xxxxyyzz_1, g_zzz_0_xxxxyzzz_0, g_zzz_0_xxxxyzzz_1, g_zzz_0_xxxxzzzz_0, g_zzz_0_xxxxzzzz_1, g_zzz_0_xxxyyyyz_0, g_zzz_0_xxxyyyyz_1, g_zzz_0_xxxyyyzz_0, g_zzz_0_xxxyyyzz_1, g_zzz_0_xxxyyzzz_0, g_zzz_0_xxxyyzzz_1, g_zzz_0_xxxyzzzz_0, g_zzz_0_xxxyzzzz_1, g_zzz_0_xxxzzzzz_0, g_zzz_0_xxxzzzzz_1, g_zzz_0_xxyyyyyz_0, g_zzz_0_xxyyyyyz_1, g_zzz_0_xxyyyyzz_0, g_zzz_0_xxyyyyzz_1, g_zzz_0_xxyyyzzz_0, g_zzz_0_xxyyyzzz_1, g_zzz_0_xxyyzzzz_0, g_zzz_0_xxyyzzzz_1, g_zzz_0_xxyzzzzz_0, g_zzz_0_xxyzzzzz_1, g_zzz_0_xxzzzzzz_0, g_zzz_0_xxzzzzzz_1, g_zzz_0_xyyyyyyz_0, g_zzz_0_xyyyyyyz_1, g_zzz_0_xyyyyyzz_0, g_zzz_0_xyyyyyzz_1, g_zzz_0_xyyyyzzz_0, g_zzz_0_xyyyyzzz_1, g_zzz_0_xyyyzzzz_0, g_zzz_0_xyyyzzzz_1, g_zzz_0_xyyzzzzz_0, g_zzz_0_xyyzzzzz_1, g_zzz_0_xyzzzzzz_0, g_zzz_0_xyzzzzzz_1, g_zzz_0_xzzzzzzz_0, g_zzz_0_xzzzzzzz_1, g_zzz_0_yyyyyyyy_0, g_zzz_0_yyyyyyyy_1, g_zzz_0_yyyyyyyz_0, g_zzz_0_yyyyyyyz_1, g_zzz_0_yyyyyyzz_0, g_zzz_0_yyyyyyzz_1, g_zzz_0_yyyyyzzz_0, g_zzz_0_yyyyyzzz_1, g_zzz_0_yyyyzzzz_0, g_zzz_0_yyyyzzzz_1, g_zzz_0_yyyzzzzz_0, g_zzz_0_yyyzzzzz_1, g_zzz_0_yyzzzzzz_0, g_zzz_0_yyzzzzzz_1, g_zzz_0_yzzzzzzz_0, g_zzz_0_yzzzzzzz_1, g_zzz_0_zzzzzzzz_0, g_zzz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzz_0_xxxxxxxx_0[i] = 2.0 * g_xxz_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxxxxx_1[i] * fz_be_0 + g_xxzz_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxzzz_0_xxxxxxxy_0[i] = 2.0 * g_xxz_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxxxxy_1[i] * fz_be_0 + g_xxzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxzzz_0_xxxxxxxz_0[i] = g_zzz_0_xxxxxxxz_0[i] * fbe_0 - g_zzz_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xzzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxxxyy_0[i] = 2.0 * g_xxz_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxxxyy_1[i] * fz_be_0 + g_xxzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxzzz_0_xxxxxxyz_0[i] = g_zzz_0_xxxxxxyz_0[i] * fbe_0 - g_zzz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxxxzz_0[i] = g_zzz_0_xxxxxxzz_0[i] * fbe_0 - g_zzz_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xzzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxxyyy_0[i] = 2.0 * g_xxz_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxxyyy_1[i] * fz_be_0 + g_xxzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxzzz_0_xxxxxyyz_0[i] = g_zzz_0_xxxxxyyz_0[i] * fbe_0 - g_zzz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxxyzz_0[i] = g_zzz_0_xxxxxyzz_0[i] * fbe_0 - g_zzz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxxzzz_0[i] = g_zzz_0_xxxxxzzz_0[i] * fbe_0 - g_zzz_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xzzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxyyyy_0[i] = 2.0 * g_xxz_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxxyyyy_1[i] * fz_be_0 + g_xxzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxzzz_0_xxxxyyyz_0[i] = g_zzz_0_xxxxyyyz_0[i] * fbe_0 - g_zzz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxyyzz_0[i] = g_zzz_0_xxxxyyzz_0[i] * fbe_0 - g_zzz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxyzzz_0[i] = g_zzz_0_xxxxyzzz_0[i] * fbe_0 - g_zzz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxxzzzz_0[i] = g_zzz_0_xxxxzzzz_0[i] * fbe_0 - g_zzz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xzzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxyyyyy_0[i] = 2.0 * g_xxz_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxxyyyyy_1[i] * fz_be_0 + g_xxzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxzzz_0_xxxyyyyz_0[i] = g_zzz_0_xxxyyyyz_0[i] * fbe_0 - g_zzz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxzzz_0_xxxyyyzz_0[i] = g_zzz_0_xxxyyyzz_0[i] * fbe_0 - g_zzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxyyzzz_0[i] = g_zzz_0_xxxyyzzz_0[i] * fbe_0 - g_zzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxyzzzz_0[i] = g_zzz_0_xxxyzzzz_0[i] * fbe_0 - g_zzz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxxzzzzz_0[i] = g_zzz_0_xxxzzzzz_0[i] * fbe_0 - g_zzz_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xzzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxyyyyyy_0[i] = 2.0 * g_xxz_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xxyyyyyy_1[i] * fz_be_0 + g_xxzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxzzz_0_xxyyyyyz_0[i] = g_zzz_0_xxyyyyyz_0[i] * fbe_0 - g_zzz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxzzz_0_xxyyyyzz_0[i] = g_zzz_0_xxyyyyzz_0[i] * fbe_0 - g_zzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxzzz_0_xxyyyzzz_0[i] = g_zzz_0_xxyyyzzz_0[i] * fbe_0 - g_zzz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxyyzzzz_0[i] = g_zzz_0_xxyyzzzz_0[i] * fbe_0 - g_zzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxyzzzzz_0[i] = g_zzz_0_xxyzzzzz_0[i] * fbe_0 - g_zzz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xxzzzzzz_0[i] = g_zzz_0_xxzzzzzz_0[i] * fbe_0 - g_zzz_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xzzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xyyyyyyy_0[i] = 2.0 * g_xxz_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_xxz_0_xyyyyyyy_1[i] * fz_be_0 + g_xxzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxzzz_0_xyyyyyyz_0[i] = g_zzz_0_xyyyyyyz_0[i] * fbe_0 - g_zzz_0_xyyyyyyz_1[i] * fz_be_0 + g_xzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxzzz_0_xyyyyyzz_0[i] = g_zzz_0_xyyyyyzz_0[i] * fbe_0 - g_zzz_0_xyyyyyzz_1[i] * fz_be_0 + g_xzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxzzz_0_xyyyyzzz_0[i] = g_zzz_0_xyyyyzzz_0[i] * fbe_0 - g_zzz_0_xyyyyzzz_1[i] * fz_be_0 + g_xzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxzzz_0_xyyyzzzz_0[i] = g_zzz_0_xyyyzzzz_0[i] * fbe_0 - g_zzz_0_xyyyzzzz_1[i] * fz_be_0 + g_xzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xyyzzzzz_0[i] = g_zzz_0_xyyzzzzz_0[i] * fbe_0 - g_zzz_0_xyyzzzzz_1[i] * fz_be_0 + g_xzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xyzzzzzz_0[i] = g_zzz_0_xyzzzzzz_0[i] * fbe_0 - g_zzz_0_xyzzzzzz_1[i] * fz_be_0 + g_xzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_xzzzzzzz_0[i] = g_zzz_0_xzzzzzzz_0[i] * fbe_0 - g_zzz_0_xzzzzzzz_1[i] * fz_be_0 + g_xzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xzzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_yyyyyyyy_0[i] = g_zzz_0_yyyyyyyy_0[i] * fbe_0 - g_zzz_0_yyyyyyyy_1[i] * fz_be_0 + g_xzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxzzz_0_yyyyyyyz_0[i] = g_zzz_0_yyyyyyyz_0[i] * fbe_0 - g_zzz_0_yyyyyyyz_1[i] * fz_be_0 + g_xzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxzzz_0_yyyyyyzz_0[i] = g_zzz_0_yyyyyyzz_0[i] * fbe_0 - g_zzz_0_yyyyyyzz_1[i] * fz_be_0 + g_xzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxzzz_0_yyyyyzzz_0[i] = g_zzz_0_yyyyyzzz_0[i] * fbe_0 - g_zzz_0_yyyyyzzz_1[i] * fz_be_0 + g_xzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxzzz_0_yyyyzzzz_0[i] = g_zzz_0_yyyyzzzz_0[i] * fbe_0 - g_zzz_0_yyyyzzzz_1[i] * fz_be_0 + g_xzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxzzz_0_yyyzzzzz_0[i] = g_zzz_0_yyyzzzzz_0[i] * fbe_0 - g_zzz_0_yyyzzzzz_1[i] * fz_be_0 + g_xzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_yyzzzzzz_0[i] = g_zzz_0_yyzzzzzz_0[i] * fbe_0 - g_zzz_0_yyzzzzzz_1[i] * fz_be_0 + g_xzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_yzzzzzzz_0[i] = g_zzz_0_yzzzzzzz_0[i] * fbe_0 - g_zzz_0_yzzzzzzz_1[i] * fz_be_0 + g_xzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxzzz_0_zzzzzzzz_0[i] = g_zzz_0_zzzzzzzz_0[i] * fbe_0 - g_zzz_0_zzzzzzzz_1[i] * fz_be_0 + g_xzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 450-495 components of targeted buffer : HSL

    auto g_xyyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 450);

    auto g_xyyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 451);

    auto g_xyyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 452);

    auto g_xyyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 453);

    auto g_xyyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 454);

    auto g_xyyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 455);

    auto g_xyyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 456);

    auto g_xyyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 457);

    auto g_xyyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 458);

    auto g_xyyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 459);

    auto g_xyyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 460);

    auto g_xyyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 461);

    auto g_xyyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 462);

    auto g_xyyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 463);

    auto g_xyyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 464);

    auto g_xyyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 465);

    auto g_xyyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 466);

    auto g_xyyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 467);

    auto g_xyyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 468);

    auto g_xyyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 469);

    auto g_xyyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 470);

    auto g_xyyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 471);

    auto g_xyyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 472);

    auto g_xyyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 473);

    auto g_xyyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 474);

    auto g_xyyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 475);

    auto g_xyyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 476);

    auto g_xyyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 477);

    auto g_xyyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 478);

    auto g_xyyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 479);

    auto g_xyyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 480);

    auto g_xyyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 481);

    auto g_xyyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 482);

    auto g_xyyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 483);

    auto g_xyyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 484);

    auto g_xyyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 485);

    auto g_xyyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 486);

    auto g_xyyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 487);

    auto g_xyyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 488);

    auto g_xyyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 489);

    auto g_xyyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 490);

    auto g_xyyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 491);

    auto g_xyyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 492);

    auto g_xyyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 493);

    auto g_xyyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 494);

    #pragma omp simd aligned(g_xyyyy_0_xxxxxxxx_0, g_xyyyy_0_xxxxxxxy_0, g_xyyyy_0_xxxxxxxz_0, g_xyyyy_0_xxxxxxyy_0, g_xyyyy_0_xxxxxxyz_0, g_xyyyy_0_xxxxxxzz_0, g_xyyyy_0_xxxxxyyy_0, g_xyyyy_0_xxxxxyyz_0, g_xyyyy_0_xxxxxyzz_0, g_xyyyy_0_xxxxxzzz_0, g_xyyyy_0_xxxxyyyy_0, g_xyyyy_0_xxxxyyyz_0, g_xyyyy_0_xxxxyyzz_0, g_xyyyy_0_xxxxyzzz_0, g_xyyyy_0_xxxxzzzz_0, g_xyyyy_0_xxxyyyyy_0, g_xyyyy_0_xxxyyyyz_0, g_xyyyy_0_xxxyyyzz_0, g_xyyyy_0_xxxyyzzz_0, g_xyyyy_0_xxxyzzzz_0, g_xyyyy_0_xxxzzzzz_0, g_xyyyy_0_xxyyyyyy_0, g_xyyyy_0_xxyyyyyz_0, g_xyyyy_0_xxyyyyzz_0, g_xyyyy_0_xxyyyzzz_0, g_xyyyy_0_xxyyzzzz_0, g_xyyyy_0_xxyzzzzz_0, g_xyyyy_0_xxzzzzzz_0, g_xyyyy_0_xyyyyyyy_0, g_xyyyy_0_xyyyyyyz_0, g_xyyyy_0_xyyyyyzz_0, g_xyyyy_0_xyyyyzzz_0, g_xyyyy_0_xyyyzzzz_0, g_xyyyy_0_xyyzzzzz_0, g_xyyyy_0_xyzzzzzz_0, g_xyyyy_0_xzzzzzzz_0, g_xyyyy_0_yyyyyyyy_0, g_xyyyy_0_yyyyyyyz_0, g_xyyyy_0_yyyyyyzz_0, g_xyyyy_0_yyyyyzzz_0, g_xyyyy_0_yyyyzzzz_0, g_xyyyy_0_yyyzzzzz_0, g_xyyyy_0_yyzzzzzz_0, g_xyyyy_0_yzzzzzzz_0, g_xyyyy_0_zzzzzzzz_0, g_yyyy_0_xxxxxxx_1, g_yyyy_0_xxxxxxxx_1, g_yyyy_0_xxxxxxxy_1, g_yyyy_0_xxxxxxxz_1, g_yyyy_0_xxxxxxy_1, g_yyyy_0_xxxxxxyy_1, g_yyyy_0_xxxxxxyz_1, g_yyyy_0_xxxxxxz_1, g_yyyy_0_xxxxxxzz_1, g_yyyy_0_xxxxxyy_1, g_yyyy_0_xxxxxyyy_1, g_yyyy_0_xxxxxyyz_1, g_yyyy_0_xxxxxyz_1, g_yyyy_0_xxxxxyzz_1, g_yyyy_0_xxxxxzz_1, g_yyyy_0_xxxxxzzz_1, g_yyyy_0_xxxxyyy_1, g_yyyy_0_xxxxyyyy_1, g_yyyy_0_xxxxyyyz_1, g_yyyy_0_xxxxyyz_1, g_yyyy_0_xxxxyyzz_1, g_yyyy_0_xxxxyzz_1, g_yyyy_0_xxxxyzzz_1, g_yyyy_0_xxxxzzz_1, g_yyyy_0_xxxxzzzz_1, g_yyyy_0_xxxyyyy_1, g_yyyy_0_xxxyyyyy_1, g_yyyy_0_xxxyyyyz_1, g_yyyy_0_xxxyyyz_1, g_yyyy_0_xxxyyyzz_1, g_yyyy_0_xxxyyzz_1, g_yyyy_0_xxxyyzzz_1, g_yyyy_0_xxxyzzz_1, g_yyyy_0_xxxyzzzz_1, g_yyyy_0_xxxzzzz_1, g_yyyy_0_xxxzzzzz_1, g_yyyy_0_xxyyyyy_1, g_yyyy_0_xxyyyyyy_1, g_yyyy_0_xxyyyyyz_1, g_yyyy_0_xxyyyyz_1, g_yyyy_0_xxyyyyzz_1, g_yyyy_0_xxyyyzz_1, g_yyyy_0_xxyyyzzz_1, g_yyyy_0_xxyyzzz_1, g_yyyy_0_xxyyzzzz_1, g_yyyy_0_xxyzzzz_1, g_yyyy_0_xxyzzzzz_1, g_yyyy_0_xxzzzzz_1, g_yyyy_0_xxzzzzzz_1, g_yyyy_0_xyyyyyy_1, g_yyyy_0_xyyyyyyy_1, g_yyyy_0_xyyyyyyz_1, g_yyyy_0_xyyyyyz_1, g_yyyy_0_xyyyyyzz_1, g_yyyy_0_xyyyyzz_1, g_yyyy_0_xyyyyzzz_1, g_yyyy_0_xyyyzzz_1, g_yyyy_0_xyyyzzzz_1, g_yyyy_0_xyyzzzz_1, g_yyyy_0_xyyzzzzz_1, g_yyyy_0_xyzzzzz_1, g_yyyy_0_xyzzzzzz_1, g_yyyy_0_xzzzzzz_1, g_yyyy_0_xzzzzzzz_1, g_yyyy_0_yyyyyyy_1, g_yyyy_0_yyyyyyyy_1, g_yyyy_0_yyyyyyyz_1, g_yyyy_0_yyyyyyz_1, g_yyyy_0_yyyyyyzz_1, g_yyyy_0_yyyyyzz_1, g_yyyy_0_yyyyyzzz_1, g_yyyy_0_yyyyzzz_1, g_yyyy_0_yyyyzzzz_1, g_yyyy_0_yyyzzzz_1, g_yyyy_0_yyyzzzzz_1, g_yyyy_0_yyzzzzz_1, g_yyyy_0_yyzzzzzz_1, g_yyyy_0_yzzzzzz_1, g_yyyy_0_yzzzzzzz_1, g_yyyy_0_zzzzzzz_1, g_yyyy_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyy_0_xxxxxxxx_0[i] = 8.0 * g_yyyy_0_xxxxxxx_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxxx_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxxxy_0[i] = 7.0 * g_yyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxxxz_0[i] = 7.0 * g_yyyy_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxxyy_0[i] = 6.0 * g_yyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxxyz_0[i] = 6.0 * g_yyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxxzz_0[i] = 6.0 * g_yyyy_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxyyy_0[i] = 5.0 * g_yyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxyyz_0[i] = 5.0 * g_yyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxyzz_0[i] = 5.0 * g_yyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxxzzz_0[i] = 5.0 * g_yyyy_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxyyyy_0[i] = 4.0 * g_yyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyyyy_0_xxxxyyyz_0[i] = 4.0 * g_yyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxyyzz_0[i] = 4.0 * g_yyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxyzzz_0[i] = 4.0 * g_yyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxxzzzz_0[i] = 4.0 * g_yyyy_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxyyyyy_0[i] = 3.0 * g_yyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyyyy_0_xxxyyyyz_0[i] = 3.0 * g_yyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyyy_0_xxxyyyzz_0[i] = 3.0 * g_yyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxyyzzz_0[i] = 3.0 * g_yyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxyzzzz_0[i] = 3.0 * g_yyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxxzzzzz_0[i] = 3.0 * g_yyyy_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxyyyyyy_0[i] = 2.0 * g_yyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyyyy_0_xxyyyyyz_0[i] = 2.0 * g_yyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyyy_0_xxyyyyzz_0[i] = 2.0 * g_yyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyyy_0_xxyyyzzz_0[i] = 2.0 * g_yyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxyyzzzz_0[i] = 2.0 * g_yyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxyzzzzz_0[i] = 2.0 * g_yyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xxzzzzzz_0[i] = 2.0 * g_yyyy_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xyyyyyyy_0[i] = g_yyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyyyy_0_xyyyyyyz_0[i] = g_yyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyyy_0_xyyyyyzz_0[i] = g_yyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyyy_0_xyyyyzzz_0[i] = g_yyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyyy_0_xyyyzzzz_0[i] = g_yyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xyyzzzzz_0[i] = g_yyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xyzzzzzz_0[i] = g_yyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_xzzzzzzz_0[i] = g_yyyy_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_yyyyyyyy_0[i] = g_yyyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyyy_0_yyyyyyyz_0[i] = g_yyyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyyy_0_yyyyyyzz_0[i] = g_yyyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyyy_0_yyyyyzzz_0[i] = g_yyyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyyy_0_yyyyzzzz_0[i] = g_yyyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyyy_0_yyyzzzzz_0[i] = g_yyyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_yyzzzzzz_0[i] = g_yyyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_yzzzzzzz_0[i] = g_yyyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyyy_0_zzzzzzzz_0[i] = g_yyyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 495-540 components of targeted buffer : HSL

    auto g_xyyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 495);

    auto g_xyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 496);

    auto g_xyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 497);

    auto g_xyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 498);

    auto g_xyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 499);

    auto g_xyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 500);

    auto g_xyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 501);

    auto g_xyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 502);

    auto g_xyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 503);

    auto g_xyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 504);

    auto g_xyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 505);

    auto g_xyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 506);

    auto g_xyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 507);

    auto g_xyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 508);

    auto g_xyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 509);

    auto g_xyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 510);

    auto g_xyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 511);

    auto g_xyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 512);

    auto g_xyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 513);

    auto g_xyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 514);

    auto g_xyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 515);

    auto g_xyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 516);

    auto g_xyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 517);

    auto g_xyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 518);

    auto g_xyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 519);

    auto g_xyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 520);

    auto g_xyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 521);

    auto g_xyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 522);

    auto g_xyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 523);

    auto g_xyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 524);

    auto g_xyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 525);

    auto g_xyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 526);

    auto g_xyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 527);

    auto g_xyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 528);

    auto g_xyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 529);

    auto g_xyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 530);

    auto g_xyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 531);

    auto g_xyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 532);

    auto g_xyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 533);

    auto g_xyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 534);

    auto g_xyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 535);

    auto g_xyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 536);

    auto g_xyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 537);

    auto g_xyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 538);

    auto g_xyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 539);

    #pragma omp simd aligned(g_xyyy_0_xxxxxxxx_1, g_xyyy_0_xxxxxxxy_1, g_xyyy_0_xxxxxxyy_1, g_xyyy_0_xxxxxyyy_1, g_xyyy_0_xxxxyyyy_1, g_xyyy_0_xxxyyyyy_1, g_xyyy_0_xxyyyyyy_1, g_xyyy_0_xyyyyyyy_1, g_xyyyz_0_xxxxxxxx_0, g_xyyyz_0_xxxxxxxy_0, g_xyyyz_0_xxxxxxxz_0, g_xyyyz_0_xxxxxxyy_0, g_xyyyz_0_xxxxxxyz_0, g_xyyyz_0_xxxxxxzz_0, g_xyyyz_0_xxxxxyyy_0, g_xyyyz_0_xxxxxyyz_0, g_xyyyz_0_xxxxxyzz_0, g_xyyyz_0_xxxxxzzz_0, g_xyyyz_0_xxxxyyyy_0, g_xyyyz_0_xxxxyyyz_0, g_xyyyz_0_xxxxyyzz_0, g_xyyyz_0_xxxxyzzz_0, g_xyyyz_0_xxxxzzzz_0, g_xyyyz_0_xxxyyyyy_0, g_xyyyz_0_xxxyyyyz_0, g_xyyyz_0_xxxyyyzz_0, g_xyyyz_0_xxxyyzzz_0, g_xyyyz_0_xxxyzzzz_0, g_xyyyz_0_xxxzzzzz_0, g_xyyyz_0_xxyyyyyy_0, g_xyyyz_0_xxyyyyyz_0, g_xyyyz_0_xxyyyyzz_0, g_xyyyz_0_xxyyyzzz_0, g_xyyyz_0_xxyyzzzz_0, g_xyyyz_0_xxyzzzzz_0, g_xyyyz_0_xxzzzzzz_0, g_xyyyz_0_xyyyyyyy_0, g_xyyyz_0_xyyyyyyz_0, g_xyyyz_0_xyyyyyzz_0, g_xyyyz_0_xyyyyzzz_0, g_xyyyz_0_xyyyzzzz_0, g_xyyyz_0_xyyzzzzz_0, g_xyyyz_0_xyzzzzzz_0, g_xyyyz_0_xzzzzzzz_0, g_xyyyz_0_yyyyyyyy_0, g_xyyyz_0_yyyyyyyz_0, g_xyyyz_0_yyyyyyzz_0, g_xyyyz_0_yyyyyzzz_0, g_xyyyz_0_yyyyzzzz_0, g_xyyyz_0_yyyzzzzz_0, g_xyyyz_0_yyzzzzzz_0, g_xyyyz_0_yzzzzzzz_0, g_xyyyz_0_zzzzzzzz_0, g_yyyz_0_xxxxxxxz_1, g_yyyz_0_xxxxxxyz_1, g_yyyz_0_xxxxxxz_1, g_yyyz_0_xxxxxxzz_1, g_yyyz_0_xxxxxyyz_1, g_yyyz_0_xxxxxyz_1, g_yyyz_0_xxxxxyzz_1, g_yyyz_0_xxxxxzz_1, g_yyyz_0_xxxxxzzz_1, g_yyyz_0_xxxxyyyz_1, g_yyyz_0_xxxxyyz_1, g_yyyz_0_xxxxyyzz_1, g_yyyz_0_xxxxyzz_1, g_yyyz_0_xxxxyzzz_1, g_yyyz_0_xxxxzzz_1, g_yyyz_0_xxxxzzzz_1, g_yyyz_0_xxxyyyyz_1, g_yyyz_0_xxxyyyz_1, g_yyyz_0_xxxyyyzz_1, g_yyyz_0_xxxyyzz_1, g_yyyz_0_xxxyyzzz_1, g_yyyz_0_xxxyzzz_1, g_yyyz_0_xxxyzzzz_1, g_yyyz_0_xxxzzzz_1, g_yyyz_0_xxxzzzzz_1, g_yyyz_0_xxyyyyyz_1, g_yyyz_0_xxyyyyz_1, g_yyyz_0_xxyyyyzz_1, g_yyyz_0_xxyyyzz_1, g_yyyz_0_xxyyyzzz_1, g_yyyz_0_xxyyzzz_1, g_yyyz_0_xxyyzzzz_1, g_yyyz_0_xxyzzzz_1, g_yyyz_0_xxyzzzzz_1, g_yyyz_0_xxzzzzz_1, g_yyyz_0_xxzzzzzz_1, g_yyyz_0_xyyyyyyz_1, g_yyyz_0_xyyyyyz_1, g_yyyz_0_xyyyyyzz_1, g_yyyz_0_xyyyyzz_1, g_yyyz_0_xyyyyzzz_1, g_yyyz_0_xyyyzzz_1, g_yyyz_0_xyyyzzzz_1, g_yyyz_0_xyyzzzz_1, g_yyyz_0_xyyzzzzz_1, g_yyyz_0_xyzzzzz_1, g_yyyz_0_xyzzzzzz_1, g_yyyz_0_xzzzzzz_1, g_yyyz_0_xzzzzzzz_1, g_yyyz_0_yyyyyyyy_1, g_yyyz_0_yyyyyyyz_1, g_yyyz_0_yyyyyyz_1, g_yyyz_0_yyyyyyzz_1, g_yyyz_0_yyyyyzz_1, g_yyyz_0_yyyyyzzz_1, g_yyyz_0_yyyyzzz_1, g_yyyz_0_yyyyzzzz_1, g_yyyz_0_yyyzzzz_1, g_yyyz_0_yyyzzzzz_1, g_yyyz_0_yyzzzzz_1, g_yyyz_0_yyzzzzzz_1, g_yyyz_0_yzzzzzz_1, g_yyyz_0_yzzzzzzz_1, g_yyyz_0_zzzzzzz_1, g_yyyz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyz_0_xxxxxxxx_0[i] = g_xyyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_xyyyz_0_xxxxxxxy_0[i] = g_xyyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xyyyz_0_xxxxxxxz_0[i] = 7.0 * g_yyyz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxxxyy_0[i] = g_xyyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xyyyz_0_xxxxxxyz_0[i] = 6.0 * g_yyyz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxxxzz_0[i] = 6.0 * g_yyyz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxxyyy_0[i] = g_xyyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xyyyz_0_xxxxxyyz_0[i] = 5.0 * g_yyyz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxxyzz_0[i] = 5.0 * g_yyyz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxxzzz_0[i] = 5.0 * g_yyyz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxyyyy_0[i] = g_xyyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xyyyz_0_xxxxyyyz_0[i] = 4.0 * g_yyyz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxyyzz_0[i] = 4.0 * g_yyyz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxyzzz_0[i] = 4.0 * g_yyyz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxxzzzz_0[i] = 4.0 * g_yyyz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxyyyyy_0[i] = g_xyyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xyyyz_0_xxxyyyyz_0[i] = 3.0 * g_yyyz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyyz_0_xxxyyyzz_0[i] = 3.0 * g_yyyz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxyyzzz_0[i] = 3.0 * g_yyyz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxyzzzz_0[i] = 3.0 * g_yyyz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxxzzzzz_0[i] = 3.0 * g_yyyz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxyyyyyy_0[i] = g_xyyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xyyyz_0_xxyyyyyz_0[i] = 2.0 * g_yyyz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyyz_0_xxyyyyzz_0[i] = 2.0 * g_yyyz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyyz_0_xxyyyzzz_0[i] = 2.0 * g_yyyz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxyyzzzz_0[i] = 2.0 * g_yyyz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxyzzzzz_0[i] = 2.0 * g_yyyz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xxzzzzzz_0[i] = 2.0 * g_yyyz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xyyyyyyy_0[i] = g_xyyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xyyyz_0_xyyyyyyz_0[i] = g_yyyz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyyz_0_xyyyyyzz_0[i] = g_yyyz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyyz_0_xyyyyzzz_0[i] = g_yyyz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyyz_0_xyyyzzzz_0[i] = g_yyyz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xyyzzzzz_0[i] = g_yyyz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xyzzzzzz_0[i] = g_yyyz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_xzzzzzzz_0[i] = g_yyyz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_yyyyyyyy_0[i] = g_yyyz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyyz_0_yyyyyyyz_0[i] = g_yyyz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyyz_0_yyyyyyzz_0[i] = g_yyyz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyyz_0_yyyyyzzz_0[i] = g_yyyz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyyz_0_yyyyzzzz_0[i] = g_yyyz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyyz_0_yyyzzzzz_0[i] = g_yyyz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_yyzzzzzz_0[i] = g_yyyz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_yzzzzzzz_0[i] = g_yyyz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyyz_0_zzzzzzzz_0[i] = g_yyyz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 540-585 components of targeted buffer : HSL

    auto g_xyyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 540);

    auto g_xyyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 541);

    auto g_xyyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 542);

    auto g_xyyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 543);

    auto g_xyyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 544);

    auto g_xyyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 545);

    auto g_xyyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 546);

    auto g_xyyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 547);

    auto g_xyyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 548);

    auto g_xyyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 549);

    auto g_xyyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 550);

    auto g_xyyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 551);

    auto g_xyyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 552);

    auto g_xyyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 553);

    auto g_xyyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 554);

    auto g_xyyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 555);

    auto g_xyyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 556);

    auto g_xyyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 557);

    auto g_xyyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 558);

    auto g_xyyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 559);

    auto g_xyyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 560);

    auto g_xyyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 561);

    auto g_xyyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 562);

    auto g_xyyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 563);

    auto g_xyyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 564);

    auto g_xyyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 565);

    auto g_xyyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 566);

    auto g_xyyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 567);

    auto g_xyyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 568);

    auto g_xyyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 569);

    auto g_xyyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 570);

    auto g_xyyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 571);

    auto g_xyyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 572);

    auto g_xyyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 573);

    auto g_xyyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 574);

    auto g_xyyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 575);

    auto g_xyyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 576);

    auto g_xyyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 577);

    auto g_xyyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 578);

    auto g_xyyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 579);

    auto g_xyyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 580);

    auto g_xyyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 581);

    auto g_xyyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 582);

    auto g_xyyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 583);

    auto g_xyyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 584);

    #pragma omp simd aligned(g_xyyzz_0_xxxxxxxx_0, g_xyyzz_0_xxxxxxxy_0, g_xyyzz_0_xxxxxxxz_0, g_xyyzz_0_xxxxxxyy_0, g_xyyzz_0_xxxxxxyz_0, g_xyyzz_0_xxxxxxzz_0, g_xyyzz_0_xxxxxyyy_0, g_xyyzz_0_xxxxxyyz_0, g_xyyzz_0_xxxxxyzz_0, g_xyyzz_0_xxxxxzzz_0, g_xyyzz_0_xxxxyyyy_0, g_xyyzz_0_xxxxyyyz_0, g_xyyzz_0_xxxxyyzz_0, g_xyyzz_0_xxxxyzzz_0, g_xyyzz_0_xxxxzzzz_0, g_xyyzz_0_xxxyyyyy_0, g_xyyzz_0_xxxyyyyz_0, g_xyyzz_0_xxxyyyzz_0, g_xyyzz_0_xxxyyzzz_0, g_xyyzz_0_xxxyzzzz_0, g_xyyzz_0_xxxzzzzz_0, g_xyyzz_0_xxyyyyyy_0, g_xyyzz_0_xxyyyyyz_0, g_xyyzz_0_xxyyyyzz_0, g_xyyzz_0_xxyyyzzz_0, g_xyyzz_0_xxyyzzzz_0, g_xyyzz_0_xxyzzzzz_0, g_xyyzz_0_xxzzzzzz_0, g_xyyzz_0_xyyyyyyy_0, g_xyyzz_0_xyyyyyyz_0, g_xyyzz_0_xyyyyyzz_0, g_xyyzz_0_xyyyyzzz_0, g_xyyzz_0_xyyyzzzz_0, g_xyyzz_0_xyyzzzzz_0, g_xyyzz_0_xyzzzzzz_0, g_xyyzz_0_xzzzzzzz_0, g_xyyzz_0_yyyyyyyy_0, g_xyyzz_0_yyyyyyyz_0, g_xyyzz_0_yyyyyyzz_0, g_xyyzz_0_yyyyyzzz_0, g_xyyzz_0_yyyyzzzz_0, g_xyyzz_0_yyyzzzzz_0, g_xyyzz_0_yyzzzzzz_0, g_xyyzz_0_yzzzzzzz_0, g_xyyzz_0_zzzzzzzz_0, g_yyzz_0_xxxxxxx_1, g_yyzz_0_xxxxxxxx_1, g_yyzz_0_xxxxxxxy_1, g_yyzz_0_xxxxxxxz_1, g_yyzz_0_xxxxxxy_1, g_yyzz_0_xxxxxxyy_1, g_yyzz_0_xxxxxxyz_1, g_yyzz_0_xxxxxxz_1, g_yyzz_0_xxxxxxzz_1, g_yyzz_0_xxxxxyy_1, g_yyzz_0_xxxxxyyy_1, g_yyzz_0_xxxxxyyz_1, g_yyzz_0_xxxxxyz_1, g_yyzz_0_xxxxxyzz_1, g_yyzz_0_xxxxxzz_1, g_yyzz_0_xxxxxzzz_1, g_yyzz_0_xxxxyyy_1, g_yyzz_0_xxxxyyyy_1, g_yyzz_0_xxxxyyyz_1, g_yyzz_0_xxxxyyz_1, g_yyzz_0_xxxxyyzz_1, g_yyzz_0_xxxxyzz_1, g_yyzz_0_xxxxyzzz_1, g_yyzz_0_xxxxzzz_1, g_yyzz_0_xxxxzzzz_1, g_yyzz_0_xxxyyyy_1, g_yyzz_0_xxxyyyyy_1, g_yyzz_0_xxxyyyyz_1, g_yyzz_0_xxxyyyz_1, g_yyzz_0_xxxyyyzz_1, g_yyzz_0_xxxyyzz_1, g_yyzz_0_xxxyyzzz_1, g_yyzz_0_xxxyzzz_1, g_yyzz_0_xxxyzzzz_1, g_yyzz_0_xxxzzzz_1, g_yyzz_0_xxxzzzzz_1, g_yyzz_0_xxyyyyy_1, g_yyzz_0_xxyyyyyy_1, g_yyzz_0_xxyyyyyz_1, g_yyzz_0_xxyyyyz_1, g_yyzz_0_xxyyyyzz_1, g_yyzz_0_xxyyyzz_1, g_yyzz_0_xxyyyzzz_1, g_yyzz_0_xxyyzzz_1, g_yyzz_0_xxyyzzzz_1, g_yyzz_0_xxyzzzz_1, g_yyzz_0_xxyzzzzz_1, g_yyzz_0_xxzzzzz_1, g_yyzz_0_xxzzzzzz_1, g_yyzz_0_xyyyyyy_1, g_yyzz_0_xyyyyyyy_1, g_yyzz_0_xyyyyyyz_1, g_yyzz_0_xyyyyyz_1, g_yyzz_0_xyyyyyzz_1, g_yyzz_0_xyyyyzz_1, g_yyzz_0_xyyyyzzz_1, g_yyzz_0_xyyyzzz_1, g_yyzz_0_xyyyzzzz_1, g_yyzz_0_xyyzzzz_1, g_yyzz_0_xyyzzzzz_1, g_yyzz_0_xyzzzzz_1, g_yyzz_0_xyzzzzzz_1, g_yyzz_0_xzzzzzz_1, g_yyzz_0_xzzzzzzz_1, g_yyzz_0_yyyyyyy_1, g_yyzz_0_yyyyyyyy_1, g_yyzz_0_yyyyyyyz_1, g_yyzz_0_yyyyyyz_1, g_yyzz_0_yyyyyyzz_1, g_yyzz_0_yyyyyzz_1, g_yyzz_0_yyyyyzzz_1, g_yyzz_0_yyyyzzz_1, g_yyzz_0_yyyyzzzz_1, g_yyzz_0_yyyzzzz_1, g_yyzz_0_yyyzzzzz_1, g_yyzz_0_yyzzzzz_1, g_yyzz_0_yyzzzzzz_1, g_yyzz_0_yzzzzzz_1, g_yyzz_0_yzzzzzzz_1, g_yyzz_0_zzzzzzz_1, g_yyzz_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzz_0_xxxxxxxx_0[i] = 8.0 * g_yyzz_0_xxxxxxx_1[i] * fi_acd_0 + g_yyzz_0_xxxxxxxx_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxxxy_0[i] = 7.0 * g_yyzz_0_xxxxxxy_1[i] * fi_acd_0 + g_yyzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxxxz_0[i] = 7.0 * g_yyzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxxyy_0[i] = 6.0 * g_yyzz_0_xxxxxyy_1[i] * fi_acd_0 + g_yyzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxxyz_0[i] = 6.0 * g_yyzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxxzz_0[i] = 6.0 * g_yyzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxyyy_0[i] = 5.0 * g_yyzz_0_xxxxyyy_1[i] * fi_acd_0 + g_yyzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxyyz_0[i] = 5.0 * g_yyzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxyzz_0[i] = 5.0 * g_yyzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxxzzz_0[i] = 5.0 * g_yyzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxyyyy_0[i] = 4.0 * g_yyzz_0_xxxyyyy_1[i] * fi_acd_0 + g_yyzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyyzz_0_xxxxyyyz_0[i] = 4.0 * g_yyzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxyyzz_0[i] = 4.0 * g_yyzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxyzzz_0[i] = 4.0 * g_yyzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxxzzzz_0[i] = 4.0 * g_yyzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxyyyyy_0[i] = 3.0 * g_yyzz_0_xxyyyyy_1[i] * fi_acd_0 + g_yyzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyyzz_0_xxxyyyyz_0[i] = 3.0 * g_yyzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyzz_0_xxxyyyzz_0[i] = 3.0 * g_yyzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxyyzzz_0[i] = 3.0 * g_yyzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxyzzzz_0[i] = 3.0 * g_yyzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxxzzzzz_0[i] = 3.0 * g_yyzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxyyyyyy_0[i] = 2.0 * g_yyzz_0_xyyyyyy_1[i] * fi_acd_0 + g_yyzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyyzz_0_xxyyyyyz_0[i] = 2.0 * g_yyzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyzz_0_xxyyyyzz_0[i] = 2.0 * g_yyzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyzz_0_xxyyyzzz_0[i] = 2.0 * g_yyzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxyyzzzz_0[i] = 2.0 * g_yyzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxyzzzzz_0[i] = 2.0 * g_yyzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xxzzzzzz_0[i] = 2.0 * g_yyzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xyyyyyyy_0[i] = g_yyzz_0_yyyyyyy_1[i] * fi_acd_0 + g_yyzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyyzz_0_xyyyyyyz_0[i] = g_yyzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyzz_0_xyyyyyzz_0[i] = g_yyzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyzz_0_xyyyyzzz_0[i] = g_yyzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyzz_0_xyyyzzzz_0[i] = g_yyzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xyyzzzzz_0[i] = g_yyzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xyzzzzzz_0[i] = g_yyzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_xzzzzzzz_0[i] = g_yyzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_yyyyyyyy_0[i] = g_yyzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyzz_0_yyyyyyyz_0[i] = g_yyzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyzz_0_yyyyyyzz_0[i] = g_yyzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyzz_0_yyyyyzzz_0[i] = g_yyzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyzz_0_yyyyzzzz_0[i] = g_yyzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyzz_0_yyyzzzzz_0[i] = g_yyzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_yyzzzzzz_0[i] = g_yyzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_yzzzzzzz_0[i] = g_yyzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyzz_0_zzzzzzzz_0[i] = g_yyzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 585-630 components of targeted buffer : HSL

    auto g_xyzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 585);

    auto g_xyzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 586);

    auto g_xyzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 587);

    auto g_xyzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 588);

    auto g_xyzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 589);

    auto g_xyzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 590);

    auto g_xyzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 591);

    auto g_xyzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 592);

    auto g_xyzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 593);

    auto g_xyzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 594);

    auto g_xyzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 595);

    auto g_xyzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 596);

    auto g_xyzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 597);

    auto g_xyzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 598);

    auto g_xyzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 599);

    auto g_xyzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 600);

    auto g_xyzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 601);

    auto g_xyzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 602);

    auto g_xyzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 603);

    auto g_xyzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 604);

    auto g_xyzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 605);

    auto g_xyzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 606);

    auto g_xyzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 607);

    auto g_xyzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 608);

    auto g_xyzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 609);

    auto g_xyzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 610);

    auto g_xyzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 611);

    auto g_xyzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 612);

    auto g_xyzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 613);

    auto g_xyzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 614);

    auto g_xyzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 615);

    auto g_xyzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 616);

    auto g_xyzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 617);

    auto g_xyzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 618);

    auto g_xyzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 619);

    auto g_xyzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 620);

    auto g_xyzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 621);

    auto g_xyzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 622);

    auto g_xyzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 623);

    auto g_xyzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 624);

    auto g_xyzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 625);

    auto g_xyzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 626);

    auto g_xyzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 627);

    auto g_xyzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 628);

    auto g_xyzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 629);

    #pragma omp simd aligned(g_xyzzz_0_xxxxxxxx_0, g_xyzzz_0_xxxxxxxy_0, g_xyzzz_0_xxxxxxxz_0, g_xyzzz_0_xxxxxxyy_0, g_xyzzz_0_xxxxxxyz_0, g_xyzzz_0_xxxxxxzz_0, g_xyzzz_0_xxxxxyyy_0, g_xyzzz_0_xxxxxyyz_0, g_xyzzz_0_xxxxxyzz_0, g_xyzzz_0_xxxxxzzz_0, g_xyzzz_0_xxxxyyyy_0, g_xyzzz_0_xxxxyyyz_0, g_xyzzz_0_xxxxyyzz_0, g_xyzzz_0_xxxxyzzz_0, g_xyzzz_0_xxxxzzzz_0, g_xyzzz_0_xxxyyyyy_0, g_xyzzz_0_xxxyyyyz_0, g_xyzzz_0_xxxyyyzz_0, g_xyzzz_0_xxxyyzzz_0, g_xyzzz_0_xxxyzzzz_0, g_xyzzz_0_xxxzzzzz_0, g_xyzzz_0_xxyyyyyy_0, g_xyzzz_0_xxyyyyyz_0, g_xyzzz_0_xxyyyyzz_0, g_xyzzz_0_xxyyyzzz_0, g_xyzzz_0_xxyyzzzz_0, g_xyzzz_0_xxyzzzzz_0, g_xyzzz_0_xxzzzzzz_0, g_xyzzz_0_xyyyyyyy_0, g_xyzzz_0_xyyyyyyz_0, g_xyzzz_0_xyyyyyzz_0, g_xyzzz_0_xyyyyzzz_0, g_xyzzz_0_xyyyzzzz_0, g_xyzzz_0_xyyzzzzz_0, g_xyzzz_0_xyzzzzzz_0, g_xyzzz_0_xzzzzzzz_0, g_xyzzz_0_yyyyyyyy_0, g_xyzzz_0_yyyyyyyz_0, g_xyzzz_0_yyyyyyzz_0, g_xyzzz_0_yyyyyzzz_0, g_xyzzz_0_yyyyzzzz_0, g_xyzzz_0_yyyzzzzz_0, g_xyzzz_0_yyzzzzzz_0, g_xyzzz_0_yzzzzzzz_0, g_xyzzz_0_zzzzzzzz_0, g_xzzz_0_xxxxxxxx_1, g_xzzz_0_xxxxxxxz_1, g_xzzz_0_xxxxxxzz_1, g_xzzz_0_xxxxxzzz_1, g_xzzz_0_xxxxzzzz_1, g_xzzz_0_xxxzzzzz_1, g_xzzz_0_xxzzzzzz_1, g_xzzz_0_xzzzzzzz_1, g_yzzz_0_xxxxxxxy_1, g_yzzz_0_xxxxxxy_1, g_yzzz_0_xxxxxxyy_1, g_yzzz_0_xxxxxxyz_1, g_yzzz_0_xxxxxyy_1, g_yzzz_0_xxxxxyyy_1, g_yzzz_0_xxxxxyyz_1, g_yzzz_0_xxxxxyz_1, g_yzzz_0_xxxxxyzz_1, g_yzzz_0_xxxxyyy_1, g_yzzz_0_xxxxyyyy_1, g_yzzz_0_xxxxyyyz_1, g_yzzz_0_xxxxyyz_1, g_yzzz_0_xxxxyyzz_1, g_yzzz_0_xxxxyzz_1, g_yzzz_0_xxxxyzzz_1, g_yzzz_0_xxxyyyy_1, g_yzzz_0_xxxyyyyy_1, g_yzzz_0_xxxyyyyz_1, g_yzzz_0_xxxyyyz_1, g_yzzz_0_xxxyyyzz_1, g_yzzz_0_xxxyyzz_1, g_yzzz_0_xxxyyzzz_1, g_yzzz_0_xxxyzzz_1, g_yzzz_0_xxxyzzzz_1, g_yzzz_0_xxyyyyy_1, g_yzzz_0_xxyyyyyy_1, g_yzzz_0_xxyyyyyz_1, g_yzzz_0_xxyyyyz_1, g_yzzz_0_xxyyyyzz_1, g_yzzz_0_xxyyyzz_1, g_yzzz_0_xxyyyzzz_1, g_yzzz_0_xxyyzzz_1, g_yzzz_0_xxyyzzzz_1, g_yzzz_0_xxyzzzz_1, g_yzzz_0_xxyzzzzz_1, g_yzzz_0_xyyyyyy_1, g_yzzz_0_xyyyyyyy_1, g_yzzz_0_xyyyyyyz_1, g_yzzz_0_xyyyyyz_1, g_yzzz_0_xyyyyyzz_1, g_yzzz_0_xyyyyzz_1, g_yzzz_0_xyyyyzzz_1, g_yzzz_0_xyyyzzz_1, g_yzzz_0_xyyyzzzz_1, g_yzzz_0_xyyzzzz_1, g_yzzz_0_xyyzzzzz_1, g_yzzz_0_xyzzzzz_1, g_yzzz_0_xyzzzzzz_1, g_yzzz_0_yyyyyyy_1, g_yzzz_0_yyyyyyyy_1, g_yzzz_0_yyyyyyyz_1, g_yzzz_0_yyyyyyz_1, g_yzzz_0_yyyyyyzz_1, g_yzzz_0_yyyyyzz_1, g_yzzz_0_yyyyyzzz_1, g_yzzz_0_yyyyzzz_1, g_yzzz_0_yyyyzzzz_1, g_yzzz_0_yyyzzzz_1, g_yzzz_0_yyyzzzzz_1, g_yzzz_0_yyzzzzz_1, g_yzzz_0_yyzzzzzz_1, g_yzzz_0_yzzzzzz_1, g_yzzz_0_yzzzzzzz_1, g_yzzz_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzz_0_xxxxxxxx_0[i] = g_xzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xyzzz_0_xxxxxxxy_0[i] = 7.0 * g_yzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_yzzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyzzz_0_xxxxxxxz_0[i] = g_xzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xyzzz_0_xxxxxxyy_0[i] = 6.0 * g_yzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_yzzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyzzz_0_xxxxxxyz_0[i] = 6.0 * g_yzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyzzz_0_xxxxxxzz_0[i] = g_xzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xyzzz_0_xxxxxyyy_0[i] = 5.0 * g_yzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_yzzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyzzz_0_xxxxxyyz_0[i] = 5.0 * g_yzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyzzz_0_xxxxxyzz_0[i] = 5.0 * g_yzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyzzz_0_xxxxxzzz_0[i] = g_xzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xyzzz_0_xxxxyyyy_0[i] = 4.0 * g_yzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_yzzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyzzz_0_xxxxyyyz_0[i] = 4.0 * g_yzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyzzz_0_xxxxyyzz_0[i] = 4.0 * g_yzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyzzz_0_xxxxyzzz_0[i] = 4.0 * g_yzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyzzz_0_xxxxzzzz_0[i] = g_xzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xyzzz_0_xxxyyyyy_0[i] = 3.0 * g_yzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_yzzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyzzz_0_xxxyyyyz_0[i] = 3.0 * g_yzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyzzz_0_xxxyyyzz_0[i] = 3.0 * g_yzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyzzz_0_xxxyyzzz_0[i] = 3.0 * g_yzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyzzz_0_xxxyzzzz_0[i] = 3.0 * g_yzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyzzz_0_xxxzzzzz_0[i] = g_xzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xyzzz_0_xxyyyyyy_0[i] = 2.0 * g_yzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_yzzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyzzz_0_xxyyyyyz_0[i] = 2.0 * g_yzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyzzz_0_xxyyyyzz_0[i] = 2.0 * g_yzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyzzz_0_xxyyyzzz_0[i] = 2.0 * g_yzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyzzz_0_xxyyzzzz_0[i] = 2.0 * g_yzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyzzz_0_xxyzzzzz_0[i] = 2.0 * g_yzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyzzz_0_xxzzzzzz_0[i] = g_xzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xyzzz_0_xyyyyyyy_0[i] = g_yzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_yzzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyzzz_0_xyyyyyyz_0[i] = g_yzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyzzz_0_xyyyyyzz_0[i] = g_yzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyzzz_0_xyyyyzzz_0[i] = g_yzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyzzz_0_xyyyzzzz_0[i] = g_yzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyzzz_0_xyyzzzzz_0[i] = g_yzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyzzz_0_xyzzzzzz_0[i] = g_yzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyzzz_0_xzzzzzzz_0[i] = g_xzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xyzzz_0_yyyyyyyy_0[i] = g_yzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyzzz_0_yyyyyyyz_0[i] = g_yzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyzzz_0_yyyyyyzz_0[i] = g_yzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyzzz_0_yyyyyzzz_0[i] = g_yzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyzzz_0_yyyyzzzz_0[i] = g_yzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyzzz_0_yyyzzzzz_0[i] = g_yzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyzzz_0_yyzzzzzz_0[i] = g_yzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyzzz_0_yzzzzzzz_0[i] = g_yzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyzzz_0_zzzzzzzz_0[i] = g_yzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 630-675 components of targeted buffer : HSL

    auto g_xzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 630);

    auto g_xzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 631);

    auto g_xzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 632);

    auto g_xzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 633);

    auto g_xzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 634);

    auto g_xzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 635);

    auto g_xzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 636);

    auto g_xzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 637);

    auto g_xzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 638);

    auto g_xzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 639);

    auto g_xzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 640);

    auto g_xzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 641);

    auto g_xzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 642);

    auto g_xzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 643);

    auto g_xzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 644);

    auto g_xzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 645);

    auto g_xzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 646);

    auto g_xzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 647);

    auto g_xzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 648);

    auto g_xzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 649);

    auto g_xzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 650);

    auto g_xzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 651);

    auto g_xzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 652);

    auto g_xzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 653);

    auto g_xzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 654);

    auto g_xzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 655);

    auto g_xzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 656);

    auto g_xzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 657);

    auto g_xzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 658);

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

    #pragma omp simd aligned(g_xzzzz_0_xxxxxxxx_0, g_xzzzz_0_xxxxxxxy_0, g_xzzzz_0_xxxxxxxz_0, g_xzzzz_0_xxxxxxyy_0, g_xzzzz_0_xxxxxxyz_0, g_xzzzz_0_xxxxxxzz_0, g_xzzzz_0_xxxxxyyy_0, g_xzzzz_0_xxxxxyyz_0, g_xzzzz_0_xxxxxyzz_0, g_xzzzz_0_xxxxxzzz_0, g_xzzzz_0_xxxxyyyy_0, g_xzzzz_0_xxxxyyyz_0, g_xzzzz_0_xxxxyyzz_0, g_xzzzz_0_xxxxyzzz_0, g_xzzzz_0_xxxxzzzz_0, g_xzzzz_0_xxxyyyyy_0, g_xzzzz_0_xxxyyyyz_0, g_xzzzz_0_xxxyyyzz_0, g_xzzzz_0_xxxyyzzz_0, g_xzzzz_0_xxxyzzzz_0, g_xzzzz_0_xxxzzzzz_0, g_xzzzz_0_xxyyyyyy_0, g_xzzzz_0_xxyyyyyz_0, g_xzzzz_0_xxyyyyzz_0, g_xzzzz_0_xxyyyzzz_0, g_xzzzz_0_xxyyzzzz_0, g_xzzzz_0_xxyzzzzz_0, g_xzzzz_0_xxzzzzzz_0, g_xzzzz_0_xyyyyyyy_0, g_xzzzz_0_xyyyyyyz_0, g_xzzzz_0_xyyyyyzz_0, g_xzzzz_0_xyyyyzzz_0, g_xzzzz_0_xyyyzzzz_0, g_xzzzz_0_xyyzzzzz_0, g_xzzzz_0_xyzzzzzz_0, g_xzzzz_0_xzzzzzzz_0, g_xzzzz_0_yyyyyyyy_0, g_xzzzz_0_yyyyyyyz_0, g_xzzzz_0_yyyyyyzz_0, g_xzzzz_0_yyyyyzzz_0, g_xzzzz_0_yyyyzzzz_0, g_xzzzz_0_yyyzzzzz_0, g_xzzzz_0_yyzzzzzz_0, g_xzzzz_0_yzzzzzzz_0, g_xzzzz_0_zzzzzzzz_0, g_zzzz_0_xxxxxxx_1, g_zzzz_0_xxxxxxxx_1, g_zzzz_0_xxxxxxxy_1, g_zzzz_0_xxxxxxxz_1, g_zzzz_0_xxxxxxy_1, g_zzzz_0_xxxxxxyy_1, g_zzzz_0_xxxxxxyz_1, g_zzzz_0_xxxxxxz_1, g_zzzz_0_xxxxxxzz_1, g_zzzz_0_xxxxxyy_1, g_zzzz_0_xxxxxyyy_1, g_zzzz_0_xxxxxyyz_1, g_zzzz_0_xxxxxyz_1, g_zzzz_0_xxxxxyzz_1, g_zzzz_0_xxxxxzz_1, g_zzzz_0_xxxxxzzz_1, g_zzzz_0_xxxxyyy_1, g_zzzz_0_xxxxyyyy_1, g_zzzz_0_xxxxyyyz_1, g_zzzz_0_xxxxyyz_1, g_zzzz_0_xxxxyyzz_1, g_zzzz_0_xxxxyzz_1, g_zzzz_0_xxxxyzzz_1, g_zzzz_0_xxxxzzz_1, g_zzzz_0_xxxxzzzz_1, g_zzzz_0_xxxyyyy_1, g_zzzz_0_xxxyyyyy_1, g_zzzz_0_xxxyyyyz_1, g_zzzz_0_xxxyyyz_1, g_zzzz_0_xxxyyyzz_1, g_zzzz_0_xxxyyzz_1, g_zzzz_0_xxxyyzzz_1, g_zzzz_0_xxxyzzz_1, g_zzzz_0_xxxyzzzz_1, g_zzzz_0_xxxzzzz_1, g_zzzz_0_xxxzzzzz_1, g_zzzz_0_xxyyyyy_1, g_zzzz_0_xxyyyyyy_1, g_zzzz_0_xxyyyyyz_1, g_zzzz_0_xxyyyyz_1, g_zzzz_0_xxyyyyzz_1, g_zzzz_0_xxyyyzz_1, g_zzzz_0_xxyyyzzz_1, g_zzzz_0_xxyyzzz_1, g_zzzz_0_xxyyzzzz_1, g_zzzz_0_xxyzzzz_1, g_zzzz_0_xxyzzzzz_1, g_zzzz_0_xxzzzzz_1, g_zzzz_0_xxzzzzzz_1, g_zzzz_0_xyyyyyy_1, g_zzzz_0_xyyyyyyy_1, g_zzzz_0_xyyyyyyz_1, g_zzzz_0_xyyyyyz_1, g_zzzz_0_xyyyyyzz_1, g_zzzz_0_xyyyyzz_1, g_zzzz_0_xyyyyzzz_1, g_zzzz_0_xyyyzzz_1, g_zzzz_0_xyyyzzzz_1, g_zzzz_0_xyyzzzz_1, g_zzzz_0_xyyzzzzz_1, g_zzzz_0_xyzzzzz_1, g_zzzz_0_xyzzzzzz_1, g_zzzz_0_xzzzzzz_1, g_zzzz_0_xzzzzzzz_1, g_zzzz_0_yyyyyyy_1, g_zzzz_0_yyyyyyyy_1, g_zzzz_0_yyyyyyyz_1, g_zzzz_0_yyyyyyz_1, g_zzzz_0_yyyyyyzz_1, g_zzzz_0_yyyyyzz_1, g_zzzz_0_yyyyyzzz_1, g_zzzz_0_yyyyzzz_1, g_zzzz_0_yyyyzzzz_1, g_zzzz_0_yyyzzzz_1, g_zzzz_0_yyyzzzzz_1, g_zzzz_0_yyzzzzz_1, g_zzzz_0_yyzzzzzz_1, g_zzzz_0_yzzzzzz_1, g_zzzz_0_yzzzzzzz_1, g_zzzz_0_zzzzzzz_1, g_zzzz_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzz_0_xxxxxxxx_0[i] = 8.0 * g_zzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxxx_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxxxy_0[i] = 7.0 * g_zzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxxxz_0[i] = 7.0 * g_zzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxxyy_0[i] = 6.0 * g_zzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxxyz_0[i] = 6.0 * g_zzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxxzz_0[i] = 6.0 * g_zzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxyyy_0[i] = 5.0 * g_zzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxyyz_0[i] = 5.0 * g_zzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxyzz_0[i] = 5.0 * g_zzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxxzzz_0[i] = 5.0 * g_zzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxyyyy_0[i] = 4.0 * g_zzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xzzzz_0_xxxxyyyz_0[i] = 4.0 * g_zzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxyyzz_0[i] = 4.0 * g_zzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxyzzz_0[i] = 4.0 * g_zzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxxzzzz_0[i] = 4.0 * g_zzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxyyyyy_0[i] = 3.0 * g_zzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xzzzz_0_xxxyyyyz_0[i] = 3.0 * g_zzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xzzzz_0_xxxyyyzz_0[i] = 3.0 * g_zzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxyyzzz_0[i] = 3.0 * g_zzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxyzzzz_0[i] = 3.0 * g_zzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxxzzzzz_0[i] = 3.0 * g_zzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxyyyyyy_0[i] = 2.0 * g_zzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xzzzz_0_xxyyyyyz_0[i] = 2.0 * g_zzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xzzzz_0_xxyyyyzz_0[i] = 2.0 * g_zzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xzzzz_0_xxyyyzzz_0[i] = 2.0 * g_zzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxyyzzzz_0[i] = 2.0 * g_zzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxyzzzzz_0[i] = 2.0 * g_zzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xxzzzzzz_0[i] = 2.0 * g_zzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xyyyyyyy_0[i] = g_zzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xzzzz_0_xyyyyyyz_0[i] = g_zzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xzzzz_0_xyyyyyzz_0[i] = g_zzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xzzzz_0_xyyyyzzz_0[i] = g_zzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xzzzz_0_xyyyzzzz_0[i] = g_zzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xyyzzzzz_0[i] = g_zzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xyzzzzzz_0[i] = g_zzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_xzzzzzzz_0[i] = g_zzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_yyyyyyyy_0[i] = g_zzzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xzzzz_0_yyyyyyyz_0[i] = g_zzzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xzzzz_0_yyyyyyzz_0[i] = g_zzzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xzzzz_0_yyyyyzzz_0[i] = g_zzzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xzzzz_0_yyyyzzzz_0[i] = g_zzzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xzzzz_0_yyyzzzzz_0[i] = g_zzzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_yyzzzzzz_0[i] = g_zzzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_yzzzzzzz_0[i] = g_zzzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xzzzz_0_zzzzzzzz_0[i] = g_zzzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 675-720 components of targeted buffer : HSL

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

    #pragma omp simd aligned(g_yyy_0_xxxxxxxx_0, g_yyy_0_xxxxxxxx_1, g_yyy_0_xxxxxxxy_0, g_yyy_0_xxxxxxxy_1, g_yyy_0_xxxxxxxz_0, g_yyy_0_xxxxxxxz_1, g_yyy_0_xxxxxxyy_0, g_yyy_0_xxxxxxyy_1, g_yyy_0_xxxxxxyz_0, g_yyy_0_xxxxxxyz_1, g_yyy_0_xxxxxxzz_0, g_yyy_0_xxxxxxzz_1, g_yyy_0_xxxxxyyy_0, g_yyy_0_xxxxxyyy_1, g_yyy_0_xxxxxyyz_0, g_yyy_0_xxxxxyyz_1, g_yyy_0_xxxxxyzz_0, g_yyy_0_xxxxxyzz_1, g_yyy_0_xxxxxzzz_0, g_yyy_0_xxxxxzzz_1, g_yyy_0_xxxxyyyy_0, g_yyy_0_xxxxyyyy_1, g_yyy_0_xxxxyyyz_0, g_yyy_0_xxxxyyyz_1, g_yyy_0_xxxxyyzz_0, g_yyy_0_xxxxyyzz_1, g_yyy_0_xxxxyzzz_0, g_yyy_0_xxxxyzzz_1, g_yyy_0_xxxxzzzz_0, g_yyy_0_xxxxzzzz_1, g_yyy_0_xxxyyyyy_0, g_yyy_0_xxxyyyyy_1, g_yyy_0_xxxyyyyz_0, g_yyy_0_xxxyyyyz_1, g_yyy_0_xxxyyyzz_0, g_yyy_0_xxxyyyzz_1, g_yyy_0_xxxyyzzz_0, g_yyy_0_xxxyyzzz_1, g_yyy_0_xxxyzzzz_0, g_yyy_0_xxxyzzzz_1, g_yyy_0_xxxzzzzz_0, g_yyy_0_xxxzzzzz_1, g_yyy_0_xxyyyyyy_0, g_yyy_0_xxyyyyyy_1, g_yyy_0_xxyyyyyz_0, g_yyy_0_xxyyyyyz_1, g_yyy_0_xxyyyyzz_0, g_yyy_0_xxyyyyzz_1, g_yyy_0_xxyyyzzz_0, g_yyy_0_xxyyyzzz_1, g_yyy_0_xxyyzzzz_0, g_yyy_0_xxyyzzzz_1, g_yyy_0_xxyzzzzz_0, g_yyy_0_xxyzzzzz_1, g_yyy_0_xxzzzzzz_0, g_yyy_0_xxzzzzzz_1, g_yyy_0_xyyyyyyy_0, g_yyy_0_xyyyyyyy_1, g_yyy_0_xyyyyyyz_0, g_yyy_0_xyyyyyyz_1, g_yyy_0_xyyyyyzz_0, g_yyy_0_xyyyyyzz_1, g_yyy_0_xyyyyzzz_0, g_yyy_0_xyyyyzzz_1, g_yyy_0_xyyyzzzz_0, g_yyy_0_xyyyzzzz_1, g_yyy_0_xyyzzzzz_0, g_yyy_0_xyyzzzzz_1, g_yyy_0_xyzzzzzz_0, g_yyy_0_xyzzzzzz_1, g_yyy_0_xzzzzzzz_0, g_yyy_0_xzzzzzzz_1, g_yyy_0_yyyyyyyy_0, g_yyy_0_yyyyyyyy_1, g_yyy_0_yyyyyyyz_0, g_yyy_0_yyyyyyyz_1, g_yyy_0_yyyyyyzz_0, g_yyy_0_yyyyyyzz_1, g_yyy_0_yyyyyzzz_0, g_yyy_0_yyyyyzzz_1, g_yyy_0_yyyyzzzz_0, g_yyy_0_yyyyzzzz_1, g_yyy_0_yyyzzzzz_0, g_yyy_0_yyyzzzzz_1, g_yyy_0_yyzzzzzz_0, g_yyy_0_yyzzzzzz_1, g_yyy_0_yzzzzzzz_0, g_yyy_0_yzzzzzzz_1, g_yyy_0_zzzzzzzz_0, g_yyy_0_zzzzzzzz_1, g_yyyy_0_xxxxxxx_1, g_yyyy_0_xxxxxxxx_1, g_yyyy_0_xxxxxxxy_1, g_yyyy_0_xxxxxxxz_1, g_yyyy_0_xxxxxxy_1, g_yyyy_0_xxxxxxyy_1, g_yyyy_0_xxxxxxyz_1, g_yyyy_0_xxxxxxz_1, g_yyyy_0_xxxxxxzz_1, g_yyyy_0_xxxxxyy_1, g_yyyy_0_xxxxxyyy_1, g_yyyy_0_xxxxxyyz_1, g_yyyy_0_xxxxxyz_1, g_yyyy_0_xxxxxyzz_1, g_yyyy_0_xxxxxzz_1, g_yyyy_0_xxxxxzzz_1, g_yyyy_0_xxxxyyy_1, g_yyyy_0_xxxxyyyy_1, g_yyyy_0_xxxxyyyz_1, g_yyyy_0_xxxxyyz_1, g_yyyy_0_xxxxyyzz_1, g_yyyy_0_xxxxyzz_1, g_yyyy_0_xxxxyzzz_1, g_yyyy_0_xxxxzzz_1, g_yyyy_0_xxxxzzzz_1, g_yyyy_0_xxxyyyy_1, g_yyyy_0_xxxyyyyy_1, g_yyyy_0_xxxyyyyz_1, g_yyyy_0_xxxyyyz_1, g_yyyy_0_xxxyyyzz_1, g_yyyy_0_xxxyyzz_1, g_yyyy_0_xxxyyzzz_1, g_yyyy_0_xxxyzzz_1, g_yyyy_0_xxxyzzzz_1, g_yyyy_0_xxxzzzz_1, g_yyyy_0_xxxzzzzz_1, g_yyyy_0_xxyyyyy_1, g_yyyy_0_xxyyyyyy_1, g_yyyy_0_xxyyyyyz_1, g_yyyy_0_xxyyyyz_1, g_yyyy_0_xxyyyyzz_1, g_yyyy_0_xxyyyzz_1, g_yyyy_0_xxyyyzzz_1, g_yyyy_0_xxyyzzz_1, g_yyyy_0_xxyyzzzz_1, g_yyyy_0_xxyzzzz_1, g_yyyy_0_xxyzzzzz_1, g_yyyy_0_xxzzzzz_1, g_yyyy_0_xxzzzzzz_1, g_yyyy_0_xyyyyyy_1, g_yyyy_0_xyyyyyyy_1, g_yyyy_0_xyyyyyyz_1, g_yyyy_0_xyyyyyz_1, g_yyyy_0_xyyyyyzz_1, g_yyyy_0_xyyyyzz_1, g_yyyy_0_xyyyyzzz_1, g_yyyy_0_xyyyzzz_1, g_yyyy_0_xyyyzzzz_1, g_yyyy_0_xyyzzzz_1, g_yyyy_0_xyyzzzzz_1, g_yyyy_0_xyzzzzz_1, g_yyyy_0_xyzzzzzz_1, g_yyyy_0_xzzzzzz_1, g_yyyy_0_xzzzzzzz_1, g_yyyy_0_yyyyyyy_1, g_yyyy_0_yyyyyyyy_1, g_yyyy_0_yyyyyyyz_1, g_yyyy_0_yyyyyyz_1, g_yyyy_0_yyyyyyzz_1, g_yyyy_0_yyyyyzz_1, g_yyyy_0_yyyyyzzz_1, g_yyyy_0_yyyyzzz_1, g_yyyy_0_yyyyzzzz_1, g_yyyy_0_yyyzzzz_1, g_yyyy_0_yyyzzzzz_1, g_yyyy_0_yyzzzzz_1, g_yyyy_0_yyzzzzzz_1, g_yyyy_0_yzzzzzz_1, g_yyyy_0_yzzzzzzz_1, g_yyyy_0_zzzzzzz_1, g_yyyy_0_zzzzzzzz_1, g_yyyyy_0_xxxxxxxx_0, g_yyyyy_0_xxxxxxxy_0, g_yyyyy_0_xxxxxxxz_0, g_yyyyy_0_xxxxxxyy_0, g_yyyyy_0_xxxxxxyz_0, g_yyyyy_0_xxxxxxzz_0, g_yyyyy_0_xxxxxyyy_0, g_yyyyy_0_xxxxxyyz_0, g_yyyyy_0_xxxxxyzz_0, g_yyyyy_0_xxxxxzzz_0, g_yyyyy_0_xxxxyyyy_0, g_yyyyy_0_xxxxyyyz_0, g_yyyyy_0_xxxxyyzz_0, g_yyyyy_0_xxxxyzzz_0, g_yyyyy_0_xxxxzzzz_0, g_yyyyy_0_xxxyyyyy_0, g_yyyyy_0_xxxyyyyz_0, g_yyyyy_0_xxxyyyzz_0, g_yyyyy_0_xxxyyzzz_0, g_yyyyy_0_xxxyzzzz_0, g_yyyyy_0_xxxzzzzz_0, g_yyyyy_0_xxyyyyyy_0, g_yyyyy_0_xxyyyyyz_0, g_yyyyy_0_xxyyyyzz_0, g_yyyyy_0_xxyyyzzz_0, g_yyyyy_0_xxyyzzzz_0, g_yyyyy_0_xxyzzzzz_0, g_yyyyy_0_xxzzzzzz_0, g_yyyyy_0_xyyyyyyy_0, g_yyyyy_0_xyyyyyyz_0, g_yyyyy_0_xyyyyyzz_0, g_yyyyy_0_xyyyyzzz_0, g_yyyyy_0_xyyyzzzz_0, g_yyyyy_0_xyyzzzzz_0, g_yyyyy_0_xyzzzzzz_0, g_yyyyy_0_xzzzzzzz_0, g_yyyyy_0_yyyyyyyy_0, g_yyyyy_0_yyyyyyyz_0, g_yyyyy_0_yyyyyyzz_0, g_yyyyy_0_yyyyyzzz_0, g_yyyyy_0_yyyyzzzz_0, g_yyyyy_0_yyyzzzzz_0, g_yyyyy_0_yyzzzzzz_0, g_yyyyy_0_yzzzzzzz_0, g_yyyyy_0_zzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyy_0_xxxxxxxx_0[i] = 4.0 * g_yyy_0_xxxxxxxx_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxxxx_1[i] * fz_be_0 + g_yyyy_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxxxy_0[i] = 4.0 * g_yyy_0_xxxxxxxy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxxxy_1[i] * fz_be_0 + g_yyyy_0_xxxxxxx_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxxy_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxxxz_0[i] = 4.0 * g_yyy_0_xxxxxxxz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxxxz_1[i] * fz_be_0 + g_yyyy_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxxyy_0[i] = 4.0 * g_yyy_0_xxxxxxyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxxyy_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxyy_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxxyz_0[i] = 4.0 * g_yyy_0_xxxxxxyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxxyz_1[i] * fz_be_0 + g_yyyy_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxxzz_0[i] = 4.0 * g_yyy_0_xxxxxxzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxxzz_1[i] * fz_be_0 + g_yyyy_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxyyy_0[i] = 4.0 * g_yyy_0_xxxxxyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxyyy_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyyy_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxyyz_0[i] = 4.0 * g_yyy_0_xxxxxyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxyzz_0[i] = 4.0 * g_yyy_0_xxxxxyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxyzz_1[i] * fz_be_0 + g_yyyy_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxxzzz_0[i] = 4.0 * g_yyy_0_xxxxxzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxxzzz_1[i] * fz_be_0 + g_yyyy_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxyyyy_0[i] = 4.0 * g_yyy_0_xxxxyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_yyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyyy_1[i] * wa_y[i];

        g_yyyyy_0_xxxxyyyz_0[i] = 4.0 * g_yyy_0_xxxxyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxyyzz_0[i] = 4.0 * g_yyy_0_xxxxyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxyzzz_0[i] = 4.0 * g_yyy_0_xxxxyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxyzzz_1[i] * fz_be_0 + g_yyyy_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxxzzzz_0[i] = 4.0 * g_yyy_0_xxxxzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxxzzzz_1[i] * fz_be_0 + g_yyyy_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxyyyyy_0[i] = 4.0 * g_yyy_0_xxxyyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyyyyy_1[i] * fz_be_0 + 5.0 * g_yyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyyy_1[i] * wa_y[i];

        g_yyyyy_0_xxxyyyyz_0[i] = 4.0 * g_yyy_0_xxxyyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyyyy_0_xxxyyyzz_0[i] = 4.0 * g_yyy_0_xxxyyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxyyzzz_0[i] = 4.0 * g_yyy_0_xxxyyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxyzzzz_0[i] = 4.0 * g_yyy_0_xxxyzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxyzzzz_1[i] * fz_be_0 + g_yyyy_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxxzzzzz_0[i] = 4.0 * g_yyy_0_xxxzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxxzzzzz_1[i] * fz_be_0 + g_yyyy_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxyyyyyy_0[i] = 4.0 * g_yyy_0_xxyyyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyyyyy_1[i] * fz_be_0 + 6.0 * g_yyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyyy_1[i] * wa_y[i];

        g_yyyyy_0_xxyyyyyz_0[i] = 4.0 * g_yyy_0_xxyyyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyyyy_0_xxyyyyzz_0[i] = 4.0 * g_yyy_0_xxyyyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyyyy_0_xxyyyzzz_0[i] = 4.0 * g_yyy_0_xxyyyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxyyzzzz_0[i] = 4.0 * g_yyy_0_xxyyzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxyzzzzz_0[i] = 4.0 * g_yyy_0_xxyzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxyzzzzz_1[i] * fz_be_0 + g_yyyy_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xxzzzzzz_0[i] = 4.0 * g_yyy_0_xxzzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xxzzzzzz_1[i] * fz_be_0 + g_yyyy_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xyyyyyyy_0[i] = 4.0 * g_yyy_0_xyyyyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyyyyy_1[i] * fz_be_0 + 7.0 * g_yyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyyy_1[i] * wa_y[i];

        g_yyyyy_0_xyyyyyyz_0[i] = 4.0 * g_yyy_0_xyyyyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyyyy_0_xyyyyyzz_0[i] = 4.0 * g_yyy_0_xyyyyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyyyy_0_xyyyyzzz_0[i] = 4.0 * g_yyy_0_xyyyyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyyyy_0_xyyyzzzz_0[i] = 4.0 * g_yyy_0_xyyyzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xyyzzzzz_0[i] = 4.0 * g_yyy_0_xyyzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xyzzzzzz_0[i] = 4.0 * g_yyy_0_xyzzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xyzzzzzz_1[i] * fz_be_0 + g_yyyy_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_xzzzzzzz_0[i] = 4.0 * g_yyy_0_xzzzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_xzzzzzzz_1[i] * fz_be_0 + g_yyyy_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_yyyyyyyy_0[i] = 4.0 * g_yyy_0_yyyyyyyy_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyyyyy_1[i] * fz_be_0 + 8.0 * g_yyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_yyyy_0_yyyyyyyy_1[i] * wa_y[i];

        g_yyyyy_0_yyyyyyyz_0[i] = 4.0 * g_yyy_0_yyyyyyyz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyy_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyyyy_0_yyyyyyzz_0[i] = 4.0 * g_yyy_0_yyyyyyzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyy_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyyyy_0_yyyyyzzz_0[i] = 4.0 * g_yyy_0_yyyyyzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyy_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyyyy_0_yyyyzzzz_0[i] = 4.0 * g_yyy_0_yyyyzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyy_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyyyy_0_yyyzzzzz_0[i] = 4.0 * g_yyy_0_yyyzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyy_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_yyzzzzzz_0[i] = 4.0 * g_yyy_0_yyzzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyy_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_yzzzzzzz_0[i] = 4.0 * g_yyy_0_yzzzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_yzzzzzzz_1[i] * fz_be_0 + g_yyyy_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyy_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyyyy_0_zzzzzzzz_0[i] = 4.0 * g_yyy_0_zzzzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_0_zzzzzzzz_1[i] * fz_be_0 + g_yyyy_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 720-765 components of targeted buffer : HSL

    auto g_yyyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 720);

    auto g_yyyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 721);

    auto g_yyyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 722);

    auto g_yyyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 723);

    auto g_yyyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 724);

    auto g_yyyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 725);

    auto g_yyyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 726);

    auto g_yyyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 727);

    auto g_yyyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 728);

    auto g_yyyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 729);

    auto g_yyyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 730);

    auto g_yyyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 731);

    auto g_yyyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 732);

    auto g_yyyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 733);

    auto g_yyyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 734);

    auto g_yyyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 735);

    auto g_yyyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 736);

    auto g_yyyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 737);

    auto g_yyyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 738);

    auto g_yyyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 739);

    auto g_yyyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 740);

    auto g_yyyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 741);

    auto g_yyyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 742);

    auto g_yyyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 743);

    auto g_yyyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 744);

    auto g_yyyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 745);

    auto g_yyyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 746);

    auto g_yyyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 747);

    auto g_yyyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 748);

    auto g_yyyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 749);

    auto g_yyyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 750);

    auto g_yyyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 751);

    auto g_yyyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 752);

    auto g_yyyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 753);

    auto g_yyyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 754);

    auto g_yyyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 755);

    auto g_yyyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 756);

    auto g_yyyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 757);

    auto g_yyyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 758);

    auto g_yyyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 759);

    auto g_yyyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 760);

    auto g_yyyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 761);

    auto g_yyyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 762);

    auto g_yyyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 763);

    auto g_yyyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 764);

    #pragma omp simd aligned(g_yyyy_0_xxxxxxx_1, g_yyyy_0_xxxxxxxx_1, g_yyyy_0_xxxxxxxy_1, g_yyyy_0_xxxxxxxz_1, g_yyyy_0_xxxxxxy_1, g_yyyy_0_xxxxxxyy_1, g_yyyy_0_xxxxxxyz_1, g_yyyy_0_xxxxxxz_1, g_yyyy_0_xxxxxxzz_1, g_yyyy_0_xxxxxyy_1, g_yyyy_0_xxxxxyyy_1, g_yyyy_0_xxxxxyyz_1, g_yyyy_0_xxxxxyz_1, g_yyyy_0_xxxxxyzz_1, g_yyyy_0_xxxxxzz_1, g_yyyy_0_xxxxxzzz_1, g_yyyy_0_xxxxyyy_1, g_yyyy_0_xxxxyyyy_1, g_yyyy_0_xxxxyyyz_1, g_yyyy_0_xxxxyyz_1, g_yyyy_0_xxxxyyzz_1, g_yyyy_0_xxxxyzz_1, g_yyyy_0_xxxxyzzz_1, g_yyyy_0_xxxxzzz_1, g_yyyy_0_xxxxzzzz_1, g_yyyy_0_xxxyyyy_1, g_yyyy_0_xxxyyyyy_1, g_yyyy_0_xxxyyyyz_1, g_yyyy_0_xxxyyyz_1, g_yyyy_0_xxxyyyzz_1, g_yyyy_0_xxxyyzz_1, g_yyyy_0_xxxyyzzz_1, g_yyyy_0_xxxyzzz_1, g_yyyy_0_xxxyzzzz_1, g_yyyy_0_xxxzzzz_1, g_yyyy_0_xxxzzzzz_1, g_yyyy_0_xxyyyyy_1, g_yyyy_0_xxyyyyyy_1, g_yyyy_0_xxyyyyyz_1, g_yyyy_0_xxyyyyz_1, g_yyyy_0_xxyyyyzz_1, g_yyyy_0_xxyyyzz_1, g_yyyy_0_xxyyyzzz_1, g_yyyy_0_xxyyzzz_1, g_yyyy_0_xxyyzzzz_1, g_yyyy_0_xxyzzzz_1, g_yyyy_0_xxyzzzzz_1, g_yyyy_0_xxzzzzz_1, g_yyyy_0_xxzzzzzz_1, g_yyyy_0_xyyyyyy_1, g_yyyy_0_xyyyyyyy_1, g_yyyy_0_xyyyyyyz_1, g_yyyy_0_xyyyyyz_1, g_yyyy_0_xyyyyyzz_1, g_yyyy_0_xyyyyzz_1, g_yyyy_0_xyyyyzzz_1, g_yyyy_0_xyyyzzz_1, g_yyyy_0_xyyyzzzz_1, g_yyyy_0_xyyzzzz_1, g_yyyy_0_xyyzzzzz_1, g_yyyy_0_xyzzzzz_1, g_yyyy_0_xyzzzzzz_1, g_yyyy_0_xzzzzzz_1, g_yyyy_0_xzzzzzzz_1, g_yyyy_0_yyyyyyy_1, g_yyyy_0_yyyyyyyy_1, g_yyyy_0_yyyyyyyz_1, g_yyyy_0_yyyyyyz_1, g_yyyy_0_yyyyyyzz_1, g_yyyy_0_yyyyyzz_1, g_yyyy_0_yyyyyzzz_1, g_yyyy_0_yyyyzzz_1, g_yyyy_0_yyyyzzzz_1, g_yyyy_0_yyyzzzz_1, g_yyyy_0_yyyzzzzz_1, g_yyyy_0_yyzzzzz_1, g_yyyy_0_yyzzzzzz_1, g_yyyy_0_yzzzzzz_1, g_yyyy_0_yzzzzzzz_1, g_yyyy_0_zzzzzzz_1, g_yyyy_0_zzzzzzzz_1, g_yyyyz_0_xxxxxxxx_0, g_yyyyz_0_xxxxxxxy_0, g_yyyyz_0_xxxxxxxz_0, g_yyyyz_0_xxxxxxyy_0, g_yyyyz_0_xxxxxxyz_0, g_yyyyz_0_xxxxxxzz_0, g_yyyyz_0_xxxxxyyy_0, g_yyyyz_0_xxxxxyyz_0, g_yyyyz_0_xxxxxyzz_0, g_yyyyz_0_xxxxxzzz_0, g_yyyyz_0_xxxxyyyy_0, g_yyyyz_0_xxxxyyyz_0, g_yyyyz_0_xxxxyyzz_0, g_yyyyz_0_xxxxyzzz_0, g_yyyyz_0_xxxxzzzz_0, g_yyyyz_0_xxxyyyyy_0, g_yyyyz_0_xxxyyyyz_0, g_yyyyz_0_xxxyyyzz_0, g_yyyyz_0_xxxyyzzz_0, g_yyyyz_0_xxxyzzzz_0, g_yyyyz_0_xxxzzzzz_0, g_yyyyz_0_xxyyyyyy_0, g_yyyyz_0_xxyyyyyz_0, g_yyyyz_0_xxyyyyzz_0, g_yyyyz_0_xxyyyzzz_0, g_yyyyz_0_xxyyzzzz_0, g_yyyyz_0_xxyzzzzz_0, g_yyyyz_0_xxzzzzzz_0, g_yyyyz_0_xyyyyyyy_0, g_yyyyz_0_xyyyyyyz_0, g_yyyyz_0_xyyyyyzz_0, g_yyyyz_0_xyyyyzzz_0, g_yyyyz_0_xyyyzzzz_0, g_yyyyz_0_xyyzzzzz_0, g_yyyyz_0_xyzzzzzz_0, g_yyyyz_0_xzzzzzzz_0, g_yyyyz_0_yyyyyyyy_0, g_yyyyz_0_yyyyyyyz_0, g_yyyyz_0_yyyyyyzz_0, g_yyyyz_0_yyyyyzzz_0, g_yyyyz_0_yyyyzzzz_0, g_yyyyz_0_yyyzzzzz_0, g_yyyyz_0_yyzzzzzz_0, g_yyyyz_0_yzzzzzzz_0, g_yyyyz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyz_0_xxxxxxxx_0[i] = g_yyyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxxxy_0[i] = g_yyyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxxxz_0[i] = g_yyyy_0_xxxxxxx_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxxz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxxyy_0[i] = g_yyyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxxyz_0[i] = g_yyyy_0_xxxxxxy_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxyz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxxzz_0[i] = 2.0 * g_yyyy_0_xxxxxxz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxxzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxyyy_0[i] = g_yyyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxyyz_0[i] = g_yyyy_0_xxxxxyy_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyyz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxyzz_0[i] = 2.0 * g_yyyy_0_xxxxxyz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxyzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxxzzz_0[i] = 3.0 * g_yyyy_0_xxxxxzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxxzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxyyyy_0[i] = g_yyyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyyyz_0_xxxxyyyz_0[i] = g_yyyy_0_xxxxyyy_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyyz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxyyzz_0[i] = 2.0 * g_yyyy_0_xxxxyyz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyyzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxyzzz_0[i] = 3.0 * g_yyyy_0_xxxxyzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxyzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxxzzzz_0[i] = 4.0 * g_yyyy_0_xxxxzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxxzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxyyyyy_0[i] = g_yyyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyyyz_0_xxxyyyyz_0[i] = g_yyyy_0_xxxyyyy_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyyz_1[i] * wa_z[i];

        g_yyyyz_0_xxxyyyzz_0[i] = 2.0 * g_yyyy_0_xxxyyyz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyyzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxyyzzz_0[i] = 3.0 * g_yyyy_0_xxxyyzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyyzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxyzzzz_0[i] = 4.0 * g_yyyy_0_xxxyzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxyzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxxzzzzz_0[i] = 5.0 * g_yyyy_0_xxxzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxxzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxyyyyyy_0[i] = g_yyyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyyyz_0_xxyyyyyz_0[i] = g_yyyy_0_xxyyyyy_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyyz_1[i] * wa_z[i];

        g_yyyyz_0_xxyyyyzz_0[i] = 2.0 * g_yyyy_0_xxyyyyz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyyzz_1[i] * wa_z[i];

        g_yyyyz_0_xxyyyzzz_0[i] = 3.0 * g_yyyy_0_xxyyyzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyyzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxyyzzzz_0[i] = 4.0 * g_yyyy_0_xxyyzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyyzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxyzzzzz_0[i] = 5.0 * g_yyyy_0_xxyzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxyzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xxzzzzzz_0[i] = 6.0 * g_yyyy_0_xxzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xxzzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xyyyyyyy_0[i] = g_yyyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyyyz_0_xyyyyyyz_0[i] = g_yyyy_0_xyyyyyy_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyyz_1[i] * wa_z[i];

        g_yyyyz_0_xyyyyyzz_0[i] = 2.0 * g_yyyy_0_xyyyyyz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyyzz_1[i] * wa_z[i];

        g_yyyyz_0_xyyyyzzz_0[i] = 3.0 * g_yyyy_0_xyyyyzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyyzzz_1[i] * wa_z[i];

        g_yyyyz_0_xyyyzzzz_0[i] = 4.0 * g_yyyy_0_xyyyzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyyzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xyyzzzzz_0[i] = 5.0 * g_yyyy_0_xyyzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyyzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xyzzzzzz_0[i] = 6.0 * g_yyyy_0_xyzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xyzzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_xzzzzzzz_0[i] = 7.0 * g_yyyy_0_xzzzzzz_1[i] * fi_acd_0 + g_yyyy_0_xzzzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_yyyyyyyy_0[i] = g_yyyy_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyyyz_0_yyyyyyyz_0[i] = g_yyyy_0_yyyyyyy_1[i] * fi_acd_0 + g_yyyy_0_yyyyyyyz_1[i] * wa_z[i];

        g_yyyyz_0_yyyyyyzz_0[i] = 2.0 * g_yyyy_0_yyyyyyz_1[i] * fi_acd_0 + g_yyyy_0_yyyyyyzz_1[i] * wa_z[i];

        g_yyyyz_0_yyyyyzzz_0[i] = 3.0 * g_yyyy_0_yyyyyzz_1[i] * fi_acd_0 + g_yyyy_0_yyyyyzzz_1[i] * wa_z[i];

        g_yyyyz_0_yyyyzzzz_0[i] = 4.0 * g_yyyy_0_yyyyzzz_1[i] * fi_acd_0 + g_yyyy_0_yyyyzzzz_1[i] * wa_z[i];

        g_yyyyz_0_yyyzzzzz_0[i] = 5.0 * g_yyyy_0_yyyzzzz_1[i] * fi_acd_0 + g_yyyy_0_yyyzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_yyzzzzzz_0[i] = 6.0 * g_yyyy_0_yyzzzzz_1[i] * fi_acd_0 + g_yyyy_0_yyzzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_yzzzzzzz_0[i] = 7.0 * g_yyyy_0_yzzzzzz_1[i] * fi_acd_0 + g_yyyy_0_yzzzzzzz_1[i] * wa_z[i];

        g_yyyyz_0_zzzzzzzz_0[i] = 8.0 * g_yyyy_0_zzzzzzz_1[i] * fi_acd_0 + g_yyyy_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 765-810 components of targeted buffer : HSL

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

    #pragma omp simd aligned(g_yyy_0_xxxxxxxy_0, g_yyy_0_xxxxxxxy_1, g_yyy_0_xxxxxxyy_0, g_yyy_0_xxxxxxyy_1, g_yyy_0_xxxxxyyy_0, g_yyy_0_xxxxxyyy_1, g_yyy_0_xxxxyyyy_0, g_yyy_0_xxxxyyyy_1, g_yyy_0_xxxyyyyy_0, g_yyy_0_xxxyyyyy_1, g_yyy_0_xxyyyyyy_0, g_yyy_0_xxyyyyyy_1, g_yyy_0_xyyyyyyy_0, g_yyy_0_xyyyyyyy_1, g_yyy_0_yyyyyyyy_0, g_yyy_0_yyyyyyyy_1, g_yyyz_0_xxxxxxxy_1, g_yyyz_0_xxxxxxyy_1, g_yyyz_0_xxxxxyyy_1, g_yyyz_0_xxxxyyyy_1, g_yyyz_0_xxxyyyyy_1, g_yyyz_0_xxyyyyyy_1, g_yyyz_0_xyyyyyyy_1, g_yyyz_0_yyyyyyyy_1, g_yyyzz_0_xxxxxxxx_0, g_yyyzz_0_xxxxxxxy_0, g_yyyzz_0_xxxxxxxz_0, g_yyyzz_0_xxxxxxyy_0, g_yyyzz_0_xxxxxxyz_0, g_yyyzz_0_xxxxxxzz_0, g_yyyzz_0_xxxxxyyy_0, g_yyyzz_0_xxxxxyyz_0, g_yyyzz_0_xxxxxyzz_0, g_yyyzz_0_xxxxxzzz_0, g_yyyzz_0_xxxxyyyy_0, g_yyyzz_0_xxxxyyyz_0, g_yyyzz_0_xxxxyyzz_0, g_yyyzz_0_xxxxyzzz_0, g_yyyzz_0_xxxxzzzz_0, g_yyyzz_0_xxxyyyyy_0, g_yyyzz_0_xxxyyyyz_0, g_yyyzz_0_xxxyyyzz_0, g_yyyzz_0_xxxyyzzz_0, g_yyyzz_0_xxxyzzzz_0, g_yyyzz_0_xxxzzzzz_0, g_yyyzz_0_xxyyyyyy_0, g_yyyzz_0_xxyyyyyz_0, g_yyyzz_0_xxyyyyzz_0, g_yyyzz_0_xxyyyzzz_0, g_yyyzz_0_xxyyzzzz_0, g_yyyzz_0_xxyzzzzz_0, g_yyyzz_0_xxzzzzzz_0, g_yyyzz_0_xyyyyyyy_0, g_yyyzz_0_xyyyyyyz_0, g_yyyzz_0_xyyyyyzz_0, g_yyyzz_0_xyyyyzzz_0, g_yyyzz_0_xyyyzzzz_0, g_yyyzz_0_xyyzzzzz_0, g_yyyzz_0_xyzzzzzz_0, g_yyyzz_0_xzzzzzzz_0, g_yyyzz_0_yyyyyyyy_0, g_yyyzz_0_yyyyyyyz_0, g_yyyzz_0_yyyyyyzz_0, g_yyyzz_0_yyyyyzzz_0, g_yyyzz_0_yyyyzzzz_0, g_yyyzz_0_yyyzzzzz_0, g_yyyzz_0_yyzzzzzz_0, g_yyyzz_0_yzzzzzzz_0, g_yyyzz_0_zzzzzzzz_0, g_yyzz_0_xxxxxxxx_1, g_yyzz_0_xxxxxxxz_1, g_yyzz_0_xxxxxxyz_1, g_yyzz_0_xxxxxxz_1, g_yyzz_0_xxxxxxzz_1, g_yyzz_0_xxxxxyyz_1, g_yyzz_0_xxxxxyz_1, g_yyzz_0_xxxxxyzz_1, g_yyzz_0_xxxxxzz_1, g_yyzz_0_xxxxxzzz_1, g_yyzz_0_xxxxyyyz_1, g_yyzz_0_xxxxyyz_1, g_yyzz_0_xxxxyyzz_1, g_yyzz_0_xxxxyzz_1, g_yyzz_0_xxxxyzzz_1, g_yyzz_0_xxxxzzz_1, g_yyzz_0_xxxxzzzz_1, g_yyzz_0_xxxyyyyz_1, g_yyzz_0_xxxyyyz_1, g_yyzz_0_xxxyyyzz_1, g_yyzz_0_xxxyyzz_1, g_yyzz_0_xxxyyzzz_1, g_yyzz_0_xxxyzzz_1, g_yyzz_0_xxxyzzzz_1, g_yyzz_0_xxxzzzz_1, g_yyzz_0_xxxzzzzz_1, g_yyzz_0_xxyyyyyz_1, g_yyzz_0_xxyyyyz_1, g_yyzz_0_xxyyyyzz_1, g_yyzz_0_xxyyyzz_1, g_yyzz_0_xxyyyzzz_1, g_yyzz_0_xxyyzzz_1, g_yyzz_0_xxyyzzzz_1, g_yyzz_0_xxyzzzz_1, g_yyzz_0_xxyzzzzz_1, g_yyzz_0_xxzzzzz_1, g_yyzz_0_xxzzzzzz_1, g_yyzz_0_xyyyyyyz_1, g_yyzz_0_xyyyyyz_1, g_yyzz_0_xyyyyyzz_1, g_yyzz_0_xyyyyzz_1, g_yyzz_0_xyyyyzzz_1, g_yyzz_0_xyyyzzz_1, g_yyzz_0_xyyyzzzz_1, g_yyzz_0_xyyzzzz_1, g_yyzz_0_xyyzzzzz_1, g_yyzz_0_xyzzzzz_1, g_yyzz_0_xyzzzzzz_1, g_yyzz_0_xzzzzzz_1, g_yyzz_0_xzzzzzzz_1, g_yyzz_0_yyyyyyyz_1, g_yyzz_0_yyyyyyz_1, g_yyzz_0_yyyyyyzz_1, g_yyzz_0_yyyyyzz_1, g_yyzz_0_yyyyyzzz_1, g_yyzz_0_yyyyzzz_1, g_yyzz_0_yyyyzzzz_1, g_yyzz_0_yyyzzzz_1, g_yyzz_0_yyyzzzzz_1, g_yyzz_0_yyzzzzz_1, g_yyzz_0_yyzzzzzz_1, g_yyzz_0_yzzzzzz_1, g_yyzz_0_yzzzzzzz_1, g_yyzz_0_zzzzzzz_1, g_yyzz_0_zzzzzzzz_1, g_yzz_0_xxxxxxxx_0, g_yzz_0_xxxxxxxx_1, g_yzz_0_xxxxxxxz_0, g_yzz_0_xxxxxxxz_1, g_yzz_0_xxxxxxyz_0, g_yzz_0_xxxxxxyz_1, g_yzz_0_xxxxxxzz_0, g_yzz_0_xxxxxxzz_1, g_yzz_0_xxxxxyyz_0, g_yzz_0_xxxxxyyz_1, g_yzz_0_xxxxxyzz_0, g_yzz_0_xxxxxyzz_1, g_yzz_0_xxxxxzzz_0, g_yzz_0_xxxxxzzz_1, g_yzz_0_xxxxyyyz_0, g_yzz_0_xxxxyyyz_1, g_yzz_0_xxxxyyzz_0, g_yzz_0_xxxxyyzz_1, g_yzz_0_xxxxyzzz_0, g_yzz_0_xxxxyzzz_1, g_yzz_0_xxxxzzzz_0, g_yzz_0_xxxxzzzz_1, g_yzz_0_xxxyyyyz_0, g_yzz_0_xxxyyyyz_1, g_yzz_0_xxxyyyzz_0, g_yzz_0_xxxyyyzz_1, g_yzz_0_xxxyyzzz_0, g_yzz_0_xxxyyzzz_1, g_yzz_0_xxxyzzzz_0, g_yzz_0_xxxyzzzz_1, g_yzz_0_xxxzzzzz_0, g_yzz_0_xxxzzzzz_1, g_yzz_0_xxyyyyyz_0, g_yzz_0_xxyyyyyz_1, g_yzz_0_xxyyyyzz_0, g_yzz_0_xxyyyyzz_1, g_yzz_0_xxyyyzzz_0, g_yzz_0_xxyyyzzz_1, g_yzz_0_xxyyzzzz_0, g_yzz_0_xxyyzzzz_1, g_yzz_0_xxyzzzzz_0, g_yzz_0_xxyzzzzz_1, g_yzz_0_xxzzzzzz_0, g_yzz_0_xxzzzzzz_1, g_yzz_0_xyyyyyyz_0, g_yzz_0_xyyyyyyz_1, g_yzz_0_xyyyyyzz_0, g_yzz_0_xyyyyyzz_1, g_yzz_0_xyyyyzzz_0, g_yzz_0_xyyyyzzz_1, g_yzz_0_xyyyzzzz_0, g_yzz_0_xyyyzzzz_1, g_yzz_0_xyyzzzzz_0, g_yzz_0_xyyzzzzz_1, g_yzz_0_xyzzzzzz_0, g_yzz_0_xyzzzzzz_1, g_yzz_0_xzzzzzzz_0, g_yzz_0_xzzzzzzz_1, g_yzz_0_yyyyyyyz_0, g_yzz_0_yyyyyyyz_1, g_yzz_0_yyyyyyzz_0, g_yzz_0_yyyyyyzz_1, g_yzz_0_yyyyyzzz_0, g_yzz_0_yyyyyzzz_1, g_yzz_0_yyyyzzzz_0, g_yzz_0_yyyyzzzz_1, g_yzz_0_yyyzzzzz_0, g_yzz_0_yyyzzzzz_1, g_yzz_0_yyzzzzzz_0, g_yzz_0_yyzzzzzz_1, g_yzz_0_yzzzzzzz_0, g_yzz_0_yzzzzzzz_1, g_yzz_0_zzzzzzzz_0, g_yzz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzz_0_xxxxxxxx_0[i] = 2.0 * g_yzz_0_xxxxxxxx_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxxxx_1[i] * fz_be_0 + g_yyzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyyzz_0_xxxxxxxy_0[i] = g_yyy_0_xxxxxxxy_0[i] * fbe_0 - g_yyy_0_xxxxxxxy_1[i] * fz_be_0 + g_yyyz_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyyzz_0_xxxxxxxz_0[i] = 2.0 * g_yzz_0_xxxxxxxz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxxxz_1[i] * fz_be_0 + g_yyzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxxxyy_0[i] = g_yyy_0_xxxxxxyy_0[i] * fbe_0 - g_yyy_0_xxxxxxyy_1[i] * fz_be_0 + g_yyyz_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyyzz_0_xxxxxxyz_0[i] = 2.0 * g_yzz_0_xxxxxxyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxxyz_1[i] * fz_be_0 + g_yyzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxxxzz_0[i] = 2.0 * g_yzz_0_xxxxxxzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxxzz_1[i] * fz_be_0 + g_yyzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxxyyy_0[i] = g_yyy_0_xxxxxyyy_0[i] * fbe_0 - g_yyy_0_xxxxxyyy_1[i] * fz_be_0 + g_yyyz_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyyzz_0_xxxxxyyz_0[i] = 2.0 * g_yzz_0_xxxxxyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxxyzz_0[i] = 2.0 * g_yzz_0_xxxxxyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxyzz_1[i] * fz_be_0 + g_yyzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxxzzz_0[i] = 2.0 * g_yzz_0_xxxxxzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxxzzz_1[i] * fz_be_0 + g_yyzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxyyyy_0[i] = g_yyy_0_xxxxyyyy_0[i] * fbe_0 - g_yyy_0_xxxxyyyy_1[i] * fz_be_0 + g_yyyz_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyyzz_0_xxxxyyyz_0[i] = 2.0 * g_yzz_0_xxxxyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxyyzz_0[i] = 2.0 * g_yzz_0_xxxxyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxyzzz_0[i] = 2.0 * g_yzz_0_xxxxyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxyzzz_1[i] * fz_be_0 + g_yyzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxxzzzz_0[i] = 2.0 * g_yzz_0_xxxxzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxxzzzz_1[i] * fz_be_0 + g_yyzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxyyyyy_0[i] = g_yyy_0_xxxyyyyy_0[i] * fbe_0 - g_yyy_0_xxxyyyyy_1[i] * fz_be_0 + g_yyyz_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyyzz_0_xxxyyyyz_0[i] = 2.0 * g_yzz_0_xxxyyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyyzz_0_xxxyyyzz_0[i] = 2.0 * g_yzz_0_xxxyyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxyyzzz_0[i] = 2.0 * g_yzz_0_xxxyyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxyzzzz_0[i] = 2.0 * g_yzz_0_xxxyzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxyzzzz_1[i] * fz_be_0 + g_yyzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxxzzzzz_0[i] = 2.0 * g_yzz_0_xxxzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxxzzzzz_1[i] * fz_be_0 + g_yyzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxyyyyyy_0[i] = g_yyy_0_xxyyyyyy_0[i] * fbe_0 - g_yyy_0_xxyyyyyy_1[i] * fz_be_0 + g_yyyz_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyyzz_0_xxyyyyyz_0[i] = 2.0 * g_yzz_0_xxyyyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyyzz_0_xxyyyyzz_0[i] = 2.0 * g_yzz_0_xxyyyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyyzz_0_xxyyyzzz_0[i] = 2.0 * g_yzz_0_xxyyyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxyyzzzz_0[i] = 2.0 * g_yzz_0_xxyyzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxyzzzzz_0[i] = 2.0 * g_yzz_0_xxyzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxyzzzzz_1[i] * fz_be_0 + g_yyzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xxzzzzzz_0[i] = 2.0 * g_yzz_0_xxzzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xxzzzzzz_1[i] * fz_be_0 + g_yyzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xyyyyyyy_0[i] = g_yyy_0_xyyyyyyy_0[i] * fbe_0 - g_yyy_0_xyyyyyyy_1[i] * fz_be_0 + g_yyyz_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyyzz_0_xyyyyyyz_0[i] = 2.0 * g_yzz_0_xyyyyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyyzz_0_xyyyyyzz_0[i] = 2.0 * g_yzz_0_xyyyyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyyzz_0_xyyyyzzz_0[i] = 2.0 * g_yzz_0_xyyyyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyyzz_0_xyyyzzzz_0[i] = 2.0 * g_yzz_0_xyyyzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xyyzzzzz_0[i] = 2.0 * g_yzz_0_xyyzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xyzzzzzz_0[i] = 2.0 * g_yzz_0_xyzzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xyzzzzzz_1[i] * fz_be_0 + g_yyzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_xzzzzzzz_0[i] = 2.0 * g_yzz_0_xzzzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_xzzzzzzz_1[i] * fz_be_0 + g_yyzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_yyyyyyyy_0[i] = g_yyy_0_yyyyyyyy_0[i] * fbe_0 - g_yyy_0_yyyyyyyy_1[i] * fz_be_0 + g_yyyz_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyyzz_0_yyyyyyyz_0[i] = 2.0 * g_yzz_0_yyyyyyyz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yyzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyyzz_0_yyyyyyzz_0[i] = 2.0 * g_yzz_0_yyyyyyzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yyzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyyzz_0_yyyyyzzz_0[i] = 2.0 * g_yzz_0_yyyyyzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yyzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyyzz_0_yyyyzzzz_0[i] = 2.0 * g_yzz_0_yyyyzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yyzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyyzz_0_yyyzzzzz_0[i] = 2.0 * g_yzz_0_yyyzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yyzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_yyzzzzzz_0[i] = 2.0 * g_yzz_0_yyzzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_yzzzzzzz_0[i] = 2.0 * g_yzz_0_yzzzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_yzzzzzzz_1[i] * fz_be_0 + g_yyzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyyzz_0_zzzzzzzz_0[i] = 2.0 * g_yzz_0_zzzzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_0_zzzzzzzz_1[i] * fz_be_0 + g_yyzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 810-855 components of targeted buffer : HSL

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

    #pragma omp simd aligned(g_yyz_0_xxxxxxxy_0, g_yyz_0_xxxxxxxy_1, g_yyz_0_xxxxxxyy_0, g_yyz_0_xxxxxxyy_1, g_yyz_0_xxxxxyyy_0, g_yyz_0_xxxxxyyy_1, g_yyz_0_xxxxyyyy_0, g_yyz_0_xxxxyyyy_1, g_yyz_0_xxxyyyyy_0, g_yyz_0_xxxyyyyy_1, g_yyz_0_xxyyyyyy_0, g_yyz_0_xxyyyyyy_1, g_yyz_0_xyyyyyyy_0, g_yyz_0_xyyyyyyy_1, g_yyz_0_yyyyyyyy_0, g_yyz_0_yyyyyyyy_1, g_yyzz_0_xxxxxxxy_1, g_yyzz_0_xxxxxxyy_1, g_yyzz_0_xxxxxyyy_1, g_yyzz_0_xxxxyyyy_1, g_yyzz_0_xxxyyyyy_1, g_yyzz_0_xxyyyyyy_1, g_yyzz_0_xyyyyyyy_1, g_yyzz_0_yyyyyyyy_1, g_yyzzz_0_xxxxxxxx_0, g_yyzzz_0_xxxxxxxy_0, g_yyzzz_0_xxxxxxxz_0, g_yyzzz_0_xxxxxxyy_0, g_yyzzz_0_xxxxxxyz_0, g_yyzzz_0_xxxxxxzz_0, g_yyzzz_0_xxxxxyyy_0, g_yyzzz_0_xxxxxyyz_0, g_yyzzz_0_xxxxxyzz_0, g_yyzzz_0_xxxxxzzz_0, g_yyzzz_0_xxxxyyyy_0, g_yyzzz_0_xxxxyyyz_0, g_yyzzz_0_xxxxyyzz_0, g_yyzzz_0_xxxxyzzz_0, g_yyzzz_0_xxxxzzzz_0, g_yyzzz_0_xxxyyyyy_0, g_yyzzz_0_xxxyyyyz_0, g_yyzzz_0_xxxyyyzz_0, g_yyzzz_0_xxxyyzzz_0, g_yyzzz_0_xxxyzzzz_0, g_yyzzz_0_xxxzzzzz_0, g_yyzzz_0_xxyyyyyy_0, g_yyzzz_0_xxyyyyyz_0, g_yyzzz_0_xxyyyyzz_0, g_yyzzz_0_xxyyyzzz_0, g_yyzzz_0_xxyyzzzz_0, g_yyzzz_0_xxyzzzzz_0, g_yyzzz_0_xxzzzzzz_0, g_yyzzz_0_xyyyyyyy_0, g_yyzzz_0_xyyyyyyz_0, g_yyzzz_0_xyyyyyzz_0, g_yyzzz_0_xyyyyzzz_0, g_yyzzz_0_xyyyzzzz_0, g_yyzzz_0_xyyzzzzz_0, g_yyzzz_0_xyzzzzzz_0, g_yyzzz_0_xzzzzzzz_0, g_yyzzz_0_yyyyyyyy_0, g_yyzzz_0_yyyyyyyz_0, g_yyzzz_0_yyyyyyzz_0, g_yyzzz_0_yyyyyzzz_0, g_yyzzz_0_yyyyzzzz_0, g_yyzzz_0_yyyzzzzz_0, g_yyzzz_0_yyzzzzzz_0, g_yyzzz_0_yzzzzzzz_0, g_yyzzz_0_zzzzzzzz_0, g_yzzz_0_xxxxxxxx_1, g_yzzz_0_xxxxxxxz_1, g_yzzz_0_xxxxxxyz_1, g_yzzz_0_xxxxxxz_1, g_yzzz_0_xxxxxxzz_1, g_yzzz_0_xxxxxyyz_1, g_yzzz_0_xxxxxyz_1, g_yzzz_0_xxxxxyzz_1, g_yzzz_0_xxxxxzz_1, g_yzzz_0_xxxxxzzz_1, g_yzzz_0_xxxxyyyz_1, g_yzzz_0_xxxxyyz_1, g_yzzz_0_xxxxyyzz_1, g_yzzz_0_xxxxyzz_1, g_yzzz_0_xxxxyzzz_1, g_yzzz_0_xxxxzzz_1, g_yzzz_0_xxxxzzzz_1, g_yzzz_0_xxxyyyyz_1, g_yzzz_0_xxxyyyz_1, g_yzzz_0_xxxyyyzz_1, g_yzzz_0_xxxyyzz_1, g_yzzz_0_xxxyyzzz_1, g_yzzz_0_xxxyzzz_1, g_yzzz_0_xxxyzzzz_1, g_yzzz_0_xxxzzzz_1, g_yzzz_0_xxxzzzzz_1, g_yzzz_0_xxyyyyyz_1, g_yzzz_0_xxyyyyz_1, g_yzzz_0_xxyyyyzz_1, g_yzzz_0_xxyyyzz_1, g_yzzz_0_xxyyyzzz_1, g_yzzz_0_xxyyzzz_1, g_yzzz_0_xxyyzzzz_1, g_yzzz_0_xxyzzzz_1, g_yzzz_0_xxyzzzzz_1, g_yzzz_0_xxzzzzz_1, g_yzzz_0_xxzzzzzz_1, g_yzzz_0_xyyyyyyz_1, g_yzzz_0_xyyyyyz_1, g_yzzz_0_xyyyyyzz_1, g_yzzz_0_xyyyyzz_1, g_yzzz_0_xyyyyzzz_1, g_yzzz_0_xyyyzzz_1, g_yzzz_0_xyyyzzzz_1, g_yzzz_0_xyyzzzz_1, g_yzzz_0_xyyzzzzz_1, g_yzzz_0_xyzzzzz_1, g_yzzz_0_xyzzzzzz_1, g_yzzz_0_xzzzzzz_1, g_yzzz_0_xzzzzzzz_1, g_yzzz_0_yyyyyyyz_1, g_yzzz_0_yyyyyyz_1, g_yzzz_0_yyyyyyzz_1, g_yzzz_0_yyyyyzz_1, g_yzzz_0_yyyyyzzz_1, g_yzzz_0_yyyyzzz_1, g_yzzz_0_yyyyzzzz_1, g_yzzz_0_yyyzzzz_1, g_yzzz_0_yyyzzzzz_1, g_yzzz_0_yyzzzzz_1, g_yzzz_0_yyzzzzzz_1, g_yzzz_0_yzzzzzz_1, g_yzzz_0_yzzzzzzz_1, g_yzzz_0_zzzzzzz_1, g_yzzz_0_zzzzzzzz_1, g_zzz_0_xxxxxxxx_0, g_zzz_0_xxxxxxxx_1, g_zzz_0_xxxxxxxz_0, g_zzz_0_xxxxxxxz_1, g_zzz_0_xxxxxxyz_0, g_zzz_0_xxxxxxyz_1, g_zzz_0_xxxxxxzz_0, g_zzz_0_xxxxxxzz_1, g_zzz_0_xxxxxyyz_0, g_zzz_0_xxxxxyyz_1, g_zzz_0_xxxxxyzz_0, g_zzz_0_xxxxxyzz_1, g_zzz_0_xxxxxzzz_0, g_zzz_0_xxxxxzzz_1, g_zzz_0_xxxxyyyz_0, g_zzz_0_xxxxyyyz_1, g_zzz_0_xxxxyyzz_0, g_zzz_0_xxxxyyzz_1, g_zzz_0_xxxxyzzz_0, g_zzz_0_xxxxyzzz_1, g_zzz_0_xxxxzzzz_0, g_zzz_0_xxxxzzzz_1, g_zzz_0_xxxyyyyz_0, g_zzz_0_xxxyyyyz_1, g_zzz_0_xxxyyyzz_0, g_zzz_0_xxxyyyzz_1, g_zzz_0_xxxyyzzz_0, g_zzz_0_xxxyyzzz_1, g_zzz_0_xxxyzzzz_0, g_zzz_0_xxxyzzzz_1, g_zzz_0_xxxzzzzz_0, g_zzz_0_xxxzzzzz_1, g_zzz_0_xxyyyyyz_0, g_zzz_0_xxyyyyyz_1, g_zzz_0_xxyyyyzz_0, g_zzz_0_xxyyyyzz_1, g_zzz_0_xxyyyzzz_0, g_zzz_0_xxyyyzzz_1, g_zzz_0_xxyyzzzz_0, g_zzz_0_xxyyzzzz_1, g_zzz_0_xxyzzzzz_0, g_zzz_0_xxyzzzzz_1, g_zzz_0_xxzzzzzz_0, g_zzz_0_xxzzzzzz_1, g_zzz_0_xyyyyyyz_0, g_zzz_0_xyyyyyyz_1, g_zzz_0_xyyyyyzz_0, g_zzz_0_xyyyyyzz_1, g_zzz_0_xyyyyzzz_0, g_zzz_0_xyyyyzzz_1, g_zzz_0_xyyyzzzz_0, g_zzz_0_xyyyzzzz_1, g_zzz_0_xyyzzzzz_0, g_zzz_0_xyyzzzzz_1, g_zzz_0_xyzzzzzz_0, g_zzz_0_xyzzzzzz_1, g_zzz_0_xzzzzzzz_0, g_zzz_0_xzzzzzzz_1, g_zzz_0_yyyyyyyz_0, g_zzz_0_yyyyyyyz_1, g_zzz_0_yyyyyyzz_0, g_zzz_0_yyyyyyzz_1, g_zzz_0_yyyyyzzz_0, g_zzz_0_yyyyyzzz_1, g_zzz_0_yyyyzzzz_0, g_zzz_0_yyyyzzzz_1, g_zzz_0_yyyzzzzz_0, g_zzz_0_yyyzzzzz_1, g_zzz_0_yyzzzzzz_0, g_zzz_0_yyzzzzzz_1, g_zzz_0_yzzzzzzz_0, g_zzz_0_yzzzzzzz_1, g_zzz_0_zzzzzzzz_0, g_zzz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzz_0_xxxxxxxx_0[i] = g_zzz_0_xxxxxxxx_0[i] * fbe_0 - g_zzz_0_xxxxxxxx_1[i] * fz_be_0 + g_yzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyzzz_0_xxxxxxxy_0[i] = 2.0 * g_yyz_0_xxxxxxxy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxxxxxy_1[i] * fz_be_0 + g_yyzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyzzz_0_xxxxxxxz_0[i] = g_zzz_0_xxxxxxxz_0[i] * fbe_0 - g_zzz_0_xxxxxxxz_1[i] * fz_be_0 + g_yzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxxxyy_0[i] = 2.0 * g_yyz_0_xxxxxxyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxxxxyy_1[i] * fz_be_0 + g_yyzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyzzz_0_xxxxxxyz_0[i] = g_zzz_0_xxxxxxyz_0[i] * fbe_0 - g_zzz_0_xxxxxxyz_1[i] * fz_be_0 + g_yzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yzzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxxxzz_0[i] = g_zzz_0_xxxxxxzz_0[i] * fbe_0 - g_zzz_0_xxxxxxzz_1[i] * fz_be_0 + g_yzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxxyyy_0[i] = 2.0 * g_yyz_0_xxxxxyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxxxyyy_1[i] * fz_be_0 + g_yyzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyzzz_0_xxxxxyyz_0[i] = g_zzz_0_xxxxxyyz_0[i] * fbe_0 - g_zzz_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yzzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxxyzz_0[i] = g_zzz_0_xxxxxyzz_0[i] * fbe_0 - g_zzz_0_xxxxxyzz_1[i] * fz_be_0 + g_yzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yzzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxxzzz_0[i] = g_zzz_0_xxxxxzzz_0[i] * fbe_0 - g_zzz_0_xxxxxzzz_1[i] * fz_be_0 + g_yzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxyyyy_0[i] = 2.0 * g_yyz_0_xxxxyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxxyyyy_1[i] * fz_be_0 + g_yyzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyzzz_0_xxxxyyyz_0[i] = g_zzz_0_xxxxyyyz_0[i] * fbe_0 - g_zzz_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yzzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxyyzz_0[i] = g_zzz_0_xxxxyyzz_0[i] * fbe_0 - g_zzz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yzzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxyzzz_0[i] = g_zzz_0_xxxxyzzz_0[i] * fbe_0 - g_zzz_0_xxxxyzzz_1[i] * fz_be_0 + g_yzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yzzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxxzzzz_0[i] = g_zzz_0_xxxxzzzz_0[i] * fbe_0 - g_zzz_0_xxxxzzzz_1[i] * fz_be_0 + g_yzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxyyyyy_0[i] = 2.0 * g_yyz_0_xxxyyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxxyyyyy_1[i] * fz_be_0 + g_yyzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyzzz_0_xxxyyyyz_0[i] = g_zzz_0_xxxyyyyz_0[i] * fbe_0 - g_zzz_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yzzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyzzz_0_xxxyyyzz_0[i] = g_zzz_0_xxxyyyzz_0[i] * fbe_0 - g_zzz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yzzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxyyzzz_0[i] = g_zzz_0_xxxyyzzz_0[i] * fbe_0 - g_zzz_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yzzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxyzzzz_0[i] = g_zzz_0_xxxyzzzz_0[i] * fbe_0 - g_zzz_0_xxxyzzzz_1[i] * fz_be_0 + g_yzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yzzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxxzzzzz_0[i] = g_zzz_0_xxxzzzzz_0[i] * fbe_0 - g_zzz_0_xxxzzzzz_1[i] * fz_be_0 + g_yzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxyyyyyy_0[i] = 2.0 * g_yyz_0_xxyyyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xxyyyyyy_1[i] * fz_be_0 + g_yyzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyzzz_0_xxyyyyyz_0[i] = g_zzz_0_xxyyyyyz_0[i] * fbe_0 - g_zzz_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yzzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyzzz_0_xxyyyyzz_0[i] = g_zzz_0_xxyyyyzz_0[i] * fbe_0 - g_zzz_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yzzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyzzz_0_xxyyyzzz_0[i] = g_zzz_0_xxyyyzzz_0[i] * fbe_0 - g_zzz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yzzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxyyzzzz_0[i] = g_zzz_0_xxyyzzzz_0[i] * fbe_0 - g_zzz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yzzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxyzzzzz_0[i] = g_zzz_0_xxyzzzzz_0[i] * fbe_0 - g_zzz_0_xxyzzzzz_1[i] * fz_be_0 + g_yzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yzzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xxzzzzzz_0[i] = g_zzz_0_xxzzzzzz_0[i] * fbe_0 - g_zzz_0_xxzzzzzz_1[i] * fz_be_0 + g_yzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xyyyyyyy_0[i] = 2.0 * g_yyz_0_xyyyyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_xyyyyyyy_1[i] * fz_be_0 + g_yyzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyzzz_0_xyyyyyyz_0[i] = g_zzz_0_xyyyyyyz_0[i] * fbe_0 - g_zzz_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yzzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyzzz_0_xyyyyyzz_0[i] = g_zzz_0_xyyyyyzz_0[i] * fbe_0 - g_zzz_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yzzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyzzz_0_xyyyyzzz_0[i] = g_zzz_0_xyyyyzzz_0[i] * fbe_0 - g_zzz_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yzzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyzzz_0_xyyyzzzz_0[i] = g_zzz_0_xyyyzzzz_0[i] * fbe_0 - g_zzz_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yzzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xyyzzzzz_0[i] = g_zzz_0_xyyzzzzz_0[i] * fbe_0 - g_zzz_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yzzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xyzzzzzz_0[i] = g_zzz_0_xyzzzzzz_0[i] * fbe_0 - g_zzz_0_xyzzzzzz_1[i] * fz_be_0 + g_yzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yzzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_xzzzzzzz_0[i] = g_zzz_0_xzzzzzzz_0[i] * fbe_0 - g_zzz_0_xzzzzzzz_1[i] * fz_be_0 + g_yzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_yyyyyyyy_0[i] = 2.0 * g_yyz_0_yyyyyyyy_0[i] * fbe_0 - 2.0 * g_yyz_0_yyyyyyyy_1[i] * fz_be_0 + g_yyzz_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyzzz_0_yyyyyyyz_0[i] = g_zzz_0_yyyyyyyz_0[i] * fbe_0 - g_zzz_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yzzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyzzz_0_yyyyyyzz_0[i] = g_zzz_0_yyyyyyzz_0[i] * fbe_0 - g_zzz_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yzzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyzzz_0_yyyyyzzz_0[i] = g_zzz_0_yyyyyzzz_0[i] * fbe_0 - g_zzz_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yzzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyzzz_0_yyyyzzzz_0[i] = g_zzz_0_yyyyzzzz_0[i] * fbe_0 - g_zzz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yzzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyzzz_0_yyyzzzzz_0[i] = g_zzz_0_yyyzzzzz_0[i] * fbe_0 - g_zzz_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yzzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_yyzzzzzz_0[i] = g_zzz_0_yyzzzzzz_0[i] * fbe_0 - g_zzz_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yzzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_yzzzzzzz_0[i] = g_zzz_0_yzzzzzzz_0[i] * fbe_0 - g_zzz_0_yzzzzzzz_1[i] * fz_be_0 + g_yzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yzzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyzzz_0_zzzzzzzz_0[i] = g_zzz_0_zzzzzzzz_0[i] * fbe_0 - g_zzz_0_zzzzzzzz_1[i] * fz_be_0 + g_yzzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 855-900 components of targeted buffer : HSL

    auto g_yzzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_hsl + 855);

    auto g_yzzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_hsl + 856);

    auto g_yzzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_hsl + 857);

    auto g_yzzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_hsl + 858);

    auto g_yzzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_hsl + 859);

    auto g_yzzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_hsl + 860);

    auto g_yzzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_hsl + 861);

    auto g_yzzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_hsl + 862);

    auto g_yzzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_hsl + 863);

    auto g_yzzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_hsl + 864);

    auto g_yzzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_hsl + 865);

    auto g_yzzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_hsl + 866);

    auto g_yzzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_hsl + 867);

    auto g_yzzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_hsl + 868);

    auto g_yzzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_hsl + 869);

    auto g_yzzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 870);

    auto g_yzzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 871);

    auto g_yzzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 872);

    auto g_yzzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 873);

    auto g_yzzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 874);

    auto g_yzzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 875);

    auto g_yzzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 876);

    auto g_yzzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 877);

    auto g_yzzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 878);

    auto g_yzzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 879);

    auto g_yzzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 880);

    auto g_yzzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 881);

    auto g_yzzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 882);

    auto g_yzzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 883);

    auto g_yzzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 884);

    auto g_yzzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 885);

    auto g_yzzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 886);

    auto g_yzzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 887);

    auto g_yzzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 888);

    auto g_yzzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 889);

    auto g_yzzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 890);

    auto g_yzzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_hsl + 891);

    auto g_yzzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_hsl + 892);

    auto g_yzzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_hsl + 893);

    auto g_yzzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_hsl + 894);

    auto g_yzzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_hsl + 895);

    auto g_yzzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 896);

    auto g_yzzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 897);

    auto g_yzzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 898);

    auto g_yzzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_hsl + 899);

    #pragma omp simd aligned(g_yzzzz_0_xxxxxxxx_0, g_yzzzz_0_xxxxxxxy_0, g_yzzzz_0_xxxxxxxz_0, g_yzzzz_0_xxxxxxyy_0, g_yzzzz_0_xxxxxxyz_0, g_yzzzz_0_xxxxxxzz_0, g_yzzzz_0_xxxxxyyy_0, g_yzzzz_0_xxxxxyyz_0, g_yzzzz_0_xxxxxyzz_0, g_yzzzz_0_xxxxxzzz_0, g_yzzzz_0_xxxxyyyy_0, g_yzzzz_0_xxxxyyyz_0, g_yzzzz_0_xxxxyyzz_0, g_yzzzz_0_xxxxyzzz_0, g_yzzzz_0_xxxxzzzz_0, g_yzzzz_0_xxxyyyyy_0, g_yzzzz_0_xxxyyyyz_0, g_yzzzz_0_xxxyyyzz_0, g_yzzzz_0_xxxyyzzz_0, g_yzzzz_0_xxxyzzzz_0, g_yzzzz_0_xxxzzzzz_0, g_yzzzz_0_xxyyyyyy_0, g_yzzzz_0_xxyyyyyz_0, g_yzzzz_0_xxyyyyzz_0, g_yzzzz_0_xxyyyzzz_0, g_yzzzz_0_xxyyzzzz_0, g_yzzzz_0_xxyzzzzz_0, g_yzzzz_0_xxzzzzzz_0, g_yzzzz_0_xyyyyyyy_0, g_yzzzz_0_xyyyyyyz_0, g_yzzzz_0_xyyyyyzz_0, g_yzzzz_0_xyyyyzzz_0, g_yzzzz_0_xyyyzzzz_0, g_yzzzz_0_xyyzzzzz_0, g_yzzzz_0_xyzzzzzz_0, g_yzzzz_0_xzzzzzzz_0, g_yzzzz_0_yyyyyyyy_0, g_yzzzz_0_yyyyyyyz_0, g_yzzzz_0_yyyyyyzz_0, g_yzzzz_0_yyyyyzzz_0, g_yzzzz_0_yyyyzzzz_0, g_yzzzz_0_yyyzzzzz_0, g_yzzzz_0_yyzzzzzz_0, g_yzzzz_0_yzzzzzzz_0, g_yzzzz_0_zzzzzzzz_0, g_zzzz_0_xxxxxxx_1, g_zzzz_0_xxxxxxxx_1, g_zzzz_0_xxxxxxxy_1, g_zzzz_0_xxxxxxxz_1, g_zzzz_0_xxxxxxy_1, g_zzzz_0_xxxxxxyy_1, g_zzzz_0_xxxxxxyz_1, g_zzzz_0_xxxxxxz_1, g_zzzz_0_xxxxxxzz_1, g_zzzz_0_xxxxxyy_1, g_zzzz_0_xxxxxyyy_1, g_zzzz_0_xxxxxyyz_1, g_zzzz_0_xxxxxyz_1, g_zzzz_0_xxxxxyzz_1, g_zzzz_0_xxxxxzz_1, g_zzzz_0_xxxxxzzz_1, g_zzzz_0_xxxxyyy_1, g_zzzz_0_xxxxyyyy_1, g_zzzz_0_xxxxyyyz_1, g_zzzz_0_xxxxyyz_1, g_zzzz_0_xxxxyyzz_1, g_zzzz_0_xxxxyzz_1, g_zzzz_0_xxxxyzzz_1, g_zzzz_0_xxxxzzz_1, g_zzzz_0_xxxxzzzz_1, g_zzzz_0_xxxyyyy_1, g_zzzz_0_xxxyyyyy_1, g_zzzz_0_xxxyyyyz_1, g_zzzz_0_xxxyyyz_1, g_zzzz_0_xxxyyyzz_1, g_zzzz_0_xxxyyzz_1, g_zzzz_0_xxxyyzzz_1, g_zzzz_0_xxxyzzz_1, g_zzzz_0_xxxyzzzz_1, g_zzzz_0_xxxzzzz_1, g_zzzz_0_xxxzzzzz_1, g_zzzz_0_xxyyyyy_1, g_zzzz_0_xxyyyyyy_1, g_zzzz_0_xxyyyyyz_1, g_zzzz_0_xxyyyyz_1, g_zzzz_0_xxyyyyzz_1, g_zzzz_0_xxyyyzz_1, g_zzzz_0_xxyyyzzz_1, g_zzzz_0_xxyyzzz_1, g_zzzz_0_xxyyzzzz_1, g_zzzz_0_xxyzzzz_1, g_zzzz_0_xxyzzzzz_1, g_zzzz_0_xxzzzzz_1, g_zzzz_0_xxzzzzzz_1, g_zzzz_0_xyyyyyy_1, g_zzzz_0_xyyyyyyy_1, g_zzzz_0_xyyyyyyz_1, g_zzzz_0_xyyyyyz_1, g_zzzz_0_xyyyyyzz_1, g_zzzz_0_xyyyyzz_1, g_zzzz_0_xyyyyzzz_1, g_zzzz_0_xyyyzzz_1, g_zzzz_0_xyyyzzzz_1, g_zzzz_0_xyyzzzz_1, g_zzzz_0_xyyzzzzz_1, g_zzzz_0_xyzzzzz_1, g_zzzz_0_xyzzzzzz_1, g_zzzz_0_xzzzzzz_1, g_zzzz_0_xzzzzzzz_1, g_zzzz_0_yyyyyyy_1, g_zzzz_0_yyyyyyyy_1, g_zzzz_0_yyyyyyyz_1, g_zzzz_0_yyyyyyz_1, g_zzzz_0_yyyyyyzz_1, g_zzzz_0_yyyyyzz_1, g_zzzz_0_yyyyyzzz_1, g_zzzz_0_yyyyzzz_1, g_zzzz_0_yyyyzzzz_1, g_zzzz_0_yyyzzzz_1, g_zzzz_0_yyyzzzzz_1, g_zzzz_0_yyzzzzz_1, g_zzzz_0_yyzzzzzz_1, g_zzzz_0_yzzzzzz_1, g_zzzz_0_yzzzzzzz_1, g_zzzz_0_zzzzzzz_1, g_zzzz_0_zzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzz_0_xxxxxxxx_0[i] = g_zzzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxxxy_0[i] = g_zzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxxy_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxxxz_0[i] = g_zzzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxxyy_0[i] = 2.0 * g_zzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxyy_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxxyz_0[i] = g_zzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxxzz_0[i] = g_zzzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxyyy_0[i] = 3.0 * g_zzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyyy_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxyyz_0[i] = 2.0 * g_zzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxyzz_0[i] = g_zzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxxzzz_0[i] = g_zzzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxyyyy_0[i] = 4.0 * g_zzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyyy_1[i] * wa_y[i];

        g_yzzzz_0_xxxxyyyz_0[i] = 3.0 * g_zzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxyyzz_0[i] = 2.0 * g_zzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxyzzz_0[i] = g_zzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxxzzzz_0[i] = g_zzzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxyyyyy_0[i] = 5.0 * g_zzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyyy_1[i] * wa_y[i];

        g_yzzzz_0_xxxyyyyz_0[i] = 4.0 * g_zzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yzzzz_0_xxxyyyzz_0[i] = 3.0 * g_zzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxyyzzz_0[i] = 2.0 * g_zzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxyzzzz_0[i] = g_zzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxxzzzzz_0[i] = g_zzzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxyyyyyy_0[i] = 6.0 * g_zzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyyy_1[i] * wa_y[i];

        g_yzzzz_0_xxyyyyyz_0[i] = 5.0 * g_zzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yzzzz_0_xxyyyyzz_0[i] = 4.0 * g_zzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yzzzz_0_xxyyyzzz_0[i] = 3.0 * g_zzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxyyzzzz_0[i] = 2.0 * g_zzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxyzzzzz_0[i] = g_zzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xxzzzzzz_0[i] = g_zzzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xyyyyyyy_0[i] = 7.0 * g_zzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyyy_1[i] * wa_y[i];

        g_yzzzz_0_xyyyyyyz_0[i] = 6.0 * g_zzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yzzzz_0_xyyyyyzz_0[i] = 5.0 * g_zzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yzzzz_0_xyyyyzzz_0[i] = 4.0 * g_zzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yzzzz_0_xyyyzzzz_0[i] = 3.0 * g_zzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xyyzzzzz_0[i] = 2.0 * g_zzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xyzzzzzz_0[i] = g_zzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_xzzzzzzz_0[i] = g_zzzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_yyyyyyyy_0[i] = 8.0 * g_zzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_zzzz_0_yyyyyyyy_1[i] * wa_y[i];

        g_yzzzz_0_yyyyyyyz_0[i] = 7.0 * g_zzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_zzzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yzzzz_0_yyyyyyzz_0[i] = 6.0 * g_zzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_zzzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yzzzz_0_yyyyyzzz_0[i] = 5.0 * g_zzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_zzzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yzzzz_0_yyyyzzzz_0[i] = 4.0 * g_zzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_zzzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yzzzz_0_yyyzzzzz_0[i] = 3.0 * g_zzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_zzzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_yyzzzzzz_0[i] = 2.0 * g_zzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_zzzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_yzzzzzzz_0[i] = g_zzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_zzzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yzzzz_0_zzzzzzzz_0[i] = g_zzzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 900-945 components of targeted buffer : HSL

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

    #pragma omp simd aligned(g_zzz_0_xxxxxxxx_0, g_zzz_0_xxxxxxxx_1, g_zzz_0_xxxxxxxy_0, g_zzz_0_xxxxxxxy_1, g_zzz_0_xxxxxxxz_0, g_zzz_0_xxxxxxxz_1, g_zzz_0_xxxxxxyy_0, g_zzz_0_xxxxxxyy_1, g_zzz_0_xxxxxxyz_0, g_zzz_0_xxxxxxyz_1, g_zzz_0_xxxxxxzz_0, g_zzz_0_xxxxxxzz_1, g_zzz_0_xxxxxyyy_0, g_zzz_0_xxxxxyyy_1, g_zzz_0_xxxxxyyz_0, g_zzz_0_xxxxxyyz_1, g_zzz_0_xxxxxyzz_0, g_zzz_0_xxxxxyzz_1, g_zzz_0_xxxxxzzz_0, g_zzz_0_xxxxxzzz_1, g_zzz_0_xxxxyyyy_0, g_zzz_0_xxxxyyyy_1, g_zzz_0_xxxxyyyz_0, g_zzz_0_xxxxyyyz_1, g_zzz_0_xxxxyyzz_0, g_zzz_0_xxxxyyzz_1, g_zzz_0_xxxxyzzz_0, g_zzz_0_xxxxyzzz_1, g_zzz_0_xxxxzzzz_0, g_zzz_0_xxxxzzzz_1, g_zzz_0_xxxyyyyy_0, g_zzz_0_xxxyyyyy_1, g_zzz_0_xxxyyyyz_0, g_zzz_0_xxxyyyyz_1, g_zzz_0_xxxyyyzz_0, g_zzz_0_xxxyyyzz_1, g_zzz_0_xxxyyzzz_0, g_zzz_0_xxxyyzzz_1, g_zzz_0_xxxyzzzz_0, g_zzz_0_xxxyzzzz_1, g_zzz_0_xxxzzzzz_0, g_zzz_0_xxxzzzzz_1, g_zzz_0_xxyyyyyy_0, g_zzz_0_xxyyyyyy_1, g_zzz_0_xxyyyyyz_0, g_zzz_0_xxyyyyyz_1, g_zzz_0_xxyyyyzz_0, g_zzz_0_xxyyyyzz_1, g_zzz_0_xxyyyzzz_0, g_zzz_0_xxyyyzzz_1, g_zzz_0_xxyyzzzz_0, g_zzz_0_xxyyzzzz_1, g_zzz_0_xxyzzzzz_0, g_zzz_0_xxyzzzzz_1, g_zzz_0_xxzzzzzz_0, g_zzz_0_xxzzzzzz_1, g_zzz_0_xyyyyyyy_0, g_zzz_0_xyyyyyyy_1, g_zzz_0_xyyyyyyz_0, g_zzz_0_xyyyyyyz_1, g_zzz_0_xyyyyyzz_0, g_zzz_0_xyyyyyzz_1, g_zzz_0_xyyyyzzz_0, g_zzz_0_xyyyyzzz_1, g_zzz_0_xyyyzzzz_0, g_zzz_0_xyyyzzzz_1, g_zzz_0_xyyzzzzz_0, g_zzz_0_xyyzzzzz_1, g_zzz_0_xyzzzzzz_0, g_zzz_0_xyzzzzzz_1, g_zzz_0_xzzzzzzz_0, g_zzz_0_xzzzzzzz_1, g_zzz_0_yyyyyyyy_0, g_zzz_0_yyyyyyyy_1, g_zzz_0_yyyyyyyz_0, g_zzz_0_yyyyyyyz_1, g_zzz_0_yyyyyyzz_0, g_zzz_0_yyyyyyzz_1, g_zzz_0_yyyyyzzz_0, g_zzz_0_yyyyyzzz_1, g_zzz_0_yyyyzzzz_0, g_zzz_0_yyyyzzzz_1, g_zzz_0_yyyzzzzz_0, g_zzz_0_yyyzzzzz_1, g_zzz_0_yyzzzzzz_0, g_zzz_0_yyzzzzzz_1, g_zzz_0_yzzzzzzz_0, g_zzz_0_yzzzzzzz_1, g_zzz_0_zzzzzzzz_0, g_zzz_0_zzzzzzzz_1, g_zzzz_0_xxxxxxx_1, g_zzzz_0_xxxxxxxx_1, g_zzzz_0_xxxxxxxy_1, g_zzzz_0_xxxxxxxz_1, g_zzzz_0_xxxxxxy_1, g_zzzz_0_xxxxxxyy_1, g_zzzz_0_xxxxxxyz_1, g_zzzz_0_xxxxxxz_1, g_zzzz_0_xxxxxxzz_1, g_zzzz_0_xxxxxyy_1, g_zzzz_0_xxxxxyyy_1, g_zzzz_0_xxxxxyyz_1, g_zzzz_0_xxxxxyz_1, g_zzzz_0_xxxxxyzz_1, g_zzzz_0_xxxxxzz_1, g_zzzz_0_xxxxxzzz_1, g_zzzz_0_xxxxyyy_1, g_zzzz_0_xxxxyyyy_1, g_zzzz_0_xxxxyyyz_1, g_zzzz_0_xxxxyyz_1, g_zzzz_0_xxxxyyzz_1, g_zzzz_0_xxxxyzz_1, g_zzzz_0_xxxxyzzz_1, g_zzzz_0_xxxxzzz_1, g_zzzz_0_xxxxzzzz_1, g_zzzz_0_xxxyyyy_1, g_zzzz_0_xxxyyyyy_1, g_zzzz_0_xxxyyyyz_1, g_zzzz_0_xxxyyyz_1, g_zzzz_0_xxxyyyzz_1, g_zzzz_0_xxxyyzz_1, g_zzzz_0_xxxyyzzz_1, g_zzzz_0_xxxyzzz_1, g_zzzz_0_xxxyzzzz_1, g_zzzz_0_xxxzzzz_1, g_zzzz_0_xxxzzzzz_1, g_zzzz_0_xxyyyyy_1, g_zzzz_0_xxyyyyyy_1, g_zzzz_0_xxyyyyyz_1, g_zzzz_0_xxyyyyz_1, g_zzzz_0_xxyyyyzz_1, g_zzzz_0_xxyyyzz_1, g_zzzz_0_xxyyyzzz_1, g_zzzz_0_xxyyzzz_1, g_zzzz_0_xxyyzzzz_1, g_zzzz_0_xxyzzzz_1, g_zzzz_0_xxyzzzzz_1, g_zzzz_0_xxzzzzz_1, g_zzzz_0_xxzzzzzz_1, g_zzzz_0_xyyyyyy_1, g_zzzz_0_xyyyyyyy_1, g_zzzz_0_xyyyyyyz_1, g_zzzz_0_xyyyyyz_1, g_zzzz_0_xyyyyyzz_1, g_zzzz_0_xyyyyzz_1, g_zzzz_0_xyyyyzzz_1, g_zzzz_0_xyyyzzz_1, g_zzzz_0_xyyyzzzz_1, g_zzzz_0_xyyzzzz_1, g_zzzz_0_xyyzzzzz_1, g_zzzz_0_xyzzzzz_1, g_zzzz_0_xyzzzzzz_1, g_zzzz_0_xzzzzzz_1, g_zzzz_0_xzzzzzzz_1, g_zzzz_0_yyyyyyy_1, g_zzzz_0_yyyyyyyy_1, g_zzzz_0_yyyyyyyz_1, g_zzzz_0_yyyyyyz_1, g_zzzz_0_yyyyyyzz_1, g_zzzz_0_yyyyyzz_1, g_zzzz_0_yyyyyzzz_1, g_zzzz_0_yyyyzzz_1, g_zzzz_0_yyyyzzzz_1, g_zzzz_0_yyyzzzz_1, g_zzzz_0_yyyzzzzz_1, g_zzzz_0_yyzzzzz_1, g_zzzz_0_yyzzzzzz_1, g_zzzz_0_yzzzzzz_1, g_zzzz_0_yzzzzzzz_1, g_zzzz_0_zzzzzzz_1, g_zzzz_0_zzzzzzzz_1, g_zzzzz_0_xxxxxxxx_0, g_zzzzz_0_xxxxxxxy_0, g_zzzzz_0_xxxxxxxz_0, g_zzzzz_0_xxxxxxyy_0, g_zzzzz_0_xxxxxxyz_0, g_zzzzz_0_xxxxxxzz_0, g_zzzzz_0_xxxxxyyy_0, g_zzzzz_0_xxxxxyyz_0, g_zzzzz_0_xxxxxyzz_0, g_zzzzz_0_xxxxxzzz_0, g_zzzzz_0_xxxxyyyy_0, g_zzzzz_0_xxxxyyyz_0, g_zzzzz_0_xxxxyyzz_0, g_zzzzz_0_xxxxyzzz_0, g_zzzzz_0_xxxxzzzz_0, g_zzzzz_0_xxxyyyyy_0, g_zzzzz_0_xxxyyyyz_0, g_zzzzz_0_xxxyyyzz_0, g_zzzzz_0_xxxyyzzz_0, g_zzzzz_0_xxxyzzzz_0, g_zzzzz_0_xxxzzzzz_0, g_zzzzz_0_xxyyyyyy_0, g_zzzzz_0_xxyyyyyz_0, g_zzzzz_0_xxyyyyzz_0, g_zzzzz_0_xxyyyzzz_0, g_zzzzz_0_xxyyzzzz_0, g_zzzzz_0_xxyzzzzz_0, g_zzzzz_0_xxzzzzzz_0, g_zzzzz_0_xyyyyyyy_0, g_zzzzz_0_xyyyyyyz_0, g_zzzzz_0_xyyyyyzz_0, g_zzzzz_0_xyyyyzzz_0, g_zzzzz_0_xyyyzzzz_0, g_zzzzz_0_xyyzzzzz_0, g_zzzzz_0_xyzzzzzz_0, g_zzzzz_0_xzzzzzzz_0, g_zzzzz_0_yyyyyyyy_0, g_zzzzz_0_yyyyyyyz_0, g_zzzzz_0_yyyyyyzz_0, g_zzzzz_0_yyyyyzzz_0, g_zzzzz_0_yyyyzzzz_0, g_zzzzz_0_yyyzzzzz_0, g_zzzzz_0_yyzzzzzz_0, g_zzzzz_0_yzzzzzzz_0, g_zzzzz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzz_0_xxxxxxxx_0[i] = 4.0 * g_zzz_0_xxxxxxxx_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxxxx_1[i] * fz_be_0 + g_zzzz_0_xxxxxxxx_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxxxy_0[i] = 4.0 * g_zzz_0_xxxxxxxy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxxxy_1[i] * fz_be_0 + g_zzzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxxxz_0[i] = 4.0 * g_zzz_0_xxxxxxxz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxxxz_1[i] * fz_be_0 + g_zzzz_0_xxxxxxx_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxxz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxxyy_0[i] = 4.0 * g_zzz_0_xxxxxxyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxxyy_1[i] * fz_be_0 + g_zzzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxxyz_0[i] = 4.0 * g_zzz_0_xxxxxxyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxxyz_1[i] * fz_be_0 + g_zzzz_0_xxxxxxy_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxyz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxxzz_0[i] = 4.0 * g_zzz_0_xxxxxxzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxxzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxxxxxz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxxzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxyyy_0[i] = 4.0 * g_zzz_0_xxxxxyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxyyy_1[i] * fz_be_0 + g_zzzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxyyz_0[i] = 4.0 * g_zzz_0_xxxxxyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxyyz_1[i] * fz_be_0 + g_zzzz_0_xxxxxyy_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyyz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxyzz_0[i] = 4.0 * g_zzz_0_xxxxxyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxxxxyz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxyzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxxzzz_0[i] = 4.0 * g_zzz_0_xxxxxzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxxzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xxxxxzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxxzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxyyyy_0[i] = 4.0 * g_zzz_0_xxxxyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxyyyy_1[i] * fz_be_0 + g_zzzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_zzzzz_0_xxxxyyyz_0[i] = 4.0 * g_zzz_0_xxxxyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxyyyz_1[i] * fz_be_0 + g_zzzz_0_xxxxyyy_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyyz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxyyzz_0[i] = 4.0 * g_zzz_0_xxxxyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxxxyyz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyyzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxyzzz_0[i] = 4.0 * g_zzz_0_xxxxyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xxxxyzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxyzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxxzzzz_0[i] = 4.0 * g_zzz_0_xxxxzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_xxxxzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxxzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxyyyyy_0[i] = 4.0 * g_zzz_0_xxxyyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyyyyy_1[i] * fz_be_0 + g_zzzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_zzzzz_0_xxxyyyyz_0[i] = 4.0 * g_zzz_0_xxxyyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyyyyz_1[i] * fz_be_0 + g_zzzz_0_xxxyyyy_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyyz_1[i] * wa_z[i];

        g_zzzzz_0_xxxyyyzz_0[i] = 4.0 * g_zzz_0_xxxyyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxxyyyz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyyzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxyyzzz_0[i] = 4.0 * g_zzz_0_xxxyyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xxxyyzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyyzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxyzzzz_0[i] = 4.0 * g_zzz_0_xxxyzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_xxxyzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxyzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxxzzzzz_0[i] = 4.0 * g_zzz_0_xxxzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxxzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_0_xxxzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxxzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxyyyyyy_0[i] = 4.0 * g_zzz_0_xxyyyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyyyyy_1[i] * fz_be_0 + g_zzzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_zzzzz_0_xxyyyyyz_0[i] = 4.0 * g_zzz_0_xxyyyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyyyyz_1[i] * fz_be_0 + g_zzzz_0_xxyyyyy_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyyz_1[i] * wa_z[i];

        g_zzzzz_0_xxyyyyzz_0[i] = 4.0 * g_zzz_0_xxyyyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xxyyyyz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyyzz_1[i] * wa_z[i];

        g_zzzzz_0_xxyyyzzz_0[i] = 4.0 * g_zzz_0_xxyyyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xxyyyzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyyzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxyyzzzz_0[i] = 4.0 * g_zzz_0_xxyyzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_xxyyzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyyzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxyzzzzz_0[i] = 4.0 * g_zzz_0_xxyzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_0_xxyzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxyzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xxzzzzzz_0[i] = 4.0 * g_zzz_0_xxzzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xxzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzz_0_xxzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xxzzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xyyyyyyy_0[i] = 4.0 * g_zzz_0_xyyyyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyyyyy_1[i] * fz_be_0 + g_zzzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_zzzzz_0_xyyyyyyz_0[i] = 4.0 * g_zzz_0_xyyyyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyyyyz_1[i] * fz_be_0 + g_zzzz_0_xyyyyyy_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyyz_1[i] * wa_z[i];

        g_zzzzz_0_xyyyyyzz_0[i] = 4.0 * g_zzz_0_xyyyyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_xyyyyyz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyyzz_1[i] * wa_z[i];

        g_zzzzz_0_xyyyyzzz_0[i] = 4.0 * g_zzz_0_xyyyyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_xyyyyzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyyzzz_1[i] * wa_z[i];

        g_zzzzz_0_xyyyzzzz_0[i] = 4.0 * g_zzz_0_xyyyzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_xyyyzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyyzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xyyzzzzz_0[i] = 4.0 * g_zzz_0_xyyzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_0_xyyzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyyzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xyzzzzzz_0[i] = 4.0 * g_zzz_0_xyzzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzz_0_xyzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xyzzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_xzzzzzzz_0[i] = 4.0 * g_zzz_0_xzzzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_xzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzzz_0_xzzzzzz_1[i] * fi_acd_0 + g_zzzz_0_xzzzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_yyyyyyyy_0[i] = 4.0 * g_zzz_0_yyyyyyyy_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyyyyy_1[i] * fz_be_0 + g_zzzz_0_yyyyyyyy_1[i] * wa_z[i];

        g_zzzzz_0_yyyyyyyz_0[i] = 4.0 * g_zzz_0_yyyyyyyz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyyyyz_1[i] * fz_be_0 + g_zzzz_0_yyyyyyy_1[i] * fi_acd_0 + g_zzzz_0_yyyyyyyz_1[i] * wa_z[i];

        g_zzzzz_0_yyyyyyzz_0[i] = 4.0 * g_zzz_0_yyyyyyzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_yyyyyyz_1[i] * fi_acd_0 + g_zzzz_0_yyyyyyzz_1[i] * wa_z[i];

        g_zzzzz_0_yyyyyzzz_0[i] = 4.0 * g_zzz_0_yyyyyzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_0_yyyyyzz_1[i] * fi_acd_0 + g_zzzz_0_yyyyyzzz_1[i] * wa_z[i];

        g_zzzzz_0_yyyyzzzz_0[i] = 4.0 * g_zzz_0_yyyyzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_0_yyyyzzz_1[i] * fi_acd_0 + g_zzzz_0_yyyyzzzz_1[i] * wa_z[i];

        g_zzzzz_0_yyyzzzzz_0[i] = 4.0 * g_zzz_0_yyyzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_0_yyyzzzz_1[i] * fi_acd_0 + g_zzzz_0_yyyzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_yyzzzzzz_0[i] = 4.0 * g_zzz_0_yyzzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzz_0_yyzzzzz_1[i] * fi_acd_0 + g_zzzz_0_yyzzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_yzzzzzzz_0[i] = 4.0 * g_zzz_0_yzzzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_yzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzzz_0_yzzzzzz_1[i] * fi_acd_0 + g_zzzz_0_yzzzzzzz_1[i] * wa_z[i];

        g_zzzzz_0_zzzzzzzz_0[i] = 4.0 * g_zzz_0_zzzzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_0_zzzzzzzz_1[i] * fz_be_0 + 8.0 * g_zzzz_0_zzzzzzz_1[i] * fi_acd_0 + g_zzzz_0_zzzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

