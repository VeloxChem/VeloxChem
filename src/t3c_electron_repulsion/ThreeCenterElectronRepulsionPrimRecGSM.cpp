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

#include "ThreeCenterElectronRepulsionPrimRecGSM.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_gsm(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsm,
                                 size_t idx_eri_0_dsm,
                                 size_t idx_eri_1_dsm,
                                 size_t idx_eri_1_fsl,
                                 size_t idx_eri_1_fsm,
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

    /// Set up components of auxilary buffer : DSM

    auto g_xx_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_dsm);

    auto g_xx_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_dsm + 1);

    auto g_xx_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_dsm + 2);

    auto g_xx_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_dsm + 3);

    auto g_xx_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_dsm + 4);

    auto g_xx_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_dsm + 5);

    auto g_xx_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_dsm + 6);

    auto g_xx_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_dsm + 7);

    auto g_xx_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_dsm + 8);

    auto g_xx_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_dsm + 9);

    auto g_xx_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_dsm + 10);

    auto g_xx_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_dsm + 11);

    auto g_xx_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_dsm + 12);

    auto g_xx_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_dsm + 13);

    auto g_xx_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_dsm + 14);

    auto g_xx_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 15);

    auto g_xx_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 16);

    auto g_xx_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 17);

    auto g_xx_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 18);

    auto g_xx_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 19);

    auto g_xx_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 20);

    auto g_xx_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 21);

    auto g_xx_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 22);

    auto g_xx_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 23);

    auto g_xx_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 24);

    auto g_xx_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 25);

    auto g_xx_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 26);

    auto g_xx_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 27);

    auto g_xx_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 28);

    auto g_xx_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 29);

    auto g_xx_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 30);

    auto g_xx_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 31);

    auto g_xx_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 32);

    auto g_xx_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 33);

    auto g_xx_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 34);

    auto g_xx_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 35);

    auto g_xx_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 36);

    auto g_xx_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 37);

    auto g_xx_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 38);

    auto g_xx_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 39);

    auto g_xx_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 40);

    auto g_xx_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 41);

    auto g_xx_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 42);

    auto g_xx_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 43);

    auto g_xx_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 44);

    auto g_xx_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 45);

    auto g_xx_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 46);

    auto g_xx_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 47);

    auto g_xx_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 48);

    auto g_xx_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 49);

    auto g_xx_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 50);

    auto g_xx_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 51);

    auto g_xx_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 52);

    auto g_xx_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 53);

    auto g_xx_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 54);

    auto g_yy_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_dsm + 165);

    auto g_yy_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_dsm + 166);

    auto g_yy_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_dsm + 167);

    auto g_yy_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_dsm + 168);

    auto g_yy_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_dsm + 169);

    auto g_yy_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_dsm + 170);

    auto g_yy_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_dsm + 171);

    auto g_yy_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_dsm + 172);

    auto g_yy_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_dsm + 173);

    auto g_yy_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_dsm + 174);

    auto g_yy_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_dsm + 175);

    auto g_yy_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_dsm + 176);

    auto g_yy_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_dsm + 177);

    auto g_yy_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_dsm + 178);

    auto g_yy_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_dsm + 179);

    auto g_yy_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 180);

    auto g_yy_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 181);

    auto g_yy_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 182);

    auto g_yy_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 183);

    auto g_yy_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 184);

    auto g_yy_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 185);

    auto g_yy_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 186);

    auto g_yy_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 187);

    auto g_yy_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 188);

    auto g_yy_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 189);

    auto g_yy_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 190);

    auto g_yy_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 191);

    auto g_yy_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 192);

    auto g_yy_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 193);

    auto g_yy_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 194);

    auto g_yy_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 195);

    auto g_yy_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 196);

    auto g_yy_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 197);

    auto g_yy_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 198);

    auto g_yy_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 199);

    auto g_yy_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 200);

    auto g_yy_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 201);

    auto g_yy_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 202);

    auto g_yy_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 203);

    auto g_yy_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 204);

    auto g_yy_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 205);

    auto g_yy_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 206);

    auto g_yy_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 207);

    auto g_yy_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 208);

    auto g_yy_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 209);

    auto g_yy_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 210);

    auto g_yy_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 211);

    auto g_yy_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 212);

    auto g_yy_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 213);

    auto g_yy_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 214);

    auto g_yy_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 215);

    auto g_yy_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 216);

    auto g_yy_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 217);

    auto g_yy_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 218);

    auto g_yy_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 219);

    auto g_zz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_dsm + 275);

    auto g_zz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_dsm + 276);

    auto g_zz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_dsm + 277);

    auto g_zz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_dsm + 278);

    auto g_zz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_dsm + 279);

    auto g_zz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_dsm + 280);

    auto g_zz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_dsm + 281);

    auto g_zz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_dsm + 282);

    auto g_zz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_dsm + 283);

    auto g_zz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_dsm + 284);

    auto g_zz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_dsm + 285);

    auto g_zz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_dsm + 286);

    auto g_zz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_dsm + 287);

    auto g_zz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_dsm + 288);

    auto g_zz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_dsm + 289);

    auto g_zz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 290);

    auto g_zz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 291);

    auto g_zz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 292);

    auto g_zz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 293);

    auto g_zz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 294);

    auto g_zz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 295);

    auto g_zz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 296);

    auto g_zz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 297);

    auto g_zz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 298);

    auto g_zz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 299);

    auto g_zz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 300);

    auto g_zz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 301);

    auto g_zz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 302);

    auto g_zz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 303);

    auto g_zz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 304);

    auto g_zz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 305);

    auto g_zz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 306);

    auto g_zz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 307);

    auto g_zz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 308);

    auto g_zz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 309);

    auto g_zz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 310);

    auto g_zz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 311);

    auto g_zz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 312);

    auto g_zz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 313);

    auto g_zz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 314);

    auto g_zz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 315);

    auto g_zz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 316);

    auto g_zz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 317);

    auto g_zz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 318);

    auto g_zz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 319);

    auto g_zz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_dsm + 320);

    auto g_zz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_dsm + 321);

    auto g_zz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_dsm + 322);

    auto g_zz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_dsm + 323);

    auto g_zz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_dsm + 324);

    auto g_zz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 325);

    auto g_zz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 326);

    auto g_zz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 327);

    auto g_zz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 328);

    auto g_zz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_dsm + 329);

    /// Set up components of auxilary buffer : DSM

    auto g_xx_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_dsm);

    auto g_xx_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_dsm + 1);

    auto g_xx_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_dsm + 2);

    auto g_xx_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_dsm + 3);

    auto g_xx_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_dsm + 4);

    auto g_xx_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_dsm + 5);

    auto g_xx_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_dsm + 6);

    auto g_xx_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_dsm + 7);

    auto g_xx_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_dsm + 8);

    auto g_xx_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_dsm + 9);

    auto g_xx_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_dsm + 10);

    auto g_xx_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_dsm + 11);

    auto g_xx_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_dsm + 12);

    auto g_xx_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_dsm + 13);

    auto g_xx_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_dsm + 14);

    auto g_xx_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 15);

    auto g_xx_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 16);

    auto g_xx_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 17);

    auto g_xx_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 18);

    auto g_xx_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 19);

    auto g_xx_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 20);

    auto g_xx_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 21);

    auto g_xx_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 22);

    auto g_xx_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 23);

    auto g_xx_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 24);

    auto g_xx_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 25);

    auto g_xx_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 26);

    auto g_xx_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 27);

    auto g_xx_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 28);

    auto g_xx_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 29);

    auto g_xx_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 30);

    auto g_xx_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 31);

    auto g_xx_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 32);

    auto g_xx_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 33);

    auto g_xx_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 34);

    auto g_xx_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 35);

    auto g_xx_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 36);

    auto g_xx_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 37);

    auto g_xx_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 38);

    auto g_xx_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 39);

    auto g_xx_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 40);

    auto g_xx_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 41);

    auto g_xx_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 42);

    auto g_xx_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 43);

    auto g_xx_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 44);

    auto g_xx_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 45);

    auto g_xx_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 46);

    auto g_xx_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 47);

    auto g_xx_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 48);

    auto g_xx_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 49);

    auto g_xx_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 50);

    auto g_xx_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 51);

    auto g_xx_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 52);

    auto g_xx_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 53);

    auto g_xx_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 54);

    auto g_yy_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_dsm + 165);

    auto g_yy_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_dsm + 166);

    auto g_yy_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_dsm + 167);

    auto g_yy_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_dsm + 168);

    auto g_yy_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_dsm + 169);

    auto g_yy_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_dsm + 170);

    auto g_yy_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_dsm + 171);

    auto g_yy_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_dsm + 172);

    auto g_yy_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_dsm + 173);

    auto g_yy_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_dsm + 174);

    auto g_yy_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_dsm + 175);

    auto g_yy_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_dsm + 176);

    auto g_yy_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_dsm + 177);

    auto g_yy_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_dsm + 178);

    auto g_yy_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_dsm + 179);

    auto g_yy_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 180);

    auto g_yy_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 181);

    auto g_yy_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 182);

    auto g_yy_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 183);

    auto g_yy_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 184);

    auto g_yy_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 185);

    auto g_yy_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 186);

    auto g_yy_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 187);

    auto g_yy_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 188);

    auto g_yy_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 189);

    auto g_yy_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 190);

    auto g_yy_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 191);

    auto g_yy_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 192);

    auto g_yy_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 193);

    auto g_yy_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 194);

    auto g_yy_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 195);

    auto g_yy_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 196);

    auto g_yy_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 197);

    auto g_yy_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 198);

    auto g_yy_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 199);

    auto g_yy_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 200);

    auto g_yy_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 201);

    auto g_yy_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 202);

    auto g_yy_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 203);

    auto g_yy_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 204);

    auto g_yy_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 205);

    auto g_yy_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 206);

    auto g_yy_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 207);

    auto g_yy_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 208);

    auto g_yy_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 209);

    auto g_yy_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 210);

    auto g_yy_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 211);

    auto g_yy_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 212);

    auto g_yy_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 213);

    auto g_yy_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 214);

    auto g_yy_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 215);

    auto g_yy_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 216);

    auto g_yy_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 217);

    auto g_yy_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 218);

    auto g_yy_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 219);

    auto g_zz_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_dsm + 275);

    auto g_zz_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_dsm + 276);

    auto g_zz_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_dsm + 277);

    auto g_zz_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_dsm + 278);

    auto g_zz_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_dsm + 279);

    auto g_zz_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_dsm + 280);

    auto g_zz_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_dsm + 281);

    auto g_zz_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_dsm + 282);

    auto g_zz_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_dsm + 283);

    auto g_zz_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_dsm + 284);

    auto g_zz_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_dsm + 285);

    auto g_zz_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_dsm + 286);

    auto g_zz_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_dsm + 287);

    auto g_zz_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_dsm + 288);

    auto g_zz_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_dsm + 289);

    auto g_zz_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 290);

    auto g_zz_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 291);

    auto g_zz_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 292);

    auto g_zz_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 293);

    auto g_zz_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 294);

    auto g_zz_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 295);

    auto g_zz_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 296);

    auto g_zz_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 297);

    auto g_zz_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 298);

    auto g_zz_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 299);

    auto g_zz_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 300);

    auto g_zz_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 301);

    auto g_zz_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 302);

    auto g_zz_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 303);

    auto g_zz_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 304);

    auto g_zz_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 305);

    auto g_zz_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 306);

    auto g_zz_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 307);

    auto g_zz_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 308);

    auto g_zz_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 309);

    auto g_zz_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 310);

    auto g_zz_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 311);

    auto g_zz_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 312);

    auto g_zz_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 313);

    auto g_zz_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 314);

    auto g_zz_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 315);

    auto g_zz_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 316);

    auto g_zz_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 317);

    auto g_zz_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 318);

    auto g_zz_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 319);

    auto g_zz_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_dsm + 320);

    auto g_zz_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_dsm + 321);

    auto g_zz_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_dsm + 322);

    auto g_zz_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_dsm + 323);

    auto g_zz_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_dsm + 324);

    auto g_zz_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 325);

    auto g_zz_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 326);

    auto g_zz_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 327);

    auto g_zz_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 328);

    auto g_zz_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_dsm + 329);

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

    auto g_xxz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 92);

    auto g_xxz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 94);

    auto g_xxz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 95);

    auto g_xxz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 97);

    auto g_xxz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 98);

    auto g_xxz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 99);

    auto g_xxz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 101);

    auto g_xxz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 102);

    auto g_xxz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 103);

    auto g_xxz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 104);

    auto g_xxz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 106);

    auto g_xxz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 107);

    auto g_xxz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 108);

    auto g_xxz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 109);

    auto g_xxz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 110);

    auto g_xxz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 112);

    auto g_xxz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 113);

    auto g_xxz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 114);

    auto g_xxz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 115);

    auto g_xxz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 116);

    auto g_xxz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 117);

    auto g_xxz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 119);

    auto g_xxz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 120);

    auto g_xxz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 121);

    auto g_xxz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 122);

    auto g_xxz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 123);

    auto g_xxz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 124);

    auto g_xxz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 125);

    auto g_xxz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 127);

    auto g_xxz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 128);

    auto g_xxz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 129);

    auto g_xxz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 130);

    auto g_xxz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 131);

    auto g_xxz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 132);

    auto g_xxz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 133);

    auto g_xxz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 134);

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

    auto g_yyz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 317);

    auto g_yyz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 319);

    auto g_yyz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 320);

    auto g_yyz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 322);

    auto g_yyz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 323);

    auto g_yyz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 324);

    auto g_yyz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 326);

    auto g_yyz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 327);

    auto g_yyz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 328);

    auto g_yyz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 329);

    auto g_yyz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 331);

    auto g_yyz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 332);

    auto g_yyz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 333);

    auto g_yyz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 334);

    auto g_yyz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 335);

    auto g_yyz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 337);

    auto g_yyz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 338);

    auto g_yyz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 339);

    auto g_yyz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 340);

    auto g_yyz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 341);

    auto g_yyz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 342);

    auto g_yyz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 344);

    auto g_yyz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 345);

    auto g_yyz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 346);

    auto g_yyz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 347);

    auto g_yyz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 348);

    auto g_yyz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 349);

    auto g_yyz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 350);

    auto g_yyz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 352);

    auto g_yyz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 353);

    auto g_yyz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 354);

    auto g_yyz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 355);

    auto g_yyz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 356);

    auto g_yyz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 357);

    auto g_yyz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 358);

    auto g_yyz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 359);

    auto g_yzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 361);

    auto g_yzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 362);

    auto g_yzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 363);

    auto g_yzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 364);

    auto g_yzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 365);

    auto g_yzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 366);

    auto g_yzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 367);

    auto g_yzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 368);

    auto g_yzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 369);

    auto g_yzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 370);

    auto g_yzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 371);

    auto g_yzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 372);

    auto g_yzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 373);

    auto g_yzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 374);

    auto g_yzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 375);

    auto g_yzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 376);

    auto g_yzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 377);

    auto g_yzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 378);

    auto g_yzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 379);

    auto g_yzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 380);

    auto g_yzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 381);

    auto g_yzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 382);

    auto g_yzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 383);

    auto g_yzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 384);

    auto g_yzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 385);

    auto g_yzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 386);

    auto g_yzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 387);

    auto g_yzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 388);

    auto g_yzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 389);

    auto g_yzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 390);

    auto g_yzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 391);

    auto g_yzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 392);

    auto g_yzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 393);

    auto g_yzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 394);

    auto g_yzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 395);

    auto g_yzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 396);

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

    /// Set up components of auxilary buffer : FSM

    auto g_xxx_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_fsm);

    auto g_xxx_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_fsm + 1);

    auto g_xxx_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_fsm + 2);

    auto g_xxx_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_fsm + 3);

    auto g_xxx_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_fsm + 4);

    auto g_xxx_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_fsm + 5);

    auto g_xxx_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_fsm + 6);

    auto g_xxx_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_fsm + 7);

    auto g_xxx_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_fsm + 8);

    auto g_xxx_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_fsm + 9);

    auto g_xxx_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_fsm + 10);

    auto g_xxx_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_fsm + 11);

    auto g_xxx_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_fsm + 12);

    auto g_xxx_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_fsm + 13);

    auto g_xxx_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_fsm + 14);

    auto g_xxx_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 15);

    auto g_xxx_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 16);

    auto g_xxx_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 17);

    auto g_xxx_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 18);

    auto g_xxx_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 19);

    auto g_xxx_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 20);

    auto g_xxx_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 21);

    auto g_xxx_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 22);

    auto g_xxx_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 23);

    auto g_xxx_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 24);

    auto g_xxx_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 25);

    auto g_xxx_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 26);

    auto g_xxx_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 27);

    auto g_xxx_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 28);

    auto g_xxx_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 29);

    auto g_xxx_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 30);

    auto g_xxx_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 31);

    auto g_xxx_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 32);

    auto g_xxx_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 33);

    auto g_xxx_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 34);

    auto g_xxx_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 35);

    auto g_xxx_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 36);

    auto g_xxx_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 37);

    auto g_xxx_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 38);

    auto g_xxx_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 39);

    auto g_xxx_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 40);

    auto g_xxx_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 41);

    auto g_xxx_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 42);

    auto g_xxx_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 43);

    auto g_xxx_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 44);

    auto g_xxx_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 45);

    auto g_xxx_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 46);

    auto g_xxx_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 47);

    auto g_xxx_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 48);

    auto g_xxx_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 49);

    auto g_xxx_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 50);

    auto g_xxx_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 51);

    auto g_xxx_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 52);

    auto g_xxx_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 53);

    auto g_xxx_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 54);

    auto g_xxy_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_fsm + 55);

    auto g_xxy_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_fsm + 56);

    auto g_xxy_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_fsm + 57);

    auto g_xxy_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_fsm + 58);

    auto g_xxy_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_fsm + 60);

    auto g_xxy_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_fsm + 61);

    auto g_xxy_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_fsm + 64);

    auto g_xxy_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_fsm + 65);

    auto g_xxy_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_fsm + 69);

    auto g_xxy_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 70);

    auto g_xxy_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 75);

    auto g_xxy_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 76);

    auto g_xxy_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 82);

    auto g_xxy_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 83);

    auto g_xxy_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 90);

    auto g_xxy_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 91);

    auto g_xxy_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 99);

    auto g_xxy_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 100);

    auto g_xxz_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_fsm + 110);

    auto g_xxz_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_fsm + 111);

    auto g_xxz_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_fsm + 112);

    auto g_xxz_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_fsm + 113);

    auto g_xxz_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_fsm + 114);

    auto g_xxz_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_fsm + 115);

    auto g_xxz_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_fsm + 116);

    auto g_xxz_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_fsm + 117);

    auto g_xxz_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_fsm + 118);

    auto g_xxz_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_fsm + 119);

    auto g_xxz_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_fsm + 120);

    auto g_xxz_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_fsm + 121);

    auto g_xxz_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_fsm + 122);

    auto g_xxz_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_fsm + 123);

    auto g_xxz_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_fsm + 124);

    auto g_xxz_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 125);

    auto g_xxz_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 126);

    auto g_xxz_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 127);

    auto g_xxz_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 128);

    auto g_xxz_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 129);

    auto g_xxz_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 130);

    auto g_xxz_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 131);

    auto g_xxz_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 132);

    auto g_xxz_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 133);

    auto g_xxz_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 134);

    auto g_xxz_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 135);

    auto g_xxz_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 136);

    auto g_xxz_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 137);

    auto g_xxz_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 138);

    auto g_xxz_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 139);

    auto g_xxz_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 140);

    auto g_xxz_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 141);

    auto g_xxz_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 142);

    auto g_xxz_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 143);

    auto g_xxz_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 144);

    auto g_xxz_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 145);

    auto g_xxz_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 146);

    auto g_xxz_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 147);

    auto g_xxz_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 148);

    auto g_xxz_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 149);

    auto g_xxz_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 150);

    auto g_xxz_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 151);

    auto g_xxz_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 152);

    auto g_xxz_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 153);

    auto g_xxz_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 154);

    auto g_xxz_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 156);

    auto g_xxz_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 157);

    auto g_xxz_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 158);

    auto g_xxz_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 159);

    auto g_xxz_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 160);

    auto g_xxz_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 161);

    auto g_xxz_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 162);

    auto g_xxz_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 163);

    auto g_xxz_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 164);

    auto g_xyy_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_fsm + 165);

    auto g_xyy_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_fsm + 166);

    auto g_xyy_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_fsm + 168);

    auto g_xyy_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_fsm + 169);

    auto g_xyy_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_fsm + 171);

    auto g_xyy_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_fsm + 172);

    auto g_xyy_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_fsm + 173);

    auto g_xyy_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_fsm + 175);

    auto g_xyy_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_fsm + 176);

    auto g_xyy_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_fsm + 177);

    auto g_xyy_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_fsm + 178);

    auto g_xyy_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 180);

    auto g_xyy_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 181);

    auto g_xyy_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 182);

    auto g_xyy_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 183);

    auto g_xyy_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 184);

    auto g_xyy_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 186);

    auto g_xyy_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 187);

    auto g_xyy_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 188);

    auto g_xyy_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 189);

    auto g_xyy_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 190);

    auto g_xyy_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 191);

    auto g_xyy_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 193);

    auto g_xyy_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 194);

    auto g_xyy_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 195);

    auto g_xyy_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 196);

    auto g_xyy_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 197);

    auto g_xyy_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 198);

    auto g_xyy_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 199);

    auto g_xyy_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 201);

    auto g_xyy_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 202);

    auto g_xyy_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 203);

    auto g_xyy_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 204);

    auto g_xyy_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 205);

    auto g_xyy_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 206);

    auto g_xyy_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 207);

    auto g_xyy_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 208);

    auto g_xyy_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 210);

    auto g_xyy_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 211);

    auto g_xyy_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 212);

    auto g_xyy_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 213);

    auto g_xyy_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 214);

    auto g_xyy_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 215);

    auto g_xyy_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 216);

    auto g_xyy_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 217);

    auto g_xyy_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 218);

    auto g_xyy_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 219);

    auto g_xzz_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_fsm + 275);

    auto g_xzz_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_fsm + 277);

    auto g_xzz_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_fsm + 279);

    auto g_xzz_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_fsm + 280);

    auto g_xzz_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_fsm + 282);

    auto g_xzz_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_fsm + 283);

    auto g_xzz_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_fsm + 284);

    auto g_xzz_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_fsm + 286);

    auto g_xzz_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_fsm + 287);

    auto g_xzz_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_fsm + 288);

    auto g_xzz_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_fsm + 289);

    auto g_xzz_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 291);

    auto g_xzz_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 292);

    auto g_xzz_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 293);

    auto g_xzz_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 294);

    auto g_xzz_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 295);

    auto g_xzz_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 297);

    auto g_xzz_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 298);

    auto g_xzz_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 299);

    auto g_xzz_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 300);

    auto g_xzz_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 301);

    auto g_xzz_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 302);

    auto g_xzz_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 304);

    auto g_xzz_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 305);

    auto g_xzz_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 306);

    auto g_xzz_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 307);

    auto g_xzz_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 308);

    auto g_xzz_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 309);

    auto g_xzz_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 310);

    auto g_xzz_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 312);

    auto g_xzz_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 313);

    auto g_xzz_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 314);

    auto g_xzz_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 315);

    auto g_xzz_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 316);

    auto g_xzz_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 317);

    auto g_xzz_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 318);

    auto g_xzz_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 319);

    auto g_xzz_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 320);

    auto g_xzz_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 321);

    auto g_xzz_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 322);

    auto g_xzz_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 323);

    auto g_xzz_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 324);

    auto g_xzz_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 325);

    auto g_xzz_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 326);

    auto g_xzz_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 327);

    auto g_xzz_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 328);

    auto g_xzz_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 329);

    auto g_yyy_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_fsm + 330);

    auto g_yyy_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_fsm + 331);

    auto g_yyy_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_fsm + 332);

    auto g_yyy_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_fsm + 333);

    auto g_yyy_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_fsm + 334);

    auto g_yyy_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_fsm + 335);

    auto g_yyy_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_fsm + 336);

    auto g_yyy_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_fsm + 337);

    auto g_yyy_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_fsm + 338);

    auto g_yyy_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_fsm + 339);

    auto g_yyy_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_fsm + 340);

    auto g_yyy_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_fsm + 341);

    auto g_yyy_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_fsm + 342);

    auto g_yyy_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_fsm + 343);

    auto g_yyy_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_fsm + 344);

    auto g_yyy_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 345);

    auto g_yyy_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 346);

    auto g_yyy_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 347);

    auto g_yyy_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 348);

    auto g_yyy_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 349);

    auto g_yyy_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 350);

    auto g_yyy_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 351);

    auto g_yyy_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 352);

    auto g_yyy_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 353);

    auto g_yyy_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 354);

    auto g_yyy_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 355);

    auto g_yyy_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 356);

    auto g_yyy_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 357);

    auto g_yyy_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 358);

    auto g_yyy_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 359);

    auto g_yyy_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 360);

    auto g_yyy_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 361);

    auto g_yyy_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 362);

    auto g_yyy_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 363);

    auto g_yyy_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 364);

    auto g_yyy_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 365);

    auto g_yyy_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 366);

    auto g_yyy_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 367);

    auto g_yyy_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 368);

    auto g_yyy_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 369);

    auto g_yyy_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 370);

    auto g_yyy_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 371);

    auto g_yyy_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 372);

    auto g_yyy_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 373);

    auto g_yyy_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 374);

    auto g_yyy_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 375);

    auto g_yyy_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 376);

    auto g_yyy_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 377);

    auto g_yyy_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 378);

    auto g_yyy_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 379);

    auto g_yyy_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 380);

    auto g_yyy_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 381);

    auto g_yyy_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 382);

    auto g_yyy_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 383);

    auto g_yyy_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 384);

    auto g_yyz_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_fsm + 386);

    auto g_yyz_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_fsm + 387);

    auto g_yyz_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_fsm + 388);

    auto g_yyz_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_fsm + 389);

    auto g_yyz_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_fsm + 390);

    auto g_yyz_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_fsm + 391);

    auto g_yyz_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_fsm + 392);

    auto g_yyz_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_fsm + 393);

    auto g_yyz_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_fsm + 394);

    auto g_yyz_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_fsm + 395);

    auto g_yyz_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_fsm + 396);

    auto g_yyz_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_fsm + 397);

    auto g_yyz_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_fsm + 398);

    auto g_yyz_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_fsm + 399);

    auto g_yyz_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 400);

    auto g_yyz_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 401);

    auto g_yyz_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 402);

    auto g_yyz_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 403);

    auto g_yyz_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 404);

    auto g_yyz_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 405);

    auto g_yyz_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 406);

    auto g_yyz_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 407);

    auto g_yyz_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 408);

    auto g_yyz_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 409);

    auto g_yyz_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 410);

    auto g_yyz_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 411);

    auto g_yyz_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 412);

    auto g_yyz_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 413);

    auto g_yyz_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 414);

    auto g_yyz_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 415);

    auto g_yyz_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 416);

    auto g_yyz_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 417);

    auto g_yyz_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 418);

    auto g_yyz_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 419);

    auto g_yyz_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 420);

    auto g_yyz_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 421);

    auto g_yyz_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 422);

    auto g_yyz_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 423);

    auto g_yyz_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 424);

    auto g_yyz_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 425);

    auto g_yyz_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 426);

    auto g_yyz_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 427);

    auto g_yyz_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 428);

    auto g_yyz_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 429);

    auto g_yyz_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 430);

    auto g_yyz_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 431);

    auto g_yyz_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 432);

    auto g_yyz_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 433);

    auto g_yyz_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 434);

    auto g_yyz_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 435);

    auto g_yyz_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 436);

    auto g_yyz_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 437);

    auto g_yyz_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 438);

    auto g_yyz_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 439);

    auto g_yzz_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_fsm + 440);

    auto g_yzz_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_fsm + 441);

    auto g_yzz_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_fsm + 442);

    auto g_yzz_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_fsm + 443);

    auto g_yzz_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_fsm + 444);

    auto g_yzz_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_fsm + 445);

    auto g_yzz_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_fsm + 446);

    auto g_yzz_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_fsm + 447);

    auto g_yzz_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_fsm + 448);

    auto g_yzz_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_fsm + 449);

    auto g_yzz_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_fsm + 450);

    auto g_yzz_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_fsm + 451);

    auto g_yzz_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_fsm + 452);

    auto g_yzz_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_fsm + 453);

    auto g_yzz_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_fsm + 454);

    auto g_yzz_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 455);

    auto g_yzz_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 456);

    auto g_yzz_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 457);

    auto g_yzz_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 458);

    auto g_yzz_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 459);

    auto g_yzz_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 460);

    auto g_yzz_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 461);

    auto g_yzz_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 462);

    auto g_yzz_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 463);

    auto g_yzz_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 464);

    auto g_yzz_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 465);

    auto g_yzz_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 466);

    auto g_yzz_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 467);

    auto g_yzz_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 468);

    auto g_yzz_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 469);

    auto g_yzz_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 470);

    auto g_yzz_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 471);

    auto g_yzz_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 472);

    auto g_yzz_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 473);

    auto g_yzz_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 474);

    auto g_yzz_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 475);

    auto g_yzz_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 476);

    auto g_yzz_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 477);

    auto g_yzz_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 478);

    auto g_yzz_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 479);

    auto g_yzz_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 480);

    auto g_yzz_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 481);

    auto g_yzz_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 482);

    auto g_yzz_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 483);

    auto g_yzz_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 484);

    auto g_yzz_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 485);

    auto g_yzz_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 486);

    auto g_yzz_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 487);

    auto g_yzz_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 488);

    auto g_yzz_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 489);

    auto g_yzz_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 490);

    auto g_yzz_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 491);

    auto g_yzz_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 492);

    auto g_yzz_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 493);

    auto g_yzz_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 494);

    auto g_zzz_0_xxxxxxxxx_1 = pbuffer.data(idx_eri_1_fsm + 495);

    auto g_zzz_0_xxxxxxxxy_1 = pbuffer.data(idx_eri_1_fsm + 496);

    auto g_zzz_0_xxxxxxxxz_1 = pbuffer.data(idx_eri_1_fsm + 497);

    auto g_zzz_0_xxxxxxxyy_1 = pbuffer.data(idx_eri_1_fsm + 498);

    auto g_zzz_0_xxxxxxxyz_1 = pbuffer.data(idx_eri_1_fsm + 499);

    auto g_zzz_0_xxxxxxxzz_1 = pbuffer.data(idx_eri_1_fsm + 500);

    auto g_zzz_0_xxxxxxyyy_1 = pbuffer.data(idx_eri_1_fsm + 501);

    auto g_zzz_0_xxxxxxyyz_1 = pbuffer.data(idx_eri_1_fsm + 502);

    auto g_zzz_0_xxxxxxyzz_1 = pbuffer.data(idx_eri_1_fsm + 503);

    auto g_zzz_0_xxxxxxzzz_1 = pbuffer.data(idx_eri_1_fsm + 504);

    auto g_zzz_0_xxxxxyyyy_1 = pbuffer.data(idx_eri_1_fsm + 505);

    auto g_zzz_0_xxxxxyyyz_1 = pbuffer.data(idx_eri_1_fsm + 506);

    auto g_zzz_0_xxxxxyyzz_1 = pbuffer.data(idx_eri_1_fsm + 507);

    auto g_zzz_0_xxxxxyzzz_1 = pbuffer.data(idx_eri_1_fsm + 508);

    auto g_zzz_0_xxxxxzzzz_1 = pbuffer.data(idx_eri_1_fsm + 509);

    auto g_zzz_0_xxxxyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 510);

    auto g_zzz_0_xxxxyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 511);

    auto g_zzz_0_xxxxyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 512);

    auto g_zzz_0_xxxxyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 513);

    auto g_zzz_0_xxxxyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 514);

    auto g_zzz_0_xxxxzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 515);

    auto g_zzz_0_xxxyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 516);

    auto g_zzz_0_xxxyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 517);

    auto g_zzz_0_xxxyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 518);

    auto g_zzz_0_xxxyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 519);

    auto g_zzz_0_xxxyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 520);

    auto g_zzz_0_xxxyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 521);

    auto g_zzz_0_xxxzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 522);

    auto g_zzz_0_xxyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 523);

    auto g_zzz_0_xxyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 524);

    auto g_zzz_0_xxyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 525);

    auto g_zzz_0_xxyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 526);

    auto g_zzz_0_xxyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 527);

    auto g_zzz_0_xxyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 528);

    auto g_zzz_0_xxyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 529);

    auto g_zzz_0_xxzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 530);

    auto g_zzz_0_xyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 531);

    auto g_zzz_0_xyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 532);

    auto g_zzz_0_xyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 533);

    auto g_zzz_0_xyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 534);

    auto g_zzz_0_xyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 535);

    auto g_zzz_0_xyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 536);

    auto g_zzz_0_xyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 537);

    auto g_zzz_0_xyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 538);

    auto g_zzz_0_xzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 539);

    auto g_zzz_0_yyyyyyyyy_1 = pbuffer.data(idx_eri_1_fsm + 540);

    auto g_zzz_0_yyyyyyyyz_1 = pbuffer.data(idx_eri_1_fsm + 541);

    auto g_zzz_0_yyyyyyyzz_1 = pbuffer.data(idx_eri_1_fsm + 542);

    auto g_zzz_0_yyyyyyzzz_1 = pbuffer.data(idx_eri_1_fsm + 543);

    auto g_zzz_0_yyyyyzzzz_1 = pbuffer.data(idx_eri_1_fsm + 544);

    auto g_zzz_0_yyyyzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 545);

    auto g_zzz_0_yyyzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 546);

    auto g_zzz_0_yyzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 547);

    auto g_zzz_0_yzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 548);

    auto g_zzz_0_zzzzzzzzz_1 = pbuffer.data(idx_eri_1_fsm + 549);

    /// Set up 0-55 components of targeted buffer : GSM

    auto g_xxxx_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm);

    auto g_xxxx_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 1);

    auto g_xxxx_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 2);

    auto g_xxxx_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 3);

    auto g_xxxx_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 4);

    auto g_xxxx_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 5);

    auto g_xxxx_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 6);

    auto g_xxxx_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 7);

    auto g_xxxx_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 8);

    auto g_xxxx_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 9);

    auto g_xxxx_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 10);

    auto g_xxxx_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 11);

    auto g_xxxx_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 12);

    auto g_xxxx_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 13);

    auto g_xxxx_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 14);

    auto g_xxxx_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 15);

    auto g_xxxx_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 16);

    auto g_xxxx_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 17);

    auto g_xxxx_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 18);

    auto g_xxxx_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 19);

    auto g_xxxx_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 20);

    auto g_xxxx_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 21);

    auto g_xxxx_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 22);

    auto g_xxxx_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 23);

    auto g_xxxx_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 24);

    auto g_xxxx_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 25);

    auto g_xxxx_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 26);

    auto g_xxxx_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 27);

    auto g_xxxx_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 28);

    auto g_xxxx_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 29);

    auto g_xxxx_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 30);

    auto g_xxxx_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 31);

    auto g_xxxx_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 32);

    auto g_xxxx_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 33);

    auto g_xxxx_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 34);

    auto g_xxxx_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 35);

    auto g_xxxx_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 36);

    auto g_xxxx_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 37);

    auto g_xxxx_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 38);

    auto g_xxxx_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 39);

    auto g_xxxx_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 40);

    auto g_xxxx_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 41);

    auto g_xxxx_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 42);

    auto g_xxxx_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 43);

    auto g_xxxx_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 44);

    auto g_xxxx_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 45);

    auto g_xxxx_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 46);

    auto g_xxxx_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 47);

    auto g_xxxx_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 48);

    auto g_xxxx_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 49);

    auto g_xxxx_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 50);

    auto g_xxxx_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 51);

    auto g_xxxx_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 52);

    auto g_xxxx_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 53);

    auto g_xxxx_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 54);

    #pragma omp simd aligned(g_xx_0_xxxxxxxxx_0, g_xx_0_xxxxxxxxx_1, g_xx_0_xxxxxxxxy_0, g_xx_0_xxxxxxxxy_1, g_xx_0_xxxxxxxxz_0, g_xx_0_xxxxxxxxz_1, g_xx_0_xxxxxxxyy_0, g_xx_0_xxxxxxxyy_1, g_xx_0_xxxxxxxyz_0, g_xx_0_xxxxxxxyz_1, g_xx_0_xxxxxxxzz_0, g_xx_0_xxxxxxxzz_1, g_xx_0_xxxxxxyyy_0, g_xx_0_xxxxxxyyy_1, g_xx_0_xxxxxxyyz_0, g_xx_0_xxxxxxyyz_1, g_xx_0_xxxxxxyzz_0, g_xx_0_xxxxxxyzz_1, g_xx_0_xxxxxxzzz_0, g_xx_0_xxxxxxzzz_1, g_xx_0_xxxxxyyyy_0, g_xx_0_xxxxxyyyy_1, g_xx_0_xxxxxyyyz_0, g_xx_0_xxxxxyyyz_1, g_xx_0_xxxxxyyzz_0, g_xx_0_xxxxxyyzz_1, g_xx_0_xxxxxyzzz_0, g_xx_0_xxxxxyzzz_1, g_xx_0_xxxxxzzzz_0, g_xx_0_xxxxxzzzz_1, g_xx_0_xxxxyyyyy_0, g_xx_0_xxxxyyyyy_1, g_xx_0_xxxxyyyyz_0, g_xx_0_xxxxyyyyz_1, g_xx_0_xxxxyyyzz_0, g_xx_0_xxxxyyyzz_1, g_xx_0_xxxxyyzzz_0, g_xx_0_xxxxyyzzz_1, g_xx_0_xxxxyzzzz_0, g_xx_0_xxxxyzzzz_1, g_xx_0_xxxxzzzzz_0, g_xx_0_xxxxzzzzz_1, g_xx_0_xxxyyyyyy_0, g_xx_0_xxxyyyyyy_1, g_xx_0_xxxyyyyyz_0, g_xx_0_xxxyyyyyz_1, g_xx_0_xxxyyyyzz_0, g_xx_0_xxxyyyyzz_1, g_xx_0_xxxyyyzzz_0, g_xx_0_xxxyyyzzz_1, g_xx_0_xxxyyzzzz_0, g_xx_0_xxxyyzzzz_1, g_xx_0_xxxyzzzzz_0, g_xx_0_xxxyzzzzz_1, g_xx_0_xxxzzzzzz_0, g_xx_0_xxxzzzzzz_1, g_xx_0_xxyyyyyyy_0, g_xx_0_xxyyyyyyy_1, g_xx_0_xxyyyyyyz_0, g_xx_0_xxyyyyyyz_1, g_xx_0_xxyyyyyzz_0, g_xx_0_xxyyyyyzz_1, g_xx_0_xxyyyyzzz_0, g_xx_0_xxyyyyzzz_1, g_xx_0_xxyyyzzzz_0, g_xx_0_xxyyyzzzz_1, g_xx_0_xxyyzzzzz_0, g_xx_0_xxyyzzzzz_1, g_xx_0_xxyzzzzzz_0, g_xx_0_xxyzzzzzz_1, g_xx_0_xxzzzzzzz_0, g_xx_0_xxzzzzzzz_1, g_xx_0_xyyyyyyyy_0, g_xx_0_xyyyyyyyy_1, g_xx_0_xyyyyyyyz_0, g_xx_0_xyyyyyyyz_1, g_xx_0_xyyyyyyzz_0, g_xx_0_xyyyyyyzz_1, g_xx_0_xyyyyyzzz_0, g_xx_0_xyyyyyzzz_1, g_xx_0_xyyyyzzzz_0, g_xx_0_xyyyyzzzz_1, g_xx_0_xyyyzzzzz_0, g_xx_0_xyyyzzzzz_1, g_xx_0_xyyzzzzzz_0, g_xx_0_xyyzzzzzz_1, g_xx_0_xyzzzzzzz_0, g_xx_0_xyzzzzzzz_1, g_xx_0_xzzzzzzzz_0, g_xx_0_xzzzzzzzz_1, g_xx_0_yyyyyyyyy_0, g_xx_0_yyyyyyyyy_1, g_xx_0_yyyyyyyyz_0, g_xx_0_yyyyyyyyz_1, g_xx_0_yyyyyyyzz_0, g_xx_0_yyyyyyyzz_1, g_xx_0_yyyyyyzzz_0, g_xx_0_yyyyyyzzz_1, g_xx_0_yyyyyzzzz_0, g_xx_0_yyyyyzzzz_1, g_xx_0_yyyyzzzzz_0, g_xx_0_yyyyzzzzz_1, g_xx_0_yyyzzzzzz_0, g_xx_0_yyyzzzzzz_1, g_xx_0_yyzzzzzzz_0, g_xx_0_yyzzzzzzz_1, g_xx_0_yzzzzzzzz_0, g_xx_0_yzzzzzzzz_1, g_xx_0_zzzzzzzzz_0, g_xx_0_zzzzzzzzz_1, g_xxx_0_xxxxxxxx_1, g_xxx_0_xxxxxxxxx_1, g_xxx_0_xxxxxxxxy_1, g_xxx_0_xxxxxxxxz_1, g_xxx_0_xxxxxxxy_1, g_xxx_0_xxxxxxxyy_1, g_xxx_0_xxxxxxxyz_1, g_xxx_0_xxxxxxxz_1, g_xxx_0_xxxxxxxzz_1, g_xxx_0_xxxxxxyy_1, g_xxx_0_xxxxxxyyy_1, g_xxx_0_xxxxxxyyz_1, g_xxx_0_xxxxxxyz_1, g_xxx_0_xxxxxxyzz_1, g_xxx_0_xxxxxxzz_1, g_xxx_0_xxxxxxzzz_1, g_xxx_0_xxxxxyyy_1, g_xxx_0_xxxxxyyyy_1, g_xxx_0_xxxxxyyyz_1, g_xxx_0_xxxxxyyz_1, g_xxx_0_xxxxxyyzz_1, g_xxx_0_xxxxxyzz_1, g_xxx_0_xxxxxyzzz_1, g_xxx_0_xxxxxzzz_1, g_xxx_0_xxxxxzzzz_1, g_xxx_0_xxxxyyyy_1, g_xxx_0_xxxxyyyyy_1, g_xxx_0_xxxxyyyyz_1, g_xxx_0_xxxxyyyz_1, g_xxx_0_xxxxyyyzz_1, g_xxx_0_xxxxyyzz_1, g_xxx_0_xxxxyyzzz_1, g_xxx_0_xxxxyzzz_1, g_xxx_0_xxxxyzzzz_1, g_xxx_0_xxxxzzzz_1, g_xxx_0_xxxxzzzzz_1, g_xxx_0_xxxyyyyy_1, g_xxx_0_xxxyyyyyy_1, g_xxx_0_xxxyyyyyz_1, g_xxx_0_xxxyyyyz_1, g_xxx_0_xxxyyyyzz_1, g_xxx_0_xxxyyyzz_1, g_xxx_0_xxxyyyzzz_1, g_xxx_0_xxxyyzzz_1, g_xxx_0_xxxyyzzzz_1, g_xxx_0_xxxyzzzz_1, g_xxx_0_xxxyzzzzz_1, g_xxx_0_xxxzzzzz_1, g_xxx_0_xxxzzzzzz_1, g_xxx_0_xxyyyyyy_1, g_xxx_0_xxyyyyyyy_1, g_xxx_0_xxyyyyyyz_1, g_xxx_0_xxyyyyyz_1, g_xxx_0_xxyyyyyzz_1, g_xxx_0_xxyyyyzz_1, g_xxx_0_xxyyyyzzz_1, g_xxx_0_xxyyyzzz_1, g_xxx_0_xxyyyzzzz_1, g_xxx_0_xxyyzzzz_1, g_xxx_0_xxyyzzzzz_1, g_xxx_0_xxyzzzzz_1, g_xxx_0_xxyzzzzzz_1, g_xxx_0_xxzzzzzz_1, g_xxx_0_xxzzzzzzz_1, g_xxx_0_xyyyyyyy_1, g_xxx_0_xyyyyyyyy_1, g_xxx_0_xyyyyyyyz_1, g_xxx_0_xyyyyyyz_1, g_xxx_0_xyyyyyyzz_1, g_xxx_0_xyyyyyzz_1, g_xxx_0_xyyyyyzzz_1, g_xxx_0_xyyyyzzz_1, g_xxx_0_xyyyyzzzz_1, g_xxx_0_xyyyzzzz_1, g_xxx_0_xyyyzzzzz_1, g_xxx_0_xyyzzzzz_1, g_xxx_0_xyyzzzzzz_1, g_xxx_0_xyzzzzzz_1, g_xxx_0_xyzzzzzzz_1, g_xxx_0_xzzzzzzz_1, g_xxx_0_xzzzzzzzz_1, g_xxx_0_yyyyyyyy_1, g_xxx_0_yyyyyyyyy_1, g_xxx_0_yyyyyyyyz_1, g_xxx_0_yyyyyyyz_1, g_xxx_0_yyyyyyyzz_1, g_xxx_0_yyyyyyzz_1, g_xxx_0_yyyyyyzzz_1, g_xxx_0_yyyyyzzz_1, g_xxx_0_yyyyyzzzz_1, g_xxx_0_yyyyzzzz_1, g_xxx_0_yyyyzzzzz_1, g_xxx_0_yyyzzzzz_1, g_xxx_0_yyyzzzzzz_1, g_xxx_0_yyzzzzzz_1, g_xxx_0_yyzzzzzzz_1, g_xxx_0_yzzzzzzz_1, g_xxx_0_yzzzzzzzz_1, g_xxx_0_zzzzzzzz_1, g_xxx_0_zzzzzzzzz_1, g_xxxx_0_xxxxxxxxx_0, g_xxxx_0_xxxxxxxxy_0, g_xxxx_0_xxxxxxxxz_0, g_xxxx_0_xxxxxxxyy_0, g_xxxx_0_xxxxxxxyz_0, g_xxxx_0_xxxxxxxzz_0, g_xxxx_0_xxxxxxyyy_0, g_xxxx_0_xxxxxxyyz_0, g_xxxx_0_xxxxxxyzz_0, g_xxxx_0_xxxxxxzzz_0, g_xxxx_0_xxxxxyyyy_0, g_xxxx_0_xxxxxyyyz_0, g_xxxx_0_xxxxxyyzz_0, g_xxxx_0_xxxxxyzzz_0, g_xxxx_0_xxxxxzzzz_0, g_xxxx_0_xxxxyyyyy_0, g_xxxx_0_xxxxyyyyz_0, g_xxxx_0_xxxxyyyzz_0, g_xxxx_0_xxxxyyzzz_0, g_xxxx_0_xxxxyzzzz_0, g_xxxx_0_xxxxzzzzz_0, g_xxxx_0_xxxyyyyyy_0, g_xxxx_0_xxxyyyyyz_0, g_xxxx_0_xxxyyyyzz_0, g_xxxx_0_xxxyyyzzz_0, g_xxxx_0_xxxyyzzzz_0, g_xxxx_0_xxxyzzzzz_0, g_xxxx_0_xxxzzzzzz_0, g_xxxx_0_xxyyyyyyy_0, g_xxxx_0_xxyyyyyyz_0, g_xxxx_0_xxyyyyyzz_0, g_xxxx_0_xxyyyyzzz_0, g_xxxx_0_xxyyyzzzz_0, g_xxxx_0_xxyyzzzzz_0, g_xxxx_0_xxyzzzzzz_0, g_xxxx_0_xxzzzzzzz_0, g_xxxx_0_xyyyyyyyy_0, g_xxxx_0_xyyyyyyyz_0, g_xxxx_0_xyyyyyyzz_0, g_xxxx_0_xyyyyyzzz_0, g_xxxx_0_xyyyyzzzz_0, g_xxxx_0_xyyyzzzzz_0, g_xxxx_0_xyyzzzzzz_0, g_xxxx_0_xyzzzzzzz_0, g_xxxx_0_xzzzzzzzz_0, g_xxxx_0_yyyyyyyyy_0, g_xxxx_0_yyyyyyyyz_0, g_xxxx_0_yyyyyyyzz_0, g_xxxx_0_yyyyyyzzz_0, g_xxxx_0_yyyyyzzzz_0, g_xxxx_0_yyyyzzzzz_0, g_xxxx_0_yyyzzzzzz_0, g_xxxx_0_yyzzzzzzz_0, g_xxxx_0_yzzzzzzzz_0, g_xxxx_0_zzzzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxx_0_xxxxxxxxx_0[i] = 3.0 * g_xx_0_xxxxxxxxx_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxxxx_1[i] * fz_be_0 + 9.0 * g_xxx_0_xxxxxxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxxx_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxxxy_0[i] = 3.0 * g_xx_0_xxxxxxxxy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxxxy_1[i] * fz_be_0 + 8.0 * g_xxx_0_xxxxxxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxxy_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxxxz_0[i] = 3.0 * g_xx_0_xxxxxxxxz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxxxz_1[i] * fz_be_0 + 8.0 * g_xxx_0_xxxxxxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxxz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxxyy_0[i] = 3.0 * g_xx_0_xxxxxxxyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxxyy_1[i] * fz_be_0 + 7.0 * g_xxx_0_xxxxxxyy_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxyy_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxxyz_0[i] = 3.0 * g_xx_0_xxxxxxxyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxxyz_1[i] * fz_be_0 + 7.0 * g_xxx_0_xxxxxxyz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxxzz_0[i] = 3.0 * g_xx_0_xxxxxxxzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxxzz_1[i] * fz_be_0 + 7.0 * g_xxx_0_xxxxxxzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxyyy_0[i] = 3.0 * g_xx_0_xxxxxxyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxyyy_1[i] * fz_be_0 + 6.0 * g_xxx_0_xxxxxyyy_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyyy_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxyyz_0[i] = 3.0 * g_xx_0_xxxxxxyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxyyz_1[i] * fz_be_0 + 6.0 * g_xxx_0_xxxxxyyz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxyzz_0[i] = 3.0 * g_xx_0_xxxxxxyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxyzz_1[i] * fz_be_0 + 6.0 * g_xxx_0_xxxxxyzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxzzz_0[i] = 3.0 * g_xx_0_xxxxxxzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxzzz_1[i] * fz_be_0 + 6.0 * g_xxx_0_xxxxxzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxyyyy_0[i] = 3.0 * g_xx_0_xxxxxyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxyyyy_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxyyyy_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyyy_1[i] * wa_x[i];

        g_xxxx_0_xxxxxyyyz_0[i] = 3.0 * g_xx_0_xxxxxyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxyyyz_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxyyyz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxyyzz_0[i] = 3.0 * g_xx_0_xxxxxyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxyyzz_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxyyzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxyzzz_0[i] = 3.0 * g_xx_0_xxxxxyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxyzzz_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxyzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxzzzz_0[i] = 3.0 * g_xx_0_xxxxxzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxzzzz_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxyyyyy_0[i] = 3.0 * g_xx_0_xxxxyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyyyyy_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyyyyy_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyyy_1[i] * wa_x[i];

        g_xxxx_0_xxxxyyyyz_0[i] = 3.0 * g_xx_0_xxxxyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyyyyz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyyyyz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xxxx_0_xxxxyyyzz_0[i] = 3.0 * g_xx_0_xxxxyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyyyzz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyyyzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxyyzzz_0[i] = 3.0 * g_xx_0_xxxxyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyyzzz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyyzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxyzzzz_0[i] = 3.0 * g_xx_0_xxxxyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyzzzz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxzzzzz_0[i] = 3.0 * g_xx_0_xxxxzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxzzzzz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxyyyyyy_0[i] = 3.0 * g_xx_0_xxxyyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyyyyy_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyyyyy_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyyy_1[i] * wa_x[i];

        g_xxxx_0_xxxyyyyyz_0[i] = 3.0 * g_xx_0_xxxyyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyyyyz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyyyyz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xxxx_0_xxxyyyyzz_0[i] = 3.0 * g_xx_0_xxxyyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyyyzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyyyzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xxxx_0_xxxyyyzzz_0[i] = 3.0 * g_xx_0_xxxyyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyyzzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyyzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxyyzzzz_0[i] = 3.0 * g_xx_0_xxxyyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyzzzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxyzzzzz_0[i] = 3.0 * g_xx_0_xxxyzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyzzzzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxzzzzzz_0[i] = 3.0 * g_xx_0_xxxzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxzzzzzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxyyyyyyy_0[i] = 3.0 * g_xx_0_xxyyyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyyyyy_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyyyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyyy_1[i] * wa_x[i];

        g_xxxx_0_xxyyyyyyz_0[i] = 3.0 * g_xx_0_xxyyyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyyyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xxxx_0_xxyyyyyzz_0[i] = 3.0 * g_xx_0_xxyyyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyyyzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xxxx_0_xxyyyyzzz_0[i] = 3.0 * g_xx_0_xxyyyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyyzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xxxx_0_xxyyyzzzz_0[i] = 3.0 * g_xx_0_xxyyyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxyyzzzzz_0[i] = 3.0 * g_xx_0_xxyyzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxyzzzzzz_0[i] = 3.0 * g_xx_0_xxyzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyzzzzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxzzzzzzz_0[i] = 3.0 * g_xx_0_xxzzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxzzzzzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xzzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxzzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xyyyyyyyy_0[i] = 3.0 * g_xx_0_xyyyyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyyyyy_1[i] * fz_be_0 + g_xxx_0_yyyyyyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyyy_1[i] * wa_x[i];

        g_xxxx_0_xyyyyyyyz_0[i] = 3.0 * g_xx_0_xyyyyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyyyyz_1[i] * fz_be_0 + g_xxx_0_yyyyyyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xxxx_0_xyyyyyyzz_0[i] = 3.0 * g_xx_0_xyyyyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyyyzz_1[i] * fz_be_0 + g_xxx_0_yyyyyyzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xxxx_0_xyyyyyzzz_0[i] = 3.0 * g_xx_0_xyyyyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyyzzz_1[i] * fz_be_0 + g_xxx_0_yyyyyzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xxxx_0_xyyyyzzzz_0[i] = 3.0 * g_xx_0_xyyyyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyzzzz_1[i] * fz_be_0 + g_xxx_0_yyyyzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xxxx_0_xyyyzzzzz_0[i] = 3.0 * g_xx_0_xyyyzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyzzzzz_1[i] * fz_be_0 + g_xxx_0_yyyzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xyyzzzzzz_0[i] = 3.0 * g_xx_0_xyyzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyzzzzzz_1[i] * fz_be_0 + g_xxx_0_yyzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xyzzzzzzz_0[i] = 3.0 * g_xx_0_xyzzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyzzzzzzz_1[i] * fz_be_0 + g_xxx_0_yzzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xzzzzzzzz_0[i] = 3.0 * g_xx_0_xzzzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xzzzzzzzz_1[i] * fz_be_0 + g_xxx_0_zzzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xzzzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_yyyyyyyyy_0[i] = 3.0 * g_xx_0_yyyyyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyyyyy_1[i] * fz_be_0 + g_xxx_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xxxx_0_yyyyyyyyz_0[i] = 3.0 * g_xx_0_yyyyyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyyyyz_1[i] * fz_be_0 + g_xxx_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xxxx_0_yyyyyyyzz_0[i] = 3.0 * g_xx_0_yyyyyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyyyzz_1[i] * fz_be_0 + g_xxx_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xxxx_0_yyyyyyzzz_0[i] = 3.0 * g_xx_0_yyyyyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyyzzz_1[i] * fz_be_0 + g_xxx_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xxxx_0_yyyyyzzzz_0[i] = 3.0 * g_xx_0_yyyyyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyzzzz_1[i] * fz_be_0 + g_xxx_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xxxx_0_yyyyzzzzz_0[i] = 3.0 * g_xx_0_yyyyzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyzzzzz_1[i] * fz_be_0 + g_xxx_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xxxx_0_yyyzzzzzz_0[i] = 3.0 * g_xx_0_yyyzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyzzzzzz_1[i] * fz_be_0 + g_xxx_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_yyzzzzzzz_0[i] = 3.0 * g_xx_0_yyzzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyzzzzzzz_1[i] * fz_be_0 + g_xxx_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_yzzzzzzzz_0[i] = 3.0 * g_xx_0_yzzzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yzzzzzzzz_1[i] * fz_be_0 + g_xxx_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_zzzzzzzzz_0[i] = 3.0 * g_xx_0_zzzzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_zzzzzzzzz_1[i] * fz_be_0 + g_xxx_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 55-110 components of targeted buffer : GSM

    auto g_xxxy_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 55);

    auto g_xxxy_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 56);

    auto g_xxxy_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 57);

    auto g_xxxy_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 58);

    auto g_xxxy_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 59);

    auto g_xxxy_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 60);

    auto g_xxxy_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 61);

    auto g_xxxy_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 62);

    auto g_xxxy_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 63);

    auto g_xxxy_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 64);

    auto g_xxxy_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 65);

    auto g_xxxy_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 66);

    auto g_xxxy_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 67);

    auto g_xxxy_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 68);

    auto g_xxxy_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 69);

    auto g_xxxy_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 70);

    auto g_xxxy_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 71);

    auto g_xxxy_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 72);

    auto g_xxxy_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 73);

    auto g_xxxy_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 74);

    auto g_xxxy_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 75);

    auto g_xxxy_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 76);

    auto g_xxxy_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 77);

    auto g_xxxy_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 78);

    auto g_xxxy_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 79);

    auto g_xxxy_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 80);

    auto g_xxxy_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 81);

    auto g_xxxy_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 82);

    auto g_xxxy_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 83);

    auto g_xxxy_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 84);

    auto g_xxxy_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 85);

    auto g_xxxy_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 86);

    auto g_xxxy_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 87);

    auto g_xxxy_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 88);

    auto g_xxxy_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 89);

    auto g_xxxy_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 90);

    auto g_xxxy_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 91);

    auto g_xxxy_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 92);

    auto g_xxxy_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 93);

    auto g_xxxy_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 94);

    auto g_xxxy_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 95);

    auto g_xxxy_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 96);

    auto g_xxxy_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 97);

    auto g_xxxy_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 98);

    auto g_xxxy_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 99);

    auto g_xxxy_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 100);

    auto g_xxxy_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 101);

    auto g_xxxy_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 102);

    auto g_xxxy_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 103);

    auto g_xxxy_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 104);

    auto g_xxxy_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 105);

    auto g_xxxy_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 106);

    auto g_xxxy_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 107);

    auto g_xxxy_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 108);

    auto g_xxxy_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 109);

    #pragma omp simd aligned(g_xxx_0_xxxxxxxx_1, g_xxx_0_xxxxxxxxx_1, g_xxx_0_xxxxxxxxy_1, g_xxx_0_xxxxxxxxz_1, g_xxx_0_xxxxxxxy_1, g_xxx_0_xxxxxxxyy_1, g_xxx_0_xxxxxxxyz_1, g_xxx_0_xxxxxxxz_1, g_xxx_0_xxxxxxxzz_1, g_xxx_0_xxxxxxyy_1, g_xxx_0_xxxxxxyyy_1, g_xxx_0_xxxxxxyyz_1, g_xxx_0_xxxxxxyz_1, g_xxx_0_xxxxxxyzz_1, g_xxx_0_xxxxxxzz_1, g_xxx_0_xxxxxxzzz_1, g_xxx_0_xxxxxyyy_1, g_xxx_0_xxxxxyyyy_1, g_xxx_0_xxxxxyyyz_1, g_xxx_0_xxxxxyyz_1, g_xxx_0_xxxxxyyzz_1, g_xxx_0_xxxxxyzz_1, g_xxx_0_xxxxxyzzz_1, g_xxx_0_xxxxxzzz_1, g_xxx_0_xxxxxzzzz_1, g_xxx_0_xxxxyyyy_1, g_xxx_0_xxxxyyyyy_1, g_xxx_0_xxxxyyyyz_1, g_xxx_0_xxxxyyyz_1, g_xxx_0_xxxxyyyzz_1, g_xxx_0_xxxxyyzz_1, g_xxx_0_xxxxyyzzz_1, g_xxx_0_xxxxyzzz_1, g_xxx_0_xxxxyzzzz_1, g_xxx_0_xxxxzzzz_1, g_xxx_0_xxxxzzzzz_1, g_xxx_0_xxxyyyyy_1, g_xxx_0_xxxyyyyyy_1, g_xxx_0_xxxyyyyyz_1, g_xxx_0_xxxyyyyz_1, g_xxx_0_xxxyyyyzz_1, g_xxx_0_xxxyyyzz_1, g_xxx_0_xxxyyyzzz_1, g_xxx_0_xxxyyzzz_1, g_xxx_0_xxxyyzzzz_1, g_xxx_0_xxxyzzzz_1, g_xxx_0_xxxyzzzzz_1, g_xxx_0_xxxzzzzz_1, g_xxx_0_xxxzzzzzz_1, g_xxx_0_xxyyyyyy_1, g_xxx_0_xxyyyyyyy_1, g_xxx_0_xxyyyyyyz_1, g_xxx_0_xxyyyyyz_1, g_xxx_0_xxyyyyyzz_1, g_xxx_0_xxyyyyzz_1, g_xxx_0_xxyyyyzzz_1, g_xxx_0_xxyyyzzz_1, g_xxx_0_xxyyyzzzz_1, g_xxx_0_xxyyzzzz_1, g_xxx_0_xxyyzzzzz_1, g_xxx_0_xxyzzzzz_1, g_xxx_0_xxyzzzzzz_1, g_xxx_0_xxzzzzzz_1, g_xxx_0_xxzzzzzzz_1, g_xxx_0_xyyyyyyy_1, g_xxx_0_xyyyyyyyy_1, g_xxx_0_xyyyyyyyz_1, g_xxx_0_xyyyyyyz_1, g_xxx_0_xyyyyyyzz_1, g_xxx_0_xyyyyyzz_1, g_xxx_0_xyyyyyzzz_1, g_xxx_0_xyyyyzzz_1, g_xxx_0_xyyyyzzzz_1, g_xxx_0_xyyyzzzz_1, g_xxx_0_xyyyzzzzz_1, g_xxx_0_xyyzzzzz_1, g_xxx_0_xyyzzzzzz_1, g_xxx_0_xyzzzzzz_1, g_xxx_0_xyzzzzzzz_1, g_xxx_0_xzzzzzzz_1, g_xxx_0_xzzzzzzzz_1, g_xxx_0_yyyyyyyy_1, g_xxx_0_yyyyyyyyy_1, g_xxx_0_yyyyyyyyz_1, g_xxx_0_yyyyyyyz_1, g_xxx_0_yyyyyyyzz_1, g_xxx_0_yyyyyyzz_1, g_xxx_0_yyyyyyzzz_1, g_xxx_0_yyyyyzzz_1, g_xxx_0_yyyyyzzzz_1, g_xxx_0_yyyyzzzz_1, g_xxx_0_yyyyzzzzz_1, g_xxx_0_yyyzzzzz_1, g_xxx_0_yyyzzzzzz_1, g_xxx_0_yyzzzzzz_1, g_xxx_0_yyzzzzzzz_1, g_xxx_0_yzzzzzzz_1, g_xxx_0_yzzzzzzzz_1, g_xxx_0_zzzzzzzz_1, g_xxx_0_zzzzzzzzz_1, g_xxxy_0_xxxxxxxxx_0, g_xxxy_0_xxxxxxxxy_0, g_xxxy_0_xxxxxxxxz_0, g_xxxy_0_xxxxxxxyy_0, g_xxxy_0_xxxxxxxyz_0, g_xxxy_0_xxxxxxxzz_0, g_xxxy_0_xxxxxxyyy_0, g_xxxy_0_xxxxxxyyz_0, g_xxxy_0_xxxxxxyzz_0, g_xxxy_0_xxxxxxzzz_0, g_xxxy_0_xxxxxyyyy_0, g_xxxy_0_xxxxxyyyz_0, g_xxxy_0_xxxxxyyzz_0, g_xxxy_0_xxxxxyzzz_0, g_xxxy_0_xxxxxzzzz_0, g_xxxy_0_xxxxyyyyy_0, g_xxxy_0_xxxxyyyyz_0, g_xxxy_0_xxxxyyyzz_0, g_xxxy_0_xxxxyyzzz_0, g_xxxy_0_xxxxyzzzz_0, g_xxxy_0_xxxxzzzzz_0, g_xxxy_0_xxxyyyyyy_0, g_xxxy_0_xxxyyyyyz_0, g_xxxy_0_xxxyyyyzz_0, g_xxxy_0_xxxyyyzzz_0, g_xxxy_0_xxxyyzzzz_0, g_xxxy_0_xxxyzzzzz_0, g_xxxy_0_xxxzzzzzz_0, g_xxxy_0_xxyyyyyyy_0, g_xxxy_0_xxyyyyyyz_0, g_xxxy_0_xxyyyyyzz_0, g_xxxy_0_xxyyyyzzz_0, g_xxxy_0_xxyyyzzzz_0, g_xxxy_0_xxyyzzzzz_0, g_xxxy_0_xxyzzzzzz_0, g_xxxy_0_xxzzzzzzz_0, g_xxxy_0_xyyyyyyyy_0, g_xxxy_0_xyyyyyyyz_0, g_xxxy_0_xyyyyyyzz_0, g_xxxy_0_xyyyyyzzz_0, g_xxxy_0_xyyyyzzzz_0, g_xxxy_0_xyyyzzzzz_0, g_xxxy_0_xyyzzzzzz_0, g_xxxy_0_xyzzzzzzz_0, g_xxxy_0_xzzzzzzzz_0, g_xxxy_0_yyyyyyyyy_0, g_xxxy_0_yyyyyyyyz_0, g_xxxy_0_yyyyyyyzz_0, g_xxxy_0_yyyyyyzzz_0, g_xxxy_0_yyyyyzzzz_0, g_xxxy_0_yyyyzzzzz_0, g_xxxy_0_yyyzzzzzz_0, g_xxxy_0_yyzzzzzzz_0, g_xxxy_0_yzzzzzzzz_0, g_xxxy_0_zzzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxy_0_xxxxxxxxx_0[i] = g_xxx_0_xxxxxxxxx_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxxxy_0[i] = g_xxx_0_xxxxxxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxxy_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxxxz_0[i] = g_xxx_0_xxxxxxxxz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxxyy_0[i] = 2.0 * g_xxx_0_xxxxxxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxyy_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxxyz_0[i] = g_xxx_0_xxxxxxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxyz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxxzz_0[i] = g_xxx_0_xxxxxxxzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxyyy_0[i] = 3.0 * g_xxx_0_xxxxxxyy_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyyy_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxyyz_0[i] = 2.0 * g_xxx_0_xxxxxxyz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyyz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxyzz_0[i] = g_xxx_0_xxxxxxzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxzzz_0[i] = g_xxx_0_xxxxxxzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxyyyy_0[i] = 4.0 * g_xxx_0_xxxxxyyy_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyyy_1[i] * wa_y[i];

        g_xxxy_0_xxxxxyyyz_0[i] = 3.0 * g_xxx_0_xxxxxyyz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyyz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxyyzz_0[i] = 2.0 * g_xxx_0_xxxxxyzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxyzzz_0[i] = g_xxx_0_xxxxxzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxzzzz_0[i] = g_xxx_0_xxxxxzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxyyyyy_0[i] = 5.0 * g_xxx_0_xxxxyyyy_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyyy_1[i] * wa_y[i];

        g_xxxy_0_xxxxyyyyz_0[i] = 4.0 * g_xxx_0_xxxxyyyz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyyz_1[i] * wa_y[i];

        g_xxxy_0_xxxxyyyzz_0[i] = 3.0 * g_xxx_0_xxxxyyzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxyyzzz_0[i] = 2.0 * g_xxx_0_xxxxyzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxyzzzz_0[i] = g_xxx_0_xxxxzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxzzzzz_0[i] = g_xxx_0_xxxxzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxyyyyyy_0[i] = 6.0 * g_xxx_0_xxxyyyyy_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyyy_1[i] * wa_y[i];

        g_xxxy_0_xxxyyyyyz_0[i] = 5.0 * g_xxx_0_xxxyyyyz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyyz_1[i] * wa_y[i];

        g_xxxy_0_xxxyyyyzz_0[i] = 4.0 * g_xxx_0_xxxyyyzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyzz_1[i] * wa_y[i];

        g_xxxy_0_xxxyyyzzz_0[i] = 3.0 * g_xxx_0_xxxyyzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxyyzzzz_0[i] = 2.0 * g_xxx_0_xxxyzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxyzzzzz_0[i] = g_xxx_0_xxxzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxzzzzzz_0[i] = g_xxx_0_xxxzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxyyyyyyy_0[i] = 7.0 * g_xxx_0_xxyyyyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyyy_1[i] * wa_y[i];

        g_xxxy_0_xxyyyyyyz_0[i] = 6.0 * g_xxx_0_xxyyyyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyyz_1[i] * wa_y[i];

        g_xxxy_0_xxyyyyyzz_0[i] = 5.0 * g_xxx_0_xxyyyyzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyzz_1[i] * wa_y[i];

        g_xxxy_0_xxyyyyzzz_0[i] = 4.0 * g_xxx_0_xxyyyzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyzzz_1[i] * wa_y[i];

        g_xxxy_0_xxyyyzzzz_0[i] = 3.0 * g_xxx_0_xxyyzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxyyzzzzz_0[i] = 2.0 * g_xxx_0_xxyzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxyzzzzzz_0[i] = g_xxx_0_xxzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxzzzzzzz_0[i] = g_xxx_0_xxzzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xyyyyyyyy_0[i] = 8.0 * g_xxx_0_xyyyyyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyyy_1[i] * wa_y[i];

        g_xxxy_0_xyyyyyyyz_0[i] = 7.0 * g_xxx_0_xyyyyyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyyz_1[i] * wa_y[i];

        g_xxxy_0_xyyyyyyzz_0[i] = 6.0 * g_xxx_0_xyyyyyzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyzz_1[i] * wa_y[i];

        g_xxxy_0_xyyyyyzzz_0[i] = 5.0 * g_xxx_0_xyyyyzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyzzz_1[i] * wa_y[i];

        g_xxxy_0_xyyyyzzzz_0[i] = 4.0 * g_xxx_0_xyyyzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyzzzz_1[i] * wa_y[i];

        g_xxxy_0_xyyyzzzzz_0[i] = 3.0 * g_xxx_0_xyyzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xyyzzzzzz_0[i] = 2.0 * g_xxx_0_xyzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xyzzzzzzz_0[i] = g_xxx_0_xzzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xzzzzzzzz_0[i] = g_xxx_0_xzzzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_yyyyyyyyy_0[i] = 9.0 * g_xxx_0_yyyyyyyy_1[i] * fi_acd_0 + g_xxx_0_yyyyyyyyy_1[i] * wa_y[i];

        g_xxxy_0_yyyyyyyyz_0[i] = 8.0 * g_xxx_0_yyyyyyyz_1[i] * fi_acd_0 + g_xxx_0_yyyyyyyyz_1[i] * wa_y[i];

        g_xxxy_0_yyyyyyyzz_0[i] = 7.0 * g_xxx_0_yyyyyyzz_1[i] * fi_acd_0 + g_xxx_0_yyyyyyyzz_1[i] * wa_y[i];

        g_xxxy_0_yyyyyyzzz_0[i] = 6.0 * g_xxx_0_yyyyyzzz_1[i] * fi_acd_0 + g_xxx_0_yyyyyyzzz_1[i] * wa_y[i];

        g_xxxy_0_yyyyyzzzz_0[i] = 5.0 * g_xxx_0_yyyyzzzz_1[i] * fi_acd_0 + g_xxx_0_yyyyyzzzz_1[i] * wa_y[i];

        g_xxxy_0_yyyyzzzzz_0[i] = 4.0 * g_xxx_0_yyyzzzzz_1[i] * fi_acd_0 + g_xxx_0_yyyyzzzzz_1[i] * wa_y[i];

        g_xxxy_0_yyyzzzzzz_0[i] = 3.0 * g_xxx_0_yyzzzzzz_1[i] * fi_acd_0 + g_xxx_0_yyyzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_yyzzzzzzz_0[i] = 2.0 * g_xxx_0_yzzzzzzz_1[i] * fi_acd_0 + g_xxx_0_yyzzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_yzzzzzzzz_0[i] = g_xxx_0_zzzzzzzz_1[i] * fi_acd_0 + g_xxx_0_yzzzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_zzzzzzzzz_0[i] = g_xxx_0_zzzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 110-165 components of targeted buffer : GSM

    auto g_xxxz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 110);

    auto g_xxxz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 111);

    auto g_xxxz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 112);

    auto g_xxxz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 113);

    auto g_xxxz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 114);

    auto g_xxxz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 115);

    auto g_xxxz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 116);

    auto g_xxxz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 117);

    auto g_xxxz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 118);

    auto g_xxxz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 119);

    auto g_xxxz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 120);

    auto g_xxxz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 121);

    auto g_xxxz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 122);

    auto g_xxxz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 123);

    auto g_xxxz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 124);

    auto g_xxxz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 125);

    auto g_xxxz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 126);

    auto g_xxxz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 127);

    auto g_xxxz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 128);

    auto g_xxxz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 129);

    auto g_xxxz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 130);

    auto g_xxxz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 131);

    auto g_xxxz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 132);

    auto g_xxxz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 133);

    auto g_xxxz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 134);

    auto g_xxxz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 135);

    auto g_xxxz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 136);

    auto g_xxxz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 137);

    auto g_xxxz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 138);

    auto g_xxxz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 139);

    auto g_xxxz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 140);

    auto g_xxxz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 141);

    auto g_xxxz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 142);

    auto g_xxxz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 143);

    auto g_xxxz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 144);

    auto g_xxxz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 145);

    auto g_xxxz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 146);

    auto g_xxxz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 147);

    auto g_xxxz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 148);

    auto g_xxxz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 149);

    auto g_xxxz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 150);

    auto g_xxxz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 151);

    auto g_xxxz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 152);

    auto g_xxxz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 153);

    auto g_xxxz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 154);

    auto g_xxxz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 155);

    auto g_xxxz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 156);

    auto g_xxxz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 157);

    auto g_xxxz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 158);

    auto g_xxxz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 159);

    auto g_xxxz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 160);

    auto g_xxxz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 161);

    auto g_xxxz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 162);

    auto g_xxxz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 163);

    auto g_xxxz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 164);

    #pragma omp simd aligned(g_xxx_0_xxxxxxxx_1, g_xxx_0_xxxxxxxxx_1, g_xxx_0_xxxxxxxxy_1, g_xxx_0_xxxxxxxxz_1, g_xxx_0_xxxxxxxy_1, g_xxx_0_xxxxxxxyy_1, g_xxx_0_xxxxxxxyz_1, g_xxx_0_xxxxxxxz_1, g_xxx_0_xxxxxxxzz_1, g_xxx_0_xxxxxxyy_1, g_xxx_0_xxxxxxyyy_1, g_xxx_0_xxxxxxyyz_1, g_xxx_0_xxxxxxyz_1, g_xxx_0_xxxxxxyzz_1, g_xxx_0_xxxxxxzz_1, g_xxx_0_xxxxxxzzz_1, g_xxx_0_xxxxxyyy_1, g_xxx_0_xxxxxyyyy_1, g_xxx_0_xxxxxyyyz_1, g_xxx_0_xxxxxyyz_1, g_xxx_0_xxxxxyyzz_1, g_xxx_0_xxxxxyzz_1, g_xxx_0_xxxxxyzzz_1, g_xxx_0_xxxxxzzz_1, g_xxx_0_xxxxxzzzz_1, g_xxx_0_xxxxyyyy_1, g_xxx_0_xxxxyyyyy_1, g_xxx_0_xxxxyyyyz_1, g_xxx_0_xxxxyyyz_1, g_xxx_0_xxxxyyyzz_1, g_xxx_0_xxxxyyzz_1, g_xxx_0_xxxxyyzzz_1, g_xxx_0_xxxxyzzz_1, g_xxx_0_xxxxyzzzz_1, g_xxx_0_xxxxzzzz_1, g_xxx_0_xxxxzzzzz_1, g_xxx_0_xxxyyyyy_1, g_xxx_0_xxxyyyyyy_1, g_xxx_0_xxxyyyyyz_1, g_xxx_0_xxxyyyyz_1, g_xxx_0_xxxyyyyzz_1, g_xxx_0_xxxyyyzz_1, g_xxx_0_xxxyyyzzz_1, g_xxx_0_xxxyyzzz_1, g_xxx_0_xxxyyzzzz_1, g_xxx_0_xxxyzzzz_1, g_xxx_0_xxxyzzzzz_1, g_xxx_0_xxxzzzzz_1, g_xxx_0_xxxzzzzzz_1, g_xxx_0_xxyyyyyy_1, g_xxx_0_xxyyyyyyy_1, g_xxx_0_xxyyyyyyz_1, g_xxx_0_xxyyyyyz_1, g_xxx_0_xxyyyyyzz_1, g_xxx_0_xxyyyyzz_1, g_xxx_0_xxyyyyzzz_1, g_xxx_0_xxyyyzzz_1, g_xxx_0_xxyyyzzzz_1, g_xxx_0_xxyyzzzz_1, g_xxx_0_xxyyzzzzz_1, g_xxx_0_xxyzzzzz_1, g_xxx_0_xxyzzzzzz_1, g_xxx_0_xxzzzzzz_1, g_xxx_0_xxzzzzzzz_1, g_xxx_0_xyyyyyyy_1, g_xxx_0_xyyyyyyyy_1, g_xxx_0_xyyyyyyyz_1, g_xxx_0_xyyyyyyz_1, g_xxx_0_xyyyyyyzz_1, g_xxx_0_xyyyyyzz_1, g_xxx_0_xyyyyyzzz_1, g_xxx_0_xyyyyzzz_1, g_xxx_0_xyyyyzzzz_1, g_xxx_0_xyyyzzzz_1, g_xxx_0_xyyyzzzzz_1, g_xxx_0_xyyzzzzz_1, g_xxx_0_xyyzzzzzz_1, g_xxx_0_xyzzzzzz_1, g_xxx_0_xyzzzzzzz_1, g_xxx_0_xzzzzzzz_1, g_xxx_0_xzzzzzzzz_1, g_xxx_0_yyyyyyyy_1, g_xxx_0_yyyyyyyyy_1, g_xxx_0_yyyyyyyyz_1, g_xxx_0_yyyyyyyz_1, g_xxx_0_yyyyyyyzz_1, g_xxx_0_yyyyyyzz_1, g_xxx_0_yyyyyyzzz_1, g_xxx_0_yyyyyzzz_1, g_xxx_0_yyyyyzzzz_1, g_xxx_0_yyyyzzzz_1, g_xxx_0_yyyyzzzzz_1, g_xxx_0_yyyzzzzz_1, g_xxx_0_yyyzzzzzz_1, g_xxx_0_yyzzzzzz_1, g_xxx_0_yyzzzzzzz_1, g_xxx_0_yzzzzzzz_1, g_xxx_0_yzzzzzzzz_1, g_xxx_0_zzzzzzzz_1, g_xxx_0_zzzzzzzzz_1, g_xxxz_0_xxxxxxxxx_0, g_xxxz_0_xxxxxxxxy_0, g_xxxz_0_xxxxxxxxz_0, g_xxxz_0_xxxxxxxyy_0, g_xxxz_0_xxxxxxxyz_0, g_xxxz_0_xxxxxxxzz_0, g_xxxz_0_xxxxxxyyy_0, g_xxxz_0_xxxxxxyyz_0, g_xxxz_0_xxxxxxyzz_0, g_xxxz_0_xxxxxxzzz_0, g_xxxz_0_xxxxxyyyy_0, g_xxxz_0_xxxxxyyyz_0, g_xxxz_0_xxxxxyyzz_0, g_xxxz_0_xxxxxyzzz_0, g_xxxz_0_xxxxxzzzz_0, g_xxxz_0_xxxxyyyyy_0, g_xxxz_0_xxxxyyyyz_0, g_xxxz_0_xxxxyyyzz_0, g_xxxz_0_xxxxyyzzz_0, g_xxxz_0_xxxxyzzzz_0, g_xxxz_0_xxxxzzzzz_0, g_xxxz_0_xxxyyyyyy_0, g_xxxz_0_xxxyyyyyz_0, g_xxxz_0_xxxyyyyzz_0, g_xxxz_0_xxxyyyzzz_0, g_xxxz_0_xxxyyzzzz_0, g_xxxz_0_xxxyzzzzz_0, g_xxxz_0_xxxzzzzzz_0, g_xxxz_0_xxyyyyyyy_0, g_xxxz_0_xxyyyyyyz_0, g_xxxz_0_xxyyyyyzz_0, g_xxxz_0_xxyyyyzzz_0, g_xxxz_0_xxyyyzzzz_0, g_xxxz_0_xxyyzzzzz_0, g_xxxz_0_xxyzzzzzz_0, g_xxxz_0_xxzzzzzzz_0, g_xxxz_0_xyyyyyyyy_0, g_xxxz_0_xyyyyyyyz_0, g_xxxz_0_xyyyyyyzz_0, g_xxxz_0_xyyyyyzzz_0, g_xxxz_0_xyyyyzzzz_0, g_xxxz_0_xyyyzzzzz_0, g_xxxz_0_xyyzzzzzz_0, g_xxxz_0_xyzzzzzzz_0, g_xxxz_0_xzzzzzzzz_0, g_xxxz_0_yyyyyyyyy_0, g_xxxz_0_yyyyyyyyz_0, g_xxxz_0_yyyyyyyzz_0, g_xxxz_0_yyyyyyzzz_0, g_xxxz_0_yyyyyzzzz_0, g_xxxz_0_yyyyzzzzz_0, g_xxxz_0_yyyzzzzzz_0, g_xxxz_0_yyzzzzzzz_0, g_xxxz_0_yzzzzzzzz_0, g_xxxz_0_zzzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxz_0_xxxxxxxxx_0[i] = g_xxx_0_xxxxxxxxx_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxxxy_0[i] = g_xxx_0_xxxxxxxxy_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxxxz_0[i] = g_xxx_0_xxxxxxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxxz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxxyy_0[i] = g_xxx_0_xxxxxxxyy_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxxyz_0[i] = g_xxx_0_xxxxxxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxyz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxxzz_0[i] = 2.0 * g_xxx_0_xxxxxxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxyyy_0[i] = g_xxx_0_xxxxxxyyy_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxyyz_0[i] = g_xxx_0_xxxxxxyy_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyyz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxyzz_0[i] = 2.0 * g_xxx_0_xxxxxxyz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxzzz_0[i] = 3.0 * g_xxx_0_xxxxxxzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxyyyy_0[i] = g_xxx_0_xxxxxyyyy_1[i] * wa_z[i];

        g_xxxz_0_xxxxxyyyz_0[i] = g_xxx_0_xxxxxyyy_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyyz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxyyzz_0[i] = 2.0 * g_xxx_0_xxxxxyyz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxyzzz_0[i] = 3.0 * g_xxx_0_xxxxxyzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxzzzz_0[i] = 4.0 * g_xxx_0_xxxxxzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxyyyyy_0[i] = g_xxx_0_xxxxyyyyy_1[i] * wa_z[i];

        g_xxxz_0_xxxxyyyyz_0[i] = g_xxx_0_xxxxyyyy_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyyz_1[i] * wa_z[i];

        g_xxxz_0_xxxxyyyzz_0[i] = 2.0 * g_xxx_0_xxxxyyyz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxyyzzz_0[i] = 3.0 * g_xxx_0_xxxxyyzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxyzzzz_0[i] = 4.0 * g_xxx_0_xxxxyzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxzzzzz_0[i] = 5.0 * g_xxx_0_xxxxzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxyyyyyy_0[i] = g_xxx_0_xxxyyyyyy_1[i] * wa_z[i];

        g_xxxz_0_xxxyyyyyz_0[i] = g_xxx_0_xxxyyyyy_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyyz_1[i] * wa_z[i];

        g_xxxz_0_xxxyyyyzz_0[i] = 2.0 * g_xxx_0_xxxyyyyz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyzz_1[i] * wa_z[i];

        g_xxxz_0_xxxyyyzzz_0[i] = 3.0 * g_xxx_0_xxxyyyzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxyyzzzz_0[i] = 4.0 * g_xxx_0_xxxyyzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxyzzzzz_0[i] = 5.0 * g_xxx_0_xxxyzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxzzzzzz_0[i] = 6.0 * g_xxx_0_xxxzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxyyyyyyy_0[i] = g_xxx_0_xxyyyyyyy_1[i] * wa_z[i];

        g_xxxz_0_xxyyyyyyz_0[i] = g_xxx_0_xxyyyyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyyz_1[i] * wa_z[i];

        g_xxxz_0_xxyyyyyzz_0[i] = 2.0 * g_xxx_0_xxyyyyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyzz_1[i] * wa_z[i];

        g_xxxz_0_xxyyyyzzz_0[i] = 3.0 * g_xxx_0_xxyyyyzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyzzz_1[i] * wa_z[i];

        g_xxxz_0_xxyyyzzzz_0[i] = 4.0 * g_xxx_0_xxyyyzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxyyzzzzz_0[i] = 5.0 * g_xxx_0_xxyyzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxyzzzzzz_0[i] = 6.0 * g_xxx_0_xxyzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxzzzzzzz_0[i] = 7.0 * g_xxx_0_xxzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxzzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xyyyyyyyy_0[i] = g_xxx_0_xyyyyyyyy_1[i] * wa_z[i];

        g_xxxz_0_xyyyyyyyz_0[i] = g_xxx_0_xyyyyyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyyz_1[i] * wa_z[i];

        g_xxxz_0_xyyyyyyzz_0[i] = 2.0 * g_xxx_0_xyyyyyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyzz_1[i] * wa_z[i];

        g_xxxz_0_xyyyyyzzz_0[i] = 3.0 * g_xxx_0_xyyyyyzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyzzz_1[i] * wa_z[i];

        g_xxxz_0_xyyyyzzzz_0[i] = 4.0 * g_xxx_0_xyyyyzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyzzzz_1[i] * wa_z[i];

        g_xxxz_0_xyyyzzzzz_0[i] = 5.0 * g_xxx_0_xyyyzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xyyzzzzzz_0[i] = 6.0 * g_xxx_0_xyyzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xyzzzzzzz_0[i] = 7.0 * g_xxx_0_xyzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xzzzzzzzz_0[i] = 8.0 * g_xxx_0_xzzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xzzzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_yyyyyyyyy_0[i] = g_xxx_0_yyyyyyyyy_1[i] * wa_z[i];

        g_xxxz_0_yyyyyyyyz_0[i] = g_xxx_0_yyyyyyyy_1[i] * fi_acd_0 + g_xxx_0_yyyyyyyyz_1[i] * wa_z[i];

        g_xxxz_0_yyyyyyyzz_0[i] = 2.0 * g_xxx_0_yyyyyyyz_1[i] * fi_acd_0 + g_xxx_0_yyyyyyyzz_1[i] * wa_z[i];

        g_xxxz_0_yyyyyyzzz_0[i] = 3.0 * g_xxx_0_yyyyyyzz_1[i] * fi_acd_0 + g_xxx_0_yyyyyyzzz_1[i] * wa_z[i];

        g_xxxz_0_yyyyyzzzz_0[i] = 4.0 * g_xxx_0_yyyyyzzz_1[i] * fi_acd_0 + g_xxx_0_yyyyyzzzz_1[i] * wa_z[i];

        g_xxxz_0_yyyyzzzzz_0[i] = 5.0 * g_xxx_0_yyyyzzzz_1[i] * fi_acd_0 + g_xxx_0_yyyyzzzzz_1[i] * wa_z[i];

        g_xxxz_0_yyyzzzzzz_0[i] = 6.0 * g_xxx_0_yyyzzzzz_1[i] * fi_acd_0 + g_xxx_0_yyyzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_yyzzzzzzz_0[i] = 7.0 * g_xxx_0_yyzzzzzz_1[i] * fi_acd_0 + g_xxx_0_yyzzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_yzzzzzzzz_0[i] = 8.0 * g_xxx_0_yzzzzzzz_1[i] * fi_acd_0 + g_xxx_0_yzzzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_zzzzzzzzz_0[i] = 9.0 * g_xxx_0_zzzzzzzz_1[i] * fi_acd_0 + g_xxx_0_zzzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 165-220 components of targeted buffer : GSM

    auto g_xxyy_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 165);

    auto g_xxyy_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 166);

    auto g_xxyy_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 167);

    auto g_xxyy_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 168);

    auto g_xxyy_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 169);

    auto g_xxyy_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 170);

    auto g_xxyy_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 171);

    auto g_xxyy_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 172);

    auto g_xxyy_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 173);

    auto g_xxyy_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 174);

    auto g_xxyy_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 175);

    auto g_xxyy_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 176);

    auto g_xxyy_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 177);

    auto g_xxyy_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 178);

    auto g_xxyy_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 179);

    auto g_xxyy_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 180);

    auto g_xxyy_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 181);

    auto g_xxyy_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 182);

    auto g_xxyy_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 183);

    auto g_xxyy_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 184);

    auto g_xxyy_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 185);

    auto g_xxyy_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 186);

    auto g_xxyy_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 187);

    auto g_xxyy_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 188);

    auto g_xxyy_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 189);

    auto g_xxyy_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 190);

    auto g_xxyy_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 191);

    auto g_xxyy_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 192);

    auto g_xxyy_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 193);

    auto g_xxyy_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 194);

    auto g_xxyy_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 195);

    auto g_xxyy_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 196);

    auto g_xxyy_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 197);

    auto g_xxyy_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 198);

    auto g_xxyy_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 199);

    auto g_xxyy_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 200);

    auto g_xxyy_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 201);

    auto g_xxyy_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 202);

    auto g_xxyy_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 203);

    auto g_xxyy_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 204);

    auto g_xxyy_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 205);

    auto g_xxyy_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 206);

    auto g_xxyy_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 207);

    auto g_xxyy_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 208);

    auto g_xxyy_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 209);

    auto g_xxyy_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 210);

    auto g_xxyy_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 211);

    auto g_xxyy_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 212);

    auto g_xxyy_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 213);

    auto g_xxyy_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 214);

    auto g_xxyy_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 215);

    auto g_xxyy_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 216);

    auto g_xxyy_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 217);

    auto g_xxyy_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 218);

    auto g_xxyy_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 219);

    #pragma omp simd aligned(g_xx_0_xxxxxxxxx_0, g_xx_0_xxxxxxxxx_1, g_xx_0_xxxxxxxxz_0, g_xx_0_xxxxxxxxz_1, g_xx_0_xxxxxxxzz_0, g_xx_0_xxxxxxxzz_1, g_xx_0_xxxxxxzzz_0, g_xx_0_xxxxxxzzz_1, g_xx_0_xxxxxzzzz_0, g_xx_0_xxxxxzzzz_1, g_xx_0_xxxxzzzzz_0, g_xx_0_xxxxzzzzz_1, g_xx_0_xxxzzzzzz_0, g_xx_0_xxxzzzzzz_1, g_xx_0_xxzzzzzzz_0, g_xx_0_xxzzzzzzz_1, g_xx_0_xzzzzzzzz_0, g_xx_0_xzzzzzzzz_1, g_xxy_0_xxxxxxxxx_1, g_xxy_0_xxxxxxxxz_1, g_xxy_0_xxxxxxxzz_1, g_xxy_0_xxxxxxzzz_1, g_xxy_0_xxxxxzzzz_1, g_xxy_0_xxxxzzzzz_1, g_xxy_0_xxxzzzzzz_1, g_xxy_0_xxzzzzzzz_1, g_xxy_0_xzzzzzzzz_1, g_xxyy_0_xxxxxxxxx_0, g_xxyy_0_xxxxxxxxy_0, g_xxyy_0_xxxxxxxxz_0, g_xxyy_0_xxxxxxxyy_0, g_xxyy_0_xxxxxxxyz_0, g_xxyy_0_xxxxxxxzz_0, g_xxyy_0_xxxxxxyyy_0, g_xxyy_0_xxxxxxyyz_0, g_xxyy_0_xxxxxxyzz_0, g_xxyy_0_xxxxxxzzz_0, g_xxyy_0_xxxxxyyyy_0, g_xxyy_0_xxxxxyyyz_0, g_xxyy_0_xxxxxyyzz_0, g_xxyy_0_xxxxxyzzz_0, g_xxyy_0_xxxxxzzzz_0, g_xxyy_0_xxxxyyyyy_0, g_xxyy_0_xxxxyyyyz_0, g_xxyy_0_xxxxyyyzz_0, g_xxyy_0_xxxxyyzzz_0, g_xxyy_0_xxxxyzzzz_0, g_xxyy_0_xxxxzzzzz_0, g_xxyy_0_xxxyyyyyy_0, g_xxyy_0_xxxyyyyyz_0, g_xxyy_0_xxxyyyyzz_0, g_xxyy_0_xxxyyyzzz_0, g_xxyy_0_xxxyyzzzz_0, g_xxyy_0_xxxyzzzzz_0, g_xxyy_0_xxxzzzzzz_0, g_xxyy_0_xxyyyyyyy_0, g_xxyy_0_xxyyyyyyz_0, g_xxyy_0_xxyyyyyzz_0, g_xxyy_0_xxyyyyzzz_0, g_xxyy_0_xxyyyzzzz_0, g_xxyy_0_xxyyzzzzz_0, g_xxyy_0_xxyzzzzzz_0, g_xxyy_0_xxzzzzzzz_0, g_xxyy_0_xyyyyyyyy_0, g_xxyy_0_xyyyyyyyz_0, g_xxyy_0_xyyyyyyzz_0, g_xxyy_0_xyyyyyzzz_0, g_xxyy_0_xyyyyzzzz_0, g_xxyy_0_xyyyzzzzz_0, g_xxyy_0_xyyzzzzzz_0, g_xxyy_0_xyzzzzzzz_0, g_xxyy_0_xzzzzzzzz_0, g_xxyy_0_yyyyyyyyy_0, g_xxyy_0_yyyyyyyyz_0, g_xxyy_0_yyyyyyyzz_0, g_xxyy_0_yyyyyyzzz_0, g_xxyy_0_yyyyyzzzz_0, g_xxyy_0_yyyyzzzzz_0, g_xxyy_0_yyyzzzzzz_0, g_xxyy_0_yyzzzzzzz_0, g_xxyy_0_yzzzzzzzz_0, g_xxyy_0_zzzzzzzzz_0, g_xyy_0_xxxxxxxxy_1, g_xyy_0_xxxxxxxy_1, g_xyy_0_xxxxxxxyy_1, g_xyy_0_xxxxxxxyz_1, g_xyy_0_xxxxxxyy_1, g_xyy_0_xxxxxxyyy_1, g_xyy_0_xxxxxxyyz_1, g_xyy_0_xxxxxxyz_1, g_xyy_0_xxxxxxyzz_1, g_xyy_0_xxxxxyyy_1, g_xyy_0_xxxxxyyyy_1, g_xyy_0_xxxxxyyyz_1, g_xyy_0_xxxxxyyz_1, g_xyy_0_xxxxxyyzz_1, g_xyy_0_xxxxxyzz_1, g_xyy_0_xxxxxyzzz_1, g_xyy_0_xxxxyyyy_1, g_xyy_0_xxxxyyyyy_1, g_xyy_0_xxxxyyyyz_1, g_xyy_0_xxxxyyyz_1, g_xyy_0_xxxxyyyzz_1, g_xyy_0_xxxxyyzz_1, g_xyy_0_xxxxyyzzz_1, g_xyy_0_xxxxyzzz_1, g_xyy_0_xxxxyzzzz_1, g_xyy_0_xxxyyyyy_1, g_xyy_0_xxxyyyyyy_1, g_xyy_0_xxxyyyyyz_1, g_xyy_0_xxxyyyyz_1, g_xyy_0_xxxyyyyzz_1, g_xyy_0_xxxyyyzz_1, g_xyy_0_xxxyyyzzz_1, g_xyy_0_xxxyyzzz_1, g_xyy_0_xxxyyzzzz_1, g_xyy_0_xxxyzzzz_1, g_xyy_0_xxxyzzzzz_1, g_xyy_0_xxyyyyyy_1, g_xyy_0_xxyyyyyyy_1, g_xyy_0_xxyyyyyyz_1, g_xyy_0_xxyyyyyz_1, g_xyy_0_xxyyyyyzz_1, g_xyy_0_xxyyyyzz_1, g_xyy_0_xxyyyyzzz_1, g_xyy_0_xxyyyzzz_1, g_xyy_0_xxyyyzzzz_1, g_xyy_0_xxyyzzzz_1, g_xyy_0_xxyyzzzzz_1, g_xyy_0_xxyzzzzz_1, g_xyy_0_xxyzzzzzz_1, g_xyy_0_xyyyyyyy_1, g_xyy_0_xyyyyyyyy_1, g_xyy_0_xyyyyyyyz_1, g_xyy_0_xyyyyyyz_1, g_xyy_0_xyyyyyyzz_1, g_xyy_0_xyyyyyzz_1, g_xyy_0_xyyyyyzzz_1, g_xyy_0_xyyyyzzz_1, g_xyy_0_xyyyyzzzz_1, g_xyy_0_xyyyzzzz_1, g_xyy_0_xyyyzzzzz_1, g_xyy_0_xyyzzzzz_1, g_xyy_0_xyyzzzzzz_1, g_xyy_0_xyzzzzzz_1, g_xyy_0_xyzzzzzzz_1, g_xyy_0_yyyyyyyy_1, g_xyy_0_yyyyyyyyy_1, g_xyy_0_yyyyyyyyz_1, g_xyy_0_yyyyyyyz_1, g_xyy_0_yyyyyyyzz_1, g_xyy_0_yyyyyyzz_1, g_xyy_0_yyyyyyzzz_1, g_xyy_0_yyyyyzzz_1, g_xyy_0_yyyyyzzzz_1, g_xyy_0_yyyyzzzz_1, g_xyy_0_yyyyzzzzz_1, g_xyy_0_yyyzzzzz_1, g_xyy_0_yyyzzzzzz_1, g_xyy_0_yyzzzzzz_1, g_xyy_0_yyzzzzzzz_1, g_xyy_0_yzzzzzzz_1, g_xyy_0_yzzzzzzzz_1, g_xyy_0_zzzzzzzzz_1, g_yy_0_xxxxxxxxy_0, g_yy_0_xxxxxxxxy_1, g_yy_0_xxxxxxxyy_0, g_yy_0_xxxxxxxyy_1, g_yy_0_xxxxxxxyz_0, g_yy_0_xxxxxxxyz_1, g_yy_0_xxxxxxyyy_0, g_yy_0_xxxxxxyyy_1, g_yy_0_xxxxxxyyz_0, g_yy_0_xxxxxxyyz_1, g_yy_0_xxxxxxyzz_0, g_yy_0_xxxxxxyzz_1, g_yy_0_xxxxxyyyy_0, g_yy_0_xxxxxyyyy_1, g_yy_0_xxxxxyyyz_0, g_yy_0_xxxxxyyyz_1, g_yy_0_xxxxxyyzz_0, g_yy_0_xxxxxyyzz_1, g_yy_0_xxxxxyzzz_0, g_yy_0_xxxxxyzzz_1, g_yy_0_xxxxyyyyy_0, g_yy_0_xxxxyyyyy_1, g_yy_0_xxxxyyyyz_0, g_yy_0_xxxxyyyyz_1, g_yy_0_xxxxyyyzz_0, g_yy_0_xxxxyyyzz_1, g_yy_0_xxxxyyzzz_0, g_yy_0_xxxxyyzzz_1, g_yy_0_xxxxyzzzz_0, g_yy_0_xxxxyzzzz_1, g_yy_0_xxxyyyyyy_0, g_yy_0_xxxyyyyyy_1, g_yy_0_xxxyyyyyz_0, g_yy_0_xxxyyyyyz_1, g_yy_0_xxxyyyyzz_0, g_yy_0_xxxyyyyzz_1, g_yy_0_xxxyyyzzz_0, g_yy_0_xxxyyyzzz_1, g_yy_0_xxxyyzzzz_0, g_yy_0_xxxyyzzzz_1, g_yy_0_xxxyzzzzz_0, g_yy_0_xxxyzzzzz_1, g_yy_0_xxyyyyyyy_0, g_yy_0_xxyyyyyyy_1, g_yy_0_xxyyyyyyz_0, g_yy_0_xxyyyyyyz_1, g_yy_0_xxyyyyyzz_0, g_yy_0_xxyyyyyzz_1, g_yy_0_xxyyyyzzz_0, g_yy_0_xxyyyyzzz_1, g_yy_0_xxyyyzzzz_0, g_yy_0_xxyyyzzzz_1, g_yy_0_xxyyzzzzz_0, g_yy_0_xxyyzzzzz_1, g_yy_0_xxyzzzzzz_0, g_yy_0_xxyzzzzzz_1, g_yy_0_xyyyyyyyy_0, g_yy_0_xyyyyyyyy_1, g_yy_0_xyyyyyyyz_0, g_yy_0_xyyyyyyyz_1, g_yy_0_xyyyyyyzz_0, g_yy_0_xyyyyyyzz_1, g_yy_0_xyyyyyzzz_0, g_yy_0_xyyyyyzzz_1, g_yy_0_xyyyyzzzz_0, g_yy_0_xyyyyzzzz_1, g_yy_0_xyyyzzzzz_0, g_yy_0_xyyyzzzzz_1, g_yy_0_xyyzzzzzz_0, g_yy_0_xyyzzzzzz_1, g_yy_0_xyzzzzzzz_0, g_yy_0_xyzzzzzzz_1, g_yy_0_yyyyyyyyy_0, g_yy_0_yyyyyyyyy_1, g_yy_0_yyyyyyyyz_0, g_yy_0_yyyyyyyyz_1, g_yy_0_yyyyyyyzz_0, g_yy_0_yyyyyyyzz_1, g_yy_0_yyyyyyzzz_0, g_yy_0_yyyyyyzzz_1, g_yy_0_yyyyyzzzz_0, g_yy_0_yyyyyzzzz_1, g_yy_0_yyyyzzzzz_0, g_yy_0_yyyyzzzzz_1, g_yy_0_yyyzzzzzz_0, g_yy_0_yyyzzzzzz_1, g_yy_0_yyzzzzzzz_0, g_yy_0_yyzzzzzzz_1, g_yy_0_yzzzzzzzz_0, g_yy_0_yzzzzzzzz_1, g_yy_0_zzzzzzzzz_0, g_yy_0_zzzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyy_0_xxxxxxxxx_0[i] = g_xx_0_xxxxxxxxx_0[i] * fbe_0 - g_xx_0_xxxxxxxxx_1[i] * fz_be_0 + g_xxy_0_xxxxxxxxx_1[i] * wa_y[i];

        g_xxyy_0_xxxxxxxxy_0[i] = g_yy_0_xxxxxxxxy_0[i] * fbe_0 - g_yy_0_xxxxxxxxy_1[i] * fz_be_0 + 8.0 * g_xyy_0_xxxxxxxy_1[i] * fi_acd_0 + g_xyy_0_xxxxxxxxy_1[i] * wa_x[i];

        g_xxyy_0_xxxxxxxxz_0[i] = g_xx_0_xxxxxxxxz_0[i] * fbe_0 - g_xx_0_xxxxxxxxz_1[i] * fz_be_0 + g_xxy_0_xxxxxxxxz_1[i] * wa_y[i];

        g_xxyy_0_xxxxxxxyy_0[i] = g_yy_0_xxxxxxxyy_0[i] * fbe_0 - g_yy_0_xxxxxxxyy_1[i] * fz_be_0 + 7.0 * g_xyy_0_xxxxxxyy_1[i] * fi_acd_0 + g_xyy_0_xxxxxxxyy_1[i] * wa_x[i];

        g_xxyy_0_xxxxxxxyz_0[i] = g_yy_0_xxxxxxxyz_0[i] * fbe_0 - g_yy_0_xxxxxxxyz_1[i] * fz_be_0 + 7.0 * g_xyy_0_xxxxxxyz_1[i] * fi_acd_0 + g_xyy_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xxyy_0_xxxxxxxzz_0[i] = g_xx_0_xxxxxxxzz_0[i] * fbe_0 - g_xx_0_xxxxxxxzz_1[i] * fz_be_0 + g_xxy_0_xxxxxxxzz_1[i] * wa_y[i];

        g_xxyy_0_xxxxxxyyy_0[i] = g_yy_0_xxxxxxyyy_0[i] * fbe_0 - g_yy_0_xxxxxxyyy_1[i] * fz_be_0 + 6.0 * g_xyy_0_xxxxxyyy_1[i] * fi_acd_0 + g_xyy_0_xxxxxxyyy_1[i] * wa_x[i];

        g_xxyy_0_xxxxxxyyz_0[i] = g_yy_0_xxxxxxyyz_0[i] * fbe_0 - g_yy_0_xxxxxxyyz_1[i] * fz_be_0 + 6.0 * g_xyy_0_xxxxxyyz_1[i] * fi_acd_0 + g_xyy_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xxyy_0_xxxxxxyzz_0[i] = g_yy_0_xxxxxxyzz_0[i] * fbe_0 - g_yy_0_xxxxxxyzz_1[i] * fz_be_0 + 6.0 * g_xyy_0_xxxxxyzz_1[i] * fi_acd_0 + g_xyy_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xxyy_0_xxxxxxzzz_0[i] = g_xx_0_xxxxxxzzz_0[i] * fbe_0 - g_xx_0_xxxxxxzzz_1[i] * fz_be_0 + g_xxy_0_xxxxxxzzz_1[i] * wa_y[i];

        g_xxyy_0_xxxxxyyyy_0[i] = g_yy_0_xxxxxyyyy_0[i] * fbe_0 - g_yy_0_xxxxxyyyy_1[i] * fz_be_0 + 5.0 * g_xyy_0_xxxxyyyy_1[i] * fi_acd_0 + g_xyy_0_xxxxxyyyy_1[i] * wa_x[i];

        g_xxyy_0_xxxxxyyyz_0[i] = g_yy_0_xxxxxyyyz_0[i] * fbe_0 - g_yy_0_xxxxxyyyz_1[i] * fz_be_0 + 5.0 * g_xyy_0_xxxxyyyz_1[i] * fi_acd_0 + g_xyy_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xxyy_0_xxxxxyyzz_0[i] = g_yy_0_xxxxxyyzz_0[i] * fbe_0 - g_yy_0_xxxxxyyzz_1[i] * fz_be_0 + 5.0 * g_xyy_0_xxxxyyzz_1[i] * fi_acd_0 + g_xyy_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xxyy_0_xxxxxyzzz_0[i] = g_yy_0_xxxxxyzzz_0[i] * fbe_0 - g_yy_0_xxxxxyzzz_1[i] * fz_be_0 + 5.0 * g_xyy_0_xxxxyzzz_1[i] * fi_acd_0 + g_xyy_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xxyy_0_xxxxxzzzz_0[i] = g_xx_0_xxxxxzzzz_0[i] * fbe_0 - g_xx_0_xxxxxzzzz_1[i] * fz_be_0 + g_xxy_0_xxxxxzzzz_1[i] * wa_y[i];

        g_xxyy_0_xxxxyyyyy_0[i] = g_yy_0_xxxxyyyyy_0[i] * fbe_0 - g_yy_0_xxxxyyyyy_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyyyyy_1[i] * fi_acd_0 + g_xyy_0_xxxxyyyyy_1[i] * wa_x[i];

        g_xxyy_0_xxxxyyyyz_0[i] = g_yy_0_xxxxyyyyz_0[i] * fbe_0 - g_yy_0_xxxxyyyyz_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyyyyz_1[i] * fi_acd_0 + g_xyy_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xxyy_0_xxxxyyyzz_0[i] = g_yy_0_xxxxyyyzz_0[i] * fbe_0 - g_yy_0_xxxxyyyzz_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyyyzz_1[i] * fi_acd_0 + g_xyy_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xxyy_0_xxxxyyzzz_0[i] = g_yy_0_xxxxyyzzz_0[i] * fbe_0 - g_yy_0_xxxxyyzzz_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyyzzz_1[i] * fi_acd_0 + g_xyy_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xxyy_0_xxxxyzzzz_0[i] = g_yy_0_xxxxyzzzz_0[i] * fbe_0 - g_yy_0_xxxxyzzzz_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyzzzz_1[i] * fi_acd_0 + g_xyy_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xxyy_0_xxxxzzzzz_0[i] = g_xx_0_xxxxzzzzz_0[i] * fbe_0 - g_xx_0_xxxxzzzzz_1[i] * fz_be_0 + g_xxy_0_xxxxzzzzz_1[i] * wa_y[i];

        g_xxyy_0_xxxyyyyyy_0[i] = g_yy_0_xxxyyyyyy_0[i] * fbe_0 - g_yy_0_xxxyyyyyy_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyyyyy_1[i] * fi_acd_0 + g_xyy_0_xxxyyyyyy_1[i] * wa_x[i];

        g_xxyy_0_xxxyyyyyz_0[i] = g_yy_0_xxxyyyyyz_0[i] * fbe_0 - g_yy_0_xxxyyyyyz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyyyyz_1[i] * fi_acd_0 + g_xyy_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xxyy_0_xxxyyyyzz_0[i] = g_yy_0_xxxyyyyzz_0[i] * fbe_0 - g_yy_0_xxxyyyyzz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyyyzz_1[i] * fi_acd_0 + g_xyy_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xxyy_0_xxxyyyzzz_0[i] = g_yy_0_xxxyyyzzz_0[i] * fbe_0 - g_yy_0_xxxyyyzzz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyyzzz_1[i] * fi_acd_0 + g_xyy_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xxyy_0_xxxyyzzzz_0[i] = g_yy_0_xxxyyzzzz_0[i] * fbe_0 - g_yy_0_xxxyyzzzz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyzzzz_1[i] * fi_acd_0 + g_xyy_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xxyy_0_xxxyzzzzz_0[i] = g_yy_0_xxxyzzzzz_0[i] * fbe_0 - g_yy_0_xxxyzzzzz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyzzzzz_1[i] * fi_acd_0 + g_xyy_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xxyy_0_xxxzzzzzz_0[i] = g_xx_0_xxxzzzzzz_0[i] * fbe_0 - g_xx_0_xxxzzzzzz_1[i] * fz_be_0 + g_xxy_0_xxxzzzzzz_1[i] * wa_y[i];

        g_xxyy_0_xxyyyyyyy_0[i] = g_yy_0_xxyyyyyyy_0[i] * fbe_0 - g_yy_0_xxyyyyyyy_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyyyyy_1[i] * fi_acd_0 + g_xyy_0_xxyyyyyyy_1[i] * wa_x[i];

        g_xxyy_0_xxyyyyyyz_0[i] = g_yy_0_xxyyyyyyz_0[i] * fbe_0 - g_yy_0_xxyyyyyyz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyyyyz_1[i] * fi_acd_0 + g_xyy_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xxyy_0_xxyyyyyzz_0[i] = g_yy_0_xxyyyyyzz_0[i] * fbe_0 - g_yy_0_xxyyyyyzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyyyzz_1[i] * fi_acd_0 + g_xyy_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xxyy_0_xxyyyyzzz_0[i] = g_yy_0_xxyyyyzzz_0[i] * fbe_0 - g_yy_0_xxyyyyzzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyyzzz_1[i] * fi_acd_0 + g_xyy_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xxyy_0_xxyyyzzzz_0[i] = g_yy_0_xxyyyzzzz_0[i] * fbe_0 - g_yy_0_xxyyyzzzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyzzzz_1[i] * fi_acd_0 + g_xyy_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xxyy_0_xxyyzzzzz_0[i] = g_yy_0_xxyyzzzzz_0[i] * fbe_0 - g_yy_0_xxyyzzzzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyzzzzz_1[i] * fi_acd_0 + g_xyy_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xxyy_0_xxyzzzzzz_0[i] = g_yy_0_xxyzzzzzz_0[i] * fbe_0 - g_yy_0_xxyzzzzzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyzzzzzz_1[i] * fi_acd_0 + g_xyy_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xxyy_0_xxzzzzzzz_0[i] = g_xx_0_xxzzzzzzz_0[i] * fbe_0 - g_xx_0_xxzzzzzzz_1[i] * fz_be_0 + g_xxy_0_xxzzzzzzz_1[i] * wa_y[i];

        g_xxyy_0_xyyyyyyyy_0[i] = g_yy_0_xyyyyyyyy_0[i] * fbe_0 - g_yy_0_xyyyyyyyy_1[i] * fz_be_0 + g_xyy_0_yyyyyyyy_1[i] * fi_acd_0 + g_xyy_0_xyyyyyyyy_1[i] * wa_x[i];

        g_xxyy_0_xyyyyyyyz_0[i] = g_yy_0_xyyyyyyyz_0[i] * fbe_0 - g_yy_0_xyyyyyyyz_1[i] * fz_be_0 + g_xyy_0_yyyyyyyz_1[i] * fi_acd_0 + g_xyy_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xxyy_0_xyyyyyyzz_0[i] = g_yy_0_xyyyyyyzz_0[i] * fbe_0 - g_yy_0_xyyyyyyzz_1[i] * fz_be_0 + g_xyy_0_yyyyyyzz_1[i] * fi_acd_0 + g_xyy_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xxyy_0_xyyyyyzzz_0[i] = g_yy_0_xyyyyyzzz_0[i] * fbe_0 - g_yy_0_xyyyyyzzz_1[i] * fz_be_0 + g_xyy_0_yyyyyzzz_1[i] * fi_acd_0 + g_xyy_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xxyy_0_xyyyyzzzz_0[i] = g_yy_0_xyyyyzzzz_0[i] * fbe_0 - g_yy_0_xyyyyzzzz_1[i] * fz_be_0 + g_xyy_0_yyyyzzzz_1[i] * fi_acd_0 + g_xyy_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xxyy_0_xyyyzzzzz_0[i] = g_yy_0_xyyyzzzzz_0[i] * fbe_0 - g_yy_0_xyyyzzzzz_1[i] * fz_be_0 + g_xyy_0_yyyzzzzz_1[i] * fi_acd_0 + g_xyy_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xxyy_0_xyyzzzzzz_0[i] = g_yy_0_xyyzzzzzz_0[i] * fbe_0 - g_yy_0_xyyzzzzzz_1[i] * fz_be_0 + g_xyy_0_yyzzzzzz_1[i] * fi_acd_0 + g_xyy_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xxyy_0_xyzzzzzzz_0[i] = g_yy_0_xyzzzzzzz_0[i] * fbe_0 - g_yy_0_xyzzzzzzz_1[i] * fz_be_0 + g_xyy_0_yzzzzzzz_1[i] * fi_acd_0 + g_xyy_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xxyy_0_xzzzzzzzz_0[i] = g_xx_0_xzzzzzzzz_0[i] * fbe_0 - g_xx_0_xzzzzzzzz_1[i] * fz_be_0 + g_xxy_0_xzzzzzzzz_1[i] * wa_y[i];

        g_xxyy_0_yyyyyyyyy_0[i] = g_yy_0_yyyyyyyyy_0[i] * fbe_0 - g_yy_0_yyyyyyyyy_1[i] * fz_be_0 + g_xyy_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xxyy_0_yyyyyyyyz_0[i] = g_yy_0_yyyyyyyyz_0[i] * fbe_0 - g_yy_0_yyyyyyyyz_1[i] * fz_be_0 + g_xyy_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xxyy_0_yyyyyyyzz_0[i] = g_yy_0_yyyyyyyzz_0[i] * fbe_0 - g_yy_0_yyyyyyyzz_1[i] * fz_be_0 + g_xyy_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xxyy_0_yyyyyyzzz_0[i] = g_yy_0_yyyyyyzzz_0[i] * fbe_0 - g_yy_0_yyyyyyzzz_1[i] * fz_be_0 + g_xyy_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xxyy_0_yyyyyzzzz_0[i] = g_yy_0_yyyyyzzzz_0[i] * fbe_0 - g_yy_0_yyyyyzzzz_1[i] * fz_be_0 + g_xyy_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xxyy_0_yyyyzzzzz_0[i] = g_yy_0_yyyyzzzzz_0[i] * fbe_0 - g_yy_0_yyyyzzzzz_1[i] * fz_be_0 + g_xyy_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xxyy_0_yyyzzzzzz_0[i] = g_yy_0_yyyzzzzzz_0[i] * fbe_0 - g_yy_0_yyyzzzzzz_1[i] * fz_be_0 + g_xyy_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xxyy_0_yyzzzzzzz_0[i] = g_yy_0_yyzzzzzzz_0[i] * fbe_0 - g_yy_0_yyzzzzzzz_1[i] * fz_be_0 + g_xyy_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xxyy_0_yzzzzzzzz_0[i] = g_yy_0_yzzzzzzzz_0[i] * fbe_0 - g_yy_0_yzzzzzzzz_1[i] * fz_be_0 + g_xyy_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xxyy_0_zzzzzzzzz_0[i] = g_yy_0_zzzzzzzzz_0[i] * fbe_0 - g_yy_0_zzzzzzzzz_1[i] * fz_be_0 + g_xyy_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 220-275 components of targeted buffer : GSM

    auto g_xxyz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 220);

    auto g_xxyz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 221);

    auto g_xxyz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 222);

    auto g_xxyz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 223);

    auto g_xxyz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 224);

    auto g_xxyz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 225);

    auto g_xxyz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 226);

    auto g_xxyz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 227);

    auto g_xxyz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 228);

    auto g_xxyz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 229);

    auto g_xxyz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 230);

    auto g_xxyz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 231);

    auto g_xxyz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 232);

    auto g_xxyz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 233);

    auto g_xxyz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 234);

    auto g_xxyz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 235);

    auto g_xxyz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 236);

    auto g_xxyz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 237);

    auto g_xxyz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 238);

    auto g_xxyz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 239);

    auto g_xxyz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 240);

    auto g_xxyz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 241);

    auto g_xxyz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 242);

    auto g_xxyz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 243);

    auto g_xxyz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 244);

    auto g_xxyz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 245);

    auto g_xxyz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 246);

    auto g_xxyz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 247);

    auto g_xxyz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 248);

    auto g_xxyz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 249);

    auto g_xxyz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 250);

    auto g_xxyz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 251);

    auto g_xxyz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 252);

    auto g_xxyz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 253);

    auto g_xxyz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 254);

    auto g_xxyz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 255);

    auto g_xxyz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 256);

    auto g_xxyz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 257);

    auto g_xxyz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 258);

    auto g_xxyz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 259);

    auto g_xxyz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 260);

    auto g_xxyz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 261);

    auto g_xxyz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 262);

    auto g_xxyz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 263);

    auto g_xxyz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 264);

    auto g_xxyz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 265);

    auto g_xxyz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 266);

    auto g_xxyz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 267);

    auto g_xxyz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 268);

    auto g_xxyz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 269);

    auto g_xxyz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 270);

    auto g_xxyz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 271);

    auto g_xxyz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 272);

    auto g_xxyz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 273);

    auto g_xxyz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 274);

    #pragma omp simd aligned(g_xxy_0_xxxxxxxxy_1, g_xxy_0_xxxxxxxyy_1, g_xxy_0_xxxxxxyyy_1, g_xxy_0_xxxxxyyyy_1, g_xxy_0_xxxxyyyyy_1, g_xxy_0_xxxyyyyyy_1, g_xxy_0_xxyyyyyyy_1, g_xxy_0_xyyyyyyyy_1, g_xxy_0_yyyyyyyyy_1, g_xxyz_0_xxxxxxxxx_0, g_xxyz_0_xxxxxxxxy_0, g_xxyz_0_xxxxxxxxz_0, g_xxyz_0_xxxxxxxyy_0, g_xxyz_0_xxxxxxxyz_0, g_xxyz_0_xxxxxxxzz_0, g_xxyz_0_xxxxxxyyy_0, g_xxyz_0_xxxxxxyyz_0, g_xxyz_0_xxxxxxyzz_0, g_xxyz_0_xxxxxxzzz_0, g_xxyz_0_xxxxxyyyy_0, g_xxyz_0_xxxxxyyyz_0, g_xxyz_0_xxxxxyyzz_0, g_xxyz_0_xxxxxyzzz_0, g_xxyz_0_xxxxxzzzz_0, g_xxyz_0_xxxxyyyyy_0, g_xxyz_0_xxxxyyyyz_0, g_xxyz_0_xxxxyyyzz_0, g_xxyz_0_xxxxyyzzz_0, g_xxyz_0_xxxxyzzzz_0, g_xxyz_0_xxxxzzzzz_0, g_xxyz_0_xxxyyyyyy_0, g_xxyz_0_xxxyyyyyz_0, g_xxyz_0_xxxyyyyzz_0, g_xxyz_0_xxxyyyzzz_0, g_xxyz_0_xxxyyzzzz_0, g_xxyz_0_xxxyzzzzz_0, g_xxyz_0_xxxzzzzzz_0, g_xxyz_0_xxyyyyyyy_0, g_xxyz_0_xxyyyyyyz_0, g_xxyz_0_xxyyyyyzz_0, g_xxyz_0_xxyyyyzzz_0, g_xxyz_0_xxyyyzzzz_0, g_xxyz_0_xxyyzzzzz_0, g_xxyz_0_xxyzzzzzz_0, g_xxyz_0_xxzzzzzzz_0, g_xxyz_0_xyyyyyyyy_0, g_xxyz_0_xyyyyyyyz_0, g_xxyz_0_xyyyyyyzz_0, g_xxyz_0_xyyyyyzzz_0, g_xxyz_0_xyyyyzzzz_0, g_xxyz_0_xyyyzzzzz_0, g_xxyz_0_xyyzzzzzz_0, g_xxyz_0_xyzzzzzzz_0, g_xxyz_0_xzzzzzzzz_0, g_xxyz_0_yyyyyyyyy_0, g_xxyz_0_yyyyyyyyz_0, g_xxyz_0_yyyyyyyzz_0, g_xxyz_0_yyyyyyzzz_0, g_xxyz_0_yyyyyzzzz_0, g_xxyz_0_yyyyzzzzz_0, g_xxyz_0_yyyzzzzzz_0, g_xxyz_0_yyzzzzzzz_0, g_xxyz_0_yzzzzzzzz_0, g_xxyz_0_zzzzzzzzz_0, g_xxz_0_xxxxxxxxx_1, g_xxz_0_xxxxxxxxz_1, g_xxz_0_xxxxxxxyz_1, g_xxz_0_xxxxxxxz_1, g_xxz_0_xxxxxxxzz_1, g_xxz_0_xxxxxxyyz_1, g_xxz_0_xxxxxxyz_1, g_xxz_0_xxxxxxyzz_1, g_xxz_0_xxxxxxzz_1, g_xxz_0_xxxxxxzzz_1, g_xxz_0_xxxxxyyyz_1, g_xxz_0_xxxxxyyz_1, g_xxz_0_xxxxxyyzz_1, g_xxz_0_xxxxxyzz_1, g_xxz_0_xxxxxyzzz_1, g_xxz_0_xxxxxzzz_1, g_xxz_0_xxxxxzzzz_1, g_xxz_0_xxxxyyyyz_1, g_xxz_0_xxxxyyyz_1, g_xxz_0_xxxxyyyzz_1, g_xxz_0_xxxxyyzz_1, g_xxz_0_xxxxyyzzz_1, g_xxz_0_xxxxyzzz_1, g_xxz_0_xxxxyzzzz_1, g_xxz_0_xxxxzzzz_1, g_xxz_0_xxxxzzzzz_1, g_xxz_0_xxxyyyyyz_1, g_xxz_0_xxxyyyyz_1, g_xxz_0_xxxyyyyzz_1, g_xxz_0_xxxyyyzz_1, g_xxz_0_xxxyyyzzz_1, g_xxz_0_xxxyyzzz_1, g_xxz_0_xxxyyzzzz_1, g_xxz_0_xxxyzzzz_1, g_xxz_0_xxxyzzzzz_1, g_xxz_0_xxxzzzzz_1, g_xxz_0_xxxzzzzzz_1, g_xxz_0_xxyyyyyyz_1, g_xxz_0_xxyyyyyz_1, g_xxz_0_xxyyyyyzz_1, g_xxz_0_xxyyyyzz_1, g_xxz_0_xxyyyyzzz_1, g_xxz_0_xxyyyzzz_1, g_xxz_0_xxyyyzzzz_1, g_xxz_0_xxyyzzzz_1, g_xxz_0_xxyyzzzzz_1, g_xxz_0_xxyzzzzz_1, g_xxz_0_xxyzzzzzz_1, g_xxz_0_xxzzzzzz_1, g_xxz_0_xxzzzzzzz_1, g_xxz_0_xyyyyyyyz_1, g_xxz_0_xyyyyyyz_1, g_xxz_0_xyyyyyyzz_1, g_xxz_0_xyyyyyzz_1, g_xxz_0_xyyyyyzzz_1, g_xxz_0_xyyyyzzz_1, g_xxz_0_xyyyyzzzz_1, g_xxz_0_xyyyzzzz_1, g_xxz_0_xyyyzzzzz_1, g_xxz_0_xyyzzzzz_1, g_xxz_0_xyyzzzzzz_1, g_xxz_0_xyzzzzzz_1, g_xxz_0_xyzzzzzzz_1, g_xxz_0_xzzzzzzz_1, g_xxz_0_xzzzzzzzz_1, g_xxz_0_yyyyyyyyz_1, g_xxz_0_yyyyyyyz_1, g_xxz_0_yyyyyyyzz_1, g_xxz_0_yyyyyyzz_1, g_xxz_0_yyyyyyzzz_1, g_xxz_0_yyyyyzzz_1, g_xxz_0_yyyyyzzzz_1, g_xxz_0_yyyyzzzz_1, g_xxz_0_yyyyzzzzz_1, g_xxz_0_yyyzzzzz_1, g_xxz_0_yyyzzzzzz_1, g_xxz_0_yyzzzzzz_1, g_xxz_0_yyzzzzzzz_1, g_xxz_0_yzzzzzzz_1, g_xxz_0_yzzzzzzzz_1, g_xxz_0_zzzzzzzz_1, g_xxz_0_zzzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyz_0_xxxxxxxxx_0[i] = g_xxz_0_xxxxxxxxx_1[i] * wa_y[i];

        g_xxyz_0_xxxxxxxxy_0[i] = g_xxy_0_xxxxxxxxy_1[i] * wa_z[i];

        g_xxyz_0_xxxxxxxxz_0[i] = g_xxz_0_xxxxxxxxz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxxxyy_0[i] = g_xxy_0_xxxxxxxyy_1[i] * wa_z[i];

        g_xxyz_0_xxxxxxxyz_0[i] = g_xxz_0_xxxxxxxz_1[i] * fi_acd_0 + g_xxz_0_xxxxxxxyz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxxxzz_0[i] = g_xxz_0_xxxxxxxzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxxyyy_0[i] = g_xxy_0_xxxxxxyyy_1[i] * wa_z[i];

        g_xxyz_0_xxxxxxyyz_0[i] = 2.0 * g_xxz_0_xxxxxxyz_1[i] * fi_acd_0 + g_xxz_0_xxxxxxyyz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxxyzz_0[i] = g_xxz_0_xxxxxxzz_1[i] * fi_acd_0 + g_xxz_0_xxxxxxyzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxxzzz_0[i] = g_xxz_0_xxxxxxzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxyyyy_0[i] = g_xxy_0_xxxxxyyyy_1[i] * wa_z[i];

        g_xxyz_0_xxxxxyyyz_0[i] = 3.0 * g_xxz_0_xxxxxyyz_1[i] * fi_acd_0 + g_xxz_0_xxxxxyyyz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxyyzz_0[i] = 2.0 * g_xxz_0_xxxxxyzz_1[i] * fi_acd_0 + g_xxz_0_xxxxxyyzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxyzzz_0[i] = g_xxz_0_xxxxxzzz_1[i] * fi_acd_0 + g_xxz_0_xxxxxyzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxzzzz_0[i] = g_xxz_0_xxxxxzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxyyyyy_0[i] = g_xxy_0_xxxxyyyyy_1[i] * wa_z[i];

        g_xxyz_0_xxxxyyyyz_0[i] = 4.0 * g_xxz_0_xxxxyyyz_1[i] * fi_acd_0 + g_xxz_0_xxxxyyyyz_1[i] * wa_y[i];

        g_xxyz_0_xxxxyyyzz_0[i] = 3.0 * g_xxz_0_xxxxyyzz_1[i] * fi_acd_0 + g_xxz_0_xxxxyyyzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxyyzzz_0[i] = 2.0 * g_xxz_0_xxxxyzzz_1[i] * fi_acd_0 + g_xxz_0_xxxxyyzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxyzzzz_0[i] = g_xxz_0_xxxxzzzz_1[i] * fi_acd_0 + g_xxz_0_xxxxyzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxzzzzz_0[i] = g_xxz_0_xxxxzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxyyyyyy_0[i] = g_xxy_0_xxxyyyyyy_1[i] * wa_z[i];

        g_xxyz_0_xxxyyyyyz_0[i] = 5.0 * g_xxz_0_xxxyyyyz_1[i] * fi_acd_0 + g_xxz_0_xxxyyyyyz_1[i] * wa_y[i];

        g_xxyz_0_xxxyyyyzz_0[i] = 4.0 * g_xxz_0_xxxyyyzz_1[i] * fi_acd_0 + g_xxz_0_xxxyyyyzz_1[i] * wa_y[i];

        g_xxyz_0_xxxyyyzzz_0[i] = 3.0 * g_xxz_0_xxxyyzzz_1[i] * fi_acd_0 + g_xxz_0_xxxyyyzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxyyzzzz_0[i] = 2.0 * g_xxz_0_xxxyzzzz_1[i] * fi_acd_0 + g_xxz_0_xxxyyzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxyzzzzz_0[i] = g_xxz_0_xxxzzzzz_1[i] * fi_acd_0 + g_xxz_0_xxxyzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxzzzzzz_0[i] = g_xxz_0_xxxzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxyyyyyyy_0[i] = g_xxy_0_xxyyyyyyy_1[i] * wa_z[i];

        g_xxyz_0_xxyyyyyyz_0[i] = 6.0 * g_xxz_0_xxyyyyyz_1[i] * fi_acd_0 + g_xxz_0_xxyyyyyyz_1[i] * wa_y[i];

        g_xxyz_0_xxyyyyyzz_0[i] = 5.0 * g_xxz_0_xxyyyyzz_1[i] * fi_acd_0 + g_xxz_0_xxyyyyyzz_1[i] * wa_y[i];

        g_xxyz_0_xxyyyyzzz_0[i] = 4.0 * g_xxz_0_xxyyyzzz_1[i] * fi_acd_0 + g_xxz_0_xxyyyyzzz_1[i] * wa_y[i];

        g_xxyz_0_xxyyyzzzz_0[i] = 3.0 * g_xxz_0_xxyyzzzz_1[i] * fi_acd_0 + g_xxz_0_xxyyyzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxyyzzzzz_0[i] = 2.0 * g_xxz_0_xxyzzzzz_1[i] * fi_acd_0 + g_xxz_0_xxyyzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxyzzzzzz_0[i] = g_xxz_0_xxzzzzzz_1[i] * fi_acd_0 + g_xxz_0_xxyzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxzzzzzzz_0[i] = g_xxz_0_xxzzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xyyyyyyyy_0[i] = g_xxy_0_xyyyyyyyy_1[i] * wa_z[i];

        g_xxyz_0_xyyyyyyyz_0[i] = 7.0 * g_xxz_0_xyyyyyyz_1[i] * fi_acd_0 + g_xxz_0_xyyyyyyyz_1[i] * wa_y[i];

        g_xxyz_0_xyyyyyyzz_0[i] = 6.0 * g_xxz_0_xyyyyyzz_1[i] * fi_acd_0 + g_xxz_0_xyyyyyyzz_1[i] * wa_y[i];

        g_xxyz_0_xyyyyyzzz_0[i] = 5.0 * g_xxz_0_xyyyyzzz_1[i] * fi_acd_0 + g_xxz_0_xyyyyyzzz_1[i] * wa_y[i];

        g_xxyz_0_xyyyyzzzz_0[i] = 4.0 * g_xxz_0_xyyyzzzz_1[i] * fi_acd_0 + g_xxz_0_xyyyyzzzz_1[i] * wa_y[i];

        g_xxyz_0_xyyyzzzzz_0[i] = 3.0 * g_xxz_0_xyyzzzzz_1[i] * fi_acd_0 + g_xxz_0_xyyyzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xyyzzzzzz_0[i] = 2.0 * g_xxz_0_xyzzzzzz_1[i] * fi_acd_0 + g_xxz_0_xyyzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xyzzzzzzz_0[i] = g_xxz_0_xzzzzzzz_1[i] * fi_acd_0 + g_xxz_0_xyzzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xzzzzzzzz_0[i] = g_xxz_0_xzzzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_yyyyyyyyy_0[i] = g_xxy_0_yyyyyyyyy_1[i] * wa_z[i];

        g_xxyz_0_yyyyyyyyz_0[i] = 8.0 * g_xxz_0_yyyyyyyz_1[i] * fi_acd_0 + g_xxz_0_yyyyyyyyz_1[i] * wa_y[i];

        g_xxyz_0_yyyyyyyzz_0[i] = 7.0 * g_xxz_0_yyyyyyzz_1[i] * fi_acd_0 + g_xxz_0_yyyyyyyzz_1[i] * wa_y[i];

        g_xxyz_0_yyyyyyzzz_0[i] = 6.0 * g_xxz_0_yyyyyzzz_1[i] * fi_acd_0 + g_xxz_0_yyyyyyzzz_1[i] * wa_y[i];

        g_xxyz_0_yyyyyzzzz_0[i] = 5.0 * g_xxz_0_yyyyzzzz_1[i] * fi_acd_0 + g_xxz_0_yyyyyzzzz_1[i] * wa_y[i];

        g_xxyz_0_yyyyzzzzz_0[i] = 4.0 * g_xxz_0_yyyzzzzz_1[i] * fi_acd_0 + g_xxz_0_yyyyzzzzz_1[i] * wa_y[i];

        g_xxyz_0_yyyzzzzzz_0[i] = 3.0 * g_xxz_0_yyzzzzzz_1[i] * fi_acd_0 + g_xxz_0_yyyzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_yyzzzzzzz_0[i] = 2.0 * g_xxz_0_yzzzzzzz_1[i] * fi_acd_0 + g_xxz_0_yyzzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_yzzzzzzzz_0[i] = g_xxz_0_zzzzzzzz_1[i] * fi_acd_0 + g_xxz_0_yzzzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_zzzzzzzzz_0[i] = g_xxz_0_zzzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 275-330 components of targeted buffer : GSM

    auto g_xxzz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 275);

    auto g_xxzz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 276);

    auto g_xxzz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 277);

    auto g_xxzz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 278);

    auto g_xxzz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 279);

    auto g_xxzz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 280);

    auto g_xxzz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 281);

    auto g_xxzz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 282);

    auto g_xxzz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 283);

    auto g_xxzz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 284);

    auto g_xxzz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 285);

    auto g_xxzz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 286);

    auto g_xxzz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 287);

    auto g_xxzz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 288);

    auto g_xxzz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 289);

    auto g_xxzz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 290);

    auto g_xxzz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 291);

    auto g_xxzz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 292);

    auto g_xxzz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 293);

    auto g_xxzz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 294);

    auto g_xxzz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 295);

    auto g_xxzz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 296);

    auto g_xxzz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 297);

    auto g_xxzz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 298);

    auto g_xxzz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 299);

    auto g_xxzz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 300);

    auto g_xxzz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 301);

    auto g_xxzz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 302);

    auto g_xxzz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 303);

    auto g_xxzz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 304);

    auto g_xxzz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 305);

    auto g_xxzz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 306);

    auto g_xxzz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 307);

    auto g_xxzz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 308);

    auto g_xxzz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 309);

    auto g_xxzz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 310);

    auto g_xxzz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 311);

    auto g_xxzz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 312);

    auto g_xxzz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 313);

    auto g_xxzz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 314);

    auto g_xxzz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 315);

    auto g_xxzz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 316);

    auto g_xxzz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 317);

    auto g_xxzz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 318);

    auto g_xxzz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 319);

    auto g_xxzz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 320);

    auto g_xxzz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 321);

    auto g_xxzz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 322);

    auto g_xxzz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 323);

    auto g_xxzz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 324);

    auto g_xxzz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 325);

    auto g_xxzz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 326);

    auto g_xxzz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 327);

    auto g_xxzz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 328);

    auto g_xxzz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 329);

    #pragma omp simd aligned(g_xx_0_xxxxxxxxx_0, g_xx_0_xxxxxxxxx_1, g_xx_0_xxxxxxxxy_0, g_xx_0_xxxxxxxxy_1, g_xx_0_xxxxxxxyy_0, g_xx_0_xxxxxxxyy_1, g_xx_0_xxxxxxyyy_0, g_xx_0_xxxxxxyyy_1, g_xx_0_xxxxxyyyy_0, g_xx_0_xxxxxyyyy_1, g_xx_0_xxxxyyyyy_0, g_xx_0_xxxxyyyyy_1, g_xx_0_xxxyyyyyy_0, g_xx_0_xxxyyyyyy_1, g_xx_0_xxyyyyyyy_0, g_xx_0_xxyyyyyyy_1, g_xx_0_xyyyyyyyy_0, g_xx_0_xyyyyyyyy_1, g_xxz_0_xxxxxxxxx_1, g_xxz_0_xxxxxxxxy_1, g_xxz_0_xxxxxxxyy_1, g_xxz_0_xxxxxxyyy_1, g_xxz_0_xxxxxyyyy_1, g_xxz_0_xxxxyyyyy_1, g_xxz_0_xxxyyyyyy_1, g_xxz_0_xxyyyyyyy_1, g_xxz_0_xyyyyyyyy_1, g_xxzz_0_xxxxxxxxx_0, g_xxzz_0_xxxxxxxxy_0, g_xxzz_0_xxxxxxxxz_0, g_xxzz_0_xxxxxxxyy_0, g_xxzz_0_xxxxxxxyz_0, g_xxzz_0_xxxxxxxzz_0, g_xxzz_0_xxxxxxyyy_0, g_xxzz_0_xxxxxxyyz_0, g_xxzz_0_xxxxxxyzz_0, g_xxzz_0_xxxxxxzzz_0, g_xxzz_0_xxxxxyyyy_0, g_xxzz_0_xxxxxyyyz_0, g_xxzz_0_xxxxxyyzz_0, g_xxzz_0_xxxxxyzzz_0, g_xxzz_0_xxxxxzzzz_0, g_xxzz_0_xxxxyyyyy_0, g_xxzz_0_xxxxyyyyz_0, g_xxzz_0_xxxxyyyzz_0, g_xxzz_0_xxxxyyzzz_0, g_xxzz_0_xxxxyzzzz_0, g_xxzz_0_xxxxzzzzz_0, g_xxzz_0_xxxyyyyyy_0, g_xxzz_0_xxxyyyyyz_0, g_xxzz_0_xxxyyyyzz_0, g_xxzz_0_xxxyyyzzz_0, g_xxzz_0_xxxyyzzzz_0, g_xxzz_0_xxxyzzzzz_0, g_xxzz_0_xxxzzzzzz_0, g_xxzz_0_xxyyyyyyy_0, g_xxzz_0_xxyyyyyyz_0, g_xxzz_0_xxyyyyyzz_0, g_xxzz_0_xxyyyyzzz_0, g_xxzz_0_xxyyyzzzz_0, g_xxzz_0_xxyyzzzzz_0, g_xxzz_0_xxyzzzzzz_0, g_xxzz_0_xxzzzzzzz_0, g_xxzz_0_xyyyyyyyy_0, g_xxzz_0_xyyyyyyyz_0, g_xxzz_0_xyyyyyyzz_0, g_xxzz_0_xyyyyyzzz_0, g_xxzz_0_xyyyyzzzz_0, g_xxzz_0_xyyyzzzzz_0, g_xxzz_0_xyyzzzzzz_0, g_xxzz_0_xyzzzzzzz_0, g_xxzz_0_xzzzzzzzz_0, g_xxzz_0_yyyyyyyyy_0, g_xxzz_0_yyyyyyyyz_0, g_xxzz_0_yyyyyyyzz_0, g_xxzz_0_yyyyyyzzz_0, g_xxzz_0_yyyyyzzzz_0, g_xxzz_0_yyyyzzzzz_0, g_xxzz_0_yyyzzzzzz_0, g_xxzz_0_yyzzzzzzz_0, g_xxzz_0_yzzzzzzzz_0, g_xxzz_0_zzzzzzzzz_0, g_xzz_0_xxxxxxxxz_1, g_xzz_0_xxxxxxxyz_1, g_xzz_0_xxxxxxxz_1, g_xzz_0_xxxxxxxzz_1, g_xzz_0_xxxxxxyyz_1, g_xzz_0_xxxxxxyz_1, g_xzz_0_xxxxxxyzz_1, g_xzz_0_xxxxxxzz_1, g_xzz_0_xxxxxxzzz_1, g_xzz_0_xxxxxyyyz_1, g_xzz_0_xxxxxyyz_1, g_xzz_0_xxxxxyyzz_1, g_xzz_0_xxxxxyzz_1, g_xzz_0_xxxxxyzzz_1, g_xzz_0_xxxxxzzz_1, g_xzz_0_xxxxxzzzz_1, g_xzz_0_xxxxyyyyz_1, g_xzz_0_xxxxyyyz_1, g_xzz_0_xxxxyyyzz_1, g_xzz_0_xxxxyyzz_1, g_xzz_0_xxxxyyzzz_1, g_xzz_0_xxxxyzzz_1, g_xzz_0_xxxxyzzzz_1, g_xzz_0_xxxxzzzz_1, g_xzz_0_xxxxzzzzz_1, g_xzz_0_xxxyyyyyz_1, g_xzz_0_xxxyyyyz_1, g_xzz_0_xxxyyyyzz_1, g_xzz_0_xxxyyyzz_1, g_xzz_0_xxxyyyzzz_1, g_xzz_0_xxxyyzzz_1, g_xzz_0_xxxyyzzzz_1, g_xzz_0_xxxyzzzz_1, g_xzz_0_xxxyzzzzz_1, g_xzz_0_xxxzzzzz_1, g_xzz_0_xxxzzzzzz_1, g_xzz_0_xxyyyyyyz_1, g_xzz_0_xxyyyyyz_1, g_xzz_0_xxyyyyyzz_1, g_xzz_0_xxyyyyzz_1, g_xzz_0_xxyyyyzzz_1, g_xzz_0_xxyyyzzz_1, g_xzz_0_xxyyyzzzz_1, g_xzz_0_xxyyzzzz_1, g_xzz_0_xxyyzzzzz_1, g_xzz_0_xxyzzzzz_1, g_xzz_0_xxyzzzzzz_1, g_xzz_0_xxzzzzzz_1, g_xzz_0_xxzzzzzzz_1, g_xzz_0_xyyyyyyyz_1, g_xzz_0_xyyyyyyz_1, g_xzz_0_xyyyyyyzz_1, g_xzz_0_xyyyyyzz_1, g_xzz_0_xyyyyyzzz_1, g_xzz_0_xyyyyzzz_1, g_xzz_0_xyyyyzzzz_1, g_xzz_0_xyyyzzzz_1, g_xzz_0_xyyyzzzzz_1, g_xzz_0_xyyzzzzz_1, g_xzz_0_xyyzzzzzz_1, g_xzz_0_xyzzzzzz_1, g_xzz_0_xyzzzzzzz_1, g_xzz_0_xzzzzzzz_1, g_xzz_0_xzzzzzzzz_1, g_xzz_0_yyyyyyyyy_1, g_xzz_0_yyyyyyyyz_1, g_xzz_0_yyyyyyyz_1, g_xzz_0_yyyyyyyzz_1, g_xzz_0_yyyyyyzz_1, g_xzz_0_yyyyyyzzz_1, g_xzz_0_yyyyyzzz_1, g_xzz_0_yyyyyzzzz_1, g_xzz_0_yyyyzzzz_1, g_xzz_0_yyyyzzzzz_1, g_xzz_0_yyyzzzzz_1, g_xzz_0_yyyzzzzzz_1, g_xzz_0_yyzzzzzz_1, g_xzz_0_yyzzzzzzz_1, g_xzz_0_yzzzzzzz_1, g_xzz_0_yzzzzzzzz_1, g_xzz_0_zzzzzzzz_1, g_xzz_0_zzzzzzzzz_1, g_zz_0_xxxxxxxxz_0, g_zz_0_xxxxxxxxz_1, g_zz_0_xxxxxxxyz_0, g_zz_0_xxxxxxxyz_1, g_zz_0_xxxxxxxzz_0, g_zz_0_xxxxxxxzz_1, g_zz_0_xxxxxxyyz_0, g_zz_0_xxxxxxyyz_1, g_zz_0_xxxxxxyzz_0, g_zz_0_xxxxxxyzz_1, g_zz_0_xxxxxxzzz_0, g_zz_0_xxxxxxzzz_1, g_zz_0_xxxxxyyyz_0, g_zz_0_xxxxxyyyz_1, g_zz_0_xxxxxyyzz_0, g_zz_0_xxxxxyyzz_1, g_zz_0_xxxxxyzzz_0, g_zz_0_xxxxxyzzz_1, g_zz_0_xxxxxzzzz_0, g_zz_0_xxxxxzzzz_1, g_zz_0_xxxxyyyyz_0, g_zz_0_xxxxyyyyz_1, g_zz_0_xxxxyyyzz_0, g_zz_0_xxxxyyyzz_1, g_zz_0_xxxxyyzzz_0, g_zz_0_xxxxyyzzz_1, g_zz_0_xxxxyzzzz_0, g_zz_0_xxxxyzzzz_1, g_zz_0_xxxxzzzzz_0, g_zz_0_xxxxzzzzz_1, g_zz_0_xxxyyyyyz_0, g_zz_0_xxxyyyyyz_1, g_zz_0_xxxyyyyzz_0, g_zz_0_xxxyyyyzz_1, g_zz_0_xxxyyyzzz_0, g_zz_0_xxxyyyzzz_1, g_zz_0_xxxyyzzzz_0, g_zz_0_xxxyyzzzz_1, g_zz_0_xxxyzzzzz_0, g_zz_0_xxxyzzzzz_1, g_zz_0_xxxzzzzzz_0, g_zz_0_xxxzzzzzz_1, g_zz_0_xxyyyyyyz_0, g_zz_0_xxyyyyyyz_1, g_zz_0_xxyyyyyzz_0, g_zz_0_xxyyyyyzz_1, g_zz_0_xxyyyyzzz_0, g_zz_0_xxyyyyzzz_1, g_zz_0_xxyyyzzzz_0, g_zz_0_xxyyyzzzz_1, g_zz_0_xxyyzzzzz_0, g_zz_0_xxyyzzzzz_1, g_zz_0_xxyzzzzzz_0, g_zz_0_xxyzzzzzz_1, g_zz_0_xxzzzzzzz_0, g_zz_0_xxzzzzzzz_1, g_zz_0_xyyyyyyyz_0, g_zz_0_xyyyyyyyz_1, g_zz_0_xyyyyyyzz_0, g_zz_0_xyyyyyyzz_1, g_zz_0_xyyyyyzzz_0, g_zz_0_xyyyyyzzz_1, g_zz_0_xyyyyzzzz_0, g_zz_0_xyyyyzzzz_1, g_zz_0_xyyyzzzzz_0, g_zz_0_xyyyzzzzz_1, g_zz_0_xyyzzzzzz_0, g_zz_0_xyyzzzzzz_1, g_zz_0_xyzzzzzzz_0, g_zz_0_xyzzzzzzz_1, g_zz_0_xzzzzzzzz_0, g_zz_0_xzzzzzzzz_1, g_zz_0_yyyyyyyyy_0, g_zz_0_yyyyyyyyy_1, g_zz_0_yyyyyyyyz_0, g_zz_0_yyyyyyyyz_1, g_zz_0_yyyyyyyzz_0, g_zz_0_yyyyyyyzz_1, g_zz_0_yyyyyyzzz_0, g_zz_0_yyyyyyzzz_1, g_zz_0_yyyyyzzzz_0, g_zz_0_yyyyyzzzz_1, g_zz_0_yyyyzzzzz_0, g_zz_0_yyyyzzzzz_1, g_zz_0_yyyzzzzzz_0, g_zz_0_yyyzzzzzz_1, g_zz_0_yyzzzzzzz_0, g_zz_0_yyzzzzzzz_1, g_zz_0_yzzzzzzzz_0, g_zz_0_yzzzzzzzz_1, g_zz_0_zzzzzzzzz_0, g_zz_0_zzzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzz_0_xxxxxxxxx_0[i] = g_xx_0_xxxxxxxxx_0[i] * fbe_0 - g_xx_0_xxxxxxxxx_1[i] * fz_be_0 + g_xxz_0_xxxxxxxxx_1[i] * wa_z[i];

        g_xxzz_0_xxxxxxxxy_0[i] = g_xx_0_xxxxxxxxy_0[i] * fbe_0 - g_xx_0_xxxxxxxxy_1[i] * fz_be_0 + g_xxz_0_xxxxxxxxy_1[i] * wa_z[i];

        g_xxzz_0_xxxxxxxxz_0[i] = g_zz_0_xxxxxxxxz_0[i] * fbe_0 - g_zz_0_xxxxxxxxz_1[i] * fz_be_0 + 8.0 * g_xzz_0_xxxxxxxz_1[i] * fi_acd_0 + g_xzz_0_xxxxxxxxz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxxxyy_0[i] = g_xx_0_xxxxxxxyy_0[i] * fbe_0 - g_xx_0_xxxxxxxyy_1[i] * fz_be_0 + g_xxz_0_xxxxxxxyy_1[i] * wa_z[i];

        g_xxzz_0_xxxxxxxyz_0[i] = g_zz_0_xxxxxxxyz_0[i] * fbe_0 - g_zz_0_xxxxxxxyz_1[i] * fz_be_0 + 7.0 * g_xzz_0_xxxxxxyz_1[i] * fi_acd_0 + g_xzz_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxxxzz_0[i] = g_zz_0_xxxxxxxzz_0[i] * fbe_0 - g_zz_0_xxxxxxxzz_1[i] * fz_be_0 + 7.0 * g_xzz_0_xxxxxxzz_1[i] * fi_acd_0 + g_xzz_0_xxxxxxxzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxxyyy_0[i] = g_xx_0_xxxxxxyyy_0[i] * fbe_0 - g_xx_0_xxxxxxyyy_1[i] * fz_be_0 + g_xxz_0_xxxxxxyyy_1[i] * wa_z[i];

        g_xxzz_0_xxxxxxyyz_0[i] = g_zz_0_xxxxxxyyz_0[i] * fbe_0 - g_zz_0_xxxxxxyyz_1[i] * fz_be_0 + 6.0 * g_xzz_0_xxxxxyyz_1[i] * fi_acd_0 + g_xzz_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxxyzz_0[i] = g_zz_0_xxxxxxyzz_0[i] * fbe_0 - g_zz_0_xxxxxxyzz_1[i] * fz_be_0 + 6.0 * g_xzz_0_xxxxxyzz_1[i] * fi_acd_0 + g_xzz_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxxzzz_0[i] = g_zz_0_xxxxxxzzz_0[i] * fbe_0 - g_zz_0_xxxxxxzzz_1[i] * fz_be_0 + 6.0 * g_xzz_0_xxxxxzzz_1[i] * fi_acd_0 + g_xzz_0_xxxxxxzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxyyyy_0[i] = g_xx_0_xxxxxyyyy_0[i] * fbe_0 - g_xx_0_xxxxxyyyy_1[i] * fz_be_0 + g_xxz_0_xxxxxyyyy_1[i] * wa_z[i];

        g_xxzz_0_xxxxxyyyz_0[i] = g_zz_0_xxxxxyyyz_0[i] * fbe_0 - g_zz_0_xxxxxyyyz_1[i] * fz_be_0 + 5.0 * g_xzz_0_xxxxyyyz_1[i] * fi_acd_0 + g_xzz_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxyyzz_0[i] = g_zz_0_xxxxxyyzz_0[i] * fbe_0 - g_zz_0_xxxxxyyzz_1[i] * fz_be_0 + 5.0 * g_xzz_0_xxxxyyzz_1[i] * fi_acd_0 + g_xzz_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxyzzz_0[i] = g_zz_0_xxxxxyzzz_0[i] * fbe_0 - g_zz_0_xxxxxyzzz_1[i] * fz_be_0 + 5.0 * g_xzz_0_xxxxyzzz_1[i] * fi_acd_0 + g_xzz_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxzzzz_0[i] = g_zz_0_xxxxxzzzz_0[i] * fbe_0 - g_zz_0_xxxxxzzzz_1[i] * fz_be_0 + 5.0 * g_xzz_0_xxxxzzzz_1[i] * fi_acd_0 + g_xzz_0_xxxxxzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxyyyyy_0[i] = g_xx_0_xxxxyyyyy_0[i] * fbe_0 - g_xx_0_xxxxyyyyy_1[i] * fz_be_0 + g_xxz_0_xxxxyyyyy_1[i] * wa_z[i];

        g_xxzz_0_xxxxyyyyz_0[i] = g_zz_0_xxxxyyyyz_0[i] * fbe_0 - g_zz_0_xxxxyyyyz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxyyyyz_1[i] * fi_acd_0 + g_xzz_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xxzz_0_xxxxyyyzz_0[i] = g_zz_0_xxxxyyyzz_0[i] * fbe_0 - g_zz_0_xxxxyyyzz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxyyyzz_1[i] * fi_acd_0 + g_xzz_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxyyzzz_0[i] = g_zz_0_xxxxyyzzz_0[i] * fbe_0 - g_zz_0_xxxxyyzzz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxyyzzz_1[i] * fi_acd_0 + g_xzz_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxyzzzz_0[i] = g_zz_0_xxxxyzzzz_0[i] * fbe_0 - g_zz_0_xxxxyzzzz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxyzzzz_1[i] * fi_acd_0 + g_xzz_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxzzzzz_0[i] = g_zz_0_xxxxzzzzz_0[i] * fbe_0 - g_zz_0_xxxxzzzzz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxzzzzz_1[i] * fi_acd_0 + g_xzz_0_xxxxzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxyyyyyy_0[i] = g_xx_0_xxxyyyyyy_0[i] * fbe_0 - g_xx_0_xxxyyyyyy_1[i] * fz_be_0 + g_xxz_0_xxxyyyyyy_1[i] * wa_z[i];

        g_xxzz_0_xxxyyyyyz_0[i] = g_zz_0_xxxyyyyyz_0[i] * fbe_0 - g_zz_0_xxxyyyyyz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyyyyyz_1[i] * fi_acd_0 + g_xzz_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xxzz_0_xxxyyyyzz_0[i] = g_zz_0_xxxyyyyzz_0[i] * fbe_0 - g_zz_0_xxxyyyyzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyyyyzz_1[i] * fi_acd_0 + g_xzz_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xxzz_0_xxxyyyzzz_0[i] = g_zz_0_xxxyyyzzz_0[i] * fbe_0 - g_zz_0_xxxyyyzzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyyyzzz_1[i] * fi_acd_0 + g_xzz_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxyyzzzz_0[i] = g_zz_0_xxxyyzzzz_0[i] * fbe_0 - g_zz_0_xxxyyzzzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyyzzzz_1[i] * fi_acd_0 + g_xzz_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxyzzzzz_0[i] = g_zz_0_xxxyzzzzz_0[i] * fbe_0 - g_zz_0_xxxyzzzzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyzzzzz_1[i] * fi_acd_0 + g_xzz_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxzzzzzz_0[i] = g_zz_0_xxxzzzzzz_0[i] * fbe_0 - g_zz_0_xxxzzzzzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxzzzzzz_1[i] * fi_acd_0 + g_xzz_0_xxxzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxyyyyyyy_0[i] = g_xx_0_xxyyyyyyy_0[i] * fbe_0 - g_xx_0_xxyyyyyyy_1[i] * fz_be_0 + g_xxz_0_xxyyyyyyy_1[i] * wa_z[i];

        g_xxzz_0_xxyyyyyyz_0[i] = g_zz_0_xxyyyyyyz_0[i] * fbe_0 - g_zz_0_xxyyyyyyz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyyyyyz_1[i] * fi_acd_0 + g_xzz_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xxzz_0_xxyyyyyzz_0[i] = g_zz_0_xxyyyyyzz_0[i] * fbe_0 - g_zz_0_xxyyyyyzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyyyyzz_1[i] * fi_acd_0 + g_xzz_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xxzz_0_xxyyyyzzz_0[i] = g_zz_0_xxyyyyzzz_0[i] * fbe_0 - g_zz_0_xxyyyyzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyyyzzz_1[i] * fi_acd_0 + g_xzz_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xxzz_0_xxyyyzzzz_0[i] = g_zz_0_xxyyyzzzz_0[i] * fbe_0 - g_zz_0_xxyyyzzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyyzzzz_1[i] * fi_acd_0 + g_xzz_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxyyzzzzz_0[i] = g_zz_0_xxyyzzzzz_0[i] * fbe_0 - g_zz_0_xxyyzzzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyzzzzz_1[i] * fi_acd_0 + g_xzz_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxyzzzzzz_0[i] = g_zz_0_xxyzzzzzz_0[i] * fbe_0 - g_zz_0_xxyzzzzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyzzzzzz_1[i] * fi_acd_0 + g_xzz_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxzzzzzzz_0[i] = g_zz_0_xxzzzzzzz_0[i] * fbe_0 - g_zz_0_xxzzzzzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xzzzzzzz_1[i] * fi_acd_0 + g_xzz_0_xxzzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xyyyyyyyy_0[i] = g_xx_0_xyyyyyyyy_0[i] * fbe_0 - g_xx_0_xyyyyyyyy_1[i] * fz_be_0 + g_xxz_0_xyyyyyyyy_1[i] * wa_z[i];

        g_xxzz_0_xyyyyyyyz_0[i] = g_zz_0_xyyyyyyyz_0[i] * fbe_0 - g_zz_0_xyyyyyyyz_1[i] * fz_be_0 + g_xzz_0_yyyyyyyz_1[i] * fi_acd_0 + g_xzz_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xxzz_0_xyyyyyyzz_0[i] = g_zz_0_xyyyyyyzz_0[i] * fbe_0 - g_zz_0_xyyyyyyzz_1[i] * fz_be_0 + g_xzz_0_yyyyyyzz_1[i] * fi_acd_0 + g_xzz_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xxzz_0_xyyyyyzzz_0[i] = g_zz_0_xyyyyyzzz_0[i] * fbe_0 - g_zz_0_xyyyyyzzz_1[i] * fz_be_0 + g_xzz_0_yyyyyzzz_1[i] * fi_acd_0 + g_xzz_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xxzz_0_xyyyyzzzz_0[i] = g_zz_0_xyyyyzzzz_0[i] * fbe_0 - g_zz_0_xyyyyzzzz_1[i] * fz_be_0 + g_xzz_0_yyyyzzzz_1[i] * fi_acd_0 + g_xzz_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xxzz_0_xyyyzzzzz_0[i] = g_zz_0_xyyyzzzzz_0[i] * fbe_0 - g_zz_0_xyyyzzzzz_1[i] * fz_be_0 + g_xzz_0_yyyzzzzz_1[i] * fi_acd_0 + g_xzz_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xyyzzzzzz_0[i] = g_zz_0_xyyzzzzzz_0[i] * fbe_0 - g_zz_0_xyyzzzzzz_1[i] * fz_be_0 + g_xzz_0_yyzzzzzz_1[i] * fi_acd_0 + g_xzz_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xyzzzzzzz_0[i] = g_zz_0_xyzzzzzzz_0[i] * fbe_0 - g_zz_0_xyzzzzzzz_1[i] * fz_be_0 + g_xzz_0_yzzzzzzz_1[i] * fi_acd_0 + g_xzz_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xzzzzzzzz_0[i] = g_zz_0_xzzzzzzzz_0[i] * fbe_0 - g_zz_0_xzzzzzzzz_1[i] * fz_be_0 + g_xzz_0_zzzzzzzz_1[i] * fi_acd_0 + g_xzz_0_xzzzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_yyyyyyyyy_0[i] = g_zz_0_yyyyyyyyy_0[i] * fbe_0 - g_zz_0_yyyyyyyyy_1[i] * fz_be_0 + g_xzz_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xxzz_0_yyyyyyyyz_0[i] = g_zz_0_yyyyyyyyz_0[i] * fbe_0 - g_zz_0_yyyyyyyyz_1[i] * fz_be_0 + g_xzz_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xxzz_0_yyyyyyyzz_0[i] = g_zz_0_yyyyyyyzz_0[i] * fbe_0 - g_zz_0_yyyyyyyzz_1[i] * fz_be_0 + g_xzz_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xxzz_0_yyyyyyzzz_0[i] = g_zz_0_yyyyyyzzz_0[i] * fbe_0 - g_zz_0_yyyyyyzzz_1[i] * fz_be_0 + g_xzz_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xxzz_0_yyyyyzzzz_0[i] = g_zz_0_yyyyyzzzz_0[i] * fbe_0 - g_zz_0_yyyyyzzzz_1[i] * fz_be_0 + g_xzz_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xxzz_0_yyyyzzzzz_0[i] = g_zz_0_yyyyzzzzz_0[i] * fbe_0 - g_zz_0_yyyyzzzzz_1[i] * fz_be_0 + g_xzz_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xxzz_0_yyyzzzzzz_0[i] = g_zz_0_yyyzzzzzz_0[i] * fbe_0 - g_zz_0_yyyzzzzzz_1[i] * fz_be_0 + g_xzz_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_yyzzzzzzz_0[i] = g_zz_0_yyzzzzzzz_0[i] * fbe_0 - g_zz_0_yyzzzzzzz_1[i] * fz_be_0 + g_xzz_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_yzzzzzzzz_0[i] = g_zz_0_yzzzzzzzz_0[i] * fbe_0 - g_zz_0_yzzzzzzzz_1[i] * fz_be_0 + g_xzz_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_zzzzzzzzz_0[i] = g_zz_0_zzzzzzzzz_0[i] * fbe_0 - g_zz_0_zzzzzzzzz_1[i] * fz_be_0 + g_xzz_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 330-385 components of targeted buffer : GSM

    auto g_xyyy_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 330);

    auto g_xyyy_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 331);

    auto g_xyyy_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 332);

    auto g_xyyy_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 333);

    auto g_xyyy_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 334);

    auto g_xyyy_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 335);

    auto g_xyyy_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 336);

    auto g_xyyy_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 337);

    auto g_xyyy_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 338);

    auto g_xyyy_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 339);

    auto g_xyyy_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 340);

    auto g_xyyy_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 341);

    auto g_xyyy_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 342);

    auto g_xyyy_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 343);

    auto g_xyyy_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 344);

    auto g_xyyy_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 345);

    auto g_xyyy_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 346);

    auto g_xyyy_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 347);

    auto g_xyyy_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 348);

    auto g_xyyy_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 349);

    auto g_xyyy_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 350);

    auto g_xyyy_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 351);

    auto g_xyyy_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 352);

    auto g_xyyy_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 353);

    auto g_xyyy_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 354);

    auto g_xyyy_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 355);

    auto g_xyyy_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 356);

    auto g_xyyy_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 357);

    auto g_xyyy_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 358);

    auto g_xyyy_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 359);

    auto g_xyyy_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 360);

    auto g_xyyy_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 361);

    auto g_xyyy_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 362);

    auto g_xyyy_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 363);

    auto g_xyyy_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 364);

    auto g_xyyy_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 365);

    auto g_xyyy_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 366);

    auto g_xyyy_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 367);

    auto g_xyyy_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 368);

    auto g_xyyy_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 369);

    auto g_xyyy_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 370);

    auto g_xyyy_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 371);

    auto g_xyyy_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 372);

    auto g_xyyy_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 373);

    auto g_xyyy_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 374);

    auto g_xyyy_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 375);

    auto g_xyyy_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 376);

    auto g_xyyy_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 377);

    auto g_xyyy_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 378);

    auto g_xyyy_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 379);

    auto g_xyyy_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 380);

    auto g_xyyy_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 381);

    auto g_xyyy_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 382);

    auto g_xyyy_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 383);

    auto g_xyyy_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 384);

    #pragma omp simd aligned(g_xyyy_0_xxxxxxxxx_0, g_xyyy_0_xxxxxxxxy_0, g_xyyy_0_xxxxxxxxz_0, g_xyyy_0_xxxxxxxyy_0, g_xyyy_0_xxxxxxxyz_0, g_xyyy_0_xxxxxxxzz_0, g_xyyy_0_xxxxxxyyy_0, g_xyyy_0_xxxxxxyyz_0, g_xyyy_0_xxxxxxyzz_0, g_xyyy_0_xxxxxxzzz_0, g_xyyy_0_xxxxxyyyy_0, g_xyyy_0_xxxxxyyyz_0, g_xyyy_0_xxxxxyyzz_0, g_xyyy_0_xxxxxyzzz_0, g_xyyy_0_xxxxxzzzz_0, g_xyyy_0_xxxxyyyyy_0, g_xyyy_0_xxxxyyyyz_0, g_xyyy_0_xxxxyyyzz_0, g_xyyy_0_xxxxyyzzz_0, g_xyyy_0_xxxxyzzzz_0, g_xyyy_0_xxxxzzzzz_0, g_xyyy_0_xxxyyyyyy_0, g_xyyy_0_xxxyyyyyz_0, g_xyyy_0_xxxyyyyzz_0, g_xyyy_0_xxxyyyzzz_0, g_xyyy_0_xxxyyzzzz_0, g_xyyy_0_xxxyzzzzz_0, g_xyyy_0_xxxzzzzzz_0, g_xyyy_0_xxyyyyyyy_0, g_xyyy_0_xxyyyyyyz_0, g_xyyy_0_xxyyyyyzz_0, g_xyyy_0_xxyyyyzzz_0, g_xyyy_0_xxyyyzzzz_0, g_xyyy_0_xxyyzzzzz_0, g_xyyy_0_xxyzzzzzz_0, g_xyyy_0_xxzzzzzzz_0, g_xyyy_0_xyyyyyyyy_0, g_xyyy_0_xyyyyyyyz_0, g_xyyy_0_xyyyyyyzz_0, g_xyyy_0_xyyyyyzzz_0, g_xyyy_0_xyyyyzzzz_0, g_xyyy_0_xyyyzzzzz_0, g_xyyy_0_xyyzzzzzz_0, g_xyyy_0_xyzzzzzzz_0, g_xyyy_0_xzzzzzzzz_0, g_xyyy_0_yyyyyyyyy_0, g_xyyy_0_yyyyyyyyz_0, g_xyyy_0_yyyyyyyzz_0, g_xyyy_0_yyyyyyzzz_0, g_xyyy_0_yyyyyzzzz_0, g_xyyy_0_yyyyzzzzz_0, g_xyyy_0_yyyzzzzzz_0, g_xyyy_0_yyzzzzzzz_0, g_xyyy_0_yzzzzzzzz_0, g_xyyy_0_zzzzzzzzz_0, g_yyy_0_xxxxxxxx_1, g_yyy_0_xxxxxxxxx_1, g_yyy_0_xxxxxxxxy_1, g_yyy_0_xxxxxxxxz_1, g_yyy_0_xxxxxxxy_1, g_yyy_0_xxxxxxxyy_1, g_yyy_0_xxxxxxxyz_1, g_yyy_0_xxxxxxxz_1, g_yyy_0_xxxxxxxzz_1, g_yyy_0_xxxxxxyy_1, g_yyy_0_xxxxxxyyy_1, g_yyy_0_xxxxxxyyz_1, g_yyy_0_xxxxxxyz_1, g_yyy_0_xxxxxxyzz_1, g_yyy_0_xxxxxxzz_1, g_yyy_0_xxxxxxzzz_1, g_yyy_0_xxxxxyyy_1, g_yyy_0_xxxxxyyyy_1, g_yyy_0_xxxxxyyyz_1, g_yyy_0_xxxxxyyz_1, g_yyy_0_xxxxxyyzz_1, g_yyy_0_xxxxxyzz_1, g_yyy_0_xxxxxyzzz_1, g_yyy_0_xxxxxzzz_1, g_yyy_0_xxxxxzzzz_1, g_yyy_0_xxxxyyyy_1, g_yyy_0_xxxxyyyyy_1, g_yyy_0_xxxxyyyyz_1, g_yyy_0_xxxxyyyz_1, g_yyy_0_xxxxyyyzz_1, g_yyy_0_xxxxyyzz_1, g_yyy_0_xxxxyyzzz_1, g_yyy_0_xxxxyzzz_1, g_yyy_0_xxxxyzzzz_1, g_yyy_0_xxxxzzzz_1, g_yyy_0_xxxxzzzzz_1, g_yyy_0_xxxyyyyy_1, g_yyy_0_xxxyyyyyy_1, g_yyy_0_xxxyyyyyz_1, g_yyy_0_xxxyyyyz_1, g_yyy_0_xxxyyyyzz_1, g_yyy_0_xxxyyyzz_1, g_yyy_0_xxxyyyzzz_1, g_yyy_0_xxxyyzzz_1, g_yyy_0_xxxyyzzzz_1, g_yyy_0_xxxyzzzz_1, g_yyy_0_xxxyzzzzz_1, g_yyy_0_xxxzzzzz_1, g_yyy_0_xxxzzzzzz_1, g_yyy_0_xxyyyyyy_1, g_yyy_0_xxyyyyyyy_1, g_yyy_0_xxyyyyyyz_1, g_yyy_0_xxyyyyyz_1, g_yyy_0_xxyyyyyzz_1, g_yyy_0_xxyyyyzz_1, g_yyy_0_xxyyyyzzz_1, g_yyy_0_xxyyyzzz_1, g_yyy_0_xxyyyzzzz_1, g_yyy_0_xxyyzzzz_1, g_yyy_0_xxyyzzzzz_1, g_yyy_0_xxyzzzzz_1, g_yyy_0_xxyzzzzzz_1, g_yyy_0_xxzzzzzz_1, g_yyy_0_xxzzzzzzz_1, g_yyy_0_xyyyyyyy_1, g_yyy_0_xyyyyyyyy_1, g_yyy_0_xyyyyyyyz_1, g_yyy_0_xyyyyyyz_1, g_yyy_0_xyyyyyyzz_1, g_yyy_0_xyyyyyzz_1, g_yyy_0_xyyyyyzzz_1, g_yyy_0_xyyyyzzz_1, g_yyy_0_xyyyyzzzz_1, g_yyy_0_xyyyzzzz_1, g_yyy_0_xyyyzzzzz_1, g_yyy_0_xyyzzzzz_1, g_yyy_0_xyyzzzzzz_1, g_yyy_0_xyzzzzzz_1, g_yyy_0_xyzzzzzzz_1, g_yyy_0_xzzzzzzz_1, g_yyy_0_xzzzzzzzz_1, g_yyy_0_yyyyyyyy_1, g_yyy_0_yyyyyyyyy_1, g_yyy_0_yyyyyyyyz_1, g_yyy_0_yyyyyyyz_1, g_yyy_0_yyyyyyyzz_1, g_yyy_0_yyyyyyzz_1, g_yyy_0_yyyyyyzzz_1, g_yyy_0_yyyyyzzz_1, g_yyy_0_yyyyyzzzz_1, g_yyy_0_yyyyzzzz_1, g_yyy_0_yyyyzzzzz_1, g_yyy_0_yyyzzzzz_1, g_yyy_0_yyyzzzzzz_1, g_yyy_0_yyzzzzzz_1, g_yyy_0_yyzzzzzzz_1, g_yyy_0_yzzzzzzz_1, g_yyy_0_yzzzzzzzz_1, g_yyy_0_zzzzzzzz_1, g_yyy_0_zzzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyy_0_xxxxxxxxx_0[i] = 9.0 * g_yyy_0_xxxxxxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxxx_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxxxy_0[i] = 8.0 * g_yyy_0_xxxxxxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxxy_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxxxz_0[i] = 8.0 * g_yyy_0_xxxxxxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxxz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxxyy_0[i] = 7.0 * g_yyy_0_xxxxxxyy_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxyy_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxxyz_0[i] = 7.0 * g_yyy_0_xxxxxxyz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxxzz_0[i] = 7.0 * g_yyy_0_xxxxxxzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxyyy_0[i] = 6.0 * g_yyy_0_xxxxxyyy_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyyy_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxyyz_0[i] = 6.0 * g_yyy_0_xxxxxyyz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxyzz_0[i] = 6.0 * g_yyy_0_xxxxxyzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxzzz_0[i] = 6.0 * g_yyy_0_xxxxxzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxyyyy_0[i] = 5.0 * g_yyy_0_xxxxyyyy_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyyy_1[i] * wa_x[i];

        g_xyyy_0_xxxxxyyyz_0[i] = 5.0 * g_yyy_0_xxxxyyyz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxyyzz_0[i] = 5.0 * g_yyy_0_xxxxyyzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxyzzz_0[i] = 5.0 * g_yyy_0_xxxxyzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxzzzz_0[i] = 5.0 * g_yyy_0_xxxxzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxyyyyy_0[i] = 4.0 * g_yyy_0_xxxyyyyy_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyyy_1[i] * wa_x[i];

        g_xyyy_0_xxxxyyyyz_0[i] = 4.0 * g_yyy_0_xxxyyyyz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xyyy_0_xxxxyyyzz_0[i] = 4.0 * g_yyy_0_xxxyyyzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxyyzzz_0[i] = 4.0 * g_yyy_0_xxxyyzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxyzzzz_0[i] = 4.0 * g_yyy_0_xxxyzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxzzzzz_0[i] = 4.0 * g_yyy_0_xxxzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxyyyyyy_0[i] = 3.0 * g_yyy_0_xxyyyyyy_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyyy_1[i] * wa_x[i];

        g_xyyy_0_xxxyyyyyz_0[i] = 3.0 * g_yyy_0_xxyyyyyz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xyyy_0_xxxyyyyzz_0[i] = 3.0 * g_yyy_0_xxyyyyzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xyyy_0_xxxyyyzzz_0[i] = 3.0 * g_yyy_0_xxyyyzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxyyzzzz_0[i] = 3.0 * g_yyy_0_xxyyzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxyzzzzz_0[i] = 3.0 * g_yyy_0_xxyzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxzzzzzz_0[i] = 3.0 * g_yyy_0_xxzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxyyyyyyy_0[i] = 2.0 * g_yyy_0_xyyyyyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyyy_1[i] * wa_x[i];

        g_xyyy_0_xxyyyyyyz_0[i] = 2.0 * g_yyy_0_xyyyyyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xyyy_0_xxyyyyyzz_0[i] = 2.0 * g_yyy_0_xyyyyyzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xyyy_0_xxyyyyzzz_0[i] = 2.0 * g_yyy_0_xyyyyzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xyyy_0_xxyyyzzzz_0[i] = 2.0 * g_yyy_0_xyyyzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxyyzzzzz_0[i] = 2.0 * g_yyy_0_xyyzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxyzzzzzz_0[i] = 2.0 * g_yyy_0_xyzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxzzzzzzz_0[i] = 2.0 * g_yyy_0_xzzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxzzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xyyyyyyyy_0[i] = g_yyy_0_yyyyyyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyyy_1[i] * wa_x[i];

        g_xyyy_0_xyyyyyyyz_0[i] = g_yyy_0_yyyyyyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xyyy_0_xyyyyyyzz_0[i] = g_yyy_0_yyyyyyzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xyyy_0_xyyyyyzzz_0[i] = g_yyy_0_yyyyyzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xyyy_0_xyyyyzzzz_0[i] = g_yyy_0_yyyyzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xyyy_0_xyyyzzzzz_0[i] = g_yyy_0_yyyzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xyyzzzzzz_0[i] = g_yyy_0_yyzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xyzzzzzzz_0[i] = g_yyy_0_yzzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xzzzzzzzz_0[i] = g_yyy_0_zzzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xzzzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_yyyyyyyyy_0[i] = g_yyy_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xyyy_0_yyyyyyyyz_0[i] = g_yyy_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xyyy_0_yyyyyyyzz_0[i] = g_yyy_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xyyy_0_yyyyyyzzz_0[i] = g_yyy_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xyyy_0_yyyyyzzzz_0[i] = g_yyy_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xyyy_0_yyyyzzzzz_0[i] = g_yyy_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xyyy_0_yyyzzzzzz_0[i] = g_yyy_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_yyzzzzzzz_0[i] = g_yyy_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_yzzzzzzzz_0[i] = g_yyy_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_zzzzzzzzz_0[i] = g_yyy_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 385-440 components of targeted buffer : GSM

    auto g_xyyz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 385);

    auto g_xyyz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 386);

    auto g_xyyz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 387);

    auto g_xyyz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 388);

    auto g_xyyz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 389);

    auto g_xyyz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 390);

    auto g_xyyz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 391);

    auto g_xyyz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 392);

    auto g_xyyz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 393);

    auto g_xyyz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 394);

    auto g_xyyz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 395);

    auto g_xyyz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 396);

    auto g_xyyz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 397);

    auto g_xyyz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 398);

    auto g_xyyz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 399);

    auto g_xyyz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 400);

    auto g_xyyz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 401);

    auto g_xyyz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 402);

    auto g_xyyz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 403);

    auto g_xyyz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 404);

    auto g_xyyz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 405);

    auto g_xyyz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 406);

    auto g_xyyz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 407);

    auto g_xyyz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 408);

    auto g_xyyz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 409);

    auto g_xyyz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 410);

    auto g_xyyz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 411);

    auto g_xyyz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 412);

    auto g_xyyz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 413);

    auto g_xyyz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 414);

    auto g_xyyz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 415);

    auto g_xyyz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 416);

    auto g_xyyz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 417);

    auto g_xyyz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 418);

    auto g_xyyz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 419);

    auto g_xyyz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 420);

    auto g_xyyz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 421);

    auto g_xyyz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 422);

    auto g_xyyz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 423);

    auto g_xyyz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 424);

    auto g_xyyz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 425);

    auto g_xyyz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 426);

    auto g_xyyz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 427);

    auto g_xyyz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 428);

    auto g_xyyz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 429);

    auto g_xyyz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 430);

    auto g_xyyz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 431);

    auto g_xyyz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 432);

    auto g_xyyz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 433);

    auto g_xyyz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 434);

    auto g_xyyz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 435);

    auto g_xyyz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 436);

    auto g_xyyz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 437);

    auto g_xyyz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 438);

    auto g_xyyz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 439);

    #pragma omp simd aligned(g_xyy_0_xxxxxxxxx_1, g_xyy_0_xxxxxxxxy_1, g_xyy_0_xxxxxxxyy_1, g_xyy_0_xxxxxxyyy_1, g_xyy_0_xxxxxyyyy_1, g_xyy_0_xxxxyyyyy_1, g_xyy_0_xxxyyyyyy_1, g_xyy_0_xxyyyyyyy_1, g_xyy_0_xyyyyyyyy_1, g_xyyz_0_xxxxxxxxx_0, g_xyyz_0_xxxxxxxxy_0, g_xyyz_0_xxxxxxxxz_0, g_xyyz_0_xxxxxxxyy_0, g_xyyz_0_xxxxxxxyz_0, g_xyyz_0_xxxxxxxzz_0, g_xyyz_0_xxxxxxyyy_0, g_xyyz_0_xxxxxxyyz_0, g_xyyz_0_xxxxxxyzz_0, g_xyyz_0_xxxxxxzzz_0, g_xyyz_0_xxxxxyyyy_0, g_xyyz_0_xxxxxyyyz_0, g_xyyz_0_xxxxxyyzz_0, g_xyyz_0_xxxxxyzzz_0, g_xyyz_0_xxxxxzzzz_0, g_xyyz_0_xxxxyyyyy_0, g_xyyz_0_xxxxyyyyz_0, g_xyyz_0_xxxxyyyzz_0, g_xyyz_0_xxxxyyzzz_0, g_xyyz_0_xxxxyzzzz_0, g_xyyz_0_xxxxzzzzz_0, g_xyyz_0_xxxyyyyyy_0, g_xyyz_0_xxxyyyyyz_0, g_xyyz_0_xxxyyyyzz_0, g_xyyz_0_xxxyyyzzz_0, g_xyyz_0_xxxyyzzzz_0, g_xyyz_0_xxxyzzzzz_0, g_xyyz_0_xxxzzzzzz_0, g_xyyz_0_xxyyyyyyy_0, g_xyyz_0_xxyyyyyyz_0, g_xyyz_0_xxyyyyyzz_0, g_xyyz_0_xxyyyyzzz_0, g_xyyz_0_xxyyyzzzz_0, g_xyyz_0_xxyyzzzzz_0, g_xyyz_0_xxyzzzzzz_0, g_xyyz_0_xxzzzzzzz_0, g_xyyz_0_xyyyyyyyy_0, g_xyyz_0_xyyyyyyyz_0, g_xyyz_0_xyyyyyyzz_0, g_xyyz_0_xyyyyyzzz_0, g_xyyz_0_xyyyyzzzz_0, g_xyyz_0_xyyyzzzzz_0, g_xyyz_0_xyyzzzzzz_0, g_xyyz_0_xyzzzzzzz_0, g_xyyz_0_xzzzzzzzz_0, g_xyyz_0_yyyyyyyyy_0, g_xyyz_0_yyyyyyyyz_0, g_xyyz_0_yyyyyyyzz_0, g_xyyz_0_yyyyyyzzz_0, g_xyyz_0_yyyyyzzzz_0, g_xyyz_0_yyyyzzzzz_0, g_xyyz_0_yyyzzzzzz_0, g_xyyz_0_yyzzzzzzz_0, g_xyyz_0_yzzzzzzzz_0, g_xyyz_0_zzzzzzzzz_0, g_yyz_0_xxxxxxxxz_1, g_yyz_0_xxxxxxxyz_1, g_yyz_0_xxxxxxxz_1, g_yyz_0_xxxxxxxzz_1, g_yyz_0_xxxxxxyyz_1, g_yyz_0_xxxxxxyz_1, g_yyz_0_xxxxxxyzz_1, g_yyz_0_xxxxxxzz_1, g_yyz_0_xxxxxxzzz_1, g_yyz_0_xxxxxyyyz_1, g_yyz_0_xxxxxyyz_1, g_yyz_0_xxxxxyyzz_1, g_yyz_0_xxxxxyzz_1, g_yyz_0_xxxxxyzzz_1, g_yyz_0_xxxxxzzz_1, g_yyz_0_xxxxxzzzz_1, g_yyz_0_xxxxyyyyz_1, g_yyz_0_xxxxyyyz_1, g_yyz_0_xxxxyyyzz_1, g_yyz_0_xxxxyyzz_1, g_yyz_0_xxxxyyzzz_1, g_yyz_0_xxxxyzzz_1, g_yyz_0_xxxxyzzzz_1, g_yyz_0_xxxxzzzz_1, g_yyz_0_xxxxzzzzz_1, g_yyz_0_xxxyyyyyz_1, g_yyz_0_xxxyyyyz_1, g_yyz_0_xxxyyyyzz_1, g_yyz_0_xxxyyyzz_1, g_yyz_0_xxxyyyzzz_1, g_yyz_0_xxxyyzzz_1, g_yyz_0_xxxyyzzzz_1, g_yyz_0_xxxyzzzz_1, g_yyz_0_xxxyzzzzz_1, g_yyz_0_xxxzzzzz_1, g_yyz_0_xxxzzzzzz_1, g_yyz_0_xxyyyyyyz_1, g_yyz_0_xxyyyyyz_1, g_yyz_0_xxyyyyyzz_1, g_yyz_0_xxyyyyzz_1, g_yyz_0_xxyyyyzzz_1, g_yyz_0_xxyyyzzz_1, g_yyz_0_xxyyyzzzz_1, g_yyz_0_xxyyzzzz_1, g_yyz_0_xxyyzzzzz_1, g_yyz_0_xxyzzzzz_1, g_yyz_0_xxyzzzzzz_1, g_yyz_0_xxzzzzzz_1, g_yyz_0_xxzzzzzzz_1, g_yyz_0_xyyyyyyyz_1, g_yyz_0_xyyyyyyz_1, g_yyz_0_xyyyyyyzz_1, g_yyz_0_xyyyyyzz_1, g_yyz_0_xyyyyyzzz_1, g_yyz_0_xyyyyzzz_1, g_yyz_0_xyyyyzzzz_1, g_yyz_0_xyyyzzzz_1, g_yyz_0_xyyyzzzzz_1, g_yyz_0_xyyzzzzz_1, g_yyz_0_xyyzzzzzz_1, g_yyz_0_xyzzzzzz_1, g_yyz_0_xyzzzzzzz_1, g_yyz_0_xzzzzzzz_1, g_yyz_0_xzzzzzzzz_1, g_yyz_0_yyyyyyyyy_1, g_yyz_0_yyyyyyyyz_1, g_yyz_0_yyyyyyyz_1, g_yyz_0_yyyyyyyzz_1, g_yyz_0_yyyyyyzz_1, g_yyz_0_yyyyyyzzz_1, g_yyz_0_yyyyyzzz_1, g_yyz_0_yyyyyzzzz_1, g_yyz_0_yyyyzzzz_1, g_yyz_0_yyyyzzzzz_1, g_yyz_0_yyyzzzzz_1, g_yyz_0_yyyzzzzzz_1, g_yyz_0_yyzzzzzz_1, g_yyz_0_yyzzzzzzz_1, g_yyz_0_yzzzzzzz_1, g_yyz_0_yzzzzzzzz_1, g_yyz_0_zzzzzzzz_1, g_yyz_0_zzzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyz_0_xxxxxxxxx_0[i] = g_xyy_0_xxxxxxxxx_1[i] * wa_z[i];

        g_xyyz_0_xxxxxxxxy_0[i] = g_xyy_0_xxxxxxxxy_1[i] * wa_z[i];

        g_xyyz_0_xxxxxxxxz_0[i] = 8.0 * g_yyz_0_xxxxxxxz_1[i] * fi_acd_0 + g_yyz_0_xxxxxxxxz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxxxyy_0[i] = g_xyy_0_xxxxxxxyy_1[i] * wa_z[i];

        g_xyyz_0_xxxxxxxyz_0[i] = 7.0 * g_yyz_0_xxxxxxyz_1[i] * fi_acd_0 + g_yyz_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxxxzz_0[i] = 7.0 * g_yyz_0_xxxxxxzz_1[i] * fi_acd_0 + g_yyz_0_xxxxxxxzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxxyyy_0[i] = g_xyy_0_xxxxxxyyy_1[i] * wa_z[i];

        g_xyyz_0_xxxxxxyyz_0[i] = 6.0 * g_yyz_0_xxxxxyyz_1[i] * fi_acd_0 + g_yyz_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxxyzz_0[i] = 6.0 * g_yyz_0_xxxxxyzz_1[i] * fi_acd_0 + g_yyz_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxxzzz_0[i] = 6.0 * g_yyz_0_xxxxxzzz_1[i] * fi_acd_0 + g_yyz_0_xxxxxxzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxyyyy_0[i] = g_xyy_0_xxxxxyyyy_1[i] * wa_z[i];

        g_xyyz_0_xxxxxyyyz_0[i] = 5.0 * g_yyz_0_xxxxyyyz_1[i] * fi_acd_0 + g_yyz_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxyyzz_0[i] = 5.0 * g_yyz_0_xxxxyyzz_1[i] * fi_acd_0 + g_yyz_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxyzzz_0[i] = 5.0 * g_yyz_0_xxxxyzzz_1[i] * fi_acd_0 + g_yyz_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxzzzz_0[i] = 5.0 * g_yyz_0_xxxxzzzz_1[i] * fi_acd_0 + g_yyz_0_xxxxxzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxyyyyy_0[i] = g_xyy_0_xxxxyyyyy_1[i] * wa_z[i];

        g_xyyz_0_xxxxyyyyz_0[i] = 4.0 * g_yyz_0_xxxyyyyz_1[i] * fi_acd_0 + g_yyz_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xyyz_0_xxxxyyyzz_0[i] = 4.0 * g_yyz_0_xxxyyyzz_1[i] * fi_acd_0 + g_yyz_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxyyzzz_0[i] = 4.0 * g_yyz_0_xxxyyzzz_1[i] * fi_acd_0 + g_yyz_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxyzzzz_0[i] = 4.0 * g_yyz_0_xxxyzzzz_1[i] * fi_acd_0 + g_yyz_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxzzzzz_0[i] = 4.0 * g_yyz_0_xxxzzzzz_1[i] * fi_acd_0 + g_yyz_0_xxxxzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxyyyyyy_0[i] = g_xyy_0_xxxyyyyyy_1[i] * wa_z[i];

        g_xyyz_0_xxxyyyyyz_0[i] = 3.0 * g_yyz_0_xxyyyyyz_1[i] * fi_acd_0 + g_yyz_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xyyz_0_xxxyyyyzz_0[i] = 3.0 * g_yyz_0_xxyyyyzz_1[i] * fi_acd_0 + g_yyz_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xyyz_0_xxxyyyzzz_0[i] = 3.0 * g_yyz_0_xxyyyzzz_1[i] * fi_acd_0 + g_yyz_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxyyzzzz_0[i] = 3.0 * g_yyz_0_xxyyzzzz_1[i] * fi_acd_0 + g_yyz_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxyzzzzz_0[i] = 3.0 * g_yyz_0_xxyzzzzz_1[i] * fi_acd_0 + g_yyz_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxzzzzzz_0[i] = 3.0 * g_yyz_0_xxzzzzzz_1[i] * fi_acd_0 + g_yyz_0_xxxzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxyyyyyyy_0[i] = g_xyy_0_xxyyyyyyy_1[i] * wa_z[i];

        g_xyyz_0_xxyyyyyyz_0[i] = 2.0 * g_yyz_0_xyyyyyyz_1[i] * fi_acd_0 + g_yyz_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xyyz_0_xxyyyyyzz_0[i] = 2.0 * g_yyz_0_xyyyyyzz_1[i] * fi_acd_0 + g_yyz_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xyyz_0_xxyyyyzzz_0[i] = 2.0 * g_yyz_0_xyyyyzzz_1[i] * fi_acd_0 + g_yyz_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xyyz_0_xxyyyzzzz_0[i] = 2.0 * g_yyz_0_xyyyzzzz_1[i] * fi_acd_0 + g_yyz_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxyyzzzzz_0[i] = 2.0 * g_yyz_0_xyyzzzzz_1[i] * fi_acd_0 + g_yyz_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxyzzzzzz_0[i] = 2.0 * g_yyz_0_xyzzzzzz_1[i] * fi_acd_0 + g_yyz_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxzzzzzzz_0[i] = 2.0 * g_yyz_0_xzzzzzzz_1[i] * fi_acd_0 + g_yyz_0_xxzzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xyyyyyyyy_0[i] = g_xyy_0_xyyyyyyyy_1[i] * wa_z[i];

        g_xyyz_0_xyyyyyyyz_0[i] = g_yyz_0_yyyyyyyz_1[i] * fi_acd_0 + g_yyz_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xyyz_0_xyyyyyyzz_0[i] = g_yyz_0_yyyyyyzz_1[i] * fi_acd_0 + g_yyz_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xyyz_0_xyyyyyzzz_0[i] = g_yyz_0_yyyyyzzz_1[i] * fi_acd_0 + g_yyz_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xyyz_0_xyyyyzzzz_0[i] = g_yyz_0_yyyyzzzz_1[i] * fi_acd_0 + g_yyz_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xyyz_0_xyyyzzzzz_0[i] = g_yyz_0_yyyzzzzz_1[i] * fi_acd_0 + g_yyz_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xyyzzzzzz_0[i] = g_yyz_0_yyzzzzzz_1[i] * fi_acd_0 + g_yyz_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xyzzzzzzz_0[i] = g_yyz_0_yzzzzzzz_1[i] * fi_acd_0 + g_yyz_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xzzzzzzzz_0[i] = g_yyz_0_zzzzzzzz_1[i] * fi_acd_0 + g_yyz_0_xzzzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_yyyyyyyyy_0[i] = g_yyz_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xyyz_0_yyyyyyyyz_0[i] = g_yyz_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xyyz_0_yyyyyyyzz_0[i] = g_yyz_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xyyz_0_yyyyyyzzz_0[i] = g_yyz_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xyyz_0_yyyyyzzzz_0[i] = g_yyz_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xyyz_0_yyyyzzzzz_0[i] = g_yyz_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xyyz_0_yyyzzzzzz_0[i] = g_yyz_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_yyzzzzzzz_0[i] = g_yyz_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_yzzzzzzzz_0[i] = g_yyz_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_zzzzzzzzz_0[i] = g_yyz_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 440-495 components of targeted buffer : GSM

    auto g_xyzz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 440);

    auto g_xyzz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 441);

    auto g_xyzz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 442);

    auto g_xyzz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 443);

    auto g_xyzz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 444);

    auto g_xyzz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 445);

    auto g_xyzz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 446);

    auto g_xyzz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 447);

    auto g_xyzz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 448);

    auto g_xyzz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 449);

    auto g_xyzz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 450);

    auto g_xyzz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 451);

    auto g_xyzz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 452);

    auto g_xyzz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 453);

    auto g_xyzz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 454);

    auto g_xyzz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 455);

    auto g_xyzz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 456);

    auto g_xyzz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 457);

    auto g_xyzz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 458);

    auto g_xyzz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 459);

    auto g_xyzz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 460);

    auto g_xyzz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 461);

    auto g_xyzz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 462);

    auto g_xyzz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 463);

    auto g_xyzz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 464);

    auto g_xyzz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 465);

    auto g_xyzz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 466);

    auto g_xyzz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 467);

    auto g_xyzz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 468);

    auto g_xyzz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 469);

    auto g_xyzz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 470);

    auto g_xyzz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 471);

    auto g_xyzz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 472);

    auto g_xyzz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 473);

    auto g_xyzz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 474);

    auto g_xyzz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 475);

    auto g_xyzz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 476);

    auto g_xyzz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 477);

    auto g_xyzz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 478);

    auto g_xyzz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 479);

    auto g_xyzz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 480);

    auto g_xyzz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 481);

    auto g_xyzz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 482);

    auto g_xyzz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 483);

    auto g_xyzz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 484);

    auto g_xyzz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 485);

    auto g_xyzz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 486);

    auto g_xyzz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 487);

    auto g_xyzz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 488);

    auto g_xyzz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 489);

    auto g_xyzz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 490);

    auto g_xyzz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 491);

    auto g_xyzz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 492);

    auto g_xyzz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 493);

    auto g_xyzz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 494);

    #pragma omp simd aligned(g_xyzz_0_xxxxxxxxx_0, g_xyzz_0_xxxxxxxxy_0, g_xyzz_0_xxxxxxxxz_0, g_xyzz_0_xxxxxxxyy_0, g_xyzz_0_xxxxxxxyz_0, g_xyzz_0_xxxxxxxzz_0, g_xyzz_0_xxxxxxyyy_0, g_xyzz_0_xxxxxxyyz_0, g_xyzz_0_xxxxxxyzz_0, g_xyzz_0_xxxxxxzzz_0, g_xyzz_0_xxxxxyyyy_0, g_xyzz_0_xxxxxyyyz_0, g_xyzz_0_xxxxxyyzz_0, g_xyzz_0_xxxxxyzzz_0, g_xyzz_0_xxxxxzzzz_0, g_xyzz_0_xxxxyyyyy_0, g_xyzz_0_xxxxyyyyz_0, g_xyzz_0_xxxxyyyzz_0, g_xyzz_0_xxxxyyzzz_0, g_xyzz_0_xxxxyzzzz_0, g_xyzz_0_xxxxzzzzz_0, g_xyzz_0_xxxyyyyyy_0, g_xyzz_0_xxxyyyyyz_0, g_xyzz_0_xxxyyyyzz_0, g_xyzz_0_xxxyyyzzz_0, g_xyzz_0_xxxyyzzzz_0, g_xyzz_0_xxxyzzzzz_0, g_xyzz_0_xxxzzzzzz_0, g_xyzz_0_xxyyyyyyy_0, g_xyzz_0_xxyyyyyyz_0, g_xyzz_0_xxyyyyyzz_0, g_xyzz_0_xxyyyyzzz_0, g_xyzz_0_xxyyyzzzz_0, g_xyzz_0_xxyyzzzzz_0, g_xyzz_0_xxyzzzzzz_0, g_xyzz_0_xxzzzzzzz_0, g_xyzz_0_xyyyyyyyy_0, g_xyzz_0_xyyyyyyyz_0, g_xyzz_0_xyyyyyyzz_0, g_xyzz_0_xyyyyyzzz_0, g_xyzz_0_xyyyyzzzz_0, g_xyzz_0_xyyyzzzzz_0, g_xyzz_0_xyyzzzzzz_0, g_xyzz_0_xyzzzzzzz_0, g_xyzz_0_xzzzzzzzz_0, g_xyzz_0_yyyyyyyyy_0, g_xyzz_0_yyyyyyyyz_0, g_xyzz_0_yyyyyyyzz_0, g_xyzz_0_yyyyyyzzz_0, g_xyzz_0_yyyyyzzzz_0, g_xyzz_0_yyyyzzzzz_0, g_xyzz_0_yyyzzzzzz_0, g_xyzz_0_yyzzzzzzz_0, g_xyzz_0_yzzzzzzzz_0, g_xyzz_0_zzzzzzzzz_0, g_xzz_0_xxxxxxxxx_1, g_xzz_0_xxxxxxxxz_1, g_xzz_0_xxxxxxxzz_1, g_xzz_0_xxxxxxzzz_1, g_xzz_0_xxxxxzzzz_1, g_xzz_0_xxxxzzzzz_1, g_xzz_0_xxxzzzzzz_1, g_xzz_0_xxzzzzzzz_1, g_xzz_0_xzzzzzzzz_1, g_yzz_0_xxxxxxxxy_1, g_yzz_0_xxxxxxxy_1, g_yzz_0_xxxxxxxyy_1, g_yzz_0_xxxxxxxyz_1, g_yzz_0_xxxxxxyy_1, g_yzz_0_xxxxxxyyy_1, g_yzz_0_xxxxxxyyz_1, g_yzz_0_xxxxxxyz_1, g_yzz_0_xxxxxxyzz_1, g_yzz_0_xxxxxyyy_1, g_yzz_0_xxxxxyyyy_1, g_yzz_0_xxxxxyyyz_1, g_yzz_0_xxxxxyyz_1, g_yzz_0_xxxxxyyzz_1, g_yzz_0_xxxxxyzz_1, g_yzz_0_xxxxxyzzz_1, g_yzz_0_xxxxyyyy_1, g_yzz_0_xxxxyyyyy_1, g_yzz_0_xxxxyyyyz_1, g_yzz_0_xxxxyyyz_1, g_yzz_0_xxxxyyyzz_1, g_yzz_0_xxxxyyzz_1, g_yzz_0_xxxxyyzzz_1, g_yzz_0_xxxxyzzz_1, g_yzz_0_xxxxyzzzz_1, g_yzz_0_xxxyyyyy_1, g_yzz_0_xxxyyyyyy_1, g_yzz_0_xxxyyyyyz_1, g_yzz_0_xxxyyyyz_1, g_yzz_0_xxxyyyyzz_1, g_yzz_0_xxxyyyzz_1, g_yzz_0_xxxyyyzzz_1, g_yzz_0_xxxyyzzz_1, g_yzz_0_xxxyyzzzz_1, g_yzz_0_xxxyzzzz_1, g_yzz_0_xxxyzzzzz_1, g_yzz_0_xxyyyyyy_1, g_yzz_0_xxyyyyyyy_1, g_yzz_0_xxyyyyyyz_1, g_yzz_0_xxyyyyyz_1, g_yzz_0_xxyyyyyzz_1, g_yzz_0_xxyyyyzz_1, g_yzz_0_xxyyyyzzz_1, g_yzz_0_xxyyyzzz_1, g_yzz_0_xxyyyzzzz_1, g_yzz_0_xxyyzzzz_1, g_yzz_0_xxyyzzzzz_1, g_yzz_0_xxyzzzzz_1, g_yzz_0_xxyzzzzzz_1, g_yzz_0_xyyyyyyy_1, g_yzz_0_xyyyyyyyy_1, g_yzz_0_xyyyyyyyz_1, g_yzz_0_xyyyyyyz_1, g_yzz_0_xyyyyyyzz_1, g_yzz_0_xyyyyyzz_1, g_yzz_0_xyyyyyzzz_1, g_yzz_0_xyyyyzzz_1, g_yzz_0_xyyyyzzzz_1, g_yzz_0_xyyyzzzz_1, g_yzz_0_xyyyzzzzz_1, g_yzz_0_xyyzzzzz_1, g_yzz_0_xyyzzzzzz_1, g_yzz_0_xyzzzzzz_1, g_yzz_0_xyzzzzzzz_1, g_yzz_0_yyyyyyyy_1, g_yzz_0_yyyyyyyyy_1, g_yzz_0_yyyyyyyyz_1, g_yzz_0_yyyyyyyz_1, g_yzz_0_yyyyyyyzz_1, g_yzz_0_yyyyyyzz_1, g_yzz_0_yyyyyyzzz_1, g_yzz_0_yyyyyzzz_1, g_yzz_0_yyyyyzzzz_1, g_yzz_0_yyyyzzzz_1, g_yzz_0_yyyyzzzzz_1, g_yzz_0_yyyzzzzz_1, g_yzz_0_yyyzzzzzz_1, g_yzz_0_yyzzzzzz_1, g_yzz_0_yyzzzzzzz_1, g_yzz_0_yzzzzzzz_1, g_yzz_0_yzzzzzzzz_1, g_yzz_0_zzzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzz_0_xxxxxxxxx_0[i] = g_xzz_0_xxxxxxxxx_1[i] * wa_y[i];

        g_xyzz_0_xxxxxxxxy_0[i] = 8.0 * g_yzz_0_xxxxxxxy_1[i] * fi_acd_0 + g_yzz_0_xxxxxxxxy_1[i] * wa_x[i];

        g_xyzz_0_xxxxxxxxz_0[i] = g_xzz_0_xxxxxxxxz_1[i] * wa_y[i];

        g_xyzz_0_xxxxxxxyy_0[i] = 7.0 * g_yzz_0_xxxxxxyy_1[i] * fi_acd_0 + g_yzz_0_xxxxxxxyy_1[i] * wa_x[i];

        g_xyzz_0_xxxxxxxyz_0[i] = 7.0 * g_yzz_0_xxxxxxyz_1[i] * fi_acd_0 + g_yzz_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xyzz_0_xxxxxxxzz_0[i] = g_xzz_0_xxxxxxxzz_1[i] * wa_y[i];

        g_xyzz_0_xxxxxxyyy_0[i] = 6.0 * g_yzz_0_xxxxxyyy_1[i] * fi_acd_0 + g_yzz_0_xxxxxxyyy_1[i] * wa_x[i];

        g_xyzz_0_xxxxxxyyz_0[i] = 6.0 * g_yzz_0_xxxxxyyz_1[i] * fi_acd_0 + g_yzz_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xyzz_0_xxxxxxyzz_0[i] = 6.0 * g_yzz_0_xxxxxyzz_1[i] * fi_acd_0 + g_yzz_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xyzz_0_xxxxxxzzz_0[i] = g_xzz_0_xxxxxxzzz_1[i] * wa_y[i];

        g_xyzz_0_xxxxxyyyy_0[i] = 5.0 * g_yzz_0_xxxxyyyy_1[i] * fi_acd_0 + g_yzz_0_xxxxxyyyy_1[i] * wa_x[i];

        g_xyzz_0_xxxxxyyyz_0[i] = 5.0 * g_yzz_0_xxxxyyyz_1[i] * fi_acd_0 + g_yzz_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xyzz_0_xxxxxyyzz_0[i] = 5.0 * g_yzz_0_xxxxyyzz_1[i] * fi_acd_0 + g_yzz_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xyzz_0_xxxxxyzzz_0[i] = 5.0 * g_yzz_0_xxxxyzzz_1[i] * fi_acd_0 + g_yzz_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xyzz_0_xxxxxzzzz_0[i] = g_xzz_0_xxxxxzzzz_1[i] * wa_y[i];

        g_xyzz_0_xxxxyyyyy_0[i] = 4.0 * g_yzz_0_xxxyyyyy_1[i] * fi_acd_0 + g_yzz_0_xxxxyyyyy_1[i] * wa_x[i];

        g_xyzz_0_xxxxyyyyz_0[i] = 4.0 * g_yzz_0_xxxyyyyz_1[i] * fi_acd_0 + g_yzz_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xyzz_0_xxxxyyyzz_0[i] = 4.0 * g_yzz_0_xxxyyyzz_1[i] * fi_acd_0 + g_yzz_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xyzz_0_xxxxyyzzz_0[i] = 4.0 * g_yzz_0_xxxyyzzz_1[i] * fi_acd_0 + g_yzz_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xyzz_0_xxxxyzzzz_0[i] = 4.0 * g_yzz_0_xxxyzzzz_1[i] * fi_acd_0 + g_yzz_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xyzz_0_xxxxzzzzz_0[i] = g_xzz_0_xxxxzzzzz_1[i] * wa_y[i];

        g_xyzz_0_xxxyyyyyy_0[i] = 3.0 * g_yzz_0_xxyyyyyy_1[i] * fi_acd_0 + g_yzz_0_xxxyyyyyy_1[i] * wa_x[i];

        g_xyzz_0_xxxyyyyyz_0[i] = 3.0 * g_yzz_0_xxyyyyyz_1[i] * fi_acd_0 + g_yzz_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xyzz_0_xxxyyyyzz_0[i] = 3.0 * g_yzz_0_xxyyyyzz_1[i] * fi_acd_0 + g_yzz_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xyzz_0_xxxyyyzzz_0[i] = 3.0 * g_yzz_0_xxyyyzzz_1[i] * fi_acd_0 + g_yzz_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xyzz_0_xxxyyzzzz_0[i] = 3.0 * g_yzz_0_xxyyzzzz_1[i] * fi_acd_0 + g_yzz_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xyzz_0_xxxyzzzzz_0[i] = 3.0 * g_yzz_0_xxyzzzzz_1[i] * fi_acd_0 + g_yzz_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xyzz_0_xxxzzzzzz_0[i] = g_xzz_0_xxxzzzzzz_1[i] * wa_y[i];

        g_xyzz_0_xxyyyyyyy_0[i] = 2.0 * g_yzz_0_xyyyyyyy_1[i] * fi_acd_0 + g_yzz_0_xxyyyyyyy_1[i] * wa_x[i];

        g_xyzz_0_xxyyyyyyz_0[i] = 2.0 * g_yzz_0_xyyyyyyz_1[i] * fi_acd_0 + g_yzz_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xyzz_0_xxyyyyyzz_0[i] = 2.0 * g_yzz_0_xyyyyyzz_1[i] * fi_acd_0 + g_yzz_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xyzz_0_xxyyyyzzz_0[i] = 2.0 * g_yzz_0_xyyyyzzz_1[i] * fi_acd_0 + g_yzz_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xyzz_0_xxyyyzzzz_0[i] = 2.0 * g_yzz_0_xyyyzzzz_1[i] * fi_acd_0 + g_yzz_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xyzz_0_xxyyzzzzz_0[i] = 2.0 * g_yzz_0_xyyzzzzz_1[i] * fi_acd_0 + g_yzz_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xyzz_0_xxyzzzzzz_0[i] = 2.0 * g_yzz_0_xyzzzzzz_1[i] * fi_acd_0 + g_yzz_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xyzz_0_xxzzzzzzz_0[i] = g_xzz_0_xxzzzzzzz_1[i] * wa_y[i];

        g_xyzz_0_xyyyyyyyy_0[i] = g_yzz_0_yyyyyyyy_1[i] * fi_acd_0 + g_yzz_0_xyyyyyyyy_1[i] * wa_x[i];

        g_xyzz_0_xyyyyyyyz_0[i] = g_yzz_0_yyyyyyyz_1[i] * fi_acd_0 + g_yzz_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xyzz_0_xyyyyyyzz_0[i] = g_yzz_0_yyyyyyzz_1[i] * fi_acd_0 + g_yzz_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xyzz_0_xyyyyyzzz_0[i] = g_yzz_0_yyyyyzzz_1[i] * fi_acd_0 + g_yzz_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xyzz_0_xyyyyzzzz_0[i] = g_yzz_0_yyyyzzzz_1[i] * fi_acd_0 + g_yzz_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xyzz_0_xyyyzzzzz_0[i] = g_yzz_0_yyyzzzzz_1[i] * fi_acd_0 + g_yzz_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xyzz_0_xyyzzzzzz_0[i] = g_yzz_0_yyzzzzzz_1[i] * fi_acd_0 + g_yzz_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xyzz_0_xyzzzzzzz_0[i] = g_yzz_0_yzzzzzzz_1[i] * fi_acd_0 + g_yzz_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xyzz_0_xzzzzzzzz_0[i] = g_xzz_0_xzzzzzzzz_1[i] * wa_y[i];

        g_xyzz_0_yyyyyyyyy_0[i] = g_yzz_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xyzz_0_yyyyyyyyz_0[i] = g_yzz_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xyzz_0_yyyyyyyzz_0[i] = g_yzz_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xyzz_0_yyyyyyzzz_0[i] = g_yzz_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xyzz_0_yyyyyzzzz_0[i] = g_yzz_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xyzz_0_yyyyzzzzz_0[i] = g_yzz_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xyzz_0_yyyzzzzzz_0[i] = g_yzz_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xyzz_0_yyzzzzzzz_0[i] = g_yzz_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xyzz_0_yzzzzzzzz_0[i] = g_yzz_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xyzz_0_zzzzzzzzz_0[i] = g_yzz_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 495-550 components of targeted buffer : GSM

    auto g_xzzz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 495);

    auto g_xzzz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 496);

    auto g_xzzz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 497);

    auto g_xzzz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 498);

    auto g_xzzz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 499);

    auto g_xzzz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 500);

    auto g_xzzz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 501);

    auto g_xzzz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 502);

    auto g_xzzz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 503);

    auto g_xzzz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 504);

    auto g_xzzz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 505);

    auto g_xzzz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 506);

    auto g_xzzz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 507);

    auto g_xzzz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 508);

    auto g_xzzz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 509);

    auto g_xzzz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 510);

    auto g_xzzz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 511);

    auto g_xzzz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 512);

    auto g_xzzz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 513);

    auto g_xzzz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 514);

    auto g_xzzz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 515);

    auto g_xzzz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 516);

    auto g_xzzz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 517);

    auto g_xzzz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 518);

    auto g_xzzz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 519);

    auto g_xzzz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 520);

    auto g_xzzz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 521);

    auto g_xzzz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 522);

    auto g_xzzz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 523);

    auto g_xzzz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 524);

    auto g_xzzz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 525);

    auto g_xzzz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 526);

    auto g_xzzz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 527);

    auto g_xzzz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 528);

    auto g_xzzz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 529);

    auto g_xzzz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 530);

    auto g_xzzz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 531);

    auto g_xzzz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 532);

    auto g_xzzz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 533);

    auto g_xzzz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 534);

    auto g_xzzz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 535);

    auto g_xzzz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 536);

    auto g_xzzz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 537);

    auto g_xzzz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 538);

    auto g_xzzz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 539);

    auto g_xzzz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 540);

    auto g_xzzz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 541);

    auto g_xzzz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 542);

    auto g_xzzz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 543);

    auto g_xzzz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 544);

    auto g_xzzz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 545);

    auto g_xzzz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 546);

    auto g_xzzz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 547);

    auto g_xzzz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 548);

    auto g_xzzz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 549);

    #pragma omp simd aligned(g_xzzz_0_xxxxxxxxx_0, g_xzzz_0_xxxxxxxxy_0, g_xzzz_0_xxxxxxxxz_0, g_xzzz_0_xxxxxxxyy_0, g_xzzz_0_xxxxxxxyz_0, g_xzzz_0_xxxxxxxzz_0, g_xzzz_0_xxxxxxyyy_0, g_xzzz_0_xxxxxxyyz_0, g_xzzz_0_xxxxxxyzz_0, g_xzzz_0_xxxxxxzzz_0, g_xzzz_0_xxxxxyyyy_0, g_xzzz_0_xxxxxyyyz_0, g_xzzz_0_xxxxxyyzz_0, g_xzzz_0_xxxxxyzzz_0, g_xzzz_0_xxxxxzzzz_0, g_xzzz_0_xxxxyyyyy_0, g_xzzz_0_xxxxyyyyz_0, g_xzzz_0_xxxxyyyzz_0, g_xzzz_0_xxxxyyzzz_0, g_xzzz_0_xxxxyzzzz_0, g_xzzz_0_xxxxzzzzz_0, g_xzzz_0_xxxyyyyyy_0, g_xzzz_0_xxxyyyyyz_0, g_xzzz_0_xxxyyyyzz_0, g_xzzz_0_xxxyyyzzz_0, g_xzzz_0_xxxyyzzzz_0, g_xzzz_0_xxxyzzzzz_0, g_xzzz_0_xxxzzzzzz_0, g_xzzz_0_xxyyyyyyy_0, g_xzzz_0_xxyyyyyyz_0, g_xzzz_0_xxyyyyyzz_0, g_xzzz_0_xxyyyyzzz_0, g_xzzz_0_xxyyyzzzz_0, g_xzzz_0_xxyyzzzzz_0, g_xzzz_0_xxyzzzzzz_0, g_xzzz_0_xxzzzzzzz_0, g_xzzz_0_xyyyyyyyy_0, g_xzzz_0_xyyyyyyyz_0, g_xzzz_0_xyyyyyyzz_0, g_xzzz_0_xyyyyyzzz_0, g_xzzz_0_xyyyyzzzz_0, g_xzzz_0_xyyyzzzzz_0, g_xzzz_0_xyyzzzzzz_0, g_xzzz_0_xyzzzzzzz_0, g_xzzz_0_xzzzzzzzz_0, g_xzzz_0_yyyyyyyyy_0, g_xzzz_0_yyyyyyyyz_0, g_xzzz_0_yyyyyyyzz_0, g_xzzz_0_yyyyyyzzz_0, g_xzzz_0_yyyyyzzzz_0, g_xzzz_0_yyyyzzzzz_0, g_xzzz_0_yyyzzzzzz_0, g_xzzz_0_yyzzzzzzz_0, g_xzzz_0_yzzzzzzzz_0, g_xzzz_0_zzzzzzzzz_0, g_zzz_0_xxxxxxxx_1, g_zzz_0_xxxxxxxxx_1, g_zzz_0_xxxxxxxxy_1, g_zzz_0_xxxxxxxxz_1, g_zzz_0_xxxxxxxy_1, g_zzz_0_xxxxxxxyy_1, g_zzz_0_xxxxxxxyz_1, g_zzz_0_xxxxxxxz_1, g_zzz_0_xxxxxxxzz_1, g_zzz_0_xxxxxxyy_1, g_zzz_0_xxxxxxyyy_1, g_zzz_0_xxxxxxyyz_1, g_zzz_0_xxxxxxyz_1, g_zzz_0_xxxxxxyzz_1, g_zzz_0_xxxxxxzz_1, g_zzz_0_xxxxxxzzz_1, g_zzz_0_xxxxxyyy_1, g_zzz_0_xxxxxyyyy_1, g_zzz_0_xxxxxyyyz_1, g_zzz_0_xxxxxyyz_1, g_zzz_0_xxxxxyyzz_1, g_zzz_0_xxxxxyzz_1, g_zzz_0_xxxxxyzzz_1, g_zzz_0_xxxxxzzz_1, g_zzz_0_xxxxxzzzz_1, g_zzz_0_xxxxyyyy_1, g_zzz_0_xxxxyyyyy_1, g_zzz_0_xxxxyyyyz_1, g_zzz_0_xxxxyyyz_1, g_zzz_0_xxxxyyyzz_1, g_zzz_0_xxxxyyzz_1, g_zzz_0_xxxxyyzzz_1, g_zzz_0_xxxxyzzz_1, g_zzz_0_xxxxyzzzz_1, g_zzz_0_xxxxzzzz_1, g_zzz_0_xxxxzzzzz_1, g_zzz_0_xxxyyyyy_1, g_zzz_0_xxxyyyyyy_1, g_zzz_0_xxxyyyyyz_1, g_zzz_0_xxxyyyyz_1, g_zzz_0_xxxyyyyzz_1, g_zzz_0_xxxyyyzz_1, g_zzz_0_xxxyyyzzz_1, g_zzz_0_xxxyyzzz_1, g_zzz_0_xxxyyzzzz_1, g_zzz_0_xxxyzzzz_1, g_zzz_0_xxxyzzzzz_1, g_zzz_0_xxxzzzzz_1, g_zzz_0_xxxzzzzzz_1, g_zzz_0_xxyyyyyy_1, g_zzz_0_xxyyyyyyy_1, g_zzz_0_xxyyyyyyz_1, g_zzz_0_xxyyyyyz_1, g_zzz_0_xxyyyyyzz_1, g_zzz_0_xxyyyyzz_1, g_zzz_0_xxyyyyzzz_1, g_zzz_0_xxyyyzzz_1, g_zzz_0_xxyyyzzzz_1, g_zzz_0_xxyyzzzz_1, g_zzz_0_xxyyzzzzz_1, g_zzz_0_xxyzzzzz_1, g_zzz_0_xxyzzzzzz_1, g_zzz_0_xxzzzzzz_1, g_zzz_0_xxzzzzzzz_1, g_zzz_0_xyyyyyyy_1, g_zzz_0_xyyyyyyyy_1, g_zzz_0_xyyyyyyyz_1, g_zzz_0_xyyyyyyz_1, g_zzz_0_xyyyyyyzz_1, g_zzz_0_xyyyyyzz_1, g_zzz_0_xyyyyyzzz_1, g_zzz_0_xyyyyzzz_1, g_zzz_0_xyyyyzzzz_1, g_zzz_0_xyyyzzzz_1, g_zzz_0_xyyyzzzzz_1, g_zzz_0_xyyzzzzz_1, g_zzz_0_xyyzzzzzz_1, g_zzz_0_xyzzzzzz_1, g_zzz_0_xyzzzzzzz_1, g_zzz_0_xzzzzzzz_1, g_zzz_0_xzzzzzzzz_1, g_zzz_0_yyyyyyyy_1, g_zzz_0_yyyyyyyyy_1, g_zzz_0_yyyyyyyyz_1, g_zzz_0_yyyyyyyz_1, g_zzz_0_yyyyyyyzz_1, g_zzz_0_yyyyyyzz_1, g_zzz_0_yyyyyyzzz_1, g_zzz_0_yyyyyzzz_1, g_zzz_0_yyyyyzzzz_1, g_zzz_0_yyyyzzzz_1, g_zzz_0_yyyyzzzzz_1, g_zzz_0_yyyzzzzz_1, g_zzz_0_yyyzzzzzz_1, g_zzz_0_yyzzzzzz_1, g_zzz_0_yyzzzzzzz_1, g_zzz_0_yzzzzzzz_1, g_zzz_0_yzzzzzzzz_1, g_zzz_0_zzzzzzzz_1, g_zzz_0_zzzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzz_0_xxxxxxxxx_0[i] = 9.0 * g_zzz_0_xxxxxxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxxx_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxxxy_0[i] = 8.0 * g_zzz_0_xxxxxxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxxy_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxxxz_0[i] = 8.0 * g_zzz_0_xxxxxxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxxz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxxyy_0[i] = 7.0 * g_zzz_0_xxxxxxyy_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxyy_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxxyz_0[i] = 7.0 * g_zzz_0_xxxxxxyz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxyz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxxzz_0[i] = 7.0 * g_zzz_0_xxxxxxzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxyyy_0[i] = 6.0 * g_zzz_0_xxxxxyyy_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyyy_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxyyz_0[i] = 6.0 * g_zzz_0_xxxxxyyz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyyz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxyzz_0[i] = 6.0 * g_zzz_0_xxxxxyzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxzzz_0[i] = 6.0 * g_zzz_0_xxxxxzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxyyyy_0[i] = 5.0 * g_zzz_0_xxxxyyyy_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyyy_1[i] * wa_x[i];

        g_xzzz_0_xxxxxyyyz_0[i] = 5.0 * g_zzz_0_xxxxyyyz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyyz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxyyzz_0[i] = 5.0 * g_zzz_0_xxxxyyzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxyzzz_0[i] = 5.0 * g_zzz_0_xxxxyzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxzzzz_0[i] = 5.0 * g_zzz_0_xxxxzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxyyyyy_0[i] = 4.0 * g_zzz_0_xxxyyyyy_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyyy_1[i] * wa_x[i];

        g_xzzz_0_xxxxyyyyz_0[i] = 4.0 * g_zzz_0_xxxyyyyz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyyz_1[i] * wa_x[i];

        g_xzzz_0_xxxxyyyzz_0[i] = 4.0 * g_zzz_0_xxxyyyzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxyyzzz_0[i] = 4.0 * g_zzz_0_xxxyyzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxyzzzz_0[i] = 4.0 * g_zzz_0_xxxyzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxzzzzz_0[i] = 4.0 * g_zzz_0_xxxzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxyyyyyy_0[i] = 3.0 * g_zzz_0_xxyyyyyy_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyyy_1[i] * wa_x[i];

        g_xzzz_0_xxxyyyyyz_0[i] = 3.0 * g_zzz_0_xxyyyyyz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyyz_1[i] * wa_x[i];

        g_xzzz_0_xxxyyyyzz_0[i] = 3.0 * g_zzz_0_xxyyyyzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyzz_1[i] * wa_x[i];

        g_xzzz_0_xxxyyyzzz_0[i] = 3.0 * g_zzz_0_xxyyyzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxyyzzzz_0[i] = 3.0 * g_zzz_0_xxyyzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxyzzzzz_0[i] = 3.0 * g_zzz_0_xxyzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxzzzzzz_0[i] = 3.0 * g_zzz_0_xxzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxyyyyyyy_0[i] = 2.0 * g_zzz_0_xyyyyyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyyy_1[i] * wa_x[i];

        g_xzzz_0_xxyyyyyyz_0[i] = 2.0 * g_zzz_0_xyyyyyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyyz_1[i] * wa_x[i];

        g_xzzz_0_xxyyyyyzz_0[i] = 2.0 * g_zzz_0_xyyyyyzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyzz_1[i] * wa_x[i];

        g_xzzz_0_xxyyyyzzz_0[i] = 2.0 * g_zzz_0_xyyyyzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyzzz_1[i] * wa_x[i];

        g_xzzz_0_xxyyyzzzz_0[i] = 2.0 * g_zzz_0_xyyyzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxyyzzzzz_0[i] = 2.0 * g_zzz_0_xyyzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxyzzzzzz_0[i] = 2.0 * g_zzz_0_xyzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxzzzzzzz_0[i] = 2.0 * g_zzz_0_xzzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxzzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xyyyyyyyy_0[i] = g_zzz_0_yyyyyyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyyy_1[i] * wa_x[i];

        g_xzzz_0_xyyyyyyyz_0[i] = g_zzz_0_yyyyyyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyyz_1[i] * wa_x[i];

        g_xzzz_0_xyyyyyyzz_0[i] = g_zzz_0_yyyyyyzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyzz_1[i] * wa_x[i];

        g_xzzz_0_xyyyyyzzz_0[i] = g_zzz_0_yyyyyzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyzzz_1[i] * wa_x[i];

        g_xzzz_0_xyyyyzzzz_0[i] = g_zzz_0_yyyyzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyzzzz_1[i] * wa_x[i];

        g_xzzz_0_xyyyzzzzz_0[i] = g_zzz_0_yyyzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xyyzzzzzz_0[i] = g_zzz_0_yyzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xyzzzzzzz_0[i] = g_zzz_0_yzzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xzzzzzzzz_0[i] = g_zzz_0_zzzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xzzzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_yyyyyyyyy_0[i] = g_zzz_0_yyyyyyyyy_1[i] * wa_x[i];

        g_xzzz_0_yyyyyyyyz_0[i] = g_zzz_0_yyyyyyyyz_1[i] * wa_x[i];

        g_xzzz_0_yyyyyyyzz_0[i] = g_zzz_0_yyyyyyyzz_1[i] * wa_x[i];

        g_xzzz_0_yyyyyyzzz_0[i] = g_zzz_0_yyyyyyzzz_1[i] * wa_x[i];

        g_xzzz_0_yyyyyzzzz_0[i] = g_zzz_0_yyyyyzzzz_1[i] * wa_x[i];

        g_xzzz_0_yyyyzzzzz_0[i] = g_zzz_0_yyyyzzzzz_1[i] * wa_x[i];

        g_xzzz_0_yyyzzzzzz_0[i] = g_zzz_0_yyyzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_yyzzzzzzz_0[i] = g_zzz_0_yyzzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_yzzzzzzzz_0[i] = g_zzz_0_yzzzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_zzzzzzzzz_0[i] = g_zzz_0_zzzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 550-605 components of targeted buffer : GSM

    auto g_yyyy_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 550);

    auto g_yyyy_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 551);

    auto g_yyyy_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 552);

    auto g_yyyy_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 553);

    auto g_yyyy_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 554);

    auto g_yyyy_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 555);

    auto g_yyyy_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 556);

    auto g_yyyy_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 557);

    auto g_yyyy_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 558);

    auto g_yyyy_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 559);

    auto g_yyyy_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 560);

    auto g_yyyy_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 561);

    auto g_yyyy_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 562);

    auto g_yyyy_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 563);

    auto g_yyyy_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 564);

    auto g_yyyy_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 565);

    auto g_yyyy_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 566);

    auto g_yyyy_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 567);

    auto g_yyyy_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 568);

    auto g_yyyy_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 569);

    auto g_yyyy_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 570);

    auto g_yyyy_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 571);

    auto g_yyyy_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 572);

    auto g_yyyy_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 573);

    auto g_yyyy_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 574);

    auto g_yyyy_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 575);

    auto g_yyyy_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 576);

    auto g_yyyy_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 577);

    auto g_yyyy_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 578);

    auto g_yyyy_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 579);

    auto g_yyyy_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 580);

    auto g_yyyy_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 581);

    auto g_yyyy_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 582);

    auto g_yyyy_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 583);

    auto g_yyyy_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 584);

    auto g_yyyy_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 585);

    auto g_yyyy_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 586);

    auto g_yyyy_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 587);

    auto g_yyyy_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 588);

    auto g_yyyy_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 589);

    auto g_yyyy_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 590);

    auto g_yyyy_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 591);

    auto g_yyyy_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 592);

    auto g_yyyy_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 593);

    auto g_yyyy_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 594);

    auto g_yyyy_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 595);

    auto g_yyyy_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 596);

    auto g_yyyy_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 597);

    auto g_yyyy_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 598);

    auto g_yyyy_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 599);

    auto g_yyyy_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 600);

    auto g_yyyy_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 601);

    auto g_yyyy_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 602);

    auto g_yyyy_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 603);

    auto g_yyyy_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 604);

    #pragma omp simd aligned(g_yy_0_xxxxxxxxx_0, g_yy_0_xxxxxxxxx_1, g_yy_0_xxxxxxxxy_0, g_yy_0_xxxxxxxxy_1, g_yy_0_xxxxxxxxz_0, g_yy_0_xxxxxxxxz_1, g_yy_0_xxxxxxxyy_0, g_yy_0_xxxxxxxyy_1, g_yy_0_xxxxxxxyz_0, g_yy_0_xxxxxxxyz_1, g_yy_0_xxxxxxxzz_0, g_yy_0_xxxxxxxzz_1, g_yy_0_xxxxxxyyy_0, g_yy_0_xxxxxxyyy_1, g_yy_0_xxxxxxyyz_0, g_yy_0_xxxxxxyyz_1, g_yy_0_xxxxxxyzz_0, g_yy_0_xxxxxxyzz_1, g_yy_0_xxxxxxzzz_0, g_yy_0_xxxxxxzzz_1, g_yy_0_xxxxxyyyy_0, g_yy_0_xxxxxyyyy_1, g_yy_0_xxxxxyyyz_0, g_yy_0_xxxxxyyyz_1, g_yy_0_xxxxxyyzz_0, g_yy_0_xxxxxyyzz_1, g_yy_0_xxxxxyzzz_0, g_yy_0_xxxxxyzzz_1, g_yy_0_xxxxxzzzz_0, g_yy_0_xxxxxzzzz_1, g_yy_0_xxxxyyyyy_0, g_yy_0_xxxxyyyyy_1, g_yy_0_xxxxyyyyz_0, g_yy_0_xxxxyyyyz_1, g_yy_0_xxxxyyyzz_0, g_yy_0_xxxxyyyzz_1, g_yy_0_xxxxyyzzz_0, g_yy_0_xxxxyyzzz_1, g_yy_0_xxxxyzzzz_0, g_yy_0_xxxxyzzzz_1, g_yy_0_xxxxzzzzz_0, g_yy_0_xxxxzzzzz_1, g_yy_0_xxxyyyyyy_0, g_yy_0_xxxyyyyyy_1, g_yy_0_xxxyyyyyz_0, g_yy_0_xxxyyyyyz_1, g_yy_0_xxxyyyyzz_0, g_yy_0_xxxyyyyzz_1, g_yy_0_xxxyyyzzz_0, g_yy_0_xxxyyyzzz_1, g_yy_0_xxxyyzzzz_0, g_yy_0_xxxyyzzzz_1, g_yy_0_xxxyzzzzz_0, g_yy_0_xxxyzzzzz_1, g_yy_0_xxxzzzzzz_0, g_yy_0_xxxzzzzzz_1, g_yy_0_xxyyyyyyy_0, g_yy_0_xxyyyyyyy_1, g_yy_0_xxyyyyyyz_0, g_yy_0_xxyyyyyyz_1, g_yy_0_xxyyyyyzz_0, g_yy_0_xxyyyyyzz_1, g_yy_0_xxyyyyzzz_0, g_yy_0_xxyyyyzzz_1, g_yy_0_xxyyyzzzz_0, g_yy_0_xxyyyzzzz_1, g_yy_0_xxyyzzzzz_0, g_yy_0_xxyyzzzzz_1, g_yy_0_xxyzzzzzz_0, g_yy_0_xxyzzzzzz_1, g_yy_0_xxzzzzzzz_0, g_yy_0_xxzzzzzzz_1, g_yy_0_xyyyyyyyy_0, g_yy_0_xyyyyyyyy_1, g_yy_0_xyyyyyyyz_0, g_yy_0_xyyyyyyyz_1, g_yy_0_xyyyyyyzz_0, g_yy_0_xyyyyyyzz_1, g_yy_0_xyyyyyzzz_0, g_yy_0_xyyyyyzzz_1, g_yy_0_xyyyyzzzz_0, g_yy_0_xyyyyzzzz_1, g_yy_0_xyyyzzzzz_0, g_yy_0_xyyyzzzzz_1, g_yy_0_xyyzzzzzz_0, g_yy_0_xyyzzzzzz_1, g_yy_0_xyzzzzzzz_0, g_yy_0_xyzzzzzzz_1, g_yy_0_xzzzzzzzz_0, g_yy_0_xzzzzzzzz_1, g_yy_0_yyyyyyyyy_0, g_yy_0_yyyyyyyyy_1, g_yy_0_yyyyyyyyz_0, g_yy_0_yyyyyyyyz_1, g_yy_0_yyyyyyyzz_0, g_yy_0_yyyyyyyzz_1, g_yy_0_yyyyyyzzz_0, g_yy_0_yyyyyyzzz_1, g_yy_0_yyyyyzzzz_0, g_yy_0_yyyyyzzzz_1, g_yy_0_yyyyzzzzz_0, g_yy_0_yyyyzzzzz_1, g_yy_0_yyyzzzzzz_0, g_yy_0_yyyzzzzzz_1, g_yy_0_yyzzzzzzz_0, g_yy_0_yyzzzzzzz_1, g_yy_0_yzzzzzzzz_0, g_yy_0_yzzzzzzzz_1, g_yy_0_zzzzzzzzz_0, g_yy_0_zzzzzzzzz_1, g_yyy_0_xxxxxxxx_1, g_yyy_0_xxxxxxxxx_1, g_yyy_0_xxxxxxxxy_1, g_yyy_0_xxxxxxxxz_1, g_yyy_0_xxxxxxxy_1, g_yyy_0_xxxxxxxyy_1, g_yyy_0_xxxxxxxyz_1, g_yyy_0_xxxxxxxz_1, g_yyy_0_xxxxxxxzz_1, g_yyy_0_xxxxxxyy_1, g_yyy_0_xxxxxxyyy_1, g_yyy_0_xxxxxxyyz_1, g_yyy_0_xxxxxxyz_1, g_yyy_0_xxxxxxyzz_1, g_yyy_0_xxxxxxzz_1, g_yyy_0_xxxxxxzzz_1, g_yyy_0_xxxxxyyy_1, g_yyy_0_xxxxxyyyy_1, g_yyy_0_xxxxxyyyz_1, g_yyy_0_xxxxxyyz_1, g_yyy_0_xxxxxyyzz_1, g_yyy_0_xxxxxyzz_1, g_yyy_0_xxxxxyzzz_1, g_yyy_0_xxxxxzzz_1, g_yyy_0_xxxxxzzzz_1, g_yyy_0_xxxxyyyy_1, g_yyy_0_xxxxyyyyy_1, g_yyy_0_xxxxyyyyz_1, g_yyy_0_xxxxyyyz_1, g_yyy_0_xxxxyyyzz_1, g_yyy_0_xxxxyyzz_1, g_yyy_0_xxxxyyzzz_1, g_yyy_0_xxxxyzzz_1, g_yyy_0_xxxxyzzzz_1, g_yyy_0_xxxxzzzz_1, g_yyy_0_xxxxzzzzz_1, g_yyy_0_xxxyyyyy_1, g_yyy_0_xxxyyyyyy_1, g_yyy_0_xxxyyyyyz_1, g_yyy_0_xxxyyyyz_1, g_yyy_0_xxxyyyyzz_1, g_yyy_0_xxxyyyzz_1, g_yyy_0_xxxyyyzzz_1, g_yyy_0_xxxyyzzz_1, g_yyy_0_xxxyyzzzz_1, g_yyy_0_xxxyzzzz_1, g_yyy_0_xxxyzzzzz_1, g_yyy_0_xxxzzzzz_1, g_yyy_0_xxxzzzzzz_1, g_yyy_0_xxyyyyyy_1, g_yyy_0_xxyyyyyyy_1, g_yyy_0_xxyyyyyyz_1, g_yyy_0_xxyyyyyz_1, g_yyy_0_xxyyyyyzz_1, g_yyy_0_xxyyyyzz_1, g_yyy_0_xxyyyyzzz_1, g_yyy_0_xxyyyzzz_1, g_yyy_0_xxyyyzzzz_1, g_yyy_0_xxyyzzzz_1, g_yyy_0_xxyyzzzzz_1, g_yyy_0_xxyzzzzz_1, g_yyy_0_xxyzzzzzz_1, g_yyy_0_xxzzzzzz_1, g_yyy_0_xxzzzzzzz_1, g_yyy_0_xyyyyyyy_1, g_yyy_0_xyyyyyyyy_1, g_yyy_0_xyyyyyyyz_1, g_yyy_0_xyyyyyyz_1, g_yyy_0_xyyyyyyzz_1, g_yyy_0_xyyyyyzz_1, g_yyy_0_xyyyyyzzz_1, g_yyy_0_xyyyyzzz_1, g_yyy_0_xyyyyzzzz_1, g_yyy_0_xyyyzzzz_1, g_yyy_0_xyyyzzzzz_1, g_yyy_0_xyyzzzzz_1, g_yyy_0_xyyzzzzzz_1, g_yyy_0_xyzzzzzz_1, g_yyy_0_xyzzzzzzz_1, g_yyy_0_xzzzzzzz_1, g_yyy_0_xzzzzzzzz_1, g_yyy_0_yyyyyyyy_1, g_yyy_0_yyyyyyyyy_1, g_yyy_0_yyyyyyyyz_1, g_yyy_0_yyyyyyyz_1, g_yyy_0_yyyyyyyzz_1, g_yyy_0_yyyyyyzz_1, g_yyy_0_yyyyyyzzz_1, g_yyy_0_yyyyyzzz_1, g_yyy_0_yyyyyzzzz_1, g_yyy_0_yyyyzzzz_1, g_yyy_0_yyyyzzzzz_1, g_yyy_0_yyyzzzzz_1, g_yyy_0_yyyzzzzzz_1, g_yyy_0_yyzzzzzz_1, g_yyy_0_yyzzzzzzz_1, g_yyy_0_yzzzzzzz_1, g_yyy_0_yzzzzzzzz_1, g_yyy_0_zzzzzzzz_1, g_yyy_0_zzzzzzzzz_1, g_yyyy_0_xxxxxxxxx_0, g_yyyy_0_xxxxxxxxy_0, g_yyyy_0_xxxxxxxxz_0, g_yyyy_0_xxxxxxxyy_0, g_yyyy_0_xxxxxxxyz_0, g_yyyy_0_xxxxxxxzz_0, g_yyyy_0_xxxxxxyyy_0, g_yyyy_0_xxxxxxyyz_0, g_yyyy_0_xxxxxxyzz_0, g_yyyy_0_xxxxxxzzz_0, g_yyyy_0_xxxxxyyyy_0, g_yyyy_0_xxxxxyyyz_0, g_yyyy_0_xxxxxyyzz_0, g_yyyy_0_xxxxxyzzz_0, g_yyyy_0_xxxxxzzzz_0, g_yyyy_0_xxxxyyyyy_0, g_yyyy_0_xxxxyyyyz_0, g_yyyy_0_xxxxyyyzz_0, g_yyyy_0_xxxxyyzzz_0, g_yyyy_0_xxxxyzzzz_0, g_yyyy_0_xxxxzzzzz_0, g_yyyy_0_xxxyyyyyy_0, g_yyyy_0_xxxyyyyyz_0, g_yyyy_0_xxxyyyyzz_0, g_yyyy_0_xxxyyyzzz_0, g_yyyy_0_xxxyyzzzz_0, g_yyyy_0_xxxyzzzzz_0, g_yyyy_0_xxxzzzzzz_0, g_yyyy_0_xxyyyyyyy_0, g_yyyy_0_xxyyyyyyz_0, g_yyyy_0_xxyyyyyzz_0, g_yyyy_0_xxyyyyzzz_0, g_yyyy_0_xxyyyzzzz_0, g_yyyy_0_xxyyzzzzz_0, g_yyyy_0_xxyzzzzzz_0, g_yyyy_0_xxzzzzzzz_0, g_yyyy_0_xyyyyyyyy_0, g_yyyy_0_xyyyyyyyz_0, g_yyyy_0_xyyyyyyzz_0, g_yyyy_0_xyyyyyzzz_0, g_yyyy_0_xyyyyzzzz_0, g_yyyy_0_xyyyzzzzz_0, g_yyyy_0_xyyzzzzzz_0, g_yyyy_0_xyzzzzzzz_0, g_yyyy_0_xzzzzzzzz_0, g_yyyy_0_yyyyyyyyy_0, g_yyyy_0_yyyyyyyyz_0, g_yyyy_0_yyyyyyyzz_0, g_yyyy_0_yyyyyyzzz_0, g_yyyy_0_yyyyyzzzz_0, g_yyyy_0_yyyyzzzzz_0, g_yyyy_0_yyyzzzzzz_0, g_yyyy_0_yyzzzzzzz_0, g_yyyy_0_yzzzzzzzz_0, g_yyyy_0_zzzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyy_0_xxxxxxxxx_0[i] = 3.0 * g_yy_0_xxxxxxxxx_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxxxx_1[i] * fz_be_0 + g_yyy_0_xxxxxxxxx_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxxxy_0[i] = 3.0 * g_yy_0_xxxxxxxxy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxxxy_1[i] * fz_be_0 + g_yyy_0_xxxxxxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxxy_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxxxz_0[i] = 3.0 * g_yy_0_xxxxxxxxz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxxxz_1[i] * fz_be_0 + g_yyy_0_xxxxxxxxz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxxyy_0[i] = 3.0 * g_yy_0_xxxxxxxyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxxyy_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxxxxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxyy_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxxyz_0[i] = 3.0 * g_yy_0_xxxxxxxyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxxyz_1[i] * fz_be_0 + g_yyy_0_xxxxxxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxyz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxxzz_0[i] = 3.0 * g_yy_0_xxxxxxxzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxxzz_1[i] * fz_be_0 + g_yyy_0_xxxxxxxzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxyyy_0[i] = 3.0 * g_yy_0_xxxxxxyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxyyy_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxxxxxyy_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyyy_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxyyz_0[i] = 3.0 * g_yy_0_xxxxxxyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxxxxyz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyyz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxyzz_0[i] = 3.0 * g_yy_0_xxxxxxyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxyzz_1[i] * fz_be_0 + g_yyy_0_xxxxxxzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxzzz_0[i] = 3.0 * g_yy_0_xxxxxxzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxzzz_1[i] * fz_be_0 + g_yyy_0_xxxxxxzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxyyyy_0[i] = 3.0 * g_yy_0_xxxxxyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxyyyy_1[i] * fz_be_0 + 4.0 * g_yyy_0_xxxxxyyy_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyyy_1[i] * wa_y[i];

        g_yyyy_0_xxxxxyyyz_0[i] = 3.0 * g_yy_0_xxxxxyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxxxxyyz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyyz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxyyzz_0[i] = 3.0 * g_yy_0_xxxxxyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxxxyzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxyzzz_0[i] = 3.0 * g_yy_0_xxxxxyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxyzzz_1[i] * fz_be_0 + g_yyy_0_xxxxxzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxzzzz_0[i] = 3.0 * g_yy_0_xxxxxzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxzzzz_1[i] * fz_be_0 + g_yyy_0_xxxxxzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxyyyyy_0[i] = 3.0 * g_yy_0_xxxxyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyyyyy_1[i] * fz_be_0 + 5.0 * g_yyy_0_xxxxyyyy_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyyy_1[i] * wa_y[i];

        g_yyyy_0_xxxxyyyyz_0[i] = 3.0 * g_yy_0_xxxxyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyy_0_xxxxyyyz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyyz_1[i] * wa_y[i];

        g_yyyy_0_xxxxyyyzz_0[i] = 3.0 * g_yy_0_xxxxyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxxxyyzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxyyzzz_0[i] = 3.0 * g_yy_0_xxxxyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxxyzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxyzzzz_0[i] = 3.0 * g_yy_0_xxxxyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyzzzz_1[i] * fz_be_0 + g_yyy_0_xxxxzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxzzzzz_0[i] = 3.0 * g_yy_0_xxxxzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxzzzzz_1[i] * fz_be_0 + g_yyy_0_xxxxzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxyyyyyy_0[i] = 3.0 * g_yy_0_xxxyyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyyyyy_1[i] * fz_be_0 + 6.0 * g_yyy_0_xxxyyyyy_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyyy_1[i] * wa_y[i];

        g_yyyy_0_xxxyyyyyz_0[i] = 3.0 * g_yy_0_xxxyyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyy_0_xxxyyyyz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyyz_1[i] * wa_y[i];

        g_yyyy_0_xxxyyyyzz_0[i] = 3.0 * g_yy_0_xxxyyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyy_0_xxxyyyzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyzz_1[i] * wa_y[i];

        g_yyyy_0_xxxyyyzzz_0[i] = 3.0 * g_yy_0_xxxyyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxxyyzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxyyzzzz_0[i] = 3.0 * g_yy_0_xxxyyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxyzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxyzzzzz_0[i] = 3.0 * g_yy_0_xxxyzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyzzzzz_1[i] * fz_be_0 + g_yyy_0_xxxzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxzzzzzz_0[i] = 3.0 * g_yy_0_xxxzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxzzzzzz_1[i] * fz_be_0 + g_yyy_0_xxxzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxyyyyyyy_0[i] = 3.0 * g_yy_0_xxyyyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyyyyy_1[i] * fz_be_0 + 7.0 * g_yyy_0_xxyyyyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyyy_1[i] * wa_y[i];

        g_yyyy_0_xxyyyyyyz_0[i] = 3.0 * g_yy_0_xxyyyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyy_0_xxyyyyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyyz_1[i] * wa_y[i];

        g_yyyy_0_xxyyyyyzz_0[i] = 3.0 * g_yy_0_xxyyyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyy_0_xxyyyyzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyzz_1[i] * wa_y[i];

        g_yyyy_0_xxyyyyzzz_0[i] = 3.0 * g_yy_0_xxyyyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyy_0_xxyyyzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyzzz_1[i] * wa_y[i];

        g_yyyy_0_xxyyyzzzz_0[i] = 3.0 * g_yy_0_xxyyyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxyyzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxyyzzzzz_0[i] = 3.0 * g_yy_0_xxyyzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxyzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxyzzzzzz_0[i] = 3.0 * g_yy_0_xxyzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyzzzzzz_1[i] * fz_be_0 + g_yyy_0_xxzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxzzzzzzz_0[i] = 3.0 * g_yy_0_xxzzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxzzzzzzz_1[i] * fz_be_0 + g_yyy_0_xxzzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xyyyyyyyy_0[i] = 3.0 * g_yy_0_xyyyyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyyyyy_1[i] * fz_be_0 + 8.0 * g_yyy_0_xyyyyyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyyy_1[i] * wa_y[i];

        g_yyyy_0_xyyyyyyyz_0[i] = 3.0 * g_yy_0_xyyyyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yyy_0_xyyyyyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyyz_1[i] * wa_y[i];

        g_yyyy_0_xyyyyyyzz_0[i] = 3.0 * g_yy_0_xyyyyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yyy_0_xyyyyyzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyzz_1[i] * wa_y[i];

        g_yyyy_0_xyyyyyzzz_0[i] = 3.0 * g_yy_0_xyyyyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yyy_0_xyyyyzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyzzz_1[i] * wa_y[i];

        g_yyyy_0_xyyyyzzzz_0[i] = 3.0 * g_yy_0_xyyyyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yyy_0_xyyyzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyzzzz_1[i] * wa_y[i];

        g_yyyy_0_xyyyzzzzz_0[i] = 3.0 * g_yy_0_xyyyzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xyyzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xyyzzzzzz_0[i] = 3.0 * g_yy_0_xyyzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xyzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xyzzzzzzz_0[i] = 3.0 * g_yy_0_xyzzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyzzzzzzz_1[i] * fz_be_0 + g_yyy_0_xzzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xzzzzzzzz_0[i] = 3.0 * g_yy_0_xzzzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xzzzzzzzz_1[i] * fz_be_0 + g_yyy_0_xzzzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_yyyyyyyyy_0[i] = 3.0 * g_yy_0_yyyyyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyyyyy_1[i] * fz_be_0 + 9.0 * g_yyy_0_yyyyyyyy_1[i] * fi_acd_0 + g_yyy_0_yyyyyyyyy_1[i] * wa_y[i];

        g_yyyy_0_yyyyyyyyz_0[i] = 3.0 * g_yy_0_yyyyyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyyyyz_1[i] * fz_be_0 + 8.0 * g_yyy_0_yyyyyyyz_1[i] * fi_acd_0 + g_yyy_0_yyyyyyyyz_1[i] * wa_y[i];

        g_yyyy_0_yyyyyyyzz_0[i] = 3.0 * g_yy_0_yyyyyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyyyzz_1[i] * fz_be_0 + 7.0 * g_yyy_0_yyyyyyzz_1[i] * fi_acd_0 + g_yyy_0_yyyyyyyzz_1[i] * wa_y[i];

        g_yyyy_0_yyyyyyzzz_0[i] = 3.0 * g_yy_0_yyyyyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyyzzz_1[i] * fz_be_0 + 6.0 * g_yyy_0_yyyyyzzz_1[i] * fi_acd_0 + g_yyy_0_yyyyyyzzz_1[i] * wa_y[i];

        g_yyyy_0_yyyyyzzzz_0[i] = 3.0 * g_yy_0_yyyyyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyzzzz_1[i] * fz_be_0 + 5.0 * g_yyy_0_yyyyzzzz_1[i] * fi_acd_0 + g_yyy_0_yyyyyzzzz_1[i] * wa_y[i];

        g_yyyy_0_yyyyzzzzz_0[i] = 3.0 * g_yy_0_yyyyzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyzzzzz_1[i] * fz_be_0 + 4.0 * g_yyy_0_yyyzzzzz_1[i] * fi_acd_0 + g_yyy_0_yyyyzzzzz_1[i] * wa_y[i];

        g_yyyy_0_yyyzzzzzz_0[i] = 3.0 * g_yy_0_yyyzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyzzzzzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_yyzzzzzz_1[i] * fi_acd_0 + g_yyy_0_yyyzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_yyzzzzzzz_0[i] = 3.0 * g_yy_0_yyzzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyzzzzzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_yzzzzzzz_1[i] * fi_acd_0 + g_yyy_0_yyzzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_yzzzzzzzz_0[i] = 3.0 * g_yy_0_yzzzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yzzzzzzzz_1[i] * fz_be_0 + g_yyy_0_zzzzzzzz_1[i] * fi_acd_0 + g_yyy_0_yzzzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_zzzzzzzzz_0[i] = 3.0 * g_yy_0_zzzzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_zzzzzzzzz_1[i] * fz_be_0 + g_yyy_0_zzzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 605-660 components of targeted buffer : GSM

    auto g_yyyz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 605);

    auto g_yyyz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 606);

    auto g_yyyz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 607);

    auto g_yyyz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 608);

    auto g_yyyz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 609);

    auto g_yyyz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 610);

    auto g_yyyz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 611);

    auto g_yyyz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 612);

    auto g_yyyz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 613);

    auto g_yyyz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 614);

    auto g_yyyz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 615);

    auto g_yyyz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 616);

    auto g_yyyz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 617);

    auto g_yyyz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 618);

    auto g_yyyz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 619);

    auto g_yyyz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 620);

    auto g_yyyz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 621);

    auto g_yyyz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 622);

    auto g_yyyz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 623);

    auto g_yyyz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 624);

    auto g_yyyz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 625);

    auto g_yyyz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 626);

    auto g_yyyz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 627);

    auto g_yyyz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 628);

    auto g_yyyz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 629);

    auto g_yyyz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 630);

    auto g_yyyz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 631);

    auto g_yyyz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 632);

    auto g_yyyz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 633);

    auto g_yyyz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 634);

    auto g_yyyz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 635);

    auto g_yyyz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 636);

    auto g_yyyz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 637);

    auto g_yyyz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 638);

    auto g_yyyz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 639);

    auto g_yyyz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 640);

    auto g_yyyz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 641);

    auto g_yyyz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 642);

    auto g_yyyz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 643);

    auto g_yyyz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 644);

    auto g_yyyz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 645);

    auto g_yyyz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 646);

    auto g_yyyz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 647);

    auto g_yyyz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 648);

    auto g_yyyz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 649);

    auto g_yyyz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 650);

    auto g_yyyz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 651);

    auto g_yyyz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 652);

    auto g_yyyz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 653);

    auto g_yyyz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 654);

    auto g_yyyz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 655);

    auto g_yyyz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 656);

    auto g_yyyz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 657);

    auto g_yyyz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 658);

    auto g_yyyz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 659);

    #pragma omp simd aligned(g_yyy_0_xxxxxxxx_1, g_yyy_0_xxxxxxxxx_1, g_yyy_0_xxxxxxxxy_1, g_yyy_0_xxxxxxxxz_1, g_yyy_0_xxxxxxxy_1, g_yyy_0_xxxxxxxyy_1, g_yyy_0_xxxxxxxyz_1, g_yyy_0_xxxxxxxz_1, g_yyy_0_xxxxxxxzz_1, g_yyy_0_xxxxxxyy_1, g_yyy_0_xxxxxxyyy_1, g_yyy_0_xxxxxxyyz_1, g_yyy_0_xxxxxxyz_1, g_yyy_0_xxxxxxyzz_1, g_yyy_0_xxxxxxzz_1, g_yyy_0_xxxxxxzzz_1, g_yyy_0_xxxxxyyy_1, g_yyy_0_xxxxxyyyy_1, g_yyy_0_xxxxxyyyz_1, g_yyy_0_xxxxxyyz_1, g_yyy_0_xxxxxyyzz_1, g_yyy_0_xxxxxyzz_1, g_yyy_0_xxxxxyzzz_1, g_yyy_0_xxxxxzzz_1, g_yyy_0_xxxxxzzzz_1, g_yyy_0_xxxxyyyy_1, g_yyy_0_xxxxyyyyy_1, g_yyy_0_xxxxyyyyz_1, g_yyy_0_xxxxyyyz_1, g_yyy_0_xxxxyyyzz_1, g_yyy_0_xxxxyyzz_1, g_yyy_0_xxxxyyzzz_1, g_yyy_0_xxxxyzzz_1, g_yyy_0_xxxxyzzzz_1, g_yyy_0_xxxxzzzz_1, g_yyy_0_xxxxzzzzz_1, g_yyy_0_xxxyyyyy_1, g_yyy_0_xxxyyyyyy_1, g_yyy_0_xxxyyyyyz_1, g_yyy_0_xxxyyyyz_1, g_yyy_0_xxxyyyyzz_1, g_yyy_0_xxxyyyzz_1, g_yyy_0_xxxyyyzzz_1, g_yyy_0_xxxyyzzz_1, g_yyy_0_xxxyyzzzz_1, g_yyy_0_xxxyzzzz_1, g_yyy_0_xxxyzzzzz_1, g_yyy_0_xxxzzzzz_1, g_yyy_0_xxxzzzzzz_1, g_yyy_0_xxyyyyyy_1, g_yyy_0_xxyyyyyyy_1, g_yyy_0_xxyyyyyyz_1, g_yyy_0_xxyyyyyz_1, g_yyy_0_xxyyyyyzz_1, g_yyy_0_xxyyyyzz_1, g_yyy_0_xxyyyyzzz_1, g_yyy_0_xxyyyzzz_1, g_yyy_0_xxyyyzzzz_1, g_yyy_0_xxyyzzzz_1, g_yyy_0_xxyyzzzzz_1, g_yyy_0_xxyzzzzz_1, g_yyy_0_xxyzzzzzz_1, g_yyy_0_xxzzzzzz_1, g_yyy_0_xxzzzzzzz_1, g_yyy_0_xyyyyyyy_1, g_yyy_0_xyyyyyyyy_1, g_yyy_0_xyyyyyyyz_1, g_yyy_0_xyyyyyyz_1, g_yyy_0_xyyyyyyzz_1, g_yyy_0_xyyyyyzz_1, g_yyy_0_xyyyyyzzz_1, g_yyy_0_xyyyyzzz_1, g_yyy_0_xyyyyzzzz_1, g_yyy_0_xyyyzzzz_1, g_yyy_0_xyyyzzzzz_1, g_yyy_0_xyyzzzzz_1, g_yyy_0_xyyzzzzzz_1, g_yyy_0_xyzzzzzz_1, g_yyy_0_xyzzzzzzz_1, g_yyy_0_xzzzzzzz_1, g_yyy_0_xzzzzzzzz_1, g_yyy_0_yyyyyyyy_1, g_yyy_0_yyyyyyyyy_1, g_yyy_0_yyyyyyyyz_1, g_yyy_0_yyyyyyyz_1, g_yyy_0_yyyyyyyzz_1, g_yyy_0_yyyyyyzz_1, g_yyy_0_yyyyyyzzz_1, g_yyy_0_yyyyyzzz_1, g_yyy_0_yyyyyzzzz_1, g_yyy_0_yyyyzzzz_1, g_yyy_0_yyyyzzzzz_1, g_yyy_0_yyyzzzzz_1, g_yyy_0_yyyzzzzzz_1, g_yyy_0_yyzzzzzz_1, g_yyy_0_yyzzzzzzz_1, g_yyy_0_yzzzzzzz_1, g_yyy_0_yzzzzzzzz_1, g_yyy_0_zzzzzzzz_1, g_yyy_0_zzzzzzzzz_1, g_yyyz_0_xxxxxxxxx_0, g_yyyz_0_xxxxxxxxy_0, g_yyyz_0_xxxxxxxxz_0, g_yyyz_0_xxxxxxxyy_0, g_yyyz_0_xxxxxxxyz_0, g_yyyz_0_xxxxxxxzz_0, g_yyyz_0_xxxxxxyyy_0, g_yyyz_0_xxxxxxyyz_0, g_yyyz_0_xxxxxxyzz_0, g_yyyz_0_xxxxxxzzz_0, g_yyyz_0_xxxxxyyyy_0, g_yyyz_0_xxxxxyyyz_0, g_yyyz_0_xxxxxyyzz_0, g_yyyz_0_xxxxxyzzz_0, g_yyyz_0_xxxxxzzzz_0, g_yyyz_0_xxxxyyyyy_0, g_yyyz_0_xxxxyyyyz_0, g_yyyz_0_xxxxyyyzz_0, g_yyyz_0_xxxxyyzzz_0, g_yyyz_0_xxxxyzzzz_0, g_yyyz_0_xxxxzzzzz_0, g_yyyz_0_xxxyyyyyy_0, g_yyyz_0_xxxyyyyyz_0, g_yyyz_0_xxxyyyyzz_0, g_yyyz_0_xxxyyyzzz_0, g_yyyz_0_xxxyyzzzz_0, g_yyyz_0_xxxyzzzzz_0, g_yyyz_0_xxxzzzzzz_0, g_yyyz_0_xxyyyyyyy_0, g_yyyz_0_xxyyyyyyz_0, g_yyyz_0_xxyyyyyzz_0, g_yyyz_0_xxyyyyzzz_0, g_yyyz_0_xxyyyzzzz_0, g_yyyz_0_xxyyzzzzz_0, g_yyyz_0_xxyzzzzzz_0, g_yyyz_0_xxzzzzzzz_0, g_yyyz_0_xyyyyyyyy_0, g_yyyz_0_xyyyyyyyz_0, g_yyyz_0_xyyyyyyzz_0, g_yyyz_0_xyyyyyzzz_0, g_yyyz_0_xyyyyzzzz_0, g_yyyz_0_xyyyzzzzz_0, g_yyyz_0_xyyzzzzzz_0, g_yyyz_0_xyzzzzzzz_0, g_yyyz_0_xzzzzzzzz_0, g_yyyz_0_yyyyyyyyy_0, g_yyyz_0_yyyyyyyyz_0, g_yyyz_0_yyyyyyyzz_0, g_yyyz_0_yyyyyyzzz_0, g_yyyz_0_yyyyyzzzz_0, g_yyyz_0_yyyyzzzzz_0, g_yyyz_0_yyyzzzzzz_0, g_yyyz_0_yyzzzzzzz_0, g_yyyz_0_yzzzzzzzz_0, g_yyyz_0_zzzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyz_0_xxxxxxxxx_0[i] = g_yyy_0_xxxxxxxxx_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxxxy_0[i] = g_yyy_0_xxxxxxxxy_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxxxz_0[i] = g_yyy_0_xxxxxxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxxz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxxyy_0[i] = g_yyy_0_xxxxxxxyy_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxxyz_0[i] = g_yyy_0_xxxxxxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxyz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxxzz_0[i] = 2.0 * g_yyy_0_xxxxxxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxyyy_0[i] = g_yyy_0_xxxxxxyyy_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxyyz_0[i] = g_yyy_0_xxxxxxyy_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyyz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxyzz_0[i] = 2.0 * g_yyy_0_xxxxxxyz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxzzz_0[i] = 3.0 * g_yyy_0_xxxxxxzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxyyyy_0[i] = g_yyy_0_xxxxxyyyy_1[i] * wa_z[i];

        g_yyyz_0_xxxxxyyyz_0[i] = g_yyy_0_xxxxxyyy_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyyz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxyyzz_0[i] = 2.0 * g_yyy_0_xxxxxyyz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxyzzz_0[i] = 3.0 * g_yyy_0_xxxxxyzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxzzzz_0[i] = 4.0 * g_yyy_0_xxxxxzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxyyyyy_0[i] = g_yyy_0_xxxxyyyyy_1[i] * wa_z[i];

        g_yyyz_0_xxxxyyyyz_0[i] = g_yyy_0_xxxxyyyy_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyyz_1[i] * wa_z[i];

        g_yyyz_0_xxxxyyyzz_0[i] = 2.0 * g_yyy_0_xxxxyyyz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxyyzzz_0[i] = 3.0 * g_yyy_0_xxxxyyzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxyzzzz_0[i] = 4.0 * g_yyy_0_xxxxyzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxzzzzz_0[i] = 5.0 * g_yyy_0_xxxxzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxyyyyyy_0[i] = g_yyy_0_xxxyyyyyy_1[i] * wa_z[i];

        g_yyyz_0_xxxyyyyyz_0[i] = g_yyy_0_xxxyyyyy_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyyz_1[i] * wa_z[i];

        g_yyyz_0_xxxyyyyzz_0[i] = 2.0 * g_yyy_0_xxxyyyyz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyzz_1[i] * wa_z[i];

        g_yyyz_0_xxxyyyzzz_0[i] = 3.0 * g_yyy_0_xxxyyyzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxyyzzzz_0[i] = 4.0 * g_yyy_0_xxxyyzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxyzzzzz_0[i] = 5.0 * g_yyy_0_xxxyzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxzzzzzz_0[i] = 6.0 * g_yyy_0_xxxzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxyyyyyyy_0[i] = g_yyy_0_xxyyyyyyy_1[i] * wa_z[i];

        g_yyyz_0_xxyyyyyyz_0[i] = g_yyy_0_xxyyyyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyyz_1[i] * wa_z[i];

        g_yyyz_0_xxyyyyyzz_0[i] = 2.0 * g_yyy_0_xxyyyyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyzz_1[i] * wa_z[i];

        g_yyyz_0_xxyyyyzzz_0[i] = 3.0 * g_yyy_0_xxyyyyzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyzzz_1[i] * wa_z[i];

        g_yyyz_0_xxyyyzzzz_0[i] = 4.0 * g_yyy_0_xxyyyzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxyyzzzzz_0[i] = 5.0 * g_yyy_0_xxyyzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxyzzzzzz_0[i] = 6.0 * g_yyy_0_xxyzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxzzzzzzz_0[i] = 7.0 * g_yyy_0_xxzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxzzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xyyyyyyyy_0[i] = g_yyy_0_xyyyyyyyy_1[i] * wa_z[i];

        g_yyyz_0_xyyyyyyyz_0[i] = g_yyy_0_xyyyyyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyyz_1[i] * wa_z[i];

        g_yyyz_0_xyyyyyyzz_0[i] = 2.0 * g_yyy_0_xyyyyyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyzz_1[i] * wa_z[i];

        g_yyyz_0_xyyyyyzzz_0[i] = 3.0 * g_yyy_0_xyyyyyzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyzzz_1[i] * wa_z[i];

        g_yyyz_0_xyyyyzzzz_0[i] = 4.0 * g_yyy_0_xyyyyzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyzzzz_1[i] * wa_z[i];

        g_yyyz_0_xyyyzzzzz_0[i] = 5.0 * g_yyy_0_xyyyzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xyyzzzzzz_0[i] = 6.0 * g_yyy_0_xyyzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xyzzzzzzz_0[i] = 7.0 * g_yyy_0_xyzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xzzzzzzzz_0[i] = 8.0 * g_yyy_0_xzzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xzzzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_yyyyyyyyy_0[i] = g_yyy_0_yyyyyyyyy_1[i] * wa_z[i];

        g_yyyz_0_yyyyyyyyz_0[i] = g_yyy_0_yyyyyyyy_1[i] * fi_acd_0 + g_yyy_0_yyyyyyyyz_1[i] * wa_z[i];

        g_yyyz_0_yyyyyyyzz_0[i] = 2.0 * g_yyy_0_yyyyyyyz_1[i] * fi_acd_0 + g_yyy_0_yyyyyyyzz_1[i] * wa_z[i];

        g_yyyz_0_yyyyyyzzz_0[i] = 3.0 * g_yyy_0_yyyyyyzz_1[i] * fi_acd_0 + g_yyy_0_yyyyyyzzz_1[i] * wa_z[i];

        g_yyyz_0_yyyyyzzzz_0[i] = 4.0 * g_yyy_0_yyyyyzzz_1[i] * fi_acd_0 + g_yyy_0_yyyyyzzzz_1[i] * wa_z[i];

        g_yyyz_0_yyyyzzzzz_0[i] = 5.0 * g_yyy_0_yyyyzzzz_1[i] * fi_acd_0 + g_yyy_0_yyyyzzzzz_1[i] * wa_z[i];

        g_yyyz_0_yyyzzzzzz_0[i] = 6.0 * g_yyy_0_yyyzzzzz_1[i] * fi_acd_0 + g_yyy_0_yyyzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_yyzzzzzzz_0[i] = 7.0 * g_yyy_0_yyzzzzzz_1[i] * fi_acd_0 + g_yyy_0_yyzzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_yzzzzzzzz_0[i] = 8.0 * g_yyy_0_yzzzzzzz_1[i] * fi_acd_0 + g_yyy_0_yzzzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_zzzzzzzzz_0[i] = 9.0 * g_yyy_0_zzzzzzzz_1[i] * fi_acd_0 + g_yyy_0_zzzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 660-715 components of targeted buffer : GSM

    auto g_yyzz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 660);

    auto g_yyzz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 661);

    auto g_yyzz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 662);

    auto g_yyzz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 663);

    auto g_yyzz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 664);

    auto g_yyzz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 665);

    auto g_yyzz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 666);

    auto g_yyzz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 667);

    auto g_yyzz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 668);

    auto g_yyzz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 669);

    auto g_yyzz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 670);

    auto g_yyzz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 671);

    auto g_yyzz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 672);

    auto g_yyzz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 673);

    auto g_yyzz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 674);

    auto g_yyzz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 675);

    auto g_yyzz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 676);

    auto g_yyzz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 677);

    auto g_yyzz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 678);

    auto g_yyzz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 679);

    auto g_yyzz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 680);

    auto g_yyzz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 681);

    auto g_yyzz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 682);

    auto g_yyzz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 683);

    auto g_yyzz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 684);

    auto g_yyzz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 685);

    auto g_yyzz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 686);

    auto g_yyzz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 687);

    auto g_yyzz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 688);

    auto g_yyzz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 689);

    auto g_yyzz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 690);

    auto g_yyzz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 691);

    auto g_yyzz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 692);

    auto g_yyzz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 693);

    auto g_yyzz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 694);

    auto g_yyzz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 695);

    auto g_yyzz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 696);

    auto g_yyzz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 697);

    auto g_yyzz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 698);

    auto g_yyzz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 699);

    auto g_yyzz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 700);

    auto g_yyzz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 701);

    auto g_yyzz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 702);

    auto g_yyzz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 703);

    auto g_yyzz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 704);

    auto g_yyzz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 705);

    auto g_yyzz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 706);

    auto g_yyzz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 707);

    auto g_yyzz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 708);

    auto g_yyzz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 709);

    auto g_yyzz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 710);

    auto g_yyzz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 711);

    auto g_yyzz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 712);

    auto g_yyzz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 713);

    auto g_yyzz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 714);

    #pragma omp simd aligned(g_yy_0_xxxxxxxxy_0, g_yy_0_xxxxxxxxy_1, g_yy_0_xxxxxxxyy_0, g_yy_0_xxxxxxxyy_1, g_yy_0_xxxxxxyyy_0, g_yy_0_xxxxxxyyy_1, g_yy_0_xxxxxyyyy_0, g_yy_0_xxxxxyyyy_1, g_yy_0_xxxxyyyyy_0, g_yy_0_xxxxyyyyy_1, g_yy_0_xxxyyyyyy_0, g_yy_0_xxxyyyyyy_1, g_yy_0_xxyyyyyyy_0, g_yy_0_xxyyyyyyy_1, g_yy_0_xyyyyyyyy_0, g_yy_0_xyyyyyyyy_1, g_yy_0_yyyyyyyyy_0, g_yy_0_yyyyyyyyy_1, g_yyz_0_xxxxxxxxy_1, g_yyz_0_xxxxxxxyy_1, g_yyz_0_xxxxxxyyy_1, g_yyz_0_xxxxxyyyy_1, g_yyz_0_xxxxyyyyy_1, g_yyz_0_xxxyyyyyy_1, g_yyz_0_xxyyyyyyy_1, g_yyz_0_xyyyyyyyy_1, g_yyz_0_yyyyyyyyy_1, g_yyzz_0_xxxxxxxxx_0, g_yyzz_0_xxxxxxxxy_0, g_yyzz_0_xxxxxxxxz_0, g_yyzz_0_xxxxxxxyy_0, g_yyzz_0_xxxxxxxyz_0, g_yyzz_0_xxxxxxxzz_0, g_yyzz_0_xxxxxxyyy_0, g_yyzz_0_xxxxxxyyz_0, g_yyzz_0_xxxxxxyzz_0, g_yyzz_0_xxxxxxzzz_0, g_yyzz_0_xxxxxyyyy_0, g_yyzz_0_xxxxxyyyz_0, g_yyzz_0_xxxxxyyzz_0, g_yyzz_0_xxxxxyzzz_0, g_yyzz_0_xxxxxzzzz_0, g_yyzz_0_xxxxyyyyy_0, g_yyzz_0_xxxxyyyyz_0, g_yyzz_0_xxxxyyyzz_0, g_yyzz_0_xxxxyyzzz_0, g_yyzz_0_xxxxyzzzz_0, g_yyzz_0_xxxxzzzzz_0, g_yyzz_0_xxxyyyyyy_0, g_yyzz_0_xxxyyyyyz_0, g_yyzz_0_xxxyyyyzz_0, g_yyzz_0_xxxyyyzzz_0, g_yyzz_0_xxxyyzzzz_0, g_yyzz_0_xxxyzzzzz_0, g_yyzz_0_xxxzzzzzz_0, g_yyzz_0_xxyyyyyyy_0, g_yyzz_0_xxyyyyyyz_0, g_yyzz_0_xxyyyyyzz_0, g_yyzz_0_xxyyyyzzz_0, g_yyzz_0_xxyyyzzzz_0, g_yyzz_0_xxyyzzzzz_0, g_yyzz_0_xxyzzzzzz_0, g_yyzz_0_xxzzzzzzz_0, g_yyzz_0_xyyyyyyyy_0, g_yyzz_0_xyyyyyyyz_0, g_yyzz_0_xyyyyyyzz_0, g_yyzz_0_xyyyyyzzz_0, g_yyzz_0_xyyyyzzzz_0, g_yyzz_0_xyyyzzzzz_0, g_yyzz_0_xyyzzzzzz_0, g_yyzz_0_xyzzzzzzz_0, g_yyzz_0_xzzzzzzzz_0, g_yyzz_0_yyyyyyyyy_0, g_yyzz_0_yyyyyyyyz_0, g_yyzz_0_yyyyyyyzz_0, g_yyzz_0_yyyyyyzzz_0, g_yyzz_0_yyyyyzzzz_0, g_yyzz_0_yyyyzzzzz_0, g_yyzz_0_yyyzzzzzz_0, g_yyzz_0_yyzzzzzzz_0, g_yyzz_0_yzzzzzzzz_0, g_yyzz_0_zzzzzzzzz_0, g_yzz_0_xxxxxxxxx_1, g_yzz_0_xxxxxxxxz_1, g_yzz_0_xxxxxxxyz_1, g_yzz_0_xxxxxxxz_1, g_yzz_0_xxxxxxxzz_1, g_yzz_0_xxxxxxyyz_1, g_yzz_0_xxxxxxyz_1, g_yzz_0_xxxxxxyzz_1, g_yzz_0_xxxxxxzz_1, g_yzz_0_xxxxxxzzz_1, g_yzz_0_xxxxxyyyz_1, g_yzz_0_xxxxxyyz_1, g_yzz_0_xxxxxyyzz_1, g_yzz_0_xxxxxyzz_1, g_yzz_0_xxxxxyzzz_1, g_yzz_0_xxxxxzzz_1, g_yzz_0_xxxxxzzzz_1, g_yzz_0_xxxxyyyyz_1, g_yzz_0_xxxxyyyz_1, g_yzz_0_xxxxyyyzz_1, g_yzz_0_xxxxyyzz_1, g_yzz_0_xxxxyyzzz_1, g_yzz_0_xxxxyzzz_1, g_yzz_0_xxxxyzzzz_1, g_yzz_0_xxxxzzzz_1, g_yzz_0_xxxxzzzzz_1, g_yzz_0_xxxyyyyyz_1, g_yzz_0_xxxyyyyz_1, g_yzz_0_xxxyyyyzz_1, g_yzz_0_xxxyyyzz_1, g_yzz_0_xxxyyyzzz_1, g_yzz_0_xxxyyzzz_1, g_yzz_0_xxxyyzzzz_1, g_yzz_0_xxxyzzzz_1, g_yzz_0_xxxyzzzzz_1, g_yzz_0_xxxzzzzz_1, g_yzz_0_xxxzzzzzz_1, g_yzz_0_xxyyyyyyz_1, g_yzz_0_xxyyyyyz_1, g_yzz_0_xxyyyyyzz_1, g_yzz_0_xxyyyyzz_1, g_yzz_0_xxyyyyzzz_1, g_yzz_0_xxyyyzzz_1, g_yzz_0_xxyyyzzzz_1, g_yzz_0_xxyyzzzz_1, g_yzz_0_xxyyzzzzz_1, g_yzz_0_xxyzzzzz_1, g_yzz_0_xxyzzzzzz_1, g_yzz_0_xxzzzzzz_1, g_yzz_0_xxzzzzzzz_1, g_yzz_0_xyyyyyyyz_1, g_yzz_0_xyyyyyyz_1, g_yzz_0_xyyyyyyzz_1, g_yzz_0_xyyyyyzz_1, g_yzz_0_xyyyyyzzz_1, g_yzz_0_xyyyyzzz_1, g_yzz_0_xyyyyzzzz_1, g_yzz_0_xyyyzzzz_1, g_yzz_0_xyyyzzzzz_1, g_yzz_0_xyyzzzzz_1, g_yzz_0_xyyzzzzzz_1, g_yzz_0_xyzzzzzz_1, g_yzz_0_xyzzzzzzz_1, g_yzz_0_xzzzzzzz_1, g_yzz_0_xzzzzzzzz_1, g_yzz_0_yyyyyyyyz_1, g_yzz_0_yyyyyyyz_1, g_yzz_0_yyyyyyyzz_1, g_yzz_0_yyyyyyzz_1, g_yzz_0_yyyyyyzzz_1, g_yzz_0_yyyyyzzz_1, g_yzz_0_yyyyyzzzz_1, g_yzz_0_yyyyzzzz_1, g_yzz_0_yyyyzzzzz_1, g_yzz_0_yyyzzzzz_1, g_yzz_0_yyyzzzzzz_1, g_yzz_0_yyzzzzzz_1, g_yzz_0_yyzzzzzzz_1, g_yzz_0_yzzzzzzz_1, g_yzz_0_yzzzzzzzz_1, g_yzz_0_zzzzzzzz_1, g_yzz_0_zzzzzzzzz_1, g_zz_0_xxxxxxxxx_0, g_zz_0_xxxxxxxxx_1, g_zz_0_xxxxxxxxz_0, g_zz_0_xxxxxxxxz_1, g_zz_0_xxxxxxxyz_0, g_zz_0_xxxxxxxyz_1, g_zz_0_xxxxxxxzz_0, g_zz_0_xxxxxxxzz_1, g_zz_0_xxxxxxyyz_0, g_zz_0_xxxxxxyyz_1, g_zz_0_xxxxxxyzz_0, g_zz_0_xxxxxxyzz_1, g_zz_0_xxxxxxzzz_0, g_zz_0_xxxxxxzzz_1, g_zz_0_xxxxxyyyz_0, g_zz_0_xxxxxyyyz_1, g_zz_0_xxxxxyyzz_0, g_zz_0_xxxxxyyzz_1, g_zz_0_xxxxxyzzz_0, g_zz_0_xxxxxyzzz_1, g_zz_0_xxxxxzzzz_0, g_zz_0_xxxxxzzzz_1, g_zz_0_xxxxyyyyz_0, g_zz_0_xxxxyyyyz_1, g_zz_0_xxxxyyyzz_0, g_zz_0_xxxxyyyzz_1, g_zz_0_xxxxyyzzz_0, g_zz_0_xxxxyyzzz_1, g_zz_0_xxxxyzzzz_0, g_zz_0_xxxxyzzzz_1, g_zz_0_xxxxzzzzz_0, g_zz_0_xxxxzzzzz_1, g_zz_0_xxxyyyyyz_0, g_zz_0_xxxyyyyyz_1, g_zz_0_xxxyyyyzz_0, g_zz_0_xxxyyyyzz_1, g_zz_0_xxxyyyzzz_0, g_zz_0_xxxyyyzzz_1, g_zz_0_xxxyyzzzz_0, g_zz_0_xxxyyzzzz_1, g_zz_0_xxxyzzzzz_0, g_zz_0_xxxyzzzzz_1, g_zz_0_xxxzzzzzz_0, g_zz_0_xxxzzzzzz_1, g_zz_0_xxyyyyyyz_0, g_zz_0_xxyyyyyyz_1, g_zz_0_xxyyyyyzz_0, g_zz_0_xxyyyyyzz_1, g_zz_0_xxyyyyzzz_0, g_zz_0_xxyyyyzzz_1, g_zz_0_xxyyyzzzz_0, g_zz_0_xxyyyzzzz_1, g_zz_0_xxyyzzzzz_0, g_zz_0_xxyyzzzzz_1, g_zz_0_xxyzzzzzz_0, g_zz_0_xxyzzzzzz_1, g_zz_0_xxzzzzzzz_0, g_zz_0_xxzzzzzzz_1, g_zz_0_xyyyyyyyz_0, g_zz_0_xyyyyyyyz_1, g_zz_0_xyyyyyyzz_0, g_zz_0_xyyyyyyzz_1, g_zz_0_xyyyyyzzz_0, g_zz_0_xyyyyyzzz_1, g_zz_0_xyyyyzzzz_0, g_zz_0_xyyyyzzzz_1, g_zz_0_xyyyzzzzz_0, g_zz_0_xyyyzzzzz_1, g_zz_0_xyyzzzzzz_0, g_zz_0_xyyzzzzzz_1, g_zz_0_xyzzzzzzz_0, g_zz_0_xyzzzzzzz_1, g_zz_0_xzzzzzzzz_0, g_zz_0_xzzzzzzzz_1, g_zz_0_yyyyyyyyz_0, g_zz_0_yyyyyyyyz_1, g_zz_0_yyyyyyyzz_0, g_zz_0_yyyyyyyzz_1, g_zz_0_yyyyyyzzz_0, g_zz_0_yyyyyyzzz_1, g_zz_0_yyyyyzzzz_0, g_zz_0_yyyyyzzzz_1, g_zz_0_yyyyzzzzz_0, g_zz_0_yyyyzzzzz_1, g_zz_0_yyyzzzzzz_0, g_zz_0_yyyzzzzzz_1, g_zz_0_yyzzzzzzz_0, g_zz_0_yyzzzzzzz_1, g_zz_0_yzzzzzzzz_0, g_zz_0_yzzzzzzzz_1, g_zz_0_zzzzzzzzz_0, g_zz_0_zzzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzz_0_xxxxxxxxx_0[i] = g_zz_0_xxxxxxxxx_0[i] * fbe_0 - g_zz_0_xxxxxxxxx_1[i] * fz_be_0 + g_yzz_0_xxxxxxxxx_1[i] * wa_y[i];

        g_yyzz_0_xxxxxxxxy_0[i] = g_yy_0_xxxxxxxxy_0[i] * fbe_0 - g_yy_0_xxxxxxxxy_1[i] * fz_be_0 + g_yyz_0_xxxxxxxxy_1[i] * wa_z[i];

        g_yyzz_0_xxxxxxxxz_0[i] = g_zz_0_xxxxxxxxz_0[i] * fbe_0 - g_zz_0_xxxxxxxxz_1[i] * fz_be_0 + g_yzz_0_xxxxxxxxz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxxxyy_0[i] = g_yy_0_xxxxxxxyy_0[i] * fbe_0 - g_yy_0_xxxxxxxyy_1[i] * fz_be_0 + g_yyz_0_xxxxxxxyy_1[i] * wa_z[i];

        g_yyzz_0_xxxxxxxyz_0[i] = g_zz_0_xxxxxxxyz_0[i] * fbe_0 - g_zz_0_xxxxxxxyz_1[i] * fz_be_0 + g_yzz_0_xxxxxxxz_1[i] * fi_acd_0 + g_yzz_0_xxxxxxxyz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxxxzz_0[i] = g_zz_0_xxxxxxxzz_0[i] * fbe_0 - g_zz_0_xxxxxxxzz_1[i] * fz_be_0 + g_yzz_0_xxxxxxxzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxxyyy_0[i] = g_yy_0_xxxxxxyyy_0[i] * fbe_0 - g_yy_0_xxxxxxyyy_1[i] * fz_be_0 + g_yyz_0_xxxxxxyyy_1[i] * wa_z[i];

        g_yyzz_0_xxxxxxyyz_0[i] = g_zz_0_xxxxxxyyz_0[i] * fbe_0 - g_zz_0_xxxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxxxxxyz_1[i] * fi_acd_0 + g_yzz_0_xxxxxxyyz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxxyzz_0[i] = g_zz_0_xxxxxxyzz_0[i] * fbe_0 - g_zz_0_xxxxxxyzz_1[i] * fz_be_0 + g_yzz_0_xxxxxxzz_1[i] * fi_acd_0 + g_yzz_0_xxxxxxyzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxxzzz_0[i] = g_zz_0_xxxxxxzzz_0[i] * fbe_0 - g_zz_0_xxxxxxzzz_1[i] * fz_be_0 + g_yzz_0_xxxxxxzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxyyyy_0[i] = g_yy_0_xxxxxyyyy_0[i] * fbe_0 - g_yy_0_xxxxxyyyy_1[i] * fz_be_0 + g_yyz_0_xxxxxyyyy_1[i] * wa_z[i];

        g_yyzz_0_xxxxxyyyz_0[i] = g_zz_0_xxxxxyyyz_0[i] * fbe_0 - g_zz_0_xxxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xxxxxyyz_1[i] * fi_acd_0 + g_yzz_0_xxxxxyyyz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxyyzz_0[i] = g_zz_0_xxxxxyyzz_0[i] * fbe_0 - g_zz_0_xxxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxxxxyzz_1[i] * fi_acd_0 + g_yzz_0_xxxxxyyzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxyzzz_0[i] = g_zz_0_xxxxxyzzz_0[i] * fbe_0 - g_zz_0_xxxxxyzzz_1[i] * fz_be_0 + g_yzz_0_xxxxxzzz_1[i] * fi_acd_0 + g_yzz_0_xxxxxyzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxzzzz_0[i] = g_zz_0_xxxxxzzzz_0[i] * fbe_0 - g_zz_0_xxxxxzzzz_1[i] * fz_be_0 + g_yzz_0_xxxxxzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxyyyyy_0[i] = g_yy_0_xxxxyyyyy_0[i] * fbe_0 - g_yy_0_xxxxyyyyy_1[i] * fz_be_0 + g_yyz_0_xxxxyyyyy_1[i] * wa_z[i];

        g_yyzz_0_xxxxyyyyz_0[i] = g_zz_0_xxxxyyyyz_0[i] * fbe_0 - g_zz_0_xxxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yzz_0_xxxxyyyz_1[i] * fi_acd_0 + g_yzz_0_xxxxyyyyz_1[i] * wa_y[i];

        g_yyzz_0_xxxxyyyzz_0[i] = g_zz_0_xxxxyyyzz_0[i] * fbe_0 - g_zz_0_xxxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xxxxyyzz_1[i] * fi_acd_0 + g_yzz_0_xxxxyyyzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxyyzzz_0[i] = g_zz_0_xxxxyyzzz_0[i] * fbe_0 - g_zz_0_xxxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxxxyzzz_1[i] * fi_acd_0 + g_yzz_0_xxxxyyzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxyzzzz_0[i] = g_zz_0_xxxxyzzzz_0[i] * fbe_0 - g_zz_0_xxxxyzzzz_1[i] * fz_be_0 + g_yzz_0_xxxxzzzz_1[i] * fi_acd_0 + g_yzz_0_xxxxyzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxzzzzz_0[i] = g_zz_0_xxxxzzzzz_0[i] * fbe_0 - g_zz_0_xxxxzzzzz_1[i] * fz_be_0 + g_yzz_0_xxxxzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxyyyyyy_0[i] = g_yy_0_xxxyyyyyy_0[i] * fbe_0 - g_yy_0_xxxyyyyyy_1[i] * fz_be_0 + g_yyz_0_xxxyyyyyy_1[i] * wa_z[i];

        g_yyzz_0_xxxyyyyyz_0[i] = g_zz_0_xxxyyyyyz_0[i] * fbe_0 - g_zz_0_xxxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yzz_0_xxxyyyyz_1[i] * fi_acd_0 + g_yzz_0_xxxyyyyyz_1[i] * wa_y[i];

        g_yyzz_0_xxxyyyyzz_0[i] = g_zz_0_xxxyyyyzz_0[i] * fbe_0 - g_zz_0_xxxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yzz_0_xxxyyyzz_1[i] * fi_acd_0 + g_yzz_0_xxxyyyyzz_1[i] * wa_y[i];

        g_yyzz_0_xxxyyyzzz_0[i] = g_zz_0_xxxyyyzzz_0[i] * fbe_0 - g_zz_0_xxxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xxxyyzzz_1[i] * fi_acd_0 + g_yzz_0_xxxyyyzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxyyzzzz_0[i] = g_zz_0_xxxyyzzzz_0[i] * fbe_0 - g_zz_0_xxxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxxyzzzz_1[i] * fi_acd_0 + g_yzz_0_xxxyyzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxyzzzzz_0[i] = g_zz_0_xxxyzzzzz_0[i] * fbe_0 - g_zz_0_xxxyzzzzz_1[i] * fz_be_0 + g_yzz_0_xxxzzzzz_1[i] * fi_acd_0 + g_yzz_0_xxxyzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxzzzzzz_0[i] = g_zz_0_xxxzzzzzz_0[i] * fbe_0 - g_zz_0_xxxzzzzzz_1[i] * fz_be_0 + g_yzz_0_xxxzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxyyyyyyy_0[i] = g_yy_0_xxyyyyyyy_0[i] * fbe_0 - g_yy_0_xxyyyyyyy_1[i] * fz_be_0 + g_yyz_0_xxyyyyyyy_1[i] * wa_z[i];

        g_yyzz_0_xxyyyyyyz_0[i] = g_zz_0_xxyyyyyyz_0[i] * fbe_0 - g_zz_0_xxyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yzz_0_xxyyyyyz_1[i] * fi_acd_0 + g_yzz_0_xxyyyyyyz_1[i] * wa_y[i];

        g_yyzz_0_xxyyyyyzz_0[i] = g_zz_0_xxyyyyyzz_0[i] * fbe_0 - g_zz_0_xxyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yzz_0_xxyyyyzz_1[i] * fi_acd_0 + g_yzz_0_xxyyyyyzz_1[i] * wa_y[i];

        g_yyzz_0_xxyyyyzzz_0[i] = g_zz_0_xxyyyyzzz_0[i] * fbe_0 - g_zz_0_xxyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yzz_0_xxyyyzzz_1[i] * fi_acd_0 + g_yzz_0_xxyyyyzzz_1[i] * wa_y[i];

        g_yyzz_0_xxyyyzzzz_0[i] = g_zz_0_xxyyyzzzz_0[i] * fbe_0 - g_zz_0_xxyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xxyyzzzz_1[i] * fi_acd_0 + g_yzz_0_xxyyyzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxyyzzzzz_0[i] = g_zz_0_xxyyzzzzz_0[i] * fbe_0 - g_zz_0_xxyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxyzzzzz_1[i] * fi_acd_0 + g_yzz_0_xxyyzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxyzzzzzz_0[i] = g_zz_0_xxyzzzzzz_0[i] * fbe_0 - g_zz_0_xxyzzzzzz_1[i] * fz_be_0 + g_yzz_0_xxzzzzzz_1[i] * fi_acd_0 + g_yzz_0_xxyzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxzzzzzzz_0[i] = g_zz_0_xxzzzzzzz_0[i] * fbe_0 - g_zz_0_xxzzzzzzz_1[i] * fz_be_0 + g_yzz_0_xxzzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xyyyyyyyy_0[i] = g_yy_0_xyyyyyyyy_0[i] * fbe_0 - g_yy_0_xyyyyyyyy_1[i] * fz_be_0 + g_yyz_0_xyyyyyyyy_1[i] * wa_z[i];

        g_yyzz_0_xyyyyyyyz_0[i] = g_zz_0_xyyyyyyyz_0[i] * fbe_0 - g_zz_0_xyyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yzz_0_xyyyyyyz_1[i] * fi_acd_0 + g_yzz_0_xyyyyyyyz_1[i] * wa_y[i];

        g_yyzz_0_xyyyyyyzz_0[i] = g_zz_0_xyyyyyyzz_0[i] * fbe_0 - g_zz_0_xyyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yzz_0_xyyyyyzz_1[i] * fi_acd_0 + g_yzz_0_xyyyyyyzz_1[i] * wa_y[i];

        g_yyzz_0_xyyyyyzzz_0[i] = g_zz_0_xyyyyyzzz_0[i] * fbe_0 - g_zz_0_xyyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yzz_0_xyyyyzzz_1[i] * fi_acd_0 + g_yzz_0_xyyyyyzzz_1[i] * wa_y[i];

        g_yyzz_0_xyyyyzzzz_0[i] = g_zz_0_xyyyyzzzz_0[i] * fbe_0 - g_zz_0_xyyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yzz_0_xyyyzzzz_1[i] * fi_acd_0 + g_yzz_0_xyyyyzzzz_1[i] * wa_y[i];

        g_yyzz_0_xyyyzzzzz_0[i] = g_zz_0_xyyyzzzzz_0[i] * fbe_0 - g_zz_0_xyyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xyyzzzzz_1[i] * fi_acd_0 + g_yzz_0_xyyyzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xyyzzzzzz_0[i] = g_zz_0_xyyzzzzzz_0[i] * fbe_0 - g_zz_0_xyyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xyzzzzzz_1[i] * fi_acd_0 + g_yzz_0_xyyzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xyzzzzzzz_0[i] = g_zz_0_xyzzzzzzz_0[i] * fbe_0 - g_zz_0_xyzzzzzzz_1[i] * fz_be_0 + g_yzz_0_xzzzzzzz_1[i] * fi_acd_0 + g_yzz_0_xyzzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xzzzzzzzz_0[i] = g_zz_0_xzzzzzzzz_0[i] * fbe_0 - g_zz_0_xzzzzzzzz_1[i] * fz_be_0 + g_yzz_0_xzzzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_yyyyyyyyy_0[i] = g_yy_0_yyyyyyyyy_0[i] * fbe_0 - g_yy_0_yyyyyyyyy_1[i] * fz_be_0 + g_yyz_0_yyyyyyyyy_1[i] * wa_z[i];

        g_yyzz_0_yyyyyyyyz_0[i] = g_zz_0_yyyyyyyyz_0[i] * fbe_0 - g_zz_0_yyyyyyyyz_1[i] * fz_be_0 + 8.0 * g_yzz_0_yyyyyyyz_1[i] * fi_acd_0 + g_yzz_0_yyyyyyyyz_1[i] * wa_y[i];

        g_yyzz_0_yyyyyyyzz_0[i] = g_zz_0_yyyyyyyzz_0[i] * fbe_0 - g_zz_0_yyyyyyyzz_1[i] * fz_be_0 + 7.0 * g_yzz_0_yyyyyyzz_1[i] * fi_acd_0 + g_yzz_0_yyyyyyyzz_1[i] * wa_y[i];

        g_yyzz_0_yyyyyyzzz_0[i] = g_zz_0_yyyyyyzzz_0[i] * fbe_0 - g_zz_0_yyyyyyzzz_1[i] * fz_be_0 + 6.0 * g_yzz_0_yyyyyzzz_1[i] * fi_acd_0 + g_yzz_0_yyyyyyzzz_1[i] * wa_y[i];

        g_yyzz_0_yyyyyzzzz_0[i] = g_zz_0_yyyyyzzzz_0[i] * fbe_0 - g_zz_0_yyyyyzzzz_1[i] * fz_be_0 + 5.0 * g_yzz_0_yyyyzzzz_1[i] * fi_acd_0 + g_yzz_0_yyyyyzzzz_1[i] * wa_y[i];

        g_yyzz_0_yyyyzzzzz_0[i] = g_zz_0_yyyyzzzzz_0[i] * fbe_0 - g_zz_0_yyyyzzzzz_1[i] * fz_be_0 + 4.0 * g_yzz_0_yyyzzzzz_1[i] * fi_acd_0 + g_yzz_0_yyyyzzzzz_1[i] * wa_y[i];

        g_yyzz_0_yyyzzzzzz_0[i] = g_zz_0_yyyzzzzzz_0[i] * fbe_0 - g_zz_0_yyyzzzzzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_yyzzzzzz_1[i] * fi_acd_0 + g_yzz_0_yyyzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_yyzzzzzzz_0[i] = g_zz_0_yyzzzzzzz_0[i] * fbe_0 - g_zz_0_yyzzzzzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_yzzzzzzz_1[i] * fi_acd_0 + g_yzz_0_yyzzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_yzzzzzzzz_0[i] = g_zz_0_yzzzzzzzz_0[i] * fbe_0 - g_zz_0_yzzzzzzzz_1[i] * fz_be_0 + g_yzz_0_zzzzzzzz_1[i] * fi_acd_0 + g_yzz_0_yzzzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_zzzzzzzzz_0[i] = g_zz_0_zzzzzzzzz_0[i] * fbe_0 - g_zz_0_zzzzzzzzz_1[i] * fz_be_0 + g_yzz_0_zzzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 715-770 components of targeted buffer : GSM

    auto g_yzzz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 715);

    auto g_yzzz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 716);

    auto g_yzzz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 717);

    auto g_yzzz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 718);

    auto g_yzzz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 719);

    auto g_yzzz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 720);

    auto g_yzzz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 721);

    auto g_yzzz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 722);

    auto g_yzzz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 723);

    auto g_yzzz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 724);

    auto g_yzzz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 725);

    auto g_yzzz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 726);

    auto g_yzzz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 727);

    auto g_yzzz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 728);

    auto g_yzzz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 729);

    auto g_yzzz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 730);

    auto g_yzzz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 731);

    auto g_yzzz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 732);

    auto g_yzzz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 733);

    auto g_yzzz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 734);

    auto g_yzzz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 735);

    auto g_yzzz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 736);

    auto g_yzzz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 737);

    auto g_yzzz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 738);

    auto g_yzzz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 739);

    auto g_yzzz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 740);

    auto g_yzzz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 741);

    auto g_yzzz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 742);

    auto g_yzzz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 743);

    auto g_yzzz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 744);

    auto g_yzzz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 745);

    auto g_yzzz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 746);

    auto g_yzzz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 747);

    auto g_yzzz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 748);

    auto g_yzzz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 749);

    auto g_yzzz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 750);

    auto g_yzzz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 751);

    auto g_yzzz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 752);

    auto g_yzzz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 753);

    auto g_yzzz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 754);

    auto g_yzzz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 755);

    auto g_yzzz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 756);

    auto g_yzzz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 757);

    auto g_yzzz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 758);

    auto g_yzzz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 759);

    auto g_yzzz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 760);

    auto g_yzzz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 761);

    auto g_yzzz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 762);

    auto g_yzzz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 763);

    auto g_yzzz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 764);

    auto g_yzzz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 765);

    auto g_yzzz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 766);

    auto g_yzzz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 767);

    auto g_yzzz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 768);

    auto g_yzzz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 769);

    #pragma omp simd aligned(g_yzzz_0_xxxxxxxxx_0, g_yzzz_0_xxxxxxxxy_0, g_yzzz_0_xxxxxxxxz_0, g_yzzz_0_xxxxxxxyy_0, g_yzzz_0_xxxxxxxyz_0, g_yzzz_0_xxxxxxxzz_0, g_yzzz_0_xxxxxxyyy_0, g_yzzz_0_xxxxxxyyz_0, g_yzzz_0_xxxxxxyzz_0, g_yzzz_0_xxxxxxzzz_0, g_yzzz_0_xxxxxyyyy_0, g_yzzz_0_xxxxxyyyz_0, g_yzzz_0_xxxxxyyzz_0, g_yzzz_0_xxxxxyzzz_0, g_yzzz_0_xxxxxzzzz_0, g_yzzz_0_xxxxyyyyy_0, g_yzzz_0_xxxxyyyyz_0, g_yzzz_0_xxxxyyyzz_0, g_yzzz_0_xxxxyyzzz_0, g_yzzz_0_xxxxyzzzz_0, g_yzzz_0_xxxxzzzzz_0, g_yzzz_0_xxxyyyyyy_0, g_yzzz_0_xxxyyyyyz_0, g_yzzz_0_xxxyyyyzz_0, g_yzzz_0_xxxyyyzzz_0, g_yzzz_0_xxxyyzzzz_0, g_yzzz_0_xxxyzzzzz_0, g_yzzz_0_xxxzzzzzz_0, g_yzzz_0_xxyyyyyyy_0, g_yzzz_0_xxyyyyyyz_0, g_yzzz_0_xxyyyyyzz_0, g_yzzz_0_xxyyyyzzz_0, g_yzzz_0_xxyyyzzzz_0, g_yzzz_0_xxyyzzzzz_0, g_yzzz_0_xxyzzzzzz_0, g_yzzz_0_xxzzzzzzz_0, g_yzzz_0_xyyyyyyyy_0, g_yzzz_0_xyyyyyyyz_0, g_yzzz_0_xyyyyyyzz_0, g_yzzz_0_xyyyyyzzz_0, g_yzzz_0_xyyyyzzzz_0, g_yzzz_0_xyyyzzzzz_0, g_yzzz_0_xyyzzzzzz_0, g_yzzz_0_xyzzzzzzz_0, g_yzzz_0_xzzzzzzzz_0, g_yzzz_0_yyyyyyyyy_0, g_yzzz_0_yyyyyyyyz_0, g_yzzz_0_yyyyyyyzz_0, g_yzzz_0_yyyyyyzzz_0, g_yzzz_0_yyyyyzzzz_0, g_yzzz_0_yyyyzzzzz_0, g_yzzz_0_yyyzzzzzz_0, g_yzzz_0_yyzzzzzzz_0, g_yzzz_0_yzzzzzzzz_0, g_yzzz_0_zzzzzzzzz_0, g_zzz_0_xxxxxxxx_1, g_zzz_0_xxxxxxxxx_1, g_zzz_0_xxxxxxxxy_1, g_zzz_0_xxxxxxxxz_1, g_zzz_0_xxxxxxxy_1, g_zzz_0_xxxxxxxyy_1, g_zzz_0_xxxxxxxyz_1, g_zzz_0_xxxxxxxz_1, g_zzz_0_xxxxxxxzz_1, g_zzz_0_xxxxxxyy_1, g_zzz_0_xxxxxxyyy_1, g_zzz_0_xxxxxxyyz_1, g_zzz_0_xxxxxxyz_1, g_zzz_0_xxxxxxyzz_1, g_zzz_0_xxxxxxzz_1, g_zzz_0_xxxxxxzzz_1, g_zzz_0_xxxxxyyy_1, g_zzz_0_xxxxxyyyy_1, g_zzz_0_xxxxxyyyz_1, g_zzz_0_xxxxxyyz_1, g_zzz_0_xxxxxyyzz_1, g_zzz_0_xxxxxyzz_1, g_zzz_0_xxxxxyzzz_1, g_zzz_0_xxxxxzzz_1, g_zzz_0_xxxxxzzzz_1, g_zzz_0_xxxxyyyy_1, g_zzz_0_xxxxyyyyy_1, g_zzz_0_xxxxyyyyz_1, g_zzz_0_xxxxyyyz_1, g_zzz_0_xxxxyyyzz_1, g_zzz_0_xxxxyyzz_1, g_zzz_0_xxxxyyzzz_1, g_zzz_0_xxxxyzzz_1, g_zzz_0_xxxxyzzzz_1, g_zzz_0_xxxxzzzz_1, g_zzz_0_xxxxzzzzz_1, g_zzz_0_xxxyyyyy_1, g_zzz_0_xxxyyyyyy_1, g_zzz_0_xxxyyyyyz_1, g_zzz_0_xxxyyyyz_1, g_zzz_0_xxxyyyyzz_1, g_zzz_0_xxxyyyzz_1, g_zzz_0_xxxyyyzzz_1, g_zzz_0_xxxyyzzz_1, g_zzz_0_xxxyyzzzz_1, g_zzz_0_xxxyzzzz_1, g_zzz_0_xxxyzzzzz_1, g_zzz_0_xxxzzzzz_1, g_zzz_0_xxxzzzzzz_1, g_zzz_0_xxyyyyyy_1, g_zzz_0_xxyyyyyyy_1, g_zzz_0_xxyyyyyyz_1, g_zzz_0_xxyyyyyz_1, g_zzz_0_xxyyyyyzz_1, g_zzz_0_xxyyyyzz_1, g_zzz_0_xxyyyyzzz_1, g_zzz_0_xxyyyzzz_1, g_zzz_0_xxyyyzzzz_1, g_zzz_0_xxyyzzzz_1, g_zzz_0_xxyyzzzzz_1, g_zzz_0_xxyzzzzz_1, g_zzz_0_xxyzzzzzz_1, g_zzz_0_xxzzzzzz_1, g_zzz_0_xxzzzzzzz_1, g_zzz_0_xyyyyyyy_1, g_zzz_0_xyyyyyyyy_1, g_zzz_0_xyyyyyyyz_1, g_zzz_0_xyyyyyyz_1, g_zzz_0_xyyyyyyzz_1, g_zzz_0_xyyyyyzz_1, g_zzz_0_xyyyyyzzz_1, g_zzz_0_xyyyyzzz_1, g_zzz_0_xyyyyzzzz_1, g_zzz_0_xyyyzzzz_1, g_zzz_0_xyyyzzzzz_1, g_zzz_0_xyyzzzzz_1, g_zzz_0_xyyzzzzzz_1, g_zzz_0_xyzzzzzz_1, g_zzz_0_xyzzzzzzz_1, g_zzz_0_xzzzzzzz_1, g_zzz_0_xzzzzzzzz_1, g_zzz_0_yyyyyyyy_1, g_zzz_0_yyyyyyyyy_1, g_zzz_0_yyyyyyyyz_1, g_zzz_0_yyyyyyyz_1, g_zzz_0_yyyyyyyzz_1, g_zzz_0_yyyyyyzz_1, g_zzz_0_yyyyyyzzz_1, g_zzz_0_yyyyyzzz_1, g_zzz_0_yyyyyzzzz_1, g_zzz_0_yyyyzzzz_1, g_zzz_0_yyyyzzzzz_1, g_zzz_0_yyyzzzzz_1, g_zzz_0_yyyzzzzzz_1, g_zzz_0_yyzzzzzz_1, g_zzz_0_yyzzzzzzz_1, g_zzz_0_yzzzzzzz_1, g_zzz_0_yzzzzzzzz_1, g_zzz_0_zzzzzzzz_1, g_zzz_0_zzzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzz_0_xxxxxxxxx_0[i] = g_zzz_0_xxxxxxxxx_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxxxy_0[i] = g_zzz_0_xxxxxxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxxy_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxxxz_0[i] = g_zzz_0_xxxxxxxxz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxxyy_0[i] = 2.0 * g_zzz_0_xxxxxxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxyy_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxxyz_0[i] = g_zzz_0_xxxxxxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxyz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxxzz_0[i] = g_zzz_0_xxxxxxxzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxyyy_0[i] = 3.0 * g_zzz_0_xxxxxxyy_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyyy_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxyyz_0[i] = 2.0 * g_zzz_0_xxxxxxyz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyyz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxyzz_0[i] = g_zzz_0_xxxxxxzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxzzz_0[i] = g_zzz_0_xxxxxxzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxyyyy_0[i] = 4.0 * g_zzz_0_xxxxxyyy_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyyy_1[i] * wa_y[i];

        g_yzzz_0_xxxxxyyyz_0[i] = 3.0 * g_zzz_0_xxxxxyyz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyyz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxyyzz_0[i] = 2.0 * g_zzz_0_xxxxxyzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxyzzz_0[i] = g_zzz_0_xxxxxzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxzzzz_0[i] = g_zzz_0_xxxxxzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxyyyyy_0[i] = 5.0 * g_zzz_0_xxxxyyyy_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyyy_1[i] * wa_y[i];

        g_yzzz_0_xxxxyyyyz_0[i] = 4.0 * g_zzz_0_xxxxyyyz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyyz_1[i] * wa_y[i];

        g_yzzz_0_xxxxyyyzz_0[i] = 3.0 * g_zzz_0_xxxxyyzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxyyzzz_0[i] = 2.0 * g_zzz_0_xxxxyzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxyzzzz_0[i] = g_zzz_0_xxxxzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxzzzzz_0[i] = g_zzz_0_xxxxzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxyyyyyy_0[i] = 6.0 * g_zzz_0_xxxyyyyy_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyyy_1[i] * wa_y[i];

        g_yzzz_0_xxxyyyyyz_0[i] = 5.0 * g_zzz_0_xxxyyyyz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyyz_1[i] * wa_y[i];

        g_yzzz_0_xxxyyyyzz_0[i] = 4.0 * g_zzz_0_xxxyyyzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyzz_1[i] * wa_y[i];

        g_yzzz_0_xxxyyyzzz_0[i] = 3.0 * g_zzz_0_xxxyyzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxyyzzzz_0[i] = 2.0 * g_zzz_0_xxxyzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxyzzzzz_0[i] = g_zzz_0_xxxzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxzzzzzz_0[i] = g_zzz_0_xxxzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxyyyyyyy_0[i] = 7.0 * g_zzz_0_xxyyyyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyyy_1[i] * wa_y[i];

        g_yzzz_0_xxyyyyyyz_0[i] = 6.0 * g_zzz_0_xxyyyyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyyz_1[i] * wa_y[i];

        g_yzzz_0_xxyyyyyzz_0[i] = 5.0 * g_zzz_0_xxyyyyzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyzz_1[i] * wa_y[i];

        g_yzzz_0_xxyyyyzzz_0[i] = 4.0 * g_zzz_0_xxyyyzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyzzz_1[i] * wa_y[i];

        g_yzzz_0_xxyyyzzzz_0[i] = 3.0 * g_zzz_0_xxyyzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxyyzzzzz_0[i] = 2.0 * g_zzz_0_xxyzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxyzzzzzz_0[i] = g_zzz_0_xxzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxzzzzzzz_0[i] = g_zzz_0_xxzzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xyyyyyyyy_0[i] = 8.0 * g_zzz_0_xyyyyyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyyy_1[i] * wa_y[i];

        g_yzzz_0_xyyyyyyyz_0[i] = 7.0 * g_zzz_0_xyyyyyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyyz_1[i] * wa_y[i];

        g_yzzz_0_xyyyyyyzz_0[i] = 6.0 * g_zzz_0_xyyyyyzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyzz_1[i] * wa_y[i];

        g_yzzz_0_xyyyyyzzz_0[i] = 5.0 * g_zzz_0_xyyyyzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyzzz_1[i] * wa_y[i];

        g_yzzz_0_xyyyyzzzz_0[i] = 4.0 * g_zzz_0_xyyyzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyzzzz_1[i] * wa_y[i];

        g_yzzz_0_xyyyzzzzz_0[i] = 3.0 * g_zzz_0_xyyzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xyyzzzzzz_0[i] = 2.0 * g_zzz_0_xyzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xyzzzzzzz_0[i] = g_zzz_0_xzzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xzzzzzzzz_0[i] = g_zzz_0_xzzzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_yyyyyyyyy_0[i] = 9.0 * g_zzz_0_yyyyyyyy_1[i] * fi_acd_0 + g_zzz_0_yyyyyyyyy_1[i] * wa_y[i];

        g_yzzz_0_yyyyyyyyz_0[i] = 8.0 * g_zzz_0_yyyyyyyz_1[i] * fi_acd_0 + g_zzz_0_yyyyyyyyz_1[i] * wa_y[i];

        g_yzzz_0_yyyyyyyzz_0[i] = 7.0 * g_zzz_0_yyyyyyzz_1[i] * fi_acd_0 + g_zzz_0_yyyyyyyzz_1[i] * wa_y[i];

        g_yzzz_0_yyyyyyzzz_0[i] = 6.0 * g_zzz_0_yyyyyzzz_1[i] * fi_acd_0 + g_zzz_0_yyyyyyzzz_1[i] * wa_y[i];

        g_yzzz_0_yyyyyzzzz_0[i] = 5.0 * g_zzz_0_yyyyzzzz_1[i] * fi_acd_0 + g_zzz_0_yyyyyzzzz_1[i] * wa_y[i];

        g_yzzz_0_yyyyzzzzz_0[i] = 4.0 * g_zzz_0_yyyzzzzz_1[i] * fi_acd_0 + g_zzz_0_yyyyzzzzz_1[i] * wa_y[i];

        g_yzzz_0_yyyzzzzzz_0[i] = 3.0 * g_zzz_0_yyzzzzzz_1[i] * fi_acd_0 + g_zzz_0_yyyzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_yyzzzzzzz_0[i] = 2.0 * g_zzz_0_yzzzzzzz_1[i] * fi_acd_0 + g_zzz_0_yyzzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_yzzzzzzzz_0[i] = g_zzz_0_zzzzzzzz_1[i] * fi_acd_0 + g_zzz_0_yzzzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_zzzzzzzzz_0[i] = g_zzz_0_zzzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 770-825 components of targeted buffer : GSM

    auto g_zzzz_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_gsm + 770);

    auto g_zzzz_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_gsm + 771);

    auto g_zzzz_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_gsm + 772);

    auto g_zzzz_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_gsm + 773);

    auto g_zzzz_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_gsm + 774);

    auto g_zzzz_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_gsm + 775);

    auto g_zzzz_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_gsm + 776);

    auto g_zzzz_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_gsm + 777);

    auto g_zzzz_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_gsm + 778);

    auto g_zzzz_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_gsm + 779);

    auto g_zzzz_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_gsm + 780);

    auto g_zzzz_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_gsm + 781);

    auto g_zzzz_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_gsm + 782);

    auto g_zzzz_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_gsm + 783);

    auto g_zzzz_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_gsm + 784);

    auto g_zzzz_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 785);

    auto g_zzzz_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 786);

    auto g_zzzz_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 787);

    auto g_zzzz_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 788);

    auto g_zzzz_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 789);

    auto g_zzzz_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 790);

    auto g_zzzz_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 791);

    auto g_zzzz_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 792);

    auto g_zzzz_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 793);

    auto g_zzzz_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 794);

    auto g_zzzz_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 795);

    auto g_zzzz_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 796);

    auto g_zzzz_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 797);

    auto g_zzzz_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 798);

    auto g_zzzz_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 799);

    auto g_zzzz_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 800);

    auto g_zzzz_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 801);

    auto g_zzzz_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 802);

    auto g_zzzz_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 803);

    auto g_zzzz_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 804);

    auto g_zzzz_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 805);

    auto g_zzzz_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 806);

    auto g_zzzz_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 807);

    auto g_zzzz_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 808);

    auto g_zzzz_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 809);

    auto g_zzzz_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 810);

    auto g_zzzz_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 811);

    auto g_zzzz_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 812);

    auto g_zzzz_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 813);

    auto g_zzzz_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 814);

    auto g_zzzz_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_gsm + 815);

    auto g_zzzz_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_gsm + 816);

    auto g_zzzz_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_gsm + 817);

    auto g_zzzz_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_gsm + 818);

    auto g_zzzz_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_gsm + 819);

    auto g_zzzz_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 820);

    auto g_zzzz_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 821);

    auto g_zzzz_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 822);

    auto g_zzzz_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 823);

    auto g_zzzz_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_gsm + 824);

    #pragma omp simd aligned(g_zz_0_xxxxxxxxx_0, g_zz_0_xxxxxxxxx_1, g_zz_0_xxxxxxxxy_0, g_zz_0_xxxxxxxxy_1, g_zz_0_xxxxxxxxz_0, g_zz_0_xxxxxxxxz_1, g_zz_0_xxxxxxxyy_0, g_zz_0_xxxxxxxyy_1, g_zz_0_xxxxxxxyz_0, g_zz_0_xxxxxxxyz_1, g_zz_0_xxxxxxxzz_0, g_zz_0_xxxxxxxzz_1, g_zz_0_xxxxxxyyy_0, g_zz_0_xxxxxxyyy_1, g_zz_0_xxxxxxyyz_0, g_zz_0_xxxxxxyyz_1, g_zz_0_xxxxxxyzz_0, g_zz_0_xxxxxxyzz_1, g_zz_0_xxxxxxzzz_0, g_zz_0_xxxxxxzzz_1, g_zz_0_xxxxxyyyy_0, g_zz_0_xxxxxyyyy_1, g_zz_0_xxxxxyyyz_0, g_zz_0_xxxxxyyyz_1, g_zz_0_xxxxxyyzz_0, g_zz_0_xxxxxyyzz_1, g_zz_0_xxxxxyzzz_0, g_zz_0_xxxxxyzzz_1, g_zz_0_xxxxxzzzz_0, g_zz_0_xxxxxzzzz_1, g_zz_0_xxxxyyyyy_0, g_zz_0_xxxxyyyyy_1, g_zz_0_xxxxyyyyz_0, g_zz_0_xxxxyyyyz_1, g_zz_0_xxxxyyyzz_0, g_zz_0_xxxxyyyzz_1, g_zz_0_xxxxyyzzz_0, g_zz_0_xxxxyyzzz_1, g_zz_0_xxxxyzzzz_0, g_zz_0_xxxxyzzzz_1, g_zz_0_xxxxzzzzz_0, g_zz_0_xxxxzzzzz_1, g_zz_0_xxxyyyyyy_0, g_zz_0_xxxyyyyyy_1, g_zz_0_xxxyyyyyz_0, g_zz_0_xxxyyyyyz_1, g_zz_0_xxxyyyyzz_0, g_zz_0_xxxyyyyzz_1, g_zz_0_xxxyyyzzz_0, g_zz_0_xxxyyyzzz_1, g_zz_0_xxxyyzzzz_0, g_zz_0_xxxyyzzzz_1, g_zz_0_xxxyzzzzz_0, g_zz_0_xxxyzzzzz_1, g_zz_0_xxxzzzzzz_0, g_zz_0_xxxzzzzzz_1, g_zz_0_xxyyyyyyy_0, g_zz_0_xxyyyyyyy_1, g_zz_0_xxyyyyyyz_0, g_zz_0_xxyyyyyyz_1, g_zz_0_xxyyyyyzz_0, g_zz_0_xxyyyyyzz_1, g_zz_0_xxyyyyzzz_0, g_zz_0_xxyyyyzzz_1, g_zz_0_xxyyyzzzz_0, g_zz_0_xxyyyzzzz_1, g_zz_0_xxyyzzzzz_0, g_zz_0_xxyyzzzzz_1, g_zz_0_xxyzzzzzz_0, g_zz_0_xxyzzzzzz_1, g_zz_0_xxzzzzzzz_0, g_zz_0_xxzzzzzzz_1, g_zz_0_xyyyyyyyy_0, g_zz_0_xyyyyyyyy_1, g_zz_0_xyyyyyyyz_0, g_zz_0_xyyyyyyyz_1, g_zz_0_xyyyyyyzz_0, g_zz_0_xyyyyyyzz_1, g_zz_0_xyyyyyzzz_0, g_zz_0_xyyyyyzzz_1, g_zz_0_xyyyyzzzz_0, g_zz_0_xyyyyzzzz_1, g_zz_0_xyyyzzzzz_0, g_zz_0_xyyyzzzzz_1, g_zz_0_xyyzzzzzz_0, g_zz_0_xyyzzzzzz_1, g_zz_0_xyzzzzzzz_0, g_zz_0_xyzzzzzzz_1, g_zz_0_xzzzzzzzz_0, g_zz_0_xzzzzzzzz_1, g_zz_0_yyyyyyyyy_0, g_zz_0_yyyyyyyyy_1, g_zz_0_yyyyyyyyz_0, g_zz_0_yyyyyyyyz_1, g_zz_0_yyyyyyyzz_0, g_zz_0_yyyyyyyzz_1, g_zz_0_yyyyyyzzz_0, g_zz_0_yyyyyyzzz_1, g_zz_0_yyyyyzzzz_0, g_zz_0_yyyyyzzzz_1, g_zz_0_yyyyzzzzz_0, g_zz_0_yyyyzzzzz_1, g_zz_0_yyyzzzzzz_0, g_zz_0_yyyzzzzzz_1, g_zz_0_yyzzzzzzz_0, g_zz_0_yyzzzzzzz_1, g_zz_0_yzzzzzzzz_0, g_zz_0_yzzzzzzzz_1, g_zz_0_zzzzzzzzz_0, g_zz_0_zzzzzzzzz_1, g_zzz_0_xxxxxxxx_1, g_zzz_0_xxxxxxxxx_1, g_zzz_0_xxxxxxxxy_1, g_zzz_0_xxxxxxxxz_1, g_zzz_0_xxxxxxxy_1, g_zzz_0_xxxxxxxyy_1, g_zzz_0_xxxxxxxyz_1, g_zzz_0_xxxxxxxz_1, g_zzz_0_xxxxxxxzz_1, g_zzz_0_xxxxxxyy_1, g_zzz_0_xxxxxxyyy_1, g_zzz_0_xxxxxxyyz_1, g_zzz_0_xxxxxxyz_1, g_zzz_0_xxxxxxyzz_1, g_zzz_0_xxxxxxzz_1, g_zzz_0_xxxxxxzzz_1, g_zzz_0_xxxxxyyy_1, g_zzz_0_xxxxxyyyy_1, g_zzz_0_xxxxxyyyz_1, g_zzz_0_xxxxxyyz_1, g_zzz_0_xxxxxyyzz_1, g_zzz_0_xxxxxyzz_1, g_zzz_0_xxxxxyzzz_1, g_zzz_0_xxxxxzzz_1, g_zzz_0_xxxxxzzzz_1, g_zzz_0_xxxxyyyy_1, g_zzz_0_xxxxyyyyy_1, g_zzz_0_xxxxyyyyz_1, g_zzz_0_xxxxyyyz_1, g_zzz_0_xxxxyyyzz_1, g_zzz_0_xxxxyyzz_1, g_zzz_0_xxxxyyzzz_1, g_zzz_0_xxxxyzzz_1, g_zzz_0_xxxxyzzzz_1, g_zzz_0_xxxxzzzz_1, g_zzz_0_xxxxzzzzz_1, g_zzz_0_xxxyyyyy_1, g_zzz_0_xxxyyyyyy_1, g_zzz_0_xxxyyyyyz_1, g_zzz_0_xxxyyyyz_1, g_zzz_0_xxxyyyyzz_1, g_zzz_0_xxxyyyzz_1, g_zzz_0_xxxyyyzzz_1, g_zzz_0_xxxyyzzz_1, g_zzz_0_xxxyyzzzz_1, g_zzz_0_xxxyzzzz_1, g_zzz_0_xxxyzzzzz_1, g_zzz_0_xxxzzzzz_1, g_zzz_0_xxxzzzzzz_1, g_zzz_0_xxyyyyyy_1, g_zzz_0_xxyyyyyyy_1, g_zzz_0_xxyyyyyyz_1, g_zzz_0_xxyyyyyz_1, g_zzz_0_xxyyyyyzz_1, g_zzz_0_xxyyyyzz_1, g_zzz_0_xxyyyyzzz_1, g_zzz_0_xxyyyzzz_1, g_zzz_0_xxyyyzzzz_1, g_zzz_0_xxyyzzzz_1, g_zzz_0_xxyyzzzzz_1, g_zzz_0_xxyzzzzz_1, g_zzz_0_xxyzzzzzz_1, g_zzz_0_xxzzzzzz_1, g_zzz_0_xxzzzzzzz_1, g_zzz_0_xyyyyyyy_1, g_zzz_0_xyyyyyyyy_1, g_zzz_0_xyyyyyyyz_1, g_zzz_0_xyyyyyyz_1, g_zzz_0_xyyyyyyzz_1, g_zzz_0_xyyyyyzz_1, g_zzz_0_xyyyyyzzz_1, g_zzz_0_xyyyyzzz_1, g_zzz_0_xyyyyzzzz_1, g_zzz_0_xyyyzzzz_1, g_zzz_0_xyyyzzzzz_1, g_zzz_0_xyyzzzzz_1, g_zzz_0_xyyzzzzzz_1, g_zzz_0_xyzzzzzz_1, g_zzz_0_xyzzzzzzz_1, g_zzz_0_xzzzzzzz_1, g_zzz_0_xzzzzzzzz_1, g_zzz_0_yyyyyyyy_1, g_zzz_0_yyyyyyyyy_1, g_zzz_0_yyyyyyyyz_1, g_zzz_0_yyyyyyyz_1, g_zzz_0_yyyyyyyzz_1, g_zzz_0_yyyyyyzz_1, g_zzz_0_yyyyyyzzz_1, g_zzz_0_yyyyyzzz_1, g_zzz_0_yyyyyzzzz_1, g_zzz_0_yyyyzzzz_1, g_zzz_0_yyyyzzzzz_1, g_zzz_0_yyyzzzzz_1, g_zzz_0_yyyzzzzzz_1, g_zzz_0_yyzzzzzz_1, g_zzz_0_yyzzzzzzz_1, g_zzz_0_yzzzzzzz_1, g_zzz_0_yzzzzzzzz_1, g_zzz_0_zzzzzzzz_1, g_zzz_0_zzzzzzzzz_1, g_zzzz_0_xxxxxxxxx_0, g_zzzz_0_xxxxxxxxy_0, g_zzzz_0_xxxxxxxxz_0, g_zzzz_0_xxxxxxxyy_0, g_zzzz_0_xxxxxxxyz_0, g_zzzz_0_xxxxxxxzz_0, g_zzzz_0_xxxxxxyyy_0, g_zzzz_0_xxxxxxyyz_0, g_zzzz_0_xxxxxxyzz_0, g_zzzz_0_xxxxxxzzz_0, g_zzzz_0_xxxxxyyyy_0, g_zzzz_0_xxxxxyyyz_0, g_zzzz_0_xxxxxyyzz_0, g_zzzz_0_xxxxxyzzz_0, g_zzzz_0_xxxxxzzzz_0, g_zzzz_0_xxxxyyyyy_0, g_zzzz_0_xxxxyyyyz_0, g_zzzz_0_xxxxyyyzz_0, g_zzzz_0_xxxxyyzzz_0, g_zzzz_0_xxxxyzzzz_0, g_zzzz_0_xxxxzzzzz_0, g_zzzz_0_xxxyyyyyy_0, g_zzzz_0_xxxyyyyyz_0, g_zzzz_0_xxxyyyyzz_0, g_zzzz_0_xxxyyyzzz_0, g_zzzz_0_xxxyyzzzz_0, g_zzzz_0_xxxyzzzzz_0, g_zzzz_0_xxxzzzzzz_0, g_zzzz_0_xxyyyyyyy_0, g_zzzz_0_xxyyyyyyz_0, g_zzzz_0_xxyyyyyzz_0, g_zzzz_0_xxyyyyzzz_0, g_zzzz_0_xxyyyzzzz_0, g_zzzz_0_xxyyzzzzz_0, g_zzzz_0_xxyzzzzzz_0, g_zzzz_0_xxzzzzzzz_0, g_zzzz_0_xyyyyyyyy_0, g_zzzz_0_xyyyyyyyz_0, g_zzzz_0_xyyyyyyzz_0, g_zzzz_0_xyyyyyzzz_0, g_zzzz_0_xyyyyzzzz_0, g_zzzz_0_xyyyzzzzz_0, g_zzzz_0_xyyzzzzzz_0, g_zzzz_0_xyzzzzzzz_0, g_zzzz_0_xzzzzzzzz_0, g_zzzz_0_yyyyyyyyy_0, g_zzzz_0_yyyyyyyyz_0, g_zzzz_0_yyyyyyyzz_0, g_zzzz_0_yyyyyyzzz_0, g_zzzz_0_yyyyyzzzz_0, g_zzzz_0_yyyyzzzzz_0, g_zzzz_0_yyyzzzzzz_0, g_zzzz_0_yyzzzzzzz_0, g_zzzz_0_yzzzzzzzz_0, g_zzzz_0_zzzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzz_0_xxxxxxxxx_0[i] = 3.0 * g_zz_0_xxxxxxxxx_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxxxx_1[i] * fz_be_0 + g_zzz_0_xxxxxxxxx_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxxxy_0[i] = 3.0 * g_zz_0_xxxxxxxxy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxxxy_1[i] * fz_be_0 + g_zzz_0_xxxxxxxxy_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxxxz_0[i] = 3.0 * g_zz_0_xxxxxxxxz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxxxz_1[i] * fz_be_0 + g_zzz_0_xxxxxxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxxz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxxyy_0[i] = 3.0 * g_zz_0_xxxxxxxyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxxyy_1[i] * fz_be_0 + g_zzz_0_xxxxxxxyy_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxxyz_0[i] = 3.0 * g_zz_0_xxxxxxxyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxxyz_1[i] * fz_be_0 + g_zzz_0_xxxxxxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxyz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxxzz_0[i] = 3.0 * g_zz_0_xxxxxxxzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxxzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxxxxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxyyy_0[i] = 3.0 * g_zz_0_xxxxxxyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxyyy_1[i] * fz_be_0 + g_zzz_0_xxxxxxyyy_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxyyz_0[i] = 3.0 * g_zz_0_xxxxxxyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxyyz_1[i] * fz_be_0 + g_zzz_0_xxxxxxyy_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyyz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxyzz_0[i] = 3.0 * g_zz_0_xxxxxxyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxxxxyz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxzzz_0[i] = 3.0 * g_zz_0_xxxxxxzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxxxxxzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxyyyy_0[i] = 3.0 * g_zz_0_xxxxxyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxyyyy_1[i] * fz_be_0 + g_zzz_0_xxxxxyyyy_1[i] * wa_z[i];

        g_zzzz_0_xxxxxyyyz_0[i] = 3.0 * g_zz_0_xxxxxyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxyyyz_1[i] * fz_be_0 + g_zzz_0_xxxxxyyy_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyyz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxyyzz_0[i] = 3.0 * g_zz_0_xxxxxyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxxxyyz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxyzzz_0[i] = 3.0 * g_zz_0_xxxxxyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxxxxyzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxzzzz_0[i] = 3.0 * g_zz_0_xxxxxzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xxxxxzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxyyyyy_0[i] = 3.0 * g_zz_0_xxxxyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyyyyy_1[i] * fz_be_0 + g_zzz_0_xxxxyyyyy_1[i] * wa_z[i];

        g_zzzz_0_xxxxyyyyz_0[i] = 3.0 * g_zz_0_xxxxyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyyyyz_1[i] * fz_be_0 + g_zzz_0_xxxxyyyy_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyyz_1[i] * wa_z[i];

        g_zzzz_0_xxxxyyyzz_0[i] = 3.0 * g_zz_0_xxxxyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxxyyyz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxyyzzz_0[i] = 3.0 * g_zz_0_xxxxyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxxxyyzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxyzzzz_0[i] = 3.0 * g_zz_0_xxxxyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xxxxyzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxzzzzz_0[i] = 3.0 * g_zz_0_xxxxzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_xxxxzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxyyyyyy_0[i] = 3.0 * g_zz_0_xxxyyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyyyyy_1[i] * fz_be_0 + g_zzz_0_xxxyyyyyy_1[i] * wa_z[i];

        g_zzzz_0_xxxyyyyyz_0[i] = 3.0 * g_zz_0_xxxyyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyyyyz_1[i] * fz_be_0 + g_zzz_0_xxxyyyyy_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyyz_1[i] * wa_z[i];

        g_zzzz_0_xxxyyyyzz_0[i] = 3.0 * g_zz_0_xxxyyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxyyyyz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyzz_1[i] * wa_z[i];

        g_zzzz_0_xxxyyyzzz_0[i] = 3.0 * g_zz_0_xxxyyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxxyyyzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxyyzzzz_0[i] = 3.0 * g_zz_0_xxxyyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xxxyyzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxyzzzzz_0[i] = 3.0 * g_zz_0_xxxyzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_xxxyzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxzzzzzz_0[i] = 3.0 * g_zz_0_xxxzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzz_0_xxxzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxyyyyyyy_0[i] = 3.0 * g_zz_0_xxyyyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyyyyy_1[i] * fz_be_0 + g_zzz_0_xxyyyyyyy_1[i] * wa_z[i];

        g_zzzz_0_xxyyyyyyz_0[i] = 3.0 * g_zz_0_xxyyyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyyyyz_1[i] * fz_be_0 + g_zzz_0_xxyyyyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyyz_1[i] * wa_z[i];

        g_zzzz_0_xxyyyyyzz_0[i] = 3.0 * g_zz_0_xxyyyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxyyyyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyzz_1[i] * wa_z[i];

        g_zzzz_0_xxyyyyzzz_0[i] = 3.0 * g_zz_0_xxyyyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxyyyyzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyzzz_1[i] * wa_z[i];

        g_zzzz_0_xxyyyzzzz_0[i] = 3.0 * g_zz_0_xxyyyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xxyyyzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxyyzzzzz_0[i] = 3.0 * g_zz_0_xxyyzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_xxyyzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxyzzzzzz_0[i] = 3.0 * g_zz_0_xxyzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzz_0_xxyzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxzzzzzzz_0[i] = 3.0 * g_zz_0_xxzzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzz_0_xxzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxzzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xyyyyyyyy_0[i] = 3.0 * g_zz_0_xyyyyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyyyyy_1[i] * fz_be_0 + g_zzz_0_xyyyyyyyy_1[i] * wa_z[i];

        g_zzzz_0_xyyyyyyyz_0[i] = 3.0 * g_zz_0_xyyyyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyyyyz_1[i] * fz_be_0 + g_zzz_0_xyyyyyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyyz_1[i] * wa_z[i];

        g_zzzz_0_xyyyyyyzz_0[i] = 3.0 * g_zz_0_xyyyyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xyyyyyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyzz_1[i] * wa_z[i];

        g_zzzz_0_xyyyyyzzz_0[i] = 3.0 * g_zz_0_xyyyyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xyyyyyzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyzzz_1[i] * wa_z[i];

        g_zzzz_0_xyyyyzzzz_0[i] = 3.0 * g_zz_0_xyyyyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xyyyyzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyzzzz_1[i] * wa_z[i];

        g_zzzz_0_xyyyzzzzz_0[i] = 3.0 * g_zz_0_xyyyzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_xyyyzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xyyzzzzzz_0[i] = 3.0 * g_zz_0_xyyzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzz_0_xyyzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xyzzzzzzz_0[i] = 3.0 * g_zz_0_xyzzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzz_0_xyzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xzzzzzzzz_0[i] = 3.0 * g_zz_0_xzzzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xzzzzzzzz_1[i] * fz_be_0 + 8.0 * g_zzz_0_xzzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xzzzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_yyyyyyyyy_0[i] = 3.0 * g_zz_0_yyyyyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyyyyy_1[i] * fz_be_0 + g_zzz_0_yyyyyyyyy_1[i] * wa_z[i];

        g_zzzz_0_yyyyyyyyz_0[i] = 3.0 * g_zz_0_yyyyyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyyyyz_1[i] * fz_be_0 + g_zzz_0_yyyyyyyy_1[i] * fi_acd_0 + g_zzz_0_yyyyyyyyz_1[i] * wa_z[i];

        g_zzzz_0_yyyyyyyzz_0[i] = 3.0 * g_zz_0_yyyyyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_yyyyyyyz_1[i] * fi_acd_0 + g_zzz_0_yyyyyyyzz_1[i] * wa_z[i];

        g_zzzz_0_yyyyyyzzz_0[i] = 3.0 * g_zz_0_yyyyyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_yyyyyyzz_1[i] * fi_acd_0 + g_zzz_0_yyyyyyzzz_1[i] * wa_z[i];

        g_zzzz_0_yyyyyzzzz_0[i] = 3.0 * g_zz_0_yyyyyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_yyyyyzzz_1[i] * fi_acd_0 + g_zzz_0_yyyyyzzzz_1[i] * wa_z[i];

        g_zzzz_0_yyyyzzzzz_0[i] = 3.0 * g_zz_0_yyyyzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_yyyyzzzz_1[i] * fi_acd_0 + g_zzz_0_yyyyzzzzz_1[i] * wa_z[i];

        g_zzzz_0_yyyzzzzzz_0[i] = 3.0 * g_zz_0_yyyzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzz_0_yyyzzzzz_1[i] * fi_acd_0 + g_zzz_0_yyyzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_yyzzzzzzz_0[i] = 3.0 * g_zz_0_yyzzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzz_0_yyzzzzzz_1[i] * fi_acd_0 + g_zzz_0_yyzzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_yzzzzzzzz_0[i] = 3.0 * g_zz_0_yzzzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yzzzzzzzz_1[i] * fz_be_0 + 8.0 * g_zzz_0_yzzzzzzz_1[i] * fi_acd_0 + g_zzz_0_yzzzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_zzzzzzzzz_0[i] = 3.0 * g_zz_0_zzzzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_zzzzzzzzz_1[i] * fz_be_0 + 9.0 * g_zzz_0_zzzzzzzz_1[i] * fi_acd_0 + g_zzz_0_zzzzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

