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

#include "ThreeCenterElectronRepulsionPrimRecISK.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_isk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_isk,
                                 size_t idx_eri_0_gsk,
                                 size_t idx_eri_1_gsk,
                                 size_t idx_eri_1_hsi,
                                 size_t idx_eri_1_hsk,
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

    /// Set up components of auxilary buffer : GSK

    auto g_xxxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk);

    auto g_xxxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 1);

    auto g_xxxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 2);

    auto g_xxxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 3);

    auto g_xxxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 4);

    auto g_xxxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 5);

    auto g_xxxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 6);

    auto g_xxxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 7);

    auto g_xxxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 8);

    auto g_xxxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 9);

    auto g_xxxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 10);

    auto g_xxxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 11);

    auto g_xxxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 12);

    auto g_xxxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 13);

    auto g_xxxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 14);

    auto g_xxxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 15);

    auto g_xxxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 16);

    auto g_xxxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 17);

    auto g_xxxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 18);

    auto g_xxxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 19);

    auto g_xxxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 20);

    auto g_xxxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 21);

    auto g_xxxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 22);

    auto g_xxxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 23);

    auto g_xxxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 24);

    auto g_xxxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 25);

    auto g_xxxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 26);

    auto g_xxxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 27);

    auto g_xxxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 28);

    auto g_xxxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 29);

    auto g_xxxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 30);

    auto g_xxxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 31);

    auto g_xxxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 32);

    auto g_xxxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 33);

    auto g_xxxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 34);

    auto g_xxxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 35);

    auto g_xxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 36);

    auto g_xxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 38);

    auto g_xxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 41);

    auto g_xxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 45);

    auto g_xxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 50);

    auto g_xxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 56);

    auto g_xxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 63);

    auto g_xxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 72);

    auto g_xxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 73);

    auto g_xxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 75);

    auto g_xxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 78);

    auto g_xxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 82);

    auto g_xxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 87);

    auto g_xxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 93);

    auto g_xxyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 108);

    auto g_xxyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 109);

    auto g_xxyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 110);

    auto g_xxyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 111);

    auto g_xxyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 112);

    auto g_xxyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 113);

    auto g_xxyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 114);

    auto g_xxyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 115);

    auto g_xxyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 116);

    auto g_xxyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 117);

    auto g_xxyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 118);

    auto g_xxyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 119);

    auto g_xxyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 120);

    auto g_xxyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 121);

    auto g_xxyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 122);

    auto g_xxyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 123);

    auto g_xxyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 124);

    auto g_xxyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 125);

    auto g_xxyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 126);

    auto g_xxyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 127);

    auto g_xxyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 128);

    auto g_xxyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 129);

    auto g_xxyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 130);

    auto g_xxyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 131);

    auto g_xxyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 132);

    auto g_xxyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 133);

    auto g_xxyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 134);

    auto g_xxyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 135);

    auto g_xxyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 136);

    auto g_xxyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 137);

    auto g_xxyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 138);

    auto g_xxyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 139);

    auto g_xxyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 140);

    auto g_xxyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 141);

    auto g_xxyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 142);

    auto g_xxyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 143);

    auto g_xxzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 180);

    auto g_xxzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 181);

    auto g_xxzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 182);

    auto g_xxzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 183);

    auto g_xxzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 184);

    auto g_xxzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 185);

    auto g_xxzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 186);

    auto g_xxzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 187);

    auto g_xxzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 188);

    auto g_xxzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 189);

    auto g_xxzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 190);

    auto g_xxzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 191);

    auto g_xxzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 192);

    auto g_xxzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 193);

    auto g_xxzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 194);

    auto g_xxzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 195);

    auto g_xxzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 196);

    auto g_xxzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 197);

    auto g_xxzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 198);

    auto g_xxzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 199);

    auto g_xxzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 200);

    auto g_xxzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 201);

    auto g_xxzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 202);

    auto g_xxzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 203);

    auto g_xxzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 204);

    auto g_xxzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 205);

    auto g_xxzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 206);

    auto g_xxzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 207);

    auto g_xxzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 208);

    auto g_xxzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 209);

    auto g_xxzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 210);

    auto g_xxzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 211);

    auto g_xxzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 212);

    auto g_xxzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 213);

    auto g_xxzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 214);

    auto g_xxzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 215);

    auto g_xyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 217);

    auto g_xyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 219);

    auto g_xyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 220);

    auto g_xyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 222);

    auto g_xyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 223);

    auto g_xyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 224);

    auto g_xyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 226);

    auto g_xyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 227);

    auto g_xyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 228);

    auto g_xyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 229);

    auto g_xyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 231);

    auto g_xyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 232);

    auto g_xyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 233);

    auto g_xyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 234);

    auto g_xyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 235);

    auto g_xyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 237);

    auto g_xyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 238);

    auto g_xyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 239);

    auto g_xyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 240);

    auto g_xyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 241);

    auto g_xyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 242);

    auto g_xyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 244);

    auto g_xyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 245);

    auto g_xyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 246);

    auto g_xyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 247);

    auto g_xyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 248);

    auto g_xyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 249);

    auto g_xyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 250);

    auto g_xyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 251);

    auto g_xzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 326);

    auto g_xzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 328);

    auto g_xzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 329);

    auto g_xzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 331);

    auto g_xzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 332);

    auto g_xzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 333);

    auto g_xzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 335);

    auto g_xzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 336);

    auto g_xzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 337);

    auto g_xzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 338);

    auto g_xzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 340);

    auto g_xzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 341);

    auto g_xzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 342);

    auto g_xzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 343);

    auto g_xzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 344);

    auto g_xzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 346);

    auto g_xzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 347);

    auto g_xzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 348);

    auto g_xzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 349);

    auto g_xzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 350);

    auto g_xzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 351);

    auto g_xzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 352);

    auto g_xzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 353);

    auto g_xzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 354);

    auto g_xzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 355);

    auto g_xzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 356);

    auto g_xzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 357);

    auto g_xzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 358);

    auto g_xzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 359);

    auto g_yyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 360);

    auto g_yyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 361);

    auto g_yyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 362);

    auto g_yyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 363);

    auto g_yyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 364);

    auto g_yyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 365);

    auto g_yyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 366);

    auto g_yyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 367);

    auto g_yyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 368);

    auto g_yyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 369);

    auto g_yyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 370);

    auto g_yyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 371);

    auto g_yyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 372);

    auto g_yyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 373);

    auto g_yyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 374);

    auto g_yyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 375);

    auto g_yyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 376);

    auto g_yyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 377);

    auto g_yyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 378);

    auto g_yyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 379);

    auto g_yyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 380);

    auto g_yyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 381);

    auto g_yyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 382);

    auto g_yyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 383);

    auto g_yyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 384);

    auto g_yyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 385);

    auto g_yyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 386);

    auto g_yyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 387);

    auto g_yyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 388);

    auto g_yyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 389);

    auto g_yyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 390);

    auto g_yyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 391);

    auto g_yyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 392);

    auto g_yyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 393);

    auto g_yyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 394);

    auto g_yyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 395);

    auto g_yyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 397);

    auto g_yyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 399);

    auto g_yyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 402);

    auto g_yyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 406);

    auto g_yyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 411);

    auto g_yyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 417);

    auto g_yyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 424);

    auto g_yyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 432);

    auto g_yyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 433);

    auto g_yyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 434);

    auto g_yyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 435);

    auto g_yyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 436);

    auto g_yyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 437);

    auto g_yyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 438);

    auto g_yyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 439);

    auto g_yyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 440);

    auto g_yyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 441);

    auto g_yyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 442);

    auto g_yyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 443);

    auto g_yyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 444);

    auto g_yyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 445);

    auto g_yyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 446);

    auto g_yyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 447);

    auto g_yyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 448);

    auto g_yyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 449);

    auto g_yyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 450);

    auto g_yyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 451);

    auto g_yyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 452);

    auto g_yyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 453);

    auto g_yyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 454);

    auto g_yyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 455);

    auto g_yyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 456);

    auto g_yyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 457);

    auto g_yyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 458);

    auto g_yyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 459);

    auto g_yyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 460);

    auto g_yyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 461);

    auto g_yyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 462);

    auto g_yyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 463);

    auto g_yyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 464);

    auto g_yyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 465);

    auto g_yyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 466);

    auto g_yyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 467);

    auto g_yzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 468);

    auto g_yzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 470);

    auto g_yzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 472);

    auto g_yzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 473);

    auto g_yzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 475);

    auto g_yzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 476);

    auto g_yzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 477);

    auto g_yzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 479);

    auto g_yzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 480);

    auto g_yzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 481);

    auto g_yzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 482);

    auto g_yzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 484);

    auto g_yzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 485);

    auto g_yzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 486);

    auto g_yzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 487);

    auto g_yzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 488);

    auto g_yzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 490);

    auto g_yzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 491);

    auto g_yzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 492);

    auto g_yzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 493);

    auto g_yzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 494);

    auto g_yzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 495);

    auto g_yzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 497);

    auto g_yzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 498);

    auto g_yzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 499);

    auto g_yzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 500);

    auto g_yzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 501);

    auto g_yzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 502);

    auto g_yzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 503);

    auto g_zzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 504);

    auto g_zzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 505);

    auto g_zzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 506);

    auto g_zzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 507);

    auto g_zzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 508);

    auto g_zzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 509);

    auto g_zzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 510);

    auto g_zzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 511);

    auto g_zzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 512);

    auto g_zzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 513);

    auto g_zzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 514);

    auto g_zzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 515);

    auto g_zzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 516);

    auto g_zzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 517);

    auto g_zzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 518);

    auto g_zzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 519);

    auto g_zzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 520);

    auto g_zzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 521);

    auto g_zzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 522);

    auto g_zzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 523);

    auto g_zzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 524);

    auto g_zzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 525);

    auto g_zzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 526);

    auto g_zzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 527);

    auto g_zzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 528);

    auto g_zzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 529);

    auto g_zzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 530);

    auto g_zzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 531);

    auto g_zzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 532);

    auto g_zzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 533);

    auto g_zzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 534);

    auto g_zzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 535);

    auto g_zzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 536);

    auto g_zzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 537);

    auto g_zzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 538);

    auto g_zzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 539);

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

    auto g_xxxy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 36);

    auto g_xxxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 38);

    auto g_xxxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 41);

    auto g_xxxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 45);

    auto g_xxxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 50);

    auto g_xxxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 56);

    auto g_xxxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 63);

    auto g_xxxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 72);

    auto g_xxxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 73);

    auto g_xxxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 75);

    auto g_xxxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 78);

    auto g_xxxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 82);

    auto g_xxxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 87);

    auto g_xxxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 93);

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

    auto g_xyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 251);

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

    auto g_xzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 352);

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

    auto g_yyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_gsk + 397);

    auto g_yyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_gsk + 399);

    auto g_yyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_gsk + 402);

    auto g_yyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_gsk + 406);

    auto g_yyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 411);

    auto g_yyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 417);

    auto g_yyyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_gsk + 424);

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

    auto g_yzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_gsk + 468);

    auto g_yzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_gsk + 470);

    auto g_yzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_gsk + 472);

    auto g_yzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_gsk + 473);

    auto g_yzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_gsk + 475);

    auto g_yzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_gsk + 476);

    auto g_yzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_gsk + 477);

    auto g_yzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_gsk + 479);

    auto g_yzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_gsk + 480);

    auto g_yzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_gsk + 481);

    auto g_yzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_gsk + 482);

    auto g_yzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 484);

    auto g_yzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 485);

    auto g_yzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 486);

    auto g_yzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 487);

    auto g_yzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 488);

    auto g_yzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_gsk + 490);

    auto g_yzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_gsk + 491);

    auto g_yzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_gsk + 492);

    auto g_yzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_gsk + 493);

    auto g_yzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 494);

    auto g_yzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_gsk + 495);

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

    auto g_xxxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 58);

    auto g_xxxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 60);

    auto g_xxxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 61);

    auto g_xxxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 63);

    auto g_xxxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 64);

    auto g_xxxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 65);

    auto g_xxxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 67);

    auto g_xxxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 68);

    auto g_xxxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 69);

    auto g_xxxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 70);

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

    auto g_xyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 358);

    auto g_xyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 359);

    auto g_xyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 360);

    auto g_xyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 361);

    auto g_xyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 362);

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

    auto g_yyyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_hsi + 450);

    auto g_yyyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_hsi + 452);

    auto g_yyyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_hsi + 453);

    auto g_yyyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_hsi + 455);

    auto g_yyyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_hsi + 456);

    auto g_yyyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_hsi + 457);

    auto g_yyyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_hsi + 459);

    auto g_yyyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_hsi + 460);

    auto g_yyyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_hsi + 461);

    auto g_yyyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_hsi + 462);

    auto g_yyyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_hsi + 464);

    auto g_yyyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_hsi + 465);

    auto g_yyyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_hsi + 466);

    auto g_yyyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_hsi + 467);

    auto g_yyyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_hsi + 468);

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

    /// Set up components of auxilary buffer : HSK

    auto g_xxxxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk);

    auto g_xxxxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 1);

    auto g_xxxxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 2);

    auto g_xxxxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 3);

    auto g_xxxxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 4);

    auto g_xxxxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 5);

    auto g_xxxxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 6);

    auto g_xxxxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 7);

    auto g_xxxxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 8);

    auto g_xxxxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 9);

    auto g_xxxxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 10);

    auto g_xxxxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 11);

    auto g_xxxxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 12);

    auto g_xxxxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 13);

    auto g_xxxxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 14);

    auto g_xxxxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 15);

    auto g_xxxxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 16);

    auto g_xxxxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 17);

    auto g_xxxxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 18);

    auto g_xxxxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 19);

    auto g_xxxxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 20);

    auto g_xxxxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 21);

    auto g_xxxxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 22);

    auto g_xxxxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 23);

    auto g_xxxxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 24);

    auto g_xxxxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 25);

    auto g_xxxxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 26);

    auto g_xxxxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 27);

    auto g_xxxxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 28);

    auto g_xxxxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 29);

    auto g_xxxxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 30);

    auto g_xxxxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 31);

    auto g_xxxxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 32);

    auto g_xxxxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 33);

    auto g_xxxxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 34);

    auto g_xxxxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 35);

    auto g_xxxxy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 36);

    auto g_xxxxy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 37);

    auto g_xxxxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 38);

    auto g_xxxxy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 39);

    auto g_xxxxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 41);

    auto g_xxxxy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 42);

    auto g_xxxxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 45);

    auto g_xxxxy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 46);

    auto g_xxxxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 50);

    auto g_xxxxy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 51);

    auto g_xxxxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 56);

    auto g_xxxxy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 57);

    auto g_xxxxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 63);

    auto g_xxxxy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 64);

    auto g_xxxxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 72);

    auto g_xxxxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 73);

    auto g_xxxxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 74);

    auto g_xxxxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 75);

    auto g_xxxxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 76);

    auto g_xxxxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 77);

    auto g_xxxxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 78);

    auto g_xxxxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 79);

    auto g_xxxxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 80);

    auto g_xxxxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 81);

    auto g_xxxxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 82);

    auto g_xxxxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 83);

    auto g_xxxxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 84);

    auto g_xxxxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 85);

    auto g_xxxxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 86);

    auto g_xxxxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 87);

    auto g_xxxxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 88);

    auto g_xxxxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 89);

    auto g_xxxxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 90);

    auto g_xxxxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 91);

    auto g_xxxxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 92);

    auto g_xxxxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 93);

    auto g_xxxxz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 94);

    auto g_xxxxz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 95);

    auto g_xxxxz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 96);

    auto g_xxxxz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 97);

    auto g_xxxxz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 98);

    auto g_xxxxz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 99);

    auto g_xxxxz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 101);

    auto g_xxxxz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 102);

    auto g_xxxxz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 103);

    auto g_xxxxz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 104);

    auto g_xxxxz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 105);

    auto g_xxxxz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 106);

    auto g_xxxxz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 107);

    auto g_xxxyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 108);

    auto g_xxxyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 109);

    auto g_xxxyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 110);

    auto g_xxxyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 111);

    auto g_xxxyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 112);

    auto g_xxxyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 113);

    auto g_xxxyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 114);

    auto g_xxxyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 115);

    auto g_xxxyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 116);

    auto g_xxxyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 117);

    auto g_xxxyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 118);

    auto g_xxxyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 119);

    auto g_xxxyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 120);

    auto g_xxxyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 121);

    auto g_xxxyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 122);

    auto g_xxxyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 123);

    auto g_xxxyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 124);

    auto g_xxxyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 125);

    auto g_xxxyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 126);

    auto g_xxxyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 127);

    auto g_xxxyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 128);

    auto g_xxxyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 129);

    auto g_xxxyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 130);

    auto g_xxxyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 131);

    auto g_xxxyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 132);

    auto g_xxxyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 133);

    auto g_xxxyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 134);

    auto g_xxxyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 135);

    auto g_xxxyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 136);

    auto g_xxxyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 137);

    auto g_xxxyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 138);

    auto g_xxxyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 139);

    auto g_xxxyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 140);

    auto g_xxxyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 141);

    auto g_xxxyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 142);

    auto g_xxxyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 143);

    auto g_xxxzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 180);

    auto g_xxxzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 181);

    auto g_xxxzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 182);

    auto g_xxxzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 183);

    auto g_xxxzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 184);

    auto g_xxxzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 185);

    auto g_xxxzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 186);

    auto g_xxxzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 187);

    auto g_xxxzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 188);

    auto g_xxxzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 189);

    auto g_xxxzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 190);

    auto g_xxxzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 191);

    auto g_xxxzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 192);

    auto g_xxxzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 193);

    auto g_xxxzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 194);

    auto g_xxxzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 195);

    auto g_xxxzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 196);

    auto g_xxxzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 197);

    auto g_xxxzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 198);

    auto g_xxxzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 199);

    auto g_xxxzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 200);

    auto g_xxxzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 201);

    auto g_xxxzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 202);

    auto g_xxxzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 203);

    auto g_xxxzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 204);

    auto g_xxxzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 205);

    auto g_xxxzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 206);

    auto g_xxxzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 207);

    auto g_xxxzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 208);

    auto g_xxxzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 209);

    auto g_xxxzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 210);

    auto g_xxxzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 211);

    auto g_xxxzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 212);

    auto g_xxxzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 213);

    auto g_xxxzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 214);

    auto g_xxxzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 215);

    auto g_xxyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 216);

    auto g_xxyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 217);

    auto g_xxyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 218);

    auto g_xxyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 219);

    auto g_xxyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 220);

    auto g_xxyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 221);

    auto g_xxyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 222);

    auto g_xxyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 223);

    auto g_xxyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 224);

    auto g_xxyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 225);

    auto g_xxyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 226);

    auto g_xxyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 227);

    auto g_xxyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 228);

    auto g_xxyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 229);

    auto g_xxyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 230);

    auto g_xxyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 231);

    auto g_xxyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 232);

    auto g_xxyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 233);

    auto g_xxyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 234);

    auto g_xxyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 235);

    auto g_xxyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 236);

    auto g_xxyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 237);

    auto g_xxyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 238);

    auto g_xxyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 239);

    auto g_xxyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 240);

    auto g_xxyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 241);

    auto g_xxyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 242);

    auto g_xxyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 243);

    auto g_xxyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 244);

    auto g_xxyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 245);

    auto g_xxyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 246);

    auto g_xxyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 247);

    auto g_xxyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 248);

    auto g_xxyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 249);

    auto g_xxyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 250);

    auto g_xxyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 251);

    auto g_xxyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 253);

    auto g_xxyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 255);

    auto g_xxyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 258);

    auto g_xxyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 262);

    auto g_xxyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 267);

    auto g_xxyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 273);

    auto g_xxyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 288);

    auto g_xxyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 290);

    auto g_xxyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 293);

    auto g_xxyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 297);

    auto g_xxyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 302);

    auto g_xxyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 308);

    auto g_xxyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 315);

    auto g_xxzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 324);

    auto g_xxzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 325);

    auto g_xxzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 326);

    auto g_xxzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 327);

    auto g_xxzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 328);

    auto g_xxzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 329);

    auto g_xxzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 330);

    auto g_xxzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 331);

    auto g_xxzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 332);

    auto g_xxzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 333);

    auto g_xxzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 334);

    auto g_xxzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 335);

    auto g_xxzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 336);

    auto g_xxzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 337);

    auto g_xxzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 338);

    auto g_xxzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 339);

    auto g_xxzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 340);

    auto g_xxzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 341);

    auto g_xxzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 342);

    auto g_xxzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 343);

    auto g_xxzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 344);

    auto g_xxzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 345);

    auto g_xxzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 346);

    auto g_xxzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 347);

    auto g_xxzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 348);

    auto g_xxzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 349);

    auto g_xxzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 350);

    auto g_xxzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 351);

    auto g_xxzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 352);

    auto g_xxzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 353);

    auto g_xxzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 354);

    auto g_xxzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 355);

    auto g_xxzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 356);

    auto g_xxzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 357);

    auto g_xxzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 358);

    auto g_xxzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 359);

    auto g_xyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 360);

    auto g_xyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 361);

    auto g_xyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 363);

    auto g_xyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 364);

    auto g_xyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 366);

    auto g_xyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 367);

    auto g_xyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 368);

    auto g_xyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 370);

    auto g_xyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 371);

    auto g_xyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 372);

    auto g_xyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 373);

    auto g_xyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 375);

    auto g_xyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 376);

    auto g_xyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 377);

    auto g_xyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 378);

    auto g_xyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 379);

    auto g_xyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 381);

    auto g_xyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 382);

    auto g_xyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 383);

    auto g_xyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 384);

    auto g_xyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 385);

    auto g_xyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 386);

    auto g_xyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 388);

    auto g_xyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 389);

    auto g_xyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 390);

    auto g_xyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 391);

    auto g_xyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 392);

    auto g_xyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 393);

    auto g_xyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 394);

    auto g_xyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 395);

    auto g_xyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 436);

    auto g_xyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 439);

    auto g_xyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 440);

    auto g_xyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 443);

    auto g_xyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 444);

    auto g_xyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 445);

    auto g_xyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 448);

    auto g_xyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 449);

    auto g_xyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 450);

    auto g_xyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 451);

    auto g_xyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 454);

    auto g_xyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 455);

    auto g_xyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 456);

    auto g_xyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 457);

    auto g_xyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 458);

    auto g_xyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 460);

    auto g_xyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 461);

    auto g_xyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 462);

    auto g_xyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 463);

    auto g_xyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 464);

    auto g_xyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 465);

    auto g_xyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 466);

    auto g_xyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 467);

    auto g_xzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 504);

    auto g_xzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 506);

    auto g_xzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 508);

    auto g_xzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 509);

    auto g_xzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 511);

    auto g_xzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 512);

    auto g_xzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 513);

    auto g_xzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 515);

    auto g_xzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 516);

    auto g_xzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 517);

    auto g_xzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 518);

    auto g_xzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 520);

    auto g_xzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 521);

    auto g_xzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 522);

    auto g_xzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 523);

    auto g_xzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 524);

    auto g_xzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 526);

    auto g_xzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 527);

    auto g_xzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 528);

    auto g_xzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 529);

    auto g_xzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 530);

    auto g_xzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 531);

    auto g_xzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 532);

    auto g_xzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 533);

    auto g_xzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 534);

    auto g_xzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 535);

    auto g_xzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 536);

    auto g_xzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 537);

    auto g_xzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 538);

    auto g_xzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 539);

    auto g_yyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 540);

    auto g_yyyyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 541);

    auto g_yyyyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 542);

    auto g_yyyyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 543);

    auto g_yyyyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 544);

    auto g_yyyyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 545);

    auto g_yyyyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 546);

    auto g_yyyyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 547);

    auto g_yyyyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 548);

    auto g_yyyyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 549);

    auto g_yyyyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 550);

    auto g_yyyyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 551);

    auto g_yyyyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 552);

    auto g_yyyyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 553);

    auto g_yyyyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 554);

    auto g_yyyyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 555);

    auto g_yyyyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 556);

    auto g_yyyyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 557);

    auto g_yyyyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 558);

    auto g_yyyyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 559);

    auto g_yyyyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 560);

    auto g_yyyyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 561);

    auto g_yyyyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 562);

    auto g_yyyyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 563);

    auto g_yyyyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 564);

    auto g_yyyyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 565);

    auto g_yyyyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 566);

    auto g_yyyyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 567);

    auto g_yyyyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 568);

    auto g_yyyyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 569);

    auto g_yyyyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 570);

    auto g_yyyyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 571);

    auto g_yyyyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 572);

    auto g_yyyyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 573);

    auto g_yyyyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 574);

    auto g_yyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 575);

    auto g_yyyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 577);

    auto g_yyyyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 578);

    auto g_yyyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 579);

    auto g_yyyyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 580);

    auto g_yyyyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 581);

    auto g_yyyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 582);

    auto g_yyyyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 583);

    auto g_yyyyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 584);

    auto g_yyyyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 585);

    auto g_yyyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 586);

    auto g_yyyyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 587);

    auto g_yyyyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 588);

    auto g_yyyyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 589);

    auto g_yyyyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 590);

    auto g_yyyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 591);

    auto g_yyyyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 592);

    auto g_yyyyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 593);

    auto g_yyyyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 594);

    auto g_yyyyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 595);

    auto g_yyyyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 596);

    auto g_yyyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 597);

    auto g_yyyyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 598);

    auto g_yyyyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 599);

    auto g_yyyyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 600);

    auto g_yyyyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 601);

    auto g_yyyyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 602);

    auto g_yyyyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 603);

    auto g_yyyyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 604);

    auto g_yyyyz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 605);

    auto g_yyyyz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 606);

    auto g_yyyyz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 607);

    auto g_yyyyz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 608);

    auto g_yyyyz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 609);

    auto g_yyyyz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 610);

    auto g_yyyyz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 611);

    auto g_yyyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 612);

    auto g_yyyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 613);

    auto g_yyyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 614);

    auto g_yyyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 615);

    auto g_yyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 616);

    auto g_yyyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 617);

    auto g_yyyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 618);

    auto g_yyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 619);

    auto g_yyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 620);

    auto g_yyyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 621);

    auto g_yyyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 622);

    auto g_yyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 623);

    auto g_yyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 624);

    auto g_yyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 625);

    auto g_yyyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 626);

    auto g_yyyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 627);

    auto g_yyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 628);

    auto g_yyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 629);

    auto g_yyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 630);

    auto g_yyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 631);

    auto g_yyyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 632);

    auto g_yyyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 633);

    auto g_yyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 634);

    auto g_yyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 635);

    auto g_yyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 636);

    auto g_yyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 637);

    auto g_yyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 638);

    auto g_yyyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 639);

    auto g_yyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 640);

    auto g_yyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 641);

    auto g_yyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 642);

    auto g_yyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 643);

    auto g_yyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 644);

    auto g_yyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 645);

    auto g_yyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 646);

    auto g_yyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 647);

    auto g_yyzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 648);

    auto g_yyzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 649);

    auto g_yyzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 650);

    auto g_yyzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 651);

    auto g_yyzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 652);

    auto g_yyzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 653);

    auto g_yyzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 654);

    auto g_yyzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 655);

    auto g_yyzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 656);

    auto g_yyzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 657);

    auto g_yyzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 658);

    auto g_yyzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 659);

    auto g_yyzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 660);

    auto g_yyzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 661);

    auto g_yyzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 662);

    auto g_yyzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 663);

    auto g_yyzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 664);

    auto g_yyzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 665);

    auto g_yyzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 666);

    auto g_yyzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 667);

    auto g_yyzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 668);

    auto g_yyzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 669);

    auto g_yyzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 670);

    auto g_yyzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 671);

    auto g_yyzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 672);

    auto g_yyzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 673);

    auto g_yyzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 674);

    auto g_yyzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 675);

    auto g_yyzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 676);

    auto g_yyzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 677);

    auto g_yyzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 678);

    auto g_yyzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 679);

    auto g_yyzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 680);

    auto g_yyzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 681);

    auto g_yyzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 682);

    auto g_yyzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 683);

    auto g_yzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 684);

    auto g_yzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 685);

    auto g_yzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 686);

    auto g_yzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 687);

    auto g_yzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 688);

    auto g_yzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 689);

    auto g_yzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 690);

    auto g_yzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 691);

    auto g_yzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 692);

    auto g_yzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 693);

    auto g_yzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 694);

    auto g_yzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 695);

    auto g_yzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 696);

    auto g_yzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 697);

    auto g_yzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 698);

    auto g_yzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 699);

    auto g_yzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 700);

    auto g_yzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 701);

    auto g_yzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 702);

    auto g_yzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 703);

    auto g_yzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 704);

    auto g_yzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 705);

    auto g_yzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 706);

    auto g_yzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 707);

    auto g_yzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 708);

    auto g_yzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 709);

    auto g_yzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 710);

    auto g_yzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 711);

    auto g_yzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 712);

    auto g_yzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 713);

    auto g_yzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 714);

    auto g_yzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 715);

    auto g_yzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 716);

    auto g_yzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 717);

    auto g_yzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 718);

    auto g_yzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 719);

    auto g_zzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 720);

    auto g_zzzzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 721);

    auto g_zzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 722);

    auto g_zzzzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 723);

    auto g_zzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 724);

    auto g_zzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 725);

    auto g_zzzzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 726);

    auto g_zzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 727);

    auto g_zzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 728);

    auto g_zzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 729);

    auto g_zzzzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 730);

    auto g_zzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 731);

    auto g_zzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 732);

    auto g_zzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 733);

    auto g_zzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 734);

    auto g_zzzzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 735);

    auto g_zzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 736);

    auto g_zzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 737);

    auto g_zzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 738);

    auto g_zzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 739);

    auto g_zzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 740);

    auto g_zzzzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 741);

    auto g_zzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 742);

    auto g_zzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 743);

    auto g_zzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 744);

    auto g_zzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 745);

    auto g_zzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 746);

    auto g_zzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 747);

    auto g_zzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 748);

    auto g_zzzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 749);

    auto g_zzzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 750);

    auto g_zzzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 751);

    auto g_zzzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 752);

    auto g_zzzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 753);

    auto g_zzzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 754);

    auto g_zzzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 755);

    /// Set up 0-36 components of targeted buffer : ISK

    auto g_xxxxxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk);

    auto g_xxxxxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 1);

    auto g_xxxxxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 2);

    auto g_xxxxxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 3);

    auto g_xxxxxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 4);

    auto g_xxxxxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 5);

    auto g_xxxxxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 6);

    auto g_xxxxxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 7);

    auto g_xxxxxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 8);

    auto g_xxxxxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 9);

    auto g_xxxxxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 10);

    auto g_xxxxxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 11);

    auto g_xxxxxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 12);

    auto g_xxxxxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 13);

    auto g_xxxxxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 14);

    auto g_xxxxxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 15);

    auto g_xxxxxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 16);

    auto g_xxxxxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 17);

    auto g_xxxxxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 18);

    auto g_xxxxxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 19);

    auto g_xxxxxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 20);

    auto g_xxxxxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 21);

    auto g_xxxxxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 22);

    auto g_xxxxxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 23);

    auto g_xxxxxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 24);

    auto g_xxxxxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 25);

    auto g_xxxxxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 26);

    auto g_xxxxxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 27);

    auto g_xxxxxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 28);

    auto g_xxxxxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 29);

    auto g_xxxxxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 30);

    auto g_xxxxxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 31);

    auto g_xxxxxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 32);

    auto g_xxxxxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 33);

    auto g_xxxxxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 34);

    auto g_xxxxxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 35);

    #pragma omp simd aligned(g_xxxx_0_xxxxxxx_0, g_xxxx_0_xxxxxxx_1, g_xxxx_0_xxxxxxy_0, g_xxxx_0_xxxxxxy_1, g_xxxx_0_xxxxxxz_0, g_xxxx_0_xxxxxxz_1, g_xxxx_0_xxxxxyy_0, g_xxxx_0_xxxxxyy_1, g_xxxx_0_xxxxxyz_0, g_xxxx_0_xxxxxyz_1, g_xxxx_0_xxxxxzz_0, g_xxxx_0_xxxxxzz_1, g_xxxx_0_xxxxyyy_0, g_xxxx_0_xxxxyyy_1, g_xxxx_0_xxxxyyz_0, g_xxxx_0_xxxxyyz_1, g_xxxx_0_xxxxyzz_0, g_xxxx_0_xxxxyzz_1, g_xxxx_0_xxxxzzz_0, g_xxxx_0_xxxxzzz_1, g_xxxx_0_xxxyyyy_0, g_xxxx_0_xxxyyyy_1, g_xxxx_0_xxxyyyz_0, g_xxxx_0_xxxyyyz_1, g_xxxx_0_xxxyyzz_0, g_xxxx_0_xxxyyzz_1, g_xxxx_0_xxxyzzz_0, g_xxxx_0_xxxyzzz_1, g_xxxx_0_xxxzzzz_0, g_xxxx_0_xxxzzzz_1, g_xxxx_0_xxyyyyy_0, g_xxxx_0_xxyyyyy_1, g_xxxx_0_xxyyyyz_0, g_xxxx_0_xxyyyyz_1, g_xxxx_0_xxyyyzz_0, g_xxxx_0_xxyyyzz_1, g_xxxx_0_xxyyzzz_0, g_xxxx_0_xxyyzzz_1, g_xxxx_0_xxyzzzz_0, g_xxxx_0_xxyzzzz_1, g_xxxx_0_xxzzzzz_0, g_xxxx_0_xxzzzzz_1, g_xxxx_0_xyyyyyy_0, g_xxxx_0_xyyyyyy_1, g_xxxx_0_xyyyyyz_0, g_xxxx_0_xyyyyyz_1, g_xxxx_0_xyyyyzz_0, g_xxxx_0_xyyyyzz_1, g_xxxx_0_xyyyzzz_0, g_xxxx_0_xyyyzzz_1, g_xxxx_0_xyyzzzz_0, g_xxxx_0_xyyzzzz_1, g_xxxx_0_xyzzzzz_0, g_xxxx_0_xyzzzzz_1, g_xxxx_0_xzzzzzz_0, g_xxxx_0_xzzzzzz_1, g_xxxx_0_yyyyyyy_0, g_xxxx_0_yyyyyyy_1, g_xxxx_0_yyyyyyz_0, g_xxxx_0_yyyyyyz_1, g_xxxx_0_yyyyyzz_0, g_xxxx_0_yyyyyzz_1, g_xxxx_0_yyyyzzz_0, g_xxxx_0_yyyyzzz_1, g_xxxx_0_yyyzzzz_0, g_xxxx_0_yyyzzzz_1, g_xxxx_0_yyzzzzz_0, g_xxxx_0_yyzzzzz_1, g_xxxx_0_yzzzzzz_0, g_xxxx_0_yzzzzzz_1, g_xxxx_0_zzzzzzz_0, g_xxxx_0_zzzzzzz_1, g_xxxxx_0_xxxxxx_1, g_xxxxx_0_xxxxxxx_1, g_xxxxx_0_xxxxxxy_1, g_xxxxx_0_xxxxxxz_1, g_xxxxx_0_xxxxxy_1, g_xxxxx_0_xxxxxyy_1, g_xxxxx_0_xxxxxyz_1, g_xxxxx_0_xxxxxz_1, g_xxxxx_0_xxxxxzz_1, g_xxxxx_0_xxxxyy_1, g_xxxxx_0_xxxxyyy_1, g_xxxxx_0_xxxxyyz_1, g_xxxxx_0_xxxxyz_1, g_xxxxx_0_xxxxyzz_1, g_xxxxx_0_xxxxzz_1, g_xxxxx_0_xxxxzzz_1, g_xxxxx_0_xxxyyy_1, g_xxxxx_0_xxxyyyy_1, g_xxxxx_0_xxxyyyz_1, g_xxxxx_0_xxxyyz_1, g_xxxxx_0_xxxyyzz_1, g_xxxxx_0_xxxyzz_1, g_xxxxx_0_xxxyzzz_1, g_xxxxx_0_xxxzzz_1, g_xxxxx_0_xxxzzzz_1, g_xxxxx_0_xxyyyy_1, g_xxxxx_0_xxyyyyy_1, g_xxxxx_0_xxyyyyz_1, g_xxxxx_0_xxyyyz_1, g_xxxxx_0_xxyyyzz_1, g_xxxxx_0_xxyyzz_1, g_xxxxx_0_xxyyzzz_1, g_xxxxx_0_xxyzzz_1, g_xxxxx_0_xxyzzzz_1, g_xxxxx_0_xxzzzz_1, g_xxxxx_0_xxzzzzz_1, g_xxxxx_0_xyyyyy_1, g_xxxxx_0_xyyyyyy_1, g_xxxxx_0_xyyyyyz_1, g_xxxxx_0_xyyyyz_1, g_xxxxx_0_xyyyyzz_1, g_xxxxx_0_xyyyzz_1, g_xxxxx_0_xyyyzzz_1, g_xxxxx_0_xyyzzz_1, g_xxxxx_0_xyyzzzz_1, g_xxxxx_0_xyzzzz_1, g_xxxxx_0_xyzzzzz_1, g_xxxxx_0_xzzzzz_1, g_xxxxx_0_xzzzzzz_1, g_xxxxx_0_yyyyyy_1, g_xxxxx_0_yyyyyyy_1, g_xxxxx_0_yyyyyyz_1, g_xxxxx_0_yyyyyz_1, g_xxxxx_0_yyyyyzz_1, g_xxxxx_0_yyyyzz_1, g_xxxxx_0_yyyyzzz_1, g_xxxxx_0_yyyzzz_1, g_xxxxx_0_yyyzzzz_1, g_xxxxx_0_yyzzzz_1, g_xxxxx_0_yyzzzzz_1, g_xxxxx_0_yzzzzz_1, g_xxxxx_0_yzzzzzz_1, g_xxxxx_0_zzzzzz_1, g_xxxxx_0_zzzzzzz_1, g_xxxxxx_0_xxxxxxx_0, g_xxxxxx_0_xxxxxxy_0, g_xxxxxx_0_xxxxxxz_0, g_xxxxxx_0_xxxxxyy_0, g_xxxxxx_0_xxxxxyz_0, g_xxxxxx_0_xxxxxzz_0, g_xxxxxx_0_xxxxyyy_0, g_xxxxxx_0_xxxxyyz_0, g_xxxxxx_0_xxxxyzz_0, g_xxxxxx_0_xxxxzzz_0, g_xxxxxx_0_xxxyyyy_0, g_xxxxxx_0_xxxyyyz_0, g_xxxxxx_0_xxxyyzz_0, g_xxxxxx_0_xxxyzzz_0, g_xxxxxx_0_xxxzzzz_0, g_xxxxxx_0_xxyyyyy_0, g_xxxxxx_0_xxyyyyz_0, g_xxxxxx_0_xxyyyzz_0, g_xxxxxx_0_xxyyzzz_0, g_xxxxxx_0_xxyzzzz_0, g_xxxxxx_0_xxzzzzz_0, g_xxxxxx_0_xyyyyyy_0, g_xxxxxx_0_xyyyyyz_0, g_xxxxxx_0_xyyyyzz_0, g_xxxxxx_0_xyyyzzz_0, g_xxxxxx_0_xyyzzzz_0, g_xxxxxx_0_xyzzzzz_0, g_xxxxxx_0_xzzzzzz_0, g_xxxxxx_0_yyyyyyy_0, g_xxxxxx_0_yyyyyyz_0, g_xxxxxx_0_yyyyyzz_0, g_xxxxxx_0_yyyyzzz_0, g_xxxxxx_0_yyyzzzz_0, g_xxxxxx_0_yyzzzzz_0, g_xxxxxx_0_yzzzzzz_0, g_xxxxxx_0_zzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxx_0_xxxxxxx_0[i] = 5.0 * g_xxxx_0_xxxxxxx_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxxx_1[i] * fz_be_0 + 7.0 * g_xxxxx_0_xxxxxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxx_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxxy_0[i] = 5.0 * g_xxxx_0_xxxxxxy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xxxxx_0_xxxxxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxxz_0[i] = 5.0 * g_xxxx_0_xxxxxxz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xxxxx_0_xxxxxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxyy_0[i] = 5.0 * g_xxxx_0_xxxxxyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xxxxx_0_xxxxyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxyz_0[i] = 5.0 * g_xxxx_0_xxxxxyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxxxx_0_xxxxyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxxzz_0[i] = 5.0 * g_xxxx_0_xxxxxzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xxxxx_0_xxxxzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxyyy_0[i] = 5.0 * g_xxxx_0_xxxxyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxyyz_0[i] = 5.0 * g_xxxx_0_xxxxyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxyzz_0[i] = 5.0 * g_xxxx_0_xxxxyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxxzzz_0[i] = 5.0 * g_xxxx_0_xxxxzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xxxxx_0_xxxzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyyyy_0[i] = 5.0 * g_xxxx_0_xxxyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyyyz_0[i] = 5.0 * g_xxxx_0_xxxyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyyzz_0[i] = 5.0 * g_xxxx_0_xxxyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxyzzz_0[i] = 5.0 * g_xxxx_0_xxxyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxxzzzz_0[i] = 5.0 * g_xxxx_0_xxxzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xxxxx_0_xxzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyyyy_0[i] = 5.0 * g_xxxx_0_xxyyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyyyz_0[i] = 5.0 * g_xxxx_0_xxyyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyyzz_0[i] = 5.0 * g_xxxx_0_xxyyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyyzzz_0[i] = 5.0 * g_xxxx_0_xxyyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxyzzzz_0[i] = 5.0 * g_xxxx_0_xxyzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xxzzzzz_0[i] = 5.0 * g_xxxx_0_xxzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxx_0_xzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyyyy_0[i] = 5.0 * g_xxxx_0_xyyyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyyyy_1[i] * fz_be_0 + g_xxxxx_0_yyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyyyz_0[i] = 5.0 * g_xxxx_0_xyyyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyyyz_1[i] * fz_be_0 + g_xxxxx_0_yyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyyzz_0[i] = 5.0 * g_xxxx_0_xyyyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyyzz_1[i] * fz_be_0 + g_xxxxx_0_yyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyyzzz_0[i] = 5.0 * g_xxxx_0_xyyyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyyzzz_1[i] * fz_be_0 + g_xxxxx_0_yyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyyzzzz_0[i] = 5.0 * g_xxxx_0_xyyzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyyzzzz_1[i] * fz_be_0 + g_xxxxx_0_yyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xyzzzzz_0[i] = 5.0 * g_xxxx_0_xyzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xyzzzzz_1[i] * fz_be_0 + g_xxxxx_0_yzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_xzzzzzz_0[i] = 5.0 * g_xxxx_0_xzzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_xzzzzzz_1[i] * fz_be_0 + g_xxxxx_0_zzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xzzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyyyy_0[i] = 5.0 * g_xxxx_0_yyyyyyy_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyyyy_1[i] * fz_be_0 + g_xxxxx_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyyyz_0[i] = 5.0 * g_xxxx_0_yyyyyyz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyyyz_1[i] * fz_be_0 + g_xxxxx_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyyzz_0[i] = 5.0 * g_xxxx_0_yyyyyzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyyzz_1[i] * fz_be_0 + g_xxxxx_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyyzzz_0[i] = 5.0 * g_xxxx_0_yyyyzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyyzzz_1[i] * fz_be_0 + g_xxxxx_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyyzzzz_0[i] = 5.0 * g_xxxx_0_yyyzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyyzzzz_1[i] * fz_be_0 + g_xxxxx_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yyzzzzz_0[i] = 5.0 * g_xxxx_0_yyzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yyzzzzz_1[i] * fz_be_0 + g_xxxxx_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_yzzzzzz_0[i] = 5.0 * g_xxxx_0_yzzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_yzzzzzz_1[i] * fz_be_0 + g_xxxxx_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxxxx_0_zzzzzzz_0[i] = 5.0 * g_xxxx_0_zzzzzzz_0[i] * fbe_0 - 5.0 * g_xxxx_0_zzzzzzz_1[i] * fz_be_0 + g_xxxxx_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 36-72 components of targeted buffer : ISK

    auto g_xxxxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 36);

    auto g_xxxxxy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 37);

    auto g_xxxxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 38);

    auto g_xxxxxy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 39);

    auto g_xxxxxy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 40);

    auto g_xxxxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 41);

    auto g_xxxxxy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 42);

    auto g_xxxxxy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 43);

    auto g_xxxxxy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 44);

    auto g_xxxxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 45);

    auto g_xxxxxy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 46);

    auto g_xxxxxy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 47);

    auto g_xxxxxy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 48);

    auto g_xxxxxy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 49);

    auto g_xxxxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 50);

    auto g_xxxxxy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 51);

    auto g_xxxxxy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 52);

    auto g_xxxxxy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 53);

    auto g_xxxxxy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 54);

    auto g_xxxxxy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 55);

    auto g_xxxxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 56);

    auto g_xxxxxy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 57);

    auto g_xxxxxy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 58);

    auto g_xxxxxy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 59);

    auto g_xxxxxy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 60);

    auto g_xxxxxy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 61);

    auto g_xxxxxy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 62);

    auto g_xxxxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 63);

    auto g_xxxxxy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 64);

    auto g_xxxxxy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 65);

    auto g_xxxxxy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 66);

    auto g_xxxxxy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 67);

    auto g_xxxxxy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 68);

    auto g_xxxxxy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 69);

    auto g_xxxxxy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 70);

    auto g_xxxxxy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 71);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxx_1, g_xxxxx_0_xxxxxxx_1, g_xxxxx_0_xxxxxxy_1, g_xxxxx_0_xxxxxxz_1, g_xxxxx_0_xxxxxy_1, g_xxxxx_0_xxxxxyy_1, g_xxxxx_0_xxxxxyz_1, g_xxxxx_0_xxxxxz_1, g_xxxxx_0_xxxxxzz_1, g_xxxxx_0_xxxxyy_1, g_xxxxx_0_xxxxyyy_1, g_xxxxx_0_xxxxyyz_1, g_xxxxx_0_xxxxyz_1, g_xxxxx_0_xxxxyzz_1, g_xxxxx_0_xxxxzz_1, g_xxxxx_0_xxxxzzz_1, g_xxxxx_0_xxxyyy_1, g_xxxxx_0_xxxyyyy_1, g_xxxxx_0_xxxyyyz_1, g_xxxxx_0_xxxyyz_1, g_xxxxx_0_xxxyyzz_1, g_xxxxx_0_xxxyzz_1, g_xxxxx_0_xxxyzzz_1, g_xxxxx_0_xxxzzz_1, g_xxxxx_0_xxxzzzz_1, g_xxxxx_0_xxyyyy_1, g_xxxxx_0_xxyyyyy_1, g_xxxxx_0_xxyyyyz_1, g_xxxxx_0_xxyyyz_1, g_xxxxx_0_xxyyyzz_1, g_xxxxx_0_xxyyzz_1, g_xxxxx_0_xxyyzzz_1, g_xxxxx_0_xxyzzz_1, g_xxxxx_0_xxyzzzz_1, g_xxxxx_0_xxzzzz_1, g_xxxxx_0_xxzzzzz_1, g_xxxxx_0_xyyyyy_1, g_xxxxx_0_xyyyyyy_1, g_xxxxx_0_xyyyyyz_1, g_xxxxx_0_xyyyyz_1, g_xxxxx_0_xyyyyzz_1, g_xxxxx_0_xyyyzz_1, g_xxxxx_0_xyyyzzz_1, g_xxxxx_0_xyyzzz_1, g_xxxxx_0_xyyzzzz_1, g_xxxxx_0_xyzzzz_1, g_xxxxx_0_xyzzzzz_1, g_xxxxx_0_xzzzzz_1, g_xxxxx_0_xzzzzzz_1, g_xxxxx_0_yyyyyy_1, g_xxxxx_0_yyyyyyy_1, g_xxxxx_0_yyyyyyz_1, g_xxxxx_0_yyyyyz_1, g_xxxxx_0_yyyyyzz_1, g_xxxxx_0_yyyyzz_1, g_xxxxx_0_yyyyzzz_1, g_xxxxx_0_yyyzzz_1, g_xxxxx_0_yyyzzzz_1, g_xxxxx_0_yyzzzz_1, g_xxxxx_0_yyzzzzz_1, g_xxxxx_0_yzzzzz_1, g_xxxxx_0_yzzzzzz_1, g_xxxxx_0_zzzzzz_1, g_xxxxx_0_zzzzzzz_1, g_xxxxxy_0_xxxxxxx_0, g_xxxxxy_0_xxxxxxy_0, g_xxxxxy_0_xxxxxxz_0, g_xxxxxy_0_xxxxxyy_0, g_xxxxxy_0_xxxxxyz_0, g_xxxxxy_0_xxxxxzz_0, g_xxxxxy_0_xxxxyyy_0, g_xxxxxy_0_xxxxyyz_0, g_xxxxxy_0_xxxxyzz_0, g_xxxxxy_0_xxxxzzz_0, g_xxxxxy_0_xxxyyyy_0, g_xxxxxy_0_xxxyyyz_0, g_xxxxxy_0_xxxyyzz_0, g_xxxxxy_0_xxxyzzz_0, g_xxxxxy_0_xxxzzzz_0, g_xxxxxy_0_xxyyyyy_0, g_xxxxxy_0_xxyyyyz_0, g_xxxxxy_0_xxyyyzz_0, g_xxxxxy_0_xxyyzzz_0, g_xxxxxy_0_xxyzzzz_0, g_xxxxxy_0_xxzzzzz_0, g_xxxxxy_0_xyyyyyy_0, g_xxxxxy_0_xyyyyyz_0, g_xxxxxy_0_xyyyyzz_0, g_xxxxxy_0_xyyyzzz_0, g_xxxxxy_0_xyyzzzz_0, g_xxxxxy_0_xyzzzzz_0, g_xxxxxy_0_xzzzzzz_0, g_xxxxxy_0_yyyyyyy_0, g_xxxxxy_0_yyyyyyz_0, g_xxxxxy_0_yyyyyzz_0, g_xxxxxy_0_yyyyzzz_0, g_xxxxxy_0_yyyzzzz_0, g_xxxxxy_0_yyzzzzz_0, g_xxxxxy_0_yzzzzzz_0, g_xxxxxy_0_zzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxy_0_xxxxxxx_0[i] = g_xxxxx_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxxy_0[i] = g_xxxxx_0_xxxxxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxxz_0[i] = g_xxxxx_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxyy_0[i] = 2.0 * g_xxxxx_0_xxxxxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxyz_0[i] = g_xxxxx_0_xxxxxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxxzz_0[i] = g_xxxxx_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxyyy_0[i] = 3.0 * g_xxxxx_0_xxxxyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxyyz_0[i] = 2.0 * g_xxxxx_0_xxxxyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxyzz_0[i] = g_xxxxx_0_xxxxzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxxzzz_0[i] = g_xxxxx_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyyyy_0[i] = 4.0 * g_xxxxx_0_xxxyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyyyz_0[i] = 3.0 * g_xxxxx_0_xxxyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyyzz_0[i] = 2.0 * g_xxxxx_0_xxxyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxyzzz_0[i] = g_xxxxx_0_xxxzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxxzzzz_0[i] = g_xxxxx_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyyyy_0[i] = 5.0 * g_xxxxx_0_xxyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyyyz_0[i] = 4.0 * g_xxxxx_0_xxyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyyzz_0[i] = 3.0 * g_xxxxx_0_xxyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyyzzz_0[i] = 2.0 * g_xxxxx_0_xxyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxyzzzz_0[i] = g_xxxxx_0_xxzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xxzzzzz_0[i] = g_xxxxx_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyyyy_0[i] = 6.0 * g_xxxxx_0_xyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyyyz_0[i] = 5.0 * g_xxxxx_0_xyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyyzz_0[i] = 4.0 * g_xxxxx_0_xyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyyzzz_0[i] = 3.0 * g_xxxxx_0_xyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyyzzzz_0[i] = 2.0 * g_xxxxx_0_xyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xyzzzzz_0[i] = g_xxxxx_0_xzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_xzzzzzz_0[i] = g_xxxxx_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyyyy_0[i] = 7.0 * g_xxxxx_0_yyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyyy_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyyyz_0[i] = 6.0 * g_xxxxx_0_yyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyyz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyyzz_0[i] = 5.0 * g_xxxxx_0_yyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyyzzz_0[i] = 4.0 * g_xxxxx_0_yyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyyzzzz_0[i] = 3.0 * g_xxxxx_0_yyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yyzzzzz_0[i] = 2.0 * g_xxxxx_0_yzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_yzzzzzz_0[i] = g_xxxxx_0_zzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yzzzzzz_1[i] * wa_y[i];

        g_xxxxxy_0_zzzzzzz_0[i] = g_xxxxx_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 72-108 components of targeted buffer : ISK

    auto g_xxxxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 72);

    auto g_xxxxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 73);

    auto g_xxxxxz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 74);

    auto g_xxxxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 75);

    auto g_xxxxxz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 76);

    auto g_xxxxxz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 77);

    auto g_xxxxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 78);

    auto g_xxxxxz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 79);

    auto g_xxxxxz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 80);

    auto g_xxxxxz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 81);

    auto g_xxxxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 82);

    auto g_xxxxxz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 83);

    auto g_xxxxxz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 84);

    auto g_xxxxxz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 85);

    auto g_xxxxxz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 86);

    auto g_xxxxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 87);

    auto g_xxxxxz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 88);

    auto g_xxxxxz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 89);

    auto g_xxxxxz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 90);

    auto g_xxxxxz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 91);

    auto g_xxxxxz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 92);

    auto g_xxxxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 93);

    auto g_xxxxxz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 94);

    auto g_xxxxxz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 95);

    auto g_xxxxxz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 96);

    auto g_xxxxxz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 97);

    auto g_xxxxxz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 98);

    auto g_xxxxxz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 99);

    auto g_xxxxxz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 100);

    auto g_xxxxxz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 101);

    auto g_xxxxxz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 102);

    auto g_xxxxxz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 103);

    auto g_xxxxxz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 104);

    auto g_xxxxxz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 105);

    auto g_xxxxxz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 106);

    auto g_xxxxxz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 107);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxx_1, g_xxxxx_0_xxxxxxx_1, g_xxxxx_0_xxxxxxy_1, g_xxxxx_0_xxxxxxz_1, g_xxxxx_0_xxxxxy_1, g_xxxxx_0_xxxxxyy_1, g_xxxxx_0_xxxxxyz_1, g_xxxxx_0_xxxxxz_1, g_xxxxx_0_xxxxxzz_1, g_xxxxx_0_xxxxyy_1, g_xxxxx_0_xxxxyyy_1, g_xxxxx_0_xxxxyyz_1, g_xxxxx_0_xxxxyz_1, g_xxxxx_0_xxxxyzz_1, g_xxxxx_0_xxxxzz_1, g_xxxxx_0_xxxxzzz_1, g_xxxxx_0_xxxyyy_1, g_xxxxx_0_xxxyyyy_1, g_xxxxx_0_xxxyyyz_1, g_xxxxx_0_xxxyyz_1, g_xxxxx_0_xxxyyzz_1, g_xxxxx_0_xxxyzz_1, g_xxxxx_0_xxxyzzz_1, g_xxxxx_0_xxxzzz_1, g_xxxxx_0_xxxzzzz_1, g_xxxxx_0_xxyyyy_1, g_xxxxx_0_xxyyyyy_1, g_xxxxx_0_xxyyyyz_1, g_xxxxx_0_xxyyyz_1, g_xxxxx_0_xxyyyzz_1, g_xxxxx_0_xxyyzz_1, g_xxxxx_0_xxyyzzz_1, g_xxxxx_0_xxyzzz_1, g_xxxxx_0_xxyzzzz_1, g_xxxxx_0_xxzzzz_1, g_xxxxx_0_xxzzzzz_1, g_xxxxx_0_xyyyyy_1, g_xxxxx_0_xyyyyyy_1, g_xxxxx_0_xyyyyyz_1, g_xxxxx_0_xyyyyz_1, g_xxxxx_0_xyyyyzz_1, g_xxxxx_0_xyyyzz_1, g_xxxxx_0_xyyyzzz_1, g_xxxxx_0_xyyzzz_1, g_xxxxx_0_xyyzzzz_1, g_xxxxx_0_xyzzzz_1, g_xxxxx_0_xyzzzzz_1, g_xxxxx_0_xzzzzz_1, g_xxxxx_0_xzzzzzz_1, g_xxxxx_0_yyyyyy_1, g_xxxxx_0_yyyyyyy_1, g_xxxxx_0_yyyyyyz_1, g_xxxxx_0_yyyyyz_1, g_xxxxx_0_yyyyyzz_1, g_xxxxx_0_yyyyzz_1, g_xxxxx_0_yyyyzzz_1, g_xxxxx_0_yyyzzz_1, g_xxxxx_0_yyyzzzz_1, g_xxxxx_0_yyzzzz_1, g_xxxxx_0_yyzzzzz_1, g_xxxxx_0_yzzzzz_1, g_xxxxx_0_yzzzzzz_1, g_xxxxx_0_zzzzzz_1, g_xxxxx_0_zzzzzzz_1, g_xxxxxz_0_xxxxxxx_0, g_xxxxxz_0_xxxxxxy_0, g_xxxxxz_0_xxxxxxz_0, g_xxxxxz_0_xxxxxyy_0, g_xxxxxz_0_xxxxxyz_0, g_xxxxxz_0_xxxxxzz_0, g_xxxxxz_0_xxxxyyy_0, g_xxxxxz_0_xxxxyyz_0, g_xxxxxz_0_xxxxyzz_0, g_xxxxxz_0_xxxxzzz_0, g_xxxxxz_0_xxxyyyy_0, g_xxxxxz_0_xxxyyyz_0, g_xxxxxz_0_xxxyyzz_0, g_xxxxxz_0_xxxyzzz_0, g_xxxxxz_0_xxxzzzz_0, g_xxxxxz_0_xxyyyyy_0, g_xxxxxz_0_xxyyyyz_0, g_xxxxxz_0_xxyyyzz_0, g_xxxxxz_0_xxyyzzz_0, g_xxxxxz_0_xxyzzzz_0, g_xxxxxz_0_xxzzzzz_0, g_xxxxxz_0_xyyyyyy_0, g_xxxxxz_0_xyyyyyz_0, g_xxxxxz_0_xyyyyzz_0, g_xxxxxz_0_xyyyzzz_0, g_xxxxxz_0_xyyzzzz_0, g_xxxxxz_0_xyzzzzz_0, g_xxxxxz_0_xzzzzzz_0, g_xxxxxz_0_yyyyyyy_0, g_xxxxxz_0_yyyyyyz_0, g_xxxxxz_0_yyyyyzz_0, g_xxxxxz_0_yyyyzzz_0, g_xxxxxz_0_yyyzzzz_0, g_xxxxxz_0_yyzzzzz_0, g_xxxxxz_0_yzzzzzz_0, g_xxxxxz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxz_0_xxxxxxx_0[i] = g_xxxxx_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxxy_0[i] = g_xxxxx_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxxz_0[i] = g_xxxxx_0_xxxxxx_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxxz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxyy_0[i] = g_xxxxx_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxyz_0[i] = g_xxxxx_0_xxxxxy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxxzz_0[i] = 2.0 * g_xxxxx_0_xxxxxz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxxzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxyyy_0[i] = g_xxxxx_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxyyz_0[i] = g_xxxxx_0_xxxxyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxyzz_0[i] = 2.0 * g_xxxxx_0_xxxxyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxxzzz_0[i] = 3.0 * g_xxxxx_0_xxxxzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxxzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyyyy_0[i] = g_xxxxx_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyyyz_0[i] = g_xxxxx_0_xxxyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyyzz_0[i] = 2.0 * g_xxxxx_0_xxxyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxyzzz_0[i] = 3.0 * g_xxxxx_0_xxxyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxxzzzz_0[i] = 4.0 * g_xxxxx_0_xxxzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxxzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyyyy_0[i] = g_xxxxx_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyyyz_0[i] = g_xxxxx_0_xxyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyyzz_0[i] = 2.0 * g_xxxxx_0_xxyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyyzzz_0[i] = 3.0 * g_xxxxx_0_xxyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxyzzzz_0[i] = 4.0 * g_xxxxx_0_xxyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxyzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xxzzzzz_0[i] = 5.0 * g_xxxxx_0_xxzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xxzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyyyy_0[i] = g_xxxxx_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyyyz_0[i] = g_xxxxx_0_xyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyyzz_0[i] = 2.0 * g_xxxxx_0_xyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyyzzz_0[i] = 3.0 * g_xxxxx_0_xyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyyzzzz_0[i] = 4.0 * g_xxxxx_0_xyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyyzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xyzzzzz_0[i] = 5.0 * g_xxxxx_0_xyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xyzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_xzzzzzz_0[i] = 6.0 * g_xxxxx_0_xzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_xzzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyyyy_0[i] = g_xxxxx_0_yyyyyyy_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyyyz_0[i] = g_xxxxx_0_yyyyyy_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyyz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyyzz_0[i] = 2.0 * g_xxxxx_0_yyyyyz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyyzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyyzzz_0[i] = 3.0 * g_xxxxx_0_yyyyzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyyzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyyzzzz_0[i] = 4.0 * g_xxxxx_0_yyyzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyyzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yyzzzzz_0[i] = 5.0 * g_xxxxx_0_yyzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yyzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_yzzzzzz_0[i] = 6.0 * g_xxxxx_0_yzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_yzzzzzz_1[i] * wa_z[i];

        g_xxxxxz_0_zzzzzzz_0[i] = 7.0 * g_xxxxx_0_zzzzzz_1[i] * fi_acd_0 + g_xxxxx_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 108-144 components of targeted buffer : ISK

    auto g_xxxxyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 108);

    auto g_xxxxyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 109);

    auto g_xxxxyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 110);

    auto g_xxxxyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 111);

    auto g_xxxxyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 112);

    auto g_xxxxyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 113);

    auto g_xxxxyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 114);

    auto g_xxxxyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 115);

    auto g_xxxxyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 116);

    auto g_xxxxyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 117);

    auto g_xxxxyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 118);

    auto g_xxxxyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 119);

    auto g_xxxxyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 120);

    auto g_xxxxyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 121);

    auto g_xxxxyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 122);

    auto g_xxxxyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 123);

    auto g_xxxxyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 124);

    auto g_xxxxyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 125);

    auto g_xxxxyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 126);

    auto g_xxxxyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 127);

    auto g_xxxxyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 128);

    auto g_xxxxyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 129);

    auto g_xxxxyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 130);

    auto g_xxxxyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 131);

    auto g_xxxxyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 132);

    auto g_xxxxyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 133);

    auto g_xxxxyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 134);

    auto g_xxxxyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 135);

    auto g_xxxxyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 136);

    auto g_xxxxyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 137);

    auto g_xxxxyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 138);

    auto g_xxxxyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 139);

    auto g_xxxxyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 140);

    auto g_xxxxyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 141);

    auto g_xxxxyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 142);

    auto g_xxxxyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 143);

    #pragma omp simd aligned(g_xxxx_0_xxxxxxx_0, g_xxxx_0_xxxxxxx_1, g_xxxx_0_xxxxxxz_0, g_xxxx_0_xxxxxxz_1, g_xxxx_0_xxxxxzz_0, g_xxxx_0_xxxxxzz_1, g_xxxx_0_xxxxzzz_0, g_xxxx_0_xxxxzzz_1, g_xxxx_0_xxxzzzz_0, g_xxxx_0_xxxzzzz_1, g_xxxx_0_xxzzzzz_0, g_xxxx_0_xxzzzzz_1, g_xxxx_0_xzzzzzz_0, g_xxxx_0_xzzzzzz_1, g_xxxxy_0_xxxxxxx_1, g_xxxxy_0_xxxxxxz_1, g_xxxxy_0_xxxxxzz_1, g_xxxxy_0_xxxxzzz_1, g_xxxxy_0_xxxzzzz_1, g_xxxxy_0_xxzzzzz_1, g_xxxxy_0_xzzzzzz_1, g_xxxxyy_0_xxxxxxx_0, g_xxxxyy_0_xxxxxxy_0, g_xxxxyy_0_xxxxxxz_0, g_xxxxyy_0_xxxxxyy_0, g_xxxxyy_0_xxxxxyz_0, g_xxxxyy_0_xxxxxzz_0, g_xxxxyy_0_xxxxyyy_0, g_xxxxyy_0_xxxxyyz_0, g_xxxxyy_0_xxxxyzz_0, g_xxxxyy_0_xxxxzzz_0, g_xxxxyy_0_xxxyyyy_0, g_xxxxyy_0_xxxyyyz_0, g_xxxxyy_0_xxxyyzz_0, g_xxxxyy_0_xxxyzzz_0, g_xxxxyy_0_xxxzzzz_0, g_xxxxyy_0_xxyyyyy_0, g_xxxxyy_0_xxyyyyz_0, g_xxxxyy_0_xxyyyzz_0, g_xxxxyy_0_xxyyzzz_0, g_xxxxyy_0_xxyzzzz_0, g_xxxxyy_0_xxzzzzz_0, g_xxxxyy_0_xyyyyyy_0, g_xxxxyy_0_xyyyyyz_0, g_xxxxyy_0_xyyyyzz_0, g_xxxxyy_0_xyyyzzz_0, g_xxxxyy_0_xyyzzzz_0, g_xxxxyy_0_xyzzzzz_0, g_xxxxyy_0_xzzzzzz_0, g_xxxxyy_0_yyyyyyy_0, g_xxxxyy_0_yyyyyyz_0, g_xxxxyy_0_yyyyyzz_0, g_xxxxyy_0_yyyyzzz_0, g_xxxxyy_0_yyyzzzz_0, g_xxxxyy_0_yyzzzzz_0, g_xxxxyy_0_yzzzzzz_0, g_xxxxyy_0_zzzzzzz_0, g_xxxyy_0_xxxxxxy_1, g_xxxyy_0_xxxxxy_1, g_xxxyy_0_xxxxxyy_1, g_xxxyy_0_xxxxxyz_1, g_xxxyy_0_xxxxyy_1, g_xxxyy_0_xxxxyyy_1, g_xxxyy_0_xxxxyyz_1, g_xxxyy_0_xxxxyz_1, g_xxxyy_0_xxxxyzz_1, g_xxxyy_0_xxxyyy_1, g_xxxyy_0_xxxyyyy_1, g_xxxyy_0_xxxyyyz_1, g_xxxyy_0_xxxyyz_1, g_xxxyy_0_xxxyyzz_1, g_xxxyy_0_xxxyzz_1, g_xxxyy_0_xxxyzzz_1, g_xxxyy_0_xxyyyy_1, g_xxxyy_0_xxyyyyy_1, g_xxxyy_0_xxyyyyz_1, g_xxxyy_0_xxyyyz_1, g_xxxyy_0_xxyyyzz_1, g_xxxyy_0_xxyyzz_1, g_xxxyy_0_xxyyzzz_1, g_xxxyy_0_xxyzzz_1, g_xxxyy_0_xxyzzzz_1, g_xxxyy_0_xyyyyy_1, g_xxxyy_0_xyyyyyy_1, g_xxxyy_0_xyyyyyz_1, g_xxxyy_0_xyyyyz_1, g_xxxyy_0_xyyyyzz_1, g_xxxyy_0_xyyyzz_1, g_xxxyy_0_xyyyzzz_1, g_xxxyy_0_xyyzzz_1, g_xxxyy_0_xyyzzzz_1, g_xxxyy_0_xyzzzz_1, g_xxxyy_0_xyzzzzz_1, g_xxxyy_0_yyyyyy_1, g_xxxyy_0_yyyyyyy_1, g_xxxyy_0_yyyyyyz_1, g_xxxyy_0_yyyyyz_1, g_xxxyy_0_yyyyyzz_1, g_xxxyy_0_yyyyzz_1, g_xxxyy_0_yyyyzzz_1, g_xxxyy_0_yyyzzz_1, g_xxxyy_0_yyyzzzz_1, g_xxxyy_0_yyzzzz_1, g_xxxyy_0_yyzzzzz_1, g_xxxyy_0_yzzzzz_1, g_xxxyy_0_yzzzzzz_1, g_xxxyy_0_zzzzzzz_1, g_xxyy_0_xxxxxxy_0, g_xxyy_0_xxxxxxy_1, g_xxyy_0_xxxxxyy_0, g_xxyy_0_xxxxxyy_1, g_xxyy_0_xxxxxyz_0, g_xxyy_0_xxxxxyz_1, g_xxyy_0_xxxxyyy_0, g_xxyy_0_xxxxyyy_1, g_xxyy_0_xxxxyyz_0, g_xxyy_0_xxxxyyz_1, g_xxyy_0_xxxxyzz_0, g_xxyy_0_xxxxyzz_1, g_xxyy_0_xxxyyyy_0, g_xxyy_0_xxxyyyy_1, g_xxyy_0_xxxyyyz_0, g_xxyy_0_xxxyyyz_1, g_xxyy_0_xxxyyzz_0, g_xxyy_0_xxxyyzz_1, g_xxyy_0_xxxyzzz_0, g_xxyy_0_xxxyzzz_1, g_xxyy_0_xxyyyyy_0, g_xxyy_0_xxyyyyy_1, g_xxyy_0_xxyyyyz_0, g_xxyy_0_xxyyyyz_1, g_xxyy_0_xxyyyzz_0, g_xxyy_0_xxyyyzz_1, g_xxyy_0_xxyyzzz_0, g_xxyy_0_xxyyzzz_1, g_xxyy_0_xxyzzzz_0, g_xxyy_0_xxyzzzz_1, g_xxyy_0_xyyyyyy_0, g_xxyy_0_xyyyyyy_1, g_xxyy_0_xyyyyyz_0, g_xxyy_0_xyyyyyz_1, g_xxyy_0_xyyyyzz_0, g_xxyy_0_xyyyyzz_1, g_xxyy_0_xyyyzzz_0, g_xxyy_0_xyyyzzz_1, g_xxyy_0_xyyzzzz_0, g_xxyy_0_xyyzzzz_1, g_xxyy_0_xyzzzzz_0, g_xxyy_0_xyzzzzz_1, g_xxyy_0_yyyyyyy_0, g_xxyy_0_yyyyyyy_1, g_xxyy_0_yyyyyyz_0, g_xxyy_0_yyyyyyz_1, g_xxyy_0_yyyyyzz_0, g_xxyy_0_yyyyyzz_1, g_xxyy_0_yyyyzzz_0, g_xxyy_0_yyyyzzz_1, g_xxyy_0_yyyzzzz_0, g_xxyy_0_yyyzzzz_1, g_xxyy_0_yyzzzzz_0, g_xxyy_0_yyzzzzz_1, g_xxyy_0_yzzzzzz_0, g_xxyy_0_yzzzzzz_1, g_xxyy_0_zzzzzzz_0, g_xxyy_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyy_0_xxxxxxx_0[i] = g_xxxx_0_xxxxxxx_0[i] * fbe_0 - g_xxxx_0_xxxxxxx_1[i] * fz_be_0 + g_xxxxy_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxxyy_0_xxxxxxy_0[i] = 3.0 * g_xxyy_0_xxxxxxy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xxxyy_0_xxxxxy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxxxz_0[i] = g_xxxx_0_xxxxxxz_0[i] * fbe_0 - g_xxxx_0_xxxxxxz_1[i] * fz_be_0 + g_xxxxy_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxxyy_0_xxxxxyy_0[i] = 3.0 * g_xxyy_0_xxxxxyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xxxyy_0_xxxxyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxxyz_0[i] = 3.0 * g_xxyy_0_xxxxxyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxxyy_0_xxxxyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxxzz_0[i] = g_xxxx_0_xxxxxzz_0[i] * fbe_0 - g_xxxx_0_xxxxxzz_1[i] * fz_be_0 + g_xxxxy_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxxyy_0_xxxxyyy_0[i] = 3.0 * g_xxyy_0_xxxxyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xxxyy_0_xxxyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxyyz_0[i] = 3.0 * g_xxyy_0_xxxxyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxxyy_0_xxxyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxyzz_0[i] = 3.0 * g_xxyy_0_xxxxyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxxyy_0_xxxyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxxzzz_0[i] = g_xxxx_0_xxxxzzz_0[i] * fbe_0 - g_xxxx_0_xxxxzzz_1[i] * fz_be_0 + g_xxxxy_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxxyy_0_xxxyyyy_0[i] = 3.0 * g_xxyy_0_xxxyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxxyyyz_0[i] = 3.0 * g_xxyy_0_xxxyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxyyzz_0[i] = 3.0 * g_xxyy_0_xxxyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxyzzz_0[i] = 3.0 * g_xxyy_0_xxxyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxxyy_0_xxyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxxzzzz_0[i] = g_xxxx_0_xxxzzzz_0[i] * fbe_0 - g_xxxx_0_xxxzzzz_1[i] * fz_be_0 + g_xxxxy_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxxyy_0_xxyyyyy_0[i] = 3.0 * g_xxyy_0_xxyyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xxyyyyz_0[i] = 3.0 * g_xxyy_0_xxyyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xxyyyzz_0[i] = 3.0 * g_xxyy_0_xxyyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxyyzzz_0[i] = 3.0 * g_xxyy_0_xxyyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxyzzzz_0[i] = 3.0 * g_xxyy_0_xxyzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxyy_0_xyzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xxzzzzz_0[i] = g_xxxx_0_xxzzzzz_0[i] * fbe_0 - g_xxxx_0_xxzzzzz_1[i] * fz_be_0 + g_xxxxy_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxxyy_0_xyyyyyy_0[i] = 3.0 * g_xxyy_0_xyyyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyyyy_1[i] * fz_be_0 + g_xxxyy_0_yyyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_xyyyyyz_0[i] = 3.0 * g_xxyy_0_xyyyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyyyz_1[i] * fz_be_0 + g_xxxyy_0_yyyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_xyyyyzz_0[i] = 3.0 * g_xxyy_0_xyyyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyyzz_1[i] * fz_be_0 + g_xxxyy_0_yyyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_xyyyzzz_0[i] = 3.0 * g_xxyy_0_xyyyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyyzzz_1[i] * fz_be_0 + g_xxxyy_0_yyyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xyyzzzz_0[i] = 3.0 * g_xxyy_0_xyyzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyyzzzz_1[i] * fz_be_0 + g_xxxyy_0_yyzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xyzzzzz_0[i] = 3.0 * g_xxyy_0_xyzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xyzzzzz_1[i] * fz_be_0 + g_xxxyy_0_yzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_xzzzzzz_0[i] = g_xxxx_0_xzzzzzz_0[i] * fbe_0 - g_xxxx_0_xzzzzzz_1[i] * fz_be_0 + g_xxxxy_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxxyy_0_yyyyyyy_0[i] = 3.0 * g_xxyy_0_yyyyyyy_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyyyy_1[i] * fz_be_0 + g_xxxyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxxyy_0_yyyyyyz_0[i] = 3.0 * g_xxyy_0_yyyyyyz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyyyz_1[i] * fz_be_0 + g_xxxyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxxyy_0_yyyyyzz_0[i] = 3.0 * g_xxyy_0_yyyyyzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyyzz_1[i] * fz_be_0 + g_xxxyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxxyy_0_yyyyzzz_0[i] = 3.0 * g_xxyy_0_yyyyzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyyzzz_1[i] * fz_be_0 + g_xxxyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxxyy_0_yyyzzzz_0[i] = 3.0 * g_xxyy_0_yyyzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyyzzzz_1[i] * fz_be_0 + g_xxxyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_yyzzzzz_0[i] = 3.0 * g_xxyy_0_yyzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yyzzzzz_1[i] * fz_be_0 + g_xxxyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_yzzzzzz_0[i] = 3.0 * g_xxyy_0_yzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_yzzzzzz_1[i] * fz_be_0 + g_xxxyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxxyy_0_zzzzzzz_0[i] = 3.0 * g_xxyy_0_zzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_zzzzzzz_1[i] * fz_be_0 + g_xxxyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 144-180 components of targeted buffer : ISK

    auto g_xxxxyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 144);

    auto g_xxxxyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 145);

    auto g_xxxxyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 146);

    auto g_xxxxyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 147);

    auto g_xxxxyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 148);

    auto g_xxxxyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 149);

    auto g_xxxxyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 150);

    auto g_xxxxyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 151);

    auto g_xxxxyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 152);

    auto g_xxxxyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 153);

    auto g_xxxxyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 154);

    auto g_xxxxyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 155);

    auto g_xxxxyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 156);

    auto g_xxxxyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 157);

    auto g_xxxxyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 158);

    auto g_xxxxyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 159);

    auto g_xxxxyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 160);

    auto g_xxxxyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 161);

    auto g_xxxxyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 162);

    auto g_xxxxyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 163);

    auto g_xxxxyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 164);

    auto g_xxxxyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 165);

    auto g_xxxxyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 166);

    auto g_xxxxyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 167);

    auto g_xxxxyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 168);

    auto g_xxxxyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 169);

    auto g_xxxxyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 170);

    auto g_xxxxyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 171);

    auto g_xxxxyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 172);

    auto g_xxxxyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 173);

    auto g_xxxxyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 174);

    auto g_xxxxyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 175);

    auto g_xxxxyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 176);

    auto g_xxxxyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 177);

    auto g_xxxxyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 178);

    auto g_xxxxyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 179);

    #pragma omp simd aligned(g_xxxxy_0_xxxxxxy_1, g_xxxxy_0_xxxxxyy_1, g_xxxxy_0_xxxxyyy_1, g_xxxxy_0_xxxyyyy_1, g_xxxxy_0_xxyyyyy_1, g_xxxxy_0_xyyyyyy_1, g_xxxxy_0_yyyyyyy_1, g_xxxxyz_0_xxxxxxx_0, g_xxxxyz_0_xxxxxxy_0, g_xxxxyz_0_xxxxxxz_0, g_xxxxyz_0_xxxxxyy_0, g_xxxxyz_0_xxxxxyz_0, g_xxxxyz_0_xxxxxzz_0, g_xxxxyz_0_xxxxyyy_0, g_xxxxyz_0_xxxxyyz_0, g_xxxxyz_0_xxxxyzz_0, g_xxxxyz_0_xxxxzzz_0, g_xxxxyz_0_xxxyyyy_0, g_xxxxyz_0_xxxyyyz_0, g_xxxxyz_0_xxxyyzz_0, g_xxxxyz_0_xxxyzzz_0, g_xxxxyz_0_xxxzzzz_0, g_xxxxyz_0_xxyyyyy_0, g_xxxxyz_0_xxyyyyz_0, g_xxxxyz_0_xxyyyzz_0, g_xxxxyz_0_xxyyzzz_0, g_xxxxyz_0_xxyzzzz_0, g_xxxxyz_0_xxzzzzz_0, g_xxxxyz_0_xyyyyyy_0, g_xxxxyz_0_xyyyyyz_0, g_xxxxyz_0_xyyyyzz_0, g_xxxxyz_0_xyyyzzz_0, g_xxxxyz_0_xyyzzzz_0, g_xxxxyz_0_xyzzzzz_0, g_xxxxyz_0_xzzzzzz_0, g_xxxxyz_0_yyyyyyy_0, g_xxxxyz_0_yyyyyyz_0, g_xxxxyz_0_yyyyyzz_0, g_xxxxyz_0_yyyyzzz_0, g_xxxxyz_0_yyyzzzz_0, g_xxxxyz_0_yyzzzzz_0, g_xxxxyz_0_yzzzzzz_0, g_xxxxyz_0_zzzzzzz_0, g_xxxxz_0_xxxxxxx_1, g_xxxxz_0_xxxxxxz_1, g_xxxxz_0_xxxxxyz_1, g_xxxxz_0_xxxxxz_1, g_xxxxz_0_xxxxxzz_1, g_xxxxz_0_xxxxyyz_1, g_xxxxz_0_xxxxyz_1, g_xxxxz_0_xxxxyzz_1, g_xxxxz_0_xxxxzz_1, g_xxxxz_0_xxxxzzz_1, g_xxxxz_0_xxxyyyz_1, g_xxxxz_0_xxxyyz_1, g_xxxxz_0_xxxyyzz_1, g_xxxxz_0_xxxyzz_1, g_xxxxz_0_xxxyzzz_1, g_xxxxz_0_xxxzzz_1, g_xxxxz_0_xxxzzzz_1, g_xxxxz_0_xxyyyyz_1, g_xxxxz_0_xxyyyz_1, g_xxxxz_0_xxyyyzz_1, g_xxxxz_0_xxyyzz_1, g_xxxxz_0_xxyyzzz_1, g_xxxxz_0_xxyzzz_1, g_xxxxz_0_xxyzzzz_1, g_xxxxz_0_xxzzzz_1, g_xxxxz_0_xxzzzzz_1, g_xxxxz_0_xyyyyyz_1, g_xxxxz_0_xyyyyz_1, g_xxxxz_0_xyyyyzz_1, g_xxxxz_0_xyyyzz_1, g_xxxxz_0_xyyyzzz_1, g_xxxxz_0_xyyzzz_1, g_xxxxz_0_xyyzzzz_1, g_xxxxz_0_xyzzzz_1, g_xxxxz_0_xyzzzzz_1, g_xxxxz_0_xzzzzz_1, g_xxxxz_0_xzzzzzz_1, g_xxxxz_0_yyyyyyz_1, g_xxxxz_0_yyyyyz_1, g_xxxxz_0_yyyyyzz_1, g_xxxxz_0_yyyyzz_1, g_xxxxz_0_yyyyzzz_1, g_xxxxz_0_yyyzzz_1, g_xxxxz_0_yyyzzzz_1, g_xxxxz_0_yyzzzz_1, g_xxxxz_0_yyzzzzz_1, g_xxxxz_0_yzzzzz_1, g_xxxxz_0_yzzzzzz_1, g_xxxxz_0_zzzzzz_1, g_xxxxz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyz_0_xxxxxxx_0[i] = g_xxxxz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxxxy_0[i] = g_xxxxy_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxxxxz_0[i] = g_xxxxz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxxyy_0[i] = g_xxxxy_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxxxyz_0[i] = g_xxxxz_0_xxxxxz_1[i] * fi_acd_0 + g_xxxxz_0_xxxxxyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxxzz_0[i] = g_xxxxz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxyyy_0[i] = g_xxxxy_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxxyyz_0[i] = 2.0 * g_xxxxz_0_xxxxyz_1[i] * fi_acd_0 + g_xxxxz_0_xxxxyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxyzz_0[i] = g_xxxxz_0_xxxxzz_1[i] * fi_acd_0 + g_xxxxz_0_xxxxyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxxzzz_0[i] = g_xxxxz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxyyyy_0[i] = g_xxxxy_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxxyyyz_0[i] = 3.0 * g_xxxxz_0_xxxyyz_1[i] * fi_acd_0 + g_xxxxz_0_xxxyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxyyzz_0[i] = 2.0 * g_xxxxz_0_xxxyzz_1[i] * fi_acd_0 + g_xxxxz_0_xxxyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxyzzz_0[i] = g_xxxxz_0_xxxzzz_1[i] * fi_acd_0 + g_xxxxz_0_xxxyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxxzzzz_0[i] = g_xxxxz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyyyyy_0[i] = g_xxxxy_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xxyyyyz_0[i] = 4.0 * g_xxxxz_0_xxyyyz_1[i] * fi_acd_0 + g_xxxxz_0_xxyyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyyyzz_0[i] = 3.0 * g_xxxxz_0_xxyyzz_1[i] * fi_acd_0 + g_xxxxz_0_xxyyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyyzzz_0[i] = 2.0 * g_xxxxz_0_xxyzzz_1[i] * fi_acd_0 + g_xxxxz_0_xxyyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxyzzzz_0[i] = g_xxxxz_0_xxzzzz_1[i] * fi_acd_0 + g_xxxxz_0_xxyzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xxzzzzz_0[i] = g_xxxxz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyyyyy_0[i] = g_xxxxy_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_xyyyyyz_0[i] = 5.0 * g_xxxxz_0_xyyyyz_1[i] * fi_acd_0 + g_xxxxz_0_xyyyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyyyzz_0[i] = 4.0 * g_xxxxz_0_xyyyzz_1[i] * fi_acd_0 + g_xxxxz_0_xyyyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyyzzz_0[i] = 3.0 * g_xxxxz_0_xyyzzz_1[i] * fi_acd_0 + g_xxxxz_0_xyyyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyyzzzz_0[i] = 2.0 * g_xxxxz_0_xyzzzz_1[i] * fi_acd_0 + g_xxxxz_0_xyyzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xyzzzzz_0[i] = g_xxxxz_0_xzzzzz_1[i] * fi_acd_0 + g_xxxxz_0_xyzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_xzzzzzz_0[i] = g_xxxxz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyyyyy_0[i] = g_xxxxy_0_yyyyyyy_1[i] * wa_z[i];

        g_xxxxyz_0_yyyyyyz_0[i] = 6.0 * g_xxxxz_0_yyyyyz_1[i] * fi_acd_0 + g_xxxxz_0_yyyyyyz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyyyzz_0[i] = 5.0 * g_xxxxz_0_yyyyzz_1[i] * fi_acd_0 + g_xxxxz_0_yyyyyzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyyzzz_0[i] = 4.0 * g_xxxxz_0_yyyzzz_1[i] * fi_acd_0 + g_xxxxz_0_yyyyzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyyzzzz_0[i] = 3.0 * g_xxxxz_0_yyzzzz_1[i] * fi_acd_0 + g_xxxxz_0_yyyzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yyzzzzz_0[i] = 2.0 * g_xxxxz_0_yzzzzz_1[i] * fi_acd_0 + g_xxxxz_0_yyzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_yzzzzzz_0[i] = g_xxxxz_0_zzzzzz_1[i] * fi_acd_0 + g_xxxxz_0_yzzzzzz_1[i] * wa_y[i];

        g_xxxxyz_0_zzzzzzz_0[i] = g_xxxxz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 180-216 components of targeted buffer : ISK

    auto g_xxxxzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 180);

    auto g_xxxxzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 181);

    auto g_xxxxzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 182);

    auto g_xxxxzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 183);

    auto g_xxxxzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 184);

    auto g_xxxxzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 185);

    auto g_xxxxzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 186);

    auto g_xxxxzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 187);

    auto g_xxxxzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 188);

    auto g_xxxxzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 189);

    auto g_xxxxzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 190);

    auto g_xxxxzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 191);

    auto g_xxxxzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 192);

    auto g_xxxxzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 193);

    auto g_xxxxzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 194);

    auto g_xxxxzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 195);

    auto g_xxxxzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 196);

    auto g_xxxxzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 197);

    auto g_xxxxzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 198);

    auto g_xxxxzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 199);

    auto g_xxxxzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 200);

    auto g_xxxxzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 201);

    auto g_xxxxzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 202);

    auto g_xxxxzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 203);

    auto g_xxxxzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 204);

    auto g_xxxxzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 205);

    auto g_xxxxzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 206);

    auto g_xxxxzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 207);

    auto g_xxxxzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 208);

    auto g_xxxxzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 209);

    auto g_xxxxzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 210);

    auto g_xxxxzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 211);

    auto g_xxxxzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 212);

    auto g_xxxxzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 213);

    auto g_xxxxzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 214);

    auto g_xxxxzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 215);

    #pragma omp simd aligned(g_xxxx_0_xxxxxxx_0, g_xxxx_0_xxxxxxx_1, g_xxxx_0_xxxxxxy_0, g_xxxx_0_xxxxxxy_1, g_xxxx_0_xxxxxyy_0, g_xxxx_0_xxxxxyy_1, g_xxxx_0_xxxxyyy_0, g_xxxx_0_xxxxyyy_1, g_xxxx_0_xxxyyyy_0, g_xxxx_0_xxxyyyy_1, g_xxxx_0_xxyyyyy_0, g_xxxx_0_xxyyyyy_1, g_xxxx_0_xyyyyyy_0, g_xxxx_0_xyyyyyy_1, g_xxxxz_0_xxxxxxx_1, g_xxxxz_0_xxxxxxy_1, g_xxxxz_0_xxxxxyy_1, g_xxxxz_0_xxxxyyy_1, g_xxxxz_0_xxxyyyy_1, g_xxxxz_0_xxyyyyy_1, g_xxxxz_0_xyyyyyy_1, g_xxxxzz_0_xxxxxxx_0, g_xxxxzz_0_xxxxxxy_0, g_xxxxzz_0_xxxxxxz_0, g_xxxxzz_0_xxxxxyy_0, g_xxxxzz_0_xxxxxyz_0, g_xxxxzz_0_xxxxxzz_0, g_xxxxzz_0_xxxxyyy_0, g_xxxxzz_0_xxxxyyz_0, g_xxxxzz_0_xxxxyzz_0, g_xxxxzz_0_xxxxzzz_0, g_xxxxzz_0_xxxyyyy_0, g_xxxxzz_0_xxxyyyz_0, g_xxxxzz_0_xxxyyzz_0, g_xxxxzz_0_xxxyzzz_0, g_xxxxzz_0_xxxzzzz_0, g_xxxxzz_0_xxyyyyy_0, g_xxxxzz_0_xxyyyyz_0, g_xxxxzz_0_xxyyyzz_0, g_xxxxzz_0_xxyyzzz_0, g_xxxxzz_0_xxyzzzz_0, g_xxxxzz_0_xxzzzzz_0, g_xxxxzz_0_xyyyyyy_0, g_xxxxzz_0_xyyyyyz_0, g_xxxxzz_0_xyyyyzz_0, g_xxxxzz_0_xyyyzzz_0, g_xxxxzz_0_xyyzzzz_0, g_xxxxzz_0_xyzzzzz_0, g_xxxxzz_0_xzzzzzz_0, g_xxxxzz_0_yyyyyyy_0, g_xxxxzz_0_yyyyyyz_0, g_xxxxzz_0_yyyyyzz_0, g_xxxxzz_0_yyyyzzz_0, g_xxxxzz_0_yyyzzzz_0, g_xxxxzz_0_yyzzzzz_0, g_xxxxzz_0_yzzzzzz_0, g_xxxxzz_0_zzzzzzz_0, g_xxxzz_0_xxxxxxz_1, g_xxxzz_0_xxxxxyz_1, g_xxxzz_0_xxxxxz_1, g_xxxzz_0_xxxxxzz_1, g_xxxzz_0_xxxxyyz_1, g_xxxzz_0_xxxxyz_1, g_xxxzz_0_xxxxyzz_1, g_xxxzz_0_xxxxzz_1, g_xxxzz_0_xxxxzzz_1, g_xxxzz_0_xxxyyyz_1, g_xxxzz_0_xxxyyz_1, g_xxxzz_0_xxxyyzz_1, g_xxxzz_0_xxxyzz_1, g_xxxzz_0_xxxyzzz_1, g_xxxzz_0_xxxzzz_1, g_xxxzz_0_xxxzzzz_1, g_xxxzz_0_xxyyyyz_1, g_xxxzz_0_xxyyyz_1, g_xxxzz_0_xxyyyzz_1, g_xxxzz_0_xxyyzz_1, g_xxxzz_0_xxyyzzz_1, g_xxxzz_0_xxyzzz_1, g_xxxzz_0_xxyzzzz_1, g_xxxzz_0_xxzzzz_1, g_xxxzz_0_xxzzzzz_1, g_xxxzz_0_xyyyyyz_1, g_xxxzz_0_xyyyyz_1, g_xxxzz_0_xyyyyzz_1, g_xxxzz_0_xyyyzz_1, g_xxxzz_0_xyyyzzz_1, g_xxxzz_0_xyyzzz_1, g_xxxzz_0_xyyzzzz_1, g_xxxzz_0_xyzzzz_1, g_xxxzz_0_xyzzzzz_1, g_xxxzz_0_xzzzzz_1, g_xxxzz_0_xzzzzzz_1, g_xxxzz_0_yyyyyyy_1, g_xxxzz_0_yyyyyyz_1, g_xxxzz_0_yyyyyz_1, g_xxxzz_0_yyyyyzz_1, g_xxxzz_0_yyyyzz_1, g_xxxzz_0_yyyyzzz_1, g_xxxzz_0_yyyzzz_1, g_xxxzz_0_yyyzzzz_1, g_xxxzz_0_yyzzzz_1, g_xxxzz_0_yyzzzzz_1, g_xxxzz_0_yzzzzz_1, g_xxxzz_0_yzzzzzz_1, g_xxxzz_0_zzzzzz_1, g_xxxzz_0_zzzzzzz_1, g_xxzz_0_xxxxxxz_0, g_xxzz_0_xxxxxxz_1, g_xxzz_0_xxxxxyz_0, g_xxzz_0_xxxxxyz_1, g_xxzz_0_xxxxxzz_0, g_xxzz_0_xxxxxzz_1, g_xxzz_0_xxxxyyz_0, g_xxzz_0_xxxxyyz_1, g_xxzz_0_xxxxyzz_0, g_xxzz_0_xxxxyzz_1, g_xxzz_0_xxxxzzz_0, g_xxzz_0_xxxxzzz_1, g_xxzz_0_xxxyyyz_0, g_xxzz_0_xxxyyyz_1, g_xxzz_0_xxxyyzz_0, g_xxzz_0_xxxyyzz_1, g_xxzz_0_xxxyzzz_0, g_xxzz_0_xxxyzzz_1, g_xxzz_0_xxxzzzz_0, g_xxzz_0_xxxzzzz_1, g_xxzz_0_xxyyyyz_0, g_xxzz_0_xxyyyyz_1, g_xxzz_0_xxyyyzz_0, g_xxzz_0_xxyyyzz_1, g_xxzz_0_xxyyzzz_0, g_xxzz_0_xxyyzzz_1, g_xxzz_0_xxyzzzz_0, g_xxzz_0_xxyzzzz_1, g_xxzz_0_xxzzzzz_0, g_xxzz_0_xxzzzzz_1, g_xxzz_0_xyyyyyz_0, g_xxzz_0_xyyyyyz_1, g_xxzz_0_xyyyyzz_0, g_xxzz_0_xyyyyzz_1, g_xxzz_0_xyyyzzz_0, g_xxzz_0_xyyyzzz_1, g_xxzz_0_xyyzzzz_0, g_xxzz_0_xyyzzzz_1, g_xxzz_0_xyzzzzz_0, g_xxzz_0_xyzzzzz_1, g_xxzz_0_xzzzzzz_0, g_xxzz_0_xzzzzzz_1, g_xxzz_0_yyyyyyy_0, g_xxzz_0_yyyyyyy_1, g_xxzz_0_yyyyyyz_0, g_xxzz_0_yyyyyyz_1, g_xxzz_0_yyyyyzz_0, g_xxzz_0_yyyyyzz_1, g_xxzz_0_yyyyzzz_0, g_xxzz_0_yyyyzzz_1, g_xxzz_0_yyyzzzz_0, g_xxzz_0_yyyzzzz_1, g_xxzz_0_yyzzzzz_0, g_xxzz_0_yyzzzzz_1, g_xxzz_0_yzzzzzz_0, g_xxzz_0_yzzzzzz_1, g_xxzz_0_zzzzzzz_0, g_xxzz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzz_0_xxxxxxx_0[i] = g_xxxx_0_xxxxxxx_0[i] * fbe_0 - g_xxxx_0_xxxxxxx_1[i] * fz_be_0 + g_xxxxz_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxxxy_0[i] = g_xxxx_0_xxxxxxy_0[i] * fbe_0 - g_xxxx_0_xxxxxxy_1[i] * fz_be_0 + g_xxxxz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxxxz_0[i] = 3.0 * g_xxzz_0_xxxxxxz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xxxzz_0_xxxxxz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxxyy_0[i] = g_xxxx_0_xxxxxyy_0[i] * fbe_0 - g_xxxx_0_xxxxxyy_1[i] * fz_be_0 + g_xxxxz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxxyz_0[i] = 3.0 * g_xxzz_0_xxxxxyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxxzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxxzz_0[i] = 3.0 * g_xxzz_0_xxxxxzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xxxzz_0_xxxxzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxyyy_0[i] = g_xxxx_0_xxxxyyy_0[i] * fbe_0 - g_xxxx_0_xxxxyyy_1[i] * fz_be_0 + g_xxxxz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxxyyz_0[i] = 3.0 * g_xxzz_0_xxxxyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxxzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxyzz_0[i] = 3.0 * g_xxzz_0_xxxxyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxxzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxxzzz_0[i] = 3.0 * g_xxzz_0_xxxxzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xxxzz_0_xxxzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxyyyy_0[i] = g_xxxx_0_xxxyyyy_0[i] * fbe_0 - g_xxxx_0_xxxyyyy_1[i] * fz_be_0 + g_xxxxz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxxyyyz_0[i] = 3.0 * g_xxzz_0_xxxyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxyyzz_0[i] = 3.0 * g_xxzz_0_xxxyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxyzzz_0[i] = 3.0 * g_xxzz_0_xxxyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxxzzzz_0[i] = 3.0 * g_xxzz_0_xxxzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xxxzz_0_xxzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyyyyy_0[i] = g_xxxx_0_xxyyyyy_0[i] * fbe_0 - g_xxxx_0_xxyyyyy_1[i] * fz_be_0 + g_xxxxz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xxyyyyz_0[i] = 3.0 * g_xxzz_0_xxyyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyyyzz_0[i] = 3.0 * g_xxzz_0_xxyyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyyzzz_0[i] = 3.0 * g_xxzz_0_xxyyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxyzzzz_0[i] = 3.0 * g_xxzz_0_xxyzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xxzzzzz_0[i] = 3.0 * g_xxzz_0_xxzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzz_0_xzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyyyyy_0[i] = g_xxxx_0_xyyyyyy_0[i] * fbe_0 - g_xxxx_0_xyyyyyy_1[i] * fz_be_0 + g_xxxxz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxxzz_0_xyyyyyz_0[i] = 3.0 * g_xxzz_0_xyyyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyyyz_1[i] * fz_be_0 + g_xxxzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyyyzz_0[i] = 3.0 * g_xxzz_0_xyyyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyyzz_1[i] * fz_be_0 + g_xxxzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyyzzz_0[i] = 3.0 * g_xxzz_0_xyyyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyzzz_1[i] * fz_be_0 + g_xxxzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyyzzzz_0[i] = 3.0 * g_xxzz_0_xyyzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyzzzz_1[i] * fz_be_0 + g_xxxzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xyzzzzz_0[i] = 3.0 * g_xxzz_0_xyzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyzzzzz_1[i] * fz_be_0 + g_xxxzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_xzzzzzz_0[i] = 3.0 * g_xxzz_0_xzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_xzzzzzz_1[i] * fz_be_0 + g_xxxzz_0_zzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyyyy_0[i] = 3.0 * g_xxzz_0_yyyyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyyyy_1[i] * fz_be_0 + g_xxxzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyyyz_0[i] = 3.0 * g_xxzz_0_yyyyyyz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyyyz_1[i] * fz_be_0 + g_xxxzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyyzz_0[i] = 3.0 * g_xxzz_0_yyyyyzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyyzz_1[i] * fz_be_0 + g_xxxzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyyzzz_0[i] = 3.0 * g_xxzz_0_yyyyzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyyzzz_1[i] * fz_be_0 + g_xxxzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyyzzzz_0[i] = 3.0 * g_xxzz_0_yyyzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyyzzzz_1[i] * fz_be_0 + g_xxxzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yyzzzzz_0[i] = 3.0 * g_xxzz_0_yyzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yyzzzzz_1[i] * fz_be_0 + g_xxxzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_yzzzzzz_0[i] = 3.0 * g_xxzz_0_yzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_yzzzzzz_1[i] * fz_be_0 + g_xxxzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxxzz_0_zzzzzzz_0[i] = 3.0 * g_xxzz_0_zzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzz_0_zzzzzzz_1[i] * fz_be_0 + g_xxxzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 216-252 components of targeted buffer : ISK

    auto g_xxxyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 216);

    auto g_xxxyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 217);

    auto g_xxxyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 218);

    auto g_xxxyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 219);

    auto g_xxxyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 220);

    auto g_xxxyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 221);

    auto g_xxxyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 222);

    auto g_xxxyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 223);

    auto g_xxxyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 224);

    auto g_xxxyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 225);

    auto g_xxxyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 226);

    auto g_xxxyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 227);

    auto g_xxxyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 228);

    auto g_xxxyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 229);

    auto g_xxxyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 230);

    auto g_xxxyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 231);

    auto g_xxxyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 232);

    auto g_xxxyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 233);

    auto g_xxxyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 234);

    auto g_xxxyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 235);

    auto g_xxxyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 236);

    auto g_xxxyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 237);

    auto g_xxxyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 238);

    auto g_xxxyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 239);

    auto g_xxxyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 240);

    auto g_xxxyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 241);

    auto g_xxxyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 242);

    auto g_xxxyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 243);

    auto g_xxxyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 244);

    auto g_xxxyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 245);

    auto g_xxxyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 246);

    auto g_xxxyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 247);

    auto g_xxxyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 248);

    auto g_xxxyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 249);

    auto g_xxxyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 250);

    auto g_xxxyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 251);

    #pragma omp simd aligned(g_xxxy_0_xxxxxxx_0, g_xxxy_0_xxxxxxx_1, g_xxxy_0_xxxxxxz_0, g_xxxy_0_xxxxxxz_1, g_xxxy_0_xxxxxzz_0, g_xxxy_0_xxxxxzz_1, g_xxxy_0_xxxxzzz_0, g_xxxy_0_xxxxzzz_1, g_xxxy_0_xxxzzzz_0, g_xxxy_0_xxxzzzz_1, g_xxxy_0_xxzzzzz_0, g_xxxy_0_xxzzzzz_1, g_xxxy_0_xzzzzzz_0, g_xxxy_0_xzzzzzz_1, g_xxxyy_0_xxxxxxx_1, g_xxxyy_0_xxxxxxz_1, g_xxxyy_0_xxxxxzz_1, g_xxxyy_0_xxxxzzz_1, g_xxxyy_0_xxxzzzz_1, g_xxxyy_0_xxzzzzz_1, g_xxxyy_0_xzzzzzz_1, g_xxxyyy_0_xxxxxxx_0, g_xxxyyy_0_xxxxxxy_0, g_xxxyyy_0_xxxxxxz_0, g_xxxyyy_0_xxxxxyy_0, g_xxxyyy_0_xxxxxyz_0, g_xxxyyy_0_xxxxxzz_0, g_xxxyyy_0_xxxxyyy_0, g_xxxyyy_0_xxxxyyz_0, g_xxxyyy_0_xxxxyzz_0, g_xxxyyy_0_xxxxzzz_0, g_xxxyyy_0_xxxyyyy_0, g_xxxyyy_0_xxxyyyz_0, g_xxxyyy_0_xxxyyzz_0, g_xxxyyy_0_xxxyzzz_0, g_xxxyyy_0_xxxzzzz_0, g_xxxyyy_0_xxyyyyy_0, g_xxxyyy_0_xxyyyyz_0, g_xxxyyy_0_xxyyyzz_0, g_xxxyyy_0_xxyyzzz_0, g_xxxyyy_0_xxyzzzz_0, g_xxxyyy_0_xxzzzzz_0, g_xxxyyy_0_xyyyyyy_0, g_xxxyyy_0_xyyyyyz_0, g_xxxyyy_0_xyyyyzz_0, g_xxxyyy_0_xyyyzzz_0, g_xxxyyy_0_xyyzzzz_0, g_xxxyyy_0_xyzzzzz_0, g_xxxyyy_0_xzzzzzz_0, g_xxxyyy_0_yyyyyyy_0, g_xxxyyy_0_yyyyyyz_0, g_xxxyyy_0_yyyyyzz_0, g_xxxyyy_0_yyyyzzz_0, g_xxxyyy_0_yyyzzzz_0, g_xxxyyy_0_yyzzzzz_0, g_xxxyyy_0_yzzzzzz_0, g_xxxyyy_0_zzzzzzz_0, g_xxyyy_0_xxxxxxy_1, g_xxyyy_0_xxxxxy_1, g_xxyyy_0_xxxxxyy_1, g_xxyyy_0_xxxxxyz_1, g_xxyyy_0_xxxxyy_1, g_xxyyy_0_xxxxyyy_1, g_xxyyy_0_xxxxyyz_1, g_xxyyy_0_xxxxyz_1, g_xxyyy_0_xxxxyzz_1, g_xxyyy_0_xxxyyy_1, g_xxyyy_0_xxxyyyy_1, g_xxyyy_0_xxxyyyz_1, g_xxyyy_0_xxxyyz_1, g_xxyyy_0_xxxyyzz_1, g_xxyyy_0_xxxyzz_1, g_xxyyy_0_xxxyzzz_1, g_xxyyy_0_xxyyyy_1, g_xxyyy_0_xxyyyyy_1, g_xxyyy_0_xxyyyyz_1, g_xxyyy_0_xxyyyz_1, g_xxyyy_0_xxyyyzz_1, g_xxyyy_0_xxyyzz_1, g_xxyyy_0_xxyyzzz_1, g_xxyyy_0_xxyzzz_1, g_xxyyy_0_xxyzzzz_1, g_xxyyy_0_xyyyyy_1, g_xxyyy_0_xyyyyyy_1, g_xxyyy_0_xyyyyyz_1, g_xxyyy_0_xyyyyz_1, g_xxyyy_0_xyyyyzz_1, g_xxyyy_0_xyyyzz_1, g_xxyyy_0_xyyyzzz_1, g_xxyyy_0_xyyzzz_1, g_xxyyy_0_xyyzzzz_1, g_xxyyy_0_xyzzzz_1, g_xxyyy_0_xyzzzzz_1, g_xxyyy_0_yyyyyy_1, g_xxyyy_0_yyyyyyy_1, g_xxyyy_0_yyyyyyz_1, g_xxyyy_0_yyyyyz_1, g_xxyyy_0_yyyyyzz_1, g_xxyyy_0_yyyyzz_1, g_xxyyy_0_yyyyzzz_1, g_xxyyy_0_yyyzzz_1, g_xxyyy_0_yyyzzzz_1, g_xxyyy_0_yyzzzz_1, g_xxyyy_0_yyzzzzz_1, g_xxyyy_0_yzzzzz_1, g_xxyyy_0_yzzzzzz_1, g_xxyyy_0_zzzzzzz_1, g_xyyy_0_xxxxxxy_0, g_xyyy_0_xxxxxxy_1, g_xyyy_0_xxxxxyy_0, g_xyyy_0_xxxxxyy_1, g_xyyy_0_xxxxxyz_0, g_xyyy_0_xxxxxyz_1, g_xyyy_0_xxxxyyy_0, g_xyyy_0_xxxxyyy_1, g_xyyy_0_xxxxyyz_0, g_xyyy_0_xxxxyyz_1, g_xyyy_0_xxxxyzz_0, g_xyyy_0_xxxxyzz_1, g_xyyy_0_xxxyyyy_0, g_xyyy_0_xxxyyyy_1, g_xyyy_0_xxxyyyz_0, g_xyyy_0_xxxyyyz_1, g_xyyy_0_xxxyyzz_0, g_xyyy_0_xxxyyzz_1, g_xyyy_0_xxxyzzz_0, g_xyyy_0_xxxyzzz_1, g_xyyy_0_xxyyyyy_0, g_xyyy_0_xxyyyyy_1, g_xyyy_0_xxyyyyz_0, g_xyyy_0_xxyyyyz_1, g_xyyy_0_xxyyyzz_0, g_xyyy_0_xxyyyzz_1, g_xyyy_0_xxyyzzz_0, g_xyyy_0_xxyyzzz_1, g_xyyy_0_xxyzzzz_0, g_xyyy_0_xxyzzzz_1, g_xyyy_0_xyyyyyy_0, g_xyyy_0_xyyyyyy_1, g_xyyy_0_xyyyyyz_0, g_xyyy_0_xyyyyyz_1, g_xyyy_0_xyyyyzz_0, g_xyyy_0_xyyyyzz_1, g_xyyy_0_xyyyzzz_0, g_xyyy_0_xyyyzzz_1, g_xyyy_0_xyyzzzz_0, g_xyyy_0_xyyzzzz_1, g_xyyy_0_xyzzzzz_0, g_xyyy_0_xyzzzzz_1, g_xyyy_0_yyyyyyy_0, g_xyyy_0_yyyyyyy_1, g_xyyy_0_yyyyyyz_0, g_xyyy_0_yyyyyyz_1, g_xyyy_0_yyyyyzz_0, g_xyyy_0_yyyyyzz_1, g_xyyy_0_yyyyzzz_0, g_xyyy_0_yyyyzzz_1, g_xyyy_0_yyyzzzz_0, g_xyyy_0_yyyzzzz_1, g_xyyy_0_yyzzzzz_0, g_xyyy_0_yyzzzzz_1, g_xyyy_0_yzzzzzz_0, g_xyyy_0_yzzzzzz_1, g_xyyy_0_zzzzzzz_0, g_xyyy_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyy_0_xxxxxxx_0[i] = 2.0 * g_xxxy_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxxxx_1[i] * fz_be_0 + g_xxxyy_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxyyy_0_xxxxxxy_0[i] = 2.0 * g_xyyy_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xxyyy_0_xxxxxy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxxxz_0[i] = 2.0 * g_xxxy_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxxxz_1[i] * fz_be_0 + g_xxxyy_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxyyy_0_xxxxxyy_0[i] = 2.0 * g_xyyy_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xxyyy_0_xxxxyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxxyz_0[i] = 2.0 * g_xyyy_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxyyy_0_xxxxyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxxzz_0[i] = 2.0 * g_xxxy_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxxzz_1[i] * fz_be_0 + g_xxxyy_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxyyy_0_xxxxyyy_0[i] = 2.0 * g_xyyy_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xxyyy_0_xxxyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxyyz_0[i] = 2.0 * g_xyyy_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxyyy_0_xxxyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxyzz_0[i] = 2.0 * g_xyyy_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxyyy_0_xxxyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxxzzz_0[i] = 2.0 * g_xxxy_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxxzzz_1[i] * fz_be_0 + g_xxxyy_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxyyy_0_xxxyyyy_0[i] = 2.0 * g_xyyy_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxxyyyz_0[i] = 2.0 * g_xyyy_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxyyzz_0[i] = 2.0 * g_xyyy_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxyzzz_0[i] = 2.0 * g_xyyy_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxyyy_0_xxyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxxzzzz_0[i] = 2.0 * g_xxxy_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxxzzzz_1[i] * fz_be_0 + g_xxxyy_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxyyy_0_xxyyyyy_0[i] = 2.0 * g_xyyy_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xxyyyyz_0[i] = 2.0 * g_xyyy_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xxyyyzz_0[i] = 2.0 * g_xyyy_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxyyzzz_0[i] = 2.0 * g_xyyy_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxyzzzz_0[i] = 2.0 * g_xyyy_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxyyy_0_xyzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xxzzzzz_0[i] = 2.0 * g_xxxy_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xxzzzzz_1[i] * fz_be_0 + g_xxxyy_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxyyy_0_xyyyyyy_0[i] = 2.0 * g_xyyy_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyyyy_1[i] * fz_be_0 + g_xxyyy_0_yyyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_xyyyyyz_0[i] = 2.0 * g_xyyy_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyyyz_1[i] * fz_be_0 + g_xxyyy_0_yyyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_xyyyyzz_0[i] = 2.0 * g_xyyy_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyyzz_1[i] * fz_be_0 + g_xxyyy_0_yyyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_xyyyzzz_0[i] = 2.0 * g_xyyy_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyyzzz_1[i] * fz_be_0 + g_xxyyy_0_yyyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xyyzzzz_0[i] = 2.0 * g_xyyy_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyyzzzz_1[i] * fz_be_0 + g_xxyyy_0_yyzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xyzzzzz_0[i] = 2.0 * g_xyyy_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_xyzzzzz_1[i] * fz_be_0 + g_xxyyy_0_yzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_xzzzzzz_0[i] = 2.0 * g_xxxy_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_xxxy_0_xzzzzzz_1[i] * fz_be_0 + g_xxxyy_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxyyy_0_yyyyyyy_0[i] = 2.0 * g_xyyy_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyyyy_1[i] * fz_be_0 + g_xxyyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxyyy_0_yyyyyyz_0[i] = 2.0 * g_xyyy_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyyyz_1[i] * fz_be_0 + g_xxyyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxyyy_0_yyyyyzz_0[i] = 2.0 * g_xyyy_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyyzz_1[i] * fz_be_0 + g_xxyyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxyyy_0_yyyyzzz_0[i] = 2.0 * g_xyyy_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyyzzz_1[i] * fz_be_0 + g_xxyyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxyyy_0_yyyzzzz_0[i] = 2.0 * g_xyyy_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyyzzzz_1[i] * fz_be_0 + g_xxyyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_yyzzzzz_0[i] = 2.0 * g_xyyy_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yyzzzzz_1[i] * fz_be_0 + g_xxyyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_yzzzzzz_0[i] = 2.0 * g_xyyy_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_yzzzzzz_1[i] * fz_be_0 + g_xxyyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxyyy_0_zzzzzzz_0[i] = 2.0 * g_xyyy_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyy_0_zzzzzzz_1[i] * fz_be_0 + g_xxyyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 252-288 components of targeted buffer : ISK

    auto g_xxxyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 252);

    auto g_xxxyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 253);

    auto g_xxxyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 254);

    auto g_xxxyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 255);

    auto g_xxxyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 256);

    auto g_xxxyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 257);

    auto g_xxxyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 258);

    auto g_xxxyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 259);

    auto g_xxxyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 260);

    auto g_xxxyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 261);

    auto g_xxxyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 262);

    auto g_xxxyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 263);

    auto g_xxxyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 264);

    auto g_xxxyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 265);

    auto g_xxxyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 266);

    auto g_xxxyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 267);

    auto g_xxxyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 268);

    auto g_xxxyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 269);

    auto g_xxxyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 270);

    auto g_xxxyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 271);

    auto g_xxxyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 272);

    auto g_xxxyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 273);

    auto g_xxxyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 274);

    auto g_xxxyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 275);

    auto g_xxxyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 276);

    auto g_xxxyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 277);

    auto g_xxxyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 278);

    auto g_xxxyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 279);

    auto g_xxxyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 280);

    auto g_xxxyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 281);

    auto g_xxxyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 282);

    auto g_xxxyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 283);

    auto g_xxxyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 284);

    auto g_xxxyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 285);

    auto g_xxxyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 286);

    auto g_xxxyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 287);

    #pragma omp simd aligned(g_xxxyy_0_xxxxxx_1, g_xxxyy_0_xxxxxxx_1, g_xxxyy_0_xxxxxxy_1, g_xxxyy_0_xxxxxxz_1, g_xxxyy_0_xxxxxy_1, g_xxxyy_0_xxxxxyy_1, g_xxxyy_0_xxxxxyz_1, g_xxxyy_0_xxxxxz_1, g_xxxyy_0_xxxxxzz_1, g_xxxyy_0_xxxxyy_1, g_xxxyy_0_xxxxyyy_1, g_xxxyy_0_xxxxyyz_1, g_xxxyy_0_xxxxyz_1, g_xxxyy_0_xxxxyzz_1, g_xxxyy_0_xxxxzz_1, g_xxxyy_0_xxxxzzz_1, g_xxxyy_0_xxxyyy_1, g_xxxyy_0_xxxyyyy_1, g_xxxyy_0_xxxyyyz_1, g_xxxyy_0_xxxyyz_1, g_xxxyy_0_xxxyyzz_1, g_xxxyy_0_xxxyzz_1, g_xxxyy_0_xxxyzzz_1, g_xxxyy_0_xxxzzz_1, g_xxxyy_0_xxxzzzz_1, g_xxxyy_0_xxyyyy_1, g_xxxyy_0_xxyyyyy_1, g_xxxyy_0_xxyyyyz_1, g_xxxyy_0_xxyyyz_1, g_xxxyy_0_xxyyyzz_1, g_xxxyy_0_xxyyzz_1, g_xxxyy_0_xxyyzzz_1, g_xxxyy_0_xxyzzz_1, g_xxxyy_0_xxyzzzz_1, g_xxxyy_0_xxzzzz_1, g_xxxyy_0_xxzzzzz_1, g_xxxyy_0_xyyyyy_1, g_xxxyy_0_xyyyyyy_1, g_xxxyy_0_xyyyyyz_1, g_xxxyy_0_xyyyyz_1, g_xxxyy_0_xyyyyzz_1, g_xxxyy_0_xyyyzz_1, g_xxxyy_0_xyyyzzz_1, g_xxxyy_0_xyyzzz_1, g_xxxyy_0_xyyzzzz_1, g_xxxyy_0_xyzzzz_1, g_xxxyy_0_xyzzzzz_1, g_xxxyy_0_xzzzzz_1, g_xxxyy_0_xzzzzzz_1, g_xxxyy_0_yyyyyy_1, g_xxxyy_0_yyyyyyy_1, g_xxxyy_0_yyyyyyz_1, g_xxxyy_0_yyyyyz_1, g_xxxyy_0_yyyyyzz_1, g_xxxyy_0_yyyyzz_1, g_xxxyy_0_yyyyzzz_1, g_xxxyy_0_yyyzzz_1, g_xxxyy_0_yyyzzzz_1, g_xxxyy_0_yyzzzz_1, g_xxxyy_0_yyzzzzz_1, g_xxxyy_0_yzzzzz_1, g_xxxyy_0_yzzzzzz_1, g_xxxyy_0_zzzzzz_1, g_xxxyy_0_zzzzzzz_1, g_xxxyyz_0_xxxxxxx_0, g_xxxyyz_0_xxxxxxy_0, g_xxxyyz_0_xxxxxxz_0, g_xxxyyz_0_xxxxxyy_0, g_xxxyyz_0_xxxxxyz_0, g_xxxyyz_0_xxxxxzz_0, g_xxxyyz_0_xxxxyyy_0, g_xxxyyz_0_xxxxyyz_0, g_xxxyyz_0_xxxxyzz_0, g_xxxyyz_0_xxxxzzz_0, g_xxxyyz_0_xxxyyyy_0, g_xxxyyz_0_xxxyyyz_0, g_xxxyyz_0_xxxyyzz_0, g_xxxyyz_0_xxxyzzz_0, g_xxxyyz_0_xxxzzzz_0, g_xxxyyz_0_xxyyyyy_0, g_xxxyyz_0_xxyyyyz_0, g_xxxyyz_0_xxyyyzz_0, g_xxxyyz_0_xxyyzzz_0, g_xxxyyz_0_xxyzzzz_0, g_xxxyyz_0_xxzzzzz_0, g_xxxyyz_0_xyyyyyy_0, g_xxxyyz_0_xyyyyyz_0, g_xxxyyz_0_xyyyyzz_0, g_xxxyyz_0_xyyyzzz_0, g_xxxyyz_0_xyyzzzz_0, g_xxxyyz_0_xyzzzzz_0, g_xxxyyz_0_xzzzzzz_0, g_xxxyyz_0_yyyyyyy_0, g_xxxyyz_0_yyyyyyz_0, g_xxxyyz_0_yyyyyzz_0, g_xxxyyz_0_yyyyzzz_0, g_xxxyyz_0_yyyzzzz_0, g_xxxyyz_0_yyzzzzz_0, g_xxxyyz_0_yzzzzzz_0, g_xxxyyz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyz_0_xxxxxxx_0[i] = g_xxxyy_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxxy_0[i] = g_xxxyy_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxxz_0[i] = g_xxxyy_0_xxxxxx_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxxz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxyy_0[i] = g_xxxyy_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxyz_0[i] = g_xxxyy_0_xxxxxy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxxzz_0[i] = 2.0 * g_xxxyy_0_xxxxxz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxxzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxyyy_0[i] = g_xxxyy_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxyyz_0[i] = g_xxxyy_0_xxxxyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxyzz_0[i] = 2.0 * g_xxxyy_0_xxxxyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxxzzz_0[i] = 3.0 * g_xxxyy_0_xxxxzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxxzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyyyy_0[i] = g_xxxyy_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyyyz_0[i] = g_xxxyy_0_xxxyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyyzz_0[i] = 2.0 * g_xxxyy_0_xxxyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxyzzz_0[i] = 3.0 * g_xxxyy_0_xxxyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxxzzzz_0[i] = 4.0 * g_xxxyy_0_xxxzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxxzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyyyy_0[i] = g_xxxyy_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyyyz_0[i] = g_xxxyy_0_xxyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyyzz_0[i] = 2.0 * g_xxxyy_0_xxyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyyzzz_0[i] = 3.0 * g_xxxyy_0_xxyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxyzzzz_0[i] = 4.0 * g_xxxyy_0_xxyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxyzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xxzzzzz_0[i] = 5.0 * g_xxxyy_0_xxzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xxzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyyyy_0[i] = g_xxxyy_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyyyz_0[i] = g_xxxyy_0_xyyyyy_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyyzz_0[i] = 2.0 * g_xxxyy_0_xyyyyz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyyzzz_0[i] = 3.0 * g_xxxyy_0_xyyyzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyyzzzz_0[i] = 4.0 * g_xxxyy_0_xyyzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyyzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xyzzzzz_0[i] = 5.0 * g_xxxyy_0_xyzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xyzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_xzzzzzz_0[i] = 6.0 * g_xxxyy_0_xzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_xzzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyyyy_0[i] = g_xxxyy_0_yyyyyyy_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyyyz_0[i] = g_xxxyy_0_yyyyyy_1[i] * fi_acd_0 + g_xxxyy_0_yyyyyyz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyyzz_0[i] = 2.0 * g_xxxyy_0_yyyyyz_1[i] * fi_acd_0 + g_xxxyy_0_yyyyyzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyyzzz_0[i] = 3.0 * g_xxxyy_0_yyyyzz_1[i] * fi_acd_0 + g_xxxyy_0_yyyyzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyyzzzz_0[i] = 4.0 * g_xxxyy_0_yyyzzz_1[i] * fi_acd_0 + g_xxxyy_0_yyyzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yyzzzzz_0[i] = 5.0 * g_xxxyy_0_yyzzzz_1[i] * fi_acd_0 + g_xxxyy_0_yyzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_yzzzzzz_0[i] = 6.0 * g_xxxyy_0_yzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_yzzzzzz_1[i] * wa_z[i];

        g_xxxyyz_0_zzzzzzz_0[i] = 7.0 * g_xxxyy_0_zzzzzz_1[i] * fi_acd_0 + g_xxxyy_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 288-324 components of targeted buffer : ISK

    auto g_xxxyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 288);

    auto g_xxxyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 289);

    auto g_xxxyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 290);

    auto g_xxxyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 291);

    auto g_xxxyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 292);

    auto g_xxxyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 293);

    auto g_xxxyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 294);

    auto g_xxxyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 295);

    auto g_xxxyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 296);

    auto g_xxxyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 297);

    auto g_xxxyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 298);

    auto g_xxxyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 299);

    auto g_xxxyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 300);

    auto g_xxxyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 301);

    auto g_xxxyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 302);

    auto g_xxxyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 303);

    auto g_xxxyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 304);

    auto g_xxxyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 305);

    auto g_xxxyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 306);

    auto g_xxxyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 307);

    auto g_xxxyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 308);

    auto g_xxxyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 309);

    auto g_xxxyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 310);

    auto g_xxxyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 311);

    auto g_xxxyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 312);

    auto g_xxxyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 313);

    auto g_xxxyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 314);

    auto g_xxxyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 315);

    auto g_xxxyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 316);

    auto g_xxxyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 317);

    auto g_xxxyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 318);

    auto g_xxxyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 319);

    auto g_xxxyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 320);

    auto g_xxxyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 321);

    auto g_xxxyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 322);

    auto g_xxxyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 323);

    #pragma omp simd aligned(g_xxxyzz_0_xxxxxxx_0, g_xxxyzz_0_xxxxxxy_0, g_xxxyzz_0_xxxxxxz_0, g_xxxyzz_0_xxxxxyy_0, g_xxxyzz_0_xxxxxyz_0, g_xxxyzz_0_xxxxxzz_0, g_xxxyzz_0_xxxxyyy_0, g_xxxyzz_0_xxxxyyz_0, g_xxxyzz_0_xxxxyzz_0, g_xxxyzz_0_xxxxzzz_0, g_xxxyzz_0_xxxyyyy_0, g_xxxyzz_0_xxxyyyz_0, g_xxxyzz_0_xxxyyzz_0, g_xxxyzz_0_xxxyzzz_0, g_xxxyzz_0_xxxzzzz_0, g_xxxyzz_0_xxyyyyy_0, g_xxxyzz_0_xxyyyyz_0, g_xxxyzz_0_xxyyyzz_0, g_xxxyzz_0_xxyyzzz_0, g_xxxyzz_0_xxyzzzz_0, g_xxxyzz_0_xxzzzzz_0, g_xxxyzz_0_xyyyyyy_0, g_xxxyzz_0_xyyyyyz_0, g_xxxyzz_0_xyyyyzz_0, g_xxxyzz_0_xyyyzzz_0, g_xxxyzz_0_xyyzzzz_0, g_xxxyzz_0_xyzzzzz_0, g_xxxyzz_0_xzzzzzz_0, g_xxxyzz_0_yyyyyyy_0, g_xxxyzz_0_yyyyyyz_0, g_xxxyzz_0_yyyyyzz_0, g_xxxyzz_0_yyyyzzz_0, g_xxxyzz_0_yyyzzzz_0, g_xxxyzz_0_yyzzzzz_0, g_xxxyzz_0_yzzzzzz_0, g_xxxyzz_0_zzzzzzz_0, g_xxxzz_0_xxxxxx_1, g_xxxzz_0_xxxxxxx_1, g_xxxzz_0_xxxxxxy_1, g_xxxzz_0_xxxxxxz_1, g_xxxzz_0_xxxxxy_1, g_xxxzz_0_xxxxxyy_1, g_xxxzz_0_xxxxxyz_1, g_xxxzz_0_xxxxxz_1, g_xxxzz_0_xxxxxzz_1, g_xxxzz_0_xxxxyy_1, g_xxxzz_0_xxxxyyy_1, g_xxxzz_0_xxxxyyz_1, g_xxxzz_0_xxxxyz_1, g_xxxzz_0_xxxxyzz_1, g_xxxzz_0_xxxxzz_1, g_xxxzz_0_xxxxzzz_1, g_xxxzz_0_xxxyyy_1, g_xxxzz_0_xxxyyyy_1, g_xxxzz_0_xxxyyyz_1, g_xxxzz_0_xxxyyz_1, g_xxxzz_0_xxxyyzz_1, g_xxxzz_0_xxxyzz_1, g_xxxzz_0_xxxyzzz_1, g_xxxzz_0_xxxzzz_1, g_xxxzz_0_xxxzzzz_1, g_xxxzz_0_xxyyyy_1, g_xxxzz_0_xxyyyyy_1, g_xxxzz_0_xxyyyyz_1, g_xxxzz_0_xxyyyz_1, g_xxxzz_0_xxyyyzz_1, g_xxxzz_0_xxyyzz_1, g_xxxzz_0_xxyyzzz_1, g_xxxzz_0_xxyzzz_1, g_xxxzz_0_xxyzzzz_1, g_xxxzz_0_xxzzzz_1, g_xxxzz_0_xxzzzzz_1, g_xxxzz_0_xyyyyy_1, g_xxxzz_0_xyyyyyy_1, g_xxxzz_0_xyyyyyz_1, g_xxxzz_0_xyyyyz_1, g_xxxzz_0_xyyyyzz_1, g_xxxzz_0_xyyyzz_1, g_xxxzz_0_xyyyzzz_1, g_xxxzz_0_xyyzzz_1, g_xxxzz_0_xyyzzzz_1, g_xxxzz_0_xyzzzz_1, g_xxxzz_0_xyzzzzz_1, g_xxxzz_0_xzzzzz_1, g_xxxzz_0_xzzzzzz_1, g_xxxzz_0_yyyyyy_1, g_xxxzz_0_yyyyyyy_1, g_xxxzz_0_yyyyyyz_1, g_xxxzz_0_yyyyyz_1, g_xxxzz_0_yyyyyzz_1, g_xxxzz_0_yyyyzz_1, g_xxxzz_0_yyyyzzz_1, g_xxxzz_0_yyyzzz_1, g_xxxzz_0_yyyzzzz_1, g_xxxzz_0_yyzzzz_1, g_xxxzz_0_yyzzzzz_1, g_xxxzz_0_yzzzzz_1, g_xxxzz_0_yzzzzzz_1, g_xxxzz_0_zzzzzz_1, g_xxxzz_0_zzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzz_0_xxxxxxx_0[i] = g_xxxzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxxy_0[i] = g_xxxzz_0_xxxxxx_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxxy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxxz_0[i] = g_xxxzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxyy_0[i] = 2.0 * g_xxxzz_0_xxxxxy_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxyz_0[i] = g_xxxzz_0_xxxxxz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxxyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxxzz_0[i] = g_xxxzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxyyy_0[i] = 3.0 * g_xxxzz_0_xxxxyy_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxyyz_0[i] = 2.0 * g_xxxzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxyzz_0[i] = g_xxxzz_0_xxxxzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxxyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxxzzz_0[i] = g_xxxzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyyyy_0[i] = 4.0 * g_xxxzz_0_xxxyyy_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyyyz_0[i] = 3.0 * g_xxxzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyyzz_0[i] = 2.0 * g_xxxzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxyzzz_0[i] = g_xxxzz_0_xxxzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxxyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxxzzzz_0[i] = g_xxxzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyyyy_0[i] = 5.0 * g_xxxzz_0_xxyyyy_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyyyz_0[i] = 4.0 * g_xxxzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyyzz_0[i] = 3.0 * g_xxxzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyyzzz_0[i] = 2.0 * g_xxxzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxyzzzz_0[i] = g_xxxzz_0_xxzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xxyzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xxzzzzz_0[i] = g_xxxzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyyyy_0[i] = 6.0 * g_xxxzz_0_xyyyyy_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyyyz_0[i] = 5.0 * g_xxxzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyyzz_0[i] = 4.0 * g_xxxzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyyzzz_0[i] = 3.0 * g_xxxzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyyzzzz_0[i] = 2.0 * g_xxxzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyyzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xyzzzzz_0[i] = g_xxxzz_0_xzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_xyzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_xzzzzzz_0[i] = g_xxxzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyyyy_0[i] = 7.0 * g_xxxzz_0_yyyyyy_1[i] * fi_acd_0 + g_xxxzz_0_yyyyyyy_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyyyz_0[i] = 6.0 * g_xxxzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxxzz_0_yyyyyyz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyyzz_0[i] = 5.0 * g_xxxzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxxzz_0_yyyyyzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyyzzz_0[i] = 4.0 * g_xxxzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxxzz_0_yyyyzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyyzzzz_0[i] = 3.0 * g_xxxzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxxzz_0_yyyzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yyzzzzz_0[i] = 2.0 * g_xxxzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_yyzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_yzzzzzz_0[i] = g_xxxzz_0_zzzzzz_1[i] * fi_acd_0 + g_xxxzz_0_yzzzzzz_1[i] * wa_y[i];

        g_xxxyzz_0_zzzzzzz_0[i] = g_xxxzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 324-360 components of targeted buffer : ISK

    auto g_xxxzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 324);

    auto g_xxxzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 325);

    auto g_xxxzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 326);

    auto g_xxxzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 327);

    auto g_xxxzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 328);

    auto g_xxxzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 329);

    auto g_xxxzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 330);

    auto g_xxxzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 331);

    auto g_xxxzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 332);

    auto g_xxxzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 333);

    auto g_xxxzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 334);

    auto g_xxxzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 335);

    auto g_xxxzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 336);

    auto g_xxxzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 337);

    auto g_xxxzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 338);

    auto g_xxxzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 339);

    auto g_xxxzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 340);

    auto g_xxxzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 341);

    auto g_xxxzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 342);

    auto g_xxxzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 343);

    auto g_xxxzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 344);

    auto g_xxxzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 345);

    auto g_xxxzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 346);

    auto g_xxxzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 347);

    auto g_xxxzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 348);

    auto g_xxxzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 349);

    auto g_xxxzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 350);

    auto g_xxxzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 351);

    auto g_xxxzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 352);

    auto g_xxxzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 353);

    auto g_xxxzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 354);

    auto g_xxxzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 355);

    auto g_xxxzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 356);

    auto g_xxxzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 357);

    auto g_xxxzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 358);

    auto g_xxxzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 359);

    #pragma omp simd aligned(g_xxxz_0_xxxxxxx_0, g_xxxz_0_xxxxxxx_1, g_xxxz_0_xxxxxxy_0, g_xxxz_0_xxxxxxy_1, g_xxxz_0_xxxxxyy_0, g_xxxz_0_xxxxxyy_1, g_xxxz_0_xxxxyyy_0, g_xxxz_0_xxxxyyy_1, g_xxxz_0_xxxyyyy_0, g_xxxz_0_xxxyyyy_1, g_xxxz_0_xxyyyyy_0, g_xxxz_0_xxyyyyy_1, g_xxxz_0_xyyyyyy_0, g_xxxz_0_xyyyyyy_1, g_xxxzz_0_xxxxxxx_1, g_xxxzz_0_xxxxxxy_1, g_xxxzz_0_xxxxxyy_1, g_xxxzz_0_xxxxyyy_1, g_xxxzz_0_xxxyyyy_1, g_xxxzz_0_xxyyyyy_1, g_xxxzz_0_xyyyyyy_1, g_xxxzzz_0_xxxxxxx_0, g_xxxzzz_0_xxxxxxy_0, g_xxxzzz_0_xxxxxxz_0, g_xxxzzz_0_xxxxxyy_0, g_xxxzzz_0_xxxxxyz_0, g_xxxzzz_0_xxxxxzz_0, g_xxxzzz_0_xxxxyyy_0, g_xxxzzz_0_xxxxyyz_0, g_xxxzzz_0_xxxxyzz_0, g_xxxzzz_0_xxxxzzz_0, g_xxxzzz_0_xxxyyyy_0, g_xxxzzz_0_xxxyyyz_0, g_xxxzzz_0_xxxyyzz_0, g_xxxzzz_0_xxxyzzz_0, g_xxxzzz_0_xxxzzzz_0, g_xxxzzz_0_xxyyyyy_0, g_xxxzzz_0_xxyyyyz_0, g_xxxzzz_0_xxyyyzz_0, g_xxxzzz_0_xxyyzzz_0, g_xxxzzz_0_xxyzzzz_0, g_xxxzzz_0_xxzzzzz_0, g_xxxzzz_0_xyyyyyy_0, g_xxxzzz_0_xyyyyyz_0, g_xxxzzz_0_xyyyyzz_0, g_xxxzzz_0_xyyyzzz_0, g_xxxzzz_0_xyyzzzz_0, g_xxxzzz_0_xyzzzzz_0, g_xxxzzz_0_xzzzzzz_0, g_xxxzzz_0_yyyyyyy_0, g_xxxzzz_0_yyyyyyz_0, g_xxxzzz_0_yyyyyzz_0, g_xxxzzz_0_yyyyzzz_0, g_xxxzzz_0_yyyzzzz_0, g_xxxzzz_0_yyzzzzz_0, g_xxxzzz_0_yzzzzzz_0, g_xxxzzz_0_zzzzzzz_0, g_xxzzz_0_xxxxxxz_1, g_xxzzz_0_xxxxxyz_1, g_xxzzz_0_xxxxxz_1, g_xxzzz_0_xxxxxzz_1, g_xxzzz_0_xxxxyyz_1, g_xxzzz_0_xxxxyz_1, g_xxzzz_0_xxxxyzz_1, g_xxzzz_0_xxxxzz_1, g_xxzzz_0_xxxxzzz_1, g_xxzzz_0_xxxyyyz_1, g_xxzzz_0_xxxyyz_1, g_xxzzz_0_xxxyyzz_1, g_xxzzz_0_xxxyzz_1, g_xxzzz_0_xxxyzzz_1, g_xxzzz_0_xxxzzz_1, g_xxzzz_0_xxxzzzz_1, g_xxzzz_0_xxyyyyz_1, g_xxzzz_0_xxyyyz_1, g_xxzzz_0_xxyyyzz_1, g_xxzzz_0_xxyyzz_1, g_xxzzz_0_xxyyzzz_1, g_xxzzz_0_xxyzzz_1, g_xxzzz_0_xxyzzzz_1, g_xxzzz_0_xxzzzz_1, g_xxzzz_0_xxzzzzz_1, g_xxzzz_0_xyyyyyz_1, g_xxzzz_0_xyyyyz_1, g_xxzzz_0_xyyyyzz_1, g_xxzzz_0_xyyyzz_1, g_xxzzz_0_xyyyzzz_1, g_xxzzz_0_xyyzzz_1, g_xxzzz_0_xyyzzzz_1, g_xxzzz_0_xyzzzz_1, g_xxzzz_0_xyzzzzz_1, g_xxzzz_0_xzzzzz_1, g_xxzzz_0_xzzzzzz_1, g_xxzzz_0_yyyyyyy_1, g_xxzzz_0_yyyyyyz_1, g_xxzzz_0_yyyyyz_1, g_xxzzz_0_yyyyyzz_1, g_xxzzz_0_yyyyzz_1, g_xxzzz_0_yyyyzzz_1, g_xxzzz_0_yyyzzz_1, g_xxzzz_0_yyyzzzz_1, g_xxzzz_0_yyzzzz_1, g_xxzzz_0_yyzzzzz_1, g_xxzzz_0_yzzzzz_1, g_xxzzz_0_yzzzzzz_1, g_xxzzz_0_zzzzzz_1, g_xxzzz_0_zzzzzzz_1, g_xzzz_0_xxxxxxz_0, g_xzzz_0_xxxxxxz_1, g_xzzz_0_xxxxxyz_0, g_xzzz_0_xxxxxyz_1, g_xzzz_0_xxxxxzz_0, g_xzzz_0_xxxxxzz_1, g_xzzz_0_xxxxyyz_0, g_xzzz_0_xxxxyyz_1, g_xzzz_0_xxxxyzz_0, g_xzzz_0_xxxxyzz_1, g_xzzz_0_xxxxzzz_0, g_xzzz_0_xxxxzzz_1, g_xzzz_0_xxxyyyz_0, g_xzzz_0_xxxyyyz_1, g_xzzz_0_xxxyyzz_0, g_xzzz_0_xxxyyzz_1, g_xzzz_0_xxxyzzz_0, g_xzzz_0_xxxyzzz_1, g_xzzz_0_xxxzzzz_0, g_xzzz_0_xxxzzzz_1, g_xzzz_0_xxyyyyz_0, g_xzzz_0_xxyyyyz_1, g_xzzz_0_xxyyyzz_0, g_xzzz_0_xxyyyzz_1, g_xzzz_0_xxyyzzz_0, g_xzzz_0_xxyyzzz_1, g_xzzz_0_xxyzzzz_0, g_xzzz_0_xxyzzzz_1, g_xzzz_0_xxzzzzz_0, g_xzzz_0_xxzzzzz_1, g_xzzz_0_xyyyyyz_0, g_xzzz_0_xyyyyyz_1, g_xzzz_0_xyyyyzz_0, g_xzzz_0_xyyyyzz_1, g_xzzz_0_xyyyzzz_0, g_xzzz_0_xyyyzzz_1, g_xzzz_0_xyyzzzz_0, g_xzzz_0_xyyzzzz_1, g_xzzz_0_xyzzzzz_0, g_xzzz_0_xyzzzzz_1, g_xzzz_0_xzzzzzz_0, g_xzzz_0_xzzzzzz_1, g_xzzz_0_yyyyyyy_0, g_xzzz_0_yyyyyyy_1, g_xzzz_0_yyyyyyz_0, g_xzzz_0_yyyyyyz_1, g_xzzz_0_yyyyyzz_0, g_xzzz_0_yyyyyzz_1, g_xzzz_0_yyyyzzz_0, g_xzzz_0_yyyyzzz_1, g_xzzz_0_yyyzzzz_0, g_xzzz_0_yyyzzzz_1, g_xzzz_0_yyzzzzz_0, g_xzzz_0_yyzzzzz_1, g_xzzz_0_yzzzzzz_0, g_xzzz_0_yzzzzzz_1, g_xzzz_0_zzzzzzz_0, g_xzzz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzz_0_xxxxxxx_0[i] = 2.0 * g_xxxz_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxxxx_1[i] * fz_be_0 + g_xxxzz_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxxxy_0[i] = 2.0 * g_xxxz_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxxxy_1[i] * fz_be_0 + g_xxxzz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxxxz_0[i] = 2.0 * g_xzzz_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xxzzz_0_xxxxxz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxxyy_0[i] = 2.0 * g_xxxz_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxxyy_1[i] * fz_be_0 + g_xxxzz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxxyz_0[i] = 2.0 * g_xzzz_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxzzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxxzz_0[i] = 2.0 * g_xzzz_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xxzzz_0_xxxxzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxyyy_0[i] = 2.0 * g_xxxz_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxxyyy_1[i] * fz_be_0 + g_xxxzz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxxyyz_0[i] = 2.0 * g_xzzz_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxzzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxyzz_0[i] = 2.0 * g_xzzz_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxzzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxxzzz_0[i] = 2.0 * g_xzzz_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xxzzz_0_xxxzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxyyyy_0[i] = 2.0 * g_xxxz_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxxyyyy_1[i] * fz_be_0 + g_xxxzz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxxyyyz_0[i] = 2.0 * g_xzzz_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxyyzz_0[i] = 2.0 * g_xzzz_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxyzzz_0[i] = 2.0 * g_xzzz_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxxzzzz_0[i] = 2.0 * g_xzzz_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xxzzz_0_xxzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyyyyy_0[i] = 2.0 * g_xxxz_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xxyyyyy_1[i] * fz_be_0 + g_xxxzz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xxyyyyz_0[i] = 2.0 * g_xzzz_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyyyzz_0[i] = 2.0 * g_xzzz_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyyzzz_0[i] = 2.0 * g_xzzz_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxyzzzz_0[i] = 2.0 * g_xzzz_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xxzzzzz_0[i] = 2.0 * g_xzzz_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzz_0_xzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyyyyy_0[i] = 2.0 * g_xxxz_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_xxxz_0_xyyyyyy_1[i] * fz_be_0 + g_xxxzz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxzzz_0_xyyyyyz_0[i] = 2.0 * g_xzzz_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyyyyz_1[i] * fz_be_0 + g_xxzzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyyyzz_0[i] = 2.0 * g_xzzz_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyyyzz_1[i] * fz_be_0 + g_xxzzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyyzzz_0[i] = 2.0 * g_xzzz_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyyzzz_1[i] * fz_be_0 + g_xxzzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyyzzzz_0[i] = 2.0 * g_xzzz_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyyzzzz_1[i] * fz_be_0 + g_xxzzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xyzzzzz_0[i] = 2.0 * g_xzzz_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xyzzzzz_1[i] * fz_be_0 + g_xxzzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_xzzzzzz_0[i] = 2.0 * g_xzzz_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_xzzzzzz_1[i] * fz_be_0 + g_xxzzz_0_zzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyyyy_0[i] = 2.0 * g_xzzz_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyyyy_1[i] * fz_be_0 + g_xxzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyyyz_0[i] = 2.0 * g_xzzz_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyyyz_1[i] * fz_be_0 + g_xxzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyyzz_0[i] = 2.0 * g_xzzz_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyyzz_1[i] * fz_be_0 + g_xxzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyyzzz_0[i] = 2.0 * g_xzzz_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyyzzz_1[i] * fz_be_0 + g_xxzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyyzzzz_0[i] = 2.0 * g_xzzz_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyyzzzz_1[i] * fz_be_0 + g_xxzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yyzzzzz_0[i] = 2.0 * g_xzzz_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yyzzzzz_1[i] * fz_be_0 + g_xxzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_yzzzzzz_0[i] = 2.0 * g_xzzz_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_yzzzzzz_1[i] * fz_be_0 + g_xxzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxzzz_0_zzzzzzz_0[i] = 2.0 * g_xzzz_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzz_0_zzzzzzz_1[i] * fz_be_0 + g_xxzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 360-396 components of targeted buffer : ISK

    auto g_xxyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 360);

    auto g_xxyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 361);

    auto g_xxyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 362);

    auto g_xxyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 363);

    auto g_xxyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 364);

    auto g_xxyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 365);

    auto g_xxyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 366);

    auto g_xxyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 367);

    auto g_xxyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 368);

    auto g_xxyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 369);

    auto g_xxyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 370);

    auto g_xxyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 371);

    auto g_xxyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 372);

    auto g_xxyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 373);

    auto g_xxyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 374);

    auto g_xxyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 375);

    auto g_xxyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 376);

    auto g_xxyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 377);

    auto g_xxyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 378);

    auto g_xxyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 379);

    auto g_xxyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 380);

    auto g_xxyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 381);

    auto g_xxyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 382);

    auto g_xxyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 383);

    auto g_xxyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 384);

    auto g_xxyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 385);

    auto g_xxyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 386);

    auto g_xxyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 387);

    auto g_xxyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 388);

    auto g_xxyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 389);

    auto g_xxyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 390);

    auto g_xxyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 391);

    auto g_xxyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 392);

    auto g_xxyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 393);

    auto g_xxyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 394);

    auto g_xxyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 395);

    #pragma omp simd aligned(g_xxyy_0_xxxxxxx_0, g_xxyy_0_xxxxxxx_1, g_xxyy_0_xxxxxxz_0, g_xxyy_0_xxxxxxz_1, g_xxyy_0_xxxxxzz_0, g_xxyy_0_xxxxxzz_1, g_xxyy_0_xxxxzzz_0, g_xxyy_0_xxxxzzz_1, g_xxyy_0_xxxzzzz_0, g_xxyy_0_xxxzzzz_1, g_xxyy_0_xxzzzzz_0, g_xxyy_0_xxzzzzz_1, g_xxyy_0_xzzzzzz_0, g_xxyy_0_xzzzzzz_1, g_xxyyy_0_xxxxxxx_1, g_xxyyy_0_xxxxxxz_1, g_xxyyy_0_xxxxxzz_1, g_xxyyy_0_xxxxzzz_1, g_xxyyy_0_xxxzzzz_1, g_xxyyy_0_xxzzzzz_1, g_xxyyy_0_xzzzzzz_1, g_xxyyyy_0_xxxxxxx_0, g_xxyyyy_0_xxxxxxy_0, g_xxyyyy_0_xxxxxxz_0, g_xxyyyy_0_xxxxxyy_0, g_xxyyyy_0_xxxxxyz_0, g_xxyyyy_0_xxxxxzz_0, g_xxyyyy_0_xxxxyyy_0, g_xxyyyy_0_xxxxyyz_0, g_xxyyyy_0_xxxxyzz_0, g_xxyyyy_0_xxxxzzz_0, g_xxyyyy_0_xxxyyyy_0, g_xxyyyy_0_xxxyyyz_0, g_xxyyyy_0_xxxyyzz_0, g_xxyyyy_0_xxxyzzz_0, g_xxyyyy_0_xxxzzzz_0, g_xxyyyy_0_xxyyyyy_0, g_xxyyyy_0_xxyyyyz_0, g_xxyyyy_0_xxyyyzz_0, g_xxyyyy_0_xxyyzzz_0, g_xxyyyy_0_xxyzzzz_0, g_xxyyyy_0_xxzzzzz_0, g_xxyyyy_0_xyyyyyy_0, g_xxyyyy_0_xyyyyyz_0, g_xxyyyy_0_xyyyyzz_0, g_xxyyyy_0_xyyyzzz_0, g_xxyyyy_0_xyyzzzz_0, g_xxyyyy_0_xyzzzzz_0, g_xxyyyy_0_xzzzzzz_0, g_xxyyyy_0_yyyyyyy_0, g_xxyyyy_0_yyyyyyz_0, g_xxyyyy_0_yyyyyzz_0, g_xxyyyy_0_yyyyzzz_0, g_xxyyyy_0_yyyzzzz_0, g_xxyyyy_0_yyzzzzz_0, g_xxyyyy_0_yzzzzzz_0, g_xxyyyy_0_zzzzzzz_0, g_xyyyy_0_xxxxxxy_1, g_xyyyy_0_xxxxxy_1, g_xyyyy_0_xxxxxyy_1, g_xyyyy_0_xxxxxyz_1, g_xyyyy_0_xxxxyy_1, g_xyyyy_0_xxxxyyy_1, g_xyyyy_0_xxxxyyz_1, g_xyyyy_0_xxxxyz_1, g_xyyyy_0_xxxxyzz_1, g_xyyyy_0_xxxyyy_1, g_xyyyy_0_xxxyyyy_1, g_xyyyy_0_xxxyyyz_1, g_xyyyy_0_xxxyyz_1, g_xyyyy_0_xxxyyzz_1, g_xyyyy_0_xxxyzz_1, g_xyyyy_0_xxxyzzz_1, g_xyyyy_0_xxyyyy_1, g_xyyyy_0_xxyyyyy_1, g_xyyyy_0_xxyyyyz_1, g_xyyyy_0_xxyyyz_1, g_xyyyy_0_xxyyyzz_1, g_xyyyy_0_xxyyzz_1, g_xyyyy_0_xxyyzzz_1, g_xyyyy_0_xxyzzz_1, g_xyyyy_0_xxyzzzz_1, g_xyyyy_0_xyyyyy_1, g_xyyyy_0_xyyyyyy_1, g_xyyyy_0_xyyyyyz_1, g_xyyyy_0_xyyyyz_1, g_xyyyy_0_xyyyyzz_1, g_xyyyy_0_xyyyzz_1, g_xyyyy_0_xyyyzzz_1, g_xyyyy_0_xyyzzz_1, g_xyyyy_0_xyyzzzz_1, g_xyyyy_0_xyzzzz_1, g_xyyyy_0_xyzzzzz_1, g_xyyyy_0_yyyyyy_1, g_xyyyy_0_yyyyyyy_1, g_xyyyy_0_yyyyyyz_1, g_xyyyy_0_yyyyyz_1, g_xyyyy_0_yyyyyzz_1, g_xyyyy_0_yyyyzz_1, g_xyyyy_0_yyyyzzz_1, g_xyyyy_0_yyyzzz_1, g_xyyyy_0_yyyzzzz_1, g_xyyyy_0_yyzzzz_1, g_xyyyy_0_yyzzzzz_1, g_xyyyy_0_yzzzzz_1, g_xyyyy_0_yzzzzzz_1, g_xyyyy_0_zzzzzzz_1, g_yyyy_0_xxxxxxy_0, g_yyyy_0_xxxxxxy_1, g_yyyy_0_xxxxxyy_0, g_yyyy_0_xxxxxyy_1, g_yyyy_0_xxxxxyz_0, g_yyyy_0_xxxxxyz_1, g_yyyy_0_xxxxyyy_0, g_yyyy_0_xxxxyyy_1, g_yyyy_0_xxxxyyz_0, g_yyyy_0_xxxxyyz_1, g_yyyy_0_xxxxyzz_0, g_yyyy_0_xxxxyzz_1, g_yyyy_0_xxxyyyy_0, g_yyyy_0_xxxyyyy_1, g_yyyy_0_xxxyyyz_0, g_yyyy_0_xxxyyyz_1, g_yyyy_0_xxxyyzz_0, g_yyyy_0_xxxyyzz_1, g_yyyy_0_xxxyzzz_0, g_yyyy_0_xxxyzzz_1, g_yyyy_0_xxyyyyy_0, g_yyyy_0_xxyyyyy_1, g_yyyy_0_xxyyyyz_0, g_yyyy_0_xxyyyyz_1, g_yyyy_0_xxyyyzz_0, g_yyyy_0_xxyyyzz_1, g_yyyy_0_xxyyzzz_0, g_yyyy_0_xxyyzzz_1, g_yyyy_0_xxyzzzz_0, g_yyyy_0_xxyzzzz_1, g_yyyy_0_xyyyyyy_0, g_yyyy_0_xyyyyyy_1, g_yyyy_0_xyyyyyz_0, g_yyyy_0_xyyyyyz_1, g_yyyy_0_xyyyyzz_0, g_yyyy_0_xyyyyzz_1, g_yyyy_0_xyyyzzz_0, g_yyyy_0_xyyyzzz_1, g_yyyy_0_xyyzzzz_0, g_yyyy_0_xyyzzzz_1, g_yyyy_0_xyzzzzz_0, g_yyyy_0_xyzzzzz_1, g_yyyy_0_yyyyyyy_0, g_yyyy_0_yyyyyyy_1, g_yyyy_0_yyyyyyz_0, g_yyyy_0_yyyyyyz_1, g_yyyy_0_yyyyyzz_0, g_yyyy_0_yyyyyzz_1, g_yyyy_0_yyyyzzz_0, g_yyyy_0_yyyyzzz_1, g_yyyy_0_yyyzzzz_0, g_yyyy_0_yyyzzzz_1, g_yyyy_0_yyzzzzz_0, g_yyyy_0_yyzzzzz_1, g_yyyy_0_yzzzzzz_0, g_yyyy_0_yzzzzzz_1, g_yyyy_0_zzzzzzz_0, g_yyyy_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyy_0_xxxxxxx_0[i] = 3.0 * g_xxyy_0_xxxxxxx_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxxx_1[i] * fz_be_0 + g_xxyyy_0_xxxxxxx_1[i] * wa_y[i];

        g_xxyyyy_0_xxxxxxy_0[i] = g_yyyy_0_xxxxxxy_0[i] * fbe_0 - g_yyyy_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xyyyy_0_xxxxxy_1[i] * fi_acd_0 + g_xyyyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxxxz_0[i] = 3.0 * g_xxyy_0_xxxxxxz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxxz_1[i] * fz_be_0 + g_xxyyy_0_xxxxxxz_1[i] * wa_y[i];

        g_xxyyyy_0_xxxxxyy_0[i] = g_yyyy_0_xxxxxyy_0[i] * fbe_0 - g_yyyy_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xyyyy_0_xxxxyy_1[i] * fi_acd_0 + g_xyyyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxxyz_0[i] = g_yyyy_0_xxxxxyz_0[i] * fbe_0 - g_yyyy_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xyyyy_0_xxxxyz_1[i] * fi_acd_0 + g_xyyyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxxzz_0[i] = 3.0 * g_xxyy_0_xxxxxzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxxzz_1[i] * fz_be_0 + g_xxyyy_0_xxxxxzz_1[i] * wa_y[i];

        g_xxyyyy_0_xxxxyyy_0[i] = g_yyyy_0_xxxxyyy_0[i] * fbe_0 - g_yyyy_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xyyyy_0_xxxyyy_1[i] * fi_acd_0 + g_xyyyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxyyz_0[i] = g_yyyy_0_xxxxyyz_0[i] * fbe_0 - g_yyyy_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xyyyy_0_xxxyyz_1[i] * fi_acd_0 + g_xyyyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxyzz_0[i] = g_yyyy_0_xxxxyzz_0[i] * fbe_0 - g_yyyy_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xyyyy_0_xxxyzz_1[i] * fi_acd_0 + g_xyyyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxxzzz_0[i] = 3.0 * g_xxyy_0_xxxxzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxxzzz_1[i] * fz_be_0 + g_xxyyy_0_xxxxzzz_1[i] * wa_y[i];

        g_xxyyyy_0_xxxyyyy_0[i] = g_yyyy_0_xxxyyyy_0[i] * fbe_0 - g_yyyy_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyyyy_1[i] * fi_acd_0 + g_xyyyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxxyyyz_0[i] = g_yyyy_0_xxxyyyz_0[i] * fbe_0 - g_yyyy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyyyz_1[i] * fi_acd_0 + g_xyyyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxyyzz_0[i] = g_yyyy_0_xxxyyzz_0[i] * fbe_0 - g_yyyy_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyyzz_1[i] * fi_acd_0 + g_xyyyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxyzzz_0[i] = g_yyyy_0_xxxyzzz_0[i] * fbe_0 - g_yyyy_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xyyyy_0_xxyzzz_1[i] * fi_acd_0 + g_xyyyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxxzzzz_0[i] = 3.0 * g_xxyy_0_xxxzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxxzzzz_1[i] * fz_be_0 + g_xxyyy_0_xxxzzzz_1[i] * wa_y[i];

        g_xxyyyy_0_xxyyyyy_0[i] = g_yyyy_0_xxyyyyy_0[i] * fbe_0 - g_yyyy_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyyyy_1[i] * fi_acd_0 + g_xyyyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xxyyyyz_0[i] = g_yyyy_0_xxyyyyz_0[i] * fbe_0 - g_yyyy_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyyyz_1[i] * fi_acd_0 + g_xyyyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xxyyyzz_0[i] = g_yyyy_0_xxyyyzz_0[i] * fbe_0 - g_yyyy_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyyzz_1[i] * fi_acd_0 + g_xyyyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxyyzzz_0[i] = g_yyyy_0_xxyyzzz_0[i] * fbe_0 - g_yyyy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyyzzz_1[i] * fi_acd_0 + g_xyyyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxyzzzz_0[i] = g_yyyy_0_xxyzzzz_0[i] * fbe_0 - g_yyyy_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xyyyy_0_xyzzzz_1[i] * fi_acd_0 + g_xyyyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xxzzzzz_0[i] = 3.0 * g_xxyy_0_xxzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xxzzzzz_1[i] * fz_be_0 + g_xxyyy_0_xxzzzzz_1[i] * wa_y[i];

        g_xxyyyy_0_xyyyyyy_0[i] = g_yyyy_0_xyyyyyy_0[i] * fbe_0 - g_yyyy_0_xyyyyyy_1[i] * fz_be_0 + g_xyyyy_0_yyyyyy_1[i] * fi_acd_0 + g_xyyyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_xyyyyyz_0[i] = g_yyyy_0_xyyyyyz_0[i] * fbe_0 - g_yyyy_0_xyyyyyz_1[i] * fz_be_0 + g_xyyyy_0_yyyyyz_1[i] * fi_acd_0 + g_xyyyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_xyyyyzz_0[i] = g_yyyy_0_xyyyyzz_0[i] * fbe_0 - g_yyyy_0_xyyyyzz_1[i] * fz_be_0 + g_xyyyy_0_yyyyzz_1[i] * fi_acd_0 + g_xyyyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_xyyyzzz_0[i] = g_yyyy_0_xyyyzzz_0[i] * fbe_0 - g_yyyy_0_xyyyzzz_1[i] * fz_be_0 + g_xyyyy_0_yyyzzz_1[i] * fi_acd_0 + g_xyyyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xyyzzzz_0[i] = g_yyyy_0_xyyzzzz_0[i] * fbe_0 - g_yyyy_0_xyyzzzz_1[i] * fz_be_0 + g_xyyyy_0_yyzzzz_1[i] * fi_acd_0 + g_xyyyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xyzzzzz_0[i] = g_yyyy_0_xyzzzzz_0[i] * fbe_0 - g_yyyy_0_xyzzzzz_1[i] * fz_be_0 + g_xyyyy_0_yzzzzz_1[i] * fi_acd_0 + g_xyyyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_xzzzzzz_0[i] = 3.0 * g_xxyy_0_xzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyy_0_xzzzzzz_1[i] * fz_be_0 + g_xxyyy_0_xzzzzzz_1[i] * wa_y[i];

        g_xxyyyy_0_yyyyyyy_0[i] = g_yyyy_0_yyyyyyy_0[i] * fbe_0 - g_yyyy_0_yyyyyyy_1[i] * fz_be_0 + g_xyyyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xxyyyy_0_yyyyyyz_0[i] = g_yyyy_0_yyyyyyz_0[i] * fbe_0 - g_yyyy_0_yyyyyyz_1[i] * fz_be_0 + g_xyyyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xxyyyy_0_yyyyyzz_0[i] = g_yyyy_0_yyyyyzz_0[i] * fbe_0 - g_yyyy_0_yyyyyzz_1[i] * fz_be_0 + g_xyyyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xxyyyy_0_yyyyzzz_0[i] = g_yyyy_0_yyyyzzz_0[i] * fbe_0 - g_yyyy_0_yyyyzzz_1[i] * fz_be_0 + g_xyyyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xxyyyy_0_yyyzzzz_0[i] = g_yyyy_0_yyyzzzz_0[i] * fbe_0 - g_yyyy_0_yyyzzzz_1[i] * fz_be_0 + g_xyyyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_yyzzzzz_0[i] = g_yyyy_0_yyzzzzz_0[i] * fbe_0 - g_yyyy_0_yyzzzzz_1[i] * fz_be_0 + g_xyyyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_yzzzzzz_0[i] = g_yyyy_0_yzzzzzz_0[i] * fbe_0 - g_yyyy_0_yzzzzzz_1[i] * fz_be_0 + g_xyyyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xxyyyy_0_zzzzzzz_0[i] = g_yyyy_0_zzzzzzz_0[i] * fbe_0 - g_yyyy_0_zzzzzzz_1[i] * fz_be_0 + g_xyyyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 396-432 components of targeted buffer : ISK

    auto g_xxyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 396);

    auto g_xxyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 397);

    auto g_xxyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 398);

    auto g_xxyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 399);

    auto g_xxyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 400);

    auto g_xxyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 401);

    auto g_xxyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 402);

    auto g_xxyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 403);

    auto g_xxyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 404);

    auto g_xxyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 405);

    auto g_xxyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 406);

    auto g_xxyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 407);

    auto g_xxyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 408);

    auto g_xxyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 409);

    auto g_xxyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 410);

    auto g_xxyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 411);

    auto g_xxyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 412);

    auto g_xxyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 413);

    auto g_xxyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 414);

    auto g_xxyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 415);

    auto g_xxyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 416);

    auto g_xxyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 417);

    auto g_xxyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 418);

    auto g_xxyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 419);

    auto g_xxyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 420);

    auto g_xxyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 421);

    auto g_xxyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 422);

    auto g_xxyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 423);

    auto g_xxyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 424);

    auto g_xxyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 425);

    auto g_xxyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 426);

    auto g_xxyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 427);

    auto g_xxyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 428);

    auto g_xxyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 429);

    auto g_xxyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 430);

    auto g_xxyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 431);

    #pragma omp simd aligned(g_xxyyy_0_xxxxxx_1, g_xxyyy_0_xxxxxxx_1, g_xxyyy_0_xxxxxxy_1, g_xxyyy_0_xxxxxxz_1, g_xxyyy_0_xxxxxy_1, g_xxyyy_0_xxxxxyy_1, g_xxyyy_0_xxxxxyz_1, g_xxyyy_0_xxxxxz_1, g_xxyyy_0_xxxxxzz_1, g_xxyyy_0_xxxxyy_1, g_xxyyy_0_xxxxyyy_1, g_xxyyy_0_xxxxyyz_1, g_xxyyy_0_xxxxyz_1, g_xxyyy_0_xxxxyzz_1, g_xxyyy_0_xxxxzz_1, g_xxyyy_0_xxxxzzz_1, g_xxyyy_0_xxxyyy_1, g_xxyyy_0_xxxyyyy_1, g_xxyyy_0_xxxyyyz_1, g_xxyyy_0_xxxyyz_1, g_xxyyy_0_xxxyyzz_1, g_xxyyy_0_xxxyzz_1, g_xxyyy_0_xxxyzzz_1, g_xxyyy_0_xxxzzz_1, g_xxyyy_0_xxxzzzz_1, g_xxyyy_0_xxyyyy_1, g_xxyyy_0_xxyyyyy_1, g_xxyyy_0_xxyyyyz_1, g_xxyyy_0_xxyyyz_1, g_xxyyy_0_xxyyyzz_1, g_xxyyy_0_xxyyzz_1, g_xxyyy_0_xxyyzzz_1, g_xxyyy_0_xxyzzz_1, g_xxyyy_0_xxyzzzz_1, g_xxyyy_0_xxzzzz_1, g_xxyyy_0_xxzzzzz_1, g_xxyyy_0_xyyyyy_1, g_xxyyy_0_xyyyyyy_1, g_xxyyy_0_xyyyyyz_1, g_xxyyy_0_xyyyyz_1, g_xxyyy_0_xyyyyzz_1, g_xxyyy_0_xyyyzz_1, g_xxyyy_0_xyyyzzz_1, g_xxyyy_0_xyyzzz_1, g_xxyyy_0_xyyzzzz_1, g_xxyyy_0_xyzzzz_1, g_xxyyy_0_xyzzzzz_1, g_xxyyy_0_xzzzzz_1, g_xxyyy_0_xzzzzzz_1, g_xxyyy_0_yyyyyy_1, g_xxyyy_0_yyyyyyy_1, g_xxyyy_0_yyyyyyz_1, g_xxyyy_0_yyyyyz_1, g_xxyyy_0_yyyyyzz_1, g_xxyyy_0_yyyyzz_1, g_xxyyy_0_yyyyzzz_1, g_xxyyy_0_yyyzzz_1, g_xxyyy_0_yyyzzzz_1, g_xxyyy_0_yyzzzz_1, g_xxyyy_0_yyzzzzz_1, g_xxyyy_0_yzzzzz_1, g_xxyyy_0_yzzzzzz_1, g_xxyyy_0_zzzzzz_1, g_xxyyy_0_zzzzzzz_1, g_xxyyyz_0_xxxxxxx_0, g_xxyyyz_0_xxxxxxy_0, g_xxyyyz_0_xxxxxxz_0, g_xxyyyz_0_xxxxxyy_0, g_xxyyyz_0_xxxxxyz_0, g_xxyyyz_0_xxxxxzz_0, g_xxyyyz_0_xxxxyyy_0, g_xxyyyz_0_xxxxyyz_0, g_xxyyyz_0_xxxxyzz_0, g_xxyyyz_0_xxxxzzz_0, g_xxyyyz_0_xxxyyyy_0, g_xxyyyz_0_xxxyyyz_0, g_xxyyyz_0_xxxyyzz_0, g_xxyyyz_0_xxxyzzz_0, g_xxyyyz_0_xxxzzzz_0, g_xxyyyz_0_xxyyyyy_0, g_xxyyyz_0_xxyyyyz_0, g_xxyyyz_0_xxyyyzz_0, g_xxyyyz_0_xxyyzzz_0, g_xxyyyz_0_xxyzzzz_0, g_xxyyyz_0_xxzzzzz_0, g_xxyyyz_0_xyyyyyy_0, g_xxyyyz_0_xyyyyyz_0, g_xxyyyz_0_xyyyyzz_0, g_xxyyyz_0_xyyyzzz_0, g_xxyyyz_0_xyyzzzz_0, g_xxyyyz_0_xyzzzzz_0, g_xxyyyz_0_xzzzzzz_0, g_xxyyyz_0_yyyyyyy_0, g_xxyyyz_0_yyyyyyz_0, g_xxyyyz_0_yyyyyzz_0, g_xxyyyz_0_yyyyzzz_0, g_xxyyyz_0_yyyzzzz_0, g_xxyyyz_0_yyzzzzz_0, g_xxyyyz_0_yzzzzzz_0, g_xxyyyz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyz_0_xxxxxxx_0[i] = g_xxyyy_0_xxxxxxx_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxxy_0[i] = g_xxyyy_0_xxxxxxy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxxz_0[i] = g_xxyyy_0_xxxxxx_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxxz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxyy_0[i] = g_xxyyy_0_xxxxxyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxyz_0[i] = g_xxyyy_0_xxxxxy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxxzz_0[i] = 2.0 * g_xxyyy_0_xxxxxz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxxzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxyyy_0[i] = g_xxyyy_0_xxxxyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxyyz_0[i] = g_xxyyy_0_xxxxyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxyzz_0[i] = 2.0 * g_xxyyy_0_xxxxyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxxzzz_0[i] = 3.0 * g_xxyyy_0_xxxxzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxxzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyyyy_0[i] = g_xxyyy_0_xxxyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyyyz_0[i] = g_xxyyy_0_xxxyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyyzz_0[i] = 2.0 * g_xxyyy_0_xxxyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxyzzz_0[i] = 3.0 * g_xxyyy_0_xxxyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxxzzzz_0[i] = 4.0 * g_xxyyy_0_xxxzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxxzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyyyy_0[i] = g_xxyyy_0_xxyyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyyyz_0[i] = g_xxyyy_0_xxyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyyzz_0[i] = 2.0 * g_xxyyy_0_xxyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyyzzz_0[i] = 3.0 * g_xxyyy_0_xxyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxyzzzz_0[i] = 4.0 * g_xxyyy_0_xxyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxyzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xxzzzzz_0[i] = 5.0 * g_xxyyy_0_xxzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xxzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyyyy_0[i] = g_xxyyy_0_xyyyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyyyz_0[i] = g_xxyyy_0_xyyyyy_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyyzz_0[i] = 2.0 * g_xxyyy_0_xyyyyz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyyzzz_0[i] = 3.0 * g_xxyyy_0_xyyyzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyyzzzz_0[i] = 4.0 * g_xxyyy_0_xyyzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyyzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xyzzzzz_0[i] = 5.0 * g_xxyyy_0_xyzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xyzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_xzzzzzz_0[i] = 6.0 * g_xxyyy_0_xzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_xzzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyyyy_0[i] = g_xxyyy_0_yyyyyyy_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyyyz_0[i] = g_xxyyy_0_yyyyyy_1[i] * fi_acd_0 + g_xxyyy_0_yyyyyyz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyyzz_0[i] = 2.0 * g_xxyyy_0_yyyyyz_1[i] * fi_acd_0 + g_xxyyy_0_yyyyyzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyyzzz_0[i] = 3.0 * g_xxyyy_0_yyyyzz_1[i] * fi_acd_0 + g_xxyyy_0_yyyyzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyyzzzz_0[i] = 4.0 * g_xxyyy_0_yyyzzz_1[i] * fi_acd_0 + g_xxyyy_0_yyyzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yyzzzzz_0[i] = 5.0 * g_xxyyy_0_yyzzzz_1[i] * fi_acd_0 + g_xxyyy_0_yyzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_yzzzzzz_0[i] = 6.0 * g_xxyyy_0_yzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_yzzzzzz_1[i] * wa_z[i];

        g_xxyyyz_0_zzzzzzz_0[i] = 7.0 * g_xxyyy_0_zzzzzz_1[i] * fi_acd_0 + g_xxyyy_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 432-468 components of targeted buffer : ISK

    auto g_xxyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 432);

    auto g_xxyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 433);

    auto g_xxyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 434);

    auto g_xxyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 435);

    auto g_xxyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 436);

    auto g_xxyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 437);

    auto g_xxyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 438);

    auto g_xxyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 439);

    auto g_xxyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 440);

    auto g_xxyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 441);

    auto g_xxyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 442);

    auto g_xxyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 443);

    auto g_xxyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 444);

    auto g_xxyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 445);

    auto g_xxyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 446);

    auto g_xxyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 447);

    auto g_xxyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 448);

    auto g_xxyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 449);

    auto g_xxyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 450);

    auto g_xxyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 451);

    auto g_xxyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 452);

    auto g_xxyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 453);

    auto g_xxyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 454);

    auto g_xxyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 455);

    auto g_xxyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 456);

    auto g_xxyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 457);

    auto g_xxyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 458);

    auto g_xxyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 459);

    auto g_xxyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 460);

    auto g_xxyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 461);

    auto g_xxyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 462);

    auto g_xxyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 463);

    auto g_xxyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 464);

    auto g_xxyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 465);

    auto g_xxyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 466);

    auto g_xxyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 467);

    #pragma omp simd aligned(g_xxyy_0_xxxxxxy_0, g_xxyy_0_xxxxxxy_1, g_xxyy_0_xxxxxyy_0, g_xxyy_0_xxxxxyy_1, g_xxyy_0_xxxxyyy_0, g_xxyy_0_xxxxyyy_1, g_xxyy_0_xxxyyyy_0, g_xxyy_0_xxxyyyy_1, g_xxyy_0_xxyyyyy_0, g_xxyy_0_xxyyyyy_1, g_xxyy_0_xyyyyyy_0, g_xxyy_0_xyyyyyy_1, g_xxyyz_0_xxxxxxy_1, g_xxyyz_0_xxxxxyy_1, g_xxyyz_0_xxxxyyy_1, g_xxyyz_0_xxxyyyy_1, g_xxyyz_0_xxyyyyy_1, g_xxyyz_0_xyyyyyy_1, g_xxyyzz_0_xxxxxxx_0, g_xxyyzz_0_xxxxxxy_0, g_xxyyzz_0_xxxxxxz_0, g_xxyyzz_0_xxxxxyy_0, g_xxyyzz_0_xxxxxyz_0, g_xxyyzz_0_xxxxxzz_0, g_xxyyzz_0_xxxxyyy_0, g_xxyyzz_0_xxxxyyz_0, g_xxyyzz_0_xxxxyzz_0, g_xxyyzz_0_xxxxzzz_0, g_xxyyzz_0_xxxyyyy_0, g_xxyyzz_0_xxxyyyz_0, g_xxyyzz_0_xxxyyzz_0, g_xxyyzz_0_xxxyzzz_0, g_xxyyzz_0_xxxzzzz_0, g_xxyyzz_0_xxyyyyy_0, g_xxyyzz_0_xxyyyyz_0, g_xxyyzz_0_xxyyyzz_0, g_xxyyzz_0_xxyyzzz_0, g_xxyyzz_0_xxyzzzz_0, g_xxyyzz_0_xxzzzzz_0, g_xxyyzz_0_xyyyyyy_0, g_xxyyzz_0_xyyyyyz_0, g_xxyyzz_0_xyyyyzz_0, g_xxyyzz_0_xyyyzzz_0, g_xxyyzz_0_xyyzzzz_0, g_xxyyzz_0_xyzzzzz_0, g_xxyyzz_0_xzzzzzz_0, g_xxyyzz_0_yyyyyyy_0, g_xxyyzz_0_yyyyyyz_0, g_xxyyzz_0_yyyyyzz_0, g_xxyyzz_0_yyyyzzz_0, g_xxyyzz_0_yyyzzzz_0, g_xxyyzz_0_yyzzzzz_0, g_xxyyzz_0_yzzzzzz_0, g_xxyyzz_0_zzzzzzz_0, g_xxyzz_0_xxxxxxx_1, g_xxyzz_0_xxxxxxz_1, g_xxyzz_0_xxxxxzz_1, g_xxyzz_0_xxxxzzz_1, g_xxyzz_0_xxxzzzz_1, g_xxyzz_0_xxzzzzz_1, g_xxyzz_0_xzzzzzz_1, g_xxzz_0_xxxxxxx_0, g_xxzz_0_xxxxxxx_1, g_xxzz_0_xxxxxxz_0, g_xxzz_0_xxxxxxz_1, g_xxzz_0_xxxxxzz_0, g_xxzz_0_xxxxxzz_1, g_xxzz_0_xxxxzzz_0, g_xxzz_0_xxxxzzz_1, g_xxzz_0_xxxzzzz_0, g_xxzz_0_xxxzzzz_1, g_xxzz_0_xxzzzzz_0, g_xxzz_0_xxzzzzz_1, g_xxzz_0_xzzzzzz_0, g_xxzz_0_xzzzzzz_1, g_xyyzz_0_xxxxxyz_1, g_xyyzz_0_xxxxyyz_1, g_xyyzz_0_xxxxyz_1, g_xyyzz_0_xxxxyzz_1, g_xyyzz_0_xxxyyyz_1, g_xyyzz_0_xxxyyz_1, g_xyyzz_0_xxxyyzz_1, g_xyyzz_0_xxxyzz_1, g_xyyzz_0_xxxyzzz_1, g_xyyzz_0_xxyyyyz_1, g_xyyzz_0_xxyyyz_1, g_xyyzz_0_xxyyyzz_1, g_xyyzz_0_xxyyzz_1, g_xyyzz_0_xxyyzzz_1, g_xyyzz_0_xxyzzz_1, g_xyyzz_0_xxyzzzz_1, g_xyyzz_0_xyyyyyz_1, g_xyyzz_0_xyyyyz_1, g_xyyzz_0_xyyyyzz_1, g_xyyzz_0_xyyyzz_1, g_xyyzz_0_xyyyzzz_1, g_xyyzz_0_xyyzzz_1, g_xyyzz_0_xyyzzzz_1, g_xyyzz_0_xyzzzz_1, g_xyyzz_0_xyzzzzz_1, g_xyyzz_0_yyyyyyy_1, g_xyyzz_0_yyyyyyz_1, g_xyyzz_0_yyyyyz_1, g_xyyzz_0_yyyyyzz_1, g_xyyzz_0_yyyyzz_1, g_xyyzz_0_yyyyzzz_1, g_xyyzz_0_yyyzzz_1, g_xyyzz_0_yyyzzzz_1, g_xyyzz_0_yyzzzz_1, g_xyyzz_0_yyzzzzz_1, g_xyyzz_0_yzzzzz_1, g_xyyzz_0_yzzzzzz_1, g_xyyzz_0_zzzzzzz_1, g_yyzz_0_xxxxxyz_0, g_yyzz_0_xxxxxyz_1, g_yyzz_0_xxxxyyz_0, g_yyzz_0_xxxxyyz_1, g_yyzz_0_xxxxyzz_0, g_yyzz_0_xxxxyzz_1, g_yyzz_0_xxxyyyz_0, g_yyzz_0_xxxyyyz_1, g_yyzz_0_xxxyyzz_0, g_yyzz_0_xxxyyzz_1, g_yyzz_0_xxxyzzz_0, g_yyzz_0_xxxyzzz_1, g_yyzz_0_xxyyyyz_0, g_yyzz_0_xxyyyyz_1, g_yyzz_0_xxyyyzz_0, g_yyzz_0_xxyyyzz_1, g_yyzz_0_xxyyzzz_0, g_yyzz_0_xxyyzzz_1, g_yyzz_0_xxyzzzz_0, g_yyzz_0_xxyzzzz_1, g_yyzz_0_xyyyyyz_0, g_yyzz_0_xyyyyyz_1, g_yyzz_0_xyyyyzz_0, g_yyzz_0_xyyyyzz_1, g_yyzz_0_xyyyzzz_0, g_yyzz_0_xyyyzzz_1, g_yyzz_0_xyyzzzz_0, g_yyzz_0_xyyzzzz_1, g_yyzz_0_xyzzzzz_0, g_yyzz_0_xyzzzzz_1, g_yyzz_0_yyyyyyy_0, g_yyzz_0_yyyyyyy_1, g_yyzz_0_yyyyyyz_0, g_yyzz_0_yyyyyyz_1, g_yyzz_0_yyyyyzz_0, g_yyzz_0_yyyyyzz_1, g_yyzz_0_yyyyzzz_0, g_yyzz_0_yyyyzzz_1, g_yyzz_0_yyyzzzz_0, g_yyzz_0_yyyzzzz_1, g_yyzz_0_yyzzzzz_0, g_yyzz_0_yyzzzzz_1, g_yyzz_0_yzzzzzz_0, g_yyzz_0_yzzzzzz_1, g_yyzz_0_zzzzzzz_0, g_yyzz_0_zzzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyzz_0_xxxxxxx_0[i] = g_xxzz_0_xxxxxxx_0[i] * fbe_0 - g_xxzz_0_xxxxxxx_1[i] * fz_be_0 + g_xxyzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxyyzz_0_xxxxxxy_0[i] = g_xxyy_0_xxxxxxy_0[i] * fbe_0 - g_xxyy_0_xxxxxxy_1[i] * fz_be_0 + g_xxyyz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxxxxz_0[i] = g_xxzz_0_xxxxxxz_0[i] * fbe_0 - g_xxzz_0_xxxxxxz_1[i] * fz_be_0 + g_xxyzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxyyzz_0_xxxxxyy_0[i] = g_xxyy_0_xxxxxyy_0[i] * fbe_0 - g_xxyy_0_xxxxxyy_1[i] * fz_be_0 + g_xxyyz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxxxyz_0[i] = g_yyzz_0_xxxxxyz_0[i] * fbe_0 - g_yyzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xyyzz_0_xxxxyz_1[i] * fi_acd_0 + g_xyyzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxxxzz_0[i] = g_xxzz_0_xxxxxzz_0[i] * fbe_0 - g_xxzz_0_xxxxxzz_1[i] * fz_be_0 + g_xxyzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxyyzz_0_xxxxyyy_0[i] = g_xxyy_0_xxxxyyy_0[i] * fbe_0 - g_xxyy_0_xxxxyyy_1[i] * fz_be_0 + g_xxyyz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxxyyz_0[i] = g_yyzz_0_xxxxyyz_0[i] * fbe_0 - g_yyzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xyyzz_0_xxxyyz_1[i] * fi_acd_0 + g_xyyzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxxyzz_0[i] = g_yyzz_0_xxxxyzz_0[i] * fbe_0 - g_yyzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xyyzz_0_xxxyzz_1[i] * fi_acd_0 + g_xyyzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxxzzz_0[i] = g_xxzz_0_xxxxzzz_0[i] * fbe_0 - g_xxzz_0_xxxxzzz_1[i] * fz_be_0 + g_xxyzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxyyzz_0_xxxyyyy_0[i] = g_xxyy_0_xxxyyyy_0[i] * fbe_0 - g_xxyy_0_xxxyyyy_1[i] * fz_be_0 + g_xxyyz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxxyyyz_0[i] = g_yyzz_0_xxxyyyz_0[i] * fbe_0 - g_yyzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xyyzz_0_xxyyyz_1[i] * fi_acd_0 + g_xyyzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxyyzz_0[i] = g_yyzz_0_xxxyyzz_0[i] * fbe_0 - g_yyzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xyyzz_0_xxyyzz_1[i] * fi_acd_0 + g_xyyzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxyzzz_0[i] = g_yyzz_0_xxxyzzz_0[i] * fbe_0 - g_yyzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xyyzz_0_xxyzzz_1[i] * fi_acd_0 + g_xyyzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxxzzzz_0[i] = g_xxzz_0_xxxzzzz_0[i] * fbe_0 - g_xxzz_0_xxxzzzz_1[i] * fz_be_0 + g_xxyzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxyyzz_0_xxyyyyy_0[i] = g_xxyy_0_xxyyyyy_0[i] * fbe_0 - g_xxyy_0_xxyyyyy_1[i] * fz_be_0 + g_xxyyz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xxyyyyz_0[i] = g_yyzz_0_xxyyyyz_0[i] * fbe_0 - g_yyzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyyyyz_1[i] * fi_acd_0 + g_xyyzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xxyyyzz_0[i] = g_yyzz_0_xxyyyzz_0[i] * fbe_0 - g_yyzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyyyzz_1[i] * fi_acd_0 + g_xyyzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxyyzzz_0[i] = g_yyzz_0_xxyyzzz_0[i] * fbe_0 - g_yyzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyyzzz_1[i] * fi_acd_0 + g_xyyzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxyzzzz_0[i] = g_yyzz_0_xxyzzzz_0[i] * fbe_0 - g_yyzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xyyzz_0_xyzzzz_1[i] * fi_acd_0 + g_xyyzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xxzzzzz_0[i] = g_xxzz_0_xxzzzzz_0[i] * fbe_0 - g_xxzz_0_xxzzzzz_1[i] * fz_be_0 + g_xxyzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxyyzz_0_xyyyyyy_0[i] = g_xxyy_0_xyyyyyy_0[i] * fbe_0 - g_xxyy_0_xyyyyyy_1[i] * fz_be_0 + g_xxyyz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxyyzz_0_xyyyyyz_0[i] = g_yyzz_0_xyyyyyz_0[i] * fbe_0 - g_yyzz_0_xyyyyyz_1[i] * fz_be_0 + g_xyyzz_0_yyyyyz_1[i] * fi_acd_0 + g_xyyzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_xyyyyzz_0[i] = g_yyzz_0_xyyyyzz_0[i] * fbe_0 - g_yyzz_0_xyyyyzz_1[i] * fz_be_0 + g_xyyzz_0_yyyyzz_1[i] * fi_acd_0 + g_xyyzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_xyyyzzz_0[i] = g_yyzz_0_xyyyzzz_0[i] * fbe_0 - g_yyzz_0_xyyyzzz_1[i] * fz_be_0 + g_xyyzz_0_yyyzzz_1[i] * fi_acd_0 + g_xyyzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xyyzzzz_0[i] = g_yyzz_0_xyyzzzz_0[i] * fbe_0 - g_yyzz_0_xyyzzzz_1[i] * fz_be_0 + g_xyyzz_0_yyzzzz_1[i] * fi_acd_0 + g_xyyzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xyzzzzz_0[i] = g_yyzz_0_xyzzzzz_0[i] * fbe_0 - g_yyzz_0_xyzzzzz_1[i] * fz_be_0 + g_xyyzz_0_yzzzzz_1[i] * fi_acd_0 + g_xyyzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_xzzzzzz_0[i] = g_xxzz_0_xzzzzzz_0[i] * fbe_0 - g_xxzz_0_xzzzzzz_1[i] * fz_be_0 + g_xxyzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxyyzz_0_yyyyyyy_0[i] = g_yyzz_0_yyyyyyy_0[i] * fbe_0 - g_yyzz_0_yyyyyyy_1[i] * fz_be_0 + g_xyyzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxyyzz_0_yyyyyyz_0[i] = g_yyzz_0_yyyyyyz_0[i] * fbe_0 - g_yyzz_0_yyyyyyz_1[i] * fz_be_0 + g_xyyzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxyyzz_0_yyyyyzz_0[i] = g_yyzz_0_yyyyyzz_0[i] * fbe_0 - g_yyzz_0_yyyyyzz_1[i] * fz_be_0 + g_xyyzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxyyzz_0_yyyyzzz_0[i] = g_yyzz_0_yyyyzzz_0[i] * fbe_0 - g_yyzz_0_yyyyzzz_1[i] * fz_be_0 + g_xyyzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxyyzz_0_yyyzzzz_0[i] = g_yyzz_0_yyyzzzz_0[i] * fbe_0 - g_yyzz_0_yyyzzzz_1[i] * fz_be_0 + g_xyyzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_yyzzzzz_0[i] = g_yyzz_0_yyzzzzz_0[i] * fbe_0 - g_yyzz_0_yyzzzzz_1[i] * fz_be_0 + g_xyyzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_yzzzzzz_0[i] = g_yyzz_0_yzzzzzz_0[i] * fbe_0 - g_yyzz_0_yzzzzzz_1[i] * fz_be_0 + g_xyyzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxyyzz_0_zzzzzzz_0[i] = g_yyzz_0_zzzzzzz_0[i] * fbe_0 - g_yyzz_0_zzzzzzz_1[i] * fz_be_0 + g_xyyzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 468-504 components of targeted buffer : ISK

    auto g_xxyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 468);

    auto g_xxyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 469);

    auto g_xxyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 470);

    auto g_xxyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 471);

    auto g_xxyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 472);

    auto g_xxyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 473);

    auto g_xxyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 474);

    auto g_xxyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 475);

    auto g_xxyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 476);

    auto g_xxyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 477);

    auto g_xxyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 478);

    auto g_xxyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 479);

    auto g_xxyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 480);

    auto g_xxyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 481);

    auto g_xxyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 482);

    auto g_xxyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 483);

    auto g_xxyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 484);

    auto g_xxyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 485);

    auto g_xxyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 486);

    auto g_xxyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 487);

    auto g_xxyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 488);

    auto g_xxyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 489);

    auto g_xxyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 490);

    auto g_xxyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 491);

    auto g_xxyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 492);

    auto g_xxyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 493);

    auto g_xxyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 494);

    auto g_xxyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 495);

    auto g_xxyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 496);

    auto g_xxyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 497);

    auto g_xxyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 498);

    auto g_xxyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 499);

    auto g_xxyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 500);

    auto g_xxyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 501);

    auto g_xxyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 502);

    auto g_xxyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 503);

    #pragma omp simd aligned(g_xxyzzz_0_xxxxxxx_0, g_xxyzzz_0_xxxxxxy_0, g_xxyzzz_0_xxxxxxz_0, g_xxyzzz_0_xxxxxyy_0, g_xxyzzz_0_xxxxxyz_0, g_xxyzzz_0_xxxxxzz_0, g_xxyzzz_0_xxxxyyy_0, g_xxyzzz_0_xxxxyyz_0, g_xxyzzz_0_xxxxyzz_0, g_xxyzzz_0_xxxxzzz_0, g_xxyzzz_0_xxxyyyy_0, g_xxyzzz_0_xxxyyyz_0, g_xxyzzz_0_xxxyyzz_0, g_xxyzzz_0_xxxyzzz_0, g_xxyzzz_0_xxxzzzz_0, g_xxyzzz_0_xxyyyyy_0, g_xxyzzz_0_xxyyyyz_0, g_xxyzzz_0_xxyyyzz_0, g_xxyzzz_0_xxyyzzz_0, g_xxyzzz_0_xxyzzzz_0, g_xxyzzz_0_xxzzzzz_0, g_xxyzzz_0_xyyyyyy_0, g_xxyzzz_0_xyyyyyz_0, g_xxyzzz_0_xyyyyzz_0, g_xxyzzz_0_xyyyzzz_0, g_xxyzzz_0_xyyzzzz_0, g_xxyzzz_0_xyzzzzz_0, g_xxyzzz_0_xzzzzzz_0, g_xxyzzz_0_yyyyyyy_0, g_xxyzzz_0_yyyyyyz_0, g_xxyzzz_0_yyyyyzz_0, g_xxyzzz_0_yyyyzzz_0, g_xxyzzz_0_yyyzzzz_0, g_xxyzzz_0_yyzzzzz_0, g_xxyzzz_0_yzzzzzz_0, g_xxyzzz_0_zzzzzzz_0, g_xxzzz_0_xxxxxx_1, g_xxzzz_0_xxxxxxx_1, g_xxzzz_0_xxxxxxy_1, g_xxzzz_0_xxxxxxz_1, g_xxzzz_0_xxxxxy_1, g_xxzzz_0_xxxxxyy_1, g_xxzzz_0_xxxxxyz_1, g_xxzzz_0_xxxxxz_1, g_xxzzz_0_xxxxxzz_1, g_xxzzz_0_xxxxyy_1, g_xxzzz_0_xxxxyyy_1, g_xxzzz_0_xxxxyyz_1, g_xxzzz_0_xxxxyz_1, g_xxzzz_0_xxxxyzz_1, g_xxzzz_0_xxxxzz_1, g_xxzzz_0_xxxxzzz_1, g_xxzzz_0_xxxyyy_1, g_xxzzz_0_xxxyyyy_1, g_xxzzz_0_xxxyyyz_1, g_xxzzz_0_xxxyyz_1, g_xxzzz_0_xxxyyzz_1, g_xxzzz_0_xxxyzz_1, g_xxzzz_0_xxxyzzz_1, g_xxzzz_0_xxxzzz_1, g_xxzzz_0_xxxzzzz_1, g_xxzzz_0_xxyyyy_1, g_xxzzz_0_xxyyyyy_1, g_xxzzz_0_xxyyyyz_1, g_xxzzz_0_xxyyyz_1, g_xxzzz_0_xxyyyzz_1, g_xxzzz_0_xxyyzz_1, g_xxzzz_0_xxyyzzz_1, g_xxzzz_0_xxyzzz_1, g_xxzzz_0_xxyzzzz_1, g_xxzzz_0_xxzzzz_1, g_xxzzz_0_xxzzzzz_1, g_xxzzz_0_xyyyyy_1, g_xxzzz_0_xyyyyyy_1, g_xxzzz_0_xyyyyyz_1, g_xxzzz_0_xyyyyz_1, g_xxzzz_0_xyyyyzz_1, g_xxzzz_0_xyyyzz_1, g_xxzzz_0_xyyyzzz_1, g_xxzzz_0_xyyzzz_1, g_xxzzz_0_xyyzzzz_1, g_xxzzz_0_xyzzzz_1, g_xxzzz_0_xyzzzzz_1, g_xxzzz_0_xzzzzz_1, g_xxzzz_0_xzzzzzz_1, g_xxzzz_0_yyyyyy_1, g_xxzzz_0_yyyyyyy_1, g_xxzzz_0_yyyyyyz_1, g_xxzzz_0_yyyyyz_1, g_xxzzz_0_yyyyyzz_1, g_xxzzz_0_yyyyzz_1, g_xxzzz_0_yyyyzzz_1, g_xxzzz_0_yyyzzz_1, g_xxzzz_0_yyyzzzz_1, g_xxzzz_0_yyzzzz_1, g_xxzzz_0_yyzzzzz_1, g_xxzzz_0_yzzzzz_1, g_xxzzz_0_yzzzzzz_1, g_xxzzz_0_zzzzzz_1, g_xxzzz_0_zzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzz_0_xxxxxxx_0[i] = g_xxzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxxy_0[i] = g_xxzzz_0_xxxxxx_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxxy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxxz_0[i] = g_xxzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxyy_0[i] = 2.0 * g_xxzzz_0_xxxxxy_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxyz_0[i] = g_xxzzz_0_xxxxxz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxxyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxxzz_0[i] = g_xxzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxyyy_0[i] = 3.0 * g_xxzzz_0_xxxxyy_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxyyz_0[i] = 2.0 * g_xxzzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxyzz_0[i] = g_xxzzz_0_xxxxzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxxyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxxzzz_0[i] = g_xxzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyyyy_0[i] = 4.0 * g_xxzzz_0_xxxyyy_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyyyz_0[i] = 3.0 * g_xxzzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyyzz_0[i] = 2.0 * g_xxzzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxyzzz_0[i] = g_xxzzz_0_xxxzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxxyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxxzzzz_0[i] = g_xxzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyyyy_0[i] = 5.0 * g_xxzzz_0_xxyyyy_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyyyz_0[i] = 4.0 * g_xxzzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyyzz_0[i] = 3.0 * g_xxzzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyyzzz_0[i] = 2.0 * g_xxzzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxyzzzz_0[i] = g_xxzzz_0_xxzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xxyzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xxzzzzz_0[i] = g_xxzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyyyy_0[i] = 6.0 * g_xxzzz_0_xyyyyy_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyyyz_0[i] = 5.0 * g_xxzzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyyzz_0[i] = 4.0 * g_xxzzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyyzzz_0[i] = 3.0 * g_xxzzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyyzzzz_0[i] = 2.0 * g_xxzzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyyzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xyzzzzz_0[i] = g_xxzzz_0_xzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_xyzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_xzzzzzz_0[i] = g_xxzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyyyy_0[i] = 7.0 * g_xxzzz_0_yyyyyy_1[i] * fi_acd_0 + g_xxzzz_0_yyyyyyy_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyyyz_0[i] = 6.0 * g_xxzzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxzzz_0_yyyyyyz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyyzz_0[i] = 5.0 * g_xxzzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxzzz_0_yyyyyzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyyzzz_0[i] = 4.0 * g_xxzzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxzzz_0_yyyyzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyyzzzz_0[i] = 3.0 * g_xxzzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxzzz_0_yyyzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yyzzzzz_0[i] = 2.0 * g_xxzzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_yyzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_yzzzzzz_0[i] = g_xxzzz_0_zzzzzz_1[i] * fi_acd_0 + g_xxzzz_0_yzzzzzz_1[i] * wa_y[i];

        g_xxyzzz_0_zzzzzzz_0[i] = g_xxzzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 504-540 components of targeted buffer : ISK

    auto g_xxzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 504);

    auto g_xxzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 505);

    auto g_xxzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 506);

    auto g_xxzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 507);

    auto g_xxzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 508);

    auto g_xxzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 509);

    auto g_xxzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 510);

    auto g_xxzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 511);

    auto g_xxzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 512);

    auto g_xxzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 513);

    auto g_xxzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 514);

    auto g_xxzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 515);

    auto g_xxzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 516);

    auto g_xxzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 517);

    auto g_xxzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 518);

    auto g_xxzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 519);

    auto g_xxzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 520);

    auto g_xxzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 521);

    auto g_xxzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 522);

    auto g_xxzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 523);

    auto g_xxzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 524);

    auto g_xxzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 525);

    auto g_xxzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 526);

    auto g_xxzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 527);

    auto g_xxzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 528);

    auto g_xxzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 529);

    auto g_xxzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 530);

    auto g_xxzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 531);

    auto g_xxzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 532);

    auto g_xxzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 533);

    auto g_xxzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 534);

    auto g_xxzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 535);

    auto g_xxzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 536);

    auto g_xxzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 537);

    auto g_xxzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 538);

    auto g_xxzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 539);

    #pragma omp simd aligned(g_xxzz_0_xxxxxxx_0, g_xxzz_0_xxxxxxx_1, g_xxzz_0_xxxxxxy_0, g_xxzz_0_xxxxxxy_1, g_xxzz_0_xxxxxyy_0, g_xxzz_0_xxxxxyy_1, g_xxzz_0_xxxxyyy_0, g_xxzz_0_xxxxyyy_1, g_xxzz_0_xxxyyyy_0, g_xxzz_0_xxxyyyy_1, g_xxzz_0_xxyyyyy_0, g_xxzz_0_xxyyyyy_1, g_xxzz_0_xyyyyyy_0, g_xxzz_0_xyyyyyy_1, g_xxzzz_0_xxxxxxx_1, g_xxzzz_0_xxxxxxy_1, g_xxzzz_0_xxxxxyy_1, g_xxzzz_0_xxxxyyy_1, g_xxzzz_0_xxxyyyy_1, g_xxzzz_0_xxyyyyy_1, g_xxzzz_0_xyyyyyy_1, g_xxzzzz_0_xxxxxxx_0, g_xxzzzz_0_xxxxxxy_0, g_xxzzzz_0_xxxxxxz_0, g_xxzzzz_0_xxxxxyy_0, g_xxzzzz_0_xxxxxyz_0, g_xxzzzz_0_xxxxxzz_0, g_xxzzzz_0_xxxxyyy_0, g_xxzzzz_0_xxxxyyz_0, g_xxzzzz_0_xxxxyzz_0, g_xxzzzz_0_xxxxzzz_0, g_xxzzzz_0_xxxyyyy_0, g_xxzzzz_0_xxxyyyz_0, g_xxzzzz_0_xxxyyzz_0, g_xxzzzz_0_xxxyzzz_0, g_xxzzzz_0_xxxzzzz_0, g_xxzzzz_0_xxyyyyy_0, g_xxzzzz_0_xxyyyyz_0, g_xxzzzz_0_xxyyyzz_0, g_xxzzzz_0_xxyyzzz_0, g_xxzzzz_0_xxyzzzz_0, g_xxzzzz_0_xxzzzzz_0, g_xxzzzz_0_xyyyyyy_0, g_xxzzzz_0_xyyyyyz_0, g_xxzzzz_0_xyyyyzz_0, g_xxzzzz_0_xyyyzzz_0, g_xxzzzz_0_xyyzzzz_0, g_xxzzzz_0_xyzzzzz_0, g_xxzzzz_0_xzzzzzz_0, g_xxzzzz_0_yyyyyyy_0, g_xxzzzz_0_yyyyyyz_0, g_xxzzzz_0_yyyyyzz_0, g_xxzzzz_0_yyyyzzz_0, g_xxzzzz_0_yyyzzzz_0, g_xxzzzz_0_yyzzzzz_0, g_xxzzzz_0_yzzzzzz_0, g_xxzzzz_0_zzzzzzz_0, g_xzzzz_0_xxxxxxz_1, g_xzzzz_0_xxxxxyz_1, g_xzzzz_0_xxxxxz_1, g_xzzzz_0_xxxxxzz_1, g_xzzzz_0_xxxxyyz_1, g_xzzzz_0_xxxxyz_1, g_xzzzz_0_xxxxyzz_1, g_xzzzz_0_xxxxzz_1, g_xzzzz_0_xxxxzzz_1, g_xzzzz_0_xxxyyyz_1, g_xzzzz_0_xxxyyz_1, g_xzzzz_0_xxxyyzz_1, g_xzzzz_0_xxxyzz_1, g_xzzzz_0_xxxyzzz_1, g_xzzzz_0_xxxzzz_1, g_xzzzz_0_xxxzzzz_1, g_xzzzz_0_xxyyyyz_1, g_xzzzz_0_xxyyyz_1, g_xzzzz_0_xxyyyzz_1, g_xzzzz_0_xxyyzz_1, g_xzzzz_0_xxyyzzz_1, g_xzzzz_0_xxyzzz_1, g_xzzzz_0_xxyzzzz_1, g_xzzzz_0_xxzzzz_1, g_xzzzz_0_xxzzzzz_1, g_xzzzz_0_xyyyyyz_1, g_xzzzz_0_xyyyyz_1, g_xzzzz_0_xyyyyzz_1, g_xzzzz_0_xyyyzz_1, g_xzzzz_0_xyyyzzz_1, g_xzzzz_0_xyyzzz_1, g_xzzzz_0_xyyzzzz_1, g_xzzzz_0_xyzzzz_1, g_xzzzz_0_xyzzzzz_1, g_xzzzz_0_xzzzzz_1, g_xzzzz_0_xzzzzzz_1, g_xzzzz_0_yyyyyyy_1, g_xzzzz_0_yyyyyyz_1, g_xzzzz_0_yyyyyz_1, g_xzzzz_0_yyyyyzz_1, g_xzzzz_0_yyyyzz_1, g_xzzzz_0_yyyyzzz_1, g_xzzzz_0_yyyzzz_1, g_xzzzz_0_yyyzzzz_1, g_xzzzz_0_yyzzzz_1, g_xzzzz_0_yyzzzzz_1, g_xzzzz_0_yzzzzz_1, g_xzzzz_0_yzzzzzz_1, g_xzzzz_0_zzzzzz_1, g_xzzzz_0_zzzzzzz_1, g_zzzz_0_xxxxxxz_0, g_zzzz_0_xxxxxxz_1, g_zzzz_0_xxxxxyz_0, g_zzzz_0_xxxxxyz_1, g_zzzz_0_xxxxxzz_0, g_zzzz_0_xxxxxzz_1, g_zzzz_0_xxxxyyz_0, g_zzzz_0_xxxxyyz_1, g_zzzz_0_xxxxyzz_0, g_zzzz_0_xxxxyzz_1, g_zzzz_0_xxxxzzz_0, g_zzzz_0_xxxxzzz_1, g_zzzz_0_xxxyyyz_0, g_zzzz_0_xxxyyyz_1, g_zzzz_0_xxxyyzz_0, g_zzzz_0_xxxyyzz_1, g_zzzz_0_xxxyzzz_0, g_zzzz_0_xxxyzzz_1, g_zzzz_0_xxxzzzz_0, g_zzzz_0_xxxzzzz_1, g_zzzz_0_xxyyyyz_0, g_zzzz_0_xxyyyyz_1, g_zzzz_0_xxyyyzz_0, g_zzzz_0_xxyyyzz_1, g_zzzz_0_xxyyzzz_0, g_zzzz_0_xxyyzzz_1, g_zzzz_0_xxyzzzz_0, g_zzzz_0_xxyzzzz_1, g_zzzz_0_xxzzzzz_0, g_zzzz_0_xxzzzzz_1, g_zzzz_0_xyyyyyz_0, g_zzzz_0_xyyyyyz_1, g_zzzz_0_xyyyyzz_0, g_zzzz_0_xyyyyzz_1, g_zzzz_0_xyyyzzz_0, g_zzzz_0_xyyyzzz_1, g_zzzz_0_xyyzzzz_0, g_zzzz_0_xyyzzzz_1, g_zzzz_0_xyzzzzz_0, g_zzzz_0_xyzzzzz_1, g_zzzz_0_xzzzzzz_0, g_zzzz_0_xzzzzzz_1, g_zzzz_0_yyyyyyy_0, g_zzzz_0_yyyyyyy_1, g_zzzz_0_yyyyyyz_0, g_zzzz_0_yyyyyyz_1, g_zzzz_0_yyyyyzz_0, g_zzzz_0_yyyyyzz_1, g_zzzz_0_yyyyzzz_0, g_zzzz_0_yyyyzzz_1, g_zzzz_0_yyyzzzz_0, g_zzzz_0_yyyzzzz_1, g_zzzz_0_yyzzzzz_0, g_zzzz_0_yyzzzzz_1, g_zzzz_0_yzzzzzz_0, g_zzzz_0_yzzzzzz_1, g_zzzz_0_zzzzzzz_0, g_zzzz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzz_0_xxxxxxx_0[i] = 3.0 * g_xxzz_0_xxxxxxx_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxxx_1[i] * fz_be_0 + g_xxzzz_0_xxxxxxx_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxxxy_0[i] = 3.0 * g_xxzz_0_xxxxxxy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxxy_1[i] * fz_be_0 + g_xxzzz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxxxz_0[i] = g_zzzz_0_xxxxxxz_0[i] * fbe_0 - g_zzzz_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxxyy_0[i] = 3.0 * g_xxzz_0_xxxxxyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxxyy_1[i] * fz_be_0 + g_xxzzz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxxyz_0[i] = g_zzzz_0_xxxxxyz_0[i] * fbe_0 - g_zzzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxxzz_0[i] = g_zzzz_0_xxxxxzz_0[i] * fbe_0 - g_zzzz_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxyyy_0[i] = 3.0 * g_xxzz_0_xxxxyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxxyyy_1[i] * fz_be_0 + g_xxzzz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxxyyz_0[i] = g_zzzz_0_xxxxyyz_0[i] * fbe_0 - g_zzzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxyzz_0[i] = g_zzzz_0_xxxxyzz_0[i] * fbe_0 - g_zzzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxxzzz_0[i] = g_zzzz_0_xxxxzzz_0[i] * fbe_0 - g_zzzz_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxyyyy_0[i] = 3.0 * g_xxzz_0_xxxyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxxyyyy_1[i] * fz_be_0 + g_xxzzz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxxyyyz_0[i] = g_zzzz_0_xxxyyyz_0[i] * fbe_0 - g_zzzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_xzzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxyyzz_0[i] = g_zzzz_0_xxxyyzz_0[i] * fbe_0 - g_zzzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxyzzz_0[i] = g_zzzz_0_xxxyzzz_0[i] * fbe_0 - g_zzzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxxzzzz_0[i] = g_zzzz_0_xxxzzzz_0[i] * fbe_0 - g_zzzz_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyyyyy_0[i] = 3.0 * g_xxzz_0_xxyyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xxyyyyy_1[i] * fz_be_0 + g_xxzzz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xxyyyyz_0[i] = g_zzzz_0_xxyyyyz_0[i] * fbe_0 - g_zzzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_xzzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyyyzz_0[i] = g_zzzz_0_xxyyyzz_0[i] * fbe_0 - g_zzzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_xzzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyyzzz_0[i] = g_zzzz_0_xxyyzzz_0[i] * fbe_0 - g_zzzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxyzzzz_0[i] = g_zzzz_0_xxyzzzz_0[i] * fbe_0 - g_zzzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xxzzzzz_0[i] = g_zzzz_0_xxzzzzz_0[i] * fbe_0 - g_zzzz_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyyyyy_0[i] = 3.0 * g_xxzz_0_xyyyyyy_0[i] * fbe_0 - 3.0 * g_xxzz_0_xyyyyyy_1[i] * fz_be_0 + g_xxzzz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxzzzz_0_xyyyyyz_0[i] = g_zzzz_0_xyyyyyz_0[i] * fbe_0 - g_zzzz_0_xyyyyyz_1[i] * fz_be_0 + g_xzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_xzzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyyyzz_0[i] = g_zzzz_0_xyyyyzz_0[i] * fbe_0 - g_zzzz_0_xyyyyzz_1[i] * fz_be_0 + g_xzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_xzzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyyzzz_0[i] = g_zzzz_0_xyyyzzz_0[i] * fbe_0 - g_zzzz_0_xyyyzzz_1[i] * fz_be_0 + g_xzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_xzzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyyzzzz_0[i] = g_zzzz_0_xyyzzzz_0[i] * fbe_0 - g_zzzz_0_xyyzzzz_1[i] * fz_be_0 + g_xzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xyzzzzz_0[i] = g_zzzz_0_xyzzzzz_0[i] * fbe_0 - g_zzzz_0_xyzzzzz_1[i] * fz_be_0 + g_xzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_xzzzzzz_0[i] = g_zzzz_0_xzzzzzz_0[i] * fbe_0 - g_zzzz_0_xzzzzzz_1[i] * fz_be_0 + g_xzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_xzzzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyyyy_0[i] = g_zzzz_0_yyyyyyy_0[i] * fbe_0 - g_zzzz_0_yyyyyyy_1[i] * fz_be_0 + g_xzzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyyyz_0[i] = g_zzzz_0_yyyyyyz_0[i] * fbe_0 - g_zzzz_0_yyyyyyz_1[i] * fz_be_0 + g_xzzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyyzz_0[i] = g_zzzz_0_yyyyyzz_0[i] * fbe_0 - g_zzzz_0_yyyyyzz_1[i] * fz_be_0 + g_xzzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyyzzz_0[i] = g_zzzz_0_yyyyzzz_0[i] * fbe_0 - g_zzzz_0_yyyyzzz_1[i] * fz_be_0 + g_xzzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyyzzzz_0[i] = g_zzzz_0_yyyzzzz_0[i] * fbe_0 - g_zzzz_0_yyyzzzz_1[i] * fz_be_0 + g_xzzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yyzzzzz_0[i] = g_zzzz_0_yyzzzzz_0[i] * fbe_0 - g_zzzz_0_yyzzzzz_1[i] * fz_be_0 + g_xzzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_yzzzzzz_0[i] = g_zzzz_0_yzzzzzz_0[i] * fbe_0 - g_zzzz_0_yzzzzzz_1[i] * fz_be_0 + g_xzzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxzzzz_0_zzzzzzz_0[i] = g_zzzz_0_zzzzzzz_0[i] * fbe_0 - g_zzzz_0_zzzzzzz_1[i] * fz_be_0 + g_xzzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 540-576 components of targeted buffer : ISK

    auto g_xyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 540);

    auto g_xyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 541);

    auto g_xyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 542);

    auto g_xyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 543);

    auto g_xyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 544);

    auto g_xyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 545);

    auto g_xyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 546);

    auto g_xyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 547);

    auto g_xyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 548);

    auto g_xyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 549);

    auto g_xyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 550);

    auto g_xyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 551);

    auto g_xyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 552);

    auto g_xyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 553);

    auto g_xyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 554);

    auto g_xyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 555);

    auto g_xyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 556);

    auto g_xyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 557);

    auto g_xyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 558);

    auto g_xyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 559);

    auto g_xyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 560);

    auto g_xyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 561);

    auto g_xyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 562);

    auto g_xyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 563);

    auto g_xyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 564);

    auto g_xyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 565);

    auto g_xyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 566);

    auto g_xyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 567);

    auto g_xyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 568);

    auto g_xyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 569);

    auto g_xyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 570);

    auto g_xyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 571);

    auto g_xyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 572);

    auto g_xyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 573);

    auto g_xyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 574);

    auto g_xyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 575);

    #pragma omp simd aligned(g_xyyyyy_0_xxxxxxx_0, g_xyyyyy_0_xxxxxxy_0, g_xyyyyy_0_xxxxxxz_0, g_xyyyyy_0_xxxxxyy_0, g_xyyyyy_0_xxxxxyz_0, g_xyyyyy_0_xxxxxzz_0, g_xyyyyy_0_xxxxyyy_0, g_xyyyyy_0_xxxxyyz_0, g_xyyyyy_0_xxxxyzz_0, g_xyyyyy_0_xxxxzzz_0, g_xyyyyy_0_xxxyyyy_0, g_xyyyyy_0_xxxyyyz_0, g_xyyyyy_0_xxxyyzz_0, g_xyyyyy_0_xxxyzzz_0, g_xyyyyy_0_xxxzzzz_0, g_xyyyyy_0_xxyyyyy_0, g_xyyyyy_0_xxyyyyz_0, g_xyyyyy_0_xxyyyzz_0, g_xyyyyy_0_xxyyzzz_0, g_xyyyyy_0_xxyzzzz_0, g_xyyyyy_0_xxzzzzz_0, g_xyyyyy_0_xyyyyyy_0, g_xyyyyy_0_xyyyyyz_0, g_xyyyyy_0_xyyyyzz_0, g_xyyyyy_0_xyyyzzz_0, g_xyyyyy_0_xyyzzzz_0, g_xyyyyy_0_xyzzzzz_0, g_xyyyyy_0_xzzzzzz_0, g_xyyyyy_0_yyyyyyy_0, g_xyyyyy_0_yyyyyyz_0, g_xyyyyy_0_yyyyyzz_0, g_xyyyyy_0_yyyyzzz_0, g_xyyyyy_0_yyyzzzz_0, g_xyyyyy_0_yyzzzzz_0, g_xyyyyy_0_yzzzzzz_0, g_xyyyyy_0_zzzzzzz_0, g_yyyyy_0_xxxxxx_1, g_yyyyy_0_xxxxxxx_1, g_yyyyy_0_xxxxxxy_1, g_yyyyy_0_xxxxxxz_1, g_yyyyy_0_xxxxxy_1, g_yyyyy_0_xxxxxyy_1, g_yyyyy_0_xxxxxyz_1, g_yyyyy_0_xxxxxz_1, g_yyyyy_0_xxxxxzz_1, g_yyyyy_0_xxxxyy_1, g_yyyyy_0_xxxxyyy_1, g_yyyyy_0_xxxxyyz_1, g_yyyyy_0_xxxxyz_1, g_yyyyy_0_xxxxyzz_1, g_yyyyy_0_xxxxzz_1, g_yyyyy_0_xxxxzzz_1, g_yyyyy_0_xxxyyy_1, g_yyyyy_0_xxxyyyy_1, g_yyyyy_0_xxxyyyz_1, g_yyyyy_0_xxxyyz_1, g_yyyyy_0_xxxyyzz_1, g_yyyyy_0_xxxyzz_1, g_yyyyy_0_xxxyzzz_1, g_yyyyy_0_xxxzzz_1, g_yyyyy_0_xxxzzzz_1, g_yyyyy_0_xxyyyy_1, g_yyyyy_0_xxyyyyy_1, g_yyyyy_0_xxyyyyz_1, g_yyyyy_0_xxyyyz_1, g_yyyyy_0_xxyyyzz_1, g_yyyyy_0_xxyyzz_1, g_yyyyy_0_xxyyzzz_1, g_yyyyy_0_xxyzzz_1, g_yyyyy_0_xxyzzzz_1, g_yyyyy_0_xxzzzz_1, g_yyyyy_0_xxzzzzz_1, g_yyyyy_0_xyyyyy_1, g_yyyyy_0_xyyyyyy_1, g_yyyyy_0_xyyyyyz_1, g_yyyyy_0_xyyyyz_1, g_yyyyy_0_xyyyyzz_1, g_yyyyy_0_xyyyzz_1, g_yyyyy_0_xyyyzzz_1, g_yyyyy_0_xyyzzz_1, g_yyyyy_0_xyyzzzz_1, g_yyyyy_0_xyzzzz_1, g_yyyyy_0_xyzzzzz_1, g_yyyyy_0_xzzzzz_1, g_yyyyy_0_xzzzzzz_1, g_yyyyy_0_yyyyyy_1, g_yyyyy_0_yyyyyyy_1, g_yyyyy_0_yyyyyyz_1, g_yyyyy_0_yyyyyz_1, g_yyyyy_0_yyyyyzz_1, g_yyyyy_0_yyyyzz_1, g_yyyyy_0_yyyyzzz_1, g_yyyyy_0_yyyzzz_1, g_yyyyy_0_yyyzzzz_1, g_yyyyy_0_yyzzzz_1, g_yyyyy_0_yyzzzzz_1, g_yyyyy_0_yzzzzz_1, g_yyyyy_0_yzzzzzz_1, g_yyyyy_0_zzzzzz_1, g_yyyyy_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyy_0_xxxxxxx_0[i] = 7.0 * g_yyyyy_0_xxxxxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxx_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxxy_0[i] = 6.0 * g_yyyyy_0_xxxxxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxxz_0[i] = 6.0 * g_yyyyy_0_xxxxxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxyy_0[i] = 5.0 * g_yyyyy_0_xxxxyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxyz_0[i] = 5.0 * g_yyyyy_0_xxxxyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxxzz_0[i] = 5.0 * g_yyyyy_0_xxxxzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxyyy_0[i] = 4.0 * g_yyyyy_0_xxxyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxyyz_0[i] = 4.0 * g_yyyyy_0_xxxyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxyzz_0[i] = 4.0 * g_yyyyy_0_xxxyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxxzzz_0[i] = 4.0 * g_yyyyy_0_xxxzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyyyy_0[i] = 3.0 * g_yyyyy_0_xxyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyyyz_0[i] = 3.0 * g_yyyyy_0_xxyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyyzz_0[i] = 3.0 * g_yyyyy_0_xxyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxyzzz_0[i] = 3.0 * g_yyyyy_0_xxyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxxzzzz_0[i] = 3.0 * g_yyyyy_0_xxzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyyyy_0[i] = 2.0 * g_yyyyy_0_xyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyyyz_0[i] = 2.0 * g_yyyyy_0_xyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyyzz_0[i] = 2.0 * g_yyyyy_0_xyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyyzzz_0[i] = 2.0 * g_yyyyy_0_xyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxyzzzz_0[i] = 2.0 * g_yyyyy_0_xyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xxzzzzz_0[i] = 2.0 * g_yyyyy_0_xzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyyyy_0[i] = g_yyyyy_0_yyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyyyz_0[i] = g_yyyyy_0_yyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyyzz_0[i] = g_yyyyy_0_yyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyyzzz_0[i] = g_yyyyy_0_yyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyyzzzz_0[i] = g_yyyyy_0_yyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xyzzzzz_0[i] = g_yyyyy_0_yzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_xzzzzzz_0[i] = g_yyyyy_0_zzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyyyy_0[i] = g_yyyyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyyyz_0[i] = g_yyyyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyyzz_0[i] = g_yyyyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyyzzz_0[i] = g_yyyyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyyzzzz_0[i] = g_yyyyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yyzzzzz_0[i] = g_yyyyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_yzzzzzz_0[i] = g_yyyyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyyyy_0_zzzzzzz_0[i] = g_yyyyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 576-612 components of targeted buffer : ISK

    auto g_xyyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 576);

    auto g_xyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 577);

    auto g_xyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 578);

    auto g_xyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 579);

    auto g_xyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 580);

    auto g_xyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 581);

    auto g_xyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 582);

    auto g_xyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 583);

    auto g_xyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 584);

    auto g_xyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 585);

    auto g_xyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 586);

    auto g_xyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 587);

    auto g_xyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 588);

    auto g_xyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 589);

    auto g_xyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 590);

    auto g_xyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 591);

    auto g_xyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 592);

    auto g_xyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 593);

    auto g_xyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 594);

    auto g_xyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 595);

    auto g_xyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 596);

    auto g_xyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 597);

    auto g_xyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 598);

    auto g_xyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 599);

    auto g_xyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 600);

    auto g_xyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 601);

    auto g_xyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 602);

    auto g_xyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 603);

    auto g_xyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 604);

    auto g_xyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 605);

    auto g_xyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 606);

    auto g_xyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 607);

    auto g_xyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 608);

    auto g_xyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 609);

    auto g_xyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 610);

    auto g_xyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 611);

    #pragma omp simd aligned(g_xyyyy_0_xxxxxxx_1, g_xyyyy_0_xxxxxxy_1, g_xyyyy_0_xxxxxyy_1, g_xyyyy_0_xxxxyyy_1, g_xyyyy_0_xxxyyyy_1, g_xyyyy_0_xxyyyyy_1, g_xyyyy_0_xyyyyyy_1, g_xyyyyz_0_xxxxxxx_0, g_xyyyyz_0_xxxxxxy_0, g_xyyyyz_0_xxxxxxz_0, g_xyyyyz_0_xxxxxyy_0, g_xyyyyz_0_xxxxxyz_0, g_xyyyyz_0_xxxxxzz_0, g_xyyyyz_0_xxxxyyy_0, g_xyyyyz_0_xxxxyyz_0, g_xyyyyz_0_xxxxyzz_0, g_xyyyyz_0_xxxxzzz_0, g_xyyyyz_0_xxxyyyy_0, g_xyyyyz_0_xxxyyyz_0, g_xyyyyz_0_xxxyyzz_0, g_xyyyyz_0_xxxyzzz_0, g_xyyyyz_0_xxxzzzz_0, g_xyyyyz_0_xxyyyyy_0, g_xyyyyz_0_xxyyyyz_0, g_xyyyyz_0_xxyyyzz_0, g_xyyyyz_0_xxyyzzz_0, g_xyyyyz_0_xxyzzzz_0, g_xyyyyz_0_xxzzzzz_0, g_xyyyyz_0_xyyyyyy_0, g_xyyyyz_0_xyyyyyz_0, g_xyyyyz_0_xyyyyzz_0, g_xyyyyz_0_xyyyzzz_0, g_xyyyyz_0_xyyzzzz_0, g_xyyyyz_0_xyzzzzz_0, g_xyyyyz_0_xzzzzzz_0, g_xyyyyz_0_yyyyyyy_0, g_xyyyyz_0_yyyyyyz_0, g_xyyyyz_0_yyyyyzz_0, g_xyyyyz_0_yyyyzzz_0, g_xyyyyz_0_yyyzzzz_0, g_xyyyyz_0_yyzzzzz_0, g_xyyyyz_0_yzzzzzz_0, g_xyyyyz_0_zzzzzzz_0, g_yyyyz_0_xxxxxxz_1, g_yyyyz_0_xxxxxyz_1, g_yyyyz_0_xxxxxz_1, g_yyyyz_0_xxxxxzz_1, g_yyyyz_0_xxxxyyz_1, g_yyyyz_0_xxxxyz_1, g_yyyyz_0_xxxxyzz_1, g_yyyyz_0_xxxxzz_1, g_yyyyz_0_xxxxzzz_1, g_yyyyz_0_xxxyyyz_1, g_yyyyz_0_xxxyyz_1, g_yyyyz_0_xxxyyzz_1, g_yyyyz_0_xxxyzz_1, g_yyyyz_0_xxxyzzz_1, g_yyyyz_0_xxxzzz_1, g_yyyyz_0_xxxzzzz_1, g_yyyyz_0_xxyyyyz_1, g_yyyyz_0_xxyyyz_1, g_yyyyz_0_xxyyyzz_1, g_yyyyz_0_xxyyzz_1, g_yyyyz_0_xxyyzzz_1, g_yyyyz_0_xxyzzz_1, g_yyyyz_0_xxyzzzz_1, g_yyyyz_0_xxzzzz_1, g_yyyyz_0_xxzzzzz_1, g_yyyyz_0_xyyyyyz_1, g_yyyyz_0_xyyyyz_1, g_yyyyz_0_xyyyyzz_1, g_yyyyz_0_xyyyzz_1, g_yyyyz_0_xyyyzzz_1, g_yyyyz_0_xyyzzz_1, g_yyyyz_0_xyyzzzz_1, g_yyyyz_0_xyzzzz_1, g_yyyyz_0_xyzzzzz_1, g_yyyyz_0_xzzzzz_1, g_yyyyz_0_xzzzzzz_1, g_yyyyz_0_yyyyyyy_1, g_yyyyz_0_yyyyyyz_1, g_yyyyz_0_yyyyyz_1, g_yyyyz_0_yyyyyzz_1, g_yyyyz_0_yyyyzz_1, g_yyyyz_0_yyyyzzz_1, g_yyyyz_0_yyyzzz_1, g_yyyyz_0_yyyzzzz_1, g_yyyyz_0_yyzzzz_1, g_yyyyz_0_yyzzzzz_1, g_yyyyz_0_yzzzzz_1, g_yyyyz_0_yzzzzzz_1, g_yyyyz_0_zzzzzz_1, g_yyyyz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyz_0_xxxxxxx_0[i] = g_xyyyy_0_xxxxxxx_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxxxy_0[i] = g_xyyyy_0_xxxxxxy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxxxz_0[i] = 6.0 * g_yyyyz_0_xxxxxz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxxyy_0[i] = g_xyyyy_0_xxxxxyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxxyz_0[i] = 5.0 * g_yyyyz_0_xxxxyz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxxzz_0[i] = 5.0 * g_yyyyz_0_xxxxzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxyyy_0[i] = g_xyyyy_0_xxxxyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxxyyz_0[i] = 4.0 * g_yyyyz_0_xxxyyz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxyzz_0[i] = 4.0 * g_yyyyz_0_xxxyzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxxzzz_0[i] = 4.0 * g_yyyyz_0_xxxzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxyyyy_0[i] = g_xyyyy_0_xxxyyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxxyyyz_0[i] = 3.0 * g_yyyyz_0_xxyyyz_1[i] * fi_acd_0 + g_yyyyz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxyyzz_0[i] = 3.0 * g_yyyyz_0_xxyyzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxyzzz_0[i] = 3.0 * g_yyyyz_0_xxyzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxxzzzz_0[i] = 3.0 * g_yyyyz_0_xxzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyyyyy_0[i] = g_xyyyy_0_xxyyyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xxyyyyz_0[i] = 2.0 * g_yyyyz_0_xyyyyz_1[i] * fi_acd_0 + g_yyyyz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyyyzz_0[i] = 2.0 * g_yyyyz_0_xyyyzz_1[i] * fi_acd_0 + g_yyyyz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyyzzz_0[i] = 2.0 * g_yyyyz_0_xyyzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxyzzzz_0[i] = 2.0 * g_yyyyz_0_xyzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xxzzzzz_0[i] = 2.0 * g_yyyyz_0_xzzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyyyyy_0[i] = g_xyyyy_0_xyyyyyy_1[i] * wa_z[i];

        g_xyyyyz_0_xyyyyyz_0[i] = g_yyyyz_0_yyyyyz_1[i] * fi_acd_0 + g_yyyyz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyyyzz_0[i] = g_yyyyz_0_yyyyzz_1[i] * fi_acd_0 + g_yyyyz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyyzzz_0[i] = g_yyyyz_0_yyyzzz_1[i] * fi_acd_0 + g_yyyyz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyyzzzz_0[i] = g_yyyyz_0_yyzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xyzzzzz_0[i] = g_yyyyz_0_yzzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_xzzzzzz_0[i] = g_yyyyz_0_zzzzzz_1[i] * fi_acd_0 + g_yyyyz_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyyyy_0[i] = g_yyyyz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyyyz_0[i] = g_yyyyz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyyzz_0[i] = g_yyyyz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyyzzz_0[i] = g_yyyyz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyyzzzz_0[i] = g_yyyyz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yyzzzzz_0[i] = g_yyyyz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_yzzzzzz_0[i] = g_yyyyz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyyyz_0_zzzzzzz_0[i] = g_yyyyz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 612-648 components of targeted buffer : ISK

    auto g_xyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 612);

    auto g_xyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 613);

    auto g_xyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 614);

    auto g_xyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 615);

    auto g_xyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 616);

    auto g_xyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 617);

    auto g_xyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 618);

    auto g_xyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 619);

    auto g_xyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 620);

    auto g_xyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 621);

    auto g_xyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 622);

    auto g_xyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 623);

    auto g_xyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 624);

    auto g_xyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 625);

    auto g_xyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 626);

    auto g_xyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 627);

    auto g_xyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 628);

    auto g_xyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 629);

    auto g_xyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 630);

    auto g_xyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 631);

    auto g_xyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 632);

    auto g_xyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 633);

    auto g_xyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 634);

    auto g_xyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 635);

    auto g_xyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 636);

    auto g_xyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 637);

    auto g_xyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 638);

    auto g_xyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 639);

    auto g_xyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 640);

    auto g_xyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 641);

    auto g_xyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 642);

    auto g_xyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 643);

    auto g_xyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 644);

    auto g_xyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 645);

    auto g_xyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 646);

    auto g_xyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 647);

    #pragma omp simd aligned(g_xyyyzz_0_xxxxxxx_0, g_xyyyzz_0_xxxxxxy_0, g_xyyyzz_0_xxxxxxz_0, g_xyyyzz_0_xxxxxyy_0, g_xyyyzz_0_xxxxxyz_0, g_xyyyzz_0_xxxxxzz_0, g_xyyyzz_0_xxxxyyy_0, g_xyyyzz_0_xxxxyyz_0, g_xyyyzz_0_xxxxyzz_0, g_xyyyzz_0_xxxxzzz_0, g_xyyyzz_0_xxxyyyy_0, g_xyyyzz_0_xxxyyyz_0, g_xyyyzz_0_xxxyyzz_0, g_xyyyzz_0_xxxyzzz_0, g_xyyyzz_0_xxxzzzz_0, g_xyyyzz_0_xxyyyyy_0, g_xyyyzz_0_xxyyyyz_0, g_xyyyzz_0_xxyyyzz_0, g_xyyyzz_0_xxyyzzz_0, g_xyyyzz_0_xxyzzzz_0, g_xyyyzz_0_xxzzzzz_0, g_xyyyzz_0_xyyyyyy_0, g_xyyyzz_0_xyyyyyz_0, g_xyyyzz_0_xyyyyzz_0, g_xyyyzz_0_xyyyzzz_0, g_xyyyzz_0_xyyzzzz_0, g_xyyyzz_0_xyzzzzz_0, g_xyyyzz_0_xzzzzzz_0, g_xyyyzz_0_yyyyyyy_0, g_xyyyzz_0_yyyyyyz_0, g_xyyyzz_0_yyyyyzz_0, g_xyyyzz_0_yyyyzzz_0, g_xyyyzz_0_yyyzzzz_0, g_xyyyzz_0_yyzzzzz_0, g_xyyyzz_0_yzzzzzz_0, g_xyyyzz_0_zzzzzzz_0, g_yyyzz_0_xxxxxx_1, g_yyyzz_0_xxxxxxx_1, g_yyyzz_0_xxxxxxy_1, g_yyyzz_0_xxxxxxz_1, g_yyyzz_0_xxxxxy_1, g_yyyzz_0_xxxxxyy_1, g_yyyzz_0_xxxxxyz_1, g_yyyzz_0_xxxxxz_1, g_yyyzz_0_xxxxxzz_1, g_yyyzz_0_xxxxyy_1, g_yyyzz_0_xxxxyyy_1, g_yyyzz_0_xxxxyyz_1, g_yyyzz_0_xxxxyz_1, g_yyyzz_0_xxxxyzz_1, g_yyyzz_0_xxxxzz_1, g_yyyzz_0_xxxxzzz_1, g_yyyzz_0_xxxyyy_1, g_yyyzz_0_xxxyyyy_1, g_yyyzz_0_xxxyyyz_1, g_yyyzz_0_xxxyyz_1, g_yyyzz_0_xxxyyzz_1, g_yyyzz_0_xxxyzz_1, g_yyyzz_0_xxxyzzz_1, g_yyyzz_0_xxxzzz_1, g_yyyzz_0_xxxzzzz_1, g_yyyzz_0_xxyyyy_1, g_yyyzz_0_xxyyyyy_1, g_yyyzz_0_xxyyyyz_1, g_yyyzz_0_xxyyyz_1, g_yyyzz_0_xxyyyzz_1, g_yyyzz_0_xxyyzz_1, g_yyyzz_0_xxyyzzz_1, g_yyyzz_0_xxyzzz_1, g_yyyzz_0_xxyzzzz_1, g_yyyzz_0_xxzzzz_1, g_yyyzz_0_xxzzzzz_1, g_yyyzz_0_xyyyyy_1, g_yyyzz_0_xyyyyyy_1, g_yyyzz_0_xyyyyyz_1, g_yyyzz_0_xyyyyz_1, g_yyyzz_0_xyyyyzz_1, g_yyyzz_0_xyyyzz_1, g_yyyzz_0_xyyyzzz_1, g_yyyzz_0_xyyzzz_1, g_yyyzz_0_xyyzzzz_1, g_yyyzz_0_xyzzzz_1, g_yyyzz_0_xyzzzzz_1, g_yyyzz_0_xzzzzz_1, g_yyyzz_0_xzzzzzz_1, g_yyyzz_0_yyyyyy_1, g_yyyzz_0_yyyyyyy_1, g_yyyzz_0_yyyyyyz_1, g_yyyzz_0_yyyyyz_1, g_yyyzz_0_yyyyyzz_1, g_yyyzz_0_yyyyzz_1, g_yyyzz_0_yyyyzzz_1, g_yyyzz_0_yyyzzz_1, g_yyyzz_0_yyyzzzz_1, g_yyyzz_0_yyzzzz_1, g_yyyzz_0_yyzzzzz_1, g_yyyzz_0_yzzzzz_1, g_yyyzz_0_yzzzzzz_1, g_yyyzz_0_zzzzzz_1, g_yyyzz_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzz_0_xxxxxxx_0[i] = 7.0 * g_yyyzz_0_xxxxxx_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxxx_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxxy_0[i] = 6.0 * g_yyyzz_0_xxxxxy_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxxz_0[i] = 6.0 * g_yyyzz_0_xxxxxz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxyy_0[i] = 5.0 * g_yyyzz_0_xxxxyy_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxyz_0[i] = 5.0 * g_yyyzz_0_xxxxyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxxzz_0[i] = 5.0 * g_yyyzz_0_xxxxzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxyyy_0[i] = 4.0 * g_yyyzz_0_xxxyyy_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxyyz_0[i] = 4.0 * g_yyyzz_0_xxxyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxyzz_0[i] = 4.0 * g_yyyzz_0_xxxyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxxzzz_0[i] = 4.0 * g_yyyzz_0_xxxzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyyyy_0[i] = 3.0 * g_yyyzz_0_xxyyyy_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyyyz_0[i] = 3.0 * g_yyyzz_0_xxyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyyzz_0[i] = 3.0 * g_yyyzz_0_xxyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxyzzz_0[i] = 3.0 * g_yyyzz_0_xxyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxxzzzz_0[i] = 3.0 * g_yyyzz_0_xxzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyyyy_0[i] = 2.0 * g_yyyzz_0_xyyyyy_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyyyz_0[i] = 2.0 * g_yyyzz_0_xyyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyyzz_0[i] = 2.0 * g_yyyzz_0_xyyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyyzzz_0[i] = 2.0 * g_yyyzz_0_xyyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxyzzzz_0[i] = 2.0 * g_yyyzz_0_xyzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xxzzzzz_0[i] = 2.0 * g_yyyzz_0_xzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyyyy_0[i] = g_yyyzz_0_yyyyyy_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyyyz_0[i] = g_yyyzz_0_yyyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyyzz_0[i] = g_yyyzz_0_yyyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyyzzz_0[i] = g_yyyzz_0_yyyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyyzzzz_0[i] = g_yyyzz_0_yyzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xyzzzzz_0[i] = g_yyyzz_0_yzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_xzzzzzz_0[i] = g_yyyzz_0_zzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyyyy_0[i] = g_yyyzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyyyz_0[i] = g_yyyzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyyzz_0[i] = g_yyyzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyyzzz_0[i] = g_yyyzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyyzzzz_0[i] = g_yyyzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yyzzzzz_0[i] = g_yyyzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_yzzzzzz_0[i] = g_yyyzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyyzz_0_zzzzzzz_0[i] = g_yyyzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 648-684 components of targeted buffer : ISK

    auto g_xyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 648);

    auto g_xyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 649);

    auto g_xyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 650);

    auto g_xyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 651);

    auto g_xyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 652);

    auto g_xyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 653);

    auto g_xyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 654);

    auto g_xyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 655);

    auto g_xyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 656);

    auto g_xyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 657);

    auto g_xyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 658);

    auto g_xyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 659);

    auto g_xyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 660);

    auto g_xyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 661);

    auto g_xyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 662);

    auto g_xyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 663);

    auto g_xyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 664);

    auto g_xyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 665);

    auto g_xyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 666);

    auto g_xyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 667);

    auto g_xyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 668);

    auto g_xyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 669);

    auto g_xyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 670);

    auto g_xyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 671);

    auto g_xyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 672);

    auto g_xyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 673);

    auto g_xyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 674);

    auto g_xyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 675);

    auto g_xyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 676);

    auto g_xyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 677);

    auto g_xyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 678);

    auto g_xyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 679);

    auto g_xyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 680);

    auto g_xyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 681);

    auto g_xyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 682);

    auto g_xyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 683);

    #pragma omp simd aligned(g_xyyzzz_0_xxxxxxx_0, g_xyyzzz_0_xxxxxxy_0, g_xyyzzz_0_xxxxxxz_0, g_xyyzzz_0_xxxxxyy_0, g_xyyzzz_0_xxxxxyz_0, g_xyyzzz_0_xxxxxzz_0, g_xyyzzz_0_xxxxyyy_0, g_xyyzzz_0_xxxxyyz_0, g_xyyzzz_0_xxxxyzz_0, g_xyyzzz_0_xxxxzzz_0, g_xyyzzz_0_xxxyyyy_0, g_xyyzzz_0_xxxyyyz_0, g_xyyzzz_0_xxxyyzz_0, g_xyyzzz_0_xxxyzzz_0, g_xyyzzz_0_xxxzzzz_0, g_xyyzzz_0_xxyyyyy_0, g_xyyzzz_0_xxyyyyz_0, g_xyyzzz_0_xxyyyzz_0, g_xyyzzz_0_xxyyzzz_0, g_xyyzzz_0_xxyzzzz_0, g_xyyzzz_0_xxzzzzz_0, g_xyyzzz_0_xyyyyyy_0, g_xyyzzz_0_xyyyyyz_0, g_xyyzzz_0_xyyyyzz_0, g_xyyzzz_0_xyyyzzz_0, g_xyyzzz_0_xyyzzzz_0, g_xyyzzz_0_xyzzzzz_0, g_xyyzzz_0_xzzzzzz_0, g_xyyzzz_0_yyyyyyy_0, g_xyyzzz_0_yyyyyyz_0, g_xyyzzz_0_yyyyyzz_0, g_xyyzzz_0_yyyyzzz_0, g_xyyzzz_0_yyyzzzz_0, g_xyyzzz_0_yyzzzzz_0, g_xyyzzz_0_yzzzzzz_0, g_xyyzzz_0_zzzzzzz_0, g_yyzzz_0_xxxxxx_1, g_yyzzz_0_xxxxxxx_1, g_yyzzz_0_xxxxxxy_1, g_yyzzz_0_xxxxxxz_1, g_yyzzz_0_xxxxxy_1, g_yyzzz_0_xxxxxyy_1, g_yyzzz_0_xxxxxyz_1, g_yyzzz_0_xxxxxz_1, g_yyzzz_0_xxxxxzz_1, g_yyzzz_0_xxxxyy_1, g_yyzzz_0_xxxxyyy_1, g_yyzzz_0_xxxxyyz_1, g_yyzzz_0_xxxxyz_1, g_yyzzz_0_xxxxyzz_1, g_yyzzz_0_xxxxzz_1, g_yyzzz_0_xxxxzzz_1, g_yyzzz_0_xxxyyy_1, g_yyzzz_0_xxxyyyy_1, g_yyzzz_0_xxxyyyz_1, g_yyzzz_0_xxxyyz_1, g_yyzzz_0_xxxyyzz_1, g_yyzzz_0_xxxyzz_1, g_yyzzz_0_xxxyzzz_1, g_yyzzz_0_xxxzzz_1, g_yyzzz_0_xxxzzzz_1, g_yyzzz_0_xxyyyy_1, g_yyzzz_0_xxyyyyy_1, g_yyzzz_0_xxyyyyz_1, g_yyzzz_0_xxyyyz_1, g_yyzzz_0_xxyyyzz_1, g_yyzzz_0_xxyyzz_1, g_yyzzz_0_xxyyzzz_1, g_yyzzz_0_xxyzzz_1, g_yyzzz_0_xxyzzzz_1, g_yyzzz_0_xxzzzz_1, g_yyzzz_0_xxzzzzz_1, g_yyzzz_0_xyyyyy_1, g_yyzzz_0_xyyyyyy_1, g_yyzzz_0_xyyyyyz_1, g_yyzzz_0_xyyyyz_1, g_yyzzz_0_xyyyyzz_1, g_yyzzz_0_xyyyzz_1, g_yyzzz_0_xyyyzzz_1, g_yyzzz_0_xyyzzz_1, g_yyzzz_0_xyyzzzz_1, g_yyzzz_0_xyzzzz_1, g_yyzzz_0_xyzzzzz_1, g_yyzzz_0_xzzzzz_1, g_yyzzz_0_xzzzzzz_1, g_yyzzz_0_yyyyyy_1, g_yyzzz_0_yyyyyyy_1, g_yyzzz_0_yyyyyyz_1, g_yyzzz_0_yyyyyz_1, g_yyzzz_0_yyyyyzz_1, g_yyzzz_0_yyyyzz_1, g_yyzzz_0_yyyyzzz_1, g_yyzzz_0_yyyzzz_1, g_yyzzz_0_yyyzzzz_1, g_yyzzz_0_yyzzzz_1, g_yyzzz_0_yyzzzzz_1, g_yyzzz_0_yzzzzz_1, g_yyzzz_0_yzzzzzz_1, g_yyzzz_0_zzzzzz_1, g_yyzzz_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzz_0_xxxxxxx_0[i] = 7.0 * g_yyzzz_0_xxxxxx_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxxx_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxxy_0[i] = 6.0 * g_yyzzz_0_xxxxxy_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxxz_0[i] = 6.0 * g_yyzzz_0_xxxxxz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxyy_0[i] = 5.0 * g_yyzzz_0_xxxxyy_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxyz_0[i] = 5.0 * g_yyzzz_0_xxxxyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxxzz_0[i] = 5.0 * g_yyzzz_0_xxxxzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxyyy_0[i] = 4.0 * g_yyzzz_0_xxxyyy_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxyyz_0[i] = 4.0 * g_yyzzz_0_xxxyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxyzz_0[i] = 4.0 * g_yyzzz_0_xxxyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxxzzz_0[i] = 4.0 * g_yyzzz_0_xxxzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyyyy_0[i] = 3.0 * g_yyzzz_0_xxyyyy_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyyyz_0[i] = 3.0 * g_yyzzz_0_xxyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyyzz_0[i] = 3.0 * g_yyzzz_0_xxyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxyzzz_0[i] = 3.0 * g_yyzzz_0_xxyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxxzzzz_0[i] = 3.0 * g_yyzzz_0_xxzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyyyy_0[i] = 2.0 * g_yyzzz_0_xyyyyy_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyyyz_0[i] = 2.0 * g_yyzzz_0_xyyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyyzz_0[i] = 2.0 * g_yyzzz_0_xyyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyyzzz_0[i] = 2.0 * g_yyzzz_0_xyyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxyzzzz_0[i] = 2.0 * g_yyzzz_0_xyzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xxzzzzz_0[i] = 2.0 * g_yyzzz_0_xzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyyyy_0[i] = g_yyzzz_0_yyyyyy_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyyyz_0[i] = g_yyzzz_0_yyyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyyzz_0[i] = g_yyzzz_0_yyyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyyzzz_0[i] = g_yyzzz_0_yyyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyyzzzz_0[i] = g_yyzzz_0_yyzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xyzzzzz_0[i] = g_yyzzz_0_yzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_xzzzzzz_0[i] = g_yyzzz_0_zzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyyyy_0[i] = g_yyzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyyyz_0[i] = g_yyzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyyzz_0[i] = g_yyzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyyzzz_0[i] = g_yyzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyyzzzz_0[i] = g_yyzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yyzzzzz_0[i] = g_yyzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_yzzzzzz_0[i] = g_yyzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyzzz_0_zzzzzzz_0[i] = g_yyzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 684-720 components of targeted buffer : ISK

    auto g_xyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 684);

    auto g_xyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 685);

    auto g_xyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 686);

    auto g_xyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 687);

    auto g_xyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 688);

    auto g_xyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 689);

    auto g_xyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 690);

    auto g_xyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 691);

    auto g_xyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 692);

    auto g_xyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 693);

    auto g_xyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 694);

    auto g_xyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 695);

    auto g_xyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 696);

    auto g_xyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 697);

    auto g_xyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 698);

    auto g_xyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 699);

    auto g_xyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 700);

    auto g_xyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 701);

    auto g_xyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 702);

    auto g_xyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 703);

    auto g_xyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 704);

    auto g_xyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 705);

    auto g_xyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 706);

    auto g_xyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 707);

    auto g_xyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 708);

    auto g_xyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 709);

    auto g_xyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 710);

    auto g_xyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 711);

    auto g_xyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 712);

    auto g_xyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 713);

    auto g_xyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 714);

    auto g_xyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 715);

    auto g_xyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 716);

    auto g_xyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 717);

    auto g_xyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 718);

    auto g_xyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 719);

    #pragma omp simd aligned(g_xyzzzz_0_xxxxxxx_0, g_xyzzzz_0_xxxxxxy_0, g_xyzzzz_0_xxxxxxz_0, g_xyzzzz_0_xxxxxyy_0, g_xyzzzz_0_xxxxxyz_0, g_xyzzzz_0_xxxxxzz_0, g_xyzzzz_0_xxxxyyy_0, g_xyzzzz_0_xxxxyyz_0, g_xyzzzz_0_xxxxyzz_0, g_xyzzzz_0_xxxxzzz_0, g_xyzzzz_0_xxxyyyy_0, g_xyzzzz_0_xxxyyyz_0, g_xyzzzz_0_xxxyyzz_0, g_xyzzzz_0_xxxyzzz_0, g_xyzzzz_0_xxxzzzz_0, g_xyzzzz_0_xxyyyyy_0, g_xyzzzz_0_xxyyyyz_0, g_xyzzzz_0_xxyyyzz_0, g_xyzzzz_0_xxyyzzz_0, g_xyzzzz_0_xxyzzzz_0, g_xyzzzz_0_xxzzzzz_0, g_xyzzzz_0_xyyyyyy_0, g_xyzzzz_0_xyyyyyz_0, g_xyzzzz_0_xyyyyzz_0, g_xyzzzz_0_xyyyzzz_0, g_xyzzzz_0_xyyzzzz_0, g_xyzzzz_0_xyzzzzz_0, g_xyzzzz_0_xzzzzzz_0, g_xyzzzz_0_yyyyyyy_0, g_xyzzzz_0_yyyyyyz_0, g_xyzzzz_0_yyyyyzz_0, g_xyzzzz_0_yyyyzzz_0, g_xyzzzz_0_yyyzzzz_0, g_xyzzzz_0_yyzzzzz_0, g_xyzzzz_0_yzzzzzz_0, g_xyzzzz_0_zzzzzzz_0, g_xzzzz_0_xxxxxxx_1, g_xzzzz_0_xxxxxxz_1, g_xzzzz_0_xxxxxzz_1, g_xzzzz_0_xxxxzzz_1, g_xzzzz_0_xxxzzzz_1, g_xzzzz_0_xxzzzzz_1, g_xzzzz_0_xzzzzzz_1, g_yzzzz_0_xxxxxxy_1, g_yzzzz_0_xxxxxy_1, g_yzzzz_0_xxxxxyy_1, g_yzzzz_0_xxxxxyz_1, g_yzzzz_0_xxxxyy_1, g_yzzzz_0_xxxxyyy_1, g_yzzzz_0_xxxxyyz_1, g_yzzzz_0_xxxxyz_1, g_yzzzz_0_xxxxyzz_1, g_yzzzz_0_xxxyyy_1, g_yzzzz_0_xxxyyyy_1, g_yzzzz_0_xxxyyyz_1, g_yzzzz_0_xxxyyz_1, g_yzzzz_0_xxxyyzz_1, g_yzzzz_0_xxxyzz_1, g_yzzzz_0_xxxyzzz_1, g_yzzzz_0_xxyyyy_1, g_yzzzz_0_xxyyyyy_1, g_yzzzz_0_xxyyyyz_1, g_yzzzz_0_xxyyyz_1, g_yzzzz_0_xxyyyzz_1, g_yzzzz_0_xxyyzz_1, g_yzzzz_0_xxyyzzz_1, g_yzzzz_0_xxyzzz_1, g_yzzzz_0_xxyzzzz_1, g_yzzzz_0_xyyyyy_1, g_yzzzz_0_xyyyyyy_1, g_yzzzz_0_xyyyyyz_1, g_yzzzz_0_xyyyyz_1, g_yzzzz_0_xyyyyzz_1, g_yzzzz_0_xyyyzz_1, g_yzzzz_0_xyyyzzz_1, g_yzzzz_0_xyyzzz_1, g_yzzzz_0_xyyzzzz_1, g_yzzzz_0_xyzzzz_1, g_yzzzz_0_xyzzzzz_1, g_yzzzz_0_yyyyyy_1, g_yzzzz_0_yyyyyyy_1, g_yzzzz_0_yyyyyyz_1, g_yzzzz_0_yyyyyz_1, g_yzzzz_0_yyyyyzz_1, g_yzzzz_0_yyyyzz_1, g_yzzzz_0_yyyyzzz_1, g_yzzzz_0_yyyzzz_1, g_yzzzz_0_yyyzzzz_1, g_yzzzz_0_yyzzzz_1, g_yzzzz_0_yyzzzzz_1, g_yzzzz_0_yzzzzz_1, g_yzzzz_0_yzzzzzz_1, g_yzzzz_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzz_0_xxxxxxx_0[i] = g_xzzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xyzzzz_0_xxxxxxy_0[i] = 6.0 * g_yzzzz_0_xxxxxy_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxxxz_0[i] = g_xzzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xyzzzz_0_xxxxxyy_0[i] = 5.0 * g_yzzzz_0_xxxxyy_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxxyz_0[i] = 5.0 * g_yzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxxzz_0[i] = g_xzzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xyzzzz_0_xxxxyyy_0[i] = 4.0 * g_yzzzz_0_xxxyyy_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxyyz_0[i] = 4.0 * g_yzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxyzz_0[i] = 4.0 * g_yzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxxzzz_0[i] = g_xzzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xyzzzz_0_xxxyyyy_0[i] = 3.0 * g_yzzzz_0_xxyyyy_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxxyyyz_0[i] = 3.0 * g_yzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxyyzz_0[i] = 3.0 * g_yzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxyzzz_0[i] = 3.0 * g_yzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxxzzzz_0[i] = g_xzzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xyzzzz_0_xxyyyyy_0[i] = 2.0 * g_yzzzz_0_xyyyyy_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xxyyyyz_0[i] = 2.0 * g_yzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xxyyyzz_0[i] = 2.0 * g_yzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxyyzzz_0[i] = 2.0 * g_yzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxyzzzz_0[i] = 2.0 * g_yzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xxzzzzz_0[i] = g_xzzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xyzzzz_0_xyyyyyy_0[i] = g_yzzzz_0_yyyyyy_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_xyyyyyz_0[i] = g_yzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_xyyyyzz_0[i] = g_yzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_xyyyzzz_0[i] = g_yzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xyyzzzz_0[i] = g_yzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xyzzzzz_0[i] = g_yzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_xzzzzzz_0[i] = g_xzzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xyzzzz_0_yyyyyyy_0[i] = g_yzzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyzzzz_0_yyyyyyz_0[i] = g_yzzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyzzzz_0_yyyyyzz_0[i] = g_yzzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyzzzz_0_yyyyzzz_0[i] = g_yzzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyzzzz_0_yyyzzzz_0[i] = g_yzzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_yyzzzzz_0[i] = g_yzzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_yzzzzzz_0[i] = g_yzzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyzzzz_0_zzzzzzz_0[i] = g_yzzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 720-756 components of targeted buffer : ISK

    auto g_xzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 720);

    auto g_xzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 721);

    auto g_xzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 722);

    auto g_xzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 723);

    auto g_xzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 724);

    auto g_xzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 725);

    auto g_xzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 726);

    auto g_xzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 727);

    auto g_xzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 728);

    auto g_xzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 729);

    auto g_xzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 730);

    auto g_xzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 731);

    auto g_xzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 732);

    auto g_xzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 733);

    auto g_xzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 734);

    auto g_xzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 735);

    auto g_xzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 736);

    auto g_xzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 737);

    auto g_xzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 738);

    auto g_xzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 739);

    auto g_xzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 740);

    auto g_xzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 741);

    auto g_xzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 742);

    auto g_xzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 743);

    auto g_xzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 744);

    auto g_xzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 745);

    auto g_xzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 746);

    auto g_xzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 747);

    auto g_xzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 748);

    auto g_xzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 749);

    auto g_xzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 750);

    auto g_xzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 751);

    auto g_xzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 752);

    auto g_xzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 753);

    auto g_xzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 754);

    auto g_xzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 755);

    #pragma omp simd aligned(g_xzzzzz_0_xxxxxxx_0, g_xzzzzz_0_xxxxxxy_0, g_xzzzzz_0_xxxxxxz_0, g_xzzzzz_0_xxxxxyy_0, g_xzzzzz_0_xxxxxyz_0, g_xzzzzz_0_xxxxxzz_0, g_xzzzzz_0_xxxxyyy_0, g_xzzzzz_0_xxxxyyz_0, g_xzzzzz_0_xxxxyzz_0, g_xzzzzz_0_xxxxzzz_0, g_xzzzzz_0_xxxyyyy_0, g_xzzzzz_0_xxxyyyz_0, g_xzzzzz_0_xxxyyzz_0, g_xzzzzz_0_xxxyzzz_0, g_xzzzzz_0_xxxzzzz_0, g_xzzzzz_0_xxyyyyy_0, g_xzzzzz_0_xxyyyyz_0, g_xzzzzz_0_xxyyyzz_0, g_xzzzzz_0_xxyyzzz_0, g_xzzzzz_0_xxyzzzz_0, g_xzzzzz_0_xxzzzzz_0, g_xzzzzz_0_xyyyyyy_0, g_xzzzzz_0_xyyyyyz_0, g_xzzzzz_0_xyyyyzz_0, g_xzzzzz_0_xyyyzzz_0, g_xzzzzz_0_xyyzzzz_0, g_xzzzzz_0_xyzzzzz_0, g_xzzzzz_0_xzzzzzz_0, g_xzzzzz_0_yyyyyyy_0, g_xzzzzz_0_yyyyyyz_0, g_xzzzzz_0_yyyyyzz_0, g_xzzzzz_0_yyyyzzz_0, g_xzzzzz_0_yyyzzzz_0, g_xzzzzz_0_yyzzzzz_0, g_xzzzzz_0_yzzzzzz_0, g_xzzzzz_0_zzzzzzz_0, g_zzzzz_0_xxxxxx_1, g_zzzzz_0_xxxxxxx_1, g_zzzzz_0_xxxxxxy_1, g_zzzzz_0_xxxxxxz_1, g_zzzzz_0_xxxxxy_1, g_zzzzz_0_xxxxxyy_1, g_zzzzz_0_xxxxxyz_1, g_zzzzz_0_xxxxxz_1, g_zzzzz_0_xxxxxzz_1, g_zzzzz_0_xxxxyy_1, g_zzzzz_0_xxxxyyy_1, g_zzzzz_0_xxxxyyz_1, g_zzzzz_0_xxxxyz_1, g_zzzzz_0_xxxxyzz_1, g_zzzzz_0_xxxxzz_1, g_zzzzz_0_xxxxzzz_1, g_zzzzz_0_xxxyyy_1, g_zzzzz_0_xxxyyyy_1, g_zzzzz_0_xxxyyyz_1, g_zzzzz_0_xxxyyz_1, g_zzzzz_0_xxxyyzz_1, g_zzzzz_0_xxxyzz_1, g_zzzzz_0_xxxyzzz_1, g_zzzzz_0_xxxzzz_1, g_zzzzz_0_xxxzzzz_1, g_zzzzz_0_xxyyyy_1, g_zzzzz_0_xxyyyyy_1, g_zzzzz_0_xxyyyyz_1, g_zzzzz_0_xxyyyz_1, g_zzzzz_0_xxyyyzz_1, g_zzzzz_0_xxyyzz_1, g_zzzzz_0_xxyyzzz_1, g_zzzzz_0_xxyzzz_1, g_zzzzz_0_xxyzzzz_1, g_zzzzz_0_xxzzzz_1, g_zzzzz_0_xxzzzzz_1, g_zzzzz_0_xyyyyy_1, g_zzzzz_0_xyyyyyy_1, g_zzzzz_0_xyyyyyz_1, g_zzzzz_0_xyyyyz_1, g_zzzzz_0_xyyyyzz_1, g_zzzzz_0_xyyyzz_1, g_zzzzz_0_xyyyzzz_1, g_zzzzz_0_xyyzzz_1, g_zzzzz_0_xyyzzzz_1, g_zzzzz_0_xyzzzz_1, g_zzzzz_0_xyzzzzz_1, g_zzzzz_0_xzzzzz_1, g_zzzzz_0_xzzzzzz_1, g_zzzzz_0_yyyyyy_1, g_zzzzz_0_yyyyyyy_1, g_zzzzz_0_yyyyyyz_1, g_zzzzz_0_yyyyyz_1, g_zzzzz_0_yyyyyzz_1, g_zzzzz_0_yyyyzz_1, g_zzzzz_0_yyyyzzz_1, g_zzzzz_0_yyyzzz_1, g_zzzzz_0_yyyzzzz_1, g_zzzzz_0_yyzzzz_1, g_zzzzz_0_yyzzzzz_1, g_zzzzz_0_yzzzzz_1, g_zzzzz_0_yzzzzzz_1, g_zzzzz_0_zzzzzz_1, g_zzzzz_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzz_0_xxxxxxx_0[i] = 7.0 * g_zzzzz_0_xxxxxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxx_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxxy_0[i] = 6.0 * g_zzzzz_0_xxxxxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxxz_0[i] = 6.0 * g_zzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxyy_0[i] = 5.0 * g_zzzzz_0_xxxxyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxyz_0[i] = 5.0 * g_zzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxxzz_0[i] = 5.0 * g_zzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxyyy_0[i] = 4.0 * g_zzzzz_0_xxxyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxyyz_0[i] = 4.0 * g_zzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxyzz_0[i] = 4.0 * g_zzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxxzzz_0[i] = 4.0 * g_zzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyyyy_0[i] = 3.0 * g_zzzzz_0_xxyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyyyz_0[i] = 3.0 * g_zzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyyzz_0[i] = 3.0 * g_zzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxyzzz_0[i] = 3.0 * g_zzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxxzzzz_0[i] = 3.0 * g_zzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyyyy_0[i] = 2.0 * g_zzzzz_0_xyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyyyz_0[i] = 2.0 * g_zzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyyzz_0[i] = 2.0 * g_zzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyyzzz_0[i] = 2.0 * g_zzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxyzzzz_0[i] = 2.0 * g_zzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xxzzzzz_0[i] = 2.0 * g_zzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyyyy_0[i] = g_zzzzz_0_yyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyyyz_0[i] = g_zzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyyzz_0[i] = g_zzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyyzzz_0[i] = g_zzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyyzzzz_0[i] = g_zzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xyzzzzz_0[i] = g_zzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_xzzzzzz_0[i] = g_zzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyyyy_0[i] = g_zzzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyyyz_0[i] = g_zzzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyyzz_0[i] = g_zzzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyyzzz_0[i] = g_zzzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyyzzzz_0[i] = g_zzzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yyzzzzz_0[i] = g_zzzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_yzzzzzz_0[i] = g_zzzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xzzzzz_0_zzzzzzz_0[i] = g_zzzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 756-792 components of targeted buffer : ISK

    auto g_yyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 756);

    auto g_yyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 757);

    auto g_yyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 758);

    auto g_yyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 759);

    auto g_yyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 760);

    auto g_yyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 761);

    auto g_yyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 762);

    auto g_yyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 763);

    auto g_yyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 764);

    auto g_yyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 765);

    auto g_yyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 766);

    auto g_yyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 767);

    auto g_yyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 768);

    auto g_yyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 769);

    auto g_yyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 770);

    auto g_yyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 771);

    auto g_yyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 772);

    auto g_yyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 773);

    auto g_yyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 774);

    auto g_yyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 775);

    auto g_yyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 776);

    auto g_yyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 777);

    auto g_yyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 778);

    auto g_yyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 779);

    auto g_yyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 780);

    auto g_yyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 781);

    auto g_yyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 782);

    auto g_yyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 783);

    auto g_yyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 784);

    auto g_yyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 785);

    auto g_yyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 786);

    auto g_yyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 787);

    auto g_yyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 788);

    auto g_yyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 789);

    auto g_yyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 790);

    auto g_yyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 791);

    #pragma omp simd aligned(g_yyyy_0_xxxxxxx_0, g_yyyy_0_xxxxxxx_1, g_yyyy_0_xxxxxxy_0, g_yyyy_0_xxxxxxy_1, g_yyyy_0_xxxxxxz_0, g_yyyy_0_xxxxxxz_1, g_yyyy_0_xxxxxyy_0, g_yyyy_0_xxxxxyy_1, g_yyyy_0_xxxxxyz_0, g_yyyy_0_xxxxxyz_1, g_yyyy_0_xxxxxzz_0, g_yyyy_0_xxxxxzz_1, g_yyyy_0_xxxxyyy_0, g_yyyy_0_xxxxyyy_1, g_yyyy_0_xxxxyyz_0, g_yyyy_0_xxxxyyz_1, g_yyyy_0_xxxxyzz_0, g_yyyy_0_xxxxyzz_1, g_yyyy_0_xxxxzzz_0, g_yyyy_0_xxxxzzz_1, g_yyyy_0_xxxyyyy_0, g_yyyy_0_xxxyyyy_1, g_yyyy_0_xxxyyyz_0, g_yyyy_0_xxxyyyz_1, g_yyyy_0_xxxyyzz_0, g_yyyy_0_xxxyyzz_1, g_yyyy_0_xxxyzzz_0, g_yyyy_0_xxxyzzz_1, g_yyyy_0_xxxzzzz_0, g_yyyy_0_xxxzzzz_1, g_yyyy_0_xxyyyyy_0, g_yyyy_0_xxyyyyy_1, g_yyyy_0_xxyyyyz_0, g_yyyy_0_xxyyyyz_1, g_yyyy_0_xxyyyzz_0, g_yyyy_0_xxyyyzz_1, g_yyyy_0_xxyyzzz_0, g_yyyy_0_xxyyzzz_1, g_yyyy_0_xxyzzzz_0, g_yyyy_0_xxyzzzz_1, g_yyyy_0_xxzzzzz_0, g_yyyy_0_xxzzzzz_1, g_yyyy_0_xyyyyyy_0, g_yyyy_0_xyyyyyy_1, g_yyyy_0_xyyyyyz_0, g_yyyy_0_xyyyyyz_1, g_yyyy_0_xyyyyzz_0, g_yyyy_0_xyyyyzz_1, g_yyyy_0_xyyyzzz_0, g_yyyy_0_xyyyzzz_1, g_yyyy_0_xyyzzzz_0, g_yyyy_0_xyyzzzz_1, g_yyyy_0_xyzzzzz_0, g_yyyy_0_xyzzzzz_1, g_yyyy_0_xzzzzzz_0, g_yyyy_0_xzzzzzz_1, g_yyyy_0_yyyyyyy_0, g_yyyy_0_yyyyyyy_1, g_yyyy_0_yyyyyyz_0, g_yyyy_0_yyyyyyz_1, g_yyyy_0_yyyyyzz_0, g_yyyy_0_yyyyyzz_1, g_yyyy_0_yyyyzzz_0, g_yyyy_0_yyyyzzz_1, g_yyyy_0_yyyzzzz_0, g_yyyy_0_yyyzzzz_1, g_yyyy_0_yyzzzzz_0, g_yyyy_0_yyzzzzz_1, g_yyyy_0_yzzzzzz_0, g_yyyy_0_yzzzzzz_1, g_yyyy_0_zzzzzzz_0, g_yyyy_0_zzzzzzz_1, g_yyyyy_0_xxxxxx_1, g_yyyyy_0_xxxxxxx_1, g_yyyyy_0_xxxxxxy_1, g_yyyyy_0_xxxxxxz_1, g_yyyyy_0_xxxxxy_1, g_yyyyy_0_xxxxxyy_1, g_yyyyy_0_xxxxxyz_1, g_yyyyy_0_xxxxxz_1, g_yyyyy_0_xxxxxzz_1, g_yyyyy_0_xxxxyy_1, g_yyyyy_0_xxxxyyy_1, g_yyyyy_0_xxxxyyz_1, g_yyyyy_0_xxxxyz_1, g_yyyyy_0_xxxxyzz_1, g_yyyyy_0_xxxxzz_1, g_yyyyy_0_xxxxzzz_1, g_yyyyy_0_xxxyyy_1, g_yyyyy_0_xxxyyyy_1, g_yyyyy_0_xxxyyyz_1, g_yyyyy_0_xxxyyz_1, g_yyyyy_0_xxxyyzz_1, g_yyyyy_0_xxxyzz_1, g_yyyyy_0_xxxyzzz_1, g_yyyyy_0_xxxzzz_1, g_yyyyy_0_xxxzzzz_1, g_yyyyy_0_xxyyyy_1, g_yyyyy_0_xxyyyyy_1, g_yyyyy_0_xxyyyyz_1, g_yyyyy_0_xxyyyz_1, g_yyyyy_0_xxyyyzz_1, g_yyyyy_0_xxyyzz_1, g_yyyyy_0_xxyyzzz_1, g_yyyyy_0_xxyzzz_1, g_yyyyy_0_xxyzzzz_1, g_yyyyy_0_xxzzzz_1, g_yyyyy_0_xxzzzzz_1, g_yyyyy_0_xyyyyy_1, g_yyyyy_0_xyyyyyy_1, g_yyyyy_0_xyyyyyz_1, g_yyyyy_0_xyyyyz_1, g_yyyyy_0_xyyyyzz_1, g_yyyyy_0_xyyyzz_1, g_yyyyy_0_xyyyzzz_1, g_yyyyy_0_xyyzzz_1, g_yyyyy_0_xyyzzzz_1, g_yyyyy_0_xyzzzz_1, g_yyyyy_0_xyzzzzz_1, g_yyyyy_0_xzzzzz_1, g_yyyyy_0_xzzzzzz_1, g_yyyyy_0_yyyyyy_1, g_yyyyy_0_yyyyyyy_1, g_yyyyy_0_yyyyyyz_1, g_yyyyy_0_yyyyyz_1, g_yyyyy_0_yyyyyzz_1, g_yyyyy_0_yyyyzz_1, g_yyyyy_0_yyyyzzz_1, g_yyyyy_0_yyyzzz_1, g_yyyyy_0_yyyzzzz_1, g_yyyyy_0_yyzzzz_1, g_yyyyy_0_yyzzzzz_1, g_yyyyy_0_yzzzzz_1, g_yyyyy_0_yzzzzzz_1, g_yyyyy_0_zzzzzz_1, g_yyyyy_0_zzzzzzz_1, g_yyyyyy_0_xxxxxxx_0, g_yyyyyy_0_xxxxxxy_0, g_yyyyyy_0_xxxxxxz_0, g_yyyyyy_0_xxxxxyy_0, g_yyyyyy_0_xxxxxyz_0, g_yyyyyy_0_xxxxxzz_0, g_yyyyyy_0_xxxxyyy_0, g_yyyyyy_0_xxxxyyz_0, g_yyyyyy_0_xxxxyzz_0, g_yyyyyy_0_xxxxzzz_0, g_yyyyyy_0_xxxyyyy_0, g_yyyyyy_0_xxxyyyz_0, g_yyyyyy_0_xxxyyzz_0, g_yyyyyy_0_xxxyzzz_0, g_yyyyyy_0_xxxzzzz_0, g_yyyyyy_0_xxyyyyy_0, g_yyyyyy_0_xxyyyyz_0, g_yyyyyy_0_xxyyyzz_0, g_yyyyyy_0_xxyyzzz_0, g_yyyyyy_0_xxyzzzz_0, g_yyyyyy_0_xxzzzzz_0, g_yyyyyy_0_xyyyyyy_0, g_yyyyyy_0_xyyyyyz_0, g_yyyyyy_0_xyyyyzz_0, g_yyyyyy_0_xyyyzzz_0, g_yyyyyy_0_xyyzzzz_0, g_yyyyyy_0_xyzzzzz_0, g_yyyyyy_0_xzzzzzz_0, g_yyyyyy_0_yyyyyyy_0, g_yyyyyy_0_yyyyyyz_0, g_yyyyyy_0_yyyyyzz_0, g_yyyyyy_0_yyyyzzz_0, g_yyyyyy_0_yyyzzzz_0, g_yyyyyy_0_yyzzzzz_0, g_yyyyyy_0_yzzzzzz_0, g_yyyyyy_0_zzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyy_0_xxxxxxx_0[i] = 5.0 * g_yyyy_0_xxxxxxx_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxxx_1[i] * fz_be_0 + g_yyyyy_0_xxxxxxx_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxxy_0[i] = 5.0 * g_yyyy_0_xxxxxxy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxxy_1[i] * fz_be_0 + g_yyyyy_0_xxxxxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxxz_0[i] = 5.0 * g_yyyy_0_xxxxxxz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxxz_1[i] * fz_be_0 + g_yyyyy_0_xxxxxxz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxyy_0[i] = 5.0 * g_yyyy_0_xxxxxyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxyy_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxxxxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxyz_0[i] = 5.0 * g_yyyy_0_xxxxxyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxyz_1[i] * fz_be_0 + g_yyyyy_0_xxxxxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxxzz_0[i] = 5.0 * g_yyyy_0_xxxxxzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxxzz_1[i] * fz_be_0 + g_yyyyy_0_xxxxxzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxyyy_0[i] = 5.0 * g_yyyy_0_xxxxyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxyyy_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xxxxyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxyyz_0[i] = 5.0 * g_yyyy_0_xxxxyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxxxyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxyzz_0[i] = 5.0 * g_yyyy_0_xxxxyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxyzz_1[i] * fz_be_0 + g_yyyyy_0_xxxxzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxxzzz_0[i] = 5.0 * g_yyyy_0_xxxxzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxxzzz_1[i] * fz_be_0 + g_yyyyy_0_xxxxzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyyyy_0[i] = 5.0 * g_yyyy_0_xxxyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyyyy_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_xxxyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyyyz_0[i] = 5.0 * g_yyyy_0_xxxyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xxxyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyyzz_0[i] = 5.0 * g_yyyy_0_xxxyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxxyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxyzzz_0[i] = 5.0 * g_yyyy_0_xxxyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxyzzz_1[i] * fz_be_0 + g_yyyyy_0_xxxzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxxzzzz_0[i] = 5.0 * g_yyyy_0_xxxzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxxzzzz_1[i] * fz_be_0 + g_yyyyy_0_xxxzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyyyy_0[i] = 5.0 * g_yyyy_0_xxyyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyyyy_1[i] * fz_be_0 + 5.0 * g_yyyyy_0_xxyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyyyz_0[i] = 5.0 * g_yyyy_0_xxyyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_xxyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyyzz_0[i] = 5.0 * g_yyyy_0_xxyyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xxyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyyzzz_0[i] = 5.0 * g_yyyy_0_xxyyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xxyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxyzzzz_0[i] = 5.0 * g_yyyy_0_xxyzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxyzzzz_1[i] * fz_be_0 + g_yyyyy_0_xxzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xxzzzzz_0[i] = 5.0 * g_yyyy_0_xxzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xxzzzzz_1[i] * fz_be_0 + g_yyyyy_0_xxzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyyyy_0[i] = 5.0 * g_yyyy_0_xyyyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyyyy_1[i] * fz_be_0 + 6.0 * g_yyyyy_0_xyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyyyz_0[i] = 5.0 * g_yyyy_0_xyyyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyyy_0_xyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyyzz_0[i] = 5.0 * g_yyyy_0_xyyyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_xyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyyzzz_0[i] = 5.0 * g_yyyy_0_xyyyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_xyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyyzzzz_0[i] = 5.0 * g_yyyy_0_xyyzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_xyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xyzzzzz_0[i] = 5.0 * g_yyyy_0_xyzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xyzzzzz_1[i] * fz_be_0 + g_yyyyy_0_xzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_xzzzzzz_0[i] = 5.0 * g_yyyy_0_xzzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_xzzzzzz_1[i] * fz_be_0 + g_yyyyy_0_xzzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyyyy_0[i] = 5.0 * g_yyyy_0_yyyyyyy_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyyyy_1[i] * fz_be_0 + 7.0 * g_yyyyy_0_yyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyyy_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyyyz_0[i] = 5.0 * g_yyyy_0_yyyyyyz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyyyy_0_yyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyyz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyyzz_0[i] = 5.0 * g_yyyy_0_yyyyyzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyyyy_0_yyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyyzzz_0[i] = 5.0 * g_yyyy_0_yyyyzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyyyy_0_yyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyyzzzz_0[i] = 5.0 * g_yyyy_0_yyyzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyyyy_0_yyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yyzzzzz_0[i] = 5.0 * g_yyyy_0_yyzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyy_0_yzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_yzzzzzz_0[i] = 5.0 * g_yyyy_0_yzzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_yzzzzzz_1[i] * fz_be_0 + g_yyyyy_0_zzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yzzzzzz_1[i] * wa_y[i];

        g_yyyyyy_0_zzzzzzz_0[i] = 5.0 * g_yyyy_0_zzzzzzz_0[i] * fbe_0 - 5.0 * g_yyyy_0_zzzzzzz_1[i] * fz_be_0 + g_yyyyy_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 792-828 components of targeted buffer : ISK

    auto g_yyyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 792);

    auto g_yyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 793);

    auto g_yyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 794);

    auto g_yyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 795);

    auto g_yyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 796);

    auto g_yyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 797);

    auto g_yyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 798);

    auto g_yyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 799);

    auto g_yyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 800);

    auto g_yyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 801);

    auto g_yyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 802);

    auto g_yyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 803);

    auto g_yyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 804);

    auto g_yyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 805);

    auto g_yyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 806);

    auto g_yyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 807);

    auto g_yyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 808);

    auto g_yyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 809);

    auto g_yyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 810);

    auto g_yyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 811);

    auto g_yyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 812);

    auto g_yyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 813);

    auto g_yyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 814);

    auto g_yyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 815);

    auto g_yyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 816);

    auto g_yyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 817);

    auto g_yyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 818);

    auto g_yyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 819);

    auto g_yyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 820);

    auto g_yyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 821);

    auto g_yyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 822);

    auto g_yyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 823);

    auto g_yyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 824);

    auto g_yyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 825);

    auto g_yyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 826);

    auto g_yyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 827);

    #pragma omp simd aligned(g_yyyyy_0_xxxxxx_1, g_yyyyy_0_xxxxxxx_1, g_yyyyy_0_xxxxxxy_1, g_yyyyy_0_xxxxxxz_1, g_yyyyy_0_xxxxxy_1, g_yyyyy_0_xxxxxyy_1, g_yyyyy_0_xxxxxyz_1, g_yyyyy_0_xxxxxz_1, g_yyyyy_0_xxxxxzz_1, g_yyyyy_0_xxxxyy_1, g_yyyyy_0_xxxxyyy_1, g_yyyyy_0_xxxxyyz_1, g_yyyyy_0_xxxxyz_1, g_yyyyy_0_xxxxyzz_1, g_yyyyy_0_xxxxzz_1, g_yyyyy_0_xxxxzzz_1, g_yyyyy_0_xxxyyy_1, g_yyyyy_0_xxxyyyy_1, g_yyyyy_0_xxxyyyz_1, g_yyyyy_0_xxxyyz_1, g_yyyyy_0_xxxyyzz_1, g_yyyyy_0_xxxyzz_1, g_yyyyy_0_xxxyzzz_1, g_yyyyy_0_xxxzzz_1, g_yyyyy_0_xxxzzzz_1, g_yyyyy_0_xxyyyy_1, g_yyyyy_0_xxyyyyy_1, g_yyyyy_0_xxyyyyz_1, g_yyyyy_0_xxyyyz_1, g_yyyyy_0_xxyyyzz_1, g_yyyyy_0_xxyyzz_1, g_yyyyy_0_xxyyzzz_1, g_yyyyy_0_xxyzzz_1, g_yyyyy_0_xxyzzzz_1, g_yyyyy_0_xxzzzz_1, g_yyyyy_0_xxzzzzz_1, g_yyyyy_0_xyyyyy_1, g_yyyyy_0_xyyyyyy_1, g_yyyyy_0_xyyyyyz_1, g_yyyyy_0_xyyyyz_1, g_yyyyy_0_xyyyyzz_1, g_yyyyy_0_xyyyzz_1, g_yyyyy_0_xyyyzzz_1, g_yyyyy_0_xyyzzz_1, g_yyyyy_0_xyyzzzz_1, g_yyyyy_0_xyzzzz_1, g_yyyyy_0_xyzzzzz_1, g_yyyyy_0_xzzzzz_1, g_yyyyy_0_xzzzzzz_1, g_yyyyy_0_yyyyyy_1, g_yyyyy_0_yyyyyyy_1, g_yyyyy_0_yyyyyyz_1, g_yyyyy_0_yyyyyz_1, g_yyyyy_0_yyyyyzz_1, g_yyyyy_0_yyyyzz_1, g_yyyyy_0_yyyyzzz_1, g_yyyyy_0_yyyzzz_1, g_yyyyy_0_yyyzzzz_1, g_yyyyy_0_yyzzzz_1, g_yyyyy_0_yyzzzzz_1, g_yyyyy_0_yzzzzz_1, g_yyyyy_0_yzzzzzz_1, g_yyyyy_0_zzzzzz_1, g_yyyyy_0_zzzzzzz_1, g_yyyyyz_0_xxxxxxx_0, g_yyyyyz_0_xxxxxxy_0, g_yyyyyz_0_xxxxxxz_0, g_yyyyyz_0_xxxxxyy_0, g_yyyyyz_0_xxxxxyz_0, g_yyyyyz_0_xxxxxzz_0, g_yyyyyz_0_xxxxyyy_0, g_yyyyyz_0_xxxxyyz_0, g_yyyyyz_0_xxxxyzz_0, g_yyyyyz_0_xxxxzzz_0, g_yyyyyz_0_xxxyyyy_0, g_yyyyyz_0_xxxyyyz_0, g_yyyyyz_0_xxxyyzz_0, g_yyyyyz_0_xxxyzzz_0, g_yyyyyz_0_xxxzzzz_0, g_yyyyyz_0_xxyyyyy_0, g_yyyyyz_0_xxyyyyz_0, g_yyyyyz_0_xxyyyzz_0, g_yyyyyz_0_xxyyzzz_0, g_yyyyyz_0_xxyzzzz_0, g_yyyyyz_0_xxzzzzz_0, g_yyyyyz_0_xyyyyyy_0, g_yyyyyz_0_xyyyyyz_0, g_yyyyyz_0_xyyyyzz_0, g_yyyyyz_0_xyyyzzz_0, g_yyyyyz_0_xyyzzzz_0, g_yyyyyz_0_xyzzzzz_0, g_yyyyyz_0_xzzzzzz_0, g_yyyyyz_0_yyyyyyy_0, g_yyyyyz_0_yyyyyyz_0, g_yyyyyz_0_yyyyyzz_0, g_yyyyyz_0_yyyyzzz_0, g_yyyyyz_0_yyyzzzz_0, g_yyyyyz_0_yyzzzzz_0, g_yyyyyz_0_yzzzzzz_0, g_yyyyyz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyz_0_xxxxxxx_0[i] = g_yyyyy_0_xxxxxxx_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxxy_0[i] = g_yyyyy_0_xxxxxxy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxxz_0[i] = g_yyyyy_0_xxxxxx_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxxz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxyy_0[i] = g_yyyyy_0_xxxxxyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxyz_0[i] = g_yyyyy_0_xxxxxy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxxzz_0[i] = 2.0 * g_yyyyy_0_xxxxxz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxxzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxyyy_0[i] = g_yyyyy_0_xxxxyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxyyz_0[i] = g_yyyyy_0_xxxxyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxyzz_0[i] = 2.0 * g_yyyyy_0_xxxxyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxxzzz_0[i] = 3.0 * g_yyyyy_0_xxxxzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxxzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyyyy_0[i] = g_yyyyy_0_xxxyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyyyz_0[i] = g_yyyyy_0_xxxyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyyzz_0[i] = 2.0 * g_yyyyy_0_xxxyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxyzzz_0[i] = 3.0 * g_yyyyy_0_xxxyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxxzzzz_0[i] = 4.0 * g_yyyyy_0_xxxzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxxzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyyyy_0[i] = g_yyyyy_0_xxyyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyyyz_0[i] = g_yyyyy_0_xxyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyyzz_0[i] = 2.0 * g_yyyyy_0_xxyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyyzzz_0[i] = 3.0 * g_yyyyy_0_xxyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxyzzzz_0[i] = 4.0 * g_yyyyy_0_xxyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxyzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xxzzzzz_0[i] = 5.0 * g_yyyyy_0_xxzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xxzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyyyy_0[i] = g_yyyyy_0_xyyyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyyyz_0[i] = g_yyyyy_0_xyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyyzz_0[i] = 2.0 * g_yyyyy_0_xyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyyzzz_0[i] = 3.0 * g_yyyyy_0_xyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyyzzzz_0[i] = 4.0 * g_yyyyy_0_xyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyyzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xyzzzzz_0[i] = 5.0 * g_yyyyy_0_xyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xyzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_xzzzzzz_0[i] = 6.0 * g_yyyyy_0_xzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_xzzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyyyy_0[i] = g_yyyyy_0_yyyyyyy_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyyyz_0[i] = g_yyyyy_0_yyyyyy_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyyz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyyzz_0[i] = 2.0 * g_yyyyy_0_yyyyyz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyyzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyyzzz_0[i] = 3.0 * g_yyyyy_0_yyyyzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyyzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyyzzzz_0[i] = 4.0 * g_yyyyy_0_yyyzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyyzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yyzzzzz_0[i] = 5.0 * g_yyyyy_0_yyzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yyzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_yzzzzzz_0[i] = 6.0 * g_yyyyy_0_yzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_yzzzzzz_1[i] * wa_z[i];

        g_yyyyyz_0_zzzzzzz_0[i] = 7.0 * g_yyyyy_0_zzzzzz_1[i] * fi_acd_0 + g_yyyyy_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 828-864 components of targeted buffer : ISK

    auto g_yyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 828);

    auto g_yyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 829);

    auto g_yyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 830);

    auto g_yyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 831);

    auto g_yyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 832);

    auto g_yyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 833);

    auto g_yyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 834);

    auto g_yyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 835);

    auto g_yyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 836);

    auto g_yyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 837);

    auto g_yyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 838);

    auto g_yyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 839);

    auto g_yyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 840);

    auto g_yyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 841);

    auto g_yyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 842);

    auto g_yyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 843);

    auto g_yyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 844);

    auto g_yyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 845);

    auto g_yyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 846);

    auto g_yyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 847);

    auto g_yyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 848);

    auto g_yyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 849);

    auto g_yyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 850);

    auto g_yyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 851);

    auto g_yyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 852);

    auto g_yyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 853);

    auto g_yyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 854);

    auto g_yyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 855);

    auto g_yyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 856);

    auto g_yyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 857);

    auto g_yyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 858);

    auto g_yyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 859);

    auto g_yyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 860);

    auto g_yyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 861);

    auto g_yyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 862);

    auto g_yyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 863);

    #pragma omp simd aligned(g_yyyy_0_xxxxxxy_0, g_yyyy_0_xxxxxxy_1, g_yyyy_0_xxxxxyy_0, g_yyyy_0_xxxxxyy_1, g_yyyy_0_xxxxyyy_0, g_yyyy_0_xxxxyyy_1, g_yyyy_0_xxxyyyy_0, g_yyyy_0_xxxyyyy_1, g_yyyy_0_xxyyyyy_0, g_yyyy_0_xxyyyyy_1, g_yyyy_0_xyyyyyy_0, g_yyyy_0_xyyyyyy_1, g_yyyy_0_yyyyyyy_0, g_yyyy_0_yyyyyyy_1, g_yyyyz_0_xxxxxxy_1, g_yyyyz_0_xxxxxyy_1, g_yyyyz_0_xxxxyyy_1, g_yyyyz_0_xxxyyyy_1, g_yyyyz_0_xxyyyyy_1, g_yyyyz_0_xyyyyyy_1, g_yyyyz_0_yyyyyyy_1, g_yyyyzz_0_xxxxxxx_0, g_yyyyzz_0_xxxxxxy_0, g_yyyyzz_0_xxxxxxz_0, g_yyyyzz_0_xxxxxyy_0, g_yyyyzz_0_xxxxxyz_0, g_yyyyzz_0_xxxxxzz_0, g_yyyyzz_0_xxxxyyy_0, g_yyyyzz_0_xxxxyyz_0, g_yyyyzz_0_xxxxyzz_0, g_yyyyzz_0_xxxxzzz_0, g_yyyyzz_0_xxxyyyy_0, g_yyyyzz_0_xxxyyyz_0, g_yyyyzz_0_xxxyyzz_0, g_yyyyzz_0_xxxyzzz_0, g_yyyyzz_0_xxxzzzz_0, g_yyyyzz_0_xxyyyyy_0, g_yyyyzz_0_xxyyyyz_0, g_yyyyzz_0_xxyyyzz_0, g_yyyyzz_0_xxyyzzz_0, g_yyyyzz_0_xxyzzzz_0, g_yyyyzz_0_xxzzzzz_0, g_yyyyzz_0_xyyyyyy_0, g_yyyyzz_0_xyyyyyz_0, g_yyyyzz_0_xyyyyzz_0, g_yyyyzz_0_xyyyzzz_0, g_yyyyzz_0_xyyzzzz_0, g_yyyyzz_0_xyzzzzz_0, g_yyyyzz_0_xzzzzzz_0, g_yyyyzz_0_yyyyyyy_0, g_yyyyzz_0_yyyyyyz_0, g_yyyyzz_0_yyyyyzz_0, g_yyyyzz_0_yyyyzzz_0, g_yyyyzz_0_yyyzzzz_0, g_yyyyzz_0_yyzzzzz_0, g_yyyyzz_0_yzzzzzz_0, g_yyyyzz_0_zzzzzzz_0, g_yyyzz_0_xxxxxxx_1, g_yyyzz_0_xxxxxxz_1, g_yyyzz_0_xxxxxyz_1, g_yyyzz_0_xxxxxz_1, g_yyyzz_0_xxxxxzz_1, g_yyyzz_0_xxxxyyz_1, g_yyyzz_0_xxxxyz_1, g_yyyzz_0_xxxxyzz_1, g_yyyzz_0_xxxxzz_1, g_yyyzz_0_xxxxzzz_1, g_yyyzz_0_xxxyyyz_1, g_yyyzz_0_xxxyyz_1, g_yyyzz_0_xxxyyzz_1, g_yyyzz_0_xxxyzz_1, g_yyyzz_0_xxxyzzz_1, g_yyyzz_0_xxxzzz_1, g_yyyzz_0_xxxzzzz_1, g_yyyzz_0_xxyyyyz_1, g_yyyzz_0_xxyyyz_1, g_yyyzz_0_xxyyyzz_1, g_yyyzz_0_xxyyzz_1, g_yyyzz_0_xxyyzzz_1, g_yyyzz_0_xxyzzz_1, g_yyyzz_0_xxyzzzz_1, g_yyyzz_0_xxzzzz_1, g_yyyzz_0_xxzzzzz_1, g_yyyzz_0_xyyyyyz_1, g_yyyzz_0_xyyyyz_1, g_yyyzz_0_xyyyyzz_1, g_yyyzz_0_xyyyzz_1, g_yyyzz_0_xyyyzzz_1, g_yyyzz_0_xyyzzz_1, g_yyyzz_0_xyyzzzz_1, g_yyyzz_0_xyzzzz_1, g_yyyzz_0_xyzzzzz_1, g_yyyzz_0_xzzzzz_1, g_yyyzz_0_xzzzzzz_1, g_yyyzz_0_yyyyyyz_1, g_yyyzz_0_yyyyyz_1, g_yyyzz_0_yyyyyzz_1, g_yyyzz_0_yyyyzz_1, g_yyyzz_0_yyyyzzz_1, g_yyyzz_0_yyyzzz_1, g_yyyzz_0_yyyzzzz_1, g_yyyzz_0_yyzzzz_1, g_yyyzz_0_yyzzzzz_1, g_yyyzz_0_yzzzzz_1, g_yyyzz_0_yzzzzzz_1, g_yyyzz_0_zzzzzz_1, g_yyyzz_0_zzzzzzz_1, g_yyzz_0_xxxxxxx_0, g_yyzz_0_xxxxxxx_1, g_yyzz_0_xxxxxxz_0, g_yyzz_0_xxxxxxz_1, g_yyzz_0_xxxxxyz_0, g_yyzz_0_xxxxxyz_1, g_yyzz_0_xxxxxzz_0, g_yyzz_0_xxxxxzz_1, g_yyzz_0_xxxxyyz_0, g_yyzz_0_xxxxyyz_1, g_yyzz_0_xxxxyzz_0, g_yyzz_0_xxxxyzz_1, g_yyzz_0_xxxxzzz_0, g_yyzz_0_xxxxzzz_1, g_yyzz_0_xxxyyyz_0, g_yyzz_0_xxxyyyz_1, g_yyzz_0_xxxyyzz_0, g_yyzz_0_xxxyyzz_1, g_yyzz_0_xxxyzzz_0, g_yyzz_0_xxxyzzz_1, g_yyzz_0_xxxzzzz_0, g_yyzz_0_xxxzzzz_1, g_yyzz_0_xxyyyyz_0, g_yyzz_0_xxyyyyz_1, g_yyzz_0_xxyyyzz_0, g_yyzz_0_xxyyyzz_1, g_yyzz_0_xxyyzzz_0, g_yyzz_0_xxyyzzz_1, g_yyzz_0_xxyzzzz_0, g_yyzz_0_xxyzzzz_1, g_yyzz_0_xxzzzzz_0, g_yyzz_0_xxzzzzz_1, g_yyzz_0_xyyyyyz_0, g_yyzz_0_xyyyyyz_1, g_yyzz_0_xyyyyzz_0, g_yyzz_0_xyyyyzz_1, g_yyzz_0_xyyyzzz_0, g_yyzz_0_xyyyzzz_1, g_yyzz_0_xyyzzzz_0, g_yyzz_0_xyyzzzz_1, g_yyzz_0_xyzzzzz_0, g_yyzz_0_xyzzzzz_1, g_yyzz_0_xzzzzzz_0, g_yyzz_0_xzzzzzz_1, g_yyzz_0_yyyyyyz_0, g_yyzz_0_yyyyyyz_1, g_yyzz_0_yyyyyzz_0, g_yyzz_0_yyyyyzz_1, g_yyzz_0_yyyyzzz_0, g_yyzz_0_yyyyzzz_1, g_yyzz_0_yyyzzzz_0, g_yyzz_0_yyyzzzz_1, g_yyzz_0_yyzzzzz_0, g_yyzz_0_yyzzzzz_1, g_yyzz_0_yzzzzzz_0, g_yyzz_0_yzzzzzz_1, g_yyzz_0_zzzzzzz_0, g_yyzz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzz_0_xxxxxxx_0[i] = 3.0 * g_yyzz_0_xxxxxxx_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxxx_1[i] * fz_be_0 + g_yyyzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxxxy_0[i] = g_yyyy_0_xxxxxxy_0[i] * fbe_0 - g_yyyy_0_xxxxxxy_1[i] * fz_be_0 + g_yyyyz_0_xxxxxxy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxxxxz_0[i] = 3.0 * g_yyzz_0_xxxxxxz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxxz_1[i] * fz_be_0 + g_yyyzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxxyy_0[i] = g_yyyy_0_xxxxxyy_0[i] * fbe_0 - g_yyyy_0_xxxxxyy_1[i] * fz_be_0 + g_yyyyz_0_xxxxxyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxxxyz_0[i] = 3.0 * g_yyzz_0_xxxxxyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxyz_1[i] * fz_be_0 + g_yyyzz_0_xxxxxz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxxzz_0[i] = 3.0 * g_yyzz_0_xxxxxzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxzz_1[i] * fz_be_0 + g_yyyzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxyyy_0[i] = g_yyyy_0_xxxxyyy_0[i] * fbe_0 - g_yyyy_0_xxxxyyy_1[i] * fz_be_0 + g_yyyyz_0_xxxxyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxxyyz_0[i] = 3.0 * g_yyzz_0_xxxxyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xxxxyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxyzz_0[i] = 3.0 * g_yyzz_0_xxxxyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxyzz_1[i] * fz_be_0 + g_yyyzz_0_xxxxzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxxzzz_0[i] = 3.0 * g_yyzz_0_xxxxzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxzzz_1[i] * fz_be_0 + g_yyyzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxyyyy_0[i] = g_yyyy_0_xxxyyyy_0[i] * fbe_0 - g_yyyy_0_xxxyyyy_1[i] * fz_be_0 + g_yyyyz_0_xxxyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxxyyyz_0[i] = 3.0 * g_yyzz_0_xxxyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_xxxyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxyyzz_0[i] = 3.0 * g_yyzz_0_xxxyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xxxyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxyzzz_0[i] = 3.0 * g_yyzz_0_xxxyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyzzz_1[i] * fz_be_0 + g_yyyzz_0_xxxzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxxzzzz_0[i] = 3.0 * g_yyzz_0_xxxzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxzzzz_1[i] * fz_be_0 + g_yyyzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyyyyy_0[i] = g_yyyy_0_xxyyyyy_0[i] * fbe_0 - g_yyyy_0_xxyyyyy_1[i] * fz_be_0 + g_yyyyz_0_xxyyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xxyyyyz_0[i] = 3.0 * g_yyzz_0_xxyyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyzz_0_xxyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyyyzz_0[i] = 3.0 * g_yyzz_0_xxyyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_xxyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyyzzz_0[i] = 3.0 * g_yyzz_0_xxyyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xxyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxyzzzz_0[i] = 3.0 * g_yyzz_0_xxyzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyzzzz_1[i] * fz_be_0 + g_yyyzz_0_xxzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xxzzzzz_0[i] = 3.0 * g_yyzz_0_xxzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxzzzzz_1[i] * fz_be_0 + g_yyyzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyyyyy_0[i] = g_yyyy_0_xyyyyyy_0[i] * fbe_0 - g_yyyy_0_xyyyyyy_1[i] * fz_be_0 + g_yyyyz_0_xyyyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_xyyyyyz_0[i] = 3.0 * g_yyzz_0_xyyyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyzz_0_xyyyyz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyyyzz_0[i] = 3.0 * g_yyzz_0_xyyyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyzz_0_xyyyzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyyzzz_0[i] = 3.0 * g_yyzz_0_xyyyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_xyyzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyyzzzz_0[i] = 3.0 * g_yyzz_0_xyyzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_xyzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xyzzzzz_0[i] = 3.0 * g_yyzz_0_xyzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyzzzzz_1[i] * fz_be_0 + g_yyyzz_0_xzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_xzzzzzz_0[i] = 3.0 * g_yyzz_0_xzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_xzzzzzz_1[i] * fz_be_0 + g_yyyzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyyyyy_0[i] = g_yyyy_0_yyyyyyy_0[i] * fbe_0 - g_yyyy_0_yyyyyyy_1[i] * fz_be_0 + g_yyyyz_0_yyyyyyy_1[i] * wa_z[i];

        g_yyyyzz_0_yyyyyyz_0[i] = 3.0 * g_yyzz_0_yyyyyyz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyyzz_0_yyyyyz_1[i] * fi_acd_0 + g_yyyzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyyyzz_0[i] = 3.0 * g_yyzz_0_yyyyyzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyyzz_0_yyyyzz_1[i] * fi_acd_0 + g_yyyzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyyzzz_0[i] = 3.0 * g_yyzz_0_yyyyzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyyzz_0_yyyzzz_1[i] * fi_acd_0 + g_yyyzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyyzzzz_0[i] = 3.0 * g_yyzz_0_yyyzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyyzz_0_yyzzzz_1[i] * fi_acd_0 + g_yyyzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yyzzzzz_0[i] = 3.0 * g_yyzz_0_yyzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzz_0_yzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_yzzzzzz_0[i] = 3.0 * g_yyzz_0_yzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_yzzzzzz_1[i] * fz_be_0 + g_yyyzz_0_zzzzzz_1[i] * fi_acd_0 + g_yyyzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yyyyzz_0_zzzzzzz_0[i] = 3.0 * g_yyzz_0_zzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzz_0_zzzzzzz_1[i] * fz_be_0 + g_yyyzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 864-900 components of targeted buffer : ISK

    auto g_yyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 864);

    auto g_yyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 865);

    auto g_yyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 866);

    auto g_yyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 867);

    auto g_yyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 868);

    auto g_yyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 869);

    auto g_yyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 870);

    auto g_yyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 871);

    auto g_yyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 872);

    auto g_yyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 873);

    auto g_yyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 874);

    auto g_yyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 875);

    auto g_yyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 876);

    auto g_yyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 877);

    auto g_yyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 878);

    auto g_yyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 879);

    auto g_yyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 880);

    auto g_yyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 881);

    auto g_yyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 882);

    auto g_yyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 883);

    auto g_yyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 884);

    auto g_yyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 885);

    auto g_yyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 886);

    auto g_yyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 887);

    auto g_yyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 888);

    auto g_yyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 889);

    auto g_yyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 890);

    auto g_yyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 891);

    auto g_yyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 892);

    auto g_yyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 893);

    auto g_yyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 894);

    auto g_yyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 895);

    auto g_yyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 896);

    auto g_yyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 897);

    auto g_yyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 898);

    auto g_yyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 899);

    #pragma omp simd aligned(g_yyyz_0_xxxxxxy_0, g_yyyz_0_xxxxxxy_1, g_yyyz_0_xxxxxyy_0, g_yyyz_0_xxxxxyy_1, g_yyyz_0_xxxxyyy_0, g_yyyz_0_xxxxyyy_1, g_yyyz_0_xxxyyyy_0, g_yyyz_0_xxxyyyy_1, g_yyyz_0_xxyyyyy_0, g_yyyz_0_xxyyyyy_1, g_yyyz_0_xyyyyyy_0, g_yyyz_0_xyyyyyy_1, g_yyyz_0_yyyyyyy_0, g_yyyz_0_yyyyyyy_1, g_yyyzz_0_xxxxxxy_1, g_yyyzz_0_xxxxxyy_1, g_yyyzz_0_xxxxyyy_1, g_yyyzz_0_xxxyyyy_1, g_yyyzz_0_xxyyyyy_1, g_yyyzz_0_xyyyyyy_1, g_yyyzz_0_yyyyyyy_1, g_yyyzzz_0_xxxxxxx_0, g_yyyzzz_0_xxxxxxy_0, g_yyyzzz_0_xxxxxxz_0, g_yyyzzz_0_xxxxxyy_0, g_yyyzzz_0_xxxxxyz_0, g_yyyzzz_0_xxxxxzz_0, g_yyyzzz_0_xxxxyyy_0, g_yyyzzz_0_xxxxyyz_0, g_yyyzzz_0_xxxxyzz_0, g_yyyzzz_0_xxxxzzz_0, g_yyyzzz_0_xxxyyyy_0, g_yyyzzz_0_xxxyyyz_0, g_yyyzzz_0_xxxyyzz_0, g_yyyzzz_0_xxxyzzz_0, g_yyyzzz_0_xxxzzzz_0, g_yyyzzz_0_xxyyyyy_0, g_yyyzzz_0_xxyyyyz_0, g_yyyzzz_0_xxyyyzz_0, g_yyyzzz_0_xxyyzzz_0, g_yyyzzz_0_xxyzzzz_0, g_yyyzzz_0_xxzzzzz_0, g_yyyzzz_0_xyyyyyy_0, g_yyyzzz_0_xyyyyyz_0, g_yyyzzz_0_xyyyyzz_0, g_yyyzzz_0_xyyyzzz_0, g_yyyzzz_0_xyyzzzz_0, g_yyyzzz_0_xyzzzzz_0, g_yyyzzz_0_xzzzzzz_0, g_yyyzzz_0_yyyyyyy_0, g_yyyzzz_0_yyyyyyz_0, g_yyyzzz_0_yyyyyzz_0, g_yyyzzz_0_yyyyzzz_0, g_yyyzzz_0_yyyzzzz_0, g_yyyzzz_0_yyzzzzz_0, g_yyyzzz_0_yzzzzzz_0, g_yyyzzz_0_zzzzzzz_0, g_yyzzz_0_xxxxxxx_1, g_yyzzz_0_xxxxxxz_1, g_yyzzz_0_xxxxxyz_1, g_yyzzz_0_xxxxxz_1, g_yyzzz_0_xxxxxzz_1, g_yyzzz_0_xxxxyyz_1, g_yyzzz_0_xxxxyz_1, g_yyzzz_0_xxxxyzz_1, g_yyzzz_0_xxxxzz_1, g_yyzzz_0_xxxxzzz_1, g_yyzzz_0_xxxyyyz_1, g_yyzzz_0_xxxyyz_1, g_yyzzz_0_xxxyyzz_1, g_yyzzz_0_xxxyzz_1, g_yyzzz_0_xxxyzzz_1, g_yyzzz_0_xxxzzz_1, g_yyzzz_0_xxxzzzz_1, g_yyzzz_0_xxyyyyz_1, g_yyzzz_0_xxyyyz_1, g_yyzzz_0_xxyyyzz_1, g_yyzzz_0_xxyyzz_1, g_yyzzz_0_xxyyzzz_1, g_yyzzz_0_xxyzzz_1, g_yyzzz_0_xxyzzzz_1, g_yyzzz_0_xxzzzz_1, g_yyzzz_0_xxzzzzz_1, g_yyzzz_0_xyyyyyz_1, g_yyzzz_0_xyyyyz_1, g_yyzzz_0_xyyyyzz_1, g_yyzzz_0_xyyyzz_1, g_yyzzz_0_xyyyzzz_1, g_yyzzz_0_xyyzzz_1, g_yyzzz_0_xyyzzzz_1, g_yyzzz_0_xyzzzz_1, g_yyzzz_0_xyzzzzz_1, g_yyzzz_0_xzzzzz_1, g_yyzzz_0_xzzzzzz_1, g_yyzzz_0_yyyyyyz_1, g_yyzzz_0_yyyyyz_1, g_yyzzz_0_yyyyyzz_1, g_yyzzz_0_yyyyzz_1, g_yyzzz_0_yyyyzzz_1, g_yyzzz_0_yyyzzz_1, g_yyzzz_0_yyyzzzz_1, g_yyzzz_0_yyzzzz_1, g_yyzzz_0_yyzzzzz_1, g_yyzzz_0_yzzzzz_1, g_yyzzz_0_yzzzzzz_1, g_yyzzz_0_zzzzzz_1, g_yyzzz_0_zzzzzzz_1, g_yzzz_0_xxxxxxx_0, g_yzzz_0_xxxxxxx_1, g_yzzz_0_xxxxxxz_0, g_yzzz_0_xxxxxxz_1, g_yzzz_0_xxxxxyz_0, g_yzzz_0_xxxxxyz_1, g_yzzz_0_xxxxxzz_0, g_yzzz_0_xxxxxzz_1, g_yzzz_0_xxxxyyz_0, g_yzzz_0_xxxxyyz_1, g_yzzz_0_xxxxyzz_0, g_yzzz_0_xxxxyzz_1, g_yzzz_0_xxxxzzz_0, g_yzzz_0_xxxxzzz_1, g_yzzz_0_xxxyyyz_0, g_yzzz_0_xxxyyyz_1, g_yzzz_0_xxxyyzz_0, g_yzzz_0_xxxyyzz_1, g_yzzz_0_xxxyzzz_0, g_yzzz_0_xxxyzzz_1, g_yzzz_0_xxxzzzz_0, g_yzzz_0_xxxzzzz_1, g_yzzz_0_xxyyyyz_0, g_yzzz_0_xxyyyyz_1, g_yzzz_0_xxyyyzz_0, g_yzzz_0_xxyyyzz_1, g_yzzz_0_xxyyzzz_0, g_yzzz_0_xxyyzzz_1, g_yzzz_0_xxyzzzz_0, g_yzzz_0_xxyzzzz_1, g_yzzz_0_xxzzzzz_0, g_yzzz_0_xxzzzzz_1, g_yzzz_0_xyyyyyz_0, g_yzzz_0_xyyyyyz_1, g_yzzz_0_xyyyyzz_0, g_yzzz_0_xyyyyzz_1, g_yzzz_0_xyyyzzz_0, g_yzzz_0_xyyyzzz_1, g_yzzz_0_xyyzzzz_0, g_yzzz_0_xyyzzzz_1, g_yzzz_0_xyzzzzz_0, g_yzzz_0_xyzzzzz_1, g_yzzz_0_xzzzzzz_0, g_yzzz_0_xzzzzzz_1, g_yzzz_0_yyyyyyz_0, g_yzzz_0_yyyyyyz_1, g_yzzz_0_yyyyyzz_0, g_yzzz_0_yyyyyzz_1, g_yzzz_0_yyyyzzz_0, g_yzzz_0_yyyyzzz_1, g_yzzz_0_yyyzzzz_0, g_yzzz_0_yyyzzzz_1, g_yzzz_0_yyzzzzz_0, g_yzzz_0_yyzzzzz_1, g_yzzz_0_yzzzzzz_0, g_yzzz_0_yzzzzzz_1, g_yzzz_0_zzzzzzz_0, g_yzzz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzz_0_xxxxxxx_0[i] = 2.0 * g_yzzz_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxxx_1[i] * fz_be_0 + g_yyzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxxxy_0[i] = 2.0 * g_yyyz_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxxxxy_1[i] * fz_be_0 + g_yyyzz_0_xxxxxxy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxxxxz_0[i] = 2.0 * g_yzzz_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxxz_1[i] * fz_be_0 + g_yyzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxxyy_0[i] = 2.0 * g_yyyz_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxxxyy_1[i] * fz_be_0 + g_yyyzz_0_xxxxxyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxxxyz_0[i] = 2.0 * g_yzzz_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxyz_1[i] * fz_be_0 + g_yyzzz_0_xxxxxz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxxzz_0[i] = 2.0 * g_yzzz_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxxzz_1[i] * fz_be_0 + g_yyzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxyyy_0[i] = 2.0 * g_yyyz_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxxyyy_1[i] * fz_be_0 + g_yyyzz_0_xxxxyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxxyyz_0[i] = 2.0 * g_yzzz_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xxxxyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxyzz_0[i] = 2.0 * g_yzzz_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxyzz_1[i] * fz_be_0 + g_yyzzz_0_xxxxzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxxzzz_0[i] = 2.0 * g_yzzz_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxxzzz_1[i] * fz_be_0 + g_yyzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxyyyy_0[i] = 2.0 * g_yyyz_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxxyyyy_1[i] * fz_be_0 + g_yyyzz_0_xxxyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxxyyyz_0[i] = 2.0 * g_yzzz_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_xxxyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxyyzz_0[i] = 2.0 * g_yzzz_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xxxyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxyzzz_0[i] = 2.0 * g_yzzz_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxyzzz_1[i] * fz_be_0 + g_yyzzz_0_xxxzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxxzzzz_0[i] = 2.0 * g_yzzz_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxxzzzz_1[i] * fz_be_0 + g_yyzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyyyyy_0[i] = 2.0 * g_yyyz_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xxyyyyy_1[i] * fz_be_0 + g_yyyzz_0_xxyyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xxyyyyz_0[i] = 2.0 * g_yzzz_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyzzz_0_xxyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyyyzz_0[i] = 2.0 * g_yzzz_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_xxyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyyzzz_0[i] = 2.0 * g_yzzz_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xxyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxyzzzz_0[i] = 2.0 * g_yzzz_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxyzzzz_1[i] * fz_be_0 + g_yyzzz_0_xxzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xxzzzzz_0[i] = 2.0 * g_yzzz_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xxzzzzz_1[i] * fz_be_0 + g_yyzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyyyyy_0[i] = 2.0 * g_yyyz_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_xyyyyyy_1[i] * fz_be_0 + g_yyyzz_0_xyyyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_xyyyyyz_0[i] = 2.0 * g_yzzz_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyzzz_0_xyyyyz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyyyzz_0[i] = 2.0 * g_yzzz_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyzzz_0_xyyyzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyyzzz_0[i] = 2.0 * g_yzzz_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_xyyzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyyzzzz_0[i] = 2.0 * g_yzzz_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_xyzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xyzzzzz_0[i] = 2.0 * g_yzzz_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xyzzzzz_1[i] * fz_be_0 + g_yyzzz_0_xzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_xzzzzzz_0[i] = 2.0 * g_yzzz_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_xzzzzzz_1[i] * fz_be_0 + g_yyzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyyyyy_0[i] = 2.0 * g_yyyz_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_yyyz_0_yyyyyyy_1[i] * fz_be_0 + g_yyyzz_0_yyyyyyy_1[i] * wa_z[i];

        g_yyyzzz_0_yyyyyyz_0[i] = 2.0 * g_yzzz_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyzzz_0_yyyyyz_1[i] * fi_acd_0 + g_yyzzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyyyzz_0[i] = 2.0 * g_yzzz_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyzzz_0_yyyyzz_1[i] * fi_acd_0 + g_yyzzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyyzzz_0[i] = 2.0 * g_yzzz_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyzzz_0_yyyzzz_1[i] * fi_acd_0 + g_yyzzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyyzzzz_0[i] = 2.0 * g_yzzz_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyzzz_0_yyzzzz_1[i] * fi_acd_0 + g_yyzzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yyzzzzz_0[i] = 2.0 * g_yzzz_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzz_0_yzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_yzzzzzz_0[i] = 2.0 * g_yzzz_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_yzzzzzz_1[i] * fz_be_0 + g_yyzzz_0_zzzzzz_1[i] * fi_acd_0 + g_yyzzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yyyzzz_0_zzzzzzz_0[i] = 2.0 * g_yzzz_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzz_0_zzzzzzz_1[i] * fz_be_0 + g_yyzzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 900-936 components of targeted buffer : ISK

    auto g_yyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 900);

    auto g_yyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 901);

    auto g_yyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 902);

    auto g_yyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 903);

    auto g_yyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 904);

    auto g_yyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 905);

    auto g_yyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 906);

    auto g_yyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 907);

    auto g_yyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 908);

    auto g_yyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 909);

    auto g_yyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 910);

    auto g_yyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 911);

    auto g_yyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 912);

    auto g_yyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 913);

    auto g_yyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 914);

    auto g_yyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 915);

    auto g_yyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 916);

    auto g_yyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 917);

    auto g_yyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 918);

    auto g_yyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 919);

    auto g_yyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 920);

    auto g_yyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 921);

    auto g_yyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 922);

    auto g_yyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 923);

    auto g_yyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 924);

    auto g_yyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 925);

    auto g_yyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 926);

    auto g_yyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 927);

    auto g_yyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 928);

    auto g_yyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 929);

    auto g_yyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 930);

    auto g_yyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 931);

    auto g_yyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 932);

    auto g_yyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 933);

    auto g_yyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 934);

    auto g_yyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 935);

    #pragma omp simd aligned(g_yyzz_0_xxxxxxy_0, g_yyzz_0_xxxxxxy_1, g_yyzz_0_xxxxxyy_0, g_yyzz_0_xxxxxyy_1, g_yyzz_0_xxxxyyy_0, g_yyzz_0_xxxxyyy_1, g_yyzz_0_xxxyyyy_0, g_yyzz_0_xxxyyyy_1, g_yyzz_0_xxyyyyy_0, g_yyzz_0_xxyyyyy_1, g_yyzz_0_xyyyyyy_0, g_yyzz_0_xyyyyyy_1, g_yyzz_0_yyyyyyy_0, g_yyzz_0_yyyyyyy_1, g_yyzzz_0_xxxxxxy_1, g_yyzzz_0_xxxxxyy_1, g_yyzzz_0_xxxxyyy_1, g_yyzzz_0_xxxyyyy_1, g_yyzzz_0_xxyyyyy_1, g_yyzzz_0_xyyyyyy_1, g_yyzzz_0_yyyyyyy_1, g_yyzzzz_0_xxxxxxx_0, g_yyzzzz_0_xxxxxxy_0, g_yyzzzz_0_xxxxxxz_0, g_yyzzzz_0_xxxxxyy_0, g_yyzzzz_0_xxxxxyz_0, g_yyzzzz_0_xxxxxzz_0, g_yyzzzz_0_xxxxyyy_0, g_yyzzzz_0_xxxxyyz_0, g_yyzzzz_0_xxxxyzz_0, g_yyzzzz_0_xxxxzzz_0, g_yyzzzz_0_xxxyyyy_0, g_yyzzzz_0_xxxyyyz_0, g_yyzzzz_0_xxxyyzz_0, g_yyzzzz_0_xxxyzzz_0, g_yyzzzz_0_xxxzzzz_0, g_yyzzzz_0_xxyyyyy_0, g_yyzzzz_0_xxyyyyz_0, g_yyzzzz_0_xxyyyzz_0, g_yyzzzz_0_xxyyzzz_0, g_yyzzzz_0_xxyzzzz_0, g_yyzzzz_0_xxzzzzz_0, g_yyzzzz_0_xyyyyyy_0, g_yyzzzz_0_xyyyyyz_0, g_yyzzzz_0_xyyyyzz_0, g_yyzzzz_0_xyyyzzz_0, g_yyzzzz_0_xyyzzzz_0, g_yyzzzz_0_xyzzzzz_0, g_yyzzzz_0_xzzzzzz_0, g_yyzzzz_0_yyyyyyy_0, g_yyzzzz_0_yyyyyyz_0, g_yyzzzz_0_yyyyyzz_0, g_yyzzzz_0_yyyyzzz_0, g_yyzzzz_0_yyyzzzz_0, g_yyzzzz_0_yyzzzzz_0, g_yyzzzz_0_yzzzzzz_0, g_yyzzzz_0_zzzzzzz_0, g_yzzzz_0_xxxxxxx_1, g_yzzzz_0_xxxxxxz_1, g_yzzzz_0_xxxxxyz_1, g_yzzzz_0_xxxxxz_1, g_yzzzz_0_xxxxxzz_1, g_yzzzz_0_xxxxyyz_1, g_yzzzz_0_xxxxyz_1, g_yzzzz_0_xxxxyzz_1, g_yzzzz_0_xxxxzz_1, g_yzzzz_0_xxxxzzz_1, g_yzzzz_0_xxxyyyz_1, g_yzzzz_0_xxxyyz_1, g_yzzzz_0_xxxyyzz_1, g_yzzzz_0_xxxyzz_1, g_yzzzz_0_xxxyzzz_1, g_yzzzz_0_xxxzzz_1, g_yzzzz_0_xxxzzzz_1, g_yzzzz_0_xxyyyyz_1, g_yzzzz_0_xxyyyz_1, g_yzzzz_0_xxyyyzz_1, g_yzzzz_0_xxyyzz_1, g_yzzzz_0_xxyyzzz_1, g_yzzzz_0_xxyzzz_1, g_yzzzz_0_xxyzzzz_1, g_yzzzz_0_xxzzzz_1, g_yzzzz_0_xxzzzzz_1, g_yzzzz_0_xyyyyyz_1, g_yzzzz_0_xyyyyz_1, g_yzzzz_0_xyyyyzz_1, g_yzzzz_0_xyyyzz_1, g_yzzzz_0_xyyyzzz_1, g_yzzzz_0_xyyzzz_1, g_yzzzz_0_xyyzzzz_1, g_yzzzz_0_xyzzzz_1, g_yzzzz_0_xyzzzzz_1, g_yzzzz_0_xzzzzz_1, g_yzzzz_0_xzzzzzz_1, g_yzzzz_0_yyyyyyz_1, g_yzzzz_0_yyyyyz_1, g_yzzzz_0_yyyyyzz_1, g_yzzzz_0_yyyyzz_1, g_yzzzz_0_yyyyzzz_1, g_yzzzz_0_yyyzzz_1, g_yzzzz_0_yyyzzzz_1, g_yzzzz_0_yyzzzz_1, g_yzzzz_0_yyzzzzz_1, g_yzzzz_0_yzzzzz_1, g_yzzzz_0_yzzzzzz_1, g_yzzzz_0_zzzzzz_1, g_yzzzz_0_zzzzzzz_1, g_zzzz_0_xxxxxxx_0, g_zzzz_0_xxxxxxx_1, g_zzzz_0_xxxxxxz_0, g_zzzz_0_xxxxxxz_1, g_zzzz_0_xxxxxyz_0, g_zzzz_0_xxxxxyz_1, g_zzzz_0_xxxxxzz_0, g_zzzz_0_xxxxxzz_1, g_zzzz_0_xxxxyyz_0, g_zzzz_0_xxxxyyz_1, g_zzzz_0_xxxxyzz_0, g_zzzz_0_xxxxyzz_1, g_zzzz_0_xxxxzzz_0, g_zzzz_0_xxxxzzz_1, g_zzzz_0_xxxyyyz_0, g_zzzz_0_xxxyyyz_1, g_zzzz_0_xxxyyzz_0, g_zzzz_0_xxxyyzz_1, g_zzzz_0_xxxyzzz_0, g_zzzz_0_xxxyzzz_1, g_zzzz_0_xxxzzzz_0, g_zzzz_0_xxxzzzz_1, g_zzzz_0_xxyyyyz_0, g_zzzz_0_xxyyyyz_1, g_zzzz_0_xxyyyzz_0, g_zzzz_0_xxyyyzz_1, g_zzzz_0_xxyyzzz_0, g_zzzz_0_xxyyzzz_1, g_zzzz_0_xxyzzzz_0, g_zzzz_0_xxyzzzz_1, g_zzzz_0_xxzzzzz_0, g_zzzz_0_xxzzzzz_1, g_zzzz_0_xyyyyyz_0, g_zzzz_0_xyyyyyz_1, g_zzzz_0_xyyyyzz_0, g_zzzz_0_xyyyyzz_1, g_zzzz_0_xyyyzzz_0, g_zzzz_0_xyyyzzz_1, g_zzzz_0_xyyzzzz_0, g_zzzz_0_xyyzzzz_1, g_zzzz_0_xyzzzzz_0, g_zzzz_0_xyzzzzz_1, g_zzzz_0_xzzzzzz_0, g_zzzz_0_xzzzzzz_1, g_zzzz_0_yyyyyyz_0, g_zzzz_0_yyyyyyz_1, g_zzzz_0_yyyyyzz_0, g_zzzz_0_yyyyyzz_1, g_zzzz_0_yyyyzzz_0, g_zzzz_0_yyyyzzz_1, g_zzzz_0_yyyzzzz_0, g_zzzz_0_yyyzzzz_1, g_zzzz_0_yyzzzzz_0, g_zzzz_0_yyzzzzz_1, g_zzzz_0_yzzzzzz_0, g_zzzz_0_yzzzzzz_1, g_zzzz_0_zzzzzzz_0, g_zzzz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzz_0_xxxxxxx_0[i] = g_zzzz_0_xxxxxxx_0[i] * fbe_0 - g_zzzz_0_xxxxxxx_1[i] * fz_be_0 + g_yzzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxxxy_0[i] = 3.0 * g_yyzz_0_xxxxxxy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxxy_1[i] * fz_be_0 + g_yyzzz_0_xxxxxxy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxxxxz_0[i] = g_zzzz_0_xxxxxxz_0[i] * fbe_0 - g_zzzz_0_xxxxxxz_1[i] * fz_be_0 + g_yzzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxxyy_0[i] = 3.0 * g_yyzz_0_xxxxxyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxxyy_1[i] * fz_be_0 + g_yyzzz_0_xxxxxyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxxxyz_0[i] = g_zzzz_0_xxxxxyz_0[i] * fbe_0 - g_zzzz_0_xxxxxyz_1[i] * fz_be_0 + g_yzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxxzz_0[i] = g_zzzz_0_xxxxxzz_0[i] * fbe_0 - g_zzzz_0_xxxxxzz_1[i] * fz_be_0 + g_yzzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxyyy_0[i] = 3.0 * g_yyzz_0_xxxxyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxxyyy_1[i] * fz_be_0 + g_yyzzz_0_xxxxyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxxyyz_0[i] = g_zzzz_0_xxxxyyz_0[i] * fbe_0 - g_zzzz_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxyzz_0[i] = g_zzzz_0_xxxxyzz_0[i] * fbe_0 - g_zzzz_0_xxxxyzz_1[i] * fz_be_0 + g_yzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxxzzz_0[i] = g_zzzz_0_xxxxzzz_0[i] * fbe_0 - g_zzzz_0_xxxxzzz_1[i] * fz_be_0 + g_yzzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxyyyy_0[i] = 3.0 * g_yyzz_0_xxxyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxxyyyy_1[i] * fz_be_0 + g_yyzzz_0_xxxyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxxyyyz_0[i] = g_zzzz_0_xxxyyyz_0[i] * fbe_0 - g_zzzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxyyzz_0[i] = g_zzzz_0_xxxyyzz_0[i] * fbe_0 - g_zzzz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxyzzz_0[i] = g_zzzz_0_xxxyzzz_0[i] * fbe_0 - g_zzzz_0_xxxyzzz_1[i] * fz_be_0 + g_yzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxxzzzz_0[i] = g_zzzz_0_xxxzzzz_0[i] * fbe_0 - g_zzzz_0_xxxzzzz_1[i] * fz_be_0 + g_yzzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyyyyy_0[i] = 3.0 * g_yyzz_0_xxyyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xxyyyyy_1[i] * fz_be_0 + g_yyzzz_0_xxyyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xxyyyyz_0[i] = g_zzzz_0_xxyyyyz_0[i] * fbe_0 - g_zzzz_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyyyzz_0[i] = g_zzzz_0_xxyyyzz_0[i] * fbe_0 - g_zzzz_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyyzzz_0[i] = g_zzzz_0_xxyyzzz_0[i] * fbe_0 - g_zzzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxyzzzz_0[i] = g_zzzz_0_xxyzzzz_0[i] * fbe_0 - g_zzzz_0_xxyzzzz_1[i] * fz_be_0 + g_yzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xxzzzzz_0[i] = g_zzzz_0_xxzzzzz_0[i] * fbe_0 - g_zzzz_0_xxzzzzz_1[i] * fz_be_0 + g_yzzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyyyyy_0[i] = 3.0 * g_yyzz_0_xyyyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_xyyyyyy_1[i] * fz_be_0 + g_yyzzz_0_xyyyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_xyyyyyz_0[i] = g_zzzz_0_xyyyyyz_0[i] * fbe_0 - g_zzzz_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyyyzz_0[i] = g_zzzz_0_xyyyyzz_0[i] * fbe_0 - g_zzzz_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyyzzz_0[i] = g_zzzz_0_xyyyzzz_0[i] * fbe_0 - g_zzzz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyyzzzz_0[i] = g_zzzz_0_xyyzzzz_0[i] * fbe_0 - g_zzzz_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xyzzzzz_0[i] = g_zzzz_0_xyzzzzz_0[i] * fbe_0 - g_zzzz_0_xyzzzzz_1[i] * fz_be_0 + g_yzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_xzzzzzz_0[i] = g_zzzz_0_xzzzzzz_0[i] * fbe_0 - g_zzzz_0_xzzzzzz_1[i] * fz_be_0 + g_yzzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyyyyy_0[i] = 3.0 * g_yyzz_0_yyyyyyy_0[i] * fbe_0 - 3.0 * g_yyzz_0_yyyyyyy_1[i] * fz_be_0 + g_yyzzz_0_yyyyyyy_1[i] * wa_z[i];

        g_yyzzzz_0_yyyyyyz_0[i] = g_zzzz_0_yyyyyyz_0[i] * fbe_0 - g_zzzz_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_yzzzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyyyzz_0[i] = g_zzzz_0_yyyyyzz_0[i] * fbe_0 - g_zzzz_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_yzzzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyyzzz_0[i] = g_zzzz_0_yyyyzzz_0[i] * fbe_0 - g_zzzz_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_yzzzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyyzzzz_0[i] = g_zzzz_0_yyyzzzz_0[i] * fbe_0 - g_zzzz_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_yzzzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yyzzzzz_0[i] = g_zzzz_0_yyzzzzz_0[i] * fbe_0 - g_zzzz_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_yzzzzzz_0[i] = g_zzzz_0_yzzzzzz_0[i] * fbe_0 - g_zzzz_0_yzzzzzz_1[i] * fz_be_0 + g_yzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_yzzzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yyzzzz_0_zzzzzzz_0[i] = g_zzzz_0_zzzzzzz_0[i] * fbe_0 - g_zzzz_0_zzzzzzz_1[i] * fz_be_0 + g_yzzzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 936-972 components of targeted buffer : ISK

    auto g_yzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 936);

    auto g_yzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 937);

    auto g_yzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 938);

    auto g_yzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 939);

    auto g_yzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 940);

    auto g_yzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 941);

    auto g_yzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 942);

    auto g_yzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 943);

    auto g_yzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 944);

    auto g_yzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 945);

    auto g_yzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 946);

    auto g_yzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 947);

    auto g_yzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 948);

    auto g_yzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 949);

    auto g_yzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 950);

    auto g_yzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 951);

    auto g_yzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 952);

    auto g_yzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 953);

    auto g_yzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 954);

    auto g_yzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 955);

    auto g_yzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 956);

    auto g_yzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 957);

    auto g_yzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 958);

    auto g_yzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 959);

    auto g_yzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 960);

    auto g_yzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 961);

    auto g_yzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 962);

    auto g_yzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 963);

    auto g_yzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 964);

    auto g_yzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 965);

    auto g_yzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 966);

    auto g_yzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 967);

    auto g_yzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 968);

    auto g_yzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 969);

    auto g_yzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 970);

    auto g_yzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 971);

    #pragma omp simd aligned(g_yzzzzz_0_xxxxxxx_0, g_yzzzzz_0_xxxxxxy_0, g_yzzzzz_0_xxxxxxz_0, g_yzzzzz_0_xxxxxyy_0, g_yzzzzz_0_xxxxxyz_0, g_yzzzzz_0_xxxxxzz_0, g_yzzzzz_0_xxxxyyy_0, g_yzzzzz_0_xxxxyyz_0, g_yzzzzz_0_xxxxyzz_0, g_yzzzzz_0_xxxxzzz_0, g_yzzzzz_0_xxxyyyy_0, g_yzzzzz_0_xxxyyyz_0, g_yzzzzz_0_xxxyyzz_0, g_yzzzzz_0_xxxyzzz_0, g_yzzzzz_0_xxxzzzz_0, g_yzzzzz_0_xxyyyyy_0, g_yzzzzz_0_xxyyyyz_0, g_yzzzzz_0_xxyyyzz_0, g_yzzzzz_0_xxyyzzz_0, g_yzzzzz_0_xxyzzzz_0, g_yzzzzz_0_xxzzzzz_0, g_yzzzzz_0_xyyyyyy_0, g_yzzzzz_0_xyyyyyz_0, g_yzzzzz_0_xyyyyzz_0, g_yzzzzz_0_xyyyzzz_0, g_yzzzzz_0_xyyzzzz_0, g_yzzzzz_0_xyzzzzz_0, g_yzzzzz_0_xzzzzzz_0, g_yzzzzz_0_yyyyyyy_0, g_yzzzzz_0_yyyyyyz_0, g_yzzzzz_0_yyyyyzz_0, g_yzzzzz_0_yyyyzzz_0, g_yzzzzz_0_yyyzzzz_0, g_yzzzzz_0_yyzzzzz_0, g_yzzzzz_0_yzzzzzz_0, g_yzzzzz_0_zzzzzzz_0, g_zzzzz_0_xxxxxx_1, g_zzzzz_0_xxxxxxx_1, g_zzzzz_0_xxxxxxy_1, g_zzzzz_0_xxxxxxz_1, g_zzzzz_0_xxxxxy_1, g_zzzzz_0_xxxxxyy_1, g_zzzzz_0_xxxxxyz_1, g_zzzzz_0_xxxxxz_1, g_zzzzz_0_xxxxxzz_1, g_zzzzz_0_xxxxyy_1, g_zzzzz_0_xxxxyyy_1, g_zzzzz_0_xxxxyyz_1, g_zzzzz_0_xxxxyz_1, g_zzzzz_0_xxxxyzz_1, g_zzzzz_0_xxxxzz_1, g_zzzzz_0_xxxxzzz_1, g_zzzzz_0_xxxyyy_1, g_zzzzz_0_xxxyyyy_1, g_zzzzz_0_xxxyyyz_1, g_zzzzz_0_xxxyyz_1, g_zzzzz_0_xxxyyzz_1, g_zzzzz_0_xxxyzz_1, g_zzzzz_0_xxxyzzz_1, g_zzzzz_0_xxxzzz_1, g_zzzzz_0_xxxzzzz_1, g_zzzzz_0_xxyyyy_1, g_zzzzz_0_xxyyyyy_1, g_zzzzz_0_xxyyyyz_1, g_zzzzz_0_xxyyyz_1, g_zzzzz_0_xxyyyzz_1, g_zzzzz_0_xxyyzz_1, g_zzzzz_0_xxyyzzz_1, g_zzzzz_0_xxyzzz_1, g_zzzzz_0_xxyzzzz_1, g_zzzzz_0_xxzzzz_1, g_zzzzz_0_xxzzzzz_1, g_zzzzz_0_xyyyyy_1, g_zzzzz_0_xyyyyyy_1, g_zzzzz_0_xyyyyyz_1, g_zzzzz_0_xyyyyz_1, g_zzzzz_0_xyyyyzz_1, g_zzzzz_0_xyyyzz_1, g_zzzzz_0_xyyyzzz_1, g_zzzzz_0_xyyzzz_1, g_zzzzz_0_xyyzzzz_1, g_zzzzz_0_xyzzzz_1, g_zzzzz_0_xyzzzzz_1, g_zzzzz_0_xzzzzz_1, g_zzzzz_0_xzzzzzz_1, g_zzzzz_0_yyyyyy_1, g_zzzzz_0_yyyyyyy_1, g_zzzzz_0_yyyyyyz_1, g_zzzzz_0_yyyyyz_1, g_zzzzz_0_yyyyyzz_1, g_zzzzz_0_yyyyzz_1, g_zzzzz_0_yyyyzzz_1, g_zzzzz_0_yyyzzz_1, g_zzzzz_0_yyyzzzz_1, g_zzzzz_0_yyzzzz_1, g_zzzzz_0_yyzzzzz_1, g_zzzzz_0_yzzzzz_1, g_zzzzz_0_yzzzzzz_1, g_zzzzz_0_zzzzzz_1, g_zzzzz_0_zzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzz_0_xxxxxxx_0[i] = g_zzzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxxy_0[i] = g_zzzzz_0_xxxxxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxxz_0[i] = g_zzzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxyy_0[i] = 2.0 * g_zzzzz_0_xxxxxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxyz_0[i] = g_zzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxxzz_0[i] = g_zzzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxyyy_0[i] = 3.0 * g_zzzzz_0_xxxxyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxyyz_0[i] = 2.0 * g_zzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxyzz_0[i] = g_zzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxxzzz_0[i] = g_zzzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyyyy_0[i] = 4.0 * g_zzzzz_0_xxxyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyyyz_0[i] = 3.0 * g_zzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyyzz_0[i] = 2.0 * g_zzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxyzzz_0[i] = g_zzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxxzzzz_0[i] = g_zzzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyyyy_0[i] = 5.0 * g_zzzzz_0_xxyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyyyz_0[i] = 4.0 * g_zzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyyzz_0[i] = 3.0 * g_zzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyyzzz_0[i] = 2.0 * g_zzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxyzzzz_0[i] = g_zzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xxzzzzz_0[i] = g_zzzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyyyy_0[i] = 6.0 * g_zzzzz_0_xyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyyyz_0[i] = 5.0 * g_zzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyyzz_0[i] = 4.0 * g_zzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyyzzz_0[i] = 3.0 * g_zzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyyzzzz_0[i] = 2.0 * g_zzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xyzzzzz_0[i] = g_zzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_xzzzzzz_0[i] = g_zzzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyyyy_0[i] = 7.0 * g_zzzzz_0_yyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyyy_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyyyz_0[i] = 6.0 * g_zzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyyzz_0[i] = 5.0 * g_zzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyyzzz_0[i] = 4.0 * g_zzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyyzzzz_0[i] = 3.0 * g_zzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yyzzzzz_0[i] = 2.0 * g_zzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_yzzzzzz_0[i] = g_zzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yzzzzz_0_zzzzzzz_0[i] = g_zzzzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 972-1008 components of targeted buffer : ISK

    auto g_zzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_isk + 972);

    auto g_zzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_isk + 973);

    auto g_zzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_isk + 974);

    auto g_zzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_isk + 975);

    auto g_zzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_isk + 976);

    auto g_zzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_isk + 977);

    auto g_zzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_isk + 978);

    auto g_zzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_isk + 979);

    auto g_zzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_isk + 980);

    auto g_zzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_isk + 981);

    auto g_zzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_isk + 982);

    auto g_zzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_isk + 983);

    auto g_zzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_isk + 984);

    auto g_zzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_isk + 985);

    auto g_zzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_isk + 986);

    auto g_zzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_isk + 987);

    auto g_zzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_isk + 988);

    auto g_zzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_isk + 989);

    auto g_zzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_isk + 990);

    auto g_zzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_isk + 991);

    auto g_zzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_isk + 992);

    auto g_zzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 993);

    auto g_zzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 994);

    auto g_zzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 995);

    auto g_zzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 996);

    auto g_zzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 997);

    auto g_zzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 998);

    auto g_zzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 999);

    auto g_zzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_isk + 1000);

    auto g_zzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_isk + 1001);

    auto g_zzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_isk + 1002);

    auto g_zzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_isk + 1003);

    auto g_zzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_isk + 1004);

    auto g_zzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_isk + 1005);

    auto g_zzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 1006);

    auto g_zzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_isk + 1007);

    #pragma omp simd aligned(g_zzzz_0_xxxxxxx_0, g_zzzz_0_xxxxxxx_1, g_zzzz_0_xxxxxxy_0, g_zzzz_0_xxxxxxy_1, g_zzzz_0_xxxxxxz_0, g_zzzz_0_xxxxxxz_1, g_zzzz_0_xxxxxyy_0, g_zzzz_0_xxxxxyy_1, g_zzzz_0_xxxxxyz_0, g_zzzz_0_xxxxxyz_1, g_zzzz_0_xxxxxzz_0, g_zzzz_0_xxxxxzz_1, g_zzzz_0_xxxxyyy_0, g_zzzz_0_xxxxyyy_1, g_zzzz_0_xxxxyyz_0, g_zzzz_0_xxxxyyz_1, g_zzzz_0_xxxxyzz_0, g_zzzz_0_xxxxyzz_1, g_zzzz_0_xxxxzzz_0, g_zzzz_0_xxxxzzz_1, g_zzzz_0_xxxyyyy_0, g_zzzz_0_xxxyyyy_1, g_zzzz_0_xxxyyyz_0, g_zzzz_0_xxxyyyz_1, g_zzzz_0_xxxyyzz_0, g_zzzz_0_xxxyyzz_1, g_zzzz_0_xxxyzzz_0, g_zzzz_0_xxxyzzz_1, g_zzzz_0_xxxzzzz_0, g_zzzz_0_xxxzzzz_1, g_zzzz_0_xxyyyyy_0, g_zzzz_0_xxyyyyy_1, g_zzzz_0_xxyyyyz_0, g_zzzz_0_xxyyyyz_1, g_zzzz_0_xxyyyzz_0, g_zzzz_0_xxyyyzz_1, g_zzzz_0_xxyyzzz_0, g_zzzz_0_xxyyzzz_1, g_zzzz_0_xxyzzzz_0, g_zzzz_0_xxyzzzz_1, g_zzzz_0_xxzzzzz_0, g_zzzz_0_xxzzzzz_1, g_zzzz_0_xyyyyyy_0, g_zzzz_0_xyyyyyy_1, g_zzzz_0_xyyyyyz_0, g_zzzz_0_xyyyyyz_1, g_zzzz_0_xyyyyzz_0, g_zzzz_0_xyyyyzz_1, g_zzzz_0_xyyyzzz_0, g_zzzz_0_xyyyzzz_1, g_zzzz_0_xyyzzzz_0, g_zzzz_0_xyyzzzz_1, g_zzzz_0_xyzzzzz_0, g_zzzz_0_xyzzzzz_1, g_zzzz_0_xzzzzzz_0, g_zzzz_0_xzzzzzz_1, g_zzzz_0_yyyyyyy_0, g_zzzz_0_yyyyyyy_1, g_zzzz_0_yyyyyyz_0, g_zzzz_0_yyyyyyz_1, g_zzzz_0_yyyyyzz_0, g_zzzz_0_yyyyyzz_1, g_zzzz_0_yyyyzzz_0, g_zzzz_0_yyyyzzz_1, g_zzzz_0_yyyzzzz_0, g_zzzz_0_yyyzzzz_1, g_zzzz_0_yyzzzzz_0, g_zzzz_0_yyzzzzz_1, g_zzzz_0_yzzzzzz_0, g_zzzz_0_yzzzzzz_1, g_zzzz_0_zzzzzzz_0, g_zzzz_0_zzzzzzz_1, g_zzzzz_0_xxxxxx_1, g_zzzzz_0_xxxxxxx_1, g_zzzzz_0_xxxxxxy_1, g_zzzzz_0_xxxxxxz_1, g_zzzzz_0_xxxxxy_1, g_zzzzz_0_xxxxxyy_1, g_zzzzz_0_xxxxxyz_1, g_zzzzz_0_xxxxxz_1, g_zzzzz_0_xxxxxzz_1, g_zzzzz_0_xxxxyy_1, g_zzzzz_0_xxxxyyy_1, g_zzzzz_0_xxxxyyz_1, g_zzzzz_0_xxxxyz_1, g_zzzzz_0_xxxxyzz_1, g_zzzzz_0_xxxxzz_1, g_zzzzz_0_xxxxzzz_1, g_zzzzz_0_xxxyyy_1, g_zzzzz_0_xxxyyyy_1, g_zzzzz_0_xxxyyyz_1, g_zzzzz_0_xxxyyz_1, g_zzzzz_0_xxxyyzz_1, g_zzzzz_0_xxxyzz_1, g_zzzzz_0_xxxyzzz_1, g_zzzzz_0_xxxzzz_1, g_zzzzz_0_xxxzzzz_1, g_zzzzz_0_xxyyyy_1, g_zzzzz_0_xxyyyyy_1, g_zzzzz_0_xxyyyyz_1, g_zzzzz_0_xxyyyz_1, g_zzzzz_0_xxyyyzz_1, g_zzzzz_0_xxyyzz_1, g_zzzzz_0_xxyyzzz_1, g_zzzzz_0_xxyzzz_1, g_zzzzz_0_xxyzzzz_1, g_zzzzz_0_xxzzzz_1, g_zzzzz_0_xxzzzzz_1, g_zzzzz_0_xyyyyy_1, g_zzzzz_0_xyyyyyy_1, g_zzzzz_0_xyyyyyz_1, g_zzzzz_0_xyyyyz_1, g_zzzzz_0_xyyyyzz_1, g_zzzzz_0_xyyyzz_1, g_zzzzz_0_xyyyzzz_1, g_zzzzz_0_xyyzzz_1, g_zzzzz_0_xyyzzzz_1, g_zzzzz_0_xyzzzz_1, g_zzzzz_0_xyzzzzz_1, g_zzzzz_0_xzzzzz_1, g_zzzzz_0_xzzzzzz_1, g_zzzzz_0_yyyyyy_1, g_zzzzz_0_yyyyyyy_1, g_zzzzz_0_yyyyyyz_1, g_zzzzz_0_yyyyyz_1, g_zzzzz_0_yyyyyzz_1, g_zzzzz_0_yyyyzz_1, g_zzzzz_0_yyyyzzz_1, g_zzzzz_0_yyyzzz_1, g_zzzzz_0_yyyzzzz_1, g_zzzzz_0_yyzzzz_1, g_zzzzz_0_yyzzzzz_1, g_zzzzz_0_yzzzzz_1, g_zzzzz_0_yzzzzzz_1, g_zzzzz_0_zzzzzz_1, g_zzzzz_0_zzzzzzz_1, g_zzzzzz_0_xxxxxxx_0, g_zzzzzz_0_xxxxxxy_0, g_zzzzzz_0_xxxxxxz_0, g_zzzzzz_0_xxxxxyy_0, g_zzzzzz_0_xxxxxyz_0, g_zzzzzz_0_xxxxxzz_0, g_zzzzzz_0_xxxxyyy_0, g_zzzzzz_0_xxxxyyz_0, g_zzzzzz_0_xxxxyzz_0, g_zzzzzz_0_xxxxzzz_0, g_zzzzzz_0_xxxyyyy_0, g_zzzzzz_0_xxxyyyz_0, g_zzzzzz_0_xxxyyzz_0, g_zzzzzz_0_xxxyzzz_0, g_zzzzzz_0_xxxzzzz_0, g_zzzzzz_0_xxyyyyy_0, g_zzzzzz_0_xxyyyyz_0, g_zzzzzz_0_xxyyyzz_0, g_zzzzzz_0_xxyyzzz_0, g_zzzzzz_0_xxyzzzz_0, g_zzzzzz_0_xxzzzzz_0, g_zzzzzz_0_xyyyyyy_0, g_zzzzzz_0_xyyyyyz_0, g_zzzzzz_0_xyyyyzz_0, g_zzzzzz_0_xyyyzzz_0, g_zzzzzz_0_xyyzzzz_0, g_zzzzzz_0_xyzzzzz_0, g_zzzzzz_0_xzzzzzz_0, g_zzzzzz_0_yyyyyyy_0, g_zzzzzz_0_yyyyyyz_0, g_zzzzzz_0_yyyyyzz_0, g_zzzzzz_0_yyyyzzz_0, g_zzzzzz_0_yyyzzzz_0, g_zzzzzz_0_yyzzzzz_0, g_zzzzzz_0_yzzzzzz_0, g_zzzzzz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzz_0_xxxxxxx_0[i] = 5.0 * g_zzzz_0_xxxxxxx_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxxx_1[i] * fz_be_0 + g_zzzzz_0_xxxxxxx_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxxy_0[i] = 5.0 * g_zzzz_0_xxxxxxy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxxy_1[i] * fz_be_0 + g_zzzzz_0_xxxxxxy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxxz_0[i] = 5.0 * g_zzzz_0_xxxxxxz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxxz_1[i] * fz_be_0 + g_zzzzz_0_xxxxxx_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxxz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxyy_0[i] = 5.0 * g_zzzz_0_xxxxxyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxyy_1[i] * fz_be_0 + g_zzzzz_0_xxxxxyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxyz_0[i] = 5.0 * g_zzzz_0_xxxxxyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxyz_1[i] * fz_be_0 + g_zzzzz_0_xxxxxy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxxzz_0[i] = 5.0 * g_zzzz_0_xxxxxzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxxzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxxzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxyyy_0[i] = 5.0 * g_zzzz_0_xxxxyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxyyy_1[i] * fz_be_0 + g_zzzzz_0_xxxxyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxyyz_0[i] = 5.0 * g_zzzz_0_xxxxyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxyyz_1[i] * fz_be_0 + g_zzzzz_0_xxxxyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxyzz_0[i] = 5.0 * g_zzzz_0_xxxxyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxxzzz_0[i] = 5.0 * g_zzzz_0_xxxxzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxxzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxxzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyyyy_0[i] = 5.0 * g_zzzz_0_xxxyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyyyy_1[i] * fz_be_0 + g_zzzzz_0_xxxyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyyyz_0[i] = 5.0 * g_zzzz_0_xxxyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyyyz_1[i] * fz_be_0 + g_zzzzz_0_xxxyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyyzz_0[i] = 5.0 * g_zzzz_0_xxxyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxyzzz_0[i] = 5.0 * g_zzzz_0_xxxyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxxzzzz_0[i] = 5.0 * g_zzzz_0_xxxzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxxzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxxzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyyyy_0[i] = 5.0 * g_zzzz_0_xxyyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyyyy_1[i] * fz_be_0 + g_zzzzz_0_xxyyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyyyz_0[i] = 5.0 * g_zzzz_0_xxyyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyyyz_1[i] * fz_be_0 + g_zzzzz_0_xxyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyyzz_0[i] = 5.0 * g_zzzz_0_xxyyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyyzzz_0[i] = 5.0 * g_zzzz_0_xxyyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxyzzzz_0[i] = 5.0 * g_zzzz_0_xxyzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxyzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xxzzzzz_0[i] = 5.0 * g_zzzz_0_xxzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xxzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xxzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyyyy_0[i] = 5.0 * g_zzzz_0_xyyyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyyyy_1[i] * fz_be_0 + g_zzzzz_0_xyyyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyyyz_0[i] = 5.0 * g_zzzz_0_xyyyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyyyz_1[i] * fz_be_0 + g_zzzzz_0_xyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyyzz_0[i] = 5.0 * g_zzzz_0_xyyyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyyzzz_0[i] = 5.0 * g_zzzz_0_xyyyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyyzzzz_0[i] = 5.0 * g_zzzz_0_xyyzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyyzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xyzzzzz_0[i] = 5.0 * g_zzzz_0_xyzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xyzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_xzzzzzz_0[i] = 5.0 * g_zzzz_0_xzzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_xzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_xzzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyyyy_0[i] = 5.0 * g_zzzz_0_yyyyyyy_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyyyy_1[i] * fz_be_0 + g_zzzzz_0_yyyyyyy_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyyyz_0[i] = 5.0 * g_zzzz_0_yyyyyyz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyyyz_1[i] * fz_be_0 + g_zzzzz_0_yyyyyy_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyyz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyyzz_0[i] = 5.0 * g_zzzz_0_yyyyyzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyyzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyyzzz_0[i] = 5.0 * g_zzzz_0_yyyyzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyyzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyyzzzz_0[i] = 5.0 * g_zzzz_0_yyyzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyyzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yyzzzzz_0[i] = 5.0 * g_zzzz_0_yyzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yyzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_yzzzzzz_0[i] = 5.0 * g_zzzz_0_yzzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_yzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_yzzzzzz_1[i] * wa_z[i];

        g_zzzzzz_0_zzzzzzz_0[i] = 5.0 * g_zzzz_0_zzzzzzz_0[i] * fbe_0 - 5.0 * g_zzzz_0_zzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_zzzzz_0_zzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

