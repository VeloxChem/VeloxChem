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

#include "ThreeCenterElectronRepulsionPrimRecKSK.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ksk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ksk,
                                 size_t idx_eri_0_hsk,
                                 size_t idx_eri_1_hsk,
                                 size_t idx_eri_1_isi,
                                 size_t idx_eri_1_isk,
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

    /// Set up components of auxilary buffer : HSK

    auto g_xxxxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk);

    auto g_xxxxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 1);

    auto g_xxxxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 2);

    auto g_xxxxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 3);

    auto g_xxxxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 4);

    auto g_xxxxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 5);

    auto g_xxxxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 6);

    auto g_xxxxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 7);

    auto g_xxxxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 8);

    auto g_xxxxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 9);

    auto g_xxxxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 10);

    auto g_xxxxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 11);

    auto g_xxxxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 12);

    auto g_xxxxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 13);

    auto g_xxxxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 14);

    auto g_xxxxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 15);

    auto g_xxxxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 16);

    auto g_xxxxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 17);

    auto g_xxxxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 18);

    auto g_xxxxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 19);

    auto g_xxxxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 20);

    auto g_xxxxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 21);

    auto g_xxxxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 22);

    auto g_xxxxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 23);

    auto g_xxxxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 24);

    auto g_xxxxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 25);

    auto g_xxxxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 26);

    auto g_xxxxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 27);

    auto g_xxxxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 28);

    auto g_xxxxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 29);

    auto g_xxxxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 30);

    auto g_xxxxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 31);

    auto g_xxxxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 32);

    auto g_xxxxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 33);

    auto g_xxxxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 34);

    auto g_xxxxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 35);

    auto g_xxxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 36);

    auto g_xxxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 38);

    auto g_xxxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 41);

    auto g_xxxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 45);

    auto g_xxxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 50);

    auto g_xxxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 56);

    auto g_xxxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 63);

    auto g_xxxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 72);

    auto g_xxxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 73);

    auto g_xxxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 75);

    auto g_xxxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 78);

    auto g_xxxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 82);

    auto g_xxxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 87);

    auto g_xxxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 93);

    auto g_xxxyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 108);

    auto g_xxxyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 109);

    auto g_xxxyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 110);

    auto g_xxxyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 111);

    auto g_xxxyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 112);

    auto g_xxxyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 113);

    auto g_xxxyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 114);

    auto g_xxxyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 115);

    auto g_xxxyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 116);

    auto g_xxxyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 117);

    auto g_xxxyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 118);

    auto g_xxxyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 119);

    auto g_xxxyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 120);

    auto g_xxxyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 121);

    auto g_xxxyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 122);

    auto g_xxxyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 123);

    auto g_xxxyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 124);

    auto g_xxxyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 125);

    auto g_xxxyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 126);

    auto g_xxxyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 127);

    auto g_xxxyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 128);

    auto g_xxxyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 129);

    auto g_xxxyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 130);

    auto g_xxxyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 131);

    auto g_xxxyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 132);

    auto g_xxxyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 133);

    auto g_xxxyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 134);

    auto g_xxxyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 135);

    auto g_xxxyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 136);

    auto g_xxxyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 137);

    auto g_xxxyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 138);

    auto g_xxxyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 139);

    auto g_xxxyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 140);

    auto g_xxxyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 141);

    auto g_xxxyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 142);

    auto g_xxxyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 143);

    auto g_xxxzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 180);

    auto g_xxxzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 181);

    auto g_xxxzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 182);

    auto g_xxxzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 183);

    auto g_xxxzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 184);

    auto g_xxxzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 185);

    auto g_xxxzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 186);

    auto g_xxxzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 187);

    auto g_xxxzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 188);

    auto g_xxxzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 189);

    auto g_xxxzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 190);

    auto g_xxxzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 191);

    auto g_xxxzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 192);

    auto g_xxxzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 193);

    auto g_xxxzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 194);

    auto g_xxxzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 195);

    auto g_xxxzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 196);

    auto g_xxxzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 197);

    auto g_xxxzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 198);

    auto g_xxxzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 199);

    auto g_xxxzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 200);

    auto g_xxxzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 201);

    auto g_xxxzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 202);

    auto g_xxxzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 203);

    auto g_xxxzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 204);

    auto g_xxxzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 205);

    auto g_xxxzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 206);

    auto g_xxxzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 207);

    auto g_xxxzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 208);

    auto g_xxxzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 209);

    auto g_xxxzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 210);

    auto g_xxxzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 211);

    auto g_xxxzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 212);

    auto g_xxxzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 213);

    auto g_xxxzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 214);

    auto g_xxxzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 215);

    auto g_xxyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 216);

    auto g_xxyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 217);

    auto g_xxyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 218);

    auto g_xxyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 219);

    auto g_xxyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 220);

    auto g_xxyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 221);

    auto g_xxyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 222);

    auto g_xxyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 223);

    auto g_xxyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 224);

    auto g_xxyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 225);

    auto g_xxyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 226);

    auto g_xxyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 227);

    auto g_xxyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 228);

    auto g_xxyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 229);

    auto g_xxyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 230);

    auto g_xxyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 231);

    auto g_xxyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 232);

    auto g_xxyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 233);

    auto g_xxyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 234);

    auto g_xxyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 235);

    auto g_xxyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 236);

    auto g_xxyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 237);

    auto g_xxyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 238);

    auto g_xxyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 239);

    auto g_xxyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 240);

    auto g_xxyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 241);

    auto g_xxyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 242);

    auto g_xxyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 243);

    auto g_xxyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 244);

    auto g_xxyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 245);

    auto g_xxyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 246);

    auto g_xxyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 247);

    auto g_xxyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 248);

    auto g_xxyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 249);

    auto g_xxyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 250);

    auto g_xxyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 251);

    auto g_xxyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 253);

    auto g_xxyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 255);

    auto g_xxyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 258);

    auto g_xxyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 262);

    auto g_xxyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 267);

    auto g_xxyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 273);

    auto g_xxyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 288);

    auto g_xxyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 290);

    auto g_xxyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 293);

    auto g_xxyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 297);

    auto g_xxyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 302);

    auto g_xxyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 308);

    auto g_xxyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 315);

    auto g_xxzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 324);

    auto g_xxzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 325);

    auto g_xxzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 326);

    auto g_xxzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 327);

    auto g_xxzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 328);

    auto g_xxzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 329);

    auto g_xxzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 330);

    auto g_xxzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 331);

    auto g_xxzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 332);

    auto g_xxzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 333);

    auto g_xxzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 334);

    auto g_xxzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 335);

    auto g_xxzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 336);

    auto g_xxzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 337);

    auto g_xxzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 338);

    auto g_xxzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 339);

    auto g_xxzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 340);

    auto g_xxzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 341);

    auto g_xxzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 342);

    auto g_xxzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 343);

    auto g_xxzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 344);

    auto g_xxzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 345);

    auto g_xxzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 346);

    auto g_xxzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 347);

    auto g_xxzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 348);

    auto g_xxzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 349);

    auto g_xxzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 350);

    auto g_xxzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 351);

    auto g_xxzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 352);

    auto g_xxzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 353);

    auto g_xxzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 354);

    auto g_xxzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 355);

    auto g_xxzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 356);

    auto g_xxzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 357);

    auto g_xxzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 358);

    auto g_xxzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 359);

    auto g_xyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 361);

    auto g_xyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 363);

    auto g_xyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 364);

    auto g_xyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 366);

    auto g_xyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 367);

    auto g_xyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 368);

    auto g_xyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 370);

    auto g_xyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 371);

    auto g_xyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 372);

    auto g_xyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 373);

    auto g_xyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 375);

    auto g_xyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 376);

    auto g_xyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 377);

    auto g_xyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 378);

    auto g_xyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 379);

    auto g_xyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 381);

    auto g_xyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 382);

    auto g_xyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 383);

    auto g_xyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 384);

    auto g_xyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 385);

    auto g_xyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 386);

    auto g_xyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 388);

    auto g_xyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 389);

    auto g_xyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 390);

    auto g_xyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 391);

    auto g_xyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 392);

    auto g_xyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 393);

    auto g_xyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 394);

    auto g_xyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 395);

    auto g_xyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 436);

    auto g_xyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 439);

    auto g_xyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 440);

    auto g_xyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 443);

    auto g_xyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 444);

    auto g_xyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 445);

    auto g_xyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 448);

    auto g_xyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 449);

    auto g_xyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 450);

    auto g_xyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 451);

    auto g_xyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 454);

    auto g_xyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 455);

    auto g_xyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 456);

    auto g_xyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 457);

    auto g_xyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 458);

    auto g_xyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 460);

    auto g_xyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 461);

    auto g_xyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 462);

    auto g_xyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 463);

    auto g_xyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 464);

    auto g_xyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 465);

    auto g_xyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 466);

    auto g_xyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 467);

    auto g_xzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 506);

    auto g_xzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 508);

    auto g_xzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 509);

    auto g_xzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 511);

    auto g_xzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 512);

    auto g_xzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 513);

    auto g_xzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 515);

    auto g_xzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 516);

    auto g_xzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 517);

    auto g_xzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 518);

    auto g_xzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 520);

    auto g_xzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 521);

    auto g_xzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 522);

    auto g_xzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 523);

    auto g_xzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 524);

    auto g_xzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 526);

    auto g_xzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 527);

    auto g_xzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 528);

    auto g_xzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 529);

    auto g_xzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 530);

    auto g_xzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 531);

    auto g_xzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 532);

    auto g_xzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 533);

    auto g_xzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 534);

    auto g_xzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 535);

    auto g_xzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 536);

    auto g_xzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 537);

    auto g_xzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 538);

    auto g_xzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 539);

    auto g_yyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 540);

    auto g_yyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 541);

    auto g_yyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 542);

    auto g_yyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 543);

    auto g_yyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 544);

    auto g_yyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 545);

    auto g_yyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 546);

    auto g_yyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 547);

    auto g_yyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 548);

    auto g_yyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 549);

    auto g_yyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 550);

    auto g_yyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 551);

    auto g_yyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 552);

    auto g_yyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 553);

    auto g_yyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 554);

    auto g_yyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 555);

    auto g_yyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 556);

    auto g_yyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 557);

    auto g_yyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 558);

    auto g_yyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 559);

    auto g_yyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 560);

    auto g_yyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 561);

    auto g_yyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 562);

    auto g_yyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 563);

    auto g_yyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 564);

    auto g_yyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 565);

    auto g_yyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 566);

    auto g_yyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 567);

    auto g_yyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 568);

    auto g_yyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 569);

    auto g_yyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 570);

    auto g_yyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 571);

    auto g_yyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 572);

    auto g_yyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 573);

    auto g_yyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 574);

    auto g_yyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 575);

    auto g_yyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 577);

    auto g_yyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 579);

    auto g_yyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 582);

    auto g_yyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 586);

    auto g_yyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 591);

    auto g_yyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 597);

    auto g_yyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 604);

    auto g_yyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 612);

    auto g_yyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 613);

    auto g_yyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 614);

    auto g_yyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 615);

    auto g_yyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 616);

    auto g_yyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 617);

    auto g_yyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 618);

    auto g_yyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 619);

    auto g_yyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 620);

    auto g_yyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 621);

    auto g_yyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 622);

    auto g_yyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 623);

    auto g_yyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 624);

    auto g_yyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 625);

    auto g_yyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 626);

    auto g_yyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 627);

    auto g_yyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 628);

    auto g_yyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 629);

    auto g_yyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 630);

    auto g_yyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 631);

    auto g_yyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 632);

    auto g_yyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 633);

    auto g_yyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 634);

    auto g_yyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 635);

    auto g_yyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 636);

    auto g_yyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 637);

    auto g_yyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 638);

    auto g_yyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 639);

    auto g_yyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 640);

    auto g_yyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 641);

    auto g_yyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 642);

    auto g_yyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 643);

    auto g_yyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 644);

    auto g_yyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 645);

    auto g_yyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 646);

    auto g_yyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 647);

    auto g_yyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 648);

    auto g_yyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 649);

    auto g_yyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 650);

    auto g_yyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 651);

    auto g_yyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 652);

    auto g_yyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 653);

    auto g_yyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 654);

    auto g_yyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 655);

    auto g_yyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 656);

    auto g_yyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 657);

    auto g_yyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 658);

    auto g_yyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 659);

    auto g_yyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 660);

    auto g_yyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 661);

    auto g_yyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 662);

    auto g_yyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 663);

    auto g_yyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 664);

    auto g_yyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 665);

    auto g_yyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 666);

    auto g_yyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 667);

    auto g_yyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 668);

    auto g_yyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 669);

    auto g_yyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 670);

    auto g_yyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 671);

    auto g_yyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 672);

    auto g_yyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 673);

    auto g_yyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 674);

    auto g_yyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 675);

    auto g_yyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 676);

    auto g_yyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 677);

    auto g_yyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 678);

    auto g_yyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 679);

    auto g_yyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 680);

    auto g_yyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 681);

    auto g_yyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 682);

    auto g_yyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 683);

    auto g_yzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 684);

    auto g_yzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 686);

    auto g_yzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 688);

    auto g_yzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 689);

    auto g_yzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 691);

    auto g_yzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 692);

    auto g_yzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 693);

    auto g_yzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 695);

    auto g_yzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 696);

    auto g_yzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 697);

    auto g_yzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 698);

    auto g_yzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 700);

    auto g_yzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 701);

    auto g_yzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 702);

    auto g_yzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 703);

    auto g_yzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 704);

    auto g_yzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 706);

    auto g_yzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 707);

    auto g_yzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 708);

    auto g_yzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 709);

    auto g_yzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 710);

    auto g_yzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 711);

    auto g_yzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 713);

    auto g_yzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 714);

    auto g_yzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 715);

    auto g_yzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 716);

    auto g_yzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 717);

    auto g_yzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 718);

    auto g_yzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 719);

    auto g_zzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_hsk + 720);

    auto g_zzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_hsk + 721);

    auto g_zzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_hsk + 722);

    auto g_zzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_hsk + 723);

    auto g_zzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_hsk + 724);

    auto g_zzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_hsk + 725);

    auto g_zzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_hsk + 726);

    auto g_zzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_hsk + 727);

    auto g_zzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_hsk + 728);

    auto g_zzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_hsk + 729);

    auto g_zzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_hsk + 730);

    auto g_zzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_hsk + 731);

    auto g_zzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_hsk + 732);

    auto g_zzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_hsk + 733);

    auto g_zzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_hsk + 734);

    auto g_zzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 735);

    auto g_zzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 736);

    auto g_zzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 737);

    auto g_zzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 738);

    auto g_zzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 739);

    auto g_zzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 740);

    auto g_zzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 741);

    auto g_zzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 742);

    auto g_zzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 743);

    auto g_zzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 744);

    auto g_zzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 745);

    auto g_zzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 746);

    auto g_zzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 747);

    auto g_zzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_hsk + 748);

    auto g_zzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_hsk + 749);

    auto g_zzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_hsk + 750);

    auto g_zzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_hsk + 751);

    auto g_zzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_hsk + 752);

    auto g_zzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 753);

    auto g_zzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 754);

    auto g_zzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_hsk + 755);

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

    auto g_xxxxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 38);

    auto g_xxxxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 41);

    auto g_xxxxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 45);

    auto g_xxxxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 50);

    auto g_xxxxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 56);

    auto g_xxxxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 63);

    auto g_xxxxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_hsk + 72);

    auto g_xxxxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_hsk + 73);

    auto g_xxxxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 75);

    auto g_xxxxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 78);

    auto g_xxxxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 82);

    auto g_xxxxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 87);

    auto g_xxxxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 93);

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

    auto g_yyyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_hsk + 579);

    auto g_yyyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_hsk + 582);

    auto g_yyyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_hsk + 586);

    auto g_yyyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 591);

    auto g_yyyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 597);

    auto g_yyyyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_hsk + 604);

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

    auto g_yzzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_hsk + 686);

    auto g_yzzzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_hsk + 688);

    auto g_yzzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_hsk + 689);

    auto g_yzzzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_hsk + 691);

    auto g_yzzzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_hsk + 692);

    auto g_yzzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_hsk + 693);

    auto g_yzzzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_hsk + 695);

    auto g_yzzzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_hsk + 696);

    auto g_yzzzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_hsk + 697);

    auto g_yzzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_hsk + 698);

    auto g_yzzzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 700);

    auto g_yzzzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 701);

    auto g_yzzzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 702);

    auto g_yzzzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 703);

    auto g_yzzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 704);

    auto g_yzzzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_hsk + 706);

    auto g_yzzzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_hsk + 707);

    auto g_yzzzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_hsk + 708);

    auto g_yzzzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_hsk + 709);

    auto g_yzzzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 710);

    auto g_yzzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_hsk + 711);

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

    auto g_xxxxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 58);

    auto g_xxxxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 60);

    auto g_xxxxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 61);

    auto g_xxxxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 63);

    auto g_xxxxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 64);

    auto g_xxxxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 65);

    auto g_xxxxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 67);

    auto g_xxxxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 68);

    auto g_xxxxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 69);

    auto g_xxxxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 70);

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

    auto g_xxyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 340);

    auto g_xxyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 343);

    auto g_xxyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 344);

    auto g_xxyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 347);

    auto g_xxyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 348);

    auto g_xxyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 349);

    auto g_xxyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 352);

    auto g_xxyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 353);

    auto g_xxyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 354);

    auto g_xxyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 355);

    auto g_xxyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 358);

    auto g_xxyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 359);

    auto g_xxyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 360);

    auto g_xxyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 361);

    auto g_xxyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 362);

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

    auto g_xyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 498);

    auto g_xyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 499);

    auto g_xyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 500);

    auto g_xyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 501);

    auto g_xyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 502);

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

    auto g_xyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_isi + 526);

    auto g_xyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_isi + 527);

    auto g_xyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_isi + 528);

    auto g_xyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_isi + 529);

    auto g_xyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_isi + 530);

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

    auto g_yyyyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_isi + 618);

    auto g_yyyyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_isi + 620);

    auto g_yyyyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_isi + 621);

    auto g_yyyyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_isi + 623);

    auto g_yyyyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_isi + 624);

    auto g_yyyyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_isi + 625);

    auto g_yyyyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_isi + 627);

    auto g_yyyyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_isi + 628);

    auto g_yyyyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_isi + 629);

    auto g_yyyyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_isi + 630);

    auto g_yyyyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_isi + 632);

    auto g_yyyyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_isi + 633);

    auto g_yyyyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_isi + 634);

    auto g_yyyyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_isi + 635);

    auto g_yyyyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_isi + 636);

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

    auto g_xxxxxy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 36);

    auto g_xxxxxy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 37);

    auto g_xxxxxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 38);

    auto g_xxxxxy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 39);

    auto g_xxxxxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 41);

    auto g_xxxxxy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 42);

    auto g_xxxxxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 45);

    auto g_xxxxxy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 46);

    auto g_xxxxxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 50);

    auto g_xxxxxy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 51);

    auto g_xxxxxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 56);

    auto g_xxxxxy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 57);

    auto g_xxxxxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 63);

    auto g_xxxxxy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 64);

    auto g_xxxxxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 72);

    auto g_xxxxxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 73);

    auto g_xxxxxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 74);

    auto g_xxxxxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 75);

    auto g_xxxxxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 76);

    auto g_xxxxxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 77);

    auto g_xxxxxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 78);

    auto g_xxxxxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 79);

    auto g_xxxxxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 80);

    auto g_xxxxxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 81);

    auto g_xxxxxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 82);

    auto g_xxxxxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 83);

    auto g_xxxxxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 84);

    auto g_xxxxxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 85);

    auto g_xxxxxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 86);

    auto g_xxxxxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 87);

    auto g_xxxxxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 88);

    auto g_xxxxxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 89);

    auto g_xxxxxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 90);

    auto g_xxxxxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 91);

    auto g_xxxxxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 92);

    auto g_xxxxxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 93);

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

    auto g_xxxyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 253);

    auto g_xxxyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 255);

    auto g_xxxyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 258);

    auto g_xxxyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 262);

    auto g_xxxyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 267);

    auto g_xxxyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 273);

    auto g_xxxyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 288);

    auto g_xxxyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 290);

    auto g_xxxyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 293);

    auto g_xxxyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 297);

    auto g_xxxyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 302);

    auto g_xxxyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 308);

    auto g_xxxyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 315);

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

    auto g_xxyyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 397);

    auto g_xxyyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 399);

    auto g_xxyyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 402);

    auto g_xxyyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 406);

    auto g_xxyyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 411);

    auto g_xxyyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 417);

    auto g_xxyyzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 432);

    auto g_xxyyzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 433);

    auto g_xxyyzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 434);

    auto g_xxyyzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 435);

    auto g_xxyyzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 436);

    auto g_xxyyzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 437);

    auto g_xxyyzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 438);

    auto g_xxyyzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 439);

    auto g_xxyyzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 440);

    auto g_xxyyzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 441);

    auto g_xxyyzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 442);

    auto g_xxyyzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 443);

    auto g_xxyyzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 444);

    auto g_xxyyzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 445);

    auto g_xxyyzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 446);

    auto g_xxyyzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 447);

    auto g_xxyyzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 448);

    auto g_xxyyzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 449);

    auto g_xxyyzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 450);

    auto g_xxyyzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 451);

    auto g_xxyyzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 452);

    auto g_xxyyzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 453);

    auto g_xxyyzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 454);

    auto g_xxyyzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 455);

    auto g_xxyyzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 456);

    auto g_xxyyzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 457);

    auto g_xxyyzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 458);

    auto g_xxyyzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 459);

    auto g_xxyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 460);

    auto g_xxyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 461);

    auto g_xxyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 462);

    auto g_xxyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 463);

    auto g_xxyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 464);

    auto g_xxyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 465);

    auto g_xxyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 466);

    auto g_xxyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 467);

    auto g_xxyzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 468);

    auto g_xxyzzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 470);

    auto g_xxyzzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 473);

    auto g_xxyzzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 477);

    auto g_xxyzzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 482);

    auto g_xxyzzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 488);

    auto g_xxyzzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 495);

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

    auto g_xyyyyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 540);

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

    auto g_xyyyyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 575);

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

    auto g_xyyyzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 640);

    auto g_xyyyzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 641);

    auto g_xyyyzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 642);

    auto g_xyyyzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 643);

    auto g_xyyyzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 644);

    auto g_xyyyzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 645);

    auto g_xyyyzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 646);

    auto g_xyyyzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 647);

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

    auto g_xyyzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 676);

    auto g_xyyzzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 677);

    auto g_xyyzzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 678);

    auto g_xyyzzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 679);

    auto g_xyyzzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 680);

    auto g_xyyzzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 681);

    auto g_xyyzzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 682);

    auto g_xyyzzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 683);

    auto g_xzzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 720);

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

    auto g_xzzzzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 748);

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

    auto g_yyyyyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_isk + 793);

    auto g_yyyyyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_isk + 794);

    auto g_yyyyyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_isk + 795);

    auto g_yyyyyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_isk + 796);

    auto g_yyyyyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_isk + 797);

    auto g_yyyyyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_isk + 798);

    auto g_yyyyyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_isk + 799);

    auto g_yyyyyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_isk + 800);

    auto g_yyyyyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_isk + 801);

    auto g_yyyyyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_isk + 802);

    auto g_yyyyyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_isk + 803);

    auto g_yyyyyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_isk + 804);

    auto g_yyyyyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_isk + 805);

    auto g_yyyyyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_isk + 806);

    auto g_yyyyyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_isk + 807);

    auto g_yyyyyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_isk + 808);

    auto g_yyyyyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_isk + 809);

    auto g_yyyyyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_isk + 810);

    auto g_yyyyyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_isk + 811);

    auto g_yyyyyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_isk + 812);

    auto g_yyyyyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 813);

    auto g_yyyyyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_isk + 814);

    auto g_yyyyyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_isk + 815);

    auto g_yyyyyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_isk + 816);

    auto g_yyyyyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_isk + 817);

    auto g_yyyyyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_isk + 818);

    auto g_yyyyyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_isk + 819);

    auto g_yyyyyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_isk + 820);

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

    auto g_yzzzzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_isk + 936);

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

    /// Set up 0-36 components of targeted buffer : KSK

    auto g_xxxxxxx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk);

    auto g_xxxxxxx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 1);

    auto g_xxxxxxx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 2);

    auto g_xxxxxxx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 3);

    auto g_xxxxxxx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 4);

    auto g_xxxxxxx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 5);

    auto g_xxxxxxx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 6);

    auto g_xxxxxxx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 7);

    auto g_xxxxxxx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 8);

    auto g_xxxxxxx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 9);

    auto g_xxxxxxx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 10);

    auto g_xxxxxxx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 11);

    auto g_xxxxxxx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 12);

    auto g_xxxxxxx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 13);

    auto g_xxxxxxx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 14);

    auto g_xxxxxxx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 15);

    auto g_xxxxxxx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 16);

    auto g_xxxxxxx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 17);

    auto g_xxxxxxx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 18);

    auto g_xxxxxxx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 19);

    auto g_xxxxxxx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 20);

    auto g_xxxxxxx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 21);

    auto g_xxxxxxx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 22);

    auto g_xxxxxxx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 23);

    auto g_xxxxxxx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 24);

    auto g_xxxxxxx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 25);

    auto g_xxxxxxx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 26);

    auto g_xxxxxxx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 27);

    auto g_xxxxxxx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 28);

    auto g_xxxxxxx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 29);

    auto g_xxxxxxx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 30);

    auto g_xxxxxxx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 31);

    auto g_xxxxxxx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 32);

    auto g_xxxxxxx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 33);

    auto g_xxxxxxx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 34);

    auto g_xxxxxxx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 35);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxxx_0, g_xxxxx_0_xxxxxxx_1, g_xxxxx_0_xxxxxxy_0, g_xxxxx_0_xxxxxxy_1, g_xxxxx_0_xxxxxxz_0, g_xxxxx_0_xxxxxxz_1, g_xxxxx_0_xxxxxyy_0, g_xxxxx_0_xxxxxyy_1, g_xxxxx_0_xxxxxyz_0, g_xxxxx_0_xxxxxyz_1, g_xxxxx_0_xxxxxzz_0, g_xxxxx_0_xxxxxzz_1, g_xxxxx_0_xxxxyyy_0, g_xxxxx_0_xxxxyyy_1, g_xxxxx_0_xxxxyyz_0, g_xxxxx_0_xxxxyyz_1, g_xxxxx_0_xxxxyzz_0, g_xxxxx_0_xxxxyzz_1, g_xxxxx_0_xxxxzzz_0, g_xxxxx_0_xxxxzzz_1, g_xxxxx_0_xxxyyyy_0, g_xxxxx_0_xxxyyyy_1, g_xxxxx_0_xxxyyyz_0, g_xxxxx_0_xxxyyyz_1, g_xxxxx_0_xxxyyzz_0, g_xxxxx_0_xxxyyzz_1, g_xxxxx_0_xxxyzzz_0, g_xxxxx_0_xxxyzzz_1, g_xxxxx_0_xxxzzzz_0, g_xxxxx_0_xxxzzzz_1, g_xxxxx_0_xxyyyyy_0, g_xxxxx_0_xxyyyyy_1, g_xxxxx_0_xxyyyyz_0, g_xxxxx_0_xxyyyyz_1, g_xxxxx_0_xxyyyzz_0, g_xxxxx_0_xxyyyzz_1, g_xxxxx_0_xxyyzzz_0, g_xxxxx_0_xxyyzzz_1, g_xxxxx_0_xxyzzzz_0, g_xxxxx_0_xxyzzzz_1, g_xxxxx_0_xxzzzzz_0, g_xxxxx_0_xxzzzzz_1, g_xxxxx_0_xyyyyyy_0, g_xxxxx_0_xyyyyyy_1, g_xxxxx_0_xyyyyyz_0, g_xxxxx_0_xyyyyyz_1, g_xxxxx_0_xyyyyzz_0, g_xxxxx_0_xyyyyzz_1, g_xxxxx_0_xyyyzzz_0, g_xxxxx_0_xyyyzzz_1, g_xxxxx_0_xyyzzzz_0, g_xxxxx_0_xyyzzzz_1, g_xxxxx_0_xyzzzzz_0, g_xxxxx_0_xyzzzzz_1, g_xxxxx_0_xzzzzzz_0, g_xxxxx_0_xzzzzzz_1, g_xxxxx_0_yyyyyyy_0, g_xxxxx_0_yyyyyyy_1, g_xxxxx_0_yyyyyyz_0, g_xxxxx_0_yyyyyyz_1, g_xxxxx_0_yyyyyzz_0, g_xxxxx_0_yyyyyzz_1, g_xxxxx_0_yyyyzzz_0, g_xxxxx_0_yyyyzzz_1, g_xxxxx_0_yyyzzzz_0, g_xxxxx_0_yyyzzzz_1, g_xxxxx_0_yyzzzzz_0, g_xxxxx_0_yyzzzzz_1, g_xxxxx_0_yzzzzzz_0, g_xxxxx_0_yzzzzzz_1, g_xxxxx_0_zzzzzzz_0, g_xxxxx_0_zzzzzzz_1, g_xxxxxx_0_xxxxxx_1, g_xxxxxx_0_xxxxxxx_1, g_xxxxxx_0_xxxxxxy_1, g_xxxxxx_0_xxxxxxz_1, g_xxxxxx_0_xxxxxy_1, g_xxxxxx_0_xxxxxyy_1, g_xxxxxx_0_xxxxxyz_1, g_xxxxxx_0_xxxxxz_1, g_xxxxxx_0_xxxxxzz_1, g_xxxxxx_0_xxxxyy_1, g_xxxxxx_0_xxxxyyy_1, g_xxxxxx_0_xxxxyyz_1, g_xxxxxx_0_xxxxyz_1, g_xxxxxx_0_xxxxyzz_1, g_xxxxxx_0_xxxxzz_1, g_xxxxxx_0_xxxxzzz_1, g_xxxxxx_0_xxxyyy_1, g_xxxxxx_0_xxxyyyy_1, g_xxxxxx_0_xxxyyyz_1, g_xxxxxx_0_xxxyyz_1, g_xxxxxx_0_xxxyyzz_1, g_xxxxxx_0_xxxyzz_1, g_xxxxxx_0_xxxyzzz_1, g_xxxxxx_0_xxxzzz_1, g_xxxxxx_0_xxxzzzz_1, g_xxxxxx_0_xxyyyy_1, g_xxxxxx_0_xxyyyyy_1, g_xxxxxx_0_xxyyyyz_1, g_xxxxxx_0_xxyyyz_1, g_xxxxxx_0_xxyyyzz_1, g_xxxxxx_0_xxyyzz_1, g_xxxxxx_0_xxyyzzz_1, g_xxxxxx_0_xxyzzz_1, g_xxxxxx_0_xxyzzzz_1, g_xxxxxx_0_xxzzzz_1, g_xxxxxx_0_xxzzzzz_1, g_xxxxxx_0_xyyyyy_1, g_xxxxxx_0_xyyyyyy_1, g_xxxxxx_0_xyyyyyz_1, g_xxxxxx_0_xyyyyz_1, g_xxxxxx_0_xyyyyzz_1, g_xxxxxx_0_xyyyzz_1, g_xxxxxx_0_xyyyzzz_1, g_xxxxxx_0_xyyzzz_1, g_xxxxxx_0_xyyzzzz_1, g_xxxxxx_0_xyzzzz_1, g_xxxxxx_0_xyzzzzz_1, g_xxxxxx_0_xzzzzz_1, g_xxxxxx_0_xzzzzzz_1, g_xxxxxx_0_yyyyyy_1, g_xxxxxx_0_yyyyyyy_1, g_xxxxxx_0_yyyyyyz_1, g_xxxxxx_0_yyyyyz_1, g_xxxxxx_0_yyyyyzz_1, g_xxxxxx_0_yyyyzz_1, g_xxxxxx_0_yyyyzzz_1, g_xxxxxx_0_yyyzzz_1, g_xxxxxx_0_yyyzzzz_1, g_xxxxxx_0_yyzzzz_1, g_xxxxxx_0_yyzzzzz_1, g_xxxxxx_0_yzzzzz_1, g_xxxxxx_0_yzzzzzz_1, g_xxxxxx_0_zzzzzz_1, g_xxxxxx_0_zzzzzzz_1, g_xxxxxxx_0_xxxxxxx_0, g_xxxxxxx_0_xxxxxxy_0, g_xxxxxxx_0_xxxxxxz_0, g_xxxxxxx_0_xxxxxyy_0, g_xxxxxxx_0_xxxxxyz_0, g_xxxxxxx_0_xxxxxzz_0, g_xxxxxxx_0_xxxxyyy_0, g_xxxxxxx_0_xxxxyyz_0, g_xxxxxxx_0_xxxxyzz_0, g_xxxxxxx_0_xxxxzzz_0, g_xxxxxxx_0_xxxyyyy_0, g_xxxxxxx_0_xxxyyyz_0, g_xxxxxxx_0_xxxyyzz_0, g_xxxxxxx_0_xxxyzzz_0, g_xxxxxxx_0_xxxzzzz_0, g_xxxxxxx_0_xxyyyyy_0, g_xxxxxxx_0_xxyyyyz_0, g_xxxxxxx_0_xxyyyzz_0, g_xxxxxxx_0_xxyyzzz_0, g_xxxxxxx_0_xxyzzzz_0, g_xxxxxxx_0_xxzzzzz_0, g_xxxxxxx_0_xyyyyyy_0, g_xxxxxxx_0_xyyyyyz_0, g_xxxxxxx_0_xyyyyzz_0, g_xxxxxxx_0_xyyyzzz_0, g_xxxxxxx_0_xyyzzzz_0, g_xxxxxxx_0_xyzzzzz_0, g_xxxxxxx_0_xzzzzzz_0, g_xxxxxxx_0_yyyyyyy_0, g_xxxxxxx_0_yyyyyyz_0, g_xxxxxxx_0_yyyyyzz_0, g_xxxxxxx_0_yyyyzzz_0, g_xxxxxxx_0_yyyzzzz_0, g_xxxxxxx_0_yyzzzzz_0, g_xxxxxxx_0_yzzzzzz_0, g_xxxxxxx_0_zzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxxx_0_xxxxxxx_0[i] = 6.0 * g_xxxxx_0_xxxxxxx_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxxx_1[i] * fz_be_0 + 7.0 * g_xxxxxx_0_xxxxxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxx_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxxy_0[i] = 6.0 * g_xxxxx_0_xxxxxxy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xxxxxx_0_xxxxxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxxz_0[i] = 6.0 * g_xxxxx_0_xxxxxxz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xxxxxx_0_xxxxxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxyy_0[i] = 6.0 * g_xxxxx_0_xxxxxyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xxxxxx_0_xxxxyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxyz_0[i] = 6.0 * g_xxxxx_0_xxxxxyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxxxxx_0_xxxxyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxxzz_0[i] = 6.0 * g_xxxxx_0_xxxxxzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xxxxxx_0_xxxxzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxyyy_0[i] = 6.0 * g_xxxxx_0_xxxxyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxyyz_0[i] = 6.0 * g_xxxxx_0_xxxxyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxyzz_0[i] = 6.0 * g_xxxxx_0_xxxxyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxxzzz_0[i] = 6.0 * g_xxxxx_0_xxxxzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xxxxxx_0_xxxzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyyyy_0[i] = 6.0 * g_xxxxx_0_xxxyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyyyz_0[i] = 6.0 * g_xxxxx_0_xxxyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyyzz_0[i] = 6.0 * g_xxxxx_0_xxxyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxyzzz_0[i] = 6.0 * g_xxxxx_0_xxxyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxxzzzz_0[i] = 6.0 * g_xxxxx_0_xxxzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xxxxxx_0_xxzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyyyy_0[i] = 6.0 * g_xxxxx_0_xxyyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyyyz_0[i] = 6.0 * g_xxxxx_0_xxyyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyyzz_0[i] = 6.0 * g_xxxxx_0_xxyyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyyzzz_0[i] = 6.0 * g_xxxxx_0_xxyyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxyzzzz_0[i] = 6.0 * g_xxxxx_0_xxyzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xxzzzzz_0[i] = 6.0 * g_xxxxx_0_xxzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_0_xzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyyyy_0[i] = 6.0 * g_xxxxx_0_xyyyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyyyy_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyyyz_0[i] = 6.0 * g_xxxxx_0_xyyyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyyyz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyyzz_0[i] = 6.0 * g_xxxxx_0_xyyyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyyzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyyzzz_0[i] = 6.0 * g_xxxxx_0_xyyyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyyzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyyzzzz_0[i] = 6.0 * g_xxxxx_0_xyyzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyyzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xyzzzzz_0[i] = 6.0 * g_xxxxx_0_xyzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xyzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_xzzzzzz_0[i] = 6.0 * g_xxxxx_0_xzzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_xzzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_zzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xzzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyyyy_0[i] = 6.0 * g_xxxxx_0_yyyyyyy_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyyyy_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyyyz_0[i] = 6.0 * g_xxxxx_0_yyyyyyz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyyyz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyyzz_0[i] = 6.0 * g_xxxxx_0_yyyyyzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyyzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyyzzz_0[i] = 6.0 * g_xxxxx_0_yyyyzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyyzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyyzzzz_0[i] = 6.0 * g_xxxxx_0_yyyzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyyzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yyzzzzz_0[i] = 6.0 * g_xxxxx_0_yyzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yyzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_yzzzzzz_0[i] = 6.0 * g_xxxxx_0_yzzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_yzzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxxxxx_0_zzzzzzz_0[i] = 6.0 * g_xxxxx_0_zzzzzzz_0[i] * fbe_0 - 6.0 * g_xxxxx_0_zzzzzzz_1[i] * fz_be_0 + g_xxxxxx_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 36-72 components of targeted buffer : KSK

    auto g_xxxxxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 36);

    auto g_xxxxxxy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 37);

    auto g_xxxxxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 38);

    auto g_xxxxxxy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 39);

    auto g_xxxxxxy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 40);

    auto g_xxxxxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 41);

    auto g_xxxxxxy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 42);

    auto g_xxxxxxy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 43);

    auto g_xxxxxxy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 44);

    auto g_xxxxxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 45);

    auto g_xxxxxxy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 46);

    auto g_xxxxxxy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 47);

    auto g_xxxxxxy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 48);

    auto g_xxxxxxy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 49);

    auto g_xxxxxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 50);

    auto g_xxxxxxy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 51);

    auto g_xxxxxxy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 52);

    auto g_xxxxxxy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 53);

    auto g_xxxxxxy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 54);

    auto g_xxxxxxy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 55);

    auto g_xxxxxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 56);

    auto g_xxxxxxy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 57);

    auto g_xxxxxxy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 58);

    auto g_xxxxxxy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 59);

    auto g_xxxxxxy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 60);

    auto g_xxxxxxy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 61);

    auto g_xxxxxxy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 62);

    auto g_xxxxxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 63);

    auto g_xxxxxxy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 64);

    auto g_xxxxxxy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 65);

    auto g_xxxxxxy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 66);

    auto g_xxxxxxy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 67);

    auto g_xxxxxxy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 68);

    auto g_xxxxxxy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 69);

    auto g_xxxxxxy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 70);

    auto g_xxxxxxy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 71);

    #pragma omp simd aligned(g_xxxxxx_0_xxxxxx_1, g_xxxxxx_0_xxxxxxx_1, g_xxxxxx_0_xxxxxxy_1, g_xxxxxx_0_xxxxxxz_1, g_xxxxxx_0_xxxxxy_1, g_xxxxxx_0_xxxxxyy_1, g_xxxxxx_0_xxxxxyz_1, g_xxxxxx_0_xxxxxz_1, g_xxxxxx_0_xxxxxzz_1, g_xxxxxx_0_xxxxyy_1, g_xxxxxx_0_xxxxyyy_1, g_xxxxxx_0_xxxxyyz_1, g_xxxxxx_0_xxxxyz_1, g_xxxxxx_0_xxxxyzz_1, g_xxxxxx_0_xxxxzz_1, g_xxxxxx_0_xxxxzzz_1, g_xxxxxx_0_xxxyyy_1, g_xxxxxx_0_xxxyyyy_1, g_xxxxxx_0_xxxyyyz_1, g_xxxxxx_0_xxxyyz_1, g_xxxxxx_0_xxxyyzz_1, g_xxxxxx_0_xxxyzz_1, g_xxxxxx_0_xxxyzzz_1, g_xxxxxx_0_xxxzzz_1, g_xxxxxx_0_xxxzzzz_1, g_xxxxxx_0_xxyyyy_1, g_xxxxxx_0_xxyyyyy_1, g_xxxxxx_0_xxyyyyz_1, g_xxxxxx_0_xxyyyz_1, g_xxxxxx_0_xxyyyzz_1, g_xxxxxx_0_xxyyzz_1, g_xxxxxx_0_xxyyzzz_1, g_xxxxxx_0_xxyzzz_1, g_xxxxxx_0_xxyzzzz_1, g_xxxxxx_0_xxzzzz_1, g_xxxxxx_0_xxzzzzz_1, g_xxxxxx_0_xyyyyy_1, g_xxxxxx_0_xyyyyyy_1, g_xxxxxx_0_xyyyyyz_1, g_xxxxxx_0_xyyyyz_1, g_xxxxxx_0_xyyyyzz_1, g_xxxxxx_0_xyyyzz_1, g_xxxxxx_0_xyyyzzz_1, g_xxxxxx_0_xyyzzz_1, g_xxxxxx_0_xyyzzzz_1, g_xxxxxx_0_xyzzzz_1, g_xxxxxx_0_xyzzzzz_1, g_xxxxxx_0_xzzzzz_1, g_xxxxxx_0_xzzzzzz_1, g_xxxxxx_0_yyyyyy_1, g_xxxxxx_0_yyyyyyy_1, g_xxxxxx_0_yyyyyyz_1, g_xxxxxx_0_yyyyyz_1, g_xxxxxx_0_yyyyyzz_1, g_xxxxxx_0_yyyyzz_1, g_xxxxxx_0_yyyyzzz_1, g_xxxxxx_0_yyyzzz_1, g_xxxxxx_0_yyyzzzz_1, g_xxxxxx_0_yyzzzz_1, g_xxxxxx_0_yyzzzzz_1, g_xxxxxx_0_yzzzzz_1, g_xxxxxx_0_yzzzzzz_1, g_xxxxxx_0_zzzzzz_1, g_xxxxxx_0_zzzzzzz_1, g_xxxxxxy_0_xxxxxxx_0, g_xxxxxxy_0_xxxxxxy_0, g_xxxxxxy_0_xxxxxxz_0, g_xxxxxxy_0_xxxxxyy_0, g_xxxxxxy_0_xxxxxyz_0, g_xxxxxxy_0_xxxxxzz_0, g_xxxxxxy_0_xxxxyyy_0, g_xxxxxxy_0_xxxxyyz_0, g_xxxxxxy_0_xxxxyzz_0, g_xxxxxxy_0_xxxxzzz_0, g_xxxxxxy_0_xxxyyyy_0, g_xxxxxxy_0_xxxyyyz_0, g_xxxxxxy_0_xxxyyzz_0, g_xxxxxxy_0_xxxyzzz_0, g_xxxxxxy_0_xxxzzzz_0, g_xxxxxxy_0_xxyyyyy_0, g_xxxxxxy_0_xxyyyyz_0, g_xxxxxxy_0_xxyyyzz_0, g_xxxxxxy_0_xxyyzzz_0, g_xxxxxxy_0_xxyzzzz_0, g_xxxxxxy_0_xxzzzzz_0, g_xxxxxxy_0_xyyyyyy_0, g_xxxxxxy_0_xyyyyyz_0, g_xxxxxxy_0_xyyyyzz_0, g_xxxxxxy_0_xyyyzzz_0, g_xxxxxxy_0_xyyzzzz_0, g_xxxxxxy_0_xyzzzzz_0, g_xxxxxxy_0_xzzzzzz_0, g_xxxxxxy_0_yyyyyyy_0, g_xxxxxxy_0_yyyyyyz_0, g_xxxxxxy_0_yyyyyzz_0, g_xxxxxxy_0_yyyyzzz_0, g_xxxxxxy_0_yyyzzzz_0, g_xxxxxxy_0_yyzzzzz_0, g_xxxxxxy_0_yzzzzzz_0, g_xxxxxxy_0_zzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxy_0_xxxxxxx_0[i] = g_xxxxxx_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxxy_0[i] = g_xxxxxx_0_xxxxxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxxz_0[i] = g_xxxxxx_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxyy_0[i] = 2.0 * g_xxxxxx_0_xxxxxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxyz_0[i] = g_xxxxxx_0_xxxxxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxxzz_0[i] = g_xxxxxx_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxyyy_0[i] = 3.0 * g_xxxxxx_0_xxxxyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxyyz_0[i] = 2.0 * g_xxxxxx_0_xxxxyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxyzz_0[i] = g_xxxxxx_0_xxxxzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxxzzz_0[i] = g_xxxxxx_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyyyy_0[i] = 4.0 * g_xxxxxx_0_xxxyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyyyz_0[i] = 3.0 * g_xxxxxx_0_xxxyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyyzz_0[i] = 2.0 * g_xxxxxx_0_xxxyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxyzzz_0[i] = g_xxxxxx_0_xxxzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxxzzzz_0[i] = g_xxxxxx_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyyyy_0[i] = 5.0 * g_xxxxxx_0_xxyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyyyz_0[i] = 4.0 * g_xxxxxx_0_xxyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyyzz_0[i] = 3.0 * g_xxxxxx_0_xxyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyyzzz_0[i] = 2.0 * g_xxxxxx_0_xxyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxyzzzz_0[i] = g_xxxxxx_0_xxzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xxzzzzz_0[i] = g_xxxxxx_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyyyy_0[i] = 6.0 * g_xxxxxx_0_xyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyyyz_0[i] = 5.0 * g_xxxxxx_0_xyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyyzz_0[i] = 4.0 * g_xxxxxx_0_xyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyyzzz_0[i] = 3.0 * g_xxxxxx_0_xyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyyzzzz_0[i] = 2.0 * g_xxxxxx_0_xyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xyzzzzz_0[i] = g_xxxxxx_0_xzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_xzzzzzz_0[i] = g_xxxxxx_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyyyy_0[i] = 7.0 * g_xxxxxx_0_yyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyyy_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyyyz_0[i] = 6.0 * g_xxxxxx_0_yyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyyz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyyzz_0[i] = 5.0 * g_xxxxxx_0_yyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyyzzz_0[i] = 4.0 * g_xxxxxx_0_yyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyyzzzz_0[i] = 3.0 * g_xxxxxx_0_yyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yyzzzzz_0[i] = 2.0 * g_xxxxxx_0_yzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_yzzzzzz_0[i] = g_xxxxxx_0_zzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yzzzzzz_1[i] * wa_y[i];

        g_xxxxxxy_0_zzzzzzz_0[i] = g_xxxxxx_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 72-108 components of targeted buffer : KSK

    auto g_xxxxxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 72);

    auto g_xxxxxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 73);

    auto g_xxxxxxz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 74);

    auto g_xxxxxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 75);

    auto g_xxxxxxz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 76);

    auto g_xxxxxxz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 77);

    auto g_xxxxxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 78);

    auto g_xxxxxxz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 79);

    auto g_xxxxxxz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 80);

    auto g_xxxxxxz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 81);

    auto g_xxxxxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 82);

    auto g_xxxxxxz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 83);

    auto g_xxxxxxz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 84);

    auto g_xxxxxxz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 85);

    auto g_xxxxxxz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 86);

    auto g_xxxxxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 87);

    auto g_xxxxxxz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 88);

    auto g_xxxxxxz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 89);

    auto g_xxxxxxz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 90);

    auto g_xxxxxxz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 91);

    auto g_xxxxxxz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 92);

    auto g_xxxxxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 93);

    auto g_xxxxxxz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 94);

    auto g_xxxxxxz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 95);

    auto g_xxxxxxz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 96);

    auto g_xxxxxxz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 97);

    auto g_xxxxxxz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 98);

    auto g_xxxxxxz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 99);

    auto g_xxxxxxz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 100);

    auto g_xxxxxxz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 101);

    auto g_xxxxxxz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 102);

    auto g_xxxxxxz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 103);

    auto g_xxxxxxz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 104);

    auto g_xxxxxxz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 105);

    auto g_xxxxxxz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 106);

    auto g_xxxxxxz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 107);

    #pragma omp simd aligned(g_xxxxxx_0_xxxxxx_1, g_xxxxxx_0_xxxxxxx_1, g_xxxxxx_0_xxxxxxy_1, g_xxxxxx_0_xxxxxxz_1, g_xxxxxx_0_xxxxxy_1, g_xxxxxx_0_xxxxxyy_1, g_xxxxxx_0_xxxxxyz_1, g_xxxxxx_0_xxxxxz_1, g_xxxxxx_0_xxxxxzz_1, g_xxxxxx_0_xxxxyy_1, g_xxxxxx_0_xxxxyyy_1, g_xxxxxx_0_xxxxyyz_1, g_xxxxxx_0_xxxxyz_1, g_xxxxxx_0_xxxxyzz_1, g_xxxxxx_0_xxxxzz_1, g_xxxxxx_0_xxxxzzz_1, g_xxxxxx_0_xxxyyy_1, g_xxxxxx_0_xxxyyyy_1, g_xxxxxx_0_xxxyyyz_1, g_xxxxxx_0_xxxyyz_1, g_xxxxxx_0_xxxyyzz_1, g_xxxxxx_0_xxxyzz_1, g_xxxxxx_0_xxxyzzz_1, g_xxxxxx_0_xxxzzz_1, g_xxxxxx_0_xxxzzzz_1, g_xxxxxx_0_xxyyyy_1, g_xxxxxx_0_xxyyyyy_1, g_xxxxxx_0_xxyyyyz_1, g_xxxxxx_0_xxyyyz_1, g_xxxxxx_0_xxyyyzz_1, g_xxxxxx_0_xxyyzz_1, g_xxxxxx_0_xxyyzzz_1, g_xxxxxx_0_xxyzzz_1, g_xxxxxx_0_xxyzzzz_1, g_xxxxxx_0_xxzzzz_1, g_xxxxxx_0_xxzzzzz_1, g_xxxxxx_0_xyyyyy_1, g_xxxxxx_0_xyyyyyy_1, g_xxxxxx_0_xyyyyyz_1, g_xxxxxx_0_xyyyyz_1, g_xxxxxx_0_xyyyyzz_1, g_xxxxxx_0_xyyyzz_1, g_xxxxxx_0_xyyyzzz_1, g_xxxxxx_0_xyyzzz_1, g_xxxxxx_0_xyyzzzz_1, g_xxxxxx_0_xyzzzz_1, g_xxxxxx_0_xyzzzzz_1, g_xxxxxx_0_xzzzzz_1, g_xxxxxx_0_xzzzzzz_1, g_xxxxxx_0_yyyyyy_1, g_xxxxxx_0_yyyyyyy_1, g_xxxxxx_0_yyyyyyz_1, g_xxxxxx_0_yyyyyz_1, g_xxxxxx_0_yyyyyzz_1, g_xxxxxx_0_yyyyzz_1, g_xxxxxx_0_yyyyzzz_1, g_xxxxxx_0_yyyzzz_1, g_xxxxxx_0_yyyzzzz_1, g_xxxxxx_0_yyzzzz_1, g_xxxxxx_0_yyzzzzz_1, g_xxxxxx_0_yzzzzz_1, g_xxxxxx_0_yzzzzzz_1, g_xxxxxx_0_zzzzzz_1, g_xxxxxx_0_zzzzzzz_1, g_xxxxxxz_0_xxxxxxx_0, g_xxxxxxz_0_xxxxxxy_0, g_xxxxxxz_0_xxxxxxz_0, g_xxxxxxz_0_xxxxxyy_0, g_xxxxxxz_0_xxxxxyz_0, g_xxxxxxz_0_xxxxxzz_0, g_xxxxxxz_0_xxxxyyy_0, g_xxxxxxz_0_xxxxyyz_0, g_xxxxxxz_0_xxxxyzz_0, g_xxxxxxz_0_xxxxzzz_0, g_xxxxxxz_0_xxxyyyy_0, g_xxxxxxz_0_xxxyyyz_0, g_xxxxxxz_0_xxxyyzz_0, g_xxxxxxz_0_xxxyzzz_0, g_xxxxxxz_0_xxxzzzz_0, g_xxxxxxz_0_xxyyyyy_0, g_xxxxxxz_0_xxyyyyz_0, g_xxxxxxz_0_xxyyyzz_0, g_xxxxxxz_0_xxyyzzz_0, g_xxxxxxz_0_xxyzzzz_0, g_xxxxxxz_0_xxzzzzz_0, g_xxxxxxz_0_xyyyyyy_0, g_xxxxxxz_0_xyyyyyz_0, g_xxxxxxz_0_xyyyyzz_0, g_xxxxxxz_0_xyyyzzz_0, g_xxxxxxz_0_xyyzzzz_0, g_xxxxxxz_0_xyzzzzz_0, g_xxxxxxz_0_xzzzzzz_0, g_xxxxxxz_0_yyyyyyy_0, g_xxxxxxz_0_yyyyyyz_0, g_xxxxxxz_0_yyyyyzz_0, g_xxxxxxz_0_yyyyzzz_0, g_xxxxxxz_0_yyyzzzz_0, g_xxxxxxz_0_yyzzzzz_0, g_xxxxxxz_0_yzzzzzz_0, g_xxxxxxz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxz_0_xxxxxxx_0[i] = g_xxxxxx_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxxy_0[i] = g_xxxxxx_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxxz_0[i] = g_xxxxxx_0_xxxxxx_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxxz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxyy_0[i] = g_xxxxxx_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxyz_0[i] = g_xxxxxx_0_xxxxxy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxxzz_0[i] = 2.0 * g_xxxxxx_0_xxxxxz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxxzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxyyy_0[i] = g_xxxxxx_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxyyz_0[i] = g_xxxxxx_0_xxxxyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxyzz_0[i] = 2.0 * g_xxxxxx_0_xxxxyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxxzzz_0[i] = 3.0 * g_xxxxxx_0_xxxxzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxxzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyyyy_0[i] = g_xxxxxx_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyyyz_0[i] = g_xxxxxx_0_xxxyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyyzz_0[i] = 2.0 * g_xxxxxx_0_xxxyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxyzzz_0[i] = 3.0 * g_xxxxxx_0_xxxyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxxzzzz_0[i] = 4.0 * g_xxxxxx_0_xxxzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxxzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyyyy_0[i] = g_xxxxxx_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyyyz_0[i] = g_xxxxxx_0_xxyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyyzz_0[i] = 2.0 * g_xxxxxx_0_xxyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyyzzz_0[i] = 3.0 * g_xxxxxx_0_xxyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxyzzzz_0[i] = 4.0 * g_xxxxxx_0_xxyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxyzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xxzzzzz_0[i] = 5.0 * g_xxxxxx_0_xxzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xxzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyyyy_0[i] = g_xxxxxx_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyyyz_0[i] = g_xxxxxx_0_xyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyyzz_0[i] = 2.0 * g_xxxxxx_0_xyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyyzzz_0[i] = 3.0 * g_xxxxxx_0_xyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyyzzzz_0[i] = 4.0 * g_xxxxxx_0_xyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyyzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xyzzzzz_0[i] = 5.0 * g_xxxxxx_0_xyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xyzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_xzzzzzz_0[i] = 6.0 * g_xxxxxx_0_xzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_xzzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyyyy_0[i] = g_xxxxxx_0_yyyyyyy_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyyyz_0[i] = g_xxxxxx_0_yyyyyy_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyyz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyyzz_0[i] = 2.0 * g_xxxxxx_0_yyyyyz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyyzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyyzzz_0[i] = 3.0 * g_xxxxxx_0_yyyyzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyyzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyyzzzz_0[i] = 4.0 * g_xxxxxx_0_yyyzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyyzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yyzzzzz_0[i] = 5.0 * g_xxxxxx_0_yyzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yyzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_yzzzzzz_0[i] = 6.0 * g_xxxxxx_0_yzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_yzzzzzz_1[i] * wa_z[i];

        g_xxxxxxz_0_zzzzzzz_0[i] = 7.0 * g_xxxxxx_0_zzzzzz_1[i] * fi_acd_0 + g_xxxxxx_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 108-144 components of targeted buffer : KSK

    auto g_xxxxxyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 108);

    auto g_xxxxxyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 109);

    auto g_xxxxxyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 110);

    auto g_xxxxxyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 111);

    auto g_xxxxxyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 112);

    auto g_xxxxxyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 113);

    auto g_xxxxxyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 114);

    auto g_xxxxxyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 115);

    auto g_xxxxxyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 116);

    auto g_xxxxxyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 117);

    auto g_xxxxxyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 118);

    auto g_xxxxxyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 119);

    auto g_xxxxxyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 120);

    auto g_xxxxxyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 121);

    auto g_xxxxxyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 122);

    auto g_xxxxxyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 123);

    auto g_xxxxxyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 124);

    auto g_xxxxxyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 125);

    auto g_xxxxxyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 126);

    auto g_xxxxxyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 127);

    auto g_xxxxxyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 128);

    auto g_xxxxxyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 129);

    auto g_xxxxxyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 130);

    auto g_xxxxxyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 131);

    auto g_xxxxxyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 132);

    auto g_xxxxxyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 133);

    auto g_xxxxxyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 134);

    auto g_xxxxxyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 135);

    auto g_xxxxxyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 136);

    auto g_xxxxxyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 137);

    auto g_xxxxxyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 138);

    auto g_xxxxxyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 139);

    auto g_xxxxxyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 140);

    auto g_xxxxxyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 141);

    auto g_xxxxxyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 142);

    auto g_xxxxxyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 143);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxxx_0, g_xxxxx_0_xxxxxxx_1, g_xxxxx_0_xxxxxxz_0, g_xxxxx_0_xxxxxxz_1, g_xxxxx_0_xxxxxzz_0, g_xxxxx_0_xxxxxzz_1, g_xxxxx_0_xxxxzzz_0, g_xxxxx_0_xxxxzzz_1, g_xxxxx_0_xxxzzzz_0, g_xxxxx_0_xxxzzzz_1, g_xxxxx_0_xxzzzzz_0, g_xxxxx_0_xxzzzzz_1, g_xxxxx_0_xzzzzzz_0, g_xxxxx_0_xzzzzzz_1, g_xxxxxy_0_xxxxxxx_1, g_xxxxxy_0_xxxxxxz_1, g_xxxxxy_0_xxxxxzz_1, g_xxxxxy_0_xxxxzzz_1, g_xxxxxy_0_xxxzzzz_1, g_xxxxxy_0_xxzzzzz_1, g_xxxxxy_0_xzzzzzz_1, g_xxxxxyy_0_xxxxxxx_0, g_xxxxxyy_0_xxxxxxy_0, g_xxxxxyy_0_xxxxxxz_0, g_xxxxxyy_0_xxxxxyy_0, g_xxxxxyy_0_xxxxxyz_0, g_xxxxxyy_0_xxxxxzz_0, g_xxxxxyy_0_xxxxyyy_0, g_xxxxxyy_0_xxxxyyz_0, g_xxxxxyy_0_xxxxyzz_0, g_xxxxxyy_0_xxxxzzz_0, g_xxxxxyy_0_xxxyyyy_0, g_xxxxxyy_0_xxxyyyz_0, g_xxxxxyy_0_xxxyyzz_0, g_xxxxxyy_0_xxxyzzz_0, g_xxxxxyy_0_xxxzzzz_0, g_xxxxxyy_0_xxyyyyy_0, g_xxxxxyy_0_xxyyyyz_0, g_xxxxxyy_0_xxyyyzz_0, g_xxxxxyy_0_xxyyzzz_0, g_xxxxxyy_0_xxyzzzz_0, g_xxxxxyy_0_xxzzzzz_0, g_xxxxxyy_0_xyyyyyy_0, g_xxxxxyy_0_xyyyyyz_0, g_xxxxxyy_0_xyyyyzz_0, g_xxxxxyy_0_xyyyzzz_0, g_xxxxxyy_0_xyyzzzz_0, g_xxxxxyy_0_xyzzzzz_0, g_xxxxxyy_0_xzzzzzz_0, g_xxxxxyy_0_yyyyyyy_0, g_xxxxxyy_0_yyyyyyz_0, g_xxxxxyy_0_yyyyyzz_0, g_xxxxxyy_0_yyyyzzz_0, g_xxxxxyy_0_yyyzzzz_0, g_xxxxxyy_0_yyzzzzz_0, g_xxxxxyy_0_yzzzzzz_0, g_xxxxxyy_0_zzzzzzz_0, g_xxxxyy_0_xxxxxxy_1, g_xxxxyy_0_xxxxxy_1, g_xxxxyy_0_xxxxxyy_1, g_xxxxyy_0_xxxxxyz_1, g_xxxxyy_0_xxxxyy_1, g_xxxxyy_0_xxxxyyy_1, g_xxxxyy_0_xxxxyyz_1, g_xxxxyy_0_xxxxyz_1, g_xxxxyy_0_xxxxyzz_1, g_xxxxyy_0_xxxyyy_1, g_xxxxyy_0_xxxyyyy_1, g_xxxxyy_0_xxxyyyz_1, g_xxxxyy_0_xxxyyz_1, g_xxxxyy_0_xxxyyzz_1, g_xxxxyy_0_xxxyzz_1, g_xxxxyy_0_xxxyzzz_1, g_xxxxyy_0_xxyyyy_1, g_xxxxyy_0_xxyyyyy_1, g_xxxxyy_0_xxyyyyz_1, g_xxxxyy_0_xxyyyz_1, g_xxxxyy_0_xxyyyzz_1, g_xxxxyy_0_xxyyzz_1, g_xxxxyy_0_xxyyzzz_1, g_xxxxyy_0_xxyzzz_1, g_xxxxyy_0_xxyzzzz_1, g_xxxxyy_0_xyyyyy_1, g_xxxxyy_0_xyyyyyy_1, g_xxxxyy_0_xyyyyyz_1, g_xxxxyy_0_xyyyyz_1, g_xxxxyy_0_xyyyyzz_1, g_xxxxyy_0_xyyyzz_1, g_xxxxyy_0_xyyyzzz_1, g_xxxxyy_0_xyyzzz_1, g_xxxxyy_0_xyyzzzz_1, g_xxxxyy_0_xyzzzz_1, g_xxxxyy_0_xyzzzzz_1, g_xxxxyy_0_yyyyyy_1, g_xxxxyy_0_yyyyyyy_1, g_xxxxyy_0_yyyyyyz_1, g_xxxxyy_0_yyyyyz_1, g_xxxxyy_0_yyyyyzz_1, g_xxxxyy_0_yyyyzz_1, g_xxxxyy_0_yyyyzzz_1, g_xxxxyy_0_yyyzzz_1, g_xxxxyy_0_yyyzzzz_1, g_xxxxyy_0_yyzzzz_1, g_xxxxyy_0_yyzzzzz_1, g_xxxxyy_0_yzzzzz_1, g_xxxxyy_0_yzzzzzz_1, g_xxxxyy_0_zzzzzzz_1, g_xxxyy_0_xxxxxxy_0, g_xxxyy_0_xxxxxxy_1, g_xxxyy_0_xxxxxyy_0, g_xxxyy_0_xxxxxyy_1, g_xxxyy_0_xxxxxyz_0, g_xxxyy_0_xxxxxyz_1, g_xxxyy_0_xxxxyyy_0, g_xxxyy_0_xxxxyyy_1, g_xxxyy_0_xxxxyyz_0, g_xxxyy_0_xxxxyyz_1, g_xxxyy_0_xxxxyzz_0, g_xxxyy_0_xxxxyzz_1, g_xxxyy_0_xxxyyyy_0, g_xxxyy_0_xxxyyyy_1, g_xxxyy_0_xxxyyyz_0, g_xxxyy_0_xxxyyyz_1, g_xxxyy_0_xxxyyzz_0, g_xxxyy_0_xxxyyzz_1, g_xxxyy_0_xxxyzzz_0, g_xxxyy_0_xxxyzzz_1, g_xxxyy_0_xxyyyyy_0, g_xxxyy_0_xxyyyyy_1, g_xxxyy_0_xxyyyyz_0, g_xxxyy_0_xxyyyyz_1, g_xxxyy_0_xxyyyzz_0, g_xxxyy_0_xxyyyzz_1, g_xxxyy_0_xxyyzzz_0, g_xxxyy_0_xxyyzzz_1, g_xxxyy_0_xxyzzzz_0, g_xxxyy_0_xxyzzzz_1, g_xxxyy_0_xyyyyyy_0, g_xxxyy_0_xyyyyyy_1, g_xxxyy_0_xyyyyyz_0, g_xxxyy_0_xyyyyyz_1, g_xxxyy_0_xyyyyzz_0, g_xxxyy_0_xyyyyzz_1, g_xxxyy_0_xyyyzzz_0, g_xxxyy_0_xyyyzzz_1, g_xxxyy_0_xyyzzzz_0, g_xxxyy_0_xyyzzzz_1, g_xxxyy_0_xyzzzzz_0, g_xxxyy_0_xyzzzzz_1, g_xxxyy_0_yyyyyyy_0, g_xxxyy_0_yyyyyyy_1, g_xxxyy_0_yyyyyyz_0, g_xxxyy_0_yyyyyyz_1, g_xxxyy_0_yyyyyzz_0, g_xxxyy_0_yyyyyzz_1, g_xxxyy_0_yyyyzzz_0, g_xxxyy_0_yyyyzzz_1, g_xxxyy_0_yyyzzzz_0, g_xxxyy_0_yyyzzzz_1, g_xxxyy_0_yyzzzzz_0, g_xxxyy_0_yyzzzzz_1, g_xxxyy_0_yzzzzzz_0, g_xxxyy_0_yzzzzzz_1, g_xxxyy_0_zzzzzzz_0, g_xxxyy_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxyy_0_xxxxxxx_0[i] = g_xxxxx_0_xxxxxxx_0[i] * fbe_0 - g_xxxxx_0_xxxxxxx_1[i] * fz_be_0 + g_xxxxxy_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxxxxy_0[i] = 4.0 * g_xxxyy_0_xxxxxxy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xxxxyy_0_xxxxxy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxxxz_0[i] = g_xxxxx_0_xxxxxxz_0[i] * fbe_0 - g_xxxxx_0_xxxxxxz_1[i] * fz_be_0 + g_xxxxxy_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxxxyy_0[i] = 4.0 * g_xxxyy_0_xxxxxyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xxxxyy_0_xxxxyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxxyz_0[i] = 4.0 * g_xxxyy_0_xxxxxyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxxxyy_0_xxxxyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxxzz_0[i] = g_xxxxx_0_xxxxxzz_0[i] * fbe_0 - g_xxxxx_0_xxxxxzz_1[i] * fz_be_0 + g_xxxxxy_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxxyyy_0[i] = 4.0 * g_xxxyy_0_xxxxyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xxxxyy_0_xxxyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxyyz_0[i] = 4.0 * g_xxxyy_0_xxxxyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxxxyy_0_xxxyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxyzz_0[i] = 4.0 * g_xxxyy_0_xxxxyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxxxyy_0_xxxyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxxzzz_0[i] = g_xxxxx_0_xxxxzzz_0[i] * fbe_0 - g_xxxxx_0_xxxxzzz_1[i] * fz_be_0 + g_xxxxxy_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxxyyyy_0[i] = 4.0 * g_xxxyy_0_xxxyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxyyyz_0[i] = 4.0 * g_xxxyy_0_xxxyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxyyzz_0[i] = 4.0 * g_xxxyy_0_xxxyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxyzzz_0[i] = 4.0 * g_xxxyy_0_xxxyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxxxyy_0_xxyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxxzzzz_0[i] = g_xxxxx_0_xxxzzzz_0[i] * fbe_0 - g_xxxxx_0_xxxzzzz_1[i] * fz_be_0 + g_xxxxxy_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xxyyyyy_0[i] = 4.0 * g_xxxyy_0_xxyyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyyyyz_0[i] = 4.0 * g_xxxyy_0_xxyyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyyyzz_0[i] = 4.0 * g_xxxyy_0_xxyyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyyzzz_0[i] = 4.0 * g_xxxyy_0_xxyyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxyzzzz_0[i] = 4.0 * g_xxxyy_0_xxyzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxyy_0_xyzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xxzzzzz_0[i] = g_xxxxx_0_xxzzzzz_0[i] * fbe_0 - g_xxxxx_0_xxzzzzz_1[i] * fz_be_0 + g_xxxxxy_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_xyyyyyy_0[i] = 4.0 * g_xxxyy_0_xyyyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyyyy_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyyyyz_0[i] = 4.0 * g_xxxyy_0_xyyyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyyyz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyyyzz_0[i] = 4.0 * g_xxxyy_0_xyyyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyyzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyyzzz_0[i] = 4.0 * g_xxxyy_0_xyyyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyyzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyyzzzz_0[i] = 4.0 * g_xxxyy_0_xyyzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyyzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xyzzzzz_0[i] = 4.0 * g_xxxyy_0_xyzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_xyzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_xzzzzzz_0[i] = g_xxxxx_0_xzzzzzz_0[i] * fbe_0 - g_xxxxx_0_xzzzzzz_1[i] * fz_be_0 + g_xxxxxy_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxxxyy_0_yyyyyyy_0[i] = 4.0 * g_xxxyy_0_yyyyyyy_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyyyy_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyyyyz_0[i] = 4.0 * g_xxxyy_0_yyyyyyz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyyyz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyyyzz_0[i] = 4.0 * g_xxxyy_0_yyyyyzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyyzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyyzzz_0[i] = 4.0 * g_xxxyy_0_yyyyzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyyzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyyzzzz_0[i] = 4.0 * g_xxxyy_0_yyyzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyyzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yyzzzzz_0[i] = 4.0 * g_xxxyy_0_yyzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yyzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_yzzzzzz_0[i] = 4.0 * g_xxxyy_0_yzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_yzzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxxxyy_0_zzzzzzz_0[i] = 4.0 * g_xxxyy_0_zzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxyy_0_zzzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 144-180 components of targeted buffer : KSK

    auto g_xxxxxyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 144);

    auto g_xxxxxyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 145);

    auto g_xxxxxyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 146);

    auto g_xxxxxyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 147);

    auto g_xxxxxyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 148);

    auto g_xxxxxyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 149);

    auto g_xxxxxyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 150);

    auto g_xxxxxyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 151);

    auto g_xxxxxyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 152);

    auto g_xxxxxyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 153);

    auto g_xxxxxyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 154);

    auto g_xxxxxyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 155);

    auto g_xxxxxyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 156);

    auto g_xxxxxyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 157);

    auto g_xxxxxyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 158);

    auto g_xxxxxyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 159);

    auto g_xxxxxyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 160);

    auto g_xxxxxyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 161);

    auto g_xxxxxyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 162);

    auto g_xxxxxyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 163);

    auto g_xxxxxyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 164);

    auto g_xxxxxyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 165);

    auto g_xxxxxyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 166);

    auto g_xxxxxyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 167);

    auto g_xxxxxyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 168);

    auto g_xxxxxyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 169);

    auto g_xxxxxyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 170);

    auto g_xxxxxyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 171);

    auto g_xxxxxyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 172);

    auto g_xxxxxyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 173);

    auto g_xxxxxyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 174);

    auto g_xxxxxyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 175);

    auto g_xxxxxyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 176);

    auto g_xxxxxyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 177);

    auto g_xxxxxyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 178);

    auto g_xxxxxyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 179);

    #pragma omp simd aligned(g_xxxxxy_0_xxxxxxy_1, g_xxxxxy_0_xxxxxyy_1, g_xxxxxy_0_xxxxyyy_1, g_xxxxxy_0_xxxyyyy_1, g_xxxxxy_0_xxyyyyy_1, g_xxxxxy_0_xyyyyyy_1, g_xxxxxy_0_yyyyyyy_1, g_xxxxxyz_0_xxxxxxx_0, g_xxxxxyz_0_xxxxxxy_0, g_xxxxxyz_0_xxxxxxz_0, g_xxxxxyz_0_xxxxxyy_0, g_xxxxxyz_0_xxxxxyz_0, g_xxxxxyz_0_xxxxxzz_0, g_xxxxxyz_0_xxxxyyy_0, g_xxxxxyz_0_xxxxyyz_0, g_xxxxxyz_0_xxxxyzz_0, g_xxxxxyz_0_xxxxzzz_0, g_xxxxxyz_0_xxxyyyy_0, g_xxxxxyz_0_xxxyyyz_0, g_xxxxxyz_0_xxxyyzz_0, g_xxxxxyz_0_xxxyzzz_0, g_xxxxxyz_0_xxxzzzz_0, g_xxxxxyz_0_xxyyyyy_0, g_xxxxxyz_0_xxyyyyz_0, g_xxxxxyz_0_xxyyyzz_0, g_xxxxxyz_0_xxyyzzz_0, g_xxxxxyz_0_xxyzzzz_0, g_xxxxxyz_0_xxzzzzz_0, g_xxxxxyz_0_xyyyyyy_0, g_xxxxxyz_0_xyyyyyz_0, g_xxxxxyz_0_xyyyyzz_0, g_xxxxxyz_0_xyyyzzz_0, g_xxxxxyz_0_xyyzzzz_0, g_xxxxxyz_0_xyzzzzz_0, g_xxxxxyz_0_xzzzzzz_0, g_xxxxxyz_0_yyyyyyy_0, g_xxxxxyz_0_yyyyyyz_0, g_xxxxxyz_0_yyyyyzz_0, g_xxxxxyz_0_yyyyzzz_0, g_xxxxxyz_0_yyyzzzz_0, g_xxxxxyz_0_yyzzzzz_0, g_xxxxxyz_0_yzzzzzz_0, g_xxxxxyz_0_zzzzzzz_0, g_xxxxxz_0_xxxxxxx_1, g_xxxxxz_0_xxxxxxz_1, g_xxxxxz_0_xxxxxyz_1, g_xxxxxz_0_xxxxxz_1, g_xxxxxz_0_xxxxxzz_1, g_xxxxxz_0_xxxxyyz_1, g_xxxxxz_0_xxxxyz_1, g_xxxxxz_0_xxxxyzz_1, g_xxxxxz_0_xxxxzz_1, g_xxxxxz_0_xxxxzzz_1, g_xxxxxz_0_xxxyyyz_1, g_xxxxxz_0_xxxyyz_1, g_xxxxxz_0_xxxyyzz_1, g_xxxxxz_0_xxxyzz_1, g_xxxxxz_0_xxxyzzz_1, g_xxxxxz_0_xxxzzz_1, g_xxxxxz_0_xxxzzzz_1, g_xxxxxz_0_xxyyyyz_1, g_xxxxxz_0_xxyyyz_1, g_xxxxxz_0_xxyyyzz_1, g_xxxxxz_0_xxyyzz_1, g_xxxxxz_0_xxyyzzz_1, g_xxxxxz_0_xxyzzz_1, g_xxxxxz_0_xxyzzzz_1, g_xxxxxz_0_xxzzzz_1, g_xxxxxz_0_xxzzzzz_1, g_xxxxxz_0_xyyyyyz_1, g_xxxxxz_0_xyyyyz_1, g_xxxxxz_0_xyyyyzz_1, g_xxxxxz_0_xyyyzz_1, g_xxxxxz_0_xyyyzzz_1, g_xxxxxz_0_xyyzzz_1, g_xxxxxz_0_xyyzzzz_1, g_xxxxxz_0_xyzzzz_1, g_xxxxxz_0_xyzzzzz_1, g_xxxxxz_0_xzzzzz_1, g_xxxxxz_0_xzzzzzz_1, g_xxxxxz_0_yyyyyyz_1, g_xxxxxz_0_yyyyyz_1, g_xxxxxz_0_yyyyyzz_1, g_xxxxxz_0_yyyyzz_1, g_xxxxxz_0_yyyyzzz_1, g_xxxxxz_0_yyyzzz_1, g_xxxxxz_0_yyyzzzz_1, g_xxxxxz_0_yyzzzz_1, g_xxxxxz_0_yyzzzzz_1, g_xxxxxz_0_yzzzzz_1, g_xxxxxz_0_yzzzzzz_1, g_xxxxxz_0_zzzzzz_1, g_xxxxxz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxyz_0_xxxxxxx_0[i] = g_xxxxxz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxxxy_0[i] = g_xxxxxy_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxxxxz_0[i] = g_xxxxxz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxxyy_0[i] = g_xxxxxy_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxxxyz_0[i] = g_xxxxxz_0_xxxxxz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxxxyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxxzz_0[i] = g_xxxxxz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxyyy_0[i] = g_xxxxxy_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxxyyz_0[i] = 2.0 * g_xxxxxz_0_xxxxyz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxxyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxyzz_0[i] = g_xxxxxz_0_xxxxzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxxyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxxzzz_0[i] = g_xxxxxz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxyyyy_0[i] = g_xxxxxy_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxxyyyz_0[i] = 3.0 * g_xxxxxz_0_xxxyyz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxyyzz_0[i] = 2.0 * g_xxxxxz_0_xxxyzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxyzzz_0[i] = g_xxxxxz_0_xxxzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxxyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxxzzzz_0[i] = g_xxxxxz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyyyyy_0[i] = g_xxxxxy_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xxyyyyz_0[i] = 4.0 * g_xxxxxz_0_xxyyyz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyyyzz_0[i] = 3.0 * g_xxxxxz_0_xxyyzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyyzzz_0[i] = 2.0 * g_xxxxxz_0_xxyzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxyzzzz_0[i] = g_xxxxxz_0_xxzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xxyzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xxzzzzz_0[i] = g_xxxxxz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyyyyy_0[i] = g_xxxxxy_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_xyyyyyz_0[i] = 5.0 * g_xxxxxz_0_xyyyyz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyyyzz_0[i] = 4.0 * g_xxxxxz_0_xyyyzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyyzzz_0[i] = 3.0 * g_xxxxxz_0_xyyzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyyzzzz_0[i] = 2.0 * g_xxxxxz_0_xyzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyyzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xyzzzzz_0[i] = g_xxxxxz_0_xzzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_xyzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_xzzzzzz_0[i] = g_xxxxxz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyyyyy_0[i] = g_xxxxxy_0_yyyyyyy_1[i] * wa_z[i];

        g_xxxxxyz_0_yyyyyyz_0[i] = 6.0 * g_xxxxxz_0_yyyyyz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyyyyz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyyyzz_0[i] = 5.0 * g_xxxxxz_0_yyyyzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyyyzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyyzzz_0[i] = 4.0 * g_xxxxxz_0_yyyzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyyzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyyzzzz_0[i] = 3.0 * g_xxxxxz_0_yyzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyyzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yyzzzzz_0[i] = 2.0 * g_xxxxxz_0_yzzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yyzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_yzzzzzz_0[i] = g_xxxxxz_0_zzzzzz_1[i] * fi_acd_0 + g_xxxxxz_0_yzzzzzz_1[i] * wa_y[i];

        g_xxxxxyz_0_zzzzzzz_0[i] = g_xxxxxz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 180-216 components of targeted buffer : KSK

    auto g_xxxxxzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 180);

    auto g_xxxxxzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 181);

    auto g_xxxxxzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 182);

    auto g_xxxxxzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 183);

    auto g_xxxxxzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 184);

    auto g_xxxxxzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 185);

    auto g_xxxxxzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 186);

    auto g_xxxxxzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 187);

    auto g_xxxxxzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 188);

    auto g_xxxxxzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 189);

    auto g_xxxxxzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 190);

    auto g_xxxxxzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 191);

    auto g_xxxxxzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 192);

    auto g_xxxxxzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 193);

    auto g_xxxxxzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 194);

    auto g_xxxxxzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 195);

    auto g_xxxxxzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 196);

    auto g_xxxxxzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 197);

    auto g_xxxxxzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 198);

    auto g_xxxxxzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 199);

    auto g_xxxxxzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 200);

    auto g_xxxxxzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 201);

    auto g_xxxxxzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 202);

    auto g_xxxxxzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 203);

    auto g_xxxxxzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 204);

    auto g_xxxxxzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 205);

    auto g_xxxxxzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 206);

    auto g_xxxxxzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 207);

    auto g_xxxxxzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 208);

    auto g_xxxxxzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 209);

    auto g_xxxxxzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 210);

    auto g_xxxxxzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 211);

    auto g_xxxxxzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 212);

    auto g_xxxxxzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 213);

    auto g_xxxxxzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 214);

    auto g_xxxxxzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 215);

    #pragma omp simd aligned(g_xxxxx_0_xxxxxxx_0, g_xxxxx_0_xxxxxxx_1, g_xxxxx_0_xxxxxxy_0, g_xxxxx_0_xxxxxxy_1, g_xxxxx_0_xxxxxyy_0, g_xxxxx_0_xxxxxyy_1, g_xxxxx_0_xxxxyyy_0, g_xxxxx_0_xxxxyyy_1, g_xxxxx_0_xxxyyyy_0, g_xxxxx_0_xxxyyyy_1, g_xxxxx_0_xxyyyyy_0, g_xxxxx_0_xxyyyyy_1, g_xxxxx_0_xyyyyyy_0, g_xxxxx_0_xyyyyyy_1, g_xxxxxz_0_xxxxxxx_1, g_xxxxxz_0_xxxxxxy_1, g_xxxxxz_0_xxxxxyy_1, g_xxxxxz_0_xxxxyyy_1, g_xxxxxz_0_xxxyyyy_1, g_xxxxxz_0_xxyyyyy_1, g_xxxxxz_0_xyyyyyy_1, g_xxxxxzz_0_xxxxxxx_0, g_xxxxxzz_0_xxxxxxy_0, g_xxxxxzz_0_xxxxxxz_0, g_xxxxxzz_0_xxxxxyy_0, g_xxxxxzz_0_xxxxxyz_0, g_xxxxxzz_0_xxxxxzz_0, g_xxxxxzz_0_xxxxyyy_0, g_xxxxxzz_0_xxxxyyz_0, g_xxxxxzz_0_xxxxyzz_0, g_xxxxxzz_0_xxxxzzz_0, g_xxxxxzz_0_xxxyyyy_0, g_xxxxxzz_0_xxxyyyz_0, g_xxxxxzz_0_xxxyyzz_0, g_xxxxxzz_0_xxxyzzz_0, g_xxxxxzz_0_xxxzzzz_0, g_xxxxxzz_0_xxyyyyy_0, g_xxxxxzz_0_xxyyyyz_0, g_xxxxxzz_0_xxyyyzz_0, g_xxxxxzz_0_xxyyzzz_0, g_xxxxxzz_0_xxyzzzz_0, g_xxxxxzz_0_xxzzzzz_0, g_xxxxxzz_0_xyyyyyy_0, g_xxxxxzz_0_xyyyyyz_0, g_xxxxxzz_0_xyyyyzz_0, g_xxxxxzz_0_xyyyzzz_0, g_xxxxxzz_0_xyyzzzz_0, g_xxxxxzz_0_xyzzzzz_0, g_xxxxxzz_0_xzzzzzz_0, g_xxxxxzz_0_yyyyyyy_0, g_xxxxxzz_0_yyyyyyz_0, g_xxxxxzz_0_yyyyyzz_0, g_xxxxxzz_0_yyyyzzz_0, g_xxxxxzz_0_yyyzzzz_0, g_xxxxxzz_0_yyzzzzz_0, g_xxxxxzz_0_yzzzzzz_0, g_xxxxxzz_0_zzzzzzz_0, g_xxxxzz_0_xxxxxxz_1, g_xxxxzz_0_xxxxxyz_1, g_xxxxzz_0_xxxxxz_1, g_xxxxzz_0_xxxxxzz_1, g_xxxxzz_0_xxxxyyz_1, g_xxxxzz_0_xxxxyz_1, g_xxxxzz_0_xxxxyzz_1, g_xxxxzz_0_xxxxzz_1, g_xxxxzz_0_xxxxzzz_1, g_xxxxzz_0_xxxyyyz_1, g_xxxxzz_0_xxxyyz_1, g_xxxxzz_0_xxxyyzz_1, g_xxxxzz_0_xxxyzz_1, g_xxxxzz_0_xxxyzzz_1, g_xxxxzz_0_xxxzzz_1, g_xxxxzz_0_xxxzzzz_1, g_xxxxzz_0_xxyyyyz_1, g_xxxxzz_0_xxyyyz_1, g_xxxxzz_0_xxyyyzz_1, g_xxxxzz_0_xxyyzz_1, g_xxxxzz_0_xxyyzzz_1, g_xxxxzz_0_xxyzzz_1, g_xxxxzz_0_xxyzzzz_1, g_xxxxzz_0_xxzzzz_1, g_xxxxzz_0_xxzzzzz_1, g_xxxxzz_0_xyyyyyz_1, g_xxxxzz_0_xyyyyz_1, g_xxxxzz_0_xyyyyzz_1, g_xxxxzz_0_xyyyzz_1, g_xxxxzz_0_xyyyzzz_1, g_xxxxzz_0_xyyzzz_1, g_xxxxzz_0_xyyzzzz_1, g_xxxxzz_0_xyzzzz_1, g_xxxxzz_0_xyzzzzz_1, g_xxxxzz_0_xzzzzz_1, g_xxxxzz_0_xzzzzzz_1, g_xxxxzz_0_yyyyyyy_1, g_xxxxzz_0_yyyyyyz_1, g_xxxxzz_0_yyyyyz_1, g_xxxxzz_0_yyyyyzz_1, g_xxxxzz_0_yyyyzz_1, g_xxxxzz_0_yyyyzzz_1, g_xxxxzz_0_yyyzzz_1, g_xxxxzz_0_yyyzzzz_1, g_xxxxzz_0_yyzzzz_1, g_xxxxzz_0_yyzzzzz_1, g_xxxxzz_0_yzzzzz_1, g_xxxxzz_0_yzzzzzz_1, g_xxxxzz_0_zzzzzz_1, g_xxxxzz_0_zzzzzzz_1, g_xxxzz_0_xxxxxxz_0, g_xxxzz_0_xxxxxxz_1, g_xxxzz_0_xxxxxyz_0, g_xxxzz_0_xxxxxyz_1, g_xxxzz_0_xxxxxzz_0, g_xxxzz_0_xxxxxzz_1, g_xxxzz_0_xxxxyyz_0, g_xxxzz_0_xxxxyyz_1, g_xxxzz_0_xxxxyzz_0, g_xxxzz_0_xxxxyzz_1, g_xxxzz_0_xxxxzzz_0, g_xxxzz_0_xxxxzzz_1, g_xxxzz_0_xxxyyyz_0, g_xxxzz_0_xxxyyyz_1, g_xxxzz_0_xxxyyzz_0, g_xxxzz_0_xxxyyzz_1, g_xxxzz_0_xxxyzzz_0, g_xxxzz_0_xxxyzzz_1, g_xxxzz_0_xxxzzzz_0, g_xxxzz_0_xxxzzzz_1, g_xxxzz_0_xxyyyyz_0, g_xxxzz_0_xxyyyyz_1, g_xxxzz_0_xxyyyzz_0, g_xxxzz_0_xxyyyzz_1, g_xxxzz_0_xxyyzzz_0, g_xxxzz_0_xxyyzzz_1, g_xxxzz_0_xxyzzzz_0, g_xxxzz_0_xxyzzzz_1, g_xxxzz_0_xxzzzzz_0, g_xxxzz_0_xxzzzzz_1, g_xxxzz_0_xyyyyyz_0, g_xxxzz_0_xyyyyyz_1, g_xxxzz_0_xyyyyzz_0, g_xxxzz_0_xyyyyzz_1, g_xxxzz_0_xyyyzzz_0, g_xxxzz_0_xyyyzzz_1, g_xxxzz_0_xyyzzzz_0, g_xxxzz_0_xyyzzzz_1, g_xxxzz_0_xyzzzzz_0, g_xxxzz_0_xyzzzzz_1, g_xxxzz_0_xzzzzzz_0, g_xxxzz_0_xzzzzzz_1, g_xxxzz_0_yyyyyyy_0, g_xxxzz_0_yyyyyyy_1, g_xxxzz_0_yyyyyyz_0, g_xxxzz_0_yyyyyyz_1, g_xxxzz_0_yyyyyzz_0, g_xxxzz_0_yyyyyzz_1, g_xxxzz_0_yyyyzzz_0, g_xxxzz_0_yyyyzzz_1, g_xxxzz_0_yyyzzzz_0, g_xxxzz_0_yyyzzzz_1, g_xxxzz_0_yyzzzzz_0, g_xxxzz_0_yyzzzzz_1, g_xxxzz_0_yzzzzzz_0, g_xxxzz_0_yzzzzzz_1, g_xxxzz_0_zzzzzzz_0, g_xxxzz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxzz_0_xxxxxxx_0[i] = g_xxxxx_0_xxxxxxx_0[i] * fbe_0 - g_xxxxx_0_xxxxxxx_1[i] * fz_be_0 + g_xxxxxz_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxxxy_0[i] = g_xxxxx_0_xxxxxxy_0[i] * fbe_0 - g_xxxxx_0_xxxxxxy_1[i] * fz_be_0 + g_xxxxxz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxxxz_0[i] = 4.0 * g_xxxzz_0_xxxxxxz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xxxxzz_0_xxxxxz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxxyy_0[i] = g_xxxxx_0_xxxxxyy_0[i] * fbe_0 - g_xxxxx_0_xxxxxyy_1[i] * fz_be_0 + g_xxxxxz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxxyz_0[i] = 4.0 * g_xxxzz_0_xxxxxyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxxxzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxxzz_0[i] = 4.0 * g_xxxzz_0_xxxxxzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xxxxzz_0_xxxxzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxyyy_0[i] = g_xxxxx_0_xxxxyyy_0[i] * fbe_0 - g_xxxxx_0_xxxxyyy_1[i] * fz_be_0 + g_xxxxxz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxxyyz_0[i] = 4.0 * g_xxxzz_0_xxxxyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxyzz_0[i] = 4.0 * g_xxxzz_0_xxxxyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxxzzz_0[i] = 4.0 * g_xxxzz_0_xxxxzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xxxxzz_0_xxxzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxyyyy_0[i] = g_xxxxx_0_xxxyyyy_0[i] * fbe_0 - g_xxxxx_0_xxxyyyy_1[i] * fz_be_0 + g_xxxxxz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxxyyyz_0[i] = 4.0 * g_xxxzz_0_xxxyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxyyzz_0[i] = 4.0 * g_xxxzz_0_xxxyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxyzzz_0[i] = 4.0 * g_xxxzz_0_xxxyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxxzzzz_0[i] = 4.0 * g_xxxzz_0_xxxzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xxxxzz_0_xxzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyyyyy_0[i] = g_xxxxx_0_xxyyyyy_0[i] * fbe_0 - g_xxxxx_0_xxyyyyy_1[i] * fz_be_0 + g_xxxxxz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xxyyyyz_0[i] = 4.0 * g_xxxzz_0_xxyyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyyyzz_0[i] = 4.0 * g_xxxzz_0_xxyyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyyzzz_0[i] = 4.0 * g_xxxzz_0_xxyyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxyzzzz_0[i] = 4.0 * g_xxxzz_0_xxyzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xxzzzzz_0[i] = 4.0 * g_xxxzz_0_xxzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_0_xzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyyyyy_0[i] = g_xxxxx_0_xyyyyyy_0[i] * fbe_0 - g_xxxxx_0_xyyyyyy_1[i] * fz_be_0 + g_xxxxxz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxxxzz_0_xyyyyyz_0[i] = 4.0 * g_xxxzz_0_xyyyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyyyyz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyyyzz_0[i] = 4.0 * g_xxxzz_0_xyyyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyyyzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyyzzz_0[i] = 4.0 * g_xxxzz_0_xyyyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyyzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyyzzzz_0[i] = 4.0 * g_xxxzz_0_xyyzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyyzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xyzzzzz_0[i] = 4.0 * g_xxxzz_0_xyzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xyzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_xzzzzzz_0[i] = 4.0 * g_xxxzz_0_xzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_xzzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_zzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyyyy_0[i] = 4.0 * g_xxxzz_0_yyyyyyy_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyyyy_1[i] * fz_be_0 + g_xxxxzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyyyz_0[i] = 4.0 * g_xxxzz_0_yyyyyyz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyyyz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyyzz_0[i] = 4.0 * g_xxxzz_0_yyyyyzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyyzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyyzzz_0[i] = 4.0 * g_xxxzz_0_yyyyzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyyzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyyzzzz_0[i] = 4.0 * g_xxxzz_0_yyyzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyyzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yyzzzzz_0[i] = 4.0 * g_xxxzz_0_yyzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yyzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_yzzzzzz_0[i] = 4.0 * g_xxxzz_0_yzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_yzzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxxxzz_0_zzzzzzz_0[i] = 4.0 * g_xxxzz_0_zzzzzzz_0[i] * fbe_0 - 4.0 * g_xxxzz_0_zzzzzzz_1[i] * fz_be_0 + g_xxxxzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 216-252 components of targeted buffer : KSK

    auto g_xxxxyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 216);

    auto g_xxxxyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 217);

    auto g_xxxxyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 218);

    auto g_xxxxyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 219);

    auto g_xxxxyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 220);

    auto g_xxxxyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 221);

    auto g_xxxxyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 222);

    auto g_xxxxyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 223);

    auto g_xxxxyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 224);

    auto g_xxxxyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 225);

    auto g_xxxxyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 226);

    auto g_xxxxyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 227);

    auto g_xxxxyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 228);

    auto g_xxxxyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 229);

    auto g_xxxxyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 230);

    auto g_xxxxyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 231);

    auto g_xxxxyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 232);

    auto g_xxxxyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 233);

    auto g_xxxxyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 234);

    auto g_xxxxyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 235);

    auto g_xxxxyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 236);

    auto g_xxxxyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 237);

    auto g_xxxxyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 238);

    auto g_xxxxyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 239);

    auto g_xxxxyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 240);

    auto g_xxxxyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 241);

    auto g_xxxxyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 242);

    auto g_xxxxyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 243);

    auto g_xxxxyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 244);

    auto g_xxxxyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 245);

    auto g_xxxxyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 246);

    auto g_xxxxyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 247);

    auto g_xxxxyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 248);

    auto g_xxxxyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 249);

    auto g_xxxxyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 250);

    auto g_xxxxyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 251);

    #pragma omp simd aligned(g_xxxxy_0_xxxxxxx_0, g_xxxxy_0_xxxxxxx_1, g_xxxxy_0_xxxxxxz_0, g_xxxxy_0_xxxxxxz_1, g_xxxxy_0_xxxxxzz_0, g_xxxxy_0_xxxxxzz_1, g_xxxxy_0_xxxxzzz_0, g_xxxxy_0_xxxxzzz_1, g_xxxxy_0_xxxzzzz_0, g_xxxxy_0_xxxzzzz_1, g_xxxxy_0_xxzzzzz_0, g_xxxxy_0_xxzzzzz_1, g_xxxxy_0_xzzzzzz_0, g_xxxxy_0_xzzzzzz_1, g_xxxxyy_0_xxxxxxx_1, g_xxxxyy_0_xxxxxxz_1, g_xxxxyy_0_xxxxxzz_1, g_xxxxyy_0_xxxxzzz_1, g_xxxxyy_0_xxxzzzz_1, g_xxxxyy_0_xxzzzzz_1, g_xxxxyy_0_xzzzzzz_1, g_xxxxyyy_0_xxxxxxx_0, g_xxxxyyy_0_xxxxxxy_0, g_xxxxyyy_0_xxxxxxz_0, g_xxxxyyy_0_xxxxxyy_0, g_xxxxyyy_0_xxxxxyz_0, g_xxxxyyy_0_xxxxxzz_0, g_xxxxyyy_0_xxxxyyy_0, g_xxxxyyy_0_xxxxyyz_0, g_xxxxyyy_0_xxxxyzz_0, g_xxxxyyy_0_xxxxzzz_0, g_xxxxyyy_0_xxxyyyy_0, g_xxxxyyy_0_xxxyyyz_0, g_xxxxyyy_0_xxxyyzz_0, g_xxxxyyy_0_xxxyzzz_0, g_xxxxyyy_0_xxxzzzz_0, g_xxxxyyy_0_xxyyyyy_0, g_xxxxyyy_0_xxyyyyz_0, g_xxxxyyy_0_xxyyyzz_0, g_xxxxyyy_0_xxyyzzz_0, g_xxxxyyy_0_xxyzzzz_0, g_xxxxyyy_0_xxzzzzz_0, g_xxxxyyy_0_xyyyyyy_0, g_xxxxyyy_0_xyyyyyz_0, g_xxxxyyy_0_xyyyyzz_0, g_xxxxyyy_0_xyyyzzz_0, g_xxxxyyy_0_xyyzzzz_0, g_xxxxyyy_0_xyzzzzz_0, g_xxxxyyy_0_xzzzzzz_0, g_xxxxyyy_0_yyyyyyy_0, g_xxxxyyy_0_yyyyyyz_0, g_xxxxyyy_0_yyyyyzz_0, g_xxxxyyy_0_yyyyzzz_0, g_xxxxyyy_0_yyyzzzz_0, g_xxxxyyy_0_yyzzzzz_0, g_xxxxyyy_0_yzzzzzz_0, g_xxxxyyy_0_zzzzzzz_0, g_xxxyyy_0_xxxxxxy_1, g_xxxyyy_0_xxxxxy_1, g_xxxyyy_0_xxxxxyy_1, g_xxxyyy_0_xxxxxyz_1, g_xxxyyy_0_xxxxyy_1, g_xxxyyy_0_xxxxyyy_1, g_xxxyyy_0_xxxxyyz_1, g_xxxyyy_0_xxxxyz_1, g_xxxyyy_0_xxxxyzz_1, g_xxxyyy_0_xxxyyy_1, g_xxxyyy_0_xxxyyyy_1, g_xxxyyy_0_xxxyyyz_1, g_xxxyyy_0_xxxyyz_1, g_xxxyyy_0_xxxyyzz_1, g_xxxyyy_0_xxxyzz_1, g_xxxyyy_0_xxxyzzz_1, g_xxxyyy_0_xxyyyy_1, g_xxxyyy_0_xxyyyyy_1, g_xxxyyy_0_xxyyyyz_1, g_xxxyyy_0_xxyyyz_1, g_xxxyyy_0_xxyyyzz_1, g_xxxyyy_0_xxyyzz_1, g_xxxyyy_0_xxyyzzz_1, g_xxxyyy_0_xxyzzz_1, g_xxxyyy_0_xxyzzzz_1, g_xxxyyy_0_xyyyyy_1, g_xxxyyy_0_xyyyyyy_1, g_xxxyyy_0_xyyyyyz_1, g_xxxyyy_0_xyyyyz_1, g_xxxyyy_0_xyyyyzz_1, g_xxxyyy_0_xyyyzz_1, g_xxxyyy_0_xyyyzzz_1, g_xxxyyy_0_xyyzzz_1, g_xxxyyy_0_xyyzzzz_1, g_xxxyyy_0_xyzzzz_1, g_xxxyyy_0_xyzzzzz_1, g_xxxyyy_0_yyyyyy_1, g_xxxyyy_0_yyyyyyy_1, g_xxxyyy_0_yyyyyyz_1, g_xxxyyy_0_yyyyyz_1, g_xxxyyy_0_yyyyyzz_1, g_xxxyyy_0_yyyyzz_1, g_xxxyyy_0_yyyyzzz_1, g_xxxyyy_0_yyyzzz_1, g_xxxyyy_0_yyyzzzz_1, g_xxxyyy_0_yyzzzz_1, g_xxxyyy_0_yyzzzzz_1, g_xxxyyy_0_yzzzzz_1, g_xxxyyy_0_yzzzzzz_1, g_xxxyyy_0_zzzzzzz_1, g_xxyyy_0_xxxxxxy_0, g_xxyyy_0_xxxxxxy_1, g_xxyyy_0_xxxxxyy_0, g_xxyyy_0_xxxxxyy_1, g_xxyyy_0_xxxxxyz_0, g_xxyyy_0_xxxxxyz_1, g_xxyyy_0_xxxxyyy_0, g_xxyyy_0_xxxxyyy_1, g_xxyyy_0_xxxxyyz_0, g_xxyyy_0_xxxxyyz_1, g_xxyyy_0_xxxxyzz_0, g_xxyyy_0_xxxxyzz_1, g_xxyyy_0_xxxyyyy_0, g_xxyyy_0_xxxyyyy_1, g_xxyyy_0_xxxyyyz_0, g_xxyyy_0_xxxyyyz_1, g_xxyyy_0_xxxyyzz_0, g_xxyyy_0_xxxyyzz_1, g_xxyyy_0_xxxyzzz_0, g_xxyyy_0_xxxyzzz_1, g_xxyyy_0_xxyyyyy_0, g_xxyyy_0_xxyyyyy_1, g_xxyyy_0_xxyyyyz_0, g_xxyyy_0_xxyyyyz_1, g_xxyyy_0_xxyyyzz_0, g_xxyyy_0_xxyyyzz_1, g_xxyyy_0_xxyyzzz_0, g_xxyyy_0_xxyyzzz_1, g_xxyyy_0_xxyzzzz_0, g_xxyyy_0_xxyzzzz_1, g_xxyyy_0_xyyyyyy_0, g_xxyyy_0_xyyyyyy_1, g_xxyyy_0_xyyyyyz_0, g_xxyyy_0_xyyyyyz_1, g_xxyyy_0_xyyyyzz_0, g_xxyyy_0_xyyyyzz_1, g_xxyyy_0_xyyyzzz_0, g_xxyyy_0_xyyyzzz_1, g_xxyyy_0_xyyzzzz_0, g_xxyyy_0_xyyzzzz_1, g_xxyyy_0_xyzzzzz_0, g_xxyyy_0_xyzzzzz_1, g_xxyyy_0_yyyyyyy_0, g_xxyyy_0_yyyyyyy_1, g_xxyyy_0_yyyyyyz_0, g_xxyyy_0_yyyyyyz_1, g_xxyyy_0_yyyyyzz_0, g_xxyyy_0_yyyyyzz_1, g_xxyyy_0_yyyyzzz_0, g_xxyyy_0_yyyyzzz_1, g_xxyyy_0_yyyzzzz_0, g_xxyyy_0_yyyzzzz_1, g_xxyyy_0_yyzzzzz_0, g_xxyyy_0_yyzzzzz_1, g_xxyyy_0_yzzzzzz_0, g_xxyyy_0_yzzzzzz_1, g_xxyyy_0_zzzzzzz_0, g_xxyyy_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxyyy_0_xxxxxxx_0[i] = 2.0 * g_xxxxy_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxxxx_1[i] * fz_be_0 + g_xxxxyy_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxxxxy_0[i] = 3.0 * g_xxyyy_0_xxxxxxy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xxxyyy_0_xxxxxy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxxxz_0[i] = 2.0 * g_xxxxy_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxxxz_1[i] * fz_be_0 + g_xxxxyy_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxxxyy_0[i] = 3.0 * g_xxyyy_0_xxxxxyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xxxyyy_0_xxxxyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxxyz_0[i] = 3.0 * g_xxyyy_0_xxxxxyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxxyyy_0_xxxxyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxxzz_0[i] = 2.0 * g_xxxxy_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxxzz_1[i] * fz_be_0 + g_xxxxyy_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxxyyy_0[i] = 3.0 * g_xxyyy_0_xxxxyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xxxyyy_0_xxxyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxyyz_0[i] = 3.0 * g_xxyyy_0_xxxxyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxxyyy_0_xxxyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxyzz_0[i] = 3.0 * g_xxyyy_0_xxxxyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxxyyy_0_xxxyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxxzzz_0[i] = 2.0 * g_xxxxy_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxxzzz_1[i] * fz_be_0 + g_xxxxyy_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxxyyyy_0[i] = 3.0 * g_xxyyy_0_xxxyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxyyyz_0[i] = 3.0 * g_xxyyy_0_xxxyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxyyzz_0[i] = 3.0 * g_xxyyy_0_xxxyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxyzzz_0[i] = 3.0 * g_xxyyy_0_xxxyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxxyyy_0_xxyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxxzzzz_0[i] = 2.0 * g_xxxxy_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxxzzzz_1[i] * fz_be_0 + g_xxxxyy_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xxyyyyy_0[i] = 3.0 * g_xxyyy_0_xxyyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyyyyz_0[i] = 3.0 * g_xxyyy_0_xxyyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyyyzz_0[i] = 3.0 * g_xxyyy_0_xxyyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyyzzz_0[i] = 3.0 * g_xxyyy_0_xxyyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxyzzzz_0[i] = 3.0 * g_xxyyy_0_xxyzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxyyy_0_xyzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xxzzzzz_0[i] = 2.0 * g_xxxxy_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xxzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_xyyyyyy_0[i] = 3.0 * g_xxyyy_0_xyyyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyyyy_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyyyyz_0[i] = 3.0 * g_xxyyy_0_xyyyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyyyz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyyyzz_0[i] = 3.0 * g_xxyyy_0_xyyyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyyzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyyzzz_0[i] = 3.0 * g_xxyyy_0_xyyyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyyzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyyzzzz_0[i] = 3.0 * g_xxyyy_0_xyyzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyyzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xyzzzzz_0[i] = 3.0 * g_xxyyy_0_xyzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_xyzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_xzzzzzz_0[i] = 2.0 * g_xxxxy_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_xxxxy_0_xzzzzzz_1[i] * fz_be_0 + g_xxxxyy_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxxyyy_0_yyyyyyy_0[i] = 3.0 * g_xxyyy_0_yyyyyyy_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyyyy_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyyyyz_0[i] = 3.0 * g_xxyyy_0_yyyyyyz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyyyz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyyyzz_0[i] = 3.0 * g_xxyyy_0_yyyyyzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyyzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyyzzz_0[i] = 3.0 * g_xxyyy_0_yyyyzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyyzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyyzzzz_0[i] = 3.0 * g_xxyyy_0_yyyzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyyzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yyzzzzz_0[i] = 3.0 * g_xxyyy_0_yyzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yyzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_yzzzzzz_0[i] = 3.0 * g_xxyyy_0_yzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_yzzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxxyyy_0_zzzzzzz_0[i] = 3.0 * g_xxyyy_0_zzzzzzz_0[i] * fbe_0 - 3.0 * g_xxyyy_0_zzzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 252-288 components of targeted buffer : KSK

    auto g_xxxxyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 252);

    auto g_xxxxyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 253);

    auto g_xxxxyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 254);

    auto g_xxxxyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 255);

    auto g_xxxxyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 256);

    auto g_xxxxyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 257);

    auto g_xxxxyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 258);

    auto g_xxxxyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 259);

    auto g_xxxxyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 260);

    auto g_xxxxyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 261);

    auto g_xxxxyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 262);

    auto g_xxxxyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 263);

    auto g_xxxxyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 264);

    auto g_xxxxyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 265);

    auto g_xxxxyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 266);

    auto g_xxxxyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 267);

    auto g_xxxxyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 268);

    auto g_xxxxyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 269);

    auto g_xxxxyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 270);

    auto g_xxxxyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 271);

    auto g_xxxxyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 272);

    auto g_xxxxyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 273);

    auto g_xxxxyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 274);

    auto g_xxxxyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 275);

    auto g_xxxxyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 276);

    auto g_xxxxyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 277);

    auto g_xxxxyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 278);

    auto g_xxxxyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 279);

    auto g_xxxxyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 280);

    auto g_xxxxyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 281);

    auto g_xxxxyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 282);

    auto g_xxxxyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 283);

    auto g_xxxxyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 284);

    auto g_xxxxyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 285);

    auto g_xxxxyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 286);

    auto g_xxxxyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 287);

    #pragma omp simd aligned(g_xxxxyy_0_xxxxxx_1, g_xxxxyy_0_xxxxxxx_1, g_xxxxyy_0_xxxxxxy_1, g_xxxxyy_0_xxxxxxz_1, g_xxxxyy_0_xxxxxy_1, g_xxxxyy_0_xxxxxyy_1, g_xxxxyy_0_xxxxxyz_1, g_xxxxyy_0_xxxxxz_1, g_xxxxyy_0_xxxxxzz_1, g_xxxxyy_0_xxxxyy_1, g_xxxxyy_0_xxxxyyy_1, g_xxxxyy_0_xxxxyyz_1, g_xxxxyy_0_xxxxyz_1, g_xxxxyy_0_xxxxyzz_1, g_xxxxyy_0_xxxxzz_1, g_xxxxyy_0_xxxxzzz_1, g_xxxxyy_0_xxxyyy_1, g_xxxxyy_0_xxxyyyy_1, g_xxxxyy_0_xxxyyyz_1, g_xxxxyy_0_xxxyyz_1, g_xxxxyy_0_xxxyyzz_1, g_xxxxyy_0_xxxyzz_1, g_xxxxyy_0_xxxyzzz_1, g_xxxxyy_0_xxxzzz_1, g_xxxxyy_0_xxxzzzz_1, g_xxxxyy_0_xxyyyy_1, g_xxxxyy_0_xxyyyyy_1, g_xxxxyy_0_xxyyyyz_1, g_xxxxyy_0_xxyyyz_1, g_xxxxyy_0_xxyyyzz_1, g_xxxxyy_0_xxyyzz_1, g_xxxxyy_0_xxyyzzz_1, g_xxxxyy_0_xxyzzz_1, g_xxxxyy_0_xxyzzzz_1, g_xxxxyy_0_xxzzzz_1, g_xxxxyy_0_xxzzzzz_1, g_xxxxyy_0_xyyyyy_1, g_xxxxyy_0_xyyyyyy_1, g_xxxxyy_0_xyyyyyz_1, g_xxxxyy_0_xyyyyz_1, g_xxxxyy_0_xyyyyzz_1, g_xxxxyy_0_xyyyzz_1, g_xxxxyy_0_xyyyzzz_1, g_xxxxyy_0_xyyzzz_1, g_xxxxyy_0_xyyzzzz_1, g_xxxxyy_0_xyzzzz_1, g_xxxxyy_0_xyzzzzz_1, g_xxxxyy_0_xzzzzz_1, g_xxxxyy_0_xzzzzzz_1, g_xxxxyy_0_yyyyyy_1, g_xxxxyy_0_yyyyyyy_1, g_xxxxyy_0_yyyyyyz_1, g_xxxxyy_0_yyyyyz_1, g_xxxxyy_0_yyyyyzz_1, g_xxxxyy_0_yyyyzz_1, g_xxxxyy_0_yyyyzzz_1, g_xxxxyy_0_yyyzzz_1, g_xxxxyy_0_yyyzzzz_1, g_xxxxyy_0_yyzzzz_1, g_xxxxyy_0_yyzzzzz_1, g_xxxxyy_0_yzzzzz_1, g_xxxxyy_0_yzzzzzz_1, g_xxxxyy_0_zzzzzz_1, g_xxxxyy_0_zzzzzzz_1, g_xxxxyyz_0_xxxxxxx_0, g_xxxxyyz_0_xxxxxxy_0, g_xxxxyyz_0_xxxxxxz_0, g_xxxxyyz_0_xxxxxyy_0, g_xxxxyyz_0_xxxxxyz_0, g_xxxxyyz_0_xxxxxzz_0, g_xxxxyyz_0_xxxxyyy_0, g_xxxxyyz_0_xxxxyyz_0, g_xxxxyyz_0_xxxxyzz_0, g_xxxxyyz_0_xxxxzzz_0, g_xxxxyyz_0_xxxyyyy_0, g_xxxxyyz_0_xxxyyyz_0, g_xxxxyyz_0_xxxyyzz_0, g_xxxxyyz_0_xxxyzzz_0, g_xxxxyyz_0_xxxzzzz_0, g_xxxxyyz_0_xxyyyyy_0, g_xxxxyyz_0_xxyyyyz_0, g_xxxxyyz_0_xxyyyzz_0, g_xxxxyyz_0_xxyyzzz_0, g_xxxxyyz_0_xxyzzzz_0, g_xxxxyyz_0_xxzzzzz_0, g_xxxxyyz_0_xyyyyyy_0, g_xxxxyyz_0_xyyyyyz_0, g_xxxxyyz_0_xyyyyzz_0, g_xxxxyyz_0_xyyyzzz_0, g_xxxxyyz_0_xyyzzzz_0, g_xxxxyyz_0_xyzzzzz_0, g_xxxxyyz_0_xzzzzzz_0, g_xxxxyyz_0_yyyyyyy_0, g_xxxxyyz_0_yyyyyyz_0, g_xxxxyyz_0_yyyyyzz_0, g_xxxxyyz_0_yyyyzzz_0, g_xxxxyyz_0_yyyzzzz_0, g_xxxxyyz_0_yyzzzzz_0, g_xxxxyyz_0_yzzzzzz_0, g_xxxxyyz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyyz_0_xxxxxxx_0[i] = g_xxxxyy_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxxy_0[i] = g_xxxxyy_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxxz_0[i] = g_xxxxyy_0_xxxxxx_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxxz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxyy_0[i] = g_xxxxyy_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxyz_0[i] = g_xxxxyy_0_xxxxxy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxxzz_0[i] = 2.0 * g_xxxxyy_0_xxxxxz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxxzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxyyy_0[i] = g_xxxxyy_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxyyz_0[i] = g_xxxxyy_0_xxxxyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxyzz_0[i] = 2.0 * g_xxxxyy_0_xxxxyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxxzzz_0[i] = 3.0 * g_xxxxyy_0_xxxxzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxxzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyyyy_0[i] = g_xxxxyy_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyyyz_0[i] = g_xxxxyy_0_xxxyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyyzz_0[i] = 2.0 * g_xxxxyy_0_xxxyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxyzzz_0[i] = 3.0 * g_xxxxyy_0_xxxyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxxzzzz_0[i] = 4.0 * g_xxxxyy_0_xxxzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxxzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyyyy_0[i] = g_xxxxyy_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyyyz_0[i] = g_xxxxyy_0_xxyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyyzz_0[i] = 2.0 * g_xxxxyy_0_xxyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyyzzz_0[i] = 3.0 * g_xxxxyy_0_xxyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxyzzzz_0[i] = 4.0 * g_xxxxyy_0_xxyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxyzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xxzzzzz_0[i] = 5.0 * g_xxxxyy_0_xxzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xxzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyyyy_0[i] = g_xxxxyy_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyyyz_0[i] = g_xxxxyy_0_xyyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyyzz_0[i] = 2.0 * g_xxxxyy_0_xyyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyyzzz_0[i] = 3.0 * g_xxxxyy_0_xyyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyyzzzz_0[i] = 4.0 * g_xxxxyy_0_xyyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyyzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xyzzzzz_0[i] = 5.0 * g_xxxxyy_0_xyzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xyzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_xzzzzzz_0[i] = 6.0 * g_xxxxyy_0_xzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_xzzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyyyy_0[i] = g_xxxxyy_0_yyyyyyy_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyyyz_0[i] = g_xxxxyy_0_yyyyyy_1[i] * fi_acd_0 + g_xxxxyy_0_yyyyyyz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyyzz_0[i] = 2.0 * g_xxxxyy_0_yyyyyz_1[i] * fi_acd_0 + g_xxxxyy_0_yyyyyzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyyzzz_0[i] = 3.0 * g_xxxxyy_0_yyyyzz_1[i] * fi_acd_0 + g_xxxxyy_0_yyyyzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyyzzzz_0[i] = 4.0 * g_xxxxyy_0_yyyzzz_1[i] * fi_acd_0 + g_xxxxyy_0_yyyzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yyzzzzz_0[i] = 5.0 * g_xxxxyy_0_yyzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_yyzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_yzzzzzz_0[i] = 6.0 * g_xxxxyy_0_yzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_yzzzzzz_1[i] * wa_z[i];

        g_xxxxyyz_0_zzzzzzz_0[i] = 7.0 * g_xxxxyy_0_zzzzzz_1[i] * fi_acd_0 + g_xxxxyy_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 288-324 components of targeted buffer : KSK

    auto g_xxxxyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 288);

    auto g_xxxxyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 289);

    auto g_xxxxyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 290);

    auto g_xxxxyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 291);

    auto g_xxxxyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 292);

    auto g_xxxxyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 293);

    auto g_xxxxyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 294);

    auto g_xxxxyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 295);

    auto g_xxxxyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 296);

    auto g_xxxxyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 297);

    auto g_xxxxyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 298);

    auto g_xxxxyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 299);

    auto g_xxxxyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 300);

    auto g_xxxxyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 301);

    auto g_xxxxyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 302);

    auto g_xxxxyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 303);

    auto g_xxxxyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 304);

    auto g_xxxxyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 305);

    auto g_xxxxyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 306);

    auto g_xxxxyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 307);

    auto g_xxxxyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 308);

    auto g_xxxxyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 309);

    auto g_xxxxyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 310);

    auto g_xxxxyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 311);

    auto g_xxxxyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 312);

    auto g_xxxxyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 313);

    auto g_xxxxyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 314);

    auto g_xxxxyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 315);

    auto g_xxxxyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 316);

    auto g_xxxxyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 317);

    auto g_xxxxyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 318);

    auto g_xxxxyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 319);

    auto g_xxxxyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 320);

    auto g_xxxxyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 321);

    auto g_xxxxyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 322);

    auto g_xxxxyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 323);

    #pragma omp simd aligned(g_xxxxyzz_0_xxxxxxx_0, g_xxxxyzz_0_xxxxxxy_0, g_xxxxyzz_0_xxxxxxz_0, g_xxxxyzz_0_xxxxxyy_0, g_xxxxyzz_0_xxxxxyz_0, g_xxxxyzz_0_xxxxxzz_0, g_xxxxyzz_0_xxxxyyy_0, g_xxxxyzz_0_xxxxyyz_0, g_xxxxyzz_0_xxxxyzz_0, g_xxxxyzz_0_xxxxzzz_0, g_xxxxyzz_0_xxxyyyy_0, g_xxxxyzz_0_xxxyyyz_0, g_xxxxyzz_0_xxxyyzz_0, g_xxxxyzz_0_xxxyzzz_0, g_xxxxyzz_0_xxxzzzz_0, g_xxxxyzz_0_xxyyyyy_0, g_xxxxyzz_0_xxyyyyz_0, g_xxxxyzz_0_xxyyyzz_0, g_xxxxyzz_0_xxyyzzz_0, g_xxxxyzz_0_xxyzzzz_0, g_xxxxyzz_0_xxzzzzz_0, g_xxxxyzz_0_xyyyyyy_0, g_xxxxyzz_0_xyyyyyz_0, g_xxxxyzz_0_xyyyyzz_0, g_xxxxyzz_0_xyyyzzz_0, g_xxxxyzz_0_xyyzzzz_0, g_xxxxyzz_0_xyzzzzz_0, g_xxxxyzz_0_xzzzzzz_0, g_xxxxyzz_0_yyyyyyy_0, g_xxxxyzz_0_yyyyyyz_0, g_xxxxyzz_0_yyyyyzz_0, g_xxxxyzz_0_yyyyzzz_0, g_xxxxyzz_0_yyyzzzz_0, g_xxxxyzz_0_yyzzzzz_0, g_xxxxyzz_0_yzzzzzz_0, g_xxxxyzz_0_zzzzzzz_0, g_xxxxzz_0_xxxxxx_1, g_xxxxzz_0_xxxxxxx_1, g_xxxxzz_0_xxxxxxy_1, g_xxxxzz_0_xxxxxxz_1, g_xxxxzz_0_xxxxxy_1, g_xxxxzz_0_xxxxxyy_1, g_xxxxzz_0_xxxxxyz_1, g_xxxxzz_0_xxxxxz_1, g_xxxxzz_0_xxxxxzz_1, g_xxxxzz_0_xxxxyy_1, g_xxxxzz_0_xxxxyyy_1, g_xxxxzz_0_xxxxyyz_1, g_xxxxzz_0_xxxxyz_1, g_xxxxzz_0_xxxxyzz_1, g_xxxxzz_0_xxxxzz_1, g_xxxxzz_0_xxxxzzz_1, g_xxxxzz_0_xxxyyy_1, g_xxxxzz_0_xxxyyyy_1, g_xxxxzz_0_xxxyyyz_1, g_xxxxzz_0_xxxyyz_1, g_xxxxzz_0_xxxyyzz_1, g_xxxxzz_0_xxxyzz_1, g_xxxxzz_0_xxxyzzz_1, g_xxxxzz_0_xxxzzz_1, g_xxxxzz_0_xxxzzzz_1, g_xxxxzz_0_xxyyyy_1, g_xxxxzz_0_xxyyyyy_1, g_xxxxzz_0_xxyyyyz_1, g_xxxxzz_0_xxyyyz_1, g_xxxxzz_0_xxyyyzz_1, g_xxxxzz_0_xxyyzz_1, g_xxxxzz_0_xxyyzzz_1, g_xxxxzz_0_xxyzzz_1, g_xxxxzz_0_xxyzzzz_1, g_xxxxzz_0_xxzzzz_1, g_xxxxzz_0_xxzzzzz_1, g_xxxxzz_0_xyyyyy_1, g_xxxxzz_0_xyyyyyy_1, g_xxxxzz_0_xyyyyyz_1, g_xxxxzz_0_xyyyyz_1, g_xxxxzz_0_xyyyyzz_1, g_xxxxzz_0_xyyyzz_1, g_xxxxzz_0_xyyyzzz_1, g_xxxxzz_0_xyyzzz_1, g_xxxxzz_0_xyyzzzz_1, g_xxxxzz_0_xyzzzz_1, g_xxxxzz_0_xyzzzzz_1, g_xxxxzz_0_xzzzzz_1, g_xxxxzz_0_xzzzzzz_1, g_xxxxzz_0_yyyyyy_1, g_xxxxzz_0_yyyyyyy_1, g_xxxxzz_0_yyyyyyz_1, g_xxxxzz_0_yyyyyz_1, g_xxxxzz_0_yyyyyzz_1, g_xxxxzz_0_yyyyzz_1, g_xxxxzz_0_yyyyzzz_1, g_xxxxzz_0_yyyzzz_1, g_xxxxzz_0_yyyzzzz_1, g_xxxxzz_0_yyzzzz_1, g_xxxxzz_0_yyzzzzz_1, g_xxxxzz_0_yzzzzz_1, g_xxxxzz_0_yzzzzzz_1, g_xxxxzz_0_zzzzzz_1, g_xxxxzz_0_zzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyzz_0_xxxxxxx_0[i] = g_xxxxzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxxy_0[i] = g_xxxxzz_0_xxxxxx_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxxy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxxz_0[i] = g_xxxxzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxyy_0[i] = 2.0 * g_xxxxzz_0_xxxxxy_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxyz_0[i] = g_xxxxzz_0_xxxxxz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxxyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxxzz_0[i] = g_xxxxzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxyyy_0[i] = 3.0 * g_xxxxzz_0_xxxxyy_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxyyz_0[i] = 2.0 * g_xxxxzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxyzz_0[i] = g_xxxxzz_0_xxxxzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxxyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxxzzz_0[i] = g_xxxxzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyyyy_0[i] = 4.0 * g_xxxxzz_0_xxxyyy_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyyyz_0[i] = 3.0 * g_xxxxzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyyzz_0[i] = 2.0 * g_xxxxzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxyzzz_0[i] = g_xxxxzz_0_xxxzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxxyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxxzzzz_0[i] = g_xxxxzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyyyy_0[i] = 5.0 * g_xxxxzz_0_xxyyyy_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyyyz_0[i] = 4.0 * g_xxxxzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyyzz_0[i] = 3.0 * g_xxxxzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyyzzz_0[i] = 2.0 * g_xxxxzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxyzzzz_0[i] = g_xxxxzz_0_xxzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xxyzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xxzzzzz_0[i] = g_xxxxzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyyyy_0[i] = 6.0 * g_xxxxzz_0_xyyyyy_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyyyz_0[i] = 5.0 * g_xxxxzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyyzz_0[i] = 4.0 * g_xxxxzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyyzzz_0[i] = 3.0 * g_xxxxzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyyzzzz_0[i] = 2.0 * g_xxxxzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyyzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xyzzzzz_0[i] = g_xxxxzz_0_xzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_xyzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_xzzzzzz_0[i] = g_xxxxzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyyyy_0[i] = 7.0 * g_xxxxzz_0_yyyyyy_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyyyy_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyyyz_0[i] = 6.0 * g_xxxxzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyyyz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyyzz_0[i] = 5.0 * g_xxxxzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyyzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyyzzz_0[i] = 4.0 * g_xxxxzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyyzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyyzzzz_0[i] = 3.0 * g_xxxxzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyyzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yyzzzzz_0[i] = 2.0 * g_xxxxzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yyzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_yzzzzzz_0[i] = g_xxxxzz_0_zzzzzz_1[i] * fi_acd_0 + g_xxxxzz_0_yzzzzzz_1[i] * wa_y[i];

        g_xxxxyzz_0_zzzzzzz_0[i] = g_xxxxzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 324-360 components of targeted buffer : KSK

    auto g_xxxxzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 324);

    auto g_xxxxzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 325);

    auto g_xxxxzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 326);

    auto g_xxxxzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 327);

    auto g_xxxxzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 328);

    auto g_xxxxzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 329);

    auto g_xxxxzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 330);

    auto g_xxxxzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 331);

    auto g_xxxxzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 332);

    auto g_xxxxzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 333);

    auto g_xxxxzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 334);

    auto g_xxxxzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 335);

    auto g_xxxxzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 336);

    auto g_xxxxzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 337);

    auto g_xxxxzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 338);

    auto g_xxxxzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 339);

    auto g_xxxxzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 340);

    auto g_xxxxzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 341);

    auto g_xxxxzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 342);

    auto g_xxxxzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 343);

    auto g_xxxxzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 344);

    auto g_xxxxzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 345);

    auto g_xxxxzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 346);

    auto g_xxxxzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 347);

    auto g_xxxxzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 348);

    auto g_xxxxzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 349);

    auto g_xxxxzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 350);

    auto g_xxxxzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 351);

    auto g_xxxxzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 352);

    auto g_xxxxzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 353);

    auto g_xxxxzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 354);

    auto g_xxxxzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 355);

    auto g_xxxxzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 356);

    auto g_xxxxzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 357);

    auto g_xxxxzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 358);

    auto g_xxxxzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 359);

    #pragma omp simd aligned(g_xxxxz_0_xxxxxxx_0, g_xxxxz_0_xxxxxxx_1, g_xxxxz_0_xxxxxxy_0, g_xxxxz_0_xxxxxxy_1, g_xxxxz_0_xxxxxyy_0, g_xxxxz_0_xxxxxyy_1, g_xxxxz_0_xxxxyyy_0, g_xxxxz_0_xxxxyyy_1, g_xxxxz_0_xxxyyyy_0, g_xxxxz_0_xxxyyyy_1, g_xxxxz_0_xxyyyyy_0, g_xxxxz_0_xxyyyyy_1, g_xxxxz_0_xyyyyyy_0, g_xxxxz_0_xyyyyyy_1, g_xxxxzz_0_xxxxxxx_1, g_xxxxzz_0_xxxxxxy_1, g_xxxxzz_0_xxxxxyy_1, g_xxxxzz_0_xxxxyyy_1, g_xxxxzz_0_xxxyyyy_1, g_xxxxzz_0_xxyyyyy_1, g_xxxxzz_0_xyyyyyy_1, g_xxxxzzz_0_xxxxxxx_0, g_xxxxzzz_0_xxxxxxy_0, g_xxxxzzz_0_xxxxxxz_0, g_xxxxzzz_0_xxxxxyy_0, g_xxxxzzz_0_xxxxxyz_0, g_xxxxzzz_0_xxxxxzz_0, g_xxxxzzz_0_xxxxyyy_0, g_xxxxzzz_0_xxxxyyz_0, g_xxxxzzz_0_xxxxyzz_0, g_xxxxzzz_0_xxxxzzz_0, g_xxxxzzz_0_xxxyyyy_0, g_xxxxzzz_0_xxxyyyz_0, g_xxxxzzz_0_xxxyyzz_0, g_xxxxzzz_0_xxxyzzz_0, g_xxxxzzz_0_xxxzzzz_0, g_xxxxzzz_0_xxyyyyy_0, g_xxxxzzz_0_xxyyyyz_0, g_xxxxzzz_0_xxyyyzz_0, g_xxxxzzz_0_xxyyzzz_0, g_xxxxzzz_0_xxyzzzz_0, g_xxxxzzz_0_xxzzzzz_0, g_xxxxzzz_0_xyyyyyy_0, g_xxxxzzz_0_xyyyyyz_0, g_xxxxzzz_0_xyyyyzz_0, g_xxxxzzz_0_xyyyzzz_0, g_xxxxzzz_0_xyyzzzz_0, g_xxxxzzz_0_xyzzzzz_0, g_xxxxzzz_0_xzzzzzz_0, g_xxxxzzz_0_yyyyyyy_0, g_xxxxzzz_0_yyyyyyz_0, g_xxxxzzz_0_yyyyyzz_0, g_xxxxzzz_0_yyyyzzz_0, g_xxxxzzz_0_yyyzzzz_0, g_xxxxzzz_0_yyzzzzz_0, g_xxxxzzz_0_yzzzzzz_0, g_xxxxzzz_0_zzzzzzz_0, g_xxxzzz_0_xxxxxxz_1, g_xxxzzz_0_xxxxxyz_1, g_xxxzzz_0_xxxxxz_1, g_xxxzzz_0_xxxxxzz_1, g_xxxzzz_0_xxxxyyz_1, g_xxxzzz_0_xxxxyz_1, g_xxxzzz_0_xxxxyzz_1, g_xxxzzz_0_xxxxzz_1, g_xxxzzz_0_xxxxzzz_1, g_xxxzzz_0_xxxyyyz_1, g_xxxzzz_0_xxxyyz_1, g_xxxzzz_0_xxxyyzz_1, g_xxxzzz_0_xxxyzz_1, g_xxxzzz_0_xxxyzzz_1, g_xxxzzz_0_xxxzzz_1, g_xxxzzz_0_xxxzzzz_1, g_xxxzzz_0_xxyyyyz_1, g_xxxzzz_0_xxyyyz_1, g_xxxzzz_0_xxyyyzz_1, g_xxxzzz_0_xxyyzz_1, g_xxxzzz_0_xxyyzzz_1, g_xxxzzz_0_xxyzzz_1, g_xxxzzz_0_xxyzzzz_1, g_xxxzzz_0_xxzzzz_1, g_xxxzzz_0_xxzzzzz_1, g_xxxzzz_0_xyyyyyz_1, g_xxxzzz_0_xyyyyz_1, g_xxxzzz_0_xyyyyzz_1, g_xxxzzz_0_xyyyzz_1, g_xxxzzz_0_xyyyzzz_1, g_xxxzzz_0_xyyzzz_1, g_xxxzzz_0_xyyzzzz_1, g_xxxzzz_0_xyzzzz_1, g_xxxzzz_0_xyzzzzz_1, g_xxxzzz_0_xzzzzz_1, g_xxxzzz_0_xzzzzzz_1, g_xxxzzz_0_yyyyyyy_1, g_xxxzzz_0_yyyyyyz_1, g_xxxzzz_0_yyyyyz_1, g_xxxzzz_0_yyyyyzz_1, g_xxxzzz_0_yyyyzz_1, g_xxxzzz_0_yyyyzzz_1, g_xxxzzz_0_yyyzzz_1, g_xxxzzz_0_yyyzzzz_1, g_xxxzzz_0_yyzzzz_1, g_xxxzzz_0_yyzzzzz_1, g_xxxzzz_0_yzzzzz_1, g_xxxzzz_0_yzzzzzz_1, g_xxxzzz_0_zzzzzz_1, g_xxxzzz_0_zzzzzzz_1, g_xxzzz_0_xxxxxxz_0, g_xxzzz_0_xxxxxxz_1, g_xxzzz_0_xxxxxyz_0, g_xxzzz_0_xxxxxyz_1, g_xxzzz_0_xxxxxzz_0, g_xxzzz_0_xxxxxzz_1, g_xxzzz_0_xxxxyyz_0, g_xxzzz_0_xxxxyyz_1, g_xxzzz_0_xxxxyzz_0, g_xxzzz_0_xxxxyzz_1, g_xxzzz_0_xxxxzzz_0, g_xxzzz_0_xxxxzzz_1, g_xxzzz_0_xxxyyyz_0, g_xxzzz_0_xxxyyyz_1, g_xxzzz_0_xxxyyzz_0, g_xxzzz_0_xxxyyzz_1, g_xxzzz_0_xxxyzzz_0, g_xxzzz_0_xxxyzzz_1, g_xxzzz_0_xxxzzzz_0, g_xxzzz_0_xxxzzzz_1, g_xxzzz_0_xxyyyyz_0, g_xxzzz_0_xxyyyyz_1, g_xxzzz_0_xxyyyzz_0, g_xxzzz_0_xxyyyzz_1, g_xxzzz_0_xxyyzzz_0, g_xxzzz_0_xxyyzzz_1, g_xxzzz_0_xxyzzzz_0, g_xxzzz_0_xxyzzzz_1, g_xxzzz_0_xxzzzzz_0, g_xxzzz_0_xxzzzzz_1, g_xxzzz_0_xyyyyyz_0, g_xxzzz_0_xyyyyyz_1, g_xxzzz_0_xyyyyzz_0, g_xxzzz_0_xyyyyzz_1, g_xxzzz_0_xyyyzzz_0, g_xxzzz_0_xyyyzzz_1, g_xxzzz_0_xyyzzzz_0, g_xxzzz_0_xyyzzzz_1, g_xxzzz_0_xyzzzzz_0, g_xxzzz_0_xyzzzzz_1, g_xxzzz_0_xzzzzzz_0, g_xxzzz_0_xzzzzzz_1, g_xxzzz_0_yyyyyyy_0, g_xxzzz_0_yyyyyyy_1, g_xxzzz_0_yyyyyyz_0, g_xxzzz_0_yyyyyyz_1, g_xxzzz_0_yyyyyzz_0, g_xxzzz_0_yyyyyzz_1, g_xxzzz_0_yyyyzzz_0, g_xxzzz_0_yyyyzzz_1, g_xxzzz_0_yyyzzzz_0, g_xxzzz_0_yyyzzzz_1, g_xxzzz_0_yyzzzzz_0, g_xxzzz_0_yyzzzzz_1, g_xxzzz_0_yzzzzzz_0, g_xxzzz_0_yzzzzzz_1, g_xxzzz_0_zzzzzzz_0, g_xxzzz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxzzz_0_xxxxxxx_0[i] = 2.0 * g_xxxxz_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxxxx_1[i] * fz_be_0 + g_xxxxzz_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxxxy_0[i] = 2.0 * g_xxxxz_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxxxy_1[i] * fz_be_0 + g_xxxxzz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxxxz_0[i] = 3.0 * g_xxzzz_0_xxxxxxz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xxxzzz_0_xxxxxz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxxyy_0[i] = 2.0 * g_xxxxz_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxxyy_1[i] * fz_be_0 + g_xxxxzz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxxyz_0[i] = 3.0 * g_xxzzz_0_xxxxxyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxxzzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxxzz_0[i] = 3.0 * g_xxzzz_0_xxxxxzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xxxzzz_0_xxxxzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxyyy_0[i] = 2.0 * g_xxxxz_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxxyyy_1[i] * fz_be_0 + g_xxxxzz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxxyyz_0[i] = 3.0 * g_xxzzz_0_xxxxyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxyzz_0[i] = 3.0 * g_xxzzz_0_xxxxyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxxzzz_0[i] = 3.0 * g_xxzzz_0_xxxxzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xxxzzz_0_xxxzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxyyyy_0[i] = 2.0 * g_xxxxz_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxxyyyy_1[i] * fz_be_0 + g_xxxxzz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxxyyyz_0[i] = 3.0 * g_xxzzz_0_xxxyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxyyzz_0[i] = 3.0 * g_xxzzz_0_xxxyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxyzzz_0[i] = 3.0 * g_xxzzz_0_xxxyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxxzzzz_0[i] = 3.0 * g_xxzzz_0_xxxzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xxxzzz_0_xxzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyyyyy_0[i] = 2.0 * g_xxxxz_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xxyyyyy_1[i] * fz_be_0 + g_xxxxzz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xxyyyyz_0[i] = 3.0 * g_xxzzz_0_xxyyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyyyzz_0[i] = 3.0 * g_xxzzz_0_xxyyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyyzzz_0[i] = 3.0 * g_xxzzz_0_xxyyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxyzzzz_0[i] = 3.0 * g_xxzzz_0_xxyzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xxzzzzz_0[i] = 3.0 * g_xxzzz_0_xxzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_0_xzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyyyyy_0[i] = 2.0 * g_xxxxz_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_xxxxz_0_xyyyyyy_1[i] * fz_be_0 + g_xxxxzz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxxzzz_0_xyyyyyz_0[i] = 3.0 * g_xxzzz_0_xyyyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyyyyz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyyyzz_0[i] = 3.0 * g_xxzzz_0_xyyyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyyyzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyyzzz_0[i] = 3.0 * g_xxzzz_0_xyyyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyyzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyyzzzz_0[i] = 3.0 * g_xxzzz_0_xyyzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyyzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xyzzzzz_0[i] = 3.0 * g_xxzzz_0_xyzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xyzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_xzzzzzz_0[i] = 3.0 * g_xxzzz_0_xzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_xzzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_zzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyyyy_0[i] = 3.0 * g_xxzzz_0_yyyyyyy_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyyyy_1[i] * fz_be_0 + g_xxxzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyyyz_0[i] = 3.0 * g_xxzzz_0_yyyyyyz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyyyz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyyzz_0[i] = 3.0 * g_xxzzz_0_yyyyyzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyyzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyyzzz_0[i] = 3.0 * g_xxzzz_0_yyyyzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyyzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyyzzzz_0[i] = 3.0 * g_xxzzz_0_yyyzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyyzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yyzzzzz_0[i] = 3.0 * g_xxzzz_0_yyzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yyzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_yzzzzzz_0[i] = 3.0 * g_xxzzz_0_yzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_yzzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxxzzz_0_zzzzzzz_0[i] = 3.0 * g_xxzzz_0_zzzzzzz_0[i] * fbe_0 - 3.0 * g_xxzzz_0_zzzzzzz_1[i] * fz_be_0 + g_xxxzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 360-396 components of targeted buffer : KSK

    auto g_xxxyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 360);

    auto g_xxxyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 361);

    auto g_xxxyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 362);

    auto g_xxxyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 363);

    auto g_xxxyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 364);

    auto g_xxxyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 365);

    auto g_xxxyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 366);

    auto g_xxxyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 367);

    auto g_xxxyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 368);

    auto g_xxxyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 369);

    auto g_xxxyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 370);

    auto g_xxxyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 371);

    auto g_xxxyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 372);

    auto g_xxxyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 373);

    auto g_xxxyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 374);

    auto g_xxxyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 375);

    auto g_xxxyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 376);

    auto g_xxxyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 377);

    auto g_xxxyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 378);

    auto g_xxxyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 379);

    auto g_xxxyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 380);

    auto g_xxxyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 381);

    auto g_xxxyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 382);

    auto g_xxxyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 383);

    auto g_xxxyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 384);

    auto g_xxxyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 385);

    auto g_xxxyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 386);

    auto g_xxxyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 387);

    auto g_xxxyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 388);

    auto g_xxxyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 389);

    auto g_xxxyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 390);

    auto g_xxxyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 391);

    auto g_xxxyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 392);

    auto g_xxxyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 393);

    auto g_xxxyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 394);

    auto g_xxxyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 395);

    #pragma omp simd aligned(g_xxxyy_0_xxxxxxx_0, g_xxxyy_0_xxxxxxx_1, g_xxxyy_0_xxxxxxz_0, g_xxxyy_0_xxxxxxz_1, g_xxxyy_0_xxxxxzz_0, g_xxxyy_0_xxxxxzz_1, g_xxxyy_0_xxxxzzz_0, g_xxxyy_0_xxxxzzz_1, g_xxxyy_0_xxxzzzz_0, g_xxxyy_0_xxxzzzz_1, g_xxxyy_0_xxzzzzz_0, g_xxxyy_0_xxzzzzz_1, g_xxxyy_0_xzzzzzz_0, g_xxxyy_0_xzzzzzz_1, g_xxxyyy_0_xxxxxxx_1, g_xxxyyy_0_xxxxxxz_1, g_xxxyyy_0_xxxxxzz_1, g_xxxyyy_0_xxxxzzz_1, g_xxxyyy_0_xxxzzzz_1, g_xxxyyy_0_xxzzzzz_1, g_xxxyyy_0_xzzzzzz_1, g_xxxyyyy_0_xxxxxxx_0, g_xxxyyyy_0_xxxxxxy_0, g_xxxyyyy_0_xxxxxxz_0, g_xxxyyyy_0_xxxxxyy_0, g_xxxyyyy_0_xxxxxyz_0, g_xxxyyyy_0_xxxxxzz_0, g_xxxyyyy_0_xxxxyyy_0, g_xxxyyyy_0_xxxxyyz_0, g_xxxyyyy_0_xxxxyzz_0, g_xxxyyyy_0_xxxxzzz_0, g_xxxyyyy_0_xxxyyyy_0, g_xxxyyyy_0_xxxyyyz_0, g_xxxyyyy_0_xxxyyzz_0, g_xxxyyyy_0_xxxyzzz_0, g_xxxyyyy_0_xxxzzzz_0, g_xxxyyyy_0_xxyyyyy_0, g_xxxyyyy_0_xxyyyyz_0, g_xxxyyyy_0_xxyyyzz_0, g_xxxyyyy_0_xxyyzzz_0, g_xxxyyyy_0_xxyzzzz_0, g_xxxyyyy_0_xxzzzzz_0, g_xxxyyyy_0_xyyyyyy_0, g_xxxyyyy_0_xyyyyyz_0, g_xxxyyyy_0_xyyyyzz_0, g_xxxyyyy_0_xyyyzzz_0, g_xxxyyyy_0_xyyzzzz_0, g_xxxyyyy_0_xyzzzzz_0, g_xxxyyyy_0_xzzzzzz_0, g_xxxyyyy_0_yyyyyyy_0, g_xxxyyyy_0_yyyyyyz_0, g_xxxyyyy_0_yyyyyzz_0, g_xxxyyyy_0_yyyyzzz_0, g_xxxyyyy_0_yyyzzzz_0, g_xxxyyyy_0_yyzzzzz_0, g_xxxyyyy_0_yzzzzzz_0, g_xxxyyyy_0_zzzzzzz_0, g_xxyyyy_0_xxxxxxy_1, g_xxyyyy_0_xxxxxy_1, g_xxyyyy_0_xxxxxyy_1, g_xxyyyy_0_xxxxxyz_1, g_xxyyyy_0_xxxxyy_1, g_xxyyyy_0_xxxxyyy_1, g_xxyyyy_0_xxxxyyz_1, g_xxyyyy_0_xxxxyz_1, g_xxyyyy_0_xxxxyzz_1, g_xxyyyy_0_xxxyyy_1, g_xxyyyy_0_xxxyyyy_1, g_xxyyyy_0_xxxyyyz_1, g_xxyyyy_0_xxxyyz_1, g_xxyyyy_0_xxxyyzz_1, g_xxyyyy_0_xxxyzz_1, g_xxyyyy_0_xxxyzzz_1, g_xxyyyy_0_xxyyyy_1, g_xxyyyy_0_xxyyyyy_1, g_xxyyyy_0_xxyyyyz_1, g_xxyyyy_0_xxyyyz_1, g_xxyyyy_0_xxyyyzz_1, g_xxyyyy_0_xxyyzz_1, g_xxyyyy_0_xxyyzzz_1, g_xxyyyy_0_xxyzzz_1, g_xxyyyy_0_xxyzzzz_1, g_xxyyyy_0_xyyyyy_1, g_xxyyyy_0_xyyyyyy_1, g_xxyyyy_0_xyyyyyz_1, g_xxyyyy_0_xyyyyz_1, g_xxyyyy_0_xyyyyzz_1, g_xxyyyy_0_xyyyzz_1, g_xxyyyy_0_xyyyzzz_1, g_xxyyyy_0_xyyzzz_1, g_xxyyyy_0_xyyzzzz_1, g_xxyyyy_0_xyzzzz_1, g_xxyyyy_0_xyzzzzz_1, g_xxyyyy_0_yyyyyy_1, g_xxyyyy_0_yyyyyyy_1, g_xxyyyy_0_yyyyyyz_1, g_xxyyyy_0_yyyyyz_1, g_xxyyyy_0_yyyyyzz_1, g_xxyyyy_0_yyyyzz_1, g_xxyyyy_0_yyyyzzz_1, g_xxyyyy_0_yyyzzz_1, g_xxyyyy_0_yyyzzzz_1, g_xxyyyy_0_yyzzzz_1, g_xxyyyy_0_yyzzzzz_1, g_xxyyyy_0_yzzzzz_1, g_xxyyyy_0_yzzzzzz_1, g_xxyyyy_0_zzzzzzz_1, g_xyyyy_0_xxxxxxy_0, g_xyyyy_0_xxxxxxy_1, g_xyyyy_0_xxxxxyy_0, g_xyyyy_0_xxxxxyy_1, g_xyyyy_0_xxxxxyz_0, g_xyyyy_0_xxxxxyz_1, g_xyyyy_0_xxxxyyy_0, g_xyyyy_0_xxxxyyy_1, g_xyyyy_0_xxxxyyz_0, g_xyyyy_0_xxxxyyz_1, g_xyyyy_0_xxxxyzz_0, g_xyyyy_0_xxxxyzz_1, g_xyyyy_0_xxxyyyy_0, g_xyyyy_0_xxxyyyy_1, g_xyyyy_0_xxxyyyz_0, g_xyyyy_0_xxxyyyz_1, g_xyyyy_0_xxxyyzz_0, g_xyyyy_0_xxxyyzz_1, g_xyyyy_0_xxxyzzz_0, g_xyyyy_0_xxxyzzz_1, g_xyyyy_0_xxyyyyy_0, g_xyyyy_0_xxyyyyy_1, g_xyyyy_0_xxyyyyz_0, g_xyyyy_0_xxyyyyz_1, g_xyyyy_0_xxyyyzz_0, g_xyyyy_0_xxyyyzz_1, g_xyyyy_0_xxyyzzz_0, g_xyyyy_0_xxyyzzz_1, g_xyyyy_0_xxyzzzz_0, g_xyyyy_0_xxyzzzz_1, g_xyyyy_0_xyyyyyy_0, g_xyyyy_0_xyyyyyy_1, g_xyyyy_0_xyyyyyz_0, g_xyyyy_0_xyyyyyz_1, g_xyyyy_0_xyyyyzz_0, g_xyyyy_0_xyyyyzz_1, g_xyyyy_0_xyyyzzz_0, g_xyyyy_0_xyyyzzz_1, g_xyyyy_0_xyyzzzz_0, g_xyyyy_0_xyyzzzz_1, g_xyyyy_0_xyzzzzz_0, g_xyyyy_0_xyzzzzz_1, g_xyyyy_0_yyyyyyy_0, g_xyyyy_0_yyyyyyy_1, g_xyyyy_0_yyyyyyz_0, g_xyyyy_0_yyyyyyz_1, g_xyyyy_0_yyyyyzz_0, g_xyyyy_0_yyyyyzz_1, g_xyyyy_0_yyyyzzz_0, g_xyyyy_0_yyyyzzz_1, g_xyyyy_0_yyyzzzz_0, g_xyyyy_0_yyyzzzz_1, g_xyyyy_0_yyzzzzz_0, g_xyyyy_0_yyzzzzz_1, g_xyyyy_0_yzzzzzz_0, g_xyyyy_0_yzzzzzz_1, g_xyyyy_0_zzzzzzz_0, g_xyyyy_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyyy_0_xxxxxxx_0[i] = 3.0 * g_xxxyy_0_xxxxxxx_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxxxx_1[i] * fz_be_0 + g_xxxyyy_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxxxxy_0[i] = 2.0 * g_xyyyy_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xxyyyy_0_xxxxxy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxxxz_0[i] = 3.0 * g_xxxyy_0_xxxxxxz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxxxz_1[i] * fz_be_0 + g_xxxyyy_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxxxyy_0[i] = 2.0 * g_xyyyy_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xxyyyy_0_xxxxyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxxyz_0[i] = 2.0 * g_xyyyy_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxyyyy_0_xxxxyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxxzz_0[i] = 3.0 * g_xxxyy_0_xxxxxzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxxzz_1[i] * fz_be_0 + g_xxxyyy_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxxyyy_0[i] = 2.0 * g_xyyyy_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xxyyyy_0_xxxyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxyyz_0[i] = 2.0 * g_xyyyy_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxyyyy_0_xxxyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxyzz_0[i] = 2.0 * g_xyyyy_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxyyyy_0_xxxyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxxzzz_0[i] = 3.0 * g_xxxyy_0_xxxxzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxxzzz_1[i] * fz_be_0 + g_xxxyyy_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxxyyyy_0[i] = 2.0 * g_xyyyy_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxyyyz_0[i] = 2.0 * g_xyyyy_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxyyzz_0[i] = 2.0 * g_xyyyy_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxyzzz_0[i] = 2.0 * g_xyyyy_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxyyyy_0_xxyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxxzzzz_0[i] = 3.0 * g_xxxyy_0_xxxzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxxzzzz_1[i] * fz_be_0 + g_xxxyyy_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xxyyyyy_0[i] = 2.0 * g_xyyyy_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyyyyz_0[i] = 2.0 * g_xyyyy_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyyyzz_0[i] = 2.0 * g_xyyyy_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyyzzz_0[i] = 2.0 * g_xyyyy_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxyzzzz_0[i] = 2.0 * g_xyyyy_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxyyyy_0_xyzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xxzzzzz_0[i] = 3.0 * g_xxxyy_0_xxzzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xxzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_xyyyyyy_0[i] = 2.0 * g_xyyyy_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyyyy_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyyyyz_0[i] = 2.0 * g_xyyyy_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyyyz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyyyzz_0[i] = 2.0 * g_xyyyy_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyyzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyyzzz_0[i] = 2.0 * g_xyyyy_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyyzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyyzzzz_0[i] = 2.0 * g_xyyyy_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyyzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xyzzzzz_0[i] = 2.0 * g_xyyyy_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_xyzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_xzzzzzz_0[i] = 3.0 * g_xxxyy_0_xzzzzzz_0[i] * fbe_0 - 3.0 * g_xxxyy_0_xzzzzzz_1[i] * fz_be_0 + g_xxxyyy_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxyyyy_0_yyyyyyy_0[i] = 2.0 * g_xyyyy_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyyyy_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyyyyz_0[i] = 2.0 * g_xyyyy_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyyyz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyyyzz_0[i] = 2.0 * g_xyyyy_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyyzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyyzzz_0[i] = 2.0 * g_xyyyy_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyyzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyyzzzz_0[i] = 2.0 * g_xyyyy_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyyzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yyzzzzz_0[i] = 2.0 * g_xyyyy_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yyzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_yzzzzzz_0[i] = 2.0 * g_xyyyy_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_yzzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxyyyy_0_zzzzzzz_0[i] = 2.0 * g_xyyyy_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyyy_0_zzzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 396-432 components of targeted buffer : KSK

    auto g_xxxyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 396);

    auto g_xxxyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 397);

    auto g_xxxyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 398);

    auto g_xxxyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 399);

    auto g_xxxyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 400);

    auto g_xxxyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 401);

    auto g_xxxyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 402);

    auto g_xxxyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 403);

    auto g_xxxyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 404);

    auto g_xxxyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 405);

    auto g_xxxyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 406);

    auto g_xxxyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 407);

    auto g_xxxyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 408);

    auto g_xxxyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 409);

    auto g_xxxyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 410);

    auto g_xxxyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 411);

    auto g_xxxyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 412);

    auto g_xxxyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 413);

    auto g_xxxyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 414);

    auto g_xxxyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 415);

    auto g_xxxyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 416);

    auto g_xxxyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 417);

    auto g_xxxyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 418);

    auto g_xxxyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 419);

    auto g_xxxyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 420);

    auto g_xxxyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 421);

    auto g_xxxyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 422);

    auto g_xxxyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 423);

    auto g_xxxyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 424);

    auto g_xxxyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 425);

    auto g_xxxyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 426);

    auto g_xxxyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 427);

    auto g_xxxyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 428);

    auto g_xxxyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 429);

    auto g_xxxyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 430);

    auto g_xxxyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 431);

    #pragma omp simd aligned(g_xxxyyy_0_xxxxxx_1, g_xxxyyy_0_xxxxxxx_1, g_xxxyyy_0_xxxxxxy_1, g_xxxyyy_0_xxxxxxz_1, g_xxxyyy_0_xxxxxy_1, g_xxxyyy_0_xxxxxyy_1, g_xxxyyy_0_xxxxxyz_1, g_xxxyyy_0_xxxxxz_1, g_xxxyyy_0_xxxxxzz_1, g_xxxyyy_0_xxxxyy_1, g_xxxyyy_0_xxxxyyy_1, g_xxxyyy_0_xxxxyyz_1, g_xxxyyy_0_xxxxyz_1, g_xxxyyy_0_xxxxyzz_1, g_xxxyyy_0_xxxxzz_1, g_xxxyyy_0_xxxxzzz_1, g_xxxyyy_0_xxxyyy_1, g_xxxyyy_0_xxxyyyy_1, g_xxxyyy_0_xxxyyyz_1, g_xxxyyy_0_xxxyyz_1, g_xxxyyy_0_xxxyyzz_1, g_xxxyyy_0_xxxyzz_1, g_xxxyyy_0_xxxyzzz_1, g_xxxyyy_0_xxxzzz_1, g_xxxyyy_0_xxxzzzz_1, g_xxxyyy_0_xxyyyy_1, g_xxxyyy_0_xxyyyyy_1, g_xxxyyy_0_xxyyyyz_1, g_xxxyyy_0_xxyyyz_1, g_xxxyyy_0_xxyyyzz_1, g_xxxyyy_0_xxyyzz_1, g_xxxyyy_0_xxyyzzz_1, g_xxxyyy_0_xxyzzz_1, g_xxxyyy_0_xxyzzzz_1, g_xxxyyy_0_xxzzzz_1, g_xxxyyy_0_xxzzzzz_1, g_xxxyyy_0_xyyyyy_1, g_xxxyyy_0_xyyyyyy_1, g_xxxyyy_0_xyyyyyz_1, g_xxxyyy_0_xyyyyz_1, g_xxxyyy_0_xyyyyzz_1, g_xxxyyy_0_xyyyzz_1, g_xxxyyy_0_xyyyzzz_1, g_xxxyyy_0_xyyzzz_1, g_xxxyyy_0_xyyzzzz_1, g_xxxyyy_0_xyzzzz_1, g_xxxyyy_0_xyzzzzz_1, g_xxxyyy_0_xzzzzz_1, g_xxxyyy_0_xzzzzzz_1, g_xxxyyy_0_yyyyyy_1, g_xxxyyy_0_yyyyyyy_1, g_xxxyyy_0_yyyyyyz_1, g_xxxyyy_0_yyyyyz_1, g_xxxyyy_0_yyyyyzz_1, g_xxxyyy_0_yyyyzz_1, g_xxxyyy_0_yyyyzzz_1, g_xxxyyy_0_yyyzzz_1, g_xxxyyy_0_yyyzzzz_1, g_xxxyyy_0_yyzzzz_1, g_xxxyyy_0_yyzzzzz_1, g_xxxyyy_0_yzzzzz_1, g_xxxyyy_0_yzzzzzz_1, g_xxxyyy_0_zzzzzz_1, g_xxxyyy_0_zzzzzzz_1, g_xxxyyyz_0_xxxxxxx_0, g_xxxyyyz_0_xxxxxxy_0, g_xxxyyyz_0_xxxxxxz_0, g_xxxyyyz_0_xxxxxyy_0, g_xxxyyyz_0_xxxxxyz_0, g_xxxyyyz_0_xxxxxzz_0, g_xxxyyyz_0_xxxxyyy_0, g_xxxyyyz_0_xxxxyyz_0, g_xxxyyyz_0_xxxxyzz_0, g_xxxyyyz_0_xxxxzzz_0, g_xxxyyyz_0_xxxyyyy_0, g_xxxyyyz_0_xxxyyyz_0, g_xxxyyyz_0_xxxyyzz_0, g_xxxyyyz_0_xxxyzzz_0, g_xxxyyyz_0_xxxzzzz_0, g_xxxyyyz_0_xxyyyyy_0, g_xxxyyyz_0_xxyyyyz_0, g_xxxyyyz_0_xxyyyzz_0, g_xxxyyyz_0_xxyyzzz_0, g_xxxyyyz_0_xxyzzzz_0, g_xxxyyyz_0_xxzzzzz_0, g_xxxyyyz_0_xyyyyyy_0, g_xxxyyyz_0_xyyyyyz_0, g_xxxyyyz_0_xyyyyzz_0, g_xxxyyyz_0_xyyyzzz_0, g_xxxyyyz_0_xyyzzzz_0, g_xxxyyyz_0_xyzzzzz_0, g_xxxyyyz_0_xzzzzzz_0, g_xxxyyyz_0_yyyyyyy_0, g_xxxyyyz_0_yyyyyyz_0, g_xxxyyyz_0_yyyyyzz_0, g_xxxyyyz_0_yyyyzzz_0, g_xxxyyyz_0_yyyzzzz_0, g_xxxyyyz_0_yyzzzzz_0, g_xxxyyyz_0_yzzzzzz_0, g_xxxyyyz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyyz_0_xxxxxxx_0[i] = g_xxxyyy_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxxy_0[i] = g_xxxyyy_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxxz_0[i] = g_xxxyyy_0_xxxxxx_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxxz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxyy_0[i] = g_xxxyyy_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxyz_0[i] = g_xxxyyy_0_xxxxxy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxxzz_0[i] = 2.0 * g_xxxyyy_0_xxxxxz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxxzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxyyy_0[i] = g_xxxyyy_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxyyz_0[i] = g_xxxyyy_0_xxxxyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxyzz_0[i] = 2.0 * g_xxxyyy_0_xxxxyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxxzzz_0[i] = 3.0 * g_xxxyyy_0_xxxxzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxxzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyyyy_0[i] = g_xxxyyy_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyyyz_0[i] = g_xxxyyy_0_xxxyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyyzz_0[i] = 2.0 * g_xxxyyy_0_xxxyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxyzzz_0[i] = 3.0 * g_xxxyyy_0_xxxyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxxzzzz_0[i] = 4.0 * g_xxxyyy_0_xxxzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxxzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyyyy_0[i] = g_xxxyyy_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyyyz_0[i] = g_xxxyyy_0_xxyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyyzz_0[i] = 2.0 * g_xxxyyy_0_xxyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyyzzz_0[i] = 3.0 * g_xxxyyy_0_xxyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxyzzzz_0[i] = 4.0 * g_xxxyyy_0_xxyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxyzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xxzzzzz_0[i] = 5.0 * g_xxxyyy_0_xxzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xxzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyyyy_0[i] = g_xxxyyy_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyyyz_0[i] = g_xxxyyy_0_xyyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyyzz_0[i] = 2.0 * g_xxxyyy_0_xyyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyyzzz_0[i] = 3.0 * g_xxxyyy_0_xyyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyyzzzz_0[i] = 4.0 * g_xxxyyy_0_xyyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyyzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xyzzzzz_0[i] = 5.0 * g_xxxyyy_0_xyzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xyzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_xzzzzzz_0[i] = 6.0 * g_xxxyyy_0_xzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_xzzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyyyy_0[i] = g_xxxyyy_0_yyyyyyy_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyyyz_0[i] = g_xxxyyy_0_yyyyyy_1[i] * fi_acd_0 + g_xxxyyy_0_yyyyyyz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyyzz_0[i] = 2.0 * g_xxxyyy_0_yyyyyz_1[i] * fi_acd_0 + g_xxxyyy_0_yyyyyzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyyzzz_0[i] = 3.0 * g_xxxyyy_0_yyyyzz_1[i] * fi_acd_0 + g_xxxyyy_0_yyyyzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyyzzzz_0[i] = 4.0 * g_xxxyyy_0_yyyzzz_1[i] * fi_acd_0 + g_xxxyyy_0_yyyzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yyzzzzz_0[i] = 5.0 * g_xxxyyy_0_yyzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_yyzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_yzzzzzz_0[i] = 6.0 * g_xxxyyy_0_yzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_yzzzzzz_1[i] * wa_z[i];

        g_xxxyyyz_0_zzzzzzz_0[i] = 7.0 * g_xxxyyy_0_zzzzzz_1[i] * fi_acd_0 + g_xxxyyy_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 432-468 components of targeted buffer : KSK

    auto g_xxxyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 432);

    auto g_xxxyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 433);

    auto g_xxxyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 434);

    auto g_xxxyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 435);

    auto g_xxxyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 436);

    auto g_xxxyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 437);

    auto g_xxxyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 438);

    auto g_xxxyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 439);

    auto g_xxxyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 440);

    auto g_xxxyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 441);

    auto g_xxxyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 442);

    auto g_xxxyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 443);

    auto g_xxxyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 444);

    auto g_xxxyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 445);

    auto g_xxxyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 446);

    auto g_xxxyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 447);

    auto g_xxxyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 448);

    auto g_xxxyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 449);

    auto g_xxxyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 450);

    auto g_xxxyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 451);

    auto g_xxxyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 452);

    auto g_xxxyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 453);

    auto g_xxxyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 454);

    auto g_xxxyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 455);

    auto g_xxxyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 456);

    auto g_xxxyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 457);

    auto g_xxxyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 458);

    auto g_xxxyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 459);

    auto g_xxxyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 460);

    auto g_xxxyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 461);

    auto g_xxxyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 462);

    auto g_xxxyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 463);

    auto g_xxxyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 464);

    auto g_xxxyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 465);

    auto g_xxxyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 466);

    auto g_xxxyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 467);

    #pragma omp simd aligned(g_xxxyy_0_xxxxxxy_0, g_xxxyy_0_xxxxxxy_1, g_xxxyy_0_xxxxxyy_0, g_xxxyy_0_xxxxxyy_1, g_xxxyy_0_xxxxyyy_0, g_xxxyy_0_xxxxyyy_1, g_xxxyy_0_xxxyyyy_0, g_xxxyy_0_xxxyyyy_1, g_xxxyy_0_xxyyyyy_0, g_xxxyy_0_xxyyyyy_1, g_xxxyy_0_xyyyyyy_0, g_xxxyy_0_xyyyyyy_1, g_xxxyyz_0_xxxxxxy_1, g_xxxyyz_0_xxxxxyy_1, g_xxxyyz_0_xxxxyyy_1, g_xxxyyz_0_xxxyyyy_1, g_xxxyyz_0_xxyyyyy_1, g_xxxyyz_0_xyyyyyy_1, g_xxxyyzz_0_xxxxxxx_0, g_xxxyyzz_0_xxxxxxy_0, g_xxxyyzz_0_xxxxxxz_0, g_xxxyyzz_0_xxxxxyy_0, g_xxxyyzz_0_xxxxxyz_0, g_xxxyyzz_0_xxxxxzz_0, g_xxxyyzz_0_xxxxyyy_0, g_xxxyyzz_0_xxxxyyz_0, g_xxxyyzz_0_xxxxyzz_0, g_xxxyyzz_0_xxxxzzz_0, g_xxxyyzz_0_xxxyyyy_0, g_xxxyyzz_0_xxxyyyz_0, g_xxxyyzz_0_xxxyyzz_0, g_xxxyyzz_0_xxxyzzz_0, g_xxxyyzz_0_xxxzzzz_0, g_xxxyyzz_0_xxyyyyy_0, g_xxxyyzz_0_xxyyyyz_0, g_xxxyyzz_0_xxyyyzz_0, g_xxxyyzz_0_xxyyzzz_0, g_xxxyyzz_0_xxyzzzz_0, g_xxxyyzz_0_xxzzzzz_0, g_xxxyyzz_0_xyyyyyy_0, g_xxxyyzz_0_xyyyyyz_0, g_xxxyyzz_0_xyyyyzz_0, g_xxxyyzz_0_xyyyzzz_0, g_xxxyyzz_0_xyyzzzz_0, g_xxxyyzz_0_xyzzzzz_0, g_xxxyyzz_0_xzzzzzz_0, g_xxxyyzz_0_yyyyyyy_0, g_xxxyyzz_0_yyyyyyz_0, g_xxxyyzz_0_yyyyyzz_0, g_xxxyyzz_0_yyyyzzz_0, g_xxxyyzz_0_yyyzzzz_0, g_xxxyyzz_0_yyzzzzz_0, g_xxxyyzz_0_yzzzzzz_0, g_xxxyyzz_0_zzzzzzz_0, g_xxxyzz_0_xxxxxxx_1, g_xxxyzz_0_xxxxxxz_1, g_xxxyzz_0_xxxxxzz_1, g_xxxyzz_0_xxxxzzz_1, g_xxxyzz_0_xxxzzzz_1, g_xxxyzz_0_xxzzzzz_1, g_xxxyzz_0_xzzzzzz_1, g_xxxzz_0_xxxxxxx_0, g_xxxzz_0_xxxxxxx_1, g_xxxzz_0_xxxxxxz_0, g_xxxzz_0_xxxxxxz_1, g_xxxzz_0_xxxxxzz_0, g_xxxzz_0_xxxxxzz_1, g_xxxzz_0_xxxxzzz_0, g_xxxzz_0_xxxxzzz_1, g_xxxzz_0_xxxzzzz_0, g_xxxzz_0_xxxzzzz_1, g_xxxzz_0_xxzzzzz_0, g_xxxzz_0_xxzzzzz_1, g_xxxzz_0_xzzzzzz_0, g_xxxzz_0_xzzzzzz_1, g_xxyyzz_0_xxxxxyz_1, g_xxyyzz_0_xxxxyyz_1, g_xxyyzz_0_xxxxyz_1, g_xxyyzz_0_xxxxyzz_1, g_xxyyzz_0_xxxyyyz_1, g_xxyyzz_0_xxxyyz_1, g_xxyyzz_0_xxxyyzz_1, g_xxyyzz_0_xxxyzz_1, g_xxyyzz_0_xxxyzzz_1, g_xxyyzz_0_xxyyyyz_1, g_xxyyzz_0_xxyyyz_1, g_xxyyzz_0_xxyyyzz_1, g_xxyyzz_0_xxyyzz_1, g_xxyyzz_0_xxyyzzz_1, g_xxyyzz_0_xxyzzz_1, g_xxyyzz_0_xxyzzzz_1, g_xxyyzz_0_xyyyyyz_1, g_xxyyzz_0_xyyyyz_1, g_xxyyzz_0_xyyyyzz_1, g_xxyyzz_0_xyyyzz_1, g_xxyyzz_0_xyyyzzz_1, g_xxyyzz_0_xyyzzz_1, g_xxyyzz_0_xyyzzzz_1, g_xxyyzz_0_xyzzzz_1, g_xxyyzz_0_xyzzzzz_1, g_xxyyzz_0_yyyyyyy_1, g_xxyyzz_0_yyyyyyz_1, g_xxyyzz_0_yyyyyz_1, g_xxyyzz_0_yyyyyzz_1, g_xxyyzz_0_yyyyzz_1, g_xxyyzz_0_yyyyzzz_1, g_xxyyzz_0_yyyzzz_1, g_xxyyzz_0_yyyzzzz_1, g_xxyyzz_0_yyzzzz_1, g_xxyyzz_0_yyzzzzz_1, g_xxyyzz_0_yzzzzz_1, g_xxyyzz_0_yzzzzzz_1, g_xxyyzz_0_zzzzzzz_1, g_xyyzz_0_xxxxxyz_0, g_xyyzz_0_xxxxxyz_1, g_xyyzz_0_xxxxyyz_0, g_xyyzz_0_xxxxyyz_1, g_xyyzz_0_xxxxyzz_0, g_xyyzz_0_xxxxyzz_1, g_xyyzz_0_xxxyyyz_0, g_xyyzz_0_xxxyyyz_1, g_xyyzz_0_xxxyyzz_0, g_xyyzz_0_xxxyyzz_1, g_xyyzz_0_xxxyzzz_0, g_xyyzz_0_xxxyzzz_1, g_xyyzz_0_xxyyyyz_0, g_xyyzz_0_xxyyyyz_1, g_xyyzz_0_xxyyyzz_0, g_xyyzz_0_xxyyyzz_1, g_xyyzz_0_xxyyzzz_0, g_xyyzz_0_xxyyzzz_1, g_xyyzz_0_xxyzzzz_0, g_xyyzz_0_xxyzzzz_1, g_xyyzz_0_xyyyyyz_0, g_xyyzz_0_xyyyyyz_1, g_xyyzz_0_xyyyyzz_0, g_xyyzz_0_xyyyyzz_1, g_xyyzz_0_xyyyzzz_0, g_xyyzz_0_xyyyzzz_1, g_xyyzz_0_xyyzzzz_0, g_xyyzz_0_xyyzzzz_1, g_xyyzz_0_xyzzzzz_0, g_xyyzz_0_xyzzzzz_1, g_xyyzz_0_yyyyyyy_0, g_xyyzz_0_yyyyyyy_1, g_xyyzz_0_yyyyyyz_0, g_xyyzz_0_yyyyyyz_1, g_xyyzz_0_yyyyyzz_0, g_xyyzz_0_yyyyyzz_1, g_xyyzz_0_yyyyzzz_0, g_xyyzz_0_yyyyzzz_1, g_xyyzz_0_yyyzzzz_0, g_xyyzz_0_yyyzzzz_1, g_xyyzz_0_yyzzzzz_0, g_xyyzz_0_yyzzzzz_1, g_xyyzz_0_yzzzzzz_0, g_xyyzz_0_yzzzzzz_1, g_xyyzz_0_zzzzzzz_0, g_xyyzz_0_zzzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyyzz_0_xxxxxxx_0[i] = g_xxxzz_0_xxxxxxx_0[i] * fbe_0 - g_xxxzz_0_xxxxxxx_1[i] * fz_be_0 + g_xxxyzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxxxxy_0[i] = g_xxxyy_0_xxxxxxy_0[i] * fbe_0 - g_xxxyy_0_xxxxxxy_1[i] * fz_be_0 + g_xxxyyz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxxxxz_0[i] = g_xxxzz_0_xxxxxxz_0[i] * fbe_0 - g_xxxzz_0_xxxxxxz_1[i] * fz_be_0 + g_xxxyzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxxxyy_0[i] = g_xxxyy_0_xxxxxyy_0[i] * fbe_0 - g_xxxyy_0_xxxxxyy_1[i] * fz_be_0 + g_xxxyyz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxxxyz_0[i] = 2.0 * g_xyyzz_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxyyzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxxxzz_0[i] = g_xxxzz_0_xxxxxzz_0[i] * fbe_0 - g_xxxzz_0_xxxxxzz_1[i] * fz_be_0 + g_xxxyzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxxyyy_0[i] = g_xxxyy_0_xxxxyyy_0[i] * fbe_0 - g_xxxyy_0_xxxxyyy_1[i] * fz_be_0 + g_xxxyyz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxxyyz_0[i] = 2.0 * g_xyyzz_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxyyzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxxyzz_0[i] = 2.0 * g_xyyzz_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxyyzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxxzzz_0[i] = g_xxxzz_0_xxxxzzz_0[i] * fbe_0 - g_xxxzz_0_xxxxzzz_1[i] * fz_be_0 + g_xxxyzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxxyyyy_0[i] = g_xxxyy_0_xxxyyyy_0[i] * fbe_0 - g_xxxyy_0_xxxyyyy_1[i] * fz_be_0 + g_xxxyyz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxxyyyz_0[i] = 2.0 * g_xyyzz_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxyyzz_0[i] = 2.0 * g_xyyzz_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxyzzz_0[i] = 2.0 * g_xyyzz_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxyyzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxxzzzz_0[i] = g_xxxzz_0_xxxzzzz_0[i] * fbe_0 - g_xxxzz_0_xxxzzzz_1[i] * fz_be_0 + g_xxxyzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xxyyyyy_0[i] = g_xxxyy_0_xxyyyyy_0[i] * fbe_0 - g_xxxyy_0_xxyyyyy_1[i] * fz_be_0 + g_xxxyyz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xxyyyyz_0[i] = 2.0 * g_xyyzz_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxyyyzz_0[i] = 2.0 * g_xyyzz_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxyyzzz_0[i] = 2.0 * g_xyyzz_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxyzzzz_0[i] = 2.0 * g_xyyzz_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxyyzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xxzzzzz_0[i] = g_xxxzz_0_xxzzzzz_0[i] * fbe_0 - g_xxxzz_0_xxzzzzz_1[i] * fz_be_0 + g_xxxyzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_xyyyyyy_0[i] = g_xxxyy_0_xyyyyyy_0[i] * fbe_0 - g_xxxyy_0_xyyyyyy_1[i] * fz_be_0 + g_xxxyyz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxyyzz_0_xyyyyyz_0[i] = 2.0 * g_xyyzz_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyyyyz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyyyyzz_0[i] = 2.0 * g_xyyzz_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyyyzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyyyzzz_0[i] = 2.0 * g_xyyzz_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyyzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyyzzzz_0[i] = 2.0 * g_xyyzz_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyyzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xyzzzzz_0[i] = 2.0 * g_xyyzz_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_xyzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxyyzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_xzzzzzz_0[i] = g_xxxzz_0_xzzzzzz_0[i] * fbe_0 - g_xxxzz_0_xzzzzzz_1[i] * fz_be_0 + g_xxxyzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxyyzz_0_yyyyyyy_0[i] = 2.0 * g_xyyzz_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyyyy_1[i] * fz_be_0 + g_xxyyzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyyyyz_0[i] = 2.0 * g_xyyzz_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyyyz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyyyzz_0[i] = 2.0 * g_xyyzz_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyyzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyyzzz_0[i] = 2.0 * g_xyyzz_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyyzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyyzzzz_0[i] = 2.0 * g_xyyzz_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyyzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yyzzzzz_0[i] = 2.0 * g_xyyzz_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yyzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_yzzzzzz_0[i] = 2.0 * g_xyyzz_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_yzzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxyyzz_0_zzzzzzz_0[i] = 2.0 * g_xyyzz_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_xyyzz_0_zzzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 468-504 components of targeted buffer : KSK

    auto g_xxxyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 468);

    auto g_xxxyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 469);

    auto g_xxxyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 470);

    auto g_xxxyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 471);

    auto g_xxxyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 472);

    auto g_xxxyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 473);

    auto g_xxxyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 474);

    auto g_xxxyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 475);

    auto g_xxxyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 476);

    auto g_xxxyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 477);

    auto g_xxxyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 478);

    auto g_xxxyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 479);

    auto g_xxxyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 480);

    auto g_xxxyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 481);

    auto g_xxxyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 482);

    auto g_xxxyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 483);

    auto g_xxxyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 484);

    auto g_xxxyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 485);

    auto g_xxxyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 486);

    auto g_xxxyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 487);

    auto g_xxxyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 488);

    auto g_xxxyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 489);

    auto g_xxxyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 490);

    auto g_xxxyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 491);

    auto g_xxxyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 492);

    auto g_xxxyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 493);

    auto g_xxxyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 494);

    auto g_xxxyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 495);

    auto g_xxxyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 496);

    auto g_xxxyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 497);

    auto g_xxxyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 498);

    auto g_xxxyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 499);

    auto g_xxxyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 500);

    auto g_xxxyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 501);

    auto g_xxxyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 502);

    auto g_xxxyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 503);

    #pragma omp simd aligned(g_xxxyzzz_0_xxxxxxx_0, g_xxxyzzz_0_xxxxxxy_0, g_xxxyzzz_0_xxxxxxz_0, g_xxxyzzz_0_xxxxxyy_0, g_xxxyzzz_0_xxxxxyz_0, g_xxxyzzz_0_xxxxxzz_0, g_xxxyzzz_0_xxxxyyy_0, g_xxxyzzz_0_xxxxyyz_0, g_xxxyzzz_0_xxxxyzz_0, g_xxxyzzz_0_xxxxzzz_0, g_xxxyzzz_0_xxxyyyy_0, g_xxxyzzz_0_xxxyyyz_0, g_xxxyzzz_0_xxxyyzz_0, g_xxxyzzz_0_xxxyzzz_0, g_xxxyzzz_0_xxxzzzz_0, g_xxxyzzz_0_xxyyyyy_0, g_xxxyzzz_0_xxyyyyz_0, g_xxxyzzz_0_xxyyyzz_0, g_xxxyzzz_0_xxyyzzz_0, g_xxxyzzz_0_xxyzzzz_0, g_xxxyzzz_0_xxzzzzz_0, g_xxxyzzz_0_xyyyyyy_0, g_xxxyzzz_0_xyyyyyz_0, g_xxxyzzz_0_xyyyyzz_0, g_xxxyzzz_0_xyyyzzz_0, g_xxxyzzz_0_xyyzzzz_0, g_xxxyzzz_0_xyzzzzz_0, g_xxxyzzz_0_xzzzzzz_0, g_xxxyzzz_0_yyyyyyy_0, g_xxxyzzz_0_yyyyyyz_0, g_xxxyzzz_0_yyyyyzz_0, g_xxxyzzz_0_yyyyzzz_0, g_xxxyzzz_0_yyyzzzz_0, g_xxxyzzz_0_yyzzzzz_0, g_xxxyzzz_0_yzzzzzz_0, g_xxxyzzz_0_zzzzzzz_0, g_xxxzzz_0_xxxxxx_1, g_xxxzzz_0_xxxxxxx_1, g_xxxzzz_0_xxxxxxy_1, g_xxxzzz_0_xxxxxxz_1, g_xxxzzz_0_xxxxxy_1, g_xxxzzz_0_xxxxxyy_1, g_xxxzzz_0_xxxxxyz_1, g_xxxzzz_0_xxxxxz_1, g_xxxzzz_0_xxxxxzz_1, g_xxxzzz_0_xxxxyy_1, g_xxxzzz_0_xxxxyyy_1, g_xxxzzz_0_xxxxyyz_1, g_xxxzzz_0_xxxxyz_1, g_xxxzzz_0_xxxxyzz_1, g_xxxzzz_0_xxxxzz_1, g_xxxzzz_0_xxxxzzz_1, g_xxxzzz_0_xxxyyy_1, g_xxxzzz_0_xxxyyyy_1, g_xxxzzz_0_xxxyyyz_1, g_xxxzzz_0_xxxyyz_1, g_xxxzzz_0_xxxyyzz_1, g_xxxzzz_0_xxxyzz_1, g_xxxzzz_0_xxxyzzz_1, g_xxxzzz_0_xxxzzz_1, g_xxxzzz_0_xxxzzzz_1, g_xxxzzz_0_xxyyyy_1, g_xxxzzz_0_xxyyyyy_1, g_xxxzzz_0_xxyyyyz_1, g_xxxzzz_0_xxyyyz_1, g_xxxzzz_0_xxyyyzz_1, g_xxxzzz_0_xxyyzz_1, g_xxxzzz_0_xxyyzzz_1, g_xxxzzz_0_xxyzzz_1, g_xxxzzz_0_xxyzzzz_1, g_xxxzzz_0_xxzzzz_1, g_xxxzzz_0_xxzzzzz_1, g_xxxzzz_0_xyyyyy_1, g_xxxzzz_0_xyyyyyy_1, g_xxxzzz_0_xyyyyyz_1, g_xxxzzz_0_xyyyyz_1, g_xxxzzz_0_xyyyyzz_1, g_xxxzzz_0_xyyyzz_1, g_xxxzzz_0_xyyyzzz_1, g_xxxzzz_0_xyyzzz_1, g_xxxzzz_0_xyyzzzz_1, g_xxxzzz_0_xyzzzz_1, g_xxxzzz_0_xyzzzzz_1, g_xxxzzz_0_xzzzzz_1, g_xxxzzz_0_xzzzzzz_1, g_xxxzzz_0_yyyyyy_1, g_xxxzzz_0_yyyyyyy_1, g_xxxzzz_0_yyyyyyz_1, g_xxxzzz_0_yyyyyz_1, g_xxxzzz_0_yyyyyzz_1, g_xxxzzz_0_yyyyzz_1, g_xxxzzz_0_yyyyzzz_1, g_xxxzzz_0_yyyzzz_1, g_xxxzzz_0_yyyzzzz_1, g_xxxzzz_0_yyzzzz_1, g_xxxzzz_0_yyzzzzz_1, g_xxxzzz_0_yzzzzz_1, g_xxxzzz_0_yzzzzzz_1, g_xxxzzz_0_zzzzzz_1, g_xxxzzz_0_zzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzzz_0_xxxxxxx_0[i] = g_xxxzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxxy_0[i] = g_xxxzzz_0_xxxxxx_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxxy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxxz_0[i] = g_xxxzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxyy_0[i] = 2.0 * g_xxxzzz_0_xxxxxy_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxyz_0[i] = g_xxxzzz_0_xxxxxz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxxyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxxzz_0[i] = g_xxxzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxyyy_0[i] = 3.0 * g_xxxzzz_0_xxxxyy_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxyyz_0[i] = 2.0 * g_xxxzzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxyzz_0[i] = g_xxxzzz_0_xxxxzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxxyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxxzzz_0[i] = g_xxxzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyyyy_0[i] = 4.0 * g_xxxzzz_0_xxxyyy_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyyyz_0[i] = 3.0 * g_xxxzzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyyzz_0[i] = 2.0 * g_xxxzzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxyzzz_0[i] = g_xxxzzz_0_xxxzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxxyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxxzzzz_0[i] = g_xxxzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyyyy_0[i] = 5.0 * g_xxxzzz_0_xxyyyy_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyyyz_0[i] = 4.0 * g_xxxzzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyyzz_0[i] = 3.0 * g_xxxzzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyyzzz_0[i] = 2.0 * g_xxxzzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxyzzzz_0[i] = g_xxxzzz_0_xxzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xxyzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xxzzzzz_0[i] = g_xxxzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyyyy_0[i] = 6.0 * g_xxxzzz_0_xyyyyy_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyyyz_0[i] = 5.0 * g_xxxzzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyyzz_0[i] = 4.0 * g_xxxzzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyyzzz_0[i] = 3.0 * g_xxxzzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyyzzzz_0[i] = 2.0 * g_xxxzzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyyzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xyzzzzz_0[i] = g_xxxzzz_0_xzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_xyzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_xzzzzzz_0[i] = g_xxxzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyyyy_0[i] = 7.0 * g_xxxzzz_0_yyyyyy_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyyyy_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyyyz_0[i] = 6.0 * g_xxxzzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyyyz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyyzz_0[i] = 5.0 * g_xxxzzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyyzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyyzzz_0[i] = 4.0 * g_xxxzzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyyzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyyzzzz_0[i] = 3.0 * g_xxxzzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyyzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yyzzzzz_0[i] = 2.0 * g_xxxzzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yyzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_yzzzzzz_0[i] = g_xxxzzz_0_zzzzzz_1[i] * fi_acd_0 + g_xxxzzz_0_yzzzzzz_1[i] * wa_y[i];

        g_xxxyzzz_0_zzzzzzz_0[i] = g_xxxzzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 504-540 components of targeted buffer : KSK

    auto g_xxxzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 504);

    auto g_xxxzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 505);

    auto g_xxxzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 506);

    auto g_xxxzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 507);

    auto g_xxxzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 508);

    auto g_xxxzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 509);

    auto g_xxxzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 510);

    auto g_xxxzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 511);

    auto g_xxxzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 512);

    auto g_xxxzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 513);

    auto g_xxxzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 514);

    auto g_xxxzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 515);

    auto g_xxxzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 516);

    auto g_xxxzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 517);

    auto g_xxxzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 518);

    auto g_xxxzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 519);

    auto g_xxxzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 520);

    auto g_xxxzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 521);

    auto g_xxxzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 522);

    auto g_xxxzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 523);

    auto g_xxxzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 524);

    auto g_xxxzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 525);

    auto g_xxxzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 526);

    auto g_xxxzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 527);

    auto g_xxxzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 528);

    auto g_xxxzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 529);

    auto g_xxxzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 530);

    auto g_xxxzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 531);

    auto g_xxxzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 532);

    auto g_xxxzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 533);

    auto g_xxxzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 534);

    auto g_xxxzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 535);

    auto g_xxxzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 536);

    auto g_xxxzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 537);

    auto g_xxxzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 538);

    auto g_xxxzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 539);

    #pragma omp simd aligned(g_xxxzz_0_xxxxxxx_0, g_xxxzz_0_xxxxxxx_1, g_xxxzz_0_xxxxxxy_0, g_xxxzz_0_xxxxxxy_1, g_xxxzz_0_xxxxxyy_0, g_xxxzz_0_xxxxxyy_1, g_xxxzz_0_xxxxyyy_0, g_xxxzz_0_xxxxyyy_1, g_xxxzz_0_xxxyyyy_0, g_xxxzz_0_xxxyyyy_1, g_xxxzz_0_xxyyyyy_0, g_xxxzz_0_xxyyyyy_1, g_xxxzz_0_xyyyyyy_0, g_xxxzz_0_xyyyyyy_1, g_xxxzzz_0_xxxxxxx_1, g_xxxzzz_0_xxxxxxy_1, g_xxxzzz_0_xxxxxyy_1, g_xxxzzz_0_xxxxyyy_1, g_xxxzzz_0_xxxyyyy_1, g_xxxzzz_0_xxyyyyy_1, g_xxxzzz_0_xyyyyyy_1, g_xxxzzzz_0_xxxxxxx_0, g_xxxzzzz_0_xxxxxxy_0, g_xxxzzzz_0_xxxxxxz_0, g_xxxzzzz_0_xxxxxyy_0, g_xxxzzzz_0_xxxxxyz_0, g_xxxzzzz_0_xxxxxzz_0, g_xxxzzzz_0_xxxxyyy_0, g_xxxzzzz_0_xxxxyyz_0, g_xxxzzzz_0_xxxxyzz_0, g_xxxzzzz_0_xxxxzzz_0, g_xxxzzzz_0_xxxyyyy_0, g_xxxzzzz_0_xxxyyyz_0, g_xxxzzzz_0_xxxyyzz_0, g_xxxzzzz_0_xxxyzzz_0, g_xxxzzzz_0_xxxzzzz_0, g_xxxzzzz_0_xxyyyyy_0, g_xxxzzzz_0_xxyyyyz_0, g_xxxzzzz_0_xxyyyzz_0, g_xxxzzzz_0_xxyyzzz_0, g_xxxzzzz_0_xxyzzzz_0, g_xxxzzzz_0_xxzzzzz_0, g_xxxzzzz_0_xyyyyyy_0, g_xxxzzzz_0_xyyyyyz_0, g_xxxzzzz_0_xyyyyzz_0, g_xxxzzzz_0_xyyyzzz_0, g_xxxzzzz_0_xyyzzzz_0, g_xxxzzzz_0_xyzzzzz_0, g_xxxzzzz_0_xzzzzzz_0, g_xxxzzzz_0_yyyyyyy_0, g_xxxzzzz_0_yyyyyyz_0, g_xxxzzzz_0_yyyyyzz_0, g_xxxzzzz_0_yyyyzzz_0, g_xxxzzzz_0_yyyzzzz_0, g_xxxzzzz_0_yyzzzzz_0, g_xxxzzzz_0_yzzzzzz_0, g_xxxzzzz_0_zzzzzzz_0, g_xxzzzz_0_xxxxxxz_1, g_xxzzzz_0_xxxxxyz_1, g_xxzzzz_0_xxxxxz_1, g_xxzzzz_0_xxxxxzz_1, g_xxzzzz_0_xxxxyyz_1, g_xxzzzz_0_xxxxyz_1, g_xxzzzz_0_xxxxyzz_1, g_xxzzzz_0_xxxxzz_1, g_xxzzzz_0_xxxxzzz_1, g_xxzzzz_0_xxxyyyz_1, g_xxzzzz_0_xxxyyz_1, g_xxzzzz_0_xxxyyzz_1, g_xxzzzz_0_xxxyzz_1, g_xxzzzz_0_xxxyzzz_1, g_xxzzzz_0_xxxzzz_1, g_xxzzzz_0_xxxzzzz_1, g_xxzzzz_0_xxyyyyz_1, g_xxzzzz_0_xxyyyz_1, g_xxzzzz_0_xxyyyzz_1, g_xxzzzz_0_xxyyzz_1, g_xxzzzz_0_xxyyzzz_1, g_xxzzzz_0_xxyzzz_1, g_xxzzzz_0_xxyzzzz_1, g_xxzzzz_0_xxzzzz_1, g_xxzzzz_0_xxzzzzz_1, g_xxzzzz_0_xyyyyyz_1, g_xxzzzz_0_xyyyyz_1, g_xxzzzz_0_xyyyyzz_1, g_xxzzzz_0_xyyyzz_1, g_xxzzzz_0_xyyyzzz_1, g_xxzzzz_0_xyyzzz_1, g_xxzzzz_0_xyyzzzz_1, g_xxzzzz_0_xyzzzz_1, g_xxzzzz_0_xyzzzzz_1, g_xxzzzz_0_xzzzzz_1, g_xxzzzz_0_xzzzzzz_1, g_xxzzzz_0_yyyyyyy_1, g_xxzzzz_0_yyyyyyz_1, g_xxzzzz_0_yyyyyz_1, g_xxzzzz_0_yyyyyzz_1, g_xxzzzz_0_yyyyzz_1, g_xxzzzz_0_yyyyzzz_1, g_xxzzzz_0_yyyzzz_1, g_xxzzzz_0_yyyzzzz_1, g_xxzzzz_0_yyzzzz_1, g_xxzzzz_0_yyzzzzz_1, g_xxzzzz_0_yzzzzz_1, g_xxzzzz_0_yzzzzzz_1, g_xxzzzz_0_zzzzzz_1, g_xxzzzz_0_zzzzzzz_1, g_xzzzz_0_xxxxxxz_0, g_xzzzz_0_xxxxxxz_1, g_xzzzz_0_xxxxxyz_0, g_xzzzz_0_xxxxxyz_1, g_xzzzz_0_xxxxxzz_0, g_xzzzz_0_xxxxxzz_1, g_xzzzz_0_xxxxyyz_0, g_xzzzz_0_xxxxyyz_1, g_xzzzz_0_xxxxyzz_0, g_xzzzz_0_xxxxyzz_1, g_xzzzz_0_xxxxzzz_0, g_xzzzz_0_xxxxzzz_1, g_xzzzz_0_xxxyyyz_0, g_xzzzz_0_xxxyyyz_1, g_xzzzz_0_xxxyyzz_0, g_xzzzz_0_xxxyyzz_1, g_xzzzz_0_xxxyzzz_0, g_xzzzz_0_xxxyzzz_1, g_xzzzz_0_xxxzzzz_0, g_xzzzz_0_xxxzzzz_1, g_xzzzz_0_xxyyyyz_0, g_xzzzz_0_xxyyyyz_1, g_xzzzz_0_xxyyyzz_0, g_xzzzz_0_xxyyyzz_1, g_xzzzz_0_xxyyzzz_0, g_xzzzz_0_xxyyzzz_1, g_xzzzz_0_xxyzzzz_0, g_xzzzz_0_xxyzzzz_1, g_xzzzz_0_xxzzzzz_0, g_xzzzz_0_xxzzzzz_1, g_xzzzz_0_xyyyyyz_0, g_xzzzz_0_xyyyyyz_1, g_xzzzz_0_xyyyyzz_0, g_xzzzz_0_xyyyyzz_1, g_xzzzz_0_xyyyzzz_0, g_xzzzz_0_xyyyzzz_1, g_xzzzz_0_xyyzzzz_0, g_xzzzz_0_xyyzzzz_1, g_xzzzz_0_xyzzzzz_0, g_xzzzz_0_xyzzzzz_1, g_xzzzz_0_xzzzzzz_0, g_xzzzz_0_xzzzzzz_1, g_xzzzz_0_yyyyyyy_0, g_xzzzz_0_yyyyyyy_1, g_xzzzz_0_yyyyyyz_0, g_xzzzz_0_yyyyyyz_1, g_xzzzz_0_yyyyyzz_0, g_xzzzz_0_yyyyyzz_1, g_xzzzz_0_yyyyzzz_0, g_xzzzz_0_yyyyzzz_1, g_xzzzz_0_yyyzzzz_0, g_xzzzz_0_yyyzzzz_1, g_xzzzz_0_yyzzzzz_0, g_xzzzz_0_yyzzzzz_1, g_xzzzz_0_yzzzzzz_0, g_xzzzz_0_yzzzzzz_1, g_xzzzz_0_zzzzzzz_0, g_xzzzz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzzzz_0_xxxxxxx_0[i] = 3.0 * g_xxxzz_0_xxxxxxx_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxxxx_1[i] * fz_be_0 + g_xxxzzz_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxxxy_0[i] = 3.0 * g_xxxzz_0_xxxxxxy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxxxy_1[i] * fz_be_0 + g_xxxzzz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxxxz_0[i] = 2.0 * g_xzzzz_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xxzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxxyy_0[i] = 3.0 * g_xxxzz_0_xxxxxyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxxyy_1[i] * fz_be_0 + g_xxxzzz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxxyz_0[i] = 2.0 * g_xzzzz_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxxzz_0[i] = 2.0 * g_xzzzz_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xxzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxyyy_0[i] = 3.0 * g_xxxzz_0_xxxxyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxxyyy_1[i] * fz_be_0 + g_xxxzzz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxxyyz_0[i] = 2.0 * g_xzzzz_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxyzz_0[i] = 2.0 * g_xzzzz_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxxzzz_0[i] = 2.0 * g_xzzzz_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xxzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxyyyy_0[i] = 3.0 * g_xxxzz_0_xxxyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxxyyyy_1[i] * fz_be_0 + g_xxxzzz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxxyyyz_0[i] = 2.0 * g_xzzzz_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxyyzz_0[i] = 2.0 * g_xzzzz_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxyzzz_0[i] = 2.0 * g_xzzzz_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxxzzzz_0[i] = 2.0 * g_xzzzz_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xxzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyyyyy_0[i] = 3.0 * g_xxxzz_0_xxyyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xxyyyyy_1[i] * fz_be_0 + g_xxxzzz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xxyyyyz_0[i] = 2.0 * g_xzzzz_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyyyzz_0[i] = 2.0 * g_xzzzz_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyyzzz_0[i] = 2.0 * g_xzzzz_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxyzzzz_0[i] = 2.0 * g_xzzzz_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xxzzzzz_0[i] = 2.0 * g_xzzzz_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyyyyy_0[i] = 3.0 * g_xxxzz_0_xyyyyyy_0[i] * fbe_0 - 3.0 * g_xxxzz_0_xyyyyyy_1[i] * fz_be_0 + g_xxxzzz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxzzzz_0_xyyyyyz_0[i] = 2.0 * g_xzzzz_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyyyyz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyyyzz_0[i] = 2.0 * g_xzzzz_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyyyzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyyzzz_0[i] = 2.0 * g_xzzzz_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyyzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyyzzzz_0[i] = 2.0 * g_xzzzz_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyyzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xyzzzzz_0[i] = 2.0 * g_xzzzz_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xyzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_xzzzzzz_0[i] = 2.0 * g_xzzzz_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_xzzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyyyy_0[i] = 2.0 * g_xzzzz_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyyyy_1[i] * fz_be_0 + g_xxzzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyyyz_0[i] = 2.0 * g_xzzzz_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyyyz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyyzz_0[i] = 2.0 * g_xzzzz_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyyzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyyzzz_0[i] = 2.0 * g_xzzzz_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyyzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyyzzzz_0[i] = 2.0 * g_xzzzz_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyyzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yyzzzzz_0[i] = 2.0 * g_xzzzz_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yyzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_yzzzzzz_0[i] = 2.0 * g_xzzzz_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_yzzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxzzzz_0_zzzzzzz_0[i] = 2.0 * g_xzzzz_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_xzzzz_0_zzzzzzz_1[i] * fz_be_0 + g_xxzzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 540-576 components of targeted buffer : KSK

    auto g_xxyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 540);

    auto g_xxyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 541);

    auto g_xxyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 542);

    auto g_xxyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 543);

    auto g_xxyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 544);

    auto g_xxyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 545);

    auto g_xxyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 546);

    auto g_xxyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 547);

    auto g_xxyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 548);

    auto g_xxyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 549);

    auto g_xxyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 550);

    auto g_xxyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 551);

    auto g_xxyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 552);

    auto g_xxyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 553);

    auto g_xxyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 554);

    auto g_xxyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 555);

    auto g_xxyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 556);

    auto g_xxyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 557);

    auto g_xxyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 558);

    auto g_xxyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 559);

    auto g_xxyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 560);

    auto g_xxyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 561);

    auto g_xxyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 562);

    auto g_xxyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 563);

    auto g_xxyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 564);

    auto g_xxyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 565);

    auto g_xxyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 566);

    auto g_xxyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 567);

    auto g_xxyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 568);

    auto g_xxyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 569);

    auto g_xxyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 570);

    auto g_xxyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 571);

    auto g_xxyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 572);

    auto g_xxyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 573);

    auto g_xxyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 574);

    auto g_xxyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 575);

    #pragma omp simd aligned(g_xxyyy_0_xxxxxxx_0, g_xxyyy_0_xxxxxxx_1, g_xxyyy_0_xxxxxxz_0, g_xxyyy_0_xxxxxxz_1, g_xxyyy_0_xxxxxzz_0, g_xxyyy_0_xxxxxzz_1, g_xxyyy_0_xxxxzzz_0, g_xxyyy_0_xxxxzzz_1, g_xxyyy_0_xxxzzzz_0, g_xxyyy_0_xxxzzzz_1, g_xxyyy_0_xxzzzzz_0, g_xxyyy_0_xxzzzzz_1, g_xxyyy_0_xzzzzzz_0, g_xxyyy_0_xzzzzzz_1, g_xxyyyy_0_xxxxxxx_1, g_xxyyyy_0_xxxxxxz_1, g_xxyyyy_0_xxxxxzz_1, g_xxyyyy_0_xxxxzzz_1, g_xxyyyy_0_xxxzzzz_1, g_xxyyyy_0_xxzzzzz_1, g_xxyyyy_0_xzzzzzz_1, g_xxyyyyy_0_xxxxxxx_0, g_xxyyyyy_0_xxxxxxy_0, g_xxyyyyy_0_xxxxxxz_0, g_xxyyyyy_0_xxxxxyy_0, g_xxyyyyy_0_xxxxxyz_0, g_xxyyyyy_0_xxxxxzz_0, g_xxyyyyy_0_xxxxyyy_0, g_xxyyyyy_0_xxxxyyz_0, g_xxyyyyy_0_xxxxyzz_0, g_xxyyyyy_0_xxxxzzz_0, g_xxyyyyy_0_xxxyyyy_0, g_xxyyyyy_0_xxxyyyz_0, g_xxyyyyy_0_xxxyyzz_0, g_xxyyyyy_0_xxxyzzz_0, g_xxyyyyy_0_xxxzzzz_0, g_xxyyyyy_0_xxyyyyy_0, g_xxyyyyy_0_xxyyyyz_0, g_xxyyyyy_0_xxyyyzz_0, g_xxyyyyy_0_xxyyzzz_0, g_xxyyyyy_0_xxyzzzz_0, g_xxyyyyy_0_xxzzzzz_0, g_xxyyyyy_0_xyyyyyy_0, g_xxyyyyy_0_xyyyyyz_0, g_xxyyyyy_0_xyyyyzz_0, g_xxyyyyy_0_xyyyzzz_0, g_xxyyyyy_0_xyyzzzz_0, g_xxyyyyy_0_xyzzzzz_0, g_xxyyyyy_0_xzzzzzz_0, g_xxyyyyy_0_yyyyyyy_0, g_xxyyyyy_0_yyyyyyz_0, g_xxyyyyy_0_yyyyyzz_0, g_xxyyyyy_0_yyyyzzz_0, g_xxyyyyy_0_yyyzzzz_0, g_xxyyyyy_0_yyzzzzz_0, g_xxyyyyy_0_yzzzzzz_0, g_xxyyyyy_0_zzzzzzz_0, g_xyyyyy_0_xxxxxxy_1, g_xyyyyy_0_xxxxxy_1, g_xyyyyy_0_xxxxxyy_1, g_xyyyyy_0_xxxxxyz_1, g_xyyyyy_0_xxxxyy_1, g_xyyyyy_0_xxxxyyy_1, g_xyyyyy_0_xxxxyyz_1, g_xyyyyy_0_xxxxyz_1, g_xyyyyy_0_xxxxyzz_1, g_xyyyyy_0_xxxyyy_1, g_xyyyyy_0_xxxyyyy_1, g_xyyyyy_0_xxxyyyz_1, g_xyyyyy_0_xxxyyz_1, g_xyyyyy_0_xxxyyzz_1, g_xyyyyy_0_xxxyzz_1, g_xyyyyy_0_xxxyzzz_1, g_xyyyyy_0_xxyyyy_1, g_xyyyyy_0_xxyyyyy_1, g_xyyyyy_0_xxyyyyz_1, g_xyyyyy_0_xxyyyz_1, g_xyyyyy_0_xxyyyzz_1, g_xyyyyy_0_xxyyzz_1, g_xyyyyy_0_xxyyzzz_1, g_xyyyyy_0_xxyzzz_1, g_xyyyyy_0_xxyzzzz_1, g_xyyyyy_0_xyyyyy_1, g_xyyyyy_0_xyyyyyy_1, g_xyyyyy_0_xyyyyyz_1, g_xyyyyy_0_xyyyyz_1, g_xyyyyy_0_xyyyyzz_1, g_xyyyyy_0_xyyyzz_1, g_xyyyyy_0_xyyyzzz_1, g_xyyyyy_0_xyyzzz_1, g_xyyyyy_0_xyyzzzz_1, g_xyyyyy_0_xyzzzz_1, g_xyyyyy_0_xyzzzzz_1, g_xyyyyy_0_yyyyyy_1, g_xyyyyy_0_yyyyyyy_1, g_xyyyyy_0_yyyyyyz_1, g_xyyyyy_0_yyyyyz_1, g_xyyyyy_0_yyyyyzz_1, g_xyyyyy_0_yyyyzz_1, g_xyyyyy_0_yyyyzzz_1, g_xyyyyy_0_yyyzzz_1, g_xyyyyy_0_yyyzzzz_1, g_xyyyyy_0_yyzzzz_1, g_xyyyyy_0_yyzzzzz_1, g_xyyyyy_0_yzzzzz_1, g_xyyyyy_0_yzzzzzz_1, g_xyyyyy_0_zzzzzzz_1, g_yyyyy_0_xxxxxxy_0, g_yyyyy_0_xxxxxxy_1, g_yyyyy_0_xxxxxyy_0, g_yyyyy_0_xxxxxyy_1, g_yyyyy_0_xxxxxyz_0, g_yyyyy_0_xxxxxyz_1, g_yyyyy_0_xxxxyyy_0, g_yyyyy_0_xxxxyyy_1, g_yyyyy_0_xxxxyyz_0, g_yyyyy_0_xxxxyyz_1, g_yyyyy_0_xxxxyzz_0, g_yyyyy_0_xxxxyzz_1, g_yyyyy_0_xxxyyyy_0, g_yyyyy_0_xxxyyyy_1, g_yyyyy_0_xxxyyyz_0, g_yyyyy_0_xxxyyyz_1, g_yyyyy_0_xxxyyzz_0, g_yyyyy_0_xxxyyzz_1, g_yyyyy_0_xxxyzzz_0, g_yyyyy_0_xxxyzzz_1, g_yyyyy_0_xxyyyyy_0, g_yyyyy_0_xxyyyyy_1, g_yyyyy_0_xxyyyyz_0, g_yyyyy_0_xxyyyyz_1, g_yyyyy_0_xxyyyzz_0, g_yyyyy_0_xxyyyzz_1, g_yyyyy_0_xxyyzzz_0, g_yyyyy_0_xxyyzzz_1, g_yyyyy_0_xxyzzzz_0, g_yyyyy_0_xxyzzzz_1, g_yyyyy_0_xyyyyyy_0, g_yyyyy_0_xyyyyyy_1, g_yyyyy_0_xyyyyyz_0, g_yyyyy_0_xyyyyyz_1, g_yyyyy_0_xyyyyzz_0, g_yyyyy_0_xyyyyzz_1, g_yyyyy_0_xyyyzzz_0, g_yyyyy_0_xyyyzzz_1, g_yyyyy_0_xyyzzzz_0, g_yyyyy_0_xyyzzzz_1, g_yyyyy_0_xyzzzzz_0, g_yyyyy_0_xyzzzzz_1, g_yyyyy_0_yyyyyyy_0, g_yyyyy_0_yyyyyyy_1, g_yyyyy_0_yyyyyyz_0, g_yyyyy_0_yyyyyyz_1, g_yyyyy_0_yyyyyzz_0, g_yyyyy_0_yyyyyzz_1, g_yyyyy_0_yyyyzzz_0, g_yyyyy_0_yyyyzzz_1, g_yyyyy_0_yyyzzzz_0, g_yyyyy_0_yyyzzzz_1, g_yyyyy_0_yyzzzzz_0, g_yyyyy_0_yyzzzzz_1, g_yyyyy_0_yzzzzzz_0, g_yyyyy_0_yzzzzzz_1, g_yyyyy_0_zzzzzzz_0, g_yyyyy_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyyy_0_xxxxxxx_0[i] = 4.0 * g_xxyyy_0_xxxxxxx_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxxxx_1[i] * fz_be_0 + g_xxyyyy_0_xxxxxxx_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxxxxy_0[i] = g_yyyyy_0_xxxxxxy_0[i] * fbe_0 - g_yyyyy_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xyyyyy_0_xxxxxy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxxxz_0[i] = 4.0 * g_xxyyy_0_xxxxxxz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxxxz_1[i] * fz_be_0 + g_xxyyyy_0_xxxxxxz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxxxyy_0[i] = g_yyyyy_0_xxxxxyy_0[i] * fbe_0 - g_yyyyy_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xyyyyy_0_xxxxyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxxyz_0[i] = g_yyyyy_0_xxxxxyz_0[i] * fbe_0 - g_yyyyy_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xyyyyy_0_xxxxyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxxzz_0[i] = 4.0 * g_xxyyy_0_xxxxxzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxxzz_1[i] * fz_be_0 + g_xxyyyy_0_xxxxxzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxxyyy_0[i] = g_yyyyy_0_xxxxyyy_0[i] * fbe_0 - g_yyyyy_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xyyyyy_0_xxxyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxyyz_0[i] = g_yyyyy_0_xxxxyyz_0[i] * fbe_0 - g_yyyyy_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xyyyyy_0_xxxyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxyzz_0[i] = g_yyyyy_0_xxxxyzz_0[i] * fbe_0 - g_yyyyy_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xyyyyy_0_xxxyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxxzzz_0[i] = 4.0 * g_xxyyy_0_xxxxzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxxzzz_1[i] * fz_be_0 + g_xxyyyy_0_xxxxzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxxyyyy_0[i] = g_yyyyy_0_xxxyyyy_0[i] * fbe_0 - g_yyyyy_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxyyyz_0[i] = g_yyyyy_0_xxxyyyz_0[i] * fbe_0 - g_yyyyy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxyyzz_0[i] = g_yyyyy_0_xxxyyzz_0[i] * fbe_0 - g_yyyyy_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxyzzz_0[i] = g_yyyyy_0_xxxyzzz_0[i] * fbe_0 - g_yyyyy_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xyyyyy_0_xxyzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxxzzzz_0[i] = 4.0 * g_xxyyy_0_xxxzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxxzzzz_1[i] * fz_be_0 + g_xxyyyy_0_xxxzzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xxyyyyy_0[i] = g_yyyyy_0_xxyyyyy_0[i] * fbe_0 - g_yyyyy_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyyyyz_0[i] = g_yyyyy_0_xxyyyyz_0[i] * fbe_0 - g_yyyyy_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyyyzz_0[i] = g_yyyyy_0_xxyyyzz_0[i] * fbe_0 - g_yyyyy_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyyzzz_0[i] = g_yyyyy_0_xxyyzzz_0[i] * fbe_0 - g_yyyyy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyyzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxyzzzz_0[i] = g_yyyyy_0_xxyzzzz_0[i] * fbe_0 - g_yyyyy_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xyyyyy_0_xyzzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xxzzzzz_0[i] = 4.0 * g_xxyyy_0_xxzzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xxzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_xxzzzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_xyyyyyy_0[i] = g_yyyyy_0_xyyyyyy_0[i] * fbe_0 - g_yyyyy_0_xyyyyyy_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyy_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyyyyz_0[i] = g_yyyyy_0_xyyyyyz_0[i] * fbe_0 - g_yyyyy_0_xyyyyyz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyyyzz_0[i] = g_yyyyy_0_xyyyyzz_0[i] * fbe_0 - g_yyyyy_0_xyyyyzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyyzzz_0[i] = g_yyyyy_0_xyyyzzz_0[i] * fbe_0 - g_yyyyy_0_xyyyzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyyzzzz_0[i] = g_yyyyy_0_xyyzzzz_0[i] * fbe_0 - g_yyyyy_0_xyyzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyzzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xyzzzzz_0[i] = g_yyyyy_0_xyzzzzz_0[i] * fbe_0 - g_yyyyy_0_xyzzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yzzzzz_1[i] * fi_acd_0 + g_xyyyyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_xzzzzzz_0[i] = 4.0 * g_xxyyy_0_xzzzzzz_0[i] * fbe_0 - 4.0 * g_xxyyy_0_xzzzzzz_1[i] * fz_be_0 + g_xxyyyy_0_xzzzzzz_1[i] * wa_y[i];

        g_xxyyyyy_0_yyyyyyy_0[i] = g_yyyyy_0_yyyyyyy_0[i] * fbe_0 - g_yyyyy_0_yyyyyyy_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyyyyz_0[i] = g_yyyyy_0_yyyyyyz_0[i] * fbe_0 - g_yyyyy_0_yyyyyyz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyyyzz_0[i] = g_yyyyy_0_yyyyyzz_0[i] * fbe_0 - g_yyyyy_0_yyyyyzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyyzzz_0[i] = g_yyyyy_0_yyyyzzz_0[i] * fbe_0 - g_yyyyy_0_yyyyzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyyzzzz_0[i] = g_yyyyy_0_yyyzzzz_0[i] * fbe_0 - g_yyyyy_0_yyyzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yyzzzzz_0[i] = g_yyyyy_0_yyzzzzz_0[i] * fbe_0 - g_yyyyy_0_yyzzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_yzzzzzz_0[i] = g_yyyyy_0_yzzzzzz_0[i] * fbe_0 - g_yyyyy_0_yzzzzzz_1[i] * fz_be_0 + g_xyyyyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xxyyyyy_0_zzzzzzz_0[i] = g_yyyyy_0_zzzzzzz_0[i] * fbe_0 - g_yyyyy_0_zzzzzzz_1[i] * fz_be_0 + g_xyyyyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 576-612 components of targeted buffer : KSK

    auto g_xxyyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 576);

    auto g_xxyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 577);

    auto g_xxyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 578);

    auto g_xxyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 579);

    auto g_xxyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 580);

    auto g_xxyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 581);

    auto g_xxyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 582);

    auto g_xxyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 583);

    auto g_xxyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 584);

    auto g_xxyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 585);

    auto g_xxyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 586);

    auto g_xxyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 587);

    auto g_xxyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 588);

    auto g_xxyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 589);

    auto g_xxyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 590);

    auto g_xxyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 591);

    auto g_xxyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 592);

    auto g_xxyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 593);

    auto g_xxyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 594);

    auto g_xxyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 595);

    auto g_xxyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 596);

    auto g_xxyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 597);

    auto g_xxyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 598);

    auto g_xxyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 599);

    auto g_xxyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 600);

    auto g_xxyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 601);

    auto g_xxyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 602);

    auto g_xxyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 603);

    auto g_xxyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 604);

    auto g_xxyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 605);

    auto g_xxyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 606);

    auto g_xxyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 607);

    auto g_xxyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 608);

    auto g_xxyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 609);

    auto g_xxyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 610);

    auto g_xxyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 611);

    #pragma omp simd aligned(g_xxyyyy_0_xxxxxx_1, g_xxyyyy_0_xxxxxxx_1, g_xxyyyy_0_xxxxxxy_1, g_xxyyyy_0_xxxxxxz_1, g_xxyyyy_0_xxxxxy_1, g_xxyyyy_0_xxxxxyy_1, g_xxyyyy_0_xxxxxyz_1, g_xxyyyy_0_xxxxxz_1, g_xxyyyy_0_xxxxxzz_1, g_xxyyyy_0_xxxxyy_1, g_xxyyyy_0_xxxxyyy_1, g_xxyyyy_0_xxxxyyz_1, g_xxyyyy_0_xxxxyz_1, g_xxyyyy_0_xxxxyzz_1, g_xxyyyy_0_xxxxzz_1, g_xxyyyy_0_xxxxzzz_1, g_xxyyyy_0_xxxyyy_1, g_xxyyyy_0_xxxyyyy_1, g_xxyyyy_0_xxxyyyz_1, g_xxyyyy_0_xxxyyz_1, g_xxyyyy_0_xxxyyzz_1, g_xxyyyy_0_xxxyzz_1, g_xxyyyy_0_xxxyzzz_1, g_xxyyyy_0_xxxzzz_1, g_xxyyyy_0_xxxzzzz_1, g_xxyyyy_0_xxyyyy_1, g_xxyyyy_0_xxyyyyy_1, g_xxyyyy_0_xxyyyyz_1, g_xxyyyy_0_xxyyyz_1, g_xxyyyy_0_xxyyyzz_1, g_xxyyyy_0_xxyyzz_1, g_xxyyyy_0_xxyyzzz_1, g_xxyyyy_0_xxyzzz_1, g_xxyyyy_0_xxyzzzz_1, g_xxyyyy_0_xxzzzz_1, g_xxyyyy_0_xxzzzzz_1, g_xxyyyy_0_xyyyyy_1, g_xxyyyy_0_xyyyyyy_1, g_xxyyyy_0_xyyyyyz_1, g_xxyyyy_0_xyyyyz_1, g_xxyyyy_0_xyyyyzz_1, g_xxyyyy_0_xyyyzz_1, g_xxyyyy_0_xyyyzzz_1, g_xxyyyy_0_xyyzzz_1, g_xxyyyy_0_xyyzzzz_1, g_xxyyyy_0_xyzzzz_1, g_xxyyyy_0_xyzzzzz_1, g_xxyyyy_0_xzzzzz_1, g_xxyyyy_0_xzzzzzz_1, g_xxyyyy_0_yyyyyy_1, g_xxyyyy_0_yyyyyyy_1, g_xxyyyy_0_yyyyyyz_1, g_xxyyyy_0_yyyyyz_1, g_xxyyyy_0_yyyyyzz_1, g_xxyyyy_0_yyyyzz_1, g_xxyyyy_0_yyyyzzz_1, g_xxyyyy_0_yyyzzz_1, g_xxyyyy_0_yyyzzzz_1, g_xxyyyy_0_yyzzzz_1, g_xxyyyy_0_yyzzzzz_1, g_xxyyyy_0_yzzzzz_1, g_xxyyyy_0_yzzzzzz_1, g_xxyyyy_0_zzzzzz_1, g_xxyyyy_0_zzzzzzz_1, g_xxyyyyz_0_xxxxxxx_0, g_xxyyyyz_0_xxxxxxy_0, g_xxyyyyz_0_xxxxxxz_0, g_xxyyyyz_0_xxxxxyy_0, g_xxyyyyz_0_xxxxxyz_0, g_xxyyyyz_0_xxxxxzz_0, g_xxyyyyz_0_xxxxyyy_0, g_xxyyyyz_0_xxxxyyz_0, g_xxyyyyz_0_xxxxyzz_0, g_xxyyyyz_0_xxxxzzz_0, g_xxyyyyz_0_xxxyyyy_0, g_xxyyyyz_0_xxxyyyz_0, g_xxyyyyz_0_xxxyyzz_0, g_xxyyyyz_0_xxxyzzz_0, g_xxyyyyz_0_xxxzzzz_0, g_xxyyyyz_0_xxyyyyy_0, g_xxyyyyz_0_xxyyyyz_0, g_xxyyyyz_0_xxyyyzz_0, g_xxyyyyz_0_xxyyzzz_0, g_xxyyyyz_0_xxyzzzz_0, g_xxyyyyz_0_xxzzzzz_0, g_xxyyyyz_0_xyyyyyy_0, g_xxyyyyz_0_xyyyyyz_0, g_xxyyyyz_0_xyyyyzz_0, g_xxyyyyz_0_xyyyzzz_0, g_xxyyyyz_0_xyyzzzz_0, g_xxyyyyz_0_xyzzzzz_0, g_xxyyyyz_0_xzzzzzz_0, g_xxyyyyz_0_yyyyyyy_0, g_xxyyyyz_0_yyyyyyz_0, g_xxyyyyz_0_yyyyyzz_0, g_xxyyyyz_0_yyyyzzz_0, g_xxyyyyz_0_yyyzzzz_0, g_xxyyyyz_0_yyzzzzz_0, g_xxyyyyz_0_yzzzzzz_0, g_xxyyyyz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyyz_0_xxxxxxx_0[i] = g_xxyyyy_0_xxxxxxx_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxxy_0[i] = g_xxyyyy_0_xxxxxxy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxxz_0[i] = g_xxyyyy_0_xxxxxx_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxxz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxyy_0[i] = g_xxyyyy_0_xxxxxyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxyz_0[i] = g_xxyyyy_0_xxxxxy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxxzz_0[i] = 2.0 * g_xxyyyy_0_xxxxxz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxxzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxyyy_0[i] = g_xxyyyy_0_xxxxyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxyyz_0[i] = g_xxyyyy_0_xxxxyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxyzz_0[i] = 2.0 * g_xxyyyy_0_xxxxyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxxzzz_0[i] = 3.0 * g_xxyyyy_0_xxxxzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxxzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyyyy_0[i] = g_xxyyyy_0_xxxyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyyyz_0[i] = g_xxyyyy_0_xxxyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyyzz_0[i] = 2.0 * g_xxyyyy_0_xxxyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxyzzz_0[i] = 3.0 * g_xxyyyy_0_xxxyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxxzzzz_0[i] = 4.0 * g_xxyyyy_0_xxxzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxxzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyyyy_0[i] = g_xxyyyy_0_xxyyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyyyz_0[i] = g_xxyyyy_0_xxyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyyzz_0[i] = 2.0 * g_xxyyyy_0_xxyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyyzzz_0[i] = 3.0 * g_xxyyyy_0_xxyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxyzzzz_0[i] = 4.0 * g_xxyyyy_0_xxyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxyzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xxzzzzz_0[i] = 5.0 * g_xxyyyy_0_xxzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xxzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyyyy_0[i] = g_xxyyyy_0_xyyyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyyyz_0[i] = g_xxyyyy_0_xyyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyyzz_0[i] = 2.0 * g_xxyyyy_0_xyyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyyzzz_0[i] = 3.0 * g_xxyyyy_0_xyyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyyzzzz_0[i] = 4.0 * g_xxyyyy_0_xyyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyyzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xyzzzzz_0[i] = 5.0 * g_xxyyyy_0_xyzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xyzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_xzzzzzz_0[i] = 6.0 * g_xxyyyy_0_xzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_xzzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyyyy_0[i] = g_xxyyyy_0_yyyyyyy_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyyyz_0[i] = g_xxyyyy_0_yyyyyy_1[i] * fi_acd_0 + g_xxyyyy_0_yyyyyyz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyyzz_0[i] = 2.0 * g_xxyyyy_0_yyyyyz_1[i] * fi_acd_0 + g_xxyyyy_0_yyyyyzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyyzzz_0[i] = 3.0 * g_xxyyyy_0_yyyyzz_1[i] * fi_acd_0 + g_xxyyyy_0_yyyyzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyyzzzz_0[i] = 4.0 * g_xxyyyy_0_yyyzzz_1[i] * fi_acd_0 + g_xxyyyy_0_yyyzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yyzzzzz_0[i] = 5.0 * g_xxyyyy_0_yyzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_yyzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_yzzzzzz_0[i] = 6.0 * g_xxyyyy_0_yzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_yzzzzzz_1[i] * wa_z[i];

        g_xxyyyyz_0_zzzzzzz_0[i] = 7.0 * g_xxyyyy_0_zzzzzz_1[i] * fi_acd_0 + g_xxyyyy_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 612-648 components of targeted buffer : KSK

    auto g_xxyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 612);

    auto g_xxyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 613);

    auto g_xxyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 614);

    auto g_xxyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 615);

    auto g_xxyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 616);

    auto g_xxyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 617);

    auto g_xxyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 618);

    auto g_xxyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 619);

    auto g_xxyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 620);

    auto g_xxyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 621);

    auto g_xxyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 622);

    auto g_xxyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 623);

    auto g_xxyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 624);

    auto g_xxyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 625);

    auto g_xxyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 626);

    auto g_xxyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 627);

    auto g_xxyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 628);

    auto g_xxyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 629);

    auto g_xxyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 630);

    auto g_xxyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 631);

    auto g_xxyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 632);

    auto g_xxyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 633);

    auto g_xxyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 634);

    auto g_xxyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 635);

    auto g_xxyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 636);

    auto g_xxyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 637);

    auto g_xxyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 638);

    auto g_xxyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 639);

    auto g_xxyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 640);

    auto g_xxyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 641);

    auto g_xxyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 642);

    auto g_xxyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 643);

    auto g_xxyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 644);

    auto g_xxyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 645);

    auto g_xxyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 646);

    auto g_xxyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 647);

    #pragma omp simd aligned(g_xxyyy_0_xxxxxxy_0, g_xxyyy_0_xxxxxxy_1, g_xxyyy_0_xxxxxyy_0, g_xxyyy_0_xxxxxyy_1, g_xxyyy_0_xxxxyyy_0, g_xxyyy_0_xxxxyyy_1, g_xxyyy_0_xxxyyyy_0, g_xxyyy_0_xxxyyyy_1, g_xxyyy_0_xxyyyyy_0, g_xxyyy_0_xxyyyyy_1, g_xxyyy_0_xyyyyyy_0, g_xxyyy_0_xyyyyyy_1, g_xxyyyz_0_xxxxxxy_1, g_xxyyyz_0_xxxxxyy_1, g_xxyyyz_0_xxxxyyy_1, g_xxyyyz_0_xxxyyyy_1, g_xxyyyz_0_xxyyyyy_1, g_xxyyyz_0_xyyyyyy_1, g_xxyyyzz_0_xxxxxxx_0, g_xxyyyzz_0_xxxxxxy_0, g_xxyyyzz_0_xxxxxxz_0, g_xxyyyzz_0_xxxxxyy_0, g_xxyyyzz_0_xxxxxyz_0, g_xxyyyzz_0_xxxxxzz_0, g_xxyyyzz_0_xxxxyyy_0, g_xxyyyzz_0_xxxxyyz_0, g_xxyyyzz_0_xxxxyzz_0, g_xxyyyzz_0_xxxxzzz_0, g_xxyyyzz_0_xxxyyyy_0, g_xxyyyzz_0_xxxyyyz_0, g_xxyyyzz_0_xxxyyzz_0, g_xxyyyzz_0_xxxyzzz_0, g_xxyyyzz_0_xxxzzzz_0, g_xxyyyzz_0_xxyyyyy_0, g_xxyyyzz_0_xxyyyyz_0, g_xxyyyzz_0_xxyyyzz_0, g_xxyyyzz_0_xxyyzzz_0, g_xxyyyzz_0_xxyzzzz_0, g_xxyyyzz_0_xxzzzzz_0, g_xxyyyzz_0_xyyyyyy_0, g_xxyyyzz_0_xyyyyyz_0, g_xxyyyzz_0_xyyyyzz_0, g_xxyyyzz_0_xyyyzzz_0, g_xxyyyzz_0_xyyzzzz_0, g_xxyyyzz_0_xyzzzzz_0, g_xxyyyzz_0_xzzzzzz_0, g_xxyyyzz_0_yyyyyyy_0, g_xxyyyzz_0_yyyyyyz_0, g_xxyyyzz_0_yyyyyzz_0, g_xxyyyzz_0_yyyyzzz_0, g_xxyyyzz_0_yyyzzzz_0, g_xxyyyzz_0_yyzzzzz_0, g_xxyyyzz_0_yzzzzzz_0, g_xxyyyzz_0_zzzzzzz_0, g_xxyyzz_0_xxxxxxx_1, g_xxyyzz_0_xxxxxxz_1, g_xxyyzz_0_xxxxxzz_1, g_xxyyzz_0_xxxxzzz_1, g_xxyyzz_0_xxxzzzz_1, g_xxyyzz_0_xxzzzzz_1, g_xxyyzz_0_xzzzzzz_1, g_xxyzz_0_xxxxxxx_0, g_xxyzz_0_xxxxxxx_1, g_xxyzz_0_xxxxxxz_0, g_xxyzz_0_xxxxxxz_1, g_xxyzz_0_xxxxxzz_0, g_xxyzz_0_xxxxxzz_1, g_xxyzz_0_xxxxzzz_0, g_xxyzz_0_xxxxzzz_1, g_xxyzz_0_xxxzzzz_0, g_xxyzz_0_xxxzzzz_1, g_xxyzz_0_xxzzzzz_0, g_xxyzz_0_xxzzzzz_1, g_xxyzz_0_xzzzzzz_0, g_xxyzz_0_xzzzzzz_1, g_xyyyzz_0_xxxxxyz_1, g_xyyyzz_0_xxxxyyz_1, g_xyyyzz_0_xxxxyz_1, g_xyyyzz_0_xxxxyzz_1, g_xyyyzz_0_xxxyyyz_1, g_xyyyzz_0_xxxyyz_1, g_xyyyzz_0_xxxyyzz_1, g_xyyyzz_0_xxxyzz_1, g_xyyyzz_0_xxxyzzz_1, g_xyyyzz_0_xxyyyyz_1, g_xyyyzz_0_xxyyyz_1, g_xyyyzz_0_xxyyyzz_1, g_xyyyzz_0_xxyyzz_1, g_xyyyzz_0_xxyyzzz_1, g_xyyyzz_0_xxyzzz_1, g_xyyyzz_0_xxyzzzz_1, g_xyyyzz_0_xyyyyyz_1, g_xyyyzz_0_xyyyyz_1, g_xyyyzz_0_xyyyyzz_1, g_xyyyzz_0_xyyyzz_1, g_xyyyzz_0_xyyyzzz_1, g_xyyyzz_0_xyyzzz_1, g_xyyyzz_0_xyyzzzz_1, g_xyyyzz_0_xyzzzz_1, g_xyyyzz_0_xyzzzzz_1, g_xyyyzz_0_yyyyyyy_1, g_xyyyzz_0_yyyyyyz_1, g_xyyyzz_0_yyyyyz_1, g_xyyyzz_0_yyyyyzz_1, g_xyyyzz_0_yyyyzz_1, g_xyyyzz_0_yyyyzzz_1, g_xyyyzz_0_yyyzzz_1, g_xyyyzz_0_yyyzzzz_1, g_xyyyzz_0_yyzzzz_1, g_xyyyzz_0_yyzzzzz_1, g_xyyyzz_0_yzzzzz_1, g_xyyyzz_0_yzzzzzz_1, g_xyyyzz_0_zzzzzzz_1, g_yyyzz_0_xxxxxyz_0, g_yyyzz_0_xxxxxyz_1, g_yyyzz_0_xxxxyyz_0, g_yyyzz_0_xxxxyyz_1, g_yyyzz_0_xxxxyzz_0, g_yyyzz_0_xxxxyzz_1, g_yyyzz_0_xxxyyyz_0, g_yyyzz_0_xxxyyyz_1, g_yyyzz_0_xxxyyzz_0, g_yyyzz_0_xxxyyzz_1, g_yyyzz_0_xxxyzzz_0, g_yyyzz_0_xxxyzzz_1, g_yyyzz_0_xxyyyyz_0, g_yyyzz_0_xxyyyyz_1, g_yyyzz_0_xxyyyzz_0, g_yyyzz_0_xxyyyzz_1, g_yyyzz_0_xxyyzzz_0, g_yyyzz_0_xxyyzzz_1, g_yyyzz_0_xxyzzzz_0, g_yyyzz_0_xxyzzzz_1, g_yyyzz_0_xyyyyyz_0, g_yyyzz_0_xyyyyyz_1, g_yyyzz_0_xyyyyzz_0, g_yyyzz_0_xyyyyzz_1, g_yyyzz_0_xyyyzzz_0, g_yyyzz_0_xyyyzzz_1, g_yyyzz_0_xyyzzzz_0, g_yyyzz_0_xyyzzzz_1, g_yyyzz_0_xyzzzzz_0, g_yyyzz_0_xyzzzzz_1, g_yyyzz_0_yyyyyyy_0, g_yyyzz_0_yyyyyyy_1, g_yyyzz_0_yyyyyyz_0, g_yyyzz_0_yyyyyyz_1, g_yyyzz_0_yyyyyzz_0, g_yyyzz_0_yyyyyzz_1, g_yyyzz_0_yyyyzzz_0, g_yyyzz_0_yyyyzzz_1, g_yyyzz_0_yyyzzzz_0, g_yyyzz_0_yyyzzzz_1, g_yyyzz_0_yyzzzzz_0, g_yyyzz_0_yyzzzzz_1, g_yyyzz_0_yzzzzzz_0, g_yyyzz_0_yzzzzzz_1, g_yyyzz_0_zzzzzzz_0, g_yyyzz_0_zzzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyyzz_0_xxxxxxx_0[i] = 2.0 * g_xxyzz_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxxxx_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxxxxy_0[i] = g_xxyyy_0_xxxxxxy_0[i] * fbe_0 - g_xxyyy_0_xxxxxxy_1[i] * fz_be_0 + g_xxyyyz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxxxxz_0[i] = 2.0 * g_xxyzz_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxxxz_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxxxyy_0[i] = g_xxyyy_0_xxxxxyy_0[i] * fbe_0 - g_xxyyy_0_xxxxxyy_1[i] * fz_be_0 + g_xxyyyz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxxxyz_0[i] = g_yyyzz_0_xxxxxyz_0[i] * fbe_0 - g_yyyzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xyyyzz_0_xxxxyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxxxzz_0[i] = 2.0 * g_xxyzz_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxxzz_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxxyyy_0[i] = g_xxyyy_0_xxxxyyy_0[i] * fbe_0 - g_xxyyy_0_xxxxyyy_1[i] * fz_be_0 + g_xxyyyz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxxyyz_0[i] = g_yyyzz_0_xxxxyyz_0[i] * fbe_0 - g_yyyzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xyyyzz_0_xxxyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxxyzz_0[i] = g_yyyzz_0_xxxxyzz_0[i] * fbe_0 - g_yyyzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xyyyzz_0_xxxyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxxzzz_0[i] = 2.0 * g_xxyzz_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxxzzz_1[i] * fz_be_0 + g_xxyyzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxxyyyy_0[i] = g_xxyyy_0_xxxyyyy_0[i] * fbe_0 - g_xxyyy_0_xxxyyyy_1[i] * fz_be_0 + g_xxyyyz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxxyyyz_0[i] = g_yyyzz_0_xxxyyyz_0[i] * fbe_0 - g_yyyzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_0_xxyyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxyyzz_0[i] = g_yyyzz_0_xxxyyzz_0[i] * fbe_0 - g_yyyzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_0_xxyyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxyzzz_0[i] = g_yyyzz_0_xxxyzzz_0[i] * fbe_0 - g_yyyzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xyyyzz_0_xxyzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxxzzzz_0[i] = 2.0 * g_xxyzz_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxxzzzz_1[i] * fz_be_0 + g_xxyyzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xxyyyyy_0[i] = g_xxyyy_0_xxyyyyy_0[i] * fbe_0 - g_xxyyy_0_xxyyyyy_1[i] * fz_be_0 + g_xxyyyz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xxyyyyz_0[i] = g_yyyzz_0_xxyyyyz_0[i] * fbe_0 - g_yyyzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyyyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxyyyzz_0[i] = g_yyyzz_0_xxyyyzz_0[i] * fbe_0 - g_yyyzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyyyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxyyzzz_0[i] = g_yyyzz_0_xxyyzzz_0[i] * fbe_0 - g_yyyzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyyzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxyzzzz_0[i] = g_yyyzz_0_xxyzzzz_0[i] * fbe_0 - g_yyyzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xyyyzz_0_xyzzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xxzzzzz_0[i] = 2.0 * g_xxyzz_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xxzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_xyyyyyy_0[i] = g_xxyyy_0_xyyyyyy_0[i] * fbe_0 - g_xxyyy_0_xyyyyyy_1[i] * fz_be_0 + g_xxyyyz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxyyyzz_0_xyyyyyz_0[i] = g_yyyzz_0_xyyyyyz_0[i] * fbe_0 - g_yyyzz_0_xyyyyyz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyyz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyyyyzz_0[i] = g_yyyzz_0_xyyyyzz_0[i] * fbe_0 - g_yyyzz_0_xyyyyzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyyyzzz_0[i] = g_yyyzz_0_xyyyzzz_0[i] * fbe_0 - g_yyyzz_0_xyyyzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyyzzzz_0[i] = g_yyyzz_0_xyyzzzz_0[i] * fbe_0 - g_yyyzz_0_xyyzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyzzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xyzzzzz_0[i] = g_yyyzz_0_xyzzzzz_0[i] * fbe_0 - g_yyyzz_0_xyzzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yzzzzz_1[i] * fi_acd_0 + g_xyyyzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_xzzzzzz_0[i] = 2.0 * g_xxyzz_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_xxyzz_0_xzzzzzz_1[i] * fz_be_0 + g_xxyyzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxyyyzz_0_yyyyyyy_0[i] = g_yyyzz_0_yyyyyyy_0[i] * fbe_0 - g_yyyzz_0_yyyyyyy_1[i] * fz_be_0 + g_xyyyzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyyyyz_0[i] = g_yyyzz_0_yyyyyyz_0[i] * fbe_0 - g_yyyzz_0_yyyyyyz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyyyzz_0[i] = g_yyyzz_0_yyyyyzz_0[i] * fbe_0 - g_yyyzz_0_yyyyyzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyyzzz_0[i] = g_yyyzz_0_yyyyzzz_0[i] * fbe_0 - g_yyyzz_0_yyyyzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyyzzzz_0[i] = g_yyyzz_0_yyyzzzz_0[i] * fbe_0 - g_yyyzz_0_yyyzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yyzzzzz_0[i] = g_yyyzz_0_yyzzzzz_0[i] * fbe_0 - g_yyyzz_0_yyzzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_yzzzzzz_0[i] = g_yyyzz_0_yzzzzzz_0[i] * fbe_0 - g_yyyzz_0_yzzzzzz_1[i] * fz_be_0 + g_xyyyzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxyyyzz_0_zzzzzzz_0[i] = g_yyyzz_0_zzzzzzz_0[i] * fbe_0 - g_yyyzz_0_zzzzzzz_1[i] * fz_be_0 + g_xyyyzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 648-684 components of targeted buffer : KSK

    auto g_xxyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 648);

    auto g_xxyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 649);

    auto g_xxyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 650);

    auto g_xxyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 651);

    auto g_xxyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 652);

    auto g_xxyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 653);

    auto g_xxyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 654);

    auto g_xxyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 655);

    auto g_xxyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 656);

    auto g_xxyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 657);

    auto g_xxyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 658);

    auto g_xxyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 659);

    auto g_xxyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 660);

    auto g_xxyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 661);

    auto g_xxyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 662);

    auto g_xxyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 663);

    auto g_xxyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 664);

    auto g_xxyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 665);

    auto g_xxyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 666);

    auto g_xxyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 667);

    auto g_xxyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 668);

    auto g_xxyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 669);

    auto g_xxyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 670);

    auto g_xxyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 671);

    auto g_xxyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 672);

    auto g_xxyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 673);

    auto g_xxyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 674);

    auto g_xxyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 675);

    auto g_xxyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 676);

    auto g_xxyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 677);

    auto g_xxyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 678);

    auto g_xxyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 679);

    auto g_xxyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 680);

    auto g_xxyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 681);

    auto g_xxyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 682);

    auto g_xxyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 683);

    #pragma omp simd aligned(g_xxyyz_0_xxxxxxy_0, g_xxyyz_0_xxxxxxy_1, g_xxyyz_0_xxxxxyy_0, g_xxyyz_0_xxxxxyy_1, g_xxyyz_0_xxxxyyy_0, g_xxyyz_0_xxxxyyy_1, g_xxyyz_0_xxxyyyy_0, g_xxyyz_0_xxxyyyy_1, g_xxyyz_0_xxyyyyy_0, g_xxyyz_0_xxyyyyy_1, g_xxyyz_0_xyyyyyy_0, g_xxyyz_0_xyyyyyy_1, g_xxyyzz_0_xxxxxxy_1, g_xxyyzz_0_xxxxxyy_1, g_xxyyzz_0_xxxxyyy_1, g_xxyyzz_0_xxxyyyy_1, g_xxyyzz_0_xxyyyyy_1, g_xxyyzz_0_xyyyyyy_1, g_xxyyzzz_0_xxxxxxx_0, g_xxyyzzz_0_xxxxxxy_0, g_xxyyzzz_0_xxxxxxz_0, g_xxyyzzz_0_xxxxxyy_0, g_xxyyzzz_0_xxxxxyz_0, g_xxyyzzz_0_xxxxxzz_0, g_xxyyzzz_0_xxxxyyy_0, g_xxyyzzz_0_xxxxyyz_0, g_xxyyzzz_0_xxxxyzz_0, g_xxyyzzz_0_xxxxzzz_0, g_xxyyzzz_0_xxxyyyy_0, g_xxyyzzz_0_xxxyyyz_0, g_xxyyzzz_0_xxxyyzz_0, g_xxyyzzz_0_xxxyzzz_0, g_xxyyzzz_0_xxxzzzz_0, g_xxyyzzz_0_xxyyyyy_0, g_xxyyzzz_0_xxyyyyz_0, g_xxyyzzz_0_xxyyyzz_0, g_xxyyzzz_0_xxyyzzz_0, g_xxyyzzz_0_xxyzzzz_0, g_xxyyzzz_0_xxzzzzz_0, g_xxyyzzz_0_xyyyyyy_0, g_xxyyzzz_0_xyyyyyz_0, g_xxyyzzz_0_xyyyyzz_0, g_xxyyzzz_0_xyyyzzz_0, g_xxyyzzz_0_xyyzzzz_0, g_xxyyzzz_0_xyzzzzz_0, g_xxyyzzz_0_xzzzzzz_0, g_xxyyzzz_0_yyyyyyy_0, g_xxyyzzz_0_yyyyyyz_0, g_xxyyzzz_0_yyyyyzz_0, g_xxyyzzz_0_yyyyzzz_0, g_xxyyzzz_0_yyyzzzz_0, g_xxyyzzz_0_yyzzzzz_0, g_xxyyzzz_0_yzzzzzz_0, g_xxyyzzz_0_zzzzzzz_0, g_xxyzzz_0_xxxxxxx_1, g_xxyzzz_0_xxxxxxz_1, g_xxyzzz_0_xxxxxzz_1, g_xxyzzz_0_xxxxzzz_1, g_xxyzzz_0_xxxzzzz_1, g_xxyzzz_0_xxzzzzz_1, g_xxyzzz_0_xzzzzzz_1, g_xxzzz_0_xxxxxxx_0, g_xxzzz_0_xxxxxxx_1, g_xxzzz_0_xxxxxxz_0, g_xxzzz_0_xxxxxxz_1, g_xxzzz_0_xxxxxzz_0, g_xxzzz_0_xxxxxzz_1, g_xxzzz_0_xxxxzzz_0, g_xxzzz_0_xxxxzzz_1, g_xxzzz_0_xxxzzzz_0, g_xxzzz_0_xxxzzzz_1, g_xxzzz_0_xxzzzzz_0, g_xxzzz_0_xxzzzzz_1, g_xxzzz_0_xzzzzzz_0, g_xxzzz_0_xzzzzzz_1, g_xyyzzz_0_xxxxxyz_1, g_xyyzzz_0_xxxxyyz_1, g_xyyzzz_0_xxxxyz_1, g_xyyzzz_0_xxxxyzz_1, g_xyyzzz_0_xxxyyyz_1, g_xyyzzz_0_xxxyyz_1, g_xyyzzz_0_xxxyyzz_1, g_xyyzzz_0_xxxyzz_1, g_xyyzzz_0_xxxyzzz_1, g_xyyzzz_0_xxyyyyz_1, g_xyyzzz_0_xxyyyz_1, g_xyyzzz_0_xxyyyzz_1, g_xyyzzz_0_xxyyzz_1, g_xyyzzz_0_xxyyzzz_1, g_xyyzzz_0_xxyzzz_1, g_xyyzzz_0_xxyzzzz_1, g_xyyzzz_0_xyyyyyz_1, g_xyyzzz_0_xyyyyz_1, g_xyyzzz_0_xyyyyzz_1, g_xyyzzz_0_xyyyzz_1, g_xyyzzz_0_xyyyzzz_1, g_xyyzzz_0_xyyzzz_1, g_xyyzzz_0_xyyzzzz_1, g_xyyzzz_0_xyzzzz_1, g_xyyzzz_0_xyzzzzz_1, g_xyyzzz_0_yyyyyyy_1, g_xyyzzz_0_yyyyyyz_1, g_xyyzzz_0_yyyyyz_1, g_xyyzzz_0_yyyyyzz_1, g_xyyzzz_0_yyyyzz_1, g_xyyzzz_0_yyyyzzz_1, g_xyyzzz_0_yyyzzz_1, g_xyyzzz_0_yyyzzzz_1, g_xyyzzz_0_yyzzzz_1, g_xyyzzz_0_yyzzzzz_1, g_xyyzzz_0_yzzzzz_1, g_xyyzzz_0_yzzzzzz_1, g_xyyzzz_0_zzzzzzz_1, g_yyzzz_0_xxxxxyz_0, g_yyzzz_0_xxxxxyz_1, g_yyzzz_0_xxxxyyz_0, g_yyzzz_0_xxxxyyz_1, g_yyzzz_0_xxxxyzz_0, g_yyzzz_0_xxxxyzz_1, g_yyzzz_0_xxxyyyz_0, g_yyzzz_0_xxxyyyz_1, g_yyzzz_0_xxxyyzz_0, g_yyzzz_0_xxxyyzz_1, g_yyzzz_0_xxxyzzz_0, g_yyzzz_0_xxxyzzz_1, g_yyzzz_0_xxyyyyz_0, g_yyzzz_0_xxyyyyz_1, g_yyzzz_0_xxyyyzz_0, g_yyzzz_0_xxyyyzz_1, g_yyzzz_0_xxyyzzz_0, g_yyzzz_0_xxyyzzz_1, g_yyzzz_0_xxyzzzz_0, g_yyzzz_0_xxyzzzz_1, g_yyzzz_0_xyyyyyz_0, g_yyzzz_0_xyyyyyz_1, g_yyzzz_0_xyyyyzz_0, g_yyzzz_0_xyyyyzz_1, g_yyzzz_0_xyyyzzz_0, g_yyzzz_0_xyyyzzz_1, g_yyzzz_0_xyyzzzz_0, g_yyzzz_0_xyyzzzz_1, g_yyzzz_0_xyzzzzz_0, g_yyzzz_0_xyzzzzz_1, g_yyzzz_0_yyyyyyy_0, g_yyzzz_0_yyyyyyy_1, g_yyzzz_0_yyyyyyz_0, g_yyzzz_0_yyyyyyz_1, g_yyzzz_0_yyyyyzz_0, g_yyzzz_0_yyyyyzz_1, g_yyzzz_0_yyyyzzz_0, g_yyzzz_0_yyyyzzz_1, g_yyzzz_0_yyyzzzz_0, g_yyzzz_0_yyyzzzz_1, g_yyzzz_0_yyzzzzz_0, g_yyzzz_0_yyzzzzz_1, g_yyzzz_0_yzzzzzz_0, g_yyzzz_0_yzzzzzz_1, g_yyzzz_0_zzzzzzz_0, g_yyzzz_0_zzzzzzz_1, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyzzz_0_xxxxxxx_0[i] = g_xxzzz_0_xxxxxxx_0[i] * fbe_0 - g_xxzzz_0_xxxxxxx_1[i] * fz_be_0 + g_xxyzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxxxxy_0[i] = 2.0 * g_xxyyz_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxxxxy_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxxxxz_0[i] = g_xxzzz_0_xxxxxxz_0[i] * fbe_0 - g_xxzzz_0_xxxxxxz_1[i] * fz_be_0 + g_xxyzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxxxyy_0[i] = 2.0 * g_xxyyz_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxxxyy_1[i] * fz_be_0 + g_xxyyzz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxxxyz_0[i] = g_yyzzz_0_xxxxxyz_0[i] * fbe_0 - g_yyzzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xyyzzz_0_xxxxyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxxxzz_0[i] = g_xxzzz_0_xxxxxzz_0[i] * fbe_0 - g_xxzzz_0_xxxxxzz_1[i] * fz_be_0 + g_xxyzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxxyyy_0[i] = 2.0 * g_xxyyz_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxxyyy_1[i] * fz_be_0 + g_xxyyzz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxxyyz_0[i] = g_yyzzz_0_xxxxyyz_0[i] * fbe_0 - g_yyzzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xyyzzz_0_xxxyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxxyzz_0[i] = g_yyzzz_0_xxxxyzz_0[i] * fbe_0 - g_yyzzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xyyzzz_0_xxxyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxxzzz_0[i] = g_xxzzz_0_xxxxzzz_0[i] * fbe_0 - g_xxzzz_0_xxxxzzz_1[i] * fz_be_0 + g_xxyzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxxyyyy_0[i] = 2.0 * g_xxyyz_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxxyyyy_1[i] * fz_be_0 + g_xxyyzz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxxyyyz_0[i] = g_yyzzz_0_xxxyyyz_0[i] * fbe_0 - g_yyzzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_0_xxyyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxyyzz_0[i] = g_yyzzz_0_xxxyyzz_0[i] * fbe_0 - g_yyzzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_0_xxyyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxyzzz_0[i] = g_yyzzz_0_xxxyzzz_0[i] * fbe_0 - g_yyzzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xyyzzz_0_xxyzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxxzzzz_0[i] = g_xxzzz_0_xxxzzzz_0[i] * fbe_0 - g_xxzzz_0_xxxzzzz_1[i] * fz_be_0 + g_xxyzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xxyyyyy_0[i] = 2.0 * g_xxyyz_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xxyyyyy_1[i] * fz_be_0 + g_xxyyzz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xxyyyyz_0[i] = g_yyzzz_0_xxyyyyz_0[i] * fbe_0 - g_yyzzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyyyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxyyyzz_0[i] = g_yyzzz_0_xxyyyzz_0[i] * fbe_0 - g_yyzzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyyyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxyyzzz_0[i] = g_yyzzz_0_xxyyzzz_0[i] * fbe_0 - g_yyzzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyyzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxyzzzz_0[i] = g_yyzzz_0_xxyzzzz_0[i] * fbe_0 - g_yyzzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xyyzzz_0_xyzzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xxzzzzz_0[i] = g_xxzzz_0_xxzzzzz_0[i] * fbe_0 - g_xxzzz_0_xxzzzzz_1[i] * fz_be_0 + g_xxyzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_xyyyyyy_0[i] = 2.0 * g_xxyyz_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_xxyyz_0_xyyyyyy_1[i] * fz_be_0 + g_xxyyzz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxyyzzz_0_xyyyyyz_0[i] = g_yyzzz_0_xyyyyyz_0[i] * fbe_0 - g_yyzzz_0_xyyyyyz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyyz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyyyyzz_0[i] = g_yyzzz_0_xyyyyzz_0[i] * fbe_0 - g_yyzzz_0_xyyyyzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyyyzzz_0[i] = g_yyzzz_0_xyyyzzz_0[i] * fbe_0 - g_yyzzz_0_xyyyzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyyzzzz_0[i] = g_yyzzz_0_xyyzzzz_0[i] * fbe_0 - g_yyzzz_0_xyyzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyzzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xyzzzzz_0[i] = g_yyzzz_0_xyzzzzz_0[i] * fbe_0 - g_yyzzz_0_xyzzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yzzzzz_1[i] * fi_acd_0 + g_xyyzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_xzzzzzz_0[i] = g_xxzzz_0_xzzzzzz_0[i] * fbe_0 - g_xxzzz_0_xzzzzzz_1[i] * fz_be_0 + g_xxyzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxyyzzz_0_yyyyyyy_0[i] = g_yyzzz_0_yyyyyyy_0[i] * fbe_0 - g_yyzzz_0_yyyyyyy_1[i] * fz_be_0 + g_xyyzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyyyyz_0[i] = g_yyzzz_0_yyyyyyz_0[i] * fbe_0 - g_yyzzz_0_yyyyyyz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyyyzz_0[i] = g_yyzzz_0_yyyyyzz_0[i] * fbe_0 - g_yyzzz_0_yyyyyzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyyzzz_0[i] = g_yyzzz_0_yyyyzzz_0[i] * fbe_0 - g_yyzzz_0_yyyyzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyyzzzz_0[i] = g_yyzzz_0_yyyzzzz_0[i] * fbe_0 - g_yyzzz_0_yyyzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yyzzzzz_0[i] = g_yyzzz_0_yyzzzzz_0[i] * fbe_0 - g_yyzzz_0_yyzzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_yzzzzzz_0[i] = g_yyzzz_0_yzzzzzz_0[i] * fbe_0 - g_yyzzz_0_yzzzzzz_1[i] * fz_be_0 + g_xyyzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxyyzzz_0_zzzzzzz_0[i] = g_yyzzz_0_zzzzzzz_0[i] * fbe_0 - g_yyzzz_0_zzzzzzz_1[i] * fz_be_0 + g_xyyzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 684-720 components of targeted buffer : KSK

    auto g_xxyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 684);

    auto g_xxyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 685);

    auto g_xxyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 686);

    auto g_xxyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 687);

    auto g_xxyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 688);

    auto g_xxyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 689);

    auto g_xxyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 690);

    auto g_xxyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 691);

    auto g_xxyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 692);

    auto g_xxyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 693);

    auto g_xxyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 694);

    auto g_xxyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 695);

    auto g_xxyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 696);

    auto g_xxyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 697);

    auto g_xxyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 698);

    auto g_xxyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 699);

    auto g_xxyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 700);

    auto g_xxyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 701);

    auto g_xxyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 702);

    auto g_xxyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 703);

    auto g_xxyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 704);

    auto g_xxyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 705);

    auto g_xxyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 706);

    auto g_xxyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 707);

    auto g_xxyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 708);

    auto g_xxyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 709);

    auto g_xxyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 710);

    auto g_xxyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 711);

    auto g_xxyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 712);

    auto g_xxyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 713);

    auto g_xxyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 714);

    auto g_xxyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 715);

    auto g_xxyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 716);

    auto g_xxyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 717);

    auto g_xxyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 718);

    auto g_xxyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 719);

    #pragma omp simd aligned(g_xxyzzzz_0_xxxxxxx_0, g_xxyzzzz_0_xxxxxxy_0, g_xxyzzzz_0_xxxxxxz_0, g_xxyzzzz_0_xxxxxyy_0, g_xxyzzzz_0_xxxxxyz_0, g_xxyzzzz_0_xxxxxzz_0, g_xxyzzzz_0_xxxxyyy_0, g_xxyzzzz_0_xxxxyyz_0, g_xxyzzzz_0_xxxxyzz_0, g_xxyzzzz_0_xxxxzzz_0, g_xxyzzzz_0_xxxyyyy_0, g_xxyzzzz_0_xxxyyyz_0, g_xxyzzzz_0_xxxyyzz_0, g_xxyzzzz_0_xxxyzzz_0, g_xxyzzzz_0_xxxzzzz_0, g_xxyzzzz_0_xxyyyyy_0, g_xxyzzzz_0_xxyyyyz_0, g_xxyzzzz_0_xxyyyzz_0, g_xxyzzzz_0_xxyyzzz_0, g_xxyzzzz_0_xxyzzzz_0, g_xxyzzzz_0_xxzzzzz_0, g_xxyzzzz_0_xyyyyyy_0, g_xxyzzzz_0_xyyyyyz_0, g_xxyzzzz_0_xyyyyzz_0, g_xxyzzzz_0_xyyyzzz_0, g_xxyzzzz_0_xyyzzzz_0, g_xxyzzzz_0_xyzzzzz_0, g_xxyzzzz_0_xzzzzzz_0, g_xxyzzzz_0_yyyyyyy_0, g_xxyzzzz_0_yyyyyyz_0, g_xxyzzzz_0_yyyyyzz_0, g_xxyzzzz_0_yyyyzzz_0, g_xxyzzzz_0_yyyzzzz_0, g_xxyzzzz_0_yyzzzzz_0, g_xxyzzzz_0_yzzzzzz_0, g_xxyzzzz_0_zzzzzzz_0, g_xxzzzz_0_xxxxxx_1, g_xxzzzz_0_xxxxxxx_1, g_xxzzzz_0_xxxxxxy_1, g_xxzzzz_0_xxxxxxz_1, g_xxzzzz_0_xxxxxy_1, g_xxzzzz_0_xxxxxyy_1, g_xxzzzz_0_xxxxxyz_1, g_xxzzzz_0_xxxxxz_1, g_xxzzzz_0_xxxxxzz_1, g_xxzzzz_0_xxxxyy_1, g_xxzzzz_0_xxxxyyy_1, g_xxzzzz_0_xxxxyyz_1, g_xxzzzz_0_xxxxyz_1, g_xxzzzz_0_xxxxyzz_1, g_xxzzzz_0_xxxxzz_1, g_xxzzzz_0_xxxxzzz_1, g_xxzzzz_0_xxxyyy_1, g_xxzzzz_0_xxxyyyy_1, g_xxzzzz_0_xxxyyyz_1, g_xxzzzz_0_xxxyyz_1, g_xxzzzz_0_xxxyyzz_1, g_xxzzzz_0_xxxyzz_1, g_xxzzzz_0_xxxyzzz_1, g_xxzzzz_0_xxxzzz_1, g_xxzzzz_0_xxxzzzz_1, g_xxzzzz_0_xxyyyy_1, g_xxzzzz_0_xxyyyyy_1, g_xxzzzz_0_xxyyyyz_1, g_xxzzzz_0_xxyyyz_1, g_xxzzzz_0_xxyyyzz_1, g_xxzzzz_0_xxyyzz_1, g_xxzzzz_0_xxyyzzz_1, g_xxzzzz_0_xxyzzz_1, g_xxzzzz_0_xxyzzzz_1, g_xxzzzz_0_xxzzzz_1, g_xxzzzz_0_xxzzzzz_1, g_xxzzzz_0_xyyyyy_1, g_xxzzzz_0_xyyyyyy_1, g_xxzzzz_0_xyyyyyz_1, g_xxzzzz_0_xyyyyz_1, g_xxzzzz_0_xyyyyzz_1, g_xxzzzz_0_xyyyzz_1, g_xxzzzz_0_xyyyzzz_1, g_xxzzzz_0_xyyzzz_1, g_xxzzzz_0_xyyzzzz_1, g_xxzzzz_0_xyzzzz_1, g_xxzzzz_0_xyzzzzz_1, g_xxzzzz_0_xzzzzz_1, g_xxzzzz_0_xzzzzzz_1, g_xxzzzz_0_yyyyyy_1, g_xxzzzz_0_yyyyyyy_1, g_xxzzzz_0_yyyyyyz_1, g_xxzzzz_0_yyyyyz_1, g_xxzzzz_0_yyyyyzz_1, g_xxzzzz_0_yyyyzz_1, g_xxzzzz_0_yyyyzzz_1, g_xxzzzz_0_yyyzzz_1, g_xxzzzz_0_yyyzzzz_1, g_xxzzzz_0_yyzzzz_1, g_xxzzzz_0_yyzzzzz_1, g_xxzzzz_0_yzzzzz_1, g_xxzzzz_0_yzzzzzz_1, g_xxzzzz_0_zzzzzz_1, g_xxzzzz_0_zzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzzz_0_xxxxxxx_0[i] = g_xxzzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxxy_0[i] = g_xxzzzz_0_xxxxxx_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxxy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxxz_0[i] = g_xxzzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxyy_0[i] = 2.0 * g_xxzzzz_0_xxxxxy_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxyz_0[i] = g_xxzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxxyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxxzz_0[i] = g_xxzzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxyyy_0[i] = 3.0 * g_xxzzzz_0_xxxxyy_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxyyz_0[i] = 2.0 * g_xxzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxyzz_0[i] = g_xxzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxxyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxxzzz_0[i] = g_xxzzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyyyy_0[i] = 4.0 * g_xxzzzz_0_xxxyyy_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyyyz_0[i] = 3.0 * g_xxzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyyzz_0[i] = 2.0 * g_xxzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxyzzz_0[i] = g_xxzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxxyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxxzzzz_0[i] = g_xxzzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyyyy_0[i] = 5.0 * g_xxzzzz_0_xxyyyy_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyyyz_0[i] = 4.0 * g_xxzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyyzz_0[i] = 3.0 * g_xxzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyyzzz_0[i] = 2.0 * g_xxzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxyzzzz_0[i] = g_xxzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xxyzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xxzzzzz_0[i] = g_xxzzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyyyy_0[i] = 6.0 * g_xxzzzz_0_xyyyyy_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyyyz_0[i] = 5.0 * g_xxzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyyzz_0[i] = 4.0 * g_xxzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyyzzz_0[i] = 3.0 * g_xxzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyyzzzz_0[i] = 2.0 * g_xxzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyyzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xyzzzzz_0[i] = g_xxzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_xyzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_xzzzzzz_0[i] = g_xxzzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyyyy_0[i] = 7.0 * g_xxzzzz_0_yyyyyy_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyyyy_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyyyz_0[i] = 6.0 * g_xxzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyyyz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyyzz_0[i] = 5.0 * g_xxzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyyzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyyzzz_0[i] = 4.0 * g_xxzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyyzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyyzzzz_0[i] = 3.0 * g_xxzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyyzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yyzzzzz_0[i] = 2.0 * g_xxzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yyzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_yzzzzzz_0[i] = g_xxzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_xxzzzz_0_yzzzzzz_1[i] * wa_y[i];

        g_xxyzzzz_0_zzzzzzz_0[i] = g_xxzzzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 720-756 components of targeted buffer : KSK

    auto g_xxzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 720);

    auto g_xxzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 721);

    auto g_xxzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 722);

    auto g_xxzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 723);

    auto g_xxzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 724);

    auto g_xxzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 725);

    auto g_xxzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 726);

    auto g_xxzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 727);

    auto g_xxzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 728);

    auto g_xxzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 729);

    auto g_xxzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 730);

    auto g_xxzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 731);

    auto g_xxzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 732);

    auto g_xxzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 733);

    auto g_xxzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 734);

    auto g_xxzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 735);

    auto g_xxzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 736);

    auto g_xxzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 737);

    auto g_xxzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 738);

    auto g_xxzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 739);

    auto g_xxzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 740);

    auto g_xxzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 741);

    auto g_xxzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 742);

    auto g_xxzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 743);

    auto g_xxzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 744);

    auto g_xxzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 745);

    auto g_xxzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 746);

    auto g_xxzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 747);

    auto g_xxzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 748);

    auto g_xxzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 749);

    auto g_xxzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 750);

    auto g_xxzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 751);

    auto g_xxzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 752);

    auto g_xxzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 753);

    auto g_xxzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 754);

    auto g_xxzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 755);

    #pragma omp simd aligned(g_xxzzz_0_xxxxxxx_0, g_xxzzz_0_xxxxxxx_1, g_xxzzz_0_xxxxxxy_0, g_xxzzz_0_xxxxxxy_1, g_xxzzz_0_xxxxxyy_0, g_xxzzz_0_xxxxxyy_1, g_xxzzz_0_xxxxyyy_0, g_xxzzz_0_xxxxyyy_1, g_xxzzz_0_xxxyyyy_0, g_xxzzz_0_xxxyyyy_1, g_xxzzz_0_xxyyyyy_0, g_xxzzz_0_xxyyyyy_1, g_xxzzz_0_xyyyyyy_0, g_xxzzz_0_xyyyyyy_1, g_xxzzzz_0_xxxxxxx_1, g_xxzzzz_0_xxxxxxy_1, g_xxzzzz_0_xxxxxyy_1, g_xxzzzz_0_xxxxyyy_1, g_xxzzzz_0_xxxyyyy_1, g_xxzzzz_0_xxyyyyy_1, g_xxzzzz_0_xyyyyyy_1, g_xxzzzzz_0_xxxxxxx_0, g_xxzzzzz_0_xxxxxxy_0, g_xxzzzzz_0_xxxxxxz_0, g_xxzzzzz_0_xxxxxyy_0, g_xxzzzzz_0_xxxxxyz_0, g_xxzzzzz_0_xxxxxzz_0, g_xxzzzzz_0_xxxxyyy_0, g_xxzzzzz_0_xxxxyyz_0, g_xxzzzzz_0_xxxxyzz_0, g_xxzzzzz_0_xxxxzzz_0, g_xxzzzzz_0_xxxyyyy_0, g_xxzzzzz_0_xxxyyyz_0, g_xxzzzzz_0_xxxyyzz_0, g_xxzzzzz_0_xxxyzzz_0, g_xxzzzzz_0_xxxzzzz_0, g_xxzzzzz_0_xxyyyyy_0, g_xxzzzzz_0_xxyyyyz_0, g_xxzzzzz_0_xxyyyzz_0, g_xxzzzzz_0_xxyyzzz_0, g_xxzzzzz_0_xxyzzzz_0, g_xxzzzzz_0_xxzzzzz_0, g_xxzzzzz_0_xyyyyyy_0, g_xxzzzzz_0_xyyyyyz_0, g_xxzzzzz_0_xyyyyzz_0, g_xxzzzzz_0_xyyyzzz_0, g_xxzzzzz_0_xyyzzzz_0, g_xxzzzzz_0_xyzzzzz_0, g_xxzzzzz_0_xzzzzzz_0, g_xxzzzzz_0_yyyyyyy_0, g_xxzzzzz_0_yyyyyyz_0, g_xxzzzzz_0_yyyyyzz_0, g_xxzzzzz_0_yyyyzzz_0, g_xxzzzzz_0_yyyzzzz_0, g_xxzzzzz_0_yyzzzzz_0, g_xxzzzzz_0_yzzzzzz_0, g_xxzzzzz_0_zzzzzzz_0, g_xzzzzz_0_xxxxxxz_1, g_xzzzzz_0_xxxxxyz_1, g_xzzzzz_0_xxxxxz_1, g_xzzzzz_0_xxxxxzz_1, g_xzzzzz_0_xxxxyyz_1, g_xzzzzz_0_xxxxyz_1, g_xzzzzz_0_xxxxyzz_1, g_xzzzzz_0_xxxxzz_1, g_xzzzzz_0_xxxxzzz_1, g_xzzzzz_0_xxxyyyz_1, g_xzzzzz_0_xxxyyz_1, g_xzzzzz_0_xxxyyzz_1, g_xzzzzz_0_xxxyzz_1, g_xzzzzz_0_xxxyzzz_1, g_xzzzzz_0_xxxzzz_1, g_xzzzzz_0_xxxzzzz_1, g_xzzzzz_0_xxyyyyz_1, g_xzzzzz_0_xxyyyz_1, g_xzzzzz_0_xxyyyzz_1, g_xzzzzz_0_xxyyzz_1, g_xzzzzz_0_xxyyzzz_1, g_xzzzzz_0_xxyzzz_1, g_xzzzzz_0_xxyzzzz_1, g_xzzzzz_0_xxzzzz_1, g_xzzzzz_0_xxzzzzz_1, g_xzzzzz_0_xyyyyyz_1, g_xzzzzz_0_xyyyyz_1, g_xzzzzz_0_xyyyyzz_1, g_xzzzzz_0_xyyyzz_1, g_xzzzzz_0_xyyyzzz_1, g_xzzzzz_0_xyyzzz_1, g_xzzzzz_0_xyyzzzz_1, g_xzzzzz_0_xyzzzz_1, g_xzzzzz_0_xyzzzzz_1, g_xzzzzz_0_xzzzzz_1, g_xzzzzz_0_xzzzzzz_1, g_xzzzzz_0_yyyyyyy_1, g_xzzzzz_0_yyyyyyz_1, g_xzzzzz_0_yyyyyz_1, g_xzzzzz_0_yyyyyzz_1, g_xzzzzz_0_yyyyzz_1, g_xzzzzz_0_yyyyzzz_1, g_xzzzzz_0_yyyzzz_1, g_xzzzzz_0_yyyzzzz_1, g_xzzzzz_0_yyzzzz_1, g_xzzzzz_0_yyzzzzz_1, g_xzzzzz_0_yzzzzz_1, g_xzzzzz_0_yzzzzzz_1, g_xzzzzz_0_zzzzzz_1, g_xzzzzz_0_zzzzzzz_1, g_zzzzz_0_xxxxxxz_0, g_zzzzz_0_xxxxxxz_1, g_zzzzz_0_xxxxxyz_0, g_zzzzz_0_xxxxxyz_1, g_zzzzz_0_xxxxxzz_0, g_zzzzz_0_xxxxxzz_1, g_zzzzz_0_xxxxyyz_0, g_zzzzz_0_xxxxyyz_1, g_zzzzz_0_xxxxyzz_0, g_zzzzz_0_xxxxyzz_1, g_zzzzz_0_xxxxzzz_0, g_zzzzz_0_xxxxzzz_1, g_zzzzz_0_xxxyyyz_0, g_zzzzz_0_xxxyyyz_1, g_zzzzz_0_xxxyyzz_0, g_zzzzz_0_xxxyyzz_1, g_zzzzz_0_xxxyzzz_0, g_zzzzz_0_xxxyzzz_1, g_zzzzz_0_xxxzzzz_0, g_zzzzz_0_xxxzzzz_1, g_zzzzz_0_xxyyyyz_0, g_zzzzz_0_xxyyyyz_1, g_zzzzz_0_xxyyyzz_0, g_zzzzz_0_xxyyyzz_1, g_zzzzz_0_xxyyzzz_0, g_zzzzz_0_xxyyzzz_1, g_zzzzz_0_xxyzzzz_0, g_zzzzz_0_xxyzzzz_1, g_zzzzz_0_xxzzzzz_0, g_zzzzz_0_xxzzzzz_1, g_zzzzz_0_xyyyyyz_0, g_zzzzz_0_xyyyyyz_1, g_zzzzz_0_xyyyyzz_0, g_zzzzz_0_xyyyyzz_1, g_zzzzz_0_xyyyzzz_0, g_zzzzz_0_xyyyzzz_1, g_zzzzz_0_xyyzzzz_0, g_zzzzz_0_xyyzzzz_1, g_zzzzz_0_xyzzzzz_0, g_zzzzz_0_xyzzzzz_1, g_zzzzz_0_xzzzzzz_0, g_zzzzz_0_xzzzzzz_1, g_zzzzz_0_yyyyyyy_0, g_zzzzz_0_yyyyyyy_1, g_zzzzz_0_yyyyyyz_0, g_zzzzz_0_yyyyyyz_1, g_zzzzz_0_yyyyyzz_0, g_zzzzz_0_yyyyyzz_1, g_zzzzz_0_yyyyzzz_0, g_zzzzz_0_yyyyzzz_1, g_zzzzz_0_yyyzzzz_0, g_zzzzz_0_yyyzzzz_1, g_zzzzz_0_yyzzzzz_0, g_zzzzz_0_yyzzzzz_1, g_zzzzz_0_yzzzzzz_0, g_zzzzz_0_yzzzzzz_1, g_zzzzz_0_zzzzzzz_0, g_zzzzz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzzzz_0_xxxxxxx_0[i] = 4.0 * g_xxzzz_0_xxxxxxx_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxxxx_1[i] * fz_be_0 + g_xxzzzz_0_xxxxxxx_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxxxy_0[i] = 4.0 * g_xxzzz_0_xxxxxxy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxxxy_1[i] * fz_be_0 + g_xxzzzz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxxxz_0[i] = g_zzzzz_0_xxxxxxz_0[i] * fbe_0 - g_zzzzz_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xzzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxxyy_0[i] = 4.0 * g_xxzzz_0_xxxxxyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxxyy_1[i] * fz_be_0 + g_xxzzzz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxxyz_0[i] = g_zzzzz_0_xxxxxyz_0[i] * fbe_0 - g_zzzzz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xzzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxxzz_0[i] = g_zzzzz_0_xxxxxzz_0[i] * fbe_0 - g_zzzzz_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xzzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxyyy_0[i] = 4.0 * g_xxzzz_0_xxxxyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxxyyy_1[i] * fz_be_0 + g_xxzzzz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxxyyz_0[i] = g_zzzzz_0_xxxxyyz_0[i] * fbe_0 - g_zzzzz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxyzz_0[i] = g_zzzzz_0_xxxxyzz_0[i] * fbe_0 - g_zzzzz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxxzzz_0[i] = g_zzzzz_0_xxxxzzz_0[i] * fbe_0 - g_zzzzz_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xzzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxyyyy_0[i] = 4.0 * g_xxzzz_0_xxxyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxxyyyy_1[i] * fz_be_0 + g_xxzzzz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxxyyyz_0[i] = g_zzzzz_0_xxxyyyz_0[i] * fbe_0 - g_zzzzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxyyzz_0[i] = g_zzzzz_0_xxxyyzz_0[i] * fbe_0 - g_zzzzz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxyzzz_0[i] = g_zzzzz_0_xxxyzzz_0[i] * fbe_0 - g_zzzzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxxzzzz_0[i] = g_zzzzz_0_xxxzzzz_0[i] * fbe_0 - g_zzzzz_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xzzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyyyyy_0[i] = 4.0 * g_xxzzz_0_xxyyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xxyyyyy_1[i] * fz_be_0 + g_xxzzzz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xxyyyyz_0[i] = g_zzzzz_0_xxyyyyz_0[i] * fbe_0 - g_zzzzz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyyyzz_0[i] = g_zzzzz_0_xxyyyzz_0[i] * fbe_0 - g_zzzzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyyzzz_0[i] = g_zzzzz_0_xxyyzzz_0[i] * fbe_0 - g_zzzzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxyzzzz_0[i] = g_zzzzz_0_xxyzzzz_0[i] * fbe_0 - g_zzzzz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xxzzzzz_0[i] = g_zzzzz_0_xxzzzzz_0[i] * fbe_0 - g_zzzzz_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyyyyy_0[i] = 4.0 * g_xxzzz_0_xyyyyyy_0[i] * fbe_0 - 4.0 * g_xxzzz_0_xyyyyyy_1[i] * fz_be_0 + g_xxzzzz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxzzzzz_0_xyyyyyz_0[i] = g_zzzzz_0_xyyyyyz_0[i] * fbe_0 - g_zzzzz_0_xyyyyyz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyyyzz_0[i] = g_zzzzz_0_xyyyyzz_0[i] * fbe_0 - g_zzzzz_0_xyyyyzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyyzzz_0[i] = g_zzzzz_0_xyyyzzz_0[i] * fbe_0 - g_zzzzz_0_xyyyzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyyzzzz_0[i] = g_zzzzz_0_xyyzzzz_0[i] * fbe_0 - g_zzzzz_0_xyyzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xyzzzzz_0[i] = g_zzzzz_0_xyzzzzz_0[i] * fbe_0 - g_zzzzz_0_xyzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_xzzzzzz_0[i] = g_zzzzz_0_xzzzzzz_0[i] * fbe_0 - g_zzzzz_0_xzzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_xzzzzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyyyy_0[i] = g_zzzzz_0_yyyyyyy_0[i] * fbe_0 - g_zzzzz_0_yyyyyyy_1[i] * fz_be_0 + g_xzzzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyyyz_0[i] = g_zzzzz_0_yyyyyyz_0[i] * fbe_0 - g_zzzzz_0_yyyyyyz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyyzz_0[i] = g_zzzzz_0_yyyyyzz_0[i] * fbe_0 - g_zzzzz_0_yyyyyzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyyzzz_0[i] = g_zzzzz_0_yyyyzzz_0[i] * fbe_0 - g_zzzzz_0_yyyyzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyyzzzz_0[i] = g_zzzzz_0_yyyzzzz_0[i] * fbe_0 - g_zzzzz_0_yyyzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yyzzzzz_0[i] = g_zzzzz_0_yyzzzzz_0[i] * fbe_0 - g_zzzzz_0_yyzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_yzzzzzz_0[i] = g_zzzzz_0_yzzzzzz_0[i] * fbe_0 - g_zzzzz_0_yzzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxzzzzz_0_zzzzzzz_0[i] = g_zzzzz_0_zzzzzzz_0[i] * fbe_0 - g_zzzzz_0_zzzzzzz_1[i] * fz_be_0 + g_xzzzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 756-792 components of targeted buffer : KSK

    auto g_xyyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 756);

    auto g_xyyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 757);

    auto g_xyyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 758);

    auto g_xyyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 759);

    auto g_xyyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 760);

    auto g_xyyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 761);

    auto g_xyyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 762);

    auto g_xyyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 763);

    auto g_xyyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 764);

    auto g_xyyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 765);

    auto g_xyyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 766);

    auto g_xyyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 767);

    auto g_xyyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 768);

    auto g_xyyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 769);

    auto g_xyyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 770);

    auto g_xyyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 771);

    auto g_xyyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 772);

    auto g_xyyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 773);

    auto g_xyyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 774);

    auto g_xyyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 775);

    auto g_xyyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 776);

    auto g_xyyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 777);

    auto g_xyyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 778);

    auto g_xyyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 779);

    auto g_xyyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 780);

    auto g_xyyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 781);

    auto g_xyyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 782);

    auto g_xyyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 783);

    auto g_xyyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 784);

    auto g_xyyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 785);

    auto g_xyyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 786);

    auto g_xyyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 787);

    auto g_xyyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 788);

    auto g_xyyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 789);

    auto g_xyyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 790);

    auto g_xyyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 791);

    #pragma omp simd aligned(g_xyyyyyy_0_xxxxxxx_0, g_xyyyyyy_0_xxxxxxy_0, g_xyyyyyy_0_xxxxxxz_0, g_xyyyyyy_0_xxxxxyy_0, g_xyyyyyy_0_xxxxxyz_0, g_xyyyyyy_0_xxxxxzz_0, g_xyyyyyy_0_xxxxyyy_0, g_xyyyyyy_0_xxxxyyz_0, g_xyyyyyy_0_xxxxyzz_0, g_xyyyyyy_0_xxxxzzz_0, g_xyyyyyy_0_xxxyyyy_0, g_xyyyyyy_0_xxxyyyz_0, g_xyyyyyy_0_xxxyyzz_0, g_xyyyyyy_0_xxxyzzz_0, g_xyyyyyy_0_xxxzzzz_0, g_xyyyyyy_0_xxyyyyy_0, g_xyyyyyy_0_xxyyyyz_0, g_xyyyyyy_0_xxyyyzz_0, g_xyyyyyy_0_xxyyzzz_0, g_xyyyyyy_0_xxyzzzz_0, g_xyyyyyy_0_xxzzzzz_0, g_xyyyyyy_0_xyyyyyy_0, g_xyyyyyy_0_xyyyyyz_0, g_xyyyyyy_0_xyyyyzz_0, g_xyyyyyy_0_xyyyzzz_0, g_xyyyyyy_0_xyyzzzz_0, g_xyyyyyy_0_xyzzzzz_0, g_xyyyyyy_0_xzzzzzz_0, g_xyyyyyy_0_yyyyyyy_0, g_xyyyyyy_0_yyyyyyz_0, g_xyyyyyy_0_yyyyyzz_0, g_xyyyyyy_0_yyyyzzz_0, g_xyyyyyy_0_yyyzzzz_0, g_xyyyyyy_0_yyzzzzz_0, g_xyyyyyy_0_yzzzzzz_0, g_xyyyyyy_0_zzzzzzz_0, g_yyyyyy_0_xxxxxx_1, g_yyyyyy_0_xxxxxxx_1, g_yyyyyy_0_xxxxxxy_1, g_yyyyyy_0_xxxxxxz_1, g_yyyyyy_0_xxxxxy_1, g_yyyyyy_0_xxxxxyy_1, g_yyyyyy_0_xxxxxyz_1, g_yyyyyy_0_xxxxxz_1, g_yyyyyy_0_xxxxxzz_1, g_yyyyyy_0_xxxxyy_1, g_yyyyyy_0_xxxxyyy_1, g_yyyyyy_0_xxxxyyz_1, g_yyyyyy_0_xxxxyz_1, g_yyyyyy_0_xxxxyzz_1, g_yyyyyy_0_xxxxzz_1, g_yyyyyy_0_xxxxzzz_1, g_yyyyyy_0_xxxyyy_1, g_yyyyyy_0_xxxyyyy_1, g_yyyyyy_0_xxxyyyz_1, g_yyyyyy_0_xxxyyz_1, g_yyyyyy_0_xxxyyzz_1, g_yyyyyy_0_xxxyzz_1, g_yyyyyy_0_xxxyzzz_1, g_yyyyyy_0_xxxzzz_1, g_yyyyyy_0_xxxzzzz_1, g_yyyyyy_0_xxyyyy_1, g_yyyyyy_0_xxyyyyy_1, g_yyyyyy_0_xxyyyyz_1, g_yyyyyy_0_xxyyyz_1, g_yyyyyy_0_xxyyyzz_1, g_yyyyyy_0_xxyyzz_1, g_yyyyyy_0_xxyyzzz_1, g_yyyyyy_0_xxyzzz_1, g_yyyyyy_0_xxyzzzz_1, g_yyyyyy_0_xxzzzz_1, g_yyyyyy_0_xxzzzzz_1, g_yyyyyy_0_xyyyyy_1, g_yyyyyy_0_xyyyyyy_1, g_yyyyyy_0_xyyyyyz_1, g_yyyyyy_0_xyyyyz_1, g_yyyyyy_0_xyyyyzz_1, g_yyyyyy_0_xyyyzz_1, g_yyyyyy_0_xyyyzzz_1, g_yyyyyy_0_xyyzzz_1, g_yyyyyy_0_xyyzzzz_1, g_yyyyyy_0_xyzzzz_1, g_yyyyyy_0_xyzzzzz_1, g_yyyyyy_0_xzzzzz_1, g_yyyyyy_0_xzzzzzz_1, g_yyyyyy_0_yyyyyy_1, g_yyyyyy_0_yyyyyyy_1, g_yyyyyy_0_yyyyyyz_1, g_yyyyyy_0_yyyyyz_1, g_yyyyyy_0_yyyyyzz_1, g_yyyyyy_0_yyyyzz_1, g_yyyyyy_0_yyyyzzz_1, g_yyyyyy_0_yyyzzz_1, g_yyyyyy_0_yyyzzzz_1, g_yyyyyy_0_yyzzzz_1, g_yyyyyy_0_yyzzzzz_1, g_yyyyyy_0_yzzzzz_1, g_yyyyyy_0_yzzzzzz_1, g_yyyyyy_0_zzzzzz_1, g_yyyyyy_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyy_0_xxxxxxx_0[i] = 7.0 * g_yyyyyy_0_xxxxxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxx_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxxy_0[i] = 6.0 * g_yyyyyy_0_xxxxxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxxz_0[i] = 6.0 * g_yyyyyy_0_xxxxxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxyy_0[i] = 5.0 * g_yyyyyy_0_xxxxyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxyz_0[i] = 5.0 * g_yyyyyy_0_xxxxyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxxzz_0[i] = 5.0 * g_yyyyyy_0_xxxxzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxyyy_0[i] = 4.0 * g_yyyyyy_0_xxxyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxyyz_0[i] = 4.0 * g_yyyyyy_0_xxxyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxyzz_0[i] = 4.0 * g_yyyyyy_0_xxxyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxxzzz_0[i] = 4.0 * g_yyyyyy_0_xxxzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyyyy_0[i] = 3.0 * g_yyyyyy_0_xxyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyyyz_0[i] = 3.0 * g_yyyyyy_0_xxyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyyzz_0[i] = 3.0 * g_yyyyyy_0_xxyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxyzzz_0[i] = 3.0 * g_yyyyyy_0_xxyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxxzzzz_0[i] = 3.0 * g_yyyyyy_0_xxzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyyyy_0[i] = 2.0 * g_yyyyyy_0_xyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyyyz_0[i] = 2.0 * g_yyyyyy_0_xyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyyzz_0[i] = 2.0 * g_yyyyyy_0_xyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyyzzz_0[i] = 2.0 * g_yyyyyy_0_xyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxyzzzz_0[i] = 2.0 * g_yyyyyy_0_xyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xxzzzzz_0[i] = 2.0 * g_yyyyyy_0_xzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyyyy_0[i] = g_yyyyyy_0_yyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyyyz_0[i] = g_yyyyyy_0_yyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyyzz_0[i] = g_yyyyyy_0_yyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyyzzz_0[i] = g_yyyyyy_0_yyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyyzzzz_0[i] = g_yyyyyy_0_yyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xyzzzzz_0[i] = g_yyyyyy_0_yzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_xzzzzzz_0[i] = g_yyyyyy_0_zzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyyyy_0[i] = g_yyyyyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyyyz_0[i] = g_yyyyyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyyzz_0[i] = g_yyyyyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyyzzz_0[i] = g_yyyyyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyyzzzz_0[i] = g_yyyyyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yyzzzzz_0[i] = g_yyyyyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_yzzzzzz_0[i] = g_yyyyyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyyyyy_0_zzzzzzz_0[i] = g_yyyyyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 792-828 components of targeted buffer : KSK

    auto g_xyyyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 792);

    auto g_xyyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 793);

    auto g_xyyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 794);

    auto g_xyyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 795);

    auto g_xyyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 796);

    auto g_xyyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 797);

    auto g_xyyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 798);

    auto g_xyyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 799);

    auto g_xyyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 800);

    auto g_xyyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 801);

    auto g_xyyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 802);

    auto g_xyyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 803);

    auto g_xyyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 804);

    auto g_xyyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 805);

    auto g_xyyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 806);

    auto g_xyyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 807);

    auto g_xyyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 808);

    auto g_xyyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 809);

    auto g_xyyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 810);

    auto g_xyyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 811);

    auto g_xyyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 812);

    auto g_xyyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 813);

    auto g_xyyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 814);

    auto g_xyyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 815);

    auto g_xyyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 816);

    auto g_xyyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 817);

    auto g_xyyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 818);

    auto g_xyyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 819);

    auto g_xyyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 820);

    auto g_xyyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 821);

    auto g_xyyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 822);

    auto g_xyyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 823);

    auto g_xyyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 824);

    auto g_xyyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 825);

    auto g_xyyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 826);

    auto g_xyyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 827);

    #pragma omp simd aligned(g_xyyyyy_0_xxxxxxx_1, g_xyyyyy_0_xxxxxxy_1, g_xyyyyy_0_xxxxxyy_1, g_xyyyyy_0_xxxxyyy_1, g_xyyyyy_0_xxxyyyy_1, g_xyyyyy_0_xxyyyyy_1, g_xyyyyy_0_xyyyyyy_1, g_xyyyyyz_0_xxxxxxx_0, g_xyyyyyz_0_xxxxxxy_0, g_xyyyyyz_0_xxxxxxz_0, g_xyyyyyz_0_xxxxxyy_0, g_xyyyyyz_0_xxxxxyz_0, g_xyyyyyz_0_xxxxxzz_0, g_xyyyyyz_0_xxxxyyy_0, g_xyyyyyz_0_xxxxyyz_0, g_xyyyyyz_0_xxxxyzz_0, g_xyyyyyz_0_xxxxzzz_0, g_xyyyyyz_0_xxxyyyy_0, g_xyyyyyz_0_xxxyyyz_0, g_xyyyyyz_0_xxxyyzz_0, g_xyyyyyz_0_xxxyzzz_0, g_xyyyyyz_0_xxxzzzz_0, g_xyyyyyz_0_xxyyyyy_0, g_xyyyyyz_0_xxyyyyz_0, g_xyyyyyz_0_xxyyyzz_0, g_xyyyyyz_0_xxyyzzz_0, g_xyyyyyz_0_xxyzzzz_0, g_xyyyyyz_0_xxzzzzz_0, g_xyyyyyz_0_xyyyyyy_0, g_xyyyyyz_0_xyyyyyz_0, g_xyyyyyz_0_xyyyyzz_0, g_xyyyyyz_0_xyyyzzz_0, g_xyyyyyz_0_xyyzzzz_0, g_xyyyyyz_0_xyzzzzz_0, g_xyyyyyz_0_xzzzzzz_0, g_xyyyyyz_0_yyyyyyy_0, g_xyyyyyz_0_yyyyyyz_0, g_xyyyyyz_0_yyyyyzz_0, g_xyyyyyz_0_yyyyzzz_0, g_xyyyyyz_0_yyyzzzz_0, g_xyyyyyz_0_yyzzzzz_0, g_xyyyyyz_0_yzzzzzz_0, g_xyyyyyz_0_zzzzzzz_0, g_yyyyyz_0_xxxxxxz_1, g_yyyyyz_0_xxxxxyz_1, g_yyyyyz_0_xxxxxz_1, g_yyyyyz_0_xxxxxzz_1, g_yyyyyz_0_xxxxyyz_1, g_yyyyyz_0_xxxxyz_1, g_yyyyyz_0_xxxxyzz_1, g_yyyyyz_0_xxxxzz_1, g_yyyyyz_0_xxxxzzz_1, g_yyyyyz_0_xxxyyyz_1, g_yyyyyz_0_xxxyyz_1, g_yyyyyz_0_xxxyyzz_1, g_yyyyyz_0_xxxyzz_1, g_yyyyyz_0_xxxyzzz_1, g_yyyyyz_0_xxxzzz_1, g_yyyyyz_0_xxxzzzz_1, g_yyyyyz_0_xxyyyyz_1, g_yyyyyz_0_xxyyyz_1, g_yyyyyz_0_xxyyyzz_1, g_yyyyyz_0_xxyyzz_1, g_yyyyyz_0_xxyyzzz_1, g_yyyyyz_0_xxyzzz_1, g_yyyyyz_0_xxyzzzz_1, g_yyyyyz_0_xxzzzz_1, g_yyyyyz_0_xxzzzzz_1, g_yyyyyz_0_xyyyyyz_1, g_yyyyyz_0_xyyyyz_1, g_yyyyyz_0_xyyyyzz_1, g_yyyyyz_0_xyyyzz_1, g_yyyyyz_0_xyyyzzz_1, g_yyyyyz_0_xyyzzz_1, g_yyyyyz_0_xyyzzzz_1, g_yyyyyz_0_xyzzzz_1, g_yyyyyz_0_xyzzzzz_1, g_yyyyyz_0_xzzzzz_1, g_yyyyyz_0_xzzzzzz_1, g_yyyyyz_0_yyyyyyy_1, g_yyyyyz_0_yyyyyyz_1, g_yyyyyz_0_yyyyyz_1, g_yyyyyz_0_yyyyyzz_1, g_yyyyyz_0_yyyyzz_1, g_yyyyyz_0_yyyyzzz_1, g_yyyyyz_0_yyyzzz_1, g_yyyyyz_0_yyyzzzz_1, g_yyyyyz_0_yyzzzz_1, g_yyyyyz_0_yyzzzzz_1, g_yyyyyz_0_yzzzzz_1, g_yyyyyz_0_yzzzzzz_1, g_yyyyyz_0_zzzzzz_1, g_yyyyyz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyyz_0_xxxxxxx_0[i] = g_xyyyyy_0_xxxxxxx_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxxxy_0[i] = g_xyyyyy_0_xxxxxxy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxxxz_0[i] = 6.0 * g_yyyyyz_0_xxxxxz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxxyy_0[i] = g_xyyyyy_0_xxxxxyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxxyz_0[i] = 5.0 * g_yyyyyz_0_xxxxyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxxzz_0[i] = 5.0 * g_yyyyyz_0_xxxxzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxyyy_0[i] = g_xyyyyy_0_xxxxyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxxyyz_0[i] = 4.0 * g_yyyyyz_0_xxxyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxyzz_0[i] = 4.0 * g_yyyyyz_0_xxxyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxxzzz_0[i] = 4.0 * g_yyyyyz_0_xxxzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxyyyy_0[i] = g_xyyyyy_0_xxxyyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxxyyyz_0[i] = 3.0 * g_yyyyyz_0_xxyyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxyyzz_0[i] = 3.0 * g_yyyyyz_0_xxyyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxyzzz_0[i] = 3.0 * g_yyyyyz_0_xxyzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxxzzzz_0[i] = 3.0 * g_yyyyyz_0_xxzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyyyyy_0[i] = g_xyyyyy_0_xxyyyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xxyyyyz_0[i] = 2.0 * g_yyyyyz_0_xyyyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyyyzz_0[i] = 2.0 * g_yyyyyz_0_xyyyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyyzzz_0[i] = 2.0 * g_yyyyyz_0_xyyzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxyzzzz_0[i] = 2.0 * g_yyyyyz_0_xyzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xxzzzzz_0[i] = 2.0 * g_yyyyyz_0_xzzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyyyyy_0[i] = g_xyyyyy_0_xyyyyyy_1[i] * wa_z[i];

        g_xyyyyyz_0_xyyyyyz_0[i] = g_yyyyyz_0_yyyyyz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyyyzz_0[i] = g_yyyyyz_0_yyyyzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyyzzz_0[i] = g_yyyyyz_0_yyyzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyyzzzz_0[i] = g_yyyyyz_0_yyzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xyzzzzz_0[i] = g_yyyyyz_0_yzzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_xzzzzzz_0[i] = g_yyyyyz_0_zzzzzz_1[i] * fi_acd_0 + g_yyyyyz_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyyyy_0[i] = g_yyyyyz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyyyz_0[i] = g_yyyyyz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyyzz_0[i] = g_yyyyyz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyyzzz_0[i] = g_yyyyyz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyyzzzz_0[i] = g_yyyyyz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yyzzzzz_0[i] = g_yyyyyz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_yzzzzzz_0[i] = g_yyyyyz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyyyyz_0_zzzzzzz_0[i] = g_yyyyyz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 828-864 components of targeted buffer : KSK

    auto g_xyyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 828);

    auto g_xyyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 829);

    auto g_xyyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 830);

    auto g_xyyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 831);

    auto g_xyyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 832);

    auto g_xyyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 833);

    auto g_xyyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 834);

    auto g_xyyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 835);

    auto g_xyyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 836);

    auto g_xyyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 837);

    auto g_xyyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 838);

    auto g_xyyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 839);

    auto g_xyyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 840);

    auto g_xyyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 841);

    auto g_xyyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 842);

    auto g_xyyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 843);

    auto g_xyyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 844);

    auto g_xyyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 845);

    auto g_xyyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 846);

    auto g_xyyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 847);

    auto g_xyyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 848);

    auto g_xyyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 849);

    auto g_xyyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 850);

    auto g_xyyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 851);

    auto g_xyyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 852);

    auto g_xyyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 853);

    auto g_xyyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 854);

    auto g_xyyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 855);

    auto g_xyyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 856);

    auto g_xyyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 857);

    auto g_xyyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 858);

    auto g_xyyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 859);

    auto g_xyyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 860);

    auto g_xyyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 861);

    auto g_xyyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 862);

    auto g_xyyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 863);

    #pragma omp simd aligned(g_xyyyyzz_0_xxxxxxx_0, g_xyyyyzz_0_xxxxxxy_0, g_xyyyyzz_0_xxxxxxz_0, g_xyyyyzz_0_xxxxxyy_0, g_xyyyyzz_0_xxxxxyz_0, g_xyyyyzz_0_xxxxxzz_0, g_xyyyyzz_0_xxxxyyy_0, g_xyyyyzz_0_xxxxyyz_0, g_xyyyyzz_0_xxxxyzz_0, g_xyyyyzz_0_xxxxzzz_0, g_xyyyyzz_0_xxxyyyy_0, g_xyyyyzz_0_xxxyyyz_0, g_xyyyyzz_0_xxxyyzz_0, g_xyyyyzz_0_xxxyzzz_0, g_xyyyyzz_0_xxxzzzz_0, g_xyyyyzz_0_xxyyyyy_0, g_xyyyyzz_0_xxyyyyz_0, g_xyyyyzz_0_xxyyyzz_0, g_xyyyyzz_0_xxyyzzz_0, g_xyyyyzz_0_xxyzzzz_0, g_xyyyyzz_0_xxzzzzz_0, g_xyyyyzz_0_xyyyyyy_0, g_xyyyyzz_0_xyyyyyz_0, g_xyyyyzz_0_xyyyyzz_0, g_xyyyyzz_0_xyyyzzz_0, g_xyyyyzz_0_xyyzzzz_0, g_xyyyyzz_0_xyzzzzz_0, g_xyyyyzz_0_xzzzzzz_0, g_xyyyyzz_0_yyyyyyy_0, g_xyyyyzz_0_yyyyyyz_0, g_xyyyyzz_0_yyyyyzz_0, g_xyyyyzz_0_yyyyzzz_0, g_xyyyyzz_0_yyyzzzz_0, g_xyyyyzz_0_yyzzzzz_0, g_xyyyyzz_0_yzzzzzz_0, g_xyyyyzz_0_zzzzzzz_0, g_yyyyzz_0_xxxxxx_1, g_yyyyzz_0_xxxxxxx_1, g_yyyyzz_0_xxxxxxy_1, g_yyyyzz_0_xxxxxxz_1, g_yyyyzz_0_xxxxxy_1, g_yyyyzz_0_xxxxxyy_1, g_yyyyzz_0_xxxxxyz_1, g_yyyyzz_0_xxxxxz_1, g_yyyyzz_0_xxxxxzz_1, g_yyyyzz_0_xxxxyy_1, g_yyyyzz_0_xxxxyyy_1, g_yyyyzz_0_xxxxyyz_1, g_yyyyzz_0_xxxxyz_1, g_yyyyzz_0_xxxxyzz_1, g_yyyyzz_0_xxxxzz_1, g_yyyyzz_0_xxxxzzz_1, g_yyyyzz_0_xxxyyy_1, g_yyyyzz_0_xxxyyyy_1, g_yyyyzz_0_xxxyyyz_1, g_yyyyzz_0_xxxyyz_1, g_yyyyzz_0_xxxyyzz_1, g_yyyyzz_0_xxxyzz_1, g_yyyyzz_0_xxxyzzz_1, g_yyyyzz_0_xxxzzz_1, g_yyyyzz_0_xxxzzzz_1, g_yyyyzz_0_xxyyyy_1, g_yyyyzz_0_xxyyyyy_1, g_yyyyzz_0_xxyyyyz_1, g_yyyyzz_0_xxyyyz_1, g_yyyyzz_0_xxyyyzz_1, g_yyyyzz_0_xxyyzz_1, g_yyyyzz_0_xxyyzzz_1, g_yyyyzz_0_xxyzzz_1, g_yyyyzz_0_xxyzzzz_1, g_yyyyzz_0_xxzzzz_1, g_yyyyzz_0_xxzzzzz_1, g_yyyyzz_0_xyyyyy_1, g_yyyyzz_0_xyyyyyy_1, g_yyyyzz_0_xyyyyyz_1, g_yyyyzz_0_xyyyyz_1, g_yyyyzz_0_xyyyyzz_1, g_yyyyzz_0_xyyyzz_1, g_yyyyzz_0_xyyyzzz_1, g_yyyyzz_0_xyyzzz_1, g_yyyyzz_0_xyyzzzz_1, g_yyyyzz_0_xyzzzz_1, g_yyyyzz_0_xyzzzzz_1, g_yyyyzz_0_xzzzzz_1, g_yyyyzz_0_xzzzzzz_1, g_yyyyzz_0_yyyyyy_1, g_yyyyzz_0_yyyyyyy_1, g_yyyyzz_0_yyyyyyz_1, g_yyyyzz_0_yyyyyz_1, g_yyyyzz_0_yyyyyzz_1, g_yyyyzz_0_yyyyzz_1, g_yyyyzz_0_yyyyzzz_1, g_yyyyzz_0_yyyzzz_1, g_yyyyzz_0_yyyzzzz_1, g_yyyyzz_0_yyzzzz_1, g_yyyyzz_0_yyzzzzz_1, g_yyyyzz_0_yzzzzz_1, g_yyyyzz_0_yzzzzzz_1, g_yyyyzz_0_zzzzzz_1, g_yyyyzz_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyzz_0_xxxxxxx_0[i] = 7.0 * g_yyyyzz_0_xxxxxx_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxxx_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxxy_0[i] = 6.0 * g_yyyyzz_0_xxxxxy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxxz_0[i] = 6.0 * g_yyyyzz_0_xxxxxz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxyy_0[i] = 5.0 * g_yyyyzz_0_xxxxyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxyz_0[i] = 5.0 * g_yyyyzz_0_xxxxyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxxzz_0[i] = 5.0 * g_yyyyzz_0_xxxxzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxyyy_0[i] = 4.0 * g_yyyyzz_0_xxxyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxyyz_0[i] = 4.0 * g_yyyyzz_0_xxxyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxyzz_0[i] = 4.0 * g_yyyyzz_0_xxxyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxxzzz_0[i] = 4.0 * g_yyyyzz_0_xxxzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyyyy_0[i] = 3.0 * g_yyyyzz_0_xxyyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyyyz_0[i] = 3.0 * g_yyyyzz_0_xxyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyyzz_0[i] = 3.0 * g_yyyyzz_0_xxyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxyzzz_0[i] = 3.0 * g_yyyyzz_0_xxyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxxzzzz_0[i] = 3.0 * g_yyyyzz_0_xxzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyyyy_0[i] = 2.0 * g_yyyyzz_0_xyyyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyyyz_0[i] = 2.0 * g_yyyyzz_0_xyyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyyzz_0[i] = 2.0 * g_yyyyzz_0_xyyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyyzzz_0[i] = 2.0 * g_yyyyzz_0_xyyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxyzzzz_0[i] = 2.0 * g_yyyyzz_0_xyzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xxzzzzz_0[i] = 2.0 * g_yyyyzz_0_xzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyyyy_0[i] = g_yyyyzz_0_yyyyyy_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyyyz_0[i] = g_yyyyzz_0_yyyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyyzz_0[i] = g_yyyyzz_0_yyyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyyzzz_0[i] = g_yyyyzz_0_yyyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyyzzzz_0[i] = g_yyyyzz_0_yyzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xyzzzzz_0[i] = g_yyyyzz_0_yzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_xzzzzzz_0[i] = g_yyyyzz_0_zzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyyyy_0[i] = g_yyyyzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyyyz_0[i] = g_yyyyzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyyzz_0[i] = g_yyyyzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyyzzz_0[i] = g_yyyyzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyyzzzz_0[i] = g_yyyyzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yyzzzzz_0[i] = g_yyyyzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_yzzzzzz_0[i] = g_yyyyzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyyyzz_0_zzzzzzz_0[i] = g_yyyyzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 864-900 components of targeted buffer : KSK

    auto g_xyyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 864);

    auto g_xyyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 865);

    auto g_xyyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 866);

    auto g_xyyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 867);

    auto g_xyyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 868);

    auto g_xyyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 869);

    auto g_xyyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 870);

    auto g_xyyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 871);

    auto g_xyyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 872);

    auto g_xyyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 873);

    auto g_xyyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 874);

    auto g_xyyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 875);

    auto g_xyyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 876);

    auto g_xyyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 877);

    auto g_xyyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 878);

    auto g_xyyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 879);

    auto g_xyyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 880);

    auto g_xyyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 881);

    auto g_xyyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 882);

    auto g_xyyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 883);

    auto g_xyyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 884);

    auto g_xyyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 885);

    auto g_xyyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 886);

    auto g_xyyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 887);

    auto g_xyyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 888);

    auto g_xyyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 889);

    auto g_xyyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 890);

    auto g_xyyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 891);

    auto g_xyyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 892);

    auto g_xyyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 893);

    auto g_xyyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 894);

    auto g_xyyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 895);

    auto g_xyyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 896);

    auto g_xyyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 897);

    auto g_xyyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 898);

    auto g_xyyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 899);

    #pragma omp simd aligned(g_xyyyzzz_0_xxxxxxx_0, g_xyyyzzz_0_xxxxxxy_0, g_xyyyzzz_0_xxxxxxz_0, g_xyyyzzz_0_xxxxxyy_0, g_xyyyzzz_0_xxxxxyz_0, g_xyyyzzz_0_xxxxxzz_0, g_xyyyzzz_0_xxxxyyy_0, g_xyyyzzz_0_xxxxyyz_0, g_xyyyzzz_0_xxxxyzz_0, g_xyyyzzz_0_xxxxzzz_0, g_xyyyzzz_0_xxxyyyy_0, g_xyyyzzz_0_xxxyyyz_0, g_xyyyzzz_0_xxxyyzz_0, g_xyyyzzz_0_xxxyzzz_0, g_xyyyzzz_0_xxxzzzz_0, g_xyyyzzz_0_xxyyyyy_0, g_xyyyzzz_0_xxyyyyz_0, g_xyyyzzz_0_xxyyyzz_0, g_xyyyzzz_0_xxyyzzz_0, g_xyyyzzz_0_xxyzzzz_0, g_xyyyzzz_0_xxzzzzz_0, g_xyyyzzz_0_xyyyyyy_0, g_xyyyzzz_0_xyyyyyz_0, g_xyyyzzz_0_xyyyyzz_0, g_xyyyzzz_0_xyyyzzz_0, g_xyyyzzz_0_xyyzzzz_0, g_xyyyzzz_0_xyzzzzz_0, g_xyyyzzz_0_xzzzzzz_0, g_xyyyzzz_0_yyyyyyy_0, g_xyyyzzz_0_yyyyyyz_0, g_xyyyzzz_0_yyyyyzz_0, g_xyyyzzz_0_yyyyzzz_0, g_xyyyzzz_0_yyyzzzz_0, g_xyyyzzz_0_yyzzzzz_0, g_xyyyzzz_0_yzzzzzz_0, g_xyyyzzz_0_zzzzzzz_0, g_yyyzzz_0_xxxxxx_1, g_yyyzzz_0_xxxxxxx_1, g_yyyzzz_0_xxxxxxy_1, g_yyyzzz_0_xxxxxxz_1, g_yyyzzz_0_xxxxxy_1, g_yyyzzz_0_xxxxxyy_1, g_yyyzzz_0_xxxxxyz_1, g_yyyzzz_0_xxxxxz_1, g_yyyzzz_0_xxxxxzz_1, g_yyyzzz_0_xxxxyy_1, g_yyyzzz_0_xxxxyyy_1, g_yyyzzz_0_xxxxyyz_1, g_yyyzzz_0_xxxxyz_1, g_yyyzzz_0_xxxxyzz_1, g_yyyzzz_0_xxxxzz_1, g_yyyzzz_0_xxxxzzz_1, g_yyyzzz_0_xxxyyy_1, g_yyyzzz_0_xxxyyyy_1, g_yyyzzz_0_xxxyyyz_1, g_yyyzzz_0_xxxyyz_1, g_yyyzzz_0_xxxyyzz_1, g_yyyzzz_0_xxxyzz_1, g_yyyzzz_0_xxxyzzz_1, g_yyyzzz_0_xxxzzz_1, g_yyyzzz_0_xxxzzzz_1, g_yyyzzz_0_xxyyyy_1, g_yyyzzz_0_xxyyyyy_1, g_yyyzzz_0_xxyyyyz_1, g_yyyzzz_0_xxyyyz_1, g_yyyzzz_0_xxyyyzz_1, g_yyyzzz_0_xxyyzz_1, g_yyyzzz_0_xxyyzzz_1, g_yyyzzz_0_xxyzzz_1, g_yyyzzz_0_xxyzzzz_1, g_yyyzzz_0_xxzzzz_1, g_yyyzzz_0_xxzzzzz_1, g_yyyzzz_0_xyyyyy_1, g_yyyzzz_0_xyyyyyy_1, g_yyyzzz_0_xyyyyyz_1, g_yyyzzz_0_xyyyyz_1, g_yyyzzz_0_xyyyyzz_1, g_yyyzzz_0_xyyyzz_1, g_yyyzzz_0_xyyyzzz_1, g_yyyzzz_0_xyyzzz_1, g_yyyzzz_0_xyyzzzz_1, g_yyyzzz_0_xyzzzz_1, g_yyyzzz_0_xyzzzzz_1, g_yyyzzz_0_xzzzzz_1, g_yyyzzz_0_xzzzzzz_1, g_yyyzzz_0_yyyyyy_1, g_yyyzzz_0_yyyyyyy_1, g_yyyzzz_0_yyyyyyz_1, g_yyyzzz_0_yyyyyz_1, g_yyyzzz_0_yyyyyzz_1, g_yyyzzz_0_yyyyzz_1, g_yyyzzz_0_yyyyzzz_1, g_yyyzzz_0_yyyzzz_1, g_yyyzzz_0_yyyzzzz_1, g_yyyzzz_0_yyzzzz_1, g_yyyzzz_0_yyzzzzz_1, g_yyyzzz_0_yzzzzz_1, g_yyyzzz_0_yzzzzzz_1, g_yyyzzz_0_zzzzzz_1, g_yyyzzz_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzzz_0_xxxxxxx_0[i] = 7.0 * g_yyyzzz_0_xxxxxx_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxxx_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxxy_0[i] = 6.0 * g_yyyzzz_0_xxxxxy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxxz_0[i] = 6.0 * g_yyyzzz_0_xxxxxz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxyy_0[i] = 5.0 * g_yyyzzz_0_xxxxyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxyz_0[i] = 5.0 * g_yyyzzz_0_xxxxyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxxzz_0[i] = 5.0 * g_yyyzzz_0_xxxxzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxyyy_0[i] = 4.0 * g_yyyzzz_0_xxxyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxyyz_0[i] = 4.0 * g_yyyzzz_0_xxxyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxyzz_0[i] = 4.0 * g_yyyzzz_0_xxxyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxxzzz_0[i] = 4.0 * g_yyyzzz_0_xxxzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyyyy_0[i] = 3.0 * g_yyyzzz_0_xxyyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyyyz_0[i] = 3.0 * g_yyyzzz_0_xxyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyyzz_0[i] = 3.0 * g_yyyzzz_0_xxyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxyzzz_0[i] = 3.0 * g_yyyzzz_0_xxyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxxzzzz_0[i] = 3.0 * g_yyyzzz_0_xxzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyyyy_0[i] = 2.0 * g_yyyzzz_0_xyyyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyyyz_0[i] = 2.0 * g_yyyzzz_0_xyyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyyzz_0[i] = 2.0 * g_yyyzzz_0_xyyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyyzzz_0[i] = 2.0 * g_yyyzzz_0_xyyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxyzzzz_0[i] = 2.0 * g_yyyzzz_0_xyzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xxzzzzz_0[i] = 2.0 * g_yyyzzz_0_xzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyyyy_0[i] = g_yyyzzz_0_yyyyyy_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyyyz_0[i] = g_yyyzzz_0_yyyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyyzz_0[i] = g_yyyzzz_0_yyyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyyzzz_0[i] = g_yyyzzz_0_yyyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyyzzzz_0[i] = g_yyyzzz_0_yyzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xyzzzzz_0[i] = g_yyyzzz_0_yzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_xzzzzzz_0[i] = g_yyyzzz_0_zzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyyyy_0[i] = g_yyyzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyyyz_0[i] = g_yyyzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyyzz_0[i] = g_yyyzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyyzzz_0[i] = g_yyyzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyyzzzz_0[i] = g_yyyzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yyzzzzz_0[i] = g_yyyzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_yzzzzzz_0[i] = g_yyyzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyyzzz_0_zzzzzzz_0[i] = g_yyyzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 900-936 components of targeted buffer : KSK

    auto g_xyyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 900);

    auto g_xyyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 901);

    auto g_xyyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 902);

    auto g_xyyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 903);

    auto g_xyyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 904);

    auto g_xyyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 905);

    auto g_xyyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 906);

    auto g_xyyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 907);

    auto g_xyyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 908);

    auto g_xyyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 909);

    auto g_xyyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 910);

    auto g_xyyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 911);

    auto g_xyyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 912);

    auto g_xyyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 913);

    auto g_xyyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 914);

    auto g_xyyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 915);

    auto g_xyyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 916);

    auto g_xyyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 917);

    auto g_xyyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 918);

    auto g_xyyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 919);

    auto g_xyyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 920);

    auto g_xyyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 921);

    auto g_xyyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 922);

    auto g_xyyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 923);

    auto g_xyyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 924);

    auto g_xyyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 925);

    auto g_xyyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 926);

    auto g_xyyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 927);

    auto g_xyyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 928);

    auto g_xyyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 929);

    auto g_xyyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 930);

    auto g_xyyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 931);

    auto g_xyyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 932);

    auto g_xyyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 933);

    auto g_xyyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 934);

    auto g_xyyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 935);

    #pragma omp simd aligned(g_xyyzzzz_0_xxxxxxx_0, g_xyyzzzz_0_xxxxxxy_0, g_xyyzzzz_0_xxxxxxz_0, g_xyyzzzz_0_xxxxxyy_0, g_xyyzzzz_0_xxxxxyz_0, g_xyyzzzz_0_xxxxxzz_0, g_xyyzzzz_0_xxxxyyy_0, g_xyyzzzz_0_xxxxyyz_0, g_xyyzzzz_0_xxxxyzz_0, g_xyyzzzz_0_xxxxzzz_0, g_xyyzzzz_0_xxxyyyy_0, g_xyyzzzz_0_xxxyyyz_0, g_xyyzzzz_0_xxxyyzz_0, g_xyyzzzz_0_xxxyzzz_0, g_xyyzzzz_0_xxxzzzz_0, g_xyyzzzz_0_xxyyyyy_0, g_xyyzzzz_0_xxyyyyz_0, g_xyyzzzz_0_xxyyyzz_0, g_xyyzzzz_0_xxyyzzz_0, g_xyyzzzz_0_xxyzzzz_0, g_xyyzzzz_0_xxzzzzz_0, g_xyyzzzz_0_xyyyyyy_0, g_xyyzzzz_0_xyyyyyz_0, g_xyyzzzz_0_xyyyyzz_0, g_xyyzzzz_0_xyyyzzz_0, g_xyyzzzz_0_xyyzzzz_0, g_xyyzzzz_0_xyzzzzz_0, g_xyyzzzz_0_xzzzzzz_0, g_xyyzzzz_0_yyyyyyy_0, g_xyyzzzz_0_yyyyyyz_0, g_xyyzzzz_0_yyyyyzz_0, g_xyyzzzz_0_yyyyzzz_0, g_xyyzzzz_0_yyyzzzz_0, g_xyyzzzz_0_yyzzzzz_0, g_xyyzzzz_0_yzzzzzz_0, g_xyyzzzz_0_zzzzzzz_0, g_yyzzzz_0_xxxxxx_1, g_yyzzzz_0_xxxxxxx_1, g_yyzzzz_0_xxxxxxy_1, g_yyzzzz_0_xxxxxxz_1, g_yyzzzz_0_xxxxxy_1, g_yyzzzz_0_xxxxxyy_1, g_yyzzzz_0_xxxxxyz_1, g_yyzzzz_0_xxxxxz_1, g_yyzzzz_0_xxxxxzz_1, g_yyzzzz_0_xxxxyy_1, g_yyzzzz_0_xxxxyyy_1, g_yyzzzz_0_xxxxyyz_1, g_yyzzzz_0_xxxxyz_1, g_yyzzzz_0_xxxxyzz_1, g_yyzzzz_0_xxxxzz_1, g_yyzzzz_0_xxxxzzz_1, g_yyzzzz_0_xxxyyy_1, g_yyzzzz_0_xxxyyyy_1, g_yyzzzz_0_xxxyyyz_1, g_yyzzzz_0_xxxyyz_1, g_yyzzzz_0_xxxyyzz_1, g_yyzzzz_0_xxxyzz_1, g_yyzzzz_0_xxxyzzz_1, g_yyzzzz_0_xxxzzz_1, g_yyzzzz_0_xxxzzzz_1, g_yyzzzz_0_xxyyyy_1, g_yyzzzz_0_xxyyyyy_1, g_yyzzzz_0_xxyyyyz_1, g_yyzzzz_0_xxyyyz_1, g_yyzzzz_0_xxyyyzz_1, g_yyzzzz_0_xxyyzz_1, g_yyzzzz_0_xxyyzzz_1, g_yyzzzz_0_xxyzzz_1, g_yyzzzz_0_xxyzzzz_1, g_yyzzzz_0_xxzzzz_1, g_yyzzzz_0_xxzzzzz_1, g_yyzzzz_0_xyyyyy_1, g_yyzzzz_0_xyyyyyy_1, g_yyzzzz_0_xyyyyyz_1, g_yyzzzz_0_xyyyyz_1, g_yyzzzz_0_xyyyyzz_1, g_yyzzzz_0_xyyyzz_1, g_yyzzzz_0_xyyyzzz_1, g_yyzzzz_0_xyyzzz_1, g_yyzzzz_0_xyyzzzz_1, g_yyzzzz_0_xyzzzz_1, g_yyzzzz_0_xyzzzzz_1, g_yyzzzz_0_xzzzzz_1, g_yyzzzz_0_xzzzzzz_1, g_yyzzzz_0_yyyyyy_1, g_yyzzzz_0_yyyyyyy_1, g_yyzzzz_0_yyyyyyz_1, g_yyzzzz_0_yyyyyz_1, g_yyzzzz_0_yyyyyzz_1, g_yyzzzz_0_yyyyzz_1, g_yyzzzz_0_yyyyzzz_1, g_yyzzzz_0_yyyzzz_1, g_yyzzzz_0_yyyzzzz_1, g_yyzzzz_0_yyzzzz_1, g_yyzzzz_0_yyzzzzz_1, g_yyzzzz_0_yzzzzz_1, g_yyzzzz_0_yzzzzzz_1, g_yyzzzz_0_zzzzzz_1, g_yyzzzz_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzzz_0_xxxxxxx_0[i] = 7.0 * g_yyzzzz_0_xxxxxx_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxxx_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxxy_0[i] = 6.0 * g_yyzzzz_0_xxxxxy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxxz_0[i] = 6.0 * g_yyzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxyy_0[i] = 5.0 * g_yyzzzz_0_xxxxyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxyz_0[i] = 5.0 * g_yyzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxxzz_0[i] = 5.0 * g_yyzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxyyy_0[i] = 4.0 * g_yyzzzz_0_xxxyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxyyz_0[i] = 4.0 * g_yyzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxyzz_0[i] = 4.0 * g_yyzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxxzzz_0[i] = 4.0 * g_yyzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyyyy_0[i] = 3.0 * g_yyzzzz_0_xxyyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyyyz_0[i] = 3.0 * g_yyzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyyzz_0[i] = 3.0 * g_yyzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxyzzz_0[i] = 3.0 * g_yyzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxxzzzz_0[i] = 3.0 * g_yyzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyyyy_0[i] = 2.0 * g_yyzzzz_0_xyyyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyyyz_0[i] = 2.0 * g_yyzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyyzz_0[i] = 2.0 * g_yyzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyyzzz_0[i] = 2.0 * g_yyzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxyzzzz_0[i] = 2.0 * g_yyzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xxzzzzz_0[i] = 2.0 * g_yyzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyyyy_0[i] = g_yyzzzz_0_yyyyyy_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyyyz_0[i] = g_yyzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyyzz_0[i] = g_yyzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyyzzz_0[i] = g_yyzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyyzzzz_0[i] = g_yyzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xyzzzzz_0[i] = g_yyzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_xzzzzzz_0[i] = g_yyzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyyyy_0[i] = g_yyzzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyyyz_0[i] = g_yyzzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyyzz_0[i] = g_yyzzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyyzzz_0[i] = g_yyzzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyyzzzz_0[i] = g_yyzzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yyzzzzz_0[i] = g_yyzzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_yzzzzzz_0[i] = g_yyzzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyzzzz_0_zzzzzzz_0[i] = g_yyzzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 936-972 components of targeted buffer : KSK

    auto g_xyzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 936);

    auto g_xyzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 937);

    auto g_xyzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 938);

    auto g_xyzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 939);

    auto g_xyzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 940);

    auto g_xyzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 941);

    auto g_xyzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 942);

    auto g_xyzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 943);

    auto g_xyzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 944);

    auto g_xyzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 945);

    auto g_xyzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 946);

    auto g_xyzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 947);

    auto g_xyzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 948);

    auto g_xyzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 949);

    auto g_xyzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 950);

    auto g_xyzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 951);

    auto g_xyzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 952);

    auto g_xyzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 953);

    auto g_xyzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 954);

    auto g_xyzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 955);

    auto g_xyzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 956);

    auto g_xyzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 957);

    auto g_xyzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 958);

    auto g_xyzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 959);

    auto g_xyzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 960);

    auto g_xyzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 961);

    auto g_xyzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 962);

    auto g_xyzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 963);

    auto g_xyzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 964);

    auto g_xyzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 965);

    auto g_xyzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 966);

    auto g_xyzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 967);

    auto g_xyzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 968);

    auto g_xyzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 969);

    auto g_xyzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 970);

    auto g_xyzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 971);

    #pragma omp simd aligned(g_xyzzzzz_0_xxxxxxx_0, g_xyzzzzz_0_xxxxxxy_0, g_xyzzzzz_0_xxxxxxz_0, g_xyzzzzz_0_xxxxxyy_0, g_xyzzzzz_0_xxxxxyz_0, g_xyzzzzz_0_xxxxxzz_0, g_xyzzzzz_0_xxxxyyy_0, g_xyzzzzz_0_xxxxyyz_0, g_xyzzzzz_0_xxxxyzz_0, g_xyzzzzz_0_xxxxzzz_0, g_xyzzzzz_0_xxxyyyy_0, g_xyzzzzz_0_xxxyyyz_0, g_xyzzzzz_0_xxxyyzz_0, g_xyzzzzz_0_xxxyzzz_0, g_xyzzzzz_0_xxxzzzz_0, g_xyzzzzz_0_xxyyyyy_0, g_xyzzzzz_0_xxyyyyz_0, g_xyzzzzz_0_xxyyyzz_0, g_xyzzzzz_0_xxyyzzz_0, g_xyzzzzz_0_xxyzzzz_0, g_xyzzzzz_0_xxzzzzz_0, g_xyzzzzz_0_xyyyyyy_0, g_xyzzzzz_0_xyyyyyz_0, g_xyzzzzz_0_xyyyyzz_0, g_xyzzzzz_0_xyyyzzz_0, g_xyzzzzz_0_xyyzzzz_0, g_xyzzzzz_0_xyzzzzz_0, g_xyzzzzz_0_xzzzzzz_0, g_xyzzzzz_0_yyyyyyy_0, g_xyzzzzz_0_yyyyyyz_0, g_xyzzzzz_0_yyyyyzz_0, g_xyzzzzz_0_yyyyzzz_0, g_xyzzzzz_0_yyyzzzz_0, g_xyzzzzz_0_yyzzzzz_0, g_xyzzzzz_0_yzzzzzz_0, g_xyzzzzz_0_zzzzzzz_0, g_xzzzzz_0_xxxxxxx_1, g_xzzzzz_0_xxxxxxz_1, g_xzzzzz_0_xxxxxzz_1, g_xzzzzz_0_xxxxzzz_1, g_xzzzzz_0_xxxzzzz_1, g_xzzzzz_0_xxzzzzz_1, g_xzzzzz_0_xzzzzzz_1, g_yzzzzz_0_xxxxxxy_1, g_yzzzzz_0_xxxxxy_1, g_yzzzzz_0_xxxxxyy_1, g_yzzzzz_0_xxxxxyz_1, g_yzzzzz_0_xxxxyy_1, g_yzzzzz_0_xxxxyyy_1, g_yzzzzz_0_xxxxyyz_1, g_yzzzzz_0_xxxxyz_1, g_yzzzzz_0_xxxxyzz_1, g_yzzzzz_0_xxxyyy_1, g_yzzzzz_0_xxxyyyy_1, g_yzzzzz_0_xxxyyyz_1, g_yzzzzz_0_xxxyyz_1, g_yzzzzz_0_xxxyyzz_1, g_yzzzzz_0_xxxyzz_1, g_yzzzzz_0_xxxyzzz_1, g_yzzzzz_0_xxyyyy_1, g_yzzzzz_0_xxyyyyy_1, g_yzzzzz_0_xxyyyyz_1, g_yzzzzz_0_xxyyyz_1, g_yzzzzz_0_xxyyyzz_1, g_yzzzzz_0_xxyyzz_1, g_yzzzzz_0_xxyyzzz_1, g_yzzzzz_0_xxyzzz_1, g_yzzzzz_0_xxyzzzz_1, g_yzzzzz_0_xyyyyy_1, g_yzzzzz_0_xyyyyyy_1, g_yzzzzz_0_xyyyyyz_1, g_yzzzzz_0_xyyyyz_1, g_yzzzzz_0_xyyyyzz_1, g_yzzzzz_0_xyyyzz_1, g_yzzzzz_0_xyyyzzz_1, g_yzzzzz_0_xyyzzz_1, g_yzzzzz_0_xyyzzzz_1, g_yzzzzz_0_xyzzzz_1, g_yzzzzz_0_xyzzzzz_1, g_yzzzzz_0_yyyyyy_1, g_yzzzzz_0_yyyyyyy_1, g_yzzzzz_0_yyyyyyz_1, g_yzzzzz_0_yyyyyz_1, g_yzzzzz_0_yyyyyzz_1, g_yzzzzz_0_yyyyzz_1, g_yzzzzz_0_yyyyzzz_1, g_yzzzzz_0_yyyzzz_1, g_yzzzzz_0_yyyzzzz_1, g_yzzzzz_0_yyzzzz_1, g_yzzzzz_0_yyzzzzz_1, g_yzzzzz_0_yzzzzz_1, g_yzzzzz_0_yzzzzzz_1, g_yzzzzz_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzzzz_0_xxxxxxx_0[i] = g_xzzzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxxxxy_0[i] = 6.0 * g_yzzzzz_0_xxxxxy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxxxz_0[i] = g_xzzzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxxxyy_0[i] = 5.0 * g_yzzzzz_0_xxxxyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxxyz_0[i] = 5.0 * g_yzzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxxzz_0[i] = g_xzzzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxxyyy_0[i] = 4.0 * g_yzzzzz_0_xxxyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxyyz_0[i] = 4.0 * g_yzzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxyzz_0[i] = 4.0 * g_yzzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxxzzz_0[i] = g_xzzzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxxyyyy_0[i] = 3.0 * g_yzzzzz_0_xxyyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxyyyz_0[i] = 3.0 * g_yzzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxyyzz_0[i] = 3.0 * g_yzzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxyzzz_0[i] = 3.0 * g_yzzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxxzzzz_0[i] = g_xzzzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xxyyyyy_0[i] = 2.0 * g_yzzzzz_0_xyyyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyyyyz_0[i] = 2.0 * g_yzzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyyyzz_0[i] = 2.0 * g_yzzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyyzzz_0[i] = 2.0 * g_yzzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxyzzzz_0[i] = 2.0 * g_yzzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xxzzzzz_0[i] = g_xzzzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_xyyyyyy_0[i] = g_yzzzzz_0_yyyyyy_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyyyyz_0[i] = g_yzzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyyyzz_0[i] = g_yzzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyyzzz_0[i] = g_yzzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyyzzzz_0[i] = g_yzzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xyzzzzz_0[i] = g_yzzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_xzzzzzz_0[i] = g_xzzzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xyzzzzz_0_yyyyyyy_0[i] = g_yzzzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyyyyz_0[i] = g_yzzzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyyyzz_0[i] = g_yzzzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyyzzz_0[i] = g_yzzzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyyzzzz_0[i] = g_yzzzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yyzzzzz_0[i] = g_yzzzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_yzzzzzz_0[i] = g_yzzzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyzzzzz_0_zzzzzzz_0[i] = g_yzzzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 972-1008 components of targeted buffer : KSK

    auto g_xzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 972);

    auto g_xzzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 973);

    auto g_xzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 974);

    auto g_xzzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 975);

    auto g_xzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 976);

    auto g_xzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 977);

    auto g_xzzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 978);

    auto g_xzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 979);

    auto g_xzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 980);

    auto g_xzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 981);

    auto g_xzzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 982);

    auto g_xzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 983);

    auto g_xzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 984);

    auto g_xzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 985);

    auto g_xzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 986);

    auto g_xzzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 987);

    auto g_xzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 988);

    auto g_xzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 989);

    auto g_xzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 990);

    auto g_xzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 991);

    auto g_xzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 992);

    auto g_xzzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 993);

    auto g_xzzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 994);

    auto g_xzzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 995);

    auto g_xzzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 996);

    auto g_xzzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 997);

    auto g_xzzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 998);

    auto g_xzzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 999);

    auto g_xzzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1000);

    auto g_xzzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1001);

    auto g_xzzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1002);

    auto g_xzzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1003);

    auto g_xzzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1004);

    auto g_xzzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1005);

    auto g_xzzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1006);

    auto g_xzzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1007);

    #pragma omp simd aligned(g_xzzzzzz_0_xxxxxxx_0, g_xzzzzzz_0_xxxxxxy_0, g_xzzzzzz_0_xxxxxxz_0, g_xzzzzzz_0_xxxxxyy_0, g_xzzzzzz_0_xxxxxyz_0, g_xzzzzzz_0_xxxxxzz_0, g_xzzzzzz_0_xxxxyyy_0, g_xzzzzzz_0_xxxxyyz_0, g_xzzzzzz_0_xxxxyzz_0, g_xzzzzzz_0_xxxxzzz_0, g_xzzzzzz_0_xxxyyyy_0, g_xzzzzzz_0_xxxyyyz_0, g_xzzzzzz_0_xxxyyzz_0, g_xzzzzzz_0_xxxyzzz_0, g_xzzzzzz_0_xxxzzzz_0, g_xzzzzzz_0_xxyyyyy_0, g_xzzzzzz_0_xxyyyyz_0, g_xzzzzzz_0_xxyyyzz_0, g_xzzzzzz_0_xxyyzzz_0, g_xzzzzzz_0_xxyzzzz_0, g_xzzzzzz_0_xxzzzzz_0, g_xzzzzzz_0_xyyyyyy_0, g_xzzzzzz_0_xyyyyyz_0, g_xzzzzzz_0_xyyyyzz_0, g_xzzzzzz_0_xyyyzzz_0, g_xzzzzzz_0_xyyzzzz_0, g_xzzzzzz_0_xyzzzzz_0, g_xzzzzzz_0_xzzzzzz_0, g_xzzzzzz_0_yyyyyyy_0, g_xzzzzzz_0_yyyyyyz_0, g_xzzzzzz_0_yyyyyzz_0, g_xzzzzzz_0_yyyyzzz_0, g_xzzzzzz_0_yyyzzzz_0, g_xzzzzzz_0_yyzzzzz_0, g_xzzzzzz_0_yzzzzzz_0, g_xzzzzzz_0_zzzzzzz_0, g_zzzzzz_0_xxxxxx_1, g_zzzzzz_0_xxxxxxx_1, g_zzzzzz_0_xxxxxxy_1, g_zzzzzz_0_xxxxxxz_1, g_zzzzzz_0_xxxxxy_1, g_zzzzzz_0_xxxxxyy_1, g_zzzzzz_0_xxxxxyz_1, g_zzzzzz_0_xxxxxz_1, g_zzzzzz_0_xxxxxzz_1, g_zzzzzz_0_xxxxyy_1, g_zzzzzz_0_xxxxyyy_1, g_zzzzzz_0_xxxxyyz_1, g_zzzzzz_0_xxxxyz_1, g_zzzzzz_0_xxxxyzz_1, g_zzzzzz_0_xxxxzz_1, g_zzzzzz_0_xxxxzzz_1, g_zzzzzz_0_xxxyyy_1, g_zzzzzz_0_xxxyyyy_1, g_zzzzzz_0_xxxyyyz_1, g_zzzzzz_0_xxxyyz_1, g_zzzzzz_0_xxxyyzz_1, g_zzzzzz_0_xxxyzz_1, g_zzzzzz_0_xxxyzzz_1, g_zzzzzz_0_xxxzzz_1, g_zzzzzz_0_xxxzzzz_1, g_zzzzzz_0_xxyyyy_1, g_zzzzzz_0_xxyyyyy_1, g_zzzzzz_0_xxyyyyz_1, g_zzzzzz_0_xxyyyz_1, g_zzzzzz_0_xxyyyzz_1, g_zzzzzz_0_xxyyzz_1, g_zzzzzz_0_xxyyzzz_1, g_zzzzzz_0_xxyzzz_1, g_zzzzzz_0_xxyzzzz_1, g_zzzzzz_0_xxzzzz_1, g_zzzzzz_0_xxzzzzz_1, g_zzzzzz_0_xyyyyy_1, g_zzzzzz_0_xyyyyyy_1, g_zzzzzz_0_xyyyyyz_1, g_zzzzzz_0_xyyyyz_1, g_zzzzzz_0_xyyyyzz_1, g_zzzzzz_0_xyyyzz_1, g_zzzzzz_0_xyyyzzz_1, g_zzzzzz_0_xyyzzz_1, g_zzzzzz_0_xyyzzzz_1, g_zzzzzz_0_xyzzzz_1, g_zzzzzz_0_xyzzzzz_1, g_zzzzzz_0_xzzzzz_1, g_zzzzzz_0_xzzzzzz_1, g_zzzzzz_0_yyyyyy_1, g_zzzzzz_0_yyyyyyy_1, g_zzzzzz_0_yyyyyyz_1, g_zzzzzz_0_yyyyyz_1, g_zzzzzz_0_yyyyyzz_1, g_zzzzzz_0_yyyyzz_1, g_zzzzzz_0_yyyyzzz_1, g_zzzzzz_0_yyyzzz_1, g_zzzzzz_0_yyyzzzz_1, g_zzzzzz_0_yyzzzz_1, g_zzzzzz_0_yyzzzzz_1, g_zzzzzz_0_yzzzzz_1, g_zzzzzz_0_yzzzzzz_1, g_zzzzzz_0_zzzzzz_1, g_zzzzzz_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzzz_0_xxxxxxx_0[i] = 7.0 * g_zzzzzz_0_xxxxxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxx_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxxy_0[i] = 6.0 * g_zzzzzz_0_xxxxxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxxz_0[i] = 6.0 * g_zzzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxyy_0[i] = 5.0 * g_zzzzzz_0_xxxxyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxyz_0[i] = 5.0 * g_zzzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxxzz_0[i] = 5.0 * g_zzzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxyyy_0[i] = 4.0 * g_zzzzzz_0_xxxyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxyyz_0[i] = 4.0 * g_zzzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxyzz_0[i] = 4.0 * g_zzzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxxzzz_0[i] = 4.0 * g_zzzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyyyy_0[i] = 3.0 * g_zzzzzz_0_xxyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyyyz_0[i] = 3.0 * g_zzzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyyzz_0[i] = 3.0 * g_zzzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxyzzz_0[i] = 3.0 * g_zzzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxxzzzz_0[i] = 3.0 * g_zzzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyyyy_0[i] = 2.0 * g_zzzzzz_0_xyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyyyz_0[i] = 2.0 * g_zzzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyyzz_0[i] = 2.0 * g_zzzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyyzzz_0[i] = 2.0 * g_zzzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxyzzzz_0[i] = 2.0 * g_zzzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xxzzzzz_0[i] = 2.0 * g_zzzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyyyy_0[i] = g_zzzzzz_0_yyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyyyz_0[i] = g_zzzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyyzz_0[i] = g_zzzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyyzzz_0[i] = g_zzzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyyzzzz_0[i] = g_zzzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xyzzzzz_0[i] = g_zzzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_xzzzzzz_0[i] = g_zzzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyyyy_0[i] = g_zzzzzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyyyz_0[i] = g_zzzzzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyyzz_0[i] = g_zzzzzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyyzzz_0[i] = g_zzzzzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyyzzzz_0[i] = g_zzzzzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yyzzzzz_0[i] = g_zzzzzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_yzzzzzz_0[i] = g_zzzzzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xzzzzzz_0_zzzzzzz_0[i] = g_zzzzzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 1008-1044 components of targeted buffer : KSK

    auto g_yyyyyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 1008);

    auto g_yyyyyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 1009);

    auto g_yyyyyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 1010);

    auto g_yyyyyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 1011);

    auto g_yyyyyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 1012);

    auto g_yyyyyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 1013);

    auto g_yyyyyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 1014);

    auto g_yyyyyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 1015);

    auto g_yyyyyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 1016);

    auto g_yyyyyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 1017);

    auto g_yyyyyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1018);

    auto g_yyyyyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1019);

    auto g_yyyyyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1020);

    auto g_yyyyyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1021);

    auto g_yyyyyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1022);

    auto g_yyyyyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1023);

    auto g_yyyyyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1024);

    auto g_yyyyyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1025);

    auto g_yyyyyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1026);

    auto g_yyyyyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1027);

    auto g_yyyyyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1028);

    auto g_yyyyyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1029);

    auto g_yyyyyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1030);

    auto g_yyyyyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1031);

    auto g_yyyyyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1032);

    auto g_yyyyyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1033);

    auto g_yyyyyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1034);

    auto g_yyyyyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1035);

    auto g_yyyyyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1036);

    auto g_yyyyyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1037);

    auto g_yyyyyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1038);

    auto g_yyyyyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1039);

    auto g_yyyyyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1040);

    auto g_yyyyyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1041);

    auto g_yyyyyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1042);

    auto g_yyyyyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1043);

    #pragma omp simd aligned(g_yyyyy_0_xxxxxxx_0, g_yyyyy_0_xxxxxxx_1, g_yyyyy_0_xxxxxxy_0, g_yyyyy_0_xxxxxxy_1, g_yyyyy_0_xxxxxxz_0, g_yyyyy_0_xxxxxxz_1, g_yyyyy_0_xxxxxyy_0, g_yyyyy_0_xxxxxyy_1, g_yyyyy_0_xxxxxyz_0, g_yyyyy_0_xxxxxyz_1, g_yyyyy_0_xxxxxzz_0, g_yyyyy_0_xxxxxzz_1, g_yyyyy_0_xxxxyyy_0, g_yyyyy_0_xxxxyyy_1, g_yyyyy_0_xxxxyyz_0, g_yyyyy_0_xxxxyyz_1, g_yyyyy_0_xxxxyzz_0, g_yyyyy_0_xxxxyzz_1, g_yyyyy_0_xxxxzzz_0, g_yyyyy_0_xxxxzzz_1, g_yyyyy_0_xxxyyyy_0, g_yyyyy_0_xxxyyyy_1, g_yyyyy_0_xxxyyyz_0, g_yyyyy_0_xxxyyyz_1, g_yyyyy_0_xxxyyzz_0, g_yyyyy_0_xxxyyzz_1, g_yyyyy_0_xxxyzzz_0, g_yyyyy_0_xxxyzzz_1, g_yyyyy_0_xxxzzzz_0, g_yyyyy_0_xxxzzzz_1, g_yyyyy_0_xxyyyyy_0, g_yyyyy_0_xxyyyyy_1, g_yyyyy_0_xxyyyyz_0, g_yyyyy_0_xxyyyyz_1, g_yyyyy_0_xxyyyzz_0, g_yyyyy_0_xxyyyzz_1, g_yyyyy_0_xxyyzzz_0, g_yyyyy_0_xxyyzzz_1, g_yyyyy_0_xxyzzzz_0, g_yyyyy_0_xxyzzzz_1, g_yyyyy_0_xxzzzzz_0, g_yyyyy_0_xxzzzzz_1, g_yyyyy_0_xyyyyyy_0, g_yyyyy_0_xyyyyyy_1, g_yyyyy_0_xyyyyyz_0, g_yyyyy_0_xyyyyyz_1, g_yyyyy_0_xyyyyzz_0, g_yyyyy_0_xyyyyzz_1, g_yyyyy_0_xyyyzzz_0, g_yyyyy_0_xyyyzzz_1, g_yyyyy_0_xyyzzzz_0, g_yyyyy_0_xyyzzzz_1, g_yyyyy_0_xyzzzzz_0, g_yyyyy_0_xyzzzzz_1, g_yyyyy_0_xzzzzzz_0, g_yyyyy_0_xzzzzzz_1, g_yyyyy_0_yyyyyyy_0, g_yyyyy_0_yyyyyyy_1, g_yyyyy_0_yyyyyyz_0, g_yyyyy_0_yyyyyyz_1, g_yyyyy_0_yyyyyzz_0, g_yyyyy_0_yyyyyzz_1, g_yyyyy_0_yyyyzzz_0, g_yyyyy_0_yyyyzzz_1, g_yyyyy_0_yyyzzzz_0, g_yyyyy_0_yyyzzzz_1, g_yyyyy_0_yyzzzzz_0, g_yyyyy_0_yyzzzzz_1, g_yyyyy_0_yzzzzzz_0, g_yyyyy_0_yzzzzzz_1, g_yyyyy_0_zzzzzzz_0, g_yyyyy_0_zzzzzzz_1, g_yyyyyy_0_xxxxxx_1, g_yyyyyy_0_xxxxxxx_1, g_yyyyyy_0_xxxxxxy_1, g_yyyyyy_0_xxxxxxz_1, g_yyyyyy_0_xxxxxy_1, g_yyyyyy_0_xxxxxyy_1, g_yyyyyy_0_xxxxxyz_1, g_yyyyyy_0_xxxxxz_1, g_yyyyyy_0_xxxxxzz_1, g_yyyyyy_0_xxxxyy_1, g_yyyyyy_0_xxxxyyy_1, g_yyyyyy_0_xxxxyyz_1, g_yyyyyy_0_xxxxyz_1, g_yyyyyy_0_xxxxyzz_1, g_yyyyyy_0_xxxxzz_1, g_yyyyyy_0_xxxxzzz_1, g_yyyyyy_0_xxxyyy_1, g_yyyyyy_0_xxxyyyy_1, g_yyyyyy_0_xxxyyyz_1, g_yyyyyy_0_xxxyyz_1, g_yyyyyy_0_xxxyyzz_1, g_yyyyyy_0_xxxyzz_1, g_yyyyyy_0_xxxyzzz_1, g_yyyyyy_0_xxxzzz_1, g_yyyyyy_0_xxxzzzz_1, g_yyyyyy_0_xxyyyy_1, g_yyyyyy_0_xxyyyyy_1, g_yyyyyy_0_xxyyyyz_1, g_yyyyyy_0_xxyyyz_1, g_yyyyyy_0_xxyyyzz_1, g_yyyyyy_0_xxyyzz_1, g_yyyyyy_0_xxyyzzz_1, g_yyyyyy_0_xxyzzz_1, g_yyyyyy_0_xxyzzzz_1, g_yyyyyy_0_xxzzzz_1, g_yyyyyy_0_xxzzzzz_1, g_yyyyyy_0_xyyyyy_1, g_yyyyyy_0_xyyyyyy_1, g_yyyyyy_0_xyyyyyz_1, g_yyyyyy_0_xyyyyz_1, g_yyyyyy_0_xyyyyzz_1, g_yyyyyy_0_xyyyzz_1, g_yyyyyy_0_xyyyzzz_1, g_yyyyyy_0_xyyzzz_1, g_yyyyyy_0_xyyzzzz_1, g_yyyyyy_0_xyzzzz_1, g_yyyyyy_0_xyzzzzz_1, g_yyyyyy_0_xzzzzz_1, g_yyyyyy_0_xzzzzzz_1, g_yyyyyy_0_yyyyyy_1, g_yyyyyy_0_yyyyyyy_1, g_yyyyyy_0_yyyyyyz_1, g_yyyyyy_0_yyyyyz_1, g_yyyyyy_0_yyyyyzz_1, g_yyyyyy_0_yyyyzz_1, g_yyyyyy_0_yyyyzzz_1, g_yyyyyy_0_yyyzzz_1, g_yyyyyy_0_yyyzzzz_1, g_yyyyyy_0_yyzzzz_1, g_yyyyyy_0_yyzzzzz_1, g_yyyyyy_0_yzzzzz_1, g_yyyyyy_0_yzzzzzz_1, g_yyyyyy_0_zzzzzz_1, g_yyyyyy_0_zzzzzzz_1, g_yyyyyyy_0_xxxxxxx_0, g_yyyyyyy_0_xxxxxxy_0, g_yyyyyyy_0_xxxxxxz_0, g_yyyyyyy_0_xxxxxyy_0, g_yyyyyyy_0_xxxxxyz_0, g_yyyyyyy_0_xxxxxzz_0, g_yyyyyyy_0_xxxxyyy_0, g_yyyyyyy_0_xxxxyyz_0, g_yyyyyyy_0_xxxxyzz_0, g_yyyyyyy_0_xxxxzzz_0, g_yyyyyyy_0_xxxyyyy_0, g_yyyyyyy_0_xxxyyyz_0, g_yyyyyyy_0_xxxyyzz_0, g_yyyyyyy_0_xxxyzzz_0, g_yyyyyyy_0_xxxzzzz_0, g_yyyyyyy_0_xxyyyyy_0, g_yyyyyyy_0_xxyyyyz_0, g_yyyyyyy_0_xxyyyzz_0, g_yyyyyyy_0_xxyyzzz_0, g_yyyyyyy_0_xxyzzzz_0, g_yyyyyyy_0_xxzzzzz_0, g_yyyyyyy_0_xyyyyyy_0, g_yyyyyyy_0_xyyyyyz_0, g_yyyyyyy_0_xyyyyzz_0, g_yyyyyyy_0_xyyyzzz_0, g_yyyyyyy_0_xyyzzzz_0, g_yyyyyyy_0_xyzzzzz_0, g_yyyyyyy_0_xzzzzzz_0, g_yyyyyyy_0_yyyyyyy_0, g_yyyyyyy_0_yyyyyyz_0, g_yyyyyyy_0_yyyyyzz_0, g_yyyyyyy_0_yyyyzzz_0, g_yyyyyyy_0_yyyzzzz_0, g_yyyyyyy_0_yyzzzzz_0, g_yyyyyyy_0_yzzzzzz_0, g_yyyyyyy_0_zzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyyy_0_xxxxxxx_0[i] = 6.0 * g_yyyyy_0_xxxxxxx_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxxx_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxxx_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxxy_0[i] = 6.0 * g_yyyyy_0_xxxxxxy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxxy_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxxz_0[i] = 6.0 * g_yyyyy_0_xxxxxxz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxxz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxxz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxyy_0[i] = 6.0 * g_yyyyy_0_xxxxxyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxyy_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxxxxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxyz_0[i] = 6.0 * g_yyyyy_0_xxxxxyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxyz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxxzz_0[i] = 6.0 * g_yyyyy_0_xxxxxzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxxzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxxzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxyyy_0[i] = 6.0 * g_yyyyy_0_xxxxyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxyyy_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xxxxyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxyyz_0[i] = 6.0 * g_yyyyy_0_xxxxyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxxxyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxyzz_0[i] = 6.0 * g_yyyyy_0_xxxxyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxyzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxxzzz_0[i] = 6.0 * g_yyyyy_0_xxxxzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxxzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxxzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyyyy_0[i] = 6.0 * g_yyyyy_0_xxxyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyyyy_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_xxxyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyyyz_0[i] = 6.0 * g_yyyyy_0_xxxyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xxxyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyyzz_0[i] = 6.0 * g_yyyyy_0_xxxyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxxyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxyzzz_0[i] = 6.0 * g_yyyyy_0_xxxyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxyzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxxzzzz_0[i] = 6.0 * g_yyyyy_0_xxxzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxxzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxxzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyyyy_0[i] = 6.0 * g_yyyyy_0_xxyyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyyyy_1[i] * fz_be_0 + 5.0 * g_yyyyyy_0_xxyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyyyz_0[i] = 6.0 * g_yyyyy_0_xxyyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_xxyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyyzz_0[i] = 6.0 * g_yyyyy_0_xxyyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xxyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyyzzz_0[i] = 6.0 * g_yyyyy_0_xxyyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xxyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxyzzzz_0[i] = 6.0 * g_yyyyy_0_xxyzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxyzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xxzzzzz_0[i] = 6.0 * g_yyyyy_0_xxzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xxzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xxzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyyyy_0[i] = 6.0 * g_yyyyy_0_xyyyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyyyy_1[i] * fz_be_0 + 6.0 * g_yyyyyy_0_xyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyyyz_0[i] = 6.0 * g_yyyyy_0_xyyyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyyyy_0_xyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyyzz_0[i] = 6.0 * g_yyyyy_0_xyyyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_xyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyyzzz_0[i] = 6.0 * g_yyyyy_0_xyyyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_xyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyyzzzz_0[i] = 6.0 * g_yyyyy_0_xyyzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_xyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xyzzzzz_0[i] = 6.0 * g_yyyyy_0_xyzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xyzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_xzzzzzz_0[i] = 6.0 * g_yyyyy_0_xzzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_xzzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_xzzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyyyy_0[i] = 6.0 * g_yyyyy_0_yyyyyyy_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyyyy_1[i] * fz_be_0 + 7.0 * g_yyyyyy_0_yyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyyy_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyyyz_0[i] = 6.0 * g_yyyyy_0_yyyyyyz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyyyyy_0_yyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyyz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyyzz_0[i] = 6.0 * g_yyyyy_0_yyyyyzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyyyyy_0_yyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyyzzz_0[i] = 6.0 * g_yyyyy_0_yyyyzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyyyyy_0_yyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyyzzzz_0[i] = 6.0 * g_yyyyy_0_yyyzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyyyyy_0_yyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yyzzzzz_0[i] = 6.0 * g_yyyyy_0_yyzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_0_yzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_yzzzzzz_0[i] = 6.0 * g_yyyyy_0_yzzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_yzzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_zzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yzzzzzz_1[i] * wa_y[i];

        g_yyyyyyy_0_zzzzzzz_0[i] = 6.0 * g_yyyyy_0_zzzzzzz_0[i] * fbe_0 - 6.0 * g_yyyyy_0_zzzzzzz_1[i] * fz_be_0 + g_yyyyyy_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1044-1080 components of targeted buffer : KSK

    auto g_yyyyyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 1044);

    auto g_yyyyyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 1045);

    auto g_yyyyyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 1046);

    auto g_yyyyyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 1047);

    auto g_yyyyyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 1048);

    auto g_yyyyyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 1049);

    auto g_yyyyyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 1050);

    auto g_yyyyyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 1051);

    auto g_yyyyyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 1052);

    auto g_yyyyyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 1053);

    auto g_yyyyyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1054);

    auto g_yyyyyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1055);

    auto g_yyyyyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1056);

    auto g_yyyyyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1057);

    auto g_yyyyyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1058);

    auto g_yyyyyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1059);

    auto g_yyyyyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1060);

    auto g_yyyyyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1061);

    auto g_yyyyyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1062);

    auto g_yyyyyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1063);

    auto g_yyyyyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1064);

    auto g_yyyyyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1065);

    auto g_yyyyyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1066);

    auto g_yyyyyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1067);

    auto g_yyyyyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1068);

    auto g_yyyyyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1069);

    auto g_yyyyyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1070);

    auto g_yyyyyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1071);

    auto g_yyyyyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1072);

    auto g_yyyyyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1073);

    auto g_yyyyyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1074);

    auto g_yyyyyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1075);

    auto g_yyyyyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1076);

    auto g_yyyyyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1077);

    auto g_yyyyyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1078);

    auto g_yyyyyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1079);

    #pragma omp simd aligned(g_yyyyyy_0_xxxxxx_1, g_yyyyyy_0_xxxxxxx_1, g_yyyyyy_0_xxxxxxy_1, g_yyyyyy_0_xxxxxxz_1, g_yyyyyy_0_xxxxxy_1, g_yyyyyy_0_xxxxxyy_1, g_yyyyyy_0_xxxxxyz_1, g_yyyyyy_0_xxxxxz_1, g_yyyyyy_0_xxxxxzz_1, g_yyyyyy_0_xxxxyy_1, g_yyyyyy_0_xxxxyyy_1, g_yyyyyy_0_xxxxyyz_1, g_yyyyyy_0_xxxxyz_1, g_yyyyyy_0_xxxxyzz_1, g_yyyyyy_0_xxxxzz_1, g_yyyyyy_0_xxxxzzz_1, g_yyyyyy_0_xxxyyy_1, g_yyyyyy_0_xxxyyyy_1, g_yyyyyy_0_xxxyyyz_1, g_yyyyyy_0_xxxyyz_1, g_yyyyyy_0_xxxyyzz_1, g_yyyyyy_0_xxxyzz_1, g_yyyyyy_0_xxxyzzz_1, g_yyyyyy_0_xxxzzz_1, g_yyyyyy_0_xxxzzzz_1, g_yyyyyy_0_xxyyyy_1, g_yyyyyy_0_xxyyyyy_1, g_yyyyyy_0_xxyyyyz_1, g_yyyyyy_0_xxyyyz_1, g_yyyyyy_0_xxyyyzz_1, g_yyyyyy_0_xxyyzz_1, g_yyyyyy_0_xxyyzzz_1, g_yyyyyy_0_xxyzzz_1, g_yyyyyy_0_xxyzzzz_1, g_yyyyyy_0_xxzzzz_1, g_yyyyyy_0_xxzzzzz_1, g_yyyyyy_0_xyyyyy_1, g_yyyyyy_0_xyyyyyy_1, g_yyyyyy_0_xyyyyyz_1, g_yyyyyy_0_xyyyyz_1, g_yyyyyy_0_xyyyyzz_1, g_yyyyyy_0_xyyyzz_1, g_yyyyyy_0_xyyyzzz_1, g_yyyyyy_0_xyyzzz_1, g_yyyyyy_0_xyyzzzz_1, g_yyyyyy_0_xyzzzz_1, g_yyyyyy_0_xyzzzzz_1, g_yyyyyy_0_xzzzzz_1, g_yyyyyy_0_xzzzzzz_1, g_yyyyyy_0_yyyyyy_1, g_yyyyyy_0_yyyyyyy_1, g_yyyyyy_0_yyyyyyz_1, g_yyyyyy_0_yyyyyz_1, g_yyyyyy_0_yyyyyzz_1, g_yyyyyy_0_yyyyzz_1, g_yyyyyy_0_yyyyzzz_1, g_yyyyyy_0_yyyzzz_1, g_yyyyyy_0_yyyzzzz_1, g_yyyyyy_0_yyzzzz_1, g_yyyyyy_0_yyzzzzz_1, g_yyyyyy_0_yzzzzz_1, g_yyyyyy_0_yzzzzzz_1, g_yyyyyy_0_zzzzzz_1, g_yyyyyy_0_zzzzzzz_1, g_yyyyyyz_0_xxxxxxx_0, g_yyyyyyz_0_xxxxxxy_0, g_yyyyyyz_0_xxxxxxz_0, g_yyyyyyz_0_xxxxxyy_0, g_yyyyyyz_0_xxxxxyz_0, g_yyyyyyz_0_xxxxxzz_0, g_yyyyyyz_0_xxxxyyy_0, g_yyyyyyz_0_xxxxyyz_0, g_yyyyyyz_0_xxxxyzz_0, g_yyyyyyz_0_xxxxzzz_0, g_yyyyyyz_0_xxxyyyy_0, g_yyyyyyz_0_xxxyyyz_0, g_yyyyyyz_0_xxxyyzz_0, g_yyyyyyz_0_xxxyzzz_0, g_yyyyyyz_0_xxxzzzz_0, g_yyyyyyz_0_xxyyyyy_0, g_yyyyyyz_0_xxyyyyz_0, g_yyyyyyz_0_xxyyyzz_0, g_yyyyyyz_0_xxyyzzz_0, g_yyyyyyz_0_xxyzzzz_0, g_yyyyyyz_0_xxzzzzz_0, g_yyyyyyz_0_xyyyyyy_0, g_yyyyyyz_0_xyyyyyz_0, g_yyyyyyz_0_xyyyyzz_0, g_yyyyyyz_0_xyyyzzz_0, g_yyyyyyz_0_xyyzzzz_0, g_yyyyyyz_0_xyzzzzz_0, g_yyyyyyz_0_xzzzzzz_0, g_yyyyyyz_0_yyyyyyy_0, g_yyyyyyz_0_yyyyyyz_0, g_yyyyyyz_0_yyyyyzz_0, g_yyyyyyz_0_yyyyzzz_0, g_yyyyyyz_0_yyyzzzz_0, g_yyyyyyz_0_yyzzzzz_0, g_yyyyyyz_0_yzzzzzz_0, g_yyyyyyz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyyz_0_xxxxxxx_0[i] = g_yyyyyy_0_xxxxxxx_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxxy_0[i] = g_yyyyyy_0_xxxxxxy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxxz_0[i] = g_yyyyyy_0_xxxxxx_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxxz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxyy_0[i] = g_yyyyyy_0_xxxxxyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxyz_0[i] = g_yyyyyy_0_xxxxxy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxxzz_0[i] = 2.0 * g_yyyyyy_0_xxxxxz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxxzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxyyy_0[i] = g_yyyyyy_0_xxxxyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxyyz_0[i] = g_yyyyyy_0_xxxxyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxyzz_0[i] = 2.0 * g_yyyyyy_0_xxxxyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxxzzz_0[i] = 3.0 * g_yyyyyy_0_xxxxzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxxzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyyyy_0[i] = g_yyyyyy_0_xxxyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyyyz_0[i] = g_yyyyyy_0_xxxyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyyzz_0[i] = 2.0 * g_yyyyyy_0_xxxyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxyzzz_0[i] = 3.0 * g_yyyyyy_0_xxxyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxxzzzz_0[i] = 4.0 * g_yyyyyy_0_xxxzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxxzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyyyy_0[i] = g_yyyyyy_0_xxyyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyyyz_0[i] = g_yyyyyy_0_xxyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyyzz_0[i] = 2.0 * g_yyyyyy_0_xxyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyyzzz_0[i] = 3.0 * g_yyyyyy_0_xxyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxyzzzz_0[i] = 4.0 * g_yyyyyy_0_xxyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxyzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xxzzzzz_0[i] = 5.0 * g_yyyyyy_0_xxzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xxzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyyyy_0[i] = g_yyyyyy_0_xyyyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyyyz_0[i] = g_yyyyyy_0_xyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyyzz_0[i] = 2.0 * g_yyyyyy_0_xyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyyzzz_0[i] = 3.0 * g_yyyyyy_0_xyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyyzzzz_0[i] = 4.0 * g_yyyyyy_0_xyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyyzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xyzzzzz_0[i] = 5.0 * g_yyyyyy_0_xyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xyzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_xzzzzzz_0[i] = 6.0 * g_yyyyyy_0_xzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_xzzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyyyy_0[i] = g_yyyyyy_0_yyyyyyy_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyyyz_0[i] = g_yyyyyy_0_yyyyyy_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyyz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyyzz_0[i] = 2.0 * g_yyyyyy_0_yyyyyz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyyzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyyzzz_0[i] = 3.0 * g_yyyyyy_0_yyyyzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyyzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyyzzzz_0[i] = 4.0 * g_yyyyyy_0_yyyzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyyzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yyzzzzz_0[i] = 5.0 * g_yyyyyy_0_yyzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yyzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_yzzzzzz_0[i] = 6.0 * g_yyyyyy_0_yzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_yzzzzzz_1[i] * wa_z[i];

        g_yyyyyyz_0_zzzzzzz_0[i] = 7.0 * g_yyyyyy_0_zzzzzz_1[i] * fi_acd_0 + g_yyyyyy_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 1080-1116 components of targeted buffer : KSK

    auto g_yyyyyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 1080);

    auto g_yyyyyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 1081);

    auto g_yyyyyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 1082);

    auto g_yyyyyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 1083);

    auto g_yyyyyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 1084);

    auto g_yyyyyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 1085);

    auto g_yyyyyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 1086);

    auto g_yyyyyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 1087);

    auto g_yyyyyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 1088);

    auto g_yyyyyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 1089);

    auto g_yyyyyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1090);

    auto g_yyyyyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1091);

    auto g_yyyyyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1092);

    auto g_yyyyyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1093);

    auto g_yyyyyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1094);

    auto g_yyyyyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1095);

    auto g_yyyyyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1096);

    auto g_yyyyyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1097);

    auto g_yyyyyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1098);

    auto g_yyyyyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1099);

    auto g_yyyyyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1100);

    auto g_yyyyyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1101);

    auto g_yyyyyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1102);

    auto g_yyyyyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1103);

    auto g_yyyyyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1104);

    auto g_yyyyyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1105);

    auto g_yyyyyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1106);

    auto g_yyyyyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1107);

    auto g_yyyyyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1108);

    auto g_yyyyyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1109);

    auto g_yyyyyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1110);

    auto g_yyyyyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1111);

    auto g_yyyyyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1112);

    auto g_yyyyyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1113);

    auto g_yyyyyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1114);

    auto g_yyyyyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1115);

    #pragma omp simd aligned(g_yyyyy_0_xxxxxxy_0, g_yyyyy_0_xxxxxxy_1, g_yyyyy_0_xxxxxyy_0, g_yyyyy_0_xxxxxyy_1, g_yyyyy_0_xxxxyyy_0, g_yyyyy_0_xxxxyyy_1, g_yyyyy_0_xxxyyyy_0, g_yyyyy_0_xxxyyyy_1, g_yyyyy_0_xxyyyyy_0, g_yyyyy_0_xxyyyyy_1, g_yyyyy_0_xyyyyyy_0, g_yyyyy_0_xyyyyyy_1, g_yyyyy_0_yyyyyyy_0, g_yyyyy_0_yyyyyyy_1, g_yyyyyz_0_xxxxxxy_1, g_yyyyyz_0_xxxxxyy_1, g_yyyyyz_0_xxxxyyy_1, g_yyyyyz_0_xxxyyyy_1, g_yyyyyz_0_xxyyyyy_1, g_yyyyyz_0_xyyyyyy_1, g_yyyyyz_0_yyyyyyy_1, g_yyyyyzz_0_xxxxxxx_0, g_yyyyyzz_0_xxxxxxy_0, g_yyyyyzz_0_xxxxxxz_0, g_yyyyyzz_0_xxxxxyy_0, g_yyyyyzz_0_xxxxxyz_0, g_yyyyyzz_0_xxxxxzz_0, g_yyyyyzz_0_xxxxyyy_0, g_yyyyyzz_0_xxxxyyz_0, g_yyyyyzz_0_xxxxyzz_0, g_yyyyyzz_0_xxxxzzz_0, g_yyyyyzz_0_xxxyyyy_0, g_yyyyyzz_0_xxxyyyz_0, g_yyyyyzz_0_xxxyyzz_0, g_yyyyyzz_0_xxxyzzz_0, g_yyyyyzz_0_xxxzzzz_0, g_yyyyyzz_0_xxyyyyy_0, g_yyyyyzz_0_xxyyyyz_0, g_yyyyyzz_0_xxyyyzz_0, g_yyyyyzz_0_xxyyzzz_0, g_yyyyyzz_0_xxyzzzz_0, g_yyyyyzz_0_xxzzzzz_0, g_yyyyyzz_0_xyyyyyy_0, g_yyyyyzz_0_xyyyyyz_0, g_yyyyyzz_0_xyyyyzz_0, g_yyyyyzz_0_xyyyzzz_0, g_yyyyyzz_0_xyyzzzz_0, g_yyyyyzz_0_xyzzzzz_0, g_yyyyyzz_0_xzzzzzz_0, g_yyyyyzz_0_yyyyyyy_0, g_yyyyyzz_0_yyyyyyz_0, g_yyyyyzz_0_yyyyyzz_0, g_yyyyyzz_0_yyyyzzz_0, g_yyyyyzz_0_yyyzzzz_0, g_yyyyyzz_0_yyzzzzz_0, g_yyyyyzz_0_yzzzzzz_0, g_yyyyyzz_0_zzzzzzz_0, g_yyyyzz_0_xxxxxxx_1, g_yyyyzz_0_xxxxxxz_1, g_yyyyzz_0_xxxxxyz_1, g_yyyyzz_0_xxxxxz_1, g_yyyyzz_0_xxxxxzz_1, g_yyyyzz_0_xxxxyyz_1, g_yyyyzz_0_xxxxyz_1, g_yyyyzz_0_xxxxyzz_1, g_yyyyzz_0_xxxxzz_1, g_yyyyzz_0_xxxxzzz_1, g_yyyyzz_0_xxxyyyz_1, g_yyyyzz_0_xxxyyz_1, g_yyyyzz_0_xxxyyzz_1, g_yyyyzz_0_xxxyzz_1, g_yyyyzz_0_xxxyzzz_1, g_yyyyzz_0_xxxzzz_1, g_yyyyzz_0_xxxzzzz_1, g_yyyyzz_0_xxyyyyz_1, g_yyyyzz_0_xxyyyz_1, g_yyyyzz_0_xxyyyzz_1, g_yyyyzz_0_xxyyzz_1, g_yyyyzz_0_xxyyzzz_1, g_yyyyzz_0_xxyzzz_1, g_yyyyzz_0_xxyzzzz_1, g_yyyyzz_0_xxzzzz_1, g_yyyyzz_0_xxzzzzz_1, g_yyyyzz_0_xyyyyyz_1, g_yyyyzz_0_xyyyyz_1, g_yyyyzz_0_xyyyyzz_1, g_yyyyzz_0_xyyyzz_1, g_yyyyzz_0_xyyyzzz_1, g_yyyyzz_0_xyyzzz_1, g_yyyyzz_0_xyyzzzz_1, g_yyyyzz_0_xyzzzz_1, g_yyyyzz_0_xyzzzzz_1, g_yyyyzz_0_xzzzzz_1, g_yyyyzz_0_xzzzzzz_1, g_yyyyzz_0_yyyyyyz_1, g_yyyyzz_0_yyyyyz_1, g_yyyyzz_0_yyyyyzz_1, g_yyyyzz_0_yyyyzz_1, g_yyyyzz_0_yyyyzzz_1, g_yyyyzz_0_yyyzzz_1, g_yyyyzz_0_yyyzzzz_1, g_yyyyzz_0_yyzzzz_1, g_yyyyzz_0_yyzzzzz_1, g_yyyyzz_0_yzzzzz_1, g_yyyyzz_0_yzzzzzz_1, g_yyyyzz_0_zzzzzz_1, g_yyyyzz_0_zzzzzzz_1, g_yyyzz_0_xxxxxxx_0, g_yyyzz_0_xxxxxxx_1, g_yyyzz_0_xxxxxxz_0, g_yyyzz_0_xxxxxxz_1, g_yyyzz_0_xxxxxyz_0, g_yyyzz_0_xxxxxyz_1, g_yyyzz_0_xxxxxzz_0, g_yyyzz_0_xxxxxzz_1, g_yyyzz_0_xxxxyyz_0, g_yyyzz_0_xxxxyyz_1, g_yyyzz_0_xxxxyzz_0, g_yyyzz_0_xxxxyzz_1, g_yyyzz_0_xxxxzzz_0, g_yyyzz_0_xxxxzzz_1, g_yyyzz_0_xxxyyyz_0, g_yyyzz_0_xxxyyyz_1, g_yyyzz_0_xxxyyzz_0, g_yyyzz_0_xxxyyzz_1, g_yyyzz_0_xxxyzzz_0, g_yyyzz_0_xxxyzzz_1, g_yyyzz_0_xxxzzzz_0, g_yyyzz_0_xxxzzzz_1, g_yyyzz_0_xxyyyyz_0, g_yyyzz_0_xxyyyyz_1, g_yyyzz_0_xxyyyzz_0, g_yyyzz_0_xxyyyzz_1, g_yyyzz_0_xxyyzzz_0, g_yyyzz_0_xxyyzzz_1, g_yyyzz_0_xxyzzzz_0, g_yyyzz_0_xxyzzzz_1, g_yyyzz_0_xxzzzzz_0, g_yyyzz_0_xxzzzzz_1, g_yyyzz_0_xyyyyyz_0, g_yyyzz_0_xyyyyyz_1, g_yyyzz_0_xyyyyzz_0, g_yyyzz_0_xyyyyzz_1, g_yyyzz_0_xyyyzzz_0, g_yyyzz_0_xyyyzzz_1, g_yyyzz_0_xyyzzzz_0, g_yyyzz_0_xyyzzzz_1, g_yyyzz_0_xyzzzzz_0, g_yyyzz_0_xyzzzzz_1, g_yyyzz_0_xzzzzzz_0, g_yyyzz_0_xzzzzzz_1, g_yyyzz_0_yyyyyyz_0, g_yyyzz_0_yyyyyyz_1, g_yyyzz_0_yyyyyzz_0, g_yyyzz_0_yyyyyzz_1, g_yyyzz_0_yyyyzzz_0, g_yyyzz_0_yyyyzzz_1, g_yyyzz_0_yyyzzzz_0, g_yyyzz_0_yyyzzzz_1, g_yyyzz_0_yyzzzzz_0, g_yyyzz_0_yyzzzzz_1, g_yyyzz_0_yzzzzzz_0, g_yyyzz_0_yzzzzzz_1, g_yyyzz_0_zzzzzzz_0, g_yyyzz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyzz_0_xxxxxxx_0[i] = 4.0 * g_yyyzz_0_xxxxxxx_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxxx_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxxxy_0[i] = g_yyyyy_0_xxxxxxy_0[i] * fbe_0 - g_yyyyy_0_xxxxxxy_1[i] * fz_be_0 + g_yyyyyz_0_xxxxxxy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxxxxz_0[i] = 4.0 * g_yyyzz_0_xxxxxxz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxxz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxxyy_0[i] = g_yyyyy_0_xxxxxyy_0[i] * fbe_0 - g_yyyyy_0_xxxxxyy_1[i] * fz_be_0 + g_yyyyyz_0_xxxxxyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxxxyz_0[i] = 4.0 * g_yyyzz_0_xxxxxyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxyz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxxzz_0[i] = 4.0 * g_yyyzz_0_xxxxxzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxxzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxyyy_0[i] = g_yyyyy_0_xxxxyyy_0[i] * fbe_0 - g_yyyyy_0_xxxxyyy_1[i] * fz_be_0 + g_yyyyyz_0_xxxxyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxxyyz_0[i] = 4.0 * g_yyyzz_0_xxxxyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xxxxyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxyzz_0[i] = 4.0 * g_yyyzz_0_xxxxyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxyzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxxzzz_0[i] = 4.0 * g_yyyzz_0_xxxxzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxxzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxyyyy_0[i] = g_yyyyy_0_xxxyyyy_0[i] * fbe_0 - g_yyyyy_0_xxxyyyy_1[i] * fz_be_0 + g_yyyyyz_0_xxxyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxxyyyz_0[i] = 4.0 * g_yyyzz_0_xxxyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_xxxyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxyyzz_0[i] = 4.0 * g_yyyzz_0_xxxyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xxxyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxyzzz_0[i] = 4.0 * g_yyyzz_0_xxxyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxyzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxxzzzz_0[i] = 4.0 * g_yyyzz_0_xxxzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxxzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyyyyy_0[i] = g_yyyyy_0_xxyyyyy_0[i] * fbe_0 - g_yyyyy_0_xxyyyyy_1[i] * fz_be_0 + g_yyyyyz_0_xxyyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xxyyyyz_0[i] = 4.0 * g_yyyzz_0_xxyyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_0_xxyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyyyzz_0[i] = 4.0 * g_yyyzz_0_xxyyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_xxyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyyzzz_0[i] = 4.0 * g_yyyzz_0_xxyyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xxyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxyzzzz_0[i] = 4.0 * g_yyyzz_0_xxyzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxyzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xxzzzzz_0[i] = 4.0 * g_yyyzz_0_xxzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xxzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyyyyy_0[i] = g_yyyyy_0_xyyyyyy_0[i] * fbe_0 - g_yyyyy_0_xyyyyyy_1[i] * fz_be_0 + g_yyyyyz_0_xyyyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_xyyyyyz_0[i] = 4.0 * g_yyyzz_0_xyyyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyyzz_0_xyyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyyyzz_0[i] = 4.0 * g_yyyzz_0_xyyyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_0_xyyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyyzzz_0[i] = 4.0 * g_yyyzz_0_xyyyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_xyyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyyzzzz_0[i] = 4.0 * g_yyyzz_0_xyyzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_xyzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xyzzzzz_0[i] = 4.0 * g_yyyzz_0_xyzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xyzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_xzzzzzz_0[i] = 4.0 * g_yyyzz_0_xzzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_xzzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyyyyy_0[i] = g_yyyyy_0_yyyyyyy_0[i] * fbe_0 - g_yyyyy_0_yyyyyyy_1[i] * fz_be_0 + g_yyyyyz_0_yyyyyyy_1[i] * wa_z[i];

        g_yyyyyzz_0_yyyyyyz_0[i] = 4.0 * g_yyyzz_0_yyyyyyz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyyyzz_0_yyyyyz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyyyzz_0[i] = 4.0 * g_yyyzz_0_yyyyyzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyyyzz_0_yyyyzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyyzzz_0[i] = 4.0 * g_yyyzz_0_yyyyzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyyyzz_0_yyyzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyyzzzz_0[i] = 4.0 * g_yyyzz_0_yyyzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyyyzz_0_yyzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yyzzzzz_0[i] = 4.0 * g_yyyzz_0_yyzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_0_yzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_yzzzzzz_0[i] = 4.0 * g_yyyzz_0_yzzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_yzzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_zzzzzz_1[i] * fi_acd_0 + g_yyyyzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yyyyyzz_0_zzzzzzz_0[i] = 4.0 * g_yyyzz_0_zzzzzzz_0[i] * fbe_0 - 4.0 * g_yyyzz_0_zzzzzzz_1[i] * fz_be_0 + g_yyyyzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1116-1152 components of targeted buffer : KSK

    auto g_yyyyzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 1116);

    auto g_yyyyzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 1117);

    auto g_yyyyzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 1118);

    auto g_yyyyzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 1119);

    auto g_yyyyzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 1120);

    auto g_yyyyzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 1121);

    auto g_yyyyzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 1122);

    auto g_yyyyzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 1123);

    auto g_yyyyzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 1124);

    auto g_yyyyzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 1125);

    auto g_yyyyzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1126);

    auto g_yyyyzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1127);

    auto g_yyyyzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1128);

    auto g_yyyyzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1129);

    auto g_yyyyzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1130);

    auto g_yyyyzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1131);

    auto g_yyyyzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1132);

    auto g_yyyyzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1133);

    auto g_yyyyzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1134);

    auto g_yyyyzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1135);

    auto g_yyyyzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1136);

    auto g_yyyyzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1137);

    auto g_yyyyzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1138);

    auto g_yyyyzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1139);

    auto g_yyyyzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1140);

    auto g_yyyyzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1141);

    auto g_yyyyzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1142);

    auto g_yyyyzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1143);

    auto g_yyyyzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1144);

    auto g_yyyyzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1145);

    auto g_yyyyzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1146);

    auto g_yyyyzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1147);

    auto g_yyyyzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1148);

    auto g_yyyyzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1149);

    auto g_yyyyzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1150);

    auto g_yyyyzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1151);

    #pragma omp simd aligned(g_yyyyz_0_xxxxxxy_0, g_yyyyz_0_xxxxxxy_1, g_yyyyz_0_xxxxxyy_0, g_yyyyz_0_xxxxxyy_1, g_yyyyz_0_xxxxyyy_0, g_yyyyz_0_xxxxyyy_1, g_yyyyz_0_xxxyyyy_0, g_yyyyz_0_xxxyyyy_1, g_yyyyz_0_xxyyyyy_0, g_yyyyz_0_xxyyyyy_1, g_yyyyz_0_xyyyyyy_0, g_yyyyz_0_xyyyyyy_1, g_yyyyz_0_yyyyyyy_0, g_yyyyz_0_yyyyyyy_1, g_yyyyzz_0_xxxxxxy_1, g_yyyyzz_0_xxxxxyy_1, g_yyyyzz_0_xxxxyyy_1, g_yyyyzz_0_xxxyyyy_1, g_yyyyzz_0_xxyyyyy_1, g_yyyyzz_0_xyyyyyy_1, g_yyyyzz_0_yyyyyyy_1, g_yyyyzzz_0_xxxxxxx_0, g_yyyyzzz_0_xxxxxxy_0, g_yyyyzzz_0_xxxxxxz_0, g_yyyyzzz_0_xxxxxyy_0, g_yyyyzzz_0_xxxxxyz_0, g_yyyyzzz_0_xxxxxzz_0, g_yyyyzzz_0_xxxxyyy_0, g_yyyyzzz_0_xxxxyyz_0, g_yyyyzzz_0_xxxxyzz_0, g_yyyyzzz_0_xxxxzzz_0, g_yyyyzzz_0_xxxyyyy_0, g_yyyyzzz_0_xxxyyyz_0, g_yyyyzzz_0_xxxyyzz_0, g_yyyyzzz_0_xxxyzzz_0, g_yyyyzzz_0_xxxzzzz_0, g_yyyyzzz_0_xxyyyyy_0, g_yyyyzzz_0_xxyyyyz_0, g_yyyyzzz_0_xxyyyzz_0, g_yyyyzzz_0_xxyyzzz_0, g_yyyyzzz_0_xxyzzzz_0, g_yyyyzzz_0_xxzzzzz_0, g_yyyyzzz_0_xyyyyyy_0, g_yyyyzzz_0_xyyyyyz_0, g_yyyyzzz_0_xyyyyzz_0, g_yyyyzzz_0_xyyyzzz_0, g_yyyyzzz_0_xyyzzzz_0, g_yyyyzzz_0_xyzzzzz_0, g_yyyyzzz_0_xzzzzzz_0, g_yyyyzzz_0_yyyyyyy_0, g_yyyyzzz_0_yyyyyyz_0, g_yyyyzzz_0_yyyyyzz_0, g_yyyyzzz_0_yyyyzzz_0, g_yyyyzzz_0_yyyzzzz_0, g_yyyyzzz_0_yyzzzzz_0, g_yyyyzzz_0_yzzzzzz_0, g_yyyyzzz_0_zzzzzzz_0, g_yyyzzz_0_xxxxxxx_1, g_yyyzzz_0_xxxxxxz_1, g_yyyzzz_0_xxxxxyz_1, g_yyyzzz_0_xxxxxz_1, g_yyyzzz_0_xxxxxzz_1, g_yyyzzz_0_xxxxyyz_1, g_yyyzzz_0_xxxxyz_1, g_yyyzzz_0_xxxxyzz_1, g_yyyzzz_0_xxxxzz_1, g_yyyzzz_0_xxxxzzz_1, g_yyyzzz_0_xxxyyyz_1, g_yyyzzz_0_xxxyyz_1, g_yyyzzz_0_xxxyyzz_1, g_yyyzzz_0_xxxyzz_1, g_yyyzzz_0_xxxyzzz_1, g_yyyzzz_0_xxxzzz_1, g_yyyzzz_0_xxxzzzz_1, g_yyyzzz_0_xxyyyyz_1, g_yyyzzz_0_xxyyyz_1, g_yyyzzz_0_xxyyyzz_1, g_yyyzzz_0_xxyyzz_1, g_yyyzzz_0_xxyyzzz_1, g_yyyzzz_0_xxyzzz_1, g_yyyzzz_0_xxyzzzz_1, g_yyyzzz_0_xxzzzz_1, g_yyyzzz_0_xxzzzzz_1, g_yyyzzz_0_xyyyyyz_1, g_yyyzzz_0_xyyyyz_1, g_yyyzzz_0_xyyyyzz_1, g_yyyzzz_0_xyyyzz_1, g_yyyzzz_0_xyyyzzz_1, g_yyyzzz_0_xyyzzz_1, g_yyyzzz_0_xyyzzzz_1, g_yyyzzz_0_xyzzzz_1, g_yyyzzz_0_xyzzzzz_1, g_yyyzzz_0_xzzzzz_1, g_yyyzzz_0_xzzzzzz_1, g_yyyzzz_0_yyyyyyz_1, g_yyyzzz_0_yyyyyz_1, g_yyyzzz_0_yyyyyzz_1, g_yyyzzz_0_yyyyzz_1, g_yyyzzz_0_yyyyzzz_1, g_yyyzzz_0_yyyzzz_1, g_yyyzzz_0_yyyzzzz_1, g_yyyzzz_0_yyzzzz_1, g_yyyzzz_0_yyzzzzz_1, g_yyyzzz_0_yzzzzz_1, g_yyyzzz_0_yzzzzzz_1, g_yyyzzz_0_zzzzzz_1, g_yyyzzz_0_zzzzzzz_1, g_yyzzz_0_xxxxxxx_0, g_yyzzz_0_xxxxxxx_1, g_yyzzz_0_xxxxxxz_0, g_yyzzz_0_xxxxxxz_1, g_yyzzz_0_xxxxxyz_0, g_yyzzz_0_xxxxxyz_1, g_yyzzz_0_xxxxxzz_0, g_yyzzz_0_xxxxxzz_1, g_yyzzz_0_xxxxyyz_0, g_yyzzz_0_xxxxyyz_1, g_yyzzz_0_xxxxyzz_0, g_yyzzz_0_xxxxyzz_1, g_yyzzz_0_xxxxzzz_0, g_yyzzz_0_xxxxzzz_1, g_yyzzz_0_xxxyyyz_0, g_yyzzz_0_xxxyyyz_1, g_yyzzz_0_xxxyyzz_0, g_yyzzz_0_xxxyyzz_1, g_yyzzz_0_xxxyzzz_0, g_yyzzz_0_xxxyzzz_1, g_yyzzz_0_xxxzzzz_0, g_yyzzz_0_xxxzzzz_1, g_yyzzz_0_xxyyyyz_0, g_yyzzz_0_xxyyyyz_1, g_yyzzz_0_xxyyyzz_0, g_yyzzz_0_xxyyyzz_1, g_yyzzz_0_xxyyzzz_0, g_yyzzz_0_xxyyzzz_1, g_yyzzz_0_xxyzzzz_0, g_yyzzz_0_xxyzzzz_1, g_yyzzz_0_xxzzzzz_0, g_yyzzz_0_xxzzzzz_1, g_yyzzz_0_xyyyyyz_0, g_yyzzz_0_xyyyyyz_1, g_yyzzz_0_xyyyyzz_0, g_yyzzz_0_xyyyyzz_1, g_yyzzz_0_xyyyzzz_0, g_yyzzz_0_xyyyzzz_1, g_yyzzz_0_xyyzzzz_0, g_yyzzz_0_xyyzzzz_1, g_yyzzz_0_xyzzzzz_0, g_yyzzz_0_xyzzzzz_1, g_yyzzz_0_xzzzzzz_0, g_yyzzz_0_xzzzzzz_1, g_yyzzz_0_yyyyyyz_0, g_yyzzz_0_yyyyyyz_1, g_yyzzz_0_yyyyyzz_0, g_yyzzz_0_yyyyyzz_1, g_yyzzz_0_yyyyzzz_0, g_yyzzz_0_yyyyzzz_1, g_yyzzz_0_yyyzzzz_0, g_yyzzz_0_yyyzzzz_1, g_yyzzz_0_yyzzzzz_0, g_yyzzz_0_yyzzzzz_1, g_yyzzz_0_yzzzzzz_0, g_yyzzz_0_yzzzzzz_1, g_yyzzz_0_zzzzzzz_0, g_yyzzz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyzzz_0_xxxxxxx_0[i] = 3.0 * g_yyzzz_0_xxxxxxx_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxxx_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxxxy_0[i] = 2.0 * g_yyyyz_0_xxxxxxy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxxxxy_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxxy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxxxxz_0[i] = 3.0 * g_yyzzz_0_xxxxxxz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxxz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxxyy_0[i] = 2.0 * g_yyyyz_0_xxxxxyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxxxyy_1[i] * fz_be_0 + g_yyyyzz_0_xxxxxyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxxxyz_0[i] = 3.0 * g_yyzzz_0_xxxxxyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxyz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxxzz_0[i] = 3.0 * g_yyzzz_0_xxxxxzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxxzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxyyy_0[i] = 2.0 * g_yyyyz_0_xxxxyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxxyyy_1[i] * fz_be_0 + g_yyyyzz_0_xxxxyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxxyyz_0[i] = 3.0 * g_yyzzz_0_xxxxyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xxxxyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxyzz_0[i] = 3.0 * g_yyzzz_0_xxxxyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxyzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxxzzz_0[i] = 3.0 * g_yyzzz_0_xxxxzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxxzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxyyyy_0[i] = 2.0 * g_yyyyz_0_xxxyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxxyyyy_1[i] * fz_be_0 + g_yyyyzz_0_xxxyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxxyyyz_0[i] = 3.0 * g_yyzzz_0_xxxyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_xxxyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxyyzz_0[i] = 3.0 * g_yyzzz_0_xxxyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xxxyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxyzzz_0[i] = 3.0 * g_yyzzz_0_xxxyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxyzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxxzzzz_0[i] = 3.0 * g_yyzzz_0_xxxzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxxzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyyyyy_0[i] = 2.0 * g_yyyyz_0_xxyyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xxyyyyy_1[i] * fz_be_0 + g_yyyyzz_0_xxyyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xxyyyyz_0[i] = 3.0 * g_yyzzz_0_xxyyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_0_xxyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyyyzz_0[i] = 3.0 * g_yyzzz_0_xxyyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_xxyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyyzzz_0[i] = 3.0 * g_yyzzz_0_xxyyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xxyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxyzzzz_0[i] = 3.0 * g_yyzzz_0_xxyzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxyzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xxzzzzz_0[i] = 3.0 * g_yyzzz_0_xxzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xxzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyyyyy_0[i] = 2.0 * g_yyyyz_0_xyyyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_xyyyyyy_1[i] * fz_be_0 + g_yyyyzz_0_xyyyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_xyyyyyz_0[i] = 3.0 * g_yyzzz_0_xyyyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyzzz_0_xyyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyyyzz_0[i] = 3.0 * g_yyzzz_0_xyyyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_0_xyyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyyzzz_0[i] = 3.0 * g_yyzzz_0_xyyyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_xyyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyyzzzz_0[i] = 3.0 * g_yyzzz_0_xyyzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_xyzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xyzzzzz_0[i] = 3.0 * g_yyzzz_0_xyzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xyzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_xzzzzzz_0[i] = 3.0 * g_yyzzz_0_xzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_xzzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyyyyy_0[i] = 2.0 * g_yyyyz_0_yyyyyyy_0[i] * fbe_0 - 2.0 * g_yyyyz_0_yyyyyyy_1[i] * fz_be_0 + g_yyyyzz_0_yyyyyyy_1[i] * wa_z[i];

        g_yyyyzzz_0_yyyyyyz_0[i] = 3.0 * g_yyzzz_0_yyyyyyz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyyzzz_0_yyyyyz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyyyzz_0[i] = 3.0 * g_yyzzz_0_yyyyyzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyyzzz_0_yyyyzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyyzzz_0[i] = 3.0 * g_yyzzz_0_yyyyzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyyzzz_0_yyyzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyyzzzz_0[i] = 3.0 * g_yyzzz_0_yyyzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyyzzz_0_yyzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yyzzzzz_0[i] = 3.0 * g_yyzzz_0_yyzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_0_yzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_yzzzzzz_0[i] = 3.0 * g_yyzzz_0_yzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_yzzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_zzzzzz_1[i] * fi_acd_0 + g_yyyzzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yyyyzzz_0_zzzzzzz_0[i] = 3.0 * g_yyzzz_0_zzzzzzz_0[i] * fbe_0 - 3.0 * g_yyzzz_0_zzzzzzz_1[i] * fz_be_0 + g_yyyzzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1152-1188 components of targeted buffer : KSK

    auto g_yyyzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 1152);

    auto g_yyyzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 1153);

    auto g_yyyzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 1154);

    auto g_yyyzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 1155);

    auto g_yyyzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 1156);

    auto g_yyyzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 1157);

    auto g_yyyzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 1158);

    auto g_yyyzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 1159);

    auto g_yyyzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 1160);

    auto g_yyyzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 1161);

    auto g_yyyzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1162);

    auto g_yyyzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1163);

    auto g_yyyzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1164);

    auto g_yyyzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1165);

    auto g_yyyzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1166);

    auto g_yyyzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1167);

    auto g_yyyzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1168);

    auto g_yyyzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1169);

    auto g_yyyzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1170);

    auto g_yyyzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1171);

    auto g_yyyzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1172);

    auto g_yyyzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1173);

    auto g_yyyzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1174);

    auto g_yyyzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1175);

    auto g_yyyzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1176);

    auto g_yyyzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1177);

    auto g_yyyzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1178);

    auto g_yyyzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1179);

    auto g_yyyzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1180);

    auto g_yyyzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1181);

    auto g_yyyzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1182);

    auto g_yyyzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1183);

    auto g_yyyzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1184);

    auto g_yyyzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1185);

    auto g_yyyzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1186);

    auto g_yyyzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1187);

    #pragma omp simd aligned(g_yyyzz_0_xxxxxxy_0, g_yyyzz_0_xxxxxxy_1, g_yyyzz_0_xxxxxyy_0, g_yyyzz_0_xxxxxyy_1, g_yyyzz_0_xxxxyyy_0, g_yyyzz_0_xxxxyyy_1, g_yyyzz_0_xxxyyyy_0, g_yyyzz_0_xxxyyyy_1, g_yyyzz_0_xxyyyyy_0, g_yyyzz_0_xxyyyyy_1, g_yyyzz_0_xyyyyyy_0, g_yyyzz_0_xyyyyyy_1, g_yyyzz_0_yyyyyyy_0, g_yyyzz_0_yyyyyyy_1, g_yyyzzz_0_xxxxxxy_1, g_yyyzzz_0_xxxxxyy_1, g_yyyzzz_0_xxxxyyy_1, g_yyyzzz_0_xxxyyyy_1, g_yyyzzz_0_xxyyyyy_1, g_yyyzzz_0_xyyyyyy_1, g_yyyzzz_0_yyyyyyy_1, g_yyyzzzz_0_xxxxxxx_0, g_yyyzzzz_0_xxxxxxy_0, g_yyyzzzz_0_xxxxxxz_0, g_yyyzzzz_0_xxxxxyy_0, g_yyyzzzz_0_xxxxxyz_0, g_yyyzzzz_0_xxxxxzz_0, g_yyyzzzz_0_xxxxyyy_0, g_yyyzzzz_0_xxxxyyz_0, g_yyyzzzz_0_xxxxyzz_0, g_yyyzzzz_0_xxxxzzz_0, g_yyyzzzz_0_xxxyyyy_0, g_yyyzzzz_0_xxxyyyz_0, g_yyyzzzz_0_xxxyyzz_0, g_yyyzzzz_0_xxxyzzz_0, g_yyyzzzz_0_xxxzzzz_0, g_yyyzzzz_0_xxyyyyy_0, g_yyyzzzz_0_xxyyyyz_0, g_yyyzzzz_0_xxyyyzz_0, g_yyyzzzz_0_xxyyzzz_0, g_yyyzzzz_0_xxyzzzz_0, g_yyyzzzz_0_xxzzzzz_0, g_yyyzzzz_0_xyyyyyy_0, g_yyyzzzz_0_xyyyyyz_0, g_yyyzzzz_0_xyyyyzz_0, g_yyyzzzz_0_xyyyzzz_0, g_yyyzzzz_0_xyyzzzz_0, g_yyyzzzz_0_xyzzzzz_0, g_yyyzzzz_0_xzzzzzz_0, g_yyyzzzz_0_yyyyyyy_0, g_yyyzzzz_0_yyyyyyz_0, g_yyyzzzz_0_yyyyyzz_0, g_yyyzzzz_0_yyyyzzz_0, g_yyyzzzz_0_yyyzzzz_0, g_yyyzzzz_0_yyzzzzz_0, g_yyyzzzz_0_yzzzzzz_0, g_yyyzzzz_0_zzzzzzz_0, g_yyzzzz_0_xxxxxxx_1, g_yyzzzz_0_xxxxxxz_1, g_yyzzzz_0_xxxxxyz_1, g_yyzzzz_0_xxxxxz_1, g_yyzzzz_0_xxxxxzz_1, g_yyzzzz_0_xxxxyyz_1, g_yyzzzz_0_xxxxyz_1, g_yyzzzz_0_xxxxyzz_1, g_yyzzzz_0_xxxxzz_1, g_yyzzzz_0_xxxxzzz_1, g_yyzzzz_0_xxxyyyz_1, g_yyzzzz_0_xxxyyz_1, g_yyzzzz_0_xxxyyzz_1, g_yyzzzz_0_xxxyzz_1, g_yyzzzz_0_xxxyzzz_1, g_yyzzzz_0_xxxzzz_1, g_yyzzzz_0_xxxzzzz_1, g_yyzzzz_0_xxyyyyz_1, g_yyzzzz_0_xxyyyz_1, g_yyzzzz_0_xxyyyzz_1, g_yyzzzz_0_xxyyzz_1, g_yyzzzz_0_xxyyzzz_1, g_yyzzzz_0_xxyzzz_1, g_yyzzzz_0_xxyzzzz_1, g_yyzzzz_0_xxzzzz_1, g_yyzzzz_0_xxzzzzz_1, g_yyzzzz_0_xyyyyyz_1, g_yyzzzz_0_xyyyyz_1, g_yyzzzz_0_xyyyyzz_1, g_yyzzzz_0_xyyyzz_1, g_yyzzzz_0_xyyyzzz_1, g_yyzzzz_0_xyyzzz_1, g_yyzzzz_0_xyyzzzz_1, g_yyzzzz_0_xyzzzz_1, g_yyzzzz_0_xyzzzzz_1, g_yyzzzz_0_xzzzzz_1, g_yyzzzz_0_xzzzzzz_1, g_yyzzzz_0_yyyyyyz_1, g_yyzzzz_0_yyyyyz_1, g_yyzzzz_0_yyyyyzz_1, g_yyzzzz_0_yyyyzz_1, g_yyzzzz_0_yyyyzzz_1, g_yyzzzz_0_yyyzzz_1, g_yyzzzz_0_yyyzzzz_1, g_yyzzzz_0_yyzzzz_1, g_yyzzzz_0_yyzzzzz_1, g_yyzzzz_0_yzzzzz_1, g_yyzzzz_0_yzzzzzz_1, g_yyzzzz_0_zzzzzz_1, g_yyzzzz_0_zzzzzzz_1, g_yzzzz_0_xxxxxxx_0, g_yzzzz_0_xxxxxxx_1, g_yzzzz_0_xxxxxxz_0, g_yzzzz_0_xxxxxxz_1, g_yzzzz_0_xxxxxyz_0, g_yzzzz_0_xxxxxyz_1, g_yzzzz_0_xxxxxzz_0, g_yzzzz_0_xxxxxzz_1, g_yzzzz_0_xxxxyyz_0, g_yzzzz_0_xxxxyyz_1, g_yzzzz_0_xxxxyzz_0, g_yzzzz_0_xxxxyzz_1, g_yzzzz_0_xxxxzzz_0, g_yzzzz_0_xxxxzzz_1, g_yzzzz_0_xxxyyyz_0, g_yzzzz_0_xxxyyyz_1, g_yzzzz_0_xxxyyzz_0, g_yzzzz_0_xxxyyzz_1, g_yzzzz_0_xxxyzzz_0, g_yzzzz_0_xxxyzzz_1, g_yzzzz_0_xxxzzzz_0, g_yzzzz_0_xxxzzzz_1, g_yzzzz_0_xxyyyyz_0, g_yzzzz_0_xxyyyyz_1, g_yzzzz_0_xxyyyzz_0, g_yzzzz_0_xxyyyzz_1, g_yzzzz_0_xxyyzzz_0, g_yzzzz_0_xxyyzzz_1, g_yzzzz_0_xxyzzzz_0, g_yzzzz_0_xxyzzzz_1, g_yzzzz_0_xxzzzzz_0, g_yzzzz_0_xxzzzzz_1, g_yzzzz_0_xyyyyyz_0, g_yzzzz_0_xyyyyyz_1, g_yzzzz_0_xyyyyzz_0, g_yzzzz_0_xyyyyzz_1, g_yzzzz_0_xyyyzzz_0, g_yzzzz_0_xyyyzzz_1, g_yzzzz_0_xyyzzzz_0, g_yzzzz_0_xyyzzzz_1, g_yzzzz_0_xyzzzzz_0, g_yzzzz_0_xyzzzzz_1, g_yzzzz_0_xzzzzzz_0, g_yzzzz_0_xzzzzzz_1, g_yzzzz_0_yyyyyyz_0, g_yzzzz_0_yyyyyyz_1, g_yzzzz_0_yyyyyzz_0, g_yzzzz_0_yyyyyzz_1, g_yzzzz_0_yyyyzzz_0, g_yzzzz_0_yyyyzzz_1, g_yzzzz_0_yyyzzzz_0, g_yzzzz_0_yyyzzzz_1, g_yzzzz_0_yyzzzzz_0, g_yzzzz_0_yyzzzzz_1, g_yzzzz_0_yzzzzzz_0, g_yzzzz_0_yzzzzzz_1, g_yzzzz_0_zzzzzzz_0, g_yzzzz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzzzz_0_xxxxxxx_0[i] = 2.0 * g_yzzzz_0_xxxxxxx_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxxx_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxxxy_0[i] = 3.0 * g_yyyzz_0_xxxxxxy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxxxxy_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxxy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxxxxz_0[i] = 2.0 * g_yzzzz_0_xxxxxxz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxxz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxxyy_0[i] = 3.0 * g_yyyzz_0_xxxxxyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxxxyy_1[i] * fz_be_0 + g_yyyzzz_0_xxxxxyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxxxyz_0[i] = 2.0 * g_yzzzz_0_xxxxxyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxyz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxxzz_0[i] = 2.0 * g_yzzzz_0_xxxxxzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxxzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxyyy_0[i] = 3.0 * g_yyyzz_0_xxxxyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxxyyy_1[i] * fz_be_0 + g_yyyzzz_0_xxxxyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxxyyz_0[i] = 2.0 * g_yzzzz_0_xxxxyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxyzz_0[i] = 2.0 * g_yzzzz_0_xxxxyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxyzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxxzzz_0[i] = 2.0 * g_yzzzz_0_xxxxzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxxzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxyyyy_0[i] = 3.0 * g_yyyzz_0_xxxyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxxyyyy_1[i] * fz_be_0 + g_yyyzzz_0_xxxyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxxyyyz_0[i] = 2.0 * g_yzzzz_0_xxxyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxyyzz_0[i] = 2.0 * g_yzzzz_0_xxxyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxyzzz_0[i] = 2.0 * g_yzzzz_0_xxxyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxyzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxxzzzz_0[i] = 2.0 * g_yzzzz_0_xxxzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxxzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyyyyy_0[i] = 3.0 * g_yyyzz_0_xxyyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xxyyyyy_1[i] * fz_be_0 + g_yyyzzz_0_xxyyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xxyyyyz_0[i] = 2.0 * g_yzzzz_0_xxyyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyyyzz_0[i] = 2.0 * g_yzzzz_0_xxyyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyyzzz_0[i] = 2.0 * g_yzzzz_0_xxyyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxyzzzz_0[i] = 2.0 * g_yzzzz_0_xxyzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxyzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xxzzzzz_0[i] = 2.0 * g_yzzzz_0_xxzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xxzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyyyyy_0[i] = 3.0 * g_yyyzz_0_xyyyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_xyyyyyy_1[i] * fz_be_0 + g_yyyzzz_0_xyyyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_xyyyyyz_0[i] = 2.0 * g_yzzzz_0_xyyyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyyyzz_0[i] = 2.0 * g_yzzzz_0_xyyyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyyzzz_0[i] = 2.0 * g_yzzzz_0_xyyyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyyzzzz_0[i] = 2.0 * g_yzzzz_0_xyyzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xyzzzzz_0[i] = 2.0 * g_yzzzz_0_xyzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xyzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_xzzzzzz_0[i] = 2.0 * g_yzzzz_0_xzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_xzzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyyyyy_0[i] = 3.0 * g_yyyzz_0_yyyyyyy_0[i] * fbe_0 - 3.0 * g_yyyzz_0_yyyyyyy_1[i] * fz_be_0 + g_yyyzzz_0_yyyyyyy_1[i] * wa_z[i];

        g_yyyzzzz_0_yyyyyyz_0[i] = 2.0 * g_yzzzz_0_yyyyyyz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyyyzz_0[i] = 2.0 * g_yzzzz_0_yyyyyzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyyzzz_0[i] = 2.0 * g_yzzzz_0_yyyyzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyyzzzz_0[i] = 2.0 * g_yzzzz_0_yyyzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yyzzzzz_0[i] = 2.0 * g_yzzzz_0_yyzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_yzzzzzz_0[i] = 2.0 * g_yzzzz_0_yzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_yzzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_yyzzzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yyyzzzz_0_zzzzzzz_0[i] = 2.0 * g_yzzzz_0_zzzzzzz_0[i] * fbe_0 - 2.0 * g_yzzzz_0_zzzzzzz_1[i] * fz_be_0 + g_yyzzzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1188-1224 components of targeted buffer : KSK

    auto g_yyzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 1188);

    auto g_yyzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 1189);

    auto g_yyzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 1190);

    auto g_yyzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 1191);

    auto g_yyzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 1192);

    auto g_yyzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 1193);

    auto g_yyzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 1194);

    auto g_yyzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 1195);

    auto g_yyzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 1196);

    auto g_yyzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 1197);

    auto g_yyzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1198);

    auto g_yyzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1199);

    auto g_yyzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1200);

    auto g_yyzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1201);

    auto g_yyzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1202);

    auto g_yyzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1203);

    auto g_yyzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1204);

    auto g_yyzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1205);

    auto g_yyzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1206);

    auto g_yyzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1207);

    auto g_yyzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1208);

    auto g_yyzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1209);

    auto g_yyzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1210);

    auto g_yyzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1211);

    auto g_yyzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1212);

    auto g_yyzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1213);

    auto g_yyzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1214);

    auto g_yyzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1215);

    auto g_yyzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1216);

    auto g_yyzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1217);

    auto g_yyzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1218);

    auto g_yyzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1219);

    auto g_yyzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1220);

    auto g_yyzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1221);

    auto g_yyzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1222);

    auto g_yyzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1223);

    #pragma omp simd aligned(g_yyzzz_0_xxxxxxy_0, g_yyzzz_0_xxxxxxy_1, g_yyzzz_0_xxxxxyy_0, g_yyzzz_0_xxxxxyy_1, g_yyzzz_0_xxxxyyy_0, g_yyzzz_0_xxxxyyy_1, g_yyzzz_0_xxxyyyy_0, g_yyzzz_0_xxxyyyy_1, g_yyzzz_0_xxyyyyy_0, g_yyzzz_0_xxyyyyy_1, g_yyzzz_0_xyyyyyy_0, g_yyzzz_0_xyyyyyy_1, g_yyzzz_0_yyyyyyy_0, g_yyzzz_0_yyyyyyy_1, g_yyzzzz_0_xxxxxxy_1, g_yyzzzz_0_xxxxxyy_1, g_yyzzzz_0_xxxxyyy_1, g_yyzzzz_0_xxxyyyy_1, g_yyzzzz_0_xxyyyyy_1, g_yyzzzz_0_xyyyyyy_1, g_yyzzzz_0_yyyyyyy_1, g_yyzzzzz_0_xxxxxxx_0, g_yyzzzzz_0_xxxxxxy_0, g_yyzzzzz_0_xxxxxxz_0, g_yyzzzzz_0_xxxxxyy_0, g_yyzzzzz_0_xxxxxyz_0, g_yyzzzzz_0_xxxxxzz_0, g_yyzzzzz_0_xxxxyyy_0, g_yyzzzzz_0_xxxxyyz_0, g_yyzzzzz_0_xxxxyzz_0, g_yyzzzzz_0_xxxxzzz_0, g_yyzzzzz_0_xxxyyyy_0, g_yyzzzzz_0_xxxyyyz_0, g_yyzzzzz_0_xxxyyzz_0, g_yyzzzzz_0_xxxyzzz_0, g_yyzzzzz_0_xxxzzzz_0, g_yyzzzzz_0_xxyyyyy_0, g_yyzzzzz_0_xxyyyyz_0, g_yyzzzzz_0_xxyyyzz_0, g_yyzzzzz_0_xxyyzzz_0, g_yyzzzzz_0_xxyzzzz_0, g_yyzzzzz_0_xxzzzzz_0, g_yyzzzzz_0_xyyyyyy_0, g_yyzzzzz_0_xyyyyyz_0, g_yyzzzzz_0_xyyyyzz_0, g_yyzzzzz_0_xyyyzzz_0, g_yyzzzzz_0_xyyzzzz_0, g_yyzzzzz_0_xyzzzzz_0, g_yyzzzzz_0_xzzzzzz_0, g_yyzzzzz_0_yyyyyyy_0, g_yyzzzzz_0_yyyyyyz_0, g_yyzzzzz_0_yyyyyzz_0, g_yyzzzzz_0_yyyyzzz_0, g_yyzzzzz_0_yyyzzzz_0, g_yyzzzzz_0_yyzzzzz_0, g_yyzzzzz_0_yzzzzzz_0, g_yyzzzzz_0_zzzzzzz_0, g_yzzzzz_0_xxxxxxx_1, g_yzzzzz_0_xxxxxxz_1, g_yzzzzz_0_xxxxxyz_1, g_yzzzzz_0_xxxxxz_1, g_yzzzzz_0_xxxxxzz_1, g_yzzzzz_0_xxxxyyz_1, g_yzzzzz_0_xxxxyz_1, g_yzzzzz_0_xxxxyzz_1, g_yzzzzz_0_xxxxzz_1, g_yzzzzz_0_xxxxzzz_1, g_yzzzzz_0_xxxyyyz_1, g_yzzzzz_0_xxxyyz_1, g_yzzzzz_0_xxxyyzz_1, g_yzzzzz_0_xxxyzz_1, g_yzzzzz_0_xxxyzzz_1, g_yzzzzz_0_xxxzzz_1, g_yzzzzz_0_xxxzzzz_1, g_yzzzzz_0_xxyyyyz_1, g_yzzzzz_0_xxyyyz_1, g_yzzzzz_0_xxyyyzz_1, g_yzzzzz_0_xxyyzz_1, g_yzzzzz_0_xxyyzzz_1, g_yzzzzz_0_xxyzzz_1, g_yzzzzz_0_xxyzzzz_1, g_yzzzzz_0_xxzzzz_1, g_yzzzzz_0_xxzzzzz_1, g_yzzzzz_0_xyyyyyz_1, g_yzzzzz_0_xyyyyz_1, g_yzzzzz_0_xyyyyzz_1, g_yzzzzz_0_xyyyzz_1, g_yzzzzz_0_xyyyzzz_1, g_yzzzzz_0_xyyzzz_1, g_yzzzzz_0_xyyzzzz_1, g_yzzzzz_0_xyzzzz_1, g_yzzzzz_0_xyzzzzz_1, g_yzzzzz_0_xzzzzz_1, g_yzzzzz_0_xzzzzzz_1, g_yzzzzz_0_yyyyyyz_1, g_yzzzzz_0_yyyyyz_1, g_yzzzzz_0_yyyyyzz_1, g_yzzzzz_0_yyyyzz_1, g_yzzzzz_0_yyyyzzz_1, g_yzzzzz_0_yyyzzz_1, g_yzzzzz_0_yyyzzzz_1, g_yzzzzz_0_yyzzzz_1, g_yzzzzz_0_yyzzzzz_1, g_yzzzzz_0_yzzzzz_1, g_yzzzzz_0_yzzzzzz_1, g_yzzzzz_0_zzzzzz_1, g_yzzzzz_0_zzzzzzz_1, g_zzzzz_0_xxxxxxx_0, g_zzzzz_0_xxxxxxx_1, g_zzzzz_0_xxxxxxz_0, g_zzzzz_0_xxxxxxz_1, g_zzzzz_0_xxxxxyz_0, g_zzzzz_0_xxxxxyz_1, g_zzzzz_0_xxxxxzz_0, g_zzzzz_0_xxxxxzz_1, g_zzzzz_0_xxxxyyz_0, g_zzzzz_0_xxxxyyz_1, g_zzzzz_0_xxxxyzz_0, g_zzzzz_0_xxxxyzz_1, g_zzzzz_0_xxxxzzz_0, g_zzzzz_0_xxxxzzz_1, g_zzzzz_0_xxxyyyz_0, g_zzzzz_0_xxxyyyz_1, g_zzzzz_0_xxxyyzz_0, g_zzzzz_0_xxxyyzz_1, g_zzzzz_0_xxxyzzz_0, g_zzzzz_0_xxxyzzz_1, g_zzzzz_0_xxxzzzz_0, g_zzzzz_0_xxxzzzz_1, g_zzzzz_0_xxyyyyz_0, g_zzzzz_0_xxyyyyz_1, g_zzzzz_0_xxyyyzz_0, g_zzzzz_0_xxyyyzz_1, g_zzzzz_0_xxyyzzz_0, g_zzzzz_0_xxyyzzz_1, g_zzzzz_0_xxyzzzz_0, g_zzzzz_0_xxyzzzz_1, g_zzzzz_0_xxzzzzz_0, g_zzzzz_0_xxzzzzz_1, g_zzzzz_0_xyyyyyz_0, g_zzzzz_0_xyyyyyz_1, g_zzzzz_0_xyyyyzz_0, g_zzzzz_0_xyyyyzz_1, g_zzzzz_0_xyyyzzz_0, g_zzzzz_0_xyyyzzz_1, g_zzzzz_0_xyyzzzz_0, g_zzzzz_0_xyyzzzz_1, g_zzzzz_0_xyzzzzz_0, g_zzzzz_0_xyzzzzz_1, g_zzzzz_0_xzzzzzz_0, g_zzzzz_0_xzzzzzz_1, g_zzzzz_0_yyyyyyz_0, g_zzzzz_0_yyyyyyz_1, g_zzzzz_0_yyyyyzz_0, g_zzzzz_0_yyyyyzz_1, g_zzzzz_0_yyyyzzz_0, g_zzzzz_0_yyyyzzz_1, g_zzzzz_0_yyyzzzz_0, g_zzzzz_0_yyyzzzz_1, g_zzzzz_0_yyzzzzz_0, g_zzzzz_0_yyzzzzz_1, g_zzzzz_0_yzzzzzz_0, g_zzzzz_0_yzzzzzz_1, g_zzzzz_0_zzzzzzz_0, g_zzzzz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzzzz_0_xxxxxxx_0[i] = g_zzzzz_0_xxxxxxx_0[i] * fbe_0 - g_zzzzz_0_xxxxxxx_1[i] * fz_be_0 + g_yzzzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxxxy_0[i] = 4.0 * g_yyzzz_0_xxxxxxy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxxxxy_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxxy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxxxxz_0[i] = g_zzzzz_0_xxxxxxz_0[i] * fbe_0 - g_zzzzz_0_xxxxxxz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxxyy_0[i] = 4.0 * g_yyzzz_0_xxxxxyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxxxyy_1[i] * fz_be_0 + g_yyzzzz_0_xxxxxyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxxxyz_0[i] = g_zzzzz_0_xxxxxyz_0[i] * fbe_0 - g_zzzzz_0_xxxxxyz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxxzz_0[i] = g_zzzzz_0_xxxxxzz_0[i] * fbe_0 - g_zzzzz_0_xxxxxzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxyyy_0[i] = 4.0 * g_yyzzz_0_xxxxyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxxyyy_1[i] * fz_be_0 + g_yyzzzz_0_xxxxyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxxyyz_0[i] = g_zzzzz_0_xxxxyyz_0[i] * fbe_0 - g_zzzzz_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxyzz_0[i] = g_zzzzz_0_xxxxyzz_0[i] * fbe_0 - g_zzzzz_0_xxxxyzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxxzzz_0[i] = g_zzzzz_0_xxxxzzz_0[i] * fbe_0 - g_zzzzz_0_xxxxzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxyyyy_0[i] = 4.0 * g_yyzzz_0_xxxyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxxyyyy_1[i] * fz_be_0 + g_yyzzzz_0_xxxyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxxyyyz_0[i] = g_zzzzz_0_xxxyyyz_0[i] * fbe_0 - g_zzzzz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxyyzz_0[i] = g_zzzzz_0_xxxyyzz_0[i] * fbe_0 - g_zzzzz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxyzzz_0[i] = g_zzzzz_0_xxxyzzz_0[i] * fbe_0 - g_zzzzz_0_xxxyzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxxzzzz_0[i] = g_zzzzz_0_xxxzzzz_0[i] * fbe_0 - g_zzzzz_0_xxxzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyyyyy_0[i] = 4.0 * g_yyzzz_0_xxyyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xxyyyyy_1[i] * fz_be_0 + g_yyzzzz_0_xxyyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xxyyyyz_0[i] = g_zzzzz_0_xxyyyyz_0[i] * fbe_0 - g_zzzzz_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyyyzz_0[i] = g_zzzzz_0_xxyyyzz_0[i] * fbe_0 - g_zzzzz_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyyzzz_0[i] = g_zzzzz_0_xxyyzzz_0[i] * fbe_0 - g_zzzzz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxyzzzz_0[i] = g_zzzzz_0_xxyzzzz_0[i] * fbe_0 - g_zzzzz_0_xxyzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xxzzzzz_0[i] = g_zzzzz_0_xxzzzzz_0[i] * fbe_0 - g_zzzzz_0_xxzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyyyyy_0[i] = 4.0 * g_yyzzz_0_xyyyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_xyyyyyy_1[i] * fz_be_0 + g_yyzzzz_0_xyyyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_xyyyyyz_0[i] = g_zzzzz_0_xyyyyyz_0[i] * fbe_0 - g_zzzzz_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yzzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyyyzz_0[i] = g_zzzzz_0_xyyyyzz_0[i] * fbe_0 - g_zzzzz_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyyzzz_0[i] = g_zzzzz_0_xyyyzzz_0[i] * fbe_0 - g_zzzzz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyyzzzz_0[i] = g_zzzzz_0_xyyzzzz_0[i] * fbe_0 - g_zzzzz_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xyzzzzz_0[i] = g_zzzzz_0_xyzzzzz_0[i] * fbe_0 - g_zzzzz_0_xyzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_xzzzzzz_0[i] = g_zzzzz_0_xzzzzzz_0[i] * fbe_0 - g_zzzzz_0_xzzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyyyyy_0[i] = 4.0 * g_yyzzz_0_yyyyyyy_0[i] * fbe_0 - 4.0 * g_yyzzz_0_yyyyyyy_1[i] * fz_be_0 + g_yyzzzz_0_yyyyyyy_1[i] * wa_z[i];

        g_yyzzzzz_0_yyyyyyz_0[i] = g_zzzzz_0_yyyyyyz_0[i] * fbe_0 - g_zzzzz_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yzzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyyyzz_0[i] = g_zzzzz_0_yyyyyzz_0[i] * fbe_0 - g_zzzzz_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yzzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyyzzz_0[i] = g_zzzzz_0_yyyyzzz_0[i] * fbe_0 - g_zzzzz_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yzzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyyzzzz_0[i] = g_zzzzz_0_yyyzzzz_0[i] * fbe_0 - g_zzzzz_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yzzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yyzzzzz_0[i] = g_zzzzz_0_yyzzzzz_0[i] * fbe_0 - g_zzzzz_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_yzzzzzz_0[i] = g_zzzzz_0_yzzzzzz_0[i] * fbe_0 - g_zzzzz_0_yzzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_yzzzzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yyzzzzz_0_zzzzzzz_0[i] = g_zzzzz_0_zzzzzzz_0[i] * fbe_0 - g_zzzzz_0_zzzzzzz_1[i] * fz_be_0 + g_yzzzzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1224-1260 components of targeted buffer : KSK

    auto g_yzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 1224);

    auto g_yzzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 1225);

    auto g_yzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 1226);

    auto g_yzzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 1227);

    auto g_yzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 1228);

    auto g_yzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 1229);

    auto g_yzzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 1230);

    auto g_yzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 1231);

    auto g_yzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 1232);

    auto g_yzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 1233);

    auto g_yzzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1234);

    auto g_yzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1235);

    auto g_yzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1236);

    auto g_yzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1237);

    auto g_yzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1238);

    auto g_yzzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1239);

    auto g_yzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1240);

    auto g_yzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1241);

    auto g_yzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1242);

    auto g_yzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1243);

    auto g_yzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1244);

    auto g_yzzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1245);

    auto g_yzzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1246);

    auto g_yzzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1247);

    auto g_yzzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1248);

    auto g_yzzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1249);

    auto g_yzzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1250);

    auto g_yzzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1251);

    auto g_yzzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1252);

    auto g_yzzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1253);

    auto g_yzzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1254);

    auto g_yzzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1255);

    auto g_yzzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1256);

    auto g_yzzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1257);

    auto g_yzzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1258);

    auto g_yzzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1259);

    #pragma omp simd aligned(g_yzzzzzz_0_xxxxxxx_0, g_yzzzzzz_0_xxxxxxy_0, g_yzzzzzz_0_xxxxxxz_0, g_yzzzzzz_0_xxxxxyy_0, g_yzzzzzz_0_xxxxxyz_0, g_yzzzzzz_0_xxxxxzz_0, g_yzzzzzz_0_xxxxyyy_0, g_yzzzzzz_0_xxxxyyz_0, g_yzzzzzz_0_xxxxyzz_0, g_yzzzzzz_0_xxxxzzz_0, g_yzzzzzz_0_xxxyyyy_0, g_yzzzzzz_0_xxxyyyz_0, g_yzzzzzz_0_xxxyyzz_0, g_yzzzzzz_0_xxxyzzz_0, g_yzzzzzz_0_xxxzzzz_0, g_yzzzzzz_0_xxyyyyy_0, g_yzzzzzz_0_xxyyyyz_0, g_yzzzzzz_0_xxyyyzz_0, g_yzzzzzz_0_xxyyzzz_0, g_yzzzzzz_0_xxyzzzz_0, g_yzzzzzz_0_xxzzzzz_0, g_yzzzzzz_0_xyyyyyy_0, g_yzzzzzz_0_xyyyyyz_0, g_yzzzzzz_0_xyyyyzz_0, g_yzzzzzz_0_xyyyzzz_0, g_yzzzzzz_0_xyyzzzz_0, g_yzzzzzz_0_xyzzzzz_0, g_yzzzzzz_0_xzzzzzz_0, g_yzzzzzz_0_yyyyyyy_0, g_yzzzzzz_0_yyyyyyz_0, g_yzzzzzz_0_yyyyyzz_0, g_yzzzzzz_0_yyyyzzz_0, g_yzzzzzz_0_yyyzzzz_0, g_yzzzzzz_0_yyzzzzz_0, g_yzzzzzz_0_yzzzzzz_0, g_yzzzzzz_0_zzzzzzz_0, g_zzzzzz_0_xxxxxx_1, g_zzzzzz_0_xxxxxxx_1, g_zzzzzz_0_xxxxxxy_1, g_zzzzzz_0_xxxxxxz_1, g_zzzzzz_0_xxxxxy_1, g_zzzzzz_0_xxxxxyy_1, g_zzzzzz_0_xxxxxyz_1, g_zzzzzz_0_xxxxxz_1, g_zzzzzz_0_xxxxxzz_1, g_zzzzzz_0_xxxxyy_1, g_zzzzzz_0_xxxxyyy_1, g_zzzzzz_0_xxxxyyz_1, g_zzzzzz_0_xxxxyz_1, g_zzzzzz_0_xxxxyzz_1, g_zzzzzz_0_xxxxzz_1, g_zzzzzz_0_xxxxzzz_1, g_zzzzzz_0_xxxyyy_1, g_zzzzzz_0_xxxyyyy_1, g_zzzzzz_0_xxxyyyz_1, g_zzzzzz_0_xxxyyz_1, g_zzzzzz_0_xxxyyzz_1, g_zzzzzz_0_xxxyzz_1, g_zzzzzz_0_xxxyzzz_1, g_zzzzzz_0_xxxzzz_1, g_zzzzzz_0_xxxzzzz_1, g_zzzzzz_0_xxyyyy_1, g_zzzzzz_0_xxyyyyy_1, g_zzzzzz_0_xxyyyyz_1, g_zzzzzz_0_xxyyyz_1, g_zzzzzz_0_xxyyyzz_1, g_zzzzzz_0_xxyyzz_1, g_zzzzzz_0_xxyyzzz_1, g_zzzzzz_0_xxyzzz_1, g_zzzzzz_0_xxyzzzz_1, g_zzzzzz_0_xxzzzz_1, g_zzzzzz_0_xxzzzzz_1, g_zzzzzz_0_xyyyyy_1, g_zzzzzz_0_xyyyyyy_1, g_zzzzzz_0_xyyyyyz_1, g_zzzzzz_0_xyyyyz_1, g_zzzzzz_0_xyyyyzz_1, g_zzzzzz_0_xyyyzz_1, g_zzzzzz_0_xyyyzzz_1, g_zzzzzz_0_xyyzzz_1, g_zzzzzz_0_xyyzzzz_1, g_zzzzzz_0_xyzzzz_1, g_zzzzzz_0_xyzzzzz_1, g_zzzzzz_0_xzzzzz_1, g_zzzzzz_0_xzzzzzz_1, g_zzzzzz_0_yyyyyy_1, g_zzzzzz_0_yyyyyyy_1, g_zzzzzz_0_yyyyyyz_1, g_zzzzzz_0_yyyyyz_1, g_zzzzzz_0_yyyyyzz_1, g_zzzzzz_0_yyyyzz_1, g_zzzzzz_0_yyyyzzz_1, g_zzzzzz_0_yyyzzz_1, g_zzzzzz_0_yyyzzzz_1, g_zzzzzz_0_yyzzzz_1, g_zzzzzz_0_yyzzzzz_1, g_zzzzzz_0_yzzzzz_1, g_zzzzzz_0_yzzzzzz_1, g_zzzzzz_0_zzzzzz_1, g_zzzzzz_0_zzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzzz_0_xxxxxxx_0[i] = g_zzzzzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxxy_0[i] = g_zzzzzz_0_xxxxxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxxz_0[i] = g_zzzzzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxyy_0[i] = 2.0 * g_zzzzzz_0_xxxxxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxyz_0[i] = g_zzzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxxzz_0[i] = g_zzzzzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxyyy_0[i] = 3.0 * g_zzzzzz_0_xxxxyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxyyz_0[i] = 2.0 * g_zzzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxyzz_0[i] = g_zzzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxxzzz_0[i] = g_zzzzzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyyyy_0[i] = 4.0 * g_zzzzzz_0_xxxyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyyyz_0[i] = 3.0 * g_zzzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyyzz_0[i] = 2.0 * g_zzzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxyzzz_0[i] = g_zzzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxxzzzz_0[i] = g_zzzzzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyyyy_0[i] = 5.0 * g_zzzzzz_0_xxyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyyyz_0[i] = 4.0 * g_zzzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyyzz_0[i] = 3.0 * g_zzzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyyzzz_0[i] = 2.0 * g_zzzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxyzzzz_0[i] = g_zzzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xxzzzzz_0[i] = g_zzzzzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyyyy_0[i] = 6.0 * g_zzzzzz_0_xyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyyyz_0[i] = 5.0 * g_zzzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyyzz_0[i] = 4.0 * g_zzzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyyzzz_0[i] = 3.0 * g_zzzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyyzzzz_0[i] = 2.0 * g_zzzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xyzzzzz_0[i] = g_zzzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_xzzzzzz_0[i] = g_zzzzzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyyyy_0[i] = 7.0 * g_zzzzzz_0_yyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyyy_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyyyz_0[i] = 6.0 * g_zzzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyyzz_0[i] = 5.0 * g_zzzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyyzzz_0[i] = 4.0 * g_zzzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyyzzzz_0[i] = 3.0 * g_zzzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yyzzzzz_0[i] = 2.0 * g_zzzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_yzzzzzz_0[i] = g_zzzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yzzzzzz_0_zzzzzzz_0[i] = g_zzzzzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 1260-1296 components of targeted buffer : KSK

    auto g_zzzzzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ksk + 1260);

    auto g_zzzzzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_ksk + 1261);

    auto g_zzzzzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_ksk + 1262);

    auto g_zzzzzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ksk + 1263);

    auto g_zzzzzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_ksk + 1264);

    auto g_zzzzzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ksk + 1265);

    auto g_zzzzzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ksk + 1266);

    auto g_zzzzzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_ksk + 1267);

    auto g_zzzzzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_ksk + 1268);

    auto g_zzzzzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ksk + 1269);

    auto g_zzzzzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1270);

    auto g_zzzzzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1271);

    auto g_zzzzzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1272);

    auto g_zzzzzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1273);

    auto g_zzzzzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1274);

    auto g_zzzzzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1275);

    auto g_zzzzzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1276);

    auto g_zzzzzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1277);

    auto g_zzzzzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1278);

    auto g_zzzzzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1279);

    auto g_zzzzzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1280);

    auto g_zzzzzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1281);

    auto g_zzzzzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1282);

    auto g_zzzzzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1283);

    auto g_zzzzzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1284);

    auto g_zzzzzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1285);

    auto g_zzzzzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1286);

    auto g_zzzzzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1287);

    auto g_zzzzzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ksk + 1288);

    auto g_zzzzzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_ksk + 1289);

    auto g_zzzzzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ksk + 1290);

    auto g_zzzzzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ksk + 1291);

    auto g_zzzzzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1292);

    auto g_zzzzzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1293);

    auto g_zzzzzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1294);

    auto g_zzzzzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ksk + 1295);

    #pragma omp simd aligned(g_zzzzz_0_xxxxxxx_0, g_zzzzz_0_xxxxxxx_1, g_zzzzz_0_xxxxxxy_0, g_zzzzz_0_xxxxxxy_1, g_zzzzz_0_xxxxxxz_0, g_zzzzz_0_xxxxxxz_1, g_zzzzz_0_xxxxxyy_0, g_zzzzz_0_xxxxxyy_1, g_zzzzz_0_xxxxxyz_0, g_zzzzz_0_xxxxxyz_1, g_zzzzz_0_xxxxxzz_0, g_zzzzz_0_xxxxxzz_1, g_zzzzz_0_xxxxyyy_0, g_zzzzz_0_xxxxyyy_1, g_zzzzz_0_xxxxyyz_0, g_zzzzz_0_xxxxyyz_1, g_zzzzz_0_xxxxyzz_0, g_zzzzz_0_xxxxyzz_1, g_zzzzz_0_xxxxzzz_0, g_zzzzz_0_xxxxzzz_1, g_zzzzz_0_xxxyyyy_0, g_zzzzz_0_xxxyyyy_1, g_zzzzz_0_xxxyyyz_0, g_zzzzz_0_xxxyyyz_1, g_zzzzz_0_xxxyyzz_0, g_zzzzz_0_xxxyyzz_1, g_zzzzz_0_xxxyzzz_0, g_zzzzz_0_xxxyzzz_1, g_zzzzz_0_xxxzzzz_0, g_zzzzz_0_xxxzzzz_1, g_zzzzz_0_xxyyyyy_0, g_zzzzz_0_xxyyyyy_1, g_zzzzz_0_xxyyyyz_0, g_zzzzz_0_xxyyyyz_1, g_zzzzz_0_xxyyyzz_0, g_zzzzz_0_xxyyyzz_1, g_zzzzz_0_xxyyzzz_0, g_zzzzz_0_xxyyzzz_1, g_zzzzz_0_xxyzzzz_0, g_zzzzz_0_xxyzzzz_1, g_zzzzz_0_xxzzzzz_0, g_zzzzz_0_xxzzzzz_1, g_zzzzz_0_xyyyyyy_0, g_zzzzz_0_xyyyyyy_1, g_zzzzz_0_xyyyyyz_0, g_zzzzz_0_xyyyyyz_1, g_zzzzz_0_xyyyyzz_0, g_zzzzz_0_xyyyyzz_1, g_zzzzz_0_xyyyzzz_0, g_zzzzz_0_xyyyzzz_1, g_zzzzz_0_xyyzzzz_0, g_zzzzz_0_xyyzzzz_1, g_zzzzz_0_xyzzzzz_0, g_zzzzz_0_xyzzzzz_1, g_zzzzz_0_xzzzzzz_0, g_zzzzz_0_xzzzzzz_1, g_zzzzz_0_yyyyyyy_0, g_zzzzz_0_yyyyyyy_1, g_zzzzz_0_yyyyyyz_0, g_zzzzz_0_yyyyyyz_1, g_zzzzz_0_yyyyyzz_0, g_zzzzz_0_yyyyyzz_1, g_zzzzz_0_yyyyzzz_0, g_zzzzz_0_yyyyzzz_1, g_zzzzz_0_yyyzzzz_0, g_zzzzz_0_yyyzzzz_1, g_zzzzz_0_yyzzzzz_0, g_zzzzz_0_yyzzzzz_1, g_zzzzz_0_yzzzzzz_0, g_zzzzz_0_yzzzzzz_1, g_zzzzz_0_zzzzzzz_0, g_zzzzz_0_zzzzzzz_1, g_zzzzzz_0_xxxxxx_1, g_zzzzzz_0_xxxxxxx_1, g_zzzzzz_0_xxxxxxy_1, g_zzzzzz_0_xxxxxxz_1, g_zzzzzz_0_xxxxxy_1, g_zzzzzz_0_xxxxxyy_1, g_zzzzzz_0_xxxxxyz_1, g_zzzzzz_0_xxxxxz_1, g_zzzzzz_0_xxxxxzz_1, g_zzzzzz_0_xxxxyy_1, g_zzzzzz_0_xxxxyyy_1, g_zzzzzz_0_xxxxyyz_1, g_zzzzzz_0_xxxxyz_1, g_zzzzzz_0_xxxxyzz_1, g_zzzzzz_0_xxxxzz_1, g_zzzzzz_0_xxxxzzz_1, g_zzzzzz_0_xxxyyy_1, g_zzzzzz_0_xxxyyyy_1, g_zzzzzz_0_xxxyyyz_1, g_zzzzzz_0_xxxyyz_1, g_zzzzzz_0_xxxyyzz_1, g_zzzzzz_0_xxxyzz_1, g_zzzzzz_0_xxxyzzz_1, g_zzzzzz_0_xxxzzz_1, g_zzzzzz_0_xxxzzzz_1, g_zzzzzz_0_xxyyyy_1, g_zzzzzz_0_xxyyyyy_1, g_zzzzzz_0_xxyyyyz_1, g_zzzzzz_0_xxyyyz_1, g_zzzzzz_0_xxyyyzz_1, g_zzzzzz_0_xxyyzz_1, g_zzzzzz_0_xxyyzzz_1, g_zzzzzz_0_xxyzzz_1, g_zzzzzz_0_xxyzzzz_1, g_zzzzzz_0_xxzzzz_1, g_zzzzzz_0_xxzzzzz_1, g_zzzzzz_0_xyyyyy_1, g_zzzzzz_0_xyyyyyy_1, g_zzzzzz_0_xyyyyyz_1, g_zzzzzz_0_xyyyyz_1, g_zzzzzz_0_xyyyyzz_1, g_zzzzzz_0_xyyyzz_1, g_zzzzzz_0_xyyyzzz_1, g_zzzzzz_0_xyyzzz_1, g_zzzzzz_0_xyyzzzz_1, g_zzzzzz_0_xyzzzz_1, g_zzzzzz_0_xyzzzzz_1, g_zzzzzz_0_xzzzzz_1, g_zzzzzz_0_xzzzzzz_1, g_zzzzzz_0_yyyyyy_1, g_zzzzzz_0_yyyyyyy_1, g_zzzzzz_0_yyyyyyz_1, g_zzzzzz_0_yyyyyz_1, g_zzzzzz_0_yyyyyzz_1, g_zzzzzz_0_yyyyzz_1, g_zzzzzz_0_yyyyzzz_1, g_zzzzzz_0_yyyzzz_1, g_zzzzzz_0_yyyzzzz_1, g_zzzzzz_0_yyzzzz_1, g_zzzzzz_0_yyzzzzz_1, g_zzzzzz_0_yzzzzz_1, g_zzzzzz_0_yzzzzzz_1, g_zzzzzz_0_zzzzzz_1, g_zzzzzz_0_zzzzzzz_1, g_zzzzzzz_0_xxxxxxx_0, g_zzzzzzz_0_xxxxxxy_0, g_zzzzzzz_0_xxxxxxz_0, g_zzzzzzz_0_xxxxxyy_0, g_zzzzzzz_0_xxxxxyz_0, g_zzzzzzz_0_xxxxxzz_0, g_zzzzzzz_0_xxxxyyy_0, g_zzzzzzz_0_xxxxyyz_0, g_zzzzzzz_0_xxxxyzz_0, g_zzzzzzz_0_xxxxzzz_0, g_zzzzzzz_0_xxxyyyy_0, g_zzzzzzz_0_xxxyyyz_0, g_zzzzzzz_0_xxxyyzz_0, g_zzzzzzz_0_xxxyzzz_0, g_zzzzzzz_0_xxxzzzz_0, g_zzzzzzz_0_xxyyyyy_0, g_zzzzzzz_0_xxyyyyz_0, g_zzzzzzz_0_xxyyyzz_0, g_zzzzzzz_0_xxyyzzz_0, g_zzzzzzz_0_xxyzzzz_0, g_zzzzzzz_0_xxzzzzz_0, g_zzzzzzz_0_xyyyyyy_0, g_zzzzzzz_0_xyyyyyz_0, g_zzzzzzz_0_xyyyyzz_0, g_zzzzzzz_0_xyyyzzz_0, g_zzzzzzz_0_xyyzzzz_0, g_zzzzzzz_0_xyzzzzz_0, g_zzzzzzz_0_xzzzzzz_0, g_zzzzzzz_0_yyyyyyy_0, g_zzzzzzz_0_yyyyyyz_0, g_zzzzzzz_0_yyyyyzz_0, g_zzzzzzz_0_yyyyzzz_0, g_zzzzzzz_0_yyyzzzz_0, g_zzzzzzz_0_yyzzzzz_0, g_zzzzzzz_0_yzzzzzz_0, g_zzzzzzz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzzz_0_xxxxxxx_0[i] = 6.0 * g_zzzzz_0_xxxxxxx_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxxx_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxxx_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxxy_0[i] = 6.0 * g_zzzzz_0_xxxxxxy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxxy_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxxy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxxz_0[i] = 6.0 * g_zzzzz_0_xxxxxxz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxxz_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxx_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxxz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxyy_0[i] = 6.0 * g_zzzzz_0_xxxxxyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxyy_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxyz_0[i] = 6.0 * g_zzzzz_0_xxxxxyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxyz_1[i] * fz_be_0 + g_zzzzzz_0_xxxxxy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxxzz_0[i] = 6.0 * g_zzzzz_0_xxxxxzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxxzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxxxxz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxxzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxyyy_0[i] = 6.0 * g_zzzzz_0_xxxxyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxyyy_1[i] * fz_be_0 + g_zzzzzz_0_xxxxyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxyyz_0[i] = 6.0 * g_zzzzz_0_xxxxyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxyyz_1[i] * fz_be_0 + g_zzzzzz_0_xxxxyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxyzz_0[i] = 6.0 * g_zzzzz_0_xxxxyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxxxyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxxzzz_0[i] = 6.0 * g_zzzzz_0_xxxxzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxxzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xxxxzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxxzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyyyy_0[i] = 6.0 * g_zzzzz_0_xxxyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyyyy_1[i] * fz_be_0 + g_zzzzzz_0_xxxyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyyyz_0[i] = 6.0 * g_zzzzz_0_xxxyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyyyz_1[i] * fz_be_0 + g_zzzzzz_0_xxxyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyyzz_0[i] = 6.0 * g_zzzzz_0_xxxyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxxyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxyzzz_0[i] = 6.0 * g_zzzzz_0_xxxyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xxxyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxxzzzz_0[i] = 6.0 * g_zzzzz_0_xxxzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxxzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_xxxzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxxzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyyyy_0[i] = 6.0 * g_zzzzz_0_xxyyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyyyy_1[i] * fz_be_0 + g_zzzzzz_0_xxyyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyyyz_0[i] = 6.0 * g_zzzzz_0_xxyyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyyyz_1[i] * fz_be_0 + g_zzzzzz_0_xxyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyyzz_0[i] = 6.0 * g_zzzzz_0_xxyyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xxyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyyzzz_0[i] = 6.0 * g_zzzzz_0_xxyyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xxyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxyzzzz_0[i] = 6.0 * g_zzzzz_0_xxyzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_xxyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxyzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xxzzzzz_0[i] = 6.0 * g_zzzzz_0_xxzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xxzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_0_xxzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xxzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyyyy_0[i] = 6.0 * g_zzzzz_0_xyyyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyyyy_1[i] * fz_be_0 + g_zzzzzz_0_xyyyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyyyz_0[i] = 6.0 * g_zzzzz_0_xyyyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyyyz_1[i] * fz_be_0 + g_zzzzzz_0_xyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyyzz_0[i] = 6.0 * g_zzzzz_0_xyyyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_xyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyyzzz_0[i] = 6.0 * g_zzzzz_0_xyyyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_xyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyyzzzz_0[i] = 6.0 * g_zzzzz_0_xyyzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_xyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyyzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xyzzzzz_0[i] = 6.0 * g_zzzzz_0_xyzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_0_xyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xyzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_xzzzzzz_0[i] = 6.0 * g_zzzzz_0_xzzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_xzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzzz_0_xzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_xzzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyyyy_0[i] = 6.0 * g_zzzzz_0_yyyyyyy_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyyyy_1[i] * fz_be_0 + g_zzzzzz_0_yyyyyyy_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyyyz_0[i] = 6.0 * g_zzzzz_0_yyyyyyz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyyyz_1[i] * fz_be_0 + g_zzzzzz_0_yyyyyy_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyyz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyyzz_0[i] = 6.0 * g_zzzzz_0_yyyyyzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_0_yyyyyz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyyzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyyzzz_0[i] = 6.0 * g_zzzzz_0_yyyyzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_0_yyyyzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyyzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyyzzzz_0[i] = 6.0 * g_zzzzz_0_yyyzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzzzz_0_yyyzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyyzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yyzzzzz_0[i] = 6.0 * g_zzzzz_0_yyzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzzzz_0_yyzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yyzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_yzzzzzz_0[i] = 6.0 * g_zzzzz_0_yzzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_yzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzzzz_0_yzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_yzzzzzz_1[i] * wa_z[i];

        g_zzzzzzz_0_zzzzzzz_0[i] = 6.0 * g_zzzzz_0_zzzzzzz_0[i] * fbe_0 - 6.0 * g_zzzzz_0_zzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzzzzz_0_zzzzzz_1[i] * fi_acd_0 + g_zzzzzz_0_zzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

