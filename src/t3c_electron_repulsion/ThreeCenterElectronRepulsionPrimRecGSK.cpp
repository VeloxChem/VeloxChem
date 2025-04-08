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

#include "ThreeCenterElectronRepulsionPrimRecGSK.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_gsk(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsk,
                                 size_t idx_eri_0_dsk,
                                 size_t idx_eri_1_dsk,
                                 size_t idx_eri_1_fsi,
                                 size_t idx_eri_1_fsk,
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

    /// Set up components of auxilary buffer : DSK

    auto g_xx_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_dsk);

    auto g_xx_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_dsk + 1);

    auto g_xx_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_dsk + 2);

    auto g_xx_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_dsk + 3);

    auto g_xx_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_dsk + 4);

    auto g_xx_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_dsk + 5);

    auto g_xx_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_dsk + 6);

    auto g_xx_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_dsk + 7);

    auto g_xx_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_dsk + 8);

    auto g_xx_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_dsk + 9);

    auto g_xx_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_dsk + 10);

    auto g_xx_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_dsk + 11);

    auto g_xx_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_dsk + 12);

    auto g_xx_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_dsk + 13);

    auto g_xx_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_dsk + 14);

    auto g_xx_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 15);

    auto g_xx_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 16);

    auto g_xx_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 17);

    auto g_xx_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 18);

    auto g_xx_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 19);

    auto g_xx_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 20);

    auto g_xx_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 21);

    auto g_xx_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 22);

    auto g_xx_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 23);

    auto g_xx_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 24);

    auto g_xx_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 25);

    auto g_xx_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 26);

    auto g_xx_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 27);

    auto g_xx_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 28);

    auto g_xx_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 29);

    auto g_xx_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 30);

    auto g_xx_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 31);

    auto g_xx_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 32);

    auto g_xx_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 33);

    auto g_xx_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 34);

    auto g_xx_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 35);

    auto g_yy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_dsk + 108);

    auto g_yy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_dsk + 109);

    auto g_yy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_dsk + 110);

    auto g_yy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_dsk + 111);

    auto g_yy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_dsk + 112);

    auto g_yy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_dsk + 113);

    auto g_yy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_dsk + 114);

    auto g_yy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_dsk + 115);

    auto g_yy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_dsk + 116);

    auto g_yy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_dsk + 117);

    auto g_yy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_dsk + 118);

    auto g_yy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_dsk + 119);

    auto g_yy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_dsk + 120);

    auto g_yy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_dsk + 121);

    auto g_yy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_dsk + 122);

    auto g_yy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 123);

    auto g_yy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 124);

    auto g_yy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 125);

    auto g_yy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 126);

    auto g_yy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 127);

    auto g_yy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 128);

    auto g_yy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 129);

    auto g_yy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 130);

    auto g_yy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 131);

    auto g_yy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 132);

    auto g_yy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 133);

    auto g_yy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 134);

    auto g_yy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 135);

    auto g_yy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 136);

    auto g_yy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 137);

    auto g_yy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 138);

    auto g_yy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 139);

    auto g_yy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 140);

    auto g_yy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 141);

    auto g_yy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 142);

    auto g_yy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 143);

    auto g_zz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_dsk + 180);

    auto g_zz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_dsk + 181);

    auto g_zz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_dsk + 182);

    auto g_zz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_dsk + 183);

    auto g_zz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_dsk + 184);

    auto g_zz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_dsk + 185);

    auto g_zz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_dsk + 186);

    auto g_zz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_dsk + 187);

    auto g_zz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_dsk + 188);

    auto g_zz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_dsk + 189);

    auto g_zz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_dsk + 190);

    auto g_zz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_dsk + 191);

    auto g_zz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_dsk + 192);

    auto g_zz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_dsk + 193);

    auto g_zz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_dsk + 194);

    auto g_zz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 195);

    auto g_zz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 196);

    auto g_zz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 197);

    auto g_zz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 198);

    auto g_zz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 199);

    auto g_zz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 200);

    auto g_zz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 201);

    auto g_zz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 202);

    auto g_zz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 203);

    auto g_zz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 204);

    auto g_zz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 205);

    auto g_zz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 206);

    auto g_zz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 207);

    auto g_zz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_dsk + 208);

    auto g_zz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_dsk + 209);

    auto g_zz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_dsk + 210);

    auto g_zz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_dsk + 211);

    auto g_zz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_dsk + 212);

    auto g_zz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 213);

    auto g_zz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 214);

    auto g_zz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_dsk + 215);

    /// Set up components of auxilary buffer : DSK

    auto g_xx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_dsk);

    auto g_xx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_dsk + 1);

    auto g_xx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_dsk + 2);

    auto g_xx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_dsk + 3);

    auto g_xx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_dsk + 4);

    auto g_xx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_dsk + 5);

    auto g_xx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_dsk + 6);

    auto g_xx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_dsk + 7);

    auto g_xx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_dsk + 8);

    auto g_xx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_dsk + 9);

    auto g_xx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_dsk + 10);

    auto g_xx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_dsk + 11);

    auto g_xx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_dsk + 12);

    auto g_xx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_dsk + 13);

    auto g_xx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_dsk + 14);

    auto g_xx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 15);

    auto g_xx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 16);

    auto g_xx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 17);

    auto g_xx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 18);

    auto g_xx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 19);

    auto g_xx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 20);

    auto g_xx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 21);

    auto g_xx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 22);

    auto g_xx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 23);

    auto g_xx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 24);

    auto g_xx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 25);

    auto g_xx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 26);

    auto g_xx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 27);

    auto g_xx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 28);

    auto g_xx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 29);

    auto g_xx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 30);

    auto g_xx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 31);

    auto g_xx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 32);

    auto g_xx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 33);

    auto g_xx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 34);

    auto g_xx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 35);

    auto g_yy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_dsk + 108);

    auto g_yy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_dsk + 109);

    auto g_yy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_dsk + 110);

    auto g_yy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_dsk + 111);

    auto g_yy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_dsk + 112);

    auto g_yy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_dsk + 113);

    auto g_yy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_dsk + 114);

    auto g_yy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_dsk + 115);

    auto g_yy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_dsk + 116);

    auto g_yy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_dsk + 117);

    auto g_yy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_dsk + 118);

    auto g_yy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_dsk + 119);

    auto g_yy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_dsk + 120);

    auto g_yy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_dsk + 121);

    auto g_yy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_dsk + 122);

    auto g_yy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 123);

    auto g_yy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 124);

    auto g_yy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 125);

    auto g_yy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 126);

    auto g_yy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 127);

    auto g_yy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 128);

    auto g_yy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 129);

    auto g_yy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 130);

    auto g_yy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 131);

    auto g_yy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 132);

    auto g_yy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 133);

    auto g_yy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 134);

    auto g_yy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 135);

    auto g_yy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 136);

    auto g_yy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 137);

    auto g_yy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 138);

    auto g_yy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 139);

    auto g_yy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 140);

    auto g_yy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 141);

    auto g_yy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 142);

    auto g_yy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 143);

    auto g_zz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_dsk + 180);

    auto g_zz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_dsk + 181);

    auto g_zz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_dsk + 182);

    auto g_zz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_dsk + 183);

    auto g_zz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_dsk + 184);

    auto g_zz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_dsk + 185);

    auto g_zz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_dsk + 186);

    auto g_zz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_dsk + 187);

    auto g_zz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_dsk + 188);

    auto g_zz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_dsk + 189);

    auto g_zz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_dsk + 190);

    auto g_zz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_dsk + 191);

    auto g_zz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_dsk + 192);

    auto g_zz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_dsk + 193);

    auto g_zz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_dsk + 194);

    auto g_zz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 195);

    auto g_zz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 196);

    auto g_zz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 197);

    auto g_zz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 198);

    auto g_zz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 199);

    auto g_zz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 200);

    auto g_zz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 201);

    auto g_zz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 202);

    auto g_zz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 203);

    auto g_zz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 204);

    auto g_zz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 205);

    auto g_zz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 206);

    auto g_zz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 207);

    auto g_zz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_dsk + 208);

    auto g_zz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_dsk + 209);

    auto g_zz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_dsk + 210);

    auto g_zz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_dsk + 211);

    auto g_zz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_dsk + 212);

    auto g_zz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 213);

    auto g_zz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 214);

    auto g_zz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_dsk + 215);

    /// Set up components of auxilary buffer : FSI

    auto g_xxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi);

    auto g_xxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 1);

    auto g_xxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 2);

    auto g_xxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 3);

    auto g_xxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 4);

    auto g_xxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 5);

    auto g_xxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 6);

    auto g_xxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 7);

    auto g_xxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 8);

    auto g_xxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 9);

    auto g_xxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 10);

    auto g_xxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 11);

    auto g_xxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 12);

    auto g_xxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 13);

    auto g_xxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 14);

    auto g_xxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 15);

    auto g_xxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 16);

    auto g_xxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 17);

    auto g_xxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 18);

    auto g_xxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 19);

    auto g_xxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 20);

    auto g_xxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 21);

    auto g_xxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 22);

    auto g_xxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 23);

    auto g_xxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 24);

    auto g_xxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 25);

    auto g_xxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 26);

    auto g_xxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 27);

    auto g_xxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 58);

    auto g_xxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 60);

    auto g_xxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 61);

    auto g_xxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 63);

    auto g_xxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 64);

    auto g_xxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 65);

    auto g_xxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 67);

    auto g_xxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 68);

    auto g_xxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 69);

    auto g_xxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 70);

    auto g_xxz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 72);

    auto g_xxz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 73);

    auto g_xxz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 74);

    auto g_xxz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 75);

    auto g_xxz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 76);

    auto g_xxz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 78);

    auto g_xxz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 79);

    auto g_xxz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 80);

    auto g_xxz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 81);

    auto g_xxz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 82);

    auto g_xxz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 83);

    auto g_xyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 85);

    auto g_xyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 87);

    auto g_xyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 88);

    auto g_xyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 90);

    auto g_xyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 91);

    auto g_xyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 92);

    auto g_xyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 94);

    auto g_xyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 95);

    auto g_xyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 96);

    auto g_xyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 97);

    auto g_xyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 99);

    auto g_xyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 100);

    auto g_xyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 101);

    auto g_xyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 102);

    auto g_xyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 103);

    auto g_xyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 105);

    auto g_xyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 106);

    auto g_xyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 107);

    auto g_xyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 108);

    auto g_xyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 109);

    auto g_xyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 110);

    auto g_xzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 142);

    auto g_xzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 144);

    auto g_xzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 145);

    auto g_xzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 147);

    auto g_xzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 148);

    auto g_xzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 149);

    auto g_xzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 151);

    auto g_xzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 152);

    auto g_xzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 153);

    auto g_xzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 154);

    auto g_xzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 156);

    auto g_xzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 157);

    auto g_xzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 158);

    auto g_xzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 159);

    auto g_xzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 160);

    auto g_xzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 162);

    auto g_xzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 163);

    auto g_xzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 164);

    auto g_xzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 165);

    auto g_xzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 166);

    auto g_xzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 167);

    auto g_yyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi + 168);

    auto g_yyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 169);

    auto g_yyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 170);

    auto g_yyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 171);

    auto g_yyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 172);

    auto g_yyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 173);

    auto g_yyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 174);

    auto g_yyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 175);

    auto g_yyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 176);

    auto g_yyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 177);

    auto g_yyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 178);

    auto g_yyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 179);

    auto g_yyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 180);

    auto g_yyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 181);

    auto g_yyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 182);

    auto g_yyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 183);

    auto g_yyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 184);

    auto g_yyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 185);

    auto g_yyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 186);

    auto g_yyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 187);

    auto g_yyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 188);

    auto g_yyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 189);

    auto g_yyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 190);

    auto g_yyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 191);

    auto g_yyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 192);

    auto g_yyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 193);

    auto g_yyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 194);

    auto g_yyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 195);

    auto g_yyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 198);

    auto g_yyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 200);

    auto g_yyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 201);

    auto g_yyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 203);

    auto g_yyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 204);

    auto g_yyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 205);

    auto g_yyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 207);

    auto g_yyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 208);

    auto g_yyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 209);

    auto g_yyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 210);

    auto g_yyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 212);

    auto g_yyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 213);

    auto g_yyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 214);

    auto g_yyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 215);

    auto g_yyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 216);

    auto g_yyz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 218);

    auto g_yyz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 219);

    auto g_yyz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 220);

    auto g_yyz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 221);

    auto g_yyz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 222);

    auto g_yyz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 223);

    auto g_yzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 225);

    auto g_yzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 226);

    auto g_yzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 227);

    auto g_yzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 228);

    auto g_yzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 229);

    auto g_yzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 230);

    auto g_yzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 231);

    auto g_yzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 232);

    auto g_yzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 233);

    auto g_yzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 234);

    auto g_yzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 235);

    auto g_yzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 236);

    auto g_yzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 237);

    auto g_yzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 238);

    auto g_yzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 239);

    auto g_yzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 240);

    auto g_yzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 241);

    auto g_yzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 242);

    auto g_yzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 243);

    auto g_yzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 244);

    auto g_yzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 245);

    auto g_yzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 246);

    auto g_yzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 247);

    auto g_yzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 248);

    auto g_yzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 249);

    auto g_yzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 250);

    auto g_yzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 251);

    auto g_zzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi + 252);

    auto g_zzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 253);

    auto g_zzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 254);

    auto g_zzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 255);

    auto g_zzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 256);

    auto g_zzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 257);

    auto g_zzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 258);

    auto g_zzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 259);

    auto g_zzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 260);

    auto g_zzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 261);

    auto g_zzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 262);

    auto g_zzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 263);

    auto g_zzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 264);

    auto g_zzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 265);

    auto g_zzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 266);

    auto g_zzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 267);

    auto g_zzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 268);

    auto g_zzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 269);

    auto g_zzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 270);

    auto g_zzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 271);

    auto g_zzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 272);

    auto g_zzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 273);

    auto g_zzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 274);

    auto g_zzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 275);

    auto g_zzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 276);

    auto g_zzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 277);

    auto g_zzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 278);

    auto g_zzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 279);

    /// Set up components of auxilary buffer : FSK

    auto g_xxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk);

    auto g_xxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 1);

    auto g_xxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 2);

    auto g_xxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 3);

    auto g_xxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 4);

    auto g_xxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 5);

    auto g_xxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 6);

    auto g_xxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 7);

    auto g_xxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 8);

    auto g_xxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 9);

    auto g_xxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 10);

    auto g_xxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 11);

    auto g_xxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 12);

    auto g_xxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 13);

    auto g_xxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 14);

    auto g_xxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 15);

    auto g_xxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 16);

    auto g_xxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 17);

    auto g_xxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 18);

    auto g_xxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 19);

    auto g_xxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 20);

    auto g_xxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 21);

    auto g_xxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 22);

    auto g_xxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 23);

    auto g_xxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 24);

    auto g_xxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 25);

    auto g_xxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 26);

    auto g_xxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 27);

    auto g_xxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 28);

    auto g_xxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 29);

    auto g_xxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 30);

    auto g_xxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 31);

    auto g_xxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 32);

    auto g_xxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 33);

    auto g_xxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 34);

    auto g_xxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 35);

    auto g_xxy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk + 36);

    auto g_xxy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 37);

    auto g_xxy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 38);

    auto g_xxy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 39);

    auto g_xxy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 41);

    auto g_xxy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 42);

    auto g_xxy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 45);

    auto g_xxy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 46);

    auto g_xxy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 50);

    auto g_xxy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 51);

    auto g_xxy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 56);

    auto g_xxy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 57);

    auto g_xxy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 63);

    auto g_xxy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 64);

    auto g_xxz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk + 72);

    auto g_xxz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 73);

    auto g_xxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 74);

    auto g_xxz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 75);

    auto g_xxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 76);

    auto g_xxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 77);

    auto g_xxz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 78);

    auto g_xxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 79);

    auto g_xxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 80);

    auto g_xxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 81);

    auto g_xxz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 82);

    auto g_xxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 83);

    auto g_xxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 84);

    auto g_xxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 85);

    auto g_xxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 86);

    auto g_xxz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 87);

    auto g_xxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 88);

    auto g_xxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 89);

    auto g_xxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 90);

    auto g_xxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 91);

    auto g_xxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 92);

    auto g_xxz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 93);

    auto g_xxz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 94);

    auto g_xxz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 95);

    auto g_xxz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 96);

    auto g_xxz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 97);

    auto g_xxz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 98);

    auto g_xxz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 99);

    auto g_xxz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 101);

    auto g_xxz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 102);

    auto g_xxz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 103);

    auto g_xxz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 104);

    auto g_xxz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 105);

    auto g_xxz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 106);

    auto g_xxz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 107);

    auto g_xyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk + 108);

    auto g_xyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 109);

    auto g_xyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 111);

    auto g_xyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 112);

    auto g_xyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 114);

    auto g_xyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 115);

    auto g_xyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 116);

    auto g_xyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 118);

    auto g_xyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 119);

    auto g_xyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 120);

    auto g_xyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 121);

    auto g_xyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 123);

    auto g_xyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 124);

    auto g_xyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 125);

    auto g_xyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 126);

    auto g_xyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 127);

    auto g_xyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 129);

    auto g_xyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 130);

    auto g_xyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 131);

    auto g_xyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 132);

    auto g_xyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 133);

    auto g_xyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 134);

    auto g_xyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 136);

    auto g_xyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 137);

    auto g_xyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 138);

    auto g_xyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 139);

    auto g_xyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 140);

    auto g_xyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 141);

    auto g_xyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 142);

    auto g_xyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 143);

    auto g_xzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk + 180);

    auto g_xzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 182);

    auto g_xzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 184);

    auto g_xzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 185);

    auto g_xzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 187);

    auto g_xzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 188);

    auto g_xzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 189);

    auto g_xzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 191);

    auto g_xzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 192);

    auto g_xzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 193);

    auto g_xzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 194);

    auto g_xzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 196);

    auto g_xzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 197);

    auto g_xzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 198);

    auto g_xzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 199);

    auto g_xzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 200);

    auto g_xzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 202);

    auto g_xzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 203);

    auto g_xzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 204);

    auto g_xzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 205);

    auto g_xzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 206);

    auto g_xzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 207);

    auto g_xzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 208);

    auto g_xzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 209);

    auto g_xzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 210);

    auto g_xzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 211);

    auto g_xzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 212);

    auto g_xzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 213);

    auto g_xzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 214);

    auto g_xzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 215);

    auto g_yyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk + 216);

    auto g_yyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 217);

    auto g_yyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 218);

    auto g_yyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 219);

    auto g_yyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 220);

    auto g_yyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 221);

    auto g_yyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 222);

    auto g_yyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 223);

    auto g_yyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 224);

    auto g_yyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 225);

    auto g_yyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 226);

    auto g_yyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 227);

    auto g_yyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 228);

    auto g_yyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 229);

    auto g_yyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 230);

    auto g_yyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 231);

    auto g_yyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 232);

    auto g_yyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 233);

    auto g_yyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 234);

    auto g_yyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 235);

    auto g_yyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 236);

    auto g_yyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 237);

    auto g_yyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 238);

    auto g_yyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 239);

    auto g_yyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 240);

    auto g_yyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 241);

    auto g_yyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 242);

    auto g_yyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 243);

    auto g_yyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 244);

    auto g_yyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 245);

    auto g_yyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 246);

    auto g_yyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 247);

    auto g_yyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 248);

    auto g_yyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 249);

    auto g_yyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 250);

    auto g_yyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 251);

    auto g_yyz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 253);

    auto g_yyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 254);

    auto g_yyz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 255);

    auto g_yyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 256);

    auto g_yyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 257);

    auto g_yyz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 258);

    auto g_yyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 259);

    auto g_yyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 260);

    auto g_yyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 261);

    auto g_yyz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 262);

    auto g_yyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 263);

    auto g_yyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 264);

    auto g_yyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 265);

    auto g_yyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 266);

    auto g_yyz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 267);

    auto g_yyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 268);

    auto g_yyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 269);

    auto g_yyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 270);

    auto g_yyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 271);

    auto g_yyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 272);

    auto g_yyz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 273);

    auto g_yyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 274);

    auto g_yyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 275);

    auto g_yyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 276);

    auto g_yyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 277);

    auto g_yyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 278);

    auto g_yyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 279);

    auto g_yyz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 280);

    auto g_yyz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 281);

    auto g_yyz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 282);

    auto g_yyz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 283);

    auto g_yyz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 284);

    auto g_yyz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 285);

    auto g_yyz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 286);

    auto g_yyz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 287);

    auto g_yzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk + 288);

    auto g_yzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 289);

    auto g_yzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 290);

    auto g_yzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 291);

    auto g_yzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 292);

    auto g_yzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 293);

    auto g_yzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 294);

    auto g_yzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 295);

    auto g_yzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 296);

    auto g_yzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 297);

    auto g_yzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 298);

    auto g_yzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 299);

    auto g_yzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 300);

    auto g_yzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 301);

    auto g_yzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 302);

    auto g_yzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 303);

    auto g_yzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 304);

    auto g_yzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 305);

    auto g_yzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 306);

    auto g_yzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 307);

    auto g_yzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 308);

    auto g_yzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 309);

    auto g_yzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 310);

    auto g_yzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 311);

    auto g_yzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 312);

    auto g_yzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 313);

    auto g_yzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 314);

    auto g_yzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 315);

    auto g_yzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 316);

    auto g_yzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 317);

    auto g_yzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 318);

    auto g_yzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 319);

    auto g_yzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 320);

    auto g_yzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 321);

    auto g_yzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 322);

    auto g_yzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 323);

    auto g_zzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk + 324);

    auto g_zzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 325);

    auto g_zzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 326);

    auto g_zzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 327);

    auto g_zzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 328);

    auto g_zzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 329);

    auto g_zzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 330);

    auto g_zzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 331);

    auto g_zzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 332);

    auto g_zzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 333);

    auto g_zzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 334);

    auto g_zzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 335);

    auto g_zzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 336);

    auto g_zzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 337);

    auto g_zzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 338);

    auto g_zzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 339);

    auto g_zzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 340);

    auto g_zzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 341);

    auto g_zzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 342);

    auto g_zzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 343);

    auto g_zzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 344);

    auto g_zzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 345);

    auto g_zzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 346);

    auto g_zzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 347);

    auto g_zzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 348);

    auto g_zzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 349);

    auto g_zzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 350);

    auto g_zzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 351);

    auto g_zzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 352);

    auto g_zzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 353);

    auto g_zzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 354);

    auto g_zzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 355);

    auto g_zzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 356);

    auto g_zzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 357);

    auto g_zzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 358);

    auto g_zzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 359);

    /// Set up 0-36 components of targeted buffer : GSK

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

    #pragma omp simd aligned(g_xx_0_xxxxxxx_0, g_xx_0_xxxxxxx_1, g_xx_0_xxxxxxy_0, g_xx_0_xxxxxxy_1, g_xx_0_xxxxxxz_0, g_xx_0_xxxxxxz_1, g_xx_0_xxxxxyy_0, g_xx_0_xxxxxyy_1, g_xx_0_xxxxxyz_0, g_xx_0_xxxxxyz_1, g_xx_0_xxxxxzz_0, g_xx_0_xxxxxzz_1, g_xx_0_xxxxyyy_0, g_xx_0_xxxxyyy_1, g_xx_0_xxxxyyz_0, g_xx_0_xxxxyyz_1, g_xx_0_xxxxyzz_0, g_xx_0_xxxxyzz_1, g_xx_0_xxxxzzz_0, g_xx_0_xxxxzzz_1, g_xx_0_xxxyyyy_0, g_xx_0_xxxyyyy_1, g_xx_0_xxxyyyz_0, g_xx_0_xxxyyyz_1, g_xx_0_xxxyyzz_0, g_xx_0_xxxyyzz_1, g_xx_0_xxxyzzz_0, g_xx_0_xxxyzzz_1, g_xx_0_xxxzzzz_0, g_xx_0_xxxzzzz_1, g_xx_0_xxyyyyy_0, g_xx_0_xxyyyyy_1, g_xx_0_xxyyyyz_0, g_xx_0_xxyyyyz_1, g_xx_0_xxyyyzz_0, g_xx_0_xxyyyzz_1, g_xx_0_xxyyzzz_0, g_xx_0_xxyyzzz_1, g_xx_0_xxyzzzz_0, g_xx_0_xxyzzzz_1, g_xx_0_xxzzzzz_0, g_xx_0_xxzzzzz_1, g_xx_0_xyyyyyy_0, g_xx_0_xyyyyyy_1, g_xx_0_xyyyyyz_0, g_xx_0_xyyyyyz_1, g_xx_0_xyyyyzz_0, g_xx_0_xyyyyzz_1, g_xx_0_xyyyzzz_0, g_xx_0_xyyyzzz_1, g_xx_0_xyyzzzz_0, g_xx_0_xyyzzzz_1, g_xx_0_xyzzzzz_0, g_xx_0_xyzzzzz_1, g_xx_0_xzzzzzz_0, g_xx_0_xzzzzzz_1, g_xx_0_yyyyyyy_0, g_xx_0_yyyyyyy_1, g_xx_0_yyyyyyz_0, g_xx_0_yyyyyyz_1, g_xx_0_yyyyyzz_0, g_xx_0_yyyyyzz_1, g_xx_0_yyyyzzz_0, g_xx_0_yyyyzzz_1, g_xx_0_yyyzzzz_0, g_xx_0_yyyzzzz_1, g_xx_0_yyzzzzz_0, g_xx_0_yyzzzzz_1, g_xx_0_yzzzzzz_0, g_xx_0_yzzzzzz_1, g_xx_0_zzzzzzz_0, g_xx_0_zzzzzzz_1, g_xxx_0_xxxxxx_1, g_xxx_0_xxxxxxx_1, g_xxx_0_xxxxxxy_1, g_xxx_0_xxxxxxz_1, g_xxx_0_xxxxxy_1, g_xxx_0_xxxxxyy_1, g_xxx_0_xxxxxyz_1, g_xxx_0_xxxxxz_1, g_xxx_0_xxxxxzz_1, g_xxx_0_xxxxyy_1, g_xxx_0_xxxxyyy_1, g_xxx_0_xxxxyyz_1, g_xxx_0_xxxxyz_1, g_xxx_0_xxxxyzz_1, g_xxx_0_xxxxzz_1, g_xxx_0_xxxxzzz_1, g_xxx_0_xxxyyy_1, g_xxx_0_xxxyyyy_1, g_xxx_0_xxxyyyz_1, g_xxx_0_xxxyyz_1, g_xxx_0_xxxyyzz_1, g_xxx_0_xxxyzz_1, g_xxx_0_xxxyzzz_1, g_xxx_0_xxxzzz_1, g_xxx_0_xxxzzzz_1, g_xxx_0_xxyyyy_1, g_xxx_0_xxyyyyy_1, g_xxx_0_xxyyyyz_1, g_xxx_0_xxyyyz_1, g_xxx_0_xxyyyzz_1, g_xxx_0_xxyyzz_1, g_xxx_0_xxyyzzz_1, g_xxx_0_xxyzzz_1, g_xxx_0_xxyzzzz_1, g_xxx_0_xxzzzz_1, g_xxx_0_xxzzzzz_1, g_xxx_0_xyyyyy_1, g_xxx_0_xyyyyyy_1, g_xxx_0_xyyyyyz_1, g_xxx_0_xyyyyz_1, g_xxx_0_xyyyyzz_1, g_xxx_0_xyyyzz_1, g_xxx_0_xyyyzzz_1, g_xxx_0_xyyzzz_1, g_xxx_0_xyyzzzz_1, g_xxx_0_xyzzzz_1, g_xxx_0_xyzzzzz_1, g_xxx_0_xzzzzz_1, g_xxx_0_xzzzzzz_1, g_xxx_0_yyyyyy_1, g_xxx_0_yyyyyyy_1, g_xxx_0_yyyyyyz_1, g_xxx_0_yyyyyz_1, g_xxx_0_yyyyyzz_1, g_xxx_0_yyyyzz_1, g_xxx_0_yyyyzzz_1, g_xxx_0_yyyzzz_1, g_xxx_0_yyyzzzz_1, g_xxx_0_yyzzzz_1, g_xxx_0_yyzzzzz_1, g_xxx_0_yzzzzz_1, g_xxx_0_yzzzzzz_1, g_xxx_0_zzzzzz_1, g_xxx_0_zzzzzzz_1, g_xxxx_0_xxxxxxx_0, g_xxxx_0_xxxxxxy_0, g_xxxx_0_xxxxxxz_0, g_xxxx_0_xxxxxyy_0, g_xxxx_0_xxxxxyz_0, g_xxxx_0_xxxxxzz_0, g_xxxx_0_xxxxyyy_0, g_xxxx_0_xxxxyyz_0, g_xxxx_0_xxxxyzz_0, g_xxxx_0_xxxxzzz_0, g_xxxx_0_xxxyyyy_0, g_xxxx_0_xxxyyyz_0, g_xxxx_0_xxxyyzz_0, g_xxxx_0_xxxyzzz_0, g_xxxx_0_xxxzzzz_0, g_xxxx_0_xxyyyyy_0, g_xxxx_0_xxyyyyz_0, g_xxxx_0_xxyyyzz_0, g_xxxx_0_xxyyzzz_0, g_xxxx_0_xxyzzzz_0, g_xxxx_0_xxzzzzz_0, g_xxxx_0_xyyyyyy_0, g_xxxx_0_xyyyyyz_0, g_xxxx_0_xyyyyzz_0, g_xxxx_0_xyyyzzz_0, g_xxxx_0_xyyzzzz_0, g_xxxx_0_xyzzzzz_0, g_xxxx_0_xzzzzzz_0, g_xxxx_0_yyyyyyy_0, g_xxxx_0_yyyyyyz_0, g_xxxx_0_yyyyyzz_0, g_xxxx_0_yyyyzzz_0, g_xxxx_0_yyyzzzz_0, g_xxxx_0_yyzzzzz_0, g_xxxx_0_yzzzzzz_0, g_xxxx_0_zzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxx_0_xxxxxxx_0[i] = 3.0 * g_xx_0_xxxxxxx_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxx_1[i] * fz_be_0 + 7.0 * g_xxx_0_xxxxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxxxx_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxy_0[i] = 3.0 * g_xx_0_xxxxxxy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xxx_0_xxxxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxxxy_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxz_0[i] = 3.0 * g_xx_0_xxxxxxz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xxx_0_xxxxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxyy_0[i] = 3.0 * g_xx_0_xxxxxyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxyy_1[i] * fi_acd_0 + g_xxx_0_xxxxxyy_1[i] * wa_x[i];

        g_xxxx_0_xxxxxyz_0[i] = 3.0 * g_xx_0_xxxxxyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxyz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxzz_0[i] = 3.0 * g_xx_0_xxxxxzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxyyy_0[i] = 3.0 * g_xx_0_xxxxyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyyy_1[i] * fi_acd_0 + g_xxx_0_xxxxyyy_1[i] * wa_x[i];

        g_xxxx_0_xxxxyyz_0[i] = 3.0 * g_xx_0_xxxxyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyyz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyz_1[i] * wa_x[i];

        g_xxxx_0_xxxxyzz_0[i] = 3.0 * g_xx_0_xxxxyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxzzz_0[i] = 3.0 * g_xx_0_xxxxzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxyyyy_0[i] = 3.0 * g_xx_0_xxxyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyyy_1[i] * fi_acd_0 + g_xxx_0_xxxyyyy_1[i] * wa_x[i];

        g_xxxx_0_xxxyyyz_0[i] = 3.0 * g_xx_0_xxxyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyyz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyz_1[i] * wa_x[i];

        g_xxxx_0_xxxyyzz_0[i] = 3.0 * g_xx_0_xxxyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyzz_1[i] * wa_x[i];

        g_xxxx_0_xxxyzzz_0[i] = 3.0 * g_xx_0_xxxyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxzzzz_0[i] = 3.0 * g_xx_0_xxxzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxyyyyy_0[i] = 3.0 * g_xx_0_xxyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyyyy_1[i] * wa_x[i];

        g_xxxx_0_xxyyyyz_0[i] = 3.0 * g_xx_0_xxyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyz_1[i] * wa_x[i];

        g_xxxx_0_xxyyyzz_0[i] = 3.0 * g_xx_0_xxyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyzz_1[i] * wa_x[i];

        g_xxxx_0_xxyyzzz_0[i] = 3.0 * g_xx_0_xxyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyzzz_1[i] * wa_x[i];

        g_xxxx_0_xxyzzzz_0[i] = 3.0 * g_xx_0_xxyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxzzzzz_0[i] = 3.0 * g_xx_0_xxzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xyyyyyy_0[i] = 3.0 * g_xx_0_xyyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyyy_1[i] * fz_be_0 + g_xxx_0_yyyyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyyyy_1[i] * wa_x[i];

        g_xxxx_0_xyyyyyz_0[i] = 3.0 * g_xx_0_xyyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyyz_1[i] * fz_be_0 + g_xxx_0_yyyyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyz_1[i] * wa_x[i];

        g_xxxx_0_xyyyyzz_0[i] = 3.0 * g_xx_0_xyyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyzz_1[i] * fz_be_0 + g_xxx_0_yyyyzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyzz_1[i] * wa_x[i];

        g_xxxx_0_xyyyzzz_0[i] = 3.0 * g_xx_0_xyyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyzzz_1[i] * fz_be_0 + g_xxx_0_yyyzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyzzz_1[i] * wa_x[i];

        g_xxxx_0_xyyzzzz_0[i] = 3.0 * g_xx_0_xyyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyzzzz_1[i] * fz_be_0 + g_xxx_0_yyzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyzzzz_1[i] * wa_x[i];

        g_xxxx_0_xyzzzzz_0[i] = 3.0 * g_xx_0_xyzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyzzzzz_1[i] * fz_be_0 + g_xxx_0_yzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xzzzzzz_0[i] = 3.0 * g_xx_0_xzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xzzzzzz_1[i] * fz_be_0 + g_xxx_0_zzzzzz_1[i] * fi_acd_0 + g_xxx_0_xzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_yyyyyyy_0[i] = 3.0 * g_xx_0_yyyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyyy_1[i] * fz_be_0 + g_xxx_0_yyyyyyy_1[i] * wa_x[i];

        g_xxxx_0_yyyyyyz_0[i] = 3.0 * g_xx_0_yyyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyyz_1[i] * fz_be_0 + g_xxx_0_yyyyyyz_1[i] * wa_x[i];

        g_xxxx_0_yyyyyzz_0[i] = 3.0 * g_xx_0_yyyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyzz_1[i] * fz_be_0 + g_xxx_0_yyyyyzz_1[i] * wa_x[i];

        g_xxxx_0_yyyyzzz_0[i] = 3.0 * g_xx_0_yyyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyzzz_1[i] * fz_be_0 + g_xxx_0_yyyyzzz_1[i] * wa_x[i];

        g_xxxx_0_yyyzzzz_0[i] = 3.0 * g_xx_0_yyyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyzzzz_1[i] * fz_be_0 + g_xxx_0_yyyzzzz_1[i] * wa_x[i];

        g_xxxx_0_yyzzzzz_0[i] = 3.0 * g_xx_0_yyzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyzzzzz_1[i] * fz_be_0 + g_xxx_0_yyzzzzz_1[i] * wa_x[i];

        g_xxxx_0_yzzzzzz_0[i] = 3.0 * g_xx_0_yzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yzzzzzz_1[i] * fz_be_0 + g_xxx_0_yzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_zzzzzzz_0[i] = 3.0 * g_xx_0_zzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_zzzzzzz_1[i] * fz_be_0 + g_xxx_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 36-72 components of targeted buffer : GSK

    auto g_xxxy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 36);

    auto g_xxxy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 37);

    auto g_xxxy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 38);

    auto g_xxxy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 39);

    auto g_xxxy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 40);

    auto g_xxxy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 41);

    auto g_xxxy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 42);

    auto g_xxxy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 43);

    auto g_xxxy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 44);

    auto g_xxxy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 45);

    auto g_xxxy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 46);

    auto g_xxxy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 47);

    auto g_xxxy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 48);

    auto g_xxxy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 49);

    auto g_xxxy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 50);

    auto g_xxxy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 51);

    auto g_xxxy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 52);

    auto g_xxxy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 53);

    auto g_xxxy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 54);

    auto g_xxxy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 55);

    auto g_xxxy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 56);

    auto g_xxxy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 57);

    auto g_xxxy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 58);

    auto g_xxxy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 59);

    auto g_xxxy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 60);

    auto g_xxxy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 61);

    auto g_xxxy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 62);

    auto g_xxxy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 63);

    auto g_xxxy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 64);

    auto g_xxxy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 65);

    auto g_xxxy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 66);

    auto g_xxxy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 67);

    auto g_xxxy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 68);

    auto g_xxxy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 69);

    auto g_xxxy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 70);

    auto g_xxxy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 71);

    #pragma omp simd aligned(g_xxx_0_xxxxxx_1, g_xxx_0_xxxxxxx_1, g_xxx_0_xxxxxxy_1, g_xxx_0_xxxxxxz_1, g_xxx_0_xxxxxy_1, g_xxx_0_xxxxxyy_1, g_xxx_0_xxxxxyz_1, g_xxx_0_xxxxxz_1, g_xxx_0_xxxxxzz_1, g_xxx_0_xxxxyy_1, g_xxx_0_xxxxyyy_1, g_xxx_0_xxxxyyz_1, g_xxx_0_xxxxyz_1, g_xxx_0_xxxxyzz_1, g_xxx_0_xxxxzz_1, g_xxx_0_xxxxzzz_1, g_xxx_0_xxxyyy_1, g_xxx_0_xxxyyyy_1, g_xxx_0_xxxyyyz_1, g_xxx_0_xxxyyz_1, g_xxx_0_xxxyyzz_1, g_xxx_0_xxxyzz_1, g_xxx_0_xxxyzzz_1, g_xxx_0_xxxzzz_1, g_xxx_0_xxxzzzz_1, g_xxx_0_xxyyyy_1, g_xxx_0_xxyyyyy_1, g_xxx_0_xxyyyyz_1, g_xxx_0_xxyyyz_1, g_xxx_0_xxyyyzz_1, g_xxx_0_xxyyzz_1, g_xxx_0_xxyyzzz_1, g_xxx_0_xxyzzz_1, g_xxx_0_xxyzzzz_1, g_xxx_0_xxzzzz_1, g_xxx_0_xxzzzzz_1, g_xxx_0_xyyyyy_1, g_xxx_0_xyyyyyy_1, g_xxx_0_xyyyyyz_1, g_xxx_0_xyyyyz_1, g_xxx_0_xyyyyzz_1, g_xxx_0_xyyyzz_1, g_xxx_0_xyyyzzz_1, g_xxx_0_xyyzzz_1, g_xxx_0_xyyzzzz_1, g_xxx_0_xyzzzz_1, g_xxx_0_xyzzzzz_1, g_xxx_0_xzzzzz_1, g_xxx_0_xzzzzzz_1, g_xxx_0_yyyyyy_1, g_xxx_0_yyyyyyy_1, g_xxx_0_yyyyyyz_1, g_xxx_0_yyyyyz_1, g_xxx_0_yyyyyzz_1, g_xxx_0_yyyyzz_1, g_xxx_0_yyyyzzz_1, g_xxx_0_yyyzzz_1, g_xxx_0_yyyzzzz_1, g_xxx_0_yyzzzz_1, g_xxx_0_yyzzzzz_1, g_xxx_0_yzzzzz_1, g_xxx_0_yzzzzzz_1, g_xxx_0_zzzzzz_1, g_xxx_0_zzzzzzz_1, g_xxxy_0_xxxxxxx_0, g_xxxy_0_xxxxxxy_0, g_xxxy_0_xxxxxxz_0, g_xxxy_0_xxxxxyy_0, g_xxxy_0_xxxxxyz_0, g_xxxy_0_xxxxxzz_0, g_xxxy_0_xxxxyyy_0, g_xxxy_0_xxxxyyz_0, g_xxxy_0_xxxxyzz_0, g_xxxy_0_xxxxzzz_0, g_xxxy_0_xxxyyyy_0, g_xxxy_0_xxxyyyz_0, g_xxxy_0_xxxyyzz_0, g_xxxy_0_xxxyzzz_0, g_xxxy_0_xxxzzzz_0, g_xxxy_0_xxyyyyy_0, g_xxxy_0_xxyyyyz_0, g_xxxy_0_xxyyyzz_0, g_xxxy_0_xxyyzzz_0, g_xxxy_0_xxyzzzz_0, g_xxxy_0_xxzzzzz_0, g_xxxy_0_xyyyyyy_0, g_xxxy_0_xyyyyyz_0, g_xxxy_0_xyyyyzz_0, g_xxxy_0_xyyyzzz_0, g_xxxy_0_xyyzzzz_0, g_xxxy_0_xyzzzzz_0, g_xxxy_0_xzzzzzz_0, g_xxxy_0_yyyyyyy_0, g_xxxy_0_yyyyyyz_0, g_xxxy_0_yyyyyzz_0, g_xxxy_0_yyyyzzz_0, g_xxxy_0_yyyzzzz_0, g_xxxy_0_yyzzzzz_0, g_xxxy_0_yzzzzzz_0, g_xxxy_0_zzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxy_0_xxxxxxx_0[i] = g_xxx_0_xxxxxxx_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxy_0[i] = g_xxx_0_xxxxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxxxy_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxz_0[i] = g_xxx_0_xxxxxxz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxyy_0[i] = 2.0 * g_xxx_0_xxxxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxxyy_1[i] * wa_y[i];

        g_xxxy_0_xxxxxyz_0[i] = g_xxx_0_xxxxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxzz_0[i] = g_xxx_0_xxxxxzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxyyy_0[i] = 3.0 * g_xxx_0_xxxxyy_1[i] * fi_acd_0 + g_xxx_0_xxxxyyy_1[i] * wa_y[i];

        g_xxxy_0_xxxxyyz_0[i] = 2.0 * g_xxx_0_xxxxyz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyz_1[i] * wa_y[i];

        g_xxxy_0_xxxxyzz_0[i] = g_xxx_0_xxxxzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxzzz_0[i] = g_xxx_0_xxxxzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxyyyy_0[i] = 4.0 * g_xxx_0_xxxyyy_1[i] * fi_acd_0 + g_xxx_0_xxxyyyy_1[i] * wa_y[i];

        g_xxxy_0_xxxyyyz_0[i] = 3.0 * g_xxx_0_xxxyyz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyz_1[i] * wa_y[i];

        g_xxxy_0_xxxyyzz_0[i] = 2.0 * g_xxx_0_xxxyzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyzz_1[i] * wa_y[i];

        g_xxxy_0_xxxyzzz_0[i] = g_xxx_0_xxxzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxzzzz_0[i] = g_xxx_0_xxxzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxyyyyy_0[i] = 5.0 * g_xxx_0_xxyyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyyyy_1[i] * wa_y[i];

        g_xxxy_0_xxyyyyz_0[i] = 4.0 * g_xxx_0_xxyyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyz_1[i] * wa_y[i];

        g_xxxy_0_xxyyyzz_0[i] = 3.0 * g_xxx_0_xxyyzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyzz_1[i] * wa_y[i];

        g_xxxy_0_xxyyzzz_0[i] = 2.0 * g_xxx_0_xxyzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyzzz_1[i] * wa_y[i];

        g_xxxy_0_xxyzzzz_0[i] = g_xxx_0_xxzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxzzzzz_0[i] = g_xxx_0_xxzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xyyyyyy_0[i] = 6.0 * g_xxx_0_xyyyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyyyy_1[i] * wa_y[i];

        g_xxxy_0_xyyyyyz_0[i] = 5.0 * g_xxx_0_xyyyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyz_1[i] * wa_y[i];

        g_xxxy_0_xyyyyzz_0[i] = 4.0 * g_xxx_0_xyyyzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyzz_1[i] * wa_y[i];

        g_xxxy_0_xyyyzzz_0[i] = 3.0 * g_xxx_0_xyyzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyzzz_1[i] * wa_y[i];

        g_xxxy_0_xyyzzzz_0[i] = 2.0 * g_xxx_0_xyzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyzzzz_1[i] * wa_y[i];

        g_xxxy_0_xyzzzzz_0[i] = g_xxx_0_xzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xzzzzzz_0[i] = g_xxx_0_xzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_yyyyyyy_0[i] = 7.0 * g_xxx_0_yyyyyy_1[i] * fi_acd_0 + g_xxx_0_yyyyyyy_1[i] * wa_y[i];

        g_xxxy_0_yyyyyyz_0[i] = 6.0 * g_xxx_0_yyyyyz_1[i] * fi_acd_0 + g_xxx_0_yyyyyyz_1[i] * wa_y[i];

        g_xxxy_0_yyyyyzz_0[i] = 5.0 * g_xxx_0_yyyyzz_1[i] * fi_acd_0 + g_xxx_0_yyyyyzz_1[i] * wa_y[i];

        g_xxxy_0_yyyyzzz_0[i] = 4.0 * g_xxx_0_yyyzzz_1[i] * fi_acd_0 + g_xxx_0_yyyyzzz_1[i] * wa_y[i];

        g_xxxy_0_yyyzzzz_0[i] = 3.0 * g_xxx_0_yyzzzz_1[i] * fi_acd_0 + g_xxx_0_yyyzzzz_1[i] * wa_y[i];

        g_xxxy_0_yyzzzzz_0[i] = 2.0 * g_xxx_0_yzzzzz_1[i] * fi_acd_0 + g_xxx_0_yyzzzzz_1[i] * wa_y[i];

        g_xxxy_0_yzzzzzz_0[i] = g_xxx_0_zzzzzz_1[i] * fi_acd_0 + g_xxx_0_yzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_zzzzzzz_0[i] = g_xxx_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 72-108 components of targeted buffer : GSK

    auto g_xxxz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 72);

    auto g_xxxz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 73);

    auto g_xxxz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 74);

    auto g_xxxz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 75);

    auto g_xxxz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 76);

    auto g_xxxz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 77);

    auto g_xxxz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 78);

    auto g_xxxz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 79);

    auto g_xxxz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 80);

    auto g_xxxz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 81);

    auto g_xxxz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 82);

    auto g_xxxz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 83);

    auto g_xxxz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 84);

    auto g_xxxz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 85);

    auto g_xxxz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 86);

    auto g_xxxz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 87);

    auto g_xxxz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 88);

    auto g_xxxz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 89);

    auto g_xxxz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 90);

    auto g_xxxz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 91);

    auto g_xxxz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 92);

    auto g_xxxz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 93);

    auto g_xxxz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 94);

    auto g_xxxz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 95);

    auto g_xxxz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 96);

    auto g_xxxz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 97);

    auto g_xxxz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 98);

    auto g_xxxz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 99);

    auto g_xxxz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 100);

    auto g_xxxz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 101);

    auto g_xxxz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 102);

    auto g_xxxz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 103);

    auto g_xxxz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 104);

    auto g_xxxz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 105);

    auto g_xxxz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 106);

    auto g_xxxz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 107);

    #pragma omp simd aligned(g_xxx_0_xxxxxx_1, g_xxx_0_xxxxxxx_1, g_xxx_0_xxxxxxy_1, g_xxx_0_xxxxxxz_1, g_xxx_0_xxxxxy_1, g_xxx_0_xxxxxyy_1, g_xxx_0_xxxxxyz_1, g_xxx_0_xxxxxz_1, g_xxx_0_xxxxxzz_1, g_xxx_0_xxxxyy_1, g_xxx_0_xxxxyyy_1, g_xxx_0_xxxxyyz_1, g_xxx_0_xxxxyz_1, g_xxx_0_xxxxyzz_1, g_xxx_0_xxxxzz_1, g_xxx_0_xxxxzzz_1, g_xxx_0_xxxyyy_1, g_xxx_0_xxxyyyy_1, g_xxx_0_xxxyyyz_1, g_xxx_0_xxxyyz_1, g_xxx_0_xxxyyzz_1, g_xxx_0_xxxyzz_1, g_xxx_0_xxxyzzz_1, g_xxx_0_xxxzzz_1, g_xxx_0_xxxzzzz_1, g_xxx_0_xxyyyy_1, g_xxx_0_xxyyyyy_1, g_xxx_0_xxyyyyz_1, g_xxx_0_xxyyyz_1, g_xxx_0_xxyyyzz_1, g_xxx_0_xxyyzz_1, g_xxx_0_xxyyzzz_1, g_xxx_0_xxyzzz_1, g_xxx_0_xxyzzzz_1, g_xxx_0_xxzzzz_1, g_xxx_0_xxzzzzz_1, g_xxx_0_xyyyyy_1, g_xxx_0_xyyyyyy_1, g_xxx_0_xyyyyyz_1, g_xxx_0_xyyyyz_1, g_xxx_0_xyyyyzz_1, g_xxx_0_xyyyzz_1, g_xxx_0_xyyyzzz_1, g_xxx_0_xyyzzz_1, g_xxx_0_xyyzzzz_1, g_xxx_0_xyzzzz_1, g_xxx_0_xyzzzzz_1, g_xxx_0_xzzzzz_1, g_xxx_0_xzzzzzz_1, g_xxx_0_yyyyyy_1, g_xxx_0_yyyyyyy_1, g_xxx_0_yyyyyyz_1, g_xxx_0_yyyyyz_1, g_xxx_0_yyyyyzz_1, g_xxx_0_yyyyzz_1, g_xxx_0_yyyyzzz_1, g_xxx_0_yyyzzz_1, g_xxx_0_yyyzzzz_1, g_xxx_0_yyzzzz_1, g_xxx_0_yyzzzzz_1, g_xxx_0_yzzzzz_1, g_xxx_0_yzzzzzz_1, g_xxx_0_zzzzzz_1, g_xxx_0_zzzzzzz_1, g_xxxz_0_xxxxxxx_0, g_xxxz_0_xxxxxxy_0, g_xxxz_0_xxxxxxz_0, g_xxxz_0_xxxxxyy_0, g_xxxz_0_xxxxxyz_0, g_xxxz_0_xxxxxzz_0, g_xxxz_0_xxxxyyy_0, g_xxxz_0_xxxxyyz_0, g_xxxz_0_xxxxyzz_0, g_xxxz_0_xxxxzzz_0, g_xxxz_0_xxxyyyy_0, g_xxxz_0_xxxyyyz_0, g_xxxz_0_xxxyyzz_0, g_xxxz_0_xxxyzzz_0, g_xxxz_0_xxxzzzz_0, g_xxxz_0_xxyyyyy_0, g_xxxz_0_xxyyyyz_0, g_xxxz_0_xxyyyzz_0, g_xxxz_0_xxyyzzz_0, g_xxxz_0_xxyzzzz_0, g_xxxz_0_xxzzzzz_0, g_xxxz_0_xyyyyyy_0, g_xxxz_0_xyyyyyz_0, g_xxxz_0_xyyyyzz_0, g_xxxz_0_xyyyzzz_0, g_xxxz_0_xyyzzzz_0, g_xxxz_0_xyzzzzz_0, g_xxxz_0_xzzzzzz_0, g_xxxz_0_yyyyyyy_0, g_xxxz_0_yyyyyyz_0, g_xxxz_0_yyyyyzz_0, g_xxxz_0_yyyyzzz_0, g_xxxz_0_yyyzzzz_0, g_xxxz_0_yyzzzzz_0, g_xxxz_0_yzzzzzz_0, g_xxxz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxz_0_xxxxxxx_0[i] = g_xxx_0_xxxxxxx_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxy_0[i] = g_xxx_0_xxxxxxy_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxz_0[i] = g_xxx_0_xxxxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxxxz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxyy_0[i] = g_xxx_0_xxxxxyy_1[i] * wa_z[i];

        g_xxxz_0_xxxxxyz_0[i] = g_xxx_0_xxxxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxxyz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxzz_0[i] = 2.0 * g_xxx_0_xxxxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxxzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxyyy_0[i] = g_xxx_0_xxxxyyy_1[i] * wa_z[i];

        g_xxxz_0_xxxxyyz_0[i] = g_xxx_0_xxxxyy_1[i] * fi_acd_0 + g_xxx_0_xxxxyyz_1[i] * wa_z[i];

        g_xxxz_0_xxxxyzz_0[i] = 2.0 * g_xxx_0_xxxxyz_1[i] * fi_acd_0 + g_xxx_0_xxxxyzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxzzz_0[i] = 3.0 * g_xxx_0_xxxxzz_1[i] * fi_acd_0 + g_xxx_0_xxxxzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxyyyy_0[i] = g_xxx_0_xxxyyyy_1[i] * wa_z[i];

        g_xxxz_0_xxxyyyz_0[i] = g_xxx_0_xxxyyy_1[i] * fi_acd_0 + g_xxx_0_xxxyyyz_1[i] * wa_z[i];

        g_xxxz_0_xxxyyzz_0[i] = 2.0 * g_xxx_0_xxxyyz_1[i] * fi_acd_0 + g_xxx_0_xxxyyzz_1[i] * wa_z[i];

        g_xxxz_0_xxxyzzz_0[i] = 3.0 * g_xxx_0_xxxyzz_1[i] * fi_acd_0 + g_xxx_0_xxxyzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxzzzz_0[i] = 4.0 * g_xxx_0_xxxzzz_1[i] * fi_acd_0 + g_xxx_0_xxxzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxyyyyy_0[i] = g_xxx_0_xxyyyyy_1[i] * wa_z[i];

        g_xxxz_0_xxyyyyz_0[i] = g_xxx_0_xxyyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyyyz_1[i] * wa_z[i];

        g_xxxz_0_xxyyyzz_0[i] = 2.0 * g_xxx_0_xxyyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyyzz_1[i] * wa_z[i];

        g_xxxz_0_xxyyzzz_0[i] = 3.0 * g_xxx_0_xxyyzz_1[i] * fi_acd_0 + g_xxx_0_xxyyzzz_1[i] * wa_z[i];

        g_xxxz_0_xxyzzzz_0[i] = 4.0 * g_xxx_0_xxyzzz_1[i] * fi_acd_0 + g_xxx_0_xxyzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxzzzzz_0[i] = 5.0 * g_xxx_0_xxzzzz_1[i] * fi_acd_0 + g_xxx_0_xxzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xyyyyyy_0[i] = g_xxx_0_xyyyyyy_1[i] * wa_z[i];

        g_xxxz_0_xyyyyyz_0[i] = g_xxx_0_xyyyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyyyz_1[i] * wa_z[i];

        g_xxxz_0_xyyyyzz_0[i] = 2.0 * g_xxx_0_xyyyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyyzz_1[i] * wa_z[i];

        g_xxxz_0_xyyyzzz_0[i] = 3.0 * g_xxx_0_xyyyzz_1[i] * fi_acd_0 + g_xxx_0_xyyyzzz_1[i] * wa_z[i];

        g_xxxz_0_xyyzzzz_0[i] = 4.0 * g_xxx_0_xyyzzz_1[i] * fi_acd_0 + g_xxx_0_xyyzzzz_1[i] * wa_z[i];

        g_xxxz_0_xyzzzzz_0[i] = 5.0 * g_xxx_0_xyzzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xzzzzzz_0[i] = 6.0 * g_xxx_0_xzzzzz_1[i] * fi_acd_0 + g_xxx_0_xzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_yyyyyyy_0[i] = g_xxx_0_yyyyyyy_1[i] * wa_z[i];

        g_xxxz_0_yyyyyyz_0[i] = g_xxx_0_yyyyyy_1[i] * fi_acd_0 + g_xxx_0_yyyyyyz_1[i] * wa_z[i];

        g_xxxz_0_yyyyyzz_0[i] = 2.0 * g_xxx_0_yyyyyz_1[i] * fi_acd_0 + g_xxx_0_yyyyyzz_1[i] * wa_z[i];

        g_xxxz_0_yyyyzzz_0[i] = 3.0 * g_xxx_0_yyyyzz_1[i] * fi_acd_0 + g_xxx_0_yyyyzzz_1[i] * wa_z[i];

        g_xxxz_0_yyyzzzz_0[i] = 4.0 * g_xxx_0_yyyzzz_1[i] * fi_acd_0 + g_xxx_0_yyyzzzz_1[i] * wa_z[i];

        g_xxxz_0_yyzzzzz_0[i] = 5.0 * g_xxx_0_yyzzzz_1[i] * fi_acd_0 + g_xxx_0_yyzzzzz_1[i] * wa_z[i];

        g_xxxz_0_yzzzzzz_0[i] = 6.0 * g_xxx_0_yzzzzz_1[i] * fi_acd_0 + g_xxx_0_yzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_zzzzzzz_0[i] = 7.0 * g_xxx_0_zzzzzz_1[i] * fi_acd_0 + g_xxx_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 108-144 components of targeted buffer : GSK

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

    #pragma omp simd aligned(g_xx_0_xxxxxxx_0, g_xx_0_xxxxxxx_1, g_xx_0_xxxxxxz_0, g_xx_0_xxxxxxz_1, g_xx_0_xxxxxzz_0, g_xx_0_xxxxxzz_1, g_xx_0_xxxxzzz_0, g_xx_0_xxxxzzz_1, g_xx_0_xxxzzzz_0, g_xx_0_xxxzzzz_1, g_xx_0_xxzzzzz_0, g_xx_0_xxzzzzz_1, g_xx_0_xzzzzzz_0, g_xx_0_xzzzzzz_1, g_xxy_0_xxxxxxx_1, g_xxy_0_xxxxxxz_1, g_xxy_0_xxxxxzz_1, g_xxy_0_xxxxzzz_1, g_xxy_0_xxxzzzz_1, g_xxy_0_xxzzzzz_1, g_xxy_0_xzzzzzz_1, g_xxyy_0_xxxxxxx_0, g_xxyy_0_xxxxxxy_0, g_xxyy_0_xxxxxxz_0, g_xxyy_0_xxxxxyy_0, g_xxyy_0_xxxxxyz_0, g_xxyy_0_xxxxxzz_0, g_xxyy_0_xxxxyyy_0, g_xxyy_0_xxxxyyz_0, g_xxyy_0_xxxxyzz_0, g_xxyy_0_xxxxzzz_0, g_xxyy_0_xxxyyyy_0, g_xxyy_0_xxxyyyz_0, g_xxyy_0_xxxyyzz_0, g_xxyy_0_xxxyzzz_0, g_xxyy_0_xxxzzzz_0, g_xxyy_0_xxyyyyy_0, g_xxyy_0_xxyyyyz_0, g_xxyy_0_xxyyyzz_0, g_xxyy_0_xxyyzzz_0, g_xxyy_0_xxyzzzz_0, g_xxyy_0_xxzzzzz_0, g_xxyy_0_xyyyyyy_0, g_xxyy_0_xyyyyyz_0, g_xxyy_0_xyyyyzz_0, g_xxyy_0_xyyyzzz_0, g_xxyy_0_xyyzzzz_0, g_xxyy_0_xyzzzzz_0, g_xxyy_0_xzzzzzz_0, g_xxyy_0_yyyyyyy_0, g_xxyy_0_yyyyyyz_0, g_xxyy_0_yyyyyzz_0, g_xxyy_0_yyyyzzz_0, g_xxyy_0_yyyzzzz_0, g_xxyy_0_yyzzzzz_0, g_xxyy_0_yzzzzzz_0, g_xxyy_0_zzzzzzz_0, g_xyy_0_xxxxxxy_1, g_xyy_0_xxxxxy_1, g_xyy_0_xxxxxyy_1, g_xyy_0_xxxxxyz_1, g_xyy_0_xxxxyy_1, g_xyy_0_xxxxyyy_1, g_xyy_0_xxxxyyz_1, g_xyy_0_xxxxyz_1, g_xyy_0_xxxxyzz_1, g_xyy_0_xxxyyy_1, g_xyy_0_xxxyyyy_1, g_xyy_0_xxxyyyz_1, g_xyy_0_xxxyyz_1, g_xyy_0_xxxyyzz_1, g_xyy_0_xxxyzz_1, g_xyy_0_xxxyzzz_1, g_xyy_0_xxyyyy_1, g_xyy_0_xxyyyyy_1, g_xyy_0_xxyyyyz_1, g_xyy_0_xxyyyz_1, g_xyy_0_xxyyyzz_1, g_xyy_0_xxyyzz_1, g_xyy_0_xxyyzzz_1, g_xyy_0_xxyzzz_1, g_xyy_0_xxyzzzz_1, g_xyy_0_xyyyyy_1, g_xyy_0_xyyyyyy_1, g_xyy_0_xyyyyyz_1, g_xyy_0_xyyyyz_1, g_xyy_0_xyyyyzz_1, g_xyy_0_xyyyzz_1, g_xyy_0_xyyyzzz_1, g_xyy_0_xyyzzz_1, g_xyy_0_xyyzzzz_1, g_xyy_0_xyzzzz_1, g_xyy_0_xyzzzzz_1, g_xyy_0_yyyyyy_1, g_xyy_0_yyyyyyy_1, g_xyy_0_yyyyyyz_1, g_xyy_0_yyyyyz_1, g_xyy_0_yyyyyzz_1, g_xyy_0_yyyyzz_1, g_xyy_0_yyyyzzz_1, g_xyy_0_yyyzzz_1, g_xyy_0_yyyzzzz_1, g_xyy_0_yyzzzz_1, g_xyy_0_yyzzzzz_1, g_xyy_0_yzzzzz_1, g_xyy_0_yzzzzzz_1, g_xyy_0_zzzzzzz_1, g_yy_0_xxxxxxy_0, g_yy_0_xxxxxxy_1, g_yy_0_xxxxxyy_0, g_yy_0_xxxxxyy_1, g_yy_0_xxxxxyz_0, g_yy_0_xxxxxyz_1, g_yy_0_xxxxyyy_0, g_yy_0_xxxxyyy_1, g_yy_0_xxxxyyz_0, g_yy_0_xxxxyyz_1, g_yy_0_xxxxyzz_0, g_yy_0_xxxxyzz_1, g_yy_0_xxxyyyy_0, g_yy_0_xxxyyyy_1, g_yy_0_xxxyyyz_0, g_yy_0_xxxyyyz_1, g_yy_0_xxxyyzz_0, g_yy_0_xxxyyzz_1, g_yy_0_xxxyzzz_0, g_yy_0_xxxyzzz_1, g_yy_0_xxyyyyy_0, g_yy_0_xxyyyyy_1, g_yy_0_xxyyyyz_0, g_yy_0_xxyyyyz_1, g_yy_0_xxyyyzz_0, g_yy_0_xxyyyzz_1, g_yy_0_xxyyzzz_0, g_yy_0_xxyyzzz_1, g_yy_0_xxyzzzz_0, g_yy_0_xxyzzzz_1, g_yy_0_xyyyyyy_0, g_yy_0_xyyyyyy_1, g_yy_0_xyyyyyz_0, g_yy_0_xyyyyyz_1, g_yy_0_xyyyyzz_0, g_yy_0_xyyyyzz_1, g_yy_0_xyyyzzz_0, g_yy_0_xyyyzzz_1, g_yy_0_xyyzzzz_0, g_yy_0_xyyzzzz_1, g_yy_0_xyzzzzz_0, g_yy_0_xyzzzzz_1, g_yy_0_yyyyyyy_0, g_yy_0_yyyyyyy_1, g_yy_0_yyyyyyz_0, g_yy_0_yyyyyyz_1, g_yy_0_yyyyyzz_0, g_yy_0_yyyyyzz_1, g_yy_0_yyyyzzz_0, g_yy_0_yyyyzzz_1, g_yy_0_yyyzzzz_0, g_yy_0_yyyzzzz_1, g_yy_0_yyzzzzz_0, g_yy_0_yyzzzzz_1, g_yy_0_yzzzzzz_0, g_yy_0_yzzzzzz_1, g_yy_0_zzzzzzz_0, g_yy_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyy_0_xxxxxxx_0[i] = g_xx_0_xxxxxxx_0[i] * fbe_0 - g_xx_0_xxxxxxx_1[i] * fz_be_0 + g_xxy_0_xxxxxxx_1[i] * wa_y[i];

        g_xxyy_0_xxxxxxy_0[i] = g_yy_0_xxxxxxy_0[i] * fbe_0 - g_yy_0_xxxxxxy_1[i] * fz_be_0 + 6.0 * g_xyy_0_xxxxxy_1[i] * fi_acd_0 + g_xyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xxyy_0_xxxxxxz_0[i] = g_xx_0_xxxxxxz_0[i] * fbe_0 - g_xx_0_xxxxxxz_1[i] * fz_be_0 + g_xxy_0_xxxxxxz_1[i] * wa_y[i];

        g_xxyy_0_xxxxxyy_0[i] = g_yy_0_xxxxxyy_0[i] * fbe_0 - g_yy_0_xxxxxyy_1[i] * fz_be_0 + 5.0 * g_xyy_0_xxxxyy_1[i] * fi_acd_0 + g_xyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xxyy_0_xxxxxyz_0[i] = g_yy_0_xxxxxyz_0[i] * fbe_0 - g_yy_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xyy_0_xxxxyz_1[i] * fi_acd_0 + g_xyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xxyy_0_xxxxxzz_0[i] = g_xx_0_xxxxxzz_0[i] * fbe_0 - g_xx_0_xxxxxzz_1[i] * fz_be_0 + g_xxy_0_xxxxxzz_1[i] * wa_y[i];

        g_xxyy_0_xxxxyyy_0[i] = g_yy_0_xxxxyyy_0[i] * fbe_0 - g_yy_0_xxxxyyy_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyyy_1[i] * fi_acd_0 + g_xyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xxyy_0_xxxxyyz_0[i] = g_yy_0_xxxxyyz_0[i] * fbe_0 - g_yy_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyyz_1[i] * fi_acd_0 + g_xyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xxyy_0_xxxxyzz_0[i] = g_yy_0_xxxxyzz_0[i] * fbe_0 - g_yy_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyzz_1[i] * fi_acd_0 + g_xyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xxyy_0_xxxxzzz_0[i] = g_xx_0_xxxxzzz_0[i] * fbe_0 - g_xx_0_xxxxzzz_1[i] * fz_be_0 + g_xxy_0_xxxxzzz_1[i] * wa_y[i];

        g_xxyy_0_xxxyyyy_0[i] = g_yy_0_xxxyyyy_0[i] * fbe_0 - g_yy_0_xxxyyyy_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyyy_1[i] * fi_acd_0 + g_xyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xxyy_0_xxxyyyz_0[i] = g_yy_0_xxxyyyz_0[i] * fbe_0 - g_yy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyyz_1[i] * fi_acd_0 + g_xyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xxyy_0_xxxyyzz_0[i] = g_yy_0_xxxyyzz_0[i] * fbe_0 - g_yy_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyzz_1[i] * fi_acd_0 + g_xyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xxyy_0_xxxyzzz_0[i] = g_yy_0_xxxyzzz_0[i] * fbe_0 - g_yy_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyzzz_1[i] * fi_acd_0 + g_xyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xxyy_0_xxxzzzz_0[i] = g_xx_0_xxxzzzz_0[i] * fbe_0 - g_xx_0_xxxzzzz_1[i] * fz_be_0 + g_xxy_0_xxxzzzz_1[i] * wa_y[i];

        g_xxyy_0_xxyyyyy_0[i] = g_yy_0_xxyyyyy_0[i] * fbe_0 - g_yy_0_xxyyyyy_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyyy_1[i] * fi_acd_0 + g_xyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xxyy_0_xxyyyyz_0[i] = g_yy_0_xxyyyyz_0[i] * fbe_0 - g_yy_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyyz_1[i] * fi_acd_0 + g_xyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xxyy_0_xxyyyzz_0[i] = g_yy_0_xxyyyzz_0[i] * fbe_0 - g_yy_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyzz_1[i] * fi_acd_0 + g_xyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xxyy_0_xxyyzzz_0[i] = g_yy_0_xxyyzzz_0[i] * fbe_0 - g_yy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyzzz_1[i] * fi_acd_0 + g_xyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xxyy_0_xxyzzzz_0[i] = g_yy_0_xxyzzzz_0[i] * fbe_0 - g_yy_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyzzzz_1[i] * fi_acd_0 + g_xyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xxyy_0_xxzzzzz_0[i] = g_xx_0_xxzzzzz_0[i] * fbe_0 - g_xx_0_xxzzzzz_1[i] * fz_be_0 + g_xxy_0_xxzzzzz_1[i] * wa_y[i];

        g_xxyy_0_xyyyyyy_0[i] = g_yy_0_xyyyyyy_0[i] * fbe_0 - g_yy_0_xyyyyyy_1[i] * fz_be_0 + g_xyy_0_yyyyyy_1[i] * fi_acd_0 + g_xyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xxyy_0_xyyyyyz_0[i] = g_yy_0_xyyyyyz_0[i] * fbe_0 - g_yy_0_xyyyyyz_1[i] * fz_be_0 + g_xyy_0_yyyyyz_1[i] * fi_acd_0 + g_xyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xxyy_0_xyyyyzz_0[i] = g_yy_0_xyyyyzz_0[i] * fbe_0 - g_yy_0_xyyyyzz_1[i] * fz_be_0 + g_xyy_0_yyyyzz_1[i] * fi_acd_0 + g_xyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xxyy_0_xyyyzzz_0[i] = g_yy_0_xyyyzzz_0[i] * fbe_0 - g_yy_0_xyyyzzz_1[i] * fz_be_0 + g_xyy_0_yyyzzz_1[i] * fi_acd_0 + g_xyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xxyy_0_xyyzzzz_0[i] = g_yy_0_xyyzzzz_0[i] * fbe_0 - g_yy_0_xyyzzzz_1[i] * fz_be_0 + g_xyy_0_yyzzzz_1[i] * fi_acd_0 + g_xyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xxyy_0_xyzzzzz_0[i] = g_yy_0_xyzzzzz_0[i] * fbe_0 - g_yy_0_xyzzzzz_1[i] * fz_be_0 + g_xyy_0_yzzzzz_1[i] * fi_acd_0 + g_xyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xxyy_0_xzzzzzz_0[i] = g_xx_0_xzzzzzz_0[i] * fbe_0 - g_xx_0_xzzzzzz_1[i] * fz_be_0 + g_xxy_0_xzzzzzz_1[i] * wa_y[i];

        g_xxyy_0_yyyyyyy_0[i] = g_yy_0_yyyyyyy_0[i] * fbe_0 - g_yy_0_yyyyyyy_1[i] * fz_be_0 + g_xyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xxyy_0_yyyyyyz_0[i] = g_yy_0_yyyyyyz_0[i] * fbe_0 - g_yy_0_yyyyyyz_1[i] * fz_be_0 + g_xyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xxyy_0_yyyyyzz_0[i] = g_yy_0_yyyyyzz_0[i] * fbe_0 - g_yy_0_yyyyyzz_1[i] * fz_be_0 + g_xyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xxyy_0_yyyyzzz_0[i] = g_yy_0_yyyyzzz_0[i] * fbe_0 - g_yy_0_yyyyzzz_1[i] * fz_be_0 + g_xyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xxyy_0_yyyzzzz_0[i] = g_yy_0_yyyzzzz_0[i] * fbe_0 - g_yy_0_yyyzzzz_1[i] * fz_be_0 + g_xyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xxyy_0_yyzzzzz_0[i] = g_yy_0_yyzzzzz_0[i] * fbe_0 - g_yy_0_yyzzzzz_1[i] * fz_be_0 + g_xyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xxyy_0_yzzzzzz_0[i] = g_yy_0_yzzzzzz_0[i] * fbe_0 - g_yy_0_yzzzzzz_1[i] * fz_be_0 + g_xyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xxyy_0_zzzzzzz_0[i] = g_yy_0_zzzzzzz_0[i] * fbe_0 - g_yy_0_zzzzzzz_1[i] * fz_be_0 + g_xyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 144-180 components of targeted buffer : GSK

    auto g_xxyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 144);

    auto g_xxyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 145);

    auto g_xxyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 146);

    auto g_xxyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 147);

    auto g_xxyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 148);

    auto g_xxyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 149);

    auto g_xxyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 150);

    auto g_xxyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 151);

    auto g_xxyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 152);

    auto g_xxyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 153);

    auto g_xxyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 154);

    auto g_xxyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 155);

    auto g_xxyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 156);

    auto g_xxyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 157);

    auto g_xxyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 158);

    auto g_xxyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 159);

    auto g_xxyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 160);

    auto g_xxyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 161);

    auto g_xxyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 162);

    auto g_xxyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 163);

    auto g_xxyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 164);

    auto g_xxyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 165);

    auto g_xxyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 166);

    auto g_xxyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 167);

    auto g_xxyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 168);

    auto g_xxyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 169);

    auto g_xxyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 170);

    auto g_xxyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 171);

    auto g_xxyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 172);

    auto g_xxyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 173);

    auto g_xxyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 174);

    auto g_xxyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 175);

    auto g_xxyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 176);

    auto g_xxyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 177);

    auto g_xxyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 178);

    auto g_xxyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 179);

    #pragma omp simd aligned(g_xxy_0_xxxxxxy_1, g_xxy_0_xxxxxyy_1, g_xxy_0_xxxxyyy_1, g_xxy_0_xxxyyyy_1, g_xxy_0_xxyyyyy_1, g_xxy_0_xyyyyyy_1, g_xxy_0_yyyyyyy_1, g_xxyz_0_xxxxxxx_0, g_xxyz_0_xxxxxxy_0, g_xxyz_0_xxxxxxz_0, g_xxyz_0_xxxxxyy_0, g_xxyz_0_xxxxxyz_0, g_xxyz_0_xxxxxzz_0, g_xxyz_0_xxxxyyy_0, g_xxyz_0_xxxxyyz_0, g_xxyz_0_xxxxyzz_0, g_xxyz_0_xxxxzzz_0, g_xxyz_0_xxxyyyy_0, g_xxyz_0_xxxyyyz_0, g_xxyz_0_xxxyyzz_0, g_xxyz_0_xxxyzzz_0, g_xxyz_0_xxxzzzz_0, g_xxyz_0_xxyyyyy_0, g_xxyz_0_xxyyyyz_0, g_xxyz_0_xxyyyzz_0, g_xxyz_0_xxyyzzz_0, g_xxyz_0_xxyzzzz_0, g_xxyz_0_xxzzzzz_0, g_xxyz_0_xyyyyyy_0, g_xxyz_0_xyyyyyz_0, g_xxyz_0_xyyyyzz_0, g_xxyz_0_xyyyzzz_0, g_xxyz_0_xyyzzzz_0, g_xxyz_0_xyzzzzz_0, g_xxyz_0_xzzzzzz_0, g_xxyz_0_yyyyyyy_0, g_xxyz_0_yyyyyyz_0, g_xxyz_0_yyyyyzz_0, g_xxyz_0_yyyyzzz_0, g_xxyz_0_yyyzzzz_0, g_xxyz_0_yyzzzzz_0, g_xxyz_0_yzzzzzz_0, g_xxyz_0_zzzzzzz_0, g_xxz_0_xxxxxxx_1, g_xxz_0_xxxxxxz_1, g_xxz_0_xxxxxyz_1, g_xxz_0_xxxxxz_1, g_xxz_0_xxxxxzz_1, g_xxz_0_xxxxyyz_1, g_xxz_0_xxxxyz_1, g_xxz_0_xxxxyzz_1, g_xxz_0_xxxxzz_1, g_xxz_0_xxxxzzz_1, g_xxz_0_xxxyyyz_1, g_xxz_0_xxxyyz_1, g_xxz_0_xxxyyzz_1, g_xxz_0_xxxyzz_1, g_xxz_0_xxxyzzz_1, g_xxz_0_xxxzzz_1, g_xxz_0_xxxzzzz_1, g_xxz_0_xxyyyyz_1, g_xxz_0_xxyyyz_1, g_xxz_0_xxyyyzz_1, g_xxz_0_xxyyzz_1, g_xxz_0_xxyyzzz_1, g_xxz_0_xxyzzz_1, g_xxz_0_xxyzzzz_1, g_xxz_0_xxzzzz_1, g_xxz_0_xxzzzzz_1, g_xxz_0_xyyyyyz_1, g_xxz_0_xyyyyz_1, g_xxz_0_xyyyyzz_1, g_xxz_0_xyyyzz_1, g_xxz_0_xyyyzzz_1, g_xxz_0_xyyzzz_1, g_xxz_0_xyyzzzz_1, g_xxz_0_xyzzzz_1, g_xxz_0_xyzzzzz_1, g_xxz_0_xzzzzz_1, g_xxz_0_xzzzzzz_1, g_xxz_0_yyyyyyz_1, g_xxz_0_yyyyyz_1, g_xxz_0_yyyyyzz_1, g_xxz_0_yyyyzz_1, g_xxz_0_yyyyzzz_1, g_xxz_0_yyyzzz_1, g_xxz_0_yyyzzzz_1, g_xxz_0_yyzzzz_1, g_xxz_0_yyzzzzz_1, g_xxz_0_yzzzzz_1, g_xxz_0_yzzzzzz_1, g_xxz_0_zzzzzz_1, g_xxz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyz_0_xxxxxxx_0[i] = g_xxz_0_xxxxxxx_1[i] * wa_y[i];

        g_xxyz_0_xxxxxxy_0[i] = g_xxy_0_xxxxxxy_1[i] * wa_z[i];

        g_xxyz_0_xxxxxxz_0[i] = g_xxz_0_xxxxxxz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxyy_0[i] = g_xxy_0_xxxxxyy_1[i] * wa_z[i];

        g_xxyz_0_xxxxxyz_0[i] = g_xxz_0_xxxxxz_1[i] * fi_acd_0 + g_xxz_0_xxxxxyz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxzz_0[i] = g_xxz_0_xxxxxzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxyyy_0[i] = g_xxy_0_xxxxyyy_1[i] * wa_z[i];

        g_xxyz_0_xxxxyyz_0[i] = 2.0 * g_xxz_0_xxxxyz_1[i] * fi_acd_0 + g_xxz_0_xxxxyyz_1[i] * wa_y[i];

        g_xxyz_0_xxxxyzz_0[i] = g_xxz_0_xxxxzz_1[i] * fi_acd_0 + g_xxz_0_xxxxyzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxzzz_0[i] = g_xxz_0_xxxxzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxyyyy_0[i] = g_xxy_0_xxxyyyy_1[i] * wa_z[i];

        g_xxyz_0_xxxyyyz_0[i] = 3.0 * g_xxz_0_xxxyyz_1[i] * fi_acd_0 + g_xxz_0_xxxyyyz_1[i] * wa_y[i];

        g_xxyz_0_xxxyyzz_0[i] = 2.0 * g_xxz_0_xxxyzz_1[i] * fi_acd_0 + g_xxz_0_xxxyyzz_1[i] * wa_y[i];

        g_xxyz_0_xxxyzzz_0[i] = g_xxz_0_xxxzzz_1[i] * fi_acd_0 + g_xxz_0_xxxyzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxzzzz_0[i] = g_xxz_0_xxxzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxyyyyy_0[i] = g_xxy_0_xxyyyyy_1[i] * wa_z[i];

        g_xxyz_0_xxyyyyz_0[i] = 4.0 * g_xxz_0_xxyyyz_1[i] * fi_acd_0 + g_xxz_0_xxyyyyz_1[i] * wa_y[i];

        g_xxyz_0_xxyyyzz_0[i] = 3.0 * g_xxz_0_xxyyzz_1[i] * fi_acd_0 + g_xxz_0_xxyyyzz_1[i] * wa_y[i];

        g_xxyz_0_xxyyzzz_0[i] = 2.0 * g_xxz_0_xxyzzz_1[i] * fi_acd_0 + g_xxz_0_xxyyzzz_1[i] * wa_y[i];

        g_xxyz_0_xxyzzzz_0[i] = g_xxz_0_xxzzzz_1[i] * fi_acd_0 + g_xxz_0_xxyzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxzzzzz_0[i] = g_xxz_0_xxzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xyyyyyy_0[i] = g_xxy_0_xyyyyyy_1[i] * wa_z[i];

        g_xxyz_0_xyyyyyz_0[i] = 5.0 * g_xxz_0_xyyyyz_1[i] * fi_acd_0 + g_xxz_0_xyyyyyz_1[i] * wa_y[i];

        g_xxyz_0_xyyyyzz_0[i] = 4.0 * g_xxz_0_xyyyzz_1[i] * fi_acd_0 + g_xxz_0_xyyyyzz_1[i] * wa_y[i];

        g_xxyz_0_xyyyzzz_0[i] = 3.0 * g_xxz_0_xyyzzz_1[i] * fi_acd_0 + g_xxz_0_xyyyzzz_1[i] * wa_y[i];

        g_xxyz_0_xyyzzzz_0[i] = 2.0 * g_xxz_0_xyzzzz_1[i] * fi_acd_0 + g_xxz_0_xyyzzzz_1[i] * wa_y[i];

        g_xxyz_0_xyzzzzz_0[i] = g_xxz_0_xzzzzz_1[i] * fi_acd_0 + g_xxz_0_xyzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xzzzzzz_0[i] = g_xxz_0_xzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_yyyyyyy_0[i] = g_xxy_0_yyyyyyy_1[i] * wa_z[i];

        g_xxyz_0_yyyyyyz_0[i] = 6.0 * g_xxz_0_yyyyyz_1[i] * fi_acd_0 + g_xxz_0_yyyyyyz_1[i] * wa_y[i];

        g_xxyz_0_yyyyyzz_0[i] = 5.0 * g_xxz_0_yyyyzz_1[i] * fi_acd_0 + g_xxz_0_yyyyyzz_1[i] * wa_y[i];

        g_xxyz_0_yyyyzzz_0[i] = 4.0 * g_xxz_0_yyyzzz_1[i] * fi_acd_0 + g_xxz_0_yyyyzzz_1[i] * wa_y[i];

        g_xxyz_0_yyyzzzz_0[i] = 3.0 * g_xxz_0_yyzzzz_1[i] * fi_acd_0 + g_xxz_0_yyyzzzz_1[i] * wa_y[i];

        g_xxyz_0_yyzzzzz_0[i] = 2.0 * g_xxz_0_yzzzzz_1[i] * fi_acd_0 + g_xxz_0_yyzzzzz_1[i] * wa_y[i];

        g_xxyz_0_yzzzzzz_0[i] = g_xxz_0_zzzzzz_1[i] * fi_acd_0 + g_xxz_0_yzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_zzzzzzz_0[i] = g_xxz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 180-216 components of targeted buffer : GSK

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

    #pragma omp simd aligned(g_xx_0_xxxxxxx_0, g_xx_0_xxxxxxx_1, g_xx_0_xxxxxxy_0, g_xx_0_xxxxxxy_1, g_xx_0_xxxxxyy_0, g_xx_0_xxxxxyy_1, g_xx_0_xxxxyyy_0, g_xx_0_xxxxyyy_1, g_xx_0_xxxyyyy_0, g_xx_0_xxxyyyy_1, g_xx_0_xxyyyyy_0, g_xx_0_xxyyyyy_1, g_xx_0_xyyyyyy_0, g_xx_0_xyyyyyy_1, g_xxz_0_xxxxxxx_1, g_xxz_0_xxxxxxy_1, g_xxz_0_xxxxxyy_1, g_xxz_0_xxxxyyy_1, g_xxz_0_xxxyyyy_1, g_xxz_0_xxyyyyy_1, g_xxz_0_xyyyyyy_1, g_xxzz_0_xxxxxxx_0, g_xxzz_0_xxxxxxy_0, g_xxzz_0_xxxxxxz_0, g_xxzz_0_xxxxxyy_0, g_xxzz_0_xxxxxyz_0, g_xxzz_0_xxxxxzz_0, g_xxzz_0_xxxxyyy_0, g_xxzz_0_xxxxyyz_0, g_xxzz_0_xxxxyzz_0, g_xxzz_0_xxxxzzz_0, g_xxzz_0_xxxyyyy_0, g_xxzz_0_xxxyyyz_0, g_xxzz_0_xxxyyzz_0, g_xxzz_0_xxxyzzz_0, g_xxzz_0_xxxzzzz_0, g_xxzz_0_xxyyyyy_0, g_xxzz_0_xxyyyyz_0, g_xxzz_0_xxyyyzz_0, g_xxzz_0_xxyyzzz_0, g_xxzz_0_xxyzzzz_0, g_xxzz_0_xxzzzzz_0, g_xxzz_0_xyyyyyy_0, g_xxzz_0_xyyyyyz_0, g_xxzz_0_xyyyyzz_0, g_xxzz_0_xyyyzzz_0, g_xxzz_0_xyyzzzz_0, g_xxzz_0_xyzzzzz_0, g_xxzz_0_xzzzzzz_0, g_xxzz_0_yyyyyyy_0, g_xxzz_0_yyyyyyz_0, g_xxzz_0_yyyyyzz_0, g_xxzz_0_yyyyzzz_0, g_xxzz_0_yyyzzzz_0, g_xxzz_0_yyzzzzz_0, g_xxzz_0_yzzzzzz_0, g_xxzz_0_zzzzzzz_0, g_xzz_0_xxxxxxz_1, g_xzz_0_xxxxxyz_1, g_xzz_0_xxxxxz_1, g_xzz_0_xxxxxzz_1, g_xzz_0_xxxxyyz_1, g_xzz_0_xxxxyz_1, g_xzz_0_xxxxyzz_1, g_xzz_0_xxxxzz_1, g_xzz_0_xxxxzzz_1, g_xzz_0_xxxyyyz_1, g_xzz_0_xxxyyz_1, g_xzz_0_xxxyyzz_1, g_xzz_0_xxxyzz_1, g_xzz_0_xxxyzzz_1, g_xzz_0_xxxzzz_1, g_xzz_0_xxxzzzz_1, g_xzz_0_xxyyyyz_1, g_xzz_0_xxyyyz_1, g_xzz_0_xxyyyzz_1, g_xzz_0_xxyyzz_1, g_xzz_0_xxyyzzz_1, g_xzz_0_xxyzzz_1, g_xzz_0_xxyzzzz_1, g_xzz_0_xxzzzz_1, g_xzz_0_xxzzzzz_1, g_xzz_0_xyyyyyz_1, g_xzz_0_xyyyyz_1, g_xzz_0_xyyyyzz_1, g_xzz_0_xyyyzz_1, g_xzz_0_xyyyzzz_1, g_xzz_0_xyyzzz_1, g_xzz_0_xyyzzzz_1, g_xzz_0_xyzzzz_1, g_xzz_0_xyzzzzz_1, g_xzz_0_xzzzzz_1, g_xzz_0_xzzzzzz_1, g_xzz_0_yyyyyyy_1, g_xzz_0_yyyyyyz_1, g_xzz_0_yyyyyz_1, g_xzz_0_yyyyyzz_1, g_xzz_0_yyyyzz_1, g_xzz_0_yyyyzzz_1, g_xzz_0_yyyzzz_1, g_xzz_0_yyyzzzz_1, g_xzz_0_yyzzzz_1, g_xzz_0_yyzzzzz_1, g_xzz_0_yzzzzz_1, g_xzz_0_yzzzzzz_1, g_xzz_0_zzzzzz_1, g_xzz_0_zzzzzzz_1, g_zz_0_xxxxxxz_0, g_zz_0_xxxxxxz_1, g_zz_0_xxxxxyz_0, g_zz_0_xxxxxyz_1, g_zz_0_xxxxxzz_0, g_zz_0_xxxxxzz_1, g_zz_0_xxxxyyz_0, g_zz_0_xxxxyyz_1, g_zz_0_xxxxyzz_0, g_zz_0_xxxxyzz_1, g_zz_0_xxxxzzz_0, g_zz_0_xxxxzzz_1, g_zz_0_xxxyyyz_0, g_zz_0_xxxyyyz_1, g_zz_0_xxxyyzz_0, g_zz_0_xxxyyzz_1, g_zz_0_xxxyzzz_0, g_zz_0_xxxyzzz_1, g_zz_0_xxxzzzz_0, g_zz_0_xxxzzzz_1, g_zz_0_xxyyyyz_0, g_zz_0_xxyyyyz_1, g_zz_0_xxyyyzz_0, g_zz_0_xxyyyzz_1, g_zz_0_xxyyzzz_0, g_zz_0_xxyyzzz_1, g_zz_0_xxyzzzz_0, g_zz_0_xxyzzzz_1, g_zz_0_xxzzzzz_0, g_zz_0_xxzzzzz_1, g_zz_0_xyyyyyz_0, g_zz_0_xyyyyyz_1, g_zz_0_xyyyyzz_0, g_zz_0_xyyyyzz_1, g_zz_0_xyyyzzz_0, g_zz_0_xyyyzzz_1, g_zz_0_xyyzzzz_0, g_zz_0_xyyzzzz_1, g_zz_0_xyzzzzz_0, g_zz_0_xyzzzzz_1, g_zz_0_xzzzzzz_0, g_zz_0_xzzzzzz_1, g_zz_0_yyyyyyy_0, g_zz_0_yyyyyyy_1, g_zz_0_yyyyyyz_0, g_zz_0_yyyyyyz_1, g_zz_0_yyyyyzz_0, g_zz_0_yyyyyzz_1, g_zz_0_yyyyzzz_0, g_zz_0_yyyyzzz_1, g_zz_0_yyyzzzz_0, g_zz_0_yyyzzzz_1, g_zz_0_yyzzzzz_0, g_zz_0_yyzzzzz_1, g_zz_0_yzzzzzz_0, g_zz_0_yzzzzzz_1, g_zz_0_zzzzzzz_0, g_zz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzz_0_xxxxxxx_0[i] = g_xx_0_xxxxxxx_0[i] * fbe_0 - g_xx_0_xxxxxxx_1[i] * fz_be_0 + g_xxz_0_xxxxxxx_1[i] * wa_z[i];

        g_xxzz_0_xxxxxxy_0[i] = g_xx_0_xxxxxxy_0[i] * fbe_0 - g_xx_0_xxxxxxy_1[i] * fz_be_0 + g_xxz_0_xxxxxxy_1[i] * wa_z[i];

        g_xxzz_0_xxxxxxz_0[i] = g_zz_0_xxxxxxz_0[i] * fbe_0 - g_zz_0_xxxxxxz_1[i] * fz_be_0 + 6.0 * g_xzz_0_xxxxxz_1[i] * fi_acd_0 + g_xzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxyy_0[i] = g_xx_0_xxxxxyy_0[i] * fbe_0 - g_xx_0_xxxxxyy_1[i] * fz_be_0 + g_xxz_0_xxxxxyy_1[i] * wa_z[i];

        g_xxzz_0_xxxxxyz_0[i] = g_zz_0_xxxxxyz_0[i] * fbe_0 - g_zz_0_xxxxxyz_1[i] * fz_be_0 + 5.0 * g_xzz_0_xxxxyz_1[i] * fi_acd_0 + g_xzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxzz_0[i] = g_zz_0_xxxxxzz_0[i] * fbe_0 - g_zz_0_xxxxxzz_1[i] * fz_be_0 + 5.0 * g_xzz_0_xxxxzz_1[i] * fi_acd_0 + g_xzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxyyy_0[i] = g_xx_0_xxxxyyy_0[i] * fbe_0 - g_xx_0_xxxxyyy_1[i] * fz_be_0 + g_xxz_0_xxxxyyy_1[i] * wa_z[i];

        g_xxzz_0_xxxxyyz_0[i] = g_zz_0_xxxxyyz_0[i] * fbe_0 - g_zz_0_xxxxyyz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxyyz_1[i] * fi_acd_0 + g_xzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xxzz_0_xxxxyzz_0[i] = g_zz_0_xxxxyzz_0[i] * fbe_0 - g_zz_0_xxxxyzz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxyzz_1[i] * fi_acd_0 + g_xzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxzzz_0[i] = g_zz_0_xxxxzzz_0[i] * fbe_0 - g_zz_0_xxxxzzz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxzzz_1[i] * fi_acd_0 + g_xzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxyyyy_0[i] = g_xx_0_xxxyyyy_0[i] * fbe_0 - g_xx_0_xxxyyyy_1[i] * fz_be_0 + g_xxz_0_xxxyyyy_1[i] * wa_z[i];

        g_xxzz_0_xxxyyyz_0[i] = g_zz_0_xxxyyyz_0[i] * fbe_0 - g_zz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyyyz_1[i] * fi_acd_0 + g_xzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xxzz_0_xxxyyzz_0[i] = g_zz_0_xxxyyzz_0[i] * fbe_0 - g_zz_0_xxxyyzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyyzz_1[i] * fi_acd_0 + g_xzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xxzz_0_xxxyzzz_0[i] = g_zz_0_xxxyzzz_0[i] * fbe_0 - g_zz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyzzz_1[i] * fi_acd_0 + g_xzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxzzzz_0[i] = g_zz_0_xxxzzzz_0[i] * fbe_0 - g_zz_0_xxxzzzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxzzzz_1[i] * fi_acd_0 + g_xzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxyyyyy_0[i] = g_xx_0_xxyyyyy_0[i] * fbe_0 - g_xx_0_xxyyyyy_1[i] * fz_be_0 + g_xxz_0_xxyyyyy_1[i] * wa_z[i];

        g_xxzz_0_xxyyyyz_0[i] = g_zz_0_xxyyyyz_0[i] * fbe_0 - g_zz_0_xxyyyyz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyyyz_1[i] * fi_acd_0 + g_xzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xxzz_0_xxyyyzz_0[i] = g_zz_0_xxyyyzz_0[i] * fbe_0 - g_zz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyyzz_1[i] * fi_acd_0 + g_xzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xxzz_0_xxyyzzz_0[i] = g_zz_0_xxyyzzz_0[i] * fbe_0 - g_zz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyzzz_1[i] * fi_acd_0 + g_xzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xxzz_0_xxyzzzz_0[i] = g_zz_0_xxyzzzz_0[i] * fbe_0 - g_zz_0_xxyzzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyzzzz_1[i] * fi_acd_0 + g_xzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxzzzzz_0[i] = g_zz_0_xxzzzzz_0[i] * fbe_0 - g_zz_0_xxzzzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xzzzzz_1[i] * fi_acd_0 + g_xzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xyyyyyy_0[i] = g_xx_0_xyyyyyy_0[i] * fbe_0 - g_xx_0_xyyyyyy_1[i] * fz_be_0 + g_xxz_0_xyyyyyy_1[i] * wa_z[i];

        g_xxzz_0_xyyyyyz_0[i] = g_zz_0_xyyyyyz_0[i] * fbe_0 - g_zz_0_xyyyyyz_1[i] * fz_be_0 + g_xzz_0_yyyyyz_1[i] * fi_acd_0 + g_xzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xxzz_0_xyyyyzz_0[i] = g_zz_0_xyyyyzz_0[i] * fbe_0 - g_zz_0_xyyyyzz_1[i] * fz_be_0 + g_xzz_0_yyyyzz_1[i] * fi_acd_0 + g_xzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xxzz_0_xyyyzzz_0[i] = g_zz_0_xyyyzzz_0[i] * fbe_0 - g_zz_0_xyyyzzz_1[i] * fz_be_0 + g_xzz_0_yyyzzz_1[i] * fi_acd_0 + g_xzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xxzz_0_xyyzzzz_0[i] = g_zz_0_xyyzzzz_0[i] * fbe_0 - g_zz_0_xyyzzzz_1[i] * fz_be_0 + g_xzz_0_yyzzzz_1[i] * fi_acd_0 + g_xzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xxzz_0_xyzzzzz_0[i] = g_zz_0_xyzzzzz_0[i] * fbe_0 - g_zz_0_xyzzzzz_1[i] * fz_be_0 + g_xzz_0_yzzzzz_1[i] * fi_acd_0 + g_xzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xzzzzzz_0[i] = g_zz_0_xzzzzzz_0[i] * fbe_0 - g_zz_0_xzzzzzz_1[i] * fz_be_0 + g_xzz_0_zzzzzz_1[i] * fi_acd_0 + g_xzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_yyyyyyy_0[i] = g_zz_0_yyyyyyy_0[i] * fbe_0 - g_zz_0_yyyyyyy_1[i] * fz_be_0 + g_xzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xxzz_0_yyyyyyz_0[i] = g_zz_0_yyyyyyz_0[i] * fbe_0 - g_zz_0_yyyyyyz_1[i] * fz_be_0 + g_xzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xxzz_0_yyyyyzz_0[i] = g_zz_0_yyyyyzz_0[i] * fbe_0 - g_zz_0_yyyyyzz_1[i] * fz_be_0 + g_xzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xxzz_0_yyyyzzz_0[i] = g_zz_0_yyyyzzz_0[i] * fbe_0 - g_zz_0_yyyyzzz_1[i] * fz_be_0 + g_xzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xxzz_0_yyyzzzz_0[i] = g_zz_0_yyyzzzz_0[i] * fbe_0 - g_zz_0_yyyzzzz_1[i] * fz_be_0 + g_xzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xxzz_0_yyzzzzz_0[i] = g_zz_0_yyzzzzz_0[i] * fbe_0 - g_zz_0_yyzzzzz_1[i] * fz_be_0 + g_xzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xxzz_0_yzzzzzz_0[i] = g_zz_0_yzzzzzz_0[i] * fbe_0 - g_zz_0_yzzzzzz_1[i] * fz_be_0 + g_xzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_zzzzzzz_0[i] = g_zz_0_zzzzzzz_0[i] * fbe_0 - g_zz_0_zzzzzzz_1[i] * fz_be_0 + g_xzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 216-252 components of targeted buffer : GSK

    auto g_xyyy_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 216);

    auto g_xyyy_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 217);

    auto g_xyyy_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 218);

    auto g_xyyy_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 219);

    auto g_xyyy_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 220);

    auto g_xyyy_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 221);

    auto g_xyyy_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 222);

    auto g_xyyy_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 223);

    auto g_xyyy_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 224);

    auto g_xyyy_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 225);

    auto g_xyyy_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 226);

    auto g_xyyy_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 227);

    auto g_xyyy_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 228);

    auto g_xyyy_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 229);

    auto g_xyyy_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 230);

    auto g_xyyy_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 231);

    auto g_xyyy_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 232);

    auto g_xyyy_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 233);

    auto g_xyyy_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 234);

    auto g_xyyy_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 235);

    auto g_xyyy_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 236);

    auto g_xyyy_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 237);

    auto g_xyyy_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 238);

    auto g_xyyy_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 239);

    auto g_xyyy_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 240);

    auto g_xyyy_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 241);

    auto g_xyyy_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 242);

    auto g_xyyy_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 243);

    auto g_xyyy_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 244);

    auto g_xyyy_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 245);

    auto g_xyyy_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 246);

    auto g_xyyy_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 247);

    auto g_xyyy_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 248);

    auto g_xyyy_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 249);

    auto g_xyyy_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 250);

    auto g_xyyy_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 251);

    #pragma omp simd aligned(g_xyyy_0_xxxxxxx_0, g_xyyy_0_xxxxxxy_0, g_xyyy_0_xxxxxxz_0, g_xyyy_0_xxxxxyy_0, g_xyyy_0_xxxxxyz_0, g_xyyy_0_xxxxxzz_0, g_xyyy_0_xxxxyyy_0, g_xyyy_0_xxxxyyz_0, g_xyyy_0_xxxxyzz_0, g_xyyy_0_xxxxzzz_0, g_xyyy_0_xxxyyyy_0, g_xyyy_0_xxxyyyz_0, g_xyyy_0_xxxyyzz_0, g_xyyy_0_xxxyzzz_0, g_xyyy_0_xxxzzzz_0, g_xyyy_0_xxyyyyy_0, g_xyyy_0_xxyyyyz_0, g_xyyy_0_xxyyyzz_0, g_xyyy_0_xxyyzzz_0, g_xyyy_0_xxyzzzz_0, g_xyyy_0_xxzzzzz_0, g_xyyy_0_xyyyyyy_0, g_xyyy_0_xyyyyyz_0, g_xyyy_0_xyyyyzz_0, g_xyyy_0_xyyyzzz_0, g_xyyy_0_xyyzzzz_0, g_xyyy_0_xyzzzzz_0, g_xyyy_0_xzzzzzz_0, g_xyyy_0_yyyyyyy_0, g_xyyy_0_yyyyyyz_0, g_xyyy_0_yyyyyzz_0, g_xyyy_0_yyyyzzz_0, g_xyyy_0_yyyzzzz_0, g_xyyy_0_yyzzzzz_0, g_xyyy_0_yzzzzzz_0, g_xyyy_0_zzzzzzz_0, g_yyy_0_xxxxxx_1, g_yyy_0_xxxxxxx_1, g_yyy_0_xxxxxxy_1, g_yyy_0_xxxxxxz_1, g_yyy_0_xxxxxy_1, g_yyy_0_xxxxxyy_1, g_yyy_0_xxxxxyz_1, g_yyy_0_xxxxxz_1, g_yyy_0_xxxxxzz_1, g_yyy_0_xxxxyy_1, g_yyy_0_xxxxyyy_1, g_yyy_0_xxxxyyz_1, g_yyy_0_xxxxyz_1, g_yyy_0_xxxxyzz_1, g_yyy_0_xxxxzz_1, g_yyy_0_xxxxzzz_1, g_yyy_0_xxxyyy_1, g_yyy_0_xxxyyyy_1, g_yyy_0_xxxyyyz_1, g_yyy_0_xxxyyz_1, g_yyy_0_xxxyyzz_1, g_yyy_0_xxxyzz_1, g_yyy_0_xxxyzzz_1, g_yyy_0_xxxzzz_1, g_yyy_0_xxxzzzz_1, g_yyy_0_xxyyyy_1, g_yyy_0_xxyyyyy_1, g_yyy_0_xxyyyyz_1, g_yyy_0_xxyyyz_1, g_yyy_0_xxyyyzz_1, g_yyy_0_xxyyzz_1, g_yyy_0_xxyyzzz_1, g_yyy_0_xxyzzz_1, g_yyy_0_xxyzzzz_1, g_yyy_0_xxzzzz_1, g_yyy_0_xxzzzzz_1, g_yyy_0_xyyyyy_1, g_yyy_0_xyyyyyy_1, g_yyy_0_xyyyyyz_1, g_yyy_0_xyyyyz_1, g_yyy_0_xyyyyzz_1, g_yyy_0_xyyyzz_1, g_yyy_0_xyyyzzz_1, g_yyy_0_xyyzzz_1, g_yyy_0_xyyzzzz_1, g_yyy_0_xyzzzz_1, g_yyy_0_xyzzzzz_1, g_yyy_0_xzzzzz_1, g_yyy_0_xzzzzzz_1, g_yyy_0_yyyyyy_1, g_yyy_0_yyyyyyy_1, g_yyy_0_yyyyyyz_1, g_yyy_0_yyyyyz_1, g_yyy_0_yyyyyzz_1, g_yyy_0_yyyyzz_1, g_yyy_0_yyyyzzz_1, g_yyy_0_yyyzzz_1, g_yyy_0_yyyzzzz_1, g_yyy_0_yyzzzz_1, g_yyy_0_yyzzzzz_1, g_yyy_0_yzzzzz_1, g_yyy_0_yzzzzzz_1, g_yyy_0_zzzzzz_1, g_yyy_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyy_0_xxxxxxx_0[i] = 7.0 * g_yyy_0_xxxxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxxxx_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxy_0[i] = 6.0 * g_yyy_0_xxxxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxxxy_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxz_0[i] = 6.0 * g_yyy_0_xxxxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxyy_0[i] = 5.0 * g_yyy_0_xxxxyy_1[i] * fi_acd_0 + g_yyy_0_xxxxxyy_1[i] * wa_x[i];

        g_xyyy_0_xxxxxyz_0[i] = 5.0 * g_yyy_0_xxxxyz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxzz_0[i] = 5.0 * g_yyy_0_xxxxzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxyyy_0[i] = 4.0 * g_yyy_0_xxxyyy_1[i] * fi_acd_0 + g_yyy_0_xxxxyyy_1[i] * wa_x[i];

        g_xyyy_0_xxxxyyz_0[i] = 4.0 * g_yyy_0_xxxyyz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyy_0_xxxxyzz_0[i] = 4.0 * g_yyy_0_xxxyzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxzzz_0[i] = 4.0 * g_yyy_0_xxxzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxyyyy_0[i] = 3.0 * g_yyy_0_xxyyyy_1[i] * fi_acd_0 + g_yyy_0_xxxyyyy_1[i] * wa_x[i];

        g_xyyy_0_xxxyyyz_0[i] = 3.0 * g_yyy_0_xxyyyz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyy_0_xxxyyzz_0[i] = 3.0 * g_yyy_0_xxyyzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyy_0_xxxyzzz_0[i] = 3.0 * g_yyy_0_xxyzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxzzzz_0[i] = 3.0 * g_yyy_0_xxzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxyyyyy_0[i] = 2.0 * g_yyy_0_xyyyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyyyy_1[i] * wa_x[i];

        g_xyyy_0_xxyyyyz_0[i] = 2.0 * g_yyy_0_xyyyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyy_0_xxyyyzz_0[i] = 2.0 * g_yyy_0_xyyyzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyy_0_xxyyzzz_0[i] = 2.0 * g_yyy_0_xyyzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyy_0_xxyzzzz_0[i] = 2.0 * g_yyy_0_xyzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxzzzzz_0[i] = 2.0 * g_yyy_0_xzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xyyyyyy_0[i] = g_yyy_0_yyyyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyyyy_1[i] * wa_x[i];

        g_xyyy_0_xyyyyyz_0[i] = g_yyy_0_yyyyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyy_0_xyyyyzz_0[i] = g_yyy_0_yyyyzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyy_0_xyyyzzz_0[i] = g_yyy_0_yyyzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyy_0_xyyzzzz_0[i] = g_yyy_0_yyzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyy_0_xyzzzzz_0[i] = g_yyy_0_yzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xzzzzzz_0[i] = g_yyy_0_zzzzzz_1[i] * fi_acd_0 + g_yyy_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_yyyyyyy_0[i] = g_yyy_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyy_0_yyyyyyz_0[i] = g_yyy_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyy_0_yyyyyzz_0[i] = g_yyy_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyy_0_yyyyzzz_0[i] = g_yyy_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyy_0_yyyzzzz_0[i] = g_yyy_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyy_0_yyzzzzz_0[i] = g_yyy_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyy_0_yzzzzzz_0[i] = g_yyy_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_zzzzzzz_0[i] = g_yyy_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 252-288 components of targeted buffer : GSK

    auto g_xyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 252);

    auto g_xyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 253);

    auto g_xyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 254);

    auto g_xyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 255);

    auto g_xyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 256);

    auto g_xyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 257);

    auto g_xyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 258);

    auto g_xyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 259);

    auto g_xyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 260);

    auto g_xyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 261);

    auto g_xyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 262);

    auto g_xyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 263);

    auto g_xyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 264);

    auto g_xyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 265);

    auto g_xyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 266);

    auto g_xyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 267);

    auto g_xyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 268);

    auto g_xyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 269);

    auto g_xyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 270);

    auto g_xyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 271);

    auto g_xyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 272);

    auto g_xyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 273);

    auto g_xyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 274);

    auto g_xyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 275);

    auto g_xyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 276);

    auto g_xyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 277);

    auto g_xyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 278);

    auto g_xyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 279);

    auto g_xyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 280);

    auto g_xyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 281);

    auto g_xyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 282);

    auto g_xyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 283);

    auto g_xyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 284);

    auto g_xyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 285);

    auto g_xyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 286);

    auto g_xyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 287);

    #pragma omp simd aligned(g_xyy_0_xxxxxxx_1, g_xyy_0_xxxxxxy_1, g_xyy_0_xxxxxyy_1, g_xyy_0_xxxxyyy_1, g_xyy_0_xxxyyyy_1, g_xyy_0_xxyyyyy_1, g_xyy_0_xyyyyyy_1, g_xyyz_0_xxxxxxx_0, g_xyyz_0_xxxxxxy_0, g_xyyz_0_xxxxxxz_0, g_xyyz_0_xxxxxyy_0, g_xyyz_0_xxxxxyz_0, g_xyyz_0_xxxxxzz_0, g_xyyz_0_xxxxyyy_0, g_xyyz_0_xxxxyyz_0, g_xyyz_0_xxxxyzz_0, g_xyyz_0_xxxxzzz_0, g_xyyz_0_xxxyyyy_0, g_xyyz_0_xxxyyyz_0, g_xyyz_0_xxxyyzz_0, g_xyyz_0_xxxyzzz_0, g_xyyz_0_xxxzzzz_0, g_xyyz_0_xxyyyyy_0, g_xyyz_0_xxyyyyz_0, g_xyyz_0_xxyyyzz_0, g_xyyz_0_xxyyzzz_0, g_xyyz_0_xxyzzzz_0, g_xyyz_0_xxzzzzz_0, g_xyyz_0_xyyyyyy_0, g_xyyz_0_xyyyyyz_0, g_xyyz_0_xyyyyzz_0, g_xyyz_0_xyyyzzz_0, g_xyyz_0_xyyzzzz_0, g_xyyz_0_xyzzzzz_0, g_xyyz_0_xzzzzzz_0, g_xyyz_0_yyyyyyy_0, g_xyyz_0_yyyyyyz_0, g_xyyz_0_yyyyyzz_0, g_xyyz_0_yyyyzzz_0, g_xyyz_0_yyyzzzz_0, g_xyyz_0_yyzzzzz_0, g_xyyz_0_yzzzzzz_0, g_xyyz_0_zzzzzzz_0, g_yyz_0_xxxxxxz_1, g_yyz_0_xxxxxyz_1, g_yyz_0_xxxxxz_1, g_yyz_0_xxxxxzz_1, g_yyz_0_xxxxyyz_1, g_yyz_0_xxxxyz_1, g_yyz_0_xxxxyzz_1, g_yyz_0_xxxxzz_1, g_yyz_0_xxxxzzz_1, g_yyz_0_xxxyyyz_1, g_yyz_0_xxxyyz_1, g_yyz_0_xxxyyzz_1, g_yyz_0_xxxyzz_1, g_yyz_0_xxxyzzz_1, g_yyz_0_xxxzzz_1, g_yyz_0_xxxzzzz_1, g_yyz_0_xxyyyyz_1, g_yyz_0_xxyyyz_1, g_yyz_0_xxyyyzz_1, g_yyz_0_xxyyzz_1, g_yyz_0_xxyyzzz_1, g_yyz_0_xxyzzz_1, g_yyz_0_xxyzzzz_1, g_yyz_0_xxzzzz_1, g_yyz_0_xxzzzzz_1, g_yyz_0_xyyyyyz_1, g_yyz_0_xyyyyz_1, g_yyz_0_xyyyyzz_1, g_yyz_0_xyyyzz_1, g_yyz_0_xyyyzzz_1, g_yyz_0_xyyzzz_1, g_yyz_0_xyyzzzz_1, g_yyz_0_xyzzzz_1, g_yyz_0_xyzzzzz_1, g_yyz_0_xzzzzz_1, g_yyz_0_xzzzzzz_1, g_yyz_0_yyyyyyy_1, g_yyz_0_yyyyyyz_1, g_yyz_0_yyyyyz_1, g_yyz_0_yyyyyzz_1, g_yyz_0_yyyyzz_1, g_yyz_0_yyyyzzz_1, g_yyz_0_yyyzzz_1, g_yyz_0_yyyzzzz_1, g_yyz_0_yyzzzz_1, g_yyz_0_yyzzzzz_1, g_yyz_0_yzzzzz_1, g_yyz_0_yzzzzzz_1, g_yyz_0_zzzzzz_1, g_yyz_0_zzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyz_0_xxxxxxx_0[i] = g_xyy_0_xxxxxxx_1[i] * wa_z[i];

        g_xyyz_0_xxxxxxy_0[i] = g_xyy_0_xxxxxxy_1[i] * wa_z[i];

        g_xyyz_0_xxxxxxz_0[i] = 6.0 * g_yyz_0_xxxxxz_1[i] * fi_acd_0 + g_yyz_0_xxxxxxz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxyy_0[i] = g_xyy_0_xxxxxyy_1[i] * wa_z[i];

        g_xyyz_0_xxxxxyz_0[i] = 5.0 * g_yyz_0_xxxxyz_1[i] * fi_acd_0 + g_yyz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxzz_0[i] = 5.0 * g_yyz_0_xxxxzz_1[i] * fi_acd_0 + g_yyz_0_xxxxxzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxyyy_0[i] = g_xyy_0_xxxxyyy_1[i] * wa_z[i];

        g_xyyz_0_xxxxyyz_0[i] = 4.0 * g_yyz_0_xxxyyz_1[i] * fi_acd_0 + g_yyz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyyz_0_xxxxyzz_0[i] = 4.0 * g_yyz_0_xxxyzz_1[i] * fi_acd_0 + g_yyz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxzzz_0[i] = 4.0 * g_yyz_0_xxxzzz_1[i] * fi_acd_0 + g_yyz_0_xxxxzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxyyyy_0[i] = g_xyy_0_xxxyyyy_1[i] * wa_z[i];

        g_xyyz_0_xxxyyyz_0[i] = 3.0 * g_yyz_0_xxyyyz_1[i] * fi_acd_0 + g_yyz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyyz_0_xxxyyzz_0[i] = 3.0 * g_yyz_0_xxyyzz_1[i] * fi_acd_0 + g_yyz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyyz_0_xxxyzzz_0[i] = 3.0 * g_yyz_0_xxyzzz_1[i] * fi_acd_0 + g_yyz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxzzzz_0[i] = 3.0 * g_yyz_0_xxzzzz_1[i] * fi_acd_0 + g_yyz_0_xxxzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxyyyyy_0[i] = g_xyy_0_xxyyyyy_1[i] * wa_z[i];

        g_xyyz_0_xxyyyyz_0[i] = 2.0 * g_yyz_0_xyyyyz_1[i] * fi_acd_0 + g_yyz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyyz_0_xxyyyzz_0[i] = 2.0 * g_yyz_0_xyyyzz_1[i] * fi_acd_0 + g_yyz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyyz_0_xxyyzzz_0[i] = 2.0 * g_yyz_0_xyyzzz_1[i] * fi_acd_0 + g_yyz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyyz_0_xxyzzzz_0[i] = 2.0 * g_yyz_0_xyzzzz_1[i] * fi_acd_0 + g_yyz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxzzzzz_0[i] = 2.0 * g_yyz_0_xzzzzz_1[i] * fi_acd_0 + g_yyz_0_xxzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xyyyyyy_0[i] = g_xyy_0_xyyyyyy_1[i] * wa_z[i];

        g_xyyz_0_xyyyyyz_0[i] = g_yyz_0_yyyyyz_1[i] * fi_acd_0 + g_yyz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyyz_0_xyyyyzz_0[i] = g_yyz_0_yyyyzz_1[i] * fi_acd_0 + g_yyz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyyz_0_xyyyzzz_0[i] = g_yyz_0_yyyzzz_1[i] * fi_acd_0 + g_yyz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyyz_0_xyyzzzz_0[i] = g_yyz_0_yyzzzz_1[i] * fi_acd_0 + g_yyz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyyz_0_xyzzzzz_0[i] = g_yyz_0_yzzzzz_1[i] * fi_acd_0 + g_yyz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xzzzzzz_0[i] = g_yyz_0_zzzzzz_1[i] * fi_acd_0 + g_yyz_0_xzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_yyyyyyy_0[i] = g_yyz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyyz_0_yyyyyyz_0[i] = g_yyz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyyz_0_yyyyyzz_0[i] = g_yyz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyyz_0_yyyyzzz_0[i] = g_yyz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyyz_0_yyyzzzz_0[i] = g_yyz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyyz_0_yyzzzzz_0[i] = g_yyz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyyz_0_yzzzzzz_0[i] = g_yyz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_zzzzzzz_0[i] = g_yyz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 288-324 components of targeted buffer : GSK

    auto g_xyzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 288);

    auto g_xyzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 289);

    auto g_xyzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 290);

    auto g_xyzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 291);

    auto g_xyzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 292);

    auto g_xyzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 293);

    auto g_xyzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 294);

    auto g_xyzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 295);

    auto g_xyzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 296);

    auto g_xyzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 297);

    auto g_xyzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 298);

    auto g_xyzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 299);

    auto g_xyzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 300);

    auto g_xyzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 301);

    auto g_xyzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 302);

    auto g_xyzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 303);

    auto g_xyzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 304);

    auto g_xyzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 305);

    auto g_xyzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 306);

    auto g_xyzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 307);

    auto g_xyzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 308);

    auto g_xyzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 309);

    auto g_xyzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 310);

    auto g_xyzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 311);

    auto g_xyzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 312);

    auto g_xyzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 313);

    auto g_xyzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 314);

    auto g_xyzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 315);

    auto g_xyzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 316);

    auto g_xyzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 317);

    auto g_xyzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 318);

    auto g_xyzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 319);

    auto g_xyzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 320);

    auto g_xyzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 321);

    auto g_xyzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 322);

    auto g_xyzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 323);

    #pragma omp simd aligned(g_xyzz_0_xxxxxxx_0, g_xyzz_0_xxxxxxy_0, g_xyzz_0_xxxxxxz_0, g_xyzz_0_xxxxxyy_0, g_xyzz_0_xxxxxyz_0, g_xyzz_0_xxxxxzz_0, g_xyzz_0_xxxxyyy_0, g_xyzz_0_xxxxyyz_0, g_xyzz_0_xxxxyzz_0, g_xyzz_0_xxxxzzz_0, g_xyzz_0_xxxyyyy_0, g_xyzz_0_xxxyyyz_0, g_xyzz_0_xxxyyzz_0, g_xyzz_0_xxxyzzz_0, g_xyzz_0_xxxzzzz_0, g_xyzz_0_xxyyyyy_0, g_xyzz_0_xxyyyyz_0, g_xyzz_0_xxyyyzz_0, g_xyzz_0_xxyyzzz_0, g_xyzz_0_xxyzzzz_0, g_xyzz_0_xxzzzzz_0, g_xyzz_0_xyyyyyy_0, g_xyzz_0_xyyyyyz_0, g_xyzz_0_xyyyyzz_0, g_xyzz_0_xyyyzzz_0, g_xyzz_0_xyyzzzz_0, g_xyzz_0_xyzzzzz_0, g_xyzz_0_xzzzzzz_0, g_xyzz_0_yyyyyyy_0, g_xyzz_0_yyyyyyz_0, g_xyzz_0_yyyyyzz_0, g_xyzz_0_yyyyzzz_0, g_xyzz_0_yyyzzzz_0, g_xyzz_0_yyzzzzz_0, g_xyzz_0_yzzzzzz_0, g_xyzz_0_zzzzzzz_0, g_xzz_0_xxxxxxx_1, g_xzz_0_xxxxxxz_1, g_xzz_0_xxxxxzz_1, g_xzz_0_xxxxzzz_1, g_xzz_0_xxxzzzz_1, g_xzz_0_xxzzzzz_1, g_xzz_0_xzzzzzz_1, g_yzz_0_xxxxxxy_1, g_yzz_0_xxxxxy_1, g_yzz_0_xxxxxyy_1, g_yzz_0_xxxxxyz_1, g_yzz_0_xxxxyy_1, g_yzz_0_xxxxyyy_1, g_yzz_0_xxxxyyz_1, g_yzz_0_xxxxyz_1, g_yzz_0_xxxxyzz_1, g_yzz_0_xxxyyy_1, g_yzz_0_xxxyyyy_1, g_yzz_0_xxxyyyz_1, g_yzz_0_xxxyyz_1, g_yzz_0_xxxyyzz_1, g_yzz_0_xxxyzz_1, g_yzz_0_xxxyzzz_1, g_yzz_0_xxyyyy_1, g_yzz_0_xxyyyyy_1, g_yzz_0_xxyyyyz_1, g_yzz_0_xxyyyz_1, g_yzz_0_xxyyyzz_1, g_yzz_0_xxyyzz_1, g_yzz_0_xxyyzzz_1, g_yzz_0_xxyzzz_1, g_yzz_0_xxyzzzz_1, g_yzz_0_xyyyyy_1, g_yzz_0_xyyyyyy_1, g_yzz_0_xyyyyyz_1, g_yzz_0_xyyyyz_1, g_yzz_0_xyyyyzz_1, g_yzz_0_xyyyzz_1, g_yzz_0_xyyyzzz_1, g_yzz_0_xyyzzz_1, g_yzz_0_xyyzzzz_1, g_yzz_0_xyzzzz_1, g_yzz_0_xyzzzzz_1, g_yzz_0_yyyyyy_1, g_yzz_0_yyyyyyy_1, g_yzz_0_yyyyyyz_1, g_yzz_0_yyyyyz_1, g_yzz_0_yyyyyzz_1, g_yzz_0_yyyyzz_1, g_yzz_0_yyyyzzz_1, g_yzz_0_yyyzzz_1, g_yzz_0_yyyzzzz_1, g_yzz_0_yyzzzz_1, g_yzz_0_yyzzzzz_1, g_yzz_0_yzzzzz_1, g_yzz_0_yzzzzzz_1, g_yzz_0_zzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzz_0_xxxxxxx_0[i] = g_xzz_0_xxxxxxx_1[i] * wa_y[i];

        g_xyzz_0_xxxxxxy_0[i] = 6.0 * g_yzz_0_xxxxxy_1[i] * fi_acd_0 + g_yzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xyzz_0_xxxxxxz_0[i] = g_xzz_0_xxxxxxz_1[i] * wa_y[i];

        g_xyzz_0_xxxxxyy_0[i] = 5.0 * g_yzz_0_xxxxyy_1[i] * fi_acd_0 + g_yzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xyzz_0_xxxxxyz_0[i] = 5.0 * g_yzz_0_xxxxyz_1[i] * fi_acd_0 + g_yzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xyzz_0_xxxxxzz_0[i] = g_xzz_0_xxxxxzz_1[i] * wa_y[i];

        g_xyzz_0_xxxxyyy_0[i] = 4.0 * g_yzz_0_xxxyyy_1[i] * fi_acd_0 + g_yzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xyzz_0_xxxxyyz_0[i] = 4.0 * g_yzz_0_xxxyyz_1[i] * fi_acd_0 + g_yzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xyzz_0_xxxxyzz_0[i] = 4.0 * g_yzz_0_xxxyzz_1[i] * fi_acd_0 + g_yzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xyzz_0_xxxxzzz_0[i] = g_xzz_0_xxxxzzz_1[i] * wa_y[i];

        g_xyzz_0_xxxyyyy_0[i] = 3.0 * g_yzz_0_xxyyyy_1[i] * fi_acd_0 + g_yzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xyzz_0_xxxyyyz_0[i] = 3.0 * g_yzz_0_xxyyyz_1[i] * fi_acd_0 + g_yzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xyzz_0_xxxyyzz_0[i] = 3.0 * g_yzz_0_xxyyzz_1[i] * fi_acd_0 + g_yzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xyzz_0_xxxyzzz_0[i] = 3.0 * g_yzz_0_xxyzzz_1[i] * fi_acd_0 + g_yzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xyzz_0_xxxzzzz_0[i] = g_xzz_0_xxxzzzz_1[i] * wa_y[i];

        g_xyzz_0_xxyyyyy_0[i] = 2.0 * g_yzz_0_xyyyyy_1[i] * fi_acd_0 + g_yzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xyzz_0_xxyyyyz_0[i] = 2.0 * g_yzz_0_xyyyyz_1[i] * fi_acd_0 + g_yzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xyzz_0_xxyyyzz_0[i] = 2.0 * g_yzz_0_xyyyzz_1[i] * fi_acd_0 + g_yzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xyzz_0_xxyyzzz_0[i] = 2.0 * g_yzz_0_xyyzzz_1[i] * fi_acd_0 + g_yzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xyzz_0_xxyzzzz_0[i] = 2.0 * g_yzz_0_xyzzzz_1[i] * fi_acd_0 + g_yzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xyzz_0_xxzzzzz_0[i] = g_xzz_0_xxzzzzz_1[i] * wa_y[i];

        g_xyzz_0_xyyyyyy_0[i] = g_yzz_0_yyyyyy_1[i] * fi_acd_0 + g_yzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xyzz_0_xyyyyyz_0[i] = g_yzz_0_yyyyyz_1[i] * fi_acd_0 + g_yzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xyzz_0_xyyyyzz_0[i] = g_yzz_0_yyyyzz_1[i] * fi_acd_0 + g_yzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xyzz_0_xyyyzzz_0[i] = g_yzz_0_yyyzzz_1[i] * fi_acd_0 + g_yzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xyzz_0_xyyzzzz_0[i] = g_yzz_0_yyzzzz_1[i] * fi_acd_0 + g_yzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xyzz_0_xyzzzzz_0[i] = g_yzz_0_yzzzzz_1[i] * fi_acd_0 + g_yzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xyzz_0_xzzzzzz_0[i] = g_xzz_0_xzzzzzz_1[i] * wa_y[i];

        g_xyzz_0_yyyyyyy_0[i] = g_yzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xyzz_0_yyyyyyz_0[i] = g_yzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xyzz_0_yyyyyzz_0[i] = g_yzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xyzz_0_yyyyzzz_0[i] = g_yzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xyzz_0_yyyzzzz_0[i] = g_yzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xyzz_0_yyzzzzz_0[i] = g_yzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xyzz_0_yzzzzzz_0[i] = g_yzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xyzz_0_zzzzzzz_0[i] = g_yzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 324-360 components of targeted buffer : GSK

    auto g_xzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 324);

    auto g_xzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 325);

    auto g_xzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 326);

    auto g_xzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 327);

    auto g_xzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 328);

    auto g_xzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 329);

    auto g_xzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 330);

    auto g_xzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 331);

    auto g_xzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 332);

    auto g_xzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 333);

    auto g_xzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 334);

    auto g_xzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 335);

    auto g_xzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 336);

    auto g_xzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 337);

    auto g_xzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 338);

    auto g_xzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 339);

    auto g_xzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 340);

    auto g_xzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 341);

    auto g_xzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 342);

    auto g_xzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 343);

    auto g_xzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 344);

    auto g_xzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 345);

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

    #pragma omp simd aligned(g_xzzz_0_xxxxxxx_0, g_xzzz_0_xxxxxxy_0, g_xzzz_0_xxxxxxz_0, g_xzzz_0_xxxxxyy_0, g_xzzz_0_xxxxxyz_0, g_xzzz_0_xxxxxzz_0, g_xzzz_0_xxxxyyy_0, g_xzzz_0_xxxxyyz_0, g_xzzz_0_xxxxyzz_0, g_xzzz_0_xxxxzzz_0, g_xzzz_0_xxxyyyy_0, g_xzzz_0_xxxyyyz_0, g_xzzz_0_xxxyyzz_0, g_xzzz_0_xxxyzzz_0, g_xzzz_0_xxxzzzz_0, g_xzzz_0_xxyyyyy_0, g_xzzz_0_xxyyyyz_0, g_xzzz_0_xxyyyzz_0, g_xzzz_0_xxyyzzz_0, g_xzzz_0_xxyzzzz_0, g_xzzz_0_xxzzzzz_0, g_xzzz_0_xyyyyyy_0, g_xzzz_0_xyyyyyz_0, g_xzzz_0_xyyyyzz_0, g_xzzz_0_xyyyzzz_0, g_xzzz_0_xyyzzzz_0, g_xzzz_0_xyzzzzz_0, g_xzzz_0_xzzzzzz_0, g_xzzz_0_yyyyyyy_0, g_xzzz_0_yyyyyyz_0, g_xzzz_0_yyyyyzz_0, g_xzzz_0_yyyyzzz_0, g_xzzz_0_yyyzzzz_0, g_xzzz_0_yyzzzzz_0, g_xzzz_0_yzzzzzz_0, g_xzzz_0_zzzzzzz_0, g_zzz_0_xxxxxx_1, g_zzz_0_xxxxxxx_1, g_zzz_0_xxxxxxy_1, g_zzz_0_xxxxxxz_1, g_zzz_0_xxxxxy_1, g_zzz_0_xxxxxyy_1, g_zzz_0_xxxxxyz_1, g_zzz_0_xxxxxz_1, g_zzz_0_xxxxxzz_1, g_zzz_0_xxxxyy_1, g_zzz_0_xxxxyyy_1, g_zzz_0_xxxxyyz_1, g_zzz_0_xxxxyz_1, g_zzz_0_xxxxyzz_1, g_zzz_0_xxxxzz_1, g_zzz_0_xxxxzzz_1, g_zzz_0_xxxyyy_1, g_zzz_0_xxxyyyy_1, g_zzz_0_xxxyyyz_1, g_zzz_0_xxxyyz_1, g_zzz_0_xxxyyzz_1, g_zzz_0_xxxyzz_1, g_zzz_0_xxxyzzz_1, g_zzz_0_xxxzzz_1, g_zzz_0_xxxzzzz_1, g_zzz_0_xxyyyy_1, g_zzz_0_xxyyyyy_1, g_zzz_0_xxyyyyz_1, g_zzz_0_xxyyyz_1, g_zzz_0_xxyyyzz_1, g_zzz_0_xxyyzz_1, g_zzz_0_xxyyzzz_1, g_zzz_0_xxyzzz_1, g_zzz_0_xxyzzzz_1, g_zzz_0_xxzzzz_1, g_zzz_0_xxzzzzz_1, g_zzz_0_xyyyyy_1, g_zzz_0_xyyyyyy_1, g_zzz_0_xyyyyyz_1, g_zzz_0_xyyyyz_1, g_zzz_0_xyyyyzz_1, g_zzz_0_xyyyzz_1, g_zzz_0_xyyyzzz_1, g_zzz_0_xyyzzz_1, g_zzz_0_xyyzzzz_1, g_zzz_0_xyzzzz_1, g_zzz_0_xyzzzzz_1, g_zzz_0_xzzzzz_1, g_zzz_0_xzzzzzz_1, g_zzz_0_yyyyyy_1, g_zzz_0_yyyyyyy_1, g_zzz_0_yyyyyyz_1, g_zzz_0_yyyyyz_1, g_zzz_0_yyyyyzz_1, g_zzz_0_yyyyzz_1, g_zzz_0_yyyyzzz_1, g_zzz_0_yyyzzz_1, g_zzz_0_yyyzzzz_1, g_zzz_0_yyzzzz_1, g_zzz_0_yyzzzzz_1, g_zzz_0_yzzzzz_1, g_zzz_0_yzzzzzz_1, g_zzz_0_zzzzzz_1, g_zzz_0_zzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzz_0_xxxxxxx_0[i] = 7.0 * g_zzz_0_xxxxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxxxx_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxy_0[i] = 6.0 * g_zzz_0_xxxxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxxxy_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxz_0[i] = 6.0 * g_zzz_0_xxxxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxyy_0[i] = 5.0 * g_zzz_0_xxxxyy_1[i] * fi_acd_0 + g_zzz_0_xxxxxyy_1[i] * wa_x[i];

        g_xzzz_0_xxxxxyz_0[i] = 5.0 * g_zzz_0_xxxxyz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxzz_0[i] = 5.0 * g_zzz_0_xxxxzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxyyy_0[i] = 4.0 * g_zzz_0_xxxyyy_1[i] * fi_acd_0 + g_zzz_0_xxxxyyy_1[i] * wa_x[i];

        g_xzzz_0_xxxxyyz_0[i] = 4.0 * g_zzz_0_xxxyyz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyz_1[i] * wa_x[i];

        g_xzzz_0_xxxxyzz_0[i] = 4.0 * g_zzz_0_xxxyzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxzzz_0[i] = 4.0 * g_zzz_0_xxxzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxyyyy_0[i] = 3.0 * g_zzz_0_xxyyyy_1[i] * fi_acd_0 + g_zzz_0_xxxyyyy_1[i] * wa_x[i];

        g_xzzz_0_xxxyyyz_0[i] = 3.0 * g_zzz_0_xxyyyz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyz_1[i] * wa_x[i];

        g_xzzz_0_xxxyyzz_0[i] = 3.0 * g_zzz_0_xxyyzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyzz_1[i] * wa_x[i];

        g_xzzz_0_xxxyzzz_0[i] = 3.0 * g_zzz_0_xxyzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxzzzz_0[i] = 3.0 * g_zzz_0_xxzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxyyyyy_0[i] = 2.0 * g_zzz_0_xyyyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyyyy_1[i] * wa_x[i];

        g_xzzz_0_xxyyyyz_0[i] = 2.0 * g_zzz_0_xyyyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyz_1[i] * wa_x[i];

        g_xzzz_0_xxyyyzz_0[i] = 2.0 * g_zzz_0_xyyyzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyzz_1[i] * wa_x[i];

        g_xzzz_0_xxyyzzz_0[i] = 2.0 * g_zzz_0_xyyzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyzzz_1[i] * wa_x[i];

        g_xzzz_0_xxyzzzz_0[i] = 2.0 * g_zzz_0_xyzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxzzzzz_0[i] = 2.0 * g_zzz_0_xzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xyyyyyy_0[i] = g_zzz_0_yyyyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyyyy_1[i] * wa_x[i];

        g_xzzz_0_xyyyyyz_0[i] = g_zzz_0_yyyyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyz_1[i] * wa_x[i];

        g_xzzz_0_xyyyyzz_0[i] = g_zzz_0_yyyyzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyzz_1[i] * wa_x[i];

        g_xzzz_0_xyyyzzz_0[i] = g_zzz_0_yyyzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyzzz_1[i] * wa_x[i];

        g_xzzz_0_xyyzzzz_0[i] = g_zzz_0_yyzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyzzzz_1[i] * wa_x[i];

        g_xzzz_0_xyzzzzz_0[i] = g_zzz_0_yzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xzzzzzz_0[i] = g_zzz_0_zzzzzz_1[i] * fi_acd_0 + g_zzz_0_xzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_yyyyyyy_0[i] = g_zzz_0_yyyyyyy_1[i] * wa_x[i];

        g_xzzz_0_yyyyyyz_0[i] = g_zzz_0_yyyyyyz_1[i] * wa_x[i];

        g_xzzz_0_yyyyyzz_0[i] = g_zzz_0_yyyyyzz_1[i] * wa_x[i];

        g_xzzz_0_yyyyzzz_0[i] = g_zzz_0_yyyyzzz_1[i] * wa_x[i];

        g_xzzz_0_yyyzzzz_0[i] = g_zzz_0_yyyzzzz_1[i] * wa_x[i];

        g_xzzz_0_yyzzzzz_0[i] = g_zzz_0_yyzzzzz_1[i] * wa_x[i];

        g_xzzz_0_yzzzzzz_0[i] = g_zzz_0_yzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_zzzzzzz_0[i] = g_zzz_0_zzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 360-396 components of targeted buffer : GSK

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

    #pragma omp simd aligned(g_yy_0_xxxxxxx_0, g_yy_0_xxxxxxx_1, g_yy_0_xxxxxxy_0, g_yy_0_xxxxxxy_1, g_yy_0_xxxxxxz_0, g_yy_0_xxxxxxz_1, g_yy_0_xxxxxyy_0, g_yy_0_xxxxxyy_1, g_yy_0_xxxxxyz_0, g_yy_0_xxxxxyz_1, g_yy_0_xxxxxzz_0, g_yy_0_xxxxxzz_1, g_yy_0_xxxxyyy_0, g_yy_0_xxxxyyy_1, g_yy_0_xxxxyyz_0, g_yy_0_xxxxyyz_1, g_yy_0_xxxxyzz_0, g_yy_0_xxxxyzz_1, g_yy_0_xxxxzzz_0, g_yy_0_xxxxzzz_1, g_yy_0_xxxyyyy_0, g_yy_0_xxxyyyy_1, g_yy_0_xxxyyyz_0, g_yy_0_xxxyyyz_1, g_yy_0_xxxyyzz_0, g_yy_0_xxxyyzz_1, g_yy_0_xxxyzzz_0, g_yy_0_xxxyzzz_1, g_yy_0_xxxzzzz_0, g_yy_0_xxxzzzz_1, g_yy_0_xxyyyyy_0, g_yy_0_xxyyyyy_1, g_yy_0_xxyyyyz_0, g_yy_0_xxyyyyz_1, g_yy_0_xxyyyzz_0, g_yy_0_xxyyyzz_1, g_yy_0_xxyyzzz_0, g_yy_0_xxyyzzz_1, g_yy_0_xxyzzzz_0, g_yy_0_xxyzzzz_1, g_yy_0_xxzzzzz_0, g_yy_0_xxzzzzz_1, g_yy_0_xyyyyyy_0, g_yy_0_xyyyyyy_1, g_yy_0_xyyyyyz_0, g_yy_0_xyyyyyz_1, g_yy_0_xyyyyzz_0, g_yy_0_xyyyyzz_1, g_yy_0_xyyyzzz_0, g_yy_0_xyyyzzz_1, g_yy_0_xyyzzzz_0, g_yy_0_xyyzzzz_1, g_yy_0_xyzzzzz_0, g_yy_0_xyzzzzz_1, g_yy_0_xzzzzzz_0, g_yy_0_xzzzzzz_1, g_yy_0_yyyyyyy_0, g_yy_0_yyyyyyy_1, g_yy_0_yyyyyyz_0, g_yy_0_yyyyyyz_1, g_yy_0_yyyyyzz_0, g_yy_0_yyyyyzz_1, g_yy_0_yyyyzzz_0, g_yy_0_yyyyzzz_1, g_yy_0_yyyzzzz_0, g_yy_0_yyyzzzz_1, g_yy_0_yyzzzzz_0, g_yy_0_yyzzzzz_1, g_yy_0_yzzzzzz_0, g_yy_0_yzzzzzz_1, g_yy_0_zzzzzzz_0, g_yy_0_zzzzzzz_1, g_yyy_0_xxxxxx_1, g_yyy_0_xxxxxxx_1, g_yyy_0_xxxxxxy_1, g_yyy_0_xxxxxxz_1, g_yyy_0_xxxxxy_1, g_yyy_0_xxxxxyy_1, g_yyy_0_xxxxxyz_1, g_yyy_0_xxxxxz_1, g_yyy_0_xxxxxzz_1, g_yyy_0_xxxxyy_1, g_yyy_0_xxxxyyy_1, g_yyy_0_xxxxyyz_1, g_yyy_0_xxxxyz_1, g_yyy_0_xxxxyzz_1, g_yyy_0_xxxxzz_1, g_yyy_0_xxxxzzz_1, g_yyy_0_xxxyyy_1, g_yyy_0_xxxyyyy_1, g_yyy_0_xxxyyyz_1, g_yyy_0_xxxyyz_1, g_yyy_0_xxxyyzz_1, g_yyy_0_xxxyzz_1, g_yyy_0_xxxyzzz_1, g_yyy_0_xxxzzz_1, g_yyy_0_xxxzzzz_1, g_yyy_0_xxyyyy_1, g_yyy_0_xxyyyyy_1, g_yyy_0_xxyyyyz_1, g_yyy_0_xxyyyz_1, g_yyy_0_xxyyyzz_1, g_yyy_0_xxyyzz_1, g_yyy_0_xxyyzzz_1, g_yyy_0_xxyzzz_1, g_yyy_0_xxyzzzz_1, g_yyy_0_xxzzzz_1, g_yyy_0_xxzzzzz_1, g_yyy_0_xyyyyy_1, g_yyy_0_xyyyyyy_1, g_yyy_0_xyyyyyz_1, g_yyy_0_xyyyyz_1, g_yyy_0_xyyyyzz_1, g_yyy_0_xyyyzz_1, g_yyy_0_xyyyzzz_1, g_yyy_0_xyyzzz_1, g_yyy_0_xyyzzzz_1, g_yyy_0_xyzzzz_1, g_yyy_0_xyzzzzz_1, g_yyy_0_xzzzzz_1, g_yyy_0_xzzzzzz_1, g_yyy_0_yyyyyy_1, g_yyy_0_yyyyyyy_1, g_yyy_0_yyyyyyz_1, g_yyy_0_yyyyyz_1, g_yyy_0_yyyyyzz_1, g_yyy_0_yyyyzz_1, g_yyy_0_yyyyzzz_1, g_yyy_0_yyyzzz_1, g_yyy_0_yyyzzzz_1, g_yyy_0_yyzzzz_1, g_yyy_0_yyzzzzz_1, g_yyy_0_yzzzzz_1, g_yyy_0_yzzzzzz_1, g_yyy_0_zzzzzz_1, g_yyy_0_zzzzzzz_1, g_yyyy_0_xxxxxxx_0, g_yyyy_0_xxxxxxy_0, g_yyyy_0_xxxxxxz_0, g_yyyy_0_xxxxxyy_0, g_yyyy_0_xxxxxyz_0, g_yyyy_0_xxxxxzz_0, g_yyyy_0_xxxxyyy_0, g_yyyy_0_xxxxyyz_0, g_yyyy_0_xxxxyzz_0, g_yyyy_0_xxxxzzz_0, g_yyyy_0_xxxyyyy_0, g_yyyy_0_xxxyyyz_0, g_yyyy_0_xxxyyzz_0, g_yyyy_0_xxxyzzz_0, g_yyyy_0_xxxzzzz_0, g_yyyy_0_xxyyyyy_0, g_yyyy_0_xxyyyyz_0, g_yyyy_0_xxyyyzz_0, g_yyyy_0_xxyyzzz_0, g_yyyy_0_xxyzzzz_0, g_yyyy_0_xxzzzzz_0, g_yyyy_0_xyyyyyy_0, g_yyyy_0_xyyyyyz_0, g_yyyy_0_xyyyyzz_0, g_yyyy_0_xyyyzzz_0, g_yyyy_0_xyyzzzz_0, g_yyyy_0_xyzzzzz_0, g_yyyy_0_xzzzzzz_0, g_yyyy_0_yyyyyyy_0, g_yyyy_0_yyyyyyz_0, g_yyyy_0_yyyyyzz_0, g_yyyy_0_yyyyzzz_0, g_yyyy_0_yyyzzzz_0, g_yyyy_0_yyzzzzz_0, g_yyyy_0_yzzzzzz_0, g_yyyy_0_zzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyy_0_xxxxxxx_0[i] = 3.0 * g_yy_0_xxxxxxx_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxx_1[i] * fz_be_0 + g_yyy_0_xxxxxxx_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxy_0[i] = 3.0 * g_yy_0_xxxxxxy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxy_1[i] * fz_be_0 + g_yyy_0_xxxxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxxxy_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxz_0[i] = 3.0 * g_yy_0_xxxxxxz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxz_1[i] * fz_be_0 + g_yyy_0_xxxxxxz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxyy_0[i] = 3.0 * g_yy_0_xxxxxyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxyy_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxxyy_1[i] * wa_y[i];

        g_yyyy_0_xxxxxyz_0[i] = 3.0 * g_yy_0_xxxxxyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxyz_1[i] * fz_be_0 + g_yyy_0_xxxxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxzz_0[i] = 3.0 * g_yy_0_xxxxxzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxzz_1[i] * fz_be_0 + g_yyy_0_xxxxxzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxyyy_0[i] = 3.0 * g_yy_0_xxxxyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyyy_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxxxyy_1[i] * fi_acd_0 + g_yyy_0_xxxxyyy_1[i] * wa_y[i];

        g_yyyy_0_xxxxyyz_0[i] = 3.0 * g_yy_0_xxxxyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxxyz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyz_1[i] * wa_y[i];

        g_yyyy_0_xxxxyzz_0[i] = 3.0 * g_yy_0_xxxxyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyzz_1[i] * fz_be_0 + g_yyy_0_xxxxzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxzzz_0[i] = 3.0 * g_yy_0_xxxxzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxzzz_1[i] * fz_be_0 + g_yyy_0_xxxxzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxyyyy_0[i] = 3.0 * g_yy_0_xxxyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyyy_1[i] * fz_be_0 + 4.0 * g_yyy_0_xxxyyy_1[i] * fi_acd_0 + g_yyy_0_xxxyyyy_1[i] * wa_y[i];

        g_yyyy_0_xxxyyyz_0[i] = 3.0 * g_yy_0_xxxyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxxyyz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyz_1[i] * wa_y[i];

        g_yyyy_0_xxxyyzz_0[i] = 3.0 * g_yy_0_xxxyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxyzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyzz_1[i] * wa_y[i];

        g_yyyy_0_xxxyzzz_0[i] = 3.0 * g_yy_0_xxxyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyzzz_1[i] * fz_be_0 + g_yyy_0_xxxzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxzzzz_0[i] = 3.0 * g_yy_0_xxxzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxzzzz_1[i] * fz_be_0 + g_yyy_0_xxxzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxyyyyy_0[i] = 3.0 * g_yy_0_xxyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyyy_1[i] * fz_be_0 + 5.0 * g_yyy_0_xxyyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyyyy_1[i] * wa_y[i];

        g_yyyy_0_xxyyyyz_0[i] = 3.0 * g_yy_0_xxyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyy_0_xxyyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyz_1[i] * wa_y[i];

        g_yyyy_0_xxyyyzz_0[i] = 3.0 * g_yy_0_xxyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxyyzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyzz_1[i] * wa_y[i];

        g_yyyy_0_xxyyzzz_0[i] = 3.0 * g_yy_0_xxyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxyzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyzzz_1[i] * wa_y[i];

        g_yyyy_0_xxyzzzz_0[i] = 3.0 * g_yy_0_xxyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyzzzz_1[i] * fz_be_0 + g_yyy_0_xxzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxzzzzz_0[i] = 3.0 * g_yy_0_xxzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxzzzzz_1[i] * fz_be_0 + g_yyy_0_xxzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xyyyyyy_0[i] = 3.0 * g_yy_0_xyyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyyy_1[i] * fz_be_0 + 6.0 * g_yyy_0_xyyyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyyyy_1[i] * wa_y[i];

        g_yyyy_0_xyyyyyz_0[i] = 3.0 * g_yy_0_xyyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyy_0_xyyyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyz_1[i] * wa_y[i];

        g_yyyy_0_xyyyyzz_0[i] = 3.0 * g_yy_0_xyyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyy_0_xyyyzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyzz_1[i] * wa_y[i];

        g_yyyy_0_xyyyzzz_0[i] = 3.0 * g_yy_0_xyyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xyyzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyzzz_1[i] * wa_y[i];

        g_yyyy_0_xyyzzzz_0[i] = 3.0 * g_yy_0_xyyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xyzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyzzzz_1[i] * wa_y[i];

        g_yyyy_0_xyzzzzz_0[i] = 3.0 * g_yy_0_xyzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyzzzzz_1[i] * fz_be_0 + g_yyy_0_xzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xzzzzzz_0[i] = 3.0 * g_yy_0_xzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xzzzzzz_1[i] * fz_be_0 + g_yyy_0_xzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_yyyyyyy_0[i] = 3.0 * g_yy_0_yyyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyyy_1[i] * fz_be_0 + 7.0 * g_yyy_0_yyyyyy_1[i] * fi_acd_0 + g_yyy_0_yyyyyyy_1[i] * wa_y[i];

        g_yyyy_0_yyyyyyz_0[i] = 3.0 * g_yy_0_yyyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyy_0_yyyyyz_1[i] * fi_acd_0 + g_yyy_0_yyyyyyz_1[i] * wa_y[i];

        g_yyyy_0_yyyyyzz_0[i] = 3.0 * g_yy_0_yyyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyy_0_yyyyzz_1[i] * fi_acd_0 + g_yyy_0_yyyyyzz_1[i] * wa_y[i];

        g_yyyy_0_yyyyzzz_0[i] = 3.0 * g_yy_0_yyyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyy_0_yyyzzz_1[i] * fi_acd_0 + g_yyy_0_yyyyzzz_1[i] * wa_y[i];

        g_yyyy_0_yyyzzzz_0[i] = 3.0 * g_yy_0_yyyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_yyzzzz_1[i] * fi_acd_0 + g_yyy_0_yyyzzzz_1[i] * wa_y[i];

        g_yyyy_0_yyzzzzz_0[i] = 3.0 * g_yy_0_yyzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_yzzzzz_1[i] * fi_acd_0 + g_yyy_0_yyzzzzz_1[i] * wa_y[i];

        g_yyyy_0_yzzzzzz_0[i] = 3.0 * g_yy_0_yzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yzzzzzz_1[i] * fz_be_0 + g_yyy_0_zzzzzz_1[i] * fi_acd_0 + g_yyy_0_yzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_zzzzzzz_0[i] = 3.0 * g_yy_0_zzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_zzzzzzz_1[i] * fz_be_0 + g_yyy_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 396-432 components of targeted buffer : GSK

    auto g_yyyz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 396);

    auto g_yyyz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 397);

    auto g_yyyz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 398);

    auto g_yyyz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 399);

    auto g_yyyz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 400);

    auto g_yyyz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 401);

    auto g_yyyz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 402);

    auto g_yyyz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 403);

    auto g_yyyz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 404);

    auto g_yyyz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 405);

    auto g_yyyz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 406);

    auto g_yyyz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 407);

    auto g_yyyz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 408);

    auto g_yyyz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 409);

    auto g_yyyz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 410);

    auto g_yyyz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 411);

    auto g_yyyz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 412);

    auto g_yyyz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 413);

    auto g_yyyz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 414);

    auto g_yyyz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 415);

    auto g_yyyz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 416);

    auto g_yyyz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 417);

    auto g_yyyz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 418);

    auto g_yyyz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 419);

    auto g_yyyz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 420);

    auto g_yyyz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 421);

    auto g_yyyz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 422);

    auto g_yyyz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 423);

    auto g_yyyz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 424);

    auto g_yyyz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 425);

    auto g_yyyz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 426);

    auto g_yyyz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 427);

    auto g_yyyz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 428);

    auto g_yyyz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 429);

    auto g_yyyz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 430);

    auto g_yyyz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 431);

    #pragma omp simd aligned(g_yyy_0_xxxxxx_1, g_yyy_0_xxxxxxx_1, g_yyy_0_xxxxxxy_1, g_yyy_0_xxxxxxz_1, g_yyy_0_xxxxxy_1, g_yyy_0_xxxxxyy_1, g_yyy_0_xxxxxyz_1, g_yyy_0_xxxxxz_1, g_yyy_0_xxxxxzz_1, g_yyy_0_xxxxyy_1, g_yyy_0_xxxxyyy_1, g_yyy_0_xxxxyyz_1, g_yyy_0_xxxxyz_1, g_yyy_0_xxxxyzz_1, g_yyy_0_xxxxzz_1, g_yyy_0_xxxxzzz_1, g_yyy_0_xxxyyy_1, g_yyy_0_xxxyyyy_1, g_yyy_0_xxxyyyz_1, g_yyy_0_xxxyyz_1, g_yyy_0_xxxyyzz_1, g_yyy_0_xxxyzz_1, g_yyy_0_xxxyzzz_1, g_yyy_0_xxxzzz_1, g_yyy_0_xxxzzzz_1, g_yyy_0_xxyyyy_1, g_yyy_0_xxyyyyy_1, g_yyy_0_xxyyyyz_1, g_yyy_0_xxyyyz_1, g_yyy_0_xxyyyzz_1, g_yyy_0_xxyyzz_1, g_yyy_0_xxyyzzz_1, g_yyy_0_xxyzzz_1, g_yyy_0_xxyzzzz_1, g_yyy_0_xxzzzz_1, g_yyy_0_xxzzzzz_1, g_yyy_0_xyyyyy_1, g_yyy_0_xyyyyyy_1, g_yyy_0_xyyyyyz_1, g_yyy_0_xyyyyz_1, g_yyy_0_xyyyyzz_1, g_yyy_0_xyyyzz_1, g_yyy_0_xyyyzzz_1, g_yyy_0_xyyzzz_1, g_yyy_0_xyyzzzz_1, g_yyy_0_xyzzzz_1, g_yyy_0_xyzzzzz_1, g_yyy_0_xzzzzz_1, g_yyy_0_xzzzzzz_1, g_yyy_0_yyyyyy_1, g_yyy_0_yyyyyyy_1, g_yyy_0_yyyyyyz_1, g_yyy_0_yyyyyz_1, g_yyy_0_yyyyyzz_1, g_yyy_0_yyyyzz_1, g_yyy_0_yyyyzzz_1, g_yyy_0_yyyzzz_1, g_yyy_0_yyyzzzz_1, g_yyy_0_yyzzzz_1, g_yyy_0_yyzzzzz_1, g_yyy_0_yzzzzz_1, g_yyy_0_yzzzzzz_1, g_yyy_0_zzzzzz_1, g_yyy_0_zzzzzzz_1, g_yyyz_0_xxxxxxx_0, g_yyyz_0_xxxxxxy_0, g_yyyz_0_xxxxxxz_0, g_yyyz_0_xxxxxyy_0, g_yyyz_0_xxxxxyz_0, g_yyyz_0_xxxxxzz_0, g_yyyz_0_xxxxyyy_0, g_yyyz_0_xxxxyyz_0, g_yyyz_0_xxxxyzz_0, g_yyyz_0_xxxxzzz_0, g_yyyz_0_xxxyyyy_0, g_yyyz_0_xxxyyyz_0, g_yyyz_0_xxxyyzz_0, g_yyyz_0_xxxyzzz_0, g_yyyz_0_xxxzzzz_0, g_yyyz_0_xxyyyyy_0, g_yyyz_0_xxyyyyz_0, g_yyyz_0_xxyyyzz_0, g_yyyz_0_xxyyzzz_0, g_yyyz_0_xxyzzzz_0, g_yyyz_0_xxzzzzz_0, g_yyyz_0_xyyyyyy_0, g_yyyz_0_xyyyyyz_0, g_yyyz_0_xyyyyzz_0, g_yyyz_0_xyyyzzz_0, g_yyyz_0_xyyzzzz_0, g_yyyz_0_xyzzzzz_0, g_yyyz_0_xzzzzzz_0, g_yyyz_0_yyyyyyy_0, g_yyyz_0_yyyyyyz_0, g_yyyz_0_yyyyyzz_0, g_yyyz_0_yyyyzzz_0, g_yyyz_0_yyyzzzz_0, g_yyyz_0_yyzzzzz_0, g_yyyz_0_yzzzzzz_0, g_yyyz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyz_0_xxxxxxx_0[i] = g_yyy_0_xxxxxxx_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxy_0[i] = g_yyy_0_xxxxxxy_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxz_0[i] = g_yyy_0_xxxxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxxxz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxyy_0[i] = g_yyy_0_xxxxxyy_1[i] * wa_z[i];

        g_yyyz_0_xxxxxyz_0[i] = g_yyy_0_xxxxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxxyz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxzz_0[i] = 2.0 * g_yyy_0_xxxxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxxzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxyyy_0[i] = g_yyy_0_xxxxyyy_1[i] * wa_z[i];

        g_yyyz_0_xxxxyyz_0[i] = g_yyy_0_xxxxyy_1[i] * fi_acd_0 + g_yyy_0_xxxxyyz_1[i] * wa_z[i];

        g_yyyz_0_xxxxyzz_0[i] = 2.0 * g_yyy_0_xxxxyz_1[i] * fi_acd_0 + g_yyy_0_xxxxyzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxzzz_0[i] = 3.0 * g_yyy_0_xxxxzz_1[i] * fi_acd_0 + g_yyy_0_xxxxzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxyyyy_0[i] = g_yyy_0_xxxyyyy_1[i] * wa_z[i];

        g_yyyz_0_xxxyyyz_0[i] = g_yyy_0_xxxyyy_1[i] * fi_acd_0 + g_yyy_0_xxxyyyz_1[i] * wa_z[i];

        g_yyyz_0_xxxyyzz_0[i] = 2.0 * g_yyy_0_xxxyyz_1[i] * fi_acd_0 + g_yyy_0_xxxyyzz_1[i] * wa_z[i];

        g_yyyz_0_xxxyzzz_0[i] = 3.0 * g_yyy_0_xxxyzz_1[i] * fi_acd_0 + g_yyy_0_xxxyzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxzzzz_0[i] = 4.0 * g_yyy_0_xxxzzz_1[i] * fi_acd_0 + g_yyy_0_xxxzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxyyyyy_0[i] = g_yyy_0_xxyyyyy_1[i] * wa_z[i];

        g_yyyz_0_xxyyyyz_0[i] = g_yyy_0_xxyyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyyyz_1[i] * wa_z[i];

        g_yyyz_0_xxyyyzz_0[i] = 2.0 * g_yyy_0_xxyyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyyzz_1[i] * wa_z[i];

        g_yyyz_0_xxyyzzz_0[i] = 3.0 * g_yyy_0_xxyyzz_1[i] * fi_acd_0 + g_yyy_0_xxyyzzz_1[i] * wa_z[i];

        g_yyyz_0_xxyzzzz_0[i] = 4.0 * g_yyy_0_xxyzzz_1[i] * fi_acd_0 + g_yyy_0_xxyzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxzzzzz_0[i] = 5.0 * g_yyy_0_xxzzzz_1[i] * fi_acd_0 + g_yyy_0_xxzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xyyyyyy_0[i] = g_yyy_0_xyyyyyy_1[i] * wa_z[i];

        g_yyyz_0_xyyyyyz_0[i] = g_yyy_0_xyyyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyyyz_1[i] * wa_z[i];

        g_yyyz_0_xyyyyzz_0[i] = 2.0 * g_yyy_0_xyyyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyyzz_1[i] * wa_z[i];

        g_yyyz_0_xyyyzzz_0[i] = 3.0 * g_yyy_0_xyyyzz_1[i] * fi_acd_0 + g_yyy_0_xyyyzzz_1[i] * wa_z[i];

        g_yyyz_0_xyyzzzz_0[i] = 4.0 * g_yyy_0_xyyzzz_1[i] * fi_acd_0 + g_yyy_0_xyyzzzz_1[i] * wa_z[i];

        g_yyyz_0_xyzzzzz_0[i] = 5.0 * g_yyy_0_xyzzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xzzzzzz_0[i] = 6.0 * g_yyy_0_xzzzzz_1[i] * fi_acd_0 + g_yyy_0_xzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_yyyyyyy_0[i] = g_yyy_0_yyyyyyy_1[i] * wa_z[i];

        g_yyyz_0_yyyyyyz_0[i] = g_yyy_0_yyyyyy_1[i] * fi_acd_0 + g_yyy_0_yyyyyyz_1[i] * wa_z[i];

        g_yyyz_0_yyyyyzz_0[i] = 2.0 * g_yyy_0_yyyyyz_1[i] * fi_acd_0 + g_yyy_0_yyyyyzz_1[i] * wa_z[i];

        g_yyyz_0_yyyyzzz_0[i] = 3.0 * g_yyy_0_yyyyzz_1[i] * fi_acd_0 + g_yyy_0_yyyyzzz_1[i] * wa_z[i];

        g_yyyz_0_yyyzzzz_0[i] = 4.0 * g_yyy_0_yyyzzz_1[i] * fi_acd_0 + g_yyy_0_yyyzzzz_1[i] * wa_z[i];

        g_yyyz_0_yyzzzzz_0[i] = 5.0 * g_yyy_0_yyzzzz_1[i] * fi_acd_0 + g_yyy_0_yyzzzzz_1[i] * wa_z[i];

        g_yyyz_0_yzzzzzz_0[i] = 6.0 * g_yyy_0_yzzzzz_1[i] * fi_acd_0 + g_yyy_0_yzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_zzzzzzz_0[i] = 7.0 * g_yyy_0_zzzzzz_1[i] * fi_acd_0 + g_yyy_0_zzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 432-468 components of targeted buffer : GSK

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

    #pragma omp simd aligned(g_yy_0_xxxxxxy_0, g_yy_0_xxxxxxy_1, g_yy_0_xxxxxyy_0, g_yy_0_xxxxxyy_1, g_yy_0_xxxxyyy_0, g_yy_0_xxxxyyy_1, g_yy_0_xxxyyyy_0, g_yy_0_xxxyyyy_1, g_yy_0_xxyyyyy_0, g_yy_0_xxyyyyy_1, g_yy_0_xyyyyyy_0, g_yy_0_xyyyyyy_1, g_yy_0_yyyyyyy_0, g_yy_0_yyyyyyy_1, g_yyz_0_xxxxxxy_1, g_yyz_0_xxxxxyy_1, g_yyz_0_xxxxyyy_1, g_yyz_0_xxxyyyy_1, g_yyz_0_xxyyyyy_1, g_yyz_0_xyyyyyy_1, g_yyz_0_yyyyyyy_1, g_yyzz_0_xxxxxxx_0, g_yyzz_0_xxxxxxy_0, g_yyzz_0_xxxxxxz_0, g_yyzz_0_xxxxxyy_0, g_yyzz_0_xxxxxyz_0, g_yyzz_0_xxxxxzz_0, g_yyzz_0_xxxxyyy_0, g_yyzz_0_xxxxyyz_0, g_yyzz_0_xxxxyzz_0, g_yyzz_0_xxxxzzz_0, g_yyzz_0_xxxyyyy_0, g_yyzz_0_xxxyyyz_0, g_yyzz_0_xxxyyzz_0, g_yyzz_0_xxxyzzz_0, g_yyzz_0_xxxzzzz_0, g_yyzz_0_xxyyyyy_0, g_yyzz_0_xxyyyyz_0, g_yyzz_0_xxyyyzz_0, g_yyzz_0_xxyyzzz_0, g_yyzz_0_xxyzzzz_0, g_yyzz_0_xxzzzzz_0, g_yyzz_0_xyyyyyy_0, g_yyzz_0_xyyyyyz_0, g_yyzz_0_xyyyyzz_0, g_yyzz_0_xyyyzzz_0, g_yyzz_0_xyyzzzz_0, g_yyzz_0_xyzzzzz_0, g_yyzz_0_xzzzzzz_0, g_yyzz_0_yyyyyyy_0, g_yyzz_0_yyyyyyz_0, g_yyzz_0_yyyyyzz_0, g_yyzz_0_yyyyzzz_0, g_yyzz_0_yyyzzzz_0, g_yyzz_0_yyzzzzz_0, g_yyzz_0_yzzzzzz_0, g_yyzz_0_zzzzzzz_0, g_yzz_0_xxxxxxx_1, g_yzz_0_xxxxxxz_1, g_yzz_0_xxxxxyz_1, g_yzz_0_xxxxxz_1, g_yzz_0_xxxxxzz_1, g_yzz_0_xxxxyyz_1, g_yzz_0_xxxxyz_1, g_yzz_0_xxxxyzz_1, g_yzz_0_xxxxzz_1, g_yzz_0_xxxxzzz_1, g_yzz_0_xxxyyyz_1, g_yzz_0_xxxyyz_1, g_yzz_0_xxxyyzz_1, g_yzz_0_xxxyzz_1, g_yzz_0_xxxyzzz_1, g_yzz_0_xxxzzz_1, g_yzz_0_xxxzzzz_1, g_yzz_0_xxyyyyz_1, g_yzz_0_xxyyyz_1, g_yzz_0_xxyyyzz_1, g_yzz_0_xxyyzz_1, g_yzz_0_xxyyzzz_1, g_yzz_0_xxyzzz_1, g_yzz_0_xxyzzzz_1, g_yzz_0_xxzzzz_1, g_yzz_0_xxzzzzz_1, g_yzz_0_xyyyyyz_1, g_yzz_0_xyyyyz_1, g_yzz_0_xyyyyzz_1, g_yzz_0_xyyyzz_1, g_yzz_0_xyyyzzz_1, g_yzz_0_xyyzzz_1, g_yzz_0_xyyzzzz_1, g_yzz_0_xyzzzz_1, g_yzz_0_xyzzzzz_1, g_yzz_0_xzzzzz_1, g_yzz_0_xzzzzzz_1, g_yzz_0_yyyyyyz_1, g_yzz_0_yyyyyz_1, g_yzz_0_yyyyyzz_1, g_yzz_0_yyyyzz_1, g_yzz_0_yyyyzzz_1, g_yzz_0_yyyzzz_1, g_yzz_0_yyyzzzz_1, g_yzz_0_yyzzzz_1, g_yzz_0_yyzzzzz_1, g_yzz_0_yzzzzz_1, g_yzz_0_yzzzzzz_1, g_yzz_0_zzzzzz_1, g_yzz_0_zzzzzzz_1, g_zz_0_xxxxxxx_0, g_zz_0_xxxxxxx_1, g_zz_0_xxxxxxz_0, g_zz_0_xxxxxxz_1, g_zz_0_xxxxxyz_0, g_zz_0_xxxxxyz_1, g_zz_0_xxxxxzz_0, g_zz_0_xxxxxzz_1, g_zz_0_xxxxyyz_0, g_zz_0_xxxxyyz_1, g_zz_0_xxxxyzz_0, g_zz_0_xxxxyzz_1, g_zz_0_xxxxzzz_0, g_zz_0_xxxxzzz_1, g_zz_0_xxxyyyz_0, g_zz_0_xxxyyyz_1, g_zz_0_xxxyyzz_0, g_zz_0_xxxyyzz_1, g_zz_0_xxxyzzz_0, g_zz_0_xxxyzzz_1, g_zz_0_xxxzzzz_0, g_zz_0_xxxzzzz_1, g_zz_0_xxyyyyz_0, g_zz_0_xxyyyyz_1, g_zz_0_xxyyyzz_0, g_zz_0_xxyyyzz_1, g_zz_0_xxyyzzz_0, g_zz_0_xxyyzzz_1, g_zz_0_xxyzzzz_0, g_zz_0_xxyzzzz_1, g_zz_0_xxzzzzz_0, g_zz_0_xxzzzzz_1, g_zz_0_xyyyyyz_0, g_zz_0_xyyyyyz_1, g_zz_0_xyyyyzz_0, g_zz_0_xyyyyzz_1, g_zz_0_xyyyzzz_0, g_zz_0_xyyyzzz_1, g_zz_0_xyyzzzz_0, g_zz_0_xyyzzzz_1, g_zz_0_xyzzzzz_0, g_zz_0_xyzzzzz_1, g_zz_0_xzzzzzz_0, g_zz_0_xzzzzzz_1, g_zz_0_yyyyyyz_0, g_zz_0_yyyyyyz_1, g_zz_0_yyyyyzz_0, g_zz_0_yyyyyzz_1, g_zz_0_yyyyzzz_0, g_zz_0_yyyyzzz_1, g_zz_0_yyyzzzz_0, g_zz_0_yyyzzzz_1, g_zz_0_yyzzzzz_0, g_zz_0_yyzzzzz_1, g_zz_0_yzzzzzz_0, g_zz_0_yzzzzzz_1, g_zz_0_zzzzzzz_0, g_zz_0_zzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzz_0_xxxxxxx_0[i] = g_zz_0_xxxxxxx_0[i] * fbe_0 - g_zz_0_xxxxxxx_1[i] * fz_be_0 + g_yzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yyzz_0_xxxxxxy_0[i] = g_yy_0_xxxxxxy_0[i] * fbe_0 - g_yy_0_xxxxxxy_1[i] * fz_be_0 + g_yyz_0_xxxxxxy_1[i] * wa_z[i];

        g_yyzz_0_xxxxxxz_0[i] = g_zz_0_xxxxxxz_0[i] * fbe_0 - g_zz_0_xxxxxxz_1[i] * fz_be_0 + g_yzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxyy_0[i] = g_yy_0_xxxxxyy_0[i] * fbe_0 - g_yy_0_xxxxxyy_1[i] * fz_be_0 + g_yyz_0_xxxxxyy_1[i] * wa_z[i];

        g_yyzz_0_xxxxxyz_0[i] = g_zz_0_xxxxxyz_0[i] * fbe_0 - g_zz_0_xxxxxyz_1[i] * fz_be_0 + g_yzz_0_xxxxxz_1[i] * fi_acd_0 + g_yzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxzz_0[i] = g_zz_0_xxxxxzz_0[i] * fbe_0 - g_zz_0_xxxxxzz_1[i] * fz_be_0 + g_yzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxyyy_0[i] = g_yy_0_xxxxyyy_0[i] * fbe_0 - g_yy_0_xxxxyyy_1[i] * fz_be_0 + g_yyz_0_xxxxyyy_1[i] * wa_z[i];

        g_yyzz_0_xxxxyyz_0[i] = g_zz_0_xxxxyyz_0[i] * fbe_0 - g_zz_0_xxxxyyz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxxxyz_1[i] * fi_acd_0 + g_yzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yyzz_0_xxxxyzz_0[i] = g_zz_0_xxxxyzz_0[i] * fbe_0 - g_zz_0_xxxxyzz_1[i] * fz_be_0 + g_yzz_0_xxxxzz_1[i] * fi_acd_0 + g_yzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxzzz_0[i] = g_zz_0_xxxxzzz_0[i] * fbe_0 - g_zz_0_xxxxzzz_1[i] * fz_be_0 + g_yzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxyyyy_0[i] = g_yy_0_xxxyyyy_0[i] * fbe_0 - g_yy_0_xxxyyyy_1[i] * fz_be_0 + g_yyz_0_xxxyyyy_1[i] * wa_z[i];

        g_yyzz_0_xxxyyyz_0[i] = g_zz_0_xxxyyyz_0[i] * fbe_0 - g_zz_0_xxxyyyz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xxxyyz_1[i] * fi_acd_0 + g_yzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yyzz_0_xxxyyzz_0[i] = g_zz_0_xxxyyzz_0[i] * fbe_0 - g_zz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxxyzz_1[i] * fi_acd_0 + g_yzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yyzz_0_xxxyzzz_0[i] = g_zz_0_xxxyzzz_0[i] * fbe_0 - g_zz_0_xxxyzzz_1[i] * fz_be_0 + g_yzz_0_xxxzzz_1[i] * fi_acd_0 + g_yzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxzzzz_0[i] = g_zz_0_xxxzzzz_0[i] * fbe_0 - g_zz_0_xxxzzzz_1[i] * fz_be_0 + g_yzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxyyyyy_0[i] = g_yy_0_xxyyyyy_0[i] * fbe_0 - g_yy_0_xxyyyyy_1[i] * fz_be_0 + g_yyz_0_xxyyyyy_1[i] * wa_z[i];

        g_yyzz_0_xxyyyyz_0[i] = g_zz_0_xxyyyyz_0[i] * fbe_0 - g_zz_0_xxyyyyz_1[i] * fz_be_0 + 4.0 * g_yzz_0_xxyyyz_1[i] * fi_acd_0 + g_yzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yyzz_0_xxyyyzz_0[i] = g_zz_0_xxyyyzz_0[i] * fbe_0 - g_zz_0_xxyyyzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xxyyzz_1[i] * fi_acd_0 + g_yzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yyzz_0_xxyyzzz_0[i] = g_zz_0_xxyyzzz_0[i] * fbe_0 - g_zz_0_xxyyzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxyzzz_1[i] * fi_acd_0 + g_yzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yyzz_0_xxyzzzz_0[i] = g_zz_0_xxyzzzz_0[i] * fbe_0 - g_zz_0_xxyzzzz_1[i] * fz_be_0 + g_yzz_0_xxzzzz_1[i] * fi_acd_0 + g_yzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxzzzzz_0[i] = g_zz_0_xxzzzzz_0[i] * fbe_0 - g_zz_0_xxzzzzz_1[i] * fz_be_0 + g_yzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xyyyyyy_0[i] = g_yy_0_xyyyyyy_0[i] * fbe_0 - g_yy_0_xyyyyyy_1[i] * fz_be_0 + g_yyz_0_xyyyyyy_1[i] * wa_z[i];

        g_yyzz_0_xyyyyyz_0[i] = g_zz_0_xyyyyyz_0[i] * fbe_0 - g_zz_0_xyyyyyz_1[i] * fz_be_0 + 5.0 * g_yzz_0_xyyyyz_1[i] * fi_acd_0 + g_yzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yyzz_0_xyyyyzz_0[i] = g_zz_0_xyyyyzz_0[i] * fbe_0 - g_zz_0_xyyyyzz_1[i] * fz_be_0 + 4.0 * g_yzz_0_xyyyzz_1[i] * fi_acd_0 + g_yzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yyzz_0_xyyyzzz_0[i] = g_zz_0_xyyyzzz_0[i] * fbe_0 - g_zz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xyyzzz_1[i] * fi_acd_0 + g_yzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yyzz_0_xyyzzzz_0[i] = g_zz_0_xyyzzzz_0[i] * fbe_0 - g_zz_0_xyyzzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xyzzzz_1[i] * fi_acd_0 + g_yzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yyzz_0_xyzzzzz_0[i] = g_zz_0_xyzzzzz_0[i] * fbe_0 - g_zz_0_xyzzzzz_1[i] * fz_be_0 + g_yzz_0_xzzzzz_1[i] * fi_acd_0 + g_yzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xzzzzzz_0[i] = g_zz_0_xzzzzzz_0[i] * fbe_0 - g_zz_0_xzzzzzz_1[i] * fz_be_0 + g_yzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_yyyyyyy_0[i] = g_yy_0_yyyyyyy_0[i] * fbe_0 - g_yy_0_yyyyyyy_1[i] * fz_be_0 + g_yyz_0_yyyyyyy_1[i] * wa_z[i];

        g_yyzz_0_yyyyyyz_0[i] = g_zz_0_yyyyyyz_0[i] * fbe_0 - g_zz_0_yyyyyyz_1[i] * fz_be_0 + 6.0 * g_yzz_0_yyyyyz_1[i] * fi_acd_0 + g_yzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yyzz_0_yyyyyzz_0[i] = g_zz_0_yyyyyzz_0[i] * fbe_0 - g_zz_0_yyyyyzz_1[i] * fz_be_0 + 5.0 * g_yzz_0_yyyyzz_1[i] * fi_acd_0 + g_yzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yyzz_0_yyyyzzz_0[i] = g_zz_0_yyyyzzz_0[i] * fbe_0 - g_zz_0_yyyyzzz_1[i] * fz_be_0 + 4.0 * g_yzz_0_yyyzzz_1[i] * fi_acd_0 + g_yzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yyzz_0_yyyzzzz_0[i] = g_zz_0_yyyzzzz_0[i] * fbe_0 - g_zz_0_yyyzzzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_yyzzzz_1[i] * fi_acd_0 + g_yzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yyzz_0_yyzzzzz_0[i] = g_zz_0_yyzzzzz_0[i] * fbe_0 - g_zz_0_yyzzzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_yzzzzz_1[i] * fi_acd_0 + g_yzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yyzz_0_yzzzzzz_0[i] = g_zz_0_yzzzzzz_0[i] * fbe_0 - g_zz_0_yzzzzzz_1[i] * fz_be_0 + g_yzz_0_zzzzzz_1[i] * fi_acd_0 + g_yzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_zzzzzzz_0[i] = g_zz_0_zzzzzzz_0[i] * fbe_0 - g_zz_0_zzzzzzz_1[i] * fz_be_0 + g_yzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 468-504 components of targeted buffer : GSK

    auto g_yzzz_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_gsk + 468);

    auto g_yzzz_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_gsk + 469);

    auto g_yzzz_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_gsk + 470);

    auto g_yzzz_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_gsk + 471);

    auto g_yzzz_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_gsk + 472);

    auto g_yzzz_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_gsk + 473);

    auto g_yzzz_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_gsk + 474);

    auto g_yzzz_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_gsk + 475);

    auto g_yzzz_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_gsk + 476);

    auto g_yzzz_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_gsk + 477);

    auto g_yzzz_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_gsk + 478);

    auto g_yzzz_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_gsk + 479);

    auto g_yzzz_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_gsk + 480);

    auto g_yzzz_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_gsk + 481);

    auto g_yzzz_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_gsk + 482);

    auto g_yzzz_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 483);

    auto g_yzzz_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 484);

    auto g_yzzz_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 485);

    auto g_yzzz_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 486);

    auto g_yzzz_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 487);

    auto g_yzzz_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 488);

    auto g_yzzz_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 489);

    auto g_yzzz_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 490);

    auto g_yzzz_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 491);

    auto g_yzzz_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 492);

    auto g_yzzz_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 493);

    auto g_yzzz_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 494);

    auto g_yzzz_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 495);

    auto g_yzzz_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_gsk + 496);

    auto g_yzzz_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_gsk + 497);

    auto g_yzzz_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_gsk + 498);

    auto g_yzzz_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_gsk + 499);

    auto g_yzzz_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_gsk + 500);

    auto g_yzzz_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 501);

    auto g_yzzz_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 502);

    auto g_yzzz_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_gsk + 503);

    #pragma omp simd aligned(g_yzzz_0_xxxxxxx_0, g_yzzz_0_xxxxxxy_0, g_yzzz_0_xxxxxxz_0, g_yzzz_0_xxxxxyy_0, g_yzzz_0_xxxxxyz_0, g_yzzz_0_xxxxxzz_0, g_yzzz_0_xxxxyyy_0, g_yzzz_0_xxxxyyz_0, g_yzzz_0_xxxxyzz_0, g_yzzz_0_xxxxzzz_0, g_yzzz_0_xxxyyyy_0, g_yzzz_0_xxxyyyz_0, g_yzzz_0_xxxyyzz_0, g_yzzz_0_xxxyzzz_0, g_yzzz_0_xxxzzzz_0, g_yzzz_0_xxyyyyy_0, g_yzzz_0_xxyyyyz_0, g_yzzz_0_xxyyyzz_0, g_yzzz_0_xxyyzzz_0, g_yzzz_0_xxyzzzz_0, g_yzzz_0_xxzzzzz_0, g_yzzz_0_xyyyyyy_0, g_yzzz_0_xyyyyyz_0, g_yzzz_0_xyyyyzz_0, g_yzzz_0_xyyyzzz_0, g_yzzz_0_xyyzzzz_0, g_yzzz_0_xyzzzzz_0, g_yzzz_0_xzzzzzz_0, g_yzzz_0_yyyyyyy_0, g_yzzz_0_yyyyyyz_0, g_yzzz_0_yyyyyzz_0, g_yzzz_0_yyyyzzz_0, g_yzzz_0_yyyzzzz_0, g_yzzz_0_yyzzzzz_0, g_yzzz_0_yzzzzzz_0, g_yzzz_0_zzzzzzz_0, g_zzz_0_xxxxxx_1, g_zzz_0_xxxxxxx_1, g_zzz_0_xxxxxxy_1, g_zzz_0_xxxxxxz_1, g_zzz_0_xxxxxy_1, g_zzz_0_xxxxxyy_1, g_zzz_0_xxxxxyz_1, g_zzz_0_xxxxxz_1, g_zzz_0_xxxxxzz_1, g_zzz_0_xxxxyy_1, g_zzz_0_xxxxyyy_1, g_zzz_0_xxxxyyz_1, g_zzz_0_xxxxyz_1, g_zzz_0_xxxxyzz_1, g_zzz_0_xxxxzz_1, g_zzz_0_xxxxzzz_1, g_zzz_0_xxxyyy_1, g_zzz_0_xxxyyyy_1, g_zzz_0_xxxyyyz_1, g_zzz_0_xxxyyz_1, g_zzz_0_xxxyyzz_1, g_zzz_0_xxxyzz_1, g_zzz_0_xxxyzzz_1, g_zzz_0_xxxzzz_1, g_zzz_0_xxxzzzz_1, g_zzz_0_xxyyyy_1, g_zzz_0_xxyyyyy_1, g_zzz_0_xxyyyyz_1, g_zzz_0_xxyyyz_1, g_zzz_0_xxyyyzz_1, g_zzz_0_xxyyzz_1, g_zzz_0_xxyyzzz_1, g_zzz_0_xxyzzz_1, g_zzz_0_xxyzzzz_1, g_zzz_0_xxzzzz_1, g_zzz_0_xxzzzzz_1, g_zzz_0_xyyyyy_1, g_zzz_0_xyyyyyy_1, g_zzz_0_xyyyyyz_1, g_zzz_0_xyyyyz_1, g_zzz_0_xyyyyzz_1, g_zzz_0_xyyyzz_1, g_zzz_0_xyyyzzz_1, g_zzz_0_xyyzzz_1, g_zzz_0_xyyzzzz_1, g_zzz_0_xyzzzz_1, g_zzz_0_xyzzzzz_1, g_zzz_0_xzzzzz_1, g_zzz_0_xzzzzzz_1, g_zzz_0_yyyyyy_1, g_zzz_0_yyyyyyy_1, g_zzz_0_yyyyyyz_1, g_zzz_0_yyyyyz_1, g_zzz_0_yyyyyzz_1, g_zzz_0_yyyyzz_1, g_zzz_0_yyyyzzz_1, g_zzz_0_yyyzzz_1, g_zzz_0_yyyzzzz_1, g_zzz_0_yyzzzz_1, g_zzz_0_yyzzzzz_1, g_zzz_0_yzzzzz_1, g_zzz_0_yzzzzzz_1, g_zzz_0_zzzzzz_1, g_zzz_0_zzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzz_0_xxxxxxx_0[i] = g_zzz_0_xxxxxxx_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxy_0[i] = g_zzz_0_xxxxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxxxy_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxz_0[i] = g_zzz_0_xxxxxxz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxyy_0[i] = 2.0 * g_zzz_0_xxxxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxxyy_1[i] * wa_y[i];

        g_yzzz_0_xxxxxyz_0[i] = g_zzz_0_xxxxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxzz_0[i] = g_zzz_0_xxxxxzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxyyy_0[i] = 3.0 * g_zzz_0_xxxxyy_1[i] * fi_acd_0 + g_zzz_0_xxxxyyy_1[i] * wa_y[i];

        g_yzzz_0_xxxxyyz_0[i] = 2.0 * g_zzz_0_xxxxyz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyz_1[i] * wa_y[i];

        g_yzzz_0_xxxxyzz_0[i] = g_zzz_0_xxxxzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxzzz_0[i] = g_zzz_0_xxxxzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxyyyy_0[i] = 4.0 * g_zzz_0_xxxyyy_1[i] * fi_acd_0 + g_zzz_0_xxxyyyy_1[i] * wa_y[i];

        g_yzzz_0_xxxyyyz_0[i] = 3.0 * g_zzz_0_xxxyyz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyz_1[i] * wa_y[i];

        g_yzzz_0_xxxyyzz_0[i] = 2.0 * g_zzz_0_xxxyzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyzz_1[i] * wa_y[i];

        g_yzzz_0_xxxyzzz_0[i] = g_zzz_0_xxxzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxzzzz_0[i] = g_zzz_0_xxxzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxyyyyy_0[i] = 5.0 * g_zzz_0_xxyyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyyyy_1[i] * wa_y[i];

        g_yzzz_0_xxyyyyz_0[i] = 4.0 * g_zzz_0_xxyyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyz_1[i] * wa_y[i];

        g_yzzz_0_xxyyyzz_0[i] = 3.0 * g_zzz_0_xxyyzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyzz_1[i] * wa_y[i];

        g_yzzz_0_xxyyzzz_0[i] = 2.0 * g_zzz_0_xxyzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyzzz_1[i] * wa_y[i];

        g_yzzz_0_xxyzzzz_0[i] = g_zzz_0_xxzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxzzzzz_0[i] = g_zzz_0_xxzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xyyyyyy_0[i] = 6.0 * g_zzz_0_xyyyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyyyy_1[i] * wa_y[i];

        g_yzzz_0_xyyyyyz_0[i] = 5.0 * g_zzz_0_xyyyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyz_1[i] * wa_y[i];

        g_yzzz_0_xyyyyzz_0[i] = 4.0 * g_zzz_0_xyyyzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyzz_1[i] * wa_y[i];

        g_yzzz_0_xyyyzzz_0[i] = 3.0 * g_zzz_0_xyyzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyzzz_1[i] * wa_y[i];

        g_yzzz_0_xyyzzzz_0[i] = 2.0 * g_zzz_0_xyzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyzzzz_1[i] * wa_y[i];

        g_yzzz_0_xyzzzzz_0[i] = g_zzz_0_xzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xzzzzzz_0[i] = g_zzz_0_xzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_yyyyyyy_0[i] = 7.0 * g_zzz_0_yyyyyy_1[i] * fi_acd_0 + g_zzz_0_yyyyyyy_1[i] * wa_y[i];

        g_yzzz_0_yyyyyyz_0[i] = 6.0 * g_zzz_0_yyyyyz_1[i] * fi_acd_0 + g_zzz_0_yyyyyyz_1[i] * wa_y[i];

        g_yzzz_0_yyyyyzz_0[i] = 5.0 * g_zzz_0_yyyyzz_1[i] * fi_acd_0 + g_zzz_0_yyyyyzz_1[i] * wa_y[i];

        g_yzzz_0_yyyyzzz_0[i] = 4.0 * g_zzz_0_yyyzzz_1[i] * fi_acd_0 + g_zzz_0_yyyyzzz_1[i] * wa_y[i];

        g_yzzz_0_yyyzzzz_0[i] = 3.0 * g_zzz_0_yyzzzz_1[i] * fi_acd_0 + g_zzz_0_yyyzzzz_1[i] * wa_y[i];

        g_yzzz_0_yyzzzzz_0[i] = 2.0 * g_zzz_0_yzzzzz_1[i] * fi_acd_0 + g_zzz_0_yyzzzzz_1[i] * wa_y[i];

        g_yzzz_0_yzzzzzz_0[i] = g_zzz_0_zzzzzz_1[i] * fi_acd_0 + g_zzz_0_yzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_zzzzzzz_0[i] = g_zzz_0_zzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 504-540 components of targeted buffer : GSK

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

    #pragma omp simd aligned(g_zz_0_xxxxxxx_0, g_zz_0_xxxxxxx_1, g_zz_0_xxxxxxy_0, g_zz_0_xxxxxxy_1, g_zz_0_xxxxxxz_0, g_zz_0_xxxxxxz_1, g_zz_0_xxxxxyy_0, g_zz_0_xxxxxyy_1, g_zz_0_xxxxxyz_0, g_zz_0_xxxxxyz_1, g_zz_0_xxxxxzz_0, g_zz_0_xxxxxzz_1, g_zz_0_xxxxyyy_0, g_zz_0_xxxxyyy_1, g_zz_0_xxxxyyz_0, g_zz_0_xxxxyyz_1, g_zz_0_xxxxyzz_0, g_zz_0_xxxxyzz_1, g_zz_0_xxxxzzz_0, g_zz_0_xxxxzzz_1, g_zz_0_xxxyyyy_0, g_zz_0_xxxyyyy_1, g_zz_0_xxxyyyz_0, g_zz_0_xxxyyyz_1, g_zz_0_xxxyyzz_0, g_zz_0_xxxyyzz_1, g_zz_0_xxxyzzz_0, g_zz_0_xxxyzzz_1, g_zz_0_xxxzzzz_0, g_zz_0_xxxzzzz_1, g_zz_0_xxyyyyy_0, g_zz_0_xxyyyyy_1, g_zz_0_xxyyyyz_0, g_zz_0_xxyyyyz_1, g_zz_0_xxyyyzz_0, g_zz_0_xxyyyzz_1, g_zz_0_xxyyzzz_0, g_zz_0_xxyyzzz_1, g_zz_0_xxyzzzz_0, g_zz_0_xxyzzzz_1, g_zz_0_xxzzzzz_0, g_zz_0_xxzzzzz_1, g_zz_0_xyyyyyy_0, g_zz_0_xyyyyyy_1, g_zz_0_xyyyyyz_0, g_zz_0_xyyyyyz_1, g_zz_0_xyyyyzz_0, g_zz_0_xyyyyzz_1, g_zz_0_xyyyzzz_0, g_zz_0_xyyyzzz_1, g_zz_0_xyyzzzz_0, g_zz_0_xyyzzzz_1, g_zz_0_xyzzzzz_0, g_zz_0_xyzzzzz_1, g_zz_0_xzzzzzz_0, g_zz_0_xzzzzzz_1, g_zz_0_yyyyyyy_0, g_zz_0_yyyyyyy_1, g_zz_0_yyyyyyz_0, g_zz_0_yyyyyyz_1, g_zz_0_yyyyyzz_0, g_zz_0_yyyyyzz_1, g_zz_0_yyyyzzz_0, g_zz_0_yyyyzzz_1, g_zz_0_yyyzzzz_0, g_zz_0_yyyzzzz_1, g_zz_0_yyzzzzz_0, g_zz_0_yyzzzzz_1, g_zz_0_yzzzzzz_0, g_zz_0_yzzzzzz_1, g_zz_0_zzzzzzz_0, g_zz_0_zzzzzzz_1, g_zzz_0_xxxxxx_1, g_zzz_0_xxxxxxx_1, g_zzz_0_xxxxxxy_1, g_zzz_0_xxxxxxz_1, g_zzz_0_xxxxxy_1, g_zzz_0_xxxxxyy_1, g_zzz_0_xxxxxyz_1, g_zzz_0_xxxxxz_1, g_zzz_0_xxxxxzz_1, g_zzz_0_xxxxyy_1, g_zzz_0_xxxxyyy_1, g_zzz_0_xxxxyyz_1, g_zzz_0_xxxxyz_1, g_zzz_0_xxxxyzz_1, g_zzz_0_xxxxzz_1, g_zzz_0_xxxxzzz_1, g_zzz_0_xxxyyy_1, g_zzz_0_xxxyyyy_1, g_zzz_0_xxxyyyz_1, g_zzz_0_xxxyyz_1, g_zzz_0_xxxyyzz_1, g_zzz_0_xxxyzz_1, g_zzz_0_xxxyzzz_1, g_zzz_0_xxxzzz_1, g_zzz_0_xxxzzzz_1, g_zzz_0_xxyyyy_1, g_zzz_0_xxyyyyy_1, g_zzz_0_xxyyyyz_1, g_zzz_0_xxyyyz_1, g_zzz_0_xxyyyzz_1, g_zzz_0_xxyyzz_1, g_zzz_0_xxyyzzz_1, g_zzz_0_xxyzzz_1, g_zzz_0_xxyzzzz_1, g_zzz_0_xxzzzz_1, g_zzz_0_xxzzzzz_1, g_zzz_0_xyyyyy_1, g_zzz_0_xyyyyyy_1, g_zzz_0_xyyyyyz_1, g_zzz_0_xyyyyz_1, g_zzz_0_xyyyyzz_1, g_zzz_0_xyyyzz_1, g_zzz_0_xyyyzzz_1, g_zzz_0_xyyzzz_1, g_zzz_0_xyyzzzz_1, g_zzz_0_xyzzzz_1, g_zzz_0_xyzzzzz_1, g_zzz_0_xzzzzz_1, g_zzz_0_xzzzzzz_1, g_zzz_0_yyyyyy_1, g_zzz_0_yyyyyyy_1, g_zzz_0_yyyyyyz_1, g_zzz_0_yyyyyz_1, g_zzz_0_yyyyyzz_1, g_zzz_0_yyyyzz_1, g_zzz_0_yyyyzzz_1, g_zzz_0_yyyzzz_1, g_zzz_0_yyyzzzz_1, g_zzz_0_yyzzzz_1, g_zzz_0_yyzzzzz_1, g_zzz_0_yzzzzz_1, g_zzz_0_yzzzzzz_1, g_zzz_0_zzzzzz_1, g_zzz_0_zzzzzzz_1, g_zzzz_0_xxxxxxx_0, g_zzzz_0_xxxxxxy_0, g_zzzz_0_xxxxxxz_0, g_zzzz_0_xxxxxyy_0, g_zzzz_0_xxxxxyz_0, g_zzzz_0_xxxxxzz_0, g_zzzz_0_xxxxyyy_0, g_zzzz_0_xxxxyyz_0, g_zzzz_0_xxxxyzz_0, g_zzzz_0_xxxxzzz_0, g_zzzz_0_xxxyyyy_0, g_zzzz_0_xxxyyyz_0, g_zzzz_0_xxxyyzz_0, g_zzzz_0_xxxyzzz_0, g_zzzz_0_xxxzzzz_0, g_zzzz_0_xxyyyyy_0, g_zzzz_0_xxyyyyz_0, g_zzzz_0_xxyyyzz_0, g_zzzz_0_xxyyzzz_0, g_zzzz_0_xxyzzzz_0, g_zzzz_0_xxzzzzz_0, g_zzzz_0_xyyyyyy_0, g_zzzz_0_xyyyyyz_0, g_zzzz_0_xyyyyzz_0, g_zzzz_0_xyyyzzz_0, g_zzzz_0_xyyzzzz_0, g_zzzz_0_xyzzzzz_0, g_zzzz_0_xzzzzzz_0, g_zzzz_0_yyyyyyy_0, g_zzzz_0_yyyyyyz_0, g_zzzz_0_yyyyyzz_0, g_zzzz_0_yyyyzzz_0, g_zzzz_0_yyyzzzz_0, g_zzzz_0_yyzzzzz_0, g_zzzz_0_yzzzzzz_0, g_zzzz_0_zzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzz_0_xxxxxxx_0[i] = 3.0 * g_zz_0_xxxxxxx_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxx_1[i] * fz_be_0 + g_zzz_0_xxxxxxx_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxy_0[i] = 3.0 * g_zz_0_xxxxxxy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxy_1[i] * fz_be_0 + g_zzz_0_xxxxxxy_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxz_0[i] = 3.0 * g_zz_0_xxxxxxz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxz_1[i] * fz_be_0 + g_zzz_0_xxxxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxxxz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxyy_0[i] = 3.0 * g_zz_0_xxxxxyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxyy_1[i] * fz_be_0 + g_zzz_0_xxxxxyy_1[i] * wa_z[i];

        g_zzzz_0_xxxxxyz_0[i] = 3.0 * g_zz_0_xxxxxyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxyz_1[i] * fz_be_0 + g_zzz_0_xxxxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxxyz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxzz_0[i] = 3.0 * g_zz_0_xxxxxzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxxzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxyyy_0[i] = 3.0 * g_zz_0_xxxxyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyyy_1[i] * fz_be_0 + g_zzz_0_xxxxyyy_1[i] * wa_z[i];

        g_zzzz_0_xxxxyyz_0[i] = 3.0 * g_zz_0_xxxxyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyyz_1[i] * fz_be_0 + g_zzz_0_xxxxyy_1[i] * fi_acd_0 + g_zzz_0_xxxxyyz_1[i] * wa_z[i];

        g_zzzz_0_xxxxyzz_0[i] = 3.0 * g_zz_0_xxxxyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxxyz_1[i] * fi_acd_0 + g_zzz_0_xxxxyzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxzzz_0[i] = 3.0 * g_zz_0_xxxxzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxxxzz_1[i] * fi_acd_0 + g_zzz_0_xxxxzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxyyyy_0[i] = 3.0 * g_zz_0_xxxyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyyy_1[i] * fz_be_0 + g_zzz_0_xxxyyyy_1[i] * wa_z[i];

        g_zzzz_0_xxxyyyz_0[i] = 3.0 * g_zz_0_xxxyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyyz_1[i] * fz_be_0 + g_zzz_0_xxxyyy_1[i] * fi_acd_0 + g_zzz_0_xxxyyyz_1[i] * wa_z[i];

        g_zzzz_0_xxxyyzz_0[i] = 3.0 * g_zz_0_xxxyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxyyz_1[i] * fi_acd_0 + g_zzz_0_xxxyyzz_1[i] * wa_z[i];

        g_zzzz_0_xxxyzzz_0[i] = 3.0 * g_zz_0_xxxyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxxyzz_1[i] * fi_acd_0 + g_zzz_0_xxxyzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxzzzz_0[i] = 3.0 * g_zz_0_xxxzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xxxzzz_1[i] * fi_acd_0 + g_zzz_0_xxxzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxyyyyy_0[i] = 3.0 * g_zz_0_xxyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyyy_1[i] * fz_be_0 + g_zzz_0_xxyyyyy_1[i] * wa_z[i];

        g_zzzz_0_xxyyyyz_0[i] = 3.0 * g_zz_0_xxyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyyz_1[i] * fz_be_0 + g_zzz_0_xxyyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyyyz_1[i] * wa_z[i];

        g_zzzz_0_xxyyyzz_0[i] = 3.0 * g_zz_0_xxyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxyyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyyzz_1[i] * wa_z[i];

        g_zzzz_0_xxyyzzz_0[i] = 3.0 * g_zz_0_xxyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxyyzz_1[i] * fi_acd_0 + g_zzz_0_xxyyzzz_1[i] * wa_z[i];

        g_zzzz_0_xxyzzzz_0[i] = 3.0 * g_zz_0_xxyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xxyzzz_1[i] * fi_acd_0 + g_zzz_0_xxyzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxzzzzz_0[i] = 3.0 * g_zz_0_xxzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_xxzzzz_1[i] * fi_acd_0 + g_zzz_0_xxzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xyyyyyy_0[i] = 3.0 * g_zz_0_xyyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyyy_1[i] * fz_be_0 + g_zzz_0_xyyyyyy_1[i] * wa_z[i];

        g_zzzz_0_xyyyyyz_0[i] = 3.0 * g_zz_0_xyyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyyz_1[i] * fz_be_0 + g_zzz_0_xyyyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyyyz_1[i] * wa_z[i];

        g_zzzz_0_xyyyyzz_0[i] = 3.0 * g_zz_0_xyyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xyyyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyyzz_1[i] * wa_z[i];

        g_zzzz_0_xyyyzzz_0[i] = 3.0 * g_zz_0_xyyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xyyyzz_1[i] * fi_acd_0 + g_zzz_0_xyyyzzz_1[i] * wa_z[i];

        g_zzzz_0_xyyzzzz_0[i] = 3.0 * g_zz_0_xyyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xyyzzz_1[i] * fi_acd_0 + g_zzz_0_xyyzzzz_1[i] * wa_z[i];

        g_zzzz_0_xyzzzzz_0[i] = 3.0 * g_zz_0_xyzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_xyzzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xzzzzzz_0[i] = 3.0 * g_zz_0_xzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzz_0_xzzzzz_1[i] * fi_acd_0 + g_zzz_0_xzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_yyyyyyy_0[i] = 3.0 * g_zz_0_yyyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyyy_1[i] * fz_be_0 + g_zzz_0_yyyyyyy_1[i] * wa_z[i];

        g_zzzz_0_yyyyyyz_0[i] = 3.0 * g_zz_0_yyyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyyz_1[i] * fz_be_0 + g_zzz_0_yyyyyy_1[i] * fi_acd_0 + g_zzz_0_yyyyyyz_1[i] * wa_z[i];

        g_zzzz_0_yyyyyzz_0[i] = 3.0 * g_zz_0_yyyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_yyyyyz_1[i] * fi_acd_0 + g_zzz_0_yyyyyzz_1[i] * wa_z[i];

        g_zzzz_0_yyyyzzz_0[i] = 3.0 * g_zz_0_yyyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_yyyyzz_1[i] * fi_acd_0 + g_zzz_0_yyyyzzz_1[i] * wa_z[i];

        g_zzzz_0_yyyzzzz_0[i] = 3.0 * g_zz_0_yyyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_yyyzzz_1[i] * fi_acd_0 + g_zzz_0_yyyzzzz_1[i] * wa_z[i];

        g_zzzz_0_yyzzzzz_0[i] = 3.0 * g_zz_0_yyzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_yyzzzz_1[i] * fi_acd_0 + g_zzz_0_yyzzzzz_1[i] * wa_z[i];

        g_zzzz_0_yzzzzzz_0[i] = 3.0 * g_zz_0_yzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzz_0_yzzzzz_1[i] * fi_acd_0 + g_zzz_0_yzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_zzzzzzz_0[i] = 3.0 * g_zz_0_zzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_zzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzz_0_zzzzzz_1[i] * fi_acd_0 + g_zzz_0_zzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

