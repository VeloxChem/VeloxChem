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

#include "TwoCenterElectronRepulsionPrimRecGI.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_gi(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_gi,
                                const size_t idx_eri_0_di,
                                const size_t idx_eri_1_di,
                                const size_t idx_eri_1_fh,
                                const size_t idx_eri_1_fi,
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

    // Set up components of auxiliary buffer : DI

    auto g_xx_xxxxxx_0 = pbuffer.data(idx_eri_0_di);

    auto g_xx_xxxxxy_0 = pbuffer.data(idx_eri_0_di + 1);

    auto g_xx_xxxxxz_0 = pbuffer.data(idx_eri_0_di + 2);

    auto g_xx_xxxxyy_0 = pbuffer.data(idx_eri_0_di + 3);

    auto g_xx_xxxxyz_0 = pbuffer.data(idx_eri_0_di + 4);

    auto g_xx_xxxxzz_0 = pbuffer.data(idx_eri_0_di + 5);

    auto g_xx_xxxyyy_0 = pbuffer.data(idx_eri_0_di + 6);

    auto g_xx_xxxyyz_0 = pbuffer.data(idx_eri_0_di + 7);

    auto g_xx_xxxyzz_0 = pbuffer.data(idx_eri_0_di + 8);

    auto g_xx_xxxzzz_0 = pbuffer.data(idx_eri_0_di + 9);

    auto g_xx_xxyyyy_0 = pbuffer.data(idx_eri_0_di + 10);

    auto g_xx_xxyyyz_0 = pbuffer.data(idx_eri_0_di + 11);

    auto g_xx_xxyyzz_0 = pbuffer.data(idx_eri_0_di + 12);

    auto g_xx_xxyzzz_0 = pbuffer.data(idx_eri_0_di + 13);

    auto g_xx_xxzzzz_0 = pbuffer.data(idx_eri_0_di + 14);

    auto g_xx_xyyyyy_0 = pbuffer.data(idx_eri_0_di + 15);

    auto g_xx_xyyyyz_0 = pbuffer.data(idx_eri_0_di + 16);

    auto g_xx_xyyyzz_0 = pbuffer.data(idx_eri_0_di + 17);

    auto g_xx_xyyzzz_0 = pbuffer.data(idx_eri_0_di + 18);

    auto g_xx_xyzzzz_0 = pbuffer.data(idx_eri_0_di + 19);

    auto g_xx_xzzzzz_0 = pbuffer.data(idx_eri_0_di + 20);

    auto g_xx_yyyyyy_0 = pbuffer.data(idx_eri_0_di + 21);

    auto g_xx_yyyyyz_0 = pbuffer.data(idx_eri_0_di + 22);

    auto g_xx_yyyyzz_0 = pbuffer.data(idx_eri_0_di + 23);

    auto g_xx_yyyzzz_0 = pbuffer.data(idx_eri_0_di + 24);

    auto g_xx_yyzzzz_0 = pbuffer.data(idx_eri_0_di + 25);

    auto g_xx_yzzzzz_0 = pbuffer.data(idx_eri_0_di + 26);

    auto g_xx_zzzzzz_0 = pbuffer.data(idx_eri_0_di + 27);

    auto g_yy_xxxxxx_0 = pbuffer.data(idx_eri_0_di + 84);

    auto g_yy_xxxxxy_0 = pbuffer.data(idx_eri_0_di + 85);

    auto g_yy_xxxxxz_0 = pbuffer.data(idx_eri_0_di + 86);

    auto g_yy_xxxxyy_0 = pbuffer.data(idx_eri_0_di + 87);

    auto g_yy_xxxxyz_0 = pbuffer.data(idx_eri_0_di + 88);

    auto g_yy_xxxxzz_0 = pbuffer.data(idx_eri_0_di + 89);

    auto g_yy_xxxyyy_0 = pbuffer.data(idx_eri_0_di + 90);

    auto g_yy_xxxyyz_0 = pbuffer.data(idx_eri_0_di + 91);

    auto g_yy_xxxyzz_0 = pbuffer.data(idx_eri_0_di + 92);

    auto g_yy_xxxzzz_0 = pbuffer.data(idx_eri_0_di + 93);

    auto g_yy_xxyyyy_0 = pbuffer.data(idx_eri_0_di + 94);

    auto g_yy_xxyyyz_0 = pbuffer.data(idx_eri_0_di + 95);

    auto g_yy_xxyyzz_0 = pbuffer.data(idx_eri_0_di + 96);

    auto g_yy_xxyzzz_0 = pbuffer.data(idx_eri_0_di + 97);

    auto g_yy_xxzzzz_0 = pbuffer.data(idx_eri_0_di + 98);

    auto g_yy_xyyyyy_0 = pbuffer.data(idx_eri_0_di + 99);

    auto g_yy_xyyyyz_0 = pbuffer.data(idx_eri_0_di + 100);

    auto g_yy_xyyyzz_0 = pbuffer.data(idx_eri_0_di + 101);

    auto g_yy_xyyzzz_0 = pbuffer.data(idx_eri_0_di + 102);

    auto g_yy_xyzzzz_0 = pbuffer.data(idx_eri_0_di + 103);

    auto g_yy_xzzzzz_0 = pbuffer.data(idx_eri_0_di + 104);

    auto g_yy_yyyyyy_0 = pbuffer.data(idx_eri_0_di + 105);

    auto g_yy_yyyyyz_0 = pbuffer.data(idx_eri_0_di + 106);

    auto g_yy_yyyyzz_0 = pbuffer.data(idx_eri_0_di + 107);

    auto g_yy_yyyzzz_0 = pbuffer.data(idx_eri_0_di + 108);

    auto g_yy_yyzzzz_0 = pbuffer.data(idx_eri_0_di + 109);

    auto g_yy_yzzzzz_0 = pbuffer.data(idx_eri_0_di + 110);

    auto g_yy_zzzzzz_0 = pbuffer.data(idx_eri_0_di + 111);

    auto g_zz_xxxxxx_0 = pbuffer.data(idx_eri_0_di + 140);

    auto g_zz_xxxxxy_0 = pbuffer.data(idx_eri_0_di + 141);

    auto g_zz_xxxxxz_0 = pbuffer.data(idx_eri_0_di + 142);

    auto g_zz_xxxxyy_0 = pbuffer.data(idx_eri_0_di + 143);

    auto g_zz_xxxxyz_0 = pbuffer.data(idx_eri_0_di + 144);

    auto g_zz_xxxxzz_0 = pbuffer.data(idx_eri_0_di + 145);

    auto g_zz_xxxyyy_0 = pbuffer.data(idx_eri_0_di + 146);

    auto g_zz_xxxyyz_0 = pbuffer.data(idx_eri_0_di + 147);

    auto g_zz_xxxyzz_0 = pbuffer.data(idx_eri_0_di + 148);

    auto g_zz_xxxzzz_0 = pbuffer.data(idx_eri_0_di + 149);

    auto g_zz_xxyyyy_0 = pbuffer.data(idx_eri_0_di + 150);

    auto g_zz_xxyyyz_0 = pbuffer.data(idx_eri_0_di + 151);

    auto g_zz_xxyyzz_0 = pbuffer.data(idx_eri_0_di + 152);

    auto g_zz_xxyzzz_0 = pbuffer.data(idx_eri_0_di + 153);

    auto g_zz_xxzzzz_0 = pbuffer.data(idx_eri_0_di + 154);

    auto g_zz_xyyyyy_0 = pbuffer.data(idx_eri_0_di + 155);

    auto g_zz_xyyyyz_0 = pbuffer.data(idx_eri_0_di + 156);

    auto g_zz_xyyyzz_0 = pbuffer.data(idx_eri_0_di + 157);

    auto g_zz_xyyzzz_0 = pbuffer.data(idx_eri_0_di + 158);

    auto g_zz_xyzzzz_0 = pbuffer.data(idx_eri_0_di + 159);

    auto g_zz_xzzzzz_0 = pbuffer.data(idx_eri_0_di + 160);

    auto g_zz_yyyyyy_0 = pbuffer.data(idx_eri_0_di + 161);

    auto g_zz_yyyyyz_0 = pbuffer.data(idx_eri_0_di + 162);

    auto g_zz_yyyyzz_0 = pbuffer.data(idx_eri_0_di + 163);

    auto g_zz_yyyzzz_0 = pbuffer.data(idx_eri_0_di + 164);

    auto g_zz_yyzzzz_0 = pbuffer.data(idx_eri_0_di + 165);

    auto g_zz_yzzzzz_0 = pbuffer.data(idx_eri_0_di + 166);

    auto g_zz_zzzzzz_0 = pbuffer.data(idx_eri_0_di + 167);

    // Set up components of auxiliary buffer : DI

    auto g_xx_xxxxxx_1 = pbuffer.data(idx_eri_1_di);

    auto g_xx_xxxxxy_1 = pbuffer.data(idx_eri_1_di + 1);

    auto g_xx_xxxxxz_1 = pbuffer.data(idx_eri_1_di + 2);

    auto g_xx_xxxxyy_1 = pbuffer.data(idx_eri_1_di + 3);

    auto g_xx_xxxxyz_1 = pbuffer.data(idx_eri_1_di + 4);

    auto g_xx_xxxxzz_1 = pbuffer.data(idx_eri_1_di + 5);

    auto g_xx_xxxyyy_1 = pbuffer.data(idx_eri_1_di + 6);

    auto g_xx_xxxyyz_1 = pbuffer.data(idx_eri_1_di + 7);

    auto g_xx_xxxyzz_1 = pbuffer.data(idx_eri_1_di + 8);

    auto g_xx_xxxzzz_1 = pbuffer.data(idx_eri_1_di + 9);

    auto g_xx_xxyyyy_1 = pbuffer.data(idx_eri_1_di + 10);

    auto g_xx_xxyyyz_1 = pbuffer.data(idx_eri_1_di + 11);

    auto g_xx_xxyyzz_1 = pbuffer.data(idx_eri_1_di + 12);

    auto g_xx_xxyzzz_1 = pbuffer.data(idx_eri_1_di + 13);

    auto g_xx_xxzzzz_1 = pbuffer.data(idx_eri_1_di + 14);

    auto g_xx_xyyyyy_1 = pbuffer.data(idx_eri_1_di + 15);

    auto g_xx_xyyyyz_1 = pbuffer.data(idx_eri_1_di + 16);

    auto g_xx_xyyyzz_1 = pbuffer.data(idx_eri_1_di + 17);

    auto g_xx_xyyzzz_1 = pbuffer.data(idx_eri_1_di + 18);

    auto g_xx_xyzzzz_1 = pbuffer.data(idx_eri_1_di + 19);

    auto g_xx_xzzzzz_1 = pbuffer.data(idx_eri_1_di + 20);

    auto g_xx_yyyyyy_1 = pbuffer.data(idx_eri_1_di + 21);

    auto g_xx_yyyyyz_1 = pbuffer.data(idx_eri_1_di + 22);

    auto g_xx_yyyyzz_1 = pbuffer.data(idx_eri_1_di + 23);

    auto g_xx_yyyzzz_1 = pbuffer.data(idx_eri_1_di + 24);

    auto g_xx_yyzzzz_1 = pbuffer.data(idx_eri_1_di + 25);

    auto g_xx_yzzzzz_1 = pbuffer.data(idx_eri_1_di + 26);

    auto g_xx_zzzzzz_1 = pbuffer.data(idx_eri_1_di + 27);

    auto g_yy_xxxxxx_1 = pbuffer.data(idx_eri_1_di + 84);

    auto g_yy_xxxxxy_1 = pbuffer.data(idx_eri_1_di + 85);

    auto g_yy_xxxxxz_1 = pbuffer.data(idx_eri_1_di + 86);

    auto g_yy_xxxxyy_1 = pbuffer.data(idx_eri_1_di + 87);

    auto g_yy_xxxxyz_1 = pbuffer.data(idx_eri_1_di + 88);

    auto g_yy_xxxxzz_1 = pbuffer.data(idx_eri_1_di + 89);

    auto g_yy_xxxyyy_1 = pbuffer.data(idx_eri_1_di + 90);

    auto g_yy_xxxyyz_1 = pbuffer.data(idx_eri_1_di + 91);

    auto g_yy_xxxyzz_1 = pbuffer.data(idx_eri_1_di + 92);

    auto g_yy_xxxzzz_1 = pbuffer.data(idx_eri_1_di + 93);

    auto g_yy_xxyyyy_1 = pbuffer.data(idx_eri_1_di + 94);

    auto g_yy_xxyyyz_1 = pbuffer.data(idx_eri_1_di + 95);

    auto g_yy_xxyyzz_1 = pbuffer.data(idx_eri_1_di + 96);

    auto g_yy_xxyzzz_1 = pbuffer.data(idx_eri_1_di + 97);

    auto g_yy_xxzzzz_1 = pbuffer.data(idx_eri_1_di + 98);

    auto g_yy_xyyyyy_1 = pbuffer.data(idx_eri_1_di + 99);

    auto g_yy_xyyyyz_1 = pbuffer.data(idx_eri_1_di + 100);

    auto g_yy_xyyyzz_1 = pbuffer.data(idx_eri_1_di + 101);

    auto g_yy_xyyzzz_1 = pbuffer.data(idx_eri_1_di + 102);

    auto g_yy_xyzzzz_1 = pbuffer.data(idx_eri_1_di + 103);

    auto g_yy_xzzzzz_1 = pbuffer.data(idx_eri_1_di + 104);

    auto g_yy_yyyyyy_1 = pbuffer.data(idx_eri_1_di + 105);

    auto g_yy_yyyyyz_1 = pbuffer.data(idx_eri_1_di + 106);

    auto g_yy_yyyyzz_1 = pbuffer.data(idx_eri_1_di + 107);

    auto g_yy_yyyzzz_1 = pbuffer.data(idx_eri_1_di + 108);

    auto g_yy_yyzzzz_1 = pbuffer.data(idx_eri_1_di + 109);

    auto g_yy_yzzzzz_1 = pbuffer.data(idx_eri_1_di + 110);

    auto g_yy_zzzzzz_1 = pbuffer.data(idx_eri_1_di + 111);

    auto g_zz_xxxxxx_1 = pbuffer.data(idx_eri_1_di + 140);

    auto g_zz_xxxxxy_1 = pbuffer.data(idx_eri_1_di + 141);

    auto g_zz_xxxxxz_1 = pbuffer.data(idx_eri_1_di + 142);

    auto g_zz_xxxxyy_1 = pbuffer.data(idx_eri_1_di + 143);

    auto g_zz_xxxxyz_1 = pbuffer.data(idx_eri_1_di + 144);

    auto g_zz_xxxxzz_1 = pbuffer.data(idx_eri_1_di + 145);

    auto g_zz_xxxyyy_1 = pbuffer.data(idx_eri_1_di + 146);

    auto g_zz_xxxyyz_1 = pbuffer.data(idx_eri_1_di + 147);

    auto g_zz_xxxyzz_1 = pbuffer.data(idx_eri_1_di + 148);

    auto g_zz_xxxzzz_1 = pbuffer.data(idx_eri_1_di + 149);

    auto g_zz_xxyyyy_1 = pbuffer.data(idx_eri_1_di + 150);

    auto g_zz_xxyyyz_1 = pbuffer.data(idx_eri_1_di + 151);

    auto g_zz_xxyyzz_1 = pbuffer.data(idx_eri_1_di + 152);

    auto g_zz_xxyzzz_1 = pbuffer.data(idx_eri_1_di + 153);

    auto g_zz_xxzzzz_1 = pbuffer.data(idx_eri_1_di + 154);

    auto g_zz_xyyyyy_1 = pbuffer.data(idx_eri_1_di + 155);

    auto g_zz_xyyyyz_1 = pbuffer.data(idx_eri_1_di + 156);

    auto g_zz_xyyyzz_1 = pbuffer.data(idx_eri_1_di + 157);

    auto g_zz_xyyzzz_1 = pbuffer.data(idx_eri_1_di + 158);

    auto g_zz_xyzzzz_1 = pbuffer.data(idx_eri_1_di + 159);

    auto g_zz_xzzzzz_1 = pbuffer.data(idx_eri_1_di + 160);

    auto g_zz_yyyyyy_1 = pbuffer.data(idx_eri_1_di + 161);

    auto g_zz_yyyyyz_1 = pbuffer.data(idx_eri_1_di + 162);

    auto g_zz_yyyyzz_1 = pbuffer.data(idx_eri_1_di + 163);

    auto g_zz_yyyzzz_1 = pbuffer.data(idx_eri_1_di + 164);

    auto g_zz_yyzzzz_1 = pbuffer.data(idx_eri_1_di + 165);

    auto g_zz_yzzzzz_1 = pbuffer.data(idx_eri_1_di + 166);

    auto g_zz_zzzzzz_1 = pbuffer.data(idx_eri_1_di + 167);

    // Set up components of auxiliary buffer : FH

    auto g_xxx_xxxxx_1 = pbuffer.data(idx_eri_1_fh);

    auto g_xxx_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 1);

    auto g_xxx_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 2);

    auto g_xxx_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 3);

    auto g_xxx_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 4);

    auto g_xxx_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 5);

    auto g_xxx_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 6);

    auto g_xxx_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 7);

    auto g_xxx_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 8);

    auto g_xxx_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 9);

    auto g_xxx_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 10);

    auto g_xxx_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 11);

    auto g_xxx_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 12);

    auto g_xxx_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 13);

    auto g_xxx_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 14);

    auto g_xxx_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 15);

    auto g_xxx_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 16);

    auto g_xxx_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 17);

    auto g_xxx_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 18);

    auto g_xxx_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 19);

    auto g_xxx_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 20);

    auto g_xxz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 44);

    auto g_xxz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 46);

    auto g_xxz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 47);

    auto g_xxz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 49);

    auto g_xxz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 50);

    auto g_xxz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 51);

    auto g_xxz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 53);

    auto g_xxz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 54);

    auto g_xxz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 55);

    auto g_xxz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 56);

    auto g_xxz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 58);

    auto g_xxz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 59);

    auto g_xxz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 60);

    auto g_xxz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 61);

    auto g_xxz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 62);

    auto g_xyy_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 64);

    auto g_xyy_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 66);

    auto g_xyy_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 67);

    auto g_xyy_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 69);

    auto g_xyy_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 70);

    auto g_xyy_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 71);

    auto g_xyy_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 73);

    auto g_xyy_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 74);

    auto g_xyy_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 75);

    auto g_xyy_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 76);

    auto g_xyy_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 78);

    auto g_xyy_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 79);

    auto g_xyy_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 80);

    auto g_xyy_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 81);

    auto g_xyy_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 82);

    auto g_xzz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 107);

    auto g_xzz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 109);

    auto g_xzz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 110);

    auto g_xzz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 112);

    auto g_xzz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 113);

    auto g_xzz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 114);

    auto g_xzz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 116);

    auto g_xzz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 117);

    auto g_xzz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 118);

    auto g_xzz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 119);

    auto g_xzz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 121);

    auto g_xzz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 122);

    auto g_xzz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 123);

    auto g_xzz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 124);

    auto g_xzz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 125);

    auto g_yyy_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 126);

    auto g_yyy_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 127);

    auto g_yyy_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 128);

    auto g_yyy_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 129);

    auto g_yyy_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 130);

    auto g_yyy_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 131);

    auto g_yyy_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 132);

    auto g_yyy_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 133);

    auto g_yyy_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 134);

    auto g_yyy_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 135);

    auto g_yyy_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 136);

    auto g_yyy_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 137);

    auto g_yyy_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 138);

    auto g_yyy_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 139);

    auto g_yyy_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 140);

    auto g_yyy_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 141);

    auto g_yyy_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 142);

    auto g_yyy_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 143);

    auto g_yyy_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 144);

    auto g_yyy_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 145);

    auto g_yyy_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 146);

    auto g_yyz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 149);

    auto g_yyz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 151);

    auto g_yyz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 152);

    auto g_yyz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 154);

    auto g_yyz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 155);

    auto g_yyz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 156);

    auto g_yyz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 158);

    auto g_yyz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 159);

    auto g_yyz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 160);

    auto g_yyz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 161);

    auto g_yyz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 163);

    auto g_yyz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 164);

    auto g_yyz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 165);

    auto g_yyz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 166);

    auto g_yyz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 167);

    auto g_yzz_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 169);

    auto g_yzz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 170);

    auto g_yzz_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 171);

    auto g_yzz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 172);

    auto g_yzz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 173);

    auto g_yzz_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 174);

    auto g_yzz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 175);

    auto g_yzz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 176);

    auto g_yzz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 177);

    auto g_yzz_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 178);

    auto g_yzz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 179);

    auto g_yzz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 180);

    auto g_yzz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 181);

    auto g_yzz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 182);

    auto g_yzz_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 183);

    auto g_yzz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 184);

    auto g_yzz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 185);

    auto g_yzz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 186);

    auto g_yzz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 187);

    auto g_yzz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 188);

    auto g_zzz_xxxxx_1 = pbuffer.data(idx_eri_1_fh + 189);

    auto g_zzz_xxxxy_1 = pbuffer.data(idx_eri_1_fh + 190);

    auto g_zzz_xxxxz_1 = pbuffer.data(idx_eri_1_fh + 191);

    auto g_zzz_xxxyy_1 = pbuffer.data(idx_eri_1_fh + 192);

    auto g_zzz_xxxyz_1 = pbuffer.data(idx_eri_1_fh + 193);

    auto g_zzz_xxxzz_1 = pbuffer.data(idx_eri_1_fh + 194);

    auto g_zzz_xxyyy_1 = pbuffer.data(idx_eri_1_fh + 195);

    auto g_zzz_xxyyz_1 = pbuffer.data(idx_eri_1_fh + 196);

    auto g_zzz_xxyzz_1 = pbuffer.data(idx_eri_1_fh + 197);

    auto g_zzz_xxzzz_1 = pbuffer.data(idx_eri_1_fh + 198);

    auto g_zzz_xyyyy_1 = pbuffer.data(idx_eri_1_fh + 199);

    auto g_zzz_xyyyz_1 = pbuffer.data(idx_eri_1_fh + 200);

    auto g_zzz_xyyzz_1 = pbuffer.data(idx_eri_1_fh + 201);

    auto g_zzz_xyzzz_1 = pbuffer.data(idx_eri_1_fh + 202);

    auto g_zzz_xzzzz_1 = pbuffer.data(idx_eri_1_fh + 203);

    auto g_zzz_yyyyy_1 = pbuffer.data(idx_eri_1_fh + 204);

    auto g_zzz_yyyyz_1 = pbuffer.data(idx_eri_1_fh + 205);

    auto g_zzz_yyyzz_1 = pbuffer.data(idx_eri_1_fh + 206);

    auto g_zzz_yyzzz_1 = pbuffer.data(idx_eri_1_fh + 207);

    auto g_zzz_yzzzz_1 = pbuffer.data(idx_eri_1_fh + 208);

    auto g_zzz_zzzzz_1 = pbuffer.data(idx_eri_1_fh + 209);

    // Set up components of auxiliary buffer : FI

    auto g_xxx_xxxxxx_1 = pbuffer.data(idx_eri_1_fi);

    auto g_xxx_xxxxxy_1 = pbuffer.data(idx_eri_1_fi + 1);

    auto g_xxx_xxxxxz_1 = pbuffer.data(idx_eri_1_fi + 2);

    auto g_xxx_xxxxyy_1 = pbuffer.data(idx_eri_1_fi + 3);

    auto g_xxx_xxxxyz_1 = pbuffer.data(idx_eri_1_fi + 4);

    auto g_xxx_xxxxzz_1 = pbuffer.data(idx_eri_1_fi + 5);

    auto g_xxx_xxxyyy_1 = pbuffer.data(idx_eri_1_fi + 6);

    auto g_xxx_xxxyyz_1 = pbuffer.data(idx_eri_1_fi + 7);

    auto g_xxx_xxxyzz_1 = pbuffer.data(idx_eri_1_fi + 8);

    auto g_xxx_xxxzzz_1 = pbuffer.data(idx_eri_1_fi + 9);

    auto g_xxx_xxyyyy_1 = pbuffer.data(idx_eri_1_fi + 10);

    auto g_xxx_xxyyyz_1 = pbuffer.data(idx_eri_1_fi + 11);

    auto g_xxx_xxyyzz_1 = pbuffer.data(idx_eri_1_fi + 12);

    auto g_xxx_xxyzzz_1 = pbuffer.data(idx_eri_1_fi + 13);

    auto g_xxx_xxzzzz_1 = pbuffer.data(idx_eri_1_fi + 14);

    auto g_xxx_xyyyyy_1 = pbuffer.data(idx_eri_1_fi + 15);

    auto g_xxx_xyyyyz_1 = pbuffer.data(idx_eri_1_fi + 16);

    auto g_xxx_xyyyzz_1 = pbuffer.data(idx_eri_1_fi + 17);

    auto g_xxx_xyyzzz_1 = pbuffer.data(idx_eri_1_fi + 18);

    auto g_xxx_xyzzzz_1 = pbuffer.data(idx_eri_1_fi + 19);

    auto g_xxx_xzzzzz_1 = pbuffer.data(idx_eri_1_fi + 20);

    auto g_xxx_yyyyyy_1 = pbuffer.data(idx_eri_1_fi + 21);

    auto g_xxx_yyyyyz_1 = pbuffer.data(idx_eri_1_fi + 22);

    auto g_xxx_yyyyzz_1 = pbuffer.data(idx_eri_1_fi + 23);

    auto g_xxx_yyyzzz_1 = pbuffer.data(idx_eri_1_fi + 24);

    auto g_xxx_yyzzzz_1 = pbuffer.data(idx_eri_1_fi + 25);

    auto g_xxx_yzzzzz_1 = pbuffer.data(idx_eri_1_fi + 26);

    auto g_xxx_zzzzzz_1 = pbuffer.data(idx_eri_1_fi + 27);

    auto g_xxy_xxxxxx_1 = pbuffer.data(idx_eri_1_fi + 28);

    auto g_xxy_xxxxxy_1 = pbuffer.data(idx_eri_1_fi + 29);

    auto g_xxy_xxxxxz_1 = pbuffer.data(idx_eri_1_fi + 30);

    auto g_xxy_xxxxyy_1 = pbuffer.data(idx_eri_1_fi + 31);

    auto g_xxy_xxxxzz_1 = pbuffer.data(idx_eri_1_fi + 33);

    auto g_xxy_xxxyyy_1 = pbuffer.data(idx_eri_1_fi + 34);

    auto g_xxy_xxxzzz_1 = pbuffer.data(idx_eri_1_fi + 37);

    auto g_xxy_xxyyyy_1 = pbuffer.data(idx_eri_1_fi + 38);

    auto g_xxy_xxzzzz_1 = pbuffer.data(idx_eri_1_fi + 42);

    auto g_xxy_xyyyyy_1 = pbuffer.data(idx_eri_1_fi + 43);

    auto g_xxy_xzzzzz_1 = pbuffer.data(idx_eri_1_fi + 48);

    auto g_xxy_yyyyyy_1 = pbuffer.data(idx_eri_1_fi + 49);

    auto g_xxz_xxxxxx_1 = pbuffer.data(idx_eri_1_fi + 56);

    auto g_xxz_xxxxxy_1 = pbuffer.data(idx_eri_1_fi + 57);

    auto g_xxz_xxxxxz_1 = pbuffer.data(idx_eri_1_fi + 58);

    auto g_xxz_xxxxyy_1 = pbuffer.data(idx_eri_1_fi + 59);

    auto g_xxz_xxxxyz_1 = pbuffer.data(idx_eri_1_fi + 60);

    auto g_xxz_xxxxzz_1 = pbuffer.data(idx_eri_1_fi + 61);

    auto g_xxz_xxxyyy_1 = pbuffer.data(idx_eri_1_fi + 62);

    auto g_xxz_xxxyyz_1 = pbuffer.data(idx_eri_1_fi + 63);

    auto g_xxz_xxxyzz_1 = pbuffer.data(idx_eri_1_fi + 64);

    auto g_xxz_xxxzzz_1 = pbuffer.data(idx_eri_1_fi + 65);

    auto g_xxz_xxyyyy_1 = pbuffer.data(idx_eri_1_fi + 66);

    auto g_xxz_xxyyyz_1 = pbuffer.data(idx_eri_1_fi + 67);

    auto g_xxz_xxyyzz_1 = pbuffer.data(idx_eri_1_fi + 68);

    auto g_xxz_xxyzzz_1 = pbuffer.data(idx_eri_1_fi + 69);

    auto g_xxz_xxzzzz_1 = pbuffer.data(idx_eri_1_fi + 70);

    auto g_xxz_xyyyyy_1 = pbuffer.data(idx_eri_1_fi + 71);

    auto g_xxz_xyyyyz_1 = pbuffer.data(idx_eri_1_fi + 72);

    auto g_xxz_xyyyzz_1 = pbuffer.data(idx_eri_1_fi + 73);

    auto g_xxz_xyyzzz_1 = pbuffer.data(idx_eri_1_fi + 74);

    auto g_xxz_xyzzzz_1 = pbuffer.data(idx_eri_1_fi + 75);

    auto g_xxz_xzzzzz_1 = pbuffer.data(idx_eri_1_fi + 76);

    auto g_xxz_yyyyyz_1 = pbuffer.data(idx_eri_1_fi + 78);

    auto g_xxz_yyyyzz_1 = pbuffer.data(idx_eri_1_fi + 79);

    auto g_xxz_yyyzzz_1 = pbuffer.data(idx_eri_1_fi + 80);

    auto g_xxz_yyzzzz_1 = pbuffer.data(idx_eri_1_fi + 81);

    auto g_xxz_yzzzzz_1 = pbuffer.data(idx_eri_1_fi + 82);

    auto g_xxz_zzzzzz_1 = pbuffer.data(idx_eri_1_fi + 83);

    auto g_xyy_xxxxxx_1 = pbuffer.data(idx_eri_1_fi + 84);

    auto g_xyy_xxxxxy_1 = pbuffer.data(idx_eri_1_fi + 85);

    auto g_xyy_xxxxyy_1 = pbuffer.data(idx_eri_1_fi + 87);

    auto g_xyy_xxxxyz_1 = pbuffer.data(idx_eri_1_fi + 88);

    auto g_xyy_xxxyyy_1 = pbuffer.data(idx_eri_1_fi + 90);

    auto g_xyy_xxxyyz_1 = pbuffer.data(idx_eri_1_fi + 91);

    auto g_xyy_xxxyzz_1 = pbuffer.data(idx_eri_1_fi + 92);

    auto g_xyy_xxyyyy_1 = pbuffer.data(idx_eri_1_fi + 94);

    auto g_xyy_xxyyyz_1 = pbuffer.data(idx_eri_1_fi + 95);

    auto g_xyy_xxyyzz_1 = pbuffer.data(idx_eri_1_fi + 96);

    auto g_xyy_xxyzzz_1 = pbuffer.data(idx_eri_1_fi + 97);

    auto g_xyy_xyyyyy_1 = pbuffer.data(idx_eri_1_fi + 99);

    auto g_xyy_xyyyyz_1 = pbuffer.data(idx_eri_1_fi + 100);

    auto g_xyy_xyyyzz_1 = pbuffer.data(idx_eri_1_fi + 101);

    auto g_xyy_xyyzzz_1 = pbuffer.data(idx_eri_1_fi + 102);

    auto g_xyy_xyzzzz_1 = pbuffer.data(idx_eri_1_fi + 103);

    auto g_xyy_yyyyyy_1 = pbuffer.data(idx_eri_1_fi + 105);

    auto g_xyy_yyyyyz_1 = pbuffer.data(idx_eri_1_fi + 106);

    auto g_xyy_yyyyzz_1 = pbuffer.data(idx_eri_1_fi + 107);

    auto g_xyy_yyyzzz_1 = pbuffer.data(idx_eri_1_fi + 108);

    auto g_xyy_yyzzzz_1 = pbuffer.data(idx_eri_1_fi + 109);

    auto g_xyy_yzzzzz_1 = pbuffer.data(idx_eri_1_fi + 110);

    auto g_xyy_zzzzzz_1 = pbuffer.data(idx_eri_1_fi + 111);

    auto g_xzz_xxxxxx_1 = pbuffer.data(idx_eri_1_fi + 140);

    auto g_xzz_xxxxxz_1 = pbuffer.data(idx_eri_1_fi + 142);

    auto g_xzz_xxxxyz_1 = pbuffer.data(idx_eri_1_fi + 144);

    auto g_xzz_xxxxzz_1 = pbuffer.data(idx_eri_1_fi + 145);

    auto g_xzz_xxxyyz_1 = pbuffer.data(idx_eri_1_fi + 147);

    auto g_xzz_xxxyzz_1 = pbuffer.data(idx_eri_1_fi + 148);

    auto g_xzz_xxxzzz_1 = pbuffer.data(idx_eri_1_fi + 149);

    auto g_xzz_xxyyyz_1 = pbuffer.data(idx_eri_1_fi + 151);

    auto g_xzz_xxyyzz_1 = pbuffer.data(idx_eri_1_fi + 152);

    auto g_xzz_xxyzzz_1 = pbuffer.data(idx_eri_1_fi + 153);

    auto g_xzz_xxzzzz_1 = pbuffer.data(idx_eri_1_fi + 154);

    auto g_xzz_xyyyyz_1 = pbuffer.data(idx_eri_1_fi + 156);

    auto g_xzz_xyyyzz_1 = pbuffer.data(idx_eri_1_fi + 157);

    auto g_xzz_xyyzzz_1 = pbuffer.data(idx_eri_1_fi + 158);

    auto g_xzz_xyzzzz_1 = pbuffer.data(idx_eri_1_fi + 159);

    auto g_xzz_xzzzzz_1 = pbuffer.data(idx_eri_1_fi + 160);

    auto g_xzz_yyyyyy_1 = pbuffer.data(idx_eri_1_fi + 161);

    auto g_xzz_yyyyyz_1 = pbuffer.data(idx_eri_1_fi + 162);

    auto g_xzz_yyyyzz_1 = pbuffer.data(idx_eri_1_fi + 163);

    auto g_xzz_yyyzzz_1 = pbuffer.data(idx_eri_1_fi + 164);

    auto g_xzz_yyzzzz_1 = pbuffer.data(idx_eri_1_fi + 165);

    auto g_xzz_yzzzzz_1 = pbuffer.data(idx_eri_1_fi + 166);

    auto g_xzz_zzzzzz_1 = pbuffer.data(idx_eri_1_fi + 167);

    auto g_yyy_xxxxxx_1 = pbuffer.data(idx_eri_1_fi + 168);

    auto g_yyy_xxxxxy_1 = pbuffer.data(idx_eri_1_fi + 169);

    auto g_yyy_xxxxxz_1 = pbuffer.data(idx_eri_1_fi + 170);

    auto g_yyy_xxxxyy_1 = pbuffer.data(idx_eri_1_fi + 171);

    auto g_yyy_xxxxyz_1 = pbuffer.data(idx_eri_1_fi + 172);

    auto g_yyy_xxxxzz_1 = pbuffer.data(idx_eri_1_fi + 173);

    auto g_yyy_xxxyyy_1 = pbuffer.data(idx_eri_1_fi + 174);

    auto g_yyy_xxxyyz_1 = pbuffer.data(idx_eri_1_fi + 175);

    auto g_yyy_xxxyzz_1 = pbuffer.data(idx_eri_1_fi + 176);

    auto g_yyy_xxxzzz_1 = pbuffer.data(idx_eri_1_fi + 177);

    auto g_yyy_xxyyyy_1 = pbuffer.data(idx_eri_1_fi + 178);

    auto g_yyy_xxyyyz_1 = pbuffer.data(idx_eri_1_fi + 179);

    auto g_yyy_xxyyzz_1 = pbuffer.data(idx_eri_1_fi + 180);

    auto g_yyy_xxyzzz_1 = pbuffer.data(idx_eri_1_fi + 181);

    auto g_yyy_xxzzzz_1 = pbuffer.data(idx_eri_1_fi + 182);

    auto g_yyy_xyyyyy_1 = pbuffer.data(idx_eri_1_fi + 183);

    auto g_yyy_xyyyyz_1 = pbuffer.data(idx_eri_1_fi + 184);

    auto g_yyy_xyyyzz_1 = pbuffer.data(idx_eri_1_fi + 185);

    auto g_yyy_xyyzzz_1 = pbuffer.data(idx_eri_1_fi + 186);

    auto g_yyy_xyzzzz_1 = pbuffer.data(idx_eri_1_fi + 187);

    auto g_yyy_xzzzzz_1 = pbuffer.data(idx_eri_1_fi + 188);

    auto g_yyy_yyyyyy_1 = pbuffer.data(idx_eri_1_fi + 189);

    auto g_yyy_yyyyyz_1 = pbuffer.data(idx_eri_1_fi + 190);

    auto g_yyy_yyyyzz_1 = pbuffer.data(idx_eri_1_fi + 191);

    auto g_yyy_yyyzzz_1 = pbuffer.data(idx_eri_1_fi + 192);

    auto g_yyy_yyzzzz_1 = pbuffer.data(idx_eri_1_fi + 193);

    auto g_yyy_yzzzzz_1 = pbuffer.data(idx_eri_1_fi + 194);

    auto g_yyy_zzzzzz_1 = pbuffer.data(idx_eri_1_fi + 195);

    auto g_yyz_xxxxxy_1 = pbuffer.data(idx_eri_1_fi + 197);

    auto g_yyz_xxxxxz_1 = pbuffer.data(idx_eri_1_fi + 198);

    auto g_yyz_xxxxyy_1 = pbuffer.data(idx_eri_1_fi + 199);

    auto g_yyz_xxxxyz_1 = pbuffer.data(idx_eri_1_fi + 200);

    auto g_yyz_xxxxzz_1 = pbuffer.data(idx_eri_1_fi + 201);

    auto g_yyz_xxxyyy_1 = pbuffer.data(idx_eri_1_fi + 202);

    auto g_yyz_xxxyyz_1 = pbuffer.data(idx_eri_1_fi + 203);

    auto g_yyz_xxxyzz_1 = pbuffer.data(idx_eri_1_fi + 204);

    auto g_yyz_xxxzzz_1 = pbuffer.data(idx_eri_1_fi + 205);

    auto g_yyz_xxyyyy_1 = pbuffer.data(idx_eri_1_fi + 206);

    auto g_yyz_xxyyyz_1 = pbuffer.data(idx_eri_1_fi + 207);

    auto g_yyz_xxyyzz_1 = pbuffer.data(idx_eri_1_fi + 208);

    auto g_yyz_xxyzzz_1 = pbuffer.data(idx_eri_1_fi + 209);

    auto g_yyz_xxzzzz_1 = pbuffer.data(idx_eri_1_fi + 210);

    auto g_yyz_xyyyyy_1 = pbuffer.data(idx_eri_1_fi + 211);

    auto g_yyz_xyyyyz_1 = pbuffer.data(idx_eri_1_fi + 212);

    auto g_yyz_xyyyzz_1 = pbuffer.data(idx_eri_1_fi + 213);

    auto g_yyz_xyyzzz_1 = pbuffer.data(idx_eri_1_fi + 214);

    auto g_yyz_xyzzzz_1 = pbuffer.data(idx_eri_1_fi + 215);

    auto g_yyz_xzzzzz_1 = pbuffer.data(idx_eri_1_fi + 216);

    auto g_yyz_yyyyyy_1 = pbuffer.data(idx_eri_1_fi + 217);

    auto g_yyz_yyyyyz_1 = pbuffer.data(idx_eri_1_fi + 218);

    auto g_yyz_yyyyzz_1 = pbuffer.data(idx_eri_1_fi + 219);

    auto g_yyz_yyyzzz_1 = pbuffer.data(idx_eri_1_fi + 220);

    auto g_yyz_yyzzzz_1 = pbuffer.data(idx_eri_1_fi + 221);

    auto g_yyz_yzzzzz_1 = pbuffer.data(idx_eri_1_fi + 222);

    auto g_yyz_zzzzzz_1 = pbuffer.data(idx_eri_1_fi + 223);

    auto g_yzz_xxxxxx_1 = pbuffer.data(idx_eri_1_fi + 224);

    auto g_yzz_xxxxxy_1 = pbuffer.data(idx_eri_1_fi + 225);

    auto g_yzz_xxxxxz_1 = pbuffer.data(idx_eri_1_fi + 226);

    auto g_yzz_xxxxyy_1 = pbuffer.data(idx_eri_1_fi + 227);

    auto g_yzz_xxxxyz_1 = pbuffer.data(idx_eri_1_fi + 228);

    auto g_yzz_xxxxzz_1 = pbuffer.data(idx_eri_1_fi + 229);

    auto g_yzz_xxxyyy_1 = pbuffer.data(idx_eri_1_fi + 230);

    auto g_yzz_xxxyyz_1 = pbuffer.data(idx_eri_1_fi + 231);

    auto g_yzz_xxxyzz_1 = pbuffer.data(idx_eri_1_fi + 232);

    auto g_yzz_xxxzzz_1 = pbuffer.data(idx_eri_1_fi + 233);

    auto g_yzz_xxyyyy_1 = pbuffer.data(idx_eri_1_fi + 234);

    auto g_yzz_xxyyyz_1 = pbuffer.data(idx_eri_1_fi + 235);

    auto g_yzz_xxyyzz_1 = pbuffer.data(idx_eri_1_fi + 236);

    auto g_yzz_xxyzzz_1 = pbuffer.data(idx_eri_1_fi + 237);

    auto g_yzz_xxzzzz_1 = pbuffer.data(idx_eri_1_fi + 238);

    auto g_yzz_xyyyyy_1 = pbuffer.data(idx_eri_1_fi + 239);

    auto g_yzz_xyyyyz_1 = pbuffer.data(idx_eri_1_fi + 240);

    auto g_yzz_xyyyzz_1 = pbuffer.data(idx_eri_1_fi + 241);

    auto g_yzz_xyyzzz_1 = pbuffer.data(idx_eri_1_fi + 242);

    auto g_yzz_xyzzzz_1 = pbuffer.data(idx_eri_1_fi + 243);

    auto g_yzz_xzzzzz_1 = pbuffer.data(idx_eri_1_fi + 244);

    auto g_yzz_yyyyyy_1 = pbuffer.data(idx_eri_1_fi + 245);

    auto g_yzz_yyyyyz_1 = pbuffer.data(idx_eri_1_fi + 246);

    auto g_yzz_yyyyzz_1 = pbuffer.data(idx_eri_1_fi + 247);

    auto g_yzz_yyyzzz_1 = pbuffer.data(idx_eri_1_fi + 248);

    auto g_yzz_yyzzzz_1 = pbuffer.data(idx_eri_1_fi + 249);

    auto g_yzz_yzzzzz_1 = pbuffer.data(idx_eri_1_fi + 250);

    auto g_yzz_zzzzzz_1 = pbuffer.data(idx_eri_1_fi + 251);

    auto g_zzz_xxxxxx_1 = pbuffer.data(idx_eri_1_fi + 252);

    auto g_zzz_xxxxxy_1 = pbuffer.data(idx_eri_1_fi + 253);

    auto g_zzz_xxxxxz_1 = pbuffer.data(idx_eri_1_fi + 254);

    auto g_zzz_xxxxyy_1 = pbuffer.data(idx_eri_1_fi + 255);

    auto g_zzz_xxxxyz_1 = pbuffer.data(idx_eri_1_fi + 256);

    auto g_zzz_xxxxzz_1 = pbuffer.data(idx_eri_1_fi + 257);

    auto g_zzz_xxxyyy_1 = pbuffer.data(idx_eri_1_fi + 258);

    auto g_zzz_xxxyyz_1 = pbuffer.data(idx_eri_1_fi + 259);

    auto g_zzz_xxxyzz_1 = pbuffer.data(idx_eri_1_fi + 260);

    auto g_zzz_xxxzzz_1 = pbuffer.data(idx_eri_1_fi + 261);

    auto g_zzz_xxyyyy_1 = pbuffer.data(idx_eri_1_fi + 262);

    auto g_zzz_xxyyyz_1 = pbuffer.data(idx_eri_1_fi + 263);

    auto g_zzz_xxyyzz_1 = pbuffer.data(idx_eri_1_fi + 264);

    auto g_zzz_xxyzzz_1 = pbuffer.data(idx_eri_1_fi + 265);

    auto g_zzz_xxzzzz_1 = pbuffer.data(idx_eri_1_fi + 266);

    auto g_zzz_xyyyyy_1 = pbuffer.data(idx_eri_1_fi + 267);

    auto g_zzz_xyyyyz_1 = pbuffer.data(idx_eri_1_fi + 268);

    auto g_zzz_xyyyzz_1 = pbuffer.data(idx_eri_1_fi + 269);

    auto g_zzz_xyyzzz_1 = pbuffer.data(idx_eri_1_fi + 270);

    auto g_zzz_xyzzzz_1 = pbuffer.data(idx_eri_1_fi + 271);

    auto g_zzz_xzzzzz_1 = pbuffer.data(idx_eri_1_fi + 272);

    auto g_zzz_yyyyyy_1 = pbuffer.data(idx_eri_1_fi + 273);

    auto g_zzz_yyyyyz_1 = pbuffer.data(idx_eri_1_fi + 274);

    auto g_zzz_yyyyzz_1 = pbuffer.data(idx_eri_1_fi + 275);

    auto g_zzz_yyyzzz_1 = pbuffer.data(idx_eri_1_fi + 276);

    auto g_zzz_yyzzzz_1 = pbuffer.data(idx_eri_1_fi + 277);

    auto g_zzz_yzzzzz_1 = pbuffer.data(idx_eri_1_fi + 278);

    auto g_zzz_zzzzzz_1 = pbuffer.data(idx_eri_1_fi + 279);

    // Set up 0-28 components of targeted buffer : GI

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

    #pragma omp simd aligned(g_xx_xxxxxx_0, g_xx_xxxxxx_1, g_xx_xxxxxy_0, g_xx_xxxxxy_1, g_xx_xxxxxz_0, g_xx_xxxxxz_1, g_xx_xxxxyy_0, g_xx_xxxxyy_1, g_xx_xxxxyz_0, g_xx_xxxxyz_1, g_xx_xxxxzz_0, g_xx_xxxxzz_1, g_xx_xxxyyy_0, g_xx_xxxyyy_1, g_xx_xxxyyz_0, g_xx_xxxyyz_1, g_xx_xxxyzz_0, g_xx_xxxyzz_1, g_xx_xxxzzz_0, g_xx_xxxzzz_1, g_xx_xxyyyy_0, g_xx_xxyyyy_1, g_xx_xxyyyz_0, g_xx_xxyyyz_1, g_xx_xxyyzz_0, g_xx_xxyyzz_1, g_xx_xxyzzz_0, g_xx_xxyzzz_1, g_xx_xxzzzz_0, g_xx_xxzzzz_1, g_xx_xyyyyy_0, g_xx_xyyyyy_1, g_xx_xyyyyz_0, g_xx_xyyyyz_1, g_xx_xyyyzz_0, g_xx_xyyyzz_1, g_xx_xyyzzz_0, g_xx_xyyzzz_1, g_xx_xyzzzz_0, g_xx_xyzzzz_1, g_xx_xzzzzz_0, g_xx_xzzzzz_1, g_xx_yyyyyy_0, g_xx_yyyyyy_1, g_xx_yyyyyz_0, g_xx_yyyyyz_1, g_xx_yyyyzz_0, g_xx_yyyyzz_1, g_xx_yyyzzz_0, g_xx_yyyzzz_1, g_xx_yyzzzz_0, g_xx_yyzzzz_1, g_xx_yzzzzz_0, g_xx_yzzzzz_1, g_xx_zzzzzz_0, g_xx_zzzzzz_1, g_xxx_xxxxx_1, g_xxx_xxxxxx_1, g_xxx_xxxxxy_1, g_xxx_xxxxxz_1, g_xxx_xxxxy_1, g_xxx_xxxxyy_1, g_xxx_xxxxyz_1, g_xxx_xxxxz_1, g_xxx_xxxxzz_1, g_xxx_xxxyy_1, g_xxx_xxxyyy_1, g_xxx_xxxyyz_1, g_xxx_xxxyz_1, g_xxx_xxxyzz_1, g_xxx_xxxzz_1, g_xxx_xxxzzz_1, g_xxx_xxyyy_1, g_xxx_xxyyyy_1, g_xxx_xxyyyz_1, g_xxx_xxyyz_1, g_xxx_xxyyzz_1, g_xxx_xxyzz_1, g_xxx_xxyzzz_1, g_xxx_xxzzz_1, g_xxx_xxzzzz_1, g_xxx_xyyyy_1, g_xxx_xyyyyy_1, g_xxx_xyyyyz_1, g_xxx_xyyyz_1, g_xxx_xyyyzz_1, g_xxx_xyyzz_1, g_xxx_xyyzzz_1, g_xxx_xyzzz_1, g_xxx_xyzzzz_1, g_xxx_xzzzz_1, g_xxx_xzzzzz_1, g_xxx_yyyyy_1, g_xxx_yyyyyy_1, g_xxx_yyyyyz_1, g_xxx_yyyyz_1, g_xxx_yyyyzz_1, g_xxx_yyyzz_1, g_xxx_yyyzzz_1, g_xxx_yyzzz_1, g_xxx_yyzzzz_1, g_xxx_yzzzz_1, g_xxx_yzzzzz_1, g_xxx_zzzzz_1, g_xxx_zzzzzz_1, g_xxxx_xxxxxx_0, g_xxxx_xxxxxy_0, g_xxxx_xxxxxz_0, g_xxxx_xxxxyy_0, g_xxxx_xxxxyz_0, g_xxxx_xxxxzz_0, g_xxxx_xxxyyy_0, g_xxxx_xxxyyz_0, g_xxxx_xxxyzz_0, g_xxxx_xxxzzz_0, g_xxxx_xxyyyy_0, g_xxxx_xxyyyz_0, g_xxxx_xxyyzz_0, g_xxxx_xxyzzz_0, g_xxxx_xxzzzz_0, g_xxxx_xyyyyy_0, g_xxxx_xyyyyz_0, g_xxxx_xyyyzz_0, g_xxxx_xyyzzz_0, g_xxxx_xyzzzz_0, g_xxxx_xzzzzz_0, g_xxxx_yyyyyy_0, g_xxxx_yyyyyz_0, g_xxxx_yyyyzz_0, g_xxxx_yyyzzz_0, g_xxxx_yyzzzz_0, g_xxxx_yzzzzz_0, g_xxxx_zzzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxx_xxxxxx_0[i] = 3.0 * g_xx_xxxxxx_0[i] * fbe_0 - 3.0 * g_xx_xxxxxx_1[i] * fz_be_0 + 6.0 * g_xxx_xxxxx_1[i] * fe_0 + g_xxx_xxxxxx_1[i] * pa_x[i];

        g_xxxx_xxxxxy_0[i] = 3.0 * g_xx_xxxxxy_0[i] * fbe_0 - 3.0 * g_xx_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxx_xxxxy_1[i] * fe_0 + g_xxx_xxxxxy_1[i] * pa_x[i];

        g_xxxx_xxxxxz_0[i] = 3.0 * g_xx_xxxxxz_0[i] * fbe_0 - 3.0 * g_xx_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxx_xxxxz_1[i] * fe_0 + g_xxx_xxxxxz_1[i] * pa_x[i];

        g_xxxx_xxxxyy_0[i] = 3.0 * g_xx_xxxxyy_0[i] * fbe_0 - 3.0 * g_xx_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxx_xxxyy_1[i] * fe_0 + g_xxx_xxxxyy_1[i] * pa_x[i];

        g_xxxx_xxxxyz_0[i] = 3.0 * g_xx_xxxxyz_0[i] * fbe_0 - 3.0 * g_xx_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxx_xxxyz_1[i] * fe_0 + g_xxx_xxxxyz_1[i] * pa_x[i];

        g_xxxx_xxxxzz_0[i] = 3.0 * g_xx_xxxxzz_0[i] * fbe_0 - 3.0 * g_xx_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxx_xxxzz_1[i] * fe_0 + g_xxx_xxxxzz_1[i] * pa_x[i];

        g_xxxx_xxxyyy_0[i] = 3.0 * g_xx_xxxyyy_0[i] * fbe_0 - 3.0 * g_xx_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxx_xxyyy_1[i] * fe_0 + g_xxx_xxxyyy_1[i] * pa_x[i];

        g_xxxx_xxxyyz_0[i] = 3.0 * g_xx_xxxyyz_0[i] * fbe_0 - 3.0 * g_xx_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxx_xxyyz_1[i] * fe_0 + g_xxx_xxxyyz_1[i] * pa_x[i];

        g_xxxx_xxxyzz_0[i] = 3.0 * g_xx_xxxyzz_0[i] * fbe_0 - 3.0 * g_xx_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxx_xxyzz_1[i] * fe_0 + g_xxx_xxxyzz_1[i] * pa_x[i];

        g_xxxx_xxxzzz_0[i] = 3.0 * g_xx_xxxzzz_0[i] * fbe_0 - 3.0 * g_xx_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxx_xxzzz_1[i] * fe_0 + g_xxx_xxxzzz_1[i] * pa_x[i];

        g_xxxx_xxyyyy_0[i] = 3.0 * g_xx_xxyyyy_0[i] * fbe_0 - 3.0 * g_xx_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxx_xyyyy_1[i] * fe_0 + g_xxx_xxyyyy_1[i] * pa_x[i];

        g_xxxx_xxyyyz_0[i] = 3.0 * g_xx_xxyyyz_0[i] * fbe_0 - 3.0 * g_xx_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxx_xyyyz_1[i] * fe_0 + g_xxx_xxyyyz_1[i] * pa_x[i];

        g_xxxx_xxyyzz_0[i] = 3.0 * g_xx_xxyyzz_0[i] * fbe_0 - 3.0 * g_xx_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxx_xyyzz_1[i] * fe_0 + g_xxx_xxyyzz_1[i] * pa_x[i];

        g_xxxx_xxyzzz_0[i] = 3.0 * g_xx_xxyzzz_0[i] * fbe_0 - 3.0 * g_xx_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxx_xyzzz_1[i] * fe_0 + g_xxx_xxyzzz_1[i] * pa_x[i];

        g_xxxx_xxzzzz_0[i] = 3.0 * g_xx_xxzzzz_0[i] * fbe_0 - 3.0 * g_xx_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxx_xzzzz_1[i] * fe_0 + g_xxx_xxzzzz_1[i] * pa_x[i];

        g_xxxx_xyyyyy_0[i] = 3.0 * g_xx_xyyyyy_0[i] * fbe_0 - 3.0 * g_xx_xyyyyy_1[i] * fz_be_0 + g_xxx_yyyyy_1[i] * fe_0 + g_xxx_xyyyyy_1[i] * pa_x[i];

        g_xxxx_xyyyyz_0[i] = 3.0 * g_xx_xyyyyz_0[i] * fbe_0 - 3.0 * g_xx_xyyyyz_1[i] * fz_be_0 + g_xxx_yyyyz_1[i] * fe_0 + g_xxx_xyyyyz_1[i] * pa_x[i];

        g_xxxx_xyyyzz_0[i] = 3.0 * g_xx_xyyyzz_0[i] * fbe_0 - 3.0 * g_xx_xyyyzz_1[i] * fz_be_0 + g_xxx_yyyzz_1[i] * fe_0 + g_xxx_xyyyzz_1[i] * pa_x[i];

        g_xxxx_xyyzzz_0[i] = 3.0 * g_xx_xyyzzz_0[i] * fbe_0 - 3.0 * g_xx_xyyzzz_1[i] * fz_be_0 + g_xxx_yyzzz_1[i] * fe_0 + g_xxx_xyyzzz_1[i] * pa_x[i];

        g_xxxx_xyzzzz_0[i] = 3.0 * g_xx_xyzzzz_0[i] * fbe_0 - 3.0 * g_xx_xyzzzz_1[i] * fz_be_0 + g_xxx_yzzzz_1[i] * fe_0 + g_xxx_xyzzzz_1[i] * pa_x[i];

        g_xxxx_xzzzzz_0[i] = 3.0 * g_xx_xzzzzz_0[i] * fbe_0 - 3.0 * g_xx_xzzzzz_1[i] * fz_be_0 + g_xxx_zzzzz_1[i] * fe_0 + g_xxx_xzzzzz_1[i] * pa_x[i];

        g_xxxx_yyyyyy_0[i] = 3.0 * g_xx_yyyyyy_0[i] * fbe_0 - 3.0 * g_xx_yyyyyy_1[i] * fz_be_0 + g_xxx_yyyyyy_1[i] * pa_x[i];

        g_xxxx_yyyyyz_0[i] = 3.0 * g_xx_yyyyyz_0[i] * fbe_0 - 3.0 * g_xx_yyyyyz_1[i] * fz_be_0 + g_xxx_yyyyyz_1[i] * pa_x[i];

        g_xxxx_yyyyzz_0[i] = 3.0 * g_xx_yyyyzz_0[i] * fbe_0 - 3.0 * g_xx_yyyyzz_1[i] * fz_be_0 + g_xxx_yyyyzz_1[i] * pa_x[i];

        g_xxxx_yyyzzz_0[i] = 3.0 * g_xx_yyyzzz_0[i] * fbe_0 - 3.0 * g_xx_yyyzzz_1[i] * fz_be_0 + g_xxx_yyyzzz_1[i] * pa_x[i];

        g_xxxx_yyzzzz_0[i] = 3.0 * g_xx_yyzzzz_0[i] * fbe_0 - 3.0 * g_xx_yyzzzz_1[i] * fz_be_0 + g_xxx_yyzzzz_1[i] * pa_x[i];

        g_xxxx_yzzzzz_0[i] = 3.0 * g_xx_yzzzzz_0[i] * fbe_0 - 3.0 * g_xx_yzzzzz_1[i] * fz_be_0 + g_xxx_yzzzzz_1[i] * pa_x[i];

        g_xxxx_zzzzzz_0[i] = 3.0 * g_xx_zzzzzz_0[i] * fbe_0 - 3.0 * g_xx_zzzzzz_1[i] * fz_be_0 + g_xxx_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : GI

    auto g_xxxy_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 28);

    auto g_xxxy_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 29);

    auto g_xxxy_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 30);

    auto g_xxxy_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 31);

    auto g_xxxy_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 32);

    auto g_xxxy_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 33);

    auto g_xxxy_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 34);

    auto g_xxxy_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 35);

    auto g_xxxy_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 36);

    auto g_xxxy_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 37);

    auto g_xxxy_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 38);

    auto g_xxxy_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 39);

    auto g_xxxy_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 40);

    auto g_xxxy_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 41);

    auto g_xxxy_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 42);

    auto g_xxxy_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 43);

    auto g_xxxy_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 44);

    auto g_xxxy_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 45);

    auto g_xxxy_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 46);

    auto g_xxxy_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 47);

    auto g_xxxy_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 48);

    auto g_xxxy_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 49);

    auto g_xxxy_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 50);

    auto g_xxxy_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 51);

    auto g_xxxy_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 52);

    auto g_xxxy_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 53);

    auto g_xxxy_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 54);

    auto g_xxxy_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 55);

    #pragma omp simd aligned(g_xxx_xxxxx_1, g_xxx_xxxxxx_1, g_xxx_xxxxxy_1, g_xxx_xxxxxz_1, g_xxx_xxxxy_1, g_xxx_xxxxyy_1, g_xxx_xxxxyz_1, g_xxx_xxxxz_1, g_xxx_xxxxzz_1, g_xxx_xxxyy_1, g_xxx_xxxyyy_1, g_xxx_xxxyyz_1, g_xxx_xxxyz_1, g_xxx_xxxyzz_1, g_xxx_xxxzz_1, g_xxx_xxxzzz_1, g_xxx_xxyyy_1, g_xxx_xxyyyy_1, g_xxx_xxyyyz_1, g_xxx_xxyyz_1, g_xxx_xxyyzz_1, g_xxx_xxyzz_1, g_xxx_xxyzzz_1, g_xxx_xxzzz_1, g_xxx_xxzzzz_1, g_xxx_xyyyy_1, g_xxx_xyyyyy_1, g_xxx_xyyyyz_1, g_xxx_xyyyz_1, g_xxx_xyyyzz_1, g_xxx_xyyzz_1, g_xxx_xyyzzz_1, g_xxx_xyzzz_1, g_xxx_xyzzzz_1, g_xxx_xzzzz_1, g_xxx_xzzzzz_1, g_xxx_yyyyy_1, g_xxx_yyyyyy_1, g_xxx_yyyyyz_1, g_xxx_yyyyz_1, g_xxx_yyyyzz_1, g_xxx_yyyzz_1, g_xxx_yyyzzz_1, g_xxx_yyzzz_1, g_xxx_yyzzzz_1, g_xxx_yzzzz_1, g_xxx_yzzzzz_1, g_xxx_zzzzz_1, g_xxx_zzzzzz_1, g_xxxy_xxxxxx_0, g_xxxy_xxxxxy_0, g_xxxy_xxxxxz_0, g_xxxy_xxxxyy_0, g_xxxy_xxxxyz_0, g_xxxy_xxxxzz_0, g_xxxy_xxxyyy_0, g_xxxy_xxxyyz_0, g_xxxy_xxxyzz_0, g_xxxy_xxxzzz_0, g_xxxy_xxyyyy_0, g_xxxy_xxyyyz_0, g_xxxy_xxyyzz_0, g_xxxy_xxyzzz_0, g_xxxy_xxzzzz_0, g_xxxy_xyyyyy_0, g_xxxy_xyyyyz_0, g_xxxy_xyyyzz_0, g_xxxy_xyyzzz_0, g_xxxy_xyzzzz_0, g_xxxy_xzzzzz_0, g_xxxy_yyyyyy_0, g_xxxy_yyyyyz_0, g_xxxy_yyyyzz_0, g_xxxy_yyyzzz_0, g_xxxy_yyzzzz_0, g_xxxy_yzzzzz_0, g_xxxy_zzzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxy_xxxxxx_0[i] = g_xxx_xxxxxx_1[i] * pa_y[i];

        g_xxxy_xxxxxy_0[i] = g_xxx_xxxxx_1[i] * fe_0 + g_xxx_xxxxxy_1[i] * pa_y[i];

        g_xxxy_xxxxxz_0[i] = g_xxx_xxxxxz_1[i] * pa_y[i];

        g_xxxy_xxxxyy_0[i] = 2.0 * g_xxx_xxxxy_1[i] * fe_0 + g_xxx_xxxxyy_1[i] * pa_y[i];

        g_xxxy_xxxxyz_0[i] = g_xxx_xxxxz_1[i] * fe_0 + g_xxx_xxxxyz_1[i] * pa_y[i];

        g_xxxy_xxxxzz_0[i] = g_xxx_xxxxzz_1[i] * pa_y[i];

        g_xxxy_xxxyyy_0[i] = 3.0 * g_xxx_xxxyy_1[i] * fe_0 + g_xxx_xxxyyy_1[i] * pa_y[i];

        g_xxxy_xxxyyz_0[i] = 2.0 * g_xxx_xxxyz_1[i] * fe_0 + g_xxx_xxxyyz_1[i] * pa_y[i];

        g_xxxy_xxxyzz_0[i] = g_xxx_xxxzz_1[i] * fe_0 + g_xxx_xxxyzz_1[i] * pa_y[i];

        g_xxxy_xxxzzz_0[i] = g_xxx_xxxzzz_1[i] * pa_y[i];

        g_xxxy_xxyyyy_0[i] = 4.0 * g_xxx_xxyyy_1[i] * fe_0 + g_xxx_xxyyyy_1[i] * pa_y[i];

        g_xxxy_xxyyyz_0[i] = 3.0 * g_xxx_xxyyz_1[i] * fe_0 + g_xxx_xxyyyz_1[i] * pa_y[i];

        g_xxxy_xxyyzz_0[i] = 2.0 * g_xxx_xxyzz_1[i] * fe_0 + g_xxx_xxyyzz_1[i] * pa_y[i];

        g_xxxy_xxyzzz_0[i] = g_xxx_xxzzz_1[i] * fe_0 + g_xxx_xxyzzz_1[i] * pa_y[i];

        g_xxxy_xxzzzz_0[i] = g_xxx_xxzzzz_1[i] * pa_y[i];

        g_xxxy_xyyyyy_0[i] = 5.0 * g_xxx_xyyyy_1[i] * fe_0 + g_xxx_xyyyyy_1[i] * pa_y[i];

        g_xxxy_xyyyyz_0[i] = 4.0 * g_xxx_xyyyz_1[i] * fe_0 + g_xxx_xyyyyz_1[i] * pa_y[i];

        g_xxxy_xyyyzz_0[i] = 3.0 * g_xxx_xyyzz_1[i] * fe_0 + g_xxx_xyyyzz_1[i] * pa_y[i];

        g_xxxy_xyyzzz_0[i] = 2.0 * g_xxx_xyzzz_1[i] * fe_0 + g_xxx_xyyzzz_1[i] * pa_y[i];

        g_xxxy_xyzzzz_0[i] = g_xxx_xzzzz_1[i] * fe_0 + g_xxx_xyzzzz_1[i] * pa_y[i];

        g_xxxy_xzzzzz_0[i] = g_xxx_xzzzzz_1[i] * pa_y[i];

        g_xxxy_yyyyyy_0[i] = 6.0 * g_xxx_yyyyy_1[i] * fe_0 + g_xxx_yyyyyy_1[i] * pa_y[i];

        g_xxxy_yyyyyz_0[i] = 5.0 * g_xxx_yyyyz_1[i] * fe_0 + g_xxx_yyyyyz_1[i] * pa_y[i];

        g_xxxy_yyyyzz_0[i] = 4.0 * g_xxx_yyyzz_1[i] * fe_0 + g_xxx_yyyyzz_1[i] * pa_y[i];

        g_xxxy_yyyzzz_0[i] = 3.0 * g_xxx_yyzzz_1[i] * fe_0 + g_xxx_yyyzzz_1[i] * pa_y[i];

        g_xxxy_yyzzzz_0[i] = 2.0 * g_xxx_yzzzz_1[i] * fe_0 + g_xxx_yyzzzz_1[i] * pa_y[i];

        g_xxxy_yzzzzz_0[i] = g_xxx_zzzzz_1[i] * fe_0 + g_xxx_yzzzzz_1[i] * pa_y[i];

        g_xxxy_zzzzzz_0[i] = g_xxx_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : GI

    auto g_xxxz_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 56);

    auto g_xxxz_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 57);

    auto g_xxxz_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 58);

    auto g_xxxz_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 59);

    auto g_xxxz_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 60);

    auto g_xxxz_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 61);

    auto g_xxxz_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 62);

    auto g_xxxz_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 63);

    auto g_xxxz_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 64);

    auto g_xxxz_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 65);

    auto g_xxxz_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 66);

    auto g_xxxz_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 67);

    auto g_xxxz_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 68);

    auto g_xxxz_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 69);

    auto g_xxxz_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 70);

    auto g_xxxz_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 71);

    auto g_xxxz_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 72);

    auto g_xxxz_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 73);

    auto g_xxxz_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 74);

    auto g_xxxz_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 75);

    auto g_xxxz_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 76);

    auto g_xxxz_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 77);

    auto g_xxxz_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 78);

    auto g_xxxz_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 79);

    auto g_xxxz_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 80);

    auto g_xxxz_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 81);

    auto g_xxxz_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 82);

    auto g_xxxz_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 83);

    #pragma omp simd aligned(g_xxx_xxxxx_1, g_xxx_xxxxxx_1, g_xxx_xxxxxy_1, g_xxx_xxxxxz_1, g_xxx_xxxxy_1, g_xxx_xxxxyy_1, g_xxx_xxxxyz_1, g_xxx_xxxxz_1, g_xxx_xxxxzz_1, g_xxx_xxxyy_1, g_xxx_xxxyyy_1, g_xxx_xxxyyz_1, g_xxx_xxxyz_1, g_xxx_xxxyzz_1, g_xxx_xxxzz_1, g_xxx_xxxzzz_1, g_xxx_xxyyy_1, g_xxx_xxyyyy_1, g_xxx_xxyyyz_1, g_xxx_xxyyz_1, g_xxx_xxyyzz_1, g_xxx_xxyzz_1, g_xxx_xxyzzz_1, g_xxx_xxzzz_1, g_xxx_xxzzzz_1, g_xxx_xyyyy_1, g_xxx_xyyyyy_1, g_xxx_xyyyyz_1, g_xxx_xyyyz_1, g_xxx_xyyyzz_1, g_xxx_xyyzz_1, g_xxx_xyyzzz_1, g_xxx_xyzzz_1, g_xxx_xyzzzz_1, g_xxx_xzzzz_1, g_xxx_xzzzzz_1, g_xxx_yyyyy_1, g_xxx_yyyyyy_1, g_xxx_yyyyyz_1, g_xxx_yyyyz_1, g_xxx_yyyyzz_1, g_xxx_yyyzz_1, g_xxx_yyyzzz_1, g_xxx_yyzzz_1, g_xxx_yyzzzz_1, g_xxx_yzzzz_1, g_xxx_yzzzzz_1, g_xxx_zzzzz_1, g_xxx_zzzzzz_1, g_xxxz_xxxxxx_0, g_xxxz_xxxxxy_0, g_xxxz_xxxxxz_0, g_xxxz_xxxxyy_0, g_xxxz_xxxxyz_0, g_xxxz_xxxxzz_0, g_xxxz_xxxyyy_0, g_xxxz_xxxyyz_0, g_xxxz_xxxyzz_0, g_xxxz_xxxzzz_0, g_xxxz_xxyyyy_0, g_xxxz_xxyyyz_0, g_xxxz_xxyyzz_0, g_xxxz_xxyzzz_0, g_xxxz_xxzzzz_0, g_xxxz_xyyyyy_0, g_xxxz_xyyyyz_0, g_xxxz_xyyyzz_0, g_xxxz_xyyzzz_0, g_xxxz_xyzzzz_0, g_xxxz_xzzzzz_0, g_xxxz_yyyyyy_0, g_xxxz_yyyyyz_0, g_xxxz_yyyyzz_0, g_xxxz_yyyzzz_0, g_xxxz_yyzzzz_0, g_xxxz_yzzzzz_0, g_xxxz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxz_xxxxxx_0[i] = g_xxx_xxxxxx_1[i] * pa_z[i];

        g_xxxz_xxxxxy_0[i] = g_xxx_xxxxxy_1[i] * pa_z[i];

        g_xxxz_xxxxxz_0[i] = g_xxx_xxxxx_1[i] * fe_0 + g_xxx_xxxxxz_1[i] * pa_z[i];

        g_xxxz_xxxxyy_0[i] = g_xxx_xxxxyy_1[i] * pa_z[i];

        g_xxxz_xxxxyz_0[i] = g_xxx_xxxxy_1[i] * fe_0 + g_xxx_xxxxyz_1[i] * pa_z[i];

        g_xxxz_xxxxzz_0[i] = 2.0 * g_xxx_xxxxz_1[i] * fe_0 + g_xxx_xxxxzz_1[i] * pa_z[i];

        g_xxxz_xxxyyy_0[i] = g_xxx_xxxyyy_1[i] * pa_z[i];

        g_xxxz_xxxyyz_0[i] = g_xxx_xxxyy_1[i] * fe_0 + g_xxx_xxxyyz_1[i] * pa_z[i];

        g_xxxz_xxxyzz_0[i] = 2.0 * g_xxx_xxxyz_1[i] * fe_0 + g_xxx_xxxyzz_1[i] * pa_z[i];

        g_xxxz_xxxzzz_0[i] = 3.0 * g_xxx_xxxzz_1[i] * fe_0 + g_xxx_xxxzzz_1[i] * pa_z[i];

        g_xxxz_xxyyyy_0[i] = g_xxx_xxyyyy_1[i] * pa_z[i];

        g_xxxz_xxyyyz_0[i] = g_xxx_xxyyy_1[i] * fe_0 + g_xxx_xxyyyz_1[i] * pa_z[i];

        g_xxxz_xxyyzz_0[i] = 2.0 * g_xxx_xxyyz_1[i] * fe_0 + g_xxx_xxyyzz_1[i] * pa_z[i];

        g_xxxz_xxyzzz_0[i] = 3.0 * g_xxx_xxyzz_1[i] * fe_0 + g_xxx_xxyzzz_1[i] * pa_z[i];

        g_xxxz_xxzzzz_0[i] = 4.0 * g_xxx_xxzzz_1[i] * fe_0 + g_xxx_xxzzzz_1[i] * pa_z[i];

        g_xxxz_xyyyyy_0[i] = g_xxx_xyyyyy_1[i] * pa_z[i];

        g_xxxz_xyyyyz_0[i] = g_xxx_xyyyy_1[i] * fe_0 + g_xxx_xyyyyz_1[i] * pa_z[i];

        g_xxxz_xyyyzz_0[i] = 2.0 * g_xxx_xyyyz_1[i] * fe_0 + g_xxx_xyyyzz_1[i] * pa_z[i];

        g_xxxz_xyyzzz_0[i] = 3.0 * g_xxx_xyyzz_1[i] * fe_0 + g_xxx_xyyzzz_1[i] * pa_z[i];

        g_xxxz_xyzzzz_0[i] = 4.0 * g_xxx_xyzzz_1[i] * fe_0 + g_xxx_xyzzzz_1[i] * pa_z[i];

        g_xxxz_xzzzzz_0[i] = 5.0 * g_xxx_xzzzz_1[i] * fe_0 + g_xxx_xzzzzz_1[i] * pa_z[i];

        g_xxxz_yyyyyy_0[i] = g_xxx_yyyyyy_1[i] * pa_z[i];

        g_xxxz_yyyyyz_0[i] = g_xxx_yyyyy_1[i] * fe_0 + g_xxx_yyyyyz_1[i] * pa_z[i];

        g_xxxz_yyyyzz_0[i] = 2.0 * g_xxx_yyyyz_1[i] * fe_0 + g_xxx_yyyyzz_1[i] * pa_z[i];

        g_xxxz_yyyzzz_0[i] = 3.0 * g_xxx_yyyzz_1[i] * fe_0 + g_xxx_yyyzzz_1[i] * pa_z[i];

        g_xxxz_yyzzzz_0[i] = 4.0 * g_xxx_yyzzz_1[i] * fe_0 + g_xxx_yyzzzz_1[i] * pa_z[i];

        g_xxxz_yzzzzz_0[i] = 5.0 * g_xxx_yzzzz_1[i] * fe_0 + g_xxx_yzzzzz_1[i] * pa_z[i];

        g_xxxz_zzzzzz_0[i] = 6.0 * g_xxx_zzzzz_1[i] * fe_0 + g_xxx_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 84-112 components of targeted buffer : GI

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

    #pragma omp simd aligned(g_xx_xxxxxx_0, g_xx_xxxxxx_1, g_xx_xxxxxz_0, g_xx_xxxxxz_1, g_xx_xxxxzz_0, g_xx_xxxxzz_1, g_xx_xxxzzz_0, g_xx_xxxzzz_1, g_xx_xxzzzz_0, g_xx_xxzzzz_1, g_xx_xzzzzz_0, g_xx_xzzzzz_1, g_xxy_xxxxxx_1, g_xxy_xxxxxz_1, g_xxy_xxxxzz_1, g_xxy_xxxzzz_1, g_xxy_xxzzzz_1, g_xxy_xzzzzz_1, g_xxyy_xxxxxx_0, g_xxyy_xxxxxy_0, g_xxyy_xxxxxz_0, g_xxyy_xxxxyy_0, g_xxyy_xxxxyz_0, g_xxyy_xxxxzz_0, g_xxyy_xxxyyy_0, g_xxyy_xxxyyz_0, g_xxyy_xxxyzz_0, g_xxyy_xxxzzz_0, g_xxyy_xxyyyy_0, g_xxyy_xxyyyz_0, g_xxyy_xxyyzz_0, g_xxyy_xxyzzz_0, g_xxyy_xxzzzz_0, g_xxyy_xyyyyy_0, g_xxyy_xyyyyz_0, g_xxyy_xyyyzz_0, g_xxyy_xyyzzz_0, g_xxyy_xyzzzz_0, g_xxyy_xzzzzz_0, g_xxyy_yyyyyy_0, g_xxyy_yyyyyz_0, g_xxyy_yyyyzz_0, g_xxyy_yyyzzz_0, g_xxyy_yyzzzz_0, g_xxyy_yzzzzz_0, g_xxyy_zzzzzz_0, g_xyy_xxxxxy_1, g_xyy_xxxxy_1, g_xyy_xxxxyy_1, g_xyy_xxxxyz_1, g_xyy_xxxyy_1, g_xyy_xxxyyy_1, g_xyy_xxxyyz_1, g_xyy_xxxyz_1, g_xyy_xxxyzz_1, g_xyy_xxyyy_1, g_xyy_xxyyyy_1, g_xyy_xxyyyz_1, g_xyy_xxyyz_1, g_xyy_xxyyzz_1, g_xyy_xxyzz_1, g_xyy_xxyzzz_1, g_xyy_xyyyy_1, g_xyy_xyyyyy_1, g_xyy_xyyyyz_1, g_xyy_xyyyz_1, g_xyy_xyyyzz_1, g_xyy_xyyzz_1, g_xyy_xyyzzz_1, g_xyy_xyzzz_1, g_xyy_xyzzzz_1, g_xyy_yyyyy_1, g_xyy_yyyyyy_1, g_xyy_yyyyyz_1, g_xyy_yyyyz_1, g_xyy_yyyyzz_1, g_xyy_yyyzz_1, g_xyy_yyyzzz_1, g_xyy_yyzzz_1, g_xyy_yyzzzz_1, g_xyy_yzzzz_1, g_xyy_yzzzzz_1, g_xyy_zzzzzz_1, g_yy_xxxxxy_0, g_yy_xxxxxy_1, g_yy_xxxxyy_0, g_yy_xxxxyy_1, g_yy_xxxxyz_0, g_yy_xxxxyz_1, g_yy_xxxyyy_0, g_yy_xxxyyy_1, g_yy_xxxyyz_0, g_yy_xxxyyz_1, g_yy_xxxyzz_0, g_yy_xxxyzz_1, g_yy_xxyyyy_0, g_yy_xxyyyy_1, g_yy_xxyyyz_0, g_yy_xxyyyz_1, g_yy_xxyyzz_0, g_yy_xxyyzz_1, g_yy_xxyzzz_0, g_yy_xxyzzz_1, g_yy_xyyyyy_0, g_yy_xyyyyy_1, g_yy_xyyyyz_0, g_yy_xyyyyz_1, g_yy_xyyyzz_0, g_yy_xyyyzz_1, g_yy_xyyzzz_0, g_yy_xyyzzz_1, g_yy_xyzzzz_0, g_yy_xyzzzz_1, g_yy_yyyyyy_0, g_yy_yyyyyy_1, g_yy_yyyyyz_0, g_yy_yyyyyz_1, g_yy_yyyyzz_0, g_yy_yyyyzz_1, g_yy_yyyzzz_0, g_yy_yyyzzz_1, g_yy_yyzzzz_0, g_yy_yyzzzz_1, g_yy_yzzzzz_0, g_yy_yzzzzz_1, g_yy_zzzzzz_0, g_yy_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyy_xxxxxx_0[i] = g_xx_xxxxxx_0[i] * fbe_0 - g_xx_xxxxxx_1[i] * fz_be_0 + g_xxy_xxxxxx_1[i] * pa_y[i];

        g_xxyy_xxxxxy_0[i] = g_yy_xxxxxy_0[i] * fbe_0 - g_yy_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xyy_xxxxy_1[i] * fe_0 + g_xyy_xxxxxy_1[i] * pa_x[i];

        g_xxyy_xxxxxz_0[i] = g_xx_xxxxxz_0[i] * fbe_0 - g_xx_xxxxxz_1[i] * fz_be_0 + g_xxy_xxxxxz_1[i] * pa_y[i];

        g_xxyy_xxxxyy_0[i] = g_yy_xxxxyy_0[i] * fbe_0 - g_yy_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xyy_xxxyy_1[i] * fe_0 + g_xyy_xxxxyy_1[i] * pa_x[i];

        g_xxyy_xxxxyz_0[i] = g_yy_xxxxyz_0[i] * fbe_0 - g_yy_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyy_xxxyz_1[i] * fe_0 + g_xyy_xxxxyz_1[i] * pa_x[i];

        g_xxyy_xxxxzz_0[i] = g_xx_xxxxzz_0[i] * fbe_0 - g_xx_xxxxzz_1[i] * fz_be_0 + g_xxy_xxxxzz_1[i] * pa_y[i];

        g_xxyy_xxxyyy_0[i] = g_yy_xxxyyy_0[i] * fbe_0 - g_yy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xyy_xxyyy_1[i] * fe_0 + g_xyy_xxxyyy_1[i] * pa_x[i];

        g_xxyy_xxxyyz_0[i] = g_yy_xxxyyz_0[i] * fbe_0 - g_yy_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyy_xxyyz_1[i] * fe_0 + g_xyy_xxxyyz_1[i] * pa_x[i];

        g_xxyy_xxxyzz_0[i] = g_yy_xxxyzz_0[i] * fbe_0 - g_yy_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyy_xxyzz_1[i] * fe_0 + g_xyy_xxxyzz_1[i] * pa_x[i];

        g_xxyy_xxxzzz_0[i] = g_xx_xxxzzz_0[i] * fbe_0 - g_xx_xxxzzz_1[i] * fz_be_0 + g_xxy_xxxzzz_1[i] * pa_y[i];

        g_xxyy_xxyyyy_0[i] = g_yy_xxyyyy_0[i] * fbe_0 - g_yy_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xyy_xyyyy_1[i] * fe_0 + g_xyy_xxyyyy_1[i] * pa_x[i];

        g_xxyy_xxyyyz_0[i] = g_yy_xxyyyz_0[i] * fbe_0 - g_yy_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyy_xyyyz_1[i] * fe_0 + g_xyy_xxyyyz_1[i] * pa_x[i];

        g_xxyy_xxyyzz_0[i] = g_yy_xxyyzz_0[i] * fbe_0 - g_yy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyy_xyyzz_1[i] * fe_0 + g_xyy_xxyyzz_1[i] * pa_x[i];

        g_xxyy_xxyzzz_0[i] = g_yy_xxyzzz_0[i] * fbe_0 - g_yy_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyy_xyzzz_1[i] * fe_0 + g_xyy_xxyzzz_1[i] * pa_x[i];

        g_xxyy_xxzzzz_0[i] = g_xx_xxzzzz_0[i] * fbe_0 - g_xx_xxzzzz_1[i] * fz_be_0 + g_xxy_xxzzzz_1[i] * pa_y[i];

        g_xxyy_xyyyyy_0[i] = g_yy_xyyyyy_0[i] * fbe_0 - g_yy_xyyyyy_1[i] * fz_be_0 + g_xyy_yyyyy_1[i] * fe_0 + g_xyy_xyyyyy_1[i] * pa_x[i];

        g_xxyy_xyyyyz_0[i] = g_yy_xyyyyz_0[i] * fbe_0 - g_yy_xyyyyz_1[i] * fz_be_0 + g_xyy_yyyyz_1[i] * fe_0 + g_xyy_xyyyyz_1[i] * pa_x[i];

        g_xxyy_xyyyzz_0[i] = g_yy_xyyyzz_0[i] * fbe_0 - g_yy_xyyyzz_1[i] * fz_be_0 + g_xyy_yyyzz_1[i] * fe_0 + g_xyy_xyyyzz_1[i] * pa_x[i];

        g_xxyy_xyyzzz_0[i] = g_yy_xyyzzz_0[i] * fbe_0 - g_yy_xyyzzz_1[i] * fz_be_0 + g_xyy_yyzzz_1[i] * fe_0 + g_xyy_xyyzzz_1[i] * pa_x[i];

        g_xxyy_xyzzzz_0[i] = g_yy_xyzzzz_0[i] * fbe_0 - g_yy_xyzzzz_1[i] * fz_be_0 + g_xyy_yzzzz_1[i] * fe_0 + g_xyy_xyzzzz_1[i] * pa_x[i];

        g_xxyy_xzzzzz_0[i] = g_xx_xzzzzz_0[i] * fbe_0 - g_xx_xzzzzz_1[i] * fz_be_0 + g_xxy_xzzzzz_1[i] * pa_y[i];

        g_xxyy_yyyyyy_0[i] = g_yy_yyyyyy_0[i] * fbe_0 - g_yy_yyyyyy_1[i] * fz_be_0 + g_xyy_yyyyyy_1[i] * pa_x[i];

        g_xxyy_yyyyyz_0[i] = g_yy_yyyyyz_0[i] * fbe_0 - g_yy_yyyyyz_1[i] * fz_be_0 + g_xyy_yyyyyz_1[i] * pa_x[i];

        g_xxyy_yyyyzz_0[i] = g_yy_yyyyzz_0[i] * fbe_0 - g_yy_yyyyzz_1[i] * fz_be_0 + g_xyy_yyyyzz_1[i] * pa_x[i];

        g_xxyy_yyyzzz_0[i] = g_yy_yyyzzz_0[i] * fbe_0 - g_yy_yyyzzz_1[i] * fz_be_0 + g_xyy_yyyzzz_1[i] * pa_x[i];

        g_xxyy_yyzzzz_0[i] = g_yy_yyzzzz_0[i] * fbe_0 - g_yy_yyzzzz_1[i] * fz_be_0 + g_xyy_yyzzzz_1[i] * pa_x[i];

        g_xxyy_yzzzzz_0[i] = g_yy_yzzzzz_0[i] * fbe_0 - g_yy_yzzzzz_1[i] * fz_be_0 + g_xyy_yzzzzz_1[i] * pa_x[i];

        g_xxyy_zzzzzz_0[i] = g_yy_zzzzzz_0[i] * fbe_0 - g_yy_zzzzzz_1[i] * fz_be_0 + g_xyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 112-140 components of targeted buffer : GI

    auto g_xxyz_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 112);

    auto g_xxyz_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 113);

    auto g_xxyz_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 114);

    auto g_xxyz_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 115);

    auto g_xxyz_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 116);

    auto g_xxyz_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 117);

    auto g_xxyz_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 118);

    auto g_xxyz_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 119);

    auto g_xxyz_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 120);

    auto g_xxyz_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 121);

    auto g_xxyz_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 122);

    auto g_xxyz_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 123);

    auto g_xxyz_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 124);

    auto g_xxyz_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 125);

    auto g_xxyz_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 126);

    auto g_xxyz_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 127);

    auto g_xxyz_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 128);

    auto g_xxyz_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 129);

    auto g_xxyz_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 130);

    auto g_xxyz_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 131);

    auto g_xxyz_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 132);

    auto g_xxyz_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 133);

    auto g_xxyz_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 134);

    auto g_xxyz_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 135);

    auto g_xxyz_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 136);

    auto g_xxyz_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 137);

    auto g_xxyz_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 138);

    auto g_xxyz_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 139);

    #pragma omp simd aligned(g_xxy_xxxxxy_1, g_xxy_xxxxyy_1, g_xxy_xxxyyy_1, g_xxy_xxyyyy_1, g_xxy_xyyyyy_1, g_xxy_yyyyyy_1, g_xxyz_xxxxxx_0, g_xxyz_xxxxxy_0, g_xxyz_xxxxxz_0, g_xxyz_xxxxyy_0, g_xxyz_xxxxyz_0, g_xxyz_xxxxzz_0, g_xxyz_xxxyyy_0, g_xxyz_xxxyyz_0, g_xxyz_xxxyzz_0, g_xxyz_xxxzzz_0, g_xxyz_xxyyyy_0, g_xxyz_xxyyyz_0, g_xxyz_xxyyzz_0, g_xxyz_xxyzzz_0, g_xxyz_xxzzzz_0, g_xxyz_xyyyyy_0, g_xxyz_xyyyyz_0, g_xxyz_xyyyzz_0, g_xxyz_xyyzzz_0, g_xxyz_xyzzzz_0, g_xxyz_xzzzzz_0, g_xxyz_yyyyyy_0, g_xxyz_yyyyyz_0, g_xxyz_yyyyzz_0, g_xxyz_yyyzzz_0, g_xxyz_yyzzzz_0, g_xxyz_yzzzzz_0, g_xxyz_zzzzzz_0, g_xxz_xxxxxx_1, g_xxz_xxxxxz_1, g_xxz_xxxxyz_1, g_xxz_xxxxz_1, g_xxz_xxxxzz_1, g_xxz_xxxyyz_1, g_xxz_xxxyz_1, g_xxz_xxxyzz_1, g_xxz_xxxzz_1, g_xxz_xxxzzz_1, g_xxz_xxyyyz_1, g_xxz_xxyyz_1, g_xxz_xxyyzz_1, g_xxz_xxyzz_1, g_xxz_xxyzzz_1, g_xxz_xxzzz_1, g_xxz_xxzzzz_1, g_xxz_xyyyyz_1, g_xxz_xyyyz_1, g_xxz_xyyyzz_1, g_xxz_xyyzz_1, g_xxz_xyyzzz_1, g_xxz_xyzzz_1, g_xxz_xyzzzz_1, g_xxz_xzzzz_1, g_xxz_xzzzzz_1, g_xxz_yyyyyz_1, g_xxz_yyyyz_1, g_xxz_yyyyzz_1, g_xxz_yyyzz_1, g_xxz_yyyzzz_1, g_xxz_yyzzz_1, g_xxz_yyzzzz_1, g_xxz_yzzzz_1, g_xxz_yzzzzz_1, g_xxz_zzzzz_1, g_xxz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyz_xxxxxx_0[i] = g_xxz_xxxxxx_1[i] * pa_y[i];

        g_xxyz_xxxxxy_0[i] = g_xxy_xxxxxy_1[i] * pa_z[i];

        g_xxyz_xxxxxz_0[i] = g_xxz_xxxxxz_1[i] * pa_y[i];

        g_xxyz_xxxxyy_0[i] = g_xxy_xxxxyy_1[i] * pa_z[i];

        g_xxyz_xxxxyz_0[i] = g_xxz_xxxxz_1[i] * fe_0 + g_xxz_xxxxyz_1[i] * pa_y[i];

        g_xxyz_xxxxzz_0[i] = g_xxz_xxxxzz_1[i] * pa_y[i];

        g_xxyz_xxxyyy_0[i] = g_xxy_xxxyyy_1[i] * pa_z[i];

        g_xxyz_xxxyyz_0[i] = 2.0 * g_xxz_xxxyz_1[i] * fe_0 + g_xxz_xxxyyz_1[i] * pa_y[i];

        g_xxyz_xxxyzz_0[i] = g_xxz_xxxzz_1[i] * fe_0 + g_xxz_xxxyzz_1[i] * pa_y[i];

        g_xxyz_xxxzzz_0[i] = g_xxz_xxxzzz_1[i] * pa_y[i];

        g_xxyz_xxyyyy_0[i] = g_xxy_xxyyyy_1[i] * pa_z[i];

        g_xxyz_xxyyyz_0[i] = 3.0 * g_xxz_xxyyz_1[i] * fe_0 + g_xxz_xxyyyz_1[i] * pa_y[i];

        g_xxyz_xxyyzz_0[i] = 2.0 * g_xxz_xxyzz_1[i] * fe_0 + g_xxz_xxyyzz_1[i] * pa_y[i];

        g_xxyz_xxyzzz_0[i] = g_xxz_xxzzz_1[i] * fe_0 + g_xxz_xxyzzz_1[i] * pa_y[i];

        g_xxyz_xxzzzz_0[i] = g_xxz_xxzzzz_1[i] * pa_y[i];

        g_xxyz_xyyyyy_0[i] = g_xxy_xyyyyy_1[i] * pa_z[i];

        g_xxyz_xyyyyz_0[i] = 4.0 * g_xxz_xyyyz_1[i] * fe_0 + g_xxz_xyyyyz_1[i] * pa_y[i];

        g_xxyz_xyyyzz_0[i] = 3.0 * g_xxz_xyyzz_1[i] * fe_0 + g_xxz_xyyyzz_1[i] * pa_y[i];

        g_xxyz_xyyzzz_0[i] = 2.0 * g_xxz_xyzzz_1[i] * fe_0 + g_xxz_xyyzzz_1[i] * pa_y[i];

        g_xxyz_xyzzzz_0[i] = g_xxz_xzzzz_1[i] * fe_0 + g_xxz_xyzzzz_1[i] * pa_y[i];

        g_xxyz_xzzzzz_0[i] = g_xxz_xzzzzz_1[i] * pa_y[i];

        g_xxyz_yyyyyy_0[i] = g_xxy_yyyyyy_1[i] * pa_z[i];

        g_xxyz_yyyyyz_0[i] = 5.0 * g_xxz_yyyyz_1[i] * fe_0 + g_xxz_yyyyyz_1[i] * pa_y[i];

        g_xxyz_yyyyzz_0[i] = 4.0 * g_xxz_yyyzz_1[i] * fe_0 + g_xxz_yyyyzz_1[i] * pa_y[i];

        g_xxyz_yyyzzz_0[i] = 3.0 * g_xxz_yyzzz_1[i] * fe_0 + g_xxz_yyyzzz_1[i] * pa_y[i];

        g_xxyz_yyzzzz_0[i] = 2.0 * g_xxz_yzzzz_1[i] * fe_0 + g_xxz_yyzzzz_1[i] * pa_y[i];

        g_xxyz_yzzzzz_0[i] = g_xxz_zzzzz_1[i] * fe_0 + g_xxz_yzzzzz_1[i] * pa_y[i];

        g_xxyz_zzzzzz_0[i] = g_xxz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : GI

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

    #pragma omp simd aligned(g_xx_xxxxxx_0, g_xx_xxxxxx_1, g_xx_xxxxxy_0, g_xx_xxxxxy_1, g_xx_xxxxyy_0, g_xx_xxxxyy_1, g_xx_xxxyyy_0, g_xx_xxxyyy_1, g_xx_xxyyyy_0, g_xx_xxyyyy_1, g_xx_xyyyyy_0, g_xx_xyyyyy_1, g_xxz_xxxxxx_1, g_xxz_xxxxxy_1, g_xxz_xxxxyy_1, g_xxz_xxxyyy_1, g_xxz_xxyyyy_1, g_xxz_xyyyyy_1, g_xxzz_xxxxxx_0, g_xxzz_xxxxxy_0, g_xxzz_xxxxxz_0, g_xxzz_xxxxyy_0, g_xxzz_xxxxyz_0, g_xxzz_xxxxzz_0, g_xxzz_xxxyyy_0, g_xxzz_xxxyyz_0, g_xxzz_xxxyzz_0, g_xxzz_xxxzzz_0, g_xxzz_xxyyyy_0, g_xxzz_xxyyyz_0, g_xxzz_xxyyzz_0, g_xxzz_xxyzzz_0, g_xxzz_xxzzzz_0, g_xxzz_xyyyyy_0, g_xxzz_xyyyyz_0, g_xxzz_xyyyzz_0, g_xxzz_xyyzzz_0, g_xxzz_xyzzzz_0, g_xxzz_xzzzzz_0, g_xxzz_yyyyyy_0, g_xxzz_yyyyyz_0, g_xxzz_yyyyzz_0, g_xxzz_yyyzzz_0, g_xxzz_yyzzzz_0, g_xxzz_yzzzzz_0, g_xxzz_zzzzzz_0, g_xzz_xxxxxz_1, g_xzz_xxxxyz_1, g_xzz_xxxxz_1, g_xzz_xxxxzz_1, g_xzz_xxxyyz_1, g_xzz_xxxyz_1, g_xzz_xxxyzz_1, g_xzz_xxxzz_1, g_xzz_xxxzzz_1, g_xzz_xxyyyz_1, g_xzz_xxyyz_1, g_xzz_xxyyzz_1, g_xzz_xxyzz_1, g_xzz_xxyzzz_1, g_xzz_xxzzz_1, g_xzz_xxzzzz_1, g_xzz_xyyyyz_1, g_xzz_xyyyz_1, g_xzz_xyyyzz_1, g_xzz_xyyzz_1, g_xzz_xyyzzz_1, g_xzz_xyzzz_1, g_xzz_xyzzzz_1, g_xzz_xzzzz_1, g_xzz_xzzzzz_1, g_xzz_yyyyyy_1, g_xzz_yyyyyz_1, g_xzz_yyyyz_1, g_xzz_yyyyzz_1, g_xzz_yyyzz_1, g_xzz_yyyzzz_1, g_xzz_yyzzz_1, g_xzz_yyzzzz_1, g_xzz_yzzzz_1, g_xzz_yzzzzz_1, g_xzz_zzzzz_1, g_xzz_zzzzzz_1, g_zz_xxxxxz_0, g_zz_xxxxxz_1, g_zz_xxxxyz_0, g_zz_xxxxyz_1, g_zz_xxxxzz_0, g_zz_xxxxzz_1, g_zz_xxxyyz_0, g_zz_xxxyyz_1, g_zz_xxxyzz_0, g_zz_xxxyzz_1, g_zz_xxxzzz_0, g_zz_xxxzzz_1, g_zz_xxyyyz_0, g_zz_xxyyyz_1, g_zz_xxyyzz_0, g_zz_xxyyzz_1, g_zz_xxyzzz_0, g_zz_xxyzzz_1, g_zz_xxzzzz_0, g_zz_xxzzzz_1, g_zz_xyyyyz_0, g_zz_xyyyyz_1, g_zz_xyyyzz_0, g_zz_xyyyzz_1, g_zz_xyyzzz_0, g_zz_xyyzzz_1, g_zz_xyzzzz_0, g_zz_xyzzzz_1, g_zz_xzzzzz_0, g_zz_xzzzzz_1, g_zz_yyyyyy_0, g_zz_yyyyyy_1, g_zz_yyyyyz_0, g_zz_yyyyyz_1, g_zz_yyyyzz_0, g_zz_yyyyzz_1, g_zz_yyyzzz_0, g_zz_yyyzzz_1, g_zz_yyzzzz_0, g_zz_yyzzzz_1, g_zz_yzzzzz_0, g_zz_yzzzzz_1, g_zz_zzzzzz_0, g_zz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzz_xxxxxx_0[i] = g_xx_xxxxxx_0[i] * fbe_0 - g_xx_xxxxxx_1[i] * fz_be_0 + g_xxz_xxxxxx_1[i] * pa_z[i];

        g_xxzz_xxxxxy_0[i] = g_xx_xxxxxy_0[i] * fbe_0 - g_xx_xxxxxy_1[i] * fz_be_0 + g_xxz_xxxxxy_1[i] * pa_z[i];

        g_xxzz_xxxxxz_0[i] = g_zz_xxxxxz_0[i] * fbe_0 - g_zz_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xzz_xxxxz_1[i] * fe_0 + g_xzz_xxxxxz_1[i] * pa_x[i];

        g_xxzz_xxxxyy_0[i] = g_xx_xxxxyy_0[i] * fbe_0 - g_xx_xxxxyy_1[i] * fz_be_0 + g_xxz_xxxxyy_1[i] * pa_z[i];

        g_xxzz_xxxxyz_0[i] = g_zz_xxxxyz_0[i] * fbe_0 - g_zz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xzz_xxxyz_1[i] * fe_0 + g_xzz_xxxxyz_1[i] * pa_x[i];

        g_xxzz_xxxxzz_0[i] = g_zz_xxxxzz_0[i] * fbe_0 - g_zz_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xzz_xxxzz_1[i] * fe_0 + g_xzz_xxxxzz_1[i] * pa_x[i];

        g_xxzz_xxxyyy_0[i] = g_xx_xxxyyy_0[i] * fbe_0 - g_xx_xxxyyy_1[i] * fz_be_0 + g_xxz_xxxyyy_1[i] * pa_z[i];

        g_xxzz_xxxyyz_0[i] = g_zz_xxxyyz_0[i] * fbe_0 - g_zz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xzz_xxyyz_1[i] * fe_0 + g_xzz_xxxyyz_1[i] * pa_x[i];

        g_xxzz_xxxyzz_0[i] = g_zz_xxxyzz_0[i] * fbe_0 - g_zz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xzz_xxyzz_1[i] * fe_0 + g_xzz_xxxyzz_1[i] * pa_x[i];

        g_xxzz_xxxzzz_0[i] = g_zz_xxxzzz_0[i] * fbe_0 - g_zz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xzz_xxzzz_1[i] * fe_0 + g_xzz_xxxzzz_1[i] * pa_x[i];

        g_xxzz_xxyyyy_0[i] = g_xx_xxyyyy_0[i] * fbe_0 - g_xx_xxyyyy_1[i] * fz_be_0 + g_xxz_xxyyyy_1[i] * pa_z[i];

        g_xxzz_xxyyyz_0[i] = g_zz_xxyyyz_0[i] * fbe_0 - g_zz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xzz_xyyyz_1[i] * fe_0 + g_xzz_xxyyyz_1[i] * pa_x[i];

        g_xxzz_xxyyzz_0[i] = g_zz_xxyyzz_0[i] * fbe_0 - g_zz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xzz_xyyzz_1[i] * fe_0 + g_xzz_xxyyzz_1[i] * pa_x[i];

        g_xxzz_xxyzzz_0[i] = g_zz_xxyzzz_0[i] * fbe_0 - g_zz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xzz_xyzzz_1[i] * fe_0 + g_xzz_xxyzzz_1[i] * pa_x[i];

        g_xxzz_xxzzzz_0[i] = g_zz_xxzzzz_0[i] * fbe_0 - g_zz_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xzz_xzzzz_1[i] * fe_0 + g_xzz_xxzzzz_1[i] * pa_x[i];

        g_xxzz_xyyyyy_0[i] = g_xx_xyyyyy_0[i] * fbe_0 - g_xx_xyyyyy_1[i] * fz_be_0 + g_xxz_xyyyyy_1[i] * pa_z[i];

        g_xxzz_xyyyyz_0[i] = g_zz_xyyyyz_0[i] * fbe_0 - g_zz_xyyyyz_1[i] * fz_be_0 + g_xzz_yyyyz_1[i] * fe_0 + g_xzz_xyyyyz_1[i] * pa_x[i];

        g_xxzz_xyyyzz_0[i] = g_zz_xyyyzz_0[i] * fbe_0 - g_zz_xyyyzz_1[i] * fz_be_0 + g_xzz_yyyzz_1[i] * fe_0 + g_xzz_xyyyzz_1[i] * pa_x[i];

        g_xxzz_xyyzzz_0[i] = g_zz_xyyzzz_0[i] * fbe_0 - g_zz_xyyzzz_1[i] * fz_be_0 + g_xzz_yyzzz_1[i] * fe_0 + g_xzz_xyyzzz_1[i] * pa_x[i];

        g_xxzz_xyzzzz_0[i] = g_zz_xyzzzz_0[i] * fbe_0 - g_zz_xyzzzz_1[i] * fz_be_0 + g_xzz_yzzzz_1[i] * fe_0 + g_xzz_xyzzzz_1[i] * pa_x[i];

        g_xxzz_xzzzzz_0[i] = g_zz_xzzzzz_0[i] * fbe_0 - g_zz_xzzzzz_1[i] * fz_be_0 + g_xzz_zzzzz_1[i] * fe_0 + g_xzz_xzzzzz_1[i] * pa_x[i];

        g_xxzz_yyyyyy_0[i] = g_zz_yyyyyy_0[i] * fbe_0 - g_zz_yyyyyy_1[i] * fz_be_0 + g_xzz_yyyyyy_1[i] * pa_x[i];

        g_xxzz_yyyyyz_0[i] = g_zz_yyyyyz_0[i] * fbe_0 - g_zz_yyyyyz_1[i] * fz_be_0 + g_xzz_yyyyyz_1[i] * pa_x[i];

        g_xxzz_yyyyzz_0[i] = g_zz_yyyyzz_0[i] * fbe_0 - g_zz_yyyyzz_1[i] * fz_be_0 + g_xzz_yyyyzz_1[i] * pa_x[i];

        g_xxzz_yyyzzz_0[i] = g_zz_yyyzzz_0[i] * fbe_0 - g_zz_yyyzzz_1[i] * fz_be_0 + g_xzz_yyyzzz_1[i] * pa_x[i];

        g_xxzz_yyzzzz_0[i] = g_zz_yyzzzz_0[i] * fbe_0 - g_zz_yyzzzz_1[i] * fz_be_0 + g_xzz_yyzzzz_1[i] * pa_x[i];

        g_xxzz_yzzzzz_0[i] = g_zz_yzzzzz_0[i] * fbe_0 - g_zz_yzzzzz_1[i] * fz_be_0 + g_xzz_yzzzzz_1[i] * pa_x[i];

        g_xxzz_zzzzzz_0[i] = g_zz_zzzzzz_0[i] * fbe_0 - g_zz_zzzzzz_1[i] * fz_be_0 + g_xzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : GI

    auto g_xyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 168);

    auto g_xyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 169);

    auto g_xyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 170);

    auto g_xyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 171);

    auto g_xyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 172);

    auto g_xyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 173);

    auto g_xyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 174);

    auto g_xyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 175);

    auto g_xyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 176);

    auto g_xyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 177);

    auto g_xyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 178);

    auto g_xyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 179);

    auto g_xyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 180);

    auto g_xyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 181);

    auto g_xyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 182);

    auto g_xyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 183);

    auto g_xyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 184);

    auto g_xyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 185);

    auto g_xyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 186);

    auto g_xyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 187);

    auto g_xyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 188);

    auto g_xyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 189);

    auto g_xyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 190);

    auto g_xyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 191);

    auto g_xyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 192);

    auto g_xyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 193);

    auto g_xyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 194);

    auto g_xyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 195);

    #pragma omp simd aligned(g_xyyy_xxxxxx_0, g_xyyy_xxxxxy_0, g_xyyy_xxxxxz_0, g_xyyy_xxxxyy_0, g_xyyy_xxxxyz_0, g_xyyy_xxxxzz_0, g_xyyy_xxxyyy_0, g_xyyy_xxxyyz_0, g_xyyy_xxxyzz_0, g_xyyy_xxxzzz_0, g_xyyy_xxyyyy_0, g_xyyy_xxyyyz_0, g_xyyy_xxyyzz_0, g_xyyy_xxyzzz_0, g_xyyy_xxzzzz_0, g_xyyy_xyyyyy_0, g_xyyy_xyyyyz_0, g_xyyy_xyyyzz_0, g_xyyy_xyyzzz_0, g_xyyy_xyzzzz_0, g_xyyy_xzzzzz_0, g_xyyy_yyyyyy_0, g_xyyy_yyyyyz_0, g_xyyy_yyyyzz_0, g_xyyy_yyyzzz_0, g_xyyy_yyzzzz_0, g_xyyy_yzzzzz_0, g_xyyy_zzzzzz_0, g_yyy_xxxxx_1, g_yyy_xxxxxx_1, g_yyy_xxxxxy_1, g_yyy_xxxxxz_1, g_yyy_xxxxy_1, g_yyy_xxxxyy_1, g_yyy_xxxxyz_1, g_yyy_xxxxz_1, g_yyy_xxxxzz_1, g_yyy_xxxyy_1, g_yyy_xxxyyy_1, g_yyy_xxxyyz_1, g_yyy_xxxyz_1, g_yyy_xxxyzz_1, g_yyy_xxxzz_1, g_yyy_xxxzzz_1, g_yyy_xxyyy_1, g_yyy_xxyyyy_1, g_yyy_xxyyyz_1, g_yyy_xxyyz_1, g_yyy_xxyyzz_1, g_yyy_xxyzz_1, g_yyy_xxyzzz_1, g_yyy_xxzzz_1, g_yyy_xxzzzz_1, g_yyy_xyyyy_1, g_yyy_xyyyyy_1, g_yyy_xyyyyz_1, g_yyy_xyyyz_1, g_yyy_xyyyzz_1, g_yyy_xyyzz_1, g_yyy_xyyzzz_1, g_yyy_xyzzz_1, g_yyy_xyzzzz_1, g_yyy_xzzzz_1, g_yyy_xzzzzz_1, g_yyy_yyyyy_1, g_yyy_yyyyyy_1, g_yyy_yyyyyz_1, g_yyy_yyyyz_1, g_yyy_yyyyzz_1, g_yyy_yyyzz_1, g_yyy_yyyzzz_1, g_yyy_yyzzz_1, g_yyy_yyzzzz_1, g_yyy_yzzzz_1, g_yyy_yzzzzz_1, g_yyy_zzzzz_1, g_yyy_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyy_xxxxxx_0[i] = 6.0 * g_yyy_xxxxx_1[i] * fe_0 + g_yyy_xxxxxx_1[i] * pa_x[i];

        g_xyyy_xxxxxy_0[i] = 5.0 * g_yyy_xxxxy_1[i] * fe_0 + g_yyy_xxxxxy_1[i] * pa_x[i];

        g_xyyy_xxxxxz_0[i] = 5.0 * g_yyy_xxxxz_1[i] * fe_0 + g_yyy_xxxxxz_1[i] * pa_x[i];

        g_xyyy_xxxxyy_0[i] = 4.0 * g_yyy_xxxyy_1[i] * fe_0 + g_yyy_xxxxyy_1[i] * pa_x[i];

        g_xyyy_xxxxyz_0[i] = 4.0 * g_yyy_xxxyz_1[i] * fe_0 + g_yyy_xxxxyz_1[i] * pa_x[i];

        g_xyyy_xxxxzz_0[i] = 4.0 * g_yyy_xxxzz_1[i] * fe_0 + g_yyy_xxxxzz_1[i] * pa_x[i];

        g_xyyy_xxxyyy_0[i] = 3.0 * g_yyy_xxyyy_1[i] * fe_0 + g_yyy_xxxyyy_1[i] * pa_x[i];

        g_xyyy_xxxyyz_0[i] = 3.0 * g_yyy_xxyyz_1[i] * fe_0 + g_yyy_xxxyyz_1[i] * pa_x[i];

        g_xyyy_xxxyzz_0[i] = 3.0 * g_yyy_xxyzz_1[i] * fe_0 + g_yyy_xxxyzz_1[i] * pa_x[i];

        g_xyyy_xxxzzz_0[i] = 3.0 * g_yyy_xxzzz_1[i] * fe_0 + g_yyy_xxxzzz_1[i] * pa_x[i];

        g_xyyy_xxyyyy_0[i] = 2.0 * g_yyy_xyyyy_1[i] * fe_0 + g_yyy_xxyyyy_1[i] * pa_x[i];

        g_xyyy_xxyyyz_0[i] = 2.0 * g_yyy_xyyyz_1[i] * fe_0 + g_yyy_xxyyyz_1[i] * pa_x[i];

        g_xyyy_xxyyzz_0[i] = 2.0 * g_yyy_xyyzz_1[i] * fe_0 + g_yyy_xxyyzz_1[i] * pa_x[i];

        g_xyyy_xxyzzz_0[i] = 2.0 * g_yyy_xyzzz_1[i] * fe_0 + g_yyy_xxyzzz_1[i] * pa_x[i];

        g_xyyy_xxzzzz_0[i] = 2.0 * g_yyy_xzzzz_1[i] * fe_0 + g_yyy_xxzzzz_1[i] * pa_x[i];

        g_xyyy_xyyyyy_0[i] = g_yyy_yyyyy_1[i] * fe_0 + g_yyy_xyyyyy_1[i] * pa_x[i];

        g_xyyy_xyyyyz_0[i] = g_yyy_yyyyz_1[i] * fe_0 + g_yyy_xyyyyz_1[i] * pa_x[i];

        g_xyyy_xyyyzz_0[i] = g_yyy_yyyzz_1[i] * fe_0 + g_yyy_xyyyzz_1[i] * pa_x[i];

        g_xyyy_xyyzzz_0[i] = g_yyy_yyzzz_1[i] * fe_0 + g_yyy_xyyzzz_1[i] * pa_x[i];

        g_xyyy_xyzzzz_0[i] = g_yyy_yzzzz_1[i] * fe_0 + g_yyy_xyzzzz_1[i] * pa_x[i];

        g_xyyy_xzzzzz_0[i] = g_yyy_zzzzz_1[i] * fe_0 + g_yyy_xzzzzz_1[i] * pa_x[i];

        g_xyyy_yyyyyy_0[i] = g_yyy_yyyyyy_1[i] * pa_x[i];

        g_xyyy_yyyyyz_0[i] = g_yyy_yyyyyz_1[i] * pa_x[i];

        g_xyyy_yyyyzz_0[i] = g_yyy_yyyyzz_1[i] * pa_x[i];

        g_xyyy_yyyzzz_0[i] = g_yyy_yyyzzz_1[i] * pa_x[i];

        g_xyyy_yyzzzz_0[i] = g_yyy_yyzzzz_1[i] * pa_x[i];

        g_xyyy_yzzzzz_0[i] = g_yyy_yzzzzz_1[i] * pa_x[i];

        g_xyyy_zzzzzz_0[i] = g_yyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 196-224 components of targeted buffer : GI

    auto g_xyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 196);

    auto g_xyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 197);

    auto g_xyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 198);

    auto g_xyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 199);

    auto g_xyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 200);

    auto g_xyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 201);

    auto g_xyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 202);

    auto g_xyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 203);

    auto g_xyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 204);

    auto g_xyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 205);

    auto g_xyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 206);

    auto g_xyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 207);

    auto g_xyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 208);

    auto g_xyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 209);

    auto g_xyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 210);

    auto g_xyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 211);

    auto g_xyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 212);

    auto g_xyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 213);

    auto g_xyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 214);

    auto g_xyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 215);

    auto g_xyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 216);

    auto g_xyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 217);

    auto g_xyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 218);

    auto g_xyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 219);

    auto g_xyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 220);

    auto g_xyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 221);

    auto g_xyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 222);

    auto g_xyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 223);

    #pragma omp simd aligned(g_xyy_xxxxxx_1, g_xyy_xxxxxy_1, g_xyy_xxxxyy_1, g_xyy_xxxyyy_1, g_xyy_xxyyyy_1, g_xyy_xyyyyy_1, g_xyyz_xxxxxx_0, g_xyyz_xxxxxy_0, g_xyyz_xxxxxz_0, g_xyyz_xxxxyy_0, g_xyyz_xxxxyz_0, g_xyyz_xxxxzz_0, g_xyyz_xxxyyy_0, g_xyyz_xxxyyz_0, g_xyyz_xxxyzz_0, g_xyyz_xxxzzz_0, g_xyyz_xxyyyy_0, g_xyyz_xxyyyz_0, g_xyyz_xxyyzz_0, g_xyyz_xxyzzz_0, g_xyyz_xxzzzz_0, g_xyyz_xyyyyy_0, g_xyyz_xyyyyz_0, g_xyyz_xyyyzz_0, g_xyyz_xyyzzz_0, g_xyyz_xyzzzz_0, g_xyyz_xzzzzz_0, g_xyyz_yyyyyy_0, g_xyyz_yyyyyz_0, g_xyyz_yyyyzz_0, g_xyyz_yyyzzz_0, g_xyyz_yyzzzz_0, g_xyyz_yzzzzz_0, g_xyyz_zzzzzz_0, g_yyz_xxxxxz_1, g_yyz_xxxxyz_1, g_yyz_xxxxz_1, g_yyz_xxxxzz_1, g_yyz_xxxyyz_1, g_yyz_xxxyz_1, g_yyz_xxxyzz_1, g_yyz_xxxzz_1, g_yyz_xxxzzz_1, g_yyz_xxyyyz_1, g_yyz_xxyyz_1, g_yyz_xxyyzz_1, g_yyz_xxyzz_1, g_yyz_xxyzzz_1, g_yyz_xxzzz_1, g_yyz_xxzzzz_1, g_yyz_xyyyyz_1, g_yyz_xyyyz_1, g_yyz_xyyyzz_1, g_yyz_xyyzz_1, g_yyz_xyyzzz_1, g_yyz_xyzzz_1, g_yyz_xyzzzz_1, g_yyz_xzzzz_1, g_yyz_xzzzzz_1, g_yyz_yyyyyy_1, g_yyz_yyyyyz_1, g_yyz_yyyyz_1, g_yyz_yyyyzz_1, g_yyz_yyyzz_1, g_yyz_yyyzzz_1, g_yyz_yyzzz_1, g_yyz_yyzzzz_1, g_yyz_yzzzz_1, g_yyz_yzzzzz_1, g_yyz_zzzzz_1, g_yyz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyz_xxxxxx_0[i] = g_xyy_xxxxxx_1[i] * pa_z[i];

        g_xyyz_xxxxxy_0[i] = g_xyy_xxxxxy_1[i] * pa_z[i];

        g_xyyz_xxxxxz_0[i] = 5.0 * g_yyz_xxxxz_1[i] * fe_0 + g_yyz_xxxxxz_1[i] * pa_x[i];

        g_xyyz_xxxxyy_0[i] = g_xyy_xxxxyy_1[i] * pa_z[i];

        g_xyyz_xxxxyz_0[i] = 4.0 * g_yyz_xxxyz_1[i] * fe_0 + g_yyz_xxxxyz_1[i] * pa_x[i];

        g_xyyz_xxxxzz_0[i] = 4.0 * g_yyz_xxxzz_1[i] * fe_0 + g_yyz_xxxxzz_1[i] * pa_x[i];

        g_xyyz_xxxyyy_0[i] = g_xyy_xxxyyy_1[i] * pa_z[i];

        g_xyyz_xxxyyz_0[i] = 3.0 * g_yyz_xxyyz_1[i] * fe_0 + g_yyz_xxxyyz_1[i] * pa_x[i];

        g_xyyz_xxxyzz_0[i] = 3.0 * g_yyz_xxyzz_1[i] * fe_0 + g_yyz_xxxyzz_1[i] * pa_x[i];

        g_xyyz_xxxzzz_0[i] = 3.0 * g_yyz_xxzzz_1[i] * fe_0 + g_yyz_xxxzzz_1[i] * pa_x[i];

        g_xyyz_xxyyyy_0[i] = g_xyy_xxyyyy_1[i] * pa_z[i];

        g_xyyz_xxyyyz_0[i] = 2.0 * g_yyz_xyyyz_1[i] * fe_0 + g_yyz_xxyyyz_1[i] * pa_x[i];

        g_xyyz_xxyyzz_0[i] = 2.0 * g_yyz_xyyzz_1[i] * fe_0 + g_yyz_xxyyzz_1[i] * pa_x[i];

        g_xyyz_xxyzzz_0[i] = 2.0 * g_yyz_xyzzz_1[i] * fe_0 + g_yyz_xxyzzz_1[i] * pa_x[i];

        g_xyyz_xxzzzz_0[i] = 2.0 * g_yyz_xzzzz_1[i] * fe_0 + g_yyz_xxzzzz_1[i] * pa_x[i];

        g_xyyz_xyyyyy_0[i] = g_xyy_xyyyyy_1[i] * pa_z[i];

        g_xyyz_xyyyyz_0[i] = g_yyz_yyyyz_1[i] * fe_0 + g_yyz_xyyyyz_1[i] * pa_x[i];

        g_xyyz_xyyyzz_0[i] = g_yyz_yyyzz_1[i] * fe_0 + g_yyz_xyyyzz_1[i] * pa_x[i];

        g_xyyz_xyyzzz_0[i] = g_yyz_yyzzz_1[i] * fe_0 + g_yyz_xyyzzz_1[i] * pa_x[i];

        g_xyyz_xyzzzz_0[i] = g_yyz_yzzzz_1[i] * fe_0 + g_yyz_xyzzzz_1[i] * pa_x[i];

        g_xyyz_xzzzzz_0[i] = g_yyz_zzzzz_1[i] * fe_0 + g_yyz_xzzzzz_1[i] * pa_x[i];

        g_xyyz_yyyyyy_0[i] = g_yyz_yyyyyy_1[i] * pa_x[i];

        g_xyyz_yyyyyz_0[i] = g_yyz_yyyyyz_1[i] * pa_x[i];

        g_xyyz_yyyyzz_0[i] = g_yyz_yyyyzz_1[i] * pa_x[i];

        g_xyyz_yyyzzz_0[i] = g_yyz_yyyzzz_1[i] * pa_x[i];

        g_xyyz_yyzzzz_0[i] = g_yyz_yyzzzz_1[i] * pa_x[i];

        g_xyyz_yzzzzz_0[i] = g_yyz_yzzzzz_1[i] * pa_x[i];

        g_xyyz_zzzzzz_0[i] = g_yyz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 224-252 components of targeted buffer : GI

    auto g_xyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 224);

    auto g_xyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 225);

    auto g_xyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 226);

    auto g_xyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 227);

    auto g_xyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 228);

    auto g_xyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 229);

    auto g_xyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 230);

    auto g_xyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 231);

    auto g_xyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 232);

    auto g_xyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 233);

    auto g_xyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 234);

    auto g_xyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 235);

    auto g_xyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 236);

    auto g_xyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 237);

    auto g_xyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 238);

    auto g_xyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 239);

    auto g_xyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 240);

    auto g_xyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 241);

    auto g_xyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 242);

    auto g_xyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 243);

    auto g_xyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 244);

    auto g_xyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 245);

    auto g_xyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 246);

    auto g_xyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 247);

    auto g_xyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 248);

    auto g_xyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 249);

    auto g_xyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 250);

    auto g_xyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 251);

    #pragma omp simd aligned(g_xyzz_xxxxxx_0, g_xyzz_xxxxxy_0, g_xyzz_xxxxxz_0, g_xyzz_xxxxyy_0, g_xyzz_xxxxyz_0, g_xyzz_xxxxzz_0, g_xyzz_xxxyyy_0, g_xyzz_xxxyyz_0, g_xyzz_xxxyzz_0, g_xyzz_xxxzzz_0, g_xyzz_xxyyyy_0, g_xyzz_xxyyyz_0, g_xyzz_xxyyzz_0, g_xyzz_xxyzzz_0, g_xyzz_xxzzzz_0, g_xyzz_xyyyyy_0, g_xyzz_xyyyyz_0, g_xyzz_xyyyzz_0, g_xyzz_xyyzzz_0, g_xyzz_xyzzzz_0, g_xyzz_xzzzzz_0, g_xyzz_yyyyyy_0, g_xyzz_yyyyyz_0, g_xyzz_yyyyzz_0, g_xyzz_yyyzzz_0, g_xyzz_yyzzzz_0, g_xyzz_yzzzzz_0, g_xyzz_zzzzzz_0, g_xzz_xxxxxx_1, g_xzz_xxxxxz_1, g_xzz_xxxxzz_1, g_xzz_xxxzzz_1, g_xzz_xxzzzz_1, g_xzz_xzzzzz_1, g_yzz_xxxxxy_1, g_yzz_xxxxy_1, g_yzz_xxxxyy_1, g_yzz_xxxxyz_1, g_yzz_xxxyy_1, g_yzz_xxxyyy_1, g_yzz_xxxyyz_1, g_yzz_xxxyz_1, g_yzz_xxxyzz_1, g_yzz_xxyyy_1, g_yzz_xxyyyy_1, g_yzz_xxyyyz_1, g_yzz_xxyyz_1, g_yzz_xxyyzz_1, g_yzz_xxyzz_1, g_yzz_xxyzzz_1, g_yzz_xyyyy_1, g_yzz_xyyyyy_1, g_yzz_xyyyyz_1, g_yzz_xyyyz_1, g_yzz_xyyyzz_1, g_yzz_xyyzz_1, g_yzz_xyyzzz_1, g_yzz_xyzzz_1, g_yzz_xyzzzz_1, g_yzz_yyyyy_1, g_yzz_yyyyyy_1, g_yzz_yyyyyz_1, g_yzz_yyyyz_1, g_yzz_yyyyzz_1, g_yzz_yyyzz_1, g_yzz_yyyzzz_1, g_yzz_yyzzz_1, g_yzz_yyzzzz_1, g_yzz_yzzzz_1, g_yzz_yzzzzz_1, g_yzz_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzz_xxxxxx_0[i] = g_xzz_xxxxxx_1[i] * pa_y[i];

        g_xyzz_xxxxxy_0[i] = 5.0 * g_yzz_xxxxy_1[i] * fe_0 + g_yzz_xxxxxy_1[i] * pa_x[i];

        g_xyzz_xxxxxz_0[i] = g_xzz_xxxxxz_1[i] * pa_y[i];

        g_xyzz_xxxxyy_0[i] = 4.0 * g_yzz_xxxyy_1[i] * fe_0 + g_yzz_xxxxyy_1[i] * pa_x[i];

        g_xyzz_xxxxyz_0[i] = 4.0 * g_yzz_xxxyz_1[i] * fe_0 + g_yzz_xxxxyz_1[i] * pa_x[i];

        g_xyzz_xxxxzz_0[i] = g_xzz_xxxxzz_1[i] * pa_y[i];

        g_xyzz_xxxyyy_0[i] = 3.0 * g_yzz_xxyyy_1[i] * fe_0 + g_yzz_xxxyyy_1[i] * pa_x[i];

        g_xyzz_xxxyyz_0[i] = 3.0 * g_yzz_xxyyz_1[i] * fe_0 + g_yzz_xxxyyz_1[i] * pa_x[i];

        g_xyzz_xxxyzz_0[i] = 3.0 * g_yzz_xxyzz_1[i] * fe_0 + g_yzz_xxxyzz_1[i] * pa_x[i];

        g_xyzz_xxxzzz_0[i] = g_xzz_xxxzzz_1[i] * pa_y[i];

        g_xyzz_xxyyyy_0[i] = 2.0 * g_yzz_xyyyy_1[i] * fe_0 + g_yzz_xxyyyy_1[i] * pa_x[i];

        g_xyzz_xxyyyz_0[i] = 2.0 * g_yzz_xyyyz_1[i] * fe_0 + g_yzz_xxyyyz_1[i] * pa_x[i];

        g_xyzz_xxyyzz_0[i] = 2.0 * g_yzz_xyyzz_1[i] * fe_0 + g_yzz_xxyyzz_1[i] * pa_x[i];

        g_xyzz_xxyzzz_0[i] = 2.0 * g_yzz_xyzzz_1[i] * fe_0 + g_yzz_xxyzzz_1[i] * pa_x[i];

        g_xyzz_xxzzzz_0[i] = g_xzz_xxzzzz_1[i] * pa_y[i];

        g_xyzz_xyyyyy_0[i] = g_yzz_yyyyy_1[i] * fe_0 + g_yzz_xyyyyy_1[i] * pa_x[i];

        g_xyzz_xyyyyz_0[i] = g_yzz_yyyyz_1[i] * fe_0 + g_yzz_xyyyyz_1[i] * pa_x[i];

        g_xyzz_xyyyzz_0[i] = g_yzz_yyyzz_1[i] * fe_0 + g_yzz_xyyyzz_1[i] * pa_x[i];

        g_xyzz_xyyzzz_0[i] = g_yzz_yyzzz_1[i] * fe_0 + g_yzz_xyyzzz_1[i] * pa_x[i];

        g_xyzz_xyzzzz_0[i] = g_yzz_yzzzz_1[i] * fe_0 + g_yzz_xyzzzz_1[i] * pa_x[i];

        g_xyzz_xzzzzz_0[i] = g_xzz_xzzzzz_1[i] * pa_y[i];

        g_xyzz_yyyyyy_0[i] = g_yzz_yyyyyy_1[i] * pa_x[i];

        g_xyzz_yyyyyz_0[i] = g_yzz_yyyyyz_1[i] * pa_x[i];

        g_xyzz_yyyyzz_0[i] = g_yzz_yyyyzz_1[i] * pa_x[i];

        g_xyzz_yyyzzz_0[i] = g_yzz_yyyzzz_1[i] * pa_x[i];

        g_xyzz_yyzzzz_0[i] = g_yzz_yyzzzz_1[i] * pa_x[i];

        g_xyzz_yzzzzz_0[i] = g_yzz_yzzzzz_1[i] * pa_x[i];

        g_xyzz_zzzzzz_0[i] = g_yzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 252-280 components of targeted buffer : GI

    auto g_xzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 252);

    auto g_xzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 253);

    auto g_xzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 254);

    auto g_xzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 255);

    auto g_xzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 256);

    auto g_xzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 257);

    auto g_xzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 258);

    auto g_xzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 259);

    auto g_xzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 260);

    auto g_xzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 261);

    auto g_xzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 262);

    auto g_xzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 263);

    auto g_xzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 264);

    auto g_xzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 265);

    auto g_xzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 266);

    auto g_xzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 267);

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

    #pragma omp simd aligned(g_xzzz_xxxxxx_0, g_xzzz_xxxxxy_0, g_xzzz_xxxxxz_0, g_xzzz_xxxxyy_0, g_xzzz_xxxxyz_0, g_xzzz_xxxxzz_0, g_xzzz_xxxyyy_0, g_xzzz_xxxyyz_0, g_xzzz_xxxyzz_0, g_xzzz_xxxzzz_0, g_xzzz_xxyyyy_0, g_xzzz_xxyyyz_0, g_xzzz_xxyyzz_0, g_xzzz_xxyzzz_0, g_xzzz_xxzzzz_0, g_xzzz_xyyyyy_0, g_xzzz_xyyyyz_0, g_xzzz_xyyyzz_0, g_xzzz_xyyzzz_0, g_xzzz_xyzzzz_0, g_xzzz_xzzzzz_0, g_xzzz_yyyyyy_0, g_xzzz_yyyyyz_0, g_xzzz_yyyyzz_0, g_xzzz_yyyzzz_0, g_xzzz_yyzzzz_0, g_xzzz_yzzzzz_0, g_xzzz_zzzzzz_0, g_zzz_xxxxx_1, g_zzz_xxxxxx_1, g_zzz_xxxxxy_1, g_zzz_xxxxxz_1, g_zzz_xxxxy_1, g_zzz_xxxxyy_1, g_zzz_xxxxyz_1, g_zzz_xxxxz_1, g_zzz_xxxxzz_1, g_zzz_xxxyy_1, g_zzz_xxxyyy_1, g_zzz_xxxyyz_1, g_zzz_xxxyz_1, g_zzz_xxxyzz_1, g_zzz_xxxzz_1, g_zzz_xxxzzz_1, g_zzz_xxyyy_1, g_zzz_xxyyyy_1, g_zzz_xxyyyz_1, g_zzz_xxyyz_1, g_zzz_xxyyzz_1, g_zzz_xxyzz_1, g_zzz_xxyzzz_1, g_zzz_xxzzz_1, g_zzz_xxzzzz_1, g_zzz_xyyyy_1, g_zzz_xyyyyy_1, g_zzz_xyyyyz_1, g_zzz_xyyyz_1, g_zzz_xyyyzz_1, g_zzz_xyyzz_1, g_zzz_xyyzzz_1, g_zzz_xyzzz_1, g_zzz_xyzzzz_1, g_zzz_xzzzz_1, g_zzz_xzzzzz_1, g_zzz_yyyyy_1, g_zzz_yyyyyy_1, g_zzz_yyyyyz_1, g_zzz_yyyyz_1, g_zzz_yyyyzz_1, g_zzz_yyyzz_1, g_zzz_yyyzzz_1, g_zzz_yyzzz_1, g_zzz_yyzzzz_1, g_zzz_yzzzz_1, g_zzz_yzzzzz_1, g_zzz_zzzzz_1, g_zzz_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzz_xxxxxx_0[i] = 6.0 * g_zzz_xxxxx_1[i] * fe_0 + g_zzz_xxxxxx_1[i] * pa_x[i];

        g_xzzz_xxxxxy_0[i] = 5.0 * g_zzz_xxxxy_1[i] * fe_0 + g_zzz_xxxxxy_1[i] * pa_x[i];

        g_xzzz_xxxxxz_0[i] = 5.0 * g_zzz_xxxxz_1[i] * fe_0 + g_zzz_xxxxxz_1[i] * pa_x[i];

        g_xzzz_xxxxyy_0[i] = 4.0 * g_zzz_xxxyy_1[i] * fe_0 + g_zzz_xxxxyy_1[i] * pa_x[i];

        g_xzzz_xxxxyz_0[i] = 4.0 * g_zzz_xxxyz_1[i] * fe_0 + g_zzz_xxxxyz_1[i] * pa_x[i];

        g_xzzz_xxxxzz_0[i] = 4.0 * g_zzz_xxxzz_1[i] * fe_0 + g_zzz_xxxxzz_1[i] * pa_x[i];

        g_xzzz_xxxyyy_0[i] = 3.0 * g_zzz_xxyyy_1[i] * fe_0 + g_zzz_xxxyyy_1[i] * pa_x[i];

        g_xzzz_xxxyyz_0[i] = 3.0 * g_zzz_xxyyz_1[i] * fe_0 + g_zzz_xxxyyz_1[i] * pa_x[i];

        g_xzzz_xxxyzz_0[i] = 3.0 * g_zzz_xxyzz_1[i] * fe_0 + g_zzz_xxxyzz_1[i] * pa_x[i];

        g_xzzz_xxxzzz_0[i] = 3.0 * g_zzz_xxzzz_1[i] * fe_0 + g_zzz_xxxzzz_1[i] * pa_x[i];

        g_xzzz_xxyyyy_0[i] = 2.0 * g_zzz_xyyyy_1[i] * fe_0 + g_zzz_xxyyyy_1[i] * pa_x[i];

        g_xzzz_xxyyyz_0[i] = 2.0 * g_zzz_xyyyz_1[i] * fe_0 + g_zzz_xxyyyz_1[i] * pa_x[i];

        g_xzzz_xxyyzz_0[i] = 2.0 * g_zzz_xyyzz_1[i] * fe_0 + g_zzz_xxyyzz_1[i] * pa_x[i];

        g_xzzz_xxyzzz_0[i] = 2.0 * g_zzz_xyzzz_1[i] * fe_0 + g_zzz_xxyzzz_1[i] * pa_x[i];

        g_xzzz_xxzzzz_0[i] = 2.0 * g_zzz_xzzzz_1[i] * fe_0 + g_zzz_xxzzzz_1[i] * pa_x[i];

        g_xzzz_xyyyyy_0[i] = g_zzz_yyyyy_1[i] * fe_0 + g_zzz_xyyyyy_1[i] * pa_x[i];

        g_xzzz_xyyyyz_0[i] = g_zzz_yyyyz_1[i] * fe_0 + g_zzz_xyyyyz_1[i] * pa_x[i];

        g_xzzz_xyyyzz_0[i] = g_zzz_yyyzz_1[i] * fe_0 + g_zzz_xyyyzz_1[i] * pa_x[i];

        g_xzzz_xyyzzz_0[i] = g_zzz_yyzzz_1[i] * fe_0 + g_zzz_xyyzzz_1[i] * pa_x[i];

        g_xzzz_xyzzzz_0[i] = g_zzz_yzzzz_1[i] * fe_0 + g_zzz_xyzzzz_1[i] * pa_x[i];

        g_xzzz_xzzzzz_0[i] = g_zzz_zzzzz_1[i] * fe_0 + g_zzz_xzzzzz_1[i] * pa_x[i];

        g_xzzz_yyyyyy_0[i] = g_zzz_yyyyyy_1[i] * pa_x[i];

        g_xzzz_yyyyyz_0[i] = g_zzz_yyyyyz_1[i] * pa_x[i];

        g_xzzz_yyyyzz_0[i] = g_zzz_yyyyzz_1[i] * pa_x[i];

        g_xzzz_yyyzzz_0[i] = g_zzz_yyyzzz_1[i] * pa_x[i];

        g_xzzz_yyzzzz_0[i] = g_zzz_yyzzzz_1[i] * pa_x[i];

        g_xzzz_yzzzzz_0[i] = g_zzz_yzzzzz_1[i] * pa_x[i];

        g_xzzz_zzzzzz_0[i] = g_zzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 280-308 components of targeted buffer : GI

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

    #pragma omp simd aligned(g_yy_xxxxxx_0, g_yy_xxxxxx_1, g_yy_xxxxxy_0, g_yy_xxxxxy_1, g_yy_xxxxxz_0, g_yy_xxxxxz_1, g_yy_xxxxyy_0, g_yy_xxxxyy_1, g_yy_xxxxyz_0, g_yy_xxxxyz_1, g_yy_xxxxzz_0, g_yy_xxxxzz_1, g_yy_xxxyyy_0, g_yy_xxxyyy_1, g_yy_xxxyyz_0, g_yy_xxxyyz_1, g_yy_xxxyzz_0, g_yy_xxxyzz_1, g_yy_xxxzzz_0, g_yy_xxxzzz_1, g_yy_xxyyyy_0, g_yy_xxyyyy_1, g_yy_xxyyyz_0, g_yy_xxyyyz_1, g_yy_xxyyzz_0, g_yy_xxyyzz_1, g_yy_xxyzzz_0, g_yy_xxyzzz_1, g_yy_xxzzzz_0, g_yy_xxzzzz_1, g_yy_xyyyyy_0, g_yy_xyyyyy_1, g_yy_xyyyyz_0, g_yy_xyyyyz_1, g_yy_xyyyzz_0, g_yy_xyyyzz_1, g_yy_xyyzzz_0, g_yy_xyyzzz_1, g_yy_xyzzzz_0, g_yy_xyzzzz_1, g_yy_xzzzzz_0, g_yy_xzzzzz_1, g_yy_yyyyyy_0, g_yy_yyyyyy_1, g_yy_yyyyyz_0, g_yy_yyyyyz_1, g_yy_yyyyzz_0, g_yy_yyyyzz_1, g_yy_yyyzzz_0, g_yy_yyyzzz_1, g_yy_yyzzzz_0, g_yy_yyzzzz_1, g_yy_yzzzzz_0, g_yy_yzzzzz_1, g_yy_zzzzzz_0, g_yy_zzzzzz_1, g_yyy_xxxxx_1, g_yyy_xxxxxx_1, g_yyy_xxxxxy_1, g_yyy_xxxxxz_1, g_yyy_xxxxy_1, g_yyy_xxxxyy_1, g_yyy_xxxxyz_1, g_yyy_xxxxz_1, g_yyy_xxxxzz_1, g_yyy_xxxyy_1, g_yyy_xxxyyy_1, g_yyy_xxxyyz_1, g_yyy_xxxyz_1, g_yyy_xxxyzz_1, g_yyy_xxxzz_1, g_yyy_xxxzzz_1, g_yyy_xxyyy_1, g_yyy_xxyyyy_1, g_yyy_xxyyyz_1, g_yyy_xxyyz_1, g_yyy_xxyyzz_1, g_yyy_xxyzz_1, g_yyy_xxyzzz_1, g_yyy_xxzzz_1, g_yyy_xxzzzz_1, g_yyy_xyyyy_1, g_yyy_xyyyyy_1, g_yyy_xyyyyz_1, g_yyy_xyyyz_1, g_yyy_xyyyzz_1, g_yyy_xyyzz_1, g_yyy_xyyzzz_1, g_yyy_xyzzz_1, g_yyy_xyzzzz_1, g_yyy_xzzzz_1, g_yyy_xzzzzz_1, g_yyy_yyyyy_1, g_yyy_yyyyyy_1, g_yyy_yyyyyz_1, g_yyy_yyyyz_1, g_yyy_yyyyzz_1, g_yyy_yyyzz_1, g_yyy_yyyzzz_1, g_yyy_yyzzz_1, g_yyy_yyzzzz_1, g_yyy_yzzzz_1, g_yyy_yzzzzz_1, g_yyy_zzzzz_1, g_yyy_zzzzzz_1, g_yyyy_xxxxxx_0, g_yyyy_xxxxxy_0, g_yyyy_xxxxxz_0, g_yyyy_xxxxyy_0, g_yyyy_xxxxyz_0, g_yyyy_xxxxzz_0, g_yyyy_xxxyyy_0, g_yyyy_xxxyyz_0, g_yyyy_xxxyzz_0, g_yyyy_xxxzzz_0, g_yyyy_xxyyyy_0, g_yyyy_xxyyyz_0, g_yyyy_xxyyzz_0, g_yyyy_xxyzzz_0, g_yyyy_xxzzzz_0, g_yyyy_xyyyyy_0, g_yyyy_xyyyyz_0, g_yyyy_xyyyzz_0, g_yyyy_xyyzzz_0, g_yyyy_xyzzzz_0, g_yyyy_xzzzzz_0, g_yyyy_yyyyyy_0, g_yyyy_yyyyyz_0, g_yyyy_yyyyzz_0, g_yyyy_yyyzzz_0, g_yyyy_yyzzzz_0, g_yyyy_yzzzzz_0, g_yyyy_zzzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyy_xxxxxx_0[i] = 3.0 * g_yy_xxxxxx_0[i] * fbe_0 - 3.0 * g_yy_xxxxxx_1[i] * fz_be_0 + g_yyy_xxxxxx_1[i] * pa_y[i];

        g_yyyy_xxxxxy_0[i] = 3.0 * g_yy_xxxxxy_0[i] * fbe_0 - 3.0 * g_yy_xxxxxy_1[i] * fz_be_0 + g_yyy_xxxxx_1[i] * fe_0 + g_yyy_xxxxxy_1[i] * pa_y[i];

        g_yyyy_xxxxxz_0[i] = 3.0 * g_yy_xxxxxz_0[i] * fbe_0 - 3.0 * g_yy_xxxxxz_1[i] * fz_be_0 + g_yyy_xxxxxz_1[i] * pa_y[i];

        g_yyyy_xxxxyy_0[i] = 3.0 * g_yy_xxxxyy_0[i] * fbe_0 - 3.0 * g_yy_xxxxyy_1[i] * fz_be_0 + 2.0 * g_yyy_xxxxy_1[i] * fe_0 + g_yyy_xxxxyy_1[i] * pa_y[i];

        g_yyyy_xxxxyz_0[i] = 3.0 * g_yy_xxxxyz_0[i] * fbe_0 - 3.0 * g_yy_xxxxyz_1[i] * fz_be_0 + g_yyy_xxxxz_1[i] * fe_0 + g_yyy_xxxxyz_1[i] * pa_y[i];

        g_yyyy_xxxxzz_0[i] = 3.0 * g_yy_xxxxzz_0[i] * fbe_0 - 3.0 * g_yy_xxxxzz_1[i] * fz_be_0 + g_yyy_xxxxzz_1[i] * pa_y[i];

        g_yyyy_xxxyyy_0[i] = 3.0 * g_yy_xxxyyy_0[i] * fbe_0 - 3.0 * g_yy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_yyy_xxxyy_1[i] * fe_0 + g_yyy_xxxyyy_1[i] * pa_y[i];

        g_yyyy_xxxyyz_0[i] = 3.0 * g_yy_xxxyyz_0[i] * fbe_0 - 3.0 * g_yy_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyy_xxxyz_1[i] * fe_0 + g_yyy_xxxyyz_1[i] * pa_y[i];

        g_yyyy_xxxyzz_0[i] = 3.0 * g_yy_xxxyzz_0[i] * fbe_0 - 3.0 * g_yy_xxxyzz_1[i] * fz_be_0 + g_yyy_xxxzz_1[i] * fe_0 + g_yyy_xxxyzz_1[i] * pa_y[i];

        g_yyyy_xxxzzz_0[i] = 3.0 * g_yy_xxxzzz_0[i] * fbe_0 - 3.0 * g_yy_xxxzzz_1[i] * fz_be_0 + g_yyy_xxxzzz_1[i] * pa_y[i];

        g_yyyy_xxyyyy_0[i] = 3.0 * g_yy_xxyyyy_0[i] * fbe_0 - 3.0 * g_yy_xxyyyy_1[i] * fz_be_0 + 4.0 * g_yyy_xxyyy_1[i] * fe_0 + g_yyy_xxyyyy_1[i] * pa_y[i];

        g_yyyy_xxyyyz_0[i] = 3.0 * g_yy_xxyyyz_0[i] * fbe_0 - 3.0 * g_yy_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyy_xxyyz_1[i] * fe_0 + g_yyy_xxyyyz_1[i] * pa_y[i];

        g_yyyy_xxyyzz_0[i] = 3.0 * g_yy_xxyyzz_0[i] * fbe_0 - 3.0 * g_yy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyy_xxyzz_1[i] * fe_0 + g_yyy_xxyyzz_1[i] * pa_y[i];

        g_yyyy_xxyzzz_0[i] = 3.0 * g_yy_xxyzzz_0[i] * fbe_0 - 3.0 * g_yy_xxyzzz_1[i] * fz_be_0 + g_yyy_xxzzz_1[i] * fe_0 + g_yyy_xxyzzz_1[i] * pa_y[i];

        g_yyyy_xxzzzz_0[i] = 3.0 * g_yy_xxzzzz_0[i] * fbe_0 - 3.0 * g_yy_xxzzzz_1[i] * fz_be_0 + g_yyy_xxzzzz_1[i] * pa_y[i];

        g_yyyy_xyyyyy_0[i] = 3.0 * g_yy_xyyyyy_0[i] * fbe_0 - 3.0 * g_yy_xyyyyy_1[i] * fz_be_0 + 5.0 * g_yyy_xyyyy_1[i] * fe_0 + g_yyy_xyyyyy_1[i] * pa_y[i];

        g_yyyy_xyyyyz_0[i] = 3.0 * g_yy_xyyyyz_0[i] * fbe_0 - 3.0 * g_yy_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyy_xyyyz_1[i] * fe_0 + g_yyy_xyyyyz_1[i] * pa_y[i];

        g_yyyy_xyyyzz_0[i] = 3.0 * g_yy_xyyyzz_0[i] * fbe_0 - 3.0 * g_yy_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyy_xyyzz_1[i] * fe_0 + g_yyy_xyyyzz_1[i] * pa_y[i];

        g_yyyy_xyyzzz_0[i] = 3.0 * g_yy_xyyzzz_0[i] * fbe_0 - 3.0 * g_yy_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyy_xyzzz_1[i] * fe_0 + g_yyy_xyyzzz_1[i] * pa_y[i];

        g_yyyy_xyzzzz_0[i] = 3.0 * g_yy_xyzzzz_0[i] * fbe_0 - 3.0 * g_yy_xyzzzz_1[i] * fz_be_0 + g_yyy_xzzzz_1[i] * fe_0 + g_yyy_xyzzzz_1[i] * pa_y[i];

        g_yyyy_xzzzzz_0[i] = 3.0 * g_yy_xzzzzz_0[i] * fbe_0 - 3.0 * g_yy_xzzzzz_1[i] * fz_be_0 + g_yyy_xzzzzz_1[i] * pa_y[i];

        g_yyyy_yyyyyy_0[i] = 3.0 * g_yy_yyyyyy_0[i] * fbe_0 - 3.0 * g_yy_yyyyyy_1[i] * fz_be_0 + 6.0 * g_yyy_yyyyy_1[i] * fe_0 + g_yyy_yyyyyy_1[i] * pa_y[i];

        g_yyyy_yyyyyz_0[i] = 3.0 * g_yy_yyyyyz_0[i] * fbe_0 - 3.0 * g_yy_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyy_yyyyz_1[i] * fe_0 + g_yyy_yyyyyz_1[i] * pa_y[i];

        g_yyyy_yyyyzz_0[i] = 3.0 * g_yy_yyyyzz_0[i] * fbe_0 - 3.0 * g_yy_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyy_yyyzz_1[i] * fe_0 + g_yyy_yyyyzz_1[i] * pa_y[i];

        g_yyyy_yyyzzz_0[i] = 3.0 * g_yy_yyyzzz_0[i] * fbe_0 - 3.0 * g_yy_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyy_yyzzz_1[i] * fe_0 + g_yyy_yyyzzz_1[i] * pa_y[i];

        g_yyyy_yyzzzz_0[i] = 3.0 * g_yy_yyzzzz_0[i] * fbe_0 - 3.0 * g_yy_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyy_yzzzz_1[i] * fe_0 + g_yyy_yyzzzz_1[i] * pa_y[i];

        g_yyyy_yzzzzz_0[i] = 3.0 * g_yy_yzzzzz_0[i] * fbe_0 - 3.0 * g_yy_yzzzzz_1[i] * fz_be_0 + g_yyy_zzzzz_1[i] * fe_0 + g_yyy_yzzzzz_1[i] * pa_y[i];

        g_yyyy_zzzzzz_0[i] = 3.0 * g_yy_zzzzzz_0[i] * fbe_0 - 3.0 * g_yy_zzzzzz_1[i] * fz_be_0 + g_yyy_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 308-336 components of targeted buffer : GI

    auto g_yyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 308);

    auto g_yyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 309);

    auto g_yyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 310);

    auto g_yyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 311);

    auto g_yyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 312);

    auto g_yyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 313);

    auto g_yyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 314);

    auto g_yyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 315);

    auto g_yyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 316);

    auto g_yyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 317);

    auto g_yyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 318);

    auto g_yyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 319);

    auto g_yyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 320);

    auto g_yyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 321);

    auto g_yyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 322);

    auto g_yyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 323);

    auto g_yyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 324);

    auto g_yyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 325);

    auto g_yyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 326);

    auto g_yyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 327);

    auto g_yyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 328);

    auto g_yyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 329);

    auto g_yyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 330);

    auto g_yyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 331);

    auto g_yyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 332);

    auto g_yyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 333);

    auto g_yyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 334);

    auto g_yyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 335);

    #pragma omp simd aligned(g_yyy_xxxxx_1, g_yyy_xxxxxx_1, g_yyy_xxxxxy_1, g_yyy_xxxxxz_1, g_yyy_xxxxy_1, g_yyy_xxxxyy_1, g_yyy_xxxxyz_1, g_yyy_xxxxz_1, g_yyy_xxxxzz_1, g_yyy_xxxyy_1, g_yyy_xxxyyy_1, g_yyy_xxxyyz_1, g_yyy_xxxyz_1, g_yyy_xxxyzz_1, g_yyy_xxxzz_1, g_yyy_xxxzzz_1, g_yyy_xxyyy_1, g_yyy_xxyyyy_1, g_yyy_xxyyyz_1, g_yyy_xxyyz_1, g_yyy_xxyyzz_1, g_yyy_xxyzz_1, g_yyy_xxyzzz_1, g_yyy_xxzzz_1, g_yyy_xxzzzz_1, g_yyy_xyyyy_1, g_yyy_xyyyyy_1, g_yyy_xyyyyz_1, g_yyy_xyyyz_1, g_yyy_xyyyzz_1, g_yyy_xyyzz_1, g_yyy_xyyzzz_1, g_yyy_xyzzz_1, g_yyy_xyzzzz_1, g_yyy_xzzzz_1, g_yyy_xzzzzz_1, g_yyy_yyyyy_1, g_yyy_yyyyyy_1, g_yyy_yyyyyz_1, g_yyy_yyyyz_1, g_yyy_yyyyzz_1, g_yyy_yyyzz_1, g_yyy_yyyzzz_1, g_yyy_yyzzz_1, g_yyy_yyzzzz_1, g_yyy_yzzzz_1, g_yyy_yzzzzz_1, g_yyy_zzzzz_1, g_yyy_zzzzzz_1, g_yyyz_xxxxxx_0, g_yyyz_xxxxxy_0, g_yyyz_xxxxxz_0, g_yyyz_xxxxyy_0, g_yyyz_xxxxyz_0, g_yyyz_xxxxzz_0, g_yyyz_xxxyyy_0, g_yyyz_xxxyyz_0, g_yyyz_xxxyzz_0, g_yyyz_xxxzzz_0, g_yyyz_xxyyyy_0, g_yyyz_xxyyyz_0, g_yyyz_xxyyzz_0, g_yyyz_xxyzzz_0, g_yyyz_xxzzzz_0, g_yyyz_xyyyyy_0, g_yyyz_xyyyyz_0, g_yyyz_xyyyzz_0, g_yyyz_xyyzzz_0, g_yyyz_xyzzzz_0, g_yyyz_xzzzzz_0, g_yyyz_yyyyyy_0, g_yyyz_yyyyyz_0, g_yyyz_yyyyzz_0, g_yyyz_yyyzzz_0, g_yyyz_yyzzzz_0, g_yyyz_yzzzzz_0, g_yyyz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyz_xxxxxx_0[i] = g_yyy_xxxxxx_1[i] * pa_z[i];

        g_yyyz_xxxxxy_0[i] = g_yyy_xxxxxy_1[i] * pa_z[i];

        g_yyyz_xxxxxz_0[i] = g_yyy_xxxxx_1[i] * fe_0 + g_yyy_xxxxxz_1[i] * pa_z[i];

        g_yyyz_xxxxyy_0[i] = g_yyy_xxxxyy_1[i] * pa_z[i];

        g_yyyz_xxxxyz_0[i] = g_yyy_xxxxy_1[i] * fe_0 + g_yyy_xxxxyz_1[i] * pa_z[i];

        g_yyyz_xxxxzz_0[i] = 2.0 * g_yyy_xxxxz_1[i] * fe_0 + g_yyy_xxxxzz_1[i] * pa_z[i];

        g_yyyz_xxxyyy_0[i] = g_yyy_xxxyyy_1[i] * pa_z[i];

        g_yyyz_xxxyyz_0[i] = g_yyy_xxxyy_1[i] * fe_0 + g_yyy_xxxyyz_1[i] * pa_z[i];

        g_yyyz_xxxyzz_0[i] = 2.0 * g_yyy_xxxyz_1[i] * fe_0 + g_yyy_xxxyzz_1[i] * pa_z[i];

        g_yyyz_xxxzzz_0[i] = 3.0 * g_yyy_xxxzz_1[i] * fe_0 + g_yyy_xxxzzz_1[i] * pa_z[i];

        g_yyyz_xxyyyy_0[i] = g_yyy_xxyyyy_1[i] * pa_z[i];

        g_yyyz_xxyyyz_0[i] = g_yyy_xxyyy_1[i] * fe_0 + g_yyy_xxyyyz_1[i] * pa_z[i];

        g_yyyz_xxyyzz_0[i] = 2.0 * g_yyy_xxyyz_1[i] * fe_0 + g_yyy_xxyyzz_1[i] * pa_z[i];

        g_yyyz_xxyzzz_0[i] = 3.0 * g_yyy_xxyzz_1[i] * fe_0 + g_yyy_xxyzzz_1[i] * pa_z[i];

        g_yyyz_xxzzzz_0[i] = 4.0 * g_yyy_xxzzz_1[i] * fe_0 + g_yyy_xxzzzz_1[i] * pa_z[i];

        g_yyyz_xyyyyy_0[i] = g_yyy_xyyyyy_1[i] * pa_z[i];

        g_yyyz_xyyyyz_0[i] = g_yyy_xyyyy_1[i] * fe_0 + g_yyy_xyyyyz_1[i] * pa_z[i];

        g_yyyz_xyyyzz_0[i] = 2.0 * g_yyy_xyyyz_1[i] * fe_0 + g_yyy_xyyyzz_1[i] * pa_z[i];

        g_yyyz_xyyzzz_0[i] = 3.0 * g_yyy_xyyzz_1[i] * fe_0 + g_yyy_xyyzzz_1[i] * pa_z[i];

        g_yyyz_xyzzzz_0[i] = 4.0 * g_yyy_xyzzz_1[i] * fe_0 + g_yyy_xyzzzz_1[i] * pa_z[i];

        g_yyyz_xzzzzz_0[i] = 5.0 * g_yyy_xzzzz_1[i] * fe_0 + g_yyy_xzzzzz_1[i] * pa_z[i];

        g_yyyz_yyyyyy_0[i] = g_yyy_yyyyyy_1[i] * pa_z[i];

        g_yyyz_yyyyyz_0[i] = g_yyy_yyyyy_1[i] * fe_0 + g_yyy_yyyyyz_1[i] * pa_z[i];

        g_yyyz_yyyyzz_0[i] = 2.0 * g_yyy_yyyyz_1[i] * fe_0 + g_yyy_yyyyzz_1[i] * pa_z[i];

        g_yyyz_yyyzzz_0[i] = 3.0 * g_yyy_yyyzz_1[i] * fe_0 + g_yyy_yyyzzz_1[i] * pa_z[i];

        g_yyyz_yyzzzz_0[i] = 4.0 * g_yyy_yyzzz_1[i] * fe_0 + g_yyy_yyzzzz_1[i] * pa_z[i];

        g_yyyz_yzzzzz_0[i] = 5.0 * g_yyy_yzzzz_1[i] * fe_0 + g_yyy_yzzzzz_1[i] * pa_z[i];

        g_yyyz_zzzzzz_0[i] = 6.0 * g_yyy_zzzzz_1[i] * fe_0 + g_yyy_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 336-364 components of targeted buffer : GI

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

    #pragma omp simd aligned(g_yy_xxxxxy_0, g_yy_xxxxxy_1, g_yy_xxxxyy_0, g_yy_xxxxyy_1, g_yy_xxxyyy_0, g_yy_xxxyyy_1, g_yy_xxyyyy_0, g_yy_xxyyyy_1, g_yy_xyyyyy_0, g_yy_xyyyyy_1, g_yy_yyyyyy_0, g_yy_yyyyyy_1, g_yyz_xxxxxy_1, g_yyz_xxxxyy_1, g_yyz_xxxyyy_1, g_yyz_xxyyyy_1, g_yyz_xyyyyy_1, g_yyz_yyyyyy_1, g_yyzz_xxxxxx_0, g_yyzz_xxxxxy_0, g_yyzz_xxxxxz_0, g_yyzz_xxxxyy_0, g_yyzz_xxxxyz_0, g_yyzz_xxxxzz_0, g_yyzz_xxxyyy_0, g_yyzz_xxxyyz_0, g_yyzz_xxxyzz_0, g_yyzz_xxxzzz_0, g_yyzz_xxyyyy_0, g_yyzz_xxyyyz_0, g_yyzz_xxyyzz_0, g_yyzz_xxyzzz_0, g_yyzz_xxzzzz_0, g_yyzz_xyyyyy_0, g_yyzz_xyyyyz_0, g_yyzz_xyyyzz_0, g_yyzz_xyyzzz_0, g_yyzz_xyzzzz_0, g_yyzz_xzzzzz_0, g_yyzz_yyyyyy_0, g_yyzz_yyyyyz_0, g_yyzz_yyyyzz_0, g_yyzz_yyyzzz_0, g_yyzz_yyzzzz_0, g_yyzz_yzzzzz_0, g_yyzz_zzzzzz_0, g_yzz_xxxxxx_1, g_yzz_xxxxxz_1, g_yzz_xxxxyz_1, g_yzz_xxxxz_1, g_yzz_xxxxzz_1, g_yzz_xxxyyz_1, g_yzz_xxxyz_1, g_yzz_xxxyzz_1, g_yzz_xxxzz_1, g_yzz_xxxzzz_1, g_yzz_xxyyyz_1, g_yzz_xxyyz_1, g_yzz_xxyyzz_1, g_yzz_xxyzz_1, g_yzz_xxyzzz_1, g_yzz_xxzzz_1, g_yzz_xxzzzz_1, g_yzz_xyyyyz_1, g_yzz_xyyyz_1, g_yzz_xyyyzz_1, g_yzz_xyyzz_1, g_yzz_xyyzzz_1, g_yzz_xyzzz_1, g_yzz_xyzzzz_1, g_yzz_xzzzz_1, g_yzz_xzzzzz_1, g_yzz_yyyyyz_1, g_yzz_yyyyz_1, g_yzz_yyyyzz_1, g_yzz_yyyzz_1, g_yzz_yyyzzz_1, g_yzz_yyzzz_1, g_yzz_yyzzzz_1, g_yzz_yzzzz_1, g_yzz_yzzzzz_1, g_yzz_zzzzz_1, g_yzz_zzzzzz_1, g_zz_xxxxxx_0, g_zz_xxxxxx_1, g_zz_xxxxxz_0, g_zz_xxxxxz_1, g_zz_xxxxyz_0, g_zz_xxxxyz_1, g_zz_xxxxzz_0, g_zz_xxxxzz_1, g_zz_xxxyyz_0, g_zz_xxxyyz_1, g_zz_xxxyzz_0, g_zz_xxxyzz_1, g_zz_xxxzzz_0, g_zz_xxxzzz_1, g_zz_xxyyyz_0, g_zz_xxyyyz_1, g_zz_xxyyzz_0, g_zz_xxyyzz_1, g_zz_xxyzzz_0, g_zz_xxyzzz_1, g_zz_xxzzzz_0, g_zz_xxzzzz_1, g_zz_xyyyyz_0, g_zz_xyyyyz_1, g_zz_xyyyzz_0, g_zz_xyyyzz_1, g_zz_xyyzzz_0, g_zz_xyyzzz_1, g_zz_xyzzzz_0, g_zz_xyzzzz_1, g_zz_xzzzzz_0, g_zz_xzzzzz_1, g_zz_yyyyyz_0, g_zz_yyyyyz_1, g_zz_yyyyzz_0, g_zz_yyyyzz_1, g_zz_yyyzzz_0, g_zz_yyyzzz_1, g_zz_yyzzzz_0, g_zz_yyzzzz_1, g_zz_yzzzzz_0, g_zz_yzzzzz_1, g_zz_zzzzzz_0, g_zz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzz_xxxxxx_0[i] = g_zz_xxxxxx_0[i] * fbe_0 - g_zz_xxxxxx_1[i] * fz_be_0 + g_yzz_xxxxxx_1[i] * pa_y[i];

        g_yyzz_xxxxxy_0[i] = g_yy_xxxxxy_0[i] * fbe_0 - g_yy_xxxxxy_1[i] * fz_be_0 + g_yyz_xxxxxy_1[i] * pa_z[i];

        g_yyzz_xxxxxz_0[i] = g_zz_xxxxxz_0[i] * fbe_0 - g_zz_xxxxxz_1[i] * fz_be_0 + g_yzz_xxxxxz_1[i] * pa_y[i];

        g_yyzz_xxxxyy_0[i] = g_yy_xxxxyy_0[i] * fbe_0 - g_yy_xxxxyy_1[i] * fz_be_0 + g_yyz_xxxxyy_1[i] * pa_z[i];

        g_yyzz_xxxxyz_0[i] = g_zz_xxxxyz_0[i] * fbe_0 - g_zz_xxxxyz_1[i] * fz_be_0 + g_yzz_xxxxz_1[i] * fe_0 + g_yzz_xxxxyz_1[i] * pa_y[i];

        g_yyzz_xxxxzz_0[i] = g_zz_xxxxzz_0[i] * fbe_0 - g_zz_xxxxzz_1[i] * fz_be_0 + g_yzz_xxxxzz_1[i] * pa_y[i];

        g_yyzz_xxxyyy_0[i] = g_yy_xxxyyy_0[i] * fbe_0 - g_yy_xxxyyy_1[i] * fz_be_0 + g_yyz_xxxyyy_1[i] * pa_z[i];

        g_yyzz_xxxyyz_0[i] = g_zz_xxxyyz_0[i] * fbe_0 - g_zz_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yzz_xxxyz_1[i] * fe_0 + g_yzz_xxxyyz_1[i] * pa_y[i];

        g_yyzz_xxxyzz_0[i] = g_zz_xxxyzz_0[i] * fbe_0 - g_zz_xxxyzz_1[i] * fz_be_0 + g_yzz_xxxzz_1[i] * fe_0 + g_yzz_xxxyzz_1[i] * pa_y[i];

        g_yyzz_xxxzzz_0[i] = g_zz_xxxzzz_0[i] * fbe_0 - g_zz_xxxzzz_1[i] * fz_be_0 + g_yzz_xxxzzz_1[i] * pa_y[i];

        g_yyzz_xxyyyy_0[i] = g_yy_xxyyyy_0[i] * fbe_0 - g_yy_xxyyyy_1[i] * fz_be_0 + g_yyz_xxyyyy_1[i] * pa_z[i];

        g_yyzz_xxyyyz_0[i] = g_zz_xxyyyz_0[i] * fbe_0 - g_zz_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yzz_xxyyz_1[i] * fe_0 + g_yzz_xxyyyz_1[i] * pa_y[i];

        g_yyzz_xxyyzz_0[i] = g_zz_xxyyzz_0[i] * fbe_0 - g_zz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yzz_xxyzz_1[i] * fe_0 + g_yzz_xxyyzz_1[i] * pa_y[i];

        g_yyzz_xxyzzz_0[i] = g_zz_xxyzzz_0[i] * fbe_0 - g_zz_xxyzzz_1[i] * fz_be_0 + g_yzz_xxzzz_1[i] * fe_0 + g_yzz_xxyzzz_1[i] * pa_y[i];

        g_yyzz_xxzzzz_0[i] = g_zz_xxzzzz_0[i] * fbe_0 - g_zz_xxzzzz_1[i] * fz_be_0 + g_yzz_xxzzzz_1[i] * pa_y[i];

        g_yyzz_xyyyyy_0[i] = g_yy_xyyyyy_0[i] * fbe_0 - g_yy_xyyyyy_1[i] * fz_be_0 + g_yyz_xyyyyy_1[i] * pa_z[i];

        g_yyzz_xyyyyz_0[i] = g_zz_xyyyyz_0[i] * fbe_0 - g_zz_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yzz_xyyyz_1[i] * fe_0 + g_yzz_xyyyyz_1[i] * pa_y[i];

        g_yyzz_xyyyzz_0[i] = g_zz_xyyyzz_0[i] * fbe_0 - g_zz_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yzz_xyyzz_1[i] * fe_0 + g_yzz_xyyyzz_1[i] * pa_y[i];

        g_yyzz_xyyzzz_0[i] = g_zz_xyyzzz_0[i] * fbe_0 - g_zz_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yzz_xyzzz_1[i] * fe_0 + g_yzz_xyyzzz_1[i] * pa_y[i];

        g_yyzz_xyzzzz_0[i] = g_zz_xyzzzz_0[i] * fbe_0 - g_zz_xyzzzz_1[i] * fz_be_0 + g_yzz_xzzzz_1[i] * fe_0 + g_yzz_xyzzzz_1[i] * pa_y[i];

        g_yyzz_xzzzzz_0[i] = g_zz_xzzzzz_0[i] * fbe_0 - g_zz_xzzzzz_1[i] * fz_be_0 + g_yzz_xzzzzz_1[i] * pa_y[i];

        g_yyzz_yyyyyy_0[i] = g_yy_yyyyyy_0[i] * fbe_0 - g_yy_yyyyyy_1[i] * fz_be_0 + g_yyz_yyyyyy_1[i] * pa_z[i];

        g_yyzz_yyyyyz_0[i] = g_zz_yyyyyz_0[i] * fbe_0 - g_zz_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yzz_yyyyz_1[i] * fe_0 + g_yzz_yyyyyz_1[i] * pa_y[i];

        g_yyzz_yyyyzz_0[i] = g_zz_yyyyzz_0[i] * fbe_0 - g_zz_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yzz_yyyzz_1[i] * fe_0 + g_yzz_yyyyzz_1[i] * pa_y[i];

        g_yyzz_yyyzzz_0[i] = g_zz_yyyzzz_0[i] * fbe_0 - g_zz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yzz_yyzzz_1[i] * fe_0 + g_yzz_yyyzzz_1[i] * pa_y[i];

        g_yyzz_yyzzzz_0[i] = g_zz_yyzzzz_0[i] * fbe_0 - g_zz_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yzz_yzzzz_1[i] * fe_0 + g_yzz_yyzzzz_1[i] * pa_y[i];

        g_yyzz_yzzzzz_0[i] = g_zz_yzzzzz_0[i] * fbe_0 - g_zz_yzzzzz_1[i] * fz_be_0 + g_yzz_zzzzz_1[i] * fe_0 + g_yzz_yzzzzz_1[i] * pa_y[i];

        g_yyzz_zzzzzz_0[i] = g_zz_zzzzzz_0[i] * fbe_0 - g_zz_zzzzzz_1[i] * fz_be_0 + g_yzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 364-392 components of targeted buffer : GI

    auto g_yzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_gi + 364);

    auto g_yzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_gi + 365);

    auto g_yzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_gi + 366);

    auto g_yzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_gi + 367);

    auto g_yzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_gi + 368);

    auto g_yzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_gi + 369);

    auto g_yzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_gi + 370);

    auto g_yzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_gi + 371);

    auto g_yzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_gi + 372);

    auto g_yzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_gi + 373);

    auto g_yzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_gi + 374);

    auto g_yzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_gi + 375);

    auto g_yzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_gi + 376);

    auto g_yzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_gi + 377);

    auto g_yzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_gi + 378);

    auto g_yzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_gi + 379);

    auto g_yzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_gi + 380);

    auto g_yzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_gi + 381);

    auto g_yzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_gi + 382);

    auto g_yzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_gi + 383);

    auto g_yzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_gi + 384);

    auto g_yzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_gi + 385);

    auto g_yzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_gi + 386);

    auto g_yzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_gi + 387);

    auto g_yzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_gi + 388);

    auto g_yzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_gi + 389);

    auto g_yzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_gi + 390);

    auto g_yzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_gi + 391);

    #pragma omp simd aligned(g_yzzz_xxxxxx_0, g_yzzz_xxxxxy_0, g_yzzz_xxxxxz_0, g_yzzz_xxxxyy_0, g_yzzz_xxxxyz_0, g_yzzz_xxxxzz_0, g_yzzz_xxxyyy_0, g_yzzz_xxxyyz_0, g_yzzz_xxxyzz_0, g_yzzz_xxxzzz_0, g_yzzz_xxyyyy_0, g_yzzz_xxyyyz_0, g_yzzz_xxyyzz_0, g_yzzz_xxyzzz_0, g_yzzz_xxzzzz_0, g_yzzz_xyyyyy_0, g_yzzz_xyyyyz_0, g_yzzz_xyyyzz_0, g_yzzz_xyyzzz_0, g_yzzz_xyzzzz_0, g_yzzz_xzzzzz_0, g_yzzz_yyyyyy_0, g_yzzz_yyyyyz_0, g_yzzz_yyyyzz_0, g_yzzz_yyyzzz_0, g_yzzz_yyzzzz_0, g_yzzz_yzzzzz_0, g_yzzz_zzzzzz_0, g_zzz_xxxxx_1, g_zzz_xxxxxx_1, g_zzz_xxxxxy_1, g_zzz_xxxxxz_1, g_zzz_xxxxy_1, g_zzz_xxxxyy_1, g_zzz_xxxxyz_1, g_zzz_xxxxz_1, g_zzz_xxxxzz_1, g_zzz_xxxyy_1, g_zzz_xxxyyy_1, g_zzz_xxxyyz_1, g_zzz_xxxyz_1, g_zzz_xxxyzz_1, g_zzz_xxxzz_1, g_zzz_xxxzzz_1, g_zzz_xxyyy_1, g_zzz_xxyyyy_1, g_zzz_xxyyyz_1, g_zzz_xxyyz_1, g_zzz_xxyyzz_1, g_zzz_xxyzz_1, g_zzz_xxyzzz_1, g_zzz_xxzzz_1, g_zzz_xxzzzz_1, g_zzz_xyyyy_1, g_zzz_xyyyyy_1, g_zzz_xyyyyz_1, g_zzz_xyyyz_1, g_zzz_xyyyzz_1, g_zzz_xyyzz_1, g_zzz_xyyzzz_1, g_zzz_xyzzz_1, g_zzz_xyzzzz_1, g_zzz_xzzzz_1, g_zzz_xzzzzz_1, g_zzz_yyyyy_1, g_zzz_yyyyyy_1, g_zzz_yyyyyz_1, g_zzz_yyyyz_1, g_zzz_yyyyzz_1, g_zzz_yyyzz_1, g_zzz_yyyzzz_1, g_zzz_yyzzz_1, g_zzz_yyzzzz_1, g_zzz_yzzzz_1, g_zzz_yzzzzz_1, g_zzz_zzzzz_1, g_zzz_zzzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzz_xxxxxx_0[i] = g_zzz_xxxxxx_1[i] * pa_y[i];

        g_yzzz_xxxxxy_0[i] = g_zzz_xxxxx_1[i] * fe_0 + g_zzz_xxxxxy_1[i] * pa_y[i];

        g_yzzz_xxxxxz_0[i] = g_zzz_xxxxxz_1[i] * pa_y[i];

        g_yzzz_xxxxyy_0[i] = 2.0 * g_zzz_xxxxy_1[i] * fe_0 + g_zzz_xxxxyy_1[i] * pa_y[i];

        g_yzzz_xxxxyz_0[i] = g_zzz_xxxxz_1[i] * fe_0 + g_zzz_xxxxyz_1[i] * pa_y[i];

        g_yzzz_xxxxzz_0[i] = g_zzz_xxxxzz_1[i] * pa_y[i];

        g_yzzz_xxxyyy_0[i] = 3.0 * g_zzz_xxxyy_1[i] * fe_0 + g_zzz_xxxyyy_1[i] * pa_y[i];

        g_yzzz_xxxyyz_0[i] = 2.0 * g_zzz_xxxyz_1[i] * fe_0 + g_zzz_xxxyyz_1[i] * pa_y[i];

        g_yzzz_xxxyzz_0[i] = g_zzz_xxxzz_1[i] * fe_0 + g_zzz_xxxyzz_1[i] * pa_y[i];

        g_yzzz_xxxzzz_0[i] = g_zzz_xxxzzz_1[i] * pa_y[i];

        g_yzzz_xxyyyy_0[i] = 4.0 * g_zzz_xxyyy_1[i] * fe_0 + g_zzz_xxyyyy_1[i] * pa_y[i];

        g_yzzz_xxyyyz_0[i] = 3.0 * g_zzz_xxyyz_1[i] * fe_0 + g_zzz_xxyyyz_1[i] * pa_y[i];

        g_yzzz_xxyyzz_0[i] = 2.0 * g_zzz_xxyzz_1[i] * fe_0 + g_zzz_xxyyzz_1[i] * pa_y[i];

        g_yzzz_xxyzzz_0[i] = g_zzz_xxzzz_1[i] * fe_0 + g_zzz_xxyzzz_1[i] * pa_y[i];

        g_yzzz_xxzzzz_0[i] = g_zzz_xxzzzz_1[i] * pa_y[i];

        g_yzzz_xyyyyy_0[i] = 5.0 * g_zzz_xyyyy_1[i] * fe_0 + g_zzz_xyyyyy_1[i] * pa_y[i];

        g_yzzz_xyyyyz_0[i] = 4.0 * g_zzz_xyyyz_1[i] * fe_0 + g_zzz_xyyyyz_1[i] * pa_y[i];

        g_yzzz_xyyyzz_0[i] = 3.0 * g_zzz_xyyzz_1[i] * fe_0 + g_zzz_xyyyzz_1[i] * pa_y[i];

        g_yzzz_xyyzzz_0[i] = 2.0 * g_zzz_xyzzz_1[i] * fe_0 + g_zzz_xyyzzz_1[i] * pa_y[i];

        g_yzzz_xyzzzz_0[i] = g_zzz_xzzzz_1[i] * fe_0 + g_zzz_xyzzzz_1[i] * pa_y[i];

        g_yzzz_xzzzzz_0[i] = g_zzz_xzzzzz_1[i] * pa_y[i];

        g_yzzz_yyyyyy_0[i] = 6.0 * g_zzz_yyyyy_1[i] * fe_0 + g_zzz_yyyyyy_1[i] * pa_y[i];

        g_yzzz_yyyyyz_0[i] = 5.0 * g_zzz_yyyyz_1[i] * fe_0 + g_zzz_yyyyyz_1[i] * pa_y[i];

        g_yzzz_yyyyzz_0[i] = 4.0 * g_zzz_yyyzz_1[i] * fe_0 + g_zzz_yyyyzz_1[i] * pa_y[i];

        g_yzzz_yyyzzz_0[i] = 3.0 * g_zzz_yyzzz_1[i] * fe_0 + g_zzz_yyyzzz_1[i] * pa_y[i];

        g_yzzz_yyzzzz_0[i] = 2.0 * g_zzz_yzzzz_1[i] * fe_0 + g_zzz_yyzzzz_1[i] * pa_y[i];

        g_yzzz_yzzzzz_0[i] = g_zzz_zzzzz_1[i] * fe_0 + g_zzz_yzzzzz_1[i] * pa_y[i];

        g_yzzz_zzzzzz_0[i] = g_zzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 392-420 components of targeted buffer : GI

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

    #pragma omp simd aligned(g_zz_xxxxxx_0, g_zz_xxxxxx_1, g_zz_xxxxxy_0, g_zz_xxxxxy_1, g_zz_xxxxxz_0, g_zz_xxxxxz_1, g_zz_xxxxyy_0, g_zz_xxxxyy_1, g_zz_xxxxyz_0, g_zz_xxxxyz_1, g_zz_xxxxzz_0, g_zz_xxxxzz_1, g_zz_xxxyyy_0, g_zz_xxxyyy_1, g_zz_xxxyyz_0, g_zz_xxxyyz_1, g_zz_xxxyzz_0, g_zz_xxxyzz_1, g_zz_xxxzzz_0, g_zz_xxxzzz_1, g_zz_xxyyyy_0, g_zz_xxyyyy_1, g_zz_xxyyyz_0, g_zz_xxyyyz_1, g_zz_xxyyzz_0, g_zz_xxyyzz_1, g_zz_xxyzzz_0, g_zz_xxyzzz_1, g_zz_xxzzzz_0, g_zz_xxzzzz_1, g_zz_xyyyyy_0, g_zz_xyyyyy_1, g_zz_xyyyyz_0, g_zz_xyyyyz_1, g_zz_xyyyzz_0, g_zz_xyyyzz_1, g_zz_xyyzzz_0, g_zz_xyyzzz_1, g_zz_xyzzzz_0, g_zz_xyzzzz_1, g_zz_xzzzzz_0, g_zz_xzzzzz_1, g_zz_yyyyyy_0, g_zz_yyyyyy_1, g_zz_yyyyyz_0, g_zz_yyyyyz_1, g_zz_yyyyzz_0, g_zz_yyyyzz_1, g_zz_yyyzzz_0, g_zz_yyyzzz_1, g_zz_yyzzzz_0, g_zz_yyzzzz_1, g_zz_yzzzzz_0, g_zz_yzzzzz_1, g_zz_zzzzzz_0, g_zz_zzzzzz_1, g_zzz_xxxxx_1, g_zzz_xxxxxx_1, g_zzz_xxxxxy_1, g_zzz_xxxxxz_1, g_zzz_xxxxy_1, g_zzz_xxxxyy_1, g_zzz_xxxxyz_1, g_zzz_xxxxz_1, g_zzz_xxxxzz_1, g_zzz_xxxyy_1, g_zzz_xxxyyy_1, g_zzz_xxxyyz_1, g_zzz_xxxyz_1, g_zzz_xxxyzz_1, g_zzz_xxxzz_1, g_zzz_xxxzzz_1, g_zzz_xxyyy_1, g_zzz_xxyyyy_1, g_zzz_xxyyyz_1, g_zzz_xxyyz_1, g_zzz_xxyyzz_1, g_zzz_xxyzz_1, g_zzz_xxyzzz_1, g_zzz_xxzzz_1, g_zzz_xxzzzz_1, g_zzz_xyyyy_1, g_zzz_xyyyyy_1, g_zzz_xyyyyz_1, g_zzz_xyyyz_1, g_zzz_xyyyzz_1, g_zzz_xyyzz_1, g_zzz_xyyzzz_1, g_zzz_xyzzz_1, g_zzz_xyzzzz_1, g_zzz_xzzzz_1, g_zzz_xzzzzz_1, g_zzz_yyyyy_1, g_zzz_yyyyyy_1, g_zzz_yyyyyz_1, g_zzz_yyyyz_1, g_zzz_yyyyzz_1, g_zzz_yyyzz_1, g_zzz_yyyzzz_1, g_zzz_yyzzz_1, g_zzz_yyzzzz_1, g_zzz_yzzzz_1, g_zzz_yzzzzz_1, g_zzz_zzzzz_1, g_zzz_zzzzzz_1, g_zzzz_xxxxxx_0, g_zzzz_xxxxxy_0, g_zzzz_xxxxxz_0, g_zzzz_xxxxyy_0, g_zzzz_xxxxyz_0, g_zzzz_xxxxzz_0, g_zzzz_xxxyyy_0, g_zzzz_xxxyyz_0, g_zzzz_xxxyzz_0, g_zzzz_xxxzzz_0, g_zzzz_xxyyyy_0, g_zzzz_xxyyyz_0, g_zzzz_xxyyzz_0, g_zzzz_xxyzzz_0, g_zzzz_xxzzzz_0, g_zzzz_xyyyyy_0, g_zzzz_xyyyyz_0, g_zzzz_xyyyzz_0, g_zzzz_xyyzzz_0, g_zzzz_xyzzzz_0, g_zzzz_xzzzzz_0, g_zzzz_yyyyyy_0, g_zzzz_yyyyyz_0, g_zzzz_yyyyzz_0, g_zzzz_yyyzzz_0, g_zzzz_yyzzzz_0, g_zzzz_yzzzzz_0, g_zzzz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzz_xxxxxx_0[i] = 3.0 * g_zz_xxxxxx_0[i] * fbe_0 - 3.0 * g_zz_xxxxxx_1[i] * fz_be_0 + g_zzz_xxxxxx_1[i] * pa_z[i];

        g_zzzz_xxxxxy_0[i] = 3.0 * g_zz_xxxxxy_0[i] * fbe_0 - 3.0 * g_zz_xxxxxy_1[i] * fz_be_0 + g_zzz_xxxxxy_1[i] * pa_z[i];

        g_zzzz_xxxxxz_0[i] = 3.0 * g_zz_xxxxxz_0[i] * fbe_0 - 3.0 * g_zz_xxxxxz_1[i] * fz_be_0 + g_zzz_xxxxx_1[i] * fe_0 + g_zzz_xxxxxz_1[i] * pa_z[i];

        g_zzzz_xxxxyy_0[i] = 3.0 * g_zz_xxxxyy_0[i] * fbe_0 - 3.0 * g_zz_xxxxyy_1[i] * fz_be_0 + g_zzz_xxxxyy_1[i] * pa_z[i];

        g_zzzz_xxxxyz_0[i] = 3.0 * g_zz_xxxxyz_0[i] * fbe_0 - 3.0 * g_zz_xxxxyz_1[i] * fz_be_0 + g_zzz_xxxxy_1[i] * fe_0 + g_zzz_xxxxyz_1[i] * pa_z[i];

        g_zzzz_xxxxzz_0[i] = 3.0 * g_zz_xxxxzz_0[i] * fbe_0 - 3.0 * g_zz_xxxxzz_1[i] * fz_be_0 + 2.0 * g_zzz_xxxxz_1[i] * fe_0 + g_zzz_xxxxzz_1[i] * pa_z[i];

        g_zzzz_xxxyyy_0[i] = 3.0 * g_zz_xxxyyy_0[i] * fbe_0 - 3.0 * g_zz_xxxyyy_1[i] * fz_be_0 + g_zzz_xxxyyy_1[i] * pa_z[i];

        g_zzzz_xxxyyz_0[i] = 3.0 * g_zz_xxxyyz_0[i] * fbe_0 - 3.0 * g_zz_xxxyyz_1[i] * fz_be_0 + g_zzz_xxxyy_1[i] * fe_0 + g_zzz_xxxyyz_1[i] * pa_z[i];

        g_zzzz_xxxyzz_0[i] = 3.0 * g_zz_xxxyzz_0[i] * fbe_0 - 3.0 * g_zz_xxxyzz_1[i] * fz_be_0 + 2.0 * g_zzz_xxxyz_1[i] * fe_0 + g_zzz_xxxyzz_1[i] * pa_z[i];

        g_zzzz_xxxzzz_0[i] = 3.0 * g_zz_xxxzzz_0[i] * fbe_0 - 3.0 * g_zz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_zzz_xxxzz_1[i] * fe_0 + g_zzz_xxxzzz_1[i] * pa_z[i];

        g_zzzz_xxyyyy_0[i] = 3.0 * g_zz_xxyyyy_0[i] * fbe_0 - 3.0 * g_zz_xxyyyy_1[i] * fz_be_0 + g_zzz_xxyyyy_1[i] * pa_z[i];

        g_zzzz_xxyyyz_0[i] = 3.0 * g_zz_xxyyyz_0[i] * fbe_0 - 3.0 * g_zz_xxyyyz_1[i] * fz_be_0 + g_zzz_xxyyy_1[i] * fe_0 + g_zzz_xxyyyz_1[i] * pa_z[i];

        g_zzzz_xxyyzz_0[i] = 3.0 * g_zz_xxyyzz_0[i] * fbe_0 - 3.0 * g_zz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_xxyyz_1[i] * fe_0 + g_zzz_xxyyzz_1[i] * pa_z[i];

        g_zzzz_xxyzzz_0[i] = 3.0 * g_zz_xxyzzz_0[i] * fbe_0 - 3.0 * g_zz_xxyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_xxyzz_1[i] * fe_0 + g_zzz_xxyzzz_1[i] * pa_z[i];

        g_zzzz_xxzzzz_0[i] = 3.0 * g_zz_xxzzzz_0[i] * fbe_0 - 3.0 * g_zz_xxzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_xxzzz_1[i] * fe_0 + g_zzz_xxzzzz_1[i] * pa_z[i];

        g_zzzz_xyyyyy_0[i] = 3.0 * g_zz_xyyyyy_0[i] * fbe_0 - 3.0 * g_zz_xyyyyy_1[i] * fz_be_0 + g_zzz_xyyyyy_1[i] * pa_z[i];

        g_zzzz_xyyyyz_0[i] = 3.0 * g_zz_xyyyyz_0[i] * fbe_0 - 3.0 * g_zz_xyyyyz_1[i] * fz_be_0 + g_zzz_xyyyy_1[i] * fe_0 + g_zzz_xyyyyz_1[i] * pa_z[i];

        g_zzzz_xyyyzz_0[i] = 3.0 * g_zz_xyyyzz_0[i] * fbe_0 - 3.0 * g_zz_xyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_xyyyz_1[i] * fe_0 + g_zzz_xyyyzz_1[i] * pa_z[i];

        g_zzzz_xyyzzz_0[i] = 3.0 * g_zz_xyyzzz_0[i] * fbe_0 - 3.0 * g_zz_xyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_xyyzz_1[i] * fe_0 + g_zzz_xyyzzz_1[i] * pa_z[i];

        g_zzzz_xyzzzz_0[i] = 3.0 * g_zz_xyzzzz_0[i] * fbe_0 - 3.0 * g_zz_xyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_xyzzz_1[i] * fe_0 + g_zzz_xyzzzz_1[i] * pa_z[i];

        g_zzzz_xzzzzz_0[i] = 3.0 * g_zz_xzzzzz_0[i] * fbe_0 - 3.0 * g_zz_xzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_xzzzz_1[i] * fe_0 + g_zzz_xzzzzz_1[i] * pa_z[i];

        g_zzzz_yyyyyy_0[i] = 3.0 * g_zz_yyyyyy_0[i] * fbe_0 - 3.0 * g_zz_yyyyyy_1[i] * fz_be_0 + g_zzz_yyyyyy_1[i] * pa_z[i];

        g_zzzz_yyyyyz_0[i] = 3.0 * g_zz_yyyyyz_0[i] * fbe_0 - 3.0 * g_zz_yyyyyz_1[i] * fz_be_0 + g_zzz_yyyyy_1[i] * fe_0 + g_zzz_yyyyyz_1[i] * pa_z[i];

        g_zzzz_yyyyzz_0[i] = 3.0 * g_zz_yyyyzz_0[i] * fbe_0 - 3.0 * g_zz_yyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_yyyyz_1[i] * fe_0 + g_zzz_yyyyzz_1[i] * pa_z[i];

        g_zzzz_yyyzzz_0[i] = 3.0 * g_zz_yyyzzz_0[i] * fbe_0 - 3.0 * g_zz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_yyyzz_1[i] * fe_0 + g_zzz_yyyzzz_1[i] * pa_z[i];

        g_zzzz_yyzzzz_0[i] = 3.0 * g_zz_yyzzzz_0[i] * fbe_0 - 3.0 * g_zz_yyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_yyzzz_1[i] * fe_0 + g_zzz_yyzzzz_1[i] * pa_z[i];

        g_zzzz_yzzzzz_0[i] = 3.0 * g_zz_yzzzzz_0[i] * fbe_0 - 3.0 * g_zz_yzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_yzzzz_1[i] * fe_0 + g_zzz_yzzzzz_1[i] * pa_z[i];

        g_zzzz_zzzzzz_0[i] = 3.0 * g_zz_zzzzzz_0[i] * fbe_0 - 3.0 * g_zz_zzzzzz_1[i] * fz_be_0 + 6.0 * g_zzz_zzzzz_1[i] * fe_0 + g_zzz_zzzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

