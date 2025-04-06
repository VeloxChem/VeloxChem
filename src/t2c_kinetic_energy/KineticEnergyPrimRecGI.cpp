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

#include "KineticEnergyPrimRecGI.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_gi(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_gi,
                            const size_t              idx_ovl_di,
                            const size_t              idx_kin_di,
                            const size_t              idx_kin_fh,
                            const size_t              idx_kin_fi,
                            const size_t              idx_ovl_gi,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpa,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : DI

    auto ts_xx_xxxxxx = pbuffer.data(idx_ovl_di);

    auto ts_xx_xxxxxy = pbuffer.data(idx_ovl_di + 1);

    auto ts_xx_xxxxxz = pbuffer.data(idx_ovl_di + 2);

    auto ts_xx_xxxxyy = pbuffer.data(idx_ovl_di + 3);

    auto ts_xx_xxxxyz = pbuffer.data(idx_ovl_di + 4);

    auto ts_xx_xxxxzz = pbuffer.data(idx_ovl_di + 5);

    auto ts_xx_xxxyyy = pbuffer.data(idx_ovl_di + 6);

    auto ts_xx_xxxyyz = pbuffer.data(idx_ovl_di + 7);

    auto ts_xx_xxxyzz = pbuffer.data(idx_ovl_di + 8);

    auto ts_xx_xxxzzz = pbuffer.data(idx_ovl_di + 9);

    auto ts_xx_xxyyyy = pbuffer.data(idx_ovl_di + 10);

    auto ts_xx_xxyyyz = pbuffer.data(idx_ovl_di + 11);

    auto ts_xx_xxyyzz = pbuffer.data(idx_ovl_di + 12);

    auto ts_xx_xxyzzz = pbuffer.data(idx_ovl_di + 13);

    auto ts_xx_xxzzzz = pbuffer.data(idx_ovl_di + 14);

    auto ts_xx_xyyyyy = pbuffer.data(idx_ovl_di + 15);

    auto ts_xx_xyyyyz = pbuffer.data(idx_ovl_di + 16);

    auto ts_xx_xyyyzz = pbuffer.data(idx_ovl_di + 17);

    auto ts_xx_xyyzzz = pbuffer.data(idx_ovl_di + 18);

    auto ts_xx_xyzzzz = pbuffer.data(idx_ovl_di + 19);

    auto ts_xx_xzzzzz = pbuffer.data(idx_ovl_di + 20);

    auto ts_xx_yyyyyy = pbuffer.data(idx_ovl_di + 21);

    auto ts_xx_yyyyyz = pbuffer.data(idx_ovl_di + 22);

    auto ts_xx_yyyyzz = pbuffer.data(idx_ovl_di + 23);

    auto ts_xx_yyyzzz = pbuffer.data(idx_ovl_di + 24);

    auto ts_xx_yyzzzz = pbuffer.data(idx_ovl_di + 25);

    auto ts_xx_yzzzzz = pbuffer.data(idx_ovl_di + 26);

    auto ts_xx_zzzzzz = pbuffer.data(idx_ovl_di + 27);

    auto ts_yy_xxxxxx = pbuffer.data(idx_ovl_di + 84);

    auto ts_yy_xxxxxy = pbuffer.data(idx_ovl_di + 85);

    auto ts_yy_xxxxxz = pbuffer.data(idx_ovl_di + 86);

    auto ts_yy_xxxxyy = pbuffer.data(idx_ovl_di + 87);

    auto ts_yy_xxxxyz = pbuffer.data(idx_ovl_di + 88);

    auto ts_yy_xxxxzz = pbuffer.data(idx_ovl_di + 89);

    auto ts_yy_xxxyyy = pbuffer.data(idx_ovl_di + 90);

    auto ts_yy_xxxyyz = pbuffer.data(idx_ovl_di + 91);

    auto ts_yy_xxxyzz = pbuffer.data(idx_ovl_di + 92);

    auto ts_yy_xxxzzz = pbuffer.data(idx_ovl_di + 93);

    auto ts_yy_xxyyyy = pbuffer.data(idx_ovl_di + 94);

    auto ts_yy_xxyyyz = pbuffer.data(idx_ovl_di + 95);

    auto ts_yy_xxyyzz = pbuffer.data(idx_ovl_di + 96);

    auto ts_yy_xxyzzz = pbuffer.data(idx_ovl_di + 97);

    auto ts_yy_xxzzzz = pbuffer.data(idx_ovl_di + 98);

    auto ts_yy_xyyyyy = pbuffer.data(idx_ovl_di + 99);

    auto ts_yy_xyyyyz = pbuffer.data(idx_ovl_di + 100);

    auto ts_yy_xyyyzz = pbuffer.data(idx_ovl_di + 101);

    auto ts_yy_xyyzzz = pbuffer.data(idx_ovl_di + 102);

    auto ts_yy_xyzzzz = pbuffer.data(idx_ovl_di + 103);

    auto ts_yy_xzzzzz = pbuffer.data(idx_ovl_di + 104);

    auto ts_yy_yyyyyy = pbuffer.data(idx_ovl_di + 105);

    auto ts_yy_yyyyyz = pbuffer.data(idx_ovl_di + 106);

    auto ts_yy_yyyyzz = pbuffer.data(idx_ovl_di + 107);

    auto ts_yy_yyyzzz = pbuffer.data(idx_ovl_di + 108);

    auto ts_yy_yyzzzz = pbuffer.data(idx_ovl_di + 109);

    auto ts_yy_yzzzzz = pbuffer.data(idx_ovl_di + 110);

    auto ts_yy_zzzzzz = pbuffer.data(idx_ovl_di + 111);

    auto ts_zz_xxxxxx = pbuffer.data(idx_ovl_di + 140);

    auto ts_zz_xxxxxy = pbuffer.data(idx_ovl_di + 141);

    auto ts_zz_xxxxxz = pbuffer.data(idx_ovl_di + 142);

    auto ts_zz_xxxxyy = pbuffer.data(idx_ovl_di + 143);

    auto ts_zz_xxxxyz = pbuffer.data(idx_ovl_di + 144);

    auto ts_zz_xxxxzz = pbuffer.data(idx_ovl_di + 145);

    auto ts_zz_xxxyyy = pbuffer.data(idx_ovl_di + 146);

    auto ts_zz_xxxyyz = pbuffer.data(idx_ovl_di + 147);

    auto ts_zz_xxxyzz = pbuffer.data(idx_ovl_di + 148);

    auto ts_zz_xxxzzz = pbuffer.data(idx_ovl_di + 149);

    auto ts_zz_xxyyyy = pbuffer.data(idx_ovl_di + 150);

    auto ts_zz_xxyyyz = pbuffer.data(idx_ovl_di + 151);

    auto ts_zz_xxyyzz = pbuffer.data(idx_ovl_di + 152);

    auto ts_zz_xxyzzz = pbuffer.data(idx_ovl_di + 153);

    auto ts_zz_xxzzzz = pbuffer.data(idx_ovl_di + 154);

    auto ts_zz_xyyyyy = pbuffer.data(idx_ovl_di + 155);

    auto ts_zz_xyyyyz = pbuffer.data(idx_ovl_di + 156);

    auto ts_zz_xyyyzz = pbuffer.data(idx_ovl_di + 157);

    auto ts_zz_xyyzzz = pbuffer.data(idx_ovl_di + 158);

    auto ts_zz_xyzzzz = pbuffer.data(idx_ovl_di + 159);

    auto ts_zz_xzzzzz = pbuffer.data(idx_ovl_di + 160);

    auto ts_zz_yyyyyy = pbuffer.data(idx_ovl_di + 161);

    auto ts_zz_yyyyyz = pbuffer.data(idx_ovl_di + 162);

    auto ts_zz_yyyyzz = pbuffer.data(idx_ovl_di + 163);

    auto ts_zz_yyyzzz = pbuffer.data(idx_ovl_di + 164);

    auto ts_zz_yyzzzz = pbuffer.data(idx_ovl_di + 165);

    auto ts_zz_yzzzzz = pbuffer.data(idx_ovl_di + 166);

    auto ts_zz_zzzzzz = pbuffer.data(idx_ovl_di + 167);

    // Set up components of auxiliary buffer : DI

    auto tk_xx_xxxxxx = pbuffer.data(idx_kin_di);

    auto tk_xx_xxxxxy = pbuffer.data(idx_kin_di + 1);

    auto tk_xx_xxxxxz = pbuffer.data(idx_kin_di + 2);

    auto tk_xx_xxxxyy = pbuffer.data(idx_kin_di + 3);

    auto tk_xx_xxxxyz = pbuffer.data(idx_kin_di + 4);

    auto tk_xx_xxxxzz = pbuffer.data(idx_kin_di + 5);

    auto tk_xx_xxxyyy = pbuffer.data(idx_kin_di + 6);

    auto tk_xx_xxxyyz = pbuffer.data(idx_kin_di + 7);

    auto tk_xx_xxxyzz = pbuffer.data(idx_kin_di + 8);

    auto tk_xx_xxxzzz = pbuffer.data(idx_kin_di + 9);

    auto tk_xx_xxyyyy = pbuffer.data(idx_kin_di + 10);

    auto tk_xx_xxyyyz = pbuffer.data(idx_kin_di + 11);

    auto tk_xx_xxyyzz = pbuffer.data(idx_kin_di + 12);

    auto tk_xx_xxyzzz = pbuffer.data(idx_kin_di + 13);

    auto tk_xx_xxzzzz = pbuffer.data(idx_kin_di + 14);

    auto tk_xx_xyyyyy = pbuffer.data(idx_kin_di + 15);

    auto tk_xx_xyyyyz = pbuffer.data(idx_kin_di + 16);

    auto tk_xx_xyyyzz = pbuffer.data(idx_kin_di + 17);

    auto tk_xx_xyyzzz = pbuffer.data(idx_kin_di + 18);

    auto tk_xx_xyzzzz = pbuffer.data(idx_kin_di + 19);

    auto tk_xx_xzzzzz = pbuffer.data(idx_kin_di + 20);

    auto tk_xx_yyyyyy = pbuffer.data(idx_kin_di + 21);

    auto tk_xx_yyyyyz = pbuffer.data(idx_kin_di + 22);

    auto tk_xx_yyyyzz = pbuffer.data(idx_kin_di + 23);

    auto tk_xx_yyyzzz = pbuffer.data(idx_kin_di + 24);

    auto tk_xx_yyzzzz = pbuffer.data(idx_kin_di + 25);

    auto tk_xx_yzzzzz = pbuffer.data(idx_kin_di + 26);

    auto tk_xx_zzzzzz = pbuffer.data(idx_kin_di + 27);

    auto tk_yy_xxxxxx = pbuffer.data(idx_kin_di + 84);

    auto tk_yy_xxxxxy = pbuffer.data(idx_kin_di + 85);

    auto tk_yy_xxxxxz = pbuffer.data(idx_kin_di + 86);

    auto tk_yy_xxxxyy = pbuffer.data(idx_kin_di + 87);

    auto tk_yy_xxxxyz = pbuffer.data(idx_kin_di + 88);

    auto tk_yy_xxxxzz = pbuffer.data(idx_kin_di + 89);

    auto tk_yy_xxxyyy = pbuffer.data(idx_kin_di + 90);

    auto tk_yy_xxxyyz = pbuffer.data(idx_kin_di + 91);

    auto tk_yy_xxxyzz = pbuffer.data(idx_kin_di + 92);

    auto tk_yy_xxxzzz = pbuffer.data(idx_kin_di + 93);

    auto tk_yy_xxyyyy = pbuffer.data(idx_kin_di + 94);

    auto tk_yy_xxyyyz = pbuffer.data(idx_kin_di + 95);

    auto tk_yy_xxyyzz = pbuffer.data(idx_kin_di + 96);

    auto tk_yy_xxyzzz = pbuffer.data(idx_kin_di + 97);

    auto tk_yy_xxzzzz = pbuffer.data(idx_kin_di + 98);

    auto tk_yy_xyyyyy = pbuffer.data(idx_kin_di + 99);

    auto tk_yy_xyyyyz = pbuffer.data(idx_kin_di + 100);

    auto tk_yy_xyyyzz = pbuffer.data(idx_kin_di + 101);

    auto tk_yy_xyyzzz = pbuffer.data(idx_kin_di + 102);

    auto tk_yy_xyzzzz = pbuffer.data(idx_kin_di + 103);

    auto tk_yy_xzzzzz = pbuffer.data(idx_kin_di + 104);

    auto tk_yy_yyyyyy = pbuffer.data(idx_kin_di + 105);

    auto tk_yy_yyyyyz = pbuffer.data(idx_kin_di + 106);

    auto tk_yy_yyyyzz = pbuffer.data(idx_kin_di + 107);

    auto tk_yy_yyyzzz = pbuffer.data(idx_kin_di + 108);

    auto tk_yy_yyzzzz = pbuffer.data(idx_kin_di + 109);

    auto tk_yy_yzzzzz = pbuffer.data(idx_kin_di + 110);

    auto tk_yy_zzzzzz = pbuffer.data(idx_kin_di + 111);

    auto tk_zz_xxxxxx = pbuffer.data(idx_kin_di + 140);

    auto tk_zz_xxxxxy = pbuffer.data(idx_kin_di + 141);

    auto tk_zz_xxxxxz = pbuffer.data(idx_kin_di + 142);

    auto tk_zz_xxxxyy = pbuffer.data(idx_kin_di + 143);

    auto tk_zz_xxxxyz = pbuffer.data(idx_kin_di + 144);

    auto tk_zz_xxxxzz = pbuffer.data(idx_kin_di + 145);

    auto tk_zz_xxxyyy = pbuffer.data(idx_kin_di + 146);

    auto tk_zz_xxxyyz = pbuffer.data(idx_kin_di + 147);

    auto tk_zz_xxxyzz = pbuffer.data(idx_kin_di + 148);

    auto tk_zz_xxxzzz = pbuffer.data(idx_kin_di + 149);

    auto tk_zz_xxyyyy = pbuffer.data(idx_kin_di + 150);

    auto tk_zz_xxyyyz = pbuffer.data(idx_kin_di + 151);

    auto tk_zz_xxyyzz = pbuffer.data(idx_kin_di + 152);

    auto tk_zz_xxyzzz = pbuffer.data(idx_kin_di + 153);

    auto tk_zz_xxzzzz = pbuffer.data(idx_kin_di + 154);

    auto tk_zz_xyyyyy = pbuffer.data(idx_kin_di + 155);

    auto tk_zz_xyyyyz = pbuffer.data(idx_kin_di + 156);

    auto tk_zz_xyyyzz = pbuffer.data(idx_kin_di + 157);

    auto tk_zz_xyyzzz = pbuffer.data(idx_kin_di + 158);

    auto tk_zz_xyzzzz = pbuffer.data(idx_kin_di + 159);

    auto tk_zz_xzzzzz = pbuffer.data(idx_kin_di + 160);

    auto tk_zz_yyyyyy = pbuffer.data(idx_kin_di + 161);

    auto tk_zz_yyyyyz = pbuffer.data(idx_kin_di + 162);

    auto tk_zz_yyyyzz = pbuffer.data(idx_kin_di + 163);

    auto tk_zz_yyyzzz = pbuffer.data(idx_kin_di + 164);

    auto tk_zz_yyzzzz = pbuffer.data(idx_kin_di + 165);

    auto tk_zz_yzzzzz = pbuffer.data(idx_kin_di + 166);

    auto tk_zz_zzzzzz = pbuffer.data(idx_kin_di + 167);

    // Set up components of auxiliary buffer : FH

    auto tk_xxx_xxxxx = pbuffer.data(idx_kin_fh);

    auto tk_xxx_xxxxy = pbuffer.data(idx_kin_fh + 1);

    auto tk_xxx_xxxxz = pbuffer.data(idx_kin_fh + 2);

    auto tk_xxx_xxxyy = pbuffer.data(idx_kin_fh + 3);

    auto tk_xxx_xxxyz = pbuffer.data(idx_kin_fh + 4);

    auto tk_xxx_xxxzz = pbuffer.data(idx_kin_fh + 5);

    auto tk_xxx_xxyyy = pbuffer.data(idx_kin_fh + 6);

    auto tk_xxx_xxyyz = pbuffer.data(idx_kin_fh + 7);

    auto tk_xxx_xxyzz = pbuffer.data(idx_kin_fh + 8);

    auto tk_xxx_xxzzz = pbuffer.data(idx_kin_fh + 9);

    auto tk_xxx_xyyyy = pbuffer.data(idx_kin_fh + 10);

    auto tk_xxx_xyyyz = pbuffer.data(idx_kin_fh + 11);

    auto tk_xxx_xyyzz = pbuffer.data(idx_kin_fh + 12);

    auto tk_xxx_xyzzz = pbuffer.data(idx_kin_fh + 13);

    auto tk_xxx_xzzzz = pbuffer.data(idx_kin_fh + 14);

    auto tk_xxx_yyyyy = pbuffer.data(idx_kin_fh + 15);

    auto tk_xxx_yyyyz = pbuffer.data(idx_kin_fh + 16);

    auto tk_xxx_yyyzz = pbuffer.data(idx_kin_fh + 17);

    auto tk_xxx_yyzzz = pbuffer.data(idx_kin_fh + 18);

    auto tk_xxx_yzzzz = pbuffer.data(idx_kin_fh + 19);

    auto tk_xxx_zzzzz = pbuffer.data(idx_kin_fh + 20);

    auto tk_xxz_xxxxz = pbuffer.data(idx_kin_fh + 44);

    auto tk_xxz_xxxyz = pbuffer.data(idx_kin_fh + 46);

    auto tk_xxz_xxxzz = pbuffer.data(idx_kin_fh + 47);

    auto tk_xxz_xxyyz = pbuffer.data(idx_kin_fh + 49);

    auto tk_xxz_xxyzz = pbuffer.data(idx_kin_fh + 50);

    auto tk_xxz_xxzzz = pbuffer.data(idx_kin_fh + 51);

    auto tk_xxz_xyyyz = pbuffer.data(idx_kin_fh + 53);

    auto tk_xxz_xyyzz = pbuffer.data(idx_kin_fh + 54);

    auto tk_xxz_xyzzz = pbuffer.data(idx_kin_fh + 55);

    auto tk_xxz_xzzzz = pbuffer.data(idx_kin_fh + 56);

    auto tk_xxz_yyyyz = pbuffer.data(idx_kin_fh + 58);

    auto tk_xxz_yyyzz = pbuffer.data(idx_kin_fh + 59);

    auto tk_xxz_yyzzz = pbuffer.data(idx_kin_fh + 60);

    auto tk_xxz_yzzzz = pbuffer.data(idx_kin_fh + 61);

    auto tk_xxz_zzzzz = pbuffer.data(idx_kin_fh + 62);

    auto tk_xyy_xxxxy = pbuffer.data(idx_kin_fh + 64);

    auto tk_xyy_xxxyy = pbuffer.data(idx_kin_fh + 66);

    auto tk_xyy_xxxyz = pbuffer.data(idx_kin_fh + 67);

    auto tk_xyy_xxyyy = pbuffer.data(idx_kin_fh + 69);

    auto tk_xyy_xxyyz = pbuffer.data(idx_kin_fh + 70);

    auto tk_xyy_xxyzz = pbuffer.data(idx_kin_fh + 71);

    auto tk_xyy_xyyyy = pbuffer.data(idx_kin_fh + 73);

    auto tk_xyy_xyyyz = pbuffer.data(idx_kin_fh + 74);

    auto tk_xyy_xyyzz = pbuffer.data(idx_kin_fh + 75);

    auto tk_xyy_xyzzz = pbuffer.data(idx_kin_fh + 76);

    auto tk_xyy_yyyyy = pbuffer.data(idx_kin_fh + 78);

    auto tk_xyy_yyyyz = pbuffer.data(idx_kin_fh + 79);

    auto tk_xyy_yyyzz = pbuffer.data(idx_kin_fh + 80);

    auto tk_xyy_yyzzz = pbuffer.data(idx_kin_fh + 81);

    auto tk_xyy_yzzzz = pbuffer.data(idx_kin_fh + 82);

    auto tk_xzz_xxxxz = pbuffer.data(idx_kin_fh + 107);

    auto tk_xzz_xxxyz = pbuffer.data(idx_kin_fh + 109);

    auto tk_xzz_xxxzz = pbuffer.data(idx_kin_fh + 110);

    auto tk_xzz_xxyyz = pbuffer.data(idx_kin_fh + 112);

    auto tk_xzz_xxyzz = pbuffer.data(idx_kin_fh + 113);

    auto tk_xzz_xxzzz = pbuffer.data(idx_kin_fh + 114);

    auto tk_xzz_xyyyz = pbuffer.data(idx_kin_fh + 116);

    auto tk_xzz_xyyzz = pbuffer.data(idx_kin_fh + 117);

    auto tk_xzz_xyzzz = pbuffer.data(idx_kin_fh + 118);

    auto tk_xzz_xzzzz = pbuffer.data(idx_kin_fh + 119);

    auto tk_xzz_yyyyz = pbuffer.data(idx_kin_fh + 121);

    auto tk_xzz_yyyzz = pbuffer.data(idx_kin_fh + 122);

    auto tk_xzz_yyzzz = pbuffer.data(idx_kin_fh + 123);

    auto tk_xzz_yzzzz = pbuffer.data(idx_kin_fh + 124);

    auto tk_xzz_zzzzz = pbuffer.data(idx_kin_fh + 125);

    auto tk_yyy_xxxxx = pbuffer.data(idx_kin_fh + 126);

    auto tk_yyy_xxxxy = pbuffer.data(idx_kin_fh + 127);

    auto tk_yyy_xxxxz = pbuffer.data(idx_kin_fh + 128);

    auto tk_yyy_xxxyy = pbuffer.data(idx_kin_fh + 129);

    auto tk_yyy_xxxyz = pbuffer.data(idx_kin_fh + 130);

    auto tk_yyy_xxxzz = pbuffer.data(idx_kin_fh + 131);

    auto tk_yyy_xxyyy = pbuffer.data(idx_kin_fh + 132);

    auto tk_yyy_xxyyz = pbuffer.data(idx_kin_fh + 133);

    auto tk_yyy_xxyzz = pbuffer.data(idx_kin_fh + 134);

    auto tk_yyy_xxzzz = pbuffer.data(idx_kin_fh + 135);

    auto tk_yyy_xyyyy = pbuffer.data(idx_kin_fh + 136);

    auto tk_yyy_xyyyz = pbuffer.data(idx_kin_fh + 137);

    auto tk_yyy_xyyzz = pbuffer.data(idx_kin_fh + 138);

    auto tk_yyy_xyzzz = pbuffer.data(idx_kin_fh + 139);

    auto tk_yyy_xzzzz = pbuffer.data(idx_kin_fh + 140);

    auto tk_yyy_yyyyy = pbuffer.data(idx_kin_fh + 141);

    auto tk_yyy_yyyyz = pbuffer.data(idx_kin_fh + 142);

    auto tk_yyy_yyyzz = pbuffer.data(idx_kin_fh + 143);

    auto tk_yyy_yyzzz = pbuffer.data(idx_kin_fh + 144);

    auto tk_yyy_yzzzz = pbuffer.data(idx_kin_fh + 145);

    auto tk_yyy_zzzzz = pbuffer.data(idx_kin_fh + 146);

    auto tk_yyz_xxxxz = pbuffer.data(idx_kin_fh + 149);

    auto tk_yyz_xxxyz = pbuffer.data(idx_kin_fh + 151);

    auto tk_yyz_xxxzz = pbuffer.data(idx_kin_fh + 152);

    auto tk_yyz_xxyyz = pbuffer.data(idx_kin_fh + 154);

    auto tk_yyz_xxyzz = pbuffer.data(idx_kin_fh + 155);

    auto tk_yyz_xxzzz = pbuffer.data(idx_kin_fh + 156);

    auto tk_yyz_xyyyz = pbuffer.data(idx_kin_fh + 158);

    auto tk_yyz_xyyzz = pbuffer.data(idx_kin_fh + 159);

    auto tk_yyz_xyzzz = pbuffer.data(idx_kin_fh + 160);

    auto tk_yyz_xzzzz = pbuffer.data(idx_kin_fh + 161);

    auto tk_yyz_yyyyz = pbuffer.data(idx_kin_fh + 163);

    auto tk_yyz_yyyzz = pbuffer.data(idx_kin_fh + 164);

    auto tk_yyz_yyzzz = pbuffer.data(idx_kin_fh + 165);

    auto tk_yyz_yzzzz = pbuffer.data(idx_kin_fh + 166);

    auto tk_yyz_zzzzz = pbuffer.data(idx_kin_fh + 167);

    auto tk_yzz_xxxxy = pbuffer.data(idx_kin_fh + 169);

    auto tk_yzz_xxxxz = pbuffer.data(idx_kin_fh + 170);

    auto tk_yzz_xxxyy = pbuffer.data(idx_kin_fh + 171);

    auto tk_yzz_xxxyz = pbuffer.data(idx_kin_fh + 172);

    auto tk_yzz_xxxzz = pbuffer.data(idx_kin_fh + 173);

    auto tk_yzz_xxyyy = pbuffer.data(idx_kin_fh + 174);

    auto tk_yzz_xxyyz = pbuffer.data(idx_kin_fh + 175);

    auto tk_yzz_xxyzz = pbuffer.data(idx_kin_fh + 176);

    auto tk_yzz_xxzzz = pbuffer.data(idx_kin_fh + 177);

    auto tk_yzz_xyyyy = pbuffer.data(idx_kin_fh + 178);

    auto tk_yzz_xyyyz = pbuffer.data(idx_kin_fh + 179);

    auto tk_yzz_xyyzz = pbuffer.data(idx_kin_fh + 180);

    auto tk_yzz_xyzzz = pbuffer.data(idx_kin_fh + 181);

    auto tk_yzz_xzzzz = pbuffer.data(idx_kin_fh + 182);

    auto tk_yzz_yyyyy = pbuffer.data(idx_kin_fh + 183);

    auto tk_yzz_yyyyz = pbuffer.data(idx_kin_fh + 184);

    auto tk_yzz_yyyzz = pbuffer.data(idx_kin_fh + 185);

    auto tk_yzz_yyzzz = pbuffer.data(idx_kin_fh + 186);

    auto tk_yzz_yzzzz = pbuffer.data(idx_kin_fh + 187);

    auto tk_yzz_zzzzz = pbuffer.data(idx_kin_fh + 188);

    auto tk_zzz_xxxxx = pbuffer.data(idx_kin_fh + 189);

    auto tk_zzz_xxxxy = pbuffer.data(idx_kin_fh + 190);

    auto tk_zzz_xxxxz = pbuffer.data(idx_kin_fh + 191);

    auto tk_zzz_xxxyy = pbuffer.data(idx_kin_fh + 192);

    auto tk_zzz_xxxyz = pbuffer.data(idx_kin_fh + 193);

    auto tk_zzz_xxxzz = pbuffer.data(idx_kin_fh + 194);

    auto tk_zzz_xxyyy = pbuffer.data(idx_kin_fh + 195);

    auto tk_zzz_xxyyz = pbuffer.data(idx_kin_fh + 196);

    auto tk_zzz_xxyzz = pbuffer.data(idx_kin_fh + 197);

    auto tk_zzz_xxzzz = pbuffer.data(idx_kin_fh + 198);

    auto tk_zzz_xyyyy = pbuffer.data(idx_kin_fh + 199);

    auto tk_zzz_xyyyz = pbuffer.data(idx_kin_fh + 200);

    auto tk_zzz_xyyzz = pbuffer.data(idx_kin_fh + 201);

    auto tk_zzz_xyzzz = pbuffer.data(idx_kin_fh + 202);

    auto tk_zzz_xzzzz = pbuffer.data(idx_kin_fh + 203);

    auto tk_zzz_yyyyy = pbuffer.data(idx_kin_fh + 204);

    auto tk_zzz_yyyyz = pbuffer.data(idx_kin_fh + 205);

    auto tk_zzz_yyyzz = pbuffer.data(idx_kin_fh + 206);

    auto tk_zzz_yyzzz = pbuffer.data(idx_kin_fh + 207);

    auto tk_zzz_yzzzz = pbuffer.data(idx_kin_fh + 208);

    auto tk_zzz_zzzzz = pbuffer.data(idx_kin_fh + 209);

    // Set up components of auxiliary buffer : FI

    auto tk_xxx_xxxxxx = pbuffer.data(idx_kin_fi);

    auto tk_xxx_xxxxxy = pbuffer.data(idx_kin_fi + 1);

    auto tk_xxx_xxxxxz = pbuffer.data(idx_kin_fi + 2);

    auto tk_xxx_xxxxyy = pbuffer.data(idx_kin_fi + 3);

    auto tk_xxx_xxxxyz = pbuffer.data(idx_kin_fi + 4);

    auto tk_xxx_xxxxzz = pbuffer.data(idx_kin_fi + 5);

    auto tk_xxx_xxxyyy = pbuffer.data(idx_kin_fi + 6);

    auto tk_xxx_xxxyyz = pbuffer.data(idx_kin_fi + 7);

    auto tk_xxx_xxxyzz = pbuffer.data(idx_kin_fi + 8);

    auto tk_xxx_xxxzzz = pbuffer.data(idx_kin_fi + 9);

    auto tk_xxx_xxyyyy = pbuffer.data(idx_kin_fi + 10);

    auto tk_xxx_xxyyyz = pbuffer.data(idx_kin_fi + 11);

    auto tk_xxx_xxyyzz = pbuffer.data(idx_kin_fi + 12);

    auto tk_xxx_xxyzzz = pbuffer.data(idx_kin_fi + 13);

    auto tk_xxx_xxzzzz = pbuffer.data(idx_kin_fi + 14);

    auto tk_xxx_xyyyyy = pbuffer.data(idx_kin_fi + 15);

    auto tk_xxx_xyyyyz = pbuffer.data(idx_kin_fi + 16);

    auto tk_xxx_xyyyzz = pbuffer.data(idx_kin_fi + 17);

    auto tk_xxx_xyyzzz = pbuffer.data(idx_kin_fi + 18);

    auto tk_xxx_xyzzzz = pbuffer.data(idx_kin_fi + 19);

    auto tk_xxx_xzzzzz = pbuffer.data(idx_kin_fi + 20);

    auto tk_xxx_yyyyyy = pbuffer.data(idx_kin_fi + 21);

    auto tk_xxx_yyyyyz = pbuffer.data(idx_kin_fi + 22);

    auto tk_xxx_yyyyzz = pbuffer.data(idx_kin_fi + 23);

    auto tk_xxx_yyyzzz = pbuffer.data(idx_kin_fi + 24);

    auto tk_xxx_yyzzzz = pbuffer.data(idx_kin_fi + 25);

    auto tk_xxx_yzzzzz = pbuffer.data(idx_kin_fi + 26);

    auto tk_xxx_zzzzzz = pbuffer.data(idx_kin_fi + 27);

    auto tk_xxy_xxxxxx = pbuffer.data(idx_kin_fi + 28);

    auto tk_xxy_xxxxxy = pbuffer.data(idx_kin_fi + 29);

    auto tk_xxy_xxxxxz = pbuffer.data(idx_kin_fi + 30);

    auto tk_xxy_xxxxyy = pbuffer.data(idx_kin_fi + 31);

    auto tk_xxy_xxxxzz = pbuffer.data(idx_kin_fi + 33);

    auto tk_xxy_xxxyyy = pbuffer.data(idx_kin_fi + 34);

    auto tk_xxy_xxxzzz = pbuffer.data(idx_kin_fi + 37);

    auto tk_xxy_xxyyyy = pbuffer.data(idx_kin_fi + 38);

    auto tk_xxy_xxzzzz = pbuffer.data(idx_kin_fi + 42);

    auto tk_xxy_xyyyyy = pbuffer.data(idx_kin_fi + 43);

    auto tk_xxy_xzzzzz = pbuffer.data(idx_kin_fi + 48);

    auto tk_xxy_yyyyyy = pbuffer.data(idx_kin_fi + 49);

    auto tk_xxz_xxxxxx = pbuffer.data(idx_kin_fi + 56);

    auto tk_xxz_xxxxxy = pbuffer.data(idx_kin_fi + 57);

    auto tk_xxz_xxxxxz = pbuffer.data(idx_kin_fi + 58);

    auto tk_xxz_xxxxyy = pbuffer.data(idx_kin_fi + 59);

    auto tk_xxz_xxxxyz = pbuffer.data(idx_kin_fi + 60);

    auto tk_xxz_xxxxzz = pbuffer.data(idx_kin_fi + 61);

    auto tk_xxz_xxxyyy = pbuffer.data(idx_kin_fi + 62);

    auto tk_xxz_xxxyyz = pbuffer.data(idx_kin_fi + 63);

    auto tk_xxz_xxxyzz = pbuffer.data(idx_kin_fi + 64);

    auto tk_xxz_xxxzzz = pbuffer.data(idx_kin_fi + 65);

    auto tk_xxz_xxyyyy = pbuffer.data(idx_kin_fi + 66);

    auto tk_xxz_xxyyyz = pbuffer.data(idx_kin_fi + 67);

    auto tk_xxz_xxyyzz = pbuffer.data(idx_kin_fi + 68);

    auto tk_xxz_xxyzzz = pbuffer.data(idx_kin_fi + 69);

    auto tk_xxz_xxzzzz = pbuffer.data(idx_kin_fi + 70);

    auto tk_xxz_xyyyyy = pbuffer.data(idx_kin_fi + 71);

    auto tk_xxz_xyyyyz = pbuffer.data(idx_kin_fi + 72);

    auto tk_xxz_xyyyzz = pbuffer.data(idx_kin_fi + 73);

    auto tk_xxz_xyyzzz = pbuffer.data(idx_kin_fi + 74);

    auto tk_xxz_xyzzzz = pbuffer.data(idx_kin_fi + 75);

    auto tk_xxz_xzzzzz = pbuffer.data(idx_kin_fi + 76);

    auto tk_xxz_yyyyyz = pbuffer.data(idx_kin_fi + 78);

    auto tk_xxz_yyyyzz = pbuffer.data(idx_kin_fi + 79);

    auto tk_xxz_yyyzzz = pbuffer.data(idx_kin_fi + 80);

    auto tk_xxz_yyzzzz = pbuffer.data(idx_kin_fi + 81);

    auto tk_xxz_yzzzzz = pbuffer.data(idx_kin_fi + 82);

    auto tk_xxz_zzzzzz = pbuffer.data(idx_kin_fi + 83);

    auto tk_xyy_xxxxxx = pbuffer.data(idx_kin_fi + 84);

    auto tk_xyy_xxxxxy = pbuffer.data(idx_kin_fi + 85);

    auto tk_xyy_xxxxyy = pbuffer.data(idx_kin_fi + 87);

    auto tk_xyy_xxxxyz = pbuffer.data(idx_kin_fi + 88);

    auto tk_xyy_xxxyyy = pbuffer.data(idx_kin_fi + 90);

    auto tk_xyy_xxxyyz = pbuffer.data(idx_kin_fi + 91);

    auto tk_xyy_xxxyzz = pbuffer.data(idx_kin_fi + 92);

    auto tk_xyy_xxyyyy = pbuffer.data(idx_kin_fi + 94);

    auto tk_xyy_xxyyyz = pbuffer.data(idx_kin_fi + 95);

    auto tk_xyy_xxyyzz = pbuffer.data(idx_kin_fi + 96);

    auto tk_xyy_xxyzzz = pbuffer.data(idx_kin_fi + 97);

    auto tk_xyy_xyyyyy = pbuffer.data(idx_kin_fi + 99);

    auto tk_xyy_xyyyyz = pbuffer.data(idx_kin_fi + 100);

    auto tk_xyy_xyyyzz = pbuffer.data(idx_kin_fi + 101);

    auto tk_xyy_xyyzzz = pbuffer.data(idx_kin_fi + 102);

    auto tk_xyy_xyzzzz = pbuffer.data(idx_kin_fi + 103);

    auto tk_xyy_yyyyyy = pbuffer.data(idx_kin_fi + 105);

    auto tk_xyy_yyyyyz = pbuffer.data(idx_kin_fi + 106);

    auto tk_xyy_yyyyzz = pbuffer.data(idx_kin_fi + 107);

    auto tk_xyy_yyyzzz = pbuffer.data(idx_kin_fi + 108);

    auto tk_xyy_yyzzzz = pbuffer.data(idx_kin_fi + 109);

    auto tk_xyy_yzzzzz = pbuffer.data(idx_kin_fi + 110);

    auto tk_xyy_zzzzzz = pbuffer.data(idx_kin_fi + 111);

    auto tk_xzz_xxxxxx = pbuffer.data(idx_kin_fi + 140);

    auto tk_xzz_xxxxxz = pbuffer.data(idx_kin_fi + 142);

    auto tk_xzz_xxxxyz = pbuffer.data(idx_kin_fi + 144);

    auto tk_xzz_xxxxzz = pbuffer.data(idx_kin_fi + 145);

    auto tk_xzz_xxxyyz = pbuffer.data(idx_kin_fi + 147);

    auto tk_xzz_xxxyzz = pbuffer.data(idx_kin_fi + 148);

    auto tk_xzz_xxxzzz = pbuffer.data(idx_kin_fi + 149);

    auto tk_xzz_xxyyyz = pbuffer.data(idx_kin_fi + 151);

    auto tk_xzz_xxyyzz = pbuffer.data(idx_kin_fi + 152);

    auto tk_xzz_xxyzzz = pbuffer.data(idx_kin_fi + 153);

    auto tk_xzz_xxzzzz = pbuffer.data(idx_kin_fi + 154);

    auto tk_xzz_xyyyyz = pbuffer.data(idx_kin_fi + 156);

    auto tk_xzz_xyyyzz = pbuffer.data(idx_kin_fi + 157);

    auto tk_xzz_xyyzzz = pbuffer.data(idx_kin_fi + 158);

    auto tk_xzz_xyzzzz = pbuffer.data(idx_kin_fi + 159);

    auto tk_xzz_xzzzzz = pbuffer.data(idx_kin_fi + 160);

    auto tk_xzz_yyyyyy = pbuffer.data(idx_kin_fi + 161);

    auto tk_xzz_yyyyyz = pbuffer.data(idx_kin_fi + 162);

    auto tk_xzz_yyyyzz = pbuffer.data(idx_kin_fi + 163);

    auto tk_xzz_yyyzzz = pbuffer.data(idx_kin_fi + 164);

    auto tk_xzz_yyzzzz = pbuffer.data(idx_kin_fi + 165);

    auto tk_xzz_yzzzzz = pbuffer.data(idx_kin_fi + 166);

    auto tk_xzz_zzzzzz = pbuffer.data(idx_kin_fi + 167);

    auto tk_yyy_xxxxxx = pbuffer.data(idx_kin_fi + 168);

    auto tk_yyy_xxxxxy = pbuffer.data(idx_kin_fi + 169);

    auto tk_yyy_xxxxxz = pbuffer.data(idx_kin_fi + 170);

    auto tk_yyy_xxxxyy = pbuffer.data(idx_kin_fi + 171);

    auto tk_yyy_xxxxyz = pbuffer.data(idx_kin_fi + 172);

    auto tk_yyy_xxxxzz = pbuffer.data(idx_kin_fi + 173);

    auto tk_yyy_xxxyyy = pbuffer.data(idx_kin_fi + 174);

    auto tk_yyy_xxxyyz = pbuffer.data(idx_kin_fi + 175);

    auto tk_yyy_xxxyzz = pbuffer.data(idx_kin_fi + 176);

    auto tk_yyy_xxxzzz = pbuffer.data(idx_kin_fi + 177);

    auto tk_yyy_xxyyyy = pbuffer.data(idx_kin_fi + 178);

    auto tk_yyy_xxyyyz = pbuffer.data(idx_kin_fi + 179);

    auto tk_yyy_xxyyzz = pbuffer.data(idx_kin_fi + 180);

    auto tk_yyy_xxyzzz = pbuffer.data(idx_kin_fi + 181);

    auto tk_yyy_xxzzzz = pbuffer.data(idx_kin_fi + 182);

    auto tk_yyy_xyyyyy = pbuffer.data(idx_kin_fi + 183);

    auto tk_yyy_xyyyyz = pbuffer.data(idx_kin_fi + 184);

    auto tk_yyy_xyyyzz = pbuffer.data(idx_kin_fi + 185);

    auto tk_yyy_xyyzzz = pbuffer.data(idx_kin_fi + 186);

    auto tk_yyy_xyzzzz = pbuffer.data(idx_kin_fi + 187);

    auto tk_yyy_xzzzzz = pbuffer.data(idx_kin_fi + 188);

    auto tk_yyy_yyyyyy = pbuffer.data(idx_kin_fi + 189);

    auto tk_yyy_yyyyyz = pbuffer.data(idx_kin_fi + 190);

    auto tk_yyy_yyyyzz = pbuffer.data(idx_kin_fi + 191);

    auto tk_yyy_yyyzzz = pbuffer.data(idx_kin_fi + 192);

    auto tk_yyy_yyzzzz = pbuffer.data(idx_kin_fi + 193);

    auto tk_yyy_yzzzzz = pbuffer.data(idx_kin_fi + 194);

    auto tk_yyy_zzzzzz = pbuffer.data(idx_kin_fi + 195);

    auto tk_yyz_xxxxxy = pbuffer.data(idx_kin_fi + 197);

    auto tk_yyz_xxxxxz = pbuffer.data(idx_kin_fi + 198);

    auto tk_yyz_xxxxyy = pbuffer.data(idx_kin_fi + 199);

    auto tk_yyz_xxxxyz = pbuffer.data(idx_kin_fi + 200);

    auto tk_yyz_xxxxzz = pbuffer.data(idx_kin_fi + 201);

    auto tk_yyz_xxxyyy = pbuffer.data(idx_kin_fi + 202);

    auto tk_yyz_xxxyyz = pbuffer.data(idx_kin_fi + 203);

    auto tk_yyz_xxxyzz = pbuffer.data(idx_kin_fi + 204);

    auto tk_yyz_xxxzzz = pbuffer.data(idx_kin_fi + 205);

    auto tk_yyz_xxyyyy = pbuffer.data(idx_kin_fi + 206);

    auto tk_yyz_xxyyyz = pbuffer.data(idx_kin_fi + 207);

    auto tk_yyz_xxyyzz = pbuffer.data(idx_kin_fi + 208);

    auto tk_yyz_xxyzzz = pbuffer.data(idx_kin_fi + 209);

    auto tk_yyz_xxzzzz = pbuffer.data(idx_kin_fi + 210);

    auto tk_yyz_xyyyyy = pbuffer.data(idx_kin_fi + 211);

    auto tk_yyz_xyyyyz = pbuffer.data(idx_kin_fi + 212);

    auto tk_yyz_xyyyzz = pbuffer.data(idx_kin_fi + 213);

    auto tk_yyz_xyyzzz = pbuffer.data(idx_kin_fi + 214);

    auto tk_yyz_xyzzzz = pbuffer.data(idx_kin_fi + 215);

    auto tk_yyz_xzzzzz = pbuffer.data(idx_kin_fi + 216);

    auto tk_yyz_yyyyyy = pbuffer.data(idx_kin_fi + 217);

    auto tk_yyz_yyyyyz = pbuffer.data(idx_kin_fi + 218);

    auto tk_yyz_yyyyzz = pbuffer.data(idx_kin_fi + 219);

    auto tk_yyz_yyyzzz = pbuffer.data(idx_kin_fi + 220);

    auto tk_yyz_yyzzzz = pbuffer.data(idx_kin_fi + 221);

    auto tk_yyz_yzzzzz = pbuffer.data(idx_kin_fi + 222);

    auto tk_yyz_zzzzzz = pbuffer.data(idx_kin_fi + 223);

    auto tk_yzz_xxxxxx = pbuffer.data(idx_kin_fi + 224);

    auto tk_yzz_xxxxxy = pbuffer.data(idx_kin_fi + 225);

    auto tk_yzz_xxxxxz = pbuffer.data(idx_kin_fi + 226);

    auto tk_yzz_xxxxyy = pbuffer.data(idx_kin_fi + 227);

    auto tk_yzz_xxxxyz = pbuffer.data(idx_kin_fi + 228);

    auto tk_yzz_xxxxzz = pbuffer.data(idx_kin_fi + 229);

    auto tk_yzz_xxxyyy = pbuffer.data(idx_kin_fi + 230);

    auto tk_yzz_xxxyyz = pbuffer.data(idx_kin_fi + 231);

    auto tk_yzz_xxxyzz = pbuffer.data(idx_kin_fi + 232);

    auto tk_yzz_xxxzzz = pbuffer.data(idx_kin_fi + 233);

    auto tk_yzz_xxyyyy = pbuffer.data(idx_kin_fi + 234);

    auto tk_yzz_xxyyyz = pbuffer.data(idx_kin_fi + 235);

    auto tk_yzz_xxyyzz = pbuffer.data(idx_kin_fi + 236);

    auto tk_yzz_xxyzzz = pbuffer.data(idx_kin_fi + 237);

    auto tk_yzz_xxzzzz = pbuffer.data(idx_kin_fi + 238);

    auto tk_yzz_xyyyyy = pbuffer.data(idx_kin_fi + 239);

    auto tk_yzz_xyyyyz = pbuffer.data(idx_kin_fi + 240);

    auto tk_yzz_xyyyzz = pbuffer.data(idx_kin_fi + 241);

    auto tk_yzz_xyyzzz = pbuffer.data(idx_kin_fi + 242);

    auto tk_yzz_xyzzzz = pbuffer.data(idx_kin_fi + 243);

    auto tk_yzz_xzzzzz = pbuffer.data(idx_kin_fi + 244);

    auto tk_yzz_yyyyyy = pbuffer.data(idx_kin_fi + 245);

    auto tk_yzz_yyyyyz = pbuffer.data(idx_kin_fi + 246);

    auto tk_yzz_yyyyzz = pbuffer.data(idx_kin_fi + 247);

    auto tk_yzz_yyyzzz = pbuffer.data(idx_kin_fi + 248);

    auto tk_yzz_yyzzzz = pbuffer.data(idx_kin_fi + 249);

    auto tk_yzz_yzzzzz = pbuffer.data(idx_kin_fi + 250);

    auto tk_yzz_zzzzzz = pbuffer.data(idx_kin_fi + 251);

    auto tk_zzz_xxxxxx = pbuffer.data(idx_kin_fi + 252);

    auto tk_zzz_xxxxxy = pbuffer.data(idx_kin_fi + 253);

    auto tk_zzz_xxxxxz = pbuffer.data(idx_kin_fi + 254);

    auto tk_zzz_xxxxyy = pbuffer.data(idx_kin_fi + 255);

    auto tk_zzz_xxxxyz = pbuffer.data(idx_kin_fi + 256);

    auto tk_zzz_xxxxzz = pbuffer.data(idx_kin_fi + 257);

    auto tk_zzz_xxxyyy = pbuffer.data(idx_kin_fi + 258);

    auto tk_zzz_xxxyyz = pbuffer.data(idx_kin_fi + 259);

    auto tk_zzz_xxxyzz = pbuffer.data(idx_kin_fi + 260);

    auto tk_zzz_xxxzzz = pbuffer.data(idx_kin_fi + 261);

    auto tk_zzz_xxyyyy = pbuffer.data(idx_kin_fi + 262);

    auto tk_zzz_xxyyyz = pbuffer.data(idx_kin_fi + 263);

    auto tk_zzz_xxyyzz = pbuffer.data(idx_kin_fi + 264);

    auto tk_zzz_xxyzzz = pbuffer.data(idx_kin_fi + 265);

    auto tk_zzz_xxzzzz = pbuffer.data(idx_kin_fi + 266);

    auto tk_zzz_xyyyyy = pbuffer.data(idx_kin_fi + 267);

    auto tk_zzz_xyyyyz = pbuffer.data(idx_kin_fi + 268);

    auto tk_zzz_xyyyzz = pbuffer.data(idx_kin_fi + 269);

    auto tk_zzz_xyyzzz = pbuffer.data(idx_kin_fi + 270);

    auto tk_zzz_xyzzzz = pbuffer.data(idx_kin_fi + 271);

    auto tk_zzz_xzzzzz = pbuffer.data(idx_kin_fi + 272);

    auto tk_zzz_yyyyyy = pbuffer.data(idx_kin_fi + 273);

    auto tk_zzz_yyyyyz = pbuffer.data(idx_kin_fi + 274);

    auto tk_zzz_yyyyzz = pbuffer.data(idx_kin_fi + 275);

    auto tk_zzz_yyyzzz = pbuffer.data(idx_kin_fi + 276);

    auto tk_zzz_yyzzzz = pbuffer.data(idx_kin_fi + 277);

    auto tk_zzz_yzzzzz = pbuffer.data(idx_kin_fi + 278);

    auto tk_zzz_zzzzzz = pbuffer.data(idx_kin_fi + 279);

    // Set up components of auxiliary buffer : GI

    auto ts_xxxx_xxxxxx = pbuffer.data(idx_ovl_gi);

    auto ts_xxxx_xxxxxy = pbuffer.data(idx_ovl_gi + 1);

    auto ts_xxxx_xxxxxz = pbuffer.data(idx_ovl_gi + 2);

    auto ts_xxxx_xxxxyy = pbuffer.data(idx_ovl_gi + 3);

    auto ts_xxxx_xxxxyz = pbuffer.data(idx_ovl_gi + 4);

    auto ts_xxxx_xxxxzz = pbuffer.data(idx_ovl_gi + 5);

    auto ts_xxxx_xxxyyy = pbuffer.data(idx_ovl_gi + 6);

    auto ts_xxxx_xxxyyz = pbuffer.data(idx_ovl_gi + 7);

    auto ts_xxxx_xxxyzz = pbuffer.data(idx_ovl_gi + 8);

    auto ts_xxxx_xxxzzz = pbuffer.data(idx_ovl_gi + 9);

    auto ts_xxxx_xxyyyy = pbuffer.data(idx_ovl_gi + 10);

    auto ts_xxxx_xxyyyz = pbuffer.data(idx_ovl_gi + 11);

    auto ts_xxxx_xxyyzz = pbuffer.data(idx_ovl_gi + 12);

    auto ts_xxxx_xxyzzz = pbuffer.data(idx_ovl_gi + 13);

    auto ts_xxxx_xxzzzz = pbuffer.data(idx_ovl_gi + 14);

    auto ts_xxxx_xyyyyy = pbuffer.data(idx_ovl_gi + 15);

    auto ts_xxxx_xyyyyz = pbuffer.data(idx_ovl_gi + 16);

    auto ts_xxxx_xyyyzz = pbuffer.data(idx_ovl_gi + 17);

    auto ts_xxxx_xyyzzz = pbuffer.data(idx_ovl_gi + 18);

    auto ts_xxxx_xyzzzz = pbuffer.data(idx_ovl_gi + 19);

    auto ts_xxxx_xzzzzz = pbuffer.data(idx_ovl_gi + 20);

    auto ts_xxxx_yyyyyy = pbuffer.data(idx_ovl_gi + 21);

    auto ts_xxxx_yyyyyz = pbuffer.data(idx_ovl_gi + 22);

    auto ts_xxxx_yyyyzz = pbuffer.data(idx_ovl_gi + 23);

    auto ts_xxxx_yyyzzz = pbuffer.data(idx_ovl_gi + 24);

    auto ts_xxxx_yyzzzz = pbuffer.data(idx_ovl_gi + 25);

    auto ts_xxxx_yzzzzz = pbuffer.data(idx_ovl_gi + 26);

    auto ts_xxxx_zzzzzz = pbuffer.data(idx_ovl_gi + 27);

    auto ts_xxxy_xxxxxx = pbuffer.data(idx_ovl_gi + 28);

    auto ts_xxxy_xxxxxy = pbuffer.data(idx_ovl_gi + 29);

    auto ts_xxxy_xxxxxz = pbuffer.data(idx_ovl_gi + 30);

    auto ts_xxxy_xxxxyy = pbuffer.data(idx_ovl_gi + 31);

    auto ts_xxxy_xxxxyz = pbuffer.data(idx_ovl_gi + 32);

    auto ts_xxxy_xxxxzz = pbuffer.data(idx_ovl_gi + 33);

    auto ts_xxxy_xxxyyy = pbuffer.data(idx_ovl_gi + 34);

    auto ts_xxxy_xxxyyz = pbuffer.data(idx_ovl_gi + 35);

    auto ts_xxxy_xxxyzz = pbuffer.data(idx_ovl_gi + 36);

    auto ts_xxxy_xxxzzz = pbuffer.data(idx_ovl_gi + 37);

    auto ts_xxxy_xxyyyy = pbuffer.data(idx_ovl_gi + 38);

    auto ts_xxxy_xxyyyz = pbuffer.data(idx_ovl_gi + 39);

    auto ts_xxxy_xxyyzz = pbuffer.data(idx_ovl_gi + 40);

    auto ts_xxxy_xxyzzz = pbuffer.data(idx_ovl_gi + 41);

    auto ts_xxxy_xxzzzz = pbuffer.data(idx_ovl_gi + 42);

    auto ts_xxxy_xyyyyy = pbuffer.data(idx_ovl_gi + 43);

    auto ts_xxxy_xyyyyz = pbuffer.data(idx_ovl_gi + 44);

    auto ts_xxxy_xyyyzz = pbuffer.data(idx_ovl_gi + 45);

    auto ts_xxxy_xyyzzz = pbuffer.data(idx_ovl_gi + 46);

    auto ts_xxxy_xyzzzz = pbuffer.data(idx_ovl_gi + 47);

    auto ts_xxxy_xzzzzz = pbuffer.data(idx_ovl_gi + 48);

    auto ts_xxxy_yyyyyy = pbuffer.data(idx_ovl_gi + 49);

    auto ts_xxxy_yyyyyz = pbuffer.data(idx_ovl_gi + 50);

    auto ts_xxxy_yyyyzz = pbuffer.data(idx_ovl_gi + 51);

    auto ts_xxxy_yyyzzz = pbuffer.data(idx_ovl_gi + 52);

    auto ts_xxxy_yyzzzz = pbuffer.data(idx_ovl_gi + 53);

    auto ts_xxxy_yzzzzz = pbuffer.data(idx_ovl_gi + 54);

    auto ts_xxxy_zzzzzz = pbuffer.data(idx_ovl_gi + 55);

    auto ts_xxxz_xxxxxx = pbuffer.data(idx_ovl_gi + 56);

    auto ts_xxxz_xxxxxy = pbuffer.data(idx_ovl_gi + 57);

    auto ts_xxxz_xxxxxz = pbuffer.data(idx_ovl_gi + 58);

    auto ts_xxxz_xxxxyy = pbuffer.data(idx_ovl_gi + 59);

    auto ts_xxxz_xxxxyz = pbuffer.data(idx_ovl_gi + 60);

    auto ts_xxxz_xxxxzz = pbuffer.data(idx_ovl_gi + 61);

    auto ts_xxxz_xxxyyy = pbuffer.data(idx_ovl_gi + 62);

    auto ts_xxxz_xxxyyz = pbuffer.data(idx_ovl_gi + 63);

    auto ts_xxxz_xxxyzz = pbuffer.data(idx_ovl_gi + 64);

    auto ts_xxxz_xxxzzz = pbuffer.data(idx_ovl_gi + 65);

    auto ts_xxxz_xxyyyy = pbuffer.data(idx_ovl_gi + 66);

    auto ts_xxxz_xxyyyz = pbuffer.data(idx_ovl_gi + 67);

    auto ts_xxxz_xxyyzz = pbuffer.data(idx_ovl_gi + 68);

    auto ts_xxxz_xxyzzz = pbuffer.data(idx_ovl_gi + 69);

    auto ts_xxxz_xxzzzz = pbuffer.data(idx_ovl_gi + 70);

    auto ts_xxxz_xyyyyy = pbuffer.data(idx_ovl_gi + 71);

    auto ts_xxxz_xyyyyz = pbuffer.data(idx_ovl_gi + 72);

    auto ts_xxxz_xyyyzz = pbuffer.data(idx_ovl_gi + 73);

    auto ts_xxxz_xyyzzz = pbuffer.data(idx_ovl_gi + 74);

    auto ts_xxxz_xyzzzz = pbuffer.data(idx_ovl_gi + 75);

    auto ts_xxxz_xzzzzz = pbuffer.data(idx_ovl_gi + 76);

    auto ts_xxxz_yyyyyy = pbuffer.data(idx_ovl_gi + 77);

    auto ts_xxxz_yyyyyz = pbuffer.data(idx_ovl_gi + 78);

    auto ts_xxxz_yyyyzz = pbuffer.data(idx_ovl_gi + 79);

    auto ts_xxxz_yyyzzz = pbuffer.data(idx_ovl_gi + 80);

    auto ts_xxxz_yyzzzz = pbuffer.data(idx_ovl_gi + 81);

    auto ts_xxxz_yzzzzz = pbuffer.data(idx_ovl_gi + 82);

    auto ts_xxxz_zzzzzz = pbuffer.data(idx_ovl_gi + 83);

    auto ts_xxyy_xxxxxx = pbuffer.data(idx_ovl_gi + 84);

    auto ts_xxyy_xxxxxy = pbuffer.data(idx_ovl_gi + 85);

    auto ts_xxyy_xxxxxz = pbuffer.data(idx_ovl_gi + 86);

    auto ts_xxyy_xxxxyy = pbuffer.data(idx_ovl_gi + 87);

    auto ts_xxyy_xxxxyz = pbuffer.data(idx_ovl_gi + 88);

    auto ts_xxyy_xxxxzz = pbuffer.data(idx_ovl_gi + 89);

    auto ts_xxyy_xxxyyy = pbuffer.data(idx_ovl_gi + 90);

    auto ts_xxyy_xxxyyz = pbuffer.data(idx_ovl_gi + 91);

    auto ts_xxyy_xxxyzz = pbuffer.data(idx_ovl_gi + 92);

    auto ts_xxyy_xxxzzz = pbuffer.data(idx_ovl_gi + 93);

    auto ts_xxyy_xxyyyy = pbuffer.data(idx_ovl_gi + 94);

    auto ts_xxyy_xxyyyz = pbuffer.data(idx_ovl_gi + 95);

    auto ts_xxyy_xxyyzz = pbuffer.data(idx_ovl_gi + 96);

    auto ts_xxyy_xxyzzz = pbuffer.data(idx_ovl_gi + 97);

    auto ts_xxyy_xxzzzz = pbuffer.data(idx_ovl_gi + 98);

    auto ts_xxyy_xyyyyy = pbuffer.data(idx_ovl_gi + 99);

    auto ts_xxyy_xyyyyz = pbuffer.data(idx_ovl_gi + 100);

    auto ts_xxyy_xyyyzz = pbuffer.data(idx_ovl_gi + 101);

    auto ts_xxyy_xyyzzz = pbuffer.data(idx_ovl_gi + 102);

    auto ts_xxyy_xyzzzz = pbuffer.data(idx_ovl_gi + 103);

    auto ts_xxyy_xzzzzz = pbuffer.data(idx_ovl_gi + 104);

    auto ts_xxyy_yyyyyy = pbuffer.data(idx_ovl_gi + 105);

    auto ts_xxyy_yyyyyz = pbuffer.data(idx_ovl_gi + 106);

    auto ts_xxyy_yyyyzz = pbuffer.data(idx_ovl_gi + 107);

    auto ts_xxyy_yyyzzz = pbuffer.data(idx_ovl_gi + 108);

    auto ts_xxyy_yyzzzz = pbuffer.data(idx_ovl_gi + 109);

    auto ts_xxyy_yzzzzz = pbuffer.data(idx_ovl_gi + 110);

    auto ts_xxyy_zzzzzz = pbuffer.data(idx_ovl_gi + 111);

    auto ts_xxyz_xxxxxx = pbuffer.data(idx_ovl_gi + 112);

    auto ts_xxyz_xxxxxy = pbuffer.data(idx_ovl_gi + 113);

    auto ts_xxyz_xxxxxz = pbuffer.data(idx_ovl_gi + 114);

    auto ts_xxyz_xxxxyy = pbuffer.data(idx_ovl_gi + 115);

    auto ts_xxyz_xxxxyz = pbuffer.data(idx_ovl_gi + 116);

    auto ts_xxyz_xxxxzz = pbuffer.data(idx_ovl_gi + 117);

    auto ts_xxyz_xxxyyy = pbuffer.data(idx_ovl_gi + 118);

    auto ts_xxyz_xxxyyz = pbuffer.data(idx_ovl_gi + 119);

    auto ts_xxyz_xxxyzz = pbuffer.data(idx_ovl_gi + 120);

    auto ts_xxyz_xxxzzz = pbuffer.data(idx_ovl_gi + 121);

    auto ts_xxyz_xxyyyy = pbuffer.data(idx_ovl_gi + 122);

    auto ts_xxyz_xxyyyz = pbuffer.data(idx_ovl_gi + 123);

    auto ts_xxyz_xxyyzz = pbuffer.data(idx_ovl_gi + 124);

    auto ts_xxyz_xxyzzz = pbuffer.data(idx_ovl_gi + 125);

    auto ts_xxyz_xxzzzz = pbuffer.data(idx_ovl_gi + 126);

    auto ts_xxyz_xyyyyy = pbuffer.data(idx_ovl_gi + 127);

    auto ts_xxyz_xyyyyz = pbuffer.data(idx_ovl_gi + 128);

    auto ts_xxyz_xyyyzz = pbuffer.data(idx_ovl_gi + 129);

    auto ts_xxyz_xyyzzz = pbuffer.data(idx_ovl_gi + 130);

    auto ts_xxyz_xyzzzz = pbuffer.data(idx_ovl_gi + 131);

    auto ts_xxyz_xzzzzz = pbuffer.data(idx_ovl_gi + 132);

    auto ts_xxyz_yyyyyy = pbuffer.data(idx_ovl_gi + 133);

    auto ts_xxyz_yyyyyz = pbuffer.data(idx_ovl_gi + 134);

    auto ts_xxyz_yyyyzz = pbuffer.data(idx_ovl_gi + 135);

    auto ts_xxyz_yyyzzz = pbuffer.data(idx_ovl_gi + 136);

    auto ts_xxyz_yyzzzz = pbuffer.data(idx_ovl_gi + 137);

    auto ts_xxyz_yzzzzz = pbuffer.data(idx_ovl_gi + 138);

    auto ts_xxyz_zzzzzz = pbuffer.data(idx_ovl_gi + 139);

    auto ts_xxzz_xxxxxx = pbuffer.data(idx_ovl_gi + 140);

    auto ts_xxzz_xxxxxy = pbuffer.data(idx_ovl_gi + 141);

    auto ts_xxzz_xxxxxz = pbuffer.data(idx_ovl_gi + 142);

    auto ts_xxzz_xxxxyy = pbuffer.data(idx_ovl_gi + 143);

    auto ts_xxzz_xxxxyz = pbuffer.data(idx_ovl_gi + 144);

    auto ts_xxzz_xxxxzz = pbuffer.data(idx_ovl_gi + 145);

    auto ts_xxzz_xxxyyy = pbuffer.data(idx_ovl_gi + 146);

    auto ts_xxzz_xxxyyz = pbuffer.data(idx_ovl_gi + 147);

    auto ts_xxzz_xxxyzz = pbuffer.data(idx_ovl_gi + 148);

    auto ts_xxzz_xxxzzz = pbuffer.data(idx_ovl_gi + 149);

    auto ts_xxzz_xxyyyy = pbuffer.data(idx_ovl_gi + 150);

    auto ts_xxzz_xxyyyz = pbuffer.data(idx_ovl_gi + 151);

    auto ts_xxzz_xxyyzz = pbuffer.data(idx_ovl_gi + 152);

    auto ts_xxzz_xxyzzz = pbuffer.data(idx_ovl_gi + 153);

    auto ts_xxzz_xxzzzz = pbuffer.data(idx_ovl_gi + 154);

    auto ts_xxzz_xyyyyy = pbuffer.data(idx_ovl_gi + 155);

    auto ts_xxzz_xyyyyz = pbuffer.data(idx_ovl_gi + 156);

    auto ts_xxzz_xyyyzz = pbuffer.data(idx_ovl_gi + 157);

    auto ts_xxzz_xyyzzz = pbuffer.data(idx_ovl_gi + 158);

    auto ts_xxzz_xyzzzz = pbuffer.data(idx_ovl_gi + 159);

    auto ts_xxzz_xzzzzz = pbuffer.data(idx_ovl_gi + 160);

    auto ts_xxzz_yyyyyy = pbuffer.data(idx_ovl_gi + 161);

    auto ts_xxzz_yyyyyz = pbuffer.data(idx_ovl_gi + 162);

    auto ts_xxzz_yyyyzz = pbuffer.data(idx_ovl_gi + 163);

    auto ts_xxzz_yyyzzz = pbuffer.data(idx_ovl_gi + 164);

    auto ts_xxzz_yyzzzz = pbuffer.data(idx_ovl_gi + 165);

    auto ts_xxzz_yzzzzz = pbuffer.data(idx_ovl_gi + 166);

    auto ts_xxzz_zzzzzz = pbuffer.data(idx_ovl_gi + 167);

    auto ts_xyyy_xxxxxx = pbuffer.data(idx_ovl_gi + 168);

    auto ts_xyyy_xxxxxy = pbuffer.data(idx_ovl_gi + 169);

    auto ts_xyyy_xxxxxz = pbuffer.data(idx_ovl_gi + 170);

    auto ts_xyyy_xxxxyy = pbuffer.data(idx_ovl_gi + 171);

    auto ts_xyyy_xxxxyz = pbuffer.data(idx_ovl_gi + 172);

    auto ts_xyyy_xxxxzz = pbuffer.data(idx_ovl_gi + 173);

    auto ts_xyyy_xxxyyy = pbuffer.data(idx_ovl_gi + 174);

    auto ts_xyyy_xxxyyz = pbuffer.data(idx_ovl_gi + 175);

    auto ts_xyyy_xxxyzz = pbuffer.data(idx_ovl_gi + 176);

    auto ts_xyyy_xxxzzz = pbuffer.data(idx_ovl_gi + 177);

    auto ts_xyyy_xxyyyy = pbuffer.data(idx_ovl_gi + 178);

    auto ts_xyyy_xxyyyz = pbuffer.data(idx_ovl_gi + 179);

    auto ts_xyyy_xxyyzz = pbuffer.data(idx_ovl_gi + 180);

    auto ts_xyyy_xxyzzz = pbuffer.data(idx_ovl_gi + 181);

    auto ts_xyyy_xxzzzz = pbuffer.data(idx_ovl_gi + 182);

    auto ts_xyyy_xyyyyy = pbuffer.data(idx_ovl_gi + 183);

    auto ts_xyyy_xyyyyz = pbuffer.data(idx_ovl_gi + 184);

    auto ts_xyyy_xyyyzz = pbuffer.data(idx_ovl_gi + 185);

    auto ts_xyyy_xyyzzz = pbuffer.data(idx_ovl_gi + 186);

    auto ts_xyyy_xyzzzz = pbuffer.data(idx_ovl_gi + 187);

    auto ts_xyyy_xzzzzz = pbuffer.data(idx_ovl_gi + 188);

    auto ts_xyyy_yyyyyy = pbuffer.data(idx_ovl_gi + 189);

    auto ts_xyyy_yyyyyz = pbuffer.data(idx_ovl_gi + 190);

    auto ts_xyyy_yyyyzz = pbuffer.data(idx_ovl_gi + 191);

    auto ts_xyyy_yyyzzz = pbuffer.data(idx_ovl_gi + 192);

    auto ts_xyyy_yyzzzz = pbuffer.data(idx_ovl_gi + 193);

    auto ts_xyyy_yzzzzz = pbuffer.data(idx_ovl_gi + 194);

    auto ts_xyyy_zzzzzz = pbuffer.data(idx_ovl_gi + 195);

    auto ts_xyyz_xxxxxx = pbuffer.data(idx_ovl_gi + 196);

    auto ts_xyyz_xxxxxy = pbuffer.data(idx_ovl_gi + 197);

    auto ts_xyyz_xxxxxz = pbuffer.data(idx_ovl_gi + 198);

    auto ts_xyyz_xxxxyy = pbuffer.data(idx_ovl_gi + 199);

    auto ts_xyyz_xxxxyz = pbuffer.data(idx_ovl_gi + 200);

    auto ts_xyyz_xxxxzz = pbuffer.data(idx_ovl_gi + 201);

    auto ts_xyyz_xxxyyy = pbuffer.data(idx_ovl_gi + 202);

    auto ts_xyyz_xxxyyz = pbuffer.data(idx_ovl_gi + 203);

    auto ts_xyyz_xxxyzz = pbuffer.data(idx_ovl_gi + 204);

    auto ts_xyyz_xxxzzz = pbuffer.data(idx_ovl_gi + 205);

    auto ts_xyyz_xxyyyy = pbuffer.data(idx_ovl_gi + 206);

    auto ts_xyyz_xxyyyz = pbuffer.data(idx_ovl_gi + 207);

    auto ts_xyyz_xxyyzz = pbuffer.data(idx_ovl_gi + 208);

    auto ts_xyyz_xxyzzz = pbuffer.data(idx_ovl_gi + 209);

    auto ts_xyyz_xxzzzz = pbuffer.data(idx_ovl_gi + 210);

    auto ts_xyyz_xyyyyy = pbuffer.data(idx_ovl_gi + 211);

    auto ts_xyyz_xyyyyz = pbuffer.data(idx_ovl_gi + 212);

    auto ts_xyyz_xyyyzz = pbuffer.data(idx_ovl_gi + 213);

    auto ts_xyyz_xyyzzz = pbuffer.data(idx_ovl_gi + 214);

    auto ts_xyyz_xyzzzz = pbuffer.data(idx_ovl_gi + 215);

    auto ts_xyyz_xzzzzz = pbuffer.data(idx_ovl_gi + 216);

    auto ts_xyyz_yyyyyy = pbuffer.data(idx_ovl_gi + 217);

    auto ts_xyyz_yyyyyz = pbuffer.data(idx_ovl_gi + 218);

    auto ts_xyyz_yyyyzz = pbuffer.data(idx_ovl_gi + 219);

    auto ts_xyyz_yyyzzz = pbuffer.data(idx_ovl_gi + 220);

    auto ts_xyyz_yyzzzz = pbuffer.data(idx_ovl_gi + 221);

    auto ts_xyyz_yzzzzz = pbuffer.data(idx_ovl_gi + 222);

    auto ts_xyyz_zzzzzz = pbuffer.data(idx_ovl_gi + 223);

    auto ts_xyzz_xxxxxx = pbuffer.data(idx_ovl_gi + 224);

    auto ts_xyzz_xxxxxy = pbuffer.data(idx_ovl_gi + 225);

    auto ts_xyzz_xxxxxz = pbuffer.data(idx_ovl_gi + 226);

    auto ts_xyzz_xxxxyy = pbuffer.data(idx_ovl_gi + 227);

    auto ts_xyzz_xxxxyz = pbuffer.data(idx_ovl_gi + 228);

    auto ts_xyzz_xxxxzz = pbuffer.data(idx_ovl_gi + 229);

    auto ts_xyzz_xxxyyy = pbuffer.data(idx_ovl_gi + 230);

    auto ts_xyzz_xxxyyz = pbuffer.data(idx_ovl_gi + 231);

    auto ts_xyzz_xxxyzz = pbuffer.data(idx_ovl_gi + 232);

    auto ts_xyzz_xxxzzz = pbuffer.data(idx_ovl_gi + 233);

    auto ts_xyzz_xxyyyy = pbuffer.data(idx_ovl_gi + 234);

    auto ts_xyzz_xxyyyz = pbuffer.data(idx_ovl_gi + 235);

    auto ts_xyzz_xxyyzz = pbuffer.data(idx_ovl_gi + 236);

    auto ts_xyzz_xxyzzz = pbuffer.data(idx_ovl_gi + 237);

    auto ts_xyzz_xxzzzz = pbuffer.data(idx_ovl_gi + 238);

    auto ts_xyzz_xyyyyy = pbuffer.data(idx_ovl_gi + 239);

    auto ts_xyzz_xyyyyz = pbuffer.data(idx_ovl_gi + 240);

    auto ts_xyzz_xyyyzz = pbuffer.data(idx_ovl_gi + 241);

    auto ts_xyzz_xyyzzz = pbuffer.data(idx_ovl_gi + 242);

    auto ts_xyzz_xyzzzz = pbuffer.data(idx_ovl_gi + 243);

    auto ts_xyzz_xzzzzz = pbuffer.data(idx_ovl_gi + 244);

    auto ts_xyzz_yyyyyy = pbuffer.data(idx_ovl_gi + 245);

    auto ts_xyzz_yyyyyz = pbuffer.data(idx_ovl_gi + 246);

    auto ts_xyzz_yyyyzz = pbuffer.data(idx_ovl_gi + 247);

    auto ts_xyzz_yyyzzz = pbuffer.data(idx_ovl_gi + 248);

    auto ts_xyzz_yyzzzz = pbuffer.data(idx_ovl_gi + 249);

    auto ts_xyzz_yzzzzz = pbuffer.data(idx_ovl_gi + 250);

    auto ts_xyzz_zzzzzz = pbuffer.data(idx_ovl_gi + 251);

    auto ts_xzzz_xxxxxx = pbuffer.data(idx_ovl_gi + 252);

    auto ts_xzzz_xxxxxy = pbuffer.data(idx_ovl_gi + 253);

    auto ts_xzzz_xxxxxz = pbuffer.data(idx_ovl_gi + 254);

    auto ts_xzzz_xxxxyy = pbuffer.data(idx_ovl_gi + 255);

    auto ts_xzzz_xxxxyz = pbuffer.data(idx_ovl_gi + 256);

    auto ts_xzzz_xxxxzz = pbuffer.data(idx_ovl_gi + 257);

    auto ts_xzzz_xxxyyy = pbuffer.data(idx_ovl_gi + 258);

    auto ts_xzzz_xxxyyz = pbuffer.data(idx_ovl_gi + 259);

    auto ts_xzzz_xxxyzz = pbuffer.data(idx_ovl_gi + 260);

    auto ts_xzzz_xxxzzz = pbuffer.data(idx_ovl_gi + 261);

    auto ts_xzzz_xxyyyy = pbuffer.data(idx_ovl_gi + 262);

    auto ts_xzzz_xxyyyz = pbuffer.data(idx_ovl_gi + 263);

    auto ts_xzzz_xxyyzz = pbuffer.data(idx_ovl_gi + 264);

    auto ts_xzzz_xxyzzz = pbuffer.data(idx_ovl_gi + 265);

    auto ts_xzzz_xxzzzz = pbuffer.data(idx_ovl_gi + 266);

    auto ts_xzzz_xyyyyy = pbuffer.data(idx_ovl_gi + 267);

    auto ts_xzzz_xyyyyz = pbuffer.data(idx_ovl_gi + 268);

    auto ts_xzzz_xyyyzz = pbuffer.data(idx_ovl_gi + 269);

    auto ts_xzzz_xyyzzz = pbuffer.data(idx_ovl_gi + 270);

    auto ts_xzzz_xyzzzz = pbuffer.data(idx_ovl_gi + 271);

    auto ts_xzzz_xzzzzz = pbuffer.data(idx_ovl_gi + 272);

    auto ts_xzzz_yyyyyy = pbuffer.data(idx_ovl_gi + 273);

    auto ts_xzzz_yyyyyz = pbuffer.data(idx_ovl_gi + 274);

    auto ts_xzzz_yyyyzz = pbuffer.data(idx_ovl_gi + 275);

    auto ts_xzzz_yyyzzz = pbuffer.data(idx_ovl_gi + 276);

    auto ts_xzzz_yyzzzz = pbuffer.data(idx_ovl_gi + 277);

    auto ts_xzzz_yzzzzz = pbuffer.data(idx_ovl_gi + 278);

    auto ts_xzzz_zzzzzz = pbuffer.data(idx_ovl_gi + 279);

    auto ts_yyyy_xxxxxx = pbuffer.data(idx_ovl_gi + 280);

    auto ts_yyyy_xxxxxy = pbuffer.data(idx_ovl_gi + 281);

    auto ts_yyyy_xxxxxz = pbuffer.data(idx_ovl_gi + 282);

    auto ts_yyyy_xxxxyy = pbuffer.data(idx_ovl_gi + 283);

    auto ts_yyyy_xxxxyz = pbuffer.data(idx_ovl_gi + 284);

    auto ts_yyyy_xxxxzz = pbuffer.data(idx_ovl_gi + 285);

    auto ts_yyyy_xxxyyy = pbuffer.data(idx_ovl_gi + 286);

    auto ts_yyyy_xxxyyz = pbuffer.data(idx_ovl_gi + 287);

    auto ts_yyyy_xxxyzz = pbuffer.data(idx_ovl_gi + 288);

    auto ts_yyyy_xxxzzz = pbuffer.data(idx_ovl_gi + 289);

    auto ts_yyyy_xxyyyy = pbuffer.data(idx_ovl_gi + 290);

    auto ts_yyyy_xxyyyz = pbuffer.data(idx_ovl_gi + 291);

    auto ts_yyyy_xxyyzz = pbuffer.data(idx_ovl_gi + 292);

    auto ts_yyyy_xxyzzz = pbuffer.data(idx_ovl_gi + 293);

    auto ts_yyyy_xxzzzz = pbuffer.data(idx_ovl_gi + 294);

    auto ts_yyyy_xyyyyy = pbuffer.data(idx_ovl_gi + 295);

    auto ts_yyyy_xyyyyz = pbuffer.data(idx_ovl_gi + 296);

    auto ts_yyyy_xyyyzz = pbuffer.data(idx_ovl_gi + 297);

    auto ts_yyyy_xyyzzz = pbuffer.data(idx_ovl_gi + 298);

    auto ts_yyyy_xyzzzz = pbuffer.data(idx_ovl_gi + 299);

    auto ts_yyyy_xzzzzz = pbuffer.data(idx_ovl_gi + 300);

    auto ts_yyyy_yyyyyy = pbuffer.data(idx_ovl_gi + 301);

    auto ts_yyyy_yyyyyz = pbuffer.data(idx_ovl_gi + 302);

    auto ts_yyyy_yyyyzz = pbuffer.data(idx_ovl_gi + 303);

    auto ts_yyyy_yyyzzz = pbuffer.data(idx_ovl_gi + 304);

    auto ts_yyyy_yyzzzz = pbuffer.data(idx_ovl_gi + 305);

    auto ts_yyyy_yzzzzz = pbuffer.data(idx_ovl_gi + 306);

    auto ts_yyyy_zzzzzz = pbuffer.data(idx_ovl_gi + 307);

    auto ts_yyyz_xxxxxx = pbuffer.data(idx_ovl_gi + 308);

    auto ts_yyyz_xxxxxy = pbuffer.data(idx_ovl_gi + 309);

    auto ts_yyyz_xxxxxz = pbuffer.data(idx_ovl_gi + 310);

    auto ts_yyyz_xxxxyy = pbuffer.data(idx_ovl_gi + 311);

    auto ts_yyyz_xxxxyz = pbuffer.data(idx_ovl_gi + 312);

    auto ts_yyyz_xxxxzz = pbuffer.data(idx_ovl_gi + 313);

    auto ts_yyyz_xxxyyy = pbuffer.data(idx_ovl_gi + 314);

    auto ts_yyyz_xxxyyz = pbuffer.data(idx_ovl_gi + 315);

    auto ts_yyyz_xxxyzz = pbuffer.data(idx_ovl_gi + 316);

    auto ts_yyyz_xxxzzz = pbuffer.data(idx_ovl_gi + 317);

    auto ts_yyyz_xxyyyy = pbuffer.data(idx_ovl_gi + 318);

    auto ts_yyyz_xxyyyz = pbuffer.data(idx_ovl_gi + 319);

    auto ts_yyyz_xxyyzz = pbuffer.data(idx_ovl_gi + 320);

    auto ts_yyyz_xxyzzz = pbuffer.data(idx_ovl_gi + 321);

    auto ts_yyyz_xxzzzz = pbuffer.data(idx_ovl_gi + 322);

    auto ts_yyyz_xyyyyy = pbuffer.data(idx_ovl_gi + 323);

    auto ts_yyyz_xyyyyz = pbuffer.data(idx_ovl_gi + 324);

    auto ts_yyyz_xyyyzz = pbuffer.data(idx_ovl_gi + 325);

    auto ts_yyyz_xyyzzz = pbuffer.data(idx_ovl_gi + 326);

    auto ts_yyyz_xyzzzz = pbuffer.data(idx_ovl_gi + 327);

    auto ts_yyyz_xzzzzz = pbuffer.data(idx_ovl_gi + 328);

    auto ts_yyyz_yyyyyy = pbuffer.data(idx_ovl_gi + 329);

    auto ts_yyyz_yyyyyz = pbuffer.data(idx_ovl_gi + 330);

    auto ts_yyyz_yyyyzz = pbuffer.data(idx_ovl_gi + 331);

    auto ts_yyyz_yyyzzz = pbuffer.data(idx_ovl_gi + 332);

    auto ts_yyyz_yyzzzz = pbuffer.data(idx_ovl_gi + 333);

    auto ts_yyyz_yzzzzz = pbuffer.data(idx_ovl_gi + 334);

    auto ts_yyyz_zzzzzz = pbuffer.data(idx_ovl_gi + 335);

    auto ts_yyzz_xxxxxx = pbuffer.data(idx_ovl_gi + 336);

    auto ts_yyzz_xxxxxy = pbuffer.data(idx_ovl_gi + 337);

    auto ts_yyzz_xxxxxz = pbuffer.data(idx_ovl_gi + 338);

    auto ts_yyzz_xxxxyy = pbuffer.data(idx_ovl_gi + 339);

    auto ts_yyzz_xxxxyz = pbuffer.data(idx_ovl_gi + 340);

    auto ts_yyzz_xxxxzz = pbuffer.data(idx_ovl_gi + 341);

    auto ts_yyzz_xxxyyy = pbuffer.data(idx_ovl_gi + 342);

    auto ts_yyzz_xxxyyz = pbuffer.data(idx_ovl_gi + 343);

    auto ts_yyzz_xxxyzz = pbuffer.data(idx_ovl_gi + 344);

    auto ts_yyzz_xxxzzz = pbuffer.data(idx_ovl_gi + 345);

    auto ts_yyzz_xxyyyy = pbuffer.data(idx_ovl_gi + 346);

    auto ts_yyzz_xxyyyz = pbuffer.data(idx_ovl_gi + 347);

    auto ts_yyzz_xxyyzz = pbuffer.data(idx_ovl_gi + 348);

    auto ts_yyzz_xxyzzz = pbuffer.data(idx_ovl_gi + 349);

    auto ts_yyzz_xxzzzz = pbuffer.data(idx_ovl_gi + 350);

    auto ts_yyzz_xyyyyy = pbuffer.data(idx_ovl_gi + 351);

    auto ts_yyzz_xyyyyz = pbuffer.data(idx_ovl_gi + 352);

    auto ts_yyzz_xyyyzz = pbuffer.data(idx_ovl_gi + 353);

    auto ts_yyzz_xyyzzz = pbuffer.data(idx_ovl_gi + 354);

    auto ts_yyzz_xyzzzz = pbuffer.data(idx_ovl_gi + 355);

    auto ts_yyzz_xzzzzz = pbuffer.data(idx_ovl_gi + 356);

    auto ts_yyzz_yyyyyy = pbuffer.data(idx_ovl_gi + 357);

    auto ts_yyzz_yyyyyz = pbuffer.data(idx_ovl_gi + 358);

    auto ts_yyzz_yyyyzz = pbuffer.data(idx_ovl_gi + 359);

    auto ts_yyzz_yyyzzz = pbuffer.data(idx_ovl_gi + 360);

    auto ts_yyzz_yyzzzz = pbuffer.data(idx_ovl_gi + 361);

    auto ts_yyzz_yzzzzz = pbuffer.data(idx_ovl_gi + 362);

    auto ts_yyzz_zzzzzz = pbuffer.data(idx_ovl_gi + 363);

    auto ts_yzzz_xxxxxx = pbuffer.data(idx_ovl_gi + 364);

    auto ts_yzzz_xxxxxy = pbuffer.data(idx_ovl_gi + 365);

    auto ts_yzzz_xxxxxz = pbuffer.data(idx_ovl_gi + 366);

    auto ts_yzzz_xxxxyy = pbuffer.data(idx_ovl_gi + 367);

    auto ts_yzzz_xxxxyz = pbuffer.data(idx_ovl_gi + 368);

    auto ts_yzzz_xxxxzz = pbuffer.data(idx_ovl_gi + 369);

    auto ts_yzzz_xxxyyy = pbuffer.data(idx_ovl_gi + 370);

    auto ts_yzzz_xxxyyz = pbuffer.data(idx_ovl_gi + 371);

    auto ts_yzzz_xxxyzz = pbuffer.data(idx_ovl_gi + 372);

    auto ts_yzzz_xxxzzz = pbuffer.data(idx_ovl_gi + 373);

    auto ts_yzzz_xxyyyy = pbuffer.data(idx_ovl_gi + 374);

    auto ts_yzzz_xxyyyz = pbuffer.data(idx_ovl_gi + 375);

    auto ts_yzzz_xxyyzz = pbuffer.data(idx_ovl_gi + 376);

    auto ts_yzzz_xxyzzz = pbuffer.data(idx_ovl_gi + 377);

    auto ts_yzzz_xxzzzz = pbuffer.data(idx_ovl_gi + 378);

    auto ts_yzzz_xyyyyy = pbuffer.data(idx_ovl_gi + 379);

    auto ts_yzzz_xyyyyz = pbuffer.data(idx_ovl_gi + 380);

    auto ts_yzzz_xyyyzz = pbuffer.data(idx_ovl_gi + 381);

    auto ts_yzzz_xyyzzz = pbuffer.data(idx_ovl_gi + 382);

    auto ts_yzzz_xyzzzz = pbuffer.data(idx_ovl_gi + 383);

    auto ts_yzzz_xzzzzz = pbuffer.data(idx_ovl_gi + 384);

    auto ts_yzzz_yyyyyy = pbuffer.data(idx_ovl_gi + 385);

    auto ts_yzzz_yyyyyz = pbuffer.data(idx_ovl_gi + 386);

    auto ts_yzzz_yyyyzz = pbuffer.data(idx_ovl_gi + 387);

    auto ts_yzzz_yyyzzz = pbuffer.data(idx_ovl_gi + 388);

    auto ts_yzzz_yyzzzz = pbuffer.data(idx_ovl_gi + 389);

    auto ts_yzzz_yzzzzz = pbuffer.data(idx_ovl_gi + 390);

    auto ts_yzzz_zzzzzz = pbuffer.data(idx_ovl_gi + 391);

    auto ts_zzzz_xxxxxx = pbuffer.data(idx_ovl_gi + 392);

    auto ts_zzzz_xxxxxy = pbuffer.data(idx_ovl_gi + 393);

    auto ts_zzzz_xxxxxz = pbuffer.data(idx_ovl_gi + 394);

    auto ts_zzzz_xxxxyy = pbuffer.data(idx_ovl_gi + 395);

    auto ts_zzzz_xxxxyz = pbuffer.data(idx_ovl_gi + 396);

    auto ts_zzzz_xxxxzz = pbuffer.data(idx_ovl_gi + 397);

    auto ts_zzzz_xxxyyy = pbuffer.data(idx_ovl_gi + 398);

    auto ts_zzzz_xxxyyz = pbuffer.data(idx_ovl_gi + 399);

    auto ts_zzzz_xxxyzz = pbuffer.data(idx_ovl_gi + 400);

    auto ts_zzzz_xxxzzz = pbuffer.data(idx_ovl_gi + 401);

    auto ts_zzzz_xxyyyy = pbuffer.data(idx_ovl_gi + 402);

    auto ts_zzzz_xxyyyz = pbuffer.data(idx_ovl_gi + 403);

    auto ts_zzzz_xxyyzz = pbuffer.data(idx_ovl_gi + 404);

    auto ts_zzzz_xxyzzz = pbuffer.data(idx_ovl_gi + 405);

    auto ts_zzzz_xxzzzz = pbuffer.data(idx_ovl_gi + 406);

    auto ts_zzzz_xyyyyy = pbuffer.data(idx_ovl_gi + 407);

    auto ts_zzzz_xyyyyz = pbuffer.data(idx_ovl_gi + 408);

    auto ts_zzzz_xyyyzz = pbuffer.data(idx_ovl_gi + 409);

    auto ts_zzzz_xyyzzz = pbuffer.data(idx_ovl_gi + 410);

    auto ts_zzzz_xyzzzz = pbuffer.data(idx_ovl_gi + 411);

    auto ts_zzzz_xzzzzz = pbuffer.data(idx_ovl_gi + 412);

    auto ts_zzzz_yyyyyy = pbuffer.data(idx_ovl_gi + 413);

    auto ts_zzzz_yyyyyz = pbuffer.data(idx_ovl_gi + 414);

    auto ts_zzzz_yyyyzz = pbuffer.data(idx_ovl_gi + 415);

    auto ts_zzzz_yyyzzz = pbuffer.data(idx_ovl_gi + 416);

    auto ts_zzzz_yyzzzz = pbuffer.data(idx_ovl_gi + 417);

    auto ts_zzzz_yzzzzz = pbuffer.data(idx_ovl_gi + 418);

    auto ts_zzzz_zzzzzz = pbuffer.data(idx_ovl_gi + 419);

    // Set up 0-28 components of targeted buffer : GI

    auto tk_xxxx_xxxxxx = pbuffer.data(idx_kin_gi);

    auto tk_xxxx_xxxxxy = pbuffer.data(idx_kin_gi + 1);

    auto tk_xxxx_xxxxxz = pbuffer.data(idx_kin_gi + 2);

    auto tk_xxxx_xxxxyy = pbuffer.data(idx_kin_gi + 3);

    auto tk_xxxx_xxxxyz = pbuffer.data(idx_kin_gi + 4);

    auto tk_xxxx_xxxxzz = pbuffer.data(idx_kin_gi + 5);

    auto tk_xxxx_xxxyyy = pbuffer.data(idx_kin_gi + 6);

    auto tk_xxxx_xxxyyz = pbuffer.data(idx_kin_gi + 7);

    auto tk_xxxx_xxxyzz = pbuffer.data(idx_kin_gi + 8);

    auto tk_xxxx_xxxzzz = pbuffer.data(idx_kin_gi + 9);

    auto tk_xxxx_xxyyyy = pbuffer.data(idx_kin_gi + 10);

    auto tk_xxxx_xxyyyz = pbuffer.data(idx_kin_gi + 11);

    auto tk_xxxx_xxyyzz = pbuffer.data(idx_kin_gi + 12);

    auto tk_xxxx_xxyzzz = pbuffer.data(idx_kin_gi + 13);

    auto tk_xxxx_xxzzzz = pbuffer.data(idx_kin_gi + 14);

    auto tk_xxxx_xyyyyy = pbuffer.data(idx_kin_gi + 15);

    auto tk_xxxx_xyyyyz = pbuffer.data(idx_kin_gi + 16);

    auto tk_xxxx_xyyyzz = pbuffer.data(idx_kin_gi + 17);

    auto tk_xxxx_xyyzzz = pbuffer.data(idx_kin_gi + 18);

    auto tk_xxxx_xyzzzz = pbuffer.data(idx_kin_gi + 19);

    auto tk_xxxx_xzzzzz = pbuffer.data(idx_kin_gi + 20);

    auto tk_xxxx_yyyyyy = pbuffer.data(idx_kin_gi + 21);

    auto tk_xxxx_yyyyyz = pbuffer.data(idx_kin_gi + 22);

    auto tk_xxxx_yyyyzz = pbuffer.data(idx_kin_gi + 23);

    auto tk_xxxx_yyyzzz = pbuffer.data(idx_kin_gi + 24);

    auto tk_xxxx_yyzzzz = pbuffer.data(idx_kin_gi + 25);

    auto tk_xxxx_yzzzzz = pbuffer.data(idx_kin_gi + 26);

    auto tk_xxxx_zzzzzz = pbuffer.data(idx_kin_gi + 27);

#pragma omp simd aligned(pa_x,               \
                             tk_xx_xxxxxx,   \
                             tk_xx_xxxxxy,   \
                             tk_xx_xxxxxz,   \
                             tk_xx_xxxxyy,   \
                             tk_xx_xxxxyz,   \
                             tk_xx_xxxxzz,   \
                             tk_xx_xxxyyy,   \
                             tk_xx_xxxyyz,   \
                             tk_xx_xxxyzz,   \
                             tk_xx_xxxzzz,   \
                             tk_xx_xxyyyy,   \
                             tk_xx_xxyyyz,   \
                             tk_xx_xxyyzz,   \
                             tk_xx_xxyzzz,   \
                             tk_xx_xxzzzz,   \
                             tk_xx_xyyyyy,   \
                             tk_xx_xyyyyz,   \
                             tk_xx_xyyyzz,   \
                             tk_xx_xyyzzz,   \
                             tk_xx_xyzzzz,   \
                             tk_xx_xzzzzz,   \
                             tk_xx_yyyyyy,   \
                             tk_xx_yyyyyz,   \
                             tk_xx_yyyyzz,   \
                             tk_xx_yyyzzz,   \
                             tk_xx_yyzzzz,   \
                             tk_xx_yzzzzz,   \
                             tk_xx_zzzzzz,   \
                             tk_xxx_xxxxx,   \
                             tk_xxx_xxxxxx,  \
                             tk_xxx_xxxxxy,  \
                             tk_xxx_xxxxxz,  \
                             tk_xxx_xxxxy,   \
                             tk_xxx_xxxxyy,  \
                             tk_xxx_xxxxyz,  \
                             tk_xxx_xxxxz,   \
                             tk_xxx_xxxxzz,  \
                             tk_xxx_xxxyy,   \
                             tk_xxx_xxxyyy,  \
                             tk_xxx_xxxyyz,  \
                             tk_xxx_xxxyz,   \
                             tk_xxx_xxxyzz,  \
                             tk_xxx_xxxzz,   \
                             tk_xxx_xxxzzz,  \
                             tk_xxx_xxyyy,   \
                             tk_xxx_xxyyyy,  \
                             tk_xxx_xxyyyz,  \
                             tk_xxx_xxyyz,   \
                             tk_xxx_xxyyzz,  \
                             tk_xxx_xxyzz,   \
                             tk_xxx_xxyzzz,  \
                             tk_xxx_xxzzz,   \
                             tk_xxx_xxzzzz,  \
                             tk_xxx_xyyyy,   \
                             tk_xxx_xyyyyy,  \
                             tk_xxx_xyyyyz,  \
                             tk_xxx_xyyyz,   \
                             tk_xxx_xyyyzz,  \
                             tk_xxx_xyyzz,   \
                             tk_xxx_xyyzzz,  \
                             tk_xxx_xyzzz,   \
                             tk_xxx_xyzzzz,  \
                             tk_xxx_xzzzz,   \
                             tk_xxx_xzzzzz,  \
                             tk_xxx_yyyyy,   \
                             tk_xxx_yyyyyy,  \
                             tk_xxx_yyyyyz,  \
                             tk_xxx_yyyyz,   \
                             tk_xxx_yyyyzz,  \
                             tk_xxx_yyyzz,   \
                             tk_xxx_yyyzzz,  \
                             tk_xxx_yyzzz,   \
                             tk_xxx_yyzzzz,  \
                             tk_xxx_yzzzz,   \
                             tk_xxx_yzzzzz,  \
                             tk_xxx_zzzzz,   \
                             tk_xxx_zzzzzz,  \
                             tk_xxxx_xxxxxx, \
                             tk_xxxx_xxxxxy, \
                             tk_xxxx_xxxxxz, \
                             tk_xxxx_xxxxyy, \
                             tk_xxxx_xxxxyz, \
                             tk_xxxx_xxxxzz, \
                             tk_xxxx_xxxyyy, \
                             tk_xxxx_xxxyyz, \
                             tk_xxxx_xxxyzz, \
                             tk_xxxx_xxxzzz, \
                             tk_xxxx_xxyyyy, \
                             tk_xxxx_xxyyyz, \
                             tk_xxxx_xxyyzz, \
                             tk_xxxx_xxyzzz, \
                             tk_xxxx_xxzzzz, \
                             tk_xxxx_xyyyyy, \
                             tk_xxxx_xyyyyz, \
                             tk_xxxx_xyyyzz, \
                             tk_xxxx_xyyzzz, \
                             tk_xxxx_xyzzzz, \
                             tk_xxxx_xzzzzz, \
                             tk_xxxx_yyyyyy, \
                             tk_xxxx_yyyyyz, \
                             tk_xxxx_yyyyzz, \
                             tk_xxxx_yyyzzz, \
                             tk_xxxx_yyzzzz, \
                             tk_xxxx_yzzzzz, \
                             tk_xxxx_zzzzzz, \
                             ts_xx_xxxxxx,   \
                             ts_xx_xxxxxy,   \
                             ts_xx_xxxxxz,   \
                             ts_xx_xxxxyy,   \
                             ts_xx_xxxxyz,   \
                             ts_xx_xxxxzz,   \
                             ts_xx_xxxyyy,   \
                             ts_xx_xxxyyz,   \
                             ts_xx_xxxyzz,   \
                             ts_xx_xxxzzz,   \
                             ts_xx_xxyyyy,   \
                             ts_xx_xxyyyz,   \
                             ts_xx_xxyyzz,   \
                             ts_xx_xxyzzz,   \
                             ts_xx_xxzzzz,   \
                             ts_xx_xyyyyy,   \
                             ts_xx_xyyyyz,   \
                             ts_xx_xyyyzz,   \
                             ts_xx_xyyzzz,   \
                             ts_xx_xyzzzz,   \
                             ts_xx_xzzzzz,   \
                             ts_xx_yyyyyy,   \
                             ts_xx_yyyyyz,   \
                             ts_xx_yyyyzz,   \
                             ts_xx_yyyzzz,   \
                             ts_xx_yyzzzz,   \
                             ts_xx_yzzzzz,   \
                             ts_xx_zzzzzz,   \
                             ts_xxxx_xxxxxx, \
                             ts_xxxx_xxxxxy, \
                             ts_xxxx_xxxxxz, \
                             ts_xxxx_xxxxyy, \
                             ts_xxxx_xxxxyz, \
                             ts_xxxx_xxxxzz, \
                             ts_xxxx_xxxyyy, \
                             ts_xxxx_xxxyyz, \
                             ts_xxxx_xxxyzz, \
                             ts_xxxx_xxxzzz, \
                             ts_xxxx_xxyyyy, \
                             ts_xxxx_xxyyyz, \
                             ts_xxxx_xxyyzz, \
                             ts_xxxx_xxyzzz, \
                             ts_xxxx_xxzzzz, \
                             ts_xxxx_xyyyyy, \
                             ts_xxxx_xyyyyz, \
                             ts_xxxx_xyyyzz, \
                             ts_xxxx_xyyzzz, \
                             ts_xxxx_xyzzzz, \
                             ts_xxxx_xzzzzz, \
                             ts_xxxx_yyyyyy, \
                             ts_xxxx_yyyyyz, \
                             ts_xxxx_yyyyzz, \
                             ts_xxxx_yyyzzz, \
                             ts_xxxx_yyzzzz, \
                             ts_xxxx_yzzzzz, \
                             ts_xxxx_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxx_xxxxxx[i] = -6.0 * ts_xx_xxxxxx[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxxxx[i] * fe_0 + 6.0 * tk_xxx_xxxxx[i] * fe_0 +
                            tk_xxx_xxxxxx[i] * pa_x[i] + 2.0 * ts_xxxx_xxxxxx[i] * fz_0;

        tk_xxxx_xxxxxy[i] = -6.0 * ts_xx_xxxxxy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxxxy[i] * fe_0 + 5.0 * tk_xxx_xxxxy[i] * fe_0 +
                            tk_xxx_xxxxxy[i] * pa_x[i] + 2.0 * ts_xxxx_xxxxxy[i] * fz_0;

        tk_xxxx_xxxxxz[i] = -6.0 * ts_xx_xxxxxz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxxxz[i] * fe_0 + 5.0 * tk_xxx_xxxxz[i] * fe_0 +
                            tk_xxx_xxxxxz[i] * pa_x[i] + 2.0 * ts_xxxx_xxxxxz[i] * fz_0;

        tk_xxxx_xxxxyy[i] = -6.0 * ts_xx_xxxxyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxxyy[i] * fe_0 + 4.0 * tk_xxx_xxxyy[i] * fe_0 +
                            tk_xxx_xxxxyy[i] * pa_x[i] + 2.0 * ts_xxxx_xxxxyy[i] * fz_0;

        tk_xxxx_xxxxyz[i] = -6.0 * ts_xx_xxxxyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxxyz[i] * fe_0 + 4.0 * tk_xxx_xxxyz[i] * fe_0 +
                            tk_xxx_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxxx_xxxxyz[i] * fz_0;

        tk_xxxx_xxxxzz[i] = -6.0 * ts_xx_xxxxzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxxzz[i] * fe_0 + 4.0 * tk_xxx_xxxzz[i] * fe_0 +
                            tk_xxx_xxxxzz[i] * pa_x[i] + 2.0 * ts_xxxx_xxxxzz[i] * fz_0;

        tk_xxxx_xxxyyy[i] = -6.0 * ts_xx_xxxyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxyyy[i] * fe_0 + 3.0 * tk_xxx_xxyyy[i] * fe_0 +
                            tk_xxx_xxxyyy[i] * pa_x[i] + 2.0 * ts_xxxx_xxxyyy[i] * fz_0;

        tk_xxxx_xxxyyz[i] = -6.0 * ts_xx_xxxyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxyyz[i] * fe_0 + 3.0 * tk_xxx_xxyyz[i] * fe_0 +
                            tk_xxx_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxxx_xxxyyz[i] * fz_0;

        tk_xxxx_xxxyzz[i] = -6.0 * ts_xx_xxxyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxyzz[i] * fe_0 + 3.0 * tk_xxx_xxyzz[i] * fe_0 +
                            tk_xxx_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxxx_xxxyzz[i] * fz_0;

        tk_xxxx_xxxzzz[i] = -6.0 * ts_xx_xxxzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxxzzz[i] * fe_0 + 3.0 * tk_xxx_xxzzz[i] * fe_0 +
                            tk_xxx_xxxzzz[i] * pa_x[i] + 2.0 * ts_xxxx_xxxzzz[i] * fz_0;

        tk_xxxx_xxyyyy[i] = -6.0 * ts_xx_xxyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxyyyy[i] * fe_0 + 2.0 * tk_xxx_xyyyy[i] * fe_0 +
                            tk_xxx_xxyyyy[i] * pa_x[i] + 2.0 * ts_xxxx_xxyyyy[i] * fz_0;

        tk_xxxx_xxyyyz[i] = -6.0 * ts_xx_xxyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxyyyz[i] * fe_0 + 2.0 * tk_xxx_xyyyz[i] * fe_0 +
                            tk_xxx_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxxx_xxyyyz[i] * fz_0;

        tk_xxxx_xxyyzz[i] = -6.0 * ts_xx_xxyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxyyzz[i] * fe_0 + 2.0 * tk_xxx_xyyzz[i] * fe_0 +
                            tk_xxx_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxxx_xxyyzz[i] * fz_0;

        tk_xxxx_xxyzzz[i] = -6.0 * ts_xx_xxyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxyzzz[i] * fe_0 + 2.0 * tk_xxx_xyzzz[i] * fe_0 +
                            tk_xxx_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxxx_xxyzzz[i] * fz_0;

        tk_xxxx_xxzzzz[i] = -6.0 * ts_xx_xxzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xxzzzz[i] * fe_0 + 2.0 * tk_xxx_xzzzz[i] * fe_0 +
                            tk_xxx_xxzzzz[i] * pa_x[i] + 2.0 * ts_xxxx_xxzzzz[i] * fz_0;

        tk_xxxx_xyyyyy[i] = -6.0 * ts_xx_xyyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyyyyy[i] * fe_0 + tk_xxx_yyyyy[i] * fe_0 +
                            tk_xxx_xyyyyy[i] * pa_x[i] + 2.0 * ts_xxxx_xyyyyy[i] * fz_0;

        tk_xxxx_xyyyyz[i] = -6.0 * ts_xx_xyyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyyyyz[i] * fe_0 + tk_xxx_yyyyz[i] * fe_0 +
                            tk_xxx_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxxx_xyyyyz[i] * fz_0;

        tk_xxxx_xyyyzz[i] = -6.0 * ts_xx_xyyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyyyzz[i] * fe_0 + tk_xxx_yyyzz[i] * fe_0 +
                            tk_xxx_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxxx_xyyyzz[i] * fz_0;

        tk_xxxx_xyyzzz[i] = -6.0 * ts_xx_xyyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyyzzz[i] * fe_0 + tk_xxx_yyzzz[i] * fe_0 +
                            tk_xxx_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxxx_xyyzzz[i] * fz_0;

        tk_xxxx_xyzzzz[i] = -6.0 * ts_xx_xyzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xyzzzz[i] * fe_0 + tk_xxx_yzzzz[i] * fe_0 +
                            tk_xxx_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxxx_xyzzzz[i] * fz_0;

        tk_xxxx_xzzzzz[i] = -6.0 * ts_xx_xzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xzzzzz[i] * fe_0 + tk_xxx_zzzzz[i] * fe_0 +
                            tk_xxx_xzzzzz[i] * pa_x[i] + 2.0 * ts_xxxx_xzzzzz[i] * fz_0;

        tk_xxxx_yyyyyy[i] =
            -6.0 * ts_xx_yyyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyyyyy[i] * fe_0 + tk_xxx_yyyyyy[i] * pa_x[i] + 2.0 * ts_xxxx_yyyyyy[i] * fz_0;

        tk_xxxx_yyyyyz[i] =
            -6.0 * ts_xx_yyyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyyyyz[i] * fe_0 + tk_xxx_yyyyyz[i] * pa_x[i] + 2.0 * ts_xxxx_yyyyyz[i] * fz_0;

        tk_xxxx_yyyyzz[i] =
            -6.0 * ts_xx_yyyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyyyzz[i] * fe_0 + tk_xxx_yyyyzz[i] * pa_x[i] + 2.0 * ts_xxxx_yyyyzz[i] * fz_0;

        tk_xxxx_yyyzzz[i] =
            -6.0 * ts_xx_yyyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyyzzz[i] * fe_0 + tk_xxx_yyyzzz[i] * pa_x[i] + 2.0 * ts_xxxx_yyyzzz[i] * fz_0;

        tk_xxxx_yyzzzz[i] =
            -6.0 * ts_xx_yyzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yyzzzz[i] * fe_0 + tk_xxx_yyzzzz[i] * pa_x[i] + 2.0 * ts_xxxx_yyzzzz[i] * fz_0;

        tk_xxxx_yzzzzz[i] =
            -6.0 * ts_xx_yzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yzzzzz[i] * fe_0 + tk_xxx_yzzzzz[i] * pa_x[i] + 2.0 * ts_xxxx_yzzzzz[i] * fz_0;

        tk_xxxx_zzzzzz[i] =
            -6.0 * ts_xx_zzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_zzzzzz[i] * fe_0 + tk_xxx_zzzzzz[i] * pa_x[i] + 2.0 * ts_xxxx_zzzzzz[i] * fz_0;
    }

    // Set up 28-56 components of targeted buffer : GI

    auto tk_xxxy_xxxxxx = pbuffer.data(idx_kin_gi + 28);

    auto tk_xxxy_xxxxxy = pbuffer.data(idx_kin_gi + 29);

    auto tk_xxxy_xxxxxz = pbuffer.data(idx_kin_gi + 30);

    auto tk_xxxy_xxxxyy = pbuffer.data(idx_kin_gi + 31);

    auto tk_xxxy_xxxxyz = pbuffer.data(idx_kin_gi + 32);

    auto tk_xxxy_xxxxzz = pbuffer.data(idx_kin_gi + 33);

    auto tk_xxxy_xxxyyy = pbuffer.data(idx_kin_gi + 34);

    auto tk_xxxy_xxxyyz = pbuffer.data(idx_kin_gi + 35);

    auto tk_xxxy_xxxyzz = pbuffer.data(idx_kin_gi + 36);

    auto tk_xxxy_xxxzzz = pbuffer.data(idx_kin_gi + 37);

    auto tk_xxxy_xxyyyy = pbuffer.data(idx_kin_gi + 38);

    auto tk_xxxy_xxyyyz = pbuffer.data(idx_kin_gi + 39);

    auto tk_xxxy_xxyyzz = pbuffer.data(idx_kin_gi + 40);

    auto tk_xxxy_xxyzzz = pbuffer.data(idx_kin_gi + 41);

    auto tk_xxxy_xxzzzz = pbuffer.data(idx_kin_gi + 42);

    auto tk_xxxy_xyyyyy = pbuffer.data(idx_kin_gi + 43);

    auto tk_xxxy_xyyyyz = pbuffer.data(idx_kin_gi + 44);

    auto tk_xxxy_xyyyzz = pbuffer.data(idx_kin_gi + 45);

    auto tk_xxxy_xyyzzz = pbuffer.data(idx_kin_gi + 46);

    auto tk_xxxy_xyzzzz = pbuffer.data(idx_kin_gi + 47);

    auto tk_xxxy_xzzzzz = pbuffer.data(idx_kin_gi + 48);

    auto tk_xxxy_yyyyyy = pbuffer.data(idx_kin_gi + 49);

    auto tk_xxxy_yyyyyz = pbuffer.data(idx_kin_gi + 50);

    auto tk_xxxy_yyyyzz = pbuffer.data(idx_kin_gi + 51);

    auto tk_xxxy_yyyzzz = pbuffer.data(idx_kin_gi + 52);

    auto tk_xxxy_yyzzzz = pbuffer.data(idx_kin_gi + 53);

    auto tk_xxxy_yzzzzz = pbuffer.data(idx_kin_gi + 54);

    auto tk_xxxy_zzzzzz = pbuffer.data(idx_kin_gi + 55);

#pragma omp simd aligned(pa_y,               \
                             tk_xxx_xxxxx,   \
                             tk_xxx_xxxxxx,  \
                             tk_xxx_xxxxxy,  \
                             tk_xxx_xxxxxz,  \
                             tk_xxx_xxxxy,   \
                             tk_xxx_xxxxyy,  \
                             tk_xxx_xxxxyz,  \
                             tk_xxx_xxxxz,   \
                             tk_xxx_xxxxzz,  \
                             tk_xxx_xxxyy,   \
                             tk_xxx_xxxyyy,  \
                             tk_xxx_xxxyyz,  \
                             tk_xxx_xxxyz,   \
                             tk_xxx_xxxyzz,  \
                             tk_xxx_xxxzz,   \
                             tk_xxx_xxxzzz,  \
                             tk_xxx_xxyyy,   \
                             tk_xxx_xxyyyy,  \
                             tk_xxx_xxyyyz,  \
                             tk_xxx_xxyyz,   \
                             tk_xxx_xxyyzz,  \
                             tk_xxx_xxyzz,   \
                             tk_xxx_xxyzzz,  \
                             tk_xxx_xxzzz,   \
                             tk_xxx_xxzzzz,  \
                             tk_xxx_xyyyy,   \
                             tk_xxx_xyyyyy,  \
                             tk_xxx_xyyyyz,  \
                             tk_xxx_xyyyz,   \
                             tk_xxx_xyyyzz,  \
                             tk_xxx_xyyzz,   \
                             tk_xxx_xyyzzz,  \
                             tk_xxx_xyzzz,   \
                             tk_xxx_xyzzzz,  \
                             tk_xxx_xzzzz,   \
                             tk_xxx_xzzzzz,  \
                             tk_xxx_yyyyy,   \
                             tk_xxx_yyyyyy,  \
                             tk_xxx_yyyyyz,  \
                             tk_xxx_yyyyz,   \
                             tk_xxx_yyyyzz,  \
                             tk_xxx_yyyzz,   \
                             tk_xxx_yyyzzz,  \
                             tk_xxx_yyzzz,   \
                             tk_xxx_yyzzzz,  \
                             tk_xxx_yzzzz,   \
                             tk_xxx_yzzzzz,  \
                             tk_xxx_zzzzz,   \
                             tk_xxx_zzzzzz,  \
                             tk_xxxy_xxxxxx, \
                             tk_xxxy_xxxxxy, \
                             tk_xxxy_xxxxxz, \
                             tk_xxxy_xxxxyy, \
                             tk_xxxy_xxxxyz, \
                             tk_xxxy_xxxxzz, \
                             tk_xxxy_xxxyyy, \
                             tk_xxxy_xxxyyz, \
                             tk_xxxy_xxxyzz, \
                             tk_xxxy_xxxzzz, \
                             tk_xxxy_xxyyyy, \
                             tk_xxxy_xxyyyz, \
                             tk_xxxy_xxyyzz, \
                             tk_xxxy_xxyzzz, \
                             tk_xxxy_xxzzzz, \
                             tk_xxxy_xyyyyy, \
                             tk_xxxy_xyyyyz, \
                             tk_xxxy_xyyyzz, \
                             tk_xxxy_xyyzzz, \
                             tk_xxxy_xyzzzz, \
                             tk_xxxy_xzzzzz, \
                             tk_xxxy_yyyyyy, \
                             tk_xxxy_yyyyyz, \
                             tk_xxxy_yyyyzz, \
                             tk_xxxy_yyyzzz, \
                             tk_xxxy_yyzzzz, \
                             tk_xxxy_yzzzzz, \
                             tk_xxxy_zzzzzz, \
                             ts_xxxy_xxxxxx, \
                             ts_xxxy_xxxxxy, \
                             ts_xxxy_xxxxxz, \
                             ts_xxxy_xxxxyy, \
                             ts_xxxy_xxxxyz, \
                             ts_xxxy_xxxxzz, \
                             ts_xxxy_xxxyyy, \
                             ts_xxxy_xxxyyz, \
                             ts_xxxy_xxxyzz, \
                             ts_xxxy_xxxzzz, \
                             ts_xxxy_xxyyyy, \
                             ts_xxxy_xxyyyz, \
                             ts_xxxy_xxyyzz, \
                             ts_xxxy_xxyzzz, \
                             ts_xxxy_xxzzzz, \
                             ts_xxxy_xyyyyy, \
                             ts_xxxy_xyyyyz, \
                             ts_xxxy_xyyyzz, \
                             ts_xxxy_xyyzzz, \
                             ts_xxxy_xyzzzz, \
                             ts_xxxy_xzzzzz, \
                             ts_xxxy_yyyyyy, \
                             ts_xxxy_yyyyyz, \
                             ts_xxxy_yyyyzz, \
                             ts_xxxy_yyyzzz, \
                             ts_xxxy_yyzzzz, \
                             ts_xxxy_yzzzzz, \
                             ts_xxxy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxy_xxxxxx[i] = tk_xxx_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxxy_xxxxxx[i] * fz_0;

        tk_xxxy_xxxxxy[i] = tk_xxx_xxxxx[i] * fe_0 + tk_xxx_xxxxxy[i] * pa_y[i] + 2.0 * ts_xxxy_xxxxxy[i] * fz_0;

        tk_xxxy_xxxxxz[i] = tk_xxx_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxxy_xxxxxz[i] * fz_0;

        tk_xxxy_xxxxyy[i] = 2.0 * tk_xxx_xxxxy[i] * fe_0 + tk_xxx_xxxxyy[i] * pa_y[i] + 2.0 * ts_xxxy_xxxxyy[i] * fz_0;

        tk_xxxy_xxxxyz[i] = tk_xxx_xxxxz[i] * fe_0 + tk_xxx_xxxxyz[i] * pa_y[i] + 2.0 * ts_xxxy_xxxxyz[i] * fz_0;

        tk_xxxy_xxxxzz[i] = tk_xxx_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxxy_xxxxzz[i] * fz_0;

        tk_xxxy_xxxyyy[i] = 3.0 * tk_xxx_xxxyy[i] * fe_0 + tk_xxx_xxxyyy[i] * pa_y[i] + 2.0 * ts_xxxy_xxxyyy[i] * fz_0;

        tk_xxxy_xxxyyz[i] = 2.0 * tk_xxx_xxxyz[i] * fe_0 + tk_xxx_xxxyyz[i] * pa_y[i] + 2.0 * ts_xxxy_xxxyyz[i] * fz_0;

        tk_xxxy_xxxyzz[i] = tk_xxx_xxxzz[i] * fe_0 + tk_xxx_xxxyzz[i] * pa_y[i] + 2.0 * ts_xxxy_xxxyzz[i] * fz_0;

        tk_xxxy_xxxzzz[i] = tk_xxx_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxxy_xxxzzz[i] * fz_0;

        tk_xxxy_xxyyyy[i] = 4.0 * tk_xxx_xxyyy[i] * fe_0 + tk_xxx_xxyyyy[i] * pa_y[i] + 2.0 * ts_xxxy_xxyyyy[i] * fz_0;

        tk_xxxy_xxyyyz[i] = 3.0 * tk_xxx_xxyyz[i] * fe_0 + tk_xxx_xxyyyz[i] * pa_y[i] + 2.0 * ts_xxxy_xxyyyz[i] * fz_0;

        tk_xxxy_xxyyzz[i] = 2.0 * tk_xxx_xxyzz[i] * fe_0 + tk_xxx_xxyyzz[i] * pa_y[i] + 2.0 * ts_xxxy_xxyyzz[i] * fz_0;

        tk_xxxy_xxyzzz[i] = tk_xxx_xxzzz[i] * fe_0 + tk_xxx_xxyzzz[i] * pa_y[i] + 2.0 * ts_xxxy_xxyzzz[i] * fz_0;

        tk_xxxy_xxzzzz[i] = tk_xxx_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxxy_xxzzzz[i] * fz_0;

        tk_xxxy_xyyyyy[i] = 5.0 * tk_xxx_xyyyy[i] * fe_0 + tk_xxx_xyyyyy[i] * pa_y[i] + 2.0 * ts_xxxy_xyyyyy[i] * fz_0;

        tk_xxxy_xyyyyz[i] = 4.0 * tk_xxx_xyyyz[i] * fe_0 + tk_xxx_xyyyyz[i] * pa_y[i] + 2.0 * ts_xxxy_xyyyyz[i] * fz_0;

        tk_xxxy_xyyyzz[i] = 3.0 * tk_xxx_xyyzz[i] * fe_0 + tk_xxx_xyyyzz[i] * pa_y[i] + 2.0 * ts_xxxy_xyyyzz[i] * fz_0;

        tk_xxxy_xyyzzz[i] = 2.0 * tk_xxx_xyzzz[i] * fe_0 + tk_xxx_xyyzzz[i] * pa_y[i] + 2.0 * ts_xxxy_xyyzzz[i] * fz_0;

        tk_xxxy_xyzzzz[i] = tk_xxx_xzzzz[i] * fe_0 + tk_xxx_xyzzzz[i] * pa_y[i] + 2.0 * ts_xxxy_xyzzzz[i] * fz_0;

        tk_xxxy_xzzzzz[i] = tk_xxx_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxxy_xzzzzz[i] * fz_0;

        tk_xxxy_yyyyyy[i] = 6.0 * tk_xxx_yyyyy[i] * fe_0 + tk_xxx_yyyyyy[i] * pa_y[i] + 2.0 * ts_xxxy_yyyyyy[i] * fz_0;

        tk_xxxy_yyyyyz[i] = 5.0 * tk_xxx_yyyyz[i] * fe_0 + tk_xxx_yyyyyz[i] * pa_y[i] + 2.0 * ts_xxxy_yyyyyz[i] * fz_0;

        tk_xxxy_yyyyzz[i] = 4.0 * tk_xxx_yyyzz[i] * fe_0 + tk_xxx_yyyyzz[i] * pa_y[i] + 2.0 * ts_xxxy_yyyyzz[i] * fz_0;

        tk_xxxy_yyyzzz[i] = 3.0 * tk_xxx_yyzzz[i] * fe_0 + tk_xxx_yyyzzz[i] * pa_y[i] + 2.0 * ts_xxxy_yyyzzz[i] * fz_0;

        tk_xxxy_yyzzzz[i] = 2.0 * tk_xxx_yzzzz[i] * fe_0 + tk_xxx_yyzzzz[i] * pa_y[i] + 2.0 * ts_xxxy_yyzzzz[i] * fz_0;

        tk_xxxy_yzzzzz[i] = tk_xxx_zzzzz[i] * fe_0 + tk_xxx_yzzzzz[i] * pa_y[i] + 2.0 * ts_xxxy_yzzzzz[i] * fz_0;

        tk_xxxy_zzzzzz[i] = tk_xxx_zzzzzz[i] * pa_y[i] + 2.0 * ts_xxxy_zzzzzz[i] * fz_0;
    }

    // Set up 56-84 components of targeted buffer : GI

    auto tk_xxxz_xxxxxx = pbuffer.data(idx_kin_gi + 56);

    auto tk_xxxz_xxxxxy = pbuffer.data(idx_kin_gi + 57);

    auto tk_xxxz_xxxxxz = pbuffer.data(idx_kin_gi + 58);

    auto tk_xxxz_xxxxyy = pbuffer.data(idx_kin_gi + 59);

    auto tk_xxxz_xxxxyz = pbuffer.data(idx_kin_gi + 60);

    auto tk_xxxz_xxxxzz = pbuffer.data(idx_kin_gi + 61);

    auto tk_xxxz_xxxyyy = pbuffer.data(idx_kin_gi + 62);

    auto tk_xxxz_xxxyyz = pbuffer.data(idx_kin_gi + 63);

    auto tk_xxxz_xxxyzz = pbuffer.data(idx_kin_gi + 64);

    auto tk_xxxz_xxxzzz = pbuffer.data(idx_kin_gi + 65);

    auto tk_xxxz_xxyyyy = pbuffer.data(idx_kin_gi + 66);

    auto tk_xxxz_xxyyyz = pbuffer.data(idx_kin_gi + 67);

    auto tk_xxxz_xxyyzz = pbuffer.data(idx_kin_gi + 68);

    auto tk_xxxz_xxyzzz = pbuffer.data(idx_kin_gi + 69);

    auto tk_xxxz_xxzzzz = pbuffer.data(idx_kin_gi + 70);

    auto tk_xxxz_xyyyyy = pbuffer.data(idx_kin_gi + 71);

    auto tk_xxxz_xyyyyz = pbuffer.data(idx_kin_gi + 72);

    auto tk_xxxz_xyyyzz = pbuffer.data(idx_kin_gi + 73);

    auto tk_xxxz_xyyzzz = pbuffer.data(idx_kin_gi + 74);

    auto tk_xxxz_xyzzzz = pbuffer.data(idx_kin_gi + 75);

    auto tk_xxxz_xzzzzz = pbuffer.data(idx_kin_gi + 76);

    auto tk_xxxz_yyyyyy = pbuffer.data(idx_kin_gi + 77);

    auto tk_xxxz_yyyyyz = pbuffer.data(idx_kin_gi + 78);

    auto tk_xxxz_yyyyzz = pbuffer.data(idx_kin_gi + 79);

    auto tk_xxxz_yyyzzz = pbuffer.data(idx_kin_gi + 80);

    auto tk_xxxz_yyzzzz = pbuffer.data(idx_kin_gi + 81);

    auto tk_xxxz_yzzzzz = pbuffer.data(idx_kin_gi + 82);

    auto tk_xxxz_zzzzzz = pbuffer.data(idx_kin_gi + 83);

#pragma omp simd aligned(pa_z,               \
                             tk_xxx_xxxxx,   \
                             tk_xxx_xxxxxx,  \
                             tk_xxx_xxxxxy,  \
                             tk_xxx_xxxxxz,  \
                             tk_xxx_xxxxy,   \
                             tk_xxx_xxxxyy,  \
                             tk_xxx_xxxxyz,  \
                             tk_xxx_xxxxz,   \
                             tk_xxx_xxxxzz,  \
                             tk_xxx_xxxyy,   \
                             tk_xxx_xxxyyy,  \
                             tk_xxx_xxxyyz,  \
                             tk_xxx_xxxyz,   \
                             tk_xxx_xxxyzz,  \
                             tk_xxx_xxxzz,   \
                             tk_xxx_xxxzzz,  \
                             tk_xxx_xxyyy,   \
                             tk_xxx_xxyyyy,  \
                             tk_xxx_xxyyyz,  \
                             tk_xxx_xxyyz,   \
                             tk_xxx_xxyyzz,  \
                             tk_xxx_xxyzz,   \
                             tk_xxx_xxyzzz,  \
                             tk_xxx_xxzzz,   \
                             tk_xxx_xxzzzz,  \
                             tk_xxx_xyyyy,   \
                             tk_xxx_xyyyyy,  \
                             tk_xxx_xyyyyz,  \
                             tk_xxx_xyyyz,   \
                             tk_xxx_xyyyzz,  \
                             tk_xxx_xyyzz,   \
                             tk_xxx_xyyzzz,  \
                             tk_xxx_xyzzz,   \
                             tk_xxx_xyzzzz,  \
                             tk_xxx_xzzzz,   \
                             tk_xxx_xzzzzz,  \
                             tk_xxx_yyyyy,   \
                             tk_xxx_yyyyyy,  \
                             tk_xxx_yyyyyz,  \
                             tk_xxx_yyyyz,   \
                             tk_xxx_yyyyzz,  \
                             tk_xxx_yyyzz,   \
                             tk_xxx_yyyzzz,  \
                             tk_xxx_yyzzz,   \
                             tk_xxx_yyzzzz,  \
                             tk_xxx_yzzzz,   \
                             tk_xxx_yzzzzz,  \
                             tk_xxx_zzzzz,   \
                             tk_xxx_zzzzzz,  \
                             tk_xxxz_xxxxxx, \
                             tk_xxxz_xxxxxy, \
                             tk_xxxz_xxxxxz, \
                             tk_xxxz_xxxxyy, \
                             tk_xxxz_xxxxyz, \
                             tk_xxxz_xxxxzz, \
                             tk_xxxz_xxxyyy, \
                             tk_xxxz_xxxyyz, \
                             tk_xxxz_xxxyzz, \
                             tk_xxxz_xxxzzz, \
                             tk_xxxz_xxyyyy, \
                             tk_xxxz_xxyyyz, \
                             tk_xxxz_xxyyzz, \
                             tk_xxxz_xxyzzz, \
                             tk_xxxz_xxzzzz, \
                             tk_xxxz_xyyyyy, \
                             tk_xxxz_xyyyyz, \
                             tk_xxxz_xyyyzz, \
                             tk_xxxz_xyyzzz, \
                             tk_xxxz_xyzzzz, \
                             tk_xxxz_xzzzzz, \
                             tk_xxxz_yyyyyy, \
                             tk_xxxz_yyyyyz, \
                             tk_xxxz_yyyyzz, \
                             tk_xxxz_yyyzzz, \
                             tk_xxxz_yyzzzz, \
                             tk_xxxz_yzzzzz, \
                             tk_xxxz_zzzzzz, \
                             ts_xxxz_xxxxxx, \
                             ts_xxxz_xxxxxy, \
                             ts_xxxz_xxxxxz, \
                             ts_xxxz_xxxxyy, \
                             ts_xxxz_xxxxyz, \
                             ts_xxxz_xxxxzz, \
                             ts_xxxz_xxxyyy, \
                             ts_xxxz_xxxyyz, \
                             ts_xxxz_xxxyzz, \
                             ts_xxxz_xxxzzz, \
                             ts_xxxz_xxyyyy, \
                             ts_xxxz_xxyyyz, \
                             ts_xxxz_xxyyzz, \
                             ts_xxxz_xxyzzz, \
                             ts_xxxz_xxzzzz, \
                             ts_xxxz_xyyyyy, \
                             ts_xxxz_xyyyyz, \
                             ts_xxxz_xyyyzz, \
                             ts_xxxz_xyyzzz, \
                             ts_xxxz_xyzzzz, \
                             ts_xxxz_xzzzzz, \
                             ts_xxxz_yyyyyy, \
                             ts_xxxz_yyyyyz, \
                             ts_xxxz_yyyyzz, \
                             ts_xxxz_yyyzzz, \
                             ts_xxxz_yyzzzz, \
                             ts_xxxz_yzzzzz, \
                             ts_xxxz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxz_xxxxxx[i] = tk_xxx_xxxxxx[i] * pa_z[i] + 2.0 * ts_xxxz_xxxxxx[i] * fz_0;

        tk_xxxz_xxxxxy[i] = tk_xxx_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxxz_xxxxxy[i] * fz_0;

        tk_xxxz_xxxxxz[i] = tk_xxx_xxxxx[i] * fe_0 + tk_xxx_xxxxxz[i] * pa_z[i] + 2.0 * ts_xxxz_xxxxxz[i] * fz_0;

        tk_xxxz_xxxxyy[i] = tk_xxx_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxxz_xxxxyy[i] * fz_0;

        tk_xxxz_xxxxyz[i] = tk_xxx_xxxxy[i] * fe_0 + tk_xxx_xxxxyz[i] * pa_z[i] + 2.0 * ts_xxxz_xxxxyz[i] * fz_0;

        tk_xxxz_xxxxzz[i] = 2.0 * tk_xxx_xxxxz[i] * fe_0 + tk_xxx_xxxxzz[i] * pa_z[i] + 2.0 * ts_xxxz_xxxxzz[i] * fz_0;

        tk_xxxz_xxxyyy[i] = tk_xxx_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxxz_xxxyyy[i] * fz_0;

        tk_xxxz_xxxyyz[i] = tk_xxx_xxxyy[i] * fe_0 + tk_xxx_xxxyyz[i] * pa_z[i] + 2.0 * ts_xxxz_xxxyyz[i] * fz_0;

        tk_xxxz_xxxyzz[i] = 2.0 * tk_xxx_xxxyz[i] * fe_0 + tk_xxx_xxxyzz[i] * pa_z[i] + 2.0 * ts_xxxz_xxxyzz[i] * fz_0;

        tk_xxxz_xxxzzz[i] = 3.0 * tk_xxx_xxxzz[i] * fe_0 + tk_xxx_xxxzzz[i] * pa_z[i] + 2.0 * ts_xxxz_xxxzzz[i] * fz_0;

        tk_xxxz_xxyyyy[i] = tk_xxx_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxxz_xxyyyy[i] * fz_0;

        tk_xxxz_xxyyyz[i] = tk_xxx_xxyyy[i] * fe_0 + tk_xxx_xxyyyz[i] * pa_z[i] + 2.0 * ts_xxxz_xxyyyz[i] * fz_0;

        tk_xxxz_xxyyzz[i] = 2.0 * tk_xxx_xxyyz[i] * fe_0 + tk_xxx_xxyyzz[i] * pa_z[i] + 2.0 * ts_xxxz_xxyyzz[i] * fz_0;

        tk_xxxz_xxyzzz[i] = 3.0 * tk_xxx_xxyzz[i] * fe_0 + tk_xxx_xxyzzz[i] * pa_z[i] + 2.0 * ts_xxxz_xxyzzz[i] * fz_0;

        tk_xxxz_xxzzzz[i] = 4.0 * tk_xxx_xxzzz[i] * fe_0 + tk_xxx_xxzzzz[i] * pa_z[i] + 2.0 * ts_xxxz_xxzzzz[i] * fz_0;

        tk_xxxz_xyyyyy[i] = tk_xxx_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxxz_xyyyyy[i] * fz_0;

        tk_xxxz_xyyyyz[i] = tk_xxx_xyyyy[i] * fe_0 + tk_xxx_xyyyyz[i] * pa_z[i] + 2.0 * ts_xxxz_xyyyyz[i] * fz_0;

        tk_xxxz_xyyyzz[i] = 2.0 * tk_xxx_xyyyz[i] * fe_0 + tk_xxx_xyyyzz[i] * pa_z[i] + 2.0 * ts_xxxz_xyyyzz[i] * fz_0;

        tk_xxxz_xyyzzz[i] = 3.0 * tk_xxx_xyyzz[i] * fe_0 + tk_xxx_xyyzzz[i] * pa_z[i] + 2.0 * ts_xxxz_xyyzzz[i] * fz_0;

        tk_xxxz_xyzzzz[i] = 4.0 * tk_xxx_xyzzz[i] * fe_0 + tk_xxx_xyzzzz[i] * pa_z[i] + 2.0 * ts_xxxz_xyzzzz[i] * fz_0;

        tk_xxxz_xzzzzz[i] = 5.0 * tk_xxx_xzzzz[i] * fe_0 + tk_xxx_xzzzzz[i] * pa_z[i] + 2.0 * ts_xxxz_xzzzzz[i] * fz_0;

        tk_xxxz_yyyyyy[i] = tk_xxx_yyyyyy[i] * pa_z[i] + 2.0 * ts_xxxz_yyyyyy[i] * fz_0;

        tk_xxxz_yyyyyz[i] = tk_xxx_yyyyy[i] * fe_0 + tk_xxx_yyyyyz[i] * pa_z[i] + 2.0 * ts_xxxz_yyyyyz[i] * fz_0;

        tk_xxxz_yyyyzz[i] = 2.0 * tk_xxx_yyyyz[i] * fe_0 + tk_xxx_yyyyzz[i] * pa_z[i] + 2.0 * ts_xxxz_yyyyzz[i] * fz_0;

        tk_xxxz_yyyzzz[i] = 3.0 * tk_xxx_yyyzz[i] * fe_0 + tk_xxx_yyyzzz[i] * pa_z[i] + 2.0 * ts_xxxz_yyyzzz[i] * fz_0;

        tk_xxxz_yyzzzz[i] = 4.0 * tk_xxx_yyzzz[i] * fe_0 + tk_xxx_yyzzzz[i] * pa_z[i] + 2.0 * ts_xxxz_yyzzzz[i] * fz_0;

        tk_xxxz_yzzzzz[i] = 5.0 * tk_xxx_yzzzz[i] * fe_0 + tk_xxx_yzzzzz[i] * pa_z[i] + 2.0 * ts_xxxz_yzzzzz[i] * fz_0;

        tk_xxxz_zzzzzz[i] = 6.0 * tk_xxx_zzzzz[i] * fe_0 + tk_xxx_zzzzzz[i] * pa_z[i] + 2.0 * ts_xxxz_zzzzzz[i] * fz_0;
    }

    // Set up 84-112 components of targeted buffer : GI

    auto tk_xxyy_xxxxxx = pbuffer.data(idx_kin_gi + 84);

    auto tk_xxyy_xxxxxy = pbuffer.data(idx_kin_gi + 85);

    auto tk_xxyy_xxxxxz = pbuffer.data(idx_kin_gi + 86);

    auto tk_xxyy_xxxxyy = pbuffer.data(idx_kin_gi + 87);

    auto tk_xxyy_xxxxyz = pbuffer.data(idx_kin_gi + 88);

    auto tk_xxyy_xxxxzz = pbuffer.data(idx_kin_gi + 89);

    auto tk_xxyy_xxxyyy = pbuffer.data(idx_kin_gi + 90);

    auto tk_xxyy_xxxyyz = pbuffer.data(idx_kin_gi + 91);

    auto tk_xxyy_xxxyzz = pbuffer.data(idx_kin_gi + 92);

    auto tk_xxyy_xxxzzz = pbuffer.data(idx_kin_gi + 93);

    auto tk_xxyy_xxyyyy = pbuffer.data(idx_kin_gi + 94);

    auto tk_xxyy_xxyyyz = pbuffer.data(idx_kin_gi + 95);

    auto tk_xxyy_xxyyzz = pbuffer.data(idx_kin_gi + 96);

    auto tk_xxyy_xxyzzz = pbuffer.data(idx_kin_gi + 97);

    auto tk_xxyy_xxzzzz = pbuffer.data(idx_kin_gi + 98);

    auto tk_xxyy_xyyyyy = pbuffer.data(idx_kin_gi + 99);

    auto tk_xxyy_xyyyyz = pbuffer.data(idx_kin_gi + 100);

    auto tk_xxyy_xyyyzz = pbuffer.data(idx_kin_gi + 101);

    auto tk_xxyy_xyyzzz = pbuffer.data(idx_kin_gi + 102);

    auto tk_xxyy_xyzzzz = pbuffer.data(idx_kin_gi + 103);

    auto tk_xxyy_xzzzzz = pbuffer.data(idx_kin_gi + 104);

    auto tk_xxyy_yyyyyy = pbuffer.data(idx_kin_gi + 105);

    auto tk_xxyy_yyyyyz = pbuffer.data(idx_kin_gi + 106);

    auto tk_xxyy_yyyyzz = pbuffer.data(idx_kin_gi + 107);

    auto tk_xxyy_yyyzzz = pbuffer.data(idx_kin_gi + 108);

    auto tk_xxyy_yyzzzz = pbuffer.data(idx_kin_gi + 109);

    auto tk_xxyy_yzzzzz = pbuffer.data(idx_kin_gi + 110);

    auto tk_xxyy_zzzzzz = pbuffer.data(idx_kin_gi + 111);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tk_xx_xxxxxx,   \
                             tk_xx_xxxxxz,   \
                             tk_xx_xxxxzz,   \
                             tk_xx_xxxzzz,   \
                             tk_xx_xxzzzz,   \
                             tk_xx_xzzzzz,   \
                             tk_xxy_xxxxxx,  \
                             tk_xxy_xxxxxz,  \
                             tk_xxy_xxxxzz,  \
                             tk_xxy_xxxzzz,  \
                             tk_xxy_xxzzzz,  \
                             tk_xxy_xzzzzz,  \
                             tk_xxyy_xxxxxx, \
                             tk_xxyy_xxxxxy, \
                             tk_xxyy_xxxxxz, \
                             tk_xxyy_xxxxyy, \
                             tk_xxyy_xxxxyz, \
                             tk_xxyy_xxxxzz, \
                             tk_xxyy_xxxyyy, \
                             tk_xxyy_xxxyyz, \
                             tk_xxyy_xxxyzz, \
                             tk_xxyy_xxxzzz, \
                             tk_xxyy_xxyyyy, \
                             tk_xxyy_xxyyyz, \
                             tk_xxyy_xxyyzz, \
                             tk_xxyy_xxyzzz, \
                             tk_xxyy_xxzzzz, \
                             tk_xxyy_xyyyyy, \
                             tk_xxyy_xyyyyz, \
                             tk_xxyy_xyyyzz, \
                             tk_xxyy_xyyzzz, \
                             tk_xxyy_xyzzzz, \
                             tk_xxyy_xzzzzz, \
                             tk_xxyy_yyyyyy, \
                             tk_xxyy_yyyyyz, \
                             tk_xxyy_yyyyzz, \
                             tk_xxyy_yyyzzz, \
                             tk_xxyy_yyzzzz, \
                             tk_xxyy_yzzzzz, \
                             tk_xxyy_zzzzzz, \
                             tk_xyy_xxxxxy,  \
                             tk_xyy_xxxxy,   \
                             tk_xyy_xxxxyy,  \
                             tk_xyy_xxxxyz,  \
                             tk_xyy_xxxyy,   \
                             tk_xyy_xxxyyy,  \
                             tk_xyy_xxxyyz,  \
                             tk_xyy_xxxyz,   \
                             tk_xyy_xxxyzz,  \
                             tk_xyy_xxyyy,   \
                             tk_xyy_xxyyyy,  \
                             tk_xyy_xxyyyz,  \
                             tk_xyy_xxyyz,   \
                             tk_xyy_xxyyzz,  \
                             tk_xyy_xxyzz,   \
                             tk_xyy_xxyzzz,  \
                             tk_xyy_xyyyy,   \
                             tk_xyy_xyyyyy,  \
                             tk_xyy_xyyyyz,  \
                             tk_xyy_xyyyz,   \
                             tk_xyy_xyyyzz,  \
                             tk_xyy_xyyzz,   \
                             tk_xyy_xyyzzz,  \
                             tk_xyy_xyzzz,   \
                             tk_xyy_xyzzzz,  \
                             tk_xyy_yyyyy,   \
                             tk_xyy_yyyyyy,  \
                             tk_xyy_yyyyyz,  \
                             tk_xyy_yyyyz,   \
                             tk_xyy_yyyyzz,  \
                             tk_xyy_yyyzz,   \
                             tk_xyy_yyyzzz,  \
                             tk_xyy_yyzzz,   \
                             tk_xyy_yyzzzz,  \
                             tk_xyy_yzzzz,   \
                             tk_xyy_yzzzzz,  \
                             tk_xyy_zzzzzz,  \
                             tk_yy_xxxxxy,   \
                             tk_yy_xxxxyy,   \
                             tk_yy_xxxxyz,   \
                             tk_yy_xxxyyy,   \
                             tk_yy_xxxyyz,   \
                             tk_yy_xxxyzz,   \
                             tk_yy_xxyyyy,   \
                             tk_yy_xxyyyz,   \
                             tk_yy_xxyyzz,   \
                             tk_yy_xxyzzz,   \
                             tk_yy_xyyyyy,   \
                             tk_yy_xyyyyz,   \
                             tk_yy_xyyyzz,   \
                             tk_yy_xyyzzz,   \
                             tk_yy_xyzzzz,   \
                             tk_yy_yyyyyy,   \
                             tk_yy_yyyyyz,   \
                             tk_yy_yyyyzz,   \
                             tk_yy_yyyzzz,   \
                             tk_yy_yyzzzz,   \
                             tk_yy_yzzzzz,   \
                             tk_yy_zzzzzz,   \
                             ts_xx_xxxxxx,   \
                             ts_xx_xxxxxz,   \
                             ts_xx_xxxxzz,   \
                             ts_xx_xxxzzz,   \
                             ts_xx_xxzzzz,   \
                             ts_xx_xzzzzz,   \
                             ts_xxyy_xxxxxx, \
                             ts_xxyy_xxxxxy, \
                             ts_xxyy_xxxxxz, \
                             ts_xxyy_xxxxyy, \
                             ts_xxyy_xxxxyz, \
                             ts_xxyy_xxxxzz, \
                             ts_xxyy_xxxyyy, \
                             ts_xxyy_xxxyyz, \
                             ts_xxyy_xxxyzz, \
                             ts_xxyy_xxxzzz, \
                             ts_xxyy_xxyyyy, \
                             ts_xxyy_xxyyyz, \
                             ts_xxyy_xxyyzz, \
                             ts_xxyy_xxyzzz, \
                             ts_xxyy_xxzzzz, \
                             ts_xxyy_xyyyyy, \
                             ts_xxyy_xyyyyz, \
                             ts_xxyy_xyyyzz, \
                             ts_xxyy_xyyzzz, \
                             ts_xxyy_xyzzzz, \
                             ts_xxyy_xzzzzz, \
                             ts_xxyy_yyyyyy, \
                             ts_xxyy_yyyyyz, \
                             ts_xxyy_yyyyzz, \
                             ts_xxyy_yyyzzz, \
                             ts_xxyy_yyzzzz, \
                             ts_xxyy_yzzzzz, \
                             ts_xxyy_zzzzzz, \
                             ts_yy_xxxxxy,   \
                             ts_yy_xxxxyy,   \
                             ts_yy_xxxxyz,   \
                             ts_yy_xxxyyy,   \
                             ts_yy_xxxyyz,   \
                             ts_yy_xxxyzz,   \
                             ts_yy_xxyyyy,   \
                             ts_yy_xxyyyz,   \
                             ts_yy_xxyyzz,   \
                             ts_yy_xxyzzz,   \
                             ts_yy_xyyyyy,   \
                             ts_yy_xyyyyz,   \
                             ts_yy_xyyyzz,   \
                             ts_yy_xyyzzz,   \
                             ts_yy_xyzzzz,   \
                             ts_yy_yyyyyy,   \
                             ts_yy_yyyyyz,   \
                             ts_yy_yyyyzz,   \
                             ts_yy_yyyzzz,   \
                             ts_yy_yyzzzz,   \
                             ts_yy_yzzzzz,   \
                             ts_yy_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyy_xxxxxx[i] =
            -2.0 * ts_xx_xxxxxx[i] * fbe_0 * fz_0 + tk_xx_xxxxxx[i] * fe_0 + tk_xxy_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxyy_xxxxxx[i] * fz_0;

        tk_xxyy_xxxxxy[i] = -2.0 * ts_yy_xxxxxy[i] * fbe_0 * fz_0 + tk_yy_xxxxxy[i] * fe_0 + 5.0 * tk_xyy_xxxxy[i] * fe_0 +
                            tk_xyy_xxxxxy[i] * pa_x[i] + 2.0 * ts_xxyy_xxxxxy[i] * fz_0;

        tk_xxyy_xxxxxz[i] =
            -2.0 * ts_xx_xxxxxz[i] * fbe_0 * fz_0 + tk_xx_xxxxxz[i] * fe_0 + tk_xxy_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxyy_xxxxxz[i] * fz_0;

        tk_xxyy_xxxxyy[i] = -2.0 * ts_yy_xxxxyy[i] * fbe_0 * fz_0 + tk_yy_xxxxyy[i] * fe_0 + 4.0 * tk_xyy_xxxyy[i] * fe_0 +
                            tk_xyy_xxxxyy[i] * pa_x[i] + 2.0 * ts_xxyy_xxxxyy[i] * fz_0;

        tk_xxyy_xxxxyz[i] = -2.0 * ts_yy_xxxxyz[i] * fbe_0 * fz_0 + tk_yy_xxxxyz[i] * fe_0 + 4.0 * tk_xyy_xxxyz[i] * fe_0 +
                            tk_xyy_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxyy_xxxxyz[i] * fz_0;

        tk_xxyy_xxxxzz[i] =
            -2.0 * ts_xx_xxxxzz[i] * fbe_0 * fz_0 + tk_xx_xxxxzz[i] * fe_0 + tk_xxy_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxyy_xxxxzz[i] * fz_0;

        tk_xxyy_xxxyyy[i] = -2.0 * ts_yy_xxxyyy[i] * fbe_0 * fz_0 + tk_yy_xxxyyy[i] * fe_0 + 3.0 * tk_xyy_xxyyy[i] * fe_0 +
                            tk_xyy_xxxyyy[i] * pa_x[i] + 2.0 * ts_xxyy_xxxyyy[i] * fz_0;

        tk_xxyy_xxxyyz[i] = -2.0 * ts_yy_xxxyyz[i] * fbe_0 * fz_0 + tk_yy_xxxyyz[i] * fe_0 + 3.0 * tk_xyy_xxyyz[i] * fe_0 +
                            tk_xyy_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxyy_xxxyyz[i] * fz_0;

        tk_xxyy_xxxyzz[i] = -2.0 * ts_yy_xxxyzz[i] * fbe_0 * fz_0 + tk_yy_xxxyzz[i] * fe_0 + 3.0 * tk_xyy_xxyzz[i] * fe_0 +
                            tk_xyy_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxyy_xxxyzz[i] * fz_0;

        tk_xxyy_xxxzzz[i] =
            -2.0 * ts_xx_xxxzzz[i] * fbe_0 * fz_0 + tk_xx_xxxzzz[i] * fe_0 + tk_xxy_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxyy_xxxzzz[i] * fz_0;

        tk_xxyy_xxyyyy[i] = -2.0 * ts_yy_xxyyyy[i] * fbe_0 * fz_0 + tk_yy_xxyyyy[i] * fe_0 + 2.0 * tk_xyy_xyyyy[i] * fe_0 +
                            tk_xyy_xxyyyy[i] * pa_x[i] + 2.0 * ts_xxyy_xxyyyy[i] * fz_0;

        tk_xxyy_xxyyyz[i] = -2.0 * ts_yy_xxyyyz[i] * fbe_0 * fz_0 + tk_yy_xxyyyz[i] * fe_0 + 2.0 * tk_xyy_xyyyz[i] * fe_0 +
                            tk_xyy_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxyy_xxyyyz[i] * fz_0;

        tk_xxyy_xxyyzz[i] = -2.0 * ts_yy_xxyyzz[i] * fbe_0 * fz_0 + tk_yy_xxyyzz[i] * fe_0 + 2.0 * tk_xyy_xyyzz[i] * fe_0 +
                            tk_xyy_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxyy_xxyyzz[i] * fz_0;

        tk_xxyy_xxyzzz[i] = -2.0 * ts_yy_xxyzzz[i] * fbe_0 * fz_0 + tk_yy_xxyzzz[i] * fe_0 + 2.0 * tk_xyy_xyzzz[i] * fe_0 +
                            tk_xyy_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxyy_xxyzzz[i] * fz_0;

        tk_xxyy_xxzzzz[i] =
            -2.0 * ts_xx_xxzzzz[i] * fbe_0 * fz_0 + tk_xx_xxzzzz[i] * fe_0 + tk_xxy_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxyy_xxzzzz[i] * fz_0;

        tk_xxyy_xyyyyy[i] = -2.0 * ts_yy_xyyyyy[i] * fbe_0 * fz_0 + tk_yy_xyyyyy[i] * fe_0 + tk_xyy_yyyyy[i] * fe_0 + tk_xyy_xyyyyy[i] * pa_x[i] +
                            2.0 * ts_xxyy_xyyyyy[i] * fz_0;

        tk_xxyy_xyyyyz[i] = -2.0 * ts_yy_xyyyyz[i] * fbe_0 * fz_0 + tk_yy_xyyyyz[i] * fe_0 + tk_xyy_yyyyz[i] * fe_0 + tk_xyy_xyyyyz[i] * pa_x[i] +
                            2.0 * ts_xxyy_xyyyyz[i] * fz_0;

        tk_xxyy_xyyyzz[i] = -2.0 * ts_yy_xyyyzz[i] * fbe_0 * fz_0 + tk_yy_xyyyzz[i] * fe_0 + tk_xyy_yyyzz[i] * fe_0 + tk_xyy_xyyyzz[i] * pa_x[i] +
                            2.0 * ts_xxyy_xyyyzz[i] * fz_0;

        tk_xxyy_xyyzzz[i] = -2.0 * ts_yy_xyyzzz[i] * fbe_0 * fz_0 + tk_yy_xyyzzz[i] * fe_0 + tk_xyy_yyzzz[i] * fe_0 + tk_xyy_xyyzzz[i] * pa_x[i] +
                            2.0 * ts_xxyy_xyyzzz[i] * fz_0;

        tk_xxyy_xyzzzz[i] = -2.0 * ts_yy_xyzzzz[i] * fbe_0 * fz_0 + tk_yy_xyzzzz[i] * fe_0 + tk_xyy_yzzzz[i] * fe_0 + tk_xyy_xyzzzz[i] * pa_x[i] +
                            2.0 * ts_xxyy_xyzzzz[i] * fz_0;

        tk_xxyy_xzzzzz[i] =
            -2.0 * ts_xx_xzzzzz[i] * fbe_0 * fz_0 + tk_xx_xzzzzz[i] * fe_0 + tk_xxy_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxyy_xzzzzz[i] * fz_0;

        tk_xxyy_yyyyyy[i] =
            -2.0 * ts_yy_yyyyyy[i] * fbe_0 * fz_0 + tk_yy_yyyyyy[i] * fe_0 + tk_xyy_yyyyyy[i] * pa_x[i] + 2.0 * ts_xxyy_yyyyyy[i] * fz_0;

        tk_xxyy_yyyyyz[i] =
            -2.0 * ts_yy_yyyyyz[i] * fbe_0 * fz_0 + tk_yy_yyyyyz[i] * fe_0 + tk_xyy_yyyyyz[i] * pa_x[i] + 2.0 * ts_xxyy_yyyyyz[i] * fz_0;

        tk_xxyy_yyyyzz[i] =
            -2.0 * ts_yy_yyyyzz[i] * fbe_0 * fz_0 + tk_yy_yyyyzz[i] * fe_0 + tk_xyy_yyyyzz[i] * pa_x[i] + 2.0 * ts_xxyy_yyyyzz[i] * fz_0;

        tk_xxyy_yyyzzz[i] =
            -2.0 * ts_yy_yyyzzz[i] * fbe_0 * fz_0 + tk_yy_yyyzzz[i] * fe_0 + tk_xyy_yyyzzz[i] * pa_x[i] + 2.0 * ts_xxyy_yyyzzz[i] * fz_0;

        tk_xxyy_yyzzzz[i] =
            -2.0 * ts_yy_yyzzzz[i] * fbe_0 * fz_0 + tk_yy_yyzzzz[i] * fe_0 + tk_xyy_yyzzzz[i] * pa_x[i] + 2.0 * ts_xxyy_yyzzzz[i] * fz_0;

        tk_xxyy_yzzzzz[i] =
            -2.0 * ts_yy_yzzzzz[i] * fbe_0 * fz_0 + tk_yy_yzzzzz[i] * fe_0 + tk_xyy_yzzzzz[i] * pa_x[i] + 2.0 * ts_xxyy_yzzzzz[i] * fz_0;

        tk_xxyy_zzzzzz[i] =
            -2.0 * ts_yy_zzzzzz[i] * fbe_0 * fz_0 + tk_yy_zzzzzz[i] * fe_0 + tk_xyy_zzzzzz[i] * pa_x[i] + 2.0 * ts_xxyy_zzzzzz[i] * fz_0;
    }

    // Set up 112-140 components of targeted buffer : GI

    auto tk_xxyz_xxxxxx = pbuffer.data(idx_kin_gi + 112);

    auto tk_xxyz_xxxxxy = pbuffer.data(idx_kin_gi + 113);

    auto tk_xxyz_xxxxxz = pbuffer.data(idx_kin_gi + 114);

    auto tk_xxyz_xxxxyy = pbuffer.data(idx_kin_gi + 115);

    auto tk_xxyz_xxxxyz = pbuffer.data(idx_kin_gi + 116);

    auto tk_xxyz_xxxxzz = pbuffer.data(idx_kin_gi + 117);

    auto tk_xxyz_xxxyyy = pbuffer.data(idx_kin_gi + 118);

    auto tk_xxyz_xxxyyz = pbuffer.data(idx_kin_gi + 119);

    auto tk_xxyz_xxxyzz = pbuffer.data(idx_kin_gi + 120);

    auto tk_xxyz_xxxzzz = pbuffer.data(idx_kin_gi + 121);

    auto tk_xxyz_xxyyyy = pbuffer.data(idx_kin_gi + 122);

    auto tk_xxyz_xxyyyz = pbuffer.data(idx_kin_gi + 123);

    auto tk_xxyz_xxyyzz = pbuffer.data(idx_kin_gi + 124);

    auto tk_xxyz_xxyzzz = pbuffer.data(idx_kin_gi + 125);

    auto tk_xxyz_xxzzzz = pbuffer.data(idx_kin_gi + 126);

    auto tk_xxyz_xyyyyy = pbuffer.data(idx_kin_gi + 127);

    auto tk_xxyz_xyyyyz = pbuffer.data(idx_kin_gi + 128);

    auto tk_xxyz_xyyyzz = pbuffer.data(idx_kin_gi + 129);

    auto tk_xxyz_xyyzzz = pbuffer.data(idx_kin_gi + 130);

    auto tk_xxyz_xyzzzz = pbuffer.data(idx_kin_gi + 131);

    auto tk_xxyz_xzzzzz = pbuffer.data(idx_kin_gi + 132);

    auto tk_xxyz_yyyyyy = pbuffer.data(idx_kin_gi + 133);

    auto tk_xxyz_yyyyyz = pbuffer.data(idx_kin_gi + 134);

    auto tk_xxyz_yyyyzz = pbuffer.data(idx_kin_gi + 135);

    auto tk_xxyz_yyyzzz = pbuffer.data(idx_kin_gi + 136);

    auto tk_xxyz_yyzzzz = pbuffer.data(idx_kin_gi + 137);

    auto tk_xxyz_yzzzzz = pbuffer.data(idx_kin_gi + 138);

    auto tk_xxyz_zzzzzz = pbuffer.data(idx_kin_gi + 139);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tk_xxy_xxxxxy,  \
                             tk_xxy_xxxxyy,  \
                             tk_xxy_xxxyyy,  \
                             tk_xxy_xxyyyy,  \
                             tk_xxy_xyyyyy,  \
                             tk_xxy_yyyyyy,  \
                             tk_xxyz_xxxxxx, \
                             tk_xxyz_xxxxxy, \
                             tk_xxyz_xxxxxz, \
                             tk_xxyz_xxxxyy, \
                             tk_xxyz_xxxxyz, \
                             tk_xxyz_xxxxzz, \
                             tk_xxyz_xxxyyy, \
                             tk_xxyz_xxxyyz, \
                             tk_xxyz_xxxyzz, \
                             tk_xxyz_xxxzzz, \
                             tk_xxyz_xxyyyy, \
                             tk_xxyz_xxyyyz, \
                             tk_xxyz_xxyyzz, \
                             tk_xxyz_xxyzzz, \
                             tk_xxyz_xxzzzz, \
                             tk_xxyz_xyyyyy, \
                             tk_xxyz_xyyyyz, \
                             tk_xxyz_xyyyzz, \
                             tk_xxyz_xyyzzz, \
                             tk_xxyz_xyzzzz, \
                             tk_xxyz_xzzzzz, \
                             tk_xxyz_yyyyyy, \
                             tk_xxyz_yyyyyz, \
                             tk_xxyz_yyyyzz, \
                             tk_xxyz_yyyzzz, \
                             tk_xxyz_yyzzzz, \
                             tk_xxyz_yzzzzz, \
                             tk_xxyz_zzzzzz, \
                             tk_xxz_xxxxxx,  \
                             tk_xxz_xxxxxz,  \
                             tk_xxz_xxxxyz,  \
                             tk_xxz_xxxxz,   \
                             tk_xxz_xxxxzz,  \
                             tk_xxz_xxxyyz,  \
                             tk_xxz_xxxyz,   \
                             tk_xxz_xxxyzz,  \
                             tk_xxz_xxxzz,   \
                             tk_xxz_xxxzzz,  \
                             tk_xxz_xxyyyz,  \
                             tk_xxz_xxyyz,   \
                             tk_xxz_xxyyzz,  \
                             tk_xxz_xxyzz,   \
                             tk_xxz_xxyzzz,  \
                             tk_xxz_xxzzz,   \
                             tk_xxz_xxzzzz,  \
                             tk_xxz_xyyyyz,  \
                             tk_xxz_xyyyz,   \
                             tk_xxz_xyyyzz,  \
                             tk_xxz_xyyzz,   \
                             tk_xxz_xyyzzz,  \
                             tk_xxz_xyzzz,   \
                             tk_xxz_xyzzzz,  \
                             tk_xxz_xzzzz,   \
                             tk_xxz_xzzzzz,  \
                             tk_xxz_yyyyyz,  \
                             tk_xxz_yyyyz,   \
                             tk_xxz_yyyyzz,  \
                             tk_xxz_yyyzz,   \
                             tk_xxz_yyyzzz,  \
                             tk_xxz_yyzzz,   \
                             tk_xxz_yyzzzz,  \
                             tk_xxz_yzzzz,   \
                             tk_xxz_yzzzzz,  \
                             tk_xxz_zzzzz,   \
                             tk_xxz_zzzzzz,  \
                             ts_xxyz_xxxxxx, \
                             ts_xxyz_xxxxxy, \
                             ts_xxyz_xxxxxz, \
                             ts_xxyz_xxxxyy, \
                             ts_xxyz_xxxxyz, \
                             ts_xxyz_xxxxzz, \
                             ts_xxyz_xxxyyy, \
                             ts_xxyz_xxxyyz, \
                             ts_xxyz_xxxyzz, \
                             ts_xxyz_xxxzzz, \
                             ts_xxyz_xxyyyy, \
                             ts_xxyz_xxyyyz, \
                             ts_xxyz_xxyyzz, \
                             ts_xxyz_xxyzzz, \
                             ts_xxyz_xxzzzz, \
                             ts_xxyz_xyyyyy, \
                             ts_xxyz_xyyyyz, \
                             ts_xxyz_xyyyzz, \
                             ts_xxyz_xyyzzz, \
                             ts_xxyz_xyzzzz, \
                             ts_xxyz_xzzzzz, \
                             ts_xxyz_yyyyyy, \
                             ts_xxyz_yyyyyz, \
                             ts_xxyz_yyyyzz, \
                             ts_xxyz_yyyzzz, \
                             ts_xxyz_yyzzzz, \
                             ts_xxyz_yzzzzz, \
                             ts_xxyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyz_xxxxxx[i] = tk_xxz_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxyz_xxxxxx[i] * fz_0;

        tk_xxyz_xxxxxy[i] = tk_xxy_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxyz_xxxxxy[i] * fz_0;

        tk_xxyz_xxxxxz[i] = tk_xxz_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxyz_xxxxxz[i] * fz_0;

        tk_xxyz_xxxxyy[i] = tk_xxy_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxyz_xxxxyy[i] * fz_0;

        tk_xxyz_xxxxyz[i] = tk_xxz_xxxxz[i] * fe_0 + tk_xxz_xxxxyz[i] * pa_y[i] + 2.0 * ts_xxyz_xxxxyz[i] * fz_0;

        tk_xxyz_xxxxzz[i] = tk_xxz_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxyz_xxxxzz[i] * fz_0;

        tk_xxyz_xxxyyy[i] = tk_xxy_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxyz_xxxyyy[i] * fz_0;

        tk_xxyz_xxxyyz[i] = 2.0 * tk_xxz_xxxyz[i] * fe_0 + tk_xxz_xxxyyz[i] * pa_y[i] + 2.0 * ts_xxyz_xxxyyz[i] * fz_0;

        tk_xxyz_xxxyzz[i] = tk_xxz_xxxzz[i] * fe_0 + tk_xxz_xxxyzz[i] * pa_y[i] + 2.0 * ts_xxyz_xxxyzz[i] * fz_0;

        tk_xxyz_xxxzzz[i] = tk_xxz_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxyz_xxxzzz[i] * fz_0;

        tk_xxyz_xxyyyy[i] = tk_xxy_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxyz_xxyyyy[i] * fz_0;

        tk_xxyz_xxyyyz[i] = 3.0 * tk_xxz_xxyyz[i] * fe_0 + tk_xxz_xxyyyz[i] * pa_y[i] + 2.0 * ts_xxyz_xxyyyz[i] * fz_0;

        tk_xxyz_xxyyzz[i] = 2.0 * tk_xxz_xxyzz[i] * fe_0 + tk_xxz_xxyyzz[i] * pa_y[i] + 2.0 * ts_xxyz_xxyyzz[i] * fz_0;

        tk_xxyz_xxyzzz[i] = tk_xxz_xxzzz[i] * fe_0 + tk_xxz_xxyzzz[i] * pa_y[i] + 2.0 * ts_xxyz_xxyzzz[i] * fz_0;

        tk_xxyz_xxzzzz[i] = tk_xxz_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxyz_xxzzzz[i] * fz_0;

        tk_xxyz_xyyyyy[i] = tk_xxy_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxyz_xyyyyy[i] * fz_0;

        tk_xxyz_xyyyyz[i] = 4.0 * tk_xxz_xyyyz[i] * fe_0 + tk_xxz_xyyyyz[i] * pa_y[i] + 2.0 * ts_xxyz_xyyyyz[i] * fz_0;

        tk_xxyz_xyyyzz[i] = 3.0 * tk_xxz_xyyzz[i] * fe_0 + tk_xxz_xyyyzz[i] * pa_y[i] + 2.0 * ts_xxyz_xyyyzz[i] * fz_0;

        tk_xxyz_xyyzzz[i] = 2.0 * tk_xxz_xyzzz[i] * fe_0 + tk_xxz_xyyzzz[i] * pa_y[i] + 2.0 * ts_xxyz_xyyzzz[i] * fz_0;

        tk_xxyz_xyzzzz[i] = tk_xxz_xzzzz[i] * fe_0 + tk_xxz_xyzzzz[i] * pa_y[i] + 2.0 * ts_xxyz_xyzzzz[i] * fz_0;

        tk_xxyz_xzzzzz[i] = tk_xxz_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxyz_xzzzzz[i] * fz_0;

        tk_xxyz_yyyyyy[i] = tk_xxy_yyyyyy[i] * pa_z[i] + 2.0 * ts_xxyz_yyyyyy[i] * fz_0;

        tk_xxyz_yyyyyz[i] = 5.0 * tk_xxz_yyyyz[i] * fe_0 + tk_xxz_yyyyyz[i] * pa_y[i] + 2.0 * ts_xxyz_yyyyyz[i] * fz_0;

        tk_xxyz_yyyyzz[i] = 4.0 * tk_xxz_yyyzz[i] * fe_0 + tk_xxz_yyyyzz[i] * pa_y[i] + 2.0 * ts_xxyz_yyyyzz[i] * fz_0;

        tk_xxyz_yyyzzz[i] = 3.0 * tk_xxz_yyzzz[i] * fe_0 + tk_xxz_yyyzzz[i] * pa_y[i] + 2.0 * ts_xxyz_yyyzzz[i] * fz_0;

        tk_xxyz_yyzzzz[i] = 2.0 * tk_xxz_yzzzz[i] * fe_0 + tk_xxz_yyzzzz[i] * pa_y[i] + 2.0 * ts_xxyz_yyzzzz[i] * fz_0;

        tk_xxyz_yzzzzz[i] = tk_xxz_zzzzz[i] * fe_0 + tk_xxz_yzzzzz[i] * pa_y[i] + 2.0 * ts_xxyz_yzzzzz[i] * fz_0;

        tk_xxyz_zzzzzz[i] = tk_xxz_zzzzzz[i] * pa_y[i] + 2.0 * ts_xxyz_zzzzzz[i] * fz_0;
    }

    // Set up 140-168 components of targeted buffer : GI

    auto tk_xxzz_xxxxxx = pbuffer.data(idx_kin_gi + 140);

    auto tk_xxzz_xxxxxy = pbuffer.data(idx_kin_gi + 141);

    auto tk_xxzz_xxxxxz = pbuffer.data(idx_kin_gi + 142);

    auto tk_xxzz_xxxxyy = pbuffer.data(idx_kin_gi + 143);

    auto tk_xxzz_xxxxyz = pbuffer.data(idx_kin_gi + 144);

    auto tk_xxzz_xxxxzz = pbuffer.data(idx_kin_gi + 145);

    auto tk_xxzz_xxxyyy = pbuffer.data(idx_kin_gi + 146);

    auto tk_xxzz_xxxyyz = pbuffer.data(idx_kin_gi + 147);

    auto tk_xxzz_xxxyzz = pbuffer.data(idx_kin_gi + 148);

    auto tk_xxzz_xxxzzz = pbuffer.data(idx_kin_gi + 149);

    auto tk_xxzz_xxyyyy = pbuffer.data(idx_kin_gi + 150);

    auto tk_xxzz_xxyyyz = pbuffer.data(idx_kin_gi + 151);

    auto tk_xxzz_xxyyzz = pbuffer.data(idx_kin_gi + 152);

    auto tk_xxzz_xxyzzz = pbuffer.data(idx_kin_gi + 153);

    auto tk_xxzz_xxzzzz = pbuffer.data(idx_kin_gi + 154);

    auto tk_xxzz_xyyyyy = pbuffer.data(idx_kin_gi + 155);

    auto tk_xxzz_xyyyyz = pbuffer.data(idx_kin_gi + 156);

    auto tk_xxzz_xyyyzz = pbuffer.data(idx_kin_gi + 157);

    auto tk_xxzz_xyyzzz = pbuffer.data(idx_kin_gi + 158);

    auto tk_xxzz_xyzzzz = pbuffer.data(idx_kin_gi + 159);

    auto tk_xxzz_xzzzzz = pbuffer.data(idx_kin_gi + 160);

    auto tk_xxzz_yyyyyy = pbuffer.data(idx_kin_gi + 161);

    auto tk_xxzz_yyyyyz = pbuffer.data(idx_kin_gi + 162);

    auto tk_xxzz_yyyyzz = pbuffer.data(idx_kin_gi + 163);

    auto tk_xxzz_yyyzzz = pbuffer.data(idx_kin_gi + 164);

    auto tk_xxzz_yyzzzz = pbuffer.data(idx_kin_gi + 165);

    auto tk_xxzz_yzzzzz = pbuffer.data(idx_kin_gi + 166);

    auto tk_xxzz_zzzzzz = pbuffer.data(idx_kin_gi + 167);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tk_xx_xxxxxx,   \
                             tk_xx_xxxxxy,   \
                             tk_xx_xxxxyy,   \
                             tk_xx_xxxyyy,   \
                             tk_xx_xxyyyy,   \
                             tk_xx_xyyyyy,   \
                             tk_xxz_xxxxxx,  \
                             tk_xxz_xxxxxy,  \
                             tk_xxz_xxxxyy,  \
                             tk_xxz_xxxyyy,  \
                             tk_xxz_xxyyyy,  \
                             tk_xxz_xyyyyy,  \
                             tk_xxzz_xxxxxx, \
                             tk_xxzz_xxxxxy, \
                             tk_xxzz_xxxxxz, \
                             tk_xxzz_xxxxyy, \
                             tk_xxzz_xxxxyz, \
                             tk_xxzz_xxxxzz, \
                             tk_xxzz_xxxyyy, \
                             tk_xxzz_xxxyyz, \
                             tk_xxzz_xxxyzz, \
                             tk_xxzz_xxxzzz, \
                             tk_xxzz_xxyyyy, \
                             tk_xxzz_xxyyyz, \
                             tk_xxzz_xxyyzz, \
                             tk_xxzz_xxyzzz, \
                             tk_xxzz_xxzzzz, \
                             tk_xxzz_xyyyyy, \
                             tk_xxzz_xyyyyz, \
                             tk_xxzz_xyyyzz, \
                             tk_xxzz_xyyzzz, \
                             tk_xxzz_xyzzzz, \
                             tk_xxzz_xzzzzz, \
                             tk_xxzz_yyyyyy, \
                             tk_xxzz_yyyyyz, \
                             tk_xxzz_yyyyzz, \
                             tk_xxzz_yyyzzz, \
                             tk_xxzz_yyzzzz, \
                             tk_xxzz_yzzzzz, \
                             tk_xxzz_zzzzzz, \
                             tk_xzz_xxxxxz,  \
                             tk_xzz_xxxxyz,  \
                             tk_xzz_xxxxz,   \
                             tk_xzz_xxxxzz,  \
                             tk_xzz_xxxyyz,  \
                             tk_xzz_xxxyz,   \
                             tk_xzz_xxxyzz,  \
                             tk_xzz_xxxzz,   \
                             tk_xzz_xxxzzz,  \
                             tk_xzz_xxyyyz,  \
                             tk_xzz_xxyyz,   \
                             tk_xzz_xxyyzz,  \
                             tk_xzz_xxyzz,   \
                             tk_xzz_xxyzzz,  \
                             tk_xzz_xxzzz,   \
                             tk_xzz_xxzzzz,  \
                             tk_xzz_xyyyyz,  \
                             tk_xzz_xyyyz,   \
                             tk_xzz_xyyyzz,  \
                             tk_xzz_xyyzz,   \
                             tk_xzz_xyyzzz,  \
                             tk_xzz_xyzzz,   \
                             tk_xzz_xyzzzz,  \
                             tk_xzz_xzzzz,   \
                             tk_xzz_xzzzzz,  \
                             tk_xzz_yyyyyy,  \
                             tk_xzz_yyyyyz,  \
                             tk_xzz_yyyyz,   \
                             tk_xzz_yyyyzz,  \
                             tk_xzz_yyyzz,   \
                             tk_xzz_yyyzzz,  \
                             tk_xzz_yyzzz,   \
                             tk_xzz_yyzzzz,  \
                             tk_xzz_yzzzz,   \
                             tk_xzz_yzzzzz,  \
                             tk_xzz_zzzzz,   \
                             tk_xzz_zzzzzz,  \
                             tk_zz_xxxxxz,   \
                             tk_zz_xxxxyz,   \
                             tk_zz_xxxxzz,   \
                             tk_zz_xxxyyz,   \
                             tk_zz_xxxyzz,   \
                             tk_zz_xxxzzz,   \
                             tk_zz_xxyyyz,   \
                             tk_zz_xxyyzz,   \
                             tk_zz_xxyzzz,   \
                             tk_zz_xxzzzz,   \
                             tk_zz_xyyyyz,   \
                             tk_zz_xyyyzz,   \
                             tk_zz_xyyzzz,   \
                             tk_zz_xyzzzz,   \
                             tk_zz_xzzzzz,   \
                             tk_zz_yyyyyy,   \
                             tk_zz_yyyyyz,   \
                             tk_zz_yyyyzz,   \
                             tk_zz_yyyzzz,   \
                             tk_zz_yyzzzz,   \
                             tk_zz_yzzzzz,   \
                             tk_zz_zzzzzz,   \
                             ts_xx_xxxxxx,   \
                             ts_xx_xxxxxy,   \
                             ts_xx_xxxxyy,   \
                             ts_xx_xxxyyy,   \
                             ts_xx_xxyyyy,   \
                             ts_xx_xyyyyy,   \
                             ts_xxzz_xxxxxx, \
                             ts_xxzz_xxxxxy, \
                             ts_xxzz_xxxxxz, \
                             ts_xxzz_xxxxyy, \
                             ts_xxzz_xxxxyz, \
                             ts_xxzz_xxxxzz, \
                             ts_xxzz_xxxyyy, \
                             ts_xxzz_xxxyyz, \
                             ts_xxzz_xxxyzz, \
                             ts_xxzz_xxxzzz, \
                             ts_xxzz_xxyyyy, \
                             ts_xxzz_xxyyyz, \
                             ts_xxzz_xxyyzz, \
                             ts_xxzz_xxyzzz, \
                             ts_xxzz_xxzzzz, \
                             ts_xxzz_xyyyyy, \
                             ts_xxzz_xyyyyz, \
                             ts_xxzz_xyyyzz, \
                             ts_xxzz_xyyzzz, \
                             ts_xxzz_xyzzzz, \
                             ts_xxzz_xzzzzz, \
                             ts_xxzz_yyyyyy, \
                             ts_xxzz_yyyyyz, \
                             ts_xxzz_yyyyzz, \
                             ts_xxzz_yyyzzz, \
                             ts_xxzz_yyzzzz, \
                             ts_xxzz_yzzzzz, \
                             ts_xxzz_zzzzzz, \
                             ts_zz_xxxxxz,   \
                             ts_zz_xxxxyz,   \
                             ts_zz_xxxxzz,   \
                             ts_zz_xxxyyz,   \
                             ts_zz_xxxyzz,   \
                             ts_zz_xxxzzz,   \
                             ts_zz_xxyyyz,   \
                             ts_zz_xxyyzz,   \
                             ts_zz_xxyzzz,   \
                             ts_zz_xxzzzz,   \
                             ts_zz_xyyyyz,   \
                             ts_zz_xyyyzz,   \
                             ts_zz_xyyzzz,   \
                             ts_zz_xyzzzz,   \
                             ts_zz_xzzzzz,   \
                             ts_zz_yyyyyy,   \
                             ts_zz_yyyyyz,   \
                             ts_zz_yyyyzz,   \
                             ts_zz_yyyzzz,   \
                             ts_zz_yyzzzz,   \
                             ts_zz_yzzzzz,   \
                             ts_zz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzz_xxxxxx[i] =
            -2.0 * ts_xx_xxxxxx[i] * fbe_0 * fz_0 + tk_xx_xxxxxx[i] * fe_0 + tk_xxz_xxxxxx[i] * pa_z[i] + 2.0 * ts_xxzz_xxxxxx[i] * fz_0;

        tk_xxzz_xxxxxy[i] =
            -2.0 * ts_xx_xxxxxy[i] * fbe_0 * fz_0 + tk_xx_xxxxxy[i] * fe_0 + tk_xxz_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxzz_xxxxxy[i] * fz_0;

        tk_xxzz_xxxxxz[i] = -2.0 * ts_zz_xxxxxz[i] * fbe_0 * fz_0 + tk_zz_xxxxxz[i] * fe_0 + 5.0 * tk_xzz_xxxxz[i] * fe_0 +
                            tk_xzz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xxzz_xxxxxz[i] * fz_0;

        tk_xxzz_xxxxyy[i] =
            -2.0 * ts_xx_xxxxyy[i] * fbe_0 * fz_0 + tk_xx_xxxxyy[i] * fe_0 + tk_xxz_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxzz_xxxxyy[i] * fz_0;

        tk_xxzz_xxxxyz[i] = -2.0 * ts_zz_xxxxyz[i] * fbe_0 * fz_0 + tk_zz_xxxxyz[i] * fe_0 + 4.0 * tk_xzz_xxxyz[i] * fe_0 +
                            tk_xzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxzz_xxxxyz[i] * fz_0;

        tk_xxzz_xxxxzz[i] = -2.0 * ts_zz_xxxxzz[i] * fbe_0 * fz_0 + tk_zz_xxxxzz[i] * fe_0 + 4.0 * tk_xzz_xxxzz[i] * fe_0 +
                            tk_xzz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xxzz_xxxxzz[i] * fz_0;

        tk_xxzz_xxxyyy[i] =
            -2.0 * ts_xx_xxxyyy[i] * fbe_0 * fz_0 + tk_xx_xxxyyy[i] * fe_0 + tk_xxz_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxzz_xxxyyy[i] * fz_0;

        tk_xxzz_xxxyyz[i] = -2.0 * ts_zz_xxxyyz[i] * fbe_0 * fz_0 + tk_zz_xxxyyz[i] * fe_0 + 3.0 * tk_xzz_xxyyz[i] * fe_0 +
                            tk_xzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxzz_xxxyyz[i] * fz_0;

        tk_xxzz_xxxyzz[i] = -2.0 * ts_zz_xxxyzz[i] * fbe_0 * fz_0 + tk_zz_xxxyzz[i] * fe_0 + 3.0 * tk_xzz_xxyzz[i] * fe_0 +
                            tk_xzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxzz_xxxyzz[i] * fz_0;

        tk_xxzz_xxxzzz[i] = -2.0 * ts_zz_xxxzzz[i] * fbe_0 * fz_0 + tk_zz_xxxzzz[i] * fe_0 + 3.0 * tk_xzz_xxzzz[i] * fe_0 +
                            tk_xzz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xxzz_xxxzzz[i] * fz_0;

        tk_xxzz_xxyyyy[i] =
            -2.0 * ts_xx_xxyyyy[i] * fbe_0 * fz_0 + tk_xx_xxyyyy[i] * fe_0 + tk_xxz_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxzz_xxyyyy[i] * fz_0;

        tk_xxzz_xxyyyz[i] = -2.0 * ts_zz_xxyyyz[i] * fbe_0 * fz_0 + tk_zz_xxyyyz[i] * fe_0 + 2.0 * tk_xzz_xyyyz[i] * fe_0 +
                            tk_xzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxzz_xxyyyz[i] * fz_0;

        tk_xxzz_xxyyzz[i] = -2.0 * ts_zz_xxyyzz[i] * fbe_0 * fz_0 + tk_zz_xxyyzz[i] * fe_0 + 2.0 * tk_xzz_xyyzz[i] * fe_0 +
                            tk_xzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxzz_xxyyzz[i] * fz_0;

        tk_xxzz_xxyzzz[i] = -2.0 * ts_zz_xxyzzz[i] * fbe_0 * fz_0 + tk_zz_xxyzzz[i] * fe_0 + 2.0 * tk_xzz_xyzzz[i] * fe_0 +
                            tk_xzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxzz_xxyzzz[i] * fz_0;

        tk_xxzz_xxzzzz[i] = -2.0 * ts_zz_xxzzzz[i] * fbe_0 * fz_0 + tk_zz_xxzzzz[i] * fe_0 + 2.0 * tk_xzz_xzzzz[i] * fe_0 +
                            tk_xzz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xxzz_xxzzzz[i] * fz_0;

        tk_xxzz_xyyyyy[i] =
            -2.0 * ts_xx_xyyyyy[i] * fbe_0 * fz_0 + tk_xx_xyyyyy[i] * fe_0 + tk_xxz_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxzz_xyyyyy[i] * fz_0;

        tk_xxzz_xyyyyz[i] = -2.0 * ts_zz_xyyyyz[i] * fbe_0 * fz_0 + tk_zz_xyyyyz[i] * fe_0 + tk_xzz_yyyyz[i] * fe_0 + tk_xzz_xyyyyz[i] * pa_x[i] +
                            2.0 * ts_xxzz_xyyyyz[i] * fz_0;

        tk_xxzz_xyyyzz[i] = -2.0 * ts_zz_xyyyzz[i] * fbe_0 * fz_0 + tk_zz_xyyyzz[i] * fe_0 + tk_xzz_yyyzz[i] * fe_0 + tk_xzz_xyyyzz[i] * pa_x[i] +
                            2.0 * ts_xxzz_xyyyzz[i] * fz_0;

        tk_xxzz_xyyzzz[i] = -2.0 * ts_zz_xyyzzz[i] * fbe_0 * fz_0 + tk_zz_xyyzzz[i] * fe_0 + tk_xzz_yyzzz[i] * fe_0 + tk_xzz_xyyzzz[i] * pa_x[i] +
                            2.0 * ts_xxzz_xyyzzz[i] * fz_0;

        tk_xxzz_xyzzzz[i] = -2.0 * ts_zz_xyzzzz[i] * fbe_0 * fz_0 + tk_zz_xyzzzz[i] * fe_0 + tk_xzz_yzzzz[i] * fe_0 + tk_xzz_xyzzzz[i] * pa_x[i] +
                            2.0 * ts_xxzz_xyzzzz[i] * fz_0;

        tk_xxzz_xzzzzz[i] = -2.0 * ts_zz_xzzzzz[i] * fbe_0 * fz_0 + tk_zz_xzzzzz[i] * fe_0 + tk_xzz_zzzzz[i] * fe_0 + tk_xzz_xzzzzz[i] * pa_x[i] +
                            2.0 * ts_xxzz_xzzzzz[i] * fz_0;

        tk_xxzz_yyyyyy[i] =
            -2.0 * ts_zz_yyyyyy[i] * fbe_0 * fz_0 + tk_zz_yyyyyy[i] * fe_0 + tk_xzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xxzz_yyyyyy[i] * fz_0;

        tk_xxzz_yyyyyz[i] =
            -2.0 * ts_zz_yyyyyz[i] * fbe_0 * fz_0 + tk_zz_yyyyyz[i] * fe_0 + tk_xzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xxzz_yyyyyz[i] * fz_0;

        tk_xxzz_yyyyzz[i] =
            -2.0 * ts_zz_yyyyzz[i] * fbe_0 * fz_0 + tk_zz_yyyyzz[i] * fe_0 + tk_xzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xxzz_yyyyzz[i] * fz_0;

        tk_xxzz_yyyzzz[i] =
            -2.0 * ts_zz_yyyzzz[i] * fbe_0 * fz_0 + tk_zz_yyyzzz[i] * fe_0 + tk_xzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xxzz_yyyzzz[i] * fz_0;

        tk_xxzz_yyzzzz[i] =
            -2.0 * ts_zz_yyzzzz[i] * fbe_0 * fz_0 + tk_zz_yyzzzz[i] * fe_0 + tk_xzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xxzz_yyzzzz[i] * fz_0;

        tk_xxzz_yzzzzz[i] =
            -2.0 * ts_zz_yzzzzz[i] * fbe_0 * fz_0 + tk_zz_yzzzzz[i] * fe_0 + tk_xzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xxzz_yzzzzz[i] * fz_0;

        tk_xxzz_zzzzzz[i] =
            -2.0 * ts_zz_zzzzzz[i] * fbe_0 * fz_0 + tk_zz_zzzzzz[i] * fe_0 + tk_xzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xxzz_zzzzzz[i] * fz_0;
    }

    // Set up 168-196 components of targeted buffer : GI

    auto tk_xyyy_xxxxxx = pbuffer.data(idx_kin_gi + 168);

    auto tk_xyyy_xxxxxy = pbuffer.data(idx_kin_gi + 169);

    auto tk_xyyy_xxxxxz = pbuffer.data(idx_kin_gi + 170);

    auto tk_xyyy_xxxxyy = pbuffer.data(idx_kin_gi + 171);

    auto tk_xyyy_xxxxyz = pbuffer.data(idx_kin_gi + 172);

    auto tk_xyyy_xxxxzz = pbuffer.data(idx_kin_gi + 173);

    auto tk_xyyy_xxxyyy = pbuffer.data(idx_kin_gi + 174);

    auto tk_xyyy_xxxyyz = pbuffer.data(idx_kin_gi + 175);

    auto tk_xyyy_xxxyzz = pbuffer.data(idx_kin_gi + 176);

    auto tk_xyyy_xxxzzz = pbuffer.data(idx_kin_gi + 177);

    auto tk_xyyy_xxyyyy = pbuffer.data(idx_kin_gi + 178);

    auto tk_xyyy_xxyyyz = pbuffer.data(idx_kin_gi + 179);

    auto tk_xyyy_xxyyzz = pbuffer.data(idx_kin_gi + 180);

    auto tk_xyyy_xxyzzz = pbuffer.data(idx_kin_gi + 181);

    auto tk_xyyy_xxzzzz = pbuffer.data(idx_kin_gi + 182);

    auto tk_xyyy_xyyyyy = pbuffer.data(idx_kin_gi + 183);

    auto tk_xyyy_xyyyyz = pbuffer.data(idx_kin_gi + 184);

    auto tk_xyyy_xyyyzz = pbuffer.data(idx_kin_gi + 185);

    auto tk_xyyy_xyyzzz = pbuffer.data(idx_kin_gi + 186);

    auto tk_xyyy_xyzzzz = pbuffer.data(idx_kin_gi + 187);

    auto tk_xyyy_xzzzzz = pbuffer.data(idx_kin_gi + 188);

    auto tk_xyyy_yyyyyy = pbuffer.data(idx_kin_gi + 189);

    auto tk_xyyy_yyyyyz = pbuffer.data(idx_kin_gi + 190);

    auto tk_xyyy_yyyyzz = pbuffer.data(idx_kin_gi + 191);

    auto tk_xyyy_yyyzzz = pbuffer.data(idx_kin_gi + 192);

    auto tk_xyyy_yyzzzz = pbuffer.data(idx_kin_gi + 193);

    auto tk_xyyy_yzzzzz = pbuffer.data(idx_kin_gi + 194);

    auto tk_xyyy_zzzzzz = pbuffer.data(idx_kin_gi + 195);

#pragma omp simd aligned(pa_x,               \
                             tk_xyyy_xxxxxx, \
                             tk_xyyy_xxxxxy, \
                             tk_xyyy_xxxxxz, \
                             tk_xyyy_xxxxyy, \
                             tk_xyyy_xxxxyz, \
                             tk_xyyy_xxxxzz, \
                             tk_xyyy_xxxyyy, \
                             tk_xyyy_xxxyyz, \
                             tk_xyyy_xxxyzz, \
                             tk_xyyy_xxxzzz, \
                             tk_xyyy_xxyyyy, \
                             tk_xyyy_xxyyyz, \
                             tk_xyyy_xxyyzz, \
                             tk_xyyy_xxyzzz, \
                             tk_xyyy_xxzzzz, \
                             tk_xyyy_xyyyyy, \
                             tk_xyyy_xyyyyz, \
                             tk_xyyy_xyyyzz, \
                             tk_xyyy_xyyzzz, \
                             tk_xyyy_xyzzzz, \
                             tk_xyyy_xzzzzz, \
                             tk_xyyy_yyyyyy, \
                             tk_xyyy_yyyyyz, \
                             tk_xyyy_yyyyzz, \
                             tk_xyyy_yyyzzz, \
                             tk_xyyy_yyzzzz, \
                             tk_xyyy_yzzzzz, \
                             tk_xyyy_zzzzzz, \
                             tk_yyy_xxxxx,   \
                             tk_yyy_xxxxxx,  \
                             tk_yyy_xxxxxy,  \
                             tk_yyy_xxxxxz,  \
                             tk_yyy_xxxxy,   \
                             tk_yyy_xxxxyy,  \
                             tk_yyy_xxxxyz,  \
                             tk_yyy_xxxxz,   \
                             tk_yyy_xxxxzz,  \
                             tk_yyy_xxxyy,   \
                             tk_yyy_xxxyyy,  \
                             tk_yyy_xxxyyz,  \
                             tk_yyy_xxxyz,   \
                             tk_yyy_xxxyzz,  \
                             tk_yyy_xxxzz,   \
                             tk_yyy_xxxzzz,  \
                             tk_yyy_xxyyy,   \
                             tk_yyy_xxyyyy,  \
                             tk_yyy_xxyyyz,  \
                             tk_yyy_xxyyz,   \
                             tk_yyy_xxyyzz,  \
                             tk_yyy_xxyzz,   \
                             tk_yyy_xxyzzz,  \
                             tk_yyy_xxzzz,   \
                             tk_yyy_xxzzzz,  \
                             tk_yyy_xyyyy,   \
                             tk_yyy_xyyyyy,  \
                             tk_yyy_xyyyyz,  \
                             tk_yyy_xyyyz,   \
                             tk_yyy_xyyyzz,  \
                             tk_yyy_xyyzz,   \
                             tk_yyy_xyyzzz,  \
                             tk_yyy_xyzzz,   \
                             tk_yyy_xyzzzz,  \
                             tk_yyy_xzzzz,   \
                             tk_yyy_xzzzzz,  \
                             tk_yyy_yyyyy,   \
                             tk_yyy_yyyyyy,  \
                             tk_yyy_yyyyyz,  \
                             tk_yyy_yyyyz,   \
                             tk_yyy_yyyyzz,  \
                             tk_yyy_yyyzz,   \
                             tk_yyy_yyyzzz,  \
                             tk_yyy_yyzzz,   \
                             tk_yyy_yyzzzz,  \
                             tk_yyy_yzzzz,   \
                             tk_yyy_yzzzzz,  \
                             tk_yyy_zzzzz,   \
                             tk_yyy_zzzzzz,  \
                             ts_xyyy_xxxxxx, \
                             ts_xyyy_xxxxxy, \
                             ts_xyyy_xxxxxz, \
                             ts_xyyy_xxxxyy, \
                             ts_xyyy_xxxxyz, \
                             ts_xyyy_xxxxzz, \
                             ts_xyyy_xxxyyy, \
                             ts_xyyy_xxxyyz, \
                             ts_xyyy_xxxyzz, \
                             ts_xyyy_xxxzzz, \
                             ts_xyyy_xxyyyy, \
                             ts_xyyy_xxyyyz, \
                             ts_xyyy_xxyyzz, \
                             ts_xyyy_xxyzzz, \
                             ts_xyyy_xxzzzz, \
                             ts_xyyy_xyyyyy, \
                             ts_xyyy_xyyyyz, \
                             ts_xyyy_xyyyzz, \
                             ts_xyyy_xyyzzz, \
                             ts_xyyy_xyzzzz, \
                             ts_xyyy_xzzzzz, \
                             ts_xyyy_yyyyyy, \
                             ts_xyyy_yyyyyz, \
                             ts_xyyy_yyyyzz, \
                             ts_xyyy_yyyzzz, \
                             ts_xyyy_yyzzzz, \
                             ts_xyyy_yzzzzz, \
                             ts_xyyy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyy_xxxxxx[i] = 6.0 * tk_yyy_xxxxx[i] * fe_0 + tk_yyy_xxxxxx[i] * pa_x[i] + 2.0 * ts_xyyy_xxxxxx[i] * fz_0;

        tk_xyyy_xxxxxy[i] = 5.0 * tk_yyy_xxxxy[i] * fe_0 + tk_yyy_xxxxxy[i] * pa_x[i] + 2.0 * ts_xyyy_xxxxxy[i] * fz_0;

        tk_xyyy_xxxxxz[i] = 5.0 * tk_yyy_xxxxz[i] * fe_0 + tk_yyy_xxxxxz[i] * pa_x[i] + 2.0 * ts_xyyy_xxxxxz[i] * fz_0;

        tk_xyyy_xxxxyy[i] = 4.0 * tk_yyy_xxxyy[i] * fe_0 + tk_yyy_xxxxyy[i] * pa_x[i] + 2.0 * ts_xyyy_xxxxyy[i] * fz_0;

        tk_xyyy_xxxxyz[i] = 4.0 * tk_yyy_xxxyz[i] * fe_0 + tk_yyy_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyyy_xxxxyz[i] * fz_0;

        tk_xyyy_xxxxzz[i] = 4.0 * tk_yyy_xxxzz[i] * fe_0 + tk_yyy_xxxxzz[i] * pa_x[i] + 2.0 * ts_xyyy_xxxxzz[i] * fz_0;

        tk_xyyy_xxxyyy[i] = 3.0 * tk_yyy_xxyyy[i] * fe_0 + tk_yyy_xxxyyy[i] * pa_x[i] + 2.0 * ts_xyyy_xxxyyy[i] * fz_0;

        tk_xyyy_xxxyyz[i] = 3.0 * tk_yyy_xxyyz[i] * fe_0 + tk_yyy_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyyy_xxxyyz[i] * fz_0;

        tk_xyyy_xxxyzz[i] = 3.0 * tk_yyy_xxyzz[i] * fe_0 + tk_yyy_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyyy_xxxyzz[i] * fz_0;

        tk_xyyy_xxxzzz[i] = 3.0 * tk_yyy_xxzzz[i] * fe_0 + tk_yyy_xxxzzz[i] * pa_x[i] + 2.0 * ts_xyyy_xxxzzz[i] * fz_0;

        tk_xyyy_xxyyyy[i] = 2.0 * tk_yyy_xyyyy[i] * fe_0 + tk_yyy_xxyyyy[i] * pa_x[i] + 2.0 * ts_xyyy_xxyyyy[i] * fz_0;

        tk_xyyy_xxyyyz[i] = 2.0 * tk_yyy_xyyyz[i] * fe_0 + tk_yyy_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyyy_xxyyyz[i] * fz_0;

        tk_xyyy_xxyyzz[i] = 2.0 * tk_yyy_xyyzz[i] * fe_0 + tk_yyy_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyyy_xxyyzz[i] * fz_0;

        tk_xyyy_xxyzzz[i] = 2.0 * tk_yyy_xyzzz[i] * fe_0 + tk_yyy_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyyy_xxyzzz[i] * fz_0;

        tk_xyyy_xxzzzz[i] = 2.0 * tk_yyy_xzzzz[i] * fe_0 + tk_yyy_xxzzzz[i] * pa_x[i] + 2.0 * ts_xyyy_xxzzzz[i] * fz_0;

        tk_xyyy_xyyyyy[i] = tk_yyy_yyyyy[i] * fe_0 + tk_yyy_xyyyyy[i] * pa_x[i] + 2.0 * ts_xyyy_xyyyyy[i] * fz_0;

        tk_xyyy_xyyyyz[i] = tk_yyy_yyyyz[i] * fe_0 + tk_yyy_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyyy_xyyyyz[i] * fz_0;

        tk_xyyy_xyyyzz[i] = tk_yyy_yyyzz[i] * fe_0 + tk_yyy_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyyy_xyyyzz[i] * fz_0;

        tk_xyyy_xyyzzz[i] = tk_yyy_yyzzz[i] * fe_0 + tk_yyy_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyyy_xyyzzz[i] * fz_0;

        tk_xyyy_xyzzzz[i] = tk_yyy_yzzzz[i] * fe_0 + tk_yyy_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyyy_xyzzzz[i] * fz_0;

        tk_xyyy_xzzzzz[i] = tk_yyy_zzzzz[i] * fe_0 + tk_yyy_xzzzzz[i] * pa_x[i] + 2.0 * ts_xyyy_xzzzzz[i] * fz_0;

        tk_xyyy_yyyyyy[i] = tk_yyy_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyyy_yyyyyy[i] * fz_0;

        tk_xyyy_yyyyyz[i] = tk_yyy_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyyy_yyyyyz[i] * fz_0;

        tk_xyyy_yyyyzz[i] = tk_yyy_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyyy_yyyyzz[i] * fz_0;

        tk_xyyy_yyyzzz[i] = tk_yyy_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyyy_yyyzzz[i] * fz_0;

        tk_xyyy_yyzzzz[i] = tk_yyy_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyyy_yyzzzz[i] * fz_0;

        tk_xyyy_yzzzzz[i] = tk_yyy_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyyy_yzzzzz[i] * fz_0;

        tk_xyyy_zzzzzz[i] = tk_yyy_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyyy_zzzzzz[i] * fz_0;
    }

    // Set up 196-224 components of targeted buffer : GI

    auto tk_xyyz_xxxxxx = pbuffer.data(idx_kin_gi + 196);

    auto tk_xyyz_xxxxxy = pbuffer.data(idx_kin_gi + 197);

    auto tk_xyyz_xxxxxz = pbuffer.data(idx_kin_gi + 198);

    auto tk_xyyz_xxxxyy = pbuffer.data(idx_kin_gi + 199);

    auto tk_xyyz_xxxxyz = pbuffer.data(idx_kin_gi + 200);

    auto tk_xyyz_xxxxzz = pbuffer.data(idx_kin_gi + 201);

    auto tk_xyyz_xxxyyy = pbuffer.data(idx_kin_gi + 202);

    auto tk_xyyz_xxxyyz = pbuffer.data(idx_kin_gi + 203);

    auto tk_xyyz_xxxyzz = pbuffer.data(idx_kin_gi + 204);

    auto tk_xyyz_xxxzzz = pbuffer.data(idx_kin_gi + 205);

    auto tk_xyyz_xxyyyy = pbuffer.data(idx_kin_gi + 206);

    auto tk_xyyz_xxyyyz = pbuffer.data(idx_kin_gi + 207);

    auto tk_xyyz_xxyyzz = pbuffer.data(idx_kin_gi + 208);

    auto tk_xyyz_xxyzzz = pbuffer.data(idx_kin_gi + 209);

    auto tk_xyyz_xxzzzz = pbuffer.data(idx_kin_gi + 210);

    auto tk_xyyz_xyyyyy = pbuffer.data(idx_kin_gi + 211);

    auto tk_xyyz_xyyyyz = pbuffer.data(idx_kin_gi + 212);

    auto tk_xyyz_xyyyzz = pbuffer.data(idx_kin_gi + 213);

    auto tk_xyyz_xyyzzz = pbuffer.data(idx_kin_gi + 214);

    auto tk_xyyz_xyzzzz = pbuffer.data(idx_kin_gi + 215);

    auto tk_xyyz_xzzzzz = pbuffer.data(idx_kin_gi + 216);

    auto tk_xyyz_yyyyyy = pbuffer.data(idx_kin_gi + 217);

    auto tk_xyyz_yyyyyz = pbuffer.data(idx_kin_gi + 218);

    auto tk_xyyz_yyyyzz = pbuffer.data(idx_kin_gi + 219);

    auto tk_xyyz_yyyzzz = pbuffer.data(idx_kin_gi + 220);

    auto tk_xyyz_yyzzzz = pbuffer.data(idx_kin_gi + 221);

    auto tk_xyyz_yzzzzz = pbuffer.data(idx_kin_gi + 222);

    auto tk_xyyz_zzzzzz = pbuffer.data(idx_kin_gi + 223);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tk_xyy_xxxxxx,  \
                             tk_xyy_xxxxxy,  \
                             tk_xyy_xxxxyy,  \
                             tk_xyy_xxxyyy,  \
                             tk_xyy_xxyyyy,  \
                             tk_xyy_xyyyyy,  \
                             tk_xyyz_xxxxxx, \
                             tk_xyyz_xxxxxy, \
                             tk_xyyz_xxxxxz, \
                             tk_xyyz_xxxxyy, \
                             tk_xyyz_xxxxyz, \
                             tk_xyyz_xxxxzz, \
                             tk_xyyz_xxxyyy, \
                             tk_xyyz_xxxyyz, \
                             tk_xyyz_xxxyzz, \
                             tk_xyyz_xxxzzz, \
                             tk_xyyz_xxyyyy, \
                             tk_xyyz_xxyyyz, \
                             tk_xyyz_xxyyzz, \
                             tk_xyyz_xxyzzz, \
                             tk_xyyz_xxzzzz, \
                             tk_xyyz_xyyyyy, \
                             tk_xyyz_xyyyyz, \
                             tk_xyyz_xyyyzz, \
                             tk_xyyz_xyyzzz, \
                             tk_xyyz_xyzzzz, \
                             tk_xyyz_xzzzzz, \
                             tk_xyyz_yyyyyy, \
                             tk_xyyz_yyyyyz, \
                             tk_xyyz_yyyyzz, \
                             tk_xyyz_yyyzzz, \
                             tk_xyyz_yyzzzz, \
                             tk_xyyz_yzzzzz, \
                             tk_xyyz_zzzzzz, \
                             tk_yyz_xxxxxz,  \
                             tk_yyz_xxxxyz,  \
                             tk_yyz_xxxxz,   \
                             tk_yyz_xxxxzz,  \
                             tk_yyz_xxxyyz,  \
                             tk_yyz_xxxyz,   \
                             tk_yyz_xxxyzz,  \
                             tk_yyz_xxxzz,   \
                             tk_yyz_xxxzzz,  \
                             tk_yyz_xxyyyz,  \
                             tk_yyz_xxyyz,   \
                             tk_yyz_xxyyzz,  \
                             tk_yyz_xxyzz,   \
                             tk_yyz_xxyzzz,  \
                             tk_yyz_xxzzz,   \
                             tk_yyz_xxzzzz,  \
                             tk_yyz_xyyyyz,  \
                             tk_yyz_xyyyz,   \
                             tk_yyz_xyyyzz,  \
                             tk_yyz_xyyzz,   \
                             tk_yyz_xyyzzz,  \
                             tk_yyz_xyzzz,   \
                             tk_yyz_xyzzzz,  \
                             tk_yyz_xzzzz,   \
                             tk_yyz_xzzzzz,  \
                             tk_yyz_yyyyyy,  \
                             tk_yyz_yyyyyz,  \
                             tk_yyz_yyyyz,   \
                             tk_yyz_yyyyzz,  \
                             tk_yyz_yyyzz,   \
                             tk_yyz_yyyzzz,  \
                             tk_yyz_yyzzz,   \
                             tk_yyz_yyzzzz,  \
                             tk_yyz_yzzzz,   \
                             tk_yyz_yzzzzz,  \
                             tk_yyz_zzzzz,   \
                             tk_yyz_zzzzzz,  \
                             ts_xyyz_xxxxxx, \
                             ts_xyyz_xxxxxy, \
                             ts_xyyz_xxxxxz, \
                             ts_xyyz_xxxxyy, \
                             ts_xyyz_xxxxyz, \
                             ts_xyyz_xxxxzz, \
                             ts_xyyz_xxxyyy, \
                             ts_xyyz_xxxyyz, \
                             ts_xyyz_xxxyzz, \
                             ts_xyyz_xxxzzz, \
                             ts_xyyz_xxyyyy, \
                             ts_xyyz_xxyyyz, \
                             ts_xyyz_xxyyzz, \
                             ts_xyyz_xxyzzz, \
                             ts_xyyz_xxzzzz, \
                             ts_xyyz_xyyyyy, \
                             ts_xyyz_xyyyyz, \
                             ts_xyyz_xyyyzz, \
                             ts_xyyz_xyyzzz, \
                             ts_xyyz_xyzzzz, \
                             ts_xyyz_xzzzzz, \
                             ts_xyyz_yyyyyy, \
                             ts_xyyz_yyyyyz, \
                             ts_xyyz_yyyyzz, \
                             ts_xyyz_yyyzzz, \
                             ts_xyyz_yyzzzz, \
                             ts_xyyz_yzzzzz, \
                             ts_xyyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyz_xxxxxx[i] = tk_xyy_xxxxxx[i] * pa_z[i] + 2.0 * ts_xyyz_xxxxxx[i] * fz_0;

        tk_xyyz_xxxxxy[i] = tk_xyy_xxxxxy[i] * pa_z[i] + 2.0 * ts_xyyz_xxxxxy[i] * fz_0;

        tk_xyyz_xxxxxz[i] = 5.0 * tk_yyz_xxxxz[i] * fe_0 + tk_yyz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xyyz_xxxxxz[i] * fz_0;

        tk_xyyz_xxxxyy[i] = tk_xyy_xxxxyy[i] * pa_z[i] + 2.0 * ts_xyyz_xxxxyy[i] * fz_0;

        tk_xyyz_xxxxyz[i] = 4.0 * tk_yyz_xxxyz[i] * fe_0 + tk_yyz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyyz_xxxxyz[i] * fz_0;

        tk_xyyz_xxxxzz[i] = 4.0 * tk_yyz_xxxzz[i] * fe_0 + tk_yyz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xyyz_xxxxzz[i] * fz_0;

        tk_xyyz_xxxyyy[i] = tk_xyy_xxxyyy[i] * pa_z[i] + 2.0 * ts_xyyz_xxxyyy[i] * fz_0;

        tk_xyyz_xxxyyz[i] = 3.0 * tk_yyz_xxyyz[i] * fe_0 + tk_yyz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyyz_xxxyyz[i] * fz_0;

        tk_xyyz_xxxyzz[i] = 3.0 * tk_yyz_xxyzz[i] * fe_0 + tk_yyz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyyz_xxxyzz[i] * fz_0;

        tk_xyyz_xxxzzz[i] = 3.0 * tk_yyz_xxzzz[i] * fe_0 + tk_yyz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xyyz_xxxzzz[i] * fz_0;

        tk_xyyz_xxyyyy[i] = tk_xyy_xxyyyy[i] * pa_z[i] + 2.0 * ts_xyyz_xxyyyy[i] * fz_0;

        tk_xyyz_xxyyyz[i] = 2.0 * tk_yyz_xyyyz[i] * fe_0 + tk_yyz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyyz_xxyyyz[i] * fz_0;

        tk_xyyz_xxyyzz[i] = 2.0 * tk_yyz_xyyzz[i] * fe_0 + tk_yyz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyyz_xxyyzz[i] * fz_0;

        tk_xyyz_xxyzzz[i] = 2.0 * tk_yyz_xyzzz[i] * fe_0 + tk_yyz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyyz_xxyzzz[i] * fz_0;

        tk_xyyz_xxzzzz[i] = 2.0 * tk_yyz_xzzzz[i] * fe_0 + tk_yyz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xyyz_xxzzzz[i] * fz_0;

        tk_xyyz_xyyyyy[i] = tk_xyy_xyyyyy[i] * pa_z[i] + 2.0 * ts_xyyz_xyyyyy[i] * fz_0;

        tk_xyyz_xyyyyz[i] = tk_yyz_yyyyz[i] * fe_0 + tk_yyz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyyz_xyyyyz[i] * fz_0;

        tk_xyyz_xyyyzz[i] = tk_yyz_yyyzz[i] * fe_0 + tk_yyz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyyz_xyyyzz[i] * fz_0;

        tk_xyyz_xyyzzz[i] = tk_yyz_yyzzz[i] * fe_0 + tk_yyz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyyz_xyyzzz[i] * fz_0;

        tk_xyyz_xyzzzz[i] = tk_yyz_yzzzz[i] * fe_0 + tk_yyz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyyz_xyzzzz[i] * fz_0;

        tk_xyyz_xzzzzz[i] = tk_yyz_zzzzz[i] * fe_0 + tk_yyz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xyyz_xzzzzz[i] * fz_0;

        tk_xyyz_yyyyyy[i] = tk_yyz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyyz_yyyyyy[i] * fz_0;

        tk_xyyz_yyyyyz[i] = tk_yyz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyyz_yyyyyz[i] * fz_0;

        tk_xyyz_yyyyzz[i] = tk_yyz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyyz_yyyyzz[i] * fz_0;

        tk_xyyz_yyyzzz[i] = tk_yyz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyyz_yyyzzz[i] * fz_0;

        tk_xyyz_yyzzzz[i] = tk_yyz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyyz_yyzzzz[i] * fz_0;

        tk_xyyz_yzzzzz[i] = tk_yyz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyyz_yzzzzz[i] * fz_0;

        tk_xyyz_zzzzzz[i] = tk_yyz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyyz_zzzzzz[i] * fz_0;
    }

    // Set up 224-252 components of targeted buffer : GI

    auto tk_xyzz_xxxxxx = pbuffer.data(idx_kin_gi + 224);

    auto tk_xyzz_xxxxxy = pbuffer.data(idx_kin_gi + 225);

    auto tk_xyzz_xxxxxz = pbuffer.data(idx_kin_gi + 226);

    auto tk_xyzz_xxxxyy = pbuffer.data(idx_kin_gi + 227);

    auto tk_xyzz_xxxxyz = pbuffer.data(idx_kin_gi + 228);

    auto tk_xyzz_xxxxzz = pbuffer.data(idx_kin_gi + 229);

    auto tk_xyzz_xxxyyy = pbuffer.data(idx_kin_gi + 230);

    auto tk_xyzz_xxxyyz = pbuffer.data(idx_kin_gi + 231);

    auto tk_xyzz_xxxyzz = pbuffer.data(idx_kin_gi + 232);

    auto tk_xyzz_xxxzzz = pbuffer.data(idx_kin_gi + 233);

    auto tk_xyzz_xxyyyy = pbuffer.data(idx_kin_gi + 234);

    auto tk_xyzz_xxyyyz = pbuffer.data(idx_kin_gi + 235);

    auto tk_xyzz_xxyyzz = pbuffer.data(idx_kin_gi + 236);

    auto tk_xyzz_xxyzzz = pbuffer.data(idx_kin_gi + 237);

    auto tk_xyzz_xxzzzz = pbuffer.data(idx_kin_gi + 238);

    auto tk_xyzz_xyyyyy = pbuffer.data(idx_kin_gi + 239);

    auto tk_xyzz_xyyyyz = pbuffer.data(idx_kin_gi + 240);

    auto tk_xyzz_xyyyzz = pbuffer.data(idx_kin_gi + 241);

    auto tk_xyzz_xyyzzz = pbuffer.data(idx_kin_gi + 242);

    auto tk_xyzz_xyzzzz = pbuffer.data(idx_kin_gi + 243);

    auto tk_xyzz_xzzzzz = pbuffer.data(idx_kin_gi + 244);

    auto tk_xyzz_yyyyyy = pbuffer.data(idx_kin_gi + 245);

    auto tk_xyzz_yyyyyz = pbuffer.data(idx_kin_gi + 246);

    auto tk_xyzz_yyyyzz = pbuffer.data(idx_kin_gi + 247);

    auto tk_xyzz_yyyzzz = pbuffer.data(idx_kin_gi + 248);

    auto tk_xyzz_yyzzzz = pbuffer.data(idx_kin_gi + 249);

    auto tk_xyzz_yzzzzz = pbuffer.data(idx_kin_gi + 250);

    auto tk_xyzz_zzzzzz = pbuffer.data(idx_kin_gi + 251);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tk_xyzz_xxxxxx, \
                             tk_xyzz_xxxxxy, \
                             tk_xyzz_xxxxxz, \
                             tk_xyzz_xxxxyy, \
                             tk_xyzz_xxxxyz, \
                             tk_xyzz_xxxxzz, \
                             tk_xyzz_xxxyyy, \
                             tk_xyzz_xxxyyz, \
                             tk_xyzz_xxxyzz, \
                             tk_xyzz_xxxzzz, \
                             tk_xyzz_xxyyyy, \
                             tk_xyzz_xxyyyz, \
                             tk_xyzz_xxyyzz, \
                             tk_xyzz_xxyzzz, \
                             tk_xyzz_xxzzzz, \
                             tk_xyzz_xyyyyy, \
                             tk_xyzz_xyyyyz, \
                             tk_xyzz_xyyyzz, \
                             tk_xyzz_xyyzzz, \
                             tk_xyzz_xyzzzz, \
                             tk_xyzz_xzzzzz, \
                             tk_xyzz_yyyyyy, \
                             tk_xyzz_yyyyyz, \
                             tk_xyzz_yyyyzz, \
                             tk_xyzz_yyyzzz, \
                             tk_xyzz_yyzzzz, \
                             tk_xyzz_yzzzzz, \
                             tk_xyzz_zzzzzz, \
                             tk_xzz_xxxxxx,  \
                             tk_xzz_xxxxxz,  \
                             tk_xzz_xxxxzz,  \
                             tk_xzz_xxxzzz,  \
                             tk_xzz_xxzzzz,  \
                             tk_xzz_xzzzzz,  \
                             tk_yzz_xxxxxy,  \
                             tk_yzz_xxxxy,   \
                             tk_yzz_xxxxyy,  \
                             tk_yzz_xxxxyz,  \
                             tk_yzz_xxxyy,   \
                             tk_yzz_xxxyyy,  \
                             tk_yzz_xxxyyz,  \
                             tk_yzz_xxxyz,   \
                             tk_yzz_xxxyzz,  \
                             tk_yzz_xxyyy,   \
                             tk_yzz_xxyyyy,  \
                             tk_yzz_xxyyyz,  \
                             tk_yzz_xxyyz,   \
                             tk_yzz_xxyyzz,  \
                             tk_yzz_xxyzz,   \
                             tk_yzz_xxyzzz,  \
                             tk_yzz_xyyyy,   \
                             tk_yzz_xyyyyy,  \
                             tk_yzz_xyyyyz,  \
                             tk_yzz_xyyyz,   \
                             tk_yzz_xyyyzz,  \
                             tk_yzz_xyyzz,   \
                             tk_yzz_xyyzzz,  \
                             tk_yzz_xyzzz,   \
                             tk_yzz_xyzzzz,  \
                             tk_yzz_yyyyy,   \
                             tk_yzz_yyyyyy,  \
                             tk_yzz_yyyyyz,  \
                             tk_yzz_yyyyz,   \
                             tk_yzz_yyyyzz,  \
                             tk_yzz_yyyzz,   \
                             tk_yzz_yyyzzz,  \
                             tk_yzz_yyzzz,   \
                             tk_yzz_yyzzzz,  \
                             tk_yzz_yzzzz,   \
                             tk_yzz_yzzzzz,  \
                             tk_yzz_zzzzzz,  \
                             ts_xyzz_xxxxxx, \
                             ts_xyzz_xxxxxy, \
                             ts_xyzz_xxxxxz, \
                             ts_xyzz_xxxxyy, \
                             ts_xyzz_xxxxyz, \
                             ts_xyzz_xxxxzz, \
                             ts_xyzz_xxxyyy, \
                             ts_xyzz_xxxyyz, \
                             ts_xyzz_xxxyzz, \
                             ts_xyzz_xxxzzz, \
                             ts_xyzz_xxyyyy, \
                             ts_xyzz_xxyyyz, \
                             ts_xyzz_xxyyzz, \
                             ts_xyzz_xxyzzz, \
                             ts_xyzz_xxzzzz, \
                             ts_xyzz_xyyyyy, \
                             ts_xyzz_xyyyyz, \
                             ts_xyzz_xyyyzz, \
                             ts_xyzz_xyyzzz, \
                             ts_xyzz_xyzzzz, \
                             ts_xyzz_xzzzzz, \
                             ts_xyzz_yyyyyy, \
                             ts_xyzz_yyyyyz, \
                             ts_xyzz_yyyyzz, \
                             ts_xyzz_yyyzzz, \
                             ts_xyzz_yyzzzz, \
                             ts_xyzz_yzzzzz, \
                             ts_xyzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzz_xxxxxx[i] = tk_xzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_xyzz_xxxxxx[i] * fz_0;

        tk_xyzz_xxxxxy[i] = 5.0 * tk_yzz_xxxxy[i] * fe_0 + tk_yzz_xxxxxy[i] * pa_x[i] + 2.0 * ts_xyzz_xxxxxy[i] * fz_0;

        tk_xyzz_xxxxxz[i] = tk_xzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_xyzz_xxxxxz[i] * fz_0;

        tk_xyzz_xxxxyy[i] = 4.0 * tk_yzz_xxxyy[i] * fe_0 + tk_yzz_xxxxyy[i] * pa_x[i] + 2.0 * ts_xyzz_xxxxyy[i] * fz_0;

        tk_xyzz_xxxxyz[i] = 4.0 * tk_yzz_xxxyz[i] * fe_0 + tk_yzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyzz_xxxxyz[i] * fz_0;

        tk_xyzz_xxxxzz[i] = tk_xzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_xyzz_xxxxzz[i] * fz_0;

        tk_xyzz_xxxyyy[i] = 3.0 * tk_yzz_xxyyy[i] * fe_0 + tk_yzz_xxxyyy[i] * pa_x[i] + 2.0 * ts_xyzz_xxxyyy[i] * fz_0;

        tk_xyzz_xxxyyz[i] = 3.0 * tk_yzz_xxyyz[i] * fe_0 + tk_yzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyzz_xxxyyz[i] * fz_0;

        tk_xyzz_xxxyzz[i] = 3.0 * tk_yzz_xxyzz[i] * fe_0 + tk_yzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyzz_xxxyzz[i] * fz_0;

        tk_xyzz_xxxzzz[i] = tk_xzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_xyzz_xxxzzz[i] * fz_0;

        tk_xyzz_xxyyyy[i] = 2.0 * tk_yzz_xyyyy[i] * fe_0 + tk_yzz_xxyyyy[i] * pa_x[i] + 2.0 * ts_xyzz_xxyyyy[i] * fz_0;

        tk_xyzz_xxyyyz[i] = 2.0 * tk_yzz_xyyyz[i] * fe_0 + tk_yzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyzz_xxyyyz[i] * fz_0;

        tk_xyzz_xxyyzz[i] = 2.0 * tk_yzz_xyyzz[i] * fe_0 + tk_yzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyzz_xxyyzz[i] * fz_0;

        tk_xyzz_xxyzzz[i] = 2.0 * tk_yzz_xyzzz[i] * fe_0 + tk_yzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyzz_xxyzzz[i] * fz_0;

        tk_xyzz_xxzzzz[i] = tk_xzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_xyzz_xxzzzz[i] * fz_0;

        tk_xyzz_xyyyyy[i] = tk_yzz_yyyyy[i] * fe_0 + tk_yzz_xyyyyy[i] * pa_x[i] + 2.0 * ts_xyzz_xyyyyy[i] * fz_0;

        tk_xyzz_xyyyyz[i] = tk_yzz_yyyyz[i] * fe_0 + tk_yzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyzz_xyyyyz[i] * fz_0;

        tk_xyzz_xyyyzz[i] = tk_yzz_yyyzz[i] * fe_0 + tk_yzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyzz_xyyyzz[i] * fz_0;

        tk_xyzz_xyyzzz[i] = tk_yzz_yyzzz[i] * fe_0 + tk_yzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyzz_xyyzzz[i] * fz_0;

        tk_xyzz_xyzzzz[i] = tk_yzz_yzzzz[i] * fe_0 + tk_yzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyzz_xyzzzz[i] * fz_0;

        tk_xyzz_xzzzzz[i] = tk_xzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_xyzz_xzzzzz[i] * fz_0;

        tk_xyzz_yyyyyy[i] = tk_yzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyzz_yyyyyy[i] * fz_0;

        tk_xyzz_yyyyyz[i] = tk_yzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyzz_yyyyyz[i] * fz_0;

        tk_xyzz_yyyyzz[i] = tk_yzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyzz_yyyyzz[i] * fz_0;

        tk_xyzz_yyyzzz[i] = tk_yzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyzz_yyyzzz[i] * fz_0;

        tk_xyzz_yyzzzz[i] = tk_yzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyzz_yyzzzz[i] * fz_0;

        tk_xyzz_yzzzzz[i] = tk_yzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyzz_yzzzzz[i] * fz_0;

        tk_xyzz_zzzzzz[i] = tk_yzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyzz_zzzzzz[i] * fz_0;
    }

    // Set up 252-280 components of targeted buffer : GI

    auto tk_xzzz_xxxxxx = pbuffer.data(idx_kin_gi + 252);

    auto tk_xzzz_xxxxxy = pbuffer.data(idx_kin_gi + 253);

    auto tk_xzzz_xxxxxz = pbuffer.data(idx_kin_gi + 254);

    auto tk_xzzz_xxxxyy = pbuffer.data(idx_kin_gi + 255);

    auto tk_xzzz_xxxxyz = pbuffer.data(idx_kin_gi + 256);

    auto tk_xzzz_xxxxzz = pbuffer.data(idx_kin_gi + 257);

    auto tk_xzzz_xxxyyy = pbuffer.data(idx_kin_gi + 258);

    auto tk_xzzz_xxxyyz = pbuffer.data(idx_kin_gi + 259);

    auto tk_xzzz_xxxyzz = pbuffer.data(idx_kin_gi + 260);

    auto tk_xzzz_xxxzzz = pbuffer.data(idx_kin_gi + 261);

    auto tk_xzzz_xxyyyy = pbuffer.data(idx_kin_gi + 262);

    auto tk_xzzz_xxyyyz = pbuffer.data(idx_kin_gi + 263);

    auto tk_xzzz_xxyyzz = pbuffer.data(idx_kin_gi + 264);

    auto tk_xzzz_xxyzzz = pbuffer.data(idx_kin_gi + 265);

    auto tk_xzzz_xxzzzz = pbuffer.data(idx_kin_gi + 266);

    auto tk_xzzz_xyyyyy = pbuffer.data(idx_kin_gi + 267);

    auto tk_xzzz_xyyyyz = pbuffer.data(idx_kin_gi + 268);

    auto tk_xzzz_xyyyzz = pbuffer.data(idx_kin_gi + 269);

    auto tk_xzzz_xyyzzz = pbuffer.data(idx_kin_gi + 270);

    auto tk_xzzz_xyzzzz = pbuffer.data(idx_kin_gi + 271);

    auto tk_xzzz_xzzzzz = pbuffer.data(idx_kin_gi + 272);

    auto tk_xzzz_yyyyyy = pbuffer.data(idx_kin_gi + 273);

    auto tk_xzzz_yyyyyz = pbuffer.data(idx_kin_gi + 274);

    auto tk_xzzz_yyyyzz = pbuffer.data(idx_kin_gi + 275);

    auto tk_xzzz_yyyzzz = pbuffer.data(idx_kin_gi + 276);

    auto tk_xzzz_yyzzzz = pbuffer.data(idx_kin_gi + 277);

    auto tk_xzzz_yzzzzz = pbuffer.data(idx_kin_gi + 278);

    auto tk_xzzz_zzzzzz = pbuffer.data(idx_kin_gi + 279);

#pragma omp simd aligned(pa_x,               \
                             tk_xzzz_xxxxxx, \
                             tk_xzzz_xxxxxy, \
                             tk_xzzz_xxxxxz, \
                             tk_xzzz_xxxxyy, \
                             tk_xzzz_xxxxyz, \
                             tk_xzzz_xxxxzz, \
                             tk_xzzz_xxxyyy, \
                             tk_xzzz_xxxyyz, \
                             tk_xzzz_xxxyzz, \
                             tk_xzzz_xxxzzz, \
                             tk_xzzz_xxyyyy, \
                             tk_xzzz_xxyyyz, \
                             tk_xzzz_xxyyzz, \
                             tk_xzzz_xxyzzz, \
                             tk_xzzz_xxzzzz, \
                             tk_xzzz_xyyyyy, \
                             tk_xzzz_xyyyyz, \
                             tk_xzzz_xyyyzz, \
                             tk_xzzz_xyyzzz, \
                             tk_xzzz_xyzzzz, \
                             tk_xzzz_xzzzzz, \
                             tk_xzzz_yyyyyy, \
                             tk_xzzz_yyyyyz, \
                             tk_xzzz_yyyyzz, \
                             tk_xzzz_yyyzzz, \
                             tk_xzzz_yyzzzz, \
                             tk_xzzz_yzzzzz, \
                             tk_xzzz_zzzzzz, \
                             tk_zzz_xxxxx,   \
                             tk_zzz_xxxxxx,  \
                             tk_zzz_xxxxxy,  \
                             tk_zzz_xxxxxz,  \
                             tk_zzz_xxxxy,   \
                             tk_zzz_xxxxyy,  \
                             tk_zzz_xxxxyz,  \
                             tk_zzz_xxxxz,   \
                             tk_zzz_xxxxzz,  \
                             tk_zzz_xxxyy,   \
                             tk_zzz_xxxyyy,  \
                             tk_zzz_xxxyyz,  \
                             tk_zzz_xxxyz,   \
                             tk_zzz_xxxyzz,  \
                             tk_zzz_xxxzz,   \
                             tk_zzz_xxxzzz,  \
                             tk_zzz_xxyyy,   \
                             tk_zzz_xxyyyy,  \
                             tk_zzz_xxyyyz,  \
                             tk_zzz_xxyyz,   \
                             tk_zzz_xxyyzz,  \
                             tk_zzz_xxyzz,   \
                             tk_zzz_xxyzzz,  \
                             tk_zzz_xxzzz,   \
                             tk_zzz_xxzzzz,  \
                             tk_zzz_xyyyy,   \
                             tk_zzz_xyyyyy,  \
                             tk_zzz_xyyyyz,  \
                             tk_zzz_xyyyz,   \
                             tk_zzz_xyyyzz,  \
                             tk_zzz_xyyzz,   \
                             tk_zzz_xyyzzz,  \
                             tk_zzz_xyzzz,   \
                             tk_zzz_xyzzzz,  \
                             tk_zzz_xzzzz,   \
                             tk_zzz_xzzzzz,  \
                             tk_zzz_yyyyy,   \
                             tk_zzz_yyyyyy,  \
                             tk_zzz_yyyyyz,  \
                             tk_zzz_yyyyz,   \
                             tk_zzz_yyyyzz,  \
                             tk_zzz_yyyzz,   \
                             tk_zzz_yyyzzz,  \
                             tk_zzz_yyzzz,   \
                             tk_zzz_yyzzzz,  \
                             tk_zzz_yzzzz,   \
                             tk_zzz_yzzzzz,  \
                             tk_zzz_zzzzz,   \
                             tk_zzz_zzzzzz,  \
                             ts_xzzz_xxxxxx, \
                             ts_xzzz_xxxxxy, \
                             ts_xzzz_xxxxxz, \
                             ts_xzzz_xxxxyy, \
                             ts_xzzz_xxxxyz, \
                             ts_xzzz_xxxxzz, \
                             ts_xzzz_xxxyyy, \
                             ts_xzzz_xxxyyz, \
                             ts_xzzz_xxxyzz, \
                             ts_xzzz_xxxzzz, \
                             ts_xzzz_xxyyyy, \
                             ts_xzzz_xxyyyz, \
                             ts_xzzz_xxyyzz, \
                             ts_xzzz_xxyzzz, \
                             ts_xzzz_xxzzzz, \
                             ts_xzzz_xyyyyy, \
                             ts_xzzz_xyyyyz, \
                             ts_xzzz_xyyyzz, \
                             ts_xzzz_xyyzzz, \
                             ts_xzzz_xyzzzz, \
                             ts_xzzz_xzzzzz, \
                             ts_xzzz_yyyyyy, \
                             ts_xzzz_yyyyyz, \
                             ts_xzzz_yyyyzz, \
                             ts_xzzz_yyyzzz, \
                             ts_xzzz_yyzzzz, \
                             ts_xzzz_yzzzzz, \
                             ts_xzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzz_xxxxxx[i] = 6.0 * tk_zzz_xxxxx[i] * fe_0 + tk_zzz_xxxxxx[i] * pa_x[i] + 2.0 * ts_xzzz_xxxxxx[i] * fz_0;

        tk_xzzz_xxxxxy[i] = 5.0 * tk_zzz_xxxxy[i] * fe_0 + tk_zzz_xxxxxy[i] * pa_x[i] + 2.0 * ts_xzzz_xxxxxy[i] * fz_0;

        tk_xzzz_xxxxxz[i] = 5.0 * tk_zzz_xxxxz[i] * fe_0 + tk_zzz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xzzz_xxxxxz[i] * fz_0;

        tk_xzzz_xxxxyy[i] = 4.0 * tk_zzz_xxxyy[i] * fe_0 + tk_zzz_xxxxyy[i] * pa_x[i] + 2.0 * ts_xzzz_xxxxyy[i] * fz_0;

        tk_xzzz_xxxxyz[i] = 4.0 * tk_zzz_xxxyz[i] * fe_0 + tk_zzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xzzz_xxxxyz[i] * fz_0;

        tk_xzzz_xxxxzz[i] = 4.0 * tk_zzz_xxxzz[i] * fe_0 + tk_zzz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xzzz_xxxxzz[i] * fz_0;

        tk_xzzz_xxxyyy[i] = 3.0 * tk_zzz_xxyyy[i] * fe_0 + tk_zzz_xxxyyy[i] * pa_x[i] + 2.0 * ts_xzzz_xxxyyy[i] * fz_0;

        tk_xzzz_xxxyyz[i] = 3.0 * tk_zzz_xxyyz[i] * fe_0 + tk_zzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xzzz_xxxyyz[i] * fz_0;

        tk_xzzz_xxxyzz[i] = 3.0 * tk_zzz_xxyzz[i] * fe_0 + tk_zzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xzzz_xxxyzz[i] * fz_0;

        tk_xzzz_xxxzzz[i] = 3.0 * tk_zzz_xxzzz[i] * fe_0 + tk_zzz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xzzz_xxxzzz[i] * fz_0;

        tk_xzzz_xxyyyy[i] = 2.0 * tk_zzz_xyyyy[i] * fe_0 + tk_zzz_xxyyyy[i] * pa_x[i] + 2.0 * ts_xzzz_xxyyyy[i] * fz_0;

        tk_xzzz_xxyyyz[i] = 2.0 * tk_zzz_xyyyz[i] * fe_0 + tk_zzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xzzz_xxyyyz[i] * fz_0;

        tk_xzzz_xxyyzz[i] = 2.0 * tk_zzz_xyyzz[i] * fe_0 + tk_zzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xzzz_xxyyzz[i] * fz_0;

        tk_xzzz_xxyzzz[i] = 2.0 * tk_zzz_xyzzz[i] * fe_0 + tk_zzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xzzz_xxyzzz[i] * fz_0;

        tk_xzzz_xxzzzz[i] = 2.0 * tk_zzz_xzzzz[i] * fe_0 + tk_zzz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xzzz_xxzzzz[i] * fz_0;

        tk_xzzz_xyyyyy[i] = tk_zzz_yyyyy[i] * fe_0 + tk_zzz_xyyyyy[i] * pa_x[i] + 2.0 * ts_xzzz_xyyyyy[i] * fz_0;

        tk_xzzz_xyyyyz[i] = tk_zzz_yyyyz[i] * fe_0 + tk_zzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xzzz_xyyyyz[i] * fz_0;

        tk_xzzz_xyyyzz[i] = tk_zzz_yyyzz[i] * fe_0 + tk_zzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xzzz_xyyyzz[i] * fz_0;

        tk_xzzz_xyyzzz[i] = tk_zzz_yyzzz[i] * fe_0 + tk_zzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xzzz_xyyzzz[i] * fz_0;

        tk_xzzz_xyzzzz[i] = tk_zzz_yzzzz[i] * fe_0 + tk_zzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xzzz_xyzzzz[i] * fz_0;

        tk_xzzz_xzzzzz[i] = tk_zzz_zzzzz[i] * fe_0 + tk_zzz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xzzz_xzzzzz[i] * fz_0;

        tk_xzzz_yyyyyy[i] = tk_zzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xzzz_yyyyyy[i] * fz_0;

        tk_xzzz_yyyyyz[i] = tk_zzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xzzz_yyyyyz[i] * fz_0;

        tk_xzzz_yyyyzz[i] = tk_zzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xzzz_yyyyzz[i] * fz_0;

        tk_xzzz_yyyzzz[i] = tk_zzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xzzz_yyyzzz[i] * fz_0;

        tk_xzzz_yyzzzz[i] = tk_zzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xzzz_yyzzzz[i] * fz_0;

        tk_xzzz_yzzzzz[i] = tk_zzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xzzz_yzzzzz[i] * fz_0;

        tk_xzzz_zzzzzz[i] = tk_zzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xzzz_zzzzzz[i] * fz_0;
    }

    // Set up 280-308 components of targeted buffer : GI

    auto tk_yyyy_xxxxxx = pbuffer.data(idx_kin_gi + 280);

    auto tk_yyyy_xxxxxy = pbuffer.data(idx_kin_gi + 281);

    auto tk_yyyy_xxxxxz = pbuffer.data(idx_kin_gi + 282);

    auto tk_yyyy_xxxxyy = pbuffer.data(idx_kin_gi + 283);

    auto tk_yyyy_xxxxyz = pbuffer.data(idx_kin_gi + 284);

    auto tk_yyyy_xxxxzz = pbuffer.data(idx_kin_gi + 285);

    auto tk_yyyy_xxxyyy = pbuffer.data(idx_kin_gi + 286);

    auto tk_yyyy_xxxyyz = pbuffer.data(idx_kin_gi + 287);

    auto tk_yyyy_xxxyzz = pbuffer.data(idx_kin_gi + 288);

    auto tk_yyyy_xxxzzz = pbuffer.data(idx_kin_gi + 289);

    auto tk_yyyy_xxyyyy = pbuffer.data(idx_kin_gi + 290);

    auto tk_yyyy_xxyyyz = pbuffer.data(idx_kin_gi + 291);

    auto tk_yyyy_xxyyzz = pbuffer.data(idx_kin_gi + 292);

    auto tk_yyyy_xxyzzz = pbuffer.data(idx_kin_gi + 293);

    auto tk_yyyy_xxzzzz = pbuffer.data(idx_kin_gi + 294);

    auto tk_yyyy_xyyyyy = pbuffer.data(idx_kin_gi + 295);

    auto tk_yyyy_xyyyyz = pbuffer.data(idx_kin_gi + 296);

    auto tk_yyyy_xyyyzz = pbuffer.data(idx_kin_gi + 297);

    auto tk_yyyy_xyyzzz = pbuffer.data(idx_kin_gi + 298);

    auto tk_yyyy_xyzzzz = pbuffer.data(idx_kin_gi + 299);

    auto tk_yyyy_xzzzzz = pbuffer.data(idx_kin_gi + 300);

    auto tk_yyyy_yyyyyy = pbuffer.data(idx_kin_gi + 301);

    auto tk_yyyy_yyyyyz = pbuffer.data(idx_kin_gi + 302);

    auto tk_yyyy_yyyyzz = pbuffer.data(idx_kin_gi + 303);

    auto tk_yyyy_yyyzzz = pbuffer.data(idx_kin_gi + 304);

    auto tk_yyyy_yyzzzz = pbuffer.data(idx_kin_gi + 305);

    auto tk_yyyy_yzzzzz = pbuffer.data(idx_kin_gi + 306);

    auto tk_yyyy_zzzzzz = pbuffer.data(idx_kin_gi + 307);

#pragma omp simd aligned(pa_y,               \
                             tk_yy_xxxxxx,   \
                             tk_yy_xxxxxy,   \
                             tk_yy_xxxxxz,   \
                             tk_yy_xxxxyy,   \
                             tk_yy_xxxxyz,   \
                             tk_yy_xxxxzz,   \
                             tk_yy_xxxyyy,   \
                             tk_yy_xxxyyz,   \
                             tk_yy_xxxyzz,   \
                             tk_yy_xxxzzz,   \
                             tk_yy_xxyyyy,   \
                             tk_yy_xxyyyz,   \
                             tk_yy_xxyyzz,   \
                             tk_yy_xxyzzz,   \
                             tk_yy_xxzzzz,   \
                             tk_yy_xyyyyy,   \
                             tk_yy_xyyyyz,   \
                             tk_yy_xyyyzz,   \
                             tk_yy_xyyzzz,   \
                             tk_yy_xyzzzz,   \
                             tk_yy_xzzzzz,   \
                             tk_yy_yyyyyy,   \
                             tk_yy_yyyyyz,   \
                             tk_yy_yyyyzz,   \
                             tk_yy_yyyzzz,   \
                             tk_yy_yyzzzz,   \
                             tk_yy_yzzzzz,   \
                             tk_yy_zzzzzz,   \
                             tk_yyy_xxxxx,   \
                             tk_yyy_xxxxxx,  \
                             tk_yyy_xxxxxy,  \
                             tk_yyy_xxxxxz,  \
                             tk_yyy_xxxxy,   \
                             tk_yyy_xxxxyy,  \
                             tk_yyy_xxxxyz,  \
                             tk_yyy_xxxxz,   \
                             tk_yyy_xxxxzz,  \
                             tk_yyy_xxxyy,   \
                             tk_yyy_xxxyyy,  \
                             tk_yyy_xxxyyz,  \
                             tk_yyy_xxxyz,   \
                             tk_yyy_xxxyzz,  \
                             tk_yyy_xxxzz,   \
                             tk_yyy_xxxzzz,  \
                             tk_yyy_xxyyy,   \
                             tk_yyy_xxyyyy,  \
                             tk_yyy_xxyyyz,  \
                             tk_yyy_xxyyz,   \
                             tk_yyy_xxyyzz,  \
                             tk_yyy_xxyzz,   \
                             tk_yyy_xxyzzz,  \
                             tk_yyy_xxzzz,   \
                             tk_yyy_xxzzzz,  \
                             tk_yyy_xyyyy,   \
                             tk_yyy_xyyyyy,  \
                             tk_yyy_xyyyyz,  \
                             tk_yyy_xyyyz,   \
                             tk_yyy_xyyyzz,  \
                             tk_yyy_xyyzz,   \
                             tk_yyy_xyyzzz,  \
                             tk_yyy_xyzzz,   \
                             tk_yyy_xyzzzz,  \
                             tk_yyy_xzzzz,   \
                             tk_yyy_xzzzzz,  \
                             tk_yyy_yyyyy,   \
                             tk_yyy_yyyyyy,  \
                             tk_yyy_yyyyyz,  \
                             tk_yyy_yyyyz,   \
                             tk_yyy_yyyyzz,  \
                             tk_yyy_yyyzz,   \
                             tk_yyy_yyyzzz,  \
                             tk_yyy_yyzzz,   \
                             tk_yyy_yyzzzz,  \
                             tk_yyy_yzzzz,   \
                             tk_yyy_yzzzzz,  \
                             tk_yyy_zzzzz,   \
                             tk_yyy_zzzzzz,  \
                             tk_yyyy_xxxxxx, \
                             tk_yyyy_xxxxxy, \
                             tk_yyyy_xxxxxz, \
                             tk_yyyy_xxxxyy, \
                             tk_yyyy_xxxxyz, \
                             tk_yyyy_xxxxzz, \
                             tk_yyyy_xxxyyy, \
                             tk_yyyy_xxxyyz, \
                             tk_yyyy_xxxyzz, \
                             tk_yyyy_xxxzzz, \
                             tk_yyyy_xxyyyy, \
                             tk_yyyy_xxyyyz, \
                             tk_yyyy_xxyyzz, \
                             tk_yyyy_xxyzzz, \
                             tk_yyyy_xxzzzz, \
                             tk_yyyy_xyyyyy, \
                             tk_yyyy_xyyyyz, \
                             tk_yyyy_xyyyzz, \
                             tk_yyyy_xyyzzz, \
                             tk_yyyy_xyzzzz, \
                             tk_yyyy_xzzzzz, \
                             tk_yyyy_yyyyyy, \
                             tk_yyyy_yyyyyz, \
                             tk_yyyy_yyyyzz, \
                             tk_yyyy_yyyzzz, \
                             tk_yyyy_yyzzzz, \
                             tk_yyyy_yzzzzz, \
                             tk_yyyy_zzzzzz, \
                             ts_yy_xxxxxx,   \
                             ts_yy_xxxxxy,   \
                             ts_yy_xxxxxz,   \
                             ts_yy_xxxxyy,   \
                             ts_yy_xxxxyz,   \
                             ts_yy_xxxxzz,   \
                             ts_yy_xxxyyy,   \
                             ts_yy_xxxyyz,   \
                             ts_yy_xxxyzz,   \
                             ts_yy_xxxzzz,   \
                             ts_yy_xxyyyy,   \
                             ts_yy_xxyyyz,   \
                             ts_yy_xxyyzz,   \
                             ts_yy_xxyzzz,   \
                             ts_yy_xxzzzz,   \
                             ts_yy_xyyyyy,   \
                             ts_yy_xyyyyz,   \
                             ts_yy_xyyyzz,   \
                             ts_yy_xyyzzz,   \
                             ts_yy_xyzzzz,   \
                             ts_yy_xzzzzz,   \
                             ts_yy_yyyyyy,   \
                             ts_yy_yyyyyz,   \
                             ts_yy_yyyyzz,   \
                             ts_yy_yyyzzz,   \
                             ts_yy_yyzzzz,   \
                             ts_yy_yzzzzz,   \
                             ts_yy_zzzzzz,   \
                             ts_yyyy_xxxxxx, \
                             ts_yyyy_xxxxxy, \
                             ts_yyyy_xxxxxz, \
                             ts_yyyy_xxxxyy, \
                             ts_yyyy_xxxxyz, \
                             ts_yyyy_xxxxzz, \
                             ts_yyyy_xxxyyy, \
                             ts_yyyy_xxxyyz, \
                             ts_yyyy_xxxyzz, \
                             ts_yyyy_xxxzzz, \
                             ts_yyyy_xxyyyy, \
                             ts_yyyy_xxyyyz, \
                             ts_yyyy_xxyyzz, \
                             ts_yyyy_xxyzzz, \
                             ts_yyyy_xxzzzz, \
                             ts_yyyy_xyyyyy, \
                             ts_yyyy_xyyyyz, \
                             ts_yyyy_xyyyzz, \
                             ts_yyyy_xyyzzz, \
                             ts_yyyy_xyzzzz, \
                             ts_yyyy_xzzzzz, \
                             ts_yyyy_yyyyyy, \
                             ts_yyyy_yyyyyz, \
                             ts_yyyy_yyyyzz, \
                             ts_yyyy_yyyzzz, \
                             ts_yyyy_yyzzzz, \
                             ts_yyyy_yzzzzz, \
                             ts_yyyy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyy_xxxxxx[i] =
            -6.0 * ts_yy_xxxxxx[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxxxx[i] * fe_0 + tk_yyy_xxxxxx[i] * pa_y[i] + 2.0 * ts_yyyy_xxxxxx[i] * fz_0;

        tk_yyyy_xxxxxy[i] = -6.0 * ts_yy_xxxxxy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxxxy[i] * fe_0 + tk_yyy_xxxxx[i] * fe_0 +
                            tk_yyy_xxxxxy[i] * pa_y[i] + 2.0 * ts_yyyy_xxxxxy[i] * fz_0;

        tk_yyyy_xxxxxz[i] =
            -6.0 * ts_yy_xxxxxz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxxxz[i] * fe_0 + tk_yyy_xxxxxz[i] * pa_y[i] + 2.0 * ts_yyyy_xxxxxz[i] * fz_0;

        tk_yyyy_xxxxyy[i] = -6.0 * ts_yy_xxxxyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxxyy[i] * fe_0 + 2.0 * tk_yyy_xxxxy[i] * fe_0 +
                            tk_yyy_xxxxyy[i] * pa_y[i] + 2.0 * ts_yyyy_xxxxyy[i] * fz_0;

        tk_yyyy_xxxxyz[i] = -6.0 * ts_yy_xxxxyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxxyz[i] * fe_0 + tk_yyy_xxxxz[i] * fe_0 +
                            tk_yyy_xxxxyz[i] * pa_y[i] + 2.0 * ts_yyyy_xxxxyz[i] * fz_0;

        tk_yyyy_xxxxzz[i] =
            -6.0 * ts_yy_xxxxzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxxzz[i] * fe_0 + tk_yyy_xxxxzz[i] * pa_y[i] + 2.0 * ts_yyyy_xxxxzz[i] * fz_0;

        tk_yyyy_xxxyyy[i] = -6.0 * ts_yy_xxxyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxyyy[i] * fe_0 + 3.0 * tk_yyy_xxxyy[i] * fe_0 +
                            tk_yyy_xxxyyy[i] * pa_y[i] + 2.0 * ts_yyyy_xxxyyy[i] * fz_0;

        tk_yyyy_xxxyyz[i] = -6.0 * ts_yy_xxxyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxyyz[i] * fe_0 + 2.0 * tk_yyy_xxxyz[i] * fe_0 +
                            tk_yyy_xxxyyz[i] * pa_y[i] + 2.0 * ts_yyyy_xxxyyz[i] * fz_0;

        tk_yyyy_xxxyzz[i] = -6.0 * ts_yy_xxxyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxyzz[i] * fe_0 + tk_yyy_xxxzz[i] * fe_0 +
                            tk_yyy_xxxyzz[i] * pa_y[i] + 2.0 * ts_yyyy_xxxyzz[i] * fz_0;

        tk_yyyy_xxxzzz[i] =
            -6.0 * ts_yy_xxxzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxxzzz[i] * fe_0 + tk_yyy_xxxzzz[i] * pa_y[i] + 2.0 * ts_yyyy_xxxzzz[i] * fz_0;

        tk_yyyy_xxyyyy[i] = -6.0 * ts_yy_xxyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxyyyy[i] * fe_0 + 4.0 * tk_yyy_xxyyy[i] * fe_0 +
                            tk_yyy_xxyyyy[i] * pa_y[i] + 2.0 * ts_yyyy_xxyyyy[i] * fz_0;

        tk_yyyy_xxyyyz[i] = -6.0 * ts_yy_xxyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxyyyz[i] * fe_0 + 3.0 * tk_yyy_xxyyz[i] * fe_0 +
                            tk_yyy_xxyyyz[i] * pa_y[i] + 2.0 * ts_yyyy_xxyyyz[i] * fz_0;

        tk_yyyy_xxyyzz[i] = -6.0 * ts_yy_xxyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxyyzz[i] * fe_0 + 2.0 * tk_yyy_xxyzz[i] * fe_0 +
                            tk_yyy_xxyyzz[i] * pa_y[i] + 2.0 * ts_yyyy_xxyyzz[i] * fz_0;

        tk_yyyy_xxyzzz[i] = -6.0 * ts_yy_xxyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxyzzz[i] * fe_0 + tk_yyy_xxzzz[i] * fe_0 +
                            tk_yyy_xxyzzz[i] * pa_y[i] + 2.0 * ts_yyyy_xxyzzz[i] * fz_0;

        tk_yyyy_xxzzzz[i] =
            -6.0 * ts_yy_xxzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xxzzzz[i] * fe_0 + tk_yyy_xxzzzz[i] * pa_y[i] + 2.0 * ts_yyyy_xxzzzz[i] * fz_0;

        tk_yyyy_xyyyyy[i] = -6.0 * ts_yy_xyyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyyyyy[i] * fe_0 + 5.0 * tk_yyy_xyyyy[i] * fe_0 +
                            tk_yyy_xyyyyy[i] * pa_y[i] + 2.0 * ts_yyyy_xyyyyy[i] * fz_0;

        tk_yyyy_xyyyyz[i] = -6.0 * ts_yy_xyyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyyyyz[i] * fe_0 + 4.0 * tk_yyy_xyyyz[i] * fe_0 +
                            tk_yyy_xyyyyz[i] * pa_y[i] + 2.0 * ts_yyyy_xyyyyz[i] * fz_0;

        tk_yyyy_xyyyzz[i] = -6.0 * ts_yy_xyyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyyyzz[i] * fe_0 + 3.0 * tk_yyy_xyyzz[i] * fe_0 +
                            tk_yyy_xyyyzz[i] * pa_y[i] + 2.0 * ts_yyyy_xyyyzz[i] * fz_0;

        tk_yyyy_xyyzzz[i] = -6.0 * ts_yy_xyyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyyzzz[i] * fe_0 + 2.0 * tk_yyy_xyzzz[i] * fe_0 +
                            tk_yyy_xyyzzz[i] * pa_y[i] + 2.0 * ts_yyyy_xyyzzz[i] * fz_0;

        tk_yyyy_xyzzzz[i] = -6.0 * ts_yy_xyzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xyzzzz[i] * fe_0 + tk_yyy_xzzzz[i] * fe_0 +
                            tk_yyy_xyzzzz[i] * pa_y[i] + 2.0 * ts_yyyy_xyzzzz[i] * fz_0;

        tk_yyyy_xzzzzz[i] =
            -6.0 * ts_yy_xzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xzzzzz[i] * fe_0 + tk_yyy_xzzzzz[i] * pa_y[i] + 2.0 * ts_yyyy_xzzzzz[i] * fz_0;

        tk_yyyy_yyyyyy[i] = -6.0 * ts_yy_yyyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyyyyy[i] * fe_0 + 6.0 * tk_yyy_yyyyy[i] * fe_0 +
                            tk_yyy_yyyyyy[i] * pa_y[i] + 2.0 * ts_yyyy_yyyyyy[i] * fz_0;

        tk_yyyy_yyyyyz[i] = -6.0 * ts_yy_yyyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyyyyz[i] * fe_0 + 5.0 * tk_yyy_yyyyz[i] * fe_0 +
                            tk_yyy_yyyyyz[i] * pa_y[i] + 2.0 * ts_yyyy_yyyyyz[i] * fz_0;

        tk_yyyy_yyyyzz[i] = -6.0 * ts_yy_yyyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyyyzz[i] * fe_0 + 4.0 * tk_yyy_yyyzz[i] * fe_0 +
                            tk_yyy_yyyyzz[i] * pa_y[i] + 2.0 * ts_yyyy_yyyyzz[i] * fz_0;

        tk_yyyy_yyyzzz[i] = -6.0 * ts_yy_yyyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyyzzz[i] * fe_0 + 3.0 * tk_yyy_yyzzz[i] * fe_0 +
                            tk_yyy_yyyzzz[i] * pa_y[i] + 2.0 * ts_yyyy_yyyzzz[i] * fz_0;

        tk_yyyy_yyzzzz[i] = -6.0 * ts_yy_yyzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yyzzzz[i] * fe_0 + 2.0 * tk_yyy_yzzzz[i] * fe_0 +
                            tk_yyy_yyzzzz[i] * pa_y[i] + 2.0 * ts_yyyy_yyzzzz[i] * fz_0;

        tk_yyyy_yzzzzz[i] = -6.0 * ts_yy_yzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yzzzzz[i] * fe_0 + tk_yyy_zzzzz[i] * fe_0 +
                            tk_yyy_yzzzzz[i] * pa_y[i] + 2.0 * ts_yyyy_yzzzzz[i] * fz_0;

        tk_yyyy_zzzzzz[i] =
            -6.0 * ts_yy_zzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_zzzzzz[i] * fe_0 + tk_yyy_zzzzzz[i] * pa_y[i] + 2.0 * ts_yyyy_zzzzzz[i] * fz_0;
    }

    // Set up 308-336 components of targeted buffer : GI

    auto tk_yyyz_xxxxxx = pbuffer.data(idx_kin_gi + 308);

    auto tk_yyyz_xxxxxy = pbuffer.data(idx_kin_gi + 309);

    auto tk_yyyz_xxxxxz = pbuffer.data(idx_kin_gi + 310);

    auto tk_yyyz_xxxxyy = pbuffer.data(idx_kin_gi + 311);

    auto tk_yyyz_xxxxyz = pbuffer.data(idx_kin_gi + 312);

    auto tk_yyyz_xxxxzz = pbuffer.data(idx_kin_gi + 313);

    auto tk_yyyz_xxxyyy = pbuffer.data(idx_kin_gi + 314);

    auto tk_yyyz_xxxyyz = pbuffer.data(idx_kin_gi + 315);

    auto tk_yyyz_xxxyzz = pbuffer.data(idx_kin_gi + 316);

    auto tk_yyyz_xxxzzz = pbuffer.data(idx_kin_gi + 317);

    auto tk_yyyz_xxyyyy = pbuffer.data(idx_kin_gi + 318);

    auto tk_yyyz_xxyyyz = pbuffer.data(idx_kin_gi + 319);

    auto tk_yyyz_xxyyzz = pbuffer.data(idx_kin_gi + 320);

    auto tk_yyyz_xxyzzz = pbuffer.data(idx_kin_gi + 321);

    auto tk_yyyz_xxzzzz = pbuffer.data(idx_kin_gi + 322);

    auto tk_yyyz_xyyyyy = pbuffer.data(idx_kin_gi + 323);

    auto tk_yyyz_xyyyyz = pbuffer.data(idx_kin_gi + 324);

    auto tk_yyyz_xyyyzz = pbuffer.data(idx_kin_gi + 325);

    auto tk_yyyz_xyyzzz = pbuffer.data(idx_kin_gi + 326);

    auto tk_yyyz_xyzzzz = pbuffer.data(idx_kin_gi + 327);

    auto tk_yyyz_xzzzzz = pbuffer.data(idx_kin_gi + 328);

    auto tk_yyyz_yyyyyy = pbuffer.data(idx_kin_gi + 329);

    auto tk_yyyz_yyyyyz = pbuffer.data(idx_kin_gi + 330);

    auto tk_yyyz_yyyyzz = pbuffer.data(idx_kin_gi + 331);

    auto tk_yyyz_yyyzzz = pbuffer.data(idx_kin_gi + 332);

    auto tk_yyyz_yyzzzz = pbuffer.data(idx_kin_gi + 333);

    auto tk_yyyz_yzzzzz = pbuffer.data(idx_kin_gi + 334);

    auto tk_yyyz_zzzzzz = pbuffer.data(idx_kin_gi + 335);

#pragma omp simd aligned(pa_z,               \
                             tk_yyy_xxxxx,   \
                             tk_yyy_xxxxxx,  \
                             tk_yyy_xxxxxy,  \
                             tk_yyy_xxxxxz,  \
                             tk_yyy_xxxxy,   \
                             tk_yyy_xxxxyy,  \
                             tk_yyy_xxxxyz,  \
                             tk_yyy_xxxxz,   \
                             tk_yyy_xxxxzz,  \
                             tk_yyy_xxxyy,   \
                             tk_yyy_xxxyyy,  \
                             tk_yyy_xxxyyz,  \
                             tk_yyy_xxxyz,   \
                             tk_yyy_xxxyzz,  \
                             tk_yyy_xxxzz,   \
                             tk_yyy_xxxzzz,  \
                             tk_yyy_xxyyy,   \
                             tk_yyy_xxyyyy,  \
                             tk_yyy_xxyyyz,  \
                             tk_yyy_xxyyz,   \
                             tk_yyy_xxyyzz,  \
                             tk_yyy_xxyzz,   \
                             tk_yyy_xxyzzz,  \
                             tk_yyy_xxzzz,   \
                             tk_yyy_xxzzzz,  \
                             tk_yyy_xyyyy,   \
                             tk_yyy_xyyyyy,  \
                             tk_yyy_xyyyyz,  \
                             tk_yyy_xyyyz,   \
                             tk_yyy_xyyyzz,  \
                             tk_yyy_xyyzz,   \
                             tk_yyy_xyyzzz,  \
                             tk_yyy_xyzzz,   \
                             tk_yyy_xyzzzz,  \
                             tk_yyy_xzzzz,   \
                             tk_yyy_xzzzzz,  \
                             tk_yyy_yyyyy,   \
                             tk_yyy_yyyyyy,  \
                             tk_yyy_yyyyyz,  \
                             tk_yyy_yyyyz,   \
                             tk_yyy_yyyyzz,  \
                             tk_yyy_yyyzz,   \
                             tk_yyy_yyyzzz,  \
                             tk_yyy_yyzzz,   \
                             tk_yyy_yyzzzz,  \
                             tk_yyy_yzzzz,   \
                             tk_yyy_yzzzzz,  \
                             tk_yyy_zzzzz,   \
                             tk_yyy_zzzzzz,  \
                             tk_yyyz_xxxxxx, \
                             tk_yyyz_xxxxxy, \
                             tk_yyyz_xxxxxz, \
                             tk_yyyz_xxxxyy, \
                             tk_yyyz_xxxxyz, \
                             tk_yyyz_xxxxzz, \
                             tk_yyyz_xxxyyy, \
                             tk_yyyz_xxxyyz, \
                             tk_yyyz_xxxyzz, \
                             tk_yyyz_xxxzzz, \
                             tk_yyyz_xxyyyy, \
                             tk_yyyz_xxyyyz, \
                             tk_yyyz_xxyyzz, \
                             tk_yyyz_xxyzzz, \
                             tk_yyyz_xxzzzz, \
                             tk_yyyz_xyyyyy, \
                             tk_yyyz_xyyyyz, \
                             tk_yyyz_xyyyzz, \
                             tk_yyyz_xyyzzz, \
                             tk_yyyz_xyzzzz, \
                             tk_yyyz_xzzzzz, \
                             tk_yyyz_yyyyyy, \
                             tk_yyyz_yyyyyz, \
                             tk_yyyz_yyyyzz, \
                             tk_yyyz_yyyzzz, \
                             tk_yyyz_yyzzzz, \
                             tk_yyyz_yzzzzz, \
                             tk_yyyz_zzzzzz, \
                             ts_yyyz_xxxxxx, \
                             ts_yyyz_xxxxxy, \
                             ts_yyyz_xxxxxz, \
                             ts_yyyz_xxxxyy, \
                             ts_yyyz_xxxxyz, \
                             ts_yyyz_xxxxzz, \
                             ts_yyyz_xxxyyy, \
                             ts_yyyz_xxxyyz, \
                             ts_yyyz_xxxyzz, \
                             ts_yyyz_xxxzzz, \
                             ts_yyyz_xxyyyy, \
                             ts_yyyz_xxyyyz, \
                             ts_yyyz_xxyyzz, \
                             ts_yyyz_xxyzzz, \
                             ts_yyyz_xxzzzz, \
                             ts_yyyz_xyyyyy, \
                             ts_yyyz_xyyyyz, \
                             ts_yyyz_xyyyzz, \
                             ts_yyyz_xyyzzz, \
                             ts_yyyz_xyzzzz, \
                             ts_yyyz_xzzzzz, \
                             ts_yyyz_yyyyyy, \
                             ts_yyyz_yyyyyz, \
                             ts_yyyz_yyyyzz, \
                             ts_yyyz_yyyzzz, \
                             ts_yyyz_yyzzzz, \
                             ts_yyyz_yzzzzz, \
                             ts_yyyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyz_xxxxxx[i] = tk_yyy_xxxxxx[i] * pa_z[i] + 2.0 * ts_yyyz_xxxxxx[i] * fz_0;

        tk_yyyz_xxxxxy[i] = tk_yyy_xxxxxy[i] * pa_z[i] + 2.0 * ts_yyyz_xxxxxy[i] * fz_0;

        tk_yyyz_xxxxxz[i] = tk_yyy_xxxxx[i] * fe_0 + tk_yyy_xxxxxz[i] * pa_z[i] + 2.0 * ts_yyyz_xxxxxz[i] * fz_0;

        tk_yyyz_xxxxyy[i] = tk_yyy_xxxxyy[i] * pa_z[i] + 2.0 * ts_yyyz_xxxxyy[i] * fz_0;

        tk_yyyz_xxxxyz[i] = tk_yyy_xxxxy[i] * fe_0 + tk_yyy_xxxxyz[i] * pa_z[i] + 2.0 * ts_yyyz_xxxxyz[i] * fz_0;

        tk_yyyz_xxxxzz[i] = 2.0 * tk_yyy_xxxxz[i] * fe_0 + tk_yyy_xxxxzz[i] * pa_z[i] + 2.0 * ts_yyyz_xxxxzz[i] * fz_0;

        tk_yyyz_xxxyyy[i] = tk_yyy_xxxyyy[i] * pa_z[i] + 2.0 * ts_yyyz_xxxyyy[i] * fz_0;

        tk_yyyz_xxxyyz[i] = tk_yyy_xxxyy[i] * fe_0 + tk_yyy_xxxyyz[i] * pa_z[i] + 2.0 * ts_yyyz_xxxyyz[i] * fz_0;

        tk_yyyz_xxxyzz[i] = 2.0 * tk_yyy_xxxyz[i] * fe_0 + tk_yyy_xxxyzz[i] * pa_z[i] + 2.0 * ts_yyyz_xxxyzz[i] * fz_0;

        tk_yyyz_xxxzzz[i] = 3.0 * tk_yyy_xxxzz[i] * fe_0 + tk_yyy_xxxzzz[i] * pa_z[i] + 2.0 * ts_yyyz_xxxzzz[i] * fz_0;

        tk_yyyz_xxyyyy[i] = tk_yyy_xxyyyy[i] * pa_z[i] + 2.0 * ts_yyyz_xxyyyy[i] * fz_0;

        tk_yyyz_xxyyyz[i] = tk_yyy_xxyyy[i] * fe_0 + tk_yyy_xxyyyz[i] * pa_z[i] + 2.0 * ts_yyyz_xxyyyz[i] * fz_0;

        tk_yyyz_xxyyzz[i] = 2.0 * tk_yyy_xxyyz[i] * fe_0 + tk_yyy_xxyyzz[i] * pa_z[i] + 2.0 * ts_yyyz_xxyyzz[i] * fz_0;

        tk_yyyz_xxyzzz[i] = 3.0 * tk_yyy_xxyzz[i] * fe_0 + tk_yyy_xxyzzz[i] * pa_z[i] + 2.0 * ts_yyyz_xxyzzz[i] * fz_0;

        tk_yyyz_xxzzzz[i] = 4.0 * tk_yyy_xxzzz[i] * fe_0 + tk_yyy_xxzzzz[i] * pa_z[i] + 2.0 * ts_yyyz_xxzzzz[i] * fz_0;

        tk_yyyz_xyyyyy[i] = tk_yyy_xyyyyy[i] * pa_z[i] + 2.0 * ts_yyyz_xyyyyy[i] * fz_0;

        tk_yyyz_xyyyyz[i] = tk_yyy_xyyyy[i] * fe_0 + tk_yyy_xyyyyz[i] * pa_z[i] + 2.0 * ts_yyyz_xyyyyz[i] * fz_0;

        tk_yyyz_xyyyzz[i] = 2.0 * tk_yyy_xyyyz[i] * fe_0 + tk_yyy_xyyyzz[i] * pa_z[i] + 2.0 * ts_yyyz_xyyyzz[i] * fz_0;

        tk_yyyz_xyyzzz[i] = 3.0 * tk_yyy_xyyzz[i] * fe_0 + tk_yyy_xyyzzz[i] * pa_z[i] + 2.0 * ts_yyyz_xyyzzz[i] * fz_0;

        tk_yyyz_xyzzzz[i] = 4.0 * tk_yyy_xyzzz[i] * fe_0 + tk_yyy_xyzzzz[i] * pa_z[i] + 2.0 * ts_yyyz_xyzzzz[i] * fz_0;

        tk_yyyz_xzzzzz[i] = 5.0 * tk_yyy_xzzzz[i] * fe_0 + tk_yyy_xzzzzz[i] * pa_z[i] + 2.0 * ts_yyyz_xzzzzz[i] * fz_0;

        tk_yyyz_yyyyyy[i] = tk_yyy_yyyyyy[i] * pa_z[i] + 2.0 * ts_yyyz_yyyyyy[i] * fz_0;

        tk_yyyz_yyyyyz[i] = tk_yyy_yyyyy[i] * fe_0 + tk_yyy_yyyyyz[i] * pa_z[i] + 2.0 * ts_yyyz_yyyyyz[i] * fz_0;

        tk_yyyz_yyyyzz[i] = 2.0 * tk_yyy_yyyyz[i] * fe_0 + tk_yyy_yyyyzz[i] * pa_z[i] + 2.0 * ts_yyyz_yyyyzz[i] * fz_0;

        tk_yyyz_yyyzzz[i] = 3.0 * tk_yyy_yyyzz[i] * fe_0 + tk_yyy_yyyzzz[i] * pa_z[i] + 2.0 * ts_yyyz_yyyzzz[i] * fz_0;

        tk_yyyz_yyzzzz[i] = 4.0 * tk_yyy_yyzzz[i] * fe_0 + tk_yyy_yyzzzz[i] * pa_z[i] + 2.0 * ts_yyyz_yyzzzz[i] * fz_0;

        tk_yyyz_yzzzzz[i] = 5.0 * tk_yyy_yzzzz[i] * fe_0 + tk_yyy_yzzzzz[i] * pa_z[i] + 2.0 * ts_yyyz_yzzzzz[i] * fz_0;

        tk_yyyz_zzzzzz[i] = 6.0 * tk_yyy_zzzzz[i] * fe_0 + tk_yyy_zzzzzz[i] * pa_z[i] + 2.0 * ts_yyyz_zzzzzz[i] * fz_0;
    }

    // Set up 336-364 components of targeted buffer : GI

    auto tk_yyzz_xxxxxx = pbuffer.data(idx_kin_gi + 336);

    auto tk_yyzz_xxxxxy = pbuffer.data(idx_kin_gi + 337);

    auto tk_yyzz_xxxxxz = pbuffer.data(idx_kin_gi + 338);

    auto tk_yyzz_xxxxyy = pbuffer.data(idx_kin_gi + 339);

    auto tk_yyzz_xxxxyz = pbuffer.data(idx_kin_gi + 340);

    auto tk_yyzz_xxxxzz = pbuffer.data(idx_kin_gi + 341);

    auto tk_yyzz_xxxyyy = pbuffer.data(idx_kin_gi + 342);

    auto tk_yyzz_xxxyyz = pbuffer.data(idx_kin_gi + 343);

    auto tk_yyzz_xxxyzz = pbuffer.data(idx_kin_gi + 344);

    auto tk_yyzz_xxxzzz = pbuffer.data(idx_kin_gi + 345);

    auto tk_yyzz_xxyyyy = pbuffer.data(idx_kin_gi + 346);

    auto tk_yyzz_xxyyyz = pbuffer.data(idx_kin_gi + 347);

    auto tk_yyzz_xxyyzz = pbuffer.data(idx_kin_gi + 348);

    auto tk_yyzz_xxyzzz = pbuffer.data(idx_kin_gi + 349);

    auto tk_yyzz_xxzzzz = pbuffer.data(idx_kin_gi + 350);

    auto tk_yyzz_xyyyyy = pbuffer.data(idx_kin_gi + 351);

    auto tk_yyzz_xyyyyz = pbuffer.data(idx_kin_gi + 352);

    auto tk_yyzz_xyyyzz = pbuffer.data(idx_kin_gi + 353);

    auto tk_yyzz_xyyzzz = pbuffer.data(idx_kin_gi + 354);

    auto tk_yyzz_xyzzzz = pbuffer.data(idx_kin_gi + 355);

    auto tk_yyzz_xzzzzz = pbuffer.data(idx_kin_gi + 356);

    auto tk_yyzz_yyyyyy = pbuffer.data(idx_kin_gi + 357);

    auto tk_yyzz_yyyyyz = pbuffer.data(idx_kin_gi + 358);

    auto tk_yyzz_yyyyzz = pbuffer.data(idx_kin_gi + 359);

    auto tk_yyzz_yyyzzz = pbuffer.data(idx_kin_gi + 360);

    auto tk_yyzz_yyzzzz = pbuffer.data(idx_kin_gi + 361);

    auto tk_yyzz_yzzzzz = pbuffer.data(idx_kin_gi + 362);

    auto tk_yyzz_zzzzzz = pbuffer.data(idx_kin_gi + 363);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tk_yy_xxxxxy,   \
                             tk_yy_xxxxyy,   \
                             tk_yy_xxxyyy,   \
                             tk_yy_xxyyyy,   \
                             tk_yy_xyyyyy,   \
                             tk_yy_yyyyyy,   \
                             tk_yyz_xxxxxy,  \
                             tk_yyz_xxxxyy,  \
                             tk_yyz_xxxyyy,  \
                             tk_yyz_xxyyyy,  \
                             tk_yyz_xyyyyy,  \
                             tk_yyz_yyyyyy,  \
                             tk_yyzz_xxxxxx, \
                             tk_yyzz_xxxxxy, \
                             tk_yyzz_xxxxxz, \
                             tk_yyzz_xxxxyy, \
                             tk_yyzz_xxxxyz, \
                             tk_yyzz_xxxxzz, \
                             tk_yyzz_xxxyyy, \
                             tk_yyzz_xxxyyz, \
                             tk_yyzz_xxxyzz, \
                             tk_yyzz_xxxzzz, \
                             tk_yyzz_xxyyyy, \
                             tk_yyzz_xxyyyz, \
                             tk_yyzz_xxyyzz, \
                             tk_yyzz_xxyzzz, \
                             tk_yyzz_xxzzzz, \
                             tk_yyzz_xyyyyy, \
                             tk_yyzz_xyyyyz, \
                             tk_yyzz_xyyyzz, \
                             tk_yyzz_xyyzzz, \
                             tk_yyzz_xyzzzz, \
                             tk_yyzz_xzzzzz, \
                             tk_yyzz_yyyyyy, \
                             tk_yyzz_yyyyyz, \
                             tk_yyzz_yyyyzz, \
                             tk_yyzz_yyyzzz, \
                             tk_yyzz_yyzzzz, \
                             tk_yyzz_yzzzzz, \
                             tk_yyzz_zzzzzz, \
                             tk_yzz_xxxxxx,  \
                             tk_yzz_xxxxxz,  \
                             tk_yzz_xxxxyz,  \
                             tk_yzz_xxxxz,   \
                             tk_yzz_xxxxzz,  \
                             tk_yzz_xxxyyz,  \
                             tk_yzz_xxxyz,   \
                             tk_yzz_xxxyzz,  \
                             tk_yzz_xxxzz,   \
                             tk_yzz_xxxzzz,  \
                             tk_yzz_xxyyyz,  \
                             tk_yzz_xxyyz,   \
                             tk_yzz_xxyyzz,  \
                             tk_yzz_xxyzz,   \
                             tk_yzz_xxyzzz,  \
                             tk_yzz_xxzzz,   \
                             tk_yzz_xxzzzz,  \
                             tk_yzz_xyyyyz,  \
                             tk_yzz_xyyyz,   \
                             tk_yzz_xyyyzz,  \
                             tk_yzz_xyyzz,   \
                             tk_yzz_xyyzzz,  \
                             tk_yzz_xyzzz,   \
                             tk_yzz_xyzzzz,  \
                             tk_yzz_xzzzz,   \
                             tk_yzz_xzzzzz,  \
                             tk_yzz_yyyyyz,  \
                             tk_yzz_yyyyz,   \
                             tk_yzz_yyyyzz,  \
                             tk_yzz_yyyzz,   \
                             tk_yzz_yyyzzz,  \
                             tk_yzz_yyzzz,   \
                             tk_yzz_yyzzzz,  \
                             tk_yzz_yzzzz,   \
                             tk_yzz_yzzzzz,  \
                             tk_yzz_zzzzz,   \
                             tk_yzz_zzzzzz,  \
                             tk_zz_xxxxxx,   \
                             tk_zz_xxxxxz,   \
                             tk_zz_xxxxyz,   \
                             tk_zz_xxxxzz,   \
                             tk_zz_xxxyyz,   \
                             tk_zz_xxxyzz,   \
                             tk_zz_xxxzzz,   \
                             tk_zz_xxyyyz,   \
                             tk_zz_xxyyzz,   \
                             tk_zz_xxyzzz,   \
                             tk_zz_xxzzzz,   \
                             tk_zz_xyyyyz,   \
                             tk_zz_xyyyzz,   \
                             tk_zz_xyyzzz,   \
                             tk_zz_xyzzzz,   \
                             tk_zz_xzzzzz,   \
                             tk_zz_yyyyyz,   \
                             tk_zz_yyyyzz,   \
                             tk_zz_yyyzzz,   \
                             tk_zz_yyzzzz,   \
                             tk_zz_yzzzzz,   \
                             tk_zz_zzzzzz,   \
                             ts_yy_xxxxxy,   \
                             ts_yy_xxxxyy,   \
                             ts_yy_xxxyyy,   \
                             ts_yy_xxyyyy,   \
                             ts_yy_xyyyyy,   \
                             ts_yy_yyyyyy,   \
                             ts_yyzz_xxxxxx, \
                             ts_yyzz_xxxxxy, \
                             ts_yyzz_xxxxxz, \
                             ts_yyzz_xxxxyy, \
                             ts_yyzz_xxxxyz, \
                             ts_yyzz_xxxxzz, \
                             ts_yyzz_xxxyyy, \
                             ts_yyzz_xxxyyz, \
                             ts_yyzz_xxxyzz, \
                             ts_yyzz_xxxzzz, \
                             ts_yyzz_xxyyyy, \
                             ts_yyzz_xxyyyz, \
                             ts_yyzz_xxyyzz, \
                             ts_yyzz_xxyzzz, \
                             ts_yyzz_xxzzzz, \
                             ts_yyzz_xyyyyy, \
                             ts_yyzz_xyyyyz, \
                             ts_yyzz_xyyyzz, \
                             ts_yyzz_xyyzzz, \
                             ts_yyzz_xyzzzz, \
                             ts_yyzz_xzzzzz, \
                             ts_yyzz_yyyyyy, \
                             ts_yyzz_yyyyyz, \
                             ts_yyzz_yyyyzz, \
                             ts_yyzz_yyyzzz, \
                             ts_yyzz_yyzzzz, \
                             ts_yyzz_yzzzzz, \
                             ts_yyzz_zzzzzz, \
                             ts_zz_xxxxxx,   \
                             ts_zz_xxxxxz,   \
                             ts_zz_xxxxyz,   \
                             ts_zz_xxxxzz,   \
                             ts_zz_xxxyyz,   \
                             ts_zz_xxxyzz,   \
                             ts_zz_xxxzzz,   \
                             ts_zz_xxyyyz,   \
                             ts_zz_xxyyzz,   \
                             ts_zz_xxyzzz,   \
                             ts_zz_xxzzzz,   \
                             ts_zz_xyyyyz,   \
                             ts_zz_xyyyzz,   \
                             ts_zz_xyyzzz,   \
                             ts_zz_xyzzzz,   \
                             ts_zz_xzzzzz,   \
                             ts_zz_yyyyyz,   \
                             ts_zz_yyyyzz,   \
                             ts_zz_yyyzzz,   \
                             ts_zz_yyzzzz,   \
                             ts_zz_yzzzzz,   \
                             ts_zz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzz_xxxxxx[i] =
            -2.0 * ts_zz_xxxxxx[i] * fbe_0 * fz_0 + tk_zz_xxxxxx[i] * fe_0 + tk_yzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_yyzz_xxxxxx[i] * fz_0;

        tk_yyzz_xxxxxy[i] =
            -2.0 * ts_yy_xxxxxy[i] * fbe_0 * fz_0 + tk_yy_xxxxxy[i] * fe_0 + tk_yyz_xxxxxy[i] * pa_z[i] + 2.0 * ts_yyzz_xxxxxy[i] * fz_0;

        tk_yyzz_xxxxxz[i] =
            -2.0 * ts_zz_xxxxxz[i] * fbe_0 * fz_0 + tk_zz_xxxxxz[i] * fe_0 + tk_yzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_yyzz_xxxxxz[i] * fz_0;

        tk_yyzz_xxxxyy[i] =
            -2.0 * ts_yy_xxxxyy[i] * fbe_0 * fz_0 + tk_yy_xxxxyy[i] * fe_0 + tk_yyz_xxxxyy[i] * pa_z[i] + 2.0 * ts_yyzz_xxxxyy[i] * fz_0;

        tk_yyzz_xxxxyz[i] = -2.0 * ts_zz_xxxxyz[i] * fbe_0 * fz_0 + tk_zz_xxxxyz[i] * fe_0 + tk_yzz_xxxxz[i] * fe_0 + tk_yzz_xxxxyz[i] * pa_y[i] +
                            2.0 * ts_yyzz_xxxxyz[i] * fz_0;

        tk_yyzz_xxxxzz[i] =
            -2.0 * ts_zz_xxxxzz[i] * fbe_0 * fz_0 + tk_zz_xxxxzz[i] * fe_0 + tk_yzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_yyzz_xxxxzz[i] * fz_0;

        tk_yyzz_xxxyyy[i] =
            -2.0 * ts_yy_xxxyyy[i] * fbe_0 * fz_0 + tk_yy_xxxyyy[i] * fe_0 + tk_yyz_xxxyyy[i] * pa_z[i] + 2.0 * ts_yyzz_xxxyyy[i] * fz_0;

        tk_yyzz_xxxyyz[i] = -2.0 * ts_zz_xxxyyz[i] * fbe_0 * fz_0 + tk_zz_xxxyyz[i] * fe_0 + 2.0 * tk_yzz_xxxyz[i] * fe_0 +
                            tk_yzz_xxxyyz[i] * pa_y[i] + 2.0 * ts_yyzz_xxxyyz[i] * fz_0;

        tk_yyzz_xxxyzz[i] = -2.0 * ts_zz_xxxyzz[i] * fbe_0 * fz_0 + tk_zz_xxxyzz[i] * fe_0 + tk_yzz_xxxzz[i] * fe_0 + tk_yzz_xxxyzz[i] * pa_y[i] +
                            2.0 * ts_yyzz_xxxyzz[i] * fz_0;

        tk_yyzz_xxxzzz[i] =
            -2.0 * ts_zz_xxxzzz[i] * fbe_0 * fz_0 + tk_zz_xxxzzz[i] * fe_0 + tk_yzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_yyzz_xxxzzz[i] * fz_0;

        tk_yyzz_xxyyyy[i] =
            -2.0 * ts_yy_xxyyyy[i] * fbe_0 * fz_0 + tk_yy_xxyyyy[i] * fe_0 + tk_yyz_xxyyyy[i] * pa_z[i] + 2.0 * ts_yyzz_xxyyyy[i] * fz_0;

        tk_yyzz_xxyyyz[i] = -2.0 * ts_zz_xxyyyz[i] * fbe_0 * fz_0 + tk_zz_xxyyyz[i] * fe_0 + 3.0 * tk_yzz_xxyyz[i] * fe_0 +
                            tk_yzz_xxyyyz[i] * pa_y[i] + 2.0 * ts_yyzz_xxyyyz[i] * fz_0;

        tk_yyzz_xxyyzz[i] = -2.0 * ts_zz_xxyyzz[i] * fbe_0 * fz_0 + tk_zz_xxyyzz[i] * fe_0 + 2.0 * tk_yzz_xxyzz[i] * fe_0 +
                            tk_yzz_xxyyzz[i] * pa_y[i] + 2.0 * ts_yyzz_xxyyzz[i] * fz_0;

        tk_yyzz_xxyzzz[i] = -2.0 * ts_zz_xxyzzz[i] * fbe_0 * fz_0 + tk_zz_xxyzzz[i] * fe_0 + tk_yzz_xxzzz[i] * fe_0 + tk_yzz_xxyzzz[i] * pa_y[i] +
                            2.0 * ts_yyzz_xxyzzz[i] * fz_0;

        tk_yyzz_xxzzzz[i] =
            -2.0 * ts_zz_xxzzzz[i] * fbe_0 * fz_0 + tk_zz_xxzzzz[i] * fe_0 + tk_yzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_yyzz_xxzzzz[i] * fz_0;

        tk_yyzz_xyyyyy[i] =
            -2.0 * ts_yy_xyyyyy[i] * fbe_0 * fz_0 + tk_yy_xyyyyy[i] * fe_0 + tk_yyz_xyyyyy[i] * pa_z[i] + 2.0 * ts_yyzz_xyyyyy[i] * fz_0;

        tk_yyzz_xyyyyz[i] = -2.0 * ts_zz_xyyyyz[i] * fbe_0 * fz_0 + tk_zz_xyyyyz[i] * fe_0 + 4.0 * tk_yzz_xyyyz[i] * fe_0 +
                            tk_yzz_xyyyyz[i] * pa_y[i] + 2.0 * ts_yyzz_xyyyyz[i] * fz_0;

        tk_yyzz_xyyyzz[i] = -2.0 * ts_zz_xyyyzz[i] * fbe_0 * fz_0 + tk_zz_xyyyzz[i] * fe_0 + 3.0 * tk_yzz_xyyzz[i] * fe_0 +
                            tk_yzz_xyyyzz[i] * pa_y[i] + 2.0 * ts_yyzz_xyyyzz[i] * fz_0;

        tk_yyzz_xyyzzz[i] = -2.0 * ts_zz_xyyzzz[i] * fbe_0 * fz_0 + tk_zz_xyyzzz[i] * fe_0 + 2.0 * tk_yzz_xyzzz[i] * fe_0 +
                            tk_yzz_xyyzzz[i] * pa_y[i] + 2.0 * ts_yyzz_xyyzzz[i] * fz_0;

        tk_yyzz_xyzzzz[i] = -2.0 * ts_zz_xyzzzz[i] * fbe_0 * fz_0 + tk_zz_xyzzzz[i] * fe_0 + tk_yzz_xzzzz[i] * fe_0 + tk_yzz_xyzzzz[i] * pa_y[i] +
                            2.0 * ts_yyzz_xyzzzz[i] * fz_0;

        tk_yyzz_xzzzzz[i] =
            -2.0 * ts_zz_xzzzzz[i] * fbe_0 * fz_0 + tk_zz_xzzzzz[i] * fe_0 + tk_yzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_yyzz_xzzzzz[i] * fz_0;

        tk_yyzz_yyyyyy[i] =
            -2.0 * ts_yy_yyyyyy[i] * fbe_0 * fz_0 + tk_yy_yyyyyy[i] * fe_0 + tk_yyz_yyyyyy[i] * pa_z[i] + 2.0 * ts_yyzz_yyyyyy[i] * fz_0;

        tk_yyzz_yyyyyz[i] = -2.0 * ts_zz_yyyyyz[i] * fbe_0 * fz_0 + tk_zz_yyyyyz[i] * fe_0 + 5.0 * tk_yzz_yyyyz[i] * fe_0 +
                            tk_yzz_yyyyyz[i] * pa_y[i] + 2.0 * ts_yyzz_yyyyyz[i] * fz_0;

        tk_yyzz_yyyyzz[i] = -2.0 * ts_zz_yyyyzz[i] * fbe_0 * fz_0 + tk_zz_yyyyzz[i] * fe_0 + 4.0 * tk_yzz_yyyzz[i] * fe_0 +
                            tk_yzz_yyyyzz[i] * pa_y[i] + 2.0 * ts_yyzz_yyyyzz[i] * fz_0;

        tk_yyzz_yyyzzz[i] = -2.0 * ts_zz_yyyzzz[i] * fbe_0 * fz_0 + tk_zz_yyyzzz[i] * fe_0 + 3.0 * tk_yzz_yyzzz[i] * fe_0 +
                            tk_yzz_yyyzzz[i] * pa_y[i] + 2.0 * ts_yyzz_yyyzzz[i] * fz_0;

        tk_yyzz_yyzzzz[i] = -2.0 * ts_zz_yyzzzz[i] * fbe_0 * fz_0 + tk_zz_yyzzzz[i] * fe_0 + 2.0 * tk_yzz_yzzzz[i] * fe_0 +
                            tk_yzz_yyzzzz[i] * pa_y[i] + 2.0 * ts_yyzz_yyzzzz[i] * fz_0;

        tk_yyzz_yzzzzz[i] = -2.0 * ts_zz_yzzzzz[i] * fbe_0 * fz_0 + tk_zz_yzzzzz[i] * fe_0 + tk_yzz_zzzzz[i] * fe_0 + tk_yzz_yzzzzz[i] * pa_y[i] +
                            2.0 * ts_yyzz_yzzzzz[i] * fz_0;

        tk_yyzz_zzzzzz[i] =
            -2.0 * ts_zz_zzzzzz[i] * fbe_0 * fz_0 + tk_zz_zzzzzz[i] * fe_0 + tk_yzz_zzzzzz[i] * pa_y[i] + 2.0 * ts_yyzz_zzzzzz[i] * fz_0;
    }

    // Set up 364-392 components of targeted buffer : GI

    auto tk_yzzz_xxxxxx = pbuffer.data(idx_kin_gi + 364);

    auto tk_yzzz_xxxxxy = pbuffer.data(idx_kin_gi + 365);

    auto tk_yzzz_xxxxxz = pbuffer.data(idx_kin_gi + 366);

    auto tk_yzzz_xxxxyy = pbuffer.data(idx_kin_gi + 367);

    auto tk_yzzz_xxxxyz = pbuffer.data(idx_kin_gi + 368);

    auto tk_yzzz_xxxxzz = pbuffer.data(idx_kin_gi + 369);

    auto tk_yzzz_xxxyyy = pbuffer.data(idx_kin_gi + 370);

    auto tk_yzzz_xxxyyz = pbuffer.data(idx_kin_gi + 371);

    auto tk_yzzz_xxxyzz = pbuffer.data(idx_kin_gi + 372);

    auto tk_yzzz_xxxzzz = pbuffer.data(idx_kin_gi + 373);

    auto tk_yzzz_xxyyyy = pbuffer.data(idx_kin_gi + 374);

    auto tk_yzzz_xxyyyz = pbuffer.data(idx_kin_gi + 375);

    auto tk_yzzz_xxyyzz = pbuffer.data(idx_kin_gi + 376);

    auto tk_yzzz_xxyzzz = pbuffer.data(idx_kin_gi + 377);

    auto tk_yzzz_xxzzzz = pbuffer.data(idx_kin_gi + 378);

    auto tk_yzzz_xyyyyy = pbuffer.data(idx_kin_gi + 379);

    auto tk_yzzz_xyyyyz = pbuffer.data(idx_kin_gi + 380);

    auto tk_yzzz_xyyyzz = pbuffer.data(idx_kin_gi + 381);

    auto tk_yzzz_xyyzzz = pbuffer.data(idx_kin_gi + 382);

    auto tk_yzzz_xyzzzz = pbuffer.data(idx_kin_gi + 383);

    auto tk_yzzz_xzzzzz = pbuffer.data(idx_kin_gi + 384);

    auto tk_yzzz_yyyyyy = pbuffer.data(idx_kin_gi + 385);

    auto tk_yzzz_yyyyyz = pbuffer.data(idx_kin_gi + 386);

    auto tk_yzzz_yyyyzz = pbuffer.data(idx_kin_gi + 387);

    auto tk_yzzz_yyyzzz = pbuffer.data(idx_kin_gi + 388);

    auto tk_yzzz_yyzzzz = pbuffer.data(idx_kin_gi + 389);

    auto tk_yzzz_yzzzzz = pbuffer.data(idx_kin_gi + 390);

    auto tk_yzzz_zzzzzz = pbuffer.data(idx_kin_gi + 391);

#pragma omp simd aligned(pa_y,               \
                             tk_yzzz_xxxxxx, \
                             tk_yzzz_xxxxxy, \
                             tk_yzzz_xxxxxz, \
                             tk_yzzz_xxxxyy, \
                             tk_yzzz_xxxxyz, \
                             tk_yzzz_xxxxzz, \
                             tk_yzzz_xxxyyy, \
                             tk_yzzz_xxxyyz, \
                             tk_yzzz_xxxyzz, \
                             tk_yzzz_xxxzzz, \
                             tk_yzzz_xxyyyy, \
                             tk_yzzz_xxyyyz, \
                             tk_yzzz_xxyyzz, \
                             tk_yzzz_xxyzzz, \
                             tk_yzzz_xxzzzz, \
                             tk_yzzz_xyyyyy, \
                             tk_yzzz_xyyyyz, \
                             tk_yzzz_xyyyzz, \
                             tk_yzzz_xyyzzz, \
                             tk_yzzz_xyzzzz, \
                             tk_yzzz_xzzzzz, \
                             tk_yzzz_yyyyyy, \
                             tk_yzzz_yyyyyz, \
                             tk_yzzz_yyyyzz, \
                             tk_yzzz_yyyzzz, \
                             tk_yzzz_yyzzzz, \
                             tk_yzzz_yzzzzz, \
                             tk_yzzz_zzzzzz, \
                             tk_zzz_xxxxx,   \
                             tk_zzz_xxxxxx,  \
                             tk_zzz_xxxxxy,  \
                             tk_zzz_xxxxxz,  \
                             tk_zzz_xxxxy,   \
                             tk_zzz_xxxxyy,  \
                             tk_zzz_xxxxyz,  \
                             tk_zzz_xxxxz,   \
                             tk_zzz_xxxxzz,  \
                             tk_zzz_xxxyy,   \
                             tk_zzz_xxxyyy,  \
                             tk_zzz_xxxyyz,  \
                             tk_zzz_xxxyz,   \
                             tk_zzz_xxxyzz,  \
                             tk_zzz_xxxzz,   \
                             tk_zzz_xxxzzz,  \
                             tk_zzz_xxyyy,   \
                             tk_zzz_xxyyyy,  \
                             tk_zzz_xxyyyz,  \
                             tk_zzz_xxyyz,   \
                             tk_zzz_xxyyzz,  \
                             tk_zzz_xxyzz,   \
                             tk_zzz_xxyzzz,  \
                             tk_zzz_xxzzz,   \
                             tk_zzz_xxzzzz,  \
                             tk_zzz_xyyyy,   \
                             tk_zzz_xyyyyy,  \
                             tk_zzz_xyyyyz,  \
                             tk_zzz_xyyyz,   \
                             tk_zzz_xyyyzz,  \
                             tk_zzz_xyyzz,   \
                             tk_zzz_xyyzzz,  \
                             tk_zzz_xyzzz,   \
                             tk_zzz_xyzzzz,  \
                             tk_zzz_xzzzz,   \
                             tk_zzz_xzzzzz,  \
                             tk_zzz_yyyyy,   \
                             tk_zzz_yyyyyy,  \
                             tk_zzz_yyyyyz,  \
                             tk_zzz_yyyyz,   \
                             tk_zzz_yyyyzz,  \
                             tk_zzz_yyyzz,   \
                             tk_zzz_yyyzzz,  \
                             tk_zzz_yyzzz,   \
                             tk_zzz_yyzzzz,  \
                             tk_zzz_yzzzz,   \
                             tk_zzz_yzzzzz,  \
                             tk_zzz_zzzzz,   \
                             tk_zzz_zzzzzz,  \
                             ts_yzzz_xxxxxx, \
                             ts_yzzz_xxxxxy, \
                             ts_yzzz_xxxxxz, \
                             ts_yzzz_xxxxyy, \
                             ts_yzzz_xxxxyz, \
                             ts_yzzz_xxxxzz, \
                             ts_yzzz_xxxyyy, \
                             ts_yzzz_xxxyyz, \
                             ts_yzzz_xxxyzz, \
                             ts_yzzz_xxxzzz, \
                             ts_yzzz_xxyyyy, \
                             ts_yzzz_xxyyyz, \
                             ts_yzzz_xxyyzz, \
                             ts_yzzz_xxyzzz, \
                             ts_yzzz_xxzzzz, \
                             ts_yzzz_xyyyyy, \
                             ts_yzzz_xyyyyz, \
                             ts_yzzz_xyyyzz, \
                             ts_yzzz_xyyzzz, \
                             ts_yzzz_xyzzzz, \
                             ts_yzzz_xzzzzz, \
                             ts_yzzz_yyyyyy, \
                             ts_yzzz_yyyyyz, \
                             ts_yzzz_yyyyzz, \
                             ts_yzzz_yyyzzz, \
                             ts_yzzz_yyzzzz, \
                             ts_yzzz_yzzzzz, \
                             ts_yzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzz_xxxxxx[i] = tk_zzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_yzzz_xxxxxx[i] * fz_0;

        tk_yzzz_xxxxxy[i] = tk_zzz_xxxxx[i] * fe_0 + tk_zzz_xxxxxy[i] * pa_y[i] + 2.0 * ts_yzzz_xxxxxy[i] * fz_0;

        tk_yzzz_xxxxxz[i] = tk_zzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_yzzz_xxxxxz[i] * fz_0;

        tk_yzzz_xxxxyy[i] = 2.0 * tk_zzz_xxxxy[i] * fe_0 + tk_zzz_xxxxyy[i] * pa_y[i] + 2.0 * ts_yzzz_xxxxyy[i] * fz_0;

        tk_yzzz_xxxxyz[i] = tk_zzz_xxxxz[i] * fe_0 + tk_zzz_xxxxyz[i] * pa_y[i] + 2.0 * ts_yzzz_xxxxyz[i] * fz_0;

        tk_yzzz_xxxxzz[i] = tk_zzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_yzzz_xxxxzz[i] * fz_0;

        tk_yzzz_xxxyyy[i] = 3.0 * tk_zzz_xxxyy[i] * fe_0 + tk_zzz_xxxyyy[i] * pa_y[i] + 2.0 * ts_yzzz_xxxyyy[i] * fz_0;

        tk_yzzz_xxxyyz[i] = 2.0 * tk_zzz_xxxyz[i] * fe_0 + tk_zzz_xxxyyz[i] * pa_y[i] + 2.0 * ts_yzzz_xxxyyz[i] * fz_0;

        tk_yzzz_xxxyzz[i] = tk_zzz_xxxzz[i] * fe_0 + tk_zzz_xxxyzz[i] * pa_y[i] + 2.0 * ts_yzzz_xxxyzz[i] * fz_0;

        tk_yzzz_xxxzzz[i] = tk_zzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_yzzz_xxxzzz[i] * fz_0;

        tk_yzzz_xxyyyy[i] = 4.0 * tk_zzz_xxyyy[i] * fe_0 + tk_zzz_xxyyyy[i] * pa_y[i] + 2.0 * ts_yzzz_xxyyyy[i] * fz_0;

        tk_yzzz_xxyyyz[i] = 3.0 * tk_zzz_xxyyz[i] * fe_0 + tk_zzz_xxyyyz[i] * pa_y[i] + 2.0 * ts_yzzz_xxyyyz[i] * fz_0;

        tk_yzzz_xxyyzz[i] = 2.0 * tk_zzz_xxyzz[i] * fe_0 + tk_zzz_xxyyzz[i] * pa_y[i] + 2.0 * ts_yzzz_xxyyzz[i] * fz_0;

        tk_yzzz_xxyzzz[i] = tk_zzz_xxzzz[i] * fe_0 + tk_zzz_xxyzzz[i] * pa_y[i] + 2.0 * ts_yzzz_xxyzzz[i] * fz_0;

        tk_yzzz_xxzzzz[i] = tk_zzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_yzzz_xxzzzz[i] * fz_0;

        tk_yzzz_xyyyyy[i] = 5.0 * tk_zzz_xyyyy[i] * fe_0 + tk_zzz_xyyyyy[i] * pa_y[i] + 2.0 * ts_yzzz_xyyyyy[i] * fz_0;

        tk_yzzz_xyyyyz[i] = 4.0 * tk_zzz_xyyyz[i] * fe_0 + tk_zzz_xyyyyz[i] * pa_y[i] + 2.0 * ts_yzzz_xyyyyz[i] * fz_0;

        tk_yzzz_xyyyzz[i] = 3.0 * tk_zzz_xyyzz[i] * fe_0 + tk_zzz_xyyyzz[i] * pa_y[i] + 2.0 * ts_yzzz_xyyyzz[i] * fz_0;

        tk_yzzz_xyyzzz[i] = 2.0 * tk_zzz_xyzzz[i] * fe_0 + tk_zzz_xyyzzz[i] * pa_y[i] + 2.0 * ts_yzzz_xyyzzz[i] * fz_0;

        tk_yzzz_xyzzzz[i] = tk_zzz_xzzzz[i] * fe_0 + tk_zzz_xyzzzz[i] * pa_y[i] + 2.0 * ts_yzzz_xyzzzz[i] * fz_0;

        tk_yzzz_xzzzzz[i] = tk_zzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_yzzz_xzzzzz[i] * fz_0;

        tk_yzzz_yyyyyy[i] = 6.0 * tk_zzz_yyyyy[i] * fe_0 + tk_zzz_yyyyyy[i] * pa_y[i] + 2.0 * ts_yzzz_yyyyyy[i] * fz_0;

        tk_yzzz_yyyyyz[i] = 5.0 * tk_zzz_yyyyz[i] * fe_0 + tk_zzz_yyyyyz[i] * pa_y[i] + 2.0 * ts_yzzz_yyyyyz[i] * fz_0;

        tk_yzzz_yyyyzz[i] = 4.0 * tk_zzz_yyyzz[i] * fe_0 + tk_zzz_yyyyzz[i] * pa_y[i] + 2.0 * ts_yzzz_yyyyzz[i] * fz_0;

        tk_yzzz_yyyzzz[i] = 3.0 * tk_zzz_yyzzz[i] * fe_0 + tk_zzz_yyyzzz[i] * pa_y[i] + 2.0 * ts_yzzz_yyyzzz[i] * fz_0;

        tk_yzzz_yyzzzz[i] = 2.0 * tk_zzz_yzzzz[i] * fe_0 + tk_zzz_yyzzzz[i] * pa_y[i] + 2.0 * ts_yzzz_yyzzzz[i] * fz_0;

        tk_yzzz_yzzzzz[i] = tk_zzz_zzzzz[i] * fe_0 + tk_zzz_yzzzzz[i] * pa_y[i] + 2.0 * ts_yzzz_yzzzzz[i] * fz_0;

        tk_yzzz_zzzzzz[i] = tk_zzz_zzzzzz[i] * pa_y[i] + 2.0 * ts_yzzz_zzzzzz[i] * fz_0;
    }

    // Set up 392-420 components of targeted buffer : GI

    auto tk_zzzz_xxxxxx = pbuffer.data(idx_kin_gi + 392);

    auto tk_zzzz_xxxxxy = pbuffer.data(idx_kin_gi + 393);

    auto tk_zzzz_xxxxxz = pbuffer.data(idx_kin_gi + 394);

    auto tk_zzzz_xxxxyy = pbuffer.data(idx_kin_gi + 395);

    auto tk_zzzz_xxxxyz = pbuffer.data(idx_kin_gi + 396);

    auto tk_zzzz_xxxxzz = pbuffer.data(idx_kin_gi + 397);

    auto tk_zzzz_xxxyyy = pbuffer.data(idx_kin_gi + 398);

    auto tk_zzzz_xxxyyz = pbuffer.data(idx_kin_gi + 399);

    auto tk_zzzz_xxxyzz = pbuffer.data(idx_kin_gi + 400);

    auto tk_zzzz_xxxzzz = pbuffer.data(idx_kin_gi + 401);

    auto tk_zzzz_xxyyyy = pbuffer.data(idx_kin_gi + 402);

    auto tk_zzzz_xxyyyz = pbuffer.data(idx_kin_gi + 403);

    auto tk_zzzz_xxyyzz = pbuffer.data(idx_kin_gi + 404);

    auto tk_zzzz_xxyzzz = pbuffer.data(idx_kin_gi + 405);

    auto tk_zzzz_xxzzzz = pbuffer.data(idx_kin_gi + 406);

    auto tk_zzzz_xyyyyy = pbuffer.data(idx_kin_gi + 407);

    auto tk_zzzz_xyyyyz = pbuffer.data(idx_kin_gi + 408);

    auto tk_zzzz_xyyyzz = pbuffer.data(idx_kin_gi + 409);

    auto tk_zzzz_xyyzzz = pbuffer.data(idx_kin_gi + 410);

    auto tk_zzzz_xyzzzz = pbuffer.data(idx_kin_gi + 411);

    auto tk_zzzz_xzzzzz = pbuffer.data(idx_kin_gi + 412);

    auto tk_zzzz_yyyyyy = pbuffer.data(idx_kin_gi + 413);

    auto tk_zzzz_yyyyyz = pbuffer.data(idx_kin_gi + 414);

    auto tk_zzzz_yyyyzz = pbuffer.data(idx_kin_gi + 415);

    auto tk_zzzz_yyyzzz = pbuffer.data(idx_kin_gi + 416);

    auto tk_zzzz_yyzzzz = pbuffer.data(idx_kin_gi + 417);

    auto tk_zzzz_yzzzzz = pbuffer.data(idx_kin_gi + 418);

    auto tk_zzzz_zzzzzz = pbuffer.data(idx_kin_gi + 419);

#pragma omp simd aligned(pa_z,               \
                             tk_zz_xxxxxx,   \
                             tk_zz_xxxxxy,   \
                             tk_zz_xxxxxz,   \
                             tk_zz_xxxxyy,   \
                             tk_zz_xxxxyz,   \
                             tk_zz_xxxxzz,   \
                             tk_zz_xxxyyy,   \
                             tk_zz_xxxyyz,   \
                             tk_zz_xxxyzz,   \
                             tk_zz_xxxzzz,   \
                             tk_zz_xxyyyy,   \
                             tk_zz_xxyyyz,   \
                             tk_zz_xxyyzz,   \
                             tk_zz_xxyzzz,   \
                             tk_zz_xxzzzz,   \
                             tk_zz_xyyyyy,   \
                             tk_zz_xyyyyz,   \
                             tk_zz_xyyyzz,   \
                             tk_zz_xyyzzz,   \
                             tk_zz_xyzzzz,   \
                             tk_zz_xzzzzz,   \
                             tk_zz_yyyyyy,   \
                             tk_zz_yyyyyz,   \
                             tk_zz_yyyyzz,   \
                             tk_zz_yyyzzz,   \
                             tk_zz_yyzzzz,   \
                             tk_zz_yzzzzz,   \
                             tk_zz_zzzzzz,   \
                             tk_zzz_xxxxx,   \
                             tk_zzz_xxxxxx,  \
                             tk_zzz_xxxxxy,  \
                             tk_zzz_xxxxxz,  \
                             tk_zzz_xxxxy,   \
                             tk_zzz_xxxxyy,  \
                             tk_zzz_xxxxyz,  \
                             tk_zzz_xxxxz,   \
                             tk_zzz_xxxxzz,  \
                             tk_zzz_xxxyy,   \
                             tk_zzz_xxxyyy,  \
                             tk_zzz_xxxyyz,  \
                             tk_zzz_xxxyz,   \
                             tk_zzz_xxxyzz,  \
                             tk_zzz_xxxzz,   \
                             tk_zzz_xxxzzz,  \
                             tk_zzz_xxyyy,   \
                             tk_zzz_xxyyyy,  \
                             tk_zzz_xxyyyz,  \
                             tk_zzz_xxyyz,   \
                             tk_zzz_xxyyzz,  \
                             tk_zzz_xxyzz,   \
                             tk_zzz_xxyzzz,  \
                             tk_zzz_xxzzz,   \
                             tk_zzz_xxzzzz,  \
                             tk_zzz_xyyyy,   \
                             tk_zzz_xyyyyy,  \
                             tk_zzz_xyyyyz,  \
                             tk_zzz_xyyyz,   \
                             tk_zzz_xyyyzz,  \
                             tk_zzz_xyyzz,   \
                             tk_zzz_xyyzzz,  \
                             tk_zzz_xyzzz,   \
                             tk_zzz_xyzzzz,  \
                             tk_zzz_xzzzz,   \
                             tk_zzz_xzzzzz,  \
                             tk_zzz_yyyyy,   \
                             tk_zzz_yyyyyy,  \
                             tk_zzz_yyyyyz,  \
                             tk_zzz_yyyyz,   \
                             tk_zzz_yyyyzz,  \
                             tk_zzz_yyyzz,   \
                             tk_zzz_yyyzzz,  \
                             tk_zzz_yyzzz,   \
                             tk_zzz_yyzzzz,  \
                             tk_zzz_yzzzz,   \
                             tk_zzz_yzzzzz,  \
                             tk_zzz_zzzzz,   \
                             tk_zzz_zzzzzz,  \
                             tk_zzzz_xxxxxx, \
                             tk_zzzz_xxxxxy, \
                             tk_zzzz_xxxxxz, \
                             tk_zzzz_xxxxyy, \
                             tk_zzzz_xxxxyz, \
                             tk_zzzz_xxxxzz, \
                             tk_zzzz_xxxyyy, \
                             tk_zzzz_xxxyyz, \
                             tk_zzzz_xxxyzz, \
                             tk_zzzz_xxxzzz, \
                             tk_zzzz_xxyyyy, \
                             tk_zzzz_xxyyyz, \
                             tk_zzzz_xxyyzz, \
                             tk_zzzz_xxyzzz, \
                             tk_zzzz_xxzzzz, \
                             tk_zzzz_xyyyyy, \
                             tk_zzzz_xyyyyz, \
                             tk_zzzz_xyyyzz, \
                             tk_zzzz_xyyzzz, \
                             tk_zzzz_xyzzzz, \
                             tk_zzzz_xzzzzz, \
                             tk_zzzz_yyyyyy, \
                             tk_zzzz_yyyyyz, \
                             tk_zzzz_yyyyzz, \
                             tk_zzzz_yyyzzz, \
                             tk_zzzz_yyzzzz, \
                             tk_zzzz_yzzzzz, \
                             tk_zzzz_zzzzzz, \
                             ts_zz_xxxxxx,   \
                             ts_zz_xxxxxy,   \
                             ts_zz_xxxxxz,   \
                             ts_zz_xxxxyy,   \
                             ts_zz_xxxxyz,   \
                             ts_zz_xxxxzz,   \
                             ts_zz_xxxyyy,   \
                             ts_zz_xxxyyz,   \
                             ts_zz_xxxyzz,   \
                             ts_zz_xxxzzz,   \
                             ts_zz_xxyyyy,   \
                             ts_zz_xxyyyz,   \
                             ts_zz_xxyyzz,   \
                             ts_zz_xxyzzz,   \
                             ts_zz_xxzzzz,   \
                             ts_zz_xyyyyy,   \
                             ts_zz_xyyyyz,   \
                             ts_zz_xyyyzz,   \
                             ts_zz_xyyzzz,   \
                             ts_zz_xyzzzz,   \
                             ts_zz_xzzzzz,   \
                             ts_zz_yyyyyy,   \
                             ts_zz_yyyyyz,   \
                             ts_zz_yyyyzz,   \
                             ts_zz_yyyzzz,   \
                             ts_zz_yyzzzz,   \
                             ts_zz_yzzzzz,   \
                             ts_zz_zzzzzz,   \
                             ts_zzzz_xxxxxx, \
                             ts_zzzz_xxxxxy, \
                             ts_zzzz_xxxxxz, \
                             ts_zzzz_xxxxyy, \
                             ts_zzzz_xxxxyz, \
                             ts_zzzz_xxxxzz, \
                             ts_zzzz_xxxyyy, \
                             ts_zzzz_xxxyyz, \
                             ts_zzzz_xxxyzz, \
                             ts_zzzz_xxxzzz, \
                             ts_zzzz_xxyyyy, \
                             ts_zzzz_xxyyyz, \
                             ts_zzzz_xxyyzz, \
                             ts_zzzz_xxyzzz, \
                             ts_zzzz_xxzzzz, \
                             ts_zzzz_xyyyyy, \
                             ts_zzzz_xyyyyz, \
                             ts_zzzz_xyyyzz, \
                             ts_zzzz_xyyzzz, \
                             ts_zzzz_xyzzzz, \
                             ts_zzzz_xzzzzz, \
                             ts_zzzz_yyyyyy, \
                             ts_zzzz_yyyyyz, \
                             ts_zzzz_yyyyzz, \
                             ts_zzzz_yyyzzz, \
                             ts_zzzz_yyzzzz, \
                             ts_zzzz_yzzzzz, \
                             ts_zzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzz_xxxxxx[i] =
            -6.0 * ts_zz_xxxxxx[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxxxx[i] * fe_0 + tk_zzz_xxxxxx[i] * pa_z[i] + 2.0 * ts_zzzz_xxxxxx[i] * fz_0;

        tk_zzzz_xxxxxy[i] =
            -6.0 * ts_zz_xxxxxy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxxxy[i] * fe_0 + tk_zzz_xxxxxy[i] * pa_z[i] + 2.0 * ts_zzzz_xxxxxy[i] * fz_0;

        tk_zzzz_xxxxxz[i] = -6.0 * ts_zz_xxxxxz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxxxz[i] * fe_0 + tk_zzz_xxxxx[i] * fe_0 +
                            tk_zzz_xxxxxz[i] * pa_z[i] + 2.0 * ts_zzzz_xxxxxz[i] * fz_0;

        tk_zzzz_xxxxyy[i] =
            -6.0 * ts_zz_xxxxyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxxyy[i] * fe_0 + tk_zzz_xxxxyy[i] * pa_z[i] + 2.0 * ts_zzzz_xxxxyy[i] * fz_0;

        tk_zzzz_xxxxyz[i] = -6.0 * ts_zz_xxxxyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxxyz[i] * fe_0 + tk_zzz_xxxxy[i] * fe_0 +
                            tk_zzz_xxxxyz[i] * pa_z[i] + 2.0 * ts_zzzz_xxxxyz[i] * fz_0;

        tk_zzzz_xxxxzz[i] = -6.0 * ts_zz_xxxxzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxxzz[i] * fe_0 + 2.0 * tk_zzz_xxxxz[i] * fe_0 +
                            tk_zzz_xxxxzz[i] * pa_z[i] + 2.0 * ts_zzzz_xxxxzz[i] * fz_0;

        tk_zzzz_xxxyyy[i] =
            -6.0 * ts_zz_xxxyyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxyyy[i] * fe_0 + tk_zzz_xxxyyy[i] * pa_z[i] + 2.0 * ts_zzzz_xxxyyy[i] * fz_0;

        tk_zzzz_xxxyyz[i] = -6.0 * ts_zz_xxxyyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxyyz[i] * fe_0 + tk_zzz_xxxyy[i] * fe_0 +
                            tk_zzz_xxxyyz[i] * pa_z[i] + 2.0 * ts_zzzz_xxxyyz[i] * fz_0;

        tk_zzzz_xxxyzz[i] = -6.0 * ts_zz_xxxyzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxyzz[i] * fe_0 + 2.0 * tk_zzz_xxxyz[i] * fe_0 +
                            tk_zzz_xxxyzz[i] * pa_z[i] + 2.0 * ts_zzzz_xxxyzz[i] * fz_0;

        tk_zzzz_xxxzzz[i] = -6.0 * ts_zz_xxxzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxxzzz[i] * fe_0 + 3.0 * tk_zzz_xxxzz[i] * fe_0 +
                            tk_zzz_xxxzzz[i] * pa_z[i] + 2.0 * ts_zzzz_xxxzzz[i] * fz_0;

        tk_zzzz_xxyyyy[i] =
            -6.0 * ts_zz_xxyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxyyyy[i] * fe_0 + tk_zzz_xxyyyy[i] * pa_z[i] + 2.0 * ts_zzzz_xxyyyy[i] * fz_0;

        tk_zzzz_xxyyyz[i] = -6.0 * ts_zz_xxyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxyyyz[i] * fe_0 + tk_zzz_xxyyy[i] * fe_0 +
                            tk_zzz_xxyyyz[i] * pa_z[i] + 2.0 * ts_zzzz_xxyyyz[i] * fz_0;

        tk_zzzz_xxyyzz[i] = -6.0 * ts_zz_xxyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxyyzz[i] * fe_0 + 2.0 * tk_zzz_xxyyz[i] * fe_0 +
                            tk_zzz_xxyyzz[i] * pa_z[i] + 2.0 * ts_zzzz_xxyyzz[i] * fz_0;

        tk_zzzz_xxyzzz[i] = -6.0 * ts_zz_xxyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxyzzz[i] * fe_0 + 3.0 * tk_zzz_xxyzz[i] * fe_0 +
                            tk_zzz_xxyzzz[i] * pa_z[i] + 2.0 * ts_zzzz_xxyzzz[i] * fz_0;

        tk_zzzz_xxzzzz[i] = -6.0 * ts_zz_xxzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xxzzzz[i] * fe_0 + 4.0 * tk_zzz_xxzzz[i] * fe_0 +
                            tk_zzz_xxzzzz[i] * pa_z[i] + 2.0 * ts_zzzz_xxzzzz[i] * fz_0;

        tk_zzzz_xyyyyy[i] =
            -6.0 * ts_zz_xyyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyyyyy[i] * fe_0 + tk_zzz_xyyyyy[i] * pa_z[i] + 2.0 * ts_zzzz_xyyyyy[i] * fz_0;

        tk_zzzz_xyyyyz[i] = -6.0 * ts_zz_xyyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyyyyz[i] * fe_0 + tk_zzz_xyyyy[i] * fe_0 +
                            tk_zzz_xyyyyz[i] * pa_z[i] + 2.0 * ts_zzzz_xyyyyz[i] * fz_0;

        tk_zzzz_xyyyzz[i] = -6.0 * ts_zz_xyyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyyyzz[i] * fe_0 + 2.0 * tk_zzz_xyyyz[i] * fe_0 +
                            tk_zzz_xyyyzz[i] * pa_z[i] + 2.0 * ts_zzzz_xyyyzz[i] * fz_0;

        tk_zzzz_xyyzzz[i] = -6.0 * ts_zz_xyyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyyzzz[i] * fe_0 + 3.0 * tk_zzz_xyyzz[i] * fe_0 +
                            tk_zzz_xyyzzz[i] * pa_z[i] + 2.0 * ts_zzzz_xyyzzz[i] * fz_0;

        tk_zzzz_xyzzzz[i] = -6.0 * ts_zz_xyzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xyzzzz[i] * fe_0 + 4.0 * tk_zzz_xyzzz[i] * fe_0 +
                            tk_zzz_xyzzzz[i] * pa_z[i] + 2.0 * ts_zzzz_xyzzzz[i] * fz_0;

        tk_zzzz_xzzzzz[i] = -6.0 * ts_zz_xzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xzzzzz[i] * fe_0 + 5.0 * tk_zzz_xzzzz[i] * fe_0 +
                            tk_zzz_xzzzzz[i] * pa_z[i] + 2.0 * ts_zzzz_xzzzzz[i] * fz_0;

        tk_zzzz_yyyyyy[i] =
            -6.0 * ts_zz_yyyyyy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyyyyy[i] * fe_0 + tk_zzz_yyyyyy[i] * pa_z[i] + 2.0 * ts_zzzz_yyyyyy[i] * fz_0;

        tk_zzzz_yyyyyz[i] = -6.0 * ts_zz_yyyyyz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyyyyz[i] * fe_0 + tk_zzz_yyyyy[i] * fe_0 +
                            tk_zzz_yyyyyz[i] * pa_z[i] + 2.0 * ts_zzzz_yyyyyz[i] * fz_0;

        tk_zzzz_yyyyzz[i] = -6.0 * ts_zz_yyyyzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyyyzz[i] * fe_0 + 2.0 * tk_zzz_yyyyz[i] * fe_0 +
                            tk_zzz_yyyyzz[i] * pa_z[i] + 2.0 * ts_zzzz_yyyyzz[i] * fz_0;

        tk_zzzz_yyyzzz[i] = -6.0 * ts_zz_yyyzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyyzzz[i] * fe_0 + 3.0 * tk_zzz_yyyzz[i] * fe_0 +
                            tk_zzz_yyyzzz[i] * pa_z[i] + 2.0 * ts_zzzz_yyyzzz[i] * fz_0;

        tk_zzzz_yyzzzz[i] = -6.0 * ts_zz_yyzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yyzzzz[i] * fe_0 + 4.0 * tk_zzz_yyzzz[i] * fe_0 +
                            tk_zzz_yyzzzz[i] * pa_z[i] + 2.0 * ts_zzzz_yyzzzz[i] * fz_0;

        tk_zzzz_yzzzzz[i] = -6.0 * ts_zz_yzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yzzzzz[i] * fe_0 + 5.0 * tk_zzz_yzzzz[i] * fe_0 +
                            tk_zzz_yzzzzz[i] * pa_z[i] + 2.0 * ts_zzzz_yzzzzz[i] * fz_0;

        tk_zzzz_zzzzzz[i] = -6.0 * ts_zz_zzzzzz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_zzzzzz[i] * fe_0 + 6.0 * tk_zzz_zzzzz[i] * fe_0 +
                            tk_zzz_zzzzzz[i] * pa_z[i] + 2.0 * ts_zzzz_zzzzzz[i] * fz_0;
    }
}

}  // namespace kinrec
