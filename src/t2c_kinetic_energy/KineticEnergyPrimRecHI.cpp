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

#include "KineticEnergyPrimRecHI.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_hi(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_hi,
                            const size_t              idx_ovl_fi,
                            const size_t              idx_kin_fi,
                            const size_t              idx_kin_gh,
                            const size_t              idx_kin_gi,
                            const size_t              idx_ovl_hi,
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

    // Set up components of auxiliary buffer : FI

    auto ts_xxx_xxxxxx = pbuffer.data(idx_ovl_fi);

    auto ts_xxx_xxxxxy = pbuffer.data(idx_ovl_fi + 1);

    auto ts_xxx_xxxxxz = pbuffer.data(idx_ovl_fi + 2);

    auto ts_xxx_xxxxyy = pbuffer.data(idx_ovl_fi + 3);

    auto ts_xxx_xxxxyz = pbuffer.data(idx_ovl_fi + 4);

    auto ts_xxx_xxxxzz = pbuffer.data(idx_ovl_fi + 5);

    auto ts_xxx_xxxyyy = pbuffer.data(idx_ovl_fi + 6);

    auto ts_xxx_xxxyyz = pbuffer.data(idx_ovl_fi + 7);

    auto ts_xxx_xxxyzz = pbuffer.data(idx_ovl_fi + 8);

    auto ts_xxx_xxxzzz = pbuffer.data(idx_ovl_fi + 9);

    auto ts_xxx_xxyyyy = pbuffer.data(idx_ovl_fi + 10);

    auto ts_xxx_xxyyyz = pbuffer.data(idx_ovl_fi + 11);

    auto ts_xxx_xxyyzz = pbuffer.data(idx_ovl_fi + 12);

    auto ts_xxx_xxyzzz = pbuffer.data(idx_ovl_fi + 13);

    auto ts_xxx_xxzzzz = pbuffer.data(idx_ovl_fi + 14);

    auto ts_xxx_xyyyyy = pbuffer.data(idx_ovl_fi + 15);

    auto ts_xxx_xyyyyz = pbuffer.data(idx_ovl_fi + 16);

    auto ts_xxx_xyyyzz = pbuffer.data(idx_ovl_fi + 17);

    auto ts_xxx_xyyzzz = pbuffer.data(idx_ovl_fi + 18);

    auto ts_xxx_xyzzzz = pbuffer.data(idx_ovl_fi + 19);

    auto ts_xxx_xzzzzz = pbuffer.data(idx_ovl_fi + 20);

    auto ts_xxx_yyyyyy = pbuffer.data(idx_ovl_fi + 21);

    auto ts_xxx_yyyyyz = pbuffer.data(idx_ovl_fi + 22);

    auto ts_xxx_yyyyzz = pbuffer.data(idx_ovl_fi + 23);

    auto ts_xxx_yyyzzz = pbuffer.data(idx_ovl_fi + 24);

    auto ts_xxx_yyzzzz = pbuffer.data(idx_ovl_fi + 25);

    auto ts_xxx_yzzzzz = pbuffer.data(idx_ovl_fi + 26);

    auto ts_xxx_zzzzzz = pbuffer.data(idx_ovl_fi + 27);

    auto ts_xxy_xxxxxx = pbuffer.data(idx_ovl_fi + 28);

    auto ts_xxy_xxxxxz = pbuffer.data(idx_ovl_fi + 30);

    auto ts_xxy_xxxxzz = pbuffer.data(idx_ovl_fi + 33);

    auto ts_xxy_xxxzzz = pbuffer.data(idx_ovl_fi + 37);

    auto ts_xxy_xxzzzz = pbuffer.data(idx_ovl_fi + 42);

    auto ts_xxy_xzzzzz = pbuffer.data(idx_ovl_fi + 48);

    auto ts_xxz_xxxxxx = pbuffer.data(idx_ovl_fi + 56);

    auto ts_xxz_xxxxxy = pbuffer.data(idx_ovl_fi + 57);

    auto ts_xxz_xxxxyy = pbuffer.data(idx_ovl_fi + 59);

    auto ts_xxz_xxxyyy = pbuffer.data(idx_ovl_fi + 62);

    auto ts_xxz_xxyyyy = pbuffer.data(idx_ovl_fi + 66);

    auto ts_xxz_xyyyyy = pbuffer.data(idx_ovl_fi + 71);

    auto ts_xyy_xxxxxy = pbuffer.data(idx_ovl_fi + 85);

    auto ts_xyy_xxxxyy = pbuffer.data(idx_ovl_fi + 87);

    auto ts_xyy_xxxxyz = pbuffer.data(idx_ovl_fi + 88);

    auto ts_xyy_xxxyyy = pbuffer.data(idx_ovl_fi + 90);

    auto ts_xyy_xxxyyz = pbuffer.data(idx_ovl_fi + 91);

    auto ts_xyy_xxxyzz = pbuffer.data(idx_ovl_fi + 92);

    auto ts_xyy_xxyyyy = pbuffer.data(idx_ovl_fi + 94);

    auto ts_xyy_xxyyyz = pbuffer.data(idx_ovl_fi + 95);

    auto ts_xyy_xxyyzz = pbuffer.data(idx_ovl_fi + 96);

    auto ts_xyy_xxyzzz = pbuffer.data(idx_ovl_fi + 97);

    auto ts_xyy_xyyyyy = pbuffer.data(idx_ovl_fi + 99);

    auto ts_xyy_xyyyyz = pbuffer.data(idx_ovl_fi + 100);

    auto ts_xyy_xyyyzz = pbuffer.data(idx_ovl_fi + 101);

    auto ts_xyy_xyyzzz = pbuffer.data(idx_ovl_fi + 102);

    auto ts_xyy_xyzzzz = pbuffer.data(idx_ovl_fi + 103);

    auto ts_xyy_yyyyyy = pbuffer.data(idx_ovl_fi + 105);

    auto ts_xyy_yyyyyz = pbuffer.data(idx_ovl_fi + 106);

    auto ts_xyy_yyyyzz = pbuffer.data(idx_ovl_fi + 107);

    auto ts_xyy_yyyzzz = pbuffer.data(idx_ovl_fi + 108);

    auto ts_xyy_yyzzzz = pbuffer.data(idx_ovl_fi + 109);

    auto ts_xyy_yzzzzz = pbuffer.data(idx_ovl_fi + 110);

    auto ts_xyy_zzzzzz = pbuffer.data(idx_ovl_fi + 111);

    auto ts_xzz_xxxxxz = pbuffer.data(idx_ovl_fi + 142);

    auto ts_xzz_xxxxyz = pbuffer.data(idx_ovl_fi + 144);

    auto ts_xzz_xxxxzz = pbuffer.data(idx_ovl_fi + 145);

    auto ts_xzz_xxxyyz = pbuffer.data(idx_ovl_fi + 147);

    auto ts_xzz_xxxyzz = pbuffer.data(idx_ovl_fi + 148);

    auto ts_xzz_xxxzzz = pbuffer.data(idx_ovl_fi + 149);

    auto ts_xzz_xxyyyz = pbuffer.data(idx_ovl_fi + 151);

    auto ts_xzz_xxyyzz = pbuffer.data(idx_ovl_fi + 152);

    auto ts_xzz_xxyzzz = pbuffer.data(idx_ovl_fi + 153);

    auto ts_xzz_xxzzzz = pbuffer.data(idx_ovl_fi + 154);

    auto ts_xzz_xyyyyz = pbuffer.data(idx_ovl_fi + 156);

    auto ts_xzz_xyyyzz = pbuffer.data(idx_ovl_fi + 157);

    auto ts_xzz_xyyzzz = pbuffer.data(idx_ovl_fi + 158);

    auto ts_xzz_xyzzzz = pbuffer.data(idx_ovl_fi + 159);

    auto ts_xzz_xzzzzz = pbuffer.data(idx_ovl_fi + 160);

    auto ts_xzz_yyyyyy = pbuffer.data(idx_ovl_fi + 161);

    auto ts_xzz_yyyyyz = pbuffer.data(idx_ovl_fi + 162);

    auto ts_xzz_yyyyzz = pbuffer.data(idx_ovl_fi + 163);

    auto ts_xzz_yyyzzz = pbuffer.data(idx_ovl_fi + 164);

    auto ts_xzz_yyzzzz = pbuffer.data(idx_ovl_fi + 165);

    auto ts_xzz_yzzzzz = pbuffer.data(idx_ovl_fi + 166);

    auto ts_xzz_zzzzzz = pbuffer.data(idx_ovl_fi + 167);

    auto ts_yyy_xxxxxx = pbuffer.data(idx_ovl_fi + 168);

    auto ts_yyy_xxxxxy = pbuffer.data(idx_ovl_fi + 169);

    auto ts_yyy_xxxxxz = pbuffer.data(idx_ovl_fi + 170);

    auto ts_yyy_xxxxyy = pbuffer.data(idx_ovl_fi + 171);

    auto ts_yyy_xxxxyz = pbuffer.data(idx_ovl_fi + 172);

    auto ts_yyy_xxxxzz = pbuffer.data(idx_ovl_fi + 173);

    auto ts_yyy_xxxyyy = pbuffer.data(idx_ovl_fi + 174);

    auto ts_yyy_xxxyyz = pbuffer.data(idx_ovl_fi + 175);

    auto ts_yyy_xxxyzz = pbuffer.data(idx_ovl_fi + 176);

    auto ts_yyy_xxxzzz = pbuffer.data(idx_ovl_fi + 177);

    auto ts_yyy_xxyyyy = pbuffer.data(idx_ovl_fi + 178);

    auto ts_yyy_xxyyyz = pbuffer.data(idx_ovl_fi + 179);

    auto ts_yyy_xxyyzz = pbuffer.data(idx_ovl_fi + 180);

    auto ts_yyy_xxyzzz = pbuffer.data(idx_ovl_fi + 181);

    auto ts_yyy_xxzzzz = pbuffer.data(idx_ovl_fi + 182);

    auto ts_yyy_xyyyyy = pbuffer.data(idx_ovl_fi + 183);

    auto ts_yyy_xyyyyz = pbuffer.data(idx_ovl_fi + 184);

    auto ts_yyy_xyyyzz = pbuffer.data(idx_ovl_fi + 185);

    auto ts_yyy_xyyzzz = pbuffer.data(idx_ovl_fi + 186);

    auto ts_yyy_xyzzzz = pbuffer.data(idx_ovl_fi + 187);

    auto ts_yyy_xzzzzz = pbuffer.data(idx_ovl_fi + 188);

    auto ts_yyy_yyyyyy = pbuffer.data(idx_ovl_fi + 189);

    auto ts_yyy_yyyyyz = pbuffer.data(idx_ovl_fi + 190);

    auto ts_yyy_yyyyzz = pbuffer.data(idx_ovl_fi + 191);

    auto ts_yyy_yyyzzz = pbuffer.data(idx_ovl_fi + 192);

    auto ts_yyy_yyzzzz = pbuffer.data(idx_ovl_fi + 193);

    auto ts_yyy_yzzzzz = pbuffer.data(idx_ovl_fi + 194);

    auto ts_yyy_zzzzzz = pbuffer.data(idx_ovl_fi + 195);

    auto ts_yyz_xxxxxy = pbuffer.data(idx_ovl_fi + 197);

    auto ts_yyz_xxxxyy = pbuffer.data(idx_ovl_fi + 199);

    auto ts_yyz_xxxyyy = pbuffer.data(idx_ovl_fi + 202);

    auto ts_yyz_xxyyyy = pbuffer.data(idx_ovl_fi + 206);

    auto ts_yyz_xyyyyy = pbuffer.data(idx_ovl_fi + 211);

    auto ts_yyz_yyyyyy = pbuffer.data(idx_ovl_fi + 217);

    auto ts_yzz_xxxxxx = pbuffer.data(idx_ovl_fi + 224);

    auto ts_yzz_xxxxxz = pbuffer.data(idx_ovl_fi + 226);

    auto ts_yzz_xxxxyz = pbuffer.data(idx_ovl_fi + 228);

    auto ts_yzz_xxxxzz = pbuffer.data(idx_ovl_fi + 229);

    auto ts_yzz_xxxyyz = pbuffer.data(idx_ovl_fi + 231);

    auto ts_yzz_xxxyzz = pbuffer.data(idx_ovl_fi + 232);

    auto ts_yzz_xxxzzz = pbuffer.data(idx_ovl_fi + 233);

    auto ts_yzz_xxyyyz = pbuffer.data(idx_ovl_fi + 235);

    auto ts_yzz_xxyyzz = pbuffer.data(idx_ovl_fi + 236);

    auto ts_yzz_xxyzzz = pbuffer.data(idx_ovl_fi + 237);

    auto ts_yzz_xxzzzz = pbuffer.data(idx_ovl_fi + 238);

    auto ts_yzz_xyyyyz = pbuffer.data(idx_ovl_fi + 240);

    auto ts_yzz_xyyyzz = pbuffer.data(idx_ovl_fi + 241);

    auto ts_yzz_xyyzzz = pbuffer.data(idx_ovl_fi + 242);

    auto ts_yzz_xyzzzz = pbuffer.data(idx_ovl_fi + 243);

    auto ts_yzz_xzzzzz = pbuffer.data(idx_ovl_fi + 244);

    auto ts_yzz_yyyyyz = pbuffer.data(idx_ovl_fi + 246);

    auto ts_yzz_yyyyzz = pbuffer.data(idx_ovl_fi + 247);

    auto ts_yzz_yyyzzz = pbuffer.data(idx_ovl_fi + 248);

    auto ts_yzz_yyzzzz = pbuffer.data(idx_ovl_fi + 249);

    auto ts_yzz_yzzzzz = pbuffer.data(idx_ovl_fi + 250);

    auto ts_yzz_zzzzzz = pbuffer.data(idx_ovl_fi + 251);

    auto ts_zzz_xxxxxx = pbuffer.data(idx_ovl_fi + 252);

    auto ts_zzz_xxxxxy = pbuffer.data(idx_ovl_fi + 253);

    auto ts_zzz_xxxxxz = pbuffer.data(idx_ovl_fi + 254);

    auto ts_zzz_xxxxyy = pbuffer.data(idx_ovl_fi + 255);

    auto ts_zzz_xxxxyz = pbuffer.data(idx_ovl_fi + 256);

    auto ts_zzz_xxxxzz = pbuffer.data(idx_ovl_fi + 257);

    auto ts_zzz_xxxyyy = pbuffer.data(idx_ovl_fi + 258);

    auto ts_zzz_xxxyyz = pbuffer.data(idx_ovl_fi + 259);

    auto ts_zzz_xxxyzz = pbuffer.data(idx_ovl_fi + 260);

    auto ts_zzz_xxxzzz = pbuffer.data(idx_ovl_fi + 261);

    auto ts_zzz_xxyyyy = pbuffer.data(idx_ovl_fi + 262);

    auto ts_zzz_xxyyyz = pbuffer.data(idx_ovl_fi + 263);

    auto ts_zzz_xxyyzz = pbuffer.data(idx_ovl_fi + 264);

    auto ts_zzz_xxyzzz = pbuffer.data(idx_ovl_fi + 265);

    auto ts_zzz_xxzzzz = pbuffer.data(idx_ovl_fi + 266);

    auto ts_zzz_xyyyyy = pbuffer.data(idx_ovl_fi + 267);

    auto ts_zzz_xyyyyz = pbuffer.data(idx_ovl_fi + 268);

    auto ts_zzz_xyyyzz = pbuffer.data(idx_ovl_fi + 269);

    auto ts_zzz_xyyzzz = pbuffer.data(idx_ovl_fi + 270);

    auto ts_zzz_xyzzzz = pbuffer.data(idx_ovl_fi + 271);

    auto ts_zzz_xzzzzz = pbuffer.data(idx_ovl_fi + 272);

    auto ts_zzz_yyyyyy = pbuffer.data(idx_ovl_fi + 273);

    auto ts_zzz_yyyyyz = pbuffer.data(idx_ovl_fi + 274);

    auto ts_zzz_yyyyzz = pbuffer.data(idx_ovl_fi + 275);

    auto ts_zzz_yyyzzz = pbuffer.data(idx_ovl_fi + 276);

    auto ts_zzz_yyzzzz = pbuffer.data(idx_ovl_fi + 277);

    auto ts_zzz_yzzzzz = pbuffer.data(idx_ovl_fi + 278);

    auto ts_zzz_zzzzzz = pbuffer.data(idx_ovl_fi + 279);

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

    auto tk_xxy_xxxxxz = pbuffer.data(idx_kin_fi + 30);

    auto tk_xxy_xxxxzz = pbuffer.data(idx_kin_fi + 33);

    auto tk_xxy_xxxzzz = pbuffer.data(idx_kin_fi + 37);

    auto tk_xxy_xxzzzz = pbuffer.data(idx_kin_fi + 42);

    auto tk_xxy_xzzzzz = pbuffer.data(idx_kin_fi + 48);

    auto tk_xxz_xxxxxx = pbuffer.data(idx_kin_fi + 56);

    auto tk_xxz_xxxxxy = pbuffer.data(idx_kin_fi + 57);

    auto tk_xxz_xxxxyy = pbuffer.data(idx_kin_fi + 59);

    auto tk_xxz_xxxyyy = pbuffer.data(idx_kin_fi + 62);

    auto tk_xxz_xxyyyy = pbuffer.data(idx_kin_fi + 66);

    auto tk_xxz_xyyyyy = pbuffer.data(idx_kin_fi + 71);

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

    auto tk_yyz_xxxxyy = pbuffer.data(idx_kin_fi + 199);

    auto tk_yyz_xxxyyy = pbuffer.data(idx_kin_fi + 202);

    auto tk_yyz_xxyyyy = pbuffer.data(idx_kin_fi + 206);

    auto tk_yyz_xyyyyy = pbuffer.data(idx_kin_fi + 211);

    auto tk_yyz_yyyyyy = pbuffer.data(idx_kin_fi + 217);

    auto tk_yzz_xxxxxx = pbuffer.data(idx_kin_fi + 224);

    auto tk_yzz_xxxxxz = pbuffer.data(idx_kin_fi + 226);

    auto tk_yzz_xxxxyz = pbuffer.data(idx_kin_fi + 228);

    auto tk_yzz_xxxxzz = pbuffer.data(idx_kin_fi + 229);

    auto tk_yzz_xxxyyz = pbuffer.data(idx_kin_fi + 231);

    auto tk_yzz_xxxyzz = pbuffer.data(idx_kin_fi + 232);

    auto tk_yzz_xxxzzz = pbuffer.data(idx_kin_fi + 233);

    auto tk_yzz_xxyyyz = pbuffer.data(idx_kin_fi + 235);

    auto tk_yzz_xxyyzz = pbuffer.data(idx_kin_fi + 236);

    auto tk_yzz_xxyzzz = pbuffer.data(idx_kin_fi + 237);

    auto tk_yzz_xxzzzz = pbuffer.data(idx_kin_fi + 238);

    auto tk_yzz_xyyyyz = pbuffer.data(idx_kin_fi + 240);

    auto tk_yzz_xyyyzz = pbuffer.data(idx_kin_fi + 241);

    auto tk_yzz_xyyzzz = pbuffer.data(idx_kin_fi + 242);

    auto tk_yzz_xyzzzz = pbuffer.data(idx_kin_fi + 243);

    auto tk_yzz_xzzzzz = pbuffer.data(idx_kin_fi + 244);

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

    // Set up components of auxiliary buffer : GH

    auto tk_xxxx_xxxxx = pbuffer.data(idx_kin_gh);

    auto tk_xxxx_xxxxy = pbuffer.data(idx_kin_gh + 1);

    auto tk_xxxx_xxxxz = pbuffer.data(idx_kin_gh + 2);

    auto tk_xxxx_xxxyy = pbuffer.data(idx_kin_gh + 3);

    auto tk_xxxx_xxxyz = pbuffer.data(idx_kin_gh + 4);

    auto tk_xxxx_xxxzz = pbuffer.data(idx_kin_gh + 5);

    auto tk_xxxx_xxyyy = pbuffer.data(idx_kin_gh + 6);

    auto tk_xxxx_xxyyz = pbuffer.data(idx_kin_gh + 7);

    auto tk_xxxx_xxyzz = pbuffer.data(idx_kin_gh + 8);

    auto tk_xxxx_xxzzz = pbuffer.data(idx_kin_gh + 9);

    auto tk_xxxx_xyyyy = pbuffer.data(idx_kin_gh + 10);

    auto tk_xxxx_xyyyz = pbuffer.data(idx_kin_gh + 11);

    auto tk_xxxx_xyyzz = pbuffer.data(idx_kin_gh + 12);

    auto tk_xxxx_xyzzz = pbuffer.data(idx_kin_gh + 13);

    auto tk_xxxx_xzzzz = pbuffer.data(idx_kin_gh + 14);

    auto tk_xxxx_yyyyy = pbuffer.data(idx_kin_gh + 15);

    auto tk_xxxx_yyyyz = pbuffer.data(idx_kin_gh + 16);

    auto tk_xxxx_yyyzz = pbuffer.data(idx_kin_gh + 17);

    auto tk_xxxx_yyzzz = pbuffer.data(idx_kin_gh + 18);

    auto tk_xxxx_yzzzz = pbuffer.data(idx_kin_gh + 19);

    auto tk_xxxx_zzzzz = pbuffer.data(idx_kin_gh + 20);

    auto tk_xxxz_xxxxz = pbuffer.data(idx_kin_gh + 44);

    auto tk_xxxz_xxxyz = pbuffer.data(idx_kin_gh + 46);

    auto tk_xxxz_xxxzz = pbuffer.data(idx_kin_gh + 47);

    auto tk_xxxz_xxyyz = pbuffer.data(idx_kin_gh + 49);

    auto tk_xxxz_xxyzz = pbuffer.data(idx_kin_gh + 50);

    auto tk_xxxz_xxzzz = pbuffer.data(idx_kin_gh + 51);

    auto tk_xxxz_xyyyz = pbuffer.data(idx_kin_gh + 53);

    auto tk_xxxz_xyyzz = pbuffer.data(idx_kin_gh + 54);

    auto tk_xxxz_xyzzz = pbuffer.data(idx_kin_gh + 55);

    auto tk_xxxz_xzzzz = pbuffer.data(idx_kin_gh + 56);

    auto tk_xxxz_yyyyz = pbuffer.data(idx_kin_gh + 58);

    auto tk_xxxz_yyyzz = pbuffer.data(idx_kin_gh + 59);

    auto tk_xxxz_yyzzz = pbuffer.data(idx_kin_gh + 60);

    auto tk_xxxz_yzzzz = pbuffer.data(idx_kin_gh + 61);

    auto tk_xxxz_zzzzz = pbuffer.data(idx_kin_gh + 62);

    auto tk_xxyy_xxxxx = pbuffer.data(idx_kin_gh + 63);

    auto tk_xxyy_xxxxy = pbuffer.data(idx_kin_gh + 64);

    auto tk_xxyy_xxxxz = pbuffer.data(idx_kin_gh + 65);

    auto tk_xxyy_xxxyy = pbuffer.data(idx_kin_gh + 66);

    auto tk_xxyy_xxxyz = pbuffer.data(idx_kin_gh + 67);

    auto tk_xxyy_xxxzz = pbuffer.data(idx_kin_gh + 68);

    auto tk_xxyy_xxyyy = pbuffer.data(idx_kin_gh + 69);

    auto tk_xxyy_xxyyz = pbuffer.data(idx_kin_gh + 70);

    auto tk_xxyy_xxyzz = pbuffer.data(idx_kin_gh + 71);

    auto tk_xxyy_xxzzz = pbuffer.data(idx_kin_gh + 72);

    auto tk_xxyy_xyyyy = pbuffer.data(idx_kin_gh + 73);

    auto tk_xxyy_xyyyz = pbuffer.data(idx_kin_gh + 74);

    auto tk_xxyy_xyyzz = pbuffer.data(idx_kin_gh + 75);

    auto tk_xxyy_xyzzz = pbuffer.data(idx_kin_gh + 76);

    auto tk_xxyy_xzzzz = pbuffer.data(idx_kin_gh + 77);

    auto tk_xxyy_yyyyy = pbuffer.data(idx_kin_gh + 78);

    auto tk_xxyy_yyyyz = pbuffer.data(idx_kin_gh + 79);

    auto tk_xxyy_yyyzz = pbuffer.data(idx_kin_gh + 80);

    auto tk_xxyy_yyzzz = pbuffer.data(idx_kin_gh + 81);

    auto tk_xxyy_yzzzz = pbuffer.data(idx_kin_gh + 82);

    auto tk_xxyy_zzzzz = pbuffer.data(idx_kin_gh + 83);

    auto tk_xxzz_xxxxx = pbuffer.data(idx_kin_gh + 105);

    auto tk_xxzz_xxxxy = pbuffer.data(idx_kin_gh + 106);

    auto tk_xxzz_xxxxz = pbuffer.data(idx_kin_gh + 107);

    auto tk_xxzz_xxxyy = pbuffer.data(idx_kin_gh + 108);

    auto tk_xxzz_xxxyz = pbuffer.data(idx_kin_gh + 109);

    auto tk_xxzz_xxxzz = pbuffer.data(idx_kin_gh + 110);

    auto tk_xxzz_xxyyy = pbuffer.data(idx_kin_gh + 111);

    auto tk_xxzz_xxyyz = pbuffer.data(idx_kin_gh + 112);

    auto tk_xxzz_xxyzz = pbuffer.data(idx_kin_gh + 113);

    auto tk_xxzz_xxzzz = pbuffer.data(idx_kin_gh + 114);

    auto tk_xxzz_xyyyy = pbuffer.data(idx_kin_gh + 115);

    auto tk_xxzz_xyyyz = pbuffer.data(idx_kin_gh + 116);

    auto tk_xxzz_xyyzz = pbuffer.data(idx_kin_gh + 117);

    auto tk_xxzz_xyzzz = pbuffer.data(idx_kin_gh + 118);

    auto tk_xxzz_xzzzz = pbuffer.data(idx_kin_gh + 119);

    auto tk_xxzz_yyyyy = pbuffer.data(idx_kin_gh + 120);

    auto tk_xxzz_yyyyz = pbuffer.data(idx_kin_gh + 121);

    auto tk_xxzz_yyyzz = pbuffer.data(idx_kin_gh + 122);

    auto tk_xxzz_yyzzz = pbuffer.data(idx_kin_gh + 123);

    auto tk_xxzz_yzzzz = pbuffer.data(idx_kin_gh + 124);

    auto tk_xxzz_zzzzz = pbuffer.data(idx_kin_gh + 125);

    auto tk_xyyy_xxxxy = pbuffer.data(idx_kin_gh + 127);

    auto tk_xyyy_xxxyy = pbuffer.data(idx_kin_gh + 129);

    auto tk_xyyy_xxxyz = pbuffer.data(idx_kin_gh + 130);

    auto tk_xyyy_xxyyy = pbuffer.data(idx_kin_gh + 132);

    auto tk_xyyy_xxyyz = pbuffer.data(idx_kin_gh + 133);

    auto tk_xyyy_xxyzz = pbuffer.data(idx_kin_gh + 134);

    auto tk_xyyy_xyyyy = pbuffer.data(idx_kin_gh + 136);

    auto tk_xyyy_xyyyz = pbuffer.data(idx_kin_gh + 137);

    auto tk_xyyy_xyyzz = pbuffer.data(idx_kin_gh + 138);

    auto tk_xyyy_xyzzz = pbuffer.data(idx_kin_gh + 139);

    auto tk_xyyy_yyyyy = pbuffer.data(idx_kin_gh + 141);

    auto tk_xyyy_yyyyz = pbuffer.data(idx_kin_gh + 142);

    auto tk_xyyy_yyyzz = pbuffer.data(idx_kin_gh + 143);

    auto tk_xyyy_yyzzz = pbuffer.data(idx_kin_gh + 144);

    auto tk_xyyy_yzzzz = pbuffer.data(idx_kin_gh + 145);

    auto tk_xzzz_xxxxz = pbuffer.data(idx_kin_gh + 191);

    auto tk_xzzz_xxxyz = pbuffer.data(idx_kin_gh + 193);

    auto tk_xzzz_xxxzz = pbuffer.data(idx_kin_gh + 194);

    auto tk_xzzz_xxyyz = pbuffer.data(idx_kin_gh + 196);

    auto tk_xzzz_xxyzz = pbuffer.data(idx_kin_gh + 197);

    auto tk_xzzz_xxzzz = pbuffer.data(idx_kin_gh + 198);

    auto tk_xzzz_xyyyz = pbuffer.data(idx_kin_gh + 200);

    auto tk_xzzz_xyyzz = pbuffer.data(idx_kin_gh + 201);

    auto tk_xzzz_xyzzz = pbuffer.data(idx_kin_gh + 202);

    auto tk_xzzz_xzzzz = pbuffer.data(idx_kin_gh + 203);

    auto tk_xzzz_yyyyz = pbuffer.data(idx_kin_gh + 205);

    auto tk_xzzz_yyyzz = pbuffer.data(idx_kin_gh + 206);

    auto tk_xzzz_yyzzz = pbuffer.data(idx_kin_gh + 207);

    auto tk_xzzz_yzzzz = pbuffer.data(idx_kin_gh + 208);

    auto tk_xzzz_zzzzz = pbuffer.data(idx_kin_gh + 209);

    auto tk_yyyy_xxxxx = pbuffer.data(idx_kin_gh + 210);

    auto tk_yyyy_xxxxy = pbuffer.data(idx_kin_gh + 211);

    auto tk_yyyy_xxxxz = pbuffer.data(idx_kin_gh + 212);

    auto tk_yyyy_xxxyy = pbuffer.data(idx_kin_gh + 213);

    auto tk_yyyy_xxxyz = pbuffer.data(idx_kin_gh + 214);

    auto tk_yyyy_xxxzz = pbuffer.data(idx_kin_gh + 215);

    auto tk_yyyy_xxyyy = pbuffer.data(idx_kin_gh + 216);

    auto tk_yyyy_xxyyz = pbuffer.data(idx_kin_gh + 217);

    auto tk_yyyy_xxyzz = pbuffer.data(idx_kin_gh + 218);

    auto tk_yyyy_xxzzz = pbuffer.data(idx_kin_gh + 219);

    auto tk_yyyy_xyyyy = pbuffer.data(idx_kin_gh + 220);

    auto tk_yyyy_xyyyz = pbuffer.data(idx_kin_gh + 221);

    auto tk_yyyy_xyyzz = pbuffer.data(idx_kin_gh + 222);

    auto tk_yyyy_xyzzz = pbuffer.data(idx_kin_gh + 223);

    auto tk_yyyy_xzzzz = pbuffer.data(idx_kin_gh + 224);

    auto tk_yyyy_yyyyy = pbuffer.data(idx_kin_gh + 225);

    auto tk_yyyy_yyyyz = pbuffer.data(idx_kin_gh + 226);

    auto tk_yyyy_yyyzz = pbuffer.data(idx_kin_gh + 227);

    auto tk_yyyy_yyzzz = pbuffer.data(idx_kin_gh + 228);

    auto tk_yyyy_yzzzz = pbuffer.data(idx_kin_gh + 229);

    auto tk_yyyy_zzzzz = pbuffer.data(idx_kin_gh + 230);

    auto tk_yyyz_xxxxz = pbuffer.data(idx_kin_gh + 233);

    auto tk_yyyz_xxxyz = pbuffer.data(idx_kin_gh + 235);

    auto tk_yyyz_xxxzz = pbuffer.data(idx_kin_gh + 236);

    auto tk_yyyz_xxyyz = pbuffer.data(idx_kin_gh + 238);

    auto tk_yyyz_xxyzz = pbuffer.data(idx_kin_gh + 239);

    auto tk_yyyz_xxzzz = pbuffer.data(idx_kin_gh + 240);

    auto tk_yyyz_xyyyz = pbuffer.data(idx_kin_gh + 242);

    auto tk_yyyz_xyyzz = pbuffer.data(idx_kin_gh + 243);

    auto tk_yyyz_xyzzz = pbuffer.data(idx_kin_gh + 244);

    auto tk_yyyz_xzzzz = pbuffer.data(idx_kin_gh + 245);

    auto tk_yyyz_yyyyz = pbuffer.data(idx_kin_gh + 247);

    auto tk_yyyz_yyyzz = pbuffer.data(idx_kin_gh + 248);

    auto tk_yyyz_yyzzz = pbuffer.data(idx_kin_gh + 249);

    auto tk_yyyz_yzzzz = pbuffer.data(idx_kin_gh + 250);

    auto tk_yyyz_zzzzz = pbuffer.data(idx_kin_gh + 251);

    auto tk_yyzz_xxxxx = pbuffer.data(idx_kin_gh + 252);

    auto tk_yyzz_xxxxy = pbuffer.data(idx_kin_gh + 253);

    auto tk_yyzz_xxxxz = pbuffer.data(idx_kin_gh + 254);

    auto tk_yyzz_xxxyy = pbuffer.data(idx_kin_gh + 255);

    auto tk_yyzz_xxxyz = pbuffer.data(idx_kin_gh + 256);

    auto tk_yyzz_xxxzz = pbuffer.data(idx_kin_gh + 257);

    auto tk_yyzz_xxyyy = pbuffer.data(idx_kin_gh + 258);

    auto tk_yyzz_xxyyz = pbuffer.data(idx_kin_gh + 259);

    auto tk_yyzz_xxyzz = pbuffer.data(idx_kin_gh + 260);

    auto tk_yyzz_xxzzz = pbuffer.data(idx_kin_gh + 261);

    auto tk_yyzz_xyyyy = pbuffer.data(idx_kin_gh + 262);

    auto tk_yyzz_xyyyz = pbuffer.data(idx_kin_gh + 263);

    auto tk_yyzz_xyyzz = pbuffer.data(idx_kin_gh + 264);

    auto tk_yyzz_xyzzz = pbuffer.data(idx_kin_gh + 265);

    auto tk_yyzz_xzzzz = pbuffer.data(idx_kin_gh + 266);

    auto tk_yyzz_yyyyy = pbuffer.data(idx_kin_gh + 267);

    auto tk_yyzz_yyyyz = pbuffer.data(idx_kin_gh + 268);

    auto tk_yyzz_yyyzz = pbuffer.data(idx_kin_gh + 269);

    auto tk_yyzz_yyzzz = pbuffer.data(idx_kin_gh + 270);

    auto tk_yyzz_yzzzz = pbuffer.data(idx_kin_gh + 271);

    auto tk_yyzz_zzzzz = pbuffer.data(idx_kin_gh + 272);

    auto tk_yzzz_xxxxy = pbuffer.data(idx_kin_gh + 274);

    auto tk_yzzz_xxxxz = pbuffer.data(idx_kin_gh + 275);

    auto tk_yzzz_xxxyy = pbuffer.data(idx_kin_gh + 276);

    auto tk_yzzz_xxxyz = pbuffer.data(idx_kin_gh + 277);

    auto tk_yzzz_xxxzz = pbuffer.data(idx_kin_gh + 278);

    auto tk_yzzz_xxyyy = pbuffer.data(idx_kin_gh + 279);

    auto tk_yzzz_xxyyz = pbuffer.data(idx_kin_gh + 280);

    auto tk_yzzz_xxyzz = pbuffer.data(idx_kin_gh + 281);

    auto tk_yzzz_xxzzz = pbuffer.data(idx_kin_gh + 282);

    auto tk_yzzz_xyyyy = pbuffer.data(idx_kin_gh + 283);

    auto tk_yzzz_xyyyz = pbuffer.data(idx_kin_gh + 284);

    auto tk_yzzz_xyyzz = pbuffer.data(idx_kin_gh + 285);

    auto tk_yzzz_xyzzz = pbuffer.data(idx_kin_gh + 286);

    auto tk_yzzz_xzzzz = pbuffer.data(idx_kin_gh + 287);

    auto tk_yzzz_yyyyy = pbuffer.data(idx_kin_gh + 288);

    auto tk_yzzz_yyyyz = pbuffer.data(idx_kin_gh + 289);

    auto tk_yzzz_yyyzz = pbuffer.data(idx_kin_gh + 290);

    auto tk_yzzz_yyzzz = pbuffer.data(idx_kin_gh + 291);

    auto tk_yzzz_yzzzz = pbuffer.data(idx_kin_gh + 292);

    auto tk_yzzz_zzzzz = pbuffer.data(idx_kin_gh + 293);

    auto tk_zzzz_xxxxx = pbuffer.data(idx_kin_gh + 294);

    auto tk_zzzz_xxxxy = pbuffer.data(idx_kin_gh + 295);

    auto tk_zzzz_xxxxz = pbuffer.data(idx_kin_gh + 296);

    auto tk_zzzz_xxxyy = pbuffer.data(idx_kin_gh + 297);

    auto tk_zzzz_xxxyz = pbuffer.data(idx_kin_gh + 298);

    auto tk_zzzz_xxxzz = pbuffer.data(idx_kin_gh + 299);

    auto tk_zzzz_xxyyy = pbuffer.data(idx_kin_gh + 300);

    auto tk_zzzz_xxyyz = pbuffer.data(idx_kin_gh + 301);

    auto tk_zzzz_xxyzz = pbuffer.data(idx_kin_gh + 302);

    auto tk_zzzz_xxzzz = pbuffer.data(idx_kin_gh + 303);

    auto tk_zzzz_xyyyy = pbuffer.data(idx_kin_gh + 304);

    auto tk_zzzz_xyyyz = pbuffer.data(idx_kin_gh + 305);

    auto tk_zzzz_xyyzz = pbuffer.data(idx_kin_gh + 306);

    auto tk_zzzz_xyzzz = pbuffer.data(idx_kin_gh + 307);

    auto tk_zzzz_xzzzz = pbuffer.data(idx_kin_gh + 308);

    auto tk_zzzz_yyyyy = pbuffer.data(idx_kin_gh + 309);

    auto tk_zzzz_yyyyz = pbuffer.data(idx_kin_gh + 310);

    auto tk_zzzz_yyyzz = pbuffer.data(idx_kin_gh + 311);

    auto tk_zzzz_yyzzz = pbuffer.data(idx_kin_gh + 312);

    auto tk_zzzz_yzzzz = pbuffer.data(idx_kin_gh + 313);

    auto tk_zzzz_zzzzz = pbuffer.data(idx_kin_gh + 314);

    // Set up components of auxiliary buffer : GI

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

    auto tk_xxxy_xxxxxx = pbuffer.data(idx_kin_gi + 28);

    auto tk_xxxy_xxxxxy = pbuffer.data(idx_kin_gi + 29);

    auto tk_xxxy_xxxxxz = pbuffer.data(idx_kin_gi + 30);

    auto tk_xxxy_xxxxyy = pbuffer.data(idx_kin_gi + 31);

    auto tk_xxxy_xxxxzz = pbuffer.data(idx_kin_gi + 33);

    auto tk_xxxy_xxxyyy = pbuffer.data(idx_kin_gi + 34);

    auto tk_xxxy_xxxzzz = pbuffer.data(idx_kin_gi + 37);

    auto tk_xxxy_xxyyyy = pbuffer.data(idx_kin_gi + 38);

    auto tk_xxxy_xxzzzz = pbuffer.data(idx_kin_gi + 42);

    auto tk_xxxy_xyyyyy = pbuffer.data(idx_kin_gi + 43);

    auto tk_xxxy_xzzzzz = pbuffer.data(idx_kin_gi + 48);

    auto tk_xxxy_yyyyyy = pbuffer.data(idx_kin_gi + 49);

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

    auto tk_xxxz_yyyyyz = pbuffer.data(idx_kin_gi + 78);

    auto tk_xxxz_yyyyzz = pbuffer.data(idx_kin_gi + 79);

    auto tk_xxxz_yyyzzz = pbuffer.data(idx_kin_gi + 80);

    auto tk_xxxz_yyzzzz = pbuffer.data(idx_kin_gi + 81);

    auto tk_xxxz_yzzzzz = pbuffer.data(idx_kin_gi + 82);

    auto tk_xxxz_zzzzzz = pbuffer.data(idx_kin_gi + 83);

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

    auto tk_xyyy_xxxxxx = pbuffer.data(idx_kin_gi + 168);

    auto tk_xyyy_xxxxxy = pbuffer.data(idx_kin_gi + 169);

    auto tk_xyyy_xxxxyy = pbuffer.data(idx_kin_gi + 171);

    auto tk_xyyy_xxxxyz = pbuffer.data(idx_kin_gi + 172);

    auto tk_xyyy_xxxyyy = pbuffer.data(idx_kin_gi + 174);

    auto tk_xyyy_xxxyyz = pbuffer.data(idx_kin_gi + 175);

    auto tk_xyyy_xxxyzz = pbuffer.data(idx_kin_gi + 176);

    auto tk_xyyy_xxyyyy = pbuffer.data(idx_kin_gi + 178);

    auto tk_xyyy_xxyyyz = pbuffer.data(idx_kin_gi + 179);

    auto tk_xyyy_xxyyzz = pbuffer.data(idx_kin_gi + 180);

    auto tk_xyyy_xxyzzz = pbuffer.data(idx_kin_gi + 181);

    auto tk_xyyy_xyyyyy = pbuffer.data(idx_kin_gi + 183);

    auto tk_xyyy_xyyyyz = pbuffer.data(idx_kin_gi + 184);

    auto tk_xyyy_xyyyzz = pbuffer.data(idx_kin_gi + 185);

    auto tk_xyyy_xyyzzz = pbuffer.data(idx_kin_gi + 186);

    auto tk_xyyy_xyzzzz = pbuffer.data(idx_kin_gi + 187);

    auto tk_xyyy_yyyyyy = pbuffer.data(idx_kin_gi + 189);

    auto tk_xyyy_yyyyyz = pbuffer.data(idx_kin_gi + 190);

    auto tk_xyyy_yyyyzz = pbuffer.data(idx_kin_gi + 191);

    auto tk_xyyy_yyyzzz = pbuffer.data(idx_kin_gi + 192);

    auto tk_xyyy_yyzzzz = pbuffer.data(idx_kin_gi + 193);

    auto tk_xyyy_yzzzzz = pbuffer.data(idx_kin_gi + 194);

    auto tk_xyyy_zzzzzz = pbuffer.data(idx_kin_gi + 195);

    auto tk_xzzz_xxxxxx = pbuffer.data(idx_kin_gi + 252);

    auto tk_xzzz_xxxxxz = pbuffer.data(idx_kin_gi + 254);

    auto tk_xzzz_xxxxyz = pbuffer.data(idx_kin_gi + 256);

    auto tk_xzzz_xxxxzz = pbuffer.data(idx_kin_gi + 257);

    auto tk_xzzz_xxxyyz = pbuffer.data(idx_kin_gi + 259);

    auto tk_xzzz_xxxyzz = pbuffer.data(idx_kin_gi + 260);

    auto tk_xzzz_xxxzzz = pbuffer.data(idx_kin_gi + 261);

    auto tk_xzzz_xxyyyz = pbuffer.data(idx_kin_gi + 263);

    auto tk_xzzz_xxyyzz = pbuffer.data(idx_kin_gi + 264);

    auto tk_xzzz_xxyzzz = pbuffer.data(idx_kin_gi + 265);

    auto tk_xzzz_xxzzzz = pbuffer.data(idx_kin_gi + 266);

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

    // Set up components of auxiliary buffer : HI

    auto ts_xxxxx_xxxxxx = pbuffer.data(idx_ovl_hi);

    auto ts_xxxxx_xxxxxy = pbuffer.data(idx_ovl_hi + 1);

    auto ts_xxxxx_xxxxxz = pbuffer.data(idx_ovl_hi + 2);

    auto ts_xxxxx_xxxxyy = pbuffer.data(idx_ovl_hi + 3);

    auto ts_xxxxx_xxxxyz = pbuffer.data(idx_ovl_hi + 4);

    auto ts_xxxxx_xxxxzz = pbuffer.data(idx_ovl_hi + 5);

    auto ts_xxxxx_xxxyyy = pbuffer.data(idx_ovl_hi + 6);

    auto ts_xxxxx_xxxyyz = pbuffer.data(idx_ovl_hi + 7);

    auto ts_xxxxx_xxxyzz = pbuffer.data(idx_ovl_hi + 8);

    auto ts_xxxxx_xxxzzz = pbuffer.data(idx_ovl_hi + 9);

    auto ts_xxxxx_xxyyyy = pbuffer.data(idx_ovl_hi + 10);

    auto ts_xxxxx_xxyyyz = pbuffer.data(idx_ovl_hi + 11);

    auto ts_xxxxx_xxyyzz = pbuffer.data(idx_ovl_hi + 12);

    auto ts_xxxxx_xxyzzz = pbuffer.data(idx_ovl_hi + 13);

    auto ts_xxxxx_xxzzzz = pbuffer.data(idx_ovl_hi + 14);

    auto ts_xxxxx_xyyyyy = pbuffer.data(idx_ovl_hi + 15);

    auto ts_xxxxx_xyyyyz = pbuffer.data(idx_ovl_hi + 16);

    auto ts_xxxxx_xyyyzz = pbuffer.data(idx_ovl_hi + 17);

    auto ts_xxxxx_xyyzzz = pbuffer.data(idx_ovl_hi + 18);

    auto ts_xxxxx_xyzzzz = pbuffer.data(idx_ovl_hi + 19);

    auto ts_xxxxx_xzzzzz = pbuffer.data(idx_ovl_hi + 20);

    auto ts_xxxxx_yyyyyy = pbuffer.data(idx_ovl_hi + 21);

    auto ts_xxxxx_yyyyyz = pbuffer.data(idx_ovl_hi + 22);

    auto ts_xxxxx_yyyyzz = pbuffer.data(idx_ovl_hi + 23);

    auto ts_xxxxx_yyyzzz = pbuffer.data(idx_ovl_hi + 24);

    auto ts_xxxxx_yyzzzz = pbuffer.data(idx_ovl_hi + 25);

    auto ts_xxxxx_yzzzzz = pbuffer.data(idx_ovl_hi + 26);

    auto ts_xxxxx_zzzzzz = pbuffer.data(idx_ovl_hi + 27);

    auto ts_xxxxy_xxxxxx = pbuffer.data(idx_ovl_hi + 28);

    auto ts_xxxxy_xxxxxy = pbuffer.data(idx_ovl_hi + 29);

    auto ts_xxxxy_xxxxxz = pbuffer.data(idx_ovl_hi + 30);

    auto ts_xxxxy_xxxxyy = pbuffer.data(idx_ovl_hi + 31);

    auto ts_xxxxy_xxxxyz = pbuffer.data(idx_ovl_hi + 32);

    auto ts_xxxxy_xxxxzz = pbuffer.data(idx_ovl_hi + 33);

    auto ts_xxxxy_xxxyyy = pbuffer.data(idx_ovl_hi + 34);

    auto ts_xxxxy_xxxyyz = pbuffer.data(idx_ovl_hi + 35);

    auto ts_xxxxy_xxxyzz = pbuffer.data(idx_ovl_hi + 36);

    auto ts_xxxxy_xxxzzz = pbuffer.data(idx_ovl_hi + 37);

    auto ts_xxxxy_xxyyyy = pbuffer.data(idx_ovl_hi + 38);

    auto ts_xxxxy_xxyyyz = pbuffer.data(idx_ovl_hi + 39);

    auto ts_xxxxy_xxyyzz = pbuffer.data(idx_ovl_hi + 40);

    auto ts_xxxxy_xxyzzz = pbuffer.data(idx_ovl_hi + 41);

    auto ts_xxxxy_xxzzzz = pbuffer.data(idx_ovl_hi + 42);

    auto ts_xxxxy_xyyyyy = pbuffer.data(idx_ovl_hi + 43);

    auto ts_xxxxy_xyyyyz = pbuffer.data(idx_ovl_hi + 44);

    auto ts_xxxxy_xyyyzz = pbuffer.data(idx_ovl_hi + 45);

    auto ts_xxxxy_xyyzzz = pbuffer.data(idx_ovl_hi + 46);

    auto ts_xxxxy_xyzzzz = pbuffer.data(idx_ovl_hi + 47);

    auto ts_xxxxy_xzzzzz = pbuffer.data(idx_ovl_hi + 48);

    auto ts_xxxxy_yyyyyy = pbuffer.data(idx_ovl_hi + 49);

    auto ts_xxxxy_yyyyyz = pbuffer.data(idx_ovl_hi + 50);

    auto ts_xxxxy_yyyyzz = pbuffer.data(idx_ovl_hi + 51);

    auto ts_xxxxy_yyyzzz = pbuffer.data(idx_ovl_hi + 52);

    auto ts_xxxxy_yyzzzz = pbuffer.data(idx_ovl_hi + 53);

    auto ts_xxxxy_yzzzzz = pbuffer.data(idx_ovl_hi + 54);

    auto ts_xxxxy_zzzzzz = pbuffer.data(idx_ovl_hi + 55);

    auto ts_xxxxz_xxxxxx = pbuffer.data(idx_ovl_hi + 56);

    auto ts_xxxxz_xxxxxy = pbuffer.data(idx_ovl_hi + 57);

    auto ts_xxxxz_xxxxxz = pbuffer.data(idx_ovl_hi + 58);

    auto ts_xxxxz_xxxxyy = pbuffer.data(idx_ovl_hi + 59);

    auto ts_xxxxz_xxxxyz = pbuffer.data(idx_ovl_hi + 60);

    auto ts_xxxxz_xxxxzz = pbuffer.data(idx_ovl_hi + 61);

    auto ts_xxxxz_xxxyyy = pbuffer.data(idx_ovl_hi + 62);

    auto ts_xxxxz_xxxyyz = pbuffer.data(idx_ovl_hi + 63);

    auto ts_xxxxz_xxxyzz = pbuffer.data(idx_ovl_hi + 64);

    auto ts_xxxxz_xxxzzz = pbuffer.data(idx_ovl_hi + 65);

    auto ts_xxxxz_xxyyyy = pbuffer.data(idx_ovl_hi + 66);

    auto ts_xxxxz_xxyyyz = pbuffer.data(idx_ovl_hi + 67);

    auto ts_xxxxz_xxyyzz = pbuffer.data(idx_ovl_hi + 68);

    auto ts_xxxxz_xxyzzz = pbuffer.data(idx_ovl_hi + 69);

    auto ts_xxxxz_xxzzzz = pbuffer.data(idx_ovl_hi + 70);

    auto ts_xxxxz_xyyyyy = pbuffer.data(idx_ovl_hi + 71);

    auto ts_xxxxz_xyyyyz = pbuffer.data(idx_ovl_hi + 72);

    auto ts_xxxxz_xyyyzz = pbuffer.data(idx_ovl_hi + 73);

    auto ts_xxxxz_xyyzzz = pbuffer.data(idx_ovl_hi + 74);

    auto ts_xxxxz_xyzzzz = pbuffer.data(idx_ovl_hi + 75);

    auto ts_xxxxz_xzzzzz = pbuffer.data(idx_ovl_hi + 76);

    auto ts_xxxxz_yyyyyy = pbuffer.data(idx_ovl_hi + 77);

    auto ts_xxxxz_yyyyyz = pbuffer.data(idx_ovl_hi + 78);

    auto ts_xxxxz_yyyyzz = pbuffer.data(idx_ovl_hi + 79);

    auto ts_xxxxz_yyyzzz = pbuffer.data(idx_ovl_hi + 80);

    auto ts_xxxxz_yyzzzz = pbuffer.data(idx_ovl_hi + 81);

    auto ts_xxxxz_yzzzzz = pbuffer.data(idx_ovl_hi + 82);

    auto ts_xxxxz_zzzzzz = pbuffer.data(idx_ovl_hi + 83);

    auto ts_xxxyy_xxxxxx = pbuffer.data(idx_ovl_hi + 84);

    auto ts_xxxyy_xxxxxy = pbuffer.data(idx_ovl_hi + 85);

    auto ts_xxxyy_xxxxxz = pbuffer.data(idx_ovl_hi + 86);

    auto ts_xxxyy_xxxxyy = pbuffer.data(idx_ovl_hi + 87);

    auto ts_xxxyy_xxxxyz = pbuffer.data(idx_ovl_hi + 88);

    auto ts_xxxyy_xxxxzz = pbuffer.data(idx_ovl_hi + 89);

    auto ts_xxxyy_xxxyyy = pbuffer.data(idx_ovl_hi + 90);

    auto ts_xxxyy_xxxyyz = pbuffer.data(idx_ovl_hi + 91);

    auto ts_xxxyy_xxxyzz = pbuffer.data(idx_ovl_hi + 92);

    auto ts_xxxyy_xxxzzz = pbuffer.data(idx_ovl_hi + 93);

    auto ts_xxxyy_xxyyyy = pbuffer.data(idx_ovl_hi + 94);

    auto ts_xxxyy_xxyyyz = pbuffer.data(idx_ovl_hi + 95);

    auto ts_xxxyy_xxyyzz = pbuffer.data(idx_ovl_hi + 96);

    auto ts_xxxyy_xxyzzz = pbuffer.data(idx_ovl_hi + 97);

    auto ts_xxxyy_xxzzzz = pbuffer.data(idx_ovl_hi + 98);

    auto ts_xxxyy_xyyyyy = pbuffer.data(idx_ovl_hi + 99);

    auto ts_xxxyy_xyyyyz = pbuffer.data(idx_ovl_hi + 100);

    auto ts_xxxyy_xyyyzz = pbuffer.data(idx_ovl_hi + 101);

    auto ts_xxxyy_xyyzzz = pbuffer.data(idx_ovl_hi + 102);

    auto ts_xxxyy_xyzzzz = pbuffer.data(idx_ovl_hi + 103);

    auto ts_xxxyy_xzzzzz = pbuffer.data(idx_ovl_hi + 104);

    auto ts_xxxyy_yyyyyy = pbuffer.data(idx_ovl_hi + 105);

    auto ts_xxxyy_yyyyyz = pbuffer.data(idx_ovl_hi + 106);

    auto ts_xxxyy_yyyyzz = pbuffer.data(idx_ovl_hi + 107);

    auto ts_xxxyy_yyyzzz = pbuffer.data(idx_ovl_hi + 108);

    auto ts_xxxyy_yyzzzz = pbuffer.data(idx_ovl_hi + 109);

    auto ts_xxxyy_yzzzzz = pbuffer.data(idx_ovl_hi + 110);

    auto ts_xxxyy_zzzzzz = pbuffer.data(idx_ovl_hi + 111);

    auto ts_xxxyz_xxxxxx = pbuffer.data(idx_ovl_hi + 112);

    auto ts_xxxyz_xxxxxy = pbuffer.data(idx_ovl_hi + 113);

    auto ts_xxxyz_xxxxxz = pbuffer.data(idx_ovl_hi + 114);

    auto ts_xxxyz_xxxxyy = pbuffer.data(idx_ovl_hi + 115);

    auto ts_xxxyz_xxxxyz = pbuffer.data(idx_ovl_hi + 116);

    auto ts_xxxyz_xxxxzz = pbuffer.data(idx_ovl_hi + 117);

    auto ts_xxxyz_xxxyyy = pbuffer.data(idx_ovl_hi + 118);

    auto ts_xxxyz_xxxyyz = pbuffer.data(idx_ovl_hi + 119);

    auto ts_xxxyz_xxxyzz = pbuffer.data(idx_ovl_hi + 120);

    auto ts_xxxyz_xxxzzz = pbuffer.data(idx_ovl_hi + 121);

    auto ts_xxxyz_xxyyyy = pbuffer.data(idx_ovl_hi + 122);

    auto ts_xxxyz_xxyyyz = pbuffer.data(idx_ovl_hi + 123);

    auto ts_xxxyz_xxyyzz = pbuffer.data(idx_ovl_hi + 124);

    auto ts_xxxyz_xxyzzz = pbuffer.data(idx_ovl_hi + 125);

    auto ts_xxxyz_xxzzzz = pbuffer.data(idx_ovl_hi + 126);

    auto ts_xxxyz_xyyyyy = pbuffer.data(idx_ovl_hi + 127);

    auto ts_xxxyz_xyyyyz = pbuffer.data(idx_ovl_hi + 128);

    auto ts_xxxyz_xyyyzz = pbuffer.data(idx_ovl_hi + 129);

    auto ts_xxxyz_xyyzzz = pbuffer.data(idx_ovl_hi + 130);

    auto ts_xxxyz_xyzzzz = pbuffer.data(idx_ovl_hi + 131);

    auto ts_xxxyz_xzzzzz = pbuffer.data(idx_ovl_hi + 132);

    auto ts_xxxyz_yyyyyy = pbuffer.data(idx_ovl_hi + 133);

    auto ts_xxxyz_yyyyyz = pbuffer.data(idx_ovl_hi + 134);

    auto ts_xxxyz_yyyyzz = pbuffer.data(idx_ovl_hi + 135);

    auto ts_xxxyz_yyyzzz = pbuffer.data(idx_ovl_hi + 136);

    auto ts_xxxyz_yyzzzz = pbuffer.data(idx_ovl_hi + 137);

    auto ts_xxxyz_yzzzzz = pbuffer.data(idx_ovl_hi + 138);

    auto ts_xxxyz_zzzzzz = pbuffer.data(idx_ovl_hi + 139);

    auto ts_xxxzz_xxxxxx = pbuffer.data(idx_ovl_hi + 140);

    auto ts_xxxzz_xxxxxy = pbuffer.data(idx_ovl_hi + 141);

    auto ts_xxxzz_xxxxxz = pbuffer.data(idx_ovl_hi + 142);

    auto ts_xxxzz_xxxxyy = pbuffer.data(idx_ovl_hi + 143);

    auto ts_xxxzz_xxxxyz = pbuffer.data(idx_ovl_hi + 144);

    auto ts_xxxzz_xxxxzz = pbuffer.data(idx_ovl_hi + 145);

    auto ts_xxxzz_xxxyyy = pbuffer.data(idx_ovl_hi + 146);

    auto ts_xxxzz_xxxyyz = pbuffer.data(idx_ovl_hi + 147);

    auto ts_xxxzz_xxxyzz = pbuffer.data(idx_ovl_hi + 148);

    auto ts_xxxzz_xxxzzz = pbuffer.data(idx_ovl_hi + 149);

    auto ts_xxxzz_xxyyyy = pbuffer.data(idx_ovl_hi + 150);

    auto ts_xxxzz_xxyyyz = pbuffer.data(idx_ovl_hi + 151);

    auto ts_xxxzz_xxyyzz = pbuffer.data(idx_ovl_hi + 152);

    auto ts_xxxzz_xxyzzz = pbuffer.data(idx_ovl_hi + 153);

    auto ts_xxxzz_xxzzzz = pbuffer.data(idx_ovl_hi + 154);

    auto ts_xxxzz_xyyyyy = pbuffer.data(idx_ovl_hi + 155);

    auto ts_xxxzz_xyyyyz = pbuffer.data(idx_ovl_hi + 156);

    auto ts_xxxzz_xyyyzz = pbuffer.data(idx_ovl_hi + 157);

    auto ts_xxxzz_xyyzzz = pbuffer.data(idx_ovl_hi + 158);

    auto ts_xxxzz_xyzzzz = pbuffer.data(idx_ovl_hi + 159);

    auto ts_xxxzz_xzzzzz = pbuffer.data(idx_ovl_hi + 160);

    auto ts_xxxzz_yyyyyy = pbuffer.data(idx_ovl_hi + 161);

    auto ts_xxxzz_yyyyyz = pbuffer.data(idx_ovl_hi + 162);

    auto ts_xxxzz_yyyyzz = pbuffer.data(idx_ovl_hi + 163);

    auto ts_xxxzz_yyyzzz = pbuffer.data(idx_ovl_hi + 164);

    auto ts_xxxzz_yyzzzz = pbuffer.data(idx_ovl_hi + 165);

    auto ts_xxxzz_yzzzzz = pbuffer.data(idx_ovl_hi + 166);

    auto ts_xxxzz_zzzzzz = pbuffer.data(idx_ovl_hi + 167);

    auto ts_xxyyy_xxxxxx = pbuffer.data(idx_ovl_hi + 168);

    auto ts_xxyyy_xxxxxy = pbuffer.data(idx_ovl_hi + 169);

    auto ts_xxyyy_xxxxxz = pbuffer.data(idx_ovl_hi + 170);

    auto ts_xxyyy_xxxxyy = pbuffer.data(idx_ovl_hi + 171);

    auto ts_xxyyy_xxxxyz = pbuffer.data(idx_ovl_hi + 172);

    auto ts_xxyyy_xxxxzz = pbuffer.data(idx_ovl_hi + 173);

    auto ts_xxyyy_xxxyyy = pbuffer.data(idx_ovl_hi + 174);

    auto ts_xxyyy_xxxyyz = pbuffer.data(idx_ovl_hi + 175);

    auto ts_xxyyy_xxxyzz = pbuffer.data(idx_ovl_hi + 176);

    auto ts_xxyyy_xxxzzz = pbuffer.data(idx_ovl_hi + 177);

    auto ts_xxyyy_xxyyyy = pbuffer.data(idx_ovl_hi + 178);

    auto ts_xxyyy_xxyyyz = pbuffer.data(idx_ovl_hi + 179);

    auto ts_xxyyy_xxyyzz = pbuffer.data(idx_ovl_hi + 180);

    auto ts_xxyyy_xxyzzz = pbuffer.data(idx_ovl_hi + 181);

    auto ts_xxyyy_xxzzzz = pbuffer.data(idx_ovl_hi + 182);

    auto ts_xxyyy_xyyyyy = pbuffer.data(idx_ovl_hi + 183);

    auto ts_xxyyy_xyyyyz = pbuffer.data(idx_ovl_hi + 184);

    auto ts_xxyyy_xyyyzz = pbuffer.data(idx_ovl_hi + 185);

    auto ts_xxyyy_xyyzzz = pbuffer.data(idx_ovl_hi + 186);

    auto ts_xxyyy_xyzzzz = pbuffer.data(idx_ovl_hi + 187);

    auto ts_xxyyy_xzzzzz = pbuffer.data(idx_ovl_hi + 188);

    auto ts_xxyyy_yyyyyy = pbuffer.data(idx_ovl_hi + 189);

    auto ts_xxyyy_yyyyyz = pbuffer.data(idx_ovl_hi + 190);

    auto ts_xxyyy_yyyyzz = pbuffer.data(idx_ovl_hi + 191);

    auto ts_xxyyy_yyyzzz = pbuffer.data(idx_ovl_hi + 192);

    auto ts_xxyyy_yyzzzz = pbuffer.data(idx_ovl_hi + 193);

    auto ts_xxyyy_yzzzzz = pbuffer.data(idx_ovl_hi + 194);

    auto ts_xxyyy_zzzzzz = pbuffer.data(idx_ovl_hi + 195);

    auto ts_xxyyz_xxxxxx = pbuffer.data(idx_ovl_hi + 196);

    auto ts_xxyyz_xxxxxy = pbuffer.data(idx_ovl_hi + 197);

    auto ts_xxyyz_xxxxxz = pbuffer.data(idx_ovl_hi + 198);

    auto ts_xxyyz_xxxxyy = pbuffer.data(idx_ovl_hi + 199);

    auto ts_xxyyz_xxxxyz = pbuffer.data(idx_ovl_hi + 200);

    auto ts_xxyyz_xxxxzz = pbuffer.data(idx_ovl_hi + 201);

    auto ts_xxyyz_xxxyyy = pbuffer.data(idx_ovl_hi + 202);

    auto ts_xxyyz_xxxyyz = pbuffer.data(idx_ovl_hi + 203);

    auto ts_xxyyz_xxxyzz = pbuffer.data(idx_ovl_hi + 204);

    auto ts_xxyyz_xxxzzz = pbuffer.data(idx_ovl_hi + 205);

    auto ts_xxyyz_xxyyyy = pbuffer.data(idx_ovl_hi + 206);

    auto ts_xxyyz_xxyyyz = pbuffer.data(idx_ovl_hi + 207);

    auto ts_xxyyz_xxyyzz = pbuffer.data(idx_ovl_hi + 208);

    auto ts_xxyyz_xxyzzz = pbuffer.data(idx_ovl_hi + 209);

    auto ts_xxyyz_xxzzzz = pbuffer.data(idx_ovl_hi + 210);

    auto ts_xxyyz_xyyyyy = pbuffer.data(idx_ovl_hi + 211);

    auto ts_xxyyz_xyyyyz = pbuffer.data(idx_ovl_hi + 212);

    auto ts_xxyyz_xyyyzz = pbuffer.data(idx_ovl_hi + 213);

    auto ts_xxyyz_xyyzzz = pbuffer.data(idx_ovl_hi + 214);

    auto ts_xxyyz_xyzzzz = pbuffer.data(idx_ovl_hi + 215);

    auto ts_xxyyz_xzzzzz = pbuffer.data(idx_ovl_hi + 216);

    auto ts_xxyyz_yyyyyy = pbuffer.data(idx_ovl_hi + 217);

    auto ts_xxyyz_yyyyyz = pbuffer.data(idx_ovl_hi + 218);

    auto ts_xxyyz_yyyyzz = pbuffer.data(idx_ovl_hi + 219);

    auto ts_xxyyz_yyyzzz = pbuffer.data(idx_ovl_hi + 220);

    auto ts_xxyyz_yyzzzz = pbuffer.data(idx_ovl_hi + 221);

    auto ts_xxyyz_yzzzzz = pbuffer.data(idx_ovl_hi + 222);

    auto ts_xxyyz_zzzzzz = pbuffer.data(idx_ovl_hi + 223);

    auto ts_xxyzz_xxxxxx = pbuffer.data(idx_ovl_hi + 224);

    auto ts_xxyzz_xxxxxy = pbuffer.data(idx_ovl_hi + 225);

    auto ts_xxyzz_xxxxxz = pbuffer.data(idx_ovl_hi + 226);

    auto ts_xxyzz_xxxxyy = pbuffer.data(idx_ovl_hi + 227);

    auto ts_xxyzz_xxxxyz = pbuffer.data(idx_ovl_hi + 228);

    auto ts_xxyzz_xxxxzz = pbuffer.data(idx_ovl_hi + 229);

    auto ts_xxyzz_xxxyyy = pbuffer.data(idx_ovl_hi + 230);

    auto ts_xxyzz_xxxyyz = pbuffer.data(idx_ovl_hi + 231);

    auto ts_xxyzz_xxxyzz = pbuffer.data(idx_ovl_hi + 232);

    auto ts_xxyzz_xxxzzz = pbuffer.data(idx_ovl_hi + 233);

    auto ts_xxyzz_xxyyyy = pbuffer.data(idx_ovl_hi + 234);

    auto ts_xxyzz_xxyyyz = pbuffer.data(idx_ovl_hi + 235);

    auto ts_xxyzz_xxyyzz = pbuffer.data(idx_ovl_hi + 236);

    auto ts_xxyzz_xxyzzz = pbuffer.data(idx_ovl_hi + 237);

    auto ts_xxyzz_xxzzzz = pbuffer.data(idx_ovl_hi + 238);

    auto ts_xxyzz_xyyyyy = pbuffer.data(idx_ovl_hi + 239);

    auto ts_xxyzz_xyyyyz = pbuffer.data(idx_ovl_hi + 240);

    auto ts_xxyzz_xyyyzz = pbuffer.data(idx_ovl_hi + 241);

    auto ts_xxyzz_xyyzzz = pbuffer.data(idx_ovl_hi + 242);

    auto ts_xxyzz_xyzzzz = pbuffer.data(idx_ovl_hi + 243);

    auto ts_xxyzz_xzzzzz = pbuffer.data(idx_ovl_hi + 244);

    auto ts_xxyzz_yyyyyy = pbuffer.data(idx_ovl_hi + 245);

    auto ts_xxyzz_yyyyyz = pbuffer.data(idx_ovl_hi + 246);

    auto ts_xxyzz_yyyyzz = pbuffer.data(idx_ovl_hi + 247);

    auto ts_xxyzz_yyyzzz = pbuffer.data(idx_ovl_hi + 248);

    auto ts_xxyzz_yyzzzz = pbuffer.data(idx_ovl_hi + 249);

    auto ts_xxyzz_yzzzzz = pbuffer.data(idx_ovl_hi + 250);

    auto ts_xxyzz_zzzzzz = pbuffer.data(idx_ovl_hi + 251);

    auto ts_xxzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 252);

    auto ts_xxzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 253);

    auto ts_xxzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 254);

    auto ts_xxzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 255);

    auto ts_xxzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 256);

    auto ts_xxzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 257);

    auto ts_xxzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 258);

    auto ts_xxzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 259);

    auto ts_xxzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 260);

    auto ts_xxzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 261);

    auto ts_xxzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 262);

    auto ts_xxzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 263);

    auto ts_xxzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 264);

    auto ts_xxzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 265);

    auto ts_xxzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 266);

    auto ts_xxzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 267);

    auto ts_xxzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 268);

    auto ts_xxzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 269);

    auto ts_xxzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 270);

    auto ts_xxzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 271);

    auto ts_xxzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 272);

    auto ts_xxzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 273);

    auto ts_xxzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 274);

    auto ts_xxzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 275);

    auto ts_xxzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 276);

    auto ts_xxzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 277);

    auto ts_xxzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 278);

    auto ts_xxzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 279);

    auto ts_xyyyy_xxxxxx = pbuffer.data(idx_ovl_hi + 280);

    auto ts_xyyyy_xxxxxy = pbuffer.data(idx_ovl_hi + 281);

    auto ts_xyyyy_xxxxxz = pbuffer.data(idx_ovl_hi + 282);

    auto ts_xyyyy_xxxxyy = pbuffer.data(idx_ovl_hi + 283);

    auto ts_xyyyy_xxxxyz = pbuffer.data(idx_ovl_hi + 284);

    auto ts_xyyyy_xxxxzz = pbuffer.data(idx_ovl_hi + 285);

    auto ts_xyyyy_xxxyyy = pbuffer.data(idx_ovl_hi + 286);

    auto ts_xyyyy_xxxyyz = pbuffer.data(idx_ovl_hi + 287);

    auto ts_xyyyy_xxxyzz = pbuffer.data(idx_ovl_hi + 288);

    auto ts_xyyyy_xxxzzz = pbuffer.data(idx_ovl_hi + 289);

    auto ts_xyyyy_xxyyyy = pbuffer.data(idx_ovl_hi + 290);

    auto ts_xyyyy_xxyyyz = pbuffer.data(idx_ovl_hi + 291);

    auto ts_xyyyy_xxyyzz = pbuffer.data(idx_ovl_hi + 292);

    auto ts_xyyyy_xxyzzz = pbuffer.data(idx_ovl_hi + 293);

    auto ts_xyyyy_xxzzzz = pbuffer.data(idx_ovl_hi + 294);

    auto ts_xyyyy_xyyyyy = pbuffer.data(idx_ovl_hi + 295);

    auto ts_xyyyy_xyyyyz = pbuffer.data(idx_ovl_hi + 296);

    auto ts_xyyyy_xyyyzz = pbuffer.data(idx_ovl_hi + 297);

    auto ts_xyyyy_xyyzzz = pbuffer.data(idx_ovl_hi + 298);

    auto ts_xyyyy_xyzzzz = pbuffer.data(idx_ovl_hi + 299);

    auto ts_xyyyy_xzzzzz = pbuffer.data(idx_ovl_hi + 300);

    auto ts_xyyyy_yyyyyy = pbuffer.data(idx_ovl_hi + 301);

    auto ts_xyyyy_yyyyyz = pbuffer.data(idx_ovl_hi + 302);

    auto ts_xyyyy_yyyyzz = pbuffer.data(idx_ovl_hi + 303);

    auto ts_xyyyy_yyyzzz = pbuffer.data(idx_ovl_hi + 304);

    auto ts_xyyyy_yyzzzz = pbuffer.data(idx_ovl_hi + 305);

    auto ts_xyyyy_yzzzzz = pbuffer.data(idx_ovl_hi + 306);

    auto ts_xyyyy_zzzzzz = pbuffer.data(idx_ovl_hi + 307);

    auto ts_xyyyz_xxxxxx = pbuffer.data(idx_ovl_hi + 308);

    auto ts_xyyyz_xxxxxy = pbuffer.data(idx_ovl_hi + 309);

    auto ts_xyyyz_xxxxxz = pbuffer.data(idx_ovl_hi + 310);

    auto ts_xyyyz_xxxxyy = pbuffer.data(idx_ovl_hi + 311);

    auto ts_xyyyz_xxxxyz = pbuffer.data(idx_ovl_hi + 312);

    auto ts_xyyyz_xxxxzz = pbuffer.data(idx_ovl_hi + 313);

    auto ts_xyyyz_xxxyyy = pbuffer.data(idx_ovl_hi + 314);

    auto ts_xyyyz_xxxyyz = pbuffer.data(idx_ovl_hi + 315);

    auto ts_xyyyz_xxxyzz = pbuffer.data(idx_ovl_hi + 316);

    auto ts_xyyyz_xxxzzz = pbuffer.data(idx_ovl_hi + 317);

    auto ts_xyyyz_xxyyyy = pbuffer.data(idx_ovl_hi + 318);

    auto ts_xyyyz_xxyyyz = pbuffer.data(idx_ovl_hi + 319);

    auto ts_xyyyz_xxyyzz = pbuffer.data(idx_ovl_hi + 320);

    auto ts_xyyyz_xxyzzz = pbuffer.data(idx_ovl_hi + 321);

    auto ts_xyyyz_xxzzzz = pbuffer.data(idx_ovl_hi + 322);

    auto ts_xyyyz_xyyyyy = pbuffer.data(idx_ovl_hi + 323);

    auto ts_xyyyz_xyyyyz = pbuffer.data(idx_ovl_hi + 324);

    auto ts_xyyyz_xyyyzz = pbuffer.data(idx_ovl_hi + 325);

    auto ts_xyyyz_xyyzzz = pbuffer.data(idx_ovl_hi + 326);

    auto ts_xyyyz_xyzzzz = pbuffer.data(idx_ovl_hi + 327);

    auto ts_xyyyz_xzzzzz = pbuffer.data(idx_ovl_hi + 328);

    auto ts_xyyyz_yyyyyy = pbuffer.data(idx_ovl_hi + 329);

    auto ts_xyyyz_yyyyyz = pbuffer.data(idx_ovl_hi + 330);

    auto ts_xyyyz_yyyyzz = pbuffer.data(idx_ovl_hi + 331);

    auto ts_xyyyz_yyyzzz = pbuffer.data(idx_ovl_hi + 332);

    auto ts_xyyyz_yyzzzz = pbuffer.data(idx_ovl_hi + 333);

    auto ts_xyyyz_yzzzzz = pbuffer.data(idx_ovl_hi + 334);

    auto ts_xyyyz_zzzzzz = pbuffer.data(idx_ovl_hi + 335);

    auto ts_xyyzz_xxxxxx = pbuffer.data(idx_ovl_hi + 336);

    auto ts_xyyzz_xxxxxy = pbuffer.data(idx_ovl_hi + 337);

    auto ts_xyyzz_xxxxxz = pbuffer.data(idx_ovl_hi + 338);

    auto ts_xyyzz_xxxxyy = pbuffer.data(idx_ovl_hi + 339);

    auto ts_xyyzz_xxxxyz = pbuffer.data(idx_ovl_hi + 340);

    auto ts_xyyzz_xxxxzz = pbuffer.data(idx_ovl_hi + 341);

    auto ts_xyyzz_xxxyyy = pbuffer.data(idx_ovl_hi + 342);

    auto ts_xyyzz_xxxyyz = pbuffer.data(idx_ovl_hi + 343);

    auto ts_xyyzz_xxxyzz = pbuffer.data(idx_ovl_hi + 344);

    auto ts_xyyzz_xxxzzz = pbuffer.data(idx_ovl_hi + 345);

    auto ts_xyyzz_xxyyyy = pbuffer.data(idx_ovl_hi + 346);

    auto ts_xyyzz_xxyyyz = pbuffer.data(idx_ovl_hi + 347);

    auto ts_xyyzz_xxyyzz = pbuffer.data(idx_ovl_hi + 348);

    auto ts_xyyzz_xxyzzz = pbuffer.data(idx_ovl_hi + 349);

    auto ts_xyyzz_xxzzzz = pbuffer.data(idx_ovl_hi + 350);

    auto ts_xyyzz_xyyyyy = pbuffer.data(idx_ovl_hi + 351);

    auto ts_xyyzz_xyyyyz = pbuffer.data(idx_ovl_hi + 352);

    auto ts_xyyzz_xyyyzz = pbuffer.data(idx_ovl_hi + 353);

    auto ts_xyyzz_xyyzzz = pbuffer.data(idx_ovl_hi + 354);

    auto ts_xyyzz_xyzzzz = pbuffer.data(idx_ovl_hi + 355);

    auto ts_xyyzz_xzzzzz = pbuffer.data(idx_ovl_hi + 356);

    auto ts_xyyzz_yyyyyy = pbuffer.data(idx_ovl_hi + 357);

    auto ts_xyyzz_yyyyyz = pbuffer.data(idx_ovl_hi + 358);

    auto ts_xyyzz_yyyyzz = pbuffer.data(idx_ovl_hi + 359);

    auto ts_xyyzz_yyyzzz = pbuffer.data(idx_ovl_hi + 360);

    auto ts_xyyzz_yyzzzz = pbuffer.data(idx_ovl_hi + 361);

    auto ts_xyyzz_yzzzzz = pbuffer.data(idx_ovl_hi + 362);

    auto ts_xyyzz_zzzzzz = pbuffer.data(idx_ovl_hi + 363);

    auto ts_xyzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 364);

    auto ts_xyzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 365);

    auto ts_xyzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 366);

    auto ts_xyzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 367);

    auto ts_xyzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 368);

    auto ts_xyzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 369);

    auto ts_xyzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 370);

    auto ts_xyzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 371);

    auto ts_xyzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 372);

    auto ts_xyzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 373);

    auto ts_xyzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 374);

    auto ts_xyzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 375);

    auto ts_xyzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 376);

    auto ts_xyzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 377);

    auto ts_xyzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 378);

    auto ts_xyzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 379);

    auto ts_xyzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 380);

    auto ts_xyzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 381);

    auto ts_xyzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 382);

    auto ts_xyzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 383);

    auto ts_xyzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 384);

    auto ts_xyzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 385);

    auto ts_xyzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 386);

    auto ts_xyzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 387);

    auto ts_xyzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 388);

    auto ts_xyzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 389);

    auto ts_xyzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 390);

    auto ts_xyzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 391);

    auto ts_xzzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 392);

    auto ts_xzzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 393);

    auto ts_xzzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 394);

    auto ts_xzzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 395);

    auto ts_xzzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 396);

    auto ts_xzzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 397);

    auto ts_xzzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 398);

    auto ts_xzzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 399);

    auto ts_xzzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 400);

    auto ts_xzzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 401);

    auto ts_xzzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 402);

    auto ts_xzzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 403);

    auto ts_xzzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 404);

    auto ts_xzzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 405);

    auto ts_xzzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 406);

    auto ts_xzzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 407);

    auto ts_xzzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 408);

    auto ts_xzzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 409);

    auto ts_xzzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 410);

    auto ts_xzzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 411);

    auto ts_xzzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 412);

    auto ts_xzzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 413);

    auto ts_xzzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 414);

    auto ts_xzzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 415);

    auto ts_xzzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 416);

    auto ts_xzzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 417);

    auto ts_xzzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 418);

    auto ts_xzzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 419);

    auto ts_yyyyy_xxxxxx = pbuffer.data(idx_ovl_hi + 420);

    auto ts_yyyyy_xxxxxy = pbuffer.data(idx_ovl_hi + 421);

    auto ts_yyyyy_xxxxxz = pbuffer.data(idx_ovl_hi + 422);

    auto ts_yyyyy_xxxxyy = pbuffer.data(idx_ovl_hi + 423);

    auto ts_yyyyy_xxxxyz = pbuffer.data(idx_ovl_hi + 424);

    auto ts_yyyyy_xxxxzz = pbuffer.data(idx_ovl_hi + 425);

    auto ts_yyyyy_xxxyyy = pbuffer.data(idx_ovl_hi + 426);

    auto ts_yyyyy_xxxyyz = pbuffer.data(idx_ovl_hi + 427);

    auto ts_yyyyy_xxxyzz = pbuffer.data(idx_ovl_hi + 428);

    auto ts_yyyyy_xxxzzz = pbuffer.data(idx_ovl_hi + 429);

    auto ts_yyyyy_xxyyyy = pbuffer.data(idx_ovl_hi + 430);

    auto ts_yyyyy_xxyyyz = pbuffer.data(idx_ovl_hi + 431);

    auto ts_yyyyy_xxyyzz = pbuffer.data(idx_ovl_hi + 432);

    auto ts_yyyyy_xxyzzz = pbuffer.data(idx_ovl_hi + 433);

    auto ts_yyyyy_xxzzzz = pbuffer.data(idx_ovl_hi + 434);

    auto ts_yyyyy_xyyyyy = pbuffer.data(idx_ovl_hi + 435);

    auto ts_yyyyy_xyyyyz = pbuffer.data(idx_ovl_hi + 436);

    auto ts_yyyyy_xyyyzz = pbuffer.data(idx_ovl_hi + 437);

    auto ts_yyyyy_xyyzzz = pbuffer.data(idx_ovl_hi + 438);

    auto ts_yyyyy_xyzzzz = pbuffer.data(idx_ovl_hi + 439);

    auto ts_yyyyy_xzzzzz = pbuffer.data(idx_ovl_hi + 440);

    auto ts_yyyyy_yyyyyy = pbuffer.data(idx_ovl_hi + 441);

    auto ts_yyyyy_yyyyyz = pbuffer.data(idx_ovl_hi + 442);

    auto ts_yyyyy_yyyyzz = pbuffer.data(idx_ovl_hi + 443);

    auto ts_yyyyy_yyyzzz = pbuffer.data(idx_ovl_hi + 444);

    auto ts_yyyyy_yyzzzz = pbuffer.data(idx_ovl_hi + 445);

    auto ts_yyyyy_yzzzzz = pbuffer.data(idx_ovl_hi + 446);

    auto ts_yyyyy_zzzzzz = pbuffer.data(idx_ovl_hi + 447);

    auto ts_yyyyz_xxxxxx = pbuffer.data(idx_ovl_hi + 448);

    auto ts_yyyyz_xxxxxy = pbuffer.data(idx_ovl_hi + 449);

    auto ts_yyyyz_xxxxxz = pbuffer.data(idx_ovl_hi + 450);

    auto ts_yyyyz_xxxxyy = pbuffer.data(idx_ovl_hi + 451);

    auto ts_yyyyz_xxxxyz = pbuffer.data(idx_ovl_hi + 452);

    auto ts_yyyyz_xxxxzz = pbuffer.data(idx_ovl_hi + 453);

    auto ts_yyyyz_xxxyyy = pbuffer.data(idx_ovl_hi + 454);

    auto ts_yyyyz_xxxyyz = pbuffer.data(idx_ovl_hi + 455);

    auto ts_yyyyz_xxxyzz = pbuffer.data(idx_ovl_hi + 456);

    auto ts_yyyyz_xxxzzz = pbuffer.data(idx_ovl_hi + 457);

    auto ts_yyyyz_xxyyyy = pbuffer.data(idx_ovl_hi + 458);

    auto ts_yyyyz_xxyyyz = pbuffer.data(idx_ovl_hi + 459);

    auto ts_yyyyz_xxyyzz = pbuffer.data(idx_ovl_hi + 460);

    auto ts_yyyyz_xxyzzz = pbuffer.data(idx_ovl_hi + 461);

    auto ts_yyyyz_xxzzzz = pbuffer.data(idx_ovl_hi + 462);

    auto ts_yyyyz_xyyyyy = pbuffer.data(idx_ovl_hi + 463);

    auto ts_yyyyz_xyyyyz = pbuffer.data(idx_ovl_hi + 464);

    auto ts_yyyyz_xyyyzz = pbuffer.data(idx_ovl_hi + 465);

    auto ts_yyyyz_xyyzzz = pbuffer.data(idx_ovl_hi + 466);

    auto ts_yyyyz_xyzzzz = pbuffer.data(idx_ovl_hi + 467);

    auto ts_yyyyz_xzzzzz = pbuffer.data(idx_ovl_hi + 468);

    auto ts_yyyyz_yyyyyy = pbuffer.data(idx_ovl_hi + 469);

    auto ts_yyyyz_yyyyyz = pbuffer.data(idx_ovl_hi + 470);

    auto ts_yyyyz_yyyyzz = pbuffer.data(idx_ovl_hi + 471);

    auto ts_yyyyz_yyyzzz = pbuffer.data(idx_ovl_hi + 472);

    auto ts_yyyyz_yyzzzz = pbuffer.data(idx_ovl_hi + 473);

    auto ts_yyyyz_yzzzzz = pbuffer.data(idx_ovl_hi + 474);

    auto ts_yyyyz_zzzzzz = pbuffer.data(idx_ovl_hi + 475);

    auto ts_yyyzz_xxxxxx = pbuffer.data(idx_ovl_hi + 476);

    auto ts_yyyzz_xxxxxy = pbuffer.data(idx_ovl_hi + 477);

    auto ts_yyyzz_xxxxxz = pbuffer.data(idx_ovl_hi + 478);

    auto ts_yyyzz_xxxxyy = pbuffer.data(idx_ovl_hi + 479);

    auto ts_yyyzz_xxxxyz = pbuffer.data(idx_ovl_hi + 480);

    auto ts_yyyzz_xxxxzz = pbuffer.data(idx_ovl_hi + 481);

    auto ts_yyyzz_xxxyyy = pbuffer.data(idx_ovl_hi + 482);

    auto ts_yyyzz_xxxyyz = pbuffer.data(idx_ovl_hi + 483);

    auto ts_yyyzz_xxxyzz = pbuffer.data(idx_ovl_hi + 484);

    auto ts_yyyzz_xxxzzz = pbuffer.data(idx_ovl_hi + 485);

    auto ts_yyyzz_xxyyyy = pbuffer.data(idx_ovl_hi + 486);

    auto ts_yyyzz_xxyyyz = pbuffer.data(idx_ovl_hi + 487);

    auto ts_yyyzz_xxyyzz = pbuffer.data(idx_ovl_hi + 488);

    auto ts_yyyzz_xxyzzz = pbuffer.data(idx_ovl_hi + 489);

    auto ts_yyyzz_xxzzzz = pbuffer.data(idx_ovl_hi + 490);

    auto ts_yyyzz_xyyyyy = pbuffer.data(idx_ovl_hi + 491);

    auto ts_yyyzz_xyyyyz = pbuffer.data(idx_ovl_hi + 492);

    auto ts_yyyzz_xyyyzz = pbuffer.data(idx_ovl_hi + 493);

    auto ts_yyyzz_xyyzzz = pbuffer.data(idx_ovl_hi + 494);

    auto ts_yyyzz_xyzzzz = pbuffer.data(idx_ovl_hi + 495);

    auto ts_yyyzz_xzzzzz = pbuffer.data(idx_ovl_hi + 496);

    auto ts_yyyzz_yyyyyy = pbuffer.data(idx_ovl_hi + 497);

    auto ts_yyyzz_yyyyyz = pbuffer.data(idx_ovl_hi + 498);

    auto ts_yyyzz_yyyyzz = pbuffer.data(idx_ovl_hi + 499);

    auto ts_yyyzz_yyyzzz = pbuffer.data(idx_ovl_hi + 500);

    auto ts_yyyzz_yyzzzz = pbuffer.data(idx_ovl_hi + 501);

    auto ts_yyyzz_yzzzzz = pbuffer.data(idx_ovl_hi + 502);

    auto ts_yyyzz_zzzzzz = pbuffer.data(idx_ovl_hi + 503);

    auto ts_yyzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 504);

    auto ts_yyzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 505);

    auto ts_yyzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 506);

    auto ts_yyzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 507);

    auto ts_yyzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 508);

    auto ts_yyzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 509);

    auto ts_yyzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 510);

    auto ts_yyzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 511);

    auto ts_yyzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 512);

    auto ts_yyzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 513);

    auto ts_yyzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 514);

    auto ts_yyzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 515);

    auto ts_yyzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 516);

    auto ts_yyzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 517);

    auto ts_yyzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 518);

    auto ts_yyzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 519);

    auto ts_yyzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 520);

    auto ts_yyzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 521);

    auto ts_yyzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 522);

    auto ts_yyzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 523);

    auto ts_yyzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 524);

    auto ts_yyzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 525);

    auto ts_yyzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 526);

    auto ts_yyzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 527);

    auto ts_yyzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 528);

    auto ts_yyzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 529);

    auto ts_yyzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 530);

    auto ts_yyzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 531);

    auto ts_yzzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 532);

    auto ts_yzzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 533);

    auto ts_yzzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 534);

    auto ts_yzzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 535);

    auto ts_yzzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 536);

    auto ts_yzzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 537);

    auto ts_yzzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 538);

    auto ts_yzzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 539);

    auto ts_yzzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 540);

    auto ts_yzzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 541);

    auto ts_yzzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 542);

    auto ts_yzzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 543);

    auto ts_yzzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 544);

    auto ts_yzzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 545);

    auto ts_yzzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 546);

    auto ts_yzzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 547);

    auto ts_yzzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 548);

    auto ts_yzzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 549);

    auto ts_yzzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 550);

    auto ts_yzzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 551);

    auto ts_yzzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 552);

    auto ts_yzzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 553);

    auto ts_yzzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 554);

    auto ts_yzzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 555);

    auto ts_yzzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 556);

    auto ts_yzzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 557);

    auto ts_yzzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 558);

    auto ts_yzzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 559);

    auto ts_zzzzz_xxxxxx = pbuffer.data(idx_ovl_hi + 560);

    auto ts_zzzzz_xxxxxy = pbuffer.data(idx_ovl_hi + 561);

    auto ts_zzzzz_xxxxxz = pbuffer.data(idx_ovl_hi + 562);

    auto ts_zzzzz_xxxxyy = pbuffer.data(idx_ovl_hi + 563);

    auto ts_zzzzz_xxxxyz = pbuffer.data(idx_ovl_hi + 564);

    auto ts_zzzzz_xxxxzz = pbuffer.data(idx_ovl_hi + 565);

    auto ts_zzzzz_xxxyyy = pbuffer.data(idx_ovl_hi + 566);

    auto ts_zzzzz_xxxyyz = pbuffer.data(idx_ovl_hi + 567);

    auto ts_zzzzz_xxxyzz = pbuffer.data(idx_ovl_hi + 568);

    auto ts_zzzzz_xxxzzz = pbuffer.data(idx_ovl_hi + 569);

    auto ts_zzzzz_xxyyyy = pbuffer.data(idx_ovl_hi + 570);

    auto ts_zzzzz_xxyyyz = pbuffer.data(idx_ovl_hi + 571);

    auto ts_zzzzz_xxyyzz = pbuffer.data(idx_ovl_hi + 572);

    auto ts_zzzzz_xxyzzz = pbuffer.data(idx_ovl_hi + 573);

    auto ts_zzzzz_xxzzzz = pbuffer.data(idx_ovl_hi + 574);

    auto ts_zzzzz_xyyyyy = pbuffer.data(idx_ovl_hi + 575);

    auto ts_zzzzz_xyyyyz = pbuffer.data(idx_ovl_hi + 576);

    auto ts_zzzzz_xyyyzz = pbuffer.data(idx_ovl_hi + 577);

    auto ts_zzzzz_xyyzzz = pbuffer.data(idx_ovl_hi + 578);

    auto ts_zzzzz_xyzzzz = pbuffer.data(idx_ovl_hi + 579);

    auto ts_zzzzz_xzzzzz = pbuffer.data(idx_ovl_hi + 580);

    auto ts_zzzzz_yyyyyy = pbuffer.data(idx_ovl_hi + 581);

    auto ts_zzzzz_yyyyyz = pbuffer.data(idx_ovl_hi + 582);

    auto ts_zzzzz_yyyyzz = pbuffer.data(idx_ovl_hi + 583);

    auto ts_zzzzz_yyyzzz = pbuffer.data(idx_ovl_hi + 584);

    auto ts_zzzzz_yyzzzz = pbuffer.data(idx_ovl_hi + 585);

    auto ts_zzzzz_yzzzzz = pbuffer.data(idx_ovl_hi + 586);

    auto ts_zzzzz_zzzzzz = pbuffer.data(idx_ovl_hi + 587);

    // Set up 0-28 components of targeted buffer : HI

    auto tk_xxxxx_xxxxxx = pbuffer.data(idx_kin_hi);

    auto tk_xxxxx_xxxxxy = pbuffer.data(idx_kin_hi + 1);

    auto tk_xxxxx_xxxxxz = pbuffer.data(idx_kin_hi + 2);

    auto tk_xxxxx_xxxxyy = pbuffer.data(idx_kin_hi + 3);

    auto tk_xxxxx_xxxxyz = pbuffer.data(idx_kin_hi + 4);

    auto tk_xxxxx_xxxxzz = pbuffer.data(idx_kin_hi + 5);

    auto tk_xxxxx_xxxyyy = pbuffer.data(idx_kin_hi + 6);

    auto tk_xxxxx_xxxyyz = pbuffer.data(idx_kin_hi + 7);

    auto tk_xxxxx_xxxyzz = pbuffer.data(idx_kin_hi + 8);

    auto tk_xxxxx_xxxzzz = pbuffer.data(idx_kin_hi + 9);

    auto tk_xxxxx_xxyyyy = pbuffer.data(idx_kin_hi + 10);

    auto tk_xxxxx_xxyyyz = pbuffer.data(idx_kin_hi + 11);

    auto tk_xxxxx_xxyyzz = pbuffer.data(idx_kin_hi + 12);

    auto tk_xxxxx_xxyzzz = pbuffer.data(idx_kin_hi + 13);

    auto tk_xxxxx_xxzzzz = pbuffer.data(idx_kin_hi + 14);

    auto tk_xxxxx_xyyyyy = pbuffer.data(idx_kin_hi + 15);

    auto tk_xxxxx_xyyyyz = pbuffer.data(idx_kin_hi + 16);

    auto tk_xxxxx_xyyyzz = pbuffer.data(idx_kin_hi + 17);

    auto tk_xxxxx_xyyzzz = pbuffer.data(idx_kin_hi + 18);

    auto tk_xxxxx_xyzzzz = pbuffer.data(idx_kin_hi + 19);

    auto tk_xxxxx_xzzzzz = pbuffer.data(idx_kin_hi + 20);

    auto tk_xxxxx_yyyyyy = pbuffer.data(idx_kin_hi + 21);

    auto tk_xxxxx_yyyyyz = pbuffer.data(idx_kin_hi + 22);

    auto tk_xxxxx_yyyyzz = pbuffer.data(idx_kin_hi + 23);

    auto tk_xxxxx_yyyzzz = pbuffer.data(idx_kin_hi + 24);

    auto tk_xxxxx_yyzzzz = pbuffer.data(idx_kin_hi + 25);

    auto tk_xxxxx_yzzzzz = pbuffer.data(idx_kin_hi + 26);

    auto tk_xxxxx_zzzzzz = pbuffer.data(idx_kin_hi + 27);

#pragma omp simd aligned(pa_x,                \
                             tk_xxx_xxxxxx,   \
                             tk_xxx_xxxxxy,   \
                             tk_xxx_xxxxxz,   \
                             tk_xxx_xxxxyy,   \
                             tk_xxx_xxxxyz,   \
                             tk_xxx_xxxxzz,   \
                             tk_xxx_xxxyyy,   \
                             tk_xxx_xxxyyz,   \
                             tk_xxx_xxxyzz,   \
                             tk_xxx_xxxzzz,   \
                             tk_xxx_xxyyyy,   \
                             tk_xxx_xxyyyz,   \
                             tk_xxx_xxyyzz,   \
                             tk_xxx_xxyzzz,   \
                             tk_xxx_xxzzzz,   \
                             tk_xxx_xyyyyy,   \
                             tk_xxx_xyyyyz,   \
                             tk_xxx_xyyyzz,   \
                             tk_xxx_xyyzzz,   \
                             tk_xxx_xyzzzz,   \
                             tk_xxx_xzzzzz,   \
                             tk_xxx_yyyyyy,   \
                             tk_xxx_yyyyyz,   \
                             tk_xxx_yyyyzz,   \
                             tk_xxx_yyyzzz,   \
                             tk_xxx_yyzzzz,   \
                             tk_xxx_yzzzzz,   \
                             tk_xxx_zzzzzz,   \
                             tk_xxxx_xxxxx,   \
                             tk_xxxx_xxxxxx,  \
                             tk_xxxx_xxxxxy,  \
                             tk_xxxx_xxxxxz,  \
                             tk_xxxx_xxxxy,   \
                             tk_xxxx_xxxxyy,  \
                             tk_xxxx_xxxxyz,  \
                             tk_xxxx_xxxxz,   \
                             tk_xxxx_xxxxzz,  \
                             tk_xxxx_xxxyy,   \
                             tk_xxxx_xxxyyy,  \
                             tk_xxxx_xxxyyz,  \
                             tk_xxxx_xxxyz,   \
                             tk_xxxx_xxxyzz,  \
                             tk_xxxx_xxxzz,   \
                             tk_xxxx_xxxzzz,  \
                             tk_xxxx_xxyyy,   \
                             tk_xxxx_xxyyyy,  \
                             tk_xxxx_xxyyyz,  \
                             tk_xxxx_xxyyz,   \
                             tk_xxxx_xxyyzz,  \
                             tk_xxxx_xxyzz,   \
                             tk_xxxx_xxyzzz,  \
                             tk_xxxx_xxzzz,   \
                             tk_xxxx_xxzzzz,  \
                             tk_xxxx_xyyyy,   \
                             tk_xxxx_xyyyyy,  \
                             tk_xxxx_xyyyyz,  \
                             tk_xxxx_xyyyz,   \
                             tk_xxxx_xyyyzz,  \
                             tk_xxxx_xyyzz,   \
                             tk_xxxx_xyyzzz,  \
                             tk_xxxx_xyzzz,   \
                             tk_xxxx_xyzzzz,  \
                             tk_xxxx_xzzzz,   \
                             tk_xxxx_xzzzzz,  \
                             tk_xxxx_yyyyy,   \
                             tk_xxxx_yyyyyy,  \
                             tk_xxxx_yyyyyz,  \
                             tk_xxxx_yyyyz,   \
                             tk_xxxx_yyyyzz,  \
                             tk_xxxx_yyyzz,   \
                             tk_xxxx_yyyzzz,  \
                             tk_xxxx_yyzzz,   \
                             tk_xxxx_yyzzzz,  \
                             tk_xxxx_yzzzz,   \
                             tk_xxxx_yzzzzz,  \
                             tk_xxxx_zzzzz,   \
                             tk_xxxx_zzzzzz,  \
                             tk_xxxxx_xxxxxx, \
                             tk_xxxxx_xxxxxy, \
                             tk_xxxxx_xxxxxz, \
                             tk_xxxxx_xxxxyy, \
                             tk_xxxxx_xxxxyz, \
                             tk_xxxxx_xxxxzz, \
                             tk_xxxxx_xxxyyy, \
                             tk_xxxxx_xxxyyz, \
                             tk_xxxxx_xxxyzz, \
                             tk_xxxxx_xxxzzz, \
                             tk_xxxxx_xxyyyy, \
                             tk_xxxxx_xxyyyz, \
                             tk_xxxxx_xxyyzz, \
                             tk_xxxxx_xxyzzz, \
                             tk_xxxxx_xxzzzz, \
                             tk_xxxxx_xyyyyy, \
                             tk_xxxxx_xyyyyz, \
                             tk_xxxxx_xyyyzz, \
                             tk_xxxxx_xyyzzz, \
                             tk_xxxxx_xyzzzz, \
                             tk_xxxxx_xzzzzz, \
                             tk_xxxxx_yyyyyy, \
                             tk_xxxxx_yyyyyz, \
                             tk_xxxxx_yyyyzz, \
                             tk_xxxxx_yyyzzz, \
                             tk_xxxxx_yyzzzz, \
                             tk_xxxxx_yzzzzz, \
                             tk_xxxxx_zzzzzz, \
                             ts_xxx_xxxxxx,   \
                             ts_xxx_xxxxxy,   \
                             ts_xxx_xxxxxz,   \
                             ts_xxx_xxxxyy,   \
                             ts_xxx_xxxxyz,   \
                             ts_xxx_xxxxzz,   \
                             ts_xxx_xxxyyy,   \
                             ts_xxx_xxxyyz,   \
                             ts_xxx_xxxyzz,   \
                             ts_xxx_xxxzzz,   \
                             ts_xxx_xxyyyy,   \
                             ts_xxx_xxyyyz,   \
                             ts_xxx_xxyyzz,   \
                             ts_xxx_xxyzzz,   \
                             ts_xxx_xxzzzz,   \
                             ts_xxx_xyyyyy,   \
                             ts_xxx_xyyyyz,   \
                             ts_xxx_xyyyzz,   \
                             ts_xxx_xyyzzz,   \
                             ts_xxx_xyzzzz,   \
                             ts_xxx_xzzzzz,   \
                             ts_xxx_yyyyyy,   \
                             ts_xxx_yyyyyz,   \
                             ts_xxx_yyyyzz,   \
                             ts_xxx_yyyzzz,   \
                             ts_xxx_yyzzzz,   \
                             ts_xxx_yzzzzz,   \
                             ts_xxx_zzzzzz,   \
                             ts_xxxxx_xxxxxx, \
                             ts_xxxxx_xxxxxy, \
                             ts_xxxxx_xxxxxz, \
                             ts_xxxxx_xxxxyy, \
                             ts_xxxxx_xxxxyz, \
                             ts_xxxxx_xxxxzz, \
                             ts_xxxxx_xxxyyy, \
                             ts_xxxxx_xxxyyz, \
                             ts_xxxxx_xxxyzz, \
                             ts_xxxxx_xxxzzz, \
                             ts_xxxxx_xxyyyy, \
                             ts_xxxxx_xxyyyz, \
                             ts_xxxxx_xxyyzz, \
                             ts_xxxxx_xxyzzz, \
                             ts_xxxxx_xxzzzz, \
                             ts_xxxxx_xyyyyy, \
                             ts_xxxxx_xyyyyz, \
                             ts_xxxxx_xyyyzz, \
                             ts_xxxxx_xyyzzz, \
                             ts_xxxxx_xyzzzz, \
                             ts_xxxxx_xzzzzz, \
                             ts_xxxxx_yyyyyy, \
                             ts_xxxxx_yyyyyz, \
                             ts_xxxxx_yyyyzz, \
                             ts_xxxxx_yyyzzz, \
                             ts_xxxxx_yyzzzz, \
                             ts_xxxxx_yzzzzz, \
                             ts_xxxxx_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxx_xxxxxx[i] = -8.0 * ts_xxx_xxxxxx[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxxxx[i] * fe_0 + 6.0 * tk_xxxx_xxxxx[i] * fe_0 +
                             tk_xxxx_xxxxxx[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxxxx[i] * fz_0;

        tk_xxxxx_xxxxxy[i] = -8.0 * ts_xxx_xxxxxy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxxxy[i] * fe_0 + 5.0 * tk_xxxx_xxxxy[i] * fe_0 +
                             tk_xxxx_xxxxxy[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxxxy[i] * fz_0;

        tk_xxxxx_xxxxxz[i] = -8.0 * ts_xxx_xxxxxz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxxxz[i] * fe_0 + 5.0 * tk_xxxx_xxxxz[i] * fe_0 +
                             tk_xxxx_xxxxxz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxxxz[i] * fz_0;

        tk_xxxxx_xxxxyy[i] = -8.0 * ts_xxx_xxxxyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxxyy[i] * fe_0 + 4.0 * tk_xxxx_xxxyy[i] * fe_0 +
                             tk_xxxx_xxxxyy[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxxyy[i] * fz_0;

        tk_xxxxx_xxxxyz[i] = -8.0 * ts_xxx_xxxxyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxxyz[i] * fe_0 + 4.0 * tk_xxxx_xxxyz[i] * fe_0 +
                             tk_xxxx_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxxyz[i] * fz_0;

        tk_xxxxx_xxxxzz[i] = -8.0 * ts_xxx_xxxxzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxxzz[i] * fe_0 + 4.0 * tk_xxxx_xxxzz[i] * fe_0 +
                             tk_xxxx_xxxxzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxxzz[i] * fz_0;

        tk_xxxxx_xxxyyy[i] = -8.0 * ts_xxx_xxxyyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxyyy[i] * fe_0 + 3.0 * tk_xxxx_xxyyy[i] * fe_0 +
                             tk_xxxx_xxxyyy[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxyyy[i] * fz_0;

        tk_xxxxx_xxxyyz[i] = -8.0 * ts_xxx_xxxyyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxyyz[i] * fe_0 + 3.0 * tk_xxxx_xxyyz[i] * fe_0 +
                             tk_xxxx_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxyyz[i] * fz_0;

        tk_xxxxx_xxxyzz[i] = -8.0 * ts_xxx_xxxyzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxyzz[i] * fe_0 + 3.0 * tk_xxxx_xxyzz[i] * fe_0 +
                             tk_xxxx_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxyzz[i] * fz_0;

        tk_xxxxx_xxxzzz[i] = -8.0 * ts_xxx_xxxzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxxzzz[i] * fe_0 + 3.0 * tk_xxxx_xxzzz[i] * fe_0 +
                             tk_xxxx_xxxzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxxzzz[i] * fz_0;

        tk_xxxxx_xxyyyy[i] = -8.0 * ts_xxx_xxyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxyyyy[i] * fe_0 + 2.0 * tk_xxxx_xyyyy[i] * fe_0 +
                             tk_xxxx_xxyyyy[i] * pa_x[i] + 2.0 * ts_xxxxx_xxyyyy[i] * fz_0;

        tk_xxxxx_xxyyyz[i] = -8.0 * ts_xxx_xxyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxyyyz[i] * fe_0 + 2.0 * tk_xxxx_xyyyz[i] * fe_0 +
                             tk_xxxx_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxyyyz[i] * fz_0;

        tk_xxxxx_xxyyzz[i] = -8.0 * ts_xxx_xxyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxyyzz[i] * fe_0 + 2.0 * tk_xxxx_xyyzz[i] * fe_0 +
                             tk_xxxx_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxyyzz[i] * fz_0;

        tk_xxxxx_xxyzzz[i] = -8.0 * ts_xxx_xxyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxyzzz[i] * fe_0 + 2.0 * tk_xxxx_xyzzz[i] * fe_0 +
                             tk_xxxx_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxyzzz[i] * fz_0;

        tk_xxxxx_xxzzzz[i] = -8.0 * ts_xxx_xxzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xxzzzz[i] * fe_0 + 2.0 * tk_xxxx_xzzzz[i] * fe_0 +
                             tk_xxxx_xxzzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xxzzzz[i] * fz_0;

        tk_xxxxx_xyyyyy[i] = -8.0 * ts_xxx_xyyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyyyyy[i] * fe_0 + tk_xxxx_yyyyy[i] * fe_0 +
                             tk_xxxx_xyyyyy[i] * pa_x[i] + 2.0 * ts_xxxxx_xyyyyy[i] * fz_0;

        tk_xxxxx_xyyyyz[i] = -8.0 * ts_xxx_xyyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyyyyz[i] * fe_0 + tk_xxxx_yyyyz[i] * fe_0 +
                             tk_xxxx_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxxxx_xyyyyz[i] * fz_0;

        tk_xxxxx_xyyyzz[i] = -8.0 * ts_xxx_xyyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyyyzz[i] * fe_0 + tk_xxxx_yyyzz[i] * fe_0 +
                             tk_xxxx_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xyyyzz[i] * fz_0;

        tk_xxxxx_xyyzzz[i] = -8.0 * ts_xxx_xyyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyyzzz[i] * fe_0 + tk_xxxx_yyzzz[i] * fe_0 +
                             tk_xxxx_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xyyzzz[i] * fz_0;

        tk_xxxxx_xyzzzz[i] = -8.0 * ts_xxx_xyzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xyzzzz[i] * fe_0 + tk_xxxx_yzzzz[i] * fe_0 +
                             tk_xxxx_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xyzzzz[i] * fz_0;

        tk_xxxxx_xzzzzz[i] = -8.0 * ts_xxx_xzzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xzzzzz[i] * fe_0 + tk_xxxx_zzzzz[i] * fe_0 +
                             tk_xxxx_xzzzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_xzzzzz[i] * fz_0;

        tk_xxxxx_yyyyyy[i] =
            -8.0 * ts_xxx_yyyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyyyyy[i] * fe_0 + tk_xxxx_yyyyyy[i] * pa_x[i] + 2.0 * ts_xxxxx_yyyyyy[i] * fz_0;

        tk_xxxxx_yyyyyz[i] =
            -8.0 * ts_xxx_yyyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyyyyz[i] * fe_0 + tk_xxxx_yyyyyz[i] * pa_x[i] + 2.0 * ts_xxxxx_yyyyyz[i] * fz_0;

        tk_xxxxx_yyyyzz[i] =
            -8.0 * ts_xxx_yyyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyyyzz[i] * fe_0 + tk_xxxx_yyyyzz[i] * pa_x[i] + 2.0 * ts_xxxxx_yyyyzz[i] * fz_0;

        tk_xxxxx_yyyzzz[i] =
            -8.0 * ts_xxx_yyyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyyzzz[i] * fe_0 + tk_xxxx_yyyzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_yyyzzz[i] * fz_0;

        tk_xxxxx_yyzzzz[i] =
            -8.0 * ts_xxx_yyzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yyzzzz[i] * fe_0 + tk_xxxx_yyzzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_yyzzzz[i] * fz_0;

        tk_xxxxx_yzzzzz[i] =
            -8.0 * ts_xxx_yzzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yzzzzz[i] * fe_0 + tk_xxxx_yzzzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_yzzzzz[i] * fz_0;

        tk_xxxxx_zzzzzz[i] =
            -8.0 * ts_xxx_zzzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_zzzzzz[i] * fe_0 + tk_xxxx_zzzzzz[i] * pa_x[i] + 2.0 * ts_xxxxx_zzzzzz[i] * fz_0;
    }

    // Set up 28-56 components of targeted buffer : HI

    auto tk_xxxxy_xxxxxx = pbuffer.data(idx_kin_hi + 28);

    auto tk_xxxxy_xxxxxy = pbuffer.data(idx_kin_hi + 29);

    auto tk_xxxxy_xxxxxz = pbuffer.data(idx_kin_hi + 30);

    auto tk_xxxxy_xxxxyy = pbuffer.data(idx_kin_hi + 31);

    auto tk_xxxxy_xxxxyz = pbuffer.data(idx_kin_hi + 32);

    auto tk_xxxxy_xxxxzz = pbuffer.data(idx_kin_hi + 33);

    auto tk_xxxxy_xxxyyy = pbuffer.data(idx_kin_hi + 34);

    auto tk_xxxxy_xxxyyz = pbuffer.data(idx_kin_hi + 35);

    auto tk_xxxxy_xxxyzz = pbuffer.data(idx_kin_hi + 36);

    auto tk_xxxxy_xxxzzz = pbuffer.data(idx_kin_hi + 37);

    auto tk_xxxxy_xxyyyy = pbuffer.data(idx_kin_hi + 38);

    auto tk_xxxxy_xxyyyz = pbuffer.data(idx_kin_hi + 39);

    auto tk_xxxxy_xxyyzz = pbuffer.data(idx_kin_hi + 40);

    auto tk_xxxxy_xxyzzz = pbuffer.data(idx_kin_hi + 41);

    auto tk_xxxxy_xxzzzz = pbuffer.data(idx_kin_hi + 42);

    auto tk_xxxxy_xyyyyy = pbuffer.data(idx_kin_hi + 43);

    auto tk_xxxxy_xyyyyz = pbuffer.data(idx_kin_hi + 44);

    auto tk_xxxxy_xyyyzz = pbuffer.data(idx_kin_hi + 45);

    auto tk_xxxxy_xyyzzz = pbuffer.data(idx_kin_hi + 46);

    auto tk_xxxxy_xyzzzz = pbuffer.data(idx_kin_hi + 47);

    auto tk_xxxxy_xzzzzz = pbuffer.data(idx_kin_hi + 48);

    auto tk_xxxxy_yyyyyy = pbuffer.data(idx_kin_hi + 49);

    auto tk_xxxxy_yyyyyz = pbuffer.data(idx_kin_hi + 50);

    auto tk_xxxxy_yyyyzz = pbuffer.data(idx_kin_hi + 51);

    auto tk_xxxxy_yyyzzz = pbuffer.data(idx_kin_hi + 52);

    auto tk_xxxxy_yyzzzz = pbuffer.data(idx_kin_hi + 53);

    auto tk_xxxxy_yzzzzz = pbuffer.data(idx_kin_hi + 54);

    auto tk_xxxxy_zzzzzz = pbuffer.data(idx_kin_hi + 55);

#pragma omp simd aligned(pa_y,                \
                             tk_xxxx_xxxxx,   \
                             tk_xxxx_xxxxxx,  \
                             tk_xxxx_xxxxxy,  \
                             tk_xxxx_xxxxxz,  \
                             tk_xxxx_xxxxy,   \
                             tk_xxxx_xxxxyy,  \
                             tk_xxxx_xxxxyz,  \
                             tk_xxxx_xxxxz,   \
                             tk_xxxx_xxxxzz,  \
                             tk_xxxx_xxxyy,   \
                             tk_xxxx_xxxyyy,  \
                             tk_xxxx_xxxyyz,  \
                             tk_xxxx_xxxyz,   \
                             tk_xxxx_xxxyzz,  \
                             tk_xxxx_xxxzz,   \
                             tk_xxxx_xxxzzz,  \
                             tk_xxxx_xxyyy,   \
                             tk_xxxx_xxyyyy,  \
                             tk_xxxx_xxyyyz,  \
                             tk_xxxx_xxyyz,   \
                             tk_xxxx_xxyyzz,  \
                             tk_xxxx_xxyzz,   \
                             tk_xxxx_xxyzzz,  \
                             tk_xxxx_xxzzz,   \
                             tk_xxxx_xxzzzz,  \
                             tk_xxxx_xyyyy,   \
                             tk_xxxx_xyyyyy,  \
                             tk_xxxx_xyyyyz,  \
                             tk_xxxx_xyyyz,   \
                             tk_xxxx_xyyyzz,  \
                             tk_xxxx_xyyzz,   \
                             tk_xxxx_xyyzzz,  \
                             tk_xxxx_xyzzz,   \
                             tk_xxxx_xyzzzz,  \
                             tk_xxxx_xzzzz,   \
                             tk_xxxx_xzzzzz,  \
                             tk_xxxx_yyyyy,   \
                             tk_xxxx_yyyyyy,  \
                             tk_xxxx_yyyyyz,  \
                             tk_xxxx_yyyyz,   \
                             tk_xxxx_yyyyzz,  \
                             tk_xxxx_yyyzz,   \
                             tk_xxxx_yyyzzz,  \
                             tk_xxxx_yyzzz,   \
                             tk_xxxx_yyzzzz,  \
                             tk_xxxx_yzzzz,   \
                             tk_xxxx_yzzzzz,  \
                             tk_xxxx_zzzzz,   \
                             tk_xxxx_zzzzzz,  \
                             tk_xxxxy_xxxxxx, \
                             tk_xxxxy_xxxxxy, \
                             tk_xxxxy_xxxxxz, \
                             tk_xxxxy_xxxxyy, \
                             tk_xxxxy_xxxxyz, \
                             tk_xxxxy_xxxxzz, \
                             tk_xxxxy_xxxyyy, \
                             tk_xxxxy_xxxyyz, \
                             tk_xxxxy_xxxyzz, \
                             tk_xxxxy_xxxzzz, \
                             tk_xxxxy_xxyyyy, \
                             tk_xxxxy_xxyyyz, \
                             tk_xxxxy_xxyyzz, \
                             tk_xxxxy_xxyzzz, \
                             tk_xxxxy_xxzzzz, \
                             tk_xxxxy_xyyyyy, \
                             tk_xxxxy_xyyyyz, \
                             tk_xxxxy_xyyyzz, \
                             tk_xxxxy_xyyzzz, \
                             tk_xxxxy_xyzzzz, \
                             tk_xxxxy_xzzzzz, \
                             tk_xxxxy_yyyyyy, \
                             tk_xxxxy_yyyyyz, \
                             tk_xxxxy_yyyyzz, \
                             tk_xxxxy_yyyzzz, \
                             tk_xxxxy_yyzzzz, \
                             tk_xxxxy_yzzzzz, \
                             tk_xxxxy_zzzzzz, \
                             ts_xxxxy_xxxxxx, \
                             ts_xxxxy_xxxxxy, \
                             ts_xxxxy_xxxxxz, \
                             ts_xxxxy_xxxxyy, \
                             ts_xxxxy_xxxxyz, \
                             ts_xxxxy_xxxxzz, \
                             ts_xxxxy_xxxyyy, \
                             ts_xxxxy_xxxyyz, \
                             ts_xxxxy_xxxyzz, \
                             ts_xxxxy_xxxzzz, \
                             ts_xxxxy_xxyyyy, \
                             ts_xxxxy_xxyyyz, \
                             ts_xxxxy_xxyyzz, \
                             ts_xxxxy_xxyzzz, \
                             ts_xxxxy_xxzzzz, \
                             ts_xxxxy_xyyyyy, \
                             ts_xxxxy_xyyyyz, \
                             ts_xxxxy_xyyyzz, \
                             ts_xxxxy_xyyzzz, \
                             ts_xxxxy_xyzzzz, \
                             ts_xxxxy_xzzzzz, \
                             ts_xxxxy_yyyyyy, \
                             ts_xxxxy_yyyyyz, \
                             ts_xxxxy_yyyyzz, \
                             ts_xxxxy_yyyzzz, \
                             ts_xxxxy_yyzzzz, \
                             ts_xxxxy_yzzzzz, \
                             ts_xxxxy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxy_xxxxxx[i] = tk_xxxx_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxxxx[i] * fz_0;

        tk_xxxxy_xxxxxy[i] = tk_xxxx_xxxxx[i] * fe_0 + tk_xxxx_xxxxxy[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxxxy[i] * fz_0;

        tk_xxxxy_xxxxxz[i] = tk_xxxx_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxxxz[i] * fz_0;

        tk_xxxxy_xxxxyy[i] = 2.0 * tk_xxxx_xxxxy[i] * fe_0 + tk_xxxx_xxxxyy[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxxyy[i] * fz_0;

        tk_xxxxy_xxxxyz[i] = tk_xxxx_xxxxz[i] * fe_0 + tk_xxxx_xxxxyz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxxyz[i] * fz_0;

        tk_xxxxy_xxxxzz[i] = tk_xxxx_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxxzz[i] * fz_0;

        tk_xxxxy_xxxyyy[i] = 3.0 * tk_xxxx_xxxyy[i] * fe_0 + tk_xxxx_xxxyyy[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxyyy[i] * fz_0;

        tk_xxxxy_xxxyyz[i] = 2.0 * tk_xxxx_xxxyz[i] * fe_0 + tk_xxxx_xxxyyz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxyyz[i] * fz_0;

        tk_xxxxy_xxxyzz[i] = tk_xxxx_xxxzz[i] * fe_0 + tk_xxxx_xxxyzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxyzz[i] * fz_0;

        tk_xxxxy_xxxzzz[i] = tk_xxxx_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxxzzz[i] * fz_0;

        tk_xxxxy_xxyyyy[i] = 4.0 * tk_xxxx_xxyyy[i] * fe_0 + tk_xxxx_xxyyyy[i] * pa_y[i] + 2.0 * ts_xxxxy_xxyyyy[i] * fz_0;

        tk_xxxxy_xxyyyz[i] = 3.0 * tk_xxxx_xxyyz[i] * fe_0 + tk_xxxx_xxyyyz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxyyyz[i] * fz_0;

        tk_xxxxy_xxyyzz[i] = 2.0 * tk_xxxx_xxyzz[i] * fe_0 + tk_xxxx_xxyyzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxyyzz[i] * fz_0;

        tk_xxxxy_xxyzzz[i] = tk_xxxx_xxzzz[i] * fe_0 + tk_xxxx_xxyzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxyzzz[i] * fz_0;

        tk_xxxxy_xxzzzz[i] = tk_xxxx_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xxzzzz[i] * fz_0;

        tk_xxxxy_xyyyyy[i] = 5.0 * tk_xxxx_xyyyy[i] * fe_0 + tk_xxxx_xyyyyy[i] * pa_y[i] + 2.0 * ts_xxxxy_xyyyyy[i] * fz_0;

        tk_xxxxy_xyyyyz[i] = 4.0 * tk_xxxx_xyyyz[i] * fe_0 + tk_xxxx_xyyyyz[i] * pa_y[i] + 2.0 * ts_xxxxy_xyyyyz[i] * fz_0;

        tk_xxxxy_xyyyzz[i] = 3.0 * tk_xxxx_xyyzz[i] * fe_0 + tk_xxxx_xyyyzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xyyyzz[i] * fz_0;

        tk_xxxxy_xyyzzz[i] = 2.0 * tk_xxxx_xyzzz[i] * fe_0 + tk_xxxx_xyyzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xyyzzz[i] * fz_0;

        tk_xxxxy_xyzzzz[i] = tk_xxxx_xzzzz[i] * fe_0 + tk_xxxx_xyzzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xyzzzz[i] * fz_0;

        tk_xxxxy_xzzzzz[i] = tk_xxxx_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_xzzzzz[i] * fz_0;

        tk_xxxxy_yyyyyy[i] = 6.0 * tk_xxxx_yyyyy[i] * fe_0 + tk_xxxx_yyyyyy[i] * pa_y[i] + 2.0 * ts_xxxxy_yyyyyy[i] * fz_0;

        tk_xxxxy_yyyyyz[i] = 5.0 * tk_xxxx_yyyyz[i] * fe_0 + tk_xxxx_yyyyyz[i] * pa_y[i] + 2.0 * ts_xxxxy_yyyyyz[i] * fz_0;

        tk_xxxxy_yyyyzz[i] = 4.0 * tk_xxxx_yyyzz[i] * fe_0 + tk_xxxx_yyyyzz[i] * pa_y[i] + 2.0 * ts_xxxxy_yyyyzz[i] * fz_0;

        tk_xxxxy_yyyzzz[i] = 3.0 * tk_xxxx_yyzzz[i] * fe_0 + tk_xxxx_yyyzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_yyyzzz[i] * fz_0;

        tk_xxxxy_yyzzzz[i] = 2.0 * tk_xxxx_yzzzz[i] * fe_0 + tk_xxxx_yyzzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_yyzzzz[i] * fz_0;

        tk_xxxxy_yzzzzz[i] = tk_xxxx_zzzzz[i] * fe_0 + tk_xxxx_yzzzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_yzzzzz[i] * fz_0;

        tk_xxxxy_zzzzzz[i] = tk_xxxx_zzzzzz[i] * pa_y[i] + 2.0 * ts_xxxxy_zzzzzz[i] * fz_0;
    }

    // Set up 56-84 components of targeted buffer : HI

    auto tk_xxxxz_xxxxxx = pbuffer.data(idx_kin_hi + 56);

    auto tk_xxxxz_xxxxxy = pbuffer.data(idx_kin_hi + 57);

    auto tk_xxxxz_xxxxxz = pbuffer.data(idx_kin_hi + 58);

    auto tk_xxxxz_xxxxyy = pbuffer.data(idx_kin_hi + 59);

    auto tk_xxxxz_xxxxyz = pbuffer.data(idx_kin_hi + 60);

    auto tk_xxxxz_xxxxzz = pbuffer.data(idx_kin_hi + 61);

    auto tk_xxxxz_xxxyyy = pbuffer.data(idx_kin_hi + 62);

    auto tk_xxxxz_xxxyyz = pbuffer.data(idx_kin_hi + 63);

    auto tk_xxxxz_xxxyzz = pbuffer.data(idx_kin_hi + 64);

    auto tk_xxxxz_xxxzzz = pbuffer.data(idx_kin_hi + 65);

    auto tk_xxxxz_xxyyyy = pbuffer.data(idx_kin_hi + 66);

    auto tk_xxxxz_xxyyyz = pbuffer.data(idx_kin_hi + 67);

    auto tk_xxxxz_xxyyzz = pbuffer.data(idx_kin_hi + 68);

    auto tk_xxxxz_xxyzzz = pbuffer.data(idx_kin_hi + 69);

    auto tk_xxxxz_xxzzzz = pbuffer.data(idx_kin_hi + 70);

    auto tk_xxxxz_xyyyyy = pbuffer.data(idx_kin_hi + 71);

    auto tk_xxxxz_xyyyyz = pbuffer.data(idx_kin_hi + 72);

    auto tk_xxxxz_xyyyzz = pbuffer.data(idx_kin_hi + 73);

    auto tk_xxxxz_xyyzzz = pbuffer.data(idx_kin_hi + 74);

    auto tk_xxxxz_xyzzzz = pbuffer.data(idx_kin_hi + 75);

    auto tk_xxxxz_xzzzzz = pbuffer.data(idx_kin_hi + 76);

    auto tk_xxxxz_yyyyyy = pbuffer.data(idx_kin_hi + 77);

    auto tk_xxxxz_yyyyyz = pbuffer.data(idx_kin_hi + 78);

    auto tk_xxxxz_yyyyzz = pbuffer.data(idx_kin_hi + 79);

    auto tk_xxxxz_yyyzzz = pbuffer.data(idx_kin_hi + 80);

    auto tk_xxxxz_yyzzzz = pbuffer.data(idx_kin_hi + 81);

    auto tk_xxxxz_yzzzzz = pbuffer.data(idx_kin_hi + 82);

    auto tk_xxxxz_zzzzzz = pbuffer.data(idx_kin_hi + 83);

#pragma omp simd aligned(pa_z,                \
                             tk_xxxx_xxxxx,   \
                             tk_xxxx_xxxxxx,  \
                             tk_xxxx_xxxxxy,  \
                             tk_xxxx_xxxxxz,  \
                             tk_xxxx_xxxxy,   \
                             tk_xxxx_xxxxyy,  \
                             tk_xxxx_xxxxyz,  \
                             tk_xxxx_xxxxz,   \
                             tk_xxxx_xxxxzz,  \
                             tk_xxxx_xxxyy,   \
                             tk_xxxx_xxxyyy,  \
                             tk_xxxx_xxxyyz,  \
                             tk_xxxx_xxxyz,   \
                             tk_xxxx_xxxyzz,  \
                             tk_xxxx_xxxzz,   \
                             tk_xxxx_xxxzzz,  \
                             tk_xxxx_xxyyy,   \
                             tk_xxxx_xxyyyy,  \
                             tk_xxxx_xxyyyz,  \
                             tk_xxxx_xxyyz,   \
                             tk_xxxx_xxyyzz,  \
                             tk_xxxx_xxyzz,   \
                             tk_xxxx_xxyzzz,  \
                             tk_xxxx_xxzzz,   \
                             tk_xxxx_xxzzzz,  \
                             tk_xxxx_xyyyy,   \
                             tk_xxxx_xyyyyy,  \
                             tk_xxxx_xyyyyz,  \
                             tk_xxxx_xyyyz,   \
                             tk_xxxx_xyyyzz,  \
                             tk_xxxx_xyyzz,   \
                             tk_xxxx_xyyzzz,  \
                             tk_xxxx_xyzzz,   \
                             tk_xxxx_xyzzzz,  \
                             tk_xxxx_xzzzz,   \
                             tk_xxxx_xzzzzz,  \
                             tk_xxxx_yyyyy,   \
                             tk_xxxx_yyyyyy,  \
                             tk_xxxx_yyyyyz,  \
                             tk_xxxx_yyyyz,   \
                             tk_xxxx_yyyyzz,  \
                             tk_xxxx_yyyzz,   \
                             tk_xxxx_yyyzzz,  \
                             tk_xxxx_yyzzz,   \
                             tk_xxxx_yyzzzz,  \
                             tk_xxxx_yzzzz,   \
                             tk_xxxx_yzzzzz,  \
                             tk_xxxx_zzzzz,   \
                             tk_xxxx_zzzzzz,  \
                             tk_xxxxz_xxxxxx, \
                             tk_xxxxz_xxxxxy, \
                             tk_xxxxz_xxxxxz, \
                             tk_xxxxz_xxxxyy, \
                             tk_xxxxz_xxxxyz, \
                             tk_xxxxz_xxxxzz, \
                             tk_xxxxz_xxxyyy, \
                             tk_xxxxz_xxxyyz, \
                             tk_xxxxz_xxxyzz, \
                             tk_xxxxz_xxxzzz, \
                             tk_xxxxz_xxyyyy, \
                             tk_xxxxz_xxyyyz, \
                             tk_xxxxz_xxyyzz, \
                             tk_xxxxz_xxyzzz, \
                             tk_xxxxz_xxzzzz, \
                             tk_xxxxz_xyyyyy, \
                             tk_xxxxz_xyyyyz, \
                             tk_xxxxz_xyyyzz, \
                             tk_xxxxz_xyyzzz, \
                             tk_xxxxz_xyzzzz, \
                             tk_xxxxz_xzzzzz, \
                             tk_xxxxz_yyyyyy, \
                             tk_xxxxz_yyyyyz, \
                             tk_xxxxz_yyyyzz, \
                             tk_xxxxz_yyyzzz, \
                             tk_xxxxz_yyzzzz, \
                             tk_xxxxz_yzzzzz, \
                             tk_xxxxz_zzzzzz, \
                             ts_xxxxz_xxxxxx, \
                             ts_xxxxz_xxxxxy, \
                             ts_xxxxz_xxxxxz, \
                             ts_xxxxz_xxxxyy, \
                             ts_xxxxz_xxxxyz, \
                             ts_xxxxz_xxxxzz, \
                             ts_xxxxz_xxxyyy, \
                             ts_xxxxz_xxxyyz, \
                             ts_xxxxz_xxxyzz, \
                             ts_xxxxz_xxxzzz, \
                             ts_xxxxz_xxyyyy, \
                             ts_xxxxz_xxyyyz, \
                             ts_xxxxz_xxyyzz, \
                             ts_xxxxz_xxyzzz, \
                             ts_xxxxz_xxzzzz, \
                             ts_xxxxz_xyyyyy, \
                             ts_xxxxz_xyyyyz, \
                             ts_xxxxz_xyyyzz, \
                             ts_xxxxz_xyyzzz, \
                             ts_xxxxz_xyzzzz, \
                             ts_xxxxz_xzzzzz, \
                             ts_xxxxz_yyyyyy, \
                             ts_xxxxz_yyyyyz, \
                             ts_xxxxz_yyyyzz, \
                             ts_xxxxz_yyyzzz, \
                             ts_xxxxz_yyzzzz, \
                             ts_xxxxz_yzzzzz, \
                             ts_xxxxz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxz_xxxxxx[i] = tk_xxxx_xxxxxx[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxxxx[i] * fz_0;

        tk_xxxxz_xxxxxy[i] = tk_xxxx_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxxxy[i] * fz_0;

        tk_xxxxz_xxxxxz[i] = tk_xxxx_xxxxx[i] * fe_0 + tk_xxxx_xxxxxz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxxxz[i] * fz_0;

        tk_xxxxz_xxxxyy[i] = tk_xxxx_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxxyy[i] * fz_0;

        tk_xxxxz_xxxxyz[i] = tk_xxxx_xxxxy[i] * fe_0 + tk_xxxx_xxxxyz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxxyz[i] * fz_0;

        tk_xxxxz_xxxxzz[i] = 2.0 * tk_xxxx_xxxxz[i] * fe_0 + tk_xxxx_xxxxzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxxzz[i] * fz_0;

        tk_xxxxz_xxxyyy[i] = tk_xxxx_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxyyy[i] * fz_0;

        tk_xxxxz_xxxyyz[i] = tk_xxxx_xxxyy[i] * fe_0 + tk_xxxx_xxxyyz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxyyz[i] * fz_0;

        tk_xxxxz_xxxyzz[i] = 2.0 * tk_xxxx_xxxyz[i] * fe_0 + tk_xxxx_xxxyzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxyzz[i] * fz_0;

        tk_xxxxz_xxxzzz[i] = 3.0 * tk_xxxx_xxxzz[i] * fe_0 + tk_xxxx_xxxzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxxzzz[i] * fz_0;

        tk_xxxxz_xxyyyy[i] = tk_xxxx_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxxxz_xxyyyy[i] * fz_0;

        tk_xxxxz_xxyyyz[i] = tk_xxxx_xxyyy[i] * fe_0 + tk_xxxx_xxyyyz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxyyyz[i] * fz_0;

        tk_xxxxz_xxyyzz[i] = 2.0 * tk_xxxx_xxyyz[i] * fe_0 + tk_xxxx_xxyyzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxyyzz[i] * fz_0;

        tk_xxxxz_xxyzzz[i] = 3.0 * tk_xxxx_xxyzz[i] * fe_0 + tk_xxxx_xxyzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxyzzz[i] * fz_0;

        tk_xxxxz_xxzzzz[i] = 4.0 * tk_xxxx_xxzzz[i] * fe_0 + tk_xxxx_xxzzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xxzzzz[i] * fz_0;

        tk_xxxxz_xyyyyy[i] = tk_xxxx_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxxxz_xyyyyy[i] * fz_0;

        tk_xxxxz_xyyyyz[i] = tk_xxxx_xyyyy[i] * fe_0 + tk_xxxx_xyyyyz[i] * pa_z[i] + 2.0 * ts_xxxxz_xyyyyz[i] * fz_0;

        tk_xxxxz_xyyyzz[i] = 2.0 * tk_xxxx_xyyyz[i] * fe_0 + tk_xxxx_xyyyzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xyyyzz[i] * fz_0;

        tk_xxxxz_xyyzzz[i] = 3.0 * tk_xxxx_xyyzz[i] * fe_0 + tk_xxxx_xyyzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xyyzzz[i] * fz_0;

        tk_xxxxz_xyzzzz[i] = 4.0 * tk_xxxx_xyzzz[i] * fe_0 + tk_xxxx_xyzzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xyzzzz[i] * fz_0;

        tk_xxxxz_xzzzzz[i] = 5.0 * tk_xxxx_xzzzz[i] * fe_0 + tk_xxxx_xzzzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_xzzzzz[i] * fz_0;

        tk_xxxxz_yyyyyy[i] = tk_xxxx_yyyyyy[i] * pa_z[i] + 2.0 * ts_xxxxz_yyyyyy[i] * fz_0;

        tk_xxxxz_yyyyyz[i] = tk_xxxx_yyyyy[i] * fe_0 + tk_xxxx_yyyyyz[i] * pa_z[i] + 2.0 * ts_xxxxz_yyyyyz[i] * fz_0;

        tk_xxxxz_yyyyzz[i] = 2.0 * tk_xxxx_yyyyz[i] * fe_0 + tk_xxxx_yyyyzz[i] * pa_z[i] + 2.0 * ts_xxxxz_yyyyzz[i] * fz_0;

        tk_xxxxz_yyyzzz[i] = 3.0 * tk_xxxx_yyyzz[i] * fe_0 + tk_xxxx_yyyzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_yyyzzz[i] * fz_0;

        tk_xxxxz_yyzzzz[i] = 4.0 * tk_xxxx_yyzzz[i] * fe_0 + tk_xxxx_yyzzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_yyzzzz[i] * fz_0;

        tk_xxxxz_yzzzzz[i] = 5.0 * tk_xxxx_yzzzz[i] * fe_0 + tk_xxxx_yzzzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_yzzzzz[i] * fz_0;

        tk_xxxxz_zzzzzz[i] = 6.0 * tk_xxxx_zzzzz[i] * fe_0 + tk_xxxx_zzzzzz[i] * pa_z[i] + 2.0 * ts_xxxxz_zzzzzz[i] * fz_0;
    }

    // Set up 84-112 components of targeted buffer : HI

    auto tk_xxxyy_xxxxxx = pbuffer.data(idx_kin_hi + 84);

    auto tk_xxxyy_xxxxxy = pbuffer.data(idx_kin_hi + 85);

    auto tk_xxxyy_xxxxxz = pbuffer.data(idx_kin_hi + 86);

    auto tk_xxxyy_xxxxyy = pbuffer.data(idx_kin_hi + 87);

    auto tk_xxxyy_xxxxyz = pbuffer.data(idx_kin_hi + 88);

    auto tk_xxxyy_xxxxzz = pbuffer.data(idx_kin_hi + 89);

    auto tk_xxxyy_xxxyyy = pbuffer.data(idx_kin_hi + 90);

    auto tk_xxxyy_xxxyyz = pbuffer.data(idx_kin_hi + 91);

    auto tk_xxxyy_xxxyzz = pbuffer.data(idx_kin_hi + 92);

    auto tk_xxxyy_xxxzzz = pbuffer.data(idx_kin_hi + 93);

    auto tk_xxxyy_xxyyyy = pbuffer.data(idx_kin_hi + 94);

    auto tk_xxxyy_xxyyyz = pbuffer.data(idx_kin_hi + 95);

    auto tk_xxxyy_xxyyzz = pbuffer.data(idx_kin_hi + 96);

    auto tk_xxxyy_xxyzzz = pbuffer.data(idx_kin_hi + 97);

    auto tk_xxxyy_xxzzzz = pbuffer.data(idx_kin_hi + 98);

    auto tk_xxxyy_xyyyyy = pbuffer.data(idx_kin_hi + 99);

    auto tk_xxxyy_xyyyyz = pbuffer.data(idx_kin_hi + 100);

    auto tk_xxxyy_xyyyzz = pbuffer.data(idx_kin_hi + 101);

    auto tk_xxxyy_xyyzzz = pbuffer.data(idx_kin_hi + 102);

    auto tk_xxxyy_xyzzzz = pbuffer.data(idx_kin_hi + 103);

    auto tk_xxxyy_xzzzzz = pbuffer.data(idx_kin_hi + 104);

    auto tk_xxxyy_yyyyyy = pbuffer.data(idx_kin_hi + 105);

    auto tk_xxxyy_yyyyyz = pbuffer.data(idx_kin_hi + 106);

    auto tk_xxxyy_yyyyzz = pbuffer.data(idx_kin_hi + 107);

    auto tk_xxxyy_yyyzzz = pbuffer.data(idx_kin_hi + 108);

    auto tk_xxxyy_yyzzzz = pbuffer.data(idx_kin_hi + 109);

    auto tk_xxxyy_yzzzzz = pbuffer.data(idx_kin_hi + 110);

    auto tk_xxxyy_zzzzzz = pbuffer.data(idx_kin_hi + 111);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tk_xxx_xxxxxx,   \
                             tk_xxx_xxxxxz,   \
                             tk_xxx_xxxxzz,   \
                             tk_xxx_xxxzzz,   \
                             tk_xxx_xxzzzz,   \
                             tk_xxx_xzzzzz,   \
                             tk_xxxy_xxxxxx,  \
                             tk_xxxy_xxxxxz,  \
                             tk_xxxy_xxxxzz,  \
                             tk_xxxy_xxxzzz,  \
                             tk_xxxy_xxzzzz,  \
                             tk_xxxy_xzzzzz,  \
                             tk_xxxyy_xxxxxx, \
                             tk_xxxyy_xxxxxy, \
                             tk_xxxyy_xxxxxz, \
                             tk_xxxyy_xxxxyy, \
                             tk_xxxyy_xxxxyz, \
                             tk_xxxyy_xxxxzz, \
                             tk_xxxyy_xxxyyy, \
                             tk_xxxyy_xxxyyz, \
                             tk_xxxyy_xxxyzz, \
                             tk_xxxyy_xxxzzz, \
                             tk_xxxyy_xxyyyy, \
                             tk_xxxyy_xxyyyz, \
                             tk_xxxyy_xxyyzz, \
                             tk_xxxyy_xxyzzz, \
                             tk_xxxyy_xxzzzz, \
                             tk_xxxyy_xyyyyy, \
                             tk_xxxyy_xyyyyz, \
                             tk_xxxyy_xyyyzz, \
                             tk_xxxyy_xyyzzz, \
                             tk_xxxyy_xyzzzz, \
                             tk_xxxyy_xzzzzz, \
                             tk_xxxyy_yyyyyy, \
                             tk_xxxyy_yyyyyz, \
                             tk_xxxyy_yyyyzz, \
                             tk_xxxyy_yyyzzz, \
                             tk_xxxyy_yyzzzz, \
                             tk_xxxyy_yzzzzz, \
                             tk_xxxyy_zzzzzz, \
                             tk_xxyy_xxxxxy,  \
                             tk_xxyy_xxxxy,   \
                             tk_xxyy_xxxxyy,  \
                             tk_xxyy_xxxxyz,  \
                             tk_xxyy_xxxyy,   \
                             tk_xxyy_xxxyyy,  \
                             tk_xxyy_xxxyyz,  \
                             tk_xxyy_xxxyz,   \
                             tk_xxyy_xxxyzz,  \
                             tk_xxyy_xxyyy,   \
                             tk_xxyy_xxyyyy,  \
                             tk_xxyy_xxyyyz,  \
                             tk_xxyy_xxyyz,   \
                             tk_xxyy_xxyyzz,  \
                             tk_xxyy_xxyzz,   \
                             tk_xxyy_xxyzzz,  \
                             tk_xxyy_xyyyy,   \
                             tk_xxyy_xyyyyy,  \
                             tk_xxyy_xyyyyz,  \
                             tk_xxyy_xyyyz,   \
                             tk_xxyy_xyyyzz,  \
                             tk_xxyy_xyyzz,   \
                             tk_xxyy_xyyzzz,  \
                             tk_xxyy_xyzzz,   \
                             tk_xxyy_xyzzzz,  \
                             tk_xxyy_yyyyy,   \
                             tk_xxyy_yyyyyy,  \
                             tk_xxyy_yyyyyz,  \
                             tk_xxyy_yyyyz,   \
                             tk_xxyy_yyyyzz,  \
                             tk_xxyy_yyyzz,   \
                             tk_xxyy_yyyzzz,  \
                             tk_xxyy_yyzzz,   \
                             tk_xxyy_yyzzzz,  \
                             tk_xxyy_yzzzz,   \
                             tk_xxyy_yzzzzz,  \
                             tk_xxyy_zzzzzz,  \
                             tk_xyy_xxxxxy,   \
                             tk_xyy_xxxxyy,   \
                             tk_xyy_xxxxyz,   \
                             tk_xyy_xxxyyy,   \
                             tk_xyy_xxxyyz,   \
                             tk_xyy_xxxyzz,   \
                             tk_xyy_xxyyyy,   \
                             tk_xyy_xxyyyz,   \
                             tk_xyy_xxyyzz,   \
                             tk_xyy_xxyzzz,   \
                             tk_xyy_xyyyyy,   \
                             tk_xyy_xyyyyz,   \
                             tk_xyy_xyyyzz,   \
                             tk_xyy_xyyzzz,   \
                             tk_xyy_xyzzzz,   \
                             tk_xyy_yyyyyy,   \
                             tk_xyy_yyyyyz,   \
                             tk_xyy_yyyyzz,   \
                             tk_xyy_yyyzzz,   \
                             tk_xyy_yyzzzz,   \
                             tk_xyy_yzzzzz,   \
                             tk_xyy_zzzzzz,   \
                             ts_xxx_xxxxxx,   \
                             ts_xxx_xxxxxz,   \
                             ts_xxx_xxxxzz,   \
                             ts_xxx_xxxzzz,   \
                             ts_xxx_xxzzzz,   \
                             ts_xxx_xzzzzz,   \
                             ts_xxxyy_xxxxxx, \
                             ts_xxxyy_xxxxxy, \
                             ts_xxxyy_xxxxxz, \
                             ts_xxxyy_xxxxyy, \
                             ts_xxxyy_xxxxyz, \
                             ts_xxxyy_xxxxzz, \
                             ts_xxxyy_xxxyyy, \
                             ts_xxxyy_xxxyyz, \
                             ts_xxxyy_xxxyzz, \
                             ts_xxxyy_xxxzzz, \
                             ts_xxxyy_xxyyyy, \
                             ts_xxxyy_xxyyyz, \
                             ts_xxxyy_xxyyzz, \
                             ts_xxxyy_xxyzzz, \
                             ts_xxxyy_xxzzzz, \
                             ts_xxxyy_xyyyyy, \
                             ts_xxxyy_xyyyyz, \
                             ts_xxxyy_xyyyzz, \
                             ts_xxxyy_xyyzzz, \
                             ts_xxxyy_xyzzzz, \
                             ts_xxxyy_xzzzzz, \
                             ts_xxxyy_yyyyyy, \
                             ts_xxxyy_yyyyyz, \
                             ts_xxxyy_yyyyzz, \
                             ts_xxxyy_yyyzzz, \
                             ts_xxxyy_yyzzzz, \
                             ts_xxxyy_yzzzzz, \
                             ts_xxxyy_zzzzzz, \
                             ts_xyy_xxxxxy,   \
                             ts_xyy_xxxxyy,   \
                             ts_xyy_xxxxyz,   \
                             ts_xyy_xxxyyy,   \
                             ts_xyy_xxxyyz,   \
                             ts_xyy_xxxyzz,   \
                             ts_xyy_xxyyyy,   \
                             ts_xyy_xxyyyz,   \
                             ts_xyy_xxyyzz,   \
                             ts_xyy_xxyzzz,   \
                             ts_xyy_xyyyyy,   \
                             ts_xyy_xyyyyz,   \
                             ts_xyy_xyyyzz,   \
                             ts_xyy_xyyzzz,   \
                             ts_xyy_xyzzzz,   \
                             ts_xyy_yyyyyy,   \
                             ts_xyy_yyyyyz,   \
                             ts_xyy_yyyyzz,   \
                             ts_xyy_yyyzzz,   \
                             ts_xyy_yyzzzz,   \
                             ts_xyy_yzzzzz,   \
                             ts_xyy_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxyy_xxxxxx[i] =
            -2.0 * ts_xxx_xxxxxx[i] * fbe_0 * fz_0 + tk_xxx_xxxxxx[i] * fe_0 + tk_xxxy_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxxyy_xxxxxx[i] * fz_0;

        tk_xxxyy_xxxxxy[i] = -4.0 * ts_xyy_xxxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxxxxy[i] * fe_0 + 5.0 * tk_xxyy_xxxxy[i] * fe_0 +
                             tk_xxyy_xxxxxy[i] * pa_x[i] + 2.0 * ts_xxxyy_xxxxxy[i] * fz_0;

        tk_xxxyy_xxxxxz[i] =
            -2.0 * ts_xxx_xxxxxz[i] * fbe_0 * fz_0 + tk_xxx_xxxxxz[i] * fe_0 + tk_xxxy_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxxyy_xxxxxz[i] * fz_0;

        tk_xxxyy_xxxxyy[i] = -4.0 * ts_xyy_xxxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxxxyy[i] * fe_0 + 4.0 * tk_xxyy_xxxyy[i] * fe_0 +
                             tk_xxyy_xxxxyy[i] * pa_x[i] + 2.0 * ts_xxxyy_xxxxyy[i] * fz_0;

        tk_xxxyy_xxxxyz[i] = -4.0 * ts_xyy_xxxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxxxyz[i] * fe_0 + 4.0 * tk_xxyy_xxxyz[i] * fe_0 +
                             tk_xxyy_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxxyy_xxxxyz[i] * fz_0;

        tk_xxxyy_xxxxzz[i] =
            -2.0 * ts_xxx_xxxxzz[i] * fbe_0 * fz_0 + tk_xxx_xxxxzz[i] * fe_0 + tk_xxxy_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxxyy_xxxxzz[i] * fz_0;

        tk_xxxyy_xxxyyy[i] = -4.0 * ts_xyy_xxxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxxyyy[i] * fe_0 + 3.0 * tk_xxyy_xxyyy[i] * fe_0 +
                             tk_xxyy_xxxyyy[i] * pa_x[i] + 2.0 * ts_xxxyy_xxxyyy[i] * fz_0;

        tk_xxxyy_xxxyyz[i] = -4.0 * ts_xyy_xxxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxxyyz[i] * fe_0 + 3.0 * tk_xxyy_xxyyz[i] * fe_0 +
                             tk_xxyy_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxxyy_xxxyyz[i] * fz_0;

        tk_xxxyy_xxxyzz[i] = -4.0 * ts_xyy_xxxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxxyzz[i] * fe_0 + 3.0 * tk_xxyy_xxyzz[i] * fe_0 +
                             tk_xxyy_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxxyy_xxxyzz[i] * fz_0;

        tk_xxxyy_xxxzzz[i] =
            -2.0 * ts_xxx_xxxzzz[i] * fbe_0 * fz_0 + tk_xxx_xxxzzz[i] * fe_0 + tk_xxxy_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxxyy_xxxzzz[i] * fz_0;

        tk_xxxyy_xxyyyy[i] = -4.0 * ts_xyy_xxyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxyyyy[i] * fe_0 + 2.0 * tk_xxyy_xyyyy[i] * fe_0 +
                             tk_xxyy_xxyyyy[i] * pa_x[i] + 2.0 * ts_xxxyy_xxyyyy[i] * fz_0;

        tk_xxxyy_xxyyyz[i] = -4.0 * ts_xyy_xxyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxyyyz[i] * fe_0 + 2.0 * tk_xxyy_xyyyz[i] * fe_0 +
                             tk_xxyy_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxxyy_xxyyyz[i] * fz_0;

        tk_xxxyy_xxyyzz[i] = -4.0 * ts_xyy_xxyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxyyzz[i] * fe_0 + 2.0 * tk_xxyy_xyyzz[i] * fe_0 +
                             tk_xxyy_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxxyy_xxyyzz[i] * fz_0;

        tk_xxxyy_xxyzzz[i] = -4.0 * ts_xyy_xxyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xxyzzz[i] * fe_0 + 2.0 * tk_xxyy_xyzzz[i] * fe_0 +
                             tk_xxyy_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_xxyzzz[i] * fz_0;

        tk_xxxyy_xxzzzz[i] =
            -2.0 * ts_xxx_xxzzzz[i] * fbe_0 * fz_0 + tk_xxx_xxzzzz[i] * fe_0 + tk_xxxy_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxxyy_xxzzzz[i] * fz_0;

        tk_xxxyy_xyyyyy[i] = -4.0 * ts_xyy_xyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyyyyy[i] * fe_0 + tk_xxyy_yyyyy[i] * fe_0 +
                             tk_xxyy_xyyyyy[i] * pa_x[i] + 2.0 * ts_xxxyy_xyyyyy[i] * fz_0;

        tk_xxxyy_xyyyyz[i] = -4.0 * ts_xyy_xyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyyyyz[i] * fe_0 + tk_xxyy_yyyyz[i] * fe_0 +
                             tk_xxyy_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxxyy_xyyyyz[i] * fz_0;

        tk_xxxyy_xyyyzz[i] = -4.0 * ts_xyy_xyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyyyzz[i] * fe_0 + tk_xxyy_yyyzz[i] * fe_0 +
                             tk_xxyy_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxxyy_xyyyzz[i] * fz_0;

        tk_xxxyy_xyyzzz[i] = -4.0 * ts_xyy_xyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyyzzz[i] * fe_0 + tk_xxyy_yyzzz[i] * fe_0 +
                             tk_xxyy_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_xyyzzz[i] * fz_0;

        tk_xxxyy_xyzzzz[i] = -4.0 * ts_xyy_xyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xyzzzz[i] * fe_0 + tk_xxyy_yzzzz[i] * fe_0 +
                             tk_xxyy_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_xyzzzz[i] * fz_0;

        tk_xxxyy_xzzzzz[i] =
            -2.0 * ts_xxx_xzzzzz[i] * fbe_0 * fz_0 + tk_xxx_xzzzzz[i] * fe_0 + tk_xxxy_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxxyy_xzzzzz[i] * fz_0;

        tk_xxxyy_yyyyyy[i] =
            -4.0 * ts_xyy_yyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyyyyy[i] * fe_0 + tk_xxyy_yyyyyy[i] * pa_x[i] + 2.0 * ts_xxxyy_yyyyyy[i] * fz_0;

        tk_xxxyy_yyyyyz[i] =
            -4.0 * ts_xyy_yyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyyyyz[i] * fe_0 + tk_xxyy_yyyyyz[i] * pa_x[i] + 2.0 * ts_xxxyy_yyyyyz[i] * fz_0;

        tk_xxxyy_yyyyzz[i] =
            -4.0 * ts_xyy_yyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyyyzz[i] * fe_0 + tk_xxyy_yyyyzz[i] * pa_x[i] + 2.0 * ts_xxxyy_yyyyzz[i] * fz_0;

        tk_xxxyy_yyyzzz[i] =
            -4.0 * ts_xyy_yyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyyzzz[i] * fe_0 + tk_xxyy_yyyzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_yyyzzz[i] * fz_0;

        tk_xxxyy_yyzzzz[i] =
            -4.0 * ts_xyy_yyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yyzzzz[i] * fe_0 + tk_xxyy_yyzzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_yyzzzz[i] * fz_0;

        tk_xxxyy_yzzzzz[i] =
            -4.0 * ts_xyy_yzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yzzzzz[i] * fe_0 + tk_xxyy_yzzzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_yzzzzz[i] * fz_0;

        tk_xxxyy_zzzzzz[i] =
            -4.0 * ts_xyy_zzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_zzzzzz[i] * fe_0 + tk_xxyy_zzzzzz[i] * pa_x[i] + 2.0 * ts_xxxyy_zzzzzz[i] * fz_0;
    }

    // Set up 112-140 components of targeted buffer : HI

    auto tk_xxxyz_xxxxxx = pbuffer.data(idx_kin_hi + 112);

    auto tk_xxxyz_xxxxxy = pbuffer.data(idx_kin_hi + 113);

    auto tk_xxxyz_xxxxxz = pbuffer.data(idx_kin_hi + 114);

    auto tk_xxxyz_xxxxyy = pbuffer.data(idx_kin_hi + 115);

    auto tk_xxxyz_xxxxyz = pbuffer.data(idx_kin_hi + 116);

    auto tk_xxxyz_xxxxzz = pbuffer.data(idx_kin_hi + 117);

    auto tk_xxxyz_xxxyyy = pbuffer.data(idx_kin_hi + 118);

    auto tk_xxxyz_xxxyyz = pbuffer.data(idx_kin_hi + 119);

    auto tk_xxxyz_xxxyzz = pbuffer.data(idx_kin_hi + 120);

    auto tk_xxxyz_xxxzzz = pbuffer.data(idx_kin_hi + 121);

    auto tk_xxxyz_xxyyyy = pbuffer.data(idx_kin_hi + 122);

    auto tk_xxxyz_xxyyyz = pbuffer.data(idx_kin_hi + 123);

    auto tk_xxxyz_xxyyzz = pbuffer.data(idx_kin_hi + 124);

    auto tk_xxxyz_xxyzzz = pbuffer.data(idx_kin_hi + 125);

    auto tk_xxxyz_xxzzzz = pbuffer.data(idx_kin_hi + 126);

    auto tk_xxxyz_xyyyyy = pbuffer.data(idx_kin_hi + 127);

    auto tk_xxxyz_xyyyyz = pbuffer.data(idx_kin_hi + 128);

    auto tk_xxxyz_xyyyzz = pbuffer.data(idx_kin_hi + 129);

    auto tk_xxxyz_xyyzzz = pbuffer.data(idx_kin_hi + 130);

    auto tk_xxxyz_xyzzzz = pbuffer.data(idx_kin_hi + 131);

    auto tk_xxxyz_xzzzzz = pbuffer.data(idx_kin_hi + 132);

    auto tk_xxxyz_yyyyyy = pbuffer.data(idx_kin_hi + 133);

    auto tk_xxxyz_yyyyyz = pbuffer.data(idx_kin_hi + 134);

    auto tk_xxxyz_yyyyzz = pbuffer.data(idx_kin_hi + 135);

    auto tk_xxxyz_yyyzzz = pbuffer.data(idx_kin_hi + 136);

    auto tk_xxxyz_yyzzzz = pbuffer.data(idx_kin_hi + 137);

    auto tk_xxxyz_yzzzzz = pbuffer.data(idx_kin_hi + 138);

    auto tk_xxxyz_zzzzzz = pbuffer.data(idx_kin_hi + 139);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tk_xxxy_xxxxxy,  \
                             tk_xxxy_xxxxyy,  \
                             tk_xxxy_xxxyyy,  \
                             tk_xxxy_xxyyyy,  \
                             tk_xxxy_xyyyyy,  \
                             tk_xxxy_yyyyyy,  \
                             tk_xxxyz_xxxxxx, \
                             tk_xxxyz_xxxxxy, \
                             tk_xxxyz_xxxxxz, \
                             tk_xxxyz_xxxxyy, \
                             tk_xxxyz_xxxxyz, \
                             tk_xxxyz_xxxxzz, \
                             tk_xxxyz_xxxyyy, \
                             tk_xxxyz_xxxyyz, \
                             tk_xxxyz_xxxyzz, \
                             tk_xxxyz_xxxzzz, \
                             tk_xxxyz_xxyyyy, \
                             tk_xxxyz_xxyyyz, \
                             tk_xxxyz_xxyyzz, \
                             tk_xxxyz_xxyzzz, \
                             tk_xxxyz_xxzzzz, \
                             tk_xxxyz_xyyyyy, \
                             tk_xxxyz_xyyyyz, \
                             tk_xxxyz_xyyyzz, \
                             tk_xxxyz_xyyzzz, \
                             tk_xxxyz_xyzzzz, \
                             tk_xxxyz_xzzzzz, \
                             tk_xxxyz_yyyyyy, \
                             tk_xxxyz_yyyyyz, \
                             tk_xxxyz_yyyyzz, \
                             tk_xxxyz_yyyzzz, \
                             tk_xxxyz_yyzzzz, \
                             tk_xxxyz_yzzzzz, \
                             tk_xxxyz_zzzzzz, \
                             tk_xxxz_xxxxxx,  \
                             tk_xxxz_xxxxxz,  \
                             tk_xxxz_xxxxyz,  \
                             tk_xxxz_xxxxz,   \
                             tk_xxxz_xxxxzz,  \
                             tk_xxxz_xxxyyz,  \
                             tk_xxxz_xxxyz,   \
                             tk_xxxz_xxxyzz,  \
                             tk_xxxz_xxxzz,   \
                             tk_xxxz_xxxzzz,  \
                             tk_xxxz_xxyyyz,  \
                             tk_xxxz_xxyyz,   \
                             tk_xxxz_xxyyzz,  \
                             tk_xxxz_xxyzz,   \
                             tk_xxxz_xxyzzz,  \
                             tk_xxxz_xxzzz,   \
                             tk_xxxz_xxzzzz,  \
                             tk_xxxz_xyyyyz,  \
                             tk_xxxz_xyyyz,   \
                             tk_xxxz_xyyyzz,  \
                             tk_xxxz_xyyzz,   \
                             tk_xxxz_xyyzzz,  \
                             tk_xxxz_xyzzz,   \
                             tk_xxxz_xyzzzz,  \
                             tk_xxxz_xzzzz,   \
                             tk_xxxz_xzzzzz,  \
                             tk_xxxz_yyyyyz,  \
                             tk_xxxz_yyyyz,   \
                             tk_xxxz_yyyyzz,  \
                             tk_xxxz_yyyzz,   \
                             tk_xxxz_yyyzzz,  \
                             tk_xxxz_yyzzz,   \
                             tk_xxxz_yyzzzz,  \
                             tk_xxxz_yzzzz,   \
                             tk_xxxz_yzzzzz,  \
                             tk_xxxz_zzzzz,   \
                             tk_xxxz_zzzzzz,  \
                             ts_xxxyz_xxxxxx, \
                             ts_xxxyz_xxxxxy, \
                             ts_xxxyz_xxxxxz, \
                             ts_xxxyz_xxxxyy, \
                             ts_xxxyz_xxxxyz, \
                             ts_xxxyz_xxxxzz, \
                             ts_xxxyz_xxxyyy, \
                             ts_xxxyz_xxxyyz, \
                             ts_xxxyz_xxxyzz, \
                             ts_xxxyz_xxxzzz, \
                             ts_xxxyz_xxyyyy, \
                             ts_xxxyz_xxyyyz, \
                             ts_xxxyz_xxyyzz, \
                             ts_xxxyz_xxyzzz, \
                             ts_xxxyz_xxzzzz, \
                             ts_xxxyz_xyyyyy, \
                             ts_xxxyz_xyyyyz, \
                             ts_xxxyz_xyyyzz, \
                             ts_xxxyz_xyyzzz, \
                             ts_xxxyz_xyzzzz, \
                             ts_xxxyz_xzzzzz, \
                             ts_xxxyz_yyyyyy, \
                             ts_xxxyz_yyyyyz, \
                             ts_xxxyz_yyyyzz, \
                             ts_xxxyz_yyyzzz, \
                             ts_xxxyz_yyzzzz, \
                             ts_xxxyz_yzzzzz, \
                             ts_xxxyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyz_xxxxxx[i] = tk_xxxz_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxxxx[i] * fz_0;

        tk_xxxyz_xxxxxy[i] = tk_xxxy_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxxyz_xxxxxy[i] * fz_0;

        tk_xxxyz_xxxxxz[i] = tk_xxxz_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxxxz[i] * fz_0;

        tk_xxxyz_xxxxyy[i] = tk_xxxy_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxxyz_xxxxyy[i] * fz_0;

        tk_xxxyz_xxxxyz[i] = tk_xxxz_xxxxz[i] * fe_0 + tk_xxxz_xxxxyz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxxyz[i] * fz_0;

        tk_xxxyz_xxxxzz[i] = tk_xxxz_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxxzz[i] * fz_0;

        tk_xxxyz_xxxyyy[i] = tk_xxxy_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxxyz_xxxyyy[i] * fz_0;

        tk_xxxyz_xxxyyz[i] = 2.0 * tk_xxxz_xxxyz[i] * fe_0 + tk_xxxz_xxxyyz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxyyz[i] * fz_0;

        tk_xxxyz_xxxyzz[i] = tk_xxxz_xxxzz[i] * fe_0 + tk_xxxz_xxxyzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxyzz[i] * fz_0;

        tk_xxxyz_xxxzzz[i] = tk_xxxz_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxxzzz[i] * fz_0;

        tk_xxxyz_xxyyyy[i] = tk_xxxy_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxxyz_xxyyyy[i] * fz_0;

        tk_xxxyz_xxyyyz[i] = 3.0 * tk_xxxz_xxyyz[i] * fe_0 + tk_xxxz_xxyyyz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxyyyz[i] * fz_0;

        tk_xxxyz_xxyyzz[i] = 2.0 * tk_xxxz_xxyzz[i] * fe_0 + tk_xxxz_xxyyzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxyyzz[i] * fz_0;

        tk_xxxyz_xxyzzz[i] = tk_xxxz_xxzzz[i] * fe_0 + tk_xxxz_xxyzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxyzzz[i] * fz_0;

        tk_xxxyz_xxzzzz[i] = tk_xxxz_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xxzzzz[i] * fz_0;

        tk_xxxyz_xyyyyy[i] = tk_xxxy_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxxyz_xyyyyy[i] * fz_0;

        tk_xxxyz_xyyyyz[i] = 4.0 * tk_xxxz_xyyyz[i] * fe_0 + tk_xxxz_xyyyyz[i] * pa_y[i] + 2.0 * ts_xxxyz_xyyyyz[i] * fz_0;

        tk_xxxyz_xyyyzz[i] = 3.0 * tk_xxxz_xyyzz[i] * fe_0 + tk_xxxz_xyyyzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xyyyzz[i] * fz_0;

        tk_xxxyz_xyyzzz[i] = 2.0 * tk_xxxz_xyzzz[i] * fe_0 + tk_xxxz_xyyzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xyyzzz[i] * fz_0;

        tk_xxxyz_xyzzzz[i] = tk_xxxz_xzzzz[i] * fe_0 + tk_xxxz_xyzzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xyzzzz[i] * fz_0;

        tk_xxxyz_xzzzzz[i] = tk_xxxz_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_xzzzzz[i] * fz_0;

        tk_xxxyz_yyyyyy[i] = tk_xxxy_yyyyyy[i] * pa_z[i] + 2.0 * ts_xxxyz_yyyyyy[i] * fz_0;

        tk_xxxyz_yyyyyz[i] = 5.0 * tk_xxxz_yyyyz[i] * fe_0 + tk_xxxz_yyyyyz[i] * pa_y[i] + 2.0 * ts_xxxyz_yyyyyz[i] * fz_0;

        tk_xxxyz_yyyyzz[i] = 4.0 * tk_xxxz_yyyzz[i] * fe_0 + tk_xxxz_yyyyzz[i] * pa_y[i] + 2.0 * ts_xxxyz_yyyyzz[i] * fz_0;

        tk_xxxyz_yyyzzz[i] = 3.0 * tk_xxxz_yyzzz[i] * fe_0 + tk_xxxz_yyyzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_yyyzzz[i] * fz_0;

        tk_xxxyz_yyzzzz[i] = 2.0 * tk_xxxz_yzzzz[i] * fe_0 + tk_xxxz_yyzzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_yyzzzz[i] * fz_0;

        tk_xxxyz_yzzzzz[i] = tk_xxxz_zzzzz[i] * fe_0 + tk_xxxz_yzzzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_yzzzzz[i] * fz_0;

        tk_xxxyz_zzzzzz[i] = tk_xxxz_zzzzzz[i] * pa_y[i] + 2.0 * ts_xxxyz_zzzzzz[i] * fz_0;
    }

    // Set up 140-168 components of targeted buffer : HI

    auto tk_xxxzz_xxxxxx = pbuffer.data(idx_kin_hi + 140);

    auto tk_xxxzz_xxxxxy = pbuffer.data(idx_kin_hi + 141);

    auto tk_xxxzz_xxxxxz = pbuffer.data(idx_kin_hi + 142);

    auto tk_xxxzz_xxxxyy = pbuffer.data(idx_kin_hi + 143);

    auto tk_xxxzz_xxxxyz = pbuffer.data(idx_kin_hi + 144);

    auto tk_xxxzz_xxxxzz = pbuffer.data(idx_kin_hi + 145);

    auto tk_xxxzz_xxxyyy = pbuffer.data(idx_kin_hi + 146);

    auto tk_xxxzz_xxxyyz = pbuffer.data(idx_kin_hi + 147);

    auto tk_xxxzz_xxxyzz = pbuffer.data(idx_kin_hi + 148);

    auto tk_xxxzz_xxxzzz = pbuffer.data(idx_kin_hi + 149);

    auto tk_xxxzz_xxyyyy = pbuffer.data(idx_kin_hi + 150);

    auto tk_xxxzz_xxyyyz = pbuffer.data(idx_kin_hi + 151);

    auto tk_xxxzz_xxyyzz = pbuffer.data(idx_kin_hi + 152);

    auto tk_xxxzz_xxyzzz = pbuffer.data(idx_kin_hi + 153);

    auto tk_xxxzz_xxzzzz = pbuffer.data(idx_kin_hi + 154);

    auto tk_xxxzz_xyyyyy = pbuffer.data(idx_kin_hi + 155);

    auto tk_xxxzz_xyyyyz = pbuffer.data(idx_kin_hi + 156);

    auto tk_xxxzz_xyyyzz = pbuffer.data(idx_kin_hi + 157);

    auto tk_xxxzz_xyyzzz = pbuffer.data(idx_kin_hi + 158);

    auto tk_xxxzz_xyzzzz = pbuffer.data(idx_kin_hi + 159);

    auto tk_xxxzz_xzzzzz = pbuffer.data(idx_kin_hi + 160);

    auto tk_xxxzz_yyyyyy = pbuffer.data(idx_kin_hi + 161);

    auto tk_xxxzz_yyyyyz = pbuffer.data(idx_kin_hi + 162);

    auto tk_xxxzz_yyyyzz = pbuffer.data(idx_kin_hi + 163);

    auto tk_xxxzz_yyyzzz = pbuffer.data(idx_kin_hi + 164);

    auto tk_xxxzz_yyzzzz = pbuffer.data(idx_kin_hi + 165);

    auto tk_xxxzz_yzzzzz = pbuffer.data(idx_kin_hi + 166);

    auto tk_xxxzz_zzzzzz = pbuffer.data(idx_kin_hi + 167);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tk_xxx_xxxxxx,   \
                             tk_xxx_xxxxxy,   \
                             tk_xxx_xxxxyy,   \
                             tk_xxx_xxxyyy,   \
                             tk_xxx_xxyyyy,   \
                             tk_xxx_xyyyyy,   \
                             tk_xxxz_xxxxxx,  \
                             tk_xxxz_xxxxxy,  \
                             tk_xxxz_xxxxyy,  \
                             tk_xxxz_xxxyyy,  \
                             tk_xxxz_xxyyyy,  \
                             tk_xxxz_xyyyyy,  \
                             tk_xxxzz_xxxxxx, \
                             tk_xxxzz_xxxxxy, \
                             tk_xxxzz_xxxxxz, \
                             tk_xxxzz_xxxxyy, \
                             tk_xxxzz_xxxxyz, \
                             tk_xxxzz_xxxxzz, \
                             tk_xxxzz_xxxyyy, \
                             tk_xxxzz_xxxyyz, \
                             tk_xxxzz_xxxyzz, \
                             tk_xxxzz_xxxzzz, \
                             tk_xxxzz_xxyyyy, \
                             tk_xxxzz_xxyyyz, \
                             tk_xxxzz_xxyyzz, \
                             tk_xxxzz_xxyzzz, \
                             tk_xxxzz_xxzzzz, \
                             tk_xxxzz_xyyyyy, \
                             tk_xxxzz_xyyyyz, \
                             tk_xxxzz_xyyyzz, \
                             tk_xxxzz_xyyzzz, \
                             tk_xxxzz_xyzzzz, \
                             tk_xxxzz_xzzzzz, \
                             tk_xxxzz_yyyyyy, \
                             tk_xxxzz_yyyyyz, \
                             tk_xxxzz_yyyyzz, \
                             tk_xxxzz_yyyzzz, \
                             tk_xxxzz_yyzzzz, \
                             tk_xxxzz_yzzzzz, \
                             tk_xxxzz_zzzzzz, \
                             tk_xxzz_xxxxxz,  \
                             tk_xxzz_xxxxyz,  \
                             tk_xxzz_xxxxz,   \
                             tk_xxzz_xxxxzz,  \
                             tk_xxzz_xxxyyz,  \
                             tk_xxzz_xxxyz,   \
                             tk_xxzz_xxxyzz,  \
                             tk_xxzz_xxxzz,   \
                             tk_xxzz_xxxzzz,  \
                             tk_xxzz_xxyyyz,  \
                             tk_xxzz_xxyyz,   \
                             tk_xxzz_xxyyzz,  \
                             tk_xxzz_xxyzz,   \
                             tk_xxzz_xxyzzz,  \
                             tk_xxzz_xxzzz,   \
                             tk_xxzz_xxzzzz,  \
                             tk_xxzz_xyyyyz,  \
                             tk_xxzz_xyyyz,   \
                             tk_xxzz_xyyyzz,  \
                             tk_xxzz_xyyzz,   \
                             tk_xxzz_xyyzzz,  \
                             tk_xxzz_xyzzz,   \
                             tk_xxzz_xyzzzz,  \
                             tk_xxzz_xzzzz,   \
                             tk_xxzz_xzzzzz,  \
                             tk_xxzz_yyyyyy,  \
                             tk_xxzz_yyyyyz,  \
                             tk_xxzz_yyyyz,   \
                             tk_xxzz_yyyyzz,  \
                             tk_xxzz_yyyzz,   \
                             tk_xxzz_yyyzzz,  \
                             tk_xxzz_yyzzz,   \
                             tk_xxzz_yyzzzz,  \
                             tk_xxzz_yzzzz,   \
                             tk_xxzz_yzzzzz,  \
                             tk_xxzz_zzzzz,   \
                             tk_xxzz_zzzzzz,  \
                             tk_xzz_xxxxxz,   \
                             tk_xzz_xxxxyz,   \
                             tk_xzz_xxxxzz,   \
                             tk_xzz_xxxyyz,   \
                             tk_xzz_xxxyzz,   \
                             tk_xzz_xxxzzz,   \
                             tk_xzz_xxyyyz,   \
                             tk_xzz_xxyyzz,   \
                             tk_xzz_xxyzzz,   \
                             tk_xzz_xxzzzz,   \
                             tk_xzz_xyyyyz,   \
                             tk_xzz_xyyyzz,   \
                             tk_xzz_xyyzzz,   \
                             tk_xzz_xyzzzz,   \
                             tk_xzz_xzzzzz,   \
                             tk_xzz_yyyyyy,   \
                             tk_xzz_yyyyyz,   \
                             tk_xzz_yyyyzz,   \
                             tk_xzz_yyyzzz,   \
                             tk_xzz_yyzzzz,   \
                             tk_xzz_yzzzzz,   \
                             tk_xzz_zzzzzz,   \
                             ts_xxx_xxxxxx,   \
                             ts_xxx_xxxxxy,   \
                             ts_xxx_xxxxyy,   \
                             ts_xxx_xxxyyy,   \
                             ts_xxx_xxyyyy,   \
                             ts_xxx_xyyyyy,   \
                             ts_xxxzz_xxxxxx, \
                             ts_xxxzz_xxxxxy, \
                             ts_xxxzz_xxxxxz, \
                             ts_xxxzz_xxxxyy, \
                             ts_xxxzz_xxxxyz, \
                             ts_xxxzz_xxxxzz, \
                             ts_xxxzz_xxxyyy, \
                             ts_xxxzz_xxxyyz, \
                             ts_xxxzz_xxxyzz, \
                             ts_xxxzz_xxxzzz, \
                             ts_xxxzz_xxyyyy, \
                             ts_xxxzz_xxyyyz, \
                             ts_xxxzz_xxyyzz, \
                             ts_xxxzz_xxyzzz, \
                             ts_xxxzz_xxzzzz, \
                             ts_xxxzz_xyyyyy, \
                             ts_xxxzz_xyyyyz, \
                             ts_xxxzz_xyyyzz, \
                             ts_xxxzz_xyyzzz, \
                             ts_xxxzz_xyzzzz, \
                             ts_xxxzz_xzzzzz, \
                             ts_xxxzz_yyyyyy, \
                             ts_xxxzz_yyyyyz, \
                             ts_xxxzz_yyyyzz, \
                             ts_xxxzz_yyyzzz, \
                             ts_xxxzz_yyzzzz, \
                             ts_xxxzz_yzzzzz, \
                             ts_xxxzz_zzzzzz, \
                             ts_xzz_xxxxxz,   \
                             ts_xzz_xxxxyz,   \
                             ts_xzz_xxxxzz,   \
                             ts_xzz_xxxyyz,   \
                             ts_xzz_xxxyzz,   \
                             ts_xzz_xxxzzz,   \
                             ts_xzz_xxyyyz,   \
                             ts_xzz_xxyyzz,   \
                             ts_xzz_xxyzzz,   \
                             ts_xzz_xxzzzz,   \
                             ts_xzz_xyyyyz,   \
                             ts_xzz_xyyyzz,   \
                             ts_xzz_xyyzzz,   \
                             ts_xzz_xyzzzz,   \
                             ts_xzz_xzzzzz,   \
                             ts_xzz_yyyyyy,   \
                             ts_xzz_yyyyyz,   \
                             ts_xzz_yyyyzz,   \
                             ts_xzz_yyyzzz,   \
                             ts_xzz_yyzzzz,   \
                             ts_xzz_yzzzzz,   \
                             ts_xzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxzz_xxxxxx[i] =
            -2.0 * ts_xxx_xxxxxx[i] * fbe_0 * fz_0 + tk_xxx_xxxxxx[i] * fe_0 + tk_xxxz_xxxxxx[i] * pa_z[i] + 2.0 * ts_xxxzz_xxxxxx[i] * fz_0;

        tk_xxxzz_xxxxxy[i] =
            -2.0 * ts_xxx_xxxxxy[i] * fbe_0 * fz_0 + tk_xxx_xxxxxy[i] * fe_0 + tk_xxxz_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxxzz_xxxxxy[i] * fz_0;

        tk_xxxzz_xxxxxz[i] = -4.0 * ts_xzz_xxxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxxxxz[i] * fe_0 + 5.0 * tk_xxzz_xxxxz[i] * fe_0 +
                             tk_xxzz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxxxxz[i] * fz_0;

        tk_xxxzz_xxxxyy[i] =
            -2.0 * ts_xxx_xxxxyy[i] * fbe_0 * fz_0 + tk_xxx_xxxxyy[i] * fe_0 + tk_xxxz_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxxzz_xxxxyy[i] * fz_0;

        tk_xxxzz_xxxxyz[i] = -4.0 * ts_xzz_xxxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxxxyz[i] * fe_0 + 4.0 * tk_xxzz_xxxyz[i] * fe_0 +
                             tk_xxzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxxxyz[i] * fz_0;

        tk_xxxzz_xxxxzz[i] = -4.0 * ts_xzz_xxxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxxxzz[i] * fe_0 + 4.0 * tk_xxzz_xxxzz[i] * fe_0 +
                             tk_xxzz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxxxzz[i] * fz_0;

        tk_xxxzz_xxxyyy[i] =
            -2.0 * ts_xxx_xxxyyy[i] * fbe_0 * fz_0 + tk_xxx_xxxyyy[i] * fe_0 + tk_xxxz_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxxzz_xxxyyy[i] * fz_0;

        tk_xxxzz_xxxyyz[i] = -4.0 * ts_xzz_xxxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxxyyz[i] * fe_0 + 3.0 * tk_xxzz_xxyyz[i] * fe_0 +
                             tk_xxzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxxyyz[i] * fz_0;

        tk_xxxzz_xxxyzz[i] = -4.0 * ts_xzz_xxxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxxyzz[i] * fe_0 + 3.0 * tk_xxzz_xxyzz[i] * fe_0 +
                             tk_xxzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxxyzz[i] * fz_0;

        tk_xxxzz_xxxzzz[i] = -4.0 * ts_xzz_xxxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxxzzz[i] * fe_0 + 3.0 * tk_xxzz_xxzzz[i] * fe_0 +
                             tk_xxzz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxxzzz[i] * fz_0;

        tk_xxxzz_xxyyyy[i] =
            -2.0 * ts_xxx_xxyyyy[i] * fbe_0 * fz_0 + tk_xxx_xxyyyy[i] * fe_0 + tk_xxxz_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxxzz_xxyyyy[i] * fz_0;

        tk_xxxzz_xxyyyz[i] = -4.0 * ts_xzz_xxyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxyyyz[i] * fe_0 + 2.0 * tk_xxzz_xyyyz[i] * fe_0 +
                             tk_xxzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxyyyz[i] * fz_0;

        tk_xxxzz_xxyyzz[i] = -4.0 * ts_xzz_xxyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxyyzz[i] * fe_0 + 2.0 * tk_xxzz_xyyzz[i] * fe_0 +
                             tk_xxzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxyyzz[i] * fz_0;

        tk_xxxzz_xxyzzz[i] = -4.0 * ts_xzz_xxyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxyzzz[i] * fe_0 + 2.0 * tk_xxzz_xyzzz[i] * fe_0 +
                             tk_xxzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxyzzz[i] * fz_0;

        tk_xxxzz_xxzzzz[i] = -4.0 * ts_xzz_xxzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xxzzzz[i] * fe_0 + 2.0 * tk_xxzz_xzzzz[i] * fe_0 +
                             tk_xxzz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xxzzzz[i] * fz_0;

        tk_xxxzz_xyyyyy[i] =
            -2.0 * ts_xxx_xyyyyy[i] * fbe_0 * fz_0 + tk_xxx_xyyyyy[i] * fe_0 + tk_xxxz_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxxzz_xyyyyy[i] * fz_0;

        tk_xxxzz_xyyyyz[i] = -4.0 * ts_xzz_xyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xyyyyz[i] * fe_0 + tk_xxzz_yyyyz[i] * fe_0 +
                             tk_xxzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxxzz_xyyyyz[i] * fz_0;

        tk_xxxzz_xyyyzz[i] = -4.0 * ts_xzz_xyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xyyyzz[i] * fe_0 + tk_xxzz_yyyzz[i] * fe_0 +
                             tk_xxzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xyyyzz[i] * fz_0;

        tk_xxxzz_xyyzzz[i] = -4.0 * ts_xzz_xyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xyyzzz[i] * fe_0 + tk_xxzz_yyzzz[i] * fe_0 +
                             tk_xxzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xyyzzz[i] * fz_0;

        tk_xxxzz_xyzzzz[i] = -4.0 * ts_xzz_xyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xyzzzz[i] * fe_0 + tk_xxzz_yzzzz[i] * fe_0 +
                             tk_xxzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xyzzzz[i] * fz_0;

        tk_xxxzz_xzzzzz[i] = -4.0 * ts_xzz_xzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xzzzzz[i] * fe_0 + tk_xxzz_zzzzz[i] * fe_0 +
                             tk_xxzz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_xzzzzz[i] * fz_0;

        tk_xxxzz_yyyyyy[i] =
            -4.0 * ts_xzz_yyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyyyyy[i] * fe_0 + tk_xxzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xxxzz_yyyyyy[i] * fz_0;

        tk_xxxzz_yyyyyz[i] =
            -4.0 * ts_xzz_yyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyyyyz[i] * fe_0 + tk_xxzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xxxzz_yyyyyz[i] * fz_0;

        tk_xxxzz_yyyyzz[i] =
            -4.0 * ts_xzz_yyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyyyzz[i] * fe_0 + tk_xxzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xxxzz_yyyyzz[i] * fz_0;

        tk_xxxzz_yyyzzz[i] =
            -4.0 * ts_xzz_yyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyyzzz[i] * fe_0 + tk_xxzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_yyyzzz[i] * fz_0;

        tk_xxxzz_yyzzzz[i] =
            -4.0 * ts_xzz_yyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yyzzzz[i] * fe_0 + tk_xxzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_yyzzzz[i] * fz_0;

        tk_xxxzz_yzzzzz[i] =
            -4.0 * ts_xzz_yzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yzzzzz[i] * fe_0 + tk_xxzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_yzzzzz[i] * fz_0;

        tk_xxxzz_zzzzzz[i] =
            -4.0 * ts_xzz_zzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_zzzzzz[i] * fe_0 + tk_xxzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xxxzz_zzzzzz[i] * fz_0;
    }

    // Set up 168-196 components of targeted buffer : HI

    auto tk_xxyyy_xxxxxx = pbuffer.data(idx_kin_hi + 168);

    auto tk_xxyyy_xxxxxy = pbuffer.data(idx_kin_hi + 169);

    auto tk_xxyyy_xxxxxz = pbuffer.data(idx_kin_hi + 170);

    auto tk_xxyyy_xxxxyy = pbuffer.data(idx_kin_hi + 171);

    auto tk_xxyyy_xxxxyz = pbuffer.data(idx_kin_hi + 172);

    auto tk_xxyyy_xxxxzz = pbuffer.data(idx_kin_hi + 173);

    auto tk_xxyyy_xxxyyy = pbuffer.data(idx_kin_hi + 174);

    auto tk_xxyyy_xxxyyz = pbuffer.data(idx_kin_hi + 175);

    auto tk_xxyyy_xxxyzz = pbuffer.data(idx_kin_hi + 176);

    auto tk_xxyyy_xxxzzz = pbuffer.data(idx_kin_hi + 177);

    auto tk_xxyyy_xxyyyy = pbuffer.data(idx_kin_hi + 178);

    auto tk_xxyyy_xxyyyz = pbuffer.data(idx_kin_hi + 179);

    auto tk_xxyyy_xxyyzz = pbuffer.data(idx_kin_hi + 180);

    auto tk_xxyyy_xxyzzz = pbuffer.data(idx_kin_hi + 181);

    auto tk_xxyyy_xxzzzz = pbuffer.data(idx_kin_hi + 182);

    auto tk_xxyyy_xyyyyy = pbuffer.data(idx_kin_hi + 183);

    auto tk_xxyyy_xyyyyz = pbuffer.data(idx_kin_hi + 184);

    auto tk_xxyyy_xyyyzz = pbuffer.data(idx_kin_hi + 185);

    auto tk_xxyyy_xyyzzz = pbuffer.data(idx_kin_hi + 186);

    auto tk_xxyyy_xyzzzz = pbuffer.data(idx_kin_hi + 187);

    auto tk_xxyyy_xzzzzz = pbuffer.data(idx_kin_hi + 188);

    auto tk_xxyyy_yyyyyy = pbuffer.data(idx_kin_hi + 189);

    auto tk_xxyyy_yyyyyz = pbuffer.data(idx_kin_hi + 190);

    auto tk_xxyyy_yyyyzz = pbuffer.data(idx_kin_hi + 191);

    auto tk_xxyyy_yyyzzz = pbuffer.data(idx_kin_hi + 192);

    auto tk_xxyyy_yyzzzz = pbuffer.data(idx_kin_hi + 193);

    auto tk_xxyyy_yzzzzz = pbuffer.data(idx_kin_hi + 194);

    auto tk_xxyyy_zzzzzz = pbuffer.data(idx_kin_hi + 195);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tk_xxy_xxxxxx,   \
                             tk_xxy_xxxxxz,   \
                             tk_xxy_xxxxzz,   \
                             tk_xxy_xxxzzz,   \
                             tk_xxy_xxzzzz,   \
                             tk_xxy_xzzzzz,   \
                             tk_xxyy_xxxxxx,  \
                             tk_xxyy_xxxxxz,  \
                             tk_xxyy_xxxxzz,  \
                             tk_xxyy_xxxzzz,  \
                             tk_xxyy_xxzzzz,  \
                             tk_xxyy_xzzzzz,  \
                             tk_xxyyy_xxxxxx, \
                             tk_xxyyy_xxxxxy, \
                             tk_xxyyy_xxxxxz, \
                             tk_xxyyy_xxxxyy, \
                             tk_xxyyy_xxxxyz, \
                             tk_xxyyy_xxxxzz, \
                             tk_xxyyy_xxxyyy, \
                             tk_xxyyy_xxxyyz, \
                             tk_xxyyy_xxxyzz, \
                             tk_xxyyy_xxxzzz, \
                             tk_xxyyy_xxyyyy, \
                             tk_xxyyy_xxyyyz, \
                             tk_xxyyy_xxyyzz, \
                             tk_xxyyy_xxyzzz, \
                             tk_xxyyy_xxzzzz, \
                             tk_xxyyy_xyyyyy, \
                             tk_xxyyy_xyyyyz, \
                             tk_xxyyy_xyyyzz, \
                             tk_xxyyy_xyyzzz, \
                             tk_xxyyy_xyzzzz, \
                             tk_xxyyy_xzzzzz, \
                             tk_xxyyy_yyyyyy, \
                             tk_xxyyy_yyyyyz, \
                             tk_xxyyy_yyyyzz, \
                             tk_xxyyy_yyyzzz, \
                             tk_xxyyy_yyzzzz, \
                             tk_xxyyy_yzzzzz, \
                             tk_xxyyy_zzzzzz, \
                             tk_xyyy_xxxxxy,  \
                             tk_xyyy_xxxxy,   \
                             tk_xyyy_xxxxyy,  \
                             tk_xyyy_xxxxyz,  \
                             tk_xyyy_xxxyy,   \
                             tk_xyyy_xxxyyy,  \
                             tk_xyyy_xxxyyz,  \
                             tk_xyyy_xxxyz,   \
                             tk_xyyy_xxxyzz,  \
                             tk_xyyy_xxyyy,   \
                             tk_xyyy_xxyyyy,  \
                             tk_xyyy_xxyyyz,  \
                             tk_xyyy_xxyyz,   \
                             tk_xyyy_xxyyzz,  \
                             tk_xyyy_xxyzz,   \
                             tk_xyyy_xxyzzz,  \
                             tk_xyyy_xyyyy,   \
                             tk_xyyy_xyyyyy,  \
                             tk_xyyy_xyyyyz,  \
                             tk_xyyy_xyyyz,   \
                             tk_xyyy_xyyyzz,  \
                             tk_xyyy_xyyzz,   \
                             tk_xyyy_xyyzzz,  \
                             tk_xyyy_xyzzz,   \
                             tk_xyyy_xyzzzz,  \
                             tk_xyyy_yyyyy,   \
                             tk_xyyy_yyyyyy,  \
                             tk_xyyy_yyyyyz,  \
                             tk_xyyy_yyyyz,   \
                             tk_xyyy_yyyyzz,  \
                             tk_xyyy_yyyzz,   \
                             tk_xyyy_yyyzzz,  \
                             tk_xyyy_yyzzz,   \
                             tk_xyyy_yyzzzz,  \
                             tk_xyyy_yzzzz,   \
                             tk_xyyy_yzzzzz,  \
                             tk_xyyy_zzzzzz,  \
                             tk_yyy_xxxxxy,   \
                             tk_yyy_xxxxyy,   \
                             tk_yyy_xxxxyz,   \
                             tk_yyy_xxxyyy,   \
                             tk_yyy_xxxyyz,   \
                             tk_yyy_xxxyzz,   \
                             tk_yyy_xxyyyy,   \
                             tk_yyy_xxyyyz,   \
                             tk_yyy_xxyyzz,   \
                             tk_yyy_xxyzzz,   \
                             tk_yyy_xyyyyy,   \
                             tk_yyy_xyyyyz,   \
                             tk_yyy_xyyyzz,   \
                             tk_yyy_xyyzzz,   \
                             tk_yyy_xyzzzz,   \
                             tk_yyy_yyyyyy,   \
                             tk_yyy_yyyyyz,   \
                             tk_yyy_yyyyzz,   \
                             tk_yyy_yyyzzz,   \
                             tk_yyy_yyzzzz,   \
                             tk_yyy_yzzzzz,   \
                             tk_yyy_zzzzzz,   \
                             ts_xxy_xxxxxx,   \
                             ts_xxy_xxxxxz,   \
                             ts_xxy_xxxxzz,   \
                             ts_xxy_xxxzzz,   \
                             ts_xxy_xxzzzz,   \
                             ts_xxy_xzzzzz,   \
                             ts_xxyyy_xxxxxx, \
                             ts_xxyyy_xxxxxy, \
                             ts_xxyyy_xxxxxz, \
                             ts_xxyyy_xxxxyy, \
                             ts_xxyyy_xxxxyz, \
                             ts_xxyyy_xxxxzz, \
                             ts_xxyyy_xxxyyy, \
                             ts_xxyyy_xxxyyz, \
                             ts_xxyyy_xxxyzz, \
                             ts_xxyyy_xxxzzz, \
                             ts_xxyyy_xxyyyy, \
                             ts_xxyyy_xxyyyz, \
                             ts_xxyyy_xxyyzz, \
                             ts_xxyyy_xxyzzz, \
                             ts_xxyyy_xxzzzz, \
                             ts_xxyyy_xyyyyy, \
                             ts_xxyyy_xyyyyz, \
                             ts_xxyyy_xyyyzz, \
                             ts_xxyyy_xyyzzz, \
                             ts_xxyyy_xyzzzz, \
                             ts_xxyyy_xzzzzz, \
                             ts_xxyyy_yyyyyy, \
                             ts_xxyyy_yyyyyz, \
                             ts_xxyyy_yyyyzz, \
                             ts_xxyyy_yyyzzz, \
                             ts_xxyyy_yyzzzz, \
                             ts_xxyyy_yzzzzz, \
                             ts_xxyyy_zzzzzz, \
                             ts_yyy_xxxxxy,   \
                             ts_yyy_xxxxyy,   \
                             ts_yyy_xxxxyz,   \
                             ts_yyy_xxxyyy,   \
                             ts_yyy_xxxyyz,   \
                             ts_yyy_xxxyzz,   \
                             ts_yyy_xxyyyy,   \
                             ts_yyy_xxyyyz,   \
                             ts_yyy_xxyyzz,   \
                             ts_yyy_xxyzzz,   \
                             ts_yyy_xyyyyy,   \
                             ts_yyy_xyyyyz,   \
                             ts_yyy_xyyyzz,   \
                             ts_yyy_xyyzzz,   \
                             ts_yyy_xyzzzz,   \
                             ts_yyy_yyyyyy,   \
                             ts_yyy_yyyyyz,   \
                             ts_yyy_yyyyzz,   \
                             ts_yyy_yyyzzz,   \
                             ts_yyy_yyzzzz,   \
                             ts_yyy_yzzzzz,   \
                             ts_yyy_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyy_xxxxxx[i] =
            -4.0 * ts_xxy_xxxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxxxxx[i] * fe_0 + tk_xxyy_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxyyy_xxxxxx[i] * fz_0;

        tk_xxyyy_xxxxxy[i] = -2.0 * ts_yyy_xxxxxy[i] * fbe_0 * fz_0 + tk_yyy_xxxxxy[i] * fe_0 + 5.0 * tk_xyyy_xxxxy[i] * fe_0 +
                             tk_xyyy_xxxxxy[i] * pa_x[i] + 2.0 * ts_xxyyy_xxxxxy[i] * fz_0;

        tk_xxyyy_xxxxxz[i] =
            -4.0 * ts_xxy_xxxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxxxxz[i] * fe_0 + tk_xxyy_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxyyy_xxxxxz[i] * fz_0;

        tk_xxyyy_xxxxyy[i] = -2.0 * ts_yyy_xxxxyy[i] * fbe_0 * fz_0 + tk_yyy_xxxxyy[i] * fe_0 + 4.0 * tk_xyyy_xxxyy[i] * fe_0 +
                             tk_xyyy_xxxxyy[i] * pa_x[i] + 2.0 * ts_xxyyy_xxxxyy[i] * fz_0;

        tk_xxyyy_xxxxyz[i] = -2.0 * ts_yyy_xxxxyz[i] * fbe_0 * fz_0 + tk_yyy_xxxxyz[i] * fe_0 + 4.0 * tk_xyyy_xxxyz[i] * fe_0 +
                             tk_xyyy_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxyyy_xxxxyz[i] * fz_0;

        tk_xxyyy_xxxxzz[i] =
            -4.0 * ts_xxy_xxxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxxxzz[i] * fe_0 + tk_xxyy_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxyyy_xxxxzz[i] * fz_0;

        tk_xxyyy_xxxyyy[i] = -2.0 * ts_yyy_xxxyyy[i] * fbe_0 * fz_0 + tk_yyy_xxxyyy[i] * fe_0 + 3.0 * tk_xyyy_xxyyy[i] * fe_0 +
                             tk_xyyy_xxxyyy[i] * pa_x[i] + 2.0 * ts_xxyyy_xxxyyy[i] * fz_0;

        tk_xxyyy_xxxyyz[i] = -2.0 * ts_yyy_xxxyyz[i] * fbe_0 * fz_0 + tk_yyy_xxxyyz[i] * fe_0 + 3.0 * tk_xyyy_xxyyz[i] * fe_0 +
                             tk_xyyy_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxyyy_xxxyyz[i] * fz_0;

        tk_xxyyy_xxxyzz[i] = -2.0 * ts_yyy_xxxyzz[i] * fbe_0 * fz_0 + tk_yyy_xxxyzz[i] * fe_0 + 3.0 * tk_xyyy_xxyzz[i] * fe_0 +
                             tk_xyyy_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxyyy_xxxyzz[i] * fz_0;

        tk_xxyyy_xxxzzz[i] =
            -4.0 * ts_xxy_xxxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxxzzz[i] * fe_0 + tk_xxyy_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxyyy_xxxzzz[i] * fz_0;

        tk_xxyyy_xxyyyy[i] = -2.0 * ts_yyy_xxyyyy[i] * fbe_0 * fz_0 + tk_yyy_xxyyyy[i] * fe_0 + 2.0 * tk_xyyy_xyyyy[i] * fe_0 +
                             tk_xyyy_xxyyyy[i] * pa_x[i] + 2.0 * ts_xxyyy_xxyyyy[i] * fz_0;

        tk_xxyyy_xxyyyz[i] = -2.0 * ts_yyy_xxyyyz[i] * fbe_0 * fz_0 + tk_yyy_xxyyyz[i] * fe_0 + 2.0 * tk_xyyy_xyyyz[i] * fe_0 +
                             tk_xyyy_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxyyy_xxyyyz[i] * fz_0;

        tk_xxyyy_xxyyzz[i] = -2.0 * ts_yyy_xxyyzz[i] * fbe_0 * fz_0 + tk_yyy_xxyyzz[i] * fe_0 + 2.0 * tk_xyyy_xyyzz[i] * fe_0 +
                             tk_xyyy_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxyyy_xxyyzz[i] * fz_0;

        tk_xxyyy_xxyzzz[i] = -2.0 * ts_yyy_xxyzzz[i] * fbe_0 * fz_0 + tk_yyy_xxyzzz[i] * fe_0 + 2.0 * tk_xyyy_xyzzz[i] * fe_0 +
                             tk_xyyy_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxyyy_xxyzzz[i] * fz_0;

        tk_xxyyy_xxzzzz[i] =
            -4.0 * ts_xxy_xxzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xxzzzz[i] * fe_0 + tk_xxyy_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxyyy_xxzzzz[i] * fz_0;

        tk_xxyyy_xyyyyy[i] = -2.0 * ts_yyy_xyyyyy[i] * fbe_0 * fz_0 + tk_yyy_xyyyyy[i] * fe_0 + tk_xyyy_yyyyy[i] * fe_0 +
                             tk_xyyy_xyyyyy[i] * pa_x[i] + 2.0 * ts_xxyyy_xyyyyy[i] * fz_0;

        tk_xxyyy_xyyyyz[i] = -2.0 * ts_yyy_xyyyyz[i] * fbe_0 * fz_0 + tk_yyy_xyyyyz[i] * fe_0 + tk_xyyy_yyyyz[i] * fe_0 +
                             tk_xyyy_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxyyy_xyyyyz[i] * fz_0;

        tk_xxyyy_xyyyzz[i] = -2.0 * ts_yyy_xyyyzz[i] * fbe_0 * fz_0 + tk_yyy_xyyyzz[i] * fe_0 + tk_xyyy_yyyzz[i] * fe_0 +
                             tk_xyyy_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxyyy_xyyyzz[i] * fz_0;

        tk_xxyyy_xyyzzz[i] = -2.0 * ts_yyy_xyyzzz[i] * fbe_0 * fz_0 + tk_yyy_xyyzzz[i] * fe_0 + tk_xyyy_yyzzz[i] * fe_0 +
                             tk_xyyy_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxyyy_xyyzzz[i] * fz_0;

        tk_xxyyy_xyzzzz[i] = -2.0 * ts_yyy_xyzzzz[i] * fbe_0 * fz_0 + tk_yyy_xyzzzz[i] * fe_0 + tk_xyyy_yzzzz[i] * fe_0 +
                             tk_xyyy_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxyyy_xyzzzz[i] * fz_0;

        tk_xxyyy_xzzzzz[i] =
            -4.0 * ts_xxy_xzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xzzzzz[i] * fe_0 + tk_xxyy_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxyyy_xzzzzz[i] * fz_0;

        tk_xxyyy_yyyyyy[i] =
            -2.0 * ts_yyy_yyyyyy[i] * fbe_0 * fz_0 + tk_yyy_yyyyyy[i] * fe_0 + tk_xyyy_yyyyyy[i] * pa_x[i] + 2.0 * ts_xxyyy_yyyyyy[i] * fz_0;

        tk_xxyyy_yyyyyz[i] =
            -2.0 * ts_yyy_yyyyyz[i] * fbe_0 * fz_0 + tk_yyy_yyyyyz[i] * fe_0 + tk_xyyy_yyyyyz[i] * pa_x[i] + 2.0 * ts_xxyyy_yyyyyz[i] * fz_0;

        tk_xxyyy_yyyyzz[i] =
            -2.0 * ts_yyy_yyyyzz[i] * fbe_0 * fz_0 + tk_yyy_yyyyzz[i] * fe_0 + tk_xyyy_yyyyzz[i] * pa_x[i] + 2.0 * ts_xxyyy_yyyyzz[i] * fz_0;

        tk_xxyyy_yyyzzz[i] =
            -2.0 * ts_yyy_yyyzzz[i] * fbe_0 * fz_0 + tk_yyy_yyyzzz[i] * fe_0 + tk_xyyy_yyyzzz[i] * pa_x[i] + 2.0 * ts_xxyyy_yyyzzz[i] * fz_0;

        tk_xxyyy_yyzzzz[i] =
            -2.0 * ts_yyy_yyzzzz[i] * fbe_0 * fz_0 + tk_yyy_yyzzzz[i] * fe_0 + tk_xyyy_yyzzzz[i] * pa_x[i] + 2.0 * ts_xxyyy_yyzzzz[i] * fz_0;

        tk_xxyyy_yzzzzz[i] =
            -2.0 * ts_yyy_yzzzzz[i] * fbe_0 * fz_0 + tk_yyy_yzzzzz[i] * fe_0 + tk_xyyy_yzzzzz[i] * pa_x[i] + 2.0 * ts_xxyyy_yzzzzz[i] * fz_0;

        tk_xxyyy_zzzzzz[i] =
            -2.0 * ts_yyy_zzzzzz[i] * fbe_0 * fz_0 + tk_yyy_zzzzzz[i] * fe_0 + tk_xyyy_zzzzzz[i] * pa_x[i] + 2.0 * ts_xxyyy_zzzzzz[i] * fz_0;
    }

    // Set up 196-224 components of targeted buffer : HI

    auto tk_xxyyz_xxxxxx = pbuffer.data(idx_kin_hi + 196);

    auto tk_xxyyz_xxxxxy = pbuffer.data(idx_kin_hi + 197);

    auto tk_xxyyz_xxxxxz = pbuffer.data(idx_kin_hi + 198);

    auto tk_xxyyz_xxxxyy = pbuffer.data(idx_kin_hi + 199);

    auto tk_xxyyz_xxxxyz = pbuffer.data(idx_kin_hi + 200);

    auto tk_xxyyz_xxxxzz = pbuffer.data(idx_kin_hi + 201);

    auto tk_xxyyz_xxxyyy = pbuffer.data(idx_kin_hi + 202);

    auto tk_xxyyz_xxxyyz = pbuffer.data(idx_kin_hi + 203);

    auto tk_xxyyz_xxxyzz = pbuffer.data(idx_kin_hi + 204);

    auto tk_xxyyz_xxxzzz = pbuffer.data(idx_kin_hi + 205);

    auto tk_xxyyz_xxyyyy = pbuffer.data(idx_kin_hi + 206);

    auto tk_xxyyz_xxyyyz = pbuffer.data(idx_kin_hi + 207);

    auto tk_xxyyz_xxyyzz = pbuffer.data(idx_kin_hi + 208);

    auto tk_xxyyz_xxyzzz = pbuffer.data(idx_kin_hi + 209);

    auto tk_xxyyz_xxzzzz = pbuffer.data(idx_kin_hi + 210);

    auto tk_xxyyz_xyyyyy = pbuffer.data(idx_kin_hi + 211);

    auto tk_xxyyz_xyyyyz = pbuffer.data(idx_kin_hi + 212);

    auto tk_xxyyz_xyyyzz = pbuffer.data(idx_kin_hi + 213);

    auto tk_xxyyz_xyyzzz = pbuffer.data(idx_kin_hi + 214);

    auto tk_xxyyz_xyzzzz = pbuffer.data(idx_kin_hi + 215);

    auto tk_xxyyz_xzzzzz = pbuffer.data(idx_kin_hi + 216);

    auto tk_xxyyz_yyyyyy = pbuffer.data(idx_kin_hi + 217);

    auto tk_xxyyz_yyyyyz = pbuffer.data(idx_kin_hi + 218);

    auto tk_xxyyz_yyyyzz = pbuffer.data(idx_kin_hi + 219);

    auto tk_xxyyz_yyyzzz = pbuffer.data(idx_kin_hi + 220);

    auto tk_xxyyz_yyzzzz = pbuffer.data(idx_kin_hi + 221);

    auto tk_xxyyz_yzzzzz = pbuffer.data(idx_kin_hi + 222);

    auto tk_xxyyz_zzzzzz = pbuffer.data(idx_kin_hi + 223);

#pragma omp simd aligned(pa_z,                \
                             tk_xxyy_xxxxx,   \
                             tk_xxyy_xxxxxx,  \
                             tk_xxyy_xxxxxy,  \
                             tk_xxyy_xxxxxz,  \
                             tk_xxyy_xxxxy,   \
                             tk_xxyy_xxxxyy,  \
                             tk_xxyy_xxxxyz,  \
                             tk_xxyy_xxxxz,   \
                             tk_xxyy_xxxxzz,  \
                             tk_xxyy_xxxyy,   \
                             tk_xxyy_xxxyyy,  \
                             tk_xxyy_xxxyyz,  \
                             tk_xxyy_xxxyz,   \
                             tk_xxyy_xxxyzz,  \
                             tk_xxyy_xxxzz,   \
                             tk_xxyy_xxxzzz,  \
                             tk_xxyy_xxyyy,   \
                             tk_xxyy_xxyyyy,  \
                             tk_xxyy_xxyyyz,  \
                             tk_xxyy_xxyyz,   \
                             tk_xxyy_xxyyzz,  \
                             tk_xxyy_xxyzz,   \
                             tk_xxyy_xxyzzz,  \
                             tk_xxyy_xxzzz,   \
                             tk_xxyy_xxzzzz,  \
                             tk_xxyy_xyyyy,   \
                             tk_xxyy_xyyyyy,  \
                             tk_xxyy_xyyyyz,  \
                             tk_xxyy_xyyyz,   \
                             tk_xxyy_xyyyzz,  \
                             tk_xxyy_xyyzz,   \
                             tk_xxyy_xyyzzz,  \
                             tk_xxyy_xyzzz,   \
                             tk_xxyy_xyzzzz,  \
                             tk_xxyy_xzzzz,   \
                             tk_xxyy_xzzzzz,  \
                             tk_xxyy_yyyyy,   \
                             tk_xxyy_yyyyyy,  \
                             tk_xxyy_yyyyyz,  \
                             tk_xxyy_yyyyz,   \
                             tk_xxyy_yyyyzz,  \
                             tk_xxyy_yyyzz,   \
                             tk_xxyy_yyyzzz,  \
                             tk_xxyy_yyzzz,   \
                             tk_xxyy_yyzzzz,  \
                             tk_xxyy_yzzzz,   \
                             tk_xxyy_yzzzzz,  \
                             tk_xxyy_zzzzz,   \
                             tk_xxyy_zzzzzz,  \
                             tk_xxyyz_xxxxxx, \
                             tk_xxyyz_xxxxxy, \
                             tk_xxyyz_xxxxxz, \
                             tk_xxyyz_xxxxyy, \
                             tk_xxyyz_xxxxyz, \
                             tk_xxyyz_xxxxzz, \
                             tk_xxyyz_xxxyyy, \
                             tk_xxyyz_xxxyyz, \
                             tk_xxyyz_xxxyzz, \
                             tk_xxyyz_xxxzzz, \
                             tk_xxyyz_xxyyyy, \
                             tk_xxyyz_xxyyyz, \
                             tk_xxyyz_xxyyzz, \
                             tk_xxyyz_xxyzzz, \
                             tk_xxyyz_xxzzzz, \
                             tk_xxyyz_xyyyyy, \
                             tk_xxyyz_xyyyyz, \
                             tk_xxyyz_xyyyzz, \
                             tk_xxyyz_xyyzzz, \
                             tk_xxyyz_xyzzzz, \
                             tk_xxyyz_xzzzzz, \
                             tk_xxyyz_yyyyyy, \
                             tk_xxyyz_yyyyyz, \
                             tk_xxyyz_yyyyzz, \
                             tk_xxyyz_yyyzzz, \
                             tk_xxyyz_yyzzzz, \
                             tk_xxyyz_yzzzzz, \
                             tk_xxyyz_zzzzzz, \
                             ts_xxyyz_xxxxxx, \
                             ts_xxyyz_xxxxxy, \
                             ts_xxyyz_xxxxxz, \
                             ts_xxyyz_xxxxyy, \
                             ts_xxyyz_xxxxyz, \
                             ts_xxyyz_xxxxzz, \
                             ts_xxyyz_xxxyyy, \
                             ts_xxyyz_xxxyyz, \
                             ts_xxyyz_xxxyzz, \
                             ts_xxyyz_xxxzzz, \
                             ts_xxyyz_xxyyyy, \
                             ts_xxyyz_xxyyyz, \
                             ts_xxyyz_xxyyzz, \
                             ts_xxyyz_xxyzzz, \
                             ts_xxyyz_xxzzzz, \
                             ts_xxyyz_xyyyyy, \
                             ts_xxyyz_xyyyyz, \
                             ts_xxyyz_xyyyzz, \
                             ts_xxyyz_xyyzzz, \
                             ts_xxyyz_xyzzzz, \
                             ts_xxyyz_xzzzzz, \
                             ts_xxyyz_yyyyyy, \
                             ts_xxyyz_yyyyyz, \
                             ts_xxyyz_yyyyzz, \
                             ts_xxyyz_yyyzzz, \
                             ts_xxyyz_yyzzzz, \
                             ts_xxyyz_yzzzzz, \
                             ts_xxyyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyyz_xxxxxx[i] = tk_xxyy_xxxxxx[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxxxx[i] * fz_0;

        tk_xxyyz_xxxxxy[i] = tk_xxyy_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxxxy[i] * fz_0;

        tk_xxyyz_xxxxxz[i] = tk_xxyy_xxxxx[i] * fe_0 + tk_xxyy_xxxxxz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxxxz[i] * fz_0;

        tk_xxyyz_xxxxyy[i] = tk_xxyy_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxxyy[i] * fz_0;

        tk_xxyyz_xxxxyz[i] = tk_xxyy_xxxxy[i] * fe_0 + tk_xxyy_xxxxyz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxxyz[i] * fz_0;

        tk_xxyyz_xxxxzz[i] = 2.0 * tk_xxyy_xxxxz[i] * fe_0 + tk_xxyy_xxxxzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxxzz[i] * fz_0;

        tk_xxyyz_xxxyyy[i] = tk_xxyy_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxyyy[i] * fz_0;

        tk_xxyyz_xxxyyz[i] = tk_xxyy_xxxyy[i] * fe_0 + tk_xxyy_xxxyyz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxyyz[i] * fz_0;

        tk_xxyyz_xxxyzz[i] = 2.0 * tk_xxyy_xxxyz[i] * fe_0 + tk_xxyy_xxxyzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxyzz[i] * fz_0;

        tk_xxyyz_xxxzzz[i] = 3.0 * tk_xxyy_xxxzz[i] * fe_0 + tk_xxyy_xxxzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxxzzz[i] * fz_0;

        tk_xxyyz_xxyyyy[i] = tk_xxyy_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxyyz_xxyyyy[i] * fz_0;

        tk_xxyyz_xxyyyz[i] = tk_xxyy_xxyyy[i] * fe_0 + tk_xxyy_xxyyyz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxyyyz[i] * fz_0;

        tk_xxyyz_xxyyzz[i] = 2.0 * tk_xxyy_xxyyz[i] * fe_0 + tk_xxyy_xxyyzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxyyzz[i] * fz_0;

        tk_xxyyz_xxyzzz[i] = 3.0 * tk_xxyy_xxyzz[i] * fe_0 + tk_xxyy_xxyzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxyzzz[i] * fz_0;

        tk_xxyyz_xxzzzz[i] = 4.0 * tk_xxyy_xxzzz[i] * fe_0 + tk_xxyy_xxzzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xxzzzz[i] * fz_0;

        tk_xxyyz_xyyyyy[i] = tk_xxyy_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxyyz_xyyyyy[i] * fz_0;

        tk_xxyyz_xyyyyz[i] = tk_xxyy_xyyyy[i] * fe_0 + tk_xxyy_xyyyyz[i] * pa_z[i] + 2.0 * ts_xxyyz_xyyyyz[i] * fz_0;

        tk_xxyyz_xyyyzz[i] = 2.0 * tk_xxyy_xyyyz[i] * fe_0 + tk_xxyy_xyyyzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xyyyzz[i] * fz_0;

        tk_xxyyz_xyyzzz[i] = 3.0 * tk_xxyy_xyyzz[i] * fe_0 + tk_xxyy_xyyzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xyyzzz[i] * fz_0;

        tk_xxyyz_xyzzzz[i] = 4.0 * tk_xxyy_xyzzz[i] * fe_0 + tk_xxyy_xyzzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xyzzzz[i] * fz_0;

        tk_xxyyz_xzzzzz[i] = 5.0 * tk_xxyy_xzzzz[i] * fe_0 + tk_xxyy_xzzzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_xzzzzz[i] * fz_0;

        tk_xxyyz_yyyyyy[i] = tk_xxyy_yyyyyy[i] * pa_z[i] + 2.0 * ts_xxyyz_yyyyyy[i] * fz_0;

        tk_xxyyz_yyyyyz[i] = tk_xxyy_yyyyy[i] * fe_0 + tk_xxyy_yyyyyz[i] * pa_z[i] + 2.0 * ts_xxyyz_yyyyyz[i] * fz_0;

        tk_xxyyz_yyyyzz[i] = 2.0 * tk_xxyy_yyyyz[i] * fe_0 + tk_xxyy_yyyyzz[i] * pa_z[i] + 2.0 * ts_xxyyz_yyyyzz[i] * fz_0;

        tk_xxyyz_yyyzzz[i] = 3.0 * tk_xxyy_yyyzz[i] * fe_0 + tk_xxyy_yyyzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_yyyzzz[i] * fz_0;

        tk_xxyyz_yyzzzz[i] = 4.0 * tk_xxyy_yyzzz[i] * fe_0 + tk_xxyy_yyzzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_yyzzzz[i] * fz_0;

        tk_xxyyz_yzzzzz[i] = 5.0 * tk_xxyy_yzzzz[i] * fe_0 + tk_xxyy_yzzzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_yzzzzz[i] * fz_0;

        tk_xxyyz_zzzzzz[i] = 6.0 * tk_xxyy_zzzzz[i] * fe_0 + tk_xxyy_zzzzzz[i] * pa_z[i] + 2.0 * ts_xxyyz_zzzzzz[i] * fz_0;
    }

    // Set up 224-252 components of targeted buffer : HI

    auto tk_xxyzz_xxxxxx = pbuffer.data(idx_kin_hi + 224);

    auto tk_xxyzz_xxxxxy = pbuffer.data(idx_kin_hi + 225);

    auto tk_xxyzz_xxxxxz = pbuffer.data(idx_kin_hi + 226);

    auto tk_xxyzz_xxxxyy = pbuffer.data(idx_kin_hi + 227);

    auto tk_xxyzz_xxxxyz = pbuffer.data(idx_kin_hi + 228);

    auto tk_xxyzz_xxxxzz = pbuffer.data(idx_kin_hi + 229);

    auto tk_xxyzz_xxxyyy = pbuffer.data(idx_kin_hi + 230);

    auto tk_xxyzz_xxxyyz = pbuffer.data(idx_kin_hi + 231);

    auto tk_xxyzz_xxxyzz = pbuffer.data(idx_kin_hi + 232);

    auto tk_xxyzz_xxxzzz = pbuffer.data(idx_kin_hi + 233);

    auto tk_xxyzz_xxyyyy = pbuffer.data(idx_kin_hi + 234);

    auto tk_xxyzz_xxyyyz = pbuffer.data(idx_kin_hi + 235);

    auto tk_xxyzz_xxyyzz = pbuffer.data(idx_kin_hi + 236);

    auto tk_xxyzz_xxyzzz = pbuffer.data(idx_kin_hi + 237);

    auto tk_xxyzz_xxzzzz = pbuffer.data(idx_kin_hi + 238);

    auto tk_xxyzz_xyyyyy = pbuffer.data(idx_kin_hi + 239);

    auto tk_xxyzz_xyyyyz = pbuffer.data(idx_kin_hi + 240);

    auto tk_xxyzz_xyyyzz = pbuffer.data(idx_kin_hi + 241);

    auto tk_xxyzz_xyyzzz = pbuffer.data(idx_kin_hi + 242);

    auto tk_xxyzz_xyzzzz = pbuffer.data(idx_kin_hi + 243);

    auto tk_xxyzz_xzzzzz = pbuffer.data(idx_kin_hi + 244);

    auto tk_xxyzz_yyyyyy = pbuffer.data(idx_kin_hi + 245);

    auto tk_xxyzz_yyyyyz = pbuffer.data(idx_kin_hi + 246);

    auto tk_xxyzz_yyyyzz = pbuffer.data(idx_kin_hi + 247);

    auto tk_xxyzz_yyyzzz = pbuffer.data(idx_kin_hi + 248);

    auto tk_xxyzz_yyzzzz = pbuffer.data(idx_kin_hi + 249);

    auto tk_xxyzz_yzzzzz = pbuffer.data(idx_kin_hi + 250);

    auto tk_xxyzz_zzzzzz = pbuffer.data(idx_kin_hi + 251);

#pragma omp simd aligned(pa_y,                \
                             tk_xxyzz_xxxxxx, \
                             tk_xxyzz_xxxxxy, \
                             tk_xxyzz_xxxxxz, \
                             tk_xxyzz_xxxxyy, \
                             tk_xxyzz_xxxxyz, \
                             tk_xxyzz_xxxxzz, \
                             tk_xxyzz_xxxyyy, \
                             tk_xxyzz_xxxyyz, \
                             tk_xxyzz_xxxyzz, \
                             tk_xxyzz_xxxzzz, \
                             tk_xxyzz_xxyyyy, \
                             tk_xxyzz_xxyyyz, \
                             tk_xxyzz_xxyyzz, \
                             tk_xxyzz_xxyzzz, \
                             tk_xxyzz_xxzzzz, \
                             tk_xxyzz_xyyyyy, \
                             tk_xxyzz_xyyyyz, \
                             tk_xxyzz_xyyyzz, \
                             tk_xxyzz_xyyzzz, \
                             tk_xxyzz_xyzzzz, \
                             tk_xxyzz_xzzzzz, \
                             tk_xxyzz_yyyyyy, \
                             tk_xxyzz_yyyyyz, \
                             tk_xxyzz_yyyyzz, \
                             tk_xxyzz_yyyzzz, \
                             tk_xxyzz_yyzzzz, \
                             tk_xxyzz_yzzzzz, \
                             tk_xxyzz_zzzzzz, \
                             tk_xxzz_xxxxx,   \
                             tk_xxzz_xxxxxx,  \
                             tk_xxzz_xxxxxy,  \
                             tk_xxzz_xxxxxz,  \
                             tk_xxzz_xxxxy,   \
                             tk_xxzz_xxxxyy,  \
                             tk_xxzz_xxxxyz,  \
                             tk_xxzz_xxxxz,   \
                             tk_xxzz_xxxxzz,  \
                             tk_xxzz_xxxyy,   \
                             tk_xxzz_xxxyyy,  \
                             tk_xxzz_xxxyyz,  \
                             tk_xxzz_xxxyz,   \
                             tk_xxzz_xxxyzz,  \
                             tk_xxzz_xxxzz,   \
                             tk_xxzz_xxxzzz,  \
                             tk_xxzz_xxyyy,   \
                             tk_xxzz_xxyyyy,  \
                             tk_xxzz_xxyyyz,  \
                             tk_xxzz_xxyyz,   \
                             tk_xxzz_xxyyzz,  \
                             tk_xxzz_xxyzz,   \
                             tk_xxzz_xxyzzz,  \
                             tk_xxzz_xxzzz,   \
                             tk_xxzz_xxzzzz,  \
                             tk_xxzz_xyyyy,   \
                             tk_xxzz_xyyyyy,  \
                             tk_xxzz_xyyyyz,  \
                             tk_xxzz_xyyyz,   \
                             tk_xxzz_xyyyzz,  \
                             tk_xxzz_xyyzz,   \
                             tk_xxzz_xyyzzz,  \
                             tk_xxzz_xyzzz,   \
                             tk_xxzz_xyzzzz,  \
                             tk_xxzz_xzzzz,   \
                             tk_xxzz_xzzzzz,  \
                             tk_xxzz_yyyyy,   \
                             tk_xxzz_yyyyyy,  \
                             tk_xxzz_yyyyyz,  \
                             tk_xxzz_yyyyz,   \
                             tk_xxzz_yyyyzz,  \
                             tk_xxzz_yyyzz,   \
                             tk_xxzz_yyyzzz,  \
                             tk_xxzz_yyzzz,   \
                             tk_xxzz_yyzzzz,  \
                             tk_xxzz_yzzzz,   \
                             tk_xxzz_yzzzzz,  \
                             tk_xxzz_zzzzz,   \
                             tk_xxzz_zzzzzz,  \
                             ts_xxyzz_xxxxxx, \
                             ts_xxyzz_xxxxxy, \
                             ts_xxyzz_xxxxxz, \
                             ts_xxyzz_xxxxyy, \
                             ts_xxyzz_xxxxyz, \
                             ts_xxyzz_xxxxzz, \
                             ts_xxyzz_xxxyyy, \
                             ts_xxyzz_xxxyyz, \
                             ts_xxyzz_xxxyzz, \
                             ts_xxyzz_xxxzzz, \
                             ts_xxyzz_xxyyyy, \
                             ts_xxyzz_xxyyyz, \
                             ts_xxyzz_xxyyzz, \
                             ts_xxyzz_xxyzzz, \
                             ts_xxyzz_xxzzzz, \
                             ts_xxyzz_xyyyyy, \
                             ts_xxyzz_xyyyyz, \
                             ts_xxyzz_xyyyzz, \
                             ts_xxyzz_xyyzzz, \
                             ts_xxyzz_xyzzzz, \
                             ts_xxyzz_xzzzzz, \
                             ts_xxyzz_yyyyyy, \
                             ts_xxyzz_yyyyyz, \
                             ts_xxyzz_yyyyzz, \
                             ts_xxyzz_yyyzzz, \
                             ts_xxyzz_yyzzzz, \
                             ts_xxyzz_yzzzzz, \
                             ts_xxyzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyzz_xxxxxx[i] = tk_xxzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxxxx[i] * fz_0;

        tk_xxyzz_xxxxxy[i] = tk_xxzz_xxxxx[i] * fe_0 + tk_xxzz_xxxxxy[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxxxy[i] * fz_0;

        tk_xxyzz_xxxxxz[i] = tk_xxzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxxxz[i] * fz_0;

        tk_xxyzz_xxxxyy[i] = 2.0 * tk_xxzz_xxxxy[i] * fe_0 + tk_xxzz_xxxxyy[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxxyy[i] * fz_0;

        tk_xxyzz_xxxxyz[i] = tk_xxzz_xxxxz[i] * fe_0 + tk_xxzz_xxxxyz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxxyz[i] * fz_0;

        tk_xxyzz_xxxxzz[i] = tk_xxzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxxzz[i] * fz_0;

        tk_xxyzz_xxxyyy[i] = 3.0 * tk_xxzz_xxxyy[i] * fe_0 + tk_xxzz_xxxyyy[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxyyy[i] * fz_0;

        tk_xxyzz_xxxyyz[i] = 2.0 * tk_xxzz_xxxyz[i] * fe_0 + tk_xxzz_xxxyyz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxyyz[i] * fz_0;

        tk_xxyzz_xxxyzz[i] = tk_xxzz_xxxzz[i] * fe_0 + tk_xxzz_xxxyzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxyzz[i] * fz_0;

        tk_xxyzz_xxxzzz[i] = tk_xxzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxxzzz[i] * fz_0;

        tk_xxyzz_xxyyyy[i] = 4.0 * tk_xxzz_xxyyy[i] * fe_0 + tk_xxzz_xxyyyy[i] * pa_y[i] + 2.0 * ts_xxyzz_xxyyyy[i] * fz_0;

        tk_xxyzz_xxyyyz[i] = 3.0 * tk_xxzz_xxyyz[i] * fe_0 + tk_xxzz_xxyyyz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxyyyz[i] * fz_0;

        tk_xxyzz_xxyyzz[i] = 2.0 * tk_xxzz_xxyzz[i] * fe_0 + tk_xxzz_xxyyzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxyyzz[i] * fz_0;

        tk_xxyzz_xxyzzz[i] = tk_xxzz_xxzzz[i] * fe_0 + tk_xxzz_xxyzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxyzzz[i] * fz_0;

        tk_xxyzz_xxzzzz[i] = tk_xxzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xxzzzz[i] * fz_0;

        tk_xxyzz_xyyyyy[i] = 5.0 * tk_xxzz_xyyyy[i] * fe_0 + tk_xxzz_xyyyyy[i] * pa_y[i] + 2.0 * ts_xxyzz_xyyyyy[i] * fz_0;

        tk_xxyzz_xyyyyz[i] = 4.0 * tk_xxzz_xyyyz[i] * fe_0 + tk_xxzz_xyyyyz[i] * pa_y[i] + 2.0 * ts_xxyzz_xyyyyz[i] * fz_0;

        tk_xxyzz_xyyyzz[i] = 3.0 * tk_xxzz_xyyzz[i] * fe_0 + tk_xxzz_xyyyzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xyyyzz[i] * fz_0;

        tk_xxyzz_xyyzzz[i] = 2.0 * tk_xxzz_xyzzz[i] * fe_0 + tk_xxzz_xyyzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xyyzzz[i] * fz_0;

        tk_xxyzz_xyzzzz[i] = tk_xxzz_xzzzz[i] * fe_0 + tk_xxzz_xyzzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xyzzzz[i] * fz_0;

        tk_xxyzz_xzzzzz[i] = tk_xxzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_xzzzzz[i] * fz_0;

        tk_xxyzz_yyyyyy[i] = 6.0 * tk_xxzz_yyyyy[i] * fe_0 + tk_xxzz_yyyyyy[i] * pa_y[i] + 2.0 * ts_xxyzz_yyyyyy[i] * fz_0;

        tk_xxyzz_yyyyyz[i] = 5.0 * tk_xxzz_yyyyz[i] * fe_0 + tk_xxzz_yyyyyz[i] * pa_y[i] + 2.0 * ts_xxyzz_yyyyyz[i] * fz_0;

        tk_xxyzz_yyyyzz[i] = 4.0 * tk_xxzz_yyyzz[i] * fe_0 + tk_xxzz_yyyyzz[i] * pa_y[i] + 2.0 * ts_xxyzz_yyyyzz[i] * fz_0;

        tk_xxyzz_yyyzzz[i] = 3.0 * tk_xxzz_yyzzz[i] * fe_0 + tk_xxzz_yyyzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_yyyzzz[i] * fz_0;

        tk_xxyzz_yyzzzz[i] = 2.0 * tk_xxzz_yzzzz[i] * fe_0 + tk_xxzz_yyzzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_yyzzzz[i] * fz_0;

        tk_xxyzz_yzzzzz[i] = tk_xxzz_zzzzz[i] * fe_0 + tk_xxzz_yzzzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_yzzzzz[i] * fz_0;

        tk_xxyzz_zzzzzz[i] = tk_xxzz_zzzzzz[i] * pa_y[i] + 2.0 * ts_xxyzz_zzzzzz[i] * fz_0;
    }

    // Set up 252-280 components of targeted buffer : HI

    auto tk_xxzzz_xxxxxx = pbuffer.data(idx_kin_hi + 252);

    auto tk_xxzzz_xxxxxy = pbuffer.data(idx_kin_hi + 253);

    auto tk_xxzzz_xxxxxz = pbuffer.data(idx_kin_hi + 254);

    auto tk_xxzzz_xxxxyy = pbuffer.data(idx_kin_hi + 255);

    auto tk_xxzzz_xxxxyz = pbuffer.data(idx_kin_hi + 256);

    auto tk_xxzzz_xxxxzz = pbuffer.data(idx_kin_hi + 257);

    auto tk_xxzzz_xxxyyy = pbuffer.data(idx_kin_hi + 258);

    auto tk_xxzzz_xxxyyz = pbuffer.data(idx_kin_hi + 259);

    auto tk_xxzzz_xxxyzz = pbuffer.data(idx_kin_hi + 260);

    auto tk_xxzzz_xxxzzz = pbuffer.data(idx_kin_hi + 261);

    auto tk_xxzzz_xxyyyy = pbuffer.data(idx_kin_hi + 262);

    auto tk_xxzzz_xxyyyz = pbuffer.data(idx_kin_hi + 263);

    auto tk_xxzzz_xxyyzz = pbuffer.data(idx_kin_hi + 264);

    auto tk_xxzzz_xxyzzz = pbuffer.data(idx_kin_hi + 265);

    auto tk_xxzzz_xxzzzz = pbuffer.data(idx_kin_hi + 266);

    auto tk_xxzzz_xyyyyy = pbuffer.data(idx_kin_hi + 267);

    auto tk_xxzzz_xyyyyz = pbuffer.data(idx_kin_hi + 268);

    auto tk_xxzzz_xyyyzz = pbuffer.data(idx_kin_hi + 269);

    auto tk_xxzzz_xyyzzz = pbuffer.data(idx_kin_hi + 270);

    auto tk_xxzzz_xyzzzz = pbuffer.data(idx_kin_hi + 271);

    auto tk_xxzzz_xzzzzz = pbuffer.data(idx_kin_hi + 272);

    auto tk_xxzzz_yyyyyy = pbuffer.data(idx_kin_hi + 273);

    auto tk_xxzzz_yyyyyz = pbuffer.data(idx_kin_hi + 274);

    auto tk_xxzzz_yyyyzz = pbuffer.data(idx_kin_hi + 275);

    auto tk_xxzzz_yyyzzz = pbuffer.data(idx_kin_hi + 276);

    auto tk_xxzzz_yyzzzz = pbuffer.data(idx_kin_hi + 277);

    auto tk_xxzzz_yzzzzz = pbuffer.data(idx_kin_hi + 278);

    auto tk_xxzzz_zzzzzz = pbuffer.data(idx_kin_hi + 279);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tk_xxz_xxxxxx,   \
                             tk_xxz_xxxxxy,   \
                             tk_xxz_xxxxyy,   \
                             tk_xxz_xxxyyy,   \
                             tk_xxz_xxyyyy,   \
                             tk_xxz_xyyyyy,   \
                             tk_xxzz_xxxxxx,  \
                             tk_xxzz_xxxxxy,  \
                             tk_xxzz_xxxxyy,  \
                             tk_xxzz_xxxyyy,  \
                             tk_xxzz_xxyyyy,  \
                             tk_xxzz_xyyyyy,  \
                             tk_xxzzz_xxxxxx, \
                             tk_xxzzz_xxxxxy, \
                             tk_xxzzz_xxxxxz, \
                             tk_xxzzz_xxxxyy, \
                             tk_xxzzz_xxxxyz, \
                             tk_xxzzz_xxxxzz, \
                             tk_xxzzz_xxxyyy, \
                             tk_xxzzz_xxxyyz, \
                             tk_xxzzz_xxxyzz, \
                             tk_xxzzz_xxxzzz, \
                             tk_xxzzz_xxyyyy, \
                             tk_xxzzz_xxyyyz, \
                             tk_xxzzz_xxyyzz, \
                             tk_xxzzz_xxyzzz, \
                             tk_xxzzz_xxzzzz, \
                             tk_xxzzz_xyyyyy, \
                             tk_xxzzz_xyyyyz, \
                             tk_xxzzz_xyyyzz, \
                             tk_xxzzz_xyyzzz, \
                             tk_xxzzz_xyzzzz, \
                             tk_xxzzz_xzzzzz, \
                             tk_xxzzz_yyyyyy, \
                             tk_xxzzz_yyyyyz, \
                             tk_xxzzz_yyyyzz, \
                             tk_xxzzz_yyyzzz, \
                             tk_xxzzz_yyzzzz, \
                             tk_xxzzz_yzzzzz, \
                             tk_xxzzz_zzzzzz, \
                             tk_xzzz_xxxxxz,  \
                             tk_xzzz_xxxxyz,  \
                             tk_xzzz_xxxxz,   \
                             tk_xzzz_xxxxzz,  \
                             tk_xzzz_xxxyyz,  \
                             tk_xzzz_xxxyz,   \
                             tk_xzzz_xxxyzz,  \
                             tk_xzzz_xxxzz,   \
                             tk_xzzz_xxxzzz,  \
                             tk_xzzz_xxyyyz,  \
                             tk_xzzz_xxyyz,   \
                             tk_xzzz_xxyyzz,  \
                             tk_xzzz_xxyzz,   \
                             tk_xzzz_xxyzzz,  \
                             tk_xzzz_xxzzz,   \
                             tk_xzzz_xxzzzz,  \
                             tk_xzzz_xyyyyz,  \
                             tk_xzzz_xyyyz,   \
                             tk_xzzz_xyyyzz,  \
                             tk_xzzz_xyyzz,   \
                             tk_xzzz_xyyzzz,  \
                             tk_xzzz_xyzzz,   \
                             tk_xzzz_xyzzzz,  \
                             tk_xzzz_xzzzz,   \
                             tk_xzzz_xzzzzz,  \
                             tk_xzzz_yyyyyy,  \
                             tk_xzzz_yyyyyz,  \
                             tk_xzzz_yyyyz,   \
                             tk_xzzz_yyyyzz,  \
                             tk_xzzz_yyyzz,   \
                             tk_xzzz_yyyzzz,  \
                             tk_xzzz_yyzzz,   \
                             tk_xzzz_yyzzzz,  \
                             tk_xzzz_yzzzz,   \
                             tk_xzzz_yzzzzz,  \
                             tk_xzzz_zzzzz,   \
                             tk_xzzz_zzzzzz,  \
                             tk_zzz_xxxxxz,   \
                             tk_zzz_xxxxyz,   \
                             tk_zzz_xxxxzz,   \
                             tk_zzz_xxxyyz,   \
                             tk_zzz_xxxyzz,   \
                             tk_zzz_xxxzzz,   \
                             tk_zzz_xxyyyz,   \
                             tk_zzz_xxyyzz,   \
                             tk_zzz_xxyzzz,   \
                             tk_zzz_xxzzzz,   \
                             tk_zzz_xyyyyz,   \
                             tk_zzz_xyyyzz,   \
                             tk_zzz_xyyzzz,   \
                             tk_zzz_xyzzzz,   \
                             tk_zzz_xzzzzz,   \
                             tk_zzz_yyyyyy,   \
                             tk_zzz_yyyyyz,   \
                             tk_zzz_yyyyzz,   \
                             tk_zzz_yyyzzz,   \
                             tk_zzz_yyzzzz,   \
                             tk_zzz_yzzzzz,   \
                             tk_zzz_zzzzzz,   \
                             ts_xxz_xxxxxx,   \
                             ts_xxz_xxxxxy,   \
                             ts_xxz_xxxxyy,   \
                             ts_xxz_xxxyyy,   \
                             ts_xxz_xxyyyy,   \
                             ts_xxz_xyyyyy,   \
                             ts_xxzzz_xxxxxx, \
                             ts_xxzzz_xxxxxy, \
                             ts_xxzzz_xxxxxz, \
                             ts_xxzzz_xxxxyy, \
                             ts_xxzzz_xxxxyz, \
                             ts_xxzzz_xxxxzz, \
                             ts_xxzzz_xxxyyy, \
                             ts_xxzzz_xxxyyz, \
                             ts_xxzzz_xxxyzz, \
                             ts_xxzzz_xxxzzz, \
                             ts_xxzzz_xxyyyy, \
                             ts_xxzzz_xxyyyz, \
                             ts_xxzzz_xxyyzz, \
                             ts_xxzzz_xxyzzz, \
                             ts_xxzzz_xxzzzz, \
                             ts_xxzzz_xyyyyy, \
                             ts_xxzzz_xyyyyz, \
                             ts_xxzzz_xyyyzz, \
                             ts_xxzzz_xyyzzz, \
                             ts_xxzzz_xyzzzz, \
                             ts_xxzzz_xzzzzz, \
                             ts_xxzzz_yyyyyy, \
                             ts_xxzzz_yyyyyz, \
                             ts_xxzzz_yyyyzz, \
                             ts_xxzzz_yyyzzz, \
                             ts_xxzzz_yyzzzz, \
                             ts_xxzzz_yzzzzz, \
                             ts_xxzzz_zzzzzz, \
                             ts_zzz_xxxxxz,   \
                             ts_zzz_xxxxyz,   \
                             ts_zzz_xxxxzz,   \
                             ts_zzz_xxxyyz,   \
                             ts_zzz_xxxyzz,   \
                             ts_zzz_xxxzzz,   \
                             ts_zzz_xxyyyz,   \
                             ts_zzz_xxyyzz,   \
                             ts_zzz_xxyzzz,   \
                             ts_zzz_xxzzzz,   \
                             ts_zzz_xyyyyz,   \
                             ts_zzz_xyyyzz,   \
                             ts_zzz_xyyzzz,   \
                             ts_zzz_xyzzzz,   \
                             ts_zzz_xzzzzz,   \
                             ts_zzz_yyyyyy,   \
                             ts_zzz_yyyyyz,   \
                             ts_zzz_yyyyzz,   \
                             ts_zzz_yyyzzz,   \
                             ts_zzz_yyzzzz,   \
                             ts_zzz_yzzzzz,   \
                             ts_zzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzzz_xxxxxx[i] =
            -4.0 * ts_xxz_xxxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxxxxx[i] * fe_0 + tk_xxzz_xxxxxx[i] * pa_z[i] + 2.0 * ts_xxzzz_xxxxxx[i] * fz_0;

        tk_xxzzz_xxxxxy[i] =
            -4.0 * ts_xxz_xxxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxxxxy[i] * fe_0 + tk_xxzz_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxzzz_xxxxxy[i] * fz_0;

        tk_xxzzz_xxxxxz[i] = -2.0 * ts_zzz_xxxxxz[i] * fbe_0 * fz_0 + tk_zzz_xxxxxz[i] * fe_0 + 5.0 * tk_xzzz_xxxxz[i] * fe_0 +
                             tk_xzzz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxxxxz[i] * fz_0;

        tk_xxzzz_xxxxyy[i] =
            -4.0 * ts_xxz_xxxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxxxyy[i] * fe_0 + tk_xxzz_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxzzz_xxxxyy[i] * fz_0;

        tk_xxzzz_xxxxyz[i] = -2.0 * ts_zzz_xxxxyz[i] * fbe_0 * fz_0 + tk_zzz_xxxxyz[i] * fe_0 + 4.0 * tk_xzzz_xxxyz[i] * fe_0 +
                             tk_xzzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxxxyz[i] * fz_0;

        tk_xxzzz_xxxxzz[i] = -2.0 * ts_zzz_xxxxzz[i] * fbe_0 * fz_0 + tk_zzz_xxxxzz[i] * fe_0 + 4.0 * tk_xzzz_xxxzz[i] * fe_0 +
                             tk_xzzz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxxxzz[i] * fz_0;

        tk_xxzzz_xxxyyy[i] =
            -4.0 * ts_xxz_xxxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxxyyy[i] * fe_0 + tk_xxzz_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxzzz_xxxyyy[i] * fz_0;

        tk_xxzzz_xxxyyz[i] = -2.0 * ts_zzz_xxxyyz[i] * fbe_0 * fz_0 + tk_zzz_xxxyyz[i] * fe_0 + 3.0 * tk_xzzz_xxyyz[i] * fe_0 +
                             tk_xzzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxxyyz[i] * fz_0;

        tk_xxzzz_xxxyzz[i] = -2.0 * ts_zzz_xxxyzz[i] * fbe_0 * fz_0 + tk_zzz_xxxyzz[i] * fe_0 + 3.0 * tk_xzzz_xxyzz[i] * fe_0 +
                             tk_xzzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxxyzz[i] * fz_0;

        tk_xxzzz_xxxzzz[i] = -2.0 * ts_zzz_xxxzzz[i] * fbe_0 * fz_0 + tk_zzz_xxxzzz[i] * fe_0 + 3.0 * tk_xzzz_xxzzz[i] * fe_0 +
                             tk_xzzz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxxzzz[i] * fz_0;

        tk_xxzzz_xxyyyy[i] =
            -4.0 * ts_xxz_xxyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xxyyyy[i] * fe_0 + tk_xxzz_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxzzz_xxyyyy[i] * fz_0;

        tk_xxzzz_xxyyyz[i] = -2.0 * ts_zzz_xxyyyz[i] * fbe_0 * fz_0 + tk_zzz_xxyyyz[i] * fe_0 + 2.0 * tk_xzzz_xyyyz[i] * fe_0 +
                             tk_xzzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxyyyz[i] * fz_0;

        tk_xxzzz_xxyyzz[i] = -2.0 * ts_zzz_xxyyzz[i] * fbe_0 * fz_0 + tk_zzz_xxyyzz[i] * fe_0 + 2.0 * tk_xzzz_xyyzz[i] * fe_0 +
                             tk_xzzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxyyzz[i] * fz_0;

        tk_xxzzz_xxyzzz[i] = -2.0 * ts_zzz_xxyzzz[i] * fbe_0 * fz_0 + tk_zzz_xxyzzz[i] * fe_0 + 2.0 * tk_xzzz_xyzzz[i] * fe_0 +
                             tk_xzzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxyzzz[i] * fz_0;

        tk_xxzzz_xxzzzz[i] = -2.0 * ts_zzz_xxzzzz[i] * fbe_0 * fz_0 + tk_zzz_xxzzzz[i] * fe_0 + 2.0 * tk_xzzz_xzzzz[i] * fe_0 +
                             tk_xzzz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xxzzzz[i] * fz_0;

        tk_xxzzz_xyyyyy[i] =
            -4.0 * ts_xxz_xyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xyyyyy[i] * fe_0 + tk_xxzz_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxzzz_xyyyyy[i] * fz_0;

        tk_xxzzz_xyyyyz[i] = -2.0 * ts_zzz_xyyyyz[i] * fbe_0 * fz_0 + tk_zzz_xyyyyz[i] * fe_0 + tk_xzzz_yyyyz[i] * fe_0 +
                             tk_xzzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xxzzz_xyyyyz[i] * fz_0;

        tk_xxzzz_xyyyzz[i] = -2.0 * ts_zzz_xyyyzz[i] * fbe_0 * fz_0 + tk_zzz_xyyyzz[i] * fe_0 + tk_xzzz_yyyzz[i] * fe_0 +
                             tk_xzzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xyyyzz[i] * fz_0;

        tk_xxzzz_xyyzzz[i] = -2.0 * ts_zzz_xyyzzz[i] * fbe_0 * fz_0 + tk_zzz_xyyzzz[i] * fe_0 + tk_xzzz_yyzzz[i] * fe_0 +
                             tk_xzzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xyyzzz[i] * fz_0;

        tk_xxzzz_xyzzzz[i] = -2.0 * ts_zzz_xyzzzz[i] * fbe_0 * fz_0 + tk_zzz_xyzzzz[i] * fe_0 + tk_xzzz_yzzzz[i] * fe_0 +
                             tk_xzzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xyzzzz[i] * fz_0;

        tk_xxzzz_xzzzzz[i] = -2.0 * ts_zzz_xzzzzz[i] * fbe_0 * fz_0 + tk_zzz_xzzzzz[i] * fe_0 + tk_xzzz_zzzzz[i] * fe_0 +
                             tk_xzzz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_xzzzzz[i] * fz_0;

        tk_xxzzz_yyyyyy[i] =
            -2.0 * ts_zzz_yyyyyy[i] * fbe_0 * fz_0 + tk_zzz_yyyyyy[i] * fe_0 + tk_xzzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xxzzz_yyyyyy[i] * fz_0;

        tk_xxzzz_yyyyyz[i] =
            -2.0 * ts_zzz_yyyyyz[i] * fbe_0 * fz_0 + tk_zzz_yyyyyz[i] * fe_0 + tk_xzzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xxzzz_yyyyyz[i] * fz_0;

        tk_xxzzz_yyyyzz[i] =
            -2.0 * ts_zzz_yyyyzz[i] * fbe_0 * fz_0 + tk_zzz_yyyyzz[i] * fe_0 + tk_xzzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xxzzz_yyyyzz[i] * fz_0;

        tk_xxzzz_yyyzzz[i] =
            -2.0 * ts_zzz_yyyzzz[i] * fbe_0 * fz_0 + tk_zzz_yyyzzz[i] * fe_0 + tk_xzzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_yyyzzz[i] * fz_0;

        tk_xxzzz_yyzzzz[i] =
            -2.0 * ts_zzz_yyzzzz[i] * fbe_0 * fz_0 + tk_zzz_yyzzzz[i] * fe_0 + tk_xzzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_yyzzzz[i] * fz_0;

        tk_xxzzz_yzzzzz[i] =
            -2.0 * ts_zzz_yzzzzz[i] * fbe_0 * fz_0 + tk_zzz_yzzzzz[i] * fe_0 + tk_xzzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_yzzzzz[i] * fz_0;

        tk_xxzzz_zzzzzz[i] =
            -2.0 * ts_zzz_zzzzzz[i] * fbe_0 * fz_0 + tk_zzz_zzzzzz[i] * fe_0 + tk_xzzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xxzzz_zzzzzz[i] * fz_0;
    }

    // Set up 280-308 components of targeted buffer : HI

    auto tk_xyyyy_xxxxxx = pbuffer.data(idx_kin_hi + 280);

    auto tk_xyyyy_xxxxxy = pbuffer.data(idx_kin_hi + 281);

    auto tk_xyyyy_xxxxxz = pbuffer.data(idx_kin_hi + 282);

    auto tk_xyyyy_xxxxyy = pbuffer.data(idx_kin_hi + 283);

    auto tk_xyyyy_xxxxyz = pbuffer.data(idx_kin_hi + 284);

    auto tk_xyyyy_xxxxzz = pbuffer.data(idx_kin_hi + 285);

    auto tk_xyyyy_xxxyyy = pbuffer.data(idx_kin_hi + 286);

    auto tk_xyyyy_xxxyyz = pbuffer.data(idx_kin_hi + 287);

    auto tk_xyyyy_xxxyzz = pbuffer.data(idx_kin_hi + 288);

    auto tk_xyyyy_xxxzzz = pbuffer.data(idx_kin_hi + 289);

    auto tk_xyyyy_xxyyyy = pbuffer.data(idx_kin_hi + 290);

    auto tk_xyyyy_xxyyyz = pbuffer.data(idx_kin_hi + 291);

    auto tk_xyyyy_xxyyzz = pbuffer.data(idx_kin_hi + 292);

    auto tk_xyyyy_xxyzzz = pbuffer.data(idx_kin_hi + 293);

    auto tk_xyyyy_xxzzzz = pbuffer.data(idx_kin_hi + 294);

    auto tk_xyyyy_xyyyyy = pbuffer.data(idx_kin_hi + 295);

    auto tk_xyyyy_xyyyyz = pbuffer.data(idx_kin_hi + 296);

    auto tk_xyyyy_xyyyzz = pbuffer.data(idx_kin_hi + 297);

    auto tk_xyyyy_xyyzzz = pbuffer.data(idx_kin_hi + 298);

    auto tk_xyyyy_xyzzzz = pbuffer.data(idx_kin_hi + 299);

    auto tk_xyyyy_xzzzzz = pbuffer.data(idx_kin_hi + 300);

    auto tk_xyyyy_yyyyyy = pbuffer.data(idx_kin_hi + 301);

    auto tk_xyyyy_yyyyyz = pbuffer.data(idx_kin_hi + 302);

    auto tk_xyyyy_yyyyzz = pbuffer.data(idx_kin_hi + 303);

    auto tk_xyyyy_yyyzzz = pbuffer.data(idx_kin_hi + 304);

    auto tk_xyyyy_yyzzzz = pbuffer.data(idx_kin_hi + 305);

    auto tk_xyyyy_yzzzzz = pbuffer.data(idx_kin_hi + 306);

    auto tk_xyyyy_zzzzzz = pbuffer.data(idx_kin_hi + 307);

#pragma omp simd aligned(pa_x,                \
                             tk_xyyyy_xxxxxx, \
                             tk_xyyyy_xxxxxy, \
                             tk_xyyyy_xxxxxz, \
                             tk_xyyyy_xxxxyy, \
                             tk_xyyyy_xxxxyz, \
                             tk_xyyyy_xxxxzz, \
                             tk_xyyyy_xxxyyy, \
                             tk_xyyyy_xxxyyz, \
                             tk_xyyyy_xxxyzz, \
                             tk_xyyyy_xxxzzz, \
                             tk_xyyyy_xxyyyy, \
                             tk_xyyyy_xxyyyz, \
                             tk_xyyyy_xxyyzz, \
                             tk_xyyyy_xxyzzz, \
                             tk_xyyyy_xxzzzz, \
                             tk_xyyyy_xyyyyy, \
                             tk_xyyyy_xyyyyz, \
                             tk_xyyyy_xyyyzz, \
                             tk_xyyyy_xyyzzz, \
                             tk_xyyyy_xyzzzz, \
                             tk_xyyyy_xzzzzz, \
                             tk_xyyyy_yyyyyy, \
                             tk_xyyyy_yyyyyz, \
                             tk_xyyyy_yyyyzz, \
                             tk_xyyyy_yyyzzz, \
                             tk_xyyyy_yyzzzz, \
                             tk_xyyyy_yzzzzz, \
                             tk_xyyyy_zzzzzz, \
                             tk_yyyy_xxxxx,   \
                             tk_yyyy_xxxxxx,  \
                             tk_yyyy_xxxxxy,  \
                             tk_yyyy_xxxxxz,  \
                             tk_yyyy_xxxxy,   \
                             tk_yyyy_xxxxyy,  \
                             tk_yyyy_xxxxyz,  \
                             tk_yyyy_xxxxz,   \
                             tk_yyyy_xxxxzz,  \
                             tk_yyyy_xxxyy,   \
                             tk_yyyy_xxxyyy,  \
                             tk_yyyy_xxxyyz,  \
                             tk_yyyy_xxxyz,   \
                             tk_yyyy_xxxyzz,  \
                             tk_yyyy_xxxzz,   \
                             tk_yyyy_xxxzzz,  \
                             tk_yyyy_xxyyy,   \
                             tk_yyyy_xxyyyy,  \
                             tk_yyyy_xxyyyz,  \
                             tk_yyyy_xxyyz,   \
                             tk_yyyy_xxyyzz,  \
                             tk_yyyy_xxyzz,   \
                             tk_yyyy_xxyzzz,  \
                             tk_yyyy_xxzzz,   \
                             tk_yyyy_xxzzzz,  \
                             tk_yyyy_xyyyy,   \
                             tk_yyyy_xyyyyy,  \
                             tk_yyyy_xyyyyz,  \
                             tk_yyyy_xyyyz,   \
                             tk_yyyy_xyyyzz,  \
                             tk_yyyy_xyyzz,   \
                             tk_yyyy_xyyzzz,  \
                             tk_yyyy_xyzzz,   \
                             tk_yyyy_xyzzzz,  \
                             tk_yyyy_xzzzz,   \
                             tk_yyyy_xzzzzz,  \
                             tk_yyyy_yyyyy,   \
                             tk_yyyy_yyyyyy,  \
                             tk_yyyy_yyyyyz,  \
                             tk_yyyy_yyyyz,   \
                             tk_yyyy_yyyyzz,  \
                             tk_yyyy_yyyzz,   \
                             tk_yyyy_yyyzzz,  \
                             tk_yyyy_yyzzz,   \
                             tk_yyyy_yyzzzz,  \
                             tk_yyyy_yzzzz,   \
                             tk_yyyy_yzzzzz,  \
                             tk_yyyy_zzzzz,   \
                             tk_yyyy_zzzzzz,  \
                             ts_xyyyy_xxxxxx, \
                             ts_xyyyy_xxxxxy, \
                             ts_xyyyy_xxxxxz, \
                             ts_xyyyy_xxxxyy, \
                             ts_xyyyy_xxxxyz, \
                             ts_xyyyy_xxxxzz, \
                             ts_xyyyy_xxxyyy, \
                             ts_xyyyy_xxxyyz, \
                             ts_xyyyy_xxxyzz, \
                             ts_xyyyy_xxxzzz, \
                             ts_xyyyy_xxyyyy, \
                             ts_xyyyy_xxyyyz, \
                             ts_xyyyy_xxyyzz, \
                             ts_xyyyy_xxyzzz, \
                             ts_xyyyy_xxzzzz, \
                             ts_xyyyy_xyyyyy, \
                             ts_xyyyy_xyyyyz, \
                             ts_xyyyy_xyyyzz, \
                             ts_xyyyy_xyyzzz, \
                             ts_xyyyy_xyzzzz, \
                             ts_xyyyy_xzzzzz, \
                             ts_xyyyy_yyyyyy, \
                             ts_xyyyy_yyyyyz, \
                             ts_xyyyy_yyyyzz, \
                             ts_xyyyy_yyyzzz, \
                             ts_xyyyy_yyzzzz, \
                             ts_xyyyy_yzzzzz, \
                             ts_xyyyy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyy_xxxxxx[i] = 6.0 * tk_yyyy_xxxxx[i] * fe_0 + tk_yyyy_xxxxxx[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxxxx[i] * fz_0;

        tk_xyyyy_xxxxxy[i] = 5.0 * tk_yyyy_xxxxy[i] * fe_0 + tk_yyyy_xxxxxy[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxxxy[i] * fz_0;

        tk_xyyyy_xxxxxz[i] = 5.0 * tk_yyyy_xxxxz[i] * fe_0 + tk_yyyy_xxxxxz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxxxz[i] * fz_0;

        tk_xyyyy_xxxxyy[i] = 4.0 * tk_yyyy_xxxyy[i] * fe_0 + tk_yyyy_xxxxyy[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxxyy[i] * fz_0;

        tk_xyyyy_xxxxyz[i] = 4.0 * tk_yyyy_xxxyz[i] * fe_0 + tk_yyyy_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxxyz[i] * fz_0;

        tk_xyyyy_xxxxzz[i] = 4.0 * tk_yyyy_xxxzz[i] * fe_0 + tk_yyyy_xxxxzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxxzz[i] * fz_0;

        tk_xyyyy_xxxyyy[i] = 3.0 * tk_yyyy_xxyyy[i] * fe_0 + tk_yyyy_xxxyyy[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxyyy[i] * fz_0;

        tk_xyyyy_xxxyyz[i] = 3.0 * tk_yyyy_xxyyz[i] * fe_0 + tk_yyyy_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxyyz[i] * fz_0;

        tk_xyyyy_xxxyzz[i] = 3.0 * tk_yyyy_xxyzz[i] * fe_0 + tk_yyyy_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxyzz[i] * fz_0;

        tk_xyyyy_xxxzzz[i] = 3.0 * tk_yyyy_xxzzz[i] * fe_0 + tk_yyyy_xxxzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxxzzz[i] * fz_0;

        tk_xyyyy_xxyyyy[i] = 2.0 * tk_yyyy_xyyyy[i] * fe_0 + tk_yyyy_xxyyyy[i] * pa_x[i] + 2.0 * ts_xyyyy_xxyyyy[i] * fz_0;

        tk_xyyyy_xxyyyz[i] = 2.0 * tk_yyyy_xyyyz[i] * fe_0 + tk_yyyy_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxyyyz[i] * fz_0;

        tk_xyyyy_xxyyzz[i] = 2.0 * tk_yyyy_xyyzz[i] * fe_0 + tk_yyyy_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxyyzz[i] * fz_0;

        tk_xyyyy_xxyzzz[i] = 2.0 * tk_yyyy_xyzzz[i] * fe_0 + tk_yyyy_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxyzzz[i] * fz_0;

        tk_xyyyy_xxzzzz[i] = 2.0 * tk_yyyy_xzzzz[i] * fe_0 + tk_yyyy_xxzzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xxzzzz[i] * fz_0;

        tk_xyyyy_xyyyyy[i] = tk_yyyy_yyyyy[i] * fe_0 + tk_yyyy_xyyyyy[i] * pa_x[i] + 2.0 * ts_xyyyy_xyyyyy[i] * fz_0;

        tk_xyyyy_xyyyyz[i] = tk_yyyy_yyyyz[i] * fe_0 + tk_yyyy_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyyyy_xyyyyz[i] * fz_0;

        tk_xyyyy_xyyyzz[i] = tk_yyyy_yyyzz[i] * fe_0 + tk_yyyy_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xyyyzz[i] * fz_0;

        tk_xyyyy_xyyzzz[i] = tk_yyyy_yyzzz[i] * fe_0 + tk_yyyy_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xyyzzz[i] * fz_0;

        tk_xyyyy_xyzzzz[i] = tk_yyyy_yzzzz[i] * fe_0 + tk_yyyy_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xyzzzz[i] * fz_0;

        tk_xyyyy_xzzzzz[i] = tk_yyyy_zzzzz[i] * fe_0 + tk_yyyy_xzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_xzzzzz[i] * fz_0;

        tk_xyyyy_yyyyyy[i] = tk_yyyy_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyyyy_yyyyyy[i] * fz_0;

        tk_xyyyy_yyyyyz[i] = tk_yyyy_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyyyy_yyyyyz[i] * fz_0;

        tk_xyyyy_yyyyzz[i] = tk_yyyy_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyyyy_yyyyzz[i] * fz_0;

        tk_xyyyy_yyyzzz[i] = tk_yyyy_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_yyyzzz[i] * fz_0;

        tk_xyyyy_yyzzzz[i] = tk_yyyy_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_yyzzzz[i] * fz_0;

        tk_xyyyy_yzzzzz[i] = tk_yyyy_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_yzzzzz[i] * fz_0;

        tk_xyyyy_zzzzzz[i] = tk_yyyy_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyy_zzzzzz[i] * fz_0;
    }

    // Set up 308-336 components of targeted buffer : HI

    auto tk_xyyyz_xxxxxx = pbuffer.data(idx_kin_hi + 308);

    auto tk_xyyyz_xxxxxy = pbuffer.data(idx_kin_hi + 309);

    auto tk_xyyyz_xxxxxz = pbuffer.data(idx_kin_hi + 310);

    auto tk_xyyyz_xxxxyy = pbuffer.data(idx_kin_hi + 311);

    auto tk_xyyyz_xxxxyz = pbuffer.data(idx_kin_hi + 312);

    auto tk_xyyyz_xxxxzz = pbuffer.data(idx_kin_hi + 313);

    auto tk_xyyyz_xxxyyy = pbuffer.data(idx_kin_hi + 314);

    auto tk_xyyyz_xxxyyz = pbuffer.data(idx_kin_hi + 315);

    auto tk_xyyyz_xxxyzz = pbuffer.data(idx_kin_hi + 316);

    auto tk_xyyyz_xxxzzz = pbuffer.data(idx_kin_hi + 317);

    auto tk_xyyyz_xxyyyy = pbuffer.data(idx_kin_hi + 318);

    auto tk_xyyyz_xxyyyz = pbuffer.data(idx_kin_hi + 319);

    auto tk_xyyyz_xxyyzz = pbuffer.data(idx_kin_hi + 320);

    auto tk_xyyyz_xxyzzz = pbuffer.data(idx_kin_hi + 321);

    auto tk_xyyyz_xxzzzz = pbuffer.data(idx_kin_hi + 322);

    auto tk_xyyyz_xyyyyy = pbuffer.data(idx_kin_hi + 323);

    auto tk_xyyyz_xyyyyz = pbuffer.data(idx_kin_hi + 324);

    auto tk_xyyyz_xyyyzz = pbuffer.data(idx_kin_hi + 325);

    auto tk_xyyyz_xyyzzz = pbuffer.data(idx_kin_hi + 326);

    auto tk_xyyyz_xyzzzz = pbuffer.data(idx_kin_hi + 327);

    auto tk_xyyyz_xzzzzz = pbuffer.data(idx_kin_hi + 328);

    auto tk_xyyyz_yyyyyy = pbuffer.data(idx_kin_hi + 329);

    auto tk_xyyyz_yyyyyz = pbuffer.data(idx_kin_hi + 330);

    auto tk_xyyyz_yyyyzz = pbuffer.data(idx_kin_hi + 331);

    auto tk_xyyyz_yyyzzz = pbuffer.data(idx_kin_hi + 332);

    auto tk_xyyyz_yyzzzz = pbuffer.data(idx_kin_hi + 333);

    auto tk_xyyyz_yzzzzz = pbuffer.data(idx_kin_hi + 334);

    auto tk_xyyyz_zzzzzz = pbuffer.data(idx_kin_hi + 335);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tk_xyyy_xxxxxx,  \
                             tk_xyyy_xxxxxy,  \
                             tk_xyyy_xxxxyy,  \
                             tk_xyyy_xxxyyy,  \
                             tk_xyyy_xxyyyy,  \
                             tk_xyyy_xyyyyy,  \
                             tk_xyyyz_xxxxxx, \
                             tk_xyyyz_xxxxxy, \
                             tk_xyyyz_xxxxxz, \
                             tk_xyyyz_xxxxyy, \
                             tk_xyyyz_xxxxyz, \
                             tk_xyyyz_xxxxzz, \
                             tk_xyyyz_xxxyyy, \
                             tk_xyyyz_xxxyyz, \
                             tk_xyyyz_xxxyzz, \
                             tk_xyyyz_xxxzzz, \
                             tk_xyyyz_xxyyyy, \
                             tk_xyyyz_xxyyyz, \
                             tk_xyyyz_xxyyzz, \
                             tk_xyyyz_xxyzzz, \
                             tk_xyyyz_xxzzzz, \
                             tk_xyyyz_xyyyyy, \
                             tk_xyyyz_xyyyyz, \
                             tk_xyyyz_xyyyzz, \
                             tk_xyyyz_xyyzzz, \
                             tk_xyyyz_xyzzzz, \
                             tk_xyyyz_xzzzzz, \
                             tk_xyyyz_yyyyyy, \
                             tk_xyyyz_yyyyyz, \
                             tk_xyyyz_yyyyzz, \
                             tk_xyyyz_yyyzzz, \
                             tk_xyyyz_yyzzzz, \
                             tk_xyyyz_yzzzzz, \
                             tk_xyyyz_zzzzzz, \
                             tk_yyyz_xxxxxz,  \
                             tk_yyyz_xxxxyz,  \
                             tk_yyyz_xxxxz,   \
                             tk_yyyz_xxxxzz,  \
                             tk_yyyz_xxxyyz,  \
                             tk_yyyz_xxxyz,   \
                             tk_yyyz_xxxyzz,  \
                             tk_yyyz_xxxzz,   \
                             tk_yyyz_xxxzzz,  \
                             tk_yyyz_xxyyyz,  \
                             tk_yyyz_xxyyz,   \
                             tk_yyyz_xxyyzz,  \
                             tk_yyyz_xxyzz,   \
                             tk_yyyz_xxyzzz,  \
                             tk_yyyz_xxzzz,   \
                             tk_yyyz_xxzzzz,  \
                             tk_yyyz_xyyyyz,  \
                             tk_yyyz_xyyyz,   \
                             tk_yyyz_xyyyzz,  \
                             tk_yyyz_xyyzz,   \
                             tk_yyyz_xyyzzz,  \
                             tk_yyyz_xyzzz,   \
                             tk_yyyz_xyzzzz,  \
                             tk_yyyz_xzzzz,   \
                             tk_yyyz_xzzzzz,  \
                             tk_yyyz_yyyyyy,  \
                             tk_yyyz_yyyyyz,  \
                             tk_yyyz_yyyyz,   \
                             tk_yyyz_yyyyzz,  \
                             tk_yyyz_yyyzz,   \
                             tk_yyyz_yyyzzz,  \
                             tk_yyyz_yyzzz,   \
                             tk_yyyz_yyzzzz,  \
                             tk_yyyz_yzzzz,   \
                             tk_yyyz_yzzzzz,  \
                             tk_yyyz_zzzzz,   \
                             tk_yyyz_zzzzzz,  \
                             ts_xyyyz_xxxxxx, \
                             ts_xyyyz_xxxxxy, \
                             ts_xyyyz_xxxxxz, \
                             ts_xyyyz_xxxxyy, \
                             ts_xyyyz_xxxxyz, \
                             ts_xyyyz_xxxxzz, \
                             ts_xyyyz_xxxyyy, \
                             ts_xyyyz_xxxyyz, \
                             ts_xyyyz_xxxyzz, \
                             ts_xyyyz_xxxzzz, \
                             ts_xyyyz_xxyyyy, \
                             ts_xyyyz_xxyyyz, \
                             ts_xyyyz_xxyyzz, \
                             ts_xyyyz_xxyzzz, \
                             ts_xyyyz_xxzzzz, \
                             ts_xyyyz_xyyyyy, \
                             ts_xyyyz_xyyyyz, \
                             ts_xyyyz_xyyyzz, \
                             ts_xyyyz_xyyzzz, \
                             ts_xyyyz_xyzzzz, \
                             ts_xyyyz_xzzzzz, \
                             ts_xyyyz_yyyyyy, \
                             ts_xyyyz_yyyyyz, \
                             ts_xyyyz_yyyyzz, \
                             ts_xyyyz_yyyzzz, \
                             ts_xyyyz_yyzzzz, \
                             ts_xyyyz_yzzzzz, \
                             ts_xyyyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyz_xxxxxx[i] = tk_xyyy_xxxxxx[i] * pa_z[i] + 2.0 * ts_xyyyz_xxxxxx[i] * fz_0;

        tk_xyyyz_xxxxxy[i] = tk_xyyy_xxxxxy[i] * pa_z[i] + 2.0 * ts_xyyyz_xxxxxy[i] * fz_0;

        tk_xyyyz_xxxxxz[i] = 5.0 * tk_yyyz_xxxxz[i] * fe_0 + tk_yyyz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxxxxz[i] * fz_0;

        tk_xyyyz_xxxxyy[i] = tk_xyyy_xxxxyy[i] * pa_z[i] + 2.0 * ts_xyyyz_xxxxyy[i] * fz_0;

        tk_xyyyz_xxxxyz[i] = 4.0 * tk_yyyz_xxxyz[i] * fe_0 + tk_yyyz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxxxyz[i] * fz_0;

        tk_xyyyz_xxxxzz[i] = 4.0 * tk_yyyz_xxxzz[i] * fe_0 + tk_yyyz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxxxzz[i] * fz_0;

        tk_xyyyz_xxxyyy[i] = tk_xyyy_xxxyyy[i] * pa_z[i] + 2.0 * ts_xyyyz_xxxyyy[i] * fz_0;

        tk_xyyyz_xxxyyz[i] = 3.0 * tk_yyyz_xxyyz[i] * fe_0 + tk_yyyz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxxyyz[i] * fz_0;

        tk_xyyyz_xxxyzz[i] = 3.0 * tk_yyyz_xxyzz[i] * fe_0 + tk_yyyz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxxyzz[i] * fz_0;

        tk_xyyyz_xxxzzz[i] = 3.0 * tk_yyyz_xxzzz[i] * fe_0 + tk_yyyz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxxzzz[i] * fz_0;

        tk_xyyyz_xxyyyy[i] = tk_xyyy_xxyyyy[i] * pa_z[i] + 2.0 * ts_xyyyz_xxyyyy[i] * fz_0;

        tk_xyyyz_xxyyyz[i] = 2.0 * tk_yyyz_xyyyz[i] * fe_0 + tk_yyyz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxyyyz[i] * fz_0;

        tk_xyyyz_xxyyzz[i] = 2.0 * tk_yyyz_xyyzz[i] * fe_0 + tk_yyyz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxyyzz[i] * fz_0;

        tk_xyyyz_xxyzzz[i] = 2.0 * tk_yyyz_xyzzz[i] * fe_0 + tk_yyyz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxyzzz[i] * fz_0;

        tk_xyyyz_xxzzzz[i] = 2.0 * tk_yyyz_xzzzz[i] * fe_0 + tk_yyyz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xxzzzz[i] * fz_0;

        tk_xyyyz_xyyyyy[i] = tk_xyyy_xyyyyy[i] * pa_z[i] + 2.0 * ts_xyyyz_xyyyyy[i] * fz_0;

        tk_xyyyz_xyyyyz[i] = tk_yyyz_yyyyz[i] * fe_0 + tk_yyyz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyyyz_xyyyyz[i] * fz_0;

        tk_xyyyz_xyyyzz[i] = tk_yyyz_yyyzz[i] * fe_0 + tk_yyyz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xyyyzz[i] * fz_0;

        tk_xyyyz_xyyzzz[i] = tk_yyyz_yyzzz[i] * fe_0 + tk_yyyz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xyyzzz[i] * fz_0;

        tk_xyyyz_xyzzzz[i] = tk_yyyz_yzzzz[i] * fe_0 + tk_yyyz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xyzzzz[i] * fz_0;

        tk_xyyyz_xzzzzz[i] = tk_yyyz_zzzzz[i] * fe_0 + tk_yyyz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_xzzzzz[i] * fz_0;

        tk_xyyyz_yyyyyy[i] = tk_yyyz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyyyz_yyyyyy[i] * fz_0;

        tk_xyyyz_yyyyyz[i] = tk_yyyz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyyyz_yyyyyz[i] * fz_0;

        tk_xyyyz_yyyyzz[i] = tk_yyyz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyyyz_yyyyzz[i] * fz_0;

        tk_xyyyz_yyyzzz[i] = tk_yyyz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_yyyzzz[i] * fz_0;

        tk_xyyyz_yyzzzz[i] = tk_yyyz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_yyzzzz[i] * fz_0;

        tk_xyyyz_yzzzzz[i] = tk_yyyz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_yzzzzz[i] * fz_0;

        tk_xyyyz_zzzzzz[i] = tk_yyyz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyyyz_zzzzzz[i] * fz_0;
    }

    // Set up 336-364 components of targeted buffer : HI

    auto tk_xyyzz_xxxxxx = pbuffer.data(idx_kin_hi + 336);

    auto tk_xyyzz_xxxxxy = pbuffer.data(idx_kin_hi + 337);

    auto tk_xyyzz_xxxxxz = pbuffer.data(idx_kin_hi + 338);

    auto tk_xyyzz_xxxxyy = pbuffer.data(idx_kin_hi + 339);

    auto tk_xyyzz_xxxxyz = pbuffer.data(idx_kin_hi + 340);

    auto tk_xyyzz_xxxxzz = pbuffer.data(idx_kin_hi + 341);

    auto tk_xyyzz_xxxyyy = pbuffer.data(idx_kin_hi + 342);

    auto tk_xyyzz_xxxyyz = pbuffer.data(idx_kin_hi + 343);

    auto tk_xyyzz_xxxyzz = pbuffer.data(idx_kin_hi + 344);

    auto tk_xyyzz_xxxzzz = pbuffer.data(idx_kin_hi + 345);

    auto tk_xyyzz_xxyyyy = pbuffer.data(idx_kin_hi + 346);

    auto tk_xyyzz_xxyyyz = pbuffer.data(idx_kin_hi + 347);

    auto tk_xyyzz_xxyyzz = pbuffer.data(idx_kin_hi + 348);

    auto tk_xyyzz_xxyzzz = pbuffer.data(idx_kin_hi + 349);

    auto tk_xyyzz_xxzzzz = pbuffer.data(idx_kin_hi + 350);

    auto tk_xyyzz_xyyyyy = pbuffer.data(idx_kin_hi + 351);

    auto tk_xyyzz_xyyyyz = pbuffer.data(idx_kin_hi + 352);

    auto tk_xyyzz_xyyyzz = pbuffer.data(idx_kin_hi + 353);

    auto tk_xyyzz_xyyzzz = pbuffer.data(idx_kin_hi + 354);

    auto tk_xyyzz_xyzzzz = pbuffer.data(idx_kin_hi + 355);

    auto tk_xyyzz_xzzzzz = pbuffer.data(idx_kin_hi + 356);

    auto tk_xyyzz_yyyyyy = pbuffer.data(idx_kin_hi + 357);

    auto tk_xyyzz_yyyyyz = pbuffer.data(idx_kin_hi + 358);

    auto tk_xyyzz_yyyyzz = pbuffer.data(idx_kin_hi + 359);

    auto tk_xyyzz_yyyzzz = pbuffer.data(idx_kin_hi + 360);

    auto tk_xyyzz_yyzzzz = pbuffer.data(idx_kin_hi + 361);

    auto tk_xyyzz_yzzzzz = pbuffer.data(idx_kin_hi + 362);

    auto tk_xyyzz_zzzzzz = pbuffer.data(idx_kin_hi + 363);

#pragma omp simd aligned(pa_x,                \
                             tk_xyyzz_xxxxxx, \
                             tk_xyyzz_xxxxxy, \
                             tk_xyyzz_xxxxxz, \
                             tk_xyyzz_xxxxyy, \
                             tk_xyyzz_xxxxyz, \
                             tk_xyyzz_xxxxzz, \
                             tk_xyyzz_xxxyyy, \
                             tk_xyyzz_xxxyyz, \
                             tk_xyyzz_xxxyzz, \
                             tk_xyyzz_xxxzzz, \
                             tk_xyyzz_xxyyyy, \
                             tk_xyyzz_xxyyyz, \
                             tk_xyyzz_xxyyzz, \
                             tk_xyyzz_xxyzzz, \
                             tk_xyyzz_xxzzzz, \
                             tk_xyyzz_xyyyyy, \
                             tk_xyyzz_xyyyyz, \
                             tk_xyyzz_xyyyzz, \
                             tk_xyyzz_xyyzzz, \
                             tk_xyyzz_xyzzzz, \
                             tk_xyyzz_xzzzzz, \
                             tk_xyyzz_yyyyyy, \
                             tk_xyyzz_yyyyyz, \
                             tk_xyyzz_yyyyzz, \
                             tk_xyyzz_yyyzzz, \
                             tk_xyyzz_yyzzzz, \
                             tk_xyyzz_yzzzzz, \
                             tk_xyyzz_zzzzzz, \
                             tk_yyzz_xxxxx,   \
                             tk_yyzz_xxxxxx,  \
                             tk_yyzz_xxxxxy,  \
                             tk_yyzz_xxxxxz,  \
                             tk_yyzz_xxxxy,   \
                             tk_yyzz_xxxxyy,  \
                             tk_yyzz_xxxxyz,  \
                             tk_yyzz_xxxxz,   \
                             tk_yyzz_xxxxzz,  \
                             tk_yyzz_xxxyy,   \
                             tk_yyzz_xxxyyy,  \
                             tk_yyzz_xxxyyz,  \
                             tk_yyzz_xxxyz,   \
                             tk_yyzz_xxxyzz,  \
                             tk_yyzz_xxxzz,   \
                             tk_yyzz_xxxzzz,  \
                             tk_yyzz_xxyyy,   \
                             tk_yyzz_xxyyyy,  \
                             tk_yyzz_xxyyyz,  \
                             tk_yyzz_xxyyz,   \
                             tk_yyzz_xxyyzz,  \
                             tk_yyzz_xxyzz,   \
                             tk_yyzz_xxyzzz,  \
                             tk_yyzz_xxzzz,   \
                             tk_yyzz_xxzzzz,  \
                             tk_yyzz_xyyyy,   \
                             tk_yyzz_xyyyyy,  \
                             tk_yyzz_xyyyyz,  \
                             tk_yyzz_xyyyz,   \
                             tk_yyzz_xyyyzz,  \
                             tk_yyzz_xyyzz,   \
                             tk_yyzz_xyyzzz,  \
                             tk_yyzz_xyzzz,   \
                             tk_yyzz_xyzzzz,  \
                             tk_yyzz_xzzzz,   \
                             tk_yyzz_xzzzzz,  \
                             tk_yyzz_yyyyy,   \
                             tk_yyzz_yyyyyy,  \
                             tk_yyzz_yyyyyz,  \
                             tk_yyzz_yyyyz,   \
                             tk_yyzz_yyyyzz,  \
                             tk_yyzz_yyyzz,   \
                             tk_yyzz_yyyzzz,  \
                             tk_yyzz_yyzzz,   \
                             tk_yyzz_yyzzzz,  \
                             tk_yyzz_yzzzz,   \
                             tk_yyzz_yzzzzz,  \
                             tk_yyzz_zzzzz,   \
                             tk_yyzz_zzzzzz,  \
                             ts_xyyzz_xxxxxx, \
                             ts_xyyzz_xxxxxy, \
                             ts_xyyzz_xxxxxz, \
                             ts_xyyzz_xxxxyy, \
                             ts_xyyzz_xxxxyz, \
                             ts_xyyzz_xxxxzz, \
                             ts_xyyzz_xxxyyy, \
                             ts_xyyzz_xxxyyz, \
                             ts_xyyzz_xxxyzz, \
                             ts_xyyzz_xxxzzz, \
                             ts_xyyzz_xxyyyy, \
                             ts_xyyzz_xxyyyz, \
                             ts_xyyzz_xxyyzz, \
                             ts_xyyzz_xxyzzz, \
                             ts_xyyzz_xxzzzz, \
                             ts_xyyzz_xyyyyy, \
                             ts_xyyzz_xyyyyz, \
                             ts_xyyzz_xyyyzz, \
                             ts_xyyzz_xyyzzz, \
                             ts_xyyzz_xyzzzz, \
                             ts_xyyzz_xzzzzz, \
                             ts_xyyzz_yyyyyy, \
                             ts_xyyzz_yyyyyz, \
                             ts_xyyzz_yyyyzz, \
                             ts_xyyzz_yyyzzz, \
                             ts_xyyzz_yyzzzz, \
                             ts_xyyzz_yzzzzz, \
                             ts_xyyzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyzz_xxxxxx[i] = 6.0 * tk_yyzz_xxxxx[i] * fe_0 + tk_yyzz_xxxxxx[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxxxx[i] * fz_0;

        tk_xyyzz_xxxxxy[i] = 5.0 * tk_yyzz_xxxxy[i] * fe_0 + tk_yyzz_xxxxxy[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxxxy[i] * fz_0;

        tk_xyyzz_xxxxxz[i] = 5.0 * tk_yyzz_xxxxz[i] * fe_0 + tk_yyzz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxxxz[i] * fz_0;

        tk_xyyzz_xxxxyy[i] = 4.0 * tk_yyzz_xxxyy[i] * fe_0 + tk_yyzz_xxxxyy[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxxyy[i] * fz_0;

        tk_xyyzz_xxxxyz[i] = 4.0 * tk_yyzz_xxxyz[i] * fe_0 + tk_yyzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxxyz[i] * fz_0;

        tk_xyyzz_xxxxzz[i] = 4.0 * tk_yyzz_xxxzz[i] * fe_0 + tk_yyzz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxxzz[i] * fz_0;

        tk_xyyzz_xxxyyy[i] = 3.0 * tk_yyzz_xxyyy[i] * fe_0 + tk_yyzz_xxxyyy[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxyyy[i] * fz_0;

        tk_xyyzz_xxxyyz[i] = 3.0 * tk_yyzz_xxyyz[i] * fe_0 + tk_yyzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxyyz[i] * fz_0;

        tk_xyyzz_xxxyzz[i] = 3.0 * tk_yyzz_xxyzz[i] * fe_0 + tk_yyzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxyzz[i] * fz_0;

        tk_xyyzz_xxxzzz[i] = 3.0 * tk_yyzz_xxzzz[i] * fe_0 + tk_yyzz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxxzzz[i] * fz_0;

        tk_xyyzz_xxyyyy[i] = 2.0 * tk_yyzz_xyyyy[i] * fe_0 + tk_yyzz_xxyyyy[i] * pa_x[i] + 2.0 * ts_xyyzz_xxyyyy[i] * fz_0;

        tk_xyyzz_xxyyyz[i] = 2.0 * tk_yyzz_xyyyz[i] * fe_0 + tk_yyzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxyyyz[i] * fz_0;

        tk_xyyzz_xxyyzz[i] = 2.0 * tk_yyzz_xyyzz[i] * fe_0 + tk_yyzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxyyzz[i] * fz_0;

        tk_xyyzz_xxyzzz[i] = 2.0 * tk_yyzz_xyzzz[i] * fe_0 + tk_yyzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxyzzz[i] * fz_0;

        tk_xyyzz_xxzzzz[i] = 2.0 * tk_yyzz_xzzzz[i] * fe_0 + tk_yyzz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xxzzzz[i] * fz_0;

        tk_xyyzz_xyyyyy[i] = tk_yyzz_yyyyy[i] * fe_0 + tk_yyzz_xyyyyy[i] * pa_x[i] + 2.0 * ts_xyyzz_xyyyyy[i] * fz_0;

        tk_xyyzz_xyyyyz[i] = tk_yyzz_yyyyz[i] * fe_0 + tk_yyzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyyzz_xyyyyz[i] * fz_0;

        tk_xyyzz_xyyyzz[i] = tk_yyzz_yyyzz[i] * fe_0 + tk_yyzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xyyyzz[i] * fz_0;

        tk_xyyzz_xyyzzz[i] = tk_yyzz_yyzzz[i] * fe_0 + tk_yyzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xyyzzz[i] * fz_0;

        tk_xyyzz_xyzzzz[i] = tk_yyzz_yzzzz[i] * fe_0 + tk_yyzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xyzzzz[i] * fz_0;

        tk_xyyzz_xzzzzz[i] = tk_yyzz_zzzzz[i] * fe_0 + tk_yyzz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_xzzzzz[i] * fz_0;

        tk_xyyzz_yyyyyy[i] = tk_yyzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyyzz_yyyyyy[i] * fz_0;

        tk_xyyzz_yyyyyz[i] = tk_yyzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyyzz_yyyyyz[i] * fz_0;

        tk_xyyzz_yyyyzz[i] = tk_yyzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyyzz_yyyyzz[i] * fz_0;

        tk_xyyzz_yyyzzz[i] = tk_yyzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_yyyzzz[i] * fz_0;

        tk_xyyzz_yyzzzz[i] = tk_yyzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_yyzzzz[i] * fz_0;

        tk_xyyzz_yzzzzz[i] = tk_yyzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_yzzzzz[i] * fz_0;

        tk_xyyzz_zzzzzz[i] = tk_yyzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyyzz_zzzzzz[i] * fz_0;
    }

    // Set up 364-392 components of targeted buffer : HI

    auto tk_xyzzz_xxxxxx = pbuffer.data(idx_kin_hi + 364);

    auto tk_xyzzz_xxxxxy = pbuffer.data(idx_kin_hi + 365);

    auto tk_xyzzz_xxxxxz = pbuffer.data(idx_kin_hi + 366);

    auto tk_xyzzz_xxxxyy = pbuffer.data(idx_kin_hi + 367);

    auto tk_xyzzz_xxxxyz = pbuffer.data(idx_kin_hi + 368);

    auto tk_xyzzz_xxxxzz = pbuffer.data(idx_kin_hi + 369);

    auto tk_xyzzz_xxxyyy = pbuffer.data(idx_kin_hi + 370);

    auto tk_xyzzz_xxxyyz = pbuffer.data(idx_kin_hi + 371);

    auto tk_xyzzz_xxxyzz = pbuffer.data(idx_kin_hi + 372);

    auto tk_xyzzz_xxxzzz = pbuffer.data(idx_kin_hi + 373);

    auto tk_xyzzz_xxyyyy = pbuffer.data(idx_kin_hi + 374);

    auto tk_xyzzz_xxyyyz = pbuffer.data(idx_kin_hi + 375);

    auto tk_xyzzz_xxyyzz = pbuffer.data(idx_kin_hi + 376);

    auto tk_xyzzz_xxyzzz = pbuffer.data(idx_kin_hi + 377);

    auto tk_xyzzz_xxzzzz = pbuffer.data(idx_kin_hi + 378);

    auto tk_xyzzz_xyyyyy = pbuffer.data(idx_kin_hi + 379);

    auto tk_xyzzz_xyyyyz = pbuffer.data(idx_kin_hi + 380);

    auto tk_xyzzz_xyyyzz = pbuffer.data(idx_kin_hi + 381);

    auto tk_xyzzz_xyyzzz = pbuffer.data(idx_kin_hi + 382);

    auto tk_xyzzz_xyzzzz = pbuffer.data(idx_kin_hi + 383);

    auto tk_xyzzz_xzzzzz = pbuffer.data(idx_kin_hi + 384);

    auto tk_xyzzz_yyyyyy = pbuffer.data(idx_kin_hi + 385);

    auto tk_xyzzz_yyyyyz = pbuffer.data(idx_kin_hi + 386);

    auto tk_xyzzz_yyyyzz = pbuffer.data(idx_kin_hi + 387);

    auto tk_xyzzz_yyyzzz = pbuffer.data(idx_kin_hi + 388);

    auto tk_xyzzz_yyzzzz = pbuffer.data(idx_kin_hi + 389);

    auto tk_xyzzz_yzzzzz = pbuffer.data(idx_kin_hi + 390);

    auto tk_xyzzz_zzzzzz = pbuffer.data(idx_kin_hi + 391);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tk_xyzzz_xxxxxx, \
                             tk_xyzzz_xxxxxy, \
                             tk_xyzzz_xxxxxz, \
                             tk_xyzzz_xxxxyy, \
                             tk_xyzzz_xxxxyz, \
                             tk_xyzzz_xxxxzz, \
                             tk_xyzzz_xxxyyy, \
                             tk_xyzzz_xxxyyz, \
                             tk_xyzzz_xxxyzz, \
                             tk_xyzzz_xxxzzz, \
                             tk_xyzzz_xxyyyy, \
                             tk_xyzzz_xxyyyz, \
                             tk_xyzzz_xxyyzz, \
                             tk_xyzzz_xxyzzz, \
                             tk_xyzzz_xxzzzz, \
                             tk_xyzzz_xyyyyy, \
                             tk_xyzzz_xyyyyz, \
                             tk_xyzzz_xyyyzz, \
                             tk_xyzzz_xyyzzz, \
                             tk_xyzzz_xyzzzz, \
                             tk_xyzzz_xzzzzz, \
                             tk_xyzzz_yyyyyy, \
                             tk_xyzzz_yyyyyz, \
                             tk_xyzzz_yyyyzz, \
                             tk_xyzzz_yyyzzz, \
                             tk_xyzzz_yyzzzz, \
                             tk_xyzzz_yzzzzz, \
                             tk_xyzzz_zzzzzz, \
                             tk_xzzz_xxxxxx,  \
                             tk_xzzz_xxxxxz,  \
                             tk_xzzz_xxxxzz,  \
                             tk_xzzz_xxxzzz,  \
                             tk_xzzz_xxzzzz,  \
                             tk_xzzz_xzzzzz,  \
                             tk_yzzz_xxxxxy,  \
                             tk_yzzz_xxxxy,   \
                             tk_yzzz_xxxxyy,  \
                             tk_yzzz_xxxxyz,  \
                             tk_yzzz_xxxyy,   \
                             tk_yzzz_xxxyyy,  \
                             tk_yzzz_xxxyyz,  \
                             tk_yzzz_xxxyz,   \
                             tk_yzzz_xxxyzz,  \
                             tk_yzzz_xxyyy,   \
                             tk_yzzz_xxyyyy,  \
                             tk_yzzz_xxyyyz,  \
                             tk_yzzz_xxyyz,   \
                             tk_yzzz_xxyyzz,  \
                             tk_yzzz_xxyzz,   \
                             tk_yzzz_xxyzzz,  \
                             tk_yzzz_xyyyy,   \
                             tk_yzzz_xyyyyy,  \
                             tk_yzzz_xyyyyz,  \
                             tk_yzzz_xyyyz,   \
                             tk_yzzz_xyyyzz,  \
                             tk_yzzz_xyyzz,   \
                             tk_yzzz_xyyzzz,  \
                             tk_yzzz_xyzzz,   \
                             tk_yzzz_xyzzzz,  \
                             tk_yzzz_yyyyy,   \
                             tk_yzzz_yyyyyy,  \
                             tk_yzzz_yyyyyz,  \
                             tk_yzzz_yyyyz,   \
                             tk_yzzz_yyyyzz,  \
                             tk_yzzz_yyyzz,   \
                             tk_yzzz_yyyzzz,  \
                             tk_yzzz_yyzzz,   \
                             tk_yzzz_yyzzzz,  \
                             tk_yzzz_yzzzz,   \
                             tk_yzzz_yzzzzz,  \
                             tk_yzzz_zzzzzz,  \
                             ts_xyzzz_xxxxxx, \
                             ts_xyzzz_xxxxxy, \
                             ts_xyzzz_xxxxxz, \
                             ts_xyzzz_xxxxyy, \
                             ts_xyzzz_xxxxyz, \
                             ts_xyzzz_xxxxzz, \
                             ts_xyzzz_xxxyyy, \
                             ts_xyzzz_xxxyyz, \
                             ts_xyzzz_xxxyzz, \
                             ts_xyzzz_xxxzzz, \
                             ts_xyzzz_xxyyyy, \
                             ts_xyzzz_xxyyyz, \
                             ts_xyzzz_xxyyzz, \
                             ts_xyzzz_xxyzzz, \
                             ts_xyzzz_xxzzzz, \
                             ts_xyzzz_xyyyyy, \
                             ts_xyzzz_xyyyyz, \
                             ts_xyzzz_xyyyzz, \
                             ts_xyzzz_xyyzzz, \
                             ts_xyzzz_xyzzzz, \
                             ts_xyzzz_xzzzzz, \
                             ts_xyzzz_yyyyyy, \
                             ts_xyzzz_yyyyyz, \
                             ts_xyzzz_yyyyzz, \
                             ts_xyzzz_yyyzzz, \
                             ts_xyzzz_yyzzzz, \
                             ts_xyzzz_yzzzzz, \
                             ts_xyzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzzz_xxxxxx[i] = tk_xzzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_xyzzz_xxxxxx[i] * fz_0;

        tk_xyzzz_xxxxxy[i] = 5.0 * tk_yzzz_xxxxy[i] * fe_0 + tk_yzzz_xxxxxy[i] * pa_x[i] + 2.0 * ts_xyzzz_xxxxxy[i] * fz_0;

        tk_xyzzz_xxxxxz[i] = tk_xzzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_xyzzz_xxxxxz[i] * fz_0;

        tk_xyzzz_xxxxyy[i] = 4.0 * tk_yzzz_xxxyy[i] * fe_0 + tk_yzzz_xxxxyy[i] * pa_x[i] + 2.0 * ts_xyzzz_xxxxyy[i] * fz_0;

        tk_xyzzz_xxxxyz[i] = 4.0 * tk_yzzz_xxxyz[i] * fe_0 + tk_yzzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyzzz_xxxxyz[i] * fz_0;

        tk_xyzzz_xxxxzz[i] = tk_xzzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_xyzzz_xxxxzz[i] * fz_0;

        tk_xyzzz_xxxyyy[i] = 3.0 * tk_yzzz_xxyyy[i] * fe_0 + tk_yzzz_xxxyyy[i] * pa_x[i] + 2.0 * ts_xyzzz_xxxyyy[i] * fz_0;

        tk_xyzzz_xxxyyz[i] = 3.0 * tk_yzzz_xxyyz[i] * fe_0 + tk_yzzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyzzz_xxxyyz[i] * fz_0;

        tk_xyzzz_xxxyzz[i] = 3.0 * tk_yzzz_xxyzz[i] * fe_0 + tk_yzzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyzzz_xxxyzz[i] * fz_0;

        tk_xyzzz_xxxzzz[i] = tk_xzzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_xyzzz_xxxzzz[i] * fz_0;

        tk_xyzzz_xxyyyy[i] = 2.0 * tk_yzzz_xyyyy[i] * fe_0 + tk_yzzz_xxyyyy[i] * pa_x[i] + 2.0 * ts_xyzzz_xxyyyy[i] * fz_0;

        tk_xyzzz_xxyyyz[i] = 2.0 * tk_yzzz_xyyyz[i] * fe_0 + tk_yzzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyzzz_xxyyyz[i] * fz_0;

        tk_xyzzz_xxyyzz[i] = 2.0 * tk_yzzz_xyyzz[i] * fe_0 + tk_yzzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyzzz_xxyyzz[i] * fz_0;

        tk_xyzzz_xxyzzz[i] = 2.0 * tk_yzzz_xyzzz[i] * fe_0 + tk_yzzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_xxyzzz[i] * fz_0;

        tk_xyzzz_xxzzzz[i] = tk_xzzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_xyzzz_xxzzzz[i] * fz_0;

        tk_xyzzz_xyyyyy[i] = tk_yzzz_yyyyy[i] * fe_0 + tk_yzzz_xyyyyy[i] * pa_x[i] + 2.0 * ts_xyzzz_xyyyyy[i] * fz_0;

        tk_xyzzz_xyyyyz[i] = tk_yzzz_yyyyz[i] * fe_0 + tk_yzzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyzzz_xyyyyz[i] * fz_0;

        tk_xyzzz_xyyyzz[i] = tk_yzzz_yyyzz[i] * fe_0 + tk_yzzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyzzz_xyyyzz[i] * fz_0;

        tk_xyzzz_xyyzzz[i] = tk_yzzz_yyzzz[i] * fe_0 + tk_yzzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_xyyzzz[i] * fz_0;

        tk_xyzzz_xyzzzz[i] = tk_yzzz_yzzzz[i] * fe_0 + tk_yzzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_xyzzzz[i] * fz_0;

        tk_xyzzz_xzzzzz[i] = tk_xzzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_xyzzz_xzzzzz[i] * fz_0;

        tk_xyzzz_yyyyyy[i] = tk_yzzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyzzz_yyyyyy[i] * fz_0;

        tk_xyzzz_yyyyyz[i] = tk_yzzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyzzz_yyyyyz[i] * fz_0;

        tk_xyzzz_yyyyzz[i] = tk_yzzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyzzz_yyyyzz[i] * fz_0;

        tk_xyzzz_yyyzzz[i] = tk_yzzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_yyyzzz[i] * fz_0;

        tk_xyzzz_yyzzzz[i] = tk_yzzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_yyzzzz[i] * fz_0;

        tk_xyzzz_yzzzzz[i] = tk_yzzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_yzzzzz[i] * fz_0;

        tk_xyzzz_zzzzzz[i] = tk_yzzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyzzz_zzzzzz[i] * fz_0;
    }

    // Set up 392-420 components of targeted buffer : HI

    auto tk_xzzzz_xxxxxx = pbuffer.data(idx_kin_hi + 392);

    auto tk_xzzzz_xxxxxy = pbuffer.data(idx_kin_hi + 393);

    auto tk_xzzzz_xxxxxz = pbuffer.data(idx_kin_hi + 394);

    auto tk_xzzzz_xxxxyy = pbuffer.data(idx_kin_hi + 395);

    auto tk_xzzzz_xxxxyz = pbuffer.data(idx_kin_hi + 396);

    auto tk_xzzzz_xxxxzz = pbuffer.data(idx_kin_hi + 397);

    auto tk_xzzzz_xxxyyy = pbuffer.data(idx_kin_hi + 398);

    auto tk_xzzzz_xxxyyz = pbuffer.data(idx_kin_hi + 399);

    auto tk_xzzzz_xxxyzz = pbuffer.data(idx_kin_hi + 400);

    auto tk_xzzzz_xxxzzz = pbuffer.data(idx_kin_hi + 401);

    auto tk_xzzzz_xxyyyy = pbuffer.data(idx_kin_hi + 402);

    auto tk_xzzzz_xxyyyz = pbuffer.data(idx_kin_hi + 403);

    auto tk_xzzzz_xxyyzz = pbuffer.data(idx_kin_hi + 404);

    auto tk_xzzzz_xxyzzz = pbuffer.data(idx_kin_hi + 405);

    auto tk_xzzzz_xxzzzz = pbuffer.data(idx_kin_hi + 406);

    auto tk_xzzzz_xyyyyy = pbuffer.data(idx_kin_hi + 407);

    auto tk_xzzzz_xyyyyz = pbuffer.data(idx_kin_hi + 408);

    auto tk_xzzzz_xyyyzz = pbuffer.data(idx_kin_hi + 409);

    auto tk_xzzzz_xyyzzz = pbuffer.data(idx_kin_hi + 410);

    auto tk_xzzzz_xyzzzz = pbuffer.data(idx_kin_hi + 411);

    auto tk_xzzzz_xzzzzz = pbuffer.data(idx_kin_hi + 412);

    auto tk_xzzzz_yyyyyy = pbuffer.data(idx_kin_hi + 413);

    auto tk_xzzzz_yyyyyz = pbuffer.data(idx_kin_hi + 414);

    auto tk_xzzzz_yyyyzz = pbuffer.data(idx_kin_hi + 415);

    auto tk_xzzzz_yyyzzz = pbuffer.data(idx_kin_hi + 416);

    auto tk_xzzzz_yyzzzz = pbuffer.data(idx_kin_hi + 417);

    auto tk_xzzzz_yzzzzz = pbuffer.data(idx_kin_hi + 418);

    auto tk_xzzzz_zzzzzz = pbuffer.data(idx_kin_hi + 419);

#pragma omp simd aligned(pa_x,                \
                             tk_xzzzz_xxxxxx, \
                             tk_xzzzz_xxxxxy, \
                             tk_xzzzz_xxxxxz, \
                             tk_xzzzz_xxxxyy, \
                             tk_xzzzz_xxxxyz, \
                             tk_xzzzz_xxxxzz, \
                             tk_xzzzz_xxxyyy, \
                             tk_xzzzz_xxxyyz, \
                             tk_xzzzz_xxxyzz, \
                             tk_xzzzz_xxxzzz, \
                             tk_xzzzz_xxyyyy, \
                             tk_xzzzz_xxyyyz, \
                             tk_xzzzz_xxyyzz, \
                             tk_xzzzz_xxyzzz, \
                             tk_xzzzz_xxzzzz, \
                             tk_xzzzz_xyyyyy, \
                             tk_xzzzz_xyyyyz, \
                             tk_xzzzz_xyyyzz, \
                             tk_xzzzz_xyyzzz, \
                             tk_xzzzz_xyzzzz, \
                             tk_xzzzz_xzzzzz, \
                             tk_xzzzz_yyyyyy, \
                             tk_xzzzz_yyyyyz, \
                             tk_xzzzz_yyyyzz, \
                             tk_xzzzz_yyyzzz, \
                             tk_xzzzz_yyzzzz, \
                             tk_xzzzz_yzzzzz, \
                             tk_xzzzz_zzzzzz, \
                             tk_zzzz_xxxxx,   \
                             tk_zzzz_xxxxxx,  \
                             tk_zzzz_xxxxxy,  \
                             tk_zzzz_xxxxxz,  \
                             tk_zzzz_xxxxy,   \
                             tk_zzzz_xxxxyy,  \
                             tk_zzzz_xxxxyz,  \
                             tk_zzzz_xxxxz,   \
                             tk_zzzz_xxxxzz,  \
                             tk_zzzz_xxxyy,   \
                             tk_zzzz_xxxyyy,  \
                             tk_zzzz_xxxyyz,  \
                             tk_zzzz_xxxyz,   \
                             tk_zzzz_xxxyzz,  \
                             tk_zzzz_xxxzz,   \
                             tk_zzzz_xxxzzz,  \
                             tk_zzzz_xxyyy,   \
                             tk_zzzz_xxyyyy,  \
                             tk_zzzz_xxyyyz,  \
                             tk_zzzz_xxyyz,   \
                             tk_zzzz_xxyyzz,  \
                             tk_zzzz_xxyzz,   \
                             tk_zzzz_xxyzzz,  \
                             tk_zzzz_xxzzz,   \
                             tk_zzzz_xxzzzz,  \
                             tk_zzzz_xyyyy,   \
                             tk_zzzz_xyyyyy,  \
                             tk_zzzz_xyyyyz,  \
                             tk_zzzz_xyyyz,   \
                             tk_zzzz_xyyyzz,  \
                             tk_zzzz_xyyzz,   \
                             tk_zzzz_xyyzzz,  \
                             tk_zzzz_xyzzz,   \
                             tk_zzzz_xyzzzz,  \
                             tk_zzzz_xzzzz,   \
                             tk_zzzz_xzzzzz,  \
                             tk_zzzz_yyyyy,   \
                             tk_zzzz_yyyyyy,  \
                             tk_zzzz_yyyyyz,  \
                             tk_zzzz_yyyyz,   \
                             tk_zzzz_yyyyzz,  \
                             tk_zzzz_yyyzz,   \
                             tk_zzzz_yyyzzz,  \
                             tk_zzzz_yyzzz,   \
                             tk_zzzz_yyzzzz,  \
                             tk_zzzz_yzzzz,   \
                             tk_zzzz_yzzzzz,  \
                             tk_zzzz_zzzzz,   \
                             tk_zzzz_zzzzzz,  \
                             ts_xzzzz_xxxxxx, \
                             ts_xzzzz_xxxxxy, \
                             ts_xzzzz_xxxxxz, \
                             ts_xzzzz_xxxxyy, \
                             ts_xzzzz_xxxxyz, \
                             ts_xzzzz_xxxxzz, \
                             ts_xzzzz_xxxyyy, \
                             ts_xzzzz_xxxyyz, \
                             ts_xzzzz_xxxyzz, \
                             ts_xzzzz_xxxzzz, \
                             ts_xzzzz_xxyyyy, \
                             ts_xzzzz_xxyyyz, \
                             ts_xzzzz_xxyyzz, \
                             ts_xzzzz_xxyzzz, \
                             ts_xzzzz_xxzzzz, \
                             ts_xzzzz_xyyyyy, \
                             ts_xzzzz_xyyyyz, \
                             ts_xzzzz_xyyyzz, \
                             ts_xzzzz_xyyzzz, \
                             ts_xzzzz_xyzzzz, \
                             ts_xzzzz_xzzzzz, \
                             ts_xzzzz_yyyyyy, \
                             ts_xzzzz_yyyyyz, \
                             ts_xzzzz_yyyyzz, \
                             ts_xzzzz_yyyzzz, \
                             ts_xzzzz_yyzzzz, \
                             ts_xzzzz_yzzzzz, \
                             ts_xzzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzzz_xxxxxx[i] = 6.0 * tk_zzzz_xxxxx[i] * fe_0 + tk_zzzz_xxxxxx[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxxxx[i] * fz_0;

        tk_xzzzz_xxxxxy[i] = 5.0 * tk_zzzz_xxxxy[i] * fe_0 + tk_zzzz_xxxxxy[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxxxy[i] * fz_0;

        tk_xzzzz_xxxxxz[i] = 5.0 * tk_zzzz_xxxxz[i] * fe_0 + tk_zzzz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxxxz[i] * fz_0;

        tk_xzzzz_xxxxyy[i] = 4.0 * tk_zzzz_xxxyy[i] * fe_0 + tk_zzzz_xxxxyy[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxxyy[i] * fz_0;

        tk_xzzzz_xxxxyz[i] = 4.0 * tk_zzzz_xxxyz[i] * fe_0 + tk_zzzz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxxyz[i] * fz_0;

        tk_xzzzz_xxxxzz[i] = 4.0 * tk_zzzz_xxxzz[i] * fe_0 + tk_zzzz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxxzz[i] * fz_0;

        tk_xzzzz_xxxyyy[i] = 3.0 * tk_zzzz_xxyyy[i] * fe_0 + tk_zzzz_xxxyyy[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxyyy[i] * fz_0;

        tk_xzzzz_xxxyyz[i] = 3.0 * tk_zzzz_xxyyz[i] * fe_0 + tk_zzzz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxyyz[i] * fz_0;

        tk_xzzzz_xxxyzz[i] = 3.0 * tk_zzzz_xxyzz[i] * fe_0 + tk_zzzz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxyzz[i] * fz_0;

        tk_xzzzz_xxxzzz[i] = 3.0 * tk_zzzz_xxzzz[i] * fe_0 + tk_zzzz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxxzzz[i] * fz_0;

        tk_xzzzz_xxyyyy[i] = 2.0 * tk_zzzz_xyyyy[i] * fe_0 + tk_zzzz_xxyyyy[i] * pa_x[i] + 2.0 * ts_xzzzz_xxyyyy[i] * fz_0;

        tk_xzzzz_xxyyyz[i] = 2.0 * tk_zzzz_xyyyz[i] * fe_0 + tk_zzzz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxyyyz[i] * fz_0;

        tk_xzzzz_xxyyzz[i] = 2.0 * tk_zzzz_xyyzz[i] * fe_0 + tk_zzzz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxyyzz[i] * fz_0;

        tk_xzzzz_xxyzzz[i] = 2.0 * tk_zzzz_xyzzz[i] * fe_0 + tk_zzzz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxyzzz[i] * fz_0;

        tk_xzzzz_xxzzzz[i] = 2.0 * tk_zzzz_xzzzz[i] * fe_0 + tk_zzzz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xxzzzz[i] * fz_0;

        tk_xzzzz_xyyyyy[i] = tk_zzzz_yyyyy[i] * fe_0 + tk_zzzz_xyyyyy[i] * pa_x[i] + 2.0 * ts_xzzzz_xyyyyy[i] * fz_0;

        tk_xzzzz_xyyyyz[i] = tk_zzzz_yyyyz[i] * fe_0 + tk_zzzz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xzzzz_xyyyyz[i] * fz_0;

        tk_xzzzz_xyyyzz[i] = tk_zzzz_yyyzz[i] * fe_0 + tk_zzzz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xyyyzz[i] * fz_0;

        tk_xzzzz_xyyzzz[i] = tk_zzzz_yyzzz[i] * fe_0 + tk_zzzz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xyyzzz[i] * fz_0;

        tk_xzzzz_xyzzzz[i] = tk_zzzz_yzzzz[i] * fe_0 + tk_zzzz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xyzzzz[i] * fz_0;

        tk_xzzzz_xzzzzz[i] = tk_zzzz_zzzzz[i] * fe_0 + tk_zzzz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_xzzzzz[i] * fz_0;

        tk_xzzzz_yyyyyy[i] = tk_zzzz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xzzzz_yyyyyy[i] * fz_0;

        tk_xzzzz_yyyyyz[i] = tk_zzzz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xzzzz_yyyyyz[i] * fz_0;

        tk_xzzzz_yyyyzz[i] = tk_zzzz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xzzzz_yyyyzz[i] * fz_0;

        tk_xzzzz_yyyzzz[i] = tk_zzzz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_yyyzzz[i] * fz_0;

        tk_xzzzz_yyzzzz[i] = tk_zzzz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_yyzzzz[i] * fz_0;

        tk_xzzzz_yzzzzz[i] = tk_zzzz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_yzzzzz[i] * fz_0;

        tk_xzzzz_zzzzzz[i] = tk_zzzz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xzzzz_zzzzzz[i] * fz_0;
    }

    // Set up 420-448 components of targeted buffer : HI

    auto tk_yyyyy_xxxxxx = pbuffer.data(idx_kin_hi + 420);

    auto tk_yyyyy_xxxxxy = pbuffer.data(idx_kin_hi + 421);

    auto tk_yyyyy_xxxxxz = pbuffer.data(idx_kin_hi + 422);

    auto tk_yyyyy_xxxxyy = pbuffer.data(idx_kin_hi + 423);

    auto tk_yyyyy_xxxxyz = pbuffer.data(idx_kin_hi + 424);

    auto tk_yyyyy_xxxxzz = pbuffer.data(idx_kin_hi + 425);

    auto tk_yyyyy_xxxyyy = pbuffer.data(idx_kin_hi + 426);

    auto tk_yyyyy_xxxyyz = pbuffer.data(idx_kin_hi + 427);

    auto tk_yyyyy_xxxyzz = pbuffer.data(idx_kin_hi + 428);

    auto tk_yyyyy_xxxzzz = pbuffer.data(idx_kin_hi + 429);

    auto tk_yyyyy_xxyyyy = pbuffer.data(idx_kin_hi + 430);

    auto tk_yyyyy_xxyyyz = pbuffer.data(idx_kin_hi + 431);

    auto tk_yyyyy_xxyyzz = pbuffer.data(idx_kin_hi + 432);

    auto tk_yyyyy_xxyzzz = pbuffer.data(idx_kin_hi + 433);

    auto tk_yyyyy_xxzzzz = pbuffer.data(idx_kin_hi + 434);

    auto tk_yyyyy_xyyyyy = pbuffer.data(idx_kin_hi + 435);

    auto tk_yyyyy_xyyyyz = pbuffer.data(idx_kin_hi + 436);

    auto tk_yyyyy_xyyyzz = pbuffer.data(idx_kin_hi + 437);

    auto tk_yyyyy_xyyzzz = pbuffer.data(idx_kin_hi + 438);

    auto tk_yyyyy_xyzzzz = pbuffer.data(idx_kin_hi + 439);

    auto tk_yyyyy_xzzzzz = pbuffer.data(idx_kin_hi + 440);

    auto tk_yyyyy_yyyyyy = pbuffer.data(idx_kin_hi + 441);

    auto tk_yyyyy_yyyyyz = pbuffer.data(idx_kin_hi + 442);

    auto tk_yyyyy_yyyyzz = pbuffer.data(idx_kin_hi + 443);

    auto tk_yyyyy_yyyzzz = pbuffer.data(idx_kin_hi + 444);

    auto tk_yyyyy_yyzzzz = pbuffer.data(idx_kin_hi + 445);

    auto tk_yyyyy_yzzzzz = pbuffer.data(idx_kin_hi + 446);

    auto tk_yyyyy_zzzzzz = pbuffer.data(idx_kin_hi + 447);

#pragma omp simd aligned(pa_y,                \
                             tk_yyy_xxxxxx,   \
                             tk_yyy_xxxxxy,   \
                             tk_yyy_xxxxxz,   \
                             tk_yyy_xxxxyy,   \
                             tk_yyy_xxxxyz,   \
                             tk_yyy_xxxxzz,   \
                             tk_yyy_xxxyyy,   \
                             tk_yyy_xxxyyz,   \
                             tk_yyy_xxxyzz,   \
                             tk_yyy_xxxzzz,   \
                             tk_yyy_xxyyyy,   \
                             tk_yyy_xxyyyz,   \
                             tk_yyy_xxyyzz,   \
                             tk_yyy_xxyzzz,   \
                             tk_yyy_xxzzzz,   \
                             tk_yyy_xyyyyy,   \
                             tk_yyy_xyyyyz,   \
                             tk_yyy_xyyyzz,   \
                             tk_yyy_xyyzzz,   \
                             tk_yyy_xyzzzz,   \
                             tk_yyy_xzzzzz,   \
                             tk_yyy_yyyyyy,   \
                             tk_yyy_yyyyyz,   \
                             tk_yyy_yyyyzz,   \
                             tk_yyy_yyyzzz,   \
                             tk_yyy_yyzzzz,   \
                             tk_yyy_yzzzzz,   \
                             tk_yyy_zzzzzz,   \
                             tk_yyyy_xxxxx,   \
                             tk_yyyy_xxxxxx,  \
                             tk_yyyy_xxxxxy,  \
                             tk_yyyy_xxxxxz,  \
                             tk_yyyy_xxxxy,   \
                             tk_yyyy_xxxxyy,  \
                             tk_yyyy_xxxxyz,  \
                             tk_yyyy_xxxxz,   \
                             tk_yyyy_xxxxzz,  \
                             tk_yyyy_xxxyy,   \
                             tk_yyyy_xxxyyy,  \
                             tk_yyyy_xxxyyz,  \
                             tk_yyyy_xxxyz,   \
                             tk_yyyy_xxxyzz,  \
                             tk_yyyy_xxxzz,   \
                             tk_yyyy_xxxzzz,  \
                             tk_yyyy_xxyyy,   \
                             tk_yyyy_xxyyyy,  \
                             tk_yyyy_xxyyyz,  \
                             tk_yyyy_xxyyz,   \
                             tk_yyyy_xxyyzz,  \
                             tk_yyyy_xxyzz,   \
                             tk_yyyy_xxyzzz,  \
                             tk_yyyy_xxzzz,   \
                             tk_yyyy_xxzzzz,  \
                             tk_yyyy_xyyyy,   \
                             tk_yyyy_xyyyyy,  \
                             tk_yyyy_xyyyyz,  \
                             tk_yyyy_xyyyz,   \
                             tk_yyyy_xyyyzz,  \
                             tk_yyyy_xyyzz,   \
                             tk_yyyy_xyyzzz,  \
                             tk_yyyy_xyzzz,   \
                             tk_yyyy_xyzzzz,  \
                             tk_yyyy_xzzzz,   \
                             tk_yyyy_xzzzzz,  \
                             tk_yyyy_yyyyy,   \
                             tk_yyyy_yyyyyy,  \
                             tk_yyyy_yyyyyz,  \
                             tk_yyyy_yyyyz,   \
                             tk_yyyy_yyyyzz,  \
                             tk_yyyy_yyyzz,   \
                             tk_yyyy_yyyzzz,  \
                             tk_yyyy_yyzzz,   \
                             tk_yyyy_yyzzzz,  \
                             tk_yyyy_yzzzz,   \
                             tk_yyyy_yzzzzz,  \
                             tk_yyyy_zzzzz,   \
                             tk_yyyy_zzzzzz,  \
                             tk_yyyyy_xxxxxx, \
                             tk_yyyyy_xxxxxy, \
                             tk_yyyyy_xxxxxz, \
                             tk_yyyyy_xxxxyy, \
                             tk_yyyyy_xxxxyz, \
                             tk_yyyyy_xxxxzz, \
                             tk_yyyyy_xxxyyy, \
                             tk_yyyyy_xxxyyz, \
                             tk_yyyyy_xxxyzz, \
                             tk_yyyyy_xxxzzz, \
                             tk_yyyyy_xxyyyy, \
                             tk_yyyyy_xxyyyz, \
                             tk_yyyyy_xxyyzz, \
                             tk_yyyyy_xxyzzz, \
                             tk_yyyyy_xxzzzz, \
                             tk_yyyyy_xyyyyy, \
                             tk_yyyyy_xyyyyz, \
                             tk_yyyyy_xyyyzz, \
                             tk_yyyyy_xyyzzz, \
                             tk_yyyyy_xyzzzz, \
                             tk_yyyyy_xzzzzz, \
                             tk_yyyyy_yyyyyy, \
                             tk_yyyyy_yyyyyz, \
                             tk_yyyyy_yyyyzz, \
                             tk_yyyyy_yyyzzz, \
                             tk_yyyyy_yyzzzz, \
                             tk_yyyyy_yzzzzz, \
                             tk_yyyyy_zzzzzz, \
                             ts_yyy_xxxxxx,   \
                             ts_yyy_xxxxxy,   \
                             ts_yyy_xxxxxz,   \
                             ts_yyy_xxxxyy,   \
                             ts_yyy_xxxxyz,   \
                             ts_yyy_xxxxzz,   \
                             ts_yyy_xxxyyy,   \
                             ts_yyy_xxxyyz,   \
                             ts_yyy_xxxyzz,   \
                             ts_yyy_xxxzzz,   \
                             ts_yyy_xxyyyy,   \
                             ts_yyy_xxyyyz,   \
                             ts_yyy_xxyyzz,   \
                             ts_yyy_xxyzzz,   \
                             ts_yyy_xxzzzz,   \
                             ts_yyy_xyyyyy,   \
                             ts_yyy_xyyyyz,   \
                             ts_yyy_xyyyzz,   \
                             ts_yyy_xyyzzz,   \
                             ts_yyy_xyzzzz,   \
                             ts_yyy_xzzzzz,   \
                             ts_yyy_yyyyyy,   \
                             ts_yyy_yyyyyz,   \
                             ts_yyy_yyyyzz,   \
                             ts_yyy_yyyzzz,   \
                             ts_yyy_yyzzzz,   \
                             ts_yyy_yzzzzz,   \
                             ts_yyy_zzzzzz,   \
                             ts_yyyyy_xxxxxx, \
                             ts_yyyyy_xxxxxy, \
                             ts_yyyyy_xxxxxz, \
                             ts_yyyyy_xxxxyy, \
                             ts_yyyyy_xxxxyz, \
                             ts_yyyyy_xxxxzz, \
                             ts_yyyyy_xxxyyy, \
                             ts_yyyyy_xxxyyz, \
                             ts_yyyyy_xxxyzz, \
                             ts_yyyyy_xxxzzz, \
                             ts_yyyyy_xxyyyy, \
                             ts_yyyyy_xxyyyz, \
                             ts_yyyyy_xxyyzz, \
                             ts_yyyyy_xxyzzz, \
                             ts_yyyyy_xxzzzz, \
                             ts_yyyyy_xyyyyy, \
                             ts_yyyyy_xyyyyz, \
                             ts_yyyyy_xyyyzz, \
                             ts_yyyyy_xyyzzz, \
                             ts_yyyyy_xyzzzz, \
                             ts_yyyyy_xzzzzz, \
                             ts_yyyyy_yyyyyy, \
                             ts_yyyyy_yyyyyz, \
                             ts_yyyyy_yyyyzz, \
                             ts_yyyyy_yyyzzz, \
                             ts_yyyyy_yyzzzz, \
                             ts_yyyyy_yzzzzz, \
                             ts_yyyyy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyy_xxxxxx[i] =
            -8.0 * ts_yyy_xxxxxx[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxxxx[i] * fe_0 + tk_yyyy_xxxxxx[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxxxx[i] * fz_0;

        tk_yyyyy_xxxxxy[i] = -8.0 * ts_yyy_xxxxxy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxxxy[i] * fe_0 + tk_yyyy_xxxxx[i] * fe_0 +
                             tk_yyyy_xxxxxy[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxxxy[i] * fz_0;

        tk_yyyyy_xxxxxz[i] =
            -8.0 * ts_yyy_xxxxxz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxxxz[i] * fe_0 + tk_yyyy_xxxxxz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxxxz[i] * fz_0;

        tk_yyyyy_xxxxyy[i] = -8.0 * ts_yyy_xxxxyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxxyy[i] * fe_0 + 2.0 * tk_yyyy_xxxxy[i] * fe_0 +
                             tk_yyyy_xxxxyy[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxxyy[i] * fz_0;

        tk_yyyyy_xxxxyz[i] = -8.0 * ts_yyy_xxxxyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxxyz[i] * fe_0 + tk_yyyy_xxxxz[i] * fe_0 +
                             tk_yyyy_xxxxyz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxxyz[i] * fz_0;

        tk_yyyyy_xxxxzz[i] =
            -8.0 * ts_yyy_xxxxzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxxzz[i] * fe_0 + tk_yyyy_xxxxzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxxzz[i] * fz_0;

        tk_yyyyy_xxxyyy[i] = -8.0 * ts_yyy_xxxyyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxyyy[i] * fe_0 + 3.0 * tk_yyyy_xxxyy[i] * fe_0 +
                             tk_yyyy_xxxyyy[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxyyy[i] * fz_0;

        tk_yyyyy_xxxyyz[i] = -8.0 * ts_yyy_xxxyyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxyyz[i] * fe_0 + 2.0 * tk_yyyy_xxxyz[i] * fe_0 +
                             tk_yyyy_xxxyyz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxyyz[i] * fz_0;

        tk_yyyyy_xxxyzz[i] = -8.0 * ts_yyy_xxxyzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxyzz[i] * fe_0 + tk_yyyy_xxxzz[i] * fe_0 +
                             tk_yyyy_xxxyzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxyzz[i] * fz_0;

        tk_yyyyy_xxxzzz[i] =
            -8.0 * ts_yyy_xxxzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxxzzz[i] * fe_0 + tk_yyyy_xxxzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxxzzz[i] * fz_0;

        tk_yyyyy_xxyyyy[i] = -8.0 * ts_yyy_xxyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxyyyy[i] * fe_0 + 4.0 * tk_yyyy_xxyyy[i] * fe_0 +
                             tk_yyyy_xxyyyy[i] * pa_y[i] + 2.0 * ts_yyyyy_xxyyyy[i] * fz_0;

        tk_yyyyy_xxyyyz[i] = -8.0 * ts_yyy_xxyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxyyyz[i] * fe_0 + 3.0 * tk_yyyy_xxyyz[i] * fe_0 +
                             tk_yyyy_xxyyyz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxyyyz[i] * fz_0;

        tk_yyyyy_xxyyzz[i] = -8.0 * ts_yyy_xxyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxyyzz[i] * fe_0 + 2.0 * tk_yyyy_xxyzz[i] * fe_0 +
                             tk_yyyy_xxyyzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxyyzz[i] * fz_0;

        tk_yyyyy_xxyzzz[i] = -8.0 * ts_yyy_xxyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxyzzz[i] * fe_0 + tk_yyyy_xxzzz[i] * fe_0 +
                             tk_yyyy_xxyzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxyzzz[i] * fz_0;

        tk_yyyyy_xxzzzz[i] =
            -8.0 * ts_yyy_xxzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xxzzzz[i] * fe_0 + tk_yyyy_xxzzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xxzzzz[i] * fz_0;

        tk_yyyyy_xyyyyy[i] = -8.0 * ts_yyy_xyyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyyyyy[i] * fe_0 + 5.0 * tk_yyyy_xyyyy[i] * fe_0 +
                             tk_yyyy_xyyyyy[i] * pa_y[i] + 2.0 * ts_yyyyy_xyyyyy[i] * fz_0;

        tk_yyyyy_xyyyyz[i] = -8.0 * ts_yyy_xyyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyyyyz[i] * fe_0 + 4.0 * tk_yyyy_xyyyz[i] * fe_0 +
                             tk_yyyy_xyyyyz[i] * pa_y[i] + 2.0 * ts_yyyyy_xyyyyz[i] * fz_0;

        tk_yyyyy_xyyyzz[i] = -8.0 * ts_yyy_xyyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyyyzz[i] * fe_0 + 3.0 * tk_yyyy_xyyzz[i] * fe_0 +
                             tk_yyyy_xyyyzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xyyyzz[i] * fz_0;

        tk_yyyyy_xyyzzz[i] = -8.0 * ts_yyy_xyyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyyzzz[i] * fe_0 + 2.0 * tk_yyyy_xyzzz[i] * fe_0 +
                             tk_yyyy_xyyzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xyyzzz[i] * fz_0;

        tk_yyyyy_xyzzzz[i] = -8.0 * ts_yyy_xyzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xyzzzz[i] * fe_0 + tk_yyyy_xzzzz[i] * fe_0 +
                             tk_yyyy_xyzzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xyzzzz[i] * fz_0;

        tk_yyyyy_xzzzzz[i] =
            -8.0 * ts_yyy_xzzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xzzzzz[i] * fe_0 + tk_yyyy_xzzzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_xzzzzz[i] * fz_0;

        tk_yyyyy_yyyyyy[i] = -8.0 * ts_yyy_yyyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyyyyy[i] * fe_0 + 6.0 * tk_yyyy_yyyyy[i] * fe_0 +
                             tk_yyyy_yyyyyy[i] * pa_y[i] + 2.0 * ts_yyyyy_yyyyyy[i] * fz_0;

        tk_yyyyy_yyyyyz[i] = -8.0 * ts_yyy_yyyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyyyyz[i] * fe_0 + 5.0 * tk_yyyy_yyyyz[i] * fe_0 +
                             tk_yyyy_yyyyyz[i] * pa_y[i] + 2.0 * ts_yyyyy_yyyyyz[i] * fz_0;

        tk_yyyyy_yyyyzz[i] = -8.0 * ts_yyy_yyyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyyyzz[i] * fe_0 + 4.0 * tk_yyyy_yyyzz[i] * fe_0 +
                             tk_yyyy_yyyyzz[i] * pa_y[i] + 2.0 * ts_yyyyy_yyyyzz[i] * fz_0;

        tk_yyyyy_yyyzzz[i] = -8.0 * ts_yyy_yyyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyyzzz[i] * fe_0 + 3.0 * tk_yyyy_yyzzz[i] * fe_0 +
                             tk_yyyy_yyyzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_yyyzzz[i] * fz_0;

        tk_yyyyy_yyzzzz[i] = -8.0 * ts_yyy_yyzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yyzzzz[i] * fe_0 + 2.0 * tk_yyyy_yzzzz[i] * fe_0 +
                             tk_yyyy_yyzzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_yyzzzz[i] * fz_0;

        tk_yyyyy_yzzzzz[i] = -8.0 * ts_yyy_yzzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yzzzzz[i] * fe_0 + tk_yyyy_zzzzz[i] * fe_0 +
                             tk_yyyy_yzzzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_yzzzzz[i] * fz_0;

        tk_yyyyy_zzzzzz[i] =
            -8.0 * ts_yyy_zzzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_zzzzzz[i] * fe_0 + tk_yyyy_zzzzzz[i] * pa_y[i] + 2.0 * ts_yyyyy_zzzzzz[i] * fz_0;
    }

    // Set up 448-476 components of targeted buffer : HI

    auto tk_yyyyz_xxxxxx = pbuffer.data(idx_kin_hi + 448);

    auto tk_yyyyz_xxxxxy = pbuffer.data(idx_kin_hi + 449);

    auto tk_yyyyz_xxxxxz = pbuffer.data(idx_kin_hi + 450);

    auto tk_yyyyz_xxxxyy = pbuffer.data(idx_kin_hi + 451);

    auto tk_yyyyz_xxxxyz = pbuffer.data(idx_kin_hi + 452);

    auto tk_yyyyz_xxxxzz = pbuffer.data(idx_kin_hi + 453);

    auto tk_yyyyz_xxxyyy = pbuffer.data(idx_kin_hi + 454);

    auto tk_yyyyz_xxxyyz = pbuffer.data(idx_kin_hi + 455);

    auto tk_yyyyz_xxxyzz = pbuffer.data(idx_kin_hi + 456);

    auto tk_yyyyz_xxxzzz = pbuffer.data(idx_kin_hi + 457);

    auto tk_yyyyz_xxyyyy = pbuffer.data(idx_kin_hi + 458);

    auto tk_yyyyz_xxyyyz = pbuffer.data(idx_kin_hi + 459);

    auto tk_yyyyz_xxyyzz = pbuffer.data(idx_kin_hi + 460);

    auto tk_yyyyz_xxyzzz = pbuffer.data(idx_kin_hi + 461);

    auto tk_yyyyz_xxzzzz = pbuffer.data(idx_kin_hi + 462);

    auto tk_yyyyz_xyyyyy = pbuffer.data(idx_kin_hi + 463);

    auto tk_yyyyz_xyyyyz = pbuffer.data(idx_kin_hi + 464);

    auto tk_yyyyz_xyyyzz = pbuffer.data(idx_kin_hi + 465);

    auto tk_yyyyz_xyyzzz = pbuffer.data(idx_kin_hi + 466);

    auto tk_yyyyz_xyzzzz = pbuffer.data(idx_kin_hi + 467);

    auto tk_yyyyz_xzzzzz = pbuffer.data(idx_kin_hi + 468);

    auto tk_yyyyz_yyyyyy = pbuffer.data(idx_kin_hi + 469);

    auto tk_yyyyz_yyyyyz = pbuffer.data(idx_kin_hi + 470);

    auto tk_yyyyz_yyyyzz = pbuffer.data(idx_kin_hi + 471);

    auto tk_yyyyz_yyyzzz = pbuffer.data(idx_kin_hi + 472);

    auto tk_yyyyz_yyzzzz = pbuffer.data(idx_kin_hi + 473);

    auto tk_yyyyz_yzzzzz = pbuffer.data(idx_kin_hi + 474);

    auto tk_yyyyz_zzzzzz = pbuffer.data(idx_kin_hi + 475);

#pragma omp simd aligned(pa_z,                \
                             tk_yyyy_xxxxx,   \
                             tk_yyyy_xxxxxx,  \
                             tk_yyyy_xxxxxy,  \
                             tk_yyyy_xxxxxz,  \
                             tk_yyyy_xxxxy,   \
                             tk_yyyy_xxxxyy,  \
                             tk_yyyy_xxxxyz,  \
                             tk_yyyy_xxxxz,   \
                             tk_yyyy_xxxxzz,  \
                             tk_yyyy_xxxyy,   \
                             tk_yyyy_xxxyyy,  \
                             tk_yyyy_xxxyyz,  \
                             tk_yyyy_xxxyz,   \
                             tk_yyyy_xxxyzz,  \
                             tk_yyyy_xxxzz,   \
                             tk_yyyy_xxxzzz,  \
                             tk_yyyy_xxyyy,   \
                             tk_yyyy_xxyyyy,  \
                             tk_yyyy_xxyyyz,  \
                             tk_yyyy_xxyyz,   \
                             tk_yyyy_xxyyzz,  \
                             tk_yyyy_xxyzz,   \
                             tk_yyyy_xxyzzz,  \
                             tk_yyyy_xxzzz,   \
                             tk_yyyy_xxzzzz,  \
                             tk_yyyy_xyyyy,   \
                             tk_yyyy_xyyyyy,  \
                             tk_yyyy_xyyyyz,  \
                             tk_yyyy_xyyyz,   \
                             tk_yyyy_xyyyzz,  \
                             tk_yyyy_xyyzz,   \
                             tk_yyyy_xyyzzz,  \
                             tk_yyyy_xyzzz,   \
                             tk_yyyy_xyzzzz,  \
                             tk_yyyy_xzzzz,   \
                             tk_yyyy_xzzzzz,  \
                             tk_yyyy_yyyyy,   \
                             tk_yyyy_yyyyyy,  \
                             tk_yyyy_yyyyyz,  \
                             tk_yyyy_yyyyz,   \
                             tk_yyyy_yyyyzz,  \
                             tk_yyyy_yyyzz,   \
                             tk_yyyy_yyyzzz,  \
                             tk_yyyy_yyzzz,   \
                             tk_yyyy_yyzzzz,  \
                             tk_yyyy_yzzzz,   \
                             tk_yyyy_yzzzzz,  \
                             tk_yyyy_zzzzz,   \
                             tk_yyyy_zzzzzz,  \
                             tk_yyyyz_xxxxxx, \
                             tk_yyyyz_xxxxxy, \
                             tk_yyyyz_xxxxxz, \
                             tk_yyyyz_xxxxyy, \
                             tk_yyyyz_xxxxyz, \
                             tk_yyyyz_xxxxzz, \
                             tk_yyyyz_xxxyyy, \
                             tk_yyyyz_xxxyyz, \
                             tk_yyyyz_xxxyzz, \
                             tk_yyyyz_xxxzzz, \
                             tk_yyyyz_xxyyyy, \
                             tk_yyyyz_xxyyyz, \
                             tk_yyyyz_xxyyzz, \
                             tk_yyyyz_xxyzzz, \
                             tk_yyyyz_xxzzzz, \
                             tk_yyyyz_xyyyyy, \
                             tk_yyyyz_xyyyyz, \
                             tk_yyyyz_xyyyzz, \
                             tk_yyyyz_xyyzzz, \
                             tk_yyyyz_xyzzzz, \
                             tk_yyyyz_xzzzzz, \
                             tk_yyyyz_yyyyyy, \
                             tk_yyyyz_yyyyyz, \
                             tk_yyyyz_yyyyzz, \
                             tk_yyyyz_yyyzzz, \
                             tk_yyyyz_yyzzzz, \
                             tk_yyyyz_yzzzzz, \
                             tk_yyyyz_zzzzzz, \
                             ts_yyyyz_xxxxxx, \
                             ts_yyyyz_xxxxxy, \
                             ts_yyyyz_xxxxxz, \
                             ts_yyyyz_xxxxyy, \
                             ts_yyyyz_xxxxyz, \
                             ts_yyyyz_xxxxzz, \
                             ts_yyyyz_xxxyyy, \
                             ts_yyyyz_xxxyyz, \
                             ts_yyyyz_xxxyzz, \
                             ts_yyyyz_xxxzzz, \
                             ts_yyyyz_xxyyyy, \
                             ts_yyyyz_xxyyyz, \
                             ts_yyyyz_xxyyzz, \
                             ts_yyyyz_xxyzzz, \
                             ts_yyyyz_xxzzzz, \
                             ts_yyyyz_xyyyyy, \
                             ts_yyyyz_xyyyyz, \
                             ts_yyyyz_xyyyzz, \
                             ts_yyyyz_xyyzzz, \
                             ts_yyyyz_xyzzzz, \
                             ts_yyyyz_xzzzzz, \
                             ts_yyyyz_yyyyyy, \
                             ts_yyyyz_yyyyyz, \
                             ts_yyyyz_yyyyzz, \
                             ts_yyyyz_yyyzzz, \
                             ts_yyyyz_yyzzzz, \
                             ts_yyyyz_yzzzzz, \
                             ts_yyyyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyyz_xxxxxx[i] = tk_yyyy_xxxxxx[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxxxx[i] * fz_0;

        tk_yyyyz_xxxxxy[i] = tk_yyyy_xxxxxy[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxxxy[i] * fz_0;

        tk_yyyyz_xxxxxz[i] = tk_yyyy_xxxxx[i] * fe_0 + tk_yyyy_xxxxxz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxxxz[i] * fz_0;

        tk_yyyyz_xxxxyy[i] = tk_yyyy_xxxxyy[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxxyy[i] * fz_0;

        tk_yyyyz_xxxxyz[i] = tk_yyyy_xxxxy[i] * fe_0 + tk_yyyy_xxxxyz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxxyz[i] * fz_0;

        tk_yyyyz_xxxxzz[i] = 2.0 * tk_yyyy_xxxxz[i] * fe_0 + tk_yyyy_xxxxzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxxzz[i] * fz_0;

        tk_yyyyz_xxxyyy[i] = tk_yyyy_xxxyyy[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxyyy[i] * fz_0;

        tk_yyyyz_xxxyyz[i] = tk_yyyy_xxxyy[i] * fe_0 + tk_yyyy_xxxyyz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxyyz[i] * fz_0;

        tk_yyyyz_xxxyzz[i] = 2.0 * tk_yyyy_xxxyz[i] * fe_0 + tk_yyyy_xxxyzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxyzz[i] * fz_0;

        tk_yyyyz_xxxzzz[i] = 3.0 * tk_yyyy_xxxzz[i] * fe_0 + tk_yyyy_xxxzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxxzzz[i] * fz_0;

        tk_yyyyz_xxyyyy[i] = tk_yyyy_xxyyyy[i] * pa_z[i] + 2.0 * ts_yyyyz_xxyyyy[i] * fz_0;

        tk_yyyyz_xxyyyz[i] = tk_yyyy_xxyyy[i] * fe_0 + tk_yyyy_xxyyyz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxyyyz[i] * fz_0;

        tk_yyyyz_xxyyzz[i] = 2.0 * tk_yyyy_xxyyz[i] * fe_0 + tk_yyyy_xxyyzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxyyzz[i] * fz_0;

        tk_yyyyz_xxyzzz[i] = 3.0 * tk_yyyy_xxyzz[i] * fe_0 + tk_yyyy_xxyzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxyzzz[i] * fz_0;

        tk_yyyyz_xxzzzz[i] = 4.0 * tk_yyyy_xxzzz[i] * fe_0 + tk_yyyy_xxzzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xxzzzz[i] * fz_0;

        tk_yyyyz_xyyyyy[i] = tk_yyyy_xyyyyy[i] * pa_z[i] + 2.0 * ts_yyyyz_xyyyyy[i] * fz_0;

        tk_yyyyz_xyyyyz[i] = tk_yyyy_xyyyy[i] * fe_0 + tk_yyyy_xyyyyz[i] * pa_z[i] + 2.0 * ts_yyyyz_xyyyyz[i] * fz_0;

        tk_yyyyz_xyyyzz[i] = 2.0 * tk_yyyy_xyyyz[i] * fe_0 + tk_yyyy_xyyyzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xyyyzz[i] * fz_0;

        tk_yyyyz_xyyzzz[i] = 3.0 * tk_yyyy_xyyzz[i] * fe_0 + tk_yyyy_xyyzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xyyzzz[i] * fz_0;

        tk_yyyyz_xyzzzz[i] = 4.0 * tk_yyyy_xyzzz[i] * fe_0 + tk_yyyy_xyzzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xyzzzz[i] * fz_0;

        tk_yyyyz_xzzzzz[i] = 5.0 * tk_yyyy_xzzzz[i] * fe_0 + tk_yyyy_xzzzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_xzzzzz[i] * fz_0;

        tk_yyyyz_yyyyyy[i] = tk_yyyy_yyyyyy[i] * pa_z[i] + 2.0 * ts_yyyyz_yyyyyy[i] * fz_0;

        tk_yyyyz_yyyyyz[i] = tk_yyyy_yyyyy[i] * fe_0 + tk_yyyy_yyyyyz[i] * pa_z[i] + 2.0 * ts_yyyyz_yyyyyz[i] * fz_0;

        tk_yyyyz_yyyyzz[i] = 2.0 * tk_yyyy_yyyyz[i] * fe_0 + tk_yyyy_yyyyzz[i] * pa_z[i] + 2.0 * ts_yyyyz_yyyyzz[i] * fz_0;

        tk_yyyyz_yyyzzz[i] = 3.0 * tk_yyyy_yyyzz[i] * fe_0 + tk_yyyy_yyyzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_yyyzzz[i] * fz_0;

        tk_yyyyz_yyzzzz[i] = 4.0 * tk_yyyy_yyzzz[i] * fe_0 + tk_yyyy_yyzzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_yyzzzz[i] * fz_0;

        tk_yyyyz_yzzzzz[i] = 5.0 * tk_yyyy_yzzzz[i] * fe_0 + tk_yyyy_yzzzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_yzzzzz[i] * fz_0;

        tk_yyyyz_zzzzzz[i] = 6.0 * tk_yyyy_zzzzz[i] * fe_0 + tk_yyyy_zzzzzz[i] * pa_z[i] + 2.0 * ts_yyyyz_zzzzzz[i] * fz_0;
    }

    // Set up 476-504 components of targeted buffer : HI

    auto tk_yyyzz_xxxxxx = pbuffer.data(idx_kin_hi + 476);

    auto tk_yyyzz_xxxxxy = pbuffer.data(idx_kin_hi + 477);

    auto tk_yyyzz_xxxxxz = pbuffer.data(idx_kin_hi + 478);

    auto tk_yyyzz_xxxxyy = pbuffer.data(idx_kin_hi + 479);

    auto tk_yyyzz_xxxxyz = pbuffer.data(idx_kin_hi + 480);

    auto tk_yyyzz_xxxxzz = pbuffer.data(idx_kin_hi + 481);

    auto tk_yyyzz_xxxyyy = pbuffer.data(idx_kin_hi + 482);

    auto tk_yyyzz_xxxyyz = pbuffer.data(idx_kin_hi + 483);

    auto tk_yyyzz_xxxyzz = pbuffer.data(idx_kin_hi + 484);

    auto tk_yyyzz_xxxzzz = pbuffer.data(idx_kin_hi + 485);

    auto tk_yyyzz_xxyyyy = pbuffer.data(idx_kin_hi + 486);

    auto tk_yyyzz_xxyyyz = pbuffer.data(idx_kin_hi + 487);

    auto tk_yyyzz_xxyyzz = pbuffer.data(idx_kin_hi + 488);

    auto tk_yyyzz_xxyzzz = pbuffer.data(idx_kin_hi + 489);

    auto tk_yyyzz_xxzzzz = pbuffer.data(idx_kin_hi + 490);

    auto tk_yyyzz_xyyyyy = pbuffer.data(idx_kin_hi + 491);

    auto tk_yyyzz_xyyyyz = pbuffer.data(idx_kin_hi + 492);

    auto tk_yyyzz_xyyyzz = pbuffer.data(idx_kin_hi + 493);

    auto tk_yyyzz_xyyzzz = pbuffer.data(idx_kin_hi + 494);

    auto tk_yyyzz_xyzzzz = pbuffer.data(idx_kin_hi + 495);

    auto tk_yyyzz_xzzzzz = pbuffer.data(idx_kin_hi + 496);

    auto tk_yyyzz_yyyyyy = pbuffer.data(idx_kin_hi + 497);

    auto tk_yyyzz_yyyyyz = pbuffer.data(idx_kin_hi + 498);

    auto tk_yyyzz_yyyyzz = pbuffer.data(idx_kin_hi + 499);

    auto tk_yyyzz_yyyzzz = pbuffer.data(idx_kin_hi + 500);

    auto tk_yyyzz_yyzzzz = pbuffer.data(idx_kin_hi + 501);

    auto tk_yyyzz_yzzzzz = pbuffer.data(idx_kin_hi + 502);

    auto tk_yyyzz_zzzzzz = pbuffer.data(idx_kin_hi + 503);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tk_yyy_xxxxxy,   \
                             tk_yyy_xxxxyy,   \
                             tk_yyy_xxxyyy,   \
                             tk_yyy_xxyyyy,   \
                             tk_yyy_xyyyyy,   \
                             tk_yyy_yyyyyy,   \
                             tk_yyyz_xxxxxy,  \
                             tk_yyyz_xxxxyy,  \
                             tk_yyyz_xxxyyy,  \
                             tk_yyyz_xxyyyy,  \
                             tk_yyyz_xyyyyy,  \
                             tk_yyyz_yyyyyy,  \
                             tk_yyyzz_xxxxxx, \
                             tk_yyyzz_xxxxxy, \
                             tk_yyyzz_xxxxxz, \
                             tk_yyyzz_xxxxyy, \
                             tk_yyyzz_xxxxyz, \
                             tk_yyyzz_xxxxzz, \
                             tk_yyyzz_xxxyyy, \
                             tk_yyyzz_xxxyyz, \
                             tk_yyyzz_xxxyzz, \
                             tk_yyyzz_xxxzzz, \
                             tk_yyyzz_xxyyyy, \
                             tk_yyyzz_xxyyyz, \
                             tk_yyyzz_xxyyzz, \
                             tk_yyyzz_xxyzzz, \
                             tk_yyyzz_xxzzzz, \
                             tk_yyyzz_xyyyyy, \
                             tk_yyyzz_xyyyyz, \
                             tk_yyyzz_xyyyzz, \
                             tk_yyyzz_xyyzzz, \
                             tk_yyyzz_xyzzzz, \
                             tk_yyyzz_xzzzzz, \
                             tk_yyyzz_yyyyyy, \
                             tk_yyyzz_yyyyyz, \
                             tk_yyyzz_yyyyzz, \
                             tk_yyyzz_yyyzzz, \
                             tk_yyyzz_yyzzzz, \
                             tk_yyyzz_yzzzzz, \
                             tk_yyyzz_zzzzzz, \
                             tk_yyzz_xxxxxx,  \
                             tk_yyzz_xxxxxz,  \
                             tk_yyzz_xxxxyz,  \
                             tk_yyzz_xxxxz,   \
                             tk_yyzz_xxxxzz,  \
                             tk_yyzz_xxxyyz,  \
                             tk_yyzz_xxxyz,   \
                             tk_yyzz_xxxyzz,  \
                             tk_yyzz_xxxzz,   \
                             tk_yyzz_xxxzzz,  \
                             tk_yyzz_xxyyyz,  \
                             tk_yyzz_xxyyz,   \
                             tk_yyzz_xxyyzz,  \
                             tk_yyzz_xxyzz,   \
                             tk_yyzz_xxyzzz,  \
                             tk_yyzz_xxzzz,   \
                             tk_yyzz_xxzzzz,  \
                             tk_yyzz_xyyyyz,  \
                             tk_yyzz_xyyyz,   \
                             tk_yyzz_xyyyzz,  \
                             tk_yyzz_xyyzz,   \
                             tk_yyzz_xyyzzz,  \
                             tk_yyzz_xyzzz,   \
                             tk_yyzz_xyzzzz,  \
                             tk_yyzz_xzzzz,   \
                             tk_yyzz_xzzzzz,  \
                             tk_yyzz_yyyyyz,  \
                             tk_yyzz_yyyyz,   \
                             tk_yyzz_yyyyzz,  \
                             tk_yyzz_yyyzz,   \
                             tk_yyzz_yyyzzz,  \
                             tk_yyzz_yyzzz,   \
                             tk_yyzz_yyzzzz,  \
                             tk_yyzz_yzzzz,   \
                             tk_yyzz_yzzzzz,  \
                             tk_yyzz_zzzzz,   \
                             tk_yyzz_zzzzzz,  \
                             tk_yzz_xxxxxx,   \
                             tk_yzz_xxxxxz,   \
                             tk_yzz_xxxxyz,   \
                             tk_yzz_xxxxzz,   \
                             tk_yzz_xxxyyz,   \
                             tk_yzz_xxxyzz,   \
                             tk_yzz_xxxzzz,   \
                             tk_yzz_xxyyyz,   \
                             tk_yzz_xxyyzz,   \
                             tk_yzz_xxyzzz,   \
                             tk_yzz_xxzzzz,   \
                             tk_yzz_xyyyyz,   \
                             tk_yzz_xyyyzz,   \
                             tk_yzz_xyyzzz,   \
                             tk_yzz_xyzzzz,   \
                             tk_yzz_xzzzzz,   \
                             tk_yzz_yyyyyz,   \
                             tk_yzz_yyyyzz,   \
                             tk_yzz_yyyzzz,   \
                             tk_yzz_yyzzzz,   \
                             tk_yzz_yzzzzz,   \
                             tk_yzz_zzzzzz,   \
                             ts_yyy_xxxxxy,   \
                             ts_yyy_xxxxyy,   \
                             ts_yyy_xxxyyy,   \
                             ts_yyy_xxyyyy,   \
                             ts_yyy_xyyyyy,   \
                             ts_yyy_yyyyyy,   \
                             ts_yyyzz_xxxxxx, \
                             ts_yyyzz_xxxxxy, \
                             ts_yyyzz_xxxxxz, \
                             ts_yyyzz_xxxxyy, \
                             ts_yyyzz_xxxxyz, \
                             ts_yyyzz_xxxxzz, \
                             ts_yyyzz_xxxyyy, \
                             ts_yyyzz_xxxyyz, \
                             ts_yyyzz_xxxyzz, \
                             ts_yyyzz_xxxzzz, \
                             ts_yyyzz_xxyyyy, \
                             ts_yyyzz_xxyyyz, \
                             ts_yyyzz_xxyyzz, \
                             ts_yyyzz_xxyzzz, \
                             ts_yyyzz_xxzzzz, \
                             ts_yyyzz_xyyyyy, \
                             ts_yyyzz_xyyyyz, \
                             ts_yyyzz_xyyyzz, \
                             ts_yyyzz_xyyzzz, \
                             ts_yyyzz_xyzzzz, \
                             ts_yyyzz_xzzzzz, \
                             ts_yyyzz_yyyyyy, \
                             ts_yyyzz_yyyyyz, \
                             ts_yyyzz_yyyyzz, \
                             ts_yyyzz_yyyzzz, \
                             ts_yyyzz_yyzzzz, \
                             ts_yyyzz_yzzzzz, \
                             ts_yyyzz_zzzzzz, \
                             ts_yzz_xxxxxx,   \
                             ts_yzz_xxxxxz,   \
                             ts_yzz_xxxxyz,   \
                             ts_yzz_xxxxzz,   \
                             ts_yzz_xxxyyz,   \
                             ts_yzz_xxxyzz,   \
                             ts_yzz_xxxzzz,   \
                             ts_yzz_xxyyyz,   \
                             ts_yzz_xxyyzz,   \
                             ts_yzz_xxyzzz,   \
                             ts_yzz_xxzzzz,   \
                             ts_yzz_xyyyyz,   \
                             ts_yzz_xyyyzz,   \
                             ts_yzz_xyyzzz,   \
                             ts_yzz_xyzzzz,   \
                             ts_yzz_xzzzzz,   \
                             ts_yzz_yyyyyz,   \
                             ts_yzz_yyyyzz,   \
                             ts_yzz_yyyzzz,   \
                             ts_yzz_yyzzzz,   \
                             ts_yzz_yzzzzz,   \
                             ts_yzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyzz_xxxxxx[i] =
            -4.0 * ts_yzz_xxxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxxxx[i] * fe_0 + tk_yyzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxxxx[i] * fz_0;

        tk_yyyzz_xxxxxy[i] =
            -2.0 * ts_yyy_xxxxxy[i] * fbe_0 * fz_0 + tk_yyy_xxxxxy[i] * fe_0 + tk_yyyz_xxxxxy[i] * pa_z[i] + 2.0 * ts_yyyzz_xxxxxy[i] * fz_0;

        tk_yyyzz_xxxxxz[i] =
            -4.0 * ts_yzz_xxxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxxxz[i] * fe_0 + tk_yyzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxxxz[i] * fz_0;

        tk_yyyzz_xxxxyy[i] =
            -2.0 * ts_yyy_xxxxyy[i] * fbe_0 * fz_0 + tk_yyy_xxxxyy[i] * fe_0 + tk_yyyz_xxxxyy[i] * pa_z[i] + 2.0 * ts_yyyzz_xxxxyy[i] * fz_0;

        tk_yyyzz_xxxxyz[i] = -4.0 * ts_yzz_xxxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxxyz[i] * fe_0 + tk_yyzz_xxxxz[i] * fe_0 +
                             tk_yyzz_xxxxyz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxxyz[i] * fz_0;

        tk_yyyzz_xxxxzz[i] =
            -4.0 * ts_yzz_xxxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxxzz[i] * fe_0 + tk_yyzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxxzz[i] * fz_0;

        tk_yyyzz_xxxyyy[i] =
            -2.0 * ts_yyy_xxxyyy[i] * fbe_0 * fz_0 + tk_yyy_xxxyyy[i] * fe_0 + tk_yyyz_xxxyyy[i] * pa_z[i] + 2.0 * ts_yyyzz_xxxyyy[i] * fz_0;

        tk_yyyzz_xxxyyz[i] = -4.0 * ts_yzz_xxxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxyyz[i] * fe_0 + 2.0 * tk_yyzz_xxxyz[i] * fe_0 +
                             tk_yyzz_xxxyyz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxyyz[i] * fz_0;

        tk_yyyzz_xxxyzz[i] = -4.0 * ts_yzz_xxxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxyzz[i] * fe_0 + tk_yyzz_xxxzz[i] * fe_0 +
                             tk_yyzz_xxxyzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxyzz[i] * fz_0;

        tk_yyyzz_xxxzzz[i] =
            -4.0 * ts_yzz_xxxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxxzzz[i] * fe_0 + tk_yyzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxxzzz[i] * fz_0;

        tk_yyyzz_xxyyyy[i] =
            -2.0 * ts_yyy_xxyyyy[i] * fbe_0 * fz_0 + tk_yyy_xxyyyy[i] * fe_0 + tk_yyyz_xxyyyy[i] * pa_z[i] + 2.0 * ts_yyyzz_xxyyyy[i] * fz_0;

        tk_yyyzz_xxyyyz[i] = -4.0 * ts_yzz_xxyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxyyyz[i] * fe_0 + 3.0 * tk_yyzz_xxyyz[i] * fe_0 +
                             tk_yyzz_xxyyyz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxyyyz[i] * fz_0;

        tk_yyyzz_xxyyzz[i] = -4.0 * ts_yzz_xxyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxyyzz[i] * fe_0 + 2.0 * tk_yyzz_xxyzz[i] * fe_0 +
                             tk_yyzz_xxyyzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxyyzz[i] * fz_0;

        tk_yyyzz_xxyzzz[i] = -4.0 * ts_yzz_xxyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxyzzz[i] * fe_0 + tk_yyzz_xxzzz[i] * fe_0 +
                             tk_yyzz_xxyzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxyzzz[i] * fz_0;

        tk_yyyzz_xxzzzz[i] =
            -4.0 * ts_yzz_xxzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xxzzzz[i] * fe_0 + tk_yyzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xxzzzz[i] * fz_0;

        tk_yyyzz_xyyyyy[i] =
            -2.0 * ts_yyy_xyyyyy[i] * fbe_0 * fz_0 + tk_yyy_xyyyyy[i] * fe_0 + tk_yyyz_xyyyyy[i] * pa_z[i] + 2.0 * ts_yyyzz_xyyyyy[i] * fz_0;

        tk_yyyzz_xyyyyz[i] = -4.0 * ts_yzz_xyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xyyyyz[i] * fe_0 + 4.0 * tk_yyzz_xyyyz[i] * fe_0 +
                             tk_yyzz_xyyyyz[i] * pa_y[i] + 2.0 * ts_yyyzz_xyyyyz[i] * fz_0;

        tk_yyyzz_xyyyzz[i] = -4.0 * ts_yzz_xyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xyyyzz[i] * fe_0 + 3.0 * tk_yyzz_xyyzz[i] * fe_0 +
                             tk_yyzz_xyyyzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xyyyzz[i] * fz_0;

        tk_yyyzz_xyyzzz[i] = -4.0 * ts_yzz_xyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xyyzzz[i] * fe_0 + 2.0 * tk_yyzz_xyzzz[i] * fe_0 +
                             tk_yyzz_xyyzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xyyzzz[i] * fz_0;

        tk_yyyzz_xyzzzz[i] = -4.0 * ts_yzz_xyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xyzzzz[i] * fe_0 + tk_yyzz_xzzzz[i] * fe_0 +
                             tk_yyzz_xyzzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xyzzzz[i] * fz_0;

        tk_yyyzz_xzzzzz[i] =
            -4.0 * ts_yzz_xzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xzzzzz[i] * fe_0 + tk_yyzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_xzzzzz[i] * fz_0;

        tk_yyyzz_yyyyyy[i] =
            -2.0 * ts_yyy_yyyyyy[i] * fbe_0 * fz_0 + tk_yyy_yyyyyy[i] * fe_0 + tk_yyyz_yyyyyy[i] * pa_z[i] + 2.0 * ts_yyyzz_yyyyyy[i] * fz_0;

        tk_yyyzz_yyyyyz[i] = -4.0 * ts_yzz_yyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yyyyyz[i] * fe_0 + 5.0 * tk_yyzz_yyyyz[i] * fe_0 +
                             tk_yyzz_yyyyyz[i] * pa_y[i] + 2.0 * ts_yyyzz_yyyyyz[i] * fz_0;

        tk_yyyzz_yyyyzz[i] = -4.0 * ts_yzz_yyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yyyyzz[i] * fe_0 + 4.0 * tk_yyzz_yyyzz[i] * fe_0 +
                             tk_yyzz_yyyyzz[i] * pa_y[i] + 2.0 * ts_yyyzz_yyyyzz[i] * fz_0;

        tk_yyyzz_yyyzzz[i] = -4.0 * ts_yzz_yyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yyyzzz[i] * fe_0 + 3.0 * tk_yyzz_yyzzz[i] * fe_0 +
                             tk_yyzz_yyyzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_yyyzzz[i] * fz_0;

        tk_yyyzz_yyzzzz[i] = -4.0 * ts_yzz_yyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yyzzzz[i] * fe_0 + 2.0 * tk_yyzz_yzzzz[i] * fe_0 +
                             tk_yyzz_yyzzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_yyzzzz[i] * fz_0;

        tk_yyyzz_yzzzzz[i] = -4.0 * ts_yzz_yzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yzzzzz[i] * fe_0 + tk_yyzz_zzzzz[i] * fe_0 +
                             tk_yyzz_yzzzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_yzzzzz[i] * fz_0;

        tk_yyyzz_zzzzzz[i] =
            -4.0 * ts_yzz_zzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_zzzzzz[i] * fe_0 + tk_yyzz_zzzzzz[i] * pa_y[i] + 2.0 * ts_yyyzz_zzzzzz[i] * fz_0;
    }

    // Set up 504-532 components of targeted buffer : HI

    auto tk_yyzzz_xxxxxx = pbuffer.data(idx_kin_hi + 504);

    auto tk_yyzzz_xxxxxy = pbuffer.data(idx_kin_hi + 505);

    auto tk_yyzzz_xxxxxz = pbuffer.data(idx_kin_hi + 506);

    auto tk_yyzzz_xxxxyy = pbuffer.data(idx_kin_hi + 507);

    auto tk_yyzzz_xxxxyz = pbuffer.data(idx_kin_hi + 508);

    auto tk_yyzzz_xxxxzz = pbuffer.data(idx_kin_hi + 509);

    auto tk_yyzzz_xxxyyy = pbuffer.data(idx_kin_hi + 510);

    auto tk_yyzzz_xxxyyz = pbuffer.data(idx_kin_hi + 511);

    auto tk_yyzzz_xxxyzz = pbuffer.data(idx_kin_hi + 512);

    auto tk_yyzzz_xxxzzz = pbuffer.data(idx_kin_hi + 513);

    auto tk_yyzzz_xxyyyy = pbuffer.data(idx_kin_hi + 514);

    auto tk_yyzzz_xxyyyz = pbuffer.data(idx_kin_hi + 515);

    auto tk_yyzzz_xxyyzz = pbuffer.data(idx_kin_hi + 516);

    auto tk_yyzzz_xxyzzz = pbuffer.data(idx_kin_hi + 517);

    auto tk_yyzzz_xxzzzz = pbuffer.data(idx_kin_hi + 518);

    auto tk_yyzzz_xyyyyy = pbuffer.data(idx_kin_hi + 519);

    auto tk_yyzzz_xyyyyz = pbuffer.data(idx_kin_hi + 520);

    auto tk_yyzzz_xyyyzz = pbuffer.data(idx_kin_hi + 521);

    auto tk_yyzzz_xyyzzz = pbuffer.data(idx_kin_hi + 522);

    auto tk_yyzzz_xyzzzz = pbuffer.data(idx_kin_hi + 523);

    auto tk_yyzzz_xzzzzz = pbuffer.data(idx_kin_hi + 524);

    auto tk_yyzzz_yyyyyy = pbuffer.data(idx_kin_hi + 525);

    auto tk_yyzzz_yyyyyz = pbuffer.data(idx_kin_hi + 526);

    auto tk_yyzzz_yyyyzz = pbuffer.data(idx_kin_hi + 527);

    auto tk_yyzzz_yyyzzz = pbuffer.data(idx_kin_hi + 528);

    auto tk_yyzzz_yyzzzz = pbuffer.data(idx_kin_hi + 529);

    auto tk_yyzzz_yzzzzz = pbuffer.data(idx_kin_hi + 530);

    auto tk_yyzzz_zzzzzz = pbuffer.data(idx_kin_hi + 531);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tk_yyz_xxxxxy,   \
                             tk_yyz_xxxxyy,   \
                             tk_yyz_xxxyyy,   \
                             tk_yyz_xxyyyy,   \
                             tk_yyz_xyyyyy,   \
                             tk_yyz_yyyyyy,   \
                             tk_yyzz_xxxxxy,  \
                             tk_yyzz_xxxxyy,  \
                             tk_yyzz_xxxyyy,  \
                             tk_yyzz_xxyyyy,  \
                             tk_yyzz_xyyyyy,  \
                             tk_yyzz_yyyyyy,  \
                             tk_yyzzz_xxxxxx, \
                             tk_yyzzz_xxxxxy, \
                             tk_yyzzz_xxxxxz, \
                             tk_yyzzz_xxxxyy, \
                             tk_yyzzz_xxxxyz, \
                             tk_yyzzz_xxxxzz, \
                             tk_yyzzz_xxxyyy, \
                             tk_yyzzz_xxxyyz, \
                             tk_yyzzz_xxxyzz, \
                             tk_yyzzz_xxxzzz, \
                             tk_yyzzz_xxyyyy, \
                             tk_yyzzz_xxyyyz, \
                             tk_yyzzz_xxyyzz, \
                             tk_yyzzz_xxyzzz, \
                             tk_yyzzz_xxzzzz, \
                             tk_yyzzz_xyyyyy, \
                             tk_yyzzz_xyyyyz, \
                             tk_yyzzz_xyyyzz, \
                             tk_yyzzz_xyyzzz, \
                             tk_yyzzz_xyzzzz, \
                             tk_yyzzz_xzzzzz, \
                             tk_yyzzz_yyyyyy, \
                             tk_yyzzz_yyyyyz, \
                             tk_yyzzz_yyyyzz, \
                             tk_yyzzz_yyyzzz, \
                             tk_yyzzz_yyzzzz, \
                             tk_yyzzz_yzzzzz, \
                             tk_yyzzz_zzzzzz, \
                             tk_yzzz_xxxxxx,  \
                             tk_yzzz_xxxxxz,  \
                             tk_yzzz_xxxxyz,  \
                             tk_yzzz_xxxxz,   \
                             tk_yzzz_xxxxzz,  \
                             tk_yzzz_xxxyyz,  \
                             tk_yzzz_xxxyz,   \
                             tk_yzzz_xxxyzz,  \
                             tk_yzzz_xxxzz,   \
                             tk_yzzz_xxxzzz,  \
                             tk_yzzz_xxyyyz,  \
                             tk_yzzz_xxyyz,   \
                             tk_yzzz_xxyyzz,  \
                             tk_yzzz_xxyzz,   \
                             tk_yzzz_xxyzzz,  \
                             tk_yzzz_xxzzz,   \
                             tk_yzzz_xxzzzz,  \
                             tk_yzzz_xyyyyz,  \
                             tk_yzzz_xyyyz,   \
                             tk_yzzz_xyyyzz,  \
                             tk_yzzz_xyyzz,   \
                             tk_yzzz_xyyzzz,  \
                             tk_yzzz_xyzzz,   \
                             tk_yzzz_xyzzzz,  \
                             tk_yzzz_xzzzz,   \
                             tk_yzzz_xzzzzz,  \
                             tk_yzzz_yyyyyz,  \
                             tk_yzzz_yyyyz,   \
                             tk_yzzz_yyyyzz,  \
                             tk_yzzz_yyyzz,   \
                             tk_yzzz_yyyzzz,  \
                             tk_yzzz_yyzzz,   \
                             tk_yzzz_yyzzzz,  \
                             tk_yzzz_yzzzz,   \
                             tk_yzzz_yzzzzz,  \
                             tk_yzzz_zzzzz,   \
                             tk_yzzz_zzzzzz,  \
                             tk_zzz_xxxxxx,   \
                             tk_zzz_xxxxxz,   \
                             tk_zzz_xxxxyz,   \
                             tk_zzz_xxxxzz,   \
                             tk_zzz_xxxyyz,   \
                             tk_zzz_xxxyzz,   \
                             tk_zzz_xxxzzz,   \
                             tk_zzz_xxyyyz,   \
                             tk_zzz_xxyyzz,   \
                             tk_zzz_xxyzzz,   \
                             tk_zzz_xxzzzz,   \
                             tk_zzz_xyyyyz,   \
                             tk_zzz_xyyyzz,   \
                             tk_zzz_xyyzzz,   \
                             tk_zzz_xyzzzz,   \
                             tk_zzz_xzzzzz,   \
                             tk_zzz_yyyyyz,   \
                             tk_zzz_yyyyzz,   \
                             tk_zzz_yyyzzz,   \
                             tk_zzz_yyzzzz,   \
                             tk_zzz_yzzzzz,   \
                             tk_zzz_zzzzzz,   \
                             ts_yyz_xxxxxy,   \
                             ts_yyz_xxxxyy,   \
                             ts_yyz_xxxyyy,   \
                             ts_yyz_xxyyyy,   \
                             ts_yyz_xyyyyy,   \
                             ts_yyz_yyyyyy,   \
                             ts_yyzzz_xxxxxx, \
                             ts_yyzzz_xxxxxy, \
                             ts_yyzzz_xxxxxz, \
                             ts_yyzzz_xxxxyy, \
                             ts_yyzzz_xxxxyz, \
                             ts_yyzzz_xxxxzz, \
                             ts_yyzzz_xxxyyy, \
                             ts_yyzzz_xxxyyz, \
                             ts_yyzzz_xxxyzz, \
                             ts_yyzzz_xxxzzz, \
                             ts_yyzzz_xxyyyy, \
                             ts_yyzzz_xxyyyz, \
                             ts_yyzzz_xxyyzz, \
                             ts_yyzzz_xxyzzz, \
                             ts_yyzzz_xxzzzz, \
                             ts_yyzzz_xyyyyy, \
                             ts_yyzzz_xyyyyz, \
                             ts_yyzzz_xyyyzz, \
                             ts_yyzzz_xyyzzz, \
                             ts_yyzzz_xyzzzz, \
                             ts_yyzzz_xzzzzz, \
                             ts_yyzzz_yyyyyy, \
                             ts_yyzzz_yyyyyz, \
                             ts_yyzzz_yyyyzz, \
                             ts_yyzzz_yyyzzz, \
                             ts_yyzzz_yyzzzz, \
                             ts_yyzzz_yzzzzz, \
                             ts_yyzzz_zzzzzz, \
                             ts_zzz_xxxxxx,   \
                             ts_zzz_xxxxxz,   \
                             ts_zzz_xxxxyz,   \
                             ts_zzz_xxxxzz,   \
                             ts_zzz_xxxyyz,   \
                             ts_zzz_xxxyzz,   \
                             ts_zzz_xxxzzz,   \
                             ts_zzz_xxyyyz,   \
                             ts_zzz_xxyyzz,   \
                             ts_zzz_xxyzzz,   \
                             ts_zzz_xxzzzz,   \
                             ts_zzz_xyyyyz,   \
                             ts_zzz_xyyyzz,   \
                             ts_zzz_xyyzzz,   \
                             ts_zzz_xyzzzz,   \
                             ts_zzz_xzzzzz,   \
                             ts_zzz_yyyyyz,   \
                             ts_zzz_yyyyzz,   \
                             ts_zzz_yyyzzz,   \
                             ts_zzz_yyzzzz,   \
                             ts_zzz_yzzzzz,   \
                             ts_zzz_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzzz_xxxxxx[i] =
            -2.0 * ts_zzz_xxxxxx[i] * fbe_0 * fz_0 + tk_zzz_xxxxxx[i] * fe_0 + tk_yzzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_yyzzz_xxxxxx[i] * fz_0;

        tk_yyzzz_xxxxxy[i] =
            -4.0 * ts_yyz_xxxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xxxxxy[i] * fe_0 + tk_yyzz_xxxxxy[i] * pa_z[i] + 2.0 * ts_yyzzz_xxxxxy[i] * fz_0;

        tk_yyzzz_xxxxxz[i] =
            -2.0 * ts_zzz_xxxxxz[i] * fbe_0 * fz_0 + tk_zzz_xxxxxz[i] * fe_0 + tk_yzzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxxxxz[i] * fz_0;

        tk_yyzzz_xxxxyy[i] =
            -4.0 * ts_yyz_xxxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xxxxyy[i] * fe_0 + tk_yyzz_xxxxyy[i] * pa_z[i] + 2.0 * ts_yyzzz_xxxxyy[i] * fz_0;

        tk_yyzzz_xxxxyz[i] = -2.0 * ts_zzz_xxxxyz[i] * fbe_0 * fz_0 + tk_zzz_xxxxyz[i] * fe_0 + tk_yzzz_xxxxz[i] * fe_0 +
                             tk_yzzz_xxxxyz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxxxyz[i] * fz_0;

        tk_yyzzz_xxxxzz[i] =
            -2.0 * ts_zzz_xxxxzz[i] * fbe_0 * fz_0 + tk_zzz_xxxxzz[i] * fe_0 + tk_yzzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxxxzz[i] * fz_0;

        tk_yyzzz_xxxyyy[i] =
            -4.0 * ts_yyz_xxxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xxxyyy[i] * fe_0 + tk_yyzz_xxxyyy[i] * pa_z[i] + 2.0 * ts_yyzzz_xxxyyy[i] * fz_0;

        tk_yyzzz_xxxyyz[i] = -2.0 * ts_zzz_xxxyyz[i] * fbe_0 * fz_0 + tk_zzz_xxxyyz[i] * fe_0 + 2.0 * tk_yzzz_xxxyz[i] * fe_0 +
                             tk_yzzz_xxxyyz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxxyyz[i] * fz_0;

        tk_yyzzz_xxxyzz[i] = -2.0 * ts_zzz_xxxyzz[i] * fbe_0 * fz_0 + tk_zzz_xxxyzz[i] * fe_0 + tk_yzzz_xxxzz[i] * fe_0 +
                             tk_yzzz_xxxyzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxxyzz[i] * fz_0;

        tk_yyzzz_xxxzzz[i] =
            -2.0 * ts_zzz_xxxzzz[i] * fbe_0 * fz_0 + tk_zzz_xxxzzz[i] * fe_0 + tk_yzzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxxzzz[i] * fz_0;

        tk_yyzzz_xxyyyy[i] =
            -4.0 * ts_yyz_xxyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xxyyyy[i] * fe_0 + tk_yyzz_xxyyyy[i] * pa_z[i] + 2.0 * ts_yyzzz_xxyyyy[i] * fz_0;

        tk_yyzzz_xxyyyz[i] = -2.0 * ts_zzz_xxyyyz[i] * fbe_0 * fz_0 + tk_zzz_xxyyyz[i] * fe_0 + 3.0 * tk_yzzz_xxyyz[i] * fe_0 +
                             tk_yzzz_xxyyyz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxyyyz[i] * fz_0;

        tk_yyzzz_xxyyzz[i] = -2.0 * ts_zzz_xxyyzz[i] * fbe_0 * fz_0 + tk_zzz_xxyyzz[i] * fe_0 + 2.0 * tk_yzzz_xxyzz[i] * fe_0 +
                             tk_yzzz_xxyyzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxyyzz[i] * fz_0;

        tk_yyzzz_xxyzzz[i] = -2.0 * ts_zzz_xxyzzz[i] * fbe_0 * fz_0 + tk_zzz_xxyzzz[i] * fe_0 + tk_yzzz_xxzzz[i] * fe_0 +
                             tk_yzzz_xxyzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxyzzz[i] * fz_0;

        tk_yyzzz_xxzzzz[i] =
            -2.0 * ts_zzz_xxzzzz[i] * fbe_0 * fz_0 + tk_zzz_xxzzzz[i] * fe_0 + tk_yzzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xxzzzz[i] * fz_0;

        tk_yyzzz_xyyyyy[i] =
            -4.0 * ts_yyz_xyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xyyyyy[i] * fe_0 + tk_yyzz_xyyyyy[i] * pa_z[i] + 2.0 * ts_yyzzz_xyyyyy[i] * fz_0;

        tk_yyzzz_xyyyyz[i] = -2.0 * ts_zzz_xyyyyz[i] * fbe_0 * fz_0 + tk_zzz_xyyyyz[i] * fe_0 + 4.0 * tk_yzzz_xyyyz[i] * fe_0 +
                             tk_yzzz_xyyyyz[i] * pa_y[i] + 2.0 * ts_yyzzz_xyyyyz[i] * fz_0;

        tk_yyzzz_xyyyzz[i] = -2.0 * ts_zzz_xyyyzz[i] * fbe_0 * fz_0 + tk_zzz_xyyyzz[i] * fe_0 + 3.0 * tk_yzzz_xyyzz[i] * fe_0 +
                             tk_yzzz_xyyyzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xyyyzz[i] * fz_0;

        tk_yyzzz_xyyzzz[i] = -2.0 * ts_zzz_xyyzzz[i] * fbe_0 * fz_0 + tk_zzz_xyyzzz[i] * fe_0 + 2.0 * tk_yzzz_xyzzz[i] * fe_0 +
                             tk_yzzz_xyyzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xyyzzz[i] * fz_0;

        tk_yyzzz_xyzzzz[i] = -2.0 * ts_zzz_xyzzzz[i] * fbe_0 * fz_0 + tk_zzz_xyzzzz[i] * fe_0 + tk_yzzz_xzzzz[i] * fe_0 +
                             tk_yzzz_xyzzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xyzzzz[i] * fz_0;

        tk_yyzzz_xzzzzz[i] =
            -2.0 * ts_zzz_xzzzzz[i] * fbe_0 * fz_0 + tk_zzz_xzzzzz[i] * fe_0 + tk_yzzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_xzzzzz[i] * fz_0;

        tk_yyzzz_yyyyyy[i] =
            -4.0 * ts_yyz_yyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_yyyyyy[i] * fe_0 + tk_yyzz_yyyyyy[i] * pa_z[i] + 2.0 * ts_yyzzz_yyyyyy[i] * fz_0;

        tk_yyzzz_yyyyyz[i] = -2.0 * ts_zzz_yyyyyz[i] * fbe_0 * fz_0 + tk_zzz_yyyyyz[i] * fe_0 + 5.0 * tk_yzzz_yyyyz[i] * fe_0 +
                             tk_yzzz_yyyyyz[i] * pa_y[i] + 2.0 * ts_yyzzz_yyyyyz[i] * fz_0;

        tk_yyzzz_yyyyzz[i] = -2.0 * ts_zzz_yyyyzz[i] * fbe_0 * fz_0 + tk_zzz_yyyyzz[i] * fe_0 + 4.0 * tk_yzzz_yyyzz[i] * fe_0 +
                             tk_yzzz_yyyyzz[i] * pa_y[i] + 2.0 * ts_yyzzz_yyyyzz[i] * fz_0;

        tk_yyzzz_yyyzzz[i] = -2.0 * ts_zzz_yyyzzz[i] * fbe_0 * fz_0 + tk_zzz_yyyzzz[i] * fe_0 + 3.0 * tk_yzzz_yyzzz[i] * fe_0 +
                             tk_yzzz_yyyzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_yyyzzz[i] * fz_0;

        tk_yyzzz_yyzzzz[i] = -2.0 * ts_zzz_yyzzzz[i] * fbe_0 * fz_0 + tk_zzz_yyzzzz[i] * fe_0 + 2.0 * tk_yzzz_yzzzz[i] * fe_0 +
                             tk_yzzz_yyzzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_yyzzzz[i] * fz_0;

        tk_yyzzz_yzzzzz[i] = -2.0 * ts_zzz_yzzzzz[i] * fbe_0 * fz_0 + tk_zzz_yzzzzz[i] * fe_0 + tk_yzzz_zzzzz[i] * fe_0 +
                             tk_yzzz_yzzzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_yzzzzz[i] * fz_0;

        tk_yyzzz_zzzzzz[i] =
            -2.0 * ts_zzz_zzzzzz[i] * fbe_0 * fz_0 + tk_zzz_zzzzzz[i] * fe_0 + tk_yzzz_zzzzzz[i] * pa_y[i] + 2.0 * ts_yyzzz_zzzzzz[i] * fz_0;
    }

    // Set up 532-560 components of targeted buffer : HI

    auto tk_yzzzz_xxxxxx = pbuffer.data(idx_kin_hi + 532);

    auto tk_yzzzz_xxxxxy = pbuffer.data(idx_kin_hi + 533);

    auto tk_yzzzz_xxxxxz = pbuffer.data(idx_kin_hi + 534);

    auto tk_yzzzz_xxxxyy = pbuffer.data(idx_kin_hi + 535);

    auto tk_yzzzz_xxxxyz = pbuffer.data(idx_kin_hi + 536);

    auto tk_yzzzz_xxxxzz = pbuffer.data(idx_kin_hi + 537);

    auto tk_yzzzz_xxxyyy = pbuffer.data(idx_kin_hi + 538);

    auto tk_yzzzz_xxxyyz = pbuffer.data(idx_kin_hi + 539);

    auto tk_yzzzz_xxxyzz = pbuffer.data(idx_kin_hi + 540);

    auto tk_yzzzz_xxxzzz = pbuffer.data(idx_kin_hi + 541);

    auto tk_yzzzz_xxyyyy = pbuffer.data(idx_kin_hi + 542);

    auto tk_yzzzz_xxyyyz = pbuffer.data(idx_kin_hi + 543);

    auto tk_yzzzz_xxyyzz = pbuffer.data(idx_kin_hi + 544);

    auto tk_yzzzz_xxyzzz = pbuffer.data(idx_kin_hi + 545);

    auto tk_yzzzz_xxzzzz = pbuffer.data(idx_kin_hi + 546);

    auto tk_yzzzz_xyyyyy = pbuffer.data(idx_kin_hi + 547);

    auto tk_yzzzz_xyyyyz = pbuffer.data(idx_kin_hi + 548);

    auto tk_yzzzz_xyyyzz = pbuffer.data(idx_kin_hi + 549);

    auto tk_yzzzz_xyyzzz = pbuffer.data(idx_kin_hi + 550);

    auto tk_yzzzz_xyzzzz = pbuffer.data(idx_kin_hi + 551);

    auto tk_yzzzz_xzzzzz = pbuffer.data(idx_kin_hi + 552);

    auto tk_yzzzz_yyyyyy = pbuffer.data(idx_kin_hi + 553);

    auto tk_yzzzz_yyyyyz = pbuffer.data(idx_kin_hi + 554);

    auto tk_yzzzz_yyyyzz = pbuffer.data(idx_kin_hi + 555);

    auto tk_yzzzz_yyyzzz = pbuffer.data(idx_kin_hi + 556);

    auto tk_yzzzz_yyzzzz = pbuffer.data(idx_kin_hi + 557);

    auto tk_yzzzz_yzzzzz = pbuffer.data(idx_kin_hi + 558);

    auto tk_yzzzz_zzzzzz = pbuffer.data(idx_kin_hi + 559);

#pragma omp simd aligned(pa_y,                \
                             tk_yzzzz_xxxxxx, \
                             tk_yzzzz_xxxxxy, \
                             tk_yzzzz_xxxxxz, \
                             tk_yzzzz_xxxxyy, \
                             tk_yzzzz_xxxxyz, \
                             tk_yzzzz_xxxxzz, \
                             tk_yzzzz_xxxyyy, \
                             tk_yzzzz_xxxyyz, \
                             tk_yzzzz_xxxyzz, \
                             tk_yzzzz_xxxzzz, \
                             tk_yzzzz_xxyyyy, \
                             tk_yzzzz_xxyyyz, \
                             tk_yzzzz_xxyyzz, \
                             tk_yzzzz_xxyzzz, \
                             tk_yzzzz_xxzzzz, \
                             tk_yzzzz_xyyyyy, \
                             tk_yzzzz_xyyyyz, \
                             tk_yzzzz_xyyyzz, \
                             tk_yzzzz_xyyzzz, \
                             tk_yzzzz_xyzzzz, \
                             tk_yzzzz_xzzzzz, \
                             tk_yzzzz_yyyyyy, \
                             tk_yzzzz_yyyyyz, \
                             tk_yzzzz_yyyyzz, \
                             tk_yzzzz_yyyzzz, \
                             tk_yzzzz_yyzzzz, \
                             tk_yzzzz_yzzzzz, \
                             tk_yzzzz_zzzzzz, \
                             tk_zzzz_xxxxx,   \
                             tk_zzzz_xxxxxx,  \
                             tk_zzzz_xxxxxy,  \
                             tk_zzzz_xxxxxz,  \
                             tk_zzzz_xxxxy,   \
                             tk_zzzz_xxxxyy,  \
                             tk_zzzz_xxxxyz,  \
                             tk_zzzz_xxxxz,   \
                             tk_zzzz_xxxxzz,  \
                             tk_zzzz_xxxyy,   \
                             tk_zzzz_xxxyyy,  \
                             tk_zzzz_xxxyyz,  \
                             tk_zzzz_xxxyz,   \
                             tk_zzzz_xxxyzz,  \
                             tk_zzzz_xxxzz,   \
                             tk_zzzz_xxxzzz,  \
                             tk_zzzz_xxyyy,   \
                             tk_zzzz_xxyyyy,  \
                             tk_zzzz_xxyyyz,  \
                             tk_zzzz_xxyyz,   \
                             tk_zzzz_xxyyzz,  \
                             tk_zzzz_xxyzz,   \
                             tk_zzzz_xxyzzz,  \
                             tk_zzzz_xxzzz,   \
                             tk_zzzz_xxzzzz,  \
                             tk_zzzz_xyyyy,   \
                             tk_zzzz_xyyyyy,  \
                             tk_zzzz_xyyyyz,  \
                             tk_zzzz_xyyyz,   \
                             tk_zzzz_xyyyzz,  \
                             tk_zzzz_xyyzz,   \
                             tk_zzzz_xyyzzz,  \
                             tk_zzzz_xyzzz,   \
                             tk_zzzz_xyzzzz,  \
                             tk_zzzz_xzzzz,   \
                             tk_zzzz_xzzzzz,  \
                             tk_zzzz_yyyyy,   \
                             tk_zzzz_yyyyyy,  \
                             tk_zzzz_yyyyyz,  \
                             tk_zzzz_yyyyz,   \
                             tk_zzzz_yyyyzz,  \
                             tk_zzzz_yyyzz,   \
                             tk_zzzz_yyyzzz,  \
                             tk_zzzz_yyzzz,   \
                             tk_zzzz_yyzzzz,  \
                             tk_zzzz_yzzzz,   \
                             tk_zzzz_yzzzzz,  \
                             tk_zzzz_zzzzz,   \
                             tk_zzzz_zzzzzz,  \
                             ts_yzzzz_xxxxxx, \
                             ts_yzzzz_xxxxxy, \
                             ts_yzzzz_xxxxxz, \
                             ts_yzzzz_xxxxyy, \
                             ts_yzzzz_xxxxyz, \
                             ts_yzzzz_xxxxzz, \
                             ts_yzzzz_xxxyyy, \
                             ts_yzzzz_xxxyyz, \
                             ts_yzzzz_xxxyzz, \
                             ts_yzzzz_xxxzzz, \
                             ts_yzzzz_xxyyyy, \
                             ts_yzzzz_xxyyyz, \
                             ts_yzzzz_xxyyzz, \
                             ts_yzzzz_xxyzzz, \
                             ts_yzzzz_xxzzzz, \
                             ts_yzzzz_xyyyyy, \
                             ts_yzzzz_xyyyyz, \
                             ts_yzzzz_xyyyzz, \
                             ts_yzzzz_xyyzzz, \
                             ts_yzzzz_xyzzzz, \
                             ts_yzzzz_xzzzzz, \
                             ts_yzzzz_yyyyyy, \
                             ts_yzzzz_yyyyyz, \
                             ts_yzzzz_yyyyzz, \
                             ts_yzzzz_yyyzzz, \
                             ts_yzzzz_yyzzzz, \
                             ts_yzzzz_yzzzzz, \
                             ts_yzzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzzz_xxxxxx[i] = tk_zzzz_xxxxxx[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxxxx[i] * fz_0;

        tk_yzzzz_xxxxxy[i] = tk_zzzz_xxxxx[i] * fe_0 + tk_zzzz_xxxxxy[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxxxy[i] * fz_0;

        tk_yzzzz_xxxxxz[i] = tk_zzzz_xxxxxz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxxxz[i] * fz_0;

        tk_yzzzz_xxxxyy[i] = 2.0 * tk_zzzz_xxxxy[i] * fe_0 + tk_zzzz_xxxxyy[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxxyy[i] * fz_0;

        tk_yzzzz_xxxxyz[i] = tk_zzzz_xxxxz[i] * fe_0 + tk_zzzz_xxxxyz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxxyz[i] * fz_0;

        tk_yzzzz_xxxxzz[i] = tk_zzzz_xxxxzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxxzz[i] * fz_0;

        tk_yzzzz_xxxyyy[i] = 3.0 * tk_zzzz_xxxyy[i] * fe_0 + tk_zzzz_xxxyyy[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxyyy[i] * fz_0;

        tk_yzzzz_xxxyyz[i] = 2.0 * tk_zzzz_xxxyz[i] * fe_0 + tk_zzzz_xxxyyz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxyyz[i] * fz_0;

        tk_yzzzz_xxxyzz[i] = tk_zzzz_xxxzz[i] * fe_0 + tk_zzzz_xxxyzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxyzz[i] * fz_0;

        tk_yzzzz_xxxzzz[i] = tk_zzzz_xxxzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxxzzz[i] * fz_0;

        tk_yzzzz_xxyyyy[i] = 4.0 * tk_zzzz_xxyyy[i] * fe_0 + tk_zzzz_xxyyyy[i] * pa_y[i] + 2.0 * ts_yzzzz_xxyyyy[i] * fz_0;

        tk_yzzzz_xxyyyz[i] = 3.0 * tk_zzzz_xxyyz[i] * fe_0 + tk_zzzz_xxyyyz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxyyyz[i] * fz_0;

        tk_yzzzz_xxyyzz[i] = 2.0 * tk_zzzz_xxyzz[i] * fe_0 + tk_zzzz_xxyyzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxyyzz[i] * fz_0;

        tk_yzzzz_xxyzzz[i] = tk_zzzz_xxzzz[i] * fe_0 + tk_zzzz_xxyzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxyzzz[i] * fz_0;

        tk_yzzzz_xxzzzz[i] = tk_zzzz_xxzzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xxzzzz[i] * fz_0;

        tk_yzzzz_xyyyyy[i] = 5.0 * tk_zzzz_xyyyy[i] * fe_0 + tk_zzzz_xyyyyy[i] * pa_y[i] + 2.0 * ts_yzzzz_xyyyyy[i] * fz_0;

        tk_yzzzz_xyyyyz[i] = 4.0 * tk_zzzz_xyyyz[i] * fe_0 + tk_zzzz_xyyyyz[i] * pa_y[i] + 2.0 * ts_yzzzz_xyyyyz[i] * fz_0;

        tk_yzzzz_xyyyzz[i] = 3.0 * tk_zzzz_xyyzz[i] * fe_0 + tk_zzzz_xyyyzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xyyyzz[i] * fz_0;

        tk_yzzzz_xyyzzz[i] = 2.0 * tk_zzzz_xyzzz[i] * fe_0 + tk_zzzz_xyyzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xyyzzz[i] * fz_0;

        tk_yzzzz_xyzzzz[i] = tk_zzzz_xzzzz[i] * fe_0 + tk_zzzz_xyzzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xyzzzz[i] * fz_0;

        tk_yzzzz_xzzzzz[i] = tk_zzzz_xzzzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_xzzzzz[i] * fz_0;

        tk_yzzzz_yyyyyy[i] = 6.0 * tk_zzzz_yyyyy[i] * fe_0 + tk_zzzz_yyyyyy[i] * pa_y[i] + 2.0 * ts_yzzzz_yyyyyy[i] * fz_0;

        tk_yzzzz_yyyyyz[i] = 5.0 * tk_zzzz_yyyyz[i] * fe_0 + tk_zzzz_yyyyyz[i] * pa_y[i] + 2.0 * ts_yzzzz_yyyyyz[i] * fz_0;

        tk_yzzzz_yyyyzz[i] = 4.0 * tk_zzzz_yyyzz[i] * fe_0 + tk_zzzz_yyyyzz[i] * pa_y[i] + 2.0 * ts_yzzzz_yyyyzz[i] * fz_0;

        tk_yzzzz_yyyzzz[i] = 3.0 * tk_zzzz_yyzzz[i] * fe_0 + tk_zzzz_yyyzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_yyyzzz[i] * fz_0;

        tk_yzzzz_yyzzzz[i] = 2.0 * tk_zzzz_yzzzz[i] * fe_0 + tk_zzzz_yyzzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_yyzzzz[i] * fz_0;

        tk_yzzzz_yzzzzz[i] = tk_zzzz_zzzzz[i] * fe_0 + tk_zzzz_yzzzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_yzzzzz[i] * fz_0;

        tk_yzzzz_zzzzzz[i] = tk_zzzz_zzzzzz[i] * pa_y[i] + 2.0 * ts_yzzzz_zzzzzz[i] * fz_0;
    }

    // Set up 560-588 components of targeted buffer : HI

    auto tk_zzzzz_xxxxxx = pbuffer.data(idx_kin_hi + 560);

    auto tk_zzzzz_xxxxxy = pbuffer.data(idx_kin_hi + 561);

    auto tk_zzzzz_xxxxxz = pbuffer.data(idx_kin_hi + 562);

    auto tk_zzzzz_xxxxyy = pbuffer.data(idx_kin_hi + 563);

    auto tk_zzzzz_xxxxyz = pbuffer.data(idx_kin_hi + 564);

    auto tk_zzzzz_xxxxzz = pbuffer.data(idx_kin_hi + 565);

    auto tk_zzzzz_xxxyyy = pbuffer.data(idx_kin_hi + 566);

    auto tk_zzzzz_xxxyyz = pbuffer.data(idx_kin_hi + 567);

    auto tk_zzzzz_xxxyzz = pbuffer.data(idx_kin_hi + 568);

    auto tk_zzzzz_xxxzzz = pbuffer.data(idx_kin_hi + 569);

    auto tk_zzzzz_xxyyyy = pbuffer.data(idx_kin_hi + 570);

    auto tk_zzzzz_xxyyyz = pbuffer.data(idx_kin_hi + 571);

    auto tk_zzzzz_xxyyzz = pbuffer.data(idx_kin_hi + 572);

    auto tk_zzzzz_xxyzzz = pbuffer.data(idx_kin_hi + 573);

    auto tk_zzzzz_xxzzzz = pbuffer.data(idx_kin_hi + 574);

    auto tk_zzzzz_xyyyyy = pbuffer.data(idx_kin_hi + 575);

    auto tk_zzzzz_xyyyyz = pbuffer.data(idx_kin_hi + 576);

    auto tk_zzzzz_xyyyzz = pbuffer.data(idx_kin_hi + 577);

    auto tk_zzzzz_xyyzzz = pbuffer.data(idx_kin_hi + 578);

    auto tk_zzzzz_xyzzzz = pbuffer.data(idx_kin_hi + 579);

    auto tk_zzzzz_xzzzzz = pbuffer.data(idx_kin_hi + 580);

    auto tk_zzzzz_yyyyyy = pbuffer.data(idx_kin_hi + 581);

    auto tk_zzzzz_yyyyyz = pbuffer.data(idx_kin_hi + 582);

    auto tk_zzzzz_yyyyzz = pbuffer.data(idx_kin_hi + 583);

    auto tk_zzzzz_yyyzzz = pbuffer.data(idx_kin_hi + 584);

    auto tk_zzzzz_yyzzzz = pbuffer.data(idx_kin_hi + 585);

    auto tk_zzzzz_yzzzzz = pbuffer.data(idx_kin_hi + 586);

    auto tk_zzzzz_zzzzzz = pbuffer.data(idx_kin_hi + 587);

#pragma omp simd aligned(pa_z,                \
                             tk_zzz_xxxxxx,   \
                             tk_zzz_xxxxxy,   \
                             tk_zzz_xxxxxz,   \
                             tk_zzz_xxxxyy,   \
                             tk_zzz_xxxxyz,   \
                             tk_zzz_xxxxzz,   \
                             tk_zzz_xxxyyy,   \
                             tk_zzz_xxxyyz,   \
                             tk_zzz_xxxyzz,   \
                             tk_zzz_xxxzzz,   \
                             tk_zzz_xxyyyy,   \
                             tk_zzz_xxyyyz,   \
                             tk_zzz_xxyyzz,   \
                             tk_zzz_xxyzzz,   \
                             tk_zzz_xxzzzz,   \
                             tk_zzz_xyyyyy,   \
                             tk_zzz_xyyyyz,   \
                             tk_zzz_xyyyzz,   \
                             tk_zzz_xyyzzz,   \
                             tk_zzz_xyzzzz,   \
                             tk_zzz_xzzzzz,   \
                             tk_zzz_yyyyyy,   \
                             tk_zzz_yyyyyz,   \
                             tk_zzz_yyyyzz,   \
                             tk_zzz_yyyzzz,   \
                             tk_zzz_yyzzzz,   \
                             tk_zzz_yzzzzz,   \
                             tk_zzz_zzzzzz,   \
                             tk_zzzz_xxxxx,   \
                             tk_zzzz_xxxxxx,  \
                             tk_zzzz_xxxxxy,  \
                             tk_zzzz_xxxxxz,  \
                             tk_zzzz_xxxxy,   \
                             tk_zzzz_xxxxyy,  \
                             tk_zzzz_xxxxyz,  \
                             tk_zzzz_xxxxz,   \
                             tk_zzzz_xxxxzz,  \
                             tk_zzzz_xxxyy,   \
                             tk_zzzz_xxxyyy,  \
                             tk_zzzz_xxxyyz,  \
                             tk_zzzz_xxxyz,   \
                             tk_zzzz_xxxyzz,  \
                             tk_zzzz_xxxzz,   \
                             tk_zzzz_xxxzzz,  \
                             tk_zzzz_xxyyy,   \
                             tk_zzzz_xxyyyy,  \
                             tk_zzzz_xxyyyz,  \
                             tk_zzzz_xxyyz,   \
                             tk_zzzz_xxyyzz,  \
                             tk_zzzz_xxyzz,   \
                             tk_zzzz_xxyzzz,  \
                             tk_zzzz_xxzzz,   \
                             tk_zzzz_xxzzzz,  \
                             tk_zzzz_xyyyy,   \
                             tk_zzzz_xyyyyy,  \
                             tk_zzzz_xyyyyz,  \
                             tk_zzzz_xyyyz,   \
                             tk_zzzz_xyyyzz,  \
                             tk_zzzz_xyyzz,   \
                             tk_zzzz_xyyzzz,  \
                             tk_zzzz_xyzzz,   \
                             tk_zzzz_xyzzzz,  \
                             tk_zzzz_xzzzz,   \
                             tk_zzzz_xzzzzz,  \
                             tk_zzzz_yyyyy,   \
                             tk_zzzz_yyyyyy,  \
                             tk_zzzz_yyyyyz,  \
                             tk_zzzz_yyyyz,   \
                             tk_zzzz_yyyyzz,  \
                             tk_zzzz_yyyzz,   \
                             tk_zzzz_yyyzzz,  \
                             tk_zzzz_yyzzz,   \
                             tk_zzzz_yyzzzz,  \
                             tk_zzzz_yzzzz,   \
                             tk_zzzz_yzzzzz,  \
                             tk_zzzz_zzzzz,   \
                             tk_zzzz_zzzzzz,  \
                             tk_zzzzz_xxxxxx, \
                             tk_zzzzz_xxxxxy, \
                             tk_zzzzz_xxxxxz, \
                             tk_zzzzz_xxxxyy, \
                             tk_zzzzz_xxxxyz, \
                             tk_zzzzz_xxxxzz, \
                             tk_zzzzz_xxxyyy, \
                             tk_zzzzz_xxxyyz, \
                             tk_zzzzz_xxxyzz, \
                             tk_zzzzz_xxxzzz, \
                             tk_zzzzz_xxyyyy, \
                             tk_zzzzz_xxyyyz, \
                             tk_zzzzz_xxyyzz, \
                             tk_zzzzz_xxyzzz, \
                             tk_zzzzz_xxzzzz, \
                             tk_zzzzz_xyyyyy, \
                             tk_zzzzz_xyyyyz, \
                             tk_zzzzz_xyyyzz, \
                             tk_zzzzz_xyyzzz, \
                             tk_zzzzz_xyzzzz, \
                             tk_zzzzz_xzzzzz, \
                             tk_zzzzz_yyyyyy, \
                             tk_zzzzz_yyyyyz, \
                             tk_zzzzz_yyyyzz, \
                             tk_zzzzz_yyyzzz, \
                             tk_zzzzz_yyzzzz, \
                             tk_zzzzz_yzzzzz, \
                             tk_zzzzz_zzzzzz, \
                             ts_zzz_xxxxxx,   \
                             ts_zzz_xxxxxy,   \
                             ts_zzz_xxxxxz,   \
                             ts_zzz_xxxxyy,   \
                             ts_zzz_xxxxyz,   \
                             ts_zzz_xxxxzz,   \
                             ts_zzz_xxxyyy,   \
                             ts_zzz_xxxyyz,   \
                             ts_zzz_xxxyzz,   \
                             ts_zzz_xxxzzz,   \
                             ts_zzz_xxyyyy,   \
                             ts_zzz_xxyyyz,   \
                             ts_zzz_xxyyzz,   \
                             ts_zzz_xxyzzz,   \
                             ts_zzz_xxzzzz,   \
                             ts_zzz_xyyyyy,   \
                             ts_zzz_xyyyyz,   \
                             ts_zzz_xyyyzz,   \
                             ts_zzz_xyyzzz,   \
                             ts_zzz_xyzzzz,   \
                             ts_zzz_xzzzzz,   \
                             ts_zzz_yyyyyy,   \
                             ts_zzz_yyyyyz,   \
                             ts_zzz_yyyyzz,   \
                             ts_zzz_yyyzzz,   \
                             ts_zzz_yyzzzz,   \
                             ts_zzz_yzzzzz,   \
                             ts_zzz_zzzzzz,   \
                             ts_zzzzz_xxxxxx, \
                             ts_zzzzz_xxxxxy, \
                             ts_zzzzz_xxxxxz, \
                             ts_zzzzz_xxxxyy, \
                             ts_zzzzz_xxxxyz, \
                             ts_zzzzz_xxxxzz, \
                             ts_zzzzz_xxxyyy, \
                             ts_zzzzz_xxxyyz, \
                             ts_zzzzz_xxxyzz, \
                             ts_zzzzz_xxxzzz, \
                             ts_zzzzz_xxyyyy, \
                             ts_zzzzz_xxyyyz, \
                             ts_zzzzz_xxyyzz, \
                             ts_zzzzz_xxyzzz, \
                             ts_zzzzz_xxzzzz, \
                             ts_zzzzz_xyyyyy, \
                             ts_zzzzz_xyyyyz, \
                             ts_zzzzz_xyyyzz, \
                             ts_zzzzz_xyyzzz, \
                             ts_zzzzz_xyzzzz, \
                             ts_zzzzz_xzzzzz, \
                             ts_zzzzz_yyyyyy, \
                             ts_zzzzz_yyyyyz, \
                             ts_zzzzz_yyyyzz, \
                             ts_zzzzz_yyyzzz, \
                             ts_zzzzz_yyzzzz, \
                             ts_zzzzz_yzzzzz, \
                             ts_zzzzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzzz_xxxxxx[i] =
            -8.0 * ts_zzz_xxxxxx[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxxxx[i] * fe_0 + tk_zzzz_xxxxxx[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxxxx[i] * fz_0;

        tk_zzzzz_xxxxxy[i] =
            -8.0 * ts_zzz_xxxxxy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxxxy[i] * fe_0 + tk_zzzz_xxxxxy[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxxxy[i] * fz_0;

        tk_zzzzz_xxxxxz[i] = -8.0 * ts_zzz_xxxxxz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxxxz[i] * fe_0 + tk_zzzz_xxxxx[i] * fe_0 +
                             tk_zzzz_xxxxxz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxxxz[i] * fz_0;

        tk_zzzzz_xxxxyy[i] =
            -8.0 * ts_zzz_xxxxyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxxyy[i] * fe_0 + tk_zzzz_xxxxyy[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxxyy[i] * fz_0;

        tk_zzzzz_xxxxyz[i] = -8.0 * ts_zzz_xxxxyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxxyz[i] * fe_0 + tk_zzzz_xxxxy[i] * fe_0 +
                             tk_zzzz_xxxxyz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxxyz[i] * fz_0;

        tk_zzzzz_xxxxzz[i] = -8.0 * ts_zzz_xxxxzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxxzz[i] * fe_0 + 2.0 * tk_zzzz_xxxxz[i] * fe_0 +
                             tk_zzzz_xxxxzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxxzz[i] * fz_0;

        tk_zzzzz_xxxyyy[i] =
            -8.0 * ts_zzz_xxxyyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxyyy[i] * fe_0 + tk_zzzz_xxxyyy[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxyyy[i] * fz_0;

        tk_zzzzz_xxxyyz[i] = -8.0 * ts_zzz_xxxyyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxyyz[i] * fe_0 + tk_zzzz_xxxyy[i] * fe_0 +
                             tk_zzzz_xxxyyz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxyyz[i] * fz_0;

        tk_zzzzz_xxxyzz[i] = -8.0 * ts_zzz_xxxyzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxyzz[i] * fe_0 + 2.0 * tk_zzzz_xxxyz[i] * fe_0 +
                             tk_zzzz_xxxyzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxyzz[i] * fz_0;

        tk_zzzzz_xxxzzz[i] = -8.0 * ts_zzz_xxxzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxxzzz[i] * fe_0 + 3.0 * tk_zzzz_xxxzz[i] * fe_0 +
                             tk_zzzz_xxxzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxxzzz[i] * fz_0;

        tk_zzzzz_xxyyyy[i] =
            -8.0 * ts_zzz_xxyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxyyyy[i] * fe_0 + tk_zzzz_xxyyyy[i] * pa_z[i] + 2.0 * ts_zzzzz_xxyyyy[i] * fz_0;

        tk_zzzzz_xxyyyz[i] = -8.0 * ts_zzz_xxyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxyyyz[i] * fe_0 + tk_zzzz_xxyyy[i] * fe_0 +
                             tk_zzzz_xxyyyz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxyyyz[i] * fz_0;

        tk_zzzzz_xxyyzz[i] = -8.0 * ts_zzz_xxyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxyyzz[i] * fe_0 + 2.0 * tk_zzzz_xxyyz[i] * fe_0 +
                             tk_zzzz_xxyyzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxyyzz[i] * fz_0;

        tk_zzzzz_xxyzzz[i] = -8.0 * ts_zzz_xxyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxyzzz[i] * fe_0 + 3.0 * tk_zzzz_xxyzz[i] * fe_0 +
                             tk_zzzz_xxyzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxyzzz[i] * fz_0;

        tk_zzzzz_xxzzzz[i] = -8.0 * ts_zzz_xxzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xxzzzz[i] * fe_0 + 4.0 * tk_zzzz_xxzzz[i] * fe_0 +
                             tk_zzzz_xxzzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xxzzzz[i] * fz_0;

        tk_zzzzz_xyyyyy[i] =
            -8.0 * ts_zzz_xyyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyyyyy[i] * fe_0 + tk_zzzz_xyyyyy[i] * pa_z[i] + 2.0 * ts_zzzzz_xyyyyy[i] * fz_0;

        tk_zzzzz_xyyyyz[i] = -8.0 * ts_zzz_xyyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyyyyz[i] * fe_0 + tk_zzzz_xyyyy[i] * fe_0 +
                             tk_zzzz_xyyyyz[i] * pa_z[i] + 2.0 * ts_zzzzz_xyyyyz[i] * fz_0;

        tk_zzzzz_xyyyzz[i] = -8.0 * ts_zzz_xyyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyyyzz[i] * fe_0 + 2.0 * tk_zzzz_xyyyz[i] * fe_0 +
                             tk_zzzz_xyyyzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xyyyzz[i] * fz_0;

        tk_zzzzz_xyyzzz[i] = -8.0 * ts_zzz_xyyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyyzzz[i] * fe_0 + 3.0 * tk_zzzz_xyyzz[i] * fe_0 +
                             tk_zzzz_xyyzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xyyzzz[i] * fz_0;

        tk_zzzzz_xyzzzz[i] = -8.0 * ts_zzz_xyzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xyzzzz[i] * fe_0 + 4.0 * tk_zzzz_xyzzz[i] * fe_0 +
                             tk_zzzz_xyzzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xyzzzz[i] * fz_0;

        tk_zzzzz_xzzzzz[i] = -8.0 * ts_zzz_xzzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xzzzzz[i] * fe_0 + 5.0 * tk_zzzz_xzzzz[i] * fe_0 +
                             tk_zzzz_xzzzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_xzzzzz[i] * fz_0;

        tk_zzzzz_yyyyyy[i] =
            -8.0 * ts_zzz_yyyyyy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyyyyy[i] * fe_0 + tk_zzzz_yyyyyy[i] * pa_z[i] + 2.0 * ts_zzzzz_yyyyyy[i] * fz_0;

        tk_zzzzz_yyyyyz[i] = -8.0 * ts_zzz_yyyyyz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyyyyz[i] * fe_0 + tk_zzzz_yyyyy[i] * fe_0 +
                             tk_zzzz_yyyyyz[i] * pa_z[i] + 2.0 * ts_zzzzz_yyyyyz[i] * fz_0;

        tk_zzzzz_yyyyzz[i] = -8.0 * ts_zzz_yyyyzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyyyzz[i] * fe_0 + 2.0 * tk_zzzz_yyyyz[i] * fe_0 +
                             tk_zzzz_yyyyzz[i] * pa_z[i] + 2.0 * ts_zzzzz_yyyyzz[i] * fz_0;

        tk_zzzzz_yyyzzz[i] = -8.0 * ts_zzz_yyyzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyyzzz[i] * fe_0 + 3.0 * tk_zzzz_yyyzz[i] * fe_0 +
                             tk_zzzz_yyyzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_yyyzzz[i] * fz_0;

        tk_zzzzz_yyzzzz[i] = -8.0 * ts_zzz_yyzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yyzzzz[i] * fe_0 + 4.0 * tk_zzzz_yyzzz[i] * fe_0 +
                             tk_zzzz_yyzzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_yyzzzz[i] * fz_0;

        tk_zzzzz_yzzzzz[i] = -8.0 * ts_zzz_yzzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yzzzzz[i] * fe_0 + 5.0 * tk_zzzz_yzzzz[i] * fe_0 +
                             tk_zzzz_yzzzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_yzzzzz[i] * fz_0;

        tk_zzzzz_zzzzzz[i] = -8.0 * ts_zzz_zzzzzz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_zzzzzz[i] * fe_0 + 6.0 * tk_zzzz_zzzzz[i] * fe_0 +
                             tk_zzzz_zzzzzz[i] * pa_z[i] + 2.0 * ts_zzzzz_zzzzzz[i] * fz_0;
    }
}

}  // namespace kinrec
