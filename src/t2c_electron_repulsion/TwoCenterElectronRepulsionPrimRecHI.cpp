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

#include "TwoCenterElectronRepulsionPrimRecHI.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_hi(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_hi,
                                const size_t idx_eri_0_fi,
                                const size_t idx_eri_1_fi,
                                const size_t idx_eri_1_gh,
                                const size_t idx_eri_1_gi,
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

    // Set up components of auxiliary buffer : FI

    auto g_xxx_xxxxxx_0 = pbuffer.data(idx_eri_0_fi);

    auto g_xxx_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 1);

    auto g_xxx_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 2);

    auto g_xxx_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 3);

    auto g_xxx_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 4);

    auto g_xxx_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 5);

    auto g_xxx_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 6);

    auto g_xxx_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 7);

    auto g_xxx_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 8);

    auto g_xxx_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 9);

    auto g_xxx_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 10);

    auto g_xxx_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 11);

    auto g_xxx_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 12);

    auto g_xxx_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 13);

    auto g_xxx_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 14);

    auto g_xxx_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 15);

    auto g_xxx_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 16);

    auto g_xxx_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 17);

    auto g_xxx_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 18);

    auto g_xxx_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 19);

    auto g_xxx_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 20);

    auto g_xxx_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 21);

    auto g_xxx_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 22);

    auto g_xxx_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 23);

    auto g_xxx_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 24);

    auto g_xxx_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 25);

    auto g_xxx_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 26);

    auto g_xxx_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 27);

    auto g_xxy_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 28);

    auto g_xxy_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 30);

    auto g_xxy_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 33);

    auto g_xxy_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 37);

    auto g_xxy_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 42);

    auto g_xxy_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 48);

    auto g_xxz_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 56);

    auto g_xxz_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 57);

    auto g_xxz_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 59);

    auto g_xxz_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 62);

    auto g_xxz_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 66);

    auto g_xxz_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 71);

    auto g_xyy_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 85);

    auto g_xyy_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 87);

    auto g_xyy_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 88);

    auto g_xyy_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 90);

    auto g_xyy_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 91);

    auto g_xyy_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 92);

    auto g_xyy_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 94);

    auto g_xyy_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 95);

    auto g_xyy_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 96);

    auto g_xyy_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 97);

    auto g_xyy_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 99);

    auto g_xyy_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 100);

    auto g_xyy_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 101);

    auto g_xyy_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 102);

    auto g_xyy_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 103);

    auto g_xyy_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 105);

    auto g_xyy_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 106);

    auto g_xyy_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 107);

    auto g_xyy_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 108);

    auto g_xyy_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 109);

    auto g_xyy_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 110);

    auto g_xyy_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 111);

    auto g_xzz_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 142);

    auto g_xzz_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 144);

    auto g_xzz_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 145);

    auto g_xzz_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 147);

    auto g_xzz_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 148);

    auto g_xzz_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 149);

    auto g_xzz_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 151);

    auto g_xzz_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 152);

    auto g_xzz_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 153);

    auto g_xzz_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 154);

    auto g_xzz_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 156);

    auto g_xzz_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 157);

    auto g_xzz_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 158);

    auto g_xzz_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 159);

    auto g_xzz_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 160);

    auto g_xzz_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 161);

    auto g_xzz_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 162);

    auto g_xzz_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 163);

    auto g_xzz_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 164);

    auto g_xzz_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 165);

    auto g_xzz_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 166);

    auto g_xzz_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 167);

    auto g_yyy_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 168);

    auto g_yyy_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 169);

    auto g_yyy_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 170);

    auto g_yyy_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 171);

    auto g_yyy_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 172);

    auto g_yyy_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 173);

    auto g_yyy_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 174);

    auto g_yyy_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 175);

    auto g_yyy_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 176);

    auto g_yyy_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 177);

    auto g_yyy_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 178);

    auto g_yyy_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 179);

    auto g_yyy_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 180);

    auto g_yyy_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 181);

    auto g_yyy_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 182);

    auto g_yyy_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 183);

    auto g_yyy_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 184);

    auto g_yyy_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 185);

    auto g_yyy_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 186);

    auto g_yyy_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 187);

    auto g_yyy_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 188);

    auto g_yyy_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 189);

    auto g_yyy_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 190);

    auto g_yyy_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 191);

    auto g_yyy_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 192);

    auto g_yyy_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 193);

    auto g_yyy_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 194);

    auto g_yyy_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 195);

    auto g_yyz_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 197);

    auto g_yyz_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 199);

    auto g_yyz_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 202);

    auto g_yyz_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 206);

    auto g_yyz_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 211);

    auto g_yyz_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 217);

    auto g_yzz_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 224);

    auto g_yzz_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 226);

    auto g_yzz_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 228);

    auto g_yzz_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 229);

    auto g_yzz_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 231);

    auto g_yzz_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 232);

    auto g_yzz_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 233);

    auto g_yzz_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 235);

    auto g_yzz_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 236);

    auto g_yzz_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 237);

    auto g_yzz_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 238);

    auto g_yzz_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 240);

    auto g_yzz_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 241);

    auto g_yzz_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 242);

    auto g_yzz_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 243);

    auto g_yzz_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 244);

    auto g_yzz_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 246);

    auto g_yzz_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 247);

    auto g_yzz_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 248);

    auto g_yzz_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 249);

    auto g_yzz_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 250);

    auto g_yzz_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 251);

    auto g_zzz_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 252);

    auto g_zzz_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 253);

    auto g_zzz_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 254);

    auto g_zzz_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 255);

    auto g_zzz_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 256);

    auto g_zzz_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 257);

    auto g_zzz_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 258);

    auto g_zzz_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 259);

    auto g_zzz_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 260);

    auto g_zzz_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 261);

    auto g_zzz_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 262);

    auto g_zzz_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 263);

    auto g_zzz_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 264);

    auto g_zzz_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 265);

    auto g_zzz_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 266);

    auto g_zzz_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 267);

    auto g_zzz_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 268);

    auto g_zzz_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 269);

    auto g_zzz_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 270);

    auto g_zzz_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 271);

    auto g_zzz_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 272);

    auto g_zzz_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 273);

    auto g_zzz_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 274);

    auto g_zzz_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 275);

    auto g_zzz_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 276);

    auto g_zzz_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 277);

    auto g_zzz_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 278);

    auto g_zzz_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 279);

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

    auto g_xxy_xxxxxz_1 = pbuffer.data(idx_eri_1_fi + 30);

    auto g_xxy_xxxxzz_1 = pbuffer.data(idx_eri_1_fi + 33);

    auto g_xxy_xxxzzz_1 = pbuffer.data(idx_eri_1_fi + 37);

    auto g_xxy_xxzzzz_1 = pbuffer.data(idx_eri_1_fi + 42);

    auto g_xxy_xzzzzz_1 = pbuffer.data(idx_eri_1_fi + 48);

    auto g_xxz_xxxxxx_1 = pbuffer.data(idx_eri_1_fi + 56);

    auto g_xxz_xxxxxy_1 = pbuffer.data(idx_eri_1_fi + 57);

    auto g_xxz_xxxxyy_1 = pbuffer.data(idx_eri_1_fi + 59);

    auto g_xxz_xxxyyy_1 = pbuffer.data(idx_eri_1_fi + 62);

    auto g_xxz_xxyyyy_1 = pbuffer.data(idx_eri_1_fi + 66);

    auto g_xxz_xyyyyy_1 = pbuffer.data(idx_eri_1_fi + 71);

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

    auto g_yyz_xxxxyy_1 = pbuffer.data(idx_eri_1_fi + 199);

    auto g_yyz_xxxyyy_1 = pbuffer.data(idx_eri_1_fi + 202);

    auto g_yyz_xxyyyy_1 = pbuffer.data(idx_eri_1_fi + 206);

    auto g_yyz_xyyyyy_1 = pbuffer.data(idx_eri_1_fi + 211);

    auto g_yyz_yyyyyy_1 = pbuffer.data(idx_eri_1_fi + 217);

    auto g_yzz_xxxxxx_1 = pbuffer.data(idx_eri_1_fi + 224);

    auto g_yzz_xxxxxz_1 = pbuffer.data(idx_eri_1_fi + 226);

    auto g_yzz_xxxxyz_1 = pbuffer.data(idx_eri_1_fi + 228);

    auto g_yzz_xxxxzz_1 = pbuffer.data(idx_eri_1_fi + 229);

    auto g_yzz_xxxyyz_1 = pbuffer.data(idx_eri_1_fi + 231);

    auto g_yzz_xxxyzz_1 = pbuffer.data(idx_eri_1_fi + 232);

    auto g_yzz_xxxzzz_1 = pbuffer.data(idx_eri_1_fi + 233);

    auto g_yzz_xxyyyz_1 = pbuffer.data(idx_eri_1_fi + 235);

    auto g_yzz_xxyyzz_1 = pbuffer.data(idx_eri_1_fi + 236);

    auto g_yzz_xxyzzz_1 = pbuffer.data(idx_eri_1_fi + 237);

    auto g_yzz_xxzzzz_1 = pbuffer.data(idx_eri_1_fi + 238);

    auto g_yzz_xyyyyz_1 = pbuffer.data(idx_eri_1_fi + 240);

    auto g_yzz_xyyyzz_1 = pbuffer.data(idx_eri_1_fi + 241);

    auto g_yzz_xyyzzz_1 = pbuffer.data(idx_eri_1_fi + 242);

    auto g_yzz_xyzzzz_1 = pbuffer.data(idx_eri_1_fi + 243);

    auto g_yzz_xzzzzz_1 = pbuffer.data(idx_eri_1_fi + 244);

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

    // Set up components of auxiliary buffer : GH

    auto g_xxxx_xxxxx_1 = pbuffer.data(idx_eri_1_gh);

    auto g_xxxx_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 1);

    auto g_xxxx_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 2);

    auto g_xxxx_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 3);

    auto g_xxxx_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 4);

    auto g_xxxx_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 5);

    auto g_xxxx_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 6);

    auto g_xxxx_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 7);

    auto g_xxxx_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 8);

    auto g_xxxx_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 9);

    auto g_xxxx_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 10);

    auto g_xxxx_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 11);

    auto g_xxxx_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 12);

    auto g_xxxx_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 13);

    auto g_xxxx_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 14);

    auto g_xxxx_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 15);

    auto g_xxxx_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 16);

    auto g_xxxx_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 17);

    auto g_xxxx_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 18);

    auto g_xxxx_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 19);

    auto g_xxxx_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 20);

    auto g_xxxz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 44);

    auto g_xxxz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 46);

    auto g_xxxz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 47);

    auto g_xxxz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 49);

    auto g_xxxz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 50);

    auto g_xxxz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 51);

    auto g_xxxz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 53);

    auto g_xxxz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 54);

    auto g_xxxz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 55);

    auto g_xxxz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 56);

    auto g_xxxz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 58);

    auto g_xxxz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 59);

    auto g_xxxz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 60);

    auto g_xxxz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 61);

    auto g_xxxz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 62);

    auto g_xxyy_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 63);

    auto g_xxyy_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 64);

    auto g_xxyy_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 65);

    auto g_xxyy_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 66);

    auto g_xxyy_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 67);

    auto g_xxyy_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 68);

    auto g_xxyy_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 69);

    auto g_xxyy_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 70);

    auto g_xxyy_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 71);

    auto g_xxyy_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 72);

    auto g_xxyy_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 73);

    auto g_xxyy_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 74);

    auto g_xxyy_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 75);

    auto g_xxyy_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 76);

    auto g_xxyy_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 77);

    auto g_xxyy_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 78);

    auto g_xxyy_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 79);

    auto g_xxyy_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 80);

    auto g_xxyy_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 81);

    auto g_xxyy_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 82);

    auto g_xxyy_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 83);

    auto g_xxzz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 105);

    auto g_xxzz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 106);

    auto g_xxzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 107);

    auto g_xxzz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 108);

    auto g_xxzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 109);

    auto g_xxzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 110);

    auto g_xxzz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 111);

    auto g_xxzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 112);

    auto g_xxzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 113);

    auto g_xxzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 114);

    auto g_xxzz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 115);

    auto g_xxzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 116);

    auto g_xxzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 117);

    auto g_xxzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 118);

    auto g_xxzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 119);

    auto g_xxzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 120);

    auto g_xxzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 121);

    auto g_xxzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 122);

    auto g_xxzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 123);

    auto g_xxzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 124);

    auto g_xxzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 125);

    auto g_xyyy_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 127);

    auto g_xyyy_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 129);

    auto g_xyyy_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 130);

    auto g_xyyy_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 132);

    auto g_xyyy_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 133);

    auto g_xyyy_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 134);

    auto g_xyyy_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 136);

    auto g_xyyy_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 137);

    auto g_xyyy_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 138);

    auto g_xyyy_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 139);

    auto g_xyyy_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 141);

    auto g_xyyy_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 142);

    auto g_xyyy_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 143);

    auto g_xyyy_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 144);

    auto g_xyyy_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 145);

    auto g_xzzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 191);

    auto g_xzzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 193);

    auto g_xzzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 194);

    auto g_xzzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 196);

    auto g_xzzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 197);

    auto g_xzzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 198);

    auto g_xzzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 200);

    auto g_xzzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 201);

    auto g_xzzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 202);

    auto g_xzzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 203);

    auto g_xzzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 205);

    auto g_xzzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 206);

    auto g_xzzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 207);

    auto g_xzzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 208);

    auto g_xzzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 209);

    auto g_yyyy_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 210);

    auto g_yyyy_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 211);

    auto g_yyyy_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 212);

    auto g_yyyy_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 213);

    auto g_yyyy_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 214);

    auto g_yyyy_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 215);

    auto g_yyyy_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 216);

    auto g_yyyy_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 217);

    auto g_yyyy_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 218);

    auto g_yyyy_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 219);

    auto g_yyyy_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 220);

    auto g_yyyy_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 221);

    auto g_yyyy_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 222);

    auto g_yyyy_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 223);

    auto g_yyyy_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 224);

    auto g_yyyy_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 225);

    auto g_yyyy_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 226);

    auto g_yyyy_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 227);

    auto g_yyyy_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 228);

    auto g_yyyy_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 229);

    auto g_yyyy_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 230);

    auto g_yyyz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 233);

    auto g_yyyz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 235);

    auto g_yyyz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 236);

    auto g_yyyz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 238);

    auto g_yyyz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 239);

    auto g_yyyz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 240);

    auto g_yyyz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 242);

    auto g_yyyz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 243);

    auto g_yyyz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 244);

    auto g_yyyz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 245);

    auto g_yyyz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 247);

    auto g_yyyz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 248);

    auto g_yyyz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 249);

    auto g_yyyz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 250);

    auto g_yyyz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 251);

    auto g_yyzz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 252);

    auto g_yyzz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 253);

    auto g_yyzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 254);

    auto g_yyzz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 255);

    auto g_yyzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 256);

    auto g_yyzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 257);

    auto g_yyzz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 258);

    auto g_yyzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 259);

    auto g_yyzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 260);

    auto g_yyzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 261);

    auto g_yyzz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 262);

    auto g_yyzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 263);

    auto g_yyzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 264);

    auto g_yyzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 265);

    auto g_yyzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 266);

    auto g_yyzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 267);

    auto g_yyzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 268);

    auto g_yyzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 269);

    auto g_yyzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 270);

    auto g_yyzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 271);

    auto g_yyzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 272);

    auto g_yzzz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 274);

    auto g_yzzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 275);

    auto g_yzzz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 276);

    auto g_yzzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 277);

    auto g_yzzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 278);

    auto g_yzzz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 279);

    auto g_yzzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 280);

    auto g_yzzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 281);

    auto g_yzzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 282);

    auto g_yzzz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 283);

    auto g_yzzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 284);

    auto g_yzzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 285);

    auto g_yzzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 286);

    auto g_yzzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 287);

    auto g_yzzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 288);

    auto g_yzzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 289);

    auto g_yzzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 290);

    auto g_yzzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 291);

    auto g_yzzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 292);

    auto g_yzzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 293);

    auto g_zzzz_xxxxx_1 = pbuffer.data(idx_eri_1_gh + 294);

    auto g_zzzz_xxxxy_1 = pbuffer.data(idx_eri_1_gh + 295);

    auto g_zzzz_xxxxz_1 = pbuffer.data(idx_eri_1_gh + 296);

    auto g_zzzz_xxxyy_1 = pbuffer.data(idx_eri_1_gh + 297);

    auto g_zzzz_xxxyz_1 = pbuffer.data(idx_eri_1_gh + 298);

    auto g_zzzz_xxxzz_1 = pbuffer.data(idx_eri_1_gh + 299);

    auto g_zzzz_xxyyy_1 = pbuffer.data(idx_eri_1_gh + 300);

    auto g_zzzz_xxyyz_1 = pbuffer.data(idx_eri_1_gh + 301);

    auto g_zzzz_xxyzz_1 = pbuffer.data(idx_eri_1_gh + 302);

    auto g_zzzz_xxzzz_1 = pbuffer.data(idx_eri_1_gh + 303);

    auto g_zzzz_xyyyy_1 = pbuffer.data(idx_eri_1_gh + 304);

    auto g_zzzz_xyyyz_1 = pbuffer.data(idx_eri_1_gh + 305);

    auto g_zzzz_xyyzz_1 = pbuffer.data(idx_eri_1_gh + 306);

    auto g_zzzz_xyzzz_1 = pbuffer.data(idx_eri_1_gh + 307);

    auto g_zzzz_xzzzz_1 = pbuffer.data(idx_eri_1_gh + 308);

    auto g_zzzz_yyyyy_1 = pbuffer.data(idx_eri_1_gh + 309);

    auto g_zzzz_yyyyz_1 = pbuffer.data(idx_eri_1_gh + 310);

    auto g_zzzz_yyyzz_1 = pbuffer.data(idx_eri_1_gh + 311);

    auto g_zzzz_yyzzz_1 = pbuffer.data(idx_eri_1_gh + 312);

    auto g_zzzz_yzzzz_1 = pbuffer.data(idx_eri_1_gh + 313);

    auto g_zzzz_zzzzz_1 = pbuffer.data(idx_eri_1_gh + 314);

    // Set up components of auxiliary buffer : GI

    auto g_xxxx_xxxxxx_1 = pbuffer.data(idx_eri_1_gi);

    auto g_xxxx_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 1);

    auto g_xxxx_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 2);

    auto g_xxxx_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 3);

    auto g_xxxx_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 4);

    auto g_xxxx_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 5);

    auto g_xxxx_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 6);

    auto g_xxxx_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 7);

    auto g_xxxx_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 8);

    auto g_xxxx_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 9);

    auto g_xxxx_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 10);

    auto g_xxxx_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 11);

    auto g_xxxx_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 12);

    auto g_xxxx_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 13);

    auto g_xxxx_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 14);

    auto g_xxxx_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 15);

    auto g_xxxx_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 16);

    auto g_xxxx_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 17);

    auto g_xxxx_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 18);

    auto g_xxxx_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 19);

    auto g_xxxx_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 20);

    auto g_xxxx_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 21);

    auto g_xxxx_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 22);

    auto g_xxxx_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 23);

    auto g_xxxx_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 24);

    auto g_xxxx_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 25);

    auto g_xxxx_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 26);

    auto g_xxxx_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 27);

    auto g_xxxy_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 28);

    auto g_xxxy_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 29);

    auto g_xxxy_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 30);

    auto g_xxxy_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 31);

    auto g_xxxy_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 33);

    auto g_xxxy_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 34);

    auto g_xxxy_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 37);

    auto g_xxxy_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 38);

    auto g_xxxy_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 42);

    auto g_xxxy_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 43);

    auto g_xxxy_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 48);

    auto g_xxxy_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 49);

    auto g_xxxz_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 56);

    auto g_xxxz_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 57);

    auto g_xxxz_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 58);

    auto g_xxxz_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 59);

    auto g_xxxz_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 60);

    auto g_xxxz_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 61);

    auto g_xxxz_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 62);

    auto g_xxxz_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 63);

    auto g_xxxz_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 64);

    auto g_xxxz_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 65);

    auto g_xxxz_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 66);

    auto g_xxxz_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 67);

    auto g_xxxz_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 68);

    auto g_xxxz_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 69);

    auto g_xxxz_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 70);

    auto g_xxxz_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 71);

    auto g_xxxz_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 72);

    auto g_xxxz_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 73);

    auto g_xxxz_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 74);

    auto g_xxxz_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 75);

    auto g_xxxz_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 76);

    auto g_xxxz_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 78);

    auto g_xxxz_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 79);

    auto g_xxxz_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 80);

    auto g_xxxz_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 81);

    auto g_xxxz_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 82);

    auto g_xxxz_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 83);

    auto g_xxyy_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 84);

    auto g_xxyy_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 85);

    auto g_xxyy_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 86);

    auto g_xxyy_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 87);

    auto g_xxyy_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 88);

    auto g_xxyy_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 89);

    auto g_xxyy_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 90);

    auto g_xxyy_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 91);

    auto g_xxyy_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 92);

    auto g_xxyy_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 93);

    auto g_xxyy_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 94);

    auto g_xxyy_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 95);

    auto g_xxyy_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 96);

    auto g_xxyy_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 97);

    auto g_xxyy_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 98);

    auto g_xxyy_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 99);

    auto g_xxyy_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 100);

    auto g_xxyy_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 101);

    auto g_xxyy_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 102);

    auto g_xxyy_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 103);

    auto g_xxyy_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 104);

    auto g_xxyy_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 105);

    auto g_xxyy_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 106);

    auto g_xxyy_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 107);

    auto g_xxyy_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 108);

    auto g_xxyy_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 109);

    auto g_xxyy_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 110);

    auto g_xxyy_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 111);

    auto g_xxzz_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 140);

    auto g_xxzz_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 141);

    auto g_xxzz_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 142);

    auto g_xxzz_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 143);

    auto g_xxzz_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 144);

    auto g_xxzz_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 145);

    auto g_xxzz_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 146);

    auto g_xxzz_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 147);

    auto g_xxzz_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 148);

    auto g_xxzz_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 149);

    auto g_xxzz_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 150);

    auto g_xxzz_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 151);

    auto g_xxzz_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 152);

    auto g_xxzz_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 153);

    auto g_xxzz_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 154);

    auto g_xxzz_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 155);

    auto g_xxzz_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 156);

    auto g_xxzz_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 157);

    auto g_xxzz_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 158);

    auto g_xxzz_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 159);

    auto g_xxzz_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 160);

    auto g_xxzz_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 161);

    auto g_xxzz_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 162);

    auto g_xxzz_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 163);

    auto g_xxzz_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 164);

    auto g_xxzz_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 165);

    auto g_xxzz_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 166);

    auto g_xxzz_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 167);

    auto g_xyyy_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 168);

    auto g_xyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 169);

    auto g_xyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 171);

    auto g_xyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 172);

    auto g_xyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 174);

    auto g_xyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 175);

    auto g_xyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 176);

    auto g_xyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 178);

    auto g_xyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 179);

    auto g_xyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 180);

    auto g_xyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 181);

    auto g_xyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 183);

    auto g_xyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 184);

    auto g_xyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 185);

    auto g_xyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 186);

    auto g_xyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 187);

    auto g_xyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 189);

    auto g_xyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 190);

    auto g_xyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 191);

    auto g_xyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 192);

    auto g_xyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 193);

    auto g_xyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 194);

    auto g_xyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 195);

    auto g_xzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 252);

    auto g_xzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 254);

    auto g_xzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 256);

    auto g_xzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 257);

    auto g_xzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 259);

    auto g_xzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 260);

    auto g_xzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 261);

    auto g_xzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 263);

    auto g_xzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 264);

    auto g_xzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 265);

    auto g_xzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 266);

    auto g_xzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 268);

    auto g_xzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 269);

    auto g_xzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 270);

    auto g_xzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 271);

    auto g_xzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 272);

    auto g_xzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 273);

    auto g_xzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 274);

    auto g_xzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 275);

    auto g_xzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 276);

    auto g_xzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 277);

    auto g_xzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 278);

    auto g_xzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 279);

    auto g_yyyy_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 280);

    auto g_yyyy_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 281);

    auto g_yyyy_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 282);

    auto g_yyyy_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 283);

    auto g_yyyy_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 284);

    auto g_yyyy_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 285);

    auto g_yyyy_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 286);

    auto g_yyyy_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 287);

    auto g_yyyy_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 288);

    auto g_yyyy_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 289);

    auto g_yyyy_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 290);

    auto g_yyyy_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 291);

    auto g_yyyy_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 292);

    auto g_yyyy_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 293);

    auto g_yyyy_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 294);

    auto g_yyyy_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 295);

    auto g_yyyy_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 296);

    auto g_yyyy_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 297);

    auto g_yyyy_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 298);

    auto g_yyyy_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 299);

    auto g_yyyy_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 300);

    auto g_yyyy_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 301);

    auto g_yyyy_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 302);

    auto g_yyyy_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 303);

    auto g_yyyy_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 304);

    auto g_yyyy_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 305);

    auto g_yyyy_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 306);

    auto g_yyyy_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 307);

    auto g_yyyz_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 309);

    auto g_yyyz_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 310);

    auto g_yyyz_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 311);

    auto g_yyyz_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 312);

    auto g_yyyz_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 313);

    auto g_yyyz_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 314);

    auto g_yyyz_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 315);

    auto g_yyyz_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 316);

    auto g_yyyz_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 317);

    auto g_yyyz_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 318);

    auto g_yyyz_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 319);

    auto g_yyyz_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 320);

    auto g_yyyz_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 321);

    auto g_yyyz_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 322);

    auto g_yyyz_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 323);

    auto g_yyyz_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 324);

    auto g_yyyz_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 325);

    auto g_yyyz_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 326);

    auto g_yyyz_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 327);

    auto g_yyyz_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 328);

    auto g_yyyz_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 329);

    auto g_yyyz_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 330);

    auto g_yyyz_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 331);

    auto g_yyyz_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 332);

    auto g_yyyz_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 333);

    auto g_yyyz_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 334);

    auto g_yyyz_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 335);

    auto g_yyzz_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 336);

    auto g_yyzz_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 337);

    auto g_yyzz_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 338);

    auto g_yyzz_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 339);

    auto g_yyzz_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 340);

    auto g_yyzz_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 341);

    auto g_yyzz_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 342);

    auto g_yyzz_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 343);

    auto g_yyzz_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 344);

    auto g_yyzz_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 345);

    auto g_yyzz_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 346);

    auto g_yyzz_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 347);

    auto g_yyzz_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 348);

    auto g_yyzz_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 349);

    auto g_yyzz_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 350);

    auto g_yyzz_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 351);

    auto g_yyzz_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 352);

    auto g_yyzz_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 353);

    auto g_yyzz_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 354);

    auto g_yyzz_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 355);

    auto g_yyzz_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 356);

    auto g_yyzz_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 357);

    auto g_yyzz_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 358);

    auto g_yyzz_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 359);

    auto g_yyzz_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 360);

    auto g_yyzz_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 361);

    auto g_yyzz_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 362);

    auto g_yyzz_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 363);

    auto g_yzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 364);

    auto g_yzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 365);

    auto g_yzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 366);

    auto g_yzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 367);

    auto g_yzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 368);

    auto g_yzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 369);

    auto g_yzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 370);

    auto g_yzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 371);

    auto g_yzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 372);

    auto g_yzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 373);

    auto g_yzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 374);

    auto g_yzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 375);

    auto g_yzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 376);

    auto g_yzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 377);

    auto g_yzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 378);

    auto g_yzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 379);

    auto g_yzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 380);

    auto g_yzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 381);

    auto g_yzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 382);

    auto g_yzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 383);

    auto g_yzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 384);

    auto g_yzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 385);

    auto g_yzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 386);

    auto g_yzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 387);

    auto g_yzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 388);

    auto g_yzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 389);

    auto g_yzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 390);

    auto g_yzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 391);

    auto g_zzzz_xxxxxx_1 = pbuffer.data(idx_eri_1_gi + 392);

    auto g_zzzz_xxxxxy_1 = pbuffer.data(idx_eri_1_gi + 393);

    auto g_zzzz_xxxxxz_1 = pbuffer.data(idx_eri_1_gi + 394);

    auto g_zzzz_xxxxyy_1 = pbuffer.data(idx_eri_1_gi + 395);

    auto g_zzzz_xxxxyz_1 = pbuffer.data(idx_eri_1_gi + 396);

    auto g_zzzz_xxxxzz_1 = pbuffer.data(idx_eri_1_gi + 397);

    auto g_zzzz_xxxyyy_1 = pbuffer.data(idx_eri_1_gi + 398);

    auto g_zzzz_xxxyyz_1 = pbuffer.data(idx_eri_1_gi + 399);

    auto g_zzzz_xxxyzz_1 = pbuffer.data(idx_eri_1_gi + 400);

    auto g_zzzz_xxxzzz_1 = pbuffer.data(idx_eri_1_gi + 401);

    auto g_zzzz_xxyyyy_1 = pbuffer.data(idx_eri_1_gi + 402);

    auto g_zzzz_xxyyyz_1 = pbuffer.data(idx_eri_1_gi + 403);

    auto g_zzzz_xxyyzz_1 = pbuffer.data(idx_eri_1_gi + 404);

    auto g_zzzz_xxyzzz_1 = pbuffer.data(idx_eri_1_gi + 405);

    auto g_zzzz_xxzzzz_1 = pbuffer.data(idx_eri_1_gi + 406);

    auto g_zzzz_xyyyyy_1 = pbuffer.data(idx_eri_1_gi + 407);

    auto g_zzzz_xyyyyz_1 = pbuffer.data(idx_eri_1_gi + 408);

    auto g_zzzz_xyyyzz_1 = pbuffer.data(idx_eri_1_gi + 409);

    auto g_zzzz_xyyzzz_1 = pbuffer.data(idx_eri_1_gi + 410);

    auto g_zzzz_xyzzzz_1 = pbuffer.data(idx_eri_1_gi + 411);

    auto g_zzzz_xzzzzz_1 = pbuffer.data(idx_eri_1_gi + 412);

    auto g_zzzz_yyyyyy_1 = pbuffer.data(idx_eri_1_gi + 413);

    auto g_zzzz_yyyyyz_1 = pbuffer.data(idx_eri_1_gi + 414);

    auto g_zzzz_yyyyzz_1 = pbuffer.data(idx_eri_1_gi + 415);

    auto g_zzzz_yyyzzz_1 = pbuffer.data(idx_eri_1_gi + 416);

    auto g_zzzz_yyzzzz_1 = pbuffer.data(idx_eri_1_gi + 417);

    auto g_zzzz_yzzzzz_1 = pbuffer.data(idx_eri_1_gi + 418);

    auto g_zzzz_zzzzzz_1 = pbuffer.data(idx_eri_1_gi + 419);

    // Set up 0-28 components of targeted buffer : HI

    auto g_xxxxx_xxxxxx_0 = pbuffer.data(idx_eri_0_hi);

    auto g_xxxxx_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 1);

    auto g_xxxxx_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 2);

    auto g_xxxxx_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 3);

    auto g_xxxxx_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 4);

    auto g_xxxxx_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 5);

    auto g_xxxxx_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 6);

    auto g_xxxxx_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 7);

    auto g_xxxxx_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 8);

    auto g_xxxxx_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 9);

    auto g_xxxxx_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 10);

    auto g_xxxxx_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 11);

    auto g_xxxxx_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 12);

    auto g_xxxxx_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 13);

    auto g_xxxxx_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 14);

    auto g_xxxxx_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 15);

    auto g_xxxxx_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 16);

    auto g_xxxxx_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 17);

    auto g_xxxxx_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 18);

    auto g_xxxxx_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 19);

    auto g_xxxxx_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 20);

    auto g_xxxxx_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 21);

    auto g_xxxxx_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 22);

    auto g_xxxxx_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 23);

    auto g_xxxxx_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 24);

    auto g_xxxxx_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 25);

    auto g_xxxxx_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 26);

    auto g_xxxxx_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 27);

    #pragma omp simd aligned(g_xxx_xxxxxx_0, g_xxx_xxxxxx_1, g_xxx_xxxxxy_0, g_xxx_xxxxxy_1, g_xxx_xxxxxz_0, g_xxx_xxxxxz_1, g_xxx_xxxxyy_0, g_xxx_xxxxyy_1, g_xxx_xxxxyz_0, g_xxx_xxxxyz_1, g_xxx_xxxxzz_0, g_xxx_xxxxzz_1, g_xxx_xxxyyy_0, g_xxx_xxxyyy_1, g_xxx_xxxyyz_0, g_xxx_xxxyyz_1, g_xxx_xxxyzz_0, g_xxx_xxxyzz_1, g_xxx_xxxzzz_0, g_xxx_xxxzzz_1, g_xxx_xxyyyy_0, g_xxx_xxyyyy_1, g_xxx_xxyyyz_0, g_xxx_xxyyyz_1, g_xxx_xxyyzz_0, g_xxx_xxyyzz_1, g_xxx_xxyzzz_0, g_xxx_xxyzzz_1, g_xxx_xxzzzz_0, g_xxx_xxzzzz_1, g_xxx_xyyyyy_0, g_xxx_xyyyyy_1, g_xxx_xyyyyz_0, g_xxx_xyyyyz_1, g_xxx_xyyyzz_0, g_xxx_xyyyzz_1, g_xxx_xyyzzz_0, g_xxx_xyyzzz_1, g_xxx_xyzzzz_0, g_xxx_xyzzzz_1, g_xxx_xzzzzz_0, g_xxx_xzzzzz_1, g_xxx_yyyyyy_0, g_xxx_yyyyyy_1, g_xxx_yyyyyz_0, g_xxx_yyyyyz_1, g_xxx_yyyyzz_0, g_xxx_yyyyzz_1, g_xxx_yyyzzz_0, g_xxx_yyyzzz_1, g_xxx_yyzzzz_0, g_xxx_yyzzzz_1, g_xxx_yzzzzz_0, g_xxx_yzzzzz_1, g_xxx_zzzzzz_0, g_xxx_zzzzzz_1, g_xxxx_xxxxx_1, g_xxxx_xxxxxx_1, g_xxxx_xxxxxy_1, g_xxxx_xxxxxz_1, g_xxxx_xxxxy_1, g_xxxx_xxxxyy_1, g_xxxx_xxxxyz_1, g_xxxx_xxxxz_1, g_xxxx_xxxxzz_1, g_xxxx_xxxyy_1, g_xxxx_xxxyyy_1, g_xxxx_xxxyyz_1, g_xxxx_xxxyz_1, g_xxxx_xxxyzz_1, g_xxxx_xxxzz_1, g_xxxx_xxxzzz_1, g_xxxx_xxyyy_1, g_xxxx_xxyyyy_1, g_xxxx_xxyyyz_1, g_xxxx_xxyyz_1, g_xxxx_xxyyzz_1, g_xxxx_xxyzz_1, g_xxxx_xxyzzz_1, g_xxxx_xxzzz_1, g_xxxx_xxzzzz_1, g_xxxx_xyyyy_1, g_xxxx_xyyyyy_1, g_xxxx_xyyyyz_1, g_xxxx_xyyyz_1, g_xxxx_xyyyzz_1, g_xxxx_xyyzz_1, g_xxxx_xyyzzz_1, g_xxxx_xyzzz_1, g_xxxx_xyzzzz_1, g_xxxx_xzzzz_1, g_xxxx_xzzzzz_1, g_xxxx_yyyyy_1, g_xxxx_yyyyyy_1, g_xxxx_yyyyyz_1, g_xxxx_yyyyz_1, g_xxxx_yyyyzz_1, g_xxxx_yyyzz_1, g_xxxx_yyyzzz_1, g_xxxx_yyzzz_1, g_xxxx_yyzzzz_1, g_xxxx_yzzzz_1, g_xxxx_yzzzzz_1, g_xxxx_zzzzz_1, g_xxxx_zzzzzz_1, g_xxxxx_xxxxxx_0, g_xxxxx_xxxxxy_0, g_xxxxx_xxxxxz_0, g_xxxxx_xxxxyy_0, g_xxxxx_xxxxyz_0, g_xxxxx_xxxxzz_0, g_xxxxx_xxxyyy_0, g_xxxxx_xxxyyz_0, g_xxxxx_xxxyzz_0, g_xxxxx_xxxzzz_0, g_xxxxx_xxyyyy_0, g_xxxxx_xxyyyz_0, g_xxxxx_xxyyzz_0, g_xxxxx_xxyzzz_0, g_xxxxx_xxzzzz_0, g_xxxxx_xyyyyy_0, g_xxxxx_xyyyyz_0, g_xxxxx_xyyyzz_0, g_xxxxx_xyyzzz_0, g_xxxxx_xyzzzz_0, g_xxxxx_xzzzzz_0, g_xxxxx_yyyyyy_0, g_xxxxx_yyyyyz_0, g_xxxxx_yyyyzz_0, g_xxxxx_yyyzzz_0, g_xxxxx_yyzzzz_0, g_xxxxx_yzzzzz_0, g_xxxxx_zzzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxx_xxxxxx_0[i] = 4.0 * g_xxx_xxxxxx_0[i] * fbe_0 - 4.0 * g_xxx_xxxxxx_1[i] * fz_be_0 + 6.0 * g_xxxx_xxxxx_1[i] * fe_0 + g_xxxx_xxxxxx_1[i] * pa_x[i];

        g_xxxxx_xxxxxy_0[i] = 4.0 * g_xxx_xxxxxy_0[i] * fbe_0 - 4.0 * g_xxx_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxxx_xxxxy_1[i] * fe_0 + g_xxxx_xxxxxy_1[i] * pa_x[i];

        g_xxxxx_xxxxxz_0[i] = 4.0 * g_xxx_xxxxxz_0[i] * fbe_0 - 4.0 * g_xxx_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxxx_xxxxz_1[i] * fe_0 + g_xxxx_xxxxxz_1[i] * pa_x[i];

        g_xxxxx_xxxxyy_0[i] = 4.0 * g_xxx_xxxxyy_0[i] * fbe_0 - 4.0 * g_xxx_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxxx_xxxyy_1[i] * fe_0 + g_xxxx_xxxxyy_1[i] * pa_x[i];

        g_xxxxx_xxxxyz_0[i] = 4.0 * g_xxx_xxxxyz_0[i] * fbe_0 - 4.0 * g_xxx_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxxx_xxxyz_1[i] * fe_0 + g_xxxx_xxxxyz_1[i] * pa_x[i];

        g_xxxxx_xxxxzz_0[i] = 4.0 * g_xxx_xxxxzz_0[i] * fbe_0 - 4.0 * g_xxx_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxxx_xxxzz_1[i] * fe_0 + g_xxxx_xxxxzz_1[i] * pa_x[i];

        g_xxxxx_xxxyyy_0[i] = 4.0 * g_xxx_xxxyyy_0[i] * fbe_0 - 4.0 * g_xxx_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxxx_xxyyy_1[i] * fe_0 + g_xxxx_xxxyyy_1[i] * pa_x[i];

        g_xxxxx_xxxyyz_0[i] = 4.0 * g_xxx_xxxyyz_0[i] * fbe_0 - 4.0 * g_xxx_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxxx_xxyyz_1[i] * fe_0 + g_xxxx_xxxyyz_1[i] * pa_x[i];

        g_xxxxx_xxxyzz_0[i] = 4.0 * g_xxx_xxxyzz_0[i] * fbe_0 - 4.0 * g_xxx_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxxx_xxyzz_1[i] * fe_0 + g_xxxx_xxxyzz_1[i] * pa_x[i];

        g_xxxxx_xxxzzz_0[i] = 4.0 * g_xxx_xxxzzz_0[i] * fbe_0 - 4.0 * g_xxx_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxxx_xxzzz_1[i] * fe_0 + g_xxxx_xxxzzz_1[i] * pa_x[i];

        g_xxxxx_xxyyyy_0[i] = 4.0 * g_xxx_xxyyyy_0[i] * fbe_0 - 4.0 * g_xxx_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxxx_xyyyy_1[i] * fe_0 + g_xxxx_xxyyyy_1[i] * pa_x[i];

        g_xxxxx_xxyyyz_0[i] = 4.0 * g_xxx_xxyyyz_0[i] * fbe_0 - 4.0 * g_xxx_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxxx_xyyyz_1[i] * fe_0 + g_xxxx_xxyyyz_1[i] * pa_x[i];

        g_xxxxx_xxyyzz_0[i] = 4.0 * g_xxx_xxyyzz_0[i] * fbe_0 - 4.0 * g_xxx_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxxx_xyyzz_1[i] * fe_0 + g_xxxx_xxyyzz_1[i] * pa_x[i];

        g_xxxxx_xxyzzz_0[i] = 4.0 * g_xxx_xxyzzz_0[i] * fbe_0 - 4.0 * g_xxx_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_xyzzz_1[i] * fe_0 + g_xxxx_xxyzzz_1[i] * pa_x[i];

        g_xxxxx_xxzzzz_0[i] = 4.0 * g_xxx_xxzzzz_0[i] * fbe_0 - 4.0 * g_xxx_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxxx_xzzzz_1[i] * fe_0 + g_xxxx_xxzzzz_1[i] * pa_x[i];

        g_xxxxx_xyyyyy_0[i] = 4.0 * g_xxx_xyyyyy_0[i] * fbe_0 - 4.0 * g_xxx_xyyyyy_1[i] * fz_be_0 + g_xxxx_yyyyy_1[i] * fe_0 + g_xxxx_xyyyyy_1[i] * pa_x[i];

        g_xxxxx_xyyyyz_0[i] = 4.0 * g_xxx_xyyyyz_0[i] * fbe_0 - 4.0 * g_xxx_xyyyyz_1[i] * fz_be_0 + g_xxxx_yyyyz_1[i] * fe_0 + g_xxxx_xyyyyz_1[i] * pa_x[i];

        g_xxxxx_xyyyzz_0[i] = 4.0 * g_xxx_xyyyzz_0[i] * fbe_0 - 4.0 * g_xxx_xyyyzz_1[i] * fz_be_0 + g_xxxx_yyyzz_1[i] * fe_0 + g_xxxx_xyyyzz_1[i] * pa_x[i];

        g_xxxxx_xyyzzz_0[i] = 4.0 * g_xxx_xyyzzz_0[i] * fbe_0 - 4.0 * g_xxx_xyyzzz_1[i] * fz_be_0 + g_xxxx_yyzzz_1[i] * fe_0 + g_xxxx_xyyzzz_1[i] * pa_x[i];

        g_xxxxx_xyzzzz_0[i] = 4.0 * g_xxx_xyzzzz_0[i] * fbe_0 - 4.0 * g_xxx_xyzzzz_1[i] * fz_be_0 + g_xxxx_yzzzz_1[i] * fe_0 + g_xxxx_xyzzzz_1[i] * pa_x[i];

        g_xxxxx_xzzzzz_0[i] = 4.0 * g_xxx_xzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_xzzzzz_1[i] * fz_be_0 + g_xxxx_zzzzz_1[i] * fe_0 + g_xxxx_xzzzzz_1[i] * pa_x[i];

        g_xxxxx_yyyyyy_0[i] = 4.0 * g_xxx_yyyyyy_0[i] * fbe_0 - 4.0 * g_xxx_yyyyyy_1[i] * fz_be_0 + g_xxxx_yyyyyy_1[i] * pa_x[i];

        g_xxxxx_yyyyyz_0[i] = 4.0 * g_xxx_yyyyyz_0[i] * fbe_0 - 4.0 * g_xxx_yyyyyz_1[i] * fz_be_0 + g_xxxx_yyyyyz_1[i] * pa_x[i];

        g_xxxxx_yyyyzz_0[i] = 4.0 * g_xxx_yyyyzz_0[i] * fbe_0 - 4.0 * g_xxx_yyyyzz_1[i] * fz_be_0 + g_xxxx_yyyyzz_1[i] * pa_x[i];

        g_xxxxx_yyyzzz_0[i] = 4.0 * g_xxx_yyyzzz_0[i] * fbe_0 - 4.0 * g_xxx_yyyzzz_1[i] * fz_be_0 + g_xxxx_yyyzzz_1[i] * pa_x[i];

        g_xxxxx_yyzzzz_0[i] = 4.0 * g_xxx_yyzzzz_0[i] * fbe_0 - 4.0 * g_xxx_yyzzzz_1[i] * fz_be_0 + g_xxxx_yyzzzz_1[i] * pa_x[i];

        g_xxxxx_yzzzzz_0[i] = 4.0 * g_xxx_yzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_yzzzzz_1[i] * fz_be_0 + g_xxxx_yzzzzz_1[i] * pa_x[i];

        g_xxxxx_zzzzzz_0[i] = 4.0 * g_xxx_zzzzzz_0[i] * fbe_0 - 4.0 * g_xxx_zzzzzz_1[i] * fz_be_0 + g_xxxx_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : HI

    auto g_xxxxy_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 28);

    auto g_xxxxy_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 29);

    auto g_xxxxy_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 30);

    auto g_xxxxy_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 31);

    auto g_xxxxy_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 32);

    auto g_xxxxy_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 33);

    auto g_xxxxy_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 34);

    auto g_xxxxy_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 35);

    auto g_xxxxy_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 36);

    auto g_xxxxy_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 37);

    auto g_xxxxy_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 38);

    auto g_xxxxy_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 39);

    auto g_xxxxy_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 40);

    auto g_xxxxy_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 41);

    auto g_xxxxy_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 42);

    auto g_xxxxy_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 43);

    auto g_xxxxy_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 44);

    auto g_xxxxy_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 45);

    auto g_xxxxy_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 46);

    auto g_xxxxy_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 47);

    auto g_xxxxy_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 48);

    auto g_xxxxy_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 49);

    auto g_xxxxy_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 50);

    auto g_xxxxy_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 51);

    auto g_xxxxy_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 52);

    auto g_xxxxy_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 53);

    auto g_xxxxy_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 54);

    auto g_xxxxy_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 55);

    #pragma omp simd aligned(g_xxxx_xxxxx_1, g_xxxx_xxxxxx_1, g_xxxx_xxxxxy_1, g_xxxx_xxxxxz_1, g_xxxx_xxxxy_1, g_xxxx_xxxxyy_1, g_xxxx_xxxxyz_1, g_xxxx_xxxxz_1, g_xxxx_xxxxzz_1, g_xxxx_xxxyy_1, g_xxxx_xxxyyy_1, g_xxxx_xxxyyz_1, g_xxxx_xxxyz_1, g_xxxx_xxxyzz_1, g_xxxx_xxxzz_1, g_xxxx_xxxzzz_1, g_xxxx_xxyyy_1, g_xxxx_xxyyyy_1, g_xxxx_xxyyyz_1, g_xxxx_xxyyz_1, g_xxxx_xxyyzz_1, g_xxxx_xxyzz_1, g_xxxx_xxyzzz_1, g_xxxx_xxzzz_1, g_xxxx_xxzzzz_1, g_xxxx_xyyyy_1, g_xxxx_xyyyyy_1, g_xxxx_xyyyyz_1, g_xxxx_xyyyz_1, g_xxxx_xyyyzz_1, g_xxxx_xyyzz_1, g_xxxx_xyyzzz_1, g_xxxx_xyzzz_1, g_xxxx_xyzzzz_1, g_xxxx_xzzzz_1, g_xxxx_xzzzzz_1, g_xxxx_yyyyy_1, g_xxxx_yyyyyy_1, g_xxxx_yyyyyz_1, g_xxxx_yyyyz_1, g_xxxx_yyyyzz_1, g_xxxx_yyyzz_1, g_xxxx_yyyzzz_1, g_xxxx_yyzzz_1, g_xxxx_yyzzzz_1, g_xxxx_yzzzz_1, g_xxxx_yzzzzz_1, g_xxxx_zzzzz_1, g_xxxx_zzzzzz_1, g_xxxxy_xxxxxx_0, g_xxxxy_xxxxxy_0, g_xxxxy_xxxxxz_0, g_xxxxy_xxxxyy_0, g_xxxxy_xxxxyz_0, g_xxxxy_xxxxzz_0, g_xxxxy_xxxyyy_0, g_xxxxy_xxxyyz_0, g_xxxxy_xxxyzz_0, g_xxxxy_xxxzzz_0, g_xxxxy_xxyyyy_0, g_xxxxy_xxyyyz_0, g_xxxxy_xxyyzz_0, g_xxxxy_xxyzzz_0, g_xxxxy_xxzzzz_0, g_xxxxy_xyyyyy_0, g_xxxxy_xyyyyz_0, g_xxxxy_xyyyzz_0, g_xxxxy_xyyzzz_0, g_xxxxy_xyzzzz_0, g_xxxxy_xzzzzz_0, g_xxxxy_yyyyyy_0, g_xxxxy_yyyyyz_0, g_xxxxy_yyyyzz_0, g_xxxxy_yyyzzz_0, g_xxxxy_yyzzzz_0, g_xxxxy_yzzzzz_0, g_xxxxy_zzzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxy_xxxxxx_0[i] = g_xxxx_xxxxxx_1[i] * pa_y[i];

        g_xxxxy_xxxxxy_0[i] = g_xxxx_xxxxx_1[i] * fe_0 + g_xxxx_xxxxxy_1[i] * pa_y[i];

        g_xxxxy_xxxxxz_0[i] = g_xxxx_xxxxxz_1[i] * pa_y[i];

        g_xxxxy_xxxxyy_0[i] = 2.0 * g_xxxx_xxxxy_1[i] * fe_0 + g_xxxx_xxxxyy_1[i] * pa_y[i];

        g_xxxxy_xxxxyz_0[i] = g_xxxx_xxxxz_1[i] * fe_0 + g_xxxx_xxxxyz_1[i] * pa_y[i];

        g_xxxxy_xxxxzz_0[i] = g_xxxx_xxxxzz_1[i] * pa_y[i];

        g_xxxxy_xxxyyy_0[i] = 3.0 * g_xxxx_xxxyy_1[i] * fe_0 + g_xxxx_xxxyyy_1[i] * pa_y[i];

        g_xxxxy_xxxyyz_0[i] = 2.0 * g_xxxx_xxxyz_1[i] * fe_0 + g_xxxx_xxxyyz_1[i] * pa_y[i];

        g_xxxxy_xxxyzz_0[i] = g_xxxx_xxxzz_1[i] * fe_0 + g_xxxx_xxxyzz_1[i] * pa_y[i];

        g_xxxxy_xxxzzz_0[i] = g_xxxx_xxxzzz_1[i] * pa_y[i];

        g_xxxxy_xxyyyy_0[i] = 4.0 * g_xxxx_xxyyy_1[i] * fe_0 + g_xxxx_xxyyyy_1[i] * pa_y[i];

        g_xxxxy_xxyyyz_0[i] = 3.0 * g_xxxx_xxyyz_1[i] * fe_0 + g_xxxx_xxyyyz_1[i] * pa_y[i];

        g_xxxxy_xxyyzz_0[i] = 2.0 * g_xxxx_xxyzz_1[i] * fe_0 + g_xxxx_xxyyzz_1[i] * pa_y[i];

        g_xxxxy_xxyzzz_0[i] = g_xxxx_xxzzz_1[i] * fe_0 + g_xxxx_xxyzzz_1[i] * pa_y[i];

        g_xxxxy_xxzzzz_0[i] = g_xxxx_xxzzzz_1[i] * pa_y[i];

        g_xxxxy_xyyyyy_0[i] = 5.0 * g_xxxx_xyyyy_1[i] * fe_0 + g_xxxx_xyyyyy_1[i] * pa_y[i];

        g_xxxxy_xyyyyz_0[i] = 4.0 * g_xxxx_xyyyz_1[i] * fe_0 + g_xxxx_xyyyyz_1[i] * pa_y[i];

        g_xxxxy_xyyyzz_0[i] = 3.0 * g_xxxx_xyyzz_1[i] * fe_0 + g_xxxx_xyyyzz_1[i] * pa_y[i];

        g_xxxxy_xyyzzz_0[i] = 2.0 * g_xxxx_xyzzz_1[i] * fe_0 + g_xxxx_xyyzzz_1[i] * pa_y[i];

        g_xxxxy_xyzzzz_0[i] = g_xxxx_xzzzz_1[i] * fe_0 + g_xxxx_xyzzzz_1[i] * pa_y[i];

        g_xxxxy_xzzzzz_0[i] = g_xxxx_xzzzzz_1[i] * pa_y[i];

        g_xxxxy_yyyyyy_0[i] = 6.0 * g_xxxx_yyyyy_1[i] * fe_0 + g_xxxx_yyyyyy_1[i] * pa_y[i];

        g_xxxxy_yyyyyz_0[i] = 5.0 * g_xxxx_yyyyz_1[i] * fe_0 + g_xxxx_yyyyyz_1[i] * pa_y[i];

        g_xxxxy_yyyyzz_0[i] = 4.0 * g_xxxx_yyyzz_1[i] * fe_0 + g_xxxx_yyyyzz_1[i] * pa_y[i];

        g_xxxxy_yyyzzz_0[i] = 3.0 * g_xxxx_yyzzz_1[i] * fe_0 + g_xxxx_yyyzzz_1[i] * pa_y[i];

        g_xxxxy_yyzzzz_0[i] = 2.0 * g_xxxx_yzzzz_1[i] * fe_0 + g_xxxx_yyzzzz_1[i] * pa_y[i];

        g_xxxxy_yzzzzz_0[i] = g_xxxx_zzzzz_1[i] * fe_0 + g_xxxx_yzzzzz_1[i] * pa_y[i];

        g_xxxxy_zzzzzz_0[i] = g_xxxx_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : HI

    auto g_xxxxz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 56);

    auto g_xxxxz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 57);

    auto g_xxxxz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 58);

    auto g_xxxxz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 59);

    auto g_xxxxz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 60);

    auto g_xxxxz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 61);

    auto g_xxxxz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 62);

    auto g_xxxxz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 63);

    auto g_xxxxz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 64);

    auto g_xxxxz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 65);

    auto g_xxxxz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 66);

    auto g_xxxxz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 67);

    auto g_xxxxz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 68);

    auto g_xxxxz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 69);

    auto g_xxxxz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 70);

    auto g_xxxxz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 71);

    auto g_xxxxz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 72);

    auto g_xxxxz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 73);

    auto g_xxxxz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 74);

    auto g_xxxxz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 75);

    auto g_xxxxz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 76);

    auto g_xxxxz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 77);

    auto g_xxxxz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 78);

    auto g_xxxxz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 79);

    auto g_xxxxz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 80);

    auto g_xxxxz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 81);

    auto g_xxxxz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 82);

    auto g_xxxxz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 83);

    #pragma omp simd aligned(g_xxxx_xxxxx_1, g_xxxx_xxxxxx_1, g_xxxx_xxxxxy_1, g_xxxx_xxxxxz_1, g_xxxx_xxxxy_1, g_xxxx_xxxxyy_1, g_xxxx_xxxxyz_1, g_xxxx_xxxxz_1, g_xxxx_xxxxzz_1, g_xxxx_xxxyy_1, g_xxxx_xxxyyy_1, g_xxxx_xxxyyz_1, g_xxxx_xxxyz_1, g_xxxx_xxxyzz_1, g_xxxx_xxxzz_1, g_xxxx_xxxzzz_1, g_xxxx_xxyyy_1, g_xxxx_xxyyyy_1, g_xxxx_xxyyyz_1, g_xxxx_xxyyz_1, g_xxxx_xxyyzz_1, g_xxxx_xxyzz_1, g_xxxx_xxyzzz_1, g_xxxx_xxzzz_1, g_xxxx_xxzzzz_1, g_xxxx_xyyyy_1, g_xxxx_xyyyyy_1, g_xxxx_xyyyyz_1, g_xxxx_xyyyz_1, g_xxxx_xyyyzz_1, g_xxxx_xyyzz_1, g_xxxx_xyyzzz_1, g_xxxx_xyzzz_1, g_xxxx_xyzzzz_1, g_xxxx_xzzzz_1, g_xxxx_xzzzzz_1, g_xxxx_yyyyy_1, g_xxxx_yyyyyy_1, g_xxxx_yyyyyz_1, g_xxxx_yyyyz_1, g_xxxx_yyyyzz_1, g_xxxx_yyyzz_1, g_xxxx_yyyzzz_1, g_xxxx_yyzzz_1, g_xxxx_yyzzzz_1, g_xxxx_yzzzz_1, g_xxxx_yzzzzz_1, g_xxxx_zzzzz_1, g_xxxx_zzzzzz_1, g_xxxxz_xxxxxx_0, g_xxxxz_xxxxxy_0, g_xxxxz_xxxxxz_0, g_xxxxz_xxxxyy_0, g_xxxxz_xxxxyz_0, g_xxxxz_xxxxzz_0, g_xxxxz_xxxyyy_0, g_xxxxz_xxxyyz_0, g_xxxxz_xxxyzz_0, g_xxxxz_xxxzzz_0, g_xxxxz_xxyyyy_0, g_xxxxz_xxyyyz_0, g_xxxxz_xxyyzz_0, g_xxxxz_xxyzzz_0, g_xxxxz_xxzzzz_0, g_xxxxz_xyyyyy_0, g_xxxxz_xyyyyz_0, g_xxxxz_xyyyzz_0, g_xxxxz_xyyzzz_0, g_xxxxz_xyzzzz_0, g_xxxxz_xzzzzz_0, g_xxxxz_yyyyyy_0, g_xxxxz_yyyyyz_0, g_xxxxz_yyyyzz_0, g_xxxxz_yyyzzz_0, g_xxxxz_yyzzzz_0, g_xxxxz_yzzzzz_0, g_xxxxz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxz_xxxxxx_0[i] = g_xxxx_xxxxxx_1[i] * pa_z[i];

        g_xxxxz_xxxxxy_0[i] = g_xxxx_xxxxxy_1[i] * pa_z[i];

        g_xxxxz_xxxxxz_0[i] = g_xxxx_xxxxx_1[i] * fe_0 + g_xxxx_xxxxxz_1[i] * pa_z[i];

        g_xxxxz_xxxxyy_0[i] = g_xxxx_xxxxyy_1[i] * pa_z[i];

        g_xxxxz_xxxxyz_0[i] = g_xxxx_xxxxy_1[i] * fe_0 + g_xxxx_xxxxyz_1[i] * pa_z[i];

        g_xxxxz_xxxxzz_0[i] = 2.0 * g_xxxx_xxxxz_1[i] * fe_0 + g_xxxx_xxxxzz_1[i] * pa_z[i];

        g_xxxxz_xxxyyy_0[i] = g_xxxx_xxxyyy_1[i] * pa_z[i];

        g_xxxxz_xxxyyz_0[i] = g_xxxx_xxxyy_1[i] * fe_0 + g_xxxx_xxxyyz_1[i] * pa_z[i];

        g_xxxxz_xxxyzz_0[i] = 2.0 * g_xxxx_xxxyz_1[i] * fe_0 + g_xxxx_xxxyzz_1[i] * pa_z[i];

        g_xxxxz_xxxzzz_0[i] = 3.0 * g_xxxx_xxxzz_1[i] * fe_0 + g_xxxx_xxxzzz_1[i] * pa_z[i];

        g_xxxxz_xxyyyy_0[i] = g_xxxx_xxyyyy_1[i] * pa_z[i];

        g_xxxxz_xxyyyz_0[i] = g_xxxx_xxyyy_1[i] * fe_0 + g_xxxx_xxyyyz_1[i] * pa_z[i];

        g_xxxxz_xxyyzz_0[i] = 2.0 * g_xxxx_xxyyz_1[i] * fe_0 + g_xxxx_xxyyzz_1[i] * pa_z[i];

        g_xxxxz_xxyzzz_0[i] = 3.0 * g_xxxx_xxyzz_1[i] * fe_0 + g_xxxx_xxyzzz_1[i] * pa_z[i];

        g_xxxxz_xxzzzz_0[i] = 4.0 * g_xxxx_xxzzz_1[i] * fe_0 + g_xxxx_xxzzzz_1[i] * pa_z[i];

        g_xxxxz_xyyyyy_0[i] = g_xxxx_xyyyyy_1[i] * pa_z[i];

        g_xxxxz_xyyyyz_0[i] = g_xxxx_xyyyy_1[i] * fe_0 + g_xxxx_xyyyyz_1[i] * pa_z[i];

        g_xxxxz_xyyyzz_0[i] = 2.0 * g_xxxx_xyyyz_1[i] * fe_0 + g_xxxx_xyyyzz_1[i] * pa_z[i];

        g_xxxxz_xyyzzz_0[i] = 3.0 * g_xxxx_xyyzz_1[i] * fe_0 + g_xxxx_xyyzzz_1[i] * pa_z[i];

        g_xxxxz_xyzzzz_0[i] = 4.0 * g_xxxx_xyzzz_1[i] * fe_0 + g_xxxx_xyzzzz_1[i] * pa_z[i];

        g_xxxxz_xzzzzz_0[i] = 5.0 * g_xxxx_xzzzz_1[i] * fe_0 + g_xxxx_xzzzzz_1[i] * pa_z[i];

        g_xxxxz_yyyyyy_0[i] = g_xxxx_yyyyyy_1[i] * pa_z[i];

        g_xxxxz_yyyyyz_0[i] = g_xxxx_yyyyy_1[i] * fe_0 + g_xxxx_yyyyyz_1[i] * pa_z[i];

        g_xxxxz_yyyyzz_0[i] = 2.0 * g_xxxx_yyyyz_1[i] * fe_0 + g_xxxx_yyyyzz_1[i] * pa_z[i];

        g_xxxxz_yyyzzz_0[i] = 3.0 * g_xxxx_yyyzz_1[i] * fe_0 + g_xxxx_yyyzzz_1[i] * pa_z[i];

        g_xxxxz_yyzzzz_0[i] = 4.0 * g_xxxx_yyzzz_1[i] * fe_0 + g_xxxx_yyzzzz_1[i] * pa_z[i];

        g_xxxxz_yzzzzz_0[i] = 5.0 * g_xxxx_yzzzz_1[i] * fe_0 + g_xxxx_yzzzzz_1[i] * pa_z[i];

        g_xxxxz_zzzzzz_0[i] = 6.0 * g_xxxx_zzzzz_1[i] * fe_0 + g_xxxx_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 84-112 components of targeted buffer : HI

    auto g_xxxyy_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 84);

    auto g_xxxyy_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 85);

    auto g_xxxyy_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 86);

    auto g_xxxyy_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 87);

    auto g_xxxyy_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 88);

    auto g_xxxyy_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 89);

    auto g_xxxyy_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 90);

    auto g_xxxyy_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 91);

    auto g_xxxyy_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 92);

    auto g_xxxyy_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 93);

    auto g_xxxyy_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 94);

    auto g_xxxyy_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 95);

    auto g_xxxyy_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 96);

    auto g_xxxyy_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 97);

    auto g_xxxyy_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 98);

    auto g_xxxyy_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 99);

    auto g_xxxyy_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 100);

    auto g_xxxyy_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 101);

    auto g_xxxyy_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 102);

    auto g_xxxyy_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 103);

    auto g_xxxyy_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 104);

    auto g_xxxyy_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 105);

    auto g_xxxyy_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 106);

    auto g_xxxyy_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 107);

    auto g_xxxyy_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 108);

    auto g_xxxyy_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 109);

    auto g_xxxyy_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 110);

    auto g_xxxyy_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 111);

    #pragma omp simd aligned(g_xxx_xxxxxx_0, g_xxx_xxxxxx_1, g_xxx_xxxxxz_0, g_xxx_xxxxxz_1, g_xxx_xxxxzz_0, g_xxx_xxxxzz_1, g_xxx_xxxzzz_0, g_xxx_xxxzzz_1, g_xxx_xxzzzz_0, g_xxx_xxzzzz_1, g_xxx_xzzzzz_0, g_xxx_xzzzzz_1, g_xxxy_xxxxxx_1, g_xxxy_xxxxxz_1, g_xxxy_xxxxzz_1, g_xxxy_xxxzzz_1, g_xxxy_xxzzzz_1, g_xxxy_xzzzzz_1, g_xxxyy_xxxxxx_0, g_xxxyy_xxxxxy_0, g_xxxyy_xxxxxz_0, g_xxxyy_xxxxyy_0, g_xxxyy_xxxxyz_0, g_xxxyy_xxxxzz_0, g_xxxyy_xxxyyy_0, g_xxxyy_xxxyyz_0, g_xxxyy_xxxyzz_0, g_xxxyy_xxxzzz_0, g_xxxyy_xxyyyy_0, g_xxxyy_xxyyyz_0, g_xxxyy_xxyyzz_0, g_xxxyy_xxyzzz_0, g_xxxyy_xxzzzz_0, g_xxxyy_xyyyyy_0, g_xxxyy_xyyyyz_0, g_xxxyy_xyyyzz_0, g_xxxyy_xyyzzz_0, g_xxxyy_xyzzzz_0, g_xxxyy_xzzzzz_0, g_xxxyy_yyyyyy_0, g_xxxyy_yyyyyz_0, g_xxxyy_yyyyzz_0, g_xxxyy_yyyzzz_0, g_xxxyy_yyzzzz_0, g_xxxyy_yzzzzz_0, g_xxxyy_zzzzzz_0, g_xxyy_xxxxxy_1, g_xxyy_xxxxy_1, g_xxyy_xxxxyy_1, g_xxyy_xxxxyz_1, g_xxyy_xxxyy_1, g_xxyy_xxxyyy_1, g_xxyy_xxxyyz_1, g_xxyy_xxxyz_1, g_xxyy_xxxyzz_1, g_xxyy_xxyyy_1, g_xxyy_xxyyyy_1, g_xxyy_xxyyyz_1, g_xxyy_xxyyz_1, g_xxyy_xxyyzz_1, g_xxyy_xxyzz_1, g_xxyy_xxyzzz_1, g_xxyy_xyyyy_1, g_xxyy_xyyyyy_1, g_xxyy_xyyyyz_1, g_xxyy_xyyyz_1, g_xxyy_xyyyzz_1, g_xxyy_xyyzz_1, g_xxyy_xyyzzz_1, g_xxyy_xyzzz_1, g_xxyy_xyzzzz_1, g_xxyy_yyyyy_1, g_xxyy_yyyyyy_1, g_xxyy_yyyyyz_1, g_xxyy_yyyyz_1, g_xxyy_yyyyzz_1, g_xxyy_yyyzz_1, g_xxyy_yyyzzz_1, g_xxyy_yyzzz_1, g_xxyy_yyzzzz_1, g_xxyy_yzzzz_1, g_xxyy_yzzzzz_1, g_xxyy_zzzzzz_1, g_xyy_xxxxxy_0, g_xyy_xxxxxy_1, g_xyy_xxxxyy_0, g_xyy_xxxxyy_1, g_xyy_xxxxyz_0, g_xyy_xxxxyz_1, g_xyy_xxxyyy_0, g_xyy_xxxyyy_1, g_xyy_xxxyyz_0, g_xyy_xxxyyz_1, g_xyy_xxxyzz_0, g_xyy_xxxyzz_1, g_xyy_xxyyyy_0, g_xyy_xxyyyy_1, g_xyy_xxyyyz_0, g_xyy_xxyyyz_1, g_xyy_xxyyzz_0, g_xyy_xxyyzz_1, g_xyy_xxyzzz_0, g_xyy_xxyzzz_1, g_xyy_xyyyyy_0, g_xyy_xyyyyy_1, g_xyy_xyyyyz_0, g_xyy_xyyyyz_1, g_xyy_xyyyzz_0, g_xyy_xyyyzz_1, g_xyy_xyyzzz_0, g_xyy_xyyzzz_1, g_xyy_xyzzzz_0, g_xyy_xyzzzz_1, g_xyy_yyyyyy_0, g_xyy_yyyyyy_1, g_xyy_yyyyyz_0, g_xyy_yyyyyz_1, g_xyy_yyyyzz_0, g_xyy_yyyyzz_1, g_xyy_yyyzzz_0, g_xyy_yyyzzz_1, g_xyy_yyzzzz_0, g_xyy_yyzzzz_1, g_xyy_yzzzzz_0, g_xyy_yzzzzz_1, g_xyy_zzzzzz_0, g_xyy_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyy_xxxxxx_0[i] = g_xxx_xxxxxx_0[i] * fbe_0 - g_xxx_xxxxxx_1[i] * fz_be_0 + g_xxxy_xxxxxx_1[i] * pa_y[i];

        g_xxxyy_xxxxxy_0[i] = 2.0 * g_xyy_xxxxxy_0[i] * fbe_0 - 2.0 * g_xyy_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxyy_xxxxy_1[i] * fe_0 + g_xxyy_xxxxxy_1[i] * pa_x[i];

        g_xxxyy_xxxxxz_0[i] = g_xxx_xxxxxz_0[i] * fbe_0 - g_xxx_xxxxxz_1[i] * fz_be_0 + g_xxxy_xxxxxz_1[i] * pa_y[i];

        g_xxxyy_xxxxyy_0[i] = 2.0 * g_xyy_xxxxyy_0[i] * fbe_0 - 2.0 * g_xyy_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxyy_xxxyy_1[i] * fe_0 + g_xxyy_xxxxyy_1[i] * pa_x[i];

        g_xxxyy_xxxxyz_0[i] = 2.0 * g_xyy_xxxxyz_0[i] * fbe_0 - 2.0 * g_xyy_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxyy_xxxyz_1[i] * fe_0 + g_xxyy_xxxxyz_1[i] * pa_x[i];

        g_xxxyy_xxxxzz_0[i] = g_xxx_xxxxzz_0[i] * fbe_0 - g_xxx_xxxxzz_1[i] * fz_be_0 + g_xxxy_xxxxzz_1[i] * pa_y[i];

        g_xxxyy_xxxyyy_0[i] = 2.0 * g_xyy_xxxyyy_0[i] * fbe_0 - 2.0 * g_xyy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxyy_xxyyy_1[i] * fe_0 + g_xxyy_xxxyyy_1[i] * pa_x[i];

        g_xxxyy_xxxyyz_0[i] = 2.0 * g_xyy_xxxyyz_0[i] * fbe_0 - 2.0 * g_xyy_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxyy_xxyyz_1[i] * fe_0 + g_xxyy_xxxyyz_1[i] * pa_x[i];

        g_xxxyy_xxxyzz_0[i] = 2.0 * g_xyy_xxxyzz_0[i] * fbe_0 - 2.0 * g_xyy_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxyy_xxyzz_1[i] * fe_0 + g_xxyy_xxxyzz_1[i] * pa_x[i];

        g_xxxyy_xxxzzz_0[i] = g_xxx_xxxzzz_0[i] * fbe_0 - g_xxx_xxxzzz_1[i] * fz_be_0 + g_xxxy_xxxzzz_1[i] * pa_y[i];

        g_xxxyy_xxyyyy_0[i] = 2.0 * g_xyy_xxyyyy_0[i] * fbe_0 - 2.0 * g_xyy_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxyy_xyyyy_1[i] * fe_0 + g_xxyy_xxyyyy_1[i] * pa_x[i];

        g_xxxyy_xxyyyz_0[i] = 2.0 * g_xyy_xxyyyz_0[i] * fbe_0 - 2.0 * g_xyy_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxyy_xyyyz_1[i] * fe_0 + g_xxyy_xxyyyz_1[i] * pa_x[i];

        g_xxxyy_xxyyzz_0[i] = 2.0 * g_xyy_xxyyzz_0[i] * fbe_0 - 2.0 * g_xyy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxyy_xyyzz_1[i] * fe_0 + g_xxyy_xxyyzz_1[i] * pa_x[i];

        g_xxxyy_xxyzzz_0[i] = 2.0 * g_xyy_xxyzzz_0[i] * fbe_0 - 2.0 * g_xyy_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxyy_xyzzz_1[i] * fe_0 + g_xxyy_xxyzzz_1[i] * pa_x[i];

        g_xxxyy_xxzzzz_0[i] = g_xxx_xxzzzz_0[i] * fbe_0 - g_xxx_xxzzzz_1[i] * fz_be_0 + g_xxxy_xxzzzz_1[i] * pa_y[i];

        g_xxxyy_xyyyyy_0[i] = 2.0 * g_xyy_xyyyyy_0[i] * fbe_0 - 2.0 * g_xyy_xyyyyy_1[i] * fz_be_0 + g_xxyy_yyyyy_1[i] * fe_0 + g_xxyy_xyyyyy_1[i] * pa_x[i];

        g_xxxyy_xyyyyz_0[i] = 2.0 * g_xyy_xyyyyz_0[i] * fbe_0 - 2.0 * g_xyy_xyyyyz_1[i] * fz_be_0 + g_xxyy_yyyyz_1[i] * fe_0 + g_xxyy_xyyyyz_1[i] * pa_x[i];

        g_xxxyy_xyyyzz_0[i] = 2.0 * g_xyy_xyyyzz_0[i] * fbe_0 - 2.0 * g_xyy_xyyyzz_1[i] * fz_be_0 + g_xxyy_yyyzz_1[i] * fe_0 + g_xxyy_xyyyzz_1[i] * pa_x[i];

        g_xxxyy_xyyzzz_0[i] = 2.0 * g_xyy_xyyzzz_0[i] * fbe_0 - 2.0 * g_xyy_xyyzzz_1[i] * fz_be_0 + g_xxyy_yyzzz_1[i] * fe_0 + g_xxyy_xyyzzz_1[i] * pa_x[i];

        g_xxxyy_xyzzzz_0[i] = 2.0 * g_xyy_xyzzzz_0[i] * fbe_0 - 2.0 * g_xyy_xyzzzz_1[i] * fz_be_0 + g_xxyy_yzzzz_1[i] * fe_0 + g_xxyy_xyzzzz_1[i] * pa_x[i];

        g_xxxyy_xzzzzz_0[i] = g_xxx_xzzzzz_0[i] * fbe_0 - g_xxx_xzzzzz_1[i] * fz_be_0 + g_xxxy_xzzzzz_1[i] * pa_y[i];

        g_xxxyy_yyyyyy_0[i] = 2.0 * g_xyy_yyyyyy_0[i] * fbe_0 - 2.0 * g_xyy_yyyyyy_1[i] * fz_be_0 + g_xxyy_yyyyyy_1[i] * pa_x[i];

        g_xxxyy_yyyyyz_0[i] = 2.0 * g_xyy_yyyyyz_0[i] * fbe_0 - 2.0 * g_xyy_yyyyyz_1[i] * fz_be_0 + g_xxyy_yyyyyz_1[i] * pa_x[i];

        g_xxxyy_yyyyzz_0[i] = 2.0 * g_xyy_yyyyzz_0[i] * fbe_0 - 2.0 * g_xyy_yyyyzz_1[i] * fz_be_0 + g_xxyy_yyyyzz_1[i] * pa_x[i];

        g_xxxyy_yyyzzz_0[i] = 2.0 * g_xyy_yyyzzz_0[i] * fbe_0 - 2.0 * g_xyy_yyyzzz_1[i] * fz_be_0 + g_xxyy_yyyzzz_1[i] * pa_x[i];

        g_xxxyy_yyzzzz_0[i] = 2.0 * g_xyy_yyzzzz_0[i] * fbe_0 - 2.0 * g_xyy_yyzzzz_1[i] * fz_be_0 + g_xxyy_yyzzzz_1[i] * pa_x[i];

        g_xxxyy_yzzzzz_0[i] = 2.0 * g_xyy_yzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_yzzzzz_1[i] * fz_be_0 + g_xxyy_yzzzzz_1[i] * pa_x[i];

        g_xxxyy_zzzzzz_0[i] = 2.0 * g_xyy_zzzzzz_0[i] * fbe_0 - 2.0 * g_xyy_zzzzzz_1[i] * fz_be_0 + g_xxyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 112-140 components of targeted buffer : HI

    auto g_xxxyz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 112);

    auto g_xxxyz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 113);

    auto g_xxxyz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 114);

    auto g_xxxyz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 115);

    auto g_xxxyz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 116);

    auto g_xxxyz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 117);

    auto g_xxxyz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 118);

    auto g_xxxyz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 119);

    auto g_xxxyz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 120);

    auto g_xxxyz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 121);

    auto g_xxxyz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 122);

    auto g_xxxyz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 123);

    auto g_xxxyz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 124);

    auto g_xxxyz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 125);

    auto g_xxxyz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 126);

    auto g_xxxyz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 127);

    auto g_xxxyz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 128);

    auto g_xxxyz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 129);

    auto g_xxxyz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 130);

    auto g_xxxyz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 131);

    auto g_xxxyz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 132);

    auto g_xxxyz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 133);

    auto g_xxxyz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 134);

    auto g_xxxyz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 135);

    auto g_xxxyz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 136);

    auto g_xxxyz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 137);

    auto g_xxxyz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 138);

    auto g_xxxyz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 139);

    #pragma omp simd aligned(g_xxxy_xxxxxy_1, g_xxxy_xxxxyy_1, g_xxxy_xxxyyy_1, g_xxxy_xxyyyy_1, g_xxxy_xyyyyy_1, g_xxxy_yyyyyy_1, g_xxxyz_xxxxxx_0, g_xxxyz_xxxxxy_0, g_xxxyz_xxxxxz_0, g_xxxyz_xxxxyy_0, g_xxxyz_xxxxyz_0, g_xxxyz_xxxxzz_0, g_xxxyz_xxxyyy_0, g_xxxyz_xxxyyz_0, g_xxxyz_xxxyzz_0, g_xxxyz_xxxzzz_0, g_xxxyz_xxyyyy_0, g_xxxyz_xxyyyz_0, g_xxxyz_xxyyzz_0, g_xxxyz_xxyzzz_0, g_xxxyz_xxzzzz_0, g_xxxyz_xyyyyy_0, g_xxxyz_xyyyyz_0, g_xxxyz_xyyyzz_0, g_xxxyz_xyyzzz_0, g_xxxyz_xyzzzz_0, g_xxxyz_xzzzzz_0, g_xxxyz_yyyyyy_0, g_xxxyz_yyyyyz_0, g_xxxyz_yyyyzz_0, g_xxxyz_yyyzzz_0, g_xxxyz_yyzzzz_0, g_xxxyz_yzzzzz_0, g_xxxyz_zzzzzz_0, g_xxxz_xxxxxx_1, g_xxxz_xxxxxz_1, g_xxxz_xxxxyz_1, g_xxxz_xxxxz_1, g_xxxz_xxxxzz_1, g_xxxz_xxxyyz_1, g_xxxz_xxxyz_1, g_xxxz_xxxyzz_1, g_xxxz_xxxzz_1, g_xxxz_xxxzzz_1, g_xxxz_xxyyyz_1, g_xxxz_xxyyz_1, g_xxxz_xxyyzz_1, g_xxxz_xxyzz_1, g_xxxz_xxyzzz_1, g_xxxz_xxzzz_1, g_xxxz_xxzzzz_1, g_xxxz_xyyyyz_1, g_xxxz_xyyyz_1, g_xxxz_xyyyzz_1, g_xxxz_xyyzz_1, g_xxxz_xyyzzz_1, g_xxxz_xyzzz_1, g_xxxz_xyzzzz_1, g_xxxz_xzzzz_1, g_xxxz_xzzzzz_1, g_xxxz_yyyyyz_1, g_xxxz_yyyyz_1, g_xxxz_yyyyzz_1, g_xxxz_yyyzz_1, g_xxxz_yyyzzz_1, g_xxxz_yyzzz_1, g_xxxz_yyzzzz_1, g_xxxz_yzzzz_1, g_xxxz_yzzzzz_1, g_xxxz_zzzzz_1, g_xxxz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyz_xxxxxx_0[i] = g_xxxz_xxxxxx_1[i] * pa_y[i];

        g_xxxyz_xxxxxy_0[i] = g_xxxy_xxxxxy_1[i] * pa_z[i];

        g_xxxyz_xxxxxz_0[i] = g_xxxz_xxxxxz_1[i] * pa_y[i];

        g_xxxyz_xxxxyy_0[i] = g_xxxy_xxxxyy_1[i] * pa_z[i];

        g_xxxyz_xxxxyz_0[i] = g_xxxz_xxxxz_1[i] * fe_0 + g_xxxz_xxxxyz_1[i] * pa_y[i];

        g_xxxyz_xxxxzz_0[i] = g_xxxz_xxxxzz_1[i] * pa_y[i];

        g_xxxyz_xxxyyy_0[i] = g_xxxy_xxxyyy_1[i] * pa_z[i];

        g_xxxyz_xxxyyz_0[i] = 2.0 * g_xxxz_xxxyz_1[i] * fe_0 + g_xxxz_xxxyyz_1[i] * pa_y[i];

        g_xxxyz_xxxyzz_0[i] = g_xxxz_xxxzz_1[i] * fe_0 + g_xxxz_xxxyzz_1[i] * pa_y[i];

        g_xxxyz_xxxzzz_0[i] = g_xxxz_xxxzzz_1[i] * pa_y[i];

        g_xxxyz_xxyyyy_0[i] = g_xxxy_xxyyyy_1[i] * pa_z[i];

        g_xxxyz_xxyyyz_0[i] = 3.0 * g_xxxz_xxyyz_1[i] * fe_0 + g_xxxz_xxyyyz_1[i] * pa_y[i];

        g_xxxyz_xxyyzz_0[i] = 2.0 * g_xxxz_xxyzz_1[i] * fe_0 + g_xxxz_xxyyzz_1[i] * pa_y[i];

        g_xxxyz_xxyzzz_0[i] = g_xxxz_xxzzz_1[i] * fe_0 + g_xxxz_xxyzzz_1[i] * pa_y[i];

        g_xxxyz_xxzzzz_0[i] = g_xxxz_xxzzzz_1[i] * pa_y[i];

        g_xxxyz_xyyyyy_0[i] = g_xxxy_xyyyyy_1[i] * pa_z[i];

        g_xxxyz_xyyyyz_0[i] = 4.0 * g_xxxz_xyyyz_1[i] * fe_0 + g_xxxz_xyyyyz_1[i] * pa_y[i];

        g_xxxyz_xyyyzz_0[i] = 3.0 * g_xxxz_xyyzz_1[i] * fe_0 + g_xxxz_xyyyzz_1[i] * pa_y[i];

        g_xxxyz_xyyzzz_0[i] = 2.0 * g_xxxz_xyzzz_1[i] * fe_0 + g_xxxz_xyyzzz_1[i] * pa_y[i];

        g_xxxyz_xyzzzz_0[i] = g_xxxz_xzzzz_1[i] * fe_0 + g_xxxz_xyzzzz_1[i] * pa_y[i];

        g_xxxyz_xzzzzz_0[i] = g_xxxz_xzzzzz_1[i] * pa_y[i];

        g_xxxyz_yyyyyy_0[i] = g_xxxy_yyyyyy_1[i] * pa_z[i];

        g_xxxyz_yyyyyz_0[i] = 5.0 * g_xxxz_yyyyz_1[i] * fe_0 + g_xxxz_yyyyyz_1[i] * pa_y[i];

        g_xxxyz_yyyyzz_0[i] = 4.0 * g_xxxz_yyyzz_1[i] * fe_0 + g_xxxz_yyyyzz_1[i] * pa_y[i];

        g_xxxyz_yyyzzz_0[i] = 3.0 * g_xxxz_yyzzz_1[i] * fe_0 + g_xxxz_yyyzzz_1[i] * pa_y[i];

        g_xxxyz_yyzzzz_0[i] = 2.0 * g_xxxz_yzzzz_1[i] * fe_0 + g_xxxz_yyzzzz_1[i] * pa_y[i];

        g_xxxyz_yzzzzz_0[i] = g_xxxz_zzzzz_1[i] * fe_0 + g_xxxz_yzzzzz_1[i] * pa_y[i];

        g_xxxyz_zzzzzz_0[i] = g_xxxz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : HI

    auto g_xxxzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 140);

    auto g_xxxzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 141);

    auto g_xxxzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 142);

    auto g_xxxzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 143);

    auto g_xxxzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 144);

    auto g_xxxzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 145);

    auto g_xxxzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 146);

    auto g_xxxzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 147);

    auto g_xxxzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 148);

    auto g_xxxzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 149);

    auto g_xxxzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 150);

    auto g_xxxzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 151);

    auto g_xxxzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 152);

    auto g_xxxzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 153);

    auto g_xxxzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 154);

    auto g_xxxzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 155);

    auto g_xxxzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 156);

    auto g_xxxzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 157);

    auto g_xxxzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 158);

    auto g_xxxzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 159);

    auto g_xxxzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 160);

    auto g_xxxzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 161);

    auto g_xxxzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 162);

    auto g_xxxzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 163);

    auto g_xxxzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 164);

    auto g_xxxzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 165);

    auto g_xxxzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 166);

    auto g_xxxzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 167);

    #pragma omp simd aligned(g_xxx_xxxxxx_0, g_xxx_xxxxxx_1, g_xxx_xxxxxy_0, g_xxx_xxxxxy_1, g_xxx_xxxxyy_0, g_xxx_xxxxyy_1, g_xxx_xxxyyy_0, g_xxx_xxxyyy_1, g_xxx_xxyyyy_0, g_xxx_xxyyyy_1, g_xxx_xyyyyy_0, g_xxx_xyyyyy_1, g_xxxz_xxxxxx_1, g_xxxz_xxxxxy_1, g_xxxz_xxxxyy_1, g_xxxz_xxxyyy_1, g_xxxz_xxyyyy_1, g_xxxz_xyyyyy_1, g_xxxzz_xxxxxx_0, g_xxxzz_xxxxxy_0, g_xxxzz_xxxxxz_0, g_xxxzz_xxxxyy_0, g_xxxzz_xxxxyz_0, g_xxxzz_xxxxzz_0, g_xxxzz_xxxyyy_0, g_xxxzz_xxxyyz_0, g_xxxzz_xxxyzz_0, g_xxxzz_xxxzzz_0, g_xxxzz_xxyyyy_0, g_xxxzz_xxyyyz_0, g_xxxzz_xxyyzz_0, g_xxxzz_xxyzzz_0, g_xxxzz_xxzzzz_0, g_xxxzz_xyyyyy_0, g_xxxzz_xyyyyz_0, g_xxxzz_xyyyzz_0, g_xxxzz_xyyzzz_0, g_xxxzz_xyzzzz_0, g_xxxzz_xzzzzz_0, g_xxxzz_yyyyyy_0, g_xxxzz_yyyyyz_0, g_xxxzz_yyyyzz_0, g_xxxzz_yyyzzz_0, g_xxxzz_yyzzzz_0, g_xxxzz_yzzzzz_0, g_xxxzz_zzzzzz_0, g_xxzz_xxxxxz_1, g_xxzz_xxxxyz_1, g_xxzz_xxxxz_1, g_xxzz_xxxxzz_1, g_xxzz_xxxyyz_1, g_xxzz_xxxyz_1, g_xxzz_xxxyzz_1, g_xxzz_xxxzz_1, g_xxzz_xxxzzz_1, g_xxzz_xxyyyz_1, g_xxzz_xxyyz_1, g_xxzz_xxyyzz_1, g_xxzz_xxyzz_1, g_xxzz_xxyzzz_1, g_xxzz_xxzzz_1, g_xxzz_xxzzzz_1, g_xxzz_xyyyyz_1, g_xxzz_xyyyz_1, g_xxzz_xyyyzz_1, g_xxzz_xyyzz_1, g_xxzz_xyyzzz_1, g_xxzz_xyzzz_1, g_xxzz_xyzzzz_1, g_xxzz_xzzzz_1, g_xxzz_xzzzzz_1, g_xxzz_yyyyyy_1, g_xxzz_yyyyyz_1, g_xxzz_yyyyz_1, g_xxzz_yyyyzz_1, g_xxzz_yyyzz_1, g_xxzz_yyyzzz_1, g_xxzz_yyzzz_1, g_xxzz_yyzzzz_1, g_xxzz_yzzzz_1, g_xxzz_yzzzzz_1, g_xxzz_zzzzz_1, g_xxzz_zzzzzz_1, g_xzz_xxxxxz_0, g_xzz_xxxxxz_1, g_xzz_xxxxyz_0, g_xzz_xxxxyz_1, g_xzz_xxxxzz_0, g_xzz_xxxxzz_1, g_xzz_xxxyyz_0, g_xzz_xxxyyz_1, g_xzz_xxxyzz_0, g_xzz_xxxyzz_1, g_xzz_xxxzzz_0, g_xzz_xxxzzz_1, g_xzz_xxyyyz_0, g_xzz_xxyyyz_1, g_xzz_xxyyzz_0, g_xzz_xxyyzz_1, g_xzz_xxyzzz_0, g_xzz_xxyzzz_1, g_xzz_xxzzzz_0, g_xzz_xxzzzz_1, g_xzz_xyyyyz_0, g_xzz_xyyyyz_1, g_xzz_xyyyzz_0, g_xzz_xyyyzz_1, g_xzz_xyyzzz_0, g_xzz_xyyzzz_1, g_xzz_xyzzzz_0, g_xzz_xyzzzz_1, g_xzz_xzzzzz_0, g_xzz_xzzzzz_1, g_xzz_yyyyyy_0, g_xzz_yyyyyy_1, g_xzz_yyyyyz_0, g_xzz_yyyyyz_1, g_xzz_yyyyzz_0, g_xzz_yyyyzz_1, g_xzz_yyyzzz_0, g_xzz_yyyzzz_1, g_xzz_yyzzzz_0, g_xzz_yyzzzz_1, g_xzz_yzzzzz_0, g_xzz_yzzzzz_1, g_xzz_zzzzzz_0, g_xzz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzz_xxxxxx_0[i] = g_xxx_xxxxxx_0[i] * fbe_0 - g_xxx_xxxxxx_1[i] * fz_be_0 + g_xxxz_xxxxxx_1[i] * pa_z[i];

        g_xxxzz_xxxxxy_0[i] = g_xxx_xxxxxy_0[i] * fbe_0 - g_xxx_xxxxxy_1[i] * fz_be_0 + g_xxxz_xxxxxy_1[i] * pa_z[i];

        g_xxxzz_xxxxxz_0[i] = 2.0 * g_xzz_xxxxxz_0[i] * fbe_0 - 2.0 * g_xzz_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxzz_xxxxz_1[i] * fe_0 + g_xxzz_xxxxxz_1[i] * pa_x[i];

        g_xxxzz_xxxxyy_0[i] = g_xxx_xxxxyy_0[i] * fbe_0 - g_xxx_xxxxyy_1[i] * fz_be_0 + g_xxxz_xxxxyy_1[i] * pa_z[i];

        g_xxxzz_xxxxyz_0[i] = 2.0 * g_xzz_xxxxyz_0[i] * fbe_0 - 2.0 * g_xzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxzz_xxxyz_1[i] * fe_0 + g_xxzz_xxxxyz_1[i] * pa_x[i];

        g_xxxzz_xxxxzz_0[i] = 2.0 * g_xzz_xxxxzz_0[i] * fbe_0 - 2.0 * g_xzz_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxzz_xxxzz_1[i] * fe_0 + g_xxzz_xxxxzz_1[i] * pa_x[i];

        g_xxxzz_xxxyyy_0[i] = g_xxx_xxxyyy_0[i] * fbe_0 - g_xxx_xxxyyy_1[i] * fz_be_0 + g_xxxz_xxxyyy_1[i] * pa_z[i];

        g_xxxzz_xxxyyz_0[i] = 2.0 * g_xzz_xxxyyz_0[i] * fbe_0 - 2.0 * g_xzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxzz_xxyyz_1[i] * fe_0 + g_xxzz_xxxyyz_1[i] * pa_x[i];

        g_xxxzz_xxxyzz_0[i] = 2.0 * g_xzz_xxxyzz_0[i] * fbe_0 - 2.0 * g_xzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxzz_xxyzz_1[i] * fe_0 + g_xxzz_xxxyzz_1[i] * pa_x[i];

        g_xxxzz_xxxzzz_0[i] = 2.0 * g_xzz_xxxzzz_0[i] * fbe_0 - 2.0 * g_xzz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxzz_xxzzz_1[i] * fe_0 + g_xxzz_xxxzzz_1[i] * pa_x[i];

        g_xxxzz_xxyyyy_0[i] = g_xxx_xxyyyy_0[i] * fbe_0 - g_xxx_xxyyyy_1[i] * fz_be_0 + g_xxxz_xxyyyy_1[i] * pa_z[i];

        g_xxxzz_xxyyyz_0[i] = 2.0 * g_xzz_xxyyyz_0[i] * fbe_0 - 2.0 * g_xzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxzz_xyyyz_1[i] * fe_0 + g_xxzz_xxyyyz_1[i] * pa_x[i];

        g_xxxzz_xxyyzz_0[i] = 2.0 * g_xzz_xxyyzz_0[i] * fbe_0 - 2.0 * g_xzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxzz_xyyzz_1[i] * fe_0 + g_xxzz_xxyyzz_1[i] * pa_x[i];

        g_xxxzz_xxyzzz_0[i] = 2.0 * g_xzz_xxyzzz_0[i] * fbe_0 - 2.0 * g_xzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_xyzzz_1[i] * fe_0 + g_xxzz_xxyzzz_1[i] * pa_x[i];

        g_xxxzz_xxzzzz_0[i] = 2.0 * g_xzz_xxzzzz_0[i] * fbe_0 - 2.0 * g_xzz_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxzz_xzzzz_1[i] * fe_0 + g_xxzz_xxzzzz_1[i] * pa_x[i];

        g_xxxzz_xyyyyy_0[i] = g_xxx_xyyyyy_0[i] * fbe_0 - g_xxx_xyyyyy_1[i] * fz_be_0 + g_xxxz_xyyyyy_1[i] * pa_z[i];

        g_xxxzz_xyyyyz_0[i] = 2.0 * g_xzz_xyyyyz_0[i] * fbe_0 - 2.0 * g_xzz_xyyyyz_1[i] * fz_be_0 + g_xxzz_yyyyz_1[i] * fe_0 + g_xxzz_xyyyyz_1[i] * pa_x[i];

        g_xxxzz_xyyyzz_0[i] = 2.0 * g_xzz_xyyyzz_0[i] * fbe_0 - 2.0 * g_xzz_xyyyzz_1[i] * fz_be_0 + g_xxzz_yyyzz_1[i] * fe_0 + g_xxzz_xyyyzz_1[i] * pa_x[i];

        g_xxxzz_xyyzzz_0[i] = 2.0 * g_xzz_xyyzzz_0[i] * fbe_0 - 2.0 * g_xzz_xyyzzz_1[i] * fz_be_0 + g_xxzz_yyzzz_1[i] * fe_0 + g_xxzz_xyyzzz_1[i] * pa_x[i];

        g_xxxzz_xyzzzz_0[i] = 2.0 * g_xzz_xyzzzz_0[i] * fbe_0 - 2.0 * g_xzz_xyzzzz_1[i] * fz_be_0 + g_xxzz_yzzzz_1[i] * fe_0 + g_xxzz_xyzzzz_1[i] * pa_x[i];

        g_xxxzz_xzzzzz_0[i] = 2.0 * g_xzz_xzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_xzzzzz_1[i] * fz_be_0 + g_xxzz_zzzzz_1[i] * fe_0 + g_xxzz_xzzzzz_1[i] * pa_x[i];

        g_xxxzz_yyyyyy_0[i] = 2.0 * g_xzz_yyyyyy_0[i] * fbe_0 - 2.0 * g_xzz_yyyyyy_1[i] * fz_be_0 + g_xxzz_yyyyyy_1[i] * pa_x[i];

        g_xxxzz_yyyyyz_0[i] = 2.0 * g_xzz_yyyyyz_0[i] * fbe_0 - 2.0 * g_xzz_yyyyyz_1[i] * fz_be_0 + g_xxzz_yyyyyz_1[i] * pa_x[i];

        g_xxxzz_yyyyzz_0[i] = 2.0 * g_xzz_yyyyzz_0[i] * fbe_0 - 2.0 * g_xzz_yyyyzz_1[i] * fz_be_0 + g_xxzz_yyyyzz_1[i] * pa_x[i];

        g_xxxzz_yyyzzz_0[i] = 2.0 * g_xzz_yyyzzz_0[i] * fbe_0 - 2.0 * g_xzz_yyyzzz_1[i] * fz_be_0 + g_xxzz_yyyzzz_1[i] * pa_x[i];

        g_xxxzz_yyzzzz_0[i] = 2.0 * g_xzz_yyzzzz_0[i] * fbe_0 - 2.0 * g_xzz_yyzzzz_1[i] * fz_be_0 + g_xxzz_yyzzzz_1[i] * pa_x[i];

        g_xxxzz_yzzzzz_0[i] = 2.0 * g_xzz_yzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_yzzzzz_1[i] * fz_be_0 + g_xxzz_yzzzzz_1[i] * pa_x[i];

        g_xxxzz_zzzzzz_0[i] = 2.0 * g_xzz_zzzzzz_0[i] * fbe_0 - 2.0 * g_xzz_zzzzzz_1[i] * fz_be_0 + g_xxzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : HI

    auto g_xxyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 168);

    auto g_xxyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 169);

    auto g_xxyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 170);

    auto g_xxyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 171);

    auto g_xxyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 172);

    auto g_xxyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 173);

    auto g_xxyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 174);

    auto g_xxyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 175);

    auto g_xxyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 176);

    auto g_xxyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 177);

    auto g_xxyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 178);

    auto g_xxyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 179);

    auto g_xxyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 180);

    auto g_xxyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 181);

    auto g_xxyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 182);

    auto g_xxyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 183);

    auto g_xxyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 184);

    auto g_xxyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 185);

    auto g_xxyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 186);

    auto g_xxyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 187);

    auto g_xxyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 188);

    auto g_xxyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 189);

    auto g_xxyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 190);

    auto g_xxyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 191);

    auto g_xxyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 192);

    auto g_xxyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 193);

    auto g_xxyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 194);

    auto g_xxyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 195);

    #pragma omp simd aligned(g_xxy_xxxxxx_0, g_xxy_xxxxxx_1, g_xxy_xxxxxz_0, g_xxy_xxxxxz_1, g_xxy_xxxxzz_0, g_xxy_xxxxzz_1, g_xxy_xxxzzz_0, g_xxy_xxxzzz_1, g_xxy_xxzzzz_0, g_xxy_xxzzzz_1, g_xxy_xzzzzz_0, g_xxy_xzzzzz_1, g_xxyy_xxxxxx_1, g_xxyy_xxxxxz_1, g_xxyy_xxxxzz_1, g_xxyy_xxxzzz_1, g_xxyy_xxzzzz_1, g_xxyy_xzzzzz_1, g_xxyyy_xxxxxx_0, g_xxyyy_xxxxxy_0, g_xxyyy_xxxxxz_0, g_xxyyy_xxxxyy_0, g_xxyyy_xxxxyz_0, g_xxyyy_xxxxzz_0, g_xxyyy_xxxyyy_0, g_xxyyy_xxxyyz_0, g_xxyyy_xxxyzz_0, g_xxyyy_xxxzzz_0, g_xxyyy_xxyyyy_0, g_xxyyy_xxyyyz_0, g_xxyyy_xxyyzz_0, g_xxyyy_xxyzzz_0, g_xxyyy_xxzzzz_0, g_xxyyy_xyyyyy_0, g_xxyyy_xyyyyz_0, g_xxyyy_xyyyzz_0, g_xxyyy_xyyzzz_0, g_xxyyy_xyzzzz_0, g_xxyyy_xzzzzz_0, g_xxyyy_yyyyyy_0, g_xxyyy_yyyyyz_0, g_xxyyy_yyyyzz_0, g_xxyyy_yyyzzz_0, g_xxyyy_yyzzzz_0, g_xxyyy_yzzzzz_0, g_xxyyy_zzzzzz_0, g_xyyy_xxxxxy_1, g_xyyy_xxxxy_1, g_xyyy_xxxxyy_1, g_xyyy_xxxxyz_1, g_xyyy_xxxyy_1, g_xyyy_xxxyyy_1, g_xyyy_xxxyyz_1, g_xyyy_xxxyz_1, g_xyyy_xxxyzz_1, g_xyyy_xxyyy_1, g_xyyy_xxyyyy_1, g_xyyy_xxyyyz_1, g_xyyy_xxyyz_1, g_xyyy_xxyyzz_1, g_xyyy_xxyzz_1, g_xyyy_xxyzzz_1, g_xyyy_xyyyy_1, g_xyyy_xyyyyy_1, g_xyyy_xyyyyz_1, g_xyyy_xyyyz_1, g_xyyy_xyyyzz_1, g_xyyy_xyyzz_1, g_xyyy_xyyzzz_1, g_xyyy_xyzzz_1, g_xyyy_xyzzzz_1, g_xyyy_yyyyy_1, g_xyyy_yyyyyy_1, g_xyyy_yyyyyz_1, g_xyyy_yyyyz_1, g_xyyy_yyyyzz_1, g_xyyy_yyyzz_1, g_xyyy_yyyzzz_1, g_xyyy_yyzzz_1, g_xyyy_yyzzzz_1, g_xyyy_yzzzz_1, g_xyyy_yzzzzz_1, g_xyyy_zzzzzz_1, g_yyy_xxxxxy_0, g_yyy_xxxxxy_1, g_yyy_xxxxyy_0, g_yyy_xxxxyy_1, g_yyy_xxxxyz_0, g_yyy_xxxxyz_1, g_yyy_xxxyyy_0, g_yyy_xxxyyy_1, g_yyy_xxxyyz_0, g_yyy_xxxyyz_1, g_yyy_xxxyzz_0, g_yyy_xxxyzz_1, g_yyy_xxyyyy_0, g_yyy_xxyyyy_1, g_yyy_xxyyyz_0, g_yyy_xxyyyz_1, g_yyy_xxyyzz_0, g_yyy_xxyyzz_1, g_yyy_xxyzzz_0, g_yyy_xxyzzz_1, g_yyy_xyyyyy_0, g_yyy_xyyyyy_1, g_yyy_xyyyyz_0, g_yyy_xyyyyz_1, g_yyy_xyyyzz_0, g_yyy_xyyyzz_1, g_yyy_xyyzzz_0, g_yyy_xyyzzz_1, g_yyy_xyzzzz_0, g_yyy_xyzzzz_1, g_yyy_yyyyyy_0, g_yyy_yyyyyy_1, g_yyy_yyyyyz_0, g_yyy_yyyyyz_1, g_yyy_yyyyzz_0, g_yyy_yyyyzz_1, g_yyy_yyyzzz_0, g_yyy_yyyzzz_1, g_yyy_yyzzzz_0, g_yyy_yyzzzz_1, g_yyy_yzzzzz_0, g_yyy_yzzzzz_1, g_yyy_zzzzzz_0, g_yyy_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyy_xxxxxx_0[i] = 2.0 * g_xxy_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxy_xxxxxx_1[i] * fz_be_0 + g_xxyy_xxxxxx_1[i] * pa_y[i];

        g_xxyyy_xxxxxy_0[i] = g_yyy_xxxxxy_0[i] * fbe_0 - g_yyy_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xyyy_xxxxy_1[i] * fe_0 + g_xyyy_xxxxxy_1[i] * pa_x[i];

        g_xxyyy_xxxxxz_0[i] = 2.0 * g_xxy_xxxxxz_0[i] * fbe_0 - 2.0 * g_xxy_xxxxxz_1[i] * fz_be_0 + g_xxyy_xxxxxz_1[i] * pa_y[i];

        g_xxyyy_xxxxyy_0[i] = g_yyy_xxxxyy_0[i] * fbe_0 - g_yyy_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xyyy_xxxyy_1[i] * fe_0 + g_xyyy_xxxxyy_1[i] * pa_x[i];

        g_xxyyy_xxxxyz_0[i] = g_yyy_xxxxyz_0[i] * fbe_0 - g_yyy_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyyy_xxxyz_1[i] * fe_0 + g_xyyy_xxxxyz_1[i] * pa_x[i];

        g_xxyyy_xxxxzz_0[i] = 2.0 * g_xxy_xxxxzz_0[i] * fbe_0 - 2.0 * g_xxy_xxxxzz_1[i] * fz_be_0 + g_xxyy_xxxxzz_1[i] * pa_y[i];

        g_xxyyy_xxxyyy_0[i] = g_yyy_xxxyyy_0[i] * fbe_0 - g_yyy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xyyy_xxyyy_1[i] * fe_0 + g_xyyy_xxxyyy_1[i] * pa_x[i];

        g_xxyyy_xxxyyz_0[i] = g_yyy_xxxyyz_0[i] * fbe_0 - g_yyy_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyyy_xxyyz_1[i] * fe_0 + g_xyyy_xxxyyz_1[i] * pa_x[i];

        g_xxyyy_xxxyzz_0[i] = g_yyy_xxxyzz_0[i] * fbe_0 - g_yyy_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyyy_xxyzz_1[i] * fe_0 + g_xyyy_xxxyzz_1[i] * pa_x[i];

        g_xxyyy_xxxzzz_0[i] = 2.0 * g_xxy_xxxzzz_0[i] * fbe_0 - 2.0 * g_xxy_xxxzzz_1[i] * fz_be_0 + g_xxyy_xxxzzz_1[i] * pa_y[i];

        g_xxyyy_xxyyyy_0[i] = g_yyy_xxyyyy_0[i] * fbe_0 - g_yyy_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xyyy_xyyyy_1[i] * fe_0 + g_xyyy_xxyyyy_1[i] * pa_x[i];

        g_xxyyy_xxyyyz_0[i] = g_yyy_xxyyyz_0[i] * fbe_0 - g_yyy_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyyy_xyyyz_1[i] * fe_0 + g_xyyy_xxyyyz_1[i] * pa_x[i];

        g_xxyyy_xxyyzz_0[i] = g_yyy_xxyyzz_0[i] * fbe_0 - g_yyy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyyy_xyyzz_1[i] * fe_0 + g_xyyy_xxyyzz_1[i] * pa_x[i];

        g_xxyyy_xxyzzz_0[i] = g_yyy_xxyzzz_0[i] * fbe_0 - g_yyy_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyyy_xyzzz_1[i] * fe_0 + g_xyyy_xxyzzz_1[i] * pa_x[i];

        g_xxyyy_xxzzzz_0[i] = 2.0 * g_xxy_xxzzzz_0[i] * fbe_0 - 2.0 * g_xxy_xxzzzz_1[i] * fz_be_0 + g_xxyy_xxzzzz_1[i] * pa_y[i];

        g_xxyyy_xyyyyy_0[i] = g_yyy_xyyyyy_0[i] * fbe_0 - g_yyy_xyyyyy_1[i] * fz_be_0 + g_xyyy_yyyyy_1[i] * fe_0 + g_xyyy_xyyyyy_1[i] * pa_x[i];

        g_xxyyy_xyyyyz_0[i] = g_yyy_xyyyyz_0[i] * fbe_0 - g_yyy_xyyyyz_1[i] * fz_be_0 + g_xyyy_yyyyz_1[i] * fe_0 + g_xyyy_xyyyyz_1[i] * pa_x[i];

        g_xxyyy_xyyyzz_0[i] = g_yyy_xyyyzz_0[i] * fbe_0 - g_yyy_xyyyzz_1[i] * fz_be_0 + g_xyyy_yyyzz_1[i] * fe_0 + g_xyyy_xyyyzz_1[i] * pa_x[i];

        g_xxyyy_xyyzzz_0[i] = g_yyy_xyyzzz_0[i] * fbe_0 - g_yyy_xyyzzz_1[i] * fz_be_0 + g_xyyy_yyzzz_1[i] * fe_0 + g_xyyy_xyyzzz_1[i] * pa_x[i];

        g_xxyyy_xyzzzz_0[i] = g_yyy_xyzzzz_0[i] * fbe_0 - g_yyy_xyzzzz_1[i] * fz_be_0 + g_xyyy_yzzzz_1[i] * fe_0 + g_xyyy_xyzzzz_1[i] * pa_x[i];

        g_xxyyy_xzzzzz_0[i] = 2.0 * g_xxy_xzzzzz_0[i] * fbe_0 - 2.0 * g_xxy_xzzzzz_1[i] * fz_be_0 + g_xxyy_xzzzzz_1[i] * pa_y[i];

        g_xxyyy_yyyyyy_0[i] = g_yyy_yyyyyy_0[i] * fbe_0 - g_yyy_yyyyyy_1[i] * fz_be_0 + g_xyyy_yyyyyy_1[i] * pa_x[i];

        g_xxyyy_yyyyyz_0[i] = g_yyy_yyyyyz_0[i] * fbe_0 - g_yyy_yyyyyz_1[i] * fz_be_0 + g_xyyy_yyyyyz_1[i] * pa_x[i];

        g_xxyyy_yyyyzz_0[i] = g_yyy_yyyyzz_0[i] * fbe_0 - g_yyy_yyyyzz_1[i] * fz_be_0 + g_xyyy_yyyyzz_1[i] * pa_x[i];

        g_xxyyy_yyyzzz_0[i] = g_yyy_yyyzzz_0[i] * fbe_0 - g_yyy_yyyzzz_1[i] * fz_be_0 + g_xyyy_yyyzzz_1[i] * pa_x[i];

        g_xxyyy_yyzzzz_0[i] = g_yyy_yyzzzz_0[i] * fbe_0 - g_yyy_yyzzzz_1[i] * fz_be_0 + g_xyyy_yyzzzz_1[i] * pa_x[i];

        g_xxyyy_yzzzzz_0[i] = g_yyy_yzzzzz_0[i] * fbe_0 - g_yyy_yzzzzz_1[i] * fz_be_0 + g_xyyy_yzzzzz_1[i] * pa_x[i];

        g_xxyyy_zzzzzz_0[i] = g_yyy_zzzzzz_0[i] * fbe_0 - g_yyy_zzzzzz_1[i] * fz_be_0 + g_xyyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 196-224 components of targeted buffer : HI

    auto g_xxyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 196);

    auto g_xxyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 197);

    auto g_xxyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 198);

    auto g_xxyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 199);

    auto g_xxyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 200);

    auto g_xxyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 201);

    auto g_xxyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 202);

    auto g_xxyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 203);

    auto g_xxyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 204);

    auto g_xxyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 205);

    auto g_xxyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 206);

    auto g_xxyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 207);

    auto g_xxyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 208);

    auto g_xxyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 209);

    auto g_xxyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 210);

    auto g_xxyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 211);

    auto g_xxyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 212);

    auto g_xxyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 213);

    auto g_xxyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 214);

    auto g_xxyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 215);

    auto g_xxyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 216);

    auto g_xxyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 217);

    auto g_xxyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 218);

    auto g_xxyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 219);

    auto g_xxyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 220);

    auto g_xxyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 221);

    auto g_xxyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 222);

    auto g_xxyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 223);

    #pragma omp simd aligned(g_xxyy_xxxxx_1, g_xxyy_xxxxxx_1, g_xxyy_xxxxxy_1, g_xxyy_xxxxxz_1, g_xxyy_xxxxy_1, g_xxyy_xxxxyy_1, g_xxyy_xxxxyz_1, g_xxyy_xxxxz_1, g_xxyy_xxxxzz_1, g_xxyy_xxxyy_1, g_xxyy_xxxyyy_1, g_xxyy_xxxyyz_1, g_xxyy_xxxyz_1, g_xxyy_xxxyzz_1, g_xxyy_xxxzz_1, g_xxyy_xxxzzz_1, g_xxyy_xxyyy_1, g_xxyy_xxyyyy_1, g_xxyy_xxyyyz_1, g_xxyy_xxyyz_1, g_xxyy_xxyyzz_1, g_xxyy_xxyzz_1, g_xxyy_xxyzzz_1, g_xxyy_xxzzz_1, g_xxyy_xxzzzz_1, g_xxyy_xyyyy_1, g_xxyy_xyyyyy_1, g_xxyy_xyyyyz_1, g_xxyy_xyyyz_1, g_xxyy_xyyyzz_1, g_xxyy_xyyzz_1, g_xxyy_xyyzzz_1, g_xxyy_xyzzz_1, g_xxyy_xyzzzz_1, g_xxyy_xzzzz_1, g_xxyy_xzzzzz_1, g_xxyy_yyyyy_1, g_xxyy_yyyyyy_1, g_xxyy_yyyyyz_1, g_xxyy_yyyyz_1, g_xxyy_yyyyzz_1, g_xxyy_yyyzz_1, g_xxyy_yyyzzz_1, g_xxyy_yyzzz_1, g_xxyy_yyzzzz_1, g_xxyy_yzzzz_1, g_xxyy_yzzzzz_1, g_xxyy_zzzzz_1, g_xxyy_zzzzzz_1, g_xxyyz_xxxxxx_0, g_xxyyz_xxxxxy_0, g_xxyyz_xxxxxz_0, g_xxyyz_xxxxyy_0, g_xxyyz_xxxxyz_0, g_xxyyz_xxxxzz_0, g_xxyyz_xxxyyy_0, g_xxyyz_xxxyyz_0, g_xxyyz_xxxyzz_0, g_xxyyz_xxxzzz_0, g_xxyyz_xxyyyy_0, g_xxyyz_xxyyyz_0, g_xxyyz_xxyyzz_0, g_xxyyz_xxyzzz_0, g_xxyyz_xxzzzz_0, g_xxyyz_xyyyyy_0, g_xxyyz_xyyyyz_0, g_xxyyz_xyyyzz_0, g_xxyyz_xyyzzz_0, g_xxyyz_xyzzzz_0, g_xxyyz_xzzzzz_0, g_xxyyz_yyyyyy_0, g_xxyyz_yyyyyz_0, g_xxyyz_yyyyzz_0, g_xxyyz_yyyzzz_0, g_xxyyz_yyzzzz_0, g_xxyyz_yzzzzz_0, g_xxyyz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyz_xxxxxx_0[i] = g_xxyy_xxxxxx_1[i] * pa_z[i];

        g_xxyyz_xxxxxy_0[i] = g_xxyy_xxxxxy_1[i] * pa_z[i];

        g_xxyyz_xxxxxz_0[i] = g_xxyy_xxxxx_1[i] * fe_0 + g_xxyy_xxxxxz_1[i] * pa_z[i];

        g_xxyyz_xxxxyy_0[i] = g_xxyy_xxxxyy_1[i] * pa_z[i];

        g_xxyyz_xxxxyz_0[i] = g_xxyy_xxxxy_1[i] * fe_0 + g_xxyy_xxxxyz_1[i] * pa_z[i];

        g_xxyyz_xxxxzz_0[i] = 2.0 * g_xxyy_xxxxz_1[i] * fe_0 + g_xxyy_xxxxzz_1[i] * pa_z[i];

        g_xxyyz_xxxyyy_0[i] = g_xxyy_xxxyyy_1[i] * pa_z[i];

        g_xxyyz_xxxyyz_0[i] = g_xxyy_xxxyy_1[i] * fe_0 + g_xxyy_xxxyyz_1[i] * pa_z[i];

        g_xxyyz_xxxyzz_0[i] = 2.0 * g_xxyy_xxxyz_1[i] * fe_0 + g_xxyy_xxxyzz_1[i] * pa_z[i];

        g_xxyyz_xxxzzz_0[i] = 3.0 * g_xxyy_xxxzz_1[i] * fe_0 + g_xxyy_xxxzzz_1[i] * pa_z[i];

        g_xxyyz_xxyyyy_0[i] = g_xxyy_xxyyyy_1[i] * pa_z[i];

        g_xxyyz_xxyyyz_0[i] = g_xxyy_xxyyy_1[i] * fe_0 + g_xxyy_xxyyyz_1[i] * pa_z[i];

        g_xxyyz_xxyyzz_0[i] = 2.0 * g_xxyy_xxyyz_1[i] * fe_0 + g_xxyy_xxyyzz_1[i] * pa_z[i];

        g_xxyyz_xxyzzz_0[i] = 3.0 * g_xxyy_xxyzz_1[i] * fe_0 + g_xxyy_xxyzzz_1[i] * pa_z[i];

        g_xxyyz_xxzzzz_0[i] = 4.0 * g_xxyy_xxzzz_1[i] * fe_0 + g_xxyy_xxzzzz_1[i] * pa_z[i];

        g_xxyyz_xyyyyy_0[i] = g_xxyy_xyyyyy_1[i] * pa_z[i];

        g_xxyyz_xyyyyz_0[i] = g_xxyy_xyyyy_1[i] * fe_0 + g_xxyy_xyyyyz_1[i] * pa_z[i];

        g_xxyyz_xyyyzz_0[i] = 2.0 * g_xxyy_xyyyz_1[i] * fe_0 + g_xxyy_xyyyzz_1[i] * pa_z[i];

        g_xxyyz_xyyzzz_0[i] = 3.0 * g_xxyy_xyyzz_1[i] * fe_0 + g_xxyy_xyyzzz_1[i] * pa_z[i];

        g_xxyyz_xyzzzz_0[i] = 4.0 * g_xxyy_xyzzz_1[i] * fe_0 + g_xxyy_xyzzzz_1[i] * pa_z[i];

        g_xxyyz_xzzzzz_0[i] = 5.0 * g_xxyy_xzzzz_1[i] * fe_0 + g_xxyy_xzzzzz_1[i] * pa_z[i];

        g_xxyyz_yyyyyy_0[i] = g_xxyy_yyyyyy_1[i] * pa_z[i];

        g_xxyyz_yyyyyz_0[i] = g_xxyy_yyyyy_1[i] * fe_0 + g_xxyy_yyyyyz_1[i] * pa_z[i];

        g_xxyyz_yyyyzz_0[i] = 2.0 * g_xxyy_yyyyz_1[i] * fe_0 + g_xxyy_yyyyzz_1[i] * pa_z[i];

        g_xxyyz_yyyzzz_0[i] = 3.0 * g_xxyy_yyyzz_1[i] * fe_0 + g_xxyy_yyyzzz_1[i] * pa_z[i];

        g_xxyyz_yyzzzz_0[i] = 4.0 * g_xxyy_yyzzz_1[i] * fe_0 + g_xxyy_yyzzzz_1[i] * pa_z[i];

        g_xxyyz_yzzzzz_0[i] = 5.0 * g_xxyy_yzzzz_1[i] * fe_0 + g_xxyy_yzzzzz_1[i] * pa_z[i];

        g_xxyyz_zzzzzz_0[i] = 6.0 * g_xxyy_zzzzz_1[i] * fe_0 + g_xxyy_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 224-252 components of targeted buffer : HI

    auto g_xxyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 224);

    auto g_xxyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 225);

    auto g_xxyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 226);

    auto g_xxyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 227);

    auto g_xxyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 228);

    auto g_xxyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 229);

    auto g_xxyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 230);

    auto g_xxyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 231);

    auto g_xxyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 232);

    auto g_xxyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 233);

    auto g_xxyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 234);

    auto g_xxyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 235);

    auto g_xxyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 236);

    auto g_xxyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 237);

    auto g_xxyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 238);

    auto g_xxyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 239);

    auto g_xxyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 240);

    auto g_xxyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 241);

    auto g_xxyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 242);

    auto g_xxyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 243);

    auto g_xxyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 244);

    auto g_xxyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 245);

    auto g_xxyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 246);

    auto g_xxyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 247);

    auto g_xxyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 248);

    auto g_xxyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 249);

    auto g_xxyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 250);

    auto g_xxyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 251);

    #pragma omp simd aligned(g_xxyzz_xxxxxx_0, g_xxyzz_xxxxxy_0, g_xxyzz_xxxxxz_0, g_xxyzz_xxxxyy_0, g_xxyzz_xxxxyz_0, g_xxyzz_xxxxzz_0, g_xxyzz_xxxyyy_0, g_xxyzz_xxxyyz_0, g_xxyzz_xxxyzz_0, g_xxyzz_xxxzzz_0, g_xxyzz_xxyyyy_0, g_xxyzz_xxyyyz_0, g_xxyzz_xxyyzz_0, g_xxyzz_xxyzzz_0, g_xxyzz_xxzzzz_0, g_xxyzz_xyyyyy_0, g_xxyzz_xyyyyz_0, g_xxyzz_xyyyzz_0, g_xxyzz_xyyzzz_0, g_xxyzz_xyzzzz_0, g_xxyzz_xzzzzz_0, g_xxyzz_yyyyyy_0, g_xxyzz_yyyyyz_0, g_xxyzz_yyyyzz_0, g_xxyzz_yyyzzz_0, g_xxyzz_yyzzzz_0, g_xxyzz_yzzzzz_0, g_xxyzz_zzzzzz_0, g_xxzz_xxxxx_1, g_xxzz_xxxxxx_1, g_xxzz_xxxxxy_1, g_xxzz_xxxxxz_1, g_xxzz_xxxxy_1, g_xxzz_xxxxyy_1, g_xxzz_xxxxyz_1, g_xxzz_xxxxz_1, g_xxzz_xxxxzz_1, g_xxzz_xxxyy_1, g_xxzz_xxxyyy_1, g_xxzz_xxxyyz_1, g_xxzz_xxxyz_1, g_xxzz_xxxyzz_1, g_xxzz_xxxzz_1, g_xxzz_xxxzzz_1, g_xxzz_xxyyy_1, g_xxzz_xxyyyy_1, g_xxzz_xxyyyz_1, g_xxzz_xxyyz_1, g_xxzz_xxyyzz_1, g_xxzz_xxyzz_1, g_xxzz_xxyzzz_1, g_xxzz_xxzzz_1, g_xxzz_xxzzzz_1, g_xxzz_xyyyy_1, g_xxzz_xyyyyy_1, g_xxzz_xyyyyz_1, g_xxzz_xyyyz_1, g_xxzz_xyyyzz_1, g_xxzz_xyyzz_1, g_xxzz_xyyzzz_1, g_xxzz_xyzzz_1, g_xxzz_xyzzzz_1, g_xxzz_xzzzz_1, g_xxzz_xzzzzz_1, g_xxzz_yyyyy_1, g_xxzz_yyyyyy_1, g_xxzz_yyyyyz_1, g_xxzz_yyyyz_1, g_xxzz_yyyyzz_1, g_xxzz_yyyzz_1, g_xxzz_yyyzzz_1, g_xxzz_yyzzz_1, g_xxzz_yyzzzz_1, g_xxzz_yzzzz_1, g_xxzz_yzzzzz_1, g_xxzz_zzzzz_1, g_xxzz_zzzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzz_xxxxxx_0[i] = g_xxzz_xxxxxx_1[i] * pa_y[i];

        g_xxyzz_xxxxxy_0[i] = g_xxzz_xxxxx_1[i] * fe_0 + g_xxzz_xxxxxy_1[i] * pa_y[i];

        g_xxyzz_xxxxxz_0[i] = g_xxzz_xxxxxz_1[i] * pa_y[i];

        g_xxyzz_xxxxyy_0[i] = 2.0 * g_xxzz_xxxxy_1[i] * fe_0 + g_xxzz_xxxxyy_1[i] * pa_y[i];

        g_xxyzz_xxxxyz_0[i] = g_xxzz_xxxxz_1[i] * fe_0 + g_xxzz_xxxxyz_1[i] * pa_y[i];

        g_xxyzz_xxxxzz_0[i] = g_xxzz_xxxxzz_1[i] * pa_y[i];

        g_xxyzz_xxxyyy_0[i] = 3.0 * g_xxzz_xxxyy_1[i] * fe_0 + g_xxzz_xxxyyy_1[i] * pa_y[i];

        g_xxyzz_xxxyyz_0[i] = 2.0 * g_xxzz_xxxyz_1[i] * fe_0 + g_xxzz_xxxyyz_1[i] * pa_y[i];

        g_xxyzz_xxxyzz_0[i] = g_xxzz_xxxzz_1[i] * fe_0 + g_xxzz_xxxyzz_1[i] * pa_y[i];

        g_xxyzz_xxxzzz_0[i] = g_xxzz_xxxzzz_1[i] * pa_y[i];

        g_xxyzz_xxyyyy_0[i] = 4.0 * g_xxzz_xxyyy_1[i] * fe_0 + g_xxzz_xxyyyy_1[i] * pa_y[i];

        g_xxyzz_xxyyyz_0[i] = 3.0 * g_xxzz_xxyyz_1[i] * fe_0 + g_xxzz_xxyyyz_1[i] * pa_y[i];

        g_xxyzz_xxyyzz_0[i] = 2.0 * g_xxzz_xxyzz_1[i] * fe_0 + g_xxzz_xxyyzz_1[i] * pa_y[i];

        g_xxyzz_xxyzzz_0[i] = g_xxzz_xxzzz_1[i] * fe_0 + g_xxzz_xxyzzz_1[i] * pa_y[i];

        g_xxyzz_xxzzzz_0[i] = g_xxzz_xxzzzz_1[i] * pa_y[i];

        g_xxyzz_xyyyyy_0[i] = 5.0 * g_xxzz_xyyyy_1[i] * fe_0 + g_xxzz_xyyyyy_1[i] * pa_y[i];

        g_xxyzz_xyyyyz_0[i] = 4.0 * g_xxzz_xyyyz_1[i] * fe_0 + g_xxzz_xyyyyz_1[i] * pa_y[i];

        g_xxyzz_xyyyzz_0[i] = 3.0 * g_xxzz_xyyzz_1[i] * fe_0 + g_xxzz_xyyyzz_1[i] * pa_y[i];

        g_xxyzz_xyyzzz_0[i] = 2.0 * g_xxzz_xyzzz_1[i] * fe_0 + g_xxzz_xyyzzz_1[i] * pa_y[i];

        g_xxyzz_xyzzzz_0[i] = g_xxzz_xzzzz_1[i] * fe_0 + g_xxzz_xyzzzz_1[i] * pa_y[i];

        g_xxyzz_xzzzzz_0[i] = g_xxzz_xzzzzz_1[i] * pa_y[i];

        g_xxyzz_yyyyyy_0[i] = 6.0 * g_xxzz_yyyyy_1[i] * fe_0 + g_xxzz_yyyyyy_1[i] * pa_y[i];

        g_xxyzz_yyyyyz_0[i] = 5.0 * g_xxzz_yyyyz_1[i] * fe_0 + g_xxzz_yyyyyz_1[i] * pa_y[i];

        g_xxyzz_yyyyzz_0[i] = 4.0 * g_xxzz_yyyzz_1[i] * fe_0 + g_xxzz_yyyyzz_1[i] * pa_y[i];

        g_xxyzz_yyyzzz_0[i] = 3.0 * g_xxzz_yyzzz_1[i] * fe_0 + g_xxzz_yyyzzz_1[i] * pa_y[i];

        g_xxyzz_yyzzzz_0[i] = 2.0 * g_xxzz_yzzzz_1[i] * fe_0 + g_xxzz_yyzzzz_1[i] * pa_y[i];

        g_xxyzz_yzzzzz_0[i] = g_xxzz_zzzzz_1[i] * fe_0 + g_xxzz_yzzzzz_1[i] * pa_y[i];

        g_xxyzz_zzzzzz_0[i] = g_xxzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 252-280 components of targeted buffer : HI

    auto g_xxzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 252);

    auto g_xxzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 253);

    auto g_xxzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 254);

    auto g_xxzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 255);

    auto g_xxzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 256);

    auto g_xxzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 257);

    auto g_xxzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 258);

    auto g_xxzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 259);

    auto g_xxzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 260);

    auto g_xxzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 261);

    auto g_xxzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 262);

    auto g_xxzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 263);

    auto g_xxzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 264);

    auto g_xxzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 265);

    auto g_xxzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 266);

    auto g_xxzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 267);

    auto g_xxzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 268);

    auto g_xxzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 269);

    auto g_xxzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 270);

    auto g_xxzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 271);

    auto g_xxzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 272);

    auto g_xxzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 273);

    auto g_xxzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 274);

    auto g_xxzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 275);

    auto g_xxzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 276);

    auto g_xxzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 277);

    auto g_xxzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 278);

    auto g_xxzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 279);

    #pragma omp simd aligned(g_xxz_xxxxxx_0, g_xxz_xxxxxx_1, g_xxz_xxxxxy_0, g_xxz_xxxxxy_1, g_xxz_xxxxyy_0, g_xxz_xxxxyy_1, g_xxz_xxxyyy_0, g_xxz_xxxyyy_1, g_xxz_xxyyyy_0, g_xxz_xxyyyy_1, g_xxz_xyyyyy_0, g_xxz_xyyyyy_1, g_xxzz_xxxxxx_1, g_xxzz_xxxxxy_1, g_xxzz_xxxxyy_1, g_xxzz_xxxyyy_1, g_xxzz_xxyyyy_1, g_xxzz_xyyyyy_1, g_xxzzz_xxxxxx_0, g_xxzzz_xxxxxy_0, g_xxzzz_xxxxxz_0, g_xxzzz_xxxxyy_0, g_xxzzz_xxxxyz_0, g_xxzzz_xxxxzz_0, g_xxzzz_xxxyyy_0, g_xxzzz_xxxyyz_0, g_xxzzz_xxxyzz_0, g_xxzzz_xxxzzz_0, g_xxzzz_xxyyyy_0, g_xxzzz_xxyyyz_0, g_xxzzz_xxyyzz_0, g_xxzzz_xxyzzz_0, g_xxzzz_xxzzzz_0, g_xxzzz_xyyyyy_0, g_xxzzz_xyyyyz_0, g_xxzzz_xyyyzz_0, g_xxzzz_xyyzzz_0, g_xxzzz_xyzzzz_0, g_xxzzz_xzzzzz_0, g_xxzzz_yyyyyy_0, g_xxzzz_yyyyyz_0, g_xxzzz_yyyyzz_0, g_xxzzz_yyyzzz_0, g_xxzzz_yyzzzz_0, g_xxzzz_yzzzzz_0, g_xxzzz_zzzzzz_0, g_xzzz_xxxxxz_1, g_xzzz_xxxxyz_1, g_xzzz_xxxxz_1, g_xzzz_xxxxzz_1, g_xzzz_xxxyyz_1, g_xzzz_xxxyz_1, g_xzzz_xxxyzz_1, g_xzzz_xxxzz_1, g_xzzz_xxxzzz_1, g_xzzz_xxyyyz_1, g_xzzz_xxyyz_1, g_xzzz_xxyyzz_1, g_xzzz_xxyzz_1, g_xzzz_xxyzzz_1, g_xzzz_xxzzz_1, g_xzzz_xxzzzz_1, g_xzzz_xyyyyz_1, g_xzzz_xyyyz_1, g_xzzz_xyyyzz_1, g_xzzz_xyyzz_1, g_xzzz_xyyzzz_1, g_xzzz_xyzzz_1, g_xzzz_xyzzzz_1, g_xzzz_xzzzz_1, g_xzzz_xzzzzz_1, g_xzzz_yyyyyy_1, g_xzzz_yyyyyz_1, g_xzzz_yyyyz_1, g_xzzz_yyyyzz_1, g_xzzz_yyyzz_1, g_xzzz_yyyzzz_1, g_xzzz_yyzzz_1, g_xzzz_yyzzzz_1, g_xzzz_yzzzz_1, g_xzzz_yzzzzz_1, g_xzzz_zzzzz_1, g_xzzz_zzzzzz_1, g_zzz_xxxxxz_0, g_zzz_xxxxxz_1, g_zzz_xxxxyz_0, g_zzz_xxxxyz_1, g_zzz_xxxxzz_0, g_zzz_xxxxzz_1, g_zzz_xxxyyz_0, g_zzz_xxxyyz_1, g_zzz_xxxyzz_0, g_zzz_xxxyzz_1, g_zzz_xxxzzz_0, g_zzz_xxxzzz_1, g_zzz_xxyyyz_0, g_zzz_xxyyyz_1, g_zzz_xxyyzz_0, g_zzz_xxyyzz_1, g_zzz_xxyzzz_0, g_zzz_xxyzzz_1, g_zzz_xxzzzz_0, g_zzz_xxzzzz_1, g_zzz_xyyyyz_0, g_zzz_xyyyyz_1, g_zzz_xyyyzz_0, g_zzz_xyyyzz_1, g_zzz_xyyzzz_0, g_zzz_xyyzzz_1, g_zzz_xyzzzz_0, g_zzz_xyzzzz_1, g_zzz_xzzzzz_0, g_zzz_xzzzzz_1, g_zzz_yyyyyy_0, g_zzz_yyyyyy_1, g_zzz_yyyyyz_0, g_zzz_yyyyyz_1, g_zzz_yyyyzz_0, g_zzz_yyyyzz_1, g_zzz_yyyzzz_0, g_zzz_yyyzzz_1, g_zzz_yyzzzz_0, g_zzz_yyzzzz_1, g_zzz_yzzzzz_0, g_zzz_yzzzzz_1, g_zzz_zzzzzz_0, g_zzz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzz_xxxxxx_0[i] = 2.0 * g_xxz_xxxxxx_0[i] * fbe_0 - 2.0 * g_xxz_xxxxxx_1[i] * fz_be_0 + g_xxzz_xxxxxx_1[i] * pa_z[i];

        g_xxzzz_xxxxxy_0[i] = 2.0 * g_xxz_xxxxxy_0[i] * fbe_0 - 2.0 * g_xxz_xxxxxy_1[i] * fz_be_0 + g_xxzz_xxxxxy_1[i] * pa_z[i];

        g_xxzzz_xxxxxz_0[i] = g_zzz_xxxxxz_0[i] * fbe_0 - g_zzz_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xzzz_xxxxz_1[i] * fe_0 + g_xzzz_xxxxxz_1[i] * pa_x[i];

        g_xxzzz_xxxxyy_0[i] = 2.0 * g_xxz_xxxxyy_0[i] * fbe_0 - 2.0 * g_xxz_xxxxyy_1[i] * fz_be_0 + g_xxzz_xxxxyy_1[i] * pa_z[i];

        g_xxzzz_xxxxyz_0[i] = g_zzz_xxxxyz_0[i] * fbe_0 - g_zzz_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xzzz_xxxyz_1[i] * fe_0 + g_xzzz_xxxxyz_1[i] * pa_x[i];

        g_xxzzz_xxxxzz_0[i] = g_zzz_xxxxzz_0[i] * fbe_0 - g_zzz_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xzzz_xxxzz_1[i] * fe_0 + g_xzzz_xxxxzz_1[i] * pa_x[i];

        g_xxzzz_xxxyyy_0[i] = 2.0 * g_xxz_xxxyyy_0[i] * fbe_0 - 2.0 * g_xxz_xxxyyy_1[i] * fz_be_0 + g_xxzz_xxxyyy_1[i] * pa_z[i];

        g_xxzzz_xxxyyz_0[i] = g_zzz_xxxyyz_0[i] * fbe_0 - g_zzz_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xzzz_xxyyz_1[i] * fe_0 + g_xzzz_xxxyyz_1[i] * pa_x[i];

        g_xxzzz_xxxyzz_0[i] = g_zzz_xxxyzz_0[i] * fbe_0 - g_zzz_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xzzz_xxyzz_1[i] * fe_0 + g_xzzz_xxxyzz_1[i] * pa_x[i];

        g_xxzzz_xxxzzz_0[i] = g_zzz_xxxzzz_0[i] * fbe_0 - g_zzz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xzzz_xxzzz_1[i] * fe_0 + g_xzzz_xxxzzz_1[i] * pa_x[i];

        g_xxzzz_xxyyyy_0[i] = 2.0 * g_xxz_xxyyyy_0[i] * fbe_0 - 2.0 * g_xxz_xxyyyy_1[i] * fz_be_0 + g_xxzz_xxyyyy_1[i] * pa_z[i];

        g_xxzzz_xxyyyz_0[i] = g_zzz_xxyyyz_0[i] * fbe_0 - g_zzz_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xzzz_xyyyz_1[i] * fe_0 + g_xzzz_xxyyyz_1[i] * pa_x[i];

        g_xxzzz_xxyyzz_0[i] = g_zzz_xxyyzz_0[i] * fbe_0 - g_zzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xzzz_xyyzz_1[i] * fe_0 + g_xzzz_xxyyzz_1[i] * pa_x[i];

        g_xxzzz_xxyzzz_0[i] = g_zzz_xxyzzz_0[i] * fbe_0 - g_zzz_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_xyzzz_1[i] * fe_0 + g_xzzz_xxyzzz_1[i] * pa_x[i];

        g_xxzzz_xxzzzz_0[i] = g_zzz_xxzzzz_0[i] * fbe_0 - g_zzz_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xzzz_xzzzz_1[i] * fe_0 + g_xzzz_xxzzzz_1[i] * pa_x[i];

        g_xxzzz_xyyyyy_0[i] = 2.0 * g_xxz_xyyyyy_0[i] * fbe_0 - 2.0 * g_xxz_xyyyyy_1[i] * fz_be_0 + g_xxzz_xyyyyy_1[i] * pa_z[i];

        g_xxzzz_xyyyyz_0[i] = g_zzz_xyyyyz_0[i] * fbe_0 - g_zzz_xyyyyz_1[i] * fz_be_0 + g_xzzz_yyyyz_1[i] * fe_0 + g_xzzz_xyyyyz_1[i] * pa_x[i];

        g_xxzzz_xyyyzz_0[i] = g_zzz_xyyyzz_0[i] * fbe_0 - g_zzz_xyyyzz_1[i] * fz_be_0 + g_xzzz_yyyzz_1[i] * fe_0 + g_xzzz_xyyyzz_1[i] * pa_x[i];

        g_xxzzz_xyyzzz_0[i] = g_zzz_xyyzzz_0[i] * fbe_0 - g_zzz_xyyzzz_1[i] * fz_be_0 + g_xzzz_yyzzz_1[i] * fe_0 + g_xzzz_xyyzzz_1[i] * pa_x[i];

        g_xxzzz_xyzzzz_0[i] = g_zzz_xyzzzz_0[i] * fbe_0 - g_zzz_xyzzzz_1[i] * fz_be_0 + g_xzzz_yzzzz_1[i] * fe_0 + g_xzzz_xyzzzz_1[i] * pa_x[i];

        g_xxzzz_xzzzzz_0[i] = g_zzz_xzzzzz_0[i] * fbe_0 - g_zzz_xzzzzz_1[i] * fz_be_0 + g_xzzz_zzzzz_1[i] * fe_0 + g_xzzz_xzzzzz_1[i] * pa_x[i];

        g_xxzzz_yyyyyy_0[i] = g_zzz_yyyyyy_0[i] * fbe_0 - g_zzz_yyyyyy_1[i] * fz_be_0 + g_xzzz_yyyyyy_1[i] * pa_x[i];

        g_xxzzz_yyyyyz_0[i] = g_zzz_yyyyyz_0[i] * fbe_0 - g_zzz_yyyyyz_1[i] * fz_be_0 + g_xzzz_yyyyyz_1[i] * pa_x[i];

        g_xxzzz_yyyyzz_0[i] = g_zzz_yyyyzz_0[i] * fbe_0 - g_zzz_yyyyzz_1[i] * fz_be_0 + g_xzzz_yyyyzz_1[i] * pa_x[i];

        g_xxzzz_yyyzzz_0[i] = g_zzz_yyyzzz_0[i] * fbe_0 - g_zzz_yyyzzz_1[i] * fz_be_0 + g_xzzz_yyyzzz_1[i] * pa_x[i];

        g_xxzzz_yyzzzz_0[i] = g_zzz_yyzzzz_0[i] * fbe_0 - g_zzz_yyzzzz_1[i] * fz_be_0 + g_xzzz_yyzzzz_1[i] * pa_x[i];

        g_xxzzz_yzzzzz_0[i] = g_zzz_yzzzzz_0[i] * fbe_0 - g_zzz_yzzzzz_1[i] * fz_be_0 + g_xzzz_yzzzzz_1[i] * pa_x[i];

        g_xxzzz_zzzzzz_0[i] = g_zzz_zzzzzz_0[i] * fbe_0 - g_zzz_zzzzzz_1[i] * fz_be_0 + g_xzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 280-308 components of targeted buffer : HI

    auto g_xyyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 280);

    auto g_xyyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 281);

    auto g_xyyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 282);

    auto g_xyyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 283);

    auto g_xyyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 284);

    auto g_xyyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 285);

    auto g_xyyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 286);

    auto g_xyyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 287);

    auto g_xyyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 288);

    auto g_xyyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 289);

    auto g_xyyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 290);

    auto g_xyyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 291);

    auto g_xyyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 292);

    auto g_xyyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 293);

    auto g_xyyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 294);

    auto g_xyyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 295);

    auto g_xyyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 296);

    auto g_xyyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 297);

    auto g_xyyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 298);

    auto g_xyyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 299);

    auto g_xyyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 300);

    auto g_xyyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 301);

    auto g_xyyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 302);

    auto g_xyyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 303);

    auto g_xyyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 304);

    auto g_xyyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 305);

    auto g_xyyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 306);

    auto g_xyyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 307);

    #pragma omp simd aligned(g_xyyyy_xxxxxx_0, g_xyyyy_xxxxxy_0, g_xyyyy_xxxxxz_0, g_xyyyy_xxxxyy_0, g_xyyyy_xxxxyz_0, g_xyyyy_xxxxzz_0, g_xyyyy_xxxyyy_0, g_xyyyy_xxxyyz_0, g_xyyyy_xxxyzz_0, g_xyyyy_xxxzzz_0, g_xyyyy_xxyyyy_0, g_xyyyy_xxyyyz_0, g_xyyyy_xxyyzz_0, g_xyyyy_xxyzzz_0, g_xyyyy_xxzzzz_0, g_xyyyy_xyyyyy_0, g_xyyyy_xyyyyz_0, g_xyyyy_xyyyzz_0, g_xyyyy_xyyzzz_0, g_xyyyy_xyzzzz_0, g_xyyyy_xzzzzz_0, g_xyyyy_yyyyyy_0, g_xyyyy_yyyyyz_0, g_xyyyy_yyyyzz_0, g_xyyyy_yyyzzz_0, g_xyyyy_yyzzzz_0, g_xyyyy_yzzzzz_0, g_xyyyy_zzzzzz_0, g_yyyy_xxxxx_1, g_yyyy_xxxxxx_1, g_yyyy_xxxxxy_1, g_yyyy_xxxxxz_1, g_yyyy_xxxxy_1, g_yyyy_xxxxyy_1, g_yyyy_xxxxyz_1, g_yyyy_xxxxz_1, g_yyyy_xxxxzz_1, g_yyyy_xxxyy_1, g_yyyy_xxxyyy_1, g_yyyy_xxxyyz_1, g_yyyy_xxxyz_1, g_yyyy_xxxyzz_1, g_yyyy_xxxzz_1, g_yyyy_xxxzzz_1, g_yyyy_xxyyy_1, g_yyyy_xxyyyy_1, g_yyyy_xxyyyz_1, g_yyyy_xxyyz_1, g_yyyy_xxyyzz_1, g_yyyy_xxyzz_1, g_yyyy_xxyzzz_1, g_yyyy_xxzzz_1, g_yyyy_xxzzzz_1, g_yyyy_xyyyy_1, g_yyyy_xyyyyy_1, g_yyyy_xyyyyz_1, g_yyyy_xyyyz_1, g_yyyy_xyyyzz_1, g_yyyy_xyyzz_1, g_yyyy_xyyzzz_1, g_yyyy_xyzzz_1, g_yyyy_xyzzzz_1, g_yyyy_xzzzz_1, g_yyyy_xzzzzz_1, g_yyyy_yyyyy_1, g_yyyy_yyyyyy_1, g_yyyy_yyyyyz_1, g_yyyy_yyyyz_1, g_yyyy_yyyyzz_1, g_yyyy_yyyzz_1, g_yyyy_yyyzzz_1, g_yyyy_yyzzz_1, g_yyyy_yyzzzz_1, g_yyyy_yzzzz_1, g_yyyy_yzzzzz_1, g_yyyy_zzzzz_1, g_yyyy_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyy_xxxxxx_0[i] = 6.0 * g_yyyy_xxxxx_1[i] * fe_0 + g_yyyy_xxxxxx_1[i] * pa_x[i];

        g_xyyyy_xxxxxy_0[i] = 5.0 * g_yyyy_xxxxy_1[i] * fe_0 + g_yyyy_xxxxxy_1[i] * pa_x[i];

        g_xyyyy_xxxxxz_0[i] = 5.0 * g_yyyy_xxxxz_1[i] * fe_0 + g_yyyy_xxxxxz_1[i] * pa_x[i];

        g_xyyyy_xxxxyy_0[i] = 4.0 * g_yyyy_xxxyy_1[i] * fe_0 + g_yyyy_xxxxyy_1[i] * pa_x[i];

        g_xyyyy_xxxxyz_0[i] = 4.0 * g_yyyy_xxxyz_1[i] * fe_0 + g_yyyy_xxxxyz_1[i] * pa_x[i];

        g_xyyyy_xxxxzz_0[i] = 4.0 * g_yyyy_xxxzz_1[i] * fe_0 + g_yyyy_xxxxzz_1[i] * pa_x[i];

        g_xyyyy_xxxyyy_0[i] = 3.0 * g_yyyy_xxyyy_1[i] * fe_0 + g_yyyy_xxxyyy_1[i] * pa_x[i];

        g_xyyyy_xxxyyz_0[i] = 3.0 * g_yyyy_xxyyz_1[i] * fe_0 + g_yyyy_xxxyyz_1[i] * pa_x[i];

        g_xyyyy_xxxyzz_0[i] = 3.0 * g_yyyy_xxyzz_1[i] * fe_0 + g_yyyy_xxxyzz_1[i] * pa_x[i];

        g_xyyyy_xxxzzz_0[i] = 3.0 * g_yyyy_xxzzz_1[i] * fe_0 + g_yyyy_xxxzzz_1[i] * pa_x[i];

        g_xyyyy_xxyyyy_0[i] = 2.0 * g_yyyy_xyyyy_1[i] * fe_0 + g_yyyy_xxyyyy_1[i] * pa_x[i];

        g_xyyyy_xxyyyz_0[i] = 2.0 * g_yyyy_xyyyz_1[i] * fe_0 + g_yyyy_xxyyyz_1[i] * pa_x[i];

        g_xyyyy_xxyyzz_0[i] = 2.0 * g_yyyy_xyyzz_1[i] * fe_0 + g_yyyy_xxyyzz_1[i] * pa_x[i];

        g_xyyyy_xxyzzz_0[i] = 2.0 * g_yyyy_xyzzz_1[i] * fe_0 + g_yyyy_xxyzzz_1[i] * pa_x[i];

        g_xyyyy_xxzzzz_0[i] = 2.0 * g_yyyy_xzzzz_1[i] * fe_0 + g_yyyy_xxzzzz_1[i] * pa_x[i];

        g_xyyyy_xyyyyy_0[i] = g_yyyy_yyyyy_1[i] * fe_0 + g_yyyy_xyyyyy_1[i] * pa_x[i];

        g_xyyyy_xyyyyz_0[i] = g_yyyy_yyyyz_1[i] * fe_0 + g_yyyy_xyyyyz_1[i] * pa_x[i];

        g_xyyyy_xyyyzz_0[i] = g_yyyy_yyyzz_1[i] * fe_0 + g_yyyy_xyyyzz_1[i] * pa_x[i];

        g_xyyyy_xyyzzz_0[i] = g_yyyy_yyzzz_1[i] * fe_0 + g_yyyy_xyyzzz_1[i] * pa_x[i];

        g_xyyyy_xyzzzz_0[i] = g_yyyy_yzzzz_1[i] * fe_0 + g_yyyy_xyzzzz_1[i] * pa_x[i];

        g_xyyyy_xzzzzz_0[i] = g_yyyy_zzzzz_1[i] * fe_0 + g_yyyy_xzzzzz_1[i] * pa_x[i];

        g_xyyyy_yyyyyy_0[i] = g_yyyy_yyyyyy_1[i] * pa_x[i];

        g_xyyyy_yyyyyz_0[i] = g_yyyy_yyyyyz_1[i] * pa_x[i];

        g_xyyyy_yyyyzz_0[i] = g_yyyy_yyyyzz_1[i] * pa_x[i];

        g_xyyyy_yyyzzz_0[i] = g_yyyy_yyyzzz_1[i] * pa_x[i];

        g_xyyyy_yyzzzz_0[i] = g_yyyy_yyzzzz_1[i] * pa_x[i];

        g_xyyyy_yzzzzz_0[i] = g_yyyy_yzzzzz_1[i] * pa_x[i];

        g_xyyyy_zzzzzz_0[i] = g_yyyy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 308-336 components of targeted buffer : HI

    auto g_xyyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 308);

    auto g_xyyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 309);

    auto g_xyyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 310);

    auto g_xyyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 311);

    auto g_xyyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 312);

    auto g_xyyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 313);

    auto g_xyyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 314);

    auto g_xyyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 315);

    auto g_xyyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 316);

    auto g_xyyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 317);

    auto g_xyyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 318);

    auto g_xyyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 319);

    auto g_xyyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 320);

    auto g_xyyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 321);

    auto g_xyyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 322);

    auto g_xyyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 323);

    auto g_xyyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 324);

    auto g_xyyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 325);

    auto g_xyyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 326);

    auto g_xyyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 327);

    auto g_xyyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 328);

    auto g_xyyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 329);

    auto g_xyyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 330);

    auto g_xyyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 331);

    auto g_xyyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 332);

    auto g_xyyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 333);

    auto g_xyyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 334);

    auto g_xyyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 335);

    #pragma omp simd aligned(g_xyyy_xxxxxx_1, g_xyyy_xxxxxy_1, g_xyyy_xxxxyy_1, g_xyyy_xxxyyy_1, g_xyyy_xxyyyy_1, g_xyyy_xyyyyy_1, g_xyyyz_xxxxxx_0, g_xyyyz_xxxxxy_0, g_xyyyz_xxxxxz_0, g_xyyyz_xxxxyy_0, g_xyyyz_xxxxyz_0, g_xyyyz_xxxxzz_0, g_xyyyz_xxxyyy_0, g_xyyyz_xxxyyz_0, g_xyyyz_xxxyzz_0, g_xyyyz_xxxzzz_0, g_xyyyz_xxyyyy_0, g_xyyyz_xxyyyz_0, g_xyyyz_xxyyzz_0, g_xyyyz_xxyzzz_0, g_xyyyz_xxzzzz_0, g_xyyyz_xyyyyy_0, g_xyyyz_xyyyyz_0, g_xyyyz_xyyyzz_0, g_xyyyz_xyyzzz_0, g_xyyyz_xyzzzz_0, g_xyyyz_xzzzzz_0, g_xyyyz_yyyyyy_0, g_xyyyz_yyyyyz_0, g_xyyyz_yyyyzz_0, g_xyyyz_yyyzzz_0, g_xyyyz_yyzzzz_0, g_xyyyz_yzzzzz_0, g_xyyyz_zzzzzz_0, g_yyyz_xxxxxz_1, g_yyyz_xxxxyz_1, g_yyyz_xxxxz_1, g_yyyz_xxxxzz_1, g_yyyz_xxxyyz_1, g_yyyz_xxxyz_1, g_yyyz_xxxyzz_1, g_yyyz_xxxzz_1, g_yyyz_xxxzzz_1, g_yyyz_xxyyyz_1, g_yyyz_xxyyz_1, g_yyyz_xxyyzz_1, g_yyyz_xxyzz_1, g_yyyz_xxyzzz_1, g_yyyz_xxzzz_1, g_yyyz_xxzzzz_1, g_yyyz_xyyyyz_1, g_yyyz_xyyyz_1, g_yyyz_xyyyzz_1, g_yyyz_xyyzz_1, g_yyyz_xyyzzz_1, g_yyyz_xyzzz_1, g_yyyz_xyzzzz_1, g_yyyz_xzzzz_1, g_yyyz_xzzzzz_1, g_yyyz_yyyyyy_1, g_yyyz_yyyyyz_1, g_yyyz_yyyyz_1, g_yyyz_yyyyzz_1, g_yyyz_yyyzz_1, g_yyyz_yyyzzz_1, g_yyyz_yyzzz_1, g_yyyz_yyzzzz_1, g_yyyz_yzzzz_1, g_yyyz_yzzzzz_1, g_yyyz_zzzzz_1, g_yyyz_zzzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyz_xxxxxx_0[i] = g_xyyy_xxxxxx_1[i] * pa_z[i];

        g_xyyyz_xxxxxy_0[i] = g_xyyy_xxxxxy_1[i] * pa_z[i];

        g_xyyyz_xxxxxz_0[i] = 5.0 * g_yyyz_xxxxz_1[i] * fe_0 + g_yyyz_xxxxxz_1[i] * pa_x[i];

        g_xyyyz_xxxxyy_0[i] = g_xyyy_xxxxyy_1[i] * pa_z[i];

        g_xyyyz_xxxxyz_0[i] = 4.0 * g_yyyz_xxxyz_1[i] * fe_0 + g_yyyz_xxxxyz_1[i] * pa_x[i];

        g_xyyyz_xxxxzz_0[i] = 4.0 * g_yyyz_xxxzz_1[i] * fe_0 + g_yyyz_xxxxzz_1[i] * pa_x[i];

        g_xyyyz_xxxyyy_0[i] = g_xyyy_xxxyyy_1[i] * pa_z[i];

        g_xyyyz_xxxyyz_0[i] = 3.0 * g_yyyz_xxyyz_1[i] * fe_0 + g_yyyz_xxxyyz_1[i] * pa_x[i];

        g_xyyyz_xxxyzz_0[i] = 3.0 * g_yyyz_xxyzz_1[i] * fe_0 + g_yyyz_xxxyzz_1[i] * pa_x[i];

        g_xyyyz_xxxzzz_0[i] = 3.0 * g_yyyz_xxzzz_1[i] * fe_0 + g_yyyz_xxxzzz_1[i] * pa_x[i];

        g_xyyyz_xxyyyy_0[i] = g_xyyy_xxyyyy_1[i] * pa_z[i];

        g_xyyyz_xxyyyz_0[i] = 2.0 * g_yyyz_xyyyz_1[i] * fe_0 + g_yyyz_xxyyyz_1[i] * pa_x[i];

        g_xyyyz_xxyyzz_0[i] = 2.0 * g_yyyz_xyyzz_1[i] * fe_0 + g_yyyz_xxyyzz_1[i] * pa_x[i];

        g_xyyyz_xxyzzz_0[i] = 2.0 * g_yyyz_xyzzz_1[i] * fe_0 + g_yyyz_xxyzzz_1[i] * pa_x[i];

        g_xyyyz_xxzzzz_0[i] = 2.0 * g_yyyz_xzzzz_1[i] * fe_0 + g_yyyz_xxzzzz_1[i] * pa_x[i];

        g_xyyyz_xyyyyy_0[i] = g_xyyy_xyyyyy_1[i] * pa_z[i];

        g_xyyyz_xyyyyz_0[i] = g_yyyz_yyyyz_1[i] * fe_0 + g_yyyz_xyyyyz_1[i] * pa_x[i];

        g_xyyyz_xyyyzz_0[i] = g_yyyz_yyyzz_1[i] * fe_0 + g_yyyz_xyyyzz_1[i] * pa_x[i];

        g_xyyyz_xyyzzz_0[i] = g_yyyz_yyzzz_1[i] * fe_0 + g_yyyz_xyyzzz_1[i] * pa_x[i];

        g_xyyyz_xyzzzz_0[i] = g_yyyz_yzzzz_1[i] * fe_0 + g_yyyz_xyzzzz_1[i] * pa_x[i];

        g_xyyyz_xzzzzz_0[i] = g_yyyz_zzzzz_1[i] * fe_0 + g_yyyz_xzzzzz_1[i] * pa_x[i];

        g_xyyyz_yyyyyy_0[i] = g_yyyz_yyyyyy_1[i] * pa_x[i];

        g_xyyyz_yyyyyz_0[i] = g_yyyz_yyyyyz_1[i] * pa_x[i];

        g_xyyyz_yyyyzz_0[i] = g_yyyz_yyyyzz_1[i] * pa_x[i];

        g_xyyyz_yyyzzz_0[i] = g_yyyz_yyyzzz_1[i] * pa_x[i];

        g_xyyyz_yyzzzz_0[i] = g_yyyz_yyzzzz_1[i] * pa_x[i];

        g_xyyyz_yzzzzz_0[i] = g_yyyz_yzzzzz_1[i] * pa_x[i];

        g_xyyyz_zzzzzz_0[i] = g_yyyz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 336-364 components of targeted buffer : HI

    auto g_xyyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 336);

    auto g_xyyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 337);

    auto g_xyyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 338);

    auto g_xyyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 339);

    auto g_xyyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 340);

    auto g_xyyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 341);

    auto g_xyyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 342);

    auto g_xyyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 343);

    auto g_xyyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 344);

    auto g_xyyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 345);

    auto g_xyyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 346);

    auto g_xyyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 347);

    auto g_xyyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 348);

    auto g_xyyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 349);

    auto g_xyyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 350);

    auto g_xyyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 351);

    auto g_xyyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 352);

    auto g_xyyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 353);

    auto g_xyyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 354);

    auto g_xyyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 355);

    auto g_xyyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 356);

    auto g_xyyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 357);

    auto g_xyyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 358);

    auto g_xyyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 359);

    auto g_xyyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 360);

    auto g_xyyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 361);

    auto g_xyyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 362);

    auto g_xyyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 363);

    #pragma omp simd aligned(g_xyyzz_xxxxxx_0, g_xyyzz_xxxxxy_0, g_xyyzz_xxxxxz_0, g_xyyzz_xxxxyy_0, g_xyyzz_xxxxyz_0, g_xyyzz_xxxxzz_0, g_xyyzz_xxxyyy_0, g_xyyzz_xxxyyz_0, g_xyyzz_xxxyzz_0, g_xyyzz_xxxzzz_0, g_xyyzz_xxyyyy_0, g_xyyzz_xxyyyz_0, g_xyyzz_xxyyzz_0, g_xyyzz_xxyzzz_0, g_xyyzz_xxzzzz_0, g_xyyzz_xyyyyy_0, g_xyyzz_xyyyyz_0, g_xyyzz_xyyyzz_0, g_xyyzz_xyyzzz_0, g_xyyzz_xyzzzz_0, g_xyyzz_xzzzzz_0, g_xyyzz_yyyyyy_0, g_xyyzz_yyyyyz_0, g_xyyzz_yyyyzz_0, g_xyyzz_yyyzzz_0, g_xyyzz_yyzzzz_0, g_xyyzz_yzzzzz_0, g_xyyzz_zzzzzz_0, g_yyzz_xxxxx_1, g_yyzz_xxxxxx_1, g_yyzz_xxxxxy_1, g_yyzz_xxxxxz_1, g_yyzz_xxxxy_1, g_yyzz_xxxxyy_1, g_yyzz_xxxxyz_1, g_yyzz_xxxxz_1, g_yyzz_xxxxzz_1, g_yyzz_xxxyy_1, g_yyzz_xxxyyy_1, g_yyzz_xxxyyz_1, g_yyzz_xxxyz_1, g_yyzz_xxxyzz_1, g_yyzz_xxxzz_1, g_yyzz_xxxzzz_1, g_yyzz_xxyyy_1, g_yyzz_xxyyyy_1, g_yyzz_xxyyyz_1, g_yyzz_xxyyz_1, g_yyzz_xxyyzz_1, g_yyzz_xxyzz_1, g_yyzz_xxyzzz_1, g_yyzz_xxzzz_1, g_yyzz_xxzzzz_1, g_yyzz_xyyyy_1, g_yyzz_xyyyyy_1, g_yyzz_xyyyyz_1, g_yyzz_xyyyz_1, g_yyzz_xyyyzz_1, g_yyzz_xyyzz_1, g_yyzz_xyyzzz_1, g_yyzz_xyzzz_1, g_yyzz_xyzzzz_1, g_yyzz_xzzzz_1, g_yyzz_xzzzzz_1, g_yyzz_yyyyy_1, g_yyzz_yyyyyy_1, g_yyzz_yyyyyz_1, g_yyzz_yyyyz_1, g_yyzz_yyyyzz_1, g_yyzz_yyyzz_1, g_yyzz_yyyzzz_1, g_yyzz_yyzzz_1, g_yyzz_yyzzzz_1, g_yyzz_yzzzz_1, g_yyzz_yzzzzz_1, g_yyzz_zzzzz_1, g_yyzz_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzz_xxxxxx_0[i] = 6.0 * g_yyzz_xxxxx_1[i] * fe_0 + g_yyzz_xxxxxx_1[i] * pa_x[i];

        g_xyyzz_xxxxxy_0[i] = 5.0 * g_yyzz_xxxxy_1[i] * fe_0 + g_yyzz_xxxxxy_1[i] * pa_x[i];

        g_xyyzz_xxxxxz_0[i] = 5.0 * g_yyzz_xxxxz_1[i] * fe_0 + g_yyzz_xxxxxz_1[i] * pa_x[i];

        g_xyyzz_xxxxyy_0[i] = 4.0 * g_yyzz_xxxyy_1[i] * fe_0 + g_yyzz_xxxxyy_1[i] * pa_x[i];

        g_xyyzz_xxxxyz_0[i] = 4.0 * g_yyzz_xxxyz_1[i] * fe_0 + g_yyzz_xxxxyz_1[i] * pa_x[i];

        g_xyyzz_xxxxzz_0[i] = 4.0 * g_yyzz_xxxzz_1[i] * fe_0 + g_yyzz_xxxxzz_1[i] * pa_x[i];

        g_xyyzz_xxxyyy_0[i] = 3.0 * g_yyzz_xxyyy_1[i] * fe_0 + g_yyzz_xxxyyy_1[i] * pa_x[i];

        g_xyyzz_xxxyyz_0[i] = 3.0 * g_yyzz_xxyyz_1[i] * fe_0 + g_yyzz_xxxyyz_1[i] * pa_x[i];

        g_xyyzz_xxxyzz_0[i] = 3.0 * g_yyzz_xxyzz_1[i] * fe_0 + g_yyzz_xxxyzz_1[i] * pa_x[i];

        g_xyyzz_xxxzzz_0[i] = 3.0 * g_yyzz_xxzzz_1[i] * fe_0 + g_yyzz_xxxzzz_1[i] * pa_x[i];

        g_xyyzz_xxyyyy_0[i] = 2.0 * g_yyzz_xyyyy_1[i] * fe_0 + g_yyzz_xxyyyy_1[i] * pa_x[i];

        g_xyyzz_xxyyyz_0[i] = 2.0 * g_yyzz_xyyyz_1[i] * fe_0 + g_yyzz_xxyyyz_1[i] * pa_x[i];

        g_xyyzz_xxyyzz_0[i] = 2.0 * g_yyzz_xyyzz_1[i] * fe_0 + g_yyzz_xxyyzz_1[i] * pa_x[i];

        g_xyyzz_xxyzzz_0[i] = 2.0 * g_yyzz_xyzzz_1[i] * fe_0 + g_yyzz_xxyzzz_1[i] * pa_x[i];

        g_xyyzz_xxzzzz_0[i] = 2.0 * g_yyzz_xzzzz_1[i] * fe_0 + g_yyzz_xxzzzz_1[i] * pa_x[i];

        g_xyyzz_xyyyyy_0[i] = g_yyzz_yyyyy_1[i] * fe_0 + g_yyzz_xyyyyy_1[i] * pa_x[i];

        g_xyyzz_xyyyyz_0[i] = g_yyzz_yyyyz_1[i] * fe_0 + g_yyzz_xyyyyz_1[i] * pa_x[i];

        g_xyyzz_xyyyzz_0[i] = g_yyzz_yyyzz_1[i] * fe_0 + g_yyzz_xyyyzz_1[i] * pa_x[i];

        g_xyyzz_xyyzzz_0[i] = g_yyzz_yyzzz_1[i] * fe_0 + g_yyzz_xyyzzz_1[i] * pa_x[i];

        g_xyyzz_xyzzzz_0[i] = g_yyzz_yzzzz_1[i] * fe_0 + g_yyzz_xyzzzz_1[i] * pa_x[i];

        g_xyyzz_xzzzzz_0[i] = g_yyzz_zzzzz_1[i] * fe_0 + g_yyzz_xzzzzz_1[i] * pa_x[i];

        g_xyyzz_yyyyyy_0[i] = g_yyzz_yyyyyy_1[i] * pa_x[i];

        g_xyyzz_yyyyyz_0[i] = g_yyzz_yyyyyz_1[i] * pa_x[i];

        g_xyyzz_yyyyzz_0[i] = g_yyzz_yyyyzz_1[i] * pa_x[i];

        g_xyyzz_yyyzzz_0[i] = g_yyzz_yyyzzz_1[i] * pa_x[i];

        g_xyyzz_yyzzzz_0[i] = g_yyzz_yyzzzz_1[i] * pa_x[i];

        g_xyyzz_yzzzzz_0[i] = g_yyzz_yzzzzz_1[i] * pa_x[i];

        g_xyyzz_zzzzzz_0[i] = g_yyzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 364-392 components of targeted buffer : HI

    auto g_xyzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 364);

    auto g_xyzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 365);

    auto g_xyzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 366);

    auto g_xyzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 367);

    auto g_xyzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 368);

    auto g_xyzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 369);

    auto g_xyzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 370);

    auto g_xyzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 371);

    auto g_xyzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 372);

    auto g_xyzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 373);

    auto g_xyzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 374);

    auto g_xyzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 375);

    auto g_xyzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 376);

    auto g_xyzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 377);

    auto g_xyzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 378);

    auto g_xyzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 379);

    auto g_xyzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 380);

    auto g_xyzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 381);

    auto g_xyzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 382);

    auto g_xyzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 383);

    auto g_xyzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 384);

    auto g_xyzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 385);

    auto g_xyzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 386);

    auto g_xyzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 387);

    auto g_xyzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 388);

    auto g_xyzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 389);

    auto g_xyzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 390);

    auto g_xyzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 391);

    #pragma omp simd aligned(g_xyzzz_xxxxxx_0, g_xyzzz_xxxxxy_0, g_xyzzz_xxxxxz_0, g_xyzzz_xxxxyy_0, g_xyzzz_xxxxyz_0, g_xyzzz_xxxxzz_0, g_xyzzz_xxxyyy_0, g_xyzzz_xxxyyz_0, g_xyzzz_xxxyzz_0, g_xyzzz_xxxzzz_0, g_xyzzz_xxyyyy_0, g_xyzzz_xxyyyz_0, g_xyzzz_xxyyzz_0, g_xyzzz_xxyzzz_0, g_xyzzz_xxzzzz_0, g_xyzzz_xyyyyy_0, g_xyzzz_xyyyyz_0, g_xyzzz_xyyyzz_0, g_xyzzz_xyyzzz_0, g_xyzzz_xyzzzz_0, g_xyzzz_xzzzzz_0, g_xyzzz_yyyyyy_0, g_xyzzz_yyyyyz_0, g_xyzzz_yyyyzz_0, g_xyzzz_yyyzzz_0, g_xyzzz_yyzzzz_0, g_xyzzz_yzzzzz_0, g_xyzzz_zzzzzz_0, g_xzzz_xxxxxx_1, g_xzzz_xxxxxz_1, g_xzzz_xxxxzz_1, g_xzzz_xxxzzz_1, g_xzzz_xxzzzz_1, g_xzzz_xzzzzz_1, g_yzzz_xxxxxy_1, g_yzzz_xxxxy_1, g_yzzz_xxxxyy_1, g_yzzz_xxxxyz_1, g_yzzz_xxxyy_1, g_yzzz_xxxyyy_1, g_yzzz_xxxyyz_1, g_yzzz_xxxyz_1, g_yzzz_xxxyzz_1, g_yzzz_xxyyy_1, g_yzzz_xxyyyy_1, g_yzzz_xxyyyz_1, g_yzzz_xxyyz_1, g_yzzz_xxyyzz_1, g_yzzz_xxyzz_1, g_yzzz_xxyzzz_1, g_yzzz_xyyyy_1, g_yzzz_xyyyyy_1, g_yzzz_xyyyyz_1, g_yzzz_xyyyz_1, g_yzzz_xyyyzz_1, g_yzzz_xyyzz_1, g_yzzz_xyyzzz_1, g_yzzz_xyzzz_1, g_yzzz_xyzzzz_1, g_yzzz_yyyyy_1, g_yzzz_yyyyyy_1, g_yzzz_yyyyyz_1, g_yzzz_yyyyz_1, g_yzzz_yyyyzz_1, g_yzzz_yyyzz_1, g_yzzz_yyyzzz_1, g_yzzz_yyzzz_1, g_yzzz_yyzzzz_1, g_yzzz_yzzzz_1, g_yzzz_yzzzzz_1, g_yzzz_zzzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzz_xxxxxx_0[i] = g_xzzz_xxxxxx_1[i] * pa_y[i];

        g_xyzzz_xxxxxy_0[i] = 5.0 * g_yzzz_xxxxy_1[i] * fe_0 + g_yzzz_xxxxxy_1[i] * pa_x[i];

        g_xyzzz_xxxxxz_0[i] = g_xzzz_xxxxxz_1[i] * pa_y[i];

        g_xyzzz_xxxxyy_0[i] = 4.0 * g_yzzz_xxxyy_1[i] * fe_0 + g_yzzz_xxxxyy_1[i] * pa_x[i];

        g_xyzzz_xxxxyz_0[i] = 4.0 * g_yzzz_xxxyz_1[i] * fe_0 + g_yzzz_xxxxyz_1[i] * pa_x[i];

        g_xyzzz_xxxxzz_0[i] = g_xzzz_xxxxzz_1[i] * pa_y[i];

        g_xyzzz_xxxyyy_0[i] = 3.0 * g_yzzz_xxyyy_1[i] * fe_0 + g_yzzz_xxxyyy_1[i] * pa_x[i];

        g_xyzzz_xxxyyz_0[i] = 3.0 * g_yzzz_xxyyz_1[i] * fe_0 + g_yzzz_xxxyyz_1[i] * pa_x[i];

        g_xyzzz_xxxyzz_0[i] = 3.0 * g_yzzz_xxyzz_1[i] * fe_0 + g_yzzz_xxxyzz_1[i] * pa_x[i];

        g_xyzzz_xxxzzz_0[i] = g_xzzz_xxxzzz_1[i] * pa_y[i];

        g_xyzzz_xxyyyy_0[i] = 2.0 * g_yzzz_xyyyy_1[i] * fe_0 + g_yzzz_xxyyyy_1[i] * pa_x[i];

        g_xyzzz_xxyyyz_0[i] = 2.0 * g_yzzz_xyyyz_1[i] * fe_0 + g_yzzz_xxyyyz_1[i] * pa_x[i];

        g_xyzzz_xxyyzz_0[i] = 2.0 * g_yzzz_xyyzz_1[i] * fe_0 + g_yzzz_xxyyzz_1[i] * pa_x[i];

        g_xyzzz_xxyzzz_0[i] = 2.0 * g_yzzz_xyzzz_1[i] * fe_0 + g_yzzz_xxyzzz_1[i] * pa_x[i];

        g_xyzzz_xxzzzz_0[i] = g_xzzz_xxzzzz_1[i] * pa_y[i];

        g_xyzzz_xyyyyy_0[i] = g_yzzz_yyyyy_1[i] * fe_0 + g_yzzz_xyyyyy_1[i] * pa_x[i];

        g_xyzzz_xyyyyz_0[i] = g_yzzz_yyyyz_1[i] * fe_0 + g_yzzz_xyyyyz_1[i] * pa_x[i];

        g_xyzzz_xyyyzz_0[i] = g_yzzz_yyyzz_1[i] * fe_0 + g_yzzz_xyyyzz_1[i] * pa_x[i];

        g_xyzzz_xyyzzz_0[i] = g_yzzz_yyzzz_1[i] * fe_0 + g_yzzz_xyyzzz_1[i] * pa_x[i];

        g_xyzzz_xyzzzz_0[i] = g_yzzz_yzzzz_1[i] * fe_0 + g_yzzz_xyzzzz_1[i] * pa_x[i];

        g_xyzzz_xzzzzz_0[i] = g_xzzz_xzzzzz_1[i] * pa_y[i];

        g_xyzzz_yyyyyy_0[i] = g_yzzz_yyyyyy_1[i] * pa_x[i];

        g_xyzzz_yyyyyz_0[i] = g_yzzz_yyyyyz_1[i] * pa_x[i];

        g_xyzzz_yyyyzz_0[i] = g_yzzz_yyyyzz_1[i] * pa_x[i];

        g_xyzzz_yyyzzz_0[i] = g_yzzz_yyyzzz_1[i] * pa_x[i];

        g_xyzzz_yyzzzz_0[i] = g_yzzz_yyzzzz_1[i] * pa_x[i];

        g_xyzzz_yzzzzz_0[i] = g_yzzz_yzzzzz_1[i] * pa_x[i];

        g_xyzzz_zzzzzz_0[i] = g_yzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 392-420 components of targeted buffer : HI

    auto g_xzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 392);

    auto g_xzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 393);

    auto g_xzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 394);

    auto g_xzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 395);

    auto g_xzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 396);

    auto g_xzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 397);

    auto g_xzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 398);

    auto g_xzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 399);

    auto g_xzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 400);

    auto g_xzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 401);

    auto g_xzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 402);

    auto g_xzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 403);

    auto g_xzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 404);

    auto g_xzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 405);

    auto g_xzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 406);

    auto g_xzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 407);

    auto g_xzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 408);

    auto g_xzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 409);

    auto g_xzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 410);

    auto g_xzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 411);

    auto g_xzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 412);

    auto g_xzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 413);

    auto g_xzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 414);

    auto g_xzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 415);

    auto g_xzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 416);

    auto g_xzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 417);

    auto g_xzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 418);

    auto g_xzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 419);

    #pragma omp simd aligned(g_xzzzz_xxxxxx_0, g_xzzzz_xxxxxy_0, g_xzzzz_xxxxxz_0, g_xzzzz_xxxxyy_0, g_xzzzz_xxxxyz_0, g_xzzzz_xxxxzz_0, g_xzzzz_xxxyyy_0, g_xzzzz_xxxyyz_0, g_xzzzz_xxxyzz_0, g_xzzzz_xxxzzz_0, g_xzzzz_xxyyyy_0, g_xzzzz_xxyyyz_0, g_xzzzz_xxyyzz_0, g_xzzzz_xxyzzz_0, g_xzzzz_xxzzzz_0, g_xzzzz_xyyyyy_0, g_xzzzz_xyyyyz_0, g_xzzzz_xyyyzz_0, g_xzzzz_xyyzzz_0, g_xzzzz_xyzzzz_0, g_xzzzz_xzzzzz_0, g_xzzzz_yyyyyy_0, g_xzzzz_yyyyyz_0, g_xzzzz_yyyyzz_0, g_xzzzz_yyyzzz_0, g_xzzzz_yyzzzz_0, g_xzzzz_yzzzzz_0, g_xzzzz_zzzzzz_0, g_zzzz_xxxxx_1, g_zzzz_xxxxxx_1, g_zzzz_xxxxxy_1, g_zzzz_xxxxxz_1, g_zzzz_xxxxy_1, g_zzzz_xxxxyy_1, g_zzzz_xxxxyz_1, g_zzzz_xxxxz_1, g_zzzz_xxxxzz_1, g_zzzz_xxxyy_1, g_zzzz_xxxyyy_1, g_zzzz_xxxyyz_1, g_zzzz_xxxyz_1, g_zzzz_xxxyzz_1, g_zzzz_xxxzz_1, g_zzzz_xxxzzz_1, g_zzzz_xxyyy_1, g_zzzz_xxyyyy_1, g_zzzz_xxyyyz_1, g_zzzz_xxyyz_1, g_zzzz_xxyyzz_1, g_zzzz_xxyzz_1, g_zzzz_xxyzzz_1, g_zzzz_xxzzz_1, g_zzzz_xxzzzz_1, g_zzzz_xyyyy_1, g_zzzz_xyyyyy_1, g_zzzz_xyyyyz_1, g_zzzz_xyyyz_1, g_zzzz_xyyyzz_1, g_zzzz_xyyzz_1, g_zzzz_xyyzzz_1, g_zzzz_xyzzz_1, g_zzzz_xyzzzz_1, g_zzzz_xzzzz_1, g_zzzz_xzzzzz_1, g_zzzz_yyyyy_1, g_zzzz_yyyyyy_1, g_zzzz_yyyyyz_1, g_zzzz_yyyyz_1, g_zzzz_yyyyzz_1, g_zzzz_yyyzz_1, g_zzzz_yyyzzz_1, g_zzzz_yyzzz_1, g_zzzz_yyzzzz_1, g_zzzz_yzzzz_1, g_zzzz_yzzzzz_1, g_zzzz_zzzzz_1, g_zzzz_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzz_xxxxxx_0[i] = 6.0 * g_zzzz_xxxxx_1[i] * fe_0 + g_zzzz_xxxxxx_1[i] * pa_x[i];

        g_xzzzz_xxxxxy_0[i] = 5.0 * g_zzzz_xxxxy_1[i] * fe_0 + g_zzzz_xxxxxy_1[i] * pa_x[i];

        g_xzzzz_xxxxxz_0[i] = 5.0 * g_zzzz_xxxxz_1[i] * fe_0 + g_zzzz_xxxxxz_1[i] * pa_x[i];

        g_xzzzz_xxxxyy_0[i] = 4.0 * g_zzzz_xxxyy_1[i] * fe_0 + g_zzzz_xxxxyy_1[i] * pa_x[i];

        g_xzzzz_xxxxyz_0[i] = 4.0 * g_zzzz_xxxyz_1[i] * fe_0 + g_zzzz_xxxxyz_1[i] * pa_x[i];

        g_xzzzz_xxxxzz_0[i] = 4.0 * g_zzzz_xxxzz_1[i] * fe_0 + g_zzzz_xxxxzz_1[i] * pa_x[i];

        g_xzzzz_xxxyyy_0[i] = 3.0 * g_zzzz_xxyyy_1[i] * fe_0 + g_zzzz_xxxyyy_1[i] * pa_x[i];

        g_xzzzz_xxxyyz_0[i] = 3.0 * g_zzzz_xxyyz_1[i] * fe_0 + g_zzzz_xxxyyz_1[i] * pa_x[i];

        g_xzzzz_xxxyzz_0[i] = 3.0 * g_zzzz_xxyzz_1[i] * fe_0 + g_zzzz_xxxyzz_1[i] * pa_x[i];

        g_xzzzz_xxxzzz_0[i] = 3.0 * g_zzzz_xxzzz_1[i] * fe_0 + g_zzzz_xxxzzz_1[i] * pa_x[i];

        g_xzzzz_xxyyyy_0[i] = 2.0 * g_zzzz_xyyyy_1[i] * fe_0 + g_zzzz_xxyyyy_1[i] * pa_x[i];

        g_xzzzz_xxyyyz_0[i] = 2.0 * g_zzzz_xyyyz_1[i] * fe_0 + g_zzzz_xxyyyz_1[i] * pa_x[i];

        g_xzzzz_xxyyzz_0[i] = 2.0 * g_zzzz_xyyzz_1[i] * fe_0 + g_zzzz_xxyyzz_1[i] * pa_x[i];

        g_xzzzz_xxyzzz_0[i] = 2.0 * g_zzzz_xyzzz_1[i] * fe_0 + g_zzzz_xxyzzz_1[i] * pa_x[i];

        g_xzzzz_xxzzzz_0[i] = 2.0 * g_zzzz_xzzzz_1[i] * fe_0 + g_zzzz_xxzzzz_1[i] * pa_x[i];

        g_xzzzz_xyyyyy_0[i] = g_zzzz_yyyyy_1[i] * fe_0 + g_zzzz_xyyyyy_1[i] * pa_x[i];

        g_xzzzz_xyyyyz_0[i] = g_zzzz_yyyyz_1[i] * fe_0 + g_zzzz_xyyyyz_1[i] * pa_x[i];

        g_xzzzz_xyyyzz_0[i] = g_zzzz_yyyzz_1[i] * fe_0 + g_zzzz_xyyyzz_1[i] * pa_x[i];

        g_xzzzz_xyyzzz_0[i] = g_zzzz_yyzzz_1[i] * fe_0 + g_zzzz_xyyzzz_1[i] * pa_x[i];

        g_xzzzz_xyzzzz_0[i] = g_zzzz_yzzzz_1[i] * fe_0 + g_zzzz_xyzzzz_1[i] * pa_x[i];

        g_xzzzz_xzzzzz_0[i] = g_zzzz_zzzzz_1[i] * fe_0 + g_zzzz_xzzzzz_1[i] * pa_x[i];

        g_xzzzz_yyyyyy_0[i] = g_zzzz_yyyyyy_1[i] * pa_x[i];

        g_xzzzz_yyyyyz_0[i] = g_zzzz_yyyyyz_1[i] * pa_x[i];

        g_xzzzz_yyyyzz_0[i] = g_zzzz_yyyyzz_1[i] * pa_x[i];

        g_xzzzz_yyyzzz_0[i] = g_zzzz_yyyzzz_1[i] * pa_x[i];

        g_xzzzz_yyzzzz_0[i] = g_zzzz_yyzzzz_1[i] * pa_x[i];

        g_xzzzz_yzzzzz_0[i] = g_zzzz_yzzzzz_1[i] * pa_x[i];

        g_xzzzz_zzzzzz_0[i] = g_zzzz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 420-448 components of targeted buffer : HI

    auto g_yyyyy_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 420);

    auto g_yyyyy_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 421);

    auto g_yyyyy_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 422);

    auto g_yyyyy_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 423);

    auto g_yyyyy_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 424);

    auto g_yyyyy_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 425);

    auto g_yyyyy_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 426);

    auto g_yyyyy_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 427);

    auto g_yyyyy_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 428);

    auto g_yyyyy_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 429);

    auto g_yyyyy_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 430);

    auto g_yyyyy_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 431);

    auto g_yyyyy_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 432);

    auto g_yyyyy_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 433);

    auto g_yyyyy_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 434);

    auto g_yyyyy_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 435);

    auto g_yyyyy_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 436);

    auto g_yyyyy_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 437);

    auto g_yyyyy_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 438);

    auto g_yyyyy_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 439);

    auto g_yyyyy_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 440);

    auto g_yyyyy_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 441);

    auto g_yyyyy_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 442);

    auto g_yyyyy_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 443);

    auto g_yyyyy_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 444);

    auto g_yyyyy_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 445);

    auto g_yyyyy_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 446);

    auto g_yyyyy_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 447);

    #pragma omp simd aligned(g_yyy_xxxxxx_0, g_yyy_xxxxxx_1, g_yyy_xxxxxy_0, g_yyy_xxxxxy_1, g_yyy_xxxxxz_0, g_yyy_xxxxxz_1, g_yyy_xxxxyy_0, g_yyy_xxxxyy_1, g_yyy_xxxxyz_0, g_yyy_xxxxyz_1, g_yyy_xxxxzz_0, g_yyy_xxxxzz_1, g_yyy_xxxyyy_0, g_yyy_xxxyyy_1, g_yyy_xxxyyz_0, g_yyy_xxxyyz_1, g_yyy_xxxyzz_0, g_yyy_xxxyzz_1, g_yyy_xxxzzz_0, g_yyy_xxxzzz_1, g_yyy_xxyyyy_0, g_yyy_xxyyyy_1, g_yyy_xxyyyz_0, g_yyy_xxyyyz_1, g_yyy_xxyyzz_0, g_yyy_xxyyzz_1, g_yyy_xxyzzz_0, g_yyy_xxyzzz_1, g_yyy_xxzzzz_0, g_yyy_xxzzzz_1, g_yyy_xyyyyy_0, g_yyy_xyyyyy_1, g_yyy_xyyyyz_0, g_yyy_xyyyyz_1, g_yyy_xyyyzz_0, g_yyy_xyyyzz_1, g_yyy_xyyzzz_0, g_yyy_xyyzzz_1, g_yyy_xyzzzz_0, g_yyy_xyzzzz_1, g_yyy_xzzzzz_0, g_yyy_xzzzzz_1, g_yyy_yyyyyy_0, g_yyy_yyyyyy_1, g_yyy_yyyyyz_0, g_yyy_yyyyyz_1, g_yyy_yyyyzz_0, g_yyy_yyyyzz_1, g_yyy_yyyzzz_0, g_yyy_yyyzzz_1, g_yyy_yyzzzz_0, g_yyy_yyzzzz_1, g_yyy_yzzzzz_0, g_yyy_yzzzzz_1, g_yyy_zzzzzz_0, g_yyy_zzzzzz_1, g_yyyy_xxxxx_1, g_yyyy_xxxxxx_1, g_yyyy_xxxxxy_1, g_yyyy_xxxxxz_1, g_yyyy_xxxxy_1, g_yyyy_xxxxyy_1, g_yyyy_xxxxyz_1, g_yyyy_xxxxz_1, g_yyyy_xxxxzz_1, g_yyyy_xxxyy_1, g_yyyy_xxxyyy_1, g_yyyy_xxxyyz_1, g_yyyy_xxxyz_1, g_yyyy_xxxyzz_1, g_yyyy_xxxzz_1, g_yyyy_xxxzzz_1, g_yyyy_xxyyy_1, g_yyyy_xxyyyy_1, g_yyyy_xxyyyz_1, g_yyyy_xxyyz_1, g_yyyy_xxyyzz_1, g_yyyy_xxyzz_1, g_yyyy_xxyzzz_1, g_yyyy_xxzzz_1, g_yyyy_xxzzzz_1, g_yyyy_xyyyy_1, g_yyyy_xyyyyy_1, g_yyyy_xyyyyz_1, g_yyyy_xyyyz_1, g_yyyy_xyyyzz_1, g_yyyy_xyyzz_1, g_yyyy_xyyzzz_1, g_yyyy_xyzzz_1, g_yyyy_xyzzzz_1, g_yyyy_xzzzz_1, g_yyyy_xzzzzz_1, g_yyyy_yyyyy_1, g_yyyy_yyyyyy_1, g_yyyy_yyyyyz_1, g_yyyy_yyyyz_1, g_yyyy_yyyyzz_1, g_yyyy_yyyzz_1, g_yyyy_yyyzzz_1, g_yyyy_yyzzz_1, g_yyyy_yyzzzz_1, g_yyyy_yzzzz_1, g_yyyy_yzzzzz_1, g_yyyy_zzzzz_1, g_yyyy_zzzzzz_1, g_yyyyy_xxxxxx_0, g_yyyyy_xxxxxy_0, g_yyyyy_xxxxxz_0, g_yyyyy_xxxxyy_0, g_yyyyy_xxxxyz_0, g_yyyyy_xxxxzz_0, g_yyyyy_xxxyyy_0, g_yyyyy_xxxyyz_0, g_yyyyy_xxxyzz_0, g_yyyyy_xxxzzz_0, g_yyyyy_xxyyyy_0, g_yyyyy_xxyyyz_0, g_yyyyy_xxyyzz_0, g_yyyyy_xxyzzz_0, g_yyyyy_xxzzzz_0, g_yyyyy_xyyyyy_0, g_yyyyy_xyyyyz_0, g_yyyyy_xyyyzz_0, g_yyyyy_xyyzzz_0, g_yyyyy_xyzzzz_0, g_yyyyy_xzzzzz_0, g_yyyyy_yyyyyy_0, g_yyyyy_yyyyyz_0, g_yyyyy_yyyyzz_0, g_yyyyy_yyyzzz_0, g_yyyyy_yyzzzz_0, g_yyyyy_yzzzzz_0, g_yyyyy_zzzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyy_xxxxxx_0[i] = 4.0 * g_yyy_xxxxxx_0[i] * fbe_0 - 4.0 * g_yyy_xxxxxx_1[i] * fz_be_0 + g_yyyy_xxxxxx_1[i] * pa_y[i];

        g_yyyyy_xxxxxy_0[i] = 4.0 * g_yyy_xxxxxy_0[i] * fbe_0 - 4.0 * g_yyy_xxxxxy_1[i] * fz_be_0 + g_yyyy_xxxxx_1[i] * fe_0 + g_yyyy_xxxxxy_1[i] * pa_y[i];

        g_yyyyy_xxxxxz_0[i] = 4.0 * g_yyy_xxxxxz_0[i] * fbe_0 - 4.0 * g_yyy_xxxxxz_1[i] * fz_be_0 + g_yyyy_xxxxxz_1[i] * pa_y[i];

        g_yyyyy_xxxxyy_0[i] = 4.0 * g_yyy_xxxxyy_0[i] * fbe_0 - 4.0 * g_yyy_xxxxyy_1[i] * fz_be_0 + 2.0 * g_yyyy_xxxxy_1[i] * fe_0 + g_yyyy_xxxxyy_1[i] * pa_y[i];

        g_yyyyy_xxxxyz_0[i] = 4.0 * g_yyy_xxxxyz_0[i] * fbe_0 - 4.0 * g_yyy_xxxxyz_1[i] * fz_be_0 + g_yyyy_xxxxz_1[i] * fe_0 + g_yyyy_xxxxyz_1[i] * pa_y[i];

        g_yyyyy_xxxxzz_0[i] = 4.0 * g_yyy_xxxxzz_0[i] * fbe_0 - 4.0 * g_yyy_xxxxzz_1[i] * fz_be_0 + g_yyyy_xxxxzz_1[i] * pa_y[i];

        g_yyyyy_xxxyyy_0[i] = 4.0 * g_yyy_xxxyyy_0[i] * fbe_0 - 4.0 * g_yyy_xxxyyy_1[i] * fz_be_0 + 3.0 * g_yyyy_xxxyy_1[i] * fe_0 + g_yyyy_xxxyyy_1[i] * pa_y[i];

        g_yyyyy_xxxyyz_0[i] = 4.0 * g_yyy_xxxyyz_0[i] * fbe_0 - 4.0 * g_yyy_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyyy_xxxyz_1[i] * fe_0 + g_yyyy_xxxyyz_1[i] * pa_y[i];

        g_yyyyy_xxxyzz_0[i] = 4.0 * g_yyy_xxxyzz_0[i] * fbe_0 - 4.0 * g_yyy_xxxyzz_1[i] * fz_be_0 + g_yyyy_xxxzz_1[i] * fe_0 + g_yyyy_xxxyzz_1[i] * pa_y[i];

        g_yyyyy_xxxzzz_0[i] = 4.0 * g_yyy_xxxzzz_0[i] * fbe_0 - 4.0 * g_yyy_xxxzzz_1[i] * fz_be_0 + g_yyyy_xxxzzz_1[i] * pa_y[i];

        g_yyyyy_xxyyyy_0[i] = 4.0 * g_yyy_xxyyyy_0[i] * fbe_0 - 4.0 * g_yyy_xxyyyy_1[i] * fz_be_0 + 4.0 * g_yyyy_xxyyy_1[i] * fe_0 + g_yyyy_xxyyyy_1[i] * pa_y[i];

        g_yyyyy_xxyyyz_0[i] = 4.0 * g_yyy_xxyyyz_0[i] * fbe_0 - 4.0 * g_yyy_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyyy_xxyyz_1[i] * fe_0 + g_yyyy_xxyyyz_1[i] * pa_y[i];

        g_yyyyy_xxyyzz_0[i] = 4.0 * g_yyy_xxyyzz_0[i] * fbe_0 - 4.0 * g_yyy_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyyy_xxyzz_1[i] * fe_0 + g_yyyy_xxyyzz_1[i] * pa_y[i];

        g_yyyyy_xxyzzz_0[i] = 4.0 * g_yyy_xxyzzz_0[i] * fbe_0 - 4.0 * g_yyy_xxyzzz_1[i] * fz_be_0 + g_yyyy_xxzzz_1[i] * fe_0 + g_yyyy_xxyzzz_1[i] * pa_y[i];

        g_yyyyy_xxzzzz_0[i] = 4.0 * g_yyy_xxzzzz_0[i] * fbe_0 - 4.0 * g_yyy_xxzzzz_1[i] * fz_be_0 + g_yyyy_xxzzzz_1[i] * pa_y[i];

        g_yyyyy_xyyyyy_0[i] = 4.0 * g_yyy_xyyyyy_0[i] * fbe_0 - 4.0 * g_yyy_xyyyyy_1[i] * fz_be_0 + 5.0 * g_yyyy_xyyyy_1[i] * fe_0 + g_yyyy_xyyyyy_1[i] * pa_y[i];

        g_yyyyy_xyyyyz_0[i] = 4.0 * g_yyy_xyyyyz_0[i] * fbe_0 - 4.0 * g_yyy_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyyy_xyyyz_1[i] * fe_0 + g_yyyy_xyyyyz_1[i] * pa_y[i];

        g_yyyyy_xyyyzz_0[i] = 4.0 * g_yyy_xyyyzz_0[i] * fbe_0 - 4.0 * g_yyy_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyyy_xyyzz_1[i] * fe_0 + g_yyyy_xyyyzz_1[i] * pa_y[i];

        g_yyyyy_xyyzzz_0[i] = 4.0 * g_yyy_xyyzzz_0[i] * fbe_0 - 4.0 * g_yyy_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_xyzzz_1[i] * fe_0 + g_yyyy_xyyzzz_1[i] * pa_y[i];

        g_yyyyy_xyzzzz_0[i] = 4.0 * g_yyy_xyzzzz_0[i] * fbe_0 - 4.0 * g_yyy_xyzzzz_1[i] * fz_be_0 + g_yyyy_xzzzz_1[i] * fe_0 + g_yyyy_xyzzzz_1[i] * pa_y[i];

        g_yyyyy_xzzzzz_0[i] = 4.0 * g_yyy_xzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_xzzzzz_1[i] * fz_be_0 + g_yyyy_xzzzzz_1[i] * pa_y[i];

        g_yyyyy_yyyyyy_0[i] = 4.0 * g_yyy_yyyyyy_0[i] * fbe_0 - 4.0 * g_yyy_yyyyyy_1[i] * fz_be_0 + 6.0 * g_yyyy_yyyyy_1[i] * fe_0 + g_yyyy_yyyyyy_1[i] * pa_y[i];

        g_yyyyy_yyyyyz_0[i] = 4.0 * g_yyy_yyyyyz_0[i] * fbe_0 - 4.0 * g_yyy_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyyy_yyyyz_1[i] * fe_0 + g_yyyy_yyyyyz_1[i] * pa_y[i];

        g_yyyyy_yyyyzz_0[i] = 4.0 * g_yyy_yyyyzz_0[i] * fbe_0 - 4.0 * g_yyy_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyyy_yyyzz_1[i] * fe_0 + g_yyyy_yyyyzz_1[i] * pa_y[i];

        g_yyyyy_yyyzzz_0[i] = 4.0 * g_yyy_yyyzzz_0[i] * fbe_0 - 4.0 * g_yyy_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyyy_yyzzz_1[i] * fe_0 + g_yyyy_yyyzzz_1[i] * pa_y[i];

        g_yyyyy_yyzzzz_0[i] = 4.0 * g_yyy_yyzzzz_0[i] * fbe_0 - 4.0 * g_yyy_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyyy_yzzzz_1[i] * fe_0 + g_yyyy_yyzzzz_1[i] * pa_y[i];

        g_yyyyy_yzzzzz_0[i] = 4.0 * g_yyy_yzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_yzzzzz_1[i] * fz_be_0 + g_yyyy_zzzzz_1[i] * fe_0 + g_yyyy_yzzzzz_1[i] * pa_y[i];

        g_yyyyy_zzzzzz_0[i] = 4.0 * g_yyy_zzzzzz_0[i] * fbe_0 - 4.0 * g_yyy_zzzzzz_1[i] * fz_be_0 + g_yyyy_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 448-476 components of targeted buffer : HI

    auto g_yyyyz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 448);

    auto g_yyyyz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 449);

    auto g_yyyyz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 450);

    auto g_yyyyz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 451);

    auto g_yyyyz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 452);

    auto g_yyyyz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 453);

    auto g_yyyyz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 454);

    auto g_yyyyz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 455);

    auto g_yyyyz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 456);

    auto g_yyyyz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 457);

    auto g_yyyyz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 458);

    auto g_yyyyz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 459);

    auto g_yyyyz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 460);

    auto g_yyyyz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 461);

    auto g_yyyyz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 462);

    auto g_yyyyz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 463);

    auto g_yyyyz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 464);

    auto g_yyyyz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 465);

    auto g_yyyyz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 466);

    auto g_yyyyz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 467);

    auto g_yyyyz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 468);

    auto g_yyyyz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 469);

    auto g_yyyyz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 470);

    auto g_yyyyz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 471);

    auto g_yyyyz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 472);

    auto g_yyyyz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 473);

    auto g_yyyyz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 474);

    auto g_yyyyz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 475);

    #pragma omp simd aligned(g_yyyy_xxxxx_1, g_yyyy_xxxxxx_1, g_yyyy_xxxxxy_1, g_yyyy_xxxxxz_1, g_yyyy_xxxxy_1, g_yyyy_xxxxyy_1, g_yyyy_xxxxyz_1, g_yyyy_xxxxz_1, g_yyyy_xxxxzz_1, g_yyyy_xxxyy_1, g_yyyy_xxxyyy_1, g_yyyy_xxxyyz_1, g_yyyy_xxxyz_1, g_yyyy_xxxyzz_1, g_yyyy_xxxzz_1, g_yyyy_xxxzzz_1, g_yyyy_xxyyy_1, g_yyyy_xxyyyy_1, g_yyyy_xxyyyz_1, g_yyyy_xxyyz_1, g_yyyy_xxyyzz_1, g_yyyy_xxyzz_1, g_yyyy_xxyzzz_1, g_yyyy_xxzzz_1, g_yyyy_xxzzzz_1, g_yyyy_xyyyy_1, g_yyyy_xyyyyy_1, g_yyyy_xyyyyz_1, g_yyyy_xyyyz_1, g_yyyy_xyyyzz_1, g_yyyy_xyyzz_1, g_yyyy_xyyzzz_1, g_yyyy_xyzzz_1, g_yyyy_xyzzzz_1, g_yyyy_xzzzz_1, g_yyyy_xzzzzz_1, g_yyyy_yyyyy_1, g_yyyy_yyyyyy_1, g_yyyy_yyyyyz_1, g_yyyy_yyyyz_1, g_yyyy_yyyyzz_1, g_yyyy_yyyzz_1, g_yyyy_yyyzzz_1, g_yyyy_yyzzz_1, g_yyyy_yyzzzz_1, g_yyyy_yzzzz_1, g_yyyy_yzzzzz_1, g_yyyy_zzzzz_1, g_yyyy_zzzzzz_1, g_yyyyz_xxxxxx_0, g_yyyyz_xxxxxy_0, g_yyyyz_xxxxxz_0, g_yyyyz_xxxxyy_0, g_yyyyz_xxxxyz_0, g_yyyyz_xxxxzz_0, g_yyyyz_xxxyyy_0, g_yyyyz_xxxyyz_0, g_yyyyz_xxxyzz_0, g_yyyyz_xxxzzz_0, g_yyyyz_xxyyyy_0, g_yyyyz_xxyyyz_0, g_yyyyz_xxyyzz_0, g_yyyyz_xxyzzz_0, g_yyyyz_xxzzzz_0, g_yyyyz_xyyyyy_0, g_yyyyz_xyyyyz_0, g_yyyyz_xyyyzz_0, g_yyyyz_xyyzzz_0, g_yyyyz_xyzzzz_0, g_yyyyz_xzzzzz_0, g_yyyyz_yyyyyy_0, g_yyyyz_yyyyyz_0, g_yyyyz_yyyyzz_0, g_yyyyz_yyyzzz_0, g_yyyyz_yyzzzz_0, g_yyyyz_yzzzzz_0, g_yyyyz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyz_xxxxxx_0[i] = g_yyyy_xxxxxx_1[i] * pa_z[i];

        g_yyyyz_xxxxxy_0[i] = g_yyyy_xxxxxy_1[i] * pa_z[i];

        g_yyyyz_xxxxxz_0[i] = g_yyyy_xxxxx_1[i] * fe_0 + g_yyyy_xxxxxz_1[i] * pa_z[i];

        g_yyyyz_xxxxyy_0[i] = g_yyyy_xxxxyy_1[i] * pa_z[i];

        g_yyyyz_xxxxyz_0[i] = g_yyyy_xxxxy_1[i] * fe_0 + g_yyyy_xxxxyz_1[i] * pa_z[i];

        g_yyyyz_xxxxzz_0[i] = 2.0 * g_yyyy_xxxxz_1[i] * fe_0 + g_yyyy_xxxxzz_1[i] * pa_z[i];

        g_yyyyz_xxxyyy_0[i] = g_yyyy_xxxyyy_1[i] * pa_z[i];

        g_yyyyz_xxxyyz_0[i] = g_yyyy_xxxyy_1[i] * fe_0 + g_yyyy_xxxyyz_1[i] * pa_z[i];

        g_yyyyz_xxxyzz_0[i] = 2.0 * g_yyyy_xxxyz_1[i] * fe_0 + g_yyyy_xxxyzz_1[i] * pa_z[i];

        g_yyyyz_xxxzzz_0[i] = 3.0 * g_yyyy_xxxzz_1[i] * fe_0 + g_yyyy_xxxzzz_1[i] * pa_z[i];

        g_yyyyz_xxyyyy_0[i] = g_yyyy_xxyyyy_1[i] * pa_z[i];

        g_yyyyz_xxyyyz_0[i] = g_yyyy_xxyyy_1[i] * fe_0 + g_yyyy_xxyyyz_1[i] * pa_z[i];

        g_yyyyz_xxyyzz_0[i] = 2.0 * g_yyyy_xxyyz_1[i] * fe_0 + g_yyyy_xxyyzz_1[i] * pa_z[i];

        g_yyyyz_xxyzzz_0[i] = 3.0 * g_yyyy_xxyzz_1[i] * fe_0 + g_yyyy_xxyzzz_1[i] * pa_z[i];

        g_yyyyz_xxzzzz_0[i] = 4.0 * g_yyyy_xxzzz_1[i] * fe_0 + g_yyyy_xxzzzz_1[i] * pa_z[i];

        g_yyyyz_xyyyyy_0[i] = g_yyyy_xyyyyy_1[i] * pa_z[i];

        g_yyyyz_xyyyyz_0[i] = g_yyyy_xyyyy_1[i] * fe_0 + g_yyyy_xyyyyz_1[i] * pa_z[i];

        g_yyyyz_xyyyzz_0[i] = 2.0 * g_yyyy_xyyyz_1[i] * fe_0 + g_yyyy_xyyyzz_1[i] * pa_z[i];

        g_yyyyz_xyyzzz_0[i] = 3.0 * g_yyyy_xyyzz_1[i] * fe_0 + g_yyyy_xyyzzz_1[i] * pa_z[i];

        g_yyyyz_xyzzzz_0[i] = 4.0 * g_yyyy_xyzzz_1[i] * fe_0 + g_yyyy_xyzzzz_1[i] * pa_z[i];

        g_yyyyz_xzzzzz_0[i] = 5.0 * g_yyyy_xzzzz_1[i] * fe_0 + g_yyyy_xzzzzz_1[i] * pa_z[i];

        g_yyyyz_yyyyyy_0[i] = g_yyyy_yyyyyy_1[i] * pa_z[i];

        g_yyyyz_yyyyyz_0[i] = g_yyyy_yyyyy_1[i] * fe_0 + g_yyyy_yyyyyz_1[i] * pa_z[i];

        g_yyyyz_yyyyzz_0[i] = 2.0 * g_yyyy_yyyyz_1[i] * fe_0 + g_yyyy_yyyyzz_1[i] * pa_z[i];

        g_yyyyz_yyyzzz_0[i] = 3.0 * g_yyyy_yyyzz_1[i] * fe_0 + g_yyyy_yyyzzz_1[i] * pa_z[i];

        g_yyyyz_yyzzzz_0[i] = 4.0 * g_yyyy_yyzzz_1[i] * fe_0 + g_yyyy_yyzzzz_1[i] * pa_z[i];

        g_yyyyz_yzzzzz_0[i] = 5.0 * g_yyyy_yzzzz_1[i] * fe_0 + g_yyyy_yzzzzz_1[i] * pa_z[i];

        g_yyyyz_zzzzzz_0[i] = 6.0 * g_yyyy_zzzzz_1[i] * fe_0 + g_yyyy_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 476-504 components of targeted buffer : HI

    auto g_yyyzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 476);

    auto g_yyyzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 477);

    auto g_yyyzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 478);

    auto g_yyyzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 479);

    auto g_yyyzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 480);

    auto g_yyyzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 481);

    auto g_yyyzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 482);

    auto g_yyyzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 483);

    auto g_yyyzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 484);

    auto g_yyyzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 485);

    auto g_yyyzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 486);

    auto g_yyyzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 487);

    auto g_yyyzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 488);

    auto g_yyyzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 489);

    auto g_yyyzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 490);

    auto g_yyyzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 491);

    auto g_yyyzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 492);

    auto g_yyyzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 493);

    auto g_yyyzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 494);

    auto g_yyyzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 495);

    auto g_yyyzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 496);

    auto g_yyyzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 497);

    auto g_yyyzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 498);

    auto g_yyyzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 499);

    auto g_yyyzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 500);

    auto g_yyyzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 501);

    auto g_yyyzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 502);

    auto g_yyyzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 503);

    #pragma omp simd aligned(g_yyy_xxxxxy_0, g_yyy_xxxxxy_1, g_yyy_xxxxyy_0, g_yyy_xxxxyy_1, g_yyy_xxxyyy_0, g_yyy_xxxyyy_1, g_yyy_xxyyyy_0, g_yyy_xxyyyy_1, g_yyy_xyyyyy_0, g_yyy_xyyyyy_1, g_yyy_yyyyyy_0, g_yyy_yyyyyy_1, g_yyyz_xxxxxy_1, g_yyyz_xxxxyy_1, g_yyyz_xxxyyy_1, g_yyyz_xxyyyy_1, g_yyyz_xyyyyy_1, g_yyyz_yyyyyy_1, g_yyyzz_xxxxxx_0, g_yyyzz_xxxxxy_0, g_yyyzz_xxxxxz_0, g_yyyzz_xxxxyy_0, g_yyyzz_xxxxyz_0, g_yyyzz_xxxxzz_0, g_yyyzz_xxxyyy_0, g_yyyzz_xxxyyz_0, g_yyyzz_xxxyzz_0, g_yyyzz_xxxzzz_0, g_yyyzz_xxyyyy_0, g_yyyzz_xxyyyz_0, g_yyyzz_xxyyzz_0, g_yyyzz_xxyzzz_0, g_yyyzz_xxzzzz_0, g_yyyzz_xyyyyy_0, g_yyyzz_xyyyyz_0, g_yyyzz_xyyyzz_0, g_yyyzz_xyyzzz_0, g_yyyzz_xyzzzz_0, g_yyyzz_xzzzzz_0, g_yyyzz_yyyyyy_0, g_yyyzz_yyyyyz_0, g_yyyzz_yyyyzz_0, g_yyyzz_yyyzzz_0, g_yyyzz_yyzzzz_0, g_yyyzz_yzzzzz_0, g_yyyzz_zzzzzz_0, g_yyzz_xxxxxx_1, g_yyzz_xxxxxz_1, g_yyzz_xxxxyz_1, g_yyzz_xxxxz_1, g_yyzz_xxxxzz_1, g_yyzz_xxxyyz_1, g_yyzz_xxxyz_1, g_yyzz_xxxyzz_1, g_yyzz_xxxzz_1, g_yyzz_xxxzzz_1, g_yyzz_xxyyyz_1, g_yyzz_xxyyz_1, g_yyzz_xxyyzz_1, g_yyzz_xxyzz_1, g_yyzz_xxyzzz_1, g_yyzz_xxzzz_1, g_yyzz_xxzzzz_1, g_yyzz_xyyyyz_1, g_yyzz_xyyyz_1, g_yyzz_xyyyzz_1, g_yyzz_xyyzz_1, g_yyzz_xyyzzz_1, g_yyzz_xyzzz_1, g_yyzz_xyzzzz_1, g_yyzz_xzzzz_1, g_yyzz_xzzzzz_1, g_yyzz_yyyyyz_1, g_yyzz_yyyyz_1, g_yyzz_yyyyzz_1, g_yyzz_yyyzz_1, g_yyzz_yyyzzz_1, g_yyzz_yyzzz_1, g_yyzz_yyzzzz_1, g_yyzz_yzzzz_1, g_yyzz_yzzzzz_1, g_yyzz_zzzzz_1, g_yyzz_zzzzzz_1, g_yzz_xxxxxx_0, g_yzz_xxxxxx_1, g_yzz_xxxxxz_0, g_yzz_xxxxxz_1, g_yzz_xxxxyz_0, g_yzz_xxxxyz_1, g_yzz_xxxxzz_0, g_yzz_xxxxzz_1, g_yzz_xxxyyz_0, g_yzz_xxxyyz_1, g_yzz_xxxyzz_0, g_yzz_xxxyzz_1, g_yzz_xxxzzz_0, g_yzz_xxxzzz_1, g_yzz_xxyyyz_0, g_yzz_xxyyyz_1, g_yzz_xxyyzz_0, g_yzz_xxyyzz_1, g_yzz_xxyzzz_0, g_yzz_xxyzzz_1, g_yzz_xxzzzz_0, g_yzz_xxzzzz_1, g_yzz_xyyyyz_0, g_yzz_xyyyyz_1, g_yzz_xyyyzz_0, g_yzz_xyyyzz_1, g_yzz_xyyzzz_0, g_yzz_xyyzzz_1, g_yzz_xyzzzz_0, g_yzz_xyzzzz_1, g_yzz_xzzzzz_0, g_yzz_xzzzzz_1, g_yzz_yyyyyz_0, g_yzz_yyyyyz_1, g_yzz_yyyyzz_0, g_yzz_yyyyzz_1, g_yzz_yyyzzz_0, g_yzz_yyyzzz_1, g_yzz_yyzzzz_0, g_yzz_yyzzzz_1, g_yzz_yzzzzz_0, g_yzz_yzzzzz_1, g_yzz_zzzzzz_0, g_yzz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzz_xxxxxx_0[i] = 2.0 * g_yzz_xxxxxx_0[i] * fbe_0 - 2.0 * g_yzz_xxxxxx_1[i] * fz_be_0 + g_yyzz_xxxxxx_1[i] * pa_y[i];

        g_yyyzz_xxxxxy_0[i] = g_yyy_xxxxxy_0[i] * fbe_0 - g_yyy_xxxxxy_1[i] * fz_be_0 + g_yyyz_xxxxxy_1[i] * pa_z[i];

        g_yyyzz_xxxxxz_0[i] = 2.0 * g_yzz_xxxxxz_0[i] * fbe_0 - 2.0 * g_yzz_xxxxxz_1[i] * fz_be_0 + g_yyzz_xxxxxz_1[i] * pa_y[i];

        g_yyyzz_xxxxyy_0[i] = g_yyy_xxxxyy_0[i] * fbe_0 - g_yyy_xxxxyy_1[i] * fz_be_0 + g_yyyz_xxxxyy_1[i] * pa_z[i];

        g_yyyzz_xxxxyz_0[i] = 2.0 * g_yzz_xxxxyz_0[i] * fbe_0 - 2.0 * g_yzz_xxxxyz_1[i] * fz_be_0 + g_yyzz_xxxxz_1[i] * fe_0 + g_yyzz_xxxxyz_1[i] * pa_y[i];

        g_yyyzz_xxxxzz_0[i] = 2.0 * g_yzz_xxxxzz_0[i] * fbe_0 - 2.0 * g_yzz_xxxxzz_1[i] * fz_be_0 + g_yyzz_xxxxzz_1[i] * pa_y[i];

        g_yyyzz_xxxyyy_0[i] = g_yyy_xxxyyy_0[i] * fbe_0 - g_yyy_xxxyyy_1[i] * fz_be_0 + g_yyyz_xxxyyy_1[i] * pa_z[i];

        g_yyyzz_xxxyyz_0[i] = 2.0 * g_yzz_xxxyyz_0[i] * fbe_0 - 2.0 * g_yzz_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyzz_xxxyz_1[i] * fe_0 + g_yyzz_xxxyyz_1[i] * pa_y[i];

        g_yyyzz_xxxyzz_0[i] = 2.0 * g_yzz_xxxyzz_0[i] * fbe_0 - 2.0 * g_yzz_xxxyzz_1[i] * fz_be_0 + g_yyzz_xxxzz_1[i] * fe_0 + g_yyzz_xxxyzz_1[i] * pa_y[i];

        g_yyyzz_xxxzzz_0[i] = 2.0 * g_yzz_xxxzzz_0[i] * fbe_0 - 2.0 * g_yzz_xxxzzz_1[i] * fz_be_0 + g_yyzz_xxxzzz_1[i] * pa_y[i];

        g_yyyzz_xxyyyy_0[i] = g_yyy_xxyyyy_0[i] * fbe_0 - g_yyy_xxyyyy_1[i] * fz_be_0 + g_yyyz_xxyyyy_1[i] * pa_z[i];

        g_yyyzz_xxyyyz_0[i] = 2.0 * g_yzz_xxyyyz_0[i] * fbe_0 - 2.0 * g_yzz_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyzz_xxyyz_1[i] * fe_0 + g_yyzz_xxyyyz_1[i] * pa_y[i];

        g_yyyzz_xxyyzz_0[i] = 2.0 * g_yzz_xxyyzz_0[i] * fbe_0 - 2.0 * g_yzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyzz_xxyzz_1[i] * fe_0 + g_yyzz_xxyyzz_1[i] * pa_y[i];

        g_yyyzz_xxyzzz_0[i] = 2.0 * g_yzz_xxyzzz_0[i] * fbe_0 - 2.0 * g_yzz_xxyzzz_1[i] * fz_be_0 + g_yyzz_xxzzz_1[i] * fe_0 + g_yyzz_xxyzzz_1[i] * pa_y[i];

        g_yyyzz_xxzzzz_0[i] = 2.0 * g_yzz_xxzzzz_0[i] * fbe_0 - 2.0 * g_yzz_xxzzzz_1[i] * fz_be_0 + g_yyzz_xxzzzz_1[i] * pa_y[i];

        g_yyyzz_xyyyyy_0[i] = g_yyy_xyyyyy_0[i] * fbe_0 - g_yyy_xyyyyy_1[i] * fz_be_0 + g_yyyz_xyyyyy_1[i] * pa_z[i];

        g_yyyzz_xyyyyz_0[i] = 2.0 * g_yzz_xyyyyz_0[i] * fbe_0 - 2.0 * g_yzz_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyzz_xyyyz_1[i] * fe_0 + g_yyzz_xyyyyz_1[i] * pa_y[i];

        g_yyyzz_xyyyzz_0[i] = 2.0 * g_yzz_xyyyzz_0[i] * fbe_0 - 2.0 * g_yzz_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyzz_xyyzz_1[i] * fe_0 + g_yyzz_xyyyzz_1[i] * pa_y[i];

        g_yyyzz_xyyzzz_0[i] = 2.0 * g_yzz_xyyzzz_0[i] * fbe_0 - 2.0 * g_yzz_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_xyzzz_1[i] * fe_0 + g_yyzz_xyyzzz_1[i] * pa_y[i];

        g_yyyzz_xyzzzz_0[i] = 2.0 * g_yzz_xyzzzz_0[i] * fbe_0 - 2.0 * g_yzz_xyzzzz_1[i] * fz_be_0 + g_yyzz_xzzzz_1[i] * fe_0 + g_yyzz_xyzzzz_1[i] * pa_y[i];

        g_yyyzz_xzzzzz_0[i] = 2.0 * g_yzz_xzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_xzzzzz_1[i] * fz_be_0 + g_yyzz_xzzzzz_1[i] * pa_y[i];

        g_yyyzz_yyyyyy_0[i] = g_yyy_yyyyyy_0[i] * fbe_0 - g_yyy_yyyyyy_1[i] * fz_be_0 + g_yyyz_yyyyyy_1[i] * pa_z[i];

        g_yyyzz_yyyyyz_0[i] = 2.0 * g_yzz_yyyyyz_0[i] * fbe_0 - 2.0 * g_yzz_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyzz_yyyyz_1[i] * fe_0 + g_yyzz_yyyyyz_1[i] * pa_y[i];

        g_yyyzz_yyyyzz_0[i] = 2.0 * g_yzz_yyyyzz_0[i] * fbe_0 - 2.0 * g_yzz_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyzz_yyyzz_1[i] * fe_0 + g_yyzz_yyyyzz_1[i] * pa_y[i];

        g_yyyzz_yyyzzz_0[i] = 2.0 * g_yzz_yyyzzz_0[i] * fbe_0 - 2.0 * g_yzz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyzz_yyzzz_1[i] * fe_0 + g_yyzz_yyyzzz_1[i] * pa_y[i];

        g_yyyzz_yyzzzz_0[i] = 2.0 * g_yzz_yyzzzz_0[i] * fbe_0 - 2.0 * g_yzz_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyzz_yzzzz_1[i] * fe_0 + g_yyzz_yyzzzz_1[i] * pa_y[i];

        g_yyyzz_yzzzzz_0[i] = 2.0 * g_yzz_yzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_yzzzzz_1[i] * fz_be_0 + g_yyzz_zzzzz_1[i] * fe_0 + g_yyzz_yzzzzz_1[i] * pa_y[i];

        g_yyyzz_zzzzzz_0[i] = 2.0 * g_yzz_zzzzzz_0[i] * fbe_0 - 2.0 * g_yzz_zzzzzz_1[i] * fz_be_0 + g_yyzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 504-532 components of targeted buffer : HI

    auto g_yyzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 504);

    auto g_yyzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 505);

    auto g_yyzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 506);

    auto g_yyzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 507);

    auto g_yyzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 508);

    auto g_yyzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 509);

    auto g_yyzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 510);

    auto g_yyzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 511);

    auto g_yyzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 512);

    auto g_yyzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 513);

    auto g_yyzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 514);

    auto g_yyzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 515);

    auto g_yyzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 516);

    auto g_yyzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 517);

    auto g_yyzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 518);

    auto g_yyzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 519);

    auto g_yyzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 520);

    auto g_yyzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 521);

    auto g_yyzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 522);

    auto g_yyzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 523);

    auto g_yyzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 524);

    auto g_yyzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 525);

    auto g_yyzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 526);

    auto g_yyzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 527);

    auto g_yyzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 528);

    auto g_yyzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 529);

    auto g_yyzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 530);

    auto g_yyzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 531);

    #pragma omp simd aligned(g_yyz_xxxxxy_0, g_yyz_xxxxxy_1, g_yyz_xxxxyy_0, g_yyz_xxxxyy_1, g_yyz_xxxyyy_0, g_yyz_xxxyyy_1, g_yyz_xxyyyy_0, g_yyz_xxyyyy_1, g_yyz_xyyyyy_0, g_yyz_xyyyyy_1, g_yyz_yyyyyy_0, g_yyz_yyyyyy_1, g_yyzz_xxxxxy_1, g_yyzz_xxxxyy_1, g_yyzz_xxxyyy_1, g_yyzz_xxyyyy_1, g_yyzz_xyyyyy_1, g_yyzz_yyyyyy_1, g_yyzzz_xxxxxx_0, g_yyzzz_xxxxxy_0, g_yyzzz_xxxxxz_0, g_yyzzz_xxxxyy_0, g_yyzzz_xxxxyz_0, g_yyzzz_xxxxzz_0, g_yyzzz_xxxyyy_0, g_yyzzz_xxxyyz_0, g_yyzzz_xxxyzz_0, g_yyzzz_xxxzzz_0, g_yyzzz_xxyyyy_0, g_yyzzz_xxyyyz_0, g_yyzzz_xxyyzz_0, g_yyzzz_xxyzzz_0, g_yyzzz_xxzzzz_0, g_yyzzz_xyyyyy_0, g_yyzzz_xyyyyz_0, g_yyzzz_xyyyzz_0, g_yyzzz_xyyzzz_0, g_yyzzz_xyzzzz_0, g_yyzzz_xzzzzz_0, g_yyzzz_yyyyyy_0, g_yyzzz_yyyyyz_0, g_yyzzz_yyyyzz_0, g_yyzzz_yyyzzz_0, g_yyzzz_yyzzzz_0, g_yyzzz_yzzzzz_0, g_yyzzz_zzzzzz_0, g_yzzz_xxxxxx_1, g_yzzz_xxxxxz_1, g_yzzz_xxxxyz_1, g_yzzz_xxxxz_1, g_yzzz_xxxxzz_1, g_yzzz_xxxyyz_1, g_yzzz_xxxyz_1, g_yzzz_xxxyzz_1, g_yzzz_xxxzz_1, g_yzzz_xxxzzz_1, g_yzzz_xxyyyz_1, g_yzzz_xxyyz_1, g_yzzz_xxyyzz_1, g_yzzz_xxyzz_1, g_yzzz_xxyzzz_1, g_yzzz_xxzzz_1, g_yzzz_xxzzzz_1, g_yzzz_xyyyyz_1, g_yzzz_xyyyz_1, g_yzzz_xyyyzz_1, g_yzzz_xyyzz_1, g_yzzz_xyyzzz_1, g_yzzz_xyzzz_1, g_yzzz_xyzzzz_1, g_yzzz_xzzzz_1, g_yzzz_xzzzzz_1, g_yzzz_yyyyyz_1, g_yzzz_yyyyz_1, g_yzzz_yyyyzz_1, g_yzzz_yyyzz_1, g_yzzz_yyyzzz_1, g_yzzz_yyzzz_1, g_yzzz_yyzzzz_1, g_yzzz_yzzzz_1, g_yzzz_yzzzzz_1, g_yzzz_zzzzz_1, g_yzzz_zzzzzz_1, g_zzz_xxxxxx_0, g_zzz_xxxxxx_1, g_zzz_xxxxxz_0, g_zzz_xxxxxz_1, g_zzz_xxxxyz_0, g_zzz_xxxxyz_1, g_zzz_xxxxzz_0, g_zzz_xxxxzz_1, g_zzz_xxxyyz_0, g_zzz_xxxyyz_1, g_zzz_xxxyzz_0, g_zzz_xxxyzz_1, g_zzz_xxxzzz_0, g_zzz_xxxzzz_1, g_zzz_xxyyyz_0, g_zzz_xxyyyz_1, g_zzz_xxyyzz_0, g_zzz_xxyyzz_1, g_zzz_xxyzzz_0, g_zzz_xxyzzz_1, g_zzz_xxzzzz_0, g_zzz_xxzzzz_1, g_zzz_xyyyyz_0, g_zzz_xyyyyz_1, g_zzz_xyyyzz_0, g_zzz_xyyyzz_1, g_zzz_xyyzzz_0, g_zzz_xyyzzz_1, g_zzz_xyzzzz_0, g_zzz_xyzzzz_1, g_zzz_xzzzzz_0, g_zzz_xzzzzz_1, g_zzz_yyyyyz_0, g_zzz_yyyyyz_1, g_zzz_yyyyzz_0, g_zzz_yyyyzz_1, g_zzz_yyyzzz_0, g_zzz_yyyzzz_1, g_zzz_yyzzzz_0, g_zzz_yyzzzz_1, g_zzz_yzzzzz_0, g_zzz_yzzzzz_1, g_zzz_zzzzzz_0, g_zzz_zzzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzz_xxxxxx_0[i] = g_zzz_xxxxxx_0[i] * fbe_0 - g_zzz_xxxxxx_1[i] * fz_be_0 + g_yzzz_xxxxxx_1[i] * pa_y[i];

        g_yyzzz_xxxxxy_0[i] = 2.0 * g_yyz_xxxxxy_0[i] * fbe_0 - 2.0 * g_yyz_xxxxxy_1[i] * fz_be_0 + g_yyzz_xxxxxy_1[i] * pa_z[i];

        g_yyzzz_xxxxxz_0[i] = g_zzz_xxxxxz_0[i] * fbe_0 - g_zzz_xxxxxz_1[i] * fz_be_0 + g_yzzz_xxxxxz_1[i] * pa_y[i];

        g_yyzzz_xxxxyy_0[i] = 2.0 * g_yyz_xxxxyy_0[i] * fbe_0 - 2.0 * g_yyz_xxxxyy_1[i] * fz_be_0 + g_yyzz_xxxxyy_1[i] * pa_z[i];

        g_yyzzz_xxxxyz_0[i] = g_zzz_xxxxyz_0[i] * fbe_0 - g_zzz_xxxxyz_1[i] * fz_be_0 + g_yzzz_xxxxz_1[i] * fe_0 + g_yzzz_xxxxyz_1[i] * pa_y[i];

        g_yyzzz_xxxxzz_0[i] = g_zzz_xxxxzz_0[i] * fbe_0 - g_zzz_xxxxzz_1[i] * fz_be_0 + g_yzzz_xxxxzz_1[i] * pa_y[i];

        g_yyzzz_xxxyyy_0[i] = 2.0 * g_yyz_xxxyyy_0[i] * fbe_0 - 2.0 * g_yyz_xxxyyy_1[i] * fz_be_0 + g_yyzz_xxxyyy_1[i] * pa_z[i];

        g_yyzzz_xxxyyz_0[i] = g_zzz_xxxyyz_0[i] * fbe_0 - g_zzz_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yzzz_xxxyz_1[i] * fe_0 + g_yzzz_xxxyyz_1[i] * pa_y[i];

        g_yyzzz_xxxyzz_0[i] = g_zzz_xxxyzz_0[i] * fbe_0 - g_zzz_xxxyzz_1[i] * fz_be_0 + g_yzzz_xxxzz_1[i] * fe_0 + g_yzzz_xxxyzz_1[i] * pa_y[i];

        g_yyzzz_xxxzzz_0[i] = g_zzz_xxxzzz_0[i] * fbe_0 - g_zzz_xxxzzz_1[i] * fz_be_0 + g_yzzz_xxxzzz_1[i] * pa_y[i];

        g_yyzzz_xxyyyy_0[i] = 2.0 * g_yyz_xxyyyy_0[i] * fbe_0 - 2.0 * g_yyz_xxyyyy_1[i] * fz_be_0 + g_yyzz_xxyyyy_1[i] * pa_z[i];

        g_yyzzz_xxyyyz_0[i] = g_zzz_xxyyyz_0[i] * fbe_0 - g_zzz_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yzzz_xxyyz_1[i] * fe_0 + g_yzzz_xxyyyz_1[i] * pa_y[i];

        g_yyzzz_xxyyzz_0[i] = g_zzz_xxyyzz_0[i] * fbe_0 - g_zzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yzzz_xxyzz_1[i] * fe_0 + g_yzzz_xxyyzz_1[i] * pa_y[i];

        g_yyzzz_xxyzzz_0[i] = g_zzz_xxyzzz_0[i] * fbe_0 - g_zzz_xxyzzz_1[i] * fz_be_0 + g_yzzz_xxzzz_1[i] * fe_0 + g_yzzz_xxyzzz_1[i] * pa_y[i];

        g_yyzzz_xxzzzz_0[i] = g_zzz_xxzzzz_0[i] * fbe_0 - g_zzz_xxzzzz_1[i] * fz_be_0 + g_yzzz_xxzzzz_1[i] * pa_y[i];

        g_yyzzz_xyyyyy_0[i] = 2.0 * g_yyz_xyyyyy_0[i] * fbe_0 - 2.0 * g_yyz_xyyyyy_1[i] * fz_be_0 + g_yyzz_xyyyyy_1[i] * pa_z[i];

        g_yyzzz_xyyyyz_0[i] = g_zzz_xyyyyz_0[i] * fbe_0 - g_zzz_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yzzz_xyyyz_1[i] * fe_0 + g_yzzz_xyyyyz_1[i] * pa_y[i];

        g_yyzzz_xyyyzz_0[i] = g_zzz_xyyyzz_0[i] * fbe_0 - g_zzz_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yzzz_xyyzz_1[i] * fe_0 + g_yzzz_xyyyzz_1[i] * pa_y[i];

        g_yyzzz_xyyzzz_0[i] = g_zzz_xyyzzz_0[i] * fbe_0 - g_zzz_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_xyzzz_1[i] * fe_0 + g_yzzz_xyyzzz_1[i] * pa_y[i];

        g_yyzzz_xyzzzz_0[i] = g_zzz_xyzzzz_0[i] * fbe_0 - g_zzz_xyzzzz_1[i] * fz_be_0 + g_yzzz_xzzzz_1[i] * fe_0 + g_yzzz_xyzzzz_1[i] * pa_y[i];

        g_yyzzz_xzzzzz_0[i] = g_zzz_xzzzzz_0[i] * fbe_0 - g_zzz_xzzzzz_1[i] * fz_be_0 + g_yzzz_xzzzzz_1[i] * pa_y[i];

        g_yyzzz_yyyyyy_0[i] = 2.0 * g_yyz_yyyyyy_0[i] * fbe_0 - 2.0 * g_yyz_yyyyyy_1[i] * fz_be_0 + g_yyzz_yyyyyy_1[i] * pa_z[i];

        g_yyzzz_yyyyyz_0[i] = g_zzz_yyyyyz_0[i] * fbe_0 - g_zzz_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yzzz_yyyyz_1[i] * fe_0 + g_yzzz_yyyyyz_1[i] * pa_y[i];

        g_yyzzz_yyyyzz_0[i] = g_zzz_yyyyzz_0[i] * fbe_0 - g_zzz_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yzzz_yyyzz_1[i] * fe_0 + g_yzzz_yyyyzz_1[i] * pa_y[i];

        g_yyzzz_yyyzzz_0[i] = g_zzz_yyyzzz_0[i] * fbe_0 - g_zzz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yzzz_yyzzz_1[i] * fe_0 + g_yzzz_yyyzzz_1[i] * pa_y[i];

        g_yyzzz_yyzzzz_0[i] = g_zzz_yyzzzz_0[i] * fbe_0 - g_zzz_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yzzz_yzzzz_1[i] * fe_0 + g_yzzz_yyzzzz_1[i] * pa_y[i];

        g_yyzzz_yzzzzz_0[i] = g_zzz_yzzzzz_0[i] * fbe_0 - g_zzz_yzzzzz_1[i] * fz_be_0 + g_yzzz_zzzzz_1[i] * fe_0 + g_yzzz_yzzzzz_1[i] * pa_y[i];

        g_yyzzz_zzzzzz_0[i] = g_zzz_zzzzzz_0[i] * fbe_0 - g_zzz_zzzzzz_1[i] * fz_be_0 + g_yzzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 532-560 components of targeted buffer : HI

    auto g_yzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 532);

    auto g_yzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 533);

    auto g_yzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 534);

    auto g_yzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 535);

    auto g_yzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 536);

    auto g_yzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 537);

    auto g_yzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 538);

    auto g_yzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 539);

    auto g_yzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 540);

    auto g_yzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 541);

    auto g_yzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 542);

    auto g_yzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 543);

    auto g_yzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 544);

    auto g_yzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 545);

    auto g_yzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 546);

    auto g_yzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 547);

    auto g_yzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 548);

    auto g_yzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 549);

    auto g_yzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 550);

    auto g_yzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 551);

    auto g_yzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 552);

    auto g_yzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 553);

    auto g_yzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 554);

    auto g_yzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 555);

    auto g_yzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 556);

    auto g_yzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 557);

    auto g_yzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 558);

    auto g_yzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 559);

    #pragma omp simd aligned(g_yzzzz_xxxxxx_0, g_yzzzz_xxxxxy_0, g_yzzzz_xxxxxz_0, g_yzzzz_xxxxyy_0, g_yzzzz_xxxxyz_0, g_yzzzz_xxxxzz_0, g_yzzzz_xxxyyy_0, g_yzzzz_xxxyyz_0, g_yzzzz_xxxyzz_0, g_yzzzz_xxxzzz_0, g_yzzzz_xxyyyy_0, g_yzzzz_xxyyyz_0, g_yzzzz_xxyyzz_0, g_yzzzz_xxyzzz_0, g_yzzzz_xxzzzz_0, g_yzzzz_xyyyyy_0, g_yzzzz_xyyyyz_0, g_yzzzz_xyyyzz_0, g_yzzzz_xyyzzz_0, g_yzzzz_xyzzzz_0, g_yzzzz_xzzzzz_0, g_yzzzz_yyyyyy_0, g_yzzzz_yyyyyz_0, g_yzzzz_yyyyzz_0, g_yzzzz_yyyzzz_0, g_yzzzz_yyzzzz_0, g_yzzzz_yzzzzz_0, g_yzzzz_zzzzzz_0, g_zzzz_xxxxx_1, g_zzzz_xxxxxx_1, g_zzzz_xxxxxy_1, g_zzzz_xxxxxz_1, g_zzzz_xxxxy_1, g_zzzz_xxxxyy_1, g_zzzz_xxxxyz_1, g_zzzz_xxxxz_1, g_zzzz_xxxxzz_1, g_zzzz_xxxyy_1, g_zzzz_xxxyyy_1, g_zzzz_xxxyyz_1, g_zzzz_xxxyz_1, g_zzzz_xxxyzz_1, g_zzzz_xxxzz_1, g_zzzz_xxxzzz_1, g_zzzz_xxyyy_1, g_zzzz_xxyyyy_1, g_zzzz_xxyyyz_1, g_zzzz_xxyyz_1, g_zzzz_xxyyzz_1, g_zzzz_xxyzz_1, g_zzzz_xxyzzz_1, g_zzzz_xxzzz_1, g_zzzz_xxzzzz_1, g_zzzz_xyyyy_1, g_zzzz_xyyyyy_1, g_zzzz_xyyyyz_1, g_zzzz_xyyyz_1, g_zzzz_xyyyzz_1, g_zzzz_xyyzz_1, g_zzzz_xyyzzz_1, g_zzzz_xyzzz_1, g_zzzz_xyzzzz_1, g_zzzz_xzzzz_1, g_zzzz_xzzzzz_1, g_zzzz_yyyyy_1, g_zzzz_yyyyyy_1, g_zzzz_yyyyyz_1, g_zzzz_yyyyz_1, g_zzzz_yyyyzz_1, g_zzzz_yyyzz_1, g_zzzz_yyyzzz_1, g_zzzz_yyzzz_1, g_zzzz_yyzzzz_1, g_zzzz_yzzzz_1, g_zzzz_yzzzzz_1, g_zzzz_zzzzz_1, g_zzzz_zzzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzz_xxxxxx_0[i] = g_zzzz_xxxxxx_1[i] * pa_y[i];

        g_yzzzz_xxxxxy_0[i] = g_zzzz_xxxxx_1[i] * fe_0 + g_zzzz_xxxxxy_1[i] * pa_y[i];

        g_yzzzz_xxxxxz_0[i] = g_zzzz_xxxxxz_1[i] * pa_y[i];

        g_yzzzz_xxxxyy_0[i] = 2.0 * g_zzzz_xxxxy_1[i] * fe_0 + g_zzzz_xxxxyy_1[i] * pa_y[i];

        g_yzzzz_xxxxyz_0[i] = g_zzzz_xxxxz_1[i] * fe_0 + g_zzzz_xxxxyz_1[i] * pa_y[i];

        g_yzzzz_xxxxzz_0[i] = g_zzzz_xxxxzz_1[i] * pa_y[i];

        g_yzzzz_xxxyyy_0[i] = 3.0 * g_zzzz_xxxyy_1[i] * fe_0 + g_zzzz_xxxyyy_1[i] * pa_y[i];

        g_yzzzz_xxxyyz_0[i] = 2.0 * g_zzzz_xxxyz_1[i] * fe_0 + g_zzzz_xxxyyz_1[i] * pa_y[i];

        g_yzzzz_xxxyzz_0[i] = g_zzzz_xxxzz_1[i] * fe_0 + g_zzzz_xxxyzz_1[i] * pa_y[i];

        g_yzzzz_xxxzzz_0[i] = g_zzzz_xxxzzz_1[i] * pa_y[i];

        g_yzzzz_xxyyyy_0[i] = 4.0 * g_zzzz_xxyyy_1[i] * fe_0 + g_zzzz_xxyyyy_1[i] * pa_y[i];

        g_yzzzz_xxyyyz_0[i] = 3.0 * g_zzzz_xxyyz_1[i] * fe_0 + g_zzzz_xxyyyz_1[i] * pa_y[i];

        g_yzzzz_xxyyzz_0[i] = 2.0 * g_zzzz_xxyzz_1[i] * fe_0 + g_zzzz_xxyyzz_1[i] * pa_y[i];

        g_yzzzz_xxyzzz_0[i] = g_zzzz_xxzzz_1[i] * fe_0 + g_zzzz_xxyzzz_1[i] * pa_y[i];

        g_yzzzz_xxzzzz_0[i] = g_zzzz_xxzzzz_1[i] * pa_y[i];

        g_yzzzz_xyyyyy_0[i] = 5.0 * g_zzzz_xyyyy_1[i] * fe_0 + g_zzzz_xyyyyy_1[i] * pa_y[i];

        g_yzzzz_xyyyyz_0[i] = 4.0 * g_zzzz_xyyyz_1[i] * fe_0 + g_zzzz_xyyyyz_1[i] * pa_y[i];

        g_yzzzz_xyyyzz_0[i] = 3.0 * g_zzzz_xyyzz_1[i] * fe_0 + g_zzzz_xyyyzz_1[i] * pa_y[i];

        g_yzzzz_xyyzzz_0[i] = 2.0 * g_zzzz_xyzzz_1[i] * fe_0 + g_zzzz_xyyzzz_1[i] * pa_y[i];

        g_yzzzz_xyzzzz_0[i] = g_zzzz_xzzzz_1[i] * fe_0 + g_zzzz_xyzzzz_1[i] * pa_y[i];

        g_yzzzz_xzzzzz_0[i] = g_zzzz_xzzzzz_1[i] * pa_y[i];

        g_yzzzz_yyyyyy_0[i] = 6.0 * g_zzzz_yyyyy_1[i] * fe_0 + g_zzzz_yyyyyy_1[i] * pa_y[i];

        g_yzzzz_yyyyyz_0[i] = 5.0 * g_zzzz_yyyyz_1[i] * fe_0 + g_zzzz_yyyyyz_1[i] * pa_y[i];

        g_yzzzz_yyyyzz_0[i] = 4.0 * g_zzzz_yyyzz_1[i] * fe_0 + g_zzzz_yyyyzz_1[i] * pa_y[i];

        g_yzzzz_yyyzzz_0[i] = 3.0 * g_zzzz_yyzzz_1[i] * fe_0 + g_zzzz_yyyzzz_1[i] * pa_y[i];

        g_yzzzz_yyzzzz_0[i] = 2.0 * g_zzzz_yzzzz_1[i] * fe_0 + g_zzzz_yyzzzz_1[i] * pa_y[i];

        g_yzzzz_yzzzzz_0[i] = g_zzzz_zzzzz_1[i] * fe_0 + g_zzzz_yzzzzz_1[i] * pa_y[i];

        g_yzzzz_zzzzzz_0[i] = g_zzzz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 560-588 components of targeted buffer : HI

    auto g_zzzzz_xxxxxx_0 = pbuffer.data(idx_eri_0_hi + 560);

    auto g_zzzzz_xxxxxy_0 = pbuffer.data(idx_eri_0_hi + 561);

    auto g_zzzzz_xxxxxz_0 = pbuffer.data(idx_eri_0_hi + 562);

    auto g_zzzzz_xxxxyy_0 = pbuffer.data(idx_eri_0_hi + 563);

    auto g_zzzzz_xxxxyz_0 = pbuffer.data(idx_eri_0_hi + 564);

    auto g_zzzzz_xxxxzz_0 = pbuffer.data(idx_eri_0_hi + 565);

    auto g_zzzzz_xxxyyy_0 = pbuffer.data(idx_eri_0_hi + 566);

    auto g_zzzzz_xxxyyz_0 = pbuffer.data(idx_eri_0_hi + 567);

    auto g_zzzzz_xxxyzz_0 = pbuffer.data(idx_eri_0_hi + 568);

    auto g_zzzzz_xxxzzz_0 = pbuffer.data(idx_eri_0_hi + 569);

    auto g_zzzzz_xxyyyy_0 = pbuffer.data(idx_eri_0_hi + 570);

    auto g_zzzzz_xxyyyz_0 = pbuffer.data(idx_eri_0_hi + 571);

    auto g_zzzzz_xxyyzz_0 = pbuffer.data(idx_eri_0_hi + 572);

    auto g_zzzzz_xxyzzz_0 = pbuffer.data(idx_eri_0_hi + 573);

    auto g_zzzzz_xxzzzz_0 = pbuffer.data(idx_eri_0_hi + 574);

    auto g_zzzzz_xyyyyy_0 = pbuffer.data(idx_eri_0_hi + 575);

    auto g_zzzzz_xyyyyz_0 = pbuffer.data(idx_eri_0_hi + 576);

    auto g_zzzzz_xyyyzz_0 = pbuffer.data(idx_eri_0_hi + 577);

    auto g_zzzzz_xyyzzz_0 = pbuffer.data(idx_eri_0_hi + 578);

    auto g_zzzzz_xyzzzz_0 = pbuffer.data(idx_eri_0_hi + 579);

    auto g_zzzzz_xzzzzz_0 = pbuffer.data(idx_eri_0_hi + 580);

    auto g_zzzzz_yyyyyy_0 = pbuffer.data(idx_eri_0_hi + 581);

    auto g_zzzzz_yyyyyz_0 = pbuffer.data(idx_eri_0_hi + 582);

    auto g_zzzzz_yyyyzz_0 = pbuffer.data(idx_eri_0_hi + 583);

    auto g_zzzzz_yyyzzz_0 = pbuffer.data(idx_eri_0_hi + 584);

    auto g_zzzzz_yyzzzz_0 = pbuffer.data(idx_eri_0_hi + 585);

    auto g_zzzzz_yzzzzz_0 = pbuffer.data(idx_eri_0_hi + 586);

    auto g_zzzzz_zzzzzz_0 = pbuffer.data(idx_eri_0_hi + 587);

    #pragma omp simd aligned(g_zzz_xxxxxx_0, g_zzz_xxxxxx_1, g_zzz_xxxxxy_0, g_zzz_xxxxxy_1, g_zzz_xxxxxz_0, g_zzz_xxxxxz_1, g_zzz_xxxxyy_0, g_zzz_xxxxyy_1, g_zzz_xxxxyz_0, g_zzz_xxxxyz_1, g_zzz_xxxxzz_0, g_zzz_xxxxzz_1, g_zzz_xxxyyy_0, g_zzz_xxxyyy_1, g_zzz_xxxyyz_0, g_zzz_xxxyyz_1, g_zzz_xxxyzz_0, g_zzz_xxxyzz_1, g_zzz_xxxzzz_0, g_zzz_xxxzzz_1, g_zzz_xxyyyy_0, g_zzz_xxyyyy_1, g_zzz_xxyyyz_0, g_zzz_xxyyyz_1, g_zzz_xxyyzz_0, g_zzz_xxyyzz_1, g_zzz_xxyzzz_0, g_zzz_xxyzzz_1, g_zzz_xxzzzz_0, g_zzz_xxzzzz_1, g_zzz_xyyyyy_0, g_zzz_xyyyyy_1, g_zzz_xyyyyz_0, g_zzz_xyyyyz_1, g_zzz_xyyyzz_0, g_zzz_xyyyzz_1, g_zzz_xyyzzz_0, g_zzz_xyyzzz_1, g_zzz_xyzzzz_0, g_zzz_xyzzzz_1, g_zzz_xzzzzz_0, g_zzz_xzzzzz_1, g_zzz_yyyyyy_0, g_zzz_yyyyyy_1, g_zzz_yyyyyz_0, g_zzz_yyyyyz_1, g_zzz_yyyyzz_0, g_zzz_yyyyzz_1, g_zzz_yyyzzz_0, g_zzz_yyyzzz_1, g_zzz_yyzzzz_0, g_zzz_yyzzzz_1, g_zzz_yzzzzz_0, g_zzz_yzzzzz_1, g_zzz_zzzzzz_0, g_zzz_zzzzzz_1, g_zzzz_xxxxx_1, g_zzzz_xxxxxx_1, g_zzzz_xxxxxy_1, g_zzzz_xxxxxz_1, g_zzzz_xxxxy_1, g_zzzz_xxxxyy_1, g_zzzz_xxxxyz_1, g_zzzz_xxxxz_1, g_zzzz_xxxxzz_1, g_zzzz_xxxyy_1, g_zzzz_xxxyyy_1, g_zzzz_xxxyyz_1, g_zzzz_xxxyz_1, g_zzzz_xxxyzz_1, g_zzzz_xxxzz_1, g_zzzz_xxxzzz_1, g_zzzz_xxyyy_1, g_zzzz_xxyyyy_1, g_zzzz_xxyyyz_1, g_zzzz_xxyyz_1, g_zzzz_xxyyzz_1, g_zzzz_xxyzz_1, g_zzzz_xxyzzz_1, g_zzzz_xxzzz_1, g_zzzz_xxzzzz_1, g_zzzz_xyyyy_1, g_zzzz_xyyyyy_1, g_zzzz_xyyyyz_1, g_zzzz_xyyyz_1, g_zzzz_xyyyzz_1, g_zzzz_xyyzz_1, g_zzzz_xyyzzz_1, g_zzzz_xyzzz_1, g_zzzz_xyzzzz_1, g_zzzz_xzzzz_1, g_zzzz_xzzzzz_1, g_zzzz_yyyyy_1, g_zzzz_yyyyyy_1, g_zzzz_yyyyyz_1, g_zzzz_yyyyz_1, g_zzzz_yyyyzz_1, g_zzzz_yyyzz_1, g_zzzz_yyyzzz_1, g_zzzz_yyzzz_1, g_zzzz_yyzzzz_1, g_zzzz_yzzzz_1, g_zzzz_yzzzzz_1, g_zzzz_zzzzz_1, g_zzzz_zzzzzz_1, g_zzzzz_xxxxxx_0, g_zzzzz_xxxxxy_0, g_zzzzz_xxxxxz_0, g_zzzzz_xxxxyy_0, g_zzzzz_xxxxyz_0, g_zzzzz_xxxxzz_0, g_zzzzz_xxxyyy_0, g_zzzzz_xxxyyz_0, g_zzzzz_xxxyzz_0, g_zzzzz_xxxzzz_0, g_zzzzz_xxyyyy_0, g_zzzzz_xxyyyz_0, g_zzzzz_xxyyzz_0, g_zzzzz_xxyzzz_0, g_zzzzz_xxzzzz_0, g_zzzzz_xyyyyy_0, g_zzzzz_xyyyyz_0, g_zzzzz_xyyyzz_0, g_zzzzz_xyyzzz_0, g_zzzzz_xyzzzz_0, g_zzzzz_xzzzzz_0, g_zzzzz_yyyyyy_0, g_zzzzz_yyyyyz_0, g_zzzzz_yyyyzz_0, g_zzzzz_yyyzzz_0, g_zzzzz_yyzzzz_0, g_zzzzz_yzzzzz_0, g_zzzzz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzz_xxxxxx_0[i] = 4.0 * g_zzz_xxxxxx_0[i] * fbe_0 - 4.0 * g_zzz_xxxxxx_1[i] * fz_be_0 + g_zzzz_xxxxxx_1[i] * pa_z[i];

        g_zzzzz_xxxxxy_0[i] = 4.0 * g_zzz_xxxxxy_0[i] * fbe_0 - 4.0 * g_zzz_xxxxxy_1[i] * fz_be_0 + g_zzzz_xxxxxy_1[i] * pa_z[i];

        g_zzzzz_xxxxxz_0[i] = 4.0 * g_zzz_xxxxxz_0[i] * fbe_0 - 4.0 * g_zzz_xxxxxz_1[i] * fz_be_0 + g_zzzz_xxxxx_1[i] * fe_0 + g_zzzz_xxxxxz_1[i] * pa_z[i];

        g_zzzzz_xxxxyy_0[i] = 4.0 * g_zzz_xxxxyy_0[i] * fbe_0 - 4.0 * g_zzz_xxxxyy_1[i] * fz_be_0 + g_zzzz_xxxxyy_1[i] * pa_z[i];

        g_zzzzz_xxxxyz_0[i] = 4.0 * g_zzz_xxxxyz_0[i] * fbe_0 - 4.0 * g_zzz_xxxxyz_1[i] * fz_be_0 + g_zzzz_xxxxy_1[i] * fe_0 + g_zzzz_xxxxyz_1[i] * pa_z[i];

        g_zzzzz_xxxxzz_0[i] = 4.0 * g_zzz_xxxxzz_0[i] * fbe_0 - 4.0 * g_zzz_xxxxzz_1[i] * fz_be_0 + 2.0 * g_zzzz_xxxxz_1[i] * fe_0 + g_zzzz_xxxxzz_1[i] * pa_z[i];

        g_zzzzz_xxxyyy_0[i] = 4.0 * g_zzz_xxxyyy_0[i] * fbe_0 - 4.0 * g_zzz_xxxyyy_1[i] * fz_be_0 + g_zzzz_xxxyyy_1[i] * pa_z[i];

        g_zzzzz_xxxyyz_0[i] = 4.0 * g_zzz_xxxyyz_0[i] * fbe_0 - 4.0 * g_zzz_xxxyyz_1[i] * fz_be_0 + g_zzzz_xxxyy_1[i] * fe_0 + g_zzzz_xxxyyz_1[i] * pa_z[i];

        g_zzzzz_xxxyzz_0[i] = 4.0 * g_zzz_xxxyzz_0[i] * fbe_0 - 4.0 * g_zzz_xxxyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_xxxyz_1[i] * fe_0 + g_zzzz_xxxyzz_1[i] * pa_z[i];

        g_zzzzz_xxxzzz_0[i] = 4.0 * g_zzz_xxxzzz_0[i] * fbe_0 - 4.0 * g_zzz_xxxzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_xxxzz_1[i] * fe_0 + g_zzzz_xxxzzz_1[i] * pa_z[i];

        g_zzzzz_xxyyyy_0[i] = 4.0 * g_zzz_xxyyyy_0[i] * fbe_0 - 4.0 * g_zzz_xxyyyy_1[i] * fz_be_0 + g_zzzz_xxyyyy_1[i] * pa_z[i];

        g_zzzzz_xxyyyz_0[i] = 4.0 * g_zzz_xxyyyz_0[i] * fbe_0 - 4.0 * g_zzz_xxyyyz_1[i] * fz_be_0 + g_zzzz_xxyyy_1[i] * fe_0 + g_zzzz_xxyyyz_1[i] * pa_z[i];

        g_zzzzz_xxyyzz_0[i] = 4.0 * g_zzz_xxyyzz_0[i] * fbe_0 - 4.0 * g_zzz_xxyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_xxyyz_1[i] * fe_0 + g_zzzz_xxyyzz_1[i] * pa_z[i];

        g_zzzzz_xxyzzz_0[i] = 4.0 * g_zzz_xxyzzz_0[i] * fbe_0 - 4.0 * g_zzz_xxyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_xxyzz_1[i] * fe_0 + g_zzzz_xxyzzz_1[i] * pa_z[i];

        g_zzzzz_xxzzzz_0[i] = 4.0 * g_zzz_xxzzzz_0[i] * fbe_0 - 4.0 * g_zzz_xxzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_xxzzz_1[i] * fe_0 + g_zzzz_xxzzzz_1[i] * pa_z[i];

        g_zzzzz_xyyyyy_0[i] = 4.0 * g_zzz_xyyyyy_0[i] * fbe_0 - 4.0 * g_zzz_xyyyyy_1[i] * fz_be_0 + g_zzzz_xyyyyy_1[i] * pa_z[i];

        g_zzzzz_xyyyyz_0[i] = 4.0 * g_zzz_xyyyyz_0[i] * fbe_0 - 4.0 * g_zzz_xyyyyz_1[i] * fz_be_0 + g_zzzz_xyyyy_1[i] * fe_0 + g_zzzz_xyyyyz_1[i] * pa_z[i];

        g_zzzzz_xyyyzz_0[i] = 4.0 * g_zzz_xyyyzz_0[i] * fbe_0 - 4.0 * g_zzz_xyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_xyyyz_1[i] * fe_0 + g_zzzz_xyyyzz_1[i] * pa_z[i];

        g_zzzzz_xyyzzz_0[i] = 4.0 * g_zzz_xyyzzz_0[i] * fbe_0 - 4.0 * g_zzz_xyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_xyyzz_1[i] * fe_0 + g_zzzz_xyyzzz_1[i] * pa_z[i];

        g_zzzzz_xyzzzz_0[i] = 4.0 * g_zzz_xyzzzz_0[i] * fbe_0 - 4.0 * g_zzz_xyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_xyzzz_1[i] * fe_0 + g_zzzz_xyzzzz_1[i] * pa_z[i];

        g_zzzzz_xzzzzz_0[i] = 4.0 * g_zzz_xzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_xzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_xzzzz_1[i] * fe_0 + g_zzzz_xzzzzz_1[i] * pa_z[i];

        g_zzzzz_yyyyyy_0[i] = 4.0 * g_zzz_yyyyyy_0[i] * fbe_0 - 4.0 * g_zzz_yyyyyy_1[i] * fz_be_0 + g_zzzz_yyyyyy_1[i] * pa_z[i];

        g_zzzzz_yyyyyz_0[i] = 4.0 * g_zzz_yyyyyz_0[i] * fbe_0 - 4.0 * g_zzz_yyyyyz_1[i] * fz_be_0 + g_zzzz_yyyyy_1[i] * fe_0 + g_zzzz_yyyyyz_1[i] * pa_z[i];

        g_zzzzz_yyyyzz_0[i] = 4.0 * g_zzz_yyyyzz_0[i] * fbe_0 - 4.0 * g_zzz_yyyyzz_1[i] * fz_be_0 + 2.0 * g_zzzz_yyyyz_1[i] * fe_0 + g_zzzz_yyyyzz_1[i] * pa_z[i];

        g_zzzzz_yyyzzz_0[i] = 4.0 * g_zzz_yyyzzz_0[i] * fbe_0 - 4.0 * g_zzz_yyyzzz_1[i] * fz_be_0 + 3.0 * g_zzzz_yyyzz_1[i] * fe_0 + g_zzzz_yyyzzz_1[i] * pa_z[i];

        g_zzzzz_yyzzzz_0[i] = 4.0 * g_zzz_yyzzzz_0[i] * fbe_0 - 4.0 * g_zzz_yyzzzz_1[i] * fz_be_0 + 4.0 * g_zzzz_yyzzz_1[i] * fe_0 + g_zzzz_yyzzzz_1[i] * pa_z[i];

        g_zzzzz_yzzzzz_0[i] = 4.0 * g_zzz_yzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_yzzzzz_1[i] * fz_be_0 + 5.0 * g_zzzz_yzzzz_1[i] * fe_0 + g_zzzz_yzzzzz_1[i] * pa_z[i];

        g_zzzzz_zzzzzz_0[i] = 4.0 * g_zzz_zzzzzz_0[i] * fbe_0 - 4.0 * g_zzz_zzzzzz_1[i] * fz_be_0 + 6.0 * g_zzzz_zzzzz_1[i] * fe_0 + g_zzzz_zzzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

